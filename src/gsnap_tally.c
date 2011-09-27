static char rcsid[] = "$Id: gsnap_tally.c 43665 2011-07-26 20:48:15Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>		/* Needed to define pthread_t on Solaris */
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>		/* For strcpy */
#include <strings.h>		/* For rindex */
#include <ctype.h>
#include <math.h>


#include "except.h"
#include "mem.h"
#include "bool.h"
#include "genomicpos.h"
#include "complement.h"
#include "list.h"
#include "iit-read.h"
#include "interval.h"
#include "genome.h"

#ifdef SAM_INPUT
#include "samflags.h"
#include "samread.h"
#elif defined(BAM_INPUT)
#include "samflags.h"		/* For flags */
#include "bamread.h"
#include "parserange.h"
#else
#include "gsnapread.h"
#endif

#include "datadir.h"
#include "getopt.h"


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* trimming */
#ifdef DEBUG8
#define debug8(x) x
#else
#define debug8(x)
#endif


#define MATCH_SCORE 1
#define MISMATCH_SCORE -3

static char *user_genomedir = NULL;
static char *dbroot = NULL;
static char *dbversion = NULL;

static char *chromosome = NULL;

static bool print_refdetails_p = false;
static bool want_genotypes_p = false;

static bool print_wig_p = false;
static bool print_allele_counts_p = true;

static bool need_concordant_p = false;
static bool uniquep = false;
static int blocksize = 1000;

static int minimum_mapq = 0;

/* SAM/BAM files should be in Sanger format */
#ifdef SAM_INPUT
static int min_mlength = 0;
static int quality_score_adj = 33;  /* For Illumina, subtract 64.  For Sanger, subtract 33 */
#define DEFAULT_QUALITY 33+40  /* h */

#elif defined(BAM_INPUT)
static int min_mlength = 0;
static int quality_score_adj = 33;  /* For Illumina, subtract 64.  For Sanger, subtract 33 */
#define DEFAULT_QUALITY 33+40  /* h */

#else
static int min_readlength = 0;
static int quality_score_adj = 64;  /* For Illumina, subtract 64.  For Sanger, subtract 33 */
#define DEFAULT_QUALITY 64+40  /* h */

#endif
static int quality_score_constant = -1;


static Genomicpos_T haltpos = -1U;

static bool trimp = false;
static bool dibasep = false;

static char desired_strand = ' ';
static char opposite_strand = ' ';

static bool want_low_p = true;
static bool want_high_p = true;
static bool allow_translocations_p = false; /* concordant translocations */


typedef enum {FIRST, SECOND} Pairend_T;
typedef enum {ALL_SCORES, MISMATCH_HIGH_SCORES, HIGH_SCORES} Mode_T;

static Mode_T mode = HIGH_SCORES;


static struct option long_options[] = {
  /* Input options */
  {"dir", required_argument, 0, 'D'},	/* user_genomedir */
  {"db", required_argument, 0, 'd'}, /* dbroot */

  {"alldetails", no_argument, 0, 'A'}, /* print_refdetails_p */
  {"genotypes", no_argument, 0, 'G'}, /* want_genotypes_p */
  {"use-quality-constant", required_argument, 0, 0}, /* quality_score_constant */

  {"chr", required_argument, 0, 'c'}, /* chromosome */
  {"fwd", no_argument, 0, 'F'}, /* desired_strand, opposite_strand */
  {"rev", no_argument, 0, 'R'}, /* desired_strand, opposite_strand */

  {"side", required_argument, 0, 's'}, /* want_low_p, want_high_p */

  {"concordant", required_argument, 0, 'C'}, /* need_concordant_p */
  {"unique", required_argument, 0, 'U'}, /* uniquep */
  {"blocksize", required_argument, 0, 'b'}, /* blocksize */

  {"mode", required_argument, 0, 'a'}, /* mode */

  {"minimum-mapq", required_argument, 0, 'Q'}, /* minimum_mapq */

  {"halt", required_argument, 0, 'h'}, /* haltpos */

  {"trim", required_argument, 0, 'T'}, /* trimp */
  {"minlength", required_argument, 0, 'l'}, /* min_readlength */
  
  {"dibase", no_argument, 0, '2'}, /* dibasep */

  /* Help options */
  {"version", no_argument, 0, 0}, /* print_program_version */
  {"help", no_argument, 0, 0}, /* print_program_usage */
  {0, 0, 0, 0}
};


static void
print_program_version () {
  fprintf(stdout,"\n");
#ifdef SAM_INPUT
  fprintf(stdout,"SAM_TALLY\n");
#elif defined(BAM_INPUT)
  fprintf(stdout,"BAM_TALLY\n");
#else
  fprintf(stdout,"GSNAP_TALLY\n");
#endif
  fprintf(stdout,"Part of GMAP package, version %s\n",PACKAGE_VERSION);
  fprintf(stdout,"Build target: %s\n",TARGET);
  fprintf(stdout,"Default gmap directory: %s\n",GMAPDB);
  fprintf(stdout,"Thomas D. Wu, Genentech, Inc.\n");
  fprintf(stdout,"Contact: twu@gene.com\n");
  fprintf(stdout,"\n");
  return;
}


static void
print_program_usage ();

typedef struct Match_T *Match_T;
struct Match_T {
  char quality;
  int shift;			/* Used to record shifts */
  int count;
};

static Match_T
Match_new (int shift, char quality) {
  Match_T new = (Match_T) MALLOC(sizeof(*new));

  new->quality = quality;
  new->shift = shift;
  new->count = 1;
  return new;
}

static void
Match_free (Match_T *old) {
  FREE(*old);
  return;
}

static Match_T
find_match_byshift (List_T matches, int shift) {
  List_T p;
  Match_T match;

  for (p = matches; p != NULL; p = List_next(p)) {
    match = (Match_T) List_head(p);
    if (match->shift == shift) {
      return match;
    }
  }
  return (Match_T) NULL;
}

static Match_T
find_match_byquality (List_T matches, char quality) {
  List_T p;
  Match_T match;

  for (p = matches; p != NULL; p = List_next(p)) {
    match = (Match_T) List_head(p);
    if (match->quality == quality) {
      return match;
    }
  }
  return (Match_T) NULL;
}

/* Go -1 to -readlength, then +readlength to +1 */
static int
Match_shift_cmp (const void *a, const void *b) {
  Match_T x = * (Match_T *) a;
  Match_T y = * (Match_T *) b;

  if (x->shift < 0 && y->shift > 0) {
    return -1;
  } else if (y->shift < 0 && x->shift > 0) {
    return +1;
  } else {
    if (x->shift > y->shift) {
      return -1;
    } else if (y->shift > x->shift) {
      return +1;
    } else {
      return 0;
    }
  }
}


static int
Match_quality_cmp (const void *a, const void *b) {
  Match_T x = * (Match_T *) a;
  Match_T y = * (Match_T *) b;

  if (x->quality > y->quality) {
    return -1;
  } else if (x->quality < y->quality) {
    return +1;
  } else {
    return 0;
  }
}




typedef struct Mismatch_T *Mismatch_T;
struct Mismatch_T {
  char nt;
  char quality;
  int shift;			/* Used to record shifts */
  int count;
  Mismatch_T next;		/* Used for linking similar mismatches together */
};


static Mismatch_T
Mismatch_new (char nt, int shift, char quality) {
  Mismatch_T new = (Mismatch_T) MALLOC(sizeof(*new));

  new->nt = nt;
  new->quality = quality;
  new->shift = shift;
  new->count = 1;
  new->next = NULL;
  return new;
}

static void
Mismatch_free (Mismatch_T *old) {
  FREE(*old);
  return;
}

static int
Mismatch_chain_length (Mismatch_T this) {
  int length = 0;

  while (this != NULL) {
    length++;
    this = this->next;
  }

  return length;
}


static int
Mismatch_count_cmp (const void *a, const void *b) {
  Mismatch_T x = * (Mismatch_T *) a;
  Mismatch_T y = * (Mismatch_T *) b;

  if (x->count > y->count) {
    return -1;
  } else if (x->count < y->count) {
    return +1;
  } else {
    return 0;
  }
}


/* Go -1 to -readlength, then +readlength to +1 */
static int
Mismatch_shift_cmp (const void *a, const void *b) {
  Mismatch_T x = * (Mismatch_T *) a;
  Mismatch_T y = * (Mismatch_T *) b;

  if (x->shift < 0 && y->shift > 0) {
    return -1;
  } else if (y->shift < 0 && x->shift > 0) {
    return +1;
  } else {
    if (x->shift > y->shift) {
      return -1;
    } else if (y->shift > x->shift) {
      return +1;
    } else {
      return 0;
    }
  }
}


static int
Mismatch_quality_cmp (const void *a, const void *b) {
  Mismatch_T x = * (Mismatch_T *) a;
  Mismatch_T y = * (Mismatch_T *) b;

  if (x->quality > y->quality) {
    return -1;
  } else if (x->quality < y->quality) {
    return +1;
  } else {
    return 0;
  }
}


static Mismatch_T
find_mismatch_byshift (List_T mismatches, char nt, int shift) {
  List_T p;
  Mismatch_T mismatch;

  for (p = mismatches; p != NULL; p = List_next(p)) {
    mismatch = (Mismatch_T) List_head(p);
    if (mismatch->nt == nt && mismatch->shift == shift) {
      return mismatch;
    }
  }
  return (Mismatch_T) NULL;
}

static Mismatch_T
find_mismatch_byquality (List_T mismatches, char nt, char quality) {
  List_T p;
  Mismatch_T mismatch;

  for (p = mismatches; p != NULL; p = List_next(p)) {
    mismatch = (Mismatch_T) List_head(p);
    if (mismatch->nt == nt && mismatch->quality == quality) {
      return mismatch;
    }
  }
  return (Mismatch_T) NULL;
}

static Mismatch_T
find_mismatch_nt (List_T mismatches, char nt) {
  List_T p;
  Mismatch_T mismatch;

  for (p = mismatches; p != NULL; p = List_next(p)) {
    mismatch = (Mismatch_T) List_head(p);
    if (mismatch->nt == nt) {
      return mismatch;
    }
  }
  return (Mismatch_T) NULL;
}


#if 0
static Genomicpos_T
translate_chromosomepos (char *chromosome, unsigned int pos, IIT_T chromosome_iit) {
  int index;
  Interval_T interval;
  
  if ((index = IIT_find_linear(chromosome_iit,chromosome)) >= 0) {
    interval = IIT_interval(chromosome_iit,index);
    return Interval_low(interval)+pos;
  } else {
    fprintf(stderr,"Cannot find chromosome %s in chromosome_iit file\n",chromosome);
    abort();
  }
}
#endif


/* querystart (typically 0) and query3 (typically querylength) are exclusive */
/* sequences may have had lower case characters marked */
static int
trim_left_end (char *genomicdir, char *query, int querystart, int queryend) {
  int bestscore, score;
  int trim5, alignlength, pos;
  char *p, *q;

  alignlength = queryend - querystart;

  p = &(query[querystart]);
  q = &(genomicdir[querystart]);

  debug8(printf("querystart = %d, queryend = %d\n",querystart,queryend));
  debug8(printf("Trim left query:   %s\n",p));
  debug8(printf("Trim left genomic: %s\n",q));

  bestscore = 0;
  score = 0;
  trim5 = 0;
  for (pos = alignlength-1; pos >= 0; pos--) {
    if (toupper(p[pos]) == toupper(q[pos])) {
      score += MATCH_SCORE;
    } else {
      score += MISMATCH_SCORE;
    }
    if (score < 0) {
      score = 0;
    }
    if (score > bestscore) {
      bestscore = score;
      trim5 = pos;
    }
    debug8(printf("Trim left pos %d (%c vs %c), score %d, trim5 %d\n",
		  pos,p[pos],q[pos],score,trim5));
  }

  debug8({
      printf("Trim left called with querystart %d and queryend %d\n",querystart,queryend);
      printf("At query ->: %.*s\n",alignlength,&(query[querystart]));
      printf("At genome->: %.*s\n",alignlength,&(genomicdir[querystart]));
      printf("trim %02d  ->: ",trim5);
      for (pos = 0; pos < trim5; pos++) {
	printf(" ");
      }
      for ( ; pos < alignlength; pos++) {
	printf("*");
      }
      printf("\n");
    });

  return trim5;
}



/* querystart (typically 0) and queryend (typically querylength) are exclusive */
/* sequences may have had lower case characters marked */
static int
trim_right_end (char *genomicdir, char *query, int querystart, int queryend) {
  int bestscore, score;
  int trim3, alignlength, pos;
  char *p, *q;

  alignlength = queryend - querystart;

  p = &(query[querystart]);
  q = &(genomicdir[querystart]);

  bestscore = 0;
  score = 0;
  trim3 = 0;
  for (pos = 0; pos < alignlength; pos++) {
    if (toupper(p[pos]) == toupper(q[pos])) {
      score += MATCH_SCORE;
    } else {
      score += MISMATCH_SCORE;
    }
    if (score < 0) {
      score = 0;
    }
    if (score > bestscore) {
      bestscore = score;
      trim3 = alignlength - pos - 1;
    }
    debug8(printf("Trim right pos %d, score %d, trim3 %d\n",pos,score,trim3));
  }

  debug8({
      printf("Trim right called with querystart %d and queryend %d\n",querystart,queryend);
      printf("At query ->: %.*s\n",alignlength,&(query[querystart]));
      printf("At genome->: %.*s\n",alignlength,&(genomicdir[querystart]));
      printf("trim %02d  ->: ",trim3);
      for (pos = 0; pos < alignlength - trim3; pos++) {
	printf("*");
      }
      for ( ; pos < alignlength; pos++) {
	printf(" ");
      }
      printf("\n");
    });

  return trim3;
}


static void
compute_trims (int *trimleft, int *trimright, char *query, char *genomic,
	       int query5, int query3) {
  *trimleft = trim_left_end(genomic,query,/*query5*/0,/*query3*/query3 - query5 + 1);
  *trimright = trim_right_end(genomic,query,/*query5*/0,/*query3*/query3 - query5 + 1);
  return;
}


static List_T
lines_gc (List_T *lines) {
  char *line;
  void *item;
  
  while (*lines != NULL) {
    *lines = List_pop(*lines,&item);
    line = (char *) item;
    FREE(line);
  }
  return NULL;
}


static void
revise_data (Genomicpos_T position, char querynt, char genomicnt,
	     int quality, int signed_shift, char *refnts, int *nmatches,
	     List_T *matches_byshift, List_T *matches_byquality,
	     List_T *mismatches_byshift, List_T *mismatches_byquality,
	     int *quality_counts_match, int *quality_counts_mismatch) {
  Match_T match;
  Mismatch_T mismatch;

  if (dibasep == true) {
    if (toupper(genomicnt) == refnts[position]) {
      nmatches[position] += 1;
      if (print_refdetails_p == true) {
	if ((match = find_match_byshift(matches_byshift[position],signed_shift)) != NULL) {
	  match->count += 1;
	} else {
	  matches_byshift[position] = List_push(matches_byshift[position],(void *) Match_new(signed_shift,quality));
	}
      }

      if ((match = find_match_byquality(matches_byquality[position],quality)) != NULL) {
	match->count += 1;
      } else {
	matches_byquality[position] = List_push(matches_byquality[position],(void *) Match_new(signed_shift,quality));
      }
      quality_counts_match[quality] += 1;

    } else {
      if ((mismatch = find_mismatch_byshift(mismatches_byshift[position],toupper(genomicnt),signed_shift)) != NULL) {
	mismatch->count += 1;
      } else {
	mismatches_byshift[position] = List_push(mismatches_byshift[position],(void *) Mismatch_new(toupper(genomicnt),signed_shift,quality));
      }
      
      if ((mismatch = find_mismatch_byquality(mismatches_byquality[position],toupper(genomicnt),quality)) != NULL) {
	mismatch->count += 1;
      } else {
	mismatches_byquality[position] = List_push(mismatches_byquality[position],(void *) Mismatch_new(toupper(genomicnt),signed_shift,quality));
      }
      quality_counts_mismatch[quality] += 1;
    }

  } else if (mode == ALL_SCORES) {
    if (toupper(querynt) == toupper(genomicnt)) {
      /* Count matches, even if query quality score was low */
      nmatches[position] += 1;
      if (print_refdetails_p == true) {
	if ((match = find_match_byshift(matches_byshift[position],signed_shift)) != NULL) {
	  match->count += 1;
	} else {
	  matches_byshift[position] = List_push(matches_byshift[position],(void *) Match_new(signed_shift,quality));
	}
      }

      if ((match = find_match_byquality(matches_byquality[position],quality)) != NULL) {
	match->count += 1;
      } else {
	matches_byquality[position] = List_push(matches_byquality[position],(void *) Match_new(signed_shift,quality));
      }
      quality_counts_match[quality] += 1;

    } else {
      if ((mismatch = find_mismatch_byshift(mismatches_byshift[position],toupper(querynt),signed_shift)) != NULL) {
	mismatch->count += 1;
      } else {
	mismatches_byshift[position] = List_push(mismatches_byshift[position],(void *) Mismatch_new(toupper(querynt),signed_shift,quality));
      }

      if ((mismatch = find_mismatch_byquality(mismatches_byquality[position],toupper(querynt),quality)) != NULL) {
	mismatch->count += 1;
      } else {
	mismatches_byquality[position] = List_push(mismatches_byquality[position],(void *) Mismatch_new(toupper(querynt),signed_shift,quality));
      }
      quality_counts_mismatch[quality] += 1;
    }

  } else if (mode == HIGH_SCORES) {
    if (islower(querynt)) {
      /* Skip.  Don't count anything where query quality score was low */

    } else if (querynt == toupper(genomicnt)) {
      nmatches[position] += 1;
      if (print_refdetails_p == true) {
	if ((match = find_match_byshift(matches_byshift[position],signed_shift)) != NULL) {
	  match->count += 1;
	} else {
	  matches_byshift[position] = List_push(matches_byshift[position],(void *) Match_new(signed_shift,quality));
	}
      }

      if ((match = find_match_byquality(matches_byquality[position],quality)) != NULL) {
	match->count += 1;
      } else {
	matches_byquality[position] = List_push(matches_byquality[position],(void *) Match_new(signed_shift,quality));
      }
      quality_counts_match[quality] += 1;

    } else {
      if ((mismatch = find_mismatch_byshift(mismatches_byshift[position],querynt,signed_shift)) != NULL) {
	mismatch->count += 1;
      } else {
	mismatches_byshift[position] = List_push(mismatches_byshift[position],(void *) Mismatch_new(querynt,signed_shift,quality));
      }

      if ((mismatch = find_mismatch_byquality(mismatches_byquality[position],querynt,quality)) != NULL) {
	mismatch->count += 1;
      } else {
	mismatches_byquality[position] = List_push(mismatches_byquality[position],(void *) Mismatch_new(querynt,signed_shift,quality));
      }
      quality_counts_mismatch[quality] += 1;
    }

  } else if (mode == MISMATCH_HIGH_SCORES) {
    /* Counts all matches, regardless of quality score, but only mismatches with high quality score */
    if (toupper(querynt) == toupper(genomicnt)) {
      /* Count matches, even if query quality score was low */
      nmatches[position] += 1;
      if (print_refdetails_p == true) {
	if ((match = find_match_byshift(matches_byshift[position],signed_shift)) != NULL) {
	  match->count += 1;
	} else {
	  matches_byshift[position] = List_push(matches_byshift[position],(void *) Match_new(signed_shift,quality));
	}
      }

      if ((match = find_match_byquality(matches_byquality[position],quality)) != NULL) {
	match->count += 1;
      } else {
	matches_byquality[position] = List_push(matches_byquality[position],(void *) Match_new(signed_shift,quality));
      }
      quality_counts_match[quality] += 1;

    } else if (islower(querynt)) {
      /* Skip.  Don't count a mismatch if query quality score was low */

    } else {
      if ((mismatch = find_mismatch_byshift(mismatches_byshift[position],querynt,signed_shift)) != NULL) {
	mismatch->count += 1;
      } else {
	mismatches_byshift[position] = List_push(mismatches_byshift[position],(void *) Mismatch_new(querynt,signed_shift,quality));
      }

      if ((mismatch = find_mismatch_byquality(mismatches_byquality[position],querynt,quality)) != NULL) {
	mismatch->count += 1;
      } else {
	mismatches_byquality[position] = List_push(mismatches_byquality[position],(void *) Mismatch_new(querynt,signed_shift,quality));
      }
      quality_counts_mismatch[quality] += 1;
    }

  }

  return;
}



#ifdef SAM_INPUT

static void
process_lines (char *refnts, int *nmatches, List_T *matches_byshift, List_T *matches_byquality,
	       List_T *mismatches_byshift, List_T *mismatches_byquality,
	       int *quality_counts_match, int *quality_counts_mismatch,
	       char *chromosome, Pairend_T pairend, List_T lines,
	       IIT_T chromosome_iit, Genome_T genome, Genomicpos_T chroffset) {
  List_T ptr, ptr2;
  char *line, *p, *q, *r;
  char genomic[1024], *acc, *chr, *cigar, *shortread, *quality_string, strand;
  int shift, signed_shift, mapq;
  Chrnum_T chrnum;
  Genomicpos_T pos;
  int trimleft = 0, trimright = 0, ntrim;
  int querylength;
  int nunknowns;
  int index;
  Genomicpos_T allocoffset = 0U, chrpos;
  bool processp;

  unsigned int flag, mlength;
  int type;
  int cigar_querylength;
  Intlist_T types, a;
  Uintlist_T npositions, b;
  
  for (ptr = lines; ptr != NULL; ptr = List_next(ptr)) {
    line = (char *) List_head(ptr);
    Samread_parse_line(&acc,&flag,&mapq,&chr,&chrpos,&cigar,&querylength,
		       &shortread,&quality_string,line);

    if (flag & QUERY_MINUSP) {
      strand = '-';
    } else {
      strand = '+';
    }

#if 0
    /* Already computed by calling procedure */
    if (flag & FIRST_READ_P) {
      pairend = FIRST;
    } else if (flag & SECOND_READ_P) {
      pairend = SECOND;
    } else {
      pairend = FIRST;
    }
#endif

    if (mapq < minimum_mapq) {
      processp = false;

    } else if (desired_strand != ' ' && ((pairend == FIRST && strand != desired_strand) ||
				  (pairend == SECOND && strand != opposite_strand))) {
      processp = false;

    } else if (chromosome == NULL) {
      if ((index = IIT_find_one(chromosome_iit,chr)) < 0) {
	processp = false;
      } else {
	allocoffset = chroffset = Interval_low(IIT_interval(chromosome_iit,index));
	processp = true;
      }

    } else if (strcmp(chr,chromosome)) {
      /* Different chromosome */
      debug(printf("Ignoring different chromosome %s.  Wanted %s\n",chr,chromosome));
      processp = false;
      
    } else {
      processp = true;
    }

    if (processp == true) {
      types = Samread_parse_cigar(&npositions,&cigar_querylength,cigar);

#if 0
      /* Doesn't hold for hard clipping */
      if (cigar_querylength != querylength) {
	fprintf(stderr,"Cigar querylength = %d, but read has length %d\n",cigar_querylength,querylength);
	exit(9);
      }
#endif

      /* Get query coordinates */
      /* Samread_get_query_coordinates(&query5,&query3,types,npositions,querylength,cigar); */

      pos = chrpos-1U;		/* SAM chrpos coordinates are 1-based */
      if (strand == '+') {
	shift = 1;
      } else {
	shift = cigar_querylength;
      }
      p = shortread;
      r = quality_string;
      for (a = types, b = npositions; a != NULL; a = Intlist_next(a), b = Uintlist_next(b)) {
	type = Intlist_head(a);
	mlength = Uintlist_head(b);
	if (type == 'S') {
	  /* pos += mlength; -- SAM assumes genome coordinates are of clipped region */
	  p += mlength;
	  r += mlength;
	  shift += ((strand == '+') ? +mlength : -mlength);
	} else if (type == 'H') {
	  /* pos += mlength; -- SAM assumes genome coordinates are of clipped region */
	  /* p += mlength; -- hard clip means query sequence is absent */
	  /* r += mlength; -- hard clip means quality string is absent */
	  shift += ((strand == '+') ? +mlength : -mlength);
	} else if (type == 'N') {
	  pos += mlength;
	} else if (type == 'I') {
	  p += mlength;
	  r += mlength;
	  shift += ((strand == '+') ? +mlength : -mlength);
	} else if (type == 'D') {
	  pos += mlength;
	} else if (type == 'M') {
	  if (mlength < min_mlength) {
	    debug(printf("mlength %d < min_mlength %d\n",mlength,min_mlength));
	    p += mlength;
	    r += mlength;
	    pos += mlength;
	    shift += ((strand == '+') ? +mlength : -mlength);

	  } else {
	    debug(printf("Genomic pos is %u + %u = %u\n",chroffset,pos,chroffset+pos));
	    Genome_fill_buffer(&chrnum,&nunknowns,genome,chroffset + pos,mlength,genomic,chromosome_iit);
	    /* printf("After (+): %s\n",genomic); */
	    q = genomic;

	    /* psave = p; qsave = q; */
	    debug(printf("Processing %.*s and %.*s\n",mlength,p,mlength,q));

#if 0
	    if (trimp == true) {
	      compute_trims(&trimleft,&trimright,p,q,mlength,strand,/*query5*/1,/*query3*/mlength,querylength);
	      debug(printf("Computed trimleft of %d and trimright of %d\n",trimleft,trimright));
	    }
#else
	    trimleft = trimright = 0;
#endif

	    ntrim = 0;
	    while (ntrim++ < trimleft) {
	      debug(printf("Trim left: Skipping %c and %c at shift %d, pos %u\n",*p,*q,shift,pos+1U));
	      p++;
	      q++;
	      r++;
	      pos++;
	      shift += ((strand == '+') ? +1 : -1);
	      mlength--;
	    }

	    while (mlength-- > trimright) {
	      debug(printf("Processing %c and %c at shift %d, pos %u+%u, mlength %u\n",
			   *p,*q,shift,allocoffset,pos+1U,mlength));
	      if (nmatches[allocoffset+pos] == 0 && mismatches_byshift[allocoffset+pos] == NULL) {
		refnts[allocoffset+pos] = toupper(*q);

	      } else if (allocoffset+pos == haltpos) {
		fprintf(stderr,"Halt requested at %u+%u:\n",allocoffset,pos);
		fprintf(stderr,"Have seen %d matches and %d types of mismatches here so far\n",
			nmatches[allocoffset+pos],List_length(mismatches_byshift[allocoffset+pos]));
		lines = List_reverse(lines);
		for (ptr2 = lines; ptr2 != NULL; ptr2 = List_next(ptr2)) {
		  fprintf(stderr,"%s",(char *) List_head(ptr2));
		}
		exit(9);

	      } else if (refnts[allocoffset+pos] != toupper(*q)) {
		fprintf(stderr,"Two different genomic chars %c and %c at position %u+%u\n",
			refnts[allocoffset+pos],*q,allocoffset,pos+1U);
		fprintf(stderr,"Have seen %d matches and %d types of mismatches here so far\n",
			nmatches[allocoffset+pos],List_length(mismatches_byshift[allocoffset+pos]));
		fprintf(stderr,"Most recent lines were:\n");
		lines = List_reverse(lines);
		for (ptr2 = lines; ptr2 != NULL; ptr2 = List_next(ptr2)) {
		  fprintf(stderr,"%s",(char *) List_head(ptr2));
		}
		fprintf(stderr,"To halt at first occurrence, re-run with -h <coord>\n");
		exit(9);

	      }

	      signed_shift = (strand == '+') ? shift : -shift;
	      revise_data(/*position*/allocoffset+pos,/*querynt*/*p,/*genomicnt*/*q,/*quality*/(int) *r,
			  signed_shift,refnts,nmatches,matches_byshift,matches_byquality,
			  mismatches_byshift,mismatches_byquality,
			  quality_counts_match,quality_counts_mismatch);
	      
	      p++;
	      q++;
	      r++;
	      pos++;
	      shift += ((strand == '+') ? +1 : -1);
	    }
	  }

	} else {
	  fprintf(stderr,"Cannot parse type '%c'\n",type);
	  exit(9);
	}
      }

      Intlist_free(&types);
      Uintlist_free(&npositions);
    }

    FREE(acc);
    FREE(chr);
    FREE(cigar);
    FREE(shortread);
    FREE(quality_string);
  }

  return;
}

static void
parse_sam (char *refnts, int *nmatches, List_T *matches_byshift, List_T *matches_byquality,
	   List_T *mismatches_byshift, List_T *mismatches_byquality,
	   int *quality_counts_match, int *quality_counts_mismatch,
	   char *chromosome, IIT_T chromosome_iit, Genome_T genome, Genomicpos_T chroffset) {
  char line[1024000], *copy;
  char *acc, *lastacc;
  unsigned int flag;
  int nhits_first = 0, nhits_second = 0;
  List_T lines_first = NULL, lines_second = NULL;
  bool concordantp;

  lastacc = (char *) CALLOC(1,sizeof(char));
  lastacc[0] = '\0';

  while (fgets(line,1024000,stdin) != NULL) {
    if (line[0] == '@') {
      /* Skip */
    } else {
      acc = Samread_get_acc(&flag,line);

      if (strcmp(acc,lastacc)) {
	if (lastacc[0] != '\0') {
	  if ((need_concordant_p == false || concordantp == true) &&
	      (uniquep == false || nhits_first == 1)) {
	    process_lines(refnts,nmatches,matches_byshift,matches_byquality,
			  mismatches_byshift,mismatches_byquality,
			  quality_counts_match,quality_counts_mismatch,chromosome,
			  /*pairend*/FIRST,lines_first,chromosome_iit,genome,chroffset);
	  }
	  if ((need_concordant_p == false || concordantp == true) && 
	      (uniquep == false || nhits_second == 1)) {
	    process_lines(refnts,nmatches,matches_byshift,matches_byquality,
			  mismatches_byshift,mismatches_byquality,
			  quality_counts_match,quality_counts_mismatch,chromosome,
			  /*pairend*/SECOND,lines_second,chromosome_iit,genome,chroffset);
	  }
	  lines_first = lines_gc(&lines_first);
	  lines_second = lines_gc(&lines_second);
	  nhits_first = 0;
	  nhits_second = 0;
	}
	FREE(lastacc);
	lastacc = acc;
      }

      if (flag & QUERY_UNMAPPED) {
	/* Skip line */
      } else {
	copy = (char *) CALLOC(strlen(line)+1,sizeof(char));
	strcpy(copy,line);
	if (!(flag & PAIRED_READ)) {
	  concordantp = false;
	  lines_first = List_push(lines_first,(void *) copy);
	  nhits_first++;
	} else if (flag & FIRST_READ_P) {
	  concordantp = (flag & PAIRED_MAPPING) ? true : false;
	  lines_first = List_push(lines_first,(void *) copy);
	  nhits_first++;
	} else if (flag & SECOND_READ_P) {
	  concordantp = (flag & PAIRED_MAPPING) ? true : false;
	  lines_second = List_push(lines_second,(void *) copy);
	  nhits_second++;
	} else {
	  fprintf(stderr,"Flag %u is paired (%u), but contains neither first_read nor second_read flag\n",
		  flag,flag & PAIRED_READ);
	  abort();
	}
      }
    }
  }

  if (lastacc[0] != '\0') {
    if ((need_concordant_p == false || concordantp == true) &&
	(uniquep == false || nhits_first == 1)) {
      process_lines(refnts,nmatches,matches_byshift,matches_byquality,
		    mismatches_byshift,mismatches_byquality,
		    quality_counts_match,quality_counts_mismatch,chromosome,
		    /*pairend*/FIRST,lines_first,chromosome_iit,genome,chroffset);
    }
    if ((need_concordant_p == false || concordantp == true) &&
	(uniquep == false || nhits_second == 1)) {
      process_lines(refnts,nmatches,matches_byshift,matches_byquality,
		    mismatches_byshift,mismatches_byquality,
		    quality_counts_match,quality_counts_mismatch,chromosome,
		    /*pairend*/SECOND,lines_second,chromosome_iit,genome,chroffset);
    }
    lines_gc(&lines_first);
    lines_gc(&lines_second);
  }
  FREE(lastacc);

  return;
}


#elif defined(BAM_INPUT)

static void
parse_bam (char *refnts, int *nmatches, List_T *matches_byshift, List_T *matches_byquality,
	   List_T *mismatches_byshift, List_T *mismatches_byquality,
	   int *quality_counts_match, int *quality_counts_mismatch,
	   Bamreader_T bamreader, char *chromosome, IIT_T chromosome_iit,
	   Genome_T genome, Genomicpos_T chroffset) {
  char *acc, *chr, *shortread, *quality_string, strand;
  unsigned int flag, mlength;
  Genomicpos_T allocoffset = 0U, chrpos;
  Intlist_T types, a;
  Uintlist_T npositions, b;
  Pairend_T pairend;
  int readlength, cigar_querylength, shift, signed_shift, mapq;
  int index;
  int type;
  bool processp;

  char genomic[1024], *p, *q, *r;
  Chrnum_T chrnum;
  Genomicpos_T pos;
  int nunknowns;
  int trimleft = 0, trimright = 0, ntrim;


  while (Bamread_next_line(bamreader,&acc,&flag,&mapq,&chr,&chrpos,&types,&npositions,&cigar_querylength,
			   &readlength,&shortread,&quality_string) > 0) {
    if (mapq < minimum_mapq) {
      processp = false;
    } else if (uniquep == true && (flag & NOT_PRIMARY)) {
      processp = false;
    } else {
      if (flag & QUERY_MINUSP) {
	strand = '-';
      } else {
	strand = '+';
      }

      if (flag & FIRST_READ_P) {
	pairend = FIRST;
      } else if (flag & SECOND_READ_P) {
	pairend = SECOND;
      } else {
	fprintf(stderr,"Flag %u for acc %s has neither first read nor second read set\n",flag,acc);
	pairend = FIRST;
      }

      if (desired_strand != ' ' && ((pairend == FIRST && strand != desired_strand) ||
				    (pairend == SECOND && strand != opposite_strand))) {
	processp = false;

      } else if (chr == NULL) {
	/* No alignment */
	processp = false;

      } else if (chromosome == NULL) {
	if ((index = IIT_find_one(chromosome_iit,chr)) < 0) {
	  processp = false;
	} else {
	  allocoffset = chroffset = Interval_low(IIT_interval(chromosome_iit,index));
	  processp = true;
	}

      } else if (strcmp(chr,chromosome)) {
	/* Different chromosome */
	debug(printf("Ignoring different chromosome %s.  Wanted %s\n",chr,chromosome));
	processp = false;
      
      } else {
	processp = true;
      }
    }

    if (processp == true) {
      pos = chrpos-1U;		/* Bamread reports chrpos as 1-based */
      if (strand == '+') {
	shift = 1;
      } else {
	shift = cigar_querylength;
      }
      p = shortread;
      r = quality_string;
      for (a = types, b = npositions; a != NULL; a = Intlist_next(a), b = Uintlist_next(b)) {
	type = Intlist_head(a);
	mlength = Uintlist_head(b);
	if (type == 'S') {
	  /* pos += mlength; -- SAM assumes genome coordinates are of clipped region */
	  p += mlength;
	  r += mlength;
	  shift += ((strand == '+') ? +mlength : -mlength);
	} else if (type == 'H') {
	  /* pos += mlength; -- SAM assumes genome coordinates are of clipped region */
	  /* p += mlength; -- hard clip means query sequence is absent */
	  /* r += mlength; -- hard clip means quality string is absent */
	  shift += ((strand == '+') ? +mlength : -mlength);
	} else if (type == 'N') {
	  pos += mlength;
	} else if (type == 'I') {
	  p += mlength;
	  r += mlength;
	  shift += ((strand == '+') ? +mlength : -mlength);
	} else if (type == 'D') {
	  pos += mlength;
	} else if (type == 'M') {
	  if (mlength < min_mlength) {
	    debug(printf("mlength %d < min_mlength %d\n",mlength,min_mlength));
	    p += mlength;
	    r += mlength;
	    pos += mlength;
	    shift += ((strand == '+') ? +mlength : -mlength);

	  } else {
	    debug(printf("Genomic pos is %u + %u = %u\n",chroffset,pos,chroffset+pos));
	    Genome_fill_buffer(&chrnum,&nunknowns,genome,chroffset + pos,mlength,genomic,chromosome_iit);
	    /* printf("After (+): %s\n",genomic); */
	    q = genomic;

	    /* psave = p; qsave = q; */
	    debug(printf("Processing %.*s and %.*s\n",mlength,p,mlength,q));
#if 0
	    if (trimp == true) {
	      compute_trims(&trimleft,&trimright,p,q,mlength,strand,/*query5*/1,/*query3*/mlength,querylength);
	      debug(printf("Computed trimleft of %d and trimright of %d\n",trimleft,trimright));
	    }
#else
	    trimleft = trimright = 0;
#endif

	    ntrim = 0;
	    while (ntrim++ < trimleft) {
	      debug(printf("Trim left: Skipping %c and %c at shift %d, pos %u\n",*p,*q,shift,pos+1U));
	      p++;
	      q++;
	      r++;
	      pos++;
	      shift += ((strand == '+') ? +1 : -1);
	      mlength--;
	    }

	    while (mlength-- > trimright) {
	      debug(printf("Processing %c and %c at shift %d, pos %u+%u, mlength %u\n",
			   *p,*q,shift,allocoffset,pos+1U,mlength));
	      if (nmatches[allocoffset+pos] == 0 && mismatches_byshift[allocoffset+pos] == NULL) {
		refnts[allocoffset+pos] = toupper(*q);

	      } else if (refnts[allocoffset+pos] != toupper(*q)) {
		fprintf(stderr,"Two different genomic chars %c and %c at position %u+%u\n",
			refnts[allocoffset+pos],*q,allocoffset,pos+1U);
		fprintf(stderr,"Have seen %d matches and %d types of mismatches here so far\n",
			nmatches[allocoffset+pos],List_length(mismatches_byshift[allocoffset+pos]));
		exit(9);

	      }

	      signed_shift = (strand == '+') ? shift : -shift;
	      revise_data(/*position*/allocoffset+pos,/*querynt*/*p,/*genomicnt*/*q,/*quality*/(int) *r,
			  signed_shift,refnts,nmatches,matches_byshift,matches_byquality,
			  mismatches_byshift,mismatches_byquality,
			  quality_counts_match,quality_counts_mismatch);
	      
	      p++;
	      q++;
	      r++;
	      pos++;
	      shift += ((strand == '+') ? +1 : -1);
	    }
	  }

	} else {
	  fprintf(stderr,"Cannot parse type '%c'\n",type);
	  exit(9);
	}
      }

      Intlist_free(&types);
      Uintlist_free(&npositions);
    }

    FREE(shortread);
  }

  return;
}



#else

static char complCode[128] = COMPLEMENT_LC;

static void
make_complement_buffered (char *complement, char *sequence, int length) {
  int i, j;

  /* complement = (char *) CALLOC(length+1,sizeof(char)); */
  for (i = length-1, j = 0; i >= 0; i--, j++) {
    complement[j] = complCode[(int) sequence[i]];
  }
  complement[length] = '\0';
  return;
}

static void
make_complement_inplace (char *sequence, int length) {
  char temp;
  int i, j;

  for (i = 0, j = length-1; i < length/2; i++, j--) {
    temp = complCode[(int) sequence[i]];
    sequence[i] = complCode[(int) sequence[j]];
    sequence[j] = temp;
  }
  if (i == j) {
    sequence[i] = complCode[(int) sequence[i]];
  }

  return;
}


static void
make_reverse (char *reverse, char *sequence, int length) {
  int i, j;

  for (i = length-1, j = 0; i >= 0; i--, j++) {
    reverse[j] = sequence[i];
  }
  reverse[length] = '\0';
  return;
}


static void
process_line (char *line, char *header, char *refnts, int *nmatches,
	      List_T *matches_byshift, List_T *matches_byquality,
	      List_T *mismatches_byshift, List_T *mismatches_byquality,
	      int *quality_counts_match, int *quality_counts_mismatch,
	      char *chromosome, IIT_T chromosome_iit, Genome_T genome, Genomicpos_T chroffset) {
  char *p, *q, *r;
  char *genomic, *chr, strand, truestrand, firstchar, secondchar;
  char *shortread, *quality_string, shortread_rc[1024], genomic_rc[1024], quality_rev[1024];
  int query5, query3, readlength, support, nmismatches, querypos, signed_querypos;
  Chrnum_T chrnum;
  Pairend_T pairend;
  Genomicpos_T pos5, pos3, pos, donorpos, acceptorpos;
  int trimleft = 0, trimright = 0, ntrim;
  int nunknowns;
  int index;
  Genomicpos_T allocoffset = 0U;
  bool sensep;

  shortread = &(header[1]);
  quality_string = Gsnapread_quality_string(header);

  if (header[0] == '>') {
    pairend = FIRST;
  } else if (header[0] == '<') {
    pairend = SECOND;
  } else {
    fprintf(stderr,"Bad header: %s\n",header);
    abort();
  }

  genomic = &(line[1]);

  Gsnapread_parse_line(line,&query5,&query3,&firstchar,&secondchar,
		       &support,&nmismatches,&strand,&chr,&pos5,&pos3,
		       &donorpos,&acceptorpos,&truestrand,&sensep);

  if (support < min_readlength) {
    debug(printf("readlength %d < min_readlength %d\n",support,min_readlength));
    readlength = 0;
  } else if (desired_strand != ' ' && ((pairend == FIRST && strand != desired_strand) ||
				       (pairend == SECOND && strand != opposite_strand))) {
    readlength = 0;
  } else if (chromosome == NULL) {
    if ((index = IIT_find_one(chromosome_iit,chr)) < 0) {
      readlength = 0;
    } else {
      allocoffset = chroffset = Interval_low(IIT_interval(chromosome_iit,index));
      readlength = support;
    }
    
  } else if (strcmp(chr,chromosome)) {
    /* Different chromosome */
    debug(printf("Ignoring different chromosome %s.  Wanted %s\n",chr,chromosome));
    readlength = 0;
    
  } else {
    readlength = support;
  }
  
  FREE(chr);

  if (readlength > 0) {
    if (strand == '+') {
      pos = pos5-1U;		/* Based pos3 and pos5 were 1-based */
      p = &(shortread[query5-1]);
      if (quality_string == NULL) {
	r = (char *) NULL;
      } else {
	r = &(quality_string[query5-1]);
      }
      if (genome == NULL) {
	/* Use genome according to GSNAP */
	q = &(genomic[query5-1]);
      } else {
	debug(printf("Looking up genomic segment\n"));
	/* Look up genomic segment.  Don't use known gene. */
	/* printf("Before (+): %s\n",&(genomic[query5-1])); */
	Genome_fill_buffer(&chrnum,&nunknowns,genome,chroffset + pos /*=pos5-1U*/,readlength,genomic,chromosome_iit);
	/* printf("After (+): %s\n",genomic); */
	q = genomic;
      }
    } else if (strand == '-') {
      pos = pos3 - 1U;
      make_complement_buffered(shortread_rc,&(shortread[query5-1]),readlength);
      p = shortread_rc;
      if (quality_string == NULL) {
	r = (char *) NULL;
      } else {
	make_reverse(quality_rev,&(quality_string[query5-1]),readlength);
	r = quality_rev;
      }
      if (genome == NULL) {
	/* Use genome according to GSNAP */
	make_complement_buffered(genomic_rc,&(genomic[query5-1]),readlength);
	q = genomic_rc;
      } else {
	debug(printf("Looking up genomic segment\n"));
	/* Look up genomic segment.  Don't use known gene. */
	/* printf("Before (-): %s\n",&(genomic[query5-1])); */
	Genome_fill_buffer(&chrnum,&nunknowns,genome,chroffset + pos /*=pos3-1U*/,readlength,genomic_rc,chromosome_iit);
	/* printf("After (-): %s\n",genomic_rc); */
	q = genomic_rc;
      }
    } else {
      fprintf(stderr,"Can't parse strand %c in %s\n",strand,line);
      abort();
    }

    debug(printf("Processing strand %c, %.*s and %.*s\n",strand,readlength,p,readlength,q));
    debug(printf("shortread: %s\n",shortread));
    debug(printf("shortread_rc: %s\n",shortread_rc));

#if 0
    if (trimp == true) {
      compute_trims(&trimleft,&trimright,p,q,query5,query3);
      debug(printf("Computed trimleft of %d and trimright of %d\n",trimleft,trimright));
    }
#endif

    if (strand == '+') {
      querypos = query5;
    } else {
      querypos = query3;
    }

    ntrim = 0;
    while (ntrim++ < trimleft) {
      debug(printf("Trim left: Skipping %c and %c at querypos %d, pos %u\n",*p,*q,querypos,pos+1U));
      p++;
      q++;
      if (r) r++;
      pos++;
      querypos += ((strand == '+') ? +1 : -1);
      readlength--;
    }

    while (--readlength >= trimright) {
      debug(printf("Processing %c and %c at querypos %d, pos %u+%u\n",*p,*q,querypos,allocoffset,pos+1U));
      if (nmatches[allocoffset+pos] == 0 && mismatches_byshift[allocoffset+pos] == NULL) {
	refnts[allocoffset+pos] = toupper(*q);

      } else if (allocoffset+pos == haltpos) {
	fprintf(stderr,"Halt requested at %u+%u:\n",allocoffset,pos);
	fprintf(stderr,"Have seen %d matches and %d types of mismatches here so far\n",
		nmatches[allocoffset+pos],List_length(mismatches_byshift[allocoffset+pos]));
	fprintf(stderr,"%s",header);
	fprintf(stderr,"%s",line);
	exit(9);

      } else if (refnts[allocoffset+pos] != toupper(*q)) {
	fprintf(stderr,"Two different genomic chars %c and %c at position %u+%u\n",
		refnts[allocoffset+pos],*q,allocoffset,pos+1U);
	fprintf(stderr,"Have seen %d matches and %d types of mismatches here so far\n",
		nmatches[allocoffset+pos],List_length(mismatches_byshift[allocoffset+pos]));
	fprintf(stderr,"Most recent lines were:\n");
	fprintf(stderr,"%s",header);
	fprintf(stderr,"%s",line);
	exit(9);

      }

      signed_querypos = (strand == '+') ? querypos : -querypos;
      revise_data(/*position*/allocoffset+pos,/*querynt*/*p,/*genomicnt*/*q,
		  /*quality*/r ? (int) *r : DEFAULT_QUALITY,
		  signed_querypos,refnts,nmatches,matches_byshift,matches_byquality,
		  mismatches_byshift,mismatches_byquality,
		  quality_counts_match,quality_counts_mismatch);

      p++;
      q++;
      if (r) r++;
      pos++;
      querypos += ((strand == '+') ? +1 : -1);
    }

    debug({
	while (readlength-- >= 0) {
	  printf("Trim right: Skipping %c and %c at querypos %d, pos %u\n",*p,*q,querypos,pos+1U);
	  p++;
	  q++;
	  if (r) r++;
	  pos++;
	  querypos += ((strand == '+') ? +1 : -1);
	}
      });

  }

  if (quality_string != NULL) {
    FREE(quality_string);
  }

  return;
}


static void
parse_single (char *refnts, int *nmatches, List_T *matches_byshift, List_T *matches_byquality,
	      List_T *mismatches_byshift, List_T *mismatches_byquality,
	      int *quality_counts_match, int *quality_counts_mismatch,
	      char *chromosome, IIT_T chromosome_iit, Genome_T genome, Genomicpos_T chroffset) {
  char line[1024000], *copy, *header = NULL;
  List_T lines = NULL, p;

  while (fgets(line,1024000,stdin) != NULL) {
    if (line[0] == '\n') {
      /* Skip blank line */

    } else if (line[0] == '>' || line[0] == '<') {

      if (header != NULL && 
	  (uniquep == false || Gsnapread_nhits(header) == 1)) {
	/* Process previous lines */
	p = lines = List_reverse(lines);
	while (p != NULL) {
	  process_line((char *) List_head(p),header,refnts,nmatches,
		       matches_byshift,matches_byquality,
		       mismatches_byshift,mismatches_byquality,
		       quality_counts_match,quality_counts_mismatch,
		       chromosome,chromosome_iit,genome,chroffset);
	  p = List_next(p);
	}
      }

      lines = lines_gc(&lines);
      FREE(header);

      /* Process this header */
      header = (char *) CALLOC(strlen(line)+1,sizeof(char));
      strcpy(header,line);

    } else {
      copy = (char *) CALLOC(strlen(line)+1,sizeof(char));
      strcpy(copy,line);
      debug(printf("Pushing %s",copy));
      lines = List_push(lines,(void *) copy);

    }

  }

  if (header != NULL && 
      (uniquep == false || Gsnapread_nhits(header) == 1)) {
    /* Process previous lines */
    p = lines = List_reverse(lines);
    while (p != NULL) {
      process_line((char *) List_head(p),header,refnts,nmatches,
		   matches_byshift,matches_byquality,
		   mismatches_byshift,mismatches_byquality,
		   quality_counts_match,quality_counts_mismatch,
		   chromosome,chromosome_iit,genome,chroffset);
      p = List_next(p);
    }
  }

  lines_gc(&lines);
  FREE(header);

  return;
}


static List_T
advance_one_hit (List_T p) {
  char *line;

  p = List_next(p);		/* Skip past first line of current hit */

  while (p != NULL) {
    line = (char *) List_head(p);
    if (line[0] == '\0') {
      return (List_T) NULL;
    } else if (line[0] == ' ') {
      return p;
    } else if (line[0] == ',') {
      p = List_next(p);
    } else {
      fprintf(stderr,"advance_one_hit: Don't recognize %s\n",line);
    }
  }

  return (List_T) NULL;
}


static void
parse_paired (char *refnts, int *nmatches, List_T *matches_byshift, List_T *matches_byquality,
	      List_T *mismatches_byshift, List_T *mismatches_byquality,
	      int *quality_counts_match, int *quality_counts_mismatch,
	      char *chromosome, IIT_T chromosome_iit, Genome_T genome, Genomicpos_T chroffset) {
  char line[1024000], *copy, *header1 = NULL, *header2 = NULL, *chr1, *chr2;
  List_T lines1 = NULL, lines2 = NULL, p, q, pp, qq;
  int query5, query3, support, nmismatches;
  char strand, truestrand, firstchar, secondchar;
  Genomicpos_T pos5_1, pos5_2, pos3, donorpos, acceptorpos;
  bool saw_pair_p = false, sensep;

  while (fgets(line,1024000,stdin) != NULL) {
    if (line[0] == '\n') {
      /* Skip blank line */

    } else if (line[0] == '>') {

      if (header1 != NULL && Gsnapread_concordantp(header1) == true &&
	  (uniquep == false || Gsnapread_nhits(header1) == 1)) {
	/* Process previous lines */
	p = lines1 = List_reverse(lines1);
	q = lines2 = List_reverse(lines2);
	while (p != NULL && q != NULL) {
	  Gsnapread_parse_line((char *) List_head(p),&query5,&query3,&firstchar,&secondchar,
			       &support,&nmismatches,&strand,
			       &chr1,&pos5_1,&pos3,&donorpos,&acceptorpos,&truestrand,&sensep);
	  Gsnapread_parse_line((char *) List_head(q),&query5,&query3,&firstchar,&secondchar,
			       &support,&nmismatches,&strand,
			       &chr2,&pos5_2,&pos3,&donorpos,&acceptorpos,&truestrand,&sensep);

	  if (strcmp(chr1,chr2) && allow_translocations_p == false) {
	      /* Skip */

	  } else if (pos5_1 < pos5_2) {
	    if (want_low_p == true) {
	      for (pp = p; pp != advance_one_hit(p); pp = List_next(pp)) {
		process_line((char *) List_head(pp),header1,refnts,nmatches,
			     matches_byshift,matches_byquality,
			     mismatches_byshift,mismatches_byquality,
			     quality_counts_match,quality_counts_mismatch,
			     chromosome,chromosome_iit,genome,chroffset);
	      }
	    }
	    if (want_high_p == true) {
	      for (qq = q; qq != advance_one_hit(q); qq = List_next(qq)) {
		process_line((char *) List_head(qq),header2,refnts,nmatches,
			     matches_byshift,matches_byquality,
			     mismatches_byshift,mismatches_byquality,
			     quality_counts_match,quality_counts_mismatch,
			     chromosome,chromosome_iit,genome,chroffset);
	      }
	    }

	  } else {
	    if (want_low_p == true) {
	      for (qq = q; qq != advance_one_hit(q); qq = List_next(qq)) {
		process_line((char *) List_head(qq),header2,refnts,nmatches,
			     matches_byshift,matches_byquality,
			     mismatches_byshift,mismatches_byquality,
			     quality_counts_match,quality_counts_mismatch,
			     chromosome,chromosome_iit,genome,chroffset);
	      }
	    }
	    if (want_high_p == true) {
	      for (pp = p; pp != advance_one_hit(p); pp = List_next(pp)) {
		process_line((char *) List_head(pp),header1,refnts,nmatches,
			     matches_byshift,matches_byquality,
			     mismatches_byshift,mismatches_byquality,
			     quality_counts_match,quality_counts_mismatch,
			     chromosome,chromosome_iit,genome,chroffset);
	      }
	    }
	  }

	  FREE(chr1);
	  FREE(chr2);

	  p = advance_one_hit(p);
	  q = advance_one_hit(q);
	}
      }
      lines_gc(&lines2);
      lines_gc(&lines1);
      FREE(header2);
      FREE(header1);
      lines1 = lines2 = (List_T) NULL;

      /* Process this header */
      header1 = (char *) CALLOC(strlen(line)+1,sizeof(char));
      strcpy(header1,line);
      saw_pair_p = false;

    } else if (line[0] == '<') {
      header2 = (char *) CALLOC(strlen(line)+1,sizeof(char));
      strcpy(header2,line);
      saw_pair_p = true;
      
    } else if (saw_pair_p == false) {
      copy = (char *) CALLOC(strlen(line)+1,sizeof(char));
      strcpy(copy,line);
      debug(printf("Pushing %s",copy));
      lines1 = List_push(lines1,(void *) copy);

    } else {
      copy = (char *) CALLOC(strlen(line)+1,sizeof(char));
      strcpy(copy,line);
      debug(printf("Pushing %s",copy));
      lines2 = List_push(lines2,(void *) copy);

    }

  }

  if (header1 != NULL && Gsnapread_concordantp(header1) == true &&
      (uniquep == false || Gsnapread_nhits(header1) == 1)) {
    /* Process previous lines */
    p = lines1 = List_reverse(lines1);
    q = lines2 = List_reverse(lines2);
    while (p != NULL && q != NULL) {
      Gsnapread_parse_line((char *) List_head(p),&query5,&query3,&firstchar,&secondchar,
			   &support,&nmismatches,&strand,
			   &chr1,&pos5_1,&pos3,&donorpos,&acceptorpos,&truestrand,&sensep);
      Gsnapread_parse_line((char *) List_head(q),&query5,&query3,&firstchar,&secondchar,
			   &support,&nmismatches,&strand,
			   &chr2,&pos5_2,&pos3,&donorpos,&acceptorpos,&truestrand,&sensep);

      if (strcmp(chr1,chr2) && allow_translocations_p == false) {
	/* Skip */

      } else if (pos5_1 < pos5_2) {
	if (want_low_p == true) {
	  for (pp = p; pp != advance_one_hit(p); pp = List_next(pp)) {
	    process_line((char *) List_head(pp),header1,refnts,nmatches,
			 matches_byshift,matches_byquality,
			 mismatches_byshift,mismatches_byquality,
			 quality_counts_match,quality_counts_mismatch,
			 chromosome,chromosome_iit,genome,chroffset);
	  }
	}
	if (want_high_p == true) {
	  for (qq = q; qq != advance_one_hit(q); qq = List_next(qq)) {
	    process_line((char *) List_head(qq),header2,refnts,nmatches,
			 matches_byshift,matches_byquality,
			 mismatches_byshift,mismatches_byquality,
			 quality_counts_match,quality_counts_mismatch,
			 chromosome,chromosome_iit,genome,chroffset);
	  }
	}

      } else {
	if (want_low_p == true) {
	  for (qq = q; qq != advance_one_hit(q); qq = List_next(qq)) {
	    process_line((char *) List_head(qq),header2,refnts,nmatches,
			 matches_byshift,matches_byquality,
			 mismatches_byshift,mismatches_byquality,
			 quality_counts_match,quality_counts_mismatch,
			 chromosome,chromosome_iit,genome,chroffset);
	  }
	}
	if (want_high_p == true) {
	  for (pp = p; pp != advance_one_hit(p); pp = List_next(pp)) {
	    process_line((char *) List_head(pp),header1,refnts,nmatches,
			 matches_byshift,matches_byquality,
			 mismatches_byshift,mismatches_byquality,
			 quality_counts_match,quality_counts_mismatch,
			 chromosome,chromosome_iit,genome,chroffset);
	  }
	}
      }

      FREE(chr1);
      FREE(chr2);

      p = advance_one_hit(p);
      q = advance_one_hit(q);
    }
  }

  lines_gc(&lines2);
  lines_gc(&lines1);
  FREE(header2);
  FREE(header1);

  return;
}

#endif	/* GSNAP_INPUT */


/************************************************************************
 *   Probability calculations
 ************************************************************************/

#define MAX_QUALITY_SCORE 40

#define NGENOTYPES 10

static double
allele_logprob[3][MAX_QUALITY_SCORE+1] = {
  /* Not in genotype: log(1/3*10^(-Q/10)) */
  {-1.098612,
   -1.328871, -1.559129, -1.789388, -2.019646, -2.249905,
   -2.480163, -2.710422, -2.940680, -3.170939, -3.401197,
   -3.631456, -3.861714, -4.091973, -4.322231, -4.552490,
   -4.782748, -5.013007, -5.243265, -5.473524, -5.703782,
   -5.934041, -6.164299, -6.394558, -6.624817, -6.855075,
   -7.085334, -7.315592, -7.545851, -7.776109, -8.006368,
   -8.236626, -8.466885, -8.697143, -8.927402, -9.157660,
   -9.387919, -9.618177, -9.848436, -10.078694, -10.308953},

  /* In heterozygote: log(1/2 - 1/3*10^(-Q/10)) */
  {-1.7917595,
   -1.4472174, -1.2389754, -1.0998002, -1.0015828, -0.9299061,
   -0.8764201, -0.8358837, -0.8048159, -0.7808079, -0.7621401,
   -0.7475561, -0.7361213, -0.7271306, -0.7200462, -0.7144544,
   -0.7100349, -0.7065382, -0.7037694, -0.7015754, -0.6998362,
   -0.6984568, -0.6973624, -0.6964940, -0.6958048, -0.6952576,
   -0.6948232, -0.6944782, -0.6942043, -0.6939868, -0.6938141,
   -0.6936769, -0.6935679, -0.6934814, -0.6934126, -0.6933580,
   -0.6933147, -0.6932802, -0.6932528, -0.6932311, -0.6932138},

  /* In homozygote: log(1 - 10^(-Q/10)) */
  {/*-Inf*/-1.58,
   -1.5814737534, -0.9968430440, -0.6955244713, -0.5076758737, -0.3801304081,
   -0.2892681872, -0.2225515160, -0.1725565729, -0.1345519603, -0.1053605157,
   -0.0827653027, -0.0651741732, -0.0514182742, -0.0406248442, -0.0321335740,
   -0.0254397275, -0.0201543648, -0.0159758692, -0.0126691702, -0.0100503359,
   -0.0079749983, -0.0063295629, -0.0050244739, -0.0039890173, -0.0031672882,
   -0.0025150465, -0.0019972555, -0.0015861505, -0.0012597185, -0.0010005003,
   -0.0007946439, -0.0006311565, -0.0005013129, -0.0003981864, -0.0003162778,
   -0.0002512202, -0.0001995461, -0.0001585019, -0.0001259005, -0.0001000050}
};


static void
make_quality_scores_constant (int quality_score_constant) {
  int i, Q;
  double value;

  for (i = 0; i < 3; i++) {
    value = allele_logprob[i][quality_score_constant];
    for (Q = 0; Q <= MAX_QUALITY_SCORE; Q++) {
      allele_logprob[i][Q] = value;
    }
  }

  return;
}


static int
genotype_compatibility[4][NGENOTYPES] = {
  
  {/*AA*/ 2, /*AC*/ 1, /*AG*/ 1, /*AT*/ 1,
   /*CC*/ 0, /*CG*/ 0, /*CT*/ 0,
   /*GG*/ 0, /*GT*/ 0,
   /*TT*/ 0},

  /* C */
  {/*AA*/ 0, /*AC*/ 1, /*AG*/ 0, /*AT*/ 0,
   /*CC*/ 2, /*CG*/ 1, /*CT*/ 1,
   /*GG*/ 0, /*GT*/ 0,
   /*TT*/ 0},

  /* G */
  {/*AA*/ 0, /*AC*/ 0, /*AG*/ 1, /*AT*/ 0,
   /*CC*/ 0, /*CG*/ 1, /*CT*/ 0,
   /*GG*/ 2, /*GT*/ 1,
   /*TT*/ 0},

  /* T */
  {/*AA*/ 0, /*AC*/ 0, /*AG*/ 0, /*AT*/ 1,
   /*CC*/ 0, /*CG*/ 0, /*CT*/ 1,
   /*GG*/ 0, /*GT*/ 1,
   /*TT*/ 2}
};


static void
compute_probs (double *probs, double *loglik, char refnt,
	       List_T matches_byquality, List_T mismatches_byquality) {
  double total, adj_loglik[NGENOTYPES], maxlik;
  List_T p;
  Match_T match;
  Mismatch_T mismatch;
  int Q;
  int g;

  for (g = 0; g < NGENOTYPES; g++) {
    loglik[g] = 0.0;
  }

  for (p = matches_byquality; p != NULL; p = List_next(p)) {
    match = (Match_T) List_head(p);
    Q = match->quality - quality_score_adj;
    if (Q < 1) {
      Q = 1;
    } else if (Q > MAX_QUALITY_SCORE) {
      Q = MAX_QUALITY_SCORE;
    }

    switch (refnt) {
    case 'A':
      for (g = 0; g < NGENOTYPES; g++) {
	loglik[g] += match->count * allele_logprob[genotype_compatibility[0][g]][Q];
      }
      break;

    case 'C':
      for (g = 0; g < NGENOTYPES; g++) {
	loglik[g] += match->count * allele_logprob[genotype_compatibility[1][g]][Q];
      }
      break;

    case 'G':
      for (g = 0; g < NGENOTYPES; g++) {
	loglik[g] += match->count * allele_logprob[genotype_compatibility[2][g]][Q];
      }
      break;

    case 'T':
      for (g = 0; g < NGENOTYPES; g++) {
	loglik[g] += match->count * allele_logprob[genotype_compatibility[3][g]][Q];
      }
      break;
    }
  }

  for (p = mismatches_byquality; p != NULL; p = List_next(p)) {
    mismatch = (Mismatch_T) List_head(p);
    Q = mismatch->quality - quality_score_adj;
    if (Q < 1) {
      Q = 1;
    } else if (Q > MAX_QUALITY_SCORE) {
      Q = MAX_QUALITY_SCORE;
    }

    switch (mismatch->nt) {
    case 'A':
      for (g = 0; g < NGENOTYPES; g++) {
	loglik[g] += mismatch->count * allele_logprob[genotype_compatibility[0][g]][Q];
      }
      break;

    case 'C':
      for (g = 0; g < NGENOTYPES; g++) {
	loglik[g] += mismatch->count * allele_logprob[genotype_compatibility[1][g]][Q];
      }
      break;

    case 'G':
      for (g = 0; g < NGENOTYPES; g++) {
	loglik[g] += mismatch->count * allele_logprob[genotype_compatibility[2][g]][Q];
      }
      break;

    case 'T':
      for (g = 0; g < NGENOTYPES; g++) {
	loglik[g] += mismatch->count * allele_logprob[genotype_compatibility[3][g]][Q];
      }
      break;
    }
  }
  
  /* Raise all loglikelihoods to maximum, to avoid underflow when taking exp() */
  maxlik = loglik[0];
  for (g = 1; g < NGENOTYPES; g++) {
    if (loglik[g] > maxlik) {
      maxlik = loglik[g];
    }
  }

  for (g = 0; g < NGENOTYPES; g++) {
    adj_loglik[g] = loglik[g] - maxlik;
  }


  total = 0.0;
  for (g = 0; g < NGENOTYPES; g++) {
    probs[g] = exp(adj_loglik[g]);
    total += probs[g];
  }

  for (g = 0; g < NGENOTYPES; g++) {
    probs[g] /= total;
  }

  return;
}




static void
print_chromosome_blocks (Genomicpos_T allocoffset, Genomicpos_T chroffset,
			 Genomicpos_T chrstart, Genomicpos_T chrend,
			 char *refnts, int *nmatches, 
			 List_T *matches_byshift, List_T *matches_byquality,
			 List_T *mismatches_byshift, List_T *mismatches_byquality,
			 Genome_T genome, IIT_T chromosome_iit, char *chr) {
  double probs[NGENOTYPES], loglik[NGENOTYPES];
  int length, i, j, g;
  Genomicpos_T position, pos, posi, firsti, lasti;
  bool hitp;
  Match_T *match_array;
  Mismatch_T mismatch, mismatch0, *array, *subarray;
  List_T unique_mismatches_byshift, unique_mismatches_byquality, ptr;
  long int total;
 
  for (pos = chrstart - 1U; pos < chrend; pos += blocksize) {
    hitp = false;
    for (posi = pos; posi < chrend && posi < pos + blocksize; posi += 1U) {
      if (nmatches[allocoffset+posi] > 0 || mismatches_byshift[allocoffset+posi] != NULL) {
	if (hitp == false) {
	  firsti = posi;
	  hitp = true;
	}
	lasti = posi;
      }
    }
      
    if (hitp == true) {
      /* old: printf(">%s_%d %u %u %s\n",chromosome,++segno,firsti+1U,lasti+1U,chromosome); */
      /* Block total */
      total = 0;
      for (posi = firsti; posi <= lasti; posi++) {
	total += nmatches[allocoffset+posi];
	for (ptr = mismatches_byshift[allocoffset+posi]; ptr != NULL; ptr = List_next(ptr)) {
	  mismatch = (Mismatch_T) List_head(ptr);
	  total += mismatch->count;
	}
      }
      printf(">%ld %s:%u..%u\n",total,chr,firsti+1U,lasti+1U);

      for (posi = firsti; posi <= lasti; posi++) {
	position = allocoffset+posi;

#if 0
	printf("%s:%u\t",chr,posi+1U);

	/* Positional total */
	total = nmatches[allocoffset+posi];
	for (ptr = mismatches_byshift[allocoffset+posi]; ptr != NULL; ptr = List_next(ptr)) {
	  mismatch = (Mismatch_T) List_head(ptr);
	  total += mismatch->count;
	}
	printf("%ld\t",total);
#endif


	/* Details */
	if (nmatches[position] == 0) {
	  printf("%c0",Genome_get_char(genome,chroffset+posi));
	} else {
	  printf("%c%d",refnts[position],nmatches[position]);

	  if (print_refdetails_p == true) {
	    printf("(");
	    length = List_length(matches_byshift[position]);
	    match_array = (Match_T *) List_to_array(matches_byshift[position],NULL);
	    qsort(match_array,length,sizeof(Match_T),Match_shift_cmp);
	    printf("%d@%d",match_array[0]->count,match_array[0]->shift);
	    for (j = 1; j < length; j++) {
	      printf(",%d@%d",match_array[j]->count,match_array[j]->shift);
	    }
	    FREE(match_array);

	    printf("|");

	    length = List_length(matches_byquality[position]);
	    match_array = (Match_T *) List_to_array(matches_byquality[position],NULL);
	    qsort(match_array,length,sizeof(Match_T),Match_quality_cmp);
	    printf("%dQ%d",match_array[0]->count,match_array[0]->quality - quality_score_adj);
	    for (j = 1; j < length; j++) {
	      printf(",%dQ%d",match_array[j]->count,match_array[j]->quality - quality_score_adj);
	    }
	    FREE(match_array);

	    printf(")");
	  }
	}

	if (mismatches_byshift[position] != NULL) {
	  unique_mismatches_byshift = NULL;
	  for (ptr = mismatches_byshift[position]; ptr != NULL; ptr = List_next(ptr)) {
	    mismatch = (Mismatch_T) List_head(ptr);
	    if ((mismatch0 = find_mismatch_nt(unique_mismatches_byshift,mismatch->nt)) != NULL) {
	      mismatch0->count += mismatch->count;
		
	      /* Insert mismatch into list */
	      mismatch->next = mismatch0->next;
	      mismatch0->next = mismatch;

	      mismatch0->shift += 1; /* Used here as nshifts */
	    } else {
	      mismatch0 = Mismatch_new(mismatch->nt,/*shift, used here as nshifts*/1,' ');
	      mismatch0->count = mismatch->count;
	      mismatch0->next = mismatch;
	      unique_mismatches_byshift = List_push(unique_mismatches_byshift,mismatch0);
	    }
	  }

	  unique_mismatches_byquality = NULL;
	  for (ptr = mismatches_byquality[position]; ptr != NULL; ptr = List_next(ptr)) {
	    mismatch = (Mismatch_T) List_head(ptr);
	    if ((mismatch0 = find_mismatch_nt(unique_mismatches_byquality,mismatch->nt)) != NULL) {
	      mismatch0->count += mismatch->count;
		
	      /* Insert mismatch into list */
	      mismatch->next = mismatch0->next;
	      mismatch0->next = mismatch;

	      mismatch0->shift += 1; /* Used here as nshifts */
	    } else {
	      mismatch0 = Mismatch_new(mismatch->nt,/*shift, used here as nshifts*/1,' ');
	      mismatch0->count = mismatch->count;
	      mismatch0->next = mismatch;
	      unique_mismatches_byquality = List_push(unique_mismatches_byquality,mismatch0);
	    }
	  }

	  array = (Mismatch_T *) List_to_array(unique_mismatches_byshift,NULL);
	  qsort(array,List_length(unique_mismatches_byshift),sizeof(Mismatch_T),Mismatch_count_cmp);

	  for (i = 0; i < List_length(unique_mismatches_byshift); i++) {
	    mismatch0 = array[i];
	    printf(" %c%d",mismatch0->nt,mismatch0->count);
	    if (print_refdetails_p == true) {
	      printf("(");
	      
	      length = Mismatch_chain_length(mismatch0->next);
	      subarray = (Mismatch_T *) CALLOC(length,sizeof(Mismatch_T));
	      for (mismatch = mismatch0->next, j = 0; mismatch != NULL; mismatch = mismatch->next, j++) {
		subarray[j] = mismatch;
	      }
	      qsort(subarray,length,sizeof(Mismatch_T),Mismatch_shift_cmp);

	      printf("%d@%d",subarray[0]->count,subarray[0]->shift);
	      for (j = 1; j < length; j++) {
		printf(",%d@%d",subarray[j]->count,subarray[j]->shift);
	      }
	      FREE(subarray);

	      printf("|");
	      
	      mismatch0 = find_mismatch_nt(unique_mismatches_byquality,mismatch0->nt);
	      length = Mismatch_chain_length(mismatch0->next);
	      subarray = (Mismatch_T *) CALLOC(length,sizeof(Mismatch_T));
	      for (mismatch = mismatch0->next, j = 0; mismatch != NULL; mismatch = mismatch->next, j++) {
		subarray[j] = mismatch;
	      }
	      qsort(subarray,length,sizeof(Mismatch_T),Mismatch_quality_cmp);
	      
	      printf("%dQ%d",subarray[0]->count,subarray[0]->quality - quality_score_adj);
	      for (j = 1; j < length; j++) {
		printf(",%dQ%d",subarray[j]->count,subarray[j]->quality - quality_score_adj);
	      }
	      FREE(subarray);

	      printf(")");
	    }
	  }

	  FREE(array);

	  for (ptr = unique_mismatches_byshift; ptr != NULL; ptr = List_next(ptr)) {
	    mismatch0 = List_head(ptr);
	    Mismatch_free(&mismatch0);
	  }
	  List_free(&unique_mismatches_byshift);

	  for (ptr = unique_mismatches_byquality; ptr != NULL; ptr = List_next(ptr)) {
	    mismatch0 = List_head(ptr);
	    Mismatch_free(&mismatch0);
	  }
	  List_free(&unique_mismatches_byquality);
	}

	if (want_genotypes_p == true) {
	  compute_probs(probs,loglik,refnts[position],matches_byquality[position],
			mismatches_byquality[position]);

	  /* Probabilities */
	  printf("\t");
	  printf("%.3g",1.0-probs[0]);
	  for (g = 1; g < NGENOTYPES; g++) {
	    printf(",%.3g",1.0-probs[g]);
	  }

	  /* Log-likelihoods */
	  printf("\t");
	  printf("%.3g",loglik[0]);
	  for (g = 1; g < NGENOTYPES; g++) {
	    printf(",%.3g",loglik[g]);
	  }
	}

	printf("\n");
      }
    }
  }

  return;
}


static void
print_chromosome_blocks_wig (Genomicpos_T allocoffset, Genomicpos_T chroffset,
			     Genomicpos_T chrstart, Genomicpos_T chrend,
			     char *refnts, int *nmatches, 
			     List_T *matches_byshift, List_T *matches_byquality,
			     List_T *mismatches_byshift, List_T *mismatches_byquality,
			     Genome_T genome, IIT_T chromosome_iit, char *chr) {
  Genomicpos_T pos, posi, firsti, lasti;
  bool hitp;
  Mismatch_T mismatch;
  List_T ptr;
  long int total;
 
  for (pos = chrstart - 1U; pos < chrend; pos += blocksize) {
    hitp = false;
    for (posi = pos; posi < chrend && posi < pos + blocksize; posi += 1U) {
      if (nmatches[allocoffset+posi] > 0 || mismatches_byshift[allocoffset+posi] != NULL) {
	if (hitp == false) {
	  firsti = posi;
	  hitp = true;
	}
	lasti = posi;
      }
    }
      
    if (hitp == true) {
      /* old: printf(">%s_%d %u %u %s\n",chromosome,++segno,firsti+1U,lasti+1U,chromosome); */
      /* Block total */
      total = 0;
      for (posi = firsti; posi <= lasti; posi++) {
	total += nmatches[allocoffset+posi];
	for (ptr = mismatches_byshift[allocoffset+posi]; ptr != NULL; ptr = List_next(ptr)) {
	  mismatch = (Mismatch_T) List_head(ptr);
	  total += mismatch->count;
	}
      }
      printf("variableStep chrom=%s start=%u\n",chr,firsti+1U);

      for (posi = firsti; posi <= lasti; posi++) {
	printf("%u",posi+1U);

	/* Positional total */
	total = nmatches[allocoffset+posi];
	for (ptr = mismatches_byshift[allocoffset+posi]; ptr != NULL; ptr = List_next(ptr)) {
	  mismatch = (Mismatch_T) List_head(ptr);
	  total += mismatch->count;
	}
	printf(" %ld",total);

	printf("\n");
      }
    }
  }

  return;
}



int
main (int argc, char *argv[]) {
  char *genomesubdir = NULL, *fileroot = NULL, *iitfile;
  IIT_T chromosome_iit;
  int index;
  int quality_counts_match[256], quality_counts_mismatch[256], i;

  Genomicpos_T allocoffset = 0U, alloclength, chroffset = 0U, chrstart, chrend, chrlength, pos;
  char *chr;
  bool allocp;

  char *refnts;

  int *nmatches;
  List_T *matches_byshift = NULL, *matches_byquality = NULL,
    *mismatches_byshift, *mismatches_byquality, p;
  Match_T match;
  Mismatch_T mismatch;

  Genome_T genome = NULL;

#ifdef SAM_INPUT
#elif defined(BAM_INPUT)
  Bamreader_T bamreader;
#endif

#ifndef SAM_INPUT
  Genomicpos_T genomicstart, genomiclength;
  bool revcomp;
#endif


  int opt;
  extern int optind;
  extern char *optarg;
  int long_option_index = 0;
  const char *long_name;

  while ((opt = getopt_long(argc,argv,
			    "D:d:2AGQ:c:FRs:C:U:b:a:h:T:l:",
			    long_options, &long_option_index)) != -1) {
    switch (opt) {
    case 0:
      long_name = long_options[long_option_index].name;
      if (!strcmp(long_name,"version")) {
	print_program_version();
	exit(0);
      } else if (!strcmp(long_name,"help")) {
	print_program_usage();
	exit(0);
      } else if (!strcmp(long_name,"use-quality-const")) {
	quality_score_constant = atoi(optarg);
	if (quality_score_constant > MAX_QUALITY_SCORE) {
	  fprintf(stderr,"quality substitution score %d is > maximum %d\n",
		  quality_score_constant,MAX_QUALITY_SCORE);
	  exit(9);
	}
      } else {
	/* Shouldn't reach here */
	fprintf(stderr,"Don't recognize option %s.  For usage, run 'gsnap --help'",long_name);
	exit(9);
      }
      break;

    case 'D': user_genomedir = optarg; break;
    case 'd': 
      dbroot = (char *) CALLOC(strlen(optarg)+1,sizeof(char));
      strcpy(dbroot,optarg);
      break;
    case '2': dibasep = true; break;

    case 'A': print_refdetails_p = true; break;
    case 'G': want_genotypes_p = true; break;
    case 'Q': minimum_mapq = atoi(optarg); break;

    case 'c': chromosome = optarg; break;

    case 'F': 
      desired_strand = '+';
      opposite_strand = '-';
      break;

    case 'R':
      desired_strand = '-';
      opposite_strand = '+';
      break;

    case 's':
      if (!strcmp(optarg,"low")) {
	want_high_p = false;
      } else if (!strcmp(optarg,"high")) {
	want_low_p = false;
      } else {
	fprintf(stderr,"Side %s not recognized.  Must be low or high\n",optarg);
	exit(9);
      }
      break;

    case 'C':
      switch (atoi(optarg)) {
      case 0: need_concordant_p = false; break;
      case 1: need_concordant_p = true; break;
      default: fprintf(stderr,"Concordant mode %s not recognized.\n",optarg); exit(9);
      }
      break;

    case 'U':
      switch (atoi(optarg)) {
      case 0: uniquep = false; break;
      case 1: uniquep = true; break;
      default: fprintf(stderr,"Unique mode %s not recognized.\n",optarg); exit(9);
      }
      break;

    case 'b': blocksize = atoi(optarg); break;

    case 'a':
      if (!strcmp(optarg,"all")) {
	mode = ALL_SCORES;
      } else if (!strcmp(optarg,"caps")) {
	mode = HIGH_SCORES;
      } else if (!strcmp(optarg,"mismatch_caps")) {
	mode = MISMATCH_HIGH_SCORES;
      } else {
	fprintf(stderr,"Mode %s not recognized.\n",optarg); exit(9);
      }
      break;

    case 'h': haltpos = strtoul(optarg,NULL,10); break;

    case 'T': 
      switch (atoi(optarg)) {
      case 0: trimp = false; break;
      case 1: trimp = true; break;
      default: fprintf(stderr,"Trim mode %s not recognized.\n",optarg); exit(9);
      }
      break;

#ifdef SAM_INPUT
    case 'l': min_mlength = atoi(optarg); break;
#elif defined(BAM_INPUT)
    case 'l': min_mlength = atoi(optarg); break;
#else
    case 'l': min_readlength = atoi(optarg); break;
#endif

    default: exit(9);
    }
  }
  argc -= optind;
  argv += optind;
      
  if (dbroot == NULL) {
    fprintf(stderr,"Need to specify the -d flag\n");
    print_program_usage();
    exit(9);
  }

  if (
#ifdef BAM_INPUT
      argc == 1 &&
#endif
      chromosome == NULL) {
    fprintf(stderr,"Processing entire genome...\n");
  }

  genomesubdir = Datadir_find_genomesubdir(&fileroot,&dbversion,user_genomedir,dbroot);

  iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			    strlen(fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
  sprintf(iitfile,"%s/%s.chromosome.iit",genomesubdir,fileroot);
  chromosome_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
			    /*divstring*/NULL,/*add_iit_p*/false,/*labels_read_p*/true);
  FREE(iitfile);

  /* Need genome to determine wild-type, because "known gene" may not match reference genome */
  genome = Genome_new(genomesubdir,fileroot,/*snps_root*/NULL,/*uncompressedp*/false,
		      /*access*/USE_MMAP_ONLY);


  if (dibasep == true && trimp == true) {
    fprintf(stderr,"Turning trimming off for 2-base encoded alignments\n");
    trimp = false;
  }


  if (chromosome != NULL) {
    if ((index = IIT_find_one(chromosome_iit,chromosome)) < 0) {
      fprintf(stderr,"Cannot find chromosome %s in chromosome IIT file\n",chromosome);
      exit(9);

    } else {
      chroffset = Interval_low(IIT_interval(chromosome_iit,index));
      alloclength = Interval_length(IIT_interval(chromosome_iit,index));
      fprintf(stderr,"Starting to allocate memory for %u positions\n",alloclength);

      chrstart = 1;
      chrend = alloclength;

      refnts = (char *) CALLOC(alloclength,sizeof(char));
      nmatches = (int *) CALLOC(alloclength,sizeof(int));
      if (print_refdetails_p == true) {
	matches_byshift = (List_T *) CALLOC(alloclength,sizeof(List_T));
      }
      matches_byquality = (List_T *) CALLOC(alloclength,sizeof(List_T));
      mismatches_byshift = (List_T *) CALLOC(alloclength,sizeof(List_T));
      mismatches_byquality = (List_T *) CALLOC(alloclength,sizeof(List_T));
    }

#ifdef BAM_INPUT
  } else if (argc > 1) {
    Parserange_universal(&chromosome,&revcomp,&genomicstart,&genomiclength,&chrstart,&chrend,
			 &chroffset,&chrlength,argv[1],genomesubdir,fileroot);

    alloclength = chrlength;
    fprintf(stderr,"Starting to allocate memory for %u positions (could be improved)\n",alloclength);
    refnts = (char *) CALLOC(alloclength,sizeof(char));
    nmatches = (int *) CALLOC(alloclength,sizeof(int));
    if (print_refdetails_p == true) {
      matches_byshift = (List_T *) CALLOC(alloclength,sizeof(List_T));
    }
    matches_byquality = (List_T *) CALLOC(alloclength,sizeof(List_T));
    mismatches_byshift = (List_T *) CALLOC(alloclength,sizeof(List_T));
    mismatches_byquality = (List_T *) CALLOC(alloclength,sizeof(List_T));
#endif

  } else {
    alloclength = IIT_totallength(chromosome_iit);
    fprintf(stderr,"Starting to allocate memory for all %u positions\n",alloclength);

    refnts = (char *) CALLOC(alloclength,sizeof(char));
    nmatches = (int *) CALLOC(alloclength,sizeof(int));
    if (print_refdetails_p == true) {
      matches_byshift = (List_T *) CALLOC(alloclength,sizeof(List_T));
    }
    matches_byquality = (List_T *) CALLOC(alloclength,sizeof(List_T));
    mismatches_byshift = (List_T *) CALLOC(alloclength,sizeof(List_T));
    mismatches_byquality = (List_T *) CALLOC(alloclength,sizeof(List_T));

  }

  fprintf(stderr,"Done allocating memory\n");

  for (i = 0; i < 256; i++) {
    quality_counts_match[i] = 0;
    quality_counts_mismatch[i] = 0;
  }


  FREE(dbroot);
  FREE(dbversion);
  FREE(fileroot);
  FREE(genomesubdir);


  /* */

  if (want_low_p == false || want_high_p == false) {
    need_concordant_p = true;
  } else {
    need_concordant_p = false;
  }

#ifdef SAM_INPUT
  parse_sam(refnts,nmatches,matches_byshift,matches_byquality,
	    mismatches_byshift,mismatches_byquality,
	    quality_counts_match,quality_counts_mismatch,
	    chromosome,chromosome_iit,genome,chroffset);

#elif defined(BAM_INPUT)

  bamreader = Bamread_new(argv[0]);
  if (argc > 1) {
    Bamread_limit_region(bamreader,chromosome,chrstart,chrend);
  }

  parse_bam(refnts,nmatches,matches_byshift,matches_byquality,
	    mismatches_byshift,mismatches_byquality,
	    quality_counts_match,quality_counts_mismatch,
	    bamreader,chromosome,chromosome_iit,genome,chroffset);
  /* Bamread_unlimit_region(bamreader); */

  /* free bamreader */

#else
  if (need_concordant_p == true) {
    parse_paired(refnts,nmatches,matches_byshift,matches_byquality,
		 mismatches_byshift,mismatches_byquality,
		 quality_counts_match,quality_counts_mismatch,
		 chromosome,chromosome_iit,genome,chroffset);
  } else {
    parse_single(refnts,nmatches,matches_byshift,matches_byquality,
		 mismatches_byshift,mismatches_byquality,
		 quality_counts_match,quality_counts_mismatch,
		 chromosome,chromosome_iit,genome,chroffset);
  }
#endif

  if (quality_score_constant > 0) {
    make_quality_scores_constant(quality_score_constant);
  }

  if (chromosome != NULL) {
    /* Print results for one chromosome */
    if (print_wig_p == true) {
      print_chromosome_blocks_wig(/*allocoffset*/0,chroffset,chrstart,chrend,refnts,
				  nmatches,matches_byshift,matches_byquality,
				  mismatches_byshift,mismatches_byquality,
				  genome,chromosome_iit,/*chr*/chromosome);
    } else if (print_allele_counts_p == true) {
      print_chromosome_blocks(/*allocoffset*/0,chroffset,chrstart,chrend,refnts,
			      nmatches,matches_byshift,matches_byquality,
			      mismatches_byshift,mismatches_byquality,
			      genome,chromosome_iit,/*chr*/chromosome);
    }

  } else {
    /* Print results for all chromosomes */
    for (index = 1; index <= IIT_total_nintervals(chromosome_iit); index++) {
      allocoffset = Interval_low(IIT_interval(chromosome_iit,index));
      chrlength = Interval_length(IIT_interval(chromosome_iit,index));
      chr = IIT_label(chromosome_iit,index,&allocp);
      fprintf(stderr,"%s\t%u\t%u\n",chr,allocoffset,chrlength);
      if (print_wig_p == true) {
	print_chromosome_blocks_wig(allocoffset,/*chroffset*/allocoffset,
				    /*chrstart*/1,/*chrend*/chrlength,refnts,
				    nmatches,matches_byshift,matches_byquality,
				    mismatches_byshift,mismatches_byquality,
				    genome,chromosome_iit,chr);
      } else if (print_allele_counts_p == true) {
	print_chromosome_blocks(allocoffset,/*chroffset*/allocoffset,
				/*chrstart*/1,/*chrend*/chrlength,refnts,
				nmatches,matches_byshift,matches_byquality,
				mismatches_byshift,mismatches_byquality,
				genome,chromosome_iit,chr);
      }
      if (allocp == true) {
	FREE(chr);
      }
    }
  }


  for (i = 0; i < 256; i++) {
    if (quality_counts_match[i] > 0 || quality_counts_mismatch[i] > 0) {
      fprintf(stderr,"Quality %c %d: %d %d\n",i,i-quality_score_adj,quality_counts_match[i],quality_counts_mismatch[i]);
    }
  }

  for (pos = 0; pos < alloclength; pos++) {
    for (p = mismatches_byshift[pos]; p != NULL; p = List_next(p)) {
      mismatch = (Mismatch_T) List_head(p);
      Mismatch_free(&mismatch);
    }
    List_free(&(mismatches_byshift[pos]));
  }
  FREE(mismatches_byshift);

  for (pos = 0; pos < alloclength; pos++) {
    for (p = mismatches_byquality[pos]; p != NULL; p = List_next(p)) {
      mismatch = (Mismatch_T) List_head(p);
      Mismatch_free(&mismatch);
    }
    List_free(&(mismatches_byquality[pos]));
  }
  FREE(mismatches_byquality);

  if (print_refdetails_p == true) {
    for (pos = 0; pos < alloclength; pos++) {
      for (p = matches_byshift[pos]; p != NULL; p = List_next(p)) {
	match = (Match_T) List_head(p);
	Match_free(&match);
      }
      List_free(&(matches_byshift[pos]));
    }
    FREE(matches_byshift);
  }

  for (pos = 0; pos < alloclength; pos++) {
    for (p = matches_byquality[pos]; p != NULL; p = List_next(p)) {
      match = (Match_T) List_head(p);
      Match_free(&match);
    }
    List_free(&(matches_byquality[pos]));
  }
  FREE(matches_byquality);

  FREE(nmatches);
  FREE(refnts);

  IIT_free(&chromosome_iit);

  if (genome != NULL) {
    Genome_free(&genome);
  }

  return 0;
}


static void
print_program_usage () {
    fprintf(stdout,"\
Usage: cat <GSNAP output> | gsnap_tally [OPTIONS...]\n\
\n\
Input options (must include -d and -c)\n\
  -D, --dir=directory            Genome directory\n\
  -d, --db=STRING                Genome database\n\
  -c, --chr=STRING               Chromosome\n\
\n\
  -G, --genotypes                Print genotype information (probability and likelihood)\n\
  -A, --alldetails               Include details about shifts and quality scores\n\
  --use-quality-const=INT        Ignore quality scores and use this value instead\n\
  -Q, --min-mapq=INT             Use only alignments with this mapping quality and higher\n\
\n\
  -F, --fwd                      Forward strand only\n\
  -R, --rev                      Reverse complement strand only\n\
  -s, --side=STRING              Side of paired-end read: low or high\n\
  -C, --concordant=INT           Concordant hits only (0=no [default], 1=yes)\n\
  -U, --unique=INT               Unique hits only (0=no [default], 1=yes)\n\
  -b, --blocksize=INT            Block size [default 1000]\n\
  -a, --mode=STRING              Mode (all, caps [default], mismatch_caps)\n\
  -l, --minlength=INT            Minimum read length (default 10).  Allows skipping of indel ends\n\
\n\
");
    return;
}
