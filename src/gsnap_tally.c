static char rcsid[] = "$Id: gsnap_tally.c,v 1.8 2010/03/08 20:21:06 twu Exp $";
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


#include "except.h"
#include "mem.h"
#include "bool.h"
#include "genomicpos.h"
#include "complement.h"
#include "list.h"
#include "iit-read.h"
#include "interval.h"
#include "genome.h"
#include "datadir.h"
#include "getopt.h"


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

#define MATCH_SCORE 1
#define MISMATCH_SCORE -4

static char *user_genomedir = NULL;
static char *dbroot = NULL;
static char *dbversion = NULL;

static char *chromosome = NULL;

static bool remappedp = false;

static bool uniquep = false;
static int blocksize = 1000;
static int min_readlength = 10;
static Genomicpos_T haltpos = -1U;

static bool trimp = true;
static bool trim_endnt_p = false;
static bool dibasep = false;

static char desired_strand = ' ';
static char opposite_strand = ' ';

static bool want_first_p = true;
static bool want_second_p = true;


typedef enum {FIVE, THREE} Pairend_T;
typedef enum {ALL_SCORES, MISMATCH_HIGH_SCORES, HIGH_SCORES, ONE_SUB, ONE_SUB_HIGH_SCORES} Mode_T;

static Mode_T mode = HIGH_SCORES;


static struct option long_options[] = {
  /* Input options */
  {"dir", required_argument, 0, 'D'},	/* user_genomedir */
  {"db", required_argument, 0, 'd'}, /* dbroot */

  {"chr", required_argument, 0, 'c'}, /* chromosome */
  {"fwd", no_argument, 0, 'F'}, /* desired_strand, opposite_strand */
  {"rev", no_argument, 0, 'R'}, /* desired_strand, opposite_strand */
  {"foo", required_argument, 0, 'P'}, /* want_first_p, want_second_p */

  {"remapped", no_argument, 0, 'M'}, /* remappedp */
  {"unique", required_argument, 0, 'u'}, /* uniquep */
  {"blocksize", required_argument, 0, 'b'}, /* blocksize */

  {"mode", required_argument, 0, 'a'}, /* mode */

  {"halt", required_argument, 0, 'h'}, /* haltpos */

  {"trim", required_argument, 0, 'T'}, /* trimp */
  {"trimendnt", required_argument, 0, 't'}, /* trim_endnt_p */
  {"minlength", required_argument, 0, 'l'}, /* min_readlength */
  
  {"dibase", no_argument, 0, '2'}, /* dibasep */

  /* Help options */
  {"version", no_argument, 0, 'V'}, /* print_program_version */
  {"help", no_argument, 0, '?'}, /* print_program_usage */
  {0, 0, 0, 0}
};


static void
print_program_version () {
  fprintf(stdout,"\n");
  fprintf(stdout,"GSNAP_TALLY\n");
  fprintf(stdout,"Part of GMAP package, version %s\n",PACKAGE_VERSION);
  fprintf(stdout,"Build target: %s\n",TARGET);
  fprintf(stdout,"Default gmap directory: %s\n",GMAPDB);
  fprintf(stdout,"Thomas D. Wu, Genentech, Inc.\n");
  fprintf(stdout,"Contact: twu@gene.com\n");
  fprintf(stdout,"\n");
  return;
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
  -F, --fwd                      Forward strand only\n\
  -R, --rev                      Reverse complement strand only\n\
  -M, --remapped                 Remapped from cDNA to genome\n\
  -u, --unique=INT               Unique hits only (0=no [default], 1=yes)\n\
  -b, --blocksize=INT            Block size [default 1000]\n\
  -a, --mode=STRING              Mode (high [default], all, mismatch_high, onesub, onesub_high)\n\
  -h, --halt=INT                 Halt at given position\n\
  -T, --trim=INT                 Trim poor matches on left and right (default 1)\n\
  -t, --trimendnt=INT            Always trim first and last nt in query (default 0)\n\
  -l, --minlength=INT            Minimum read length (default 10).  Allows skipping of indel ends\n\
\n\
");
    return;
}

typedef struct Mismatch_T *Mismatch_T;
struct Mismatch_T {
  char nt;
  int querypos;			/* Used to record shifts */
  int count;
};


static Mismatch_T
Mismatch_new (char nt, int querypos) {
  Mismatch_T new = (Mismatch_T) MALLOC(sizeof(*new));

  new->nt = nt;
  new->querypos = querypos;
  new->count = 1;
  return new;
}

static void
Mismatch_free (Mismatch_T *old) {
  FREE(*old);
  return;
}

static int
Mismatch_cmp (const void *a, const void *b) {
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

static Mismatch_T
find_mismatch (List_T mismatches, char nt, int querypos) {
  List_T p;
  Mismatch_T mismatch;

  for (p = mismatches; p != NULL; p = List_next(p)) {
    mismatch = (Mismatch_T) List_head(p);
    if (mismatch->nt == nt && mismatch->querypos == querypos) {
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
compute_trims (int *trimleft, int *trimright, char *query, char *genomic, int readlength,
	       char strand, int query5, int query3, int querylength) {
  int score, pos;

  score = 0;
  *trimleft = 0;
  for (pos = 0; pos < readlength; pos++) {
    if (toupper(query[pos]) == toupper(genomic[pos])) {
      score += MATCH_SCORE;
    } else {
      score += MISMATCH_SCORE;
      if (score < 0) {
	score = 0;
	*trimleft = pos + 1;
      }
    }
  }

  score = 0;
  *trimright = 0;
  for (pos = readlength-1; pos >= 0; pos--) {
    if (toupper(query[pos]) == toupper(genomic[pos])) {
      score += MATCH_SCORE;
    } else {
      score += MISMATCH_SCORE;
      if (score < 0) {
	score = 0;
	*trimright = readlength - pos;
      }
    }
  }

  if (trim_endnt_p) {
    if (strand == '+') {
      if (query5 == 1 && *trimleft == 0) {
	debug(printf("query5 is 1 and trimleft is 0, so trimming 1 on left\n"));
	*trimleft = 1;
      }
      if (query3 == querylength && *trimright == 0) {
	debug(printf("query3 is %d and trimright is 0, so trimming 1 on right\n",querylength));
	*trimright = 1;
      }
    } else if (strand == '-') {
      if (query3 == querylength && *trimleft == 0) {
	debug(printf("query3 is %d and trimleft is 0, so trimming 1 on left\n",querylength));
	*trimleft = 1;
      }
      if (query5 == 1 && *trimright == 0) {
	debug(printf("query5 is 1 and trimright is 0, so trimming 1 on right\n"));
	*trimright = 1;
      }
    }
  }

  return;
}



static void
process_lines (char *refnts, int *nmatches, List_T *mismatches, int querylength, char *chromosome,
	       Pairend_T pairend, char *header, List_T lines, char *shortread, 
	       IIT_T chromosome_iit, Genome_T genome, Genomicpos_T chrstart) {
  List_T ptr, ptr2;
  char *line, *p, *q, *psave, *qsave, lastsub;
  char genomic[1024], type[1024], chr[1024], strand;
  char shortread_rc[1024], genomic_rc[1024];
  int query5, query3, querypos;
  Chrnum_T chrnum;
  Genomicpos_T pos5, pos3, pos, lastpos;
  Mismatch_T mismatch;
  int trimleft = 0, trimright = 0, ntrim;
  int readlength, len, nsubs, npos;
  int nunknowns;
  int index;
  Genomicpos_T allocoffset = 0U;
  bool processp;

  
  for (ptr = lines; ptr != NULL; ptr = List_next(ptr)) {
    nsubs = npos = 0;
    line = (char *) List_head(ptr);
    line++;			/* Skip first character */
    if (sscanf(line,"%s",genomic) != 1) {
      fprintf(stderr,"Can't parse genomic part of %s\n",line);
      abort();
    }

    /* Get query coordinates */
    p = &(line[1]);      
    while (*p != '\0' && *p != '\t') p++;
    if (*p != '\0' && *p == '\t') p++;
    if (*p == '\0') {
      fprintf(stderr,"Can't parse query coordinates part of %s\n",line);
      abort();
    }

    if (sscanf(p,"%d",&query5) != 1) {
      fprintf(stderr,"Can't parse first query coordinate in %s\n",line);
      abort();
    }
    while (*p != '\0' && isdigit(*p)) p++;
    while (*p != '\0' && !isdigit(*p)) p++;
    if (sscanf(p,"%d",&query3) != 1) {
      fprintf(stderr,"Can't parse second query coordinate in %s\n",line);
      abort();
    }

    readlength = query3 - query5 + 1;
    if (readlength < min_readlength) {
      /* Skip rest of line */
      while (*p != '\0') p++;
    } else {
      debug(printf("Got query %d to %d, with length %d\n",query5,query3,readlength));

      /* Advance to strand */
      while (*p != '\0' && *p != '\t') p++;
      if (*p != '\0' && *p == '\t') p++;
      if (*p == '\0') {
	fprintf(stderr,"Can't parse strand part of %s\n",line);
	abort();
      } else {
	strand = *p++;
      }

      /* Get chr part */
      q = p;
      len = 0;
      while (*q != '\0' && *q != ':') {
	q++;
	len++;
      }
      if (*q == '\0') {
	fprintf(stderr,"Can't parse chr part of %s\n",line);
	abort();
      } else {
	strncpy(chr,p,len);
	chr[len] = '\0';
      }

      if (desired_strand != ' ' && ((pairend == FIVE && strand != desired_strand) ||
				    (pairend == THREE && strand != opposite_strand))) {
	processp = false;

      } else if (chromosome == NULL) {
	if ((index = IIT_find_one(chromosome_iit,chr)) < 0) {
	  processp = false;
	} else {
	  allocoffset = chrstart = Interval_low(IIT_interval(chromosome_iit,index));
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
	/* Desired chromosome */
	p = ++q;
	if (sscanf(p,"%u",&pos5) != 1) {
	  fprintf(stderr,"Can't parse first chrpos in %s\n",line);
	  abort();
	}

	/* Advance past first chrpos */
	while (*p != '\0' && isdigit(*p)) p++;
	while (*p != '\0' && !isdigit(*p)) p++;
	if (sscanf(p,"%u",&pos3) != 1) {
	  fprintf(stderr,"Can't parse second chrpos in %s\n",line);
	  abort();
	}

	if (strand == '+') {
	  pos = pos5-1U;		/* Based pos3 and pos5 were 1-based */
	  p = &(shortread[query5-1]);
	  if (genome == NULL) {
	    /* Use genome according to GSNAP */
	    q = &(genomic[query5-1]);
	  } else {
	    /* Look up genomic segment.  Don't use known gene. */
	    /* printf("Before (+): %s\n",&(genomic[query5-1])); */
	    Genome_fill_buffer(&chrnum,&nunknowns,genome,chrstart + pos5 - 1U,readlength,genomic,chromosome_iit);
	    /* printf("After (+): %s\n",genomic); */
	    q = genomic;
	  }
	} else if (strand == '-') {
	  pos = pos3 - 1U;
	  make_complement_buffered(shortread_rc,&(shortread[query5-1]),readlength);
	  p = shortread_rc;
	  if (genome == NULL) {
	    /* Use genome according to GSNAP */
	    make_complement_buffered(genomic_rc,&(genomic[query5-1]),readlength);
	    q = genomic_rc;
	  } else {
	    /* Look up genomic segment.  Don't use known gene. */
	    /* printf("Before (-): %s\n",&(genomic[query5-1])); */
	    Genome_fill_buffer(&chrnum,&nunknowns,genome,chrstart + pos3 - 1U,readlength,genomic_rc,chromosome_iit);
	    /* printf("After (-): %s\n",genomic_rc); */
	    q = genomic_rc;
	  }
	} else {
	  fprintf(stderr,"Can't parse strand %c in %s\n",strand,line);
	  abort();
	}

	psave = p; qsave = q;
	debug(printf("Processing %.*s and %.*s\n",readlength,p,readlength,q));

	if (trimp == true) {
	  compute_trims(&trimleft,&trimright,p,q,readlength,strand,query5,query3,querylength);
	  debug(printf("Computed trimleft of %d and trimright of %d\n",trimleft,trimright));
	}

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
	  pos++;
	  querypos += ((strand == '+') ? +1 : -1);
	  readlength--;
	}

	while (--readlength >= trimright) {
	  debug(printf("Processing %c and %c at querypos %d, pos %u+%u\n",*p,*q,querypos,allocoffset,pos+1U));
	  if (nmatches[allocoffset+pos] == 0 && mismatches[allocoffset+pos] == NULL) {
	    refnts[allocoffset+pos] = toupper(*q);

	  } else if (allocoffset+pos == haltpos) {
	    fprintf(stderr,"Halt requested at %u+%u:\n",allocoffset,pos);
	    fprintf(stderr,"Have seen %d matches and %d types of mismatches here so far\n",
		    nmatches[allocoffset+pos],List_length(mismatches[allocoffset+pos]));
	    fprintf(stderr,"%s",header);
	    lines = List_reverse(lines);
	    for (ptr2 = lines; ptr2 != NULL; ptr2 = List_next(ptr2)) {
	      fprintf(stderr,"%s",(char *) List_head(ptr2));
	    }
	    exit(9);

	  } else if (refnts[allocoffset+pos] != toupper(*q)) {
	    fprintf(stderr,"Two different genomic chars %c and %c at position %u+%u\n",
		    refnts[allocoffset+pos],*q,allocoffset,pos+1U);
	    fprintf(stderr,"Have seen %d matches and %d types of mismatches here so far\n",
		    nmatches[allocoffset+pos],List_length(mismatches[allocoffset+pos]));
	    fprintf(stderr,"Most recent lines were:\n");
	    fprintf(stderr,"%s",header);
	    lines = List_reverse(lines);
	    for (ptr2 = lines; ptr2 != NULL; ptr2 = List_next(ptr2)) {
	      fprintf(stderr,"%s",(char *) List_head(ptr2));
	    }
	    fprintf(stderr,"To halt at first occurrence, re-run with -h <coord>\n");
	    exit(9);

	  }

	  if (dibasep == true) {
	    if (toupper(*q) == refnts[allocoffset+pos]) {
	      nmatches[allocoffset+pos] += 1;
	    } else if ((mismatch = find_mismatch(mismatches[allocoffset+pos],toupper(*q),querypos)) != NULL) {
	      mismatch->count += 1;
	    } else {
	      mismatches[allocoffset+pos] = List_push(mismatches[allocoffset+pos],(void *) Mismatch_new(toupper(*q),querypos));
	    }

	  } else if (mode == ALL_SCORES) {
	    if (toupper(*p) == toupper(*q)) {
	      /* Count matches, even if query quality score was low */
	      debug(printf("Incrementing nmatches at pos %u\n",pos+1U));
	      nmatches[allocoffset+pos] += 1;
	    } else if ((mismatch = find_mismatch(mismatches[allocoffset+pos],toupper(*p),querypos)) != NULL) {
	      mismatch->count += 1;
	    } else {
	      mismatches[allocoffset+pos] = List_push(mismatches[allocoffset+pos],(void *) Mismatch_new(toupper(*p),querypos));
	    }
	  } else if (mode == HIGH_SCORES) {
	    if (islower(*p)) {
	      /* Skip.  Don't count anything where query quality score was low */
	    } else if (*p == toupper(*q)) {
	      debug(printf("Incrementing nmatches at pos %u+%u\n",allocoffset,pos+1U));
	      nmatches[allocoffset+pos] += 1;
	    } else if ((mismatch = find_mismatch(mismatches[allocoffset+pos],*p,querypos)) != NULL) {
	      debug(printf("Incrementing mismatch at pos %u+%u\n",allocoffset,pos+1U));
	      mismatch->count += 1;
	    } else {
	      debug(printf("Adding new mismatch at pos %u+%u\n",allocoffset,pos+1U));
	      mismatches[allocoffset+pos] = List_push(mismatches[allocoffset+pos],(void *) Mismatch_new(*p,querypos));
	    }
	  } else if (mode == MISMATCH_HIGH_SCORES) {
	    /* Counts all matches, regardless of quality score, but only mismatches with high quality score */
	    if (toupper(*p) == toupper(*q)) {
	      /* Count matches, even if query quality score was low */
	      debug(printf("Incrementing nmatches at pos %u+%u\n",allocoffset,pos+1U));
	      nmatches[allocoffset+pos] += 1;
	    } else if (islower(*p)) {
	      /* Skip.  Don't count a mismatch if query quality score was low */
	    } else if ((mismatch = find_mismatch(mismatches[allocoffset+pos],*p,querypos)) != NULL) {
	      mismatch->count += 1;
	    } else {
	      mismatches[allocoffset+pos] = List_push(mismatches[allocoffset+pos],(void *) Mismatch_new(*p,querypos));
	    }
	  } else if (mode == ONE_SUB || mode == ONE_SUB_HIGH_SCORES) {
	    /* Counts only sequences where zero or one subs relative to genome */
	    if (toupper(*p) != toupper(*q)) {
	      nsubs++;
	      lastpos = pos;
	      lastsub = toupper(*p);
	    }
	    npos++;
	  }

	  p++;
	  q++;
	  pos++;
	  querypos += ((strand == '+') ? +1 : -1);
	}

	debug({
	    while (readlength-- >= 0) {
	      printf("Trim right: Skipping %c and %c at querypos %d, pos %u\n",*p,*q,querypos,pos+1U);
	      p++;
	      q++;
	      pos++;
	      querypos += ((strand == '+') ? +1 : -1);
	    }
	  });

      }
    }
  }

  return;
}


static void
print_chromosome_blocks (Genomicpos_T allocoffset, Genomicpos_T chrlength,
			 char *refnts, List_T *mismatches, int *nmatches,
			 IIT_T chromosome_iit, char *chr) {
  int segno;
  int n, i;
  Genomicpos_T pos, posi, firsti, lasti;
  bool hitp;
  Mismatch_T mismatch, mismatch0, *array;
  List_T unique_mismatches, ptr;
 
  segno = 0;
  for (pos = 0; pos < chrlength; pos += blocksize) {
    hitp = false;
    for (posi = pos; posi < chrlength && posi < pos + blocksize; posi += 1U) {
      if (nmatches[allocoffset+posi] > 0 || mismatches[allocoffset+posi] != NULL) {
	if (hitp == false) {
	  firsti = posi;
	  hitp = true;
	}
	lasti = posi;
      }
    }
      
    if (hitp == true) {
      /* old: printf(">%s_%d %u %u %s\n",chromosome,++segno,firsti+1U,lasti+1U,chromosome); */
      printf(">%s_%d %s:%u..%u\n",chr,++segno,chr,firsti+1U,lasti+1U);
      for (posi = firsti; posi <= lasti; posi++) {
	if (nmatches[allocoffset+posi] > 0 || mismatches[allocoffset+posi] != NULL) {
	  printf("%c%d",refnts[allocoffset+posi],nmatches[allocoffset+posi]);
	  if (mismatches[allocoffset+posi] != NULL) {
	    if ((n = List_length(mismatches[allocoffset+posi])) == 1) {
	      mismatch = (Mismatch_T) List_head(mismatches[allocoffset+posi]);
	      printf(" %c%d/1",mismatch->nt,mismatch->count);
	    } else {
	      unique_mismatches = NULL;
	      for (ptr = mismatches[allocoffset+posi]; ptr != NULL; ptr = List_next(ptr)) {
		mismatch0 = List_head(ptr);
		if ((mismatch = find_mismatch_nt(unique_mismatches,mismatch0->nt)) != NULL) {
		  mismatch->count += mismatch0->count;
		  mismatch->querypos += 1; /* Using querypos here as number of unique querypos's */
		} else {
		  mismatch = Mismatch_new(mismatch0->nt,1);
		  mismatch->count = mismatch0->count;
		  unique_mismatches = List_push(unique_mismatches,mismatch);
		}
	      }

	      array = (Mismatch_T *) List_to_array(unique_mismatches,NULL);
	      qsort(array,List_length(unique_mismatches),sizeof(Mismatch_T),Mismatch_cmp);

	      for (i = 0; i < List_length(unique_mismatches); i++) {
		mismatch = array[i];
		printf(" %c%d/%d",mismatch->nt,mismatch->count,mismatch->querypos);
	      }
	      FREE(array);

	      for (ptr = unique_mismatches; ptr != NULL; ptr = List_next(ptr)) {
		mismatch = List_head(ptr);
		Mismatch_free(&mismatch);
	      }
	      List_free(&unique_mismatches);

	    }
	  }
	}
	printf("\n");
      }
    }
  }

  return;
}



static void
lines_gc (List_T *lines) {
  char *line;
  void *item;
  
  while (*lines != NULL) {
    *lines = List_pop(*lines,&item);
    line = (char *) item;
    FREE(line);
  }
  return;
}




int
main (int argc, char *argv[]) {
  char line[1024000], *copy, *header = NULL;
  char *genomesubdir = NULL, *fileroot = NULL, *iitfile;
  IIT_T chromosome_iit;
  int querylength;
  int index;
  Pairend_T pairend = FIVE;

  char shortread[1000];
  List_T lines;
  Genomicpos_T allocoffset = 0U, alloclength, chrstart = 0U, chrlength;
  char *chr;
  bool allocp;

  char *refnts;
  int *nmatches, nhits = 0;
  List_T *mismatches;

  Genome_T genome = NULL;

  int opt;
  extern int optind;
  extern char *optarg;
  int long_option_index = 0;

  while ((opt = getopt_long(argc,argv,
			    "D:d:2c:FRP:Mu:b:a:h:T:t:l:V?",
			    long_options, &long_option_index)) != -1) {
    switch (opt) {
    case 'D': user_genomedir = optarg; break;
    case 'd': 
      dbroot = (char *) CALLOC(strlen(optarg)+1,sizeof(char));
      strcpy(dbroot,optarg);
      break;
    case '2': dibasep = true; break;

    case 'c': chromosome = optarg; break;
    case 'F': desired_strand = '+'; opposite_strand = '-'; break;
    case 'R': desired_strand = '-'; opposite_strand = '+'; break;
    case 'P':
      switch (atoi(optarg)) {
      case 1: want_second_p = false; break;
      case 2: want_first_p = false; break;
      default: fprintf(stderr,"Flag -P mode %s not recognized.\n",optarg); exit(9);
      }
      break;

    case 'M': remappedp = true; break;

    case 'u':
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
      } else if (!strcmp(optarg,"high")) {
	mode = HIGH_SCORES;
      } else if (!strcmp(optarg,"mismatch_high")) {
	mode = MISMATCH_HIGH_SCORES;
      } else if (!strcmp(optarg,"onesub")) {
	mode = ONE_SUB;
      } else if (!strcmp(optarg,"onesub_high")) {
	mode = ONE_SUB_HIGH_SCORES;
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

    case 't': 
      switch (atoi(optarg)) {
      case 0: trim_endnt_p = false; break;
      case 1: trim_endnt_p = true; break;
      default: fprintf(stderr,"Trimendnt mode %s not recognized.\n",optarg); exit(9);
      }
      break;

    case 'l': min_readlength = atoi(optarg); break;

    case 'V': print_program_version(); exit(0);
    case '?': print_program_usage(); exit(0);
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

  genomesubdir = Datadir_find_genomesubdir(&fileroot,&dbversion,user_genomedir,dbroot);

  iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			    strlen(fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
  sprintf(iitfile,"%s/%s.chromosome.iit",genomesubdir,fileroot);
  chromosome_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
			    /*divstring*/NULL,/*add_iit_p*/false,/*labels_read_p*/true);
  FREE(iitfile);

  if (remappedp == true || dibasep == true) {
    /* Need genome to determine wild-type, because "known gene" may not match reference genome */
    genome = Genome_new(genomesubdir,fileroot,/*snps_root*/NULL,/*uncompressedp*/false,
			/*batchp*/true);
  }

  if (dibasep == true && trimp == true) {
    fprintf(stderr,"Turning trimming off for 2-base encoded alignments\n");
    trimp = false;
  }

  if (chromosome == NULL) {
    alloclength = IIT_totallength(chromosome_iit);
    fprintf(stderr,"Starting to allocate memory for %u positions\n",alloclength);

    refnts = (char *) CALLOC(alloclength,sizeof(char));
    nmatches = (int *) CALLOC(alloclength,sizeof(int));
    mismatches = (List_T *) CALLOC(alloclength,sizeof(List_T));

  } else if ((index = IIT_find_one(chromosome_iit,chromosome)) < 0) {
    fprintf(stderr,"Cannot find chromosome %s in chromosome IIT file\n",chromosome);
    exit(9);

  } else {
    chrstart = Interval_low(IIT_interval(chromosome_iit,index));
    alloclength = Interval_length(IIT_interval(chromosome_iit,index));
    fprintf(stderr,"Starting to allocate memory for %u positions\n",alloclength);

    refnts = (char *) CALLOC(alloclength,sizeof(char));
    nmatches = (int *) CALLOC(alloclength,sizeof(int));
    mismatches = (List_T *) CALLOC(alloclength,sizeof(List_T));
  }

  fprintf(stderr,"Done allocating memory\n");

  while (fgets(line,1024000,stdin) != NULL) {
    if (line[0] == '>') {
      pairend = FIVE;
      sscanf(&(line[1]),"%s",shortread);
      querylength = strlen(shortread);
      header = (char *) CALLOC(strlen(line)+1,sizeof(char));
      strcpy(header,line);
      lines = (List_T) NULL;
      nhits = 0;

    } else if (line[0] == '<') {
      pairend = THREE;
      sscanf(&(line[1]),"%s",shortread);
      querylength = strlen(shortread);
      header = (char *) CALLOC(strlen(line)+1,sizeof(char));
      strcpy(header,line);
      lines = (List_T) NULL;
      nhits = 0;

    } else if (line[0] == '\n') {
      if (uniquep == false || nhits == 1) {
	if ((pairend == FIVE && want_first_p == true) ||
	    (pairend == THREE && want_second_p == true)) {
	  process_lines(refnts,nmatches,mismatches,querylength,chromosome,
			pairend,header,lines,shortread,chromosome_iit,genome,chrstart);
	}
      }
      lines_gc(&lines);
      FREE(header);

    } else {
      if (line[0] != ',') {
	nhits++;
      }
      copy = (char *) CALLOC(strlen(line)+1,sizeof(char));
      strcpy(copy,line);
      debug(printf("Pushing %s",copy));
      lines = List_push(lines,(void *) copy);
    }
  }

  if (header != NULL) {
    if (uniquep == false || nhits == 1) {
      if ((pairend == FIVE && want_first_p == true) ||
	  (pairend == THREE && want_second_p == true)) {
	process_lines(refnts,nmatches,mismatches,querylength,chromosome,
		      pairend,header,lines,shortread,chromosome_iit,genome,chrstart);
      }
    }
    lines_gc(&lines);
    FREE(header);
  }


  /* Print results */
  if (chromosome != NULL) {
    print_chromosome_blocks(/*allocoffset*/0,/*chrlength*/alloclength,refnts,mismatches,nmatches,chromosome_iit,/*chr*/chromosome);
  } else {
    for (index = 0; index < IIT_total_nintervals(chromosome_iit); index++) {
      allocoffset = Interval_low(IIT_interval(chromosome_iit,index));
      chrlength = Interval_length(IIT_interval(chromosome_iit,index));
      chr = IIT_label(chromosome_iit,index,&allocp);
      print_chromosome_blocks(allocoffset,chrlength,refnts,mismatches,nmatches,chromosome_iit,chr);
      if (allocp == true) {
	FREE(chr);
      }
    }
  }

  if (genome != NULL) {
    Genome_free(&genome);
  }

  return 0;
}

