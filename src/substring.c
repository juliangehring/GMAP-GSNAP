static char rcsid[] = "$Id: substring.c 37257 2011-03-28 16:39:14Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "substring.h"
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "assert.h"
#include "mem.h"
#include "stage1hr.h"		/* For MAX_QUERYLENGTH */
#include "maxent_hr.h"
#include "listdef.h"
#include "list.h"
#include "complement.h"
#include "genome_hr.h"
#include "mapq.h"


#define TRIM_MATCH_SCORE 1

#define CONSISTENT_TEXT "consistent"
#define TRANSLOCATION_TEXT "translocation"
#define INVERSION_TEXT "inversion"
#define SCRAMBLE_TEXT "scramble"


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* mark_mismatches */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

/* Substring_new */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

/* Substring_overlap_p and Substring_insert_length */
#ifdef DEBUG3
#define debug3(x) x
#else
#define debug3(x)
#endif


/* splice site probs */
#ifdef DEBUG4
#define debug4(x) x
#else
#define debug4(x)
#endif


/* trimming */
#ifdef DEBUG8
#define debug8(x) x
#else
#define debug8(x)
#endif



static bool print_nsnpdiffs_p;
static bool print_ncolordiffs_p;
static bool print_snplabels_p;
static bool show_refdiff_p;

static IIT_T snps_iit;
static int *snps_divint_crosstable;

static IIT_T splicesites_iit;
static int *splicesites_divint_crosstable;

static int donor_typeint;
static int acceptor_typeint;

static int trim_mismatch_score;
static bool output_sam_p;


void
Substring_setup (bool print_ncolordiffs_p_in, bool print_nsnpdiffs_p_in, bool print_snplabels_p_in,
		 bool show_refdiff_p_in, IIT_T snps_iit_in, int *snps_divint_crosstable_in,
		 IIT_T splicesites_iit_in, int *splicesites_divint_crosstable_in,
		 int donor_typeint_in, int acceptor_typeint_in, int trim_mismatch_score_in,
		 bool output_sam_p_in) {
  print_ncolordiffs_p = print_ncolordiffs_p_in;
  print_nsnpdiffs_p = print_nsnpdiffs_p_in;
  print_snplabels_p = print_snplabels_p_in;
  show_refdiff_p = show_refdiff_p_in;

  snps_iit = snps_iit_in;
  snps_divint_crosstable = snps_divint_crosstable_in;

  splicesites_iit = splicesites_iit_in;
  splicesites_divint_crosstable = splicesites_divint_crosstable_in;

  donor_typeint = donor_typeint_in;
  acceptor_typeint = acceptor_typeint_in;

  trim_mismatch_score = trim_mismatch_score_in;
  output_sam_p = output_sam_p_in;

  return;
}


static char *
endtype_string (Endtype_T endtype) {
  switch (endtype) {
  case END: return "end";
  case INS: return "ins";
  case DEL: return "del";
  case DON: return "don";
  case ACC: return "acc";
  case AMB_DON: return "amb_don";
  case AMB_ACC: return "amb_acc";
  case TERM: return "term";
  default:
    fprintf(stderr,"Unexpected endtype %d\n",endtype);
    abort();
  }
  return "";
}


static char complCode[128] = COMPLEMENT_LC;

static char *
make_complement_buffered (char *complement, char *sequence, unsigned int length) {
  int i, j;

  /* complement = (char *) CALLOC(length+1,sizeof(char)); */
  for (i = length-1, j = 0; i >= 0; i--, j++) {
    complement[j] = complCode[(int) sequence[i]];
  }
  complement[length] = '\0';
  return complement;
}


#define T Substring_T

struct T {
  bool exactp;
  int nmismatches_whole;        /* Over total substring.  Comes from calling procedure */
  int nmismatches_bothdiff;	/* Over region left after trimming */
  int nmismatches_refdiff;	/* Over region left after trimming */
  /* nsnpdiffs = nmismatches_bothdiff - nmismatches_refdiff */
  int nmatches;			/* Over region left after trimming */

  int ncolordiffs;
  int trim_left;
  int trim_right;
  Chrnum_T chrnum;
  Genomicpos_T chroffset;

  Genomicpos_T left_genomicseg;	/* left needed to retrieve genomicseg */
  Genomicpos_T left;	 /* adjusted by Substring_new for aligndiff */
  Genomicpos_T genomicstart;	/* For region corresponding to querylength */
  Genomicpos_T genomicend;

  Endtype_T start_endtype;
  Endtype_T end_endtype;

  int querystart_orig;
  int queryend_orig;
  int querystart;		/* For part that aligns to genome */
  int queryend;
  int querylength;

  Genomicpos_T alignstart;	/* For part that aligns to genome, including part that is trimmed */
  Genomicpos_T alignend;

  Genomicpos_T alignstart_trim;	/* For part that aligns to genome, excluding part that is trimmed */
  Genomicpos_T alignend_trim;

  int alignoffset;
  int extraleft;
  int extraright;

  int genomiclength;
  bool plusp;

  char *genomic_bothdiff; /* In same direction as query.  NULL if same
  			     as query.  Has dashes outside of
  			     querystart..(queryend-1) for indels and
  			     splices.  Has lowercase to indicate
  			     trimmed regions, terminal regions,
  			     deletions, splice dinucleotides, and
  			     mismatches from ref and alt.  Use for
  			     MAPQ computations and scoring. */

  char *genomic_refdiff;  /* Same as above, but lowercase for
			     mismatches from ref only.  For
			     non-SNP-tolerant alignment, this is just
			     a pointer to genomic_bothdiff.  Use for
			     NM and MD computations.  */

  double mapq_loglik;

  /* for splices */
  bool chimera_sensep;

  int splicesites_i;		/* Keep this for finding short-end
				   splicing later.  If < 0, then it's
				   a novel splice site. */
  bool chimera_knownp;
  Genomicpos_T chimera_modelpos;
  int chimera_pos;
  double chimera_prob;

  /* for shortexon (always use *_1 for acceptor and *_2 for donor) */
  /* for donor/acceptor: the ambiguous site */
  int splicesites_i_2;
  bool chimera_knownp_2;
  Genomicpos_T chimera_modelpos_2;
  int chimera_pos_2;
  double chimera_prob_2;
};


static void
fill_w_dashes (char *string, int start, int end) {
  int i;

  for (i = start; i < end; i++) {
    string[i] = '-';
  }
  return;
}


static bool
dibase_mismatch_p (char leftnt, char leftcolor, char nt, char rightcolor, char rightnt) {
  char left_predict, right_predict;

  switch (toupper(leftnt)) {
  case 'A':
    switch (leftcolor) {
    case '0': left_predict = 'A'; break;
    case '1': left_predict = 'C'; break;
    case '2': left_predict = 'G'; break;
    case '3': left_predict = 'T'; break;
    }
    break;
  case 'C':
    switch (leftcolor) {
    case '0': left_predict = 'C'; break;
    case '1': left_predict = 'A'; break;
    case '2': left_predict = 'T'; break;
    case '3': left_predict = 'G'; break;
    }
    break;
  case 'G':
    switch (leftcolor) {
    case '0': left_predict = 'G'; break;
    case '1': left_predict = 'T'; break;
    case '2': left_predict = 'A'; break;
    case '3': left_predict = 'C'; break;
    }
    break;
  case 'T':
    switch (leftcolor) {
    case '0': left_predict = 'T'; break;
    case '1': left_predict = 'G'; break;
    case '2': left_predict = 'C'; break;
    case '3': left_predict = 'A'; break;
    }
    break;
  }

  switch (toupper(rightnt)) {
  case 'A':
    switch (rightcolor) {
    case '0': right_predict = 'A'; break;
    case '1': right_predict = 'C'; break;
    case '2': right_predict = 'G'; break;
    case '3': right_predict = 'T'; break;
    }
    break;
  case 'C':
    switch (rightcolor) {
    case '0': right_predict = 'C'; break;
    case '1': right_predict = 'A'; break;
    case '2': right_predict = 'T'; break;
    case '3': right_predict = 'G'; break;
    }
    break;
  case 'G':
    switch (rightcolor) {
    case '0': right_predict = 'G'; break;
    case '1': right_predict = 'T'; break;
    case '2': right_predict = 'A'; break;
    case '3': right_predict = 'C'; break;
    }
    break;
  case 'T':
    switch (rightcolor) {
    case '0': right_predict = 'T'; break;
    case '1': right_predict = 'G'; break;
    case '2': right_predict = 'C'; break;
    case '3': right_predict = 'A'; break;
    }
    break;
  }

  if (left_predict != right_predict) {
    /* printf("left_predict %c != right_predict %c\n",left_predict,right_predict); */
    return false;
  } else if (left_predict == toupper(nt)) {
    return false;
  } else {
    return true;
  }
}


/* query is numbers */
static void
mark_mismatches_dibase (char *genome, char *query, int start, int end, bool plusp) {
  int i;

  debug1(printf("query:  %s\n",query));
  debug1(printf("genome: %s\n",genome));
  debug1(printf("mark:   "));

  if (plusp) {
    /* querypos start */
    debug1(printf("*"));
    if (dibase_mismatch_p(/*leftnt*/genome[start],/*leftcolor*/'0',genome[start],query[start],genome[start+1])) {
      genome[start] = (char) tolower(genome[start]);
    }

    for (i = start+1; i < end-1; i++) {
      debug1(printf("*"));

      if (dibase_mismatch_p(genome[i-1],query[i-1],genome[i],query[i],genome[i+1])) {
	/* printf("Got a dibase mismatch at %d with %c %c %c %c %c\n",
	   i,genome[i-1],query[i-1],genome[i],query[i],genome[i+1]); */
	genome[i] = (char) tolower(genome[i]);
      }
    }

    /* querypos (end-1) */
    if (dibase_mismatch_p(genome[end-2],query[end-2],genome[end-1],
			  /*rightcolor*/'0',/*rightnt*/genome[end-1])) {
      genome[end-1] = (char) tolower(genome[end-1]);
    }

  } else {
    /* querypos start */
    debug1(printf("*"));
    if (dibase_mismatch_p(/*leftnt*/genome[start],/*leftcolor*/'0',genome[start],query[start+1],genome[start+1])) {
      genome[start] = (char) tolower(genome[start]);
    }

    for (i = start+1; i < end-1; i++) {
      debug1(printf("*"));

      if (dibase_mismatch_p(genome[i-1],query[i],genome[i],query[i+1],genome[i+1])) {
	/* printf("Got a dibase mismatch at %d with %c %c %c %c %c\n",
	   i,genome[i-1],query[i],genome[i],query[i+1],genome[i+1]); */
	genome[i] = (char) tolower(genome[i]);
      }
    }

    /* querypos (end-1) */
    if (dibase_mismatch_p(genome[end-2],query[end-1],genome[end-1],
			  /*rightcolor*/'0',/*rightnt*/genome[end-1])) {
      genome[end-1] = (char) tolower(genome[end-1]);
    }
  }

  debug1(printf("\n"));

  return;
}


static int
count_mismatches_dibase (char *genome, char *query, int start, int end, bool plusp) {
  int nmismatches = 0, i;

  debug1(printf("query:  %s\n",query));
  debug1(printf("genome: %s\n",genome));
  debug1(printf("mark:   "));

  if (plusp) {
    /* querypos start */
    debug1(printf("*"));
    if (dibase_mismatch_p(/*leftnt*/genome[start],/*leftcolor*/'0',genome[start],query[start],genome[start+1])) {
      nmismatches++;
    }

    for (i = start+1; i < end-1; i++) {
      debug1(printf("*"));

      if (dibase_mismatch_p(genome[i-1],query[i-1],genome[i],query[i],genome[i+1])) {
	/* printf("Got a dibase mismatch at %d with %c %c %c %c %c\n",
	   i,genome[i-1],query[i-1],genome[i],query[i],genome[i+1]); */
	nmismatches++;
      }
    }

    /* querypos (end-1) */
    if (dibase_mismatch_p(genome[end-2],query[end-2],genome[end-1],
			  /*rightcolor*/'0',/*rightnt*/genome[end-1])) {
      nmismatches++;
    }

  } else {
    /* querypos start */
    debug1(printf("*"));
    if (dibase_mismatch_p(/*leftnt*/genome[start],/*leftcolor*/'0',genome[start],query[start+1],genome[start+1])) {
      nmismatches++;
    }

    for (i = start+1; i < end-1; i++) {
      debug1(printf("*"));

      if (dibase_mismatch_p(genome[i-1],query[i],genome[i],query[i+1],genome[i+1])) {
	/* printf("Got a dibase mismatch at %d with %c %c %c %c %c\n",
	   i,genome[i-1],query[i],genome[i],query[i+1],genome[i+1]); */
	nmismatches++;
      }
    }

    /* querypos (end-1) */
    if (dibase_mismatch_p(genome[end-2],query[end-1],genome[end-1],
			  /*rightcolor*/'0',/*rightnt*/genome[end-1])) {
      nmismatches++;
    }
  }

  debug1(printf("\n"));
  debug1(printf("%d mismatches\n\n",nmismatches));

  return nmismatches;
}


static int
trim_left_end (char *query, Compress_T query_compress, UINT4 *genome_blocks, UINT4 *snp_blocks,
	       Genomicpos_T left, int querystart, int queryend,
	       int querylength, bool dibasep, bool cmetp, bool plusp) {
  int bestscore, score;
  int trim5, alignlength, pos, prevpos, i;

  int mismatch_positions[MAX_QUERYLENGTH];
  int nmismatches, ncolordiffs;

  debug8(printf("Entered trim_left_end with querystart %d, queryend %d\n",querystart,queryend));

  alignlength = queryend - querystart;
  bestscore = 0;
  score = 0;
  trim5 = 0;

  if (plusp == true) {
    nmismatches = Genome_mismatches_right(mismatch_positions,&ncolordiffs,/*max_mismatches*/alignlength,
					  query,query_compress,genome_blocks,snp_blocks,
					  left,/*pos5*/querystart,/*pos3*/queryend,dibasep,cmetp,plusp);
    debug8(printf("%d mismatches:",nmismatches));
    debug8(
	   for (i = 0; i < nmismatches; i++) {
	     printf(" %d",mismatch_positions[i]);
	   }
	   printf("\n");
	   );

    prevpos = queryend;
    for (i = 0; i < nmismatches; i++) {
      pos = mismatch_positions[i] + 1; /* the position just after the mismatch */
      score += (prevpos - pos)*TRIM_MATCH_SCORE;
      if (score > bestscore) {
	bestscore = score;
	trim5 = pos;
      }
      debug8(printf("Trim left pos %d, score %d, trim5 %d\n",pos,score,trim5));

      score += trim_mismatch_score; /* For the mismatch */
      if (score < 0) {
	score = 0;
      }
      prevpos = pos - 1;		/* On the mismatch */
    }

  } else {
    nmismatches = Genome_mismatches_left(mismatch_positions,&ncolordiffs,/*max_mismatches*/alignlength,
					 query,query_compress,genome_blocks,snp_blocks,
					 left,/*pos5*/querylength - queryend,/*pos3*/querylength - querystart,
					 dibasep,cmetp,plusp);

    debug8(printf("%d mismatches:",nmismatches));
    debug8(
	   for (i = 0; i < nmismatches; i++) {
	     printf(" %d",querylength - 1 - mismatch_positions[i]);
	   }
	   printf("\n");
	   );

    prevpos = queryend;
    for (i = 0; i < nmismatches; i++) {
      pos = (querylength - 1 - mismatch_positions[i]) + 1; /* the position just after the mismatch */
      score += (prevpos - pos)*TRIM_MATCH_SCORE;
      if (score > bestscore) {
	bestscore = score;
	trim5 = pos;
      }
      debug8(printf("Trim left pos %d, score %d, trim5 %d\n",pos,score,trim5));

      score += trim_mismatch_score; /* For the mismatch */
      if (score < 0) {
	score = 0;
      }
      prevpos = pos - 1;		/* On the mismatch */
    }
  }

  score += prevpos*TRIM_MATCH_SCORE;
  if (score > bestscore) {
    bestscore = score;
    trim5 = 0;
  }

  debug8(printf("Trim left pos 0, score %d, trim5 %d\n",score,trim5));
  debug8(printf("\n"));

  return trim5;
}



/* querystart (typically 0) and queryend (typically querylength) are exclusive */
/* sequences may have had lower case characters marked */
static int
trim_right_end (char *query, Compress_T query_compress, UINT4 *genome_blocks, UINT4 *snp_blocks,
		Genomicpos_T left, int querystart, int queryend,
		int querylength, bool dibasep, bool cmetp, bool plusp) {
  int bestscore, score;
  int trim3, alignlength, pos, prevpos, i;

  int mismatch_positions[MAX_QUERYLENGTH];
  int nmismatches, ncolordiffs;


  debug8(printf("Entered trim_right_end with querystart %d, queryend %d\n",querystart,queryend));

  alignlength = queryend - querystart;
  bestscore = 0;
  score = 0;
  trim3 = 0;

  if (plusp == true) {
    nmismatches = Genome_mismatches_left(mismatch_positions,&ncolordiffs,/*max_mismatches*/alignlength,
					 query,query_compress,genome_blocks,snp_blocks,
					 left,/*pos5*/querystart,/*pos3*/queryend,dibasep,cmetp,plusp);

    debug8(printf("%d mismatches:",nmismatches));
    debug8(
	   for (i = 0; i < nmismatches; i++) {
	     printf(" %d",mismatch_positions[i]);
	   }
	   printf("\n");
	   );

    prevpos = querystart - 1;
    for (i = 0; i < nmismatches; i++) {
      pos = mismatch_positions[i] - 1; /* the position just before the mismatch */
      score += (pos - prevpos)*TRIM_MATCH_SCORE;
      if (score > bestscore) {
	bestscore = score;
	trim3 = querylength - pos - 1;
      }
      debug8(printf("Trim right pos %d, score %d, trim3 %d\n",pos,score,trim3));
      
      score += trim_mismatch_score; /* For the mismatch */
      if (score < 0) {
	score = 0;
      }
      prevpos = pos + 1;		/* On the mismatch */
    }

  } else {
    nmismatches = Genome_mismatches_right(mismatch_positions,&ncolordiffs,/*max_mismatches*/alignlength,
					  query,query_compress,genome_blocks,snp_blocks,
					  left,/*pos5*/querylength - queryend,/*pos3*/querylength - querystart,
					  dibasep,cmetp,plusp);

    debug8(printf("%d mismatches:",nmismatches));
    debug8(
	   for (i = 0; i < nmismatches; i++) {
	     printf(" %d",querylength - 1 - mismatch_positions[i]);
	   }
	   printf("\n");
	   );

    prevpos = querystart - 1;
    for (i = 0; i < nmismatches; i++) {
      pos = (querylength - 1 - mismatch_positions[i]) - 1; /* the position just before the mismatch */
      score += (pos - prevpos)*TRIM_MATCH_SCORE;
      if (score > bestscore) {
	bestscore = score;
	trim3 = querylength - pos - 1;
      }
      debug8(printf("Trim right pos %d, score %d, trim3 %d\n",pos,score,trim3));
      
      score += trim_mismatch_score; /* For the mismatch */
      if (score < 0) {
	score = 0;
      }
      prevpos = pos + 1;		/* On the mismatch */
    }
  }

  score += (queryend - 1 - prevpos)*TRIM_MATCH_SCORE;
  if (score > bestscore) {
    bestscore = score;
    trim3 = 0;
  }

  debug8(printf("Trim right pos %d, score %d, trim3 %d\n",queryend-1,score,trim3));
  debug8(printf("\n"));

  return trim3;
}





void
Substring_free (T *old) {
  if ((*old)->genomic_bothdiff != NULL) {
    if ((*old)->genomic_refdiff != (*old)->genomic_bothdiff) {
      FREE((*old)->genomic_refdiff);
    }
    FREE((*old)->genomic_bothdiff);
  }

  FREE(*old);
  return;
}


bool
Substring_overlap_p (T substring1, T substring2) {
  Genomicpos_T low1, high1, low2, high2;

  if (substring1->plusp == true) {
    low1 = substring1->alignstart;
    high1 = substring1->alignend;
  } else {
    low1 = substring1->alignend;
    high1 = substring1->alignstart;
  }

  if (substring2->plusp == true) {
    low2 = substring2->alignstart;
    high2 = substring2->alignend;
  } else {
    low2 = substring2->alignend;
    high2 = substring2->alignstart;
  }

  debug3(printf("Checking overlap between %u..%u and %u..%u",low1,high1,low2,high2));

  if (high2 < low1) {
    debug3(printf(" => no because %u < %u\n",high2,low1));
    return false;
  } else if (low2 > high1) {
    debug3(printf(" => no because %u > %u\n",low2,high1));
    return false;
  } else {
    debug3(printf(" => yes\n"));
    return true;
  }
}


Genomicpos_T
Substring_insert_length (T substring5, T substring3) {
  Genomicpos_T pos5, pos3;

  pos5 = substring5->genomicstart;
  debug3(printf("pos5 %u\n",substring5->genomicstart));

  pos3 = substring3->genomicend;
  debug3(printf("pos3 %u\n",substring3->genomicend));

  if (pos5 > pos3) {
    return pos5 - pos3;
  } else {
    return pos3 - pos5;
  }
}



static void
mark_mismatches_met_gsnap (char *gbuffer, char *query, int start, int end, bool plusp) {
  int i;
  
  debug1(printf("query:  %s\n",query));
  debug1(printf("genome: %s\n",gbuffer));
  debug1(printf("count:  "));

#if 0
  if (plusp) {
#endif
    for (i = start; i < end; i++) {
      if (gbuffer[i] == 'C' && query[i] == 'T') {
	debug1(printf("."));
	gbuffer[i] = '.';
      } else if (query[i] != gbuffer[i]) {
	debug1(printf("x"));
	assert(gbuffer[i] != OUTOFBOUNDS);
	gbuffer[i] = (char) tolower(gbuffer[i]);
      } else {
	debug1(printf("*"));
      }
    }

#if 0
  } else {
    for (i = start; i < end; i++) {
      if (gbuffer[i] == 'G' && query[i] == 'A') {
	debug1(printf("."));
	gbuffer[i] = '.';
      } else if (query[i] != gbuffer[i]) {
	debug1(printf("x"));
	assert(gbuffer[i] != OUTOFBOUNDS);
	gbuffer[i] = (char) tolower(gbuffer[i]);
      } else {
	debug1(printf("*"));
      }
    }
  }
#endif

  debug1(printf("\n"));

  debug1(printf("query:  %s\n",query));
  debug1(printf("genome: %s\n",gbuffer));
  debug1(printf("\n"));

  return;
}

static void
mark_mismatches_met_sam (char *gbuffer, char *query, int start, int end) {
  int i;
  
  debug1(printf("query:  %s\n",query));
  debug1(printf("genome: %s\n",gbuffer));
  debug1(printf("count:  "));

  for (i = start; i < end; i++) {
    if (gbuffer[i] == 'C' && query[i] == 'T') {
      debug1(printf("."));
      gbuffer[i] = 'T';		/* Avoids showing mismatches in MD and NM strings */
    } else if (query[i] != gbuffer[i]) {
      debug1(printf("x"));
      assert(gbuffer[i] != OUTOFBOUNDS);
      gbuffer[i] = (char) tolower(gbuffer[i]);
    } else {
      debug1(printf("*"));
    }
  }

  debug1(printf("\n"));

  debug1(printf("query:  %s\n",query));
  debug1(printf("genome: %s\n",gbuffer));
  debug1(printf("\n"));

  return;
}



static char *
embellish_genomic (char *genomic_diff, char *query, int querystart, int queryend, int querylength,
		   int alignoffset, int extraleft, int extraright, bool cmetp, bool plusp) {
  char *result;
  int i, j, k;

  result = (char *) CALLOC(querylength+1,sizeof(char));
  result[querylength] = '\0';

  /* Add aligned region with lower-case diffs, surrounded by dashes */
  fill_w_dashes(result,0,querystart);
  debug2(printf("g1: %s\n",result));

  strncpy(&(result[querystart]),&(genomic_diff[alignoffset]),queryend-querystart);
  debug2(printf("g1: %s\n",result));

  if (cmetp == true) {
    mark_mismatches_met_gsnap(result,query,querystart,queryend,plusp);
  }

  fill_w_dashes(result,queryend,querylength);
  debug2(printf("g1: %s\n",result));

  /* Add terminal ends as lower-case */
  for (k = 0, i = querystart-1, j = alignoffset-1; k < extraleft && i >= 0; k++, i--, j--) {
    result[i] = (char) tolower(genomic_diff[j]);
#if 0
    printf("k=%d i=%d result[i]=%c\n",k,i,result[i]);
#endif
    assert(result[i] == 'a' || result[i] == 'c' || result[i] == 'g' || result[i] == 't' || result[i] == 'n');
    }
  for (k = 0, i = queryend, j = alignoffset+queryend-querystart; k < extraright && i < querylength; k++, i++, j++) {
    result[i] = (char) tolower(genomic_diff[j]);
#if 0
    printf("k=%d i=%d result[i]=%c\n",k,i,result[i]);
#endif
    assert(result[i] == 'a' || result[i] == 'c' || result[i] == 'g' || result[i] == 't' || result[i] == 'n');
  }
  debug2(printf("g1: %s\n",result));

  return result;
}


static char *
embellish_genomic_sam (char *genomic_diff, char *query, int querystart, int queryend, int querylength,
		       int genomiclength, int alignoffset, bool cmetp, bool plusp) {
  char *result;
  int i, j, k;

  result = (char *) CALLOC(querylength+1,sizeof(char));
  result[querylength] = '\0';

  strncpy(&(result[querystart]),&(genomic_diff[alignoffset]),queryend-querystart);

  if (cmetp == true) {
    mark_mismatches_met_sam(result,query,querystart,queryend);
  }

  /* Add terminal ends as lower-case */
  for (k = 0, i = querystart-1, j = alignoffset-1; i >= 0 && j >= 0; k++, i--, j--) {
    if (query[i] == genomic_diff[j]) {
      result[i] = genomic_diff[j];
    } else {
      result[i] = (char) tolower(genomic_diff[j]);
    }
#if 0
    printf("k=%d i=%d j=%d result[i]=%c\n",k,i,j,result[i]);
#endif
  }

  for (k = 0, i = queryend, j = alignoffset+queryend-querystart; i < querylength && j < genomiclength; k++, i++, j++) {
    if (query[i] == genomic_diff[j]) {
      result[i] = genomic_diff[j];
    } else {
      result[i] = (char) tolower(genomic_diff[j]);
    }
#if 0
    printf("k=%d i=%d j=%d result[i]=%c\n",k,i,j,result[i]);
#endif
  }

  return result;
}




T
Substring_new (int nmismatches_whole, int ncolordiffs, Chrnum_T chrnum, Genomicpos_T chroffset,
	       Genomicpos_T left, Genomicpos_T genomicstart, Genomicpos_T genomicend,
	       Compress_T query_compress, UINT4 *genome_blocks, UINT4 *snp_blocks,
	       Endtype_T start_endtype, Endtype_T end_endtype,
	       int querystart, int queryend, int querylength,
	       Genomicpos_T alignstart, Genomicpos_T alignend, int genomiclength,
	       int extraleft, int extraright, bool exactp, char *query,
	       bool plusp, bool trim_left_p, bool trim_right_p,
	       int minlength, bool dibasep, bool cmetp) {
  T new = (T) MALLOC(sizeof(*new));
  int aligndiff;

  new->exactp = exactp;
  new->nmismatches_whole = nmismatches_whole;
  new->ncolordiffs = ncolordiffs;
  new->chrnum = chrnum;
  new->chroffset = chroffset;

  new->left_genomicseg = left;
  new->genomicstart = genomicstart;
  new->genomicend = genomicend;

  new->start_endtype = start_endtype;
  new->end_endtype = end_endtype;

  new->querystart_orig = new->querystart = querystart;
  new->queryend_orig = new->queryend = queryend;
  new->querylength = querylength;

  new->alignstart = new->alignstart_trim = alignstart;
  new->alignend = new->alignend_trim = alignend;

  new->extraleft = extraleft;
  new->extraright = extraright;

  new->genomiclength = genomiclength;
  new->plusp = plusp;

  new->trim_left = 0;
  new->trim_right = 0;


  /* Compute coordinates */
  if (plusp == true) {
    new->alignoffset = alignstart - genomicstart;
    aligndiff = /* (alignstart - genomicstart) - querystart = */ new->alignoffset - querystart;
    left += aligndiff;
    new->left = left;

    debug2(printf("\n"));
    debug2(printf("querylength is %d, genomiclength is %d, alignstart is %u, alignend is %u, genomicstart is %u, genomicend is %u, alignoffset is %d\n",
		  querylength,genomiclength,alignstart,alignend,genomicstart,genomicend,new->alignoffset));
    debug2(printf("query:  %s\n",query));
    
  } else {
    new->alignoffset = genomicstart - alignstart;
    aligndiff = (alignend - genomicend) - (querylength - queryend);
    left += aligndiff;
    new->left = left;

    debug2(printf("\n"));
    debug2(printf("querylength is %d, genomiclength is %d, alignstart is %u, alignend is %u, genomicstart is %u, genomicend is %u, alignoffset is %d\n",
		  querylength,genomiclength,alignstart,alignend,genomicstart,genomicend,new->alignoffset));
    debug2(printf("query:  %s\n",query));
  }


  /* Do trimming */
  debug8(printf("trim_left_p %d, trim_right_p %d\n",trim_left_p,trim_right_p));

  if (trim_mismatch_score < 0 && trim_left_p == true) {
    new->trim_left = trim_left_end(query,query_compress,genome_blocks,snp_blocks,left,
				   querystart,queryend,querylength,dibasep,cmetp,plusp);
    new->querystart += new->trim_left;
    if (plusp == true) {
      new->alignstart_trim += new->trim_left;
    } else {
      new->alignstart_trim -= new->trim_left;
    }
  }

  if (trim_mismatch_score < 0 && trim_right_p == true) {
    new->trim_right = trim_right_end(query,query_compress,genome_blocks,snp_blocks,left,
				     querystart,queryend,querylength,dibasep,cmetp,plusp);
    new->queryend -= new->trim_right;
    if (plusp == true) {
      new->alignend_trim -= new->trim_right;
    } else {
      new->alignend_trim += new->trim_right;
    }
  }

  /* Initialize these so Substring_free knows what to do */
  new->genomic_bothdiff = (char *) NULL;
  new->genomic_refdiff = (char *) NULL;


  /* Check for minlength.  Needed to avoid nonsensical terminal alignments */
  if (new->queryend - new->querystart <= minlength) {
    debug8(printf("queryend %d - querystart %d <= minlength %d, so returning NULL\n",
		  new->queryend,new->querystart,minlength));
    Substring_free(&new);

    return (T) NULL;
  }



  /* nmismatches_bothdiff */
  if (new->trim_left == 0 && new->trim_right == 0) {
    new->nmismatches_bothdiff = nmismatches_whole;

  } else if (plusp == true) {
    new->nmismatches_bothdiff = 
      Genome_count_mismatches_substring(&ncolordiffs,query,query_compress,genome_blocks,snp_blocks,left,
					/*pos5*/new->alignstart_trim-left,/*pos3*/new->alignend_trim-left,
					dibasep,cmetp,plusp);

    /* Re-compute nmismatches_whole for terminal alignments */
    if (start_endtype == TERM) {
      new->nmismatches_whole = 
	Genome_count_mismatches_substring(&ncolordiffs,query,query_compress,genome_blocks,snp_blocks,left,
					  /*pos5*/new->alignstart_trim-left,/*pos3*/new->alignend-left,
					  dibasep,cmetp,plusp);

    } else if (end_endtype == TERM) {
      new->nmismatches_whole = 
	Genome_count_mismatches_substring(&ncolordiffs,query,query_compress,genome_blocks,snp_blocks,left,
					  /*pos5*/new->alignstart-left,/*pos3*/new->alignend_trim-left,
					  dibasep,cmetp,plusp);
    }

  } else {
    new->nmismatches_bothdiff = 
      Genome_count_mismatches_substring(&ncolordiffs,query,query_compress,genome_blocks,snp_blocks,left,
					/*pos5*/new->alignend_trim-left,/*pos3*/new->alignstart_trim-left,
					dibasep,cmetp,plusp);

    /* Re-compute nmismatches_whole for terminal alignments */
    if (start_endtype == TERM) {
      new->nmismatches_whole = 
	Genome_count_mismatches_substring(&ncolordiffs,query,query_compress,genome_blocks,snp_blocks,left,
					  /*pos5*/new->alignend-left,/*pos3*/new->alignstart_trim-left,
					  dibasep,cmetp,plusp);
    } else if (end_endtype == TERM) {
      new->nmismatches_whole = 
	Genome_count_mismatches_substring(&ncolordiffs,query,query_compress,genome_blocks,snp_blocks,left,
					  /*pos5*/new->alignend_trim-left,/*pos3*/new->alignstart-left,
					  dibasep,cmetp,plusp);
    }
  }

  /* Counts matches only within trimmed region */
  new->nmatches = (new->queryend - new->querystart) - new->nmismatches_bothdiff;

  return new;
}

double
Substring_compute_mapq (T this, char *query, Compress_T query_compress,
			UINT4 *genome_blocks, UINT4 *snp_blocks, char *quality_string,
			bool dibasep, bool cmetp) {
  int mapq_start, mapq_end;

  /* mapq */
  mapq_start = this->querystart_orig;
  mapq_end = this->queryend_orig;

  if (this->start_endtype == TERM) {
    mapq_start += this->trim_left;
    /* mapq_end -= this->trim_right; */ /* Don't allow trimming found on right */

  } else if (this->end_endtype == TERM) {
    /* mapq_start += this->trim_left; */ /* Don't allow trimming found on left */
    mapq_end -= this->trim_right;

  } else {
    mapq_start += this->trim_left;
    mapq_end -= this->trim_right;
  }


  if (this->exactp == true) {
    /* this->mapq_loglik = MAPQ_loglik_exact(quality_string,0,querylength); */
    this->mapq_loglik = 0.0;
  } else {
    debug2(printf("trim_left %d, trim_right %d, mapq_start = %d, mapq_end = %d\n",
		  this->trim_left,this->trim_right,mapq_start,mapq_end));
    this->mapq_loglik = MAPQ_loglik(query,query_compress,genome_blocks,snp_blocks,
				    this->left,mapq_start,mapq_end,this->querylength,
				    quality_string,dibasep,cmetp,this->plusp);
  }

  return this->mapq_loglik;
}


/* Note: query needed only for dibase */
int
Substring_display_prep (char **deletion, T this, char *query, Compress_T query_compress_fwd, Compress_T query_compress_rev,
			UINT4 *genome_blocks, UINT4 *snp_blocks, Genome_T genome, int deletion_pos, int deletion_length,
			bool dibasep, bool cmetp) {
  char *genomic_diff, gbuffer_alloc[MAX_QUERYLENGTH+1], genomicsegrc_alloc[MAX_QUERYLENGTH+1];
  char *gbuffer, *genomicsegrc;
  int mismatch_offset;
  int ncolordiffs;
  bool allocp;


  mismatch_offset = this->alignoffset - this->querystart_orig;

  /* genomic_bothdiff, genomic_refdiff, and nmismatches_refdiff */
  if (this->exactp == true) {
    this->genomic_bothdiff = (char *) NULL;
    this->genomic_refdiff = (char *) NULL;
    this->nmismatches_refdiff = this->nmismatches_whole;

  } else if (this->plusp == true) {
    if (this->genomiclength < MAX_QUERYLENGTH) {
      gbuffer = gbuffer_alloc;
      allocp = false;
    } else {
      gbuffer = (char *) CALLOC(this->genomiclength+1,sizeof(char));
      allocp = true;
    }

    Genome_fill_buffer_simple(genome,this->left_genomicseg,this->genomiclength,gbuffer);
    if (deletion_pos >= 0) {
      *deletion = (char *) CALLOC(deletion_length+1,sizeof(char));
      strncpy(*deletion,&(gbuffer[deletion_pos]),deletion_length);
    }
    genomic_diff = gbuffer;

    Genome_mark_mismatches(genomic_diff,this->querylength,query_compress_fwd,genome_blocks,snp_blocks,
			   this->left,/*pos5*/this->querystart_orig,/*pos3*/this->queryend_orig,
			   mismatch_offset,dibasep,cmetp,/*plusp*/true);

    this->genomic_bothdiff = embellish_genomic(genomic_diff,query,this->querystart_orig,this->queryend_orig,
					       this->querylength,this->alignoffset,this->extraleft,this->extraright,
					       cmetp,/*plusp*/true);

    if (snp_blocks == NULL) {
      this->genomic_refdiff = this->genomic_bothdiff;
      this->nmismatches_refdiff = this->nmismatches_bothdiff;

    } else {
      this->nmismatches_refdiff =
	Genome_count_mismatches_substring(&ncolordiffs,query,query_compress_fwd,genome_blocks,/*snp_blocks*/NULL,
					  this->left,/*pos5*/this->alignstart_trim - this->left,
					  /*pos3*/this->alignend_trim - this->left,
					  dibasep,cmetp,/*plusp*/true);

      Genome_mark_mismatches(genomic_diff,this->querylength,query_compress_fwd,
			     genome_blocks,/*snp_blocks*/NULL,this->left,
			     /*pos5*/this->querystart_orig,/*pos3*/this->queryend_orig,
			     mismatch_offset,dibasep,cmetp,/*plusp*/true);
      if (output_sam_p == false) {
	this->genomic_refdiff = embellish_genomic(genomic_diff,query,this->querystart_orig,this->queryend_orig,
						  this->querylength,this->alignoffset,this->extraleft,this->extraright,
						  cmetp,/*plusp*/true);
      }
    }

    if (output_sam_p == true) {
      this->genomic_refdiff = embellish_genomic_sam(genomic_diff,query,this->querystart_orig,this->queryend_orig,
						    this->querylength,this->genomiclength,this->alignoffset,
						    cmetp,/*plusp*/true);
    }

    if (allocp == true) {
      FREE(gbuffer);
    }

  } else {
    if (this->genomiclength < MAX_QUERYLENGTH) {
      gbuffer = gbuffer_alloc;
      genomicsegrc = genomicsegrc_alloc;
      allocp = false;
    } else {
      gbuffer = (char *) CALLOC(this->genomiclength+1,sizeof(char));
      genomicsegrc = (char *) CALLOC(this->genomiclength+1,sizeof(char));
      allocp = true;
    }

    Genome_fill_buffer_simple(genome,this->left_genomicseg,this->genomiclength,gbuffer);
    genomic_diff = make_complement_buffered(genomicsegrc,gbuffer,this->genomiclength);
    if (deletion_pos >= 0) {
      *deletion = (char *) CALLOC(deletion_length+1,sizeof(char));
      strncpy(*deletion,&(genomicsegrc[deletion_pos]),deletion_length);
    }

    Genome_mark_mismatches(genomic_diff,this->querylength,query_compress_rev,genome_blocks,snp_blocks,
			   this->left,/*pos5*/this->querylength - this->queryend_orig,
			   /*pos3*/this->querylength - this->querystart_orig,
			   mismatch_offset,dibasep,cmetp,/*plusp*/false);

    this->genomic_bothdiff = embellish_genomic(genomic_diff,query,this->querystart_orig,this->queryend_orig,
					       this->querylength,this->alignoffset,this->extraleft,this->extraright,
					       cmetp,/*plusp*/false);

    if (snp_blocks == NULL) {
      this->genomic_refdiff = this->genomic_bothdiff;
      this->nmismatches_refdiff = this->nmismatches_bothdiff;

    } else {
      this->nmismatches_refdiff = 
	Genome_count_mismatches_substring(&ncolordiffs,query,query_compress_rev,genome_blocks,/*snp_blocks*/NULL,
					  this->left,/*pos5*/this->alignend_trim - this->left,
					  /*pos3*/this->alignstart_trim - this->left,dibasep,cmetp,/*plusp*/false);

      Genome_mark_mismatches(genomic_diff,this->querylength,query_compress_rev,
			     genome_blocks,/*snp_blocks*/NULL,this->left,
			     /*pos5*/this->querylength - this->queryend_orig,
			     /*pos3*/this->querylength - this->querystart_orig,
			     mismatch_offset,dibasep,cmetp,/*plusp*/false);

      if (output_sam_p == false) {
	this->genomic_refdiff = embellish_genomic(genomic_diff,query,this->querystart_orig,this->queryend_orig,
						  this->querylength,this->alignoffset,this->extraleft,this->extraright,
						  cmetp,/*plusp*/false);
      }
    }

    if (output_sam_p == true) {
      this->genomic_refdiff = embellish_genomic_sam(genomic_diff,query,this->querystart_orig,this->queryend_orig,
						    this->querylength,this->genomiclength,this->alignoffset,
						    cmetp,/*plusp*/false);
    }

    if (allocp == true) {
      FREE(genomicsegrc);
      FREE(gbuffer);
    }
  }

  return this->nmismatches_refdiff;
}


int
Substring_splicesites_i (T this) {
  return this->splicesites_i;
}

int
Substring_splicesites_i_A (T this) {
  return this->splicesites_i;
}

int
Substring_splicesites_i_D (T this) {
  return this->splicesites_i_2;
}

bool
Substring_plusp (T this) {
  return this->plusp;
}

char *
Substring_genomic_bothdiff (T this) {
  return this->genomic_bothdiff;
}

char *
Substring_genomic_refdiff (T this) {
  return this->genomic_refdiff;
}


int
Substring_nmismatches_whole (T this) {
  return this->nmismatches_whole;
}


int
Substring_nmismatches_bothdiff (T this) {
  return this->nmismatches_bothdiff;
}


int
Substring_nmismatches_refdiff (T this) {
  return this->nmismatches_refdiff;
}

int
Substring_nmatches (T this) {
  return this->nmatches;
}

int
Substring_ncolordiffs (T this) {
  return this->ncolordiffs;
}

Endtype_T
Substring_start_endtype (T this) {
  return this->start_endtype;
}

Endtype_T
Substring_end_endtype (T this) {
  return this->end_endtype;
}

void
Substring_set_start_endtype (T this, Endtype_T start_endtype) {
  this->start_endtype = start_endtype;
  return;
}

void
Substring_set_end_endtype (T this, Endtype_T end_endtype) {
  this->end_endtype = end_endtype;
  return;
}

double
Substring_mapq_loglik (T this) {
  return this->mapq_loglik;
}

int
Substring_trim_left (T this) {
  return this->trim_left;
}

int
Substring_trim_right (T this) {
  return this->trim_right;
}


int
Substring_querystart (T this) {
  return this->querystart;
}

int
Substring_querystart_orig (T this) {
  return this->querystart_orig;
}

int
Substring_queryend (T this) {
  return this->queryend;
}

int
Substring_querylength (T this) {
  return this->querylength;
}

int
Substring_match_length (T this) {
  if (this == NULL) {
    return 0;
  } else {
    return this->queryend - this->querystart;
  }
}

/* Before trimming */
int
Substring_match_length_orig (T this) {
  if (this == NULL) {
    return 0;
  } else {
    return this->queryend_orig - this->querystart_orig;
  }
}

/* Used only by Goby */
Genomicpos_T
Substring_genomic_alignment_length (T this) {
  if (this == NULL) {
    return 0U;
  } else if (this->alignend_trim > this->alignstart_trim) {
    return this->alignend_trim - this->alignstart_trim;
  } else {
    return this->alignstart_trim - this->alignend_trim;
  }
}


Chrnum_T
Substring_chrnum (T this) {
  return this->chrnum;
}

Genomicpos_T
Substring_chroffset (T this) {
  return this->chroffset;
}

Genomicpos_T
Substring_alignstart_trim (T this) {
  return this->alignstart_trim;
}

Genomicpos_T
Substring_alignend_trim (T this) {
  return this->alignend_trim;
}


Genomicpos_T
Substring_genomicstart (T this) {
  return this->genomicstart;
}

Genomicpos_T
Substring_genomicend (T this) {
  return this->genomicend;
}

Genomicpos_T
Substring_genomiclength (T this) {
  return this->genomiclength;
}

double
Substring_chimera_prob (T this) {
  return this->chimera_prob;
}

int
Substring_chimera_pos (T this) {
  return this->chimera_pos;
}

/* For shortexon */
int
Substring_chimera_pos_A (T this) {
  return this->chimera_pos;
}

/* For shortexon */
int
Substring_chimera_pos_D (T this) {
  return this->chimera_pos_2;
}

bool
Substring_chimera_knownp (T this) {
  return this->chimera_knownp;
}

bool
Substring_chimera_sensep (T this) {
  return this->chimera_sensep;
}



T
Substring_copy (T old) {
  T new;

  if (old == NULL) {
    return NULL;
  } else {
    new = (T) MALLOC(sizeof(*new));

    new->exactp = old->exactp;
    new->nmismatches_whole = old->nmismatches_whole;
    new->nmismatches_bothdiff = old->nmismatches_bothdiff;
    new->nmismatches_refdiff = old->nmismatches_refdiff;
    new->nmatches = old->nmatches;

    new->ncolordiffs = old->ncolordiffs;
    new->trim_left = old->trim_left;
    new->trim_right = old->trim_right;
    new->chrnum = old->chrnum;
    new->chroffset = old->chroffset;

    new->left_genomicseg = old->left_genomicseg;
    new->left = old->left;
    new->genomicstart = old->genomicstart;
    new->genomicend = old->genomicend;

    new->start_endtype = old->start_endtype;
    new->end_endtype = old->end_endtype;

    new->querystart_orig = old->querystart_orig;
    new->queryend_orig = old->queryend_orig;
    new->querystart = old->querystart;
    new->queryend = old->queryend;
    new->querylength = old->querylength;

    new->alignstart = old->alignstart;
    new->alignend = old->alignend;

    new->alignstart_trim = old->alignstart_trim;
    new->alignend_trim = old->alignend_trim;

    new->alignoffset = old->alignoffset;
    new->extraleft = old->extraleft;
    new->extraright = old->extraright;

    new->genomiclength = old->genomiclength;
    new->plusp = old->plusp;

    if (old->genomic_bothdiff == NULL) {
      new->genomic_bothdiff = (char *) NULL;
      new->genomic_refdiff = (char *) NULL;
    } else if (old->genomic_refdiff == old->genomic_bothdiff) {
      new->genomic_bothdiff = (char *) CALLOC(strlen(old->genomic_bothdiff)+1,sizeof(char));
      strcpy(new->genomic_bothdiff,old->genomic_bothdiff);
      new->genomic_refdiff = new->genomic_bothdiff;
    } else {
      new->genomic_bothdiff = (char *) CALLOC(strlen(old->genomic_bothdiff)+1,sizeof(char));
      strcpy(new->genomic_bothdiff,old->genomic_bothdiff);
      new->genomic_refdiff = (char *) CALLOC(strlen(old->genomic_refdiff)+1,sizeof(char));
      strcpy(new->genomic_refdiff,old->genomic_refdiff);
    }

    new->mapq_loglik = old->mapq_loglik;

    new->chimera_sensep = old->chimera_sensep;

    new->splicesites_i = old->splicesites_i;
    new->chimera_knownp = old->chimera_knownp;
    new->chimera_modelpos = old->chimera_modelpos;
    new->chimera_pos = old->chimera_pos;
    new->chimera_prob = old->chimera_prob;

    new->splicesites_i_2 = old->splicesites_i_2;
    new->chimera_knownp_2 = old->chimera_knownp_2;
    new->chimera_modelpos_2 = old->chimera_modelpos_2;
    new->chimera_pos_2 = old->chimera_pos_2;
    new->chimera_prob_2 = old->chimera_prob_2;

    return new;
  }
}



T
Substring_new_donor (int splicesites_i, int splicesites_offset, int donor_pos, int donor_nmismatches, int donor_ncolordiffs,
		     double donor_prob, Genomicpos_T left, Compress_T query_compress,
		     UINT4 *genome_blocks, UINT4 *snp_blocks, int querylength, bool plusp, bool sensep,
		     char *query, Chrnum_T chrnum, Genomicpos_T chroffset, bool dibasep, bool cmetp) {
  T new;
  int querystart, queryend, extraleft, extraright;
  Genomicpos_T genomicstart, genomicend, alignstart, alignend;
  Endtype_T start_endtype, end_endtype;
  bool trim_left_p, trim_right_p;

  if (dibasep) {
  } else {
    donor_ncolordiffs = 0;
  }

  if (plusp == true) {
    genomicstart = left;
    genomicend = left + querylength;
    if (sensep == true) {
      start_endtype = END;
      end_endtype = DON;

      querystart = 0;
      queryend = donor_pos;
      extraleft = 0;
      extraright = 2;
      alignstart = genomicstart;
      alignend = genomicstart + donor_pos;
      trim_left_p = true;	/* querystart == 0 */
      trim_right_p = false;

    } else {
      extraleft = 2;
      extraright = 0;
      start_endtype = DON;
      end_endtype = END;

      querystart = donor_pos;
      queryend = querylength;
      alignstart = genomicstart + donor_pos;
      alignend = genomicend;
      trim_left_p = false;
      trim_right_p = true;	/* queryend == querylength */
    }

  } else {
    genomicstart = left + querylength;
    genomicend = left;
    if (sensep == true) {
      start_endtype = END;
      end_endtype = DON;

      querystart = 0;
      queryend = querylength - donor_pos;
      extraleft = 0;
      extraright = 2;
      alignstart = genomicstart;
      alignend = genomicstart - (querylength - donor_pos);
      trim_left_p = true;	/* querystart == 0 */
      trim_right_p = false;

    } else {
      extraleft = 2;
      extraright = 0;
      start_endtype = DON;
      end_endtype = END;

      querystart = querylength - donor_pos;
      queryend = querylength;
      alignstart = genomicstart - (querylength - donor_pos);
      alignend = genomicend;
      trim_left_p = false;
      trim_right_p = true;	/* queryend == querylength */
    }
  }

  if ((new = Substring_new(donor_nmismatches,donor_ncolordiffs,chrnum,chroffset,
			   left,genomicstart,genomicend,query_compress,genome_blocks,snp_blocks,
			   start_endtype,end_endtype,querystart,queryend,querylength,
			   alignstart,alignend,/*genomiclength*/querylength,
			   extraleft,extraright,/*exactp*/false,query,plusp,
			   trim_left_p,trim_right_p,/*minlength*/0,dibasep,cmetp)) == NULL) {
    return (T) NULL;
  }

  debug2(printf("Making new donor with splicesites_i %d and left %u, plusp %d\n",splicesites_i,left,plusp));
  new->chimera_modelpos = left + donor_pos;
  new->chimera_sensep = sensep;
  if (splicesites_i >= 0) {
    new->chimera_knownp = true;
    new->splicesites_i = splicesites_offset + splicesites_i;
  } else {
    new->chimera_knownp = false;
    new->splicesites_i = -1;
  }

  if (plusp == true) {
    new->chimera_pos = donor_pos;
  } else {
    new->chimera_pos = querylength - donor_pos;
  }
  new->chimera_prob = donor_prob;

  return new;
}


T
Substring_new_acceptor (int splicesites_i, int splicesites_offset, int acceptor_pos, int acceptor_nmismatches, int acceptor_ncolordiffs,
			double acceptor_prob, Genomicpos_T left, Compress_T query_compress,
			UINT4 *genome_blocks, UINT4 *snp_blocks, int querylength, bool plusp, bool sensep,
			char *query, Chrnum_T chrnum, Genomicpos_T chroffset, bool dibasep, bool cmetp) {
  T new;
  int querystart, queryend, extraleft, extraright;
  Genomicpos_T genomicstart, genomicend, alignstart, alignend;
  Endtype_T start_endtype, end_endtype;
  bool trim_left_p, trim_right_p;

  if (dibasep) {
  } else {
    acceptor_ncolordiffs = 0;
  }

  if (plusp == true) {
    genomicstart = left;
    genomicend = left + querylength;
    if (sensep == true) {
      extraleft = 2;
      extraright = 0;
      start_endtype = ACC;
      end_endtype = END;

      querystart = acceptor_pos;
      queryend = querylength;
      alignstart = genomicstart + acceptor_pos;
      alignend = genomicend;
      trim_left_p = false;
      trim_right_p = true;	/* queryend == querylength */

    } else {
      extraleft = 0;
      extraright = 2;
      start_endtype = END;
      end_endtype = ACC;

      querystart = 0;
      queryend = acceptor_pos;
      alignstart = genomicstart;
      alignend = genomicstart + acceptor_pos;
      trim_left_p = true;	/* querystart == 0 */
      trim_right_p = false;
    }

  } else {
    genomicstart = left + querylength;
    genomicend = left;
    if (sensep == true) {
      extraleft = 2;
      extraright = 0;
      start_endtype = ACC;
      end_endtype = END;

      querystart = querylength - acceptor_pos;
      queryend = querylength;
      alignstart = genomicstart - (querylength - acceptor_pos);
      alignend = genomicend;
      trim_left_p = false;
      trim_right_p = true;	/* queryend == querylength */

    } else {
      extraleft = 0;
      extraright = 2;
      start_endtype = END;
      end_endtype = ACC;

      querystart = 0;
      queryend = querylength - acceptor_pos;
      alignstart = genomicstart;
      alignend = genomicstart - (querylength - acceptor_pos);
      trim_left_p = true;	/* querystart == 0 */
      trim_right_p = false;
    }
  }

  if ((new = Substring_new(acceptor_nmismatches,acceptor_ncolordiffs,chrnum,chroffset,
			   left,genomicstart,genomicend,query_compress,genome_blocks,snp_blocks,
			   start_endtype,end_endtype,querystart,queryend,querylength,
			   alignstart,alignend,/*genomiclength*/querylength,
			   extraleft,extraright,/*exactp*/false,query,plusp,
			   trim_left_p,trim_right_p,/*minlength*/0,dibasep,cmetp)) == NULL) {
    return (T) NULL;
  }

  debug2(printf("Making new acceptor with splicesites_i %d and left %u, plusp %d\n",splicesites_i,left,plusp));
  new->chimera_modelpos = left + acceptor_pos;
  new->chimera_sensep = sensep;
  if (splicesites_i >= 0) {
    new->chimera_knownp = true;
    new->splicesites_i = splicesites_i + splicesites_offset;
  } else {
    new->chimera_knownp = false;
    new->splicesites_i = -1;
  }

  if (plusp == true) {
    new->chimera_pos = acceptor_pos;
  } else {
    new->chimera_pos = querylength - acceptor_pos;
  }
  new->chimera_prob = acceptor_prob;

  return new;
}



T
Substring_new_shortexon (int acceptor_splicesites_i, int donor_splicesites_i, int splicesites_offset,
			 int acceptor_pos, int donor_pos, int nmismatches, int ncolordiffs,
			 double acceptor_prob, double donor_prob, Genomicpos_T left,
			 Compress_T query_compress, UINT4 *genome_blocks, UINT4 *snp_blocks,
			 int querylength, bool plusp, bool sensep,
			 bool acceptor_ambp, bool donor_ambp, char *query,
			 Chrnum_T chrnum, Genomicpos_T chroffset,
			 bool dibasep, bool cmetp) {
  T new;
  int querystart, queryend;
  Genomicpos_T genomicstart, genomicend, alignstart, alignend;
  Endtype_T start_endtype, end_endtype;

  if (dibasep) {
  } else {
    ncolordiffs = 0;
  }

  if (plusp == true) {
    genomicstart = left;
    genomicend = left + querylength;
    if (sensep == true) {
      start_endtype = (acceptor_ambp == true) ? AMB_ACC : ACC;
      end_endtype = (donor_ambp == true) ? AMB_DON : DON;
      querystart = acceptor_pos;
      queryend = donor_pos;
      alignstart = genomicstart + acceptor_pos;
      alignend = genomicstart + donor_pos;
    } else {
      start_endtype = (donor_ambp == true) ? AMB_DON : DON;
      end_endtype = (acceptor_ambp == true) ? AMB_ACC : ACC;
      querystart = donor_pos;
      queryend = acceptor_pos;
      alignstart = genomicstart + donor_pos;
      alignend = genomicstart + acceptor_pos;
    }

  } else {
    genomicstart = left + querylength;
    genomicend = left;
    if (sensep == true) {
      start_endtype = (acceptor_ambp == true) ? AMB_ACC : ACC;
      end_endtype = (donor_ambp == true) ? AMB_DON : DON;
      querystart = querylength - acceptor_pos;
      queryend = querylength - donor_pos;
      alignstart = genomicstart - (querylength - acceptor_pos);
      alignend = genomicstart - (querylength - donor_pos);
    } else {
      start_endtype = (donor_ambp == true) ? AMB_DON : DON;
      end_endtype = (acceptor_ambp == true) ? AMB_ACC : ACC;
      querystart = querylength - donor_pos;
      queryend = querylength - acceptor_pos;
      alignstart = genomicstart - (querylength - donor_pos);
      alignend = genomicstart - (querylength - acceptor_pos);
    }
  }

  if ((new = Substring_new(nmismatches,ncolordiffs,chrnum,chroffset,
			   left,genomicstart,genomicend,query_compress,genome_blocks,snp_blocks,
			   start_endtype,end_endtype,querystart,queryend,querylength,
			   alignstart,alignend,/*genomiclength*/querylength,
			   /*extraleft*/2,/*extraright*/2,/*exactp*/false,query,plusp,
			   /*trim_left_p*/false,/*trim_right_p*/false,
			   /*minlength*/0,dibasep,cmetp)) == NULL) {
    return (T) NULL;
  }

  debug2(printf("Making new middle with left %u, plusp %d\n",left,plusp));
  new->chimera_modelpos = left + acceptor_pos;
  new->chimera_modelpos_2 = left + donor_pos;
  new->chimera_sensep = sensep;

  if (acceptor_splicesites_i >= 0) {
    new->chimera_knownp = true;
    new->splicesites_i = acceptor_splicesites_i + splicesites_offset;
  } else {
    new->chimera_knownp = false;
    new->splicesites_i = -1;
  }

  if (donor_splicesites_i >= 0) {
    new->chimera_knownp_2 = true;
    new->splicesites_i_2 = donor_splicesites_i + splicesites_offset;
  } else {
    new->chimera_knownp_2 = false;
    new->splicesites_i_2 = -1;
  }

  if (plusp == true) {
    new->chimera_pos = acceptor_pos;
    new->chimera_pos_2 = donor_pos;
  } else {
    new->chimera_pos = querylength - acceptor_pos;
    new->chimera_pos_2 = querylength - donor_pos;
  }

  new->chimera_prob = acceptor_prob;
  new->chimera_prob_2 = donor_prob;

  return new;
}



void
Substring_assign_donor_prob (T donor, UINT4 *genome_blocks) {

  if (donor == NULL) {
    return;
  }

  if (donor->chimera_knownp == false) {
    /* Prob already assigned */
  } else if (donor->plusp == donor->chimera_sensep) {
    donor->chimera_prob = Maxent_hr_donor_prob(donor->chimera_modelpos,genome_blocks);
  } else {
    donor->chimera_prob = Maxent_hr_antidonor_prob(donor->chimera_modelpos,genome_blocks);
  }

  return;
}

void
Substring_assign_acceptor_prob (T acceptor, UINT4 *genome_blocks) {

  if (acceptor == NULL) {
    return;
  }

  if (acceptor->chimera_knownp == false) {
    /* Prob already assigned */
  } else if (acceptor->plusp == acceptor->chimera_sensep) {
    acceptor->chimera_prob = Maxent_hr_acceptor_prob(acceptor->chimera_modelpos,genome_blocks);
  } else {
    acceptor->chimera_prob = Maxent_hr_antiacceptor_prob(acceptor->chimera_modelpos,genome_blocks);
  }

  return;
}


void
Substring_assign_shortexon_prob (T shortexon, UINT4 *genome_blocks) {

  if (shortexon->chimera_knownp == false) {
    /* Prob1 already assigned */
  } else if (shortexon->plusp == shortexon->chimera_sensep) {
    shortexon->chimera_prob = Maxent_hr_acceptor_prob(shortexon->chimera_modelpos,genome_blocks);
  } else {
    shortexon->chimera_prob = Maxent_hr_antiacceptor_prob(shortexon->chimera_modelpos,genome_blocks);
  }


  if (shortexon->chimera_knownp_2 == false) {
    /* Prob2 already assigned */
  } else if (shortexon->plusp == shortexon->chimera_sensep) {
    shortexon->chimera_prob_2 = Maxent_hr_donor_prob(shortexon->chimera_modelpos_2,genome_blocks);
  } else {
    shortexon->chimera_prob_2 = Maxent_hr_antidonor_prob(shortexon->chimera_modelpos_2,genome_blocks);
  }

  return;
}



static int
ascending_pos_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->chimera_pos < y->chimera_pos) {
    return -1;
  } else if (x->chimera_pos > y->chimera_pos) {
    return +1;
  } else if (x->genomicstart < y->genomicstart) {
    return -1;
  } else if (x->genomicstart > y->genomicstart) {
    return +1;
  } else if (x->chimera_knownp == true && y->chimera_knownp == false) {
    return -1;
  } else if (y->chimera_knownp == true && x->chimera_knownp == false) {
    return +1;
  } else {
    return 0;
  }
}

static int
descending_pos_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->chimera_pos < y->chimera_pos) {
    return -1;
  } else if (x->chimera_pos > y->chimera_pos) {
    return +1;
  } else if (x->genomicstart > y->genomicstart) {
    return -1;
  } else if (x->genomicstart < y->genomicstart) {
    return +1;
  } else if (x->chimera_knownp == true && y->chimera_knownp == false) {
    return -1;
  } else if (y->chimera_knownp == true && x->chimera_knownp == false) {
    return +1;
  } else {
    return 0;
  }
}

List_T
Substring_sort_chimera_halves (List_T hitlist, bool ascendingp) {
  List_T sorted = NULL, p;
  T x, *hits;
  int n, i, j;
  bool *eliminate;

  n = List_length(hitlist);
  debug(printf("Checking %d spliceends for duplicates...",n));
  if (n == 0) {
    debug(printf("\n"));
    return NULL;
  }

  hits = (T *) CALLOC(n,sizeof(T));
  for (p = hitlist, i = 0; p != NULL; p = p->rest) {
    hits[i++] = (T) p->first;
  }
  List_free(&hitlist);

  if (ascendingp == true) {
    qsort(hits,n,sizeof(T),ascending_pos_cmp);
  } else {
    qsort(hits,n,sizeof(T),descending_pos_cmp);
  }

  /* Check for duplicates */
  eliminate = (bool *) CALLOC(n,sizeof(bool));
  for (i = 0; i < n; i++) {
    x = hits[i];
    j = i+1;
    while (j < n && hits[j]->chimera_pos == x->chimera_pos && hits[j]->genomicstart == x->genomicstart) {
      eliminate[j] = true;
      j++;
    }
  }

  debug(j = 0);
  for (i = n-1; i >= 0; i--) {
    x = hits[i];
    if (eliminate[i] == false) {
      sorted = List_push(sorted,x);
    } else {
      Substring_free(&x);
      debug(j++);
    }
  }
  debug(printf("%d eliminated\n",j));

  FREE(hits);
  FREE(eliminate);

  return sorted;
}


static void
print_snp_labels (FILE *fp, T this, Shortread_T queryseq) {
  int *snps, nsnps, querypos, i;
  char *label, *seq1, *seq2;
  bool allocp, printp = false;
  Interval_T interval;
  Genomicpos_T position;

  if (this->plusp == true) {
    snps = IIT_get_with_divno(&nsnps,snps_iit,
			      snps_divint_crosstable[this->chrnum],
			      this->alignstart - this->chroffset + 1,this->alignend - this->chroffset,
			      /*sortp*/false);
  } else {
    snps = IIT_get_with_divno(&nsnps,snps_iit,
			      snps_divint_crosstable[this->chrnum],
			      this->alignend - this->chroffset + 1,this->alignstart - this->chroffset,
			      /*sortp*/false);
  }

  fprintf(fp,",snps:");

  seq1 = Shortread_fullpointer_uc(queryseq);
  if (this->genomic_bothdiff == NULL) {
    seq2 = Shortread_fullpointer(queryseq);
  } else {
    seq2 = this->genomic_bothdiff;
  }

  if (this->plusp) {

#if 0
    for (i = 0; i < nsnps; i++) {
      interval = IIT_interval(snps_iit,snps[i]);
      position = Interval_low(interval);
      querypos = position - (this->genomicstart-this->chroffset) - 1;
      printf("%d ",querypos);
    }
    printf("\n");
#endif

    for (i = 0; i < nsnps; i++) {
      label = IIT_label(snps_iit,snps[i],&allocp);
    
      interval = IIT_interval(snps_iit,snps[i]);
      position = Interval_low(interval);
      querypos = position - (this->genomicstart-this->chroffset) - 1;
      if (querypos < 0 || querypos >= this->genomiclength) {
	abort();
      }

#if 0
      alleles = IIT_typestring(snps_iit,Interval_type(interval));
      if (c == alleles[0] || c == alleles[1]) ;
#endif

      if (1 || seq1[querypos] != toupper(seq2[querypos])) {
	if (printp) {
	  fprintf(fp,"|");
	}
	fprintf(fp,"%d@",querypos+1);
	fprintf(fp,"%s",label);
	printp = true;
      }

      if (allocp) FREE(label);
    }

  } else {

    for (i = nsnps-1; i >= 0; i--) {
      label = IIT_label(snps_iit,snps[i],&allocp);
    
      interval = IIT_interval(snps_iit,snps[i]);
      position = Interval_low(interval);
      querypos = (this->genomicstart-this->chroffset) - position;
      if (querypos < 0 || querypos >= this->genomiclength) {
	abort();
      }

#if 0
      /* printf("\n%d%c\n",querypos,c); */
      alleles = IIT_typestring(snps_iit,Interval_type(interval));
      if (c == alleles[0] || c == alleles[1]) ;
#endif

      if (1 || seq1[querypos] != toupper(seq2[querypos])) {
	if (printp) {
	  fprintf(fp,"|");
	}
	fprintf(fp,"%d@",querypos+1);
	fprintf(fp,"%s",label);
	printp = true;
      }

      if (allocp) FREE(label);
    }
  }

  FREE(snps);

  return;
}


static void
print_splicesite_labels (FILE *fp, T this, int typeint, int chimera_pos, char *tag) {
  Genomicpos_T splicesitepos;
  int *splicesites, nsplicesites, i;
  char *label;
  bool allocp;

  if (this->plusp == true) {
    splicesitepos = this->genomicstart - this->chroffset + chimera_pos;
  } else {
    splicesitepos = this->genomicstart - this->chroffset - chimera_pos;
  }

  if (this->chimera_knownp == true) {
    /* Note: IIT_get_typed_signed_with_divno does not work here */
    splicesites = IIT_get_exact_multiple_with_divno(&nsplicesites,splicesites_iit,
						    splicesites_divint_crosstable[this->chrnum],
						    splicesitepos,splicesitepos+1U,typeint);
    if (nsplicesites == 0) {
      fprintf(stderr,"Supposed to have a splicesite at chrnum %d, %u..%u, type %d\n",
	      this->chrnum,splicesitepos,splicesitepos+1U,typeint);
    }
    assert(nsplicesites > 0);

    fprintf(fp,",%s:",tag);
    label = IIT_label(splicesites_iit,splicesites[0],&allocp);
    fprintf(fp,"%s",label);
    if (allocp) FREE(label);

    for (i = 1; i < nsplicesites; i++) {
      label = IIT_label(splicesites_iit,splicesites[i],&allocp);
      fprintf(fp,"|%s",label);
      if (allocp) FREE(label);
    }
    FREE(splicesites);
  }

  return;
}


static void
print_forward (FILE *fp, char *string, int n) {
  int i;
  
  for (i = 0; i < n; i++) {
    fprintf(fp,"%c",string[i]);
  }
  return;
}




static void
print_lc (FILE *fp, char *string, int n) {
  int i;
  
  for (i = 0; i < n; i++) {
    fprintf(fp,"%c",(char) tolower(string[i]));
  }
  return;
}



static void
print_revcomp (FILE *fp, char *nt, int len) {
  int i;

  for (i = len-1; i >= 0; --i) {
    fprintf(fp,"%c",complCode[(int) nt[i]]);
  }
  return;
}

static void
print_revcomp_lc (FILE *fp, char *nt, int len) {
  int i;

  for (i = len-1; i >= 0; --i) {
    fprintf(fp,"%c",(char) tolower(complCode[(int) nt[i]]));
  }
  return;
}



static void
print_genomic (FILE *fp, T substring, char *deletion, int deletionlength, bool invertp,
	       Shortread_T queryseq) {
  int i;

  if (invertp == false) {
    if (substring->genomic_bothdiff == NULL) {
      /* Exact match */
      Shortread_print_oneline_uc(fp,queryseq);

    } else if (show_refdiff_p == true) {
      print_forward(fp,substring->genomic_refdiff,substring->queryend);
      if (deletion != NULL) {
	print_lc(fp,deletion,deletionlength);
      }
      print_forward(fp,&(substring->genomic_refdiff[substring->queryend]),substring->querylength - substring->queryend);
    } else {
      print_forward(fp,substring->genomic_bothdiff,substring->queryend);
      if (deletion != NULL) {
	print_lc(fp,deletion,deletionlength);
      }
      print_forward(fp,&(substring->genomic_bothdiff[substring->queryend]),substring->querylength - substring->queryend);
    }

    for (i = 0; i < Shortread_choplength(queryseq); i++) {
      fprintf(fp,"*");
    }
    fprintf(fp,"\t");
    fprintf(fp,"%d..%d",1 + substring->querystart,substring->queryend);

  } else {
    if (substring->genomic_bothdiff == NULL) {
      /* Exact match */
      Shortread_print_oneline_revcomp_uc(fp,queryseq);

    } else if (show_refdiff_p == true) {
      print_revcomp(fp,&(substring->genomic_refdiff[substring->querystart]),substring->querylength - substring->querystart);
      if (deletion != NULL) {
	print_revcomp_lc(fp,deletion,deletionlength);
      }
      print_revcomp(fp,substring->genomic_refdiff,substring->querystart);

    } else {
      print_revcomp(fp,&(substring->genomic_bothdiff[substring->querystart]),substring->querylength - substring->querystart);
      if (deletion != NULL) {
	print_revcomp_lc(fp,deletion,deletionlength);
      }
      print_revcomp(fp,substring->genomic_bothdiff,substring->querystart);
    }
    for (i = 0; i < Shortread_choplength(queryseq); i++) {
      fprintf(fp,"*");
    }
    fprintf(fp,"\t");
    fprintf(fp,"%d..%d",1 + substring->querylength - substring->queryend,
	   substring->querylength - substring->querystart);
  }
  return;
}


static void
print_coordinates (FILE *fp, T substring, char *chr, bool invertp) {

  if (substring->plusp == true) {
    if (invertp == false) {
      fprintf(fp,"+%s:%u..%u",chr,substring->alignstart_trim - substring->chroffset + 1U,
	     substring->alignend_trim - substring->chroffset);
    } else {
      fprintf(fp,"-%s:%u..%u",chr,substring->alignend_trim - substring->chroffset,
	     substring->alignstart_trim - substring->chroffset + 1U);
    }
  } else {
    if (invertp == false) {
      fprintf(fp,"-%s:%u..%u",chr,substring->alignstart_trim - substring->chroffset,
	     substring->alignend_trim - substring->chroffset + 1U);
    } else {
      fprintf(fp,"+%s:%u..%u",chr,substring->alignend_trim - substring->chroffset + 1U,
	     substring->alignstart_trim - substring->chroffset);
    }
  }

  return;
}


void
Substring_print_single (FILE *fp, T substring, Shortread_T queryseq,
			char *chr, int querylength, bool invertp) {

  print_genomic(fp,substring,/*deletion*/(char *) NULL,/*deletionlength*/0,invertp,queryseq);
  fprintf(fp,"\t");
  print_coordinates(fp,substring,chr,invertp);

  fprintf(fp,"\t");
  if (invertp == false) {
    switch (substring->start_endtype) {
    case END: fprintf(fp,"start:%d",substring->trim_left); break;
    case TERM: fprintf(fp,"term:%d",substring->trim_left); break;
    default: fprintf(stderr,"start_endtype is %d\n",substring->start_endtype); abort(); break;
    }
  } else {
    switch (substring->end_endtype) {
    case END: fprintf(fp,"start:%d",substring->trim_right); break;
    case TERM: fprintf(fp,"term:%d",substring->trim_right); break;
    default: fprintf(stderr,"end_endtype is %d\n",substring->end_endtype); abort(); break;
    }
  }

  fprintf(fp,"..");

  if (invertp == false) {
    switch (substring->end_endtype) {
    case END: fprintf(fp,"end:%d",substring->trim_right); break;
    case TERM: fprintf(fp,"term:%d",substring->trim_right); break;
    default: fprintf(stderr,"end_endtype is %d\n",substring->end_endtype); abort(); break;
    }
  } else {
    switch (substring->start_endtype) {
    case END: fprintf(fp,"end:%d",substring->trim_left); break;
    case TERM: fprintf(fp,"term:%d",substring->trim_left); break;
    default: fprintf(stderr,"start_endtype is %d\n",substring->start_endtype); abort(); break;
    }
  }

  fprintf(fp,",sub:%d",substring->nmismatches_bothdiff);
  if (print_ncolordiffs_p) {
    fprintf(fp,"+%d=%d",substring->ncolordiffs,substring->nmismatches_bothdiff+substring->ncolordiffs);
  } else if (print_nsnpdiffs_p) {
    fprintf(fp,"+%d=%d",substring->nmismatches_refdiff - substring->nmismatches_bothdiff,substring->nmismatches_refdiff);
    if (print_snplabels_p && substring->nmismatches_refdiff > substring->nmismatches_bothdiff) {
      print_snp_labels(fp,substring,queryseq);
    }
  }

  return;
}


void
Substring_print_insertion_1 (FILE *fp, T substring1, T substring2, int nindels, 
			     Shortread_T queryseq, char *chr, bool invertp) {
  T substring;

  if (invertp == false) {
    substring = substring1;
    print_genomic(fp,substring1,/*deletion*/NULL,/*deletionlength*/0,/*invertp*/false,
		  queryseq);

  } else {
    substring = substring2;
    print_genomic(fp,substring2,/*deletion*/NULL,/*deletionlength*/0,/*invertp*/true,
		  queryseq);
  }

  fprintf(fp,"\t");

  print_coordinates(fp,substring,chr,invertp);


  fprintf(fp,"\t");
  if (invertp == false) {
    fprintf(fp,"start:%d..ins:%d,sub:%d",substring->trim_left,nindels,substring->nmismatches_bothdiff);
  } else {
    fprintf(fp,"start:%d..ins:%d,sub:%d",substring->trim_right,nindels,substring->nmismatches_bothdiff);
  }
  if (print_ncolordiffs_p) {
    fprintf(fp,"+%d=%d",substring->ncolordiffs,substring->nmismatches_bothdiff+substring->ncolordiffs);
  } else if (print_nsnpdiffs_p) {
    fprintf(fp,"+%d=%d",substring->nmismatches_refdiff - substring->nmismatches_bothdiff,substring->nmismatches_refdiff);
    if (print_snplabels_p && substring->nmismatches_refdiff > substring->nmismatches_bothdiff) {
      print_snp_labels(fp,substring,queryseq);
    }
  }


  return;
}

void
Substring_print_insertion_2 (FILE *fp, T substring1, T substring2, int nindels,
			     Shortread_T queryseq, char *chr, bool invertp) {
  T substring;

  if (invertp == false) {
    substring = substring2;
    print_genomic(fp,substring2,/*deletion*/NULL,/*deletionlength*/0,/*invertp*/false,
		  queryseq);

  } else {
    substring = substring1;
    print_genomic(fp,substring1,/*deletion*/NULL,/*deletionlength*/0,/*invertp*/true,
		  queryseq);
  }

  fprintf(fp,"\t");

  print_coordinates(fp,substring,chr,invertp);


  fprintf(fp,"\t");
  if (invertp == false) {
    fprintf(fp,"ins:%d..end:%d,sub:%d",nindels,substring->trim_right,substring->nmismatches_bothdiff);
  } else {
    fprintf(fp,"ins:%d..end:%d,sub:%d",nindels,substring->trim_left,substring->nmismatches_bothdiff);
  }
  if (print_ncolordiffs_p) {
    fprintf(fp,"+%d=%d",substring->ncolordiffs,substring->nmismatches_bothdiff+substring->ncolordiffs);
  } else if (print_nsnpdiffs_p) {
    fprintf(fp,"+%d=%d",substring->nmismatches_refdiff - substring->nmismatches_bothdiff,substring->nmismatches_refdiff);
    if (print_snplabels_p && substring->nmismatches_refdiff > substring->nmismatches_bothdiff) {
      print_snp_labels(fp,substring,queryseq);
    }
  }

  return;
}


void
Substring_print_deletion_1 (FILE *fp, T substring1, T substring2, int nindels, 
			    char *deletion, Shortread_T queryseq, char *chr, 
			    bool invertp) {
  T substring;

  if (invertp == false) {
    substring = substring1;
    print_genomic(fp,substring1,deletion,nindels,/*invertp*/false,queryseq);
  } else {
    substring = substring2;
    print_genomic(fp,substring2,deletion,nindels,/*invertp*/true,queryseq);
  }

  fprintf(fp,"\t");

  print_coordinates(fp,substring,chr,invertp);


  fprintf(fp,"\t");
  if (invertp == false) {
    fprintf(fp,"start:%d..del:%d,sub:%d",substring->trim_left,nindels,substring->nmismatches_bothdiff);
  } else {
    fprintf(fp,"start:%d..del:%d,sub:%d",substring->trim_right,nindels,substring->nmismatches_bothdiff);
  }
  if (print_ncolordiffs_p) {
    fprintf(fp,"+%d=%d",substring->ncolordiffs,substring->nmismatches_bothdiff+substring->ncolordiffs);
  } else if (print_nsnpdiffs_p) {
    fprintf(fp,"+%d=%d",substring->nmismatches_refdiff - substring->nmismatches_bothdiff,substring->nmismatches_refdiff);
    if (print_snplabels_p && substring->nmismatches_refdiff > substring->nmismatches_bothdiff) {
      print_snp_labels(fp,substring,queryseq);
    }
  }

  return;
}


void
Substring_print_deletion_2 (FILE *fp, T substring1, T substring2, int nindels, 
			    char *deletion, Shortread_T queryseq, char *chr, 
			    bool invertp) {
  T substring;

  if (invertp == false) {
    substring = substring2;
    print_genomic(fp,substring2,/*deletion*/NULL,/*deletionlength*/0,/*invertp*/false,
		  queryseq);

  } else {
    substring = substring1;
    print_genomic(fp,substring1,/*deletion*/NULL,/*deletionlength*/0,/*invertp*/true,
		  queryseq);
  }

  fprintf(fp,"\t");

  print_coordinates(fp,substring,chr,invertp);

  fprintf(fp,"\t");
  if (invertp == false) {
    fprintf(fp,"del:%d..end:%d,sub:%d",nindels,substring->trim_right,substring->nmismatches_bothdiff);
  } else {
    fprintf(fp,"del:%d..end:%d,sub:%d",nindels,substring->trim_left,substring->nmismatches_bothdiff);
  }
  if (print_ncolordiffs_p) {
    fprintf(fp,"+%d=%d",substring->ncolordiffs,substring->nmismatches_bothdiff+substring->ncolordiffs);
  } else if (print_nsnpdiffs_p) {
    fprintf(fp,"+%d=%d",substring->nmismatches_refdiff - substring->nmismatches_bothdiff,substring->nmismatches_refdiff);
    if (print_snplabels_p && substring->nmismatches_refdiff > substring->nmismatches_bothdiff) {
      print_snp_labels(fp,substring,queryseq);
    }
  }
  
  return;
}


/* This logic used in splice part of SAM_print */
static void
print_splice_distance (FILE *fp, T donor, T acceptor, Genomicpos_T distance, bool sensep) {
  bool normalp = true;

  if (donor == NULL || acceptor == NULL) {
    /* Don't print anything */
  } else if (distance == 0U) {
    fprintf(fp,",splice_type:%s",TRANSLOCATION_TEXT);
  } else {
    if (donor->plusp != acceptor->plusp) {
      fprintf(fp,",splice_type:%s",INVERSION_TEXT);
      normalp = false;
    } else if (donor->plusp == true) {
      if (sensep == true) {
	if (acceptor->genomicstart < donor->genomicstart) {
	  fprintf(fp,",splice_type:%s",SCRAMBLE_TEXT);
	  normalp = false;
	}
      } else {
	if (donor->genomicstart < acceptor->genomicstart) {
	  fprintf(fp,",splice_type:%s",SCRAMBLE_TEXT);
	  normalp = false;
	}
      }
    } else {
      if (sensep == true) {
	if (donor->genomicstart < acceptor->genomicstart) {
	  fprintf(fp,",splice_type:%s",SCRAMBLE_TEXT);
	  normalp = false;
	}
      } else {
	if (acceptor->genomicstart < donor->genomicstart) {
	  fprintf(fp,",splice_type:%s",SCRAMBLE_TEXT);
	  normalp = false;
	}
      }
    }
    if (normalp == true) {
      fprintf(fp,",splice_type:%s",CONSISTENT_TEXT);
    }
    fprintf(fp,",splice_dist:%u",distance);
  }

  return;
}

static void
print_shortexon_splice_distances (FILE *fp, Genomicpos_T distance1, Genomicpos_T distance2) {
  if (distance1 == 0U || distance2 == 0U) {
    /* Skip */
  } else {
    fprintf(fp,",splice_type:%s",CONSISTENT_TEXT);
    fprintf(fp,",splice_dist_1:%u",distance1);
    fprintf(fp,",splice_dist_2:%u",distance2);
  }

  return;
}



void
Substring_print_donor (FILE *fp, T donor, bool sensep, bool invertp, Shortread_T queryseq,
		       IIT_T chromosome_iit, T acceptor, Genomicpos_T chimera_distance) {
  char *chr;
  bool allocp;

#ifdef CHECK_KNOWNI
  Genomicpos_T splicesitepos;
  int *splicesites;
  int nsplicesites;
#endif

  print_genomic(fp,donor,/*deletion*/NULL,/*deletionlength*/0,invertp,queryseq);

  fprintf(fp,"\t");
  chr = IIT_label(chromosome_iit,donor->chrnum,&allocp);
  print_coordinates(fp,donor,chr,invertp);

  /* printf("donor chimera_pos is %d\n",donor->chimera_pos); */
  fprintf(fp,"\t");
  if (sensep == true && invertp == false) {
    fprintf(fp,"start:%d..donor:%.2f",donor->trim_left,donor->chimera_prob);
  } else if (sensep == true && invertp == true) {
    fprintf(fp,"donor:%.2f..end:%d",donor->chimera_prob,donor->trim_left);
  } else if (sensep == false && invertp == false) {
    fprintf(fp,"donor:%.2f..end:%d",donor->chimera_prob,donor->trim_right);
  } else if (sensep == false && invertp == true) {
    fprintf(fp,"start:%d..donor:%.2f",donor->trim_right,donor->chimera_prob);
  }

  fprintf(fp,",sub:%d",donor->nmismatches_bothdiff);
  if (print_ncolordiffs_p) {
    fprintf(fp,"+%d=%d",donor->ncolordiffs,donor->nmismatches_bothdiff+donor->ncolordiffs);
  } else if (print_nsnpdiffs_p) {
    fprintf(fp,"+%d=%d",donor->nmismatches_refdiff - donor->nmismatches_bothdiff,donor->nmismatches_refdiff);
    if (print_snplabels_p && donor->nmismatches_refdiff > donor->nmismatches_bothdiff) {
      print_snp_labels(fp,donor,queryseq);
    }
  }

  if (sensep == true && invertp == false) {
    fprintf(fp,",dir:sense");
  } else if (sensep == true && invertp == true) {
    fprintf(fp,",dir:antisense");
  } else if (sensep == false && invertp == false) {
    fprintf(fp,",dir:antisense");
  } else if (sensep == false && invertp == true) {
    fprintf(fp,",dir:sense");
  }

  if (acceptor != NULL) {
    print_splice_distance(fp,donor,acceptor,chimera_distance,sensep);
  }

#ifdef CHECK_KNOWNI
  if (donor->chimera_knownp == false && splicesites_iit) {
    if (donor->plusp == true) {
      splicesitepos = donor->genomicstart - donor->chroffset + donor->chimera_pos;
    } else {
      splicesitepos = donor->genomicstart - donor->chroffset - donor->chimera_pos;
    }
    splicesites = IIT_get_exact_multiple_with_divno(&nsplicesites,splicesites_iit,
						    splicesites_divint_crosstable[donor->chrnum],
						    splicesitepos,splicesitepos+1U,donor_typeint);
    assert(nsplicesites == 0);
  }
#endif

  if (donor->chimera_knownp && splicesites_iit) {
    print_splicesite_labels(fp,donor,donor_typeint,donor->chimera_pos,/*tag*/"label");
  }

  if (allocp == true) {
    FREE(chr);
  }

  return;
}

void 
Substring_print_acceptor (FILE *fp, T acceptor, bool sensep, bool invertp, Shortread_T queryseq, 
			  IIT_T chromosome_iit, T donor, Genomicpos_T chimera_distance) {
  char *chr;
  bool allocp;

#ifdef CHECK_KNOWNI
  Genomicpos_T splicesitepos;
  int *splicesites;
  int nsplicesites;
#endif

  print_genomic(fp,acceptor,/*deletion*/NULL,/*deletionlength*/0,invertp,queryseq);

  fprintf(fp,"\t");
  chr = IIT_label(chromosome_iit,acceptor->chrnum,&allocp);
  print_coordinates(fp,acceptor,chr,invertp);

  /* printf("acceptor chimera_pos is %d\n",acceptor->chimera_pos); */
  fprintf(fp,"\t");
  if (sensep == true && invertp == false) {
    fprintf(fp,"acceptor:%.2f..end:%d",acceptor->chimera_prob,acceptor->trim_right);
  } else if (sensep == true && invertp == true) {
    fprintf(fp,"start:%d..acceptor:%.2f",acceptor->trim_right,acceptor->chimera_prob);
  } else if (sensep == false && invertp == false) {
    fprintf(fp,"start:%d..acceptor:%.2f",acceptor->trim_left,acceptor->chimera_prob);
  } else if (sensep == false && invertp == true) {
    fprintf(fp,"acceptor:%.2f..end:%d",acceptor->chimera_prob,acceptor->trim_left);
  }

  fprintf(fp,",sub:%d",acceptor->nmismatches_bothdiff);
  if (print_ncolordiffs_p) {
    fprintf(fp,"+%d=%d",acceptor->ncolordiffs,acceptor->nmismatches_bothdiff+acceptor->ncolordiffs);
  } else if (print_nsnpdiffs_p) {
    fprintf(fp,"+%d=%d",acceptor->nmismatches_refdiff - acceptor->nmismatches_bothdiff,acceptor->nmismatches_refdiff);
    if (print_snplabels_p && acceptor->nmismatches_refdiff > acceptor->nmismatches_bothdiff) {
      print_snp_labels(fp,acceptor,queryseq);
    }
  }

  if (sensep == true && invertp == false) {
    fprintf(fp,",dir:sense");
  } else if (sensep == true && invertp == true) {
    fprintf(fp,",dir:antisense");
  } else if (sensep == false && invertp == false) {
    fprintf(fp,",dir:antisense");
  } else if (sensep == false && invertp == true) {
    fprintf(fp,",dir:sense");
  }

  if (donor != NULL) {
    print_splice_distance(fp,donor,acceptor,chimera_distance,sensep);
  }

#ifdef CHECK_KNOWNI
  if (acceptor->chimera_knownp == false && splicesites_iit) {
    if (acceptor->plusp == true) {
      splicesitepos = acceptor->genomicstart - acceptor->chroffset + acceptor->chimera_pos;
    } else {
      splicesitepos = acceptor->genomicstart - acceptor->chroffset - acceptor->chimera_pos;
    }
    splicesites = IIT_get_exact_multiple_with_divno(&nsplicesites,splicesites_iit,
						    splicesites_divint_crosstable[acceptor->chrnum],
						    splicesitepos,splicesitepos+1U,acceptor_typeint);
    assert(nsplicesites == 0);
  }
#endif


  if (acceptor->chimera_knownp && splicesites_iit) {
    print_splicesite_labels(fp,acceptor,acceptor_typeint,
			    acceptor->chimera_pos,/*tag*/"label");
  }


  if (allocp == true) {
    FREE(chr);
  }

  return;
}


void
Substring_print_shortexon (FILE *fp, T shortexon, bool sensep, bool invertp, Shortread_T queryseq,
			   IIT_T chromosome_iit, Genomicpos_T distance1, Genomicpos_T distance2) {
  char *chr;
  bool allocp;

#ifdef CHECK_KNOWNI
  Genomicpos_T splicesitepos;
  int *splicesites;
  int nsplicesites;
#endif

#ifdef CHECK_KNOWNI
    if (shortexon->chimera_knownp == false && splicesites_iit) {
      if (shortexon->plusp == true) {
	splicesitepos = shortexon->genomicstart - shortexon->chroffset + shortexon->chimera_pos;
      } else {
	splicesitepos = shortexon->genomicstart - shortexon->chroffset - shortexon->chimera_pos;
      }
      splicesites = IIT_get_exact_multiple_with_divno(&nsplicesites,splicesites_iit,
						      splicesites_divint_crosstable[shortexon->chrnum],
						      splicesitepos,splicesitepos+1U,acceptor_typeint);
      assert(nsplicesites == 0);
    }

    if (shortexon->chimera_knownp_2 == false && splicesites_iit) {
      if (shortexon->plusp == true) {
	splicesitepos = shortexon->genomicstart - shortexon->chroffset + shortexon->chimera_pos_2;
      } else {
	splicesitepos = shortexon->genomicstart - shortexon->chroffset - shortexon->chimera_pos_2;
      }
      splicesites = IIT_get_exact_multiple_with_divno(&nsplicesites,splicesites_iit,
						      splicesites_divint_crosstable[shortexon->chrnum],
						      splicesitepos,splicesitepos+1U,donor_typeint);
      assert(nsplicesites == 0);
    }
#endif

  print_genomic(fp,shortexon,/*deletion*/NULL,/*deletionlength*/0,invertp,queryseq);

  fprintf(fp,"\t");

  chr = IIT_label(chromosome_iit,shortexon->chrnum,&allocp);
  print_coordinates(fp,shortexon,chr,invertp);

  fprintf(fp,"\t");
  if (sensep == true && invertp == false) {
    fprintf(fp,"acceptor:%.2f..donor:%.2f",shortexon->chimera_prob,shortexon->chimera_prob_2);
  } else if (sensep == true && invertp == true) {
    fprintf(fp,"donor:%.2f..acceptor:%.2f",shortexon->chimera_prob_2,shortexon->chimera_prob);
  } else if (sensep == false && invertp == false) {
    fprintf(fp,"donor:%.2f..acceptor:%.2f",shortexon->chimera_prob_2,shortexon->chimera_prob);
  } else if (sensep == false && invertp == true) {
    fprintf(fp,"acceptor:%.2f..donor:%.2f",shortexon->chimera_prob,shortexon->chimera_prob_2);
  }

  fprintf(fp,",sub:%d",shortexon->nmismatches_bothdiff);
  if (print_ncolordiffs_p) {
    fprintf(fp,"+%d=%d",shortexon->ncolordiffs,shortexon->nmismatches_bothdiff+shortexon->ncolordiffs);
  } else if (print_nsnpdiffs_p) {
    fprintf(fp,"+%d=%d",shortexon->nmismatches_refdiff - shortexon->nmismatches_bothdiff,shortexon->nmismatches_refdiff);
    if (print_snplabels_p && shortexon->nmismatches_refdiff > shortexon->nmismatches_bothdiff) {
      print_snp_labels(fp,shortexon,queryseq);
    }
  }

  if (sensep == true && invertp == false) {
    fprintf(fp,",dir:sense");
    print_shortexon_splice_distances(fp,distance1,distance2);

    if (shortexon->chimera_knownp && splicesites_iit) {
      print_splicesite_labels(fp,shortexon,acceptor_typeint,
			      shortexon->chimera_pos,/*tag*/"label1");
    }
    if (shortexon->chimera_knownp_2 && splicesites_iit) {
      print_splicesite_labels(fp,shortexon,donor_typeint,
			      shortexon->chimera_pos_2,/*tag*/"label2");
    }

  } else if (sensep == true && invertp == true) {
    fprintf(fp,",dir:antisense");
    print_shortexon_splice_distances(fp,distance1,distance2);

    if (shortexon->chimera_knownp_2 && splicesites_iit) {
      print_splicesite_labels(fp,shortexon,donor_typeint,
			      shortexon->chimera_pos_2,/*tag*/"label1");
    }
    if (shortexon->chimera_knownp && splicesites_iit) {
      print_splicesite_labels(fp,shortexon,acceptor_typeint,
			      shortexon->chimera_pos,/*tag*/"label2");
    }

  } else if (sensep == false && invertp == false) {
    fprintf(fp,",dir:antisense");
    print_shortexon_splice_distances(fp,distance1,distance2);



    if (shortexon->chimera_knownp_2 && splicesites_iit) {
      print_splicesite_labels(fp,shortexon,donor_typeint,
			      shortexon->chimera_pos_2,/*tag*/"label1");
    }

    if (shortexon->chimera_knownp && splicesites_iit) {
      print_splicesite_labels(fp,shortexon,acceptor_typeint,
			      shortexon->chimera_pos,/*tag*/"label2");
    }

  } else if (sensep == false && invertp == true) {
    fprintf(fp,",dir:sense");
    print_shortexon_splice_distances(fp,distance1,distance2);
    if (shortexon->chimera_knownp && splicesites_iit) {
      print_splicesite_labels(fp,shortexon,acceptor_typeint,
			      shortexon->chimera_pos,/*tag*/"label1");
    }
    if (shortexon->chimera_knownp_2 && splicesites_iit) {
      print_splicesite_labels(fp,shortexon,donor_typeint,
			      shortexon->chimera_pos_2,/*tag*/"label2");
    }
  }

  if (allocp == true) {
    FREE(chr);
  }

  return;
}


