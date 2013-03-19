static char rcsid[] = "$Id: stage2.c 89122 2013-03-13 22:21:01Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "stage2.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "assert.h"
#include "mem.h"
#include "comp.h"
#include "pair.h"
#include "pairdef.h"
#include "listdef.h"
#include "intlist.h"
#include "diag.h"
#include "oligoindex_hr.h"
#include "genome_hr.h"
#include "complement.h"


/* Tests whether genomicseg == query in convert_to_nucleotides, and
   whether oligoindex_hr gives same results as oligoindex */
/* #define EXTRACT_GENOMICSEG 1 */


/* #define SQUARE 1 */

#define SUFF_PCTCOVERAGE_OLIGOINDEX 0.90

/* #define SUFF_PCTCOVERAGE_STAGE2 0.10 */
#define SUFF_NCOVERED 200
#define SUFF_MAXNCONSECUTIVE 20

#define MAX_NACTIVE 100	/* 100 previously considered too low, but may
			   be okay in conjunction with
			   diagonalization */
#define MAX_GRAND_LOOKBACK 200

/* Penalty for genomic distances */

#define INTRON_PENALTY_UNKNOWN 8
#define INTRON_PENALTY_INCONSISTENT 16

#define NINTRON_PENALTY_MISMATCH 8

#define ENOUGH_CONSECUTIVE 32

#define INFINITE 1000000

/* EQUAL_DISTANCE used to be 3 for PMAP and 6 for GMAP, but that
   allowed indels in repetitive regions.  Now have separate
   variables. */
#ifdef PMAP
#define EQUAL_DISTANCE_FOR_CONSECUTIVE 0
#define EQUAL_DISTANCE_NOT_SPLICING 3
#else
#define EQUAL_DISTANCE_FOR_CONSECUTIVE 0
#define EQUAL_DISTANCE_NOT_SPLICING 9
#endif


#define INTRON_DEFN 9		/* Cannot exceed 9 */
#define NEAR_END_LENGTH 20	/* Determines whether to ignore EXON_DEFN at ends */
#define EXON_DEFN 30
#define MAX_SKIPPED 3

#define SCORE_FOR_RESTRICT 10
/* #define SUFFICIENT_ROOTNLINKS 10  */ /* Setting this too low can slow down program considerably */


#ifdef PMAP
#define SAMPLE_INTERVAL 1
#define NT_PER_MATCH 3
#define CONSEC_POINTS_PER_MATCH 3 /* Possible increase to reward consecutiveness */
#define NONCODON_INDEL_PENALTY 15
#else
#define SAMPLE_INTERVAL 2	/* For cases where adjacentp == false.
				   Means that we can find islands of
				   9-mers */
#define NT_PER_MATCH 1
#define NT_PER_CODON 3
#define CONSEC_POINTS_PER_MATCH 1 /* Possible increase to reward consecutiveness */
#define CONSEC_POINTS_PER_CODON 3 /* Possible increase to reward consecutiveness */
#endif

#define SHIFT_EXTRA 15

static bool splicingp;
static int suboptimal_score_end;
static int suboptimal_score_start;
static Mode_T mode;
static bool snps_p;

void
Stage2_setup (bool splicingp_in, int suboptimal_score_start_in, int suboptimal_score_end_in,
	      Mode_T mode_in, bool snps_p_in) {
  splicingp = splicingp_in;
  suboptimal_score_start = suboptimal_score_start_in;
  suboptimal_score_end = suboptimal_score_end_in;
  mode = mode_in;
  snps_p = snps_p_in;
  return;
}


/* General */
#ifdef DEBUG
#define debug(x) x
#else 
#define debug(x)
#endif

/* Final results of stage 2 */
#ifdef DEBUG0
#define debug0(x) x
#else 
#define debug0(x)
#endif

/* Print all links */
#ifdef DEBUG1
#define debug1(x) x
#else 
#define debug1(x)
#endif

/* For generating a graph */
#ifdef DEBUG3
#define debug3(x) x
#else 
#define debug3(x)
#endif

/* Converting to nucleotides */
#ifdef DEBUG5
#define debug5(x) x
#else
#define debug5(x)
#endif

/* revise_active */
#ifdef DEBUG6
#define debug6(x) x
#else 
#define debug6(x)
#endif

/* Shifted canonical */
#ifdef DEBUG7
#define debug7(x) x
#else 
#define debug7(x)
#endif

/* find_canonical_dinucleotides */
#ifdef DEBUG8
#define debug8(x) x
#else 
#define debug8(x)
#endif

/* Dynamic programming */
/* Can also define debug9(x) as: if (querypos == XX) {x;} */
#ifdef DEBUG9
#define debug9(x) x
#else 
#define debug9(x)
#endif

/* Grand winner */
#ifdef DEBUG10
#define debug10(x) x
#else 
#define debug10(x)
#endif

/* Multiple alignments */
#ifdef DEBUG11
#define debug11(x) x
#else 
#define debug11(x)
#endif



#define T Stage2_T
struct T {
  List_T path;
  int nfwdintrons;
  int ncanonical;
  int nnoncanonical;
  int cdna_direction;

  double diag_runtime;
  double align_runtime;
  int indexsize;
};

List_T
Stage2_path (T this) {
  return this->path;
}

int
Stage2_ncanonical (T this) {
  return this->ncanonical;
}

int
Stage2_nnoncanonical (T this) {
  return this->nnoncanonical;
}

int
Stage2_cdna_direction (T this) {
  return this->cdna_direction;
}

double
Stage2_diag_runtime (T this) {
  return this->diag_runtime;
}

double
Stage2_align_runtime (T this) {
  return this->align_runtime;
}

int
Stage2_indexsize (T this) {
  return this->indexsize;
}

int
Stage2_pathlength (T this) {
  return List_length(this->path);
}

#if 0
static T
Stage2_new (List_T path, int nfwdintrons, int nrevintrons, int nunkintrons,
	    int cdna_direction, int indexsize_nt) {
  T new = (T) MALLOC(sizeof(*new));

  new->path = path;
  if (cdna_direction > 0) {
    new->ncanonical = nfwdintrons;
    new->nnoncanonical = nrevintrons + nunkintrons;
  } else if (cdna_direction < 0) {
    new->ncanonical = nrevintrons;
    new->nnoncanonical = nfwdintrons + nunkintrons;
  } else {
    abort();
  }
  new->cdna_direction = cdna_direction;
  new->indexsize = indexsize_nt;

  return new;
}
#endif	    

void
Stage2_free (T *old) {
  if (*old) {
    /* These cells are taken from Pairpool_T and should not be freed. */
    /* List_free(&(*old)->path); */
    FREE(*old);
  }
  return;
}


/************************************************************************/

typedef struct Link_T *Link_T;
struct Link_T {
  int fwd_consecutive;
  int fwd_rootposition;		/* Last branch, or equivalently beginning of consecutive run */
  /*int fwd_rootnlinks;*/		/* Number of links in last branch */
  int fwd_score;

  int fwd_pos;
  int fwd_hit;
#ifdef USE_SUBOPTIMAL_STARTS
  int fwd_initposition;
  Intlist_T fwd_pos_ties;
  Intlist_T fwd_hit_ties;
#endif

  int fwd_intronnfwd;
  int fwd_intronnrev;
  int fwd_intronnunk;
#ifndef PMAP
  int rev_consecutive;
  int rev_rootposition;		/* Last branch, or equivalently beginning of consecutive run */
  /*int rev_rootnlinks;*/		/* Number of links in last branch */
  int rev_score;

  int rev_pos;
  int rev_hit;
#ifdef USE_SUBOPTIMAL_STARTS
  int rev_initposition;
  Intlist_T rev_pos_ties;
  Intlist_T rev_hit_ties;
#endif

  int rev_intronnfwd;
  int rev_intronnrev;
  int rev_intronnunk;
#endif
};


/* lengths2 is has length1 entries.  Note that lengths2 may have
   negative entries */
static struct Link_T **
Linkmatrix_1d_new (int length1, int *lengths2, int totallength) {
  struct Link_T **links;
  int i;
  
  links = (struct Link_T **) CALLOC(length1,sizeof(struct Link_T *));
  links[0] = (struct Link_T *) CALLOC(totallength,sizeof(struct Link_T));
  for (i = 1; i < length1; i++) {
    if (lengths2[i-1] < 0) {
      links[i] = links[i-1];
    } else {
      links[i] = &(links[i-1][lengths2[i-1]]);
    }
  }
  return links;
}

#ifdef USE_SUBOPTIMAL_STARTS
static void
Linkmatrix_gc (struct Link_T **links, int length1, int *lengths2) {
  int i, j;

  for (i = 0; i < length1; i++) {
    for (j = 0; j < lengths2[i]; j++) {
      if (links[i][j].fwd_pos_ties) {
	Intlist_free(&(links[i][j].fwd_pos_ties));
	Intlist_free(&(links[i][j].fwd_hit_ties));
      }
#ifndef PMAP
      if (links[i][j].rev_pos_ties) {
	Intlist_free(&(links[i][j].rev_pos_ties));
	Intlist_free(&(links[i][j].rev_hit_ties));
      }
#endif
    }
  }

  return;
}
#endif


static void
Linkmatrix_1d_free (struct Link_T ***links) {
  FREE((*links)[0]);
  FREE(*links);
  return;
}


static struct Link_T **
Linkmatrix_2d_new (int length1, int *lengths2) {
  struct Link_T **links;
  int i;
  
  links = (struct Link_T **) CALLOC(length1,sizeof(struct Link_T *));
  for (i = 0; i < length1; i++) {
    if (lengths2[i] <= 0) {
      links[i] = (struct Link_T *) NULL;
    } else {
      links[i] = (struct Link_T *) CALLOC(lengths2[i],sizeof(struct Link_T));
    }
  }
  return links;
}

static void
Linkmatrix_2d_free (struct Link_T ***links, int length1) {
  int i;

  for (i = 0; i < length1; i++) {
    if ((*links)[i]) {
      FREE((*links)[i]);
    }
  }
  FREE(*links);
  return;
}



#ifdef DEBUG1
#ifdef PMAP
/* For PMAP, indexsize is in aa */
static void
Linkmatrix_print_fwd (struct Link_T **links, unsigned int **mappings, int *npositions,
		      char *queryseq_ptr, int indexsize) {
  int i, j, lastpos;
  char *oligo;
  Intlist_T p, q;

  oligo = (char *) CALLOC(indexsize+1,sizeof(char));
  lastpos = length1 - indexsize;

  for (i = 0; i <= lastpos; i++) {
    strncpy(oligo,&(queryseq_ptr[i]),indexsize);
    printf("Querypos %d (%s, %d positions):",i,oligo,npositions[i]);
    for (j = 0; j < npositions[i]; j++) {
#ifdef USE_SUBOPTIMAL_STARTS
      printf(" %d.%u:%d(%d,%d)[%u..%u]",
	     j,mappings[i][j],links[i][j].fwd_score,
	     links[i][j].fwd_pos,links[i][j].fwd_hit,links[i][j].fwd_initposition,links[i][j].fwd_rootposition);
      for (p = links[i][j].fwd_pos_ties, q = links[i][j].fwd_hit_ties; p != NULL; p = Intlist_next(p), q = Intlist_next(q)) {
	printf(" TIE at %d,%d",Intlist_head(p),Intlist_head(q));
      }
#else
      printf(" %d.%u:%d(%d,%d)[%u]",
	     j,mappings[i][j],links[i][j].fwd_score,
	     links[i][j].fwd_pos,links[i][j].fwd_hit,links[i][j].fwd_rootposition);
#endif
    }
    printf("\n");
  }
  printf("\n");

  FREE(oligo);
  return;
}

#else

static void
Linkmatrix_print_both (struct Link_T **links, unsigned int **mappings, int length1,
		       int *npositions, char *queryseq_ptr, int indexsize) {
  int i, j;
  char *oligo;
#ifdef USE_SUBOPTIMAL_STARTS
  Intlist_T p, q;
#endif

  oligo = (char *) CALLOC(indexsize+1,sizeof(char));
  for (i = 0; i <= length1-indexsize; i++) {
    strncpy(oligo,&(queryseq_ptr[i]),indexsize);
    printf("Querypos %d (%s, %d positions):",i,oligo,npositions[i]);
    for (j = 0; j < npositions[i]; j++) {
#ifdef USE_SUBOPTIMAL_STARTS
      printf(" %d.%u:%d(%d,%d)[%u..%u]-%d(%d,%d)[%u..%u]",
	     j,mappings[i][j],links[i][j].fwd_score,
	     links[i][j].fwd_pos,links[i][j].fwd_hit,links[i][j].fwd_initposition,links[i][j].fwd_rootposition,
	     links[i][j].rev_score,
	     links[i][j].rev_pos,links[i][j].rev_hit,links[i][j].rev_initposition,links[i][j].rev_rootposition);
      for (p = links[i][j].fwd_pos_ties, q = links[i][j].fwd_hit_ties; p != NULL; p = Intlist_next(p), q = Intlist_next(q)) {
	printf(" FWD_TIE at %d,%d",Intlist_head(p),Intlist_head(q));
      }
      for (p = links[i][j].rev_pos_ties, q = links[i][j].rev_hit_ties; p != NULL; p = Intlist_next(p), q = Intlist_next(q)) {
	printf(" REV_TIE at %d,%d",Intlist_head(p),Intlist_head(q));
      }
#else
      printf(" %d.%u:%d(%d,%d)[%u]-%d(%d,%d)[%u]",
	     j,mappings[i][j],links[i][j].fwd_score,
	     links[i][j].fwd_pos,links[i][j].fwd_hit,links[i][j].fwd_rootposition,
	     links[i][j].rev_score,
	     links[i][j].rev_pos,links[i][j].rev_hit,links[i][j].rev_rootposition);
#endif
    }
    printf("\n");
  }
  printf("\n");

  FREE(oligo);
  return;
}
#endif
#endif

static void
mappings_dump_R (unsigned int **mappings, int *npositions, int length1,
		 int **active, int *firstactive, int indexsize, char *varname) {
  int querypos;
  int lastpos, hit;
  bool printp = false;

  lastpos = length1 - indexsize;
  printf("%s <- matrix(c(\n",varname);
  for (querypos = 0; querypos < lastpos; querypos++) {
    if (firstactive) {
      if (mappings[querypos] != NULL) {
	hit = firstactive[querypos];
	while (hit != -1) {
	  /* Last elt is for score */
	  if (printp == false) {
	    printp = true;
	  } else {
	    printf(",\n");
	  }
	  printf("%d,%d,%d,%d",querypos,mappings[querypos][hit],
		 hit,active[querypos][hit]);
	  hit = active[querypos][hit];
	}
      }
    } else {
      for (hit = 0; hit < npositions[querypos]; hit++) {
	if (printp == false) {
	  printp = true;
	} else {
	  printf(",\n");
	}
	printf("%d,%d,%d",querypos,mappings[querypos][hit],hit);
      }
    }
  }
  printf("),ncol=2,byrow=T)\n");

  return;
}
    

static void
best_path_dump_R (struct Link_T **links, unsigned int **mappings,
		  int querypos, int hit, bool fwdp, char *varname) {
  unsigned int position;
  int prev_querypos, prevhit, save_querypos, savehit;
  bool printp = false;

  save_querypos = querypos;
  savehit = hit;

  printf("%s <- matrix(c(\n",varname);
  prev_querypos = querypos+1;
  while (querypos >= 0) {
    position = mappings[querypos][hit];

    if (printp == false) {
      printp = true;
    } else {
      printf(",\n");
    }
    printf("%d,%d",querypos,position);

    prev_querypos = querypos;
    prevhit = hit;
    if (fwdp) {
      querypos = links[prev_querypos][prevhit].fwd_pos;
      hit = links[prev_querypos][prevhit].fwd_hit;
#ifndef PMAP
    } else {
      querypos = links[prev_querypos][prevhit].rev_pos;
      hit = links[prev_querypos][prevhit].rev_hit;
#endif
    }
  }
  printf("),ncol=2,byrow=T)\n");

  querypos = save_querypos;
  hit = savehit;

  printp = false;
  printf("%s <- matrix(c(\n","scores");
  prev_querypos = querypos+1;
  while (querypos >= 0) {
    position = mappings[querypos][hit];

    if (printp == false) {
      printp = true;
    } else {
      printf(",\n");
    }
    if (fwdp == true) {
      printf("%d,%d",querypos,links[querypos][hit].fwd_score);
#ifndef PMAP
    } else {
      printf("%d,%d",querypos,links[querypos][hit].rev_score);
#endif
    }

    prev_querypos = querypos;
    prevhit = hit;
    if (fwdp) {
      querypos = links[prev_querypos][prevhit].fwd_pos;
      hit = links[prev_querypos][prevhit].fwd_hit;
#ifndef PMAP
    } else {
      querypos = links[prev_querypos][prevhit].rev_pos;
      hit = links[prev_querypos][prevhit].rev_hit;
#endif
    }
  }
  printf("),ncol=2,byrow=T)\n");

  return;
}

static void
active_bounds_dump_R (unsigned int *minactive, unsigned int *maxactive,
		      int querylength) {
  int querypos;
  bool printp = false;

  printf("querypos <- 0:%d\n",querylength-1);
  printf("%s <- c(\n","minactive");
  for (querypos = 0; querypos < querylength; querypos++) {
    if (printp == false) {
      printp = true;
    } else {
      printf(",\n");
    }
    printf("%d",minactive[querypos]);
  }
  printf(")\n");

  printp = false;
  printf("%s <- c(\n","maxactive");
  for (querypos = 0; querypos < querylength; querypos++) {
    if (printp == false) {
      printp = true;
    } else {
      printf(",\n");
    }
    printf("%d",maxactive[querypos]);
  }
  printf(")\n");

  return;
}


#define TEN_THOUSAND 10000.0
#define HUNDRED_THOUSAND 100000.0
#define ONE_MILLION 1000000.0


#if 0
/* Not used anymore */
static int
gendist_penalty (unsigned int diffdistance, int querydistance) {
  int penalty;

  /* Add 1 to get effect of querydistance multiplier below */
  penalty = diffdistance/HUNDRED_THOUSAND + 1;
  if (diffdistance <= EQUAL_DISTANCE_NOT_SPLICING) {
    return 0;			/* was penalty */
#if 0
  } else if (querydistance < 8) {
    return penalty;
#endif
  } else {
    penalty *= (querydistance/8) * (querydistance/8);
    return penalty;
  }
}
#endif


#define diffdist_penalty_nosplicing(x) (x + 1)
#define diffdist_penalty_splicing(x) (x/TEN_THOUSAND + 1)


/* querydistance already has indexsize_nt subtracted */
static int
querydist_penalty (int querydistance) {
#ifdef PMAP
  return querydistance/2;
#else
  return querydistance/8;
#endif
}


#if 0
static int
diffdist_mismatch_penalty (int diffdistance, bool crossspeciesp) {
  int penalty;

  if (crossspeciesp == true) {
    if (diffdistance <= 24) {
      return (diffdistance+7)/2;
    } else if (diffdistance <= 40) {
      return (diffdistance+7)/2;
    } else if (diffdistance <= 56) {
      return (diffdistance+7);
    } else {
      return (diffdistance+7)*3/2;
    }
  } else {
#ifdef PMAP
    if (diffdistance <= 24) {
      penalty = (diffdistance+7)/4;
    } else if (diffdistance <= 40) {
      penalty = (diffdistance+7)/4;
    } else if (diffdistance <= 56) {
      penalty = (diffdistance+7)/2;
    } else {
      penalty = (diffdistance+7)*3/4;
    }
    return penalty;
#else
    if (diffdistance <= 24) {
      return (diffdistance+7)/8;
    } else if (diffdistance <= 40) {
      return (diffdistance+7)/4;
    } else if (diffdistance <= 56) {
      return (diffdistance+7)/2;
    } else {
      return (diffdistance+7)*3/4;
    }
#endif
  }
}
#endif


/************************************************************************
 *  Procedures for finding canonical introns quickly
 ************************************************************************/

#ifdef DEBUG8

static void
print_last_dinucl (int *last_dinucl, int genomiclength) {
  int pos;

  for (pos = 0; pos < genomiclength - 3 + SHIFT_EXTRA; pos++) {
    printf("%d %d\n",pos,last_dinucl[pos]);
  }
  printf("\n");

  return;
}

#endif


#if 0
/* A specialized Boyer-Moore algorithm */
static void
find_canonical_dinucleotides (int *lastGT, int *lastAG, 
#ifndef PMAP
			      int *lastCT, int *lastAC,
#endif
			      char *genomicuc_ptr, int genomiclength) {
  int pos, GT = -1, AG = -1;
#ifndef PMAP
  int CT = -1, AC = -1;
#endif
  char c1, c2;

  /* printf("%.*s\n",100,genomicuc_ptr); */

  lastGT[0] = lastGT[1] = lastGT[2] = -1;
  lastAG[0] = lastAG[1] = lastAG[2] = -1;
#ifndef PMAP
  lastCT[0] = lastCT[1] = lastCT[2] = -1;
  lastAC[0] = lastAC[1] = lastAC[2] = -1;
#endif

  pos = 1;			/* Not 0, so we don't overwrite initial -1's */
  c1 = genomicuc_ptr[pos+1];
  while (pos <= genomiclength-4) {
    if ((c2 = genomicuc_ptr[pos+2]) == 'T') {
      if (c1 == 'G') {
	GT = pos;
#ifndef PMAP
      } else if (c1 == 'C') {
	CT = pos;
#endif
      }
      
      lastGT[pos] = lastGT[pos+1] = GT;
      lastAG[pos+3] = lastAG[pos+4] = AG;
#ifndef PMAP
      lastCT[pos] = lastCT[pos+1] = CT;
      lastAC[pos+3] = lastAC[pos+4] = AC;
#endif

      pos += 2;
      c1 = genomicuc_ptr[pos+1];

    } else {
      if (c2 == 'G') {
	if (c1 == 'A') {
	  AG = pos + 3;
	}
#ifndef PMAP
      } else if (c2 == 'C') {
	if (c1 == 'A') {
	  AC = pos + 3;
	}
#endif
      }
	
      lastGT[pos] = GT;
      lastAG[pos+3] = AG;
#ifndef PMAP
      lastCT[pos] = CT;
      lastAC[pos+3] = AC;
#endif

      pos++;
      c1 = c2;
    }
  }

  /* pos is now either genomiclength-3 (one more dinucleotide) or genomiclength-2 (none) */
  if (pos == genomiclength-3) {
    if ((c2 = genomicuc_ptr[pos+2]) == 'T') {
      if (c1 == 'G') {
	GT = pos;
#ifndef PMAP
      } else if (c1 == 'C') {
	CT = pos;
#endif
      }
    } else {
      if (c2 == 'G') {
	if (c1 == 'A') {
	  AG = pos + 3;
	}
#ifndef PMAP
      } else if (c2 == 'C') {
	if (c1 == 'A') {
	  AC = pos + 3;
	}
#endif
      }
    }
    
    lastGT[pos] = GT;
    lastAG[pos+3] = AG;
#ifndef PMAP
    lastCT[pos] = CT;
    lastAC[pos+3] = AC;
#endif
  }

  /* Fill in rest */
  while (pos < genomiclength - 3 + SHIFT_EXTRA) {
    lastGT[pos] = GT;
    lastAG[pos+3] = AG;
#ifndef PMAP
    lastCT[pos] = CT;
    lastAC[pos+3] = AC;
#endif
    pos++;
  }


  debug8(printf("Old GT\n"));
  debug8(print_last_dinucl(lastGT,genomiclength));
  debug8(printf("Old AG\n"));
  debug8(print_last_dinucl(lastAG,genomiclength));
  debug8(printf("Old AC\n"));
  debug8(print_last_dinucl(lastAC,genomiclength));
  debug8(printf("Old CT\n"));
  debug8(print_last_dinucl(lastCT,genomiclength));

  return;
}
#endif


#if 0
static void
find_canonical_dinucleotides_hr (int *lastGT, int *lastAG, 
#ifndef PMAP
				 int *lastCT, int *lastAC,
#endif
				 Genomicpos_T genomicstart, bool plusp,
				 int genomiclength) {

  Genome_last_donor_positions(lastGT,genomicstart,/*margin5*/3,/*margin3*/3,
			      genomiclength,plusp);
  Genome_last_acceptor_positions(lastAG,genomicstart,/*margin5*/3,/*margin3*/3,
				 genomiclength,plusp);
#ifndef PMAP
  Genome_last_antidonor_positions(lastAC,genomicstart,/*margin5*/3,/*margin3*/3,
				  genomiclength,plusp);
  Genome_last_antiacceptor_positions(lastCT,genomicstart,/*margin5*/3,/*margin3*/3,
				     genomiclength,plusp);
#endif


  debug8(printf("New GT\n"));
  debug8(print_last_dinucl(lastGT,genomiclength));
  debug8(printf("New AG\n"));
  debug8(print_last_dinucl(lastAG,genomiclength));
  debug8(printf("New AC\n"));
  debug8(print_last_dinucl(lastAC,genomiclength));
  debug8(printf("New CT\n"));
  debug8(print_last_dinucl(lastCT,genomiclength));

  return;
}
#endif


#if 0
static void
check_canonical_dinucleotides_hr (int *lastGT, int *lastAG, 
#ifndef PMAP
				  int *lastCT, int *lastAC,
#endif
				  Genomicpos_T genomicstart, Genomicpos_T genomicend, bool plusp,
				  int genomiclength) {

  int pos;
  int lastpos1, lastpos2;

  /* printf("genomicstart %u, genomicend %u, genomiclength %d\n",genomicstart,genomicend,genomiclength); */

#if 1
  pos = 6;
  while (pos <= genomiclength-4) {
    lastpos1 = lastGT[pos];
    lastpos2 = Genome_prev_donor_position(pos,/*pos5 (wrong)*/1,chroffset,chrhigh,plusp);
    /* printf("at pos %d, lastpos1 %d, lastpos2 %d\n",pos,lastpos1,lastpos2); */
    if (lastpos1 != lastpos2 && lastpos1 != -1) {
      printf("donor plusp %d: at pos %d, lastpos1 %d != lastpos2 %d\n",plusp,pos,lastpos1,lastpos2);
      abort();
    }
    pos++;
  }
#endif

#if 1
  pos = 6;
  while (pos <= genomiclength-4) {
    lastpos1 = lastAG[pos];
    lastpos2 = Genome_prev_acceptor_position(pos,/*pos5 (wrong)*/1,chroffset,chrhigh,plusp);
    /* printf("at pos %d, lastpos1 %d, lastpos2 %d\n",pos,lastpos1,lastpos2); */
    if (lastpos1 != lastpos2 && lastpos1 != -1) {
      printf("acceptor plusp %d: at pos %d, lastpos1 %d != lastpos2 %d\n",plusp,pos,lastpos1,lastpos2);
      abort();
    }
    pos++;
  }
#endif

#ifndef PMAP
#if 1
  pos = 6;
  while (pos <= genomiclength-4) {
    lastpos1 = lastAC[pos];
    lastpos2 = Genome_prev_antidonor_position(pos,/*pos5 (wrong)*/1,chroffset,chrhigh,plusp);
    /* printf("at pos %d, lastpos1 %d, lastpos2 %d\n",pos,lastpos1,lastpos2); */
    if (lastpos1 != lastpos2 && lastpos1 != -1) {
      printf("antidonor plusp %d: at pos %d, lastpos1 %d != lastpos2 %d\n",plusp,pos,lastpos1,lastpos2);
      abort();
    }
    pos++;
  }
#endif

#if 1
  pos = 6;
  while (pos <= genomiclength-4) {
    lastpos1 = lastCT[pos];
    lastpos2 = Genome_prev_antiacceptor_position(pos,/*pos5 (wrong)*/1,chroffset,chrhigh,plusp);
    /* printf("at pos %d, lastpos1 %d, lastpos2 %d\n",pos,lastpos1,lastpos2); */
    if (lastpos1 != lastpos2 && lastpos1 != -1) {
      printf("antiacceptor plusp %d: at pos %d, lastpos1 %d != lastpos2 %d\n",plusp,pos,lastpos1,lastpos2);
      abort();
    }
    pos++;
  }
#endif
#endif

  return;
}
#endif



#if 0
static int
get_last (int *last_dinucl, int pos) {
  int lastpos;
  int i;

  if (pos < 0) {
    return -1;
  }

#if 0
  i = pos;
  while (i >= 0 && last_dinucl[i] == 0) {
    i--;
  }
  if (i < 0) {
    lastpos = -1;
  } else {
    lastpos = last_dinucl[i];
  }
#else
  i = pos;
  while (last_dinucl[i] == 0) {
    i--;
  }
  lastpos = last_dinucl[i];
#endif


#if 0
  /* Fill in region for future use */
  while (i <= pos) {
    last_dinucl[i] = lastpos;
    i++;
  }
#else
  /* Fill in position future use */
  last_dinucl[pos] = lastpos;
#endif

  return lastpos;
}
#endif


#if 0
/* Need this procedure because we are skipping some oligomers */
static bool
find_shifted_canonical_memoize (int leftpos, int rightpos,
				int querydistance, int **last_leftdi, int **last_rightdi,
				int (*genome_left_position)(int, Genomicpos_T, Genomicpos_T, int, bool),
				int (*genome_right_position)(int, Genomicpos_T, Genomicpos_T, int, bool),
				Genomicpos_T genomicstart, Genomicpos_T genomicend, bool plusp,
				bool skip_repetitive_p) {
  int leftdi, rightdi;
  int shift, leftmiss, rightmiss;

  /* leftpos = prevposition + querydistance + indexsize_nt - 1; */
  /* rightpos = position; */

  debug7(printf("Looking for shifted canonical at leftpos %d to rightpos %d at genomic %u..%u\n",
		leftpos,rightpos,genomicstart,genomicend));

#if 0
  /* previously checked against genomiclength */
  if (leftpos > genomiclength || rightpos > genomiclength) {
    return false;
  }
#else
  /* Checking just before call to genome_right_position */
#endif

  assert(leftpos < rightpos);
  if (skip_repetitive_p == false) {
    if ((Genomicpos_T) rightpos >= genomicend) {
      return false;
    }

    if (*last_leftdi == NULL) {
      debug7(printf("Allocating %d ints for leftdi\n",genomicend-genomicstart+SHIFT_EXTRA+1));
      *last_leftdi = (int *) CALLOC(genomicend-genomicstart+SHIFT_EXTRA+1,sizeof(int));
    }
    if ((*last_leftdi)[leftpos] == 0) {
      (*last_leftdi)[leftpos] = (*genome_left_position)(leftpos,genomicstart,genomicend,/*pos5*/3,plusp);
    }

    if (*last_rightdi == NULL) {
      debug7(printf("Allocating %d ints for rightdi\n",genomicend-genomicstart+SHIFT_EXTRA+1));
      *last_rightdi = (int *) CALLOC(genomicend-genomicstart+SHIFT_EXTRA+1,sizeof(int));
    }
    if ((*last_rightdi)[rightpos] == 0) {
      (*last_rightdi)[rightpos] = (*genome_right_position)(rightpos,genomicstart,genomicend,/*pos5*/3,plusp);
    }

    assert((*last_leftdi)[leftpos] != 0);
    assert((*last_rightdi)[rightpos] != 0);

    return (leftpos == (*last_leftdi)[leftpos] && rightpos == (*last_rightdi)[rightpos]);
  }

  /* Allow canonical to be to right of match */
  leftpos += SHIFT_EXTRA;
  rightpos += SHIFT_EXTRA;
  debug7(printf("after shift, leftpos = %u, rightpos = %u\n",leftpos,rightpos));

  shift = 0;
  while (shift <= querydistance + SHIFT_EXTRA + SHIFT_EXTRA) {

    if (leftpos < 0) {
      return false;
#if 0
    } else if (rightpos < 0) {
      /* Shouldn't need to check if leftpos >= 0 and rightpos >= leftpos, in the other two conditions) */
      return false;
#endif
    } else if ((Genomicpos_T) rightpos >= genomicend) {
      return false;
    } else if (leftpos > rightpos) {
      return false;
    }
    assert(rightpos >= 0);

    
    if ((*last_leftdi)[leftpos] == 0) {
      (*last_leftdi)[leftpos] = (*genome_left_position)(leftpos,genomicstart,genomicend,/*pos5*/3,plusp);
    }
    assert((*last_leftdi)[leftpos] != 0);
    if ((leftdi = (*last_leftdi)[leftpos]) < 0) {
      debug7(printf("\n"));
      return false;
    } else {
      leftmiss = leftpos - leftdi;
    }

    if ((*last_rightdi)[rightpos] == 0) {
      (*last_rightdi)[rightpos] = (*genome_right_position)(rightpos,genomicstart,genomicend,/*pos5*/3,plusp);
    }
    assert((*last_rightdi)[rightpos] != 0);
    if ((rightdi = (*last_rightdi)[rightpos]) < 0) {
      debug7(printf("\n"));
      return false;
    } else {
      rightmiss = rightpos - rightdi;
    }

    debug7(printf("shift %d/left %d (miss %d)/right %d (miss %d)\n",shift,leftpos,leftmiss,rightpos,rightmiss));
    if (leftmiss == rightmiss) {  /* was leftmiss == 0 && rightmiss == 0, which doesn't allow for a shift */
      debug7(printf(" => Success\n\n"));
      return true;
    } else if (leftmiss >= rightmiss) {
      shift += leftmiss;
      leftpos -= leftmiss;
      rightpos -= leftmiss;
    } else {
      shift += rightmiss;
      leftpos -= rightmiss;
      rightpos -= rightmiss;
    }
  }

  debug7(printf("\n"));
  return false;
}
#endif


/* assert(chrstart < chrend) */
/* For plus, chrinit = chrstart, chrterm = chrend.  For minus, chrinit = (chrhigh - chroffset) - chrend, chrterm = (chrhigh - chroffset) - chrstart. */

/* Need this procedure because we are skipping some oligomers */
static bool
find_shifted_canonical (Genomicpos_T leftpos, Genomicpos_T rightpos, int querydistance, 
			Genomicpos_T (*genome_left_position)(Genomicpos_T, Genomicpos_T, Genomicpos_T, Genomicpos_T, bool),
			Genomicpos_T (*genome_right_position)(Genomicpos_T, Genomicpos_T, Genomicpos_T, Genomicpos_T, bool),
			Genomicpos_T chroffset, Genomicpos_T chrhigh, bool plusp, bool skip_repetitive_p) {
  Genomicpos_T leftdi, rightdi;
  Genomicpos_T last_leftpos, last_rightpos;
  int shift, leftmiss, rightmiss;
  Genomicpos_T left_chrbound, right_chrbound;
  
  /* leftpos = prevposition + querydistance + indexsize_nt - 1; */
  /* rightpos = position; */

  debug7(printf("Looking for shifted canonical at leftpos %u to rightpos %u, chrhigh %u\n",leftpos,rightpos,chrhigh));

#if 0
  /* previously checked against genomiclength */
  if (leftpos > genomiclength || rightpos > genomiclength) {
    return false;
  }
#else
  /* Checking just before call to genome_right_position */
#endif

  if (leftpos >= rightpos) {
    debug7(printf("leftpos %u >= rightpos %u, so returning false\n",leftpos,rightpos));
    return false;
  }

  if (leftpos < 100) {
    left_chrbound = 0;
  } else {
    left_chrbound = leftpos - 100;
  }

  if (rightpos < 100) {
    right_chrbound = 0;
  } else {
    right_chrbound = rightpos - 100;
  }

#if 0
  if (skip_repetitive_p == false) {
#if 0    
    if ((Genomicpos_T) rightpos >= chrterm) {
      return false;
    }
#endif

    last_leftpos = (*genome_left_position)(leftpos,left_chrbound,chroffset,chrhigh,plusp);
    last_rightpos = (*genome_right_position)(rightpos,right_chrbound,chroffset,chrhigh,plusp);
    debug7(printf("last_leftpos %u, last_rightpos %u\n",last_leftpos,last_rightpos));

    debug7(printf("skip_repetitive_p == false, so returning %u == %u && %u == %u\n",
		  leftpos,last_leftpos,rightpos,last_rightpos));
    return (leftpos == last_leftpos && rightpos == last_rightpos);
  }
#endif

  /* Allow canonical to be to right of match */
  leftpos += SHIFT_EXTRA;
  rightpos += SHIFT_EXTRA;
  debug7(printf("after shift, leftpos = %u, rightpos = %u\n",leftpos,rightpos));

  shift = 0;
  while (shift <= querydistance + SHIFT_EXTRA + SHIFT_EXTRA) {

#if 0
    if (leftpos < chrinit) {
      return false;
    } else if (rightpos < 0) {
      /* Shouldn't need to check if leftpos >= 0 and rightpos >= leftpos, in the other two conditions) */
      return false;
    } else if ((Genomicpos_T) rightpos >= chrterm) {
      return false;
    }
#endif
    if (leftpos > rightpos) {
      return false;
    }

    last_leftpos = (*genome_left_position)(leftpos,left_chrbound,chroffset,chrhigh,plusp);
    debug7(printf("last_leftpos %u\n",last_leftpos));
    assert(last_leftpos != 0U);
    if ((leftdi = last_leftpos) == -1U) {
      debug7(printf("\n"));
      return false;
    } else {
      leftmiss = (int) (leftpos - leftdi);
    }

    last_rightpos = (*genome_right_position)(rightpos,right_chrbound,chroffset,chrhigh,plusp);
    debug7(printf("last_rightpos %u\n",last_rightpos));
    assert(last_rightpos != 0U);
    if ((rightdi = last_rightpos) == -1U) {
      debug7(printf("\n"));
      return false;
    } else {
      rightmiss = (int) (rightpos - rightdi);
    }

    debug7(printf("shift %d/left %d (miss %d)/right %d (miss %d)\n",shift,leftpos,leftmiss,rightpos,rightmiss));
    if (leftmiss == rightmiss) {  /* was leftmiss == 0 && rightmiss == 0, which doesn't allow for a shift */
      debug7(printf(" => Success at %u..%u (fwd) or %u..%u (rev)\n\n",
		    leftpos-leftmiss+/*onebasedp*/1U,rightpos-rightmiss+/*onebasedp*/1U,
		    chrhigh-chroffset-(leftpos-leftmiss),chrhigh-chroffset-(rightpos-rightmiss)));
      return true;
    } else if (leftmiss >= rightmiss) {
      shift += leftmiss;
      leftpos -= leftmiss;
      rightpos -= leftmiss;
    } else {
      shift += rightmiss;
      leftpos -= rightmiss;
      rightpos -= rightmiss;
    }
  }

  debug7(printf("\n"));
  return false;
}


/* For PMAP, indexsize is in aa. */
static void
score_querypos_general (Link_T currlink, int querypos,
			int querystart, int queryend, unsigned int position,
			struct Link_T **links, unsigned int **mappings,
			int **active, int *firstactive,
#ifdef USE_SUBOPTIMAL_STARTS
			int *fwd_initposition_bestscore, int *fwd_initposition_bestpos, int *fwd_initposition_besthit,
			int *rev_initposition_bestscore, int *rev_initposition_bestpos, int *rev_initposition_besthit,
#endif
#ifdef BLUETARP
			int *nactive,
			int grand_fwd_querypos, int grand_fwd_hit,
#ifndef PMAP
			int grand_rev_querypos, int grand_rev_hit,
#endif
#endif
#if 0
			int *lastGT, int *lastAG,
#ifndef PMAP
			int *lastCT, int *lastAC, 
#endif
#endif
			Genomicpos_T chroffset, Genomicpos_T chrhigh, bool plusp,
			int indexsize, Intlist_T processed, int sufflookback, int nsufflookback, int maxintronlen, 
			bool localp, bool skip_repetitive_p, bool use_shifted_canonical_p) {
  Link_T prevlink;
#ifdef USE_SUBOPTIMAL_STARTS
  Intlist_T fwd_prevposes = NULL, fwd_prevhits = NULL, fwd_scores = NULL,
    rev_prevposes = NULL, rev_prevhits = NULL, rev_scores = NULL, p, q, r;
  int best_fwd_initposition = -1;
#else
  Intlist_T p;
#endif

  int best_fwd_consecutive = indexsize*NT_PER_MATCH;
  int best_fwd_score = 0, fwd_score;
  int best_fwd_prevpos = -1, best_fwd_prevhit = -1;
  int best_fwd_intronnfwd = 0, best_fwd_intronnrev = 0, best_fwd_intronnunk = 0;
  int best_fwd_rootposition = position;
  /* int best_fwd_rootnlinks = 1; */
#ifndef PMAP
  int best_rev_consecutive = indexsize*NT_PER_MATCH;
  int best_rev_score = 0, rev_score;
  int best_rev_prevpos = -1, best_rev_prevhit = -1;
  int best_rev_intronnfwd = 0, best_rev_intronnrev = 0, best_rev_intronnunk = 0;
  int best_rev_rootposition = position;
#ifdef USE_SUBOPTIMAL_STARTS
  int best_rev_initposition = -1;
#endif
  /* int best_rev_rootnlinks = 1; */
#endif
  bool adjacentp = false, donep;
  int prev_querypos, prevhit;
  unsigned int prevposition, gendistance;
  int querydistance, diffdistance, lookback, nlookback, nseen, indexsize_nt;
  int leftpos, last_leftpos, rightpos;
  bool canonicalp, last_canonicalp = false;
  int canonicalsgn;
  int enough_consecutive;
  bool near_end_p;

#ifdef PMAP
  indexsize_nt = indexsize*3;
#else
  indexsize_nt = indexsize;
#endif

  enough_consecutive = 32;

  if (querypos < querystart + NEAR_END_LENGTH) {
    near_end_p = true;
  } else if (querypos > queryend - NEAR_END_LENGTH) {
    near_end_p = true;
  } else {
    near_end_p = false;
  }


  /* A. Evaluate adjacent position (at last one processed) */
  if (processed != NULL) {
    prev_querypos = Intlist_head(processed);
#ifdef PMAP
    querydistance = (querypos - prev_querypos)*3;
#else
    querydistance = querypos - prev_querypos;
#endif
    prevhit = firstactive[prev_querypos];
    while (prevhit != -1 && (prevposition = mappings[prev_querypos][prevhit]) + querydistance <= position) {
      if (prevposition + querydistance == position) {
	prevlink = &(links[prev_querypos][prevhit]);
	best_fwd_consecutive = prevlink->fwd_consecutive + querydistance;
	best_fwd_rootposition = prevlink->fwd_rootposition;
#ifdef USE_SUBOPTIMAL_STARTS
	best_fwd_initposition = prevlink->fwd_initposition;
#endif
	/* best_fwd_rootnlinks = prevlink->fwd_rootnlinks + 1; */
	best_fwd_score = prevlink->fwd_score + querydistance;
#ifndef PMAP
	best_rev_consecutive = prevlink->rev_consecutive + querydistance;
	best_rev_rootposition = prevlink->rev_rootposition;
#ifdef USE_SUBOPTIMAL_STARTS
	best_rev_initposition = prevlink->rev_initposition;
#endif
	/* best_rev_rootnlinks = prevlink->rev_rootnlinks + 1; */
	best_rev_score = prevlink->rev_score + querydistance;
#endif

	best_fwd_prevpos = 
#ifndef PMAP
	  best_rev_prevpos = 
#endif
	  prev_querypos;
	best_fwd_prevhit = 
#ifndef PMAP
	  best_rev_prevhit = 
#endif
	  prevhit;
	best_fwd_intronnfwd = prevlink->fwd_intronnfwd;
	best_fwd_intronnrev = prevlink->fwd_intronnrev;
	best_fwd_intronnunk = prevlink->fwd_intronnunk;
#ifndef PMAP
	best_rev_intronnfwd = prevlink->rev_intronnfwd;
	best_rev_intronnrev = prevlink->rev_intronnrev;
	best_rev_intronnunk = prevlink->rev_intronnunk;
#endif
	adjacentp = true;
#ifdef PMAP
	debug9(printf("\tA. Adjacent qpos %d,%d at %ux%d (scores = %d -> %d, consec = %d (from %d), intr = %d-%d-%d)\n",
		      prev_querypos,prevhit,prevposition,active[prev_querypos][prevhit],prevlink->fwd_score,
		      best_fwd_score,best_fwd_consecutive,best_fwd_rootposition,
		      best_fwd_intronnfwd,best_fwd_intronnrev,best_fwd_intronnunk));
#else
	debug9(printf("\tA. Adjacent qpos %d,%d at %ux%d (scores = %d/%d -> %d/%d, consec = %d/%d (from %d/%d), intr = %d-%d-%d/%d-%d-%d)\n",
		      prev_querypos,prevhit,prevposition,active[prev_querypos][prevhit],prevlink->fwd_score,prevlink->rev_score,
		      best_fwd_score,best_rev_score,best_fwd_consecutive,best_rev_consecutive,best_fwd_rootposition,best_rev_rootposition,
		      best_fwd_intronnfwd,best_fwd_intronnrev,best_fwd_intronnunk,
		      best_rev_intronnfwd,best_rev_intronnrev,best_rev_intronnunk));
#endif
	prevhit = -1;		/* Exit loop */
      } else {
	prevhit = active[prev_querypos][prevhit];
      }
    }

#ifdef USE_SUBOPTIMAL_STARTS
    /* Just push one for all adjacent positions */
    if (best_fwd_prevhit >= 0) {
      fwd_prevposes = Intlist_push(fwd_prevposes,best_fwd_prevpos);
      fwd_prevhits = Intlist_push(fwd_prevhits,best_fwd_prevhit);
      fwd_scores = Intlist_push(fwd_scores,best_fwd_score);
    }
#ifndef PMAP
    if (best_rev_prevhit >= 0) {
      rev_prevposes = Intlist_push(rev_prevposes,best_rev_prevpos);
      rev_prevhits = Intlist_push(rev_prevhits,best_rev_prevhit);
      rev_scores = Intlist_push(rev_scores,best_rev_score);
    }
#endif
#endif
  }


#if 0
  /* AA. Evaluate querypos - 3 (for wobble) */
  if (adjacentp == false && querypos >= 3) {
    prevhit = firstactive[querypos-3];
    while (prevhit != -1 && (prevposition = mappings[querypos-3][prevhit]) + NT_PER_CODON <= position) {
      if (prevposition + NT_PER_CODON == position) {
	prevlink = &(links[querypos-3][prevhit]);
	best_fwd_consecutive = prevlink->fwd_consecutive + NT_PER_CODON;
	best_fwd_rootposition = prevlink->fwd_rootposition;
#ifdef USE_SUBOPTIMAL_STARTS
	best_fwd_initposition = prevlink->fwd_initposition;
#endif
	/* best_fwd_rootnlinks = prevlink->fwd_rootnlinks + 1; */
	best_fwd_score = prevlink->fwd_score + CONSEC_POINTS_PER_CODON;
	best_rev_consecutive = prevlink->rev_consecutive + NT_PER_CODON;
	best_rev_rootposition = prevlink->rev_rootposition;
#ifdef USE_SUBOPTIMAL_STARTS
	best_rev_initposition = prevlink->rev_initposition;
#endif
	/* best_rev_rootnlinks = prevlink->rev_rootnlinks + 1; */
	best_rev_score = prevlink->rev_score + CONSEC_POINTS_PER_CODON;

	best_fwd_prevpos = best_rev_prevpos = querypos-3;
	best_fwd_prevhit = best_rev_prevhit = prevhit;
	best_fwd_intronnfwd = prevlink->fwd_intronnfwd;
	best_fwd_intronnrev = prevlink->fwd_intronnrev;
	best_fwd_intronnunk = prevlink->fwd_intronnunk;
	best_rev_intronnfwd = prevlink->rev_intronnfwd;
	best_rev_intronnrev = prevlink->rev_intronnrev;
	best_rev_intronnunk = prevlink->rev_intronnunk;
	adjacentp = true;

	debug9(printf("\tAA. Codon-adjacent qpos %d,%d at %u (scores = %d/%d -> %d/%d, consec = %d/%d (from %d/%d), intr = %d-%d-%d/%d-%d-%d)\n",
		      querypos-3,prevhit,position - NT_PER_CODON,prevlink->fwd_score,prevlink->rev_score,
		      best_fwd_score,best_rev_score,best_fwd_consecutive,best_rev_consecutive,best_fwd_rootposition,best_rev_rootposition,
		      best_fwd_intronnfwd,best_fwd_intronnrev,best_fwd_intronnunk,
		      best_rev_intronnfwd,best_rev_intronnrev,best_rev_intronnunk));
	prevhit = -1; /* Exit loop */
      } else {
	prevhit = active[querypos-3][prevhit];
      }
    }
  }
#endif


  if (best_fwd_consecutive < enough_consecutive
#ifndef PMAP
      || best_rev_consecutive < enough_consecutive
#endif
      ) {

#ifdef BLUETARP
    /* It looks like case D can cover cases B and C */

    /* B. Evaluate for intron (at querypos - indexsize).  
       Don't update best_consecutive, because a reasonable alternative 
       might be a mismatch, regardless of the quality of the intron. */

#ifdef PMAP
    /* Use indexsize, not indexsize_nt */
    intronpos_lower = querypos - indexsize - 1;
    intronpos_upper = querypos - indexsize;
#else
    intronpos_lower = intronpos_upper = querypos - indexsize_nt;
#endif
    if (intronpos_lower < 0) {
      intronpos_lower = 0;
    }

    for (prev_querypos = intronpos_upper; prev_querypos >= intronpos_lower; --prev_querypos) {
      if (0 && skip_repetitive_p && nactive[prev_querypos] > MAX_NACTIVE) { /* Turn off.  If too many active, then firstactive will be -1 */
	debug9(printf("Not evaluating for intron because nactive at prev_querypos %d is %d > max %d\n",
		      prev_querypos,nactive[prev_querypos],MAX_NACTIVE));
	prevhit = -1;
      } else {
#ifdef PMAP
	querydistance = (querypos - prev_querypos)*3 - indexsize_nt
#else
	querydistance = (querypos - prev_querypos) - indexsize_nt;
#endif
	prevhit = firstactive[prev_querypos];
      }

      debug9(printf("Evaluating for intron at prev_querypos %d, prevhit %d\n",prev_querypos,prevhit));
      while (prevhit != -1 && (prevposition = mappings[prev_querypos][prevhit]) + indexsize_nt <= position) {
	prevlink = &(links[prev_querypos][prevhit]);
	if (1 /* want to find canonicalsgn for short exons */ || prevlink->fwd_consecutive > EXON_DEFN
#ifndef PMAP
	    || prevlink->rev_consecutive > EXON_DEFN
#endif
	    ) {
	  gendistance = position - prevposition - indexsize_nt;
	  diffdistance = abs(gendistance - querydistance);
	
	  if (splicingp == false) {
	    canonicalsgn = 0;
	  } else if (diffdistance > maxintronlen) {
	    canonicalsgn = 0;
	  } else if (find_shifted_canonical(/*leftpos*/prevposition + querydistance + indexsize_nt - 1,
					    /*rightpos*/position,querydistance,
					    /* &lastGT,&lastAG, */
					    Genome_prev_donor_position,Genome_prev_acceptor_position,
					    chroffset,chrhigh,plusp,skip_repetitive_p) == true) {
	    canonicalsgn = +1;
#ifndef PMAP
	  } else if (find_shifted_canonical(/*leftpos*/prevposition + querydistance + indexsize_nt - 1,
					    /*rightpos*/position,querydistance,
					    /* &lastCT,&lastAC, */
					    Genome_prev_antiacceptor_position,Genome_prev_antidonor_position,
					    chroffset,chrhigh,plusp,skip_repetitive_p) == true) {
	    canonicalsgn = -1;
#endif
	  } else {
	    canonicalsgn = 0;
	  }

	  /* Need to both compensate for the oligomer lost across
	     the intron and penalize for each additional intron */
#ifdef PMAP
	  fwd_score = prevlink->fwd_score + querydistance + indexsize_nt;
#else
	  fwd_score = prevlink->fwd_score + querydistance + indexsize_nt;
	  rev_score = prevlink->rev_score + querydistance + indexsize_nt;
#endif
	
	  /* Penalty for introns. */
	  if (diffdistance > maxintronlen) {
	    fwd_gendistance_penalty = 
#ifndef PMAP
	      rev_gendistance_penalty = 
#endif
	      INFINITE;
	  } else {
	    if (canonicalsgn == 0) {
	      /* Extra penalty for known non-canonical intron */
	      fwd_gendistance_penalty = 
#ifndef PMAP
		rev_gendistance_penalty = 
#endif
		gendist_penalty(diffdistance,querydistance/*-indexsize_nt*/) + INTRON_PENALTY_UNKNOWN;
	    } else if (canonicalsgn == +1) {
	      /* Don't penalize for intron */
	      fwd_gendistance_penalty = gendist_penalty(diffdistance,querydistance/*-indexsize_nt*/) /* - INTRON_REWARD_CONSISTENT*/;
#ifndef PMAP
	      rev_gendistance_penalty = gendist_penalty(diffdistance,querydistance/*-indexsize_nt*/) + INTRON_PENALTY_INCONSISTENT;
#endif
#ifndef PMAP
	    } else if (canonicalsgn == -1) {
	      fwd_gendistance_penalty = gendist_penalty(diffdistance,querydistance) + INTRON_PENALTY_INCONSISTENT;
	      /* Don't penalize for intron */
	      rev_gendistance_penalty = gendist_penalty(diffdistance,querydistance) /* - INTRON_REWARD_CONSISTENT*/;
#endif
	    }
	  }
	  fwd_score -= fwd_gendistance_penalty;
#ifndef PMAP
	  rev_score -= rev_gendistance_penalty;
#endif

	  /* Allow ties, which should favor shorter intron */
	  if (fwd_score >= best_fwd_score) {
	    best_fwd_consecutive = indexsize_nt;
	    best_fwd_rootposition = position;
	    /* best_fwd_rootnlinks = 1; */
	    /* best_fwd_score = fwd_score; -- Put below so debug9 works */
	    best_fwd_prevpos = prev_querypos;
	    best_fwd_prevhit = prevhit;
	    best_fwd_intronnfwd = prevlink->fwd_intronnfwd;
	    best_fwd_intronnrev = prevlink->fwd_intronnrev;
	    best_fwd_intronnunk = prevlink->fwd_intronnunk;
	    switch (canonicalsgn) {
	    case 1: best_fwd_intronnfwd++; break;
	    case -1: best_fwd_intronnrev++; break;
	    case 0: best_fwd_intronnunk++; break;
	    }
	  }

#ifndef PMAP
	  /* Allow ties, which should favor shorter intron */
	  if (rev_score >= best_rev_score) {
	    best_rev_consecutive = indexsize_nt;
	    best_rev_rootposition = position;
	    /* best_rev_rootnlinks = 1; */
	    /* best_rev_score = rev_score; -- Put below so debug9 works */
	    best_rev_prevpos = prev_querypos;
	    best_rev_prevhit = prevhit;
	    best_rev_intronnfwd = prevlink->rev_intronnfwd;
	    best_rev_intronnrev = prevlink->rev_intronnrev;
	    best_rev_intronnunk = prevlink->rev_intronnunk;
	    switch (canonicalsgn) {
	    case 1: best_rev_intronnfwd++; break;
	    case -1: best_rev_intronnrev++; break;
	    case 0: best_rev_intronnunk++; break;
	    }
	  }
#endif

#ifdef PMAP
	  debug9(printf("\tB. Intron qpos %d,%d at %ux%d (scores = %d -> %d, consec = %d (from %d), intr = %d-%d-%d, gendist %u, pen %d)",
			prev_querypos,prevhit,prevposition,active[prev_querypos][prevhit],
			prevlink->fwd_score,fwd_score,prevlink->fwd_consecutive,prevlink->fwd_rootposition,
			best_fwd_intronnfwd,best_fwd_intronnrev,best_fwd_intronnunk,
			gendistance,fwd_gendistance_penalty));
#else
	  debug9(printf("\tB. Intron qpos %d,%d at %ux%d (scores = %d/%d -> %d/%d, consec = %d/%d (from %d/%d), intr = %d-%d-%d/%d-%d-%d, gendist %u, pen %d/%d)",
			prev_querypos,prevhit,prevposition,active[prev_querypos][prevhit],
			prevlink->fwd_score,prevlink->rev_score,fwd_score,rev_score,
			prevlink->fwd_consecutive,prevlink->rev_consecutive,prevlink->fwd_rootposition,prevlink->rev_rootposition,
			best_fwd_intronnfwd,best_fwd_intronnrev,best_fwd_intronnunk,
			best_rev_intronnfwd,best_rev_intronnrev,best_rev_intronnunk,
			gendistance,fwd_gendistance_penalty,rev_gendistance_penalty));
#endif
	  debug9({
	    printf(" =>");
	    if (fwd_score > best_fwd_score) {
	      printf(" Best fwd at %d",fwd_score);
	    } else {
	      printf(" Fwd loses to %d",best_fwd_score);
	    }
#ifndef PMAP
	    if (rev_score > best_rev_score) {
	      printf(" Best rev at %d",rev_score);
	    } else {
	      printf(" Rev loses to %d\n",best_rev_score);
	    }
#endif
	    printf("\n");
	  });
	}

	/* Put these down here so debug9 works */
	if (fwd_score > best_fwd_score) {
	  best_fwd_score = fwd_score;
	}
#ifndef PMAP
	if (rev_score > best_rev_score) {
	  best_rev_score = rev_score;
	}
#endif
	prevhit = active[prev_querypos][prevhit];
      }
    }

#endif /* bluetarp */


#ifdef BLUETARP
    /* It looks like case D can cover cases B and C */

    /* C (fwd). Look back at grand querypos and hit */
    debug9(printf("\tGrand fwd querypos is %d,%d at position %u\n",
		  grand_fwd_querypos,grand_fwd_hit,grand_fwd_hit >= 0 ? mappings[grand_fwd_querypos][grand_fwd_hit] : 0));
    if (grand_fwd_querypos > 0 && grand_fwd_querypos + indexsize_nt <= querypos &&
	(prevposition = mappings[grand_fwd_querypos][grand_fwd_hit]) + indexsize_nt <= position) {
      prev_querypos = grand_fwd_querypos;
      prevhit = grand_fwd_hit;
#ifdef PMAP
      querydistance = (querypos - prev_querypos) - indexsize_nt;
#else
      querydistance = querypos - prev_querypos - indexsize_nt;
#endif
      if (querydistance < MAX_GRAND_LOOKBACK) {
	prevlink = &(links[prev_querypos][prevhit]);
	if (1 /* want to find canonicalsgn for short exons */ || prevlink->fwd_consecutive > EXON_DEFN) {
	  gendistance = position - prevposition - indexsize_nt;
	  diffdistance = abs(gendistance - querydistance);

	  if (splicingp == false) {
	    canonicalsgn = 0;
	  } else if (diffdistance <= EQUAL_DISTANCE_NOT_SPLICING) {
	    canonicalsgn = 9;
	  } else if (diffdistance > maxintronlen) {
	    canonicalsgn = 0;
	  } else if (find_shifted_canonical(/*leftpos*/prevposition + querydistance + indexsize_nt - 1,
					    /*rightpos*/position,querydistance,
					    /* &lastGT,&lastAG, */
					    Genome_prev_donor_position,Genome_prev_acceptor_position,
					    chroffset,chrhigh,plusp,skip_repetitive_p) == true) {
	    canonicalsgn = +1;
	  } else {
	    canonicalsgn = 0;
	  }

	  fwd_score = prevlink->fwd_score;
	  if (diffdistance > maxintronlen) {
	    fwd_score -= INFINITE;
	  } else if (diffdistance <= EQUAL_DISTANCE_FOR_CONSECUTIVE) {
	    fwd_score = prevlink->fwd_score + CONSEC_POINTS_PER_MATCH;
	  } else if (diffdistance < INTRON_DEFN) {
	    /* fwd_score -= querydist_mismatch_penalty(querydistance); */
#ifdef PMAP
	    if (diffdistance % 3 != 0) {
	      fwd_score -= NONCODON_INDEL_PENALTY;
	    }
#endif
	  } else if (canonicalsgn == +1) {
	    /* Don't penalize for intron.  Don't penalize for querydistance to grand querypos. */
	    /* fwd_score -= gendist_penalty(diffdistance,querydistance); */
	    if (splicingp == true) {
	      fwd_score -= diffdist_penalty_splicing(diffdistance);
	    } else {
	      fwd_score -= diffdist_penalty_nosplicing(diffdistance);
	    }

	  } else {
	    /* fwd_score -= querydist_mismatch_penalty(querydistance); */
	    /* fwd_gendistance_penalty = gendist_penalty(diffdistance,querydistance); */
	    if (splicingp == true) {
	      fwd_gendistance_penalty = diffdist_penalty_splicing(diffdistance);
	    } else {
	      fwd_gendistance_penalty = diffdist_penalty_nosplicing(diffdistance);
	    }
	    fwd_score -= fwd_gendistance_penalty;
	    fwd_score -= NINTRON_PENALTY_MISMATCH;
	  }

	  debug9(printf("\tC+. Fwd grand qpos %d,%d at %ux%d (score = %d -> %d, consec = %d (from %d), intr = %d-%d-%d, gendist %u, querydist %d, canonicalsgn %d)",
			prev_querypos,prevhit,prevposition,active[prev_querypos][prevhit],
			prevlink->fwd_score,fwd_score,prevlink->fwd_consecutive,prevlink->fwd_rootposition,
			best_fwd_intronnfwd,best_fwd_intronnrev,best_fwd_intronnunk,
			gendistance,querydistance,canonicalsgn));
	    
	  /* Allow ties, which should favor shorter intron */
	  if (fwd_score >= best_fwd_score) {
	    if (diffdistance <= EQUAL_DISTANCE_FOR_CONSECUTIVE) {
	      best_fwd_consecutive = prevlink->fwd_consecutive + NT_PER_MATCH;
	      best_fwd_rootposition = prevlink->fwd_rootposition;
	      /* best_fwd_rootnlinks = prevlink->rootnlinks + 1; */
	    } else {
	      best_fwd_consecutive = indexsize_nt;
	      best_fwd_rootposition = position;
	      /* best_fwd_rootnlinks = 1; */
	    }
	    best_fwd_score = fwd_score;
	    best_fwd_prevpos = prev_querypos;
	    best_fwd_prevhit = prevhit;
	    best_fwd_intronnfwd = prevlink->fwd_intronnfwd;
	    best_fwd_intronnrev = prevlink->fwd_intronnrev;
	    best_fwd_intronnunk = prevlink->fwd_intronnunk;
	    switch (canonicalsgn) {
	    case 1: best_fwd_intronnfwd++; break;
	    case -1: best_fwd_intronnrev++; break;
	    case 0: best_fwd_intronnunk++; break;
	    }
	    debug9(printf(" => Best fwd at %d (consec = %d)\n",fwd_score,best_fwd_consecutive));
	  } else {
	    debug9(printf(" => Loses to %d\n",best_fwd_score));
	  }
	}
      }
    }

#ifndef PMAP
    /* C (rev). Look back at grand querypos and hit */
    debug9(printf("\tGrand rev querypos is %d,%d at position %u\n",
		  grand_rev_querypos,grand_rev_hit,grand_rev_hit >= 0 ? mappings[grand_rev_querypos][grand_rev_hit] : 0));
    if (grand_rev_querypos > 0 && grand_rev_querypos + indexsize_nt <= querypos &&
	(prevposition = mappings[grand_rev_querypos][grand_rev_hit]) + indexsize_nt <= position) {
      prev_querypos = grand_rev_querypos;
      prevhit = grand_rev_hit;
#ifdef PMAP
      querydistance = (querypos - prev_querypos)*3 - indexsize_nt;
#else
      querydistance = querypos - prev_querypos - indexsize_nt;
#endif
      if (querydistance < MAX_GRAND_LOOKBACK) {
	prevlink = &(links[prev_querypos][prevhit]);
	if (1 /* want to find canonicalsgn for short exons */ || prevlink->rev_consecutive > EXON_DEFN) {
	  gendistance = position - prevposition - indexsize_nt;
	  diffdistance = abs(gendistance - querydistance);

	  if (splicingp == false) {
	    canonicalsgn = 0;
	  } else if (diffdistance <= EQUAL_DISTANCE_NOT_SPLICING) {
	    canonicalsgn = 9;
	  } else if (diffdistance > maxintronlen) {
	    canonicalsgn = 0;
	  } else if (find_shifted_canonical(/*leftpos*/prevposition + querydistance + indexsize_nt - 1,
					    /*rightpos*/position,querydistance,
					    /* &lastCT,&lastAC, */
					    Genome_prev_antiacceptor_position,Genome_prev_antidonor_position,
					    chroffset,chrhigh,plusp,skip_repetitive_p) == true) {
	    canonicalsgn = -1;
	  } else {
	    canonicalsgn = 0;
	  }

	  rev_score = prevlink->rev_score;
	  if (diffdistance > maxintronlen) {
	    rev_score -= INFINITE;
	  } else if (diffdistance <= EQUAL_DISTANCE_FOR_CONSECUTIVE) {
	    rev_score = prevlink->rev_score + CONSEC_POINTS_PER_MATCH;
	  } else if (diffdistance < INTRON_DEFN) {
	    /* rev_score -= querydist_mismatch_penalty(querydistance); */
	  } else if (canonicalsgn == -1) {
	    /* Don't penalize for intron.  Don't penalize for querydistance back to grand querypos. */
	    /* rev_score -= gendist_penalty(diffdistance,querydistance); */
	    if (splicingp == true) {
	      rev_score -= diffdist_penalty_splicing(diffdistance);
	    } else {
	      rev_score -= diffdist_penalty_nosplicing(diffdistance);
	    }

	  } else {
	    /* rev_score -= querydist_mismatch_penalty(querydistance); */
	    /* rev_gendistance_penalty = gendist_penalty(diffdistance,querydistance); */
	    if (splicingp == true) {
	      rev_gendistance_penalty = diffdist_penalty_splicing(diffdistance);
	    } else {
	      rev_gendistance_penalty = diffdist_penalty_nosplicing(diffdistance);
	    }
	    rev_score -= rev_gendistance_penalty;
	    rev_score -= NINTRON_PENALTY_MISMATCH;
	  }

	  debug9(printf("\tC-. Rev grand qpos %d,%d at %ux%d (score = %d -> %d, consec = %d (from %d), intr = %d-%d-%d, gendist %u, querydist %d, canonicalsgn %d)",
			prev_querypos,prevhit,prevposition,active[prev_querypos][prevhit],
			prevlink->rev_score,rev_score,prevlink->rev_consecutive,prevlink->rev_rootposition,
			best_rev_intronnrev,best_rev_intronnrev,best_rev_intronnunk,
			gendistance,querydistance,canonicalsgn));
	    
	  /* Allow ties, which should favor shorter intron */
	  if (rev_score >= best_rev_score) {
	    if (diffdistance <= EQUAL_DISTANCE_FOR_CONSECUTIVE) {
	      best_rev_consecutive = prevlink->rev_consecutive + NT_PER_MATCH;
	      best_rev_rootposition = prevlink->rev_rootposition;
	      /* best_rev_rootnlinks = prevlink->rev_rootnlinks + 1; */
	    } else {
	      best_rev_consecutive = indexsize_nt;
	      best_rev_rootposition = position;
	      /* best_rev_rootnlinks = 1; */
	    }
	    best_rev_score = rev_score;
	    best_rev_prevpos = prev_querypos;
	    best_rev_prevhit = prevhit;
	    best_rev_intronnfwd = prevlink->rev_intronnfwd;
	    best_rev_intronnrev = prevlink->rev_intronnrev;
	    best_rev_intronnunk = prevlink->rev_intronnunk;
	    switch (canonicalsgn) {
	    case 1: best_rev_intronnfwd++; break;
	    case -1: best_rev_intronnrev++; break;
	    case 0: best_rev_intronnunk++; break;
	    }
	    debug9(printf(" => Best rev at %d (consec = %d)\n",rev_score,best_rev_consecutive));
	  } else {
	    debug9(printf(" => Loses to %d\n",best_rev_score));
	  }
	}
      }
    }
#endif

#endif /* bluetarp */


    /* D (fwd). Evaluate for mismatches (all other previous querypos) */
    /* Set parameters */
    if (adjacentp == true) {
      /* Look not so far back */
      nlookback = 1;
      lookback = sufflookback/2;
    } else {
      /* Look farther back */
      nlookback = nsufflookback;
      lookback = sufflookback;
    }

    last_leftpos = 0;		/* if use_shifted_canonical_p is true */
    rightpos = position;	/* if use_shifted_canonical_p is true */

    donep = false;
    for (p = processed, nseen = 0;
	 p != NULL && best_fwd_consecutive < enough_consecutive && donep == false;
	 p = Intlist_next(p), nseen++) {
				  
      prev_querypos = Intlist_head(p);
#ifdef PMAP
      querydistance = (querypos - prev_querypos)*3 - indexsize_nt;
#else
      querydistance = querypos - prev_querypos - indexsize_nt;
#endif

      if (nseen > nlookback && querydistance > lookback) {
	donep = true;
      }

      prevhit = firstactive[prev_querypos];
      while (prevhit != -1 && (prevposition = mappings[prev_querypos][prevhit]) + indexsize_nt <= position) {
	/* printf("fwd: prevposition %u, prevhit %d\n",prevposition,prevhit); */
	prevlink = &(links[prev_querypos][prevhit]);

	gendistance = position - prevposition - indexsize_nt;
	diffdistance = abs(gendistance - querydistance);

	if (diffdistance < maxintronlen) {
	  if (diffdistance <= EQUAL_DISTANCE_NOT_SPLICING) {
	    canonicalsgn = 9;
	    fwd_score = prevlink->fwd_score + CONSEC_POINTS_PER_MATCH;
#ifdef PMAP
	    if (diffdistance % 3 != 0) {
	      fwd_score -= NONCODON_INDEL_PENALTY;
	    }
#endif
	  } else if (near_end_p == false && prevlink->fwd_consecutive < EXON_DEFN) {
	    canonicalsgn = 0;
	    if (splicingp == true) {
	      fwd_score = prevlink->fwd_score - diffdist_penalty_splicing(diffdistance) - querydist_penalty(querydistance) - NINTRON_PENALTY_MISMATCH;
	    } else {
	      fwd_score = prevlink->fwd_score - diffdist_penalty_nosplicing(diffdistance) - querydist_penalty(querydistance) - NINTRON_PENALTY_MISMATCH;
	    }

	  } else {
	    /* Called only near end */
	    if (splicingp == false) {
	      canonicalsgn = 0;
	      fwd_score = prevlink->fwd_score - diffdist_penalty_nosplicing(diffdistance) - querydist_penalty(querydistance);

	    } else if (use_shifted_canonical_p == true) {
	      leftpos = prevposition + querydistance + indexsize_nt - 1;
	      /* printf("leftpos %d, last_leftpos %d, rightpos %d\n",leftpos,last_leftpos,rightpos); */
	      if (leftpos == last_leftpos) {
		canonicalp = last_canonicalp;
	      } else {
		debug7(printf("Calling find_shift_canonical fwd\n"));
		canonicalp = find_shifted_canonical(leftpos,rightpos,querydistance,
						    /* &lastGT,&lastAG, */
						    Genome_prev_donor_position,Genome_prev_acceptor_position,
						    chroffset,chrhigh,plusp,skip_repetitive_p);
		last_leftpos = leftpos;
		last_canonicalp = canonicalp;
	      }
	      if (canonicalp == true) {
		canonicalsgn = +1;
		fwd_score = prevlink->fwd_score - diffdist_penalty_splicing(diffdistance) - querydist_penalty(querydistance);
	      } else {
		canonicalsgn = 0;
		fwd_score = prevlink->fwd_score - diffdist_penalty_splicing(diffdistance) - querydist_penalty(querydistance) - NINTRON_PENALTY_MISMATCH;
	      }

	    } else {
	      canonicalsgn = +1;
	      fwd_score = prevlink->fwd_score - diffdist_penalty_splicing(diffdistance) - querydist_penalty(querydistance);
	    }
#if 0
	    if (prevlink->fwd_consecutive < 20) {
	      /* Intended to avoid too many jumps, but misses short introns */
	      fwd_score -= diffdist_penalty(diffdistance);
	    }
#endif
	  }


	  debug9(printf("\tD+. Fwd mismatch qpos %d,%d at %ux%d (score = %d -> %d, consec = %d (from %d), intr = %d-%d-%d, gendist %u, querydist %d, canonicalsgn %d)",
			prev_querypos,prevhit,prevposition,active[prev_querypos][prevhit],
			prevlink->fwd_score,fwd_score,prevlink->fwd_consecutive,prevlink->fwd_rootposition,
			best_fwd_intronnfwd,best_fwd_intronnrev,best_fwd_intronnunk,
			gendistance,querydistance,canonicalsgn));
	    
#ifdef USE_SUBOPTIMAL_STARTS
	  fwd_prevposes = Intlist_push(fwd_prevposes,prev_querypos);
	  fwd_prevhits = Intlist_push(fwd_prevhits,prevhit);
	  fwd_scores = Intlist_push(fwd_scores,fwd_score);
#endif

	  /* Allow ties, which should favor shorter intron */
	  if (fwd_score >= best_fwd_score) {
	    if (diffdistance <= EQUAL_DISTANCE_FOR_CONSECUTIVE) {
	      best_fwd_consecutive = prevlink->fwd_consecutive + (querydistance + indexsize_nt);
	      best_fwd_rootposition = prevlink->fwd_rootposition;
	      /* best_fwd_rootnlinks = prevlink->fwd_rootnlinks + 1; */
	    } else {
	      best_fwd_consecutive = 0;
	      best_fwd_rootposition = position;
	      /* best_fwd_rootnlinks = 1; */
	    }
#ifdef USE_SUBOPTIMAL_STARTS
	    best_fwd_initposition = prevlink->fwd_initposition;
#endif
	    best_fwd_score = fwd_score;
	    best_fwd_prevpos = prev_querypos;
	    best_fwd_prevhit = prevhit;
	    best_fwd_intronnfwd = prevlink->fwd_intronnfwd;
	    best_fwd_intronnrev = prevlink->fwd_intronnrev;
	    best_fwd_intronnunk = prevlink->fwd_intronnunk;
	    switch (canonicalsgn) {
	    case 1: best_fwd_intronnfwd++; break;
	    case 0: best_fwd_intronnunk++; break;
	    }
	    debug9(printf(" => Best fwd at %d (consec = %d)\n",fwd_score,best_fwd_consecutive));
	  } else {
	    debug9(printf(" => Loses to %d\n",best_fwd_score));
	  }
	}

	prevhit = active[prev_querypos][prevhit];
      }
    }


#ifndef PMAP
    /* D (rev). Evaluate for mismatches (all other previous querypos) */
    /* Parameters already set in D (fwd). */

    if (splicingp == false || use_shifted_canonical_p == false) {
      /* The rev scores are identical to the fwd scores */

    } else {
      last_leftpos = 0;
      rightpos = position;
      donep = false;
      for (p = processed, nseen = 0;
	   p != NULL && best_rev_consecutive < enough_consecutive && donep == false;
	   p = Intlist_next(p), nseen++) {

	prev_querypos = Intlist_head(p);
#ifdef PMAP
	querydistance = (querypos - prev_querypos)*3 - indexsize_nt;
#else
	querydistance = querypos - prev_querypos - indexsize_nt;
#endif

	if (nseen > nlookback && querydistance > lookback) {
	  donep = true;
	}

	prevhit = firstactive[prev_querypos];
	while (prevhit != -1 && (prevposition = mappings[prev_querypos][prevhit]) + indexsize_nt <= position) {
	  /* printf("rev: prevposition %u, prevhit %d\n",prevposition,prevhit); */
	  prevlink = &(links[prev_querypos][prevhit]);

	  gendistance = position - prevposition - indexsize_nt;
	  diffdistance = abs(gendistance - querydistance);

	  if (diffdistance < maxintronlen) {
	    if (diffdistance <= EQUAL_DISTANCE_NOT_SPLICING) {
	      canonicalsgn = 9;
	      rev_score = prevlink->rev_score + CONSEC_POINTS_PER_MATCH;

	    } else if (near_end_p == false && prevlink->rev_consecutive < EXON_DEFN) {
	      canonicalsgn = 0;
	      rev_score = prevlink->rev_score - diffdist_penalty_splicing(diffdistance) - querydist_penalty(querydistance) - NINTRON_PENALTY_MISMATCH;

	    } else {
	      /* Called only near end */
#if 0
	      if (splicingp == false) {
		/* Not possible because splicingp == true && use_shifted_canonical_p == true */
		canonicalsgn = 0;
		rev_score = prevlink->rev_score - diffdist_penalty_nosplicing(diffdistance) - querydist_penalty(querydistance);
	      } else {
#endif

		leftpos = prevposition + querydistance + indexsize_nt - 1;
		/* printf("leftpos %d, last_leftpos %d, rightpos %d\n",leftpos,last_leftpos,rightpos); */
		if (leftpos == last_leftpos) {
		  canonicalp = last_canonicalp;
		} else {
		  debug7(printf("Calling find_shift_canonical rev\n"));
		  canonicalp = find_shifted_canonical(leftpos,rightpos,querydistance,
						      /* &lastCT,&lastAC, */
						      Genome_prev_antiacceptor_position,Genome_prev_antidonor_position,
						      chroffset,chrhigh,plusp,skip_repetitive_p);
		  last_leftpos = leftpos;
		  last_canonicalp = canonicalp;
		}
		if (canonicalp == true) {
		  canonicalsgn = -1;
		  rev_score = prevlink->rev_score - diffdist_penalty_splicing(diffdistance) - querydist_penalty(querydistance);
		} else {
		  canonicalsgn = 0;
		  rev_score = prevlink->rev_score - diffdist_penalty_splicing(diffdistance) - querydist_penalty(querydistance) - NINTRON_PENALTY_MISMATCH;
		}
#if 0
	      }
#endif

#if 0
	      if (prevlink->rev_consecutive < 20) {
		/* Intended to avoid too many jumps, but misses short introns */
		rev_score -= diffdist_penalty(diffdistance);
	      }
#endif
	    }

	  
	    debug9(printf("\tD-. Rev mismatch qpos %d,%d at %ux%d (score = %d -> %d, consec = %d (from %d), intr = %d-%d-%d, gendist %u, querydist %d, canonicalsgn %d)",
			  prev_querypos,prevhit,prevposition,active[prev_querypos][prevhit],
			  prevlink->rev_score,rev_score,prevlink->rev_consecutive,prevlink->rev_rootposition,
			  best_rev_intronnfwd,best_rev_intronnrev,best_rev_intronnunk,
			  gendistance,querydistance,canonicalsgn));
	    
#ifdef USE_SUBOPTIMAL_STARTS
	    rev_prevposes = Intlist_push(rev_prevposes,prev_querypos);
	    rev_prevhits = Intlist_push(rev_prevhits,prevhit);
	    rev_scores = Intlist_push(rev_scores,rev_score);
#endif

	    /* Allow ties, which should favor shorter intron */
	    if (rev_score >= best_rev_score) {
	      if (diffdistance <= EQUAL_DISTANCE_FOR_CONSECUTIVE) {
		best_rev_consecutive = prevlink->rev_consecutive + (querydistance + indexsize_nt);
		best_rev_rootposition = prevlink->rev_rootposition;
		/* best_rev_rootnlinks = prevlink->rev_rootnlinks + 1; */
	      } else {
		best_rev_consecutive = 0;
		best_rev_rootposition = position;
		/* best_rev_rootnlinks = 1; */
	      }
#ifdef USE_SUBOPTIMAL_STARTS
	      best_rev_initposition = prevlink->rev_initposition;
#endif
	      best_rev_score = rev_score;
	      best_rev_prevpos = prev_querypos;
	      best_rev_prevhit = prevhit;
	      best_rev_intronnfwd = prevlink->rev_intronnfwd;
	      best_rev_intronnrev = prevlink->rev_intronnrev;
	      best_rev_intronnunk = prevlink->rev_intronnunk;
	      switch (canonicalsgn) {
	      case 1: best_rev_intronnfwd++; break;
	      case -1: best_rev_intronnrev++; break;
	      case 0: best_rev_intronnunk++; break;
	      }
	      debug9(printf(" => Best rev at %d (consec = %d)\n",rev_score,best_rev_consecutive));
	    } else {
	      debug9(printf(" => Loses to %d\n",best_rev_score));
	    }
	  }

	  prevhit = active[prev_querypos][prevhit];
	}
      }
    }
#endif
  }

  /* Best_score needs to beat something positive to prevent a
     small local extension from beating a good canonical intron.
     If querypos is too small, don't insert an intron.  */
  /* linksconsecutive already assigned above */
  currlink->fwd_consecutive = best_fwd_consecutive;
  currlink->fwd_rootposition = best_fwd_rootposition;
  /* currlink->fwd_rootnlinks = best_fwd_rootnlinks; */
  currlink->fwd_pos = best_fwd_prevpos;
  currlink->fwd_hit = best_fwd_prevhit;
  if (localp && currlink->fwd_pos < 0) {
    currlink->fwd_score = indexsize_nt;
#ifdef USE_SUBOPTIMAL_STARTS
    currlink->fwd_initposition = position;
#endif
  } else {
    currlink->fwd_score = best_fwd_score;
#ifdef USE_SUBOPTIMAL_STARTS
    currlink->fwd_initposition = best_fwd_initposition;
#endif
  }
  currlink->fwd_intronnfwd = best_fwd_intronnfwd;
  currlink->fwd_intronnrev = best_fwd_intronnrev;
  currlink->fwd_intronnunk = best_fwd_intronnunk;

#ifndef PMAP
  if (splicingp == true) {
    /* linksconsecutive already assigned above */
    currlink->rev_consecutive = best_rev_consecutive;
    currlink->rev_rootposition = best_rev_rootposition;
    /* currlink->rev_rootnlinks = best_rev_rootnlinks; */
    currlink->rev_pos = best_rev_prevpos;
    currlink->rev_hit = best_rev_prevhit;
    if (localp && currlink->rev_pos < 0) {
      currlink->rev_score = indexsize_nt;
#ifdef USE_SUBOPTIMAL_STARTS
      currlink->rev_initposition = position;
#endif
    } else {
      currlink->rev_score = best_rev_score;
#ifdef USE_SUBOPTIMAL_STARTS
      currlink->rev_initposition = best_rev_initposition;
#endif
    }
    currlink->rev_intronnfwd = best_rev_intronnfwd;
    currlink->rev_intronnrev = best_rev_intronnrev;
    currlink->rev_intronnunk = best_rev_intronnunk;
  }
#endif

#ifdef PMAP
  debug9(printf("\tChose %d,%d with score %d (fwd)\n",
		currlink->fwd_pos,currlink->fwd_hit,currlink->fwd_score));
#else
  debug9(
	 if (splicingp == false) {
	   printf("\tChose %d,%d with score %d (fwd)\n",
		  currlink->fwd_pos,currlink->fwd_hit,currlink->fwd_score);
	 } else {
	   printf("\tChose %d,%d with score %d (fwd) and %d,%d with score %d (rev)\n",
		  currlink->fwd_pos,currlink->fwd_hit,currlink->fwd_score,currlink->rev_pos,currlink->rev_hit,currlink->rev_score);
	 });
#endif
  debug3(printf("%d %d  %d %d  1\n",querypos,hit,best_prevpos,best_prevhit));

#ifdef USE_SUBOPTIMAL_STARTS
  /* Add ties */
  for (p = fwd_prevposes, q = fwd_prevhits, r = fwd_scores; p != NULL; p = Intlist_next(p), q = Intlist_next(q), r = Intlist_next(r)) {
    prev_querypos = Intlist_head(p);
    prevhit = Intlist_head(q);
    fwd_score = Intlist_head(r);
    if (fwd_score < best_fwd_score - suboptimal_score_start) {
      /* Skip */
    } else if (prev_querypos == best_fwd_prevpos && prevhit == best_fwd_prevhit) {
      /* Skip */
    } else if (links[prev_querypos][prevhit].fwd_initposition == best_fwd_initposition) {
      /* Skip */
    } else {
      currlink->fwd_pos_ties = Intlist_push(currlink->fwd_pos_ties,prev_querypos);
      currlink->fwd_hit_ties = Intlist_push(currlink->fwd_hit_ties,prevhit);
    }
  }
  Intlist_free(&fwd_scores);
  Intlist_free(&fwd_prevhits);
  Intlist_free(&fwd_prevposes);

#ifndef PMAP
  for (p = rev_prevposes, q = rev_prevhits, r = rev_scores; p != NULL; p = Intlist_next(p), q = Intlist_next(q), r = Intlist_next(r)) {
    prev_querypos = Intlist_head(p);
    prevhit = Intlist_head(q);
    rev_score = Intlist_head(r);
    if (rev_score < best_rev_score - suboptimal_score_start) {
      /* Skip */
    } else if (prev_querypos == best_rev_prevpos && prevhit == best_rev_prevhit) {
      /* Skip */
    } else if (links[prev_querypos][prevhit].rev_initposition == best_rev_initposition) {
      /* Skip */
    } else {
      currlink->rev_pos_ties = Intlist_push(currlink->rev_pos_ties,prev_querypos);
      currlink->rev_hit_ties = Intlist_push(currlink->rev_hit_ties,prevhit);
    }
  }
  Intlist_free(&rev_scores);
  Intlist_free(&rev_prevhits);
  Intlist_free(&rev_prevposes);
#endif

#endif

  return;
}


/* For PMAP, indexsize is in aa. */
static void
score_querypos_splicing_no_shifted (Link_T currlink, int querypos, int hit,
				    int querystart, int queryend, unsigned int position,
				    struct Link_T **links, unsigned int **mappings,
				    int **active, int *firstactive,
#ifdef USE_SUBOPTIMAL_STARTS
				    int *fwd_initposition_bestscore, int *fwd_initposition_bestpos, int *fwd_initposition_besthit,
#endif
				    int indexsize, Intlist_T processed, int sufflookback, int nsufflookback, int maxintronlen, 
				    bool localp) {
  Link_T prevlink;
#ifdef USE_SUBOPTIMAL_STARTS
  Intlist_T fwd_prevposes = NULL, fwd_prevhits = NULL, fwd_scores = NULL, p, q, r;
  int initposition;
  int best_fwd_initposition = -1;
#else
  Intlist_T p;
#endif

  int best_fwd_consecutive = indexsize*NT_PER_MATCH;
  int best_fwd_rootposition = position;
  /* int best_fwd_rootnlinks = 1; */
  int best_fwd_score = 0, fwd_score;
  int best_fwd_prevpos = -1, best_fwd_prevhit = -1;
  int best_fwd_intronnfwd = 0, best_fwd_intronnrev = 0, best_fwd_intronnunk = 0;
  bool adjacentp = false, donep;
  int prev_querypos, prevhit;
  unsigned int prevposition, gendistance;
  int querydistance, diffdistance, lookback, nlookback, nseen, indexsize_nt;
  int canonicalsgn;
  int enough_consecutive;
  bool near_end_p;

#ifdef PMAP
  indexsize_nt = indexsize*3;
#else
  indexsize_nt = indexsize;
#endif

  enough_consecutive = 32;

  if (querypos < querystart + NEAR_END_LENGTH) {
    near_end_p = true;
  } else if (querypos > queryend - NEAR_END_LENGTH) {
    near_end_p = true;
  } else {
    near_end_p = false;
  }


  /* A. Evaluate adjacent position (at last one processed) */
  if (processed != NULL) {
    prev_querypos = Intlist_head(processed);
#ifdef PMAP
    querydistance = (querypos - prev_querypos)*3;
#else
    querydistance = querypos - prev_querypos;
#endif
    prevhit = firstactive[prev_querypos];
    while (prevhit != -1 && (prevposition = mappings[prev_querypos][prevhit]) + querydistance <= position) {
      if (prevposition + querydistance == position) {
	prevlink = &(links[prev_querypos][prevhit]);
	best_fwd_consecutive = prevlink->fwd_consecutive + querydistance;
	best_fwd_rootposition = prevlink->fwd_rootposition;
#ifdef USE_SUBOPTIMAL_STARTS
	best_fwd_initposition = prevlink->fwd_initposition;
#endif
	/* best_fwd_rootnlinks = prevlink->fwd_rootnlinks + 1; */
	best_fwd_score = prevlink->fwd_score + querydistance;

	best_fwd_prevpos = prev_querypos;
	best_fwd_prevhit = prevhit;
	best_fwd_intronnfwd = prevlink->fwd_intronnfwd;
	best_fwd_intronnrev = prevlink->fwd_intronnrev;
	best_fwd_intronnunk = prevlink->fwd_intronnunk;
	adjacentp = true;

	debug9(printf("\tA. Adjacent qpos %d,%d at %ux%d (scores = %d -> %d, consec = %d (from %d), intr = %d-%d-%d)\n",
		      prev_querypos,prevhit,prevposition,active[prev_querypos][prevhit],prevlink->fwd_score,
		      best_fwd_score,best_fwd_consecutive,best_fwd_rootposition,
		      best_fwd_intronnfwd,best_fwd_intronnrev,best_fwd_intronnunk));
	prevhit = -1;		/* Exit loop */
      } else {
	prevhit = active[prev_querypos][prevhit];
      }
    }

#ifdef USE_SUBOPTIMAL_STARTS
    /* Just push one for all adjacent positions */
    if (best_fwd_prevhit >= 0) {
      fwd_prevposes = Intlist_push(fwd_prevposes,best_fwd_prevpos);
      fwd_prevhits = Intlist_push(fwd_prevhits,best_fwd_prevhit);
      fwd_scores = Intlist_push(fwd_scores,best_fwd_score);
    }
#endif
  }

  if (best_fwd_consecutive < enough_consecutive) {

    /* D. Evaluate for mismatches (all other previous querypos) */
    /* Set parameters */
    if (adjacentp == true) {
      /* Look not so far back */
      nlookback = 1;
      lookback = sufflookback/2;
    } else {
      /* Look farther back */
      nlookback = nsufflookback;
      lookback = sufflookback;
    }

    donep = false;
    for (p = processed, nseen = 0;
	 p != NULL && best_fwd_consecutive < enough_consecutive && donep == false;
	 p = Intlist_next(p), nseen++) {
				  
      prev_querypos = Intlist_head(p);
#ifdef PMAP
      querydistance = (querypos - prev_querypos)*3 - indexsize_nt;
#else
      querydistance = querypos - prev_querypos - indexsize_nt;
#endif

      if (nseen > nlookback && querydistance > lookback) {
	donep = true;
      }

      prevhit = firstactive[prev_querypos];
      while (prevhit != -1 && (prevposition = mappings[prev_querypos][prevhit]) + indexsize_nt <= position) {
	/* printf("fwd: prevposition %u, prevhit %d\n",prevposition,prevhit); */
	prevlink = &(links[prev_querypos][prevhit]);

	gendistance = position - prevposition - indexsize_nt;
	diffdistance = abs(gendistance - querydistance);

	if (diffdistance < maxintronlen) {
	  if (diffdistance <= EQUAL_DISTANCE_NOT_SPLICING) {
	    canonicalsgn = 9;
	    fwd_score = prevlink->fwd_score + CONSEC_POINTS_PER_MATCH;
#ifdef PMAP
	    if (diffdistance % 3 != 0) {
	      fwd_score -= NONCODON_INDEL_PENALTY;
	    }
#endif
	  } else if (near_end_p == false && prevlink->fwd_consecutive < EXON_DEFN) {
	    canonicalsgn = 0;
	    if (splicingp == true) {
	      fwd_score = prevlink->fwd_score - diffdist_penalty_splicing(diffdistance) - querydist_penalty(querydistance) - NINTRON_PENALTY_MISMATCH;
	    } else {
	      fwd_score = prevlink->fwd_score - diffdist_penalty_nosplicing(diffdistance) - querydist_penalty(querydistance) - NINTRON_PENALTY_MISMATCH;
	    }

	  } else if (splicingp == false) {
	    canonicalsgn = 0;
	    fwd_score = prevlink->fwd_score - diffdist_penalty_nosplicing(diffdistance) - querydist_penalty(querydistance);
	  } else {
	    canonicalsgn = +1;
	    fwd_score = prevlink->fwd_score - diffdist_penalty_splicing(diffdistance) - querydist_penalty(querydistance);
	  }

	  debug9(printf("\tD. Fwd mismatch qpos %d,%d at %ux%d (score = %d -> %d, consec = %d (from %d), intr = %d-%d-%d, gendist %u, querydist %d, canonicalsgn %d)",
			prev_querypos,prevhit,prevposition,active[prev_querypos][prevhit],
			prevlink->fwd_score,fwd_score,prevlink->fwd_consecutive,prevlink->fwd_rootposition,
			best_fwd_intronnfwd,best_fwd_intronnrev,best_fwd_intronnunk,
			gendistance,querydistance,canonicalsgn));
	    
#ifdef USE_SUBOPTIMAL_STARTS
	  fwd_prevposes = Intlist_push(fwd_prevposes,prev_querypos);
	  fwd_prevhits = Intlist_push(fwd_prevhits,prevhit);
	  fwd_scores = Intlist_push(fwd_scores,fwd_score);
#endif

	  /* Allow ties, which should favor shorter intron */
	  if (fwd_score >= best_fwd_score) {
	    if (diffdistance <= EQUAL_DISTANCE_FOR_CONSECUTIVE) {
	      best_fwd_consecutive = prevlink->fwd_consecutive + (querydistance + indexsize_nt);
	      best_fwd_rootposition = prevlink->fwd_rootposition;
	      /* best_fwd_rootnlinks = prevlink->fwd_rootnlinks + 1; */
	    } else {
	      best_fwd_consecutive = 0;
	      best_fwd_rootposition = position;
	      /* best_fwd_rootnlinks = 1; */
	    }
#ifdef USE_SUBOPTIMAL_STARTS
	    best_fwd_initposition = prevlink->fwd_initposition;
#endif
	    best_fwd_score = fwd_score;
	    best_fwd_prevpos = prev_querypos;
	    best_fwd_prevhit = prevhit;
	    best_fwd_intronnfwd = prevlink->fwd_intronnfwd;
	    best_fwd_intronnrev = prevlink->fwd_intronnrev;
	    best_fwd_intronnunk = prevlink->fwd_intronnunk;
	    switch (canonicalsgn) {
	    case 1: best_fwd_intronnfwd++; break;
	    case 0: best_fwd_intronnunk++; break;
	    }
	    debug9(printf(" => Best fwd at %d (consec = %d)\n",fwd_score,best_fwd_consecutive));
	  } else {
	    debug9(printf(" => Loses to %d\n",best_fwd_score));
	  }
	}

	prevhit = active[prev_querypos][prevhit];
      }
    }
  }

  /* Best_score needs to beat something positive to prevent a
     small local extension from beating a good canonical intron.
     If querypos is too small, don't insert an intron.  */
  /* linksconsecutive already assigned above */
  currlink->fwd_consecutive = best_fwd_consecutive;
  currlink->fwd_rootposition = best_fwd_rootposition;
  /* currlink->fwd_rootnlinks = best_fwd_rootnlinks; */
  currlink->fwd_pos = best_fwd_prevpos;
  currlink->fwd_hit = best_fwd_prevhit;
  if (localp && currlink->fwd_pos < 0) {
    currlink->fwd_score = indexsize_nt;
#ifdef USE_SUBOPTIMAL_STARTS
    currlink->fwd_initposition = position;
#endif
  } else {
    currlink->fwd_score = best_fwd_score;
#ifdef USE_SUBOPTIMAL_STARTS
    currlink->fwd_initposition = best_fwd_initposition;
#endif
  }
  currlink->fwd_intronnfwd = best_fwd_intronnfwd;
  currlink->fwd_intronnrev = best_fwd_intronnrev;
  currlink->fwd_intronnunk = best_fwd_intronnunk;

  debug9(printf("\tChose %d,%d with score %d (fwd)\n",
	       currlink->fwd_pos,currlink->fwd_hit,currlink->fwd_score));
  debug3(printf("%d %d  %d %d  1\n",querypos,hit,best_prevpos,best_prevhit));

#ifdef USE_SUBOPTIMAL_STARTS
  /* Add ties */
  for (p = fwd_prevposes, q = fwd_prevhits, r = fwd_scores; p != NULL; p = Intlist_next(p), q = Intlist_next(q), r = Intlist_next(r)) {
    prev_querypos = Intlist_head(p);
    prevhit = Intlist_head(q);
    fwd_score = Intlist_head(r);
    if (fwd_score < best_fwd_score - suboptimal_score_start) {
      /* Skip */
    } else if (prev_querypos == best_fwd_prevpos && prevhit == best_fwd_prevhit) {
      /* Skip */
    } else if ((initposition = links[prev_querypos][prevhit].fwd_initposition) == best_fwd_initposition) {
      /* Skip */
    } else if (fwd_score <= fwd_initposition_bestscore[initposition]) {
      /* Skip */
    } else {
      currlink->fwd_pos_ties = Intlist_push(currlink->fwd_pos_ties,prev_querypos);
      currlink->fwd_hit_ties = Intlist_push(currlink->fwd_hit_ties,prevhit);
      fwd_initposition_bestscore[initposition] = fwd_score;
      fwd_initposition_bestpos[initposition] = querypos;
      fwd_initposition_besthit[initposition] = hit;
    }
  }
  Intlist_free(&fwd_scores);
  Intlist_free(&fwd_prevhits);
  Intlist_free(&fwd_prevposes);
#endif

  return;
}



static void
revise_active (int **active, int *firstactive, int *nactive, 
	       int low_hit, int high_hit, struct Link_T **links, int querypos) {
  int best_score, threshold, score;
  int hit, *ptr;

  debug6(printf("Revising querypos %d from low_hit %d to high_hit %d.  Scores:\n",querypos,low_hit,high_hit));
  if ((hit = low_hit) >= high_hit) {
    firstactive[querypos] = -1;
  } else {
    debug6(printf("At hit %d, fwd_score is %d",hit,links[querypos][hit].fwd_score));
    best_score = links[querypos][hit].fwd_score;
#ifndef PMAP
    debug6(printf(" and rev_score is %d",links[querypos][hit].rev_score));
    if ((score = links[querypos][hit].rev_score) > best_score) {
      best_score = score;
    }
#endif
    debug6(printf("\n"));

    for (hit++; hit < high_hit; hit++) {
      debug6(printf("At hit %d, fwd_score is %d",hit,links[querypos][hit].fwd_score));
      if ((score = links[querypos][hit].fwd_score) > best_score) {
	best_score = score;
      }
#ifndef PMAP
      debug6(printf(" and rev_score is %d",links[querypos][hit].rev_score));
      if ((score = links[querypos][hit].rev_score) > best_score) {
	best_score = score;
      }
#endif
      debug6(printf("\n"));
    }

    threshold = best_score - SCORE_FOR_RESTRICT;
    if (threshold < 0) {
      threshold = 0;
    }

    nactive[querypos] = 0;
    ptr = &(firstactive[querypos]);
    hit = low_hit;
    while (hit < high_hit) {
      while (hit < high_hit && links[querypos][hit].fwd_score <= threshold
#ifndef PMAP
	     && links[querypos][hit].rev_score <= threshold
#endif
	     ) {
	hit++;
      }
      *ptr = hit;
      if (hit < high_hit) {
	nactive[querypos] += 1;
	ptr = &(active[querypos][hit]);
	hit++;
      }
    }
    *ptr = -1;
  }

  debug6(
	 printf("Valid hits (%d) at querypos %d:",nactive[querypos],querypos);
	 hit = firstactive[querypos];
	 while (hit != -1) {
	   printf(" %d",hit);
	   hit = active[querypos][hit];
	 }
	 printf("\n");
	 );

  return;
}

static int **
intmatrix_1d_new (int length1, int *lengths2, int totallength) {
  int **matrix;
  int i;
  
  matrix = (int **) CALLOC(length1,sizeof(int *));
  matrix[0] = (int *) CALLOC(totallength,sizeof(int));
  for (i = 1; i < length1; i++) {
    if (lengths2[i-1] <= 0) {
      matrix[i] = matrix[i-1];
    } else {
      matrix[i] = &(matrix[i-1][lengths2[i-1]]);
    }
  }
  return matrix;
}

static void
intmatrix_1d_free (int ***matrix) {
  FREE((*matrix)[0]);
  FREE(*matrix);
  return;
}


static int **
intmatrix_2d_new (int length1, int *lengths2) {
  int **matrix;
  int i;
  
  matrix = (int **) CALLOC(length1,sizeof(int *));
  for (i = 0; i < length1; i++) {
    if (lengths2[i] <= 0) {
      matrix[i] = (int *) NULL;
    } else {
      matrix[i] = (int *) CALLOC(lengths2[i],sizeof(int));
    }
  }
  return matrix;
}

static void
intmatrix_2d_free (int ***matrix, int length1) {
  int i;

  for (i = 0; i < length1; i++) {
    if ((*matrix)[i]) {
      FREE((*matrix)[i]);
    }
  }
  FREE(*matrix);
  return;
}


/************************************************************************
 *   Cells used for ranking hits
 ************************************************************************/

typedef struct Cell_T *Cell_T;
struct Cell_T {
  int rootposition;
  int querypos;
  int hit;
  bool fwdp;
  int score;
};


static void
Cell_free (Cell_T *old) {
  FREE(*old);
  return;
}


static Cell_T
Cell_new (int rootposition, int querypos, int hit, bool fwdp, int score) {
  Cell_T new = (Cell_T) MALLOC(sizeof(*new));

  new->rootposition = rootposition;
  new->querypos = querypos;
  new->hit = hit;
  new->fwdp = fwdp;
  new->score = score;
  return new;
}


static int
Cell_rootposition_left_cmp (const void *a, const void *b) {
  Cell_T x = * (Cell_T *) a;
  Cell_T y = * (Cell_T *) b;

  if (x->rootposition < y->rootposition) {
    return -1;
  } else if (y->rootposition < x->rootposition) {
    return +1;
  } else if (x->score > y->score) {
    return -1;
  } else if (y->score > x->score) {
    return +1;
  } else if (x->querypos > y->querypos) {
    return -1;
  } else if (y->querypos > x->querypos) {
    return +1;
  } else if (x->hit < y->hit) {
    return -1;
  } else if (y->hit < x->hit) {
    return +1;
  } else if (x->fwdp == true && y->fwdp == false) {
    return -1;
  } else if (y->fwdp == true && x->fwdp == false) {
    return +1;
  } else {
    return 0;
  }
}

static int
Cell_rootposition_right_cmp (const void *a, const void *b) {
  Cell_T x = * (Cell_T *) a;
  Cell_T y = * (Cell_T *) b;

  if (x->rootposition < y->rootposition) {
    return -1;
  } else if (y->rootposition < x->rootposition) {
    return +1;
  } else if (x->score > y->score) {
    return -1;
  } else if (y->score > x->score) {
    return +1;
  } else if (x->querypos > y->querypos) {
    return -1;
  } else if (y->querypos > x->querypos) {
    return +1;
  } else if (x->hit > y->hit) {
    return -1;
  } else if (y->hit > x->hit) {
    return +1;
  } else if (x->fwdp == true && y->fwdp == false) {
    return -1;
  } else if (y->fwdp == true && x->fwdp == false) {
    return +1;
  } else {
    return 0;
  }
}


static int
Cell_score_cmp (const void *a, const void *b) {
  Cell_T x = * (Cell_T *) a;
  Cell_T y = * (Cell_T *) b;

  if (x->score > y->score) {
    return -1;
  } else if (y->score > x->score) {
    return +1;
  } else {
    return 0;
  }
}



static Cell_T *
Linkmatrix_get_cells_fwd (int *nunique, struct Link_T **links, int length1, int *npositions,
			  int indexsize, int bestscore, bool favor_right_p) {
  Cell_T *sorted, *cells;
  List_T celllist = NULL;
  int querypos, hit, lastpos;
  int rootposition, last_rootposition;
  int threshold_score;
  int ngood, ncells, i, k;

  lastpos = length1 - indexsize;

  if (bestscore > 2*suboptimal_score_end) {
    threshold_score = bestscore - suboptimal_score_end;
  } else {
    threshold_score = bestscore/2;
  }
  if (threshold_score <= indexsize) {
    threshold_score = indexsize + 1;
  }

  ncells = 0;
  for (querypos = 0; querypos <= lastpos; querypos++) {
    ngood = 0;
    for (hit = 0; hit < npositions[querypos]; hit++) {
      if (links[querypos][hit].fwd_score >= threshold_score) {
	ngood++;
      }
    }
    if (ngood > 0 && ngood <= 10) {
      for (hit = 0; hit < npositions[querypos]; hit++) {
	if (links[querypos][hit].fwd_score >= threshold_score) {
	  rootposition = links[querypos][hit].fwd_rootposition;
	  celllist = List_push(celllist,(void *) Cell_new(rootposition,querypos,hit,/*fwdp*/true,
							  links[querypos][hit].fwd_score));
	  ncells++;
	}
      }
    }
  }

  if (ncells == 0) {
    *nunique = 0;
    return (Cell_T *) NULL;

  } else {
    /* Take best result for each rootposition */
    cells = (Cell_T *) List_to_array(celllist,NULL);
    List_free(&celllist);

    if (favor_right_p == true) {
      qsort(cells,ncells,sizeof(Cell_T),Cell_rootposition_right_cmp);
    } else {
      qsort(cells,ncells,sizeof(Cell_T),Cell_rootposition_left_cmp);
    }

    sorted = (Cell_T *) CALLOC(ncells,sizeof(Cell_T));
    k = 0;

    last_rootposition = -1;
    for (i = 0; i < ncells; i++) {
      if (cells[i]->rootposition == last_rootposition) {
	Cell_free(&(cells[i]));
      } else {
	debug11(printf("Pushing position %d, score %d, pos %d, hit %d\n",
		       cells[i]->rootposition,cells[i]->score,cells[i]->querypos,cells[i]->hit));
	sorted[k++] = cells[i];
	last_rootposition = cells[i]->rootposition;
      }
    }
    debug11(printf("\n"));
    FREE(cells);
  
    *nunique = k;
    qsort(sorted,*nunique,sizeof(Cell_T),Cell_score_cmp);

    return sorted;
  }
}


static Cell_T *
Linkmatrix_get_cells_both (int *nunique, struct Link_T **links, int length1, int *npositions,
			   int indexsize, int bestscore, bool favor_right_p) {
  Cell_T *sorted, *cells;
  List_T celllist = NULL;
  int querypos, hit, lastpos;
  int rootposition, last_rootposition;
  int threshold_score;
  int ngood, ncells, i, k;

  lastpos = length1 - indexsize;

  if (bestscore > 2*suboptimal_score_end) {
    threshold_score = bestscore - suboptimal_score_end;
  } else {
    threshold_score = bestscore/2;
  }
  if (threshold_score <= indexsize) {
    threshold_score = indexsize + 1;
  }

  ncells = 0;
  for (querypos = 0; querypos <= lastpos; querypos++) {
    ngood = 0;
    for (hit = 0; hit < npositions[querypos]; hit++) {
      if (links[querypos][hit].fwd_score >= threshold_score) {
	ngood++;
      }
#ifndef PMAP
      if (links[querypos][hit].rev_score >= threshold_score) {
	ngood++;
      }
#endif
    }
    if (ngood > 0 && ngood <= 10) {
      for (hit = 0; hit < npositions[querypos]; hit++) {
	if (links[querypos][hit].fwd_score >= threshold_score) {
	  rootposition = links[querypos][hit].fwd_rootposition;
	  celllist = List_push(celllist,(void *) Cell_new(rootposition,querypos,hit,/*fwdp*/true,
							  links[querypos][hit].fwd_score));
	  ncells++;
	}
#ifndef PMAP
	if (links[querypos][hit].rev_score >= threshold_score) {
	  rootposition = links[querypos][hit].rev_rootposition;
	  celllist = List_push(celllist,(void *) Cell_new(rootposition,querypos,hit,/*fwdp*/false,
							  links[querypos][hit].rev_score));
	  ncells++;
	}
#endif
      }
    }
  }

  if (ncells == 0) {
    *nunique = 0;
    return (Cell_T *) NULL;

  } else {
    /* Take best result for each rootposition */
    cells = (Cell_T *) List_to_array(celllist,NULL);
    List_free(&celllist);

    if (favor_right_p == true) {
      qsort(cells,ncells,sizeof(Cell_T),Cell_rootposition_right_cmp);
    } else {
      qsort(cells,ncells,sizeof(Cell_T),Cell_rootposition_left_cmp);
    }

    sorted = (Cell_T *) CALLOC(ncells,sizeof(Cell_T));
    k = 0;

    last_rootposition = -1;
    for (i = 0; i < ncells; i++) {
      if (cells[i]->rootposition == last_rootposition) {
	Cell_free(&(cells[i]));
      } else {
	debug11(printf("position %d, score %d, pos %d, hit %d\n",
		       cells[i]->rootposition,cells[i]->score,cells[i]->querypos,cells[i]->hit));
	sorted[k++] = cells[i];
	last_rootposition = cells[i]->rootposition;
      }
    }
    debug11(printf("\n"));
    FREE(cells);

    *nunique = k;
    qsort(sorted,*nunique,sizeof(Cell_T),Cell_score_cmp);

    return sorted;
  }
}




/* Returns celllist */
/* For PMAP, indexsize is in aa. */
static Cell_T *
align_compute_scores (int *ncells, struct Link_T **links, unsigned int **mappings, int *npositions, int totalpositions,
		      bool oned_matrix_p, unsigned int *minactive, unsigned int *maxactive,
#ifdef USE_SUBOPTIMAL_STARTS
		      int *fwd_initposition_bestpos, int *fwd_initposition_besthit,
		      int *rev_initposition_bestpos, int *rev_initposition_besthit,
#endif
		      int querystart, int queryend, int querylength,

		      Genomicpos_T chroffset, Genomicpos_T chrhigh, bool plusp,

		      int indexsize, int sufflookback, int nsufflookback, int maxintronlen,
#ifdef DEBUG9
		      char *queryseq_ptr,
#endif
		      bool localp, bool skip_repetitive_p, 
		      bool use_shifted_canonical_p, bool debug_graphic_p, bool favor_right_p) {
  Cell_T *cells;
  Link_T currlink, prevlink;
  int querypos, indexsize_nt, hit, low_hit, high_hit;
  int nskipped, min_hits, specific_querypos, specific_low_hit, specific_high_hit, next_querypos;
  Intlist_T processed = NULL;
  int best_overall_score = 0;
  int grand_fwd_score, grand_fwd_querypos, grand_fwd_hit, best_fwd_hit, best_fwd_score;
#ifndef PMAP
  int grand_rev_score, grand_rev_querypos, grand_rev_hit, best_rev_hit, best_rev_score;
#endif
  int **active, *firstactive, *nactive;
#ifdef USE_SUBOPTIMAL_STARTS
  int *fwd_initposition_bestscore, *rev_initposition_bestscore;
#endif
  unsigned int position, prevposition;
#if 0
  int *lastGT, *lastAG;
#ifndef PMAP
  int *lastCT, *lastAC;
#endif
#endif
#ifdef DEBUG9
  char *oligo;
#endif
#ifdef DEBUG10
  Link_T termlink = NULL;
#endif

#ifdef PMAP
  indexsize_nt = indexsize*3;
#else
  indexsize_nt = indexsize;
#endif

#ifdef DEBUG9
  oligo = (char *) CALLOC(indexsize+1,sizeof(char));
#endif
  debug0(printf("Querystart = %d, queryend = %d, indexsize = %d\n",querystart,queryend,indexsize));

  if (oned_matrix_p == true) {
    active = intmatrix_1d_new(querylength,npositions,totalpositions);
  } else {
    active = intmatrix_2d_new(querylength,npositions);
  }

  firstactive = (int *) CALLOC(querylength,sizeof(int));
  nactive = (int *) CALLOC(querylength,sizeof(int));
#ifdef USE_SUBOPTIMAL_STARTS
  fwd_initposition_bestscore = (int *) CALLOC(genomiclength,sizeof(int));
  rev_initposition_bestscore = (int *) CALLOC(genomiclength,sizeof(int));
#endif

  /* Initialize */
  for (querypos = 0; querypos < querystart; querypos++) {
    firstactive[querypos] = -1;
    nactive[querypos] = 0;
  }
  while (querypos <= queryend - indexsize && npositions[querypos] <= 0) {
    firstactive[querypos] = -1;
    querypos++;
  }
  if (querypos <= queryend - indexsize) {
    for (hit = 0; hit < npositions[querypos]; hit++) {
      currlink = &(links[querypos][hit]);
#ifdef PMAP    
      currlink->fwd_pos = currlink->fwd_hit = -1;
      currlink->fwd_score = indexsize_nt;
      currlink->fwd_consecutive = indexsize_nt;
      currlink->fwd_rootposition = mappings[querypos][hit];
#ifdef USE_SUBOPTIMAL_STARTS
      currlink->fwd_initposition = mappings[querypos][hit];
      currlink->fwd_pos_ties = currlink->fwd_hit_ties = (Intlist_T) NULL;
#endif
      /* currlink->fwd_rootnlinks = 1; */
#else
      currlink->fwd_pos = currlink->fwd_hit = -1;
      currlink->fwd_score = indexsize_nt;
      currlink->fwd_consecutive = indexsize_nt;
      currlink->fwd_rootposition = mappings[querypos][hit];
#ifdef USE_SUBOPTIMAL_STARTS
      currlink->fwd_initposition = mappings[querypos][hit];
      currlink->fwd_pos_ties = currlink->fwd_hit_ties = (Intlist_T) NULL;
#endif
      /* currlink->fwd_rootnlinks = 1; */
      if (splicingp == true) {
	currlink->rev_pos = currlink->rev_hit = -1;
	currlink->rev_score = indexsize_nt;
	currlink->rev_consecutive = indexsize_nt;
	currlink->rev_rootposition = mappings[querypos][hit];
#ifdef USE_SUBOPTIMAL_STARTS
	currlink->rev_initposition = mappings[querypos][hit];
	currlink->rev_pos_ties = currlink->rev_hit_ties = (Intlist_T) NULL;
#endif
	/* currlink->rev_rootnlinks = 1; */
      }
#endif
    }
    revise_active(active,firstactive,nactive,0,npositions[querypos],links,querypos);
  }


  grand_fwd_score = 0;
  grand_fwd_querypos = -1;
  grand_fwd_hit = -1;
#ifndef PMAP
  if (splicingp == true) {
    grand_rev_score = 0;
    grand_rev_querypos = -1;
    grand_rev_hit = -1;
  }
#endif

  nskipped = 0;
  min_hits = 1000000;
  specific_querypos = -1;

  /* querypos += 1; -- this causes querypos 0 to be ignored */
  while (querypos <= queryend - indexsize) {
    best_fwd_score = 0;
    best_fwd_hit = -1;
#ifndef PMAP
    best_rev_score = 0;
    best_rev_hit = -1;
#endif
    
    hit = 0;
    while (hit < npositions[querypos] && mappings[querypos][hit] < minactive[querypos]) {
      hit++;
    }
    low_hit = hit;
    while (hit < npositions[querypos] && mappings[querypos][hit] < maxactive[querypos]) {
      hit++;
    }
    high_hit = hit;
    debug9(printf("Querypos %d has hit %d..%d out of %d (minactive = %u, maxactive = %u)\n",
		 querypos,low_hit,high_hit-1,npositions[querypos],minactive[querypos],maxactive[querypos]));

    /* Can't use nactive yet, so use high_hit - low_hit */
    if (skip_repetitive_p && high_hit - low_hit >= MAX_NACTIVE && nskipped <= MAX_SKIPPED) { /* Previously turned off */
      debug6(printf("Too many active (%d - %d) at querypos %d.  Setting firstactive to be -1\n",high_hit,low_hit,querypos));
      firstactive[querypos] = -1;
      nskipped++;
      debug9(printf("  %d skipped because of %d hits\n",nskipped,high_hit - low_hit + 1));

      /* Store most specific querypos in section of skipped */
      if (high_hit - low_hit < min_hits) {
	min_hits = high_hit - low_hit;
	specific_querypos = querypos;
	specific_low_hit = low_hit;
	specific_high_hit = high_hit;
      }
      querypos++;

    } else {
      if (nskipped > MAX_SKIPPED) {
	debug9(printf("Too many skipped.  Going back to specific querypos %d\n",specific_querypos));
	next_querypos = querypos;
	querypos = specific_querypos;
	low_hit = specific_low_hit;
	high_hit = specific_high_hit;
      } else {
	next_querypos = querypos + 1;
      }

      if (use_shifted_canonical_p == true) {
	for (hit = low_hit; hit < high_hit; hit++) {
	  currlink = &(links[querypos][hit]);
	  position = mappings[querypos][hit];
	  
	  debug9(strncpy(oligo,&(queryseq_ptr[querypos]),indexsize));
	  debug9(printf("Finding link at querypos %d,%d at %ux%d (%s).  prev_querypos was %d\n",
			querypos,hit,position,active[querypos][hit],oligo,processed ? Intlist_head(processed) : -1));

	  score_querypos_general(currlink,querypos,querystart,queryend,position,
				 links,mappings,active,firstactive,
#ifdef USE_SUBOPTIMAL_STARTS
				 fwd_initposition_bestscore,fwd_initposition_bestpos,fwd_initposition_besthit,
				 rev_initposition_bestscore,rev_initposition_bestpos,rev_initposition_besthit,
#endif
#if 0
				 lastGT,lastAG,
#ifndef PMAP
				 lastCT,lastAC,
#endif
#endif
				 chroffset,chrhigh,plusp,
				 indexsize,processed,sufflookback,nsufflookback,maxintronlen,
				 localp,skip_repetitive_p,use_shifted_canonical_p);

	  if (currlink->fwd_score > best_fwd_score) {
	    best_fwd_score = currlink->fwd_score;
	    best_fwd_hit = hit;
	  }
#ifndef PMAP
	  if (currlink->rev_score > best_rev_score) {
	    best_rev_score = currlink->rev_score;
	    best_rev_hit = hit;
	  }
#endif
	}

      } else {
	for (hit = low_hit; hit < high_hit; hit++) {
	  currlink = &(links[querypos][hit]);
	  position = mappings[querypos][hit];
	  
	  debug9(strncpy(oligo,&(queryseq_ptr[querypos]),indexsize));
	  debug9(printf("Finding link at querypos %d,%d at %ux%d (%s).  prev_querypos was %d\n",
			querypos,hit,position,active[querypos][hit],oligo,processed ? Intlist_head(processed) : -1));

	  score_querypos_splicing_no_shifted(currlink,querypos,hit,querystart,queryend,position,
					     links,mappings,active,firstactive,
#ifdef USE_SUBOPTIMAL_STARTS
					     fwd_initposition_bestscore,fwd_initposition_bestpos,fwd_initposition_besthit,
#endif
					     indexsize,processed,sufflookback,nsufflookback,maxintronlen,
					     localp);

	  if (currlink->fwd_score > best_fwd_score) {
	    best_fwd_score = currlink->fwd_score;
	    best_fwd_hit = hit;
	  }
	}
      }

      if (best_fwd_score > best_overall_score) {
	best_overall_score = best_fwd_score;
      }

      nskipped = 0;
      min_hits = 1000000;
      specific_querypos = -1;
      
#ifdef PMAP
      debug9(printf("Overall result at querypos %d yields best_fwd_hit %d\n",
		   querypos,best_fwd_hit));
#else
      debug9(printf("Overall result at querypos %d yields best_fwd_hit %d and best_rev_hit %d\n",
		   querypos,best_fwd_hit,best_rev_hit));
#endif

      if (splicingp == true && best_fwd_hit >= 0 && links[querypos][best_fwd_hit].fwd_hit < 0 && 
	  grand_fwd_querypos >= 0 && querypos >= grand_fwd_querypos + indexsize_nt) {
	prevlink = &(links[grand_fwd_querypos][grand_fwd_hit]);
	if ((best_fwd_score = prevlink->fwd_score - (querypos - grand_fwd_querypos)) > 0) {
	  prevposition = mappings[grand_fwd_querypos][grand_fwd_hit];
	  for (hit = low_hit; hit < high_hit; hit++) {
	    currlink = &(links[querypos][hit]);
	    if ((position = mappings[querypos][hit]) >= prevposition + indexsize_nt) {
	      currlink->fwd_consecutive = indexsize_nt;
	      currlink->fwd_rootposition = position;
#ifdef USE_SUBOPTIMAL_STARTS
	      currlink->fwd_initposition = prevlink->fwd_initposition;
#endif
	      /* currlink->fwd_rootnlinks = 1; */
	      currlink->fwd_pos = grand_fwd_querypos;
	      currlink->fwd_hit = grand_fwd_hit;
	      currlink->fwd_score = best_fwd_score;
	      currlink->fwd_intronnfwd = prevlink->fwd_intronnfwd;
	      currlink->fwd_intronnrev = prevlink->fwd_intronnrev;
	      currlink->fwd_intronnunk = prevlink->fwd_intronnunk + 1;
	    }
	  }
	  debug10(printf("At querypos %d, setting all fwd hits to point back to grand_fwd %d,%d with a score of %d\n",
		       querypos,grand_fwd_querypos,grand_fwd_hit,prevlink->fwd_score));
	}
      }

      /* Use >= to favor longer path in case of ties */
      if (best_fwd_hit >= 0 && best_fwd_score >= grand_fwd_score && 
	  links[querypos][best_fwd_hit].fwd_consecutive > EXON_DEFN) {
	grand_fwd_score = best_fwd_score;
	grand_fwd_querypos = querypos;
	grand_fwd_hit = best_fwd_hit;
	debug10(termlink = &(links[querypos][best_fwd_hit]));
	debug10(printf("At querypos %d, revising grand fwd to be hit %d with score of %d (pointing back to %d,%d)\n",
		     querypos,best_fwd_hit,best_fwd_score,termlink->fwd_pos,termlink->fwd_hit));
      }

#ifndef PMAP
      if (best_rev_score > best_overall_score) {
	best_overall_score = best_rev_score;
      }

      if (splicingp == false || use_shifted_canonical_p == false) {
	/* rev scores should be the same as the fwd scores */
      } else {
	if (best_rev_hit >= 0 && links[querypos][best_rev_hit].rev_hit < 0 && 
	    grand_rev_querypos >= 0 && querypos >= grand_rev_querypos + indexsize_nt) {
	  prevlink = &(links[grand_rev_querypos][grand_rev_hit]);
	  if ((best_rev_score = prevlink->rev_score - (querypos - grand_rev_querypos)) > 0) {
	    prevposition = mappings[grand_rev_querypos][grand_rev_hit];
	    for (hit = low_hit; hit < high_hit; hit++) {
	      currlink = &(links[querypos][hit]);
	      if ((position = mappings[querypos][hit]) >= prevposition + indexsize_nt) {
		currlink->rev_consecutive = indexsize_nt;
		currlink->rev_rootposition = position;
#ifdef USE_SUBOPTIMAL_STARTS
		currlink->rev_initposition = prevlink->rev_initposition;
#endif
		/* currlink->rev_rootnlinks = 1; */
		currlink->rev_pos = grand_rev_querypos;
		currlink->rev_hit = grand_rev_hit;
		currlink->rev_score = best_rev_score;
		currlink->rev_intronnrev = prevlink->rev_intronnfwd;
		currlink->rev_intronnrev = prevlink->rev_intronnrev;
		currlink->rev_intronnunk = prevlink->rev_intronnunk + 1;
	      }
	    }
	    debug10(printf("At querypos %d, setting all rev hits to point back to grand_rev %d,%d with a score of %d\n",
			  querypos,grand_rev_querypos,grand_rev_hit,prevlink->rev_score));
	  }
	}

	/* Use >= to favor longer path in case of ties */
	if (best_rev_hit >= 0 && best_rev_score >= grand_rev_score &&
	    links[querypos][best_rev_hit].rev_consecutive > EXON_DEFN) {
	  grand_rev_score = best_rev_score;
	  grand_rev_querypos = querypos;
	  grand_rev_hit = best_rev_hit;
	}
      }
#endif

      revise_active(active,firstactive,nactive,low_hit,high_hit,links,querypos);
      debug10(printf("Pushing querypos %d onto processed\n",querypos));
      processed = Intlist_push(processed,querypos);
      querypos = next_querypos;
    }
  }

  Intlist_free(&processed);

#if 0
  if (use_shifted_canonical_p == true) {
#ifndef PMAP
    if (lastAC != NULL) {
      FREE(lastAC);
    }
    if (lastCT != NULL) {
      FREE(lastCT);
    }
#endif
    if (lastAG != NULL) {
      FREE(lastAG);
    }
    if (lastGT != NULL) {
      FREE(lastGT);
    }
  }
#endif

  /* These are the final active oligomers, after pruning by score */
  if (debug_graphic_p == true) {
    mappings_dump_R(mappings,npositions,querylength,active,firstactive,indexsize,"active.mers");
  }

#ifdef USE_SUBOPTIMAL_STARTS
  FREE(rev_initposition_bestscore);
  FREE(fwd_initposition_bestscore);
#endif
  FREE(nactive);
  FREE(firstactive);

  if (oned_matrix_p == true) {
    intmatrix_1d_free(&active);
  } else {
    intmatrix_2d_free(&active,querylength);
  }


  /* Grand winners */
  debug10(printf("Finding grand winners, using root position method\n"));
  if (splicingp == false || use_shifted_canonical_p == false) {
    cells = Linkmatrix_get_cells_fwd(&(*ncells),links,querylength,npositions,indexsize,best_overall_score,favor_right_p);
  } else {
    cells = Linkmatrix_get_cells_both(&(*ncells),links,querylength,npositions,indexsize,best_overall_score,favor_right_p);
  }

  debug9(FREE(oligo));

  return cells;
}


static char complCode[128] = COMPLEMENT_LC;

/* genomicstart == chroffset + chrpos */
/* arguments were genomicpos, genomicstart, genomiclength */

static char
get_genomic_nt (char *g_alt, Genomicpos_T chrpos, Genomicpos_T chroffset,
		Genomicpos_T chrhigh, bool watsonp) {
  char c2, c2_alt;

  if (watsonp) {
    return Genome_get_char_blocks(&(*g_alt),chroffset + chrpos);

  } else {
    c2 = Genome_get_char_blocks(&c2_alt,chrhigh - chrpos);
    *g_alt = complCode[(int) c2_alt];
    return complCode[(int) c2];
  }
}


static List_T
traceback_one (int querypos, int hit, struct Link_T **links, Genomicpos_T **mappings,
	       char *queryseq_ptr, char *queryuc_ptr, 
#ifdef PMAP
	       Genomicpos_T chroffset, Genomicpos_T chrhigh, bool watsonp,
#endif
#ifdef DEBUG0
	       int indexsize,
#endif	       
	       Pairpool_T pairpool, bool fwdp) {
  List_T path = NULL;
  Genomicpos_T position;
  int prev_querypos, prevhit;
  char c2;
#ifdef PMAP
  char c2_alt;
#endif

#ifdef DEBUG0
  char *oligo;
#endif


  while (querypos >= 0) {
    position = mappings[querypos][hit];

#ifdef PMAP
    /* Change querypos positions from protein to nucleotide */
    c2 = get_genomic_nt(&c2_alt,position+2,chroffset,chrhigh,watsonp);
    path = Pairpool_push(path,pairpool,querypos*3+2,position+2,/*cdna*/c2,MATCH_COMP,c2,c2_alt,
			 /*dynprogindex*/0);
    c2 = get_genomic_nt(&c2_alt,position+1,chroffset,chrhigh,watsonp);
    path = Pairpool_push(path,pairpool,querypos*3+1,position+1,/*cdna*/c2,MATCH_COMP,c2,c2_alt,
			 /*dynprogindex*/0);
    c2 = get_genomic_nt(&c2_alt,position,chroffset,chrhigh,watsonp);
    path = Pairpool_push(path,pairpool,querypos*3,position,/*cdna*/c2,MATCH_COMP,c2,c2_alt,
			 /*dynprogindex*/0);
#else
    /* genomic nucleotide same as queryseq */
    c2 = queryuc_ptr[querypos];
    path = Pairpool_push(path,pairpool,querypos,position,queryseq_ptr[querypos],MATCH_COMP,
			 c2,/*genomealt*/c2,/*dynprogindex*/0);
#endif


#ifdef DEBUG0
    debug0(oligo = (char *) CALLOC(indexsize+1,sizeof(char)));
    debug0(strncpy(oligo,&(queryseq_ptr[querypos]),indexsize));
    if (fwdp == true) {
      debug0(printf("Pushing %d,%d (%s) at %u, score = %d, consec = %d (from %d), intr = %d(+)/%d(-)/%d(?)\n",
		    querypos,hit,oligo,position,
		    links[querypos][hit].fwd_score,links[querypos][hit].fwd_consecutive,links[querypos][hit].fwd_rootposition,
		    links[querypos][hit].fwd_intronnfwd,links[querypos][hit].fwd_intronnrev,
		    links[querypos][hit].fwd_intronnunk));
#ifndef PMAP
    } else {
      debug0(printf("Pushing %d,%d (%s) at %u, score = %d, consec = %d (from %d), intr = %d(+)/%d(-)/%d(?)\n",
		    querypos,hit,oligo,position,
		    links[querypos][hit].rev_score,links[querypos][hit].rev_consecutive,links[querypos][hit].rev_rootposition,
		    links[querypos][hit].rev_intronnfwd,links[querypos][hit].rev_intronnrev,
		    links[querypos][hit].rev_intronnunk));
#endif
    }
#endif
    debug0(FREE(oligo));

    /* prevposition = position; */
    prev_querypos = querypos;
    prevhit = hit;
    if (fwdp == true) {
      querypos = links[prev_querypos][prevhit].fwd_pos;
      hit = links[prev_querypos][prevhit].fwd_hit;
#ifndef PMAP
    } else {
      querypos = links[prev_querypos][prevhit].rev_pos;
      hit = links[prev_querypos][prevhit].rev_hit;
#endif
    }
    debug3(printf("%d %d  %d %d  3\n",prev_querypos,prevhit,querypos,hit));
  }

  return path;
}


static List_T
traceback_one_snps (int querypos, int hit, struct Link_T **links, Genomicpos_T **mappings,
		    char *queryseq_ptr, char *queryuc_ptr, 

		    Genomicpos_T chroffset, Genomicpos_T chrhigh, bool watsonp,
#ifdef DEBUG0
		    int indexsize,
#endif
		    Pairpool_T pairpool, bool fwdp) {
  List_T path = NULL;
  Genomicpos_T position;
  int prev_querypos, prevhit;
  char c2, c2_alt;

#ifdef DEBUG0
  char *oligo;
#endif


  while (querypos >= 0) {
    position = mappings[querypos][hit];

#ifdef PMAP
    /* Change querypos positions from protein to nucleotide */
    c2 = get_genomic_nt(&c2_alt,position+2,chroffset,chrhigh,watsonp);
    path = Pairpool_push(path,pairpool,querypos*3+2,position+2,/*cdna*/c2,MATCH_COMP,c2,c2_alt,
			 /*dynprogindex*/0);
    c2 = get_genomic_nt(&c2_alt,position+1,chroffset,chrhigh,watsonp);
    path = Pairpool_push(path,pairpool,querypos*3+1,position+1,/*cdna*/c2,MATCH_COMP,c2,c2_alt,
			 /*dynprogindex*/0);
    c2 = get_genomic_nt(&c2_alt,position,chroffset,chrhigh,watsonp);
    path = Pairpool_push(path,pairpool,querypos*3,position,/*cdna*/c2,MATCH_COMP,c2,c2_alt,
			 /*dynprogindex*/0);
#else
    /* genomic nucleotide same as queryseq */
    c2 = get_genomic_nt(&c2_alt,position,chroffset,chrhigh,watsonp);
    path = Pairpool_push(path,pairpool,querypos,position,queryseq_ptr[querypos],MATCH_COMP,c2,c2_alt,
			 /*dynprogindex*/0);
#endif


#ifdef DEBUG0
    debug0(oligo = (char *) CALLOC(indexsize+1,sizeof(char)));
    debug0(strncpy(oligo,&(queryseq_ptr[querypos]),indexsize));
    if (fwdp == true) {
      debug0(printf("Pushing %d,%d (%s) at %u, score = %d, consec = %d (from %d), intr = %d(+)/%d(-)/%d(?)\n",
		    querypos,hit,oligo,position,
		    links[querypos][hit].fwd_score,links[querypos][hit].fwd_consecutive,links[querypos][hit].fwd_rootposition,
		    links[querypos][hit].fwd_intronnfwd,links[querypos][hit].fwd_intronnrev,
		    links[querypos][hit].fwd_intronnunk));
#ifndef PMAP
    } else {
      debug0(printf("Pushing %d,%d (%s) at %u, score = %d, consec = %d (from %d), intr = %d(+)/%d(-)/%d(?)\n",
		    querypos,hit,oligo,position,
		    links[querypos][hit].rev_score,links[querypos][hit].rev_consecutive,links[querypos][hit].rev_rootposition,
		    links[querypos][hit].rev_intronnfwd,links[querypos][hit].rev_intronnrev,
		    links[querypos][hit].rev_intronnunk));
#endif
    }
#endif
    debug0(FREE(oligo));

    /* prevposition = position; */
    prev_querypos = querypos;
    prevhit = hit;
    if (fwdp == true) {
      querypos = links[prev_querypos][prevhit].fwd_pos;
      hit = links[prev_querypos][prevhit].fwd_hit;
#ifndef PMAP
    } else {
      querypos = links[prev_querypos][prevhit].rev_pos;
      hit = links[prev_querypos][prevhit].rev_hit;
#endif
    }
    debug3(printf("%d %d  %d %d  3\n",prev_querypos,prevhit,querypos,hit));
  }

  return path;
}



#ifdef USE_SUBOPTIMAL_STARTS
static List_T
traceback_ties (List_T path, int querypos, int hit, struct Link_T **links, Genomicpos_T **mappings,
		int *fwd_initposition_bestpos, int *fwd_initposition_besthit,
		int *rev_initposition_bestpos, int *rev_initposition_besthit,
		char *queryseq_ptr, Genomicpos_T chroffset, Genomicpos_T chrhigh, bool watsonp,
		Pairpool_T pairpool, int indexsize, bool fwdp) {
  List_T newpaths = NULL;
  List_T copy;
  Genomicpos_T position, initposition;
  Link_T currlink;
  List_T p;
  Intlist_T q, r;
  char c2, c2_alt;

#ifdef DEBUG0
  char *oligo;
#endif


  if (querypos < 0) {
    return List_push(NULL,(void *) path);

  } else {
    position = mappings[querypos][hit];
    currlink = &(links[querypos][hit]);

#ifdef PMAP
    /* Change querypos positions from protein to nucleotide */
    c2 = get_genomic_nt(&c2_alt,position+2,chroffset,chrhigh,watsonp);
    path = Pairpool_push(path,pairpool,querypos*3+2,position+2,/*cdna*/c2,MATCH_COMP,c2,c2_alt,
			 /*dynprogindex*/0);
    c2 = get_genomic_nt(&c2_alt,position+1,chroffset,chrhigh,watsonp);
    path = Pairpool_push(path,pairpool,querypos*3+1,position+1,/*cdna*/c2,MATCH_COMP,c2,c2_alt,
			 /*dynprogindex*/0);
    c2 = get_genomic_nt(&c2_alt,position,chroffset,chrhigh,watsonp);
    path = Pairpool_push(path,pairpool,querypos*3,position,/*cdna*/c2,MATCH_COMP,c2,c2_alt,
			 /*dynprogindex*/0);
#else
    /* genomic nucleotide same as queryseq */
    path = Pairpool_push(path,pairpool,querypos,position,queryseq_ptr[querypos],MATCH_COMP,
			 /*genome*/queryuc_ptr[querypos],/*genomealt*/GENOMEALT_DEFERRED,
			 /*dynprogindex*/0);
#endif

#ifdef DEBUG0
    debug0(oligo = (char *) CALLOC(indexsize+1,sizeof(char)));
    debug0(strncpy(oligo,&(queryseq_ptr[querypos]),indexsize));
    if (fwdp == true) {
      debug0(printf("Pushing %d,%d (%s) at %u, score = %d, consec = %d (from %d), intr = %d(+)/%d(-)/%d(?)\n",
		    querypos,hit,oligo,position,
		    links[querypos][hit].fwd_score,links[querypos][hit].fwd_consecutive,links[querypos][hit].fwd_rootposition,
		    links[querypos][hit].fwd_intronnfwd,links[querypos][hit].fwd_intronnrev,
		    links[querypos][hit].fwd_intronnunk));
#ifndef PMAP
    } else {
      debug0(printf("Pushing %d,%d (%s) at %u, score = %d, consec = %d (from %d), intr = %d(+)/%d(-)/%d(?)\n",
		    querypos,hit,oligo,position,
		    links[querypos][hit].rev_score,links[querypos][hit].rev_consecutive,links[querypos][hit].rev_rootposition,
		    links[querypos][hit].rev_intronnfwd,links[querypos][hit].rev_intronnrev,
		    links[querypos][hit].rev_intronnunk));
#endif
    }
#endif
    debug0(FREE(oligo));

    if (fwdp == true) {
      newpaths = traceback_ties(path,/*querypos*/currlink->fwd_pos,/*hit*/currlink->fwd_hit,links,mappings,
				fwd_initposition_bestpos,fwd_initposition_besthit,
				rev_initposition_bestpos,rev_initposition_besthit,
				queryseq_ptr,chroffset,chrhigh,watsonp,
				pairpool,indexsize,fwdp);
      for (q = currlink->fwd_pos_ties, r = currlink->fwd_hit_ties; q != NULL; q = Intlist_next(q), r = Intlist_next(r)) {
	initposition = links[Intlist_head(q)][Intlist_head(r)].fwd_initposition;
	if (fwd_initposition_bestpos[initposition] == querypos && fwd_initposition_besthit[initposition] == hit) {
	  debug0(printf("Recursing from %d,%d to %d,%d\n",querypos,hit,Intlist_head(q),Intlist_head(r)));
	  copy = Pairpool_copy(path,pairpool);
	  newpaths = List_append(newpaths,traceback_ties(copy,/*querypos*/Intlist_head(q),/*hit*/Intlist_head(r),links,mappings,
							 fwd_initposition_bestpos,fwd_initposition_besthit,
							 rev_initposition_bestpos,rev_initposition_besthit,
							 queryseq_ptr,chroffset,chrhigh,watsonp,
							 pairpool,indexsize,fwdp));
	  debug0(printf("Returning from recursion\n"));
	}
      }

#ifndef PMAP
    } else {
      newpaths = traceback_ties(path,/*querypos*/currlink->rev_pos,/*hit*/currlink->rev_hit,links,mappings,
				fwd_initposition_bestpos,fwd_initposition_besthit,
				rev_initposition_bestpos,rev_initposition_besthit,
				queryseq_ptr,chroffset,chrhigh,watsonp,
				pairpool,indexsize,fwdp);
      for (q = currlink->rev_pos_ties, r = currlink->rev_hit_ties; q != NULL; q = Intlist_next(q), r = Intlist_next(r)) {
	initposition = links[Intlist_head(q)][Intlist_head(r)].rev_initposition;
	if (rev_initposition_bestpos[initposition] == querypos && rev_initposition_besthit[initposition] == hit) {
	  debug0(printf("Recursing from %d,%d to %d,%d\n",querypos,hit,Intlist_head(q),Intlist_head(r)));
	  copy = Pairpool_copy(path,pairpool);
	  newpaths = List_append(newpaths,traceback_ties(copy,/*querypos*/Intlist_head(q),/*hit*/Intlist_head(r),links,mappings,
							 fwd_initposition_bestpos,fwd_initposition_besthit,
							 rev_initposition_bestpos,rev_initposition_besthit,
							 queryseq_ptr,chroffset,chrhigh,watsonp,
							 pairpool,indexsize,fwdp));
	  debug0(printf("Returning from recursion\n"));
	}
      }
#endif
    }

    return newpaths;
  }
}
#endif


/* Performs dynamic programming.  For PMAP, indexsize is in aa. */
static List_T
align_compute (Genomicpos_T **mappings, int *npositions, int totalpositions,
	       bool oned_matrix_p, unsigned int *minactive, unsigned int *maxactive,
	       char *queryseq_ptr, char *queryuc_ptr, int querylength, int queryseq_trim_start, int queryseq_trim_end,

	       Genomicpos_T chroffset, Genomicpos_T chrhigh, bool plusp,
	       int indexsize, int sufflookback, int nsufflookback, int maxintronlen, Pairpool_T pairpool,
	       bool localp, bool skip_repetitive_p, bool use_shifted_canonical_p,
	       bool favor_right_p, int max_nalignments, bool debug_graphic_p) {
  List_T all_paths = NULL;
  struct Link_T **links;
#ifdef USE_SUBOPTIMAL_STARTS
  int *fwd_initposition_bestpos, *fwd_initposition_besthit, *rev_initposition_bestpos, *rev_initposition_besthit;
  unsigned int position;
  int prev_querypos, prevhit;
#endif

  Cell_T *cells, cell;
  int ncells, i;

  bool fwdp;
  int querypos, hit;
  int querystart, queryend;
  int bestscore;


  querystart = queryseq_trim_start;
  queryend = queryseq_trim_end;

#ifdef USE_SUBOPTIMAL_STARTS
  fwd_initposition_bestpos = (int *) CALLOC(genomiclength,sizeof(int));
  fwd_initposition_besthit = (int *) CALLOC(genomiclength,sizeof(int));
  rev_initposition_bestpos = (int *) CALLOC(genomiclength,sizeof(int));
  rev_initposition_besthit = (int *) CALLOC(genomiclength,sizeof(int));
#endif

  if (oned_matrix_p == true) {
    links = Linkmatrix_1d_new(querylength,npositions,totalpositions);
  } else {
    links = Linkmatrix_2d_new(querylength,npositions);
  }

  /* These are all oligomers */
  if (debug_graphic_p == true) {
    mappings_dump_R(mappings,npositions,querylength,/*active*/NULL,/*firstactive*/NULL,indexsize,"all.mers");
  }
  
  cells = align_compute_scores(&ncells,links,mappings,npositions,totalpositions,
			       oned_matrix_p,minactive,maxactive,
#ifdef USE_SUBOPTIMAL_STARTS
			       fwd_initposition_bestpos,fwd_initposition_besthit,
			       rev_initposition_bestpos,rev_initposition_besthit,
#endif
			       querystart,queryend,querylength,
			       
			       chroffset,chrhigh,plusp,

			       indexsize,sufflookback,nsufflookback,maxintronlen,
#ifdef DEBUG9
			       queryseq_ptr,
#endif
			       localp,skip_repetitive_p,use_shifted_canonical_p,debug_graphic_p,
			       favor_right_p);

#ifdef PMAP
  debug1(Linkmatrix_print_fwd(links,mappings,querylength,npositions,queryseq_ptr,indexsize));
#else
  debug1(Linkmatrix_print_both(links,mappings,querylength,npositions,queryseq_ptr,indexsize));
#endif

  if (ncells == 0) {
    all_paths = (List_T) NULL;

  } else {
    bestscore = cells[0]->score;

    debug11(printf("Looping on %d cells, allowing up to %d alignments, plus any with best score %d\n",
		   ncells,max_nalignments,bestscore));

    for (i = 0; i < ncells && (i < max_nalignments || cells[i]->score == bestscore); i++) {
      cell = cells[i];
      querypos = cell->querypos;
      hit = cell->hit;
      fwdp = cell->fwdp;
      debug11(printf("Starting subpath %d at %d with score %d, querypos %d, hit %d\n",
		     i,mappings[querypos][hit],cell->score,querypos,hit));


      if (debug_graphic_p == true) {
	best_path_dump_R(links,mappings,querypos,hit,fwdp,"best.path");
	printf("plot(all.mers,col=\"black\",pch=\".\",xlab=\"Query\",ylab=\"Genomic\")\n");
	printf("points(active.mers,col=\"red\",pch=\".\")\n");
	printf("points(best.path,col=\"green\",pch=\".\")\n");
	printf("lines(querypos,minactive,col=\"blue\")\n");
	printf("lines(querypos,maxactive,col=\"blue\")\n");
      }


#ifdef USE_SUBOPTIMAL_STARTS
      all_paths = List_append(all_paths,traceback_ties(/*path*/(List_T) NULL,querypos,hit,links,mappings,
						       fwd_initposition_bestpos,fwd_initposition_besthit,
						       rev_initposition_bestpos,rev_initposition_besthit,
						       queryseq_ptr,chroffset,chrhigh,/*watsonp*/plusp,
						       pairpool,indexsize,fwdp));
#else
      if (snps_p == true) {
	all_paths = List_push(all_paths,(void *) traceback_one_snps(querypos,hit,links,mappings,queryseq_ptr,queryuc_ptr,	
								    chroffset,chrhigh,/*watsonp*/plusp,
#ifdef DEBUG0
								    indexsize,
#endif
								    pairpool,fwdp));
      } else {
	all_paths = List_push(all_paths,(void *) traceback_one(querypos,hit,links,mappings,queryseq_ptr,queryuc_ptr,	
#ifdef PMAP
							       chroffset,chrhigh,/*watsonp*/plusp,
#endif
#ifdef DEBUG0
							       indexsize,
#endif
							       pairpool,fwdp));
      }
#endif

    }
    debug11(printf("\n"));

    for (i = 0; i < ncells; i++) {
      cell = cells[i];
      Cell_free(&cell);
    }
    FREE(cells);
  }


#ifdef USE_SUBOPTIMAL_STARTS
  Linkmatrix_gc(links,querylength,npositions);
#endif
  if (oned_matrix_p == true) {
    Linkmatrix_1d_free(&links);
  } else {
    Linkmatrix_2d_free(&links,querylength);
  }

#ifdef USE_SUBOPTIMAL_STARTS
  FREE(rev_initposition_besthit);
  FREE(rev_initposition_bestpos);
  FREE(fwd_initposition_besthit);
  FREE(fwd_initposition_bestpos);
#endif

#if 0
  for (p = all_paths; p != NULL; p = List_next(p)) {
    Pair_dump_list(List_head(p),/*zerobasedp*/true);
    printf("\n");
  }
#endif

  return all_paths;
}


/* queryseq_ptr is NULL for PMAP.  querypos here is in nt. */
static List_T
convert_to_nucleotides (List_T path,
#ifndef PMAP
			char *queryseq_ptr, char *queryuc_ptr, 
#endif
			Genomicpos_T chroffset, Genomicpos_T chrhigh, bool watsonp,
			int query_offset, Pairpool_T pairpool, int indexsize_nt) {
  List_T pairs = NULL, pairptr;
  Pair_T pair;
  int querypos, lastquerypos, queryjump, genomejump, fill, default_fill;
  Genomicpos_T genomepos, lastgenomepos;
  char c, c_alt;

  debug5(printf("Beginning convert_to_nucleotides with %d pairs\n",List_length(path)));

  pairptr = path;
  path = Pairpool_pop(path,&pair);
  querypos = pair->querypos;
  genomepos = pair->genomepos;

#ifdef PMAP
  default_fill = indexsize_nt - 3;
#else
  default_fill = indexsize_nt - 1;
#endif

  lastquerypos = querypos + default_fill;
  lastgenomepos = genomepos + default_fill;
  while (lastquerypos > querypos) {
    debug5(printf("lastquerypos %d, lastgenomepos %d\n",
		  lastquerypos,lastgenomepos));

#ifdef PMAP
    c = get_genomic_nt(&c_alt,lastgenomepos,chroffset,chrhigh,watsonp);
    pairs = Pairpool_push(pairs,pairpool,lastquerypos+query_offset,lastgenomepos,/*cdna*/c,MATCH_COMP,c,c_alt,
			  /*dynprogindex*/0);
    debug5(printf("Pushing %c | %c at %d,%d\n",c,c,lastquerypos,lastgenomepos));
#elif defined(EXTRACT_GENOMICSEG)
    if (queryuc_ptr[lastquerypos] == genomicuc_ptr[lastgenomepos]) {
      pairs = Pairpool_push(pairs,pairpool,lastquerypos+query_offset,lastgenomepos,
			    queryseq_ptr[lastquerypos],MISMATCH_COMP,
			    genomicseg_ptr[lastgenomepos],/*genomealt*/GENOMEALT_DEFERRED,
			    /*dynprogindex*/0);
      debug5(printf("Pushing %c | %c at %d,%d\n",queryseq_ptr[lastquerypos],queryuc_ptr[lastquerypos],
		    lastquerypos+query_offset,lastgenomepos));
    } else {
      abort();
      pairs = Pairpool_push(pairs,pairpool,lastquerypos+query_offset,lastgenomepos,
			    queryseq_ptr[lastquerypos],MISMATCH_COMP,
			    genomicseg_ptr[lastgenomepos],/*genomealt*/GENOMEALT_DEFERRED,
			    /*dynprogindex*/0);
      debug5(printf("Pushing %c   %c at %d,%d\n",queryseq_ptr[lastquerypos],genomicseg_ptr[lastgenomepos],
		    lastquerypos+query_offset,lastgenomepos));
    }
#else
    if (mode == STANDARD) {
      c = queryuc_ptr[lastquerypos];
      pairs = Pairpool_push(pairs,pairpool,lastquerypos+query_offset,lastgenomepos,
			    queryseq_ptr[lastquerypos],MATCH_COMP,c,/*genomealt*/c,
			    /*dynprogindex*/0);
      debug5(printf("Pushing %c | %c at %d,%d\n",queryseq_ptr[lastquerypos],queryuc_ptr[lastquerypos],
		    lastquerypos+query_offset,lastgenomepos));
    } else {
      c = get_genomic_nt(&c_alt,lastgenomepos,chroffset,chrhigh,watsonp);
      if (queryuc_ptr[lastquerypos] == c) {
	pairs = Pairpool_push(pairs,pairpool,lastquerypos+query_offset,lastgenomepos,
			      queryseq_ptr[lastquerypos],MATCH_COMP,c,c_alt,/*dynprogindex*/0);
	debug5(printf("Pushing %c | %c at %d,%d\n",queryseq_ptr[lastquerypos],c,
		      lastquerypos+query_offset,lastgenomepos));
      } else {
	pairs = Pairpool_push(pairs,pairpool,lastquerypos+query_offset,lastgenomepos,
			      queryseq_ptr[lastquerypos],AMBIGUOUS_COMP,c,c_alt,/*dynprogindex*/0);
	debug5(printf("Pushing %c : %c at %d,%d\n",queryseq_ptr[lastquerypos],c,
		      lastquerypos+query_offset,lastgenomepos));
      }
    }
#endif
    --lastquerypos;
    --lastgenomepos;
  }

  /* Take care of first pair */
  if (mode == STANDARD) {
    pair->querypos += query_offset; /* Revise coordinates */
    /*pair->genomepos += genomic_offset;*/ /* Revise coordinates */
#ifdef WASTE
    pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
    pairs = List_push_existing(pairs,pairptr);
#endif
  } else {
    c = get_genomic_nt(&c_alt,pair->genomepos,chroffset,chrhigh,watsonp);
    if (pair->cdna == c) {
#ifdef WASTE
      pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
      pairs = List_push_existing(pairs,pairptr);
#endif
    } else {
      pairs = Pairpool_push(pairs,pairpool,pair->querypos+query_offset,pair->genomepos,
			    pair->cdna,AMBIGUOUS_COMP,c,c_alt,/*dynprogindex*/0);
      debug5(printf("Pushing %c : %c at %d,%d (first pair)\n",pair->cdna,c,
		    pair->querypos+query_offset,pair->genomepos));
    }
  }

  lastquerypos = querypos;
  lastgenomepos = genomepos;

  while (path != NULL) {
    pairptr = path;
    path = Pairpool_pop(path,&pair);
    querypos = pair->querypos;
    genomepos = pair->genomepos;
    
    queryjump = lastquerypos - 1 - querypos;
    genomejump = lastgenomepos - 1 - genomepos;

    if (queryjump == 0 && genomejump == 0) {
      /* Do nothing */
    } else {
      debug5(printf("At querypos %d, saw queryjump of %d and genomejump of %d\n",querypos,queryjump,genomejump));

      if (querypos + default_fill >= lastquerypos || genomepos + default_fill >= lastgenomepos) {
	if (lastquerypos - querypos < (int) (lastgenomepos - genomepos)) {
#if 0
	  /* This can occur with wobble mask */
	  fprintf(stderr,"Partial fill from querypos %d to %d (genomepos goes from %u to %u)\n",
		  querypos,lastquerypos,genomepos,lastgenomepos);
	  abort();
#endif
	  fill = lastquerypos - querypos - 1;
	} else {
#if 0
	  /* This can occur with wobble mask */
	  fprintf(stderr,"Partial fill from genomepos %u to %u (querypos goes from %d to %d)\n",
		  genomepos,lastgenomepos,querypos,lastquerypos);
	  abort();
#endif
	  fill = lastgenomepos - genomepos - 1;
	}
      } else {
	fill = default_fill;
      }

      lastquerypos = querypos + fill;
      lastgenomepos = genomepos + fill;
      debug5(printf("  Fill from querypos %d down to %d\n",lastquerypos,querypos));
      while (lastquerypos > querypos) {
#ifdef PMAP
	c = get_genomic_nt(&c_alt,lastgenomepos,chroffset,chrhigh,watsonp);
	pairs = Pairpool_push(pairs,pairpool,lastquerypos+query_offset,lastgenomepos,/*cdna*/c,MATCH_COMP,c,c_alt,
			      /*dynprogindex*/0);
	debug5(printf("Pushing %c | %c at %d,%d\n",c,c,lastquerypos+query_offset,lastgenomepos));
#elif defined(EXTRACT_GENOMICSEG)
	if (queryuc_ptr[lastquerypos] == genomicuc_ptr[lastgenomepos]) {
	  pairs = Pairpool_push(pairs,pairpool,lastquerypos+query_offset,lastgenomepos,
				queryseq_ptr[lastquerypos],MATCH_COMP,
				queryuc_ptr[lastquerypos],/*genomealt*/GENOMEALT_DEFERRED,
				/*dynprogindex*/0);
	  debug5(printf("Pushing %c | %c at %d,%d\n",queryseq_ptr[lastquerypos],genomicseg_ptr[lastgenomepos],
			lastquerypos+query_offset,lastgenomepos));
	} else {
	  abort();
	  pairs = Pairpool_push(pairs,pairpool,lastquerypos+query_offset,lastgenomepos,
				queryseq_ptr[lastquerypos],MISMATCH_COMP,
				genomicseg_ptr[lastgenomepos],/*genomealt*/GENOMEALT_DEFERRED,
				/*dynprogindex*/0);
	  debug5(printf("Pushing %c   %c at %d,%d\n",queryseq_ptr[lastquerypos],genomicseg_ptr[lastgenomepos],
			lastquerypos+query_offset,lastgenomepos));
	}
#else
	if (mode == STANDARD) {
	  /* assert(queryuc_ptr[lastquerypos] == get_genomic_nt(&c_alt,lastgenomepos,genomicstart,genomiclength,watsonp)); */
	  c = queryuc_ptr[lastquerypos];
	  pairs = Pairpool_push(pairs,pairpool,lastquerypos+query_offset,lastgenomepos,
				queryseq_ptr[lastquerypos],MATCH_COMP,c,/*genomealt*/c,
				/*dynprogindex*/0);
	  debug5(printf("Pushing %c | %c at %d,%d\n",queryseq_ptr[lastquerypos],queryuc_ptr[lastquerypos],
			lastquerypos+query_offset,lastgenomepos));
	} else {
	  c = get_genomic_nt(&c_alt,lastgenomepos,chroffset,chrhigh,watsonp);
	  if (queryuc_ptr[lastquerypos] == c) {
	    pairs = Pairpool_push(pairs,pairpool,lastquerypos+query_offset,lastgenomepos,
				  queryseq_ptr[lastquerypos],MATCH_COMP,c,c_alt,/*dynprogindex*/0);
	    debug5(printf("Pushing %c | %c at %d,%d\n",queryseq_ptr[lastquerypos],c,
			  lastquerypos+query_offset,lastgenomepos));
	  } else {
	    pairs = Pairpool_push(pairs,pairpool,lastquerypos+query_offset,lastgenomepos,
				  queryseq_ptr[lastquerypos],AMBIGUOUS_COMP,c,c_alt,/*dynprogindex*/0);
	    debug5(printf("Pushing %c : %c at %d,%d\n",queryseq_ptr[lastquerypos],c,
			  lastquerypos+query_offset,lastgenomepos));
	  }
	}
#endif
	--lastquerypos;
	--lastgenomepos;
      }
    }

    /* Take care of observed match */
    if (mode == STANDARD) {
      pair->querypos += query_offset; /* Revise coordinates */
      /*pair->genomepos += genomic_offset;*/ /* Revise coordinates */
#ifdef WASTE
      pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
      pairs = List_push_existing(pairs,pairptr);
#endif
    } else {
      c = get_genomic_nt(&c_alt,pair->genomepos,chroffset,chrhigh,watsonp);
      if (pair->cdna == c) {
#ifdef WASTE
	pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
	pairs = List_push_existing(pairs,pairptr);
#endif
      } else {
	pairs = Pairpool_push(pairs,pairpool,pair->querypos+query_offset,pair->genomepos,
			      pair->cdna,AMBIGUOUS_COMP,c,c_alt,/*dynprogindex*/0);
	debug5(printf("Pushing %c : %c at %d,%d (observed)\n",pair->cdna,c,
		      pair->querypos+query_offset,pair->genomepos));
      }
    }

    lastquerypos = querypos;
    lastgenomepos = genomepos;
  }

  debug5(Pair_dump_list(pairs,true));
  return List_reverse(pairs);
}


/* queryseq_ptr is NULL for PMAP.  querypos here is in nt. */
static List_T
convert_to_nucleotides_snps (List_T path,
#ifndef PMAP
			     char *queryseq_ptr, char *queryuc_ptr, 
#endif
			     Genomicpos_T chroffset, Genomicpos_T chrhigh, bool watsonp,
			     int query_offset, Pairpool_T pairpool, int indexsize_nt) {
  List_T pairs = NULL, pairptr;
  Pair_T pair;
  int querypos, genomepos, lastquerypos, lastgenomepos, queryjump, genomejump, fill, default_fill;
  char c, c_alt;

  debug5(printf("Beginning convert_to_nucleotides with %d pairs\n",List_length(path)));

  pairptr = path;
  path = Pairpool_pop(path,&pair);
  querypos = pair->querypos;
  genomepos = pair->genomepos;

#ifdef PMAP
  default_fill = indexsize_nt - 3;
#else
  default_fill = indexsize_nt - 1;
#endif

  lastquerypos = querypos + default_fill;
  lastgenomepos = genomepos + default_fill;
  while (lastquerypos > querypos) {
    debug5(printf("lastquerypos %d, lastgenomepos %d\n",
		  lastquerypos,lastgenomepos));

#ifdef PMAP
    c = get_genomic_nt(&c_alt,lastgenomepos,chroffset,chrhigh,watsonp);
    pairs = Pairpool_push(pairs,pairpool,lastquerypos+query_offset,lastgenomepos,/*cdna*/c,MATCH_COMP,c,c_alt,
			  /*dynprogindex*/0);
    debug5(printf("Pushing %c | %c at %d,%d\n",c,c,lastquerypos,lastgenomepos));
#elif defined(EXTRACT_GENOMICSEG)
    if (queryuc_ptr[lastquerypos] == genomicuc_ptr[lastgenomepos]) {
      pairs = Pairpool_push(pairs,pairpool,lastquerypos+query_offset,lastgenomepos,
			    queryseq_ptr[lastquerypos],MISMATCH_COMP,
			    genomicseg_ptr[lastgenomepos],/*genomealt*/GENOMEALT_DEFERRED,
			    /*dynprogindex*/0);
      debug5(printf("Pushing %c | %c at %d,%d\n",queryseq_ptr[lastquerypos],queryuc_ptr[lastquerypos],
		    lastquerypos+query_offset,lastgenomepos));
    } else {
      abort();
      pairs = Pairpool_push(pairs,pairpool,lastquerypos+query_offset,lastgenomepos,
			    queryseq_ptr[lastquerypos],MISMATCH_COMP,
			    genomicseg_ptr[lastgenomepos],/*genomealt*/GENOMEALT_DEFERRED,
			    /*dynprogindex*/0);
      debug5(printf("Pushing %c   %c at %d,%d\n",queryseq_ptr[lastquerypos],genomicseg_ptr[lastgenomepos],
		    lastquerypos+query_offset,lastgenomepos));
    }
#else
    if (mode == STANDARD) {
      /* assert(queryuc_ptr[lastquerypos] == get_genomic_nt(&c_alt,lastgenomepos,chroffset,chrhigh,watsonp)); */
      c = get_genomic_nt(&c_alt,lastgenomepos,chroffset,chrhigh,watsonp);
      pairs = Pairpool_push(pairs,pairpool,lastquerypos+query_offset,lastgenomepos,
			    queryseq_ptr[lastquerypos],MATCH_COMP,c,c_alt,
			    /*dynprogindex*/0);
      debug5(printf("Pushing %c | %c at %d,%d\n",queryseq_ptr[lastquerypos],queryuc_ptr[lastquerypos],
		    lastquerypos+query_offset,lastgenomepos));
    } else {
      c = get_genomic_nt(&c_alt,lastgenomepos,chroffset,chrhigh,watsonp);
      if (queryuc_ptr[lastquerypos] == c) {
	pairs = Pairpool_push(pairs,pairpool,lastquerypos+query_offset,lastgenomepos,
			      queryseq_ptr[lastquerypos],MATCH_COMP,c,c_alt,/*dynprogindex*/0);
	debug5(printf("Pushing %c | %c at %d,%d\n",queryseq_ptr[lastquerypos],c,
		      lastquerypos+query_offset,lastgenomepos));
      } else {
	pairs = Pairpool_push(pairs,pairpool,lastquerypos+query_offset,lastgenomepos,
			      queryseq_ptr[lastquerypos],AMBIGUOUS_COMP,c,c_alt,/*dynprogindex*/0);
	debug5(printf("Pushing %c : %c at %d,%d\n",queryseq_ptr[lastquerypos],c,
		      lastquerypos+query_offset,lastgenomepos));
      }
    }
#endif
    --lastquerypos;
    --lastgenomepos;
  }

  /* Take care of first pair */
  if (mode == STANDARD) {
    pair->querypos += query_offset; /* Revise coordinates */
    /*pair->genomepos += genomic_offset;*/ /* Revise coordinates */
#ifdef WASTE
    pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
    pairs = List_push_existing(pairs,pairptr);
#endif
  } else {
    c = get_genomic_nt(&c_alt,pair->genomepos,chroffset,chrhigh,watsonp);
    if (pair->cdna == c) {
#ifdef WASTE
      pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
      pairs = List_push_existing(pairs,pairptr);
#endif
    } else {
      pairs = Pairpool_push(pairs,pairpool,pair->querypos+query_offset,pair->genomepos,
			    pair->cdna,AMBIGUOUS_COMP,c,c_alt,/*dynprogindex*/0);
      debug5(printf("Pushing %c : %c at %d,%d (first pair)\n",pair->cdna,c,
		    pair->querypos+query_offset,pair->genomepos));
    }
  }

  lastquerypos = querypos;
  lastgenomepos = genomepos;

  while (path != NULL) {
    pairptr = path;
    path = Pairpool_pop(path,&pair);
    querypos = pair->querypos;
    genomepos = pair->genomepos;
    
    queryjump = lastquerypos - 1 - querypos;
    genomejump = lastgenomepos - 1 - genomepos;

    if (queryjump == 0 && genomejump == 0) {
      /* Do nothing */
    } else {
      debug5(printf("At querypos %d, saw queryjump of %d and genomejump of %d\n",querypos,queryjump,genomejump));

      if (querypos + default_fill >= lastquerypos || genomepos + default_fill >= lastgenomepos) {
	if (lastquerypos - querypos < lastgenomepos - genomepos) {
#if 0
	  /* This can occur with wobble mask */
	  fprintf(stderr,"Partial fill from querypos %d to %d (genomepos goes from %u to %u)\n",
		  querypos,lastquerypos,genomepos,lastgenomepos);
	  abort();
#endif
	  fill = lastquerypos - querypos - 1;
	} else {
#if 0
	  /* This can occur with wobble mask */
	  fprintf(stderr,"Partial fill from genomepos %u to %u (querypos goes from %d to %d)\n",
		  genomepos,lastgenomepos,querypos,lastquerypos);
	  abort();
#endif
	  fill = lastgenomepos - genomepos - 1;
	}
      } else {
	fill = default_fill;
      }

      lastquerypos = querypos + fill;
      lastgenomepos = genomepos + fill;
      debug5(printf("  Fill from querypos %d down to %d\n",lastquerypos,querypos));
      while (lastquerypos > querypos) {
#ifdef PMAP
	c = get_genomic_nt(&c_alt,lastgenomepos,chroffset,chrhigh,watsonp);
	pairs = Pairpool_push(pairs,pairpool,lastquerypos+query_offset,lastgenomepos,/*cdna*/c,MATCH_COMP,c,c_alt,
			      /*dynprogindex*/0);
	debug5(printf("Pushing %c | %c at %d,%d\n",c,c,lastquerypos+query_offset,lastgenomepos));
#elif defined(EXTRACT_GENOMICSEG)
	if (queryuc_ptr[lastquerypos] == genomicuc_ptr[lastgenomepos]) {
	  pairs = Pairpool_push(pairs,pairpool,lastquerypos+query_offset,lastgenomepos,
				queryseq_ptr[lastquerypos],MATCH_COMP,
				queryuc_ptr[lastquerypos],/*genomealt*/GENOMEALT_DEFERRED,
				/*dynprogindex*/0);
	  debug5(printf("Pushing %c | %c at %d,%d\n",queryseq_ptr[lastquerypos],genomicseg_ptr[lastgenomepos],
			lastquerypos+query_offset,lastgenomepos));
	} else {
	  abort();
	  pairs = Pairpool_push(pairs,pairpool,lastquerypos+query_offset,lastgenomepos,
				queryseq_ptr[lastquerypos],MISMATCH_COMP,
				genomicseg_ptr[lastgenomepos],/*genomealt*/GENOMEALT_DEFERRED,
				/*dynprogindex*/0);
	  debug5(printf("Pushing %c   %c at %d,%d\n",queryseq_ptr[lastquerypos],genomicseg_ptr[lastgenomepos],
			lastquerypos+query_offset,lastgenomepos));
	}
#else
	if (mode == STANDARD) {
	  c = get_genomic_nt(&c_alt,lastgenomepos,chroffset,chrhigh,watsonp);
	  pairs = Pairpool_push(pairs,pairpool,lastquerypos+query_offset,lastgenomepos,
				queryseq_ptr[lastquerypos],MATCH_COMP,c,c_alt,
				/*dynprogindex*/0);
	  debug5(printf("Pushing %c | %c at %d,%d\n",queryseq_ptr[lastquerypos],queryuc_ptr[lastquerypos],
			lastquerypos+query_offset,lastgenomepos));
	} else {
	  c = get_genomic_nt(&c_alt,lastgenomepos,chroffset,chrhigh,watsonp);
	  if (queryuc_ptr[lastquerypos] == c) {
	    pairs = Pairpool_push(pairs,pairpool,lastquerypos+query_offset,lastgenomepos,
				  queryseq_ptr[lastquerypos],MATCH_COMP,c,c_alt,/*dynprogindex*/0);
	    debug5(printf("Pushing %c | %c at %d,%d\n",queryseq_ptr[lastquerypos],c,
			  lastquerypos+query_offset,lastgenomepos));
	  } else {
	    pairs = Pairpool_push(pairs,pairpool,lastquerypos+query_offset,lastgenomepos,
				  queryseq_ptr[lastquerypos],AMBIGUOUS_COMP,c,c_alt,/*dynprogindex*/0);
	    debug5(printf("Pushing %c : %c at %d,%d\n",queryseq_ptr[lastquerypos],c,
			  lastquerypos+query_offset,lastgenomepos));
	  }
	}
#endif
	--lastquerypos;
	--lastgenomepos;
      }
    }

    /* Take care of observed match */
    if (mode == STANDARD) {
      pair->querypos += query_offset; /* Revise coordinates */
      /*pair->genomepos += genomic_offset;*/ /* Revise coordinates */
#ifdef WASTE
      pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
      pairs = List_push_existing(pairs,pairptr);
#endif
    } else {
      c = get_genomic_nt(&c_alt,pair->genomepos,chroffset,chrhigh,watsonp);
      if (pair->cdna == c) {
#ifdef WASTE
	pairs = Pairpool_push_existing(pairs,pairpool,pair);
#else
	pairs = List_push_existing(pairs,pairptr);
#endif
      } else {
	pairs = Pairpool_push(pairs,pairpool,pair->querypos+query_offset,pair->genomepos,
			      pair->cdna,AMBIGUOUS_COMP,c,c_alt,/*dynprogindex*/0);
	debug5(printf("Pushing %c : %c at %d,%d (observed)\n",pair->cdna,c,
		      pair->querypos+query_offset,pair->genomepos));
      }
    }

    lastquerypos = querypos;
    lastgenomepos = genomepos;
  }

  debug5(Pair_dump_list(pairs,true));
  return List_reverse(pairs);
}



/* Returns ncovered */
int
Stage2_scan (int *stage2_source, char *queryuc_ptr, int querylength,
#ifdef PMAP
	     char *genomicuc_ptr,
#endif
	     Genomicpos_T chrstart, Genomicpos_T chrend,
	     Genomicpos_T chroffset, Genomicpos_T chrhigh, bool plusp,
	     int genestrand, Oligoindex_T *oligoindices, int noligoindices,
	     Diagpool_T diagpool, bool debug_graphic_p, bool diagnosticp) {
  int ncovered;
  int source;
  int indexsize;
  Oligoindex_T oligoindex;
  unsigned int **mappings;
  bool *coveredp, oned_matrix_p;
  int *npositions, totalpositions;
  double pct_coverage;
  int maxnconsecutive;
  /* double diag_runtime; */
  List_T diagonals;
#ifndef USE_DIAGPOOL
  List_p;
  Diag_T diag;
#endif
#ifdef EXTRACT_GENOMICSEG
  int *counts;
#endif

  if (debug_graphic_p == true) {
    /* printf("par(mfrow=c(1,2),cex=0.2)\n"); */
    printf("par(cex=0.3)\n");
    printf("layout(matrix(c(1,2),1,2),widths=c(0.5,0.5),heights=c(1))\n");
  }

  coveredp = (bool *) CALLOC(querylength,sizeof(bool));
  mappings = (unsigned int **) CALLOC(querylength,sizeof(unsigned int *));
  npositions = (int *) CALLOC(querylength,sizeof(int));
  totalpositions = 0;
  maxnconsecutive = 0;

  source = 0;
  pct_coverage = 0.0;
  Diagpool_reset(diagpool);
  diagonals = (List_T) NULL;
  while (source < noligoindices && pct_coverage < SUFF_PCTCOVERAGE_OLIGOINDEX) {
    oligoindex = oligoindices[source];
    indexsize = Oligoindex_indexsize(oligoindex); /* Different sources can have different indexsizes */
#ifdef PMAP
    Oligoindex_tally(oligoindex,genomicuc_ptr,/*genomiclength*/chrend-chrstart,queryuc_ptr,querylength,
		     /*sequencepos*/0);
#else

#ifdef EXTRACT_GENOMICSEG
    Oligoindex_tally(oligoindex,genomicuc_ptr,/*genomiclength*/chrend-chrstart,queryuc_ptr,querylength,
		     /*sequencepos*/0);
    counts = Oligoindex_counts_copy(oligoindex);
#endif

    if (plusp == true) {
      Oligoindex_hr_tally(oligoindex,/*mappingstart*/chroffset+chrstart,
			  /*mappingend*/chroffset+chrend,/*plusp*/true,
			  queryuc_ptr,querylength,/*chrpos*/chrstart,genestrand);
    } else {
      Oligoindex_hr_tally(oligoindex,/*mappingstart*/chroffset+chrstart,
			  /*mappingend*/chroffset+chrend+1,/*plusp*/false,
			  queryuc_ptr,querylength,/*chrpos*/(chrhigh-chroffset)-chrend,genestrand);
    }

#ifdef EXTRACT_GENOMICSEG
    assert(Oligoindex_counts_equal(oligoindex,counts));
    /* Oligoindex_counts_dump(oligoindex,counts); */
    FREE(counts);
#endif

#endif

    diagonals = Oligoindex_get_mappings(diagonals,coveredp,mappings,npositions,&totalpositions,
					&oned_matrix_p,&maxnconsecutive,oligoindex,queryuc_ptr,
					querylength,chrstart,chrend,chroffset,chrhigh,plusp,diagpool);
    pct_coverage = Diag_update_coverage(coveredp,&ncovered,diagonals,querylength);
    debug(printf("Stage2_scan: source = %d, ncovered = %d, pct_coverage = %f\n",source,ncovered,pct_coverage));

    source++;
  }
  *stage2_source = source;

#ifdef USE_DIAGPOOL
  /* No need to free diagonals */
#else
    for (p = diagonals; p != NULL; p = List_next(p)) {
      diag = (Diag_T) List_head(p);
      Diag_free(&diag);
    }
    List_free(&diagonals);
#endif

  FREE(npositions);
  FREE(coveredp);
  FREE(mappings);		/* Don't need to free contents of mappings */

  for (source = 0; source < noligoindices; source++) {
    oligoindex = oligoindices[source];
    Oligoindex_untally(oligoindex);
  }

  return ncovered;
}


List_T
Stage2_compute (int *stage2_source, int *stage2_indexsize,
		char *queryseq_ptr, char *queryuc_ptr, int querylength, int query_offset,	
#ifdef PMAP
		char *genomicuc_ptr,
#endif
		Genomicpos_T chrstart, Genomicpos_T chrend,
		Genomicpos_T chroffset, Genomicpos_T chrhigh, bool plusp, int genestrand,
		Oligoindex_T *oligoindices, int noligoindices, double proceed_pctcoverage,
		Pairpool_T pairpool, Diagpool_T diagpool, int sufflookback, int nsufflookback,
		int maxintronlen, bool localp, bool skip_repetitive_p, bool use_shifted_canonical_p,
		bool favor_right_p, int max_nalignments, bool debug_graphic_p, bool diagnosticp,
		Stopwatch_T stopwatch, bool diag_debug) {
  List_T all_pairs = NULL, all_paths, path, pairs, p;
  int indexsize, indexsize_nt;
  Oligoindex_T oligoindex;
  Genomicpos_T **mappings;
  bool *coveredp, oned_matrix_p;
  int source;
  int *npositions, totalpositions;
  unsigned int *minactive, *maxactive;
  int ncovered;
  double pct_coverage;
  int maxnconsecutive;
  /* double diag_runtime; */
  List_T diagonals;
#ifndef USE_DIAGPOOL
  List_T p;
  Diag_T diag;
#endif
#ifdef DEBUG
  int nunique;
#endif

#ifdef EXTRACT_GENOMICSEG
  int *counts;
#endif

  debug(printf("Entered Stage2_compute with chrstart %u and chrend %u\n",chrstart,chrend));

  Stopwatch_start(stopwatch);

  if (debug_graphic_p == true) {
    /* printf("par(mfrow=c(1,2),cex=0.2)\n"); */
    printf("par(cex=0.3)\n");
    printf("layout(matrix(c(1,2),1,2),widths=c(0.5,0.5),heights=c(1))\n");
  }

  coveredp = (bool *) CALLOC(querylength,sizeof(bool));
  mappings = (Genomicpos_T **) CALLOC(querylength,sizeof(Genomicpos_T *));
  npositions = (int *) CALLOC(querylength,sizeof(int));
  totalpositions = 0;
  maxnconsecutive = 0;

  source = 0;
  pct_coverage = 0.0;
#ifdef USE_DIAGPOOL
  Diagpool_reset(diagpool);
#endif
  diagonals = (List_T) NULL;
  while (source < noligoindices && pct_coverage < SUFF_PCTCOVERAGE_OLIGOINDEX) {
    oligoindex = oligoindices[source];
    indexsize = Oligoindex_indexsize(oligoindex); /* Different sources can have different indexsizes */

#ifdef PMAP
    Oligoindex_tally(oligoindex,genomicuc_ptr,/*genomiclength*/chrend-chrstart,queryuc_ptr,querylength,
		     /*sequencepos*/0);
#else

#if 0
    /* Previously used this for user_genomicseg, but now creating genome_blocks on the fly */
    Oligoindex_tally(oligoindex,genomicuc_ptr,/*genomiclength*/chrend-chrstart,queryuc_ptr,querylength,
		     /*sequencepos*/0);
#endif

#ifdef EXTRACT_GENOMICSEG
    /* printf("indexsize = %d\n",indexsize); */
    /* printf("Query:  %.*s\n",querylength,queryuc_ptr); */
    /* printf("Genome: %s\n",genomicuc_ptr); */
    Oligoindex_tally(oligoindex,genomicuc_ptr,/*genomiclength*/mappingend-mappingstart,
		     queryuc_ptr,querylength,sequencepos);
    counts = Oligoindex_counts_copy(oligoindex);

    /* printf("plusp %d\n",plusp); */
    /* printf("genomicstart %u, genomicend %u, genomiclength %d\n",genomicstart,genomicend,genomiclength); */
    /* printf("mappingstart %u, mappingend %u\n",mappingstart,mappingend); */
#endif

    if (plusp == true) {
      Oligoindex_hr_tally(oligoindex,/*mappingstart*/chroffset+chrstart,
			  /*mappingend*/chroffset+chrend,/*plusp*/true,
			  queryuc_ptr,querylength,/*chrpos*/chrstart,genestrand);
    } else {
      Oligoindex_hr_tally(oligoindex,/*mappingstart*/chroffset+chrstart,
			  /*mappingend*/chroffset+chrend+1,/*plusp*/false,
			  queryuc_ptr,querylength,/*chrpos*/(chrhigh-chroffset)-chrend,genestrand);
    }

#ifdef EXTRACT_GENOMICSEG
    assert(Oligoindex_counts_equal(oligoindex,counts));
    /* Oligoindex_counts_dump(oligoindex,counts); */

    FREE(counts);
#endif

#endif

    diagonals = Oligoindex_get_mappings(diagonals,coveredp,mappings,npositions,&totalpositions,
					&oned_matrix_p,&maxnconsecutive,oligoindex,queryuc_ptr,
					querylength,chrstart,chrend,chroffset,chrhigh,plusp,diagpool);
    pct_coverage = Diag_update_coverage(coveredp,&ncovered,diagonals,querylength);
    debug(printf("Stage2_compute: source = %d, ncovered = %d, pct_coverage = %f\n",source,ncovered,pct_coverage));

    source++;
  }
  *stage2_source = source;
  *stage2_indexsize = indexsize;

  /* diag_runtime = */ Stopwatch_stop(stopwatch);


  minactive = (unsigned int *) CALLOC(querylength,sizeof(unsigned int));
  maxactive = (unsigned int *) CALLOC(querylength,sizeof(unsigned int));


  Stopwatch_start(stopwatch);

  if (diag_debug == true) {
    /* Do nothing */
  } else if (totalpositions == 0) {
    debug(printf("Quitting because totalpositions is zero\n"));

  } else if (querylength > 150 && pct_coverage < proceed_pctcoverage && ncovered < SUFF_NCOVERED) {
    /* Filter only on long queries */
    debug(printf("Quitting because querylength %d > 150, and pct_coverage is only %f < %f, and ncovered is only %d < %d, maxnconsecutive = %d\n",
		 querylength,pct_coverage,proceed_pctcoverage,ncovered,SUFF_NCOVERED,maxnconsecutive));

  } else {
    debug(printf("Proceeding because maxnconsecutive is %d and pct_coverage is %f > %f or ncovered = %d > %d\n",
		 maxnconsecutive,pct_coverage,proceed_pctcoverage,ncovered,SUFF_NCOVERED));

    debug(printf("Performing diag on genomiclength %u\n",chrend-chrstart));
    Diag_compute_bounds(minactive,maxactive,diagonals,querylength,
			debug_graphic_p,chrstart,chrend,chroffset,chrhigh,plusp);
    
    debug(
	  nunique = Diag_compute_bounds(minactive,maxactive,diagonals,querylength,
					debug_graphic_p,chrstart,chrend,chroffset,chrhigh,plusp);
	  fprintf(stderr,"%d diagonals (%d not dominated), maxnconsecutive = %d\n",
		  List_length(diagonals),nunique,maxnconsecutive);
	  );

    if (debug_graphic_p == true) {
      active_bounds_dump_R(minactive,maxactive,querylength);
      printf("lines(querypos,minactive,col=\"blue\")\n");
      printf("lines(querypos,maxactive,col=\"blue\")\n");
    }

    all_paths = align_compute(mappings,npositions,totalpositions,
			      oned_matrix_p,minactive,maxactive,
			      queryseq_ptr,queryuc_ptr,querylength,
			      /*query_trim_start*/0,/*query_trim_end*/querylength,
			      chroffset,chrhigh,plusp,
			      indexsize,sufflookback,nsufflookback,maxintronlen,pairpool,
			      localp,skip_repetitive_p,use_shifted_canonical_p,
			      favor_right_p,max_nalignments,debug_graphic_p);

#ifdef PMAP
    indexsize_nt = 3*indexsize;
#else
    indexsize_nt = indexsize;
#endif

    for (p = all_paths; p != NULL; p = List_next(p)) {
      path = (List_T) List_head(p);
      if (path != NULL) {
	if (snps_p == true) {
	  pairs = convert_to_nucleotides_snps(List_reverse(path),
#ifndef PMAP
					      queryseq_ptr,queryuc_ptr,
#endif
					      chroffset,chrhigh,/*watsonp*/plusp,
					      query_offset,pairpool,indexsize_nt);
	} else {
	  pairs = convert_to_nucleotides(List_reverse(path),
#ifndef PMAP
					 queryseq_ptr,queryuc_ptr,
#endif
					 chroffset,chrhigh,/*watsonp*/plusp,
					 query_offset,pairpool,indexsize_nt);
	}
	/* Don't need to free path, because its memory belongs to pairpool */
	all_pairs = List_push(all_pairs,(void *) pairs);
      }
    }

    List_free(&all_paths);
  }



  FREE(maxactive);
  FREE(minactive);
  FREE(npositions);
  FREE(coveredp);
  FREE(mappings);		/* Don't need to free contents of mappings */

  for (source = 0; source < noligoindices; source++) {
    oligoindex = oligoindices[source];
    Oligoindex_untally(oligoindex);
  }

  Stopwatch_stop(stopwatch);

  if (diag_debug == true) {
    return diagonals;
  } else {

#ifdef USE_DIAGPOOL
  /* No need to free diagonals */
#else
    for (p = diagonals; p != NULL; p = List_next(p)) {
      diag = (Diag_T) List_head(p);
      Diag_free(&diag);
    }
    List_free(&diagonals);
#endif
  }

  return all_pairs;
}



List_T
Stage2_compute_one (int *stage2_source, int *stage2_indexsize,
		    char *queryseq_ptr, char *queryuc_ptr, int querylength, int query_offset,	
#ifdef PMAP
		    char *genomicuc_ptr,
#endif
		    Genomicpos_T chrstart, Genomicpos_T chrend,
		    Genomicpos_T chroffset, Genomicpos_T chrhigh, bool plusp, int genestrand,

		    Oligoindex_T *oligoindices, int noligoindices, double proceed_pctcoverage,
		    Pairpool_T pairpool, Diagpool_T diagpool, int sufflookback, int nsufflookback,
		    int maxintronlen, bool localp, bool skip_repetitive_p, bool use_shifted_canonical_p,
		    bool favor_right_p, bool debug_graphic_p, bool diagnosticp,
		    Stopwatch_T stopwatch, bool diag_debug) {

  List_T pairs, all_pairs;

  all_pairs = Stage2_compute(&(*stage2_source),&(*stage2_indexsize),
			     queryseq_ptr,queryuc_ptr,querylength,query_offset,
#ifdef PMAP
			     genomicuc_ptr,
#endif
			     chrstart,chrend,chroffset,chrhigh,plusp,
			     genestrand,oligoindices,noligoindices,proceed_pctcoverage,
			     pairpool,diagpool,sufflookback,nsufflookback,
			     maxintronlen,localp,skip_repetitive_p,use_shifted_canonical_p,
			     favor_right_p,/*max_nalignments*/1,debug_graphic_p,
			     diagnosticp,stopwatch,diag_debug);
  if (all_pairs == NULL) {
    return (List_T) NULL;
  } else {
    pairs = (List_T) List_head(all_pairs);
    List_free(&all_pairs);
    return pairs;
  }
}


