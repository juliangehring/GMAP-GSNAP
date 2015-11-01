static char rcsid[] = "$Id: sarray-read.c 132776 2014-04-09 01:03:33Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifndef HAVE_MEMCPY
#define memcpy(d,s,n) bcopy((s),(d),(n))
#endif

#include "sarray-read.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>		/* For munmap */
#include "mem.h"
#include "bool.h"
#include "assert.h"
#include "access.h"
#include "listdef.h"
#include "list.h"
#include "genome_hr.h"
#include "splice.h"
#include "indel.h"
#include "stage3hr.h"
#include "bitpack64-read.h"
#include "bitpack64-access.h"

#ifdef USE_CHILD_BP
#include "bp.h"
#include "bp-read.h"
#endif


#ifdef HAVE_SSE2
#include <emmintrin.h>
#endif


/* A value of 10000 misses various splices, although they are caught by GSNAP algorithm */
#define EXCESS_SARRAY_HITS 100000
#define GUESS_ALLOCATION 10

#define get_bit(bitvector,i) ((bitvector)[(i)/32] & (1 << ((i)%32)))

/* Results of each suffix array search */
#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Details of suffix array search */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

/* Search through saindex */
#ifdef DEBUG1A
#define debug1a(x) x
#else
#define debug1a(x)
#endif

/* get_child */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

/* known splicing */
#ifdef DEBUG4S
#define debug4s(x) x
#else
#define debug4s(x)
#endif

/* find_multimiss_iter */
#ifdef DEBUG7
#define debug7(x) x
#else
#define debug7(x)
#endif

/* find_multimiss_iter details */
#ifdef DEBUG7A
#define debug7a(x) x
#else
#define debug7a(x)
#endif

/* Comparing SIMD with non-SIMD */
#ifdef DEBUG8
#define debug8(x) x
#else
#define debug8(x)
#endif

/* binary_search */
#ifdef DEBUG10
#define debug10(x) x
#else
#define debug10(x)
#endif

/* Compare sarray_search with sarray_search_simple */
#ifdef DEBUG14
#define debug14(x) x
#else
#define debug14(x)
#endif


#define T Sarray_T
struct T {
  Univcoord_T n;
  Univcoord_T n_plus_one;

  Univcoord_T *array;

  UINT4 *plcp_ptrs, *plcp_comp;		/* permuted LCP */

#ifdef USE_CHILD_BP
  UINT4 *childbp;		    /* child balanced parentheses */
  UINT4 *childfc;		    /* first child (array B in CST++) */
  UINT4 *childs_pages, *childs_ptrs, *childs_comp; /* child select */
  UINT4 *childr_ptrs, *childr_comp; /* child rank */
  UINT4 *childx_ptrs, *childx_comp; /* child block excess */

  UINT4 *pioneerbp;		/* pioneer balanced parentheses */
  UINT4 *pior_ptrs, *pior_comp;	/* pioneer rank */
  UINT4 *piom_ptrs, *piom_comp;	/* pioneer matches */
#else
  UINT4 *child_ptrs;
  UINT4 *child_comp;
  UINT4 *nextp;
#endif

  Sarrayptr_T initindexi[4];	/* For A, C, G, T */
  Sarrayptr_T initindexj[4];	/* For A, C, G, T */

  int indexsize;
  UINT4 indexspace;		/* 4^indexsize.  Used by sarray_search to detect when we have a poly-T oligo shorter than indexsize */
  UINT4 *indexi_ptrs, *indexi_comp; /* oligomer lookup into suffix array */
  UINT4 *indexj_ptrs, *indexj_comp;

  Access_T access;

  int array_fd; size_t array_len;
  int plcp_ptrs_fd; size_t plcp_ptrs_len; int plcp_comp_fd; size_t plcp_comp_len;

#ifdef USE_CHILD_BP
  int childbp_fd; size_t childbp_len;
  int childfc_fd; size_t childfc_len;
  int childs_pages_fd; size_t childs_pages_len;
  int childs_ptrs_fd; size_t childs_ptrs_len; int childs_comp_fd; size_t childs_comp_len;
  int childr_ptrs_fd; size_t childr_ptrs_len; int childr_comp_fd; size_t childr_comp_len;
  int childx_ptrs_fd; size_t childx_ptrs_len; int childx_comp_fd; size_t childx_comp_len;

  int pioneerbp_fd; size_t pioneerbp_len;
  int pior_ptrs_fd; size_t pior_ptrs_len; int pior_comp_fd; size_t pior_comp_len;
  int piom_ptrs_fd; size_t piom_ptrs_len; int piom_comp_fd; size_t piom_comp_len;
#else
  int child_ptrs_fd; size_t child_ptrs_len; int child_comp_fd; size_t child_comp_len;
  int nextp_fd; size_t nextp_len;
#endif

  int indexi_ptrs_fd; size_t indexi_ptrs_len; int indexi_comp_fd; size_t indexi_comp_len;
  int indexj_ptrs_fd; size_t indexj_ptrs_len; int indexj_comp_fd; size_t indexj_comp_len;
};


/* For benchmarking */
UINT4 *
Sarray_plcp_ptrs (Sarray_T this) {
  return this->plcp_ptrs;
}

/* For benchmarking */
UINT4 *
Sarray_plcp_comp (Sarray_T this) {
  return this->plcp_comp;
}

/* For benchmarking */
Univcoord_T
Sarray_size (Sarray_T this) {
  return this->n_plus_one;
}


static Genome_T genome;
static Sarray_T sarray;
static Univ_IIT_T chromosome_iit;
static int circular_typeint;
static int splicing_penalty;

static Chrpos_T overall_max_distance;
static Chrpos_T shortsplicedist;
static Chrpos_T max_deletionlen;
static Chrpos_T max_insertionlen;
static Chrpos_T max_end_deletions;

/* Splicing */
static Univcoord_T *splicesites;
static Splicetype_T *splicetypes;
static Chrpos_T *splicedists;
static int nsplicesites;


/* Simplified from sarray_search_simple in sarray-write.c */
static void
sarray_search_char (Sarrayptr_T *initptr, Sarrayptr_T *finalptr, char desired_char,
		    UINT4 *SA, UINT4 n) {
  Sarrayptr_T low, high, mid;
  Univcoord_T pos;
  char c;


  low = 1;
  high = n + 1;

  while (low < high) {
    /* Compute mid for unsigned ints.  Want floor((low+high)/2). */
    mid = low/2 + high/2;
    if (low % 2 == 1 && high % 2 == 1) {
      mid += 1;
    }
    pos = SA[mid];
    c = Genome_get_char_lex(genome,pos,n);
    if (desired_char > c) {
      low = mid + 1;
    } else {
      high = mid;
    }
  }

  *initptr = low;

  low--;
  high = n;
  while (low < high) {
    /* Compute mid for unsigned ints.  Want ceil((low+high)/2). */
    mid = low/2 + high/2;
    if (low % 2 == 1 || high % 2 == 1) {
      mid += 1;
    }
    pos = SA[mid];
    c = Genome_get_char_lex(genome,pos,n);
    if (desired_char >= c) {
      low = mid;
    } else {
      high = mid - 1;
    }
  }

  *finalptr = high;
  return;
}


void
Sarray_setup (T sarray_in, Genome_T genome_in, Univ_IIT_T chromosome_iit_in, int circular_typeint_in,
	      Chrpos_T shortsplicedist_in, int splicing_penalty_in,
	      int max_deletionlength, int max_end_deletions_in,
	      int max_middle_insertions, int max_end_insertions,
	      Univcoord_T *splicesites_in, Splicetype_T *splicetypes_in,
	      Chrpos_T *splicedists_in, int nsplicesites_in) {

  sarray = sarray_in;
  genome = genome_in;
  chromosome_iit = chromosome_iit_in;
  circular_typeint = circular_typeint_in;
  shortsplicedist = shortsplicedist_in;
  splicing_penalty = splicing_penalty_in;

  max_deletionlen = max_deletionlength;
  max_end_deletions = max_end_deletions_in;
  if (max_middle_insertions > max_end_insertions) {
    max_insertionlen = max_middle_insertions;
  } else {
    max_insertionlen = max_end_insertions;
  }

  if (shortsplicedist > max_deletionlen) {
    overall_max_distance = shortsplicedist;
  } else {
    overall_max_distance = max_deletionlen;
  }

  splicesites = splicesites_in;
  splicetypes = splicetypes_in;
  splicedists = splicedists_in;
  nsplicesites = nsplicesites_in;

  Bitpack64_read_setup();

  sarray_search_char(&(sarray->initindexi[0]),&(sarray->initindexj[0]),/*desired_char*/'A',sarray->array,sarray->n);
  sarray_search_char(&(sarray->initindexi[1]),&(sarray->initindexj[1]),/*desired_char*/'C',sarray->array,sarray->n);
  sarray_search_char(&(sarray->initindexi[2]),&(sarray->initindexj[2]),/*desired_char*/'G',sarray->array,sarray->n);
  sarray_search_char(&(sarray->initindexi[3]),&(sarray->initindexj[3]),/*desired_char*/'T',sarray->array,sarray->n);

#if 0
  printf("A => %u %u\n",sarray->initindexi[0],sarray->initindexj[0]);
  printf("C => %u %u\n",sarray->initindexi[1],sarray->initindexj[1]);
  printf("G => %u %u\n",sarray->initindexi[2],sarray->initindexj[2]);
  printf("T => %u %u\n",sarray->initindexi[3],sarray->initindexj[3]);
#endif

  return;
}


static int
log4 (int result) {
  int exponent = 0;

  while (result > 1) {
    result /= 4;
    exponent++;
  }

  return exponent;
}

static UINT4
power (int base, int exponent) {
  UINT4 result = 1;
  int i;

  for (i = 0; i < exponent; i++) {
    result *= base;
  }

  return result;
}




/* Ignores snps_root */
T
Sarray_new (char *dir, char *fileroot, char *snps_root, Access_mode_T access) {
  T new;
  char *comma1;
  double seconds;
  int npages;

  char *sarrayfile;
  char *plcp_ptrsfile, *plcp_compfile;
#ifdef USE_CHILD_BP
  char *childbpfile, *childfcfile, *childs_pagesfile, *childs_ptrsfile, *childs_compfile, *childr_ptrsfile, *childr_compfile;
  char *childx_ptrsfile, *childx_compfile;
  char *pioneerbpfile, *pior_ptrsfile, *pior_compfile, *piom_ptrsfile, *piom_compfile;
#else
  char *child_ptrsfile, *child_compfile, *nextpfile;
#endif
  char *indexi_ptrsfile, *indexi_compfile, *indexj_ptrsfile, *indexj_compfile;


  sarrayfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(".sarray")+1,sizeof(char));
  sprintf(sarrayfile,"%s/%s.sarray",dir,fileroot);

  plcp_ptrsfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(".saplcpptrs")+1,sizeof(char));
  sprintf(plcp_ptrsfile,"%s/%s.saplcpptrs",dir,fileroot);
  plcp_compfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(".saplcpcomp")+1,sizeof(char));
  sprintf(plcp_compfile,"%s/%s.saplcpcomp",dir,fileroot);

#ifdef USE_CHILD_BP
  childbpfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(".sachildbp")+1,sizeof(char));
  sprintf(childbpfile,"%s/%s.sachildbp",dir,fileroot);
  childfcfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(".sachildfc")+1,sizeof(char));
  sprintf(childfcfile,"%s/%s.sachildfc",dir,fileroot);

  childs_pagesfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(".sachildspages")+1,sizeof(char));
  sprintf(childs_pagesfile,"%s/%s.sachildspages",dir,fileroot);
  childs_ptrsfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(".sachildsptrs")+1,sizeof(char));
  sprintf(childs_ptrsfile,"%s/%s.sachildsptrs",dir,fileroot);
  childs_compfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(".sachildscomp")+1,sizeof(char));
  sprintf(childs_compfile,"%s/%s.sachildscomp",dir,fileroot);

  childr_ptrsfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(".sachildrptrs")+1,sizeof(char));
  sprintf(childr_ptrsfile,"%s/%s.sachildrptrs",dir,fileroot);
  childr_compfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(".sachildrcomp")+1,sizeof(char));
  sprintf(childr_compfile,"%s/%s.sachildrcomp",dir,fileroot);

  childx_ptrsfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(".sachildxptrs")+1,sizeof(char));
  sprintf(childx_ptrsfile,"%s/%s.sachildxptrs",dir,fileroot);
  childx_compfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(".sachildxcomp")+1,sizeof(char));
  sprintf(childx_compfile,"%s/%s.sachildxcomp",dir,fileroot);

  pioneerbpfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(".sapiobp")+1,sizeof(char));
  sprintf(pioneerbpfile,"%s/%s.sapiobp",dir,fileroot);
  pior_ptrsfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(".sapiorptrs")+1,sizeof(char));
  sprintf(pior_ptrsfile,"%s/%s.sapiorptrs",dir,fileroot);
  pior_compfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(".sapiorcomp")+1,sizeof(char));
  sprintf(pior_compfile,"%s/%s.sapiorcomp",dir,fileroot);
  piom_ptrsfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(".sapiomptrs")+1,sizeof(char));
  sprintf(piom_ptrsfile,"%s/%s.sapiomptrs",dir,fileroot);
  piom_compfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(".sapiomcomp")+1,sizeof(char));
  sprintf(piom_compfile,"%s/%s.sapiomcomp",dir,fileroot);
#else
  child_ptrsfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(".sachildptrs")+1,sizeof(char));
  sprintf(child_ptrsfile,"%s/%s.sachildptrs",dir,fileroot);
  child_compfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(".sachildcomp")+1,sizeof(char));
  sprintf(child_compfile,"%s/%s.sachildcomp",dir,fileroot);
  nextpfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(".sanextp")+1,sizeof(char));
  sprintf(nextpfile,"%s/%s.sanextp",dir,fileroot);
#endif

  indexi_ptrsfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(".saindexiptrs")+1,sizeof(char));
  sprintf(indexi_ptrsfile,"%s/%s.saindexiptrs",dir,fileroot);
  indexi_compfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(".saindexicomp")+1,sizeof(char));
  sprintf(indexi_compfile,"%s/%s.saindexicomp",dir,fileroot);
  indexj_ptrsfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(".saindexjptrs")+1,sizeof(char));
  sprintf(indexj_ptrsfile,"%s/%s.saindexjptrs",dir,fileroot);
  indexj_compfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(".saindexjcomp")+1,sizeof(char));
  sprintf(indexj_compfile,"%s/%s.saindexjcomp",dir,fileroot);


  if (Access_file_exists_p(sarrayfile) == false) {
    fprintf(stderr,"Suffix array index file %s does not exist\n",sarrayfile);
    new = (T) NULL;

  } else if (Access_file_exists_p(plcp_ptrsfile) == false || Access_file_exists_p(plcp_compfile) == false) {
    fprintf(stderr,"Suffix array plcp file %s does not exist.  The genome may be built using an older suffix array format\n",sarrayfile);
    new = (T) NULL;

  } else {
    new = (T) MALLOC(sizeof(*new));

    if (access == USE_MMAP_PRELOAD) {
      fprintf(stderr,"Pre-loading suffix array...");
      new->array = (UINT4 *) Access_mmap_and_preload(&new->array_fd,&new->array_len,&npages,&seconds,sarrayfile,
						     sizeof(UINT4));
      new->access = MMAPPED;
      comma1 = Genomicpos_commafmt(new->array_len);
      fprintf(stderr,"done (%s bytes)\n",comma1);
      FREE(comma1);
    } else if (access == USE_MMAP_ONLY) {
      new->array = (UINT4 *) Access_mmap(&new->array_fd,&new->array_len,sarrayfile,sizeof(UINT4),/*randomp*/true);
      new->access = MMAPPED;
    } else if (access == USE_ALLOCATE) {
      /* Always memory map suffix array, regardless of access */
      new->array = (UINT4 *) Access_mmap(&new->array_fd,&new->array_len,sarrayfile,sizeof(UINT4),/*randomp*/true);
      new->access = ALLOCATED;
    }

    new->n_plus_one = new->array_len/sizeof(UINT4); /* Should be genomiclength + 1*/
    new->n = new->n_plus_one - 1;

    new->indexi_ptrs = (UINT4 *) Access_allocated(&new->indexi_ptrs_len,&seconds,indexi_ptrsfile,sizeof(UINT4));
    /* 8 is for two DIFFERENTIAL_METAINFO_SIZE words */
    new->indexsize = 3 + log4(((new->indexi_ptrs_len - 8)/sizeof(UINT4))/ /*DIFFERENTIAL_METAINFO_SIZE*/2);
    new->indexspace = power(4,new->indexsize);

    new->indexi_comp = (UINT4 *) Access_allocated(&new->indexi_comp_len,&seconds,indexi_compfile,sizeof(UINT4));
    new->indexj_ptrs = (UINT4 *) Access_allocated(&new->indexj_ptrs_len,&seconds,indexj_ptrsfile,sizeof(UINT4));
    new->indexj_comp = (UINT4 *) Access_allocated(&new->indexj_comp_len,&seconds,indexj_compfile,sizeof(UINT4));


    new->plcp_ptrs = (UINT4 *) Access_allocated(&new->plcp_ptrs_len,&seconds,plcp_ptrsfile,sizeof(UINT4));
    if (access == USE_ALLOCATE) {
      new->plcp_comp = (UINT4 *) Access_allocated(&new->plcp_comp_len,&seconds,plcp_compfile,sizeof(UINT4));
    } else {
      new->plcp_comp = (UINT4 *) Access_mmap(&new->plcp_comp_fd,&new->plcp_comp_len,plcp_compfile,sizeof(UINT4),/*randomp*/true);
    }

#ifdef USE_CHILD_BP
    if (access == USE_ALLOCATE) {
      new->childbp = (UINT4 *) Access_allocated(&new->childbp_len,&seconds,childbpfile,sizeof(UINT4));
    } else {
      new->childbp = (UINT4 *) Access_mmap(&new->childbp_fd,&new->childbp_len,childbpfile,sizeof(UINT4),/*randomp*/true);
    }
    
    /* Always allocate childfc */
    new->childfc = (UINT4 *) Access_allocated(&new->childfc_len,&seconds,childfcfile,sizeof(UINT4));

    if (Access_file_exists_p(childs_pagesfile) == false) {
      new->childs_pages = (UINT4 *) NULL;
    } else {
      new->childs_pages = (UINT4 *) Access_allocated(&new->childs_pages_len,&seconds,childs_pagesfile,sizeof(UINT4));
    }
    new->childs_ptrs = (UINT4 *) Access_allocated(&new->childs_ptrs_len,&seconds,childs_ptrsfile,sizeof(UINT4));
    if (access == USE_ALLOCATE) {
      new->childs_comp = (UINT4 *) Access_allocated(&new->childs_comp_len,&seconds,childs_compfile,sizeof(UINT4));
    } else {
      new->childs_comp = (UINT4 *) Access_mmap(&new->childs_comp_fd,&new->childs_comp_len,childs_compfile,sizeof(UINT4),/*randomp*/true);
    }

    new->childr_ptrs = (UINT4 *) Access_allocated(&new->childr_ptrs_len,&seconds,childr_ptrsfile,sizeof(UINT4));
    if (access == USE_ALLOCATE) {
      new->childr_comp = (UINT4 *) Access_allocated(&new->childr_comp_len,&seconds,childr_compfile,sizeof(UINT4));
    } else {
      new->childr_comp = (UINT4 *) Access_mmap(&new->childr_comp_fd,&new->childr_comp_len,childr_compfile,sizeof(UINT4),/*randomp*/true);
    }

    new->childx_ptrs = (UINT4 *) Access_allocated(&new->childx_ptrs_len,&seconds,childx_ptrsfile,sizeof(UINT4));
    if (access == USE_ALLOCATE) {
      new->childx_comp = (UINT4 *) Access_allocated(&new->childx_comp_len,&seconds,childx_compfile,sizeof(UINT4));
    } else {
      new->childx_comp = (UINT4 *) Access_mmap(&new->childx_comp_fd,&new->childx_comp_len,childx_compfile,sizeof(UINT4),/*randomp*/true);
    }


    if (access == USE_ALLOCATE) {
      new->pioneerbp = (UINT4 *) Access_allocated(&new->pioneerbp_len,&seconds,pioneerbpfile,sizeof(UINT4));
    } else {
      new->pioneerbp = (UINT4 *) Access_mmap(&new->pioneerbp_fd,&new->pioneerbp_len,pioneerbpfile,sizeof(UINT4),/*randomp*/true);
    }

    new->pior_ptrs = (UINT4 *) Access_allocated(&new->pior_ptrs_len,&seconds,pior_ptrsfile,sizeof(UINT4));
    if (access == USE_ALLOCATE) {
      new->pior_comp = (UINT4 *) Access_allocated(&new->pior_comp_len,&seconds,pior_compfile,sizeof(UINT4));
    } else {
      new->pior_comp = (UINT4 *) Access_mmap(&new->pior_comp_fd,&new->pior_comp_len,pior_compfile,sizeof(UINT4),/*randomp*/true);
    }

    new->piom_ptrs = (UINT4 *) Access_allocated(&new->piom_ptrs_len,&seconds,piom_ptrsfile,sizeof(UINT4));
    if (access == USE_ALLOCATE) {
      new->piom_comp = (UINT4 *) Access_allocated(&new->piom_comp_len,&seconds,piom_compfile,sizeof(UINT4));
    } else {
      new->piom_comp = (UINT4 *) Access_mmap(&new->piom_comp_fd,&new->piom_comp_len,piom_compfile,sizeof(UINT4),/*randomp*/true);
    }

    BP_read_setup(BLOCKSIZE,SELECT_SAMPLING_INTERVAL);

#else
    new->nextp = (UINT4 *) Access_allocated(&new->nextp_len,&seconds,nextpfile,sizeof(UINT4));
    
    new->child_ptrs = (UINT4 *) Access_allocated(&new->child_ptrs_len,&seconds,child_ptrsfile,sizeof(UINT4));
    if (access == USE_ALLOCATE) {
      new->child_comp = (UINT4 *) Access_allocated(&new->child_comp_len,&seconds,child_compfile,sizeof(UINT4));
    } else {
      new->child_comp = (UINT4 *) Access_mmap(&new->child_comp_fd,&new->child_comp_len,child_compfile,sizeof(UINT4),/*randomp*/true);
    }
#endif
  }


  FREE(indexj_compfile);
  FREE(indexj_ptrsfile);
  FREE(indexi_compfile);
  FREE(indexi_ptrsfile);

#ifdef USE_CHILD_BP
  FREE(piom_compfile);
  FREE(piom_ptrsfile);
  FREE(pior_compfile);
  FREE(pior_ptrsfile);
  FREE(pioneerbpfile);

  FREE(childx_compfile);
  FREE(childx_ptrsfile);
  FREE(childr_compfile);
  FREE(childr_ptrsfile);
  FREE(childs_compfile);
  FREE(childs_ptrsfile);
  FREE(childs_pagesfile);
  FREE(childfcfile);
  FREE(childbpfile);
#else
  FREE(nextpfile);
  FREE(child_compfile);
  FREE(child_ptrsfile);
#endif

  FREE(plcp_compfile);
  FREE(plcp_ptrsfile);
  FREE(sarrayfile);

  return new;
}


void
Sarray_free (T *old) {
  if (*old) {
    FREE((*old)->indexi_ptrs);
    FREE((*old)->indexi_comp);
    FREE((*old)->indexj_ptrs);
    FREE((*old)->indexj_comp);
    FREE((*old)->plcp_ptrs);
#ifdef USE_CHILD_BP
    FREE((*old)->childfc);
    FREE((*old)->childs_pages);
    FREE((*old)->childs_ptrs);
    FREE((*old)->childr_ptrs);
    FREE((*old)->childx_ptrs);
    FREE((*old)->pior_ptrs);
    FREE((*old)->piom_ptrs);
#else
    FREE((*old)->child_ptrs);
    FREE((*old)->nextp);
#endif

    munmap((void *) (*old)->array,(*old)->array_len);
    close((*old)->array_fd);

    if ((*old)->access == ALLOCATED) {
      FREE((*old)->plcp_comp);
#ifdef USE_CHILD_BP
      FREE((*old)->childbp);
      FREE((*old)->childs_comp);
      FREE((*old)->childr_comp);
      FREE((*old)->childx_comp);
      FREE((*old)->pior_comp);
      FREE((*old)->piom_comp);
#endif

    } else if ((*old)->access == MMAPPED) {
      munmap((void *) (*old)->plcp_comp,(*old)->plcp_comp_len);
      close((*old)->plcp_comp_fd);
#ifdef USE_CHILD_BP
      munmap((void *) (*old)->childs_comp,(*old)->childs_comp_len);
      munmap((void *) (*old)->childr_comp,(*old)->childr_comp_len);
      munmap((void *) (*old)->childx_comp,(*old)->childx_comp_len);
      munmap((void *) (*old)->pior_comp,(*old)->pior_comp_len);
      munmap((void *) (*old)->piom_comp,(*old)->piom_comp_len);
      close((*old)->childs_comp_fd);
      close((*old)->childr_comp_fd);
      close((*old)->childx_comp_fd);
      close((*old)->pior_comp_fd);
      close((*old)->piom_comp_fd);
#endif
    }

    FREE(*old);
  }

  return;
}



#if 0
/* Old search method.  O(m*(log n)), where m is the querylength and n
   is the size of the suffix array searched */
static Sarrayptr_T
sarray_search_init (char *query, int querylength, int queryoffset, Compress_T query_compress, bool plusp,
		    Sarrayptr_T low, Sarrayptr_T high, Univcoord_T nmatches_low, Univcoord_T nmatches_high) {
  Sarrayptr_T mid;
  Univcoord_T pos;
  Univcoord_T nmatches_mid, fasti;
  char c;
  UINT4 sa_low, sa_mid;
  UINT4 lcp_low, lcp_mid;

  assert(querylength > 0);

  debug1(printf("sarray_search_init on querylength %d with low %u, high %u\n",querylength,low,high));
  while (low + 1 < high) {
    /* Compute mid for unsigned ints */
    mid = low/2 + high/2;
    if (low % 2 == 1 && high % 2 == 1) {
      mid += 1;
    }
    debug1(printf("low %u, high %u => mid %u\n",low,high,mid));
    nmatches_mid =  (nmatches_low < nmatches_high) ? nmatches_low : nmatches_high;

    fasti = nmatches_mid +
      (Univcoord_T) Genome_consecutive_matches_rightward(query_compress,/*left*/sarray->array[mid]-queryoffset,
							 /*pos5*/queryoffset+nmatches_mid,
							 /*pos3*/queryoffset+querylength,plusp,/*genestrand*/0);
    pos = sarray->array[mid] + fasti;
    c = Genome_get_char_lex(genome,pos,sarray->n);

    if (fasti == (Univcoord_T) querylength || c > query[fasti]) {
      high = mid;
      /* nmatches_high = (sarray->lcp[mid] < nmatches_mid) ? sarray->lcp[mid] : nmatches_mid; */
      sa_mid = sarray->array[mid];
      lcp_mid = Bitpack64_offsetptr_only(sa_mid,sarray->plcp_ptrs,sarray->plcp_comp) - sa_mid;
#ifdef USE_LCP
      if (lcp_mid != sarray->lcp[mid]) {
	fprintf(stderr,"LCP compression error at %u\n",mid);
      }
#endif
      nmatches_high = (lcp_mid < nmatches_mid) ? lcp_mid : nmatches_mid;
    } else {
      low = mid;
      /* nmatches_low = (sarray->lcp[low] < nmatches_mid) ? sarray->lcp[low] : nmatches_mid; */
      sa_low = sarray->array[low];
      lcp_low = Bitpack64_offsetptr_only(sa_low,sarray->plcp_ptrs,sarray->plcp_comp) - sa_low;
#ifdef USE_LCP
      if (lcp_low != sarray->lcp[low]) {
	fprintf(stderr,"LCP compression error at %u\n",mid);
      }
#endif
      nmatches_low = (lcp_low < nmatches_mid) ? lcp_low : nmatches_mid;
    }

    debug1(printf("sarray_search_init with low %u, high %u\n",low,high));
  }

  debug1(printf("sarray_search_init ended.  Returning low %u+1\n\n",low));
  return low + 1;
}
#endif


#if 0
/* Old search method.  O(m*(log n)), where m is the querylength and n
   is the size of the suffix array searched */
static Sarrayptr_T
sarray_search_final (char *query, int querylength, int queryoffset, Compress_T query_compress, bool plusp,
		     Sarrayptr_T low, Sarrayptr_T high, Univcoord_T nmatches_low, Univcoord_T nmatches_high) {
  Sarrayptr_T mid;
  Univcoord_T pos;
  Univcoord_T nmatches_mid, fasti;
  UINT4 sa_low, sa_mid;
  UINT4 lcp_low, lcp_mid;
  char c;

  assert(querylength > 0);

  debug1(printf("sarray_search_final on querylength %d with low %u, high %u\n",querylength,low,high));
  while (low + 1 < high) {
    /* Compute mid for unsigned ints */
    mid = low/2 + high/2;
    if (low % 2 == 1 && high % 2 == 1) {
      mid += 1;
    }
    debug1(printf("low %u, high %u => mid %u\n",low,high,mid));
    nmatches_mid =  (nmatches_low < nmatches_high) ? nmatches_low : nmatches_high;

    fasti = nmatches_mid +
      (Univcoord_T) Genome_consecutive_matches_rightward(query_compress,/*left*/sarray->array[mid]-queryoffset,
							 /*pos5*/queryoffset+nmatches_mid,
							 /*pos3*/queryoffset+querylength,plusp,/*genestrand*/0);
    pos = sarray->array[mid] + fasti;
    c = Genome_get_char_lex(genome,pos,sarray->n);

    if (fasti == (Univcoord_T) querylength || c < query[fasti]) {
      low = mid;
      /* nmatches_low = (sarray->lcp[low] < nmatches_mid) ? sarray->lcp[low] : nmatches_mid; */
      sa_low = sarray->array[low];
      lcp_low = Bitpack64_offsetptr_only(sa_low,sarray->plcp_ptrs,sarray->plcp_comp) - sa_low;
#ifdef USE_LCP
      if (lcp_low != sarray->lcp[low]) {
	fprintf(stderr,"LCP compression error at %u\n",mid);
      }
#endif
      nmatches_low = (lcp_low < nmatches_mid) ? lcp_low : nmatches_mid;
    } else {
      high = mid;
      /* nmatches_high = (sarray->lcp[mid] < nmatches_mid) ? sarray->lcp[mid] : nmatches_mid; */
      sa_mid = sarray->array[mid];
      lcp_mid = Bitpack64_offsetptr_only(sa_mid,sarray->plcp_ptrs,sarray->plcp_comp) - sa_mid;
#ifdef USE_LCP
      if (lcp_mid != sarray->lcp[mid]) {
	fprintf(stderr,"LCP compression error at %u\n",mid);
      }
#endif
      nmatches_high = (lcp_mid < nmatches_mid) ? lcp_mid : nmatches_mid;
    }

    debug1(printf("sarray_search_final with low %u, high %u\n",low,high));
  }

  debug1(printf("sarray_search_final ended.  Returning high %u-1\n\n",high-1));
  return high - 1;
}
#endif


int
nt_querylength (char *query, int querylength) {
  int i;
  char c;

  i = 0;
  while (i < querylength && ((c = query[i]) == 'A' || c == 'C' || c == 'G' || c == 'T')) {
    i++;
  }

  return i;
}


Storedoligomer_T
nt_oligo (char *query, int indexsize) {
  Storedoligomer_T oligo = 0U;
  int i;

  for (i = 0; i < indexsize; i++) {
    oligo *= 4;
    
    switch (query[i]) {
    case 'A': break;
    case 'C': oligo += 1; break;
    case 'G': oligo += 2; break;
    case 'T': oligo += 3; break;
    default:
      fprintf(stderr,"Saw N in nt_oligo\n");
      abort();
    }
  }

  return oligo;
}

Storedoligomer_T
nt_oligo_truncate (char *query, int truncsize, int indexsize, int subst_value) {
  Storedoligomer_T oligo = 0U;
  int i;

  for (i = 0; i < truncsize; i++) {
    oligo *= 4;
    
    switch (query[i]) {
    case 'A': break;
    case 'C': oligo += 1; break;
    case 'G': oligo += 2; break;
    case 'T': oligo += 3; break;
    default:
      fprintf(stderr,"Saw N in nt_oligo\n");
      abort();
    }
  }

  for ( ; i < indexsize; i++) {
    oligo *= 4;
    oligo += subst_value;
  }

  return oligo;
}


#ifndef USE_CHILD_BP

/* For child[index+1].up, just calling child[index] */
#define decode_up(index,child_ptrs,child_comp) index - Bitpack64_access(index,child_ptrs,child_comp)
#define decode_down(index,child_ptrs,child_comp) Bitpack64_access(index,child_ptrs,child_comp) + index + 1
#define decode_next(index,child_ptrs,child_comp) Bitpack64_access(index,child_ptrs,child_comp) + index + 1


#if 0
/* For benchmarking */
void
Sarray_traverse_children (Sarrayptr_T i, Sarrayptr_T j, T sarray) {
  UINT4 up, nextl;

  /* LCP interval */
  debug1(printf("lcp-interval %u..%u\n",i,j));
  up = decode_up(j,sarray->child_ptrs,sarray->child_comp);
  if (i < up && up <= j) {
    nextl = up;
    debug1(printf("nextl is up: %d\n",nextl));
  } else {
    nextl = decode_down(i,sarray->child_ptrs,sarray->child_comp); /* down */
    debug1(printf("nextl is down: %d\n",nextl));
  }

  /* Test for child[i] being down: lcp[child[i]] > lcp[i] */
  /* Test for child[i] being next_lindex: lcp[child[i]] == lcp[i] */
  while (get_bit(sarray->nextp,nextl) != 0) {
    debug2(printf("Child: %u to %u, char %c\n",nextl,decode_next(nextl,child_ptrs,child_comp)-1,c));
    nextl = decode_next(nextl,sarray->child_ptrs,sarray->child_comp); /* child[nextl] */
  }

  return;
}
#endif


/* Previously did not use this code */
static bool
get_child (Sarrayptr_T *l, Sarrayptr_T *r, Sarrayptr_T i, Sarrayptr_T j, char desired_char,
	   UINT4 *child_ptrs, UINT4 *child_comp, INT4 *nextp, UINT4 *SA, UINT4 *plcp_ptrs, UINT4 *plcp_comp) {
  UINT4 up, nextl;
  Sarrayptr_T sa_nextl;
  UINT4 lcp_whole;
  UINT4 pos;
  char c;

  debug2(printf("Getting children for l-interval from %u to %u, char %c\n",i,j,desired_char));

  /* Test for child[j] being up: lcp[j] > lcp[j+1] */

  up = decode_up(j,child_ptrs,child_comp); /* up: child[j] = childtab[j+1].up */
  if (i < up && up <= j) {
    nextl = up;
  } else {
    nextl = decode_down(i,child_ptrs,child_comp);	/* down: child[i] */
  }
  sa_nextl = SA[nextl];
  lcp_whole = Bitpack64_offsetptr_only(sa_nextl,plcp_ptrs,plcp_comp) - sa_nextl;
  debug2(printf("LCP of whole is %u\n",lcp_whole));

  pos = SA[i] + lcp_whole;
  c = Genome_get_char_lex(genome,pos,sarray->n);
  if (c > desired_char) {
    debug2(printf("Returning false\n"));
    return false;
  } else if (c == desired_char) {
    *l = i;
    *r = nextl - 1;
    debug2(printf("Child: %u to %u, char %c\n",*l,*r,c));
    debug2(printf("Returning true\n\n"));
    return true;
  } else {
    /* Skip */
  }
  
  /* Test for child[i] being down: lcp[child[i]] > lcp[i] */
  /* Test for child[i] being next_lindex: lcp[child[i]] == lcp[i] */
  while (get_bit(nextp,nextl) != 0) {
    pos = SA[nextl] + lcp_whole;
    c = Genome_get_char_lex(genome,pos,sarray->n);
    if (c > desired_char) {
      debug2(printf("Returning false\n"));
      return false;
    } else if (c == desired_char) {
      *l = nextl;
      nextl = decode_next(nextl,child_ptrs,child_comp);  /* child[nextl]; */
      *r = nextl - 1;
      debug2(printf("Child: %u to %u, char %c\n",*l,*r,c));
      debug2(printf("Returning true\n\n"));
      return true;
    } else {
      nextl = decode_next(nextl,child_ptrs,child_comp);  /* child[nextl]; */
    }
  }

  pos = SA[nextl] + lcp_whole;
  c = Genome_get_char_lex(genome,pos,sarray->n);
  if (c == desired_char) {
    *l = nextl;
    *r = j;
    debug2(printf("Child: %u to %u, char %c\n",*l,*r,c));
    debug2(printf("Returning true\n\n"));
    return true;
  } else {
    debug2(printf("Returning false\n"));
    return false;
  }
}


/* Previously used this code.  Avoids recomputing lcp_whole */
static bool
get_child_given_first (Sarrayptr_T *l, Sarrayptr_T *r, Sarrayptr_T i, Sarrayptr_T j, char desired_char,
		       UINT4 *child_ptrs, UINT4 *child_comp, UINT4 *nextp, UINT4 *SA, UINT4 lcp_whole, UINT4 nextl) {
  UINT4 pos;
  char c;

  debug2(printf("Getting children for l-interval from %u to %u, char %c\n",i,j,desired_char));

  /* First child already given */
  pos = SA[i] + lcp_whole;
  c = Genome_get_char_lex(genome,pos,sarray->n);
  if (c > desired_char) {
    debug2(printf("Child: %u to %u, char %c\n",i,nextl-1,c));
    debug2(printf("1.  Returning false, because %c > desired %c\n",c,desired_char));
    return false;
  } else if (c == desired_char) {
    *l = i;
    *r = nextl - 1;
    debug2(printf("Child: %u to %u, char %c\n",i,nextl-1,c));
    debug2(printf("Returning true\n\n"));
    return true;
  } else {
    /* Skip */
    debug2(printf("Child: %u to %u, char %c\n",i,nextl-1,c));
  }
  
  /* Test for child[i] being down: lcp[child[i]] > lcp[i] */
  /* Test for child[i] being next_lindex: lcp[child[i]] == lcp[i] */
  while (get_bit(nextp,nextl) != 0) {
    pos = SA[nextl] + lcp_whole;
    c = Genome_get_char_lex(genome,pos,sarray->n);
    if (c > desired_char) {
      debug2(printf("Child: %u to %u, char %c\n",nextl,decode_next(nextl,child_ptrs,child_comp)-1,c));
      debug2(printf("2.  Returning false, because %c > desired %c\n",c,desired_char));
      return false;
    } else if (c == desired_char) {
      *l = nextl;
      *r = decode_next(nextl,child_ptrs,child_comp) - 1; /* child[nextl] - 1 */
      debug2(printf("Child: %u to %u, char %c\n",nextl,decode_next(nextl,child_ptrs,child_comp)-1,c));
      debug2(printf("Returning true\n\n"));
      return true;
    } else {
      debug2(printf("Child: %u to %u, char %c\n",nextl,decode_next(nextl,child_ptrs,child_comp)-1,c));
      nextl = decode_next(nextl,child_ptrs,child_comp); /* child[nextl] */
    }
  }

  debug2(printf("Processing last interval\n"));
  pos = SA[nextl] + lcp_whole;
  c = Genome_get_char_lex(genome,pos,sarray->n);
  if (c == desired_char) {
    *l = nextl;
    *r = j;
    debug2(printf("Child: %u to %u, char %c\n",nextl,j,c));
    debug2(printf("Returning true\n\n"));
    return true;
  } else {
    debug2(printf("Child: %u to %u, char %c\n",nextl,j,c));
    debug2(printf("3.  Returning false, because %c != desired %c\n",c,desired_char));
    return false;
  }
}


static UINT4
find_longest_match (UINT4 nmatches, Sarrayptr_T *initptr, Sarrayptr_T *finalptr,
		    Sarrayptr_T i, Sarrayptr_T j, char *query, UINT4 querylength,
		    int queryoffset, Compress_T query_compress, T sarray, bool plusp) {
  UINT4 lcp_whole, nextl, up, sa_nextl;
  UINT4 minlength;
  UINT4 l, r;

  while (nmatches < querylength) {
    if (i == j) {
      /* Singleton interval */
      debug1(printf("Singleton interval %u..%u\n",i,j));
      nmatches +=
	Genome_consecutive_matches_rightward(query_compress,/*left*/sarray->array[i]-queryoffset,
					     /*pos5*/queryoffset+nmatches,/*pos3*/queryoffset+querylength,
					     plusp,/*genestrand*/0);
      *initptr = i;
      *finalptr = j;
      return nmatches;

    } else {
      /* LCP interval */
      debug1(printf("lcp-interval %u..%u\n",i,j));
      up = decode_up(j,sarray->child_ptrs,sarray->child_comp);
      if (i < up && up <= j) {
	nextl = up;
	debug1(printf("nextl is up: %u\n",nextl));
      } else {
	nextl = decode_down(i,sarray->child_ptrs,sarray->child_comp); /* down */
	debug1(printf("nextl is down: %u\n",nextl));
      }
      sa_nextl = sarray->array[nextl];
      lcp_whole = Bitpack64_offsetptr_only(sa_nextl,sarray->plcp_ptrs,sarray->plcp_comp) - sa_nextl; /* lcp(i,j) */
      debug1(printf("lcp_whole for %u..%u is %d, compared with nmatches %d\n",i,j,lcp_whole,nmatches));

      if (lcp_whole > nmatches) {
	/* Check only up to minlength, so we validate the entire interval */
	minlength = (lcp_whole < querylength) ? lcp_whole : querylength;
	debug1(printf("Looking up genome for query from %d .. %d - 1\n",nmatches,minlength));
	nmatches +=
	  Genome_consecutive_matches_rightward(query_compress,/*left*/sarray->array[i]-queryoffset,
					       /*pos5*/queryoffset+nmatches,/*pos3*/queryoffset+minlength,
					       plusp,/*genestrand*/0);
	if (nmatches < minlength) {
	  *initptr = i;
	  *finalptr = j;
	  return nmatches;

	} else if (nmatches >= querylength) {
	  debug1(printf("nmatches is now %d >= querylength %d => success\n",nmatches,querylength));
	  *initptr = i;
	  *finalptr = j;
	  return nmatches;
	}
      }
	
      debug1(printf("nmatches is now %d => desired_char is %c\n",nmatches,query[nmatches]));
      if (get_child_given_first(&l,&r,i,j,/*desired_char*/query[nmatches],
				sarray->child_ptrs,sarray->child_comp,sarray->nextp,
				sarray->array,lcp_whole,nextl) == false) {
	*initptr = i;
	*finalptr = j;
	return nmatches;
      } else {
	nmatches += 1;
	i = l;
	j = r;
      }
    }
  }

  *initptr = i;
  *finalptr = j;
  return nmatches;
}




/* Searches using LCP and child arrays.  Should be O(m * |Sigma|),
   where m wis the querylength and |Sigma| is the size of the alphabet
   (4 for DNA) */
/* query is a substring of the original, starting with queryoffset */
static void
sarray_search (Sarrayptr_T *initptr, Sarrayptr_T *finalptr, bool *successp,
	       UINT4 *nmatches, char *query, UINT4 querylength, int queryoffset,
	       Compress_T query_compress, bool plusp) {
  int effective_querylength;	/* length to first N */
  Storedoligomer_T oligo;
  UINT4 l, r;


#ifdef DEBUG
  int k = 0;
  UINT4 recount;
  char Buffer[1000];
  Univcoord_T hit;
  bool failp;
#elif defined(DEBUG1)
  char Buffer[1000];
#endif

  debug(printf("sarray_search on %s, querylength %d\n",query,querylength));

  /* Find initial lcp-interval */
  effective_querylength = nt_querylength(query,querylength);

  *nmatches = 0;
  if (effective_querylength == 0) {
    *initptr = *finalptr = 0;
    *successp = false;
    return;

  } else if (effective_querylength < sarray->indexsize) {
    debug1(printf("string %.*s with effective querylength %d is shorter than indexsize",
		  querylength,query,effective_querylength));

#if 1
    l = 1;
    r = sarray->n;

#else
    /* Try to infer from 12-mer index, but can be tricky when N's are present */
    oligo = nt_oligo_truncate(query,effective_querylength,sarray->indexsize,/*subst_value for A*/0);
    l = Bitpack64_offsetptr_only(oligo,sarray->indexi_ptrs,sarray->indexi_comp);
    debug1(printf(" => oligo %08X",oligo));

    /* Because $ < A, we need to check for this case.  Need to back up just 1. */
    if (l > 1 && sarray->array[l-1] + effective_querylength == sarray->n) {
      debug1(printf(" (backing up one position for l, because at end of genome)"));
      l--;
    }

    /* Add 1 to rollover to next oligo, to handle Ns in genome */
    oligo = nt_oligo_truncate(query,effective_querylength,sarray->indexsize,/*subst_value for T*/3) + 1;
    r = Bitpack64_read_one(oligo*2,sarray->indexij_ptrs,sarray->indexij_comp) - 1;
    debug1(printf(" => ending oligo %08X => r-value %u\n",oligo,r));

    /* Potential bug here if N is present.  Could read dc and if .X, back up a certain number */

    if (oligo == sarray->indexspace) {
      /* We have a poly-T, so we cannot determine r.  For example,
	 TTTTTN has a different r value than TTN. */
      debug1(printf(" but poly-T => 1-letter for T: %u..%u\n",l,r));
      l = sarray->initindexi[3];
      r = sarray->initindexj[3];
      /* Keep nmatches = 0, because there may not be a T in the genome */

    } else {
      r = Bitpack64_offsetptr_only(oligo,sarray->indexi_ptrs,sarray->indexi_comp) - 1;

      /* Because $ < A, we need to check for this case.  Need to back up just 1. */
      debug1(printf(" (checking %u + %d >= %u)",sarray->array[r],effective_querylength,sarray->n));
      if (r > 0 && sarray->array[r] + effective_querylength >= sarray->n) {
	debug1(printf(" (backing up one position for r, because at end of genome)"));
	r--;
      }

      debug1(printf(" and %08X => interval %u..%u (effective_querylength %d)",
		    oligo,l,r,effective_querylength));

      if (l <= r) {
	/* Keep nmatches = 0, since we don't know the value yet */
	debug1(printf(" (good)\n"));
      } else {
#if 0
	/* Did not find a match using saindex, so resort to one letter */
	switch (query[0]) {
	case 'A': l = sarray->initindexi[0]; r = sarray->initindexj[0]; break;
	case 'C': l = sarray->initindexi[1]; r = sarray->initindexj[1]; break;
	case 'G': l = sarray->initindexi[2]; r = sarray->initindexj[2]; break;
	case 'T': l = sarray->initindexi[3]; r = sarray->initindexj[3]; break;
	default: l = 1; r = 0;
	}
	debug1(printf(" (bad) => 1-letter from %c: %u..%u\n",query[0],l,r));
#else
	/* The entire lcp-interval [1,sarray->n] should also work without initindex */
	l = 1;
	r = sarray->n;
	debug1(printf(" (bad) => entire lcp-interval: %u..%u\n",l,r));
#endif
      }
    }
    /* End of code to infer from 12-mers */
#endif

  } else {
    oligo = nt_oligo(query,sarray->indexsize);
    l = Bitpack64_offsetptr_only(oligo,sarray->indexi_ptrs,sarray->indexi_comp);
    r = Bitpack64_offsetptr_only(oligo,sarray->indexj_ptrs,sarray->indexj_comp);
    debug1(printf("string %.*s is equal/longer than indexsize %d => oligo %u => interval %u..%u",
		  querylength,query,sarray->indexsize,oligo,l,r));
    if (l <= r) {
      debug1(printf(" (good)\n"));
      *nmatches = sarray->indexsize;
      /* i = l; */
      /* j = r; */
    } else {
#if 0
      /* Did not find a match using saindex, so resort to one letter */
      switch (query[0]) {
      case 'A': l = sarray->initindexi[0]; r = sarray->initindexj[0]; break;
      case 'C': l = sarray->initindexi[1]; r = sarray->initindexj[1]; break;
      case 'G': l = sarray->initindexi[2]; r = sarray->initindexj[2]; break;
      case 'T': l = sarray->initindexi[3]; r = sarray->initindexj[3]; break;
      default: l = 1; r = 0;
      }
      debug1(printf(" (bad) => 1-letter from %c: %u..%u\n",query[0],l,r));
#else
      /* The entire lcp-interval [1,sarray->n] should also work without initindex */
      l = 1;
      r = sarray->n;
      debug1(printf(" (bad) => entire lcp-interval: %u..%u\n",l,r));
#endif
    }
  }

  if (l > r) {
    /* Did not find a match using saindex or one letter */
    *initptr = l;
    *finalptr = r;
  } else {
    *nmatches = find_longest_match(*nmatches,&(*initptr),&(*finalptr),/*i*/l,/*j*/r,
				   query,querylength,queryoffset,query_compress,sarray,plusp);
  }

  /* Search through suffix tree */
  debug(printf("initptr gets %u, finalptr gets %u\n",*initptr,*finalptr));

  if (*nmatches < querylength) {
    *successp = false;
    debug(printf("%s fail at %d: got %d hits with %d matches:\n",
		 plusp ? "plus" : "minus",queryoffset,(*finalptr - *initptr + 1),*nmatches));
  } else {
    *successp = true;
    debug(printf("%s success at %d: got %d hits with %d matches:\n",
		 plusp ? "plus" : "minus",queryoffset,(*finalptr - *initptr + 1),*nmatches));
  }

#ifdef DEBUG
  failp = false;

  /* Before */
  if (*nmatches > 0 && *initptr > 0U) {
    recount = Genome_consecutive_matches_rightward(query_compress,/*left*/sarray->array[(*initptr)-1]-queryoffset,
						   /*pos5*/queryoffset,/*pos3*/queryoffset+querylength,
						   plusp,/*genestrand*/0);
    printf("%d\t%u\t%u\t",recount,(*initptr)-1,sarray->array[(*initptr)-1] /*+ 1U*/);
    Genome_fill_buffer_simple(genome,sarray->array[(*initptr)-1],recount+1,Buffer);
    printf("%s\n",Buffer);
    if (recount >= *nmatches) {
      printf("querylength is %d\n",querylength);
      printf("false negative: recount %d at %u before init does equal expected nmatches %d\n",
	     recount,sarray->array[(*initptr)-1],*nmatches);
      failp = true;
    }
  }
  printf("\n");


  /* Hits */
  for (k = 0; k < (int) (*finalptr - *initptr + 1) && k < 100; k++) {
    hit = sarray->array[(*initptr)+k];
    recount = Genome_consecutive_matches_rightward(query_compress,/*left*/hit-queryoffset,
						   /*pos5*/queryoffset,/*pos3*/queryoffset+querylength,
						   plusp,/*genestrand*/0);
    printf("%d\t%u\t%u\t",recount,(*initptr)+k,hit /*+ 1U*/);
    Genome_fill_buffer_simple(genome,sarray->array[(*initptr)+k],recount+1,Buffer);
    
    printf("%s\n",Buffer);
    if (recount != *nmatches) {
      printf("querylength is %d\n",querylength);
      printf("false positive: recount %d at %u does not equal expected nmatches %d\n",
	     recount,sarray->array[(*initptr)],*nmatches);
      failp = true;
    }
    /* hits[k] = sarray->array[(*initptr)++]; */
  }


  /* After */
  if (*nmatches > 0) {
    printf("\n");
    recount = Genome_consecutive_matches_rightward(query_compress,/*left*/sarray->array[(*finalptr)+1]-queryoffset,
						   /*pos5*/queryoffset,/*pos3*/queryoffset+querylength,
						   plusp,/*genestrand*/0);
    printf("%d\t%u\t%u\t",recount,(*finalptr)+1,sarray->array[(*finalptr)+1] /*+ 1U*/);
    Genome_fill_buffer_simple(genome,sarray->array[(*finalptr)+1],recount+1,Buffer);
    printf("%s\n",Buffer);
    if (recount >= *nmatches) {
      printf("querylength is %d\n",querylength);
      printf("false negative: recount %d at %u after (*finalptr) does equal expected nmatches %d\n",
	     recount,sarray->array[(*finalptr)+1],*nmatches);
      failp = true;
    }
  }

  if (failp == true) {
    /* Can happen because $ ranks below 0 */
    abort();
  }
#endif

  return;
}

#else

#if 0
/* For benchmarking */
void
Sarray_traverse_children (Sarrayptr_T i, Sarrayptr_T j, T sarray) {
  UINT4 nextl, nsv_i;
  BP_size_t selecti, initw;
  BP_size_t w, b, x;

  /* First child */
  /* Compute NSV(i) (next smallest value) */
  selecti = BP_select(i,sarray->childbp,sarray->childs_pages,sarray->childs_ptrs,sarray->childs_comp,
		      HALF_BLOCKSIZE,BLOCKSIZE,SELECT_SAMPLING_INTERVAL);
  w = BP_find_closeparen(selecti,sarray->childbp,sarray->childx_ptrs,sarray->childx_comp,sarray->pioneerbp,
			 sarray->pior_ptrs,sarray->pior_comp,sarray->piom_ptrs,sarray->piom_comp,BLOCKSIZE);
  nsv_i = BP_rank_open(w,sarray->childbp,sarray->childr_ptrs,sarray->childr_comp,BLOCKSIZE) /*+1*/;
  debug2(printf("NSV(i)-1 = %u, compared with j %u\n",nsv_i,j));

  /* Check for <= and not ==, because of $ at end of genome */
  /* Compare against j, and not j+1, because we computed NSV(i) - 1 */
  if (nsv_i <= j) {
    /* First child: case 1 */
    initw = w-1;
    w = BP_find_openparen(w-1,sarray->childbp,sarray->childx_ptrs,sarray->childx_comp,sarray->pioneerbp,
			  sarray->pior_ptrs,sarray->pior_comp,sarray->piom_ptrs,sarray->piom_comp,BLOCKSIZE);
    nextl = BP_rank_open(w,sarray->childbp,sarray->childr_ptrs,sarray->childr_comp,BLOCKSIZE);
    debug2(printf("First child, case 1: %u\n",nextl));
  } else {
    /* First child: case 2 */
    w = BP_select(j+1,sarray->childbp,sarray->childs_pages,sarray->childs_ptrs,sarray->childs_comp,
		  HALF_BLOCKSIZE,BLOCKSIZE,SELECT_SAMPLING_INTERVAL);
    initw = w-1;
    w = BP_find_openparen(w-1,sarray->childbp,sarray->childx_ptrs,sarray->childx_comp,sarray->pioneerbp,
			  sarray->pior_ptrs,sarray->pior_comp,sarray->piom_ptrs,sarray->piom_comp,BLOCKSIZE);
    nextl = BP_rank_open(w,sarray->childbp,sarray->childr_ptrs,sarray->childr_comp,BLOCKSIZE);
    debug2(printf("First child, case 2: %u\n",nextl));
  }


  debug2(printf("Requesting rank_close of %llu",w));
  b = BP_rank_close(w,sarray->childbp,sarray->childr_ptrs,sarray->childr_comp,BLOCKSIZE);
  debug2(printf("  Got %u\n",b));

  debug2(printf("close paren w is %u => %08X, and b is %u => %08X\n",
		w-1,get_bit(sarray->childbp,w-1),b-1,get_bit(sarray->childfc,b-1)));
  w--;
  b--;

  while (get_bit(sarray->childbp,w) == close_paren && get_bit(sarray->childfc,b) == 0) {
    debug2(printf("Child: %u",nextl));
    x = BP_find_openparen(w,sarray->childbp,sarray->childx_ptrs,sarray->childx_comp,sarray->pioneerbp,
			  sarray->pior_ptrs,sarray->pior_comp,sarray->piom_ptrs,sarray->piom_comp,BLOCKSIZE);
    nextl = BP_rank_open(x,sarray->childbp,sarray->childr_ptrs,sarray->childr_comp,BLOCKSIZE);
    debug2(printf(" to %u, char %c\n",nextl-1,c));

    w--;
    b--;
    debug2(printf("close paren w is %u => %08X, and b is %u => %08X\n",
		  w,get_bit(sarray->childbp,w),b,get_bit(sarray->childfc,b)));
  }

  return;
}
#endif


static bool
get_child_given_first (Sarrayptr_T *l, Sarrayptr_T *r, Sarrayptr_T i, Sarrayptr_T j, char desired_char,
		       T sarray, UINT4 lcp_whole, UINT4 nextl, BP_size_t initw) {
  UINT4 pos;
  char c;
  BP_size_t w, b, x;
#if 0
  BP_size_t v;
#endif


  debug2(printf("Getting children for l-interval from %u to %u, char %c, lcp_whole = %d\n",
		i,j,desired_char,lcp_whole));

  /* First interval already given */
  pos = sarray->array[i] + lcp_whole;
  c = Genome_get_char_lex(genome,pos,sarray->n);
  if (c > desired_char) {
    debug2(printf("Child: %u to %u, char %c\n",i,nextl-1,c));
    debug2(printf("1.  Returning false, because %c > desired %c\n",c,desired_char));
    return false;
  } else if (c == desired_char) {
    *l = i;
    *r = nextl - 1;
    debug2(printf("Child: %u to %u, char %c\n",i,nextl-1,c));
    debug2(printf("Returning true\n\n"));
    return true;
  } else {
    /* Skip */
    debug2(printf("Child: %u to %u, char %c\n",i,nextl-1,c));
  }

  /* Use first child info */
#if 0
  debug2(printf("Requesting select of nextl %u\n",nextl));
  v = BP_select(nextl,sarray->childbp,sarray->childs_pages,sarray->childs_ptrs,sarray->childs_comp,
		HALF_BLOCKSIZE,BLOCKSIZE,SELECT_SAMPLING_INTERVAL);
  debug2(printf("  Got %llu\n",v));

  debug2(printf("Requesting close paren of %llu",v));
  w = BP_find_closeparen(v,sarray->childbp,sarray->childx_ptrs,sarray->childx_comp,sarray->pioneerbp,
			 sarray->pior_ptrs,sarray->pior_comp,sarray->piom_ptrs,sarray->piom_comp,BLOCKSIZE);
  debug2(printf("  Got %llu\n",w));
  assert(w == initw);
#else
  w = initw;
#endif

  debug2(printf("Requesting rank_close of %llu",w));
  b = BP_rank_close(w,sarray->childbp,sarray->childr_ptrs,sarray->childr_comp,BLOCKSIZE);
  debug2(printf("  Got %u\n",b));

  debug2(printf("close paren w is %u => %08X, and b is %u => %08X\n",
		w-1,get_bit(sarray->childbp,w-1),b-1,get_bit(sarray->childfc,b-1)));
  w--;
  b--;

  while (get_bit(sarray->childbp,w) == close_paren && get_bit(sarray->childfc,b) == 0) {
    /* Test l-interval */
    pos = sarray->array[nextl] + lcp_whole;
    c = Genome_get_char_lex(genome,pos,sarray->n);
    if (c > desired_char) {
      debug2(printf("Child: %u to (nextl), char %c\n",nextl,c));
      debug2(printf("2.  Returning false, because %c > desired %c\n",c,desired_char));
      return false;
    } else if (c == desired_char) {
      *l = nextl;
      x = BP_find_openparen(w,sarray->childbp,sarray->childx_ptrs,sarray->childx_comp,sarray->pioneerbp,
			    sarray->pior_ptrs,sarray->pior_comp,sarray->piom_ptrs,sarray->piom_comp,BLOCKSIZE);
      *r = BP_rank_open(x,sarray->childbp,sarray->childr_ptrs,sarray->childr_comp,BLOCKSIZE) - 1; /* nextl - 1 */
      debug2(printf("Child: %u to %u, char %c\n",nextl,*r,c));
      debug2(printf("Returning true\n\n"));
      return true;
    } else {
      debug2(printf("Child: %u",nextl));
      x = BP_find_openparen(w,sarray->childbp,sarray->childx_ptrs,sarray->childx_comp,sarray->pioneerbp,
			    sarray->pior_ptrs,sarray->pior_comp,sarray->piom_ptrs,sarray->piom_comp,BLOCKSIZE);
      nextl = BP_rank_open(x,sarray->childbp,sarray->childr_ptrs,sarray->childr_comp,BLOCKSIZE);
      debug2(printf(" to %u, char %c\n",nextl-1,c));

      w--;
      b--;
      debug2(printf("close paren w is %u => %08X, and b is %u => %08X\n",
		    w,get_bit(sarray->childbp,w),b,get_bit(sarray->childfc,b)));
    }
  }

  /* Last interval */
  debug2(printf("Processing last interval\n"));
  pos = sarray->array[nextl] + lcp_whole;
  c = Genome_get_char_lex(genome,pos,sarray->n);
  if (c == desired_char) {
    *l = nextl;
    *r = j;
    debug2(printf("Child: %u to %u, char %c\n",nextl,j,c));
    debug2(printf("Returning true\n\n"));
    return true;
  } else {
    debug2(printf("Child: %u to %u, char %c\n",nextl,j,c));
    debug2(printf("3.  Returning false, because %c != desired %c\n",c,desired_char));
    return false;
  }


}


/* Searches using LCP and BP arrays.  Should be O(m * |Sigma|),
   where m wis the querylength and |Sigma| is the size of the alphabet
   (4 for DNA) */
/* query is a substring of the original, starting with queryoffset */
static void
sarray_search (Sarrayptr_T *initptr, Sarrayptr_T *finalptr, bool *successp,
	       UINT4 *nmatches, char *query, UINT4 querylength, int queryoffset,
	       Compress_T query_compress, bool plusp) {
  int effective_querylength;	/* length to first N */
  Storedoligomer_T oligo;

  UINT4 sa_i, sa_j_plus_one, sa_nextl;
  int lcp_i, lcp_j_plus_one, lcp_whole;
  UINT4 minlength;
  bool next_interval_p;
  UINT4 i, j;
  UINT4 l, r, nextl, nsv_i;
  BP_size_t selecti, v, w, initw;


#ifdef DEBUG
  int k = 0;
  UINT4 recount;
  char Buffer[1000];
  Univcoord_T hit;
  bool failp;
#elif defined(DEBUG1)
  char Buffer[1000];
#endif

  debug(printf("sarray_search on %s, querylength %d\n",query,querylength));

  /* Find initial lcp-interval */
  effective_querylength = nt_querylength(query,querylength);

  *nmatches = 0;
  if (effective_querylength == 0) {
    *initptr = *finalptr = 0;
    *successp = false;
    return;

  } else if (effective_querylength < sarray->indexsize) {
    debug1(printf("string %.*s is shorter than indexsize",querylength,query));

    oligo = nt_oligo_truncate(query,effective_querylength,sarray->indexsize,/*subst_value for A*/0);
    l = Bitpack64_offsetptr_only(oligo,sarray->indexi_ptrs,sarray->indexi_comp);
    debug1(printf(" => oligo %08X",oligo));

    /* Because $ < A, we need to check for this case.  Need to back up just 1. */
    if (l > 1 && sarray->array[l-1] + effective_querylength == sarray->n) {
      debug1(printf(" (backing up one position, because at end of genome)"));
      l--;
    }

    /* Add 1 to rollover to next oligo, to handle Ns in genome */
    oligo = nt_oligo_truncate(query,effective_querylength,sarray->indexsize,/*subst_value for T*/3) + 1;
    if (oligo == sarray->indexspace) {
      /* We have a poly-T, so we cannot determine r.  For example,
	 TTTTTN has a different r value than TTN. */
      debug1(printf(" but poly-T => 1-letter for T: %u..%u\n",l,r));
      l = sarray->initindexi[3];
      r = sarray->initindexj[3];
      /* Keep nmatches = 0, because there may not be a T in the genome */

    } else {
      r = Bitpack64_offsetptr_only(oligo,sarray->indexi_ptrs,sarray->indexi_comp) - 1;

      /* Because $ < A, we need to check for this case.  Need to back up just 1. */
      debug1(printf(" (checking %u + %d >= %u)",sarray->array[r],effective_querylength,sarray->n));
      if (r > 0 && sarray->array[r] + effective_querylength >= sarray->n) {
	debug1(printf(" (backing up one position, because at end of genome)"));
	r--;
      }

      debug1(printf(" and %08X => interval %u..%u (effective_querylength %d)",
		    oligo,l,r,effective_querylength));

      if (l <= r) {
	/* Keep nmatches = 0, since we don't know the value yet */
	debug1(printf(" (good)\n"));
      } else {
#if 0
	/* Did not find a match using saindex, so resort to one letter */
	switch (query[0]) {
	case 'A': l = sarray->initindexi[0]; r = sarray->initindexj[0]; break;
	case 'C': l = sarray->initindexi[1]; r = sarray->initindexj[1]; break;
	case 'G': l = sarray->initindexi[2]; r = sarray->initindexj[2]; break;
	case 'T': l = sarray->initindexi[3]; r = sarray->initindexj[3]; break;
	default: l = 1; r = 0;
	}
	debug1(printf(" (bad) => 1-letter from %c: %u..%u\n",query[0],l,r));
#else
	/* The entire lcp-interval [1,sarray->n] should also work without initindex */
	l = 1;
	r = sarray->n;
	debug1(printf(" (bad) => entire lcp-interval: %u..%u\n",l,r));
#endif
      }
    }

  } else {
    oligo = nt_oligo(query,sarray->indexsize);
    l = Bitpack64_offsetptr_only(oligo,sarray->indexi_ptrs,sarray->indexi_comp);
    r = Bitpack64_offsetptr_only(oligo,sarray->indexj_ptrs,sarray->indexj_comp);
    debug1(printf("string %.*s is equal/longer than indexsize %d => oligo %u => interval %u..%u",
		  querylength,query,sarray->indexsize,oligo,l,r));
    if (l <= r) {
      debug1(printf(" (good)\n"));
      *nmatches = sarray->indexsize;
      i = l;
      j = r;
    } else {
#if 0
      /* Did not find a match using saindex, so resort to one letter */
      switch (query[0]) {
      case 'A': l = sarray->initindexi[0]; r = sarray->initindexj[0]; break;
      case 'C': l = sarray->initindexi[1]; r = sarray->initindexj[1]; break;
      case 'G': l = sarray->initindexi[2]; r = sarray->initindexj[2]; break;
      case 'T': l = sarray->initindexi[3]; r = sarray->initindexj[3]; break;
      default: l = 1; r = 0;
      }
      debug1(printf(" (bad) => 1-letter from %c: %u..%u\n",query[0],l,r));
#else
      /* The entire lcp-interval [1,sarray->n] should also work without initindex */
      l = 1;
      r = sarray->n;
      debug1(printf(" (bad) => entire lcp-interval: %u..%u\n",l,r));
#endif
    }
  }

  if (l > r) {
    /* Did not find a match using saindex or one letter */
    i = l;
    j = r;
    next_interval_p = false;
  } else {
    next_interval_p = true;
  }

  /* Search through suffix tree */
  while (*nmatches < querylength && next_interval_p == true) {
    i = l;
    j = r;
    if (i == j) {
      /* Singleton interval */
      *nmatches +=
	Genome_consecutive_matches_rightward(query_compress,/*left*/sarray->array[i]-queryoffset,
					     /*pos5*/queryoffset+(*nmatches),/*pos3*/queryoffset+querylength,
					     plusp,/*genestrand*/0);
      next_interval_p = false;

    } else {
      /* LCP interval */
      debug2(printf("Initial i..j is %u..%u\n",i,j));

      /* First child */
      /* Compute NSV(i) (next smallest value) */
      selecti = BP_select(i,sarray->childbp,sarray->childs_pages,sarray->childs_ptrs,sarray->childs_comp,
			  HALF_BLOCKSIZE,BLOCKSIZE,SELECT_SAMPLING_INTERVAL);
      w = BP_find_closeparen(selecti,sarray->childbp,sarray->childx_ptrs,sarray->childx_comp,sarray->pioneerbp,
			     sarray->pior_ptrs,sarray->pior_comp,sarray->piom_ptrs,sarray->piom_comp,BLOCKSIZE);
      nsv_i = BP_rank_open(w,sarray->childbp,sarray->childr_ptrs,sarray->childr_comp,BLOCKSIZE) /*+1*/;
      debug2(printf("NSV(i)-1 = %u, compared with j %u\n",nsv_i,j));

      /* Check for <= and not ==, because of $ at end of genome */
      /* Compare against j, and not j+1, because we computed NSV(i) - 1 */
      if (nsv_i <= j) {
	/* First child: case 1 */
	initw = w-1;
	w = BP_find_openparen(w-1,sarray->childbp,sarray->childx_ptrs,sarray->childx_comp,sarray->pioneerbp,
			      sarray->pior_ptrs,sarray->pior_comp,sarray->piom_ptrs,sarray->piom_comp,BLOCKSIZE);
	nextl = BP_rank_open(w,sarray->childbp,sarray->childr_ptrs,sarray->childr_comp,BLOCKSIZE);
	debug2(printf("First child, case 1: %u\n",nextl));
      } else {
	/* First child: case 2 */
	w = BP_select(j+1,sarray->childbp,sarray->childs_pages,sarray->childs_ptrs,sarray->childs_comp,
		      HALF_BLOCKSIZE,BLOCKSIZE,SELECT_SAMPLING_INTERVAL);
	initw = w-1;
	w = BP_find_openparen(w-1,sarray->childbp,sarray->childx_ptrs,sarray->childx_comp,sarray->pioneerbp,
			      sarray->pior_ptrs,sarray->pior_comp,sarray->piom_ptrs,sarray->piom_comp,BLOCKSIZE);
	nextl = BP_rank_open(w,sarray->childbp,sarray->childr_ptrs,sarray->childr_comp,BLOCKSIZE);
	debug2(printf("First child, case 2: %u\n",nextl));
      }

      sa_nextl = sarray->array[nextl];
      lcp_whole = Bitpack64_offsetptr_only(sa_nextl,sarray->plcp_ptrs,sarray->plcp_comp) - sa_nextl; /* lcp(i,j) */
      /* printf("lcp_whole for %u..%u is %d\n",i,j,lcp_whole); */

      /* Check only up to minlength, so we validate the entire interval */
      minlength = (lcp_whole < querylength) ? lcp_whole : querylength;
      *nmatches +=
	Genome_consecutive_matches_rightward(query_compress,/*left*/sarray->array[i]-queryoffset,
					     /*pos5*/queryoffset+(*nmatches),/*pos3*/queryoffset+minlength,
					     plusp,/*genestrand*/0);
      if (*nmatches < minlength) {
	next_interval_p = false;

      } else if (*nmatches >= querylength) {
	debug1(printf("nmatches is now %d >= querylength %d => success\n",*nmatches,querylength));
	
      } else {
	debug1(printf("nmatches is now %d => desired_char is %c\n",*nmatches,query[minlength]));
	next_interval_p =
	  get_child_given_first(&l,&r,i,j,/*desired_char*/query[minlength],sarray,lcp_whole,nextl,initw);
      }
    }
  }

  debug(printf("initptr gets %u, finalptr gets %u\n",i,j));

  *initptr = i;
  *finalptr = j;
  if (*nmatches < querylength) {
    *successp = false;
    debug(printf("%s fail at %d: got %d hits with %d matches:\n",
		 plusp ? "plus" : "minus",queryoffset,(*finalptr - *initptr + 1),*nmatches));
  } else {
    *successp = true;
    debug(printf("%s success at %d: got %d hits with %d matches:\n",
		 plusp ? "plus" : "minus",queryoffset,(*finalptr - *initptr + 1),*nmatches));
  }

#ifdef DEBUG
  failp = false;

  /* Before */
  if (*nmatches > 0 && *initptr > 0U) {
    recount = Genome_consecutive_matches_rightward(query_compress,/*left*/sarray->array[(*initptr)-1]-queryoffset,
						   /*pos5*/queryoffset,/*pos3*/queryoffset+querylength,
						   plusp,/*genestrand*/0);
    printf("%d\t%u\t%u\t",recount,(*initptr)-1,sarray->array[(*initptr)-1] /*+ 1U*/);
    Genome_fill_buffer_simple(genome,sarray->array[(*initptr)-1],recount+1,Buffer);
    printf("%s\n",Buffer);
    if (recount >= *nmatches) {
      printf("querylength is %d\n",querylength);
      printf("false negative: recount %d at %u before init does equal expected nmatches %d\n",
	     recount,sarray->array[(*initptr)-1],*nmatches);
      failp = true;
    }
  }
  printf("\n");


  /* Hits */
  for (k = 0; k < (int) (*finalptr - *initptr + 1) && k < 100; k++) {
    hit = sarray->array[(*initptr)+k];
    recount = Genome_consecutive_matches_rightward(query_compress,/*left*/hit-queryoffset,
						   /*pos5*/queryoffset,/*pos3*/queryoffset+querylength,
						   plusp,/*genestrand*/0);
    printf("%d\t%u\t%u\t",recount,(*initptr)+k,hit /*+ 1U*/);
    Genome_fill_buffer_simple(genome,sarray->array[(*initptr)+k],recount+1,Buffer);
    
    printf("%s\n",Buffer);
    if (recount != *nmatches) {
      printf("querylength is %d\n",querylength);
      printf("false positive: recount %d at %u does not equal expected nmatches %d\n",
	     recount,sarray->array[(*initptr)],*nmatches);
      failp = true;
    }
    /* hits[k] = sarray->array[(*initptr)++]; */
  }


  /* After */
  if (*nmatches > 0) {
    printf("\n");
    recount = Genome_consecutive_matches_rightward(query_compress,/*left*/sarray->array[(*finalptr)+1]-queryoffset,
						   /*pos5*/queryoffset,/*pos3*/queryoffset+querylength,
						   plusp,/*genestrand*/0);
    printf("%d\t%u\t%u\t",recount,(*finalptr)+1,sarray->array[(*finalptr)+1] /*+ 1U*/);
    Genome_fill_buffer_simple(genome,sarray->array[(*finalptr)+1],recount+1,Buffer);
    printf("%s\n",Buffer);
    if (recount >= *nmatches) {
      printf("querylength is %d\n",querylength);
      printf("false negative: recount %d at %u after (*finalptr) does equal expected nmatches %d\n",
	     recount,sarray->array[(*finalptr)+1],*nmatches);
      failp = true;
    }
  }

  if (failp == true) {
    /* Can happen because $ ranks below 0 */
    abort();
  }
#endif

  return;
}

#endif




/* Simplified version of Spanningelt_T */
typedef struct Elt_T *Elt_T;
struct Elt_T {
  int querystart;
  int queryend;
  Univcoord_T nmatches;

  Sarrayptr_T initptr;			/* in sarray */
  Sarrayptr_T finalptr;

  Univcoord_T *positions_allocated;
  Univcoord_T *positions;
  int npositions;		/* from goal to high */
  bool filledp;			/* for development purposes */
};


static Elt_T
Elt_new (int querypos, int nmatches, Sarrayptr_T initptr, Sarrayptr_T finalptr) {
  Elt_T new = (Elt_T) MALLOC(sizeof(*new));

  new->querystart = querypos;
  new->queryend = querypos + nmatches - 1;
  new->nmatches = nmatches;

  new->initptr = initptr;
  new->finalptr = finalptr;

  new->positions_allocated = new->positions = (Univcoord_T *) NULL;
  new->npositions = 0;

  new->filledp = false;

  return new;
}

static void
Elt_replace (Elt_T this, int querypos, int nmatches, Sarrayptr_T initptr, Sarrayptr_T finalptr) {
  this->querystart = querypos;
  this->queryend = querypos + nmatches - 1;
  this->nmatches = nmatches;

  this->initptr = initptr;
  this->finalptr = finalptr;

  if (this->positions_allocated != NULL) {
    FREE(this->positions_allocated);
  }
  this->positions_allocated = this->positions = (Univcoord_T *) NULL;
  this->npositions = 0;

  this->filledp = false;

  return;
}


static void
Elt_free (Elt_T *old) {

  if ((*old)->positions_allocated != NULL) {
    FREE((*old)->positions_allocated);
  }
  FREE(*old);
  return;
}


#if 0
static int
Elt_nmatches_cmp (const void *a, const void *b) {
  Elt_T x = * (Elt_T *) a;
  Elt_T y = * (Elt_T *) b;

  if (x->nmatches > y->nmatches) {
    return -1;
  } else if (y->nmatches > x->nmatches) {
    return +1;
  } else {
    return 0;
  }
}
#endif

static int
Elt_querypos_ascending_cmp (const void *a, const void *b) {
  Elt_T x = * (Elt_T *) a;
  Elt_T y = * (Elt_T *) b;

  if (x->querystart < y->querystart) {
    return -1;
  } else if (y->querystart < x->querystart) {
    return +1;
  } else {
    return 0;
  }
}

static int
Elt_querypos_descending_cmp (const void *a, const void *b) {
  Elt_T x = * (Elt_T *) a;
  Elt_T y = * (Elt_T *) b;

  if (x->querystart > y->querystart) {
    return -1;
  } else if (y->querystart > x->querystart) {
    return +1;
  } else {
    return 0;
  }
}


static void
Elt_fill_positions_all (Elt_T this, T sarray) {
  Sarrayptr_T ptr;
  Univcoord_T pos;
  int i;

  if (this->positions_allocated == NULL) {
    this->npositions = this->finalptr - this->initptr + 1;
    if (this->nmatches == 0 || this->npositions > EXCESS_SARRAY_HITS) {
      this->positions_allocated = this->positions = (Univcoord_T *) NULL;
      this->npositions = 0;
    } else {
      this->positions_allocated = this->positions = (Univcoord_T *) CALLOC(this->npositions,sizeof(Univcoord_T));
      i = 0;
      ptr = this->initptr;
      while (ptr <= this->finalptr) {
	if ((pos = sarray->array[ptr++]) >= (Univcoord_T) this->querystart) {
	  this->positions[i++] = pos - this->querystart;
	}
      }
      this->npositions = i;
      qsort(this->positions,this->npositions,sizeof(Univcoord_T),Univcoord_compare);
    }

    this->filledp = true;
  }

  return;
}


#ifdef DEBUG7
static void
print_vector (__m128i x, char *label) {
  __m128i a[1];
  unsigned int *s = a;

  _mm_store_si128(a,x);
  printf("%s: %u %u %u %u\n",label,s[0],s[1],s[2],s[3]);
  return;
}

static void
print_vector_looking (__m128i x, Univcoord_T low, Univcoord_T high) {
  __m128i a[1];
  unsigned int *s = a;

  _mm_store_si128(a,x);
  printf("Looking at value %u, relative to low %u and high %u\n",s[0],low,high);
  printf("Looking at value %u, relative to low %u and high %u\n",s[1],low,high);
  printf("Looking at value %u, relative to low %u and high %u\n",s[2],low,high);
  printf("Looking at value %u, relative to low %u and high %u\n",s[3],low,high);
  return;
}
#endif


#ifdef DEBUG8
/* Non-SIMD methods for comparison */
static void
positions_compare (Univcoord_T *positions, int npositions,
		   Univcoord_T *positions_std, int npositions_std) {
  int i;

  if (npositions != npositions_std) {
    fprintf(stderr,"npositions %d != npositions_std %d\n",npositions,npositions_std);
    abort();
  } else {
    for (i = 0; i < npositions; i++) {
      if (positions[i] != positions_std[i]) {
	fprintf(stderr,"At %d, positions %d != positions_std %d\n",i,positions[i],positions_std[i]);
	abort();
      }
    }
  }

  return;
}

static Univcoord_T *
fill_positions_std (int *npositions, Univcoord_T low_adj, Univcoord_T high_adj,
		    Sarrayptr_T initptr, Sarrayptr_T finalptr,
		    int querystart, Univcoord_T *array) {
  Univcoord_T *more_positions;
  Univcoord_T *positions, value;
  Sarrayptr_T ptr, lastptr;
  int i;

  positions = (Univcoord_T *) CALLOC(GUESS_ALLOCATION,sizeof(Univcoord_T));

  *npositions = 0;
  ptr = initptr;      

  while (ptr <= finalptr) {
    debug7a(printf("Std: Looking at value %u, relative to low %u and high %u\n",array[ptr],low_adj,high_adj));
    if ((value = array[ptr++]) < low_adj) {
      /* Skip */
    } else if (value > high_adj) {
      /* Skip */
    } else if (*npositions < GUESS_ALLOCATION) {
      debug7(printf("Std: Found position %u between low %u and high %u, and within allocation\n",value,low_adj,high_adj));
      positions[(*npositions)++] = value - querystart;
    } else {
      debug7(printf("Std: Found position %u between low %u and high %u, but exceeds allocation\n",value,low_adj,high_adj));
      (*npositions)++;
      lastptr = ptr;		/* saves us from going through the entire sarray below */
    }
  }

  debug7(printf("Std method found %d positions\n",*npositions));
  if (*npositions > GUESS_ALLOCATION) {
    /* Copy the positions we have stored so far */
    more_positions = (Univcoord_T *) CALLOC(*npositions,sizeof(Univcoord_T));
    memcpy(more_positions,positions,GUESS_ALLOCATION*sizeof(Univcoord_T));
    FREE(positions);
    positions = more_positions;
    
    i = GUESS_ALLOCATION;	/* Start count with the number stored */
    ptr = lastptr;	/* One past the last ptr with a result */

    while (i < *npositions) {
      if ((value = array[--ptr]) < low_adj) {
	/* Skip */
      } else if (value > high_adj) {
	/* Skip */
      } else {
	positions[i++] = value - querystart;
      }
    }
  }

  return positions;
}
#endif



static void
Elt_fill_positions_filtered (Elt_T this, T sarray, Univcoord_T goal, Univcoord_T low, Univcoord_T high,
			     Compress_T query_compress, bool plusp, int genestrand) {
  Sarrayptr_T ptr, lastptr;
  int nmatches;
  int i;
  Univcoord_T *array = sarray->array, low_adj, high_adj, value;
  Univcoord_T *more_positions;
#ifdef HAVE_SSE2
#ifdef HAVE_64_BIT
  UINT8 pointer;
#else
  UINT4 pointer;
#endif
  __m128i base, floor, ceiling, values, compare;
  int n_prealign, k;
#endif
#ifdef DEBUG8
  Univcoord_T *positions_std;
  int npositions_std;
#endif


  debug7(printf("Entered Elt_fill_positions_filtered with goal %u, low %u and high %u, initptr %u and finalptr %u (n = %d), nmatches %d\n",
		goal,low,high,this->initptr,this->finalptr,this->finalptr - this->initptr + 1,this->nmatches));
  
  if (this->positions_allocated != NULL) {
    /* Filled from a previous call */
    FREE(this->positions_allocated);
  }

  if (this->nmatches == 0 || this->finalptr - this->initptr + 1 > EXCESS_SARRAY_HITS) {
    nmatches = Genome_consecutive_matches_rightward(query_compress,/*left*/goal,/*pos5*/this->querystart,
						    /*pos3*/this->queryend + 1,plusp,genestrand);
    debug7(printf("rightward at goal %u from %d to %d shows %d matches (want %d)\n",goal,this->querystart,this->queryend,
		  nmatches,this->queryend - this->querystart + 1));
    if (nmatches == this->queryend - this->querystart + 1) {
      /* Create a position that works */
      this->positions_allocated = this->positions = (Univcoord_T *) CALLOC(1,sizeof(Univcoord_T));
      this->positions[0] = goal;
      this->npositions = 1;
    } else {
      this->positions_allocated = this->positions = (Univcoord_T *) NULL;
      this->npositions = 0;
    }
  } else {

#ifdef DEBUG8
    positions_std = fill_positions_std(&npositions_std,/*low_adj*/low + this->querystart,
				       /*high_adj*/high + this->querystart,
				       this->initptr,this->finalptr,this->querystart,array);
#endif


#ifdef HAVE_SSE2
    base = _mm_set1_epi32(2147483648); /* 2^31 */
#endif

    /* Guess at allocation size */
    this->positions_allocated = this->positions = (Univcoord_T *) CALLOC(GUESS_ALLOCATION,sizeof(Univcoord_T));

    low_adj = low + this->querystart;
    high_adj = high + this->querystart;

    this->npositions = 0;
    ptr = this->initptr;      
#ifdef HAVE_SSE2
    if (ptr + 3 > this->finalptr) { /* ptr + 4 > (this->finalptr + 1) */
      /* Handle in normal manner */
      debug7(printf("Small batch, because %u + 3 <= %u\n",ptr,this->finalptr));
      while (ptr <= this->finalptr) {
	debug7a(printf("Looking at value %u, relative to low %u and high %u\n",array[ptr],low_adj,high_adj));
	if ((value = array[ptr++]) < low_adj) {
	  /* Skip */
	} else if (value > high_adj) {
	  /* Skip */
	} else if (this->npositions < GUESS_ALLOCATION) {
	  debug7(printf("Found position %u between low %u and high %u, and within allocation\n",value,low_adj,high_adj));
	  this->positions[this->npositions++] = value - this->querystart;
	} else {
	  debug7(printf("Found position %u between low %u and high %u, but exceeds allocation\n",value,low_adj,high_adj));
	  this->npositions++;
	  lastptr = ptr;		/* saves us from going through the entire sarray below */
	}
      }

    } else {
#ifdef HAVE_64_BIT
      pointer = (UINT8) &(array[ptr]);
#else
      pointer = (UINT4) &(array[ptr]);
#endif
      n_prealign = ((16 - (pointer & 0xF))/4) & 0x3;
      debug7(printf("Initial ptr is at location %p.  Need %d to get to 128-bit boundary\n",
		    &(array[ptr]),n_prealign));

      /* Initial part */
      debug7(printf("Initial part:\n"));
      for (k = 0; k < n_prealign; k++) {
	debug7a(printf("Looking at value %u, relative to low %u and high %u\n",array[ptr],low_adj,high_adj));
	if ((value = array[ptr++]) < low_adj) {
	  /* Skip */
	} else if (value > high_adj) {
	  /* Skip */
	} else if (this->npositions < GUESS_ALLOCATION) {
	  debug7(printf("Found position %u between low %u and high %u, and within allocation\n",value,low_adj,high_adj));
	  this->positions[this->npositions++] = value - this->querystart;
	} else {
	  debug7(printf("Found position %u between low %u and high %u, but exceeds allocation\n",value,low_adj,high_adj));
	  this->npositions++;
	  lastptr = ptr;		/* saves us from going through the entire sarray below */
	}
      }

      /* Aligned part */
      debug7(printf("\nAligned part:\n"));
      /* Since compare operations not available for unsigned ints, using the fact that
	 unsigned_gt(a,b) is equivalent to signed_gt(a - 2^31, b - 2^31) */
      floor = _mm_set1_epi32(low_adj - 1 - 2147483648);
      ceiling = _mm_set1_epi32(high_adj + 1 - 2147483648);
      while (ptr + 3 <= this->finalptr) { /* ptr + 4 < this->finalptr + 1 */
	values = _mm_load_si128((__m128i *) &(array[ptr]));
	debug7a(print_vector_looking(values,low_adj,high_adj));
	values = _mm_sub_epi32(values,base);
	compare = _mm_and_si128(_mm_cmpgt_epi32(values,floor),_mm_cmplt_epi32(values,ceiling));
	if (/*cmp*/_mm_movemask_epi8(compare) == 0x0000) {
	  /* All results are false, indicating no values between low_adj and high_adj (most common case) */
	  ptr += 4;
	} else {
	  for (k = 0; k < 4; k++) {
	    if ((value = array[ptr++]) < low_adj) {
	      /* Skip */
	      debug7(printf("Skipping position %u < low %u\n",value,low_adj));
	    } else if (value > high_adj) {
	      /* Skip */
	      debug7(printf("Skipping position %u > high %u\n",value,high_adj));
	    } else if (this->npositions < GUESS_ALLOCATION) {
	      debug7(printf("Found position %u between low %u and high %u, and within allocation\n",value,low_adj,high_adj));
	      this->positions[this->npositions++] = value - this->querystart;
	    } else {
	      debug7(printf("Found position %u between low %u and high %u, but exceeds allocation\n",value,low_adj,high_adj));
	      this->npositions++;
	      lastptr = ptr;		/* saves us from going through the entire sarray below */
	    }
	  }
	}
      }

      /* Final part */
      debug7(printf("\nFinal part:\n"));
      while (ptr <= this->finalptr) {
	debug7a(printf("Looking at value %u, relative to low %u and high %u\n",array[ptr],low_adj,high_adj));
	if ((value = array[ptr++]) < low_adj) {
	  /* Skip */
	} else if (value > high_adj) {
	  /* Skip */
	} else if (this->npositions < GUESS_ALLOCATION) {
	  debug7(printf("Found position %u between low %u and high %u, and within allocation\n",value,low_adj,high_adj));
	  this->positions[this->npositions++] = value - this->querystart;
	} else {
	  debug7(printf("Found position %u between low %u and high %u, but exceeds allocation\n",value,low_adj,high_adj));
	  this->npositions++;
	  lastptr = ptr;		/* saves us from going through the entire sarray below */
	}
      }
    }

#else

    while (ptr <= this->finalptr) {
      debug7a(printf("Looking at value %u, relative to low %u and high %u\n",array[ptr],low_adj,high_adj));
      if ((value = array[ptr++]) < low_adj) {
	/* Skip */
      } else if (value > high_adj) {
	/* Skip */
      } else if (this->npositions < GUESS_ALLOCATION) {
	debug7(printf("Found position %u between low %u and high %u, and within allocation\n",value,low_adj,high_adj));
	this->positions[this->npositions++] = value - this->querystart;
      } else {
	debug7(printf("Found position %u between low %u and high %u, but exceeds allocation\n",value,low_adj,high_adj));
	this->npositions++;
	lastptr = ptr;		/* saves us from going through the entire sarray below */
      }
    }
#endif

    debug7(printf("SIMD method found %d positions\n",this->npositions));
    if (this->npositions > GUESS_ALLOCATION) {
      /* Handle the case if we exceeded GUESS_ALLOCATION */

      /* Copy the positions we have stored so far */
      more_positions = (Univcoord_T *) CALLOC(this->npositions,sizeof(Univcoord_T));
      memcpy(more_positions,this->positions,GUESS_ALLOCATION*sizeof(Univcoord_T));
      FREE(this->positions_allocated);
      this->positions_allocated = this->positions = more_positions;

      i = GUESS_ALLOCATION;	/* Start count with the number stored */
      ptr = lastptr;		/* One past the last ptr with a result */
#ifdef HAVE_SSE2
      if (this->initptr + 4 < ptr) {
	while (i < this->npositions) {
	  if ((value = array[--ptr]) < low_adj) {
	    /* Skip */
	  } else if (value > high_adj) {
	    /* Skip */
	  } else {
	    this->positions[i++] = value - this->querystart;
	  }
	}

      } else {
#ifdef HAVE_64_BIT
	pointer = (UINT8) &(array[ptr]);
#else
	pointer = (UINT4) &(array[ptr]);
#endif
	n_prealign = ((pointer & 0xF)/4) & 0x3;
	debug7(printf("Initial ptr is at location %p.  Need %d to get to 128-bit boundary\n",
		      &(array[ptr]),n_prealign));

	/* Initial part */
	while (i < this->npositions) {
	  if ((value = array[--ptr]) < low_adj) {
	    /* Skip */
	  } else if (value > high_adj) {
	    /* Skip */
	  } else {
	    this->positions[i++] = value - this->querystart;
	  }
	}

	/* Aligned part */
	while (i < this->npositions && this->initptr + 4 < ptr) {
	  values = _mm_load_si128((__m128i *) &(array[ptr-4]));
	  values = _mm_sub_epi32(values,base);
	  compare = _mm_and_si128(_mm_cmpgt_epi32(values,floor),_mm_cmplt_epi32(values,ceiling));
	  if (/*cmp*/_mm_movemask_epi8(compare) == 0x0000) {
	    /* All results are false, indicating no values between low_adj and high_adj (most common case) */
	    ptr -= 4;
	  } else {
	    for (k = 0; k < 4; k++) {
	      if ((value = array[--ptr]) < low_adj) {
		/* Skip */
	      } else if (value > high_adj) {
		/* Skip */
	      } else {
		this->positions[i++] = value - this->querystart;
	      }
	    }
	  }
	}
	  
	/* Last part */
	while (i < this->npositions) {
	  if ((value = array[--ptr]) < low_adj) {
	    /* Skip */
	  } else if (value > high_adj) {
	    /* Skip */
	  } else {
	    this->positions[i++] = value - this->querystart;
	  }
	}
      }

#else

      while (i < this->npositions) {
	if ((value = array[--ptr]) < low_adj) {
	  /* Skip */
	} else if (value > high_adj) {
	  /* Skip */
	} else {
	  this->positions[i++] = value - this->querystart;
	}
      }
#endif
    }

#ifdef DEBUG8
    positions_compare(this->positions,this->npositions,positions_std,npositions_std);
    FREE(positions_std);
#endif

    qsort(this->positions,this->npositions,sizeof(Univcoord_T),Univcoord_compare);
    debug7(printf("Sorting %d positions\n",this->npositions));

    /* Advance pointer to goal (note: do not want goal_adj, since we have already subtracted this->querystart) */
    /* Have tested positions[i] <= goal, but want positions[-1] to be < goal, or positions[0] >= goal */
    i = 0;
    while (i < this->npositions && this->positions[i] < goal) {
      debug7(printf("Skipping position %u < goal %u\n",this->positions[i],goal));
      i++;
    }
    this->positions += i;
    this->npositions -= i;
    debug7(printf("Remaining: %d positions\n",this->npositions));
  }

  this->filledp = true;

  return;
}
  

static void
Elt_dump_list (List_T list) {
  List_T p;
  Elt_T elt;
  int maxn = 0, k;

  for (p = list; p != NULL; p = p->rest) {
    elt = (Elt_T) p->first;
    if (elt->npositions > maxn) {
      maxn = elt->npositions;
    }
  }

  for (k = 0; k < maxn /* && k < 100 */; k++) {
    for (p = list; p != NULL; p = p->rest) {
      elt = (Elt_T) p->first;
      if (k >= elt->npositions) {
	printf("\t");
      } else {
	printf("%d..%d:%u\t",elt->querystart,elt->queryend,elt->positions[k]);
      }
    }
    printf("\n");
  }
  printf("\n");

  return;
}

static void
Elt_dump (Elt_T elt) {
  int k;

  for (k = 0; k < elt->npositions; k++) {
    printf("%d..%d:%u\n",elt->querystart,elt->queryend,elt->positions[k]);
  }
  printf("\n");

  return;
}



static int
binary_search (int lowi, int highi, Univcoord_T *positions, Univcoord_T goal) {
  int middlei;

  debug10(printf("entered binary search with lowi=%d, highi=%d, goal=%u\n",lowi,highi,goal));

  while (lowi < highi) {
    middlei = lowi + ((highi - lowi) / 2);
    debug10(printf("  binary: %d:%u %d:%u %d:%u   vs. %u\n",
		   lowi,positions[lowi],middlei,positions[middlei],
		   highi,positions[highi],goal));
    if (goal < positions[middlei]) {
      highi = middlei;
    } else if (goal > positions[middlei]) {
      lowi = middlei + 1;
    } else {
      debug10(printf("binary search returns %d\n",middlei));
      return middlei;
    }
  }

  debug10(printf("binary search returns %d\n",highi));
  return highi;
}


#define add_bounded(x,plusterm,highbound) ((x + (plusterm) >= highbound) ? (highbound - 1) : x + (plusterm))
#define subtract_bounded(x,minusterm,lowbound) ((x < lowbound + (minusterm)) ? lowbound : x - (minusterm))


/* Taken from stage1hr.c identify_multimiss_iter */
static bool
extend_rightward (Univcoord_T goal, Univcoord_T chroffset, Univcoord_T chrhigh,
		  List_T set, T sarray, Compress_T query_compress,
		  bool plusp, int genestrand, int best_queryend) {
  Elt_T elt;
  Univcoord_T low, high;

  debug7(printf("extend_rightward, with goal %u\n",goal));

  for ( ; set /* != NULL */; set = set->rest) {
    debug7(Elt_dump_list(set));
    elt = (Elt_T) set->first;

    debug7(printf("remaining elts %d: ",List_length(set)));
    debug7(printf("%d..%d\n",elt->querystart,elt->queryend));
    if (elt->querystart > best_queryend) {
      /* Allow for deletion with higher goal */
      low = subtract_bounded(goal,/*minusterm*/max_insertionlen,chroffset);
      high = add_bounded(goal,/*plusterm*/overall_max_distance,chrhigh);
      Elt_fill_positions_filtered(elt,sarray,goal,low,high,query_compress,plusp,genestrand);
      debug7(printf("Allow for deletion with higher goal: %d positions\n",elt->npositions));

      if (elt->npositions <= 0) {
	/* List is empty, so one more miss seen. */
	debug7(printf(" positions empty, so not spanning\n"));
	return false;
	
      } else if (*elt->positions > high) {
	/* Already advanced past goal, so one more miss seen. */
	debug7(printf(" %u advanced past goal %u + %d, so not spanning\n",*elt->positions,goal,overall_max_distance));
	return false;

      } else {
	/* Found goal.  Advance past goal and continue with loop. */
	debug7(printf(" advancing\n"));
	++elt->positions;
	--elt->npositions;
	/* continue */
      }
    } else {
      /* Allow for deletion with lower goal */
      low = subtract_bounded(goal,/*minusterm*/overall_max_distance,chroffset);
      high = add_bounded(goal,/*plusterm*/max_insertionlen,chrhigh);
      Elt_fill_positions_filtered(elt,sarray,goal,low,high,query_compress,plusp,genestrand);
      debug7(printf("Allow for deletion with lower goal: %d positions\n",elt->npositions));

      if (elt->npositions <= 0) {
	/* List is empty, so test previous one only, which must exist
	   since positions had at least one entry. */
	if (elt->positions[-1] >= low) {
	  /* Found goal with deletion */
	  debug7(printf(" possible deletion, continuing\n"));
	  /* continue */
	} else {
	  debug7(printf(" previous %u before goal %u - %d, so not spanning\n",elt->positions[-1],goal,shortsplicedist));
	  return false;
	}
	
      } else if (elt->positions == elt->positions_allocated) {
	/* List is at beginning, so test current one only, not the previous one */
	if (*elt->positions > goal) {
	  /* Already advanced past goal, so one more miss seen. */
	  debug7(printf(" %u advanced past goal_high %u, so not spanning\n",*elt->positions,goal));
	  return false;
	} else {
	  /* Found goal.  Advance past goal and continue with loop. */
	  debug7(printf(" advancing\n"));
	  ++elt->positions;
	  --elt->npositions;
	  /* continue */
	}

      } else {
	/* Test both current one (for goal) and previous one (for deletion) */
	if (*elt->positions == goal) {
	  /* Found goal.  Advance past goal and continue with loop. */
	  debug7(printf(" advancing\n"));
	  ++elt->positions;
	  --elt->npositions;
	  /* continue */

	} else if (elt->positions[-1] >= low) {
	  /* Found goal with deletion */
	  debug7(printf(" possible deletion, continuing\n"));
	  /* continue */

	} else {
	  debug7(printf(" %u advanced past goal %u, and previous %u before goal %u - %d, so not spanning\n",
			*elt->positions,goal,elt->positions[-1],goal,overall_max_distance));
	  return false;
	}
      }
    }
  }

  debug7(printf("Returning true\n"));
  return true;
}


/* Taken from stage1hr.c identify_multimiss_iter */
static bool
extend_leftward (Univcoord_T goal, Univcoord_T chroffset, Univcoord_T chrhigh,
		 List_T set, T sarray, char *queryptr, Compress_T query_compress,
		 bool plusp, int genestrand, int best_querystart, int best_queryend) {
  Elt_T elt;
  UINT4 nmatches;
  Sarrayptr_T initptr, finalptr;
  bool successp;
  UINT4 queryend, querypos;
  Univcoord_T low, high;


  debug7(printf("extend_leftward, plusp %d, with goal %u, querystart..queryend %d..%d\n",
		plusp,goal,best_querystart,best_queryend));

  queryend = best_querystart - 2;

  for ( ; set /* != NULL */; set = set->rest) {
    debug7(Elt_dump_list(set));
    elt = (Elt_T) set->first;
    debug7(printf("remaining elts %d: ",List_length(set)));
    debug7(printf("%d..%d\n",elt->querystart,elt->queryend));
    
    debug7(printf("Checking for re-compute of left region: elt->queryend %d vs queryend %d\n",elt->queryend,queryend));
    if (/* elt->queryend != queryend && */ elt->queryend > queryend) {
      debug7(printf("Re-computing left region\n"));
      querypos = elt->querystart;

      sarray_search(&initptr,&finalptr,&successp,&nmatches,&(queryptr[querypos]),
		    /*querylength*/(queryend + 1) - querypos,/*queryoffset*/querypos,
		    query_compress,plusp);
      Elt_replace(elt,querypos,nmatches,initptr,finalptr);
      /* set->first = (void *) elt; */
    }
    queryend = elt->querystart - 2;

    debug7(printf("remaining elts %d: ",List_length(set)));
    debug7(printf("%d..%d\n",elt->querystart,elt->queryend));
    if (elt->querystart > best_queryend) {
      /* Allow for deletion with higher goal */
      debug7(printf("Allow for deletion with higher goal: %d positions\n",elt->npositions));
      low = subtract_bounded(goal,/*minusterm*/max_insertionlen,chroffset);
      high = add_bounded(goal,/*plusterm*/overall_max_distance,chrhigh);
      Elt_fill_positions_filtered(elt,sarray,goal,low,high,query_compress,plusp,genestrand);
      if (elt->npositions <= 0) {
	/* List is empty, so one more miss seen. */
	debug7(printf(" positions empty, so not spanning\n"));
	return false;
	
      } else if (*elt->positions > high) {
	/* Already advanced past goal, so one more miss seen. */
	debug7(printf(" %u advanced past goal %u + %d, so not spanning\n",*elt->positions,goal,overall_max_distance));
	return false;

      } else {
	/* Found goal.  Advance past goal and continue with loop. */
	debug7(printf(" advancing\n"));
	if ((nmatches = Genome_consecutive_matches_leftward(query_compress,/*left*/*elt->positions,
							    /*pos5*/0,/*pos3*/elt->querystart,
							    plusp,genestrand)) > 0) {
	  debug7(printf(" extending querystart %d leftward by %d matches\n",elt->querystart,nmatches));
	  elt->querystart -= nmatches;
	  queryend = elt->querystart - 2;
	}
	++elt->positions;
	--elt->npositions;
	/* continue */
      }
    } else {
      /* Allow for deletion with lower goal */
      debug7(printf("Allow for deletion with lower goal: %d positions\n",elt->npositions));
      low = subtract_bounded(goal,/*minusterm*/overall_max_distance,chroffset);
      high = add_bounded(goal,/*plusterm*/max_insertionlen,chrhigh);
      Elt_fill_positions_filtered(elt,sarray,goal,low,high,query_compress,plusp,genestrand);
      if (elt->npositions <= 0 && elt->positions == elt->positions_allocated) {
	/* List is empty, and no previous one exists */
	debug7(printf(" list is empty and no previous, so not spanning\n"));
	return false;

      } else if (elt->npositions <= 0) {
	/* List is empty, but previous one exists */
	if (elt->positions[-1] >= low) {
	  /* Found goal with deletion */
	  debug7(printf(" possible deletion, continuing\n"));
	  if ((nmatches = Genome_consecutive_matches_leftward(query_compress,/*left*/elt->positions[-1],
							      /*pos5*/0,/*pos3*/elt->querystart,
							      plusp,genestrand)) > 0) {
	    debug7(printf(" extending querystart %d leftward by %d matches\n",elt->querystart,nmatches));
	    elt->querystart -= nmatches;
	    queryend = elt->querystart - 2;
	  }
	  /* continue */
	} else {
	  debug7(printf(" previous %u before goal %u - %d, so not spanning\n",elt->positions[-1],goal,overall_max_distance));
	  return false;
	}
	
      } else if (elt->positions == elt->positions_allocated) {
	/* List is at beginning, but current one exists */
	if (*elt->positions > goal) {
	  /* Already advanced past goal, so one more miss seen. */
	  debug7(printf(" %u advanced past goal_high %u, so not spanning\n",*elt->positions,goal));
	  return false;
	} else {
	  /* Found goal.  Advance past goal and continue with loop. */
	  debug7(printf(" advancing\n"));
	  if ((nmatches = Genome_consecutive_matches_leftward(query_compress,/*left*/*elt->positions,
							      /*pos5*/0,/*pos3*/elt->querystart,
							      plusp,genestrand)) > 0) {
	    debug7(printf(" extending querystart %d leftward by %d matches\n",elt->querystart,nmatches));
	    elt->querystart -= nmatches;
	    queryend = elt->querystart - 2;
	  }
	  ++elt->positions;
	  --elt->npositions;
	  /* continue */
	}

      } else {
	/* Test both current one (for goal) and previous one (for deletion) */
	if (*elt->positions == goal) {
	  /* Found goal.  Advance past goal and continue with loop. */
	  debug7(printf(" advancing\n"));
	  if ((nmatches = Genome_consecutive_matches_leftward(query_compress,/*left*/*elt->positions,
							      /*pos5*/0,/*pos3*/elt->querystart,
							      plusp,genestrand)) > 0) {
	    debug7(printf(" extending querystart %d leftward by %d matches\n",elt->querystart,nmatches));
	    elt->querystart -= nmatches;
	    queryend = elt->querystart - 2;
	  }
	  ++elt->positions;
	  --elt->npositions;
	  /* continue */

	} else if (elt->positions[-1] >= low) {
	  /* Found goal with deletion */
	  debug7(printf(" possible deletion, continuing\n"));
	  if ((nmatches = Genome_consecutive_matches_leftward(query_compress,/*left*/elt->positions[-1],
							      /*pos5*/0,/*pos3*/elt->querystart,
							      plusp,genestrand)) > 0) {
	    debug7(printf(" extending querystart %d leftward by %d matches\n",elt->querystart,nmatches));
	    elt->querystart -= nmatches;
	    queryend = elt->querystart - 2;
	  }
	  /* continue */

	} else {
	  debug7(printf(" %u advanced past goal %u, and previous %u before goal %u - %d, so not spanning\n",
			*elt->positions,goal,elt->positions[-1],goal,overall_max_distance));
	  return false;
	}
      }
    }
  }

  debug7(printf("Returning true\n"));
  return true;
}



static void
collect_elt_matches (int *found_score, List_T *subs, List_T *indels, List_T *singlesplicing,
		     List_T *doublesplicing, int querystart_same, int queryend_same,
		     Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
		     Chrpos_T chrlength, Univcoord_T goal, 
		     List_T rightward_set, List_T leftward_set, int querylength, Compress_T query_compress,
		     bool plusp, int genestrand, int nmisses_allowed, bool first_read_p) {
  List_T set, p;
  Stage3end_T hit;
  Elt_T elt;
  Univcoord_T left, left1, left2, *array;
  Uintlist_T difflist = NULL;	/* Won't work with LARGE_GENOMES */
  int nmismatches, nindels;
  int nsame, ndiff;
  int querystart_diff, queryend_diff, indel_pos;
#if 0
  int nmismatches1, nmismatches2;
  int query_indel_pos;
#endif

  List_T lowprob;
  int nhits;
  int segmenti_donor_knownpos[MAX_READLENGTH+1], segmentj_acceptor_knownpos[MAX_READLENGTH+1],
    segmentj_antidonor_knownpos[MAX_READLENGTH+1], segmenti_antiacceptor_knownpos[MAX_READLENGTH+1];
  int segmenti_donor_knowni[MAX_READLENGTH+1], segmentj_acceptor_knowni[MAX_READLENGTH+1],
    segmentj_antidonor_knowni[MAX_READLENGTH+1], segmenti_antiacceptor_knowni[MAX_READLENGTH+1];
  int segmenti_donor_nknown, segmentj_acceptor_nknown,
    segmentj_antidonor_nknown, segmenti_antiacceptor_nknown;
  int j, i, n;
  bool segmenti_usedp, segmentj_usedp;
  bool foundp;

  /* Potential success */
  debug7(printf("  successful candidate found\n"));
  if (goal < (Univcoord_T) querylength) {
    debug7(printf("  Goes over beginning of chromosome\n"));
    return;
  } else if (goal + querylength > chrhigh) {
    debug7(printf("  Goes over end of chromosome\n"));
    return;
  } else {
    left = goal /* - querylength */;
  }

  nsame = ndiff = 0;
  querystart_diff = querylength;
  queryend_diff = 0;
  for (set = rightward_set; set /* != NULL */; set = set->rest) {
    elt = (Elt_T) set->first;
    debug7(printf("%d..%d:%u vs %u: ",elt->querystart,elt->queryend,elt->positions[-1],goal));
    assert(elt->filledp == true);
    if (elt->positions[-1] == goal) {
      debug7(printf("same\n"));
      if (elt->querystart < querystart_same) {
	querystart_same = elt->querystart;
      }
      if (elt->queryend > queryend_same) {
	queryend_same = elt->queryend;
      }
      nsame++;

    } else {
      debug7(printf("diff (npositions %d)\n",elt->npositions));
      debug7(printf("Pushing position %u\n",elt->positions[-1]));
      difflist = Uintlist_push(difflist,elt->positions[-1]);
      for (i = 0; i < elt->npositions; i++) {
	debug7(printf("Pushing position %u\n",elt->positions[i]));
	difflist = Uintlist_push(difflist,elt->positions[i]);
      }
      if (elt->querystart < querystart_diff) {
	querystart_diff = elt->querystart;
      }
      if (elt->queryend > queryend_diff) {
	queryend_diff = elt->queryend;
      }
      ndiff++;
    }
  }

  for (set = leftward_set; set /* != NULL */; set = set->rest) {
    elt = (Elt_T) set->first;
    debug7(printf("%d..%d:%u vs %u: ",elt->querystart,elt->queryend,elt->positions[-1],goal));
    assert(elt->filledp == true);
    if (elt->positions[-1] == goal) {
      debug7(printf("same\n"));
      if (elt->querystart < querystart_same) {
	querystart_same = elt->querystart;
      }
      if (elt->queryend > queryend_same) {
	queryend_same = elt->queryend;
      }
      nsame++;

    } else {
      debug7(printf("diff (npositions %d)\n",elt->npositions));
      debug7(printf("Pushing position %u\n",elt->positions[-1]));
      difflist = Uintlist_push(difflist,elt->positions[-1]);
      for (i = 0; i < elt->npositions; i++) {
	debug7(printf("Pushing position %u\n",elt->positions[i]));
	difflist = Uintlist_push(difflist,elt->positions[i]);
      }
      if (elt->querystart < querystart_diff) {
	querystart_diff = elt->querystart;
      }
      if (elt->queryend > queryend_diff) {
	queryend_diff = elt->queryend;
      }
      ndiff++;
    }
  }

  debug7(printf("Got %d same, %d diff\n",nsame,ndiff));

  if (ndiff == 0) {
    /* sub */
    debug7(printf("  Testing in entire query\n"));
    nmismatches = Genome_count_mismatches_substring(query_compress,left,/*pos5*/0,/*pos3*/querylength,
						    plusp,genestrand);
    debug7(printf("nmismatches = %d (vs %d misses allowed)\n",nmismatches,nmisses_allowed));

    if (nmismatches > nmisses_allowed) {
      debug7(printf("Result: too many mismatches\n"));

    } else {
      debug7(printf("Result: successful hit saved\n"));
      debug(printf("Reporting hit with %d mismatches\n",nmismatches));
      if ((hit = Stage3end_new_substitution(&(*found_score),nmismatches,
					    left,/*genomiclength*/querylength,
					    query_compress,plusp,genestrand,
					    chrnum,chroffset,chrhigh,chrlength,
					    /*sarrayp*/true)) != NULL) {
	*subs = List_push(*subs,(void *) hit);
      }
    }
    assert(difflist == NULL);

  } else if (querystart_same == 0 && queryend_diff == querylength - 1) {
    left1 = left;
    indel_pos = queryend_same + 1;
    debug7(printf("same is at %u from %d to %d\n",left,querystart_same,queryend_same));

    array = Uintlist_to_array(&n,difflist);
    qsort(array,n,sizeof(Univcoord_T),Univcoord_compare);
    Uintlist_free(&difflist);

    for (i = 0; i < n; i++) {
      left2 = array[i];
      debug7(printf("diff is at %u, from %d to %d\n",left2,querystart_diff - 1,queryend_diff));

      if (i > 0 && left2 == array[i-1]) {
	/* Already processed */

      } else if (left2 + querylength >= chrhigh) {
	/* Splice or deletion would extend to next chromosome */

      } else if (left2 > left1 + max_deletionlen) {
	debug7(printf("A splice..."));

	segmenti_donor_nknown = segmenti_antiacceptor_nknown = 0;
	if (nsplicesites > 0 &&
	    Splicetrie_splicesite_p(left1,/*pos5*/1,/*pos3*/querylength) == true) {
	  j = binary_search(0,nsplicesites,splicesites,left1);
	  while (j < nsplicesites && splicesites[j] < left1 + querylength) {
	    if (splicetypes[j] == DONOR) {
	      debug4s(printf("Setting known donor %d for segmenti at %u\n",j,splicesites[j]));
	      segmenti_donor_knownpos[segmenti_donor_nknown] = splicesites[j] - left1;
	      segmenti_donor_knowni[segmenti_donor_nknown++] = j;
	    } else if (splicetypes[j] == ANTIACCEPTOR) {
	      debug4s(printf("Setting known antiacceptor %d for segmenti at %u\n",j,splicesites[j]));
	      segmenti_antiacceptor_knownpos[segmenti_antiacceptor_nknown] = splicesites[j] - left1;
	      segmenti_antiacceptor_knowni[segmenti_antiacceptor_nknown++] = j;
	    }
	    j++;
	  }
	}
	segmenti_donor_knownpos[segmenti_donor_nknown] = MAX_READLENGTH;
	segmenti_antiacceptor_knownpos[segmenti_antiacceptor_nknown] = MAX_READLENGTH;
	  
	segmentj_acceptor_nknown = segmentj_antidonor_nknown = 0;
	if (nsplicesites > 0 &&
	    Splicetrie_splicesite_p(left2,/*pos5*/1,/*pos3*/querylength) == true) {
	  j = binary_search(0,nsplicesites,splicesites,left2);
	  while (j < nsplicesites && splicesites[j] < left2 + querylength) {
	    if (splicetypes[j] == ACCEPTOR) {
	      debug4s(printf("Setting known acceptor %d for segmentj at %u\n",j,splicesites[j]));
	      segmentj_acceptor_knownpos[segmentj_acceptor_nknown] = splicesites[j] - left2;
	      segmentj_acceptor_knowni[segmentj_acceptor_nknown++] = j;
	    } else if (splicetypes[j] == ANTIDONOR) {
	      debug4s(printf("Setting known antidonor %d for segmentj at %u\n",j,splicesites[j]));
	      segmentj_antidonor_knownpos[segmentj_antidonor_nknown] = splicesites[j] - left2;
	      segmentj_antidonor_knowni[segmentj_antidonor_nknown++] = j;
	    }
	    j++;
	  }
	}
	segmentj_acceptor_knownpos[segmentj_acceptor_nknown] = MAX_READLENGTH;
	segmentj_antidonor_knownpos[segmentj_antidonor_nknown] = MAX_READLENGTH;

	lowprob = (List_T) NULL;
	*singlesplicing = Splice_solve_single(&(*found_score),&nhits,*singlesplicing,&lowprob,
					      &segmenti_usedp,&segmentj_usedp,
					      /*segmenti_left*/left1,/*segmentj_left*/left2,
					      chrnum,chroffset,chrhigh,chrlength,
					      chrnum,chroffset,chrhigh,chrlength,
					      querylength,query_compress,
					      segmenti_donor_knownpos,segmentj_acceptor_knownpos,
					      segmentj_antidonor_knownpos,segmenti_antiacceptor_knownpos,
					      segmenti_donor_knowni,segmentj_acceptor_knowni,
					      segmentj_antidonor_knowni,segmenti_antiacceptor_knowni,
					      segmenti_donor_nknown,segmentj_acceptor_nknown,
					      segmentj_antidonor_nknown,segmenti_antiacceptor_nknown,
					      splicing_penalty,/*max_mismatches_allowed*/1000,
					      first_read_p,plusp,genestrand,/*subs_or_indels_p*/false,
					      /*sarrayp*/true);
	for (p = lowprob; p != NULL; p = List_next(p)) {
	  debug7(printf("freeing lowprob..."));
	  hit = (Stage3end_T) List_head(p);
	  Stage3end_free(&hit);
	}
	List_free(&lowprob);
	debug7(printf("\n"));

      } else if (left2 > left1) {
	nindels = left2 - left1;
	debug7(printf("B deletion of %d bp relative to max_deletionlen %d (nmisses allowed %d)...",
		      nindels,max_deletionlen,nmisses_allowed));
	if ((indel_pos < 17 || querylength - indel_pos < 17) && nindels > max_end_deletions) {
	  /* Allow regular GSNAP algorithm to find this */
	  debug7(printf("too long for end deletion"));
	} else {
#if 0
	  nmismatches1 = Genome_count_mismatches_substring(query_compress,left1,/*pos5*/0,/*pos3*/indel_pos,
							   plusp,genestrand);
	  nmismatches2 = Genome_count_mismatches_substring(query_compress,left2,/*pos5*/indel_pos,
							   /*pos3*/querylength,plusp,genestrand);
	  if (plusp == true) {
	    query_indel_pos = indel_pos;
	  } else {
	    query_indel_pos = querylength - indel_pos;
	  }
	  if ((hit = Stage3end_new_deletion(&(*found_score),nindels,query_indel_pos,
					    nmismatches1,nmismatches2,
					    left1,/*genomiclength*/querylength+nindels,
					    query_compress,querylength,plusp,genestrand,
					    chrnum,chroffset,chrhigh,chrlength,
					    /*indel_penalty*/2,/*sarrayp*/true)) != NULL) {
	    debug7(printf("successful"));
	    *indels = List_push(*indels,(void *) hit);
	  }
#else
	  *indels = Indel_solve_middle_deletion(&foundp,&(*found_score),&nhits,*indels,
						/*left*/left1,chrnum,chroffset,chrhigh,chrlength,
						/*indels*/-nindels,query_compress,querylength,nmisses_allowed,
						plusp,genestrand,/*sarray*/true);
	  debug7(
		 if (foundp == true) {
		   printf("successful");
		 }
		 );
#endif
	}
	debug7(printf("\n"));
      
      } else if (left2 < left1) {
	nindels = left1 - left2;
	if (nindels >= indel_pos || indel_pos + nindels >= querylength) {
	  debug7(printf("X insertion of %d bp too long\n",nindels));
	} else {
	  debug7(printf("C insertion of %d bp (nmisses allowed %d)...",nindels,nmisses_allowed));
#if 0
	  nmismatches1 = Genome_count_mismatches_substring(query_compress,left1,/*pos5*/0,/*pos3*/indel_pos-nindels,
							   plusp,genestrand);
	  nmismatches2 = Genome_count_mismatches_substring(query_compress,left2,/*pos5*/indel_pos+nindels,
							   /*pos3*/querylength,plusp,genestrand);
	  if (plusp == true) {
	    query_indel_pos = indel_pos;
	  } else {
	    query_indel_pos = querylength - indel_pos - nindels;
	  }
	  if ((hit = Stage3end_new_insertion(&(*found_score),nindels,query_indel_pos,
					     nmismatches1,nmismatches2,
					     left1,/*genomiclength*/querylength-nindels,
					     query_compress,querylength,plusp,genestrand,
					     chrnum,chroffset,chrhigh,chrlength,
					     /*indel_penalty*/2,/*sarrayp*/true)) != NULL) {
	    debug7(printf("successful"));
	    *indels = List_push(*indels,(void *) hit);
	  }
#else
	  *indels = Indel_solve_middle_insertion(&foundp,&(*found_score),&nhits,*indels,
						 /*left*/left1,chrnum,chroffset,chrhigh,chrlength,
						 /*indels*/+nindels,query_compress,querylength,nmisses_allowed,
						 plusp,genestrand,/*sarrayp*/true);
	  debug7(
		 if (foundp == true) {
		   printf("successful");
		 }
		 );
#endif
	  debug7(printf("\n"));
	}
      }
    }

    FREE(array);

  } else if (querystart_diff == 0 && queryend_same == querylength - 1) {
    left2 = left;
    indel_pos = querystart_same;
    debug7(printf("same is at %u from %d to %d\n",left,querystart_same,queryend_same));
    
    array = Uintlist_to_array(&n,difflist);
    qsort(array,n,sizeof(Univcoord_T),Univcoord_compare);
    Uintlist_free(&difflist);

    for (i = 0; i < n; i++) {
      left1 = array[i];
      debug7(printf("diff is at %u, from %d to %d\n",left1,querystart_diff,queryend_diff));

      if (i > 0 && left1 == array[i-1]) {
	/* Already processed */

      } else if (left2 + querylength >= chrhigh) {
	/* Splice or deletion would extend to next chromosome */

      } else if (left2 > left1 + max_deletionlen) {
	debug7(printf("A splice..."));

	segmenti_donor_nknown = segmenti_antiacceptor_nknown = 0;
	if (nsplicesites > 0 &&
	    Splicetrie_splicesite_p(left1,/*pos5*/1,/*pos3*/querylength) == true) {
	  j = binary_search(0,nsplicesites,splicesites,left1);
	  while (j < nsplicesites && splicesites[j] < left1 + querylength) {
	    if (splicetypes[j] == DONOR) {
	      debug4s(printf("Setting known donor %d for segmenti at %u\n",j,splicesites[j]));
	      segmenti_donor_knownpos[segmenti_donor_nknown] = splicesites[j] - left1;
	      segmenti_donor_knowni[segmenti_donor_nknown++] = j;
	    } else if (splicetypes[j] == ANTIACCEPTOR) {
	      debug4s(printf("Setting known antiacceptor %d for segmenti at %u\n",j,splicesites[j]));
	      segmenti_antiacceptor_knownpos[segmenti_antiacceptor_nknown] = splicesites[j] - left1;
	      segmenti_antiacceptor_knowni[segmenti_antiacceptor_nknown++] = j;
	    }
	    j++;
	  }
	}
	segmenti_donor_knownpos[segmenti_donor_nknown] = MAX_READLENGTH;
	segmenti_antiacceptor_knownpos[segmenti_antiacceptor_nknown] = MAX_READLENGTH;
	  
	segmentj_acceptor_nknown = segmentj_antidonor_nknown = 0;
	if (nsplicesites > 0 &&
	    Splicetrie_splicesite_p(left2,/*pos5*/1,/*pos3*/querylength) == true) {
	  j = binary_search(0,nsplicesites,splicesites,left2);
	  while (j < nsplicesites && splicesites[j] < left2 + querylength) {
	    if (splicetypes[j] == ACCEPTOR) {
	      debug4s(printf("Setting known acceptor %d for segmentj at %u\n",j,splicesites[j]));
	      segmentj_acceptor_knownpos[segmentj_acceptor_nknown] = splicesites[j] - left2;
	      segmentj_acceptor_knowni[segmentj_acceptor_nknown++] = j;
	    } else if (splicetypes[j] == ANTIDONOR) {
	      debug4s(printf("Setting known antidonor %d for segmentj at %u\n",j,splicesites[j]));
	      segmentj_antidonor_knownpos[segmentj_antidonor_nknown] = splicesites[j] - left2;
	      segmentj_antidonor_knowni[segmentj_antidonor_nknown++] = j;
	    }
	    j++;
	  }
	}
	segmentj_acceptor_knownpos[segmentj_acceptor_nknown] = MAX_READLENGTH;
	segmentj_antidonor_knownpos[segmentj_antidonor_nknown] = MAX_READLENGTH;

	lowprob = (List_T) NULL;
	*singlesplicing = Splice_solve_single(&(*found_score),&nhits,*singlesplicing,&lowprob,
					      &segmenti_usedp,&segmentj_usedp,
					      /*segmenti_left*/left1,/*segmentj_left*/left2,
					      chrnum,chroffset,chrhigh,chrlength,
					      chrnum,chroffset,chrhigh,chrlength,
					      querylength,query_compress,
					      segmenti_donor_knownpos,segmentj_acceptor_knownpos,
					      segmentj_antidonor_knownpos,segmenti_antiacceptor_knownpos,
					      segmenti_donor_knowni,segmentj_acceptor_knowni,
					      segmentj_antidonor_knowni,segmenti_antiacceptor_knowni,
					      segmenti_donor_nknown,segmentj_acceptor_nknown,
					      segmentj_antidonor_nknown,segmenti_antiacceptor_nknown,
					      splicing_penalty,/*max_mismatches_allowed*/1000,
					      first_read_p,plusp,genestrand,/*subs_or_indels_p*/false,
					      /*sarrayp*/true);
	for (p = lowprob; p != NULL; p = List_next(p)) {
	  debug7(printf("freeing lowprob..."));
	  hit = (Stage3end_T) List_head(p);
	  Stage3end_free(&hit);
	}
	List_free(&lowprob);
	debug7(printf("\n"));

      } else if (left2 > left1) {
	nindels = left2 - left1;
	debug7(printf("B deletion of %d bp relative to max_deletionlen %d (nmisses allowed %d)...",
		      nindels,max_deletionlen,nmisses_allowed));
	if ((indel_pos < 17 || querylength - indel_pos < 17) && nindels > max_end_deletions) {
	  /* Allow regular GSNAP algorithm to find this */
	  debug7(printf("too long for end deletion"));
	} else {
#if 0
	  nmismatches1 = Genome_count_mismatches_substring(query_compress,left1,/*pos5*/0,/*pos3*/indel_pos,
							   plusp,genestrand);
	  nmismatches2 = Genome_count_mismatches_substring(query_compress,left2,/*pos5*/indel_pos,
							   /*pos3*/querylength,plusp,genestrand);
	  if (plusp == true) {
	    query_indel_pos = indel_pos;
	  } else {
	    query_indel_pos = querylength - indel_pos;
	  }
	  if ((hit = Stage3end_new_deletion(&(*found_score),nindels,query_indel_pos,
					    nmismatches1,nmismatches2,
					    left1,/*genomiclength*/querylength+nindels,
					    query_compress,querylength,plusp,genestrand,
					    chrnum,chroffset,chrhigh,chrlength,
					    /*indel_penalty*/2,/*sarrayp*/true)) != NULL) {
	    debug7(printf("successful"));
	    *indels = List_push(*indels,(void *) hit);
	  }
#else
	  *indels = Indel_solve_middle_deletion(&foundp,&(*found_score),&nhits,*indels,
						/*left*/left1,chrnum,chroffset,chrhigh,chrlength,
						/*indels*/-nindels,query_compress,querylength,nmisses_allowed,
						plusp,genestrand,/*sarray*/true);
	  debug7(
		 if (foundp == true) {
		   printf("successful");
		 }
		 );
#endif
	}
	debug7(printf("\n"));
      
      } else if (left2 < left1) {
	nindels = left1 - left2;
	if (nindels >= indel_pos || indel_pos + nindels >= querylength) {
	  debug7(printf("X insertion of %d bp too long\n",nindels));
	} else {
	  debug7(printf("C insertion of %d bp (nmisses allowed %d)...",nindels,nmisses_allowed));
#if 0      
	  nmismatches1 = Genome_count_mismatches_substring(query_compress,left1,/*pos5*/0,/*pos3*/indel_pos-nindels,
							   plusp,genestrand);
	  nmismatches2 = Genome_count_mismatches_substring(query_compress,left2,/*pos5*/indel_pos+nindels,
							   /*pos3*/querylength,plusp,genestrand);
	  if (plusp == true) {
	    query_indel_pos = indel_pos;
	  } else {
	    query_indel_pos = querylength - indel_pos - nindels;
	  }
	  if ((hit = Stage3end_new_insertion(&(*found_score),nindels,query_indel_pos,
					     nmismatches1,nmismatches2,
					     left1,/*genomiclength*/querylength-nindels,
					     query_compress,querylength,plusp,genestrand,
					     chrnum,chroffset,chrhigh,chrlength,
					     /*indel_penalty*/2,/*sarrayp*/true)) != NULL) {
	    debug7(printf("successful"));
	    *indels = List_push(*indels,(void *) hit);
	  }
#else
	  *indels = Indel_solve_middle_insertion(&foundp,&(*found_score),&nhits,*indels,
						 /*left*/left1,chrnum,chroffset,chrhigh,chrlength,
						 /*indels*/+nindels,query_compress,querylength,nmisses_allowed,
						 plusp,genestrand,/*sarrayp*/true);
	  debug7(
		 if (foundp == true) {
		   printf("successful");
		 }
		 );
#endif
	  debug7(printf("\n"));
	}
      }
    }

    FREE(array);

  } else {
    Uintlist_free(&difflist);
  }

  return;
}


void
Sarray_search_greedy (int *found_score, List_T *subs, List_T *indels, List_T *singlesplicing,
		      List_T *doublesplicing, char *queryuc_ptr, char *queryrc, int querylength,
		      Compress_T query_compress_fwd, Compress_T query_compress_rev,
		      int nmisses_allowed, bool first_read_p) {
  List_T plus_set, minus_set, p;
  List_T rightward_set, leftward_set;
  Elt_T best_plus_elt, best_minus_elt, elt, *array;
  UINT4 best_plus_nmatches, best_minus_nmatches, nmatches;
  Sarrayptr_T initptr, finalptr;
  int plus_niter, minus_niter;
  bool successp;
  int plus_querypos, minus_querypos, halfwaypos;
  int nelts, i;
  Chrnum_T chrnum;
  Univcoord_T chroffset, chrhigh, left;
  Chrpos_T chrlength;
  Stage3end_T hit;
  int nmismatches;


  if (nmisses_allowed < 0) {
    nmisses_allowed = 0;
  }
  debug(printf("\nStarting Sarray_search_greedy with querylength %d and indexsize %d and nmisses_allowed %d\n",
	       querylength,sarray->indexsize,nmisses_allowed));

  *found_score = querylength;

  /* Do one plus round */
  plus_querypos = 0;
  sarray_search(&initptr,&finalptr,&successp,&best_plus_nmatches,&(queryuc_ptr[plus_querypos]),
		querylength - plus_querypos,/*queryoffset*/plus_querypos,
		query_compress_fwd,/*plusp*/true);
  best_plus_elt = Elt_new(plus_querypos,best_plus_nmatches,initptr,finalptr);
  plus_querypos += (int) best_plus_nmatches;
  plus_querypos += 1;		/* To skip the presumed mismatch */


  /* Do one minus round */
  minus_querypos = 0;
  sarray_search(&initptr,&finalptr,&successp,&best_minus_nmatches,&(queryrc[minus_querypos]),
		querylength - minus_querypos,/*queryoffset*/minus_querypos,
		query_compress_rev,/*plusp*/false);
  best_minus_elt = Elt_new(minus_querypos,best_minus_nmatches,initptr,finalptr);
  minus_querypos += (int) best_minus_nmatches;
  minus_querypos += 1;		/* To skip the presumed mismatch */


  if (best_plus_nmatches >= querylength/2) {
    /* See if we have a winner */
    debug(printf("best_plus_nmatches = %d > %d/2, so checking mismatches against %d allowed\n",
		 best_plus_nmatches,querylength,nmisses_allowed));
    Elt_fill_positions_all(best_plus_elt,sarray);
    for (i = 0; i < best_plus_elt->npositions; i++) {
      left = best_plus_elt->positions[i];
      /* Should return max_mismatches + 1 if it exceeds the limit */
      if ((nmismatches = Genome_count_mismatches_limit(query_compress_fwd,left,/*pos5*/0,/*pos3*/querylength,
						       /*max_mismatches*/nmisses_allowed,/*plusp*/true,/*genestrand*/0)) <= nmisses_allowed) {
	chrnum = Univ_IIT_get_one(chromosome_iit,left,left);
	Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
	debug(printf("Case 1: New substitution from beginning\n"));
	if ((hit = Stage3end_new_substitution(&(*found_score),nmismatches,
					      left,/*genomiclength*/querylength,
					      query_compress_fwd,/*plusp*/true,/*genestrand*/0,
					      chrnum,chroffset,chrhigh,chrlength,
					      /*sarrayp*/true)) != NULL) {
	  *subs = List_push(*subs,(void *) hit);
	}
      }
      debug(printf("Looking at plus position %u => %d mismatches\n",left,nmismatches));
    }

  } else {
    /* Try starting from middle of read */
    halfwaypos = querylength/2;
    debug(printf("Starting from halfway point on plus\n"));
    sarray_search(&initptr,&finalptr,&successp,&nmatches,&(queryuc_ptr[halfwaypos]),
		  querylength - halfwaypos,/*queryoffset*/halfwaypos,
		  query_compress_fwd,/*plusp*/true);
    if (nmatches >= querylength - halfwaypos) {
      elt = Elt_new(halfwaypos,nmatches,initptr,finalptr);
      Elt_fill_positions_all(elt,sarray);
      for (i = 0; i < elt->npositions; i++) {
	left = elt->positions[i];
	/* Should return max_mismatches + 1 if it exceeds the limit */
	if ((nmismatches = Genome_count_mismatches_limit(query_compress_fwd,left,/*pos5*/0,/*pos3*/querylength,
							 /*max_mismatches*/nmisses_allowed,/*plusp*/true,/*genestrand*/0)) <= nmisses_allowed) {
	  chrnum = Univ_IIT_get_one(chromosome_iit,left,left);
	  Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
	  debug(printf("Case 1: New substitution from middle\n"));
	  if ((hit = Stage3end_new_substitution(&(*found_score),nmismatches,
						left,/*genomiclength*/querylength,
						query_compress_fwd,/*plusp*/true,/*genestrand*/0,
						chrnum,chroffset,chrhigh,chrlength,
						/*sarrayp*/true)) != NULL) {
	    *subs = List_push(*subs,(void *) hit);
	  }
	}
	debug(printf("Looking at plus position %u => %d mismatches\n",left,nmismatches));
      }
      Elt_free(&elt);
    }
  }


  if (best_minus_nmatches >= querylength/2) {
    /* See if we have a winner */
    debug(printf("best_minus_nmatches = %d > %d/2, so checking mismatches against %d allowed\n",
		 best_minus_nmatches,querylength,nmisses_allowed));
    Elt_fill_positions_all(best_minus_elt,sarray);
    for (i = 0; i < best_minus_elt->npositions; i++) {
      left = best_minus_elt->positions[i];
      /* Should return max_mismatches + 1 if it exceeds the limit */
      if ((nmismatches = Genome_count_mismatches_limit(query_compress_rev,left,/*pos5*/0,/*pos3*/querylength,
						       /*max_mismatches*/nmisses_allowed,/*plusp*/false,/*genestrand*/0)) <= nmisses_allowed) {
	chrnum = Univ_IIT_get_one(chromosome_iit,left,left);
	Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
	debug(printf("Case 2: New substitution from beginning\n"));
	if ((hit = Stage3end_new_substitution(&(*found_score),nmismatches,
					      left,/*genomiclength*/querylength,
					      query_compress_rev,/*plusp*/false,/*genestrand*/0,
					      chrnum,chroffset,chrhigh,chrlength,
					      /*sarrayp*/true)) != NULL) {
	  *subs = List_push(*subs,(void *) hit);
	}
      }
      debug(printf("Looking at minus position %u => %d mismatches\n",left,nmismatches));
    }
  } else {
    /* Try starting from middle of read */
    halfwaypos = querylength/2;
    debug(printf("Starting from halfway point on minus\n"));
    sarray_search(&initptr,&finalptr,&successp,&nmatches,&(queryrc[halfwaypos]),
		  querylength - halfwaypos,/*queryoffset*/halfwaypos,
		  query_compress_rev,/*plusp*/false);
    if (nmatches >= querylength - halfwaypos) {
      elt = Elt_new(halfwaypos,nmatches,initptr,finalptr);
      Elt_fill_positions_all(elt,sarray);
      for (i = 0; i < elt->npositions; i++) {
	left = elt->positions[i];
	/* Should return max_mismatches + 1 if it exceeds the limit */
	if ((nmismatches = Genome_count_mismatches_limit(query_compress_rev,left,/*pos5*/0,/*pos3*/querylength,
							 /*max_mismatches*/nmisses_allowed,/*plusp*/false,/*genestrand*/0)) <= nmisses_allowed) {
	  chrnum = Univ_IIT_get_one(chromosome_iit,left,left);
	  Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
	  debug(printf("Case 2: New substitution from middle\n"));
	  if ((hit = Stage3end_new_substitution(&(*found_score),nmismatches,
						left,/*genomiclength*/querylength,
						query_compress_rev,/*plusp*/false,/*genestrand*/0,
						chrnum,chroffset,chrhigh,chrlength,
						/*sarrayp*/true)) != NULL) {
	    *subs = List_push(*subs,(void *) hit);
	  }
	}
	debug(printf("Looking at minus position %u => %d mismatches\n",left,nmismatches));
      }
      Elt_free(&elt);
    }
  }

    
  debug(printf("Found %d subs\n",List_length(*subs)));
#if 0
  /* Allow identification of splicing, even if substitutions are found */
  if (*subs != NULL) {
    /* Be satisfied with 1-mismatch results */
    Elt_free(&best_plus_elt);
    Elt_free(&best_minus_elt);
    return;
  }
#endif

  if (plus_querypos >= querylength) {
    plus_set = (List_T) NULL;
  } else {
    /* Extend plus side a second time */
    sarray_search(&initptr,&finalptr,&successp,&nmatches,&(queryuc_ptr[plus_querypos]),
		  querylength - plus_querypos,/*queryoffset*/plus_querypos,
		  query_compress_fwd,/*plusp*/true);
    elt = Elt_new(plus_querypos,nmatches,initptr,finalptr);
    plus_querypos += nmatches;
    plus_querypos += 1;		/* To skip the presumed mismatch */

    debug(printf("plus_querypos %d vs querylength %d\n",plus_querypos,querylength));
    if (nmatches <= best_plus_nmatches) {
      /* Initial (left) elt was best */
      debug(printf("Initial elt was best\n"));
      plus_set = List_push(NULL,elt);
      if (plus_querypos >= querylength) {
	chrhigh = 0U;
	Elt_fill_positions_all(best_plus_elt,sarray);
	for (i = 0; i < best_plus_elt->npositions; i++) {
	  left = best_plus_elt->positions[i];
	  if (left > chrhigh) {
	    chrnum = Univ_IIT_get_one(chromosome_iit,left,left);
	    Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
	    /* *chrhigh += 1U; */
	  }
	  if (extend_rightward(/*goal*/left,chroffset,chrhigh,/*rightward_set*/plus_set,sarray,
			       query_compress_fwd,/*plusp*/true,/*genestrand*/0,
			       best_plus_elt->queryend) == true) {
	    collect_elt_matches(&(*found_score),&(*subs),&(*indels),&(*singlesplicing),&(*doublesplicing),
				best_plus_elt->querystart,best_plus_elt->queryend,
				chrnum,chroffset,chrhigh,chrlength,
				/*goal*/left,/*rightward_set*/plus_set,/*leftward_set*/NULL,
				querylength,query_compress_fwd,/*plusp*/true,/*genestrand*/0,
				nmisses_allowed,first_read_p);
	  }
	}
      }
    } else {
      /* Second (right) elt is best */
      debug(printf("Second elt is best\n"));
      plus_set = List_push(NULL,best_plus_elt);
      best_plus_elt = elt;
      best_plus_nmatches = nmatches;
      if (plus_querypos >= querylength) {
	chrhigh = 0U;
	Elt_fill_positions_all(best_plus_elt,sarray);
	for (i = 0; i < best_plus_elt->npositions; i++) {
	  left = best_plus_elt->positions[i];
	  if (left > chrhigh) {
	    chrnum = Univ_IIT_get_one(chromosome_iit,left,left);
	    Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
	    /* *chrhigh += 1U; */
	  }
	  nmatches = Genome_consecutive_matches_leftward(query_compress_fwd,left,
							 /*pos5*/0,/*pos3*/best_plus_elt->querystart,
							 /*plusp*/true,/*genestrand*/0);
	  debug(printf("Looking at position %u => %d matches leftward\n",left,nmatches));
	  best_plus_elt->querystart -= nmatches;
	  if (extend_leftward(/*goal*/left,chroffset,chrhigh,/*leftward_set*/plus_set,sarray,
			      /*queryptr*/queryuc_ptr,query_compress_fwd,
			      /*plusp*/true,/*genestrand*/0,
			      best_plus_elt->querystart,best_plus_elt->queryend) == true) {
	    collect_elt_matches(&(*found_score),&(*subs),&(*indels),&(*singlesplicing),&(*doublesplicing),
				best_plus_elt->querystart,best_plus_elt->queryend,
				chrnum,chroffset,chrhigh,chrlength,
				/*goal*/left,/*rightward_set*/NULL,/*leftward_set*/plus_set,
				querylength,query_compress_fwd,/*plusp*/true,/*genestrand*/0,
				nmisses_allowed,first_read_p);
	  }
	  best_plus_elt->querystart += nmatches;
	}
      }
    }
  }
    
  if (minus_querypos >= querylength) {
    minus_set = (List_T) NULL;
  } else {
    /* Extend minus side a second time */
    sarray_search(&initptr,&finalptr,&successp,&nmatches,&(queryrc[minus_querypos]),
		  querylength - minus_querypos,/*queryoffset*/minus_querypos,
		  query_compress_rev,/*plusp*/false);
    elt = Elt_new(minus_querypos,nmatches,initptr,finalptr);
    minus_querypos += nmatches;
    minus_querypos += 1;		/* To skip the presumed mismatch */

    debug(printf("minus_querypos %d vs querylength %d\n",minus_querypos,querylength));
    if (nmatches <= best_minus_nmatches) {
      /* Initial (left) elt was best */
      debug(printf("Initial elt was best\n"));
      minus_set = List_push(NULL,elt);
      if (minus_querypos >= querylength) {
	chrhigh = 0U;
	Elt_fill_positions_all(best_minus_elt,sarray);
	for (i = 0; i < best_minus_elt->npositions; i++) {
	  left = best_minus_elt->positions[i];
	  if (left > chrhigh) {
	    chrnum = Univ_IIT_get_one(chromosome_iit,left,left);
	    Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
	    /* *chrhigh += 1U; */
	  }
	  if (extend_rightward(/*goal*/left,chroffset,chrhigh,/*rightward_set*/minus_set,sarray,
			       query_compress_rev,/*plusp*/false,/*genestrand*/0,
			       best_minus_elt->queryend) == true) {
	    collect_elt_matches(&(*found_score),&(*subs),&(*indels),&(*singlesplicing),&(*doublesplicing),
				best_minus_elt->querystart,best_minus_elt->queryend,
				chrnum,chroffset,chrhigh,chrlength,
				/*goal*/left,/*rightward_set*/minus_set,/*leftward_set*/NULL,
				querylength,query_compress_rev,/*plusp*/false,/*genestrand*/0,
				nmisses_allowed,first_read_p);
	  }
	}
      }
    } else {
      /* Second (right) elt is best */
      debug(printf("Second elt is best\n"));
      minus_set = List_push(NULL,best_minus_elt);
      best_minus_elt = elt;
      best_minus_nmatches = nmatches;
      if (minus_querypos >= querylength) {
	chrhigh = 0U;
	Elt_fill_positions_all(best_minus_elt,sarray);
	for (i = 0; i < best_minus_elt->npositions; i++) {
	  left = best_minus_elt->positions[i];
	  if (left > chrhigh) {
	    chrnum = Univ_IIT_get_one(chromosome_iit,left,left);
	    Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
	    /* *chrhigh += 1U; */
	  }
	  nmatches = Genome_consecutive_matches_leftward(query_compress_rev,left,
							 /*pos5*/0,/*pos3*/best_minus_elt->querystart,
							 /*plusp*/false,/*genestrand*/0);
	  debug(printf(" extending bestelt querystart %d leftward by %d matches\n",best_minus_elt->querystart,nmatches));
	  best_minus_elt->querystart -= nmatches;
	  if (extend_leftward(/*goal*/left,chroffset,chrhigh,/*leftward_set*/minus_set,sarray,
			      /*queryptr*/queryrc,query_compress_rev,
			      /*plusp*/false,/*genestrand*/0,
			      best_minus_elt->querystart,best_minus_elt->queryend) == true) {
	    collect_elt_matches(&(*found_score),&(*subs),&(*indels),&(*singlesplicing),&(*doublesplicing),
				best_minus_elt->querystart,best_minus_elt->queryend,
				chrnum,chroffset,chrhigh,chrlength,
				/*goal*/left,/*rightward_set*/NULL,/*leftward_set*/minus_set,
				querylength,query_compress_rev,/*plusp*/false,/*genestrand*/0,
				nmisses_allowed,first_read_p);
	  }
	  best_minus_elt->querystart += nmatches;
	}
      }
    }
  }

  debug(printf("Found %d subs, %d indels, %d singlesplices, %d doublesplices\n",
	       List_length(*subs),List_length(*indels),List_length(*singlesplicing),List_length(*doublesplicing)));
  debug(printf("found_score %d vs querylength %d\n",*found_score,querylength));

  if (*found_score < querylength) {
    /* Be satisfied with a two-part alignment */
    if (plus_set != NULL) {
      elt = List_head(plus_set);
      Elt_free(&elt);
      List_free(&plus_set);
    }
    Elt_free(&best_plus_elt);
    if (minus_set != NULL) {
      elt = List_head(minus_set);
      Elt_free(&elt);
      List_free(&minus_set);
    }
    Elt_free(&best_minus_elt);
    return;
  } else {
    plus_set = List_push(plus_set,best_plus_elt);
    minus_set = List_push(minus_set,best_minus_elt);

#if 0
    /* Checking middle of read above */
    halfwaypos = querylength/2;
    if (best_plus_nmatches < halfwaypos) {
      /* Start from middle of read */
      debug(printf("Starting from halfway point on plus\n"));
      sarray_search(&initptr,&finalptr,&successp,&nmatches,&(queryuc_ptr[halfwaypos]),
		    querylength - halfwaypos,/*queryoffset*/halfwaypos,
		    query_compress_fwd,/*plusp*/true);
      elt = Elt_new(halfwaypos,nmatches,initptr,finalptr);
      if (nmatches > best_plus_nmatches) {
	best_plus_elt = elt;
	best_plus_nmatches = nmatches;
      }
      plus_set = List_push(plus_set,elt);
    }

    if (best_minus_nmatches < halfwaypos) {
      /* Start from middle of read */
      debug(printf("Starting from halfway point on minus\n"));
      sarray_search(&initptr,&finalptr,&successp,&nmatches,&(queryrc[halfwaypos]),
		    querylength - halfwaypos,/*queryoffset*/halfwaypos,
		    query_compress_rev,/*plusp*/false);
      elt = Elt_new(halfwaypos,nmatches,initptr,finalptr);
      if (nmatches > best_minus_nmatches) {
	best_minus_elt = elt;
	best_minus_nmatches = nmatches;
      }
      minus_set = List_push(minus_set,elt);
    }
#endif
  }


  plus_niter = minus_niter = 2;
  /* Both sides have failed and we don't have a good best hit.  Use up given allotment of attempts. */
  while (plus_querypos < querylength && plus_niter < nmisses_allowed) {
    sarray_search(&initptr,&finalptr,&successp,&nmatches,&(queryuc_ptr[plus_querypos]),
		  querylength - plus_querypos,/*queryoffset*/plus_querypos,
		  query_compress_fwd,/*plusp*/true);
    elt = Elt_new(plus_querypos,nmatches,initptr,finalptr);
    plus_set = List_push(plus_set,(void *) elt);
    if (nmatches > best_plus_nmatches) {
      best_plus_elt = elt;
      best_plus_nmatches = nmatches;

      /* See if we have a substitution winner */
      Elt_fill_positions_all(best_plus_elt,sarray);
      for (i = 0; i < best_plus_elt->npositions; i++) {
	left = best_plus_elt->positions[i];
	/* Should return max_mismatches + 1 if it exceeds the limit */
	if ((nmismatches = Genome_count_mismatches_limit(query_compress_fwd,left,/*pos5*/0,/*pos3*/querylength,
							 /*max_mismatches*/nmisses_allowed,/*plusp*/true,/*genestrand*/0)) <= nmisses_allowed) {
	  chrnum = Univ_IIT_get_one(chromosome_iit,left,left);
	  Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
	  if ((hit = Stage3end_new_substitution(&(*found_score),nmismatches,
						left,/*genomiclength*/querylength,
						query_compress_fwd,/*plusp*/true,/*genestrand*/0,
						chrnum,chroffset,chrhigh,chrlength,
						/*sarrayp*/true)) != NULL) {
	    *subs = List_push(*subs,(void *) hit);
	  }
	}
	debug(printf("Looking at plus position %u => %d mismatches\n",left,nmismatches));
      }
    }
    
    plus_querypos += nmatches;
    plus_querypos += 1;		/* To skip the presumed mismatch */
    plus_niter++;
  }

  while (minus_querypos < querylength && minus_niter < nmisses_allowed) {
    sarray_search(&initptr,&finalptr,&successp,&nmatches,&(queryrc[minus_querypos]),
		  querylength - minus_querypos,/*queryoffset*/minus_querypos,
		  query_compress_rev,/*plusp*/false);
    elt = Elt_new(minus_querypos,nmatches,initptr,finalptr);
    minus_set = List_push(minus_set,(void *) elt);
    if (nmatches > best_minus_nmatches) {
      best_minus_elt = elt;
      best_minus_nmatches = nmatches;

      /* See if we have a substitution winner */
      Elt_fill_positions_all(best_minus_elt,sarray);
      for (i = 0; i < best_minus_elt->npositions; i++) {
	left = best_minus_elt->positions[i];
	/* Should return max_mismatches + 1 if it exceeds the limit */
	if ((nmismatches = Genome_count_mismatches_limit(query_compress_rev,left,/*pos5*/0,/*pos3*/querylength,
							 /*max_mismatches*/nmisses_allowed,/*plusp*/false,/*genestrand*/0)) <= nmisses_allowed) {
	  chrnum = Univ_IIT_get_one(chromosome_iit,left,left);
	  Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
	  if ((hit = Stage3end_new_substitution(&(*found_score),nmismatches,
						left,/*genomiclength*/querylength,
						query_compress_rev,/*plusp*/false,/*genestrand*/0,
						chrnum,chroffset,chrhigh,chrlength,
						/*sarrayp*/true)) != NULL) {
	    *subs = List_push(*subs,(void *) hit);
	  }
	}
	debug(printf("Looking at minus position %u => %d mismatches\n",left,nmismatches));
      }
    }

    minus_querypos += nmatches;
    minus_querypos += 1;		/* To skip the presumed mismatch */
    minus_niter++;
  }

  debug(printf("Ended with %d plus iterations and %d minus iterations\n",plus_niter,minus_niter));

  if (plus_querypos >= querylength) {
    /* Handle plus extensions around best elt */
    debug(printf("BEST PLUS:\n"));
    debug(Elt_dump(best_plus_elt));

    leftward_set = rightward_set = (List_T) NULL;
    for (p = plus_set; p != NULL; p = p->rest) {
      elt = (Elt_T) p->first;
      if (elt == best_plus_elt) {
	/* Skip */

      } else if (elt->queryend < best_plus_elt->querystart) {
	leftward_set = List_push(leftward_set,(void *) elt);

      } else if (elt->querystart > best_plus_elt->queryend) {
	rightward_set = List_push(rightward_set,(void *) elt);

      } else {
	/* Duplicate -- skip */
      }
    }

    if ((nelts = List_length(rightward_set)) > 0) {
      array = (Elt_T *) List_to_array(rightward_set,NULL);
      List_free(&rightward_set);
      rightward_set = (List_T) NULL;
    
      qsort(array,nelts,sizeof(Elt_T),Elt_querypos_ascending_cmp);
      for (i = nelts-1; i >= 0; --i) {
	rightward_set = List_push(rightward_set,(void *) array[i]);
      }
      FREE(array);
    }

    if ((nelts = List_length(leftward_set)) > 0) {
      array = (Elt_T *) List_to_array(leftward_set,NULL);
      List_free(&leftward_set);
      leftward_set = (List_T) NULL;
    
      qsort(array,nelts,sizeof(Elt_T),Elt_querypos_descending_cmp);
      for (i = nelts-1; i >= 0; --i) {
	leftward_set = List_push(leftward_set,(void *) array[i]);
      }
      FREE(array);
    }


    chrhigh = 0U;
    Elt_fill_positions_all(best_plus_elt,sarray);
    for (i = 0; i < best_plus_elt->npositions; i++) {
      left = best_plus_elt->positions[i];
      if (left > chrhigh) {
	chrnum = Univ_IIT_get_one(chromosome_iit,left,left);
	Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
	/* *chrhigh += 1U; */
      }
      if (extend_rightward(/*goal*/left,chroffset,chrhigh,rightward_set,sarray,
			   query_compress_fwd,/*plusp*/true,/*genestrand*/0,
			   best_plus_elt->queryend) == true) {
	nmatches = Genome_consecutive_matches_leftward(query_compress_fwd,left,
						       /*pos5*/0,/*pos3*/best_plus_elt->querystart,
						       /*plusp*/true,/*genestrand*/0);
	debug(printf(" extending bestelt querystart %d leftward by %d matches\n",best_plus_elt->querystart,nmatches));
	best_plus_elt->querystart -= nmatches;
	if (extend_leftward(/*goal*/left,chroffset,chrhigh,leftward_set,sarray,
			    /*queryptr*/queryuc_ptr,query_compress_fwd,
			    /*plusp*/true,/*genestrand*/0,
			    best_plus_elt->querystart,best_plus_elt->queryend) == true) {
	  collect_elt_matches(&(*found_score),&(*subs),&(*indels),&(*singlesplicing),&(*doublesplicing),
			      best_plus_elt->querystart,best_plus_elt->queryend,
			      chrnum,chroffset,chrhigh,chrlength,
			      /*goal*/left,rightward_set,leftward_set,
			      querylength,query_compress_fwd,/*plusp*/true,/*genestrand*/0,
			      nmisses_allowed,first_read_p);
	}
	best_plus_elt->querystart += nmatches;
      }
    }

    List_free(&rightward_set);
    List_free(&leftward_set);
  }

  if (minus_querypos >= querylength) {
    /* Handle minus extensions around best elt */
    debug(printf("BEST MINUS:\n"));
    debug(Elt_dump(best_minus_elt));

    leftward_set = rightward_set = (List_T) NULL;
    for (p = minus_set; p != NULL; p = p->rest) {
      elt = (Elt_T) p->first;
      if (elt == best_minus_elt) {
	/* Skip */

      } else if (elt->queryend < best_minus_elt->querystart) {
	leftward_set = List_push(leftward_set,(void *) elt);

      } else if (elt->querystart > best_minus_elt->queryend) {
	rightward_set = List_push(rightward_set,(void *) elt);

      } else {
	/* Duplicate -- skip */
      }
    }

    if ((nelts = List_length(rightward_set)) > 0) {
      array = (Elt_T *) List_to_array(rightward_set,NULL);
      List_free(&rightward_set);
      rightward_set = (List_T) NULL;
    
      qsort(array,nelts,sizeof(Elt_T),Elt_querypos_ascending_cmp);
      for (i = nelts-1; i >= 0; --i) {
	rightward_set = List_push(rightward_set,(void *) array[i]);
      }
      FREE(array);
    }

    if ((nelts = List_length(leftward_set)) > 0) {
      array = (Elt_T *) List_to_array(leftward_set,NULL);
      List_free(&leftward_set);
      leftward_set = (List_T) NULL;
    
      qsort(array,nelts,sizeof(Elt_T),Elt_querypos_descending_cmp);
      for (i = nelts-1; i >= 0; --i) {
	leftward_set = List_push(leftward_set,(void *) array[i]);
      }
      FREE(array);
    }

    chrhigh = 0U;
    Elt_fill_positions_all(best_minus_elt,sarray);
    for (i = 0; i < best_minus_elt->npositions; i++) {
      left = best_minus_elt->positions[i];
      if (left > chrhigh) {
	chrnum = Univ_IIT_get_one(chromosome_iit,left,left);
	Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
	/* *chrhigh += 1U; */
      }
      if (extend_rightward(/*goal*/left,chroffset,chrhigh,rightward_set,sarray,
			   query_compress_rev,/*plusp*/false,/*genestrand*/0,
			   best_minus_elt->queryend) == true) {
	nmatches = Genome_consecutive_matches_leftward(query_compress_rev,left,
						       /*pos5*/0,/*pos3*/best_minus_elt->querystart,
						       /*plusp*/false,/*genestrand*/0);
	debug(printf(" extending bestelt querystart %d leftward by %d matches\n",best_minus_elt->querystart,nmatches));
	best_minus_elt->querystart -= nmatches;
	if (extend_leftward(/*goal*/left,chroffset,chrhigh,leftward_set,sarray,
			    /*queryptr*/queryrc,query_compress_rev,
			    /*plusp*/false,/*genestrand*/0,
			    best_minus_elt->querystart,best_minus_elt->queryend) == true) {
	  collect_elt_matches(&(*found_score),&(*subs),&(*indels),&(*singlesplicing),&(*doublesplicing),
			      best_minus_elt->querystart,best_minus_elt->queryend,
			      chrnum,chroffset,chrhigh,chrlength,
			      /*goal*/left,rightward_set,leftward_set,
			      querylength,query_compress_rev,/*plusp*/false,/*genestrand*/0,
			      nmisses_allowed,first_read_p);
	}
	best_minus_elt->querystart += nmatches;
      }
    }

    List_free(&rightward_set);
    List_free(&leftward_set);
  }

  for (p = plus_set; p != NULL; p = p->rest) {
    elt = (Elt_T) p->first;
    Elt_free(&elt);
  }
  List_free(&plus_set);

  for (p = minus_set; p != NULL; p = p->rest) {
    elt = (Elt_T) p->first;
    Elt_free(&elt);
  }
  List_free(&minus_set);

  debug(printf("Found %d subs, %d indels, %d singlesplices, %d doublesplices\n",
	       List_length(*subs),List_length(*indels),List_length(*singlesplicing),List_length(*doublesplicing)));

  return;
}

