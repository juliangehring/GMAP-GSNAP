static char rcsid[] = "$Id: sarray-read.c 151053 2014-10-16 19:57:23Z twu $";
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
#include "types.h"
#include "listdef.h"
#include "list.h"
#include "genome128_hr.h"
#include "splice.h"
#include "indel.h"
#include "stage3hr.h"
#include "bytecoding.h"
#include "bitpack64-read.h"
#include "bitpack64-readtwo.h"
#include "bitpack64-access.h"


#ifdef HAVE_SSE2
#include <emmintrin.h>
#endif
#ifdef HAVE_SSSE3
#include <tmmintrin.h>
#endif
#ifdef HAVE_POPCNT
#include <immintrin.h>
#elif defined(HAVE_MM_POPCNT)
#include <nmmintrin.h>
#endif


/* A value of 10000 misses various splices, although they are caught by GSNAP algorithm */
#define EXCESS_SARRAY_HITS 100000
#define LOCALSPLICING_SLOP 0.05

#define USE_SHUFFLE_MASK 1	/* Alternative requires AVX, and that part of the code isn't called much */

#define GUESS_ALLOCATION 10

/* #define USE_SEPARATE_BUCKETS 1 */

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

/* SIMD new filtering */
#ifdef DEBUG7B
#define debug7b(x) x
#else
#define debug7b(x)
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

/* Compare separate buckets with a single one */
#ifdef DEBUG15
#define debug15(x) x
#else
#define debug15(x)
#endif


#ifdef DEBUG7B
static void
print_vector_hex (__m128i x) {
  UINT4 *s = (UINT4 *) &x;

  /* printf("%08X %08X %08X %08X\n",s[0],s[1],s[2],s[3]); */
  printf("%08X %08X %08X %08X\n",s[3],s[2],s[1],s[0]);
  return;
}

static void
print_vector_uint (__m128i x) {
  UINT4 *s = (UINT4 *) &x;

  /* printf("%d %d %d %d\n",s[0],s[1],s[2],s[3]); */
  printf("%u %u %u %u\n",s[3],s[2],s[1],s[0]);
  return;
}
#endif



#define T Sarray_T
struct T {
  Univcoord_T n;
  Univcoord_T n_plus_one;

  Univcoord_T *array;

  unsigned char *lcpchilddc;

  UINT4 *lcp_guide;
  UINT4 *lcp_exceptions;
  int n_lcp_exceptions;		/* Won't be necessary if we change lcpchilddc to use guide array */
  /* int lcp_guide_interval; -- Always use 1024 */
  
  UINT4 *child_guide;
  UINT4 *child_exceptions;
  /* int n_child_exceptions; */
  int child_guide_interval; /* Always use 1024 */

#if 0
  Sarrayptr_T initindexi[4];	/* For A, C, G, T */
  Sarrayptr_T initindexj[4];	/* For A, C, G, T */
#endif

  int indexsize;
  UINT4 indexspace;		/* 4^indexsize.  Used by sarray_search to detect when we have a poly-T oligo shorter than indexsize */
#ifdef DEBUG15
  UINT4 *indexi_ptrs, *indexi_comp, *indexj_ptrs, *indexj_comp; /* bucket array: oligomer lookup into suffix array */
  UINT4 *indexij_ptrs, *indexij_comp;
#elif defined(USE_SEPARATE_BUCKETS)
  UINT4 *indexi_ptrs, *indexi_comp, *indexj_ptrs, *indexj_comp; /* bucket array: oligomer lookup into suffix array */
#else
  UINT4 *indexij_ptrs, *indexij_comp;
#endif

  Access_T sarray_access;
  Access_T aux_access;

  int array_fd; size_t array_len;
#ifdef DEBUG15
  int indexi_ptrs_fd; size_t indexi_ptrs_len; int indexi_comp_fd; size_t indexi_comp_len;
  int indexj_ptrs_fd; size_t indexj_ptrs_len; int indexj_comp_fd; size_t indexj_comp_len;
  int indexij_ptrs_fd; size_t indexij_ptrs_len; int indexij_comp_fd; size_t indexij_comp_len;
#elif defined(USE_SEPARATE_BUCKETS)
  int indexi_ptrs_fd; size_t indexi_ptrs_len; int indexi_comp_fd; size_t indexi_comp_len;
  int indexj_ptrs_fd; size_t indexj_ptrs_len; int indexj_comp_fd; size_t indexj_comp_len;
#else
  int indexij_ptrs_fd; size_t indexij_ptrs_len; int indexij_comp_fd; size_t indexij_comp_len;
#endif

  int lcpchilddc_fd; size_t lcpchilddc_len;

  int lcp_guide_fd; size_t lcp_guide_len;
  int lcp_exceptions_fd; size_t lcp_exceptions_len;

  int child_guide_fd; size_t child_guide_len;
  int child_exceptions_fd; size_t child_exceptions_len;

};


/* For benchmarking */
Univcoord_T
Sarray_size (Sarray_T this) {
  return this->n_plus_one;
}


static Sarray_T sarray_fwd;
static Sarray_T sarray_rev;
static Genome_T genome;

static char conversion_fwd[128];
static char conversion_rev[128];

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

#if defined(HAVE_SSE2) && defined(USE_SHUFFLE_MASK)
static __m128i shuffle_mask16[16];
#endif


#if 0
/* Simplified from sarray_search_simple in sarray-write.c */
static void
sarray_search_char (Sarrayptr_T *initptr, Sarrayptr_T *finalptr, char desired_char,
		    UINT4 *SA, UINT4 n, char *chartable) {
  Sarrayptr_T low, high, mid;
  Univcoord_T pos;
  char c;

  low = 1;
  high = n + 1;

  while (low < high) {
#if 0
    /* Compute mid for unsigned ints.  Want floor((low+high)/2). */
    mid = low/2 + high/2;
    if (low % 2 == 1 && high % 2 == 1) {
      mid += 1;
    }
#else
    mid = low + ((high - low) / 2);
#endif
    pos = SA[mid];
    c = Genome_get_char_lex(genome,pos,n,chartable);
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
#if 1
    /* Compute mid for unsigned ints.  Want ceil((low+high)/2). */
    mid = low/2 + high/2;
    if (low % 2 == 1 || high % 2 == 1) {
      mid += 1;
    }
#else
    /* This does not work for ceiling */
    mid = low + ((high - low) / 2);
#endif
    pos = SA[mid];
    c = Genome_get_char_lex(genome,pos,n,chartable);
    if (desired_char >= c) {
      low = mid;
    } else {
      high = mid - 1;
    }
  }

  *finalptr = high;
  return;
}
#endif


void
Sarray_setup (T sarray_fwd_in, T sarray_rev_in, Genome_T genome_in, Mode_T mode,
	      Univ_IIT_T chromosome_iit_in, int circular_typeint_in,
	      Chrpos_T shortsplicedist_in, int splicing_penalty_in,
	      int max_deletionlength, int max_end_deletions_in,
	      int max_middle_insertions, int max_end_insertions,
	      Univcoord_T *splicesites_in, Splicetype_T *splicetypes_in,
	      Chrpos_T *splicedists_in, int nsplicesites_in) {
  int i;

  sarray_fwd = sarray_fwd_in;
  sarray_rev = sarray_rev_in;
  genome = genome_in;

  for (i = 0; i < 128; i++) {
    conversion_fwd[i] = i;
    conversion_rev[i] = i;
  }
  if (mode == STANDARD) {
    /* Don't change conversion */
  } else if (mode == CMET_STRANDED || mode == CMET_NONSTRANDED) {
    conversion_fwd['C'] = 'T';	/* CT */
    conversion_rev['G'] = 'A';	/* GA */
  } else if (mode == ATOI_STRANDED || mode == ATOI_NONSTRANDED) {
    conversion_fwd['A'] = 'G';	/* AG */
    conversion_rev['T'] = 'C';	/* TC */
  }

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

#if 0
  sarray_search_char(&(sarray->initindexi[0]),&(sarray->initindexj[0]),/*desired_char*/'A',sarray->array,sarray->n);
  sarray_search_char(&(sarray->initindexi[1]),&(sarray->initindexj[1]),/*desired_char*/'C',sarray->array,sarray->n);
  sarray_search_char(&(sarray->initindexi[2]),&(sarray->initindexj[2]),/*desired_char*/'G',sarray->array,sarray->n);
  sarray_search_char(&(sarray->initindexi[3]),&(sarray->initindexj[3]),/*desired_char*/'T',sarray->array,sarray->n);
#endif

#if 0
  printf("A => %u %u\n",sarray->initindexi[0],sarray->initindexj[0]);
  printf("C => %u %u\n",sarray->initindexi[1],sarray->initindexj[1]);
  printf("G => %u %u\n",sarray->initindexi[2],sarray->initindexj[2]);
  printf("T => %u %u\n",sarray->initindexi[3],sarray->initindexj[3]);
#endif

#if defined(HAVE_SSE2) && defined(USE_SHUFFLE_MASK)
  /* Used by Elt_fill_positions_filtered */
  shuffle_mask16[0] =  _mm_set_epi8(-1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1);
  shuffle_mask16[1] =  _mm_set_epi8(-1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1,  3, 2, 1, 0);
  shuffle_mask16[2] =  _mm_set_epi8(-1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1,  7, 6, 5, 4);
  shuffle_mask16[3] =  _mm_set_epi8(-1,-1,-1,-1, -1,-1,-1,-1,  7, 6, 5, 4,  3, 2, 1, 0);
  shuffle_mask16[4] =  _mm_set_epi8(-1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1, 11,10, 9, 8);
  shuffle_mask16[5] =  _mm_set_epi8(-1,-1,-1,-1, -1,-1,-1,-1, 11,10, 9, 8,  3, 2, 1, 0);
  shuffle_mask16[6] =  _mm_set_epi8(-1,-1,-1,-1, -1,-1,-1,-1, 11,10, 9, 8,  7, 6, 5, 4);
  shuffle_mask16[7] =  _mm_set_epi8(-1,-1,-1,-1, 11,10, 9, 8,  7, 6, 5, 4,  3, 2, 1, 0);
  shuffle_mask16[8] =  _mm_set_epi8(-1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1, 15,14,13,12);
  shuffle_mask16[9] =  _mm_set_epi8(-1,-1,-1,-1, -1,-1,-1,-1, 15,14,13,12,  3, 2, 1, 0);
  shuffle_mask16[10] = _mm_set_epi8(-1,-1,-1,-1, -1,-1,-1,-1, 15,14,13,12,  7, 6, 5, 4);
  shuffle_mask16[11] = _mm_set_epi8(-1,-1,-1,-1, 15,14,13,12,  7, 6, 5, 4,  3, 2, 1, 0);
  shuffle_mask16[12] = _mm_set_epi8(-1,-1,-1,-1, -1,-1,-1,-1, 15,14,13,12, 11,10, 9, 8);
  shuffle_mask16[13] = _mm_set_epi8(-1,-1,-1,-1, 15,14,13,12, 11,10, 9, 8,  3, 2, 1, 0);
  shuffle_mask16[14] = _mm_set_epi8(-1,-1,-1,-1, 15,14,13,12, 11,10, 9, 8,  7, 6, 5, 4);
  shuffle_mask16[15] = _mm_set_epi8(15,14,13,12, 11,10, 9, 8,  7, 6, 5, 4,  3, 2, 1, 0);
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
Sarray_new (char *dir, char *fileroot, char *snps_root, Access_mode_T sarray_access, Access_mode_T aux_access,
	    Mode_T mode, bool fwdp) {
  T new;
  char *comma1;
  double seconds;
  int npages;

  char *sarrayfile;
  char *lcpchilddcfile;
  char *lcp_guidefile, *lcp_exceptionsfile;
  char *child_guidefile, *child_exceptionsfile;
#ifdef DEBUG15
  char *indexi_ptrsfile, *indexi_compfile;
  char *indexj_ptrsfile, *indexj_compfile;
  char *indexij_ptrsfile, *indexij_compfile;
#elif defined(USE_SEPARATE_BUCKETS)
  char *indexi_ptrsfile, *indexi_compfile;
  char *indexj_ptrsfile, *indexj_compfile;
#else
  char *indexij_ptrsfile, *indexij_compfile;
#endif

  char *mode_prefix;

  if (mode == STANDARD) {
    mode_prefix = ".";
  } else if (mode == CMET_STRANDED || mode == CMET_NONSTRANDED) {
    if (fwdp == true) {
      mode_prefix = ".metct.";
    } else {
      mode_prefix = ".metga.";
    }
  } else if (mode == ATOI_STRANDED || mode == ATOI_NONSTRANDED) {
    if (fwdp == true) {
      mode_prefix = ".a2iag.";
    } else {
      mode_prefix = ".a2itc.";
    }
  }

  sarrayfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(mode_prefix)+strlen("sarray")+1,sizeof(char));
  sprintf(sarrayfile,"%s/%s%ssarray",dir,fileroot,mode_prefix);

  lcpchilddcfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(mode_prefix)+strlen("salcpchilddc")+1,sizeof(char));
  sprintf(lcpchilddcfile,"%s/%s%ssalcpchilddc",dir,fileroot,mode_prefix);

  lcp_guidefile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(mode_prefix)+strlen("salcpguide1024")+1,sizeof(char));
  sprintf(lcp_guidefile,"%s/%s%ssalcpguide1024",dir,fileroot,mode_prefix);
  lcp_exceptionsfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(mode_prefix)+strlen("salcpexc")+1,sizeof(char));
  sprintf(lcp_exceptionsfile,"%s/%s%ssalcpexc",dir,fileroot,mode_prefix);

  child_guidefile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(mode_prefix)+strlen("sachildguide1024")+1,sizeof(char));
  sprintf(child_guidefile,"%s/%s%ssachildguide1024",dir,fileroot,mode_prefix);
  child_exceptionsfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(mode_prefix)+strlen("sachildexc")+1,sizeof(char));
  sprintf(child_exceptionsfile,"%s/%s%ssachildexc",dir,fileroot,mode_prefix);

#ifdef DEBUG15
  indexi_ptrsfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(".saindexi64meta")+1,sizeof(char));
  sprintf(indexi_ptrsfile,"%s/%s%ssaindexi64meta",dir,fileroot,mode_prefix);
  indexi_compfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(".saindexi64strm")+1,sizeof(char));
  sprintf(indexi_compfile,"%s/%s%ssaindexi64strm",dir,fileroot,mode_prefix);
  indexj_ptrsfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(".saindexj64meta")+1,sizeof(char));
  sprintf(indexj_ptrsfile,"%s/%s%ssaindexj64meta",dir,fileroot,mode_prefix);
  indexj_compfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(".saindexj64strm")+1,sizeof(char));
  sprintf(indexj_compfile,"%s/%s%ssaindexj64strm",dir,fileroot,mode_prefix);
  indexij_ptrsfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(".saindex64meta")+1,sizeof(char));
  sprintf(indexij_ptrsfile,"%s/%s%ssaindex64meta",dir,fileroot,mode_prefix);
  indexij_compfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(".saindex64strm")+1,sizeof(char));
  sprintf(indexij_compfile,"%s/%s%ssaindex64strm",dir,fileroot,mode_prefix);
#elif defined(USE_SEPARATE_BUCKETS)
  indexi_ptrsfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(".saindexi64meta")+1,sizeof(char));
  sprintf(indexi_ptrsfile,"%s/%s%ssaindexi64meta",dir,fileroot,mode_prefix);
  indexi_compfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(".saindexi64strm")+1,sizeof(char));
  sprintf(indexi_compfile,"%s/%s%ssaindexi64strm",dir,fileroot,mode_prefix);
  indexj_ptrsfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(".saindexj64meta")+1,sizeof(char));
  sprintf(indexj_ptrsfile,"%s/%s%ssaindexj64meta",dir,fileroot,mode_prefix);
  indexj_compfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(".saindexj64strm")+1,sizeof(char));
  sprintf(indexj_compfile,"%s/%s%ssaindexj64strm",dir,fileroot,mode_prefix);
#else
  indexij_ptrsfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(mode_prefix)+strlen("saindex64meta")+1,sizeof(char));
  sprintf(indexij_ptrsfile,"%s/%s%ssaindex64meta",dir,fileroot,mode_prefix);
  indexij_compfile = (char *) CALLOC(strlen(dir)+strlen("/")+strlen(fileroot)+strlen(mode_prefix)+strlen("saindex64strm")+1,sizeof(char));
  sprintf(indexij_compfile,"%s/%s%ssaindex64strm",dir,fileroot,mode_prefix);
#endif

  if (Access_file_exists_p(sarrayfile) == false) {
    fprintf(stderr,"Suffix array index file %s does not exist\n",sarrayfile);
    new = (T) NULL;

  } else if (Access_file_exists_p(lcpchilddcfile) == false) {
    fprintf(stderr,"Enhanced suffix array file %s does not exist.  The genome was built using an obsolete version\n",
	    lcpchilddcfile);
    new = (T) NULL;
    exit(9);

  } else {
    new = (T) MALLOC(sizeof(*new));

    if (sarray_access == USE_MMAP_PRELOAD) {
      fprintf(stderr,"Pre-loading suffix array...");
      new->array = (UINT4 *) Access_mmap_and_preload(&new->array_fd,&new->array_len,&npages,&seconds,sarrayfile,
						     sizeof(UINT4));
      new->sarray_access = MMAPPED;
      comma1 = Genomicpos_commafmt(new->array_len);
      fprintf(stderr,"done (%s bytes)\n",comma1);
      FREE(comma1);
    } else if (sarray_access == USE_MMAP_ONLY) {
      new->array = (UINT4 *) Access_mmap(&new->array_fd,&new->array_len,sarrayfile,sizeof(UINT4),/*randomp*/true);
      new->sarray_access = MMAPPED;
    } else if (sarray_access == USE_ALLOCATE) {
      new->array = (UINT4 *) Access_allocated(&new->array_len,&seconds,sarrayfile,sizeof(UINT4));
      new->sarray_access = ALLOCATED;
    }

    new->n_plus_one = new->array_len/sizeof(UINT4); /* Should be genomiclength + 1*/
    new->n = new->n_plus_one - 1;

#ifdef DEBUG15
    /* 8 is for two DIFFERENTIAL_METAINFO_SIZE words */
    new->indexi_ptrs = (UINT4 *) Access_allocated(&new->indexi_ptrs_len,&seconds,indexi_ptrsfile,sizeof(UINT4));
    new->indexi_comp = (UINT4 *) Access_allocated(&new->indexi_comp_len,&seconds,indexi_compfile,sizeof(UINT4));
    new->indexj_ptrs = (UINT4 *) Access_allocated(&new->indexj_ptrs_len,&seconds,indexj_ptrsfile,sizeof(UINT4));
    new->indexj_comp = (UINT4 *) Access_allocated(&new->indexj_comp_len,&seconds,indexj_compfile,sizeof(UINT4));
    new->indexij_ptrs = (UINT4 *) Access_allocated(&new->indexij_ptrs_len,&seconds,indexij_ptrsfile,sizeof(UINT4));
    new->indexij_comp = (UINT4 *) Access_allocated(&new->indexij_comp_len,&seconds,indexij_compfile,sizeof(UINT4));
    new->indexsize = 3 + log4(((new->indexij_ptrs_len - 8)/sizeof(UINT4)/2)/ /*DIFFERENTIAL_METAINFO_SIZE*/2);
#elif defined(USE_SEPARATE_BUCKETS)
    /* 8 is for two DIFFERENTIAL_METAINFO_SIZE words */
    new->indexi_ptrs = (UINT4 *) Access_allocated(&new->indexi_ptrs_len,&seconds,indexi_ptrsfile,sizeof(UINT4));
    new->indexi_comp = (UINT4 *) Access_allocated(&new->indexi_comp_len,&seconds,indexi_compfile,sizeof(UINT4));
    new->indexj_ptrs = (UINT4 *) Access_allocated(&new->indexj_ptrs_len,&seconds,indexj_ptrsfile,sizeof(UINT4));
    new->indexj_comp = (UINT4 *) Access_allocated(&new->indexj_comp_len,&seconds,indexj_compfile,sizeof(UINT4));
    new->indexsize = 3 + log4(((new->indexi_ptrs_len - 8)/sizeof(UINT4))/ /*DIFFERENTIAL_METAINFO_SIZE*/2);
#else
    /* 8 is for two DIFFERENTIAL_METAINFO_SIZE words */
    new->indexij_ptrs = (UINT4 *) Access_allocated(&new->indexij_ptrs_len,&seconds,indexij_ptrsfile,sizeof(UINT4));
    new->indexij_comp = (UINT4 *) Access_allocated(&new->indexij_comp_len,&seconds,indexij_compfile,sizeof(UINT4));
    new->indexsize = 3 + log4(((new->indexij_ptrs_len - 8)/sizeof(UINT4)/2)/ /*DIFFERENTIAL_METAINFO_SIZE*/2);
#endif
    new->indexspace = power(4,new->indexsize);

    if (aux_access == USE_MMAP_PRELOAD) {
      fprintf(stderr,"Pre-loading LCP/child/DC arrays...");
      new->lcpchilddc = (unsigned char *) Access_mmap_and_preload(&new->lcpchilddc_fd,&new->lcpchilddc_len,&npages,&seconds,
								  lcpchilddcfile,sizeof(unsigned char));
      new->aux_access = MMAPPED;
      comma1 = Genomicpos_commafmt(new->lcpchilddc_len);
      fprintf(stderr,"done (%s bytes)\n",comma1);
      FREE(comma1);
    } else if (aux_access == USE_MMAP_ONLY) {
      new->lcpchilddc = (unsigned char *) Access_mmap(&new->lcpchilddc_fd,&new->lcpchilddc_len,lcpchilddcfile,
						      sizeof(unsigned char),/*randomp*/true);
      new->aux_access = MMAPPED;
    } else if (aux_access == USE_ALLOCATE) {
      new->lcpchilddc = (unsigned char *) Access_allocated(&new->lcpchilddc_len,&seconds,lcpchilddcfile,sizeof(unsigned char));
      new->aux_access = ALLOCATED;
    }

    new->lcp_guide = (UINT4 *) Access_allocated(&new->lcp_guide_len,&seconds,lcp_guidefile,sizeof(UINT4));
    new->lcp_exceptions = (UINT4 *) Access_allocated(&new->lcp_exceptions_len,&seconds,lcp_exceptionsfile,sizeof(UINT4));
    new->n_lcp_exceptions = new->lcp_exceptions_len/(sizeof(UINT4) + sizeof(UINT4));

    new->child_guide = (UINT4 *) Access_allocated(&new->child_guide_len,&seconds,child_guidefile,sizeof(UINT4));
    new->child_exceptions = (UINT4 *) Access_allocated(&new->child_exceptions_len,&seconds,child_exceptionsfile,sizeof(UINT4));
    new->child_guide_interval = 1024;
  }


  FREE(child_exceptionsfile);
  FREE(child_guidefile);

  FREE(lcp_exceptionsfile);
  FREE(lcp_guidefile);

  FREE(lcpchilddcfile);

#ifdef DEBUG15
  FREE(indexi_compfile);
  FREE(indexi_ptrsfile);
  FREE(indexj_compfile);
  FREE(indexj_ptrsfile);
  FREE(indexij_compfile);
  FREE(indexij_ptrsfile);
#elif defined(USE_SEPARATE_BUCKETS)
  FREE(indexi_compfile);
  FREE(indexi_ptrsfile);
  FREE(indexj_compfile);
  FREE(indexj_ptrsfile);
#else
  FREE(indexij_compfile);
  FREE(indexij_ptrsfile);
#endif

  FREE(sarrayfile);

  return new;
}


void
Sarray_free (T *old) {
  if (*old) {
#ifdef DEBUG15
    FREE((*old)->indexi_ptrs);
    FREE((*old)->indexi_comp);
    FREE((*old)->indexj_ptrs);
    FREE((*old)->indexj_comp);
    FREE((*old)->indexij_ptrs);
    FREE((*old)->indexij_comp);
#elif defined(USE_SEPARATE_BUCKETS)
    FREE((*old)->indexi_ptrs);
    FREE((*old)->indexi_comp);
    FREE((*old)->indexj_ptrs);
    FREE((*old)->indexj_comp);
#else
    FREE((*old)->indexij_ptrs);
    FREE((*old)->indexij_comp);
#endif
    FREE((*old)->lcp_exceptions);
    FREE((*old)->lcp_guide);
    FREE((*old)->child_exceptions);
    FREE((*old)->child_guide);

    if ((*old)->aux_access == MMAPPED) {
      munmap((void *) (*old)->lcpchilddc,(*old)->lcpchilddc_len);
      close((*old)->lcpchilddc_fd);
    } else {
      FREE((*old)->lcpchilddc);
    }

    if ((*old)->sarray_access == MMAPPED) {
      munmap((void *) (*old)->array,(*old)->array_len);
      close((*old)->array_fd);
    } else {
      FREE((*old)->array);
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
#if 0
    /* Compute mid for unsigned ints */
    mid = low/2 + high/2;
    if (low % 2 == 1 && high % 2 == 1) {
      mid += 1;
    }
#else
    mid = low + ((high - low) / 2);
#endif

    debug1(printf("low %u, high %u => mid %u\n",low,high,mid));
    nmatches_mid =  (nmatches_low < nmatches_high) ? nmatches_low : nmatches_high;

    fasti = nmatches_mid +
      (Univcoord_T) Genome_consecutive_matches_rightward(query_compress,/*left*/sarray->array[mid]-queryoffset,
							 /*pos5*/queryoffset+nmatches_mid,
							 /*pos3*/queryoffset+querylength,plusp,genestrand,first_read_p);
    pos = sarray->array[mid] + fasti;
    c = Genome_get_char_lex(genome,pos,sarray->n,chartable);

    if (fasti == (Univcoord_T) querylength || c > query[fasti]) {
      high = mid;
      /* nmatches_high = (sarray->lcp[mid] < nmatches_mid) ? sarray->lcp[mid] : nmatches_mid; */
      sa_mid = sarray->array[mid];
      lcp_mid = Bitpack64_read_one(sa_mid,sarray->plcp_ptrs,sarray->plcp_comp) - sa_mid;
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
      lcp_low = Bitpack64_read_one(sa_low,sarray->plcp_ptrs,sarray->plcp_comp) - sa_low;
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
#if 0
    /* Compute mid for unsigned ints */
    mid = low/2 + high/2;
    if (low % 2 == 1 && high % 2 == 1) {
      mid += 1;
    }
#else
    mid = low + ((high - low) / 2);
#endif
    debug1(printf("low %u, high %u => mid %u\n",low,high,mid));
    nmatches_mid =  (nmatches_low < nmatches_high) ? nmatches_low : nmatches_high;

    fasti = nmatches_mid +
      (Univcoord_T) Genome_consecutive_matches_rightward(query_compress,/*left*/sarray->array[mid]-queryoffset,
							 /*pos5*/queryoffset+nmatches_mid,
							 /*pos3*/queryoffset+querylength,plusp,genestrand,first_read_p);
    pos = sarray->array[mid] + fasti;
    c = Genome_get_char_lex(genome,pos,sarray->n,chartable);

    if (fasti == (Univcoord_T) querylength || c < query[fasti]) {
      low = mid;
      /* nmatches_low = (sarray->lcp[low] < nmatches_mid) ? sarray->lcp[low] : nmatches_mid; */
      sa_low = sarray->array[low];
      lcp_low = Bitpack64_read_one(sa_low,sarray->plcp_ptrs,sarray->plcp_comp) - sa_low;
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
      lcp_mid = Bitpack64_read_one(sa_mid,sarray->plcp_ptrs,sarray->plcp_comp) - sa_mid;
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



/* For child[index+1].up, just calling child[index] */
#define decode_up(index,child_bytes,child_guide,child_exceptions,child_guide_interval) index - Bytecoding_read_wguide(index,child_bytes,child_guide,child_exceptions,child_guide_interval)
#define decode_down(index,child_bytes,child_guide,child_exceptions,child_guide_interval) Bytecoding_read_wguide(index,child_bytes,child_guide,child_exceptions,child_guide_interval) + index + 1
#define decode_next(index,child_bytes,child_guide,child_exceptions,child_guide_interval) Bytecoding_read_wguide(index,child_bytes,child_guide,child_exceptions,child_guide_interval) + index + 1

/*                                      0   1   2   3   4   5   6   7   8   9   A   B   C   D   E   F */
static char discrim_char_before[16] = {'?','$','$','$','$','$','A','A','A','A','C','C','C','G','G','T'};
static char discrim_char_after[16]  = {'?','A','C','G','T','X','C','G','T','X','G','T','X','T','X','X'};

static bool
get_child_given_first (Sarrayptr_T *l, Sarrayptr_T *r, Sarrayptr_T i, Sarrayptr_T j, char desired_char,
		       T sarray, unsigned char *lcpchilddc, UINT4 lcp_whole, UINT4 nextl) {
  char c1, c2;

  debug2(printf("Getting children for l-interval from %u to %u, char %c\n",i,j,desired_char));

#if 0
  /* First child already given */
  debug1(printf("lcp-interval %u..%u\n",i,j));
  up = decode_up(j,sarray->child_bytes,sarray->child_guide,sarray->child_exceptions,sarray->child_guide_interval);
  if (i < up && up <= j) {
    nextl = up;
    debug2(printf("nextl is up: %u\n",nextl));
  } else {
    nextl = decode_down(i,sarray->child_bytes,sarray->child_guide,sarray->child_exceptions,sarray->child_guide_interval); /* down */
    debug2(printf("nextl is down: %u\n",nextl));
  }
#endif

  /* Test first child: Use discrim_chars, rather than looking up S[SA[i] + lcp_whole] */
  c2 = Bytecoding_lcpchilddc_dc(&c1,nextl,lcpchilddc);
  debug2(printf("First child: %u to %u, discrim chars %c and %c\n",i,nextl-1,c1,c2));

  if (desired_char < c1) {
    debug2(printf("1.  Returning false, because desired %c < c1 %c\n",desired_char,c1));
    return false;
  } else if (desired_char == c1) {
    *l = i;
    *r = nextl - 1;
    debug2(printf("Returning true\n\n"));
    return true;
  } else if (desired_char < c2) {
    debug2(printf("1.  Returning false, because desired %c < c2 %c\n",desired_char,c2));
    return false;
  } else {
    /* Advance to middle children or final child */
    debug2(printf("1.  Advancing\n"));
  }

  /* Test for child[i] being down: lcp[child[i]] > lcp[i] */
  /* Test for child[i] being next_lindex: lcp[child[i]] == lcp[i] */
  /* Test middle children */
  while (nextl < j && Bytecoding_lcpchilddc_lcp_next(nextl,/*bytes*/lcpchilddc,sarray->child_guide,sarray->child_exceptions,
						     sarray->child_guide_interval,sarray->lcp_exceptions,sarray->n_lcp_exceptions) == lcp_whole) {
    /* Already tested for desired_char < c2 */
    if (desired_char == c2) {
      *l = nextl;
      *r = Bytecoding_lcpchilddc_child_next(nextl,lcpchilddc,sarray->child_guide,sarray->child_exceptions,
					    sarray->child_guide_interval) - 1; /* child[nextl] - 1 */
      debug2(printf("Child: %u to %u, c2 %c\n",nextl,*r,c2));
      debug2(printf("Returning true\n\n"));
      return true;
    } else {
      debug2(printf("Child: %u",nextl));
      nextl = Bytecoding_lcpchilddc_child_next(nextl,lcpchilddc,sarray->child_guide,sarray->child_exceptions,
					       sarray->child_guide_interval); /* child[nextl] */
      c2 = Bytecoding_lcpchilddc_dc(&c1,nextl,lcpchilddc);
      debug2(printf(" to %u, discrim chars %c and %c\n",nextl-1,c1,c2));

      if (desired_char < c2) {
	debug2(printf("M.  Returning false, because desired %c < c2 %c\n",desired_char,c2));
	return false;
      } else {
	debug2(printf("M.  Advancing\n"));
      }
    }
  }

  /* Test last child */
  /* Already tested for desired_char < c2 */
  debug2(printf("Final child: %u to %u, c2 %c\n",nextl,j,c2));
  if (desired_char == c2) {
    *l = nextl;
    *r = j;
    debug2(printf("Returning true\n\n"));
    return true;
  } else {
    debug2(printf("3.  Returning false, because desired %c != c2 %c\n",desired_char,c2));
    return false;
  }
}


static UINT4
find_longest_match (UINT4 nmatches, Sarrayptr_T *initptr, Sarrayptr_T *finalptr,
		    Sarrayptr_T i, Sarrayptr_T j, char *query, UINT4 querylength,
		    int queryoffset, Compress_T query_compress, T sarray, bool plusp,
		    int genestrand, bool first_read_p, char conversion[]) {
  UINT4 lcp_whole, nextl, up;
  UINT4 minlength;
  UINT4 l, r;

  while (nmatches < querylength) {
    if (i == j) {
      /* Singleton interval */
      debug1(printf("Singleton interval %u..%u\n",i,j));
      nmatches +=
	Genome_consecutive_matches_rightward(query_compress,/*left*/sarray->array[i]-queryoffset,
					     /*pos5*/queryoffset+nmatches,/*pos3*/queryoffset+querylength,
					     plusp,genestrand,first_read_p);
      *initptr = i;
      *finalptr = j;
      return nmatches;

    } else {
      /* First child */
      debug1(printf("lcp-interval %u..%u\n",i,j));
      up = Bytecoding_lcpchilddc_child_up(j,sarray->lcpchilddc,sarray->child_guide,sarray->child_exceptions,
					  sarray->child_guide_interval);
      if (i < up && up <= j) {
	nextl = up;
	debug2(printf("nextl is up: %u\n",nextl));
      } else {
	nextl = Bytecoding_lcpchilddc_child_next(i,sarray->lcpchilddc,sarray->child_guide,sarray->child_exceptions,
						 sarray->child_guide_interval); /* really down */
	debug2(printf("nextl is down: %u\n",nextl));
      }

      lcp_whole = Bytecoding_lcpchilddc_lcp(nextl,sarray->lcpchilddc,sarray->lcp_exceptions,
					    sarray->n_lcp_exceptions); /* lcp(i,j) */
      debug1(printf("lcp_whole for %u..%u is %d, compared with nmatches %d\n",i,j,lcp_whole,nmatches));

      if (lcp_whole > nmatches) {
	/* Check only up to minlength, so we validate the entire interval */
	minlength = (lcp_whole < querylength) ? lcp_whole : querylength;
	debug1(printf("Looking up genome for query from %d .. %d - 1\n",nmatches,minlength));
	nmatches +=
	  Genome_consecutive_matches_rightward(query_compress,/*left*/sarray->array[i]-queryoffset,
					       /*pos5*/queryoffset+nmatches,/*pos3*/queryoffset+minlength,
					       plusp,genestrand,first_read_p);
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
	
      debug1(printf("nmatches is now %d => desired_char is %c => %c\n",
		    nmatches,query[nmatches],conversion[query[nmatches]]));
      if (get_child_given_first(&l,&r,i,j,/*desired_char*/conversion[query[nmatches]],
				sarray,sarray->lcpchilddc,lcp_whole,nextl) == false) {
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
	       Compress_T query_compress, T sarray, bool plusp, int genestrand,
	       bool first_read_p, char conversion[]) {
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

  debug(printf("sarray_search on %s, querylength %d, plusp %d\n",query,querylength,plusp));

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
#ifdef DEBUG15
      if ((l = Bitpack64_read_one(oligo*2,sarray->indexij_ptrs,sarray->indexij_comp)) !=
	  Bitpack64_read_one(oligo,sarray->indexi_ptrs,sarray->indexi_comp)) {
	abort();
      }
#elif defined(USE_SEPARATE_BUCKETS)
      l = Bitpack64_read_one(oligo,sarray->indexi_ptrs,sarray->indexi_comp);
#else
      l = Bitpack64_read_one(oligo*2,sarray->indexij_ptrs,sarray->indexij_comp);
#endif
      debug1(printf(" => oligo %08X",oligo));
    }

    /* Because $ < A, we need to check for this case.  Need to back up just 1. */
    /* Test is SA[l-1] + indexsize - 1 >= n, or SA[l-1] + indexsize > n */
    debug1(printf("Comparing SA %u + indexsize %d with n %u\n",sarray->array[l-1],sarray->indexsize,sarray->n));
    if (l > 1 && sarray->array[l-1] + sarray->indexsize > sarray->n) {
      debug1(printf(" (backing up one position for l, because at end of genome)"));
      l--;
    }

    /* Add 1 to rollover to next oligo, to handle Ns in genome */
    oligo = nt_oligo_truncate(query,effective_querylength,sarray->indexsize,/*subst_value for T*/3) + 1;
    oligo_prev = oligo - 1;

    r = Bitpack64_read_one(oligo*2,sarray->indexij_ptrs,sarray->indexij_comp) - 1;
    r_prev = Bitpack64_offsetptr_only(oligo_prev*2+1,sarray->indexij_ptrs,sarray->indexij_comp);

    if (r != r_prev) {
      debug1(printf("r is %u - 1, but r_prev is %u, indicating the presence of N's => Starting from root\n",
		    r+1,r_prev));
      l = 1;
      r = sarray->n;

    } else if (oligo == sarray->indexspace) {
      /* We have a poly-T, so we cannot determine r.  For example,
	 TTTTTN has a different r value than TTN. */
      debug1(printf(" but poly-T => 1-letter for T: %u..%u\n",l,r));
#if 0
      l = sarray->initindexi[3];
      r = sarray->initindexj[3];
      /* Keep nmatches = 0, because there may not be a T in the genome */
#else
      l = 1;
      r = sarray->n;
#endif

    } else {
#if 0
      /* Already computed above */
#ifdef DEBUG15
      if ((r = Bitpack64_read_one(oligo*2,sarray->indexij_ptrs,sarray->indexij_comp) - 1) !=
	  Bitpack64_read_one(oligo,sarray->indexi_ptrs,sarray->indexi_comp) - 1) {
	abort();
      }
#elif defined(USE_SEPARATE_BUCKETS)
      r = Bitpack64_read_one(oligo,sarray->indexi_ptrs,sarray->indexi_comp) - 1;
#else
      r = Bitpack64_read_one(oligo*2,sarray->indexij_ptrs,sarray->indexij_comp) - 1;
#endif
#endif

      /* Because $ < A, we need to check for this case.  Need to back up just 1. */
      /* Test is SA[r] + indexsize - 1 >= n, or SA[r] + indexsize > n */
      debug1(printf(" (checking %u + %d >= %u)",sarray->array[r],effective_querylength,sarray->n));
      if (r > 0 && sarray->array[r] + sarray->indexsize > sarray->n) {
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
#ifdef DEBUG15
    if ((l = Bitpack64_read_two(&r,oligo*2,sarray->indexij_ptrs,sarray->indexij_comp)) !=
	Bitpack64_read_one(oligo,sarray->indexi_ptrs,sarray->indexi_comp)) {
      abort();
    } else if (r - 1 != Bitpack64_read_one(oligo,sarray->indexj_ptrs,sarray->indexj_comp)) {
      printf("For oligo %u, separate buckets give %u and %u, while single bucket gives %u and %u\n",
	     oligo,
	     Bitpack64_read_one(oligo,sarray->indexi_ptrs,sarray->indexi_comp),
	     Bitpack64_read_one(oligo,sarray->indexj_ptrs,sarray->indexj_comp),
	     l,r);
      abort();
    }
    r--;			/* Because interleaved writes r+1 to maintain monotonicity */
#elif defined(USE_SEPARATE_BUCKETS)
    l = Bitpack64_read_one(oligo,sarray->indexi_ptrs,sarray->indexi_comp);
    r = Bitpack64_read_one(oligo,sarray->indexj_ptrs,sarray->indexj_comp);
#else
    l = Bitpack64_read_two(&r,oligo*2,sarray->indexij_ptrs,sarray->indexij_comp);
    r--;			/* Because interleaved writes r+1 to maintain monotonicity */
#endif
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
				   query,querylength,queryoffset,query_compress,sarray,
				   plusp,genestrand,first_read_p,conversion);
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
						   plusp,genestrand,first_read_p);
    printf("%d\t%u\t%u\t",recount,(*initptr)-1,sarray->array[(*initptr)-1] /*+ 1U*/);
    if (genestrand == +2) {
      if (plusp) {
	Genome_fill_buffer_convert_rev(sarray->array[(*initptr)-1],recount+1,Buffer);
      } else {
	Genome_fill_buffer_convert_fwd(sarray->array[(*initptr)-1],recount+1,Buffer);
      }
    } else {
      if (plusp) {
	Genome_fill_buffer_convert_fwd(sarray->array[(*initptr)-1],recount+1,Buffer);
      } else {
	Genome_fill_buffer_convert_rev(sarray->array[(*initptr)-1],recount+1,Buffer);
      }
    }
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
						   plusp,genestrand,first_read_p);
    printf("%d\t%u\t%u\t",recount,(*initptr)+k,hit /*+ 1U*/);
    if (genestrand == +2) {
      if (plusp) {
	Genome_fill_buffer_convert_rev(sarray->array[(*initptr)+k],recount+1,Buffer);
      } else {
	Genome_fill_buffer_convert_fwd(sarray->array[(*initptr)+k],recount+1,Buffer);
      }
    } else {
      if (plusp) {
	Genome_fill_buffer_convert_fwd(sarray->array[(*initptr)+k],recount+1,Buffer);
      } else {
	Genome_fill_buffer_convert_rev(sarray->array[(*initptr)+k],recount+1,Buffer);
      }
    }
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
  if (*nmatches > 0 && sarray->array[(*finalptr)+1] > 0U) {
    printf("\n");
    recount = Genome_consecutive_matches_rightward(query_compress,/*left*/sarray->array[(*finalptr)+1]-queryoffset,
						   /*pos5*/queryoffset,/*pos3*/queryoffset+querylength,
						   plusp,genestrand,first_read_p);
    printf("%d\t%u\t%u\t",recount,(*finalptr)+1,sarray->array[(*finalptr)+1] /*+ 1U*/);
    if (genestrand == +2) {
      if (plusp) {
	Genome_fill_buffer_convert_rev(sarray->array[(*finalptr)+1],recount+1,Buffer);
      } else {
	Genome_fill_buffer_convert_fwd(sarray->array[(*finalptr)+1],recount+1,Buffer);
      }
    } else {
      if (plusp) {
	Genome_fill_buffer_convert_fwd(sarray->array[(*finalptr)+1],recount+1,Buffer);
      } else {
	Genome_fill_buffer_convert_rev(sarray->array[(*finalptr)+1],recount+1,Buffer);
      }
    }
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
    /* Can also happen with CMET or ATOI, since genome128_hr procedures find genome-to-query mismatches */
    /* abort(); */
  }
#endif

  return;
}




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

  debug(printf("Entering Elt_fill_positions_all on %p\n",this));
  if (this->positions_allocated != NULL) {
    debug(printf("  positions_allocated is already non-NULL, so skipping\n"));

  } else {
    this->npositions = this->finalptr - this->initptr + 1;
    debug(printf("  filling %d positions\n",this->npositions));

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
    for (i = 0; i < npositions; i++) {
      printf("%u\n",positions[i]);
    }
    printf("\n");

    for (i = 0; i < npositions_std; i++) {
      printf("%u\n",positions_std[i]);
    }
    printf("\n");
    abort();

  } else {
    qsort(positions,npositions,sizeof(Univcoord_T),Univcoord_compare);
    qsort(positions_std,npositions,sizeof(Univcoord_T),Univcoord_compare);
    for (i = 0; i < npositions; i++) {
      if (positions[i] != positions_std[i]) {
	fprintf(stderr,"At %d, positions %u != positions_std %u\n",i,positions[i],positions_std[i]);
	abort();
      }
    }
  }

  return;
}
#endif


#ifdef DEBUG8
static Univcoord_T *
fill_positions_std (int *npositions, Univcoord_T low_adj, Univcoord_T high_adj,
		    Sarrayptr_T initptr, Sarrayptr_T finalptr,
		    int querystart, Univcoord_T *array) {
  Univcoord_T *more_positions;
  Univcoord_T *positions, value;
  Sarrayptr_T ptr, lastptr;
  int i;

  positions = (Univcoord_T *) MALLOC(GUESS_ALLOCATION * sizeof(Univcoord_T)); /* Return value, so cannot use alloca */

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



#ifdef HAVE_ALLOCA

#if defined(HAVE_SSSE3) && defined(HAVE_SSE2)
/* SSSE3 needed for _mm_shuffle_epi8 */
static void
Elt_fill_positions_filtered (Elt_T this, T sarray, Univcoord_T goal, Univcoord_T low, Univcoord_T high,
			     Compress_T query_compress, bool plusp, int genestrand, bool first_read_p) {
  int nmatches;
  int i;
  Univcoord_T *array = sarray->array, low_adj, high_adj, value;
  Univcoord_T *positions_temp;
#ifdef HAVE_64_BIT
  UINT8 pointer;
#else
  UINT4 pointer;
#endif
  Sarrayptr_T *array_stop, *array_end, *array_ptr;
  Univcoord_T *out;
  __m128i converted, adjusted, match;
  __m128i base, floor, ceiling, values, adj, p;
  int matchbits;
#ifdef REQUIRE_ALIGNMENT
  int n_prealign, k;
#endif
#ifndef USE_SHUFFLE_MASK
  __m128i MASTER_CONTROL;
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
						    /*pos3*/this->queryend + 1,plusp,genestrand,first_read_p);
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


    positions_temp = out = (Univcoord_T *) MALLOCA((this->finalptr - this->initptr + 1) * sizeof(Univcoord_T));

    low_adj = low + this->querystart;
    high_adj = high + this->querystart;

    base = _mm_set1_epi32(2147483648); /* 2^31 */
    floor = _mm_set1_epi32(low_adj - 1 - 2147483648);
    ceiling = _mm_set1_epi32(high_adj + 1 - 2147483648);
    adj = _mm_set1_epi32(this->querystart);

    this->npositions = 0;
    array_ptr = &(array[this->initptr]);
    
#ifdef REQUIRE_ALIGNMENT
    /* Initial part */
#ifdef HAVE_64_BIT
    n_prealign = ((16 - ((UINT8) array_ptr & 0xF))/4) & 0x3;
#else
    n_prealign = ((16 - ((UINT4) array_ptr & 0xF))/4) & 0x3;
#endif
    debug7(printf("Initial ptr is at location %p.  Need %d to get to 128-bit boundary\n",pointer,n_prealign));

    debug7(printf("Initial part:\n"));
    if (n_prealign > this->finalptr - this->initptr + 1) {
      n_prealign = this->finalptr - this->initptr + 1;
    }
    for (k = 0; k < n_prealign; k++) {
      debug7a(printf("Looking at value %u, relative to low %u and high %u\n",array[ptr],low_adj,high_adj));
      if ((value = *array_ptr++) >= low_adj && value <= high_adj) {
	*out++ = value - this->querystart;
      }
    }
#endif

    /* Aligned part */
    if (this->finalptr < 4) {
      array_stop = &(array[0]);
    } else {
      array_stop = &(array[this->finalptr - 4]);
    }
    array_end = &(array[this->finalptr]);

#ifndef USE_SHUFFLE_MASK
    MASTER_CONTROL = _mm_setr_epi8(0x10, 0x12, 0x13, 0x12, 0x40, 0x68, 0x7C, 0x6B,
                                   0x00, 0x80, 0xC0, 0xBC, 0x00, 0x00, 0x00, 0xC0);
#endif

    while (array_ptr < array_stop) {
#ifdef REQUIRE_ALIGNMENT
      values = _mm_load_si128((__m128i *) array_ptr);
#else
      /* It looks like loadu is just as fast as load */
      values = _mm_loadu_si128((__m128i *) array_ptr);
#endif
      debug7b(print_vector_uint(values));

      converted = _mm_sub_epi32(values,base);
      /* match = _mm_andnot_si128(_mm_cmpgt_epi32(floor,converted),_mm_cmpgt_epi32(ceiling,converted)); -- This is off by 1 at floor */
      match = _mm_and_si128(_mm_cmpgt_epi32(converted,floor),_mm_cmplt_epi32(converted,ceiling));
      debug7b(print_vector_hex(match));

      matchbits = _mm_movemask_ps(_mm_castsi128_ps(match));
      if (matchbits) {
	adjusted = _mm_sub_epi32(values,adj);
#ifdef USE_SHUFFLE_MASK
	p = _mm_shuffle_epi8(adjusted, shuffle_mask16[matchbits]);
#else
	p = _mm_castps_si128(_mm_permutevar_ps(_mm_castsi128_ps(adjusted),_mm_srli_epi32(MASTER_CONTROL,matchbits*2)));
#endif
	_mm_storeu_si128((__m128i *) out, p);

#ifdef HAVE_POPCNT
	out += _popcnt32(matchbits);
	debug7b(printf("matchbits: %08X (%d ones)\n",matchbits,_popcnt32(matchbits)));
#elif defined HAVE_MM_POPCNT
	out += _mm_popcnt_u32(matchbits);
	debug7b(printf("matchbits: %08X (%d ones)\n",matchbits,_popcnt32(matchbits)));
#else
	out += __builtin_popcount(matchbits);
	debug7b(printf("matchbits: %08X (%d ones)\n",matchbits,__builtin_popcount(matchbits)));
#endif
      }

      array_ptr += 4;
    }

    /* Partial block at end; do scalar */
    debug7(printf("\nFinal part:\n"));
    while (array_ptr <= array_end) {
      if ((value = *array_ptr++) >= low_adj && value <= high_adj) {
	*out++ = value - this->querystart;
      }
    }

    this->npositions = out - positions_temp;

    debug7(printf("SIMD method found %d positions\n",this->npositions));
#ifdef DEBUG8
    positions_compare(positions_temp,this->npositions,positions_std,npositions_std);
    FREE(positions_std);
#endif

    /* Copy the positions from temp */
    if (this->npositions == 0) {
      this->positions_allocated = this->positions = (Univcoord_T *) NULL;
    } else {
      debug7(printf("Sorting %d positions\n",this->npositions));
      qsort(positions_temp,this->npositions,sizeof(Univcoord_T),Univcoord_compare);

      /* Need to copy positions before the goal */
      this->positions_allocated = this->positions = MALLOC(this->npositions * sizeof(Univcoord_T));
      memcpy(this->positions,positions_temp,this->npositions * sizeof(Univcoord_T));

      /* Advance pointer to goal (note: do not want goal_adj, since we have already subtracted this->querystart) */
      /* Have tested positions[i] <= goal, but want positions[-1] to be < goal, or positions[0] >= goal */
      /* ? Replace with a binary search */
      i = 0;
      while (i < this->npositions && positions_temp[i] < goal) {
	debug7(printf("Skipping position %u < goal %u\n",positions_temp[i],goal));
	i++;
      }
      this->positions += i;
      this->npositions -= i;
      debug7(printf("Remaining: %d positions\n",this->npositions));
    }
    
    FREEA(positions_temp);
  }

  this->filledp = true;

  return;
}

#else

static void
Elt_fill_positions_filtered (Elt_T this, T sarray, Univcoord_T goal, Univcoord_T low, Univcoord_T high,
			     Compress_T query_compress, bool plusp, int genestrand, bool first_read_p) {
  Sarrayptr_T ptr, lastptr;
  int nmatches;
  int i;
  Univcoord_T *array = sarray->array, low_adj, high_adj, value;
  Univcoord_T *positions_temp;
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
						    /*pos3*/this->queryend + 1,plusp,genestrand,first_read_p);
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

    positions_temp = (Univcoord_T *) MALLOCA((this->finalptr - this->initptr + 1) * sizeof(Univcoord_T));

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
	} else {
	  debug7(printf("Found position %u between low %u and high %u, and within allocation\n",value,low_adj,high_adj));
	  positions_temp[this->npositions++] = value - this->querystart;
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
	} else {
	  debug7(printf("Found position %u between low %u and high %u, and within allocation\n",value,low_adj,high_adj));
	  positions_temp[this->npositions++] = value - this->querystart;
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
	    } else {
	      debug7(printf("Found position %u between low %u and high %u, and within allocation\n",value,low_adj,high_adj));
	      positions_temp[this->npositions++] = value - this->querystart;
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
	} else {
	  debug7(printf("Found position %u between low %u and high %u, and within allocation\n",value,low_adj,high_adj));
	  positions_temp[this->npositions++] = value - this->querystart;
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
      } else {
	debug7(printf("Found position %u between low %u and high %u, and within allocation\n",value,low_adj,high_adj));
	positions_temp[this->npositions++] = value - this->querystart;
      }
    }
#endif

    debug7(printf("SIMD method found %d positions\n",this->npositions));
#ifdef DEBUG8
    positions_compare(positions_temp,this->npositions,positions_std,npositions_std);
    FREE(positions_std);
#endif

    /* Copy the positions from temp */
    if (this->npositions == 0) {
      this->positions_allocated = this->positions = (Univcoord_T *) NULL;
    } else {
      debug7(printf("Sorting %d positions\n",this->npositions));
      qsort(positions_temp,this->npositions,sizeof(Univcoord_T),Univcoord_compare);

      /* Need to copy positions before the goal */
      this->positions_allocated = this->positions = MALLOC(this->npositions * sizeof(Univcoord_T));
      memcpy(this->positions,positions_temp,this->npositions * sizeof(Univcoord_T));

      /* Advance pointer to goal (note: do not want goal_adj, since we have already subtracted this->querystart) */
      /* Have tested positions[i] <= goal, but want positions[-1] to be < goal, or positions[0] >= goal */
      /* ? Replace with a binary search */
      i = 0;
      while (i < this->npositions && positions_temp[i] < goal) {
	debug7(printf("Skipping position %u < goal %u\n",positions_temp[i],goal));
	i++;
      }
      this->positions += i;
      this->npositions -= i;
      debug7(printf("Remaining: %d positions\n",this->npositions));
    }
    
    FREEA(positions_temp);
  }

  this->filledp = true;

  return;
}
#endif


#else

static void
Elt_fill_positions_filtered (Elt_T this, T sarray, Univcoord_T goal, Univcoord_T low, Univcoord_T high,
			     Compress_T query_compress, bool plusp, int genestrand, bool first_read_p) {
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
						    /*pos3*/this->queryend + 1,plusp,genestrand,first_read_p);
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
  
#endif



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

  printf("Elt with %d positions:\n",elt->npositions);
  for (k = 0; k < elt->npositions; k++) {
    printf("  %d..%d:%u\n",elt->querystart,elt->queryend,elt->positions[k]);
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
		  List_T set, Compress_T query_compress,
		  T sarray, bool plusp, int genestrand, bool first_read_p, int best_queryend) {
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
      Elt_fill_positions_filtered(elt,sarray,goal,low,high,query_compress,plusp,genestrand,first_read_p);
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
      Elt_fill_positions_filtered(elt,sarray,goal,low,high,query_compress,plusp,genestrand,first_read_p);
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
		 List_T set, char *queryptr, Compress_T query_compress,
		 T sarray, bool plusp, int genestrand, bool first_read_p, char conversion[],
		 int best_querystart, int best_queryend) {
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
		    query_compress,sarray,plusp,genestrand,first_read_p,conversion);
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
      Elt_fill_positions_filtered(elt,sarray,goal,low,high,query_compress,plusp,genestrand,first_read_p);
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
							    plusp,genestrand,first_read_p)) > 0) {
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
      Elt_fill_positions_filtered(elt,sarray,goal,low,high,query_compress,plusp,genestrand,first_read_p);
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
							      plusp,genestrand,first_read_p)) > 0) {
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
							      plusp,genestrand,first_read_p)) > 0) {
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
							      plusp,genestrand,first_read_p)) > 0) {
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
							      plusp,genestrand,first_read_p)) > 0) {
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
collect_elt_matches (int *found_score, List_T *subs, List_T *indels, List_T *ambiguous, List_T *singlesplicing,
		     List_T *doublesplicing, int querystart_same, int queryend_same,
		     Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
		     Chrpos_T chrlength, Univcoord_T goal, 
		     List_T rightward_set, List_T leftward_set, int querylength, Compress_T query_compress,
		     bool plusp, int genestrand, bool first_read_p, int nmisses_allowed) {
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

  List_T spliceends_sense, spliceends_antisense, lowprob;
  int nhits, nspliceends_sense, nspliceends_antisense, n_good_spliceends;
  int best_nmismatches, nmismatches_donor, nmismatches_acceptor;
  double best_prob, prob;
  Substring_T donor, acceptor;

  int sensedir;
  Uintlist_T ambcoords, ambcoords_left, ambcoords_right;
  Intlist_T amb_knowni, amb_nmismatches;

  int segmenti_donor_nknown, segmentj_acceptor_nknown,
    segmentj_antidonor_nknown, segmenti_antiacceptor_nknown;
  int j, i, n;
  bool segmenti_usedp, segmentj_usedp;
  bool foundp;

#ifdef HAVE_ALLOCA
  int *segmenti_donor_knownpos = (int *) ALLOCA((querylength+1)*sizeof(int));
  int *segmentj_acceptor_knownpos = (int *) ALLOCA((querylength+1)*sizeof(int));
  int *segmentj_antidonor_knownpos = (int *) ALLOCA((querylength+1)*sizeof(int));
  int *segmenti_antiacceptor_knownpos = (int *) ALLOCA((querylength+1)*sizeof(int));
  int *segmenti_donor_knowni = (int *) ALLOCA((querylength+1)*sizeof(int));
  int *segmentj_acceptor_knowni = (int *) ALLOCA((querylength+1)*sizeof(int));
  int *segmentj_antidonor_knowni = (int *) ALLOCA((querylength+1)*sizeof(int));
  int *segmenti_antiacceptor_knowni = (int *) ALLOCA((querylength+1)*sizeof(int));
#else
  int segmenti_donor_knownpos[MAX_READLENGTH+1], segmentj_acceptor_knownpos[MAX_READLENGTH+1],
    segmentj_antidonor_knownpos[MAX_READLENGTH+1], segmenti_antiacceptor_knownpos[MAX_READLENGTH+1];
  int segmenti_donor_knowni[MAX_READLENGTH+1], segmentj_acceptor_knowni[MAX_READLENGTH+1],
    segmentj_antidonor_knowni[MAX_READLENGTH+1], segmenti_antiacceptor_knowni[MAX_READLENGTH+1];
#endif


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
						    plusp,genestrand,first_read_p);
    debug7(printf("nmismatches = %d (vs %d misses allowed)\n",nmismatches,nmisses_allowed));

    if (nmismatches > nmisses_allowed) {
      debug7(printf("Result: too many mismatches\n"));

    } else {
      debug7(printf("Result: successful hit saved\n"));
      debug(printf("Reporting hit with %d mismatches\n",nmismatches));
      if ((hit = Stage3end_new_substitution(&(*found_score),nmismatches,
					    left,/*genomiclength*/querylength,
					    query_compress,plusp,genestrand,first_read_p,
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

    n = Uintlist_length(difflist);
    array = (UINT4 *) MALLOCA(n * sizeof(UINT4));
    Uintlist_fill_array_and_free(array,&difflist);
    qsort(array,n,sizeof(Univcoord_T),Univcoord_compare);
    debug7(printf("Have %d matching diffs\n",n));

    spliceends_sense = spliceends_antisense = (List_T) NULL;
    lowprob = (List_T) NULL;
    for (i = 0; i < n; i++) {
      left2 = array[i];
      debug7(printf("diff %d/%d is at %u, from %d to %d\n",i,n,left2,querystart_diff - 1,queryend_diff));

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
	segmenti_donor_knownpos[segmenti_donor_nknown] = querylength;
	segmenti_antiacceptor_knownpos[segmenti_antiacceptor_nknown] = querylength;
	  
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
	segmentj_acceptor_knownpos[segmentj_acceptor_nknown] = querylength;
	segmentj_antidonor_knownpos[segmentj_antidonor_nknown] = querylength;

	/* nspliceends = 0; */
	spliceends_sense =
	  Splice_solve_single_sense(&(*found_score),&nspliceends_sense,spliceends_sense,&lowprob,
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
				    plusp,genestrand,first_read_p,/*subs_or_indels_p*/false,
				    /*sarrayp*/true);
	spliceends_antisense =
	  Splice_solve_single_antisense(&(*found_score),&nspliceends_antisense,spliceends_antisense,&lowprob,
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
				    plusp,genestrand,first_read_p,/*subs_or_indels_p*/false,
				    /*sarrayp*/true);

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
							   plusp,genestrand,first_read_p);
	  nmismatches2 = Genome_count_mismatches_substring(query_compress,left2,/*pos5*/indel_pos,
							   /*pos3*/querylength,plusp,genestrand,first_read_p);
	  if (plusp == true) {
	    query_indel_pos = indel_pos;
	  } else {
	    query_indel_pos = querylength - indel_pos;
	  }
	  if ((hit = Stage3end_new_deletion(&(*found_score),nindels,query_indel_pos,
					    nmismatches1,nmismatches2,
					    left1,/*genomiclength*/querylength+nindels,
					    query_compress,querylength,plusp,genestrand,first_read_p,
					    chrnum,chroffset,chrhigh,chrlength,
					    /*indel_penalty*/2,/*sarrayp*/true)) != NULL) {
	    debug7(printf("successful"));
	    *indels = List_push(*indels,(void *) hit);
	  }
#else
	  *indels = Indel_solve_middle_deletion(&foundp,&(*found_score),&nhits,*indels,
						/*left*/left1,chrnum,chroffset,chrhigh,chrlength,
						/*indels*/-nindels,query_compress,querylength,nmisses_allowed,
						plusp,genestrand,first_read_p,/*sarray*/true);
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
							   plusp,genestrand,first_read_p);
	  nmismatches2 = Genome_count_mismatches_substring(query_compress,left2,/*pos5*/indel_pos+nindels,
							   /*pos3*/querylength,plusp,genestrand,first_read_p);
	  if (plusp == true) {
	    query_indel_pos = indel_pos;
	  } else {
	    query_indel_pos = querylength - indel_pos - nindels;
	  }
	  if ((hit = Stage3end_new_insertion(&(*found_score),nindels,query_indel_pos,
					     nmismatches1,nmismatches2,
					     left1,/*genomiclength*/querylength-nindels,
					     query_compress,querylength,plusp,genestrand,first_read_p,
					     chrnum,chroffset,chrhigh,chrlength,
					     /*indel_penalty*/2,/*sarrayp*/true)) != NULL) {
	    debug7(printf("successful"));
	    *indels = List_push(*indels,(void *) hit);
	  }
#else
	  *indels = Indel_solve_middle_insertion(&foundp,&(*found_score),&nhits,*indels,
						 /*left*/left1,chrnum,chroffset,chrhigh,chrlength,
						 /*indels*/+nindels,query_compress,querylength,nmisses_allowed,
						 plusp,genestrand,first_read_p,/*sarrayp*/true);
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

    if (spliceends_sense != NULL) {
      /* nmismatches should be the same for all spliceends, so pick based on prob */
      hit = (Stage3end_T) List_head(spliceends_sense);
      best_nmismatches = Stage3end_nmismatches_whole(hit);

      best_prob = 0.0;
      for (p = spliceends_sense; p != NULL; p = List_next(p)) {
	hit = (Stage3end_T) List_head(p);
	debug7(printf("analyzing distance %d, probabilities %f and %f\n",
		      Stage3end_distance(hit),Substring_chimera_prob(Stage3end_substring_donor(hit)),
		      Substring_chimera_prob(Stage3end_substring_acceptor(hit))));
	if ((prob = Stage3end_chimera_prob(hit)) > best_prob) {
	  best_prob = prob;
	}
      }

      n_good_spliceends = 0;
      for (p = spliceends_sense; p != NULL; p = List_next(p)) {
	hit = (Stage3end_T) List_head(p);
	if (Stage3end_chimera_prob(hit) > best_prob - LOCALSPLICING_SLOP) {
	  debug7(printf("accepting distance %d, probabilities %f and %f\n",
			Stage3end_distance(hit),Substring_chimera_prob(Stage3end_substring_donor(hit)),
			Substring_chimera_prob(Stage3end_substring_acceptor(hit))));
	  n_good_spliceends += 1;
	}
      }

      if (n_good_spliceends == 1) {
	for (p = spliceends_sense; p != NULL; p = List_next(p)) {
	  hit = (Stage3end_T) List_head(p);
	  if (Stage3end_chimera_prob(hit) == best_prob) {
	    debug7(printf("pushing distance %d, probabilities %f and %f\n",
			  Stage3end_distance(hit),Substring_chimera_prob(Stage3end_substring_donor(hit)),
			  Substring_chimera_prob(Stage3end_substring_acceptor(hit))));
	    *singlesplicing = List_push(*singlesplicing,(void *) hit);
	    nhits += 1;
	  } else {
	    Stage3end_free(&hit);
	  }
	}
	List_free(&spliceends_sense);

      } else {
	/* Create ambiguous, sense */
	hit = (Stage3end_T) List_head(spliceends_sense);
	donor = Stage3end_substring_donor(hit);
	acceptor = Stage3end_substring_acceptor(hit);
	sensedir = Stage3end_sensedir(hit);

	ambcoords = (Uintlist_T) NULL;
	amb_knowni = (Intlist_T) NULL;
	amb_nmismatches = (Intlist_T) NULL;

	if (Substring_left_genomicseg(donor) == left1) {
	  for (p = spliceends_sense; p != NULL; p = List_next(p)) {
	    hit = (Stage3end_T) List_head(p);
	    acceptor = Stage3end_substring_acceptor(hit);
#ifdef LARGE_GENOMES
	    ambcoords = Uint8list_push(ambcoords,Substring_splicecoord(acceptor));
#else
	    ambcoords = Uintlist_push(ambcoords,Substring_splicecoord(acceptor));
#endif
	    amb_knowni = Intlist_push(amb_knowni,-1);
	    amb_nmismatches = Intlist_push(amb_nmismatches,Substring_nmismatches_whole(acceptor));
	  }

	  nmismatches_acceptor = best_nmismatches - Substring_nmismatches_whole(donor);
	  *ambiguous = List_push(*ambiguous,
				 (void *) Stage3end_new_splice(&(*found_score),
							       /*nmismatches_donor*/Substring_nmismatches_whole(donor),nmismatches_acceptor,
							       donor,/*acceptor*/NULL,/*distance*/0U,
							       /*shortdistancep*/false,/*penalty*/0,querylength,/*amb_nmatches*/Substring_nmatches_posttrim(acceptor),
							       /*ambcoords_donor*/NULL,ambcoords,
							       /*amb_knowni_donor*/NULL,amb_knowni,
							       /*amb_nmismatches_donor*/NULL,amb_nmismatches,
							       /*copy_donor_p*/true,/*copy_acceptor_p*/false,first_read_p,
							       sensedir,/*sarrayp*/true));
	  Intlist_free(&amb_nmismatches);
	  Intlist_free(&amb_knowni);
	  Uintlist_free(&ambcoords); /* LARGE_GENOMES not possible with suffix array */

	} else if (Substring_left_genomicseg(acceptor) == left1) {
	  for (p = spliceends_sense; p != NULL; p = List_next(p)) {
	    hit = (Stage3end_T) List_head(p);
	    donor = Stage3end_substring_donor(hit);
#ifdef LARGE_GENOMES
	    ambcoords = Uint8list_push(ambcoords,Substring_splicecoord(donor));
#else
	    ambcoords = Uintlist_push(ambcoords,Substring_splicecoord(donor));
#endif
	    amb_knowni = Intlist_push(amb_knowni,-1);
	    amb_nmismatches = Intlist_push(amb_nmismatches,Substring_nmismatches_whole(donor));
	  }

	  nmismatches_donor = best_nmismatches - Substring_nmismatches_whole(acceptor);
	  *ambiguous = List_push(*ambiguous,
				 (void *) Stage3end_new_splice(&(*found_score),
							       nmismatches_donor,/*nmismatches_acceptor*/Substring_nmismatches_whole(acceptor),
							       /*donor*/NULL,acceptor,/*distance*/0U,
							       /*shortdistancep*/false,/*penalty*/0,querylength,/*amb_nmatches*/Substring_nmatches_posttrim(donor),
							       ambcoords,/*ambcoords_acceptor*/NULL,
							       amb_knowni,/*amb_knowni_acceptor*/NULL,
							       amb_nmismatches,/*amb_nmismatches_acceptor*/NULL,
							       /*copy_donor_p*/false,/*copy_acceptor_p*/true,first_read_p,
							       sensedir,/*sarrayp*/true));
	  Intlist_free(&amb_nmismatches);
	  Intlist_free(&amb_knowni);
	  Uintlist_free(&ambcoords); /* LARGE_GENOMES not possible with suffix array */

	} else {
	  fprintf(stderr,"Unexpected: Neither donor left %u nor acceptor left %u equals left1 %u\n",
		  Substring_left_genomicseg(donor),Substring_left_genomicseg(acceptor),left1);
	  abort();
	}

	for (p = spliceends_sense; p != NULL; p = List_next(p)) {
	  hit = (Stage3end_T) List_head(p);
	  Stage3end_free(&hit);
	}
	List_free(&spliceends_sense);
      }
    }

    if (spliceends_antisense != NULL) {
      /* nmismatches should be the same for all spliceends, so pick based on prob */
      hit = (Stage3end_T) List_head(spliceends_antisense);
      best_nmismatches = Stage3end_nmismatches_whole(hit);

      best_prob = 0.0;
      for (p = spliceends_antisense; p != NULL; p = List_next(p)) {
	hit = (Stage3end_T) List_head(p);
	debug7(printf("analyzing distance %d, probabilities %f and %f\n",
		      Stage3end_distance(hit),Substring_chimera_prob(Stage3end_substring_donor(hit)),
		      Substring_chimera_prob(Stage3end_substring_acceptor(hit))));
	if ((prob = Stage3end_chimera_prob(hit)) > best_prob) {
	  best_prob = prob;
	}
      }

      n_good_spliceends = 0;
      for (p = spliceends_antisense; p != NULL; p = List_next(p)) {
	hit = (Stage3end_T) List_head(p);
	if (Stage3end_chimera_prob(hit) > best_prob - LOCALSPLICING_SLOP) {
	  debug7(printf("accepting distance %d, probabilities %f and %f\n",
			Stage3end_distance(hit),Substring_chimera_prob(Stage3end_substring_donor(hit)),
			Substring_chimera_prob(Stage3end_substring_acceptor(hit))));
	  n_good_spliceends += 1;
	}
      }

      if (n_good_spliceends == 1) {
	for (p = spliceends_antisense; p != NULL; p = List_next(p)) {
	  hit = (Stage3end_T) List_head(p);
	  if (Stage3end_chimera_prob(hit) == best_prob) {
	    debug7(printf("pushing distance %d, probabilities %f and %f\n",
			  Stage3end_distance(hit),Substring_chimera_prob(Stage3end_substring_donor(hit)),
			  Substring_chimera_prob(Stage3end_substring_acceptor(hit))));
	    *singlesplicing = List_push(*singlesplicing,(void *) hit);
	    nhits += 1;
	  } else {
	    Stage3end_free(&hit);
	  }
	}
	List_free(&spliceends_antisense);

      } else {
	/* Create ambiguous, antisense */
	hit = (Stage3end_T) List_head(spliceends_antisense);
	donor = Stage3end_substring_donor(hit);
	acceptor = Stage3end_substring_acceptor(hit);
	sensedir = Stage3end_sensedir(hit);

	ambcoords = (Uintlist_T) NULL;
	amb_knowni = (Intlist_T) NULL;
	amb_nmismatches = (Intlist_T) NULL;

	if (Substring_left_genomicseg(donor) == left1) {
	  for (p = spliceends_antisense; p != NULL; p = List_next(p)) {
	    hit = (Stage3end_T) List_head(p);
	    acceptor = Stage3end_substring_acceptor(hit);
#ifdef LARGE_GENOMES
	    ambcoords = Uint8list_push(ambcoords,Substring_splicecoord(acceptor));
#else
	    ambcoords = Uintlist_push(ambcoords,Substring_splicecoord(acceptor));
#endif
	    amb_knowni = Intlist_push(amb_knowni,-1);
	    amb_nmismatches = Intlist_push(amb_nmismatches,Substring_nmismatches_whole(acceptor));
	  }

	  nmismatches_acceptor = best_nmismatches - Substring_nmismatches_whole(donor);
	  *ambiguous = List_push(*ambiguous,
				 (void *) Stage3end_new_splice(&(*found_score),
							       /*nmismatches_donor*/Substring_nmismatches_whole(donor),nmismatches_acceptor,
							       donor,/*acceptor*/NULL,/*distance*/0U,
							       /*shortdistancep*/false,/*penalty*/0,querylength,/*amb_nmatches*/Substring_nmatches_posttrim(acceptor),
							       /*ambcoords_donor*/NULL,ambcoords,
							       /*amb_knowni_donor*/NULL,amb_knowni,
							       /*amb_nmismatches_donor*/NULL,amb_nmismatches,
							       /*copy_donor_p*/true,/*copy_acceptor_p*/false,first_read_p,
							       sensedir,/*sarrayp*/true));
	  Intlist_free(&amb_nmismatches);
	  Intlist_free(&amb_knowni);
	  Uintlist_free(&ambcoords); /* LARGE_GENOMES not possible with suffix array */

	} else if (Substring_left_genomicseg(acceptor) == left1) {
	  for (p = spliceends_antisense; p != NULL; p = List_next(p)) {
	    hit = (Stage3end_T) List_head(p);
	    donor = Stage3end_substring_donor(hit);
#ifdef LARGE_GENOMES
	    ambcoords = Uint8list_push(ambcoords,Substring_splicecoord(donor));
#else
	    ambcoords = Uintlist_push(ambcoords,Substring_splicecoord(donor));
#endif
	    amb_knowni = Intlist_push(amb_knowni,-1);
	    amb_nmismatches = Intlist_push(amb_nmismatches,Substring_nmismatches_whole(donor));
	  }

	  nmismatches_donor = best_nmismatches - Substring_nmismatches_whole(acceptor);
	  *ambiguous = List_push(*ambiguous,
				 (void *) Stage3end_new_splice(&(*found_score),
							       nmismatches_donor,/*nmismatches_acceptor*/Substring_nmismatches_whole(acceptor),
							       /*donor*/NULL,acceptor,/*distance*/0U,
							       /*shortdistancep*/false,/*penalty*/0,querylength,/*amb_nmatches*/Substring_nmatches_posttrim(donor),
							       ambcoords,/*ambcoords_acceptor*/NULL,
							       amb_knowni,/*amb_knowni_acceptor*/NULL,
							       amb_nmismatches,/*amb_nmismatches_acceptor*/NULL,
							       /*copy_donor_p*/false,/*copy_acceptor_p*/true,first_read_p,
							       sensedir,/*sarrayp*/true));
	  Intlist_free(&amb_nmismatches);
	  Intlist_free(&amb_knowni);
	  Uintlist_free(&ambcoords); /* LARGE_GENOMES not possible with suffix array */

	} else {
	  fprintf(stderr,"Unexpected: Neither donor left %u nor acceptor left %u equals left1 %u\n",
		  Substring_left_genomicseg(donor),Substring_left_genomicseg(acceptor),left1);
	  abort();
	}

	for (p = spliceends_antisense; p != NULL; p = List_next(p)) {
	  hit = (Stage3end_T) List_head(p);
	  Stage3end_free(&hit);
	}
	List_free(&spliceends_antisense);
      }
    }

    /* Don't use lowprob in suffix array stage */
    debug7(printf("freeing lowprobs\n"));
    for (p = lowprob; p != NULL; p = List_next(p)) {
      hit = (Stage3end_T) List_head(p);
      Stage3end_free(&hit);
    }
    List_free(&lowprob);

    FREEA(array);

  } else if (querystart_diff == 0 && queryend_same == querylength - 1) {
    left2 = left;
    indel_pos = querystart_same;
    debug7(printf("same is at %u from %d to %d\n",left,querystart_same,queryend_same));
    
    n = Uintlist_length(difflist);
    array = (UINT4 *) MALLOCA(n * sizeof(UINT4));
    Uintlist_fill_array_and_free(array,&difflist);
    qsort(array,n,sizeof(Univcoord_T),Univcoord_compare);
    debug7(printf("Have %d matching diffs\n",n));

    spliceends_sense = spliceends_antisense = (List_T) NULL;
    lowprob = (List_T) NULL;
    for (i = 0; i < n; i++) {
      left1 = array[i];
      debug7(printf("diff %d/%d is at %u, from %d to %d\n",i,n,left1,querystart_diff,queryend_diff));

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
	segmenti_donor_knownpos[segmenti_donor_nknown] = querylength;
	segmenti_antiacceptor_knownpos[segmenti_antiacceptor_nknown] = querylength;
	  
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
	segmentj_acceptor_knownpos[segmentj_acceptor_nknown] = querylength;
	segmentj_antidonor_knownpos[segmentj_antidonor_nknown] = querylength;

	/* nspliceends = 0; */
	spliceends_sense =
	  Splice_solve_single_sense(&(*found_score),&nspliceends_sense,spliceends_sense,&lowprob,
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
				    plusp,genestrand,first_read_p,/*subs_or_indels_p*/false,
				    /*sarrayp*/true);
	spliceends_antisense =
	  Splice_solve_single_antisense(&(*found_score),&nspliceends_antisense,spliceends_antisense,&lowprob,
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
				    plusp,genestrand,first_read_p,/*subs_or_indels_p*/false,
				    /*sarrayp*/true);

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
							   plusp,genestrand,first_read_p);
	  nmismatches2 = Genome_count_mismatches_substring(query_compress,left2,/*pos5*/indel_pos,
							   /*pos3*/querylength,plusp,genestrand,first_read_p);
	  if (plusp == true) {
	    query_indel_pos = indel_pos;
	  } else {
	    query_indel_pos = querylength - indel_pos;
	  }
	  if ((hit = Stage3end_new_deletion(&(*found_score),nindels,query_indel_pos,
					    nmismatches1,nmismatches2,
					    left1,/*genomiclength*/querylength+nindels,
					    query_compress,querylength,plusp,genestrand,first_read_p,
					    chrnum,chroffset,chrhigh,chrlength,
					    /*indel_penalty*/2,/*sarrayp*/true)) != NULL) {
	    debug7(printf("successful"));
	    *indels = List_push(*indels,(void *) hit);
	  }
#else
	  *indels = Indel_solve_middle_deletion(&foundp,&(*found_score),&nhits,*indels,
						/*left*/left1,chrnum,chroffset,chrhigh,chrlength,
						/*indels*/-nindels,query_compress,querylength,nmisses_allowed,
						plusp,genestrand,first_read_p,/*sarray*/true);
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
							   plusp,genestrand,first_read_p);
	  nmismatches2 = Genome_count_mismatches_substring(query_compress,left2,/*pos5*/indel_pos+nindels,
							   /*pos3*/querylength,plusp,genestrand,first_read_p);
	  if (plusp == true) {
	    query_indel_pos = indel_pos;
	  } else {
	    query_indel_pos = querylength - indel_pos - nindels;
	  }
	  if ((hit = Stage3end_new_insertion(&(*found_score),nindels,query_indel_pos,
					     nmismatches1,nmismatches2,
					     left1,/*genomiclength*/querylength-nindels,
					     query_compress,querylength,plusp,genestrand,first_read_p,
					     chrnum,chroffset,chrhigh,chrlength,
					     /*indel_penalty*/2,/*sarrayp*/true)) != NULL) {
	    debug7(printf("successful"));
	    *indels = List_push(*indels,(void *) hit);
	  }
#else
	  *indels = Indel_solve_middle_insertion(&foundp,&(*found_score),&nhits,*indels,
						 /*left*/left1,chrnum,chroffset,chrhigh,chrlength,
						 /*indels*/+nindels,query_compress,querylength,nmisses_allowed,
						 plusp,genestrand,first_read_p,/*sarrayp*/true);
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

    if (spliceends_sense != NULL) {
      /* nmismatches should be the same for all spliceends, so pick based on prob */
      hit = (Stage3end_T) List_head(spliceends_sense);
      best_nmismatches = Stage3end_nmismatches_whole(hit);

      best_prob = 0.0;
      for (p = spliceends_sense; p != NULL; p = List_next(p)) {
	hit = (Stage3end_T) List_head(p);
	debug7(printf("analyzing distance %d, probabilities %f and %f\n",
		      Stage3end_distance(hit),Substring_chimera_prob(Stage3end_substring_donor(hit)),
		      Substring_chimera_prob(Stage3end_substring_acceptor(hit))));
	if ((prob = Stage3end_chimera_prob(hit)) > best_prob) {
	  best_prob = prob;
	}
      }

      n_good_spliceends = 0;
      for (p = spliceends_sense; p != NULL; p = List_next(p)) {
	hit = (Stage3end_T) List_head(p);
	if (Stage3end_chimera_prob(hit) > best_prob - LOCALSPLICING_SLOP) {
	  debug7(printf("accepting distance %d, probabilities %f and %f\n",
			Stage3end_distance(hit),Substring_chimera_prob(Stage3end_substring_donor(hit)),
			Substring_chimera_prob(Stage3end_substring_acceptor(hit))));
	  n_good_spliceends += 1;
	}
      }
      
      if (n_good_spliceends == 1) {
	for (p = spliceends_sense; p != NULL; p = List_next(p)) {
	  hit = (Stage3end_T) List_head(p);
	  if (Stage3end_chimera_prob(hit) == best_prob) {
	    debug7(printf("pushing distance %d, probabilities %f and %f\n",
			  Stage3end_distance(hit),Substring_chimera_prob(Stage3end_substring_donor(hit)),
			  Substring_chimera_prob(Stage3end_substring_acceptor(hit))));
	    *singlesplicing = List_push(*singlesplicing,(void *) hit);
	    nhits += 1;
	  } else {
	    Stage3end_free(&hit);
	  }
	}
	List_free(&spliceends_sense);

      } else {
	/* Create ambiguous, sense */
	hit = (Stage3end_T) List_head(spliceends_sense);
	donor = Stage3end_substring_donor(hit);
	acceptor = Stage3end_substring_acceptor(hit);
	sensedir = Stage3end_sensedir(hit);

	ambcoords = (Uintlist_T) NULL;
	amb_knowni = (Intlist_T) NULL;
	amb_nmismatches = (Intlist_T) NULL;

	if (Substring_left_genomicseg(donor) == left2) {
	  for (p = spliceends_sense; p != NULL; p = List_next(p)) {
	    hit = (Stage3end_T) List_head(p);
	    acceptor = Stage3end_substring_acceptor(hit);
#ifdef LARGE_GENOMES
	    ambcoords = Uint8list_push(ambcoords,Substring_splicecoord(acceptor));
#else
	    ambcoords = Uintlist_push(ambcoords,Substring_splicecoord(acceptor));
#endif
	    amb_knowni = Intlist_push(amb_knowni,-1);
	    amb_nmismatches = Intlist_push(amb_nmismatches,Substring_nmismatches_whole(acceptor));
	  }

	  nmismatches_acceptor = best_nmismatches - Substring_nmismatches_whole(donor);
	  *ambiguous = List_push(*ambiguous,
				 (void *) Stage3end_new_splice(&(*found_score),
							       /*nmismatches_donor*/Substring_nmismatches_whole(donor),nmismatches_acceptor,
							       donor,/*acceptor*/NULL,/*distance*/0U,
							       /*shortdistancep*/false,/*penalty*/0,querylength,/*amb_nmatches*/Substring_nmatches_posttrim(acceptor),
							       /*ambcoords_donor*/NULL,ambcoords,
							       /*amb_knowni_donor*/NULL,amb_knowni,
							       /*amb_nmismatches_donor*/NULL,amb_nmismatches,
							       /*copy_donor_p*/true,/*copy_acceptor_p*/false,first_read_p,
							       sensedir,/*sarrayp*/true));
	  Intlist_free(&amb_nmismatches);
	  Intlist_free(&amb_knowni);
	  Uintlist_free(&ambcoords); /* LARGE_GENOMES not possible with suffix array */
	  
	} else if (Substring_left_genomicseg(acceptor) == left2) {
	  for (p = spliceends_sense; p != NULL; p = List_next(p)) {
	    hit = (Stage3end_T) List_head(p);
	    donor = Stage3end_substring_donor(hit);
#ifdef LARGE_GENOMES
	    ambcoords = Uint8list_push(ambcoords,Substring_splicecoord(donor));
#else
	    ambcoords = Uintlist_push(ambcoords,Substring_splicecoord(donor));
#endif
	    amb_knowni = Intlist_push(amb_knowni,-1);
	    amb_nmismatches = Intlist_push(amb_nmismatches,Substring_nmismatches_whole(donor));
	  }

	  nmismatches_donor = best_nmismatches - Substring_nmismatches_whole(acceptor);
	  *ambiguous = List_push(*ambiguous,
				 (void *) Stage3end_new_splice(&(*found_score),
							       nmismatches_donor,/*nmismatches_acceptor*/Substring_nmismatches_whole(acceptor),
							       /*donor*/NULL,acceptor,/*distance*/0U,
							       /*shortdistancep*/false,/*penalty*/0,querylength,/*amb_nmatches*/Substring_nmatches_posttrim(donor),
							       ambcoords,/*ambcoords_acceptor*/NULL,
							       amb_knowni,/*amb_known_acceptor*/NULL,
							       amb_nmismatches,/*amb_nmismatches_acceptor*/NULL,
							       /*copy_donor_p*/false,/*copy_acceptor_p*/true,first_read_p,
							       sensedir,/*sarrayp*/true));
	  Intlist_free(&amb_nmismatches);
	  Intlist_free(&amb_knowni);
	  Uintlist_free(&ambcoords); /* LARGE_GENOMES not possible with suffix array */

	} else {
	  fprintf(stderr,"Unexpected: Neither donor left %u nor acceptor left %u equals left2 %u\n",
		  Substring_left_genomicseg(donor),Substring_left_genomicseg(acceptor),left2);
	  abort();
	}

	for (p = spliceends_sense; p != NULL; p = List_next(p)) {
	  hit = (Stage3end_T) List_head(p);
	  Stage3end_free(&hit);
	}
	List_free(&spliceends_sense);
      }
    }

    if (spliceends_antisense != NULL) {
      /* nmismatches should be the same for all spliceends, so pick based on prob */
      hit = (Stage3end_T) List_head(spliceends_antisense);
      best_nmismatches = Stage3end_nmismatches_whole(hit);

      best_prob = 0.0;
      for (p = spliceends_antisense; p != NULL; p = List_next(p)) {
	hit = (Stage3end_T) List_head(p);
	debug7(printf("analyzing distance %d, probabilities %f and %f\n",
		      Stage3end_distance(hit),Substring_chimera_prob(Stage3end_substring_donor(hit)),
		      Substring_chimera_prob(Stage3end_substring_acceptor(hit))));
	if ((prob = Stage3end_chimera_prob(hit)) > best_prob) {
	  best_prob = prob;
	}
      }

      n_good_spliceends = 0;
      for (p = spliceends_antisense; p != NULL; p = List_next(p)) {
	hit = (Stage3end_T) List_head(p);
	if (Stage3end_chimera_prob(hit) > best_prob - LOCALSPLICING_SLOP) {
	  debug7(printf("accepting distance %d, probabilities %f and %f\n",
			Stage3end_distance(hit),Substring_chimera_prob(Stage3end_substring_donor(hit)),
			Substring_chimera_prob(Stage3end_substring_acceptor(hit))));
	  n_good_spliceends += 1;
	}
      }
      
      if (n_good_spliceends == 1) {
	for (p = spliceends_antisense; p != NULL; p = List_next(p)) {
	  hit = (Stage3end_T) List_head(p);
	  if (Stage3end_chimera_prob(hit) == best_prob) {
	    debug7(printf("pushing distance %d, probabilities %f and %f\n",
			  Stage3end_distance(hit),Substring_chimera_prob(Stage3end_substring_donor(hit)),
			  Substring_chimera_prob(Stage3end_substring_acceptor(hit))));
	    *singlesplicing = List_push(*singlesplicing,(void *) hit);
	    nhits += 1;
	  } else {
	    Stage3end_free(&hit);
	  }
	}
	List_free(&spliceends_antisense);

      } else {
	/* Create ambiguous, antisense */
	hit = (Stage3end_T) List_head(spliceends_antisense);
	donor = Stage3end_substring_donor(hit);
	acceptor = Stage3end_substring_acceptor(hit);
	sensedir = Stage3end_sensedir(hit);

	ambcoords = (Uintlist_T) NULL;
	amb_knowni = (Intlist_T) NULL;
	amb_nmismatches = (Intlist_T) NULL;

	if (Substring_left_genomicseg(donor) == left2) {
	  for (p = spliceends_antisense; p != NULL; p = List_next(p)) {
	    hit = (Stage3end_T) List_head(p);
	    acceptor = Stage3end_substring_acceptor(hit);
#ifdef LARGE_GENOMES
	    ambcoords = Uint8list_push(ambcoords,Substring_splicecoord(acceptor));
#else
	    ambcoords = Uintlist_push(ambcoords,Substring_splicecoord(acceptor));
#endif
	    amb_knowni = Intlist_push(amb_knowni,-1);
	    amb_nmismatches = Intlist_push(amb_nmismatches,Substring_nmismatches_whole(acceptor));
	  }

	  nmismatches_acceptor = best_nmismatches - Substring_nmismatches_whole(donor);
	  *ambiguous = List_push(*ambiguous,
				 (void *) Stage3end_new_splice(&(*found_score),
							       /*nmismatches_donor*/Substring_nmismatches_whole(donor),nmismatches_acceptor,
							       donor,/*acceptor*/NULL,/*distance*/0U,
							       /*shortdistancep*/false,/*penalty*/0,querylength,/*amb_nmatches*/Substring_nmatches_posttrim(acceptor),
							       /*ambcoords_donor*/NULL,ambcoords,
							       /*amb_knowni_donor*/NULL,amb_knowni,
							       /*amb_nmismatches_donor*/NULL,amb_nmismatches,
							       /*copy_donor_p*/true,/*copy_acceptor_p*/false,first_read_p,
							       sensedir,/*sarrayp*/true));
	  Intlist_free(&amb_nmismatches);
	  Intlist_free(&amb_knowni);
	  Uintlist_free(&ambcoords); /* LARGE_GENOMES not possible with suffix array */
	  
	} else if (Substring_left_genomicseg(acceptor) == left2) {
	  for (p = spliceends_antisense; p != NULL; p = List_next(p)) {
	    hit = (Stage3end_T) List_head(p);
	    donor = Stage3end_substring_donor(hit);
#ifdef LARGE_GENOMES
	    ambcoords = Uint8list_push(ambcoords,Substring_splicecoord(donor));
#else
	    ambcoords = Uintlist_push(ambcoords,Substring_splicecoord(donor));
#endif
	    amb_knowni = Intlist_push(amb_knowni,-1);
	    amb_nmismatches = Intlist_push(amb_nmismatches,Substring_nmismatches_whole(donor));
	  }

	  nmismatches_donor = best_nmismatches - Substring_nmismatches_whole(acceptor);
	  *ambiguous = List_push(*ambiguous,
				 (void *) Stage3end_new_splice(&(*found_score),
							       nmismatches_donor,/*nmismatches_acceptor*/Substring_nmismatches_whole(acceptor),
							       /*donor*/NULL,acceptor,/*distance*/0U,
							       /*shortdistancep*/false,/*penalty*/0,querylength,/*amb_nmatches*/Substring_nmatches_posttrim(donor),
							       ambcoords,/*ambcoords_acceptor*/NULL,
							       amb_knowni,/*amb_known_acceptor*/NULL,
							       amb_nmismatches,/*amb_nmismatches_acceptor*/NULL,
							       /*copy_donor_p*/false,/*copy_acceptor_p*/true,first_read_p,
							       sensedir,/*sarrayp*/true));
	  Intlist_free(&amb_nmismatches);
	  Intlist_free(&amb_knowni);
	  Uintlist_free(&ambcoords); /* LARGE_GENOMES not possible with suffix array */

	} else {
	  fprintf(stderr,"Unexpected: Neither donor left %u nor acceptor left %u equals left2 %u\n",
		  Substring_left_genomicseg(donor),Substring_left_genomicseg(acceptor),left2);
	  abort();
	}

	for (p = spliceends_antisense; p != NULL; p = List_next(p)) {
	  hit = (Stage3end_T) List_head(p);
	  Stage3end_free(&hit);
	}
	List_free(&spliceends_antisense);
      }
    }


    /* Don't use lowprob in suffix array stage */
    debug7(printf("freeing lowprobs\n"));
    for (p = lowprob; p != NULL; p = List_next(p)) {
      hit = (Stage3end_T) List_head(p);
      Stage3end_free(&hit);
    }
    List_free(&lowprob);

    FREEA(array);

  } else {
    Uintlist_free(&difflist);
  }

  return;
}


void
Sarray_search_greedy (int *found_score, List_T *subs, List_T *indels, List_T *ambiguous, List_T *singlesplicing,
		      List_T *doublesplicing, char *queryuc_ptr, char *queryrc, int querylength,
		      Compress_T query_compress_fwd, Compress_T query_compress_rev,
		      int nmisses_allowed, int genestrand, bool first_read_p) {
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
  T plus_sarray, minus_sarray;
  char *plus_conversion, *minus_conversion;


  if (nmisses_allowed < 0) {
    nmisses_allowed = 0;
  }
  debug(printf("\nStarting Sarray_search_greedy with querylength %d and indexsize %d and nmisses_allowed %d\n",
	       querylength,sarray_fwd->indexsize,nmisses_allowed));
  debug(printf("genestrand = %d\n",genestrand));

  *found_score = querylength;

  if (genestrand == +2) {
    plus_conversion = conversion_rev;
    minus_conversion = conversion_fwd;
    plus_sarray = sarray_rev;
    minus_sarray = sarray_fwd;
  } else {
    plus_conversion = conversion_fwd;
    minus_conversion = conversion_rev;
    plus_sarray = sarray_fwd;
    minus_sarray = sarray_rev;
  }

  /* Do one plus round */
  plus_querypos = 0;
  sarray_search(&initptr,&finalptr,&successp,&best_plus_nmatches,&(queryuc_ptr[plus_querypos]),
		querylength - plus_querypos,/*queryoffset*/plus_querypos,
		query_compress_fwd,plus_sarray,/*plusp*/true,genestrand,first_read_p,plus_conversion);
  best_plus_elt = Elt_new(plus_querypos,best_plus_nmatches,initptr,finalptr);
  plus_querypos += (int) best_plus_nmatches;
  plus_querypos += 1;		/* To skip the presumed mismatch */


  /* Do one minus round */
  minus_querypos = 0;
  sarray_search(&initptr,&finalptr,&successp,&best_minus_nmatches,&(queryrc[minus_querypos]),
		querylength - minus_querypos,/*queryoffset*/minus_querypos,
		query_compress_rev,minus_sarray,/*plusp*/false,genestrand,first_read_p,minus_conversion);
  best_minus_elt = Elt_new(minus_querypos,best_minus_nmatches,initptr,finalptr);
  minus_querypos += (int) best_minus_nmatches;
  minus_querypos += 1;		/* To skip the presumed mismatch */


  if (best_plus_nmatches >= querylength/2) {
    /* See if we have a winner */
    debug(printf("best_plus_nmatches = %d > %d/2, so checking mismatches against %d allowed\n",
		 best_plus_nmatches,querylength,nmisses_allowed));
    Elt_fill_positions_all(best_plus_elt,plus_sarray);
    for (i = 0; i < best_plus_elt->npositions; i++) {
      left = best_plus_elt->positions[i];
      /* Should return max_mismatches + 1 if it exceeds the limit */
      if ((nmismatches = Genome_count_mismatches_limit(query_compress_fwd,left,/*pos5*/0,/*pos3*/querylength,
						       /*max_mismatches*/nmisses_allowed,
						       /*plusp*/true,genestrand,first_read_p)) <= nmisses_allowed) {
	chrnum = Univ_IIT_get_one(chromosome_iit,left,left);
	Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
	debug(printf("Case 1: New substitution from beginning\n"));
	if ((hit = Stage3end_new_substitution(&(*found_score),nmismatches,
					      left,/*genomiclength*/querylength,
					      query_compress_fwd,/*plusp*/true,genestrand,first_read_p,
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
		  query_compress_fwd,plus_sarray,/*plusp*/true,genestrand,first_read_p,plus_conversion);
    /* Don't want to limit based on nmatches */
    if (1 || nmatches >= querylength - halfwaypos) {
      elt = Elt_new(halfwaypos,nmatches,initptr,finalptr);
      Elt_fill_positions_all(elt,plus_sarray);
      for (i = 0; i < elt->npositions; i++) {
	left = elt->positions[i];
	/* Should return max_mismatches + 1 if it exceeds the limit */
	if ((nmismatches = Genome_count_mismatches_limit(query_compress_fwd,left,/*pos5*/0,/*pos3*/querylength,
							 /*max_mismatches*/nmisses_allowed,
							 /*plusp*/true,genestrand,first_read_p)) <= nmisses_allowed) {
	  chrnum = Univ_IIT_get_one(chromosome_iit,left,left);
	  Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
	  debug(printf("Case 1: New substitution from middle\n"));
	  if ((hit = Stage3end_new_substitution(&(*found_score),nmismatches,
						left,/*genomiclength*/querylength,
						query_compress_fwd,/*plusp*/true,genestrand,first_read_p,
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
    Elt_fill_positions_all(best_minus_elt,minus_sarray);
    for (i = 0; i < best_minus_elt->npositions; i++) {
      left = best_minus_elt->positions[i];
      /* Should return max_mismatches + 1 if it exceeds the limit */
      if ((nmismatches = Genome_count_mismatches_limit(query_compress_rev,left,/*pos5*/0,/*pos3*/querylength,
						       /*max_mismatches*/nmisses_allowed,
						       /*plusp*/false,genestrand,first_read_p)) <= nmisses_allowed) {
	chrnum = Univ_IIT_get_one(chromosome_iit,left,left);
	Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
	debug(printf("Case 2: New substitution from beginning\n"));
	if ((hit = Stage3end_new_substitution(&(*found_score),nmismatches,
					      left,/*genomiclength*/querylength,
					      query_compress_rev,/*plusp*/false,genestrand,first_read_p,
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
		  query_compress_rev,minus_sarray,/*plusp*/false,genestrand,first_read_p,minus_conversion);
    /* Don't want to limit based on nmatches */
    if (1 || nmatches >= querylength - halfwaypos) {
      elt = Elt_new(halfwaypos,nmatches,initptr,finalptr);
      Elt_fill_positions_all(elt,minus_sarray);
      for (i = 0; i < elt->npositions; i++) {
	left = elt->positions[i];
	/* Should return max_mismatches + 1 if it exceeds the limit */
	if ((nmismatches = Genome_count_mismatches_limit(query_compress_rev,left,/*pos5*/0,/*pos3*/querylength,
							 /*max_mismatches*/nmisses_allowed,
							 /*plusp*/false,genestrand,first_read_p)) <= nmisses_allowed) {
	  chrnum = Univ_IIT_get_one(chromosome_iit,left,left);
	  Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
	  debug(printf("Case 2: New substitution from middle\n"));
	  if ((hit = Stage3end_new_substitution(&(*found_score),nmismatches,
						left,/*genomiclength*/querylength,
						query_compress_rev,/*plusp*/false,genestrand,first_read_p,
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
		  query_compress_fwd,plus_sarray,/*plusp*/true,genestrand,first_read_p,plus_conversion);
    elt = Elt_new(plus_querypos,nmatches,initptr,finalptr);
    plus_querypos += nmatches;
    plus_querypos += 1;		/* To skip the presumed mismatch */

    debug(printf("plus_querypos %d vs querylength %d\n",plus_querypos,querylength));
    if (nmatches <= best_plus_nmatches) {
      /* Initial (left) elt was best */
      debug(printf("Initial elt %p was best:\n",best_plus_elt));
      plus_set = List_push(NULL,elt);
      if (plus_querypos >= querylength) {
	chrhigh = 0U;
	Elt_fill_positions_all(best_plus_elt,plus_sarray);
	debug(Elt_dump(best_plus_elt));
	for (i = 0; i < best_plus_elt->npositions; i++) {
	  left = best_plus_elt->positions[i];
	  if (left > chrhigh) {
	    chrnum = Univ_IIT_get_one(chromosome_iit,left,left);
	    Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
	    /* *chrhigh += 1U; */
	  }
	  if (extend_rightward(/*goal*/left,chroffset,chrhigh,/*rightward_set*/plus_set,
			       query_compress_fwd,plus_sarray,/*plusp*/true,genestrand,first_read_p,
			       best_plus_elt->queryend) == true) {
	    collect_elt_matches(&(*found_score),&(*subs),&(*indels),&(*ambiguous),&(*singlesplicing),&(*doublesplicing),
				best_plus_elt->querystart,best_plus_elt->queryend,
				chrnum,chroffset,chrhigh,chrlength,
				/*goal*/left,/*rightward_set*/plus_set,/*leftward_set*/NULL,
				querylength,query_compress_fwd,/*plusp*/true,genestrand,first_read_p,
				nmisses_allowed);
	  }
	}
      }
    } else {
      /* Second (right) elt is best */
      debug(printf("Second elt %p is best:\n",elt));
      plus_set = List_push(NULL,best_plus_elt);
      best_plus_elt = elt;
      best_plus_nmatches = nmatches;
      if (plus_querypos >= querylength) {
	chrhigh = 0U;
	Elt_fill_positions_all(best_plus_elt,plus_sarray);
	debug(Elt_dump(best_plus_elt));
	for (i = 0; i < best_plus_elt->npositions; i++) {
	  left = best_plus_elt->positions[i];
	  if (left > chrhigh) {
	    chrnum = Univ_IIT_get_one(chromosome_iit,left,left);
	    Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
	    /* *chrhigh += 1U; */
	  }
	  nmatches = Genome_consecutive_matches_leftward(query_compress_fwd,left,
							 /*pos5*/0,/*pos3*/best_plus_elt->querystart,
							 /*plusp*/true,genestrand,first_read_p);
	  debug(printf("Looking at position %u => %d matches leftward\n",left,nmatches));
	  best_plus_elt->querystart -= nmatches;
	  if (extend_leftward(/*goal*/left,chroffset,chrhigh,/*leftward_set*/plus_set,
			      /*queryptr*/queryuc_ptr,query_compress_fwd,
			      plus_sarray,/*plusp*/true,genestrand,first_read_p,plus_conversion,
			      best_plus_elt->querystart,best_plus_elt->queryend) == true) {
	    collect_elt_matches(&(*found_score),&(*subs),&(*indels),&(*ambiguous),&(*singlesplicing),&(*doublesplicing),
				best_plus_elt->querystart,best_plus_elt->queryend,
				chrnum,chroffset,chrhigh,chrlength,
				/*goal*/left,/*rightward_set*/NULL,/*leftward_set*/plus_set,
				querylength,query_compress_fwd,/*plusp*/true,genestrand,first_read_p,
				nmisses_allowed);
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
		  query_compress_rev,minus_sarray,/*plusp*/false,genestrand,first_read_p,minus_conversion);
    elt = Elt_new(minus_querypos,nmatches,initptr,finalptr);
    minus_querypos += nmatches;
    minus_querypos += 1;		/* To skip the presumed mismatch */

    debug(printf("minus_querypos %d vs querylength %d\n",minus_querypos,querylength));
    if (nmatches <= best_minus_nmatches) {
      /* Initial (left) elt was best */
      debug(printf("Initial elt %p was best:\n",best_minus_elt));
      minus_set = List_push(NULL,elt);
      if (minus_querypos >= querylength) {
	chrhigh = 0U;
	Elt_fill_positions_all(best_minus_elt,minus_sarray);
	debug(Elt_dump(best_minus_elt));
	for (i = 0; i < best_minus_elt->npositions; i++) {
	  left = best_minus_elt->positions[i];
	  if (left > chrhigh) {
	    chrnum = Univ_IIT_get_one(chromosome_iit,left,left);
	    Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
	    /* *chrhigh += 1U; */
	  }
	  if (extend_rightward(/*goal*/left,chroffset,chrhigh,/*rightward_set*/minus_set,
			       query_compress_rev,minus_sarray,/*plusp*/false,genestrand,first_read_p,
			       best_minus_elt->queryend) == true) {
	    collect_elt_matches(&(*found_score),&(*subs),&(*indels),&(*ambiguous),&(*singlesplicing),&(*doublesplicing),
				best_minus_elt->querystart,best_minus_elt->queryend,
				chrnum,chroffset,chrhigh,chrlength,
				/*goal*/left,/*rightward_set*/minus_set,/*leftward_set*/NULL,
				querylength,query_compress_rev,/*plusp*/false,genestrand,first_read_p,
				nmisses_allowed);
	  }
	}
      }
    } else {
      /* Second (right) elt is best */
      debug(printf("Second elt %p is best:\n",elt));
      minus_set = List_push(NULL,best_minus_elt);
      best_minus_elt = elt;
      best_minus_nmatches = nmatches;
      if (minus_querypos >= querylength) {
	chrhigh = 0U;
	Elt_fill_positions_all(best_minus_elt,minus_sarray);
	debug(Elt_dump(best_minus_elt));
	for (i = 0; i < best_minus_elt->npositions; i++) {
	  left = best_minus_elt->positions[i];
	  if (left > chrhigh) {
	    chrnum = Univ_IIT_get_one(chromosome_iit,left,left);
	    Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
	    /* *chrhigh += 1U; */
	  }
	  nmatches = Genome_consecutive_matches_leftward(query_compress_rev,left,
							 /*pos5*/0,/*pos3*/best_minus_elt->querystart,
							 /*plusp*/false,genestrand,first_read_p);
	  debug(printf(" extending bestelt querystart %d leftward by %d matches\n",best_minus_elt->querystart,nmatches));
	  best_minus_elt->querystart -= nmatches;
	  if (extend_leftward(/*goal*/left,chroffset,chrhigh,/*leftward_set*/minus_set,
			      /*queryptr*/queryrc,query_compress_rev,
			      minus_sarray,/*plusp*/false,genestrand,first_read_p,minus_conversion,
			      best_minus_elt->querystart,best_minus_elt->queryend) == true) {
	    collect_elt_matches(&(*found_score),&(*subs),&(*indels),&(*ambiguous),&(*singlesplicing),&(*doublesplicing),
				best_minus_elt->querystart,best_minus_elt->queryend,
				chrnum,chroffset,chrhigh,chrlength,
				/*goal*/left,/*rightward_set*/NULL,/*leftward_set*/minus_set,
				querylength,query_compress_rev,/*plusp*/false,genestrand,first_read_p,
				nmisses_allowed);
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
		    query_compress_fwd,plus_sarray,/*plusp*/true,genestrand,first_read_p,plus_conversion);
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
		    query_compress_rev,minus_sarray,/*plusp*/false,genestrand,first_read_p,minus_conversion);
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
		  query_compress_fwd,plus_sarray,/*plusp*/true,genestrand,first_read_p,plus_conversion);
    elt = Elt_new(plus_querypos,nmatches,initptr,finalptr);
    plus_set = List_push(plus_set,(void *) elt);
    if (nmatches > best_plus_nmatches) {
      best_plus_elt = elt;
      best_plus_nmatches = nmatches;

      /* See if we have a substitution winner */
      Elt_fill_positions_all(best_plus_elt,plus_sarray);
      for (i = 0; i < best_plus_elt->npositions; i++) {
	left = best_plus_elt->positions[i];
	/* Should return max_mismatches + 1 if it exceeds the limit */
	if ((nmismatches = Genome_count_mismatches_limit(query_compress_fwd,left,/*pos5*/0,/*pos3*/querylength,
							 /*max_mismatches*/nmisses_allowed,
							 /*plusp*/true,genestrand,first_read_p)) <= nmisses_allowed) {
	  chrnum = Univ_IIT_get_one(chromosome_iit,left,left);
	  Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
	  if ((hit = Stage3end_new_substitution(&(*found_score),nmismatches,
						left,/*genomiclength*/querylength,
						query_compress_fwd,/*plusp*/true,genestrand,first_read_p,
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
		  query_compress_rev,minus_sarray,/*plusp*/false,genestrand,first_read_p,minus_conversion);
    elt = Elt_new(minus_querypos,nmatches,initptr,finalptr);
    minus_set = List_push(minus_set,(void *) elt);
    if (nmatches > best_minus_nmatches) {
      best_minus_elt = elt;
      best_minus_nmatches = nmatches;

      /* See if we have a substitution winner */
      Elt_fill_positions_all(best_minus_elt,minus_sarray);
      for (i = 0; i < best_minus_elt->npositions; i++) {
	left = best_minus_elt->positions[i];
	/* Should return max_mismatches + 1 if it exceeds the limit */
	if ((nmismatches = Genome_count_mismatches_limit(query_compress_rev,left,/*pos5*/0,/*pos3*/querylength,
							 /*max_mismatches*/nmisses_allowed,
							 /*plusp*/false,genestrand,first_read_p)) <= nmisses_allowed) {
	  chrnum = Univ_IIT_get_one(chromosome_iit,left,left);
	  Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
	  if ((hit = Stage3end_new_substitution(&(*found_score),nmismatches,
						left,/*genomiclength*/querylength,
						query_compress_rev,/*plusp*/false,genestrand,first_read_p,
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
      array = (Elt_T *) MALLOCA(nelts * sizeof(Elt_T));
      List_fill_array_and_free((void **) array,&rightward_set);
      rightward_set = (List_T) NULL;
    
      qsort(array,nelts,sizeof(Elt_T),Elt_querypos_ascending_cmp);
      for (i = nelts-1; i >= 0; --i) {
	rightward_set = List_push(rightward_set,(void *) array[i]);
      }
      FREEA(array);
    }

    if ((nelts = List_length(leftward_set)) > 0) {
      array = (Elt_T *) MALLOCA(nelts * sizeof(Elt_T));
      List_fill_array_and_free((void **) array,&leftward_set);
      leftward_set = (List_T) NULL;
    
      qsort(array,nelts,sizeof(Elt_T),Elt_querypos_descending_cmp);
      for (i = nelts-1; i >= 0; --i) {
	leftward_set = List_push(leftward_set,(void *) array[i]);
      }
      FREEA(array);
    }


    chrhigh = 0U;
    Elt_fill_positions_all(best_plus_elt,plus_sarray);
    for (i = 0; i < best_plus_elt->npositions; i++) {
      left = best_plus_elt->positions[i];
      if (left > chrhigh) {
	chrnum = Univ_IIT_get_one(chromosome_iit,left,left);
	Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
	/* *chrhigh += 1U; */
      }
      if (extend_rightward(/*goal*/left,chroffset,chrhigh,rightward_set,
			   query_compress_fwd,plus_sarray,/*plusp*/true,genestrand,first_read_p,
			   best_plus_elt->queryend) == true) {
	nmatches = Genome_consecutive_matches_leftward(query_compress_fwd,left,
						       /*pos5*/0,/*pos3*/best_plus_elt->querystart,
						       /*plusp*/true,genestrand,first_read_p);
	debug(printf(" extending bestelt querystart %d leftward by %d matches\n",best_plus_elt->querystart,nmatches));
	best_plus_elt->querystart -= nmatches;
	if (extend_leftward(/*goal*/left,chroffset,chrhigh,leftward_set,
			    /*queryptr*/queryuc_ptr,query_compress_fwd,
			    plus_sarray,/*plusp*/true,genestrand,first_read_p,plus_conversion,
			    best_plus_elt->querystart,best_plus_elt->queryend) == true) {
	  collect_elt_matches(&(*found_score),&(*subs),&(*indels),&(*ambiguous),&(*singlesplicing),&(*doublesplicing),
			      best_plus_elt->querystart,best_plus_elt->queryend,
			      chrnum,chroffset,chrhigh,chrlength,
			      /*goal*/left,rightward_set,leftward_set,
			      querylength,query_compress_fwd,/*plusp*/true,genestrand,first_read_p,
			      nmisses_allowed);
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
      array = (Elt_T *) MALLOCA(nelts * sizeof(Elt_T));
      List_fill_array_and_free((void **) array,&rightward_set);
      rightward_set = (List_T) NULL;
    
      qsort(array,nelts,sizeof(Elt_T),Elt_querypos_ascending_cmp);
      for (i = nelts-1; i >= 0; --i) {
	rightward_set = List_push(rightward_set,(void *) array[i]);
      }
      FREEA(array);
    }

    if ((nelts = List_length(leftward_set)) > 0) {
      array = (Elt_T *) MALLOCA(nelts * sizeof(Elt_T));
      List_fill_array_and_free((void **) array,&leftward_set);
      leftward_set = (List_T) NULL;
    
      qsort(array,nelts,sizeof(Elt_T),Elt_querypos_descending_cmp);
      for (i = nelts-1; i >= 0; --i) {
	leftward_set = List_push(leftward_set,(void *) array[i]);
      }
      FREEA(array);
    }

    chrhigh = 0U;
    Elt_fill_positions_all(best_minus_elt,minus_sarray);
    for (i = 0; i < best_minus_elt->npositions; i++) {
      left = best_minus_elt->positions[i];
      if (left > chrhigh) {
	chrnum = Univ_IIT_get_one(chromosome_iit,left,left);
	Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,circular_typeint);
	/* *chrhigh += 1U; */
      }
      if (extend_rightward(/*goal*/left,chroffset,chrhigh,rightward_set,
			   query_compress_rev,minus_sarray,/*plusp*/false,genestrand,first_read_p,
			   best_minus_elt->queryend) == true) {
	nmatches = Genome_consecutive_matches_leftward(query_compress_rev,left,
						       /*pos5*/0,/*pos3*/best_minus_elt->querystart,
						       /*plusp*/false,genestrand,first_read_p);
	debug(printf(" extending bestelt querystart %d leftward by %d matches\n",best_minus_elt->querystart,nmatches));
	best_minus_elt->querystart -= nmatches;
	if (extend_leftward(/*goal*/left,chroffset,chrhigh,leftward_set,
			    /*queryptr*/queryrc,query_compress_rev,
			    minus_sarray,/*plusp*/false,genestrand,first_read_p,minus_conversion,
			    best_minus_elt->querystart,best_minus_elt->queryend) == true) {
	  collect_elt_matches(&(*found_score),&(*subs),&(*indels),&(*ambiguous),&(*singlesplicing),&(*doublesplicing),
			      best_minus_elt->querystart,best_minus_elt->queryend,
			      chrnum,chroffset,chrhigh,chrlength,
			      /*goal*/left,rightward_set,leftward_set,
			      querylength,query_compress_rev,/*plusp*/false,genestrand,first_read_p,
			      nmisses_allowed);
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

