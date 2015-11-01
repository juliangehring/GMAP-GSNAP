static char rcsid[] = "$Id: sarray-write.c 122381 2013-12-24 01:22:55Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "sarray-write.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/mman.h>		/* For munmap */
#include "bool.h"
#include "access.h"
#include "mem.h"
#include "genomicpos.h"
#include "assert.h"
#include "compress.h"
#include "bitpack64-write.h"
#include "bitpack64-read.h"
#include "bitpack64-access.h"	/* For Sarray_plcp_compare */
#include "fopen.h"
#include "saca-k.h"
#include "genome_hr.h"
#include "uintlist.h"
#include "intlist.h"
#include "popcount.h"


#ifdef USE_CHILD_BP
#include "bp.h"
#include "bp-write.h"
#include "bp-read.h"		/* For Sarray_child_test() */
#endif


#ifdef WORDS_BIGENDIAN
#include "bigendian.h"
#else
#include "littleendian.h"
#endif


/* make_index */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

/* make_child */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

/* child_bp */
#ifdef DEBUG3
#define debug3(x) x
#else
#define debug3(x)
#endif

/* correctness of using Genome_consecutive_matches_pair */
#ifdef DEBUG14
#define debug14(x) x
#else
#define debug14(x)
#endif

/* correctness of make_index_incremental */
#ifdef DEBUG15
#define debug15(x) x
#else
#define debug15(x)
#endif



/* For computing LCP.  Comment out because mmap is faster than fread */
/* #define READ_SA_FROM_FILE 1 */

#define MONITOR_INTERVAL 100000000 /* 100 million nt */

#define set_bit(bitvector,i) (bitvector)[(i)/32] |= (1 << ((i)%32))
#define get_bit(bitvector,i) ((bitvector)[(i)/32] & (1 << ((i)%32)))


#if 0
static void
compute_lcp (UINT4 *lcp,
#ifdef DEBUG14
	     unsigned char *s,
#endif
	     UINT4 *SA, UINT4 n) {
  UINT4 *rank, h;
  UINT4 i, j;
  char *comma;
#ifdef DEBUG14
  UINT4 horig;
#endif

  rank = (UINT4 *) CALLOC(n+1,sizeof(UINT4));
  for (i = 0; i <= n; i++) {
    rank[SA[i]] = i;
  }

  lcp[0] = 0;			/* -1 ? */
  h = 0;
  for (i = 0; i <= n; i++) {
    if (rank[i] > 0) {
      j = SA[rank[i] - 1];
#ifdef DEBUG14
      horig = h;
      while (i + h < n && j + h < n && s[i+h] == s[j+h]) {
	h++;
      }
      if ((h - horig) != Genome_consecutive_matches_pair(i+horig,j+horig,/*genomelength*/n)) {
	abort();
      }
#else
      h += Genome_consecutive_matches_pair(i+h,j+h,/*genomelength*/n);
#endif
      lcp[rank[i]] = h;
      if (h > 0) {
	h--;
      }
    }
    if (i % MONITOR_INTERVAL == 0) {
      comma = Genomicpos_commafmt(i);
      fprintf(stderr,"Computing lcp index %s\n",comma);
      FREE(comma);
    }
  }

  FREE(rank);

  return;
}
#endif


/* Computes permuted lcp (Karkkainen, CPM 2009) */
/* Two methods for storing plcp as a cumulative sum: (1) Cum value at k is
   sum_{i=1}^k (plcp[i] + 1 - plcp[i-1]), or more simply (2) plcp[k] + k */
static void
compute_plcp (UINT4 *plcp,
#ifdef READ_SA_FROM_FILE
	      FILE *fp,
#else
	      UINT4 *SA,
#endif
	      UINT4 n) {
  UINT4 *phi, h;
  UINT4 i, j;
  char *comma;
#ifdef READ_SA_FROM_FILE
  UINT4 sa, prev_sa;
#endif

  phi = plcp;			/* Use space allocated for plcp */

#ifdef READ_SA_FROM_FILE
  fprintf(stderr,"Inverting suffix array (via file)...");
  FREAD_UINT(&prev_sa,fp);
  for (i = 1; i <= n; i++) {
    FREAD_UINT(&sa,fp);
    phi[sa] = prev_sa;
    prev_sa = sa;
  }
#else
  fprintf(stderr,"Inverting suffix array (via mmap)...");
  for (i = 1; i <= n; i++) {
    phi[SA[i]] = SA[i-1];
  }
#endif
  fprintf(stderr,"done\n");
  /* Note that phi[n] is not assigned, because SA[i] == n for i == 0, and we don't look up i == 0 */


  h = 0;
  for (i = 0; i < n; i++) {
    j = phi[i];			/* To be overwritten by plcp[i] */
    h += Genome_consecutive_matches_pair(i+h,j+h,/*genomelength*/n);
    plcp[i] = h;				 /* overwrites phi[i] */
    if (h > 0) {
      h--;
    }

    if (i % MONITOR_INTERVAL == 0) {
      comma = Genomicpos_commafmt(i);
      fprintf(stderr,"Computing permuted lcp index %s\n",comma);
      FREE(comma);
    }
  }

  /* This makes lcp[0] = -1, because lcp[0] = plcp[SA[0]] = plcp[n] = -1 */
  plcp[n] = -1;

  return;
}



static UINT4
power (int base, int exponent) {
  UINT4 result = 1U;
  int i;

  for (i = 0; i < exponent; i++) {
    result *= base;
  }
  return result;
}


#define LOW_TWO_BITS 0x3
#define RIGHT_A 0
#define RIGHT_C 1
#define RIGHT_G 2
#define RIGHT_T 3


static void
oligo_nt (char *nt, UINT4 oligo, int oligosize) {
  int i, j;
  UINT4 lowbits;

  j = oligosize-1;
  for (i = 0; i < oligosize; i++) {
    lowbits = oligo & LOW_TWO_BITS;
    switch (lowbits) {
    case RIGHT_A: nt[j] = 'A'; break;
    case RIGHT_C: nt[j] = 'C'; break;
    case RIGHT_G: nt[j] = 'G'; break;
    case RIGHT_T: nt[j] = 'T'; break;
    }
    oligo >>= 2;
    j--;
  }

  return;
}


/* Taken from Johannes Fischer, Advanced Text Indexing Techniques, Algorithm 1 */
/* Does not use LCP, so time is O(m * log(n)) */
/* Result should be within [i..j] */
static void
sarray_search_simple (Sarrayptr_T *initptr, Sarrayptr_T *finalptr, char *query,
		      int querylength, Genome_T genomecomp, UINT4 *SA,
		      UINT4 i, UINT4 j, UINT4 n) {
  Sarrayptr_T low, high, mid;
  Univcoord_T pos;
  int nmatches;
  char c;


  low = i;
  high = j+1;

  while (low < high) {
    /* Compute mid for unsigned ints.  Want floor((low+high)/2). */
    mid = low/2 + high/2;
    if (low % 2 == 1 && high % 2 == 1) {
      mid += 1;
    }

    nmatches = 0;
    pos = SA[mid];

    while (nmatches < querylength && (c = Genome_get_char_lex(genomecomp,pos,n)) == query[nmatches]) {
      nmatches++;
      pos++;
    }
    if (nmatches == querylength || c > query[nmatches]) {
      high = mid;
    } else {
      low = mid + 1;
    }
  }

  *initptr = low;

  low--;
  high = j;

  while (low < high) {
    /* Compute mid for unsigned ints.  Want ceil((low+high)/2). */
    mid = low/2 + high/2;
    if (low % 2 == 1 || high % 2 == 1) {
      mid += 1;
    }

    nmatches = 0;
    pos = SA[mid];

    while (nmatches < querylength && (c = Genome_get_char_lex(genomecomp,pos,n)) == query[nmatches]) {
      nmatches++;
      pos++;
    }
    if (nmatches == querylength || c < query[nmatches]) {
      low = mid;
    } else {
      high = mid - 1;
    }
  }

  *finalptr = high;

  return;
}


#if 0
static int
get_all_children (bool *filledp, Sarrayptr_T *l, Sarrayptr_T *r, Sarrayptr_T i, Sarrayptr_T j,
		  UINT4 *child, UINT4 *nextp, Genome_T genomecomp, UINT4 *SA, UINT4 *plcpptrs, UINT4 *plcpcomp,
		  int indexsize, UINT4 n) {
  int noccupied = 0;
  UINT4 up, nextl;
  Sarrayptr_T sa_nextl;
  UINT4 lcp_whole;
  UINT4 pos;
  char c;

  /* Test for child[j] being up: lcp[j] > lcp[j+1] */
  up = child[j];		/* childtab[j+1].up */
  if (i < up && up <= j) {
    nextl = up;
  } else {
    nextl = child[i];	/* down */
  }
  sa_nextl = SA[nextl];
  lcp_whole = Bitpack64_offsetptr_only(sa_nextl,plcpptrs,plcpcomp) - sa_nextl;

  if (lcp_whole != (UINT4) (indexsize - 1)) {
    /* Not at desired level, so exit this procedure */
    return 0;
  } else {
    debug1(printf("Filling children for lcp-interval %u..%u with lcp_whole %d\n",i,j,lcp_whole));
  }

  pos = SA[i] + lcp_whole;
  c = Genome_get_char_lex(genomecomp,pos,n);

  debug1(printf("For char %c, creating interval %u..%u\n",c,i,nextl-1));
  switch (c) {
  case 'A': l[0] = i; r[0] = nextl - 1; filledp[0] = true; noccupied++; break;
  case 'C': l[1] = i; r[1] = nextl - 1; filledp[1] = true; noccupied++; break;
  case 'G': l[2] = i; r[2] = nextl - 1; filledp[2] = true; noccupied++; break;
  case 'T': l[3] = i; r[3] = nextl - 1; filledp[3] = true; noccupied++; break;
  }

  
  /* Test for child[i] being down: lcp[child[i]] > lcp[i] */
  /* Test for child[i] being next_lindex: lcp[child[i]] == lcp[i] */
  while (get_bit(nextp,nextl) != 0) {
    pos = SA[nextl] + lcp_whole;
    c = Genome_get_char_lex(genomecomp,pos,n);

    debug1(printf("For char %c, creating interval %u..%u\n",c,nextl,child[nextl]-1));
    switch (c) {
    case 'A': l[0] = nextl; r[0] = child[nextl] - 1; filledp[0] = true; noccupied++; break;
    case 'C': l[1] = nextl; r[1] = child[nextl] - 1; filledp[1] = true; noccupied++; break;
    case 'G': l[2] = nextl; r[2] = child[nextl] - 1; filledp[2] = true; noccupied++; break;
    case 'T': l[3] = nextl; r[3] = child[nextl] - 1; filledp[3] = true; noccupied++; break;
    }
    
    nextl = child[nextl];
  }

  pos = SA[nextl] + lcp_whole;
  c = Genome_get_char_lex(genomecomp,pos,n);

  debug1(printf("For char %c, creating interval %u..%u\n",c,nextl,j));
  switch (c) {
  case 'A': l[0] = nextl; r[0] = j; filledp[0] = true; noccupied++; break;
  case 'C': l[1] = nextl; r[1] = j; filledp[1] = true; noccupied++; break;
  case 'G': l[2] = nextl; r[2] = j; filledp[2] = true; noccupied++; break;
  case 'T': l[3] = nextl; r[3] = j; filledp[3] = true; noccupied++; break;
  }

  return noccupied;
}
#endif



/* oligo is based on old indexsize. indexsize is new indexsize. */
static int
sarray_search_incremental (Sarrayptr_T *initptrs, Sarrayptr_T *finalptrs,
			   Oligospace_T oligo, char *query, Genome_T genomecomp, UINT4 *SA,
			   Sarrayptr_T *saindexi, Sarrayptr_T *saindexj,
			   int indexsize, Oligospace_T prev_oligospace, UINT4 n) {
  int noccupied = 0, k;
  Sarrayptr_T i, j;
  Oligospace_T prev_oligo;

  prev_oligo = oligo/4;
  if (saindexj[prev_oligo] < saindexi[prev_oligo]) {
    /* Oligo from 0..(indexsize)-1 does not match, so lengthening also will not match */
    initptrs[0] = initptrs[1] = initptrs[2] = initptrs[3] = saindexi[prev_oligo];
    finalptrs[0] = finalptrs[1] = finalptrs[2] = finalptrs[3] = saindexj[prev_oligo];

  } else {
    if (oligo == 0) {
      i = 1;
      j = saindexj[prev_oligo + 1];
    } else if (prev_oligo + 1 == prev_oligospace) {
      i = saindexi[prev_oligo - 1];
      j = n;
    } else {
      i = saindexi[prev_oligo - 1];
      j = saindexj[prev_oligo + 1];
    }

    for (k = 0; k < 4; k++) {
      oligo_nt(query,oligo+k,indexsize);
      sarray_search_simple(&(initptrs[k]),&(finalptrs[k]),query,/*querylength*/indexsize,
			   genomecomp,SA,i,j,n);
      if (initptrs[k] <= finalptrs[k]) {
	noccupied++;
      }
    }
  }

  return noccupied;
}


#define INDEX_MONITOR_INTERVAL 100000


static UINT4
make_index (Sarrayptr_T *saindexi, Sarrayptr_T *saindexj,
	    UINT4 oligospace, int querylength, Genome_T genomecomp, UINT4 *SA, UINT4 n) {
  UINT4 noccupied = 0;
  char *queryuc_ptr;
  UINT4 oligo;
  char *comma;

  queryuc_ptr = (char *) CALLOC(querylength+1,sizeof(char));

  for (oligo = 0; oligo < oligospace; oligo++) {
    oligo_nt(queryuc_ptr,oligo,querylength);
    sarray_search_simple(&(saindexi[oligo]),&(saindexj[oligo]),queryuc_ptr,querylength,
			 genomecomp,SA,/*i*/1,/*j*/n,n);
    if (saindexi[oligo] <= saindexj[oligo]) {
      debug1(printf("%u\t%s\t%u\t%u\n",oligo,queryuc_ptr,saindexi[oligo],saindexj[oligo]));
      noccupied++;
    }
#if 0
    if (oligo % INDEX_MONITOR_INTERVAL == 0) {
      comma = Genomicpos_commafmt(oligo);
      fprintf(stderr,"Computing index %s\n",comma);
      FREE(comma);
    }
#endif
  }
  FREE(queryuc_ptr);

  return noccupied;
}


static UINT4
make_index_incremental (Sarrayptr_T *initptrs, Sarrayptr_T *finalptrs,
			Genome_T genomecomp, UINT4 *SA, Sarrayptr_T *saindexi, Sarrayptr_T *saindexj, 
			int indexsize, Oligospace_T oligospace, Oligospace_T prev_oligospace, UINT4 n) {
  int noccupied = 0;
  char *query;
  Oligospace_T oligo;
  char *comma;

  query = (char *) CALLOC(indexsize+1,sizeof(char));

  for (oligo = 0; oligo < oligospace; oligo += 4) {
    noccupied += sarray_search_incremental(&(initptrs[oligo]),&(finalptrs[oligo]),oligo,query,
					   genomecomp,SA,saindexi,saindexj,indexsize,prev_oligospace,n);
#if 0
    if (oligo % INDEX_MONITOR_INTERVAL == 0) {
      comma = Genomicpos_commafmt(oligo);
      fprintf(stderr,"Computing index %s\n",comma);
      FREE(comma);
    }
#endif
  }

  FREE(query);

  return noccupied;
}


void
Sarray_write_array (char *sarrayfile, Genome_T genomecomp, UINT4 genomelength) {
  UINT4 *SA;
  UINT4 n = genomelength;
  unsigned char *gbuffer;
  FILE *fp;


  SA = (UINT4 *) MALLOC((n+1)*sizeof(UINT4));
  gbuffer = (unsigned char *) CALLOC(n+1,sizeof(unsigned char));
  Genome_fill_buffer_int_string(genomecomp,/*left*/0,/*length*/n,gbuffer);
  gbuffer[n] = 0;		       /* '\0', terminator */
  SACA_K(gbuffer,SA,n+/*virtual sentinel*/1,/*K, alphabet_size*/5,/*m*/n+1,/*level*/0);

  if ((fp = FOPEN_WRITE_BINARY(sarrayfile)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",sarrayfile);
    exit(9);
  } else {
    FWRITE_UINTS(SA,n+1,fp);
    fclose(fp);
  }

  FREE(gbuffer);
  FREE(SA);

  return;
}


/* Uses permuted lcp for speed and reduced memory usage */
void
Sarray_write_plcp (char *plcpptrsfile, char *plcpcompfile, char *sarrayfile,
		   Genome_T genomecomp, UINT4 genomelength) {
  UINT4 *plcp;
  UINT4 *ramp;

  UINT4 n = genomelength, i;
#ifdef READ_SA_FROM_FILE
  FILE *fp;
#else
  UINT4 *SA;
  int fd;
  size_t len;
#endif

  plcp = (UINT4 *) MALLOC((n+1)*sizeof(UINT4));
  ramp = plcp;

#ifdef READ_SA_FROM_FILE
  if ((fp = FOPEN_READ_BINARY(sarrayfile)) == NULL) {
    fprintf(stderr,"Can't read from file %s\n",sarrayfile);
    exit(9);
  }
  compute_plcp(plcp,fp,n);
#else
  SA = (UINT4 *) Access_mmap(&fd,&len,sarrayfile,sizeof(UINT4),/*randomp*/false);
  compute_plcp(plcp,SA,n);
#endif

#ifdef USE_CUMDEV
  /* Compute deviations from downward ramp */
  for (i = n; i >= 1; i--) {
    dev[i] = (plcp[i] + 1) - plcp[i-1];
  }
  /* dev[0] = plcp[0]; */
#else
  for (i = 0; i <= n; i++) {
    ramp[i] = plcp[i] + i;
  }
#endif

  fprintf(stderr,"Writing permuted lcp file...");
  /* Provide n to write values [0..n] */
  Bitpack64_write_differential(plcpptrsfile,plcpcompfile,ramp,n);
  fprintf(stderr,"done\n");

  FREE(plcp);

#ifdef READ_SA_FROM_FILE
  fclose(fp);
#else
  munmap((void *) SA,len);
  close(fd);
#endif
  
  return;
}




/* Without encoding, would be child[index] = value */
#define encode_up(child,index,value) child[index-1] = (index - 1) - value
#define encode_down(child,index,value) child[index] = value - 1 - index
#define encode_next(child,index,value) child[index] = value - 1 - index

/* For child[index+1].up, just calling child[index] */
#define decode_up(child_i,index) index - child_i
#define decode_down(child_i,index) child_i + index + 1
#define decode_next(child_i,index) child_i + index + 1


#ifdef USE_CHILD_BP
static bp *
make_child_bp (bool **first_child_p, UINT4 *SA, UINT4 *plcpptrs, UINT4 *plcpcomp, UINT4 n) {
  bp *childbp;
  UINT4 i, j = 0;
  BP_size_t k = 0, bplength;
  Uintlist_T lcpstack, right_sibling_stack;
  UINT4 sa_i, lcp_i, lcp_lastindex;
  UINT4 right_sibling_p;
  char *comma;

  bplength = 2 * (BP_size_t) n;
  fprintf(stderr,"Allocating %lu bytes for child bp\n",bplength*sizeof(bp));
  childbp = (bp *) MALLOC(bplength * sizeof(bp));
  *first_child_p = (bool *) MALLOC((n+3) * sizeof(bool));

  *first_child_p[j++] = true;	/* Initial value.  Helps rank work properly. */

  i = 1;
  debug3(printf("("));
  childbp[k++] = open_paren;
      

  lcp_i = 0;			/* Just as good as -1 and avoids comparison with -1U */
  lcpstack = Uintlist_push(NULL,lcp_i);
  right_sibling_stack = Uintlist_push(NULL,(UINT4) false);

  for (i = 2; i <= n; i++) {
    sa_i = SA[i];
    lcp_i = Bitpack64_offsetptr_only(sa_i,plcpptrs,plcpcomp) - sa_i;

    while (lcp_i < Uintlist_head(lcpstack)) {
      lcpstack = Uintlist_pop(lcpstack,&lcp_lastindex);
      right_sibling_stack = Uintlist_pop(right_sibling_stack,&right_sibling_p);
      if (right_sibling_p == true) {
	(*first_child_p)[j++] = false;
	debug3(printf("]"));
      } else {
	(*first_child_p)[j++] = true;
	debug3(printf(")"));
      }
      childbp[k++] = close_paren;
    }

    if (lcp_i == Uintlist_head(lcpstack)) {
      /* This is a right sibling */
      debug3(printf("["));
      childbp[k++] = open_paren;
      lcpstack = Uintlist_push(lcpstack,lcp_i);
      right_sibling_stack = Uintlist_push(right_sibling_stack,(UINT4) true);
    } else {
      /* This is a right child */
      debug3(printf("("));
      childbp[k++] = open_paren;
      lcpstack = Uintlist_push(lcpstack,lcp_i);
      right_sibling_stack = Uintlist_push(right_sibling_stack,(UINT4) false);
    }

    if (i % MONITOR_INTERVAL == 0) {
      comma = Genomicpos_commafmt(i);
      fprintf(stderr,"Computing child bp %s\n",comma);
      FREE(comma);
    }
  }

  /* i == n + 1.  lcp[n+1] = -1 */
  while (lcpstack != NULL) {
    lcpstack = Uintlist_pop(lcpstack,&lcp_lastindex);
    right_sibling_stack = Uintlist_pop(right_sibling_stack,&right_sibling_p);
    if (right_sibling_p == true) {
      (*first_child_p)[j++] = false;
      debug3(printf("]"));
    } else {
      (*first_child_p)[j++] = true;
      debug3(printf(")"));
    }
    childbp[k++] = close_paren;
  }
      
  /* rank0 = rank; */
  /* subrank_stack = (Intlist_T) NULL; */

  /* close of node n, which is the right sibling of node 0 */
  (*first_child_p)[j++] = false;
  /* close of node 0 */
  (*first_child_p)[j++] = true;

  assert(k == bplength);

  return childbp;
}

#else

static UINT4 *
make_child_twopass (UINT4 **nextp, UINT4 *nbytes, UINT4 *SA, UINT4 *plcpptrs, UINT4 *plcpcomp, UINT4 n) {
  UINT4 *child;
  UINT4 lastindex, i;
  Uintlist_T lcpstack, indexstack;
  UINT4 sa_i, lcp_i, lcp_lastindex;
  char *comma;

  *nbytes = ((n+1) + 31)/32;
  *nextp = (UINT4 *) CALLOC(*nbytes,sizeof(UINT4));

#if 0
  child = (UINT4 *) MALLOC((n+1) * sizeof(UINT4));
  for (i = 0; i <= n; i++) {
    child[i] = -1U;
  }
#else
  child = (UINT4 *) CALLOC(n+1,sizeof(UINT4));
#endif

  
  /* Because we sort suffixes with $ < rest of alphabet, we never use
     the entry at 0, where SA[0] = n and lcp[0] = 0 */

  /* Algorithm 6.2: Compute up and down values */
  lastindex = 0;

  fprintf(stderr,"Computing child up/down index 0\n");
  i = 1;
  indexstack = Uintlist_push(NULL,i);

  sa_i = SA[i];
  lcp_i = Bitpack64_offsetptr_only(sa_i,plcpptrs,plcpcomp) - sa_i;
  lcpstack = Uintlist_push(NULL,lcp_i);

  for (i = 2; i <= n; i++) {
    sa_i = SA[i];
    lcp_i = Bitpack64_offsetptr_only(sa_i,plcpptrs,plcpcomp) - sa_i;
    while (lcp_i < Uintlist_head(lcpstack)) {
      lcpstack = Uintlist_pop(lcpstack,&lcp_lastindex);
      indexstack = Uintlist_pop(indexstack,&lastindex);
      if (lcp_i <= Uintlist_head(lcpstack) && Uintlist_head(lcpstack) != lcp_lastindex) {
#ifdef NO_ENCODING
	child[Uintlist_head(indexstack)] = lastindex; /* down */
#else
	encode_down(child,Uintlist_head(indexstack),lastindex);
#endif
      }
    }
    /* Now lcp[i] >= lcp[stack->first] holds */
    if (lastindex != 0) {
#ifdef NO_ENCODING
      child[i-1] = lastindex;	/* up */
#else
      encode_up(child,i,lastindex);
#endif
      lastindex = 0;
    }
    indexstack = Uintlist_push(indexstack,i);
    lcpstack = Uintlist_push(lcpstack,lcp_i);

    if (i % MONITOR_INTERVAL == 0) {
      comma = Genomicpos_commafmt(i);
      fprintf(stderr,"Computing child up/down index %s\n",comma);
      FREE(comma);
    }
  }

  /* Handle end of suffix array */
  lcp_i = 0;
  while (lcp_i < Uintlist_head(lcpstack)) {
    lcpstack = Uintlist_pop(lcpstack,&lcp_lastindex);
    indexstack = Uintlist_pop(indexstack,&lastindex);
    if (lcp_i <= Uintlist_head(lcpstack) && Uintlist_head(lcpstack) != lcp_lastindex) {
#ifdef NO_ENCODING
      child[Uintlist_head(indexstack)] = lastindex; /* down */
#else
      encode_down(child,Uintlist_head(indexstack),lastindex);
#endif
    }
  }
  if (lastindex != 0) {
#ifdef NO_ENCODING
    child[i-1] = lastindex;	/* up */
#else
    encode_up(child,i,lastindex);
#endif
    /* lastindex = 0; */
  }

  Uintlist_free(&lcpstack);
  Uintlist_free(&indexstack);


  /* Algorithm 6.5: Compute next l-index values */
  fprintf(stderr,"Computing child next index 0\n");
  i = 1;
  indexstack = Uintlist_push(NULL,i);

  sa_i = SA[i];
  lcp_i = Bitpack64_offsetptr_only(sa_i,plcpptrs,plcpcomp) - sa_i;
  lcpstack = Uintlist_push(NULL,lcp_i);

  for (i = 2; i <= n; i++) {
    sa_i = SA[i];
    lcp_i = Bitpack64_offsetptr_only(sa_i,plcpptrs,plcpcomp) - sa_i;
    while (lcp_i < Uintlist_head(lcpstack)) {
      lcpstack = Uintlist_pop(lcpstack,&lcp_lastindex);
      indexstack = Uintlist_pop(indexstack,&lastindex);
    }
    if (lcp_i == Uintlist_head(lcpstack)) {
      lcpstack = Uintlist_pop(lcpstack,&lcp_lastindex);
      indexstack = Uintlist_pop(indexstack,&lastindex);
#ifdef NO_ENCODING
      child[lastindex] = i;
#else
      encode_next(child,lastindex,i);
#endif
      set_bit(*nextp,lastindex);
    }
    indexstack = Uintlist_push(indexstack,i);
    lcpstack = Uintlist_push(lcpstack,lcp_i);

    if (i % MONITOR_INTERVAL == 0) {
      comma = Genomicpos_commafmt(i);
      fprintf(stderr,"Computing child next index %s\n",comma);
      FREE(comma);
    }
  }

  Uintlist_free(&lcpstack);
  Uintlist_free(&indexstack);

  return child;
}


#if 0
/* Appears to be buggy */
static UINT4 *
make_child_onepass_buggy (UINT4 **nextp, UINT4 *nbytes, UINT4 *SA, UINT4 *plcpptrs, UINT4 *plcpcomp, UINT4 n) {
  UINT4 *child;
  UINT4 lastindex, i;
  Uintlist_T lcpstack = NULL, indexstack = NULL;
  UINT4 sa_i, lcp_i, lcp_lastindex;
  char *comma;

  *nbytes = ((n+1) + 31)/32;
  *nextp = (UINT4 *) CALLOC(*nbytes,sizeof(UINT4));

#if 0
  child = (UINT4 *) MALLOC((n+1) * sizeof(UINT4));
  for (i = 0; i <= n; i++) {
    child[i] = -1U;
  }
#else
  child = (UINT4 *) CALLOC(n+1,sizeof(UINT4));
#endif

  
  /* Because we sort suffixes with $ < rest of alphabet, we never use
     the entry at 0, where SA[0] = n and lcp[0] = -1 */

  fprintf(stderr,"Computing child index 0\n");
  lastindex = -1U;
  i = 1;
  indexstack = Uintlist_push(NULL,i);

  sa_i = SA[i];
  lcp_i = Bitpack64_offsetptr_only(sa_i,plcpptrs,plcpcomp) - sa_i;
  lcpstack = Uintlist_push(NULL,lcp_i);

  for (i = 2; i <= n; i++) {
    sa_i = SA[i];
    lcp_i = Bitpack64_offsetptr_only(sa_i,plcpptrs,plcpcomp) - sa_i;
    while (lcp_i < Uintlist_head(lcpstack)) {
      lcpstack = Uintlist_pop(lcpstack,&lcp_lastindex);
      indexstack = Uintlist_pop(indexstack,&lastindex);
      debug2(printf("Popping index %u, lcp %u\n",lastindex,lcp_lastindex));
      if (lcp_i < Uintlist_head(lcpstack) && Uintlist_head(lcpstack) != lcp_lastindex) {
	debug2(printf("Storing down: child[%u] = %u\n",Uintlist_head(indexstack),lastindex));
	encode_down(child,Uintlist_head(indexstack),lastindex);
      }
    }

    if (lastindex != -1U) {
      debug2(printf("Storing up: child[%u] = %u\n",i-1,lastindex));
      encode_up(child,i,lastindex);
      lastindex = -1U;
    }
    if (lcp_i == Uintlist_head(lcpstack)) {
      debug2(printf("Storing next: child[%u] = %u\n",Uintlist_head(indexstack),i));
      encode_next(child,Uintlist_head(indexstack),i);
      set_bit(*nextp,Uintlist_head(indexstack));
    }
    debug2(printf("Pushing index %u, lcp %u\n",i,lcp_i));
    indexstack = Uintlist_push(indexstack,i);
    lcpstack = Uintlist_push(lcpstack,lcp_i);

    if (i % MONITOR_INTERVAL == 0) {
      comma = Genomicpos_commafmt(i);
      fprintf(stderr,"Computing child index %s\n",comma);
      FREE(comma);
    }
  }

  /* Handle end of suffix array */
  lcp_i = 0;
  while (lcp_i < Uintlist_head(lcpstack)) {
    lcpstack = Uintlist_pop(lcpstack,&lcp_lastindex);
    indexstack = Uintlist_pop(indexstack,&lastindex);
    debug2(printf("Popping index %u, lcp %u (final)\n",lastindex,lcp_lastindex));
    if (lcp_i < Uintlist_head(lcpstack) && Uintlist_head(lcpstack) != lcp_lastindex) {
      debug2(printf("Storing down: child[%u] = %u (final)\n",Uintlist_head(indexstack),lastindex));
      encode_down(child,Uintlist_head(indexstack),lastindex);
    }
  }
  
  if (lastindex != -1U) {
    debug2(printf("Storing up: child[%u] = %u (final)\n",i-1,lastindex));
    encode_up(child,i,lastindex);
    /* lastindex = -1U; */
  }
  if (lcp_i == Uintlist_head(lcpstack)) {
    debug2(printf("Storing next: child[%u] = %u (final)\n",Uintlist_head(indexstack),i));
    encode_next(child,Uintlist_head(indexstack),i);
    set_bit(*nextp,Uintlist_head(indexstack));
  }


  Uintlist_free(&lcpstack);
  Uintlist_free(&indexstack);

  return child;
}
#endif

/* Modeled after make_child_bp */
static UINT4 *
make_child_onepass (UINT4 **nextp, UINT4 *nbytes, UINT4 *SA, UINT4 *plcpptrs, UINT4 *plcpcomp, UINT4 n) {
  UINT4 *child;
  UINT4 lastindex, i;
  Uintlist_T lcpstack, indexstack;
  UINT4 sa_i, lcp_i, lcp_lastindex;
  char *comma;

  *nbytes = ((n+1) + 31)/32;
  *nextp = (UINT4 *) CALLOC(*nbytes,sizeof(UINT4));
  child = (UINT4 *) CALLOC(n+1,sizeof(UINT4));

  i = 1;
  lcp_i = 0;			/* Just as good as -1 and avoids comparison with -1U */
  indexstack = Uintlist_push(NULL,i);
  lcpstack = Uintlist_push(NULL,lcp_i);

  for (i = 2; i <= n; i++) {
    sa_i = SA[i];
    lcp_i = Bitpack64_offsetptr_only(sa_i,plcpptrs,plcpcomp) - sa_i;

    while (lcp_i < Uintlist_head(lcpstack)) {
      indexstack = Uintlist_pop(indexstack,&lastindex);
      lcpstack = Uintlist_pop(lcpstack,&lcp_lastindex);
      /* Mark as a right child */
#ifdef NO_ENCODING
      child[Uintlist_head(indexstack)] = lastindex;
#else
      encode_down(child,Uintlist_head(indexstack),lastindex);
#endif
    }
#ifdef NO_ENCODING
    child[i-1] = lastindex;	/* up */
#else
    encode_up(child,i,lastindex);
#endif

    if (lcp_i == Uintlist_head(lcpstack)) {
      /* This is a right sibling, so mark previous index as having a right sibling */
#ifdef NO_ENCODING
      child[Uintlist_head(indexstack)] = i;
#else
      encode_next(child,Uintlist_head(indexstack),i);
#endif
      set_bit(*nextp,Uintlist_head(indexstack));
    }
    indexstack = Uintlist_push(indexstack,i);
    lcpstack = Uintlist_push(lcpstack,lcp_i);

    if (i % MONITOR_INTERVAL == 0) {
      comma = Genomicpos_commafmt(i);
      fprintf(stderr,"Computing child index %s\n",comma);
      FREE(comma);
    }
  }

  /* i == n + 1.  lcp[n+1] = -1 */
  lcp_i = 0;
  while (lcp_i < Uintlist_head(lcpstack)) {
    indexstack = Uintlist_pop(indexstack,&lastindex);
    lcpstack = Uintlist_pop(lcpstack,&lcp_lastindex);
    /* Mark as a right child */
#ifdef NO_ENCODING
    child[Uintlist_head(indexstack)] = lastindex; /* down */
#else
    encode_down(child,Uintlist_head(indexstack),lastindex);
#endif
  }

  /* Final value for up needs to fail: 1 < up && up <= n, so can use
     0, 1, or n+1.  The value 1 makes sense, since the left pointer
     from node n+1 goes to 1. */
#ifdef NO_ENCODING
  child[/*i-1*/n] = 1;	/* up: start with node 1 */
#else
  encode_up(child,/*i*/n+1,1);
#endif

      
  Uintlist_free(&lcpstack);
  Uintlist_free(&indexstack);

  return child;
}

#endif


void
Sarray_write_child (
#ifdef USE_CHILD_BP
		    char *childbpfile, char *childfcfile, char *childs_pagesfile, char *childs_ptrsfile, char *childs_compfile,
		    char *childr_ptrsfile,  char *childr_compfile, char *childx_ptrsfile, char *childx_compfile,
		    char *pioneerbpfile, char *pior_ptrsfile, char *pior_compfile,
		    char *piom_ptrsfile, char *piom_compfile,
#else
		    char *childptrsfile, char *childcompfile, char *sanextpfile,
#endif
		    char *sarrayfile, char *plcpptrsfile, char *plcpcompfile, UINT4 genomelength) {
  UINT4 n = genomelength, nbytes;
  UINT4 *SA, *plcpptrs, *plcpcomp;
  FILE *fp;
  int sa_fd, plcpcomp_fd;
  size_t sa_len, plcpptrs_len, plcpcomp_len;
  double seconds;

#ifdef USE_CHILD_BP
  BP_size_t bplength;

  UINT4 nblocks;

  bp *childbp;
  bool *first_child_p;
  UINT4 *child_blocks, *first_child_blocks;
  UINT8 *child_select;
  UINT4 *block_excess;

  UINT4 *cumrank;

  bool *pioneerbp;
  UINT4 *pioneer_blocks;
  UINT4 *pioneer_matches;
  UINT4 n_pioneer_pairs;

  int i;
#else
  UINT4 *child, *nextp;
#endif


  SA = (UINT4 *) Access_mmap(&sa_fd,&sa_len,sarrayfile,sizeof(UINT4),/*randomp*/true);

  plcpptrs = (UINT4 *) Access_allocated(&plcpptrs_len,&seconds,plcpptrsfile,sizeof(UINT4));
  plcpcomp = (UINT4 *) Access_mmap(&plcpcomp_fd,&plcpcomp_len,plcpcompfile,sizeof(UINT4),
				  /*randomp*/true);

  Bitpack64_read_setup();

#ifdef USE_CHILD_BP
  bplength = 2 * (BP_size_t) n;

  /* Perhaps we can create and then read blocks instead of childbp */
  childbp = make_child_bp(&first_child_p,SA,plcpptrs,plcpcomp,n); /* childbp valid from [0..2*n-1] */

  fprintf(stderr,"Creating child fc blocks\n");
  first_child_blocks = BP_make_blocks(first_child_p,&nblocks,n+2,BLOCKSIZE);   /* valid from [0..nblocks-1] */
  if ((fp = FOPEN_WRITE_BINARY(childfcfile)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",childfcfile);
    exit(9);
  } else {
    FWRITE_UINTS(first_child_blocks,nblocks,fp);
    fclose(fp);
  }
  FREE(first_child_blocks);
  FREE(first_child_p);


  fprintf(stderr,"Creating child bp select\n");
  child_select = BP_make_select_sampled(&nblocks,childbp,bplength,/*n_ones*/n,SELECT_SAMPLING_INTERVAL);  /* valid from [0..nblocks] */
  Bitpack64_write_differential_huge(childs_pagesfile,childs_ptrsfile,childs_compfile,child_select,nblocks); /* writes values [0..nblocks] */
  FREE(child_select);

  fprintf(stderr,"Creating child bp blocks\n");
  child_blocks = BP_make_blocks(childbp,&nblocks,bplength,BLOCKSIZE);   /* valid from [0..nblocks-1] */
  if ((fp = FOPEN_WRITE_BINARY(childbpfile)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",childbpfile);
    exit(9);
  } else {
    FWRITE_UINTS(child_blocks,nblocks,fp);
    fclose(fp);
  }
  
  fprintf(stderr,"Creating child bp excess\n");
  block_excess = BP_compute_block_excess(&nblocks,childbp,bplength,BLOCKSIZE); /* valid from [0..nblocks-1] */
  Bitpack64_write_direct(childx_ptrsfile,childx_compfile,block_excess,nblocks); /* writes [0..nblocks-1] */
  FREE(block_excess);

  fprintf(stderr,"Creating child bp rank\n");
  cumrank = BP_make_cumrank(child_blocks,nblocks); /* valid from [0..nblocks] */
  FREE(child_blocks);

  Bitpack64_write_differential(childr_ptrsfile,childr_compfile,cumrank,nblocks); /* writes values [0..nblocks] */
  FREE(cumrank);


  /* If pseudop is true, this step modifies childbp so BP_compute_matches works */
  fprintf(stderr,"Computing pioneer bp\n");
  pioneerbp = BP_compute_pioneers(&n_pioneer_pairs,childbp,bplength,BLOCKSIZE,/*pseudop*/true);  /* valid from [0..2*npairs-1] */

  fprintf(stderr,"Creating pioneer blocks\n");
  pioneer_blocks = BP_make_blocks(pioneerbp,&nblocks,bplength,BLOCKSIZE);
  if ((fp = FOPEN_WRITE_BINARY(pioneerbpfile)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",pioneerbpfile);
    exit(9);
  } else {
    FWRITE_UINTS(pioneer_blocks,nblocks,fp);
    fclose(fp);
  }

  fprintf(stderr,"Creating pioneer bp rank\n");
  cumrank = BP_make_cumrank(pioneer_blocks,nblocks); /* valid from [0..nblocks] */
  FREE(pioneer_blocks);

  Bitpack64_write_differential(pior_ptrsfile,pior_compfile,cumrank,nblocks); /* writes values [0..nblocks] */
  FREE(cumrank);


  fprintf(stderr,"Creating pioneer matches\n");
  pioneer_matches = BP_compute_matches(/*markp*/pioneerbp,/*parenbp*/childbp,bplength,
				       /*n_ones*/n_pioneer_pairs,BLOCKSIZE);
  FREE(pioneerbp);
  FREE(childbp);

  Bitpack64_write_direct(piom_ptrsfile,piom_compfile,pioneer_matches,2*n_pioneer_pairs+1); /* writes 0..2*n_pioneer_pairs */
  FREE(pioneer_matches);

#else

  /* child = make_child_twopass(&nextp,&nbytes,SA,plcpptrs,plcpcomp,n); */
  child = make_child_onepass(&nextp,&nbytes,SA,plcpptrs,plcpcomp,n);

  fprintf(stderr,"Writing child nextp file...");
  if ((fp = FOPEN_WRITE_BINARY(sanextpfile)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",sanextpfile);
    exit(9);
  } else {
    FWRITE_UINTS(nextp,nbytes,fp);
    fclose(fp);
  }
  fprintf(stderr,"done\n");

  FREE(nextp);

  fprintf(stderr,"Writing child array file...");
  /* Provide n+1 to write values [0..n] */
  Bitpack64_write_direct(childptrsfile,childcompfile,child,n+1);
  fprintf(stderr,"done\n");

  FREE(child);
#endif

  FREE(plcpptrs);
  munmap((void *) plcpcomp,plcpcomp_len);
  close(plcpcomp_fd);
  munmap((void *) SA,sa_len);
  close(sa_fd);


  return;
}


static void
sarray_search_char (Sarrayptr_T *initptr, Sarrayptr_T *finalptr, char desired_char,
		    Genome_T genomecomp, UINT4 *SA, int n) {
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
    c = Genome_get_char_lex(genomecomp,pos,n);
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
    c = Genome_get_char_lex(genomecomp,pos,n);
    if (desired_char >= c) {
      low = mid;
    } else {
      high = mid - 1;
    }
  }

  *finalptr = high;
  return;
}


#define MIN_INDEXSIZE 12
#define MAX_INDEXSIZE 12

void
Sarray_write_index (char *indexiptrsfile, char *indexicompfile, char *indexjptrsfile,char *indexjcompfile,
		    char *sarrayfile, Genome_T genomecomp, UINT4 genomelength) {
  UINT4 n = genomelength;
  Oligospace_T oligospace, prev_oligospace, noccupied, prev_noccupied;
  Sarrayptr_T *saindexi_new, *saindexj_new, *saindexi_old, *saindexj_old;
  Sarrayptr_T initindexi[4], initindexj[4];
  UINT4 *SA;
  int sa_fd;
  size_t sa_len;
  int indexsize;


  SA = (UINT4 *) Access_mmap(&sa_fd,&sa_len,sarrayfile,sizeof(UINT4),/*randomp*/true);

  sarray_search_char(&(initindexi[0]),&(initindexj[0]),/*desired_char*/'A',genomecomp,SA,n);
  sarray_search_char(&(initindexi[1]),&(initindexj[1]),/*desired_char*/'C',genomecomp,SA,n);
  sarray_search_char(&(initindexi[2]),&(initindexj[2]),/*desired_char*/'G',genomecomp,SA,n);
  sarray_search_char(&(initindexi[3]),&(initindexj[3]),/*desired_char*/'T',genomecomp,SA,n);


  indexsize = MIN_INDEXSIZE;
  oligospace = power(4,/*querylength*/indexsize);
  saindexi_old = (Sarrayptr_T *) CALLOC(oligospace,sizeof(Sarrayptr_T));
  saindexj_old = (Sarrayptr_T *) CALLOC(oligospace,sizeof(Sarrayptr_T));
  prev_noccupied = 0;
  noccupied = make_index(saindexi_old,saindexj_old,oligospace,
			 /*querylength*/indexsize,genomecomp,SA,n);
  fprintf(stderr,"For indexsize %d, occupied %u/%u\n",indexsize,noccupied,oligospace);

  indexsize++;
  while (indexsize <= MAX_INDEXSIZE && noccupied > prev_noccupied) {
    prev_noccupied = noccupied;
    prev_oligospace = oligospace;
    oligospace = power(4,/*querylength*/indexsize);
    saindexi_new = (Sarrayptr_T *) CALLOC(oligospace,sizeof(Sarrayptr_T));
    saindexj_new = (Sarrayptr_T *) CALLOC(oligospace,sizeof(Sarrayptr_T));
    noccupied = make_index_incremental(saindexi_new,saindexj_new,genomecomp,SA,
				       saindexi_old,saindexj_old,indexsize,
				       oligospace,prev_oligospace,n);
    fprintf(stderr,"For indexsize %d, occupied %u/%u\n",indexsize,noccupied,oligospace);
    if (noccupied > prev_noccupied) {
      FREE(saindexj_old);
      FREE(saindexi_old);
      saindexi_old = saindexi_new;
      saindexj_old = saindexj_new;
      indexsize++;
    } else {
      FREE(saindexj_new);
      FREE(saindexi_new);
    }
  }

  indexsize--;
  oligospace = power(4,/*querylength*/indexsize);
  fprintf(stderr,"Optimal indexsize = %d\n",indexsize);

#ifdef DEBUG15
  /* For comparison */
  fprintf(stderr,"Checking...");
  query = (char *) CALLOC(indexsize+1,sizeof(char));
  saindexi_new = (Sarrayptr_T *) CALLOC(oligospace,sizeof(Sarrayptr_T));
  saindexj_new = (Sarrayptr_T *) CALLOC(oligospace,sizeof(Sarrayptr_T));
  noccupied = make_index(saindexi_new,saindexj_new,oligospace,
			 /*querylength*/indexsize,genomecomp,SA,n);
  
  for (oligo = 0; oligo < oligospace; oligo++) {
    if (saindexi_old[oligo] != saindexi_new[oligo]) {
      oligo_nt(query,oligo,indexsize);
      printf("%u\t%s\t%u\t%u\t%u\t%u\n",oligo,query,saindexi_old[oligo],saindexj_old[oligo],saindexi_new[oligo],saindexj_new[oligo]);
      abort();
    } else if (saindexj_old[oligo] != saindexj_new[oligo]) {
      oligo_nt(query,oligo,indexsize);
      printf("%u\t%s\t%u\t%u\t%u\t%u\n",oligo,query,saindexi_old[oligo],saindexj_old[oligo],saindexi_new[oligo],saindexj_new[oligo]);
      abort();
    }
  }
  FREE(query);
  FREE(saindexj_new);
  FREE(saindexi_new);
  fprintf(stderr,"done\n");
#endif


  Bitpack64_write_differential(indexiptrsfile,indexicompfile,saindexi_old,oligospace-1);
  Bitpack64_write_differential(indexjptrsfile,indexjcompfile,saindexj_old,oligospace-1);

  FREE(saindexj_old);
  FREE(saindexi_old);

  munmap((void *) SA,sa_len);
  close(sa_fd);

  return;
}


void
Sarray_array_uncompress (Genome_T genomecomp, char *sarrayfile, char *plcpptrsfile, char *plcpcompfile,
			 UINT4 genomelength, UINT4 start, UINT4 end) {
  UINT4 n = genomelength, pos, match, h;
  unsigned char *gbuffer;

  UINT4 *SA, *plcpptrs, *plcpcomp;

  int sa_fd, plcpcomp_fd;
  size_t sa_len, plcpptrs_len, plcpcomp_len;

  double seconds;
  UINT4 sa_i, lcp_i, sa_nexti, lcp_nexti;


  gbuffer = (unsigned char *) CALLOC(n+1,sizeof(unsigned char));
  Genome_fill_buffer_simple(genomecomp,/*left*/0,/*length*/n,gbuffer);
  gbuffer[n] = 0;		       /* '\0', terminator */

  if (end == 0) {
    start = 0;
    end = n;
  }

  SA = (UINT4 *) Access_mmap(&sa_fd,&sa_len,sarrayfile,sizeof(UINT4),/*randomp*/false);
  plcpptrs = (UINT4 *) Access_allocated(&plcpptrs_len,&seconds,plcpptrsfile,sizeof(UINT4));
  plcpcomp = (UINT4 *) Access_mmap(&plcpcomp_fd,&plcpcomp_len,plcpcompfile,sizeof(UINT4),
				  /*randomp*/true);
  plcpcomp = (UINT4 *) Access_mmap(&plcpcomp_fd,&plcpcomp_len,plcpcompfile,sizeof(UINT4),
				  /*randomp*/true);

  Bitpack64_read_setup();


  printf("i\tSA\tLCP\n");

  pos = start;
  sa_i = SA[pos];
  lcp_i = Bitpack64_offsetptr_only(sa_i,plcpptrs,plcpcomp) - sa_i;

  sa_nexti = SA[pos+1];
  lcp_nexti = Bitpack64_offsetptr_only(sa_nexti,plcpptrs,plcpcomp) - sa_nexti;

  if (pos == 0) {
    /* lcp_i is -1 */
    printf("%u\t%u\t-1\t",pos,sa_i);
    printf("R");
    printf("\n");
  }

  for (pos = start + 1; pos < end; pos++) {
    sa_i = sa_nexti;
    lcp_i = lcp_nexti;

    sa_nexti = SA[pos+1];
    lcp_nexti = Bitpack64_offsetptr_only(sa_nexti,plcpptrs,plcpcomp) - sa_nexti;

    printf("%u\t%u\t%u\t",pos,sa_i,lcp_i);
    if (lcp_i > lcp_nexti) {
      printf("L");
    } else {
      printf("R");
    }

    printf("\t%c",gbuffer[sa_i]);
    for (h = 1; (h <= lcp_i || h <= lcp_nexti) && h < 25; h++) {
      if (gbuffer[sa_i+h] == '\0') {
	printf("$");
      } else {
	printf("%c",gbuffer[sa_i+h]);
      }
    }
    if (h <= lcp_i || h <= lcp_nexti) {
      printf("...");
    }
    
    printf("\n");
  }

  sa_i = sa_nexti;
  lcp_i = lcp_nexti;
  if (pos == n) {
    printf("%u\t%u\t%u\t",pos,sa_i,lcp_i);
    printf("L");

    printf("\t%c",gbuffer[sa_i]);
    for (h = 1; h <= lcp_i && h < 25; h++) {
      if (gbuffer[sa_i+h] == '\0') {
	printf("$");
      } else {
	printf("%c",gbuffer[sa_i+h]);
      }
    }
    if (h <= lcp_i) {
      printf("...");
    }

    printf("\n");
  }

  FREE(plcpptrs);
  munmap((void *) plcpcomp,plcpcomp_len);
  close(plcpcomp_fd);
  munmap((void *) SA,sa_len);
  close(sa_fd);

  return;
}


#ifdef USE_CHILD_BP
void
Sarray_child_uncompress (char *childbpfile, char *childfcfile, char *childs_pagesfile, char *childs_ptrsfile, char *childs_compfile,
			 char *childr_ptrsfile, char *childr_compfile, char *childx_ptrsfile, char *childx_compfile,
			 char *pioneerbpfile, char *pior_ptrsfile, char *pior_compfile,
			 char *piom_ptrsfile, char *piom_compfile, UINT4 genomelength, BP_size_t startl, BP_size_t endl) {
  UINT4 n = genomelength, ranki;
  BP_size_t bplength, selecti, match;

  UINT4 *childbp, *childfc, *pioneerbp;
  UINT4 *childs_pages, *childs_ptrs, *childs_comp, *childr_ptrs, *childr_comp, *childx_ptrs, *childx_comp;
  UINT4 *pior_ptrs, *pior_comp, *piom_ptrs, *piom_comp;

  int childbp_fd, childfc_fd, childs_comp_fd, childr_comp_fd, childx_comp_fd;
  int pioneerbp_fd, pior_comp_fd, piom_comp_fd;
  size_t childbp_len, childfc_len, childs_pages_len, childs_ptrs_len, childs_comp_len, childr_ptrs_len, childr_comp_len,
    childx_ptrs_len, childx_comp_len;
  size_t pioneerbp_len, pior_ptrs_len, pior_comp_len, piom_ptrs_len, piom_comp_len;

  double seconds;


  childbp = (UINT4 *) Access_mmap(&childbp_fd,&childbp_len,childbpfile,sizeof(UINT4),/*randomp*/true);
  childfc = (UINT4 *) Access_mmap(&childfc_fd,&childfc_len,childfcfile,sizeof(UINT4),/*randomp*/true);

  if (Access_file_exists_p(childs_pagesfile) == false) {
    childs_pages = (UINT4 *) NULL;
  } else {
    childs_pages = (UINT4 *) Access_allocated(&childs_pages_len,&seconds,childs_pagesfile,sizeof(UINT4));
  }
  childs_ptrs = (UINT4 *) Access_allocated(&childs_ptrs_len,&seconds,childs_ptrsfile,sizeof(UINT4));
  childs_comp = (UINT4 *) Access_mmap(&childs_comp_fd,&childs_comp_len,childs_compfile,sizeof(UINT4),
				  /*randomp*/true);

  childr_ptrs = (UINT4 *) Access_allocated(&childr_ptrs_len,&seconds,childr_ptrsfile,sizeof(UINT4));
  childr_comp = (UINT4 *) Access_mmap(&childr_comp_fd,&childr_comp_len,childr_compfile,sizeof(UINT4),
				  /*randomp*/true);

  childx_ptrs = (UINT4 *) Access_allocated(&childx_ptrs_len,&seconds,childx_ptrsfile,sizeof(UINT4));
  childx_comp = (UINT4 *) Access_mmap(&childx_comp_fd,&childx_comp_len,childx_compfile,sizeof(UINT4),
				  /*randomp*/true);

  pioneerbp = (UINT4 *) Access_mmap(&pioneerbp_fd,&pioneerbp_len,pioneerbpfile,sizeof(UINT4),
				    /*randomp*/true);
  pior_ptrs = (UINT4 *) Access_allocated(&pior_ptrs_len,&seconds,pior_ptrsfile,sizeof(UINT4));
  pior_comp = (UINT4 *) Access_mmap(&pior_comp_fd,&pior_comp_len,pior_compfile,sizeof(UINT4),
				  /*randomp*/true);

  piom_ptrs = (UINT4 *) Access_allocated(&piom_ptrs_len,&seconds,piom_ptrsfile,sizeof(UINT4));
  piom_comp = (UINT4 *) Access_mmap(&piom_comp_fd,&piom_comp_len,piom_compfile,sizeof(UINT4),
				  /*randomp*/true);

  BP_read_setup(BLOCKSIZE,SELECT_SAMPLING_INTERVAL);
  Bitpack64_read_setup();


  printf("selecti\tchildbp\tpioneerbp\tmatch\tranki\n");

  bplength = 2 * (BP_size_t) n;
  if (endl == 0) {
    startl = 0;
    endl = bplength;
  }

  for (selecti = startl; selecti < endl; selecti++) {
    if (selecti % BLOCKSIZE == 0) {
      printf(">Block %u: excess %u\n",
	     selecti/BLOCKSIZE,Bitpack64_access(selecti/BLOCKSIZE,childx_ptrs,childx_comp));
    }

    printf("%llu",selecti);
    if (get_bit(childbp,selecti) == close_paren) {
      printf("\t  )");
    } else {
      printf("\t(  ");
    }
    if (get_bit(pioneerbp,selecti) == close_paren) {
      printf("\t");
    } else {
      printf("\tP");
      /* printf("\t%u",Bitpack64_access(rank_q,piom_ptrs,piom_comp)); */
    }

    if (get_bit(childbp,selecti) == close_paren) {
      match = BP_find_openparen(selecti,childbp,childx_ptrs,childx_comp,pioneerbp,pior_ptrs,pior_comp,
				piom_ptrs,piom_comp,BLOCKSIZE);
      printf("\t%llu",match);
      if (match/BLOCKSIZE != selecti/BLOCKSIZE) {
	printf("*");
      }
      ranki = BP_rank_close(selecti,childbp,childr_ptrs,childr_comp,BLOCKSIZE);
      printf("\t  %u",ranki);

      if (get_bit(childfc,ranki) == 0) {
	printf("\tright-sib");
      } else {
	printf("\tfirst-sib");
      }

    } else {
      match = BP_find_closeparen(selecti,childbp,childx_ptrs,childx_comp,pioneerbp,pior_ptrs,pior_comp,
				piom_ptrs,piom_comp,BLOCKSIZE);
      printf("\t%llu",match);
      if (match/BLOCKSIZE != selecti/BLOCKSIZE) {
	printf("*");
      }
      printf("\t%u  ",BP_rank_open(selecti,childbp,childr_ptrs,childr_comp,BLOCKSIZE));
    }
    
    printf("\n");
  }

  return;
}

#else

void
Sarray_child_uncompress (Genome_T genomecomp, char *plcpptrsfile, char *plcpcompfile, char *childptrsfile, char *childcompfile,
			 char *sanextpfile, char *sarrayfile, UINT4 genomelength, UINT4 start, UINT4 end) {
  UINT4 n = genomelength, pos, h;
  unsigned char *gbuffer;

  UINT4 *SA, *plcpptrs, *plcpcomp, *childptrs, *childcomp, *nextp;
  int sa_fd, childptrs_fd, childcomp_fd, nextp_fd, plcpcomp_fd;
  size_t sa_len, childptrs_len, childcomp_len, nextp_len, plcpptrs_len, plcpcomp_len;
  double seconds;
  UINT4 sa_i, lcp_i, child_i, sa_nexti, lcp_nexti, child_nexti;


  gbuffer = (unsigned char *) CALLOC(n+1,sizeof(unsigned char));
  Genome_fill_buffer_simple(genomecomp,/*left*/0,/*length*/n,gbuffer);
  gbuffer[n] = 0;		       /* '\0', terminator */

  if (end == 0) {
    start = 0;
    end = n;
  }

  SA = (UINT4 *) Access_mmap(&sa_fd,&sa_len,sarrayfile,sizeof(UINT4),/*randomp*/false);
  plcpptrs = (UINT4 *) Access_allocated(&plcpptrs_len,&seconds,plcpptrsfile,sizeof(UINT4));
  plcpcomp = (UINT4 *) Access_mmap(&plcpcomp_fd,&plcpcomp_len,plcpcompfile,sizeof(UINT4),
				  /*randomp*/true);
  plcpcomp = (UINT4 *) Access_mmap(&plcpcomp_fd,&plcpcomp_len,plcpcompfile,sizeof(UINT4),
				  /*randomp*/true);
  childptrs = (UINT4 *) Access_mmap(&childptrs_fd,&childptrs_len,childptrsfile,sizeof(UINT4),/*randomp*/false);
  childcomp = (UINT4 *) Access_mmap(&childcomp_fd,&childcomp_len,childcompfile,sizeof(UINT4),/*randomp*/false);
  nextp = (UINT4 *) Access_mmap(&nextp_fd,&nextp_len,sanextpfile,sizeof(UINT4),/*randomp*/false);

  Bitpack64_read_setup();

  printf("i\tSA\tLCP\tChild\n");

  pos = start;
  sa_i = SA[pos];
  lcp_i = Bitpack64_offsetptr_only(sa_i,plcpptrs,plcpcomp) - sa_i;
  child_i = Bitpack64_access(pos,childptrs,childcomp);

  sa_nexti = SA[pos+1];
  lcp_nexti = Bitpack64_offsetptr_only(sa_nexti,plcpptrs,plcpcomp) - sa_nexti;
  child_nexti = Bitpack64_access(pos+1,childptrs,childcomp);

  if (pos == 0) {
    /* Print lcp[0] as -1 */
    printf("%u\t%u\t-1\t",pos,sa_i);
#ifdef NO_ENCODING
    printf("%u\tdown",child_i);
#else
    printf("0");
#endif
    printf("\n");
  }

  for (pos = start + 1; pos < end; pos++) {
    sa_i = sa_nexti;
    lcp_i = lcp_nexti;
    child_i = child_nexti;

    sa_nexti = SA[pos+1];
    lcp_nexti = Bitpack64_offsetptr_only(sa_nexti,plcpptrs,plcpcomp) - sa_nexti;
    child_nexti = Bitpack64_access(pos+1,childptrs,childcomp);

    printf("%u\t%u\t%u\t",pos,sa_i,lcp_i);
    if (get_bit(nextp,pos) != 0) {
#ifdef NO_ENCODING
      printf("%u\tR-next",child_i);
#else
      printf("%u\tR-next",decode_next(child_i,pos));
#endif
    } else if (lcp_i > lcp_nexti) {
#ifdef NO_ENCODING
      printf("%u\tL",child_i);
#else
      printf("%u\tL",decode_up(child_i,pos));
#endif
    } else {
#ifdef NO_ENCODING
      printf("%u\tR-down",child_i);
#else
      printf("%u\tR-down",decode_down(child_i,pos));
#endif
    }

    printf("\t%c",gbuffer[sa_i]);
    for (h = 1; (h <= lcp_i || h <= lcp_nexti) && h < 25; h++) {
      if (gbuffer[sa_i+h] == '\0') {
	printf("$");
      } else {
	printf("%c",gbuffer[sa_i+h]);
      }
    }
    if (h <= lcp_i || h <= lcp_nexti) {
      printf("...");
    }
    
    printf("\n");
  }

  sa_i = sa_nexti;
  lcp_i = lcp_nexti;
  child_i = child_nexti;
  if (pos == n) {
    printf("%u\t%u\t%u\t",pos,sa_i,lcp_i);
#ifdef NO_ENCODING
    printf("%u\tL",child_i);
#else
    printf("%u\tL",decode_up(child_i,pos));
#endif

    printf("\t%c",gbuffer[sa_i]);
    for (h = 1; h <= lcp_i && h < 25; h++) {
      if (gbuffer[sa_i+h] == '\0') {
	printf("$");
      } else {
	printf("%c",gbuffer[sa_i+h]);
      }
    }
    if (h <= lcp_i) {
      printf("...");
    }

    printf("\n");
  }


  munmap((void *) nextp,nextp_len);
  close(nextp_fd);
  munmap((void *) childcomp,childcomp_len);
  close(childcomp_fd);
  munmap((void *) childptrs,childptrs_len);
  close(childptrs_fd);
  FREE(plcpptrs);
  munmap((void *) plcpcomp,plcpcomp_len);
  close(plcpcomp_fd);
  munmap((void *) SA,sa_len);
  close(sa_fd);

  return;
}

#endif


#ifdef USE_CHILD_BP

#define clear_end(diff,enddiscard) (diff & ~(~0U << (enddiscard)))

void
Sarray_child_test (char *sarrayfile, char *plcpptrsfile, char *plcpcompfile,
		   char *childbpfile, char *childs_pagesfile, char *childs_ptrsfile, char *childs_compfile,
		   char *childr_ptrsfile, char *childr_compfile, char *childx_ptrsfile, char *childx_compfile,
		   char *pioneerbpfile, char *pior_ptrsfile, char *pior_compfile,
		   char *piom_ptrsfile, char *piom_compfile, UINT4 genomelength) {
  UINT8 i, j;
  UINT4 ranki, rankj;
  UINT4 n = genomelength;
  UINT4 *SA, *plcpptrs, *plcpcomp;
  UINT4 *childbp, *pioneerbp;
  UINT4 *childs_pages, *childs_ptrs, *childs_comp, *childr_ptrs, *childr_comp, *childx_ptrs, *childx_comp;
  UINT4 *pior_ptrs, *pior_comp, *piom_ptrs, *piom_comp;

  double seconds;
  int sa_fd, plcpcomp_fd;
  int childbp_fd, childs_comp_fd, childr_comp_fd, childx_comp_fd;
  int pioneerbp_fd, pior_comp_fd, piom_comp_fd;
  size_t sa_len, plcpptrs_len, plcpcomp_len;
  size_t childbp_len, childs_pages_len, childs_ptrs_len, childs_comp_len, childr_ptrs_len, childr_comp_len,
    childx_ptrs_len, childx_comp_len;
  size_t pioneerbp_len, pior_ptrs_len, pior_comp_len, piom_ptrs_len, piom_comp_len;


  SA = (UINT4 *) Access_mmap(&sa_fd,&sa_len,sarrayfile,sizeof(UINT4),/*randomp*/true);

  plcpptrs = (UINT4 *) Access_allocated(&plcpptrs_len,&seconds,plcpptrsfile,sizeof(UINT4));
  plcpcomp = (UINT4 *) Access_mmap(&plcpcomp_fd,&plcpcomp_len,plcpcompfile,sizeof(UINT4),
				  /*randomp*/true);
  
  childbp = (UINT4 *) Access_mmap(&childbp_fd,&childbp_len,childbpfile,sizeof(UINT4),
				  /*randomp*/true);

  if (Access_file_exists_p(childs_pagesfile) == false) {
    childs_pages = (UINT4 *) NULL;
  } else {
    childs_pages = (UINT4 *) Access_allocated(&childs_pages_len,&seconds,childs_pagesfile,sizeof(UINT4));
  }
  childs_ptrs = (UINT4 *) Access_allocated(&childs_ptrs_len,&seconds,childs_ptrsfile,sizeof(UINT4));
  childs_comp = (UINT4 *) Access_mmap(&childs_comp_fd,&childs_comp_len,childs_compfile,sizeof(UINT4),
				  /*randomp*/true);

  childr_ptrs = (UINT4 *) Access_allocated(&childr_ptrs_len,&seconds,childr_ptrsfile,sizeof(UINT4));
  childr_comp = (UINT4 *) Access_mmap(&childr_comp_fd,&childr_comp_len,childr_compfile,sizeof(UINT4),
				  /*randomp*/true);

  childx_ptrs = (UINT4 *) Access_allocated(&childx_ptrs_len,&seconds,childx_ptrsfile,sizeof(UINT4));
  childx_comp = (UINT4 *) Access_mmap(&childx_comp_fd,&childx_comp_len,childx_compfile,sizeof(UINT4),
				  /*randomp*/true);

  pioneerbp = (UINT4 *) Access_mmap(&pioneerbp_fd,&pioneerbp_len,pioneerbpfile,sizeof(UINT4),
				    /*randomp*/true);
  pior_ptrs = (UINT4 *) Access_allocated(&pior_ptrs_len,&seconds,pior_ptrsfile,sizeof(UINT4));
  pior_comp = (UINT4 *) Access_mmap(&pior_comp_fd,&pior_comp_len,pior_compfile,sizeof(UINT4),
				  /*randomp*/true);

  piom_ptrs = (UINT4 *) Access_allocated(&piom_ptrs_len,&seconds,piom_ptrsfile,sizeof(UINT4));
  piom_comp = (UINT4 *) Access_mmap(&piom_comp_fd,&piom_comp_len,piom_compfile,sizeof(UINT4),
				  /*randomp*/true);

  BP_read_setup(BLOCKSIZE,SELECT_SAMPLING_INTERVAL);
  Bitpack64_read_setup();

  for (ranki = 1; ranki <= n; ranki++) {
    j = BP_select(ranki,childbp,childs_pages,childs_ptrs,childs_comp,
		  BLOCKSIZE/2,BLOCKSIZE,SELECT_SAMPLING_INTERVAL);
    printf("Select of %u is %llu\n",ranki,j);
    rankj = BP_rank_open(j,childbp,childr_ptrs,childr_comp,BLOCKSIZE);
    printf("Rank of %llu is %u\n",j,rankj);
    if (rankj != ranki) {
      abort();
    }
  }
  exit(0);

  for (i = 0; i < 2*n; i++) {
    if ((childbp[i/BLOCKSIZE] & (1 << (i % BLOCKSIZE))) != 0) {
      j = BP_find_closeparen(i,childbp,childx_ptrs,childx_comp,pioneerbp,pior_ptrs,pior_comp,
			     piom_ptrs,piom_comp,BLOCKSIZE);
      printf("Close paren of %d is %u\n",i,j);
      j = BP_find_openparen(j,childbp,childx_ptrs,childx_comp,pioneerbp,pior_ptrs,pior_comp,
			    piom_ptrs,piom_comp,BLOCKSIZE);
      if (j != i) {
	printf("But open paren of j yields %u\n",j);
	abort();
      }
    } else {
      j = BP_find_openparen(i,childbp,childx_ptrs,childx_comp,pioneerbp,pior_ptrs,pior_comp,
			    piom_ptrs,piom_comp,BLOCKSIZE);
      printf("Open paren of %d is %u\n",i,j);
      j = BP_find_closeparen(j,childbp,childx_ptrs,childx_comp,pioneerbp,pior_ptrs,pior_comp,
			     piom_ptrs,piom_comp,BLOCKSIZE);
      if (j != i) {
	printf("But close paren of j yields %u\n",j);
	abort();
      }
    }
  }

  return;
}

#endif


