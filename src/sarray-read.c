static char rcsid[] = "$Id: sarray-read.c 110160 2013-10-05 01:45:16Z twu $";
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
#include "stage3hr.h"
#include "bitpack64-access.h"


#ifdef HAVE_SSE2
#include <emmintrin.h>
#endif


/* A value of 10000 misses various splices, although they are caught by GSNAP algorithm */
#define EXCESS_SARRAY_HITS 100000
#define GUESS_ALLOCATION 10


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



#define T Sarray_T
struct T {
  Univcoord_T *array;
#ifdef USE_LCP
  Univcoord_T *lcp;
#endif
  UINT4 *lcpptrs;
  UINT4 *lcpcomp;
  Univcoord_T n_plus_one;

  Sarrayptr_T *saindex;
  int indexsize;

  Access_T access;
  int array_fd;
  size_t array_len;
#ifdef USE_LCP
  int lcp_fd;
  size_t lcp_len;
#endif
  int lcpptrs_fd;
  size_t lcpptrs_len;
  int lcpcomp_fd;
  size_t lcpcomp_len;

  int saindex_fd;
  size_t saindex_len;
};


/* For benchmarking */
UINT4 *
Sarray_lcpptrs (Sarray_T this) {
  return this->lcpptrs;
}

/* For benchmarking */
UINT4 *
Sarray_lcpcomp (Sarray_T this) {
  return this->lcpcomp;
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


void
Sarray_setup (T sarray_in, Genome_T genome_in, Univ_IIT_T chromosome_iit_in, int circular_typeint_in,
	      Chrpos_T shortsplicedist_in, int splicing_penalty_in,
	      int max_deletionlength, int max_end_deletions_in,
	      int max_middle_insertions, int max_end_insertions,
	      Univcoord_T *splicesites_in, Splicetype_T *splicetypes_in,
	      Chrpos_T *splicedists_in, int nsplicesites_in) {
#ifdef USE_LCP
  Univcoord_T position;
#endif

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

  Bitpack64_access_setup(sarray->lcpptrs,sarray->lcpcomp);
#ifdef USE_LCP
  fprintf(stderr,"Checking validity of compressed LCP file\n");
  for (position = 0; position < sarray->n_plus_one; position++) {
    if (sarray->lcp[position] != Bitpack64_access(position)) {
      abort();
    }
  }
#endif

  return;
}


T
Sarray_new (char *directory, char *fileroot, char *snps_root, Access_mode_T access) {
  T new = (T) MALLOC(sizeof(*new));
  char *filename;
  char *comma1, *comma2;
  double seconds1, seconds2, seconds;
  int npages;

  if (snps_root != NULL) {
    /* Always allocate saindex and lcpptrs */
    filename = (char *) CALLOC(strlen(directory)+strlen("/")+strlen(fileroot)+strlen(".saindex.")+
			       strlen(snps_root)+1,sizeof(char));
    sprintf(filename,"%s/%s.saindex.%s",directory,fileroot,snps_root);
    if (Access_file_exists_p(filename) == false) {
      fprintf(stderr,"Suffix array index file %s does not exist\n",filename);
      FREE(filename);
      FREE(new);
      return (T) NULL;
    } else {
      new->saindex = (unsigned int *) Access_allocated(&new->saindex_len,&seconds,filename,sizeof(unsigned int));
      new->indexsize = 12;
      FREE(filename);
    }

    filename = (char *) CALLOC(strlen(directory)+strlen("/")+strlen(fileroot)+strlen(".salcpptrs.")+
			       strlen(snps_root)+1,sizeof(char));
    sprintf(filename,"%s/%s.salcpptrs.%s",directory,fileroot,snps_root);
    if (Access_file_exists_p(filename) == false) {
      fprintf(stderr,"Suffix array lcp file %s does not exist\n",filename);
      FREE(filename);
      FREE(new->saindex);
      FREE(new);
      return (T) NULL;
    } else {
      new->lcpptrs = (UINT4 *) Access_allocated(&new->lcpptrs_len,&seconds2,filename,sizeof(UINT4));
      FREE(filename);
    }


    if (access == USE_ALLOCATE) {
      fprintf(stderr,"Allocating memory for suffix array...");

      filename = (char *) CALLOC(strlen(directory)+strlen("/")+strlen(fileroot)+strlen(".sarray.")+
				 strlen(snps_root)+1,sizeof(char));
      sprintf(filename,"%s/%s.sarray.%s",directory,fileroot,snps_root);
      if (Access_file_exists_p(filename) == false) {
	fprintf(stderr,"Suffix array file %s does not exist\n",filename);
	FREE(filename);
	FREE(new->saindex);
	FREE(new->lcpptrs);
	FREE(new);
	return (T) NULL;
      } else {
	new->array = (Univcoord_T *) Access_allocated(&new->array_len,&seconds1,filename,sizeof(Univcoord_T));
	FREE(filename);
      }

#ifdef USE_LCP
      filename = (char *) CALLOC(strlen(directory)+strlen("/")+strlen(fileroot)+strlen(".salcp.")+
				 strlen(snps_root)+1,sizeof(char));
      sprintf(filename,"%s/%s.salcp.%s",directory,fileroot,snps_root);
      if (Access_file_exists_p(filename) == false) {
	fprintf(stderr,"Suffix array lcp file %s does not exist\n",filename);
	FREE(filename);
	FREE(new->saindex);
	FREE(new->lcpptrs);
	FREE(new->array);
	FREE(new);
	return (T) NULL;
      } else {
	new->lcp = (Univcoord_T *) Access_allocated(&new->lcp_len,&seconds2,filename,sizeof(Univcoord_T));
	FREE(filename);
      }
#endif

      filename = (char *) CALLOC(strlen(directory)+strlen("/")+strlen(fileroot)+strlen(".salcpcomp.")+
				 strlen(snps_root)+1,sizeof(char));
      sprintf(filename,"%s/%s.salcpcomp.%s",directory,fileroot,snps_root);
      if (Access_file_exists_p(filename) == false) {
	fprintf(stderr,"Suffix array lcp file %s does not exist\n",filename);
	FREE(filename);
	FREE(new->saindex);
	FREE(new->lcpptrs);
	FREE(new->array);
	FREE(new);
	return (T) NULL;
      } else {
	new->lcpcomp = (Univcoord_T *) Access_allocated(&new->lcpcomp_len,&seconds2,filename,sizeof(Univcoord_T));
	FREE(filename);
      }

      if (new->array == NULL || new->lcpcomp == NULL) {
	fprintf(stderr,"insufficient memory (need to use memory mapping)\n");
	exit(9);
      } else {
	comma1 = Genomicpos_commafmt(new->array_len);
	comma2 = Genomicpos_commafmt(new->lcpcomp_len);
	fprintf(stderr,"done (%s + %s bytes, %.2f sec)\n",comma1,comma2,seconds1+seconds2);
	FREE(comma2);
	FREE(comma1);
	new->access = ALLOCATED;
      }

    } else if (access == USE_MMAP_PRELOAD) {
      fprintf(stderr,"Memory mapping suffix array...");

      filename = (char *) CALLOC(strlen(directory)+strlen("/")+strlen(fileroot)+strlen(".sarray.")+
				 strlen(snps_root)+1,sizeof(char));
      sprintf(filename,"%s/%s.sarray.%s",directory,fileroot,snps_root);
      if (Access_file_exists_p(filename) == false) {
	fprintf(stderr,"Suffix array file %s does not exist\n",filename);
	FREE(filename);
	FREE(new->saindex);
	FREE(new->lcpptrs);
	FREE(new);
	return (T) NULL;
      } else {
	new->array = (Univcoord_T *) Access_mmap_and_preload(&new->array_fd,&new->array_len,&npages,&seconds1,
							     filename,sizeof(Univcoord_T));
	FREE(filename);
      }

#ifdef USE_LCP
      filename = (char *) CALLOC(strlen(directory)+strlen("/")+strlen(fileroot)+strlen(".salcp.")+
				 strlen(snps_root)+1,sizeof(char));
      sprintf(filename,"%s/%s.salcp.%s",directory,fileroot,snps_root);
      if (Access_file_exists_p(filename) == false) {
	fprintf(stderr,"Suffix array lcp file %s does not exist\n",filename);
	FREE(filename);
	FREE(new->saindex);
	FREE(new->lcpptrs);
	munmap((void *) new->array,new->array_len);
	FREE(new);
	return (T) NULL;
      } else {
	new->lcp = (Univcoord_T *) Access_mmap_and_preload(&new->lcp_fd,&new->lcp_len,&npages,&seconds,
							   filename,sizeof(Univcoord_T));
	FREE(filename);
      }
#endif

      filename = (char *) CALLOC(strlen(directory)+strlen("/")+strlen(fileroot)+strlen(".salcpcomp.")+
				 strlen(snps_root)+1,sizeof(char));
      sprintf(filename,"%s/%s.salcpcomp.%s",directory,fileroot,snps_root);
      if (Access_file_exists_p(filename) == false) {
	fprintf(stderr,"Suffix array lcpcomp file %s does not exist\n",filename);
	FREE(filename);
	FREE(new->saindex);
	FREE(new->lcpptrs);
	munmap((void *) new->array,new->array_len);
	FREE(new);
	return (T) NULL;
      } else {
	new->lcpcomp = (Univcoord_T *) Access_mmap_and_preload(&new->lcpcomp_fd,&new->lcpcomp_len,&npages,&seconds2,
							       filename,sizeof(Univcoord_T));
	FREE(filename);
      }

      if (new->array == NULL || new->lcpcomp == NULL) {
	fprintf(stderr,"insufficient memory (need more virtual memory or run without suffix array)\n");
	exit(9);
      } else {
	comma1 = Genomicpos_commafmt(new->array_len);
	comma2 = Genomicpos_commafmt(new->lcpcomp_len);
	fprintf(stderr,"done (%s + %s bytes, %.2f sec)\n",comma1,comma2,seconds1+seconds2);
	FREE(comma2);
	FREE(comma1);
	new->access = MMAPPED;
      }

    } else if (access == USE_MMAP_ONLY) {
      fprintf(stderr,"Memory mapping suffix array...");

      filename = (char *) CALLOC(strlen(directory)+strlen("/")+strlen(fileroot)+strlen(".sarray.")+
				 strlen(snps_root)+1,sizeof(char));
      sprintf(filename,"%s/%s.sarray.%s",directory,fileroot,snps_root);
      if (Access_file_exists_p(filename) == false) {
	fprintf(stderr,"Suffix array file %s does not exist\n",filename);
	FREE(filename);
	FREE(new->saindex);
	FREE(new->lcpptrs);
	FREE(new);
	return (T) NULL;
      } else {
	new->array = (Univcoord_T *) Access_mmap(&new->array_fd,&new->array_len,filename,sizeof(Univcoord_T),
						 /*randomp*/true);
	FREE(filename);
      }

#ifdef USE_LCP
      filename = (char *) CALLOC(strlen(directory)+strlen("/")+strlen(fileroot)+strlen(".salcp.")+
				 strlen(snps_root)+1,sizeof(char));
      sprintf(filename,"%s/%s.salcp.%s",directory,fileroot,snps_root);
      if (Access_file_exists_p(filename) == false) {
	fprintf(stderr,"Suffix array lcp file %s does not exist\n",filename);
	FREE(filename);
	FREE(new->saindex);
	FREE(new->lcpptrs);
	munmap((void *) new->array,new->array_len);
	FREE(new);
	return (T) NULL;
      } else {
	new->lcp = (Univcoord_T *) Access_mmap(&new->lcp_fd,&new->lcp_len,filename,sizeof(Univcoord_T),
					       /*randomp*/true);
	FREE(filename);
      }
#endif

      filename = (char *) CALLOC(strlen(directory)+strlen("/")+strlen(fileroot)+strlen(".salcpcomp.")+
				 strlen(snps_root)+1,sizeof(char));
      sprintf(filename,"%s/%s.salcpcomp.%s",directory,fileroot,snps_root);
      if (Access_file_exists_p(filename) == false) {
	fprintf(stderr,"Suffix array file %s does not exist\n",filename);
	FREE(filename);
	FREE(new->saindex);
	FREE(new->lcpptrs);
	FREE(new);
	return (T) NULL;
      } else {
	new->lcpcomp = (Univcoord_T *) Access_mmap(&new->lcpcomp_fd,&new->lcpcomp_len,filename,sizeof(Univcoord_T),
						 /*randomp*/true);
	FREE(filename);
      }

      if (new->array == NULL || new->lcpcomp == NULL) {
	fprintf(stderr,"insufficient memory (need more virtual memory or run without suffix array)\n");
	exit(9);
      } else {
	comma1 = Genomicpos_commafmt(new->array_len);
	comma2 = Genomicpos_commafmt(new->lcpcomp_len);
	fprintf(stderr,"done (%s + %s bytes)\n",comma1,comma2);
	FREE(comma2);
	FREE(comma1);
	new->access = MMAPPED;
      }
    }

  } else {
    /* Always allocate saindex and salcpptrs */
    filename = (char *) CALLOC(strlen(directory)+strlen("/")+strlen(fileroot)+strlen(".saindex")+1,sizeof(char));
    sprintf(filename,"%s/%s.saindex",directory,fileroot);
    if (Access_file_exists_p(filename) == false) {
      fprintf(stderr,"Suffix array index file %s does not exist\n",filename);
      FREE(filename);
      FREE(new);
      return (T) NULL;
    } else {
      new->saindex = (Sarrayptr_T *) Access_allocated(&new->saindex_len,&seconds,filename,sizeof(Sarrayptr_T));
      new->indexsize = 12;
      FREE(filename);
    }

    filename = (char *) CALLOC(strlen(directory)+strlen("/")+strlen(fileroot)+strlen(".salcpptrs")+1,sizeof(char));
    sprintf(filename,"%s/%s.salcpptrs",directory,fileroot);
    if (Access_file_exists_p(filename) == false) {
      fprintf(stderr,"Suffix array lcp file %s does not exist\n",filename);
      FREE(filename);
      FREE(new->saindex);
      FREE(new);
      return (T) NULL;
    } else {
      new->lcpptrs = (UINT4 *) Access_allocated(&new->lcpptrs_len,&seconds,filename,sizeof(UINT4));
      FREE(filename);
    }


    if (access == USE_ALLOCATE) {
      fprintf(stderr,"Allocating memory for suffix array...");

      filename = (char *) CALLOC(strlen(directory)+strlen("/")+strlen(fileroot)+strlen(".sarray")+1,sizeof(char));
      sprintf(filename,"%s/%s.sarray",directory,fileroot);
      if (Access_file_exists_p(filename) == false) {
	fprintf(stderr,"Suffix array file %s does not exist\n",filename);
	FREE(filename);
	FREE(new->saindex);
	FREE(new->lcpptrs);
	FREE(new);
	return (T) NULL;
      } else {
	new->array = (Univcoord_T *) Access_allocated(&new->array_len,&seconds1,filename,sizeof(Univcoord_T));
	FREE(filename);
      }

#ifdef USE_LCP
      filename = (char *) CALLOC(strlen(directory)+strlen("/")+strlen(fileroot)+strlen(".salcp")+1,sizeof(char));
      sprintf(filename,"%s/%s.salcp",directory,fileroot);
      if (Access_file_exists_p(filename) == false) {
	fprintf(stderr,"Suffix array lcp file %s does not exist\n",filename);
	FREE(filename);
	FREE(new->saindex);
	FREE(new->lcpptrs);
	FREE(new->array);
	FREE(new);
	return (T) NULL;
      } else {
	new->lcp = (Univcoord_T *) Access_allocated(&new->lcp_len,&seconds2,filename,sizeof(Univcoord_T));
	FREE(filename);
      }
#endif

      filename = (char *) CALLOC(strlen(directory)+strlen("/")+strlen(fileroot)+strlen(".salcpcomp")+1,sizeof(char));
      sprintf(filename,"%s/%s.salcpcomp",directory,fileroot);
      if (Access_file_exists_p(filename) == false) {
	fprintf(stderr,"Suffix array lcp file %s does not exist\n",filename);
	FREE(filename);
	FREE(new->saindex);
	FREE(new->lcpptrs);
	FREE(new->array);
	FREE(new);
	return (T) NULL;
      } else {
	new->lcpcomp = (Univcoord_T *) Access_allocated(&new->lcpcomp_len,&seconds2,filename,sizeof(Univcoord_T));
	FREE(filename);
      }

      if (new->array == NULL || new->lcpcomp == NULL) {
	fprintf(stderr,"insufficient memory (need to use memory mapping)\n");
	exit(9);
      } else {
	comma1 = Genomicpos_commafmt(new->array_len);
	comma2 = Genomicpos_commafmt(new->lcpcomp_len);
	fprintf(stderr,"done (%s + %s bytes, %.2f sec)\n",comma1,comma2,seconds1+seconds2);
	FREE(comma2);
	FREE(comma1);
	new->access = ALLOCATED;
      }

    } else if (access == USE_MMAP_PRELOAD) {
      fprintf(stderr,"Memory mapping suffix array...");

      filename = (char *) CALLOC(strlen(directory)+strlen("/")+strlen(fileroot)+strlen(".sarray")+1,sizeof(char));
      sprintf(filename,"%s/%s.sarray",directory,fileroot);
      if (Access_file_exists_p(filename) == false) {
	fprintf(stderr,"Suffix array file %s does not exist\n",filename);
	FREE(filename);
	FREE(new->saindex);
	FREE(new->lcpptrs);
	FREE(new);
	return (T) NULL;
      } else {
	new->array = (Univcoord_T *) Access_mmap_and_preload(&new->array_fd,&new->array_len,&npages,&seconds,
							     filename,sizeof(Univcoord_T));
	FREE(filename);
      }

#ifdef USE_LCP
      filename = (char *) CALLOC(strlen(directory)+strlen("/")+strlen(fileroot)+strlen(".salcp")+1,sizeof(char));
      sprintf(filename,"%s/%s.salcp",directory,fileroot);
      if (Access_file_exists_p(filename) == false) {
	fprintf(stderr,"Suffix array lcp file %s does not exist\n",filename);
	FREE(filename);
	FREE(new->saindex);
	FREE(new->lcpptrs);
	munmap((void *) new->array,new->array_len);
	FREE(new);
	return (T) NULL;
      } else {
	new->lcp = (Univcoord_T *) Access_mmap_and_preload(&new->lcp_fd,&new->lcp_len,&npages,&seconds,
							   filename,sizeof(Univcoord_T));
	FREE(filename);
      }
#endif

      filename = (char *) CALLOC(strlen(directory)+strlen("/")+strlen(fileroot)+strlen(".salcpcomp")+1,sizeof(char));
      sprintf(filename,"%s/%s.salcpcomp",directory,fileroot);
      if (Access_file_exists_p(filename) == false) {
	fprintf(stderr,"Suffix array lcpcomp file %s does not exist\n",filename);
	FREE(filename);
	FREE(new->saindex);
	FREE(new->lcpptrs);
	munmap((void *) new->array,new->array_len);
	FREE(new);
	return (T) NULL;
      } else {
	new->lcpcomp = (UINT4 *) Access_mmap_and_preload(&new->lcpcomp_fd,&new->lcpcomp_len,&npages,&seconds,
							 filename,sizeof(UINT4));
	FREE(filename);
      }

      if (new->array == NULL || new->lcpcomp == NULL) {
	fprintf(stderr,"insufficient memory (need more virtual memory or run without suffix array)\n");
	exit(9);
      } else {
	comma1 = Genomicpos_commafmt(new->array_len);
	comma2 = Genomicpos_commafmt(new->lcpcomp_len);
	fprintf(stderr,"done (%s + %s bytes)\n",comma1,comma2);
	FREE(comma2);
	FREE(comma1);
	new->access = MMAPPED;
      }

    } else if (access == USE_MMAP_ONLY) {
      fprintf(stderr,"Memory mapping suffix array...");

      filename = (char *) CALLOC(strlen(directory)+strlen("/")+strlen(fileroot)+strlen(".sarray")+1,sizeof(char));
      sprintf(filename,"%s/%s.sarray",directory,fileroot);
      if (Access_file_exists_p(filename) == false) {
	fprintf(stderr,"Suffix array file %s does not exist\n",filename);
	FREE(filename);
	FREE(new->saindex);
	FREE(new->lcpptrs);
	FREE(new);
	return (T) NULL;
      } else {
	new->array = (Univcoord_T *) Access_mmap(&new->array_fd,&new->array_len,
						 filename,sizeof(Univcoord_T),/*randomp*/true);
	FREE(filename);
      }

#ifdef USE_LCP
      filename = (char *) CALLOC(strlen(directory)+strlen("/")+strlen(fileroot)+strlen(".salcp")+1,sizeof(char));
      sprintf(filename,"%s/%s.salcp",directory,fileroot);
      if (Access_file_exists_p(filename) == false) {
	fprintf(stderr,"Suffix array lcp file %s does not exist\n",filename);
	FREE(filename);
	FREE(new->saindex);
	FREE(new->lcpptrs);
	munmap((void *) new->array,new->array_len);
	FREE(new);
	return (T) NULL;
      } else {
	new->lcp = (Univcoord_T *) Access_mmap(&new->lcp_fd,&new->lcp_len,
					       filename,sizeof(Univcoord_T),/*randomp*/true);
	FREE(filename);
      }
#endif

      filename = (char *) CALLOC(strlen(directory)+strlen("/")+strlen(fileroot)+strlen(".salcpcomp")+1,sizeof(char));
      sprintf(filename,"%s/%s.salcpcomp",directory,fileroot);
      if (Access_file_exists_p(filename) == false) {
	fprintf(stderr,"Suffix array lcp file %s does not exist\n",filename);
	FREE(filename);
	FREE(new->saindex);
	FREE(new->lcpptrs);
	munmap((void *) new->array,new->array_len);
	FREE(new);
	return (T) NULL;
      } else {
	new->lcpcomp = (Univcoord_T *) Access_mmap(&new->lcpcomp_fd,&new->lcpcomp_len,
						   filename,sizeof(Univcoord_T),/*randomp*/true);
	FREE(filename);
      }

      if (new->array == NULL || new->lcpcomp == NULL) {
	fprintf(stderr,"insufficient memory (need more virtual memory or run without suffix array)\n");
	exit(9);
      } else {
	comma1 = Genomicpos_commafmt(new->array_len);
	comma2 = Genomicpos_commafmt(new->lcpcomp_len);
	fprintf(stderr,"done (%s + %s bytes)\n",comma1,comma2);
	FREE(comma2);
	FREE(comma1);
	new->access = MMAPPED;
      }
    }
  }

  /* Should be genomiclength + 1*/
  new->n_plus_one = new->array_len/sizeof(Univcoord_T);

  return new;
}


void
Sarray_free (T *old) {
  if (*old) {
    FREE((*old)->saindex);
    FREE((*old)->lcpptrs);
    if ((*old)->access == ALLOCATED) {
      FREE((*old)->array);
      FREE((*old)->lcpcomp);
#ifdef USE_LCP
      FREE((*old)->lcp);
#endif

#ifdef HAVE_MMAP
    } else if ((*old)->access == MMAPPED) {
#ifdef USE_LCP
      munmap((void *) (*old)->lcp,(*old)->lcp_len);
      close((*old)->lcp_fd);
#endif
      munmap((void *) (*old)->lcpcomp,(*old)->lcpcomp_len);
      munmap((void *) (*old)->array,(*old)->array_len);
      close((*old)->lcpcomp_fd);
      close((*old)->array_fd);
#endif
    }

    FREE(*old);
  }

  return;
}



static Sarrayptr_T
sarray_search_init (char *query, int querylength, int queryoffset, Compress_T query_compress, bool plusp,
		    Sarrayptr_T low, Sarrayptr_T high, Univcoord_T nmatches_low, Univcoord_T nmatches_high) {
  Sarrayptr_T mid;
  Univcoord_T pos;
  Univcoord_T nmatches_mid, fasti;
  char c;
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
    if ((c = Genome_get_char(genome,pos)) == 'N') {
      c = 'X';
    }

    if (fasti == (Univcoord_T) querylength || c > query[fasti]) {
      high = mid;
      /* nmatches_high = (sarray->lcp[mid] < nmatches_mid) ? sarray->lcp[mid] : nmatches_mid; */
      lcp_mid = Bitpack64_access(mid);
#ifdef USE_LCP
      if (lcp_mid != sarray->lcp[mid]) {
	fprintf(stderr,"LCP compression error at %u\n",mid);
      }
#endif
      nmatches_high = (lcp_mid < nmatches_mid) ? lcp_mid : nmatches_mid;
    } else {
      low = mid;
      /* nmatches_low = (sarray->lcp[low] < nmatches_mid) ? sarray->lcp[low] : nmatches_mid; */
      lcp_low = Bitpack64_access(low);
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


static Sarrayptr_T
sarray_search_final (char *query, int querylength, int queryoffset, Compress_T query_compress, bool plusp,
		     Sarrayptr_T low, Sarrayptr_T high, Univcoord_T nmatches_low, Univcoord_T nmatches_high) {
  Sarrayptr_T mid;
  Univcoord_T pos;
  Univcoord_T nmatches_mid, fasti;
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
    if ((c = Genome_get_char(genome,pos)) == 'N') {
      c = 'X';
    }

    if (fasti == (Univcoord_T) querylength || c < query[fasti]) {
      low = mid;
      /* nmatches_low = (sarray->lcp[low] < nmatches_mid) ? sarray->lcp[low] : nmatches_mid; */
      lcp_low = Bitpack64_access(low);
#ifdef USE_LCP
      if (lcp_low != sarray->lcp[low]) {
	fprintf(stderr,"LCP compression error at %u\n",mid);
      }
#endif
      nmatches_low = (lcp_low < nmatches_mid) ? lcp_low : nmatches_mid;
    } else {
      high = mid;
      /* nmatches_high = (sarray->lcp[mid] < nmatches_mid) ? sarray->lcp[mid] : nmatches_mid; */
      lcp_mid = Bitpack64_access(mid);
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


static void
sarray_search (Sarrayptr_T *initptr, Sarrayptr_T *finalptr, bool *successp,
	       int *nmatches, char *query, int querylength,
	       int queryoffset, Compress_T query_compress, bool plusp) {
  Sarrayptr_T low, high, mid;
  Univcoord_T pos;
  Univcoord_T nmatches_low, nmatches_high, nmatches_mid, fasti;
  UINT4 lcp_low, lcp_mid;

  Sarrayptr_T prevlow, prevhigh;
  Univcoord_T nmatches_prevlow, nmatches_prevhigh, nmatches_best = 0;

  int effective_querylength;	/* length to first N */
  Storedoligomer_T oligo;
  char c;

#ifdef DEBUG
  int i = 0;
  int recount;
  char Buffer[1000];
  Univcoord_T hit;
  bool failp;
#elif defined(DEBUG1)
  char Buffer[1000];
#endif


  *successp = false;
  effective_querylength = nt_querylength(query,querylength);

  if (effective_querylength == 0) {
    *initptr = *finalptr = 0;
    *nmatches = 0;
    return;

  } else if (effective_querylength < sarray->indexsize) {
    low = prevlow = 0;
    high = prevhigh = sarray->n_plus_one;
    
  } else {
    oligo = nt_oligo(query,sarray->indexsize);
    if (sarray->saindex[oligo] == -1U) {
      low = prevlow = 0;
    } else {
      low = prevlow = sarray->saindex[oligo] - 1;
    }

    if (sarray->saindex[oligo+1] == -1U) {
      high = prevhigh = sarray->n_plus_one;
    } else {
      high = prevhigh = sarray->saindex[oligo+1];
    }
  }

  debug1(printf("sarray_search on %s, querylength %d, with low %u, high %u\n",
		query,querylength,low,high));

  nmatches_low = nmatches_high = 0;
  while (low + 1 < high && *successp == false) {
    /* Compute mid for unsigned ints */
    mid = low/2 + high/2;
    if (low % 2 == 1 && high % 2 == 1) {
      mid += 1;
    }
    debug1(printf("low %u, high %u => mid %u\n",low,high,mid));
    nmatches_mid = (nmatches_low < nmatches_high) ? nmatches_low : nmatches_high;

    fasti = nmatches_mid +
      (Univcoord_T) Genome_consecutive_matches_rightward(query_compress,/*left*/sarray->array[mid]-queryoffset,
							 /*pos5*/queryoffset+nmatches_mid,/*pos3*/queryoffset+querylength,
							 plusp,/*genestrand*/0);
    debug1(Genome_fill_buffer_simple(genome,sarray->array[mid],querylength,Buffer));
    debug1(printf("fasti at %u is %d: %s\n",sarray->array[mid],fasti,Buffer));

    pos = sarray->array[mid] + fasti;
    if ((c = Genome_get_char(genome,pos)) == 'N') {
      c = 'X';
    }

    if (fasti > nmatches_best) {
      debug1(printf("fasti %d > nmatches_best %d.  Saving prevlow %u and prevhigh %u.\n",
		   fasti,nmatches_best,low,high));
      prevlow = low;
      prevhigh = high;
      nmatches_prevlow = nmatches_low;
      nmatches_prevhigh = nmatches_high;
      nmatches_best = fasti;
    }

    if (fasti == (Univcoord_T) querylength) {
      *successp = true;

    } else if (c < query[fasti]) {
      low = mid;
      /* nmatches_low = (sarray->lcp[low] < nmatches_mid) ? sarray->lcp[low] : nmatches_mid; */
      lcp_low = Bitpack64_access(low);
#ifdef USE_LCP
      if (lcp_low != sarray->lcp[low]) {
	fprintf(stderr,"LCP compression error at %u\n",mid);
      }
#endif
      nmatches_low = (lcp_low < nmatches_mid) ? lcp_low : nmatches_mid;
      debug1(printf("genome %c < query (%c) => low gets %u @ %u\n",c,query[fasti],low,sarray->array[low]));

    } else if (c > query[fasti]) {
      high = mid;
      /* nmatches_high = (sarray->lcp[mid] < nmatches_mid) ? sarray->lcp[mid] : nmatches_mid; */
      lcp_mid = Bitpack64_access(mid);
#ifdef USE_LCP
      if (lcp_mid != sarray->lcp[mid]) {
	fprintf(stderr,"LCP compression error at %u\n",mid);
      }
#endif
      nmatches_high = (lcp_mid < nmatches_mid) ? lcp_mid : nmatches_mid;
      debug1(printf("genome %c > query (%c) => high gets %u @ %u\n",c,query[fasti],high,sarray->array[high]));

    } else {
      debug1(printf("genome %c == query (%c) => should not happen after Genome_consecutive_matches\n",
		   c,query[fasti]));
      abort();
    }

    debug1(printf("sarray_search with low %u @ %u, high %u @ %u\n",low,sarray->array[low],high,sarray->array[high]));
  }
  debug1(printf("\n"));

  if ((*nmatches = (int) nmatches_best) == 0) {
    debug(printf("Got no matches at all\n"));
    return;
  } else if (*successp == false) {
    /* Search only on part of string that does match.  Back up to prevlow and prevhigh. */
    debug(printf("%s fail at %d: calling init/final on prevlow %u, prevhigh %u\n",
		 plusp ? "plus" : "minus",queryoffset,prevlow,prevhigh));

    *initptr = sarray_search_init(query,/*querylength*/*nmatches,queryoffset,query_compress,plusp,
				  prevlow,prevhigh,nmatches_prevlow,nmatches_prevhigh);
    *finalptr = sarray_search_final(query,/*querylength*/*nmatches,queryoffset,query_compress,plusp,
				    prevlow,prevhigh,nmatches_prevlow,nmatches_prevhigh);
    debug(printf("%s fail at %d: got %d hits with %d matches:\n",
		 plusp ? "plus" : "minus",queryoffset,(*finalptr - *initptr + 1),*nmatches));
  } else {
    debug(printf("%s success at %d: calling init/final on low %u, high %u\n",
		 plusp ? "plus" : "minus",queryoffset,low,high));

    *initptr = sarray_search_init(query,querylength,queryoffset,query_compress,plusp,
				  low,mid,nmatches_low,nmatches_mid);
    *finalptr = sarray_search_final(query,querylength,queryoffset,query_compress,plusp,
				    mid,high,nmatches_mid,nmatches_high);
    debug(printf("%s success at %d: got %d hits with %d matches:\n",
		 plusp ? "plus" : "minus",queryoffset,(*finalptr - *initptr + 1),*nmatches));
  }

  if ((int) (*finalptr - *initptr + 1) < 0) {
    abort();
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
  for (i = 0; i < (int) (*finalptr - *initptr + 1) && i < 100; i++) {
    hit = sarray->array[(*initptr)+i];
    recount = Genome_consecutive_matches_rightward(query_compress,/*left*/hit-queryoffset,
						   /*pos5*/queryoffset,/*pos3*/queryoffset+querylength,
						   plusp,/*genestrand*/0);
    printf("%d\t%u\t%u\t",recount,(*initptr)+i,hit /*+ 1U*/);
    Genome_fill_buffer_simple(genome,sarray->array[(*initptr)+i],recount+1,Buffer);
    
    printf("%s\n",Buffer);
    if (recount != *nmatches) {
      printf("querylength is %d\n",querylength);
      printf("false positive: recount %d at %u does not equal expected nmatches %d\n",
	     recount,sarray->array[(*initptr)],*nmatches);
      failp = true;
    }
    /* hits[i] = sarray->array[(*initptr)++]; */
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
  int i, j;
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
      n_prealign = (16 - (pointer & 0xF))/4;
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
	n_prealign = (pointer & 0xF)/4;
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
    middlei = (lowi+highi)/2;
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
		  List_T set, Sarray_T sarray, Compress_T query_compress,
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
		 List_T set, Sarray_T sarray, char *queryptr, Compress_T query_compress,
		 bool plusp, int genestrand, int best_querystart, int best_queryend) {
  Elt_T elt;
  int nmatches;
  Sarrayptr_T initptr, finalptr;
  bool successp;
  int queryend, querypos;
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
  int nmismatches, nmismatches1, nmismatches2, nindels;
  int nsame, ndiff;
  int querystart_diff, queryend_diff, query_indel_pos, indel_pos;

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
	debug7(printf("B deletion of %d bp relative to max_deletionlen %d...",nindels,max_deletionlen));
	if ((indel_pos < 17 || querylength - indel_pos < 17) && nindels > max_end_deletions) {
	  debug7(printf("too long for end deletion"));
	} else {
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
	}
	debug7(printf("\n"));
      
      } else if (left2 < left1) {
	nindels = left1 - left2;
	if (nindels >= indel_pos || indel_pos + nindels >= querylength) {
	  debug7(printf("X insertion of %d bp too long\n",nindels));
	} else {
	  debug7(printf("C insertion of %d bp...",nindels));
      
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
	debug7(printf("B deletion of %d bp relative to max_deletionlen %d...",nindels,max_deletionlen));
	if ((indel_pos < 17 || querylength - indel_pos < 17) && nindels > max_end_deletions) {
	  debug7(printf("too long for end deletion"));
	} else {
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
	}
	debug7(printf("\n"));
      
      } else if (left2 < left1) {
	nindels = left1 - left2;
	if (nindels >= indel_pos || indel_pos + nindels >= querylength) {
	  debug7(printf("X insertion of %d bp too long\n",nindels));
	} else {
	  debug7(printf("C insertion of %d bp...",nindels));
      
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
  int best_plus_nmatches, best_minus_nmatches, nmatches;
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


  debug(printf("\nStarting Sarray_search_greedy with querylength %d and indexsize %d\n",querylength,sarray->indexsize));
  if (nmisses_allowed < 0) {
    nmisses_allowed = 0;
  }

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

