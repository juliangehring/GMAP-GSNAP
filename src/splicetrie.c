static char rcsid[] = "$Id: splicetrie.c 36298 2011-03-09 05:26:09Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "splicetrie.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "assert.h"
#include "mem.h"
#include "iitdef.h"
#include "interval.h"
#include "genome_hr.h"


#define MININTRONLEN 9

#define SPLICEDIST_EXTRA 10 /* Reported intron length does not include dinucleotides */


#define INTERNAL_NODE -1U
#define EMPTY_POINTER 0U

#define MAX_DUPLICATES 1000
#define DUPLICATE_NODE -1000U	/* Needs to be -MAX_DUPLICATES */

#if 0
#define MAX_BEST_NMISMATCHES 1
#endif

#define single_leaf_p(x) (x) < DUPLICATE_NODE
#define multiple_leaf_p(x) (x) < INTERNAL_NODE

/* #define USE_2BYTE_RELOFFSETS 1 */
#ifdef USE_2BYTE_RELOFFSETS
#define NULL_POINTER 65535
#else
#define NULL_POINTER -1U    /* Note: 0 does not work */
#endif


/* Generating trie */
#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Generating trie, details */
#ifdef DEBUG0
#define debug0(x) x
#else
#define debug0(x)
#endif


/* Finding short-overlap splicing.  Also may want to turn on DEBUG4H in stage1hr.c. */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif


/* Finding short exons.  Also may want to turn on DEBUG4H in stage1hr.c. */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif




/*               87654321 */
#define LEFT_A 0x00000000
#define LEFT_C 0x40000000
#define LEFT_G 0x80000000
#define LEFT_T 0xC0000000
#define HIGH2  0xC0000000

#define RIGHT_A 0x00
#define RIGHT_C 0x01
#define RIGHT_G 0x02
#define RIGHT_T 0x03


/* Puts leftmost character into lowest bits */
/* For right splicestrings, we want the leftmost character in the highest bits */

static UINT4
compress16 (bool *saw_n_p, char *buffer) {
  UINT4 low = 0U;
  int c;
  int i;

  /* *saw_n_p = false; -- Want to check both ref and alt, so rely on caller to set */
  for (i = 0; i < 16; i++) {
    c = buffer[i];
    low >>= 2;
    switch (c) {
    case 'A': break;
    case 'C': low |= LEFT_C; break;
    case 'G': low |= LEFT_G; break;
    case 'T': low |= LEFT_T; break;
    default: *saw_n_p = true; break;
    }
  }

  return low;
}

static UINT4
uint4_reverse (UINT4 forward) {
  UINT4 reverse = 0U;
  int c;
  int i;

  for (i = 0; i < 16; i++) {
    c = forward & 0x03;
    reverse <<= 2;
    reverse |= c;
    forward >>= 2;
  }

  return reverse;
}



typedef struct Splicestring_T *Splicestring_T;
struct Splicestring_T {
  UINT4 string;
  UINT4 splicesite;
  UINT4 splicesite_i;
};


static void
Splicestring_free (Splicestring_T *old) {
  FREE(*old);
  return;
}


void
Splicestring_gc (List_T *splicestrings, int nsplicesites) {
  int i;
  List_T list, p;
  Splicestring_T splicestring;

  for (i = 0; i < nsplicesites; i++) {
    list = splicestrings[i];
    for (p = list; p != NULL; p = List_next(p)) {
      splicestring = (Splicestring_T) List_head(p);
      Splicestring_free(&splicestring);
    }
    List_free(&list);
  }
  FREE(splicestrings);
  return;
}


static Splicestring_T
Splicestring_new (UINT4 string, UINT4 splicesite, int splicesite_i) {
  Splicestring_T new = (Splicestring_T) MALLOC(sizeof(*new));

  new->string = string;
  new->splicesite = splicesite;
  new->splicesite_i = (UINT4) splicesite_i;
  return new;
}

static int
Splicestring_cmp (const void *a, const void *b) {
  Splicestring_T x = * (Splicestring_T *) a;
  Splicestring_T y = * (Splicestring_T *) b;

#ifdef USE_STRINGS
  return strcmp(x->string,y->string);
#else
  if (x->string < y->string) {
    return -1;
  } else if (x->string > y->string) {
    return +1;
  } else {
    return 0;
  }
#endif
}



static List_T
allelic_combinations (UINT4 refstring, UINT4 altstring, UINT4 splicesite, int splicesite_i) {
  List_T splicestrings = NULL;
  Uintlist_T combinations, newcombinations, temp, p;
  int refc, altc;
  int i;

  combinations = Uintlist_push(NULL,0U);
  for (i = 0; i < 16; i++) {
    newcombinations = NULL;

    refc = (refstring & HIGH2) >> 30;
    altc = (altstring & HIGH2) >> 30;
    refstring <<= 2;
    altstring <<= 2;
    if (refc != altc) {
      for (p = combinations; p != NULL; p = Uintlist_next(p)) {
	newcombinations = Uintlist_push(newcombinations,(Uintlist_head(p) << 2) | refc);
	newcombinations = Uintlist_push(newcombinations,(Uintlist_head(p) << 2) | altc);
      }
    } else {
      for (p = combinations; p != NULL; p = Uintlist_next(p)) {
	newcombinations = Uintlist_push(newcombinations,(Uintlist_head(p) << 2) | refc);
      }
    }

    temp = combinations;
    combinations = Uintlist_reverse(newcombinations);
    Uintlist_free(&temp);
  }

  for (p = combinations; p != NULL; p = Uintlist_next(p)) {
    splicestrings = List_push(splicestrings,(void *) Splicestring_new(Uintlist_head(p),splicesite,splicesite_i));
  }

  Uintlist_free(&combinations);
  return splicestrings;
}




/* If distances are provided, splicedists are the observed distances */
Genomicpos_T *
Splicetrie_retrieve_splicesites (bool *distances_observed_p, Splicetype_T **splicetypes, Genomicpos_T **splicedists,
				 List_T **splicestrings, UINT4 **splicefrags_ref, UINT4 **splicefrags_alt,
				 int *nsplicesites, IIT_T splicesites_iit, int *splicesites_divint_crosstable,
				 int donor_typeint, int acceptor_typeint, IIT_T chromosome_iit,
				 Genome_T genome, Genome_T genomealt, Genomicpos_T shortsplicedist) {
  Genomicpos_T *splicesites, chrlength, chroffset, position, chrpos;
  Genomicpos_T last_donor, last_antidonor, last_acceptor, last_antiacceptor;
  int last_donor_k, last_antidonor_k, last_acceptor_k, last_antiacceptor_k;
  Genomicpos_T distance;
  UINT4 refstring, altstring;
  int *splicesites1;
  int divno, nsplicesites1, i, k;
  Chrnum_T chrnum;
  Interval_T *intervals, interval;
  char gbuffer_ref[17], gbuffer_alt[17], *chr;
  char *restofheader, *annot;
  bool firstp = true, saw_n_p, allocp, alloc_header_p;

  k = 0;
  for (chrnum = 1; chrnum <= IIT_total_nintervals(chromosome_iit); chrnum++) {
    if ((divno = splicesites_divint_crosstable[chrnum]) > 0) {
      chroffset = IIT_interval_low(chromosome_iit,chrnum);
      chrlength = IIT_length(chromosome_iit,chrnum);
      splicesites1 = IIT_get_with_divno(&nsplicesites1,splicesites_iit,divno,
					0U,chrlength-1U,/*sortp*/false);
      if (nsplicesites1 > 0) {
	if (firstp == true) {
	  annot = IIT_annotation(&restofheader,splicesites_iit,splicesites1[0],&alloc_header_p);
	  if (restofheader[0] == '\0') {
	    fprintf(stderr,"splice distances absent...");
	    *distances_observed_p = false;
	  } else if (sscanf(restofheader,"%u",&distance) < 1) {
	    fprintf(stderr,"splice distances absent...");
	    *distances_observed_p = false;
	  } else {
	    fprintf(stderr,"splice distances present...");
	    *distances_observed_p = true;
	  }
	  if (alloc_header_p == true) {
	    FREE(restofheader);
	  }
	  firstp = false;
	}

	intervals = (Interval_T *) CALLOC(nsplicesites1,sizeof(Interval_T));
	for (i = 0; i < nsplicesites1; i++) {
	  intervals[i] = &(splicesites_iit->intervals[divno][i]);
	}
	qsort(intervals,nsplicesites1,sizeof(Interval_T),Interval_cmp_low);

	last_donor = last_antidonor = last_acceptor = last_antiacceptor = 0U;
	for (i = 0; i < nsplicesites1; i++) {
	  interval = intervals[i];
	  position = Interval_low(intervals[i]) + chroffset;
	  if (Interval_type(interval) == donor_typeint) {
	    if (Interval_sign(interval) > 0) {
	      if (position != last_donor) {
		last_donor = position;
		k++;
	      }
	    } else {
	      if (position != last_antidonor) {
		last_antidonor = position;
		k++;
	      }
	    }
	  } else if (Interval_type(interval) == acceptor_typeint) {
	    if (Interval_sign(interval) > 0) {
	      if (position != last_acceptor) {
		last_acceptor = position;
		k++;
	      }
	    } else {
	      if (position != last_antiacceptor) {
		last_antiacceptor = position;
		k++;
	      }
	    }
	  }
	}
	FREE(intervals);
	FREE(splicesites1);
      }
    }
  }

  *nsplicesites = k;
  debug(printf("total unique splicesites: %d\n",*nsplicesites));
  fprintf(stderr,"%d unique splicesites...",*nsplicesites);

  /* The above procedure determines number of unique splicesites */

  if (*nsplicesites == 0) {
    splicesites = (Genomicpos_T *) NULL;
    *splicetypes = (Splicetype_T *) NULL;
    *splicedists = (Genomicpos_T *) NULL;
    *splicestrings = (List_T *) NULL;
    *splicefrags_ref = (UINT4 *) NULL;
    *splicefrags_alt = (UINT4 *) NULL;
  } else {
    splicesites = (Genomicpos_T *) CALLOC(*nsplicesites,sizeof(Genomicpos_T));
    *splicetypes = (Splicetype_T *) CALLOC(*nsplicesites,sizeof(Splicetype_T));
    *splicedists = (Genomicpos_T *) CALLOC(*nsplicesites,sizeof(Genomicpos_T));
    *splicestrings = (List_T *) CALLOC(*nsplicesites,sizeof(List_T));
    *splicefrags_ref = (UINT4 *) CALLOC(*nsplicesites,sizeof(UINT4));
    if (genomealt == NULL) {
      *splicefrags_alt = *splicefrags_ref;
    } else {
      *splicefrags_alt = (UINT4 *) CALLOC(*nsplicesites,sizeof(UINT4));
    }
  }

  k = 0;
  for (chrnum = 1; chrnum <= IIT_total_nintervals(chromosome_iit); chrnum++) {
    if ((divno = splicesites_divint_crosstable[chrnum]) > 0) {
      chroffset = IIT_interval_low(chromosome_iit,chrnum);
      chrlength = IIT_length(chromosome_iit,chrnum);
      splicesites1 = IIT_get_with_divno(&nsplicesites1,splicesites_iit,divno,
					0U,chrlength-1U,/*sortp*/false);
      if (nsplicesites1 > 0) {
	chr = IIT_label(chromosome_iit,chrnum,&allocp);
	intervals = (Interval_T *) CALLOC(nsplicesites1,sizeof(Interval_T));
	for (i = 0; i < nsplicesites1; i++) {
	  /* intervals[i] = &(splicesites_iit->intervals[divno][i]); */
	  /* Copy so we can store distance information in Interval_high */
	  intervals[i] = Interval_copy(&(splicesites_iit->intervals[divno][i]));
	  if (*distances_observed_p == false) {
	    /* No, want to have essentially zero distance */
	    /* Interval_store_length(intervals[i],shortsplicedist); */
	  } else {
	    annot = IIT_annotation(&restofheader,splicesites_iit,splicesites1[i],&alloc_header_p);
	    if (sscanf(restofheader,"%u",&distance) != 1) {
	      fprintf(stderr,"splicesites file missing distance in an entry\n");
	      exit(9);
	    } else {
	      Interval_store_length(intervals[i],distance + SPLICEDIST_EXTRA);
	    }
	    if (alloc_header_p == true) {
	      FREE(restofheader);
	    }
	  }
	}

	qsort(intervals,nsplicesites1,sizeof(Interval_T),Interval_cmp_low);

	last_donor = last_antidonor = last_acceptor = last_antiacceptor = 0U;
	for (i = 0; i < nsplicesites1; i++) {
	  interval = intervals[i];
	  chrpos = Interval_low(intervals[i]);
	  position = chrpos + chroffset;

	  if (Interval_type(interval) == donor_typeint) {
	    if (Interval_sign(interval) > 0) {
	      if (position == last_donor) {
		if (Interval_length(interval) > (*splicedists)[last_donor_k]) {
		  (*splicedists)[last_donor_k] = Interval_length(interval);
		}

	      } else {
		last_donor_k = k;
		last_donor = splicesites[k] = position;
		(*splicetypes)[k] = DONOR;
		(*splicedists)[k] = Interval_length(interval);

		saw_n_p = false;
		Genome_fill_buffer_simple(genome,position-16,16,gbuffer_ref);
		refstring = (*splicefrags_ref)[k] = compress16(&saw_n_p,gbuffer_ref);
		if (genomealt) {
		  Genome_fill_buffer_simple_alt(genomealt,position-16,16,gbuffer_alt);
		  altstring = (*splicefrags_alt)[k] = compress16(&saw_n_p,gbuffer_alt);
		}

		if (saw_n_p == true) {
		  fprintf(stderr,"\nNon-standard nucleotide N near splice site %s:%u.  Discarding...\n",
			  chr,chrpos);
		} else if (genomealt) {
		  (*splicestrings)[k] = allelic_combinations(refstring,altstring,position,k);
		  k++;
		} else {
		  (*splicestrings)[k] = List_push(NULL,(void *) Splicestring_new(refstring,position,k));
		  k++;
		}
	      }

	    } else {
	      if (position == last_antidonor) {
		if (Interval_length(interval) > (*splicedists)[last_antidonor_k]) {
		  (*splicedists)[last_antidonor_k] = Interval_length(interval);
		}

	      } else {
		last_antidonor_k = k;
		last_antidonor = splicesites[k] = position;
		(*splicetypes)[k] = ANTIDONOR;
		(*splicedists)[k] = Interval_length(interval);

		saw_n_p = false;
		Genome_fill_buffer_simple(genome,position,16,gbuffer_ref);
		refstring = (*splicefrags_ref)[k] = compress16(&saw_n_p,gbuffer_ref);
		refstring = uint4_reverse(refstring);
		if (genomealt) {
		  Genome_fill_buffer_simple_alt(genomealt,position,16,gbuffer_alt);
		  altstring = (*splicefrags_alt)[k] = compress16(&saw_n_p,gbuffer_alt);
		  altstring = uint4_reverse(altstring);
		}

		if (saw_n_p == true) {
		  fprintf(stderr,"\nNon-standard nucleotide N near splice site %s:%u.  Discarding...\n",
			  chr,chrpos);
		} else if (genomealt) {
		  (*splicestrings)[k] = allelic_combinations(refstring,altstring,position,k);
		  k++;
		} else {
		  (*splicestrings)[k] = List_push(NULL,(void *) Splicestring_new(refstring,position,k));
		  k++;
		}
	      }
	    }

	  } else if (Interval_type(interval) == acceptor_typeint) {
	    if (Interval_sign(interval) > 0) {
	      if (position == last_acceptor) {
		if (Interval_length(interval) > (*splicedists)[last_acceptor_k]) {
		  (*splicedists)[last_acceptor_k] = Interval_length(interval);
		}

	      } else {
		last_acceptor_k = k;
		last_acceptor = splicesites[k] = position;
		(*splicetypes)[k] = ACCEPTOR;
		(*splicedists)[k] = Interval_length(interval);

		saw_n_p = false;
		Genome_fill_buffer_simple(genome,position,16,gbuffer_ref);
		refstring = (*splicefrags_ref)[k] = compress16(&saw_n_p,gbuffer_ref);
		refstring = uint4_reverse(refstring);
		if (genomealt) {
		  Genome_fill_buffer_simple_alt(genomealt,position,16,gbuffer_alt);
		  altstring = (*splicefrags_alt)[k] = compress16(&saw_n_p,gbuffer_alt);
		  altstring = uint4_reverse(altstring);
		}

		if (saw_n_p == true) {
		  fprintf(stderr,"\nNon-standard nucleotide N near splice site %s:%u.  Discarding...\n",
			  chr,chrpos);
		} else if (genomealt) {
		  (*splicestrings)[k] = allelic_combinations(refstring,altstring,position,k);
		  k++;
		} else {
		  (*splicestrings)[k] = List_push(NULL,(void *) Splicestring_new(refstring,position,k));
		  k++;
		}
	      }

	    } else {
	      if (position == last_antiacceptor) {
		if (Interval_length(interval) > (*splicedists)[last_antiacceptor_k]) {
		  (*splicedists)[last_antiacceptor_k] = Interval_length(interval);
		}

	      } else {
		last_antiacceptor_k = k;
		last_antiacceptor = splicesites[k] = position;
		(*splicetypes)[k] = ANTIACCEPTOR;
		(*splicedists)[k] = Interval_length(interval);

		saw_n_p = false;
		Genome_fill_buffer_simple(genome,position-16,16,gbuffer_ref);
		refstring = (*splicefrags_ref)[k] = compress16(&saw_n_p,gbuffer_ref);
		if (genomealt) {
		  Genome_fill_buffer_simple_alt(genomealt,position-16,16,gbuffer_alt);
		  altstring = (*splicefrags_alt)[k] = compress16(&saw_n_p,gbuffer_alt);
		}

		if (saw_n_p == true) {
		  fprintf(stderr,"\nNon-standard nucleotide N near splice site %s:%u.  Discarding...\n",
			  chr,chrpos);
		} else if (genomealt) {
		  (*splicestrings)[k] = allelic_combinations(refstring,altstring,position,k);
		  k++;
		} else {
		  (*splicestrings)[k] = List_push(NULL,(void *) Splicestring_new(refstring,position,k));
		  k++;
		}
	      }

	    }
	  }
	}

	for (i = 0; i < nsplicesites1; i++) {
	  Interval_free(&(intervals[i]));
	}
	FREE(intervals);
	FREE(splicesites1);
	if (allocp == true) {
	  FREE(chr);
	}
      }
    }
  }

  *nsplicesites = k;
  fprintf(stderr,"%d splicesites without Ns...",*nsplicesites);

  return splicesites;
}

/************************************************************************/



typedef struct Trie_T *Trie_T;
struct Trie_T {
  Splicestring_T leaf;

  int nsites;
  int na;
  int nc;
  int ng;
  int nt;

  Trie_T triea;
  Trie_T triec;
  Trie_T trieg;
  Trie_T triet;
};


#if 0
static void
Trie_free (Trie_T *old) {
  if ((*old) == NULL) {
    return;
  } else if ((*old)->leaf != NULL) {
    FREE(*old);
    return;
  } else {
    Trie_free(&(*old)->triea);
    Trie_free(&(*old)->triec);
    Trie_free(&(*old)->trieg);
    Trie_free(&(*old)->triet);
    FREE(*old);
    return;
  }
}
#endif

/************************************************************************
 *   Building splicetries
 ************************************************************************/

#if 0
static Trie_T
Trie_new (Splicestring_T *sites, int nsites, int charpos) {
  Trie_T trie;
  int i, ptr;

  if (nsites == 0) {
    return (Trie_T) NULL;

  } else {

    if (nsites == 1) {
      trie = (Trie_T) MALLOC(sizeof(*trie));
      trie->nsites = 1;
      trie->na = trie->nc = trie->ng = trie->nt = 0;
      trie->triea = trie->triec = trie->trieg = trie->triet = (Trie_T) NULL;
      trie->leaf = sites[0];
      /* trie->remainder = &(sites[0]->string[charpos]); */
      return trie;

    } else if (charpos >= 16) {
      return (Trie_T) NULL;

    } else {
      trie = (Trie_T) MALLOC(sizeof(*trie));
      trie->nsites = nsites;
      trie->na = trie->nc = trie->ng = trie->nt = 0;
      trie->leaf = (Splicestring_T) NULL;
      /* trie->remainder = (char *) NULL; */

      for (i = 0; i < nsites; i++) {
	switch ((sites[i]->string >> (30 - 2*charpos)) & 0x03) {
	case RIGHT_A: trie->na++; break;
	case RIGHT_C: trie->nc++; break;
	case RIGHT_G: trie->ng++; break;
	case RIGHT_T: trie->nt++; break;
	default: abort();
	}
      }

      ptr = 0;
      trie->triea = Trie_new(&(sites[ptr]),trie->na,charpos+1);
      ptr += trie->na;

      trie->triec = Trie_new(&(sites[ptr]),trie->nc,charpos+1);
      ptr += trie->nc;

      trie->trieg = Trie_new(&(sites[ptr]),trie->ng,charpos+1);
      ptr += trie->ng;

      trie->triet = Trie_new(&(sites[ptr]),trie->nt,charpos+1);
      /* ptr += trie->nt; */

      return trie;
    }
  }
}
#endif


/* Combination of Trie_new and Trie_output */
static Uintlist_T
Trie_output_new (unsigned int *ptr, int *nprinted, Uintlist_T triecontents_list,
		 Splicestring_T *sites, int nsites, int charpos) {
  int i, k;
  unsigned int posa, posc, posg, post;
  int na, nc, ng, nt;

  if (nsites == 0) {
    debug0(printf("nsites == 0, so NULL\n"));
    *ptr = NULL_POINTER;

  } else if (nsites == 1) {
    debug0(printf("nsites == 1, so pushing %d\n",sites[0]->splicesite_i));
    *ptr = *nprinted;
    triecontents_list = Uintlist_push(triecontents_list,sites[0]->splicesite_i);
    *nprinted += 1;
    
  } else if (charpos >= 16) {
    if (nsites > MAX_DUPLICATES) {
      fprintf(stderr,"Warning: Splicetrie exceeded max duplicates value of %d\n",MAX_DUPLICATES);
      *ptr = NULL_POINTER;
    } else {
      *ptr = *nprinted;

      triecontents_list = Uintlist_push(triecontents_list,(unsigned int) (-nsites));
      *nprinted += 1;
      for (i = 0; i < nsites; i++) {
	debug0(printf(" %d",sites[i]->splicesite_i));
	triecontents_list = Uintlist_push(triecontents_list,sites[i]->splicesite_i);
	*nprinted += 1;
      }
      debug0(printf("\n"));
    }
      
  } else {
    na = nc = ng = nt = 0;
    for (i = 0; i < nsites; i++) {
      switch ((sites[i]->string >> (30 - 2*charpos)) & 0x03) {
      case RIGHT_A: na++; break;
      case RIGHT_C: nc++; break;
      case RIGHT_G: ng++; break;
      case RIGHT_T: nt++; break;
      default: abort();
      }
    }
    debug0(printf("%d A, %d C, %d G, %d T\n",na,nc,ng,nt));

    k = 0;
    triecontents_list = Trie_output_new(&posa,&(*nprinted),triecontents_list,
					&(sites[k]),na,charpos+1);
    k += na;
    triecontents_list = Trie_output_new(&posc,&(*nprinted),triecontents_list,
					&(sites[k]),nc,charpos+1);
    k += nc;
    triecontents_list = Trie_output_new(&posg,&(*nprinted),triecontents_list,
					&(sites[k]),ng,charpos+1);
    k += ng;
    triecontents_list = Trie_output_new(&post,&(*nprinted),triecontents_list,
					&(sites[k]),nt,charpos+1);
    /* k += nt; */

    *ptr = *nprinted;
    triecontents_list = Uintlist_push(triecontents_list,INTERNAL_NODE);

    if (posa == NULL_POINTER) {
      triecontents_list = Uintlist_push(triecontents_list,EMPTY_POINTER);
    } else {
      triecontents_list = Uintlist_push(triecontents_list,*ptr - posa);
    }
    if (posc == NULL_POINTER) {
      triecontents_list = Uintlist_push(triecontents_list,EMPTY_POINTER);
    } else {
      triecontents_list = Uintlist_push(triecontents_list,*ptr - posc);
    }
    if (posg == NULL_POINTER) {
      triecontents_list = Uintlist_push(triecontents_list,EMPTY_POINTER);
    } else {
      triecontents_list = Uintlist_push(triecontents_list,*ptr - posg);
    }
    if (post == NULL_POINTER) {
      triecontents_list = Uintlist_push(triecontents_list,EMPTY_POINTER);
    } else {
      triecontents_list = Uintlist_push(triecontents_list,*ptr - post);
    }
    *nprinted += 5;
  }

  return triecontents_list;
}


#if 0
static void
Trie_print (Trie_T trie, int level) {
  int i;

  if (trie == NULL) {
    printf("\n");

  } else {
    if (trie->leaf != NULL) {
      printf(" %d@%u\n",(int) trie->leaf->splicesite_i,trie->leaf->splicesite);

    } else {
      printf("\n");
      for (i = 0; i < level; i++) {
	printf("\t");
      }
      printf("A (%d)",trie->na);
      Trie_print(trie->triea,level+1);

      for (i = 0; i < level; i++) {
	printf("\t");
      }
      printf("C (%d)",trie->nc);
      Trie_print(trie->triec,level+1);

      for (i = 0; i < level; i++) {
	printf("\t");
      }
      printf("G (%d)",trie->ng);
      Trie_print(trie->trieg,level+1);

      for (i = 0; i < level; i++) {
	printf("\t");
      }
      printf("T (%d)",trie->nt);
      Trie_print(trie->triet,level+1);
    }
  }

  return;
}
#endif

#if 0
static int
Trie_print_compact (int *ptr, FILE *fp, Trie_T trie, int nprinted) {
  int posa, posc, posg, post;

  if (trie == NULL) {
    *ptr = NULL_POINTER;
    return nprinted;

  } else {
    if (trie->leaf != NULL) {
      printf("%d. %d@%u\n",nprinted,(int) trie->leaf->splicesite_i,trie->leaf->splicesite);
      *ptr = nprinted;
      return nprinted + 1;

    } else {
      nprinted = Trie_print_compact(&posa,fp,trie->triea,nprinted);
      nprinted = Trie_print_compact(&posc,fp,trie->triec,nprinted);
      nprinted = Trie_print_compact(&posg,fp,trie->trieg,nprinted);
      nprinted = Trie_print_compact(&post,fp,trie->triet,nprinted);

      *ptr = nprinted;

      printf("%d. %u\n",nprinted,INTERNAL_NODE);

      printf("relative: %d %d %d %d\n",
	     *ptr - posa,*ptr - posc,*ptr - posg,*ptr - post);
      printf("absolute: %d %d %d %d\n",
	     posa,posc,posg,post);

      return nprinted + 3;
    }
  }
}
#endif


static Uintlist_T
Trie_output_empty (unsigned int *ptr, int *nprinted, Uintlist_T triecontents_list) {
  *ptr = (unsigned int) *nprinted;
  triecontents_list = Uintlist_push(triecontents_list,INTERNAL_NODE);
#ifdef USE_2BYTE_RELOFFSETS
  triecontents_list = Uintlist_push(triecontents_list,EMPTY_POINTER);
  triecontents_list = Uintlist_push(triecontents_list,EMPTY_POINTER);
  *nprinted += 3;
#else
  triecontents_list = Uintlist_push(triecontents_list,EMPTY_POINTER);
  triecontents_list = Uintlist_push(triecontents_list,EMPTY_POINTER);
  triecontents_list = Uintlist_push(triecontents_list,EMPTY_POINTER);
  triecontents_list = Uintlist_push(triecontents_list,EMPTY_POINTER);
  *nprinted += 5;
#endif
  return triecontents_list;
}


#if 0
static Uintlist_T
Trie_output (int *ptr, int *nprinted, Uintlist_T triecontents_list, Trie_T trie) {
  int posa, posc, posg, post;
#ifdef USE_2BYTE_RELOFFSETS
  unsigned int pointersAC, pointersGT;
  int reloffset;
#endif

  if (trie == NULL) {
    *ptr = NULL_POINTER;

  } else if (trie->leaf != NULL) {
    *ptr = *nprinted;
    /* Entering splicesite_i, rather than position */
    triecontents_list = Uintlist_push(triecontents_list,trie->leaf->splicesite_i);
    *nprinted += 1;

  } else {
    triecontents_list = Trie_output(&posa,&(*nprinted),triecontents_list,trie->triea);
    triecontents_list = Trie_output(&posc,&(*nprinted),triecontents_list,trie->triec);
    triecontents_list = Trie_output(&posg,&(*nprinted),triecontents_list,trie->trieg);
    triecontents_list = Trie_output(&post,&(*nprinted),triecontents_list,trie->triet);

    *ptr = *nprinted;
    triecontents_list = Uintlist_push(triecontents_list,INTERNAL_NODE);


#ifdef USE_2BYTE_RELOFFSETS
    if (posa == NULL_POINTER) {
      pointersAC = 0;
    } else if ((reloffset = *ptr - posa) >= 65536) {
      fprintf(stderr,"Reloffset A %d = %d - %d is too big for 2 bytes\n",reloffset,*ptr,posa);
      abort();
    } else {
      pointersAC = reloffset;
    }

    if (posc == NULL_POINTER) {
      pointersAC = (pointersAC << 16);
    } else if ((reloffset = *ptr - posc) >= 65536) {
      fprintf(stderr,"Reloffset C %d = %d - %d is too big for 2 bytes\n",reloffset,*ptr,posc);
      abort();
    } else {
      pointersAC = (pointersAC << 16) + reloffset;
    }

    if (posg == NULL_POINTER) {
      pointersGT = 0;
    } else if ((reloffset = *ptr - posg) >= 65536) {
      fprintf(stderr,"Reloffset G %d = %d - %d is too big for 2 bytes\n",reloffset,*ptr,posg);
      abort();
    } else {
      pointersGT = reloffset;
    }

    if (post == NULL_POINTER) {
      pointersGT = (pointersGT << 16);
    } else if ((reloffset = *ptr - post) >= 65536) {
      fprintf(stderr,"Reloffset T %d = %d - %d is too big for 2 bytes\n",reloffset,*ptr,post);
      abort();
    } else {
      pointersGT = (pointersGT << 16) + reloffset;
    }

    triecontents_list = Uintlist_push(triecontents_list,pointersAC);
    triecontents_list = Uintlist_push(triecontents_list,pointersGT);

    *nprinted += 3;
#else
    if (posa == NULL_POINTER) {
      triecontents_list = Uintlist_push(triecontents_list,EMPTY_POINTER);
    } else {
      triecontents_list = Uintlist_push(triecontents_list,*ptr - posa);
    }
    if (posc == NULL_POINTER) {
      triecontents_list = Uintlist_push(triecontents_list,EMPTY_POINTER);
    } else {
      triecontents_list = Uintlist_push(triecontents_list,*ptr - posc);
    }
    if (posg == NULL_POINTER) {
      triecontents_list = Uintlist_push(triecontents_list,EMPTY_POINTER);
    } else {
      triecontents_list = Uintlist_push(triecontents_list,*ptr - posg);
    }
    if (post == NULL_POINTER) {
      triecontents_list = Uintlist_push(triecontents_list,EMPTY_POINTER);
    } else {
      triecontents_list = Uintlist_push(triecontents_list,*ptr - post);
    }
    *nprinted += 5;
#endif

  }
  return triecontents_list;
}
#endif


void
Splicetrie_npartners (int **nsplicepartners_skip, int **nsplicepartners_obs, int **nsplicepartners_max,
		      Genomicpos_T *splicesites, Splicetype_T *splicetypes,
		      Genomicpos_T *splicedists, List_T *splicestrings, int nsplicesites,
		      IIT_T chromosome_iit, Genomicpos_T max_distance,
		      bool distances_observed_p) {
  int nsites_skip, nsites_obs, nsites_max, j, j1;
  Genomicpos_T chroffset, chrhigh, leftbound_obs, leftbound_max, rightbound_obs, rightbound_max;
  Genomicpos_T leftbound_min, rightbound_min;
  Chrnum_T chrnum;

  (*nsplicepartners_skip) = (int *) CALLOC(nsplicesites,sizeof(int));
  if (distances_observed_p == true) {
    (*nsplicepartners_obs) = (int *) CALLOC(nsplicesites,sizeof(int));
  } else {
    (*nsplicepartners_obs) = (int *) NULL;
  }
  (*nsplicepartners_max) = (int *) CALLOC(nsplicesites,sizeof(int));

  chrhigh = 0U;
  for (j = 0; j < nsplicesites; j++) {
    if (splicesites[j] > chrhigh) {
      chrnum = IIT_get_one(chromosome_iit,/*divstring*/NULL,splicesites[j],splicesites[j]);
      IIT_interval_bounds(&chroffset,&chrhigh,chromosome_iit,chrnum);
      chrhigh += 1U;
    }

    switch (splicetypes[j]) {
    case DONOR:
      nsites_skip = nsites_obs = nsites_max = 0;
      j1 = j + 1;

      if ((rightbound_min = splicesites[j] + MININTRONLEN) > chrhigh) {
	rightbound_min = chrhigh;
      }
      while (j1 < nsplicesites && splicesites[j1] < rightbound_min) {
	nsites_skip++;
	j1++;
      }

      if (distances_observed_p) {
	if ((rightbound_obs = splicesites[j] + splicedists[j]) > chrhigh) {
	  rightbound_obs = chrhigh;
	}
	while (j1 < nsplicesites && splicesites[j1] < rightbound_obs) {
	  if (splicetypes[j1] == ACCEPTOR) {
	    nsites_obs += List_length(splicestrings[j1]);
	  }
	  j1++;
	}
      }

      if ((rightbound_max = splicesites[j] + max_distance) > chrhigh) {
	rightbound_max = chrhigh;
      }
      while (j1 < nsplicesites && splicesites[j1] < rightbound_max) {
	if (splicetypes[j1] == ACCEPTOR) {
	  nsites_max += List_length(splicestrings[j1]);
	}
	j1++;
      }

      break;

    case ACCEPTOR:
      nsites_skip = nsites_obs = nsites_max = 0;
      j1 = j - 1;

      if (splicesites[j] < chroffset + MININTRONLEN) {
	leftbound_min = chroffset;
      } else {
	leftbound_min = splicesites[j] - MININTRONLEN;
      }
      while (j1 >= 0 && splicesites[j1] > leftbound_min) {
	nsites_skip++;
	j1--;
      }

      if (distances_observed_p) {
	if (splicesites[j] < chroffset + splicedists[j]) {
	  leftbound_obs = chroffset;
	} else {
	  leftbound_obs = splicesites[j] - splicedists[j];
	}
	while (j1 >= 0 && splicesites[j1] > leftbound_obs) {
	  if (splicetypes[j1] == DONOR) {
	    nsites_obs += List_length(splicestrings[j1]);
	  }
	  j1--;
	}
      }

      if (splicesites[j] < chroffset + max_distance) {
	leftbound_max = chroffset;
      } else {
	leftbound_max = splicesites[j] - max_distance;
      }
      while (j1 >= 0 && splicesites[j1] > leftbound_max) {
	if (splicetypes[j1] == DONOR) {
	  nsites_max += List_length(splicestrings[j1]);
	}
	j1--;
      }

      break;

    case ANTIDONOR:
      nsites_skip = nsites_obs = nsites_max = 0;
      j1 = j - 1;

      if (splicesites[j] < chroffset + MININTRONLEN) {
	leftbound_min = chroffset;
      } else {
	leftbound_min = splicesites[j] - MININTRONLEN;
      }
      while (j1 >= 0 && splicesites[j1] > leftbound_min) {
	nsites_skip++;
	j1--;
      }

      if (distances_observed_p) {
	if (splicesites[j] < chroffset + splicedists[j]) {
	  leftbound_obs = chroffset;
	} else {
	  leftbound_obs = splicesites[j] - splicedists[j];
	}
	while (j1 >= 0 && splicesites[j1] > leftbound_obs) {
	  if (splicetypes[j1] == ANTIACCEPTOR) {
	    nsites_obs += List_length(splicestrings[j1]);
	  }
	  j1--;
	}
      }

      if (splicesites[j] < chroffset + max_distance) {
	leftbound_max = chroffset;
      } else {
	leftbound_max = splicesites[j] - max_distance;
      }
      while (j1 >= 0 && splicesites[j1] > leftbound_max) {
	if (splicetypes[j1] == ANTIACCEPTOR) {
	  nsites_max += List_length(splicestrings[j1]);
	}
	j1--;
      }

      break;

    case ANTIACCEPTOR:
      nsites_skip = nsites_obs = nsites_max = 0;
      j1 = j + 1;

      if ((rightbound_min = splicesites[j] + MININTRONLEN) > chrhigh) {
	rightbound_min = chrhigh;
      }
      while (j1 < nsplicesites && splicesites[j1] < rightbound_min) {
	nsites_skip++;
	j1++;
      }

      if (distances_observed_p) {
	if ((rightbound_obs = splicesites[j] + splicedists[j]) > chrhigh) {
	  rightbound_obs = chrhigh;
	}
	while (j1 < nsplicesites && splicesites[j1] < rightbound_obs) {
	  if (splicetypes[j1] == ANTIDONOR) {
	    nsites_obs += List_length(splicestrings[j1]);
	  }
	  j1++;
	}
      }

      if ((rightbound_max = splicesites[j] + max_distance) > chrhigh) {
	rightbound_max = chrhigh;
      }
      while (j1 < nsplicesites && splicesites[j1] < rightbound_max) {
	if (splicetypes[j1] == ANTIDONOR) {
	  nsites_max += List_length(splicestrings[j1]);
	}
	j1++;
      }

      break;

    default: 
      fprintf(stderr,"Unexpected splicetype %d\n",splicetypes[j]);
      abort();
    }

    (*nsplicepartners_skip)[j] = nsites_skip;
    if (distances_observed_p) {
      (*nsplicepartners_obs)[j] = nsites_obs;
    }
    (*nsplicepartners_max)[j] = nsites_max;
  }

  return;
}



/************************************************************************
 *   Using splicetries
 ************************************************************************/

#ifdef USE_2BYTE_RELOFFSETS
static void
get_offsets (int *offseta, int *offsetc, int *offsetg, int *offsett,
	     unsigned int offsets1, unsigned int offsets2) {

  *offsetc = (int) (offsets1 & 0xffff);
  *offseta = (int) ((offsets1 >>= 16) & 0xffff);

  *offsett = (int) (offsets2 & 0xffff);
  *offsetg = (int) ((offsets2 >>= 16) & 0xffff);

  return;
}
#endif


#if 0
static void
make_reverse_buffered (char *reverse, char *sequence, unsigned int length) {
  int i, j;

  for (i = length-1, j = 0; i >= 0; i--, j++) {
    reverse[j] = sequence[i];
  }
  reverse[length] = '\0';
  return;
}
#endif


/*                      87654321 */
#define LOW_TWO_BITS  0x00000003

static void
splicefrag_nt_leftward (char *nt, UINT4 splicefrag) {
  int i, j;
  UINT4 lowbits;

  j = 15;
  for (i = 0; i < 16; i++) {
    lowbits = splicefrag & LOW_TWO_BITS;
    switch (lowbits) {
    case RIGHT_A: nt[j] = 'A'; break;
    case RIGHT_C: nt[j] = 'C'; break;
    case RIGHT_G: nt[j] = 'G'; break;
    case RIGHT_T: nt[j] = 'T'; break;
    }
    splicefrag >>= 2;
    j--;
  }

  return;
}

static void
splicefrag_nt_rightward (char *nt, UINT4 splicefrag) {
  int i, j;
  UINT4 lowbits;

  j = 0;
  for (i = 0; i < 16; i++) {
    lowbits = splicefrag & LOW_TWO_BITS;
    switch (lowbits) {
    case RIGHT_A: nt[j] = 'A'; break;
    case RIGHT_C: nt[j] = 'C'; break;
    case RIGHT_G: nt[j] = 'G'; break;
    case RIGHT_T: nt[j] = 'T'; break;
    }
    splicefrag >>= 2;
    j++;
  }

  return;
}



#ifdef DEBUG2
/* Splicetype is for the anchor splice */
static void
Splicetrie_dump (unsigned int *triestart, Genomicpos_T *splicesites,
		 Splicetype_T splicetype, UINT4 *splicefrags_ref) {
  unsigned int leaf;
  int nleaves, i;
  Genomicpos_T position;
  int offseta, offsetc, offsetg, offsett;
  char gbuffer[17];

  gbuffer[16] = '\0';

  if (single_leaf_p(leaf = triestart[0])) {
    position = splicesites[leaf];
    printf("%d %u",(int) leaf,position);
    if (splicetype == DONOR || splicetype == ANTIACCEPTOR) {
      splicefrag_nt_rightward(gbuffer,splicefrags_ref[leaf]);
      printf(" %s (rightward)\n",gbuffer);
    } else {
      splicefrag_nt_leftward(gbuffer,splicefrags_ref[leaf]);
      printf(" %s (leftward)\n",gbuffer);
    }

  } else if (multiple_leaf_p(leaf)) {
    nleaves = (int) (-leaf);
    for (i = 1; i <= nleaves; i++) {
      leaf = triestart[i];
      position = splicesites[leaf];
      printf("%d %u",(int) leaf,position);
      if (splicetype == DONOR || splicetype == ANTIACCEPTOR) {
	splicefrag_nt_rightward(gbuffer,splicefrags_ref[leaf]);
	printf(" %s (rightward)\n",gbuffer);
      } else {
	splicefrag_nt_leftward(gbuffer,splicefrags_ref[leaf]);
	printf(" %s (leftward)\n",gbuffer);
      }
    }

  } else {
#ifdef USE_2BYTE_RELOFFSETS
    get_offsets(&offseta,&offsetc,&offsetg,&offsett,triestart[1],triestart[2]);
#else
    offseta = (int) triestart[1];
    offsetc = (int) triestart[2];
    offsetg = (int) triestart[3];
    offsett = (int) triestart[4];
#endif

    if (offseta > 0) {
      Splicetrie_dump(&(triestart[-offseta]),splicesites,splicetype,splicefrags_ref);
    }
    if (offsetc > 0) {
      Splicetrie_dump(&(triestart[-offsetc]),splicesites,splicetype,splicefrags_ref);
    }
    if (offsetg > 0) {
      Splicetrie_dump(&(triestart[-offsetg]),splicesites,splicetype,splicefrags_ref);
    }
    if (offsett > 0) {
      Splicetrie_dump(&(triestart[-offsett]),splicesites,splicetype,splicefrags_ref);
    }
  }

  return;
}
#endif


#ifdef DEBUG
/* Splicetype is for the anchor splice */
static void
dump_sites (Splicestring_T *sites, int nsites, Splicetype_T splicetype) {
  int i;
  char gbuffer[17];

  gbuffer[16] = '\0';

  for (i = 0; i < nsites; i++) {
    printf("%d %u",sites[i]->splicesite_i,sites[i]->splicesite);
    if (splicetype == DONOR || splicetype == ANTIACCEPTOR) {
      splicefrag_nt_rightward(gbuffer,sites[i]->string);
      printf(" %s (rightward)\n",gbuffer);
    } else {
      splicefrag_nt_leftward(gbuffer,sites[i]->string);
      printf(" %s (leftward)\n",gbuffer);
    }
  }
  return;
}
#endif


void
Splicetrie_build (unsigned int **triecontents_obs, unsigned int **trieoffsets_obs,
		  unsigned int **triecontents_max, unsigned int **trieoffsets_max,
		  int *nsplicepartners_skip, int *nsplicepartners_obs, int *nsplicepartners_max,
		  Splicetype_T *splicetypes, List_T *splicestrings, int nsplicesites) {
  Uintlist_T triecontents_obs_list = NULL, triecontents_max_list = NULL;
  List_T p;
  int nsites, j, j1;
  Splicestring_T *sites;
  int nprinted_obs = 0, nprinted_max = 0;
  bool distances_observed_p;

  if (nsplicepartners_obs == NULL) {
    distances_observed_p = false;
  } else {
    distances_observed_p = true;
  }

  if (distances_observed_p == true) {
    *trieoffsets_obs = (unsigned int *) CALLOC(nsplicesites,sizeof(unsigned int));
  } else {
    *trieoffsets_obs = (unsigned int *) NULL;
  }
  *trieoffsets_max = (unsigned int *) CALLOC(nsplicesites,sizeof(unsigned int));

  for (j = 0; j < nsplicesites; j++) {
    switch (splicetypes[j]) {
    case DONOR:
      debug(
	    if (distances_observed_p == true) {
	      printf("donor #%d (%d partners obs, %d partners max):",
		     j,nsplicepartners_obs[j],nsplicepartners_max[j]);
	    } else {
	      printf("donor #%d (%d partners max):",
		     j,nsplicepartners_max[j]);
	    });

      j1 = j + 1 + nsplicepartners_skip[j];

      if (distances_observed_p) {
	if (nsplicepartners_obs[j] == 0) {
	  triecontents_obs_list = Trie_output_empty(&((*trieoffsets_obs)[j]),&nprinted_obs,triecontents_obs_list);

	} else {
	  sites = (Splicestring_T *) CALLOC(nsplicepartners_obs[j],sizeof(Splicestring_T));
	  nsites = 0;
	  while (nsites < nsplicepartners_obs[j]) {
	    if (splicetypes[j1] == ACCEPTOR) {
	      for (p = splicestrings[j1]; p != NULL; p = List_next(p)) {
		debug(printf(" %d",j1));
		sites[nsites++] = (Splicestring_T) List_head(p);
	      }
	    }
	    j1++;
	  }
	  assert(nsites == nsplicepartners_obs[j]);
	  qsort(sites,nsites,sizeof(Splicestring_T),Splicestring_cmp);
	  triecontents_obs_list = Trie_output_new(&((*trieoffsets_obs)[j]),&nprinted_obs,triecontents_obs_list,
						  sites,nsites,/*charpos*/0);
	  FREE(sites);
	}
      }

      if (nsplicepartners_max[j] == 0) {
	triecontents_max_list = Trie_output_empty(&((*trieoffsets_max)[j]),&nprinted_max,triecontents_max_list);

      } else {
	sites = (Splicestring_T *) CALLOC(nsplicepartners_max[j],sizeof(Splicestring_T));
	nsites = 0;
	while (nsites < nsplicepartners_max[j]) {
	  if (splicetypes[j1] == ACCEPTOR) {
	    for (p = splicestrings[j1]; p != NULL; p = List_next(p)) {
	      debug(printf(" %d",j1));
	      sites[nsites++] = (Splicestring_T) List_head(p);
	    }
	  }
	  j1++;
	}
	assert(nsites == nsplicepartners_max[j]);
	qsort(sites,nsites,sizeof(Splicestring_T),Splicestring_cmp);
	triecontents_max_list = Trie_output_new(&((*trieoffsets_max)[j]),&nprinted_max,triecontents_max_list,
						sites,nsites,/*charpos*/0);
	FREE(sites);
      }

      debug(printf("\n"));
      break;

    case ACCEPTOR:
      debug(
	    if (distances_observed_p == true) {
	      printf("acceptor #%d (%d partners obs, %d partners max):",
		     j,nsplicepartners_obs[j],nsplicepartners_max[j]);
	    } else {
	      printf("acceptor #%d (%d partners max):",
		     j,nsplicepartners_max[j]);
	    });

      j1 = j - 1 - nsplicepartners_skip[j];

      if (distances_observed_p) {
	if (nsplicepartners_obs[j] == 0) {
	  triecontents_obs_list = Trie_output_empty(&((*trieoffsets_obs)[j]),&nprinted_obs,triecontents_obs_list);

	} else {
	  sites = (Splicestring_T *) CALLOC(nsplicepartners_obs[j],sizeof(Splicestring_T));
	  nsites = 0;
	  while (nsites < nsplicepartners_obs[j]) {
	    if (splicetypes[j1] == DONOR) {
	      for (p = splicestrings[j1]; p != NULL; p = List_next(p)) {
		debug(printf(" %d",j1));
		sites[nsites++] = (Splicestring_T) List_head(p);
	      }
	    }
	    j1--;
	  }
	  assert(nsites == nsplicepartners_obs[j]);
	  qsort(sites,nsites,sizeof(Splicestring_T),Splicestring_cmp);
	  triecontents_obs_list = Trie_output_new(&((*trieoffsets_obs)[j]),&nprinted_obs,triecontents_obs_list,
						  sites,nsites,/*charpos*/0);
	  FREE(sites);
	}
      }

      if (nsplicepartners_max[j] == 0) {
	triecontents_max_list = Trie_output_empty(&((*trieoffsets_max)[j]),&nprinted_max,triecontents_max_list);

      } else {
	sites = (Splicestring_T *) CALLOC(nsplicepartners_max[j],sizeof(Splicestring_T));
	nsites = 0;
	while (nsites < nsplicepartners_max[j]) {
	  if (splicetypes[j1] == DONOR) {
	    for (p = splicestrings[j1]; p != NULL; p = List_next(p)) {
	      debug(printf(" %d",j1));
	      sites[nsites++] = (Splicestring_T) List_head(p);
	    }
	  }
	  j1--;
	}
	assert(nsites == nsplicepartners_max[j]);
	qsort(sites,nsites,sizeof(Splicestring_T),Splicestring_cmp);
	triecontents_max_list = Trie_output_new(&((*trieoffsets_max)[j]),&nprinted_max,triecontents_max_list,
						sites,nsites,/*charpos*/0);
	FREE(sites);
      }

      debug(printf("\n"));
      break;

    case ANTIDONOR:
      debug(
	    if (distances_observed_p == true) {
	      printf("antidonor #%d (%d partners obs, %d partners max):",
		     j,nsplicepartners_obs[j],nsplicepartners_max[j]);
	    } else {
	      printf("antidonor #%d (%d partners max):",
		     j,nsplicepartners_max[j]);
	    });

      j1 = j - 1 - nsplicepartners_skip[j];

      if (distances_observed_p) {
	if (nsplicepartners_obs[j] == 0) {
	  triecontents_obs_list = Trie_output_empty(&((*trieoffsets_obs)[j]),&nprinted_obs,triecontents_obs_list);

	} else {
	  sites = (Splicestring_T *) CALLOC(nsplicepartners_obs[j],sizeof(Splicestring_T));
	  nsites = 0;
	  while (nsites < nsplicepartners_obs[j]) {
	    if (splicetypes[j1] == ANTIACCEPTOR) {
	      for (p = splicestrings[j1]; p != NULL; p = List_next(p)) {
		debug(printf(" %d",j1));
		sites[nsites++] = (Splicestring_T) List_head(p);
	      }
	    }
	    j1--;
	  }
	  assert(nsites == nsplicepartners_obs[j]);
	  qsort(sites,nsites,sizeof(Splicestring_T),Splicestring_cmp);
	  triecontents_obs_list = Trie_output_new(&((*trieoffsets_obs)[j]),&nprinted_obs,triecontents_obs_list,
						  sites,nsites,/*charpos*/0);
	  FREE(sites);
	}
      }

      if (nsplicepartners_max[j] == 0) {
	triecontents_max_list = Trie_output_empty(&((*trieoffsets_max)[j]),&nprinted_max,triecontents_max_list);

      } else {
	sites = (Splicestring_T *) CALLOC(nsplicepartners_max[j],sizeof(Splicestring_T));
	nsites = 0;
	while (nsites < nsplicepartners_max[j]) {
	  if (splicetypes[j1] == ANTIACCEPTOR) {
	    for (p = splicestrings[j1]; p != NULL; p = List_next(p)) {
	      debug(printf(" %d",j1));
	      sites[nsites++] = (Splicestring_T) List_head(p);
	    }
	  }
	  j1--;
	}
	assert(nsites == nsplicepartners_max[j]);
	qsort(sites,nsites,sizeof(Splicestring_T),Splicestring_cmp);
	triecontents_max_list = Trie_output_new(&((*trieoffsets_max)[j]),&nprinted_max,triecontents_max_list,
						sites,nsites,/*charpos*/0);
	FREE(sites);
      }

      debug(printf("\n"));
      break;

    case ANTIACCEPTOR:
      debug(
	    if (distances_observed_p == true) {
	      printf("antiacceptor #%d (%d partners obs, %d partners max):",
		     j,nsplicepartners_obs[j],nsplicepartners_max[j]);
	    } else {
	      printf("antiacceptor #%d (%d partners max):",
		     j,nsplicepartners_max[j]);
	    });

      j1 = j + 1 + nsplicepartners_skip[j];

      if (distances_observed_p) {
	if (nsplicepartners_obs[j] == 0) {
	  triecontents_obs_list = Trie_output_empty(&((*trieoffsets_obs)[j]),&nprinted_obs,triecontents_obs_list);

	} else {
	  sites = (Splicestring_T *) CALLOC(nsplicepartners_obs[j],sizeof(Splicestring_T));
	  nsites = 0;
	  while (nsites < nsplicepartners_obs[j]) {
	    if (splicetypes[j1] == ANTIDONOR) {
	      for (p = splicestrings[j1]; p != NULL; p = List_next(p)) {
		debug(printf(" %d",j1));
		sites[nsites++] = (Splicestring_T) List_head(p);
	      }
	    }
	    j1++;
	  }
	  assert(nsites == nsplicepartners_obs[j]);
	  qsort(sites,nsites,sizeof(Splicestring_T),Splicestring_cmp);
	  triecontents_obs_list = Trie_output_new(&((*trieoffsets_obs)[j]),&nprinted_obs,triecontents_obs_list,
						  sites,nsites,/*charpos*/0);
	  FREE(sites);
	}
      }

      if (nsplicepartners_max[j] == 0) {
	triecontents_max_list = Trie_output_empty(&((*trieoffsets_max)[j]),&nprinted_max,triecontents_max_list);

      } else {
	sites = (Splicestring_T *) CALLOC(nsplicepartners_max[j],sizeof(Splicestring_T));
	nsites = 0;
	while (nsites < nsplicepartners_max[j]) {
	  if (splicetypes[j1] == ANTIDONOR) {
	    for (p = splicestrings[j1]; p != NULL; p = List_next(p)) {
	      debug(printf(" %d",j1));
	      sites[nsites++] = (Splicestring_T) List_head(p);
	    }
	  }
	  j1++;
	}
	assert(nsites == nsplicepartners_max[j]);
	qsort(sites,nsites,sizeof(Splicestring_T),Splicestring_cmp);
	triecontents_max_list = Trie_output_new(&((*trieoffsets_max)[j]),&nprinted_max,triecontents_max_list,
						sites,nsites,/*charpos*/0);
	FREE(sites);
      }

      debug(printf("\n"));
      break;

    default:
      fprintf(stderr,"Unexpected splicetype %d\n",splicetypes[j]);
      abort();
    }
  }

  if (distances_observed_p) {
    fprintf(stderr,"splicetrie_obs has %d entries...",nprinted_obs);
    triecontents_obs_list = Uintlist_reverse(triecontents_obs_list);
    *triecontents_obs = Uintlist_to_array(&nprinted_obs,triecontents_obs_list);
    Uintlist_free(&triecontents_obs_list);
  } else {
    *triecontents_obs = (unsigned int *) NULL;
  }

  fprintf(stderr,"splicetrie_max has %d entries...",nprinted_max);
  triecontents_max_list = Uintlist_reverse(triecontents_max_list);
  *triecontents_max = Uintlist_to_array(&nprinted_max,triecontents_max_list);
  Uintlist_free(&triecontents_max_list);

  return;
}


static void
Splicetrie_build_one (unsigned int **triecontents_obs, unsigned int **triestart_obs,
		      unsigned int **triecontents_max, unsigned int **triestart_max,
		      int *nsplicepartners_skip, int *nsplicepartners_obs, int *nsplicepartners_max, 
		      int j, Splicetype_T *splicetypes, List_T *splicestrings) {
  Uintlist_T triecontents_obs_list = NULL, triecontents_max_list = NULL;
  List_T p;
  int nsites, j1;
  Splicestring_T *sites;
  int nprinted_obs = 0, nprinted_max = 0;
  unsigned int ptr_obs, ptr_max;
  bool distances_observed_p;
#ifdef DEBUG
  Splicestring_T splicestring;
  char gbuffer[17];
#endif

  if (nsplicepartners_obs == NULL) {
    distances_observed_p = false;
  } else {
    distances_observed_p = true;
  }


  debug(gbuffer[16] = '\0');

  switch (splicetypes[j]) {
  case DONOR:
    debug(
	  if (distances_observed_p == true) {
	    printf("donor #%d (%d partners obs, %d partners max):",
		   j,nsplicepartners_obs[j],nsplicepartners_max[j]);
	  } else {
	    printf("donor #%d (%d partners max):",
		   j,nsplicepartners_max[j]);
	  });

    j1 = j + 1 + nsplicepartners_skip[j];

    if (distances_observed_p) {
      if (nsplicepartners_obs[j] == 0) {
	triecontents_obs_list = Trie_output_empty(&ptr_obs,&nprinted_obs,triecontents_obs_list);
      
      } else {
	sites = (Splicestring_T *) CALLOC(nsplicepartners_obs[j],sizeof(Splicestring_T));
	nsites = 0;
	while (nsites < nsplicepartners_obs[j]) {
	  if (splicetypes[j1] == ACCEPTOR) {
	    for (p = splicestrings[j1]; p != NULL; p = List_next(p)) {
	      debug(splicestring = (Splicestring_T) List_head(p));
	      debug(splicefrag_nt_rightward(gbuffer,splicestring->string));
	      debug(printf(" %d (%s)",j1,gbuffer));
	      sites[nsites++] = (Splicestring_T) List_head(p);
	    }
	  }
	  j1++;
	}
	assert(nsites == nsplicepartners_obs[j]);
	qsort(sites,nsites,sizeof(Splicestring_T),Splicestring_cmp);
	debug(printf("\n"));
	debug(dump_sites(sites,nsites,splicetypes[j]));

	triecontents_obs_list = Trie_output_new(&ptr_obs,&nprinted_obs,triecontents_obs_list,
						sites,nsites,/*charpos*/0);
	FREE(sites);
      }
    }

    if (nsplicepartners_max[j] == 0) {
      triecontents_max_list = Trie_output_empty(&ptr_max,&nprinted_max,triecontents_max_list);
      
    } else {
      sites = (Splicestring_T *) CALLOC(nsplicepartners_max[j],sizeof(Splicestring_T));
      nsites = 0;
      while (nsites < nsplicepartners_max[j]) {
	if (splicetypes[j1] == ACCEPTOR) {
	  for (p = splicestrings[j1]; p != NULL; p = List_next(p)) {
	    debug(splicestring = (Splicestring_T) List_head(p));
	    debug(splicefrag_nt_rightward(gbuffer,splicestring->string));
	    debug(printf(" %d (%s)",j1,gbuffer));
	    sites[nsites++] = (Splicestring_T) List_head(p);
	  }
	}
	j1++;
      }
      assert(nsites == nsplicepartners_max[j]);
      qsort(sites,nsites,sizeof(Splicestring_T),Splicestring_cmp);
      debug(printf("\n"));
      debug(dump_sites(sites,nsites,splicetypes[j]));

      triecontents_max_list = Trie_output_new(&ptr_max,&nprinted_max,triecontents_max_list,
					      sites,nsites,/*charpos*/0);
      FREE(sites);
    }

    debug(printf("\n"));
    break;

  case ACCEPTOR:
    debug(
	  if (distances_observed_p == true) {
	    printf("acceptor #%d (%d partners obs, %d partners max):",
		   j,nsplicepartners_obs[j],nsplicepartners_max[j]);
	  } else {
	    printf("acceptor #%d (%d partners max):",
		   j,nsplicepartners_max[j]);
	  });

    j1 = j - 1 - nsplicepartners_skip[j];

    if (distances_observed_p) {
      if (nsplicepartners_obs[j] == 0) {
	triecontents_obs_list = Trie_output_empty(&ptr_obs,&nprinted_obs,triecontents_obs_list);

      } else {
	sites = (Splicestring_T *) CALLOC(nsplicepartners_obs[j],sizeof(Splicestring_T));
	nsites = 0;
	while (nsites < nsplicepartners_obs[j]) {
	  if (splicetypes[j1] == DONOR) {
	    for (p = splicestrings[j1]; p != NULL; p = List_next(p)) {
	      debug(splicestring = (Splicestring_T) List_head(p));
	      debug(splicefrag_nt_leftward(gbuffer,splicestring->string));
	      debug(printf(" %d (%s)",j1,gbuffer));
	      sites[nsites++] = (Splicestring_T) List_head(p);
	    }
	  }
	  j1--;
	}
	assert(nsites == nsplicepartners_obs[j]);
	qsort(sites,nsites,sizeof(Splicestring_T),Splicestring_cmp);
	debug(printf("\n"));
	debug(dump_sites(sites,nsites,splicetypes[j]));

	triecontents_obs_list = Trie_output_new(&ptr_obs,&nprinted_obs,triecontents_obs_list,
						sites,nsites,/*charpos*/0);
	FREE(sites);
      }
    }

    if (nsplicepartners_max[j] == 0) {
      triecontents_max_list = Trie_output_empty(&ptr_max,&nprinted_max,triecontents_max_list);

    } else {
      sites = (Splicestring_T *) CALLOC(nsplicepartners_max[j],sizeof(Splicestring_T));
      nsites = 0;
      while (nsites < nsplicepartners_max[j]) {
	if (splicetypes[j1] == DONOR) {
	  for (p = splicestrings[j1]; p != NULL; p = List_next(p)) {
	    debug(splicestring = (Splicestring_T) List_head(p));
	    debug(splicefrag_nt_leftward(gbuffer,splicestring->string));
	    debug(printf(" %d (%s)",j1,gbuffer));
	    sites[nsites++] = (Splicestring_T) List_head(p);
	  }
	}
	j1--;
      }
      assert(nsites == nsplicepartners_max[j]);
      qsort(sites,nsites,sizeof(Splicestring_T),Splicestring_cmp);
      debug(printf("\n"));
      debug(dump_sites(sites,nsites,splicetypes[j]));

      triecontents_max_list = Trie_output_new(&ptr_max,&nprinted_max,triecontents_max_list,
					      sites,nsites,/*charpos*/0);
      FREE(sites);
    }

    debug(printf("\n"));
    break;

  case ANTIDONOR:
    debug(
	  if (distances_observed_p == true) {
	    printf("antidonor #%d (%d partners obs, %d partners max):",
		   j,nsplicepartners_obs[j],nsplicepartners_max[j]);
	  } else {
	    printf("antidonor #%d (%d partners max):",
		   j,nsplicepartners_max[j]);
	  });

    j1 = j - 1 - nsplicepartners_skip[j];

    if (distances_observed_p) {
      if (nsplicepartners_obs[j] == 0) {
	triecontents_obs_list = Trie_output_empty(&ptr_obs,&nprinted_obs,triecontents_obs_list);

      } else {
	sites = (Splicestring_T *) CALLOC(nsplicepartners_obs[j],sizeof(Splicestring_T));
	nsites = 0;
	while (nsites < nsplicepartners_obs[j]) {
	  if (splicetypes[j1] == ANTIACCEPTOR) {
	    for (p = splicestrings[j1]; p != NULL; p = List_next(p)) {
	      debug(splicestring = (Splicestring_T) List_head(p));
	      debug(splicefrag_nt_leftward(gbuffer,splicestring->string));
	      debug(printf(" %d (%s)",j1,gbuffer));
	      sites[nsites++] = (Splicestring_T) List_head(p);
	    }
	  }
	  j1--;
	}
	assert(nsites = nsplicepartners_obs[j]);
	qsort(sites,nsites,sizeof(Splicestring_T),Splicestring_cmp);
	debug(printf("\n"));
	debug(dump_sites(sites,nsites,splicetypes[j]));

	triecontents_obs_list = Trie_output_new(&ptr_obs,&nprinted_obs,triecontents_obs_list,
						sites,nsites,/*charpos*/0);
	FREE(sites);
      }
    }

    if (nsplicepartners_max[j] == 0) {
      triecontents_max_list = Trie_output_empty(&ptr_max,&nprinted_max,triecontents_max_list);

    } else {
      sites = (Splicestring_T *) CALLOC(nsplicepartners_max[j],sizeof(Splicestring_T));
      nsites = 0;
      while (nsites < nsplicepartners_max[j]) {
	if (splicetypes[j1] == ANTIACCEPTOR) {
	  for (p = splicestrings[j1]; p != NULL; p = List_next(p)) {
	    debug(splicestring = (Splicestring_T) List_head(p));
	    debug(splicefrag_nt_leftward(gbuffer,splicestring->string));
	    debug(printf(" %d (%s)",j1,gbuffer));
	    sites[nsites++] = (Splicestring_T) List_head(p);
	  }
	}
	j1--;
      }
      assert(nsites = nsplicepartners_max[j]);
      qsort(sites,nsites,sizeof(Splicestring_T),Splicestring_cmp);
      debug(printf("\n"));
      debug(dump_sites(sites,nsites,splicetypes[j]));

      triecontents_max_list = Trie_output_new(&ptr_max,&nprinted_max,triecontents_max_list,
					      sites,nsites,/*charpos*/0);
      FREE(sites);
    }

    debug(printf("\n"));
    break;

  case ANTIACCEPTOR:
    debug(
	  if (distances_observed_p == true) {
	    printf("antiacceptor #%d (%d partners obs, %d partners max):",
		   j,nsplicepartners_obs[j],nsplicepartners_max[j]);
	  } else {
	    printf("antiacceptor #%d (%d partners max):",
		   j,nsplicepartners_max[j]);
	  });

    j1 = j + 1 + nsplicepartners_skip[j];

    if (distances_observed_p) {
      if (nsplicepartners_obs[j] == 0) {
	triecontents_obs_list = Trie_output_empty(&ptr_obs,&nprinted_obs,triecontents_obs_list);

      } else {
	sites = (Splicestring_T *) CALLOC(nsplicepartners_obs[j],sizeof(Splicestring_T));
	nsites = 0;
	while (nsites < nsplicepartners_obs[j]) {
	  if (splicetypes[j1] == ANTIDONOR) {
	    for (p = splicestrings[j1]; p != NULL; p = List_next(p)) {
	      debug(splicestring = (Splicestring_T) List_head(p));
	      debug(splicefrag_nt_rightward(gbuffer,splicestring->string));
	      debug(printf(" %d (%s)",j1,gbuffer));
	      sites[nsites++] = (Splicestring_T) List_head(p);
	    }
	  }
	  j1++;
	}
	assert(nsites = nsplicepartners_obs[j]);

	qsort(sites,nsites,sizeof(Splicestring_T),Splicestring_cmp);
	debug(printf("\n"));
	debug(dump_sites(sites,nsites,splicetypes[j]));

	triecontents_obs_list = Trie_output_new(&ptr_obs,&nprinted_obs,triecontents_obs_list,
						sites,nsites,/*charpos*/0);
	FREE(sites);
      }
    }

    if (nsplicepartners_max[j] == 0) {
      triecontents_max_list = Trie_output_empty(&ptr_max,&nprinted_max,triecontents_max_list);

    } else {
      sites = (Splicestring_T *) CALLOC(nsplicepartners_max[j],sizeof(Splicestring_T));
      nsites = 0;
      while (nsites < nsplicepartners_max[j]) {
	if (splicetypes[j1] == ANTIDONOR) {
	  for (p = splicestrings[j1]; p != NULL; p = List_next(p)) {
	    debug(splicestring = (Splicestring_T) List_head(p));
	    debug(splicefrag_nt_rightward(gbuffer,splicestring->string));
	    debug(printf(" %d (%s)",j1,gbuffer));
	    sites[nsites++] = (Splicestring_T) List_head(p);
	  }
	}
	j1++;
      }
      assert(nsites = nsplicepartners_max[j]);

      qsort(sites,nsites,sizeof(Splicestring_T),Splicestring_cmp);
      debug(printf("\n"));
      debug(dump_sites(sites,nsites,splicetypes[j]));

      triecontents_max_list = Trie_output_new(&ptr_max,&nprinted_max,triecontents_max_list,
					      sites,nsites,/*charpos*/0);
      FREE(sites);
    }

    debug(printf("\n"));
    break;

  default:
    fprintf(stderr,"Unexpected splicetype %d at splicesites_i %d\n",splicetypes[j],j);
    abort();

  }

  if (distances_observed_p) {
    triecontents_obs_list = Uintlist_reverse(triecontents_obs_list);
    *triecontents_obs = Uintlist_to_array(&nprinted_obs,triecontents_obs_list);
    Uintlist_free(&triecontents_obs_list);
    *triestart_obs = &((*triecontents_obs)[ptr_obs]);
  } else {
    *triecontents_obs = (unsigned int *) NULL;
    *triestart_obs = (unsigned int *) NULL;
  }

  triecontents_max_list = Uintlist_reverse(triecontents_max_list);
  *triecontents_max = Uintlist_to_array(&nprinted_max,triecontents_max_list);
  Uintlist_free(&triecontents_max_list);
  *triestart_max = &((*triecontents_max)[ptr_max]);

  return;
}



#ifdef OLD_CODE

/* Used for short-overlap splicing */
static int
count_subtree_left (int *splicesite_i, unsigned int *triecontents,
		    UINT4 *splicefrags_ref, UINT4 *splicefrags_alt,
		    int max_mismatches, int maxhits, Compress_T query_compress,
		    bool dibasep, bool cmetp, int pos5, int pos3) {
  unsigned int leaf;
  int nleaves, i;
  int offseta, offsetc, offsetg, offsett;
  UINT4 query_shifted, flags, mask;
  int nmismatches, n;

  assert(pos3 - pos5 <= 16);

  if (single_leaf_p(leaf = triecontents[0])) {
    *splicesite_i = (int) leaf;

    /* Recomputes entire segment to determine mismatches (necessary because of splicefrags) */
    query_shifted = Genome_query_shift_fragment_left(&flags,&mask,query_compress,pos5,pos3);
    nmismatches = Genome_count_mismatches_fragment(query_shifted,flags,mask,
						   splicefrags_ref[leaf],splicefrags_alt[leaf]);

    debug1(printf("Found leaf %u => %d mismatches\n",leaf,nmismatches));
    if (nmismatches <= max_mismatches) {
      return 1;
    } else {
      return 0;
    }

  } else if (multiple_leaf_p(leaf)) {
    nleaves = (int) (-leaf);

    n = 0;
    for (i = 1; i <= nleaves; i++) {
      leaf = triecontents[i];
      /* Recomputes entire segment to determine mismatches (necessary because of splicefrags) */
      query_shifted = Genome_query_shift_fragment_left(&flags,&mask,query_compress,pos5,pos3);
      nmismatches = Genome_count_mismatches_fragment(query_shifted,flags,mask,
						     splicefrags_ref[leaf],splicefrags_alt[leaf]);
      
      debug1(printf("Found leaf %u => %d mismatches\n",leaf,nmismatches));
      if (nmismatches <= max_mismatches) {
	n++;
      }
    }

    return n;

  } else {
#ifdef USE_2BYTE_RELOFFSETS
    get_offsets(&offseta,&offsetc,&offsetg,&offsett,triecontents[1],triecontents[2]);
#else
    offseta = (int) triecontents[1];
    offsetc = (int) triecontents[2];
    offsetg = (int) triecontents[3];
    offsett = (int) triecontents[4];
#endif

    n = 0;
    if (offseta > 0) {
      n += count_subtree_left(&(*splicesite_i),&(triecontents[-offseta]),
			      splicefrags_ref,splicefrags_alt,max_mismatches,maxhits,
			      query_compress,dibasep,cmetp,pos5,pos3);
#ifdef LIMIT_MAXHITS
      if (n > maxhits) {
	return n;
      }
#endif
    }
    if (offsetc > 0) {
      n += count_subtree_left(&(*splicesite_i),&(triecontents[-offsetc]),
			      splicefrags_ref,splicefrags_alt,max_mismatches,maxhits,
			      query_compress,dibasep,cmetp,pos5,pos3);
#ifdef LIMIT_MAXHITS
      if (n > maxhits) {
	return n;
      }
#endif
    }
    if (offsetg > 0) {
      n += count_subtree_left(&(*splicesite_i),&(triecontents[-offsetg]),
			      splicefrags_ref,splicefrags_alt,max_mismatches,maxhits,
			      query_compress,dibasep,cmetp,pos5,pos3);
#ifdef LIMIT_MAXHITS
      if (n > maxhits) {
	return n;
      }
#endif
    }
    if (offsett > 0) {
      n += count_subtree_left(&(*splicesite_i),&(triecontents[-offsett]),
			      splicefrags_ref,splicefrags_alt,max_mismatches,maxhits,
			      query_compress,dibasep,cmetp,pos5,pos3);
#ifdef LIMIT_MAXHITS
      if (n > maxhits) {
	return n;
      }
#endif
    }

    return n;
  }
}

#endif

#ifdef OLD_CODE

static int
count_subtree_right (int *splicesite_i, unsigned int *triecontents,
		     UINT4 *splicefrags_ref, UINT4 *splicefrags_alt,
		     int max_mismatches, int maxhits, Compress_T query_compress,
		     bool dibasep, bool cmetp, int pos5, int pos3) {
  unsigned int leaf;
  int nleaves, i;
  int offseta, offsetc, offsetg, offsett;
  UINT4 query_shifted, flags, mask;
  int nmismatches, n;

  assert(pos3 - pos5 <= 16);

  if (single_leaf_p(leaf = triecontents[0])) {
    *splicesite_i = (int) leaf;

    /* Recomputes entire segment to determine mismatches (necessary because of splicefrags) */
    query_shifted = Genome_query_shift_fragment_right(&flags,&mask,query_compress,pos5,pos3);
    nmismatches = Genome_count_mismatches_fragment(query_shifted,flags,mask,
						   splicefrags_ref[leaf],splicefrags_alt[leaf]);

    debug1(printf("Found leaf %u => %d mismatches\n",leaf,nmismatches));
    if (nmismatches <= max_mismatches) {
      return 1;
    } else {
      return 0;
    }

  } else if (multiple_leaf_p(leaf)) {
    nleaves = (int) (-leaf);

    n = 0;
    for (i = 1; i <= nleaves; i++) {
      leaf = triecontents[i];

      /* Recomputes entire segment to determine mismatches (necessary because of splicefrags) */
      query_shifted = Genome_query_shift_fragment_right(&flags,&mask,query_compress,pos5,pos3);
      nmismatches = Genome_count_mismatches_fragment(query_shifted,flags,mask,
						     splicefrags_ref[leaf],splicefrags_alt[leaf]);
      
      debug1(printf("Found leaf %u => %d mismatches\n",leaf,nmismatches));
      if (nmismatches <= max_mismatches) {
	n++;
      }
    }

    return n;

  } else {
#ifdef USE_2BYTE_RELOFFSETS
    get_offsets(&offseta,&offsetc,&offsetg,&offsett,triecontents[1],triecontents[2]);
#else
    offseta = (int) triecontents[1];
    offsetc = (int) triecontents[2];
    offsetg = (int) triecontents[3];
    offsett = (int) triecontents[4];
#endif

    n = 0;
    if (offseta > 0) {
      n += count_subtree_right(&(*splicesite_i),&(triecontents[-offseta]),
			       splicefrags_ref,splicefrags_alt,max_mismatches,maxhits,
			       query_compress,dibasep,cmetp,pos5,pos3);
#ifdef LIMIT_MAXHITS
      if (n > maxhits) {
	return n;
      }
#endif
    }
    if (offsetc > 0) {
      n += count_subtree_right(&(*splicesite_i),&(triecontents[-offsetc]),
			       splicefrags_ref,splicefrags_alt,max_mismatches,maxhits,
			       query_compress,dibasep,cmetp,pos5,pos3);
#ifdef LIMIT_MAXHITS
      if (n > maxhits) {
	return n;
      }
#endif
    }
    if (offsetg > 0) {
      n += count_subtree_right(&(*splicesite_i),&(triecontents[-offsetg]),
			       splicefrags_ref,splicefrags_alt,max_mismatches,maxhits,
			       query_compress,dibasep,cmetp,pos5,pos3);
#ifdef LIMIT_MAXHITS
      if (n > maxhits) {
	return n;
      }
#endif
    }
    if (offsett > 0) {
      n += count_subtree_right(&(*splicesite_i),&(triecontents[-offsett]),
			       splicefrags_ref,splicefrags_alt,max_mismatches,maxhits,
			       query_compress,dibasep,cmetp,pos5,pos3);
#ifdef LIMIT_MAXHITS
      if (n > maxhits) {
	return n;
      }
#endif
    }

    return n;
  }
}

#endif

#ifdef OLD_CODE

static void
search_subtree_left (int *splicesite_i, int *best_nmismatches, int *nbest,
		     unsigned int *triecontents, Genomicpos_T *splicesites,
		     UINT4 *splicefrags_ref, UINT4 *splicefrags_alt, char *query,
		     Compress_T query_compress, UINT4 *genome_blocks, UINT4 *snp_blocks,
		     bool dibasep, bool cmetp, int pos5, int pos3, bool plusp) {
  unsigned int leaf;
  int nleaves, i;
  int offseta, offsetc, offsetg, offsett;
  UINT4 query_shifted, flags, mask;
  int nmismatches, ncolordiffs;
  Genomicpos_T segment_left;

  if (single_leaf_p(leaf = triecontents[0])) {
    if (pos3 - pos5 <= 16) {
      /* Recomputes entire segment to determine mismatches (necessary because of splicefrags) */
      query_shifted = Genome_query_shift_fragment_left(&flags,&mask,query_compress,pos5,pos3);
      nmismatches = Genome_count_mismatches_fragment(query_shifted,flags,mask,
						     splicefrags_ref[leaf],splicefrags_alt[leaf]);
    } else {
      /* Can happen in search for short middle exon */
      segment_left = splicesites[leaf] - pos3;
      nmismatches =
	Genome_count_mismatches_substring(&ncolordiffs,query,query_compress,genome_blocks,snp_blocks,
					  segment_left,pos5,pos3,dibasep,cmetp,plusp);
    }

    debug2(printf("Found leaf %u => %d mismatches\n",leaf,nmismatches));
    if (nmismatches < *best_nmismatches) {
      *best_nmismatches = nmismatches;
      *splicesite_i = (int) leaf;
      *nbest = 1;
    } else if (nmismatches == *best_nmismatches) {
      *nbest += 1;
    }

  } else if (multiple_leaf_p(leaf)) {
    nleaves = (int) (-leaf);

    for (i = 1; i <= nleaves; i++) {
      leaf = triecontents[i];

      if (pos3 - pos5 <= 16) {
	/* Recomputes entire segment to determine mismatches (necessary because of splicefrags) */
	query_shifted = Genome_query_shift_fragment_left(&flags,&mask,query_compress,pos5,pos3);
	nmismatches = Genome_count_mismatches_fragment(query_shifted,flags,mask,
						       splicefrags_ref[leaf],splicefrags_alt[leaf]);
      } else {
	/* Can happen in search for short middle exon */
	segment_left = splicesites[leaf] - pos3;
	nmismatches =
	  Genome_count_mismatches_substring(&ncolordiffs,query,query_compress,genome_blocks,snp_blocks,
					    segment_left,pos5,pos3,dibasep,cmetp,plusp);
      }

      debug2(printf("Found leaf %u => %d mismatches\n",leaf,nmismatches));
      if (nmismatches < *best_nmismatches) {
	*best_nmismatches = nmismatches;
	*splicesite_i = (int) leaf;
	*nbest = 1;
      } else if (nmismatches == *best_nmismatches) {
	*nbest += 1;
      }
    }

  } else {
#ifdef USE_2BYTE_RELOFFSETS
    get_offsets(&offseta,&offsetc,&offsetg,&offsett,triecontents[1],triecontents[2]);
#else
    offseta = (int) triecontents[1];
    offsetc = (int) triecontents[2];
    offsetg = (int) triecontents[3];
    offsett = (int) triecontents[4];
#endif

    if (offseta > 0) {
      search_subtree_left(&(*splicesite_i),&(*best_nmismatches),&(*nbest),
			  &(triecontents[-offseta]),splicesites,
			  splicefrags_ref,splicefrags_alt,
			  query,query_compress,genome_blocks,snp_blocks,
			  dibasep,cmetp,pos5,pos3,plusp);
    }
    if (offsetc > 0) {
      search_subtree_left(&(*splicesite_i),&(*best_nmismatches),&(*nbest),
			  &(triecontents[-offsetc]),splicesites,
			  splicefrags_ref,splicefrags_alt,
			  query,query_compress,genome_blocks,snp_blocks,
			  dibasep,cmetp,pos5,pos3,plusp);
    }
    if (offsetg > 0) {
      search_subtree_left(&(*splicesite_i),&(*best_nmismatches),&(*nbest),
			  &(triecontents[-offsetg]),splicesites,
			  splicefrags_ref,splicefrags_alt,
			  query,query_compress,genome_blocks,snp_blocks,
			  dibasep,cmetp,pos5,pos3,plusp);
    }
    if (offsett > 0) {
      search_subtree_left(&(*splicesite_i),&(*best_nmismatches),&(*nbest),
			  &(triecontents[-offsett]),splicesites,
			  splicefrags_ref,splicefrags_alt,
			  query,query_compress,genome_blocks,snp_blocks,
			  dibasep,cmetp,pos5,pos3,plusp);
    }
  }

  return;
}

#endif

#ifdef OLD_CODE

static void
search_subtree_right (int *splicesite_i, int *best_nmismatches, int *nbest,
		      unsigned int *triecontents, Genomicpos_T *splicesites,
		      UINT4 *splicefrags_ref, UINT4 *splicefrags_alt, char *query,
		      Compress_T query_compress, UINT4 *genome_blocks, UINT4 *snp_blocks,
		      bool dibasep, bool cmetp, int pos5, int pos3, bool plusp) {
  unsigned int leaf;
  int nleaves, i;
  int offseta, offsetc, offsetg, offsett;
  UINT4 query_shifted, flags, mask;
  int nmismatches, ncolordiffs;
  Genomicpos_T segment_left;

  if (single_leaf_p(leaf = triecontents[0])) {
    if (pos3 - pos5 <= 16) {
      /* Recomputes entire segment to determine mismatches (necessary because of splicefrags) */
      query_shifted = Genome_query_shift_fragment_right(&flags,&mask,query_compress,pos5,pos3);
      nmismatches = Genome_count_mismatches_fragment(query_shifted,flags,mask,
						     splicefrags_ref[leaf],splicefrags_alt[leaf]);
    } else {
      /* Can happen in search for short middle exon */
      segment_left = splicesites[leaf] - pos5;
      nmismatches =
	Genome_count_mismatches_substring(&ncolordiffs,query,query_compress,genome_blocks,snp_blocks,
					  segment_left,pos5,pos3,dibasep,cmetp,plusp);
    }

    debug2(printf("Found leaf %u => %d mismatches\n",leaf,nmismatches));
    if (nmismatches < *best_nmismatches) {
      *best_nmismatches = nmismatches;
      *splicesite_i = (int) leaf;
      *nbest = 1;
    } else if (nmismatches == *best_nmismatches) {
      *nbest += 1;
    }

  } else if (multiple_leaf_p(leaf)) {
    nleaves = (int) (-leaf);

    for (i = 1; i <= nleaves; i++) {
      leaf = triecontents[i];

      if (pos3 - pos5 <= 16) {
	/* Recomputes entire segment to determine mismatches (necessary because of splicefrags) */
	query_shifted = Genome_query_shift_fragment_right(&flags,&mask,query_compress,pos5,pos3);
	nmismatches = Genome_count_mismatches_fragment(query_shifted,flags,mask,
						       splicefrags_ref[leaf],splicefrags_alt[leaf]);
      } else {
	/* Can happen in search for short middle exon */
	segment_left = splicesites[leaf] - pos5;
	nmismatches =
	  Genome_count_mismatches_substring(&ncolordiffs,query,query_compress,genome_blocks,snp_blocks,
					    segment_left,pos5,pos3,dibasep,cmetp,plusp);
      }
      
      debug2(printf("Found leaf %u => %d mismatches\n",leaf,nmismatches));
      if (nmismatches < *best_nmismatches) {
	*best_nmismatches = nmismatches;
	*splicesite_i = (int) leaf;
	*nbest = 1;
      } else if (nmismatches == *best_nmismatches) {
	*nbest += 1;
      }
    }

  } else {
#ifdef USE_2BYTE_RELOFFSETS
    get_offsets(&offseta,&offsetc,&offsetg,&offsett,triecontents[1],triecontents[2]);
#else
    offseta = (int) triecontents[1];
    offsetc = (int) triecontents[2];
    offsetg = (int) triecontents[3];
    offsett = (int) triecontents[4];
#endif

    if (offseta > 0) {
      search_subtree_right(&(*splicesite_i),&(*best_nmismatches),&(*nbest),
			   &(triecontents[-offseta]),splicesites,
			   splicefrags_ref,splicefrags_alt,
			   query,query_compress,genome_blocks,snp_blocks,
			   dibasep,cmetp,pos5,pos3,plusp);
    }
    if (offsetc > 0) {
      search_subtree_right(&(*splicesite_i),&(*best_nmismatches),&(*nbest),
			   &(triecontents[-offsetc]),splicesites,
			   splicefrags_ref,splicefrags_alt,
			   query,query_compress,genome_blocks,snp_blocks,
			   dibasep,cmetp,pos5,pos3,plusp);
    }
    if (offsetg > 0) {
      search_subtree_right(&(*splicesite_i),&(*best_nmismatches),&(*nbest),
			   &(triecontents[-offsetg]),splicesites,
			   splicefrags_ref,splicefrags_alt,
			   query,query_compress,genome_blocks,snp_blocks,
			   dibasep,cmetp,pos5,pos3,plusp);
    }
    if (offsett > 0) {
      search_subtree_right(&(*splicesite_i),&(*best_nmismatches),&(*nbest),
			   &(triecontents[-offsett]),splicesites,
			   splicefrags_ref,splicefrags_alt,
			   query,query_compress,genome_blocks,snp_blocks,
			   dibasep,cmetp,pos5,pos3,plusp);
    }
  }

  return;
}

#endif


/**********************************************************************************/

#ifdef OLD_CODE

/* exclude is designed for apparent similarity due to SNP */
/* Call initially with charpos = pos3 */
static int
count_left (int *splicesite_i, int exclude_i, unsigned int *triecontents, 
	    UINT4 *splicefrags_ref, UINT4 *splicefrags_alt, char *queryptr,
	    int max_mismatches, int maxhits, Compress_T query_compress,
	    bool dibasep, bool cmetp, int pos5, int pos3, int nmismatches, int charpos) {
  Genomicpos_T leaf;
  int nleaves, i;
  UINT4 query_shifted, flags, mask;
  int offseta, offsetc, offsetg, offsett, jump;
  int ncolordiffs;
  char c;
  int n;
  
  assert(pos3 - pos5 <= 16);

  debug1(printf("Entered Trie_count_left with nmismatches %d, charpos %d\n",nmismatches,charpos));

  if (single_leaf_p(leaf = triecontents[0])) {

    if (charpos - 1 >= pos5) {
      if ((int) leaf == exclude_i) {
	debug1(printf("Found leaf %u, but excluded\n",leaf));
	return 0;
      } else {
	/* Recomputes entire segment to determine mismatches (necessary because of splicefrags) */
	query_shifted = Genome_query_shift_fragment_left(&flags,&mask,query_compress,pos5,pos3);
	nmismatches = Genome_count_mismatches_fragment(query_shifted,flags,mask,
						       splicefrags_ref[leaf],splicefrags_alt[leaf]);

	debug1(printf("Found leaf %u, but still have characters to check against Genome: %.*s => %d mismatches\n",
		      leaf,charpos + 1,&(queryptr[0]),nmismatches));
	if (nmismatches <= max_mismatches) {
	  *splicesite_i = (int) leaf;
	  return 1;
	} else {
	  return 0;
	}
      }
      
    } else {
      debug1(printf("Found leaf %u, and completely finished query, so unique\n",leaf));
      *splicesite_i = (int) leaf;
      return 1;
    }

  } else if (multiple_leaf_p(leaf)) {
    nleaves = (int) (-leaf);

    n = 0;
    for (i = 1; i <= nleaves; i++) {
      leaf = triecontents[i];

      if (charpos - 1 >= pos5) {
	if ((int) leaf == exclude_i) {
	  debug1(printf("Found leaf %u, but excluded\n",leaf));
	} else {
	  /* Recomputes entire segment to determine mismatches (necessary because of splicefrags) */
	  query_shifted = Genome_query_shift_fragment_left(&flags,&mask,query_compress,pos5,pos3);
	  nmismatches = Genome_count_mismatches_fragment(query_shifted,flags,mask,
							 splicefrags_ref[leaf],splicefrags_alt[leaf]);

	  debug1(printf("Found leaf %u, but still have characters to check against Genome: %.*s => %d mismatches\n",
			leaf,charpos + 1,&(queryptr[0]),nmismatches));
	  if (nmismatches <= max_mismatches) {
	    *splicesite_i = (int) leaf;
	    n++;
	  }
	}
      
      } else {
	debug1(printf("Found leaf %u, and completely finished query, so unique\n",leaf));
	*splicesite_i = (int) leaf;
	n++;
      }
    }

    return n;

  } else if (charpos - 1 < pos5) {
    if (splicefrags_alt == splicefrags_ref) {
      /* Non-leaf, and no more characters left, so have multiple hits */
      debug1(printf("Non-leaf, but no more characters left\n"));
#ifdef LIMIT_MAXHITS
      return 2;
#else
      return count_subtree_left(&(*splicesite_i),triecontents,
				splicefrags_ref,splicefrags_alt,
				max_mismatches,maxhits,query_compress,
				dibasep,cmetp,pos5,pos3);
#endif

    } else {
      /* Alternate genome.  Need to re-check entire subtree, because nmismatches is just a lower bound */
      debug1(printf("Non-leaf, no more characters left, so check all splicefrags\n"));
      return count_subtree_left(&(*splicesite_i),triecontents,
				splicefrags_ref,splicefrags_alt,
				max_mismatches,maxhits,query_compress,
				dibasep,cmetp,pos5,pos3);
    }

  } else {
    /* Non-leaf, and characters left, so recurse */
    c = queryptr[charpos - 1];

#ifdef USE_2BYTE_RELOFFSETS
    get_offsets(&offseta,&offsetc,&offsetg,&offsett,triecontents[1],triecontents[2]);
#else
    offseta = (int) triecontents[1];
    offsetc = (int) triecontents[2];
    offsetg = (int) triecontents[3];
    offsett = (int) triecontents[4];
#endif

    debug1(printf("Got offsets %d %d %d %d",offseta,offsetc,offsetg,offsett));
    if (nmismatches == max_mismatches) {
      /* Continue only with a match */
      switch (c) {
      case 'A': jump = offseta; break;
      case 'C': jump = offsetc; break;
      case 'G': jump = offsetg; break;
      case 'T': jump = offsett; break;
      default: jump = 0;
      }

      if (jump == 0) {
	debug1(printf(" => %c: No matches\n",c));
	return 0;
      } else {
	debug1(printf(" => %c: Going to offset %d\n",c,jump));
	return count_left(&(*splicesite_i),exclude_i,&(triecontents[-jump]),
			  splicefrags_ref,splicefrags_alt,queryptr,
			  max_mismatches,maxhits,query_compress,
			  dibasep,cmetp,pos5,pos3,nmismatches,charpos-1);
      }

    } else {
      debug1(printf("\n"));
      n = 0;
      if (offseta > 0) {
	debug1(printf("Branching to A: offset %d\n",offseta));
	n += count_left(&(*splicesite_i),exclude_i,&(triecontents[-offseta]),
			splicefrags_ref,splicefrags_alt,queryptr,
			max_mismatches,maxhits,query_compress,
			dibasep,cmetp,pos5,pos3,nmismatches+(c != 'A'),charpos-1);
#ifdef LIMIT_MAXHITS
	if (n > maxhits) {
	  debug1(printf("Returning %d\n",n));
	  return n;
	}
#endif
      }
      
      if (offsetc > 0) {
	debug1(printf("Branching to C: offset %d\n",offsetc));
	n += count_left(&(*splicesite_i),exclude_i,&(triecontents[-offsetc]),
			splicefrags_ref,splicefrags_alt,queryptr,
			max_mismatches,maxhits,query_compress,
			dibasep,cmetp,pos5,pos3,nmismatches+(c != 'C'),charpos-1);
#ifdef LIMIT_MAXHITS
	if (n > maxhits) {
	  debug1(printf("Returning %d\n",n));
	  return n;
	}
#endif
      }

      if (offsetg > 0) {
	debug1(printf("Branching to G: offset %d\n",offsetg));
	n += count_left(&(*splicesite_i),exclude_i,&(triecontents[-offsetg]),
			splicefrags_ref,splicefrags_alt,queryptr,
			max_mismatches,maxhits,query_compress,
			dibasep,cmetp,pos5,pos3,nmismatches+(c != 'G'),charpos-1);
#ifdef LIMIT_MAXHITS
	if (n > maxhits) {
	  debug1(printf("Returning %d\n",n));
	  return n;
	}
#endif
      }

      if (offsett > 0) {
	debug1(printf("Branching to T: offset %d\n",offsett));
	n += count_left(&(*splicesite_i),exclude_i,&(triecontents[-offsett]),
			splicefrags_ref,splicefrags_alt,queryptr,
			max_mismatches,maxhits,query_compress,
			dibasep,cmetp,pos5,pos3,nmismatches+(c != 'T'),charpos-1);
#ifdef LIMIT_MAXHITS
	if (n > maxhits) {
	  debug1(printf("Returning %d\n",n));
	  return n;
	}
#endif
      }

      debug1(printf("Returning %d\n",n));
      return n;
    }
  }
}

#endif

#ifdef OLD_CODE

/* exclude is designed for apparent similarity due to SNP */
/* Call initially with charpos = pos5 - 1 */
static int
count_right (int *splicesite_i, int exclude_i, unsigned int *triecontents,
	     UINT4 *splicefrags_ref, UINT4 *splicefrags_alt, char *queryptr,
	     int max_mismatches, int maxhits, Compress_T query_compress,
	     bool dibasep, bool cmetp, int pos5, int pos3, int nmismatches, int charpos) {

  Genomicpos_T leaf;
  int nleaves, i;
  UINT4 query_shifted, flags, mask;
  int offseta, offsetc, offsetg, offsett, jump;
  int ncolordiffs;
  char c;
  int n;
  
  assert(pos3 - pos5 <= 16);

  debug1(printf("Entered Trie_count_right with nmismatches %d, charpos %d\n",nmismatches,charpos));

  if (single_leaf_p(leaf = triecontents[0])) {

    if (charpos + 1 < pos3) {
      if ((int) leaf == exclude_i) {
	debug1(printf("Found leaf %u, but excluded\n",leaf));
	return 0;
      } else {
	/* Recomputes entire segment to determine mismatches (necessary because of splicefrags) */
	query_shifted = Genome_query_shift_fragment_right(&flags,&mask,query_compress,pos5,pos3);
	nmismatches = Genome_count_mismatches_fragment(query_shifted,flags,mask,
						       splicefrags_ref[leaf],splicefrags_alt[leaf]);

	debug1(printf("Found leaf %u, but still have characters to check against Genome: %s => %d mismatches\n",
		      leaf,&(queryptr[charpos + 1]),nmismatches));
	if (nmismatches <= max_mismatches) {
	  *splicesite_i = (int) leaf;
	  return 1;
	} else {
	  return 0;
	}
      }
      
    } else {
      debug1(printf("Found leaf %u, and completely finished query, so unique\n",leaf));
      *splicesite_i = (int) leaf;
      return 1;
    }

  } else if (multiple_leaf_p(leaf)) {
    nleaves = (int) (-leaf);

    n = 0;
    for (i = 1; i <= nleaves; i++) {
      leaf = triecontents[i];

      if (charpos + 1 < pos3) {
	if ((int) leaf == exclude_i) {
	  debug1(printf("Found leaf %u, but excluded\n",leaf));
	} else {
	  /* Recomputes entire segment to determine mismatches (necessary because of splicefrags) */
	  query_shifted = Genome_query_shift_fragment_right(&flags,&mask,query_compress,pos5,pos3);
	  nmismatches = Genome_count_mismatches_fragment(query_shifted,flags,mask,
							 splicefrags_ref[leaf],splicefrags_alt[leaf]);

	  debug1(printf("Found leaf %u, but still have characters to check against Genome: %s => %d mismatches\n",
			leaf,&(queryptr[charpos + 1]),nmismatches));
	  if (nmismatches <= max_mismatches) {
	    *splicesite_i = (int) leaf;
	    n++;
	  }
	}
      
      } else {
	debug1(printf("Found leaf %u, and completely finished query, so unique\n",leaf));
	*splicesite_i = (int) leaf;
	n++;
      }
    }

    return n;

  } else if (charpos + 1 >= pos3) {
    if (splicefrags_alt == splicefrags_ref) {
      /* Non-leaf, and no more characters left, so have multiple hits */
      debug1(printf("Non-leaf, but no more characters left\n"));
#ifdef LIMIT_MAXHITS
      return 2;
#else
      return count_subtree_right(&(*splicesite_i),triecontents,
				 splicefrags_ref,splicefrags_alt,
				 max_mismatches,maxhits,query_compress,
				 dibasep,cmetp,pos5,pos3);
#endif

    } else {
      /* Alternate genome.  Need to re-check entire subtree, because nmismatches is just a lower bound */
      debug1(printf("Non-leaf, no more characters left, so check all splicefrags\n"));
      return count_subtree_right(&(*splicesite_i),triecontents,
				 splicefrags_ref,splicefrags_alt,
				 max_mismatches,maxhits,query_compress,
				 dibasep,cmetp,pos5,pos3);
    }

  } else {
    /* Non-leaf, and characters left, so recurse */
    c = queryptr[charpos + 1];

#ifdef USE_2BYTE_RELOFFSETS
    get_offsets(&offseta,&offsetc,&offsetg,&offsett,triecontents[1],triecontents[2]);
#else
    offseta = (int) triecontents[1];
    offsetc = (int) triecontents[2];
    offsetg = (int) triecontents[3];
    offsett = (int) triecontents[4];
#endif

    debug1(printf("Got offsets %d %d %d %d\n",offseta,offsetc,offsetg,offsett));
    if (nmismatches == max_mismatches) {
      /* Continue only with a match */
      switch (c) {
      case 'A': jump = offseta; break;
      case 'C': jump = offsetc; break;
      case 'G': jump = offsetg; break;
      case 'T': jump = offsett; break;
      default: jump = 0;
      }

      if (jump == 0) {
	debug1(printf("No matches\n"));
	return 0;
      } else {
	debug1(printf("Going to offset %d\n",jump));
	return count_right(&(*splicesite_i),exclude_i,&(triecontents[-jump]),
			   splicefrags_ref,splicefrags_alt,queryptr,
			   max_mismatches,maxhits,query_compress,
			   dibasep,cmetp,pos5,pos3,nmismatches,charpos+1);
      }

    } else {
      n = 0;
      if (offseta > 0) {
	n += count_right(&(*splicesite_i),exclude_i,&(triecontents[-offseta]),
			 splicefrags_ref,splicefrags_alt,queryptr,
			 max_mismatches,maxhits,query_compress,
			 dibasep,cmetp,pos5,pos3,nmismatches+(c != 'A'),charpos+1);
#ifdef LIMIT_MAXHITS
	if (n > maxhits) {
	  debug1(printf("Returning %d\n",n));
	  return n;
	}
#endif
      }
      
      if (offsetc > 0) {
	n += count_right(&(*splicesite_i),exclude_i,&(triecontents[-offsetc]),
			 splicefrags_ref,splicefrags_alt,queryptr,
			 max_mismatches,maxhits,query_compress,
			 dibasep,cmetp,pos5,pos3,nmismatches+(c != 'C'),charpos+1);
#ifdef LIMIT_MAXHITS
	if (n > maxhits) {
	  debug1(printf("Returning %d\n",n));
	  return n;
	}
#endif
      }

      if (offsetg > 0) {
	n += count_right(&(*splicesite_i),exclude_i,&(triecontents[-offsetg]),
			 splicefrags_ref,splicefrags_alt,queryptr,
			 max_mismatches,maxhits,query_compress,
			 dibasep,cmetp,pos5,pos3,nmismatches+(c != 'G'),charpos+1);
#ifdef LIMIT_MAXHITS
	if (n > maxhits) {
	  debug1(printf("Returning %d\n",n));
	  return n;
	}
#endif
      }

      if (offsett > 0) {
	n += count_right(&(*splicesite_i),exclude_i,&(triecontents[-offsett]),
			 splicefrags_ref,splicefrags_alt,queryptr,
			 max_mismatches,maxhits,query_compress,
			 dibasep,cmetp,pos5,pos3,nmismatches+(c != 'T'),charpos+1);
#ifdef LIMIT_MAXHITS
	if (n > maxhits) {
	  debug1(printf("Returning %d\n",n));
	  return n;
	}
#endif
      }

      debug1(printf("Returning %d\n",n));
      return n;
    }
  }
}

#endif

#ifdef OLD_CODE

static void
search_left_single (int *splicesite_i, int *best_nmismatches, int *nbest,
		    unsigned int *triecontents, Genomicpos_T *splicesites,
		    UINT4 *splicefrags_ref, UINT4 *splicefrags_alt, char *query, char *queryptr,
		    Compress_T query_compress, UINT4 *genome_blocks, UINT4 *snp_blocks,
		    bool dibasep, bool cmetp, int pos5, int pos3, bool plusp,
		    int max_mismatches, int nmismatches, int charpos) {
  Genomicpos_T leaf, segment_left;
  int nleaves, i;
  UINT4 query_shifted, flags, mask;
  int offseta, offsetc, offsetg, offsett, jump;
  int ncolordiffs;
  char c;
  
  debug2(printf("Entered search_left_single with nmismatches %d, charpos %d, best_nmismatches %d, nbest %d, bestj %d\n",
		nmismatches,charpos,*best_nmismatches,*nbest,*splicesite_i));

  if (single_leaf_p(leaf = triecontents[0])) {

    if (charpos - 1 >= pos5) {
      if (pos3 - pos5 <= 16) {
	/* Recomputes entire segment to determine mismatches (necessary because of splicefrags) */
	query_shifted = Genome_query_shift_fragment_left(&flags,&mask,query_compress,pos5,pos3);
	nmismatches = Genome_count_mismatches_fragment(query_shifted,flags,mask,
						       splicefrags_ref[leaf],splicefrags_alt[leaf]);
      } else {
	/* Can happen in search for short middle exon */
	segment_left = splicesites[leaf] - pos3;
	nmismatches =
	  Genome_count_mismatches_substring(&ncolordiffs,query,query_compress,genome_blocks,snp_blocks,
					    segment_left,pos5,pos3,dibasep,cmetp,plusp);
      }

      debug2(printf("Found leaf %u, but still have characters to check against Genome: %.*s => %d mismatches\n",
		    leaf,charpos + 1,&(queryptr[0]),nmismatches));
      if (nmismatches < *best_nmismatches) {
	debug2(printf("  nmismatches %d < best_nmismatches %d => setting bestj to be %d (new best_nmismatches %d)\n",
		      nmismatches,*best_nmismatches,(int) leaf,nmismatches));
	*best_nmismatches = nmismatches;
	*splicesite_i = (int) leaf;
	*nbest = 1;
      } else if (nmismatches == *best_nmismatches) {
	*nbest += 1;
	debug2(printf("  nmismatches %d == best_nmismatches %d => incrementing nbest to be %d (same best_nmismatches %d)\n",
		      nmismatches,*best_nmismatches,*nbest,nmismatches));
      }
      
    } else {
      debug2(printf("Found leaf %u, and completely finished query, so unique\n",leaf));
      if (nmismatches < *best_nmismatches) {
	debug2(printf("  nmismatches %d < best_nmismatches %d => setting bestj to be %d (new best_nmismatches %d)\n",
		      nmismatches,*best_nmismatches,(int) leaf,nmismatches));
	*best_nmismatches = nmismatches;
	*splicesite_i = (int) leaf;
	*nbest = 1;
      } else if (nmismatches == *best_nmismatches) {
	*nbest += 1;
	debug2(printf("  nmismatches %d == best_nmismatches %d => incrementing nbest to be %d (same best_nmismatches %d)\n",
		      nmismatches,*best_nmismatches,*nbest,nmismatches));
      }
    }

    return;

  } else if (multiple_leaf_p(leaf)) {
    nleaves = (int) (-leaf);

    for (i = 1; i <= nleaves; i++) {
      leaf = triecontents[i];

      if (charpos - 1 >= pos5) {
	if (pos3 - pos5 <= 16) {
	  /* Recomputes entire segment to determine mismatches (necessary because of splicefrags) */
	  query_shifted = Genome_query_shift_fragment_left(&flags,&mask,query_compress,pos5,pos3);
	  nmismatches = Genome_count_mismatches_fragment(query_shifted,flags,mask,
							 splicefrags_ref[leaf],splicefrags_alt[leaf]);
	} else {
	  /* Can happen in search for short middle exon */
	  segment_left = splicesites[leaf] - pos3;
	  nmismatches =
	    Genome_count_mismatches_substring(&ncolordiffs,query,query_compress,genome_blocks,snp_blocks,
					      segment_left,pos5,pos3,dibasep,cmetp,plusp);
	}

	debug2(printf("Found leaf %u, but still have characters to check against Genome: %.*s => %d mismatches\n",
		      leaf,charpos + 1,&(queryptr[0]),nmismatches));
	if (nmismatches < *best_nmismatches) {
	  debug2(printf("  nmismatches %d < best_nmismatches %d => setting bestj to be %d (new best_nmismatches %d)\n",
			nmismatches,*best_nmismatches,(int) leaf,nmismatches));
	  *best_nmismatches = nmismatches;
	  *splicesite_i = (int) leaf;
	  *nbest = 1;
	} else if (nmismatches == *best_nmismatches) {
	  *nbest += 1;
	  debug2(printf("  nmismatches %d == best_nmismatches %d => incrementing nbest to be %d (same best_nmismatches %d)\n",
			nmismatches,*best_nmismatches,*nbest,nmismatches));
	}
      
      } else {
	debug2(printf("Found leaf %u, and completely finished query, so unique\n",leaf));
	if (nmismatches < *best_nmismatches) {
	  debug2(printf("  nmismatches %d < best_nmismatches %d => setting bestj to be %d (new best_nmismatches %d)\n",
			nmismatches,*best_nmismatches,(int) leaf,nmismatches));
	  *best_nmismatches = nmismatches;
	  *splicesite_i = (int) leaf;
	  *nbest = 1;
	} else if (nmismatches == *best_nmismatches) {
	  *nbest += 1;
	  debug2(printf("  nmismatches %d == best_nmismatches %d => incrementing nbest to be %d (same best_nmismatches %d)\n",
			nmismatches,*best_nmismatches,*nbest,nmismatches));
	}
      }
    }

    return;

  } else if (charpos - 1 < pos5) {
    if (splicefrags_alt == splicefrags_ref) {
      /* Non-leaf, and no more characters left, so have multiple hits */
      debug2(printf("Non-leaf, but no more characters left: nmismatches %d vs best_nmismatches %d\n",
		    nmismatches,*best_nmismatches));
      if (nmismatches < *best_nmismatches) {
	debug2(printf("  nmismatches %d < best_nmismatches %d => setting nbest to be %d (new best_nmismatches %d)\n",
		      nmismatches,*best_nmismatches,2,nmismatches));
	*best_nmismatches = nmismatches;
	/* *splicesite_i = (int) leaf; -- What should we set splicesite_i to be? */
	*nbest = 2;
      } else if (nmismatches == *best_nmismatches) {
	*nbest += 2;
	debug2(printf("  nmismatches %d == best_nmismatches %d => incrementing nbest to be %d (same best_nmismatches %d)\n",
		      nmismatches,*best_nmismatches,*nbest,nmismatches));
      }
      return;

    } else {
      /* Alternate genome.  Need to re-check entire subtree, because nmismatches is just a lower bound */
      debug2(printf("Non-leaf, no more characters left, so check all splicefrags\n"));
      search_subtree_left(&(*splicesite_i),&(*best_nmismatches),&(*nbest),
			  triecontents,splicesites,splicefrags_ref,splicefrags_alt,
			  query,query_compress,genome_blocks,snp_blocks,
			  dibasep,cmetp,pos5,pos3,plusp);
      return;
    }

  } else {
    /* Non-leaf, and characters left, so recurse */
    c = queryptr[charpos - 1];

#ifdef USE_2BYTE_RELOFFSETS
    get_offsets(&offseta,&offsetc,&offsetg,&offsett,triecontents[1],triecontents[2]);
#else
    offseta = (int) triecontents[1];
    offsetc = (int) triecontents[2];
    offsetg = (int) triecontents[3];
    offsett = (int) triecontents[4];
#endif

    debug2(printf("Got offsets %d %d %d %d",offseta,offsetc,offsetg,offsett));
    if (nmismatches == max_mismatches) {
      /* Continue only with a match */
      switch (c) {
      case 'A': jump = offseta; break;
      case 'C': jump = offsetc; break;
      case 'G': jump = offsetg; break;
      case 'T': jump = offsett; break;
      default: jump = 0;
      }

      if (jump == 0) {
	debug2(printf(" => %c: No matches\n",c));
	return;
      } else {
	debug2(printf(" => %c: Going to offset %d\n",c,jump));
	search_left_single(&(*splicesite_i),&(*best_nmismatches),&(*nbest),
			   &(triecontents[-jump]),splicesites,
			   splicefrags_ref,splicefrags_alt,query,queryptr,
			   query_compress,genome_blocks,snp_blocks,
			   dibasep,cmetp,pos5,pos3,plusp,
			   max_mismatches,nmismatches,charpos-1);
	return;
      }

    } else {
      debug2(printf("\n"));
      if (offseta > 0) {
	debug2(printf("Branching to A: offset %d\n",offseta));
	search_left_single(&(*splicesite_i),&(*best_nmismatches),&(*nbest),
			   &(triecontents[-offseta]),splicesites,
			   splicefrags_ref,splicefrags_alt,query,queryptr,
			   query_compress,genome_blocks,snp_blocks,
			   dibasep,cmetp,pos5,pos3,plusp,
			   max_mismatches,nmismatches+(c != 'A'),charpos-1);
      }
      
      if (offsetc > 0) {
	debug2(printf("Branching to C: offset %d\n",offsetc));
	search_left_single(&(*splicesite_i),&(*best_nmismatches),&(*nbest),
			   &(triecontents[-offsetc]),splicesites,
			   splicefrags_ref,splicefrags_alt,query,queryptr,
			   query_compress,genome_blocks,snp_blocks,
			   dibasep,cmetp,pos5,pos3,plusp,
			   max_mismatches,nmismatches+(c != 'C'),charpos-1);
      }

      if (offsetg > 0) {
	debug2(printf("Branching to G: offset %d\n",offsetg));
	search_left_single(&(*splicesite_i),&(*best_nmismatches),&(*nbest),
			   &(triecontents[-offsetg]),splicesites,
			   splicefrags_ref,splicefrags_alt,query,queryptr,
			   query_compress,genome_blocks,snp_blocks,
			   dibasep,cmetp,pos5,pos3,plusp,
			   max_mismatches,nmismatches+(c != 'G'),charpos-1);
      }

      if (offsett > 0) {
	debug2(printf("Branching to T: offset %d\n",offsett));
	search_left_single(&(*splicesite_i),&(*best_nmismatches),&(*nbest),
			   &(triecontents[-offsett]),splicesites,
			   splicefrags_ref,splicefrags_alt,query,queryptr,
			   query_compress,genome_blocks,snp_blocks,
			   dibasep,cmetp,pos5,pos3,plusp,
			   max_mismatches,nmismatches+(c != 'T'),charpos-1);
      }

      return;
    }
  }
}

#endif

#ifdef OLD_CODE

static void
search_right_single (int *splicesite_i, int *best_nmismatches, int *nbest,
		     unsigned int *triecontents, Genomicpos_T *splicesites,
		     UINT4 *splicefrags_ref, UINT4 *splicefrags_alt, char *query, char *queryptr,
		     Compress_T query_compress, UINT4 *genome_blocks, UINT4 *snp_blocks,
		     bool dibasep, bool cmetp, int pos5, int pos3, bool plusp,
		     int max_mismatches, int nmismatches, int charpos) {
  Genomicpos_T leaf, segment_left;
  int nleaves, i;
  UINT4 query_shifted, flags, mask;
  int offseta, offsetc, offsetg, offsett, jump;
  int ncolordiffs;
  char c;
  
  debug2(printf("Entered search_right_single with nmismatches %d, charpos %d, best_nmismatches %d, nbest %d, bestj %d\n",
		nmismatches,charpos,*best_nmismatches,*nbest,*splicesite_i));

  if (single_leaf_p(leaf = triecontents[0])) {

    if (charpos + 1 < pos3) {
      if (pos3 - pos5 <= 16) {
	/* Recomputes entire segment to determine mismatches (necessary because of splicefrags) */
	query_shifted = Genome_query_shift_fragment_right(&flags,&mask,query_compress,pos5,pos3);
	nmismatches = Genome_count_mismatches_fragment(query_shifted,flags,mask,
						       splicefrags_ref[leaf],splicefrags_alt[leaf]);
      } else {
	/* Can happen in search for short middle exon */
	segment_left = splicesites[leaf] - pos5;
	nmismatches =
	  Genome_count_mismatches_substring(&ncolordiffs,query,query_compress,genome_blocks,snp_blocks,
					    segment_left,pos5,pos3,dibasep,cmetp,plusp);
      }

      debug2(printf("Found leaf %u => %d mismatches\n",leaf,nmismatches));
      if (nmismatches < *best_nmismatches) {
	debug2(printf("  nmismatches %d < *best_nmismatches %d => setting bestj to be %d (new best_nmismatches %d)\n",
		      nmismatches,*best_nmismatches,(int) leaf,nmismatches));
	*best_nmismatches = nmismatches;
	*splicesite_i = (int) leaf;
	*nbest = 1;
      } else if (nmismatches == *best_nmismatches) {
	*nbest += 1;
	debug2(printf("  nmismatches %d == *best_nmismatches %d => incrementing nbest to be %d (same best_nmismatches %d)\n",
		      nmismatches,*best_nmismatches,*nbest,nmismatches));
      }
      
    } else {
      debug2(printf("Found leaf %u, and completely finished query, so unique\n",leaf));
      if (nmismatches < *best_nmismatches) {
	debug2(printf("  nmismatches %d < *best_nmismatches %d => setting bestj to be %d (new best_nmismatches %d)\n",
		      nmismatches,*best_nmismatches,(int) leaf,nmismatches));
	*best_nmismatches = nmismatches;
	*splicesite_i = (int) leaf;
	*nbest = 1;
      } else if (nmismatches == *best_nmismatches) {
	*nbest += 1;
	debug2(printf("  nmismatches %d == *best_nmismatches %d => incrementing nbest to be %d (same best_nmismatches %d)\n",
		      nmismatches,*best_nmismatches,*nbest,nmismatches));
      }
    }

    return;

  } else if (multiple_leaf_p(leaf)) {
    nleaves = (int) (-leaf);

    for (i = 1; i <= nleaves; i++) {
      leaf = triecontents[i];

      if (charpos + 1 < pos3) {
	if (pos3 - pos5 <= 16) {
	  /* Recomputes entire segment to determine mismatches (necessary because of splicefrags) */
	  query_shifted = Genome_query_shift_fragment_right(&flags,&mask,query_compress,pos5,pos3);
	  nmismatches = Genome_count_mismatches_fragment(query_shifted,flags,mask,
							 splicefrags_ref[leaf],splicefrags_alt[leaf]);
	} else {
	  /* Can happen in search for short middle exon */
	  segment_left = splicesites[leaf] - pos5;
	  nmismatches =
	    Genome_count_mismatches_substring(&ncolordiffs,query,query_compress,genome_blocks,snp_blocks,
					      segment_left,pos5,pos3,dibasep,cmetp,plusp);
	}

	debug2(printf("Found leaf %u => %d mismatches\n",leaf,nmismatches));
	if (nmismatches < *best_nmismatches) {
	  debug2(printf("  nmismatches %d < *best_nmismatches %d => setting bestj to be %d (new best_nmismatches %d)\n",
			nmismatches,*best_nmismatches,(int) leaf,nmismatches));
	  *best_nmismatches = nmismatches;
	  *splicesite_i = (int) leaf;
	  *nbest = 1;
	} else if (nmismatches == *best_nmismatches) {
	  *nbest += 1;
	  debug2(printf("  nmismatches %d == *best_nmismatches %d => incrementing nbest to be %d (same best_nmismatches %d)\n",
			nmismatches,*best_nmismatches,*nbest,nmismatches));
	}
      
      } else {
	debug2(printf("Found leaf %u, and completely finished query, so unique\n",leaf));
	if (nmismatches < *best_nmismatches) {
	  debug2(printf("  nmismatches %d < *best_nmismatches %d => setting bestj to be %d (new best_nmismatches %d)\n",
			nmismatches,*best_nmismatches,(int) leaf,nmismatches));
	  *best_nmismatches = nmismatches;
	  *splicesite_i = (int) leaf;
	  *nbest = 1;
	} else if (nmismatches == *best_nmismatches) {
	  *nbest += 1;
	  debug2(printf("  nmismatches %d == *best_nmismatches %d => incrementing nbest to be %d (same best_nmismatches %d)\n",
			nmismatches,*best_nmismatches,*nbest,nmismatches));
	}
      }
    }

    return;

  } else if (charpos + 1 >= pos3) {
    if (splicefrags_alt == splicefrags_ref) {
      /* Non-leaf, and no more characters left, so have multiple hits */
      debug2(printf("Non-leaf, but no more characters left: nmismatches %d vs best_nmismatches %d\n",
		    nmismatches,*best_nmismatches));
      if (nmismatches < *best_nmismatches) {
	debug2(printf("  nmismatches %d < *best_nmismatches %d => setting nbest to be %d (new best_nmismatches %d)\n",
		      nmismatches,*best_nmismatches,2,nmismatches));
	*best_nmismatches = nmismatches;
	/* *splicesite_i = (int) leaf; -- What should we set splicesite_i to be? */
	*nbest = 2;
      } else if (nmismatches == *best_nmismatches) {
	*nbest += 2;
	debug2(printf("  nmismatches %d == *best_nmismatches %d => incrementing nbest to be %d (same best_nmismatches %d)\n",
		      nmismatches,*best_nmismatches,*nbest,nmismatches));
      }
      return;

    } else {
      /* Alternate genome.  Need to re-check entire subtree, because nmismatches is just a lower bound */
      debug2(printf("Non-leaf, no more characters left, so check all splicefrags\n"));
      search_subtree_right(&(*splicesite_i),&(*best_nmismatches),&(*nbest),
			   triecontents,splicesites,splicefrags_ref,splicefrags_alt,
			   query,query_compress,genome_blocks,snp_blocks,
			   dibasep,cmetp,pos5,pos3,plusp);
      return;
    }

  } else {
    /* Non-leaf, and characters left, so recurse */
    c = queryptr[charpos + 1];

#ifdef USE_2BYTE_RELOFFSETS
    get_offsets(&offseta,&offsetc,&offsetg,&offsett,triecontents[1],triecontents[2]);
#else
    offseta = (int) triecontents[1];
    offsetc = (int) triecontents[2];
    offsetg = (int) triecontents[3];
    offsett = (int) triecontents[4];
#endif

    debug2(printf("Got offsets %d %d %d %d\n",offseta,offsetc,offsetg,offsett));
    if (nmismatches == max_mismatches) {
      /* Continue only with a match */
      switch (c) {
      case 'A': jump = offseta; break;
      case 'C': jump = offsetc; break;
      case 'G': jump = offsetg; break;
      case 'T': jump = offsett; break;
      default: jump = 0;
      }

      if (jump == 0) {
	debug2(printf("No matches\n"));
	return;
      } else {
	debug2(printf("Going to offset %d\n",jump));
	search_right_single(&(*splicesite_i),&(*best_nmismatches),&(*nbest),
			    &(triecontents[-jump]),splicesites,
			    splicefrags_ref,splicefrags_alt,query,queryptr,
			    query_compress,genome_blocks,snp_blocks,
			    dibasep,cmetp,pos5,pos3,plusp,
			    max_mismatches,nmismatches,charpos+1);
	return;
      }

    } else {
      if (offseta > 0) {
	debug2(printf("Branching to A: offset %d\n",offseta));
	search_right_single(&(*splicesite_i),&(*best_nmismatches),&(*nbest),
			    &(triecontents[-offseta]),splicesites,
			    splicefrags_ref,splicefrags_alt,query,queryptr,
			    query_compress,genome_blocks,snp_blocks,
			    dibasep,cmetp,pos5,pos3,plusp,
			    max_mismatches,nmismatches+(c != 'A'),charpos+1);
      }
      
      if (offsetc > 0) {
	debug2(printf("Branching to C: offset %d\n",offsetc));
	search_right_single(&(*splicesite_i),&(*best_nmismatches),&(*nbest),
			    &(triecontents[-offsetc]),splicesites,
			    splicefrags_ref,splicefrags_alt,query,queryptr,
			    query_compress,genome_blocks,snp_blocks,
			    dibasep,cmetp,pos5,pos3,plusp,
			    max_mismatches,nmismatches+(c != 'C'),charpos+1);
      }

      if (offsetg > 0) {
	debug2(printf("Branching to G: offset %d\n",offsetg));
	search_right_single(&(*splicesite_i),&(*best_nmismatches),&(*nbest),
			    &(triecontents[-offsetg]),splicesites,
			    splicefrags_ref,splicefrags_alt,query,queryptr,
			    query_compress,genome_blocks,snp_blocks,
			    dibasep,cmetp,pos5,pos3,plusp,
			    max_mismatches,nmismatches+(c != 'G'),charpos+1);
      }

      if (offsett > 0) {
	debug2(printf("Branching to T: offset %d\n",offsett));
	search_right_single(&(*splicesite_i),&(*best_nmismatches),&(*nbest),
			    &(triecontents[-offsett]),splicesites,
			    splicefrags_ref,splicefrags_alt,query,queryptr,
			    query_compress,genome_blocks,snp_blocks,
			    dibasep,cmetp,pos5,pos3,plusp,
			    max_mismatches,nmismatches+(c != 'T'),charpos+1);
      }

      return;
    }
  }
}

#endif


/************************************************************************
 *   General search procedures
 ************************************************************************/

static Intlist_T
search_left (int *best_nmismatches, Intlist_T *nmismatches_list, Intlist_T splicesites_i,
	     unsigned int *triecontents, Genomicpos_T *splicesites,
	     UINT4 *splicefrags_ref, UINT4 *splicefrags_alt, char *query, char *queryptr,
	     Compress_T query_compress, UINT4 *genome_blocks, UINT4 *snp_blocks,
	     bool dibasep, bool cmetp, int pos5, int pos3, bool plusp,
	     int nmismatches, int charpos) {
  Genomicpos_T leaf, segment_left;
  int nleaves, i;
  UINT4 query_shifted, flags, mask;
  int offseta, offsetc, offsetg, offsett;
  int ncolordiffs;
  char c;
  
  debug2(printf("Entered search_left with nmismatches %d, charpos %d, best_nmismatches %d\n",
		nmismatches,charpos,*best_nmismatches));

  if (nmismatches > *best_nmismatches) {
    return splicesites_i;

  } else if (single_leaf_p(leaf = triecontents[0])) {

    if (charpos - 1 >= pos5) {
      if (pos3 - pos5 <= 16) {
	/* Recomputes entire segment to determine mismatches (necessary because of splicefrags) */
	query_shifted = Genome_query_shift_fragment_left(&flags,&mask,query_compress,pos5,pos3);
	nmismatches = Genome_count_mismatches_fragment(query_shifted,flags,mask,
						       splicefrags_ref[leaf],splicefrags_alt[leaf]);
      } else {
	/* Can happen in search for short middle exon */
	segment_left = splicesites[leaf] - pos3;
	nmismatches =
	  Genome_count_mismatches_substring(&ncolordiffs,query,query_compress,genome_blocks,snp_blocks,
					    segment_left,pos5,pos3,dibasep,cmetp,plusp);
      }

      debug2(printf("Found leaf %u, but still have characters to check against Genome: %.*s => %d mismatches\n",
		    leaf,charpos + 1,&(queryptr[0]),nmismatches));
      if (nmismatches < *best_nmismatches) {
	debug2(printf("  nmismatches %d < best_nmismatches %d => setting splicesites_i to be %d (new best_nmismatches %d)\n",
		      nmismatches,*best_nmismatches,(int) leaf,nmismatches));
	*best_nmismatches = nmismatches;
	Intlist_free(&(*nmismatches_list));
	*nmismatches_list = Intlist_push(NULL,nmismatches);
	Intlist_free(&splicesites_i);
	return Intlist_push(NULL,(int) leaf);
      } else if (nmismatches == *best_nmismatches) {
	debug2(printf("  nmismatches %d == best_nmismatches %d => pushing leaf %d (same best_nmismatches %d)\n",
		      nmismatches,*best_nmismatches,(int) leaf,nmismatches));
	*nmismatches_list = Intlist_push(*nmismatches_list,nmismatches);
	return Intlist_push(splicesites_i,(int) leaf);
      } else {
	return splicesites_i;
      }
      
    } else {
      debug2(printf("Found leaf %u, and completely finished query, so unique\n",leaf));
      if (nmismatches < *best_nmismatches) {
	debug2(printf("  nmismatches %d < best_nmismatches %d => setting splicesites_i to be %d (new best_nmismatches %d)\n",
		      nmismatches,*best_nmismatches,(int) leaf,nmismatches));
	*best_nmismatches = nmismatches;
	Intlist_free(&(*nmismatches_list));
	*nmismatches_list = Intlist_push(NULL,nmismatches);
	Intlist_free(&splicesites_i);
	return Intlist_push(NULL,(int) leaf);
      } else if (nmismatches == *best_nmismatches) {
	debug2(printf("  nmismatches %d == best_nmismatches %d => pushing leaf %d (same best_nmismatches %d)\n",
		      nmismatches,*best_nmismatches,(int) leaf,nmismatches));
	*nmismatches_list = Intlist_push(*nmismatches_list,nmismatches);
	return Intlist_push(splicesites_i,(int) leaf);
      } else {
	return splicesites_i;
      }
    }

  } else if (multiple_leaf_p(leaf)) {
    nleaves = (int) (-leaf);

    /* Need to compute nmismatches multiple times if we have SNPs */
    for (i = 1; i <= nleaves; i++) {
      leaf = triecontents[i];

      if (charpos - 1 >= pos5) {
	if (pos3 - pos5 <= 16) {
	  /* Recomputes entire segment to determine mismatches (necessary because of splicefrags) */
	  query_shifted = Genome_query_shift_fragment_left(&flags,&mask,query_compress,pos5,pos3);
	  nmismatches = Genome_count_mismatches_fragment(query_shifted,flags,mask,
							 splicefrags_ref[leaf],splicefrags_alt[leaf]);
	} else {
	  /* Can happen in search for short middle exon */
	  segment_left = splicesites[leaf] - pos3;
	  nmismatches =
	    Genome_count_mismatches_substring(&ncolordiffs,query,query_compress,genome_blocks,snp_blocks,
					      segment_left,pos5,pos3,dibasep,cmetp,plusp);
	}

	debug2(printf("Found leaf %u => %d mismatches\n",leaf,nmismatches));
	if (nmismatches < *best_nmismatches) {
	  debug2(printf("  nmismatches %d < best_nmismatches %d => setting splicesites_i to be %d (new best_nmismatches %d)\n",
			nmismatches,*best_nmismatches,(int) leaf,nmismatches));
	  *best_nmismatches = nmismatches;
	  Intlist_free(&(*nmismatches_list));
	  *nmismatches_list = Intlist_push(NULL,nmismatches);
	  Intlist_free(&splicesites_i);
	  splicesites_i = Intlist_push(NULL,(int) leaf);
	} else if (nmismatches == *best_nmismatches) {
	  debug2(printf("  nmismatches %d == best_nmismatches %d => pushing leaf %d (same best_nmismatches %d)\n",
			nmismatches,*best_nmismatches,(int) leaf,nmismatches));
	  *nmismatches_list = Intlist_push(*nmismatches_list,nmismatches);
	  splicesites_i = Intlist_push(splicesites_i,(int) leaf);
	}
      
      } else {
	debug2(printf("Found leaf %u, and completely finished query, so unique\n",leaf));
	if (nmismatches < *best_nmismatches) {
	  debug2(printf("  nmismatches %d < best_nmismatches %d => setting splicesites_i to be %d (new best_nmismatches %d)\n",
			nmismatches,*best_nmismatches,(int) leaf,nmismatches));
	  *best_nmismatches = nmismatches;
	  Intlist_free(&(*nmismatches_list));
	  *nmismatches_list = Intlist_push(NULL,nmismatches);
	  Intlist_free(&splicesites_i);
	  splicesites_i = Intlist_push(NULL,(int) leaf);
	} else if (nmismatches == *best_nmismatches) {
	  debug2(printf("  nmismatches %d == best_nmismatches %d => pushing leaf %d (same best_nmismatches %d)\n",
			nmismatches,*best_nmismatches,(int) leaf,nmismatches));
	  *nmismatches_list = Intlist_push(*nmismatches_list,nmismatches);
	  splicesites_i = Intlist_push(splicesites_i,(int) leaf);
	}
      }
    }

    return splicesites_i;

  } else {
    /* Non-leaf, and characters left, so recurse */
#ifdef USE_2BYTE_RELOFFSETS
    get_offsets(&offseta,&offsetc,&offsetg,&offsett,triecontents[1],triecontents[2]);
#else
    offseta = (int) triecontents[1];
    offsetc = (int) triecontents[2];
    offsetg = (int) triecontents[3];
    offsett = (int) triecontents[4];
#endif

    debug2(printf("Got offsets %d %d %d %d\n",offseta,offsetc,offsetg,offsett));

    if (charpos - 1 >= pos5) {
      c = queryptr[charpos - 1];
    } else {
      c = '\0';
    }

    if (offseta > 0) {
      debug2(printf("Branching to A: offset %d\n",offseta));
      splicesites_i = search_left(&(*best_nmismatches),&(*nmismatches_list),splicesites_i,
				  &(triecontents[-offseta]),splicesites,
				  splicefrags_ref,splicefrags_alt,query,queryptr,
				  query_compress,genome_blocks,snp_blocks,
				  dibasep,cmetp,pos5,pos3,plusp,
				  nmismatches+(c != '\0' && c != 'A') ,charpos-1);
    }
      
    if (offsetc > 0) {
      debug2(printf("Branching to C: offset %d\n",offsetc));
      splicesites_i = search_left(&(*best_nmismatches),&(*nmismatches_list),splicesites_i,
				  &(triecontents[-offsetc]),splicesites,
				  splicefrags_ref,splicefrags_alt,query,queryptr,
				  query_compress,genome_blocks,snp_blocks,
				  dibasep,cmetp,pos5,pos3,plusp,
				  nmismatches+(c != '\0' && c != 'C'),charpos-1);
    }

    if (offsetg > 0) {
      debug2(printf("Branching to G: offset %d\n",offsetg));
      splicesites_i = search_left(&(*best_nmismatches),&(*nmismatches_list),splicesites_i,
				  &(triecontents[-offsetg]),splicesites,
				  splicefrags_ref,splicefrags_alt,query,queryptr,
				  query_compress,genome_blocks,snp_blocks,
				  dibasep,cmetp,pos5,pos3,plusp,
				  nmismatches+(c != '\0' && c != 'G'),charpos-1);
    }

    if (offsett > 0) {
      debug2(printf("Branching to T: offset %d\n",offsett));
      splicesites_i = search_left(&(*best_nmismatches),&(*nmismatches_list),splicesites_i,
				  &(triecontents[-offsett]),splicesites,
				  splicefrags_ref,splicefrags_alt,query,queryptr,
				  query_compress,genome_blocks,snp_blocks,
				  dibasep,cmetp,pos5,pos3,plusp,
				  nmismatches+(c != '\0' && c != 'T'),charpos-1);
    }

    return splicesites_i;
  }
}



static Intlist_T
search_right (int *best_nmismatches, Intlist_T *nmismatches_list, Intlist_T splicesites_i,
	      unsigned int *triecontents, Genomicpos_T *splicesites,
	      UINT4 *splicefrags_ref, UINT4 *splicefrags_alt, char *query, char *queryptr,
	      Compress_T query_compress, UINT4 *genome_blocks, UINT4 *snp_blocks,
	      bool dibasep, bool cmetp, int pos5, int pos3, bool plusp,
	      int nmismatches, int charpos) {
  Genomicpos_T leaf, segment_left;
  int nleaves, i;
  UINT4 query_shifted, flags, mask;
  int offseta, offsetc, offsetg, offsett;
  int ncolordiffs;
  char c;
  
  debug2(printf("Entered search_right with nmismatches %d, charpos %d, best_nmismatches %d\n",
		nmismatches,charpos,*best_nmismatches));

  if (nmismatches > *best_nmismatches) {
    return splicesites_i;

  } else if (single_leaf_p(leaf = triecontents[0])) {

    if (charpos + 1 < pos3) {
      if (pos3 - pos5 <= 16) {
	/* Recomputes entire segment to determine mismatches (necessary because of splicefrags) */
	query_shifted = Genome_query_shift_fragment_right(&flags,&mask,query_compress,pos5,pos3);
	nmismatches = Genome_count_mismatches_fragment(query_shifted,flags,mask,
						       splicefrags_ref[leaf],splicefrags_alt[leaf]);
      } else {
	/* Can happen in search for short middle exon */
	segment_left = splicesites[leaf] - pos5;
	nmismatches =
	  Genome_count_mismatches_substring(&ncolordiffs,query,query_compress,genome_blocks,snp_blocks,
					    segment_left,pos5,pos3,dibasep,cmetp,plusp);
      }

      debug2(printf("Found leaf %u => %d mismatches\n",leaf,nmismatches));
      if (nmismatches < *best_nmismatches) {
	debug2(printf("  nmismatches %d < *best_nmismatches %d => setting splicesites_i to be %d (new best_nmismatches %d)\n",
		      nmismatches,*best_nmismatches,(int) leaf,nmismatches));
	*best_nmismatches = nmismatches;
	Intlist_free(&(*nmismatches_list));
	*nmismatches_list = Intlist_push(NULL,nmismatches);
	Intlist_free(&splicesites_i);
	return Intlist_push(NULL,(int) leaf);
      } else if (nmismatches == *best_nmismatches) {
	debug2(printf("  nmismatches %d == *best_nmismatches %d => pushing leaf %d (same best_nmismatches %d)\n",
		      nmismatches,*best_nmismatches,(int) leaf,nmismatches));
	*nmismatches_list = Intlist_push(*nmismatches_list,nmismatches);
	return Intlist_push(splicesites_i,(int) leaf);
      } else {
	return splicesites_i;
      }
      
    } else {
      debug2(printf("Found leaf %u, and completely finished query, so unique\n",leaf));
      if (nmismatches < *best_nmismatches) {
	debug2(printf("  nmismatches %d < *best_nmismatches %d => setting splicesites_i to be %d (new best_nmismatches %d)\n",
		      nmismatches,*best_nmismatches,(int) leaf,nmismatches));
	*best_nmismatches = nmismatches;
	Intlist_free(&(*nmismatches_list));
	*nmismatches_list = Intlist_push(NULL,nmismatches);
	Intlist_free(&splicesites_i);
	return Intlist_push(NULL,(int) leaf);
      } else if (nmismatches == *best_nmismatches) {
	debug2(printf("  nmismatches %d == *best_nmismatches %d => pushing leaf %d (same best_nmismatches %d)\n",
		      nmismatches,*best_nmismatches,(int) leaf,nmismatches));
	*nmismatches_list = Intlist_push(*nmismatches_list,nmismatches);
	return Intlist_push(splicesites_i,(int) leaf);
      } else {
	return splicesites_i;
      }
    }

  } else if (multiple_leaf_p(leaf)) {
    nleaves = (int) (-leaf);

    /* Need to compute nmismatches multiple times if we have SNPs */
    for (i = 1; i <= nleaves; i++) {
      leaf = triecontents[i];

      if (charpos + 1 < pos3) {
	if (pos3 - pos5 <= 16) {
	  /* Recomputes entire segment to determine mismatches (necessary because of splicefrags) */
	  query_shifted = Genome_query_shift_fragment_right(&flags,&mask,query_compress,pos5,pos3);
	  nmismatches = Genome_count_mismatches_fragment(query_shifted,flags,mask,
							 splicefrags_ref[leaf],splicefrags_alt[leaf]);
	} else {
	  /* Can happen in search for short middle exon */
	  segment_left = splicesites[leaf] - pos5;
	  nmismatches =
	    Genome_count_mismatches_substring(&ncolordiffs,query,query_compress,genome_blocks,snp_blocks,
					      segment_left,pos5,pos3,dibasep,cmetp,plusp);
	}

	debug2(printf("Found leaf %u => %d mismatches\n",leaf,nmismatches));
	if (nmismatches < *best_nmismatches) {
	  debug2(printf("  nmismatches %d < *best_nmismatches %d => setting splicesites_i to be %d (new best_nmismatches %d)\n",
			nmismatches,*best_nmismatches,(int) leaf,nmismatches));
	  *best_nmismatches = nmismatches;
	  Intlist_free(&(*nmismatches_list));
	  *nmismatches_list = Intlist_push(NULL,nmismatches);
	  Intlist_free(&splicesites_i);
	  splicesites_i =  Intlist_push(NULL,(int) leaf);
	} else if (nmismatches == *best_nmismatches) {
	  debug2(printf("  nmismatches %d == *best_nmismatches %d => pushing leaf %d (same best_nmismatches %d)\n",
			nmismatches,*best_nmismatches,(int) leaf,nmismatches));
	  *nmismatches_list = Intlist_push(*nmismatches_list,nmismatches);
	  splicesites_i = Intlist_push(splicesites_i,(int) leaf);
	}
      
      } else {
	debug2(printf("Found leaf %u, and completely finished query, so unique\n",leaf));
	if (nmismatches < *best_nmismatches) {
	  debug2(printf("  nmismatches %d < *best_nmismatches %d => setting splicesites_i to be %d (new best_nmismatches %d)\n",
			nmismatches,*best_nmismatches,(int) leaf,nmismatches));
	  *best_nmismatches = nmismatches;
	  Intlist_free(&(*nmismatches_list));
	  *nmismatches_list = Intlist_push(NULL,nmismatches);
	  Intlist_free(&splicesites_i);
	  splicesites_i = Intlist_push(NULL,(int) leaf);
	} else if (nmismatches == *best_nmismatches) {
	  debug2(printf("  nmismatches %d == *best_nmismatches %d => pushing leaf %d (same best_nmismatches %d)\n",
			nmismatches,*best_nmismatches,(int) leaf,nmismatches));
	  *nmismatches_list = Intlist_push(*nmismatches_list,nmismatches);
	  splicesites_i = Intlist_push(splicesites_i,(int) leaf);
	}
      }
    }

    return splicesites_i;

  } else {
    /* Non-leaf, and characters left, so recurse */
#ifdef USE_2BYTE_RELOFFSETS
    get_offsets(&offseta,&offsetc,&offsetg,&offsett,triecontents[1],triecontents[2]);
#else
    offseta = (int) triecontents[1];
    offsetc = (int) triecontents[2];
    offsetg = (int) triecontents[3];
    offsett = (int) triecontents[4];
#endif

    debug2(printf("Got offsets %d %d %d %d\n",offseta,offsetc,offsetg,offsett));

    if (charpos + 1 < pos3) {
      c = queryptr[charpos + 1];
    } else {
      c = '\0';
    }

    if (offseta > 0) {
      debug2(printf("Branching to A: offset %d\n",offseta));
      splicesites_i = search_right(&(*best_nmismatches),&(*nmismatches_list),splicesites_i,
				   &(triecontents[-offseta]),splicesites,
				   splicefrags_ref,splicefrags_alt,query,queryptr,
				   query_compress,genome_blocks,snp_blocks,
				   dibasep,cmetp,pos5,pos3,plusp,
				   nmismatches+(c != '\0' && c != 'A'),charpos+1);
    }
      
    if (offsetc > 0) {
      debug2(printf("Branching to C: offset %d\n",offsetc));
      splicesites_i = search_right(&(*best_nmismatches),&(*nmismatches_list),splicesites_i,
				   &(triecontents[-offsetc]),splicesites,
				   splicefrags_ref,splicefrags_alt,query,queryptr,
				   query_compress,genome_blocks,snp_blocks,
				   dibasep,cmetp,pos5,pos3,plusp,
				   nmismatches+(c != '\0' && c != 'C'),charpos+1);
    }

    if (offsetg > 0) {
      debug2(printf("Branching to G: offset %d\n",offsetg));
      splicesites_i = search_right(&(*best_nmismatches),&(*nmismatches_list),splicesites_i,
				   &(triecontents[-offsetg]),splicesites,
				   splicefrags_ref,splicefrags_alt,query,queryptr,
				   query_compress,genome_blocks,snp_blocks,
				   dibasep,cmetp,pos5,pos3,plusp,
				   nmismatches+(c != '\0' && c != 'G'),charpos+1);
    }

    if (offsett > 0) {
      debug2(printf("Branching to T: offset %d\n",offsett));
      splicesites_i = search_right(&(*best_nmismatches),&(*nmismatches_list),splicesites_i,
				   &(triecontents[-offsett]),splicesites,
				   splicefrags_ref,splicefrags_alt,query,queryptr,
				   query_compress,genome_blocks,snp_blocks,
				   dibasep,cmetp,pos5,pos3,plusp,
				   nmismatches+(c != '\0' && c != 'T'),charpos+1);
    }

    return splicesites_i;
  }
}


#ifdef SUBOPTIMAL_SEPARATION

bool
Splicetrie_find_left_short (bool *ambp, int *best_nmismatches, int *bestj,
			    int i, int pos5, int pos3, 
			    Genomicpos_T *splicesites, UINT4 *splicefrags_ref, UINT4 *splicefrags_alt,
			    unsigned int *trieoffsets, unsigned int *triecontents,
			    Compress_T query_compress, char *queryptr, bool dibasep, bool cmetp) {
  int nhits0, nhits1, nhits2;
  int bestj0, bestj1, bestj2;

  debug1(printf("Queryptr is %s from %d to %d: %.*s\n",queryptr,pos5,pos3,pos3-pos5,&(queryptr[pos5])));

  debug1(printf("**Searching trie with max_mismatches 0, maxhits 1\n"));
  nhits0 = Trie_count_left(&bestj0,/*exclude_i*/-1,&(triecontents[trieoffsets[i]]),
			   splicefrags_ref,splicefrags_alt,queryptr,
			   /*max_mismatches*/0,/*maxhits*/1,query_compress,
			   dibasep,cmetp,pos5,pos3,/*nmismatches*/0,/*charpos*/pos3);
  if (nhits0 > 1) {
    /* nhits0 > 1 => amb */
    *ambp = true;
    *best_nmismatches = 0;
    return true;

  } else if (nhits0 == 1) {
    debug1(printf("**Searching trie with max_mismatches 1, maxhits 0\n"));
    nhits1 = count_left(&bestj1,/*exclude_i*/bestj0,&(triecontents[trieoffsets[i]]),
			splicefrags_ref,splicefrags_alt,queryptr,
			/*max_mismatches*/1,/*maxhits*/0,query_compress,
			dibasep,cmetp,pos5,pos3,/*nmismatches*/0,/*charpos*/pos3);
    if (nhits1 > 0) {
      /* nhits0 == 1, nhits1 > 0 * => amb */
      *ambp = true;
      *best_nmismatches = 0;
      return true;
    } else {
      /* nhits0 == 1, nhits1 == 0 => unique */
      *ambp = false;
      *best_nmismatches = 0;
      *bestj = bestj0;
      return true;
    }

  } else {
    debug1(printf("**Searching trie with max_mismatches 1, maxhits 1\n"));
    nhits1 = count_left(&bestj1,/*exclude_i*/-1,&(triecontents[trieoffsets[i]]),
			splicefrags_ref,splicefrags_alt,queryptr,
			/*max_mismatches*/1,/*maxhits*/1,query_compress,
			dibasep,cmetp,pos5,pos3,/*nmismatches*/0,/*charpos*/pos3);
    if (nhits1 > 1) {
      /* nhits0 == 0, nhits1 > 1 => amb */
      *ambp = true;
      *best_nmismatches = 1;
      return true;

    } else if (nhits1 == 0) {
      /* nhits0 == 0, nhits1 == 0 => fail */
      *bestj = -1;
      return false;

    } else {
      debug1(printf("**Searching trie with max_mismatches 2, maxhits 0\n"));
      nhits2 = count_left(&bestj2,/*exclude_i*/bestj1,&(triecontents[trieoffsets[i]]),
			  splicefrags_ref,splicefrags_alt,queryptr,
			  /*max_mismatches*/2,/*maxhits*/0,query_compress,
			  dibasep,cmetp,pos5,pos3,/*nmismatches*/0,/*charpos*/pos3);
      if (nhits2 > 0) {
	/* nhits0 == 0, nhits1 == 1, nhits2 > 0 => amb */
	*ambp = true;
	*best_nmismatches = 1;
	return true;
      } else {
	/* nhits0 == 0, nhits1 == 1, nhits2 == 0 => unique */
	*ambp = false;
	*best_nmismatches = 1;
	*bestj = bestj1;
	return true;
      }
    }
  }
}

bool
Splicetrie_find_right_short (bool *ambp, int *best_nmismatches, int *bestj,
			     int i, int pos5, int pos3, 
			     Genomicpos_T *splicesites, UINT4 *splicefrags_ref, UINT4 *splicefrags_alt,
			     unsigned int *trieoffsets, unsigned int *triecontents,
			     Compress_T query_compress, char *queryptr, bool dibasep, bool cmetp) {
  int nhits0, nhits1, nhits2;
  int bestj0, bestj1, bestj2;
  
  debug1(printf("Queryptr is %s from %d to %d: %.*s\n",queryptr,pos5,pos3,pos3-pos5,&(queryptr[pos5])));

  debug1(printf("**Searching trie with max_mismatches 0, maxhits 1\n"));
  nhits0 = count_right(&bestj0,/*exclude_i*/-1,&(triecontents[trieoffsets[i]]),
		       splicefrags_ref,splicefrags_alt,queryptr,
		       /*max_mismatches*/0,/*maxhits*/1,query_compress,
		       dibasep,cmetp,pos5,pos3,/*nmismatches*/0,/*charpos*/pos5-1);
  if (nhits0 > 1) {
    /* nhits0 > 1 => amb */
    *ambp = true;
    *best_nmismatches = 0;
    return true;

  } else if (nhits0 == 1) {
    debug1(printf("**Searching trie with max_mismatches 1, maxhits 0\n"));
    nhits1 = count_right(&bestj1,/*exclude_i*/bestj0,&(triecontents[trieoffsets[i]]),
			 splicefrags_ref,splicefrags_alt,queryptr,
			 /*max_mismatches*/1,/*maxhits*/0,query_compress,
			 dibasep,cmetp,pos5,pos3,/*nmismatches*/0,/*charpos*/pos5-1);
    if (nhits1 > 0) {
      /* nhits0 == 1, nhits1 > 0 * => amb */
      *ambp = true;
      *best_nmismatches = 0;
      return true;
    } else {
      /* nhits0 == 1, nhits1 == 0 => unique */
      *ambp = false;
      *best_nmismatches = 0;
      *bestj = bestj0;
      return true;
    }

  } else {
    debug1(printf("**Searching trie with max_mismatches 1, maxhits 1\n"));
    nhits1 = Trie_count_right(&bestj1,/*exclude_i*/-1,&(triecontents[trieoffsets[i]]),
			      splicefrags_ref,splicefrags_alt,queryptr,
			      /*max_mismatches*/1,/*maxhits*/1,query_compress,
			      dibasep,cmetp,pos5,pos3,/*nmismatches*/0,/*charpos*/pos5-1);
    if (nhits1 > 1) {
      /* nhits0 == 0, nhits1 > 1 => amb */
      *ambp = true;
      *best_nmismatches = 1;
      return true;

    } else if (nhits1 == 0) {
      /* nhits0 == 0, nhits1 == 0 => fail */
      *bestj = -1;
      return false;

    } else {
      debug1(printf("**Searching trie with max_mismatches 2, maxhits 0\n"));
      nhits2 = count_right(&bestj2,/*exclude_i*/bestj1,&(triecontents[trieoffsets[i]]),
			   splicefrags_ref,splicefrags_alt,queryptr,
			   /*max_mismatches*/2,/*maxhits*/0,query_compress,
			   dibasep,cmetp,pos5,pos3,/*nmismatches*/0,/*charpos*/pos5-1);
      if (nhits2 > 0) {
	/* nhits0 == 0, nhits1 == 1, nhits2 > 0 => amb */
	*ambp = true;
	*best_nmismatches = 1;
	return true;
      } else {
	/* nhits0 == 0, nhits1 == 1, nhits2 == 0 => unique */
	*ambp = false;
	*best_nmismatches = 1;
	*bestj = bestj1;
	return true;
      }
    }
  }
}


#else

#ifdef OLD_CODE

bool
Splicetrie_find_left_single (bool *ambp, int *best_nmismatches, int *bestj,
			     int i, Genomicpos_T origleft, int pos5, int pos3,
			     Genomicpos_T *splicesites, UINT4 *splicefrags_ref, UINT4 *splicefrags_alt,
			     unsigned int *trieoffsets, unsigned int *triecontents,
			     Splicetype_T *splicetypes, int *nsplicepartners, List_T *splicestrings,
			     Compress_T query_compress, UINT4 *genome_blocks, UINT4 *snp_blocks,
			     char *query, char *queryptr, int max_mismatches_allowed,
			     bool dibasep, bool cmetp, bool plusp) {
  int nbest = 0;
  int ncolordiffs, nmatches;
  unsigned int *triebranch = NULL, *triestart;

  debug2(printf("Splicetrie_find_left called with #%d at origleft %u, pos5 %d, pos3 %d\n",i,origleft,pos5,pos3));

  if (pos5 == pos3) {
    debug2(printf("Splicetrie_find_left returning false, because pos5 == pos3\n"));
    return false;
  }

#if 0
  *best_nmismatches = pos3 - pos5;
#else
  /* Extension.  Cannot use splicefrags. */
  *best_nmismatches = 
    Genome_count_mismatches_substring(&ncolordiffs,query,query_compress,genome_blocks,snp_blocks,
				      origleft,pos5,pos3,dibasep,cmetp,plusp);
  *bestj = -1;
  nbest = 1;
  debug2(printf("  extension at %u => %d nmismatches\n",origleft,*best_nmismatches));
#endif

  if (triecontents) {
    triestart = &(triecontents[trieoffsets[i]]);
  } else {
    triebranch = Splicetrie_build_one(&triestart,i,splicetypes,nsplicepartners,splicestrings);
  }

  debug2(Splicetrie_dump(triestart,splicesites,splicetypes[i],splicefrags_ref));

  search_left_single(&(*bestj),&(*best_nmismatches),&nbest,triestart,splicesites,
		     splicefrags_ref,splicefrags_alt,query,queryptr,
		     query_compress,genome_blocks,snp_blocks,
		     dibasep,cmetp,pos5,pos3,plusp,max_mismatches_allowed,
		     /*nmismatches*/0,/*charpos*/pos3);

  if (triecontents == NULL) {
    FREE(triebranch);
  }

  debug2(printf("nbest = %d, *bestj = %d, *best_nmismatches = %d => ",
		nbest,*bestj,*best_nmismatches));
  if (nbest > 1) {
    *ambp = true;
    debug2(printf("Splicetrie_find_left returning true, with ambp true\n"));
    return true;
  } else if (*bestj < 0) {
    debug2(printf("Splicetrie_find_left returning false\n"));
    return false;
  } else if (nbest == 1) {

#if 0
    if (pos3 - pos5 > 16) {
      *ambp = false;
      debug2(printf("Splicetrie_find_left returning true, with ambp false\n"));
      return true;
    } else if (*best_nmismatches > MAX_BEST_NMISMATCHES) {
      debug2(printf("Splicetrie_find_left returning false\n"));
      return false;
    } else {
      *ambp = false;
      debug2(printf("Splicetrie_find_left returning true, with ambp false\n"));
      return true;
    }
#else
    nmatches = pos3 - pos5 - (*best_nmismatches);
    debug2(printf("nmatches %d => ",nmatches));
    if (3 * (*best_nmismatches) >= nmatches) {
      debug2(printf("Splicetrie_find_left returning false\n"));
      return false;
    } else {
      *ambp = false;
      debug2(printf("Splicetrie_find_left returning true, with ambp false\n"));
      return true;
    }
#endif

  } else {
    debug2(printf("Splicetrie_find_left returning false\n"));
    return false;
  }
}

#endif

#ifdef OLD_CODE

bool
Splicetrie_find_right_single (bool *ambp, int *best_nmismatches, int *bestj,
			      int i, Genomicpos_T origleft, int pos5, int pos3,
			      Genomicpos_T *splicesites, UINT4 *splicefrags_ref, UINT4 *splicefrags_alt,
			      unsigned int *trieoffsets, unsigned int *triecontents,
			      Splicetype_T *splicetypes, int *nsplicepartners, List_T *splicestrings,
			      Compress_T query_compress, UINT4 *genome_blocks, UINT4 *snp_blocks,
			      char *query, char *queryptr, int max_mismatches_allowed,
			      bool dibasep, bool cmetp, bool plusp) {
  int nbest = 0;
  int ncolordiffs, nmatches;
  unsigned int *triebranch = NULL, *triestart;

  debug2(printf("Splicetrie_find_right called with #%d at origleft %u, pos5 %d, pos3 %d\n",i,origleft,pos5,pos3));

  if (pos5 == pos3) {
    debug2(printf("Splicetrie_find_right returning false, because pos5 == pos3\n"));
    return false;
  }

#if 0
  *best_nmismatches = pos3 - pos5;
#else
  *best_nmismatches = 
    Genome_count_mismatches_substring(&ncolordiffs,query,query_compress,genome_blocks,snp_blocks,
				      origleft,pos5,pos3,dibasep,cmetp,plusp);
  *bestj = -1;
  nbest = 1;
  debug2(printf("  extension at %u => %d nmismatches\n",origleft,*best_nmismatches));
#endif

  if (triecontents) {
    triestart = &(triecontents[trieoffsets[i]]);
  } else {
    triebranch = Splicetrie_build_one(&triestart,i,splicetypes,nsplicepartners,splicestrings);
  }

  debug2(Splicetrie_dump(triestart,splicesites,splicetypes[i],splicefrags_ref));

  search_right_single(&(*bestj),&(*best_nmismatches),&nbest,triestart,splicesites,
		      splicefrags_ref,splicefrags_alt,query,queryptr,
		      query_compress,genome_blocks,snp_blocks,
		      dibasep,cmetp,pos5,pos3,plusp,max_mismatches_allowed,
		      /*nmismatches*/0,/*charpos*/pos5-1);

  if (triecontents == NULL) {
    FREE(triebranch);
  }

  debug2(printf("nbest = %d, *bestj = %d, *best_nmismatches = %d => ",
		nbest,*bestj,*best_nmismatches));
  if (nbest > 1) {
    *ambp = true;
    debug2(printf("Splicetrie_find_right returning true, with ambp true\n"));
    return true;
  } else if (*bestj < 0) {
    debug2(printf("Splicetrie_find_right returning false\n"));
    return false;
  } else if (nbest == 1) {

#if 0
    if (pos3 - pos5 > 16) {
      *ambp = false;
      debug2(printf("Splicetrie_find_right returning true, with ambp false\n"));
      return true;
    } else if (*best_nmismatches > MAX_BEST_NMISMATCHES) {
      debug2(printf("Splicetrie_find_right returning false\n"));
      return false;
    } else {
      *ambp = false;
      debug2(printf("Splicetrie_find_right returning true, with ambp false\n"));
      return true;
    }
#else
    nmatches = pos3 - pos5 - (*best_nmismatches);
    debug2(printf("nmatches %d => ",nmatches));
    if (3 * (*best_nmismatches) >= nmatches) {
      debug2(printf("Splicetrie_find_right returning false\n"));
      return false;
    } else {
      *ambp = false;
      debug2(printf("Splicetrie_find_right returning true, with ambp false\n"));
      return true;
    }
#endif

  } else {
    debug2(printf("Splicetrie_find_right returning false\n"));
    return false;
  }
}

#endif

#ifdef OLD_CODE

bool
Splicetrie_find_left_short (bool *ambp, int *best_nmismatches, int *bestj,
			    int i, Genomicpos_T origleft, int pos5, int pos3, 
			    Genomicpos_T *splicesites, UINT4 *splicefrags_ref, UINT4 *splicefrags_alt,
			    unsigned int *trieoffsets, unsigned int *triecontents,
			    Splicetype_T *splicetypes, int *nsplicepartners, List_T *splicestrings,
			    Compress_T query_compress, UINT4 *genome_blocks, UINT4 *snp_blocks,
			    char *query, char *queryptr, bool dibasep, bool cmetp, bool plusp) {
  int ext_nmismatches, ncolordiffs;
  int nhits0, nhits1;
  int bestj0, bestj1;
  unsigned int *triebranch = NULL, *triestart;

  debug1(printf("Queryptr is %s from %d to %d: %.*s\n",queryptr,pos5,pos3,pos3-pos5,&(queryptr[pos5])));

  assert(pos5 != pos3);

  if (triecontents) {
    triestart = &(triecontents[trieoffsets[i]]);
  } else {
    triebranch = Splicetrie_build_one(&triestart,i,splicetypes,nsplicepartners,splicestrings);
  }

#if 0
  debug2(Splicetrie_dump(triestart,splicesites,splicetypes[i],splicefrags_ref));
#endif

  debug1(printf("**Searching trie with max_mismatches 0, maxhits 1\n"));
  nhits0 = count_left(&bestj0,/*exclude_i*/-1,triestart,
		      splicefrags_ref,splicefrags_alt,queryptr,
		      /*max_mismatches*/0,/*maxhits*/1,query_compress,
		      dibasep,cmetp,pos5,pos3,/*nmismatches*/0,/*charpos*/pos3);
  debug1(printf("**Result is nhits0 = %d, bestj0 = %d\n",nhits0,bestj0));

  if (nhits0 > 1) {
    /* nhits0 > 1 => amb */
    *ambp = true;
    *best_nmismatches = 0;
    FREE(triebranch);
    return true;

  } else if (nhits0 == 1) {
    ext_nmismatches = Genome_count_mismatches_substring(&ncolordiffs,query,query_compress,genome_blocks,snp_blocks,
							origleft,pos5,pos3,dibasep,cmetp,plusp);
    if (ext_nmismatches == 0) {
      /* nhits0 == 1, but extension also has 0 mismatches => amb */
      *ambp = true;
      *best_nmismatches = 0;
      FREE(triebranch);
      return true;

    } else {
      /* nhits0 == 1, and extension has > 0 mismatches => unique */
      *ambp = false;
      *best_nmismatches = 0;
      *bestj = bestj0;
      FREE(triebranch);
      return true;
    }

  } else {
    /* nhits0 == 0.  Check extension for 0 mismatches */
    ext_nmismatches = Genome_count_mismatches_substring(&ncolordiffs,query,query_compress,genome_blocks,snp_blocks,
							origleft,pos5,pos3,dibasep,cmetp,plusp);
    if (ext_nmismatches == 0) {
      /* nhits0 == 0, but extension has 0 mismatches => false */
      *bestj = -1;
      FREE(triebranch);
      return false;

    } else {
      debug1(printf("**Searching trie with max_mismatches 1, maxhits 1\n"));
      nhits1 = count_left(&bestj1,/*exclude_i*/-1,triestart,
			  splicefrags_ref,splicefrags_alt,queryptr,
			  /*max_mismatches*/1,/*maxhits*/1,query_compress,
			  dibasep,cmetp,pos5,pos3,/*nmismatches*/0,/*charpos*/pos3);
      debug1(printf("**Result is nhits1 = %d, bestj1 = %d\n",nhits1,bestj1));

      if (nhits1 > 1) {
	/* nhits0 == 0, nhits1 > 1 => amb */
	*ambp = true;
	*best_nmismatches = 1;
	FREE(triebranch);
	return true;
	
      } else if (nhits1 == 1) {
	if (ext_nmismatches == 1) {
	  /* nhits0 == 0, nhits1 == 1, but extension also has 1 mismatch => amb */
	  *ambp = true;
	  *best_nmismatches = 1;
	  FREE(triebranch);
	  return true;
	  
	} else {
	  /* nhits0 == 0, nhits1 == 1, and extension has > 1 mismatch => unique */
	  *ambp = false;
	  *best_nmismatches = 1;
	  *bestj = bestj1;
	  FREE(triebranch);
	  return true;
	}

      } else {
	/* nhits0 == 0, nhits1 == 0 => fail */
	*bestj = -1;
	FREE(triebranch);
	return false;
      }
    }
  }
}

#endif

#ifdef OLD_CODE

bool
Splicetrie_find_right_short (bool *ambp, int *best_nmismatches, int *bestj,
			     int i, Genomicpos_T origleft, int pos5, int pos3, 
			     Genomicpos_T *splicesites, UINT4 *splicefrags_ref, UINT4 *splicefrags_alt,
			     unsigned int *trieoffsets, unsigned int *triecontents,
			     Splicetype_T *splicetypes, int *nsplicepartners, List_T *splicestrings,
			     Compress_T query_compress, UINT4 *genome_blocks, UINT4 *snp_blocks,
			     char *query, char *queryptr, bool dibasep, bool cmetp, bool plusp) {
  int ext_nmismatches, ncolordiffs;
  int nhits0, nhits1;
  int bestj0, bestj1;
  unsigned int *triebranch = NULL, *triestart;
  
  debug1(printf("Queryptr is %s from %d to %d: %.*s\n",queryptr,pos5,pos3,pos3-pos5,&(queryptr[pos5])));

  assert(pos5 != pos3);

  if (triecontents) {
    triestart = &(triecontents[trieoffsets[i]]);
  } else {
    triebranch = Splicetrie_build_one(&triestart,i,splicetypes,nsplicepartners,splicestrings);
  }

#if 0
  debug2(Splicetrie_dump(triestart,splicesites,splicetypes[i],splicefrags_ref));
#endif

  debug1(printf("**Searching trie with max_mismatches 0, maxhits 1\n"));
  nhits0 = count_right(&bestj0,/*exclude_i*/-1,triestart,
		       splicefrags_ref,splicefrags_alt,queryptr,
		       /*max_mismatches*/0,/*maxhits*/1,query_compress,
		       dibasep,cmetp,pos5,pos3,/*nmismatches*/0,/*charpos*/pos5-1);
  debug1(printf("**Result is nhits0 = %d, bestj0 = %d\n",nhits0,bestj0));

  if (nhits0 > 1) {
    /* nhits0 > 1 => amb */
    *ambp = true;
    *best_nmismatches = 0;
    FREE(triebranch);
    return true;

  } else if (nhits0 == 1) {
    ext_nmismatches = Genome_count_mismatches_substring(&ncolordiffs,query,query_compress,genome_blocks,snp_blocks,
							origleft,pos5,pos3,dibasep,cmetp,plusp);
    if (ext_nmismatches == 0) {
      /* nhits0 == 1, but extension also has 0 mismatches => amb */
      *ambp = true;
      *best_nmismatches = 0;
      FREE(triebranch);
      return true;

    } else {
      /* nhits0 == 1, and extension also has 0 mismatches => unique */
      *ambp = false;
      *best_nmismatches = 0;
      *bestj = bestj0;
      FREE(triebranch);
      return true;
    }

  } else {
    /* nhits0 == 0.  Check extension for 0 mismatches */
    ext_nmismatches = Genome_count_mismatches_substring(&ncolordiffs,query,query_compress,genome_blocks,snp_blocks,
							origleft,pos5,pos3,dibasep,cmetp,plusp);
    if (ext_nmismatches == 0) {
      /* nhits0 == 0, but extension has 0 mismatches => false */
      *bestj = -1;
      FREE(triebranch);
      return false;

    } else {
      debug1(printf("**Searching trie with max_mismatches 1, maxhits 1\n"));
      nhits1 = count_right(&bestj1,/*exclude_i*/-1,triestart,
			   splicefrags_ref,splicefrags_alt,queryptr,
			   /*max_mismatches*/1,/*maxhits*/1,query_compress,
			   dibasep,cmetp,pos5,pos3,/*nmismatches*/0,/*charpos*/pos5-1);
      debug1(printf("**Result is nhits1 = %d, bestj1 = %d\n",nhits1,bestj1));

      if (nhits1 > 1) {
	/* nhits0 == 0, nhits1 > 1 => amb */
	*ambp = true;
	*best_nmismatches = 1;
	FREE(triebranch);
	return true;

      } else if (nhits1 == 1) {
	if (ext_nmismatches == 1) {
	  /* nhits0 == 0, nhits1 == 1, but extension also has 1 mismatch => amb */
	  *ambp = true;
	  *best_nmismatches = 1;
	  FREE(triebranch);
	  return true;

	} else {
	  /* nhits0 == 0, nhits1 == 1, and extension has > 1 mismatch => unique */
	  *ambp = false;
	  *best_nmismatches = 1;
	  *bestj = bestj1;
	  FREE(triebranch);
	  return true;
	}

      } else {
	/* nhits0 == 0, nhits1 == 0 => fail */
	*bestj = -1;
	FREE(triebranch);
	return false;
      }
    }
  }
}

#endif


Intlist_T
Splicetrie_find_left (int *best_nmismatches, Intlist_T *nmismatches_list, int i,
		      Genomicpos_T origleft, int pos5, int pos3,
		      Genomicpos_T *splicesites, UINT4 *splicefrags_ref, UINT4 *splicefrags_alt, int *nsplicepartners_skip,
		      unsigned int *trieoffsets_obs, unsigned int *triecontents_obs, int *nsplicepartners_obs,
		      unsigned int *trieoffsets_max, unsigned int *triecontents_max, int *nsplicepartners_max,
		      Splicetype_T *splicetypes, List_T *splicestrings,
		      Compress_T query_compress, UINT4 *genome_blocks, UINT4 *snp_blocks,
		      char *query, char *queryptr, int max_mismatches_allowed,
		      bool dibasep, bool cmetp, bool plusp) {
  Intlist_T splicesites_i;
  int ncolordiffs;
  unsigned int *triebranch_obs = NULL, *triebranch_max = NULL, *triestart_obs, *triestart_max;
  int best_nmismatches_1, best_nmismatches_2;

  debug1(printf("Splicetrie_find_left called with #%d at origleft %u, pos5 %d, pos3 %d\n",i,origleft,pos5,pos3));

  if (pos5 == pos3) {
    debug1(printf("Splicetrie_find_left returning NULL, because pos5 == pos3\n"));
    return (Intlist_T) NULL;
  }


  /* Look at extension.  Cannot use splicefrags. */
  best_nmismatches_1 = 
    Genome_count_mismatches_substring(&ncolordiffs,query,query_compress,genome_blocks,snp_blocks,
				      origleft,pos5,pos3,dibasep,cmetp,plusp);
  debug1(printf("  extension at %u => %d nmismatches",origleft,best_nmismatches_1));
  if (best_nmismatches_1 > max_mismatches_allowed) {
    best_nmismatches_1 = max_mismatches_allowed;
  }
  if (best_nmismatches_1 > (pos3 - pos5)/3) {
    best_nmismatches_1 = (pos3 - pos5)/3;
    debug1(printf(" => limited by %d mismatches for length of %d\n",
		  best_nmismatches_1,pos3-pos5));
  }
  debug1(printf("\n"));

  if (triecontents_max == NULL) {
    Splicetrie_build_one(&triebranch_obs,&triestart_obs,
			 &triebranch_max,&triestart_max,
			 nsplicepartners_skip,nsplicepartners_obs,nsplicepartners_max,
			 i,splicetypes,splicestrings);
  } else if (trieoffsets_obs == NULL) {
    triestart_obs = (unsigned int *) NULL;
    triestart_max = &(triecontents_max[trieoffsets_max[i]]);
  } else {
    triestart_obs = &(triecontents_obs[trieoffsets_obs[i]]);
    triestart_max = &(triecontents_max[trieoffsets_max[i]]);
  }

  *nmismatches_list = NULL;
  if (triestart_obs == NULL) {
    debug2(printf("Trie for maximum allowed splice lengths:\n"));
    debug2(Splicetrie_dump(triestart_max,splicesites,splicetypes[i],splicefrags_ref));

    splicesites_i = search_left(&best_nmismatches_1,&(*nmismatches_list),/*splicesites_i*/NULL,
				triestart_max,splicesites,splicefrags_ref,splicefrags_alt,
				query,queryptr,query_compress,genome_blocks,snp_blocks,
				dibasep,cmetp,pos5,pos3,plusp,
				/*nmismatches*/0,/*charpos*/pos3);
    best_nmismatches_2 = best_nmismatches_1;

  } else {
    debug2(printf("Trie for observed splice lengths:\n"));
    debug2(Splicetrie_dump(triestart_obs,splicesites,splicetypes[i],splicefrags_ref));

    debug2(printf("Trie for maximum allowed splice lengths:\n"));
    debug2(Splicetrie_dump(triestart_max,splicesites,splicetypes[i],splicefrags_ref));

    splicesites_i = search_left(&best_nmismatches_1,&(*nmismatches_list),/*splicesites_i*/NULL,
				triestart_obs,splicesites,splicefrags_ref,splicefrags_alt,
				query,queryptr,query_compress,genome_blocks,snp_blocks,
				dibasep,cmetp,pos5,pos3,plusp,
				/*nmismatches*/0,/*charpos*/pos3);

    if ((best_nmismatches_2 = best_nmismatches_1 - 1) < 0) {
      best_nmismatches_2 = best_nmismatches_1;
    } else {
      splicesites_i = search_left(&best_nmismatches_2,&(*nmismatches_list),splicesites_i,
				  triestart_max,splicesites,splicefrags_ref,splicefrags_alt,
				  query,queryptr,query_compress,genome_blocks,snp_blocks,
				  dibasep,cmetp,pos5,pos3,plusp,
				  /*nmismatches*/0,/*charpos*/pos3);
    }
  }

  if (triecontents_max == NULL) {
    FREE(triebranch_max);
    FREE(triebranch_obs);
  }

  if (best_nmismatches_1 <= best_nmismatches_2) {
    *best_nmismatches = best_nmismatches_1;
  } else {
    *best_nmismatches = best_nmismatches_2;
  }
  assert(*best_nmismatches >= 0);

  debug1(printf("Splicetrie_find_left returning %s\n",Intlist_to_string(splicesites_i)));

  return splicesites_i;
}


Intlist_T
Splicetrie_find_right (int *best_nmismatches, Intlist_T *nmismatches_list, int i,
		       Genomicpos_T origleft, int pos5, int pos3,	
		       Genomicpos_T *splicesites, UINT4 *splicefrags_ref, UINT4 *splicefrags_alt, int *nsplicepartners_skip,
		       unsigned int *trieoffsets_obs, unsigned int *triecontents_obs, int *nsplicepartners_obs,
		       unsigned int *trieoffsets_max, unsigned int *triecontents_max, int *nsplicepartners_max,
		       Splicetype_T *splicetypes, List_T *splicestrings,
		       Compress_T query_compress, UINT4 *genome_blocks, UINT4 *snp_blocks,
		       char *query, char *queryptr, int max_mismatches_allowed,
		       bool dibasep, bool cmetp, bool plusp) {
  Intlist_T splicesites_i;
  int ncolordiffs;
  unsigned int *triebranch_obs = NULL, *triebranch_max = NULL, *triestart_obs, *triestart_max;
  int best_nmismatches_1, best_nmismatches_2;

  debug1(printf("Splicetrie_find_right called with #%d at origleft %u, pos5 %d, pos3 %d\n",i,origleft,pos5,pos3));

  if (pos5 == pos3) {
    debug1(printf("Splicetrie_find_right returning NULL, because pos5 == pos3\n"));
    return (Intlist_T) NULL;
  }

  /* Look at extension.  Cannot use splicefrags. */
  best_nmismatches_1 = 
    Genome_count_mismatches_substring(&ncolordiffs,query,query_compress,genome_blocks,snp_blocks,
				      origleft,pos5,pos3,dibasep,cmetp,plusp);
  debug1(printf("  extension at %u => %d nmismatches",origleft,best_nmismatches_1));
  if (best_nmismatches_1 > max_mismatches_allowed) {
    best_nmismatches_1 = max_mismatches_allowed;
  }
  if (best_nmismatches_1 > (pos3 - pos5)/3) {
    best_nmismatches_1 = (pos3 - pos5)/3;
    debug1(printf(" => limited by %d mismatches for length of %d\n",
		  best_nmismatches_1,pos3-pos5));
  }
  debug1(printf("\n"));


  if (triecontents_max == NULL) {
    Splicetrie_build_one(&triebranch_obs,&triestart_obs,
			 &triebranch_max,&triestart_max,
			 nsplicepartners_skip,nsplicepartners_obs,nsplicepartners_max,
			 i,splicetypes,splicestrings);
  } else if (trieoffsets_obs == NULL) {
    triestart_obs = (unsigned int *) NULL;
    triestart_max = &(triecontents_max[trieoffsets_max[i]]);
  } else {
    triestart_obs = &(triecontents_obs[trieoffsets_obs[i]]);
    triestart_max = &(triecontents_max[trieoffsets_max[i]]);
  }

  *nmismatches_list = NULL;
  if (triestart_obs == NULL) {
    debug2(printf("Trie for maximum allowed splice lengths:\n"));
    debug2(Splicetrie_dump(triestart_max,splicesites,splicetypes[i],splicefrags_ref));

    splicesites_i = search_right(&best_nmismatches_1,&(*nmismatches_list),/*splicesites_i*/NULL,
				 triestart_max,splicesites,splicefrags_ref,splicefrags_alt,
				 query,queryptr,query_compress,genome_blocks,snp_blocks,
				 dibasep,cmetp,pos5,pos3,plusp,
				 /*nmismatches*/0,/*charpos*/pos5-1);
    best_nmismatches_2 = best_nmismatches_1;

  } else {
    debug2(printf("Trie for observed splice lengths:\n"));
    debug2(Splicetrie_dump(triestart_obs,splicesites,splicetypes[i],splicefrags_ref));

    debug2(printf("Trie for maximum allowed splice lengths:\n"));
    debug2(Splicetrie_dump(triestart_max,splicesites,splicetypes[i],splicefrags_ref));
    
    splicesites_i = search_right(&best_nmismatches_1,&(*nmismatches_list),/*splicesites_i*/NULL,
				 triestart_obs,splicesites,splicefrags_ref,splicefrags_alt,
				 query,queryptr,query_compress,genome_blocks,snp_blocks,
				 dibasep,cmetp,pos5,pos3,plusp,
				 /*nmismatches*/0,/*charpos*/pos5-1);

    if ((best_nmismatches_2 = best_nmismatches_1 - 1) < 0) {
      best_nmismatches_2 = best_nmismatches_1;
    } else {
      splicesites_i = search_right(&best_nmismatches_2,&(*nmismatches_list),splicesites_i,
				   triestart_max,splicesites,splicefrags_ref,splicefrags_alt,
				   query,queryptr,query_compress,genome_blocks,snp_blocks,
				   dibasep,cmetp,pos5,pos3,plusp,
				   /*nmismatches*/0,/*charpos*/pos5-1);
    }
  }

  if (triecontents_max == NULL) {
    FREE(triebranch_max);
    FREE(triebranch_obs);
  }

  if (best_nmismatches_1 <= best_nmismatches_2) {
    *best_nmismatches = best_nmismatches_1;
  } else {
    *best_nmismatches = best_nmismatches_2;
  }
  assert(*best_nmismatches >= 0);

  debug1(printf("Splicetrie_find_right returning %s\n",Intlist_to_string(splicesites_i)));

  return splicesites_i;
}


#endif

