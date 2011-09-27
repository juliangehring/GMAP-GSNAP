static char rcsid[] = "$Id: snpindex.c,v 1.11 2010-07-27 01:10:20 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifndef HAVE_MEMCPY
# define memcpy(d,s,n) bcopy((s),(d),(n))
#endif
#ifndef HAVE_MEMMOVE
# define memmove(d,s,n) bcopy((s),(d),(n))
#endif

#ifdef WORDS_BIGENDIAN
#include "bigendian.h"
#else
#include "littleendian.h"
#endif

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>		/* For memset */
#include <ctype.h>		/* For toupper */
#include <sys/mman.h>		/* For munmap */
#include <math.h>		/* For qsort */
#ifdef HAVE_UNISTD_H
#include <unistd.h>		/* For lseek and close */
#endif
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>		/* For off_t */
#endif
#if HAVE_DIRENT_H
# include <dirent.h>
# define NAMLEN(dirent) strlen((dirent)->d_name)
#else
# define dirent direct
# define NAMLEN(dirent) (dirent)->d_namlen
# if HAVE_SYS_NDIR_H
#  include <sys/ndir.h>
# endif
# if HAVE_SYS_DIR_H
#  include <sys/dir.h>
# endif
# if HAVE_NDIR_H
#  include <ndir.h>
# endif
#endif

#include "mem.h"
#include "fopen.h"
#include "access.h"
#include "types.h"

#include "compress.h"
#include "interval.h"
#include "complement.h"

#include "bool.h"
#include "genomicpos.h"
#include "iitdef.h"
#include "uintlist.h"
#include "chrnum.h"
#include "genome.h"
#include "datadir.h"
#include "iit-read.h"
#include "indexdb.h"
#include "indexdbdef.h"


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* SNP blocks and writing of positions */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

#ifdef CHECK
#define check(x) x
#else
#define check(x)
#endif



static char *genomesubdir = ".";
static char *dbroot = NULL;
static char *fileroot = NULL;

static char *destdir = ".";
static char *dbversion = NULL;
static char *user_genomedir = NULL;

static char *snps_root = NULL;
static bool show_warnings_p = true;


static int
power (int base, int exponent) {
  int result = 1, i;

  for (i = 0; i < exponent; i++) {
    result *= base;
  }
  return result;
}


/*                87654321 */
#define RIGHT_A 0x00000000
#define RIGHT_C 0x00000001
#define RIGHT_G 0x00000002
#define RIGHT_T 0x00000003

/*                      87654321 */
#define LOW_TWO_BITS  0x00000003

static char *
shortoligo_nt (Storedoligomer_T oligo, int oligosize) {
  char *nt;
  int i, j;
  Storedoligomer_T lowbits;

  nt = (char *) CALLOC(oligosize+1,sizeof(char));
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

  return nt;
}

static char
check_acgt (char nt, IIT_T snps_iit, int divno, int i) {
  int index;
  char *label;
  bool allocp;

  if (nt == 'N') {
    return nt;
  } else if (nt != 'A' && nt != 'C' && nt != 'G' && nt != 'T') {
    index = IIT_index(snps_iit,divno,i);
    label = IIT_label(snps_iit,index,&allocp);
    fprintf(stderr,"\nFor %s, alternate allele %c is not ACGT, using N as alternate allele",
	    label,nt);
    if (allocp) {
      FREE(label);
    }
    return 'N';
  } else {
    return nt;
  }
}


static Storedoligomer_T
revise_oligomer (Storedoligomer_T oligo, char nt) {
  switch (nt) {
  case 'A': return (oligo << 2);
  case 'C': return (oligo << 2) | 1U;
  case 'G': return (oligo << 2) | 2U;
  case 'T': return (oligo << 2) | 3U;
  }
  fprintf(stderr,"altstring char is %c\n",nt);
  exit(9);
  return 0;
}


#if 0
static bool
samep (char *string1, char *string2, int length) {
  int i;

  for (i = 0; i < length; i++) {
    if (string1[i] != string2[i]) {
      return false;
    }
  }
  return true;
}
#endif

static int
process_snp_block (Positionsptr_T *offsets, Genomicpos_T *positions,
		   Interval_T *intervals, int nintervals,
		   Genomicpos_T chroffset, Genome_T genome,
		   UINT4 *genome_blocks, int divno, int intervali,
		   IIT_T snps_iit, IIT_T chromosome_iit) {
  int nerrors = 0;
  bool *snpp;
  char *refstring;
  char *altstring;

  int index, length;
  int nsnps, stringi, starti, shift, i, k;
  char *snptype, *label, refnt, altnt;
  unsigned int ptr;
  Genomicpos_T snpposition, startposition, endposition, first_snppos, last_snppos, position, chrpos;
  Chrnum_T chrnum;
  Interval_T interval;
  int nunknowns;
  bool badcharp, allocp;
  
  Uintlist_T oligomers, newoligomers, p;
  Storedoligomer_T oligo;

  /* Subtract 1 because snps_iit is 1-based */
  first_snppos = Interval_low(intervals[0]) - 1U;
  last_snppos = Interval_low(intervals[nintervals-1]) - 1U;
  debug1(
	 if (nintervals == 1) {
	   printf("Processing snp at chrpos %u\n",first_snppos);
	 } else {
	   printf("Processing block of %d snps from chrpos %u to %u\n",
		  nintervals,first_snppos,last_snppos);
	 }
	 );

  if (first_snppos < (INDEX1PART - 1)) {
    startposition = chroffset;
  } else {
    startposition = chroffset + first_snppos - (INDEX1PART - 1);
  }
  endposition = chroffset + last_snppos + (INDEX1PART - 1);
  length = endposition - startposition + 1;

  snpp = (bool *) CALLOC(length,sizeof(bool));
  refstring = (char *) CALLOC(length + 1,sizeof(char));
  altstring = (char *) CALLOC(length,sizeof(char));

  Genome_fill_buffer(&chrnum,&nunknowns,genome,/*left*/startposition,length,
		     refstring,chromosome_iit);

  for (i = 0; i < nintervals; i++) {
    interval = intervals[i];
    snpposition = chroffset + Interval_low(interval) - 1U; /* Subtract 1 because snps_iit is 1-based */
    debug1(printf("  Neighbor at %u\n",snpposition));
    stringi = snpposition - startposition;

    snptype = IIT_typestring(snps_iit,Interval_type(interval));
    if (strlen(snptype) != 2) {
      fprintf(stderr,"Unrecognized snptype %s\n",snptype);
      abort();
    } else {
      refnt = refstring[stringi];
      if (refnt == snptype[0]) {
	altnt = check_acgt(snptype[1],snps_iit,divno,intervali);
      } else if (refnt == snptype[1]) {
	altnt = check_acgt(snptype[0],snps_iit,divno,intervali);
      } else {
	nerrors++;
	if (show_warnings_p == true) {
	  index = IIT_index(snps_iit,divno,intervali);
	  label = IIT_label(snps_iit,index,&allocp);
	  fprintf(stderr,"\nFor %s at %u, snptype %s not consistent with reference allele %c, using N as alternate allele",
		  label,snpposition,snptype,refstring[stringi]);
	  if (allocp == true) {
	    FREE(label);
	  }
	}
	altnt = 'N';
      }

      if (altnt != 'A' && altnt != 'C' && altnt != 'G' && altnt != 'T') {
	/* fprintf(stderr,"\nAt %u, alternate allele %c is not ACGT\n",snpposition,altnt); */
      } else if (altnt == refnt) {
	fprintf(stderr,"\nAt %u, alternate allele %c is same as reference allele\n",snpposition,altnt);
      } else {
	altstring[stringi] = altnt;
	snpp[stringi] = true;
	if (genome_blocks != NULL) {
	  /* revising genome */
	  ptr = snpposition/32U*3;
	  shift = snpposition % 32U;

	  genome_blocks[ptr+2] |= (1 << shift);	/* flags */
	  if (shift >= 16) {
	    /* high */
	    shift -= 16;
	    genome_blocks[ptr] &= ~(0x3 << (2*shift)); /* clear bits */
	    switch (altnt) {	/* set bits */
	    case 'A': break;
	    case 'C': genome_blocks[ptr] |= (0x1 << (2*shift)); break;
	    case 'G': genome_blocks[ptr] |= (0x2 << (2*shift)); break;
	    case 'T': genome_blocks[ptr] |= (0x3 << (2*shift)); break;
	    default: abort();
	    }
	  } else {
	    /* low */
	    genome_blocks[ptr+1] &= ~(0x3 << (2*shift)); /* clear bits */
	    switch (altnt) {	/* set bits */
	    case 'A': break;
	    case 'C': genome_blocks[ptr+1] |= (0x1 << (2*shift)); break;
	    case 'G': genome_blocks[ptr+1] |= (0x2 << (2*shift)); break;
	    case 'T': genome_blocks[ptr+1] |= (0x3 << (2*shift)); break;
	    default: abort();
	    }
	  }
	}
      }
    }
  }
  
  for (starti = 0, position = startposition, chrpos = startposition - chroffset; 
       position <= endposition-INDEX1PART+1U; starti++, position++, chrpos++) {
    if (chrpos % 3 == 0) {
      /* chrpos % 3 == 0 is same as the condition in indexdb.c for storing a position */
      nsnps = 0;
      badcharp = false;
      for (k = starti; k < starti + INDEX1PART; k++) {
	if (snpp[k] == true) {
	  nsnps++;
	}
	if (refstring[k] != 'A' && refstring[k] != 'C' && refstring[k] != 'G' && refstring[k] != 'T') {
	  badcharp = true;
	}
      }
      if (nsnps == 0) {
	/* no snps */
	/* fprintf(stderr,"\nNo snps at position %u",position); */
      } else if (nsnps > 4) {
	/* too many snps */
      } else if (badcharp == true) {
	/* bad reference char */
      } else {
	oligomers = Uintlist_push(NULL,0U);
	for (k = starti, i = 0; k < starti + INDEX1PART; k++, i++) {
	  newoligomers = NULL;
	  if (snpp[k] == false) {
	    for (p = oligomers; p != NULL; p = Uintlist_next(p)) {
	      newoligomers = Uintlist_push(newoligomers,revise_oligomer(Uintlist_head(p),refstring[k]));
	    }
	    
	  } else if (altstring[k] != 'N') {
	    for (p = oligomers; p != NULL; p = Uintlist_next(p)) {
	      newoligomers = Uintlist_push(newoligomers,revise_oligomer(Uintlist_head(p),refstring[k]));
	      newoligomers = Uintlist_push(newoligomers,revise_oligomer(Uintlist_head(p),altstring[k]));
	    }

	  } else {
	    for (p = oligomers; p != NULL; p = Uintlist_next(p)) {
	      newoligomers = Uintlist_push(newoligomers,revise_oligomer(Uintlist_head(p),refstring[k]));
	      if (refstring[k] != 'A') {
		newoligomers = Uintlist_push(newoligomers,revise_oligomer(Uintlist_head(p),'A'));
	      }
	      if (refstring[k] != 'C') {
		newoligomers = Uintlist_push(newoligomers,revise_oligomer(Uintlist_head(p),'C'));
	      }
	      if (refstring[k] != 'G') {
		newoligomers = Uintlist_push(newoligomers,revise_oligomer(Uintlist_head(p),'G'));
	      }
	      if (refstring[k] != 'T') {
		newoligomers = Uintlist_push(newoligomers,revise_oligomer(Uintlist_head(p),'T'));
	      }
	    }
	  }

	  Uintlist_free(&oligomers);
	  oligomers = Uintlist_reverse(newoligomers);
	}

#ifdef CHECK
	for (p = Uintlist_next(oligomers); p != NULL; p = Uintlist_next(p)) {
	  oligo = Uintlist_head(p);
	  nt = shortoligo_nt(oligo,INDEX1PART);
	  if (samep(nt,&(refstring[starti]),INDEX1PART) == true) {
	    fprintf(stderr,"Storing oligomer %s that is the same as the reference at %u\n",
		    nt,position);
	    abort();
	  }
	  FREE(nt);
	}
#endif	

	/* Ignore the first element in oligomers, which is all reference */
	if (positions == NULL) {
	  /* writing offsets */
	  for (p = Uintlist_next(oligomers); p != NULL; p = Uintlist_next(p)) {
	    oligo = Uintlist_head(p);
	    offsets[oligo + 1U] += 1;
	    debug1(nt = shortoligo_nt(oligo,INDEX1PART);
		   printf("Storing %s at %u, chrpos %u\n",nt,position,chrpos);
		   FREE(nt));
	  }
	} else {
	  /* writing positions */
	  for (p = Uintlist_next(oligomers); p != NULL; p = Uintlist_next(p)) {
	    oligo = Uintlist_head(p);
	    positions[offsets[oligo]++] = position;
	    debug1(nt = shortoligo_nt(oligo,INDEX1PART);
		   printf("Storing %s at %u, chrpos %u\n",nt,position,chrpos);
		   FREE(nt));
	  }
	}
	Uintlist_free(&oligomers);
      }
    }
  }

  FREE(altstring);
  FREE(refstring);
  FREE(snpp);

  return nerrors;
}


static Positionsptr_T *
compute_offsets (IIT_T snps_iit, IIT_T chromosome_iit, Genome_T genome, UINT4 *genome_blocks,
		 int oligospace) {
  Positionsptr_T *offsets;

  Interval_T *intervals;
  int nintervals, nerrors, divno, i, j;
  char *divstring;
  Chrnum_T chrnum;
  Genomicpos_T chroffset;


  offsets = (Positionsptr_T *) CALLOC(oligospace+1,sizeof(Positionsptr_T));

  for (divno = 1; divno < snps_iit->ndivs; divno++) {
    divstring = IIT_divstring(snps_iit,divno);
    fprintf(stderr,"Processing offsets for chromosome %s...",divstring);
    if ((chrnum = IIT_find_one(chromosome_iit,divstring)) <= 0) {
      fprintf(stderr,"not found in chromosome iit\n");
    } else {
      nerrors = 0;
      fprintf(stderr,"has %d snps...",IIT_nintervals(snps_iit,divno));
      chroffset = IIT_interval_low(chromosome_iit,chrnum);

      nintervals = IIT_nintervals(snps_iit,divno);
      intervals = (Interval_T *) CALLOC(nintervals,sizeof(Interval_T));
      for (i = 0; i < nintervals; i++) {
	intervals[i] = &(snps_iit->intervals[divno][i]);
      }
      qsort(intervals,nintervals,sizeof(Interval_T),Interval_cmp);

      i = 0;
      while (i < nintervals) {
	j = i + 1;
	while (j < nintervals && Interval_low(intervals[j]) < Interval_low(intervals[j-1]) + INDEX1PART) {
	  j++;
	}
	nerrors += process_snp_block(offsets,/*positions*/NULL,&(intervals[i]),j-i,
				     chroffset,genome,genome_blocks,
				     divno,/*intervali*/i,snps_iit,chromosome_iit);
	i = j;
      }

      FREE(intervals);
      fprintf(stderr,"done (%d snps inconsistent with reference genome)\n",nerrors);
    }
  }

  /* Do not add extra for sentinel */

  for (i = 1; i <= oligospace; i++) {
    offsets[i] = offsets[i] + offsets[i-1];
    debug(if (offsets[i] != offsets[i-1]) {
	    printf("Offset for %06X: %u\n",i,offsets[i]);
	  });
  }

  return offsets;
}


static Genomicpos_T *
compute_positions (Positionsptr_T *offsets, IIT_T snps_iit, IIT_T chromosome_iit,
		   Genome_T genome, int oligospace) {
  Genomicpos_T *positions;

  Interval_T *intervals;
  int nintervals, divno, i, j;
  char *divstring;
  Chrnum_T chrnum;
  Genomicpos_T chroffset;
  Positionsptr_T *pointers, totalcounts, block_start, block_end, npositions;


  totalcounts = offsets[oligospace];
  if (totalcounts == 0) {
    fprintf(stderr,"Something is wrong with the offsets.  Total counts is zero.\n");
    exit(9);
  } else {
    fprintf(stderr,"Trying to allocate %u*%d bytes of memory...",totalcounts,(int) sizeof(Genomicpos_T));
    positions = (Genomicpos_T *) CALLOC_NO_EXCEPTION(totalcounts,sizeof(Genomicpos_T));
    if (positions == NULL) {
      fprintf(stderr,"failed.  Need a computer with sufficient memory.\n");
      exit(9);
    } else {
      fprintf(stderr,"done\n");
    }
  }

  /* Copy offsets */
  pointers = (Positionsptr_T *) CALLOC(oligospace+1,sizeof(Positionsptr_T));
  for (i = 0; i <= oligospace; i++) {
    pointers[i] = offsets[i];
  }

  for (divno = 1; divno < snps_iit->ndivs; divno++) {
    divstring = IIT_divstring(snps_iit,divno);
    fprintf(stderr,"Processing positions for chromosome %s...",divstring);
    if ((chrnum = IIT_find_one(chromosome_iit,divstring)) <= 0) {
      fprintf(stderr,"not found in chromosome iit\n");
    } else {
      fprintf(stderr,"has %d snps...",IIT_nintervals(snps_iit,divno));
      chroffset = IIT_interval_low(chromosome_iit,chrnum);

      nintervals = IIT_nintervals(snps_iit,divno);
      intervals = (Interval_T *) CALLOC(nintervals,sizeof(Interval_T));
      for (i = 0; i < nintervals; i++) {
	intervals[i] = &(snps_iit->intervals[divno][i]);
      }
      qsort(intervals,nintervals,sizeof(Interval_T),Interval_cmp);

      i = 0;
      while (i < nintervals) {
	j = i + 1;
	while (j < nintervals && Interval_low(intervals[j]) < Interval_low(intervals[j-1]) + INDEX1PART) {
	  j++;
	}
	process_snp_block(pointers,positions,&(intervals[i]),j-i,
			  chroffset,genome,/*genome_blocks*/NULL,
			  divno,/*intervali*/i,snps_iit,chromosome_iit);
	i = j;
      }

      FREE(intervals);
      fprintf(stderr,"done\n");
    }
  }

  FREE(pointers);

  /* Sort positions in each block */
  for (i = 0; i < oligospace; i++) {
    block_start = offsets[i];
    block_end = offsets[i+1];
    if ((npositions = block_end - block_start) > 1) {
      qsort(&(positions[block_start]),npositions,sizeof(Genomicpos_T),Genomicpos_compare);
    }
  }

  return positions;
}


static void
merge_positions (FILE *positions_fp, Genomicpos_T *start1, Genomicpos_T *end1,
		 Genomicpos_T *start2, Genomicpos_T *end2, Positionsptr_T oligo) {
  Genomicpos_T *ptr1 = start1, *ptr2 = start2;
  char *nt;

  while (ptr1 < end1 && ptr2 < end2) {
    if (*ptr1 < *ptr2) {
      FWRITE_UINT(*ptr1,positions_fp);
      ptr1++;
    } else if (*ptr2 < *ptr1) {
      FWRITE_UINT(*ptr2,positions_fp);
      ptr2++;
    } else {
      nt = shortoligo_nt(oligo,INDEX1PART);
      fprintf(stderr,"Problem: saw duplicate positions %u in oligo %s\n",*ptr1,nt);
      FREE(nt);
      abort();
      /*
      FWRITE_UINT(*ptr1,positions_fp);
      ptr1++;
      ptr2++;
      */
    }
  }

  while (ptr1 < end1) {
    FWRITE_UINT(*ptr1,positions_fp);
    ptr1++;
  }

  while (ptr2 < end2) {
    FWRITE_UINT(*ptr2,positions_fp);
    ptr2++;
  }

  return;
}


/* Usage: snpindex -d <genome> -V <snps_root> <iitfile> */


/* Note: Depends on having gmapindex sampled on mod 3 bounds */
int
main (int argc, char *argv[]) {
  char *iitfile;
  IIT_T chromosome_iit, snps_iit;
  Genome_T genome;
  Positionsptr_T *offsets, *snp_offsets, *ref_offsets, npositions;
  Genomicpos_T *snp_positions, *ref_positions, nblocks;
  UINT4 *genome_blocks;
  int oligospace, i;

  char *filename;
  FILE *genome_fp, *offsets_fp, *positions_fp, *ref_offsets_fp, *ref_positions_fp;
  int ref_offsets_fd, ref_positions_fd;
  size_t ref_offsets_len, ref_positions_len;

  int c;
  extern int optind;
  extern char *optarg;

  while ((c = getopt(argc,argv,"D:d:V:")) != -1) {
    switch (c) {
    case 'D': user_genomedir = optarg; break;
    case 'd': dbroot = optarg; break;
    case 'V': snps_root = optarg; break;
    }
  }
  argc -= (optind - 1);
  argv += (optind - 1);

  if (dbroot == NULL) {
    fprintf(stderr,"Missing name of genome database.  Must specify with -d flag.\n");
    fprintf(stderr,"Usage: snpindex -d <genome> -V <snpversion> <snp_iitfile>\n");
    exit(9);
  } else if (snps_root == NULL) {
    fprintf(stderr,"Missing name of snp version.  Must specify with -V flag.\n");
    fprintf(stderr,"Usage: snpindex -d <genome> -V <snpversion> <snp_iitfile>\n");
    exit(9);
  }


  if (argc <= 1) {
    fprintf(stderr,"Must specify an iit file\n");
    exit(9);
  } else if (Access_file_exists_p(argv[1]) == false) {
    fprintf(stderr,"SNP IIT file %s not found\n",argv[1]);
    exit(9);
  } else {
    fprintf(stderr,"Reading SNPs IIT file %s...",argv[1]);
    if ((snps_iit = IIT_read(argv[1],/*name*/NULL,/*readonlyp*/true,
			     /*divread*/READ_ALL,/*divstring*/NULL,/*add_iit_p*/false,
			     /*labels_read_p*/true)) == NULL) {
      fprintf(stderr,"SNP IIT file %s is not valid\n",argv[1]);
      exit(9);
    }
  }


  /* Chromosome IIT file */
  iitfile = (char *) CALLOC(strlen(dbroot)+strlen(".chromosome.iit")+1,sizeof(char));
  sprintf(iitfile,"%s.chromosome.iit",dbroot);
  if (Access_file_exists_p(iitfile) == true) {
    genomesubdir = ".";
    fileroot = dbroot;
    if ((chromosome_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,
				   /*divread*/READ_ALL,/*divstring*/NULL,/*add_iit_p*/false,
				   /*labels_read_p*/true)) == NULL) {
      fprintf(stderr,"Found %s in this directory, but not valid\n",iitfile);
      exit(9);
    } else {
      FREE(iitfile);
    }

  } else {
    genomesubdir = Datadir_find_genomesubdir(&fileroot,&dbversion,user_genomedir,dbroot);
    iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			      strlen(fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
    sprintf(iitfile,"%s/%s.chromosome.iit",genomesubdir,fileroot);
    if ((chromosome_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,
				   /*divread*/READ_ALL,/*divstring*/NULL,/*add_iit_p*/false,
				   /*labels_read_p*/true)) == NULL) {
      fprintf(stderr,"Could not find valid IIT in %s\n",iitfile);
      exit(9);
    } else {
      FREE(iitfile);
    }
  }


  genome = Genome_new(genomesubdir,fileroot,/*snps_root*/NULL,/*uncompressedp*/false,
		      /*batch_genome_p*/false);

  /* Copy genome */
  nblocks = Genome_totallength(genome)/32U;
  genome_blocks = (UINT4 *) CALLOC(nblocks*3,sizeof(UINT4));
  fprintf(stderr,"Allocating %u*3*%lu bytes for compressed genome\n",nblocks,sizeof(UINT4));
  memcpy(genome_blocks,Genome_blocks(genome),nblocks*3*sizeof(UINT4));


  /* Compute offsets and write genome */
  oligospace = power(4,INDEX1PART);
  snp_offsets = compute_offsets(snps_iit,chromosome_iit,genome,genome_blocks,oligospace);
  filename = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
			     strlen(".genomecomp.")+strlen(snps_root)+1,sizeof(char));
  sprintf(filename,"%s/%s.genomecomp.%s",destdir,fileroot,snps_root);
  if ((genome_fp = FOPEN_WRITE_BINARY(filename)) == NULL) {
    fprintf(stderr,"Can't open file %s for writing genome\n",filename);
    exit(9);
  }
  FREE(filename);
  FWRITE_UINTS(genome_blocks,nblocks*3,genome_fp);
  fclose(genome_fp);


  /* Compute positions */
  snp_positions = compute_positions(snp_offsets,snps_iit,chromosome_iit,genome,oligospace);


  /* Read reference offsets and update */
  filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+strlen(fileroot)+
			     strlen(".")+strlen(IDX_FILESUFFIX)+strlen("3")+
			     strlen(OFFSETS_FILESUFFIX)+1,sizeof(char));
  sprintf(filename,"%s/%s.%s%s%s",genomesubdir,fileroot,IDX_FILESUFFIX,"3",OFFSETS_FILESUFFIX);
  if ((ref_offsets_fp = FOPEN_READ_BINARY(filename)) == NULL) {
    fprintf(stderr,"Can't open ref_offsets file %s\n",filename);
    exit(9);
  } else {
    fclose(ref_offsets_fp);
  }
  ref_offsets = (Positionsptr_T *) Access_mmap(&ref_offsets_fd,&ref_offsets_len,
					       filename,sizeof(Positionsptr_T),/*randomp*/false);
  FREE(filename);

  offsets = (Positionsptr_T *) CALLOC(oligospace+1,sizeof(Positionsptr_T));
  offsets[0] = 0U;
  for (i = 1; i <= oligospace; i++) {
    npositions = (ref_offsets[i] - ref_offsets[i-1]) + (snp_offsets[i] - snp_offsets[i-1]);
    offsets[i] = offsets[i-1] + npositions;
  }


  /* Write offsets */
  filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+strlen(fileroot)+
			     strlen(".")+strlen("ref")+strlen("3")+
			     strlen(OFFSETS_FILESUFFIX)+strlen(".")+strlen(snps_root)+1,sizeof(char));
  sprintf(filename,"%s/%s.%s%s%s.%s",genomesubdir,fileroot,"ref","3",OFFSETS_FILESUFFIX,snps_root);
  if ((offsets_fp = FOPEN_WRITE_BINARY(filename)) == NULL) {
    fprintf(stderr,"Can't open file %s for writing\n",filename);
    exit(9);
  }
  FREE(filename);

  fprintf(stderr,"Writing %d offsets to file with %u total positions\n",
	  oligospace+1,offsets[oligospace]);
  FWRITE_UINTS(offsets,oligospace+1,offsets_fp);
  FREE(offsets);


  /* Read reference positions and merge */
  filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+strlen(fileroot)+
			     strlen(".")+strlen(IDX_FILESUFFIX)+strlen("3")+
			     strlen(POSITIONS_FILESUFFIX)+1,sizeof(char));
  sprintf(filename,"%s/%s.%s%s%s",genomesubdir,fileroot,IDX_FILESUFFIX,"3",POSITIONS_FILESUFFIX);
  if ((ref_positions_fp = FOPEN_READ_BINARY(filename)) == NULL) {
    fprintf(stderr,"Can't open file %s\n",filename);
    exit(9);
  } else {
    fclose(ref_positions_fp);
  }
  ref_positions = (Genomicpos_T *) Access_mmap(&ref_positions_fd,&ref_positions_len,
					       filename,sizeof(Genomicpos_T),/*randomp*/false);
  FREE(filename);

  filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+strlen(fileroot)+
			     strlen(".")+strlen("ref")+strlen("3")+
			     strlen(POSITIONS_FILESUFFIX)+strlen(".")+strlen(snps_root)+1,sizeof(char));
  sprintf(filename,"%s/%s.%s%s%s.%s",genomesubdir,fileroot,"ref","3",POSITIONS_FILESUFFIX,snps_root);
  if ((positions_fp = FOPEN_WRITE_BINARY(filename)) == NULL) {
    fprintf(stderr,"Can't open file %s for writing\n",filename);
    exit(9);
  }
  FREE(filename);

  fprintf(stderr,"Merging snp positions with reference positions\n");
  for (i = 0; i < oligospace; i++) {
    merge_positions(positions_fp,&(snp_positions[snp_offsets[i]]),&(snp_positions[snp_offsets[i+1]]),
		    &(ref_positions[ref_offsets[i]]),&(ref_positions[ref_offsets[i+1]]),i);
  }
  fclose(positions_fp);


  /* Clean up */
  munmap((void *) ref_positions,ref_positions_len);
  close(ref_positions_fd);

  munmap((void *) ref_offsets,ref_offsets_len);
  close(ref_offsets_fd);

  FREE(snp_positions);
  FREE(snp_offsets);

  Genome_free(&genome);
  IIT_free(&chromosome_iit);
  IIT_free(&snps_iit);

  return 0;
}
