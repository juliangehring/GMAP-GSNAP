static char rcsid[] = "$Id: indexdb-write.c 127628 2014-02-18 19:45:01Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef HAVE_MEMCPY
# define memcpy(d,s,n) bcopy((s),(d),(n))
#endif
#ifndef HAVE_MEMMOVE
# define memmove(d,s,n) bcopy((s),(d),(n))
#endif

#include "indexdb-write.h"

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
#include "types.h"		/* For Oligospace_T */

#include "compress-write.h"	/* For Compress_get_char */
#include "interval.h"
#include "complement.h"
#include "access.h"
#include "genomicpos.h"
#include "bool.h"
#include "indexdbdef.h"
#include "genome_hr.h"		/* For read_gammas procedures */
#include "iit-read-univ.h"
#include "indexdb.h"
#include "popcount.h"

#include "bitpack64-read.h"
#include "bitpack64-write.h"


/* Another MONITOR_INTERVAL is in compress.c */
#define MONITOR_INTERVAL 100000000 /* 100 million nt */
#define OFFSETS_BUFFER_SIZE 1000000
#define WRITE_CHUNK 1000000
#define POSITIONS8_HIGH_SHIFT 32
#define POSITIONS8_LOW_MASK 0xFFFFFFFF

/* #define ALLOW_ODD_PACKSIZES 1 */

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Writing of positions */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

/* Gammas */
#ifdef DEBUG3
#define debug3(x) x
#else
#define debug3(x)
#endif




/************************************************************************
 *   Write procedures -- called by gmapindex/pmapindex
 ************************************************************************/

static Oligospace_T
power (int base, int exponent) {
#ifdef OLIGOSPACE_NOT_LONG
  Oligospace_T result = 1U;
#else
  Oligospace_T result = 1UL;
#endif
  int i;

  for (i = 0; i < exponent; i++) {
    result *= base;
  }
  return result;
}


/************************************************************************
 *   Elias gamma representation
 ************************************************************************/

/* ctr is 32 at high bit and 1 at low bit */
static int
write_gamma (FILE *offsetscomp_fp, Offsetscomp_T *offsets_buffer, int offsets_buffer_size, int *offsets_buffer_i,
	     unsigned int *nwritten, Offsetscomp_T *buffer, int ctr, Offsetscomp_T gamma) {
  int length;
  Positionsptr_T nn;
  
  debug(printf("Entering write_gamma with gamma %u, ctr %d\n",gamma,ctr));

  gamma += 1;			/* To allow 0 to be represented */

  /* Compute length */
  length = 1;
  nn = 2;
  while (nn <= gamma) {
    length += 2;
    nn += nn;
  }
  debug(printf("gamma is %u (%08X), length is %u\n",gamma,gamma,length));


  /* Update buffer and write */
  while (length > ctr) {
    if (length - ctr < 32) {
      *buffer |= (gamma >> (length - ctr));
    }
    debug(printf("writing gamma %08X\n",*buffer));
    offsets_buffer[(*offsets_buffer_i)++] = *buffer;
    if (*offsets_buffer_i == offsets_buffer_size) {
      FWRITE_UINTS(offsets_buffer,offsets_buffer_size,offsetscomp_fp);
      *offsets_buffer_i = 0;
    }
    *nwritten += 1;
    length -= ctr;
    ctr = 32;
    *buffer = 0U;
  }
  
  debug(printf("  shifting gamma left by %d\n",ctr - length));
  *buffer |= (gamma << (ctr - length));
  debug(printf("  buffer is %08X\n",*buffer));
  ctr -= length;

  debug(printf("  returning ctr %d\n",ctr));
  return ctr;
}


void
Indexdb_write_gammaptrs (char *gammaptrsfile, char *offsetsfile, Positionsptr_T *offsets,
			 Oligospace_T oligospace, Blocksize_T blocksize) {
  FILE *offsets_fp;
  FILE *gammaptrs_fp, *offsetscomp_fp;
  Gammaptr_T *gammaptrs;
  int gammaptri;
  Blocksize_T j;
  Oligospace_T oligoi;

  Offsetscomp_T *offsets_buffer;
  int offsets_buffer_size = OFFSETS_BUFFER_SIZE;
  int offsets_buffer_i;

  Offsetscomp_T buffer;			/* gammabuffer */
  int ctr;
  Gammaptr_T nwritten;


  if (blocksize == 1) {
    /* Don't write gammaptrs.  Write offsetscomp in a single command. */
    if ((offsets_fp = FOPEN_WRITE_BINARY(offsetsfile)) == NULL) {
      fprintf(stderr,"Can't write to file %s\n",offsetsfile);
      exit(9);
    } else {
      FWRITE_UINTS(offsets,oligospace+1,offsets_fp);
    }

  } else {
    gammaptrs = (Gammaptr_T *) CALLOC(oligospace/blocksize+1,sizeof(Gammaptr_T));
    gammaptri = 0;

    if ((offsetscomp_fp = FOPEN_WRITE_BINARY(offsetsfile)) == NULL) {
      fprintf(stderr,"Can't write to file %s\n",offsetsfile);
      exit(9);
    }
    offsets_buffer = (Offsetscomp_T *) CALLOC(offsets_buffer_size,sizeof(Offsetscomp_T));
    offsets_buffer_i = 0;

    nwritten = 0U;
    for (oligoi = 0; oligoi < oligospace; oligoi += blocksize) {
      gammaptrs[gammaptri++] = nwritten;

      offsets_buffer[offsets_buffer_i++] = offsets[oligoi];
      if (offsets_buffer_i == offsets_buffer_size) {
	FWRITE_UINTS(offsets_buffer,offsets_buffer_size,offsetscomp_fp);
	offsets_buffer_i = 0;
      }
      nwritten += 1;

      if (blocksize > 1) {
	buffer = 0U;
	ctr = 32;
	for (j = 1; j < blocksize; j++) {
	  ctr = write_gamma(offsetscomp_fp,offsets_buffer,offsets_buffer_size,&offsets_buffer_i,
			    &nwritten,&buffer,ctr,offsets[oligoi+j]-offsets[oligoi+j-1]);
	}
	debug(printf("writing gamma %08X\n",buffer));
	offsets_buffer[offsets_buffer_i++] = buffer;
	if (offsets_buffer_i == offsets_buffer_size) {
	  FWRITE_UINTS(offsets_buffer,offsets_buffer_size,offsetscomp_fp);
	  offsets_buffer_i = 0;
	}
	nwritten += 1;
      }
    }


    /* Final entries for i == oligospace */
    gammaptrs[gammaptri++] = nwritten;
    if ((gammaptrs_fp = FOPEN_WRITE_BINARY(gammaptrsfile)) == NULL) {
      fprintf(stderr,"Can't write to file %s\n",gammaptrsfile);
      exit(9);
    } else {
      FWRITE_UINTS(gammaptrs,gammaptri,gammaptrs_fp);
      FREE(gammaptrs);
      fclose(gammaptrs_fp);
    }
    
    offsets_buffer[offsets_buffer_i++] = offsets[oligoi];
    if (offsets_buffer_i > 0) {
      FWRITE_UINTS(offsets_buffer,offsets_buffer_i,offsetscomp_fp);
      offsets_buffer_i = 0;
    }
    FREE(offsets_buffer);
    fclose(offsetscomp_fp);
    nwritten += 1;

  }

  return;
}


#if 0
static void
check_bitpackptrs (char *bitpackptrsfile, char *offsetscompfile,
		   Oligospace_T oligospace, Blocksize_T offsetscomp_blocksize) {
  UINT4 *bitpackptrs, *offsetscomp, *info;
  int bitpackptrs_fd, offsetscomp_fd;
  size_t bitpackptrs_len, offsetscomp_len;
  UINT4 *blockptr;
  Oligospace_T oligo;
  int packsize;
#ifndef HAVE_MMAP
  double seconds;
#endif

#ifdef HAVE_MMAP
  bitpackptrs = (UINT4 *) Access_mmap(&bitpackptrs_fd,&bitpackptrs_len,bitpackptrsfile,sizeof(UINT4),/*randomp*/false);
  offsetscomp = (UINT4 *) Access_mmap(&offsetscomp_fd,&offsetscomp_len,offsetscompfile,sizeof(UINT4),/*randomp*/false);
#else
  bitpackptrs = (UINT4 *) Access_allocated(&bitpackptrs_len,&seconds,bitpackptrsfile,sizeof(UINT4));
  offsetscomp = (UINT4 *) Access_allocated(&offsetscomp_len,&seconds,offsetscompfile,sizeof(UINT4));
#endif

  for (oligo = 0; oligo < oligospace; oligo += offsetscomp_blocksize) {
    info = &(bitpackptrs[oligo/offsetscomp_blocksize * METAINFO_SIZE]);
    blockptr = &(offsetscomp[info[0]]);
    packsize = (info[2] - info[0])/2;
    if (packsize > 32) {
      fprintf(stderr,"Error:\n");
      fprintf(stderr,"oligo: %08X\n",oligo);
      fprintf(stderr,"nwritten: %u\n",info[0]);
      fprintf(stderr,"packsize: %d\n",packsize);
      fprintf(stderr,"offset: %u\n",info[1]);
      abort();
    }
  }

#ifdef HAVE_MMAP
  munmap((void *) offsetscomp,offsetscomp_len);
  munmap((void *) bitpackptrs,bitpackptrs_len);
#else
  FREE(offsetscomp);
  FREE(bitpackptrs);
#endif

  return;
}
#endif


#ifndef PMAP
static void
check_offsets_from_bitpack (char *bitpackptrsfile, char *bitpackcompfile, Positionsptr_T *offsets,
			    Oligospace_T oligospace, Blocksize_T blocksize) {
  UINT4 *bitpackptrs;
  UINT4 *bitpackcomp;
  Positionsptr_T offsets64[65];
  int bitpackptrs_fd, bitpackcomp_fd;
  size_t bitpackptrs_len, bitpackcomp_len;
  Oligospace_T oligoi, i;
#ifndef HAVE_MMAP
  double seconds;
#endif


#ifdef HAVE_MMAP
  bitpackptrs = (UINT4 *) Access_mmap(&bitpackptrs_fd,&bitpackptrs_len,bitpackptrsfile,sizeof(UINT4),/*randomp*/false);
  bitpackcomp = (UINT4 *) Access_mmap(&bitpackcomp_fd,&bitpackcomp_len,bitpackcompfile,sizeof(UINT4),/*randomp*/false);
#else
  bitpackptrs = (UINT4 *) Access_allocated(&bitpackptrs_len,&seconds,bitpackptrsfile,sizeof(UINT4));
  bitpackcomp = (UINT4 *) Access_allocated(&bitpackcomp_len,&seconds,bitpackcompfile,sizeof(UINT4));
#endif

  Bitpack64_read_setup();

  for (oligoi = 0UL; oligoi < oligospace; oligoi += blocksize) {
    Bitpack64_block_offsets(offsets64,oligoi,bitpackptrs,bitpackcomp);
    for (i = 0; i <= 64; i++) {
      if (offsets64[i] != offsets[oligoi+i]) {
#ifdef OLIGOSPACE_NOT_LONG
	fprintf(stderr,"\nProblem with bitpack at oligo %u+%u = %u: uncompressed %u != expected %u.  Your compiler may be defective.  Please inform twu@gene.com",
		oligoi,i,oligoi+i,offsets64[i],offsets[oligoi+i]);
#else
	fprintf(stderr,"\nProblem with bitpack at oligo %lu+%lu = %lu: uncompressed %u != expected %u.  Your compiler may be defective.  Please inform twu@gene.com",
		oligoi,i,oligoi+i,offsets64[i],offsets[oligoi+i]);
#endif
      }
    }
  }

#ifdef HAVE_MMAP
  munmap((void *) bitpackcomp,bitpackcomp_len);
  munmap((void *) bitpackptrs,bitpackptrs_len);
#else
  FREE(bitpackcomp);
  FREE(bitpackptrs);
#endif

  return;
}
#endif


#ifndef PMAP
static void
check_offsets_from_bitpack_huge (char *bitpackpagesfile, char *bitpackptrsfile, char *bitpackcompfile,
				 Hugepositionsptr_T *offsets, Oligospace_T oligospace, Blocksize_T blocksize) {
  UINT4 *bitpackpages;
  UINT4 *bitpackptrs;
  UINT4 *bitpackcomp;
  Hugepositionsptr_T offsets64[65];
  int bitpackptrs_fd, bitpackcomp_fd;
  size_t bitpackpages_len, bitpackptrs_len, bitpackcomp_len;
  Oligospace_T oligoi, i;
  double seconds;


  bitpackpages = (UINT4 *) Access_allocated(&bitpackpages_len,&seconds,bitpackpagesfile,sizeof(UINT4));
#ifdef HAVE_MMAP
  bitpackptrs = (UINT4 *) Access_mmap(&bitpackptrs_fd,&bitpackptrs_len,bitpackptrsfile,sizeof(UINT4),/*randomp*/false);
  bitpackcomp = (UINT4 *) Access_mmap(&bitpackcomp_fd,&bitpackcomp_len,bitpackcompfile,sizeof(UINT4),/*randomp*/false);
#else
  bitpackptrs = (UINT4 *) Access_allocated(&bitpackptrs_len,&seconds,bitpackptrsfile,sizeof(UINT4));
  bitpackcomp = (UINT4 *) Access_allocated(&bitpackcomp_len,&seconds,bitpackcompfile,sizeof(UINT4));
#endif

  Bitpack64_read_setup();

  for (oligoi = 0UL; oligoi < oligospace; oligoi += blocksize) {
    Bitpack64_block_offsets_huge(offsets64,oligoi,bitpackpages,bitpackptrs,bitpackcomp);
    for (i = 0; i <= 64; i++) {
      if (offsets64[i] != offsets[oligoi+i]) {
#ifdef OLIGOSPACE_NOT_LONG
	fprintf(stderr,"\nProblem with bitpack at oligo %u+%u = %u: uncompressed %lu != expected %lu.  Your compiler may be defective.  Please inform twu@gene.com\n",
		oligoi,i,oligoi+i,offsets64[i],offsets[oligoi+i]);
#else
	fprintf(stderr,"\nProblem with bitpack at oligo %lu+%lu = %lu: uncompressed %lu != expected %lu.  Your compiler may be defective.  Please inform twu@gene.com\n",
		oligoi,i,oligoi+i,offsets64[i],offsets[oligoi+i]);
#endif
      }
    }
  }

#ifdef HAVE_MMAP
  munmap((void *) bitpackcomp,bitpackcomp_len);
  munmap((void *) bitpackptrs,bitpackptrs_len);
#else
  FREE(bitpackcomp);
  FREE(bitpackptrs);
#endif
  FREE(bitpackpages);

  return;
}
#endif



#if 0
void
Indexdb_convert_gammas (char *gammaptrsfile, char *offsetscompfile, FILE *offsets_fp,
#ifdef PMAP
			int alphabet_size, Width_T index1part_aa,
#else
			Width_T index1part,
#endif
			Blocksize_T blocksize) {
  int gammaptrs_fd, offsetscomp_fd;
  Positionsptr_T *offsets = NULL, totalcounts;
  Blocksize_T j;
  Oligospace_T oligospace, oligoi;

  UINT4 buffer;
  int ctr;
  UINT4 nwritten;

#ifdef PMAP
  oligospace = power(alphabet_size,index1part_aa);
#else
  oligospace = power(4,index1part);
#endif

  offsets = (Positionsptr_T *) CALLOC(oligospace+1,sizeof(Positionsptr_T));
  FREAD_UINTS(offsets,oligospace+1,offsets_fp);
  totalcounts = offsets[oligospace];
  if (totalcounts == 0) {
    fprintf(stderr,"Something is wrong with the offsets file.  Total counts is zero.\n");
    exit(9);
  }

  gammaptrs_fd = Access_fileio_rw(gammaptrsfile);
  offsetscomp_fd = Access_fileio_rw(offsetscompfile);

  nwritten = 0U;
  for (oligoi = 0; oligoi < oligospace; oligoi += blocksize) {
    WRITE_UINT(nwritten,gammaptrs_fd);

    WRITE_UINT(offsets[oligoi],offsetscomp_fd);
    nwritten += 1;

    buffer = 0U;
    ctr = 32;
    for (j = 1; j < blocksize; j++) {
      ctr = write_gamma(offsetscomp_fp,offsets_buffer,offsets_buffer_size,&offsets_buffer_i,
			&nwritten,&buffer,ctr,offsets[oligoi+j]-offsets[oligoi+j-1]);
    }
    debug3(printf("writing gamma %08X\n",buffer));
    WRITE_UINT(buffer,offsetscomp_fd);
    nwritten += 1;
  }

  WRITE_UINT(offsets[oligoi],offsetscomp_fd);
  nwritten += 1;
  WRITE_UINT(nwritten,gammaptrs_fd);

  close(offsetscomp_fd);
  close(gammaptrs_fd);

  FREE(offsets);
  return;
}
#endif



static void
check_offsets_from_gammas (char *gammaptrsfile, char *offsetscompfile, Positionsptr_T *offsets,
			   Oligospace_T oligospace, Blocksize_T blocksize) {
  Gammaptr_T *gammaptrs;
  Offsetscomp_T *offsetscomp;
  int gammaptrs_fd, offsetscomp_fd;
  size_t gammaptrs_len, offsetscomp_len;
  Oligospace_T oligoi, oligok;
  int j, p;
#ifndef HAVE_MMAP
  double seconds;
#endif

  Positionsptr_T *ptr, cum;
  int ctr;


#ifdef HAVE_MMAP
  gammaptrs = (Gammaptr_T *) Access_mmap(&gammaptrs_fd,&gammaptrs_len,gammaptrsfile,sizeof(Gammaptr_T),/*randomp*/false);
  offsetscomp = (Offsetscomp_T *) Access_mmap(&offsetscomp_fd,&offsetscomp_len,offsetscompfile,sizeof(Offsetscomp_T),/*randomp*/false);
#else
  gammaptrs = (Gammaptr_T *) Access_allocated(&gammaptrs_len,&seconds,gammaptrsfile,sizeof(Gammaptr_T));
  offsetscomp = (Offsetscomp_T *) Access_allocated(&offsetscomp_len,&seconds,offsetscompfile,sizeof(Offsetscomp_T));
#endif

  ptr = offsetscomp;
  oligok = 0UL;
  p = 0;

  for (oligoi = 0UL; oligoi < oligospace; oligoi += blocksize) {
#ifdef HAVE_MMAP
#ifdef WORDS_BIGENDIAN
    cum = Bigendian_convert_uint(*ptr++);
#else
    cum = *ptr++;
#endif
#else
    cum = *ptr++;
#endif

    if (offsetscomp[gammaptrs[p++]] != cum) {
#ifdef OLIGOSPACE_NOT_LONG
      fprintf(stderr,"Problem with gammaptrs at oligo %u: %u != %u.  Please inform twu@gene.com\n",
	      oligok,offsetscomp[gammaptrs[p-1]],cum);
#else
      fprintf(stderr,"Problem with gammaptrs at oligo %lu: %u != %u.  Please inform twu@gene.com\n",
	      oligok,offsetscomp[gammaptrs[p-1]],cum);
#endif
      exit(9);
    }

    if (offsets[oligok++] != cum) {
#ifdef OLIGOSPACE_NOT_LONG
      fprintf(stderr,"Problem with offsetscomp at oligo %u: %u != %u.  Please inform twu@gene.com\n",
	      oligok-1U,offsets[oligok-1],cum);
#else
      fprintf(stderr,"Problem with offsetscomp at oligo %lu: %u != %u.  Please inform twu@gene.com\n",
	      oligok-1UL,offsets[oligok-1],cum);
#endif
      exit(9);
    }

    ctr = 0;
    for (j = 1; j < blocksize; j++) {
#ifdef HAVE_MMAP
#ifdef WORDS_BIGENDIAN
      ctr = Genome_read_gamma_bigendian(&ptr,ctr,&cum);
#else
      ctr = Genome_read_gamma(&ptr,ctr,&cum);
#endif
#else
      ctr = Genome_read_gamma(&ptr,ctr,&cum);
#endif

      if (offsets[oligok++] != cum) {
#ifdef OLIGOSPACE_NOT_LONG
	fprintf(stderr,"Problem with offsetscomp at oligo %u: %u != %u.  Please inform twu@gene.com\n",
		oligok-1U,offsets[oligok-1],cum);
#else
	fprintf(stderr,"Problem with offsetscomp at oligo %lu: %u != %u.  Please inform twu@gene.com\n",
		oligok-1UL,offsets[oligok-1],cum);
#endif
	exit(9);
      }
    }
    if (ctr > 0) {
      ptr++;			/* Done with last gamma byte */
    }
  }


#ifdef HAVE_MMAP
#ifdef WORDS_BIGENDIAN
  cum = Bigendian_convert_uint(*ptr++);
#else
  cum = *ptr++;
#endif
#else
  cum = *ptr++;
#endif

  if (offsetscomp[gammaptrs[p]] != cum) {
#ifdef OLIGOSPACE_NOT_LONG
    fprintf(stderr,"Problem with gammaptrs at oligo %u: %u != %u.  Please inform twu@gene.com\n",
	    oligok,offsetscomp[gammaptrs[p-1]],cum);
#else
    fprintf(stderr,"Problem with gammaptrs at oligo %lu: %u != %u.  Please inform twu@gene.com\n",
	    oligok,offsetscomp[gammaptrs[p-1]],cum);
#endif
    exit(9);
  }
    
  if (offsets[oligok] != cum) {
#ifdef OLIGOSPACE_NOT_LONG
    fprintf(stderr,"Problem with offsetscomp at oligo %u: %u != %u.  Please inform twu@gene.com\n",
	    oligok,offsets[oligok],*ptr);
#else
    fprintf(stderr,"Problem with offsetscomp at oligo %lu: %u != %u.  Please inform twu@gene.com\n",
	    oligok,offsets[oligok],*ptr);
#endif
    exit(9);
  }

#ifdef HAVE_MMAP
  munmap((void *) offsetscomp,offsetscomp_len);
  munmap((void *) gammaptrs,gammaptrs_len);
#else
  FREE(offsetscomp);
  FREE(gammaptrs);
#endif

  return;
}


#ifdef HAVE_64_BIT
UINT8
Indexdb_count_offsets (FILE *sequence_fp, Univ_IIT_T chromosome_iit,
#ifdef PMAP
		       Width_T index1part_aa, bool watsonp,
#else
		       Width_T index1part,
#endif
		       Width_T index1interval, bool genome_lc_p, char *fileroot, bool mask_lowercase_p) {
  char *uppercaseCode;
  Univcoord_T position = 0, next_chrbound;
  Chrpos_T chrpos = 0U;
  UINT8 noffsets = 0;
  int c;
  char *comma;
#ifdef PMAP
  int frame = -1, between_counter[3], in_counter[3];
  Storedoligomer_T low = 0U;
#else
  int between_counter = 0, in_counter = 0;
#endif
  int circular_typeint;
  int nchrs, chrnum;


  if (mask_lowercase_p == false) {
    uppercaseCode = UPPERCASE_U2T; /* We are reading DNA sequence */
  } else {
    uppercaseCode = NO_UPPERCASE;
  }

#ifdef PMAP
  between_counter[0] = between_counter[1] = between_counter[2] = 0;
  in_counter[0] = in_counter[1] = in_counter[2] = 0;
#endif

  circular_typeint = Univ_IIT_typeint(chromosome_iit,"circular");
  chrnum = 1;
  nchrs = Univ_IIT_total_nintervals(chromosome_iit);
  next_chrbound = Univ_IIT_next_chrbound(chromosome_iit,chrnum,circular_typeint);

  while ((c = Compress_get_char(sequence_fp,position,genome_lc_p)) != EOF) {
#ifdef PMAP
    if (++frame == 3) {
      frame = 0;
    }
    between_counter[frame] += 1;
    in_counter[frame] += 1;
#else
    between_counter++;
    in_counter++;
#endif

    if (position % MONITOR_INTERVAL == 0) {
      comma = Genomicpos_commafmt(position);
#ifdef PMAP
      fprintf(stderr,"Counting positions in genome %s (%d aa every %d aa), position %s",
	      fileroot,index1part_aa,index1interval,comma);
#else
      fprintf(stderr,"Counting positions in genome %s (%d bp every %d bp), position %s",
	      fileroot,index1part,index1interval,comma);
#endif
      FREE(comma);
#ifdef PMAP
      if (watsonp == true) {
	fprintf(stderr," (fwd)");
      } else {
	fprintf(stderr," (rev)");
      }
#endif
      fprintf(stderr,"\n");
    }

#ifdef PMAP
    switch (uppercaseCode[c]) {
    case 'A': low = (low << 2); break;
    case 'C': low = (low << 2) | 1U; break;
    case 'G': low = (low << 2) | 2U; break;
    case 'T': low = (low << 2) | 3U; break;
    case 'X': case 'N': 
      in_counter[0] = in_counter[1] = in_counter[2] = 0;
      low = 0U; 
      break;
    default: 
      if (genome_lc_p == true) {
	low = 0U;
	in_counter[0] = in_counter[1] = in_counter[2] = 0;
      } else {
	fprintf(stderr,"Bad character %c at position %u\n",c,position);
	abort();
      }
    }
#else
    switch (uppercaseCode[c]) {
    case 'A': break;
    case 'C': break;
    case 'G': break;
    case 'T': break;
    case 'X': case 'N': in_counter = 0; break;
    default: 
      if (genome_lc_p == true) {
	in_counter = 0;
      } else {
	fprintf(stderr,"Bad character %c at position %lu\n",c,position);
	abort();
      }
    }
#endif

#ifdef PMAP
    debug(printf("frame=%d char=%c bc=%d ic=%d high=%08X low=%08X\n",
		 frame,c,between_counter[frame],in_counter[frame],high,low));

    if (in_counter[frame] > 0) {
      if (watsonp == true) {
	if (Alphabet_get_codon_fwd(low) == AA_STOP) {
	  debug(printf("Resetting in_counter for frame %d to 0\n",frame));
	  in_counter[frame] = 0; 
	}
      } else {
	if (Alphabet_get_codon_rev(low) == AA_STOP) {
	  debug(printf("Resetting in_counter for frame %d to 0\n",frame));
	  in_counter[frame] = 0; 
	}
      }
    }
    if (in_counter[frame] == index1part_aa + 1) {
      if (between_counter[frame] >= index1interval) {
	noffsets += 1;
	between_counter[frame] = 0;
      }
      in_counter[frame] -= 1;
    }
#else
    if (in_counter == index1part) {
      if (
#ifdef NONMODULAR
	  between_counter >= index1interval
#else
	  (chrpos-index1part+1U) % index1interval == 0
#endif
	  ) {
	noffsets += 1;
	between_counter = 0;
      }
      in_counter--;
    }
#endif

    chrpos++;			/* Needs to go here, before we reset chrpos to 0 */
    if (position >= next_chrbound) {
#ifdef PMAP
      low = 0U;
      in_counter[0] = in_counter[1] = in_counter[2] = 0;
#else
      in_counter = 0;
#endif
      chrpos = 0U;
      chrnum++;
      while (chrnum <= nchrs && (next_chrbound = Univ_IIT_next_chrbound(chromosome_iit,chrnum,circular_typeint)) < position) {
	chrnum++;
      }
    }
    position++;
  }

  return noffsets;
}
#endif



void
Indexdb_write_offsets (char *destdir, char interval_char, FILE *sequence_fp, Univ_IIT_T chromosome_iit,
		       Width_T offsetscomp_basesize,
#ifdef PMAP
		       Alphabet_T alphabet, Width_T index1part_aa, bool watsonp,
#else
		       Width_T index1part,
#endif
		       Width_T index1interval, bool genome_lc_p, char *fileroot, bool mask_lowercase_p,
		       int compression_types) {
  char *uppercaseCode;
  char *pointersfile, *offsetsfile;
  char *filesuffix;
  FILE *offsets_fp;

  /* If offsets[oligospace] > 2^32, then will will want to allocate and write 8-mers for offsets file */
  Positionsptr_T *offsets;

  char *comma;
  int c, nchrs, chrnum;
  Oligospace_T oligospace, oligoi;
  Univcoord_T position = 0, next_chrbound;
  Chrpos_T chrpos = 0U;
#ifdef PMAP
  int frame = -1, between_counter[3], in_counter[3];
  Storedoligomer_T high = 0U, low = 0U, carry;
  Storedoligomer_T aaindex;
  int index1part_nt = 3*index1part_aa;
  int alphabet_size;
#else
  int between_counter = 0, in_counter = 0;
  Storedoligomer_T oligo = 0U, masked, mask;
  int bitpack_basesize;
#endif
#ifdef DEBUG1
  char *aa;
#endif

  int gamma_basesize;
  int offsetscomp_blocksize;
  int circular_typeint;

  if (mask_lowercase_p == false) {
    uppercaseCode = UPPERCASE_U2T; /* We are reading DNA sequence */
  } else {
    uppercaseCode = NO_UPPERCASE;
  }

#ifdef PMAP
  alphabet_size = Alphabet_get_size(alphabet);
  oligospace = power(alphabet_size,index1part_aa);
  between_counter[0] = between_counter[1] = between_counter[2] = 0;
  in_counter[0] = in_counter[1] = in_counter[2] = 0;
#ifdef OLIGOSPACE_NOT_LONG
  fprintf(stderr,"Allocating %u*%lu bytes for offsets\n",oligospace+1U,sizeof(Positionsptr_T));
#else
  fprintf(stderr,"Allocating %lu*%lu bytes for offsets\n",oligospace+1UL,sizeof(Positionsptr_T));
#endif
  offsets = (Positionsptr_T *) CALLOC_NO_EXCEPTION(oligospace+1,sizeof(Positionsptr_T));
  if (offsets == NULL) {
#ifdef OLIGOSPACE_NOT_LONG
    fprintf(stderr,"Unable to allocate %u bytes of memory, needed to build offsets with %d-mers\n",oligospace+1U,index1part_aa);
#else
    fprintf(stderr,"Unable to allocate %lu bytes of memory, needed to build offsets with %d-mers\n",oligospace+1UL,index1part_aa);
#endif
    fprintf(stderr,"Either find a computer with more RAM, or lower your value for the k-mer size\n");
    exit(9);
  }
#else
  mask = ~(~0UL << 2*index1part);
  oligospace = power(4,index1part);
#ifdef OLIGOSPACE_NOT_LONG
  fprintf(stderr,"Allocating %u*%lu bytes for offsets\n",oligospace+1U,sizeof(Positionsptr_T));
#else
  fprintf(stderr,"Allocating %lu*%lu bytes for offsets\n",oligospace+1UL,sizeof(Positionsptr_T));
#endif
  offsets = (Positionsptr_T *) CALLOC_NO_EXCEPTION(oligospace+1,sizeof(Positionsptr_T));
  if (offsets == NULL) {
#ifdef OLIGOSPACE_NOT_LONG
    fprintf(stderr,"Unable to allocate %u bytes of memory, needed to build offsets with %d-mers\n",oligospace+1U,index1part);
#else
    fprintf(stderr,"Unable to allocate %lu bytes of memory, needed to build offsets with %d-mers\n",oligospace+1UL,index1part);
#endif
    fprintf(stderr,"Either find a computer with more RAM, or lower your value for the k-mer size\n");
    exit(9);
  }
#endif

  /* Handle reference strain */
  circular_typeint = Univ_IIT_typeint(chromosome_iit,"circular");
  chrnum = 1;
  nchrs = Univ_IIT_total_nintervals(chromosome_iit);
  next_chrbound = Univ_IIT_next_chrbound(chromosome_iit,chrnum,circular_typeint);

  while ((c = Compress_get_char(sequence_fp,position,genome_lc_p)) != EOF) {

#ifdef PMAP
    if (++frame == 3) {
      frame = 0;
    }
    between_counter[frame] += 1;
    in_counter[frame] += 1;
#else
    between_counter++;
    in_counter++;
#endif

    if (position % MONITOR_INTERVAL == 0) {
      comma = Genomicpos_commafmt(position);
#ifdef PMAP
      fprintf(stderr,"Indexing offsets of oligomers in genome %s (%d aa every %d aa), position %s",
	      fileroot,index1part_aa,index1interval,comma);
#else
      fprintf(stderr,"Indexing offsets of oligomers in genome %s (%d bp every %d bp), position %s",
	      fileroot,index1part,index1interval,comma);
#endif
      FREE(comma);
#ifdef PMAP
      if (watsonp == true) {
	fprintf(stderr," (fwd)");
      } else {
	fprintf(stderr," (rev)");
      }
#endif
      fprintf(stderr,"\n");
    }

#ifdef PMAP
    carry = (low >> 30);
    switch (uppercaseCode[c]) {
    case 'A': low = (low << 2); break;
    case 'C': low = (low << 2) | 1U; break;
    case 'G': low = (low << 2) | 2U; break;
    case 'T': low = (low << 2) | 3U; break;
    case 'X': case 'N': 
      high = low = carry = 0U; 
      in_counter[0] = in_counter[1] = in_counter[2] = 0;
      break;
    default: 
      if (genome_lc_p == true) {
	high = low = carry = 0U;
	in_counter[0] = in_counter[1] = in_counter[2] = 0;
      } else {
	fprintf(stderr,"Bad character %c at position %u\n",c,position);
	abort();
      }
    }
    high = (high << 2) | carry; 
#else
    switch (uppercaseCode[c]) {
    case 'A': oligo = (oligo << 2); break;
    case 'C': oligo = (oligo << 2) | 1U; break;
    case 'G': oligo = (oligo << 2) | 2U; break;
    case 'T': oligo = (oligo << 2) | 3U; break;
    case 'X': case 'N': oligo = 0U; in_counter = 0; break;
    default: 
      if (genome_lc_p == true) {
	oligo = 0U; in_counter = 0;
      } else {
	fprintf(stderr,"Bad character %c at position %lu\n",c,position);
	abort();
      }
    }
#endif

#ifdef PMAP
    debug(printf("frame=%d char=%c bc=%d ic=%d high=%08X low=%08X\n",
		 frame,c,between_counter[frame],in_counter[frame],high,low));

    if (in_counter[frame] > 0) {
      if (watsonp == true) {
	if (Alphabet_get_codon_fwd(low) == AA_STOP) {
	  debug(printf("Resetting in_counter for frame %d to 0\n",frame));
	  in_counter[frame] = 0; 
	}
      } else {
	if (Alphabet_get_codon_rev(low) == AA_STOP) {
	  debug(printf("Resetting in_counter for frame %d to 0\n",frame));
	  in_counter[frame] = 0; 
	}
      }
    }
    if (in_counter[frame] == index1part_aa + 1) {
      if (between_counter[frame] >= index1interval) {
	aaindex = Alphabet_get_aa_index(high,low,watsonp,index1part_nt);
#ifdef OLIGOSPACE_NOT_LONG
	oligoi = (Oligospace_T) aaindex + 1U;
#else
	oligoi = (Oligospace_T) aaindex + 1UL;
#endif
	offsets[oligoi] += 1;
	debug1(
	       aa = aaindex_aa(aaindex);
	       if (watsonp == true) {
		 printf("Storing %s (%u) at %u\n",aa,aaindex,position-index1part_nt+1U);
	       } else {
		 printf("Storing %s (%u) at %u\n",aa,aaindex,position);
	       }
	       FREE(aa);
	       );
	between_counter[frame] = 0;
      }
      in_counter[frame] -= 1;
    }
#else
    if (in_counter == index1part) {
      if (
#ifdef NONMODULAR
	  between_counter >= index1interval
#else
	  (chrpos-index1part+1U) % index1interval == 0
#endif
	  ) {
	masked = oligo & mask;
#ifdef OLIGOSPACE_NOT_LONG
	oligoi = (Oligospace_T) masked + 1U;
#else
	oligoi = (Oligospace_T) masked + 1UL;
#endif
	offsets[oligoi] += 1;
	debug(printf("Found oligo %06X.  Incremented offsets for %lu to be %u\n",
		     masked,oligoi,offsets[oligoi]));
	between_counter = 0;
      }
      in_counter--;
    }
#endif

    chrpos++;			/* Needs to go here, before we reset chrpos to 0 */
    if (position >= next_chrbound) {
#ifdef PMAP
      high = low = carry = 0U;
      in_counter[0] = in_counter[1] = in_counter[2] = 0;
#else
      oligo = 0U; in_counter = 0;
#endif
      chrpos = 0U;
      chrnum++;
      while (chrnum <= nchrs && (next_chrbound = Univ_IIT_next_chrbound(chromosome_iit,chrnum,circular_typeint)) < position) {
	chrnum++;
      }
    }
    position++;
  }


#ifdef ADDSENTINEL
  for (oligoi = 1; oligoi <= oligospace; oligoi++) {
    offsets[oligoi] = offsets[oligoi] + offsets[oligoi-1] + 1U;
    debug(if (offsets[oligoi] != offsets[oligoi-1]) {
	    printf("Offset for %06X: %u\n",oligoi,offsets[oligoi]);
	  });
  }
#else
  for (oligoi = 1; oligoi <= oligospace; oligoi++) {
    offsets[oligoi] = offsets[oligoi] + offsets[oligoi-1];
    debug(if (offsets[oligoi] != offsets[oligoi-1]) {
	    printf("Offset for %06X: %u\n",oligoi,offsets[oligoi]);
	  });
  }
#endif

  /*
  fprintf(stderr,"Offset for A...A is %u to %u\n",offsets[0],offsets[1]);
  fprintf(stderr,"Offset for T...T is %u to %u\n",offsets[oligospace-1],offsets[oligospace]);
  */

#ifdef PMAP
  if (watsonp == true) {
    filesuffix = FWD_FILESUFFIX;
  } else {
    filesuffix = REV_FILESUFFIX;
  }
#else
  filesuffix = IDX_FILESUFFIX;
#endif

  if (compression_types == NO_COMPRESSION) {
    /* New convention.  Previously this was called offsetscomp.  Will
       need to create a symbolic link for backward compatibility */
    offsetsfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
				      strlen(".")+strlen(filesuffix)+
				      /*for kmer*/2+/*for interval char*/1+
				      strlen("offsets")+1,sizeof(char));
#ifdef PMAP
    sprintf(offsetsfile,"%s/%s.%s.%s%d%c%s",
	    destdir,fileroot,Alphabet_string(alphabet),filesuffix,index1part_aa,interval_char,"offsets");
#else
    sprintf(offsetsfile,"%s/%s.%s%02d%c%s",
	    destdir,fileroot,filesuffix,index1part,interval_char,"offsets");
#endif

    if ((offsets_fp = FOPEN_WRITE_BINARY(offsetsfile)) == NULL) {
      fprintf(stderr,"Can't write to file %s\n",offsetsfile);
      exit(9);
    } else {
      fprintf(stderr,"Writing %lu offsets to file with total of %u k-mers...",oligospace+1,offsets[oligospace]);
      FWRITE_UINTS(offsets,oligospace+1,offsets_fp);
      fprintf(stderr,"done\n");
    }
    fclose(offsets_fp);
    FREE(offsetsfile);

  } else {
    if (compression_types & GAMMA_COMPRESSION) {
      gamma_basesize = offsetscomp_basesize;
#ifdef PMAP
      offsetscomp_blocksize = power(alphabet_size,index1part_aa - gamma_basesize);
#else
      offsetscomp_blocksize = power(4,index1part - gamma_basesize);
#endif

      pointersfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
				     strlen(".")+strlen(filesuffix)+
				     /*for basesize*/2+/*for kmer*/2+/*for interval char*/1+
				     strlen("gammaptrs")+1,sizeof(char));
#ifdef PMAP
      sprintf(pointersfile,"%s/%s.%s.%s%d%d%c%s",
	      destdir,fileroot,Alphabet_string(alphabet),filesuffix,
	      gamma_basesize,index1part_aa,interval_char,"gammaptrs");
#else
      sprintf(pointersfile,"%s/%s.%s%02d%02d%c%s",
	      destdir,fileroot,filesuffix,gamma_basesize,index1part,interval_char,"gammaptrs");
#endif
	
      offsetsfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
					strlen(".")+strlen(filesuffix)+
					/*for basesize*/2+/*for kmer*/2+/*for interval char*/1+
					strlen("offsetscomp")+1,sizeof(char));
#ifdef PMAP
      sprintf(offsetsfile,"%s/%s.%s.%s%d%d%c%s",
	      destdir,fileroot,Alphabet_string(alphabet),filesuffix,
	      gamma_basesize,index1part_aa,interval_char,"offsetscomp");
#else
      sprintf(offsetsfile,"%s/%s.%s%02d%02d%c%s",
	      destdir,fileroot,filesuffix,gamma_basesize,index1part,interval_char,"offsetscomp");
#endif


#ifdef OLIGOSPACE_NOT_LONG
      fprintf(stderr,"Writing %u offsets compressed via gammas to file with total of %u k-mers...",oligospace+1U,offsets[oligospace]);
#else
      fprintf(stderr,"Writing %lu offsets compressed via gammas to file with total of %u k-mers...",oligospace+1UL,offsets[oligospace]);
#endif

      Indexdb_write_gammaptrs(pointersfile,offsetsfile,offsets,oligospace,offsetscomp_blocksize);
      if (offsetscomp_blocksize > 1) {
	fprintf(stderr,"Checking gammas...");
	check_offsets_from_gammas(pointersfile,offsetsfile,offsets,oligospace,offsetscomp_blocksize);
	fprintf(stderr,"done\n");
      }

      FREE(offsetsfile);
      FREE(pointersfile);
    }

    if (compression_types & BITPACK64_COMPRESSION) {
#ifdef PMAP
#else
      bitpack_basesize = index1part - 3;
      offsetscomp_blocksize = power(4,index1part - offsetscomp_basesize);
      if (offsetscomp_basesize != bitpack_basesize) {
	fprintf(stderr,"For bitpack compression, basesize must be 3 less than blocksize.  Changing requested basesize to %d\n",
		bitpack_basesize);
      }
#endif

      pointersfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
				     strlen(".")+strlen(filesuffix)+
				     /*for basesize*/2+/*for kmer*/2+/*for interval char*/1+
				     strlen("bitpackptrs")+1,sizeof(char));
#ifdef PMAP
      sprintf(pointersfile,"%s/%s.%s.%s%d%c%s",
	      destdir,fileroot,Alphabet_string(alphabet),filesuffix,index1part_aa,interval_char,"bitpackptrs");
#else
      sprintf(pointersfile,"%s/%s.%s%02d%02d%c%s",
	      destdir,fileroot,filesuffix,bitpack_basesize,index1part,interval_char,"bitpackptrs");
#endif
	
      offsetsfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
					strlen(".")+strlen(filesuffix)+
					/*for basesize*/2+/*for kmer*/2+/*for interval char*/1+
					strlen("bitpackcomp")+1,sizeof(char));
#ifdef PMAP
      sprintf(offsetsfile,"%s/%s.%s.%s%d%c%s",
	      destdir,fileroot,Alphabet_string(alphabet),filesuffix,index1part_aa,interval_char,"bitpackcomp");
#else
      sprintf(offsetsfile,"%s/%s.%s%02d%02d%c%s",
	      destdir,fileroot,filesuffix,bitpack_basesize,index1part,interval_char,"bitpackcomp");
#endif

#ifdef OLIGOSPACE_NOT_LONG
      fprintf(stderr,"Writing %u offsets compressed via bitpack to file with total of %u k-mers...",oligospace+1U,offsets[oligospace]);
#else
      fprintf(stderr,"Writing %lu offsets compressed via bitpack to file with total of %u k-mers...",oligospace+1UL,offsets[oligospace]);
#endif

      Bitpack64_write_differential(pointersfile,offsetsfile,offsets,oligospace);
#ifdef PMAP
      /* 20-mers do not coincide with 64-element blocks */
#else
      if (offsets[oligospace] == 0) {
	/* Don't check bitpack, which will crash on an empty file */
      } else if (offsetscomp_blocksize > 1) {
	fprintf(stderr,"Checking bitpack...");
	check_offsets_from_bitpack(pointersfile,offsetsfile,offsets,oligospace,offsetscomp_blocksize);
      }
      fprintf(stderr,"done\n");
#endif

      FREE(offsetsfile);
      FREE(pointersfile);
    }
  }
  
  FREE(offsets);

  return;
}


#ifdef HAVE_64_BIT
void
Indexdb_write_offsets_huge (char *destdir, char interval_char, FILE *sequence_fp, Univ_IIT_T chromosome_iit,
			    Width_T offsetscomp_basesize,
#ifdef PMAP
			    Alphabet_T alphabet, Width_T index1part_aa, bool watsonp,
#else
			    Width_T index1part,
#endif
			    Width_T index1interval, bool genome_lc_p, char *fileroot, bool mask_lowercase_p,
			    int compression_types) {
  char *uppercaseCode;
  char *pagesfile, *pointersfile, *offsetsfile;
  char *filesuffix;

  /* If offsets[oligospace] > 2^32, then we allocate and write 8-mers for offsets file */
  Hugepositionsptr_T *offsets;

  char *comma;
  int c, nchrs, chrnum;
  Oligospace_T oligospace, oligoi;
  Univcoord_T position = 0, next_chrbound;
  Chrpos_T chrpos = 0U;
#ifdef PMAP
  int alphabet_size;
  int frame = -1, between_counter[3], in_counter[3];
  Storedoligomer_T high = 0U, low = 0U, carry;
  Storedoligomer_T aaindex;
  int index1part_nt = 3*index1part_aa;
#else
  int between_counter = 0, in_counter = 0;
  Storedoligomer_T oligo = 0U, masked, mask;
#endif
#ifdef DEBUG1
  char *aa;
#endif

  int bitpack_basesize;
  int offsetscomp_blocksize;
  int circular_typeint;

  if (mask_lowercase_p == false) {
    uppercaseCode = UPPERCASE_U2T; /* We are reading DNA sequence */
  } else {
    uppercaseCode = NO_UPPERCASE;
  }

#ifdef PMAP
  alphabet_size = Alphabet_get_size(alphabet);
  oligospace = power(alphabet_size,index1part_aa);
  between_counter[0] = between_counter[1] = between_counter[2] = 0;
  in_counter[0] = in_counter[1] = in_counter[2] = 0;
#ifdef OLIGOSPACE_NOT_LONG
  fprintf(stderr,"Allocating %u*%lu bytes for offsets\n",oligospace+1U,sizeof(Hugepositionsptr_T));
#else
  fprintf(stderr,"Allocating %lu*%lu bytes for offsets\n",oligospace+1UL,sizeof(Hugepositionsptr_T));
#endif
  offsets = (Hugepositionsptr_T *) CALLOC_NO_EXCEPTION(oligospace+1,sizeof(Hugepositionsptr_T));
  if (offsets == NULL) {
#ifdef OLIGOSPACE_NOT_LONG
    fprintf(stderr,"Unable to allocate %u bytes of memory, needed to build offsets with %d-mers\n",oligospace+1U,index1part_aa);
#else
    fprintf(stderr,"Unable to allocate %lu bytes of memory, needed to build offsets with %d-mers\n",oligospace+1UL,index1part_aa);
#endif
    fprintf(stderr,"Either find a computer with more RAM, or lower your value for the k-mer size\n");
    exit(9);
  }
#else
  mask = ~(~0UL << 2*index1part);
  oligospace = power(4,index1part);
#ifdef OLIGOSPACE_NOT_LONG
  fprintf(stderr,"Allocating %u*%lu bytes for offsets\n",oligospace+1U,sizeof(Hugepositionsptr_T));
#else
  fprintf(stderr,"Allocating %lu*%lu bytes for offsets\n",oligospace+1UL,sizeof(Hugepositionsptr_T));
#endif
  offsets = (Hugepositionsptr_T *) CALLOC_NO_EXCEPTION(oligospace+1,sizeof(Hugepositionsptr_T));
  if (offsets == NULL) {
#ifdef OLIGOSPACE_NOT_LONG
    fprintf(stderr,"Unable to allocate %u bytes of memory, needed to build offsets with %d-mers\n",oligospace+1U,index1part);
#else
    fprintf(stderr,"Unable to allocate %lu bytes of memory, needed to build offsets with %d-mers\n",oligospace+1UL,index1part);
#endif
    fprintf(stderr,"Either find a computer with more RAM, or lower your value for the k-mer size\n");
    exit(9);
  }
#endif

  /* Handle reference strain */
  circular_typeint = Univ_IIT_typeint(chromosome_iit,"circular");
  chrnum = 1;
  nchrs = Univ_IIT_total_nintervals(chromosome_iit);
  next_chrbound = Univ_IIT_next_chrbound(chromosome_iit,chrnum,circular_typeint);

  while ((c = Compress_get_char(sequence_fp,position,genome_lc_p)) != EOF) {
#ifdef PMAP
    if (++frame == 3) {
      frame = 0;
    }
    between_counter[frame] += 1;
    in_counter[frame] += 1;
#else
    between_counter++;
    in_counter++;
#endif

    if (position % MONITOR_INTERVAL == 0) {
      comma = Genomicpos_commafmt(position);
#ifdef PMAP
      fprintf(stderr,"Indexing offsets of oligomers in genome %s (%d aa every %d aa), position %s",
	      fileroot,index1part_aa,index1interval,comma);
#else
      fprintf(stderr,"Indexing offsets of oligomers in genome %s (%d bp every %d bp), position %s",
	      fileroot,index1part,index1interval,comma);
#endif
      FREE(comma);
#ifdef PMAP
      if (watsonp == true) {
	fprintf(stderr," (fwd)");
      } else {
	fprintf(stderr," (rev)");
      }
#endif
      fprintf(stderr,"\n");
    }

#ifdef PMAP
    carry = (low >> 30);
    switch (uppercaseCode[c]) {
    case 'A': low = (low << 2); break;
    case 'C': low = (low << 2) | 1U; break;
    case 'G': low = (low << 2) | 2U; break;
    case 'T': low = (low << 2) | 3U; break;
    case 'X': case 'N': 
      high = low = carry = 0U; 
      in_counter[0] = in_counter[1] = in_counter[2] = 0;
      break;
    default: 
      if (genome_lc_p == true) {
	high = low = carry = 0U;
	in_counter[0] = in_counter[1] = in_counter[2] = 0;
      } else {
	fprintf(stderr,"Bad character %c at position %u\n",c,position);
	abort();
      }
    }
    high = (high << 2) | carry; 
#else
    switch (uppercaseCode[c]) {
    case 'A': oligo = (oligo << 2); break;
    case 'C': oligo = (oligo << 2) | 1U; break;
    case 'G': oligo = (oligo << 2) | 2U; break;
    case 'T': oligo = (oligo << 2) | 3U; break;
    case 'X': case 'N': oligo = 0U; in_counter = 0; break;
    default: 
      if (genome_lc_p == true) {
	oligo = 0U; in_counter = 0;
      } else {
	fprintf(stderr,"Bad character %c at position %lu\n",c,position);
	abort();
      }
    }
#endif

#ifdef PMAP
    debug(printf("frame=%d char=%c bc=%d ic=%d high=%08X low=%08X\n",
		 frame,c,between_counter[frame],in_counter[frame],high,low));

    if (in_counter[frame] > 0) {
      if (watsonp == true) {
	if (Alphabet_get_codon_fwd(low) == AA_STOP) {
	  debug(printf("Resetting in_counter for frame %d to 0\n",frame));
	  in_counter[frame] = 0; 
	}
      } else {
	if (Alphabet_get_codon_rev(low) == AA_STOP) {
	  debug(printf("Resetting in_counter for frame %d to 0\n",frame));
	  in_counter[frame] = 0; 
	}
      }
    }
    if (in_counter[frame] == index1part_aa + 1) {
      if (between_counter[frame] >= index1interval) {
	aaindex = Alphabet_get_aa_index(high,low,watsonp,index1part_nt);
#ifdef OLIGOSPACE_NOT_LONG
	oligoi = (Oligospace_T) aaindex + 1U;
#else
	oligoi = (Oligospace_T) aaindex + 1UL;
#endif
	offsets[oligoi] += 1;
	debug1(
	       aa = aaindex_aa(aaindex);
	       if (watsonp == true) {
		 printf("Storing %s (%u) at %u\n",aa,aaindex,position-index1part_nt+1U);
	       } else {
		 printf("Storing %s (%u) at %u\n",aa,aaindex,position);
	       }
	       FREE(aa);
	       );
	between_counter[frame] = 0;
      }
      in_counter[frame] -= 1;
    }
#else
    if (in_counter == index1part) {
      if (
#ifdef NONMODULAR
	  between_counter >= index1interval
#else
	  (chrpos-index1part+1U) % index1interval == 0
#endif
	  ) {
	masked = oligo & mask;
#ifdef OLIGOSPACE_NOT_LONG
	oligoi = (Oligospace_T) masked + 1U;
#else
	oligoi = (Oligospace_T) masked + 1UL;
#endif
	offsets[oligoi] += 1;
	debug(printf("Found oligo %06X.  Incremented offsets for %lu to be %u\n",
		     masked,oligoi,offsets[oligoi]));
	between_counter = 0;
      }
      in_counter--;
    }
#endif

    chrpos++;			/* Needs to go here, before we reset chrpos to 0 */
    if (position >= next_chrbound) {
#ifdef PMAP
      high = low = carry = 0U;
      in_counter[0] = in_counter[1] = in_counter[2] = 0;
#else
      oligo = 0U; in_counter = 0;
#endif
      chrpos = 0U;
      chrnum++;
      while (chrnum <= nchrs && (next_chrbound = Univ_IIT_next_chrbound(chromosome_iit,chrnum,circular_typeint)) < position) {
	chrnum++;
      }
    }
    position++;
  }


#ifdef ADDSENTINEL
  for (oligoi = 1; oligoi <= oligospace; oligoi++) {
    offsets[oligoi] = offsets[oligoi] + offsets[oligoi-1] + 1U;
    debug(if (offsets[oligoi] != offsets[oligoi-1]) {
	    printf("Offset for %06X: %u\n",oligoi,offsets[oligoi]);
	  });
  }
#else
  for (oligoi = 1; oligoi <= oligospace; oligoi++) {
    offsets[oligoi] = offsets[oligoi] + offsets[oligoi-1];
    debug(if (offsets[oligoi] != offsets[oligoi-1]) {
	    printf("Offset for %06X: %u\n",oligoi,offsets[oligoi]);
	  });
  }
#endif

  /*
  fprintf(stderr,"Offset for A...A is %u to %u\n",offsets[0],offsets[1]);
  fprintf(stderr,"Offset for T...T is %u to %u\n",offsets[oligospace-1],offsets[oligospace]);
  */

#ifdef PMAP
  if (watsonp == true) {
    filesuffix = FWD_FILESUFFIX;
  } else {
    filesuffix = REV_FILESUFFIX;
  }
#else
  filesuffix = IDX_FILESUFFIX;
#endif

  if (compression_types == NO_COMPRESSION) {
    fprintf(stderr,"Not supporting non-compression on huge genomes\n");
    abort();

  } else {
    if (compression_types & GAMMA_COMPRESSION) {
      /* Don't need to make huge genomes backwards compatible */
      fprintf(stderr,"Not supporting gamma compression on huge genomes\n");
      abort();
    }

    if (compression_types & BITPACK64_COMPRESSION) {
#ifdef PMAP
      bitpack_basesize = index1part_aa - 3;
      offsetscomp_blocksize = power(alphabet_size,index1part_aa - offsetscomp_basesize);
#else
      bitpack_basesize = index1part - 3;
      offsetscomp_blocksize = power(4,index1part - offsetscomp_basesize);
#endif
      if (offsetscomp_basesize != bitpack_basesize) {
	fprintf(stderr,"For bitpack compression, basesize must be 3 less than blocksize.  Changing requested basesize to %d\n",
		bitpack_basesize);
      }

      pagesfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
				  strlen(".")+strlen(filesuffix)+
				  /*for basesize*/2+/*for kmer*/2+/*for interval char*/1+
				  strlen("bitpackpages")+1,sizeof(char));
#ifdef PMAP
      sprintf(pagesfile,"%s/%s.%s.%s%d%c%s",
	      destdir,fileroot,Alphabet_string(alphabet),filesuffix,index1part_aa,interval_char,"bitpackpages");
#else
      sprintf(pagesfile,"%s/%s.%s%02d%02d%c%s",
	      destdir,fileroot,filesuffix,bitpack_basesize,index1part,interval_char,"bitpackpages");
#endif

      pointersfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
				     strlen(".")+strlen(filesuffix)+
				     /*for basesize*/2+/*for kmer*/2+/*for interval char*/1+
				     strlen("bitpackptrs")+1,sizeof(char));
#ifdef PMAP
      sprintf(pointersfile,"%s/%s.%s.%s%d%c%s",
	      destdir,fileroot,Alphabet_string(alphabet),filesuffix,index1part_aa,interval_char,"bitpackptrs");
#else
      sprintf(pointersfile,"%s/%s.%s%02d%02d%c%s",
	      destdir,fileroot,filesuffix,bitpack_basesize,index1part,interval_char,"bitpackptrs");
#endif
	
      offsetsfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
					strlen(".")+strlen(filesuffix)+
					/*for basesize*/2+/*for kmer*/2+/*for interval char*/1+
					strlen("bitpackcomp")+1,sizeof(char));
#ifdef PMAP
      sprintf(offsetsfile,"%s/%s.%s.%s%d%c%s",
	      destdir,fileroot,Alphabet_string(alphabet),filesuffix,index1part_aa,interval_char,"bitpackcomp");
#else
      sprintf(offsetsfile,"%s/%s.%s%02d%02d%c%s",
	      destdir,fileroot,filesuffix,bitpack_basesize,index1part,interval_char,"bitpackcomp");
#endif

#ifdef OLIGOSPACE_NOT_LONG
      fprintf(stderr,"Writing %u offsets compressed via bitpack to file with total of %lu k-mers...",oligospace+1U,offsets[oligospace]);
#else
      fprintf(stderr,"Writing %lu offsets compressed via bitpack to file with total of %lu k-mers...",oligospace+1UL,offsets[oligospace]);
#endif

      Bitpack64_write_differential_huge(pagesfile,pointersfile,offsetsfile,offsets,oligospace);
#ifndef PMAP
      if (offsets[oligospace] == 0) {
	/* Don't check bitpack, which will crash on an empty file */
      } else if (offsetscomp_blocksize > 1) {
	fprintf(stderr,"Checking bitpack...");
	check_offsets_from_bitpack_huge(pagesfile,pointersfile,offsetsfile,offsets,oligospace,offsetscomp_blocksize);
      }
      fprintf(stderr,"done\n");
#endif

      FREE(offsetsfile);
      FREE(pointersfile);
      FREE(pagesfile);
    }
  }
  
  FREE(offsets);

  return;
}
#endif


/* FILE *fp is preferable to int fd, because former is buffered.  No
   need for fseeko, because offsets file is < 2 Gigabytes */
#if 0
static void
offsetsfile_move_absolute (FILE *fp, int ptr) {
  long int offset = ptr*((long int) sizeof(Positionsptr_T));

  if (fseek(fp,offset,SEEK_SET) < 0) {
    fprintf(stderr,"Attempted to do fseek on offset %u*%lu=%lu\n",ptr,sizeof(Positionsptr_T),offset);
    perror("Error in indexdb.c, offsetsfile_move_absolute");
    exit(9);
  }

  return;
}
#endif


#if 0
static bool
need_to_sort_p (Univcoord_T *positions, int length) {
  Univcoord_T prevpos;
  int j;

  prevpos = positions[0];
  for (j = 1; j < length; j++) {
    if (positions[j] <= prevpos) {
      return true;
    } else {
      prevpos = positions[j];
    }
  }
  return false;
}
#endif


static void
positions_move_absolute_1 (int positions_fd, Positionsptr_T ptr) {
  off_t offset = ptr*((off_t) sizeof(unsigned char));

  if (lseek(positions_fd,offset,SEEK_SET) < 0) {
    fprintf(stderr,"Attempted to do lseek on offset %u*%lu=%lu\n",ptr,sizeof(unsigned char),offset);
    perror("Error in indexdb.c, positions_move_absolute_1");
    exit(9);
  }
  return;
}

static void
positions_move_absolute_4 (int positions_fd, Positionsptr_T ptr) {
  off_t offset = ptr*((off_t) sizeof(UINT4));

  if (lseek(positions_fd,offset,SEEK_SET) < 0) {
    fprintf(stderr,"Attempted to do lseek on offset %u*%lu=%lu\n",ptr,sizeof(UINT4),offset);
    perror("Error in indexdb.c, positions_move_absolute_4");
    exit(9);
  }
  return;
}

#if 0
/* Replaced by positions_move_absolute_1 and positions_move_absolute_4 */
static void
positions_move_absolute_8 (int positions_fd, Positionsptr_T ptr) {
  off_t offset = ptr*((off_t) sizeof(UINT8));

  if (lseek(positions_fd,offset,SEEK_SET) < 0) {
    fprintf(stderr,"Attempted to do lseek on offset %u*%lu=%lu\n",ptr,sizeof(UINT8),offset);
    perror("Error in indexdb.c, positions_move_absolute_8");
    exit(9);
  }
  return;
}
#endif


/* Works directly in file, so we don't need to allocate memory */
static void
compute_positions_in_file (int positions_high_fd, int positions_low_fd, Positionsptr_T *offsets,
			   FILE *sequence_fp, Univ_IIT_T chromosome_iit,
#ifdef PMAP
			   Width_T index1part_aa, bool watsonp,
#else
			   Width_T index1part,
#endif
			   Width_T index1interval, bool genome_lc_p, char *fileroot,
			   bool mask_lowercase_p, bool coord_values_8p) {
  char *uppercaseCode;
  Univcoord_T position = 0, next_chrbound;
  UINT8 adjposition8;
  UINT4 adjposition4;
  Chrpos_T chrpos = 0U;
  char *comma;
  int c, nchrs, chrnum;
#ifdef PMAP
  int frame = -1, between_counter[3], in_counter[3];
  Storedoligomer_T high = 0U, low = 0U, carry;
  Oligospace_T aaindex;
  Width_T index1part_nt = 3*index1part_aa;
#else
  int between_counter = 0, in_counter = 0;
  Storedoligomer_T oligo = 0U, masked, mask;
#endif
  int circular_typeint;

#ifdef ADDSENTINEL
  Oligospace_T oligospace, oligoi;
  oligospace = power(4,index1part);
#endif

  if (mask_lowercase_p == false) {
    uppercaseCode = UPPERCASE_U2T; /* We are reading DNA sequence */
  } else {
    uppercaseCode = NO_UPPERCASE;
  }

#ifdef PMAP
  /* oligospace = power(alphabet_size,index1part_aa); */
#else
  mask = ~(~0UL << 2*index1part);
  /* oligospace = power(4,index1part); */
#endif

  /* Handle reference strain */
  circular_typeint = Univ_IIT_typeint(chromosome_iit,"circular");
  chrnum = 1;
  nchrs = Univ_IIT_total_nintervals(chromosome_iit);
  next_chrbound = Univ_IIT_next_chrbound(chromosome_iit,chrnum,circular_typeint);

  while ((c = Compress_get_char(sequence_fp,position,genome_lc_p)) != EOF) {
#ifdef PMAP
    if (++frame == 3) {
      frame = 0;
    }
    between_counter[frame] += 1;
    in_counter[frame] += 1;
#else
    between_counter++;
    in_counter++;
#endif

    if (position % MONITOR_INTERVAL == 0) {
      comma = Genomicpos_commafmt(position);
#ifdef PMAP
      fprintf(stderr,"Indexing positions of oligomers in genome %s (%d aa every %d aa), position %s",
	      fileroot,index1part_aa,index1interval,comma);
#else
      fprintf(stderr,"Indexing positions of oligomers in genome %s (%d bp every %d bp), position %s",
	      fileroot,index1part,index1interval,comma);
#endif
      FREE(comma);
#ifdef PMAP
      if (watsonp == true) {
	fprintf(stderr," (fwd)");
      } else {
	fprintf(stderr," (rev)");
      }
#endif
      fprintf(stderr,"\n");
    }

#ifdef PMAP
    carry = (low >> 30);
    switch (uppercaseCode[c]) {
    case 'A': low = (low << 2); break;
    case 'C': low = (low << 2) | 1; break;
    case 'G': low = (low << 2) | 2; break;
    case 'T': low = (low << 2) | 3; break;
    case 'X': case 'N': 
      high = low = carry = 0U;
      in_counter[0] = in_counter[1] = in_counter[2] = 0;
      break;
    default: 
      if (genome_lc_p == true) {
	high = low = carry = 0U;
	in_counter[0] = in_counter[1] = in_counter[2] = 0;
      } else {
	fprintf(stderr,"Bad character %c at position %u\n",c,position);
	abort();
      }
    }
    high = (high << 2) | carry; 
#else
    switch (uppercaseCode[c]) {
    case 'A': oligo = (oligo << 2); break;
    case 'C': oligo = (oligo << 2) | 1; break;
    case 'G': oligo = (oligo << 2) | 2; break;
    case 'T': oligo = (oligo << 2) | 3; break;
    case 'X': case 'N': oligo = 0U; in_counter = 0; break;
    default: 
      if (genome_lc_p == true) {
	oligo = 0U; in_counter = 0;
      } else {
	fprintf(stderr,"Bad character %c at position %lu\n",c,position);
	abort();
      }
    }
#endif

    /*
    debug(printf("char=%c bc=%d ic=%d oligo=%016lX\n",
		 c,between_counter,in_counter,oligo));
    */
    
#ifdef PMAP
    if (in_counter[frame] > 0) {
      if (watsonp == true) {
	if (Alphabet_get_codon_fwd(low) == AA_STOP) {
	  in_counter[frame] = 0; 
	}
      } else {
	if (Alphabet_get_codon_rev(low) == AA_STOP) {
	  in_counter[frame] = 0; 
	}
      }
    }
    if (in_counter[frame] == index1part_aa + 1) {
      if (between_counter[frame] >= index1interval) {
	aaindex = Alphabet_get_aa_index(high,low,watsonp,index1part_nt);
	if (coord_values_8p == true) {
	  positions_move_absolute_1(positions_high_fd,offsets[aaindex]);
	  positions_move_absolute_4(positions_low_fd,offsets[aaindex]);
	  offsets[aaindex] += 1;
	  if (watsonp == true) {
	    adjposition8 =  position - index1part_nt + 1;
	  } else {
	    adjposition8 = position;
	  }
	  WRITE_CHAR((unsigned char) (adjposition8 >> POSITIONS8_HIGH_SHIFT),positions_high_fd);
	  WRITE_UINT((UINT4) (adjposition8 & POSITIONS8_LOW_MASK),positions_low_fd);
	} else {
	  positions_move_absolute_4(positions_low_fd,offsets[aaindex]);
	  offsets[aaindex] += 1;
	  if (watsonp == true) {
	    adjposition4 = (UINT4) (position - index1part_nt + 1);
	  } else {
	    adjposition4 = (UINT4) position;
	  }
	  WRITE_UINT(adjposition4,positions_low_fd);
	}
	between_counter[frame] = 0;
      }
      in_counter[frame] -= 1;
    }
#else
    if (in_counter == index1part) {
      if (
#ifdef NONMODULAR
	  between_counter >= index1interval
#else
	  (chrpos-index1part+1U) % index1interval == 0
#endif
	  ) {
	masked = oligo & mask;
	if (coord_values_8p == true) {
	  positions_move_absolute_1(positions_high_fd,offsets[masked]);
	  positions_move_absolute_4(positions_low_fd,offsets[masked]);
	  offsets[masked] += 1;
	  adjposition8 = position - index1part + 1;
	  WRITE_CHAR((unsigned char) (adjposition8 >> POSITIONS8_HIGH_SHIFT),positions_high_fd);
	  WRITE_UINT((UINT4) (adjposition8 & POSITIONS8_LOW_MASK),positions_low_fd);
	} else {
	  positions_move_absolute_4(positions_low_fd,offsets[masked]);
	  offsets[masked] += 1;
	  adjposition4 = (UINT4) (position - index1part + 1);
	  WRITE_UINT(adjposition4,positions_low_fd);
	}
	between_counter = 0;
      }
      in_counter--;
    }
#endif
    
    chrpos++;			/* Needs to go here, before we reset chrpos to 0 */
    if (position >= next_chrbound) {
#ifdef PMAP
      high = low = carry = 0U;
      in_counter[0] = in_counter[1] = in_counter[2] = 0;
#else
      oligo = 0U; in_counter = 0;
#endif
      chrpos = 0U;
      chrnum++;
      while (chrnum <= nchrs && (next_chrbound = Univ_IIT_next_chrbound(chromosome_iit,chrnum,circular_typeint)) < position) {
	chrnum++;
      }
    }
    position++;
  }

#ifdef ADDSENTINEL
  if (coord_values_8p == true) {
    for (oligoi = 0; oligoi < oligospace; oligoi++) {
      positions_move_absolute_1(positions_high_fd,offsets[oligoi]);
      positions_move_absolute_4(positions_low_fd,offsets[oligoi]);
      WRITE_CHAR((unsigned char) (-1UL >> POSITIONS8_HIGH_SHIFT),positions_high_fd);
      WRITE_UINT((UINT4) (-1UL & POSITIONS8_LOW_SHIFT),positions_low_fd);
    }
  } else {
    for (oligoi = 0; oligoi < oligospace; oligoi++) {
      positions_move_absolute_4(positions_low_fd,offsets[oligoi]);
      WRITE_UINT(-1U,positions_low_fd);
    }
  }
#endif

  return;
}

/* Requires sufficient memory to hold all positions */
static void
compute_positions_in_memory (UINT4 *positions4, unsigned char *positions8_high, UINT4 *positions8_low,
			     Positionsptr_T *offsets, FILE *sequence_fp, Univ_IIT_T chromosome_iit,
#ifdef PMAP
			     Width_T index1part_aa, bool watsonp,
#else
			     Width_T index1part,
#endif
			     Width_T index1interval, bool genome_lc_p, char *fileroot, bool mask_lowercase_p,
			     bool coord_values_8p) {
  char *uppercaseCode;
  Univcoord_T position = 0, next_chrbound;
  Chrpos_T chrpos = 0U;
  char *comma;
  int c, nchrs, chrnum;
  UINT8 adjposition8;
#ifdef PMAP
  int frame = -1, between_counter[3], in_counter[3];
  Storedoligomer_T high = 0U, low = 0U, carry;
  Storedoligomer_T aaindex;
  Width_T index1part_nt = 3*index1part_aa;
  debug1(char *aa);
#else
  int between_counter = 0, in_counter = 0;
  Storedoligomer_T oligo = 0U, masked, mask;
  debug1(char *nt);
#endif
  int circular_typeint;

#ifdef ADDSENTINEL
  Oligospace_T oligospace, oligoi;
  oligospace = power(4,index1part);
#endif

  if (mask_lowercase_p == false) {
    uppercaseCode = UPPERCASE_U2T; /* We are reading DNA sequence */
  } else {
    uppercaseCode = NO_UPPERCASE;
  }

#ifdef PMAP
  between_counter[0] = between_counter[1] = between_counter[2] = 0;
  in_counter[0] = in_counter[1] = in_counter[2] = 0;
#else
  mask = ~(~0UL << 2*index1part);
#endif

  /* Handle reference strain */
  circular_typeint = Univ_IIT_typeint(chromosome_iit,"circular");
  chrnum = 1;
  nchrs = Univ_IIT_total_nintervals(chromosome_iit);
  next_chrbound = Univ_IIT_next_chrbound(chromosome_iit,chrnum,circular_typeint);

  while ((c = Compress_get_char(sequence_fp,position,genome_lc_p)) != EOF) {
#ifdef PMAP
    if (++frame == 3) {
      frame = 0;
    }
    between_counter[frame] += 1;
    in_counter[frame] += 1;
#else
    between_counter++;
    in_counter++;
#endif

    if (position % MONITOR_INTERVAL == 0) {
      comma = Genomicpos_commafmt(position);
#ifdef PMAP
      fprintf(stderr,"Indexing positions of oligomers in genome %s (%d aa every %d aa), position %s",
	      fileroot,index1part_aa,index1interval,comma);
#else
      fprintf(stderr,"Indexing positions of oligomers in genome %s (%d bp every %d bp), position %s",
	      fileroot,index1part,index1interval,comma);
#endif
      FREE(comma);
#ifdef PMAP
      if (watsonp == true) {
	fprintf(stderr," (fwd)");
      } else {
	fprintf(stderr," (rev)");
      }
#endif
      fprintf(stderr,"\n");
    }

#ifdef PMAP
    carry = (low >> 30);
    switch (uppercaseCode[c]) {
    case 'A': low = (low << 2); break;
    case 'C': low = (low << 2) | 1U; break;
    case 'G': low = (low << 2) | 2U; break;
    case 'T': low = (low << 2) | 3U; break;
    case 'X': case 'N':
      high = low = carry = 0U;
      in_counter[0] = in_counter[1] = in_counter[2] = 0;
      break;
    default: 
      if (genome_lc_p == true) {
	high = low = carry = 0U;
	in_counter[0] = in_counter[1] = in_counter[2] = 0;
      } else {
	fprintf(stderr,"Bad character %c at position %u\n",c,position);
	abort();
      }
    }
    high = (high << 2) | carry; 

    debug(printf("frame=%d char=%c bc=%d ic=%d high=%08X low=%08X\n",
		 frame,c,between_counter[frame],in_counter[frame],high,low));
#else
    switch (uppercaseCode[c]) {
    case 'A': oligo = (oligo << 2); break;
    case 'C': oligo = (oligo << 2) | 1U; break;
    case 'G': oligo = (oligo << 2) | 2U; break;
    case 'T': oligo = (oligo << 2) | 3U; break;
    case 'X': case 'N': oligo = 0U; in_counter = 0; break;
    default: 
      if (genome_lc_p == true) {
	oligo = 0U; in_counter = 0;
      } else {
	fprintf(stderr,"Bad character %c at position %lu\n",c,position);
	abort();
      }
    }

    debug(printf("char=%c bc=%d ic=%d oligo=%08X\n",
		 c,between_counter,in_counter,oligo));
#endif

    
#ifdef PMAP
    if (in_counter[frame] > 0) {
      if (watsonp == true) {
	if (Alphabet_get_codon_fwd(low) == AA_STOP) {
	  in_counter[frame] = 0; 
	}
      } else {
	if (Alphabet_get_codon_rev(low) == AA_STOP) {
	  in_counter[frame] = 0; 
	}
      }
    }
    if (in_counter[frame] == index1part_aa + 1) {
      if (between_counter[frame] >= index1interval) {
	aaindex = Alphabet_get_aa_index(high,low,watsonp,index1part_nt);
	if (coord_values_8p == true) {
	  if (watsonp == true) {
	    adjposition8 = position - index1part_nt + 1;
	    positions8_high[offsets[aaindex]/*++*/] = (unsigned char) (adjposition8 >> POSITIONS8_HIGH_SHIFT);
	    positions8_low[offsets[aaindex]++] = (UINT4) (adjposition8 & POSITIONS8_LOW_MASK);
	    debug1(adjposition = position-index1part_nt+1UL);
	  } else {
	    positions8_high[offsets[aaindex]/*++*/] = (unsigned char) (position >> POSITIONS8_HIGH_SHIFT);
	    positions8_low[offsets[aaindex]++] = (UINT4) (position & POSITIONS8_LOW_MASK);
	    debug1(adjposition = position);
	  }
	} else {
	  if (watsonp == true) {
	    positions4[offsets[aaindex]++] = (UINT4) (position - index1part_nt + 1);
	    debug1(adjposition = position - index1part_nt + 1);
	  } else {
	    positions4[offsets[aaindex]++] = (UINT4) position;
	    debug1(adjposition = position);
	  }
	}
	debug1(
	       aa = aaindex_aa(aaindex);
	       printf("Storing %s (%u) at %u\n",aa,aaindex,adjposition);
	       FREE(aa);
	       );
	between_counter[frame] = 0;
      }
      in_counter[frame] -= 1;
    }
#else
    if (in_counter == index1part) {
      if (
#ifdef NONMODULAR
	  between_counter >= index1interval
#else
	  (chrpos-index1part+1U) % index1interval == 0
#endif
	  ) {
	masked = oligo & mask;
#if 0
	if (offsets[masked] >= totalcounts) {
	  fprintf(stderr,"masked = %u, offsets[masked] = %u >= totalcounts %u\n",
		  masked,offsets[masked],totalcounts);
	  exit(9);
	}
#endif
	if (coord_values_8p == true) {
	  adjposition8 = position - index1part + 1;
	  positions8_high[offsets[masked]/*++*/] = (unsigned char) (adjposition8 >> POSITIONS8_HIGH_SHIFT);
	  positions8_low[offsets[masked]++] = (UINT4) (adjposition8 & POSITIONS8_LOW_MASK);
	} else {
	  positions4[offsets[masked]++] = (UINT4) (position - index1part + 1);
	}
	debug1(nt = shortoligo_nt(masked,index1part);
	       printf("Storing %s at %lu, chrpos %u\n",nt,position-index1part+1U,chrpos-index1part+1U);
	       FREE(nt));
	between_counter = 0;
      }
      in_counter--;
    }
#endif
    
    chrpos++;			/* Needs to go here, before we reset chrpos to 0 */
    if (position >= next_chrbound) {
      debug1(printf("Skipping because position %u is at chrbound\n",position));
#ifdef PMAP
      high = low = carry = 0U;
      in_counter[0] = in_counter[1] = in_counter[2] = 0;
#else
      oligo = 0U; in_counter = 0;
#endif
      chrpos = 0U;
      chrnum++;
      while (chrnum <= nchrs && (next_chrbound = Univ_IIT_next_chrbound(chromosome_iit,chrnum,circular_typeint)) < position) {
	chrnum++;
      }
    }
    position++;
  }

#ifdef ADDSENTINEL
  if (coord_values_8p == true) {
    for (oligoi = 0; oligoi < oligospace; oligoi++) {
      positions8_high[offsets[oligoi]] = (unsigned char) (-1UL >> POSITIONS8_HIGH_SHIFT);
      positions8_low[offsets[oligoi]] = (UINT4) (-1UL & POSITIONS8_LOW_MASK);
    }
  } else {
    for (oligoi = 0; oligoi < oligospace; oligoi++) {
      positions4[offsets[oligoi]] = (UINT4) -1U;
    }
  }
  /*
  fprintf(stderr,"Put a sentinel at %u\n",offsets[0]);
  fprintf(stderr,"Put a sentinel at %u\n",offsets[oligospace-2]);
  fprintf(stderr,"Put a sentinel at %u\n",offsets[oligospace-1]);
  */
#endif

  return;
}


#ifdef HAVE_64_BIT
/* Requires sufficient memory to hold all positions */
/* Main difference is that we use Hugepositionsptr_T * for offsets.
   Also, huge genome implies large genome, so need only positions8. */
static void
compute_positions_in_memory_huge (UINT4 *positions4, unsigned char *positions8_high, UINT4 *positions8_low,
				  Hugepositionsptr_T *offsets, FILE *sequence_fp, Univ_IIT_T chromosome_iit,
#ifdef PMAP
				  Width_T index1part_aa, bool watsonp,
#else
				  Width_T index1part,
#endif
				  Width_T index1interval, bool genome_lc_p, char *fileroot, bool mask_lowercase_p,
				  bool coord_values_8p) {
  char *uppercaseCode;
  Univcoord_T position = 0, next_chrbound;
  Chrpos_T chrpos = 0U;
  char *comma;
  int c, nchrs, chrnum;
  UINT8 adjposition8;
#ifdef PMAP
  int frame = -1, between_counter[3], in_counter[3];
  Storedoligomer_T high = 0U, low = 0U, carry;
  Storedoligomer_T aaindex;
  Width_T index1part_nt = 3*index1part_aa;
  debug1(char *aa);
#else
  int between_counter = 0, in_counter = 0;
  Storedoligomer_T oligo = 0U, masked, mask;
  debug1(char *nt);
#endif
  int circular_typeint;

#ifdef ADDSENTINEL
  Oligospace_T oligospace, oligoi;
  oligospace = power(4,index1part);
#endif

  if (mask_lowercase_p == false) {
    uppercaseCode = UPPERCASE_U2T; /* We are reading DNA sequence */
  } else {
    uppercaseCode = NO_UPPERCASE;
  }

#ifdef PMAP
  between_counter[0] = between_counter[1] = between_counter[2] = 0;
  in_counter[0] = in_counter[1] = in_counter[2] = 0;
#else
  mask = ~(~0UL << 2*index1part);
#endif

  /* Handle reference strain */
  circular_typeint = Univ_IIT_typeint(chromosome_iit,"circular");
  chrnum = 1;
  nchrs = Univ_IIT_total_nintervals(chromosome_iit);
  next_chrbound = Univ_IIT_next_chrbound(chromosome_iit,chrnum,circular_typeint);

  while ((c = Compress_get_char(sequence_fp,position,genome_lc_p)) != EOF) {
#ifdef PMAP
    if (++frame == 3) {
      frame = 0;
    }
    between_counter[frame] += 1;
    in_counter[frame] += 1;
#else
    between_counter++;
    in_counter++;
#endif

    if (position % MONITOR_INTERVAL == 0) {
      comma = Genomicpos_commafmt(position);
#ifdef PMAP
      fprintf(stderr,"Indexing positions of oligomers in genome %s (%d aa every %d aa), position %s",
	      fileroot,index1part_aa,index1interval,comma);
#else
      fprintf(stderr,"Indexing positions of oligomers in genome %s (%d bp every %d bp), position %s",
	      fileroot,index1part,index1interval,comma);
#endif
      FREE(comma);
#ifdef PMAP
      if (watsonp == true) {
	fprintf(stderr," (fwd)");
      } else {
	fprintf(stderr," (rev)");
      }
#endif
      fprintf(stderr,"\n");
    }

#ifdef PMAP
    carry = (low >> 30);
    switch (uppercaseCode[c]) {
    case 'A': low = (low << 2); break;
    case 'C': low = (low << 2) | 1U; break;
    case 'G': low = (low << 2) | 2U; break;
    case 'T': low = (low << 2) | 3U; break;
    case 'X': case 'N':
      high = low = carry = 0U;
      in_counter[0] = in_counter[1] = in_counter[2] = 0;
      break;
    default: 
      if (genome_lc_p == true) {
	high = low = carry = 0U;
	in_counter[0] = in_counter[1] = in_counter[2] = 0;
      } else {
	fprintf(stderr,"Bad character %c at position %u\n",c,position);
	abort();
      }
    }
    high = (high << 2) | carry; 

    debug(printf("frame=%d char=%c bc=%d ic=%d high=%08X low=%08X\n",
		 frame,c,between_counter[frame],in_counter[frame],high,low));
#else
    switch (uppercaseCode[c]) {
    case 'A': oligo = (oligo << 2); break;
    case 'C': oligo = (oligo << 2) | 1U; break;
    case 'G': oligo = (oligo << 2) | 2U; break;
    case 'T': oligo = (oligo << 2) | 3U; break;
    case 'X': case 'N': oligo = 0U; in_counter = 0; break;
    default: 
      if (genome_lc_p == true) {
	oligo = 0U; in_counter = 0;
      } else {
	fprintf(stderr,"Bad character %c at position %lu\n",c,position);
	abort();
      }
    }

    debug(printf("char=%c bc=%d ic=%d oligo=%08X\n",
		 c,between_counter,in_counter,oligo));
#endif

    
#ifdef PMAP
    if (in_counter[frame] > 0) {
      if (watsonp == true) {
	if (Alphabet_get_codon_fwd(low) == AA_STOP) {
	  in_counter[frame] = 0; 
	}
      } else {
	if (Alphabet_get_codon_rev(low) == AA_STOP) {
	  in_counter[frame] = 0; 
	}
      }
    }
    if (in_counter[frame] == index1part_aa + 1) {
      if (between_counter[frame] >= index1interval) {
	aaindex = Alphabet_get_aa_index(high,low,watsonp,index1part_nt);
	if (coord_values_8p == true) {
	  if (watsonp == true) {
	    adjposition8 = position - index1part_nt + 1;
	    positions8_high[offsets[aaindex]/*++*/] = (unsigned char) (adjposition8 >> POSITIONS8_HIGH_SHIFT);
	    positions8_low[offsets[aaindex]++] = (UINT4) (adjposition8 & POSITIONS8_LOW_MASK);
	    debug1(adjposition = position-index1part_nt+1);
	  } else {
	    positions8_high[offsets[aaindex]/*++*/] = (unsigned char) (position >> POSITIONS8_HIGH_SHIFT);
	    positions8_low[offsets[aaindex]++] = (UINT4) (position & POSITIONS8_LOW_MASK);
	    debug1(adjposition = position);
	  }
	} else {
	  if (watsonp == true) {
	    positions4[offsets[aaindex]++] = (UINT4) (position-index1part_nt+1U);
	    debug1(adjposition = position-index1part_nt+1U);
	  } else {
	    positions4[offsets[aaindex]++] = (UINT4) position;
	    debug1(adjposition = position);
	  }
	}
	debug1(
	       aa = aaindex_aa(aaindex);
	       printf("Storing %s (%u) at %u\n",aa,aaindex,adjposition);
	       FREE(aa);
	       );
	between_counter[frame] = 0;
      }
      in_counter[frame] -= 1;
    }
#else
    if (in_counter == index1part) {
      if (
#ifdef NONMODULAR
	  between_counter >= index1interval
#else
	  (chrpos-index1part+1U) % index1interval == 0
#endif
	  ) {
	masked = oligo & mask;
#if 0
	if (offsets[masked] >= totalcounts) {
	  fprintf(stderr,"masked = %u, offsets[masked] = %u >= totalcounts %u\n",
		  masked,offsets[masked],totalcounts);
	  exit(9);
	}
#endif
	if (coord_values_8p == true) {
	  adjposition8 = position - index1part + 1;
	  positions8_high[offsets[masked]/*++*/] = (unsigned char) (adjposition8 >> POSITIONS8_HIGH_SHIFT);
	  positions8_low[offsets[masked]++] = (UINT4) (adjposition8 & POSITIONS8_LOW_MASK);
	} else {
	  positions4[offsets[masked]++] = (UINT4) (position - index1part + 1);
	}
	debug1(nt = shortoligo_nt(masked,index1part);
	       printf("Storing %s at %lu, chrpos %u\n",nt,position-index1part+1U,chrpos-index1part+1U);
	       FREE(nt));
	between_counter = 0;
      }
      in_counter--;
    }
#endif
    
    chrpos++;			/* Needs to go here, before we reset chrpos to 0 */
    if (position >= next_chrbound) {
      debug1(printf("Skipping because position %u is at chrbound\n",position));
#ifdef PMAP
      high = low = carry = 0U;
      in_counter[0] = in_counter[1] = in_counter[2] = 0;
#else
      oligo = 0U; in_counter = 0;
#endif
      chrpos = 0U;
      chrnum++;
      while (chrnum <= nchrs && (next_chrbound = Univ_IIT_next_chrbound(chromosome_iit,chrnum,circular_typeint)) < position) {
	chrnum++;
      }
    }
    position++;
  }

#ifdef ADDSENTINEL
  if (coord_values_8p == true) {
    for (oligoi = 0; oligoi < oligospace; oligoi++) {
      positions8_high[offsets[oligoi]] = (unsigned char) (-1UL >> POSITIONS8_HIGH_SHIFT);
      positions8_low[offsets[oligoi]] = (UINT4) (-1UL & POSITIONS8_LOW_MASK);
    }
  } else {
    for (oligoi = 0; oligoi < oligospace; oligoi++) {
      positions4[offsets[oligoi]] = (UINT4) -1U;
    }
  }
  /*
  fprintf(stderr,"Put a sentinel at %u\n",offsets[0]);
  fprintf(stderr,"Put a sentinel at %u\n",offsets[oligospace-2]);
  fprintf(stderr,"Put a sentinel at %u\n",offsets[oligospace-1]);
  */
#endif

  return;
}
#endif



#define WRITE_CHUNK 1000000

void
Indexdb_write_positions (char *positionsfile_high, char *positionsfile_low, char *pointersfile, char *offsetsfile,
			 FILE *sequence_fp, Univ_IIT_T chromosome_iit, int offsetscomp_basesize,
#ifdef PMAP
			 Alphabet_T alphabet, int index1part_aa, bool watsonp,
#else
			 int index1part,
#endif
			 int index1interval, bool genome_lc_p, bool writefilep,
			 char *fileroot, bool mask_lowercase_p, int compression_type,
			 bool coord_values_8p) {
  FILE *positions_high_fp, *positions_low_fp; /* For building positions in memory */
  int positions_high_fd, positions_low_fd; /* For building positions in file */
  Positionsptr_T *offsets = NULL, totalcounts, count;
  UINT4 *positions4;
  unsigned char *positions8_high;
  UINT4 *positions8_low;
  Oligospace_T oligospace;
  off_t filesize;
  size_t offsetscomp_len;
  double seconds;
#ifdef PMAP
  int alphabet_size;
  int offsetscomp_fd;
#endif



#ifdef PMAP
  alphabet_size = Alphabet_get_size(alphabet);
  if (compression_type == BITPACK64_COMPRESSION) {
    offsets = Indexdb_offsets_from_bitpack(pointersfile,offsetsfile,alphabet_size,index1part_aa);
  } else if (compression_type == GAMMA_COMPRESSION) {
    offsets = Indexdb_offsets_from_gammas(pointersfile,offsetsfile,offsetscomp_basesize,
					  alphabet_size,index1part_aa);
  } else {
    offsets = (UINT4 *) Access_allocated(&offsetscomp_len,&seconds,offsetsfile,sizeof(UINT4));
  }
  oligospace = power(alphabet_size,index1part_aa);
#else
  if (compression_type == BITPACK64_COMPRESSION) {
    offsets = Indexdb_offsets_from_bitpack(pointersfile,offsetsfile,offsetscomp_basesize,index1part);
  } else if (compression_type == GAMMA_COMPRESSION) {
    offsets = Indexdb_offsets_from_gammas(pointersfile,offsetsfile,offsetscomp_basesize,index1part);
  } else {
    offsets = (UINT4 *) Access_allocated(&offsetscomp_len,&seconds,offsetsfile,sizeof(UINT4));
  }
  oligospace = power(4,index1part);
#endif
  totalcounts = offsets[oligospace];
  if (totalcounts == 0) {
    fprintf(stderr,"Something is wrong with the offsets file.  Total counts is zero.\n");
    exit(9);
  }

  if (writefilep == true) {
    fprintf(stderr,"User requested build of positions in file\n");
    positions_high_fd = Access_fileio_rw(positionsfile_high);
    positions_low_fd = Access_fileio_rw(positionsfile_low);
#ifdef PMAP
    compute_positions_in_file(positions_high_fd,positions_low_fd,offsets,sequence_fp,chromosome_iit,
			      index1part_aa,watsonp,index1interval,genome_lc_p,fileroot,mask_lowercase_p,
			      coord_values_8p);
#else
    compute_positions_in_file(positions_high_fd,positions_low_fd,offsets,sequence_fp,chromosome_iit,
			      index1part,index1interval,genome_lc_p,fileroot,mask_lowercase_p,
			      coord_values_8p);
#endif
    close(positions_high_fd);
    close(positions_low_fd);

  } else if (coord_values_8p == true) {
    fprintf(stderr,"Trying to allocate %u*(%d+%d) bytes of memory for positions...",totalcounts,(int) sizeof(unsigned char),(int) sizeof(UINT4));
    positions8_high = (unsigned char *) CALLOC_NO_EXCEPTION(totalcounts,sizeof(unsigned char));
    positions8_low = (UINT4 *) CALLOC_NO_EXCEPTION(totalcounts,sizeof(UINT4));
    if (positions8_high == NULL || positions8_low == NULL) {
      fprintf(stderr,"failed.  Building positions in file.\n");
      positions_high_fd = Access_fileio_rw(positionsfile_high);
      positions_low_fd = Access_fileio_rw(positionsfile_low);
#ifdef PMAP
      compute_positions_in_file(positions_high_fd,positions_low_fd,offsets,sequence_fp,chromosome_iit,
				index1part_aa,watsonp,index1interval,genome_lc_p,fileroot,
				mask_lowercase_p,/*coord_values_8p*/true);
#else
      compute_positions_in_file(positions_high_fd,positions_low_fd,offsets,sequence_fp,chromosome_iit,
				index1part,index1interval,genome_lc_p,fileroot,
				mask_lowercase_p,/*coord_values_8p*/true);
#endif
      close(positions_high_fd);
      close(positions_low_fd);

      if ((filesize = Access_filesize(positionsfile_high)) != totalcounts * (off_t) sizeof(unsigned char)) {
	fprintf(stderr,"Error after file-based build: expected file size for %s is %lu, but observed only %lu.  Please notify twu@gene.com of this error.\n",
		positionsfile_high,totalcounts*sizeof(unsigned char),filesize);
	abort();
      } else if ((filesize = Access_filesize(positionsfile_low)) != totalcounts * (off_t) sizeof(UINT4)) {
	fprintf(stderr,"Error after file-based build: expected file size for %s is %lu, but observed only %lu.  Please notify twu@gene.com of this error.\n",
		positionsfile_low,totalcounts*sizeof(UINT4),filesize);
	abort();
      }

    } else {
      fprintf(stderr,"succeeded.  Building positions in memory.\n");
      if ((positions_high_fp = FOPEN_WRITE_BINARY(positionsfile_high)) == NULL) {
	fprintf(stderr,"Can't open file %s\n",positionsfile_high);
	exit(9);
      } else if ((positions_low_fp = FOPEN_WRITE_BINARY(positionsfile_low)) == NULL) {
	fprintf(stderr,"Can't open file %s\n",positionsfile_low);
	exit(9);
      }
#ifdef PMAP
      compute_positions_in_memory(/*positions4*/NULL,positions8_high,positions8_low,offsets,
				  sequence_fp,chromosome_iit,index1part_aa,watsonp,
				  index1interval,genome_lc_p,fileroot,
				  mask_lowercase_p,/*coord_values_8p*/true);
#else
      compute_positions_in_memory(/*positions4*/NULL,positions8_high,positions8_low,offsets,
				  sequence_fp,chromosome_iit,index1part,
				  index1interval,genome_lc_p,fileroot,
				  mask_lowercase_p,/*coord_values_8p*/true);
#endif
      fprintf(stderr,"Writing %u genomic positions to files %s and %s ...\n",
	      totalcounts,positionsfile_high,positionsfile_low);
      FWRITE_CHARS(positions8_high,totalcounts,positions_high_fp);
      FWRITE_UINTS(positions8_low,totalcounts,positions_low_fp);

      fclose(positions_high_fp);
      fclose(positions_low_fp);

      if ((filesize = Access_filesize(positionsfile_high)) != totalcounts * (off_t) sizeof(unsigned char)) {
	fprintf(stderr,"Error: expected file size for %s is %lu, but observed only %lu.  Trying now to write with smaller chunks.\n",
		positionsfile_high,totalcounts*sizeof(unsigned char),filesize);
	if ((positions_high_fp = FOPEN_WRITE_BINARY(positionsfile_high)) == NULL) {
	  fprintf(stderr,"Can't open file %s\n",positionsfile_high);
	  exit(9);
	}
	for (count = 0; count + WRITE_CHUNK < totalcounts; count += WRITE_CHUNK) {
	  FWRITE_CHARS(&(positions8_high[count]),count,positions_high_fp);
	}
	if (count < totalcounts) {
	  FWRITE_CHARS(&(positions8_high[count]),totalcounts - count,positions_high_fp);
	}
	fclose(positions_high_fp);

	if ((filesize = Access_filesize(positionsfile_high)) != totalcounts * (off_t) sizeof(unsigned char)) {
	  fprintf(stderr,"Error persists: expected file size for %s is %lu, but observed only %lu.  Please notify twu@gene.com of this error.\n",
		  positionsfile_high,totalcounts*sizeof(unsigned char),filesize);
	  abort();
	}
      }

      if ((filesize = Access_filesize(positionsfile_low)) != totalcounts * (off_t) sizeof(UINT4)) {
	fprintf(stderr,"Error: expected file size for %s is %lu, but observed only %lu.  Trying now to write with smaller chunks.\n",
		positionsfile_low,totalcounts*sizeof(UINT4),filesize);
	if ((positions_low_fp = FOPEN_WRITE_BINARY(positionsfile_low)) == NULL) {
	  fprintf(stderr,"Can't open file %s\n",positionsfile_low);
	  exit(9);
	}
	for (count = 0; count + WRITE_CHUNK < totalcounts; count += WRITE_CHUNK) {
	  FWRITE_UINTS(&(positions8_low[count]),count,positions_low_fp);
	}
	if (count < totalcounts) {
	  FWRITE_UINTS(&(positions8_low[count]),totalcounts - count,positions_low_fp);
	}
	fclose(positions_low_fp);

	if ((filesize = Access_filesize(positionsfile_low)) != totalcounts * (off_t) sizeof(UINT4)) {
	  fprintf(stderr,"Error persists: expected file size for %s is %lu, but observed only %lu.  Please notify twu@gene.com of this error.\n",
		  positionsfile_low,totalcounts*sizeof(UINT4),filesize);
	  abort();
	}
      }

      FREE(positions8_high);
      FREE(positions8_low);
    }

  } else {
    fprintf(stderr,"Trying to allocate %u*%d bytes of memory for positions...",totalcounts,(int) sizeof(UINT4));
    positions4 = (UINT4 *) CALLOC_NO_EXCEPTION(totalcounts,sizeof(UINT4));
    if (positions4 == NULL) {
      fprintf(stderr,"failed.  Building positions in file.\n");
      positions_low_fd = Access_fileio_rw(positionsfile_low);
#ifdef PMAP
      compute_positions_in_file(/*positions_high_fd*/0,positions_low_fd,offsets,sequence_fp,chromosome_iit,
				index1part_aa,watsonp,index1interval,genome_lc_p,fileroot,mask_lowercase_p,
				/*coord_values_8*/false);
#else
      compute_positions_in_file(/*positions_high_fd*/0,positions_low_fd,offsets,sequence_fp,chromosome_iit,
				index1part,index1interval,genome_lc_p,fileroot,mask_lowercase_p,
				/*coord_values_8p*/false);
#endif
      close(positions_low_fd);

      if ((filesize = Access_filesize(positionsfile_low)) != totalcounts * (off_t) sizeof(UINT4)) {
	fprintf(stderr,"Error after file-based build: expected file size for %s is %lu, but observed only %lu.  Please notify twu@gene.com of this error.\n",
		positionsfile_low,totalcounts*sizeof(UINT4),filesize);
	abort();
      }

    } else {
      fprintf(stderr,"succeeded.  Building positions in memory.\n");
      if ((positions_low_fp = FOPEN_WRITE_BINARY(positionsfile_low)) == NULL) {
	fprintf(stderr,"Can't open file %s\n",positionsfile_low);
	exit(9);
      }
#ifdef PMAP
      compute_positions_in_memory(positions4,/*positions8_high*/NULL,/*positions8_low*/NULL,
				  offsets,sequence_fp,chromosome_iit,
				  index1part_aa,watsonp,index1interval,genome_lc_p,fileroot,
				  mask_lowercase_p,/*coord_values_8p*/false);
#else
      compute_positions_in_memory(positions4,/*positions8_high*/NULL,/*positions8_low*/NULL,
				  offsets,sequence_fp,chromosome_iit,
				  index1part,index1interval,genome_lc_p,fileroot,
				  mask_lowercase_p,/*coord_values_8p*/false);
#endif
      fprintf(stderr,"Writing %u genomic positions to file %s ...\n",
	      totalcounts,positionsfile_low);
      FWRITE_UINTS(positions4,totalcounts,positions_low_fp);

      fclose(positions_low_fp);

      if ((filesize = Access_filesize(positionsfile_low)) != totalcounts * (off_t) sizeof(UINT4)) {
	fprintf(stderr,"Error: expected file size for %s is %lu, but observed only %lu.  Trying now to write with smaller chunks.\n",
		positionsfile_low,totalcounts*sizeof(UINT4),filesize);
	if ((positions_low_fp = FOPEN_WRITE_BINARY(positionsfile_low)) == NULL) {
	  fprintf(stderr,"Can't open file %s\n",positionsfile_low);
	  exit(9);
	}
	for (count = 0; count + WRITE_CHUNK < totalcounts; count += WRITE_CHUNK) {
	  FWRITE_UINTS(&(positions4[count]),count,positions_low_fp);
	}
	if (count < totalcounts) {
	  FWRITE_UINTS(&(positions4[count]),totalcounts - count,positions_low_fp);
	}
	fclose(positions_low_fp);

	if ((filesize = Access_filesize(positionsfile_low)) != totalcounts * (off_t) sizeof(UINT4)) {
	  fprintf(stderr,"Error persists: expected file size for %s is %lu, but observed only %lu.  Please notify twu@gene.com of this error.\n",
		  positionsfile_low,totalcounts*sizeof(UINT4),filesize);
	  abort();
	}
      }

      FREE(positions4);
    }
  }

  FREE(offsets);

  return;
}


#ifdef HAVE_64_BIT
void
Indexdb_write_positions_huge (char *positionsfile_high, char *positionsfile_low,
			      char *pagesfile, char *pointersfile, char *offsetsfile,
			      FILE *sequence_fp, Univ_IIT_T chromosome_iit, int offsetscomp_basesize,
#ifdef PMAP
			      Alphabet_T alphabet, int index1part_aa, bool watsonp,
#else
			      int index1part,
#endif
			      int index1interval, bool genome_lc_p, bool writefilep,
			      char *fileroot, bool mask_lowercase_p, int compression_type,
			      bool coord_values_8p) {
  FILE *positions_high_fp, *positions_low_fp; /* For building positions in memory */
  Hugepositionsptr_T *offsets = NULL, totalcounts, count;
  UINT4 *positions4;
  unsigned char *positions8_high;
  UINT4 *positions8_low;
  Oligospace_T oligospace;
  off_t filesize;
#ifdef PMAP
  int alphabet_size;
#endif


#ifdef PMAP
  alphabet_size = Alphabet_get_size(alphabet);
  offsets = Indexdb_offsets_from_gammas(pointersfile,offsetsfile,offsetscomp_basesize,
					alphabet_size,index1part_aa);
  oligospace = power(alphabet_size,index1part_aa);
#else
  if (compression_type == BITPACK64_COMPRESSION) {
    offsets = Indexdb_offsets_from_bitpack_huge(pagesfile,pointersfile,offsetsfile,offsetscomp_basesize,index1part);
  } else if (compression_type == GAMMA_COMPRESSION) {
    fprintf(stderr,"Not supporting gamma compression on huge genomes\n");
    /* offsets = Indexdb_offsets_from_gammas(pointersfile,offsetsfile,offsetscomp_basesize,index1part); */
    return;
  } else {
    fprintf(stderr,"Not supporting non-compression on huge genomes\n");
    /* offsets = Indexdb_offsets_from_gammas(pointersfile,offsetsfile,offsetscomp_basesize,index1part); */
    return;
  }
  oligospace = power(4,index1part);
#endif
  totalcounts = offsets[oligospace];
  if (totalcounts == 0) {
    fprintf(stderr,"Something is wrong with the offsets file.  Total counts is zero.\n");
    exit(9);
  }

  if (writefilep == true) {
    fprintf(stderr,"User requested build of positions in file, but not supported.\n");
    abort();

  } else if (coord_values_8p == true) {
    fprintf(stderr,"Trying to allocate %lu*(%d+%d) bytes of memory for positions...",
	    totalcounts,(int) sizeof(unsigned char),(int) sizeof(UINT8));
    positions8_high = (unsigned char *) CALLOC_NO_EXCEPTION(totalcounts,sizeof(unsigned char));
    positions8_low = (UINT4 *) CALLOC_NO_EXCEPTION(totalcounts,sizeof(UINT4));
    if (positions8_high == NULL || positions8_low == NULL) {
      fprintf(stderr,"failed.  Not able to proceed.\n");
      abort();

    } else {
      fprintf(stderr,"succeeded.  Building positions in memory.\n");
      if ((positions_high_fp = FOPEN_WRITE_BINARY(positionsfile_high)) == NULL) {
	fprintf(stderr,"Can't open file %s\n",positionsfile_high);
	exit(9);
      } else if ((positions_low_fp = FOPEN_WRITE_BINARY(positionsfile_low)) == NULL) {
	fprintf(stderr,"Can't open file %s\n",positionsfile_low);
	exit(9);
      }
#ifdef PMAP
      compute_positions_in_memory_huge(/*positions4*/NULL,positions8_high,positions8_low,offsets,
				       sequence_fp,chromosome_iit,index1part_aa,watsonp,
				       index1interval,genome_lc_p,fileroot,
				       mask_lowercase_p,/*coord_values_8p*/true);
#else
      compute_positions_in_memory_huge(/*positions4*/NULL,positions8_high,positions8_low,offsets,
				       sequence_fp,chromosome_iit,index1part,index1interval,
				       genome_lc_p,fileroot,mask_lowercase_p,/*coord_values_8p*/true);
#endif
      fprintf(stderr,"Writing %lu genomic positions to files %s and %s...\n",
	      totalcounts,positionsfile_high,positionsfile_low);
      FWRITE_CHARS(positions8_high,totalcounts,positions_high_fp);
      FWRITE_UINTS(positions8_low,totalcounts,positions_low_fp);

      fclose(positions_high_fp);
      fclose(positions_low_fp);

      if ((filesize = Access_filesize(positionsfile_high)) != totalcounts * (off_t) sizeof(unsigned char)) {
	fprintf(stderr,"Error: expected file size for %s is %lu, but observed only %lu.  Trying now to write with smaller chunks.\n",
		positionsfile_high,totalcounts*sizeof(unsigned char),filesize);
	if ((positions_high_fp = FOPEN_WRITE_BINARY(positionsfile_high)) == NULL) {
	  fprintf(stderr,"Can't open file %s\n",positionsfile_high);
	  exit(9);
	}
	for (count = 0; count + WRITE_CHUNK < totalcounts; count += WRITE_CHUNK) {
	  FWRITE_CHARS(&(positions8_high[count]),count,positions_high_fp);
	}
	if (count < totalcounts) {
	  FWRITE_CHARS(&(positions8_high[count]),totalcounts - count,positions_high_fp);
	}
	fclose(positions_high_fp);

	if ((filesize = Access_filesize(positionsfile_high)) != totalcounts * (off_t) sizeof(unsigned char)) {
	  fprintf(stderr,"Error persists: expected file size for %s is %lu, but observed only %lu.  Please notify twu@gene.com of this error.\n",
		  positionsfile_high,totalcounts*sizeof(unsigned char),filesize);
	  abort();
	}
      }

      if ((filesize = Access_filesize(positionsfile_low)) != totalcounts * (off_t) sizeof(UINT4)) {
	fprintf(stderr,"Error: expected file size for %s is %lu, but observed only %lu.  Trying now to write with smaller chunks.\n",
		positionsfile_low,totalcounts*sizeof(UINT4),filesize);
	if ((positions_low_fp = FOPEN_WRITE_BINARY(positionsfile_low)) == NULL) {
	  fprintf(stderr,"Can't open file %s\n",positionsfile_low);
	  exit(9);
	}
	for (count = 0; count + WRITE_CHUNK < totalcounts; count += WRITE_CHUNK) {
	  FWRITE_UINTS(&(positions8_low[count]),count,positions_low_fp);
	}
	if (count < totalcounts) {
	  FWRITE_UINTS(&(positions8_low[count]),totalcounts - count,positions_low_fp);
	}
	fclose(positions_low_fp);

	if ((filesize = Access_filesize(positionsfile_low)) != totalcounts * (off_t) sizeof(UINT4)) {
	  fprintf(stderr,"Error persists: expected file size for %s is %lu, but observed only %lu.  Please notify twu@gene.com of this error.\n",
		  positionsfile_low,totalcounts*sizeof(UINT4),filesize);
	  abort();
	}
      }

      FREE(positions8_high);
      FREE(positions8_low);
    }

  } else {
    fprintf(stderr,"Something is wrong.  We have a huge genome (offsets need 8 bytes), but positions can fit in 4 bytes\n");
    fprintf(stderr,"Please report this bug to twu@gene.com\n");
    abort();

    fprintf(stderr,"Trying to allocate %lu*%d bytes of memory for positions...",totalcounts,(int) sizeof(UINT4));
    positions4 = (UINT4 *) CALLOC_NO_EXCEPTION(totalcounts,sizeof(UINT4));
    if (positions4 == NULL) {
      fprintf(stderr,"failed.  Not able to proceed.\n");
      abort();

    } else {
      fprintf(stderr,"succeeded.  Building positions in memory.\n");
      if ((positions_low_fp = FOPEN_WRITE_BINARY(positionsfile_low)) == NULL) {
	fprintf(stderr,"Can't open file %s\n",positionsfile_low);
	exit(9);
      }
#ifdef PMAP
      compute_positions_in_memory_huge(positions4,/*positions8_high*/NULL,/*positions8_low*/NULL,offsets,
				       sequence_fp,chromosome_iit,index1part_aa,watsonp,index1interval,
				       genome_lc_p,fileroot,mask_lowercase_p,/*coord_values_8p*/false);
#else
      compute_positions_in_memory_huge(positions4,/*positions8_high*/NULL,/*positions8_low*/NULL,offsets,
				       sequence_fp,chromosome_iit,index1part,index1interval,
				       genome_lc_p,fileroot,mask_lowercase_p,/*coord_values_8p*/false);
#endif
      fprintf(stderr,"Writing %lu genomic positions to file %s ...\n",
	      totalcounts,positionsfile_low);
      FWRITE_UINTS(positions4,totalcounts,positions_low_fp);

      fclose(positions_low_fp);

      if ((filesize = Access_filesize(positionsfile_low)) != totalcounts * (off_t) sizeof(UINT4)) {
	fprintf(stderr,"Error: expected file size for %s is %lu, but observed only %lu.  Trying now to write with smaller chunks.\n",
		positionsfile_low,totalcounts*sizeof(UINT4),filesize);
	if ((positions_low_fp = FOPEN_WRITE_BINARY(positionsfile_low)) == NULL) {
	  fprintf(stderr,"Can't open file %s\n",positionsfile_low);
	  exit(9);
	}
	for (count = 0; count + WRITE_CHUNK < totalcounts; count += WRITE_CHUNK) {
	  FWRITE_UINTS(&(positions4[count]),count,positions_low_fp);
	}
	if (count < totalcounts) {
	  FWRITE_UINTS(&(positions4[count]),totalcounts - count,positions_low_fp);
	}
	fclose(positions_low_fp);

	if ((filesize = Access_filesize(positionsfile_low)) != totalcounts * (off_t) sizeof(UINT4)) {
	  fprintf(stderr,"Error persists: expected file size for %s is %lu, but observed only %lu.  Please notify twu@gene.com of this error.\n",
		  positionsfile_low,totalcounts*sizeof(UINT4),filesize);
	  abort();
	}
      }

      FREE(positions4);
    }
  }

  FREE(offsets);

  return;
}
#endif


