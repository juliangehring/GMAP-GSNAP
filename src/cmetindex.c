static char rcsid[] = "$Id: cmetindex.c 37290 2011-03-29 01:28:47Z twu $";
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

#include "cmet.h"

#include "bool.h"
#include "genomicpos.h"
#include "indexdb.h"		/* For Storedoligomer_T */
#include "indexdbdef.h"		/* For Positionsptr_T */
#include "datadir.h"
#include "getopt.h"


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


static char *user_sourcedir = NULL;
static char *user_destdir = NULL;
static char *dbroot = NULL;
static char *dbversion = NULL;

static char *snps_root = NULL;


static struct option long_options[] = {
  /* Input options */
  {"sourcedir", required_argument, 0, 'F'},	/* user_sourcedir */
  {"destdir", required_argument, 0, 'D'},	/* user_destdir */
  {"db", required_argument, 0, 'd'}, /* dbroot */
  {"usesnps", required_argument, 0, 'v'}, /* snps_root */

  /* Help options */
  {"version", no_argument, 0, 0}, /* print_program_version */
  {"help", no_argument, 0, 0}, /* print_program_usage */
  {0, 0, 0, 0}
};


static void
print_program_version () {
  fprintf(stdout,"\n");
  fprintf(stdout,"CMETINDEX: Builds GMAP index files for methylated genomes\n");
  fprintf(stdout,"Part of GMAP package, version %s\n",PACKAGE_VERSION);
  fprintf(stdout,"Default gmap directory: %s\n",GMAPDB);
  fprintf(stdout,"Thomas D. Wu, Genentech, Inc.\n");
  fprintf(stdout,"Contact: twu@gene.com\n");
  fprintf(stdout,"\n");
  return;
}

static void
print_program_usage ();


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

/*               87654321 */
#define LEFT_A 0x00000000
#define LEFT_C 0x40000000
#define LEFT_G 0x80000000
#define LEFT_T 0xC0000000

/*                      87654321 */
#define LOW_TWO_BITS  0x00000003


static Storedoligomer_T
reduce_oligo_old (Storedoligomer_T oligo, int oligosize) {
  Storedoligomer_T reduced = 0U, lowbits;
  int i;

  for (i = 0; i < oligosize; i++) {
    lowbits = oligo & LOW_TWO_BITS;
    reduced >>= 2;
    switch (lowbits) {
    case RIGHT_A: break;
    case RIGHT_C: reduced |= LEFT_T; break;
    case RIGHT_G: reduced |= LEFT_G; break;
    case RIGHT_T: reduced |= LEFT_T; break;
    }
    oligo >>= 2;
  }

  for ( ; i < 16; i++) {
    reduced >>= 2;
  }

  return reduced;
}


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


static Positionsptr_T *
compute_offsets_ct (Positionsptr_T *oldoffsets, int oligospace) {
  Positionsptr_T *offsets, *sizes;
  Storedoligomer_T reduced;
  int i;
#ifdef DEBUG
  char *nt1, *nt2;
#endif

  /* Fill with sizes */
  offsets = (Positionsptr_T *) CALLOC(oligospace+1,sizeof(Positionsptr_T));
  sizes = (Positionsptr_T *) CALLOC(oligospace,sizeof(Positionsptr_T));

  for (i = 0; i < oligospace; i++) {
#if 0
    if (reduce_oligo(i) != reduce_oligo_old(i,INDEX1PART)) {
      abort();
    }
#endif
    reduced = Cmet_reduce_ct(i);
    debug(
	  nt1 = shortoligo_nt(i,INDEX1PART);
	  nt2 = shortoligo_nt(reduced,INDEX1PART);
	  printf("For oligo %s, updating sizes for %s from %d",nt1,nt2,sizes[reduced]);
	  );
#ifdef WORDS_BIGENDIAN
    sizes[reduced] += (Bigendian_convert_uint(oldoffsets[i+1]) - Bigendian_convert_uint(oldoffsets[i]));
#else
    sizes[reduced] += (oldoffsets[i+1] - oldoffsets[i]);
#endif
    debug(
	  printf(" to %d\n",sizes[reduced]);
	  FREE(nt2);
	  FREE(nt1);
	  );
  }

  offsets[0] = 0U;
  for (i = 1; i <= oligospace; i++) {
    offsets[i] = sizes[i-1] + offsets[i-1];
    debug(if (offsets[i] != offsets[i-1]) {
	    printf("Offset for %06X: %u\n",i,offsets[i]);
	  });
  }

  FREE(sizes);
  return offsets;
}

static Positionsptr_T *
compute_offsets_ga (Positionsptr_T *oldoffsets, int oligospace) {
  Positionsptr_T *offsets, *sizes;
  Storedoligomer_T reduced;
  int i;
  char *nt1, *nt2;

  /* Fill with sizes */
  offsets = (Positionsptr_T *) CALLOC(oligospace+1,sizeof(Positionsptr_T));
  sizes = (Positionsptr_T *) CALLOC(oligospace,sizeof(Positionsptr_T));

  for (i = 0; i < oligospace; i++) {
    reduced = Cmet_reduce_ga(i);
    debug(
	  nt1 = shortoligo_nt(i,INDEX1PART);
	  nt2 = shortoligo_nt(reduced,INDEX1PART);
	  printf("For oligo %s, updating sizes for %s from %d",nt1,nt2,sizes[reduced]);
	  );
#ifdef WORDS_BIGENDIAN
    sizes[reduced] += (Bigendian_convert_uint(oldoffsets[i+1]) - Bigendian_convert_uint(oldoffsets[i]));
#else
    sizes[reduced] += (oldoffsets[i+1] - oldoffsets[i]);
#endif
    debug(
	  printf(" to %d\n",sizes[reduced]);
	  FREE(nt2);
	  FREE(nt1);
	  );
  }

  offsets[0] = 0U;
  for (i = 1; i <= oligospace; i++) {
    offsets[i] = sizes[i-1] + offsets[i-1];
    debug(if (offsets[i] != offsets[i-1]) {
	    printf("Offset for %06X: %u\n",i,offsets[i]);
	  });
  }

  FREE(sizes);
  return offsets;
}


static void
compute_positions_ct (FILE *offsets_fp, FILE *positions_fp, Positionsptr_T *offsets, Positionsptr_T *oldoffsets,
		      Genomicpos_T *oldpositions, int oligospace) {
  Genomicpos_T *positions;
  Storedoligomer_T reduced;
  int i, j;
  Positionsptr_T *pointers, preunique_totalcounts, block_start, block_end, npositions, offsets_ptr;
#ifdef DEBUG
  char *nt1, *nt2;
#endif

  preunique_totalcounts = offsets[oligospace];
  if (preunique_totalcounts == 0) {
    fprintf(stderr,"Something is wrong with the offsets.  Total counts is zero.\n");
    exit(9);
  } else {
    fprintf(stderr,"Trying to allocate %u*%d bytes of memory...",preunique_totalcounts,(int) sizeof(Genomicpos_T));
    positions = (Genomicpos_T *) CALLOC_NO_EXCEPTION(preunique_totalcounts,sizeof(Genomicpos_T));
    if (positions == NULL) {
      fprintf(stderr,"failed.  Need a computer with sufficient memory.\n");
      exit(9);
    } else {
      fprintf(stderr,"done\n");
    }
  }

  /* Copy offsets */
  pointers = (Positionsptr_T *) CALLOC(oligospace,sizeof(Positionsptr_T));
  for (i = 0; i < oligospace; i++) {
    pointers[i] = offsets[i];
  }

  fprintf(stderr,"Rearranging CT positions");
  for (i = 0; i < oligospace; i++) {
    if (i % 100000 == 0) {
      fprintf(stderr,".");
    }
    reduced = Cmet_reduce_ct(i);
#ifdef WORDS_BIGENDIAN
    for (j = Bigendian_convert_uint(oldoffsets[i]); j < Bigendian_convert_uint(oldoffsets[i+1]); j++) {
      debug(nt1 = shortoligo_nt(i,INDEX1PART);
	    nt2 = shortoligo_nt(reduced,INDEX1PART);
	    printf("Oligo %s => %s: copying position %u to location %u\n",
		   nt1,nt2,oldpositions[j],pointers[i]);
	    FREE(nt2);
	    FREE(nt1);
	    );
      positions[pointers[reduced]++] = Bigendian_convert_uint(oldpositions[j]);
    }
#else
    for (j = oldoffsets[i]; j < oldoffsets[i+1]; j++) {
      debug(nt1 = shortoligo_nt(i,INDEX1PART);
	    nt2 = shortoligo_nt(reduced,INDEX1PART);
	    printf("Oligo %s => %s: copying position %u to location %u\n",
		   nt1,nt2,oldpositions[j],pointers[i]);
	    FREE(nt2);
	    FREE(nt1);
	    );
      positions[pointers[reduced]++] = oldpositions[j];
    }
#endif
  }
  fprintf(stderr,"\n");


  fprintf(stderr,"Sorting CT positions");
  /* Sort positions in each block */
  if (snps_root) {
    offsets_ptr = 0U;
    FWRITE_UINT(offsets_ptr,offsets_fp);
  }
  for (i = 0; i < oligospace; i++) {
    if (i % 100000 == 0) {
      fprintf(stderr,".");
    }
    block_start = offsets[i];
    block_end = offsets[i+1];
    if ((npositions = block_end - block_start) > 0) {
      qsort(&(positions[block_start]),npositions,sizeof(Genomicpos_T),Genomicpos_compare);
      if (snps_root == NULL) {
	FWRITE_UINTS(&(positions[block_start]),npositions,positions_fp);
      } else {
	FWRITE_UINT(positions[block_start],positions_fp);
	for (j = block_start+1; j < block_end; j++) {
	  if (positions[j] == positions[j-1]) {
	    npositions--;
	  } else {
	    FWRITE_UINT(positions[j],positions_fp);
	  }
	}
	offsets_ptr += npositions;
      }
    }
    if (snps_root) {
      FWRITE_UINT(offsets_ptr,offsets_fp);
    }
  }
  fprintf(stderr,"\n");

  if (snps_root == NULL) {
    FWRITE_UINTS(offsets,oligospace+1,offsets_fp);
  }

  FREE(positions);

  return;
}

static void
compute_positions_ga (FILE *offsets_fp, FILE *positions_fp, Positionsptr_T *offsets, Positionsptr_T *oldoffsets,
		      Genomicpos_T *oldpositions, int oligospace) {
  Genomicpos_T *positions;
  Storedoligomer_T reduced;
  int size, i, j;
  Positionsptr_T *pointers, preunique_totalcounts, block_start, block_end, npositions, offsets_ptr;
#ifdef DEBUG
  char *nt1, *nt2;
#endif

  preunique_totalcounts = offsets[oligospace];
  if (preunique_totalcounts == 0) {
    fprintf(stderr,"Something is wrong with the offsets.  Total counts is zero.\n");
    exit(9);
  } else {
    fprintf(stderr,"Trying to allocate %u*%d bytes of memory...",preunique_totalcounts,(int) sizeof(Genomicpos_T));
    positions = (Genomicpos_T *) CALLOC_NO_EXCEPTION(preunique_totalcounts,sizeof(Genomicpos_T));
    if (positions == NULL) {
      fprintf(stderr,"failed.  Need a computer with sufficient memory.\n");
      exit(9);
    } else {
      fprintf(stderr,"done\n");
    }
  }

  /* Copy offsets */
  pointers = (Positionsptr_T *) CALLOC(oligospace,sizeof(Positionsptr_T));
  for (i = 0; i < oligospace; i++) {
    pointers[i] = offsets[i];
  }

  fprintf(stderr,"Rearranging GA positions");
  for (i = 0; i < oligospace; i++) {
    if (i % 100000 == 0) {
      fprintf(stderr,".");
    }

    reduced = Cmet_reduce_ga(i);
#ifdef WORDS_BIGENDIAN
    for (j = Bigendian_convert_uint(oldoffsets[i]); j < Bigendian_convert_uint(oldoffsets[i+1]); j++) {
      debug(nt1 = shortoligo_nt(i,INDEX1PART);
	    nt2 = shortoligo_nt(reduced,INDEX1PART);
	    printf("Oligo %s => %s: copying position %u to location %u\n",
		   nt1,nt2,oldpositions[j],pointers[i]);
	    FREE(nt2);
	    FREE(nt1);
	    );
      positions[pointers[reduced]++] = Bigendian_convert_uint(oldpositions[j]);
    }
#else
    for (j = oldoffsets[i]; j < oldoffsets[i+1]; j++) {
      debug(nt1 = shortoligo_nt(i,INDEX1PART);
	    nt2 = shortoligo_nt(reduced,INDEX1PART);
	    printf("Oligo %s => %s: copying position %u to location %u\n",
		   nt1,nt2,oldpositions[j],pointers[i]);
	    FREE(nt2);
	    FREE(nt1);
	    );
      positions[pointers[reduced]++] = oldpositions[j];
    }
#endif
  }
  fprintf(stderr,"\n");


  fprintf(stderr,"Sorting GA positions");
  /* Sort positions in each block */
  if (snps_root) {
    offsets_ptr = 0U;
    FWRITE_UINT(offsets_ptr,offsets_fp);
  }
  for (i = 0; i < oligospace; i++) {
    if (i % 100000 == 0) {
      fprintf(stderr,".");
    }
    block_start = offsets[i];
    block_end = offsets[i+1];
    if ((npositions = block_end - block_start) > 0) {
      qsort(&(positions[block_start]),npositions,sizeof(Genomicpos_T),Genomicpos_compare);
      if (snps_root == NULL) {
	FWRITE_UINTS(&(positions[block_start]),npositions,positions_fp);
      } else {
	FWRITE_UINT(positions[block_start],positions_fp);
	for (j = block_start+1; j < block_end; j++) {
	  if (positions[j] == positions[j-1]) {
	    npositions--;
	  } else {
	    FWRITE_UINT(positions[j],positions_fp);
	  }
	}
	offsets_ptr += npositions;
      }
    }
    if (snps_root) {
      FWRITE_UINT(offsets_ptr,offsets_fp);
    }
  }
  fprintf(stderr,"\n");

  if (snps_root == NULL) {
    FWRITE_UINTS(offsets,oligospace+1,offsets_fp);
  }

  FREE(positions);
  return;
}


/* Usage: cmetindex -d <genome> */


/* Note: Depends on having gmapindex sampled on mod 3 bounds */
int
main (int argc, char *argv[]) {
  char *genomesubdir = NULL, *sourcedir = NULL, *destdir = NULL, *filename, *fileroot;
  FILE *input;
  Positionsptr_T *offsets_ct, *offsets_ga, *ref_offsets, totalcounts;
  Genomicpos_T *positions, *ref_positions, nblocks;
  int oligospace, i;

  FILE *offsets_fp, *positions_fp, *ref_offsets_fp, *ref_positions_fp;
  int ref_offsets_fd, ref_positions_fd;
  size_t ref_offsets_len, ref_positions_len;

  int opt, len, c;
  extern int optind;
  extern char *optarg;
  int long_option_index = 0;
  const char *long_name;

  while ((opt = getopt_long(argc,argv,"F:D:d:v:",
			    long_options,&long_option_index)) != -1) {
    switch (opt) {
    case 0: 
      long_name = long_options[long_option_index].name;
      if (!strcmp(long_name,"version")) {
	print_program_version();
	exit(0);
      } else if (!strcmp(long_name,"help")) {
	print_program_usage();
	exit(0);
      } else {
	/* Shouldn't reach here */
	fprintf(stderr,"Don't recognize option %s.  For usage, run 'cmetindex --help'",long_name);
	exit(9);
      }
      break;

    case 'F': user_sourcedir = optarg; break;
    case 'D': user_destdir = optarg; break;
    case 'd': dbroot = optarg; break;
    case 'v': snps_root = optarg; break;
    default: fprintf(stderr,"Do not recognize flag %c\n",opt); exit(9);
    }
  }
  argc -= (optind - 1);
  argv += (optind - 1);

  if (dbroot == NULL) {
    fprintf(stderr,"Missing name of genome database.  Must specify with -d flag.\n");
    fprintf(stderr,"Usage: cmetindex -d <genome>\n");
    exit(9);
  } else {
    sourcedir = Datadir_find_genomesubdir(&fileroot,&dbversion,user_sourcedir,dbroot);
    fprintf(stderr,"Reading source files from %s\n",sourcedir);
  }


  oligospace = power(4,INDEX1PART);

  /* Read offsets */
  if (snps_root == NULL) {
    filename = (char *) CALLOC(strlen(sourcedir)+strlen("/")+strlen(fileroot)+
			       strlen(".")+strlen(IDX_FILESUFFIX)+strlen("3")+
			       strlen(OFFSETS_FILESUFFIX)+1,sizeof(char));
    sprintf(filename,"%s/%s.%s%s%s",sourcedir,fileroot,IDX_FILESUFFIX,"3",OFFSETS_FILESUFFIX);

  } else {
    filename = (char *) CALLOC(strlen(sourcedir)+strlen("/")+strlen(fileroot)+
			       strlen(".")+strlen(IDX_FILESUFFIX)+strlen("3")+
			       strlen(OFFSETS_FILESUFFIX)+strlen(".")+strlen(snps_root)+1,sizeof(char));
    sprintf(filename,"%s/%s.%s%s%s.%s",sourcedir,fileroot,IDX_FILESUFFIX,"3",OFFSETS_FILESUFFIX,snps_root);
  }
  if ((ref_offsets_fp = FOPEN_READ_BINARY(filename)) == NULL) {
    fprintf(stderr,"Can't open file %s\n",filename);
    exit(9);
  }
  ref_offsets = (Positionsptr_T *) Access_mmap(&ref_offsets_fd,&ref_offsets_len,
					       filename,sizeof(Positionsptr_T),/*randomp*/false);
  FREE(filename);

  /* Read positions */
  if (snps_root == NULL) {
    filename = (char *) CALLOC(strlen(sourcedir)+strlen("/")+strlen(fileroot)+
			       strlen(".")+strlen(IDX_FILESUFFIX)+strlen("3")+
			       strlen(POSITIONS_FILESUFFIX)+1,sizeof(char));
    sprintf(filename,"%s/%s.%s%s%s",sourcedir,fileroot,IDX_FILESUFFIX,"3",POSITIONS_FILESUFFIX);
  } else {
    filename = (char *) CALLOC(strlen(sourcedir)+strlen("/")+strlen(fileroot)+
			       strlen(".")+strlen(IDX_FILESUFFIX)+strlen("3")+
			       strlen(POSITIONS_FILESUFFIX)+strlen(".")+strlen(snps_root)+1,sizeof(char));
    sprintf(filename,"%s/%s.%s%s%s.%s",sourcedir,fileroot,IDX_FILESUFFIX,"3",POSITIONS_FILESUFFIX,snps_root);
  }
  if ((ref_positions_fp = FOPEN_READ_BINARY(filename)) == NULL) {
    fprintf(stderr,"Can't open file %s\n",filename);
    exit(9);
  }
  ref_positions = (Genomicpos_T *) Access_mmap(&ref_positions_fd,&ref_positions_len,
					       filename,sizeof(Genomicpos_T),/*randomp*/false);
  FREE(filename);


  /* Open CT output files */
  if (user_destdir == NULL) {
    destdir = sourcedir;
  } else {
    destdir = user_destdir;
  }
  fprintf(stderr,"Reading cmet index files to %s\n",destdir);

  if (snps_root == NULL) {
    filename = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
			       strlen(".")+strlen("metct")+strlen("3")+
			       strlen(OFFSETS_FILESUFFIX)+1,sizeof(char));
    sprintf(filename,"%s/%s.%s%s%s",destdir,fileroot,"metct","3",OFFSETS_FILESUFFIX);
  } else {
    filename = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
			       strlen(".")+strlen("metct")+strlen("3")+
			       strlen(OFFSETS_FILESUFFIX)+strlen(".")+strlen(snps_root)+1,sizeof(char));
    sprintf(filename,"%s/%s.%s%s%s.%s",destdir,fileroot,"metct","3",OFFSETS_FILESUFFIX,snps_root);
  }
  if ((offsets_fp = FOPEN_WRITE_BINARY(filename)) == NULL) {
    fprintf(stderr,"Can't open file %s for writing\n",filename);
    exit(9);
  }
  FREE(filename);

  if (snps_root == NULL) {
    filename = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
			       strlen(".")+strlen("metct")+strlen("3")+
			       strlen(POSITIONS_FILESUFFIX)+1,sizeof(char));
    sprintf(filename,"%s/%s.%s%s%s",destdir,fileroot,"metct","3",POSITIONS_FILESUFFIX);
  } else {
    filename = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
			       strlen(".")+strlen("metct")+strlen("3")+
			       strlen(POSITIONS_FILESUFFIX)+strlen(".")+strlen(snps_root)+1,sizeof(char));
    sprintf(filename,"%s/%s.%s%s%s.%s",destdir,fileroot,"metct","3",POSITIONS_FILESUFFIX,snps_root);
  }
  if ((positions_fp = FOPEN_WRITE_BINARY(filename)) == NULL) {
    fprintf(stderr,"Can't open file %s for writing\n",filename);
    exit(9);
  }
  FREE(filename);

  /* Compute and write CT files */
  offsets_ct = compute_offsets_ct(ref_offsets,oligospace);
  compute_positions_ct(offsets_fp,positions_fp,offsets_ct,ref_offsets,ref_positions,oligospace);
  FREE(offsets_ct);
  fclose(positions_fp);
  fclose(offsets_fp);



  /* Open GA output files */
  if (snps_root == NULL) {
    filename = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
			       strlen(".")+strlen("metga")+strlen("3")+
			       strlen(OFFSETS_FILESUFFIX)+1,sizeof(char));
    sprintf(filename,"%s/%s.%s%s%s",destdir,fileroot,"metga","3",OFFSETS_FILESUFFIX);
  } else {
    filename = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
			       strlen(".")+strlen("metga")+strlen("3")+
			       strlen(OFFSETS_FILESUFFIX)+strlen(".")+strlen(snps_root)+1,sizeof(char));
    sprintf(filename,"%s/%s.%s%s%s.%s",destdir,fileroot,"metga","3",OFFSETS_FILESUFFIX,snps_root);
  }
  if ((offsets_fp = FOPEN_WRITE_BINARY(filename)) == NULL) {
    fprintf(stderr,"Can't open file %s for writing\n",filename);
    exit(9);
  }
  FREE(filename);

  if (snps_root == NULL) {
    filename = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
			       strlen(".")+strlen("metga")+strlen("3")+
			       strlen(POSITIONS_FILESUFFIX)+1,sizeof(char));
    sprintf(filename,"%s/%s.%s%s%s",destdir,fileroot,"metga","3",POSITIONS_FILESUFFIX);
  } else {
    filename = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
			       strlen(".")+strlen("metga")+strlen("3")+
			       strlen(POSITIONS_FILESUFFIX)+strlen(".")+strlen(snps_root)+1,sizeof(char));
    sprintf(filename,"%s/%s.%s%s%s.%s",destdir,fileroot,"metga","3",POSITIONS_FILESUFFIX,snps_root);
  }
  if ((positions_fp = FOPEN_WRITE_BINARY(filename)) == NULL) {
    fprintf(stderr,"Can't open file %s for writing\n",filename);
    exit(9);
  }
  FREE(filename);

  /* Compute and write GA files */
  offsets_ga = compute_offsets_ga(ref_offsets,oligospace);
  compute_positions_ga(offsets_fp,positions_fp,offsets_ga,ref_offsets,ref_positions,oligospace);
  FREE(offsets_ga);
  fclose(positions_fp);
  fclose(offsets_fp);


  /* Clean up */
  munmap((void *) ref_positions,ref_positions_len);
  close(ref_positions_fd);

  munmap((void *) ref_offsets,ref_offsets_len);
  close(ref_offsets_fd);

  return 0;
}



static void
print_program_usage () {
  fprintf(stdout,"\
Usage: cmetindex [OPTIONS...] -d <genome>\n\
\n\
");

  /* Input options */
  fprintf(stdout,"Options (must include -d)\n");
  fprintf(stdout,"\
  -F, --sourcedir=directory      Directory where to read cmet index files (default is\n\
                                   GMAP genome directory specified at compile time)\n\
  -D, --destdir=directory        Directory where to write cmet index files (default is\n\
                                   GMAP genome directory specified at compile time)\n\
  -d, --db=STRING                Genome database\n\
  -v, --use-snps=STRING          Use database containing known SNPs (in <STRING>.iit, built\n\
                                   previously using snpindex) for tolerance to SNPs\n\
\n\
  --version                      Show version\n\
  --help                         Show this help message\n\
");
  return;
}

