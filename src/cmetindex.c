static char rcsid[] = "$Id: cmetindex.c 60849 2012-03-30 19:25:41Z twu $";
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
#include "types.h"		/* For Positionsptr_T, Oligospace_T, and Storedoligomer_T */

#include "cmet.h"

#include "bool.h"
#include "genomicpos.h"
#include "indexdb.h"
#include "datadir.h"
#include "getopt.h"


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

#define MONITOR_INTERVAL 10000000


static char *user_sourcedir = NULL;
static char *user_destdir = NULL;
static char *dbroot = NULL;
static char *dbversion = NULL;
static int offsetscomp_basesize = 12;
static int required_basesize = 0;
static int index1part = 15;
static int required_index1part = 0;
static int index1interval;
static int required_interval = 0;

static char *snps_root = NULL;


static struct option long_options[] = {
  /* Input options */
  {"sourcedir", required_argument, 0, 'F'},	/* user_sourcedir */
  {"destdir", required_argument, 0, 'D'},	/* user_destdir */
  {"basesize", required_argument, 0, 'b'},	/* required_basesize */
  {"kmer", required_argument, 0, 'k'}, /* required_index1part */
  {"sampling", required_argument, 0, 'q'}, /* required_interval */
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


static Oligospace_T
power (int base, int exponent) {
  Oligospace_T result = 1UL;
  int i;

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
compute_offsets_ct (Positionsptr_T *oldoffsets, Oligospace_T oligospace, Storedoligomer_T mask) {
  Positionsptr_T *offsets;
  Oligospace_T oligoi, reduced;
#ifdef DEBUG
  char *nt1, *nt2;
#endif

  /* Fill with sizes */
  offsets = (Positionsptr_T *) CALLOC(oligospace+1,sizeof(Positionsptr_T));

  for (oligoi = 0; oligoi < oligospace; oligoi++) {
#if 0
    if (reduce_oligo(oligoi) != reduce_oligo_old(oligoi,index1part)) {
      abort();
    }
#endif
    reduced = Cmet_reduce_ct(oligoi) & mask;
    debug(
	  nt1 = shortoligo_nt(oligoi,index1part);
	  nt2 = shortoligo_nt(reduced,index1part);
	  printf("For oligo %s, updating sizes for %s from %u",nt1,nt2,offsets[reduced+1]);
	  );
#ifdef WORDS_BIGENDIAN
    /*size*/offsets[reduced+1] += (Bigendian_convert_uint(oldoffsets[oligoi+1]) - Bigendian_convert_uint(oldoffsets[oligoi]));
#else
    /*size*/offsets[reduced+1] += (oldoffsets[oligoi+1] - oldoffsets[oligoi]);
#endif
    debug(
	  printf(" to %u\n",offsets[reduced+1]);
	  FREE(nt2);
	  FREE(nt1);
	  );
  }

  offsets[0] = 0U;
  for (oligoi = 1; oligoi <= oligospace; oligoi++) {
    offsets[oligoi] = offsets[oligoi-1] + /*size*/offsets[oligoi];
    debug(if (offsets[oligoi] != offsets[oligoi-1]) {
	    printf("Offset for %X: %u\n",oligoi,offsets[oligoi]);
	  });
  }

  return offsets;
}


static Positionsptr_T *
compute_offsets_ga (Positionsptr_T *oldoffsets, Oligospace_T oligospace, Storedoligomer_T mask) {
  Positionsptr_T *offsets;
  Oligospace_T oligoi, reduced;
#ifdef DEBUG
  char *nt1, *nt2;
#endif

  /* Fill with sizes */
  offsets = (Positionsptr_T *) CALLOC(oligospace+1,sizeof(Positionsptr_T));

  for (oligoi = 0; oligoi < oligospace; oligoi++) {
    reduced = Cmet_reduce_ga(oligoi) & mask;
    debug(
	  nt1 = shortoligo_nt(oligoi,index1part);
	  nt2 = shortoligo_nt(reduced,index1part);
	  printf("For oligo %s, updating sizes for %s from %u",nt1,nt2,offsets[reduced+1]);
	  );
#ifdef WORDS_BIGENDIAN
    /*size*/offsets[reduced+1] += (Bigendian_convert_uint(oldoffsets[oligoi+1]) - Bigendian_convert_uint(oldoffsets[oligoi]));
#else
    /*size*/offsets[reduced+1] += (oldoffsets[oligoi+1] - oldoffsets[oligoi]);
#endif
    debug(
	  printf(" to %d\n",offsets[reduced+1]);
	  FREE(nt2);
	  FREE(nt1);
	  );
  }

  offsets[0] = 0U;
  for (oligoi = 1; oligoi <= oligospace; oligoi++) {
    offsets[oligoi] = offsets[oligoi-1] + /*size*/offsets[oligoi];
    debug(if (offsets[oligoi] != offsets[oligoi-1]) {
	    printf("Offset for %X: %u\n",oligoi,offsets[oligoi]);
	  });
  }

  return offsets;
}


static void
compute_ct (char *gammaptrs_filename, char *offsetscomp_filename,
	    FILE *positions_fp, Positionsptr_T *oldoffsets, Genomicpos_T *oldpositions,
	    Oligospace_T oligospace, Storedoligomer_T mask) {
  Genomicpos_T *positions;
  Positionsptr_T *snpoffsets, j;
  Oligospace_T oligoi, oligok, reduced;
  Positionsptr_T *pointers, *offsets, preunique_totalcounts, block_start, block_end, npositions, offsets_ptr;
#ifdef DEBUG
  char *nt1, *nt2;
#endif

  offsets = compute_offsets_ct(oldoffsets,oligospace,mask);

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

  /* Point to offsets and revise (previously copied) */
  pointers = offsets;
  fprintf(stderr,"Rearranging CT positions...");

  for (oligoi = 0; oligoi < oligospace; oligoi++) {
    if (oligoi % MONITOR_INTERVAL == 0) {
      fprintf(stderr,".");
    }
    reduced = Cmet_reduce_ct(oligoi) & mask;
#ifdef WORDS_BIGENDIAN
    for (j = Bigendian_convert_uint(oldoffsets[oligoi]); j < Bigendian_convert_uint(oldoffsets[oligoi+1]); j++) {
      debug(nt1 = shortoligo_nt(oligoi,index1part);
	    nt2 = shortoligo_nt(reduced,index1part);
	    printf("Oligo %s => %s: copying position %u to location %u\n",
		   nt1,nt2,oldpositions[j],pointers[oligoi]);
	    FREE(nt2);
	    FREE(nt1);
	    );
      positions[pointers[reduced]++] = Bigendian_convert_uint(oldpositions[j]);
    }
#else
    for (j = oldoffsets[oligoi]; j < oldoffsets[oligoi+1]; j++) {
      debug(nt1 = shortoligo_nt(oligoi,index1part);
	    nt2 = shortoligo_nt(reduced,index1part);
	    printf("Oligo %s => %s: copying position %u to location %u\n",
		   nt1,nt2,oldpositions[j],pointers[oligoi]);
	    FREE(nt2);
	    FREE(nt1);
	    );
      positions[pointers[reduced]++] = oldpositions[j];
    }
#endif
  }
  FREE(offsets);
  fprintf(stderr,"done\n");


  fprintf(stderr,"Sorting CT positions...");
  offsets = compute_offsets_ct(oldoffsets,oligospace,mask);

  /* Sort positions in each block */
  if (snps_root) {
    oligok = 0;
    snpoffsets = (Positionsptr_T *) CALLOC(oligospace+1,sizeof(Positionsptr_T));
    offsets_ptr = 0U;
    snpoffsets[oligok++] = offsets_ptr;
  }
  for (oligoi = 0; oligoi < oligospace; oligoi++) {
    if (oligoi % MONITOR_INTERVAL == 0) {
      fprintf(stderr,".");
    }
    block_start = offsets[oligoi];
    block_end = offsets[oligoi+1];
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
      snpoffsets[oligok++] = offsets_ptr;
    }
  }
  fprintf(stderr,"done\n");

  if (snps_root == NULL) {
#ifdef PRE_GAMMA
    FWRITE_UINTS(offsets,oligospace+1,offsets_fp);
#else
    Indexdb_write_gammaptrs(gammaptrs_filename,offsetscomp_filename,offsets,oligospace,
			    /*blocksize*/power(4,index1part - offsetscomp_basesize));
#endif
  } else {
    Indexdb_write_gammaptrs(gammaptrs_filename,offsetscomp_filename,snpoffsets,oligospace,
			    /*blocksize*/power(4,index1part - offsetscomp_basesize));
    FREE(snpoffsets);
  }

  FREE(offsets);
  FREE(positions);

  return;
}

static void
compute_ga (char *gammaptrs_filename, char *offsetscomp_filename,
	    FILE *positions_fp, Positionsptr_T *oldoffsets, Genomicpos_T *oldpositions,
	    Oligospace_T oligospace, Storedoligomer_T mask) {
  Genomicpos_T *positions;
  Positionsptr_T *snpoffsets, j;
  Oligospace_T oligoi, oligok, reduced;
  Positionsptr_T *pointers, *offsets, preunique_totalcounts, block_start, block_end, npositions, offsets_ptr;
#ifdef DEBUG
  char *nt1, *nt2;
#endif

  offsets = compute_offsets_ga(oldoffsets,oligospace,mask);

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

  /* Point to offsets and revise (previously copied) */
  pointers = offsets;
  fprintf(stderr,"Rearranging GA positions...");

  for (oligoi = 0; oligoi < oligospace; oligoi++) {
    if (oligoi % MONITOR_INTERVAL == 0) {
      fprintf(stderr,".");
    }

    reduced = Cmet_reduce_ga(oligoi) & mask;
#ifdef WORDS_BIGENDIAN
    for (j = Bigendian_convert_uint(oldoffsets[oligoi]); j < Bigendian_convert_uint(oldoffsets[oligoi+1]); j++) {
      debug(nt1 = shortoligo_nt(oligoi,index1part);
	    nt2 = shortoligo_nt(reduced,index1part);
	    printf("Oligo %s => %s: copying position %u to location %u\n",
		   nt1,nt2,oldpositions[j],pointers[oligoi]);
	    FREE(nt2);
	    FREE(nt1);
	    );
      positions[pointers[reduced]++] = Bigendian_convert_uint(oldpositions[j]);
    }
#else
    for (j = oldoffsets[oligoi]; j < oldoffsets[oligoi+1]; j++) {
      debug(nt1 = shortoligo_nt(oligoi,index1part);
	    nt2 = shortoligo_nt(reduced,index1part);
	    printf("Oligo %s => %s: copying position %u to location %u\n",
		   nt1,nt2,oldpositions[j],pointers[oligoi]);
	    FREE(nt2);
	    FREE(nt1);
	    );
      positions[pointers[reduced]++] = oldpositions[j];
    }
#endif
  }
  FREE(offsets);
  fprintf(stderr,"done\n");


  fprintf(stderr,"Sorting GA positions...");
  offsets = compute_offsets_ga(oldoffsets,oligospace,mask);

  /* Sort positions in each block */
  if (snps_root) {
    oligok = 0;
    snpoffsets = (Positionsptr_T *) CALLOC(oligospace+1,sizeof(Positionsptr_T));
    offsets_ptr = 0U;
    snpoffsets[oligok++] = offsets_ptr;
  }
  for (oligoi = 0; oligoi < oligospace; oligoi++) {
    if (oligoi % MONITOR_INTERVAL == 0) {
      fprintf(stderr,".");
    }
    block_start = offsets[oligoi];
    block_end = offsets[oligoi+1];
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
      snpoffsets[oligok++] = offsets_ptr;
    }
  }
  fprintf(stderr,"done\n");

  if (snps_root == NULL) {
#ifdef PRE_GAMMAS
    FWRITE_UINTS(offsets,oligospace+1,offsets_fp);
#else
    Indexdb_write_gammaptrs(gammaptrs_filename,offsetscomp_filename,offsets,oligospace,
			    /*blocksize*/power(4,index1part - offsetscomp_basesize));
#endif
  } else {
    Indexdb_write_gammaptrs(gammaptrs_filename,offsetscomp_filename,snpoffsets,oligospace,
			    /*blocksize*/power(4,index1part - offsetscomp_basesize));
    FREE(snpoffsets);
  }

  FREE(offsets);
  FREE(positions);
  return;
}


/* Usage: cmetindex -d <genome> */


/* Note: Depends on having gmapindex sampled on mod 3 bounds */
int
main (int argc, char *argv[]) {
  char *sourcedir = NULL, *destdir = NULL, *filename, *fileroot;
  char *gammaptrs_filename, *offsetscomp_filename, *positions_filename,
    *gammaptrs_basename_ptr, *offsetscomp_basename_ptr, *positions_basename_ptr,
    *gammaptrs_index1info_ptr, *offsetscomp_index1info_ptr, *positions_index1info_ptr;
  char *new_gammaptrs_filename, *new_offsetscomp_filename;
  Positionsptr_T *ref_offsets;
  Storedoligomer_T mask;
  Genomicpos_T *ref_positions;
  Oligospace_T oligospace;

  FILE *positions_fp, *ref_positions_fp;
  int ref_positions_fd;
  size_t ref_positions_len;
#ifndef HAVE_MMAP
  double seconds;
#endif

  int opt;
  extern int optind;
  extern char *optarg;
  int long_option_index = 0;
  const char *long_name;

  while ((opt = getopt_long(argc,argv,"F:D:d:b:k:q:v:",
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
    case 'b': required_basesize = atoi(optarg); break;
    case 'k': required_index1part = atoi(optarg); break;
    case 'q': required_interval = atoi(optarg); break;
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


  Indexdb_get_filenames(&gammaptrs_filename,&offsetscomp_filename,&positions_filename,
			&gammaptrs_basename_ptr,&offsetscomp_basename_ptr,&positions_basename_ptr,
			&gammaptrs_index1info_ptr,&offsetscomp_index1info_ptr,&positions_index1info_ptr,
			&offsetscomp_basesize,&index1part,&index1interval,
			sourcedir,fileroot,IDX_FILESUFFIX,snps_root,
			required_basesize,required_index1part,required_interval);

  mask = ~(~0UL << 2*index1part);
  oligospace = power(4,index1part);

  /* Read offsets */
  ref_offsets = Indexdb_offsets_from_gammas(gammaptrs_filename,offsetscomp_filename,offsetscomp_basesize,index1part);

  /* Read positions */
  if ((ref_positions_fp = FOPEN_READ_BINARY(positions_filename)) == NULL) {
    fprintf(stderr,"Can't open file %s\n",positions_filename);
    exit(9);
  }

#ifdef HAVE_MMAP
  ref_positions = (Genomicpos_T *) Access_mmap(&ref_positions_fd,&ref_positions_len,
					       positions_filename,sizeof(Genomicpos_T),/*randomp*/false);
#else
  ref_positions = (Genomicpos_T *) Access_allocated(&ref_positions_len,&seconds,
						    positions_filename,sizeof(Genomicpos_T));
#endif


  /* Open CT output files */
  if (user_destdir == NULL) {
    destdir = sourcedir;
  } else {
    destdir = user_destdir;
  }
  fprintf(stderr,"Writing cmet index files to %s\n",destdir);

  if (index1part == offsetscomp_basesize) {
    new_gammaptrs_filename = (char *) NULL;
  } else {
    new_gammaptrs_filename = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
					     strlen(".")+strlen("metct")+strlen(gammaptrs_index1info_ptr)+1,sizeof(char));
    sprintf(new_gammaptrs_filename,"%s/%s.%s%s",destdir,fileroot,"metct",gammaptrs_index1info_ptr);
  }

  new_offsetscomp_filename = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
			     strlen(".")+strlen("metct")+strlen(offsetscomp_index1info_ptr)+1,sizeof(char));
  sprintf(new_offsetscomp_filename,"%s/%s.%s%s",destdir,fileroot,"metct",offsetscomp_index1info_ptr);


  filename = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
			     strlen(".")+strlen("metct")+strlen(positions_index1info_ptr)+1,sizeof(char));
  sprintf(filename,"%s/%s.%s%s",destdir,fileroot,"metct",positions_index1info_ptr);

  if ((positions_fp = FOPEN_WRITE_BINARY(filename)) == NULL) {
    fprintf(stderr,"Can't open file %s for writing\n",filename);
    exit(9);
  }
  FREE(filename);

  /* Compute and write CT files */
  compute_ct(new_gammaptrs_filename,new_offsetscomp_filename,
	     positions_fp,ref_offsets,ref_positions,oligospace,mask);
  fclose(positions_fp);
  FREE(new_offsetscomp_filename);
  if (index1part != offsetscomp_basesize) {
    FREE(new_gammaptrs_filename);
  }


  /* Open GA output files */
  if (index1part == offsetscomp_basesize) {
    new_gammaptrs_filename = (char *) NULL;
  } else {
    new_gammaptrs_filename = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
					     strlen(".")+strlen("metga")+strlen(gammaptrs_index1info_ptr)+1,sizeof(char));
    sprintf(new_gammaptrs_filename,"%s/%s.%s%s",destdir,fileroot,"metga",gammaptrs_index1info_ptr);
  }

  new_offsetscomp_filename = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
			     strlen(".")+strlen("metga")+strlen(offsetscomp_index1info_ptr)+1,sizeof(char));
  sprintf(new_offsetscomp_filename,"%s/%s.%s%s",destdir,fileroot,"metga",offsetscomp_index1info_ptr);


  filename = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
			     strlen(".")+strlen("metga")+strlen(positions_index1info_ptr)+1,sizeof(char));
  sprintf(filename,"%s/%s.%s%s",destdir,fileroot,"metga",positions_index1info_ptr);

  if ((positions_fp = FOPEN_WRITE_BINARY(filename)) == NULL) {
    fprintf(stderr,"Can't open file %s for writing\n",filename);
    exit(9);
  }
  FREE(filename);

  /* Compute and write GA files */
  compute_ga(new_gammaptrs_filename,new_offsetscomp_filename,
	     positions_fp,ref_offsets,ref_positions,oligospace,mask);
  fclose(positions_fp);
  FREE(new_offsetscomp_filename);
  if (index1part != offsetscomp_basesize) {
    FREE(new_gammaptrs_filename);
  }



  /* Clean up */
  FREE(ref_offsets);

#ifdef HAVE_MMAP
  munmap((void *) ref_positions,ref_positions_len);
  close(ref_positions_fd);
#else
  FREE(ref_positions);
#endif

  FREE(positions_filename);
  FREE(offsetscomp_filename);
  FREE(gammaptrs_filename);

  FREE(dbversion);
  FREE(fileroot);
  FREE(sourcedir);

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
                                   value of -F, if provided; otherwise the value of the\n\
                                   GMAP genome directory specified at compile time)\n\
  -d, --db=STRING                Genome database\n\
  -k, --kmer=INT                 kmer size to use in genome database (allowed values: 16 or less).\n\
                                   If not specified, the program will find the highest available\n\
                                   kmer size in the genome database\n\
  -b, --basesize=INT             Base size to use in genome database.  If not specified, the program\n\
                                   will find the highest available base size in the genome database\n\
                                   within selected k-mer size\n\
  -q, --sampling=INT             Sampling to use in genome database.  If not specified, the program\n\
                                   will find the smallest available sampling value in the genome database\n\
                                   within selected basesize and k-mer size\n\
  -v, --use-snps=STRING          Use database containing known SNPs (in <STRING>.iit, built\n\
                                   previously using snpindex) for tolerance to SNPs\n\
\n\
  --version                      Show version\n\
  --help                         Show this help message\n\
");
  return;
}

