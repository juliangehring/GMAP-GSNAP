static char rcsid[] = "$Id: indexdb.c 131819 2014-03-28 23:21:14Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifndef HAVE_MEMCPY
# define memcpy(d,s,n) bcopy((s),(d),(n))
#endif
#ifndef HAVE_MEMMOVE
# define memmove(d,s,n) bcopy((s),(d),(n))
#endif
#ifdef HAVE_SSE2
#include <emmintrin.h>
#endif

#include "indexdb.h"
#include "indexdbdef.h"
#include "genome_hr.h"		/* For read_gammas procedures */
#include "bitpack64-read.h"


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

#include "compress.h"
#include "interval.h"
#include "complement.h"
#include "bitpack64-read.h"	/* For Bitpack64_read_setup() */


#ifdef HAVE_PTHREAD
#include <pthread.h>		/* sys/types.h already included above */
#endif

#define MAXENTRIES 20

/* Note: NONMODULAR is the old behavior.  Now we store only when
   startposition % index1interval == 0 */


/* Low-level codon hacking */
#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Calls to Indexdb_read */
#ifdef DEBUG0
#define debug0(x) x
#else
#define debug0(x)
#endif

/* Shifting of high word to low word for PMAP */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif


#ifdef PMAP

#if (defined(DEBUG) || defined(DEBUG0))
static Width_T index1part_aa;
#endif

void
Indexdb_setup (Width_T index1part_aa_in) {
#if (defined(DEBUG) || defined(DEBUG0))
  index1part_aa = index1part_aa_in;
#endif
  return;
}

#else

#define poly_A 0U
static Storedoligomer_T poly_T;  /* Was LOW12MER 0x00FFFFFF */

#if (defined(DEBUG) || defined(DEBUG0))
static Width_T index1part;
#endif

void
Indexdb_setup (Width_T index1part_in) {
#if (defined(DEBUG) || defined(DEBUG0))
  index1part = index1part_in;
#endif

  poly_T = ~(~0UL << 2*index1part_in);
  return;
}
#endif


#define T Indexdb_T


void
Indexdb_free (T *old) {
  if (*old) {
    if ((*old)->positions_access == ALLOCATED) {
#ifdef LARGE_GENOMES
      FREE((*old)->positions_high);
      FREE((*old)->positions_low);
#else
      FREE((*old)->positions);
#endif
#ifdef HAVE_MMAP
    } else if ((*old)->positions_access == MMAPPED) {
#ifdef LARGE_GENOMES
      munmap((void *) (*old)->positions_high,(*old)->positions_high_len);
      close((*old)->positions_high_fd);
      munmap((void *) (*old)->positions_high,(*old)->positions_low_len);
      close((*old)->positions_low_fd);
#else
      munmap((void *) (*old)->positions,(*old)->positions_len);
      close((*old)->positions_fd);
#endif
#endif
    } else if ((*old)->positions_access == FILEIO) {
#ifdef HAVE_PTHREAD
      pthread_mutex_destroy(&(*old)->positions_read_mutex);
#endif
#ifdef LARGE_GENOMES
      close((*old)->positions_high_fd);
      close((*old)->positions_low_fd);
#else
      close((*old)->positions_fd);
#endif
    }

    if ((*old)->offsetscomp_access == ALLOCATED) {
      FREE((*old)->offsetscomp);
#ifdef HAVE_MMAP
    } else if ((*old)->offsetscomp_access == MMAPPED) {
      munmap((void *) (*old)->offsetscomp,(*old)->offsetscomp_len);
      close((*old)->offsetscomp_fd);
#endif
    }
      
    FREE((*old)->gammaptrs);	/* Always ALLOCATED */

    FREE(*old);
  }
  return;
}


Width_T
Indexdb_interval (T this) {
  return this->index1interval;
}


bool
Indexdb_positions_fileio_p (T this) {
  if (this->positions_access == FILEIO) {
    return true;
  } else {
    return false;
  }
}

static Oligospace_T
power (int base, Width_T exponent) {
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

double
Indexdb_mean_size (T this, Mode_T mode, Width_T index1part) {
  Oligospace_T oligospace, n;

#ifdef PMAP
  /* index1part should be in aa */
  n = oligospace = power(this->alphabet_size,index1part);
#else
  n = oligospace = power(4,index1part);
  if (mode == CMET_STRANDED || mode == CMET_NONSTRANDED || mode == ATOI_STRANDED || mode == ATOI_NONSTRANDED) {
    n = power(3,index1part);
  }
#endif

#ifdef WORDS_BIGENDIAN
  if (this->offsetscomp_access == ALLOCATED) {
    return (double) this->offsetscomp[this->gammaptrs[oligospace/this->offsetscomp_blocksize]]/(double) n;
  } else {
    return (double) Bigendian_convert_uint(this->offsetscomp[this->gammaptrs[oligospace/this->offsetscomp_blocksize]])/(double) n;
  }
#else
  return (double) this->offsetscomp[this->gammaptrs[oligospace/this->offsetscomp_blocksize]]/(double) n;
#endif
}



static Filenames_T
Filenames_new (char *pages_filename, char *pointers_filename, char *offsets_filename,
	       char *positions_high_filename, char *positions_low_filename,
	       char *pointers_basename_ptr, char *offsets_basename_ptr,
	       char *positions_high_basename_ptr, char *positions_low_basename_ptr,
	       char *pointers_index1info_ptr, char *offsets_index1info_ptr,
	       char *positions_high_index1info_ptr, char *positions_low_index1info_ptr) {
  Filenames_T new = (Filenames_T) MALLOC(sizeof(*new));

  new->pages_filename = pages_filename;
  new->pointers_filename = pointers_filename;
  new->offsets_filename = offsets_filename;
  new->positions_high_filename = positions_high_filename;
  new->positions_low_filename = positions_low_filename;

  new->pointers_basename_ptr = pointers_basename_ptr;
  new->offsets_basename_ptr = offsets_basename_ptr;
  new->positions_high_basename_ptr = positions_high_basename_ptr;
  new->positions_low_basename_ptr = positions_low_basename_ptr;

  new->pointers_index1info_ptr = pointers_index1info_ptr;
  new->offsets_index1info_ptr = offsets_index1info_ptr;
  new->positions_high_index1info_ptr = positions_high_index1info_ptr;
  new->positions_low_index1info_ptr = positions_low_index1info_ptr;

  return new;
}

void
Filenames_free (Filenames_T *old) {

  FREE((*old)->pages_filename);
  FREE((*old)->pointers_filename);
  FREE((*old)->offsets_filename);
  FREE((*old)->positions_high_filename);
  FREE((*old)->positions_low_filename);

  FREE(*old);

  return;
}


Filenames_T
Indexdb_get_filenames_no_compression (Width_T *index1part, Width_T *index1interval,
				      char *genomesubdir, char *fileroot, char *idx_filesuffix, char *snps_root,
				      Width_T required_interval, bool offsets_only_p) {
  char *offsets_filename, *positions_high_filename, *positions_low_filename,
    *offsets_basename_ptr, *positions_high_basename_ptr, *positions_low_basename_ptr,
    *offsets_index1info_ptr, *positions_high_index1info_ptr, *positions_low_index1info_ptr;

  char *base_filename, *filename;
  char *pattern, interval_char, digit_string[2], *p, *q;
  char tens, ones, ones0;
  Width_T found_index1part, found_interval;
  int rootlength, patternlength;

  char *offsets_suffix, *positions_high_suffix, *positions_low_suffix;
  struct dirent *entry;
  DIR *dp;


  if (snps_root == NULL) {
    offsets_suffix = "offsets";
    positions_high_suffix = POSITIONS_HIGH_FILESUFFIX;
    positions_low_suffix = POSITIONS_LOW_FILESUFFIX;
  } else {
    offsets_suffix = (char *) CALLOC(strlen("offsets")+strlen(".")+strlen(snps_root)+1,sizeof(char));
    sprintf(offsets_suffix,"%s.%s","offsets",snps_root);
    positions_high_suffix = (char *) CALLOC(strlen(POSITIONS_HIGH_FILESUFFIX)+strlen(".")+strlen(snps_root)+1,sizeof(char));
    sprintf(positions_high_suffix,"%s.%s",POSITIONS_HIGH_FILESUFFIX,snps_root);
    positions_low_suffix = (char *) CALLOC(strlen(POSITIONS_LOW_FILESUFFIX)+strlen(".")+strlen(snps_root)+1,sizeof(char));
    sprintf(positions_low_suffix,"%s.%s",POSITIONS_LOW_FILESUFFIX,snps_root);
  }

  *index1part = 0;
  *index1interval = 1000;
  base_filename = (char *) NULL;

  if ((dp = opendir(genomesubdir)) == NULL) {
    fprintf(stderr,"Unable to open directory %s\n",genomesubdir);
    exit(9);
  }

  pattern = (char *) CALLOC(strlen(fileroot)+strlen(".")+strlen(idx_filesuffix)+1,sizeof(char));
  sprintf(pattern,"%s.%s",fileroot,idx_filesuffix);
  patternlength = strlen(pattern);

  digit_string[1] = '\0';	/* Needed for atoi */
  while ((entry = readdir(dp)) != NULL) {
    filename = entry->d_name;
    if (!strncmp(filename,pattern,patternlength)) {
      p = &(filename[strlen(pattern)]); /* Points after idx_filesuffix, e.g., "ref" */
      if ((q = strstr(p,offsets_suffix)) != NULL && !strcmp(q,offsets_suffix)) {
	if (q - p == 3) {
#ifdef PMAP
	  /* e.g., pf677 */
	  if (sscanf(p,"%c%c%c",&ones0,&ones,&interval_char) == 3) {
	    /* digit_string[0] = ones0; */
	    /* found_basesize = atoi(digit_string); */

	    digit_string[0] = ones;
	    found_index1part = atoi(digit_string);
	    
	    digit_string[0] = interval_char;
	    found_interval = atoi(digit_string);
	  } else {
	    abort();
	  }
#else
	  /* e.g., ref123offsets */
	  if (sscanf(p,"%c%c%c",&tens,&ones,&interval_char) == 3) {
	    digit_string[0] = tens;
	    found_index1part = 10*atoi(digit_string);
	    digit_string[0] = ones;
	    found_index1part += atoi(digit_string);

	    digit_string[0] = interval_char;
	    found_interval = atoi(digit_string);
	  } else {
	    abort();
	  }
#endif

	} else if (q - p == 1) {
	  /* Old style, e.g, idx or ref3 */
	  if (sscanf(p,"%c",&interval_char) == 1) {
	    if (interval_char == 'x') {
	      found_interval = 6;
	    } else {
	      digit_string[0] = interval_char;
	      found_interval = atoi(digit_string);
	    }
	  } else {
	    abort();
	  }
#ifdef PMAP
	  found_index1part = found_interval;
#else
	  found_index1part = 12;
#endif

	} else {
	  fprintf(stderr,"Cannot parse part between %s and offsets in filename %s\n",idx_filesuffix,filename);
	  if (snps_root != NULL) {
	    FREE(offsets_suffix);
	    FREE(positions_high_suffix);
	    FREE(positions_low_suffix);
	  }
	  return (Filenames_T) NULL;
	}

	if (required_interval != 0) {
	  if (found_interval == required_interval) {
	    *index1part = found_index1part;
	    *index1interval = found_interval;
	    FREE(base_filename);
	    base_filename = (char *) CALLOC(strlen(filename)+1,sizeof(char));
	    strcpy(base_filename,filename);
	  }
	} else {
	  if (found_interval < *index1interval) {
	    *index1part = found_index1part;
	    *index1interval = found_interval;
	    FREE(base_filename);
	    base_filename = (char *) CALLOC(strlen(filename)+1,sizeof(char));
	    strcpy(base_filename,filename);
	  }
	}
      }
    }
  }

  FREE(pattern);

  if (closedir(dp) < 0) {
    fprintf(stderr,"Unable to close directory %s\n",genomesubdir);
  }

  /* Construct full filenames */
  if (base_filename == NULL) {
#if 0
    fprintf(stderr,"Cannot find offsets file containing %s and %s",idx_filesuffix,offsets_suffix);
    if (required_interval > 0) {
      fprintf(stderr," and having sampling interval of %d",required_interval);
    }
    fprintf(stderr,"\n");
#endif

    /* offsets_filename = (char *) NULL; */
    /* positions_high_filename = (char *) NULL; */
    /* positions_low_filename = (char *) NULL; */
    if (snps_root != NULL) {
      FREE(offsets_suffix);
      FREE(positions_high_suffix);
      FREE(positions_low_suffix);
    }
    return (Filenames_T) NULL;

  } else {
    offsets_filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+strlen(base_filename)+1,sizeof(char));
    offsets_basename_ptr = &(offsets_filename[strlen(genomesubdir)+strlen("/")]);
    offsets_index1info_ptr = &(offsets_basename_ptr[patternlength]);

    sprintf(offsets_filename,"%s/%s",genomesubdir,base_filename);
    if (Access_file_exists_p(offsets_filename) == false) {
      fprintf(stderr,"Offsets filename %s does not exist\n",offsets_filename);
      FREE(offsets_filename);
      /* offsets_filename = (char *) NULL; */
      /* positions_high_filename = (char *) NULL; */
      /* positions_low_filename = (char *) NULL; */
      FREE(base_filename);
      if (snps_root != NULL) {
	FREE(offsets_suffix);
	FREE(positions_high_suffix);
	FREE(positions_low_suffix);
      }
      return (Filenames_T) NULL;
    }


    if ((q = strstr(base_filename,offsets_suffix)) == NULL) {
      abort();
    } else {
      rootlength = q - base_filename;
    }

    positions_high_filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+rootlength+strlen(positions_high_suffix)+1,sizeof(char));
    positions_high_basename_ptr = &(positions_high_filename[strlen(genomesubdir)+strlen("/")]);
    positions_high_index1info_ptr = &(positions_high_basename_ptr[patternlength]);

    positions_low_filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+rootlength+strlen(positions_low_suffix)+1,sizeof(char));
    positions_low_basename_ptr = &(positions_low_filename[strlen(genomesubdir)+strlen("/")]);
    positions_low_index1info_ptr = &(positions_low_basename_ptr[patternlength]);

    sprintf(positions_high_filename,"%s/",genomesubdir);
    strncpy(positions_high_basename_ptr,base_filename,rootlength);
    strcpy(&(positions_high_basename_ptr[rootlength]),positions_high_suffix);

    sprintf(positions_low_filename,"%s/",genomesubdir);
    strncpy(positions_low_basename_ptr,base_filename,rootlength);
    strcpy(&(positions_low_basename_ptr[rootlength]),positions_low_suffix);

    if (offsets_only_p == true) {
      /* Do not look for a positions file */
    } else if (Access_file_exists_p(positions_low_filename) == false) {
      fprintf(stderr,"Positions filename %s does not exist\n",positions_low_filename);
      FREE(offsets_filename);
      FREE(positions_high_filename);
      FREE(positions_low_filename);
      /* offsets_filename = (char *) NULL; */
      /* positions_high_filename = (char *) NULL; */
      /* positions_low_filename = (char *) NULL; */
      FREE(base_filename);
      if (snps_root != NULL) {
	FREE(offsets_suffix);
	FREE(positions_high_suffix);
	FREE(positions_low_suffix);
      }
      return (Filenames_T) NULL;
    }

    if (Access_file_exists_p(positions_high_filename) == false) {
      /* Not a large genome */
      FREE(positions_high_filename);
      positions_high_filename = (char *) NULL;
    }

    if (snps_root != NULL) {
      FREE(offsets_suffix);
      FREE(positions_high_suffix);
      FREE(positions_low_suffix);
    }

    FREE(base_filename);
    fprintf(stderr,"Looking for index files in directory %s (offsets not compressed)\n",genomesubdir);
    fprintf(stderr,"  Offsets file is %s\n",offsets_basename_ptr);
    if (positions_high_filename == NULL) {
      fprintf(stderr,"  Positions file is %s\n",positions_low_basename_ptr);
    } else {
      fprintf(stderr,"  Positions files are %s and %s\n",positions_high_basename_ptr,positions_low_basename_ptr);
    }
    return Filenames_new(/*pages_filename*/NULL,/*pointers_filename*/NULL,offsets_filename,
			 positions_high_filename,positions_low_filename,
			 /*pointers_basename_ptr*/NULL,offsets_basename_ptr,
			 positions_high_basename_ptr,positions_low_basename_ptr,
			 /*pointers_index1info_ptr*/NULL,offsets_index1info_ptr,
			 positions_high_index1info_ptr,positions_low_index1info_ptr);
  }
}



#ifdef PMAP
#define BASE_KMER_SAMPLING 3   /* e.g., 677 */
#define KMER_SAMPLING 2   /* e.g., 77 */
#else
#define BASE_KMER_SAMPLING 5   /* e.g., 12153 */
#define KMER_SAMPLING 3   /* e.g., 153 */
#endif


Filenames_T
Indexdb_get_filenames_bitpack64 (
#ifdef PMAP
				 Alphabet_T *alphabet, Alphabet_T required_alphabet,
#endif
				 Width_T *basesize, Width_T *index1part, Width_T *index1interval,
				 char *genomesubdir, char *fileroot, char *idx_filesuffix, char *snps_root,
				 Width_T required_basesize, Width_T required_index1part, Width_T required_interval,
				 bool offsets_only_p) {
  char *bitpackpages_filename, *bitpackptrs_filename, *offsetscomp_filename,
    *positions_high_filename, *positions_low_filename,
    *bitpackpages_basename_ptr, *bitpackptrs_basename_ptr, *offsetscomp_basename_ptr,
    *positions_high_basename_ptr, *positions_low_basename_ptr,
    *bitpackpages_index1info_ptr, *bitpackptrs_index1info_ptr, *offsetscomp_index1info_ptr,
    *positions_high_index1info_ptr, *positions_low_index1info_ptr;

  char *base_filename, *filename;
#ifdef PMAP
  char *pattern1, *pattern2, *a;
  int patternlength1, patternlength2, alphabet_strlen;
  Alphabet_T found_alphabet;
#else
  char *pattern;
  char tens0, tens;
#endif
  char interval_char, digit_string[2], *p, *q;
  Width_T found_basesize = 0, found_index1part = 0, found_interval = 0;
  int rootlength, patternlength;

  char ones0, ones;
  char *bitpackpages_suffix, *bitpackptrs_suffix, *offsetscomp_suffix,
    *positions_high_suffix, *positions_low_suffix;
  struct dirent *entry;
  DIR *dp;


  if (snps_root == NULL) {
    bitpackpages_suffix = "bitpackpages";
    bitpackptrs_suffix = "bitpackptrs";
    offsetscomp_suffix = "bitpackcomp";
    positions_high_suffix = POSITIONS_HIGH_FILESUFFIX;
    positions_low_suffix = POSITIONS_LOW_FILESUFFIX;
  } else {
    bitpackpages_suffix = (char *) CALLOC(strlen("bitpackpages")+strlen(".")+strlen(snps_root)+1,sizeof(char));
    sprintf(bitpackpages_suffix,"%s.%s","bitpackpages",snps_root);
    bitpackptrs_suffix = (char *) CALLOC(strlen("bitpackptrs")+strlen(".")+strlen(snps_root)+1,sizeof(char));
    sprintf(bitpackptrs_suffix,"%s.%s","bitpackptrs",snps_root);
    offsetscomp_suffix = (char *) CALLOC(strlen("bitpackcomp")+strlen(".")+strlen(snps_root)+1,sizeof(char));
    sprintf(offsetscomp_suffix,"%s.%s","bitpackcomp",snps_root);
    positions_high_suffix = (char *) CALLOC(strlen(POSITIONS_HIGH_FILESUFFIX)+strlen(".")+strlen(snps_root)+1,sizeof(char));
    sprintf(positions_high_suffix,"%s.%s",POSITIONS_HIGH_FILESUFFIX,snps_root);
    positions_low_suffix = (char *) CALLOC(strlen(POSITIONS_LOW_FILESUFFIX)+strlen(".")+strlen(snps_root)+1,sizeof(char));
    sprintf(positions_low_suffix,"%s.%s",POSITIONS_LOW_FILESUFFIX,snps_root);
  }

#ifdef PMAP
  *alphabet = NALPHABETS + 1;
#endif
  *basesize = 0;
  *index1part = 0;
  *index1interval = 1000;
  base_filename = (char *) NULL;

  if ((dp = opendir(genomesubdir)) == NULL) {
    fprintf(stderr,"Unable to open directory %s\n",genomesubdir);
    exit(9);
  }

#ifdef PMAP
  pattern1 = (char *) CALLOC(strlen(fileroot)+strlen(".")+1,sizeof(char)); /* e.g., "hg19." */
  sprintf(pattern1,"%s.",fileroot);
  patternlength1 = strlen(pattern1);

  pattern2 = (char *) CALLOC(strlen(".")+strlen(idx_filesuffix)+1,sizeof(char)); /* e.g., ".pr" */
  sprintf(pattern2,".%s",idx_filesuffix);
  patternlength2 = strlen(pattern2);

  digit_string[1] = '\0';	/* Needed for atoi */
  while ((entry = readdir(dp)) != NULL) {
    filename = entry->d_name;
    if (!strncmp(filename,pattern1,patternlength1)) {
      a = &(filename[strlen(pattern1)]); /* Points after fileroot, e.g., "hg19." */
      if ((p = strstr(a,pattern2)) != NULL && (q = strstr(p,offsetscomp_suffix)) != NULL && !strcmp(q,offsetscomp_suffix)) {
	if ((found_alphabet = Alphabet_find(a)) != AA0) {
	  alphabet_strlen = p - a;
	  p += patternlength2;

	  if (q - p == KMER_SAMPLING) {
	    /* Latest style, e.g., pf77 */
	    if (sscanf(p,"%c%c",&ones,&interval_char) == 2) {
	      digit_string[0] = ones;
	      found_index1part = atoi(digit_string);

	      digit_string[0] = interval_char;
	      found_interval = atoi(digit_string);
	    } else {
	      abort();
	    }

	  } else if (q - p == BASE_KMER_SAMPLING) {
	    /* Previous style, e.g., pf677 */
	    if (sscanf(p,"%c%c%c",&ones0,&ones,&interval_char) == 3) {
	      digit_string[0] = ones0;
	      found_basesize = atoi(digit_string);

	      digit_string[0] = ones;
	      found_index1part = atoi(digit_string);

	      digit_string[0] = interval_char;
	      found_interval = atoi(digit_string);
	    } else {
	      abort();
	    }

	  } else {
	    /* fprintf(stderr,"Cannot parse part between %s and offsets in filename %s\n",idx_filesuffix,filename); */
	    if (snps_root != NULL) {
	      FREE(offsetscomp_suffix);
	      FREE(bitpackptrs_suffix);
	      FREE(bitpackpages_suffix);
	      FREE(positions_high_suffix);
	      FREE(positions_low_suffix);
	    }
	    return (Filenames_T) NULL;
	  }

	  if ((required_alphabet == AA0 || found_alphabet == required_alphabet) &&
	      (required_index1part == 0 || found_index1part == required_index1part) &&
	      (required_interval == 0 || found_interval == required_interval)) {
	    if (required_alphabet == AA0 && found_alphabet > *alphabet) {
	      /* Skip, since we have already found an earlier alphabet */
	    } else if (required_index1part == 0 && found_index1part < *index1part) {
	      /* Skip, since we have already found a larger index1part */
	    } else if (required_interval == 0 && found_interval > *index1interval) {
	      /* Skip, since we have already found a smaller interval */
	    } else {
	      patternlength = patternlength1 + alphabet_strlen + patternlength2;
	      *basesize = found_basesize;
	      *index1part = found_index1part;
	      *index1interval = found_interval;
	      *alphabet = found_alphabet;
	      FREE(base_filename);
	      base_filename = (char *) CALLOC(strlen(filename)+1,sizeof(char));
	      strcpy(base_filename,filename);
	    }
	  }
	}
      }
    }
  }

  FREE(pattern2);
  FREE(pattern1);

#else

  pattern = (char *) CALLOC(strlen(fileroot)+strlen(".")+strlen(idx_filesuffix)+1,sizeof(char));
  sprintf(pattern,"%s.%s",fileroot,idx_filesuffix);
  patternlength = strlen(pattern);

  digit_string[1] = '\0';	/* Needed for atoi */
  while ((entry = readdir(dp)) != NULL) {
    filename = entry->d_name;
    if (!strncmp(filename,pattern,patternlength)) {
      p = &(filename[strlen(pattern)]); /* Points after idx_filesuffix, e.g., "ref" */
      if ((q = strstr(p,offsetscomp_suffix)) != NULL && !strcmp(q,offsetscomp_suffix)) {

	if (q - p == BASE_KMER_SAMPLING) {
	  /* New style, e.g., ref12153 */
	  if (sscanf(p,"%c%c%c%c%c",&tens0,&ones0,&tens,&ones,&interval_char) == 5) {
	    digit_string[0] = tens0;
	    found_basesize = 10*atoi(digit_string);
	    digit_string[0] = ones0;
	    found_basesize += atoi(digit_string);

	    digit_string[0] = tens;
	    found_index1part = 10*atoi(digit_string);
	    digit_string[0] = ones;
	    found_index1part += atoi(digit_string);

	    digit_string[0] = interval_char;
	    found_interval = atoi(digit_string);
	  } else {
	    abort();
	  }

	} else {
	  fprintf(stderr,"Cannot parse part between %s and offsets in filename %s: found %ld characters, expecting %d\n",
		  idx_filesuffix,filename,q-p,BASE_KMER_SAMPLING);
	  if (snps_root != NULL) {
	    FREE(offsetscomp_suffix);
	    FREE(bitpackptrs_suffix);
	    FREE(bitpackpages_suffix);
	    FREE(positions_high_suffix);
	    FREE(positions_low_suffix);
	  }
	  return (Filenames_T) NULL;
	}

	if ((required_index1part == 0 || found_index1part == required_index1part) &&
	    (required_basesize == 0 || found_basesize == required_basesize) &&
	    (required_interval == 0 || found_interval == required_interval)) {
	  if (required_index1part == 0 && found_index1part < *index1part) {
	    /* Skip, since we have already found a larger index1part */
	  } else if (required_basesize == 0 && found_basesize < *basesize) {
	    /* Skip, since we have already found a larger basesize */
	  } else if (required_interval == 0 && found_interval > *index1interval) {
	    /* Skip, since we have already found a smaller interval */
	  } else {
#ifdef PMAP
	    *basesize = 0;
#else
	    *basesize = found_basesize;
#endif
	    *index1part = found_index1part;
	    *index1interval = found_interval;
	    FREE(base_filename);
	    base_filename = (char *) CALLOC(strlen(filename)+1,sizeof(char));
	    strcpy(base_filename,filename);
	  }
	}
      }
    }
  }

  FREE(pattern);
#endif


  if (closedir(dp) < 0) {
    fprintf(stderr,"Unable to close directory %s\n",genomesubdir);
  }

  /* Construct full filenames */
  if (base_filename == NULL) {
#if 0
    fprintf(stderr,"Cannot find offsets file containing %s and %s",idx_filesuffix,offsetscomp_suffix);
#ifdef PMAP
    if (required_alphabet > AA0) {
      fprintf(stderr," and having alphabet %s",Alphabet_string(required_alphabet));
    }
#endif
    if (required_index1part > 0) {
      fprintf(stderr," and having k-mer of %d",required_index1part);
    }
    if (required_interval > 0) {
      fprintf(stderr," and having sampling interval of %d",required_interval);
    }
    fprintf(stderr,"\n");
#endif

    /* bitpackpages_filename = (char *) NULL; */
    /* bitpackptrs_filename = (char *) NULL; */
    /* offsetscomp_filename = (char *) NULL; */
    /* positions_high_filename = (char *) NULL; */
    /* positions_low_filename = (char *) NULL; */
    if (snps_root != NULL) {
      FREE(offsetscomp_suffix);
      FREE(bitpackptrs_suffix);
      FREE(bitpackpages_suffix);
      FREE(positions_high_suffix);
      FREE(positions_low_suffix);
    }
    return (Filenames_T) NULL;

  } else {
    offsetscomp_filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+strlen(base_filename)+1,sizeof(char));
    offsetscomp_basename_ptr = &(offsetscomp_filename[strlen(genomesubdir)+strlen("/")]);
    offsetscomp_index1info_ptr = &(offsetscomp_basename_ptr[patternlength]);

    sprintf(offsetscomp_filename,"%s/%s",genomesubdir,base_filename);
    if (Access_file_exists_p(offsetscomp_filename) == false) {
      fprintf(stderr,"Offsets filename %s does not exist\n",offsetscomp_filename);
      FREE(offsetscomp_filename);
      /* offsetscomp_filename = (char *) NULL; */
      /* positions_high_filename = (char *) NULL; */
      /* positions_low_filename = (char *) NULL; */
      FREE(base_filename);
      if (snps_root != NULL) {
	FREE(offsetscomp_suffix);
	FREE(bitpackptrs_suffix);
	FREE(bitpackpages_suffix);
	FREE(positions_high_suffix);
	FREE(positions_low_suffix);
      }
      return (Filenames_T) NULL;
    }


    if ((q = strstr(base_filename,offsetscomp_suffix)) == NULL) {
      abort();
    } else {
      rootlength = q - base_filename;
    }

    if (*index1part == *basesize) {
      /* bitpackpages_filename = (char *) NULL; */
      /* bitpackptrs_filename = (char *) NULL; */
      /* bitpackptrs_basename_ptr = (char *) NULL; */
      /* bitpackptrs_index1info_ptr = (char *) NULL; */
    } else {
      bitpackptrs_filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+rootlength+strlen(bitpackptrs_suffix)+1,sizeof(char));
      bitpackptrs_basename_ptr = &(bitpackptrs_filename[strlen(genomesubdir)+strlen("/")]);
      bitpackptrs_index1info_ptr = &(bitpackptrs_basename_ptr[patternlength]);

      sprintf(bitpackptrs_filename,"%s/",genomesubdir);
      strncpy(bitpackptrs_basename_ptr,base_filename,rootlength);
      strcpy(&(bitpackptrs_basename_ptr[rootlength]),bitpackptrs_suffix);

      if (Access_file_exists_p(bitpackptrs_filename) == false) {
	fprintf(stderr,"Bitpackptrs filename %s does not exist\n",bitpackptrs_filename);
	FREE(offsetscomp_filename);
	/* offsetscomp_filename = (char *) NULL; */
	FREE(base_filename);
	if (snps_root != NULL) {
	  FREE(offsetscomp_suffix);
	  FREE(bitpackptrs_suffix);
	  FREE(positions_high_suffix);
	  FREE(positions_low_suffix);
	}
	return (Filenames_T) NULL;
      }
    }

    bitpackpages_filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+rootlength+strlen(bitpackpages_suffix)+1,sizeof(char));
    bitpackpages_basename_ptr = &(bitpackpages_filename[strlen(genomesubdir)+strlen("/")]);
    bitpackpages_index1info_ptr = &(bitpackpages_basename_ptr[patternlength]);

    sprintf(bitpackpages_filename,"%s/",genomesubdir);
    strncpy(bitpackpages_basename_ptr,base_filename,rootlength);
    strcpy(&(bitpackpages_basename_ptr[rootlength]),bitpackpages_suffix);
    if (Access_file_exists_p(bitpackpages_filename) == false) {
      /* Not a huge genome */
      FREE(bitpackpages_filename);
      bitpackpages_filename = (char *) NULL;
    }

    positions_high_filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+rootlength+strlen(positions_high_suffix)+1,sizeof(char));
    positions_high_basename_ptr = &(positions_high_filename[strlen(genomesubdir)+strlen("/")]);
    positions_high_index1info_ptr = &(positions_high_basename_ptr[patternlength]);

    positions_low_filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+rootlength+strlen(positions_low_suffix)+1,sizeof(char));
    positions_low_basename_ptr = &(positions_low_filename[strlen(genomesubdir)+strlen("/")]);
    positions_low_index1info_ptr = &(positions_low_basename_ptr[patternlength]);

    sprintf(positions_high_filename,"%s/",genomesubdir);
    strncpy(positions_high_basename_ptr,base_filename,rootlength);
    strcpy(&(positions_high_basename_ptr[rootlength]),positions_high_suffix);

    sprintf(positions_low_filename,"%s/",genomesubdir);
    strncpy(positions_low_basename_ptr,base_filename,rootlength);
    strcpy(&(positions_low_basename_ptr[rootlength]),positions_low_suffix);

    if (offsets_only_p == true) {
      /* Do not look for a positions file */
    } else if (Access_file_exists_p(positions_low_filename) == false) {
      /* Try newer naming scheme: ref153positions instead of ref12153positions */
      sprintf(positions_high_filename,"%s/",genomesubdir);
      strncpy(positions_high_basename_ptr,base_filename,rootlength-BASE_KMER_SAMPLING); /* e.g., skip "12153" */
      strncpy(&(positions_high_basename_ptr[rootlength-BASE_KMER_SAMPLING]),&(base_filename[rootlength-KMER_SAMPLING]),KMER_SAMPLING);
      strcpy(&(positions_high_basename_ptr[rootlength+KMER_SAMPLING-BASE_KMER_SAMPLING]),positions_high_suffix);

      sprintf(positions_low_filename,"%s/",genomesubdir);
      strncpy(positions_low_basename_ptr,base_filename,rootlength-BASE_KMER_SAMPLING); /* e.g., skip "12153" */
      strncpy(&(positions_low_basename_ptr[rootlength-BASE_KMER_SAMPLING]),&(base_filename[rootlength-KMER_SAMPLING]),KMER_SAMPLING);
      strcpy(&(positions_low_basename_ptr[rootlength+KMER_SAMPLING-BASE_KMER_SAMPLING]),positions_low_suffix);

      if (Access_file_exists_p(positions_low_filename) == false) {
	fprintf(stderr,"Positions filename %s does not exist\n",positions_low_filename);
	FREE(bitpackpages_filename);
	FREE(bitpackptrs_filename);
	FREE(offsetscomp_filename);
	FREE(positions_high_filename);
	FREE(positions_low_filename);
	/* bitpackptrs_filename = (char *) NULL; */
	/* offsetscomp_filename = (char *) NULL; */
	/* positions_high_filename = (char *) NULL; */
	/* positions_low_filename = (char *) NULL; */
	FREE(base_filename);
	if (snps_root != NULL) {
	  FREE(offsetscomp_suffix);
	  FREE(bitpackptrs_suffix);
	  FREE(bitpackpages_suffix);
	  FREE(positions_high_suffix);
	  FREE(positions_low_suffix);
	}
	return (Filenames_T) NULL;
      }
    }

    if (Access_file_exists_p(positions_high_filename) == false) {
      /* Not a large genome */
      FREE(positions_high_filename);
      positions_high_filename = (char *) NULL;
    }

    if (snps_root != NULL) {
      FREE(offsetscomp_suffix);
      FREE(bitpackptrs_suffix);
      FREE(bitpackpages_suffix);
      FREE(positions_high_suffix);
      FREE(positions_low_suffix);
    }

    FREE(base_filename);

    fprintf(stderr,"Looking for index files in directory %s\n",genomesubdir);
    if (*basesize == *index1part) {
      /* No pointers file */
      bitpackpages_filename = bitpackptrs_filename = bitpackptrs_basename_ptr = bitpackptrs_index1info_ptr = (char *) NULL;
    } else {
      if (bitpackpages_filename != NULL)  {
	fprintf(stderr,"  Pages file is %s\n",bitpackpages_basename_ptr);
      }
      fprintf(stderr,"  Pointers file is %s\n",bitpackptrs_basename_ptr);
    }
    fprintf(stderr,"  Offsets file is %s\n",offsetscomp_basename_ptr);
    if (positions_high_filename == NULL) {
      fprintf(stderr,"  Positions file is %s\n",positions_low_basename_ptr);
    } else {
      fprintf(stderr,"  Positions files are %s and %s\n",positions_high_basename_ptr,positions_low_basename_ptr);
    }
    return Filenames_new(bitpackpages_filename,bitpackptrs_filename,offsetscomp_filename,
			 positions_high_filename,positions_low_filename,
			 bitpackptrs_basename_ptr,offsetscomp_basename_ptr,
			 positions_high_basename_ptr,positions_low_basename_ptr,
			 bitpackptrs_index1info_ptr,offsetscomp_index1info_ptr,
			 positions_high_index1info_ptr,positions_low_index1info_ptr);

  }
}



Filenames_T
Indexdb_get_filenames_gamma (
#ifdef PMAP
			     Alphabet_T *alphabet, Alphabet_T required_alphabet,
#endif
			     Width_T *basesize, Width_T *index1part, Width_T *index1interval,
			     char *genomesubdir, char *fileroot, char *idx_filesuffix, char *snps_root,
			     Width_T required_basesize, Width_T required_index1part, Width_T required_interval,
			     bool offsets_only_p) {
  char *gammaptrs_filename, *offsetscomp_filename, *positions_high_filename, *positions_low_filename,
    *gammaptrs_basename_ptr, *offsetscomp_basename_ptr,
    *positions_high_basename_ptr, *positions_low_basename_ptr,
    *gammaptrs_index1info_ptr, *offsetscomp_index1info_ptr,
    *positions_high_index1info_ptr, *positions_low_index1info_ptr;

  char *base_filename, *filename;
#ifdef PMAP
  char *pattern1, *pattern2, *a;
  int patternlength1, patternlength2, alphabet_strlen;
  Alphabet_T found_alphabet;
#else
  char *pattern;
  char tens0, tens;
#endif
  char interval_char, digit_string[2], *p, *q;
  Width_T found_basesize = 0, found_index1part = 0, found_interval = 0;
  int rootlength, patternlength;

  char ones0, ones;
  char *gammaptrs_suffix, *offsetscomp_suffix, *positions_high_suffix, *positions_low_suffix;
  struct dirent *entry;
  DIR *dp;


  if (snps_root == NULL) {
    gammaptrs_suffix = "gammaptrs";
    offsetscomp_suffix = "offsetscomp";
    positions_high_suffix = POSITIONS_HIGH_FILESUFFIX;
    positions_low_suffix = POSITIONS_LOW_FILESUFFIX;
  } else {
    gammaptrs_suffix = (char *) CALLOC(strlen("gammaptrs")+strlen(".")+strlen(snps_root)+1,sizeof(char));
    sprintf(gammaptrs_suffix,"%s.%s","gammaptrs",snps_root);
    offsetscomp_suffix = (char *) CALLOC(strlen("offsetscomp")+strlen(".")+strlen(snps_root)+1,sizeof(char));
    sprintf(offsetscomp_suffix,"%s.%s","offsetscomp",snps_root);
    positions_high_suffix = (char *) CALLOC(strlen(POSITIONS_HIGH_FILESUFFIX)+strlen(".")+strlen(snps_root)+1,sizeof(char));
    sprintf(positions_high_suffix,"%s.%s",POSITIONS_HIGH_FILESUFFIX,snps_root);
    positions_low_suffix = (char *) CALLOC(strlen(POSITIONS_LOW_FILESUFFIX)+strlen(".")+strlen(snps_root)+1,sizeof(char));
    sprintf(positions_low_suffix,"%s.%s",POSITIONS_LOW_FILESUFFIX,snps_root);
  }

#ifdef PMAP
  *alphabet = NALPHABETS + 1;
#endif
  *basesize = 0;
  *index1part = 0;
  *index1interval = 1000;
  base_filename = (char *) NULL;

  if ((dp = opendir(genomesubdir)) == NULL) {
    fprintf(stderr,"Unable to open directory %s\n",genomesubdir);
    exit(9);
  }

#ifdef PMAP
  pattern1 = (char *) CALLOC(strlen(fileroot)+strlen(".")+1,sizeof(char)); /* e.g., "hg19." */
  sprintf(pattern1,"%s.",fileroot);
  patternlength1 = strlen(pattern1);

  pattern2 = (char *) CALLOC(strlen(".")+strlen(idx_filesuffix)+1,sizeof(char)); /* e.g., ".pr" */
  sprintf(pattern2,".%s",idx_filesuffix);
  patternlength2 = strlen(pattern2);

  digit_string[1] = '\0';	/* Needed for atoi */
  while ((entry = readdir(dp)) != NULL) {
    filename = entry->d_name;
    if (!strncmp(filename,pattern1,patternlength1)) {
      a = &(filename[strlen(pattern1)]); /* Points after fileroot, e.g., "hg19." */
      if ((p = strstr(a,pattern2)) != NULL && (q = strstr(p,offsetscomp_suffix)) != NULL && !strcmp(q,offsetscomp_suffix)) {
	if ((found_alphabet = Alphabet_find(a)) != AA0) {
	  alphabet_strlen = p - a;
	  p += patternlength2;

	  if (q - p == BASE_KMER_SAMPLING) {
	    /* New style, e.g., pf677 */
	    if (sscanf(p,"%c%c%c",&ones0,&ones,&interval_char) == 3) {
	      digit_string[0] = ones0;
	      found_basesize = atoi(digit_string);

	      digit_string[0] = ones;
	      found_index1part = atoi(digit_string);

	      digit_string[0] = interval_char;
	      found_interval = atoi(digit_string);
	    } else {
	      abort();
	    }

	  } else {
	    /* fprintf(stderr,"Cannot parse part between %s and offsets in filename %s\n",idx_filesuffix,filename); */
	    if (snps_root != NULL) {
	      FREE(offsetscomp_suffix);
	      FREE(gammaptrs_suffix);
	      FREE(positions_high_suffix);
	      FREE(positions_low_suffix);
	    }
	    return (Filenames_T) NULL;
	  }

	  if ((required_alphabet == AA0 || found_alphabet == required_alphabet) &&
	      (required_index1part == 0 || found_index1part == required_index1part) &&
	      (required_basesize == 0 || found_basesize == required_basesize) &&
	      (required_interval == 0 || found_interval == required_interval)) {
	    if (required_alphabet == AA0 && found_alphabet > *alphabet) {
	      /* Skip, since we have already found an earlier alphabet */
	    } else if (required_index1part == 0 && found_index1part < *index1part) {
	      /* Skip, since we have already found a larger index1part */
	    } else if (required_basesize == 0 && found_basesize < *basesize) {
	      /* Skip, since we have already found a larger basesize */
	    } else if (required_interval == 0 && found_interval > *index1interval) {
	      /* Skip, since we have already found a smaller interval */
	    } else {
	      patternlength = patternlength1 + alphabet_strlen + patternlength2;
	      *basesize = found_basesize;
	      *index1part = found_index1part;
	      *index1interval = found_interval;
	      *alphabet = found_alphabet;
	      FREE(base_filename);
	      base_filename = (char *) CALLOC(strlen(filename)+1,sizeof(char));
	      strcpy(base_filename,filename);
	    }
	  }
	}
      }
    }
  }

  FREE(pattern2);
  FREE(pattern1);

#else

  pattern = (char *) CALLOC(strlen(fileroot)+strlen(".")+strlen(idx_filesuffix)+1,sizeof(char));
  sprintf(pattern,"%s.%s",fileroot,idx_filesuffix);
  patternlength = strlen(pattern);

  digit_string[1] = '\0';	/* Needed for atoi */
  while ((entry = readdir(dp)) != NULL) {
    filename = entry->d_name;
    if (!strncmp(filename,pattern,patternlength)) {
      p = &(filename[strlen(pattern)]); /* Points after idx_filesuffix, e.g., "ref" */
      if ((q = strstr(p,offsetscomp_suffix)) != NULL && !strcmp(q,offsetscomp_suffix)) {

	if (q - p == BASE_KMER_SAMPLING) {
	  /* New style, e.g., ref12153 */
	  if (sscanf(p,"%c%c%c%c%c",&tens0,&ones0,&tens,&ones,&interval_char) == 5) {
	    digit_string[0] = tens0;
	    found_basesize = 10*atoi(digit_string);
	    digit_string[0] = ones0;
	    found_basesize += atoi(digit_string);

	    digit_string[0] = tens;
	    found_index1part = 10*atoi(digit_string);
	    digit_string[0] = ones;
	    found_index1part += atoi(digit_string);

	    digit_string[0] = interval_char;
	    found_interval = atoi(digit_string);
	  } else {
	    abort();
	  }

	} else {
	  fprintf(stderr,"Cannot parse part between %s and offsets in filename %s: found %ld characters, expecting %d\n",
		  idx_filesuffix,filename,q-p,BASE_KMER_SAMPLING);
	  if (snps_root != NULL) {
	    FREE(offsetscomp_suffix);
	    FREE(gammaptrs_suffix);
	    FREE(positions_high_suffix);
	    FREE(positions_low_suffix);
	  }
	  return (Filenames_T) NULL;
	}

	if ((required_index1part == 0 || found_index1part == required_index1part) &&
	    (required_basesize == 0 || found_basesize == required_basesize) &&
	    (required_interval == 0 || found_interval == required_interval)) {
	  if (required_index1part == 0 && found_index1part < *index1part) {
	    /* Skip, since we have already found a larger index1part */
	  } else if (required_basesize == 0 && found_basesize < *basesize) {
	    /* Skip, since we have already found a larger basesize */
	  } else if (required_interval == 0 && found_interval > *index1interval) {
	    /* Skip, since we have already found a smaller interval */
	  } else {
	    *basesize = found_basesize;
	    *index1part = found_index1part;
	    *index1interval = found_interval;
	    FREE(base_filename);
	    base_filename = (char *) CALLOC(strlen(filename)+1,sizeof(char));
	    strcpy(base_filename,filename);
	  }
	}
      }
    }
  }

  FREE(pattern);
#endif


  if (closedir(dp) < 0) {
    fprintf(stderr,"Unable to close directory %s\n",genomesubdir);
  }

  /* Construct full filenames */
  if (base_filename == NULL) {
#if 0
    fprintf(stderr,"Cannot find offsets file containing %s and %s",idx_filesuffix,offsetscomp_suffix);
#ifdef PMAP
    if (required_alphabet > AA0) {
      fprintf(stderr," and having alphabet %s",Alphabet_string(required_alphabet));
    }
#endif
    if (required_index1part > 0) {
      fprintf(stderr," and having k-mer of %d",required_index1part);
    }
    if (required_interval > 0) {
      fprintf(stderr," and having sampling interval of %d",required_interval);
    }
    fprintf(stderr,"\n");
#endif

    /* gammaptrs_filename = (char *) NULL; */
    /* offsetscomp_filename = (char *) NULL; */
    /* positions_high_filename = (char *) NULL; */
    /* positions_low_filename = (char *) NULL; */
    if (snps_root != NULL) {
      FREE(offsetscomp_suffix);
      FREE(gammaptrs_suffix);
      FREE(positions_high_suffix);
      FREE(positions_low_suffix);
    }
    return (Filenames_T) NULL;

  } else {
    offsetscomp_filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+strlen(base_filename)+1,sizeof(char));
    offsetscomp_basename_ptr = &(offsetscomp_filename[strlen(genomesubdir)+strlen("/")]);
    offsetscomp_index1info_ptr = &(offsetscomp_basename_ptr[patternlength]);

    sprintf(offsetscomp_filename,"%s/%s",genomesubdir,base_filename);
    if (Access_file_exists_p(offsetscomp_filename) == false) {
      fprintf(stderr,"Offsets filename %s does not exist\n",offsetscomp_filename);
      FREE(offsetscomp_filename);
      /* offsetscomp_filename = (char *) NULL; */
      /* positions_high_filename = (char *) NULL; */
      /* positions_low_filename = (char *) NULL; */
      FREE(base_filename);
      if (snps_root != NULL) {
	FREE(offsetscomp_suffix);
	FREE(gammaptrs_suffix);
	FREE(positions_high_suffix);
	FREE(positions_low_suffix);
      }
      return (Filenames_T) NULL;
    }


    if ((q = strstr(base_filename,offsetscomp_suffix)) == NULL) {
      abort();
    } else {
      rootlength = q - base_filename;
    }

    if (*index1part == *basesize) {
      /* gammaptrs_filename = (char *) NULL; */
      /* gammaptrs_basename_ptr = (char *) NULL; */
      /* gammaptrs_index1info_ptr = (char *) NULL; */
    } else {
      gammaptrs_filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+rootlength+strlen(gammaptrs_suffix)+1,sizeof(char));
      gammaptrs_basename_ptr = &(gammaptrs_filename[strlen(genomesubdir)+strlen("/")]);
      gammaptrs_index1info_ptr = &(gammaptrs_basename_ptr[patternlength]);

      sprintf(gammaptrs_filename,"%s/",genomesubdir);
      strncpy(gammaptrs_basename_ptr,base_filename,rootlength);
      strcpy(&(gammaptrs_basename_ptr[rootlength]),gammaptrs_suffix);

      if (Access_file_exists_p(gammaptrs_filename) == false) {
	fprintf(stderr,"Gammaptrs filename %s does not exist\n",gammaptrs_filename);
	FREE(offsetscomp_filename);
	/* offsetscomp_filename = (char *) NULL; */
	FREE(base_filename);
	if (snps_root != NULL) {
	  FREE(offsetscomp_suffix);
	  FREE(gammaptrs_suffix);
	  FREE(positions_high_suffix);
	  FREE(positions_low_suffix);
	}
	return (Filenames_T) NULL;
      }
    }


    positions_high_filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+rootlength+strlen(positions_high_suffix)+1,sizeof(char));
    positions_high_basename_ptr = &(positions_high_filename[strlen(genomesubdir)+strlen("/")]);
    positions_high_index1info_ptr = &(positions_high_basename_ptr[patternlength]);

    positions_low_filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+rootlength+strlen(positions_low_suffix)+1,sizeof(char));
    positions_low_basename_ptr = &(positions_low_filename[strlen(genomesubdir)+strlen("/")]);
    positions_low_index1info_ptr = &(positions_low_basename_ptr[patternlength]);

    sprintf(positions_high_filename,"%s/",genomesubdir);
    strncpy(positions_high_basename_ptr,base_filename,rootlength);
    strcpy(&(positions_high_basename_ptr[rootlength]),positions_high_suffix);

    sprintf(positions_low_filename,"%s/",genomesubdir);
    strncpy(positions_low_basename_ptr,base_filename,rootlength);
    strcpy(&(positions_low_basename_ptr[rootlength]),positions_low_suffix);

    if (offsets_only_p == true) {
      /* Do not look for a positions file */
    } else if (Access_file_exists_p(positions_low_filename) == false) {
      /* Try newer naming scheme: ref153positions instead of ref12153positions */
      sprintf(positions_high_filename,"%s/",genomesubdir);
      strncpy(positions_high_basename_ptr,base_filename,rootlength-BASE_KMER_SAMPLING); /* e.g., skip "12153" */
      strncpy(&(positions_high_basename_ptr[rootlength-BASE_KMER_SAMPLING]),&(base_filename[rootlength-KMER_SAMPLING]),KMER_SAMPLING);
      strcpy(&(positions_high_basename_ptr[rootlength+KMER_SAMPLING-BASE_KMER_SAMPLING]),positions_high_suffix);

      sprintf(positions_low_filename,"%s/",genomesubdir);
      strncpy(positions_low_basename_ptr,base_filename,rootlength-BASE_KMER_SAMPLING); /* e.g., skip "12153" */
      strncpy(&(positions_low_basename_ptr[rootlength-BASE_KMER_SAMPLING]),&(base_filename[rootlength-KMER_SAMPLING]),KMER_SAMPLING);
      strcpy(&(positions_low_basename_ptr[rootlength+KMER_SAMPLING-BASE_KMER_SAMPLING]),positions_low_suffix);

      if (Access_file_exists_p(positions_low_filename) == false) {
	fprintf(stderr,"Positions filename %s does not exist\n",positions_low_filename);
	FREE(gammaptrs_filename);
	FREE(offsetscomp_filename);
	FREE(positions_high_filename);
	FREE(positions_low_filename);
	/* gammaptrs_filename = (char *) NULL; */
	/* offsetscomp_filename = (char *) NULL; */
	/* positions_high_filename = (char *) NULL; */
	/* positions_low_filename = (char *) NULL; */
	FREE(base_filename);
	if (snps_root != NULL) {
	  FREE(offsetscomp_suffix);
	  FREE(gammaptrs_suffix);
	  FREE(positions_high_suffix);
	  FREE(positions_low_suffix);
	}
	return (Filenames_T) NULL;
      }
    }

    if (Access_file_exists_p(positions_high_filename) == false) {
      /* Not a large genome */
      FREE(positions_high_filename);
      positions_high_filename = (char *) NULL;
    }

    if (snps_root != NULL) {
      FREE(offsetscomp_suffix);
      FREE(gammaptrs_suffix);
      FREE(positions_high_suffix);
      FREE(positions_low_suffix);
    }

    FREE(base_filename);

    fprintf(stderr,"Looking for index files in directory %s\n",genomesubdir);
    if (*basesize == *index1part) {
      /* No pointers file */
      gammaptrs_filename = gammaptrs_basename_ptr = gammaptrs_index1info_ptr = (char *) NULL;
    } else {
      fprintf(stderr,"  Pointers file is %s\n",gammaptrs_basename_ptr);
    }
    fprintf(stderr,"  Offsets file is %s\n",offsetscomp_basename_ptr);
    if (positions_high_filename == NULL) {
      fprintf(stderr,"  Positions file is %s\n",positions_low_basename_ptr);
    } else {
      fprintf(stderr,"  Positions files are %s and %s\n",positions_high_basename_ptr,positions_low_basename_ptr);
    }
    return Filenames_new(/*pages_filename*/NULL,gammaptrs_filename,offsetscomp_filename,
			 positions_high_filename,positions_low_filename,
			 gammaptrs_basename_ptr,offsetscomp_basename_ptr,
			 positions_high_basename_ptr,positions_low_basename_ptr,
			 gammaptrs_index1info_ptr,offsetscomp_index1info_ptr,
			 positions_high_index1info_ptr,positions_low_index1info_ptr);

  }
}



Filenames_T
Indexdb_get_filenames (int *compression_type,
#ifdef PMAP
		       Alphabet_T *alphabet, Alphabet_T required_alphabet,
#endif
		       Width_T *basesize, Width_T *index1part, Width_T *index1interval, char *genomesubdir,
		       char *fileroot, char *idx_filesuffix, char *snps_root,
		       Width_T required_basesize, Width_T required_index1part, Width_T required_interval,
		       bool offsets_only_p) {
  Filenames_T filenames;

  if ((filenames = Indexdb_get_filenames_no_compression(&(*index1part),&(*index1interval),
							genomesubdir,fileroot,idx_filesuffix,snps_root,
							required_interval,offsets_only_p)) != NULL) {
    *compression_type = NO_COMPRESSION;
    *basesize = *index1part;
    return filenames;
    

  } else if ((filenames = Indexdb_get_filenames_bitpack64(
#ifdef PMAP
							  &(*alphabet),required_alphabet,
#endif
							  &(*basesize),&(*index1part),&(*index1interval),
							  genomesubdir,fileroot,idx_filesuffix,snps_root,
							  required_basesize,required_index1part,required_interval,
							  offsets_only_p)) != NULL) {
    *compression_type = BITPACK64_COMPRESSION;
    return filenames;
    
  } else if ((filenames = Indexdb_get_filenames_gamma(
#ifdef PMAP
						      &(*alphabet),required_alphabet,
#endif
						      &(*basesize),&(*index1part),&(*index1interval),
						      genomesubdir,fileroot,idx_filesuffix,snps_root,
						      required_basesize,required_index1part,required_interval,
						      offsets_only_p)) != NULL) {
    *compression_type = GAMMA_COMPRESSION;
    return filenames;

  } else {
    return (Filenames_T) NULL;
  }
}


T
Indexdb_new_genome (Width_T *basesize, Width_T *index1part, Width_T *index1interval,
		    char *genomesubdir, char *fileroot, char *idx_filesuffix, char *snps_root,
#ifdef PMAP
		    Alphabet_T *alphabet, int *alphabet_size, Alphabet_T required_alphabet,
#endif
		    Width_T required_basesize, Width_T required_index1part, Width_T required_interval,
		    bool expand_offsets_p, Access_mode_T offsetscomp_access, Access_mode_T positions_access) {
  T new = (T) MALLOC(sizeof(*new));
  Filenames_T filenames;
  Oligospace_T basespace, base;

  unsigned int poly_T;
  Positionsptr_T ptr0, end0;
  off_t filesize;

#ifdef LARGE_GENOMES
  size_t bitpackpages_len;
#endif

  char *comma;
  double seconds;
#ifdef HAVE_MMAP
  int npages;
#endif


  if ((filenames = Indexdb_get_filenames_no_compression(&new->index1part,&new->index1interval,
							genomesubdir,fileroot,idx_filesuffix,snps_root,
							required_interval,/*offsets_only_p*/false)) != NULL) {
    /* Try non-compressed files */
    fprintf(stderr,"Offsets compression type: none\n");
    new->compression_type = NO_COMPRESSION;
    new->offsetscomp_basesize = new->index1part;
    new->offsetscomp_blocksize = 1;

    *basesize = new->index1part;
    *index1part = new->index1part;
    *index1interval = new->index1interval;

#ifdef PMAP
    basespace = power(*alphabet_size,new->index1part);
#else
    basespace = power(4,new->index1part);
#endif
    new->gammaptrs = (Gammaptr_T *) CALLOC(basespace+1,sizeof(Gammaptr_T));
    for (base = 0; base <= basespace; base++) {
      new->gammaptrs[base] = base;
    }


    if (offsetscomp_access == USE_ALLOCATE) {
      if (snps_root) {
	fprintf(stderr,"Allocating memory for %s (%s) offsets, kmer %d, interval %d...",
		idx_filesuffix,snps_root,new->index1part,new->index1interval);
      } else {
	fprintf(stderr,"Allocating memory for %s offsets, kmer %d, interval %d...",
		idx_filesuffix,new->index1part,new->index1interval);
      }
      new->offsetscomp = (Offsetscomp_T *) Access_allocated(&new->offsetscomp_len,&seconds,
							    filenames->offsets_filename,sizeof(Offsetscomp_T));
      if (new->offsetscomp == NULL) {
	fprintf(stderr,"insufficient memory (need to use a lower batch mode (-B))\n");
	exit(9);
      } else {
	comma = Genomicpos_commafmt(new->offsetscomp_len);
	fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma,seconds);
	FREE(comma);
	new->offsetscomp_access = ALLOCATED;
      }

#ifdef HAVE_MMAP
    } else if (offsetscomp_access == USE_MMAP_PRELOAD) {
      if (snps_root) {
	fprintf(stderr,"Pre-loading %s (%s) offsets, kmer %d, interval %d...",
		idx_filesuffix,snps_root,new->index1part,new->index1interval);
      } else {
	fprintf(stderr,"Pre-loading %s offsets, kmer %d, interval %d...",
		idx_filesuffix,new->index1part,new->index1interval);
      }
      new->offsetscomp = (Offsetscomp_T *) Access_mmap_and_preload(&new->offsetscomp_fd,&new->offsetscomp_len,&npages,&seconds,
								   filenames->offsets_filename,sizeof(Offsetscomp_T));
      if (new->offsetscomp == NULL) {
	fprintf(stderr,"insufficient memory (will use disk file instead, but program may not run)\n");
#ifdef PMAP
	new->offsetscomp_access = FILEIO;
#else
	exit(9);
#endif
      } else {
	comma = Genomicpos_commafmt(new->offsetscomp_len);
	fprintf(stderr,"done (%s bytes, %d pages, %.2f sec)\n",comma,npages,seconds);
	FREE(comma);
	new->offsetscomp_access = MMAPPED;
      }

    } else if (offsetscomp_access == USE_MMAP_ONLY) {
      new->offsetscomp = (Offsetscomp_T *) Access_mmap(&new->offsetscomp_fd,&new->offsetscomp_len,
						       filenames->offsets_filename,sizeof(Offsetscomp_T),/*randomp*/false);
      if (new->offsetscomp == NULL) {
	fprintf(stderr,"Insufficient memory for mmap of %s (will use disk file instead, but program may not run)\n",
		filenames->offsets_filename);
#ifdef PMAP
	new->offsetscomp_access = FILEIO;
#else
	exit(9);
#endif
      } else {
	new->offsetscomp_access = MMAPPED;
      }
#endif

    } else if (offsetscomp_access == USE_FILEIO) {
      fprintf(stderr,"Offsets file I/O access of %s not allowed\n",filenames->offsets_filename);
      exit(9);

    } else {
      fprintf(stderr,"Don't recognize offsetscomp_access type %d\n",offsetscomp_access);
      abort();
    }

 } else if ((filenames = Indexdb_get_filenames_bitpack64(
#ifdef PMAP
							 &(*alphabet),required_alphabet,
#endif
							 &new->offsetscomp_basesize,&new->index1part,&new->index1interval,
							 genomesubdir,fileroot,idx_filesuffix,snps_root,
							 required_basesize,required_index1part,required_interval,
							 /*offsets_only_p*/false)) != NULL) {
    /* Try bitpack compression  */
    *index1part = new->index1part;
    *index1interval = new->index1interval;

    if (expand_offsets_p == true) {
#ifdef LARGE_GENOMES
      fprintf(stderr,"Expansion of bitpack offsets not supported for large genomes\n");
      abort();
#else      
      fprintf(stderr,"Offsets compression type: none (bitpack expanded)\n");
      new->compression_type = NO_COMPRESSION;

      *basesize = *index1part;
#ifdef PMAP
      *alphabet_size = Alphabet_get_size(*alphabet);
      new->offsetscomp_blocksize = 1;
      basespace = power(*alphabet_size,new->index1part);
#else
      new->offsetscomp_blocksize = 1;
      basespace = power(4,new->index1part);
#endif
      new->gammaptrs = (Gammaptr_T *) CALLOC(basespace+1,sizeof(Gammaptr_T));
      for (base = 0; base <= basespace; base++) {
	new->gammaptrs[base] = base;
      }

#ifdef PMAP
      new->offsetscomp = Indexdb_offsets_from_bitpack(filenames->pointers_filename,filenames->offsets_filename,
						      *alphabet_size,new->index1part);
#else
      new->offsetscomp = Indexdb_offsets_from_bitpack(filenames->pointers_filename,filenames->offsets_filename,
						      new->offsetscomp_basesize,new->index1part);
#endif
      new->offsetscomp_access = ALLOCATED;

#endif

#ifdef PMAP
#else
      /* Sanity check on positions filesize */
      poly_T = ~(~0UL << 2*new->index1part);
      end0 = new->offsetscomp[poly_T+1];
#ifdef LARGE_GENOMES
      if ((filesize = Access_filesize(filenames->positions_high_filename)) != end0 * (off_t) sizeof(unsigned char)) {
	fprintf(stderr,"Something is wrong with the genomic index: expected file size for %s is %lu, but observed %lu.\n",
		filenames->positions_high_filename,end0*sizeof(unsigned char),filesize);
	abort();
      }
#endif
      if ((filesize = Access_filesize(filenames->positions_low_filename)) != end0 * (off_t) sizeof(UINT4)) {
	fprintf(stderr,"Something is wrong with the genomic index: expected file size for %s is %lu, but observed %lu.\n",
		filenames->positions_low_filename,end0*sizeof(UINT4),filesize);
	abort();
      }
#endif	/* PMAP */

    } else {
      fprintf(stderr,"Offsets compression type: bitpack\n");
      new->compression_type = BITPACK64_COMPRESSION;

      *basesize = new->offsetscomp_basesize;
#ifdef PMAP
      *alphabet_size = Alphabet_get_size(*alphabet);
      new->offsetscomp_blocksize = power(*alphabet_size,(*index1part) - new->offsetscomp_basesize);
#else
      new->offsetscomp_blocksize = power(4,(*index1part) - new->offsetscomp_basesize);
#endif

      if (new->index1part == new->offsetscomp_basesize) {
#ifdef PMAP
	basespace = power(*alphabet_size,new->offsetscomp_basesize);
#else
	basespace = power(4,new->offsetscomp_basesize);
#endif
	new->gammaptrs = (Gammaptr_T *) CALLOC(basespace+1,sizeof(Gammaptr_T));
	for (base = 0; base <= basespace; base++) {
	  new->gammaptrs[base] = base;
	}

      } else {
	/* bitpackptrs and gammaptrs always ALLOCATED */
	if (snps_root) {
	  fprintf(stderr,"Allocating memory for %s (%s) offset pointers, kmer %d, interval %d...",
		  idx_filesuffix,snps_root,new->index1part,new->index1interval);
	} else {
	  fprintf(stderr,"Allocating memory for %s offset pointers, kmer %d, interval %d...",
		  idx_filesuffix,new->index1part,new->index1interval);
	}
	new->gammaptrs = (Gammaptr_T *) Access_allocated(&new->gammaptrs_len,&seconds,
							 filenames->pointers_filename,sizeof(Gammaptr_T));

	comma = Genomicpos_commafmt(new->gammaptrs_len);
	fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma,seconds);
	FREE(comma);
      }

      /* offsetscomp could be ALLOCATED or MMAPPED +/- PRELOAD */
      if (offsetscomp_access == USE_ALLOCATE) {
	if (snps_root) {
	  fprintf(stderr,"Allocating memory for %s (%s) offsets, kmer %d, interval %d...",
		  idx_filesuffix,snps_root,new->index1part,new->index1interval);
	} else {
	  fprintf(stderr,"Allocating memory for %s offsets, kmer %d, interval %d...",
		  idx_filesuffix,new->index1part,new->index1interval);
	}
	new->offsetscomp = (Offsetscomp_T *) Access_allocated(&new->offsetscomp_len,&seconds,
							      filenames->offsets_filename,sizeof(Offsetscomp_T));
	if (new->offsetscomp == NULL) {
	  fprintf(stderr,"insufficient memory (need to use a lower batch mode (-B))\n");
	  exit(9);
	} else {
	  comma = Genomicpos_commafmt(new->offsetscomp_len);
	  fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma,seconds);
	  FREE(comma);
	  new->offsetscomp_access = ALLOCATED;
	}

#ifdef HAVE_MMAP
      } else if (offsetscomp_access == USE_MMAP_PRELOAD) {
	if (snps_root) {
	  fprintf(stderr,"Pre-loading %s (%s) offsets, kmer %d, interval %d...",
		  idx_filesuffix,snps_root,new->index1part,new->index1interval);
	} else {
	  fprintf(stderr,"Pre-loading %s offsets, kmer %d, interval %d...",
		  idx_filesuffix,new->index1part,new->index1interval);
	}
	new->offsetscomp = (Offsetscomp_T *) Access_mmap_and_preload(&new->offsetscomp_fd,&new->offsetscomp_len,&npages,&seconds,
								     filenames->offsets_filename,sizeof(Offsetscomp_T));
	if (new->offsetscomp == NULL) {
	  fprintf(stderr,"insufficient memory (will use disk file instead, but program may not run)\n");
#ifdef PMAP
	  new->offsetscomp_access = FILEIO;
#else
	  exit(9);
#endif
	} else {
	  comma = Genomicpos_commafmt(new->offsetscomp_len);
	  fprintf(stderr,"done (%s bytes, %d pages, %.2f sec)\n",comma,npages,seconds);
	  FREE(comma);
	  new->offsetscomp_access = MMAPPED;
	}

      } else if (offsetscomp_access == USE_MMAP_ONLY) {
	new->offsetscomp = (Offsetscomp_T *) Access_mmap(&new->offsetscomp_fd,&new->offsetscomp_len,
							 filenames->offsets_filename,sizeof(Offsetscomp_T),/*randomp*/false);
	if (new->offsetscomp == NULL) {
	  fprintf(stderr,"Insufficient memory for mmap of %s (will use disk file instead, but program may not run)\n",
		  filenames->offsets_filename);
#ifdef PMAP
	  new->offsetscomp_access = FILEIO;
#else
	  exit(9);
#endif
	} else {
	  new->offsetscomp_access = MMAPPED;
	}
#endif

      } else if (offsetscomp_access == USE_FILEIO) {
	fprintf(stderr,"Offsetscomp file I/O access of %s not allowed\n",filenames->offsets_filename);
	exit(9);

      } else {
	fprintf(stderr,"Don't recognize offsetscomp_access type %d\n",offsetscomp_access);
	abort();
      }

      Bitpack64_read_setup();
#ifdef PMAP
#else
      /* Sanity check on positions filesize */
#ifdef LARGE_GENOMES
      if (filenames->pages_filename != NULL) {
	new->offsetspages = (UINT4 *) Access_allocated(&bitpackpages_len,&seconds,filenames->pages_filename,sizeof(UINT4));
      } else {
	new->offsetspages = (UINT4 *) MALLOC(1*sizeof(UINT4));
	new->offsetspages[0] = -1U;
      }
#endif
      poly_T = ~(~0UL << 2*new->index1part);
#ifdef LARGE_GENOMES
      ptr0 = Bitpack64_offsetptr_huge(&end0,poly_T,new->offsetspages,new->gammaptrs,new->offsetscomp);
#else
      ptr0 = Bitpack64_offsetptr(&end0,poly_T,new->gammaptrs,new->offsetscomp);
#endif

#ifdef LARGE_GENOMES
      ptr0 = Bitpack64_offsetptr_huge(&end0,poly_T,new->offsetspages,new->gammaptrs,new->offsetscomp);
      if ((filesize = Access_filesize(filenames->positions_high_filename)) != end0 * (off_t) sizeof(unsigned char)) {
	fprintf(stderr,"Something is wrong with the genomic index: expected file size for %s is %lu, but observed %lu.\n",
		filenames->positions_high_filename,end0*sizeof(unsigned char),filesize);
	abort();
      }
#endif
      if ((filesize = Access_filesize(filenames->positions_low_filename)) != end0 * (off_t) sizeof(UINT4)) {
	fprintf(stderr,"Something is wrong with the genomic index: expected file size for %s is %lu, but observed %lu.\n",
		filenames->positions_low_filename,end0*sizeof(UINT4),filesize);
	abort();
      }
#endif	/* PMAP */

    }

  } else if ((filenames = Indexdb_get_filenames_gamma(
#ifdef PMAP
						      &(*alphabet),required_alphabet,
#endif
						      &new->offsetscomp_basesize,&new->index1part,&new->index1interval,
						      genomesubdir,fileroot,idx_filesuffix,snps_root,
						      required_basesize,required_index1part,required_interval,
						      /*offsets_only_p*/false)) != NULL) {
    /* Try gamma compression */
    *index1part = new->index1part;
    *index1interval = new->index1interval;

    if (expand_offsets_p == true) {
#ifdef LARGE_GENOMES
      fprintf(stderr,"Expansion of gamma offsets not supported for large genomes\n");
      abort();
#else
      fprintf(stderr,"Offsets compression type: none (gamma expanded)\n");
      new->compression_type = NO_COMPRESSION;

      *basesize = *index1part;
#ifdef PMAP
      *alphabet_size = Alphabet_get_size(*alphabet);
      new->offsetscomp_blocksize = 1;
      basespace = power(*alphabet_size,new->index1part);
#else
      new->offsetscomp_blocksize = 1;
      basespace = power(4,new->index1part);
#endif
      new->gammaptrs = (Gammaptr_T *) CALLOC(basespace+1,sizeof(Gammaptr_T));
      for (base = 0; base <= basespace; base++) {
	new->gammaptrs[base] = base;
      }

#ifdef PMAP
      new->offsetscomp = Indexdb_offsets_from_gammas(filenames->pointers_filename,filenames->offsets_filename,
						     new->offsetscomp_basesize,*alphabet_size,new->index1part);
#else
      new->offsetscomp = Indexdb_offsets_from_gammas(filenames->pointers_filename,filenames->offsets_filename,
						     new->offsetscomp_basesize,new->index1part);
#endif
      new->offsetscomp_access = ALLOCATED;

#endif

    } else {
      fprintf(stderr,"Offsets compression type: gamma\n");
      new->compression_type = GAMMA_COMPRESSION;

      *basesize = new->offsetscomp_basesize;
#ifdef PMAP
      *alphabet_size = Alphabet_get_size(*alphabet);
      new->offsetscomp_blocksize = power(*alphabet_size,(*index1part) - new->offsetscomp_basesize);
#else
      new->offsetscomp_blocksize = power(4,(*index1part) - new->offsetscomp_basesize);
#endif

      if (new->index1part == new->offsetscomp_basesize) {
#ifdef PMAP
	basespace = power(*alphabet_size,new->offsetscomp_basesize);
#else
	basespace = power(4,new->offsetscomp_basesize);
#endif
	new->gammaptrs = (Gammaptr_T *) CALLOC(basespace+1,sizeof(Gammaptr_T));
	for (base = 0; base <= basespace; base++) {
	  new->gammaptrs[base] = base;
	}

      } else {
	/* gammaptrs always ALLOCATED */
	if (snps_root) {
	  fprintf(stderr,"Allocating memory for %s (%s) gammaptrs, kmer %d, interval %d...",
		  idx_filesuffix,snps_root,new->index1part,new->index1interval);
	} else {
	  fprintf(stderr,"Allocating memory for %s gammaptrs, kmer %d, interval %d...",
		  idx_filesuffix,new->index1part,new->index1interval);
	}
	new->gammaptrs = (Gammaptr_T *) Access_allocated(&new->gammaptrs_len,&seconds,
							 filenames->pointers_filename,sizeof(Gammaptr_T));
	comma = Genomicpos_commafmt(new->gammaptrs_len);
	fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma,seconds);
	FREE(comma);
      }

      /* offsetscomp could be ALLOCATED or MMAPPED +/- PRELOAD */
      if (offsetscomp_access == USE_ALLOCATE) {
	if (snps_root) {
	  fprintf(stderr,"Allocating memory for %s (%s) offsets, kmer %d, interval %d...",
		  idx_filesuffix,snps_root,new->index1part,new->index1interval);
	} else {
	  fprintf(stderr,"Allocating memory for %s offsets, kmer %d, interval %d...",
		  idx_filesuffix,new->index1part,new->index1interval);
	}
	new->offsetscomp = (Offsetscomp_T *) Access_allocated(&new->offsetscomp_len,&seconds,
							      filenames->offsets_filename,sizeof(Offsetscomp_T));
	if (new->offsetscomp == NULL) {
	  fprintf(stderr,"insufficient memory (need to use a lower batch mode (-B))\n");
	  exit(9);
	} else {
	  comma = Genomicpos_commafmt(new->offsetscomp_len);
	  fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma,seconds);
	  FREE(comma);
	  new->offsetscomp_access = ALLOCATED;
	}

#ifdef HAVE_MMAP
      } else if (offsetscomp_access == USE_MMAP_PRELOAD) {
	if (snps_root) {
	  fprintf(stderr,"Pre-loading %s (%s) offsets, kmer %d, interval %d...",
		  idx_filesuffix,snps_root,new->index1part,new->index1interval);
	} else {
	  fprintf(stderr,"Pre-loading %s offsets, kmer %d, interval %d...",
		  idx_filesuffix,new->index1part,new->index1interval);
	}
	new->offsetscomp = (Offsetscomp_T *) Access_mmap_and_preload(&new->offsetscomp_fd,&new->offsetscomp_len,&npages,&seconds,
								     filenames->offsets_filename,sizeof(Offsetscomp_T));
	if (new->offsetscomp == NULL) {
	  fprintf(stderr,"insufficient memory (will use disk file instead, but program may not run)\n");
#ifdef PMAP
	  new->offsetscomp_access = FILEIO;
#else
	  exit(9);
#endif
	} else {
	  comma = Genomicpos_commafmt(new->offsetscomp_len);
	  fprintf(stderr,"done (%s bytes, %d pages, %.2f sec)\n",comma,npages,seconds);
	  FREE(comma);
	  new->offsetscomp_access = MMAPPED;
	}

      } else if (offsetscomp_access == USE_MMAP_ONLY) {
	new->offsetscomp = (Offsetscomp_T *) Access_mmap(&new->offsetscomp_fd,&new->offsetscomp_len,
							 filenames->offsets_filename,sizeof(Offsetscomp_T),/*randomp*/false);
	if (new->offsetscomp == NULL) {
	  fprintf(stderr,"Insufficient memory for mmap of %s (will use disk file instead, but program may not run)\n",
		  filenames->offsets_filename);
#ifdef PMAP
	  new->offsetscomp_access = FILEIO;
#else
	  exit(9);
#endif
	} else {
	  new->offsetscomp_access = MMAPPED;
	}
#endif

      } else if (offsetscomp_access == USE_FILEIO) {
	fprintf(stderr,"Offsetscomp file I/O access of %s not allowed\n",filenames->offsets_filename);
	exit(9);

      } else {
	fprintf(stderr,"Don't recognize offsetscomp_access type %d\n",offsetscomp_access);
	abort();
      }
    }

#ifdef PMAP
#else
    /* Sanity check on positions filesize */
    poly_T = ~(~0UL << 2*new->index1part);
#ifdef WORDS_BIGENDIAN
    if (offsetscomp_access == ALLOCATED) {
      ptr0 = Genome_offsetptr_from_gammas(&end0,new->gammaptrs,new->offsetscomp,new->offsetscomp_blocksize,poly_T);
    } else {
      ptr0 = Genome_offsetptr_from_gammas_bigendian(&end0,new->gammaptrs,new->offsetscomp,new->offsetscomp_blocksize,poly_T);
    }
#else
    ptr0 = Genome_offsetptr_from_gammas(&end0,new->gammaptrs,new->offsetscomp,new->offsetscomp_blocksize,poly_T);
#endif

    if ((filesize = Access_filesize(filenames->positions_low_filename)) != end0 * (off_t) sizeof(UINT4)) {
      fprintf(stderr,"Something is wrong with the genomic index: expected file size for %s is %lu, but observed %lu.\n",
	      filenames->positions_low_filename,end0*sizeof(UINT4),filesize);
      abort();
    }
#endif	/* PMAP */

  } else {
    fprintf(stderr,"Cannot find genomic index files in either current or old format\n");
    exit(9);
  }


  /* Read or memory map positions file */
  if (positions_access == USE_ALLOCATE) {
    if (snps_root) {
      fprintf(stderr,"Allocating memory for %s (%s) positions, kmer %d, interval %d...",
	      idx_filesuffix,snps_root,new->index1part,new->index1interval);
    } else {
      fprintf(stderr,"Allocating memory for %s positions, kmer %d, interval %d...",
	      idx_filesuffix,new->index1part,new->index1interval);
    }
#ifdef LARGE_GENOMES
    new->positions_high = (unsigned char *) Access_allocated(&new->positions_high_len,&seconds,
							     filenames->positions_high_filename,sizeof(unsigned char));
    if (new->positions_high == NULL) {
      fprintf(stderr,"insufficient memory (need to use a lower batch mode (-B)\n");
      exit(9);
    } else {
      comma = Genomicpos_commafmt(new->positions_high_len);
      fprintf(stderr,"done (%s bytes, %.2f sec), ",comma,seconds);
      FREE(comma);

      new->positions_low = (UINT4 *) Access_allocated(&new->positions_low_len,&seconds,
						      filenames->positions_low_filename,sizeof(UINT4));
      if (new->positions_low == NULL) {
	fprintf(stderr,"insufficient memory (need to use a lower batch mode (-B)\n");
	exit(9);
      } else {
	comma = Genomicpos_commafmt(new->positions_low_len);
	fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma,seconds);
	FREE(comma);

	new->positions_access = ALLOCATED;
      }
    }
#else
    new->positions = (UINT4 *) Access_allocated(&new->positions_len,&seconds,
						filenames->positions_low_filename,sizeof(UINT4));
    if (new->positions == NULL) {
      fprintf(stderr,"insufficient memory (need to use a lower batch mode (-B)\n");
      exit(9);
    } else {
      comma = Genomicpos_commafmt(new->positions_len);
      fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma,seconds);
      FREE(comma);
      new->positions_access = ALLOCATED;
    }
#endif
    

#ifdef HAVE_MMAP
  } else if (positions_access == USE_MMAP_PRELOAD) {
    if (snps_root) {
      fprintf(stderr,"Pre-loading %s (%s) positions, kmer %d, interval %d...",
	      idx_filesuffix,snps_root,new->index1part,new->index1interval);
    } else {
      fprintf(stderr,"Pre-loading %s positions, kmer %d, interval %d...",
	      idx_filesuffix,new->index1part,new->index1interval);
    }
#ifdef LARGE_GENOMES
    new->positions_high = (unsigned char *) Access_mmap_and_preload(&new->positions_high_fd,&new->positions_high_len,&npages,&seconds,
								  filenames->positions_high_filename,sizeof(unsigned char));
    if (new->positions_high == NULL) {
      fprintf(stderr,"insufficient memory (will use disk file instead, but program will be slow)\n");
      new->positions_access = FILEIO;
    } else {
      comma = Genomicpos_commafmt(new->positions_high_len);
      fprintf(stderr,"done (%s bytes, %d pages, %.2f sec), ",comma,npages,seconds);
      FREE(comma);

      new->positions_low = (UINT4 *) Access_mmap_and_preload(&new->positions_low_fd,&new->positions_low_len,&npages,&seconds,
							     filenames->positions_low_filename,sizeof(UINT4));
      if (new->positions_low == NULL) {
	fprintf(stderr,"insufficient memory (will use disk file instead, but program will be slow)\n");
	new->positions_access = FILEIO;
      } else {
	comma = Genomicpos_commafmt(new->positions_low_len);
	fprintf(stderr,"done (%s bytes, %d pages, %.2f sec)\n",comma,npages,seconds);
	FREE(comma);

	new->positions_access = MMAPPED;
      }
    }
#else
    new->positions = (UINT4 *) Access_mmap_and_preload(&new->positions_fd,&new->positions_len,&npages,&seconds,
						       filenames->positions_low_filename,sizeof(UINT4));
    if (new->positions == NULL) {
      fprintf(stderr,"insufficient memory (will use disk file instead, but program will be slow)\n");
      new->positions_access = FILEIO;
    } else {
      comma = Genomicpos_commafmt(new->positions_len);
      fprintf(stderr,"done (%s bytes, %d pages, %.2f sec)\n",comma,npages,seconds);
      FREE(comma);
      new->positions_access = MMAPPED;
    }
#endif


  } else if (positions_access == USE_MMAP_ONLY) {
#ifdef LARGE_GENOMES
    new->positions_high = (unsigned char *) Access_mmap(&new->positions_high_fd,&new->positions_high_len,
							filenames->positions_high_filename,sizeof(unsigned char),/*randomp*/true);
    if (new->positions_high == NULL) {
      fprintf(stderr,"Insufficient memory for mmap of %s (will use disk file instead, but program will be slow)\n",
	      filenames->positions_high_filename);
      new->positions_access = FILEIO;
    } else {
      new->positions_low = (UINT4 *) Access_mmap(&new->positions_low_fd,&new->positions_low_len,
						 filenames->positions_low_filename,sizeof(UINT4),/*randomp*/true);
      if (new->positions_low == NULL) {
	fprintf(stderr,"Insufficient memory for mmap of %s (will use disk file instead, but program will be slow)\n",
		filenames->positions_low_filename);
	new->positions_access = FILEIO;
      } else {
	new->positions_access = MMAPPED;
      }
    }
#else
    new->positions = (UINT4 *) Access_mmap(&new->positions_fd,&new->positions_len,
					   filenames->positions_low_filename,sizeof(UINT4),/*randomp*/true);

    if (new->positions == NULL) {
      fprintf(stderr,"Insufficient memory for mmap of %s (will use disk file instead, but program will be slow)\n",
	      filenames->positions_low_filename);
      new->positions_access = FILEIO;
    } else {
      new->positions_access = MMAPPED;
    }
#endif

#endif

  } else if (positions_access == USE_FILEIO) {
    new->positions_access = FILEIO;
  } else {
    fprintf(stderr,"Don't recognize positions_access %d\n",positions_access);
    abort();
  }

#ifdef HAVE_PTHREAD
  if (new->positions_access == FILEIO) {
    pthread_mutex_init(&new->positions_read_mutex,NULL);
  }
#endif

  Filenames_free(&filenames);

  return new;
}


/************************************************************************
 *   Debugging procedures
 ************************************************************************/

#ifndef PMAP

/*                87654321 */
#define RIGHT_A 0x00000000
#define RIGHT_C 0x00000001
#define RIGHT_G 0x00000002
#define RIGHT_T 0x00000003

/*                      87654321 */
#define LOW_TWO_BITS  0x00000003

#if (defined(DEBUG0) || defined(DEBUG1) || defined(DEBUG2))
static char *
shortoligo_nt (Storedoligomer_T oligo, Width_T oligosize) {
  char *nt;
  Width_T i, j;
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
#endif

#endif


/************************************************************************
 *   Read procedures
 ************************************************************************/


static void
positions_move_absolute_1 (int positions_fd, off_t ptr) {
  off_t offset = ptr*((off_t) sizeof(unsigned char));

  if (lseek(positions_fd,offset,SEEK_SET) < 0) {
    fprintf(stderr,"Attempted to do lseek on offset %ld*%lu=%ld\n",
	    ptr,sizeof(unsigned char),offset);
    perror("Error in indexdb.c, positions_move_absolute");
    exit(9);
  }
  return;
}

static void
positions_move_absolute_4 (int positions_fd, off_t ptr) {
  off_t offset = ptr*((off_t) sizeof(UINT4));

  if (lseek(positions_fd,offset,SEEK_SET) < 0) {
    fprintf(stderr,"Attempted to do lseek on offset %ld*%lu=%ld\n",
	    ptr,sizeof(UINT4),offset);
    perror("Error in indexdb.c, positions_move_absolute");
    exit(9);
  }
  return;
}

#if 0
static Univcoord_T
positions_read_forward (int positions_fd) {
  Univcoord_T value;
  char buffer[4];

  read(positions_fd,buffer,4);

  value = (buffer[3] & 0xff);
  value <<= 8;
  value |= (buffer[2] & 0xff);
  value <<= 8;
  value |= (buffer[1] & 0xff);
  value <<= 8;
  value |= (buffer[0] & 0xff);

  return value;
}
#endif

#ifdef LARGE_GENOMES
static void
positions_read_multiple_large (int positions_high_fd, int positions_low_fd, Univcoord_T *values, int n) {
  int i;
  Univcoord_T value;
  unsigned char buffer[4];

#ifdef WORDS_BIGENDIAN
  /* Need to keep in bigendian format */
  for (i = 0; i < n; i++) {
    read(positions_high_fd,buffer,1);
    value = (buffer[0] & 0xff);
    value <<= 8;

    read(positions_low_fd,buffer,4);
    value |= (buffer[0] & 0xff);
    value <<= 8;
    value |= (buffer[1] & 0xff);
    value <<= 8;
    value |= (buffer[2] & 0xff);
    value <<= 8;
    value |= (buffer[3] & 0xff);

    values[i] = value;
  }
#else
  for (i = 0; i < n; i++) {
    read(positions_high_fd,buffer,1);
    value = (buffer[0] & 0xff);
    value <<= 8;

    read(positions_low_fd,buffer,4);
    value |= (buffer[3] & 0xff);
    value <<= 8;
    value |= (buffer[2] & 0xff);
    value <<= 8;
    value |= (buffer[1] & 0xff);
    value <<= 8;
    value |= (buffer[0] & 0xff);

    values[i] = value;
  }
#endif

  return;
}

static void
positions_copy_multiple_large (Univcoord_T *positions, unsigned char *positions_high, UINT4 *positions_low, int n) {
  int i;

  for (i = 0; i < n; i++) {
    positions[i] = ((Univcoord_T) positions_high[i] << 32) + positions_low[i];
  }

  return;
}

#else

static void
positions_read_multiple (int positions_fd, Univcoord_T *values, int n) {
  int i;
  Univcoord_T value;
  unsigned char buffer[4];

#ifdef WORDS_BIGENDIAN
  /* Need to keep in bigendian format */
  for (i = 0; i < n; i++) {
    read(positions_fd,buffer,4);

    value = (buffer[0] & 0xff);
    value <<= 8;
    value |= (buffer[1] & 0xff);
    value <<= 8;
    value |= (buffer[2] & 0xff);
    value <<= 8;
    value |= (buffer[3] & 0xff);

    values[i] = value;
  }
#else
  for (i = 0; i < n; i++) {
    read(positions_fd,buffer,4);

    value = (buffer[3] & 0xff);
    value <<= 8;
    value |= (buffer[2] & 0xff);
    value <<= 8;
    value |= (buffer[1] & 0xff);
    value <<= 8;
    value |= (buffer[0] & 0xff);

    values[i] = value;
  }
#endif

  return;
}
#endif




#if 0
static Univcoord_T
positions_read_backward (int positions_fd) {
  Univcoord_T value;
  char buffer[4];
  off_t reloffset = -2*((off_t) sizeof(Univcoord_T)); /* 1 to undo the effect of read */

  read(positions_fd,buffer,4);

  value = (buffer[3] & 0xff);
  value <<= 8;
  value |= (buffer[2] & 0xff);
  value <<= 8;
  value |= (buffer[1] & 0xff);
  value <<= 8;
  value |= (buffer[0] & 0xff);

  if (lseek(positions_fd,reloffset,SEEK_CUR) < 0) {
    fprintf(stderr,"Attempted to do lseek on relative offset %ld\n",(long int) reloffset);
    perror("Error in indexdb.c, positions_read_backward");
    exit(9);
  }
    
  return value;
}
#endif


Positionsptr_T *
Indexdb_offsets_from_bitpack (char *bitpackptrsfile, char *offsetscompfile, 
#ifdef PMAP
			      int alphabet_size, Width_T index1part_aa
#else
			      Width_T offsetscomp_basesize , Width_T index1part
#endif
			      ) {
  UINT4 *bitpackptrs;
  Offsetscomp_T *offsetscomp;
  int bitpackptrs_fd, offsetscomp_fd;
  size_t bitpackptrs_len, offsetscomp_len;
  UINT4 *offsets = NULL;
  Oligospace_T oligospace, oligoi;
#ifndef PMAP
  Blocksize_T blocksize;
#endif
  double seconds;

#ifdef PMAP
  oligospace = power(alphabet_size,index1part_aa);
#else
  oligospace = power(4,index1part);
  blocksize = power(4,index1part - offsetscomp_basesize);
#endif


#ifdef HAVE_MMAP
  bitpackptrs = (UINT4 *) Access_mmap(&bitpackptrs_fd,&bitpackptrs_len,bitpackptrsfile,sizeof(UINT4),/*randomp*/false);
  offsetscomp = (Offsetscomp_T *) Access_mmap(&offsetscomp_fd,&offsetscomp_len,offsetscompfile,sizeof(Offsetscomp_T),/*randomp*/false);
#else
  bitpackptrs = (UINT4 *) Access_allocated(&bitpackptrs_len,&seconds,bitpackptrsfile,sizeof(UINT4));
  offsetscomp = (Offsetscomp_T *) Access_allocated(&offsetscomp_len,&seconds,offsetscompfile,sizeof(Offsetscomp_T));
#endif

#ifdef OLIGOSPACE_NOT_LONG
  fprintf(stderr,"Allocating memory (%u 4-byte words) for offsets, kmer %d...",oligospace+1U,
#ifdef PMAP
	  index1part_aa
#else
	  index1part
#endif
	  );
#else
  fprintf(stderr,"Allocating memory (%lu 4-byte words) for offsets, kmer %d...",oligospace+1UL,
#ifdef PMAP
	  index1part_aa
#else
	  index1part
#endif
	  );
#endif

  /* Bitpack procedures start from offsets[1], so we need to print offsets[0] as a special case */
  offsets = (Positionsptr_T *) MALLOC((oligospace+1) * sizeof(Positionsptr_T));

  if (offsets == NULL) {
    fprintf(stderr,"cannot allocated requested memory.  Cannot run expand offsets on this machine.\n");
    exit(9);
  } else {
    fprintf(stderr,"done\n");
  }

  fprintf(stderr,"Expanding bitpackcomp into offsets...");
  Bitpack64_read_setup();
#ifdef PMAP
  for (oligoi = 0UL; oligoi <= oligospace; oligoi += 1) {
    offsets[oligoi] = Bitpack64_offsetptr_only(oligoi,bitpackptrs,offsetscomp);
  }
#else
  for (oligoi = 0UL; oligoi < oligospace; oligoi += blocksize) {
    Bitpack64_block_offsets(&(offsets[oligoi]),oligoi,bitpackptrs,offsetscomp);
  }
#endif

  fprintf(stderr,"done\n");

#ifdef HAVE_MMAP
  munmap((void *) offsetscomp,offsetscomp_len);
  munmap((void *) bitpackptrs,bitpackptrs_len);
#else
  FREE(offsetscomp);
  FREE(bitpackptrs);
#endif

  return offsets;
}


#if defined(HAVE_64_BIT) && defined(UTILITYP)
/* Used by utility programs */
Hugepositionsptr_T *
Indexdb_offsets_from_bitpack_huge (char *bitpackpagesfile, char *bitpackptrsfile, char *offsetscompfile, Width_T offsetscomp_basesize
#ifdef PMAP
				   , int alphabet_size, Width_T index1part_aa
#else
				   , Width_T index1part
#endif
				   ) {
  UINT4 *bitpackpages;
  UINT4 *bitpackptrs;
  Offsetscomp_T *offsetscomp;
  int bitpackptrs_fd, offsetscomp_fd;
  size_t bitpackpages_len, bitpackptrs_len, offsetscomp_len;
  Hugepositionsptr_T *offsets = NULL;
  Oligospace_T oligospace, oligoi;
  Blocksize_T blocksize;
  double seconds;

#ifdef PMAP
  oligospace = power(alphabet_size,index1part_aa);
  blocksize = power(alphabet_size,index1part_aa - offsetscomp_basesize);
#else
  oligospace = power(4,index1part);
  blocksize = power(4,index1part - offsetscomp_basesize);
#endif

  if (blocksize == 1) {
    return (Hugepositionsptr_T *) Access_allocated(&offsetscomp_len,&seconds,offsetscompfile,sizeof(Hugepositionsptr_T));

  } else {

    if (bitpackpagesfile == NULL) {
      bitpackpages = (UINT4 *) MALLOC(1*sizeof(UINT4));
      bitpackpages[0] = -1U;
    } else {
      bitpackpages = (UINT4 *) Access_allocated(&bitpackpages_len,&seconds,bitpackpagesfile,sizeof(UINT4));
    }
#ifdef HAVE_MMAP
    bitpackptrs = (UINT4 *) Access_mmap(&bitpackptrs_fd,&bitpackptrs_len,bitpackptrsfile,sizeof(UINT4),/*randomp*/false);
    offsetscomp = (Offsetscomp_T *) Access_mmap(&offsetscomp_fd,&offsetscomp_len,offsetscompfile,sizeof(Offsetscomp_T),/*randomp*/false);
#else
    bitpackptrs = (UINT4 *) Access_allocated(&bitpackptrs_len,&seconds,bitpackptrsfile,sizeof(UINT4));
    offsetscomp = (Offsetscomp_T *) Access_allocated(&offsetscomp_len,&seconds,offsetscompfile,sizeof(Offsetscomp_T));
#endif

#ifdef OLIGOSPACE_NOT_LONG
    fprintf(stderr,"Allocating memory (%u 8-byte words) for offsets, kmer %d...",oligospace+1U,
#ifdef PMAP
	    index1part_aa
#else
	    index1part
#endif
	    );
#else
    fprintf(stderr,"Allocating memory (%lu 8-byte words) for offsets, kmer %d...",oligospace+1UL,
#ifdef PMAP
	    index1part_aa
#else
	    index1part
#endif
	    );
#endif

    /* Bitpack procedures start from offsets[1], so we need to print offsets[0] as a special case */
    offsets = (Hugepositionsptr_T *) MALLOC((oligospace+1) * sizeof(Hugepositionsptr_T));

    if (offsets == NULL) {
      fprintf(stderr,"cannot allocated requested memory.  Cannot run expand offsets on this machine.\n");
      exit(9);
    } else {
      fprintf(stderr,"done\n");
    }

    fprintf(stderr,"Expanding bitpackcomp into offsets...");
    Bitpack64_read_setup();
#ifdef PMAP
    for (oligoi = 0UL; oligoi <= oligospace; oligoi += 1) {
      offsets[oligoi] = Bitpack64_offsetptr_only_huge(oligoi,bitpackpages,bitpackptrs,offsetscomp);
    }
#else
    for (oligoi = 0UL; oligoi < oligospace; oligoi += blocksize) {
      Bitpack64_block_offsets_huge(&(offsets[oligoi]),oligoi,bitpackpages,bitpackptrs,offsetscomp);
    }
#endif

    fprintf(stderr,"done\n");

#ifdef HAVE_MMAP
    munmap((void *) offsetscomp,offsetscomp_len);
    munmap((void *) bitpackptrs,bitpackptrs_len);
#else
    FREE(offsetscomp);
    FREE(bitpackptrs);
#endif
    FREE(bitpackpages);

    return offsets;
  }
}
#endif



Positionsptr_T *
Indexdb_offsets_from_gammas (char *gammaptrsfile, char *offsetscompfile, Width_T offsetscomp_basesize
#ifdef PMAP
			     , int alphabet_size, Width_T index1part_aa
#else
			     , Width_T index1part
#endif
			     ) {
  Gammaptr_T *gammaptrs;
  Offsetscomp_T *offsetscomp, *ptr;
  int gammaptrs_fd, offsetscomp_fd;
  size_t gammaptrs_len, offsetscomp_len;
  Positionsptr_T *offsets = NULL;
  Oligospace_T oligospace, oligoi, oligok;
  Blocksize_T blocksize, j;
  double seconds;

  Positionsptr_T cum;
  int ctr;


#ifdef PMAP
  oligospace = power(alphabet_size,index1part_aa);
  blocksize = power(alphabet_size,index1part_aa - offsetscomp_basesize);
#else
  oligospace = power(4,index1part);
  blocksize = power(4,index1part - offsetscomp_basesize);
#endif

  if (blocksize == 1) {
    return (Positionsptr_T *) Access_allocated(&offsetscomp_len,&seconds,offsetscompfile,sizeof(Positionsptr_T));

  } else {

#ifdef HAVE_MMAP
    gammaptrs = (Gammaptr_T *) Access_mmap(&gammaptrs_fd,&gammaptrs_len,gammaptrsfile,sizeof(Gammaptr_T),/*randomp*/false);
    offsetscomp = (Offsetscomp_T *) Access_mmap(&offsetscomp_fd,&offsetscomp_len,offsetscompfile,sizeof(Offsetscomp_T),/*randomp*/false);
#else
    gammaptrs = (Gammaptr_T *) Access_allocated(&gammaptrs_len,&seconds,gammaptrsfile,sizeof(Gammaptr_T));
    offsetscomp = (Offsetscomp_T *) Access_allocated(&offsetscomp_len,&seconds,offsetscompfile,sizeof(Offsetscomp_T));
#endif

#ifdef OLIGOSPACE_NOT_LONG
    fprintf(stderr,"Allocating memory (%u words) for offsets, kmer %d...",oligospace+1U,
#ifdef PMAP
	    index1part_aa
#else
	    index1part
#endif
	    );
#else
    fprintf(stderr,"Allocating memory (%lu words) for offsets, kmer %d...",oligospace+1UL,
#ifdef PMAP
	    index1part_aa
#else
	    index1part
#endif
	    );
#endif

    offsets = (Positionsptr_T *) CALLOC(oligospace+1,sizeof(Positionsptr_T));
    if (offsets == NULL) {
      fprintf(stderr,"cannot allocated requested memory.  Cannot run expand offsets on this machine.\n");
      exit(9);
    } else {
      fprintf(stderr,"done\n");
    }


    fprintf(stderr,"Expanding offsetscomp into offsets...");

    ptr = offsetscomp;
    oligok = 0UL;


    for (oligoi = 0UL; oligoi < oligospace; oligoi += blocksize) {
#ifdef HAVE_MMAP
#ifdef WORDS_BIGENDIAN
      cum = offsets[oligok++] = Bigendian_convert_uint(*ptr++);
#else
      cum = offsets[oligok++] = *ptr++;
#endif
#else
      cum = offsets[oligok++] = *ptr++;
#endif

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
	offsets[oligok++] = cum;
      }
      if (ctr > 0) {
	ptr++;			/* Done with last gamma byte */
      }
    }

#ifdef HAVE_MMAP
#ifdef WORDS_BIGENDIAN
    offsets[oligok++] = Bigendian_convert_uint(*ptr++);
#else
    offsets[oligok++] = *ptr++;
#endif
#else
    offsets[oligok++] = *ptr++;
#endif

    fprintf(stderr,"done\n");

#ifdef HAVE_MMAP
    munmap((void *) offsetscomp,offsetscomp_len);
    munmap((void *) gammaptrs,gammaptrs_len);
#else
    FREE(offsetscomp);
    FREE(gammaptrs);
#endif

    return offsets;
  }
}


#ifdef PMAP

/* PMAP version.  Doesn't mask bottom 12 nt. */
Univcoord_T *
Indexdb_read (int *nentries, T this, Storedoligomer_T aaindex) {
  Univcoord_T *positions, bigendian, littleendian;
  Positionsptr_T ptr, ptr0, end0;
  int i;
  char byte1, byte2, byte3;

  debug0(printf("%u (%s)\n",aaindex,Alphabet_aaindex_aa(aaindex,this->alphabet)));

  if (this->compression_type == NO_COMPRESSION) {
#ifdef WORDS_BIGENDIAN
    if (this->offsetscomp_access == ALLOCATED) {
      ptr0 = this->offsetscomp[aaindex];
      end0 = this->offsetscomp[aaindex+1];
    } else {
      ptr0 = Bigendian_convert_uint(this->offsetscomp[aaindex]);
      end0 = Bigendian_convert_uint(this->offsetscomp[aaindex+1]);
    }
#else
    ptr0 = this->offsetscomp[aaindex];
    end0 = this->offsetscomp[aaindex+1];
#endif

  } else if (this->compression_type == BITPACK64_COMPRESSION) {
    ptr0 = Bitpack64_offsetptr(&end0,aaindex,this->gammaptrs,this->offsetscomp);

  } else if (this->compression_type == GAMMA_COMPRESSION) {
#ifdef WORDS_BIGENDIAN
    if (this->offsetscomp_access == ALLOCATED) {
      ptr0 = Genome_offsetptr_from_gammas(&end0,this->gammaptrs,this->offsetscomp,this->offsetscomp_blocksize,aaindex);
    } else {
      ptr0 = Genome_offsetptr_from_gammas_bigendian(&end0,this->gammaptrs,this->offsetscomp,this->offsetscomp_blocksize,aaindex);
    }
#else
    ptr0 = Genome_offsetptr_from_gammas(&end0,this->gammaptrs,this->offsetscomp,this->offsetscomp_blocksize,aaindex);
#endif

  } else {
    abort();
  }


  debug0(printf("offset pointers are %u and %u\n",ptr0,end0));

#ifdef ALLOW_DUPLICATES
  /* Not used */
  /* Skip backward over bad values, due to duplicates */
  if (this->positions_access == FILEIO) {
#ifdef HAVE_PTHREAD
    pthread_mutex_lock(&this->positions_read_mutex);
#endif
    positions_move_absolute(this->positions_fd,end0-1);
    while (end0 > ptr0 && positions_read_backward(this->positions_fd) == BADVAL) {
      end0--;
    }
#ifdef HAVE_PTHREAD
    pthread_mutex_unlock(&this->positions_read_mutex);
#endif
  } else {
    while (end0 > ptr0 && this->positions[end0-1] == BADVAL) {
      end0--;
    }
  }
#endif	/* ALLOW_DUPLICATES */

  if ((*nentries = end0 - ptr0) == 0) {
    return NULL;
  } else {
    positions = (Univcoord_T *) CALLOC(*nentries,sizeof(Univcoord_T));
    if (this->positions_access == FILEIO) {
#ifdef HAVE_PTHREAD
      pthread_mutex_lock(&this->positions_read_mutex);
#endif
#ifdef LARGE_GENOMES
      positions_move_absolute_1(this->positions_high_fd,ptr0);
      positions_move_absolute_4(this->positions_low_fd,ptr0);
      positions_read_multiple_large(this->positions_high_fd,this->positions_low_fd,positions,*nentries);
#else
      positions_move_absolute_4(this->positions_fd,ptr0);
      positions_read_multiple(this->positions_fd,positions,*nentries);
#endif
#ifdef HAVE_PTHREAD
      pthread_mutex_unlock(&this->positions_read_mutex);
#endif

    } else if (this->positions_access == ALLOCATED) {
#ifdef LARGE_GENOMES
      positions_copy_multiple_large(positions,&(this->positions_high[ptr0]),&(this->positions_low[ptr0]),*nentries);
#else
      memcpy(positions,&(this->positions[ptr0]),(*nentries)*sizeof(Univcoord_T));
#endif

    } else {
#ifdef WORDS_BIGENDIAN
#ifdef LARGE_GENOMES
      for (ptr = ptr0, i = 0; ptr < end0; ptr++, i++) {
	bigendian = (Univcoord_T) this->positions_high[ptr];
	bigendian <<= 8;

	littleendian = this->positions_low[ptr];
	bigendian |= littleendian & 0xff; /* 0 */
	bigendian <<= 8;
	bigendian |= ((littleendian >>= 8) & 0xff); /* 1 */
	bigendian <<= 8;
	bigendian |= ((littleendian >>= 8) & 0xff); /* 2 */
	bigendian <<= 8;
	bigendian |= ((littleendian >>= 8) & 0xff); /* 3 */
	positions[i] = bigendian;
      }
#else
      for (ptr = ptr0, i = 0; ptr < end0; ptr++, i++) {
	littleendian = this->positions[ptr];
	bigendian = littleendian & 0xff; /* 0 */
	bigendian <<= 8;
	bigendian |= ((littleendian >>= 8) & 0xff); /* 1 */
	bigendian <<= 8;
	bigendian |= ((littleendian >>= 8) & 0xff); /* 2 */
	bigendian <<= 8;
	bigendian |= ((littleendian >>= 8) & 0xff); /* 3 */
	positions[i] = bigendian;
      }
#endif

#else  /* littleendian */

#ifdef LARGE_GENOMES
      positions_copy_multiple_large(positions,&(this->positions_high[ptr0]),&(this->positions_low[ptr0]),*nentries);
#else
      memcpy(positions,&(this->positions[ptr0]),(*nentries)*sizeof(Univcoord_T));
#endif
#endif
    }
    debug0(
	   printf("%d entries:",*nentries);
	   for (i = 0; i < *nentries; i++) {
	     printf(" %u",positions[i]);
	   }
	   printf("\n");
	   );

    return positions;
  }
}

#else

/* GMAP version: Allocates memory */
Univcoord_T *
Indexdb_read (int *nentries, T this, Storedoligomer_T oligo) {
  Univcoord_T *positions;
  Positionsptr_T ptr0, end0;
  Storedoligomer_T part0;
#ifdef WORDS_BIGENDIAN
  int i;
  Positionsptr_T ptr;
  Univcoord_T bigendian, littleendian;
  char byte1, byte2, byte3;
#endif
#ifdef DEBUG0
  int j;
#endif


#if 0
  debug0(printf("%06X (%s)\n",oligo,shortoligo_nt(oligo,index1part)));
#endif
  part0 = oligo & poly_T;

  /* Ignore poly A and poly T on stage 1 */
  /* Was commented out */
  if (part0 == poly_A || part0 == poly_T) {
    *nentries = 0;
    return NULL;
  }

  if (this->compression_type == NO_COMPRESSION) {
#ifdef WORDS_BIGENDIAN
    if (this->offsetscomp_access == ALLOCATED) {
      ptr0 = this->offsetscomp[part0];
      end0 = this->offsetscomp[part0+1];
    } else {
      ptr0 = Bigendian_convert_uint(this->offsetscomp[part0]);
      end0 = Bigendian_convert_uint(this->offsetscomp[part0+1]);
    }
#else
    ptr0 = this->offsetscomp[part0];
    end0 = this->offsetscomp[part0+1];
#endif

  } else if (this->compression_type == BITPACK64_COMPRESSION) {
#ifdef LARGE_GENOMES
    ptr0 = Bitpack64_offsetptr_huge(&end0,part0,this->offsetspages,this->gammaptrs,this->offsetscomp);
#else
    ptr0 = Bitpack64_offsetptr(&end0,part0,this->gammaptrs,this->offsetscomp);
#endif

  } else if (this->compression_type == GAMMA_COMPRESSION) {
#ifdef WORDS_BIGENDIAN
    if (this->offsetscomp_access == ALLOCATED) {
      ptr0 = Genome_offsetptr_from_gammas(&end0,this->gammaptrs,this->offsetscomp,this->offsetscomp_blocksize,part0);
    } else {
      ptr0 = Genome_offsetptr_from_gammas_bigendian(&end0,this->gammaptrs,this->offsetscomp,this->offsetscomp_blocksize,part0);
    }
#else
    ptr0 = Genome_offsetptr_from_gammas(&end0,this->gammaptrs,this->offsetscomp,this->offsetscomp_blocksize,part0);
#endif

  } else {
    abort();
  }


#ifdef ALLOW_DUPLICATES
  /* Not used */
  /* Skip backward over bad values, due to duplicates */
  if (this->positions_access == FILEIO) {
#ifdef HAVE_PTHREAD
    pthread_mutex_lock(&this->positions_read_mutex);
#endif
    positions_move_absolute(this->positions_fd,end0-1);
    while (end0 > ptr0 && positions_read_backward(this->positions_fd) == BADVAL) {
      end0--;
    }
#ifdef HAVE_PTHREAD
    pthread_mutex_unlock(&this->positions_read_mutex);
#endif
  } else {
    while (end0 > ptr0 && this->positions[end0-1] == BADVAL) {
      end0--;
    }
  }
#endif

  if ((*nentries = end0 - ptr0) == 0) {
    return NULL;
  } else {
    debug0(printf("Indexdb_read: offset pointers are %u and %u\n",ptr0,end0));
  
    positions = (Univcoord_T *) CALLOC(*nentries,sizeof(Univcoord_T));
    if (this->positions_access == FILEIO) {
#ifdef HAVE_PTHREAD
      pthread_mutex_lock(&this->positions_read_mutex);
#endif
#ifdef LARGE_GENOMES
      positions_move_absolute_1(this->positions_high_fd,ptr0);
      positions_move_absolute_4(this->positions_low_fd,ptr0);
      positions_read_multiple_large(this->positions_high_fd,this->positions_low_fd,positions,*nentries);
#else
      positions_move_absolute_4(this->positions_fd,ptr0);
      positions_read_multiple(this->positions_fd,positions,*nentries);
#endif
#ifdef HAVE_PTHREAD
      pthread_mutex_unlock(&this->positions_read_mutex);
#endif
    } else if (this->positions_access == ALLOCATED) {
#ifdef LARGE_GENOMES
      positions_copy_multiple_large(positions,&(this->positions_high[ptr0]),&(this->positions_low[ptr0]),*nentries);
#else
      memcpy(positions,&(this->positions[ptr0]),(*nentries)*sizeof(Univcoord_T));
#endif

    } else {
#ifdef WORDS_BIGENDIAN
#ifdef LARGE_GENOMES
      for (ptr = ptr0, i = 0; ptr < end0; ptr++, i++) {
	bigendian = (Univcoord_T) this->positions_high[ptr];
	bigendian <<= 8;

	littleendian = this->positions_low[ptr];
	bigendian |= littleendian & 0xff; /* 0 */
	bigendian <<= 8;
	bigendian |= ((littleendian >>= 8) & 0xff); /* 1 */
	bigendian <<= 8;
	bigendian |= ((littleendian >>= 8) & 0xff); /* 2 */
	bigendian <<= 8;
	bigendian |= ((littleendian >>= 8) & 0xff); /* 3 */
	positions[i] = bigendian;
      }
#else
      for (ptr = ptr0, i = 0; ptr < end0; ptr++, i++) {
	littleendian = this->positions[ptr];
	bigendian = littleendian & 0xff; /* 0 */
	bigendian <<= 8;
	bigendian |= ((littleendian >>= 8) & 0xff); /* 1 */
	bigendian <<= 8;
	bigendian |= ((littleendian >>= 8) & 0xff); /* 2 */
	bigendian <<= 8;
	bigendian |= ((littleendian >>= 8) & 0xff); /* 3 */
	positions[i] = bigendian;
      }
#endif

#else  /* littleendian */

#ifdef LARGE_GENOMES
      positions_copy_multiple_large(positions,&(this->positions_high[ptr0]),&(this->positions_low[ptr0]),*nentries);
#else
      memcpy(positions,&(this->positions[ptr0]),(*nentries)*sizeof(Univcoord_T));
#endif
#endif
    }

    debug0(
	   printf("%d entries:",*nentries);
	   for (j = 0; j < *nentries; j++) {
	     printf(" %u",positions[j]);
	   }
	   printf("\n");
	   );

    return positions;
  }
}


/* GSNAP version.  Expects calling procedure to handle bigendian conversion. */
UINT4 *
Indexdb_read_inplace (int *nentries,
#ifdef LARGE_GENOMES
		      unsigned char **positions_high,
#endif
		      T this, Storedoligomer_T oligo) {
  Positionsptr_T ptr0, end0;
  Storedoligomer_T part0;
#ifdef DEBUG0
  Positionsptr_T ptr;
#endif

  debug0(printf("%08X (%s)\n",oligo,shortoligo_nt(oligo,index1part)));
  part0 = oligo & poly_T;	/* Probably unnecessary, since stage1 procedure already masks oligo */

  /* Needed to avoid overflow on 15-mers */
  if (part0 == poly_A || part0 == poly_T) {
    *nentries = 0;
    return NULL;
  }

  if (this->compression_type == NO_COMPRESSION) {
#ifdef WORDS_BIGENDIAN
    if (this->offsetscomp_access == ALLOCATED) {
      ptr0 = this->offsetscomp[part0];
      end0 = this->offsetscomp[part0+1];
    } else {
      ptr0 = Bigendian_convert_uint(this->offsetscomp[part0]);
      end0 = Bigendian_convert_uint(this->offsetscomp[part0+1]);
    }
#else
    ptr0 = this->offsetscomp[part0];
    end0 = this->offsetscomp[part0+1];
#endif

  } else if (this->compression_type == BITPACK64_COMPRESSION) {
#ifdef LARGE_GENOMES
    ptr0 = Bitpack64_offsetptr_huge(&end0,oligo,this->offsetspages,this->gammaptrs,this->offsetscomp);
#else
    ptr0 = Bitpack64_offsetptr(&end0,oligo,this->gammaptrs,this->offsetscomp);
#endif

  } else if (this->compression_type == GAMMA_COMPRESSION) {
#ifdef WORDS_BIGENDIAN
    if (this->offsetscomp_access == ALLOCATED) {
      ptr0 = Genome_offsetptr_from_gammas(&end0,this->gammaptrs,this->offsetscomp,this->offsetscomp_blocksize,oligo);
    } else {
      ptr0 = Genome_offsetptr_from_gammas_bigendian(&end0,this->gammaptrs,this->offsetscomp,this->offsetscomp_blocksize,oligo);
    }
#else
    ptr0 = Genome_offsetptr_from_gammas(&end0,this->gammaptrs,this->offsetscomp,this->offsetscomp_blocksize,oligo);
#endif

  } else {
    abort();
  }

  debug0(printf("Indexdb_read_inplace: offset pointers are %u and %u\n",ptr0,end0));

  *nentries = end0 - ptr0;

  if (*nentries == 0) {
    return NULL;
  } else if (this->positions_access == FILEIO) {
    abort();
  } else {
    debug0(
	   printf("%d entries:",*nentries);
	   for (ptr = ptr0; ptr < end0; ptr++) {
	     printf(" %lu",this->positions[ptr]);
	   }
	   printf("\n");
	   );
#ifdef LARGE_GENOMES
    *positions_high = &(this->positions_high[ptr0]);
    return &(this->positions_low[ptr0]);
#else
    return &(this->positions[ptr0]);
#endif
  }
}

#endif	/* ifdef PMAP */


/* Analogous to Indexdb_read, except this includes diagterm.  Always allocates memory. */
Univcoord_T *
Indexdb_read_with_diagterm (int *nentries, T this, Storedoligomer_T oligo, int diagterm) {
  Univcoord_T *positions;
  Positionsptr_T ptr0, end0, ptr;
  int i;

  if (this->compression_type == NO_COMPRESSION) {
#ifdef WORDS_BIGENDIAN
    if (this->offsetscomp_access == ALLOCATED) {
      ptr0 = this->offsetscomp[oligo];
      end0 = this->offsetscomp[oligo+1];
    } else {
      ptr0 = Bigendian_convert_uint(this->offsetscomp[oligo]);
      end0 = Bigendian_convert_uint(this->offsetscomp[oligo+1]);
    }
#else
    ptr0 = this->offsetscomp[oligo];
    end0 = this->offsetscomp[oligo+1];
#endif

  } else if (this->compression_type == BITPACK64_COMPRESSION) {
#ifdef LARGE_GENOMES
    ptr0 = Bitpack64_offsetptr_huge(&end0,oligo,this->offsetspages,this->gammaptrs,this->offsetscomp);
#else
    ptr0 = Bitpack64_offsetptr(&end0,oligo,this->gammaptrs,this->offsetscomp);
#endif

  } else if (this->compression_type == GAMMA_COMPRESSION) {
#ifdef WORDS_BIGENDIAN
    if (this->offsetscomp_access == ALLOCATED) {
      ptr0 = Genome_offsetptr_from_gammas(&end0,this->gammaptrs,this->offsetscomp,this->offsetscomp_blocksize,oligo);
    } else {
      ptr0 = Genome_offsetptr_from_gammas_bigendian(&end0,this->gammaptrs,this->offsetscomp,this->offsetscomp_blocksize,oligo);
    }
#else
    ptr0 = Genome_offsetptr_from_gammas(&end0,this->gammaptrs,this->offsetscomp,this->offsetscomp_blocksize,oligo);
#endif

  } else {
    abort();
  }

  debug0(printf("read_zero_shift: oligo = %06X, offset pointers are %u and %u\n",oligo,ptr0,end0));

  if ((*nentries = end0 - ptr0) == 0) {
    return (Univcoord_T *) NULL;
  } else {
    positions = (Univcoord_T *) CALLOC(*nentries,sizeof(Univcoord_T));
    if (this->positions_access == FILEIO) {
#ifdef HAVE_PTHREAD
      pthread_mutex_lock(&this->positions_read_mutex);
#endif
#ifdef LARGE_GENOMES
      positions_move_absolute_1(this->positions_high_fd,ptr0);
      positions_move_absolute_4(this->positions_low_fd,ptr0);
      positions_read_multiple_large(this->positions_high_fd,this->positions_low_fd,positions,*nentries);
#else
      positions_move_absolute_4(this->positions_fd,ptr0);
      positions_read_multiple(this->positions_fd,positions,*nentries);
#endif
#ifdef HAVE_PTHREAD
      pthread_mutex_unlock(&this->positions_read_mutex);
#endif

    } else if (this->positions_access == ALLOCATED) {
#ifdef LARGE_GENOMES
      for (ptr = ptr0, i = 0; ptr < end0; ptr++) {
	positions[i++] = ((Univcoord_T) this->positions_high[ptr] << 32) + this->positions_low[ptr] + diagterm;
      }
#else
      for (ptr = ptr0, i = 0; ptr < end0; ptr++) {
	positions[i++] = this->positions[ptr] + diagterm;
      }
#endif

    } else {

#ifdef WORDS_BIGENDIAN
#ifdef LARGE_GENOMES
      for (ptr = ptr0, i = 0; ptr < end0; ptr++) {
	positions[i++] = ((Univcoord_T) this->positions_high[ptr] << 32) + Bigendian_convert_uint(this->positions_low[ptr]) + diagterm;
      }
#else
      for (ptr = ptr0, i = 0; ptr < end0; ptr++) {
	positions[i++] = Bigendian_convert_univcoord(this->positions[ptr]) + diagterm;
      }
#endif
#else
#ifdef LARGE_GENOMES
      for (ptr = ptr0, i = 0; ptr < end0; ptr++) {
	positions[i++] = ((Univcoord_T) this->positions_high[ptr] << 32) + this->positions_low[ptr] + diagterm;
      }
#else
      for (ptr = ptr0, i = 0; ptr < end0; ptr++) {
	positions[i++] = this->positions[ptr] + diagterm;
      }
#endif
#endif
    }
  }
      
  debug0(
	printf("%d entries:",*nentries);
	for (i = 0; i < *nentries; i++) {
	  printf(" %u",positions[i]);
	}
	printf("\n");
	);
  
  return positions;
}


/* Analogous to Indexdb_read, except this includes diagterm.  Always allocates memory. */
Univcoord_T *
Indexdb_read_with_diagterm_sizelimit (int *nentries, T this, Storedoligomer_T oligo, int diagterm,
				      int size_threshold) {
  Univcoord_T *positions;
  Positionsptr_T ptr0, end0, ptr;
  int i;

  if (this->compression_type == NO_COMPRESSION) {
#ifdef WORDS_BIGENDIAN
    if (this->offsetscomp_access == ALLOCATED) {
      ptr0 = this->offsetscomp[oligo];
      end0 = this->offsetscomp[oligo+1];
    } else {
      ptr0 = Bigendian_convert_uint(this->offsetscomp[oligo]);
      end0 = Bigendian_convert_uint(this->offsetscomp[oligo+1]);
    }
#else
    ptr0 = this->offsetscomp[oligo];
    end0 = this->offsetscomp[oligo+1];
#endif

  } else if (this->compression_type == BITPACK64_COMPRESSION) {
#ifdef LARGE_GENOMES
    ptr0 = Bitpack64_offsetptr_huge(&end0,oligo,this->offsetspages,this->gammaptrs,this->offsetscomp);
#else
    ptr0 = Bitpack64_offsetptr(&end0,oligo,this->gammaptrs,this->offsetscomp);
#endif

  } else if (this->compression_type == GAMMA_COMPRESSION) {
#ifdef WORDS_BIGENDIAN
    if (this->offsetscomp_access == ALLOCATED) {
      ptr0 = Genome_offsetptr_from_gammas(&end0,this->gammaptrs,this->offsetscomp,this->offsetscomp_blocksize,oligo);
    } else {
      ptr0 = Genome_offsetptr_from_gammas_bigendian(&end0,this->gammaptrs,this->offsetscomp,this->offsetscomp_blocksize,oligo);
    }
#else
    ptr0 = Genome_offsetptr_from_gammas(&end0,this->gammaptrs,this->offsetscomp,this->offsetscomp_blocksize,oligo);
#endif

  } else {
    abort();
  }


  debug0(printf("read_zero_shift: oligo = %06X, offset pointers are %u and %u\n",oligo,ptr0,end0));

  if ((*nentries = end0 - ptr0) == 0) {
    return (Univcoord_T *) NULL;

  } else if (*nentries > size_threshold) {
    *nentries = 0;
    return (Univcoord_T *) NULL;

  } else {
    positions = (Univcoord_T *) CALLOC(*nentries,sizeof(Univcoord_T));
    if (this->positions_access == FILEIO) {
#ifdef HAVE_PTHREAD
      pthread_mutex_lock(&this->positions_read_mutex);
#endif
#ifdef LARGE_GENOMES
      positions_move_absolute_1(this->positions_high_fd,ptr0);
      positions_move_absolute_4(this->positions_low_fd,ptr0);
      positions_read_multiple_large(this->positions_high_fd,this->positions_low_fd,positions,*nentries);
#else
      positions_move_absolute_4(this->positions_fd,ptr0);
      positions_read_multiple(this->positions_fd,positions,*nentries);
#endif

#ifdef HAVE_PTHREAD
      pthread_mutex_unlock(&this->positions_read_mutex);
#endif

    } else if (this->positions_access == ALLOCATED) {
#ifdef LARGE_GENOMES
      for (ptr = ptr0, i = 0; ptr < end0; ptr++) {
	positions[i++] = ((Univcoord_T) this->positions_high[ptr] << 32) + this->positions_low[ptr] + diagterm;
      }
#else
      for (ptr = ptr0, i = 0; ptr < end0; ptr++) {
	positions[i++] = this->positions[ptr] + diagterm;
      }
#endif

    } else {

#ifdef WORDS_BIGENDIAN
#ifdef LARGE_GENOMES
      for (ptr = ptr0, i = 0; ptr < end0; ptr++) {
	positions[i++] = ((Univcoord_T) this->positions_high[ptr] << 32) + Bigendian_convert_uint(this->positions_low[ptr]) + diagterm;
      }
#else
      for (ptr = ptr0, i = 0; ptr < end0; ptr++) {
	positions[i++] = Bigendian_convert_univcoord(this->positions[ptr]) + diagterm;
      }
#endif
#else
#ifdef LARGE_GENOMES
      for (ptr = ptr0, i = 0; ptr < end0; ptr++) {
	positions[i++] = ((Univcoord_T) this->positions_high[ptr] << 32) + this->positions_low[ptr] + diagterm;
      }
#else
      for (ptr = ptr0, i = 0; ptr < end0; ptr++) {
	positions[i++] = this->positions[ptr] + diagterm;
      }
#endif
#endif
    }
  }
      
  debug0(
	printf("%d entries:",*nentries);
	for (i = 0; i < *nentries; i++) {
	  printf(" %u",positions[i]);
	}
	printf("\n");
	);
  
  return positions;
}


/************************************************************************
 *   Create procedure -- for user-provided genomic segment
 ************************************************************************/

#if defined(UTILITYP) || defined(LARGE_GENOMES)
#else
T
Indexdb_new_segment (char *genomicseg,
#ifdef PMAP
		     int alphabet_size, Width_T index1part_aa, bool watsonp,
#else
		     Width_T index1part,
#endif
		     Width_T index1interval) {
  T new = (T) MALLOC(sizeof(*new));
  char *uppercaseCode;
  Positionsptr_T *work_offsets;	/* Working set for use in calculating positions */
  int totalcounts = 0;
  int c;
  Oligospace_T oligospace, oligoi;
  char *p;
  Univcoord_T position = 0UL;
#ifdef PMAP
  int frame = -1, between_counter[3], in_counter[3];
  Storedoligomer_T high = 0U, low = 0U, carry;
  Storedoligomer_T aaindex;
  int index1part_nt = 3*index1part_aa;
#else
  int between_counter = 0, in_counter = 0;
  Storedoligomer_T oligo = 0U, masked, mask;
#endif

  uppercaseCode = UPPERCASE_U2T;

#ifdef PMAP
  oligospace = power(alphabet_size,index1part_aa);
  between_counter[0] = between_counter[1] = between_counter[2] = 0;
  in_counter[0] = in_counter[1] = in_counter[2] = 0;
  new->index1part = index1part_aa;
#else
  mask = ~(~0UL << 2*index1part);
  oligospace = power(4,index1part);
  new->index1part = index1part;
#endif
  new->index1interval = 1;

  new->compression_type = NO_COMPRESSION;

  new->gammaptrs = (Gammaptr_T *) CALLOC(oligospace+1,sizeof(Gammaptr_T));
  for (oligoi = 0; oligoi <= oligospace; oligoi++) {
    new->gammaptrs[oligoi] = oligoi;
  }

  new->offsetscomp = (Offsetscomp_T *) CALLOC(oligospace+1,sizeof(Offsetscomp_T));

  p = genomicseg;
  while ((c = *(p++)) != '\0') {
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

#ifdef PMAP
    carry = (low >> 30);
    switch (uppercaseCode[c]) {
    case 'A': low = (low << 2); break;
    case 'C': low = (low << 2) | 1U; break;
    case 'G': low = (low << 2) | 2U; break;
    case 'T': low = (low << 2) | 3U; break;
    default:
      high = low = carry = 0U;
      in_counter[0] = in_counter[1] = in_counter[2] = 0;
      break;
    }
    high = (high << 2) | carry; 
#else
    switch (uppercaseCode[c]) {
    case 'A': oligo = (oligo << 2); break;
    case 'C': oligo = (oligo << 2) | 1U; break;
    case 'G': oligo = (oligo << 2) | 2U; break;
    case 'T': oligo = (oligo << 2) | 3U; break;
    default: oligo = 0U; in_counter = 0; break;
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
#ifdef OLIGOSPACE_NOT_LONG
	oligoi = (Oligospace_T) aaindex + 1U;
#else
	oligoi = (Oligospace_T) aaindex + 1UL;
#endif
	new->offsetscomp[oligoi] += 1;
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
	  /* Actually, modular condition not needed for user-supplied genomic segment */
	  (position-index1part+1U) % index1interval == 0
#endif
	  ) {
	masked = oligo & mask;
#ifdef OLIGOSPACE_NOT_LONG
	oligoi = (Oligospace_T) masked + 1U;
#else
	oligoi = (Oligospace_T) masked + 1UL;
#endif
	new->offsetscomp[oligoi] += 1;
	debug(printf("Found oligo %06X.  Incremented offsets for %lu to be %u\n",
		     masked,oligoi,new->offsetscomp[oligoi]));
	between_counter = 0;
      }
      in_counter--;
    }
#endif

    position++;
  }

#ifdef ADDSENTINEL
  for (oligoi = 1; oligoi <= oligospace; oligoi++) {
    new->offsetscomp[oligoi] = new->offsetscomp[oligoi] + new->offsetscomp[oligoi-1] + 1;
    debug(printf("Offset for %06X: %u\n",oligoi,new->offsetscomp[oligoi]));
  }
#else
  for (oligoi = 1; oligoi <= oligospace; oligoi++) {
    new->offsetscomp[oligoi] = new->offsetscomp[oligoi] + new->offsetscomp[oligoi-1];
    debug(printf("Offset for %06X: %u\n",oligoi,new->offsetscomp[oligoi]));
  }
#endif


  /* Create positions */
  position = 0U;
#ifdef PMAP
  frame = -1;
  between_counter[0] = between_counter[1] = between_counter[2] = 0;
  in_counter[0] = in_counter[1] = in_counter[2] = 0;
  high = low = 0U;
#else
  between_counter = in_counter = 0;
  oligo = 0UL;
#endif

  work_offsets = (Positionsptr_T *) CALLOC(oligospace+1,sizeof(Positionsptr_T));
  for (oligoi = 0; oligoi <= oligospace; oligoi++) {
    work_offsets[oligoi] = new->offsetscomp[oligoi];
  }

  totalcounts = new->offsetscomp[oligospace];
  if (totalcounts == 0) {
#ifdef PMAP
    fprintf(stderr,"Error: user-provided genomic segment has no valid oligomers of size %d\n",index1part_nt);
#else
    fprintf(stderr,"Error: user-provided genomic segment has no valid oligomers of size %d\n",index1part);
#endif
    exit(9);
  }
  new->positions = (Univcoord_T *) CALLOC(totalcounts,sizeof(Univcoord_T));
  new->positions_access = ALLOCATED;

  p = genomicseg;
  while ((c = *(p++)) != '\0') {
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

#ifdef PMAP
    carry = (low >> 30);
    switch (uppercaseCode[c]) {
    case 'A': low = (low << 2); break;
    case 'C': low = (low << 2) | 1; break;
    case 'G': low = (low << 2) | 2; break;
    case 'T': low = (low << 2) | 3; break;
    default:
      high = low = carry = 0U;
      in_counter[0] = in_counter[1] = in_counter[2] = 0;
      break;
    }
    high = (high << 2) | carry; 
#else
    switch (uppercaseCode[c]) {
    case 'A': oligo = (oligo << 2); break;
    case 'C': oligo = (oligo << 2) | 1; break;
    case 'G': oligo = (oligo << 2) | 2; break;
    case 'T': oligo = (oligo << 2) | 3; break;
    default: oligo = 0U; in_counter = 0; break;
    }
#endif

    /*
    debug(printf("char=%c bc=%d ic=%d oligo=%06X\n",
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
	if (watsonp == true) {
	  new->positions[work_offsets[aaindex]++] = position-index1part_nt+1;
	} else {
	  new->positions[work_offsets[aaindex]++] = position;
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
	  /* Actually, modular condition not needed for user-supplied genomic segment */
	  (position-index1part+1) % index1interval == 0
#endif
	  ) {
	masked = oligo & mask;
	new->positions[work_offsets[masked]++] = position-index1part+1;
	between_counter = 0;
      }
      in_counter--;
    }
#endif

    position++;
  }

#ifdef ADDSENTINEL
  for (oligoi = 0; oligoi < oligospace; oligoi++) {
    new->positions[work_offsets[oligoi]] = (Univcoord_T) -1;
  }
#endif

  FREE(work_offsets);

  return new;
}
#endif


int
Storedoligomer_compare (const void *a, const void *b) {
  Storedoligomer_T x = * (Storedoligomer_T *) a;
  Storedoligomer_T y = * (Storedoligomer_T *) b;

  if (x < y) {
    return -1;
  } else if (y < x) {
    return 1;
  } else {
    return 0;
  }
}

