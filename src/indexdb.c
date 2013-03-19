static char rcsid[] = "$Id: indexdb.c 88195 2013-03-06 21:42:45Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifndef HAVE_MEMCPY
# define memcpy(d,s,n) bcopy((s),(d),(n))
#endif
#ifndef HAVE_MEMMOVE
# define memmove(d,s,n) bcopy((s),(d),(n))
#endif

#include "indexdb.h"
#include "indexdbdef.h"
#include "genome_hr.h"		/* For read_gammas procedures */


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

/* Writing of positions */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

/* Shifting of high word to low word for PMAP */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

/* Gammas */
#ifdef DEBUG3
#define debug3(x) x
#else
#define debug3(x)
#endif


#ifdef PMAP

#if (defined(DEBUG) || defined(DEBUG0) || defined(DEBUG1))
static int index1part_aa;
#endif

void
Indexdb_setup (int index1part_aa_in) {
#if (defined(DEBUG) || defined(DEBUG0) || defined(DEBUG1))
  index1part_aa = index1part_aa_in;
#endif
  return;
}

#else

#define poly_A 0U
static unsigned int poly_T;  /* Was LOW12MER 0x00FFFFFF */

#if (defined(DEBUG) || defined(DEBUG0) || defined(DEBUG1))
static int index1part;
#endif

void
Indexdb_setup (int index1part_in) {
#if (defined(DEBUG) || defined(DEBUG0) || defined(DEBUG1))
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
      FREE((*old)->positions);
#ifdef HAVE_MMAP
    } else if ((*old)->positions_access == MMAPPED) {
      munmap((void *) (*old)->positions,(*old)->positions_len);
      close((*old)->positions_fd);
#endif
    } else if ((*old)->positions_access == FILEIO) {
#ifdef HAVE_PTHREAD
      pthread_mutex_destroy(&(*old)->positions_read_mutex);
#endif
      close((*old)->positions_fd);
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


int
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

double
Indexdb_mean_size (T this, Mode_T mode, int index1part) {
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



#define LARGEVALUE 1000000

bool
Indexdb_get_filenames_pregamma (char **offsets_filename, char **positions_filename,
				char **offsets_basename_ptr, char **positions_basename_ptr,
				char **offsets_index1info_ptr, char **positions_index1info_ptr,
				int *index1part, int *index1interval, char *genomesubdir,
				char *fileroot, char *idx_filesuffix, char *snps_root,
				int required_interval) {
  char *base_filename, *filename;
  char *pattern, interval_char, digit_string[2], *p, *q;
  int found_index1part, found_interval;
  int rootlength, patternlength;

  char *offsets_suffix, *positions_suffix;
  struct dirent *entry;
  DIR *dp;


  if (snps_root == NULL) {
    offsets_suffix = "offsets";
    positions_suffix = POSITIONS_FILESUFFIX;
  } else {
    offsets_suffix = (char *) CALLOC(strlen("offsets")+strlen(".")+strlen(snps_root)+1,sizeof(char));
    sprintf(offsets_suffix,"%s.%s","offsets",snps_root);
    positions_suffix = (char *) CALLOC(strlen(POSITIONS_FILESUFFIX)+strlen(".")+strlen(snps_root)+1,sizeof(char));
    sprintf(positions_suffix,"%s.%s",POSITIONS_FILESUFFIX,snps_root);
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
	if (q - p == 1) {
	  /* Old style, e.g, idx or ref3 */
	  if (sscanf(p,"%c",&interval_char) == 1) {
	    if (interval_char == 'x') {
	      found_interval = 6;
	    } else {
	      digit_string[0] = interval_char;
	      found_interval = atoi(digit_string);
	    }
	  }
#ifdef PMAP
	  found_index1part = found_interval;
#else
	  found_index1part = 12;
#endif

	} else {
	  fprintf(stderr,"Cannot parse part between %s and offsets in filename %s\n",idx_filesuffix,filename);
	  return false;
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
    fprintf(stderr,"Cannot find offsets file containing %s and %s",
	    idx_filesuffix,offsets_suffix);
    if (required_interval > 0) {
      fprintf(stderr," and having sampling interval of %d",required_interval);
    }
    fprintf(stderr,"\n");

    *offsets_filename = (char *) NULL;
    *positions_filename = (char *) NULL;
    return false;

  } else {
    *offsets_filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+strlen(base_filename)+1,sizeof(char));
    *offsets_basename_ptr = &((*offsets_filename)[strlen(genomesubdir)+strlen("/")]);
    *offsets_index1info_ptr = &((*offsets_basename_ptr)[patternlength]);

    sprintf(*offsets_filename,"%s/%s",genomesubdir,base_filename);
    if (Access_file_exists_p(*offsets_filename) == false) {
      fprintf(stderr,"Offsets filename %s does not exist\n",*offsets_filename);
      FREE(*offsets_filename);
      *offsets_filename = (char *) NULL;
      *positions_filename = (char *) NULL;
      FREE(base_filename);
      return false;
    }


    if ((q = strstr(base_filename,offsets_suffix)) == NULL) {
      abort();
    } else {
      rootlength = q - base_filename;
    }

    *positions_filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+rootlength+strlen(positions_suffix)+1,sizeof(char));
    *positions_basename_ptr = &((*positions_filename)[strlen(genomesubdir)+strlen("/")]);
    *positions_index1info_ptr = &((*positions_basename_ptr)[patternlength]);

    sprintf(*positions_filename,"%s/",genomesubdir);
    strncpy(*positions_basename_ptr,base_filename,rootlength);
    strcpy(&((*positions_basename_ptr)[rootlength]),positions_suffix);

    if (Access_file_exists_p(*positions_filename) == false) {
      fprintf(stderr,"Positions filename %s does not exist\n",*positions_filename);
      FREE(*offsets_filename);
      FREE(*positions_filename);
      *offsets_filename = (char *) NULL;
      *positions_filename = (char *) NULL;
      FREE(base_filename);
      return false;
    }

    if (snps_root != NULL) {
      FREE(offsets_suffix);
      FREE(positions_suffix);
    }

    FREE(base_filename);
    fprintf(stderr,"Looking for index files in directory %s (offsets not compressed)\n",genomesubdir);
    fprintf(stderr,"  Offsets file is %s\n",*offsets_basename_ptr);
    fprintf(stderr,"  Positions file is %s\n",*positions_basename_ptr);
    return true;
  }
}



#ifdef PMAP
#define BASE_KMER_SAMPLING 3   /* e.g., 677 */
#define KMER_SAMPLING 2   /* e.g., 77 */
#else
#define BASE_KMER_SAMPLING 5   /* e.g., 12153 */
#define KMER_SAMPLING 3   /* e.g., 153 */
#endif


bool
Indexdb_get_filenames (char **gammaptrs_filename, char **offsetscomp_filename, char **positions_filename,
		       char **gammaptrs_basename_ptr, char **offsetscomp_basename_ptr, char **positions_basename_ptr,
		       char **gammaptrs_index1info_ptr, char **offsetscomp_index1info_ptr, char **positions_index1info_ptr,
#ifdef PMAP
		       Alphabet_T *alphabet, Alphabet_T required_alphabet,
#endif
		       int *basesize, int *index1part, int *index1interval, char *genomesubdir,
		       char *fileroot, char *idx_filesuffix, char *snps_root,
		       int required_basesize, int required_index1part, int required_interval) {
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
  int found_basesize, found_index1part, found_interval;
  int rootlength, patternlength;

  char ones0, ones;
  char *gammaptrs_suffix, *offsetscomp_suffix, *positions_suffix;
  struct dirent *entry;
  DIR *dp;


  if (snps_root == NULL) {
    gammaptrs_suffix = "gammaptrs";
    offsetscomp_suffix = "offsetscomp";
    positions_suffix = POSITIONS_FILESUFFIX;
  } else {
    gammaptrs_suffix = (char *) CALLOC(strlen("gammaptrs")+strlen(".")+strlen(snps_root)+1,sizeof(char));
    sprintf(gammaptrs_suffix,"%s.%s","gammaptrs",snps_root);
    offsetscomp_suffix = (char *) CALLOC(strlen("offsetscomp")+strlen(".")+strlen(snps_root)+1,sizeof(char));
    sprintf(offsetscomp_suffix,"%s.%s","offsetscomp",snps_root);
    positions_suffix = (char *) CALLOC(strlen(POSITIONS_FILESUFFIX)+strlen(".")+strlen(snps_root)+1,sizeof(char));
    sprintf(positions_suffix,"%s.%s",POSITIONS_FILESUFFIX,snps_root);
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
	    }
	  } else {
	    /* fprintf(stderr,"Cannot parse part between %s and offsets in filename %s\n",idx_filesuffix,filename); */
	    return false;
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
	  }
	} else {
	  fprintf(stderr,"Cannot parse part between %s and offsets in filename %s: found %ld characters, expecting %d\n",
		  idx_filesuffix,filename,q-p,BASE_KMER_SAMPLING);
	  return false;
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
    fprintf(stderr,"Cannot find offsetscomp file containing %s and %s",
	    idx_filesuffix,offsetscomp_suffix);
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

    *gammaptrs_filename = (char *) NULL;
    *offsetscomp_filename = (char *) NULL;
    *positions_filename = (char *) NULL;
    return false;

  } else {
    *offsetscomp_filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+strlen(base_filename)+1,sizeof(char));
    *offsetscomp_basename_ptr = &((*offsetscomp_filename)[strlen(genomesubdir)+strlen("/")]);
    *offsetscomp_index1info_ptr = &((*offsetscomp_basename_ptr)[patternlength]);

    sprintf(*offsetscomp_filename,"%s/%s",genomesubdir,base_filename);
    if (Access_file_exists_p(*offsetscomp_filename) == false) {
      fprintf(stderr,"Offsets filename %s does not exist\n",*offsetscomp_filename);
      FREE(*offsetscomp_filename);
      *offsetscomp_filename = (char *) NULL;
      *positions_filename = (char *) NULL;
      FREE(base_filename);
      return false;
    }


    if ((q = strstr(base_filename,offsetscomp_suffix)) == NULL) {
      abort();
    } else {
      rootlength = q - base_filename;
    }

    if (*index1part == *basesize) {
      *gammaptrs_filename = (char *) NULL;
      *gammaptrs_basename_ptr = (char *) NULL;
      *gammaptrs_index1info_ptr = (char *) NULL;
    } else {
      *gammaptrs_filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+rootlength+strlen(gammaptrs_suffix)+1,sizeof(char));
      *gammaptrs_basename_ptr = &((*gammaptrs_filename)[strlen(genomesubdir)+strlen("/")]);
      *gammaptrs_index1info_ptr = &((*gammaptrs_basename_ptr)[patternlength]);

      sprintf(*gammaptrs_filename,"%s/",genomesubdir);
      strncpy(*gammaptrs_basename_ptr,base_filename,rootlength);
      strcpy(&((*gammaptrs_basename_ptr)[rootlength]),gammaptrs_suffix);

      if (Access_file_exists_p(*gammaptrs_filename) == false) {
	fprintf(stderr,"Gammaptrs filename %s does not exist\n",*gammaptrs_filename);
	FREE(*offsetscomp_filename);
	*offsetscomp_filename = (char *) NULL;
	FREE(base_filename);
	return false;
      }
    }


    *positions_filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+rootlength+strlen(positions_suffix)+1,sizeof(char));
    *positions_basename_ptr = &((*positions_filename)[strlen(genomesubdir)+strlen("/")]);
    *positions_index1info_ptr = &((*positions_basename_ptr)[patternlength]);

    sprintf(*positions_filename,"%s/",genomesubdir);
    strncpy(*positions_basename_ptr,base_filename,rootlength);
    strcpy(&((*positions_basename_ptr)[rootlength]),positions_suffix);

    if (Access_file_exists_p(*positions_filename) == false) {
      /* Try newer naming scheme: ref153positions instead of ref12153positions */
      sprintf(*positions_filename,"%s/",genomesubdir);
      strncpy(*positions_basename_ptr,base_filename,rootlength-BASE_KMER_SAMPLING); /* e.g., skip "12153" */
      strncpy(&((*positions_basename_ptr)[rootlength-BASE_KMER_SAMPLING]),&(base_filename[rootlength-KMER_SAMPLING]),KMER_SAMPLING);
      strcpy(&((*positions_basename_ptr)[rootlength+KMER_SAMPLING-BASE_KMER_SAMPLING]),positions_suffix);

      if (Access_file_exists_p(*positions_filename) == false) {
	fprintf(stderr,"Positions filename %s does not exist\n",*positions_filename);
	FREE(*gammaptrs_filename);
	FREE(*offsetscomp_filename);
	FREE(*positions_filename);
	*gammaptrs_filename = (char *) NULL;
	*offsetscomp_filename = (char *) NULL;
	*positions_filename = (char *) NULL;
	FREE(base_filename);
	return false;
      }
    }

    if (snps_root != NULL) {
      FREE(offsetscomp_suffix);
      FREE(gammaptrs_suffix);
      FREE(positions_suffix);
    }

    FREE(base_filename);

    fprintf(stderr,"Looking for index files in directory %s\n",genomesubdir);
    if (*gammaptrs_filename == NULL) {
      fprintf(stderr,"  No gammaptrs file, because kmersize %d == basesize %d\n",
	      *index1part,*basesize);
    } else {
      fprintf(stderr,"  Gammaptrs file is %s\n",*gammaptrs_basename_ptr);
    }
    fprintf(stderr,"  Offsetscomp file is %s\n",*offsetscomp_basename_ptr);
    fprintf(stderr,"  Positions file is %s\n",*positions_basename_ptr);
    return true;
  }
}



T
Indexdb_new_genome (int *basesize, int *index1part, int *index1interval,
		    char *genomesubdir, char *fileroot, char *idx_filesuffix, char *snps_root,
#ifdef PMAP
		    Alphabet_T *alphabet, int *alphabet_size, Alphabet_T required_alphabet,
#endif
		    int required_basesize, int required_index1part, int required_interval, bool expand_offsets_p,
		    Access_mode_T offsetscomp_access, Access_mode_T positions_access) {
  T new = (T) MALLOC(sizeof(*new));
  char *gammaptrs_filename, *offsetscomp_filename, *positions_filename,
    *gammaptrs_basename_ptr, *offsetscomp_basename_ptr, *positions_basename_ptr,
    *gammaptrs_index1info_ptr, *offsetscomp_index1info_ptr, *positions_index1info_ptr;
  char *offsets_filename, *offsets_basename_ptr, *offsets_index1info_ptr;
  Oligospace_T basespace, base;

  char *comma;
  double seconds;
#ifdef HAVE_MMAP
  int npages;
#endif

  /* Read offsets file */
  if (Indexdb_get_filenames(&gammaptrs_filename,&offsetscomp_filename,&positions_filename,
			    &gammaptrs_basename_ptr,&offsetscomp_basename_ptr,&positions_basename_ptr,
			    &gammaptrs_index1info_ptr,&offsetscomp_index1info_ptr,&positions_index1info_ptr,
#ifdef PMAP
			    &(*alphabet),required_alphabet,
#endif
			    &new->offsetscomp_basesize,&new->index1part,&new->index1interval,
			    genomesubdir,fileroot,idx_filesuffix,snps_root,
			    required_basesize,required_index1part,required_interval) == true) {
    *index1part = new->index1part;
    *index1interval = new->index1interval;

    if (expand_offsets_p == true) {
      *basesize = *index1part;
#ifdef PMAP
      *alphabet_size = Alphabet_get_size(*alphabet);
      new->offsetscomp_blocksize = 1;
      basespace = power(*alphabet_size,new->index1part);
#else
      new->offsetscomp_blocksize = 1;
      basespace = power(4,new->index1part);
#endif
      new->gammaptrs = (UINT4 *) CALLOC(basespace+1,sizeof(UINT4));
      for (base = 0; base <= basespace; base++) {
	new->gammaptrs[base] = base;
      }

#ifdef PMAP
      new->offsetscomp = Indexdb_offsets_from_gammas(gammaptrs_filename,offsetscomp_filename,
						     new->offsetscomp_basesize,*alphabet_size,new->index1part);
#else
      new->offsetscomp = Indexdb_offsets_from_gammas(gammaptrs_filename,offsetscomp_filename,
						     new->offsetscomp_basesize,new->index1part);
#endif
      new->offsetscomp_access = ALLOCATED;

    } else {
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
	new->gammaptrs = (UINT4 *) CALLOC(basespace+1,sizeof(UINT4));
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
	new->gammaptrs = (UINT4 *) Access_allocated(&new->gammaptrs_len,&seconds,
						    gammaptrs_filename,sizeof(UINT4));
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
	new->offsetscomp = (UINT4 *) Access_allocated(&new->offsetscomp_len,&seconds,
						      offsetscomp_filename,sizeof(Positionsptr_T));
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
	new->offsetscomp = (UINT4 *) Access_mmap_and_preload(&new->offsetscomp_fd,&new->offsetscomp_len,&npages,&seconds,
							     offsetscomp_filename,sizeof(Positionsptr_T));
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
	new->offsetscomp = (UINT4 *) Access_mmap(&new->offsetscomp_fd,&new->offsetscomp_len,
						 offsetscomp_filename,sizeof(Positionsptr_T),/*randomp*/false);
	if (new->offsetscomp == NULL) {
	  fprintf(stderr,"Insufficient memory for mmap of %s (will use disk file instead, but program may not run)\n",
		  offsetscomp_filename);
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
	fprintf(stderr,"Offsetscomp file I/O access of %s not allowed\n",offsetscomp_filename);
	exit(9);

      } else {
	fprintf(stderr,"Don't recognize offsetscomp_access type %d\n",offsetscomp_access);
	abort();
      }
    }


    FREE(offsetscomp_filename);
    FREE(gammaptrs_filename);


  } else if (Indexdb_get_filenames_pregamma(&offsets_filename,&positions_filename,
					    &offsets_basename_ptr,&positions_basename_ptr,
					    &offsets_index1info_ptr,&positions_index1info_ptr,
					    &new->index1part,&new->index1interval,
					    genomesubdir,fileroot,idx_filesuffix,snps_root,
					    required_interval) == true) {

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
    new->gammaptrs = (UINT4 *) CALLOC(basespace+1,sizeof(UINT4));
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
      new->offsetscomp = (UINT4 *) Access_allocated(&new->offsetscomp_len,&seconds,
						    offsets_filename,sizeof(Positionsptr_T));
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
      new->offsetscomp = (UINT4 *) Access_mmap_and_preload(&new->offsetscomp_fd,&new->offsetscomp_len,&npages,&seconds,
							   offsets_filename,sizeof(Positionsptr_T));
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
      new->offsetscomp = (UINT4 *) Access_mmap(&new->offsetscomp_fd,&new->offsetscomp_len,
					       offsets_filename,sizeof(Positionsptr_T),/*randomp*/false);
      if (new->offsetscomp == NULL) {
	fprintf(stderr,"Insufficient memory for mmap of %s (will use disk file instead, but program may not run)\n",
		offsets_filename);
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
      fprintf(stderr,"Offsets file I/O access of %s not allowed\n",offsets_filename);
      exit(9);

    } else {
      fprintf(stderr,"Don't recognize offsetscomp_access type %d\n",offsetscomp_access);
      abort();
    }

    FREE(offsets_filename);

  } else {
    fprintf(stderr,"Cannot find genomic index files in either current or old format\n");
    exit(9);
  }


  /* Positions */

  if (positions_access == USE_ALLOCATE) {
    if (snps_root) {
      fprintf(stderr,"Allocating memory for %s (%s) positions, kmer %d, interval %d...",
	      idx_filesuffix,snps_root,new->index1part,new->index1interval);
    } else {
      fprintf(stderr,"Allocating memory for %s positions, kmer %d, interval %d...",
	      idx_filesuffix,new->index1part,new->index1interval);
    }
    new->positions = (Genomicpos_T *) Access_allocated(&new->positions_len,&seconds,
						       positions_filename,sizeof(Genomicpos_T));
    if (new->positions == NULL) {
      fprintf(stderr,"insufficient memory (need to use a lower batch mode (-B)\n");
      exit(9);
    } else {
      comma = Genomicpos_commafmt(new->positions_len);
      fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma,seconds);
      FREE(comma);
      new->positions_access = ALLOCATED;
    }

#ifdef HAVE_MMAP
  } else if (positions_access == USE_MMAP_PRELOAD) {
    if (snps_root) {
      fprintf(stderr,"Pre-loading %s (%s) positions, kmer %d, interval %d...",
	      idx_filesuffix,snps_root,new->index1part,new->index1interval);
    } else {
      fprintf(stderr,"Pre-loading %s positions, kmer %d, interval %d...",
	      idx_filesuffix,new->index1part,new->index1interval);
    }
    new->positions = (Genomicpos_T *) Access_mmap_and_preload(&new->positions_fd,&new->positions_len,&npages,&seconds,
							    positions_filename,sizeof(Genomicpos_T));
    if (new->positions == NULL) {
      fprintf(stderr,"insufficient memory (will use disk file instead, but program will be slow)\n");
      new->positions_access = FILEIO;
    } else {
      comma = Genomicpos_commafmt(new->positions_len);
      fprintf(stderr,"done (%s bytes, %d pages, %.2f sec)\n",comma,npages,seconds);
      FREE(comma);
      new->positions_access = MMAPPED;
    }

  } else if (positions_access == USE_MMAP_ONLY) {
    new->positions = (Genomicpos_T *) Access_mmap(&new->positions_fd,&new->positions_len,
						  positions_filename,sizeof(Genomicpos_T),/*randomp*/true);
    if (new->positions == NULL) {
      fprintf(stderr,"Insufficient memory for mmap of %s (will use disk file instead, but program will be slow)\n",
	      positions_filename);
      new->positions_access = FILEIO;
    } else {
      new->positions_access = MMAPPED;
    }
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

  FREE(positions_filename);

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
#endif

#endif


/************************************************************************
 *   Read procedures
 ************************************************************************/


static void
positions_move_absolute (int positions_fd, Positionsptr_T ptr) {
  off_t offset = ptr*((off_t) sizeof(Genomicpos_T));

  if (lseek(positions_fd,offset,SEEK_SET) < 0) {
    fprintf(stderr,"Attempted to do lseek on offset %u*%lu=%lu\n",
	    ptr,sizeof(Genomicpos_T),(long unsigned int) offset);
    perror("Error in indexdb.c, positions_move_absolute");
    exit(9);
  }
  return;
}

#if 0
static Genomicpos_T
positions_read_forward (int positions_fd) {
  Genomicpos_T value;
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

static void
positions_read_multiple (int positions_fd, Genomicpos_T *values, int n) {
  int i;
  Genomicpos_T value;
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

#if 0
static Genomicpos_T
positions_read_backward (int positions_fd) {
  Genomicpos_T value;
  char buffer[4];
  off_t reloffset = -2*((off_t) sizeof(Genomicpos_T)); /* 1 to undo the effect of read */

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


/************************************************************************
 *   Elias gamma representation
 ************************************************************************/

/* ctr is 32 at high bit and 1 at low bit */
static int
write_gamma (FILE *offsetscomp_fp, UINT4 *offsets_buffer, int offsets_buffer_size, int *offsets_buffer_i,
	     unsigned int *nwritten, unsigned int *buffer, int ctr, unsigned int gamma) {
  int length;
  unsigned int nn;
  
  debug3(printf("Entering write_gamma with gamma %u, ctr %d\n",gamma,ctr));

  gamma += 1;			/* To allow 0 to be represented */

  /* Compute length */
  length = 1;
  nn = 2;
  while (nn <= gamma) {
    length += 2;
    nn += nn;
  }
  debug3(printf("gamma is %u (%08X), length is %u\n",gamma,gamma,length));


  /* Update buffer and write */
  while (length > ctr) {
    if (length - ctr < 32) {
      *buffer |= (gamma >> (length - ctr));
    }
    debug3(printf("writing gamma %08X\n",*buffer));
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
  
  debug3(printf("  shifting gamma left by %d\n",ctr - length));
  *buffer |= (gamma << (ctr - length));
  debug3(printf("  buffer is %08X\n",*buffer));
  ctr -= length;

  debug3(printf("  returning ctr %d\n",ctr));
  return ctr;
}



#if 0
void
Indexdb_convert_gammas (char *gammaptrsfile, char *offsetscompfile, FILE *offsets_fp,
#ifdef PMAP
			int alphabet_size, int index1part_aa,
#else
			int index1part,
#endif
			int blocksize) {
  int gammaptrs_fd, offsetscomp_fd;
  Positionsptr_T *offsets = NULL, totalcounts;
  int j;
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




Positionsptr_T *
Indexdb_offsets_from_gammas (char *gammaptrsfile, char *offsetscompfile, int offsetscomp_basesize
#ifdef PMAP
			     , int alphabet_size, int index1part_aa
#else
			     , int index1part
#endif
			     ) {
  UINT4 *gammaptrs, *offsetscomp;
  int gammaptrs_fd, offsetscomp_fd;
  size_t gammaptrs_len, offsetscomp_len;
  Positionsptr_T *offsets = NULL;
  int j;
  Oligospace_T oligospace, oligoi, oligok;
  int blocksize;
  double seconds;

  UINT4 *ptr, cum;
  int ctr;


#ifdef PMAP
  oligospace = power(alphabet_size,index1part_aa);
  blocksize = (int) power(alphabet_size,index1part_aa - offsetscomp_basesize);
#else
  oligospace = power(4,index1part);
  blocksize = (int) power(4,index1part - offsetscomp_basesize);
#endif

  if (blocksize == 1) {
    return (UINT4 *) Access_allocated(&offsetscomp_len,&seconds,offsetscompfile,sizeof(UINT4));

  } else {

#ifdef HAVE_MMAP
    gammaptrs = (UINT4 *) Access_mmap(&gammaptrs_fd,&gammaptrs_len,gammaptrsfile,sizeof(UINT4),/*randomp*/false);
    offsetscomp = (UINT4 *) Access_mmap(&offsetscomp_fd,&offsetscomp_len,offsetscompfile,sizeof(UINT4),/*randomp*/false);
#else
    gammaptrs = (UINT4 *) Access_allocated(&gammaptrs_len,&seconds,gammaptrsfile,sizeof(UINT4));
    offsetscomp = (UINT4 *) Access_allocated(&offsetscomp_len,&seconds,offsetscompfile,sizeof(UINT4));
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
      fprintf(stderr,"cannot allocated requested memory.  Cannot run -B 5 mode on this machine.\n");
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



static void
check_offsets_from_gammas (char *gammaptrsfile, char *offsetscompfile, Positionsptr_T *offsets,
			   Oligospace_T oligospace, int blocksize) {
  UINT4 *gammaptrs, *offsetscomp;
  int gammaptrs_fd, offsetscomp_fd;
  size_t gammaptrs_len, offsetscomp_len;
  Oligospace_T oligoi, oligok;
  int j, p;
#ifndef HAVE_MMAP
  double seconds;
#endif

  UINT4 *ptr, cum;
  int ctr;


#ifdef HAVE_MMAP
  gammaptrs = (UINT4 *) Access_mmap(&gammaptrs_fd,&gammaptrs_len,gammaptrsfile,sizeof(UINT4),/*randomp*/false);
  offsetscomp = (UINT4 *) Access_mmap(&offsetscomp_fd,&offsetscomp_len,offsetscompfile,sizeof(UINT4),/*randomp*/false);
#else
  gammaptrs = (UINT4 *) Access_allocated(&gammaptrs_len,&seconds,gammaptrsfile,sizeof(UINT4));
  offsetscomp = (UINT4 *) Access_allocated(&offsetscomp_len,&seconds,offsetscompfile,sizeof(UINT4));
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



#ifdef PMAP

/* PMAP version.  Doesn't mask bottom 12 nt. */
Genomicpos_T *
Indexdb_read (int *nentries, T this, unsigned int aaindex) {
  Genomicpos_T *positions;
  Positionsptr_T ptr, ptr0, end0;
  int i;
  char byte1, byte2, byte3;
  int bigendian, littleendian;

  debug0(printf("%u (%s)\n",aaindex,Alphabet_aaindex_aa(aaindex,this->alphabet)));

#ifdef WORDS_BIGENDIAN
  if (this->offsetscomp_access == ALLOCATED) {
    ptr0 = Genome_offsetptr_from_gammas(&end0,this->gammaptrs,this->offsetscomp,this->offsetscomp_blocksize,aaindex);
  } else {
    ptr0 = Genome_offsetptr_from_gammas_bigendian(&end0,this->gammaptrs,this->offsetscomp,this->offsetscomp_blocksize,aaindex);
  }
#else
  ptr0 = Genome_offsetptr_from_gammas(&end0,this->gammaptrs,this->offsetscomp,this->offsetscomp_blocksize,aaindex);
#endif

  debug0(printf("offset pointers are %u and %u\n",ptr0,end0));

#ifdef ALLOW_DUPLICATES
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
    positions = (Genomicpos_T *) CALLOC(*nentries,sizeof(Genomicpos_T));
    if (this->positions_access == FILEIO) {
#ifdef HAVE_PTHREAD
      pthread_mutex_lock(&this->positions_read_mutex);
#endif
      positions_move_absolute(this->positions_fd,ptr0);
      positions_read_multiple(this->positions_fd,positions,*nentries);
#ifdef HAVE_PTHREAD
      pthread_mutex_unlock(&this->positions_read_mutex);
#endif

    } else if (this->positions_access == ALLOCATED) {
      memcpy(positions,&(this->positions[ptr0]),(*nentries)*sizeof(Genomicpos_T));

    } else {
#ifdef WORDS_BIGENDIAN
      for (ptr = ptr0, i = 0; ptr < end0; ptr++, i++) {
	littleendian = this->positions[ptr];
	bigendian = littleendian & 0xff;
	byte1 = (littleendian >>= 8);
	byte2 = (littleendian >>= 8);
	byte3 = (littleendian >>= 8);
	  
	bigendian <<= 8;
	bigendian |= (byte1 & 0xff);
	bigendian <<= 8;
	bigendian |= (byte2 & 0xff);
	bigendian <<= 8;
	bigendian |= (byte3 & 0xff);
	  
	positions[i] = bigendian;
      }
#else
      memcpy(positions,&(this->positions[ptr0]),(*nentries)*sizeof(Genomicpos_T));
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
Genomicpos_T *
Indexdb_read (int *nentries, T this, Storedoligomer_T oligo) {
  Genomicpos_T *positions;
  Positionsptr_T ptr0, end0;
  Storedoligomer_T part0;
#ifdef WORDS_BIGENDIAN
  int i;
  Positionsptr_T ptr;
  int bigendian, littleendian;
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

#ifdef WORDS_BIGENDIAN
  if (this->offsetscomp_access == ALLOCATED) {
    ptr0 = Genome_offsetptr_from_gammas(&end0,this->gammaptrs,this->offsetscomp,this->offsetscomp_blocksize,part0);
  } else {
    ptr0 = Genome_offsetptr_from_gammas_bigendian(&end0,this->gammaptrs,this->offsetscomp,this->offsetscomp_blocksize,part0);
  }
#else
  ptr0 = Genome_offsetptr_from_gammas(&end0,this->gammaptrs,this->offsetscomp,this->offsetscomp_blocksize,part0);
#endif

#ifdef ALLOW_DUPLICATES
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
  
    positions = (Genomicpos_T *) CALLOC(*nentries,sizeof(Genomicpos_T));
    if (this->positions_access == FILEIO) {
#ifdef HAVE_PTHREAD
      pthread_mutex_lock(&this->positions_read_mutex);
#endif
      positions_move_absolute(this->positions_fd,ptr0);
      positions_read_multiple(this->positions_fd,positions,*nentries);
#ifdef HAVE_PTHREAD
      pthread_mutex_unlock(&this->positions_read_mutex);
#endif
    } else if (this->positions_access == ALLOCATED) {
      memcpy(positions,&(this->positions[ptr0]),(*nentries)*sizeof(Genomicpos_T));

    } else {
#ifdef WORDS_BIGENDIAN
      for (ptr = ptr0, i = 0; ptr < end0; ptr++, i++) {
	littleendian = this->positions[ptr];
	bigendian = littleendian & 0xff;
	byte1 = (littleendian >>= 8);
	byte2 = (littleendian >>= 8);
	byte3 = (littleendian >>= 8);
	
	bigendian <<= 8;
	bigendian |= (byte1 & 0xff);
	bigendian <<= 8;
	bigendian |= (byte2 & 0xff);
	bigendian <<= 8;
	bigendian |= (byte3 & 0xff);
	
	positions[i] = bigendian;
      }
#else
      memcpy(positions,&(this->positions[ptr0]),(*nentries)*sizeof(Genomicpos_T));
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
Genomicpos_T *
Indexdb_read_inplace (int *nentries, T this, Storedoligomer_T oligo) {
  Genomicpos_T *positions;
  Positionsptr_T ptr0, end0;
  Storedoligomer_T part0;
#ifdef DEBUG0
  Positionsptr_T ptr;
#endif

  debug0(printf("%08X (%s)\n",oligo,shortoligo_nt(oligo,index1part)));
  part0 = oligo & poly_T;

  /* Needed to avoid overflow on 15-mers */
  if (part0 == poly_A || part0 == poly_T) {
    *nentries = 0;
    return NULL;
  }

#ifdef WORDS_BIGENDIAN
  if (this->offsetscomp_access == ALLOCATED) {
    ptr0 = Genome_offsetptr_from_gammas(&end0,this->gammaptrs,this->offsetscomp,this->offsetscomp_blocksize,oligo);
  } else {
    ptr0 = Genome_offsetptr_from_gammas_bigendian(&end0,this->gammaptrs,this->offsetscomp,this->offsetscomp_blocksize,oligo);
  }
#else
  ptr0 = Genome_offsetptr_from_gammas(&end0,this->gammaptrs,this->offsetscomp,this->offsetscomp_blocksize,oligo);
#endif

  debug0(printf("Indexdb_read_inplace: offset pointers are %u and %u\n",ptr0,end0));

  *nentries = end0 - ptr0;

  if (*nentries == 0) {
    return NULL;
  } else if (this->positions_access == FILEIO) {
    positions = (Genomicpos_T *) CALLOC(*nentries,sizeof(Genomicpos_T));
#ifdef HAVE_PTHREAD
    pthread_mutex_lock(&this->positions_read_mutex);
#endif
    positions_move_absolute(this->positions_fd,ptr0);
    positions_read_multiple(this->positions_fd,positions,*nentries);
#ifdef HAVE_PTHREAD
    pthread_mutex_unlock(&this->positions_read_mutex);
#endif
    return positions;
  } else {
    debug0(
	   printf("%d entries:",*nentries);
	   for (ptr = ptr0; ptr < end0; ptr++) {
	     printf(" %u",this->positions[ptr]);
	   }
	   printf("\n");
	   );

    return &(this->positions[ptr0]);
  }
}

#endif	/* ifdef PMAP */


/* Analogous to Indexdb_read, except this includes diagterm.  Always allocates memory. */
Genomicpos_T *
Indexdb_read_with_diagterm (int *nentries, T this, Storedoligomer_T oligo, int diagterm) {
  Genomicpos_T *positions;
  Positionsptr_T ptr0, end0, ptr;
  int i;

#ifdef WORDS_BIGENDIAN
  if (this->offsetscomp_access == ALLOCATED) {
    ptr0 = Genome_offsetptr_from_gammas(&end0,this->gammaptrs,this->offsetscomp,this->offsetscomp_blocksize,oligo);
  } else {
    ptr0 = Genome_offsetptr_from_gammas_bigendian(&end0,this->gammaptrs,this->offsetscomp,this->offsetscomp_blocksize,oligo);
  }
#else
  ptr0 = Genome_offsetptr_from_gammas(&end0,this->gammaptrs,this->offsetscomp,this->offsetscomp_blocksize,oligo);
#endif

  debug0(printf("read_zero_shift: oligo = %06X, offset pointers are %u and %u\n",oligo,ptr0,end0));

  if ((*nentries = end0 - ptr0) == 0) {
    return (Genomicpos_T *) NULL;
  } else {
    positions = (Genomicpos_T *) CALLOC(*nentries,sizeof(Genomicpos_T));
    if (this->positions_access == FILEIO) {
#ifdef HAVE_PTHREAD
      pthread_mutex_lock(&this->positions_read_mutex);
#endif
      positions_move_absolute(this->positions_fd,ptr0);
      positions_read_multiple(this->positions_fd,positions,*nentries);

#ifdef HAVE_PTHREAD
      pthread_mutex_unlock(&this->positions_read_mutex);
#endif

    } else if (this->positions_access == ALLOCATED) {
      for (ptr = ptr0, i = 0; ptr < end0; ptr++) {
	positions[i++] = this->positions[ptr] + diagterm;
      }

    } else {

#ifdef WORDS_BIGENDIAN
      for (ptr = ptr0, i = 0; ptr < end0; ptr++) {
	positions[i++] = Bigendian_convert_uint(this->positions[ptr]) + diagterm;
      }
#else
      for (ptr = ptr0, i = 0; ptr < end0; ptr++) {
	positions[i++] = this->positions[ptr] + diagterm;
      }
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
Genomicpos_T *
Indexdb_read_with_diagterm_sizelimit (int *nentries, T this, Storedoligomer_T oligo, int diagterm,
				      int size_threshold) {
  Genomicpos_T *positions;
  Positionsptr_T ptr0, end0, ptr;
  int i;

#ifdef WORDS_BIGENDIAN
  if (this->offsetscomp_access == ALLOCATED) {
    ptr0 = Genome_offsetptr_from_gammas(&end0,this->gammaptrs,this->offsetscomp,this->offsetscomp_blocksize,oligo);
  } else {
    ptr0 = Genome_offsetptr_from_gammas_bigendian(&end0,this->gammaptrs,this->offsetscomp,this->offsetscomp_blocksize,oligo);
  }
#else
  ptr0 = Genome_offsetptr_from_gammas(&end0,this->gammaptrs,this->offsetscomp,this->offsetscomp_blocksize,oligo);
#endif

  debug0(printf("read_zero_shift: oligo = %06X, offset pointers are %u and %u\n",oligo,ptr0,end0));

  if ((*nentries = end0 - ptr0) == 0) {
    return (Genomicpos_T *) NULL;

  } else if (*nentries > size_threshold) {
    *nentries = 0;
    return (Genomicpos_T *) NULL;

  } else {
    positions = (Genomicpos_T *) CALLOC(*nentries,sizeof(Genomicpos_T));
    if (this->positions_access == FILEIO) {
#ifdef HAVE_PTHREAD
      pthread_mutex_lock(&this->positions_read_mutex);
#endif
      positions_move_absolute(this->positions_fd,ptr0);
      positions_read_multiple(this->positions_fd,positions,*nentries);

#ifdef HAVE_PTHREAD
      pthread_mutex_unlock(&this->positions_read_mutex);
#endif

    } else if (this->positions_access == ALLOCATED) {
      for (ptr = ptr0, i = 0; ptr < end0; ptr++) {
	positions[i++] = this->positions[ptr] + diagterm;
      }

    } else {

#ifdef WORDS_BIGENDIAN
      for (ptr = ptr0, i = 0; ptr < end0; ptr++) {
	positions[i++] = Bigendian_convert_uint(this->positions[ptr]) + diagterm;
      }
#else
      for (ptr = ptr0, i = 0; ptr < end0; ptr++) {
	positions[i++] = this->positions[ptr] + diagterm;
      }
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
 *   Write procedures -- called by gmapindex/pmapindex
 ************************************************************************/


/* Another MONITOR_INTERVAL is in compress.c */
#define MONITOR_INTERVAL 10000000 /* 10 million nt */
#define OFFSETS_BUFFER_SIZE 1000000

void
Indexdb_write_gammaptrs (char *gammaptrsfile, char *offsetsfile, Positionsptr_T *offsets,
			 Oligospace_T oligospace, int blocksize) {
  FILE *offsets_fp;
  FILE *gammaptrs_fp, *offsetscomp_fp;
  UINT4 *gammaptrs;
  int gammaptri;
  int j;
  Oligospace_T oligoi;

  UINT4 *offsets_buffer;
  int offsets_buffer_size = OFFSETS_BUFFER_SIZE;
  int offsets_buffer_i;

  UINT4 buffer;			/* gammabuffer */
  int ctr;
  UINT4 nwritten;


  if (blocksize == 1) {
    /* Don't write gammaptrs.  Write offsetscomp in a single command. */
    if ((offsets_fp = FOPEN_WRITE_BINARY(offsetsfile)) == NULL) {
      fprintf(stderr,"Can't write to file %s\n",offsetsfile);
      exit(9);
    } else {
      FWRITE_UINTS(offsets,oligospace+1,offsets_fp);
    }

  } else {
    gammaptrs = (UINT4 *) CALLOC(oligospace/blocksize+1,sizeof(UINT4));
    gammaptri = 0;

    if ((offsetscomp_fp = FOPEN_WRITE_BINARY(offsetsfile)) == NULL) {
      fprintf(stderr,"Can't write to file %s\n",offsetsfile);
      exit(9);
    }
    offsets_buffer = (UINT4 *) CALLOC(offsets_buffer_size,sizeof(UINT4));
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
	debug3(printf("writing gamma %08X\n",buffer));
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

static int frame
process_position (int *between_counter, int *in_counter,
		  Storedoligomer_T *high, Storedoligomer_T *low, Storedoligomer_T *carry,
		  Positionsptr_T *offsets, Genomicpos_T *positions,
		  int frame, int c, Genomicpos_T position, Genomicpos_T chrpos,
		  int index1part_aa, int index1interval, Storedoligomer_T mask, bool genome_lc_p,
		  bool watsonp) {
  Oligospace_T oligoi;
  Storedoligomer_T masked;
  unsigned int aaindex;

  if (++frame == 3) {
    frame = 0;
  }
  between_counter[frame] += 1;
  in_counter[frame] += 1;

  *carry = ((*low) >> 30);
  switch (uppercaseCode[c]) {
  case 'A': *low = ((*low) << 2); break;
  case 'C': *low = ((*low) << 2) | 1U; break;
  case 'G': *low = ((*low) << 2) | 2U; break;
  case 'T': *low = ((*low) << 2) | 3U; break;
  case 'X': case 'N': 
    *high = *low = *carry = 0U; 
    in_counter[0] = in_counter[1] = in_counter[2] = 0;
    break;
  default: 
    if (genome_lc_p == true) {
      *high = *low = *carry = 0U;
      in_counter[0] = in_counter[1] = in_counter[2] = 0;
    } else {
      fprintf(stderr,"Bad character %c at position %u\n",c,position);
      abort();
    }
  }
  *high = ((*high) << 2) | (*carry); 

  debug(printf("frame=%d char=%c bc=%d ic=%d high=%08X low=%08X\n",
	       frame,c,between_counter[frame],in_counter[frame],*high,*low));

  if (in_counter[frame] > 0) {
    if (watsonp == true) {
      if (Alphabet_get_codon_fwd(*low) == AA_STOP) {
	debug(printf("Resetting in_counter for frame %d to 0\n",frame));
	in_counter[frame] = 0; 
      }
    } else {
      if (Alphabet_get_codon_rev(*low) == AA_STOP) {
	debug(printf("Resetting in_counter for frame %d to 0\n",frame));
	in_counter[frame] = 0; 
      }
    }
  }
  if (in_counter[frame] == index1part_aa + 1) {
    if (between_counter[frame] >= index1interval) {
      aaindex = Alphabet_get_aa_index(*high,*low,watsonp,index1part_nt);
      if (positions == NULL) {
#ifdef OLIGOSPACE_NOT_LONG
	oligoi = (Oligospace_T) aaindex + 1U;
#else
	oligoi = (Oligospace_T) aaindex + 1UL;
#endif
	offsets[oligoi] += 1;
      } else if (watsonp == true) {
	positions[offsets[aaindex]++] = position-index1part_nt+1U;
	debug1(adjposition = position-index1part_nt+1U);
      } else {
	positions[offsets[aaindex]++] = position;
	debug1(adjposition = position);
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

  return frame;
}

#endif

#if 0

static Oligospace_T
process_position (int *between_counter, int *in_counter, Positionsptr_T *offsets, Genomicpos_T *positions,
		  Oligospace_T oligo, int c, Genomicpos_T position, Genomicpos_T chrpos,
		  int index1part, int index1interval, Storedoligomer_T mask, bool genome_lc_p) {
  Oligospace_T oligoi;
  Storedoligomer_T masked;

  (*between_counter)++;
  (*in_counter)++;

  switch (uppercaseCode[c]) {
  case 'A': oligo = (oligo << 2); break;
  case 'C': oligo = (oligo << 2) | 1U; break;
  case 'G': oligo = (oligo << 2) | 2U; break;
  case 'T': oligo = (oligo << 2) | 3U; break;
  case 'X': case 'N': oligo = 0U; *in_counter = 0; break;
  default: 
    if (genome_lc_p == true) {
      oligo = 0U; *in_counter = 0;
    } else {
      fprintf(stderr,"Bad character %c at position %u\n",c,position);
      abort();
    }
  }

  if (*in_counter == index1part) {
    if ((chrpos-index1part+1U) % index1interval == 0) { /* previous was between_counter >= index1interval */
      masked = oligo & mask;
      if (positions == NULL) {
#ifdef OLIGOSPACE_NOT_LONG
	oligoi = (Oligospace_T) masked + 1U;
#else
	oligoi = (Oligospace_T) masked + 1UL;
#endif
	offsets[oligoi] += 1;
	debug(printf("Found oligo %06X.  Incremented offsets for %lu to be %u\n",
		     masked,oligoi,offsets[oligoi]));
      } else {
	positions[offsets[masked]++] = position-index1part+1U;
	debug1(nt = shortoligo_nt(masked,index1part);
	       printf("Storing %s at %u, chrpos %u\n",nt,position-index1part+1U,chrpos-index1part+1U);
	       FREE(nt));
      }
      *between_counter = 0;
    }
    (*in_counter)--;
  }

  return oligo;
}

#endif



void
Indexdb_write_offsets (char *gammaptrsfile, char *offsetscompfile, FILE *sequence_fp, IIT_T chromosome_iit,
		       int offsetscomp_basesize,
#ifdef PMAP
		       int alphabet_size, int index1part_aa, bool watsonp,
#else
		       int index1part,
#endif
		       int index1interval, bool genome_lc_p, char *fileroot, bool mask_lowercase_p) {
  char *uppercaseCode;
  Positionsptr_T *offsets;
  char *comma;
  int c, nchrs, chrnum;
  Oligospace_T oligospace, oligoi;
  Genomicpos_T position = 0U, chrpos = 0U, next_chrbound;
#ifdef PMAP
  int frame = -1, between_counter[3], in_counter[3];
  Storedoligomer_T high = 0U, low = 0U, carry;
  unsigned int aaindex;
  int index1part_nt = 3*index1part_aa;
#else
  int between_counter = 0, in_counter = 0;
  Storedoligomer_T oligo = 0U, masked, mask;
#endif
#ifdef DEBUG1
  char *aa;
#endif

  int offsetscomp_blocksize;
  int circular_typeint;

  if (mask_lowercase_p == false) {
    uppercaseCode = UPPERCASE_U2T; /* We are reading DNA sequence */
  } else {
    uppercaseCode = NO_UPPERCASE;
  }

#ifdef PMAP
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
  offsetscomp_blocksize = power(alphabet_size,index1part_aa - offsetscomp_basesize);
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
  offsetscomp_blocksize = power(4,index1part - offsetscomp_basesize);
#endif

  /* Handle reference strain */
  circular_typeint = IIT_typeint(chromosome_iit,"circular");
  chrnum = 1;
  nchrs = IIT_total_nintervals(chromosome_iit);
  next_chrbound = IIT_next_chrbound(chromosome_iit,chrnum,circular_typeint);

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
	fprintf(stderr,"Bad character %c at position %u\n",c,position);
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
      while (chrnum <= nchrs && (next_chrbound = IIT_next_chrbound(chromosome_iit,chrnum,circular_typeint)) < position) {
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

#ifdef PRE_GAMMAS
  fprintf(stderr,"Writing %lu offsets to file with total of %u positions...",oligospace+1,offsets[oligospace]);
  FWRITE_UINTS(offsets,oligospace+1,offsets_fp);
#else
#ifdef OLIGOSPACE_NOT_LONG
  fprintf(stderr,"Writing %u offsets compressed to file with total of %u positions...",oligospace+1U,offsets[oligospace]);
#else
  fprintf(stderr,"Writing %lu offsets compressed to file with total of %u positions...",oligospace+1UL,offsets[oligospace]);
#endif
  Indexdb_write_gammaptrs(gammaptrsfile,offsetscompfile,offsets,oligospace,offsetscomp_blocksize);
#endif
  fprintf(stderr,"done\n");

  
  if (offsetscomp_blocksize > 1) {
    fprintf(stderr,"Checking gammas...");
    check_offsets_from_gammas(gammaptrsfile,offsetscompfile,offsets,oligospace,offsetscomp_blocksize);
    fprintf(stderr,"done\n");
  }

  FREE(offsets);

  return;
}


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
need_to_sort_p (Genomicpos_T *positions, int length) {
  Genomicpos_T prevpos;
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


/* Works directly in file, so we don't need to allocate memory */
static void
compute_positions_in_file (int positions_fd, Positionsptr_T *offsets,
			   FILE *sequence_fp, IIT_T chromosome_iit,
#ifdef PMAP
			   int index1part_aa, bool watsonp,
#else
			   int index1part,
#endif
			   int index1interval, bool genome_lc_p, char *fileroot, bool mask_lowercase_p) {
  char *uppercaseCode;
  Genomicpos_T position = 0U, chrpos = 0U, next_chrbound;
  char *comma;
  int c, nchrs, chrnum;
#ifdef PMAP
  int frame = -1, between_counter[3], in_counter[3];
  Storedoligomer_T high = 0U, low = 0U, carry;
  Oligospace_T aaindex;
  Genomicpos_T adjposition;
  int index1part_nt = 3*index1part_aa;
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
  circular_typeint = IIT_typeint(chromosome_iit,"circular");
  chrnum = 1;
  nchrs = IIT_total_nintervals(chromosome_iit);
  next_chrbound = IIT_next_chrbound(chromosome_iit,chrnum,circular_typeint);

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
	fprintf(stderr,"Bad character %c at position %u\n",c,position);
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
	positions_move_absolute(positions_fd,offsets[aaindex]);
	offsets[aaindex] += 1;
	if (watsonp == true) {
	  adjposition = position-index1part_nt+1U;
	} else {
	  adjposition = position;
	}
	WRITE_UINT(adjposition,positions_fd);
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
	positions_move_absolute(positions_fd,offsets[masked]);
	offsets[masked] += 1;
	WRITE_UINT(position-index1part+1U,positions_fd);
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
      while (chrnum <= nchrs && (next_chrbound = IIT_next_chrbound(chromosome_iit,chrnum,circular_typeint)) < position) {
	chrnum++;
      }
    }
    position++;
  }

#ifdef ADDSENTINEL
  for (oligoi = 0; oligoi < oligospace; oligoi++) {
    positions_move_absolute(positions_fd,offsets[oligoi]);
    WRITE_UINT(-1U,positions_fd);
  }
#endif

  return;
}

/* Requires sufficient memory to hold all positions */
static void
compute_positions_in_memory (Genomicpos_T *positions, Positionsptr_T *offsets,
			     FILE *sequence_fp, IIT_T chromosome_iit,
#ifdef PMAP
			     int index1part_aa, bool watsonp,
#else
			     int index1part,
#endif
			     int index1interval, bool genome_lc_p, char *fileroot, bool mask_lowercase_p) {
  char *uppercaseCode;
  Genomicpos_T position = 0U, chrpos = 0U, next_chrbound;
  char *comma;
  int c, nchrs, chrnum;
#ifdef PMAP
  int frame = -1, between_counter[3], in_counter[3];
  Storedoligomer_T high = 0U, low = 0U, carry;
  unsigned int aaindex;
  int index1part_nt = 3*index1part_aa;
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
  circular_typeint = IIT_typeint(chromosome_iit,"circular");
  chrnum = 1;
  nchrs = IIT_total_nintervals(chromosome_iit);
  next_chrbound = IIT_next_chrbound(chromosome_iit,chrnum,circular_typeint);

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
	fprintf(stderr,"Bad character %c at position %u\n",c,position);
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
	if (watsonp == true) {
	  positions[offsets[aaindex]++] = position-index1part_nt+1U;
	  debug1(adjposition = position-index1part_nt+1U);
	} else {
	  positions[offsets[aaindex]++] = position;
	  debug1(adjposition = position);
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
	positions[offsets[masked]++] = position-index1part+1U;
	debug1(nt = shortoligo_nt(masked,index1part);
	       printf("Storing %s at %u, chrpos %u\n",nt,position-index1part+1U,chrpos-index1part+1U);
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
      while (chrnum <= nchrs && (next_chrbound = IIT_next_chrbound(chromosome_iit,chrnum,circular_typeint)) < position) {
	chrnum++;
      }
    }
    position++;
  }

#ifdef ADDSENTINEL
  for (oligoi = 0; oligoi < oligospace; oligoi++) {
    positions[offsets[oligoi]] = -1U;
  }
  /*
  fprintf(stderr,"Put a sentinel at %u\n",offsets[0]);
  fprintf(stderr,"Put a sentinel at %u\n",offsets[oligospace-2]);
  fprintf(stderr,"Put a sentinel at %u\n",offsets[oligospace-1]);
  */
#endif

  return;
}


void
Indexdb_write_positions (char *positionsfile, char *gammaptrsfile, char *offsetscompfile,
			 FILE *sequence_fp, IIT_T chromosome_iit, int offsetscomp_basesize,
#ifdef PMAP
			 int alphabet_size, int index1part_aa, bool watsonp,
#else
			 int index1part,
#endif
			 int index1interval, bool genome_lc_p, bool writefilep,
			 char *fileroot, bool mask_lowercase_p) {
  FILE *positions_fp;		/* For building positions in memory */
  int positions_fd;		/* For building positions in file */
  Positionsptr_T *offsets = NULL, totalcounts;
  Genomicpos_T *positions;
  Oligospace_T oligospace;


#ifdef PMAP
  offsets = Indexdb_offsets_from_gammas(gammaptrsfile,offsetscompfile,offsetscomp_basesize,
					alphabet_size,index1part_aa);
  oligospace = power(alphabet_size,index1part_aa);
#else
  offsets = Indexdb_offsets_from_gammas(gammaptrsfile,offsetscompfile,offsetscomp_basesize,index1part);
  oligospace = power(4,index1part);
#endif
  totalcounts = offsets[oligospace];
  if (totalcounts == 0) {
    fprintf(stderr,"Something is wrong with the offsets file.  Total counts is zero.\n");
    exit(9);
  }


  if (writefilep == true) {
    fprintf(stderr,"User requested build of positions in file\n");
    positions_fd = Access_fileio_rw(positionsfile);
#ifdef PMAP
    compute_positions_in_file(positions_fd,offsets,sequence_fp,chromosome_iit,
			      index1part_aa,watsonp,index1interval,genome_lc_p,fileroot,mask_lowercase_p);
#else
    compute_positions_in_file(positions_fd,offsets,sequence_fp,chromosome_iit,
			      index1part,index1interval,genome_lc_p,fileroot,mask_lowercase_p);
#endif
    close(positions_fd);

  } else {
    fprintf(stderr,"Trying to allocate %u*%d bytes of memory...",totalcounts,(int) sizeof(Genomicpos_T));
    positions = (Genomicpos_T *) CALLOC_NO_EXCEPTION(totalcounts,sizeof(Genomicpos_T));
    if (positions == NULL) {
      fprintf(stderr,"failed.  Building positions in file.\n");
      positions_fd = Access_fileio_rw(positionsfile);
#ifdef PMAP
      compute_positions_in_file(positions_fd,offsets,sequence_fp,chromosome_iit,
				index1part_aa,watsonp,index1interval,genome_lc_p,fileroot,mask_lowercase_p);
#else
      compute_positions_in_file(positions_fd,offsets,sequence_fp,chromosome_iit,
				index1part,index1interval,genome_lc_p,fileroot,mask_lowercase_p);
#endif
      close(positions_fd);

    } else {
      fprintf(stderr,"succeeded.  Building positions in memory.\n");
      if ((positions_fp = FOPEN_WRITE_BINARY(positionsfile)) == NULL) {
	fprintf(stderr,"Can't open file %s\n",positionsfile);
	exit(9);
      }
#ifdef PMAP
      compute_positions_in_memory(positions,offsets,sequence_fp,chromosome_iit,
				  index1part_aa,watsonp,index1interval,genome_lc_p,fileroot,mask_lowercase_p);
#else
      compute_positions_in_memory(positions,offsets,sequence_fp,chromosome_iit,
				  index1part,index1interval,genome_lc_p,fileroot,mask_lowercase_p);
#endif
      fprintf(stderr,"Writing %u genomic positions to file %s ...\n",
	      totalcounts,positionsfile);
      FWRITE_UINTS(positions,totalcounts,positions_fp);

      fclose(positions_fp);
      FREE(positions);
    }
  }

  FREE(offsets);

  return;
}



/************************************************************************
 *   Create procedure -- for user-provided genomic segment
 ************************************************************************/

T
Indexdb_new_segment (char *genomicseg,
#ifdef PMAP
		     int alphabet_size, int index1part_aa, bool watsonp,
#else
		     int index1part,
#endif
		     int index1interval) {
  T new = (T) MALLOC(sizeof(*new));
  char *uppercaseCode;
  Positionsptr_T *work_offsets;	/* Working set for use in calculating positions */
  int totalcounts = 0;
  int c;
  Oligospace_T oligospace, oligoi;
  char *p;
  Genomicpos_T position = 0U;
#ifdef PMAP
  int frame = -1, between_counter[3], in_counter[3];
  Storedoligomer_T high = 0U, low = 0U, carry;
  unsigned int aaindex;
  int index1part_nt = 3*index1part_aa;
#else
  int between_counter = 0, in_counter = 0;
  Storedoligomer_T oligo = 0U, masked, mask;
#endif

  uppercaseCode = UPPERCASE_U2T;

#ifdef PMAP
  oligospace = power(alphabet_size,index1part_aa);
#else
  mask = ~(~0UL << 2*index1part);
  oligospace = power(4,index1part);
  new->index1interval = 1;
#endif


  new->gammaptrs = (UINT4 *) CALLOC(oligospace+1,sizeof(UINT4));
  for (oligoi = 0; oligoi <= oligospace; oligoi++) {
    new->gammaptrs[oligoi] = oligoi;
  }

  new->offsetscomp = (Positionsptr_T *) CALLOC(oligospace+1,sizeof(Positionsptr_T));

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
    new->offsetscomp[oligoi] = new->offsetscomp[oligoi] + new->offsetscomp[oligoi-1] + 1U;
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
  high = low = 0UL;
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
  new->positions = (Genomicpos_T *) CALLOC(totalcounts,sizeof(Genomicpos_T));
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
	  new->positions[work_offsets[aaindex]++] = position-index1part_nt+1U;
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
	  (position-index1part+1U) % index1interval == 0
#endif
	  ) {
	masked = oligo & mask;
	new->positions[work_offsets[masked]++] = position-index1part+1U;
	between_counter = 0;
      }
      in_counter--;
    }
#endif

    position++;
  }

#ifdef ADDSENTINEL
  for (oligoi = 0; oligoi < oligospace; oligoi++) {
    new->positions[work_offsets[oligoi]] = -1U;
  }
#endif

  FREE(work_offsets);

  return new;
}


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

