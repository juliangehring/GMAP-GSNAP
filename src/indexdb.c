static char rcsid[] = "$Id: indexdb.c 49487 2011-10-10 20:58:00Z twu $";
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

  poly_T = ~(~0U << 2*index1part_in);
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

    if ((*old)->offsets != NULL) {
      if ((*old)->offsets_access == ALLOCATED) {
	/* Could be expansion or backward compatibility with pregamma indices */
	FREE((*old)->offsets);
#ifdef HAVE_MMAP
      } else if ((*old)->offsets_access == MMAPPED) {
	/* Backward compatibility with pregamma indices */
	munmap((void *) (*old)->offsets,(*old)->offsets_len);
	close((*old)->offsets_fd);
#endif
      }

    } else {
      if ((*old)->offsetscomp_access == ALLOCATED) {
	FREE((*old)->offsetscomp);
#ifdef HAVE_MMAP
      } else if ((*old)->offsetscomp_access == MMAPPED) {
	munmap((void *) (*old)->offsetscomp,(*old)->offsetscomp_len);
	close((*old)->offsetscomp_fd);
#endif
      }
      
      FREE((*old)->gammaptrs);	/* Always ALLOCATED */
    }

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

static int
power (int base, int exponent) {
  int result = 1, i;

  for (i = 0; i < exponent; i++) {
    result *= base;
  }
  return result;
}

double
Indexdb_mean_size (T this, Mode_T mode, int index1part) {
  int oligospace, n;

#ifdef PMAP
  /* index1part should be in aa */
  n = oligospace = power(NAMINOACIDS,index1part);
#else
  n = oligospace = power(4,index1part);
  if (mode == CMET_STRANDED || mode == CMET_NONSTRANDED || mode == ATOI_STRANDED || mode == ATOI_NONSTRANDED) {
    n = power(3,index1part);
  }
#endif

  if (this->offsets) {
#ifdef WORDS_BIGENDIAN
    if (this->offsets_access == ALLOCATED) {
      return (double) this->offsets[oligospace]/(double) n;
    } else {
      return (double) Bigendian_convert_uint(this->offsets[oligospace])/(double) n;
    }
#else
    return (double) this->offsets[oligospace]/(double) n;
#endif
  } else {
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
      if ((q = strstr(p,offsets_suffix)) != NULL) {
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




bool
Indexdb_get_filenames (char **gammaptrs_filename, char **offsetscomp_filename, char **positions_filename,
		       char **gammaptrs_basename_ptr, char **offsetscomp_basename_ptr, char **positions_basename_ptr,
		       char **gammaptrs_index1info_ptr, char **offsetscomp_index1info_ptr, char **positions_index1info_ptr,
		       int *basesize, int *index1part, int *index1interval, char *genomesubdir,
		       char *fileroot, char *idx_filesuffix, char *snps_root,
		       int required_index1part, int required_interval) {
  char *base_filename, *filename;
  char *pattern, interval_char, digit_string[2], *p, *q;
  int found_basesize, found_index1part, found_interval;
  int rootlength, patternlength;

  char tens0, ones0, tens, ones;
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
      if ((q = strstr(p,offsetscomp_suffix)) != NULL) {

	if (q - p == 5) {
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
	  fprintf(stderr,"Cannot parse part between %s and offsets in filename %s\n",idx_filesuffix,filename);
	  return false;
	}

	if (required_index1part != 0 && required_interval != 0) {
	  if (found_index1part == required_index1part && found_interval == required_interval) {
	    *basesize = found_basesize;
	    *index1part = found_index1part;
	    *index1interval = found_interval;
	    FREE(base_filename);
	    base_filename = (char *) CALLOC(strlen(filename)+1,sizeof(char));
	    strcpy(base_filename,filename);
	  }
	} else if (required_index1part != 0 && required_interval == 0) {
	  if (found_index1part == required_index1part && found_interval < *index1interval) {
	    *basesize = found_basesize;
	    *index1part = found_index1part;
	    *index1interval = found_interval;
	    FREE(base_filename);
	    base_filename = (char *) CALLOC(strlen(filename)+1,sizeof(char));
	    strcpy(base_filename,filename);
	  }
	} else if (required_index1part == 0 && required_interval != 0) {
	  if (found_index1part > *index1part && found_interval == required_interval) {
	    *basesize = found_basesize;
	    *index1part = found_index1part;
	    *index1interval = found_interval;
	    FREE(base_filename);
	    base_filename = (char *) CALLOC(strlen(filename)+1,sizeof(char));
	    strcpy(base_filename,filename);
	  }
	} else {
	  if (found_index1part > *index1part) {
	    *basesize = found_basesize;
	    *index1part = found_index1part;
	    *index1interval = found_interval;
	    FREE(base_filename);
	    base_filename = (char *) CALLOC(strlen(filename)+1,sizeof(char));
	    strcpy(base_filename,filename);
	  } else if (found_index1part == *index1part && found_interval < *index1interval) {
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

  if (closedir(dp) < 0) {
    fprintf(stderr,"Unable to close directory %s\n",genomesubdir);
  }

  /* Construct full filenames */
  if (base_filename == NULL) {
    fprintf(stderr,"Cannot find offsetscomp file containing %s and %s",
	    idx_filesuffix,offsetscomp_suffix);
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
Indexdb_new_genome (int *index1part, char *genomesubdir, char *fileroot, char *idx_filesuffix, char *snps_root,
		    int required_index1part, int required_interval, bool expand_offsets_p,
		    Access_mode_T offsetscomp_access, Access_mode_T positions_access) {
  T new = (T) MALLOC(sizeof(*new));
  char *gammaptrs_filename, *offsetscomp_filename, *positions_filename,
    *gammaptrs_basename_ptr, *offsetscomp_basename_ptr, *positions_basename_ptr,
    *gammaptrs_index1info_ptr, *offsetscomp_index1info_ptr, *positions_index1info_ptr;
  char *offsets_filename, *offsets_basename_ptr, *offsets_index1info_ptr;
  Access_mode_T offsets_access;

  char *comma;
  double seconds;
#ifdef HAVE_MMAP
  int npages;
#endif

  /* Read offsets file */
  if (Indexdb_get_filenames(&gammaptrs_filename,&offsetscomp_filename,&positions_filename,
			    &gammaptrs_basename_ptr,&offsetscomp_basename_ptr,&positions_basename_ptr,
			    &gammaptrs_index1info_ptr,&offsetscomp_index1info_ptr,&positions_index1info_ptr,
			    &new->offsetscomp_basesize,&new->index1part,&new->index1interval,
			    genomesubdir,fileroot,idx_filesuffix,snps_root,
			    required_index1part,required_interval) == true) {
    *index1part = new->index1part;
    new->offsetscomp_blocksize = power(4,(*index1part) - new->offsetscomp_basesize);

    if (new->index1part == new->offsetscomp_basesize || expand_offsets_p == true) {
      new->offsets = Indexdb_offsets_from_gammas(gammaptrs_filename,offsetscomp_filename,
						 new->offsetscomp_basesize,new->index1part);
      new->offsets_access = ALLOCATED;

      new->gammaptrs = (UINT4 *) NULL;
      new->offsetscomp = (UINT4 *) NULL;

    } else {
      new->offsets = (Positionsptr_T *) NULL;

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
    *index1part = new->index1part;

    new->gammaptrs = (UINT4 *) NULL;
    new->offsetscomp = (UINT4 *) NULL;

    /* Interpret offsetscomp_access to be offsets_access */
    offsets_access = offsetscomp_access;

    if (offsets_access == USE_ALLOCATE) {
      if (snps_root) {
	fprintf(stderr,"Allocating memory for %s (%s) offsets, kmer %d, interval %d...",
		idx_filesuffix,snps_root,new->index1part,new->index1interval);
      } else {
	fprintf(stderr,"Allocating memory for %s offsets, kmer %d, interval %d...",
		idx_filesuffix,new->index1part,new->index1interval);
      }
      new->offsets = (UINT4 *) Access_allocated(&new->offsets_len,&seconds,
						offsets_filename,sizeof(Positionsptr_T));
      if (new->offsets == NULL) {
	fprintf(stderr,"insufficient memory (need to use a lower batch mode (-B))\n");
	exit(9);
      } else {
	comma = Genomicpos_commafmt(new->offsets_len);
	fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma,seconds);
	FREE(comma);
	new->offsets_access = ALLOCATED;
      }

#ifdef HAVE_MMAP
    } else if (offsets_access == USE_MMAP_PRELOAD) {
      if (snps_root) {
	fprintf(stderr,"Pre-loading %s (%s) offsets, kmer %d, interval %d...",
		idx_filesuffix,snps_root,new->index1part,new->index1interval);
      } else {
	fprintf(stderr,"Pre-loading %s offsets, kmer %d, interval %d...",
		idx_filesuffix,new->index1part,new->index1interval);
      }
      new->offsets = (UINT4 *) Access_mmap_and_preload(&new->offsets_fd,&new->offsets_len,&npages,&seconds,
						       offsets_filename,sizeof(Positionsptr_T));
      if (new->offsets == NULL) {
	fprintf(stderr,"insufficient memory (will use disk file instead, but program may not run)\n");
#ifdef PMAP
	new->offsets_access = FILEIO;
#else
	exit(9);
#endif
      } else {
	comma = Genomicpos_commafmt(new->offsets_len);
	fprintf(stderr,"done (%s bytes, %d pages, %.2f sec)\n",comma,npages,seconds);
	FREE(comma);
	new->offsets_access = MMAPPED;
      }

    } else if (offsets_access == USE_MMAP_ONLY) {
      new->offsets = (UINT4 *) Access_mmap(&new->offsets_fd,&new->offsets_len,
					   offsets_filename,sizeof(Positionsptr_T),/*randomp*/false);
      if (new->offsets == NULL) {
	fprintf(stderr,"Insufficient memory for mmap of %s (will use disk file instead, but program may not run)\n",
		offsets_filename);
#ifdef PMAP
	new->offsets_access = FILEIO;
#else
	exit(9);
#endif
      } else {
	new->offsets_access = MMAPPED;
      }
#endif

    } else if (offsets_access == USE_FILEIO) {
      fprintf(stderr,"Offsets file I/O access of %s not allowed\n",offsets_filename);
      exit(9);

    } else {
      fprintf(stderr,"Don't recognize offsets_access type %d\n",offsets_access);
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

/*                87654321 */
#define RIGHT_A 0x00000000
#define RIGHT_C 0x00000001
#define RIGHT_G 0x00000002
#define RIGHT_T 0x00000003

/*                      87654321 */
#define LOW_TWO_BITS  0x00000003

#ifdef PMAP

#if NAMINOACIDS == 20
static char aa_table[NAMINOACIDS] = "ACDEFGHIKLMNPQRSTVWY";
#elif NAMINOACIDS == 15
static char aa_table[NAMINOACIDS] = "ACDFGHIKMNPRSWY";
#elif NAMINOACIDS == 12
static char aa_table[NAMINOACIDS] = "ACDEFGHIKPSW";
#endif

#if (defined(DEBUG) || defined(DEBUG0) || defined(DEBUG1))
static char *
aaindex_aa (unsigned int aaindex) {
  char *aa;
  int i, j;
  int aasize = index1part_aa;

  aa = (char *) CALLOC(aasize+1,sizeof(char));
  j = aasize-1;
  for (i = 0; i < aasize; i++) {
    aa[j] = aa_table[aaindex % NAMINOACIDS];
    aaindex /= NAMINOACIDS;
    j--;
  }

  return aa;
}
#endif

/*               87654321 */
#define LEFT_A 0x00000000
#define LEFT_C 0x40000000
#define LEFT_G 0x80000000
#define LEFT_T 0xC0000000

#ifdef DEBUG
static char *
highlow_nt (Storedoligomer_T high, Storedoligomer_T low) {
  char *nt;
  int i, j;
  Storedoligomer_T lowbits;
  int oligosize = index1part_aa*3;

  nt = (char *) CALLOC(oligosize+1,sizeof(char));
  j = oligosize-1;
  for (i = 0; i < 16; i++) {
    lowbits = low & LOW_TWO_BITS;
    switch (lowbits) {
    case RIGHT_A: nt[j] = 'A'; break;
    case RIGHT_C: nt[j] = 'C'; break;
    case RIGHT_G: nt[j] = 'G'; break;
    case RIGHT_T: nt[j] = 'T'; break;
    }
    low >>= 2;
    j--;
  }
  for ( ; i < oligosize; i++) {
    lowbits = high & LOW_TWO_BITS;
    switch (lowbits) {
    case RIGHT_A: nt[j] = 'A'; break;
    case RIGHT_C: nt[j] = 'C'; break;
    case RIGHT_G: nt[j] = 'G'; break;
    case RIGHT_T: nt[j] = 'T'; break;
    }
    high >>= 2;
    j--;
  }

  return nt;
}
#endif


#else

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
write_gamma (int fd, unsigned int *nwritten, unsigned int *buffer, int ctr, unsigned int gamma) {
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
    WRITE_UINT(*buffer,fd);
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
			int index1part_aa,
#else
			int index1part,
#endif
			int blocksize) {
  int gammaptrs_fd, offsetscomp_fd;
  Positionsptr_T *offsets = NULL, totalcounts;
  int oligospace, i, j;

  UINT4 buffer;
  int ctr;
  UINT4 nwritten;

#ifdef PMAP
  oligospace = power(NAMINOACIDS,index1part_aa);
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
  for (i = 0; i < oligospace; i += blocksize) {
    WRITE_UINT(nwritten,gammaptrs_fd);

    WRITE_UINT(offsets[i],offsetscomp_fd);
    nwritten += 1;

    buffer = 0U;
    ctr = 32;
    for (j = 1; j < blocksize; j++) {
#ifdef ABSOLUTE_GAMMAS
      ctr = write_gamma(offsetscomp_fd,&nwritten,&buffer,ctr,offsets[i+j]-offsets[i]);
#else
      ctr = write_gamma(offsetscomp_fd,&nwritten,&buffer,ctr,offsets[i+j]-offsets[i+j-1]);
#endif
    }
    debug3(printf("writing gamma %08X\n",buffer));
    WRITE_UINT(buffer,offsetscomp_fd);
    nwritten += 1;
  }

  WRITE_UINT(offsets[i],offsetscomp_fd);
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
			     , int index1part_aa
#else
			     , int index1part
#endif
			     ) {
  UINT4 *gammaptrs, *offsetscomp;
  int gammaptrs_fd, offsetscomp_fd;
  size_t gammaptrs_len, offsetscomp_len;
  Positionsptr_T *offsets = NULL;
  int oligospace, i, j, k;
  int blocksize;
  double seconds;

  UINT4 *ptr, cum;
  int ctr;
#ifdef ABSOLUTE_GAMMAS
  UINT4 value;
#endif


#ifdef PMAP
  oligospace = power(NAMINOACIDS,index1part_aa);
  blocksize = power(4,index1part_aa - offsetscomp_basesize);
#else
  oligospace = power(4,index1part);
  blocksize = power(4,index1part - offsetscomp_basesize);
#endif

  if (blocksize == 1) {
    return (UINT4 *) Access_allocated(&offsetscomp_len,&seconds,offsetscompfile,sizeof(UINT4));

  } else {

#ifdef HAVE_MMAP
    gammaptrs = (UINT4 *) Access_mmap(&gammaptrs_fd,&gammaptrs_len,gammaptrsfile,sizeof(UINT4),/*randomp*/false);
    offsetscomp = (UINT4 *) Access_mmap(&offsetscomp_fd,&offsetscomp_len,offsetscompfile,sizeof(UINT4),/*randomp*/false);
#else
    gammaptrs = (UINT4 *) Access_allocated(&gammaptrs_len,&seconds,gammaptrsfile,sizeof(UINT4));
    offsetscomp = (UINT4 *) Access_allocated(&offsetscomp_len,&second,offsetscompfile,sizeof(UINT4));
#endif

#ifdef PMAP
    fprintf(stderr,"Allocating memory (%u words) for offsets, kmer %d...",oligospace+1,index1part_aa);
#else
    fprintf(stderr,"Allocating memory (%u words) for offsets, kmer %d...",oligospace+1,index1part);
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
    k = 0;


    for (i = 0; i < oligospace; i += blocksize) {
#ifdef HAVE_MMAP
#ifdef WORDS_BIGENDIAN
      cum = offsets[k++] = Bigendian_convert_uint(*ptr++);
#else
      cum = offsets[k++] = *ptr++;
#endif
#else
      cum = offsets[k++] = *ptr++;
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
	offsets[k++] = cum;
      }
      if (ctr > 0) {
	ptr++;			/* Done with last gamma byte */
      }
    }

#ifdef HAVE_MMAP
#ifdef WORDS_BIGENDIAN
    offsets[k++] = Bigendian_convert_uint(*ptr++);
#else
    offsets[k++] = *ptr++;
#endif
#else
    offsets[k++] = *ptr++;
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
			   int oligospace, int blocksize) {
  UINT4 *gammaptrs, *offsetscomp;
  int gammaptrs_fd, offsetscomp_fd;
  size_t gammaptrs_len, offsetscomp_len;
  int i, j, k, p;
#ifndef HAVE_MMAP
  double seconds;
#endif

  UINT4 *ptr, cum;
  int ctr;
#ifdef ABSOLUTE_GAMMAS
  UINT4 value;
#endif


#ifdef HAVE_MMAP
  gammaptrs = (UINT4 *) Access_mmap(&gammaptrs_fd,&gammaptrs_len,gammaptrsfile,sizeof(UINT4),/*randomp*/false);
  offsetscomp = (UINT4 *) Access_mmap(&offsetscomp_fd,&offsetscomp_len,offsetscompfile,sizeof(UINT4),/*randomp*/false);
#else
  gammaptrs = (UINT4 *) Access_allocated(&gammaptrs_len,&seconds,gammaptrsfile,sizeof(UINT4));
  offsetscomp = (UINT4 *) Access_allocated(&offsetscomp_len,&second,offsetscompfile,sizeof(UINT4));
#endif

  ptr = offsetscomp;
  k = 0;
  p = 0;

  for (i = 0; i < oligospace; i += blocksize) {
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
      fprintf(stderr,"Problem with gammaptrs at oligo %d: %u != %u.  Please inform twu@gene.com\n",
	      k,offsetscomp[gammaptrs[p-1]],cum);
      exit(9);
    }

    if (offsets[k++] != cum) {
      fprintf(stderr,"Problem with offsetscomp at oligo %d: %u != %u.  Please inform twu@gene.com\n",k-1,offsets[k-1],cum);
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

      if (offsets[k++] != cum) {
	fprintf(stderr,"Problem with offsetscomp at oligo %d: %u != %u.  Please inform twu@gene.com\n",k-1,offsets[k-1],cum);
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
    fprintf(stderr,"Problem with gammaptrs at oligo %d: %u != %u.  Please inform twu@gene.com\n",
	    k,offsetscomp[gammaptrs[p-1]],cum);
    exit(9);
  }
    
  if (offsets[k] != cum) {
    fprintf(stderr,"Problem with offsetscomp at oligo %d: %u != %u.  Please inform twu@gene.com\n",k,offsets[k],*ptr);
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

  debug0(printf("%u (%s)\n",aaindex,aaindex_aa(aaindex)));

  if (this->offsets) {
#ifdef WORDS_BIGENDIAN
    if (this->offsets_access == ALLOCATED) {
      ptr0 = this->offsets[aaindex];
      end0 = this->offsets[aaindex+1];
    } else {
      ptr0 = Bigendian_convert_uint(this->offsets[aaindex]);
      end0 = Bigendian_convert_uint(this->offsets[aaindex+1]);
    }
#else
    ptr0 = this->offsets[aaindex];
    end0 = this->offsets[aaindex+1];
#endif

  } else {
#ifdef WORDS_BIGENDIAN
    if (this->offsetscomp_access == ALLOCATED) {
      ptr0 = Genome_offsetptr_from_gammas(&end0,this->gammaptrs,this->offsetscomp,this->offsetscomp_blocksize,aaindex);
    } else {
      ptr0 = Genome_offsetptr_from_gammas_bigendian(&end0,this->gammaptrs,this->offsetscomp,this->offsetscomp_blocksize,aaindex);
    }
#else
    ptr0 = Genome_offsetptr_from_gammas(&end0,this->gammaptrs,this->offsetscomp,this->offsetscomp_blocksize,aaindex);
#endif
  }

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

  if (this->offsets) {
#ifdef WORDS_BIGENDIAN
    if (this->offsets_access == ALLOCATED) {
      ptr0 = this->offsets[part0];
      end0 = this->offsets[part0+1];
    } else {
      ptr0 = Bigendian_convert_uint(this->offsets[part0]);
      end0 = Bigendian_convert_uint(this->offsets[part0+1]);
    }
#else
    ptr0 = this->offsets[part0];
    end0 = this->offsets[part0+1];
#endif

  } else {
#ifdef WORDS_BIGENDIAN
    if (this->offsetscomp_access == ALLOCATED) {
      ptr0 = Genome_offsetptr_from_gammas(&end0,this->gammaptrs,this->offsetscomp,this->offsetscomp_blocksize,part0);
    } else {
      ptr0 = Genome_offsetptr_from_gammas_bigendian(&end0,this->gammaptrs,this->offsetscomp,this->offsetscomp_blocksize,part0);
    }
#else
    ptr0 = Genome_offsetptr_from_gammas(&end0,this->gammaptrs,this->offsetscomp,this->offsetscomp_blocksize,part0);
#endif
  }

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

  debug0(printf("%06X (%s)\n",oligo,shortoligo_nt(oligo,index1part)));
  part0 = oligo & poly_T;

  /* Needed to avoid overflow on 15-mers */
  if (part0 == poly_A || part0 == poly_T) {
    *nentries = 0;
    return NULL;
  }

  if (this->offsets) {
#ifdef WORDS_BIGENDIAN
    if (this->offsets_access == ALLOCATED) {
      ptr0 = this->offsets[oligo];
      end0 = this->offsets[oligo+1];
    } else {
      ptr0 = Bigendian_convert_uint(this->offsets[oligo]);
      end0 = Bigendian_convert_uint(this->offsets[oligo+1]);
    }
#else
    ptr0 = this->offsets[oligo];
    end0 = this->offsets[oligo+1];
#endif

  } else {
#ifdef WORDS_BIGENDIAN
    if (this->offsetscomp_access == ALLOCATED) {
      ptr0 = Genome_offsetptr_from_gammas(&end0,this->gammaptrs,this->offsetscomp,this->offsetscomp_blocksize,oligo);
    } else {
      ptr0 = Genome_offsetptr_from_gammas_bigendian(&end0,this->gammaptrs,this->offsetscomp,this->offsetscomp_blocksize,oligo);
    }
#else
    ptr0 = Genome_offsetptr_from_gammas(&end0,this->gammaptrs,this->offsetscomp,this->offsetscomp_blocksize,oligo);
#endif
  }

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

  if (this->offsets) {
#ifdef WORDS_BIGENDIAN
    if (this->offsets_access == ALLOCATED) {
      ptr0 = this->offsets[oligo];
      end0 = this->offsets[oligo+1];
    } else {
      ptr0 = Bigendian_convert_uint(this->offsets[oligo]);
      end0 = Bigendian_convert_uint(this->offsets[oligo+1]);
    }
#else
    ptr0 = this->offsets[oligo];
    end0 = this->offsets[oligo+1];
#endif

  } else {
#ifdef WORDS_BIGENDIAN
    if (this->offsetscomp_access == ALLOCATED) {
      ptr0 = Genome_offsetptr_from_gammas(&end0,this->gammaptrs,this->offsetscomp,this->offsetscomp_blocksize,oligo);
    } else {
      ptr0 = Genome_offsetptr_from_gammas_bigendian(&end0,this->gammaptrs,this->offsetscomp,this->offsetscomp_blocksize,oligo);
    }
#else
    ptr0 = Genome_offsetptr_from_gammas(&end0,this->gammaptrs,this->offsetscomp,this->offsetscomp_blocksize,oligo);
#endif
  }

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

  if (this->offsets) {
#ifdef WORDS_BIGENDIAN
    if (this->offsets_access == ALLOCATED) {
      ptr0 = this->offsets[oligo];
      end0 = this->offsets[oligo+1];
    } else {
      ptr0 = Bigendian_convert_uint(this->offsets[oligo]);
      end0 = Bigendian_convert_uint(this->offsets[oligo+1]);
    }
#else
    ptr0 = this->offsets[oligo];
    end0 = this->offsets[oligo+1];
#endif

  } else {
#ifdef WORDS_BIGENDIAN
    if (this->offsetscomp_access == ALLOCATED) {
      ptr0 = Genome_offsetptr_from_gammas(&end0,this->gammaptrs,this->offsetscomp,this->offsetscomp_blocksize,oligo);
    } else {
      ptr0 = Genome_offsetptr_from_gammas_bigendian(&end0,this->gammaptrs,this->offsetscomp,this->offsetscomp_blocksize,oligo);
    }
#else
    ptr0 = Genome_offsetptr_from_gammas(&end0,this->gammaptrs,this->offsetscomp,this->offsetscomp_blocksize,oligo);
#endif
  }

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

#ifdef PMAP
/*              87654321 */
#define CHAR1 0x00000030
#define CHAR2 0x0000000C
#define CHAR3 0x00000003

#define A1    0x00000000
#define C1    0x00000010
#define G1    0x00000020
#define T1    0x00000030

#define A2    0x00000000
#define C2    0x00000004
#define G2    0x00000008
#define T2    0x0000000C

#define A3    0x00000000
#define C3    0x00000001
#define G3    0x00000002
#define T3    0x00000003


#if NAMINOACIDS == 20
typedef enum {AA_A, AA_C, AA_D, AA_E, AA_F, 
	      AA_G, AA_H, AA_I, AA_K, AA_L,
	      AA_M, AA_N, AA_P, AA_Q, AA_R,
	      AA_S, AA_T, AA_V, AA_W, AA_Y,
	      AA_STOP} Aminoacid_T;
#elif NAMINOACIDS == 15
/* Amino acids are collapsed from 20 to 15, because 20^7*4 is greater
   than 2 GB, and 15^7*4 is less than 2 GB.  Now have DE, NQ, ST, and
   ILV. */
typedef enum {AA_A, AA_C, AA_DE, AA_F, 
	      AA_G, AA_H, AA_ILV, AA_K,
	      AA_M, AA_NQ, AA_P, AA_R,
	      AA_ST, AA_W, AA_Y, AA_STOP} Aminoacid_T;
#elif NAMINOACIDS == 12
/* 12-letter alphabet taken from Murphy, Wallqvist, and Levy.  Protein
   Engineering 13, 149-152, 2000. */
typedef enum {AA_A, AA_C, AA_DN, AA_EQ, AA_FY,
	      AA_G, AA_H, AA_ILMV, AA_KR, AA_P,
	      AA_ST, AA_W, AA_STOP} Aminoacid_T;
#else
#error The given value of NAMINOACIDS is not supported by indexdb.c
#endif

static int
get_codon_fwd (Storedoligomer_T codon) {
  switch (codon & CHAR2) {
  case T2:
    switch (codon & CHAR1) {
    case T1:
      switch (codon & CHAR3) {
      case T3: case C3:
#if NAMINOACIDS == 20
	return AA_F;
#elif NAMINOACIDS == 15
	return AA_F;
#elif NAMINOACIDS == 12
	return AA_FY;
#endif
      default: /* case A3: case G3: */
#if NAMINOACIDS == 20
	return AA_L;
#elif NAMINOACIDS == 15
	return AA_ILV;
#elif NAMINOACIDS == 12
	return AA_ILMV;
#endif
      } /* CHAR3 */

#if NAMINOACIDS == 20
    case A1:
      switch (codon & CHAR3) {
      case G3: return AA_M;
      default: /* case T3: case A3: case C3: */ return AA_I;
      }
    case C1: return AA_L;
    default: /* case G1: */ return AA_V;
#elif NAMINOACIDS == 15
    case A1:
      switch (codon & CHAR3) {
      case G3: return AA_M;
      default: /* case T3: case A3: case C3: */ return AA_ILV;
      }
    default: /* case C1: case G1: */ return AA_ILV;
#elif NAMINOACIDS == 12
    default: /* case A1: case C1: case G1: */ return AA_ILMV;
#endif
    } /* CHAR1 */

  case C2:
    switch (codon & CHAR1) {
    case C1: return AA_P;
#if NAMINOACIDS == 20
    case T1: return AA_S;
    case A1: return AA_T;
#elif NAMINOACIDS == 15
    case T1: case A1: return AA_ST;
#elif NAMINOACIDS == 12
    case T1: case A1: return AA_ST;
#endif
    default: /* case G1: */ return AA_A;
    } /* CHAR1 */
  case A2:
    switch (codon & CHAR1) {
    case T1:
      switch (codon & CHAR3) {
      case T3: case C3:
#if NAMINOACIDS == 20
	return AA_Y;
#elif NAMINOACIDS == 15
	return AA_Y;
#elif NAMINOACIDS == 12
	return AA_FY;
#endif
      default: /* case A3: case G3: */ return AA_STOP;
      }	/* CHAR3 */
    case C1:
      switch (codon & CHAR3) {
      case T3: case C3: return AA_H;
      default: /* case A3: case G3: */
#if NAMINOACIDS == 20
	return AA_Q;
#elif NAMINOACIDS == 15
	return AA_NQ;
#elif NAMINOACIDS == 12
	return AA_EQ;
#endif
      }	/* CHAR3 */
    case A1:
      switch (codon & CHAR3) {
      case T3: case C3: 
#if NAMINOACIDS == 20
	return AA_N;
#elif NAMINOACIDS == 15
	return AA_NQ;
#elif NAMINOACIDS == 12
	return AA_DN;
#endif
      default: /* case A3: case G3: */
#if NAMINOACIDS == 20
	return AA_K;
#elif NAMINOACIDS == 15
	return AA_K;
#elif NAMINOACIDS == 12
	return AA_KR;
#endif      
      }	/* CHAR3 */

    default: /* case G1: */
#if NAMINOACIDS == 20
      switch (codon & CHAR3) {
      case T3: case C3: return AA_D;
      default: /* case A3: case G3: */ return AA_E;
      }
#elif NAMINOACIDS == 15
      return AA_DE;
#elif NAMINOACIDS == 12
      switch (codon & CHAR3) {
      case T3: case C3: return AA_DN;
      default: /* case A3: case G3: */ return AA_EQ;
      }
#endif
    } /* CHAR1 */

  default: /* case G2: */
    switch (codon & CHAR1) {
    case T1:
      switch (codon & CHAR3) {
      case T3: case C3: return AA_C;
      case A3: return AA_STOP;
      default: /* case G3: */ return AA_W;
      }
    case C1:
#if NAMINOACIDS == 20
      return AA_R;
#elif NAMINOACIDS == 15
      return AA_R;
#elif NAMINOACIDS == 12
      return AA_KR;
#endif
    case A1:
      switch (codon & CHAR3) {
      case T3: case C3:
#if NAMINOACIDS == 20
	return AA_S;
#elif NAMINOACIDS == 15
	return AA_ST;
#elif NAMINOACIDS == 12
	return AA_ST;
#endif
      default: /* case A3: case G3: */
#if NAMINOACIDS == 20
	return AA_R;
#elif NAMINOACIDS == 15
	return AA_R;
#elif NAMINOACIDS == 12
	return AA_KR;
#endif
      }	/* CHAR3 */
    default: /* case G1: */ return AA_G;
    } /* CHAR1 */
  } /* CHAR2 */

  abort();
  return -1;
}


static int
get_codon_rev (Storedoligomer_T codon) {
  switch (codon & CHAR2) {
  case A2:
    switch (codon & CHAR3) {
    case A3:
      switch (codon & CHAR1) {
      case A1: case G1: 
#if NAMINOACIDS == 20
	return AA_F;
#elif NAMINOACIDS == 15
	return AA_F;
#elif NAMINOACIDS == 12
	return AA_FY;
#endif
      default: /* case T1: case C1: */
#if NAMINOACIDS == 20
	return AA_L;
#elif NAMINOACIDS == 15
	return AA_ILV;
#elif NAMINOACIDS == 12
	return AA_ILMV;
#endif
      } /* CHAR1 */

#if NAMINOACIDS == 20
    case T3: 
      switch (codon & CHAR1) {
      case C1: return AA_M;
      default: /* case A1: case T1: case G1: */ return AA_I;
      }
    case G3: return AA_L;
    default: /* case C3: */ return AA_V;
#elif NAMINOACIDS == 15
    case T3: 
      switch (codon & CHAR1) {
      case C1: return AA_M;
      default: /* case A1: case T1: case G1: */ return AA_ILV;
      }
    default: /* case G3: case C3: */ return AA_ILV;
#elif NAMINOACIDS == 12
    default: /* case G3: case T3: case C3: */ return AA_ILMV;
#endif
    }
  case G2:
    switch (codon & CHAR3) {
    case G3: return AA_P;
#if NAMINOACIDS == 20
    case A3: return AA_S;
    case T3: return AA_T;
#elif NAMINOACIDS == 15
    case A3: case T3: return AA_ST;
#else
    case A3: case T3: return AA_ST;
#endif
    default: /* case C3: */ return AA_A;
    }
  case T2:
    switch (codon & CHAR3) {
    case A3:
      switch (codon & CHAR1) {
      case A1: case G1: 
#if NAMINOACIDS == 20
	return AA_Y;
#elif NAMINOACIDS == 15
	return AA_Y;
#elif NAMINOACIDS == 12
	return AA_FY;
#endif
      default: /* case T1: case C1: */ return AA_STOP;
      }
    case G3:
      switch (codon & CHAR1) {
      case A1: case G1: return AA_H;
      default: /* case T1: case C1: */
#if NAMINOACIDS == 20
	return AA_Q;
#elif NAMINOACIDS == 15
	return AA_NQ;
#elif NAMINOACIDS == 12
	return AA_EQ;
#endif
      }	/* CHAR1 */
    case T3:
      switch (codon & CHAR1) {
      case A1: case G1: 
#if NAMINOACIDS == 20
	return AA_N;
#elif NAMINOACIDS == 15
	return AA_NQ;
#elif NAMINOACIDS == 12
	return AA_DN;
#endif
      default: /* case T1: case C1: */ 
#if NAMINOACIDS == 20
	return AA_K;
#elif NAMINOACIDS == 15
	return AA_K;
#elif NAMINOACIDS == 12
	return AA_KR;
#endif      
      }	/* CHAR1 */
    default: /* case C3: */
#if NAMINOACIDS == 20
      switch (codon & CHAR1) {
      case A1: case G1: return AA_D;
      default: /* case T1: case C1: */ return AA_E;
      }
#elif NAMINOACIDS == 15
      return AA_DE;
#elif NAMINOACIDS == 12
      switch (codon & CHAR1) {
      case A1: case G1: return AA_DN;
      default: /* case T1: case C1: */ return AA_EQ;
      }
#endif
    } /* CHAR3 */

  default: /* case C2: */
    switch (codon & CHAR3) {
    case A3:
      switch (codon & CHAR1) {
      case A1: case G1: return AA_C;
      case T1: return AA_STOP;
      default: /* case C1: */ return AA_W;
      }
    case G3: 
#if NAMINOACIDS == 20
      return AA_R;
#elif NAMINOACIDS == 15
      return AA_R;
#elif NAMINOACIDS == 12
      return AA_KR;
#endif
    case T3:
      switch (codon & CHAR1) {
      case A1: case G1:
#if NAMINOACIDS == 20
	return AA_S;
#elif NAMINOACIDS == 15
	return AA_ST;
#elif NAMINOACIDS == 12
	return AA_ST;
#endif
      default: /* case T1: case C1: */ 
#if NAMINOACIDS == 20
	return AA_R;
#elif NAMINOACIDS == 15
	return AA_R;
#elif NAMINOACIDS == 12
	return AA_KR;
#endif
      }
    default: /* case C3: */ return AA_G;
    }
  }

  abort();
  return -1;
}



static Storedoligomer_T
offset_codon (Storedoligomer_T high, Storedoligomer_T low, int offset) {
  Storedoligomer_T shifted;
  char *nt;

  if (offset < 14) {
    shifted = (low >> (offset + offset));
  } else if (offset == 14) {
    shifted = (high << 4) | (low >> 28);
  } else if (offset == 15) {
    shifted = (high << 2) | (low >> 30);
  } else {
    shifted = (high >> (offset + offset - 32));
  }

  debug2(
	nt = shortoligo_nt(shifted,3);
	printf("Offset = %d.  Shifted = %s.  %2X %2X %2X\n",
	       offset,nt,shifted & CHAR1,shifted & CHAR2,shifted & CHAR3);
	FREE(nt);
	);

  return shifted;
}


static unsigned int
get_aa_index (Storedoligomer_T high, Storedoligomer_T low, bool watsonp, int index1part_nt) {
  unsigned int aaindex = 0U;
  Storedoligomer_T shifted;
  int i, codonindex;
  char *nt, *aa;

  if (watsonp == true) {
    for (i = index1part_nt-3; i >= 0; i -= 3) {
      shifted = offset_codon(high,low,i);
      if ((codonindex = get_codon_fwd(shifted)) == AA_STOP) {
	fprintf(stderr,"Unexpected stop codon in get_aa_index\n");
	abort();
      } else {
	aaindex = aaindex * NAMINOACIDS + codonindex;
      }
    }
  } else {
    for (i = 0; i < index1part_nt; i += 3) {
      shifted = offset_codon(high,low,i);
      if ((codonindex = get_codon_rev(shifted)) == AA_STOP) {
	fprintf(stderr,"Unexpected stop codon in get_aa_index\n");
	abort();
      } else {
	aaindex = aaindex * NAMINOACIDS + codonindex;
      }
    }
  }

  debug(
	nt = highlow_nt(high,low);
	aa = aaindex_aa(aaindex);
	printf("For oligo %s, aa is %s\n",nt,aa);
	FREE(aa);
	FREE(nt);
	);

  return aaindex;
}
#endif


/* Another MONITOR_INTERVAL is in compress.c */
#define MONITOR_INTERVAL 10000000 /* 10 million nt */


void
Indexdb_write_gammaptrs (char *gammaptrsfile, char *offsetsfile, Positionsptr_T *offsets,
			 int oligospace, int blocksize) {
  int gammaptrs_fd, offsets_fd;
  int i, j;

  UINT4 buffer;
  int ctr;
  UINT4 nwritten;


  if (blocksize == 1) {
    /* Don't write gammaptrs */
    offsets_fd = Access_fileio_rw(offsetsfile);
    for (i = 0; i <= oligospace; i++) {
      WRITE_UINT(offsets[i],offsets_fd);
    }
    close(offsets_fd);

  } else {
    gammaptrs_fd = Access_fileio_rw(gammaptrsfile);
    offsets_fd = Access_fileio_rw(offsetsfile);

    nwritten = 0U;
    for (i = 0; i < oligospace; i += blocksize) {
      WRITE_UINT(nwritten,gammaptrs_fd);

      WRITE_UINT(offsets[i],offsets_fd);
      nwritten += 1;

      if (blocksize > 1) {
	buffer = 0U;
	ctr = 32;
	for (j = 1; j < blocksize; j++) {
#ifdef ABSOLUTE_GAMMAS
	  ctr = write_gamma(offsets_fd,&nwritten,&buffer,ctr,offsets[i+j]-offsets[i]);
#else
	  ctr = write_gamma(offsets_fd,&nwritten,&buffer,ctr,offsets[i+j]-offsets[i+j-1]);
#endif
	}
	debug3(printf("writing gamma %08X\n",buffer));
	WRITE_UINT(buffer,offsets_fd);
	nwritten += 1;
      }
    }


    /* Final entries for i == oligospace */
    WRITE_UINT(nwritten,gammaptrs_fd);
    WRITE_UINT(offsets[i],offsets_fd);
    nwritten += 1;

    close(offsets_fd);
    close(gammaptrs_fd);
  }

  return;
}


void
Indexdb_write_offsets (char *gammaptrsfile, char *offsetscompfile, FILE *sequence_fp, IIT_T chromosome_iit,
		       int offsetscomp_basesize,
#ifdef PMAP
		       int index1part_aa, bool watsonp,
#else
		       int index1part,
#endif
		       int index1interval, bool genome_lc_p, char *fileroot, bool mask_lowercase_p) {
  char *uppercaseCode;
  Positionsptr_T *offsets;
  char *comma;
  int c, nchrs, chrnum, oligospace, i;
  Genomicpos_T position = 0U, chrpos = 0U, next_chrbound;
#ifdef PMAP
  int frame = -1, between_counter[3], in_counter[3];
  Storedoligomer_T high = 0U, low = 0U, carry;
  unsigned int aaindex;
  char *aa;
  int index1part_nt = 3*index1part_aa;
#else
  int between_counter = 0, in_counter = 0;
  Storedoligomer_T oligo = 0U, masked, mask;
#endif

  int offsetscomp_blocksize;


  if (mask_lowercase_p == false) {
    uppercaseCode = UPPERCASE_U2T; /* We are reading DNA sequence */
  } else {
    uppercaseCode = NO_UPPERCASE;
  }

#ifdef PMAP
  oligospace = power(NAMINOACIDS,index1part_aa);
  between_counter[0] = between_counter[1] = between_counter[2] = 0;
  in_counter[0] = in_counter[1] = in_counter[2] = 0;
  offsets = (Positionsptr_T *) CALLOC_NO_EXCEPTION(oligospace+1,sizeof(Positionsptr_T));
  if (offsets == NULL) {
    fprintf(stderr,"Unable to allocate %d bytes of memory, needed to build offsets with %d-mers\n",oligospace+1,index1part_aa);
    fprintf(stderr,"Either find a computer with more RAM, or lower your value for the k-mer size\n");
    exit(9);
  }
  offsetscomp_blocksize = power(4,index1part_aa - offsetscomp_basesize);
#else
  mask = ~(~0UL << 2*index1part);
  oligospace = power(4,index1part);
  offsets = (Positionsptr_T *) CALLOC_NO_EXCEPTION(oligospace+1,sizeof(Positionsptr_T));
  if (offsets == NULL) {
    fprintf(stderr,"Unable to allocate %d bytes of memory, needed to build offsets with %d-mers\n",oligospace+1,index1part);
    fprintf(stderr,"Either find a computer with more RAM, or lower your value for the k-mer size\n");
    exit(9);
  }
  offsetscomp_blocksize = power(4,index1part - offsetscomp_basesize);
#endif

  /* Handle reference strain */
  chrnum = 1;
  nchrs = IIT_total_nintervals(chromosome_iit);
  next_chrbound = IIT_interval_high(chromosome_iit,chrnum);

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
	if (get_codon_fwd(low) == AA_STOP) {
	  debug(printf("Resetting in_counter for frame %d to 0\n",frame));
	  in_counter[frame] = 0; 
	}
      } else {
	if (get_codon_rev(low) == AA_STOP) {
	  debug(printf("Resetting in_counter for frame %d to 0\n",frame));
	  in_counter[frame] = 0; 
	}
      }
    }
    if (in_counter[frame] == index1part_aa + 1) {
      if (between_counter[frame] >= index1interval) {
	aaindex = get_aa_index(high,low,watsonp,index1part_nt);
	offsets[aaindex + 1U] += 1;
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
	offsets[masked + 1U] += 1;
	debug(printf("Found oligo %06X.  Incremented offsets for %d to be %d\n",
		     masked,masked+1U,offsets[masked+1U]));
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
      while (chrnum <= nchrs && (next_chrbound = IIT_interval_high(chromosome_iit,chrnum)) < position) {
	chrnum++;
      }
    }
    position++;
  }


#ifdef ADDSENTINEL
  for (i = 1; i <= oligospace; i++) {
    offsets[i] = offsets[i] + offsets[i-1] + 1U;
    debug(if (offsets[i] != offsets[i-1]) {
	    printf("Offset for %06X: %u\n",i,offsets[i]);
	  });
  }
#else
  for (i = 1; i <= oligospace; i++) {
    offsets[i] = offsets[i] + offsets[i-1];
    debug(if (offsets[i] != offsets[i-1]) {
	    printf("Offset for %06X: %u\n",i,offsets[i]);
	  });
  }
#endif

  /*
  fprintf(stderr,"Offset for A...A is %u to %u\n",offsets[0],offsets[1]);
  fprintf(stderr,"Offset for T...T is %u to %u\n",offsets[oligospace-1],offsets[oligospace]);
  */

  fprintf(stderr,"Writing %d offsets to file with total of %u positions...",oligospace+1,offsets[oligospace]);
#ifdef PRE_GAMMAS
  FWRITE_UINTS(offsets,oligospace+1,offsets_fp);
#else
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
  unsigned int aaindex;
  Genomicpos_T adjposition;
  int index1part_nt = 3*index1part_aa;
#else
  int between_counter = 0, in_counter = 0;
  Storedoligomer_T oligo = 0U, masked, mask;
#endif

#ifdef ADDSENTINEL
  int oligospace;
  oligospace = power(4,index1part);
#endif

  if (mask_lowercase_p == false) {
    uppercaseCode = UPPERCASE_U2T; /* We are reading DNA sequence */
  } else {
    uppercaseCode = NO_UPPERCASE;
  }

#ifdef PMAP
  /* oligospace = power(NAMINOACIDS,index1part_aa); */
#else
  mask = ~(~0UL << 2*index1part);
  /* oligospace = power(4,index1part); */
#endif

  /* Handle reference strain */
  chrnum = 1;
  nchrs = IIT_total_nintervals(chromosome_iit);
  next_chrbound = IIT_interval_high(chromosome_iit,chrnum);

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
	if (get_codon_fwd(low) == AA_STOP) {
	  in_counter[frame] = 0; 
	}
      } else {
	if (get_codon_rev(low) == AA_STOP) {
	  in_counter[frame] = 0; 
	}
      }
    }
    if (in_counter[frame] == index1part_aa + 1) {
      if (between_counter[frame] >= index1interval) {
	aaindex = get_aa_index(high,low,watsonp,index1part_nt);
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
      while (chrnum <= nchrs && (next_chrbound = IIT_interval_high(chromosome_iit,chrnum)) < position) {
	chrnum++;
      }
    }
    position++;
  }

#ifdef ADDSENTINEL
  for (i = 0; i < oligospace; i++) {
    positions_move_absolute(positions_fd,offsets[i]);
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

#ifdef ADDSENTINEL
  int oligospace;
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
  chrnum = 1;
  nchrs = IIT_total_nintervals(chromosome_iit);
  next_chrbound = IIT_interval_high(chromosome_iit,chrnum);

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
	if (get_codon_fwd(low) == AA_STOP) {
	  in_counter[frame] = 0; 
	}
      } else {
	if (get_codon_rev(low) == AA_STOP) {
	  in_counter[frame] = 0; 
	}
      }
    }
    if (in_counter[frame] == index1part_aa + 1) {
      if (between_counter[frame] >= index1interval) {
	aaindex = get_aa_index(high,low,watsonp,index1part_nt);
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
      while (chrnum <= nchrs && (next_chrbound = IIT_interval_high(chromosome_iit,chrnum)) < position) {
	chrnum++;
      }
    }
    position++;
  }

#ifdef ADDSENTINEL
  for (i = 0; i < oligospace; i++) {
    positions[offsets[i]] = -1U;
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
			 int index1part_aa, bool watsonp,
#else
			 int index1part,
#endif
			 int index1interval, bool genome_lc_p, bool writefilep,
			 char *fileroot, bool mask_lowercase_p) {
  FILE *positions_fp;		/* For building positions in memory */
  int positions_fd;		/* For building positions in file */
  Positionsptr_T *offsets = NULL, totalcounts;
  Genomicpos_T *positions;
  int oligospace;


#ifdef PMAP
  offsets = Indexdb_offsets_from_gammas(gammaptrsfile,offsetscompfile,offsetscomp_basesize,index1part_aa);
  oligospace = power(NAMINOACIDS,index1part_aa);
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
		     int index1part_aa, bool watsonp,
#else
		     int index1part,
#endif
		     int index1interval) {
  T new = (T) MALLOC(sizeof(*new));
  char *uppercaseCode;
  Positionsptr_T *work_offsets;	/* Working set for use in calculating positions */
  int totalcounts = 0;
  int c, oligospace, i;
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
  oligospace = power(NAMINOACIDS,index1part_aa);
#else
  mask = ~(~0UL << 2*index1part);
  oligospace = power(4,index1part);
  new->index1interval = 1;
#endif


  new->offsets = (Positionsptr_T *) CALLOC(oligospace+1,sizeof(Positionsptr_T));


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
	if (get_codon_fwd(low) == AA_STOP) {
	  in_counter[frame] = 0; 
	}
      } else {
	if (get_codon_rev(low) == AA_STOP) {
	  in_counter[frame] = 0; 
	}
      }
    }
    if (in_counter[frame] == index1part_aa + 1) {
      if (between_counter[frame] >= index1interval) {
	aaindex = get_aa_index(high,low,watsonp,index1part_nt);
	new->offsets[aaindex + 1U] += 1;
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
	new->offsets[masked + 1U] += 1;
	debug(printf("Found oligo %06X.  Incremented offsets for %d to be %d\n",
		     masked,masked+1U,new->offsets[masked+1U]));
	between_counter = 0;
      }
      in_counter--;
    }
#endif

    position++;
  }

#ifdef ADDSENTINEL
  for (i = 1; i <= oligospace; i++) {
    new->offsets[i] = new->offsets[i] + new->offsets[i-1] + 1U;
    debug(printf("Offset for %06X: %u\n",i,new->offsets[i]));
  }
#else
  for (i = 1; i <= oligospace; i++) {
    new->offsets[i] = new->offsets[i] + new->offsets[i-1];
    debug(printf("Offset for %06X: %u\n",i,new->offsets[i]));
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
  for (i = 0; i <= oligospace; i++) {
    work_offsets[i] = new->offsets[i];
  }

  totalcounts = new->offsets[oligospace];
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
	if (get_codon_fwd(low) == AA_STOP) {
	  in_counter[frame] = 0; 
	}
      } else {
	if (get_codon_rev(low) == AA_STOP) {
	  in_counter[frame] = 0; 
	}
      }
    }

    if (in_counter[frame] == index1part_aa + 1) {
      if (between_counter[frame] >= index1interval) {
	aaindex = get_aa_index(high,low,watsonp,index1part_nt);
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
  for (i = 0; i < oligospace; i++) {
    new->positions[work_offsets[i]] = -1U;
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

