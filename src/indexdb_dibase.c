static char rcsid[] = "$Id: indexdb.c,v 1.136 2009/08/14 15:27:12 twu Exp $";
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

#define EXTRA_FOR_DIBASE 1	/* To get 12 colors, need in_counter to go to 13 */
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

    if ((*old)->offsets_access == ALLOCATED) {
      FREE((*old)->offsets);
#ifdef HAVE_MMAP
    } else if ((*old)->offsets_access == MMAPPED) {
      munmap((void *) (*old)->offsets,(*old)->offsets_len);
      close((*old)->offsets_fd);
#endif
    } else if ((*old)->offsets_access == FILEIO) {
      close((*old)->offsets_fd);
    }

    FREE(*old);
  }
  return;
}


int
Indexdb_interval (T this) {
  return this->index1interval;
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
Indexdb_mean_size (T this, bool cmetp) {
  int oligospace, n;

  n = oligospace = power(4,INDEX1PART);
  if (cmetp == true) {
    n = power(3,INDEX1PART);
  }

  return (double) this->offsets[oligospace]/(double) n;
}



/*                 87654321 */
#define LOW12MER 0x00FFFFFF	/* Applies if we have index1interval = 1 */


#define LARGEVALUE 1000000

T
Indexdb_new_genome (char *genomesubdir, char *fileroot, char *idx_filesuffix, char *snps_root,
		    int required_interval, bool batch_offsets_p, bool batch_positions_p) {
  T new = (T) MALLOC(sizeof(*new));
  char *filename, *pattern, interval_char, interval_string[2], best_interval_char, *p;
  char *offsets_suffix, *positions_suffix;
  struct dirent *entry;
  DIR *dp;
  size_t len;
  double seconds;
  int interval;
#ifdef HAVE_MMAP
  int npages;
#endif

  /* Read offsets file */
  if (snps_root == NULL) {
    offsets_suffix = OFFSETS_FILESUFFIX;
  } else {
    offsets_suffix = (char *) CALLOC(strlen(OFFSETS_FILESUFFIX)+strlen(".")+strlen(snps_root)+1,sizeof(char));
    sprintf(offsets_suffix,"%s.%s",OFFSETS_FILESUFFIX,snps_root);
  }

  if (required_interval > 0) {
    new->index1interval = required_interval;
    sprintf(interval_string,"%d",required_interval);
    best_interval_char = interval_string[0];

  } else {
    new->index1interval = 0;
    best_interval_char = ' ';

    if ((dp = opendir(genomesubdir)) == NULL) {
      fprintf(stderr,"Unable to open directory %s\n",genomesubdir);
      exit(9);
    }

    pattern = (char *) CALLOC(strlen(fileroot)+strlen(".")+strlen(idx_filesuffix)+1,sizeof(char));
    sprintf(pattern,"%s.%s",fileroot,idx_filesuffix);

    interval_string[1] = '\0';
    while ((entry = readdir(dp)) != NULL) {
      filename = entry->d_name;
      if (!strncmp(filename,pattern,strlen(pattern))) {
	p = &(filename[strlen(pattern)]);
	p += 1;			/* Advance past character after "id" */
	if (!strcmp(p,"offsets")) {
	  p -= 1;
	  if (sscanf(p,"%c",&interval_char) == 1) {
	    if (interval_char == 'x') {
	      interval = 6;
	    } else {
	      interval_string[0] = interval_char;
	      interval = atoi(interval_string);
	    }

	    if (new->index1interval == 0 || interval < new->index1interval) {
	      new->index1interval = interval;
	      best_interval_char = interval_char;
	    }
	  }
	}
      }
    }
    FREE(pattern);

    if (closedir(dp) < 0) {
      fprintf(stderr,"Unable to close directory %s\n",genomesubdir);
    }
  }

  /* Offsets */
  filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			     strlen(fileroot)+strlen(".")+strlen(idx_filesuffix)+
			     /*for interval_char*/1+strlen(offsets_suffix)+1,sizeof(char));
  sprintf(filename,"%s/%s.%s%c%s",genomesubdir,fileroot,idx_filesuffix,best_interval_char,offsets_suffix);
  if (Access_file_exists_p(filename) == false) {
    FREE(filename);
    FREE(new);
    return NULL;
  }
  if (snps_root != NULL) {
    FREE(offsets_suffix);
  }

#ifndef HAVE_MMAP
  new->offsets = (Positionsptr_T *) Access_allocated(&len,&seconds,filename,sizeof(Positionsptr_T));
  new->offsets_access = ALLOCATED;
#else
  if (batch_offsets_p == true) {
    if (snps_root) {
      fprintf(stderr,"Pre-reading %s%c (%s) offsets db...",idx_filesuffix,best_interval_char,snps_root);
    } else {
      fprintf(stderr,"Pre-reading %s%c offsets db...",idx_filesuffix,best_interval_char);
    }
    new->offsets = (Positionsptr_T *) Access_allocated(&len,&seconds,filename,sizeof(Positionsptr_T));
    fprintf(stderr,"done (%u bytes, %.2f sec)\n",(unsigned int) len,seconds);
    new->offsets_access = ALLOCATED;
  } else {
    new->offsets = (Positionsptr_T *) Access_mmap(&new->offsets_fd,&new->offsets_len,
						  filename,sizeof(Positionsptr_T),/*randomp*/false);
    if (new->offsets == NULL) {
      new->offsets = (Positionsptr_T *) Access_allocated(&len,&seconds,filename,sizeof(Positionsptr_T));
      close(new->offsets_fd);	/* Needed here because ALLOCATED implies fd is closed */
      new->offsets_access = ALLOCATED;
    } else {
      new->offsets_access = MMAPPED;
    }
  }
#endif

  FREE(filename);


  /* Positions */
  if (snps_root == NULL) {
    positions_suffix = POSITIONS_FILESUFFIX;
  } else {
    positions_suffix = (char *) CALLOC(strlen(POSITIONS_FILESUFFIX)+strlen(".")+strlen(snps_root)+1,sizeof(char));
    sprintf(positions_suffix,"%s.%s",POSITIONS_FILESUFFIX,snps_root);
  }

  filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			     strlen(fileroot)+strlen(".")+strlen(idx_filesuffix)+
			     /*for interval_char*/1+strlen(positions_suffix)+1,sizeof(char));
  sprintf(filename,"%s/%s.%s%c%s",genomesubdir,fileroot,idx_filesuffix,best_interval_char,positions_suffix);

#ifndef HAVE_MMAP
  new->positions = (Genomicpos_T *) NULL;
  new->positions_fd = Access_fileio(filename);
  new->positions_access = FILEIO;
#else
  if (batch_positions_p == true) {
    if (snps_root) {
      fprintf(stderr,"Pre-loading %s%c (%s) positions db...",idx_filesuffix,best_interval_char,snps_root);
    } else {
      fprintf(stderr,"Pre-loading %s%c positions db...",idx_filesuffix,best_interval_char);
    }
    new->positions = (Genomicpos_T *) Access_mmap_and_preload(&new->positions_fd,&new->positions_len,&npages,&seconds,
							    filename,sizeof(Genomicpos_T));
    if (new->positions == NULL) {
      fprintf(stderr,"insufficient memory (will use disk file instead)\n");
      new->positions_access = FILEIO;
    } else {
      fprintf(stderr,"done (%lu bytes, %d pages, %.2f sec)\n",
	      (long unsigned int) new->positions_len,npages,seconds);
      new->positions_access = MMAPPED;
    }
  } else {
    new->positions = (Genomicpos_T *) Access_mmap(&new->positions_fd,&new->positions_len,
						  filename,sizeof(Genomicpos_T),/*randomp*/true);
    if (new->positions == NULL) {
      new->positions_access = FILEIO;
    } else {
      new->positions_access = MMAPPED;
    }
  }
#endif

#ifdef HAVE_PTHREAD
  if (new->positions_access == FILEIO) {
    pthread_mutex_init(&new->positions_read_mutex,NULL);
  }
#endif

  if (snps_root != NULL) {
    FREE(positions_suffix);
  }

  FREE(filename);

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

static void
positions_write (int fd, Genomicpos_T value) {
  char buffer[4];

  buffer[0] = value & 0xff;
  buffer[1] = (value >>= 8) & 0xff;
  buffer[2] = (value >>= 8) & 0xff;
  buffer[3] = (value >>= 8) & 0xff;

  write(fd,buffer,4);
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
  char buffer[4];

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

  debug0(printf("%06X (%s)\n",oligo,shortoligo_nt(oligo,INDEX1PART)));
  part0 = oligo & LOW12MER;

#if 0
  /* Ignore poly A and poly T on stage 1 */
  if (part0 == 0U || part0 == LOW12MER) {
    *nentries = 0;
    return NULL;
  }
#endif

  switch (this->offsets_access) {
  case ALLOCATED:
    ptr0 = this->offsets[part0];
    end0 = this->offsets[part0+1];
    break;

  case MMAPPED:
#ifdef WORDS_BIGENDIAN
    ptr0 = Bigendian_convert_uint(this->offsets[part0]);
    end0 = Bigendian_convert_uint(this->offsets[part0+1]);
#else
    ptr0 = this->offsets[part0];
    end0 = this->offsets[part0+1];
#endif
    break;

  case FILEIO:
    abort();
  }

  debug0(printf("Indexdb_read: offset pointers are %u and %u\n",ptr0,end0));

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


/* GSNAP version.  Expects calling procedure to handle bigendian conversion. */
Genomicpos_T *
Indexdb_read_inplace (int *nentries, T this, Storedoligomer_T oligo) {
  Positionsptr_T ptr0, end0;
  Storedoligomer_T part0;

  debug0(printf("%06X (%s)\n",oligo,shortoligo_nt(oligo,INDEX1PART)));
  part0 = oligo & LOW12MER;

  switch (this->offsets_access) {
  case ALLOCATED:
    ptr0 = this->offsets[part0];
    end0 = this->offsets[part0+1];
    break;
    
  case MMAPPED:
#ifdef WORDS_BIGENDIAN
    ptr0 = Bigendian_convert_uint(this->offsets[part0]);
    end0 = Bigendian_convert_uint(this->offsets[part0+1]);
#else
    ptr0 = this->offsets[part0];
    end0 = this->offsets[part0+1];
#endif
    break;

  case FILEIO:
    abort();
  }

  debug0(printf("Indexdb_read_inplace: offset pointers are %u and %u\n",ptr0,end0));

  *nentries = end0 - ptr0;

  debug0(
	 printf("%d entries:",*nentries);
	 for (ptr = ptr0; ptr < end0; ptr++) {
	   printf(" %u",this->positions[ptr]);
	 }
	 printf("\n");
	 );

  if (*nentries == 0) {
    return NULL;
  } else {
    return &(this->positions[ptr0]);
  }
}


/* Analogous to Indexdb_read, except this includes diagterm.  Always allocates memory. */
Genomicpos_T *
Indexdb_read_with_diagterm (int *nentries, T this, Storedoligomer_T oligo, int diagterm) {
  Genomicpos_T *positions;
  Positionsptr_T ptr0, end0, ptr;
  int i;

  switch (this->offsets_access) {
  case ALLOCATED:
    ptr0 = this->offsets[oligo];
    end0 = this->offsets[oligo+1];
    break;

  case MMAPPED:
#ifdef WORDS_BIGENDIAN
    ptr0 = Bigendian_convert_uint(this->offsets[oligo]);
    end0 = Bigendian_convert_uint(this->offsets[oligo+1]);
#else
    ptr0 = this->offsets[oligo];
    end0 = this->offsets[oligo+1];
#endif
    break;

  case FILEIO:
    abort();
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

  switch (this->offsets_access) {
  case ALLOCATED:
    ptr0 = this->offsets[oligo];
    end0 = this->offsets[oligo+1];
    break;

  case MMAPPED:
#ifdef WORDS_BIGENDIAN
    ptr0 = Bigendian_convert_uint(this->offsets[oligo]);
    end0 = Bigendian_convert_uint(this->offsets[oligo+1]);
#else
    ptr0 = this->offsets[oligo];
    end0 = this->offsets[oligo+1];
#endif
    break;

  case FILEIO:
    abort();
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

/* Another MONITOR_INTERVAL is in compress.c */
#define MONITOR_INTERVAL 10000000 /* 10 million nt */

/*               87654321 */
#define LEFT_0 0x00000000
#define LEFT_1 0x40000000
#define LEFT_2 0x80000000
#define LEFT_3 0xC0000000
#define LEFT_BIT 0x80000000

#if 0
void
Indexdb_write_genomedibase (FILE *out, FILE *fp) {
  UINT4 low = 0U, high = 0U, flags = 0U, carry;
  Genomicpos_T position = 0U;
  int c, lastc;
  int in_counter = 0;


  lastc = Compress_get_char(fp,position++,/*genome_lc_p*/false);
  while ((c = Compress_get_char(fp,position,/*genome_lc_p*/false)) != EOF) {
    if (isalpha(c)) {
      in_counter++;

      carry = high & 3U;
      high >>= 2;
      low >>= 2;
      flags >>= 1;
      switch (carry) {
      case 0U: break;
      case 1U: low |= LEFT_1; break;
      case 2U: low |= LEFT_2; break;
      case 3U: low |= LEFT_3; break;
      default: abort();
      }

      switch (c) {
      case 'A': 
      switch (lastc) {
      case 'A': break;
      case 'C': high |= LEFT_1; break;
      case 'G': high |= LEFT_2; break;
      case 'T': high |= LEFT_3; break;
      default: flags |= LEFT_BIT; break;
      }
      break;
      case 'C':
	switch (lastc) {
	case 'A': high |= LEFT_1; break;
	case 'C': break;
	case 'G': high |= LEFT_3; break;
	case 'T': high |= LEFT_2; break;
	default: flags |= LEFT_BIT; break;
      }
      break;
      case 'G':
	switch (lastc) {
	case 'A': high |= LEFT_2; break;
	case 'C': high |= LEFT_3; break;
	case 'G': break;
	case 'T': high |= LEFT_1; break;
	default: flags |= LEFT_BIT; break;
      }
      break;
      case 'T':
	switch (lastc) {
	case 'A': high |= LEFT_3; break;
	case 'C': high |= LEFT_2; break;
	case 'G': high |= LEFT_1; break;
	case 'T': break;
	default: flags |= LEFT_BIT; break;
	}
	break;

      default: 
	/* fprintf(stderr,"Non-standard nucleotide %c at position %u.  Using N instead\n",c,position); */
	flags |= LEFT_BIT;
	break;
      }
      
      if (in_counter == 8*sizeof(Genomicpos_T)) {
	FWRITE_UINT(high,out);
	FWRITE_UINT(low,out);
	FWRITE_UINT(flags,out);

	low = high = flags = 0U;
	in_counter = 0;
      }
    }
    position++;
    if (position % MONITOR_INTERVAL == 0) {
      fprintf(stderr,"Compressing position %u\n",position);
    }

    lastc = c;
  }

  if (in_counter > 0) {
    while (in_counter < 8*sizeof(Genomicpos_T)) {
      carry = high & 3U;
      high >>= 2;
      low >>= 2;
      flags >>= 1;
      switch (carry) {
      case 0U: break;
      case 1U: low |= LEFT_1; break;
      case 2U: low |= LEFT_2; break;
      case 3U: low |= LEFT_3; break;
      default: abort();
      }
      high |= LEFT_3; flags |= LEFT_BIT;
      in_counter++;
    }

    FWRITE_UINT(high,out);
    FWRITE_UINT(low,out);
    FWRITE_UINT(flags,out);
  }

  return;
}
#endif


void
Indexdb_write_offsets (FILE *offsets_fp, FILE *sequence_fp, IIT_T chromosome_iit,
		       int index1interval, bool genome_lc_p, char *fileroot,
		       bool mask_lowercase_p) {
  char *uppercaseCode;
  Positionsptr_T *offsets;
  char *comma;
  int c, lastc = 'X', nchrs, chrnum, oligospace, i;
  Genomicpos_T position = 0U, chrpos = 0U, next_chrbound;
  int between_counter = 0, in_counter = 0;
  Storedoligomer_T oligo = 0U, masked, mask;

  if (mask_lowercase_p == false) {
    uppercaseCode = UPPERCASE_U2T; /* We are reading DNA sequence */
  } else {
    uppercaseCode = NO_UPPERCASE;
  }

  mask = ~(~0UL << 2*INDEX1PART);
  oligospace = power(4,INDEX1PART);
  offsets = (Positionsptr_T *) CALLOC(oligospace+1,sizeof(Positionsptr_T));

  /* Handle reference strain */
  chrnum = 1;
  nchrs = IIT_total_nintervals(chromosome_iit);
  next_chrbound = IIT_interval_high(chromosome_iit,chrnum);

  while ((c = Compress_get_char(sequence_fp,position,genome_lc_p)) != EOF) {
    between_counter++;
    in_counter++;

    if (position % MONITOR_INTERVAL == 0) {
      comma = Genomicpos_commafmt(position);
      fprintf(stderr,"Indexing offsets of oligomers in genome %s (every %d bp), position %s",
	      fileroot,index1interval,comma);
      FREE(comma);
      fprintf(stderr,"\n");
    }

    switch (uppercaseCode[c]) {
    case 'A': 
      switch (lastc) {
      case 'A': oligo = (oligo << 2)/*|0U*/; break;
      case 'C': oligo = (oligo << 2) | 1U; break;
      case 'G': oligo = (oligo << 2) | 2U; break;
      case 'T': oligo = (oligo << 2) | 3U; break;
      default: 	oligo = 0U; in_counter = 0;
      }
      break;
    case 'C':
      switch (lastc) {
      case 'A': oligo = (oligo << 2) | 1U; break;
      case 'C': oligo = (oligo << 2)/*|0U*/; break;
      case 'G': oligo = (oligo << 2) | 3U; break;
      case 'T': oligo = (oligo << 2) | 2U; break;
      default: 	oligo = 0U; in_counter = 0;
      }
      break;
    case 'G':
      switch (lastc) {
      case 'A': oligo = (oligo << 2) | 2U; break;
      case 'C': oligo = (oligo << 2) | 3U; break;
      case 'G': oligo = (oligo << 2)/*|0U*/; break;
      case 'T': oligo = (oligo << 2) | 1U; break;
      default: 	oligo = 0U; in_counter = 0;
      }
      break;
    case 'T':
      switch (lastc) {
      case 'A': oligo = (oligo << 2) | 3U; break;
      case 'C': oligo = (oligo << 2) | 2U; break;
      case 'G': oligo = (oligo << 2) | 1U; break;
      case 'T': oligo = (oligo << 2)/*|0U*/; break;
      default: 	oligo = 0U; in_counter = 0;
      }
      break;
    default: oligo = 0U; in_counter = 0; break;
    }

    if (in_counter == INDEX1PART + EXTRA_FOR_DIBASE) {
      if ((chrpos-INDEX1PART-EXTRA_FOR_DIBASE+1U) % index1interval == 0) {
	masked = oligo & mask;
	offsets[masked + 1U] += 1;
	debug(printf("Found oligo %06X.  Incremented offsets for %d to be %d\n",
		     masked,masked+1U,offsets[masked+1U]));
	between_counter = 0;
      }
      in_counter--;
    }

    chrpos++;			/* Needs to go here, before we reset chrpos to 0 */
    if (position >= next_chrbound) {
      oligo = 0U; in_counter = 0;
      chrpos = 0U;
      chrnum++;
      while (chrnum <= nchrs && (next_chrbound = IIT_interval_high(chromosome_iit,chrnum)) < position) {
	chrnum++;
      }
    }
    position++;
    lastc = c;
  }

  for (i = 1; i <= oligospace; i++) {
    offsets[i] = offsets[i] + offsets[i-1];
    debug(if (offsets[i] != offsets[i-1]) {
	    printf("Offset for %06X: %u\n",i,offsets[i]);
	  });
  }

  /*
  fprintf(stderr,"Offset for A...A is %u to %u\n",offsets[0],offsets[1]);
  fprintf(stderr,"Offset for T...T is %u to %u\n",offsets[oligospace-1],offsets[oligospace]);
  */

  fprintf(stderr,"Writing %u offsets to file with total of %d positions\n",oligospace+1,offsets[oligospace]);
  FWRITE_UINTS(offsets,oligospace+1,offsets_fp);
  FREE(offsets);

  return;
}


/* Works directly in file, so we don't need to allocate memory */
static void
compute_positions_in_file (int positions_fd, FILE *offsets_fp, Positionsptr_T *offsets,
			   FILE *sequence_fp, IIT_T chromosome_iit, int index1interval,
			   bool genome_lc_p, char *fileroot, bool mask_lowercase_p) {
  char *uppercaseCode;
  Genomicpos_T position = 0U, chrpos = 0U, next_chrbound;
  char *comma;
  int c, lastc = 'X', nchrs, chrnum;
  int oligospace;
  int between_counter = 0, in_counter = 0;
  Storedoligomer_T oligo = 0U, masked, mask;

  oligospace = power(4,INDEX1PART);

  if (mask_lowercase_p == false) {
    uppercaseCode = UPPERCASE_U2T; /* We are reading DNA sequence */
  } else {
    uppercaseCode = NO_UPPERCASE;
  }

  mask = ~(~0UL << 2*INDEX1PART);

  /* Handle reference strain */
  chrnum = 1;
  nchrs = IIT_total_nintervals(chromosome_iit);
  next_chrbound = IIT_interval_high(chromosome_iit,chrnum);

  while ((c = Compress_get_char(sequence_fp,position,genome_lc_p)) != EOF) {
    between_counter++;
    in_counter++;

    if (position % MONITOR_INTERVAL == 0) {
      comma = Genomicpos_commafmt(position);
      fprintf(stderr,"Indexing positions of oligomers in genome %s (every %d bp), position %s",
	      fileroot,index1interval,comma);
      FREE(comma);
      fprintf(stderr,"\n");
    }

    switch (uppercaseCode[c]) {
    case 'A': 
      switch (lastc) {
      case 'A': oligo = (oligo << 2)/*|0U*/; break;
      case 'C': oligo = (oligo << 2) | 1U; break;
      case 'G': oligo = (oligo << 2) | 2U; break;
      case 'T': oligo = (oligo << 2) | 3U; break;
      default: 	oligo = 0U; in_counter = 0;
      }
      break;
    case 'C':
      switch (lastc) {
      case 'A': oligo = (oligo << 2) | 1U; break;
      case 'C': oligo = (oligo << 2)/*|0U*/; break;
      case 'G': oligo = (oligo << 2) | 3U; break;
      case 'T': oligo = (oligo << 2) | 2U; break;
      default: 	oligo = 0U; in_counter = 0;
      }
      break;
    case 'G':
      switch (lastc) {
      case 'A': oligo = (oligo << 2) | 2U; break;
      case 'C': oligo = (oligo << 2) | 3U; break;
      case 'G': oligo = (oligo << 2)/*|0U*/; break;
      case 'T': oligo = (oligo << 2) | 1U; break;
      default: 	oligo = 0U; in_counter = 0;
      }
      break;
    case 'T':
      switch (lastc) {
      case 'A': oligo = (oligo << 2) | 3U; break;
      case 'C': oligo = (oligo << 2) | 2U; break;
      case 'G': oligo = (oligo << 2) | 1U; break;
      case 'T': oligo = (oligo << 2)/*|0U*/; break;
      default: 	oligo = 0U; in_counter = 0;
      }
      break;
    default: oligo = 0U; in_counter = 0; break;
    }

    /*
    debug(printf("char=%c bc=%d ic=%d oligo=%016lX\n",
		 c,between_counter,in_counter,oligo));
    */
    
    if (in_counter == INDEX1PART + EXTRA_FOR_DIBASE) {
      if ((chrpos-INDEX1PART-EXTRA_FOR_DIBASE+1U) % index1interval == 0) {
	masked = oligo & mask;
	positions_move_absolute(positions_fd,offsets[masked]);
	offsets[masked] += 1;
	positions_write(positions_fd,position-INDEX1PART-EXTRA_FOR_DIBASE+1U);
	between_counter = 0;
      }
      in_counter--;
    }
    
    chrpos++;			/* Needs to go here, before we reset chrpos to 0 */
    if (position >= next_chrbound) {
      oligo = 0U; in_counter = 0;
      chrpos = 0U;
      chrnum++;
      while (chrnum <= nchrs && (next_chrbound = IIT_interval_high(chromosome_iit,chrnum)) < position) {
	chrnum++;
      }
    }
    position++;
    lastc = c;
  }

  return;
}

/* Requires sufficient memory to hold all positions */
static void
compute_positions_in_memory (Genomicpos_T *positions, FILE *offsets_fp, Positionsptr_T *offsets,
			     FILE *sequence_fp, IIT_T chromosome_iit, int index1interval,
			     bool genome_lc_p, char *fileroot, bool mask_lowercase_p) {
  char *uppercaseCode;
  Genomicpos_T position = 0U, chrpos = 0U, next_chrbound;
  char *comma;
  int c, lastc = 'X', nchrs, chrnum;
  int oligospace;
  int between_counter = 0, in_counter = 0;
  Storedoligomer_T oligo = 0U, masked, mask;
  debug1(char *nt);

  oligospace = power(4,INDEX1PART);

  if (mask_lowercase_p == false) {
    uppercaseCode = UPPERCASE_U2T; /* We are reading DNA sequence */
  } else {
    uppercaseCode = NO_UPPERCASE;
  }

  mask = ~(~0UL << 2*INDEX1PART);

  /* Handle reference strain */
  chrnum = 1;
  nchrs = IIT_total_nintervals(chromosome_iit);
  next_chrbound = IIT_interval_high(chromosome_iit,chrnum);

  while ((c = Compress_get_char(sequence_fp,position,genome_lc_p)) != EOF) {
    between_counter++;
    in_counter++;

    if (position % MONITOR_INTERVAL == 0) {
      comma = Genomicpos_commafmt(position);
      fprintf(stderr,"Indexing positions of oligomers in genome %s (every %d bp), position %s",
	      fileroot,index1interval,comma);
      FREE(comma);
      fprintf(stderr,"\n");
    }

    switch (uppercaseCode[c]) {
    case 'A': 
      switch (lastc) {
      case 'A': oligo = (oligo << 2)/*|0U*/; break;
      case 'C': oligo = (oligo << 2) | 1U; break;
      case 'G': oligo = (oligo << 2) | 2U; break;
      case 'T': oligo = (oligo << 2) | 3U; break;
      default: 	oligo = 0U; in_counter = 0;
      }
      break;
    case 'C':
      switch (lastc) {
      case 'A': oligo = (oligo << 2) | 1U; break;
      case 'C': oligo = (oligo << 2)/*|0U*/; break;
      case 'G': oligo = (oligo << 2) | 3U; break;
      case 'T': oligo = (oligo << 2) | 2U; break;
      default: 	oligo = 0U; in_counter = 0;
      }
      break;
    case 'G':
      switch (lastc) {
      case 'A': oligo = (oligo << 2) | 2U; break;
      case 'C': oligo = (oligo << 2) | 3U; break;
      case 'G': oligo = (oligo << 2)/*|0U*/; break;
      case 'T': oligo = (oligo << 2) | 1U; break;
      default: 	oligo = 0U; in_counter = 0;
      }
      break;
    case 'T':
      switch (lastc) {
      case 'A': oligo = (oligo << 2) | 3U; break;
      case 'C': oligo = (oligo << 2) | 2U; break;
      case 'G': oligo = (oligo << 2) | 1U; break;
      case 'T': oligo = (oligo << 2)/*|0U*/; break;
      default: 	oligo = 0U; in_counter = 0;
      }
      break;
    default: oligo = 0U; in_counter = 0; break;
    }

    debug(printf("char=%c bc=%d ic=%d oligo=%08X\n",
		 c,between_counter,in_counter,oligo));

    if (in_counter == INDEX1PART + EXTRA_FOR_DIBASE) {
      if ((chrpos-INDEX1PART-EXTRA_FOR_DIBASE+1U) % index1interval == 0) {
	masked = oligo & mask;
	positions[offsets[masked]++] = position-INDEX1PART-EXTRA_FOR_DIBASE+1U;
	debug1(nt = shortoligo_nt(masked,INDEX1PART);
	       printf("Storing %s at %u, chrpos %u\n",
		      nt,position-INDEX1PART-EXTRA_FOR_DIBASE+1U,chrpos-INDEX1PART-EXTRA_FOR_DIBASE+1U);
	       FREE(nt));
	between_counter = 0;
      }
      in_counter--;
    }
    
    chrpos++;			/* Needs to go here, before we reset chrpos to 0 */
    if (position >= next_chrbound) {
      debug1(printf("Skipping because position %u is at chrbound\n",position));
      oligo = 0U; in_counter = 0;
      chrpos = 0U;
      chrnum++;
      while (chrnum <= nchrs && (next_chrbound = IIT_interval_high(chromosome_iit,chrnum)) < position) {
	chrnum++;
      }
    }
    position++;
    lastc = c;
  }

  return;
}


void
Indexdb_write_positions (char *positionsfile, FILE *offsets_fp, FILE *sequence_fp,
			 IIT_T chromosome_iit, int index1interval,
			 bool genome_lc_p, bool writefilep, char *fileroot, bool mask_lowercase_p) {
  FILE *positions_fp;		/* For building positions in memory */
  int positions_fd;		/* For building positions in file */
  Positionsptr_T *offsets = NULL, totalcounts;
  Genomicpos_T *positions;
  int oligospace;

  oligospace = power(4,INDEX1PART);

  offsets = (Positionsptr_T *) CALLOC(oligospace+1,sizeof(Positionsptr_T));
  FREAD_UINTS(offsets,oligospace+1,offsets_fp);
  totalcounts = offsets[oligospace];
  if (totalcounts == 0) {
    fprintf(stderr,"Something is wrong with the offsets file.  Total counts is zero.\n");
    exit(9);
  }

  if (writefilep == true) {
    fprintf(stderr,"User requested build of positions in file\n");
    positions_fd = Access_fileio_rw(positionsfile);
    compute_positions_in_file(positions_fd,offsets_fp,offsets,sequence_fp,chromosome_iit,
			      index1interval,genome_lc_p,fileroot,mask_lowercase_p);
    close(positions_fd);

  } else {
    fprintf(stderr,"Trying to allocate %u*%d bytes of memory...",totalcounts,(int) sizeof(Genomicpos_T));
    positions = (Genomicpos_T *) CALLOC_NO_EXCEPTION(totalcounts,sizeof(Genomicpos_T));
    if (positions == NULL) {
      fprintf(stderr,"failed.  Building positions in file.\n");
      positions_fd = Access_fileio_rw(positionsfile);
      compute_positions_in_file(positions_fd,offsets_fp,offsets,sequence_fp,chromosome_iit,
				index1interval,genome_lc_p,fileroot,mask_lowercase_p);
      close(positions_fd);

    } else {
      fprintf(stderr,"succeeded.  Building positions in memory.\n");
      if ((positions_fp = FOPEN_WRITE_BINARY(positionsfile)) == NULL) {
	fprintf(stderr,"Can't open file %s\n",positionsfile);
	exit(9);
      }
      compute_positions_in_memory(positions,offsets_fp,offsets,sequence_fp,chromosome_iit,
				  index1interval,genome_lc_p,fileroot,mask_lowercase_p);
      fprintf(stderr,"Writing %u genomic positions to file %s ...\n",
	      totalcounts,positionsfile);
      FWRITE_UINTS(positions,totalcounts,positions_fp);

      fclose(positions_fp);
      FREE(positions);
    }
  }

  if (offsets != NULL) {
    FREE(offsets);
  }

  return;
}


/************************************************************************
 *   Create procedure -- for user-provided genomic segment
 ************************************************************************/

T
Indexdb_new_segment (char *genomicseg, int index1interval) {
  T new = (T) MALLOC(sizeof(*new));
  char *uppercaseCode;
  Positionsptr_T *work_offsets;	/* Working set for use in calculating positions */
  int totalcounts = 0;
  int c, lastc = 'X', oligospace, i;
  char *p;
  Genomicpos_T position = 0U;
  int between_counter = 0, in_counter = 0;
  Storedoligomer_T oligo = 0U, masked, mask;

  uppercaseCode = UPPERCASE_U2T;

  mask = ~(~0UL << 2*INDEX1PART);
  oligospace = power(4,INDEX1PART);
  new->index1interval = 1;

  /* Create offsets */
  new->offsets = (Positionsptr_T *) CALLOC(oligospace+1,sizeof(Positionsptr_T));
  new->offsets_access = ALLOCATED;

  p = genomicseg;
  while ((c = *(p++)) != '\0') {
    between_counter++;
    in_counter++;

    switch (uppercaseCode[c]) {
    case 'A': 
      switch (lastc) {
      case 'A': oligo = (oligo << 2)/*|0U*/; break;
      case 'C': oligo = (oligo << 2) | 1U; break;
      case 'G': oligo = (oligo << 2) | 2U; break;
      case 'T': oligo = (oligo << 2) | 3U; break;
      default: 	oligo = 0U; in_counter = 0;
      }
      break;
    case 'C':
      switch (lastc) {
      case 'A': oligo = (oligo << 2) | 1U; break;
      case 'C': oligo = (oligo << 2)/*|0U*/; break;
      case 'G': oligo = (oligo << 2) | 3U; break;
      case 'T': oligo = (oligo << 2) | 2U; break;
      default: 	oligo = 0U; in_counter = 0;
      }
      break;
    case 'G':
      switch (lastc) {
      case 'A': oligo = (oligo << 2) | 2U; break;
      case 'C': oligo = (oligo << 2) | 3U; break;
      case 'G': oligo = (oligo << 2)/*|0U*/; break;
      case 'T': oligo = (oligo << 2) | 1U; break;
      default: 	oligo = 0U; in_counter = 0;
      }
      break;
    case 'T':
      switch (lastc) {
      case 'A': oligo = (oligo << 2) | 3U; break;
      case 'C': oligo = (oligo << 2) | 2U; break;
      case 'G': oligo = (oligo << 2) | 1U; break;
      case 'T': oligo = (oligo << 2)/*|0U*/; break;
      default: 	oligo = 0U; in_counter = 0;
      }
      break;
    default: oligo = 0U; in_counter = 0; break;
    }

    /*
    debug(printf("char=%c bc=%d ic=%d oligo=%016lX\n",
		 c,between_counter,in_counter,oligo));
    */

    if (in_counter == INDEX1PART + EXTRA_FOR_DIBASE) {
      if ((position-INDEX1PART-EXTRA_FOR_DIBASE+1U) % index1interval == 0) {
	masked = oligo & mask;
	new->offsets[masked + 1U] += 1;
	debug(printf("Found oligo %06X.  Incremented offsets for %d to be %d\n",
		     masked,masked+1U,new->offsets[masked+1U]));
	between_counter = 0;
      }
      in_counter--;
    }

    position++;
    lastc = c;
  }

  for (i = 1; i <= oligospace; i++) {
    new->offsets[i] = new->offsets[i] + new->offsets[i-1];
    debug(printf("Offset for %06X: %u\n",i,new->offsets[i]));
  }


  /* Create positions */
  position = 0U;
  between_counter = in_counter = 0;
  oligo = 0UL;
  lastc = 'X';

  work_offsets = (Positionsptr_T *) CALLOC(oligospace+1,sizeof(Positionsptr_T));
  for (i = 0; i <= oligospace; i++) {
    work_offsets[i] = new->offsets[i];
  }

  totalcounts = new->offsets[oligospace];
  if (totalcounts == 0) {
    fprintf(stderr,"Error: user-provided genomic segment has no valid oligomers of size %d\n",INDEX1PART);
    exit(9);
  }
  new->positions = (Genomicpos_T *) CALLOC(totalcounts,sizeof(Genomicpos_T));
  new->positions_access = ALLOCATED;

  p = genomicseg;
  while ((c = *(p++)) != '\0') {
    between_counter++;
    in_counter++;

    switch (uppercaseCode[c]) {
    case 'A': 
      switch (lastc) {
      case 'A': oligo = (oligo << 2)/*|0U*/; break;
      case 'C': oligo = (oligo << 2) | 1U; break;
      case 'G': oligo = (oligo << 2) | 2U; break;
      case 'T': oligo = (oligo << 2) | 3U; break;
      default: 	oligo = 0U; in_counter = 0;
      }
      break;
    case 'C':
      switch (lastc) {
      case 'A': oligo = (oligo << 2) | 1U; break;
      case 'C': oligo = (oligo << 2)/*|0U*/; break;
      case 'G': oligo = (oligo << 2) | 3U; break;
      case 'T': oligo = (oligo << 2) | 2U; break;
      default: 	oligo = 0U; in_counter = 0;
      }
      break;
    case 'G':
      switch (lastc) {
      case 'A': oligo = (oligo << 2) | 2U; break;
      case 'C': oligo = (oligo << 2) | 3U; break;
      case 'G': oligo = (oligo << 2)/*|0U*/; break;
      case 'T': oligo = (oligo << 2) | 1U; break;
      default: 	oligo = 0U; in_counter = 0;
      }
      break;
    case 'T':
      switch (lastc) {
      case 'A': oligo = (oligo << 2) | 3U; break;
      case 'C': oligo = (oligo << 2) | 2U; break;
      case 'G': oligo = (oligo << 2) | 1U; break;
      case 'T': oligo = (oligo << 2)/*|0U*/; break;
      default: 	oligo = 0U; in_counter = 0;
      }
      break;
    default: oligo = 0U; in_counter = 0; break;
    }

    /*
    debug(printf("char=%c bc=%d ic=%d oligo=%06X\n",
		 c,between_counter,in_counter,oligo));
    */
    
    if (in_counter == INDEX1PART + EXTRA_FOR_DIBASE) {
      if ((position-INDEX1PART-EXTRA_FOR_DIBASE+1U) % index1interval == 0) {
	masked = oligo & mask;
	new->positions[work_offsets[masked]++] = position-INDEX1PART-EXTRA_FOR_DIBASE+1U;
	between_counter = 0;
      }
      in_counter--;
    }

    position++;
    lastc = c;
  }

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

