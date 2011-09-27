static char rcsid[] = "$Id: indexdb.c,v 1.138 2010/02/03 19:07:14 twu Exp $";
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
#ifdef HAVE_PTHREAD
#ifdef PMAP
      pthread_mutex_destroy(&(*old)->offsets_read_mutex);
#endif
#endif
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

#ifdef PMAP
  n = oligospace = power(NAMINOACIDS,INDEX1PART_AA);
#else
  n = oligospace = power(4,INDEX1PART);
  if (cmetp == true) {
    n = power(3,INDEX1PART);
  }
#endif

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


#ifdef PMAP
#ifndef HAVE_MMAP
  new->offsets = (Positionsptr_T *) NULL;
  new->offsets_fd = Access_fileio(filename);
  new->offsets_access = FILEIO;
#else
  if (batch_offsets_p == true) {
    if (snps_root) {
      fprintf(stderr,"Pre-loading %s%c (%s) offsets db...",idx_filesuffix,best_interval_char,snps_root);
    } else {
      fprintf(stderr,"Pre-loading %s%c offsets db...",idx_filesuffix,best_interval_char);
    }
    new->offsets = (Positionsptr_T *) Access_mmap_and_preload(&new->offsets_fd,&new->offsets_len,&npages,&seconds,
							    filename,sizeof(Positionsptr_T));
    if (new->offsets == NULL) {
      fprintf(stderr,"insufficient memory (will use disk file instead)\n");
      new->offsets_access = FILEIO;
    } else {
      fprintf(stderr,"done (%lu bytes, %d pages, %.2f sec)\n",
	      new->offsets_len,npages,seconds);
      new->offsets_access = MMAPPED;
    }
  } else {
    new->offsets = (Positionsptr_T *) Access_mmap(&new->offsets_fd,&new->offsets_len,
						  filename,sizeof(Positionsptr_T),/*randomp*/false);
    if (new->offsets == NULL) {
      new->offsets_access = FILEIO;
    } else {
      new->offsets_access = MMAPPED;
    }
  }

#ifdef HAVE_PTHREAD
  if (new->offsets_access == FILEIO) {
    pthread_mutex_init(&new->offsets_read_mutex,NULL);
  }
#endif

#endif

#else  /* PMAP */

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

#ifdef PMAP

#if NAMINOACIDS == 20
static char aa_table[NAMINOACIDS] = "ACDEFGHIKLMNPQRSTVWY";
#elif NAMINOACIDS == 15
static char aa_table[NAMINOACIDS] = "ACDFGHIKMNPRSWY";
#elif NAMINOACIDS == 12
static char aa_table[NAMINOACIDS] = "ACDEFGHIKPSW";
#endif

static char *
aaindex_aa (unsigned int aaindex) {
  char *aa;
  int i, j;
  int aasize = INDEX1PART_AA;

  aa = (char *) CALLOC(aasize+1,sizeof(char));
  j = aasize-1;
  for (i = 0; i < aasize; i++) {
    aa[j] = aa_table[aaindex % NAMINOACIDS];
    aaindex /= NAMINOACIDS;
    j--;
  }

  return aa;
}

/*               87654321 */
#define LEFT_A 0x00000000
#define LEFT_C 0x40000000
#define LEFT_G 0x80000000
#define LEFT_T 0xC0000000

static char *
highlow_nt (Storedoligomer_T high, Storedoligomer_T low) {
  char *nt;
  int i, j;
  Storedoligomer_T lowbits;
  int oligosize = INDEX1PART_NT;

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


#else

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


/************************************************************************
 *   Read procedures
 ************************************************************************/

/* PMAP only, because GMAP does allocation rather than fileio */
#ifdef PMAP
static void
offsets_move_absolute (int offsets_fd, unsigned int aaindex) {
  off_t offset = aaindex*((off_t) sizeof(Positionsptr_T));

  if (lseek(offsets_fd,offset,SEEK_SET) < 0) {
    fprintf(stderr,"Attempted to do lseek on offset %u*%lu=%lu\n",aaindex,sizeof(Positionsptr_T),offset);
    perror("Error in indexdb.c, offsets_move_absolute");
    exit(9);
  }
  return;
}

static Positionsptr_T
offsets_read_forward (int offsets_fd) {
  Positionsptr_T value;
  char buffer[4];

  read(offsets_fd,buffer,4);

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

  switch (this->offsets_access) {
  case FILEIO:
#ifdef HAVE_PTHREAD
    pthread_mutex_lock(&this->offsets_read_mutex);
#endif
    offsets_move_absolute(this->offsets_fd,aaindex);
    ptr0 = offsets_read_forward(this->offsets_fd);
    end0 = offsets_read_forward(this->offsets_fd);
#ifdef HAVE_PTHREAD
    pthread_mutex_unlock(&this->offsets_read_mutex);
#endif
    break;

  case ALLOCATED:
    ptr0 = this->offsets[aaindex];
    end0 = this->offsets[aaindex+1];
    break;

  case MMAPPED:
#ifdef WORDS_BIGENDIAN
    ptr0 = Bigendian_convert_uint(this->offsets[aaindex]);
    end0 = Bigendian_convert_uint(this->offsets[aaindex+1]);
#else
    ptr0 = this->offsets[aaindex];
    end0 = this->offsets[aaindex+1];
#endif
    break;
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

#if 0
  debug0(printf("%06X (%s)\n",oligo,shortoligo_nt(oligo,INDEX1PART)));
#endif
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
#ifdef DEBUG0
  Positionsptr_T ptr;
#endif

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

#endif	/* ifdef PMAP */


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
get_aa_index (Storedoligomer_T high, Storedoligomer_T low, bool watsonp) {
  unsigned int aaindex = 0U;
  Storedoligomer_T shifted;
  int i, codonindex;
  char *nt, *aa;

  if (watsonp == true) {
    for (i = INDEX1PART_NT-3; i >= 0; i -= 3) {
      shifted = offset_codon(high,low,i);
      if ((codonindex = get_codon_fwd(shifted)) == AA_STOP) {
	fprintf(stderr,"Unexpected stop codon in get_aa_index\n");
	abort();
      } else {
	aaindex = aaindex * NAMINOACIDS + codonindex;
      }
    }
  } else {
    for (i = 0; i < INDEX1PART_NT; i += 3) {
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
Indexdb_write_offsets (FILE *offsets_fp, FILE *sequence_fp, IIT_T chromosome_iit,
		       int index1interval,
#ifdef PMAP
		       bool watsonp,
#endif
		       bool genome_lc_p, char *fileroot, bool mask_lowercase_p) {
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
#else
  int between_counter = 0, in_counter = 0;
  Storedoligomer_T oligo = 0U, masked, mask;
#endif

  if (mask_lowercase_p == false) {
    uppercaseCode = UPPERCASE_U2T; /* We are reading DNA sequence */
  } else {
    uppercaseCode = NO_UPPERCASE;
  }

#ifdef PMAP
  oligospace = power(NAMINOACIDS,INDEX1PART_AA);
  between_counter[0] = between_counter[1] = between_counter[2] = 0;
  in_counter[0] = in_counter[1] = in_counter[2] = 0;
#else
  mask = ~(~0UL << 2*INDEX1PART);
  oligospace = power(4,INDEX1PART);
#endif
  offsets = (Positionsptr_T *) CALLOC(oligospace+1,sizeof(Positionsptr_T));

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
      fprintf(stderr,"Indexing offsets of oligomers in genome %s (every %d aa), position %s",
	      fileroot,index1interval,comma);
#else
      fprintf(stderr,"Indexing offsets of oligomers in genome %s (every %d bp), position %s",
	      fileroot,index1interval,comma);
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
    if (in_counter[frame] == INDEX1PART_AA + 1) {
      if (between_counter[frame] >= index1interval) {
	aaindex = get_aa_index(high,low,watsonp);
	offsets[aaindex + 1U] += 1;
	debug1(
	       aa = aaindex_aa(aaindex);
	       if (watsonp == true) {
		 printf("Storing %s (%u) at %u\n",aa,aaindex,position-INDEX1PART_NT+1U);
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
    if (in_counter == INDEX1PART) {
      if (
#ifdef NONMODULAR
	  between_counter >= index1interval
#else
	  (chrpos-INDEX1PART+1U) % index1interval == 0
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

  fprintf(stderr,"Writing %u offsets to file with total of %d positions\n",oligospace+1,offsets[oligospace]);
  FWRITE_UINTS(offsets,oligospace+1,offsets_fp);
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
compute_positions_in_file (int positions_fd, FILE *offsets_fp, Positionsptr_T *offsets,
			   FILE *sequence_fp, IIT_T chromosome_iit, int index1interval,
#ifdef PMAP
			   bool watsonp,
#endif
			   bool genome_lc_p, char *fileroot, bool mask_lowercase_p) {
  char *uppercaseCode;
  Genomicpos_T position = 0U, chrpos = 0U, next_chrbound;
  char *comma;
  int c, nchrs, chrnum;
#ifdef PMAP
  int frame = -1, between_counter[3], in_counter[3];
  Storedoligomer_T high = 0U, low = 0U, carry;
  unsigned int aaindex;
  Genomicpos_T adjposition;
#else
  int oligospace;
  int between_counter = 0, in_counter = 0;
  Storedoligomer_T oligo = 0U, masked, mask;
#endif

#ifndef PMAP
  oligospace = power(4,INDEX1PART);
#endif

  if (mask_lowercase_p == false) {
    uppercaseCode = UPPERCASE_U2T; /* We are reading DNA sequence */
  } else {
    uppercaseCode = NO_UPPERCASE;
  }

#ifdef PMAP
  /* oligospace = power(NAMINOACIDS,INDEX1PART_AA); */
#else
  mask = ~(~0UL << 2*INDEX1PART);
  /* oligospace = power(4,INDEX1PART); */
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
      fprintf(stderr,"Indexing positions of oligomers in genome %s (every %d aa), position %s",
	      fileroot,index1interval,comma);
#else
      fprintf(stderr,"Indexing positions of oligomers in genome %s (every %d bp), position %s",
	      fileroot,index1interval,comma);
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
    if (in_counter[frame] == INDEX1PART_AA + 1) {
      if (between_counter[frame] >= index1interval) {
	aaindex = get_aa_index(high,low,watsonp);
	positions_move_absolute(positions_fd,offsets[aaindex]);
	offsets[aaindex] += 1;
	if (watsonp == true) {
	  adjposition = position-INDEX1PART_NT+1U;
	} else {
	  adjposition = position;
	}
	WRITE_UINT(adjposition,positions_fd);
	between_counter[frame] = 0;
      }
      in_counter[frame] -= 1;
    }
#else
    if (in_counter == INDEX1PART) {
      if (
#ifdef NONMODULAR
	  between_counter >= index1interval
#else
	  (chrpos-INDEX1PART+1U) % index1interval == 0
#endif
	  ) {
	masked = oligo & mask;
	positions_move_absolute(positions_fd,offsets[masked]);
	offsets[masked] += 1;
	WRITE_UINT(position-INDEX1PART+1U,positions_fd);
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
compute_positions_in_memory (Genomicpos_T *positions, FILE *offsets_fp, Positionsptr_T *offsets,
			     FILE *sequence_fp, IIT_T chromosome_iit, int index1interval,
#ifdef PMAP
			     bool watsonp,
#endif
			     bool genome_lc_p, char *fileroot, bool mask_lowercase_p) {
  char *uppercaseCode;
  Genomicpos_T position = 0U, chrpos = 0U, next_chrbound;
  char *comma;
  int c, nchrs, chrnum;
#ifdef PMAP
  int frame = -1, between_counter[3], in_counter[3];
  Storedoligomer_T high = 0U, low = 0U, carry;
  unsigned int aaindex;
  debug1(char *aa);
#else
  int between_counter = 0, in_counter = 0;
  Storedoligomer_T oligo = 0U, masked, mask;
  debug1(char *nt);
#endif

#ifdef ADDSENTINEL
  int oligospace;
  oligospace = power(4,INDEX1PART);
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
  mask = ~(~0UL << 2*INDEX1PART);
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
      fprintf(stderr,"Indexing positions of oligomers in genome %s (every %d aa), position %s",
	      fileroot,index1interval,comma);
#else
      fprintf(stderr,"Indexing positions of oligomers in genome %s (every %d bp), position %s",
	      fileroot,index1interval,comma);
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
    if (in_counter[frame] == INDEX1PART_AA + 1) {
      if (between_counter[frame] >= index1interval) {
	aaindex = get_aa_index(high,low,watsonp);
	if (watsonp == true) {
	  positions[offsets[aaindex]++] = position-INDEX1PART_NT+1U;
	  debug1(adjposition = position-INDEX1PART_NT+1U);
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
    if (in_counter == INDEX1PART) {
      if (
#ifdef NONMODULAR
	  between_counter >= index1interval
#else
	  (chrpos-INDEX1PART+1U) % index1interval == 0
#endif
	  ) {
	masked = oligo & mask;
	positions[offsets[masked]++] = position-INDEX1PART+1U;
	debug1(nt = shortoligo_nt(masked,INDEX1PART);
	       printf("Storing %s at %u, chrpos %u\n",nt,position-INDEX1PART+1U,chrpos-INDEX1PART+1U);
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
Indexdb_write_positions (char *positionsfile, FILE *offsets_fp, FILE *sequence_fp,
			 IIT_T chromosome_iit, int index1interval,
#ifdef PMAP
			 bool watsonp,
#endif
			 bool genome_lc_p, bool writefilep, char *fileroot, bool mask_lowercase_p) {
  FILE *positions_fp;		/* For building positions in memory */
  int positions_fd;		/* For building positions in file */
  Positionsptr_T *offsets = NULL, totalcounts;
  Genomicpos_T *positions;
  int oligospace;

#ifdef PMAP
  oligospace = power(NAMINOACIDS,INDEX1PART_AA);
#else
  oligospace = power(4,INDEX1PART);
#endif

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
#ifdef PMAP
    compute_positions_in_file(positions_fd,offsets_fp,offsets,sequence_fp,chromosome_iit,
			      index1interval,watsonp,genome_lc_p,fileroot,mask_lowercase_p);
#else
    compute_positions_in_file(positions_fd,offsets_fp,offsets,sequence_fp,chromosome_iit,
			      index1interval,genome_lc_p,fileroot,mask_lowercase_p);
#endif
    close(positions_fd);

  } else {
    fprintf(stderr,"Trying to allocate %u*%d bytes of memory...",totalcounts,(int) sizeof(Genomicpos_T));
    positions = (Genomicpos_T *) CALLOC_NO_EXCEPTION(totalcounts,sizeof(Genomicpos_T));
    if (positions == NULL) {
      fprintf(stderr,"failed.  Building positions in file.\n");
      positions_fd = Access_fileio_rw(positionsfile);
#ifdef PMAP
      compute_positions_in_file(positions_fd,offsets_fp,offsets,sequence_fp,chromosome_iit,
				index1interval,watsonp,genome_lc_p,fileroot,mask_lowercase_p);
#else
      compute_positions_in_file(positions_fd,offsets_fp,offsets,sequence_fp,chromosome_iit,
				index1interval,genome_lc_p,fileroot,mask_lowercase_p);
#endif
      close(positions_fd);

    } else {
      fprintf(stderr,"succeeded.  Building positions in memory.\n");
      if ((positions_fp = FOPEN_WRITE_BINARY(positionsfile)) == NULL) {
	fprintf(stderr,"Can't open file %s\n",positionsfile);
	exit(9);
      }
#ifdef PMAP
      compute_positions_in_memory(positions,offsets_fp,offsets,sequence_fp,chromosome_iit,
				  index1interval,watsonp,genome_lc_p,fileroot,mask_lowercase_p);
#else
      compute_positions_in_memory(positions,offsets_fp,offsets,sequence_fp,chromosome_iit,
				  index1interval,genome_lc_p,fileroot,mask_lowercase_p);
#endif
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
Indexdb_new_segment (char *genomicseg, int index1interval
#ifdef PMAP
		     , bool watsonp
#endif
		     ) {
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
#else
  int between_counter = 0, in_counter = 0;
  Storedoligomer_T oligo = 0U, masked, mask;
#endif

  uppercaseCode = UPPERCASE_U2T;

#ifdef PMAP
  oligospace = power(NAMINOACIDS,INDEX1PART_AA);
#else
  mask = ~(~0UL << 2*INDEX1PART);
  oligospace = power(4,INDEX1PART);
  new->index1interval = 1;
#endif

  /* Create offsets */
  new->offsets = (Positionsptr_T *) CALLOC(oligospace+1,sizeof(Positionsptr_T));
  new->offsets_access = ALLOCATED;

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
    if (in_counter[frame] == INDEX1PART_AA + 1) {
      if (between_counter[frame] >= index1interval) {
	aaindex = get_aa_index(high,low,watsonp);
	new->offsets[aaindex + 1U] += 1;
	between_counter[frame] = 0;
      }
      in_counter[frame] -= 1;
    }
#else
    if (in_counter == INDEX1PART) {
      if (
#ifdef NONMODULAR
	  between_counter >= index1interval
#else
	  /* Actually, modular condition not needed for user-supplied genomic segment */
	  (position-INDEX1PART+1U) % index1interval == 0
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
    fprintf(stderr,"Error: user-provided genomic segment has no valid oligomers of size %d\n",INDEX1PART_NT);
#else
    fprintf(stderr,"Error: user-provided genomic segment has no valid oligomers of size %d\n",INDEX1PART);
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

    if (in_counter[frame] == INDEX1PART_AA + 1) {
      if (between_counter[frame] >= index1interval) {
	aaindex = get_aa_index(high,low,watsonp);
	if (watsonp == true) {
	  new->positions[work_offsets[aaindex]++] = position-INDEX1PART_NT+1U;
	} else {
	  new->positions[work_offsets[aaindex]++] = position;
	}
	between_counter[frame] = 0;
      }
      in_counter[frame] -= 1;
    }
#else
    if (in_counter == INDEX1PART) {
      if (
#ifdef NONMODULAR
	  between_counter >= index1interval
#else
	  /* Actually, modular condition not needed for user-supplied genomic segment */
	  (position-INDEX1PART+1U) % index1interval == 0
#endif
	  ) {
	masked = oligo & mask;
	new->positions[work_offsets[masked]++] = position-INDEX1PART+1U;
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

