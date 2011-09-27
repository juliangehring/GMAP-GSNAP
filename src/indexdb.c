static char rcsid[] = "$Id: indexdb.c,v 1.97 2006/11/30 17:17:00 twu Exp $";
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

#include "mem.h"
#include "fopen.h"

#include "compress.h"
#include "access.h"
#include "interval.h"
#include "complement.h"

#ifdef HAVE_PTHREAD
#include <pthread.h>		/* sys/types.h already included above */
#endif

/* An offset into the positions file of an IndexDB.  Typically, 3
   billion divided by sampling interval, requiring a maximum of 32
   bits or 4 bytes */
typedef UINT4 Positionsptr_T;


#define MAXENTRIES 20
#define BADVAL (Genomicpos_T) -1

/* Low-level codon hacking */
#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Output of positions */
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

/* get_codon_rev */
#ifdef DEBUG3
#define debug3(x) x
#else
#define debug3(x)
#endif


#define T Indexdb_T
struct T {
  Access_T offsets_access;
  int offsets_fd;
  size_t offsets_len;
  Positionsptr_T *offsets;
#ifdef HAVE_PTHREAD
#ifdef PMAP
  pthread_mutex_t offsets_read_mutex;
#endif
#endif

  Access_T positions_access;
  int positions_fd;
  size_t positions_len;
  Genomicpos_T *positions;
#ifdef HAVE_PTHREAD
  pthread_mutex_t positions_read_mutex;
#endif
};


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


T
Indexdb_new_genome (char *genomesubdir, char *fileroot,
#ifdef PMAP
		    bool watsonp,
#endif
		    bool batch_offsets_p, bool batch_positions_p) {
  T new = (T) MALLOC(sizeof(*new));
  char *filename;
  size_t len;
  int npages;
  double seconds;

  /* Read offsets file */
#ifdef PMAP
  if (watsonp == true) {
    filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			       strlen(fileroot)+strlen(".")+strlen(PFXOFFSETS)+1,sizeof(char));
    sprintf(filename,"%s/%s.%s",genomesubdir,fileroot,PFXOFFSETS);
  } else {
    filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			       strlen(fileroot)+strlen(".")+strlen(PRXOFFSETS)+1,sizeof(char));
    sprintf(filename,"%s/%s.%s",genomesubdir,fileroot,PRXOFFSETS);
  }
#else
  filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			     strlen(fileroot)+strlen(".")+strlen(IDXOFFSETS)+1,sizeof(char));
  sprintf(filename,"%s/%s.%s",genomesubdir,fileroot,IDXOFFSETS);
#endif

#ifdef PMAP
#ifndef HAVE_MMAP
  new->offsets = (Positionsptr_T *) NULL;
  new->offsets_fd = Access_fileio(filename);
  new->offsets_access = FILEIO;
#else
  if (batch_offsets_p == true) {
    if (watsonp == true) {
      fprintf(stderr,"Pre-loading forward offsets db...");
    } else {
      fprintf(stderr,"Pre-loading reverse offsets db...");
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

#else

#ifndef HAVE_MMAP
  new->offsets = (Positionsptr_T *) Access_allocated(&len,&seconds,filename,sizeof(Positionsptr_T));
  new->offsets_access = ALLOCATED;
#else
  if (batch_offsets_p == true) {
    fprintf(stderr,"Pre-reading offsets db...");
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


  /* Read positions file */
#ifdef PMAP
  if (watsonp == true) {
    filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			       strlen(fileroot)+strlen(".")+strlen(PFXPOSITIONS)+1,sizeof(char));
    sprintf(filename,"%s/%s.%s",genomesubdir,fileroot,PFXPOSITIONS);
  } else {
    filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			       strlen(fileroot)+strlen(".")+strlen(PRXPOSITIONS)+1,sizeof(char));
    sprintf(filename,"%s/%s.%s",genomesubdir,fileroot,PRXPOSITIONS);
  }
#else
  filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			     strlen(fileroot)+strlen(".")+strlen(IDXPOSITIONS)+1,sizeof(char));
  sprintf(filename,"%s/%s.%s",genomesubdir,fileroot,IDXPOSITIONS);
#endif

#ifndef HAVE_MMAP
  new->positions = (Genomicpos_T *) NULL;
  new->positions_fd = Access_fileio(filename);
  new->positions_access = FILEIO;
#else
  if (batch_positions_p == true) {
#ifdef PMAP
    if (watsonp == true) {
      fprintf(stderr,"Pre-loading forward positions db...");
    } else {
      fprintf(stderr,"Pre-loading reverse positions db...");
    }
#else
    fprintf(stderr,"Pre-loading positions db...");
#endif
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
    fprintf(stderr,"Attempted to do lseek on offset %u*%d=%lu\n",aaindex,sizeof(Positionsptr_T),offset);
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
    fprintf(stderr,"Attempted to do lseek on offset %u*%d=%lu\n",
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

static void
positions_write_multiple (int fd, Genomicpos_T *values, int n) {
  int i;
  Genomicpos_T value;
  char buffer[4];

  for (i = 0; i < n; i++) {
    value = values[i];
    buffer[0] = value & 0xff;
    buffer[1] = (value >>= 8) & 0xff;
    buffer[2] = (value >>= 8) & 0xff;
    buffer[3] = (value >>= 8) & 0xff;

    write(fd,buffer,4);
  }

  return;
}


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

#ifdef PMAP

/* PMAP version */
Genomicpos_T *
Indexdb_read (int *nentries, T this, unsigned int aaindex) {
  Genomicpos_T *positions, position0, expected1, position1;
  Positionsptr_T ptr, ptr0, end0;
  int i;
  char byte1, byte2, byte3;
  int bigendian, littleendian;

  debug(printf("%u (%s)\n",aaindex,aaindex_aa(aaindex)));

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
  debug(printf("offset pointers are %u and %u\n",ptr0,end0));

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

  *nentries = end0 - ptr0;

  if (*nentries == 0) {
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
    debug(
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
/*                 87654321 */
#define LOW12MER 0x00FFFFFF

/* GMAP version */
Genomicpos_T *
Indexdb_read (int *nentries, T this, Storedoligomer_T oligo) {
  Genomicpos_T *positions;
  Positionsptr_T ptr0, end0;
  Storedoligomer_T part0;
  int i;
#ifdef WORDS_BIGENDIAN
  Positionsptr_T ptr;
  int bigendian, littleendian;
  char byte1, byte2, byte3;
#endif

  debug(printf("%06X (%s)\n",oligo,shortoligo_nt(oligo,INDEX1PART)));
  part0 = oligo & LOW12MER;

  /* Ignore poly A and poly T on stage 1 */
  if (part0 == 0U || part0 == LOW12MER) {
    return NULL;
  }

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

  debug(printf("offset pointers are %u and %u\n",ptr0,end0));

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

  *nentries = end0 - ptr0;

  if (*nentries == 0) {
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
    debug(
	  printf("%d entries:",*nentries);
	  for (i = 0; i < *nentries; i++) {
	    printf(" %u",positions[i]);
	  }
	  printf("\n");
	  );

    return positions;
  }
}

#endif

/************************************************************************
 *   Write procedures -- called by gmapindex/pmapindex
 ************************************************************************/

static int
power (int base, int exponent) {
  int result = 1, i;

  for (i = 0; i < exponent; i++) {
    result *= base;
  }
  return result;
}

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


static char uppercaseCode[128] = UPPERCASE_U2T;

/* Another MONITOR_INTERVAL is in compress.c */
#define MONITOR_INTERVAL 10000000 /* 10 million nt */

void
Indexdb_write_offsets (FILE *offsets_fp, FILE *sequence_fp, IIT_T altstrain_iit,
		       int index1interval,
#ifdef PMAP
		       bool watsonp,
#endif
		       bool uncompressedp, char *fileroot) {
  Positionsptr_T *offsets;
  char *comma, *p, b;
  bool allocp;
  int c, oligospace, index, i;
  Genomicpos_T position = 0U;
#ifdef PMAP
  int frame = -1, between_counter[3], in_counter[3];
  Storedoligomer_T high = 0U, low = 0U, carry;
  unsigned int aaindex;
  char *aa;
#else
  int between_counter = 0, in_counter = 0;
  Storedoligomer_T oligo = 0U, masked, mask;
#endif
  Interval_T interval;

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
  while ((c = Compress_get_char(sequence_fp,position,uncompressedp)) != EOF) {
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
      in_counter[0] = 0; in_counter[1] = in_counter[2] = 0;
      break;
    default: 
      if (uncompressedp == true) {
	high = low = carry = 0U;
	in_counter[0] = 0; in_counter[1] = in_counter[2] = 0;
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
      if (uncompressedp == true) {
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
      if (between_counter >= index1interval) {
	masked = oligo & mask;
	offsets[masked + 1U] += 1;
	debug(printf("Found oligo %06X.  Incremented offsets for %d to be %d\n",
		     masked,masked+1U,offsets[masked+1U]));
	between_counter = 0;
      }
      in_counter--;
    }
#endif

    position++;
  }

  /* Handle alternate strains */
  if (altstrain_iit != NULL) {
    for (index = 1; index <= IIT_nintervals(altstrain_iit); index++) {
      interval = IIT_interval(altstrain_iit,index);

      position = Interval_low(interval);
      p = IIT_annotation(altstrain_iit,index,&allocp); /* Holds the sequence */
      
#ifdef PMAP
      frame = -1;
      between_counter[0] = between_counter[1] = between_counter[2] = 0;
      in_counter[0] = in_counter[1] = in_counter[2] = 0;
      high = low = 0U;
#else
      between_counter = in_counter = 0;
      oligo = 0U;
#endif
      while ((b = *p++) != '\0') {
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
	switch (uppercaseCode[b]) {
	case 'A': low = (low << 2); break;
	case 'C': low = (low << 2) | 1U; break;
	case 'G': low = (low << 2) | 2U; break;
	case 'T': low = (low << 2) | 3U; break;
	default:
	  high = low = carry = 0U;
	  in_counter[0] = 0; in_counter[1] = in_counter[2] = 0;
	  break;
	}
	high = (high << 2) | carry; 
#else
	switch (uppercaseCode[b]) {
	case 'A': oligo = (oligo << 2); break;
	case 'C': oligo = (oligo << 2) | 1; break;
	case 'G': oligo = (oligo << 2) | 2; break;
	case 'T': oligo = (oligo << 2) | 3; break;
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
	    offsets[aaindex + 1U] += 1;
	    between_counter[frame] = 0;
	  }
	  in_counter[frame] -= 1;
	}
#else
	if (in_counter == INDEX1PART) {
	  if (between_counter >= index1interval) {
	    masked = oligo & mask;
	    offsets[masked + 1U] += 1;
	    debug(printf("Found oligo %06X.  Incremented offsets for %d to be %d\n",
			 masked,masked+1U,offsets[masked+1U]));
	    between_counter = 0;
	  }
	  in_counter--;
	}
#endif

	position++;
      }

      if (allocp == true) {
	FREE(p);
      }
    }
  }

  for (i = 1; i <= oligospace; i++) {
    offsets[i] = offsets[i] + offsets[i-1];
    debug(if (offsets[i] != offsets[i-1]) {
	    printf("Offset for %06X: %u\n",i,offsets[i]);
	  });
  }

  FWRITE_UINTS(offsets,oligospace+1,offsets_fp);
  FREE(offsets);

  return;
}


/* FILE *fp is preferable to int fd, because former is buffered.  No
   need for fseeko, because offsets file is < 2 Gigabytes */
static void
offsetsfile_move_absolute (FILE *fp, int ptr) {
  long int offset = ptr*((long int) sizeof(Positionsptr_T));

  if (fseek(fp,offset,SEEK_SET) < 0) {
    fprintf(stderr,"Attempted to do fseek on offset %u*%d=%lu\n",ptr,sizeof(Positionsptr_T),offset);
    perror("Error in indexdb.c, offsetsfile_move_absolute");
    exit(9);
  }

  return;
}


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



/* Works directly in file, so we don't need to allocate memory */
static void
compute_positions_in_file (int positions_fd, FILE *offsets_fp, Positionsptr_T *offsets,
			   FILE *sequence_fp, IIT_T altstrain_iit, int index1interval,
#ifdef PMAP
			   bool watsonp,
#endif
			   bool uncompressedp, char *fileroot) {
  Positionsptr_T block_start, block_end;
  Genomicpos_T *positions_for_block, position = 0U, adjposition, prevpos;
  char *comma, *p, b;
  bool allocp;
  int c, oligospace, index, i, npositions, nbadvals, j, k;
#ifdef PMAP
  int frame = -1, between_counter[3], in_counter[3];
  Storedoligomer_T high = 0U, low = 0U, carry;
  unsigned int aaindex;
#else
  int between_counter = 0, in_counter = 0;
  Storedoligomer_T oligo = 0U, masked, mask;
#endif
  Interval_T interval;

#ifdef PMAP
  oligospace = power(NAMINOACIDS,INDEX1PART_AA);
#else
  mask = ~(~0UL << 2*INDEX1PART);
  oligospace = power(4,INDEX1PART);
#endif

  /* Handle reference strain */
  while ((c = Compress_get_char(sequence_fp,position,uncompressedp)) != EOF) {
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
    switch (c) {
    case 'A': low = (low << 2); break;
    case 'C': low = (low << 2) | 1; break;
    case 'G': low = (low << 2) | 2; break;
    case 'T': low = (low << 2) | 3; break;
    case 'X': case 'N': 
      high = low = carry = 0U;
      in_counter[0] = 0; in_counter[1] = in_counter[2] = 0;
      break;
    default: 
      if (uncompressedp == true) {
	high = low = carry = 0U;
	in_counter[0] = 0; in_counter[1] = in_counter[2] = 0;
      } else {
	fprintf(stderr,"Bad character %c at position %u\n",c,position);
	abort();
      }
    }
    high = (high << 2) | carry; 
#else
    switch (c) {
    case 'A': oligo = (oligo << 2); break;
    case 'C': oligo = (oligo << 2) | 1; break;
    case 'G': oligo = (oligo << 2) | 2; break;
    case 'T': oligo = (oligo << 2) | 3; break;
    case 'X': case 'N': oligo = 0U; in_counter = 0; break;
    default: 
      if (uncompressedp == true) {
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
	positions_write(positions_fd,adjposition);
	between_counter[frame] = 0;
      }
      in_counter[frame] -= 1;
    }
#else
    if (in_counter == INDEX1PART) {
      if (between_counter >= index1interval) {
	masked = oligo & mask;
	positions_move_absolute(positions_fd,offsets[masked]);
	offsets[masked] += 1;
	adjposition = position-INDEX1PART+1U;
	positions_write(positions_fd,adjposition);
	between_counter = 0;
      }
      in_counter--;
    }
#endif
    
    position++;
  }

  /* Handle alternate strains */
  if (altstrain_iit != NULL) {
    for (index = 1; index <= IIT_nintervals(altstrain_iit); index++) {
      interval = IIT_interval(altstrain_iit,index);
      position = Interval_low(interval);
      p = IIT_annotation(altstrain_iit,index,&allocp); /* Holds the sequence */
      
#ifdef PMAP
      frame = -1;
      between_counter[0] = between_counter[1] = between_counter[2] = 0;
      in_counter[0] = in_counter[1] = in_counter[2] = 0;
      high = low = 0U;
#else
      between_counter = in_counter = 0;
      oligo = 0U;
#endif
      while ((b = *p++) != '\0') {
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
	switch (uppercaseCode[b]) {
	case 'A': low = (low << 2); break;
	case 'C': low = (low << 2) | 1; break;
	case 'G': low = (low << 2) | 2; break;
	case 'T': low = (low << 2) | 3; break;
	default:
	  high = low = carry = 0U;
	  in_counter[0] = 0; in_counter[1] = in_counter[2] = 0;
	  break;
	}
	high = (high << 2) | carry; 
#else
	switch (uppercaseCode[b]) {
	case 'A': oligo = (oligo << 2); break;
	case 'C': oligo = (oligo << 2) | 1; break;
	case 'G': oligo = (oligo << 2) | 2; break;
	case 'T': oligo = (oligo << 2) | 3; break;
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
	    positions_move_absolute(positions_fd,offsets[aaindex]);
	    offsets[aaindex] += 1;
	    if (watsonp == true) {
	      adjposition = position-INDEX1PART_NT+1U;
	    } else {
	      adjposition = position;
	    }
	    positions_write(positions_fd,adjposition);
	    between_counter[frame] = 0;
	  }
	  in_counter[frame] -= 1;
	}
#else
	if (in_counter == INDEX1PART) {
	  if (between_counter >= index1interval) {
	    masked = oligo & mask;
	    positions_move_absolute(positions_fd,offsets[masked]);
	    offsets[masked] += 1;
	    adjposition = position-INDEX1PART+1U;
	    positions_write(positions_fd,adjposition);
	    between_counter = 0;
	  }
	  in_counter--;
	}
#endif

	position++;
      }
      if (allocp == true) {
	FREE(p);
      }
    }

    /* Warning: at this point, offsets is not the original offsets, so
       need to re-read them. */
    offsetsfile_move_absolute(offsets_fp,0);
    FREAD_UINTS(offsets,oligospace+1,offsets_fp);

    /* Sorting positions (necessary if alternate strains are present) */
    /* Also removing duplicates.  ?Done for efficiency, since stage1.c
       can handle them. */
    for (i = 0; i < oligospace; i++) {
      block_start = offsets[i];
      block_end = offsets[i+1];
      if ((npositions = block_end - block_start) > 1) {
	positions_for_block = (Genomicpos_T *) CALLOC(npositions,sizeof(Genomicpos_T));
	positions_move_absolute(positions_fd,block_start);
	positions_read_multiple(positions_fd,positions_for_block,npositions);

	if (need_to_sort_p(positions_for_block,npositions)) {	
	  qsort(positions_for_block,npositions,sizeof(Genomicpos_T),Genomicpos_compare);

	  /* Mark duplicates */
	  prevpos = positions_for_block[0];
	  nbadvals = 0;
	  for (j = 1; j < npositions; j++) {
	    if (positions_for_block[j] == prevpos) {
	      /* fprintf(stderr,"Found duplicate position %u in offset %d\n",prevpos,i); */
	      positions_for_block[j] = BADVAL;
	      nbadvals++;
	    } else {
	      prevpos = positions_for_block[j];
	    }
	  }

	  /* Move unique positions forward and bad values to end */
	  if (nbadvals > 0) {
	    k = 1;
	    for (j = 1; j < npositions; j++) {
	      if (positions_for_block[j] != BADVAL) {
		positions_for_block[k++] = positions_for_block[j];
	      }
	    }
	    for ( ; k < npositions; k++) {
	      positions_for_block[k] = BADVAL;
	    }
	  }

	  positions_move_absolute(positions_fd,block_start);
	  positions_write_multiple(positions_fd,positions_for_block,npositions);
	}

	FREE(positions_for_block);
      }
    }
  }

  return;
}

/* Requires sufficient memory to hold all positions */
static void
compute_positions_in_memory (Genomicpos_T *positions, FILE *offsets_fp, Positionsptr_T *offsets,
			     FILE *sequence_fp, IIT_T altstrain_iit, int index1interval,
#ifdef PMAP
			     bool watsonp,
#endif
			     bool uncompressedp, char *fileroot) {
  Positionsptr_T block_start, block_end, j, k;
  Genomicpos_T position = 0U, prevpos;
  char *comma, *p, b;
  bool allocp;
  int c, oligospace, index, i, npositions, nbadvals;
#ifdef PMAP
  int frame = -1, between_counter[3], in_counter[3];
  Storedoligomer_T high = 0U, low = 0U, carry;
  unsigned int aaindex;
  debug1(char *aa);
#else
  int between_counter = 0, in_counter = 0;
  Storedoligomer_T oligo = 0U, masked, mask;
#endif
  Interval_T interval;

#ifdef PMAP
  oligospace = power(NAMINOACIDS,INDEX1PART_AA);
  between_counter[0] = between_counter[1] = between_counter[2] = 0;
  in_counter[0] = in_counter[1] = in_counter[2] = 0;
#else
  mask = ~(~0UL << 2*INDEX1PART);
  oligospace = power(4,INDEX1PART);
#endif

  /* Handle reference strain */
  while ((c = Compress_get_char(sequence_fp,position,uncompressedp)) != EOF) {
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
    switch (c) {
    case 'A': low = (low << 2); break;
    case 'C': low = (low << 2) | 1U; break;
    case 'G': low = (low << 2) | 2U; break;
    case 'T': low = (low << 2) | 3U; break;
    case 'X': case 'N':
      high = low = carry = 0U;
      in_counter[0] = 0; in_counter[1] = in_counter[2] = 0;
      break;
    default: 
      if (uncompressedp == true) {
	high = low = carry = 0U;
	in_counter[0] = 0; in_counter[1] = in_counter[2] = 0;
      } else {
	fprintf(stderr,"Bad character %c at position %u\n",c,position);
	abort();
      }
    }
    high = (high << 2) | carry; 

    debug(printf("frame=%d char=%c bc=%d ic=%d high=%08X low=%08X\n",
		 frame,c,between_counter[frame],in_counter[frame],high,low));
#else
    switch (c) {
    case 'A': oligo = (oligo << 2); break;
    case 'C': oligo = (oligo << 2) | 1; break;
    case 'G': oligo = (oligo << 2) | 2; break;
    case 'T': oligo = (oligo << 2) | 3; break;
    case 'X': case 'N': oligo = 0U; in_counter = 0; break;
    default: 
      if (uncompressedp == true) {
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
	  debug1(
		 aa = aaindex_aa(aaindex);
		 printf("Storing %s (%u) at %u\n",aa,aaindex,position-INDEX1PART_NT+1U);
		 FREE(aa);
		 );
	} else {
	  positions[offsets[aaindex]++] = position;
	  debug1(
		 aa = aaindex_aa(aaindex);
		 printf("Storing %s (%u) at %u\n",aa,aaindex,position);
		 FREE(aa);
		 );
	}
	between_counter[frame] = 0;
      }
      in_counter[frame] -= 1;
    }
#else
    if (in_counter == INDEX1PART) {
      if (between_counter >= index1interval) {
	masked = oligo & mask;
	positions[offsets[masked]++] = position-INDEX1PART+1U;
	between_counter = 0;
      }
      in_counter--;
    }
#endif
    
    position++;
  }

  /* Handle alternate strains */
  if (altstrain_iit != NULL) {
    for (index = 1; index <= IIT_nintervals(altstrain_iit); index++) {
      interval = IIT_interval(altstrain_iit,index);
      position = Interval_low(interval);
      p = IIT_annotation(altstrain_iit,index,&allocp); /* Holds the sequence */
      
#ifdef PMAP
      frame = -1;
      between_counter[0] = between_counter[1] = between_counter[2] = 0;
      in_counter[0] = in_counter[1] = in_counter[2] = 0;
      high = low = 0U;
#else
      between_counter = in_counter = 0;
      oligo = 0U;
#endif
      while ((b = *p++) != '\0') {
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
	switch (uppercaseCode[b]) {
	case 'A': low = (low << 2); break;
	case 'C': low = (low << 2) | 1U; break;
	case 'G': low = (low << 2) | 2U; break;
	case 'T': low = (low << 2) | 3U; break;
	default:
	  high = low = carry = 0U;
	  in_counter[0] = 0; in_counter[1] = in_counter[2] = 0;
	  break;
	}
	high = (high << 2) | carry; 
#else
	switch (uppercaseCode[b]) {
	case 'A': oligo = (oligo << 2); break;
	case 'C': oligo = (oligo << 2) | 1; break;
	case 'G': oligo = (oligo << 2) | 2; break;
	case 'T': oligo = (oligo << 2) | 3; break;
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
	    if (watsonp == true) {
	      positions[offsets[aaindex]++] = position-INDEX1PART_NT+1U;
	    } else {
	      positions[offsets[aaindex]++] = position;
	    }
	    between_counter[frame] = 0;
	  }
	  in_counter[frame] -= 1;
	}
#else
	if (in_counter == INDEX1PART) {
	  if (between_counter >= index1interval) {
	    masked = oligo & mask;
	    positions[offsets[masked]++] = position-INDEX1PART+1U;
	    between_counter = 0;
	  }
	  in_counter--;
	}
#endif
	position++;
      }

      if (allocp == true) {
	FREE(p);
      }
    }

    /* Warning: at this point, offsets is not the original offsets, so
       need to re-read them. */
    offsetsfile_move_absolute(offsets_fp,0);
    FREAD_UINTS(offsets,oligospace+1,offsets_fp);

    /* Sorting positions (necessary if alternate strains are present) */
    /* Also removing duplicates.  ?Done for efficiency, since stage1.c
       can handle them. */
    for (i = 0; i < oligospace; i++) {
      block_start = offsets[i];
      block_end = offsets[i+1];
      if ((npositions = block_end - block_start) > 1) {

	if (need_to_sort_p(&(positions[block_start]),npositions)) {
	  qsort(&(positions[block_start]),npositions,sizeof(Genomicpos_T),Genomicpos_compare);

	  /* Mark duplicates */
	  prevpos = positions[block_start];
	  nbadvals = 0;
	  for (j = block_start+1U; j < block_end; j++) {
	    if (positions[j] == prevpos) {
	      /* fprintf(stderr,"Found duplicate position %u in offset %d\n",prevpos,i); */
	      positions[j] = BADVAL;
	      nbadvals++;
	    } else {
	      prevpos = positions[j];
	    }
	  }

	  /* Move unique positions forward and bad values to end */
	  if (nbadvals > 0) {
	    k = block_start+1U;
	    for (j = block_start+1U; j < block_end; j++) {
	      if (positions[j] != BADVAL) {
		positions[k++] = positions[j];
	      }
	    }
	    for ( ; k < block_end; k++) {
	      positions[k] = BADVAL;
	    }
	  }
	}
      }
    }
  }

  return;
}


void
Indexdb_write_positions (char *positionsfile, FILE *offsets_fp, FILE *sequence_fp, 
			 IIT_T altstrain_iit, int index1interval, bool uncompressedp, 
#ifdef PMAP
			 bool watsonp,
#endif
			 bool writefilep, char *fileroot) {
  FILE *positions_fp;		/* For building positions in memory */
  int positions_fd;		/* For building positions in file */
  Positionsptr_T *offsets, totalcounts;
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
    compute_positions_in_file(positions_fd,offsets_fp,offsets,sequence_fp,altstrain_iit,
			      index1interval,watsonp,uncompressedp,fileroot);
#else
    compute_positions_in_file(positions_fd,offsets_fp,offsets,sequence_fp,altstrain_iit,
			      index1interval,uncompressedp,fileroot);
#endif
    close(positions_fd);

  } else {
    fprintf(stderr,"Trying to allocate %u*%d bytes of memory...",totalcounts,(int) sizeof(Genomicpos_T));
    positions = (Genomicpos_T *) CALLOC_NO_EXCEPTION(totalcounts,sizeof(Genomicpos_T));
    if (positions == NULL) {
      fprintf(stderr,"failed.  Building positions in file.\n");
      positions_fd = Access_fileio_rw(positionsfile);
#ifdef PMAP
      compute_positions_in_file(positions_fd,offsets_fp,offsets,sequence_fp,altstrain_iit,
				index1interval,watsonp,uncompressedp,fileroot);
#else
      compute_positions_in_file(positions_fd,offsets_fp,offsets,sequence_fp,altstrain_iit,
				index1interval,uncompressedp,fileroot);
#endif
      close(positions_fd);

    } else {
      fprintf(stderr,"succeeded.  Building positions in memory.\n");
      if ((positions_fp = FOPEN_WRITE_BINARY(positionsfile)) == NULL) {
	fprintf(stderr,"Can't open file %s\n",positionsfile);
	exit(9);
      }
#ifdef PMAP
      compute_positions_in_memory(positions,offsets_fp,offsets,sequence_fp,altstrain_iit,
				  index1interval,watsonp,uncompressedp,fileroot);
#else
      compute_positions_in_memory(positions,offsets_fp,offsets,sequence_fp,altstrain_iit,
				  index1interval,uncompressedp,fileroot);
#endif
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
Indexdb_new_segment (char *genomicseg, int index1interval
#ifdef PMAP
		     , bool watsonp
#endif
		     ) {
  T new = (T) MALLOC(sizeof(*new));
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

#ifdef PMAP
  oligospace = power(NAMINOACIDS,INDEX1PART_AA);
#else
  mask = ~(~0UL << 2*INDEX1PART);
  oligospace = power(4,INDEX1PART);
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
    case 'C': oligo = (oligo << 2) | 1; break;
    case 'G': oligo = (oligo << 2) | 2; break;
    case 'T': oligo = (oligo << 2) | 3; break;
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
      if (between_counter >= index1interval) {
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

  for (i = 1; i <= oligospace; i++) {
    new->offsets[i] = new->offsets[i] + new->offsets[i-1];
    debug(printf("Offset for %06X: %u\n",i,new->offsets[i]));
  }

  /* Create positions */
  position = 0;
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
      in_counter[0] = 0; in_counter[1] = in_counter[2] = 0;
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
      if (between_counter >= index1interval) {
	masked = oligo & mask;
	new->positions[work_offsets[masked]++] = position-INDEX1PART+1U;
	between_counter = 0;
      }
      in_counter--;
    }
#endif

    position++;
  }

  FREE(work_offsets);

  return new;
}


