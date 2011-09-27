static char rcsid[] = "$Id: genome.c 46773 2011-09-07 22:05:18Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifndef HAVE_MEMCPY
# define memcpy(d,s,n) bcopy((s),(d),(n))
#endif
#ifndef HAVE_MEMMOVE
# define memmove(d,s,n) bcopy((s),(d),(n))
#endif

#include "genome.h"

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>		/* For munmap */
#ifdef HAVE_UNISTD_H
#include <unistd.h>		/* For lseek and close */
#endif
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>		/* For off_t */
#endif

#ifdef HAVE_PTHREAD
#include <pthread.h>		/* sys/types.h already included above */
#endif

#include "assert.h"
#include "except.h"
#include "mem.h"
#include "complement.h"
#include "interval.h"
#include "genomicpos.h"		/* For Genomicpos_commafmt */


#ifdef WORDS_BIGENDIAN
#include "bigendian.h"
#else
#include "littleendian.h"
#endif


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Print genomic segment */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

/* Patching with strain information */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

/* Nucleotides */
#ifdef DEBUG3
#define debug3(x) x
#else
#define debug3(x)
#endif


#define T Genome_T
struct T {
  Access_T access;
  int fd;
  size_t len;

  char *chars;
  UINT4 *blocks;
  bool compressedp;

  char *ptr;
  unsigned int left;
#ifdef HAVE_PTHREAD
  pthread_mutex_t read_mutex;
#endif
};


UINT4 *
Genome_blocks (T this) {
  return this->blocks;
}


Genomicpos_T
Genome_totallength (T this) {
  if (this->compressedp == false) {
    return (Genomicpos_T) this->len;
  } else {
    return (Genomicpos_T) ((this->len/3) * 8);
  }
}

void
Genome_free (T *old) {
  if (*old) {
    if ((*old)->access == ALLOCATED) {
      FREE((*old)->blocks);
#ifdef HAVE_MMAP
    } else if ((*old)->access == MMAPPED) {
      if ((*old)->compressedp == true) {
	munmap((void *) (*old)->blocks,(*old)->len);
      } else {
	munmap((void *) (*old)->chars,(*old)->len);
      }
      close((*old)->fd);
#endif
    } else if ((*old)->access == FILEIO) {
#ifdef HAVE_PTHREAD
      pthread_mutex_destroy(&(*old)->read_mutex);
#endif
      close((*old)->fd);
    }

    FREE(*old);
  }
  return;
}


T
Genome_new (char *genomesubdir, char *fileroot, char *snps_root, bool genome_lc_p, Access_mode_T access) {
  T new = (T) MALLOC(sizeof(*new));
  char *filename;
  bool compressedp = !genome_lc_p;
  char *comma;
  int npages;
  double seconds;

  new->compressedp = compressedp;

  if (compressedp == true) {
    if (snps_root != NULL) {
      filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+strlen(fileroot)+
				 strlen(".genomecomp.")+strlen(snps_root)+1,sizeof(char));
      sprintf(filename,"%s/%s.genomecomp.%s",genomesubdir,fileroot,snps_root);
    } else {
      filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+strlen(fileroot)+
				 strlen(".genomecomp")+1,sizeof(char));
      sprintf(filename,"%s/%s.genomecomp",genomesubdir,fileroot);
    }
  } else {
    filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+strlen(fileroot)+
			       strlen(".genome")+1,sizeof(char));
    sprintf(filename,"%s/%s.genome",genomesubdir,fileroot);
  }

  if (compressedp == true) {
    new->chars = (char *) NULL;

    if (access == USE_ALLOCATE) {
      fprintf(stderr,"Allocating memory for compressed genome...");
      new->blocks = (UINT4 *) Access_allocated(&new->len,&seconds,filename,sizeof(UINT4));
      if (new->blocks == NULL) {
	fprintf(stderr,"insufficient memory (need to use a lower batch mode (-B))\n");
	exit(9);
      } else {
	comma = Genomicpos_commafmt(new->len);
	fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma,seconds);
	FREE(comma);
	new->access = ALLOCATED;
      }

#ifdef HAVE_MMAP
    } else if (access == USE_MMAP_PRELOAD) {
      fprintf(stderr,"Pre-loading compressed genome...");
      new->blocks = (UINT4 *) Access_mmap_and_preload(&new->fd,&new->len,&npages,&seconds,
						      filename,sizeof(UINT4));
      if (new->blocks == NULL) {
	fprintf(stderr,"insufficient memory (will use disk file instead)\n");
	new->access = FILEIO;
      } else {
	comma = Genomicpos_commafmt(new->len);
	fprintf(stderr,"done (%s bytes, %d pages, %.2f sec)\n",comma,npages,seconds);
	FREE(comma);
	new->access = MMAPPED;
      }
    } else if (access == USE_MMAP_ONLY) {
      new->blocks = (UINT4 *) Access_mmap(&new->fd,&new->len,filename,sizeof(UINT4),/*randomp*/false);
      if (new->blocks == NULL) {
	fprintf(stderr,"Insufficient memory for genome mmap (will use disk file instead)\n");
	new->access = FILEIO;
      } else {
	new->access = MMAPPED;
      }
#endif

    } else if (access == USE_FILEIO) {
      new->blocks = (UINT4 *) NULL;
      new->fd = Access_fileio(filename);
      new->access = FILEIO;
    }

  } else {
    new->blocks = (UINT4 *) NULL;

    if (access == USE_ALLOCATE) {
      fprintf(stderr,"Allocating memory for uncompressed genome...");
      new->chars = (char *) Access_allocated(&new->len,&seconds,filename,sizeof(char));
      if (new->chars == NULL) {
	fprintf(stderr,"insufficient memory (need to use a lower batch mode (-B))\n");
	exit(9);
      } else {
	comma = Genomicpos_commafmt(new->len);
	fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma,seconds);
	FREE(comma);
	new->access = ALLOCATED;
      }
      
#ifdef HAVE_MMAP
    } else if (access == USE_MMAP_PRELOAD) {
      fprintf(stderr,"Pre-loading uncompressed genome...");
      new->chars = (char *) Access_mmap_and_preload(&new->fd,&new->len,&npages,&seconds,
						    filename,sizeof(char));
      if (new->chars == NULL) {
	fprintf(stderr,"insufficient memory (will use disk file instead)\n");
	new->access = FILEIO;
      } else {
	comma = Genomicpos_commafmt(new->len);
	fprintf(stderr,"done (%s bytes, %d pages, %.2f sec)\n",comma,npages,seconds);
	FREE(comma);
	new->access = MMAPPED;
      }

    } else if (access == USE_MMAP_ONLY) {
      new->chars = (char *) Access_mmap(&new->fd,&new->len,filename,sizeof(char),/*randomp*/false);
      if (new->chars == NULL) {
	fprintf(stderr,"Insufficient memory for genome mmap (will use disk file instead)\n");
	new->access = FILEIO;
      } else {
	new->access = MMAPPED;
      }
#endif

    } else if (access == USE_FILEIO) {
      new->chars = (char *) NULL;
      new->fd = Access_fileio(filename);
      new->access = FILEIO;
    }
  }

#ifdef HAVE_PTHREAD
  if (new->access == FILEIO) {
    pthread_mutex_init(&new->read_mutex,NULL);
  }
#endif

  FREE(filename);

  /* Initialize for Genome_next_char */
  new->ptr = new->chars;
  new->left = 0U;
  if (new->compressedp == false) {
    if (new->access == FILEIO) {
#ifdef HAVE_PTHREAD
      pthread_mutex_lock(&new->read_mutex);
#endif
      if (lseek(new->fd,/*left*/0U,SEEK_SET) < 0) {
	perror("Error in Genome_new");
	exit(9);
      }
#ifdef HAVE_PTHREAD
      pthread_mutex_unlock(&new->read_mutex);
#endif
    }
  }

  return new;
}


static char complCode[128] = COMPLEMENT_LC;

static void
make_complement_buffered (char *complement, char *sequence, Genomicpos_T length) {
  int i, j;

  for (i = length-1, j = 0; i >= 0; i--, j++) {
    complement[j] = complCode[(int) sequence[i]];
  }
  complement[length] = '\0';
  return;
}

static void
make_complement_inplace (char *sequence, Genomicpos_T length) {
  char temp;
  unsigned int i, j;

  for (i = 0, j = length-1; i < length/2; i++, j--) {
    temp = complCode[(int) sequence[i]];
    sequence[i] = complCode[(int) sequence[j]];
    sequence[j] = temp;
  }
  if (i == j) {
    sequence[i] = complCode[(int) sequence[i]];
  }

  return;
}


/************************************************************************
 *   Read procedures
 ************************************************************************/

static void
genomecomp_move_absolute (T this, Genomicpos_T ptr) {
  off_t offset = ptr*((off_t) sizeof(UINT4));

  if (lseek(this->fd,offset,SEEK_SET) < 0) {
    perror("Error in gmap, genomecomp_move_absolute");
    exit(9);
  }
  return;
}

static UINT4
genomecomp_read_current (T this) {
  UINT4 value;
  char buffer[4];

  read(this->fd,buffer,4);

  value = (buffer[3] & 0xff);
  value <<= 8;
  value |= (buffer[2] & 0xff);
  value <<= 8;
  value |= (buffer[1] & 0xff);
  value <<= 8;
  value |= (buffer[0] & 0xff);

  return value;
}


static char DEFAULT_CHARS[4] = {'A','C','G','T'};
static char DEFAULT_FLAGS[4] = {'N','N','N','N'};

static char SNP_CHARS[4] = {' ',' ',' ',' '};
static char SNP_FLAGS[4] = {'A','C','G','T'};

static char *global_chars = DEFAULT_CHARS;
static char *global_flags = DEFAULT_FLAGS;


#if 0
void
Genome_replace_x (void) {
  non_acgt[3] = 'N';
}
#endif


static void
uncompress_fileio (char *gbuffer1, T this, Genomicpos_T startpos, 
		   Genomicpos_T endpos, const char defaultchars[],
		   const char flagchars[]) {
  /* Genomicpos_T length = endpos - startpos; */
  Genomicpos_T startblock, endblock, ptr;
  UINT4 high, low, flags;
  char Buffer[32];
  int startdiscard, enddiscard, i, k = 0;

  /* sequence = (char *) CALLOC(length+1,sizeof(char)); */

  ptr = startblock = startpos/32U*3;
  endblock = endpos/32U*3;
  startdiscard = startpos % 32;
  enddiscard = endpos % 32;
  
  if (endblock == startblock) {
    /* Special case */
    genomecomp_move_absolute(this,ptr);
    high = genomecomp_read_current(this);
    low = genomecomp_read_current(this);
    flags = genomecomp_read_current(this);

    for (i = 0; i < 16; i++) {
      Buffer[i] = (char) ((flags & 1U) ? flagchars[low & 3U] : defaultchars[low & 3U]);
      low >>= 2;
      flags >>= 1;
    }
    for ( ; i < 32; i++) {
      Buffer[i] = (char) ((flags & 1U) ? flagchars[high & 3U] : defaultchars[high & 3U]);
      high >>= 2;
      flags >>= 1;
    }
    for (i = startdiscard; i < enddiscard; i++) {
      gbuffer1[k++] = Buffer[i];
    }
  } else {
    genomecomp_move_absolute(this,ptr);
    high = genomecomp_read_current(this);
    low = genomecomp_read_current(this);
    flags = genomecomp_read_current(this);

    for (i = 0; i < 16; i++) {
      Buffer[i] = (char) ((flags & 1U) ? flagchars[low & 3U] : defaultchars[low & 3U]);
      low >>= 2;
      flags >>= 1;
    }
    for ( ; i < 32; i++) {
      Buffer[i] = (char) ((flags & 1U) ? flagchars[high & 3U] : defaultchars[high & 3U]);
      high >>= 2;
      flags >>= 1;
    }
    for (i = startdiscard; i < 32; i++) {
      gbuffer1[k++] = Buffer[i];
    }
    ptr += 3;
      
    while (ptr < endblock) {
      high = genomecomp_read_current(this);
      low = genomecomp_read_current(this);
      flags = genomecomp_read_current(this);

      for (i = 0; i < 16; i++) {
	gbuffer1[k++] = (char) ((flags & 1U) ? flagchars[low & 3U] : defaultchars[low & 3U]);
	low >>= 2;
	flags >>= 1;
      }
      for ( ; i < 32; i++) {
	gbuffer1[k++] = (char) ((flags & 1U) ? flagchars[high & 3U] : defaultchars[high & 3U]);
	high >>= 2;
	flags >>= 1;
      }
      ptr += 3;
    }

    if (enddiscard > 0) {
      high = genomecomp_read_current(this);
      low = genomecomp_read_current(this);
      flags = genomecomp_read_current(this);

      for (i = 0; i < 16; i++) {
	Buffer[i] = (char) ((flags & 1U) ? flagchars[low & 3U] : defaultchars[low & 3U]);
	low >>= 2;
	flags >>= 1;
      }
      for ( ; i < 32; i++) {
	Buffer[i] = (char) ((flags & 1U) ? flagchars[high & 3U] : defaultchars[high & 3U]);
	high >>= 2;
	flags >>= 1;
      }
      for (i = 0; i < enddiscard; i++) {
	gbuffer1[k++] = Buffer[i];
      }
    }
  }

  return;
}

static void
ntcounts_fileio (int *na, int *nc, int *ng, int *nt,
		 T this, Genomicpos_T startpos, Genomicpos_T endpos,
		 const char defaultchars[], const char flagchars[]) {
  /* Genomicpos_T length = endpos - startpos; */
  Genomicpos_T startblock, endblock, ptr;
  UINT4 high, low, flags;
  char Buffer[32], c;
  int startdiscard, enddiscard, i;

  /* sequence = (char *) CALLOC(length+1,sizeof(char)); */

  ptr = startblock = startpos/32U*3;
  endblock = endpos/32U*3;
  startdiscard = startpos % 32;
  enddiscard = endpos % 32;
  
  if (endblock == startblock) {
    /* Special case */
    genomecomp_move_absolute(this,ptr);
    high = genomecomp_read_current(this);
    low = genomecomp_read_current(this);
    flags = genomecomp_read_current(this);

    for (i = 0; i < 16; i++) {
      Buffer[i] = (char) ((flags & 1U) ? flagchars[low & 3U] : defaultchars[low & 3U]);
      low >>= 2;
      flags >>= 1;
    }
    for ( ; i < 32; i++) {
      Buffer[i] = (char) ((flags & 1U) ? flagchars[high & 3U] : defaultchars[high & 3U]);
      high >>= 2;
      flags >>= 1;
    }
    for (i = startdiscard; i < enddiscard; i++) {
      switch (Buffer[i]) {
      case 'A': case 'a': (*na)++; break;
      case 'C': case 'c': (*nc)++; break;
      case 'G': case 'g': (*ng)++; break;
      case 'T': case 't': (*nt)++; break;
      }
    }
  } else {
    genomecomp_move_absolute(this,ptr);
    high = genomecomp_read_current(this);
    low = genomecomp_read_current(this);
    flags = genomecomp_read_current(this);

    for (i = 0; i < 16; i++) {
      Buffer[i] = (char) ((flags & 1U) ? flagchars[low & 3U] : defaultchars[low & 3U]);
      low >>= 2;
      flags >>= 1;
    }
    for ( ; i < 32; i++) {
      Buffer[i] = (char) ((flags & 1U) ? flagchars[high & 3U] : defaultchars[high & 3U]);
      high >>= 2;
      flags >>= 1;
    }
    for (i = startdiscard; i < 32; i++) {
      switch (Buffer[i]) {
      case 'A': case 'a': (*na)++; break;
      case 'C': case 'c': (*nc)++; break;
      case 'G': case 'g': (*ng)++; break;
      case 'T': case 't': (*nt)++; break;
      }
    }
    ptr += 3;
      
    while (ptr < endblock) {
      high = genomecomp_read_current(this);
      low = genomecomp_read_current(this);
      flags = genomecomp_read_current(this);

      for (i = 0; i < 16; i++) {
	c = (char) ((flags & 1U) ? flagchars[low & 3U] : defaultchars[low & 3U]);
	switch (c) {
	case 'A': case 'a': (*na)++; break;
	case 'C': case 'c': (*nc)++; break;
	case 'G': case 'g': (*ng)++; break;
	case 'T': case 't': (*nt)++; break;
	}
	low >>= 2;
	flags >>= 1;
      }
      for ( ; i < 32; i++) {
	c = (char) ((flags & 1U) ? flagchars[high & 3U] : defaultchars[high & 3U]);
	switch (c) {
	case 'A': case 'a': (*na)++; break;
	case 'C': case 'c': (*nc)++; break;
	case 'G': case 'g': (*ng)++; break;
	case 'T': case 't': (*nt)++; break;
	}
	high >>= 2;
	flags >>= 1;
      }
      ptr += 3;
    }

    if (enddiscard > 0) {
      high = genomecomp_read_current(this);
      low = genomecomp_read_current(this);
      flags = genomecomp_read_current(this);

      for (i = 0; i < 16; i++) {
	Buffer[i] = (char) ((flags & 1U) ? flagchars[low & 3U] : defaultchars[low & 3U]);
	low >>= 2;
	flags >>= 1;
      }
      for ( ; i < 32; i++) {
	Buffer[i] = (char) ((flags & 1U) ? flagchars[high & 3U] : defaultchars[high & 3U]);
	high >>= 2;
	flags >>= 1;
      }
      for (i = 0; i < enddiscard; i++) {
	switch (Buffer[i]) {
	case 'A': case 'a': (*na)++; break;
	case 'C': case 'c': (*nc)++; break;
	case 'G': case 'g': (*ng)++; break;
	case 'T': case 't': (*nt)++; break;
	}
      }
    }
  }

  return;
}


static void
uncompress_mmap (char *gbuffer1, UINT4 *blocks, Genomicpos_T startpos, 
		 Genomicpos_T endpos, const char defaultchars[],
		 const char flagchars[]) {
  /* Genomicpos_T length = endpos - startpos; */
  Genomicpos_T startblock, endblock, ptr;
  UINT4 high, low, flags;
  char Buffer[32];
  int startdiscard, enddiscard, i, k = 0;

  /* sequence = (char *) CALLOC(length+1,sizeof(char)); */

  ptr = startblock = startpos/32U*3;
  endblock = endpos/32U*3;
  startdiscard = startpos % 32;
  enddiscard = endpos % 32;
  
  if (endblock == startblock) {
    /* Special case */
#ifdef WORDS_BIGENDIAN
    high = Bigendian_convert_uint(blocks[ptr]);
    low = Bigendian_convert_uint(blocks[ptr+1]);
    flags = Bigendian_convert_uint(blocks[ptr+2]);
#else
    high = blocks[ptr]; low = blocks[ptr+1]; flags = blocks[ptr+2];
#endif
    for (i = 0; i < 16; i++) {
      Buffer[i] = (char) ((flags & 1U) ? flagchars[low & 3U] : defaultchars[low & 3U]);
      low >>= 2;
      flags >>= 1;
    }
    for ( ; i < 32; i++) {
      Buffer[i] = (char) ((flags & 1U) ? flagchars[high & 3U] : defaultchars[high & 3U]);
      high >>= 2;
      flags >>= 1;
    }
    for (i = startdiscard; i < enddiscard; i++) {
      gbuffer1[k++] = Buffer[i];
    }
  } else {
#ifdef WORDS_BIGENDIAN
    high = Bigendian_convert_uint(blocks[ptr]);
    low = Bigendian_convert_uint(blocks[ptr+1]);
    flags = Bigendian_convert_uint(blocks[ptr+2]);
#else
    high = blocks[ptr]; low = blocks[ptr+1]; flags = blocks[ptr+2];
#endif
    for (i = 0; i < 16; i++) {
      Buffer[i] = (char) ((flags & 1U) ? flagchars[low & 3U] : defaultchars[low & 3U]);
      low >>= 2;
      flags >>= 1;
    }
    for ( ; i < 32; i++) {
      Buffer[i] = (char) ((flags & 1U) ? flagchars[high & 3U] : defaultchars[high & 3U]);
      high >>= 2;
      flags >>= 1;
    }
    for (i = startdiscard; i < 32; i++) {
      gbuffer1[k++] = Buffer[i];
    }
    ptr += 3;
      
    while (ptr < endblock) {
#ifdef WORDS_BIGENDIAN
      high = Bigendian_convert_uint(blocks[ptr]);
      low = Bigendian_convert_uint(blocks[ptr+1]);
      flags = Bigendian_convert_uint(blocks[ptr+2]);
#else
      high = blocks[ptr]; low = blocks[ptr+1]; flags = blocks[ptr+2];
#endif
      for (i = 0; i < 16; i++) {
	gbuffer1[k++] = (char) ((flags & 1U) ? flagchars[low & 3U] : defaultchars[low & 3U]);
	low >>= 2;
	flags >>= 1;
      }
      for ( ; i < 32; i++) {
	gbuffer1[k++] = (char) ((flags & 1U) ? flagchars[high & 3U] : defaultchars[high & 3U]);
	high >>= 2;
	flags >>= 1;
      }
      ptr += 3;
    }

    if (enddiscard > 0) {
#ifdef WORDS_BIGENDIAN
      high = Bigendian_convert_uint(blocks[ptr]);
      low = Bigendian_convert_uint(blocks[ptr+1]);
      flags = Bigendian_convert_uint(blocks[ptr+2]);
#else
      high = blocks[ptr]; low = blocks[ptr+1]; flags = blocks[ptr+2];
#endif
      for (i = 0; i < 16; i++) {
	Buffer[i] = (char) ((flags & 1U) ? flagchars[low & 3U] : defaultchars[low & 3U]);
	low >>= 2;
	flags >>= 1;
      }
      for ( ; i < 32; i++) {
	Buffer[i] = (char) ((flags & 1U) ? flagchars[high & 3U] : defaultchars[high & 3U]);
	high >>= 2;
	flags >>= 1;
      }
      for (i = 0; i < enddiscard; i++) {
	gbuffer1[k++] = Buffer[i];
      }
    }
  }

  return;
}


static void
Genome_ntcounts_mmap (int *na, int *nc, int *ng, int *nt, UINT4 *blocks,
		      Genomicpos_T startpos, Genomicpos_T endpos, const char defaultchars[],
		      const char flagchars[]) {
  /* Genomicpos_T length = endpos - startpos; */
  Genomicpos_T startblock, endblock, ptr;
  UINT4 high, low, flags;
  char Buffer[32], c;
  int startdiscard, enddiscard, i;

  /* *na = *nc = *ng = *nt = 0; */

  /* sequence = (char *) CALLOC(length+1,sizeof(char)); */

  ptr = startblock = startpos/32U*3;
  endblock = endpos/32U*3;
  startdiscard = startpos % 32;
  enddiscard = endpos % 32;
  
  if (endblock == startblock) {
    /* Special case */
#ifdef WORDS_BIGENDIAN
    high = Bigendian_convert_uint(blocks[ptr]);
    low = Bigendian_convert_uint(blocks[ptr+1]);
    flags = Bigendian_convert_uint(blocks[ptr+2]);
#else
    high = blocks[ptr]; low = blocks[ptr+1]; flags = blocks[ptr+2];
#endif
    for (i = 0; i < 16; i++) {
      Buffer[i] = (char) ((flags & 1U) ? flagchars[low & 3U] : defaultchars[low & 3U]);
      low >>= 2;
      flags >>= 1;
    }
    for ( ; i < 32; i++) {
      Buffer[i] = (char) ((flags & 1U) ? flagchars[high & 3U] : defaultchars[high & 3U]);
      high >>= 2;
      flags >>= 1;
    }
    for (i = startdiscard; i < enddiscard; i++) {
      switch (Buffer[i]) {
      case 'A': case 'a': (*na)++; break;
      case 'C': case 'c': (*nc)++; break;
      case 'G': case 'g': (*ng)++; break;
      case 'T': case 't': (*nt)++; break;
      }
    }
  } else {
#ifdef WORDS_BIGENDIAN
    high = Bigendian_convert_uint(blocks[ptr]);
    low = Bigendian_convert_uint(blocks[ptr+1]);
    flags = Bigendian_convert_uint(blocks[ptr+2]);
#else
    high = blocks[ptr]; low = blocks[ptr+1]; flags = blocks[ptr+2];
#endif
    for (i = 0; i < 16; i++) {
      Buffer[i] = (char) ((flags & 1U) ? flagchars[low & 3U] : defaultchars[low & 3U]);
      low >>= 2;
      flags >>= 1;
    }
    for ( ; i < 32; i++) {
      Buffer[i] = (char) ((flags & 1U) ? flagchars[high & 3U] : defaultchars[high & 3U]);
      high >>= 2;
      flags >>= 1;
    }
    for (i = startdiscard; i < 32; i++) {
      switch (Buffer[i]) {
      case 'A': case 'a': (*na)++; break;
      case 'C': case 'c': (*nc)++; break;
      case 'G': case 'g': (*ng)++; break;
      case 'T': case 't': (*nt)++; break;
      }
    }
    ptr += 3;
      
    while (ptr < endblock) {
#ifdef WORDS_BIGENDIAN
      high = Bigendian_convert_uint(blocks[ptr]);
      low = Bigendian_convert_uint(blocks[ptr+1]);
      flags = Bigendian_convert_uint(blocks[ptr+2]);
#else
      high = blocks[ptr]; low = blocks[ptr+1]; flags = blocks[ptr+2];
#endif
      for (i = 0; i < 16; i++) {
	c = (char) ((flags & 1U) ? flagchars[low & 3U] : defaultchars[low & 3U]);
	switch (c) {
	case 'A': case 'a': (*na)++; break;
	case 'C': case 'c': (*nc)++; break;
	case 'G': case 'g': (*ng)++; break;
	case 'T': case 't': (*nt)++; break;
	}
	low >>= 2;
	flags >>= 1;
      }
      for ( ; i < 32; i++) {
	c = (char) ((flags & 1U) ? flagchars[high & 3U] : defaultchars[high & 3U]);
	switch (c) {
	case 'A': case 'a': (*na)++; break;
	case 'C': case 'c': (*nc)++; break;
	case 'G': case 'g': (*ng)++; break;
	case 'T': case 't': (*nt)++; break;
	}
	high >>= 2;
	flags >>= 1;
      }
      ptr += 3;
    }

    if (enddiscard > 0) {
#ifdef WORDS_BIGENDIAN
      high = Bigendian_convert_uint(blocks[ptr]);
      low = Bigendian_convert_uint(blocks[ptr+1]);
      flags = Bigendian_convert_uint(blocks[ptr+2]);
#else
      high = blocks[ptr]; low = blocks[ptr+1]; flags = blocks[ptr+2];
#endif
      for (i = 0; i < 16; i++) {
	Buffer[i] = (char) ((flags & 1U) ? flagchars[low & 3U] : defaultchars[low & 3U]);
	low >>= 2;
	flags >>= 1;
      }
      for ( ; i < 32; i++) {
	Buffer[i] = (char) ((flags & 1U) ? flagchars[high & 3U] : defaultchars[high & 3U]);
	high >>= 2;
	flags >>= 1;
      }
      for (i = 0; i < enddiscard; i++) {
	switch (Buffer[i]) {
	case 'A': case 'a': (*na)++; break;
	case 'C': case 'c': (*nc)++; break;
	case 'G': case 'g': (*ng)++; break;
	case 'T': case 't': (*nt)++; break;
	}
      }
    }
  }

  return;
}


#ifdef DEBUG3
static bool
check_nucleotide (unsigned char nucleotide, char sequence) {
  printf("%d %c\n",nucleotide,sequence);
  switch (sequence) {
  case 'A': if (nucleotide != 0) return true; break;
  case 'C': if (nucleotide != 1) return true; break;
  case 'G': if (nucleotide != 2) return true; break;
  case 'T': if (nucleotide != 3) return true; break;
  }
  return false;
}
#endif

#if 0
static bool
uncompress_mmap_nucleotides_wstatus (unsigned char *gbuffer, UINT4 *blocks, Genomicpos_T startpos, 
				     Genomicpos_T endpos) {
  /* Genomicpos_T length = endpos - startpos; */
  Genomicpos_T startblock, endblock, startdiscard, enddiscard, ptr;
  UINT4 high, low, flags;
  int i, k = 0;
#ifdef DEBUG3
  char gbuffer_debug[1024];
#endif

  /* sequence = (char *) CALLOC(length+1,sizeof(char)); */

#ifdef DEBUG3
  uncompress_mmap(gbuffer_debug,blocks,startpos,endpos,DEFAULT_CHARS,DEFAULT_FLAGS);
#endif

  ptr = startblock = startpos/32U*3;
  endblock = endpos/32U*3;
  startdiscard = startpos % 32;
  enddiscard = endpos % 32;

  
  if (endblock == startblock) {
    /* Special case */
#ifdef WORDS_BIGENDIAN
    high = Bigendian_convert_uint(blocks[ptr]);
    low = Bigendian_convert_uint(blocks[ptr+1]);
    flags = Bigendian_convert_uint(blocks[ptr+2]);
#else
    high = blocks[ptr]; low = blocks[ptr+1]; flags = blocks[ptr+2];
#endif
    flags >>= startdiscard;
    if (startdiscard < 16) {
      low >>= (startdiscard+startdiscard);
      while (startdiscard < enddiscard && startdiscard < 16) {
	if (flags & 0x01) {
	  /* return false; */
	}
	gbuffer[k++] = (unsigned char) (low & 0x03);
#ifdef DEBUG3
	if (check_nucleotide(gbuffer[k-1],gbuffer_debug[k-1])) {
	  fprintf(stderr,"Case 1, startdiscard %d, enddiscard %d\n",startdiscard,enddiscard);
	  abort();
	}
#endif
	flags >>= 1;
	low >>= 2;
	startdiscard++;
      }
    }
    if (enddiscard >= 16) {
      startdiscard -= 16;
      enddiscard -= 16;
      high >>= (startdiscard+startdiscard);
      while (startdiscard < enddiscard) {
	if (flags & 0x01) {
	  /* return false; */
	}
	gbuffer[k++] = (unsigned char) (high & 0x03);
#ifdef DEBUG3
	if (check_nucleotide(gbuffer[k-1],gbuffer_debug[k-1])) {
	  fprintf(stderr,"Case 2\n");
	  abort();
	}
#endif
	flags >>= 1;
	high >>= 2;
	startdiscard++;
      }
    }
    return true;

  } else {
#ifdef WORDS_BIGENDIAN
    high = Bigendian_convert_uint(blocks[ptr]);
    low = Bigendian_convert_uint(blocks[ptr+1]);
    flags = Bigendian_convert_uint(blocks[ptr+2]);
#else
    high = blocks[ptr]; low = blocks[ptr+1]; flags = blocks[ptr+2];
#endif
    flags >>= startdiscard;
    if (startdiscard < 16) {
      low >>= (startdiscard+startdiscard);
      while (startdiscard < 16) {
	if (flags & 0x01) {
	  /* return false; */
	}
	gbuffer[k++] = (unsigned char) (low & 0x03);
#ifdef DEBUG3
	if (check_nucleotide(gbuffer[k-1],gbuffer_debug[k-1])) {
	  fprintf(stderr,"Case 3\n");
	  abort();
	}
#endif
	flags >>= 1;
	low >>= 2;
	startdiscard++;
      }
    }
    startdiscard -= 16;
    high >>= (startdiscard+startdiscard);
    while (startdiscard < 16) {
      if (flags & 0x01) {
	/* return false; */
      }
      gbuffer[k++] = (unsigned char) (high & 0x03);
#ifdef DEBUG3
      if (check_nucleotide(gbuffer[k-1],gbuffer_debug[k-1])) {
	fprintf(stderr,"Case 4, startdiscard %d (after subtracting 16)\n",startdiscard);
	abort();
      }
#endif
      flags >>= 1;
      high >>= 2;
      startdiscard++;
    }
    ptr += 3;
      
    while (ptr < endblock) {
#ifdef WORDS_BIGENDIAN
      high = Bigendian_convert_uint(blocks[ptr]);
      low = Bigendian_convert_uint(blocks[ptr+1]);
      flags = Bigendian_convert_uint(blocks[ptr+2]);
#else
      high = blocks[ptr]; low = blocks[ptr+1]; flags = blocks[ptr+2];
#endif
      if (flags) {
	/* return false; */
      }
      for (i = 0; i < 16; i++) {
	gbuffer[k++] = (unsigned char) (low & 0x03);
#ifdef DEBUG3
	if (check_nucleotide(gbuffer[k-1],gbuffer_debug[k-1])) {
	  fprintf(stderr,"Case 5\n");
	  abort();
	}
#endif
	low >>= 2;
      }
      for ( ; i < 32; i++) {
	gbuffer[k++] = (unsigned char) (high & 0x03);
#ifdef DEBUG3
	if (check_nucleotide(gbuffer[k-1],gbuffer_debug[k-1])) {
	  fprintf(stderr,"Case 6\n");
	  abort();
	}
#endif
	high >>= 2;
      }

      ptr += 3;
    }

    if (enddiscard > 0) {
#ifdef WORDS_BIGENDIAN
      high = Bigendian_convert_uint(blocks[ptr]);
      low = Bigendian_convert_uint(blocks[ptr+1]);
      flags = Bigendian_convert_uint(blocks[ptr+2]);
#else
      high = blocks[ptr]; low = blocks[ptr+1]; flags = blocks[ptr+2];
#endif
      i = 0;
      while (i < enddiscard && i < 16) {
	if (flags & 0x01) {
	  /* return false; */
	}
	gbuffer[k++] = (unsigned char) (low & 0x03);
#ifdef DEBUG3
	if (check_nucleotide(gbuffer[k-1],gbuffer_debug[k-1])) {
	  fprintf(stderr,"Case 7\n");
	  abort();
	}
#endif
	flags >>= 1;
	low >>= 2;
	i++;
      }
      while (i < enddiscard) {
	if (flags & 0x01) {
	  /* return false; */
	}
	gbuffer[k++] = (unsigned char) (high & 0x03);
#ifdef DEBUG3
	if (check_nucleotide(gbuffer[k-1],gbuffer_debug[k-1])) {
	  fprintf(stderr,"Case 8\n");
	  abort();
	}
#endif
	flags >>= 1;
	high >>= 2;
	i++;
      }
    }
  }

  return true;
}
#endif


static void
uncompress_mmap_nucleotides (unsigned char *gbuffer, UINT4 *blocks, Genomicpos_T startpos, 
			     Genomicpos_T endpos) {
  /* Genomicpos_T length = endpos - startpos; */
  Genomicpos_T startblock, endblock, ptr;
  UINT4 high, low;
  int startdiscard, enddiscard, i, k = 0;

  ptr = startblock = startpos/32U*3;
  endblock = endpos/32U*3;
  startdiscard = startpos % 32;
  enddiscard = endpos % 32;

  /* printf("startdiscard = %d, enddiscard = %d\n",startdiscard,enddiscard); */
  if (endblock == startblock) {
    /* Special case */
    if (startdiscard < 16) {
#ifdef WORDS_BIGENDIAN
      low = Bigendian_convert_uint(blocks[ptr+1]);
#else
      low = blocks[ptr+1];
#endif
      low >>= (startdiscard+startdiscard);
      while (startdiscard < enddiscard && startdiscard < 16) {
	gbuffer[k++] = (unsigned char) (low & 0x03);
	low >>= 2;
	startdiscard++;
      }
    }
    if (enddiscard >= 16) {
      startdiscard -= 16;
      enddiscard -= 16;
#ifdef WORDS_BIGENDIAN
      high = Bigendian_convert_uint(blocks[ptr]);
#else
      high = blocks[ptr];
#endif
      high >>= (startdiscard+startdiscard);
      while (startdiscard < enddiscard) {
	gbuffer[k++] = (unsigned char) (high & 0x03);
	high >>= 2;
	startdiscard++;
      }
    }

  } else {
#ifdef WORDS_BIGENDIAN
    high = Bigendian_convert_uint(blocks[ptr]);
    low = Bigendian_convert_uint(blocks[ptr+1]);
#else
    high = blocks[ptr]; low = blocks[ptr+1];
#endif
    if (startdiscard < 16) {
      low >>= (startdiscard+startdiscard);
      while (startdiscard < 16) {
	gbuffer[k++] = (unsigned char) (low & 0x03);
	low >>= 2;
	startdiscard++;
      }
    }
    startdiscard -= 16;
    high >>= (startdiscard+startdiscard);
    while (startdiscard < 16) {
      gbuffer[k++] = (unsigned char) (high & 0x03);
      high >>= 2;
      startdiscard++;
    }
    /* printf("Block %d assigned up to %d\n",ptr,k); */
    ptr += 3;
      
    while (ptr < endblock) {
#ifdef WORDS_BIGENDIAN
      high = Bigendian_convert_uint(blocks[ptr]);
      low = Bigendian_convert_uint(blocks[ptr+1]);
#else
      high = blocks[ptr]; low = blocks[ptr+1];
#endif
      for (i = 0; i < 16; i++) {
	gbuffer[k++] = (unsigned char) (low & 0x03);
	low >>= 2;
      }
      for ( ; i < 32; i++) {
	gbuffer[k++] = (unsigned char) (high & 0x03);
	high >>= 2;
      }

      /* printf("Block %d assigned up to %d\n",ptr,k); */
      ptr += 3;
    }

#ifdef WORDS_BIGENDIAN
    high = Bigendian_convert_uint(blocks[ptr]);
    low = Bigendian_convert_uint(blocks[ptr+1]);
#else
    high = blocks[ptr]; low = blocks[ptr+1];
#endif
    i = 0;
    while (i < enddiscard && i < 16) {
      gbuffer[k++] = (unsigned char) (low & 0x03);
      low >>= 2;
      i++;
    }
    while (i < enddiscard) {
      gbuffer[k++] = (unsigned char) (high & 0x03);
      high >>= 2;
      i++;
    }
    /* printf("Block %d assigned up to %d\n",ptr,k); */
  }

#if 0
  for (ptr = startpos, k = 0; ptr < endpos; ptr++, k++) {
    if (gbuffer[k] > 3) {
      printf("startpos = %u, endpos = %u, k = %d\n",startpos,endpos,k);
      printf("startblock %u, endblock %u\n",startblock,endblock);
      for (ptr = startpos, k = 0; ptr < endpos; ptr++, k++) {
	printf("%u %d %d\n",ptr,k,gbuffer[k]);
      }
      abort();
    }
  }
#endif

  return;
}



static const Except_T gbufferlen_error = { "Insufficient allocation" };

static bool
fill_buffer (Chrnum_T *chrnum, int *nunknowns, T this, Genomicpos_T left, Genomicpos_T length, char *gbuffer1,
	     IIT_T chromosome_iit, const char defaultchars[], const char flagchars[]) {
  Chrnum_T chrnum_left, chrnum_right, chrnumi;
  Genomicpos_T inbounds_low, inbounds_high, pos, low, high, maxoverlap, overlap;
  
  *nunknowns = 0;
  if (length == 0) {
    *chrnum = 0;
    return false;
  }

  /* Fix out of bounds resulting from negative numbers */
  if (left + length < left) {
    debug(printf("Got negative left\n"));
    while (left != 0U) {
      *(gbuffer1++) = OUTOFBOUNDS;
      *nunknowns += 1;
      length--;
      left++;
    }
  }

  if (this->compressedp == false) {
    if (this->access == FILEIO) {
#ifdef HAVE_PTHREAD
      pthread_mutex_lock(&this->read_mutex);
#endif
      if (lseek(this->fd,left,SEEK_SET) < 0) {
	perror("Error in fill_buffer");
	exit(9);
      }
      read(this->fd,gbuffer1,length);
#ifdef HAVE_PTHREAD
      pthread_mutex_unlock(&this->read_mutex);
#endif

    } else {
      memcpy(gbuffer1,&(this->chars[left]),length*sizeof(char));
    }

  } else {
    if (this->access == FILEIO) {
#ifdef HAVE_PTHREAD
      pthread_mutex_lock(&this->read_mutex);
#endif
      uncompress_fileio(gbuffer1,this,left,left+length,defaultchars,flagchars);
#ifdef HAVE_PTHREAD
      pthread_mutex_unlock(&this->read_mutex);
#endif
    } else {
      uncompress_mmap(gbuffer1,this->blocks,left,left+length,defaultchars,flagchars);
    }
  }
  gbuffer1[length] = '\0';

  debug(printf("Got sequence at %u with length %u, forward\n",left,length));

  /* Fix out of bounds resulting from crossing chromosomes */
  if (chromosome_iit == NULL) {
    debug(printf("No chr info because chromosome_iit is null\n"));
    /* None provided, perhaps because aligning to user segment */
  } else if ((chrnum_left = IIT_get_one(chromosome_iit,/*divstring*/NULL,left,left)) == 
	     (chrnum_right = IIT_get_one(chromosome_iit,/*divstring*/NULL,left+length-1U,left+length-1U))) {
    debug(printf("Chr at beginning = %d, at end = %d\n",chrnum_left,chrnum_right));
    *chrnum = chrnum_left;

    /* Fix out of bounds resulting from going past last chromosome */
    IIT_interval_bounds(&low,&high,chromosome_iit,chrnum_left);
    if (left > high) {
      for (pos = 0; pos < length; pos++) {
	gbuffer1[pos] = OUTOFBOUNDS;
	*nunknowns += 1;
      }
    } else if (left+length > high) {
      for (pos = high-left+1; pos < length; pos++) {
	gbuffer1[pos] = OUTOFBOUNDS;
	*nunknowns += 1;
      }
    }

  } else {
    debug(printf("Chr at beginning = %d, at end = %d\n",chrnum_left,chrnum_right));
    maxoverlap = 0;
    for (chrnumi = chrnum_left; chrnumi <= chrnum_right; chrnumi++) {
      IIT_interval_bounds(&low,&high,chromosome_iit,chrnumi);
      if (left > low) {
	low = left;
      }
      if (left+length < high) {
	high = left+length-1U;
      }
      if ((overlap = high - low) > maxoverlap) {
	*chrnum = chrnumi;
	maxoverlap = overlap;
	inbounds_low = low - left;
	inbounds_high = high - left;
      }
    }
    debug(printf("in-bounds at %d..%d\n",inbounds_low,inbounds_high));
    for (pos = 0; pos < inbounds_low; pos++) {
      gbuffer1[pos] = OUTOFBOUNDS;
      *nunknowns += 1;
    }
    for (pos = inbounds_high+1; pos < length; pos++) {
      gbuffer1[pos] = OUTOFBOUNDS;
      *nunknowns += 1;
    }
    debug(printf("%s\n",gbuffer1));
  }

  return true;
}



/************************************************************************
 *   Usage procedures
 ************************************************************************/

static T genome;
static UINT4 *genome_blocks;

void
Genome_setup (T genome_in) {
  genome = genome_in;
  genome_blocks = genome->blocks;
  return;
}


bool
Genome_fill_buffer (Chrnum_T *chrnum, int *nunknowns, T this, Genomicpos_T left, Genomicpos_T length, char *gbuffer1,
		    IIT_T chromosome_iit) {
  return fill_buffer(&(*chrnum),&(*nunknowns),this,left,length,gbuffer1,chromosome_iit,DEFAULT_CHARS,DEFAULT_FLAGS);
}

#if 0
/* See get_segment_alt and get_segment_cmp below */
bool
Genome_fill_buffer_alt (Chrnum_T *chrnum, int *nunknowns, T this, Genomicpos_T left, Genomicpos_T length, char *gbuffer1,
			IIT_T chromosome_iit) {
  return fill_buffer(&(*chrnum),&(*nunknowns),this,left,length,gbuffer1,chromosome_iit,DEFAULTCHARS,DEFAULTCHARS);
}

bool
Genome_fill_buffer_snp (Chrnum_T *chrnum, int *nunknowns, T this, Genomicpos_T left, Genomicpos_T length, char *gbuffer1,
			IIT_T chromosome_iit) {
  return fill_buffer(&(*chrnum),&(*nunknowns),this,left,length,gbuffer1,chromosome_iit,SNP_DEFAULTCHARS,SNP_FLAGCHARS);
}
#endif


void
Genome_fill_buffer_simple (T this, Genomicpos_T left, Genomicpos_T length, char *gbuffer1) {
  int delta, i;
  
#if 0
  assert(left + length >= left);
#endif

  /* Fix out of bounds resulting from negative numbers */
  if (left + length < left) {
    fprintf(stderr,"left %u + length %u < left %u\n",left,length,left);
    delta = -left;
    length -= delta;
    for (i = 0; i < delta; i++) {
      gbuffer1[i] = 'N';
    }
    gbuffer1[i] = '\0';
    gbuffer1 += delta;
    left = 0U;
  }

  if (length == 0) {
    return;
  }


  if (this->compressedp == false) {
    if (this->access == FILEIO) {
#ifdef HAVE_PTHREAD
      pthread_mutex_lock(&this->read_mutex);
#endif
      if (lseek(this->fd,left,SEEK_SET) < 0) {
	perror("Error in Genome_fill_buffer_simple");
	exit(9);
      }
      read(this->fd,gbuffer1,length);
#ifdef HAVE_PTHREAD
      pthread_mutex_unlock(&this->read_mutex);
#endif

    } else {
      memcpy(gbuffer1,&(this->chars[left]),length*sizeof(char));
    }

  } else {
    if (this->access == FILEIO) {
#ifdef HAVE_PTHREAD
      pthread_mutex_lock(&this->read_mutex);
#endif
      uncompress_fileio(gbuffer1,this,left,left+length,global_chars,global_flags);
#ifdef HAVE_PTHREAD
      pthread_mutex_unlock(&this->read_mutex);
#endif
    } else {
      uncompress_mmap(gbuffer1,this->blocks,left,left+length,global_chars,global_flags);
    }
  }
  gbuffer1[length] = '\0';

  debug(printf("Got sequence at %u with length %u, forward\n",left,length));

  return;
}


void
Genome_fill_buffer_blocks (Genomicpos_T left, Genomicpos_T length, char *gbuffer1) {
  
  if (length > 0) {
    assert(left + length >= left);
    uncompress_mmap(gbuffer1,genome_blocks,left,left+length,global_chars,global_flags);
  }
  gbuffer1[length] = '\0';
  return;
}

void
Genome_fill_buffer_blocks_noterm (Genomicpos_T left, Genomicpos_T length, char *gbuffer1) {
  
  if (length > 0) {
    assert(left + length >= left);
    uncompress_mmap(gbuffer1,genome_blocks,left,left+length,global_chars,global_flags);
  }
  /* gbuffer1[length] = '\0'; */
  return;
}


void
Genome_fill_buffer_simple_alt (T this, Genomicpos_T left, Genomicpos_T length, char *gbuffer1) {
  int delta, i;
  
#if 0
  assert(left + length >= left);
#endif

  /* Fix out of bounds resulting from negative numbers */
  if (left + length < left) {
    fprintf(stderr,"left %u + length %u < left %u\n",left,length,left);
    delta = -left;
    length -= delta;
    for (i = 0; i < delta; i++) {
      gbuffer1[i] = 'N';
    }
    gbuffer1[i] = '\0';
    gbuffer1 += delta;
    left = 0U;
  }

  if (length == 0) {
    return;
  }


  if (this->compressedp == false) {
    if (this->access == FILEIO) {
#ifdef HAVE_PTHREAD
      pthread_mutex_lock(&this->read_mutex);
#endif
      if (lseek(this->fd,left,SEEK_SET) < 0) {
	perror("Error in gmap, Genome_get_segment");
	exit(9);
      }
      read(this->fd,gbuffer1,length);
#ifdef HAVE_PTHREAD
      pthread_mutex_unlock(&this->read_mutex);
#endif

    } else {
      memcpy(gbuffer1,&(this->chars[left]),length*sizeof(char));
    }

  } else {
    if (this->access == FILEIO) {
#ifdef HAVE_PTHREAD
      pthread_mutex_lock(&this->read_mutex);
#endif
      uncompress_fileio(gbuffer1,this,left,left+length,DEFAULT_CHARS,SNP_FLAGS);
#ifdef HAVE_PTHREAD
      pthread_mutex_unlock(&this->read_mutex);
#endif
    } else {
      uncompress_mmap(gbuffer1,this->blocks,left,left+length,DEFAULT_CHARS,SNP_FLAGS);
    }
  }
  gbuffer1[length] = '\0';

  debug(printf("Got sequence at %u with length %u, forward\n",left,length));

  return;
}


void
Genome_fill_buffer_nucleotides (T this, Genomicpos_T left, Genomicpos_T length, unsigned char *gbuffer) {
  
  if (length == 0) {
    return;
  }

  /* Fix out of bounds resulting from negative numbers */
#if 0
  if (left + length < left) {
    fprintf(stderr,"left %u + length %u < left %u\n",left,length,left);
    abort();
  }
#else
  assert(left + length >= left);
#endif

  if (this->compressedp == false) {
    fprintf(stderr,"Procedure Genome_fill_buffer_nucleotides not designed to work for non-compressed genomes\n");
    exit(9);
    if (this->access == FILEIO) {
#ifdef HAVE_PTHREAD
      pthread_mutex_lock(&this->read_mutex);
#endif
#ifdef HAVE_PTHREAD
      pthread_mutex_unlock(&this->read_mutex);
#endif

    } else {
      memcpy(gbuffer,&(this->chars[left]),length*sizeof(char));
    }

  } else {
    if (this->access == FILEIO) {
#ifdef HAVE_PTHREAD
      pthread_mutex_lock(&this->read_mutex);
#endif
      fprintf(stderr,"Procedure Genome_fill_buffer_nucleotides not designed to work under FILEIO access\n");
      exit(9);
#ifdef HAVE_PTHREAD
      pthread_mutex_unlock(&this->read_mutex);
#endif
    } else {
      uncompress_mmap_nucleotides(gbuffer,this->blocks,left,left+length);
    }
  }

  gbuffer[length] = 0xFF;

#if 0
  for (i = 0; i <= length; i++) {
    printf("%d ",gbuffer[i]);
  }
  printf("\n");
#endif

  return;
}


char
Genome_get_char (T this, Genomicpos_T left) {
  char gbuffer1[1];
  
  if (this->compressedp == false) {
    if (this->access == FILEIO) {
#ifdef HAVE_PTHREAD
      pthread_mutex_lock(&this->read_mutex);
#endif
      if (lseek(this->fd,left,SEEK_SET) < 0) {
	perror("Error in Genome_get_char");
	exit(9);
      }
      read(this->fd,gbuffer1,/*length*/1);
#ifdef HAVE_PTHREAD
      pthread_mutex_unlock(&this->read_mutex);
#endif

    } else {
      memcpy(gbuffer1,&(this->chars[left]),sizeof(char));
    }

  } else {
    if (this->access == FILEIO) {
#ifdef HAVE_PTHREAD
      pthread_mutex_lock(&this->read_mutex);
#endif
      uncompress_fileio(gbuffer1,this,left,left+1,global_chars,global_flags);
#ifdef HAVE_PTHREAD
      pthread_mutex_unlock(&this->read_mutex);
#endif
    } else {
      uncompress_mmap(gbuffer1,this->blocks,left,left+1,global_chars,global_flags);
    }
  }

  return gbuffer1[0];
}


char
Genome_get_char_blocks (Genomicpos_T left) {
  char gbuffer1[1];
  
  uncompress_mmap(gbuffer1,genome_blocks,left,left+1,global_chars,global_flags);
  return gbuffer1[0];
}


bool
buffers_diff_p (char *buffer1, char *buffer2, int length) {
  int i;

  for (i = 0; i < length; i++) {
    if (buffer1[i] != buffer2[i]) {
      return true;
    }
  }
  return false;
}


Sequence_T
Genome_get_segment (T this, Genomicpos_T left, Genomicpos_T length, IIT_T chromosome_iit,
		    bool revcomp) {
  Chrnum_T chrnum;
  int nunknowns;
  char *gbuffer;
  
  gbuffer = (char *) CALLOC(length+1,sizeof(char));

  fill_buffer(&chrnum,&nunknowns,this,left,length,gbuffer,chromosome_iit,DEFAULT_CHARS,DEFAULT_FLAGS);

  if (revcomp == true) {
    /* make_complement_buffered(gbuffer2,gbuffer1,length);*/
    make_complement_inplace(gbuffer,length);
    debug(printf("Got sequence at %u with length %u, revcomp\n",left,length));
    debug1(Sequence_print(stdout,Sequence_genomic_new(gbuffer,length,/*copyp*/false),false,60,true));
    return Sequence_genomic_new(gbuffer,length,/*copyp*/false);
  } else {
    debug(printf("Got sequence at %u with length %u, forward\n",left,length));
    debug1(Sequence_print(stdout,Sequence_genomic_new(gbuffer,length,/*copyp*/false),false,60,true));
    return Sequence_genomic_new(gbuffer,length,/*copyp*/false);
  }


}

Sequence_T
Genome_get_segment_alt (T this, Genomicpos_T left, Genomicpos_T length, IIT_T chromosome_iit,
			bool revcomp) {
  Chrnum_T chrnum;
  int nunknowns;
  char *gbuffer;
  
  gbuffer = (char *) CALLOC(length+1,sizeof(char));
  
  fill_buffer(&chrnum,&nunknowns,this,left,length,gbuffer,chromosome_iit,DEFAULT_CHARS,SNP_FLAGS);

  if (revcomp == true) {
    /* make_complement_buffered(gbuffer2,gbuffer1,length); */
    make_complement_inplace(gbuffer,length);
    debug(printf("Got sequence at %u with length %u, revcomp\n",left,length));
    debug1(Sequence_print(stdout,Sequence_genomic_new(gbuffer,length,/*copyp*/false),false,60,true));
    return Sequence_genomic_new(gbuffer,length,/*copyp*/false);
  } else {
    debug(printf("Got sequence at %u with length %u, forward\n",left,length));
    debug1(Sequence_print(stdout,Sequence_genomic_new(gbuffer1,length,/*copyp*/false),false,60,true));
    return Sequence_genomic_new(gbuffer,length,/*copyp*/false);
  }
}

Sequence_T
Genome_get_segment_snp (T this, Genomicpos_T left, Genomicpos_T length, IIT_T chromosome_iit,
			bool revcomp) {
  Chrnum_T chrnum;
  int nunknowns;
  char *gbuffer;
  
  gbuffer = (char *) CALLOC(length+1,sizeof(char));

  fill_buffer(&chrnum,&nunknowns,this,left,length,gbuffer,chromosome_iit,SNP_CHARS,SNP_FLAGS);

  if (revcomp == true) {
    /* make_complement_buffered(gbuffer2,gbuffer1,length); */
    make_complement_inplace(gbuffer,length);
    debug(printf("Got sequence at %u with length %u, revcomp\n",left,length));
    debug1(Sequence_print(stdout,Sequence_genomic_new(gbuffer,length,/*copyp*/false),false,60,true));
    return Sequence_genomic_new(gbuffer,length,/*copyp*/false);
  } else {
    debug(printf("Got sequence at %u with length %u, forward\n",left,length));
    debug1(Sequence_print(stdout,Sequence_genomic_new(gbuffer,length,/*copyp*/false),false,60,true));
    return Sequence_genomic_new(gbuffer,length,/*copyp*/false);
  }
}


int
Genome_ntcounts (int *na, int *nc, int *ng, int *nt,
		 T this, Genomicpos_T left, Genomicpos_T length) {
  char *gbuffer, *p;
  unsigned int i;
  
  *na = *nc = *ng = *nt = 0;

  if (length == 0) {
    return 0;
  }

  /* Fix out of bounds resulting from negative numbers */
#if 0
  if (left + length < left) {
    fprintf(stderr,"left %u + length %u < left %u\n",left,length,left);
    abort();
  }
#else
  assert(left + length >= left);
#endif

  if (this->compressedp == false) {
    if (this->access == FILEIO) {
      gbuffer = (char *) CALLOC(length,sizeof(char));
#ifdef HAVE_PTHREAD
      pthread_mutex_lock(&this->read_mutex);
#endif
      if (lseek(this->fd,left,SEEK_SET) < 0) {
	perror("Error in gmap, Genome_get_segment");
	exit(9);
      }
      read(this->fd,gbuffer,length);
#ifdef HAVE_PTHREAD
      pthread_mutex_unlock(&this->read_mutex);
#endif
      p = &(gbuffer[0]);
      for (i = 0; i < length; i++) {
	switch (*p) {
	case 'A': case 'a': (*na)++; break;
	case 'C': case 'c': (*nc)++; break;
	case 'G': case 'g': (*ng)++; break;
	case 'T': case 't': (*nt)++; break;
	}
	p++;
      }
      FREE(gbuffer);

    } else {
      p = &(this->chars[left]);
      for (i = 0; i < length; i++) {
	switch (*p) {
	case 'A': case 'a': (*na)++; break;
	case 'C': case 'c': (*nc)++; break;
	case 'G': case 'g': (*ng)++; break;
	case 'T': case 't': (*nt)++; break;
	}
	p++;
      }
    }

  } else {
    if (this->access == FILEIO) {
#ifdef HAVE_PTHREAD
      pthread_mutex_lock(&this->read_mutex);
#endif
      ntcounts_fileio(&(*na),&(*nc),&(*ng),&(*nt),
		      this,left,left+length,global_chars,global_flags);
#ifdef HAVE_PTHREAD
      pthread_mutex_unlock(&this->read_mutex);
#endif
    } else {
      Genome_ntcounts_mmap(&(*na),&(*nc),&(*ng),&(*nt),
			   this->blocks,left,left+length,global_chars,global_flags);
    }
  }

  debug(printf("Got sequence at %u with length %u, forward\n",left,length));

  return (*na) + (*nc) + (*ng) + (*nt);
}




#if 0
int
Genome_next_char (T this) {
  char gbuffer[2];
  
  if (this->compressedp == false) {
    if (this->access == FILEIO) {
#ifdef HAVE_PTHREAD
      pthread_mutex_lock(&this->read_mutex);
#endif
      read(this->fd,gbuffer,1);
#ifdef HAVE_PTHREAD
      pthread_mutex_unlock(&this->read_mutex);
#endif
      return (int) gbuffer[0];

    } else {
      return *(this->ptr++);
    }

  } else {
    if (this->access == FILEIO) {
#ifdef HAVE_PTHREAD
      pthread_mutex_lock(&this->read_mutex);
#endif
      uncompress_fileio(gbuffer,this,this->left,this->left+1U,DEFAULT_CHARS,DEFAULT_FLAGS);
#ifdef HAVE_PTHREAD
      pthread_mutex_unlock(&this->read_mutex);
#endif
    } else {
      uncompress_mmap(gbuffer,this->blocks,this->left,this->left+1U,DEFAULT_CHARS,DEFAULT_FLAGS);
    }
    this->left += 1U;
    return (int) gbuffer[0];
  }
}
#endif



/* sbuffer is the strain buffer, with coordinates low--high.  gbuffer1
   has the reference sequence in the forward direction.  The reference
   sequence has coordinates left + length.  The calling procedure has
   ordered the indices in descending genomic position, so that right
   shifting works. */
Sequence_T
Genome_patch_strain (int *indices, int nindices, IIT_T altstrain_iit, 
		     Genomicpos_T refL, Genomicpos_T reflen,
		     bool revcomp, char *gbuffer1, char *gbuffer2, char *gbuffer3,
		     int gbuffer3len) {
  Genomicpos_T refR, srcL, srcR, matR;
  Interval_T interval;
  char *dest, *src, *matbuffer, *shiftdest, *shiftsrc, *restofheader;
  int index, i, matlen, patchlen, shiftlen, expansion;
  bool allocp;
  
  assert(reflen <= gbuffer3len);
  refR = refL + reflen;
  debug2(printf("refL=%u refR=%u reflen=%u\n",refL,refR,reflen));

  /* Work in gbuffer3 */
  memcpy(gbuffer3,gbuffer1,reflen);

  for (i = 0; i < nindices; i++) {
    index = indices[i];
    interval = IIT_interval(altstrain_iit,index);
    srcL = Interval_low(interval);
    srcR = Interval_high(interval) + 1;	/* Intervals are inclusive */
    matbuffer = IIT_annotation(&restofheader,altstrain_iit,index,&allocp); /* Holds the sequence */
    matlen = IIT_annotation_strlen(altstrain_iit,index);
    matR = srcL + matlen;

    /* Truncate srcR and matR */
    if (srcR > refR) {
      srcR = refR;
    }
    if (matR > refR) {
      matR = refR;
    }
    expansion = matR - srcR;

    /* Find dest and src */
    if (srcL < refL) {
      dest = &(gbuffer3[0]);
      src = &(matbuffer[refL-srcL]);
      patchlen = matR - refL;
    } else {
      dest = &(gbuffer3[srcL-refL]);
      src = &(matbuffer[0]);
      patchlen = matR - srcL;
    }
    if (allocp == true) {
      FREE(restofheader);
    }
    debug2(printf("srcL=%u matR=%u srcR=%u matlen=%u patchlen=%d expansion=%d\n",
		  srcL,matR,srcR,matlen,patchlen,expansion));

    /* If patchlen < 0, then matR is to left of refL and we are done */
    if (patchlen > 0) {
      memcpy(dest,src,patchlen*sizeof(char));
      if (expansion < 0) {
	/* Contraction: shouldn't occur because gmapindex will fill in with x's */
	dest += patchlen;
	while (expansion < 0) {
	  *dest++ = 'x';
	  expansion++;
	}
      } else if (expansion > 0) {
	dest += patchlen;
	src += patchlen;
	shiftlen = refR - matR - expansion;
	shiftsrc = dest;
	shiftdest = shiftsrc + expansion;
	memmove(shiftdest,shiftsrc,shiftlen*sizeof(char));
	memcpy(dest,src,expansion*sizeof(char));
	debug2(printf("  shifted %d\n",shiftlen));
      }
    }
  }
  debug2(printf("\n"));

  if (revcomp == true) {
    make_complement_buffered(gbuffer2,gbuffer3,reflen);
    /* FREE(sequence); */
    debug(printf("Got sequence at %u with length %u, revcomp\n",refL,reflen));
    debug1(Sequence_print(stdout,Sequence_genomic_new(gbuffer2,reflen,/*copyp*/false),false,60,true));
    return Sequence_genomic_new(gbuffer2,reflen,/*copyp*/false);
  } else {
    debug(printf("Got sequence at %u with length %u, forward\n",refL,reflen));
    debug1(Sequence_print(stdout,Sequence_genomic_new(gbuffer3,reflen,/*copyp*/false),false,60,true));
    return Sequence_genomic_new(gbuffer3,reflen,/*copyp*/false);
  }
}


