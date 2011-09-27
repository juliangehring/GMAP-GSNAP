static char rcsid[] = "$Id: genome.c,v 1.81 2005/10/21 16:42:55 twu Exp $";
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
#include "mem.h"
#include "types.h"
#include "access.h"
#include "complement.h"
#include "interval.h"

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

/* Patching with strain information */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif



#define T Genome_T
struct T {
  Access_T access;
  int fd;
  size_t len;

  char *chars;
  UINT4 *blocks;
  bool compressedp;

#ifdef HAVE_PTHREAD
  pthread_mutex_t read_mutex;
#endif
};


void
Genome_free (T *old) {
  if (*old) {
    if ((*old)->access == ALLOCATED) {
      abort();
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
Genome_new (char *genomesubdir, char *fileroot, bool uncompressedp, bool batchp) {
  T new = (T) MALLOC(sizeof(*new));
  char *filename;
  bool compressedp = !uncompressedp;
  int npages;
  double seconds;

  new->compressedp = compressedp;

  if (compressedp == true) {
    filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+strlen(fileroot)+
			       strlen(".genomecomp")+1,sizeof(char));
    sprintf(filename,"%s/%s.genomecomp",genomesubdir,fileroot);
  } else {
    filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+strlen(fileroot)+
			       strlen(".genome")+1,sizeof(char));
    sprintf(filename,"%s/%s.genome",genomesubdir,fileroot);
  }

  if (compressedp == true) {
    new->chars = (char *) NULL;
#ifndef HAVE_MMAP
    new->blocks = (UINT4 *) NULL;
    new->fd = Access_fileio(filename);
    new->access = FILEIO;
#else
    if (batchp == true) {
      fprintf(stderr,"Pre-loading uncompressed genome...");
      new->blocks = (UINT4 *) Access_mmap_and_preload(&new->fd,&new->len,&npages,&seconds,
						      filename,sizeof(UINT4));
      if (new->blocks == NULL) {
	fprintf(stderr,"insufficient memory (will use disk file instead)\n");
	new->access = FILEIO;
      } else {
	fprintf(stderr,"done (%lu bytes, %d pages, %.2f sec)\n",
		new->len,npages,seconds);
	new->access = MMAPPED;
      }
    } else {
      new->blocks = (UINT4 *) Access_mmap(&new->fd,&new->len,filename,sizeof(UINT4),/*randomp*/false);
      if (new->blocks == NULL) {
	new->access = FILEIO;
      } else {
	new->access = MMAPPED;
      }
    }
#endif /* HAVE_MMAP */

  } else {
    new->blocks = (UINT4 *) NULL;

#ifndef HAVE_MMAP
    new->chars = (char *) NULL;
    new->fd = Access_fileio(filename);
    new->access = FILEIO;
#else
    if (batchp == true) {
      fprintf(stderr,"Pre-loading compressed genome...");
      new->chars = (char *) Access_mmap_and_preload(&new->fd,&new->len,&npages,&seconds,
						    filename,sizeof(char));
      if (new->chars == NULL) {
	fprintf(stderr,"insufficient memory (will use disk file instead)\n");
	new->access = FILEIO;
      } else {
	fprintf(stderr,"done (%lu bytes, %d pages, %.2f sec)\n",
		new->len,npages,seconds);
	new->access = MMAPPED;
      }
    } else {
      new->chars = (char *) Access_mmap(&new->fd,&new->len,filename,sizeof(char),/*randomp*/false);
      if (new->chars == NULL) {
	new->access = FILEIO;
      } else {
	new->access = MMAPPED;
      }
    }
#endif /* HAVE_MMAP */
  }

#ifdef HAVE_PTHREAD
  if (new->access == FILEIO) {
    pthread_mutex_init(&new->read_mutex,NULL);
  }
#endif

  FREE(filename);

  return new;
}


static void
make_complement_buffered (char *complement, char *sequence, Genomicpos_T length) {
  char complCode[128] = COMPLEMENT;
  int i, j;

  for (i = length-1, j = 0; i >= 0; i--, j++) {
    complement[j] = complCode[(int) sequence[i]];
  }
  complement[length] = '\0';
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


static char translate[8] = {'A','C','G','T','N','?','?','X'};

void
Genome_replace_x (void) {
  translate[7] = 'N';
}


static void
uncompress_fileio (char *gbuffer1, T this, Genomicpos_T startpos, 
		   Genomicpos_T endpos) {
  /* Genomicpos_T length = endpos - startpos; */
  Genomicpos_T startblock, endblock, startdiscard, enddiscard, ptr;
  UINT4 high, low, flags;
  char Buffer[32];
  int i, k = 0;

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
      Buffer[i] = translate[(flags & 1U) ? (low & 3U) + 4 : low & 3U];
      low >>= 2;
      flags >>= 1;
    }
    for ( ; i < 32; i++) {
      Buffer[i] = translate[(flags & 1U) ? (high & 3U) + 4 : high & 3U];
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
      Buffer[i] = translate[(flags & 1U) ? (low & 3U) + 4 : low & 3U];
      low >>= 2;
      flags >>= 1;
    }
    for ( ; i < 32; i++) {
      Buffer[i] = translate[(flags & 1U) ? (high & 3U) + 4 : high & 3U];
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
	gbuffer1[k++] = translate[(flags & 1U) ? (low & 3U) + 4 : low & 3U];
	low >>= 2;
	flags >>= 1;
      }
      for ( ; i < 32; i++) {
	gbuffer1[k++] = translate[(flags & 1U) ? (high & 3U) + 4 : high & 3U];
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
	Buffer[i] = translate[(flags & 1U) ? (low & 3U) + 4 : low & 3U];
	low >>= 2;
	flags >>= 1;
      }
      for ( ; i < 32; i++) {
	Buffer[i] = translate[(flags & 1U) ? (high & 3U) + 4 : high & 3U];
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
uncompress_mmap (char *gbuffer1, UINT4 *blocks, Genomicpos_T startpos, 
		 Genomicpos_T endpos) {
  /* Genomicpos_T length = endpos - startpos; */
  Genomicpos_T startblock, endblock, startdiscard, enddiscard, ptr;
  UINT4 high, low, flags;
  char Buffer[32];
  int i, k = 0;

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
      Buffer[i] = translate[(flags & 1U) ? (low & 3U) + 4 : low & 3U];
      low >>= 2;
      flags >>= 1;
    }
    for ( ; i < 32; i++) {
      Buffer[i] = translate[(flags & 1U) ? (high & 3U) + 4 : high & 3U];
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
      Buffer[i] = translate[(flags & 1U) ? (low & 3U) + 4 : low & 3U];
      low >>= 2;
      flags >>= 1;
    }
    for ( ; i < 32; i++) {
      Buffer[i] = translate[(flags & 1U) ? (high & 3U) + 4 : high & 3U];
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
	gbuffer1[k++] = translate[(flags & 1U) ? (low & 3U) + 4 : low & 3U];
	low >>= 2;
	flags >>= 1;
      }
      for ( ; i < 32; i++) {
	gbuffer1[k++] = translate[(flags & 1U) ? (high & 3U) + 4 : high & 3U];
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
	Buffer[i] = translate[(flags & 1U) ? (low & 3U) + 4 : low & 3U];
	low >>= 2;
	flags >>= 1;
      }
      for ( ; i < 32; i++) {
	Buffer[i] = translate[(flags & 1U) ? (high & 3U) + 4 : high & 3U];
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

static const Except_T gbufferlen_error = { "Insufficient allocation" };

Sequence_T
Genome_get_segment (T this, Genomicpos_T left, Genomicpos_T length, 
		    bool revcomp, char *gbuffer1, char *gbuffer2, int gbufferlen) {
  
  if (length > gbufferlen) {
    fprintf(stderr,"Didn't allocate enough space for gbufferlen (%d < %u)\n",
	    gbufferlen,length);
    Except_raise(&gbufferlen_error,__FILE__,__LINE__);
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
      uncompress_fileio(gbuffer1,this,left,left+length);
#ifdef HAVE_PTHREAD
      pthread_mutex_unlock(&this->read_mutex);
#endif
    } else {
      uncompress_mmap(gbuffer1,this->blocks,left,left+length);
    }
  }
  gbuffer1[length] = '\0';

  if (revcomp == true) {
    debug(printf("Got sequence at %u with length %u, revcomp\n",left,length));
    make_complement_buffered(gbuffer2,gbuffer1,length);
    return Sequence_genomic_new(gbuffer2,length);
  } else {
    debug(printf("Got sequence at %u with length %u, forward\n",left,length));
    return Sequence_genomic_new(gbuffer1,length);
  }
}

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
  char *dest, *src, *matbuffer, *shiftdest, *shiftsrc;
  int index, i, matlen, patchlen, shiftlen, expansion;
  bool allocp;
  
  assert(reflen <= gbuffer3len);
  refR = refL + reflen;
  debug1(printf("refL=%u refR=%u reflen=%u\n",refL,refR,reflen));

  /* Work in gbuffer3 */
  memcpy(gbuffer3,gbuffer1,reflen);

  for (i = 0; i < nindices; i++) {
    index = indices[i];
    interval = IIT_interval(altstrain_iit,index);
    srcL = Interval_low(interval);
    srcR = Interval_high(interval) + 1;	/* Intervals are inclusive */
    matbuffer = IIT_annotation(altstrain_iit,index,&allocp);
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
      FREE(matbuffer);
    }
    debug1(printf("srcL=%u matR=%u srcR=%u matlen=%u patchlen=%d expansion=%d\n",
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
	debug1(printf("  shifted %d\n",shiftlen));
      }
    }
  }
  debug1(printf("\n"));

  if (revcomp == true) {
    make_complement_buffered(gbuffer2,gbuffer3,reflen);
    /* FREE(sequence); */
    debug(printf("Got sequence at %u with length %u, revcomp\n",refL,reflen));
    return Sequence_genomic_new(gbuffer2,reflen);
  } else {
    debug(printf("Got sequence at %u with length %u, forward\n",refL,reflen));
    return Sequence_genomic_new(gbuffer3,reflen);
  }
}

