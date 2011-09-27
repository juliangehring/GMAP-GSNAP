static char rcsid[] = "$Id: genome.c,v 1.70 2005/05/10 02:14:32 twu Exp $";
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
#ifdef HAVE_UNISTD_H
#include <unistd.h>		/* For mmap on Linux and for lseek */
#endif
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>		/* For mmap and off_t */
#endif
#include <sys/mman.h>		/* For mmap */
#ifdef HAVE_FCNTL_H
#include <fcntl.h>		/* For open */
#endif
#ifdef HAVE_SYS_STAT_H
#include <sys/stat.h>		/* For open and fstat */
#endif
/* Not sure why this was included
#include <errno.h>
*/
#ifdef PAGESIZE_VIA_SYSCTL
#include <sys/sysctl.h>
#endif

#ifdef HAVE_PTHREAD
#include <pthread.h>		/* sys/types.h already included above */
#endif

#include "assert.h"
#include "mem.h"
#include "types.h"
#include "stopwatch.h"
#include "complement.h"
#include "interval.h"

#define PAGESIZE 1024*4

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
  int fd;
  size_t len;
  char *chars;
  UINT4 *blocks;
  bool compressedp;
};


void
Genome_free (T *old) {
  if (*old) {
    if ((*old)->compressedp == true) {
#ifdef HAVE_MMAP
      if ((*old)->blocks != NULL) {
	munmap((void *) (*old)->blocks,(*old)->len);
      }
#endif
      close((*old)->fd);
    } else {
#ifdef HAVE_MMAP
      if ((*old)->chars != NULL) {
	munmap((void *) (*old)->chars,(*old)->len);
      }
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
  FILE *fp;
  struct stat sb;
  caddr_t region;
  int pagesize, indicesperpage, totalindices;
  size_t len, i;
  int nzero = 0, npos = 0;
  int mib[2];
  bool compressedp = !uncompressedp;

  /* batchp = false; */

  if (compressedp == true) {
    filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+strlen(fileroot)+
			       strlen(".genomecomp")+1,sizeof(char));
    sprintf(filename,"%s/%s.genomecomp",genomesubdir,fileroot);
    if ((new->fd = open(filename,O_RDONLY,0764)) < 0) {
      fprintf(stderr,"Warning: can't open file %s.  Trying uncompressed genome.\n",filename);
      compressedp = false;
      FREE(filename);
      filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+strlen(fileroot)+
				 strlen(".genome")+1,sizeof(char));
      sprintf(filename,"%s/%s.genome",genomesubdir,fileroot);
      if ((new->fd = open(filename,O_RDONLY,0764)) < 0) {
	fprintf(stderr,"Error: can't open file %s.\n",filename);
	exit(9);
      }
    }
    FREE(filename);
  } else {
    filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+strlen(fileroot)+
			       strlen(".genome")+1,sizeof(char));
    sprintf(filename,"%s/%s.genome",genomesubdir,fileroot);
    if ((new->fd = open(filename,O_RDONLY,0764)) < 0) {
      fprintf(stderr,"Warning: can't open file %s.  Trying compressed genome.\n",filename);
      compressedp = true;
      FREE(filename);
      filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+strlen(fileroot)+
				 strlen(".genomecomp")+1,sizeof(char));
      sprintf(filename,"%s/%s.genomecomp",genomesubdir,fileroot);
      if ((new->fd = open(filename,O_RDONLY,0764)) < 0) {
	fprintf(stderr,"Error: can't open file %s.\n",filename);
	exit(9);
      }
    }
    FREE(filename);
  }
  fstat(new->fd,&sb);
  new->len = sb.st_size;

  new->compressedp = compressedp;
  if (batchp == true) {
#ifdef HAVE_GETPAGESIZE
    pagesize = getpagesize();
#elif PAGESIZE_VIA_SYSCONF
    pagesize = sysconf(_SC_PAGESIZE);
#elif PAGESIZE_VIA_SYSCTL
    len = sizeof(pagesize);
    mib[0] = CTL_HW;
    mib[1] = HW_PAGESIZE;
    sysctl(mib,2,&pagesize,&len,NULL,0);
#else
    pagesize = PAGESIZE;
#endif
  }

  if (compressedp == true) {
    new->chars = (char *) NULL;
#ifndef HAVE_MMAP
    new->blocks = NULL;
#else
    new->blocks = (UINT4 *) mmap(NULL,sb.st_size,PROT_READ,0
#ifdef HAVE_MMAP_MAP_SHARED
					|MAP_SHARED
#endif
#ifdef HAVE_MMAP_MAP_FILE
					|MAP_FILE
#endif
#ifdef HAVE_MMAP_MAP_VARIABLE
					|MAP_VARIABLE
#endif
					,new->fd,0);
    if (new->blocks == MAP_FAILED) {
      new->blocks = NULL;
    }
    region = (caddr_t) new->blocks;

#endif /* HAVE_MMAP */

  } else {
    new->blocks = (UINT4 *) NULL;
#ifndef HAVE_MMAP
    new->chars = NULL;
#else
    new->chars = (char *) mmap(NULL,sb.st_size,PROT_READ,0
#ifdef HAVE_MMAP_MAP_SHARED
			       |MAP_SHARED
#endif
#ifdef HAVE_MMAP_MAP_FILE
			       |MAP_FILE
#endif
#ifdef HAVE_MMAP_MAP_VARIABLE
			       |MAP_VARIABLE
#endif
			       ,new->fd,0);
    if (new->chars == MAP_FAILED) {
       new->chars = NULL;
    }
    region = (caddr_t) new->chars;

#endif /* HAVE_MMAP */

  }

#ifdef HAVE_MMAP
  if (batchp == false) {
    if (region != NULL) {
#ifdef HAVE_MADVISE
      madvise(region,new->len,MADV_DONTNEED);
#endif
    }
  } else if (region == NULL) {
    if (compressedp == true) {
      indicesperpage = pagesize/sizeof(UINT4);
      fprintf(stderr,"Pre-loading genome (compressed) every %d indices/page...",
	      indicesperpage);
      fprintf(stderr,"insufficient memory (will use disk file instead)\n");
    } else {
      indicesperpage = pagesize/sizeof(char);
      fprintf(stderr,"Pre-loading genome (full) every %d indices/page...",
	      indicesperpage);
      fprintf(stderr,"insufficient memory (will use disk file instead)\n");
    }
  } else {
    /* Touch all pages */

    if (compressedp == true) {
      indicesperpage = pagesize/sizeof(UINT4);
      totalindices = new->len/sizeof(UINT4);
      fprintf(stderr,"Pre-loading genome (compressed) every %d indices/page...",
	      indicesperpage);
      Stopwatch_start();

#ifdef HAVE_MADVISE
      madvise(region,new->len,MADV_WILLNEED);
#endif
      for (i = 0; i < totalindices; i += indicesperpage) {
	/* memcpy(temp,region + i,pagesize); */
	if (new->blocks[i] == 0U) {
	  nzero++;
	} else {
	  npos++;
	}
	if (i % 10000 == 0) {
	  fprintf(stderr,".");
	}
      }

    } else {
      indicesperpage = pagesize/sizeof(char);
      totalindices = new->len/sizeof(char);
      fprintf(stderr,"Pre-loading genome (full) every %d indices/page...",
	      indicesperpage);
      Stopwatch_start();

#ifdef HAVE_MADVISE
      madvise(region,new->len,MADV_WILLNEED);
#endif
      for (i = 0; i < totalindices; i += indicesperpage) {
	/* memcpy(temp,region + i,pagesize); */
	if (new->chars[i] == '0') {
	  nzero++;
	} else {
	  npos++;
	}
	if (i % 10000 == 0) {
	  fprintf(stderr,".");
	}
      }
    }

    fprintf(stderr,"done (%d pages, %.2f sec)\n",nzero+npos, Stopwatch_stop());
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
    complement[j] = complCode[sequence[i]];
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
Genome_replace_x () {
  translate[7] = 'N';
}


static void
uncompress_without_mmap (char *gbuffer1, T this, Genomicpos_T startpos, 
			 Genomicpos_T endpos) {
  char *sequence;
  Genomicpos_T length;
  Genomicpos_T startblock, endblock, startdiscard, enddiscard, ptr;
  UINT4 high, low, flags;
  char Buffer[32];
  int i, j, k = 0;

  length = endpos - startpos;
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
  char *sequence;
  Genomicpos_T length;
  Genomicpos_T startblock, endblock, startdiscard, enddiscard, ptr;
  UINT4 high, low, flags;
  char Buffer[32];
  int i, j, k = 0;

  length = endpos - startpos;
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

#ifdef HAVE_PTHREAD
static pthread_mutex_t genome_read_mutex = PTHREAD_MUTEX_INITIALIZER;
#endif


Sequence_T
Genome_get_segment (T this, Genomicpos_T left, Genomicpos_T length, 
		    bool revcomp, char *gbuffer1, char *gbuffer2, int gbufferlen) {
  char *sequence, *complement;
  
  if (length > gbufferlen) {
    fprintf(stderr,"Didn't allocate enough space for gbufferlen (%d < %d)\n",
	    gbufferlen,length);
    abort();
  }

  if (this->compressedp == false) {
    /* sequence = (char *) CALLOC(length+1,sizeof(char)); */
    if (this->chars == NULL) {
      /* non-mmap procedure, uncompressed genome */
#ifdef HAVE_PTHREAD
      pthread_mutex_lock(&genome_read_mutex);
#endif
      if (lseek(this->fd,left,SEEK_SET) < 0) {
	perror("Error in gmap, Genome_get_segment");
	exit(9);
      }
      read(this->fd,gbuffer1,length);
#ifdef HAVE_PTHREAD
      pthread_mutex_unlock(&genome_read_mutex);
#endif

    } else {
      /* mmap procedure, uncompressed genome */
      memcpy(gbuffer1,&(this->chars[left]),length*sizeof(char));
    }
  } else if (this->blocks == NULL) {
    /* non-mmap procedure, compressed genome */
#ifdef HAVE_PTHREAD
    pthread_mutex_lock(&genome_read_mutex);
#endif
    uncompress_without_mmap(gbuffer1,this,left,left+length);
#ifdef HAVE_PTHREAD
    pthread_mutex_unlock(&genome_read_mutex);
#endif
  } else {
    /* mmap procedure, compressed genome */
    uncompress_mmap(gbuffer1,this->blocks,left,left+length);
  }
  gbuffer1[length] = '\0';

  if (revcomp == true) {
    make_complement_buffered(gbuffer2,gbuffer1,length);
    /* FREE(sequence); */

    return Sequence_genomic_new(gbuffer2,length);
  } else {
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
    matbuffer = IIT_annotation(altstrain_iit,index);
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

