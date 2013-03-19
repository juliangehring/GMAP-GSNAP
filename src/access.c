static char rcsid[] = "$Id: access.c 78503 2012-11-06 21:31:24Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "access.h"

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>		/* For strerror */
#include <errno.h>

/* <unistd.h> and <sys/types.h> included in access.h */
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
#ifdef PAGESIZE_VIA_SYSCONF
#include <unistd.h>
#endif
#ifdef PAGESIZE_VIA_SYSCTL
#include <sys/sysctl.h>
#endif

#include "assert.h"
#include "mem.h"
#include "types.h"
#include "fopen.h"
#include "stopwatch.h"

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


bool
Access_file_exists_p (char *filename) {
#ifdef HAVE_STRUCT_STAT64
  struct stat64 sb;
#else
  struct stat sb;
#endif

#ifdef HAVE_STAT64
  if (stat64(filename,&sb) == 0) {
    return true;
  } else {
    return false;
  }
#else
  if (stat(filename,&sb) == 0) {
    return true;
  } else {
    return false;
  }
#endif
}


off_t
Access_filesize (char *filename) {
#ifdef HAVE_STRUCT_STAT64
  struct stat64 sb;
#else
  struct stat sb;
#endif

#ifdef HAVE_STAT64
  stat64(filename,&sb);
#else
  stat(filename,&sb);
#endif
  debug(printf("filesize is %lu\n",sb.st_size));
  return sb.st_size;
}


size_t
Access_file_copy (char *dest_file, char *source_file) {
  size_t nbytes = 0;
  FILE *dest, *source;
  int c;

  if ((source = FOPEN_READ_BINARY(source_file)) == NULL) {
    fprintf(stderr,"Cannot open source file %s\n",source_file);
    return 0;

  } else if ((dest = FOPEN_WRITE_BINARY(dest_file)) == NULL) {
    fprintf(stderr,"Cannot open destination file %s\n",dest_file);
    fclose(source);
    return 0;

  } else {
    while ((c = fgetc(source)) != EOF) {
      fputc(c,dest);
      nbytes++;
    }
    fclose(dest);
    fclose(source);
    return nbytes;
  }
}


bool
Access_file_equal (char *file1, char *file2) {
  FILE *fp1, *fp2;
  int c1, c2;

  if ((fp1 = FOPEN_READ_BINARY(file1)) == NULL) {
    fprintf(stderr,"Cannot open file %s\n",file1);
    exit(9);

  } else if ((fp2 = FOPEN_READ_BINARY(file2)) == NULL) {
    fprintf(stderr,"Cannot open file %s\n",file2);
    fclose(fp1);
    exit(9);

  } else {
    c1 = fgetc(fp1);
    c2 = fgetc(fp2);
    while (c1 != EOF && c2 != EOF) {
      if (c1 != c2) {
	fclose(fp2);
	fclose(fp1);
	return false;
      }
      c1 = fgetc(fp1);
      c2 = fgetc(fp2);
    }
    fclose(fp2);
    fclose(fp1);

    if (c1 == EOF && c2 == EOF) {
      return true;
    } else {
      return false;
    }
  }
}


int
Access_fileio (char *filename) {
  int fd;

  if ((fd = open(filename,O_RDONLY,0764)) < 0) {
    fprintf(stderr,"Error: can't open file %s with open for reading\n",filename);
    exit(9);
  }
  return fd;
}

int
Access_fileio_rw (char *filename) {
  int fd;

  if ((fd = open(filename,O_RDWR | O_CREAT | O_TRUNC,0764)) < 0) {
    fprintf(stderr,"Error: can't open file %s with open for reading/writing\n",filename);
    exit(9);
  }
  return fd;
}


#ifndef WORDS_BIGENDIAN
static unsigned int
first_nonzero (int *i, char *filename) {
  size_t len;
  FILE *fp;
  unsigned int value = 0;
  void *p;

  len = (size_t) Access_filesize(filename);

  if ((fp = FOPEN_READ_BINARY(filename)) == NULL) {
    fprintf(stderr,"Error: can't open file %s with fopen\n",filename);
    exit(9);
  } else {
    *i = 0;
    p = (void *) &value;
    while ((size_t) *i < len && fread(p,sizeof(unsigned int),1,fp) > 0 &&
	   value == 0) {
      *i += 1;
    }
    if (value == 0) {
      *i = -1;
    }
    fclose(fp);
    return value;
  }
}
#endif


#define FREAD_BATCH 10000000	/* 10 million at a time */

/* Bigendian conversion not needed after this */
void *
Access_allocated (size_t *len, double *seconds, char *filename, size_t eltsize) {
  void *memory;
  FILE *fp;
  Stopwatch_T stopwatch;
  unsigned int value;
  void *p;
  int i;

  *len = (size_t) Access_filesize(filename);

  if ((fp = FOPEN_READ_BINARY(filename)) == NULL) {
    fprintf(stderr,"Error: can't open file %s with fopen\n",filename);
    exit(9);
  }

  Stopwatch_start(stopwatch = Stopwatch_new());
  memory = (void *) MALLOC(*len);
  FREAD_UINTS(memory,(*len)/eltsize,fp);
  fclose(fp);

#ifndef WORDS_BIGENDIAN
  value = first_nonzero(&i,filename);
  if (i >= 0 && ((unsigned int *) memory)[i] != value) {
    fprintf(stderr,"single fread command failed (observed on Macs with -B 3 or greater on large genomes)...reading file in smaller batches...");
    fp = FOPEN_READ_BINARY(filename);

    for (i = 0; i < (int) ((*len)/eltsize); i += FREAD_BATCH) {
      p = (void *) &(((unsigned int *) memory)[i]);
      fread(p,sizeof(unsigned int),FREAD_BATCH,fp);
    }

    if (i < (int) (*len)/eltsize) {
      p = (void *) &(((unsigned int *) memory)[i]);
      fread(p,sizeof(unsigned int),(*len)/eltsize - i,fp);
    }

    fclose(fp);
  }
#endif

  /* Note: the following (old non-batch mode) requires conversion to bigendian later, as needed */
  /* fread(new->offsets,eltsize,sb.st_size/eltsize,fp); */

  *seconds = Stopwatch_stop(stopwatch);
  Stopwatch_free(&stopwatch);

  return memory;
}


#define PAGESIZE 1024*4

static int
get_pagesize () {

#ifdef PAGESIZE_VIA_SYSCTL
  int pagesize;
  size_t pagelen;
  int mib[2];
#endif

#ifdef __STRICT_ANSI__
  return PAGESIZE;
#elif defined(HAVE_GETPAGESIZE)
  return getpagesize();
#elif defined(PAGESIZE_VIA_SYSCONF)
  return (int) sysconf(_SC_PAGESIZE);
#elif defined(PAGESIZE_VIA_SYSCTL)
  pagelen = sizeof(pagesize);
  mib[0] = CTL_HW;
  mib[1] = HW_PAGESIZE;
  sysctl(mib,2,&pagesize,&pagelen,NULL,0);
  return pagesize;
#else
  return PAGESIZE;
#endif

}


#ifdef HAVE_MMAP
/* Returns NULL if mmap fails.  Bigendian conversion required */
#ifdef HAVE_CADDR_T
caddr_t
#else
void *
#endif
Access_mmap (int *fd, size_t *len, char *filename, size_t eltsize, bool randomp) {
  off_t length;
#ifdef HAVE_CADDR_T
  caddr_t memory;
#else
  void *memory;
#endif

  if ((*len = length = Access_filesize(filename)) == 0U) {
    fprintf(stderr,"Error: file %s is empty\n",filename);
    exit(9);
  } 

  if ((*fd = open(filename,O_RDONLY,0764)) < 0) {
    fprintf(stderr,"Error: can't open file %s with open for reading\n",filename);
    exit(9);
  }

  if (sizeof(size_t) <= 4 && length > MAX32BIT) {
    debug(printf("Too big to mmap\n"));
    *len = 0;
    memory = NULL;
  } else {
    *len = (size_t) length;
    memory = mmap(NULL,length,PROT_READ,0
#ifdef HAVE_MMAP_MAP_SHARED
		  |MAP_SHARED
#endif
#ifdef HAVE_MMAP_MAP_FILE
		  |MAP_FILE
#endif
#ifdef HAVE_MMAP_MAP_VARIABLE
		  |MAP_VARIABLE
#endif
		  ,*fd,0);
    if (memory == MAP_FAILED) {
      fprintf(stderr,"Got mmap failure on len %ju from length %ju.  Error %d: %s\n",
	      length,length,errno,strerror(errno));
      debug(printf("Got MAP_FAILED on len %lu from length %lu\n",length,length));
      memory = NULL;
    } else if (randomp == true) {
      debug(printf("Got mmap of %lu bytes at %p to %p\n",length,memory,memory+length-1));
#ifdef HAVE_MADVISE
#ifdef HAVE_MADVISE_MADV_RANDOM
      madvise(memory,*len,MADV_RANDOM);
#endif
#endif
    } else {
      debug(printf("Got mmap of %lu bytes at %p to %p\n",length,memory,memory+length-1));
#ifdef HAVE_MADVISE
#ifdef HAVE_MADVISE_MADV_DONTNEED
      madvise(memory,*len,MADV_DONTNEED);
#endif
#endif
    }
  }

  return memory;
}
#endif



#ifdef HAVE_MMAP
/* Returns NULL if mmap fails.  Bigendian conversion required */
#ifdef HAVE_CADDR_T
caddr_t
#else
void *
#endif
Access_mmap_offset (int *remainder, int fd, off_t offset, size_t length, size_t eltsize, bool randomp) {
#ifdef HAVE_CADDR_T
  caddr_t memory;
#else
  void *memory;
#endif

  if (length == 0) {
    abort();
  }

  *remainder = offset % get_pagesize();
  offset -= (off_t) *remainder;
  length += (size_t) *remainder;

  if (sizeof(size_t) <= 4 && length > MAX32BIT) {
    debug(printf("Too big to mmap\n"));
    memory = NULL;
  } else {
    memory = mmap(NULL,length,PROT_READ,0
#ifdef HAVE_MMAP_MAP_SHARED
		  |MAP_SHARED
#endif
#ifdef HAVE_MMAP_MAP_FILE
		  |MAP_FILE
#endif
#ifdef HAVE_MMAP_MAP_VARIABLE
		  |MAP_VARIABLE
#endif
		  ,fd,offset);
    if (memory == MAP_FAILED) {
      fprintf(stderr,"Got mmap failure on fd %d, offset %ju, length %ju.  Error %d: %s\n",
	      fd,offset,length,errno,strerror(errno));
      debug(printf("Got MAP_FAILED on fd %d, offset %lu, length %lu\n",fd,offset,length));
      memory = NULL;
    } else if (randomp == true) {
      debug(printf("Got mmap of %lu bytes at %p to %p\n",length,memory,memory+length-1));
#ifdef HAVE_MADVISE
#ifdef HAVE_MADVISE_MADV_RANDOM
      madvise(memory,length,MADV_RANDOM);
#endif
#endif
    } else {
      debug(printf("Got mmap of %lu bytes at %p to %p\n",length,memory,memory+length-1));
#ifdef HAVE_MADVISE
#ifdef HAVE_MADVISE_MADV_DONTNEED
      madvise(memory,length,MADV_DONTNEED);
#endif
#endif
    }
  }

  return memory;
}
#endif



#ifdef HAVE_MMAP
/* Returns NULL if mmap fails.  Bigendian conversion required */
#ifdef HAVE_CADDR_T
caddr_t
#else
void *
#endif
Access_mmap_rw (int *fd, size_t *len, char *filename, size_t eltsize, bool randomp) {
  off_t length;
#ifdef HAVE_CADDR_T
  caddr_t memory;
#else
  void *memory;
#endif

  if ((*len = length = Access_filesize(filename)) == 0U) {
    fprintf(stderr,"Error: file %s is empty\n",filename);
    exit(9);
  }

  if ((*fd = open(filename,O_RDWR,0764)) < 0) {
    fprintf(stderr,"Error: can't open file %s with open for reading/writing\n",filename);
    exit(9);
  }

  if (sizeof(size_t) <= 4 && length > MAX32BIT) {
    debug(printf("Too big to mmap\n"));
    *len = 0;
    memory = NULL;
  } else {
    *len = (size_t) length;
    memory = mmap(NULL,length,PROT_READ|PROT_WRITE,0
#ifdef HAVE_MMAP_MAP_SHARED
		  |MAP_SHARED
#endif
#ifdef HAVE_MMAP_MAP_FILE
		  |MAP_FILE
#endif
#ifdef HAVE_MMAP_MAP_VARIABLE
		  |MAP_VARIABLE
#endif
		  ,*fd,0);
    if (memory == MAP_FAILED) {
      fprintf(stderr,"Got mmap failure on len %ju from length %ju.  Error %d: %s\n",
	      *len,length,errno,strerror(errno));
      debug(printf("Got MAP_FAILED on len %lu from length %lu\n",*len,length));
      memory = NULL;
    } else if (randomp == true) {
      debug(printf("Got mmap of %lu bytes at %p to %p\n",length,memory,memory+length-1));
#ifdef HAVE_MADVISE
#ifdef HAVE_MADVISE_MADV_RANDOM
      madvise(memory,*len,MADV_RANDOM);
#endif
#endif
    } else {
      debug(printf("Got mmap of %lu bytes at %p to %p\n",length,memory,memory+length-1));
#ifdef HAVE_MADVISE
#ifdef HAVE_MADVISE_MADV_DONTNEED
      madvise(memory,*len,MADV_DONTNEED);
#endif
#endif
    }
  }

  return memory;
}
#endif

#ifdef HAVE_MMAP
/* Returns NULL if mmap fails.  Bigendian conversion required */
#ifdef HAVE_CADDR_T
caddr_t
#else
void *
#endif
Access_mmap_offset_rw (int *remainder, int fd, off_t offset, size_t length, size_t eltsize, bool randomp) {
#ifdef HAVE_CADDR_T
  caddr_t memory;
#else
  void *memory;
#endif

  if (length == 0) {
    abort();
  }

  *remainder = offset % get_pagesize();
  offset -= (off_t) *remainder;
  length += (size_t) *remainder;

  if (sizeof(size_t) <= 4 && length > MAX32BIT) {
    debug(printf("Too big to mmap\n"));
    memory = NULL;
  } else {
    memory = mmap(NULL,length,PROT_READ|PROT_WRITE,0
#ifdef HAVE_MMAP_MAP_SHARED
		  |MAP_SHARED
#endif
#ifdef HAVE_MMAP_MAP_FILE
		  |MAP_FILE
#endif
#ifdef HAVE_MMAP_MAP_VARIABLE
		  |MAP_VARIABLE
#endif
		  ,fd,offset);
    if (memory == MAP_FAILED) {
      fprintf(stderr,"Got mmap failure on offset %ju, length %ju.  Error %d: %s\n",
	      offset,length,errno,strerror(errno));
      debug(printf("Got MAP_FAILED on offset %lu, length %lu\n",offset,length));
      memory = NULL;
    } else if (randomp == true) {
      debug(printf("Got mmap of %lu bytes at %p to %p\n",length,memory,memory+length-1));
#ifdef HAVE_MADVISE
#ifdef HAVE_MADVISE_MADV_RANDOM
      madvise(memory,length,MADV_RANDOM);
#endif
#endif
    } else {
      debug(printf("Got mmap of %lu bytes at %p to %p\n",length,memory,memory+length-1));
#ifdef HAVE_MADVISE
#ifdef HAVE_MADVISE_MADV_DONTNEED
      madvise(memory,length,MADV_DONTNEED);
#endif
#endif
    }
  }

  return memory;
}
#endif




#ifdef HAVE_MMAP

#ifdef HAVE_CADDR_T
caddr_t
#else
void *
#endif
Access_mmap_and_preload (int *fd, size_t *len, int *npages, double *seconds, char *filename, size_t eltsize) {
  off_t length;
#ifdef HAVE_CADDR_T
  caddr_t memory;
#else
  void *memory;
#endif
  int pagesize, indicesperpage;
  size_t totalindices, i;	/* Needs to handle uncompressed genomes > 2 gigabytes */
  int nzero = 0, npos = 0;
  Stopwatch_T stopwatch;


  if ((*len = length = Access_filesize(filename)) == 0U) {
    fprintf(stderr,"Error: file %s is empty\n",filename);
    exit(9);
  }

  if ((*fd = open(filename,O_RDONLY,0764)) < 0) {
    fprintf(stderr,"Error: can't open file %s with open for reading\n",filename);
    exit(9);
  }

  if (sizeof(size_t) <= 4 && *len > MAX32BIT) {
    debug(printf("Too big to mmap\n"));
    *len = 0;
    *npages = 0;
    *seconds = 0.0;
    memory = NULL;

  } else {

    pagesize = get_pagesize();

    indicesperpage = pagesize/eltsize;
    
    Stopwatch_start(stopwatch = Stopwatch_new());

    memory = mmap(NULL,length,PROT_READ,0
#ifdef HAVE_MMAP_MAP_SHARED
		  |MAP_SHARED
#endif
#ifdef HAVE_MMAP_MAP_FILE
		  |MAP_FILE
#endif
#ifdef HAVE_MMAP_MAP_VARIABLE
		  |MAP_VARIABLE
#endif
		  ,*fd,0);
    if (memory == MAP_FAILED) {
      fprintf(stderr,"Got mmap failure on len %ju from length %ju.  Error %d: %s\n",
	      *len,length,errno,strerror(errno));
      debug(printf("Got MAP_FAILED on len %lu from length %lu\n",*len,length));
      memory = NULL;
      Stopwatch_stop(stopwatch);
      Stopwatch_free(&stopwatch);
    } else {
      /* Touch all pages */
      debug(printf("Got mmap of %lu bytes at %p to %p\n",length,memory,memory+length-1));
#ifdef HAVE_MADVISE
#ifdef HAVE_MADVISE_MADV_WILLNEED
      madvise(memory,*len,MADV_WILLNEED);
#endif
#endif
      totalindices = (*len)/eltsize;
      for (i = 0; i < totalindices; i += indicesperpage) {
	if (((char *) memory)[i] == 0) {
	  nzero++;
	} else {
	  npos++;
	}
	if (i % 10000 == 0) {
	  fprintf(stderr,".");
	}
      }
      *npages = nzero + npos;
      *seconds = Stopwatch_stop(stopwatch);
      Stopwatch_free(&stopwatch);
    }
  }

  return memory;
}
#endif


