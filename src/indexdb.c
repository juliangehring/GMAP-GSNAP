static char rcsid[] = "$Id: indexdb.c,v 1.63 2005/05/10 02:14:52 twu Exp $";
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
#ifdef HAVE_UNISTD_H
#include <unistd.h>		/* For mmap on Linux and for lseek */
#endif
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>		/* For mmap and off_t */
#endif
#ifdef HAVE_SYS_STAT_H
#include <sys/stat.h>
#endif
#ifdef HAVE_FCNTL_H
#include <fcntl.h>		/* For open */
#endif
#include <sys/mman.h>		/* For mmap */
/* Not sure why this was included
#include <errno.h>
*/
#ifdef PAGESIZE_VIA_SYSCTL
#include <sys/sysctl.h>
#endif
#include "mem.h"
#include "stopwatch.h"
#include "compress.h"

#ifdef HAVE_PTHREAD
#include <pthread.h>		/* sys/types.h already included above */
#endif

#define PAGESIZE 1024*4

/* An offset into the positions file of an IndexDB.  Typically, 3
   billion divided by sampling interval, requiring a maximum of 32
   bits or 4 bytes */
typedef UINT4 Positionsptr_T;


#define MAXENTRIES 20
#define BADVAL (Genomicpos_T) -1

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

#define T Indexdb_T
struct T {
  Positionsptr_T *offsets;

  int positions_fd;
  size_t positions_len;
  Genomicpos_T *positions;

  bool user_provided_p;
  bool batchp;
};


void
Indexdb_free (T *old) {
  if (*old) {
    if ((*old)->user_provided_p == true) {
      FREE((*old)->positions);
    } else {
#ifdef HAVE_MMAP
      if ((*old)->positions != NULL) {
	munmap((void *) (*old)->positions,(*old)->positions_len);
      }
#endif
      close((*old)->positions_fd);
    }
    FREE((*old)->offsets);
    FREE(*old);
  }
  return;
}


T
Indexdb_new_genome (char *genomesubdir, char *fileroot, bool batchp) {
  T new = (T) MALLOC(sizeof(*new));
  char *filename;
  FILE *fp;
  int offsets_fd;
  struct stat sb;
  int oligospace;
  int pagesize, indicesperpage, totalindices;
  size_t len, i;
  int nzero = 0, npos = 0;
  int mib[2];

  new->user_provided_p = false;
  new->batchp = batchp;

  /* Read offsets file */
  filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			     strlen(fileroot)+strlen(".idxoffsets")+1,sizeof(char));
  sprintf(filename,"%s/%s.idxoffsets",genomesubdir,fileroot);
  if ((offsets_fd = open(filename,O_RDONLY,0764)) < 0) {
    fprintf(stderr,"Error: can't open file %s\n",filename);
    exit(9);
  }
  fstat(offsets_fd,&sb);
  close(offsets_fd);

  if ((fp = fopen(filename,"r")) == NULL) {
    fprintf(stderr,"Error: can't open file %s\n",filename);
    exit(9);
  }
  new->offsets = (Positionsptr_T *) MALLOC(sb.st_size);
  if (batchp == true) {
    FREAD_UINTS(new->offsets,sb.st_size/sizeof(Positionsptr_T),fp);
  } else {
    /* Note: Need to convert to bigendian later, as needed */
    fread(new->offsets,sizeof(Positionsptr_T),sb.st_size/sizeof(Positionsptr_T),fp);
  }
  fclose(fp);
  FREE(filename);


  /* Read positions file */
  filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			     strlen(fileroot)+strlen(".idxpositions")+1,sizeof(char));
  sprintf(filename,"%s/%s.idxpositions",genomesubdir,fileroot);
  if ((new->positions_fd = open(filename,O_RDONLY,0764)) < 0) {
    fprintf(stderr,"Error: can't open file %s\n",filename);
    exit(9);
  }
  fstat(new->positions_fd,&sb);
  new->positions_len = sb.st_size;
#ifndef HAVE_MMAP
  new->positions = NULL;
#else

  if (batchp == true) {

#ifdef HAVE_GETPAGESIZE
    pagesize = getpagesize();
#elif defined(PAGESIZE_VIA_SYSCONF)
    pagesize = (int) sysconf(_SC_PAGESIZE);
#elif defined(PAGESIZE_VIA_SYSCTL)
    len = sizeof(pagesize);
    mib[0] = CTL_HW;
    mib[1] = HW_PAGESIZE;
    sysctl(mib,2,&pagesize,&len,NULL,0);
#else
    pagesize = PAGESIZE;
#endif

    indicesperpage = pagesize/sizeof(Genomicpos_T);
    fprintf(stderr,"Pre-loading index db every %d indices/page...",
	    indicesperpage);
    Stopwatch_start();
  }

  new->positions = (Genomicpos_T *) mmap(NULL,sb.st_size,PROT_READ,0
#ifdef HAVE_MMAP_MAP_SHARED
					 |MAP_SHARED
#endif
#ifdef HAVE_MMAP_MAP_FILE
					 |MAP_FILE
#endif
#ifdef HAVE_MMAP_MAP_VARIABLE
					 |MAP_VARIABLE
#endif
					 ,new->positions_fd,0);
  if (new->positions == MAP_FAILED) {
    new->positions = NULL;
    if (batchp == true) {
      Stopwatch_stop();
      fprintf(stderr,"insufficient memory (will use disk file instead)\n");
    }
  } else if (batchp == false) {
#ifdef HAVE_MADVISE
    madvise((caddr_t) new->positions,new->positions_len,MADV_RANDOM);
#endif
  } else {
    /* Touch all pages */
#ifdef HAVE_MADVISE
    madvise((caddr_t) new->positions,new->positions_len,MADV_WILLNEED);
#endif
    totalindices = new->positions_len/sizeof(Genomicpos_T);
    for (i = 0; i < totalindices; i += indicesperpage) {
      /* memcpy(temp,(caddr_t) new->positions + i,pagesize); */
      if (new->positions[i] == 0U) {
	nzero++;
      } else {
	npos++;
      }
      if (i % 10000 == 0) {
	fprintf(stderr,".");
      }
    }
    fprintf(stderr,"done (%d pages, %.2f sec)\n",nzero+npos,Stopwatch_stop());
  }
#endif /* HAVE_MMAP */
  FREE(filename);

  return new;
}


/************************************************************************
 *   Read procedures
 ************************************************************************/

static void
positions_move_absolute (T this, Positionsptr_T ptr) {
  off_t offset = ptr*((off_t) sizeof(Genomicpos_T));

  if (lseek(this->positions_fd,offset,SEEK_SET) < 0) {
    perror("Error in gmap, positions_move_absolute");
    exit(9);
  }
  return;
}

static Genomicpos_T
positions_read_current (T this) {
  Genomicpos_T value;
  char buffer[4];

  read(this->positions_fd,buffer,4);

  value = (buffer[3] & 0xff);
  value <<= 8;
  value |= (buffer[2] & 0xff);
  value <<= 8;
  value |= (buffer[1] & 0xff);
  value <<= 8;
  value |= (buffer[0] & 0xff);

  return value;
}

static Genomicpos_T
positions_read_backwards (T this) {
  Genomicpos_T value;
  char buffer[4];
  off_t reloffset = -2*((off_t) sizeof(Genomicpos_T)); /* 1 to undo the effect of read */

  read(this->positions_fd,buffer,4);

  value = (buffer[3] & 0xff);
  value <<= 8;
  value |= (buffer[2] & 0xff);
  value <<= 8;
  value |= (buffer[1] & 0xff);
  value <<= 8;
  value |= (buffer[0] & 0xff);

  if (lseek(this->positions_fd,reloffset,SEEK_CUR) < 0) {
    perror("Error in gmap, positions_read_backwards");
    exit(9);
  }
    
  return value;
}

#ifdef HAVE_PTHREAD
static pthread_mutex_t indexdb_read_mutex = PTHREAD_MUTEX_INITIALIZER;
#endif

/*                 87654321 */
#define LOW12MER 0x00FFFFFF

Genomicpos_T *
Indexdb_read (int *nentries, T this, Storedoligomer_T oligo) {
  Genomicpos_T *positions, position0, expected1, position1;
  Positionsptr_T ptr, ptr0, end0;
  int i;
  Storedoligomer_T part0;
  char byte1, byte2, byte3;
  int bigendian, littleendian;

  debug(printf("%06X\n",oligo));
  part0 = oligo & LOW12MER;

  /* Ignore poly A and poly T on stage 1 */
  if (part0 == 0U || part0 == LOW12MER) {
    return NULL;
  }

#ifdef WORDS_BIGENDIAN
  if (this->user_provided_p == true || this->batchp == true) {
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
  debug(printf("offset pointers are %u and %u\n",ptr0,end0));

  /* Skip backward over bad values, due to duplicates */
  if (this->positions == NULL) {
    positions_move_absolute(this,end0-1);
    while (end0 > ptr0 && positions_read_backwards(this) == BADVAL) {
      end0--;
    }
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
    if (this->positions == NULL) {
      /* non-mmap procedures */
#ifdef HAVE_PTHREAD
      pthread_mutex_lock(&indexdb_read_mutex);
#endif
      positions_move_absolute(this,ptr0);
      for (ptr = ptr0, i = 0; ptr < end0; ptr++, i++) {
	positions[i] = positions_read_current(this);
      }
#ifdef HAVE_PTHREAD
      pthread_mutex_unlock(&indexdb_read_mutex);
#endif
    } else {
      /* mmap procedures */
#ifdef WORDS_BIGENDIAN
      if (this->user_provided_p == true) {
	memcpy(positions,&(this->positions[ptr0]),(*nentries)*sizeof(Genomicpos_T));
      } else {
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


/************************************************************************
 *   Write procedures -- called by gmapindex
 ************************************************************************/

static int
power (int base, int exponent) {
  int result = 1, i;

  for (i = 0; i < exponent; i++) {
    result *= base;
  }
  return result;
}


#define MONITOR_INTERVAL 1000000

void
Indexdb_write_offsets (FILE *offsets_fp, FILE *sequence_fp, IIT_T altstrain_iit,
		       int index1interval, bool uncompressedp) {
  Positionsptr_T *offsets;
  char *p, b;
  int c, oligospace, index, i;
  Genomicpos_T position = 0;
  int between_counter = 0, in_counter = 0;
  Storedoligomer_T oligo = 0, masked, mask;
  Interval_T interval;

  mask = ~(~0UL << 2*INDEX1PART);
  oligospace = power(4,INDEX1PART);
  offsets = (Positionsptr_T *) CALLOC(oligospace+1,sizeof(Positionsptr_T));

  /* Handle reference strain */
  while ((c = Compress_get_char(sequence_fp,position,uncompressedp)) != EOF) {
    between_counter++;
    in_counter++;
    position++;

    if (position % MONITOR_INTERVAL == 0) {
      fprintf(stderr,"Indexing offsets of oligomers in genome (every %d), position %u\n",index1interval,position);
    }

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
	fprintf(stderr,"Bad character %c at position %u\n",c,position-1);
	abort();
      }
    }

    /*
    debug(printf("char=%c bc=%d ic=%d oligo=%016lX\n",
		 c,between_counter,in_counter,oligo));
    */

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
  }

  /* Handle alternate strains */
  if (altstrain_iit != NULL) {
    for (index = 1; index <= IIT_nintervals(altstrain_iit); index++) {
      interval = IIT_interval(altstrain_iit,index);

      position = Interval_low(interval);
      p = IIT_annotation(altstrain_iit,index); /* Holds the sequence */
      
      between_counter = in_counter = 0;
      oligo = 0U;
      while ((b = *p++) != '\0') {
	between_counter++;
	in_counter++;
	position++;

	switch (toupper(b)) {
	case 'A': oligo = (oligo << 2); break;
	case 'C': oligo = (oligo << 2) | 1; break;
	case 'G': oligo = (oligo << 2) | 2; break;
	case 'T': oligo = (oligo << 2) | 3; break;
	default: oligo = 0U; in_counter = 0; break;
	}
	
	/*
	  debug(printf("char=%c bc=%d ic=%d oligo=%016lX\n",
	  c,between_counter,in_counter,oligo));
	*/

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
      }
    }
  }

  for (i = 1; i <= oligospace; i++) {
    offsets[i] = offsets[i] + offsets[i-1];
    debug(printf("Offset for %06X: %u\n",i,offsets[i]));
  }

  FWRITE_UINTS(offsets,oligospace+1,offsets_fp);
  FREE(offsets);

  return;
}


static void
offsetsfile_move_absolute (FILE *fp, int ptr) {
  long int offset = (long int) (ptr*sizeof(Positionsptr_T));

  if (fseek(fp,offset,SEEK_SET) < 0) {
    perror("Error in gmapindex, offsetsfile_move_absolute");
    exit(9);
  }
  return;
}


static void
positionsfile_move_absolute (FILE *fp, Positionsptr_T ptr) {
  long int offset = (long int) (ptr*sizeof(Genomicpos_T));

  if (fseek(fp,offset,SEEK_SET) < 0) {
    perror("Error in gmapindex, positionsfile_move_absolute");
    exit(9);
  }
  return;
}

static Genomicpos_T
positionsfile_read_current (FILE *fp) {
  Genomicpos_T value;
  char buffer[4];

  fread(buffer,sizeof(char),4,fp);

  value = (buffer[3] & 0xff);
  value <<= 8;
  value |= (buffer[2] & 0xff);
  value <<= 8;
  value |= (buffer[1] & 0xff);
  value <<= 8;
  value |= (buffer[0] & 0xff);

  return value;
}

/* This doesn't seem to work */
/*
static void
positionsfile_move_relative (FILE *fp, int dist) {
  long int offset = (long int) (dist*sizeof(Genomicpos_T));

  if (fseek(fp,offset,SEEK_CUR) < 0) {
    fprintf(stderr,"Input to positionsfile_move_relative: %d\n",dist);
    perror("Error in gmapindex, positionsfile_move_relative");
    exit(9);
  }
  return;
}
*/


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
compute_positions_in_file (FILE *positions_fp, FILE *offsets_fp, Positionsptr_T *offsets,
			   FILE *sequence_fp, IIT_T altstrain_iit, int index1interval,
			   bool uncompressedp) {
  Positionsptr_T block_start, block_end, totalcounts;
  Genomicpos_T *positions_for_block, position = 0, adjposition, prevpos;
  char *p, b;
  int c, oligospace, index, i, npositions, nbadvals, j, k;
  int between_counter = 0, in_counter = 0;
  Storedoligomer_T oligo = 0U, masked, mask;
  Interval_T interval;

  mask = ~(~0UL << 2*INDEX1PART);
  oligospace = power(4,INDEX1PART);

  /* Handle reference strain */
  while ((c = Compress_get_char(sequence_fp,position,uncompressedp)) != EOF) {
    between_counter++;
    in_counter++;
    position++;

    if (position % MONITOR_INTERVAL == 0) {
      fprintf(stderr,"Indexing positions of oligomers in genome (every %d), position %u\n",index1interval,position);
    }

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
	fprintf(stderr,"Bad character %c at position %u\n",c,position-1U);
	abort();
      }
    }

    /*
    debug(printf("char=%c bc=%d ic=%d oligo=%016lX\n",
		 c,between_counter,in_counter,oligo));
    */
    
    if (in_counter == INDEX1PART) {
      if (between_counter >= index1interval) {
	masked = oligo & mask;
	positionsfile_move_absolute(positions_fp,offsets[masked]);
	offsets[masked] += 1;
	adjposition = position-INDEX1PART;
	FWRITE_UINT(adjposition,positions_fp);
	between_counter = 0;
      }
      in_counter--;
    }
    
  }

  /* Handle alternate strains */
  if (altstrain_iit != NULL) {
    for (index = 1; index <= IIT_nintervals(altstrain_iit); index++) {
      interval = IIT_interval(altstrain_iit,index);
      position = Interval_low(interval);
      p = IIT_annotation(altstrain_iit,index); /* Holds the sequence */
      
      between_counter = in_counter = 0;
      oligo = 0U;
      while ((b = *p++) != '\0') {
	between_counter++;
	in_counter++;
	position++;

	switch (toupper(b)) {
	case 'A': oligo = (oligo << 2); break;
	case 'C': oligo = (oligo << 2) | 1; break;
	case 'G': oligo = (oligo << 2) | 2; break;
	case 'T': oligo = (oligo << 2) | 3; break;
	default: oligo = 0U; in_counter = 0; break;
	}

	/*
	  debug(printf("char=%c bc=%d ic=%d oligo=%016lX\n",
	  c,between_counter,in_counter,oligo));
	*/
    
	if (in_counter == INDEX1PART) {
	  if (between_counter >= index1interval) {
	    masked = oligo & mask;
	    positionsfile_move_absolute(positions_fp,offsets[masked]);
	    offsets[masked] += 1;
	    adjposition = position-INDEX1PART;
	    FWRITE_UINT(adjposition,positions_fp);
	    between_counter = 0;
	  }
	  in_counter--;
	}
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
	positionsfile_move_absolute(positions_fp,block_start);
	FREAD_UINTS(positions_for_block,npositions,positions_fp);

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

	  positionsfile_move_absolute(positions_fp,block_start);
	  FWRITE_UINTS(positions_for_block,npositions,positions_fp);
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
			     bool uncompressedp) {
  Positionsptr_T block_start, block_end, j, k;
  Genomicpos_T position = 0, prevpos;
  char *p, b;
  int c, oligospace, index, i, npositions, nbadvals;
  int between_counter = 0, in_counter = 0;
  Storedoligomer_T oligo = 0U, masked, mask;
  Interval_T interval;

  mask = ~(~0UL << 2*INDEX1PART);
  oligospace = power(4,INDEX1PART);

  /* Handle reference strain */
  while ((c = Compress_get_char(sequence_fp,position,uncompressedp)) != EOF) {
    between_counter++;
    in_counter++;
    position++;

    if (position % MONITOR_INTERVAL == 0) {
      fprintf(stderr,"Indexing positions of oligomers in genome (every %d), position %u\n",index1interval,position);
    }

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
	fprintf(stderr,"Bad character %c at position %u\n",c,position-1U);
	abort();
      }
    }

    /*
    debug(printf("char=%c bc=%d ic=%d oligo=%016lX\n",
		 c,between_counter,in_counter,oligo));
    */
    
    if (in_counter == INDEX1PART) {
      if (between_counter >= index1interval) {
	masked = oligo & mask;
	positions[offsets[masked]++] = position-INDEX1PART;
	between_counter = 0;
      }
      in_counter--;
    }
    
  }

  /* Handle alternate strains */
  if (altstrain_iit != NULL) {
    for (index = 1; index <= IIT_nintervals(altstrain_iit); index++) {
      interval = IIT_interval(altstrain_iit,index);
      position = Interval_low(interval);
      p = IIT_annotation(altstrain_iit,index); /* Holds the sequence */
      
      between_counter = in_counter = 0;
      oligo = 0U;
      while ((b = *p++) != '\0') {
	between_counter++;
	in_counter++;
	position++;

	switch (toupper(b)) {
	case 'A': oligo = (oligo << 2); break;
	case 'C': oligo = (oligo << 2) | 1; break;
	case 'G': oligo = (oligo << 2) | 2; break;
	case 'T': oligo = (oligo << 2) | 3; break;
	default: oligo = 0U; in_counter = 0; break;
	}

	/*
	  debug(printf("char=%c bc=%d ic=%d oligo=%016lX\n",
	  c,between_counter,in_counter,oligo));
	*/
    
	if (in_counter == INDEX1PART) {
	  if (between_counter >= index1interval) {
	    masked = oligo & mask;
	    positions[offsets[masked]++] = position-INDEX1PART;
	    between_counter = 0;
	  }
	  in_counter--;
	}
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
			 bool writefilep) {
  FILE *positions_fp;
  Positionsptr_T *offsets, totalcounts;
  Genomicpos_T *positions;
  int oligospace;

  oligospace = power(4,INDEX1PART);
  offsets = (Positionsptr_T *) CALLOC(oligospace+1,sizeof(Positionsptr_T));
  FREAD_UINTS(offsets,oligospace+1,offsets_fp);
  totalcounts = offsets[oligospace];

  if (writefilep == true) {
    fprintf(stderr,"User requested build of positions in file\n");
    if ((positions_fp = fopen(positionsfile,"w+")) == NULL) {
      fprintf(stderr,"Can't open file %s\n",positionsfile);
      exit(9);
    }
    compute_positions_in_file(positions_fp,offsets_fp,offsets,sequence_fp,altstrain_iit,
			      index1interval,uncompressedp);
    fclose(positions_fp);

  } else {
    fprintf(stderr,"Trying to allocate %d*%d bytes of memory...",totalcounts,sizeof(Genomicpos_T));
    positions = (Genomicpos_T *) CALLOC_NO_EXCEPTION(totalcounts,sizeof(Genomicpos_T));
    if (positions == NULL) {
      fprintf(stderr,"failed.  Building positions in file.\n");
      if ((positions_fp = fopen(positionsfile,"w+")) == NULL) {
	fprintf(stderr,"Can't open file %s\n",positionsfile);
	exit(9);
      }
      compute_positions_in_file(positions_fp,offsets_fp,offsets,sequence_fp,altstrain_iit,
				index1interval,uncompressedp);
      fclose(positions_fp);

    } else {
      fprintf(stderr,"succeeded.  Building positions in memory.\n");
      if ((positions_fp = fopen(positionsfile,"w")) == NULL) {
	fprintf(stderr,"Can't open file %s\n",positionsfile);
	exit(9);
      }
      compute_positions_in_memory(positions,offsets_fp,offsets,sequence_fp,altstrain_iit,
				  index1interval,uncompressedp);
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
Indexdb_new_segment (char *genomicseg, int index1interval) {
  T new = (T) MALLOC(sizeof(*new));
  Positionsptr_T *work_offsets;	/* Working set for use in calculating positions */
  int totalcounts = 0;
  int c, oligospace, i;
  char *p;
  Genomicpos_T position = 0;
  int between_counter = 0, in_counter = 0;
  Storedoligomer_T oligo = 0U, masked, mask;

  new->user_provided_p = true;
  new->batchp = false;

  mask = ~(~0UL << 2*INDEX1PART);
  oligospace = power(4,INDEX1PART);

  /* Create offsets */
  new->offsets = (Positionsptr_T *) CALLOC(oligospace+1,sizeof(Positionsptr_T));

  p = genomicseg;
  while ((c = *(p++)) != '\0') {
    between_counter++;
    in_counter++;
    position++;

    switch (toupper(c)) {
    case 'A': oligo = (oligo << 2); break;
    case 'C': oligo = (oligo << 2) | 1; break;
    case 'G': oligo = (oligo << 2) | 2; break;
    case 'T': oligo = (oligo << 2) | 3; break;
    default: oligo = 0U; in_counter = 0; break;
    }

    /*
    debug(printf("char=%c bc=%d ic=%d oligo=%016lX\n",
		 c,between_counter,in_counter,oligo));
    */

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
  }

  for (i = 1; i <= oligospace; i++) {
    new->offsets[i] = new->offsets[i] + new->offsets[i-1];
    debug(printf("Offset for %06X: %u\n",i,new->offsets[i]));
  }

  /* Create positions */
  position = 0;
  between_counter = in_counter = 0;
  oligo = 0UL;

  work_offsets = (Positionsptr_T *) CALLOC(oligospace+1,sizeof(Positionsptr_T));
  for (i = 0; i <= oligospace; i++) {
    work_offsets[i] = new->offsets[i];
  }

  totalcounts = new->offsets[oligospace];
  new->positions = (Genomicpos_T *) CALLOC(totalcounts,sizeof(Genomicpos_T));

  p = genomicseg;
  while ((c = *(p++)) != '\0') {
    between_counter++;
    in_counter++;
    position++;

    switch (toupper(c)) {
    case 'A': oligo = (oligo << 2); break;
    case 'C': oligo = (oligo << 2) | 1; break;
    case 'G': oligo = (oligo << 2) | 2; break;
    case 'T': oligo = (oligo << 2) | 3; break;
    default: oligo = 0U; in_counter = 0; break;
    }

    /*
    debug(printf("char=%c bc=%d ic=%d oligo=%06X\n",
		 c,between_counter,in_counter,oligo));
    */
    
    if (in_counter == INDEX1PART) {
      if (between_counter >= index1interval) {
	masked = oligo & mask;
	new->positions[work_offsets[masked]++] = position-INDEX1PART;
	between_counter = 0;
      }
      in_counter--;
    }
    
  }

  FREE(work_offsets);

  return new;
}


