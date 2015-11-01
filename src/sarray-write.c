static char rcsid[] = "$Id: sarray-write.c 109825 2013-10-02 22:31:23Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "sarray-write.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/mman.h>		/* For munmap */
#include "bool.h"
#include "mem.h"
#include "genomicpos.h"
#include "assert.h"
#include "compress.h"
#include "bitpack64-write.h"
#include "bitpack64-read.h"
#include "fopen.h"
#include "saca-k.h"

#ifdef WORDS_BIGENDIAN
#include "bigendian.h"
#else
#include "littleendian.h"
#endif

/* #define WRITE_LCP 1 */

/* make_index */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif


static void
compute_lcp (UINT4 *lcp, unsigned char *s, UINT4 *SA, UINT4 n) {
  UINT4 *rank, h;
  UINT4 i, j;

  rank = (UINT4 *) CALLOC(n+1,sizeof(UINT4));
  
  for (i = 0; i <= n; i++) {
    rank[SA[i]] = i;
  }

  lcp[0] = 0;
  h = 0;
  for (i = 0; i <= n; i++) {
    if (rank[i] > 0) {
      j = SA[rank[i] - 1];
      while (i + h <= n && j + h <= n && s[i+h] == s[j+h]) {
	h++;
      }
      lcp[rank[i]] = h;
      if (h > 0) {
	h--;
      }
    }
  }

  FREE(rank);

  return;
}



static UINT4
power (int base, int exponent) {
  UINT4 result = 1U;
  int i;

  for (i = 0; i < exponent; i++) {
    result *= base;
  }
  return result;
}


#define LOW_TWO_BITS 0x3
#define RIGHT_A 0
#define RIGHT_C 1
#define RIGHT_G 2
#define RIGHT_T 3


static void
oligo_nt (char *nt, UINT4 oligo, int oligosize) {
  int i, j;
  UINT4 lowbits;

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

  return;
}



static UINT4
sarray_search_init (char *query, int querylength, Genome_T genome, UINT4 *SA, UINT4 *lcp,
		    UINT4 low, UINT4 high, UINT4 nmatches_low, UINT4 nmatches_high) {
  UINT4 mid, pos;
  UINT4 nmatches_mid;
  int fasti;
  char c;

  debug1(printf("sarray_search_init on querylength %d with low %u, high %u\n",querylength,low,high));
  while (low + 1 < high) {
    /* Compute mid for UINT4s */
    mid = low/2 + high/2;
    if (low % 2 == 1 && high % 2 == 1) {
      mid += 1;
    }
    nmatches_mid =  (nmatches_low < nmatches_high) ? nmatches_low : nmatches_high;

    fasti = nmatches_mid;
    pos = SA[mid] + nmatches_mid;
    while (fasti < querylength && (c = Genome_get_char(genome,pos)) == query[fasti]) {
      debug1(printf("Comparing query %d with pos %u\n",fasti,pos));
      fasti++;
      pos++;
    }
    if (c == 'N') {
      c = 'X';
    }

    if (fasti == querylength || c > query[fasti]) {
      high = mid;
      nmatches_high = (lcp[mid] < nmatches_mid) ? lcp[mid] : nmatches_mid;
    } else {
      low = mid;
      nmatches_low = (lcp[low] < nmatches_mid) ? lcp[low] : nmatches_mid;
    }

    debug1(printf("sarray_search_init with low %u, high %u\n",low,high));
  }

  debug1(printf("sarray_search_init ended.  Returning low %u+1\n\n",low));
  return low + 1;
}


static UINT4
sarray_find_index (bool *successp, char *query, int querylength, Genome_T genome,
		   UINT4 *SA, UINT4 *lcp, UINT4 n) {
  UINT4 low, high, mid, pos;
  UINT4 nmatches_low, nmatches_high, nmatches_mid;

  UINT4 prevlow, prevhigh;
  int nmatches_prevlow, nmatches_prevhigh, nmatches_best = 0;

  int fasti;
  char c;

  *successp = false;
  low = prevlow = 0;
  high = prevhigh = n;
  nmatches_low = nmatches_high = 0;

  debug1(printf("sarray_search on %s, querylength %d, with low %u, high %u\n",
	       query,querylength,low,high));
  while (low + 1 < high && *successp == false) {
    /* Compute mid for UINT4s */
    mid = low/2 + high/2;
    if (low % 2 == 1 && high % 2 == 1) {
      mid += 1;
    }
    nmatches_mid = (nmatches_low < nmatches_high) ? nmatches_low : nmatches_high;

    fasti = nmatches_mid;
    pos = SA[mid] + nmatches_mid;
    while (fasti < querylength && (c = Genome_get_char(genome,pos)) == query[fasti]) {
      debug1(printf("Comparing query %d with pos %u\n",fasti,pos));
      fasti++;
      pos++;
    }
    if (c == 'N') {
      c = 'X';
    }

    if (fasti > nmatches_best) {
      debug1(printf("fasti %d > nmatches_best %d.  Saving prevlow %u and prevhigh %u.\n",
		    fasti,nmatches_best,low,high));
      prevlow = low;
      prevhigh = high;
      nmatches_prevlow = nmatches_low;
      nmatches_prevhigh = nmatches_high;
      nmatches_best = fasti;
    }

    if (fasti == querylength) {
      *successp = true;

    } else if (c < query[fasti]) {
      low = mid;
      nmatches_low = (lcp[low] < nmatches_mid) ? lcp[low] : nmatches_mid;
      debug1(printf("genome %c < query (%c) => low gets %u @ %u\n",c,query[fasti],low,SA[low]));

    } else if (c > query[fasti]) {
      high = mid;
      nmatches_high = (lcp[mid] < nmatches_mid) ? lcp[mid] : nmatches_mid;
      debug1(printf("genome %c > query (%c) => high gets %u @ %u\n",c,query[fasti],high,SA[high]));

    } else {
      debug1(printf("genome %c == query (%c) => should not happen after Genome_consecutive_matches\n",
		   c,query[fasti]));
      abort();
    }

    debug1(printf("sarray_search with low %u @ %u, high %u @ %u\n",low,SA[low],high,SA[high]));
  }
  debug1(printf("\n"));

  if (*successp == false) {
    return -1U;
  } else {
    return sarray_search_init(query,querylength,genome,SA,lcp,
			      low,mid,nmatches_low,nmatches_mid);
  }
}



static UINT4 *
make_index (UINT4 *oligospace, int querylength, Genome_T genome, UINT4 *SA,
	    UINT4 *lcp, UINT4 n) {
  UINT4 *saindex;
  char *queryuc_ptr;
  UINT4 i;
  bool successp;

  *oligospace = power(4,querylength);
  saindex = (UINT4 *) CALLOC((*oligospace)+1,sizeof(UINT4));

  queryuc_ptr = (char *) CALLOC(querylength+1,sizeof(char));

  saindex[0] = 1U;
  for (i = 1; i < *oligospace; i++) {
    oligo_nt(queryuc_ptr,i,querylength);
#if 0
    if ((ptr = sarray_find_index(&successp,queryuc_ptr,querylength,genome,SA,lcp,n)) == -1) {
      saindex[i] = saindex[i-1];
    } else {
      saindex[i] = ptr;
    }
#else
    saindex[i] = sarray_find_index(&successp,queryuc_ptr,querylength,genome,SA,lcp,n);
#endif
    /* printf("%s\t%u\t%u\n",queryuc_ptr,i,saindex[i]); */
  }
  FREE(queryuc_ptr);

  saindex[*oligospace] = n;

  return saindex;
}


static int
compute_packsize (UINT4 *values) {
  UINT4 packsize;
  UINT4 maxvalue = 0, top;
  int i;
  int firstbit, msb;

  for (i = 0; i < 64; i++) {
    maxvalue |= values[i];
  }

#ifdef HAVE_BUILTIN_CLZ
  firstbit = __builtin_clz(maxvalue);
  packsize = 32 - firstbit;
#elif defined(HAVE_ASM_BSR)
  asm("bsr %1,%0" : "=r"(msb) : "r"(maxvalue));
  packsize = msb + 1;
#else
  firstbit = ((top = maxvalue >> 16) ? clz_table[top] : 16 + clz_table[maxvalue]);
  packsize = 32 - firstbit;
#endif

#ifdef ALLOW_ODD_PACKSIZES
  return packsize;
#else
  return (packsize + 1) & ~1;	/* Converts packsizes to the next multiple of 2 */
#endif
}

#define LCP_BUFFER_SIZE 1000000
#define METAINFO_SIZE 1
#define BLOCKSIZE 64

void
Sarray_write_lcpptrs (char *lcpptrsfile, char *lcpcompfile, UINT4 *lcp,
		      UINT4 /*n*/genomelength) {
  FILE *lcpptrs_fp, *lcpcomp_fp;
  UINT4 *lcpptrs;
  int lcpptri, i;
  UINT4 positioni;

  UINT4 *lcp_buffer;
  int lcp_buffer_size = LCP_BUFFER_SIZE;
  int lcp_buffer_i;

  UINT4 values[64];

  UINT4 nwritten;
  int packsize;


  Bitpack64_write_setup();

  /* 1 metavalue: nwritten (pointer).  Packsize can be
     computed from difference between successive pointers, if only
     even packsizes are allowed */
  lcpptrs = (UINT4 *) CALLOC(((genomelength + 1) + BLOCKSIZE - 1)/BLOCKSIZE * METAINFO_SIZE + 1,sizeof(UINT4));
  lcpptri = 0;

  if ((lcpcomp_fp = FOPEN_WRITE_BINARY(lcpcompfile)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",lcpcompfile);
    exit(9);
  }
  lcp_buffer = (UINT4 *) CALLOC(lcp_buffer_size,sizeof(UINT4));
  lcp_buffer_i = 0;

  nwritten = 0U;

  /* Last value of lcp is lcp[genomelength], with lcp[0] = 0. */
  for (positioni = 0; positioni + BLOCKSIZE <= genomelength; positioni += BLOCKSIZE) {
    /* Pointer */
    lcpptrs[lcpptri++] = nwritten;

    /* Pack block of 64 diffs */
    packsize = compute_packsize(&(lcp[positioni]));
    lcp_buffer_i = Bitpack64_write_vert(lcpcomp_fp,lcp_buffer,lcp_buffer_size,
					lcp_buffer_i,&(lcp[positioni]),packsize);

#ifdef ALLOW_ODD_PACKSIZES
    nwritten += 2 * ((packsize + 1) & ~1);
#else
    nwritten += 2 * packsize;
#endif
  }

  if (positioni <= genomelength) {
    /* Finish last block of 64 */
    lcpptrs[lcpptri++] = nwritten;
    
    i = 0;
    while (positioni <= genomelength) {
      values[i++] = lcp[positioni++];
    }
    while (i < BLOCKSIZE) {
      values[i++] = 0;
    }

    packsize = compute_packsize(values);
    lcp_buffer_i = Bitpack64_write_vert(lcpcomp_fp,lcp_buffer,lcp_buffer_size,
					lcp_buffer_i,values,packsize);

#ifdef ALLOW_ODD_PACKSIZES
    nwritten += 2 * ((packsize + 1) & ~1);
#else
    nwritten += 2 * packsize;
#endif
  }

  /* Write the final pointer, which will point after the end of the
     file */
  lcpptrs[lcpptri++] = nwritten;

  if ((lcpptrs_fp = FOPEN_WRITE_BINARY(lcpptrsfile)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",lcpptrsfile);
    exit(9);
  } else {
    FWRITE_UINTS(lcpptrs,lcpptri,lcpptrs_fp);
    FREE(lcpptrs);
    fclose(lcpptrs_fp);
  }
    
  /* Empty lcp buffer */
  if (lcp_buffer_i > 0) {
    FWRITE_UINTS(lcp_buffer,lcp_buffer_i,lcpcomp_fp);	
    lcp_buffer_i = 0;
  }
  FREE(lcp_buffer);
  fclose(lcpcomp_fp);

  return;
}


#if 0
/* Does not work, since Bitpack64_block_offsets is expecting cumulative sums */
static void
check_lcp_from_bitpack (char *lcpptrsfile, char *lcpcompfile, UINT4 *lcp, UINT4 n) {
  UINT4 *lcpptrs;
  UINT4 *lcpcomp;
  UINT4 lcp64[65];
  int lcpptrs_fd, lcpcomp_fd;
  size_t lcpptrs_len, lcpcomp_len;
  UINT4 pos, i;
#ifndef HAVE_MMAP
  double seconds;
#endif


#ifdef HAVE_MMAP
  lcpptrs = (UINT4 *) Access_mmap(&lcpptrs_fd,&lcpptrs_len,lcpptrsfile,sizeof(UINT4),/*randomp*/false);
  lcpcomp = (UINT4 *) Access_mmap(&lcpcomp_fd,&lcpcomp_len,lcpcompfile,sizeof(UINT4),/*randomp*/false);
#else
  lcpptrs = (UINT4 *) Access_allocated(&lcpptrs_len,&seconds,lcpptrsfile,sizeof(UINT4));
  lcpcomp = (UINT4 *) Access_allocated(&lcpcomp_len,&seconds,lcpcompfile,sizeof(UINT4));
#endif

  Bitpack64_read_setup(lcpptrs,lcpcomp,/*blocksize*/64);

  for (pos = 0UL; pos + /*blocksize*/64 <= n; pos += /*blocksize*/64) {
    Bitpack64_block_offsets(lcp64,lcpptrs,lcpcomp,/*blocksize*/64,pos);
    for (i = 0; i <= 64; i++) {
      if (lcp64[i] != lcp[pos+i]) {
#ifdef OLIGOSPACE_NOT_LONG
	fprintf(stderr,"Problem with lcp bitpack at oligo %u: %u != %u.  Please inform twu@gene.com\n",
		pos+i,lcp64[i],lcp[pos+i]);
#else
	fprintf(stderr,"Problem with lcp bitpack at oligo %lu: %u != %u.  Please inform twu@gene.com\n",
		pos+i,lcp64[i],lcp[pos+i]);
#endif
	exit(9);
      }
    }
  }

  if (pos <= n) {
    Bitpack64_block_offsets(lcp64,lcpptrs,lcpcomp,/*blocksize*/64,pos);
    i = 0;
    while (pos <= n) {
      if (lcp64[i] != lcp[pos]) {
#ifdef OLIGOSPACE_NOT_LONG
	fprintf(stderr,"Problem with lcp bitpack at oligo %u: %u != %u.  Please inform twu@gene.com\n",
		pos+i,lcp64[i],lcp[pos+i]);
#else
	fprintf(stderr,"Problem with lcp bitpack at oligo %lu: %u != %u.  Please inform twu@gene.com\n",
		pos+i,lcp64[i],lcp[pos+i]);
#endif
	exit(9);
      }
      pos++;
      i++;
    }
  }

#ifdef HAVE_MMAP
  munmap((void *) lcpcomp,lcpcomp_len);
  munmap((void *) lcpptrs,lcpptrs_len);
#else
  FREE(lcpcomp);
  FREE(lcpptrs);
#endif

  return;
}
#endif


void
Sarray_write_lcp (char *lcpfile, char *lcpptrsfile, char *lcpcompfile, char *saindexfile,
		  char *sarrayfile, Genome_T genome, UINT4 genomelength) {
  UINT4 *SA, *lcp, *saindex;
  UINT4 n;
  unsigned char *gbuffer;
  FILE *fp;
  UINT4 oligospace;
  int fd;
  size_t len;

  n = genomelength;
  gbuffer = (unsigned char *) CALLOC(n+1,sizeof(unsigned char));
  Genome_fill_buffer_int_string(genome,/*left*/0,/*length*/n,gbuffer);
  gbuffer[n] = 0;		       /* '\0', terminator */

  SA = (UINT4 *) Access_mmap(&fd,&len,sarrayfile,sizeof(UINT4),/*randomp*/true);

  fprintf(stderr,"Computing lcp...");
  lcp = (UINT4 *) MALLOC((n+1)*sizeof(UINT4));
  compute_lcp(lcp,gbuffer,SA,n);
  fprintf(stderr,"done\n");
  FREE(gbuffer);

#ifdef WRITE_LCP
  if ((fp = FOPEN_WRITE_BINARY(lcpfile)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",lcpfile);
    exit(9);
  } else {
    FWRITE_UINTS(lcp,n+1,fp);
    fclose(fp);
  }
#endif

  Sarray_write_lcpptrs(lcpptrsfile,lcpcompfile,lcp,n);
#if 0
  /* Check does not work */
  fprintf(stderr,"Checking lcpptrs...");
  check_lcp_from_bitpack(lcpptrsfile,lcpcompfile,lcp,n);
  fprintf(stderr,"done\n");
#endif

  fprintf(stderr,"Computing saindex...");
  saindex = make_index(&oligospace,/*querylength*/12,genome,SA,lcp,n);
  fprintf(stderr,"done\n");

  if ((fp = FOPEN_WRITE_BINARY(saindexfile)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",saindexfile);
    exit(9);
  } else {
    FWRITE_UINTS(saindex,oligospace+1,fp);
    fclose(fp);
  }

  FREE(saindex);
  FREE(lcp);

  munmap((void *) SA,len);
  close(fd);

  return;
}



void
Sarray_write_all (char *genomesubdir, char *fileroot, Genome_T genome, UINT4 genomelength) {
  UINT4 *SA, *lcp, *saindex;
  UINT4 n;
  unsigned char *gbuffer;
  FILE *fp;
  char *filename, *filename1;
  UINT4 oligospace;


  n = genomelength;
  SA = (UINT4 *) MALLOC((n+1)*sizeof(UINT4));
  gbuffer = (unsigned char *) CALLOC(n+1,sizeof(unsigned char));
  Genome_fill_buffer_int_string(genome,/*left*/0,/*length*/n,gbuffer);
  gbuffer[n] = 0;		       /* '\0', terminator */
  SACA_K(gbuffer,SA,n+/*virtual sentinel*/1,/*K, alphabet_size*/5,/*m*/n+1,/*level*/0);

  filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			     strlen(fileroot)+strlen(".sarray")+1,sizeof(char));
  sprintf(filename,"%s/%s.sarray",genomesubdir,fileroot);
  if ((fp = FOPEN_WRITE_BINARY(filename)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",filename);
    exit(9);
  } else {
    FWRITE_UINTS(SA,n+1,fp);
    fclose(fp);
    FREE(filename);
  }

  fprintf(stderr,"Computing lcp...");
  lcp = (UINT4 *) MALLOC((n+1)*sizeof(UINT4));
  compute_lcp(lcp,gbuffer,SA,n);
  fprintf(stderr,"done\n");
  FREE(gbuffer);

#ifdef WRITE_LCP
  filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			     strlen(fileroot)+strlen(".salcp")+1,sizeof(char));
  sprintf(filename,"%s/%s.salcp",genomesubdir,fileroot);
  if ((fp = FOPEN_WRITE_BINARY(filename)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",filename);
    exit(9);
  } else {
    FWRITE_UINTS(lcp,n+1,fp);
    fclose(fp);
    FREE(filename);
  }
#endif

  filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			     strlen(fileroot)+strlen(".salcpptrs")+1,sizeof(char));
  sprintf(filename,"%s/%s.salcpptrs",genomesubdir,fileroot);

  filename1 = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			     strlen(fileroot)+strlen(".salcpcomp")+1,sizeof(char));
  sprintf(filename1,"%s/%s.salcpcomp",genomesubdir,fileroot);

  Sarray_write_lcpptrs(filename,filename1,lcp,n);
  FREE(filename1);
  FREE(filename);


  fprintf(stderr,"Computing saindex...");
  saindex = make_index(&oligospace,/*querylength*/12,genome,SA,lcp,n);
  fprintf(stderr,"done\n");

  filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			     strlen(fileroot)+strlen(".saindex")+1,sizeof(char));
  sprintf(filename,"%s/%s.saindex",genomesubdir,fileroot);
  if ((fp = FOPEN_WRITE_BINARY(filename)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",filename);
    exit(9);
  } else {
    FWRITE_UINTS(saindex,oligospace+1,fp);
    fclose(fp);
  }

  FREE(saindex);

  FREE(lcp);
  FREE(SA);

  return;
}


