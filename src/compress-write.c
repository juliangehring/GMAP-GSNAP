static char rcsid[] = "$Id: compress-write.c 101503 2013-07-15 18:26:20Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifndef HAVE_MEMCPY
# define memcpy(d,s,n) bcopy((s),(d),(n))
#endif
#ifndef HAVE_MEMMOVE
# define memmove(d,s,n) bcopy((s),(d),(n))
#endif

#include "compress-write.h"

#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <ctype.h>		/* For isalpha, toupper */
#ifdef WORDS_BIGENDIAN
#include "bigendian.h"
#else
#include "littleendian.h"
#endif
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>		/* For off_t */
#endif
#include "complement.h"
#include "mem.h"		/* For Compress_new */


/* Another MONITOR_INTERVAL is in indexdb.c */
#define MONITOR_INTERVAL 10000000 /* 10 million nt */

#define MAX_BADCHAR_MESSAGES 10
#define BADCHAR_INTERVAL 1000000


static char uppercaseCode[128] = UPPERCASE_U2T;

/* We use int *, rather than char *, because we eventually return an int,
   and we see problems converting from char to int */
static void
fill_buffer (int *Buffer, Genomecomp_T high, Genomecomp_T low, Genomecomp_T flags, Univcoord_T position) {
  int i;

  /* printf("%08X %08X %08X => ",high,low,flags); */
  for (i = 0; i < 16; i++) {
    switch (low & 3U) {
    case 0U: Buffer[i] = 'A'; break;
    case 1U: Buffer[i] = 'C'; break;
    case 2U: Buffer[i] = 'G'; break;
    case 3U: Buffer[i] = 'T'; break;
    default: abort();
    }
    low >>= 2;
  }
  for ( ; i < 32; i++) {
    switch (high & 3U) {
    case 0U: Buffer[i] = 'A'; break;
    case 1U: Buffer[i] = 'C'; break;
    case 2U: Buffer[i] = 'G'; break;
    case 3U: Buffer[i] = 'T'; break;
    default: abort();
    }
    high >>= 2;
  }
  for (i = 0; i < 32; i++) {
    if ((flags & 1U) == 1U) {
      if (Buffer[i] == 'A') {
	Buffer[i] = 'N';
      } else if (Buffer[i] == 'T') {
	Buffer[i] = 'X';
      } else {
	printf("Parsing error; saw non-ACGT flag plus %c at position %lu\n",Buffer[i],position+i);
	exit(9);
      }
    }
    flags >>= 1;
  }

  return;
}


/* Based on genomecomp */
int
Compress_get_char (FILE *sequence_fp, Univcoord_T position, bool uncompressedp) {
  Genomecomp_T high, low, flags;
  static int SAVEBUFFER[32];
  int ptr, c;

  if (uncompressedp == true) {
    while ((c = fgetc(sequence_fp)) != EOF && isspace(c)) {
    }
    if (c == EOF) {
      return EOF;
    } else {
      return c;
    }
  } else if ((ptr = position % 32) == 0) {
    if (FREAD_UINT(&high,sequence_fp) <= 0 ||
	FREAD_UINT(&low,sequence_fp) <= 0 ||
	FREAD_UINT(&flags,sequence_fp) <= 0) {
      return EOF;
    } else {
      fill_buffer(SAVEBUFFER,high,low,flags,position);
      return SAVEBUFFER[0];
    }
  } else {
    return SAVEBUFFER[ptr];
  }
}


/************************************************************************
 *   Compression and uncompression of the genome
 ************************************************************************/

/*               87654321 */
#define LEFT_A 0x00000000
#define LEFT_C 0x40000000
#define LEFT_G 0x80000000
#define LEFT_T 0xC0000000
#define LEFT_BIT 0x80000000

/*               87654321 */
#define LEFT_0 0x00000000
#define LEFT_1 0x40000000
#define LEFT_2 0x80000000
#define LEFT_3 0xC0000000


/* Genomecomp format */
/* A = 000  Stored as 00 in first two bytes, 0 in flag byte
   C = 001            01                     0
   G = 010            10                     0
   T = 011            11                     0
   N = 100            00                     1
   X = 111            11                     1
*/


/*                   87654321 */
#define LEFT_SET   0x80000000
#define LEFT_CLEAR 0x00000000


/* Genome128 format */
/*          High bit   Low bit   Flag bit
   A = 000     0          0         0
   C = 001     0          1         0
   G = 010     1          0         0
   T = 011     1          1         1
   N = 100     0          0         1
   X = 111     1          1         1
*/




/************************************************************************/

static void
genomecomp_move_absolute (FILE *fp, Univcoord_T ptr) {
#ifdef HAVE_FSEEKO
  off_t offset = ptr*((off_t) sizeof(Genomecomp_T));

  if (fseeko(fp,offset,SEEK_SET) < 0) {
    perror("Error in gmapindex, genomecomp_move_absolute");
    exit(9);
  }
#else
  long int offset = ptr*((long int) sizeof(Genomecomp_T));

  if (fseek(fp,offset,SEEK_SET) < 0) {
    perror("Error in gmapindex, genomecomp_move_absolute");
    exit(9);
  }
#endif

  return;
}

static void
genomecomp_read_current (Genomecomp_T *high, Genomecomp_T *low, Genomecomp_T *flags, FILE *fp,
			 int index1part) {
  char section[15];

  if (fread(section,sizeof(char),index1part,fp) < (unsigned int) index1part) {
    *high = 0xFFFFFFFF;
    *low = 0xFFFFFFFF;
    *flags = 0xFFFFFFFF;
    return;
  }

  *high = (section[3] & 0xff);
  *high <<= 8;
  *high |= (section[2] & 0xff);
  *high <<= 8;
  *high |= (section[1] & 0xff);
  *high <<= 8;
  *high |= (section[0] & 0xff);

  *low = (section[7] & 0xff);
  *low <<= 8;
  *low |= (section[6] & 0xff);
  *low <<= 8;
  *low |= (section[5] & 0xff);
  *low <<= 8;
  *low |= (section[4] & 0xff);

  *flags = (section[11] & 0xff);
  *flags <<= 8;
  *flags |= (section[10] & 0xff);
  *flags <<= 8;
  *flags |= (section[9] & 0xff);
  *flags <<= 8;
  *flags |= (section[8] & 0xff);

  return;
}


static void
write_compressed_one (FILE *fp, int *nbadchars, char Buffer[], Univcoord_T position) {
  Genomecomp_T high = 0U, low = 0U, flags = 0U, carry;
  int i;

  for (i = 0; i < 32; i++) {
    carry = high & 3U;
    high >>= 2;
    low >>= 2;
    flags >>= 1;
    switch (carry) {
    case 0U: break;
    case 1U: low |= LEFT_C; break;
    case 2U: low |= LEFT_G; break;
    case 3U: low |= LEFT_T; break;
    default: abort();
    }

    switch (uppercaseCode[(int) Buffer[i]]) {
    case 'A': break;
    case 'C': high |= LEFT_C; break;
    case 'G': high |= LEFT_G; break;
    case 'T': high |= LEFT_T; break;
    case 'N': flags |= LEFT_BIT; break;
    case 'X': high |= LEFT_T; flags |= LEFT_BIT; break;
    default: 
      (*nbadchars) += 1;
      if (*nbadchars < MAX_BADCHAR_MESSAGES) {
	fprintf(stderr,"Don't recognize character %c at position %lu.  Using N instead\n",
		Buffer[i],position+i);
      } else if (*nbadchars == MAX_BADCHAR_MESSAGES) {
	fprintf(stderr,"Too many non-recognizable characters.  Not reporting each individual occurrence anymore.\n");
      } else if ((*nbadchars) % BADCHAR_INTERVAL == 0) {
	fprintf(stderr,"A total of %d non-ACGTNX characters seen so far.\n",*nbadchars);
      }
      flags |= LEFT_BIT;
      break;
    }
  }
  
  FWRITE_UINT(high,fp);
  FWRITE_UINT(low,fp);
  FWRITE_UINT(flags,fp);
  
  return;
}


static void
put_compressed_one (Genomecomp_T *sectioncomp, int *nbadchars, char Buffer[], Univcoord_T position) {
  Genomecomp_T high = 0U, low = 0U, flags = 0U, carry;
  int i;

  for (i = 0; i < 32; i++) {
    carry = high & 3U;
    high >>= 2;
    low >>= 2;
    flags >>= 1;
    switch (carry) {
    case 0U: break;
    case 1U: low |= LEFT_C; break;
    case 2U: low |= LEFT_G; break;
    case 3U: low |= LEFT_T; break;
    default: abort();
    }

    switch (uppercaseCode[(int) Buffer[i]]) {
    case 'A': break;
    case 'C': high |= LEFT_C; break;
    case 'G': high |= LEFT_G; break;
    case 'T': high |= LEFT_T; break;
    case 'N': flags |= LEFT_BIT; break;
    case 'X': high |= LEFT_T; flags |= LEFT_BIT; break;
    default: 
      (*nbadchars) += 1;
      if (*nbadchars < MAX_BADCHAR_MESSAGES) {
	fprintf(stderr,"Don't recognize character %c at position %lu.  Using N instead\n",
		Buffer[i],position+i);
      } else if (*nbadchars == MAX_BADCHAR_MESSAGES) {
	fprintf(stderr,"Too many non-recognizable characters.  Not reporting each individual occurrence anymore.\n");
      } else if ((*nbadchars) % BADCHAR_INTERVAL == 0) {
	fprintf(stderr,"A total of %d non-ACGTNX characters seen so far.\n",*nbadchars);
      }
      flags |= LEFT_BIT;
      break;
    }
  }
  
  sectioncomp[0] = high;
  sectioncomp[1] = low;
  sectioncomp[2] = flags;
  
  return;
}


static char acgt[4] = {'A','C','G','T'};
static char non_acgt[4] = {'N','?','?','X'};

/* if gbuffer is NULL, then we fill with X's */
/* Based on genomecomp.  Version for genome128 not implemented yet */
int
Compress_update_file (int nbadchars, FILE *fp, char *gbuffer, Univcoord_T startpos,
		      Univcoord_T endpos, int index1part) {
  /* Chrpos_T length = endpos - startpos; */
  Univcoord_T startblock, endblock, ptr;
  unsigned int startdiscard, enddiscard, i;
  Genomecomp_T high, low, flags;
  char Buffer[32];


  ptr = startblock = startpos/32U*3;
  endblock = endpos/32U*3;
  startdiscard = startpos % 32;
  enddiscard = endpos % 32;
  
  if (endblock == startblock) {
    /* Special case */
    genomecomp_move_absolute(fp,ptr);
    genomecomp_read_current(&high,&low,&flags,fp,index1part);

    for (i = 0; i < 16; i++) {
      Buffer[i] = flags & 1U ? non_acgt[low & 3U] : acgt[low & 3U];
      low >>= 2;
      flags >>= 1;
    }
    for ( ; i < 32; i++) {
      Buffer[i] = flags & 1U ? non_acgt[high & 3U] : acgt[high & 3U];
      high >>= 2;
      flags >>= 1;
    }
    for (i = startdiscard; i < enddiscard; i++) {
      Buffer[i] = gbuffer ? *gbuffer++ : 'X';
    }
    genomecomp_move_absolute(fp,ptr);
    write_compressed_one(fp,&nbadchars,Buffer,startpos);

  } else {

    genomecomp_move_absolute(fp,ptr);
    genomecomp_read_current(&high,&low,&flags,fp,index1part);

    for (i = 0; i < 16; i++) {
      Buffer[i] = flags & 1U ? non_acgt[low & 3U] : acgt[low & 3U];
      low >>= 2;
      flags >>= 1;
    }
    for ( ; i < 32; i++) {
      Buffer[i] = flags & 1U ? non_acgt[high & 3U] : acgt[high & 3U];
      high >>= 2;
      flags >>= 1;
    }
    for (i = startdiscard; i < 32; i++) {
      Buffer[i] = gbuffer ? *gbuffer++ : 'X';
    }
    genomecomp_move_absolute(fp,ptr);
    write_compressed_one(fp,&nbadchars,Buffer,startpos);
    ptr += 3;
      
    while (ptr < endblock) {
      for (i = 0; i < 32; i++) {
	Buffer[i] = gbuffer ? *gbuffer++ : 'X';
      }
      write_compressed_one(fp,&nbadchars,Buffer,ptr/3*32U);
      ptr += 3;
    }

    if (enddiscard > 0) {
      genomecomp_move_absolute(fp,ptr);
      genomecomp_read_current(&high,&low,&flags,fp,index1part);

      for (i = 0; i < 16; i++) {
	Buffer[i] = flags & 1U ? non_acgt[low & 3U] : acgt[low & 3U];
	low >>= 2;
	flags >>= 1;
      }
      for ( ; i < 32; i++) {
	Buffer[i] = flags & 1U ? non_acgt[high & 3U] : acgt[high & 3U];
	high >>= 2;
	flags >>= 1;
      }
      for (i = 0; i < enddiscard; i++) {
	Buffer[i] = gbuffer ? *gbuffer++ : 'X';
      }
      genomecomp_move_absolute(fp,ptr);
      write_compressed_one(fp,&nbadchars,Buffer,ptr/3*32U);
    }
  }

  return nbadchars;
}


int
Compress_update_memory (int nbadchars, Genomecomp_T *genomecomp, char *gbuffer, Univcoord_T startpos,
			Univcoord_T endpos) {
  /* Chrpos_T length = endpos - startpos; */
  Univcoord_T startblock, endblock, ptr;
  Genomecomp_T high, low, flags;
  char Buffer[32];
  unsigned int startdiscard, enddiscard, i;


  ptr = startblock = startpos/32U*3;
  endblock = endpos/32U*3;
  startdiscard = startpos % 32;
  enddiscard = endpos % 32;
  
  if (endblock == startblock) {
    /* Special case */
    high = genomecomp[ptr];
    low = genomecomp[ptr+1];
    flags = genomecomp[ptr+2];

    /* Fill Buffer with original contents */
    for (i = 0; i < 16; i++) {
      Buffer[i] = flags & 1U ? non_acgt[low & 3U] : acgt[low & 3U];
      low >>= 2;
      flags >>= 1;
    }
    for ( ; i < 32; i++) {
      Buffer[i] = flags & 1U ? non_acgt[high & 3U] : acgt[high & 3U];
      high >>= 2;
      flags >>= 1;
    }
    for (i = startdiscard; i < enddiscard; i++) {
      Buffer[i] = gbuffer ? *gbuffer++ : 'X';
    }
    put_compressed_one(&(genomecomp[ptr]),&nbadchars,Buffer,startpos);

  } else {

    high = genomecomp[ptr];
    low = genomecomp[ptr+1];
    flags = genomecomp[ptr+2];

    /* Fill Buffer with original contents */
    for (i = 0; i < 16; i++) {
      Buffer[i] = flags & 1U ? non_acgt[low & 3U] : acgt[low & 3U];
      low >>= 2;
      flags >>= 1;
    }
    for ( ; i < 32; i++) {
      Buffer[i] = flags & 1U ? non_acgt[high & 3U] : acgt[high & 3U];
      high >>= 2;
      flags >>= 1;
    }
    for (i = startdiscard; i < 32; i++) {
      Buffer[i] = gbuffer ? *gbuffer++ : 'X';
    }
    put_compressed_one(&(genomecomp[ptr]),&nbadchars,Buffer,startpos);
    ptr += 3;
      
    while (ptr < endblock) {
      for (i = 0; i < 32; i++) {
	Buffer[i] = gbuffer ? *gbuffer++ : 'X';
      }
      put_compressed_one(&(genomecomp[ptr]),&nbadchars,Buffer,ptr/3*32U);
      ptr += 3;
    }

    if (enddiscard > 0) {
      high = genomecomp[ptr];
      low = genomecomp[ptr+1];
      flags = genomecomp[ptr+2];

      /* Fill Buffer with original contents */
      for (i = 0; i < 16; i++) {
	Buffer[i] = flags & 1U ? non_acgt[low & 3U] : acgt[low & 3U];
	low >>= 2;
	flags >>= 1;
      }
      for ( ; i < 32; i++) {
	Buffer[i] = flags & 1U ? non_acgt[high & 3U] : acgt[high & 3U];
	high >>= 2;
	flags >>= 1;
      }
      for (i = 0; i < enddiscard; i++) {
	Buffer[i] = gbuffer ? *gbuffer++ : 'X';
      }
      put_compressed_one(&(genomecomp[ptr]),&nbadchars,Buffer,ptr/3*32U);
    }
  }

  return nbadchars;
}


#ifdef HAVE_64_BIT
static inline void
nt_unshuffle (UINT4 *highbits, UINT4 *lowbits, UINT4 high, UINT4 low) {
  UINT8 x, t;

  x = (UINT8) high;
  x <<= 32;
  x |= low;

  t = (x ^ (x >> 1))  & 0x2222222222222222;  x = x ^ t ^ (t << 1);
  t = (x ^ (x >> 2))  & 0x0C0C0C0C0C0C0C0C;  x = x ^ t ^ (t << 2);
  t = (x ^ (x >> 4))  & 0x00F000F000F000F0;  x = x ^ t ^ (t << 4);
  t = (x ^ (x >> 8))  & 0x0000FF000000FF00;  x = x ^ t ^ (t << 8);
  t = (x ^ (x >> 16)) & 0x00000000FFFF0000;  x = x ^ t ^ (t << 16);

  *highbits = (UINT4) (x >> 32);
  *lowbits = (UINT4) x;

  return;
}

#else

static inline void
nt_unshuffle (UINT4 *highbits, UINT4 *lowbits, UINT4 high, UINT4 low) {
  UINT4 t;

  /* unshuffle high */
  t = (high ^ (high >> 1)) & 0x22222222;  high = high ^ t ^ (t << 1);
  t = (high ^ (high >> 2)) & 0x0C0C0C0C;  high = high ^ t ^ (t << 2);
  t = (high ^ (high >> 4)) & 0x00F000F0;  high = high ^ t ^ (t << 4);
  t = (high ^ (high >> 8)) & 0x0000FF00;  high = high ^ t ^ (t << 8);

  /* unshuffle low */
  t = (low ^ (low >> 1)) & 0x22222222;  low = low ^ t ^ (t << 1);
  t = (low ^ (low >> 2)) & 0x0C0C0C0C;  low = low ^ t ^ (t << 2);
  t = (low ^ (low >> 4)) & 0x00F000F0;  low = low ^ t ^ (t << 4);
  t = (low ^ (low >> 8)) & 0x0000FF00;  low = low ^ t ^ (t << 8);

  *highbits = (high & 0xFFFF0000) | (low >> 16);
  *lowbits = (high << 16) | (low & 0x0000FFFF);

  return;
}
#endif


void
Compress_unshuffle (FILE *out, FILE *in) {
  Genomecomp_T high, low, flags;
  Genomecomp_T highbits, lowbits;

  while (FREAD_UINT(&high,in) > 0 &&
	 FREAD_UINT(&low,in) > 0 &&
	 FREAD_UINT(&flags,in) > 0) {
    nt_unshuffle(&highbits,&lowbits,high,low);
    FWRITE_UINT(highbits,out);
    FWRITE_UINT(lowbits,out);
    FWRITE_UINT(flags,out);
  }

  return;
}


/* Needed for user-provided segment in GMAP */
Genomecomp_T *
Compress_create_blocks_comp (char *genomicseg, Univcoord_T genomelength) {
  Genomecomp_T *genomecomp;
  size_t nuint4;

  nuint4 = ((genomelength + 31)/32U)*3;
  genomecomp = (Genomecomp_T *) CALLOC(nuint4+4,sizeof(Genomecomp_T));
  /* Add 4 because Oligoindex_hr procedures point to nextlow as ptr+4 */

  /* Creates X's at end */
  genomecomp[nuint4-3] = 0xFFFFFFFF;
  genomecomp[nuint4-2] = 0xFFFFFFFF;
  genomecomp[nuint4-1] = 0xFFFFFFFF;

  /* Plus extra 4 */
  genomecomp[nuint4]   = 0xFFFFFFFF;
  genomecomp[nuint4+1] = 0xFFFFFFFF;
  genomecomp[nuint4+2] = 0xFFFFFFFF;
  genomecomp[nuint4+3] = 0xFFFFFFFF;

  Compress_update_memory(/*nbadchars*/0,genomecomp,genomicseg,/*currposition*/0,genomelength);

  return genomecomp;
}


/* Needed for user-provided segment in GMAP */
Genomecomp_T *
Compress_create_blocks_bits (Genomecomp_T *genomecomp, Univcoord_T genomelength) {
  Genomecomp_T *genomebits, highbits, lowbits, high, low, flags;
  size_t nuint4, ptr;

  nuint4 = ((genomelength + 31)/32U)*3;
  genomebits = (Genomecomp_T *) CALLOC(nuint4+4,sizeof(Genomecomp_T));

  for (ptr = 0; ptr < nuint4; ptr += 3) {
    high = genomecomp[ptr];
    low = genomecomp[ptr+1];
    flags = genomecomp[ptr+2];
    
    nt_unshuffle(&highbits,&lowbits,high,low);
    genomebits[ptr] = highbits;
    genomebits[ptr+1] = lowbits;
    genomebits[ptr+2] = flags;
  }

  return genomebits;
}
