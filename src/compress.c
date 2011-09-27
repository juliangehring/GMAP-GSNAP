static char rcsid[] = "$Id: compress.c,v 1.13 2005/10/27 22:53:09 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "compress.h"

#include <stdlib.h>
#include <stddef.h>
#include <ctype.h>		/* For isalpha, toupper */
#ifdef WORDS_BIGENDIAN
#include "bigendian.h"
#else
#include "littleendian.h"
#endif

/* Another MONITOR_INTERVAL is in indexdb.c */
#define MONITOR_INTERVAL 10000000 /* 10 million nt */


/* We use int *, rather than char *, because we eventually return an int,
   and we see problems converting from char to int */
static void
fill_buffer (int *Buffer, UINT4 high, UINT4 low, UINT4 flags) {
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
	fprintf(stderr,"Parsing error\n");
	exit(9);
      }
    }
    flags >>= 1;
  }

  return;
}

int
Compress_get_char (FILE *sequence_fp, Genomicpos_T position, bool uncompressedp) {
  UINT4 high, low, flags;
  static int SAVEBUFFER[32];
  int ptr;

  if (uncompressedp == true) {
    return toupper(fgetc(sequence_fp));
  } else if ((ptr = position % 32) == 0) {
    if (FREAD_UINT(&high,sequence_fp) <= 0 ||
	FREAD_UINT(&low,sequence_fp) <= 0 ||
	FREAD_UINT(&flags,sequence_fp) <= 0) {
      return EOF;
    } else {
      fill_buffer(SAVEBUFFER,high,low,flags);
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

/* A = 000  Stored as 00 in first two bytes, 0 in flag byte
   C = 001            01                     0
   G = 010            10                     0
   T = 011            11                     0
   N = 100            00                     1
   X = 111            11                     1
*/

void
Compress_compress (FILE *fp) {
  UINT4 low = 0U, high = 0U, flags = 0U, carry;
  Genomicpos_T position = 0U;
  int c;
  int in_counter = 0;

  while ((c = fgetc(fp)) != EOF) {
    if (isalpha(c)) {
      in_counter++;

      carry = high & 3U;
      high >>= 2;
      low >>= 2;
      flags >>= 1;
      switch (carry) {
      case 0U: break;
      case 1U: low |= LEFT_C; break;
      case 2U: low |= LEFT_G; break;
      case 3U: low |= LEFT_T; break;
      }

      switch (toupper(c)) {
      case 'A': break;
      case 'C': high |= LEFT_C; break;
      case 'G': high |= LEFT_G; break;
      case 'T': high |= LEFT_T; break;
      case 'N': flags |= LEFT_BIT; break;
      case 'X': high |= LEFT_T; flags |= LEFT_BIT; break;
      default: 
	fprintf(stderr,"Non-standard nucleotide %c at position %u.  Using N instead\n",c,position);
	flags |= LEFT_BIT;
	break;
      }
      
      if (in_counter == 8*sizeof(Genomicpos_T)) {
	FWRITE_UINT(high,stdout);
	FWRITE_UINT(low,stdout);
	FWRITE_UINT(flags,stdout);

	low = high = flags = 0U;
	in_counter = 0;
      }
    }
    position++;
    if (position % MONITOR_INTERVAL == 0) {
      fprintf(stderr,"Compressing position %u\n",position);
    }
  }

  if (in_counter > 0) {
    while (in_counter < 8*sizeof(Genomicpos_T)) {
      carry = high & 3U;
      high >>= 2;
      low >>= 2;
      flags >>= 1;
      switch (carry) {
      case 0U: break;
      case 1U: low |= LEFT_C; break;
      case 2U: low |= LEFT_G; break;
      case 3U: low |= LEFT_T; break;
      }
      high |= LEFT_T; flags |= LEFT_BIT;
      in_counter++;
    }

    FWRITE_UINT(high,stdout);
    FWRITE_UINT(low,stdout);
    FWRITE_UINT(flags,stdout);
  }

  return;
}

void
Compress_uncompress (FILE *fp, int wraplength) {
  int c;
  Genomicpos_T position = 0U;

  if (wraplength <= 0) {
    while ((c = Compress_get_char(fp,position,/*uncompressedp*/false)) != EOF) {
      printf("%c",c);
      position++;

      if (position % MONITOR_INTERVAL == 0) {
	fprintf(stderr,"Uncompressing position %u\n",position);
      }
    }
  } else {
    while ((c = Compress_get_char(fp,position,/*uncompressedp*/false)) != EOF) {
      printf("%c",c);
      position++;
      if (position % wraplength == 0) {
	printf("\n");
      }
      if (position % MONITOR_INTERVAL == 0) {
	fprintf(stderr,"Uncompressing position %u\n",position);
      }
    }
    if (position % wraplength != 0) {
      printf("\n");
    }
  }

  return;
}


/************************************************************************/

static void
genomecomp_move_absolute (FILE *fp, Genomicpos_T ptr) {
#ifdef HAVE_FSEEKO
  off_t offset = ptr*((off_t) sizeof(UINT4));

  if (fseeko(fp,offset,SEEK_SET) < 0) {
    perror("Error in gmapindex, genomecomp_move_absolute");
    exit(9);
  }
#else
  long int offset = ptr*((long int) sizeof(UINT4));

  if (fseek(fp,offset,SEEK_SET) < 0) {
    perror("Error in gmapindex, genomecomp_move_absolute");
    exit(9);
  }
#endif

  return;
}

static void
genomecomp_read_current (UINT4 *high, UINT4 *low, UINT4 *flags, FILE *fp) {
  char section[12];

  if (fread(section,sizeof(char),12,fp) < 12) {
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
write_compressed_one (FILE *fp, char Buffer[], Genomicpos_T position) {
  UINT4 high = 0U, low = 0U, flags = 0U, carry;
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
    }

    switch (toupper((int) Buffer[i])) {
    case 'A': break;
    case 'C': high |= LEFT_C; break;
    case 'G': high |= LEFT_G; break;
    case 'T': high |= LEFT_T; break;
    case 'N': flags |= LEFT_BIT; break;
    case 'X': high |= LEFT_T; flags |= LEFT_BIT; break;
    default: 
      fprintf(stderr,"Non-standard nucleotide %c at position %u.  Using N instead\n",
	      Buffer[i],position+i);
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
put_compressed_one (UINT4 *sectioncomp, char Buffer[], Genomicpos_T position) {
  UINT4 high = 0U, low = 0U, flags = 0U, carry;
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
    }

    switch (toupper((int) Buffer[i])) {
    case 'A': break;
    case 'C': high |= LEFT_C; break;
    case 'G': high |= LEFT_G; break;
    case 'T': high |= LEFT_T; break;
    case 'N': flags |= LEFT_BIT; break;
    case 'X': high |= LEFT_T; flags |= LEFT_BIT; break;
    default: 
      fprintf(stderr,"Don't recognize character %c at position %u.  Using N instead\n",
	      Buffer[i],position+i);
      flags |= LEFT_BIT;
      break;
    }
  }
  
  sectioncomp[0] = high;
  sectioncomp[1] = low;
  sectioncomp[2] = flags;
  
  return;
}


static char translate[8] = {'A','C','G','T','N','?','?','X'};

/* if gbuffer is NULL, then we fill with X's */
void
Compress_update_file (FILE *fp, char *gbuffer, Genomicpos_T startpos,
		      Genomicpos_T endpos) {
  /* Genomicpos_T length = endpos - startpos; */
  Genomicpos_T startblock, endblock, startdiscard, enddiscard, ptr;
  UINT4 high, low, flags;
  char Buffer[32];
  int i, k = 0;

  ptr = startblock = startpos/32U*3;
  endblock = endpos/32U*3;
  startdiscard = startpos % 32;
  enddiscard = endpos % 32;
  
  if (endblock == startblock) {
    /* Special case */
    genomecomp_move_absolute(fp,ptr);
    genomecomp_read_current(&high,&low,&flags,fp);

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
      Buffer[i] = gbuffer ? gbuffer[k++] : 'X';
    }
    genomecomp_move_absolute(fp,ptr);
    write_compressed_one(fp,Buffer,startpos);

  } else {

    genomecomp_move_absolute(fp,ptr);
    genomecomp_read_current(&high,&low,&flags,fp);

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
      Buffer[i] = gbuffer ? gbuffer[k++] : 'X';
    }
    genomecomp_move_absolute(fp,ptr);
    write_compressed_one(fp,Buffer,startpos);
    ptr += 3;
      
    while (ptr < endblock) {
      for (i = 0; i < 32; i++) {
	Buffer[i] = gbuffer ? gbuffer[k++] : 'X';
      }
      write_compressed_one(fp,Buffer,ptr/3*32U);
      ptr += 3;
    }

    if (enddiscard > 0) {
      genomecomp_move_absolute(fp,ptr);
      genomecomp_read_current(&high,&low,&flags,fp);

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
	Buffer[i] = gbuffer ? gbuffer[k++] : 'X';
      }
      genomecomp_move_absolute(fp,ptr);
      write_compressed_one(fp,Buffer,ptr/3*32U);
    }
  }

  return;
}


void
Compress_update_memory (UINT4 *genomecomp, char *gbuffer, Genomicpos_T startpos,
			Genomicpos_T endpos) {
  /* Genomicpos_T length = endpos - startpos; */
  Genomicpos_T startblock, endblock, startdiscard, enddiscard, ptr;
  UINT4 high, low, flags;
  char Buffer[32];
  int i, k = 0;

  ptr = startblock = startpos/32U*3;
  endblock = endpos/32U*3;
  startdiscard = startpos % 32;
  enddiscard = endpos % 32;
  
  if (endblock == startblock) {
    /* Special case */
    high = genomecomp[ptr];
    low = genomecomp[ptr+1];
    flags = genomecomp[ptr+2];

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
      Buffer[i] = gbuffer ? gbuffer[k++] : 'X';
    }
    put_compressed_one(&(genomecomp[ptr]),Buffer,startpos);

  } else {

    high = genomecomp[ptr];
    low = genomecomp[ptr+1];
    flags = genomecomp[ptr+2];

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
      Buffer[i] = gbuffer ? gbuffer[k++] : 'X';
    }
    put_compressed_one(&(genomecomp[ptr]),Buffer,startpos);
    ptr += 3;
      
    while (ptr < endblock) {
      for (i = 0; i < 32; i++) {
	Buffer[i] = gbuffer ? gbuffer[k++] : 'X';
      }
      put_compressed_one(&(genomecomp[ptr]),Buffer,ptr/3*32U);
      ptr += 3;
    }

    if (enddiscard > 0) {
      high = genomecomp[ptr];
      low = genomecomp[ptr+1];
      flags = genomecomp[ptr+2];

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
	Buffer[i] = gbuffer ? gbuffer[k++] : 'X';
      }
      put_compressed_one(&(genomecomp[ptr]),Buffer,ptr/3*32U);
    }
  }

  return;
}


