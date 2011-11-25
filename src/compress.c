static char rcsid[] = "$Id: compress.c 48791 2011-09-30 18:39:38Z twu $";
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
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>		/* For off_t */
#endif
#include "complement.h"
#include "mem.h"		/* For Compress_new */


/* Another MONITOR_INTERVAL is in indexdb.c */
#define MONITOR_INTERVAL 10000000 /* 10 million nt */

#define MAX_BADCHAR_MESSAGES 10
#define BADCHAR_INTERVAL 1000000


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* print_blocks */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif


static char uppercaseCode[128] = UPPERCASE_U2T;

/* We use int *, rather than char *, because we eventually return an int,
   and we see problems converting from char to int */
static void
fill_buffer (int *Buffer, UINT4 high, UINT4 low, UINT4 flags, Genomicpos_T position) {
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
	printf("Parsing error; saw non-ACGT flag plus %c at position %u\n",Buffer[i],position+i);
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
      default: abort();
      }

      switch (uppercaseCode[c]) {
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
      
      if (in_counter == 8 * (int) sizeof(Genomicpos_T)) {
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
    while (in_counter < 8 * (int) sizeof(Genomicpos_T)) {
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
genomecomp_read_current (UINT4 *high, UINT4 *low, UINT4 *flags, FILE *fp,
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
write_compressed_one (FILE *fp, int *nbadchars, char Buffer[], Genomicpos_T position) {
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
	fprintf(stderr,"Don't recognize character %c at position %u.  Using N instead\n",
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
put_compressed_one (UINT4 *sectioncomp, int *nbadchars, char Buffer[], Genomicpos_T position) {
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
	fprintf(stderr,"Don't recognize character %c at position %u.  Using N instead\n",
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
int
Compress_update_file (int nbadchars, FILE *fp, char *gbuffer, Genomicpos_T startpos,
		      Genomicpos_T endpos, int index1part) {
  /* Genomicpos_T length = endpos - startpos; */
  Genomicpos_T startblock, endblock, startdiscard, enddiscard, ptr;
  UINT4 high, low, flags;
  char Buffer[32];
  unsigned int i;
  int k = 0;

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
      Buffer[i] = gbuffer ? gbuffer[k++] : 'X';
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
      Buffer[i] = gbuffer ? gbuffer[k++] : 'X';
    }
    genomecomp_move_absolute(fp,ptr);
    write_compressed_one(fp,&nbadchars,Buffer,startpos);
    ptr += 3;
      
    while (ptr < endblock) {
      for (i = 0; i < 32; i++) {
	Buffer[i] = gbuffer ? gbuffer[k++] : 'X';
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
	Buffer[i] = gbuffer ? gbuffer[k++] : 'X';
      }
      genomecomp_move_absolute(fp,ptr);
      write_compressed_one(fp,&nbadchars,Buffer,ptr/3*32U);
    }
  }

  return nbadchars;
}


int
Compress_update_memory (int nbadchars, UINT4 *genomecomp, char *gbuffer, Genomicpos_T startpos,
			Genomicpos_T endpos) {
  /* Genomicpos_T length = endpos - startpos; */
  Genomicpos_T startblock, endblock, startdiscard, enddiscard, ptr;
  UINT4 high, low, flags;
  char Buffer[32];
  unsigned int i;
  int k = 0;

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
      Buffer[i] = gbuffer ? gbuffer[k++] : 'X';
    }
    put_compressed_one(&(genomecomp[ptr]),&nbadchars,Buffer,startpos);

  } else {

    high = genomecomp[ptr];
    low = genomecomp[ptr+1];
    flags = genomecomp[ptr+2];

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
      Buffer[i] = gbuffer ? gbuffer[k++] : 'X';
    }
    put_compressed_one(&(genomecomp[ptr]),&nbadchars,Buffer,startpos);
    ptr += 3;
      
    while (ptr < endblock) {
      for (i = 0; i < 32; i++) {
	Buffer[i] = gbuffer ? gbuffer[k++] : 'X';
      }
      put_compressed_one(&(genomecomp[ptr]),&nbadchars,Buffer,ptr/3*32U);
      ptr += 3;
    }

    if (enddiscard > 0) {
      high = genomecomp[ptr];
      low = genomecomp[ptr+1];
      flags = genomecomp[ptr+2];

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
	Buffer[i] = gbuffer ? gbuffer[k++] : 'X';
      }
      put_compressed_one(&(genomecomp[ptr]),&nbadchars,Buffer,ptr/3*32U);
    }
  }

  return nbadchars;
}



#define T Compress_T
struct T {
  UINT4 *blocks;
  int nblocks;
  UINT4 **shift_array;
  bool availp[32];
};



void
Compress_free (T *old) {
  if (*old) {
    FREE((*old)->shift_array[0]);
    FREE((*old)->shift_array);
    FREE((*old)->blocks);
    FREE(*old);
  }
  return;
}

void
Compress_print (T this) {
  int ptr = 0;

  while (ptr < this->nblocks*3) {
    printf("high: %08X  low: %08X  flags: %08X\n",
	   this->blocks[ptr],this->blocks[ptr+1],this->blocks[ptr+2]);
    ptr += 3;
  }
  printf("\n");
  return;
}


int
Compress_nblocks (T this) {
  return this->nblocks;
}


#ifdef DEBUG1
static void
print_blocks (UINT4 *blocks, int nblocks) {
  int ptr = 0;

  while (ptr < nblocks*3) {
    printf("high: %08X  low: %08X  flags: %08X\n",
	   blocks[ptr],blocks[ptr+1],blocks[ptr+2]);
    ptr += 3;
  }
  printf("\n");
  return;
}
#endif



T
Compress_new (char *gbuffer, int length, bool plusp) {
  T new = (T) MALLOC(sizeof(*new));
  UINT4 low = 0U, high = 0U, flags = 0U, carry;
  Genomicpos_T ptr;
  int position;
  int c, i;
  int in_counter = 0;

  new->nblocks = (length+31)/32U;

  new->blocks = (UINT4 *) CALLOC((new->nblocks+1)*3,sizeof(UINT4));

  /* Note that elements of shift_array do not have extra block at beginning */
  new->shift_array = (UINT4 **) CALLOC(32,sizeof(UINT4 *));
  new->shift_array[0] = (UINT4 *) CALLOC(32*(new->nblocks+1)*3,sizeof(UINT4));
  new->availp[0] = false;
  for (i = 1; i < 32; i++) {
    new->shift_array[i] = &(new->shift_array[i-1][(new->nblocks+1)*3]);
    new->availp[i] = false;
  }

  ptr = 0;
  if (plusp == true) {
    for (position = 0U; position < length; position++) {
      c = gbuffer[position];
      /* printf("char: %c\n",c); */
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
      default: abort();
      }

      /* Assume that gbuffer is upper case */
      switch /*(uppercaseCode[c])*/ (c) {
      case 'A': break;
      case 'C': high |= LEFT_C; break;
      case 'G': high |= LEFT_G; break;
      case 'T': high |= LEFT_T; break;
      default: flags |= LEFT_BIT; break;
      }
      /* printf("high: %08X  low: %08X  flags: %08X\n",high,low,flags); */
      
      if (in_counter == 8 * (int) sizeof(Genomicpos_T)) {
	new->blocks[ptr] = high;
	new->blocks[ptr+1] = low;
	new->blocks[ptr+2] = flags;
	ptr += 3;

	low = high = flags = 0U;
	in_counter = 0;
      }
    }
  } else {

    for (position = length-1; position >= 0; position--) {
      c = gbuffer[position];
      /* printf("char: %c\n",c); */
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
      default: abort();
      }

      /* Assume that gbuffer is upper case */
      switch /*(uppercaseCode[c])*/ (c) {
      case 'A': high |= LEFT_T; break;
      case 'C': high |= LEFT_G; break;
      case 'G': high |= LEFT_C; break;
      case 'T': break;
      default: flags |= LEFT_BIT; break;
      }
      /* printf("high: %08X  low: %08X  flags: %08X\n",high,low,flags); */
      
      if (in_counter == 8 * (int) sizeof(Genomicpos_T)) {
	new->blocks[ptr] = high;
	new->blocks[ptr+1] = low;
	new->blocks[ptr+2] = flags;
	ptr += 3;

	low = high = flags = 0U;
	in_counter = 0;
      }
    }
  }

  if (in_counter > 0) {
    while (in_counter < 8 * (int) sizeof(Genomicpos_T)) {
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
#if 0
      /* Don't put in X's */
      high |= LEFT_T; flags |= LEFT_BIT;
#endif
      in_counter++;
    }

    new->blocks[ptr] = high;
    new->blocks[ptr+1] = low;
    new->blocks[ptr+2] = flags;
  }

  debug1(printf("Compress_new, plusp %d\n",plusp));
  debug1(print_blocks(new->blocks,new->nblocks));
  debug1(printf("\n"));

  return new;
}

#if 0
T
Compress_dibase_new (char *gbuffer, Genomicpos_T length, bool plusp) {
  T new = (T) MALLOC(sizeof(*new));
  UINT4 low = 0U, high = 0U, flags = 0U, carry;
  Genomicpos_T ptr;
  int position;
  int c, i;
  int in_counter = 0;

  new->nblocks = (length+31)/32U;

  new->blocks = (UINT4 *) CALLOC((new->nblocks+1)*3,sizeof(UINT4));

  /* Note that elements of shift_array do not have extra block at beginning */
  new->shift_array = (UINT4 **) CALLOC(32,sizeof(UINT4 *));
  new->shift_array[0] = (UINT4 *) CALLOC(32*(new->nblocks+1)*3,sizeof(UINT4));
  new->availp[0] = false;
  for (i = 1; i < 32; i++) {
    new->shift_array[i] = &(new->shift_array[i-1][(new->nblocks+1)*3]);
    new->availp[i] = false;
  }

  ptr = 0;
  if (plusp == true) {
    for (position = 0U; position < length; position++) {
      c = gbuffer[position];
      /* printf("char: %c\n",c); */
      in_counter++;

      carry = high & 3U;
      high >>= 2;
      low >>= 2;
      flags >>= 1;
      switch (carry) {
      case 0U: break;
      case 1U: low |= LEFT_1; break;
      case 2U: low |= LEFT_2; break;
      case 3U: low |= LEFT_3; break;
      default: abort();
      }

      switch (c) {
      case '0': break;
      case '1': high |= LEFT_1; break;
      case '2': high |= LEFT_2; break;
      case '3': high |= LEFT_3; break;
      default: flags |= LEFT_BIT; break;
      }
      /* printf("high: %08X  low: %08X  flags: %08X\n",high,low,flags); */
      
      if (in_counter == 8 * (int) sizeof(Genomicpos_T)) {
	new->blocks[ptr] = high;
	new->blocks[ptr+1] = low;
	new->blocks[ptr+2] = flags;
	ptr += 3;

	low = high = flags = 0U;
	in_counter = 0;
      }
    }
  } else {

    for (position = length-1; position >= 0; position--) {
      c = gbuffer[position];
      /* printf("char: %c\n",c); */
      in_counter++;

      carry = high & 3U;
      high >>= 2;
      low >>= 2;
      flags >>= 1;
      switch (carry) {
      case 0U: break;
      case 1U: low |= LEFT_1; break;
      case 2U: low |= LEFT_2; break;
      case 3U: low |= LEFT_3; break;
      default: abort();
      }

      switch (c) {
      case '0': break;
      case '1': high |= LEFT_1; break;
      case '2': high |= LEFT_2; break;
      case '3': high |= LEFT_3; break;
      default: flags |= LEFT_BIT; break;
      }
      /* printf("high: %08X  low: %08X  flags: %08X\n",high,low,flags); */
      
      if (in_counter == 8 * (int) sizeof(Genomicpos_T)) {
	new->blocks[ptr] = high;
	new->blocks[ptr+1] = low;
	new->blocks[ptr+2] = flags;
	ptr += 3;

	low = high = flags = 0U;
	in_counter = 0;
      }
    }
  }

  if (in_counter > 0) {
    while (in_counter < 8 * (int) sizeof(Genomicpos_T)) {
      carry = high & 3U;
      high >>= 2;
      low >>= 2;
      flags >>= 1;
      switch (carry) {
      case 0U: break;
      case 1U: low |= LEFT_1; break;
      case 2U: low |= LEFT_2; break;
      case 3U: low |= LEFT_3; break;
      default: abort();
      }
#if 0
      /* Don't put in X's */
      high |= LEFT_T; flags |= LEFT_BIT;
#endif
      in_counter++;
    }

    new->blocks[ptr] = high;
    new->blocks[ptr+1] = low;
    new->blocks[ptr+2] = flags;
  }

  debug1(printf("Compress_new, plusp %d\n",plusp));
  debug1(print_blocks(new->blocks,new->nblocks));
  debug1(printf("\n"));

  return new;
}


static char start[4] = "ACGT";

T *
Compress_dibase_array_new (char *gbuffer, Genomicpos_T dibase_length, bool plusp) {
  T *array = (T *) CALLOC(4,sizeof(T)), new;
  UINT4 low, high, flags, carry;
  Genomicpos_T ptr;
  int position;
  int d, c, i, nt, lastnt;
  int in_counter;
  int nt_length = dibase_length + 1;

  for (d = 0; d < 4; d ++) {
    new = array[d] = (T) MALLOC(sizeof(*new));

    new->nblocks = (nt_length+31)/32U;

    new->blocks = (UINT4 *) CALLOC((new->nblocks+1)*3,sizeof(UINT4));

    /* Note that elements of shift_array do not have extra block at beginning */
    new->shift_array = (UINT4 **) CALLOC(32,sizeof(UINT4 *));
    new->shift_array[0] = (UINT4 *) CALLOC(32*(new->nblocks+1)*3,sizeof(UINT4));
    new->availp[0] = false;
    for (i = 1; i < 32; i++) {
      new->shift_array[i] = &(new->shift_array[i-1][(new->nblocks+1)*3]);
      new->availp[i] = false;
    }

    ptr = 0;
    in_counter = 0;
    low = high = flags = 0U;
    lastnt = start[d];

    /* Handle starting character */
    in_counter++;
    switch (lastnt) {
    case 'A': break;
    case 'C': high |= LEFT_C; break;
    case 'G': high |= LEFT_G; break;
    case 'T': high |= LEFT_T; break;
    }

    if (plusp == true) {
      for (position = 0U; position < dibase_length; position++) {
	c = gbuffer[position];
	/* printf("char: %c\n",c); */
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
	default: abort();
	}

	/* Assume that gbuffer is upper case */
	switch /*(uppercaseCode[c])*/ (c) {
	case '0':
	  switch (lastnt) {
	  case 'A': nt = 'A'; break;
	  case 'C': nt = 'C'; break;
	  case 'G': nt = 'G'; break;
	  case 'T': nt = 'T'; break;
	  }
	  break;
	case '1':
	  switch (lastnt) {
	  case 'A': nt = 'C'; break;
	  case 'C': nt = 'A'; break;
	  case 'G': nt = 'T'; break;
	  case 'T': nt = 'G'; break;
	  }
	  break;
	case '2':
	  switch (lastnt) {
	  case 'A': nt = 'G'; break;
	  case 'C': nt = 'T'; break;
	  case 'G': nt = 'A'; break;
	  case 'T': nt = 'C'; break;
	  }
	  break;
	case '3':
	  switch (lastnt) {
	  case 'A': nt = 'T'; break;
	  case 'C': nt = 'G'; break;
	  case 'G': nt = 'C'; break;
	  case 'T': nt = 'A'; break;
	  }
	  break;
	default: abort();
	}

	/* debug1(printf("lastnt %c + dibase %d => nt %c\n",lastnt,c,nt)); */

	switch (nt) {
	case 'A': break;
	case 'C': high |= LEFT_C; break;
	case 'G': high |= LEFT_G; break;
	case 'T': high |= LEFT_T; break;
	}

	/* printf("high: %08X  low: %08X  flags: %08X\n",high,low,flags); */
      
	if (in_counter == 8 * (int) sizeof(Genomicpos_T)) {
	  new->blocks[ptr] = high;
	  new->blocks[ptr+1] = low;
	  new->blocks[ptr+2] = flags;
	  ptr += 3;

	  low = high = flags = 0U;
	  in_counter = 0;
	}

	lastnt = nt;
      }

    } else {

      for (position = dibase_length-1; position >= 0; position--) {
	c = gbuffer[position];
	       
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
	default: abort();
	}

	/* Assume that gbuffer is upper case */
	switch /*(uppercaseCode[c])*/ (c) {
	case '0':
	  switch (lastnt) {
	  case 'A': nt = 'A'; break;
	  case 'C': nt = 'C'; break;
	  case 'G': nt = 'G'; break;
	  case 'T': nt = 'T'; break;
	  }
	  break;
	case '1':
	  switch (lastnt) {
	  case 'A': nt = 'C'; break;
	  case 'C': nt = 'A'; break;
	  case 'G': nt = 'T'; break;
	  case 'T': nt = 'G'; break;
	  }
	  break;
	case '2':
	  switch (lastnt) {
	  case 'A': nt = 'G'; break;
	  case 'C': nt = 'T'; break;
	  case 'G': nt = 'A'; break;
	  case 'T': nt = 'C'; break;
	  }
	  break;
	case '3':
	  switch (lastnt) {
	  case 'A': nt = 'T'; break;
	  case 'C': nt = 'G'; break;
	  case 'G': nt = 'C'; break;
	  case 'T': nt = 'A'; break;
	  }
	  break;
	default: abort();
	}

	/* debug1(printf("lastnt %c + dibase %d => nt %c\n",lastnt,c,nt)); */

	switch (nt) {
	case 'A': break;
	case 'C': high |= LEFT_C; break;
	case 'G': high |= LEFT_G; break;
	case 'T': high |= LEFT_T; break;
	}

	/* printf("high: %08X  low: %08X  flags: %08X\n",high,low,flags); */
      
	if (in_counter == 8 * (int) sizeof(Genomicpos_T)) {
	  new->blocks[ptr] = high;
	  new->blocks[ptr+1] = low;
	  new->blocks[ptr+2] = flags;
	  ptr += 3;

	  low = high = flags = 0U;
	  in_counter = 0;
	}

	lastnt = nt;
      }
    }

    if (in_counter > 0) {
      while (in_counter < 8 * (int) sizeof(Genomicpos_T)) {
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
#if 0
	/* Don't put in X's */
	high |= LEFT_T; flags |= LEFT_BIT;
#endif
	in_counter++;
      }

      new->blocks[ptr] = high;
      new->blocks[ptr+1] = low;
      new->blocks[ptr+2] = flags;
    }

    debug1(printf("Compress_dibase_new, plusp %d, startc %c\n",plusp,start[d]));
    debug1(print_blocks(new->blocks,new->nblocks));
    debug1(printf("\n"));
  }

  return array;
}
#endif


UINT4 *
Compress_shift (T this, int nshift) {
  UINT4 *shifted;
  int leftshift, rightshift, rightshift_flags;
  int ptr;

  debug(printf("nshift is %d\n",nshift));
  if (this->availp[nshift] == false) {

    shifted = this->shift_array[nshift];

    /* Shift flags */
    if (nshift == 0) {
      ptr = this->nblocks*3;
      while (ptr >= 0) {
	shifted[ptr+2] = this->blocks[ptr+2];
	shifted[ptr+1] = this->blocks[ptr+1];
	shifted[ptr] = this->blocks[ptr];
	ptr -= 3;
      }

    } else if (nshift < 16) {

      leftshift = nshift*2;
      rightshift = 32 - leftshift;
      rightshift_flags = 32 - nshift;

      ptr = this->nblocks*3;
      while (ptr > 0) {
	shifted[ptr+2] = (this->blocks[ptr+2] << nshift) | (this->blocks[ptr-1] >> rightshift_flags);
	debug(printf("Making flag at %d by combining %08X at %d => %08X and %08X at %d => %08X\n",
		     ptr+2,this->blocks[ptr+2],ptr+2,this->blocks[ptr+2] << nshift,
		     this->blocks[ptr-1],ptr-1,this->blocks[ptr-1] >> rightshift_flags));
	shifted[ptr+1] = (this->blocks[ptr+1] << leftshift) | (this->blocks[ptr-3] >> rightshift); 
	debug(printf("Making low at %d by combining %08X at %d => %08X and %08X at %d => %08X\n",
		     ptr+1,this->blocks[ptr+1],ptr+1,this->blocks[ptr+1] << leftshift,
		     this->blocks[ptr-3],ptr-3,this->blocks[ptr-3] >> rightshift));
	shifted[ptr] = (this->blocks[ptr] << leftshift) | (this->blocks[ptr+1] >> rightshift);
	debug(printf("Making high at %d by combining %08X at %d => %08X and %08X at %d => %08X\n",
		     ptr,this->blocks[ptr],ptr,this->blocks[ptr] << leftshift,
		     this->blocks[ptr+1],ptr+1,this->blocks[ptr+1] >> rightshift));
	ptr -= 3;
      }
      shifted[2] = this->blocks[2] << nshift;
      debug(printf("Making flag at %d by combining %08X at %d => %08X\n",
		   ptr+2,this->blocks[ptr+2],ptr+2,this->blocks[ptr+2] << nshift));
      shifted[1] = this->blocks[1] << leftshift;
      debug(printf("Making low at %d by combining %08X at %d => %08X\n",
		   ptr+1,this->blocks[ptr+1],ptr+1,this->blocks[ptr+1] << leftshift));
      shifted[0] = (this->blocks[0] << leftshift) | (this->blocks[1] >> rightshift);
      debug(printf("Making high at %d by combining %08X at %d => %08X and %08X at %d => %08X\n",
		   ptr,this->blocks[ptr],ptr,this->blocks[ptr] << leftshift,
		   this->blocks[ptr+1],ptr+1,this->blocks[ptr+1] >> rightshift));
      
    } else if (nshift == 16) {
      ptr = this->nblocks*3;
      while (ptr > 0) {
	shifted[ptr+2] = (this->blocks[ptr+2] << 16) | (this->blocks[ptr-1] >> 16);
	debug(printf("Making flag at %d by combining %08X at %d => %08X and %08X at %d => %08X\n",
		     ptr+2,this->blocks[ptr+2],ptr+2,this->blocks[ptr+2] << 16,
		     this->blocks[ptr-1],ptr-1,this->blocks[ptr-1] >> 16));
	shifted[ptr+1] = this->blocks[ptr-3];
	debug(printf("Making low at %d by copying %08X at %d\n",ptr+1,this->blocks[ptr-3],ptr-3));
	shifted[ptr] = this->blocks[ptr+1];
	debug(printf("Making high at %d by copying %08X at %d\n",ptr,this->blocks[ptr+1],ptr+1));
	ptr -= 3;
      }
      shifted[2] = this->blocks[2] << 16;
      debug(printf("Making flag at %d by combining %08X at %d => %08X\n",
		   ptr+2,this->blocks[ptr+2],ptr+2,this->blocks[ptr+2] << 16));
      shifted[1] = 0U;
      debug(printf("Making low at 1 by setting 0U\n"));
      shifted[0] = this->blocks[1];
      debug(printf("Making high at %d by copying %08X at %d\n",0,this->blocks[ptr+1],ptr+1));
      
    } else {
      leftshift = (nshift - 16) * 2;
      rightshift = 32 - leftshift;
      rightshift_flags = 32 - nshift;
      
      ptr = this->nblocks*3;
      while (ptr > 0) {
	shifted[ptr+2] = (this->blocks[ptr+2] << nshift) | (this->blocks[ptr-1] >> rightshift_flags);
	debug(printf("Making flag at %d by combining %08X at %d => %08X and %08X at %d => %08X\n",
		     ptr+2,this->blocks[ptr+2],ptr+2,this->blocks[ptr+2] << nshift,
		     this->blocks[ptr-1],ptr-1,this->blocks[ptr-1] >> rightshift_flags));
	shifted[ptr+1] = (this->blocks[ptr-3] << leftshift) | (this->blocks[ptr-2] >> rightshift);
	debug(printf("Making low at %d by combining %08X at %d => %08X and %08X at %d => %08X\n",
		     ptr+1,this->blocks[ptr-3],ptr-3,this->blocks[ptr-3] << leftshift,
		     this->blocks[ptr-2],ptr-2,this->blocks[ptr-2] >> rightshift));
	shifted[ptr] = (this->blocks[ptr+1] << leftshift) | (this->blocks[ptr-3] >> rightshift);
	debug(printf("Making high at %d by combining %08X at %d => %08X and %08X at %d => %08X\n",
		     ptr,this->blocks[ptr+1],ptr+1,this->blocks[ptr+1] << leftshift,
		     this->blocks[ptr-3],ptr-3,this->blocks[ptr-3] >> rightshift));
	ptr -= 3;
      }
      shifted[2] = this->blocks[2] << nshift;
      debug(printf("Making flag at %d by combining %08X at %d => %08X\n",
		   ptr+2,this->blocks[ptr+2],ptr+2,this->blocks[ptr+2] << nshift));
      shifted[1] = 0U;
      debug(printf("Making low at 1 by setting 0U\n"));
      shifted[0] = this->blocks[1] << leftshift;
      debug(printf("Making high at %d by combining %08X at %d => %08X\n",
		   0,this->blocks[ptr+1],ptr+1,this->blocks[ptr+1] << leftshift));
    }
    this->availp[nshift] = true;
  }
			     
  return this->shift_array[nshift];
}

