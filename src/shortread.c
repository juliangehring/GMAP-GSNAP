static char rcsid[] = "$Id: shortread.c 108220 2013-09-17 18:51:00Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifndef HAVE_MEMCPY
# define memcpy(d,s,n) bcopy((s),(d),(n))
#endif
#ifndef HAVE_MEMMOVE
# define memmove(d,s,n) bcopy((s),(d),(n))
#endif

#include "shortread.h"

#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>		/* For rindex */
#include <ctype.h>		/* For iscntrl and isspace */

#ifdef HAVE_ZLIB
#include <zlib.h>
#define GZBUFFER_SIZE 131072
#endif

#define PAIRED_ADAPTER_NMISMATCHES_ALLOWED 1
#define PAIRED_ADAPTER_MINLENGTH 20 /* Must exceed 14 or stage1 will complain */

#define OVERLAP_NMISMATCHES_ALLOWED 1
#define OVERLAP_MINLENGTH 10

#include "assert.h"
#include "mem.h"
#include "complement.h"
#include "intlist.h"
#include "fopen.h"


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Pointers for first half and second half */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

/* Adapter stripping */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif



/***********************************************************************
 *    Definitions:
 *
 *   XXXXXXXXXX  TTTTTT ACGT ...... ACGT AAAAAA
 *   barcode-->  <-------- fulllength -------->
 *               ^contents
 *               ^quality
 *
 ************************************************************************/


#define T Shortread_T
struct T {
  char *acc;			/* Accession */
  char *restofheader;		/* Rest of header */
  bool filterp;			/* If true, then sequence should be skipped */
  bool invertedp;

  char *contents;		/* Original sequence, ends with '\0' */
  char *contents_alloc;		/* Allocation */
  int fulllength;		/* Full length (not including chopped sequence) */

  /* GSNAP-specific fields */
  char *contents_uc;		/* Original sequence, ends with '\0' */
  char *contents_uc_alloc;	/* Allocation */

  char *barcode;
  int barcode_length;

  int overlap;

  char *chop;
  char *chop_quality;
  int choplength;

  char *quality;		/* For Illumina short reads read via FASTQ */
  char *quality_alloc;		/* Allocation */

  /* bool free_contents_p; */
};

char *
Shortread_accession (T this) {
  return this->acc;
}

char *
Shortread_header (T this) {
  return this->restofheader;
}

bool
Shortread_filterp (T this) {
  return this->filterp;
}

bool
Shortread_invertedp (T this) {
  return this->invertedp;
}

char *
Shortread_fullpointer (T this) {
  return this->contents;
}

char *
Shortread_fullpointer_uc (T this) {
  return this->contents_uc;
}

char *
Shortread_contents_uc (T this) {
  return this->contents_uc;
}

int
Shortread_barcode_length (T this) {
  return this->barcode_length;
}

char *
Shortread_barcode (T this) {
  return this->barcode;
}


int
Shortread_choplength (T this) {
  return this->choplength;
}

char *
Shortread_quality_string (T this) {
  return this->quality;
}

int
Shortread_fulllength (T this) {
  return this->fulllength;
}


void
Shortread_free (T *old) {

  if (*old) {
    if ((*old)->restofheader != NULL) {
      FREE_IN((*old)->restofheader);
    }
    if ((*old)->acc != NULL) {
      FREE_IN((*old)->acc);
    }

    FREE_IN((*old)->contents_alloc);

    FREE_IN((*old)->barcode);
    FREE_IN((*old)->chop_quality);
    FREE_IN((*old)->chop);
    FREE_IN((*old)->quality_alloc);
    FREE_IN((*old)->contents_uc_alloc);

    FREE_IN(*old);
  }

  return;
}


#define HEADERLEN 512
#define DISCARDLEN 8192

static char Header[HEADERLEN];
static char Discard[DISCARDLEN];

static char Read1[MAX_READLENGTH+1];
static char Read2[MAX_READLENGTH+1];
static char Quality[MAX_READLENGTH+1];


/* The first element of Sequence is always the null character, to mark
   the end of the string */

/* Skipping of dashes might still be buggy */
/*
#define DASH '-'
*/


/* Returns '>' if FASTA file, first sequence char if not */
int
Shortread_input_init (FILE *fp) {
  int c;
  bool okayp = false;

  Header[0] = '\0';

  while (okayp == false && (c = fgetc(fp)) != EOF) {
    debug(printf("Read character %c\n",c));
    if (iscntrl(c)) {
#ifdef DASH
    } else if (c == DASH) {
#endif
    } else if (isspace(c)) {
    } else {
      okayp = true;
    }
  }

  debug(printf("Returning initial character %c\n",c));
  return c;
}


#ifdef HAVE_ZLIB
/* Returns '>' if FASTA file, first sequence char if not */
int
Shortread_input_init_gzip (gzFile fp) {
  int c;
  bool okayp = false;

  Header[0] = '\0';

  while (okayp == false && (c = gzgetc(fp)) != EOF) {
    debug(printf("Read character %c\n",c));
    if (iscntrl(c)) {
#ifdef DASH
    } else if (c == DASH) {
#endif
    } else if (isspace(c)) {
    } else {
      okayp = true;
    }
  }

  debug(printf("Returning initial character %c\n",c));
  return c;
}
#endif


#ifdef HAVE_BZLIB
/* Returns '>' if FASTA file, first sequence char if not */
int
Shortread_input_init_bzip2 (Bzip2_T fp) {
  int c;
  bool okayp = false;

  Header[0] = '\0';

  while (okayp == false && (c = bzgetc(fp)) != EOF) {
    debug(printf("Read character %c\n",c));
    if (iscntrl(c)) {
#ifdef DASH
    } else if (c == DASH) {
#endif
    } else if (isspace(c)) {
    } else {
      okayp = true;
    }
  }

  debug(printf("Returning initial character %c\n",c));
  return c;
}
#endif

static int acc_fieldi_start = 0;
static int acc_fieldi_end = 0;
static bool force_single_end_p = false;
static bool filter_chastity_p = true;
static bool allow_paired_end_mismatch_p = false;

void
Shortread_setup (int acc_fieldi_start_in, int acc_fieldi_end_in,
		 bool force_single_end_p_in, bool filter_chastity_p_in,
		 bool allow_paired_end_mismatch_p_in) {
  acc_fieldi_start = acc_fieldi_start_in;
  acc_fieldi_end = acc_fieldi_end_in;
  force_single_end_p = force_single_end_p_in;
  filter_chastity_p = filter_chastity_p_in;
  allow_paired_end_mismatch_p = allow_paired_end_mismatch_p_in;

  return;
}


static char *skipped_acc = "Skipped";

static char *
input_header (bool *filterp, char **restofheader, FILE *fp, bool skipp) {
  char *acc = NULL, *p, *q;
  size_t length;

  *filterp = false;

  if (feof(fp)) {
    return (char *) NULL;
  } else if (fgets(&(Header[0]),HEADERLEN,fp) == NULL) {
    /* File must terminate after > */
    return (char *) NULL;
  }

  if (Header[0] == '\n') {
    Header[0] = '\0';
  } else if ((p = rindex(&(Header[0]),'\n')) != NULL) {
    if (p[-1] == '\r') {
      p--;
    }
    *p = '\0';
  } else {
    /* Eliminate rest of header from input */
    while (fgets(&(Discard[0]),DISCARDLEN,fp) != NULL &&
	   rindex(&(Discard[0]),'\n') == NULL) ;
  }

  if (skipp == true) {
    return (char *) skipped_acc;
  } else {
    p = &(Header[0]);
    while (*p != '\0' && !isspace((int) *p)) {
      p++;
    }

    if (filter_chastity_p == true) {
      q = p;
      /* Expecting <read>:<is filtered>:<control number>:<index sequence>, e.g., 1:Y:0:CTTGTA */
      while (*q != '\0' && *q != ':') {
	q++;
      }
      if (*q != '\0') {
	q++;
	if (*q == 'Y') {
	  *filterp = true;
	}
      }
    }

    if (*p == '\0') {
      /* Accession only */
      length = (p - &(Header[0]))/sizeof(char);
      acc = (char *) CALLOC_IN(length+1,sizeof(char));
      strcpy(acc,Header);
      (*restofheader) = (char *) CALLOC_IN(1,sizeof(char));
      (*restofheader)[0] = '\0';
    } else {
      *p = '\0';
      length = (p - &(Header[0]))/sizeof(char);
      acc = (char *) CALLOC_IN(length+1,sizeof(char));
      strcpy(acc,Header);
      p++;
      *restofheader = (char *) CALLOC_IN(strlen(p)+1,sizeof(char));
      strcpy(*restofheader,p);
    }
  
    return acc;
  }
} 


#ifdef HAVE_ZLIB
static char *
input_header_gzip (bool *filterp, char **restofheader, gzFile fp, bool skipp) {
  char *acc = NULL, *p, *q;
  size_t length;

  *filterp = false;

  if (gzeof(fp)) {
    return NULL;
  } else if (gzgets(fp,&(Header[0]),HEADERLEN) == NULL) {
    /* File must terminate after > */
    return NULL;
  }

  if (Header[0] == '\n') {
    Header[0] = '\0';
  } else if ((p = rindex(&(Header[0]),'\n')) != NULL) {
    if (p[-1] == '\r') {
      p--;
    }
    *p = '\0';
  } else {
    /* Eliminate rest of header from input */
    while (gzgets(fp,&(Discard[0]),DISCARDLEN) != NULL &&
	   rindex(&(Discard[0]),'\n') == NULL) ;
  }

  if (skipp) {
    return (char *) skipped_acc;
  } else {
    p = &(Header[0]);
    while (*p != '\0' && !isspace((int) *p)) {
      p++;
    }

    if (filter_chastity_p == true) {
      q = p;
      /* Expecting <read>:<is filtered>:<control number>:<index sequence>, e.g., 1:Y:0:CTTGTA */
      while (*q != '\0' && *q != ':') {
	q++;
      }
      if (*q != '\0') {
	q++;
	if (*q == 'Y') {
	  *filterp = true;
	}
      }
    }

    if (*p == '\0') {
      /* Accession only */
      length = (p - &(Header[0]))/sizeof(char);
      acc = (char *) CALLOC_IN(length+1,sizeof(char));
      strcpy(acc,Header);
      *restofheader = (char *) CALLOC_IN(1,sizeof(char));
      (*restofheader)[0] = '\0';
    } else {
      *p = '\0';
      length = (p - &(Header[0]))/sizeof(char);
      acc = (char *) CALLOC_IN(length+1,sizeof(char));
      strcpy(acc,Header);
      p++;
      *restofheader = (char *) CALLOC_IN(strlen(p)+1,sizeof(char));
      strcpy(*restofheader,p);
    }

    return acc;
  }
} 
#endif


#ifdef HAVE_BZLIB
static char *
input_header_bzip2 (bool *filterp, char **restofheader, Bzip2_T fp, bool skipp) {
  char *acc = NULL, *p, *q;
  size_t length;

  *filterp = false;

  if (bzeof(fp)) {
    return NULL;
  } else if (bzgets(fp,&(Header[0]),HEADERLEN) == NULL) {
    /* File must terminate after > */
    return NULL;
  }

  if (Header[0] == '\n') {
    Header[0] = '\0';
  } else if ((p = rindex(&(Header[0]),'\n')) != NULL) {
    if (p[-1] == '\r') {
      p--;
    }
    *p = '\0';
  } else {
    /* Eliminate rest of header from input */
    while (bzgets(fp,&(Discard[0]),DISCARDLEN) != NULL &&
	   rindex(&(Discard[0]),'\n') == NULL) ;
  }

  if (skipp) {
    return (char *) skipped_acc;
  } else {
    p = &(Header[0]);
    while (*p != '\0' && !isspace((int) *p)) {
      p++;
    }

    if (filter_chastity_p == true) {
      q = p;
      /* Expecting <read>:<is filtered>:<control number>:<index sequence>, e.g., 1:Y:0:CTTGTA */
      while (*q != '\0' && *q != ':') {
	q++;
      }
      if (*q != '\0') {
	q++;
	if (*q == 'Y') {
	  *filterp = true;
	}
      }
    }

    if (*p == '\0') {
      /* Accession only */
      length = (p - &(Header[0]))/sizeof(char);
      acc = (char *) CALLOC_IN(length+1,sizeof(char));
      strcpy(acc,Header);
      *restofheader = (char *) CALLOC_IN(1,sizeof(char));
      (*restofheader)[0] = '\0';
    } else {
      *p = '\0';
      length = (p - &(Header[0]))/sizeof(char);
      acc = (char *) CALLOC_IN(length+1,sizeof(char));
      strcpy(acc,Header);
      p++;
      *restofheader = (char *) CALLOC_IN(strlen(p)+1,sizeof(char));
      strcpy(*restofheader,p);
    }

    return acc;
  }
} 
#endif


#if FILE_CONSISTENT
static bool stripp = true;
#endif

/* Deletes /\D1/ or /\D2 or 3/ endings. */
static void
strip_illumina_acc_ending (char *acc1, char *acc2) {
  char *p, *q;
  char slash1, slash2;

#if FILE_CONSISTENT
  if (stripp == true) {
#endif
    p = acc1;
    while (*p != '\0') {
      p++;
    }

    q = acc2;
    while (*q != '\0') {
      q++;
    }

    /* Handle old style Illumina data that ends in ":0" or ":1" */
    if (p[-2] == ':' &&	q[-2] == ':' && p[-1] == q[-1]) {
      p -= 2;
      q -= 2;
    }

    /* Delete "/1" or "/2 or /3" endings */
    slash1 = p[-2];
    slash2 = q[-2];

    if (slash1 != slash2 || isdigit(slash1)) {
#if FILE_CONSISTENT
      /* fprintf(stderr,"Do not see /1 or /2 endings in header fields %s and %s.  Will no longer look for them.\n",acc1,acc2); */
      stripp = false;
#endif

    } else if (p[-1] != '1' || (q[-1] != '2' && q[-1] != '3')) {
#if FILE_CONSISTENT
      /* fprintf(stderr,"Do not see /1 or /2 endings in header fields %s and %s.  Will no longer look for them.\n",acc1,acc2); */
      stripp = false;
#endif

    } else {
      p[-2] = '\0';
      q[-2] = '\0';
    }
#if FILE_CONSISTENT
  }
#endif

  return;
}


static char *
input_header_fastq (bool *filterp, char **restofheader, FILE *fp, bool skipp) {
  char *acc, *p, *q, *start;
  size_t length;
  int fieldi = 0;

  *filterp = false;

  if (feof(fp)) {
    return NULL;
  } else if (fgets(&(Header[0]),HEADERLEN,fp) == NULL) {
    /* File must terminate after > */
    return NULL;
  }

  if (Header[0] == '\n') {
    Header[0] = '\0';
  } else if ((p = rindex(&(Header[0]),'\n')) != NULL) {
    if (p[-1] == '\r') {
      p--;
    }
    *p = '\0';
  } else {
    /* Eliminate rest of header from input */
    while (fgets(&(Discard[0]),DISCARDLEN,fp) != NULL &&
	   rindex(&(Discard[0]),'\n') == NULL) ;
  }

  if (skipp == true) {
    return (char *) skipped_acc;
  } else {
    p = start = &(Header[0]);
    while (fieldi < acc_fieldi_start) {
      while (*p != '\0' && !isspace((int) *p)) {
	p++;
      }
      if (*p != '\0') {
	p++;
      }
      start = p;
      fieldi++;
    }

    while (fieldi < acc_fieldi_end) {
      while (*p != '\0' && !isspace((int) *p)) {
	p++;
      }
      if (*p != '\0') {
	p++;
      }
      fieldi++;
    }

    while (*p != '\0' && !isspace((int) *p)) {
      p++;
    }

    if (filter_chastity_p == true) {
      q = p;
      /* Expecting <read>:<is filtered>:<control number>:<index sequence>, e.g., 1:Y:0:CTTGTA */
      while (*q != '\0' && *q != ':') {
	q++;
      }
      if (*q != '\0') {
	q++;
	if (*q == 'Y') {
	  *filterp = true;
	}
      }
    }

    *p = '\0';

    length = (p - start)/sizeof(char);
    acc = (char *) CALLOC_IN(length+1,sizeof(char));
    strcpy(acc,start);
    *restofheader = (char *) CALLOC_IN(1,sizeof(char));
    (*restofheader)[0] = '\0';

    return acc;
  }
} 


#ifdef HAVE_ZLIB
static char *
input_header_fastq_gzip (bool *filterp, char **restofheader, gzFile fp, bool skipp) {
  char *acc, *p, *q, *start;
  size_t length;
  int fieldi = 0;

  *filterp = false;

  if (gzeof(fp)) {
    return NULL;
  } else if (gzgets(fp,&(Header[0]),HEADERLEN) == NULL) {
    /* File must terminate after > */
    return NULL;
  }

  if (Header[0] == '\n') {
    Header[0] = '\0';
  } else if ((p = rindex(&(Header[0]),'\n')) != NULL) {
    if (p[-1] == '\r') {
      p--;
    }
    *p = '\0';
  } else {
    /* Eliminate rest of header from input */
    while (gzgets(fp,&(Discard[0]),DISCARDLEN) != NULL &&
	   rindex(&(Discard[0]),'\n') == NULL) ;
  }

  if (skipp) {
    return (char *) skipped_acc;
  } else {
    p = start = &(Header[0]);
    while (fieldi < acc_fieldi_start) {
      while (*p != '\0' && !isspace((int) *p)) {
	p++;
      }
      if (*p != '\0') {
	p++;
      }
      start = p;
      fieldi++;
    }

    while (fieldi < acc_fieldi_end) {
      while (*p != '\0' && !isspace((int) *p)) {
	p++;
      }
      if (*p != '\0') {
	p++;
      }
      fieldi++;
    }

    while (*p != '\0' && !isspace((int) *p)) {
      p++;
    }

    if (filter_chastity_p == true) {
      q = p;
      /* Expecting <read>:<is filtered>:<control number>:<index sequence>, e.g., 1:Y:0:CTTGTA */
      while (*q != '\0' && *q != ':') {
	q++;
      }
      if (*q != '\0') {
	q++;
	if (*q == 'Y') {
	  *filterp = true;
	}
      }
    }

    *p = '\0';

    length = (p - start)/sizeof(char);
    acc = (char *) CALLOC_IN(length+1,sizeof(char));
    strcpy(acc,start);
    *restofheader = (char *) CALLOC_IN(1,sizeof(char));
    (*restofheader)[0] = '\0';

    return acc;
  }
} 
#endif


#ifdef HAVE_BZLIB
static char *
input_header_fastq_bzip2 (bool *filterp, char **restofheader, Bzip2_T fp, bool skipp) {
  char *acc, *p, *q, *start;
  size_t length;
  int fieldi = 0;

  *filterp = false;

  if (bzeof(fp)) {
    return NULL;
  } else if (bzgets(fp,&(Header[0]),HEADERLEN) == NULL) {
    /* File must terminate after > */
    return NULL;
  }

  if (Header[0] == '\n') {
    Header[0] = '\0';
  } else if ((p = rindex(&(Header[0]),'\n')) != NULL) {
    if (p[-1] == '\r') {
      p--;
    }
    *p = '\0';
  } else {
    /* Eliminate rest of header from input */
    while (bzgets(fp,&(Discard[0]),DISCARDLEN) != NULL &&
	   rindex(&(Discard[0]),'\n') == NULL) ;
  }

  if (skipp) {
    return (char *) skipped_acc;
  } else {
    p = start = &(Header[0]);
    while (fieldi < acc_fieldi_start) {
      while (*p != '\0' && !isspace((int) *p)) {
	p++;
      }
      if (*p != '\0') {
	p++;
      }
      start = p;
      fieldi++;
    }

    while (fieldi < acc_fieldi_end) {
      while (*p != '\0' && !isspace((int) *p)) {
	p++;
      }
      if (*p != '\0') {
	p++;
      }
      fieldi++;
    }

    while (*p != '\0' && !isspace((int) *p)) {
      p++;
    }

    if (filter_chastity_p == true) {
      q = p;
      /* Expecting <read>:<is filtered>:<control number>:<index sequence>, e.g., 1:Y:0:CTTGTA */
      while (*q != '\0' && *q != ':') {
	q++;
      }
      if (*q != '\0') {
	q++;
	if (*q == 'Y') {
	  *filterp = true;
	}
      }
    }

    *p = '\0';

    length = (p - start)/sizeof(char);
    acc = (char *) CALLOC_IN(length+1,sizeof(char));
    strcpy(acc,start);
    *restofheader = (char *) CALLOC_IN(1,sizeof(char));
    (*restofheader)[0] = '\0';

    return acc;
  }
} 
#endif


static bool
skip_header (FILE *fp) {

  if (feof(fp)) {
    return false;
  } else if (fgets(&(Header[0]),HEADERLEN,fp) == NULL) {
    /* File must terminate after > */
    return false;
  }

  if (rindex(&(Header[0]),'\n') == NULL) {
    /* Eliminate rest of header from input */
    while (fgets(&(Discard[0]),DISCARDLEN,fp) != NULL &&
	   rindex(&(Discard[0]),'\n') == NULL) ;
  }

  return true;
} 

#ifdef HAVE_ZLIB
static bool
skip_header_gzip (gzFile fp) {

  if (gzeof(fp)) {
    return false;
  } else if (gzgets(fp,&(Header[0]),HEADERLEN) == NULL) {
    /* File must terminate after > */
    return false;
  }

  if (rindex(&(Header[0]),'\n') == NULL) {
    /* Eliminate rest of header from input */
    while (gzgets(fp,&(Discard[0]),DISCARDLEN) != NULL &&
	   rindex(&(Discard[0]),'\n') == NULL) ;
  }

  return true;
} 
#endif

#ifdef HAVE_BZLIB
static bool
skip_header_bzip2 (Bzip2_T fp) {

  if (bzeof(fp)) {
    return false;
  } else if (bzgets(fp,&(Header[0]),HEADERLEN) == NULL) {
    /* File must terminate after > */
    return false;
  }

  if (rindex(&(Header[0]),'\n') == NULL) {
    /* Eliminate rest of header from input */
    while (bzgets(fp,&(Discard[0]),DISCARDLEN) != NULL &&
	   rindex(&(Discard[0]),'\n') == NULL) ;
  }

  return true;
} 
#endif



#if 0
#define CONTROLM 13		/* From PC */
#define SPACE 32

static char *
find_bad_char (char *line) {
  char *first, *p1, *p2;
#ifdef DASH
  char *p3;
#endif

  p1 = index(line,CONTROLM);
  p2 = index(line,SPACE);
  /* p3 = index(line,DASH); */

  if (p1 == NULL && p2 == NULL /* && p3 == NULL*/) {
    return NULL;
  } else {
    if (p1) {
      first = p1;
    }
    if (p2) {
      first = p2;
    }
    /*
    if (p3) {
      first = p3;
    }
    */
    if (p1 && p1 < first) {
      first = p1;
    }
    if (p2 && p2 < first) {
      first = p2;
    }
    /*
    if (p3 && p3 < first) {
      first = p3;
    }
    */
    return first;
  }
}
#endif



static int
input_oneline (int *nextchar, char *Start, FILE *fp, char *acc, bool possible_fasta_header_p) {
  int remainder;
  char *ptr, *p = NULL;

  debug(printf("Entering input_oneline with nextchar = %c\n",*nextchar));

  ptr = &(Start[0]);
  remainder = (&(Start[MAX_READLENGTH]) - ptr)/sizeof(char);
  if (*nextchar == EOF || (possible_fasta_header_p == true && (*nextchar == '>' || *nextchar == '+'))) {
    debug(printf("Returning 0\n"));
    return 0;
  } else {
    *ptr++ = (char) *nextchar;
    if ((p = fgets(ptr,remainder+1,fp)) == NULL) {
      /* NULL if file ends with a blank line */
      printf("Blank line. read %s.\n",ptr);
    } else {
      debug(printf("Read %s.\n",ptr));
#if 0
      if (pc_linefeeds_p == true) {
	/* Should not expect to see PC line feeds */
	while ((p = find_bad_char(ptr)) != NULL) {
	  /* Handle PC line feed ^M */
	  ptr = p++;
	  strlenp = strlen(p);
	  memmove(ptr,p,strlenp);
	  ptr[strlenp] = '\0';
	  debug(printf("Found control-M/space.  Did memmove of %d chars at %p to %p\n",strlenp,p,ptr));
	}
      }
#endif

      if (*ptr == '\n') {
	*ptr = '\0';
	debug(printf("Now string is %s.\n",ptr));
      } else if ((p = index(ptr,'\n')) != NULL) {
	if (p[-1] == '\r') {
	  p--;
	}
	*p = '\0';
	debug(printf("Now string is %s.\n",ptr));
      } else if (feof(fp)) {
	/* No line feed, but end of file.  Handle below. */
	debug(printf("End of file seen\n"));
      } else {
	/* No line feed, but not end of file */
	fprintf(stderr,"Line %s is too long for allocated buffer size of %d.\n",&(Start[0]),MAX_READLENGTH);
	fprintf(stderr,"Problem occurred at accession %s.  Aborting.\n",acc);
	exit(9);
      }
    }
    ptr += strlen(ptr);

    /* Peek at character after eoln */
    if (feof(fp)) {
      *nextchar = EOF;
    } else {
      while ((*nextchar = fgetc(fp)) != EOF && (*nextchar == '\r' || *nextchar == '\n' || isspace(*nextchar))) {
      }
    }

    debug(printf("Returning %ld\n",(ptr - &(Start[0]))/sizeof(char)));
    return (ptr - &(Start[0]))/sizeof(char);
  }
}



#ifdef HAVE_ZLIB
static int
input_oneline_gzip (int *nextchar, char *Start, gzFile fp, char *acc, bool possible_fasta_header_p) {
  int remainder;
  char *ptr, *p = NULL;

  debug(printf("Entering input_oneline with nextchar = %c\n",*nextchar));

  ptr = &(Start[0]);
  remainder = (&(Start[MAX_READLENGTH]) - ptr)/sizeof(char);
  if (*nextchar == EOF || (possible_fasta_header_p == true && (*nextchar == '>' || *nextchar == '+'))) {
    debug(printf("Returning 0\n"));
    return 0;
  } else {
    *ptr++ = (char) *nextchar;
    if ((p = gzgets(fp,ptr,remainder+1)) == NULL) {
      /* NULL if file ends with a blank line */
      printf("Blank line. read %s.\n",ptr);
    } else {
      debug(printf("Read %s.\n",ptr));
#if 0
      if (pc_linefeeds_p == true) {
	/* Should not expect to see PC line feeds */
	while ((p = find_bad_char(ptr)) != NULL) {
	  /* Handle PC line feed ^M */
	  ptr = p++;
	  strlenp = strlen(p);
	  memmove(ptr,p,strlenp);
	  ptr[strlenp] = '\0';
	  debug(printf("Found control-M/space.  Did memmove of %d chars at %p to %p\n",strlenp,p,ptr));
	}
      }
#endif

      if (*ptr == '\n') {
	*ptr = '\0';
	debug(printf("Now string is %s.\n",ptr));
      } else if ((p = index(ptr,'\n')) != NULL) {
	if (p[-1] == '\r') {
	  p--;
	}
	*p = '\0';
	debug(printf("Now string is %s.\n",ptr));
      } else if (gzeof(fp)) {
	/* No line feed, but end of file.  Handle below. */
	debug(printf("End of file seen\n"));
      } else {
	/* No line feed, but not end of file */
	fprintf(stderr,"Line %s is too long for allocated buffer size of %d.\n",&(Start[0]),MAX_READLENGTH);
	fprintf(stderr,"Problem occurred at accession %s.  Aborting.\n",acc);
	exit(9);
      }
    }
    ptr += strlen(ptr);

    /* Peek at character after eoln */
    if (gzeof(fp)) {
      *nextchar = EOF;
    } else {
      while ((*nextchar = gzgetc(fp)) != EOF && (*nextchar == '\r' || *nextchar == '\n' || isspace(*nextchar))) {
      }
    }

    debug(printf("Returning %ld\n",(ptr - &(Start[0]))/sizeof(char)));
    return (ptr - &(Start[0]))/sizeof(char);
  }
}
#endif

#ifdef HAVE_BZLIB
static int
input_oneline_bzip2 (int *nextchar, char *Start, Bzip2_T fp, char *acc, bool possible_fasta_header_p) {
  int remainder;
  char *ptr, *p = NULL;

  debug(printf("Entering input_oneline with nextchar = %c\n",*nextchar));

  ptr = &(Start[0]);
  remainder = (&(Start[MAX_READLENGTH]) - ptr)/sizeof(char);
  if (*nextchar == EOF || (possible_fasta_header_p == true && (*nextchar == '>' || *nextchar == '+'))) {
    debug(printf("Returning 0\n"));
    return 0;
  } else {
    *ptr++ = (char) *nextchar;
    if ((p = bzgets(fp,ptr,remainder+1)) == NULL) {
      /* NULL if file ends with a blank line */
      printf("Blank line. read %s.\n",ptr);
    } else {
      debug(printf("Read %s.\n",ptr));
#if 0
      if (pc_linefeeds_p == true) {
	/* Should not expect to see PC line feeds */
	while ((p = find_bad_char(ptr)) != NULL) {
	  /* Handle PC line feed ^M */
	  ptr = p++;
	  strlenp = strlen(p);
	  memmove(ptr,p,strlenp);
	  ptr[strlenp] = '\0';
	  debug(printf("Found control-M/space.  Did memmove of %d chars at %p to %p\n",strlenp,p,ptr));
	}
      }
#endif

      if (*ptr == '\n') {
	*ptr = '\0';
	debug(printf("Now string is %s.\n",ptr));
      } else if ((p = index(ptr,'\n')) != NULL) {
	if (p[-1] == '\r') {
	  p--;
	}
	*p = '\0';
	debug(printf("Now string is %s.\n",ptr));
      } else if (bzeof(fp)) {
	/* No line feed, but end of file.  Handle below. */
	debug(printf("End of file seen\n"));
      } else {
	/* No line feed, but not end of file */
	fprintf(stderr,"Line %s is too long for allocated buffer size of %d.\n",&(Start[0]),MAX_READLENGTH);
	fprintf(stderr,"Problem occurred at accession %s.  Aborting.\n",acc);
	exit(9);
      }
    }
    ptr += strlen(ptr);

    /* Peek at character after eoln */
    if (bzeof(fp)) {
      *nextchar = EOF;
    } else {
      while ((*nextchar = bzgetc(fp)) != EOF && (*nextchar == '\r' || *nextchar == '\n' || isspace(*nextchar))) {
      }
    }

    debug(printf("Returning %ld\n",(ptr - &(Start[0]))/sizeof(char)));
    return (ptr - &(Start[0]))/sizeof(char);
  }
}
#endif


static char complCode[128] = COMPLEMENT_LC;

static char *
make_complement (char *sequence, unsigned int length) {
  char *complement;
  int i, j;

  complement = (char *) CALLOC_IN(length+1,sizeof(char));
  for (i = length-1, j = 0; i >= 0; i--, j++) {
    complement[j] = complCode[(int) sequence[i]];
  }
  complement[length] = '\0';
  return complement;
}


#if 0
static char *
make_complement_uppercase (char *sequence, unsigned int length) {
  char *complement;
  char uppercaseCode[128] = UPPERCASE_U2T;
  int i, j;

  complement = (char *) CALLOC_IN(length+1,sizeof(char));
  for (i = length-1, j = 0; i >= 0; i--, j++) {
    complement[j] = uppercaseCode[(int) complCode[(int) sequence[i]]];
  }
  complement[length] = '\0';
  return complement;
}
#endif


static char *
make_reverse (char *sequence, unsigned int length) {
  char *reverse;
  int i, j;

  if (sequence == NULL) {
    return (char *) NULL;
  } else {
    reverse = (char *) CALLOC_IN(length+1,sizeof(char));
    for (i = length-1, j = 0; i >= 0; i--, j++) {
      reverse[j] = sequence[i];
    }
    reverse[length] = '\0';
    return reverse;
  }
}

static char *
make_uppercase (char *sequence, unsigned int length) {
  char *uppercase;
#ifdef PMAP
  char uppercaseCode[128] = UPPERCASE_STD;
#else
  char uppercaseCode[128] = UPPERCASE_U2T;
#endif
  unsigned int i;

  uppercase = (char *) CALLOC_IN(length+1,sizeof(char));
  for (i = 0; i < length; i++) {
    uppercase[i] = uppercaseCode[(int) sequence[i]];
  }
  uppercase[length] = '\0';
  return uppercase;
}


/************************************************************************
 *   Dynamic programming for paired-end short reads
 ************************************************************************/

#ifdef USE_DYNAMIC_PROGRAMMING

typedef char Direction_T;
#define VERT 4
#define HORIZ 2
#define DIAG 1
#define STOP 0

static int **matrix_ptrs;
static int *matrix_space;
static Direction_T **directions_ptrs;
static Direction_T *directions_space;
static int **jump_ptrs;
static int *jump_space;


#define MATCH 1
#define MISMATCH -3
#define END_OPEN -14
#define OPEN -10
#define EXTEND -3

static int **pairdistance_array;
static int *jump_penalty_array;

static void
pairdistance_init () {
  int j, ptr;
  int c, c1, c2;

  pairdistance_array = (int **) CALLOC(128,sizeof(int *));
  pairdistance_array[0] = (int *) CALLOC(128*128,sizeof(int));
  ptr = 0;
  for (j = 1; j < 128; j++) {
    ptr += 128;
    pairdistance_array[j] = &(pairdistance_array[0][ptr]);
  }

  for (c1 = 'A'; c1 <= 'Z'; c1++) {
    for (c2 = 'A'; c2 < 'Z'; c2++) {
      pairdistance_array[c1][c2] = MISMATCH;
    }
  }

  for (c = 'A'; c < 'Z'; c++) {
    pairdistance_array[c][c] = MATCH;
  }

  return;
}

static void
jump_penalty_init (int maxlength) {
  int length, penalty;

  jump_penalty_array = (int *) CALLOC(maxlength+1,sizeof(int));

  penalty = OPEN;
  for (length = 0; length <= maxlength; length++) {
    jump_penalty_array[length] = penalty;
    penalty += EXTEND;
  }

  return;
}



#if 0
/* Called by reads_store.c */
void
Shortread_dynprog_init (int maxlength) {

  matrix_ptrs = (int **) CALLOC(maxlength+1,sizeof(int *));
  matrix_space = (int *) CALLOC((maxlength+1)*(maxlength+1),sizeof(int));
  directions_ptrs = (Direction_T **) CALLOC(maxlength+1,sizeof(Direction_T *));
  directions_space = (Direction_T *) CALLOC((maxlength+1)*(maxlength+1),sizeof(Direction_T));
  jump_ptrs = (int **) CALLOC(maxlength+1,sizeof(int *));
  jump_space = (int *) CALLOC((maxlength+1)*(maxlength+1),sizeof(int));

  pairdistance_init();
  jump_penalty_init(maxlength);

  return;
}

void
Shortread_dynprog_term () {
  FREE(jump_penalty_array);
  FREE(pairdistance_array[0]);
  FREE(pairdistance_array);

  FREE(matrix_ptrs);
  FREE(matrix_space);
  FREE(directions_ptrs);
  FREE(directions_space);
  FREE(jump_ptrs);
  FREE(jump_space);
  return;
}
#endif


/* Makes a matrix of dimensions 0..length1 x 0..length2 inclusive */
static int **
Matrix_alloc (int length1, int length2, int **ptrs, int *space) {
  int **matrix, i;

  if (length1 < 0 || length2 < 0) {
    fprintf(stderr,"shortread: lengths are negative: %d %d\n",length1,length2);
    abort();
  }

  matrix = ptrs;
  matrix[0] = space;
  for (i = 1; i <= length1; i++) {
    matrix[i] = &(matrix[i-1][length2 + 1]);
  }

  /* Clear memory only around the band, but doesn't work with Gotoh P1
     and Q1 matrices */
  /*
  for (r = 0; r <= length1; r++) {
    if ((clo = r - lband - 1) >= 0) {
      matrix[r][clo] = 0;
    }
    if ((chigh = r + rband + 1) <= length2) {
      matrix[r][chigh] = 0;
    }
  }
  */
  memset((void *) space,0,(length1+1)*(length2+1)*sizeof(int));

  return matrix;
}


#ifdef DEBUG2
static void
Matrix_print (int **matrix, int length1, int length2, char *sequence1, char *sequence2,
	      bool revp) {
  int i, j;

  printf("  ");
  for (j = 0; j <= length2; ++j) {
    if (j == 0) {
      printf("    ");
    } else {
      printf("  %c ",revp ? sequence2[-j+1] : sequence2[j-1]);
    }
  }
  printf("\n");

  for (i = 0; i <= length1; ++i) {
    if (i == 0) {
      printf("  ");
    } else {
      printf("%c ",revp ? sequence1[-i+1] : sequence1[i-1]);
    }
    for (j = 0; j <= length2; ++j) {
      printf("%3d ",matrix[i][j]);
    }
    printf("\n");
  }
  printf("\n");

  return;
}
#endif


/* Makes a matrix of dimensions 0..length1 x 0..length2 inclusive */
static Direction_T **
Directions_alloc (int length1, int length2, Direction_T **ptrs, Direction_T *space) {
  Direction_T **directions;
  int i;

  directions = ptrs;
  directions[0] = space;
  for (i = 1; i <= length1; i++) {
    directions[i] = &(directions[i-1][length2 + 1]);
  }

  /* Clear memory only around the band, but may not work with Gotoh
     method */
  /*
    for (r = 0; r <= length1; r++) {
    if ((clo = r - lband - 1) >= 0) {
    directions[r][clo] = STOP;
    }
    if ((chigh = r + rband + 1) <= length2) {
    directions[r][chigh] = STOP;
    }
    }
  */
  memset((void *) space,0,(length1+1)*(length2+1)*sizeof(Direction_T));

  return directions;
}


#ifdef DEBUG2
static void
Directions_print (Direction_T **directions, int **jump, int length1, int length2, 
		  char *sequence1, char *sequence2, bool revp) {
  int i, j;
  char buffer[4];

  printf("  ");
  for (j = 0; j <= length2; ++j) {
    if (j == 0) {
      printf("    ");
    } else {
      printf("  %c ",revp ? sequence2[-j+1] : sequence2[j-1]);
    }
  }
  printf("\n");

  for (i = 0; i <= length1; ++i) {
    if (i == 0) {
      printf("  ");
    } else {
      printf("%c ",revp ? sequence1[-i+1] : sequence1[i-1]);
    }
    for (j = 0; j <= length2; ++j) {
      if (directions[i][j] == DIAG) {
	sprintf(buffer,"D%d",jump[i][j]);
      } else if (directions[i][j] == HORIZ) {
	sprintf(buffer,"H%d",jump[i][j]);
      } else if (directions[i][j] == VERT) {
	sprintf(buffer,"V%d",jump[i][j]);
      } else {
	sprintf(buffer,"S%d",0);
      }
      printf("%3s ",buffer);
    }
    printf("\n");
  }
  printf("\n");
  return;
}
#endif


static int **
compute_scores_lookup (Direction_T ***directions, int ***jump, char *sequence1, char *sequence2, 
		       int length1, int length2) {
  int **matrix;
  int r, c, r1, c1;
  int bestscore, score, bestjump, j;
  Direction_T bestdir;

  matrix = Matrix_alloc(length1,length2,matrix_ptrs,matrix_space);
  *directions = Directions_alloc(length1,length2,directions_ptrs,directions_space);
  *jump = Matrix_alloc(length1,length2,jump_ptrs,jump_space);

  matrix[0][0] = 0;
  (*directions)[0][0] = STOP;
  (*jump)[0][0] = 0;

  /* Row 0 initialization */
  for (c = 1; c <= length2; c++) {
    matrix[0][c] = END_OPEN;
    (*directions)[0][c] = HORIZ;
    (*jump)[0][c] = c;
  }

  /* Column 0 initialization */
  for (r = 1; r <= length1; r++) {
    matrix[r][0] = END_OPEN;
    (*directions)[r][0] = VERT;
    (*jump)[r][0] = r;
  }

  for (r = 1; r <= length1; r++) {
    /* na1 = sequence1[r-1]; */

    for (c = 1; c <= length2; c++) {
      /* na2 = sequence2[c-1]; */

      /* Diagonal case */
      bestscore = matrix[r-1][c-1] + pairdistance_array[(int) sequence1[r-1]][(int) sequence2[c-1]];
      bestdir = DIAG;
      bestjump = 1;
      
      /* Horizontal case */
      for (c1 = c-1, j = 1; c1 >= 1; c1--, j++) {
	if ((*directions)[r][c1] == DIAG) {
	  score = matrix[r][c1] + jump_penalty_array[j];
	  if (score > bestscore) {
	    bestscore = score;
	    bestdir = HORIZ;
	    bestjump = j;
	  }
	}
      }

      /* Vertical case */
      for (r1 = r-1, j = 1; r1 >= 1; r1--, j++) {
	if ((*directions)[r1][c] == DIAG) {
	  score = matrix[r1][c] + jump_penalty_array[j];
	  if (score > bestscore) {
	    bestscore = score;
	    bestdir = VERT;
	    bestjump = j;
	  }
	}
      }

      /*
	debug(printf("At %d,%d, scoreV = %d, scoreH = %d, scoreD = %d\n",
	r,c,scoreV,scoreH,scoreD));
      */
      
      /* Update */
      matrix[r][c] = bestscore;
      (*directions)[r][c] = bestdir;
      (*jump)[r][c] = bestjump;
    }
  }

  /* Final cell */
  bestscore = matrix[length1-1][length2-1] + pairdistance_array[(int) sequence1[length1-1]][(int) sequence2[length2-1]];
  bestdir = DIAG;
  bestjump = 1;

  for (c = 1; c < length2; c++) {
    if ((score = matrix[length1][c] + END_OPEN) > bestscore) {
      bestscore = score;
      bestdir = HORIZ;
      bestjump = length2 - c;
    }
  }

  for (r = 1; r < length1; r++) {
    if ((score = matrix[r][length2] + END_OPEN) > bestscore) {
      bestscore = score;
      bestdir = VERT;
      bestjump = length1 - r;
    }
  }

  /* Update */
  matrix[length1][length2] = bestscore;
  (*directions)[length1][length2] = bestdir;
  (*jump)[length1][length2] = bestjump;

  /*
  debug2(Matrix_print(P0,length1,length2));
  debug2(Matrix_print(P1,length1,length2));
  debug2(Matrix_print(Q0,length1,length2));
  debug2(Matrix_print(Q1,length1,length2));
  */
  debug2(Matrix_print(matrix,length1,length2,sequence1,sequence2,/*revp*/false));
  debug2(Directions_print(*directions,*jump,length1,length2,
			  sequence1,sequence2,/*revp*/false));

  return matrix;
}


static void
traceback (int *nmatches, int *nmismatches, int *chop1, int *chop2,
	   Direction_T **directions, int **jump, char *sequence1, char *sequence2,
	   int length1, int length2) {
  int lastjump;
  int r, c;
  Direction_T lastdir;

  if (directions[length1][length2] != VERT) {
    *chop1 = *chop2 = 0;
  } else {
    *chop1 = jump[length1][length2];

    r = length1;
    c = length2;
    while (directions[r][c] != STOP) {
      if ((lastdir = directions[r][c]) == DIAG) {
	if (sequence1[r-1] == sequence2[c-1]) {
	  (*nmatches)++;
	} else {
	  (*nmismatches)++;
	}

	lastjump = 1;
	r--; c--;
	
      } else if (lastdir == HORIZ) {
	lastjump = jump[r][c];
	c -= lastjump;
	
      } else if (directions[r][c] == VERT) {
	lastjump = jump[r][c];
	r -= lastjump;
	
      } else {
	abort();
      }
    }

    if (lastdir != HORIZ) {
      *chop1 = *chop2 = 0;
    } else {
      *chop2 = lastjump;
    }
  }

  return;
}



/* Modeled after Dynprog_single_gap */
void
Shortread_chop_primers_slow (T queryseq1, T queryseq2) {
  int finalscore;
  int **matrix, **jump;
  Direction_T **directions;
  int chop1, chop2;
  int nmatches, nmismatches;

  matrix = compute_scores_lookup(&directions,&jump,queryseq1->contents_uc,queryseq2->contents_uc,
				 queryseq1->fulllength,queryseq2->fulllength);
  finalscore = matrix[queryseq1->fulllength][queryseq2->fulllength];

  if (finalscore <= 0) {
    chop1 = chop2 = 0;
    debug2(nmatches = 0);
    debug2(nmismatches = 0);
  } else {
    nmatches = nmismatches = 0;
    traceback(&nmatches,&nmismatches,&chop1,&chop2,directions,jump,
	      queryseq1->contents_uc,queryseq2->contents_uc,
	      queryseq1->fulllength,queryseq2->fulllength);
    if (nmismatches > PAIRED_ADAPTER_NMISMATCHES_ALLOWED) {
      chop1 = chop2 = 0;
    } else {
      if (chop1 > 0) {
	queryseq1->choplength = chop1;

	queryseq1->chop = (char *) CALLOC_IN(chop1+1,sizeof(char));
	strncpy(queryseq1->chop,&(queryseq1->contents[queryseq1->fulllength - chop1]),chop1);
	queryseq1->contents[queryseq1->fulllength - chop1] = '\0';
	queryseq1->contents_uc[queryseq1->fulllength - chop1] = '\0';

	if (queryseq1->quality != NULL) {
	  queryseq1->chop_quality = (char *) CALLOC_IN(chop1+1,sizeof(char));
	  strncpy(queryseq1->chop_quality,&(queryseq1->quality[queryseq1->fulllength - chop1]),chop1);
	  queryseq1->quality[queryseq1->fulllength - chop1] = '\0';
	}

	queryseq1->fulllength -= chop1;
      }
      if (chop2 > 0) {
	queryseq2->choplength = chop2;

	queryseq2->chop = (char *) CALLOC_IN(chop2+1,sizeof(char));
	strncpy(queryseq2->chop,queryseq2->contents,chop2);
	queryseq2->contents = &(queryseq2->contents[chop2]);
	queryseq2->contents_uc = &(queryseq2->contents_uc[chop2]);

	if (queryseq2->quality != NULL) {
	  queryseq2->chop_quality = (char *) CALLOC_IN(chop2+1,sizeof(char));
	  strncpy(queryseq2->chop_quality,queryseq2->quality,chop2);
	  queryseq2->quality = &(queryseq2->quality[chop2]);
	}

	queryseq2->fulllength -= chop2;
      }
    }
  }

  debug2(printf("finalscore = %d, nmatches = %d, nmismatches = %d, chop1 = %d, chop2 = %d\n",
		finalscore,nmatches,nmismatches,chop1,chop2));

  return;
}

#endif


bool
Shortread_chop_primers (T queryseq1, T queryseq2) {
  bool choppedp = false;
  int chop1 = 0, chop2 = 0;
  int nmatches, nmismatches;
  int fulllength1 = queryseq1->fulllength;
  int fulllength2 = queryseq2->fulllength;
  int minlength;
  char *contents1 = queryseq1->contents_uc;
  char *contents2 = queryseq2->contents_uc;

  int best_score = 0, score;
  int jstart, j, i;

  if (fulllength1 < fulllength2) {
    minlength = fulllength1;
  } else {
    minlength = fulllength2;
  }

  debug2(printf("jstart must be < %d - %d and < %d - %d\n",
		fulllength2,PAIRED_ADAPTER_MINLENGTH,fulllength1,PAIRED_ADAPTER_MINLENGTH));
  for (jstart = 0; jstart < minlength - PAIRED_ADAPTER_MINLENGTH; jstart++) {
    debug2(printf("jstart = %d\n",jstart));
    i = 0;
    j = jstart;
    nmismatches = 0;
    while (nmismatches <= PAIRED_ADAPTER_NMISMATCHES_ALLOWED && i < minlength - jstart) {
      debug2(printf("At i = %d, j = %d, comparing %c and %c\n",i,j,contents1[i],contents2[j]));
      assert(i < fulllength1);
      assert(j < fulllength2);
      if (contents1[i] != contents2[j]) {
	nmismatches++;
      }
      i++;
      j++;
    }
    if (nmismatches <= PAIRED_ADAPTER_NMISMATCHES_ALLOWED) {
      nmatches = j - nmismatches;
      if ((score = nmatches*3 - nmismatches) > best_score) {
	best_score = score;
	chop1 = jstart;
	chop2 = jstart;
      }
    }
  }

  debug2(printf("chop1 = %d, chop2 = %d\n",chop1,chop2));

  if (chop1 > 0) {
    choppedp = true;
    queryseq1->choplength = chop1;
    queryseq1->chop = (char *) CALLOC_IN(chop1+1,sizeof(char));
    queryseq1->fulllength -= chop1;
    fulllength1 = queryseq1->fulllength;

    strncpy(queryseq1->chop,&(queryseq1->contents[fulllength1]),chop1);
    queryseq1->contents[fulllength1] = '\0';
    queryseq1->contents_uc[fulllength1] = '\0';
    
    if (queryseq1->quality != NULL) {
      queryseq1->chop_quality = (char *) CALLOC_IN(chop1+1,sizeof(char));
      strncpy(queryseq1->chop_quality,&(queryseq1->quality[fulllength1]),chop1);
      queryseq1->quality[fulllength1] = '\0';
    }

  }

  if (chop2 > 0) {
    choppedp = true;
    queryseq2->choplength = chop2;
    queryseq2->chop = (char *) CALLOC_IN(chop2+1,sizeof(char));
    strncpy(queryseq2->chop,queryseq2->contents,chop2);
    queryseq2->contents = &(queryseq2->contents[chop2]);
    queryseq2->contents_uc = &(queryseq2->contents_uc[chop2]);
    
    if (queryseq2->quality != NULL) {
      queryseq2->chop_quality = (char *) CALLOC_IN(chop2+1,sizeof(char));
      strncpy(queryseq2->chop_quality,queryseq2->quality,chop2);
      queryseq2->quality = &(queryseq2->quality[chop2]);
    }
    queryseq2->fulllength -= chop2;

  }

  return choppedp;
}


bool
Shortread_find_primers (T queryseq1, T queryseq2) {
  int chop1 = 0, chop2 = 0;
  int nmatches, nmismatches;
  int fulllength1 = queryseq1->fulllength;
  int fulllength2 = queryseq2->fulllength;
  int minlength;
  char *contents1 = queryseq1->contents_uc;
  char *contents2 = queryseq2->contents_uc;

  int best_score = 0, score;
  int jstart, j, i;

  if (fulllength1 < fulllength2) {
    minlength = fulllength1;
  } else {
    minlength = fulllength2;
  }

  debug2(printf("jstart must be < %d - %d and < %d - %d\n",
		fulllength2,PAIRED_ADAPTER_MINLENGTH,fulllength1,PAIRED_ADAPTER_MINLENGTH));
  for (jstart = 0; jstart < minlength - PAIRED_ADAPTER_MINLENGTH; jstart++) {
    debug2(printf("jstart = %d\n",jstart));
    i = 0;
    j = jstart;
    nmismatches = 0;
    while (nmismatches <= PAIRED_ADAPTER_NMISMATCHES_ALLOWED && i < minlength - jstart) {
      debug2(printf("At i = %d, j = %d, comparing %c and %c\n",i,j,contents1[i],contents2[j]));
      assert(i < fulllength1);
      assert(j < fulllength2);
      if (contents1[i] != contents2[j]) {
	nmismatches++;
      }
      i++;
      j++;
    }
    if (nmismatches <= PAIRED_ADAPTER_NMISMATCHES_ALLOWED) {
      nmatches = j - nmismatches;
      if ((score = nmatches*3 - nmismatches) > best_score) {
	best_score = score;
	chop1 = jstart;
	chop2 = jstart;
      }
    }
  }

  debug2(printf("chop1 = %d, chop2 = %d\n",chop1,chop2));

  if (chop1 > 0 || chop2 > 0) {
    return true;
  } else {
    return false;
  }
}



int
Shortread_max_overlap (T queryseq1, T queryseq2) {
  if (queryseq1->fulllength < queryseq2->fulllength) {
    return queryseq1->fulllength;
  } else {
    return queryseq2->fulllength;
  }
}


int
Shortread_find_overlap (T queryseq1, T queryseq2) {
  int overlap = 0;
  int nmatches, nmismatches;
  int fulllength1 = queryseq1->fulllength;
  int fulllength2 = queryseq2->fulllength;
  char *contents1 = queryseq1->contents_uc;
  char *contents2 = queryseq2->contents_uc;

  int best_score = 0, score;
  int istart, i, j;

  if (queryseq1->overlap >= 0) {
    return queryseq1->overlap;
  } else {
    debug2(printf("istart must be < %d - %d and < %d - %d\n",
		  fulllength1,OVERLAP_MINLENGTH,fulllength2,OVERLAP_MINLENGTH));
    for (istart = 0; istart < fulllength1 - OVERLAP_MINLENGTH && istart < fulllength2 - OVERLAP_MINLENGTH; istart++) {
      debug2(printf("istart = %d\n",istart));
      i = istart;
      j = 0;
      nmismatches = 0;
      while (nmismatches <= OVERLAP_NMISMATCHES_ALLOWED && i < fulllength1 - istart) {
	debug2(printf("At i = %d, j = %d, comparing %c and %c\n",i,j,contents1[i],contents2[j]));
	if (contents1[i] != contents2[j]) {
	  nmismatches++;
	}
	i++;
	j++;
      }
      if (nmismatches <= OVERLAP_NMISMATCHES_ALLOWED) {
	nmatches = j - nmismatches;
	if ((score = nmatches*3 - nmismatches) > best_score) {
	  best_score = score;
	  overlap = istart;
	}
      }
    }

    debug2(printf("overlap = %d\n",overlap));

    queryseq1->overlap = queryseq2->overlap = overlap;
    return overlap;
  }
}


/************************************************************************
 *   Short reads
 ************************************************************************/

#define SKIPPED 1

T
Shortread_new (char *acc, char *restofheader, bool filterp,
	       char *sequence, int sequence_length, char *quality, int quality_length,
	       int barcode_length, bool invertp, bool copy_acc_p, bool skipp) {
  T new;

  if (sequence_length == 0) {
    return (T) NULL;
  } else if (skipp == true) {
    return (T) SKIPPED;
  }

  new = (T) MALLOC_IN(sizeof(*new));

  if (acc == NULL) {
    new->acc = (char *) NULL;
  } else if (copy_acc_p == false) {
    new->acc = acc;
  } else {
    new->acc = (char *) CALLOC_IN(strlen(acc)+1,sizeof(char));
    strcpy(new->acc,acc);
  }

  if (restofheader == NULL) {
    new->restofheader = (char *) NULL;
  } else if (copy_acc_p == false) {
    new->restofheader = restofheader;
  } else {
    new->restofheader = (char *) CALLOC_IN(strlen(restofheader)+1,sizeof(char));
    strcpy(new->restofheader,restofheader);
  }

  new->filterp = filterp;
  new->invertedp = invertp;

  /* printf("Before: sequence %s, quality %s\n",sequence,quality); */

  new->barcode_length = barcode_length;
  if (barcode_length == 0) {
    new->barcode = (char *) NULL;
  } else {
    new->barcode = (char *) CALLOC_IN(barcode_length+1,sizeof(char));
    strncpy(new->barcode,&(sequence[0]),barcode_length);
  }

  sequence_length -= barcode_length;
  quality_length -= barcode_length;
  new->fulllength = sequence_length; /* After barcode_length was removed */

  sequence = &(sequence[barcode_length]);
  quality = &(quality[barcode_length]);

  /* printf("After: barcode %s, sequence %s, quality %s\n",new->barcode,sequence,quality); */

  if (invertp == true) {
    new->contents = new->contents_alloc =  make_complement(sequence,sequence_length);
    new->contents_uc = new->contents_uc_alloc = 
      make_uppercase(new->contents,sequence_length);
    if (quality_length == 0) {
      new->quality = new->quality_alloc = (char *) NULL;
    } else {
      new->quality = new->quality_alloc = make_reverse(quality,quality_length);
    }

  } else {
    new->contents = new->contents_alloc = (char *) CALLOC_IN(sequence_length+1,sizeof(char));
    strncpy(new->contents,sequence,sequence_length);
    new->contents_uc = new->contents_uc_alloc =
      make_uppercase(new->contents,sequence_length);
    if (quality_length == 0) {
      new->quality = new->quality_alloc = (char *) NULL;
    } else {
      new->quality = new->quality_alloc = (char *) CALLOC_IN(quality_length+1,sizeof(char));
      strncpy(new->quality,quality,quality_length);
    }
  }

  new->overlap = -1;		/* Indicates not computed yet */

  new->chop = (char *) NULL;
  new->chop_quality = (char *) NULL;
  new->choplength = 0;

  return new;
}


T
Shortread_read_fasta_shortreads (int *nextchar, T *queryseq2, FILE **input1, FILE **input2,
				 char ***files, int *nfiles, bool skipp,
				 int barcode_length, bool invert_first_p, bool invert_second_p) {
  T queryseq1;
  char *acc, *restofheader, *acc2, *restofheader2;
  int nextchar2;
  int fulllength1, fulllength2, quality_length;
  bool filterp;

  while (1) {
    if (*input1 == NULL || feof(*input1)) {
      if (*input1 != NULL) {
	fclose(*input1);
	*input1 = NULL;
      }
      if (*input2 != NULL) {
	fclose(*input2);
	*input2 = NULL;
      }

      if (*nfiles == 0) {
	*nextchar = EOF;
	return (T) NULL;

      } else if (*nfiles == 1 || force_single_end_p == true) {
	if ((*input1 = FOPEN_READ_TEXT((*files)[0])) == NULL) {
	  fprintf(stderr,"Can't open file %s => skipping it.\n",(*files)[0]);
	  (*files) += 1;
	  (*nfiles) -= 1;
	  *nextchar = EOF;
	  return (T) NULL;
	} else {
	  *input2 = NULL;
	  (*files) += 1;
	  (*nfiles) -= 1;
	  *nextchar = '\0';
	}

      } else {
	while (*nfiles > 0 && (*input1 = FOPEN_READ_TEXT((*files)[0])) == NULL) {
	  fprintf(stderr,"Can't open file %s => skipping it.\n",(*files)[0]);
	  (*files) += 1;
	  (*nfiles) -= 1;
	}
	if (*input1 == NULL) {
	  *nextchar = EOF;
	  return (T) NULL;
	} else {
	  (*files) += 1;
	  (*nfiles) -= 1;
	  *nextchar = '\0';
	}
      }
    }

    if (*nextchar == '\0') {
      if ((*nextchar = Shortread_input_init(*input1)) == EOF) {
	*nextchar = EOF;
	return NULL;
      }
    }

    debug(printf("** Getting header\n"));
    if ((acc = input_header(&filterp,&restofheader,*input1,skipp)) == NULL) {
      /* fprintf(stderr,"No header\n"); */
      /* File ends after >.  Don't process, but loop again */
      *nextchar = EOF;
    } else if ((*nextchar = fgetc(*input1)) == '\r' || *nextchar == '\n') {
      /* Process blank lines and loop again */
      while (*nextchar != EOF && ((*nextchar = fgetc(*input1)) != '>')) {
      }
    } else if ((fulllength1 = input_oneline(&(*nextchar),&(Read1[0]),*input1,acc,
					    /*possible_fasta_header_p*/true)) == 0) {
      /* fprintf(stderr,"length is zero\n"); */
      /* No sequence1.  Don't process, but loop again */
      /* *nextchar = EOF; */
    } else {
      /* queryseq1 is in Read1 */
      /* See what is in next line */
      if ((fulllength2 = input_oneline(&(*nextchar),&(Read2[0]),*input1,acc,
				       /*possible_fasta_header_p*/true)) > 0) {
	/* Paired-end, single file.  queryseq1 is in Read1 and queryseq2 is in Read2 */
	if (*nextchar == '+') {
	  /* Paired-end with quality strings */
	  skip_header(*input1);
	  *nextchar = fgetc(*input1);
	  quality_length = input_oneline(&(*nextchar),&(Quality[0]),*input1,acc,
					 /*possible_fasta_header_p*/false);
	  if (quality_length != fulllength1) {
	    fprintf(stderr,"Length %d of quality score differs from length %d of nucleotides in sequence %s\n",
		    quality_length,fulllength1,acc);
	    abort();
	  }

	  queryseq1 = Shortread_new(acc,restofheader,filterp,Read1,fulllength1,
				    Quality,quality_length,barcode_length,
				    invert_first_p,/*copy_acc_p*/false,skipp);
	  quality_length = input_oneline(&(*nextchar),&(Quality[0]),*input1,acc,
					 /*possible_fasta_header_p*/false);
	  if (quality_length != fulllength2) {
	    fprintf(stderr,"Length %d of quality score differs from length %d of nucleotides in sequence %s\n",
		    quality_length,fulllength2,acc);
	    abort();
	  }
	      
	  (*queryseq2) = Shortread_new(/*acc*/NULL,/*restofheader*/NULL,filterp,Read2,fulllength2,
				       Quality,quality_length,barcode_length,
				       invert_second_p,/*copy_acc_p*/false,skipp);
	} else {
	  queryseq1 = Shortread_new(acc,restofheader,filterp,Read1,fulllength1,
				    /*quality*/NULL,/*quality_length*/0,barcode_length,
				    invert_first_p,/*copy_acc_p*/false,skipp);
	  (*queryseq2) = Shortread_new(/*acc*/NULL,/*restofheader*/NULL,filterp,Read2,fulllength2,
				       /*quality*/NULL,/*quality_length*/0,barcode_length,
				       invert_second_p,/*copy_acc_p*/false,skipp);
	}

      } else {
	if (*input2 == NULL && *nfiles > 0 && force_single_end_p == false &&
	    (*input2 = FOPEN_READ_TEXT((*files)[0])) != NULL) {
	  (*files) += 1;
	  (*nfiles) -= 1;
	  nextchar2 = '\0';
	}

	if (*input2 != NULL) {
	  /* Paired-end in two files */
	  if ((acc2 = input_header(&filterp,&restofheader2,*input2,skipp)) == NULL) {
	    /* fprintf(stderr,"No header\n"); */
	    /* File ends after >.  Don't process, but loop again */
	    (*queryseq2) = (T) NULL;
	    nextchar2 = EOF;
	  } else if ((nextchar2 = fgetc(*input2)) == '\r' || nextchar2 == '\n') {
	    /* Process blank lines and loop again */
	    while (nextchar2 != EOF && ((nextchar2 = fgetc(*input2)) != '>')) {
	    }
	    (*queryseq2) = (T) NULL;
	  } else if ((fulllength2 = input_oneline(&nextchar2,&(Read2[0]),*input2,acc2,
						  /*possible_fasta_header_p*/true)) == 0) {
	    /* fprintf(stderr,"length is zero\n"); */
	    /* No sequence1.  Don't process, but loop again */
	    /* nextchar2= EOF; */
	    (*queryseq2) = (T) NULL;
	  } else {
	    if (*nextchar == '+') {
	      /* End 1 with a quality string */
	      skip_header(*input1);
	      *nextchar = fgetc(*input1);
	      quality_length = input_oneline(&(*nextchar),&(Quality[0]),*input1,acc,
					     /*possible_fasta_header_p*/false);
	      if (quality_length != fulllength1) {
		fprintf(stderr,"Length %d of quality score differs from length %d of nucleotides in sequence %s\n",
			quality_length,fulllength1,acc);
		abort();
	      } else {
		queryseq1 = Shortread_new(acc,restofheader,filterp,Read1,fulllength1,
					  Quality,quality_length,barcode_length,
					  invert_first_p,/*copy_acc_p*/false,skipp);
	      }
	    } else {
	      /* End 1 without quality string */
	      queryseq1 = Shortread_new(acc,restofheader,filterp,Read1,fulllength1,
					/*quality*/NULL,/*quality_length*/0,barcode_length,
					invert_first_p,/*copy_acc_p*/false,skipp);
	    }

	    if (nextchar2 == '+') {
	      /* End 2 with a quality string */
	      skip_header(*input2);
	      nextchar2 = fgetc(*input2);
	      quality_length = input_oneline(&nextchar2,&(Quality[0]),*input2,acc2,
					     /*possible_fasta_header_p*/false);
	      if (quality_length != fulllength2) {
		fprintf(stderr,"Length %d of quality score differs from length %d of nucleotides in sequence %s\n",
			quality_length,fulllength2,acc2);
		abort();
	      } else {
		/* For FASTA, drop second accession */
		(*queryseq2) = Shortread_new(/*acc2*/NULL,/*restofheader2*/NULL,filterp,Read2,fulllength2,
					     Quality,quality_length,barcode_length,
					     invert_second_p,/*copy_acc_p*/false,skipp);
		FREE_IN(acc2);
		FREE_IN(restofheader2);
	      }
	    } else {
	      /* End 2 without quality string */
	      (*queryseq2) = Shortread_new(/*acc2*/NULL,/*restofheader2*/NULL,filterp,Read2,fulllength2,
					   /*quality*/NULL,/*quality_length*/0,barcode_length,
					   invert_second_p,/*copy_acc_p*/false,skipp);
	      FREE_IN(acc2);
	      FREE_IN(restofheader2);
	    }
	  }

	} else {
	  /* Single-end: Either EOF, '>', or '+' */
	  if (*nextchar == '+') {
	    /* Single-end with a quality string */
	    skip_header(*input1);
	    *nextchar = fgetc(*input1);
	    quality_length = input_oneline(&(*nextchar),&(Quality[0]),*input1,acc,
					   /*possible_fasta_header_p*/false);
	    if (quality_length != fulllength1) {
	      fprintf(stderr,"Length %d of quality score differs from length %d of nucleotides in sequence %s\n",
		      quality_length,fulllength1,acc);
	      abort();
	    } else {
	      queryseq1 = Shortread_new(acc,restofheader,filterp,Read1,fulllength1,
					Quality,quality_length,barcode_length,
					invert_first_p,/*copy_acc_p*/false,skipp);
	    }
	  } else {
	    /* Single-end without quality string */
	    queryseq1 = Shortread_new(acc,restofheader,filterp,Read1,fulllength1,
				      /*quality*/NULL,/*quality_length*/0,barcode_length,
				      invert_first_p,/*copy_acc_p*/false,skipp);
	  }
	  (*queryseq2) = (T) NULL;
	}
      }

      debug(printf("Returning queryseq with contents %s\n",queryseq1->contents));

      return queryseq1;
    }
  }
}




#ifdef HAVE_ZLIB
T
Shortread_read_fasta_shortreads_gzip (int *nextchar, T *queryseq2, gzFile *input1, gzFile*input2,
				      char ***files, int *nfiles, bool skipp,
				      int barcode_length, bool invert_first_p, bool invert_second_p) {
  T queryseq1;
  char *acc, *restofheader, *acc2, *restofheader2;
  int nextchar2;
  int fulllength1, fulllength2;
  bool filterp;

  while (1) {
    if (*input1 == NULL || gzeof(*input1)) {
      if (*input1 != NULL) {
	gzclose(*input1);
	*input1 = NULL;
      }
      if (*input2 != NULL) {
	gzclose(*input2);
	*input2 = NULL;
      }

      if (*nfiles == 0) {
	*nextchar = EOF;
	return (T) NULL;

      } else if (*nfiles == 1 || force_single_end_p == true) {
	if ((*input1 = gzopen((*files)[0],"rb")) == NULL) {
	  fprintf(stderr,"Can't open file %s => skipping it.\n",(*files)[0]);
	  (*files) += 1;
	  (*nfiles) -= 1;
	  *nextchar = EOF;
	  return (T) NULL;
	} else {
	  *input2 = NULL;
	  (*files) += 1;
	  (*nfiles) -= 1;
	  *nextchar = '\0';
	}

      } else {
	while (*nfiles > 0 && (*input1 = gzopen((*files)[0],"rb")) == NULL) {
	  fprintf(stderr,"Can't open file %s => skipping it.\n",(*files)[0]);
	  (*files)++;
	  (*nfiles)--;
	}
	if (*input1 == NULL) {
	  *nextchar = EOF;
	  return (T) NULL;
	} else {
#ifdef HAVE_ZLIB_GZBUFFER
	  gzbuffer(*input1,GZBUFFER_SIZE);
#endif
	  (*files)++;
	  (*nfiles)--;
	  *nextchar = '\0';
	}
      }
    }

    if (*nextchar == '\0') {
      if ((*nextchar = Shortread_input_init_gzip(*input1)) == EOF) {
	*nextchar = EOF;
	return NULL;
      }
    }

    debug(printf("** Getting header\n"));
    if ((acc = input_header_gzip(&filterp,&restofheader,*input1,skipp)) == NULL) {
      /* fprintf(stderr,"No header\n"); */
      /* File ends after >.  Don't process, but loop again */
      *nextchar = EOF;
    } else if ((*nextchar = gzgetc(*input1)) == '\r' || *nextchar == '\n') {
      /* Process blank lines and loop again */
      while (*nextchar != EOF && ((*nextchar = gzgetc(*input1)) != '>')) {
      }
    } else if ((fulllength1 = input_oneline_gzip(&(*nextchar),&(Read1[0]),*input1,acc,
						 /*possible_fasta_header_p*/true)) == 0) {
      /* fprintf(stderr,"length is zero\n"); */
      /* No sequence1.  Don't process, but loop again */
      /* *nextchar = EOF; */
    } else {
      /* queryseq1 is in Read1 */
      /* See what is in next line */
      if ((fulllength2 = input_oneline_gzip(&(*nextchar),&(Read2[0]),*input1,acc,
					    /*possible_fasta_header_p*/true)) > 0) {
	/* Paired-end, single file.  queryseq1 is in Read1 and queryseq2 is in Read2 */
	queryseq1 = Shortread_new(acc,restofheader,filterp,Read1,fulllength1,
				  /*quality*/NULL,/*quality_length*/0,barcode_length,
				  invert_first_p,/*copy_acc_p*/false,skipp);
	(*queryseq2) = Shortread_new(/*acc*/NULL,/*restofheader*/NULL,filterp,Read2,fulllength2,
				     /*quality*/NULL,/*quality_length*/0,barcode_length,
				     invert_second_p,/*copy_acc_p*/false,skipp);

      } else {
	if (*input2 == NULL && *nfiles > 0 && force_single_end_p == false &&
	    (*input2 = gzopen((*files)[0],"rb")) != NULL) {
#ifdef HAVE_ZLIB_GZBUFFER
	  gzbuffer(*input2,GZBUFFER_SIZE);
#endif
	  (*files) += 1;
	  (*nfiles) -= 1;
	  nextchar2 = '\0';
	}

	if (*input2 != NULL) {
	  /* Paired-end in two files */
	  if ((acc2 = input_header_gzip(&filterp,&restofheader2,*input2,skipp)) == NULL) {
	    /* fprintf(stderr,"No header\n"); */
	    /* File ends after >.  Don't process, but loop again */
	    (*queryseq2) = (T) NULL;
	    nextchar2 = EOF;
	  } else if ((nextchar2 = gzgetc(*input2)) == '\r' || nextchar2 == '\n') {
	    /* Process blank lines and loop again */
	    while (nextchar2 != EOF && ((nextchar2 = gzgetc(*input2)) != '>')) {
	    }
	    (*queryseq2) = (T) NULL;
	  } else if ((fulllength2 = input_oneline_gzip(&nextchar2,&(Read2[0]),*input2,acc2,
						       /*possible_fasta_header_p*/true)) == 0) {
	    /* fprintf(stderr,"length is zero\n"); */
	    /* No sequence1.  Don't process, but loop again */
	    /* *nextchar = EOF; */
	    (*queryseq2) = (T) NULL;
	  } else {
	    queryseq1 = Shortread_new(acc,restofheader,filterp,Read1,fulllength1,
				      /*quality*/NULL,/*quality_length*/0,barcode_length,
				      invert_first_p,/*copy_acc_p*/false,skipp);
	    (*queryseq2) = Shortread_new(/*acc2*/NULL,/*restofheader2*/NULL,filterp,Read2,fulllength2,
					 /*quality*/NULL,/*quality_length*/0,barcode_length,
					 invert_second_p,/*copy_acc_p*/false,skipp);
	    FREE_IN(acc2);
	    FREE_IN(restofheader2);
	  }

	} else {
	  /* Single-end: Either EOF, '>', or '+' */
	  queryseq1 = Shortread_new(acc,restofheader,filterp,Read1,fulllength1,
				    /*quality*/NULL,/*quality_length*/0,barcode_length,
				    invert_first_p,/*copy_acc_p*/false,skipp);
	  (*queryseq2) = (T) NULL;
	}
      }

      debug(printf("Returning queryseq with contents %s\n",queryseq1->contents));

      return queryseq1;
    }

  }
}
#endif


#ifdef HAVE_BZLIB
T
Shortread_read_fasta_shortreads_bzip2 (int *nextchar, T *queryseq2, Bzip2_T *input1, Bzip2_T *input2,
				       char ***files, int *nfiles, bool skipp,
				       int barcode_length, bool invert_first_p, bool invert_second_p) {
  T queryseq1;
  char *acc, *restofheader, *acc2, *restofheader2;
  int nextchar2;
  int fulllength1, fulllength2;
  bool filterp;

  while (1) {
    if (*input1 == NULL || bzeof(*input1)) {
      if (*input1 != NULL) {
	Bzip2_free(&(*input1));
	*input1 = NULL;
      }
      if (*input2 != NULL) {
	Bzip2_free(&(*input2));
	*input2 = NULL;
      }

      if (*nfiles == 0) {
	*nextchar = EOF;
	return (T) NULL;

      } else if (*nfiles == 1 || force_single_end_p == true) {
	if ((*input1 = Bzip2_new((*files)[0])) == NULL) {
	  fprintf(stderr,"Can't open file %s => skipping it.\n",(*files)[0]);
	  (*files) += 1;
	  (*nfiles) -= 1;
	  *nextchar = EOF;
	  return (T) NULL;
	} else {
	  *input2 = NULL;
	  (*files) += 1;
	  (*nfiles) -= 1;
	  *nextchar = '\0';
	}

      } else {
	while (*nfiles > 0 && (*input1 = Bzip2_new((*files)[0])) == NULL) {
	  fprintf(stderr,"Can't open file %s => skipping it.\n",(*files)[0]);
	  (*files)++;
	  (*nfiles)--;
	}
	if (*input1 == NULL) {
	  *nextchar = EOF;
	  return (T) NULL;
	} else {
	  (*files)++;
	  (*nfiles)--;
	  *nextchar = '\0';
	}
      }
    }

    if (*nextchar == '\0') {
      if ((*nextchar = Shortread_input_init_bzip2(*input1)) == EOF) {
	*nextchar = EOF;
	return NULL;
      }
    }

    debug(printf("** Getting header\n"));
    if ((acc = input_header_bzip2(&filterp,&restofheader,*input1,skipp)) == NULL) {
      /* fprintf(stderr,"No header\n"); */
      /* File ends after >.  Don't process, but loop again */
      *nextchar = EOF;
    } else if ((*nextchar = bzgetc(*input1)) == '\r' || *nextchar == '\n') {
      /* Process blank lines and loop again */
      while (*nextchar != EOF && ((*nextchar = bzgetc(*input1)) != '>')) {
      }
    } else if ((fulllength1 = input_oneline_bzip2(&(*nextchar),&(Read1[0]),*input1,acc,
						 /*possible_fasta_header_p*/true)) == 0) {
      /* fprintf(stderr,"length is zero\n"); */
      /* No sequence1.  Don't process, but loop again */
      /* *nextchar = EOF; */
    } else {
      /* queryseq1 is in Read1 */
      /* See what is in next line */
      if ((fulllength2 = input_oneline_bzip2(&(*nextchar),&(Read2[0]),*input1,acc,
					     /*possible_fasta_header_p*/true)) > 0) {
	/* Paired-end, single file.  queryseq1 is in Read1 and queryseq2 is in Read2 */
	queryseq1 = Shortread_new(acc,restofheader,filterp,Read1,fulllength1,
				  /*quality*/NULL,/*quality_length*/0,barcode_length,
				  invert_first_p,/*copy_acc_p*/false,skipp);
	(*queryseq2) = Shortread_new(/*acc*/NULL,/*restofheader*/NULL,filterp,Read2,fulllength2,
				     /*quality*/NULL,/*quality_length*/0,barcode_length,
				     invert_second_p,/*copy_acc_p*/false,skipp);
      } else {
	if (*input2 == NULL && *nfiles > 0 && force_single_end_p == false &&
	    (*input2 = Bzip2_new((*files)[0])) != NULL) {
	  (*files) += 1;
	  (*nfiles) -= 1;
	  nextchar2 = '\0';
	}

	if (*input2 != NULL) {
	  /* Paired-end in two files */
	  if ((acc2 = input_header_bzip2(&filterp,&restofheader2,*input2,skipp)) == NULL) {
	    /* fprintf(stderr,"No header\n"); */
	    /* File ends after >.  Don't process, but loop again */
	    (*queryseq2) = (T) NULL;
	    nextchar2 = EOF;
	  } else if ((nextchar2 = bzgetc(*input2)) == '\r' || nextchar2 == '\n') {
	    /* Process blank lines and loop again */
	    while (nextchar2 != EOF && ((nextchar2 = bzgetc(*input2)) != '>')) {
	    }
	    (*queryseq2) = (T) NULL;
	  } else if ((fulllength2 = input_oneline_bzip2(&nextchar2,&(Read2[0]),*input2,acc2,
							/*possible_fasta_header_p*/true)) == 0) {
	    /* fprintf(stderr,"length is zero\n"); */
	    /* No sequence1.  Don't process, but loop again */
	    /* *nextchar = EOF; */
	    (*queryseq2) = (T) NULL;
	  } else {
	    queryseq1 = Shortread_new(acc,restofheader,filterp,Read1,fulllength1,
				      /*quality*/NULL,/*quality_length*/0,barcode_length,
				      invert_first_p,/*copy_acc_p*/false,skipp);
	    (*queryseq2) = Shortread_new(/*acc2*/NULL,/*restofheader2*/NULL,filterp,Read2,fulllength2,
					 /*quality*/NULL,/*quality_length*/0,barcode_length,
					 invert_second_p,/*copy_acc_p*/false,skipp);
	    FREE_IN(acc2);
	    FREE_IN(restofheader2);
	  }

	} else {
	  /* Single-end: Either EOF, '>', or '+' */
	  queryseq1 = Shortread_new(acc,restofheader,filterp,Read1,fulllength1,
				    /*quality*/NULL,/*quality_length*/0,barcode_length,
				    invert_first_p,/*copy_acc_p*/false,skipp);
	  (*queryseq2) = (T) NULL;
	}
      }

      debug(printf("Returning queryseq with contents %s\n",queryseq1->contents));

      return queryseq1;
    }

  }
}
#endif


T
Shortread_read_fastq_shortreads (int *nextchar, T *queryseq2, FILE **input1, FILE **input2,
				 char ***files, int *nfiles, bool skipp,
				 int barcode_length, bool invert_first_p, bool invert_second_p) {
  T queryseq1;
  char *acc, *restofheader;
  int nextchar2;
  int fulllength, quality_length;
  bool filterp;

  while (1) {
    if (*input1 == NULL || feof(*input1)) {
      if (*input1 != NULL) {
	fclose(*input1);
	*input1 = NULL;
      }
      if (*input2 != NULL) {
	fclose(*input2);
	*input2 = NULL;
      }

      if (*nfiles == 0) {
	*nextchar = EOF;
	return (T) NULL;

      } else if (*nfiles == 1 || force_single_end_p == true) {
	if ((*input1 = FOPEN_READ_TEXT((*files)[0])) == NULL) {
	  fprintf(stderr,"Can't open file %s => skipping.\n",(*files)[0]);
	  (*files) += 1;
	  (*nfiles) -= 1;
	  *nextchar = EOF;
	  return (T) NULL;
	} else {
	  *input2 = NULL;
	  (*files) += 1;
	  (*nfiles) -= 1;
	  *nextchar = '\0';
	}

      } else {
	while (*nfiles > 0 &&
	       ((*input1 = FOPEN_READ_TEXT((*files)[0])) == NULL ||
		(*input2 = FOPEN_READ_TEXT((*files)[1])) == NULL)) {
	  fprintf(stderr,"Can't open file %s or %s => skipping both.\n",
		  (*files)[0],(*files)[1]);
	  (*files) += 2;
	  (*nfiles) -= 2;
	}
	if (*input1 == NULL) {
	  *nextchar = EOF;
	  return (T) NULL;
	} else {
	  (*files) += 2;
	  (*nfiles) -= 2;
	  *nextchar = '\0';
	}
      }
    }

    debug(printf("** Getting header\n"));
    if ((acc = input_header_fastq(&filterp,&restofheader,*input1,skipp)) == NULL) {
      /* fprintf(stderr,"No header\n"); */
      /* File ends after >.  Don't process, but loop again */
      *nextchar = EOF;
    } else {
      *nextchar = fgetc(*input1);
      if ((fulllength = input_oneline(&(*nextchar),&(Read1[0]),*input1,acc,
				      /*possible_fasta_header_p*/true)) == 0) {
	/* fprintf(stderr,"length is zero\n"); */
	/* No sequence1.  Don't process, but loop again */
	/* *nextchar = EOF; */

      } else if (*nextchar != '+') {
	/* No quality */
	queryseq1 = Shortread_new(acc,restofheader,filterp,Read1,fulllength,
				  /*quality*/NULL,/*quality_length*/0,barcode_length,
				  invert_first_p,/*copy_acc_p*/false,skipp);
      } else {
	skip_header(*input1);
	*nextchar = fgetc(*input1);
	quality_length = input_oneline(&(*nextchar),&(Quality[0]),*input1,acc,
				       /*possible_fasta_header_p*/false);
	if (quality_length != fulllength) {
	  fprintf(stderr,"Length %d of quality score differs from length %d of nucleotides in sequence %s\n",
		  quality_length,fulllength,acc);
	  abort();
	} else {
	  /* Has quality */
	  queryseq1 = Shortread_new(acc,restofheader,filterp,Read1,fulllength,
				    Quality,quality_length,barcode_length,
				    invert_first_p,/*copy_acc_p*/false,skipp);
	}
      }
    }

    if (*input2 == NULL) {
      *queryseq2 = (T) NULL;
    } else {
      if ((acc = input_header_fastq(&filterp,&restofheader,*input2,skipp)) == NULL) {
	/* fprintf(stderr,"No header\n"); */
	/* File ends after >.  Don't process, but loop again */
	nextchar2 = EOF;
      } else {
	if (allow_paired_end_mismatch_p == true) {
	  /* Do not strip endings, and keep second accession */
	  FREE_IN(restofheader);
	} else if (skipp == true) {
	  /* Do not check endings */
	} else {
	  strip_illumina_acc_ending(queryseq1->acc,acc);
	  if (strcmp(queryseq1->acc,acc)) {
	    fprintf(stderr,"Paired-end accessions %s and %s do not match\n",queryseq1->acc,acc);
	    exit(9);
	  } else {
	    FREE_IN(acc);
	    FREE_IN(restofheader);
	    acc = (char *) NULL;
	  }
	}
	nextchar2 = fgetc(*input2);
	if ((fulllength = input_oneline(&nextchar2,&(Read2[0]),*input2,acc,
					/*possible_fasta_header_p*/true)) == 0) {
	  /* fprintf(stderr,"length is zero\n"); */
	  /* No sequence2.  Don't process, but loop again */
	  /* *nextchar = EOF; */

	} else if (nextchar2 != '+') {
	  /* No quality */
	  (*queryseq2) = Shortread_new(acc,/*restofheader*/NULL,filterp,Read2,fulllength,
				       /*quality*/NULL,/*quality_length*/0,barcode_length,
				       invert_second_p,/*copy_acc_p*/false,skipp);
	} else {
	  skip_header(*input2);
	  nextchar2 = fgetc(*input2);
	  quality_length = input_oneline(&nextchar2,&(Quality[0]),*input2,acc,
					 /*possible_fasta_header_p*/false);
	  if (quality_length != fulllength) {
	    fprintf(stderr,"Length %d of quality score differs from length %d of nucleotides in sequence %s\n",
		    quality_length,fulllength,acc);
	    abort();
	  } else {
	    /* Has quality */
	    (*queryseq2) = Shortread_new(acc,/*restofheader*/NULL,filterp,Read2,fulllength,
					 Quality,quality_length,barcode_length,
					 invert_second_p,/*copy_acc_p*/false,skipp);

	  }
	}
      }
    }

    return queryseq1;
  }
}


#ifdef HAVE_ZLIB
T
Shortread_read_fastq_shortreads_gzip (int *nextchar, T *queryseq2, gzFile *input1, gzFile *input2,
				      char ***files, int *nfiles, bool skipp,
				      int barcode_length, bool invert_first_p, bool invert_second_p) {
  T queryseq1;
  char *acc, *restofheader;
  int nextchar2;
  int fulllength, quality_length;
  bool filterp;

  while (1) {
    if (*input1 == NULL || gzeof(*input1)) {
      if (*input1 != NULL) {
	gzclose(*input1);
	*input1 = NULL;
      }
      if (*input2 != NULL) {
	gzclose(*input2);
	*input2 = NULL;
      }

      if (*nfiles == 0) {
	*nextchar = EOF;
	return (T) NULL;

      } else if (*nfiles == 1 || force_single_end_p == true) {
	if ((*input1 = gzopen((*files)[0],"rb")) == NULL) {
	  fprintf(stderr,"Cannot open gzipped file %s\n",(*files)[0]);
	  exit(9);
	} else {
#ifdef HAVE_ZLIB_GZBUFFER
	  gzbuffer(*input1,GZBUFFER_SIZE);
#endif
	}
	*input2 = NULL;
	(*files) += 1;
	(*nfiles) -= 1;
	*nextchar = '\0';
	
      } else {
	if ((*input1 = gzopen((*files)[0],"rb")) == NULL) {
	  fprintf(stderr,"Cannot open gzipped file %s\n",(*files)[0]);
	  exit(9);
	} else {
#ifdef HAVE_ZLIB_GZBUFFER
	  gzbuffer(*input1,GZBUFFER_SIZE);
#endif
	}

	if ((*input2 = gzopen((*files)[1],"rb")) == NULL) {
	  fprintf(stderr,"Cannot open gzipped file %s\n",(*files)[1]);
	  exit(9);
	} else {
#ifdef HAVE_ZLIB_GZBUFFER
	  gzbuffer(*input2,GZBUFFER_SIZE);
#endif
	}

	(*files) += 2;
	(*nfiles) -= 2;
	*nextchar = '\0';
      }
    }

    debug(printf("** Getting header\n"));
    if ((acc = input_header_fastq_gzip(&filterp,&restofheader,*input1,skipp)) == NULL) {
      /* fprintf(stderr,"No header\n"); */
      /* File ends after >.  Don't process. */
      *nextchar = EOF;
    } else {
      *nextchar = gzgetc(*input1);
      if ((fulllength = input_oneline_gzip(&(*nextchar),&(Read1[0]),*input1,acc,
					   /*possible_fasta_header_p*/true)) == 0) {
	/* fprintf(stderr,"length is zero\n"); */
	/* No sequence1.  Don't process, but loop again */
	/* *nextchar = EOF; */

      } else if (*nextchar != '+') {
	/* No quality */
	queryseq1 = Shortread_new(acc,restofheader,filterp,Read1,fulllength,
				  /*quality*/NULL,/*quality_length*/0,barcode_length,
				  invert_first_p,/*copy_acc_p*/false,skipp);
      } else {
	skip_header_gzip(*input1);
	*nextchar = gzgetc(*input1);
	quality_length = input_oneline_gzip(&(*nextchar),&(Quality[0]),*input1,acc,
					    /*possible_fasta_header_p*/false);
	if (quality_length != fulllength) {
	  fprintf(stderr,"Length %d of quality score differs from length %d of nucleotides in sequence %s\n",
		  quality_length,fulllength,acc);
	  abort();
	} else {
	  /* Has quality */
	  queryseq1 = Shortread_new(acc,restofheader,filterp,Read1,fulllength,
				    Quality,quality_length,barcode_length,
				    invert_first_p,/*copy_acc_p*/false,skipp);
	}
      }
    }

    if (*input2 == NULL) {
      *queryseq2 = (T) NULL;
    } else {
      if ((acc = input_header_fastq_gzip(&filterp,&restofheader,*input2,skipp)) == NULL) {
	/* fprintf(stderr,"No header\n"); */
	/* File ends after >.  Don't process, but loop again */
	nextchar2 = EOF;
      } else {
	if (allow_paired_end_mismatch_p == true) {
	  /* Do not strip endings, and keep second accession */
	  FREE_IN(restofheader);
	} else if (skipp == true) {
	  /* Do not check endings */
	} else {
	  strip_illumina_acc_ending(queryseq1->acc,acc);
	  if (strcmp(queryseq1->acc,acc)) {
	    fprintf(stderr,"Paired-end accessions %s and %s do not match\n",queryseq1->acc,acc);
	    exit(9);
	  } else {
	    FREE_IN(acc);
	    FREE_IN(restofheader);
	    acc = (char *) NULL;
	  }
	}
	nextchar2 = gzgetc(*input2);
	if ((fulllength = input_oneline_gzip(&nextchar2,&(Read2[0]),*input2,acc,
					     /*possible_fasta_header_p*/true)) == 0) {
	  /* fprintf(stderr,"length is zero\n"); */
	  /* No sequence2.  Don't process, but loop again */
	  /* *nextchar = EOF; */

	} else if (nextchar2 != '+') {
	  /* No quality */
	  (*queryseq2) = Shortread_new(acc,/*restofheader*/NULL,filterp,Read2,fulllength,
				       /*quality*/NULL,/*quality_length*/0,barcode_length,
				       invert_second_p,/*copy_acc_p*/false,skipp);
	} else {
	  skip_header_gzip(*input2);
	  nextchar2 = gzgetc(*input2);
	  quality_length = input_oneline_gzip(&nextchar2,&(Quality[0]),*input2,acc,
					      /*possible_fasta_header_p*/false);
	  if (quality_length != fulllength) {
	    fprintf(stderr,"Length %d of quality score differs from length %d of nucleotides in sequence %s\n",
		    quality_length,fulllength,acc);
	    abort();
	  } else {
	    /* Has quality */
	    (*queryseq2) = Shortread_new(acc,/*restofheader*/NULL,filterp,Read2,fulllength,
					 Quality,quality_length,barcode_length,
					 invert_second_p,/*copy_acc_p*/false,skipp);
	  }
	}
      }
    }

    return queryseq1;
  }
}
#endif


#ifdef HAVE_BZLIB
T
Shortread_read_fastq_shortreads_bzip2 (int *nextchar, T *queryseq2, Bzip2_T *input1, Bzip2_T *input2,
				       char ***files, int *nfiles, bool skipp,
				       int barcode_length, bool invert_first_p, bool invert_second_p) {
  T queryseq1;
  char *acc, *restofheader;
  int nextchar2;
  int fulllength, quality_length;
  bool filterp;

  while (1) {
    if (*input1 == NULL || bzeof(*input1)) {
      if (*input1 != NULL) {
	Bzip2_free(&(*input1));
	*input1 = NULL;
      }
      if (*input2 != NULL) {
	Bzip2_free(&(*input2));
	*input2 = NULL;
      }

      if (*nfiles == 0) {
	*nextchar = EOF;
	return (T) NULL;

      } else if (*nfiles == 1 || force_single_end_p == true) {
	if ((*input1 = Bzip2_new((*files)[0])) == NULL) {
	  fprintf(stderr,"Cannot open bzip2 file %s\n",(*files)[0]);
	  exit(9);
	}
	*input2 = NULL;
	(*files) += 1;
	(*nfiles) -= 1;
	*nextchar = '\0';
	
      } else {
	if ((*input1 = Bzip2_new((*files)[0])) == NULL) {
	  fprintf(stderr,"Cannot open bzip2 file %s\n",(*files)[0]);
	  exit(9);
	}

	if ((*input2 = Bzip2_new((*files)[1])) == NULL) {
	  fprintf(stderr,"Cannot open bzip2 file %s\n",(*files)[1]);
	  exit(9);
	}

	(*files) += 2;
	(*nfiles) -= 2;
	*nextchar = '\0';
      }
    }

    debug(printf("** Getting header\n"));
    if ((acc = input_header_fastq_bzip2(&filterp,&restofheader,*input1,skipp)) == NULL) {
      /* fprintf(stderr,"No header\n"); */
      /* File ends after >.  Don't process. */
      *nextchar = EOF;
    } else {
      *nextchar = bzgetc(*input1);
      if ((fulllength = input_oneline_bzip2(&(*nextchar),&(Read1[0]),*input1,acc,
					    /*possible_fasta_header_p*/true)) == 0) {
	/* fprintf(stderr,"length is zero\n"); */
	/* No sequence1.  Don't process, but loop again */
	/* *nextchar = EOF; */

      } else if (*nextchar != '+') {
	/* No quality */
	queryseq1 = Shortread_new(acc,restofheader,filterp,Read1,fulllength,
				  /*quality*/NULL,/*quality_length*/0,barcode_length,
				  invert_first_p,/*copy_acc_p*/false,skipp);
      } else {
	skip_header_bzip2(*input1);
	*nextchar = bzgetc(*input1);
	quality_length = input_oneline_bzip2(&(*nextchar),&(Quality[0]),*input1,acc,
					     /*possible_fasta_header_p*/false);
	if (quality_length != fulllength) {
	  fprintf(stderr,"Length %d of quality score differs from length %d of nucleotides in sequence %s\n",
		  quality_length,fulllength,acc);
	  abort();
	} else {
	  /* Has quality */
	  queryseq1 = Shortread_new(acc,restofheader,filterp,Read1,fulllength,
				    Quality,quality_length,barcode_length,
				    invert_first_p,/*copy_acc_p*/false,skipp);
	}
      }
    }

    if (*input2 == NULL) {
      *queryseq2 = (T) NULL;
    } else {
      if ((acc = input_header_fastq_bzip2(&filterp,&restofheader,*input2,skipp)) == NULL) {
	/* fprintf(stderr,"No header\n"); */
	/* File ends after >.  Don't process, but loop again */
	nextchar2 = EOF;
      } else {
	if (allow_paired_end_mismatch_p == true) {
	  /* Do not strip endings, and keep second accession */
	  FREE_IN(restofheader);
	} else if (skipp == true) {
	  /* Do not check endings */
	} else {
	  strip_illumina_acc_ending(queryseq1->acc,acc);
	  if (strcmp(queryseq1->acc,acc)) {
	    fprintf(stderr,"Paired-end accessions %s and %s do not match\n",queryseq1->acc,acc);
	    exit(9);
	  } else {
	    FREE_IN(acc);
	    FREE_IN(restofheader);
	    acc = (char *) NULL;
	  }
	}
	nextchar2 = bzgetc(*input2);
	if ((fulllength = input_oneline_bzip2(&nextchar2,&(Read2[0]),*input2,acc,
					      /*possible_fasta_header_p*/true)) == 0) {
	  /* fprintf(stderr,"length is zero\n"); */
	  /* No sequence2.  Don't process, but loop again */
	  /* *nextchar = EOF; */

	} else if (nextchar2 != '+') {
	  /* No quality */
	  (*queryseq2) = Shortread_new(acc,/*restofheader*/NULL,filterp,Read2,fulllength,
				       /*quality*/NULL,/*quality_length*/0,barcode_length,
				       invert_second_p,/*copy_acc_p*/false,skipp);
	} else {
	  skip_header_bzip2(*input2);
	  nextchar2 = bzgetc(*input2);
	  quality_length = input_oneline_bzip2(&nextchar2,&(Quality[0]),*input2,acc,
					       /*possible_fasta_header_p*/false);
	  if (quality_length != fulllength) {
	    fprintf(stderr,"Length %d of quality score differs from length %d of nucleotides in sequence %s\n",
		    quality_length,fulllength,acc);
	    abort();
	  } else {
	    /* Has quality */
	    (*queryseq2) = Shortread_new(acc,/*restofheader*/NULL,filterp,Read2,fulllength,
					 Quality,quality_length,barcode_length,
					 invert_second_p,/*copy_acc_p*/false,skipp);
	  }
	}
      }
    }

    return queryseq1;
  }
}
#endif


/* Calling procedure needs to print the initial ">", if desired */
void
Shortread_print_header (FILE *fp, T queryseq1, T queryseq2) {

  if (queryseq2 == NULL || queryseq2->acc == NULL) {
    fprintf(fp,"%s",queryseq1->acc);
  } else {
    fprintf(fp,"%s,%s",queryseq1->acc,queryseq2->acc);
  }

  if (queryseq1->restofheader == NULL || queryseq1->restofheader[0] == '\0') {
    /* Don't print restofheader */
  } else {
    fprintf(fp," %s",queryseq1->restofheader);
  }

  fprintf(fp,"\n");

  return;
}


void
Shortread_print_query_singleend_fasta (FILE *fp, T queryseq) {
  fprintf(fp,">");
  Shortread_print_header(fp,queryseq,/*queryseq2*/NULL);
  /* fprintf(fp,"\n"); -- included in header */
  Shortread_print_oneline(fp,queryseq);
  fprintf(fp,"\n");

  return;
}

void
Shortread_print_query_singleend_fastq (FILE *fp, T queryseq) {
  fprintf(fp,"@");
  Shortread_print_header(fp,queryseq,/*queryseq2*/NULL);
  /* fprintf(fp,"\n"); -- included in header */
  Shortread_print_oneline(fp,queryseq);
  fprintf(fp,"\n");

  if (queryseq->quality != NULL) {
    fprintf(fp,"+\n");
    Shortread_print_quality(fp,queryseq,/*hardclip_low*/0,/*hardclip_high*/0,
			    /*shift*/0,/*choppedp*/false);
    fprintf(fp,"\n");
  }

  return;
}

void
Shortread_print_query_pairedend_fasta (FILE *fp, T queryseq1, T queryseq2,
				       bool invert_first_p, bool invert_second_p) {
  fprintf(fp,">");
  Shortread_print_header(fp,queryseq1,queryseq2);
  /* fprintf(fp,"\n"); -- included in header */

  if (invert_first_p == true) {
    Shortread_print_oneline_revcomp(fp,queryseq1);
    fprintf(fp,"\n");
  } else {
    Shortread_print_oneline(fp,queryseq1);
    fprintf(fp,"\n");
  }

  if (invert_second_p == true) {
    Shortread_print_oneline_revcomp(fp,queryseq2);
    fprintf(fp,"\n");
  } else {
    Shortread_print_oneline(fp,queryseq2);
    fprintf(fp,"\n");
  }

  return;
}


void
Shortread_print_query_pairedend_fastq (FILE *fp1, FILE *fp2, T queryseq1, T queryseq2,
				      bool invert_first_p, bool invert_second_p) {
  /* First end */
  if (queryseq2->acc == NULL) {
    fprintf(fp1,"@%s/1\n",queryseq1->acc);
  } else {
    fprintf(fp2,"@%s\n",queryseq1->acc); /* Allowing paired-end name mismatch */
  }

  if (invert_first_p == true) {
    Shortread_print_oneline_revcomp(fp1,queryseq1);
    fprintf(fp1,"\n");
    if (queryseq1->quality != NULL) {
      fprintf(fp1,"+\n");
      Shortread_print_quality_revcomp(fp1,queryseq1,/*hardclip_low*/0,/*hardclip_high*/0,
				      /*shift*/0,/*choppedp*/false);
      fprintf(fp1,"\n");
    }
  } else {
    Shortread_print_oneline(fp1,queryseq1);
    fprintf(fp1,"\n");
    if (queryseq1->quality != NULL) {
      fprintf(fp1,"+\n");
      Shortread_print_quality(fp1,queryseq1,/*hardclip_low*/0,/*hardclip_high*/0,
			      /*shift*/0,/*choppedp*/false);
      fprintf(fp1,"\n");
    }
  }

  /* Second end */
  if (queryseq2->acc == NULL) {
    fprintf(fp2,"@%s/2\n",queryseq1->acc); /* Acc stored only for first end, not second end */
  } else {
    fprintf(fp2,"@%s\n",queryseq2->acc); /* Allowing paired-end name mismatch */
  }

  if (invert_second_p == true) {
    Shortread_print_oneline_revcomp(fp2,queryseq2);
    fprintf(fp2,"\n");
    if (queryseq2->quality != NULL) {
      fprintf(fp2,"+\n");
      Shortread_print_quality_revcomp(fp2,queryseq2,/*hardclip_low*/0,/*hardclip_high*/0,
				      /*shift*/0,/*chopped*/false);
      fprintf(fp2,"\n");
    }
  } else {
    Shortread_print_oneline(fp2,queryseq2);
    fprintf(fp2,"\n");
    if (queryseq2->quality != NULL) {
      fprintf(fp2,"+\n");
      Shortread_print_quality(fp2,queryseq2,/*hardclip_low*/0,/*hardclip_high*/0,
			      /*shift*/0,/*choppedp*/false);
      fprintf(fp2,"\n");
    }
  }

  return;
}


void
Shortread_print_oneline (FILE *fp, T this) {
  int i = 0;

  if (this->fulllength == 0 || isspace(this->contents[0])) {
    fprintf(fp,"(null)");
  } else {
    for (i = 0; i < this->fulllength; i++) {
      fprintf(fp,"%c",this->contents[i]);
    }
    for (i = 0; i < this->choplength; i++) {
      fprintf(fp,"%c",this->chop[i]);
    }
  }
  return;
}

void
Shortread_print_oneline_revcomp (FILE *fp, T this) {
  int i = 0;

  for (i = this->fulllength-1; i >= 0; --i) {
    fprintf(fp,"%c",complCode[(int) this->contents[i]]);
  }
  for (i = this->choplength-1; i >= 0; --i) {
    fprintf(fp,"%c",complCode[(int) this->chop[i]]);
  }

  return;
}


void
Shortread_print_chopped (FILE *fp, T this, int hardclip_low, int hardclip_high) {
  int i;

  if (this->fulllength == 0 || isspace(this->contents[0])) {
    fprintf(fp,"(null)");
  } else {
    for (i = hardclip_low; i < this->fulllength - hardclip_high; i++) {
      fprintf(fp,"%c",this->contents[i]);
    }
  }
  return;
}

void
Shortread_print_chopped_revcomp (FILE *fp, T this, int hardclip_low, int hardclip_high) {
  int i;

  for (i = this->fulllength - 1 - hardclip_low; i >= hardclip_high; --i) {
    fprintf(fp,"%c",complCode[(int) this->contents[i]]);
  }

  return;
}


void
Shortread_print_barcode (FILE *fp, T this) {

  if (this->barcode != NULL) {
    fprintf(fp,"\tXB:Z:%s",this->barcode);
  }
    
  return;
}

void
Shortread_print_chop (FILE *fp, T this, bool invertp) {
  int i;

  if (this->chop != NULL) {
    fprintf(fp,"\tXP:Z:");
    if (invertp == false) {
      fprintf(fp,"%s",this->chop);
    } else {
      for (i = this->choplength - 1; i >= 0; i--) {
	fprintf(fp,"%c",complCode[(int) this->chop[i]]);
      }
    }
  }
    
  return;
}


void
Shortread_print_chop_symbols (FILE *fp, T this) {
  int i;

  for (i = 0; i < this->choplength; i++) {
    fprintf(fp,"*");
  }
  return;
}

void
Shortread_print_quality (FILE *fp, T this, int hardclip_low, int hardclip_high,
			int shift, bool show_chopped_p) {
  int i;
  int c;

  if (this->quality == NULL) {
    fprintf(fp,"*");
  } else {
    for (i = hardclip_low; i < this->fulllength - hardclip_high; i++) {
      if ((c = this->quality[i] + shift) <= 32) {
	fprintf(stderr,"Warning: With a quality-print-shift of %d, QC score %c becomes non-printable.  May need to specify --quality-protocol or --quality-print-shift\n",
		shift,this->quality[i]);
	abort();
      } else {
	fprintf(fp,"%c",c);
      }
    }

    if (show_chopped_p == true) {
      assert(hardclip_high == 0);
      for (i = 0; i < this->choplength; i++) {
	if ((c = this->chop_quality[i] + shift) <= 32) {
	  fprintf(stderr,"Warning: With a quality-print-shift of %d, QC score %c becomes non-printable.  May need to specify --quality-protocol or --quality-print-shift\n",
		  shift,this->chop_quality[i]);
	  abort();
	} else {
	  fprintf(fp,"%c",c);
	}
      }
    }

  }
  return;
}

void
Shortread_print_quality_revcomp (FILE *fp, T this, int hardclip_low, int hardclip_high,
				int shift, bool show_chopped_p) {
  int i;
  int c;

  if (this->quality == NULL) {
    fprintf(fp,"*");
  } else {
    for (i = this->fulllength - 1 - hardclip_low; i >= hardclip_high; --i) {
      if ((c = this->quality[i] + shift) <= 32) {
	fprintf(stderr,"Warning: With a quality-print-shift of %d, QC score %c becomes non-printable.  May need to specify --quality-protocol or --quality-print-shift\n",
		shift,this->quality[i]);
	abort();
      } else {
	fprintf(fp,"%c",c);
      }
    }

    if (show_chopped_p == true) {
      /* assert(hardclip_low == 0); */
      for (i = this->choplength - 1; i >= 0; i--) {
	if ((c = this->chop_quality[i] + shift) <= 32) {
	  fprintf(stderr,"Warning: With a quality-print-shift of %d, QC score %c becomes non-printable.  May need to specify --quality-protocol or --quality-print-shift\n",
		  shift,this->chop_quality[i]);
	  abort();
	} else {
	  fprintf(fp,"%c",c);
	}
      }
    }
  }

  return;
}

void
Shortread_print_oneline_uc (FILE *fp, T this) {
  int i = 0;

  for (i = 0; i < this->fulllength; i++) {
    fprintf(fp,"%c",this->contents_uc[i]);
  }
  return;
}

void
Shortread_print_oneline_revcomp_uc (FILE *fp, T this) {
  int i = 0;

  for (i = this->fulllength-1; i >= 0; --i) {
    fprintf(fp,"%c",complCode[(int) this->contents_uc[i]]);
  }
  return;
}
