static char rcsid[] = "$Id: sequence.c,v 1.100 2010-07-22 00:11:06 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifndef HAVE_MEMCPY
# define memcpy(d,s,n) bcopy((s),(d),(n))
#endif
#ifndef HAVE_MEMMOVE
# define memmove(d,s,n) bcopy((s),(d),(n))
#endif

#include "sequence.h"

#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>		/* For rindex */
#include <ctype.h>		/* For iscntrl and isspace */
#include "mem.h"
#include "complement.h"
#include "intlist.h"
#include "fopen.h"
#include "md5.h"

/* Before setting DEBUG, may want to reduce MAXSEQLEN in sequence.h */
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

/* Dynamic programming */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif




/***********************************************************************
 *    Definitions:
 *
 *   TTTTTT ACGT ...... ACGT AAAAAA
 *          <- trimlength ->
 *   <-------- fulllength -------->
 *          ^trimstart
 *   ^contents
 *
 *   Trimming is determined by Oligoindex_set_inquery(), based on
 *   finding unique 8-mers on each end.
 ************************************************************************/


#define T Sequence_T
struct T {
  char *acc;			/* Accession */
  char *restofheader;		/* Rest of header */
  char *contents;		/* Original sequence, ends with '\0' */
  char *contents_alloc;		/* Allocation */
  int fulllength;		/* Full length (not including chopped sequence) */

#ifdef GSNAP
  char *contents_uc;		/* Original sequence, ends with '\0' */
  char *contents_uc_alloc;	/* Allocation */
  char *chop;
  int choplength;
  char *quality;		/* For Illumina short reads read via FASTQ */
  char *quality_alloc;		/* Allocation */
#endif

  int trimstart;		/* Start of trim */
  int trimend;			/* End of trim */
#ifdef PMAP
  int fulllength_given;		/* Full length minus implicit stop codon at end */
#endif
  int subseq_offset;		/* Used only for subsequences */
  int skiplength;		/* Used only for sequences longer than MAXSEQLEN */
  /* bool free_contents_p; */
};

char *
Sequence_accession (T this) {
  return this->acc;
}

char *
Sequence_fullpointer (T this) {
  return this->contents;
}

#ifdef GSNAP
char *
Sequence_fullpointer_uc (T this) {
  return this->contents_uc;
}
#endif

char *
Sequence_trimpointer (T this) {
  return &(this->contents[this->trimstart]);
}

#ifdef GSNAP
int
Sequence_choplength (T this) {
  return this->choplength;
}

char *
Sequence_quality (T this) {
  return this->quality;
}
#endif


int
Sequence_ntlength (T this) {
#ifdef PMAP
  return 3*this->fulllength;
#else
  return this->fulllength;
#endif
}

int
Sequence_fulllength (T this) {
  return this->fulllength;
}

int
Sequence_fulllength_given (T this) {
#ifdef PMAP
  return this->fulllength_given;
#else
  return this->fulllength;
#endif
}

int
Sequence_trimlength (T this) {
  return this->trimend - this->trimstart + 1;
}

void
Sequence_trim (T this, int trim_start, int trim_end) {
  this->trimstart = trim_start;
  this->trimend = trim_end;
  return;
}

int
Sequence_trim_start (T this) {
  return this->trimstart;
}

int
Sequence_trim_end (T this) {
  return this->trimend;
}

int
Sequence_subseq_offset (T this) {
  return this->subseq_offset;
}

int
Sequence_skiplength (T this) {
  return this->skiplength;
}

void
Sequence_free (T *old) {
  if (*old) {
    if ((*old)->restofheader != NULL) {
      FREE((*old)->restofheader);
    }
    if ((*old)->acc != NULL) {
      FREE((*old)->acc);
    }

    FREE((*old)->contents_alloc);

#ifdef GSNAP
    FREE((*old)->chop);
    FREE((*old)->quality_alloc);
    FREE((*old)->contents_uc_alloc);
#endif

    FREE(*old);
  }
  return;
}

#ifdef PMAP
T
Sequence_convert_to_nucleotides (T this) {
  T new = (T) MALLOC(sizeof(*new));
  int i;

  new->acc = (char *) NULL;
  new->restofheader = (char *) NULL;
  new->fulllength = this->fulllength*3;
  new->fulllength_given = this->fulllength_given*3;
  new->contents = new->contents_alloc = (char *) CALLOC(new->fulllength+1,sizeof(char));
  for (i = 0; i < new->fulllength; i++) {
    new->contents[i] = '?';
  }

#ifdef GSNAP
  new->contents_uc = new->contents_uc_alloc = (char *) NULL;
  new->chop = (char *) NULL;
  new->choplength = 0;
  new->quality = new->quality_alloc = (char *) NULL;
#endif

  new->trimstart = 0;
  new->trimend = new->fulllength_given;
  /* new->free_contents_p = true; */
  new->subseq_offset = 0;
  new->skiplength = this->skiplength;

  return new;
}
#endif

#if 0
int
Sequence_count_bad (T this, int pos, int max, int direction) {
  int nbad = 0;

  if (direction > 0) {
    while (--max >= 0 && pos < this->fulllength) {
      if (this->contents[pos] == 'X') {
	nbad++;
      }
      pos++;
    }
  } else {
    while (--max >= 0 && pos >= 0) {
      if (this->contents[pos] == 'X') {
	nbad++;
      }
      pos--;
    }
  }

  return nbad;
}
#endif


#define HEADERLEN 512
#define DISCARDLEN 8192

static char Header[HEADERLEN];
static char Discard[DISCARDLEN];

static char Sequence[1+MAXSEQLEN+1]; /* Used by Sequence_read_unlimited */
static char Sequence1[HALFLEN+1]; /* 1 at end for '\0' */
static char Sequence2[HALFLEN+3]; /* 1 at end for '\0' and 2 extra in cyclic part for '\n' and '\0' */

static char *Firsthalf;
static char *Secondhalf;
static int Initc = '\0';


/* The first element of Sequence is always the null character, to mark
   the end of the string */

/* Skipping of dashes might still be buggy */
/*
#define DASH '-'
*/


/* Returns '>' if FASTA file, first sequence char if not */
int
Sequence_input_init (FILE *fp) {
  int c;
  bool okayp = false;

  Header[0] = '\0';
  Sequence[0] = '\0';
  Firsthalf = &(Sequence1[0]);
  Secondhalf = &(Sequence2[0]);

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

static void
blank_header (T this) {
  this->acc = (char *) CALLOC(strlen("NO_HEADER")+1,sizeof(char));
  strcpy(this->acc,"NO_HEADER");
  this->restofheader = (char *) CALLOC(1,sizeof(char));
  this->restofheader[0] = '\0';
  return;
}

static char *
input_header (FILE *fp, T this) {
  char *p;
  size_t length;

  if (feof(fp)) {
    return NULL;
  } else if (fgets(&(Header[0]),HEADERLEN,fp) == NULL) {
    /* File must terminate after > */
    return NULL;
  }

  if ((p = rindex(&(Header[0]),'\n')) != NULL) {
    *p = '\0';
  } else {
    /* Eliminate rest of header from input */
    while (fgets(&(Discard[0]),DISCARDLEN,fp) != NULL &&
	   rindex(&(Discard[0]),'\n') == NULL) ;
  }

  p = &(Header[0]);
  while (*p != '\0' && !isspace((int) *p)) {
    p++;
  }
  if (*p == '\0') {
    /* Accession only */
    length = (p - &(Header[0]))/sizeof(char);
    this->acc = (char *) CALLOC(length+1,sizeof(char));
    strcpy(this->acc,Header);
    this->restofheader = (char *) CALLOC(1,sizeof(char));
    this->restofheader[0] = '\0';
  } else {
    *p = '\0';
    length = (p - &(Header[0]))/sizeof(char);
    this->acc = (char *) CALLOC(length+1,sizeof(char));
    strcpy(this->acc,Header);
    p++;
    this->restofheader = (char *) CALLOC(strlen(p)+1,sizeof(char));
    strcpy(this->restofheader,p);
  }

  return this->acc;
} 

#ifdef GSNAP
/* Deletes "/1" or "/2" endings */
static void
strip_illumina_acc_ending (char *acc) {
  char *p;

  p = acc;
  while (*p != '\0') {
    p++;
  }

  /* Delete "/1" or "/2" endings */
  if (p[-2] == '/' && (p[-1] == '1' || p[-1] == '2')) {
    p -= 2;
    *p = '\0';
  }

  return;
}

static char *
input_header_fastq (FILE *fp, T this) {
  char *p;
  size_t length;

  if (feof(fp)) {
    return NULL;
  } else if (fgets(&(Header[0]),HEADERLEN,fp) == NULL) {
    /* File must terminate after > */
    return NULL;
  }

  if ((p = rindex(&(Header[0]),'\n')) != NULL) {
    *p = '\0';
  } else {
    /* Eliminate rest of header from input */
    while (fgets(&(Discard[0]),DISCARDLEN,fp) != NULL &&
	   rindex(&(Discard[0]),'\n') == NULL) ;
  }

  p = &(Header[0]);
  while (*p != '\0' && !isspace((int) *p)) {
    p++;
  }

  *p = '\0';

  length = (p - &(Header[0]))/sizeof(char);
  this->acc = (char *) CALLOC(length+1,sizeof(char));
  strcpy(this->acc,Header);
  this->restofheader = (char *) CALLOC(1,sizeof(char));
  this->restofheader[0] = '\0';

  return this->acc;
} 
#endif

#ifdef GSNAP
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
#endif

#ifdef DEBUG
static void
print_contents (char *p, int length) {
  int i;
  FILE *fp = stdout;
	
  fprintf(fp,"\"");
  for (i = 0; i < length; i++) {
    if (*p == '\0') {
      fprintf(fp,"_");
    } else {
      fprintf(fp,"%c",*p);
    }
    p++;
  }
  fprintf(fp,"\"\n");
  return;
}
#endif


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


static int
read_first_half (int *nextchar, bool *eolnp, FILE *fp) {
  int remainder, strlenp;
  char *ptr, *p = NULL;
  int c;

  ptr = &(Firsthalf[0]);
  if (*nextchar != '>') {
    *ptr++ = (char) *nextchar;
  }
  remainder = (&(Firsthalf[HALFLEN]) - ptr)/sizeof(char);

  while (1) {
    if (remainder <= 0) {
      debug(printf("remainder <= 0.  Returning false\n"));
      *nextchar = EOF;
      debug1(printf("read_first_half returning length1 of %d\n",(ptr - &(Firsthalf[0]))/sizeof(char)));
      return (ptr - &(Firsthalf[0]))/sizeof(char);

    } else if (feof(fp)) {
      /* EOF in middle of line */
      debug(printf("EOF.  Returning true\n"));
      *nextchar = EOF;
      debug1(printf("read_first_half returning length1 of %d\n",(ptr - &(Firsthalf[0]))/sizeof(char)));
      return (ptr - &(Firsthalf[0]))/sizeof(char);

    } else if (*eolnp == true) {
      /* Peek at character after eoln */
      if ((c = fgetc(fp)) == EOF || c == '>') {
	debug(printf("c == EOF or >.  Returning true\n"));
	*nextchar = c;
	return (ptr - &(Firsthalf[0]))/sizeof(char);
      } else if (iscntrl(c)) {
	debug(printf("c == control char.  Continuing\n"));
#ifdef DASH
      } else if (c == DASH) {
	debug(printf("c == dash.  Continuing\n"));
#endif
      } else if (isspace(c)) {
	*eolnp = true;
	debug(printf("c == NULL.  Continuing\n"));
      } else {
	*ptr++ = (char) c;
	remainder--;
	*eolnp = false;
	p = NULL;
	debug(printf("c == sth.  Continuing\n"));
      }

    } else {
      debug(printf("Trying to read remainder of %d\n",remainder));
      if (p != NULL) {
	strlenp = strlen(p);
	memmove(ptr,p,strlenp);
        ptr[strlenp] = '\0';
	debug(printf("Did memmove of %d chars at %p to %p\n",strlenp,p,ptr));
      } else {
	p = fgets(ptr,remainder+1,fp);
      }
      if (p == NULL) {
	debug(printf("line == NULL.  Returning true\n"));
	*nextchar = EOF;
	debug1(printf("read_first_half returning length1 of %d\n",(ptr - &(Firsthalf[0]))/sizeof(char)));
	return (ptr - &(Firsthalf[0]))/sizeof(char);
      } else {
	debug(printf("Read %s.\n",ptr));
	while ((p = find_bad_char(ptr)) != NULL) {
	  /* Handle PC line feed ^M */
	  ptr = p++;
	  strlenp = strlen(p);
	  memmove(ptr,p,strlenp);
	  ptr[strlenp] = '\0';
	  debug(printf("Found control-M/space.  Did memmove of %d chars at %p to %p\n",strlenp,p,ptr));
	}
	if ((p = index(ptr,'\n')) != NULL) {
	  ptr = p;
	  *eolnp = true;
	  debug(printf("line == EOLN.  Continuing\n"));
	} else {
	  ptr += strlen(ptr);
	  *eolnp = false;
	  p = NULL;
	  debug(printf("line != EOLN.  Continuing\n"));
	}
	remainder = (&(Firsthalf[HALFLEN]) - ptr)/sizeof(char);
      }
    }

    debug(print_contents(&(Firsthalf[0]),HALFLEN+1));
  }
}

/* returns skip length */
static int
read_second_half (int *nextchar, char **pointer2a, int *length2a, char **pointer2b, int *length2b,
		  bool eolnp, FILE *fp) {
  int skiplength, ncycles = 0, remainder, terminator, strlenp;
  char *ptr;
  char *p = NULL;
  int c;
  
  ptr = &(Secondhalf[0]);
  remainder = (&(Secondhalf[HALFLEN+2]) - ptr)/sizeof(char);

  while (1) {
    debug(printf("\nEnd: %d\n",remainder));

    if (feof(fp)) {
      debug(printf("EOF.  Returning\n"));
      *nextchar = EOF;
      break;

    } else if (remainder <= 0) {
      ptr = &(Secondhalf[0]);
      remainder = (&(Secondhalf[HALFLEN+2]) - ptr)/sizeof(char);
      ncycles++;
      debug(printf("remainder <= 0.  Cycling\n"));

    } else if (eolnp == true) {
      /* Peek at character after eoln */
      if ((c = fgetc(fp)) == EOF || c == '>') {
	debug(printf("c == EOF or >.  Returning\n"));
	*nextchar = c;
	break;
      } else if (iscntrl(c)) {
	debug(printf("c == control char.  Continuing\n"));
#ifdef DASH
      } else if (c == DASH) {
	debug(printf("c == dash.  Continuing\n"));
#endif
      } else if (isspace(c)) {
	debug(printf("c == NULL.  Continuing\n"));
      } else {
	*ptr++ = (char) c;
	remainder--;
	eolnp = false;
	p = NULL;
	debug(printf("c == sth.  Continuing\n"));
      }
      
    } else {
      if (p != NULL) {
        strlenp = strlen(p);
	memmove(ptr,p,strlenp);
        ptr[strlenp] = '\0';
	debug(printf("Did memmove of %d chars at %p to %p\n",strlenp,p,ptr));
      } else {
	p = fgets(ptr,remainder+1,fp);
      }
      if (p == NULL) {
	debug(printf("line == NULL.  Returning\n"));
	*nextchar = EOF;
	break;
      } else {
	debug(printf("Read %s.\n",ptr));
	while ((p = find_bad_char(ptr)) != NULL) {
	  /* Handle PC line feed ^M */
	  ptr = p++;
	  strlenp = strlen(p);
	  memmove(ptr,p,strlenp);
	  ptr[strlenp] = '\0';
	  debug(printf("Found control-M/space.  Did memmove of %d chars at %p to %p\n",strlenp,p,ptr));
	} 
	if ((p = index(ptr,'\n')) != NULL) {
	  ptr = p;
	  eolnp = true;
	  debug(printf("line == EOLN.  Continuing\n"));
	} else {
	  ptr += strlen(ptr);
	  eolnp = false;
	  p = NULL;
	  debug(printf("line != EOLN.  Continuing\n"));
	}
	remainder = (&(Secondhalf[HALFLEN+2]) - ptr)/sizeof(char);
      }
    }

    debug(print_contents(&(Secondhalf[0]),HALFLEN+3));
  }

  terminator = (ptr - &(Secondhalf[0]))/sizeof(char);
  debug(printf("ncycles = %d, terminator is %d\n",ncycles,terminator));
  if (ncycles == 0) {
    *length2a = 0;
    if (terminator < HALFLEN) {
      skiplength = 0;
    } else {
      skiplength = terminator-HALFLEN;
    }
  } else {
    *length2a = HALFLEN-terminator;
    skiplength = ncycles*(HALFLEN+2) + terminator-HALFLEN;
  }
  if (*length2a <= 0) {
    *length2a = 0;
    *pointer2a = (char *) NULL;
  } else {
    *pointer2a = &(Secondhalf[HALFLEN+2-(*length2a)]);
  }
  if (terminator == 0) {
    *length2b = 0;
    *pointer2b = (char *) NULL;
  } else if (terminator > HALFLEN) {
    *length2b = HALFLEN;
    *pointer2b = &(ptr[-(*length2b)]);
  } else {
    *length2b = terminator;
    *pointer2b = &(Secondhalf[0]);
  }

  return skiplength;
}


#ifdef GSNAP
static int
input_oneline (int *nextchar, FILE *fp, bool possible_fasta_header_p) {
  int remainder, strlenp;
  char *ptr, *p = NULL;

  debug(printf("Entering input_oneline with nextchar = %c\n",*nextchar));

  ptr = &(Firsthalf[0]);
  remainder = (&(Firsthalf[HALFLEN]) - ptr)/sizeof(char);
  if (*nextchar == EOF || (possible_fasta_header_p == true && *nextchar == '>')) {
    debug(printf("Returning 0\n"));
    return 0;
  } else {
    *ptr++ = (char) *nextchar;
    if ((p = fgets(ptr,remainder+1,fp)) == NULL) {
      /* NULL if file ends with a blank line */
      printf("Blank line. read %s.\n",ptr);
    } else {
      debug(printf("Read %s.\n",ptr));
      while ((p = find_bad_char(ptr)) != NULL) {
	/* Handle PC line feed ^M */
	ptr = p++;
	strlenp = strlen(p);
	memmove(ptr,p,strlenp);
	ptr[strlenp] = '\0';
	debug(printf("Found control-M/space.  Did memmove of %d chars at %p to %p\n",strlenp,p,ptr));
      }

      if ((p = index(ptr,'\n')) != NULL) {
	*p = '\0';
	debug(printf("Now string is %s.\n",ptr));
      } else if (feof(fp)) {
	/* No line feed, but end of file.  Handle below. */
	debug(printf("End of file seen\n"));
      } else {
	/* No line feed, but not end of file */
	fprintf(stderr,"Line %s is too long for allocated buffer size of %d.  Aborting.\n",&(Firsthalf[0]),HALFLEN);
	exit(9);
      }
    }
    ptr += strlen(ptr);

    /* Peek at character after eoln */
    if (feof(fp)) {
      *nextchar = EOF;
    } else {
      while ((*nextchar = fgetc(fp)) != EOF && (*nextchar == '\n' || isspace(*nextchar))) {
      }
    }

    debug(printf("Returning %ld\n",(ptr - &(Firsthalf[0]))/sizeof(char)));
    return (ptr - &(Firsthalf[0]))/sizeof(char);
  }
}
#endif



/* Returns sequence length */
static int
input_sequence (int *nextchar, char **pointer1, int *length1, char **pointer2a, int *length2a,
		char **pointer2b, int *length2b, int *skiplength, FILE *fp) {
  bool eolnp = true;

  *pointer1 = &(Firsthalf[0]);
  *pointer2a = (char *) NULL;
  *length2a = 0;
  *pointer2b = (char *) NULL;
  *length2b = 0;

  if ((*length1 = read_first_half(&(*nextchar),&eolnp,fp)) == 0) {
    *pointer1 = (char *) NULL;
    *skiplength = 0;
  } else if (*length1 < HALFLEN) {
    *skiplength = 0;
  } else {
    *skiplength = read_second_half(&(*nextchar),&(*pointer2a),&(*length2a),
				   &(*pointer2b),&(*length2b),eolnp,fp);
    debug1(printf("read_second_half returns skiplength of %d, length2a=%d, length2b=%d\n",
		  *skiplength,*length2a,*length2b));
  }

  debug1(printf("length1 = %d, length2a = %d, length2b = %d\n",
		*length1,*length2a,*length2b));

  return (*length1) + (*length2a) + (*length2b);
}

/* Used only by extern procedures (outside of this file).  Internal
   procedures have their own specialized creators. */
T
Sequence_genomic_new (char *contents, int length) {
  T new = (T) MALLOC(sizeof(*new));

  new->acc = (char *) NULL;
  new->restofheader = (char *) NULL;
  new->contents = contents;

#ifdef GSNAP
  new->contents_uc = (char *) NULL;
  new->chop = (char *) NULL;
  new->choplength = 0;
  new->quality = new->quality_alloc = (char *) NULL;
#endif

  new->trimstart = 0;
  new->trimend = new->fulllength = length;
#ifdef PMAP
  new->fulllength_given = length;
#endif

  /* Called only by Genome_get_segment, which provides
     its own buffer */
  /* new->free_contents_p = false; */
  new->contents_alloc = (char *) NULL;
#ifdef GSNAP
  new->contents_uc_alloc = (char *) NULL;
#endif

  new->subseq_offset = 0;
  new->skiplength = 0;
  return new;
}

static char complCode[128] = COMPLEMENT_LC;

static char *
make_complement (char *sequence, unsigned int length) {
  char *complement;
  int i, j;

  complement = (char *) CALLOC(length+1,sizeof(char));
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

  complement = (char *) CALLOC(length+1,sizeof(char));
  for (i = length-1, j = 0; i >= 0; i--, j++) {
    complement[j] = uppercaseCode[(int) complCode[(int) sequence[i]]];
  }
  complement[length] = '\0';
  return complement;
}
#endif


static void
make_complement_buffered (char *complement, char *sequence, unsigned int length) {
  int i, j;

  /* complement = (char *) CALLOC(length+1,sizeof(char)); */
  for (i = length-1, j = 0; i >= 0; i--, j++) {
    complement[j] = complCode[(int) sequence[i]];
  }
  complement[length] = '\0';
  return;
}

#ifdef GSNAP
static char *
make_reverse (char *sequence, unsigned int length) {
  char *reverse;
  int i, j;

  if (sequence == NULL) {
    return (char *) NULL;
  } else {
    reverse = (char *) CALLOC(length+1,sizeof(char));
    for (i = length-1, j = 0; i >= 0; i--, j++) {
      reverse[j] = sequence[i];
    }
    reverse[length] = '\0';
    return reverse;
  }
}
#endif


/************************************************************************
 *  Original:
 *   TTTTTT ACGT ...... ACGT AAAAAA
 *          ^trimstart     ^trimend
 *   ^contents
 ************************************************************************
 *  Subsequence:
 *       ^start                ^end
 *          ^trimstart     ^trimend
 *       ^contents
 ************************************************************************/

T
Sequence_subsequence (T this, int start, int end) {
  T new;

#ifdef PMAP
  start /= 3;
  end /= 3;
#endif

  if (start < 0) {
    start = 0;
  }
  if (end > this->fulllength) {
    end = this->fulllength;
  }

  if (end <= start) {
    return NULL;
  } else {
    new = (T) MALLOC(sizeof(*new));

    new->acc = (char *) NULL;
    new->restofheader = (char *) NULL;
    new->contents = &(this->contents[start]); 

#ifdef GSNAP
    new->contents_uc = (char *) NULL;
    new->chop = (char *) NULL;
    new->choplength = 0;
    new->quality = new->quality_alloc = (char *) NULL;
#endif

    new->fulllength = end - start;
#ifdef PMAP
    new->fulllength_given = new->fulllength;
#endif
    if ((new->trimstart = this->trimstart - start) < 0) {
      new->trimstart = 0;
    }
    if ((new->trimend = this->trimend - start) > new->fulllength) {
      new->trimend = new->fulllength;
    }

    /* new->free_contents_p = false; */
    new->contents_alloc = (char *) NULL;
#ifdef GSNAP
    new->contents_uc_alloc = (char *) NULL;
#endif

    new->subseq_offset = start;
    new->skiplength = this->skiplength;
    return new;
  }
}


T
Sequence_revcomp (T this) {
  T new = (T) MALLOC(sizeof(*new));

  new->acc = (char *) NULL;
  new->restofheader = (char *) NULL;
  new->contents = new->contents_alloc = make_complement(this->contents,this->fulllength);

#ifdef GSNAP
  new->contents_uc = new->contents_uc_alloc = (char *) NULL;
  new->chop = (char *) NULL;
  new->choplength = 0;
  new->quality = new->quality_alloc = make_reverse(this->quality,this->fulllength);
#endif

  new->fulllength = this->fulllength;
#ifdef PMAP
  new->fulllength_given = this->fulllength_given;
#endif
  new->trimstart = this->trimstart;
  new->trimend = this->trimend;
  /* new->free_contents_p = true; */
  new->subseq_offset = 0;	/* Not sure if this is right */
  new->skiplength = this->skiplength;
  return new;
}

static char *
make_uppercase (char *sequence, unsigned int length) {
  char *uppercase;
#ifdef PMAP
  char uppercaseCode[128] = UPPERCASE_STD;
#else
  char uppercaseCode[128] = UPPERCASE_U2T;
#endif
  int i;

  uppercase = (char *) CALLOC(length+1,sizeof(char));
  for (i = 0; i < length; i++) {
    uppercase[i] = uppercaseCode[(int) sequence[i]];
  }
  uppercase[length] = '\0';
  return uppercase;
}

#if 0
static void
make_uppercase_inplace (char *sequence, unsigned int length) {
#ifdef PMAP
  char uppercaseCode[128] = UPPERCASE_STD;
#else
  char uppercaseCode[128] = UPPERCASE_U2T;
#endif
  int i;

  for (i = 0; i < length; i++) {
    sequence[i] = uppercaseCode[(int) sequence[i]];
  }
  return;
}
#endif


T
Sequence_uppercase (T this) {
  T new = (T) MALLOC(sizeof(*new));

  new->acc = (char *) NULL;
  new->restofheader = (char *) NULL;
  new->contents = new->contents_alloc = make_uppercase(this->contents,this->fulllength);

#ifdef GSNAP
  new->contents_uc = new->contents_uc_alloc = (char *) NULL;
  new->chop = (char *) NULL;
  new->choplength = 0;
  new->quality = new->quality_alloc = (char *) NULL;
#endif

  new->fulllength = this->fulllength;
#ifdef PMAP
  new->fulllength_given = this->fulllength_given;
#endif
  new->trimstart = this->trimstart;
  new->trimend = this->trimend;
  /* new->free_contents_p = true; */
  new->subseq_offset = this->subseq_offset;
  new->skiplength = this->skiplength;
  return new;
}

T
Sequence_alias (T this) {
  T new = (T) MALLOC(sizeof(*new));

  new->acc = (char *) NULL;
  new->restofheader = (char *) NULL;
  new->contents = this->contents;

#ifdef GSNAP
  new->contents_uc = (char *) NULL;
  new->chop = (char *) NULL;
  new->choplength = 0;
  new->quality = new->quality_alloc = (char *) NULL;
#endif

  new->fulllength = this->fulllength;
#ifdef PMAP
  new->fulllength_given = this->fulllength_given;
#endif
  new->trimstart = this->trimstart;
  new->trimend = this->trimend;

  /* new->free_contents_p = false; */
  new->contents_alloc = (char *) NULL;
#ifdef GSNAP
  new->contents_uc_alloc = (char *) NULL;
#endif

  new->subseq_offset = this->subseq_offset;
  return new;
}


/*
void
Sequence_endstream () {
  Initc = '\0';
  return;
}
*/

T
Sequence_read (int *nextchar, FILE *input, bool maponlyp) {
  T new;
  int fulllength, skiplength;
  char *pointer1, *pointer2a, *pointer2b;
  int length1, length2a, length2b;
#ifdef PMAP
  char lastchar = '*';
#endif

  if (feof(input)) {
    *nextchar = EOF;
    return NULL;
  }

  new = (T) MALLOC(sizeof(*new));

  if (*nextchar == '\0') {
    if ((*nextchar = Sequence_input_init(input)) == EOF) {
      *nextchar = EOF;
      return NULL;
    }
  }
  if (*nextchar != '>') {
    blank_header(new);
  } else if (input_header(input,new) == NULL) {
    /* File ends after >.  Don't process. */
    *nextchar = EOF;
    return NULL;
  } 
  if ((fulllength = input_sequence(&(*nextchar),&pointer1,&length1,&pointer2a,&length2a,
				   &pointer2b,&length2b,&skiplength,input)) == 0) {
    /* File ends during header.  Continue with a sequence of length 0. */
    /* fprintf(stderr,"File ends after header\n"); */
  }

  if (skiplength > 0) {
    if (maponlyp == false) {
      fprintf(stderr,"Warning: cDNA sequence length of %d exceeds maximum length of %d.  Truncating %d chars in middle.\n",
	      fulllength+skiplength,MAXSEQLEN,skiplength);
      fprintf(stderr,"  (For long sequences, perhaps you want maponly mode, by providing the '-1' flag.)\n");
    }
  }

#ifdef PMAP
  if (length1 > 0) {
    lastchar = pointer1[length1-1];
    if (length2a > 0) {
      lastchar = pointer2a[length2a-1];
    }
    if (length2b > 0) {
      lastchar = pointer2b[length2b-1];
    }
  }

  new->fulllength_given = fulllength;
  if (lastchar != '*') {
    debug(printf("Sequence does not end with *, so adding it\n"));
    fulllength++;
  }
#endif

  debug(printf("fulllength = %d\n",fulllength));
  new->fulllength = fulllength;
  new->skiplength = skiplength;

  new->trimstart = 0;
#ifdef PMAP
  new->trimend = new->fulllength_given;
#else
  new->trimend = fulllength;
#endif

  new->contents = new->contents_alloc = (char *) CALLOC(fulllength+1,sizeof(char));
  if (length1 > 0) {
    strncpy(new->contents,pointer1,length1);
    if (length2a > 0) {
      strncpy(&(new->contents[length1]),pointer2a,length2a);
    }
    if (length2b > 0) {
      strncpy(&(new->contents[length1+length2a]),pointer2b,length2b);
    }
  }

#ifdef GSNAP
  new->contents_uc = new->contents_uc_alloc = (char *) NULL;
  new->chop = (char *) NULL;
  new->choplength = 0;
  new->quality = new->quality_alloc = (char *) NULL;
#endif

#ifdef PMAP
  if (lastchar != '*') {
    new->contents[fulllength-1] = '*';
  }
#endif
  /* new->free_contents_p = true; */
  new->subseq_offset = 0;

  debug(printf("Final query sequence is:\n"));
  debug(Sequence_print(new,/*uppercasep*/false,/*wraplength*/60,/*trimmedp*/false));
  return new;
}


#if 0
static bool
Sequence_skip (int *nextchar, FILE *input, bool maponlyp) {
  int skiplength;
  char *pointer1, *pointer2a, *pointer2b;
  int length1, length2a, length2b;
#ifdef PMAP
  char lastchar = '*';
#endif

  if (feof(input)) {
    *nextchar = EOF;
    return false;
  }

  if (*nextchar == '\0') {
    if ((*nextchar = Sequence_input_init(input)) == EOF) {
      *nextchar = EOF;
      return false;
    }
  }
  if (*nextchar != '>') {
    /* Header is blank */
  } else if (skip_header(input) == false) {
    /* File ends after >.  Don't process. */
    *nextchar = EOF;
    return false;
  } 
  if (input_sequence(&(*nextchar),&pointer1,&length1,&pointer2a,&length2a,
		     &pointer2b,&length2b,&skiplength,input) == 0) {
    /* File ends during header.  Continue with a sequence of length 0. */
    /* fprintf(stderr,"File ends after header\n"); */
  }

  return true;
}
#endif


T
Sequence_read_multifile (int *nextchar, FILE **input, char ***files, int *nfiles, bool maponlyp) {
  T queryseq;

#ifdef LEAKCHECK
  Mem_leak_check_start(__FILE__,__LINE__);
#endif

  while (1) {
    if (*input == NULL || feof(*input)) {
      if (*nfiles == 0) {
	*nextchar = EOF;
	return NULL;
      } else {
	if ((*input = FOPEN_READ_TEXT((*files)[0])) == NULL) {
	  fprintf(stderr,"Can't open file %s => skipping it.\n",(*files)[0]);
	}
	(*files)++;
	(*nfiles)--;
	*nextchar = '\0';
      }
    }
    if (*input != NULL) {
      if ((queryseq = Sequence_read(&(*nextchar),*input,maponlyp)) != NULL) {
	return queryseq;
      }
    }
  }
}

#if 0
static bool
Sequence_skip_multifile (int *nextchar, FILE **input, char ***files, int *nfiles, bool maponlyp) {

  while (1) {
    if (*input == NULL || feof(*input)) {
      if (*nfiles == 0) {
	*nextchar = EOF;
	return false;
      } else {
	if ((*input = FOPEN_READ_TEXT((*files)[0])) == NULL) {
	  fprintf(stderr,"Can't open file %s => skipping it.\n",(*files)[0]);
	}
	(*files)++;
	(*nfiles)--;
	*nextchar = '\0';
      }
    }
    if (*input != NULL) {
      return Sequence_skip(&(*nextchar),*input,maponlyp);
    }
  }
}
#endif


#ifdef GSNAP

/************************************************************************
 *   Dynamic programming for paired-end short reads
 ************************************************************************/

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



void
Sequence_dynprog_init (int maxlength) {

#ifdef LEAKCHECK
  Mem_leak_check_deactivate(__FILE__,__LINE__);
#endif

  matrix_ptrs = (int **) CALLOC(maxlength+1,sizeof(int *));
  matrix_space = (int *) CALLOC((maxlength+1)*(maxlength+1),sizeof(int));
  directions_ptrs = (Direction_T **) CALLOC(maxlength+1,sizeof(Direction_T *));
  directions_space = (Direction_T *) CALLOC((maxlength+1)*(maxlength+1),sizeof(Direction_T));
  jump_ptrs = (int **) CALLOC(maxlength+1,sizeof(int *));
  jump_space = (int *) CALLOC((maxlength+1)*(maxlength+1),sizeof(int));

  pairdistance_init();
  jump_penalty_init(maxlength);

#ifdef LEAKCHECK
  Mem_leak_check_activate(__FILE__,__LINE__);
#endif

  return;
}

void
Sequence_dynprog_term () {
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


/* Makes a matrix of dimensions 0..length1 x 0..length2 inclusive */
static int **
Matrix_alloc (int length1, int length2, int **ptrs, int *space) {
  int **matrix, i;

  if (length1 < 0 || length2 < 0) {
    fprintf(stderr,"lengths are negative: %d %d\n",length1,length2);
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
      bestscore = matrix[r-1][c-1] + pairdistance_array[sequence1[r-1]][sequence2[c-1]];
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
  bestscore = matrix[length1-1][length2-1] + pairdistance_array[sequence1[length1-1]][sequence2[length2-1]];
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
traceback (int *chop1, int *chop2, Direction_T **directions, int **jump, int length1, int length2) {
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
Sequence_chop_primers (T queryseq1, T queryseq2) {
  int finalscore;
  int **matrix, **jump;
  Direction_T **directions;
  int chop1, chop2;

  matrix = compute_scores_lookup(&directions,&jump,queryseq1->contents_uc,queryseq2->contents_uc,
				 queryseq1->fulllength,queryseq2->fulllength);
  finalscore = matrix[queryseq1->fulllength][queryseq2->fulllength];

  if (finalscore <= 0) {
    chop1 = chop2 = 0;
  } else {
    traceback(&chop1,&chop2,directions,jump,queryseq1->fulllength,queryseq2->fulllength);
    if (chop1 > 0) {
      queryseq1->chop = (char *) CALLOC(chop1+1,sizeof(char));
      strncpy(queryseq1->chop,&(queryseq1->contents[queryseq1->fulllength - chop1]),chop1);
      queryseq1->choplength = chop1;

      queryseq1->contents[queryseq1->fulllength - chop1] = '\0';
      queryseq1->contents_uc[queryseq1->fulllength - chop1] = '\0';
      if (queryseq1->quality != NULL) {
	queryseq1->quality[queryseq1->fulllength - chop1] = '\0';
      }

      queryseq1->fulllength -= chop1;
    }
    if (chop2 > 0) {
      queryseq2->chop = (char *) CALLOC(chop2+1,sizeof(char));
      strncpy(queryseq2->chop,queryseq2->contents,chop2);
      queryseq2->choplength = chop2;

      queryseq2->contents = &(queryseq2->contents[chop2]);
      queryseq2->contents_uc = &(queryseq2->contents_uc[chop2]);
      if (queryseq2->quality != NULL) {
	queryseq2->quality = &(queryseq2->quality[chop2]);
      }

      queryseq2->fulllength -= chop2;
    }
  }

  debug2(printf("finalscore = %d, chop1 = %d, chop2 = %d\n",finalscore,chop1,chop2));

  return;
}

#endif


/************************************************************************
 *   Short reads
 ************************************************************************/

#ifdef GSNAP
T
Sequence_read_multifile_shortreads (int *nextchar, T *queryseq2, FILE **input, char ***files, int *nfiles,
				    bool circularp) {
  T queryseq1, temp;
  int fulllength;

#ifdef LEAKCHECK
  Mem_leak_check_start(__FILE__,__LINE__);
#endif

  while (1) {
    if (*input == NULL || feof(*input)) {
      if (*nfiles == 0) {
	*nextchar = EOF;
	return (T) NULL;
      } else {
	if ((*input = FOPEN_READ_TEXT((*files)[0])) == NULL) {
	  fprintf(stderr,"Can't open file %s => skipping it.\n",(*files)[0]);
	}
	(*files)++;
	(*nfiles)--;
	*nextchar = '\0';
      }
    }
    if (*input != NULL) {
      if (*nextchar == '\0') {
	if ((*nextchar = Sequence_input_init(*input)) == EOF) {
	  *nextchar = EOF;
	  return NULL;
	}
      }

      queryseq1 = (T) MALLOC(sizeof(*queryseq1));
      (*queryseq2) = (T) NULL;
      debug(printf("** Getting header\n"));
      if (input_header(*input,queryseq1) == NULL) {
	/* fprintf(stderr,"No header\n"); */
	/* File ends after >.  Don't process. */
	*nextchar = EOF;
      } else {
	queryseq1->chop = (char *) NULL;
	queryseq1->choplength = 0;
	queryseq1->quality = queryseq1->quality_alloc = (char *) NULL;
	if ((*nextchar = fgetc(*input)) == '\n') {
	  queryseq1->contents = (char *) NULL;
	  queryseq1->fulllength = 0;
	  /* queryseq1->free_contents_p = false; */
	  queryseq1->contents_alloc = (char *) NULL;
	  queryseq1->contents_uc_alloc = (char *) NULL;
	  while (*nextchar != EOF && ((*nextchar = fgetc(*input)) != '>')) {
	  }
	} else if ((queryseq1->fulllength = input_oneline(&(*nextchar),*input,/*possible_fasta_header_p*/true)) == 0) {
	  /* fprintf(stderr,"length is zero\n"); */
	  /* No sequence1.  Don't process. */
	  /* *nextchar = EOF; */
	  queryseq1->contents = (char *) NULL;
	  /* queryseq1->free_contents_p = false; */
	  queryseq1->contents_alloc = (char *) NULL;
	  queryseq1->contents_uc_alloc = (char *) NULL;
	} else {
	  queryseq1->contents = queryseq1->contents_alloc = (char *) CALLOC(queryseq1->fulllength+1,sizeof(char));
	  strncpy(queryseq1->contents,&(Firsthalf[0]),queryseq1->fulllength);
	  queryseq1->contents_uc = queryseq1->contents_uc_alloc = make_uppercase(queryseq1->contents,queryseq1->fulllength);
	  /* queryseq1->free_contents_p = true; */
	  if ((fulllength = input_oneline(&(*nextchar),*input,/*possible_fasta_header_p*/true)) > 0) {
	    (*queryseq2) = (T) MALLOC(sizeof(*(*queryseq2)));
	    (*queryseq2)->acc = (char *) NULL;
	    (*queryseq2)->restofheader = (char *) NULL;
	    (*queryseq2)->fulllength = fulllength;
	    (*queryseq2)->contents = (*queryseq2)->contents_alloc = make_complement(&(Firsthalf[0]),fulllength);
	    (*queryseq2)->contents_uc = (*queryseq2)->contents_uc_alloc = make_uppercase((*queryseq2)->contents,fulllength);
	    (*queryseq2)->chop = (char *) NULL;
	    (*queryseq2)->choplength = 0;
	    (*queryseq2)->quality = (*queryseq2)->quality_alloc = (char *) NULL;
	    /* (*queryseq2)->free_contents_p = true; */
	    if (circularp == true) {
	      /* Circular-end read.  Swap queryseq1 and queryseq2. */
	      temp = queryseq1;
	      queryseq1 = *queryseq2;
	      *queryseq2 = temp;
	    }
	  }
	}
	debug(printf("Returning queryseq with contents %s\n",queryseq1->contents));
	return queryseq1;
      }
    }
  }
}
#endif

#ifdef GSNAP
/* Does not do reverse complement */
T
Sequence_read_multifile_shortreads_simple (int *nextchar, T *queryseq2, FILE **input, char ***files, int *nfiles) {
  T queryseq1;
  int fulllength;

  while (1) {
    if (*input == NULL || feof(*input)) {
      if (*nfiles == 0) {
	*nextchar = EOF;
	return (T) NULL;
      } else {
	if ((*input = FOPEN_READ_TEXT((*files)[0])) == NULL) {
	  fprintf(stderr,"Can't open file %s => skipping it.\n",(*files)[0]);
	}
	(*files)++;
	(*nfiles)--;
	*nextchar = '\0';
      }
    }
    if (*input != NULL) {
      if (*nextchar == '\0') {
	if ((*nextchar = Sequence_input_init(*input)) == EOF) {
	  *nextchar = EOF;
	  return NULL;
	}
      }

      queryseq1 = (T) MALLOC(sizeof(*queryseq1));
      (*queryseq2) = (T) NULL;
      debug(printf("** Getting header\n"));
      if (input_header(*input,queryseq1) == NULL) {
	/* fprintf(stderr,"No header\n"); */
	/* File ends after >.  Don't process. */
	*nextchar = EOF;
      } else {
	queryseq1->chop = (char *) NULL;
	queryseq1->choplength = 0;
	queryseq1->quality = queryseq1->quality_alloc = (char *) NULL;

	if ((*nextchar = fgetc(*input)) == '\n') {
	  queryseq1->contents = (char *) NULL;
	  queryseq1->fulllength = 0;
	  /* queryseq1->free_contents_p = false; */
	  queryseq1->contents_alloc = (char *) NULL;
	  queryseq1->contents_uc_alloc = (char *) NULL;
	  while (*nextchar != EOF && ((*nextchar = fgetc(*input)) != '>')) {
	  }
	} else if ((queryseq1->fulllength = input_oneline(&(*nextchar),*input,/*possible_fasta_header_p*/true)) == 0) {
	  /* fprintf(stderr,"length is zero\n"); */
	  /* No sequence1.  Don't process. */
	  /* *nextchar = EOF; */
	  queryseq1->contents = (char *) NULL;
	  /* queryseq1->free_contents_p = false; */
	  queryseq1->contents_alloc = (char *) NULL;
	  queryseq1->contents_uc_alloc = (char *) NULL;
	} else {
	  queryseq1->contents = queryseq1->contents_alloc = (char *) CALLOC(queryseq1->fulllength+1,sizeof(char));
	  strncpy(queryseq1->contents,&(Firsthalf[0]),queryseq1->fulllength);
	  queryseq1->contents_uc = queryseq1->contents_uc_alloc = make_uppercase(queryseq1->contents,queryseq1->fulllength);
	  /* queryseq1->free_contents_p = true; */
	  if ((fulllength = input_oneline(&(*nextchar),*input,/*possible_fasta_header_p*/true)) > 0) {
	    (*queryseq2) = (T) MALLOC(sizeof(*(*queryseq2)));
	    (*queryseq2)->acc = (char *) NULL;
	    (*queryseq2)->restofheader = (char *) NULL;
	    (*queryseq2)->fulllength = fulllength;
	    (*queryseq2)->contents = (*queryseq2)->contents_alloc = (char *) CALLOC((*queryseq2)->fulllength+1,sizeof(char));
	    strncpy((*queryseq2)->contents,&(Firsthalf[0]),(*queryseq2)->fulllength);
	    (*queryseq2)->contents_uc = (*queryseq2)->contents_uc_alloc = make_uppercase((*queryseq2)->contents,(*queryseq2)->fulllength);
	    (*queryseq2)->chop = (char *) NULL;
	    (*queryseq2)->choplength = 0;
	    (*queryseq2)->quality = (*queryseq2)->quality_alloc = (char *) NULL;
	    /* (*queryseq2)->free_contents_p = true; */
	  }
	}
	debug(printf("Returning queryseq with contents %s\n",queryseq1->contents));
	return queryseq1;
      }
    }
  }
}
#endif


#ifdef GSNAP
T
Sequence_read_fastq_shortreads (int *nextchar, T *queryseq2, FILE *input1, FILE *input2,
				bool circularp) {
  T queryseq1, temp;
  int nextchar2;
  int qualitylength;

#ifdef LEAKCHECK
  Mem_leak_check_start(__FILE__,__LINE__);
#endif

  while (1) {
    if (input1 == NULL || feof(input1)) {
      *nextchar = EOF;
      return (T) NULL;
    }

    queryseq1 = (T) MALLOC(sizeof(*queryseq1));
    debug(printf("** Getting header\n"));
    if (input_header_fastq(input1,queryseq1) == NULL) {
      /* fprintf(stderr,"No header\n"); */
      /* File ends after >.  Don't process. */
      *nextchar = EOF;
    } else {
      *nextchar = fgetc(input1);
      if ((queryseq1->fulllength = input_oneline(&(*nextchar),input1,/*possible_fasta_header_p*/true)) == 0) {
	/* fprintf(stderr,"length is zero\n"); */
	/* No sequence1.  Don't process. */
	/* *nextchar = EOF; */
	queryseq1->contents = (char *) NULL;
	/* queryseq1->free_contents_p = false; */
	queryseq1->contents_alloc = (char *) NULL;
	queryseq1->contents_uc_alloc = (char *) NULL;

      } else {
	queryseq1->contents = queryseq1->contents_alloc = (char *) CALLOC(queryseq1->fulllength+1,sizeof(char));
	strncpy(queryseq1->contents,&(Firsthalf[0]),queryseq1->fulllength);
	queryseq1->contents_uc = queryseq1->contents_uc_alloc = make_uppercase(queryseq1->contents,queryseq1->fulllength);
	/* queryseq1->free_contents_p = true; */
      }

      queryseq1->chop = (char *) NULL;
      queryseq1->choplength = 0;

      if (*nextchar != '+') {
	queryseq1->quality = queryseq1->quality_alloc = (char *) NULL;
      } else {
	skip_header(input1);
	*nextchar = fgetc(input1);
	qualitylength = input_oneline(&(*nextchar),input1,/*possible_fasta_header_p*/false);
	if (qualitylength != queryseq1->fulllength) {
	  fprintf(stderr,"Length %d of quality score differs from length %d of nucleotides in sequence %s\n",
		  qualitylength,queryseq1->fulllength,queryseq1->acc);
	  abort();
	} else {
	  queryseq1->quality = queryseq1->quality_alloc = (char *) CALLOC(qualitylength+1,sizeof(char));
	  strncpy(queryseq1->quality,&(Firsthalf[0]),qualitylength);
	}
      }
    }

    if (input2 == NULL) {
      *queryseq2 = (T) NULL;
    } else {
      *queryseq2 = (T) MALLOC(sizeof(*(*queryseq2)));
      if (input_header_fastq(input2,*queryseq2) == NULL) {
	/* fprintf(stderr,"No header\n"); */
	/* File ends after >.  Don't process. */
	nextchar2 = EOF;
      } else {
	strip_illumina_acc_ending(queryseq1->acc);
	strip_illumina_acc_ending((*queryseq2)->acc);

	if (strcmp(queryseq1->acc,(*queryseq2)->acc)) {
	  fprintf(stderr,"Paired-end accessions %s and %s do not match\n",queryseq1->acc,(*queryseq2)->acc);
	  exit(9);
	}
	nextchar2 = fgetc(input2);
	if (((*queryseq2)->fulllength = input_oneline(&nextchar2,input2,/*possible_fasta_header_p*/true)) == 0) {
	  /* fprintf(stderr,"length is zero\n"); */
	  /* No sequence1.  Don't process. */
	  /* *nextchar = EOF; */
	  (*queryseq2)->contents = (char *) NULL;
	  /* (*queryseq2)->free_contents_p = false; */
	  (*queryseq2)->contents_alloc = (char *) NULL;
	  (*queryseq2)->contents_uc_alloc = (char *) NULL;

	} else {
	  (*queryseq2)->contents = (*queryseq2)->contents_alloc = make_complement(&(Firsthalf[0]),(*queryseq2)->fulllength);
	  (*queryseq2)->contents_uc = (*queryseq2)->contents_uc_alloc = make_uppercase((*queryseq2)->contents,(*queryseq2)->fulllength);
	  /* (*queryseq2)->free_contents_p = true; */
	  if (circularp == true) {
	    /* Circular-end read.  Swap queryseq1 and queryseq2. */
	    temp = queryseq1;
	    queryseq1 = *queryseq2;
	    *queryseq2 = temp;
	  }
	}

	(*queryseq2)->chop = (char *) NULL;
	(*queryseq2)->choplength = 0;

	if (nextchar2 != '+') {
	  (*queryseq2)->quality = (*queryseq2)->quality_alloc = (char *) NULL;
	} else {
	  skip_header(input2);
	  nextchar2 = fgetc(input2);
	  qualitylength = input_oneline(&nextchar2,input2,/*possible_fasta_header_p*/false);
	  if (qualitylength != (*queryseq2)->fulllength) {
	    fprintf(stderr,"Length %d of quality score differs from length %d of nucleotides in sequence %s\n",
		    qualitylength,(*queryseq2)->fulllength,(*queryseq2)->acc);
	    abort();
	  } else {
	    (*queryseq2)->quality = (*queryseq2)->quality_alloc = make_reverse(&(Firsthalf[0]),qualitylength);
	  }
	}
      }
    }

    return queryseq1;
  }
}
#endif




#if 0
static bool
Sequence_skip_multifile_shortreads (int *nextchar, FILE **input, char ***files, int *nfiles) {

  while (1) {
    if (*input == NULL || feof(*input)) {
      if (*nfiles == 0) {
	*nextchar = EOF;
	return false;
      } else {
	if ((*input = FOPEN_READ_TEXT((*files)[0])) == NULL) {
	  fprintf(stderr,"Can't open file %s => skipping it.\n",(*files)[0]);
	}
	(*files)++;
	(*nfiles)--;
	*nextchar = '\0';
      }
    }
    if (*input != NULL) {
      if (*nextchar == '\0') {
	if ((*nextchar = Sequence_input_init(*input)) == EOF) {
	  *nextchar = EOF;
	  return false;
	}
      }

      if (skip_header(*input) == false) {
	/* fprintf(stderr,"No header\n"); */
	/* File ends after >.  Don't process. */
	*nextchar = EOF;
      } else {
	*nextchar = fgetc(*input);
	if (input_oneline(&(*nextchar),*input,/*possible_fasta_header_p*/true) == 0) {
	  /* fprintf(stderr,"length is zero\n"); */
	  /* No sequence1.  Don't process. */
	  *nextchar = EOF;
	} else {
	  /* queryseq2, if present */
	  input_oneline(&(*nextchar),*input,/*possible_fasta_header_p*/true);
	  return true;
	}
      }
    }
  }
}
#endif

T
Sequence_read_unlimited (FILE *input) {
  T new;
  Intlist_T intlist = NULL;
  char *p;
  int length, startpos = 1, maxseqlen = MAXSEQLEN;
  bool eolnp;

  if (feof(input)) {
    return NULL;
  }

  new = (T) MALLOC(sizeof(*new));

  if (Initc == '\0') {
    if ((Initc = Sequence_input_init(input)) == EOF) {
      return NULL;
    }
  }
  if (Initc != '>') {
    blank_header(new);
    Sequence[startpos++] = Initc;
    maxseqlen--;
  } else if (input_header(input,new) == NULL) {
    /* File ends after >.  Don't process. */
     return NULL;
  } 
  /* Don't touch Sequence[0], because subsequent calls to
     Sequence_read depend on it being '\0'. */
  eolnp = true;
  while (fgets(&(Sequence[startpos]),maxseqlen,input) != NULL &&
	 (eolnp == false || Sequence[1] != '>')) {
    for (p = &(Sequence[1]); *p != '\n' && *p != '\0'; p++) {
      if (!iscntrl((int) *p)
#ifdef DASH
	  && *p != DASH
#endif
	  ) {
	intlist = Intlist_push(intlist,(int) *p);
      }
    }
    if (*p == '\n') {
      eolnp = true;
    } else {
      eolnp = false;
    }
    startpos = 1;
    maxseqlen = MAXSEQLEN;
  }
  intlist = Intlist_reverse(intlist);
  new->contents = new->contents_alloc = Intlist_to_char_array(&length,intlist);

#ifdef GSNAP
  new->contents_uc = new->contents_uc_alloc = (char *) NULL;
  new->chop = (char *) NULL;
  new->choplength = 0;
  new->quality = new->quality_alloc = (char *) NULL;
#endif

  Intlist_free(&intlist);

  if (length == 0) {
    return NULL;
  } else {
    new->fulllength = new->trimend = length;
#ifdef PMAP
    new->fulllength_given = length;
#endif
    new->trimstart = 0;

    /* new->free_contents_p = true; */
    new->subseq_offset = 0;
    new->skiplength = 0;

    /* Important to initialize for subsequent cDNA reads */
    Initc = '\0';

    return new;
  }
}

void
Sequence_print_digest (T this) {
  unsigned char *digest;

  digest = MD5_compute((unsigned char *) this->contents,this->fulllength);
  MD5_print(digest);
  FREE(digest);
  return;
}

/* Calling procedure needs to print the initial ">", if desired */
void
Sequence_print_header (T this, bool checksump) {

#if 0
#ifdef PMAP
  printf("%s (%d aa) %s",this->acc,this->fulllength_given+this->skiplength,this->restofheader);
#else
  printf("%s (%d bp) %s",this->acc,this->fulllength+this->skiplength,this->restofheader);
#endif
#else
  printf("%s %s",this->acc,this->restofheader);
#endif
  if (checksump == true) {
    printf(" md5:");
    Sequence_print_digest(this);
  }
  printf("\n");
  return;
}

void
Sequence_print_header_revcomp (T this) {
  printf(">%s %s",this->acc,this->restofheader);
  printf(" REVCOMP");
  printf("\n");
  return;
}


void
Sequence_print (T this, bool uppercasep, int wraplength, bool trimmedp) {
  int i = 0, pos, start, end;
  char uppercaseCode[128] = UPPERCASE_STD;

  if (trimmedp == true) {
    start = this->trimstart;
    end = this->trimend;
  } else {
    start = 0;
    end = this->fulllength;
  }

  if (uppercasep == true) {
    for (pos = start; pos < end; pos++, i++) {
      printf("%c",uppercaseCode[(int) this->contents[i]]);
      if ((i+1) % wraplength == 0) {
	printf("\n");
      }
    }
  } else {
    for (pos = start; pos < end; pos++, i++) {
      printf("%c",this->contents[i]);
      if ((i+1) % wraplength == 0) {
	printf("\n");
      }
    }
  }
  if (i % wraplength != 0) {
    printf("\n");
  }
  return;
}

void
Sequence_print_two (T this, T alt, bool uppercasep, int wraplength) {
  int i = 0, pos, pos2, startpos, end;
  char uppercaseCode[128] = UPPERCASE_STD;

  end = this->fulllength;

  pos = 0;
  i = 0;
  if (uppercasep == true) {
    printf("ref\t");
    startpos = pos;
    while (pos < end) {
      printf("%c",uppercaseCode[(int) this->contents[pos]]);
      if (++i % wraplength == 0) {
	printf("\n");
	printf("alt\t");
	for (pos2 = startpos, i = 0; i < wraplength; pos2++, i++) {
	  printf("%c",uppercaseCode[(int) alt->contents[pos2]]);
	}
	printf("\n\n");
	if (pos+1 < end) {
	  printf("ref\t");
	}
	startpos = pos+1;
      }
      pos++;
    }
    if (i % wraplength != 0) {
      printf("\n");
      printf("alt\t");
      for (pos2 = startpos; pos2 < end; pos2++) {
	printf("%c",uppercaseCode[(int) alt->contents[pos2]]);
      }
      printf("\n\n");
    }

  } else {
    printf("ref\t");
    startpos = pos;
    while (pos < end) {
      printf("%c",this->contents[pos]);
      if (++i % wraplength == 0) {
	printf("\n");
	printf("alt\t");
	for (pos2 = startpos, i = 0; i < wraplength; pos2++, i++) {
	  printf("%c",alt->contents[pos2]);
	}
	printf("\n\n");
	if (pos+1 < end) {
	  printf("ref\t");
	}
	startpos = pos+1;
      }
      pos++;
    }
    if (i % wraplength != 0) {
      printf("\n");
      printf("alt\t");
      for (pos2 = startpos; pos2 < end; pos2++) {
	printf("%c",alt->contents[pos2]);
      }
      printf("\n\n");
    }
  }

  return;
}



#ifdef GSNAP
void
Sequence_print_chop_symbols (FILE *fp, T this) {
  int i;

  for (i = 0; i < this->choplength; i++) {
    fprintf(fp,"*");
  }
  return;
}


void
Sequence_print_oneline (FILE *fp, T this) {
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
Sequence_print_oneline_revcomp (FILE *fp, T this) {
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
Sequence_print_chopped (FILE *fp, T this, int hardclip_low, int hardclip_high) {
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
Sequence_print_chopped_revcomp (FILE *fp, T this, int hardclip_low, int hardclip_high) {
  int i;

  for (i = this->fulllength - 1 - hardclip_low; i >= hardclip_high; --i) {
    fprintf(fp,"%c",complCode[(int) this->contents[i]]);
  }

  return;
}

void
Sequence_print_quality (FILE *fp, T this, int hardclip_low, int hardclip_high,
			int shift) {
  int i;

  if (this->quality == NULL) {
    printf("*");
  } else {
    for (i = hardclip_low; i < this->fulllength - hardclip_high; i++) {
      fprintf(fp,"%c",this->quality[i] + shift);
    }
  }
  return;
}

void
Sequence_print_quality_revcomp (FILE *fp, T this, int hardclip_low, int hardclip_high,
				int shift) {
  int i;

  if (this->quality == NULL) {
    printf("*");
  } else {
    for (i = this->fulllength - 1 - hardclip_low; i >= hardclip_high; --i) {
      fprintf(fp,"%c",this->quality[i] + shift);
    }
  }

  return;
}

void
Sequence_print_oneline_uc (FILE *fp, T this) {
  int i = 0;

  for (i = 0; i < this->fulllength; i++) {
    fprintf(fp,"%c",this->contents_uc[i]);
  }
  return;
}

void
Sequence_print_oneline_revcomp_uc (FILE *fp, T this) {
  int i = 0;

  for (i = this->fulllength-1; i >= 0; --i) {
    fprintf(fp,"%c",complCode[(int) this->contents_uc[i]]);
  }
  return;
}
#endif

void
Sequence_print_raw (T this) {
  int i = 0, pos, start, end;

  start = 0;
  end = this->fulllength;

  for (pos = start; pos < end; pos++, i++) {
    printf("%d\n",(int) this->contents[i]);
  }
  return;
}


T
Sequence_substring (T usersegment, unsigned int left, unsigned int length, 
		    bool revcomp, char *gbuffer1, char *gbuffer2, int gbufferlen) {
  if (length > gbufferlen) {
    fprintf(stderr,"Didn't allocate enough space for gbufferlen (%d < %u)\n",
	    gbufferlen,length);
    abort();
  }

  memcpy(gbuffer1,&(usersegment->contents[left]),length*sizeof(char));
  gbuffer1[length] = '\0';

  if (revcomp == true) {
    make_complement_buffered(gbuffer2,gbuffer1,length);
    debug(fprintf(stderr,"Got sequence at %u with length %u, revcomp\n",left,length));
    return Sequence_genomic_new(gbuffer2,length);
  } else {
    debug(fprintf(stderr,"Got sequence at %u with length %u, forward\n",left,length));
    return Sequence_genomic_new(gbuffer1,length);
  }
}

