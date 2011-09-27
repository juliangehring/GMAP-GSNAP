static char rcsid[] = "$Id: sequence.c,v 1.52 2005/07/08 14:42:23 twu Exp $";
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
#include <ctype.h>		/* For toupper */
#include "mem.h"
#include "complement.h"
#include "intlist.h"
#include "md5.h"

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
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
  int trimstart;		/* Start of trim */
  int trimend;			/* End of trim */
  int fulllength;		/* Full length */
  bool free_contents_p;
};

char *
Sequence_accession (T this) {
  return this->acc;
}

char *
Sequence_fullpointer (T this) {
  return this->contents;
}

char *
Sequence_trimpointer (T this) {
  return &(this->contents[this->trimstart]);
}

int
Sequence_fulllength (T this) {
  return this->fulllength;
}

int
Sequence_trimlength (T this) {
  return this->trimend - this->trimstart;
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


void
Sequence_free (T *old) {
  if (*old) {
    if ((*old)->restofheader != NULL) {
      FREE((*old)->restofheader);
    }
    if ((*old)->acc != NULL) {
      FREE((*old)->acc);
    }
    if ((*old)->free_contents_p == true) {
      FREE((*old)->contents);
    }
    FREE(*old);
  }
  return;
}


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


#define HEADERLEN 512
#define DISCARDLEN 8192
#define SEQUENCELEN MAXSEQLEN+2	/* extra spaces for beginning and end */

static char Header[HEADERLEN];
static char Sequence[SEQUENCELEN];
static char Discard[DISCARDLEN];

static char *startinit;
static char *endinit;
static int Initc = '\0';


/* The first element of Sequence is always the null character, to mark
   the end of the string */

/* Returns '>' if FASTA file, first sequence char if not */
static int
input_init (FILE *fp) {
  Header[0] = '\0';
  Sequence[0] = '\0';
  startinit = &(Sequence[1]);

  return fgetc(fp);
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

static void
print_contents (char *p, int length) {
  int i;
	
  fprintf(stderr,"\"");
  for (i = 0; i < length; i++) {
    if (*p == '\0') {
      fprintf(stderr,"_");
    } else {
      fprintf(stderr,"%c",*p);
    }
    p++;
  }
  fprintf(stderr,"\"\n");
  return;
}


/* Returns 1 if done reading sequence, 0 if not */
static bool
read_one_sequence (int *nextchar, bool *eolnp, FILE *fp) {
  int remainder;
  char *sequence;
  char *p;
  int c;

  sequence = &(Sequence[1]);
  if (Initc != '>') {
    *sequence++ = (char) Initc;
  }
  remainder = (&(Sequence[SEQUENCELEN]) - sequence)/sizeof(char);

  while (1) {
    if (remainder <= 1) {
      debug(fprintf(stderr,"remainder <= 1.  Returning false\n"));
      *nextchar = EOF;
      return false;

    } else if (feof(fp)) {
      /* EOF in middle of line */
      *sequence++ = '\0';
      debug(fprintf(stderr,"EOF.  Returning true\n"));
      *nextchar = EOF;
      return true;

    } else if (*eolnp == true) {
      if ((c = fgetc(fp)) == EOF || c == '>') {
	*sequence++ = '\0';
	debug(fprintf(stderr,"c == EOF or >.  Returning true\n"));
	*nextchar = c;
	return true;
      } else if (iscntrl(c)) {
	debug(fprintf(stderr,"c == control char.  Continuing\n"));
      } else if (isspace(c)) {
	*eolnp = true;
	debug(fprintf(stderr,"c == NULL.  Continuing\n"));
      } else {
	*sequence++ = (char) c;
	remainder--;
	*eolnp = false;
	debug(fprintf(stderr,"c == sth.  Continuing\n"));
      }

    } else {
      if (fgets(sequence,remainder,fp) == NULL) {
	*sequence++ = '\0';
	debug(fprintf(stderr,"line == NULL.  Returning true\n"));
	*nextchar = EOF;
	return true;
      } else {
	if ((p = rindex(sequence,13)) != NULL) {
	  /* Handle PC line feed ^M */
	  *p = '\0';
	  *eolnp = true;
	  debug(fprintf(stderr,"line == EOLN.  Continuing\n"));
	} else if ((p = rindex(sequence,'\n')) != NULL) {
	  *p = '\0';
	  *eolnp = true;
	  debug(fprintf(stderr,"line == EOLN.  Continuing\n"));
	} else {
	  p = rindex(sequence,'\0');
	  *eolnp = false;
	  debug(fprintf(stderr,"line != EOLN.  Continuing\n"));
	}
	sequence = p;
	remainder = (&(Sequence[SEQUENCELEN]) - sequence)/sizeof(char);
      }
    }

    debug(print_contents(&(Sequence[0]),SEQUENCELEN));
  }
}

static int
discard_one_sequence (int *nextchar, bool eolnp, FILE *fp) {
  int ncycles = 0;
  int remainder;
  char *discard;
  char *p;
  int c;
  
  discard = &(Discard[0]);
  remainder = (&(Discard[DISCARDLEN]) - discard)/sizeof(char);

  while (1) {
    debug(fprintf(stderr,"\nEnd: %d\n",remainder));

    if (feof(fp)) {
      debug(fprintf(stderr,"EOF.  Returning\n"));
      *nextchar = EOF;
      return ncycles*DISCARDLEN + (discard - &(Discard[0]))/sizeof(char);

    } else if (remainder <= 1) {
      discard = &(Discard[0]);
      remainder = (&(Discard[DISCARDLEN]) - discard)/sizeof(char);
      debug(fprintf(stderr,"remainder <= 1.  Cycling\n"));

    } else if (eolnp == true) {
      if ((c = fgetc(fp)) == EOF || c == '>') {
	debug(fprintf(stderr,"c == EOF or >.  Returning\n"));
	*nextchar = c;
	return ncycles*DISCARDLEN + (discard - &(Discard[0]))/sizeof(char);
      } else if (iscntrl(c)) {
	debug(fprintf(stderr,"c == control char.  Continuing\n"));
      } else if (isspace(c)) {
	debug(fprintf(stderr,"c == NULL.  Continuing\n"));
      } else {
	*discard++ = (char) c;
	remainder--;
	eolnp = false;
	debug(fprintf(stderr,"c == sth.  Continuing\n"));
      }
      
    } else {

      if (fgets(discard,remainder,fp) == NULL) {
	debug(fprintf(stderr,"line == NULL.  Returning\n"));
	*nextchar = EOF;
	return ncycles*DISCARDLEN + (discard - &(Discard[0]))/sizeof(char);
      } else {
	if ((p = rindex(discard,'\n')) != NULL) {
	  *p = '\0';
	  eolnp = true;
	  debug(fprintf(stderr,"line == EOLN.  Continuing\n"));
	} else {
	  p = rindex(discard,'\0');
	  eolnp = false;
	  debug(fprintf(stderr,"line != EOLN.  Continuing\n"));
	}
	discard = p;
	remainder = (&(Discard[DISCARDLEN]) - discard)/sizeof(char);
      }
    }

    debug(print_contents(&(Discard[0]),DISCARDLEN));
  }
}

/* Returns sequence length */
static int
input_sequence (int *nextchar, int *length, int *skiplength, FILE *fp) {
  bool eolnp = true;

  if (read_one_sequence(&(*nextchar),&eolnp,fp) == true) {
    *skiplength = 0;
  } else {
    *skiplength = discard_one_sequence(&(*nextchar),eolnp,fp);
  }
  endinit = rindex(startinit,'\0') - 1;
  *length = (endinit - startinit + 1)/sizeof(char);
  return *length;
}

/* Used only by extern procedures (outside of this file).  Internal
   procedures have their own specialized creators. */
T
Sequence_genomic_new (char *contents, int length) {
  T new = (T) MALLOC(sizeof(*new));

  new->acc = (char *) NULL;
  new->restofheader = (char *) NULL;
  new->contents = contents;
  new->trimstart = 0;
  new->trimend = new->fulllength = length;
  new->free_contents_p = false;	/* Called only by Genome_get_segment, which provides
				   its own buffer */
  return new;
}

static char *
make_complement (char *sequence, unsigned int length) {
  char *complement;
  char complCode[128] = COMPLEMENT;
  int i, j;

  complement = (char *) CALLOC(length+1,sizeof(char));
  for (i = length-1, j = 0; i >= 0; i--, j++) {
    complement[j] = complCode[(int) sequence[i]];
  }
  complement[length] = '\0';
  return complement;
}

static void
make_complement_buffered (char *complement, char *sequence, unsigned int length) {
  char complCode[128] = COMPLEMENT;
  int i, j;

  /* complement = (char *) CALLOC(length+1,sizeof(char)); */
  for (i = length-1, j = 0; i >= 0; i--, j++) {
    complement[j] = complCode[(int) sequence[i]];
  }
  complement[length] = '\0';
  return;
}

T
Sequence_subsequence (T this, int start, int end) {
  T new;

  if (start < this->trimstart) {
    start = this->trimstart;
  }
  if (end > this->trimend) {
    end = this->trimend;
  }

  if (end <= start) {
    return NULL;
  } else {
    new = (T) MALLOC(sizeof(*new));

    new->acc = (char *) NULL;
    new->restofheader = (char *) NULL;
    new->contents = this->contents;
    new->fulllength = this->fulllength;
    new->trimstart = start;
    new->trimend = end;
    new->free_contents_p = false;
    return new;
  }
}


T
Sequence_revcomp (T this) {
  T new = (T) MALLOC(sizeof(*new));

  new->acc = (char *) NULL;
  new->restofheader = (char *) NULL;
  new->contents = make_complement(this->contents,this->fulllength);
  new->fulllength = this->fulllength;
  new->trimstart = this->trimstart;
  new->trimend = this->trimend;
  new->free_contents_p = true;
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
Sequence_read (int *nextchar, FILE *input) {
  T new;
  int savelength, skiplength;

  if (feof(input)) {
    *nextchar = EOF;
    return NULL;
  }

  new = (T) MALLOC(sizeof(*new));

  if (Initc == '\0') {
    if ((Initc = input_init(input)) == EOF) {
      *nextchar = EOF;
      return NULL;
    }
  }
  if (Initc != '>') {
    blank_header(new);
  } else if (input_header(input,new) == NULL) {
    /* File ends after >.  Don't process. */
      *nextchar = EOF;
     return NULL;
  } 
  if (input_sequence(&(*nextchar),&savelength,&skiplength,input) == 0) {
    /* File ends during header.  Continue with a sequence of length 0. */
    /* fprintf(stderr,"File ends after header\n"); */
  }

  if (skiplength > 0) {
    fprintf(stderr,"Warning: sequence exceeds maximum length of %d.  Truncating middle.\n",MAXSEQLEN);
  }

  new->fulllength = savelength + skiplength;

  new->trimstart = 0;
  new->trimend = new->fulllength;

  new->contents = (char *) CALLOC(new->fulllength+1,sizeof(char));
  strncpy(new->contents,startinit,new->fulllength);
  new->free_contents_p = true;

  return new;
}

T
Sequence_read_unlimited (FILE *input) {
  T new;
  Intlist_T intlist = NULL;
  char *p;
  int length;
  bool eolnp;

  if (feof(input)) {
    return NULL;
  }

  new = (T) MALLOC(sizeof(*new));

  if (Initc == '\0') {
    if ((Initc = input_init(input)) == EOF) {
      return NULL;
    }
  }
  if (Initc != '>') {
    blank_header(new);
  } else if (input_header(input,new) == NULL) {
    /* File ends after >.  Don't process. */
     return NULL;
  } 
  /* Don't touch Sequence[0], because subsequent calls to
     Sequence_read depend on it being '\0'. */
  eolnp = true;
  while (fgets(&(Sequence[1]),MAXSEQLEN,input) != NULL &&
	 (eolnp == false || Sequence[1] != '>')) {
    for (p = &(Sequence[1]); *p != '\n' && *p != '\0'; p++) {
      if (!iscntrl((int) *p)) {
	intlist = Intlist_push(intlist,(int) *p);
      }
    }
    if (*p == '\n') {
      eolnp = true;
    } else {
      eolnp = false;
    }
  }
  intlist = Intlist_reverse(intlist);
  new->contents = Intlist_to_char_array(&length,intlist);
  Intlist_free(&intlist);

  if (length == 0) {
    return NULL;
  } else {
    new->fulllength = new->trimend = length;
    new->trimstart = 0;

    new->free_contents_p = true;

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

void
Sequence_print_header (T this, bool checksump) {

  printf(">%s (%d bp) %s",this->acc,this->fulllength,this->restofheader);
  if (checksump == true) {
    printf(" md5:");
    Sequence_print_digest(this);
  }
  printf("\n");
  return;
}

void
Sequence_print (T this, bool uppercasep, int wraplength, bool trimmedp) {
  int i = 0, pos, start, end;

  if (trimmedp == true) {
    start = this->trimstart;
    end = this->trimend;
  } else {
    start = 0;
    end = this->fulllength;
  }

  if (uppercasep == true) {
    for (pos = start; pos < end; pos++, i++) {
      printf("%c",toupper((int) (this->contents[i])));
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

