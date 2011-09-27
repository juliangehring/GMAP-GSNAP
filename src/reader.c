static char rcsid[] = "$Id: reader.c,v 1.13 2005/02/15 01:58:50 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "reader.h"
#include "mem.h"

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


#define T Reader_T
struct T {
  char *startinit;
  char *startptr;
  char *endptr;
};

int
Reader_startpos (T this) {
  return (this->startptr - this->startinit);
}

int
Reader_endpos (T this) {
  return (this->endptr - this->startinit);
}


void
Reader_reset_ends (T this, int querystart, int queryend) {
  char *sequence;

  sequence = this->startinit;
  this->startptr = &(sequence[querystart]);
  this->endptr = &(sequence[queryend]);
  return;
}



T
Reader_new (char *sequence, int querystart, int queryend) {
  T new = (T) MALLOC(sizeof(*new));

  new->startinit = sequence;
  new->startptr = &(sequence[querystart]);
  new->endptr = &(sequence[queryend]);
  return new;
}

void
Reader_free (T *old) {
  if (*old) {
    FREE(*old);
  }
  return;
}

char
Reader_getc (T this, cDNAEnd_T cdnaend) {
  debug(fprintf(stderr,"Read_getc has startptr %d and endptr %d\n",
		this->startptr-this->startinit,this->endptr-this->startinit));
  if (this->startptr >= this->endptr) {
    return '\0';
  } else if (cdnaend == FIVE) {
    return *(this->startptr++);
  } else { 
    return *(this->endptr--);
  }
}


/* For testing */
/*
static void
process_input (FILE *input) {
  bool ssfile;
  int queryno = 0, seqlength, skiplength;
  char *query, initc;

  if ((initc = Reader_input_init(input)) == '>') {
    ssfile = false;
    query = Reader_input_header(input);
  } else {
    ssfile = true;
    query = (char *) CALLOC(strlen("NO_HEADER")+1,sizeof(char));
    strcpy(query,"NO_HEADER");
  }
  fprintf(stderr,"Read sequence %s\n",query);
  Reader_input_sequence(&seqlength,&skiplength,input,initc);
  print_contents(&(Sequence[0]),SEQUENCELEN);
  fprintf(stderr,"\nSeqlength = %d\n",seqlength);
  FREE(query);
  queryno++;

  while ((query = Reader_input_header(input)) != NULL) {
    fprintf(stderr,"Read sequence %s\n",query);
    Reader_input_sequence(&seqlength,&skiplength,input,initc);
    print_contents(&(Sequence[0]),SEQUENCELEN);
    fprintf(stderr,"\nSeqlength = %d\n",seqlength);
    FREE(query);
    queryno++;
  }

  return;
}

int
main (int argc, char *argv[]) {
  process_input(stdin);
  return 0;
}
*/
