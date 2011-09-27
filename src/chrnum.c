static char rcsid[] = "$Id: chrnum.c,v 1.20 2005/07/25 17:58:09 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "chrnum.h"
#include <stdio.h>
#include <string.h>
#include <ctype.h>		/* for isalpha */
#include "mem.h"
#include "interval.h"

char *
Chrnum_to_string (Chrnum_T chrnum, IIT_T chromosome_iit) {
  char *string, *label;

  label = IIT_label(chromosome_iit,chrnum);
  string = (char *) CALLOC(strlen(label)+1,sizeof(char));
  strcpy(string,label);

  return string;
}

char *
Chrnum_to_string_signed (Chrnum_T chrnum, IIT_T chromosome_iit, bool watsonp) {
  char *string, *label;

  label = IIT_label(chromosome_iit,chrnum);
  string = (char *) CALLOC(strlen(label)+2,sizeof(char));
  if (watsonp) {
    string[0] = '+';
  } else {
    string[0] = '-';
  }
  strcpy(&(string[1]),label);
  return string;
}

unsigned int
Chrnum_length (Chrnum_T chrnum, IIT_T chromosome_iit) {
  return IIT_length(chromosome_iit,chrnum);
}

/* Can use Chrom_string_from_position instead */
unsigned int
Chrnum_offset (Chrnum_T chrnum, IIT_T chromosome_iit) {
  Interval_T interval;

  interval = IIT_interval(chromosome_iit,chrnum);
  return Interval_low(interval);
}

