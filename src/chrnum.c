static char rcsid[] = "$Id: chrnum.c,v 1.18 2005/02/07 23:56:55 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "chrnum.h"
#include <stdio.h>
#include <string.h>
#include <ctype.h>		/* for isalpha */
#include "mem.h"

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

