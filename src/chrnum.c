static char rcsid[] = "$Id: chrnum.c 40271 2011-05-28 02:29:18Z twu $";
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
  bool allocp;

  label = IIT_label(chromosome_iit,chrnum,&allocp);
#if 0
  if (strip_spaces_p) {
    while (*label != '\0' && isspace(*label)) {
      label++;
    }
  }
#endif
  if (allocp == true) {
    return label;
  } else {
    string = (char *) CALLOC(strlen(label)+1,sizeof(char));
    strcpy(string,label);
    return string;
  }
}

char *
Chrnum_to_string_signed (Chrnum_T chrnum, IIT_T chromosome_iit, bool watsonp) {
  char *string, *label;
  bool allocp;

  label = IIT_label(chromosome_iit,chrnum,&allocp);
  string = (char *) CALLOC(strlen(label)+2,sizeof(char));
  if (watsonp) {
    string[0] = '+';
  } else {
    string[0] = '-';
  }
  strcpy(&(string[1]),label);

  if (allocp == true) {
    FREE(label);
  }
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

void
Chrnum_print_position (Genomicpos_T position, IIT_T chromosome_iit) {
  Chrnum_T chrnum;
  Genomicpos_T chrpos;
  int index;
  char *label;
  bool allocp;

  if (chromosome_iit == NULL) {
    chrnum = 0;
    chrpos = position;
  } else {
    index = IIT_get_one(chromosome_iit,/*divstring*/NULL,position,position);
    chrpos = position - Interval_low(IIT_interval(chromosome_iit,index));
    chrnum = index;
  }
  label = IIT_label(chromosome_iit,chrnum,&allocp);
  printf("#%d (chr%s):%u ",chrnum,label,chrpos);
  if (allocp == true) {
    FREE(label);
  }
  return;
}




