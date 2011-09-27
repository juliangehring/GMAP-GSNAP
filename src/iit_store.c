static char rcsid[] = "$Id: iit_store.c,v 1.23 2005/11/23 17:29:33 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif
#include <string.h>		/* For strlen */
#include <strings.h>		/* For rindex */
#include <ctype.h>
#include "assert.h"
#include "mem.h"
#include "fopen.h"

#include "list.h"
#include "interval.h"
#include "tableint.h"
#include "iit-write.h"


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


/* Empties contents of lines */
static char *
concatenate_lines (List_T lines, int content_size) {
  char *string, *temp;
  List_T l;

  string = (char *) CALLOC(content_size+1,sizeof(char));
  for (l = lines; l; l = List_next(l)) {
    temp = (char *) List_head(l);
    strcat(string,temp);
    FREE(temp);
  }
  
  /* Keep last return
  if (string[content_size-1] == '\n') {
    string[content_size-1] = '\0';
  }
  */

  return string;
}


static int
string_compare (const void *x, const void *y) {
  char *a = (char *) x;
  char *b = (char *) y;

  return strcmp(a,b);
}

/* This is the X31 hash function */
static unsigned int
string_hash (const void *x) {
  unsigned int h = 0U;
  const char *p;
  
  for (p = x; *p != '\0'; p++) {
    h = (h << 5) - h + *p;
  }
  return h;
}


static List_T
scan_header (List_T typelist, Tableint_T typetable, 
	     char **label, unsigned int *start, unsigned int *end, int *type,
	     char *header) {
  char Buffer[1024], *typestring, *p, *ptr;

  /* Example: >A 1 10 FWD */
  if (sscanf(header,">%s %u %u\n",Buffer,&(*start),&(*end)) < 3) {
    fprintf(stderr,"Error parsing %s\n",header);
    exit(9);
  } else {
    *label = (char *) CALLOC(strlen(Buffer)+1,sizeof(char));
    strcpy(*label,Buffer);

    p = header;
    while (!isspace((int) *p)) { p++; } /* First word */
    while (isspace((int) *p)) { p++; } /* First space */
    while (!isspace((int) *p)) { p++; } /* Second word */
    while (isspace((int) *p)) { p++; } /* Second space */
    while (!isspace((int) *p)) { p++; } /* Third word */
    while (*p != '\0' && isspace((int) *p)) { p++; } /* Third space */
    
    if (*p == '\0') {
      *type = 0;		/* Empty type string */
    } else {
      if ((ptr = rindex(p,'\n')) != NULL) {
	while (isspace((int) *ptr)) { ptr--; } /* Erase empty space */
	ptr++;
	*ptr = '\0';
      }
      if ((*type = Tableint_get(typetable,(void *) p)) == 0) {
	/* Store types as 1-based */
	*type = Tableint_length(typetable) + 1;
	typestring = (char *) CALLOC(strlen(p)+1,sizeof(char));
	strcpy(typestring,p);
	Tableint_put(typetable,typestring,*type);
	typelist = List_push(typelist,typestring);
	/* debug(printf("Entering new type %s.\n",typestring)); */
      }
    }
  }
  return typelist;
}

#ifdef __STRICT_ANSI__
int getopt (int argc, char *const argv[], const char *optstring);
#endif

int 
main (int argc, char *argv[]) {
  char *inputfile = NULL, *outputfile = NULL, *iitfile, 
    *tempstring, *typestring, *label, *p;
  char Buffer[8192];
  List_T lines = NULL, l, intervallist = NULL, typelist = NULL, labellist = NULL, annotlist = NULL;
  FILE *fp;
  unsigned int start, end;
  int content_size;
  Interval_T interval;
  Tableint_T typetable;
  int type;

  int c;
  extern int optind;
  extern char *optarg;

  while ((c = getopt(argc,argv,"o:")) != -1) {
    switch (c) {
    case 'o': outputfile = optarg; break;
    }
  }
  argc -= (optind - 1);
  argv += (optind - 1);

  if (outputfile == NULL) {
    fprintf(stderr,"Need to specify an output file with the -o flag\n");
    exit(9);
  }

  if (argc < 2) {
    fp = stdin;
  } else {
    inputfile = argv[1];
    fp = FOPEN_READ_TEXT(inputfile);
    if (!fp) {
      fprintf(stderr,"Can't open file %s\n",inputfile);
      exit(9);
    }
  }

  typetable = Tableint_new(1000,string_compare,string_hash);
  /* The zeroth type is empty */
  typestring = (char *) CALLOC(1,sizeof(char));
  typestring[0] = '\0';
  typelist = List_push(typelist,typestring);

  fgets(Buffer,8192,fp);
  typelist = scan_header(typelist,typetable,&label,&start,&end,&type,Buffer);
  labellist = List_push(NULL,label);
  lines = NULL;
  content_size = 0;

  while (fgets(Buffer,8192,fp) != NULL) {
    if (Buffer[0] == '>') {

      intervallist = List_push(intervallist,(void *) Interval_new(start,end,type));
      lines = List_reverse(lines);
      annotlist = List_push(annotlist,concatenate_lines(lines,content_size));
      List_free(&lines);

      typelist = scan_header(typelist,typetable,&label,&start,&end,&type,Buffer);
      labellist = List_push(labellist,label);
      lines = NULL;
      content_size = 0;

    } else {
      p = Buffer;
      tempstring = (char *) CALLOC(strlen(p)+1,sizeof(char));
      strcpy(tempstring,p);
      lines = List_push(lines,(void *) tempstring);
      content_size += strlen(p);
    }
  }

  intervallist = List_push(intervallist,(void *) Interval_new(start,end,type));
  lines = List_reverse(lines);
  annotlist = List_push(annotlist,concatenate_lines(lines,content_size));
  List_free(&lines);
  
  intervallist = List_reverse(intervallist);
  typelist = List_reverse(typelist);
  labellist = List_reverse(labellist);
  annotlist = List_reverse(annotlist);

  if (inputfile != NULL) {
    fclose(fp);
  }

  /* Figure out name of iit file */
  if (strlen(outputfile) < 4) {
    iitfile = (char *) CALLOC(strlen(outputfile)+strlen(".iit")+1,sizeof(char));
    sprintf(iitfile,"%s.iit",outputfile);
  } else {
    p = &(outputfile[strlen(outputfile)]);
    p -= 4;
    if (!strcmp(p,".iit")) {
      iitfile = (char *) CALLOC(strlen(outputfile)+1,sizeof(char));
      strcpy(iitfile,outputfile);
    } else {
      iitfile = (char *) CALLOC(strlen(outputfile)+strlen(".iit")+1,sizeof(char));
      sprintf(iitfile,"%s.iit",outputfile);
    }
  }

  IIT_write(iitfile,intervallist,typelist,labellist,annotlist,NULL);
  FREE(iitfile);

  for (l = annotlist; l != NULL; l = List_next(l)) {
    tempstring = (char *) List_head(l);
    FREE(tempstring);
  }
  List_free(&annotlist);

  for (l = labellist; l != NULL; l = List_next(l)) {
    tempstring = (char *) List_head(l);
    FREE(tempstring);
  }
  List_free(&labellist);

  for (l = typelist; l != NULL; l = List_next(l)) {
    tempstring = (char *) List_head(l);
    FREE(tempstring);
  }
  List_free(&typelist);

  Tableint_free(&typetable);

  for (l = intervallist; l != NULL; l = List_next(l)) {
    interval = (Interval_T) List_head(l);
    Interval_free(&interval);
  }
  List_free(&intervallist);

  return 0;
}
