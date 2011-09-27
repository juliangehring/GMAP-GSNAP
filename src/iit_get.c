static char rcsid[] = "$Id: iit_get.c,v 1.28 2005/02/15 01:55:14 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>		/* For rindex */
#include <ctype.h>
#include "bool.h"
#include "mem.h"
#include "intlist.h"
#include "interval.h"
#include "iit-read.h"

#define BUFLEN 1024

static bool annotationonlyp = false;

static bool
isnumberp (char *string) {
  char *p = string;

  while (*p != '\0') {
    if (!isdigit(*p)) {
      return false;
    }
    p++;
  }
  return true;
}

int 
main (int argc, char *argv[]) {
  char *iitfile, *ptr;
  char Buffer[BUFLEN], nocomment[BUFLEN], label[BUFLEN], typestring[BUFLEN], *annotation;
  unsigned int query1, query2;
  int typeint;
  int index, nargs;
  Intlist_T matches, p;
  IIT_T iit;
  Interval_T interval;
  bool skipp;
  
  int c;
  extern int optind;
  extern char *optarg;

  while ((c = getopt(argc,argv,"")) != -1) {
    switch (c) {
    case 'A': annotationonlyp = true; break;
    }
  }

  argc -= optind;
  argv += optind;

  if (argc == 0) {
    fprintf(stderr,"Need to specify an iit file.  Type \"iit_get --help\" for help.\n");
    exit(9);
  }

  iitfile = argv[0];
  if ((iit = IIT_read(iitfile,NULL,true)) == NULL) {
    iitfile = (char *) CALLOC(strlen(argv[0])+strlen(".iit")+1,sizeof(char));
    sprintf(iitfile,"%s.iit",argv[0]);
    if ((iit = IIT_read(iitfile,NULL,true)) == NULL) {
      fprintf(stderr,"Can't open IIT file %s\n",argv[0]);
      exit(9);
    } else {
      FREE(iitfile);
    }
  }

  if (argc == 1) {
    /* stdin */
    while (fgets(Buffer,BUFLEN,stdin) != NULL) {
      if ((ptr = rindex(Buffer,'\n')) != NULL) {
	*ptr = '\0';
      }
      strcpy(nocomment,Buffer);
      if ((ptr = rindex(nocomment,'#')) != NULL) {
	*ptr = '\0';
      }

      skipp = false;
      if ((nargs = sscanf(nocomment,"%u %u %s",&query1,&query2,typestring)) == 3) {
	if ((typeint = IIT_typeint(iit,typestring)) < 0) {
	  fprintf(stderr,"For %s, no such type as %s.  Ignoring the type.\n",nocomment,typestring);
	  matches = IIT_get(iit,query1,query2);
	} else {
	  matches = IIT_get_typed(iit,query1,query2,typeint);
	}
      } else if (nargs == 2) {
	matches = IIT_get(iit,query1,query2);
      } else if (sscanf(nocomment,"%s",label) == 1) {
	matches = IIT_find(iit,label);
      } else {
	fprintf(stderr,"Can't parse line %s.  Ignoring.\n",nocomment);
	skipp = true;
      }
	
      if (skipp == false) {
	printf(">>%s\n",Buffer);
	for (p = matches; p != NULL; p = Intlist_next(p)) {
	  index = Intlist_head(p);
      
	  if (annotationonlyp == false) {
	    printf(">%s",IIT_label(iit,index));

	    interval = IIT_interval(iit,index);
	    printf(" %u %u",Interval_low(interval),Interval_high(interval));
	    if (Interval_type(interval) > 0) {
	      printf(" %s",IIT_typestring(iit,Interval_type(interval)));
	    }
	    printf("\n");
	  }
	  annotation = IIT_annotation(iit,index);
	  if (strlen(annotation) == 0) {
	  } else if (annotation[strlen(annotation)-1] == '\n') {
	    printf("%s",annotation);
	  } else {
	    printf("%s\n",annotation);
	  }
	}
      }
      Intlist_free(&matches);
    }

  } else {
    if (argc == 2) {
      /* Try as iitfile label */
      matches = IIT_find(iit,argv[1]);
      if (matches == NULL && isnumberp(argv[1])) {
	query1 = (unsigned int) atol(argv[1]);
	matches = IIT_get(iit,query1,query1);
      }
    } else if (argc == 3) {
      /* iitfile start end */
      if (!isnumberp(argv[1]) || !isnumberp(argv[2])) {
	fprintf(stderr,"Need to specify a numeric start and end\n");
	exit(9);
      } else {
	query1 = (unsigned int) atol(argv[1]);
	query2 = (unsigned int) atol(argv[2]);
	matches = IIT_get(iit,query1,query2);
      }
    } else {
      query1 = (unsigned int) atol(argv[1]);
      query2 = (unsigned int) atol(argv[2]);
      if ((typeint = IIT_typeint(iit,argv[3])) < 0) {
	fprintf(stderr,"For %s %s %s, no such type as %s.  Ignoring the type.\n",
		argv[1],argv[2],argv[3],argv[3]);
	matches = IIT_get(iit,query1,query2);
      } else {
	matches = IIT_get_typed(iit,query1,query2,typeint);
      }
    }

    for (p = matches; p != NULL; p = Intlist_next(p)) {
      index = Intlist_head(p);
      
      if (annotationonlyp == false) {
	printf(">%s",IIT_label(iit,index));
      
	interval = IIT_interval(iit,index);
	printf(" %u %u",Interval_low(interval),Interval_high(interval));
	if (Interval_type(interval) > 0) {
	  printf(" %s",IIT_typestring(iit,Interval_type(interval)));
	}
	printf("\n");
      }
      annotation = IIT_annotation(iit,index);
      if (strlen(annotation) == 0) {
      } else if (annotation[strlen(annotation)-1] == '\n') {
	printf("%s",annotation);
      } else {
	printf("%s\n",annotation);
      }
    }
    Intlist_free(&matches);
  }

  IIT_free_mmapped(&iit);

  return 0;
}
