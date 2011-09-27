static char rcsid[] = "$Id: iit_get.c,v 1.33 2005/10/19 03:53:07 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>		/* For getopt */
#endif
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
    if (!isdigit((int) *p)) {
      return false;
    }
    p++;
  }
  return true;
}

#ifdef __STRICT_ANSI__
int getopt (int argc, char *const argv[], const char *optstring);
#endif

int 
main (int argc, char *argv[]) {
  char *iitfile, *ptr;
  char Buffer[BUFLEN], nocomment[BUFLEN], label[BUFLEN], typestring[BUFLEN], *annotation;
  unsigned int query1, query2;
  int typeint;
  int index, nargs;
  int *matches, nmatches, i;
  IIT_T iit;
  Interval_T interval;
  bool skipp, allocp;
  
  int c;
  extern int optind;
  extern char *optarg;

  while ((c = getopt(argc,argv,"A")) != -1) {
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
	  matches = IIT_get(&nmatches,iit,query1,query2);
	} else {
	  matches = IIT_get_typed(&nmatches,iit,query1,query2,typeint);
	}
      } else if (nargs == 2) {
	matches = IIT_get(&nmatches,iit,query1,query2);
      } else if (sscanf(nocomment,"%s",label) == 1) {
	matches = IIT_find(&nmatches,iit,label);
      } else {
	fprintf(stderr,"Can't parse line %s.  Ignoring.\n",nocomment);
	skipp = true;
      }
	
      if (skipp == false) {
	fprintf(stdout,"# Query: %s\n",Buffer);
	for (i = 0; i < nmatches; i++) {
	  index = matches[i];
      
	  if (annotationonlyp == false) {
	    fprintf(stdout,">%s",IIT_label(iit,index));

	    interval = IIT_interval(iit,index);
	    fprintf(stdout," %u %u",Interval_low(interval),Interval_high(interval));
	    if (Interval_type(interval) > 0) {
	      fprintf(stdout," %s",IIT_typestring(iit,Interval_type(interval)));
	    }
	    fprintf(stdout,"\n");
	  }
	  annotation = IIT_annotation(iit,index,&allocp);
	  if (strlen(annotation) == 0) {
	  } else if (annotation[strlen(annotation)-1] == '\n') {
	    fprintf(stdout,"%s",annotation);
	  } else {
	    fprintf(stdout,"%s\n",annotation);
	  }
	  if (allocp == true) {
	    FREE(annotation);
	  }
	}
      }
      FREE(matches);
      fprintf(stdout,"# End\n");
      fflush(stdout);
    }

  } else {
    if (argc == 2) {
      /* Try as iitfile label */
      matches = IIT_find(&nmatches,iit,argv[1]);
      if (matches == NULL && isnumberp(argv[1])) {
	query1 = (unsigned int) atol(argv[1]);
	matches = IIT_get(&nmatches,iit,query1,query1);
      }
    } else if (argc == 3) {
      /* iitfile start end */
      if (!isnumberp(argv[1]) || !isnumberp(argv[2])) {
	fprintf(stderr,"Need to specify a numeric start and end\n");
	exit(9);
      } else {
	query1 = (unsigned int) atol(argv[1]);
	query2 = (unsigned int) atol(argv[2]);
	matches = IIT_get(&nmatches,iit,query1,query2);
      }
    } else {
      query1 = (unsigned int) atol(argv[1]);
      query2 = (unsigned int) atol(argv[2]);
      if ((typeint = IIT_typeint(iit,argv[3])) < 0) {
	fprintf(stderr,"For %s %s %s, no such type as %s.  Ignoring the type.\n",
		argv[1],argv[2],argv[3],argv[3]);
	matches = IIT_get(&nmatches,iit,query1,query2);
      } else {
	matches = IIT_get_typed(&nmatches,iit,query1,query2,typeint);
      }
    }

    for (i = 0; i < nmatches; i++) {
      index = matches[i];
      
      if (annotationonlyp == false) {
	printf(">%s",IIT_label(iit,index));
      
	interval = IIT_interval(iit,index);
	printf(" %u %u",Interval_low(interval),Interval_high(interval));
	if (Interval_type(interval) > 0) {
	  printf(" %s",IIT_typestring(iit,Interval_type(interval)));
	}
	printf("\n");
      }
      annotation = IIT_annotation(iit,index,&allocp);
      if (strlen(annotation) == 0) {
      } else if (annotation[strlen(annotation)-1] == '\n') {
	printf("%s",annotation);
      } else {
	printf("%s\n",annotation);
      }
      if (allocp == true) {
	FREE(annotation);
      }
    }
    FREE(matches);
  }

  IIT_free(&iit);

  return 0;
}
