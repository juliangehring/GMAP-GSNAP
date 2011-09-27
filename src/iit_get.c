static char rcsid[] = "$Id: iit_get.c,v 1.38 2006/11/02 02:19:10 twu Exp $";
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
#include "getopt.h"

#define BUFLEN 1024

/************************************************************************
 *   Program options
 ************************************************************************/

static bool annotationonlyp = false;
static int nflanking = 0;

static struct option long_options[] = {
  /* Input options */
  {"annotonly", no_argument, 0, 'A'},	/* annotationonlyp */

  /* Help options */
  {"version", no_argument, 0, 'V'}, /* print_program_version */
  {"help", no_argument, 0, '?'}, /* print_program_usage */
  {0, 0, 0, 0}
};

static void
print_program_version () {
  fprintf(stdout,"\n");
  fprintf(stdout,"iit_get: retrieval utility for Interval Index Trees\n");
  fprintf(stdout,"Part of GMAP package, version %s\n",PACKAGE_VERSION);
  fprintf(stdout,"Thomas D. Wu, Genentech, Inc.\n");
  fprintf(stdout,"Contact: twu@gene.com\n");
  fprintf(stdout,"\n");
  return;
}

static void
print_program_usage () {
  fprintf(stdout,"\
Usage: iit_get [OPTIONS...] iitfile start, or\n\
       iit_get [OPTIONS...] iitfile start end, or\n\
       iit_get [OPTIONS...] iitfile start end tags..., or\n\
       iit_get [OPTIONS...] iitfile label, or\n\
       iit_get [OPTIONS...] iitfile\n\
\n\
Options\n\
  -A, --annotonly         Show annotation lines only (no headers)\n\
  -u, --flanking=INT      Show flanking segments on left and right\n\
\n\
  -V, --version           Show version\n\
  -?, --help              Show this help message\n\
\n\
\n\
The iit_get program retrieves segments from an iit file that overlap a\n\
given coordinate or pair of coordinates.  Retrieval is done in\n\
logarithmic time (for more details on IIT files, see Wu and Watanabe,\n\
Bioinformatics 21:1859-1875, 2005).  The start coordinate should be\n\
less than or equal to the end coordinate.  If only a single start\n\
coordinate is provided, this is equivalent to providing the same\n\
number for the start and end coordinate.\n\
\n\
The given iit file may contain tags (which can be displayed by using\n\
the -T flag of the iit_dump program).  These tags may be used to\n\
filter the output of the iit file.  Multiple tags may be specified,\n\
which indicates a disjunctive query, such that iit_get returns entries\n\
that match any one of the given tags.\n\
\n\
In the last usage, the program will expect one or more sets of coordinates\n\
from stdin, one per line.\n\
\n\
See also: iit_store, iit_dump\n\
");
  return;
}


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

static void
print_interval (int index, IIT_T iit) {
  Interval_T interval;
  char *annotation;
  bool allocp;

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
  return;
}


int 
main (int argc, char *argv[]) {
  char *iitfile, *ptr;
  char Buffer[BUFLEN], nocomment[BUFLEN], label[BUFLEN], typestring[BUFLEN], *annotation;
  unsigned int query1, query2;
  int typeint, *types, c;
  int nargs, ntypes;
  int *matches, nmatches, i, *leftflanks, *rightflanks, nleftflanks = 0, nrightflanks = 0;
  IIT_T iit;
  bool skipp;
  
  int opt;
  extern int optind;
  extern char *optarg;
  int long_option_index = 0;

  while ((opt = getopt_long(argc,argv,"Au:",long_options,&long_option_index)) != -1) {
    switch (opt) {
    case 'A': annotationonlyp = true; break;
    case 'u': nflanking = atoi(optarg); break;

    case 'V': print_program_version(); exit(0);
    case '?': print_program_usage(); exit(0);
    default: exit(9);
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
#ifdef ROBUST
	  fprintf(stderr,"For %s, no such type as %s.  Ignoring the type.\n",nocomment,typestring);
	  matches = IIT_get(&nmatches,iit,query1,query2);
	  if (nflanking > 0) {
	    IIT_get_flanking(&leftflanks,&nleftflanks,&rightflanks,&nrightflanks,iit,query1,query2,nflanking);
	  }
#else
	  matches = NULL;
	  nmatches = 0;
#endif
	} else {
	  matches = IIT_get_typed(&nmatches,iit,query1,query2,typeint);
	  if (nflanking > 0) {
	    IIT_get_flanking_typed(&leftflanks,&nleftflanks,&rightflanks,&nrightflanks,iit,query1,query2,nflanking,typeint);
	  }
	}
      } else if (nargs == 2) {
	matches = IIT_get(&nmatches,iit,query1,query2);
	if (nflanking > 0) {
	  IIT_get_flanking(&leftflanks,&nleftflanks,&rightflanks,&nrightflanks,iit,query1,query2,nflanking);
	}
      } else if (sscanf(nocomment,"%s",label) == 1) {
	matches = IIT_find(&nmatches,iit,label);
      } else {
	fprintf(stderr,"Can't parse line %s.  Ignoring.\n",nocomment);
	skipp = true;
      }
	
      if (skipp == false) {
	fprintf(stdout,"# Query: %s\n",Buffer);
	for (i = 0; i < nmatches; i++) {
	  print_interval(matches[i],iit);
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
	query1 = (unsigned int) strtoul(argv[1],(char **) NULL,10);
	matches = IIT_get(&nmatches,iit,query1,query1);
	if (nflanking > 0) {
	  IIT_get_flanking(&leftflanks,&nleftflanks,&rightflanks,&nrightflanks,iit,query1,query1,nflanking);
	}
      }
    } else if (argc == 3) {
      /* iitfile start end */
      if (!isnumberp(argv[1]) || !isnumberp(argv[2])) {
	fprintf(stderr,"Need to specify a numeric start and end\n");
	exit(9);
      } else {
	query1 = (unsigned int) strtoul(argv[1],(char **) NULL,10);
	query2 = (unsigned int) strtoul(argv[2],(char **) NULL,10);
	matches = IIT_get(&nmatches,iit,query1,query2);
	if (nflanking > 0) {
	  IIT_get_flanking(&leftflanks,&nleftflanks,&rightflanks,&nrightflanks,iit,query1,query2,nflanking);
	}
      }
    } else if (argc == 4) {
      query1 = (unsigned int) strtoul(argv[1],(char **) NULL,10);
      query2 = (unsigned int) strtoul(argv[2],(char **) NULL,10);
      if ((typeint = IIT_typeint(iit,argv[3])) < 0) {
#ifdef ROBUST
	fprintf(stderr,"For %s %s %s, no such type as %s.  Ignoring the type.\n",
		argv[1],argv[2],argv[3],argv[3]);
	matches = IIT_get(&nmatches,iit,query1,query2);
	if (nflanking > 0) {
	  IIT_get_flanking(&leftflanks,&nleftflanks,&rightflanks,&nrightflanks,iit,query1,query2,nflanking);
	}
#else
	matches = NULL;
	nmatches = 0;
#endif
      } else {
	matches = IIT_get_typed(&nmatches,iit,query1,query2,typeint);
	if (nflanking > 0) {
	  IIT_get_flanking_typed(&leftflanks,&nleftflanks,&rightflanks,&nrightflanks,iit,query1,query2,nflanking,typeint);
	}
      }
    } else {
      query1 = (unsigned int) strtoul(argv[1],(char **) NULL,10);
      query2 = (unsigned int) strtoul(argv[2],(char **) NULL,10);
      types = (int *) CALLOC(argc-3,sizeof(int));
      for (c = 3, ntypes = 0; c < argc; c++) {
	if ((typeint = IIT_typeint(iit,argv[c])) < 0) {
	  fprintf(stderr,"No such type as %s.  Ignoring the type.\n",argv[c]);
	} else {
	  types[ntypes++] = typeint;
	}
      }
      if (ntypes == 0) {
#ifdef ROBUST
	matches = IIT_get(&nmatches,iit,query1,query2);
	if (nflanking > 0) {
	  IIT_get_flanking(&leftflanks,&nleftflanks,&rightflanks,&nrightflanks,iit,query1,query2,nflanking);
	}
#else
	matches = NULL;
	nmatches = 0;
#endif
      } else if (ntypes == 1) {
	matches = IIT_get_typed(&nmatches,iit,query1,query2,types[0]);
	if (nflanking > 0) {
	  IIT_get_flanking_typed(&leftflanks,&nleftflanks,&rightflanks,&nrightflanks,iit,query1,query2,nflanking,typeint);
	}
      } else {
	matches = IIT_get_multiple_typed(&nmatches,iit,query1,query2,types,ntypes);
	if (nflanking > 0) {
	  IIT_get_flanking_multiple_typed(&leftflanks,&nleftflanks,&rightflanks,&nrightflanks,iit,query1,query2,nflanking,types,ntypes);
	}
      }
    }

    if (nflanking > 0) {
      for (i = nleftflanks-1; i >= 0; i--) {
	print_interval(leftflanks[i],iit);
      }
      printf("====================\n");
      FREE(leftflanks);
    }

    for (i = 0; i < nmatches; i++) {
      print_interval(matches[i],iit);
    }

    if (nflanking > 0) {
      printf("====================\n");
      for (i = 0; i < nrightflanks; i++) {
	print_interval(rightflanks[i],iit);
      }
      FREE(rightflanks);
    }

    FREE(matches);
  }

  IIT_free(&iit);

  return 0;
}
