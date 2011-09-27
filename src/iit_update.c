static char rcsid[] = "$Id: iit_update.c,v 1.1 2007/09/18 20:55:56 twu Exp $";
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



/************************************************************************
 *   Program options
 ************************************************************************/

static char *outputfile = NULL;
static char iit_version = 0;	/* Means latest version */

static struct option long_options[] = {
  /* Input options */
  {"output", required_argument, 0, 'o'}, /* outputfile */
  {"iitversion", required_argument, 0, 'v'}, /* iit_version */

  /* Help options */
  {"version", no_argument, 0, 'V'}, /* print_program_version */
  {"help", no_argument, 0, '?'}, /* print_program_usage */
  {0, 0, 0, 0}
};


static void
print_program_version () {
  fprintf(stdout,"\n");
  fprintf(stdout,"iit_update: indexing utility for Interval Index Trees\n");
  fprintf(stdout,"Part of GMAP package, version %s\n",PACKAGE_VERSION);
  fprintf(stdout,"Thomas D. Wu, Genentech, Inc.\n");
  fprintf(stdout,"Contact: twu@gene.com\n");
  fprintf(stdout,"\n");
  return;
}


static void
print_program_usage () {
  fprintf(stdout,"\
Usage: iit_update [OPTIONS...] -o output_iit input_iit\n\
where\n\
   output_iit is the desired filename for the new iit file\n\
       (.iit will be added as a suffix if necessary), and\n\
   input_iit is the input iit file\n\
\n\
If no output_iit file is specified, the existing input_iit file is overwritten.\n\
\n\
Options\n\
  -o, --output=STRING       Name of output iit file (required)\n\
  -v, --iitversion=STRING   Desired iit version for output iit\n\
                            (default = 0, which means latest version)\n\
\n\
  -V, --version             Show version\n\
  -?, --help                Show this help message\n\
\n\
");
  return;
}


#ifdef __STRICT_ANSI__
int getopt (int argc, char *const argv[], const char *optstring);
#endif

int 
main (int argc, char *argv[]) {
  char *iitfile1, *iitfile2, *p;
  IIT_T iit;
  
  int opt;
  extern int optind;
  extern char *optarg;
  int long_option_index = 0;

  while ((opt = getopt_long(argc,argv,"o:v:",long_options,&long_option_index)) != -1) {
    switch (opt) {
    case 'o': outputfile = optarg; break;
    case 'v': iit_version = atoi(optarg); break;

    case 'V': print_program_version(); exit(0);
    case '?': print_program_usage(); exit(0);
    default: exit(9);
    }
  }

  argc -= optind;
  argv += optind;

  if (argc == 0) {
    fprintf(stderr,"Need to specify an iit file.  Type \"iit_update --help\" for help.\n");
    exit(9);
  } else if (outputfile == NULL) {
    fprintf(stderr,"Need to specify an output iit file with the -o flag.  Type \"iit_update --help\" for help.\n");
    exit(9);
  }

  iitfile1 = argv[0];
  if ((iit = IIT_read(iitfile1,NULL,true)) == NULL) {
    iitfile1 = (char *) CALLOC(strlen(argv[0])+strlen(".iit")+1,sizeof(char));
    sprintf(iitfile1,"%s.iit",argv[0]);
    if ((iit = IIT_read(iitfile1,NULL,true)) == NULL) {
      fprintf(stderr,"Can't open file %s or file is not an IIT file\n",argv[0]);
      exit(9);
    }
  }

  if (strlen(outputfile) < 4) {
    iitfile2 = (char *) CALLOC(strlen(outputfile)+strlen(".iit")+1,sizeof(char));
    sprintf(iitfile2,"%s.iit",outputfile);
  } else {
    p = &(outputfile[strlen(outputfile)]);
    p -= 4;
    if (!strcmp(p,".iit")) {
      iitfile2 = (char *) CALLOC(strlen(outputfile)+1,sizeof(char));
      strcpy(iitfile2,outputfile);
    } else {
      iitfile2 = (char *) CALLOC(strlen(outputfile)+strlen(".iit")+1,sizeof(char));
      sprintf(iitfile2,"%s.iit",outputfile);
    }
  }

  if (!strcmp(iitfile1,iitfile2)) {
    fprintf(stderr,"Output file %s must have a different name from input file\n",iitfile2);
    exit(9);
  } else {
    IIT_output_direct(iitfile2,iit,iit_version);
    FREE(iitfile2);
    FREE(iitfile1);
  }

  return 0;
}

