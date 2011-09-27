static char rcsid[] = "$Id: iit_dump.c,v 1.15 2009-01-21 21:22:33 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>		/* For getopt */
#endif
#include "bool.h"
#include "iit-read.h"
#include "list.h"
#include "getopt.h"

/************************************************************************
 *   Program options
 ************************************************************************/

static bool debugp = false;
static bool tagsp = false;
static bool countsp = false;
static bool annotationonlyp = false;

static struct option long_options[] = {
  /* Input options */
  {"debug", no_argument, 0, '9'}, /* debugp */
  {"tags", no_argument, 0, 'T'}, /* tagsp */
  {"counts", no_argument, 0, 'C'}, /* countsp */
  {"annotonly", no_argument, 0, 'A'},	/* annotationonlyp */

  /* Help options */
  {"version", no_argument, 0, 'V'}, /* print_program_version */
  {"help", no_argument, 0, '?'}, /* print_program_usage */
  {0, 0, 0, 0}
};

static void
print_program_version () {
  fprintf(stdout,"\n");
  fprintf(stdout,"iit_dump: debugging utility for Interval Index Trees\n");
  fprintf(stdout,"Part of GMAP package, version %s\n",PACKAGE_VERSION);
  fprintf(stdout,"Thomas D. Wu, Genentech, Inc.\n");
  fprintf(stdout,"Contact: twu@gene.com\n");
  fprintf(stdout,"\n");
  return;
}

static void
print_program_usage () {
  fprintf(stdout,"\
Usage: iit_dump [OPTIONS...] iitfile\n\
\n\
Options\n\
  -T, --tags              Show tags present in iit file\n\
  -C, --counts            Show counts for every boundary in iit file\n\
  -9, --debug             Provide debugging information\n\
  -A, --annotonly         Dump annotation lines only (no headers)\n\
\n\
  -V, --version           Show version\n\
  -?, --help              Show this help message\n\
\n\
The iit_dump program shows the entire contents of a given iit file.\n\
The default behavior is generate FASTA-type output, with both headers\n\
and annotations.  If only the annotations are desired, the -A flag\n\
may be used.  This flag may be useful for iit files created using the -G\n\
flag to iit_store, which stores the original gff3-formatted lines as\n\
the annotation.\n\
\n\
See also: iit_store, iit_get\n\
");

  return;
}

/*
static void
show_types (IIT_T iit) {
  List_T typelist, p;

  typelist = IIT_typelist(iit);
  for (p = typelist; p != NULL; p = List_next(p)) {
    printf("%s.\n",(char *) List_head(p));
  }
  return;
}
*/

int
main (int argc, char *argv[]) {
  char *iitfile;
  IIT_T iit;
  int type;
  
  int opt;
  extern int optind;
  extern char *optarg;
  int long_option_index = 0;

  while ((opt = getopt_long(argc,argv,"9TCA",long_options,&long_option_index)) != -1) {
    switch (opt) {
    case '9': debugp = true; break;
    case 'T': tagsp = true; break;
    case 'C': countsp = true; break;
    case 'A': annotationonlyp = true; break;

    case 'V': print_program_version(); exit(0);
    case '?': print_program_usage(); exit(0);
    default: exit(9);
    }
  }

  argc -= optind;
  argv += optind;

  iitfile = argv[0];
  if (debugp == true) {
    IIT_debug(iitfile);
    return 0;
  } else if ((iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
			     /*divstring*/NULL,/*add_iit_p*/true,/*labels_read_p*/true)) == NULL) {
    fprintf(stderr,"Unable to open or parse IIT file %s\n",iitfile);
    exit(9);
  }

  if (tagsp == true) {
    for (type = 1; type < IIT_ntypes(iit); type++) {
      printf("%s\n",IIT_typestring(iit,type));
    }
  } else if (countsp == true) {
    fprintf(stderr,"Flag -C not implemented\n");
#if 0
    IIT_dump_counts(iit,/*alphabetizep*/true);
#endif

  } else {
    IIT_dump(iit,annotationonlyp);
  }

  IIT_free(&iit);

  return 0;
}
