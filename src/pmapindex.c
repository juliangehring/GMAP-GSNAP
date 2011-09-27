static char rcsid[] = "$Id: pmapindex.c,v 1.6 2005/11/23 17:29:34 twu Exp $";
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
#include <ctype.h>
#include <string.h>
#include <strings.h>		/* For rindex */
#include "bool.h"
#include "assert.h"
#include "mem.h"
#include "fopen.h"

#include "genomicpos.h"
#include "iit-write.h"
#include "iit-read.h"
#include "indexdb.h"

#define BUFFERSIZE 8192

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


/* Program variables */
typedef enum{NONE,OFFSETS,POSITIONS} Action_T;
static Action_T action = NONE;
static char *sourcedir = ".";
static char *destdir = ".";
static char *fileroot = NULL;
static int index1interval = INDEX1INTERVAL;
static bool uncompressedp = false;
static bool writefilep = false;
static bool watsonp = true;


int 
main (int argc, char *argv[]) {
  IIT_T altstrain_iit;
  int c, i;
  unsigned int genomelength;
  char *chromosomefile, *iitfile, *iittypefile, *offsetsfile, *positionsfile,
    *typestring, *copy;
  FILE *offsets_fp, *fp;

  extern int optind;
  extern char *optarg;

  while ((c = getopt(argc,argv,"F:D:d:OoPpRs:Wq:")) != -1) {
    switch (c) {
    case 'F': sourcedir = optarg; break;
    case 'D': destdir = optarg; break;
    case 'd': fileroot = optarg; break;
    case 'O': action = OFFSETS; uncompressedp = false; break;
    case 'o': action = OFFSETS; uncompressedp = true; break;
    case 'P': action = POSITIONS; uncompressedp = false; break;
    case 'p': action = POSITIONS; uncompressedp = true; break;
    case 'R': watsonp = false; break;
    case 'W': writefilep = true; break;
    case 'q': 
      index1interval = atoi(optarg);
      break;
    }
  }
  argc -= (optind - 1);
  argv += (optind - 1);

  if (action == OFFSETS) {
    /* Usage: cat <genomefile> | pmapindex [-F <sourcedir>] [-D <destdir>] -d <dbname> -O
       If alternate strains are present, requires <sourcedir>/<dbname>.altstrain.iit
       Creates <destdir>/<dbname>.pfxoffsets and .prxoffsets */

    /* Reference strain */
    if (watsonp == true) {
      offsetsfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
				    strlen(".pfxoffsets")+1,sizeof(char));
      sprintf(offsetsfile,"%s/%s.pfxoffsets",destdir,fileroot);
    } else {
      offsetsfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
				    strlen(".prxoffsets")+1,sizeof(char));
      sprintf(offsetsfile,"%s/%s.prxoffsets",destdir,fileroot);
    }
    if ((offsets_fp = FOPEN_WRITE_BINARY(offsetsfile)) == NULL) {
      fprintf(stderr,"Can't write to file %s\n",offsetsfile);
      exit(9);
    }

    iitfile = (char *) CALLOC(strlen(sourcedir)+strlen("/")+
			      strlen(fileroot)+strlen(".altstrain.iit")+1,sizeof(char));
    sprintf(iitfile,"%s/%s.altstrain.iit",sourcedir,fileroot);
    altstrain_iit = IIT_read(iitfile,NULL,true);

    Indexdb_write_offsets(offsets_fp,stdin,altstrain_iit,index1interval,watsonp,uncompressedp);

    IIT_free(&altstrain_iit);
    FREE(iitfile);
    fclose(offsets_fp);
    FREE(offsetsfile);


  } else if (action == POSITIONS) {
    /* Usage: cat <genomefile> | pmapindex [-F <sourcedir>] [-D <destdir>] -d <dbname> -P
       Requires <sourcedir>/<dbname>.pfxoffsets and .prxoffsets, and possibly <sourcedir>/<dbname>.altstrain.iit
       Creates <destdir>/<dbname>.pfxpositions and .prxpositions */

    /* Reference strain */
    if (watsonp == true) {
      offsetsfile = (char *) CALLOC(strlen(sourcedir)+strlen("/")+strlen(fileroot)+
				    strlen(".pfxoffsets")+1,sizeof(char));
      sprintf(offsetsfile,"%s/%s.pfxoffsets",sourcedir,fileroot);
    } else {
      offsetsfile = (char *) CALLOC(strlen(sourcedir)+strlen("/")+strlen(fileroot)+
				    strlen(".prxoffsets")+1,sizeof(char));
      sprintf(offsetsfile,"%s/%s.prxoffsets",sourcedir,fileroot);
    }
    if ((offsets_fp = FOPEN_READ_BINARY(offsetsfile)) == NULL) {
      fprintf(stderr,"Can't open file %s\n",offsetsfile);
      exit(9);
    }

    if (watsonp == true) {
      positionsfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
				      strlen(".pfxpositions")+1,sizeof(char));
      sprintf(positionsfile,"%s/%s.pfxpositions",destdir,fileroot);
    } else {
      positionsfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
				      strlen(".prxpositions")+1,sizeof(char));
      sprintf(positionsfile,"%s/%s.prxpositions",destdir,fileroot);
    }

    iitfile = (char *) CALLOC(strlen(sourcedir)+strlen("/")+
			      strlen(fileroot)+strlen(".altstrain.iit")+1,sizeof(char));
    sprintf(iitfile,"%s/%s.altstrain.iit",sourcedir,fileroot);
    altstrain_iit = IIT_read(iitfile,NULL,true);

    Indexdb_write_positions(positionsfile,offsets_fp,stdin,altstrain_iit,index1interval,
			    uncompressedp,watsonp,writefilep);

    IIT_free(&altstrain_iit);
    FREE(iitfile);
    fclose(offsets_fp);
    FREE(positionsfile);
    FREE(offsetsfile);
  }

  return 0;
}

