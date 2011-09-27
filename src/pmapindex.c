static char rcsid[] = "$Id: pmapindex.c,v 1.9 2006/11/30 17:17:07 twu Exp $";
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


static void
write_offsets (bool watsonp, FILE *input, char *fileroot) {
  char *offsetsfile, *iitfile;
  FILE *offsets_fp;
  IIT_T altstrain_iit;

  /* Reference strain */
  if (watsonp == true) {
    offsetsfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
				  strlen(".")+strlen(PFXOFFSETS)+1,sizeof(char));
    sprintf(offsetsfile,"%s/%s.%s",destdir,fileroot,PFXOFFSETS);
  } else {
    offsetsfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
				  strlen(".")+strlen(PRXOFFSETS)+1,sizeof(char));
    sprintf(offsetsfile,"%s/%s.%s",destdir,fileroot,PRXOFFSETS);
  }
  if ((offsets_fp = FOPEN_WRITE_BINARY(offsetsfile)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",offsetsfile);
    exit(9);
  }
  
  iitfile = (char *) CALLOC(strlen(sourcedir)+strlen("/")+
			    strlen(fileroot)+strlen(".altstrain.iit")+1,sizeof(char));
  sprintf(iitfile,"%s/%s.altstrain.iit",sourcedir,fileroot);
  altstrain_iit = IIT_read(iitfile,NULL,true);
  
  Indexdb_write_offsets(offsets_fp,input,altstrain_iit,index1interval,watsonp,uncompressedp,fileroot);
  
  IIT_free(&altstrain_iit);
  FREE(iitfile);
  fclose(offsets_fp);
  FREE(offsetsfile);

  return;
}

static void
write_positions (bool watsonp, FILE *input, char *fileroot) {
  char *positionsfile, *offsetsfile, *iitfile;
  FILE *offsets_fp;
  IIT_T altstrain_iit;

  /* Reference strain */
  if (watsonp == true) {
    offsetsfile = (char *) CALLOC(strlen(sourcedir)+strlen("/")+strlen(fileroot)+
				  strlen(".")+strlen(PFXOFFSETS)+1,sizeof(char));
    sprintf(offsetsfile,"%s/%s.%s",sourcedir,fileroot,PFXOFFSETS);
  } else {
    offsetsfile = (char *) CALLOC(strlen(sourcedir)+strlen("/")+strlen(fileroot)+
				  strlen(".")+strlen(PRXOFFSETS)+1,sizeof(char));
    sprintf(offsetsfile,"%s/%s.%s",sourcedir,fileroot,PRXOFFSETS);
  }
  if ((offsets_fp = FOPEN_READ_BINARY(offsetsfile)) == NULL) {
    fprintf(stderr,"Can't open file %s\n",offsetsfile);
    exit(9);
  }

  if (watsonp == true) {
    positionsfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
				    strlen(".")+strlen(PFXPOSITIONS)+1,sizeof(char));
    sprintf(positionsfile,"%s/%s.%s",destdir,fileroot,PFXPOSITIONS);
  } else {
    positionsfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
				    strlen(".")+strlen(PRXPOSITIONS)+1,sizeof(char));
    sprintf(positionsfile,"%s/%s.%s",destdir,fileroot,PRXPOSITIONS);
  }

  iitfile = (char *) CALLOC(strlen(sourcedir)+strlen("/")+
			    strlen(fileroot)+strlen(".altstrain.iit")+1,sizeof(char));
  sprintf(iitfile,"%s/%s.altstrain.iit",sourcedir,fileroot);
  altstrain_iit = IIT_read(iitfile,NULL,true);
  
  Indexdb_write_positions(positionsfile,offsets_fp,input,altstrain_iit,index1interval,
			  uncompressedp,watsonp,writefilep,fileroot);
  
  IIT_free(&altstrain_iit);
  FREE(iitfile);
  fclose(offsets_fp);
  FREE(positionsfile);
  FREE(offsetsfile);

  return;
}

int 
main (int argc, char *argv[]) {
  char *genomefile;
  FILE *input;

  int c;
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

  if (fileroot == NULL) {
    fprintf(stderr,"Missing name of genome database.  Must specify with -d flag.\n");
    exit(9);
  }

  if (action == OFFSETS) {
    /* Usage: cat <genomefile> | pmapindex [-F <sourcedir>] [-D <destdir>] -d <dbname> -O
       If alternate strains are present, requires <sourcedir>/<dbname>.altstrain.iit
       Creates <destdir>/<dbname>.pfxoffsets and .prxoffsets */

    write_offsets(watsonp,stdin,fileroot);

  } else if (action == POSITIONS) {
    /* Usage: cat <genomefile> | pmapindex [-F <sourcedir>] [-D <destdir>] -d <dbname> -P
       Requires <sourcedir>/<dbname>.pfxoffsets and .prxoffsets, and possibly <sourcedir>/<dbname>.altstrain.iit
       Creates <destdir>/<dbname>.pfxpositions and .prxpositions */

    write_positions(watsonp,stdin,fileroot);

  } else {
    /* Perform all four procedures.  Assume that current directory has the input file.  */

    genomefile = (char *) CALLOC(strlen(sourcedir)+strlen("/")+strlen(fileroot)+strlen(".genomecomp")+1,sizeof(char));
    sprintf(genomefile,"%s/%s.genomecomp",sourcedir,fileroot);
    if ((input = FOPEN_READ_BINARY(genomefile)) != NULL) {
      uncompressedp = false;
    } else {
      FREE(genomefile);
      genomefile = (char *) CALLOC(strlen(sourcedir)+strlen("/")+strlen(fileroot)+strlen(".genome")+1,sizeof(char));
      sprintf(genomefile,"%s/%s.genome",sourcedir,fileroot);
      if ((input = FOPEN_READ_TEXT(genomefile)) != NULL) {
	uncompressedp = true;
      } else {
	fprintf(stderr,"Cannot find file %s.genomecomp or %s.genome in source directory %s.\n",fileroot,fileroot,sourcedir);
	fprintf(stderr,"You may need to connect first to the directory that contains these files,\n");
	fprintf(stderr,"  or specify the source directory with the -F flag.\n");
	exit(9);
      }
    }

    fprintf(stderr,"Writing forward offsets file\n");
    write_offsets(/*watsonp*/true,input,fileroot);
    fclose(input);

    fprintf(stderr,"Writing reverse offsets file\n");
    write_offsets(/*watsonp*/false,input=FOPEN_READ_BINARY(genomefile),fileroot);
    fclose(input);

    fprintf(stderr,"Writing forward positions file\n");
    write_positions(/*watsonp*/true,input=FOPEN_READ_BINARY(genomefile),fileroot);
    fclose(input);

    fprintf(stderr,"Writing reverse positions file\n");
    write_positions(/*watsonp*/false,input=FOPEN_READ_BINARY(genomefile),fileroot);
    fclose(input);
      
    FREE(genomefile);

  }

  return 0;
}

