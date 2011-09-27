static char rcsid[] = "$Id: pmapindex.c,v 1.17 2009/02/02 03:25:40 twu Exp $";
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
#include "indexdb_dibase.h"

#define BUFFERSIZE 8192

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


/* Program variables */
typedef enum{NONE, COMPRESS, OFFSETS, POSITIONS} Action_T;
static Action_T action = NONE;
static char *sourcedir = ".";
static char *destdir = ".";
static char *fileroot = NULL;
static int index1interval = 3;
static bool genome_lc_p = false;
static bool writefilep = false;
static bool mask_lowercase_p = false;


static void
write_offsets (FILE *input, char *fileroot, char interval_char) {
  char *offsetsfile, *iitfile;
  FILE *offsets_fp;
  IIT_T chromosome_iit;

  /* Reference strain */
  offsetsfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
				strlen(".")+strlen(DIBASE_FILESUFFIX)+/*for interval_char*/1+
				strlen(OFFSETS_FILESUFFIX)+1,sizeof(char));
  sprintf(offsetsfile,"%s/%s.%s%c%s",destdir,fileroot,DIBASE_FILESUFFIX,interval_char,OFFSETS_FILESUFFIX);
  if ((offsets_fp = FOPEN_WRITE_BINARY(offsetsfile)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",offsetsfile);
    exit(9);
  }
  
  iitfile = (char *) CALLOC(strlen(sourcedir)+strlen("/")+
			    strlen(fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
  sprintf(iitfile,"%s/%s.chromosome.iit",sourcedir,fileroot);
  chromosome_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,
			    /*divread*/READ_ALL,/*divstring*/NULL,/*add_iit_p*/false,
			    /*labels_read_p*/false);
  FREE(iitfile);

  Indexdb_write_offsets(offsets_fp,input,chromosome_iit,index1interval,
			genome_lc_p,fileroot,mask_lowercase_p);
  
  IIT_free(&chromosome_iit);
  fclose(offsets_fp);
  FREE(offsetsfile);

  return;
}

static void
write_positions (FILE *input, char *fileroot, char interval_char) {
  char *positionsfile, *offsetsfile, *iitfile;
  FILE *offsets_fp;
  IIT_T chromosome_iit;

  /* Reference strain */
  offsetsfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
				strlen(".")+strlen(DIBASE_FILESUFFIX)+/*for interval_char*/1+
				strlen(OFFSETS_FILESUFFIX)+1,sizeof(char));
  sprintf(offsetsfile,"%s/%s.%s%c%s",destdir,fileroot,DIBASE_FILESUFFIX,interval_char,OFFSETS_FILESUFFIX);
  if ((offsets_fp = FOPEN_READ_BINARY(offsetsfile)) == NULL) {
    fprintf(stderr,"Can't open file %s\n",offsetsfile);
    exit(9);
  }

  positionsfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
				  strlen(".")+strlen(DIBASE_FILESUFFIX)+/*for interval char*/1+
				  strlen(POSITIONS_FILESUFFIX)+1,sizeof(char));
  sprintf(positionsfile,"%s/%s.%s%c%s",destdir,fileroot,DIBASE_FILESUFFIX,interval_char,POSITIONS_FILESUFFIX);

  iitfile = (char *) CALLOC(strlen(sourcedir)+strlen("/")+
			    strlen(fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
  sprintf(iitfile,"%s/%s.chromosome.iit",sourcedir,fileroot);
  chromosome_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,
			    /*divread*/READ_ALL,/*divstring*/NULL,/*add_iit_p*/false,
			    /*labels_read_p*/false);
  FREE(iitfile);

  Indexdb_write_positions(positionsfile,offsets_fp,input,chromosome_iit,index1interval,
			  genome_lc_p,writefilep,fileroot,mask_lowercase_p);
  
  IIT_free(&chromosome_iit);
  fclose(offsets_fp);
  FREE(positionsfile);
  FREE(offsetsfile);

  return;
}

int 
main (int argc, char *argv[]) {
  char *genomefile, *outfile;
  FILE *input, *output;
  int interval_char;

  int c;
  extern int optind;
  extern char *optarg;

  while ((c = getopt(argc,argv,"F:D:d:COPs:Wq:m")) != -1) {
    switch (c) {
    case 'F': sourcedir = optarg; break;
    case 'D': destdir = optarg; break;
    case 'd': fileroot = optarg; break;
    case 'C': action = COMPRESS; break;
    case 'O': action = OFFSETS; break;
    case 'P': action = POSITIONS; break;
    case 'W': writefilep = true; break;
    case 'q': 
      index1interval = atoi(optarg);
      break;
    case 'm': mask_lowercase_p = true; break;
    }
  }
  argc -= (optind - 1);
  argv += (optind - 1);

  if (index1interval == 6) {
    interval_char = '6';
  } else if (index1interval == 3) {
    interval_char = '3';
  } else if (index1interval == 1) {
    interval_char = '1';
  } else {
    fprintf(stderr,"Selected indexing interval %d is not allowed.  Only values allowed are 6, 3, or 1\n",index1interval);
    exit(9);
  }

  if (fileroot == NULL) {
    fprintf(stderr,"Missing name of genome database.  Must specify with -d flag.\n");
    exit(9);
  }

  if (action == COMPRESS) {
    /* Usage: cat <genomefile> | dibaseindex -C > <genomecompfile>, or
              dibaseindex -C <genomefile> > <genomecompfile> */
    fprintf(stderr,"No need to compress color genome\n");
    exit(0);

#if 0
    if (argc > 1) {
      input = FOPEN_READ_BINARY(argv[1]);
      Indexdb_write_genomedibase(stdout,input=FOPEN_READ_BINARY(genomefile));
      fclose(input);
    } else {
      Indexdb_write_genomedibase(stdout,stdin);
    }
#endif

  } else if (action == OFFSETS) {
    /* Usage: cat <genomefile> | dibaseindex [-F <sourcedir>] [-D <destdir>] -d <dbname> -O
       Creates <destdir>/<dbname>.pfxoffsets and .prxoffsets */

    write_offsets(stdin,fileroot,interval_char);

  } else if (action == POSITIONS) {
    /* Usage: cat <genomefile> | dibaseindex [-F <sourcedir>] [-D <destdir>] -d <dbname> -P
       Requires <sourcedir>/<dbname>.pfxoffsets and .prxoffsets.
       Creates <destdir>/<dbname>.pfxpositions and .prxpositions */

    write_positions(stdin,fileroot,interval_char);

  } else {
    /* Perform all procedures.  Assume that current directory has the input file.  */

    genomefile = (char *) CALLOC(strlen(sourcedir)+strlen("/")+strlen(fileroot)+strlen(".genome")+1,sizeof(char));
    sprintf(genomefile,"%s/%s.genome",sourcedir,fileroot);
    if ((input = FOPEN_READ_TEXT(genomefile)) != NULL) {
      genome_lc_p = true;
    } else {
      FREE(genomefile);
      genomefile = (char *) CALLOC(strlen(sourcedir)+strlen("/")+strlen(fileroot)+strlen(".genomecomp")+1,sizeof(char));
      sprintf(genomefile,"%s/%s.genomecomp",sourcedir,fileroot);
      if ((input = FOPEN_READ_BINARY(genomefile)) != NULL) {
	genome_lc_p = false;
      } else {
	fprintf(stderr,"Cannot find file %s.genomecomp or %s.genome in source directory %s.\n",fileroot,fileroot,sourcedir);
	fprintf(stderr,"You may need to connect first to the directory that contains these files,\n");
	fprintf(stderr,"  or specify the source directory with the -F flag.\n");
	exit(9);
      }
    }

    fprintf(stderr,"Writing offsets file\n");
    write_offsets(input,fileroot,interval_char);
    fclose(input);

    fprintf(stderr,"Writing positions file\n");
    write_positions(input=FOPEN_READ_BINARY(genomefile),fileroot,interval_char);
    fclose(input);

#if 0
    fprintf(stderr,"Writing genomedibase file from %s\n",genomefile);
    outfile = (char *) CALLOC(strlen(sourcedir)+strlen("/")+strlen(fileroot)+strlen(".genomedibase")+1,sizeof(char));
    sprintf(outfile,"%s/%s.genomedibase",sourcedir,fileroot);
    output = FOPEN_WRITE_BINARY(outfile);
    Indexdb_write_genomedibase(output,input=FOPEN_READ_BINARY(genomefile));
    fclose(input);
    fclose(output);
    FREE(outfile);
#endif

    FREE(genomefile);

  }

  return 0;
}

