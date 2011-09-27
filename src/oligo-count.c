static char rcsid[] = "$Id: oligo-count.c,v 1.10 2005/07/08 07:58:33 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif
#include "mem.h"
#include "indexdb.h"
#include "sequence.h"
#include "oligo.h"
#include "block.h"
#include "reader.h"


#define DEFAULT_DATADIR "/usr/snap/data/genomes"

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/*               8765432187654321 */
#define LEFT_A 0x0000000000000000
#define LEFT_C 0x4000000000000000
#define LEFT_G 0x8000000000000000
#define LEFT_T 0xC000000000000000


/* Adapted from stage1.c */
static int
count_matches (int querypos0, Genomicpos_T *positions0, int npositions0, 
	       int querypos1, Genomicpos_T *positions1, int npositions1,
	       int expecteddist, Block_T block, bool forwardp) {
  int nentries = 0, i = 0, j = 0;
  Genomicpos_T position0, expected1, position1, lastmatch0;
  bool donep = false;

  if (npositions0 == 0) {
    return 0;
  } else {
    position0 = positions0[0];
    expected1 = position0 + expecteddist;
    if (npositions1 == 0) {
      return 0;
    } else {
      position1 = positions1[0];
    }
  }

  lastmatch0 = position0 - 1U;	/* Guaranteed not to be a match */
  while (!donep) {
    debug(printf("  %d:%u %d:%u\n",i,position0,j,position1));
    debug(
	  if (abs(position1-position0) < 100) {
	    printf("Close: %u %u\n",position0,position1);
	  }
	  );

    if (expected1 < position1) {
      /* Advance position0 */
      if (++i >= npositions0) {
	donep = true;
      } else {
	position0 = positions0[i];
	expected1 = position0 + expecteddist;
      }
    } else if (expected1 > position1) {
      /* Advance position1 */
      if (++j >= npositions1) {
	donep = true;
      } else {
	position1 = positions1[j];
      }

    } else {
      if (position0 != lastmatch0) {
	nentries++;
      }
      lastmatch0 = position0;

      /* Advance both */
      if (++i >= npositions0) {
	donep = true;
      } else {
	position0 = positions0[i];
	expected1 = position0 + expecteddist;
	if (++j >= npositions1) {
	  donep = true;
	} else {
	  position1 = positions1[j];
	}
      }
    }
  }

  return nentries;
}

/* Adapted from stage1.c */
static void
find_5prime_pairs (int *nplus, int *nminus,
		   Genomicpos_T **plus_positions, int *plus_npositions,
		   Genomicpos_T **minus_positions, int *minus_npositions,
		   Block_T block5, int querystart, int interval) {
  int prevpos;

  if ((prevpos = querystart - interval) > 0) {
    *nplus = count_matches(prevpos,plus_positions[prevpos],plus_npositions[prevpos],
			   querystart,plus_positions[querystart],plus_npositions[querystart],
			   interval,block5,true);
    *nminus = count_matches(querystart,minus_positions[querystart],minus_npositions[querystart],
			    prevpos,minus_positions[prevpos],minus_npositions[prevpos],
			    interval,block5,false);
  } else {
    *nplus = *nminus = 0;
  }

  return;
}


#ifdef __STRICT_ANSI__
int getopt (int argc, char *const argv[], const char *optstring);
#endif

int
main (int argc, char *argv[]) {
  Indexdb_T indexdb;
  Sequence_T queryseq;
  Reader_T reader;
  Block_T block5;
  FILE *input;
  char *datadir = NULL, *releasestring = NULL, *dbroot = NULL;
  char *oligo;
  int nextchar;
  int oligosize = 24, interval, nplus, nminus;
  int c, seqlength, pos = 0;
  Genomicpos_T **plus_positions, **minus_positions;
  int *plus_npositions, *minus_npositions;
  int querystart, prevpos;
  Storedoligomer_T *oligos;

  extern int optind;
  extern char *optarg;

  while ((c = getopt(argc,argv,"d:R:s:")) != -1) {
    switch (c) {
    case 'd': dbroot = optarg; break;
    case 'R': releasestring = optarg; break;
    case 's': oligosize = atoi(optarg); break;
    }
  }
  argc -= (optind - 1);
  argv += (optind - 1);

  if (releasestring == NULL || !strcmp(releasestring,"0")) {
    datadir = (char *) calloc(strlen(DEFAULT_DATADIR)+strlen("/")+strlen(dbroot)+1,
			      sizeof(char));
    sprintf(datadir,"%s/%s",DEFAULT_DATADIR,dbroot);
  } else {
    datadir = (char *) calloc(strlen(DEFAULT_DATADIR)+strlen("/")+strlen(dbroot)+strlen("_R")+
			      strlen(releasestring)+1,sizeof(char));
    sprintf(datadir,"%s/%s_R%s",DEFAULT_DATADIR,dbroot,releasestring);
  }

  if (argc > 1) {
    if ((input = fopen(argv[1],"r")) == NULL) {
      fprintf(stderr,"Can't open file %s\n",argv[1]);
      exit(9);
    }
  } else {
    input = stdin;
  }

  indexdb = Indexdb_new_genome(datadir,dbroot,false);

  queryseq = Sequence_read(&nextchar,input);
  seqlength = Sequence_fulllength(queryseq);
  reader = Reader_new(Sequence_fullpointer(queryseq),0,seqlength);
  block5 = Block_new(FIVE,reader);
  oligos = (Storedoligomer_T *) CALLOC(seqlength,sizeof(Storedoligomer_T));
  interval = oligosize - INDEX1PART;

  plus_positions = (Genomicpos_T **) CALLOC(seqlength,sizeof(Genomicpos_T *));
  minus_positions = (Genomicpos_T **) CALLOC(seqlength,sizeof(Genomicpos_T *));
  plus_npositions = (int *) CALLOC(seqlength,sizeof(int));
  minus_npositions = (int *) CALLOC(seqlength,sizeof(int));

  while (Block_next(block5) == true) {
    querystart = Block_querypos(block5);
    prevpos = querystart - interval;
    oligos[querystart] = Block_forward(block5);
    Block_process_oligo(&(plus_positions[querystart]),&(plus_npositions[querystart]),
			&(minus_positions[querystart]),&(minus_npositions[querystart]),
			block5,indexdb);
    find_5prime_pairs(&nplus,&nminus,plus_positions,plus_npositions,
		      minus_positions,minus_npositions,
		      block5,querystart,interval);

    if (prevpos >= 0) {
      oligo = Oligo_nt(oligos[prevpos],oligos[querystart],oligosize);
      printf("%d\t%s\t%d\t%d\t%d\n",prevpos,oligo,nplus+nminus,nplus,nminus);
      FREE(oligo);
    }
  }

  FREE(plus_positions);
  FREE(minus_positions);
  Block_free(&block5);
  Reader_free(&reader);
  Sequence_free(&queryseq);
  Indexdb_free(&indexdb);

  return 0;
}

