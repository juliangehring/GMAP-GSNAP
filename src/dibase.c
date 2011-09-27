static char rcsid[] = "$Id: dibase.c,v 1.1 2009/08/29 00:28:47 twu Exp $";

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include "dibase.h"
#include <stdio.h>
#include <stdlib.h>		/* For abort() */

#include "stage1hr.h"		/* for MAX_QUERYLENGTH */
#include "genome.h"

#define COLOR_MISS 1
#define NT_MISS 1

/* mismatches_limit and mismatches_substring */
#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* mismatches_left */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

/* mismatches_right */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif


static char GENOME_CHARS[4] = {'A','C','G','T'};
static char GENOME_FLAGS[4] = {'N','N','N','N'};

typedef enum {A, C, G, T} Nucleotide_T;

struct Node_T {
  int colorscore[4];
  int ntscore[4];
  int totalscore[4];
  Nucleotide_T prev[4];
};



int
Dibase_count_mismatches_limit (int *color_nmismatches, char *queryuc_ptr, int pos5, int pos3, UINT4 *genome_blocks,
			       Genomicpos_T startpos, Genomicpos_T endpos, int max_mismatches) {
  char gbuffer[MAX_QUERYLENGTH+1], g, d;
  struct Node_T nodes[MAX_QUERYLENGTH];
  int bestscore, score, colorscore, ntscore;
  int seqlength, i, j;
  Nucleotide_T curr;

  Genome_uncompress_mmap(gbuffer,genome_blocks,startpos,endpos,GENOME_CHARS,GENOME_FLAGS);

  seqlength = endpos - startpos;

  debug({
    printf("Dibase_count_mismatches_limit at %u..%u with max_mismatches %d\n",
	   startpos,endpos,max_mismatches);
    printf("gbuffer: ");
    for (i = 0; i < pos5; i++) {
      printf(" ");
    }
    printf("%.*s\n",seqlength,gbuffer);
    printf("query:   %s\n",queryuc_ptr);
  });

  g = gbuffer[0];

  nodes[0].colorscore[A] = 0;
  nodes[0].colorscore[C] = 0;
  nodes[0].colorscore[G] = 0;
  nodes[0].colorscore[T] = 0;

  nodes[0].totalscore[A] = nodes[0].ntscore[A] = (g == 'A' ? 0 : NT_MISS);
  nodes[0].totalscore[C] = nodes[0].ntscore[C] = (g == 'C' ? 0 : NT_MISS);
  nodes[0].totalscore[G] = nodes[0].ntscore[G] = (g == 'G' ? 0 : NT_MISS);
  nodes[0].totalscore[T] = nodes[0].ntscore[T] = (g == 'T' ? 0 : NT_MISS);

  ntscore = 0;
  for (i = 1, j = pos5; i < seqlength && ntscore <= max_mismatches * NT_MISS; i++, j++) {
    g = gbuffer[i];
    d = queryuc_ptr[j];

    debug(printf("q[%d]:%c g[%d]:%c",j,d,i,g));


    /* A */
    bestscore = nodes[i-1].totalscore[A] + (d == '0' ? 0 : COLOR_MISS);
    colorscore = nodes[i-1].colorscore[A] + (d == '0' ? 0 : COLOR_MISS);
    ntscore = nodes[i-1].ntscore[A];

    if ((score = nodes[i-1].totalscore[C] + (d == '1' ? 0 : COLOR_MISS)) < bestscore) {
      bestscore = score;
      colorscore = nodes[i-1].colorscore[C] + (d == '1' ? 0 : COLOR_MISS);
      ntscore = nodes[i-1].ntscore[C];
    }

    if ((score = nodes[i-1].totalscore[G] + (d == '2' ? 0 : COLOR_MISS)) < bestscore) {
      bestscore = score;
      colorscore = nodes[i-1].colorscore[G] + (d == '2' ? 0 : COLOR_MISS);
      ntscore = nodes[i-1].ntscore[G];
    }

    if ((score = nodes[i-1].totalscore[T] + (d == '3' ? 0 : COLOR_MISS)) < bestscore) {
      bestscore = score;
      colorscore = nodes[i-1].colorscore[T] + (d == '3' ? 0 : COLOR_MISS);
      ntscore = nodes[i-1].ntscore[T];
    }

    nodes[i].totalscore[A] = bestscore + (g == 'A' ? 0 : NT_MISS);
    nodes[i].colorscore[A] = colorscore;
    nodes[i].ntscore[A] = ntscore + (g == 'A' ? 0 : NT_MISS);

    debug(printf("\t| A: %d=%d+%d",
		 nodes[i].totalscore[A],nodes[i].colorscore[A],nodes[i].ntscore[A]));

    /* C */
    bestscore = nodes[i-1].totalscore[A] + (d == '1' ? 0 : COLOR_MISS);
    colorscore = nodes[i-1].colorscore[A] + (d == '1' ? 0 : COLOR_MISS);
    ntscore = nodes[i-1].ntscore[A];

    if ((score = nodes[i-1].totalscore[C] + (d == '0' ? 0 : COLOR_MISS)) < bestscore) {
      bestscore = score;
      colorscore = nodes[i-1].colorscore[C] + (d == '0' ? 0 : COLOR_MISS);
      ntscore = nodes[i-1].ntscore[C];
    }

    if ((score = nodes[i-1].totalscore[G] + (d == '3' ? 0 : COLOR_MISS)) < bestscore) {
      bestscore = score;
      colorscore = nodes[i-1].colorscore[G] + (d == '3' ? 0 : COLOR_MISS);
      ntscore = nodes[i-1].ntscore[G];
    }

    if ((score = nodes[i-1].totalscore[T] + (d == '2' ? 0 : COLOR_MISS)) < bestscore) {
      bestscore = score;
      colorscore = nodes[i-1].colorscore[T] + (d == '2' ? 0 : COLOR_MISS);
      ntscore = nodes[i-1].ntscore[T];
    }

    nodes[i].totalscore[C] = bestscore + (g == 'C' ? 0 : NT_MISS);
    nodes[i].colorscore[C] = colorscore;
    nodes[i].ntscore[C] = ntscore + (g == 'C' ? 0 : NT_MISS);

    debug(printf("\t| C: %d=%d+%d",
		 nodes[i].totalscore[C],nodes[i].colorscore[C],nodes[i].ntscore[C]));


    /* G */
    bestscore = nodes[i-1].totalscore[A] + (d == '2' ? 0 : COLOR_MISS);
    colorscore = nodes[i-1].colorscore[A] + (d == '2' ? 0 : COLOR_MISS);
    ntscore = nodes[i-1].ntscore[A];

    if ((score = nodes[i-1].totalscore[C] + (d == '3' ? 0 : COLOR_MISS)) < bestscore) {
      bestscore = score;
      colorscore = nodes[i-1].colorscore[C] + (d == '3' ? 0 : COLOR_MISS);
      ntscore = nodes[i-1].ntscore[C];
    }

    if ((score = nodes[i-1].totalscore[G] + (d == '0' ? 0 : COLOR_MISS)) < bestscore) {
      bestscore = score;
      colorscore = nodes[i-1].colorscore[G] + (d == '0' ? 0 : COLOR_MISS);
      ntscore = nodes[i-1].ntscore[G];
    }

    if ((score = nodes[i-1].totalscore[T] + (d == '1' ? 0 : COLOR_MISS)) < bestscore) {
      bestscore = score;
      colorscore = nodes[i-1].colorscore[T] + (d == '1' ? 0 : COLOR_MISS);
      ntscore = nodes[i-1].ntscore[T];
    }

    nodes[i].totalscore[G] = bestscore + (g == 'G' ? 0 : NT_MISS);
    nodes[i].colorscore[G] = colorscore;
    nodes[i].ntscore[G] = ntscore + (g == 'G' ? 0 : NT_MISS);

    debug(printf("\t| G: %d=%d+%d",
		 nodes[i].totalscore[G],nodes[i].colorscore[G],nodes[i].ntscore[G]));


    /* T */
    bestscore = nodes[i-1].totalscore[A] + (d == '3' ? 0 : COLOR_MISS);
    colorscore = nodes[i-1].colorscore[A] + (d == '3' ? 0 : COLOR_MISS);
    ntscore = nodes[i-1].ntscore[A];

    if ((score = nodes[i-1].totalscore[C] + (d == '2' ? 0 : COLOR_MISS)) < bestscore) {
      bestscore = score;
      colorscore = nodes[i-1].colorscore[C] + (d == '2' ? 0 : COLOR_MISS);
      ntscore = nodes[i-1].ntscore[C];
    }

    if ((score = nodes[i-1].totalscore[G] + (d == '1' ? 0 : COLOR_MISS)) < bestscore) {
      bestscore = score;
      colorscore = nodes[i-1].colorscore[G] + (d == '1' ? 0 : COLOR_MISS);
      ntscore = nodes[i-1].ntscore[G];
    }

    if ((score = nodes[i-1].totalscore[T] + (d == '0' ? 0 : COLOR_MISS)) < bestscore) {
      bestscore = score;
      colorscore = nodes[i-1].colorscore[T] + (d == '0' ? 0 : COLOR_MISS);
      ntscore = nodes[i-1].ntscore[T];
    }

    nodes[i].totalscore[T] = bestscore + (g == 'T' ? 0 : NT_MISS);
    nodes[i].colorscore[T] = colorscore;
    nodes[i].ntscore[T] = ntscore + (g == 'T' ? 0 : NT_MISS);
    
    debug(printf("\t| T: %d=%d+%d",
		 nodes[i].totalscore[T],nodes[i].colorscore[T],nodes[i].ntscore[T]));

    /* Find overall best at position i */
    bestscore = nodes[i].totalscore[A];
    ntscore = nodes[i].ntscore[A];
    curr = A;
    if ((score = nodes[i].totalscore[C]) < bestscore) {
      bestscore = score;
      ntscore = nodes[i].ntscore[C];
      curr = C;
    } else if (score == bestscore && nodes[i].ntscore[C] < ntscore) {
      bestscore = score;
      ntscore = nodes[i].ntscore[C];
      curr = C;
    }
    if ((score = nodes[i].totalscore[G]) < bestscore) {
      bestscore = score;
      ntscore = nodes[i].ntscore[G];
      curr = G;
    } else if (score == bestscore && nodes[i].ntscore[G] < ntscore) {
      bestscore = score;
      ntscore = nodes[i].ntscore[G];
      curr = G;
    }
    if ((score = nodes[i].totalscore[T]) < bestscore) {
      bestscore = score;
      ntscore = nodes[i].ntscore[T];
      curr = T;
    } else if (score == bestscore && nodes[i].ntscore[T] < ntscore) {
      bestscore = score;
      ntscore = nodes[i].ntscore[T];
      curr = T;
    }

    debug(printf("\t| %c: %d=x+%d\n",
		 GENOME_CHARS[curr],bestscore,ntscore));
  }

  debug(printf("\n"));

  colorscore = nodes[i-1].colorscore[curr];
  *color_nmismatches = colorscore / COLOR_MISS;
  return ntscore / NT_MISS;
}


int
Dibase_count_mismatches_substring (int *color_nmismatches, char *queryuc_ptr, int pos5, int pos3, UINT4 *genome_blocks,
				   Genomicpos_T startpos, Genomicpos_T endpos) {
  char gbuffer[MAX_QUERYLENGTH+1], g, d;
  struct Node_T nodes[MAX_QUERYLENGTH];
  int bestscore, score, colorscore, ntscore;
  int seqlength, i, j;
  Nucleotide_T curr;

  Genome_uncompress_mmap(gbuffer,genome_blocks,startpos,endpos,GENOME_CHARS,GENOME_FLAGS);

  seqlength = endpos - startpos;

  debug({
    printf("Dibase_count_mismatches_substring at %u..%u\n",startpos,endpos);
    printf("gbuffer: ");
    for (i = 0; i < pos5; i++) {
      printf(" ");
    }
    printf("%.*s\n",seqlength,gbuffer);
    printf("query:   %s\n",queryuc_ptr);
  });

  g = gbuffer[0];

  nodes[0].colorscore[A] = 0;
  nodes[0].colorscore[C] = 0;
  nodes[0].colorscore[G] = 0;
  nodes[0].colorscore[T] = 0;

  nodes[0].totalscore[A] = nodes[0].ntscore[A] = (g == 'A' ? 0 : NT_MISS);
  nodes[0].totalscore[C] = nodes[0].ntscore[C] = (g == 'C' ? 0 : NT_MISS);
  nodes[0].totalscore[G] = nodes[0].ntscore[G] = (g == 'G' ? 0 : NT_MISS);
  nodes[0].totalscore[T] = nodes[0].ntscore[T] = (g == 'T' ? 0 : NT_MISS);

  for (i = 1, j = pos5; i < seqlength; i++, j++) {
    g = gbuffer[i];
    d = queryuc_ptr[j];

    debug(printf("q[%d]:%c g[%d]:%c",j,d,i,g));


    /* A */
    bestscore = nodes[i-1].totalscore[A] + (d == '0' ? 0 : COLOR_MISS);
    colorscore = nodes[i-1].colorscore[A] + (d == '0' ? 0 : COLOR_MISS);
    ntscore = nodes[i-1].ntscore[A];

    if ((score = nodes[i-1].totalscore[C] + (d == '1' ? 0 : COLOR_MISS)) < bestscore) {
      bestscore = score;
      colorscore = nodes[i-1].colorscore[C] + (d == '1' ? 0 : COLOR_MISS);
      ntscore = nodes[i-1].ntscore[C];
    }

    if ((score = nodes[i-1].totalscore[G] + (d == '2' ? 0 : COLOR_MISS)) < bestscore) {
      bestscore = score;
      colorscore = nodes[i-1].colorscore[G] + (d == '2' ? 0 : COLOR_MISS);
      ntscore = nodes[i-1].ntscore[G];
    }

    if ((score = nodes[i-1].totalscore[T] + (d == '3' ? 0 : COLOR_MISS)) < bestscore) {
      bestscore = score;
      colorscore = nodes[i-1].colorscore[T] + (d == '3' ? 0 : COLOR_MISS);
      ntscore = nodes[i-1].ntscore[T];
    }

    nodes[i].totalscore[A] = bestscore + (g == 'A' ? 0 : NT_MISS);
    nodes[i].colorscore[A] = colorscore;
    nodes[i].ntscore[A] = ntscore + (g == 'A' ? 0 : NT_MISS);

    debug(printf("\t| A: %d=%d+%d",
		 nodes[i].totalscore[A],nodes[i].colorscore[A],nodes[i].ntscore[A]));

    /* C */
    bestscore = nodes[i-1].totalscore[A] + (d == '1' ? 0 : COLOR_MISS);
    colorscore = nodes[i-1].colorscore[A] + (d == '1' ? 0 : COLOR_MISS);
    ntscore = nodes[i-1].ntscore[A];

    if ((score = nodes[i-1].totalscore[C] + (d == '0' ? 0 : COLOR_MISS)) < bestscore) {
      bestscore = score;
      colorscore = nodes[i-1].colorscore[C] + (d == '0' ? 0 : COLOR_MISS);
      ntscore = nodes[i-1].ntscore[C];
    }

    if ((score = nodes[i-1].totalscore[G] + (d == '3' ? 0 : COLOR_MISS)) < bestscore) {
      bestscore = score;
      colorscore = nodes[i-1].colorscore[G] + (d == '3' ? 0 : COLOR_MISS);
      ntscore = nodes[i-1].ntscore[G];
    }

    if ((score = nodes[i-1].totalscore[T] + (d == '2' ? 0 : COLOR_MISS)) < bestscore) {
      bestscore = score;
      colorscore = nodes[i-1].colorscore[T] + (d == '2' ? 0 : COLOR_MISS);
      ntscore = nodes[i-1].ntscore[T];
    }

    nodes[i].totalscore[C] = bestscore + (g == 'C' ? 0 : NT_MISS);
    nodes[i].colorscore[C] = colorscore;
    nodes[i].ntscore[C] = ntscore + (g == 'C' ? 0 : NT_MISS);

    debug(printf("\t| C: %d=%d+%d",
		 nodes[i].totalscore[C],nodes[i].colorscore[C],nodes[i].ntscore[C]));


    /* G */
    bestscore = nodes[i-1].totalscore[A] + (d == '2' ? 0 : COLOR_MISS);
    colorscore = nodes[i-1].colorscore[A] + (d == '2' ? 0 : COLOR_MISS);
    ntscore = nodes[i-1].ntscore[A];

    if ((score = nodes[i-1].totalscore[C] + (d == '3' ? 0 : COLOR_MISS)) < bestscore) {
      bestscore = score;
      colorscore = nodes[i-1].colorscore[C] + (d == '3' ? 0 : COLOR_MISS);
      ntscore = nodes[i-1].ntscore[C];
    }

    if ((score = nodes[i-1].totalscore[G] + (d == '0' ? 0 : COLOR_MISS)) < bestscore) {
      bestscore = score;
      colorscore = nodes[i-1].colorscore[G] + (d == '0' ? 0 : COLOR_MISS);
      ntscore = nodes[i-1].ntscore[G];
    }

    if ((score = nodes[i-1].totalscore[T] + (d == '1' ? 0 : COLOR_MISS)) < bestscore) {
      bestscore = score;
      colorscore = nodes[i-1].colorscore[T] + (d == '1' ? 0 : COLOR_MISS);
      ntscore = nodes[i-1].ntscore[T];
    }

    nodes[i].totalscore[G] = bestscore + (g == 'G' ? 0 : NT_MISS);
    nodes[i].colorscore[G] = colorscore;
    nodes[i].ntscore[G] = ntscore + (g == 'G' ? 0 : NT_MISS);

    debug(printf("\t| G: %d=%d+%d",
		 nodes[i].totalscore[G],nodes[i].colorscore[G],nodes[i].ntscore[G]));


    /* T */
    bestscore = nodes[i-1].totalscore[A] + (d == '3' ? 0 : COLOR_MISS);
    colorscore = nodes[i-1].colorscore[A] + (d == '3' ? 0 : COLOR_MISS);
    ntscore = nodes[i-1].ntscore[A];

    if ((score = nodes[i-1].totalscore[C] + (d == '2' ? 0 : COLOR_MISS)) < bestscore) {
      bestscore = score;
      colorscore = nodes[i-1].colorscore[C] + (d == '2' ? 0 : COLOR_MISS);
      ntscore = nodes[i-1].ntscore[C];
    }

    if ((score = nodes[i-1].totalscore[G] + (d == '1' ? 0 : COLOR_MISS)) < bestscore) {
      bestscore = score;
      colorscore = nodes[i-1].colorscore[G] + (d == '1' ? 0 : COLOR_MISS);
      ntscore = nodes[i-1].ntscore[G];
    }

    if ((score = nodes[i-1].totalscore[T] + (d == '0' ? 0 : COLOR_MISS)) < bestscore) {
      bestscore = score;
      colorscore = nodes[i-1].colorscore[T] + (d == '0' ? 0 : COLOR_MISS);
      ntscore = nodes[i-1].ntscore[T];
    }

    nodes[i].totalscore[T] = bestscore + (g == 'T' ? 0 : NT_MISS);
    nodes[i].colorscore[T] = colorscore;
    nodes[i].ntscore[T] = ntscore + (g == 'T' ? 0 : NT_MISS);
    
    debug(printf("\t| T: %d=%d+%d",
		 nodes[i].totalscore[T],nodes[i].colorscore[T],nodes[i].ntscore[T]));

    debug(printf("\n"));
  }

  /* Find overall best at position seqlength-1 */
  i = i - 1; /*seqlength-1*/
  bestscore = nodes[i].totalscore[A];
  ntscore = nodes[i].ntscore[A];
  curr = A;
  if ((score = nodes[i].totalscore[C]) < bestscore) {
    bestscore = score;
    ntscore = nodes[i].ntscore[C];
    curr = C;
  } else if (score == bestscore && nodes[i].ntscore[C] < ntscore) {
    bestscore = score;
    ntscore = nodes[i].ntscore[C];
    curr = C;
  }
  if ((score = nodes[i].totalscore[G]) < bestscore) {
    bestscore = score;
    ntscore = nodes[i].ntscore[G];
    curr = G;
  } else if (score == bestscore && nodes[i].ntscore[G] < ntscore) {
    bestscore = score;
    ntscore = nodes[i].ntscore[G];
    curr = G;
  }
  if ((score = nodes[i].totalscore[T]) < bestscore) {
    bestscore = score;
    ntscore = nodes[i].ntscore[T];
    curr = T;
  } else if (score == bestscore && nodes[i].ntscore[T] < ntscore) {
    bestscore = score;
    ntscore = nodes[i].ntscore[T];
    curr = T;
  }

  colorscore = nodes[i].colorscore[curr];
  *color_nmismatches = colorscore / COLOR_MISS;
  return ntscore / NT_MISS;
}


int
Dibase_mismatches_left (int *mismatch_positions, int *colordiffs, int max_mismatches, char *queryuc_ptr,
			int pos5, int pos3, UINT4 *genome_blocks, Genomicpos_T startpos, Genomicpos_T endpos) {
  int nmismatches;
  char gbuffer[MAX_QUERYLENGTH+1], g, d;
  struct Node_T nodes[MAX_QUERYLENGTH];
  int bestscore, score, prevscore, colorscore, ntscore;
  int seqlength, i, j, k;
  Nucleotide_T prev, curr;

  Genome_uncompress_mmap(gbuffer,genome_blocks,startpos,endpos,GENOME_CHARS,GENOME_FLAGS);

  seqlength = endpos - startpos;

  debug1({
    printf("Dibase_mismatches_left at %u..%u, max_mismatches %d\n",startpos,endpos,max_mismatches);
    printf("gbuffer: ");
    for (i = 0; i < pos5; i++) {
      printf(" ");
    }
    printf("%.*s\n",seqlength,gbuffer);
    printf("query:   %s\n",queryuc_ptr);
  });

  g = gbuffer[0];

  nodes[0].colorscore[A] = 0;
  nodes[0].colorscore[C] = 0;
  nodes[0].colorscore[G] = 0;
  nodes[0].colorscore[T] = 0;

  nodes[0].totalscore[A] = nodes[0].ntscore[A] = (g == 'A' ? 0 : NT_MISS);
  nodes[0].totalscore[C] = nodes[0].ntscore[C] = (g == 'C' ? 0 : NT_MISS);
  nodes[0].totalscore[G] = nodes[0].ntscore[G] = (g == 'G' ? 0 : NT_MISS);
  nodes[0].totalscore[T] = nodes[0].ntscore[T] = (g == 'T' ? 0 : NT_MISS);

  bestscore = 0;
  for (i = 1, j = pos5; i < seqlength && bestscore <= max_mismatches; i++, j++) {
    g = gbuffer[i];
    d = queryuc_ptr[j];

    debug1(printf("q[%d]:%c g[%d]:%c",j,d,i,g));


    /* A */
    bestscore = nodes[i-1].totalscore[A] + (d == '0' ? 0 : COLOR_MISS);
    colorscore = nodes[i-1].colorscore[A] + (d == '0' ? 0 : COLOR_MISS);
    ntscore = nodes[i-1].ntscore[A];
    prev = A;

    if ((score = nodes[i-1].totalscore[C] + (d == '1' ? 0 : COLOR_MISS)) < bestscore) {
      bestscore = score;
      colorscore = nodes[i-1].colorscore[C] + (d == '1' ? 0 : COLOR_MISS);
      ntscore = nodes[i-1].ntscore[C];
      prev = C;
    } else if (score == bestscore && nodes[i-1].ntscore[C] < ntscore) {
      bestscore = score;
      colorscore = nodes[i-1].colorscore[C] + (d == '1' ? 0 : COLOR_MISS);
      ntscore = nodes[i-1].ntscore[C];
      prev = C;
    }

    if ((score = nodes[i-1].totalscore[G] + (d == '2' ? 0 : COLOR_MISS)) < bestscore) {
      bestscore = score;
      colorscore = nodes[i-1].colorscore[G] + (d == '2' ? 0 : COLOR_MISS);
      ntscore = nodes[i-1].ntscore[G];
      prev = G;
    } else if (score == bestscore && nodes[i-1].ntscore[G] < ntscore) {
      bestscore = score;
      colorscore = nodes[i-1].colorscore[G] + (d == '2' ? 0 : COLOR_MISS);
      ntscore = nodes[i-1].ntscore[G];
      prev = G;
    }

    if ((score = nodes[i-1].totalscore[T] + (d == '3' ? 0 : COLOR_MISS)) < bestscore) {
      bestscore = score;
      colorscore = nodes[i-1].colorscore[T] + (d == '3' ? 0 : COLOR_MISS);
      ntscore = nodes[i-1].ntscore[T];
      prev = T;
    } else if (score == bestscore && nodes[i-1].ntscore[T] < ntscore) {
      bestscore = score;
      colorscore = nodes[i-1].colorscore[T] + (d == '3' ? 0 : COLOR_MISS);
      ntscore = nodes[i-1].ntscore[T];
      prev = T;
    }

    nodes[i].totalscore[A] = bestscore + (g == 'A' ? 0 : NT_MISS);
    nodes[i].colorscore[A] = colorscore;
    nodes[i].ntscore[A] = ntscore + (g == 'A' ? 0 : NT_MISS);
    nodes[i].prev[A] = prev;

    debug1(printf("\t| A: %d=%d+%d prev:%c",
		 nodes[i].totalscore[A],nodes[i].colorscore[A],nodes[i].ntscore[A],GENOME_CHARS[prev]));


    /* C */
    bestscore = nodes[i-1].totalscore[A] + (d == '1' ? 0 : COLOR_MISS);
    colorscore = nodes[i-1].colorscore[A] + (d == '1' ? 0 : COLOR_MISS);
    ntscore = nodes[i-1].ntscore[A];
    prev = A;

    if ((score = nodes[i-1].totalscore[C] + (d == '0' ? 0 : COLOR_MISS)) < bestscore) {
      bestscore = score;
      colorscore = nodes[i-1].colorscore[C] + (d == '0' ? 0 : COLOR_MISS);
      ntscore = nodes[i-1].ntscore[C];
      prev = C;
    } else if (score == bestscore && nodes[i-1].ntscore[C] < ntscore) {
      bestscore = score;
      colorscore = nodes[i-1].colorscore[C] + (d == '0' ? 0 : COLOR_MISS);
      ntscore = nodes[i-1].ntscore[C];
      prev = C;
    }

    if ((score = nodes[i-1].totalscore[G] + (d == '3' ? 0 : COLOR_MISS)) < bestscore) {
      bestscore = score;
      colorscore = nodes[i-1].colorscore[G] + (d == '3' ? 0 : COLOR_MISS);
      ntscore = nodes[i-1].ntscore[G];
      prev = G;
    } else if (score == bestscore && nodes[i-1].ntscore[G] < ntscore) {
      bestscore = score;
      colorscore = nodes[i-1].colorscore[G] + (d == '3' ? 0 : COLOR_MISS);
      ntscore = nodes[i-1].ntscore[G];
      prev = G;
    }

    if ((score = nodes[i-1].totalscore[T] + (d == '2' ? 0 : COLOR_MISS)) < bestscore) {
      bestscore = score;
      colorscore = nodes[i-1].colorscore[T] + (d == '2' ? 0 : COLOR_MISS);
      ntscore = nodes[i-1].ntscore[T];
      prev = T;
    } else if (score == bestscore && nodes[i-1].ntscore[T] < ntscore) {
      bestscore = score;
      colorscore = nodes[i-1].colorscore[T] + (d == '2' ? 0 : COLOR_MISS);
      ntscore = nodes[i-1].ntscore[T];
      prev = T;
    }

    nodes[i].totalscore[C] = bestscore + (g == 'C' ? 0 : NT_MISS);
    nodes[i].colorscore[C] = colorscore;
    nodes[i].ntscore[C] = ntscore + (g == 'C' ? 0 : NT_MISS);
    nodes[i].prev[C] = prev;

    debug1(printf("\t| C: %d=%d+%d prev:%c",
		 nodes[i].totalscore[C],nodes[i].colorscore[C],nodes[i].ntscore[C],GENOME_CHARS[prev]));

    /* G */
    bestscore = nodes[i-1].totalscore[A] + (d == '2' ? 0 : COLOR_MISS);
    colorscore = nodes[i-1].colorscore[A] + (d == '2' ? 0 : COLOR_MISS);
    ntscore = nodes[i-1].ntscore[A];
    prev = A;

    if ((score = nodes[i-1].totalscore[C] + (d == '3' ? 0 : COLOR_MISS)) < bestscore) {
      bestscore = score;
      colorscore = nodes[i-1].colorscore[C] + (d == '3' ? 0 : COLOR_MISS);
      ntscore = nodes[i-1].ntscore[C];
      prev = C;
    } else if (score == bestscore && nodes[i-1].ntscore[C] < ntscore) {
      bestscore = score;
      colorscore = nodes[i-1].colorscore[C] + (d == '3' ? 0 : COLOR_MISS);
      ntscore = nodes[i-1].ntscore[C];
      prev = C;
    }

    if ((score = nodes[i-1].totalscore[G] + (d == '0' ? 0 : COLOR_MISS)) < bestscore) {
      bestscore = score;
      colorscore = nodes[i-1].colorscore[G] + (d == '0' ? 0 : COLOR_MISS);
      ntscore = nodes[i-1].ntscore[G];
      prev = G;
    } else if (score == bestscore && nodes[i-1].ntscore[G] < ntscore) {
      bestscore = score;
      colorscore = nodes[i-1].colorscore[G] + (d == '0' ? 0 : COLOR_MISS);
      ntscore = nodes[i-1].ntscore[G];
      prev = G;
    }

    if ((score = nodes[i-1].totalscore[T] + (d == '1' ? 0 : COLOR_MISS)) < bestscore) {
      bestscore = score;
      colorscore = nodes[i-1].colorscore[T] + (d == '1' ? 0 : COLOR_MISS);
      ntscore = nodes[i-1].ntscore[T];
      prev = T;
    } else if (score == bestscore && nodes[i-1].ntscore[T] < ntscore) {
      bestscore = score;
      colorscore = nodes[i-1].colorscore[T] + (d == '1' ? 0 : COLOR_MISS);
      ntscore = nodes[i-1].ntscore[T];
      prev = T;
    }

    nodes[i].totalscore[G] = bestscore + (g == 'G' ? 0 : NT_MISS);
    nodes[i].colorscore[G] = colorscore;
    nodes[i].ntscore[G] = ntscore + (g == 'G' ? 0 : NT_MISS);
    nodes[i].prev[G] = prev;

    debug1(printf("\t| G: %d=%d+%d prev:%c",
		 nodes[i].totalscore[G],nodes[i].colorscore[G],nodes[i].ntscore[G],GENOME_CHARS[prev]));

    /* T */
    bestscore = nodes[i-1].totalscore[A] + (d == '3' ? 0 : COLOR_MISS);
    colorscore = nodes[i-1].colorscore[A] + (d == '3' ? 0 : COLOR_MISS);
    ntscore = nodes[i-1].ntscore[A];
    prev = A;

    if ((score = nodes[i-1].totalscore[C] + (d == '2' ? 0 : COLOR_MISS)) < bestscore) {
      bestscore = score;
      colorscore = nodes[i-1].colorscore[C] + (d == '2' ? 0 : COLOR_MISS);
      ntscore = nodes[i-1].ntscore[C];
      prev = C;
    } else if (score == bestscore && nodes[i-1].ntscore[C] < ntscore) {
      bestscore = score;
      colorscore = nodes[i-1].colorscore[C] + (d == '2' ? 0 : COLOR_MISS);
      ntscore = nodes[i-1].ntscore[C];
      prev = C;
    }

    if ((score = nodes[i-1].totalscore[G] + (d == '1' ? 0 : COLOR_MISS)) < bestscore) {
      bestscore = score;
      colorscore = nodes[i-1].colorscore[G] + (d == '1' ? 0 : COLOR_MISS);
      ntscore = nodes[i-1].ntscore[G];
      prev = G;
    } else if (score == bestscore && nodes[i-1].ntscore[G] < ntscore) {
      bestscore = score;
      colorscore = nodes[i-1].colorscore[G] + (d == '1' ? 0 : COLOR_MISS);
      ntscore = nodes[i-1].ntscore[G];
      prev = G;
    }

    if ((score = nodes[i-1].totalscore[T] + (d == '0' ? 0 : COLOR_MISS)) < bestscore) {
      bestscore = score;
      colorscore = nodes[i-1].colorscore[T] + (d == '0' ? 0 : COLOR_MISS);
      ntscore = nodes[i-1].ntscore[T];
      prev = T;
    } else if (score == bestscore && nodes[i-1].ntscore[T] < ntscore) {
      bestscore = score;
      colorscore = nodes[i-1].colorscore[T] + (d == '0' ? 0 : COLOR_MISS);
      ntscore = nodes[i-1].ntscore[T];
      prev = T;
    }

    nodes[i].totalscore[T] = bestscore + (g == 'T' ? 0 : NT_MISS);
    nodes[i].colorscore[T] = colorscore;
    nodes[i].ntscore[T] = ntscore + (g == 'T' ? 0 : NT_MISS);
    nodes[i].prev[T] = prev;

    debug1(printf("\t| T: %d=%d+%d prev:%c",
		 nodes[i].totalscore[T],nodes[i].colorscore[T],nodes[i].ntscore[T],GENOME_CHARS[prev]));

    /* Find overall best at position i */
    bestscore = nodes[i].totalscore[A];
    ntscore = nodes[i].ntscore[A];
    curr = A;
    if ((score = nodes[i].totalscore[C]) < bestscore) {
      bestscore = score;
      ntscore = nodes[i].ntscore[C];
      curr = C;
    } else if (score == bestscore && nodes[i].ntscore[C] < ntscore) {
      bestscore = score;
      ntscore = nodes[i].ntscore[C];
      curr = C;
    }
    if ((score = nodes[i].totalscore[G]) < bestscore) {
      bestscore = score;
      ntscore = nodes[i].ntscore[G];
      curr = G;
    } else if (score == bestscore && nodes[i].ntscore[G] < ntscore) {
      bestscore = score;
      ntscore = nodes[i].ntscore[G];
      curr = G;
    }
    if ((score = nodes[i].totalscore[T]) < bestscore) {
      bestscore = score;
      ntscore = nodes[i].ntscore[T];
      curr = T;
    } else if (score == bestscore && nodes[i].ntscore[T] < ntscore) {
      bestscore = score;
      ntscore = nodes[i].ntscore[T];
      curr = T;
    }

    debug1(printf("\t| %c: %d=x+%d\n",
		 GENOME_CHARS[curr],bestscore,ntscore));
  }

  debug1(printf("\n"));

  k = nmismatches = prevscore = bestscore;
  debug1(printf("%d mismatches:\n",k));
  for (i = i - 1, j = j - 1; i >= 1; i--, j--) {
    g = gbuffer[i];
    d = queryuc_ptr[j];

    debug1(printf("q[%d]:%c g[%d]:%c",j,d,i,g));

    prev = nodes[i].prev[curr];
    debug1(printf("\t| %c: %d=%d+%d prev:%c",
		 GENOME_CHARS[curr],nodes[i].totalscore[curr],nodes[i].colorscore[curr],nodes[i].ntscore[curr],
		 GENOME_CHARS[prev]));

    score = nodes[i].totalscore[curr];
    if (score < prevscore) {
      mismatch_positions[--k] = j+1; /* Need +1 to be consistent with
					  nt sequence, which starts with
					  an assumed nt */
      debug1(printf("  Mismatch: Putting %d+1 in index %d",j,k));
    }
#if 0
    debug1(printf("  Putting %d in colordiffs for %d",colordiffs[k],k));
#endif
    debug1(printf("\n"));

    curr = prev;
    prevscore = score;
  }

  if (k == 1) {
    debug1(printf("Putting %d in index 0\n",pos5));
    mismatch_positions[0] = pos5;
  }

  return nmismatches;
}


int
Dibase_mismatches_right (int *mismatch_positions, int *colordiffs, int max_mismatches, char *queryuc_ptr,
			 int pos5, int pos3, UINT4 *genome_blocks, Genomicpos_T startpos, Genomicpos_T endpos) {
  int nmismatches;
  char gbuffer[MAX_QUERYLENGTH+1], g, d;
  struct Node_T nodes[MAX_QUERYLENGTH];
  int bestscore, score, prevscore, colorscore, ntscore;
  int seqlength, i, j, k;
  Nucleotide_T prev, curr;

  Genome_uncompress_mmap(gbuffer,genome_blocks,startpos,endpos,GENOME_CHARS,GENOME_FLAGS);

  seqlength = endpos - startpos;

  debug2({
    printf("Dibase_mismatches_right at %u..%u, max_mismatches %d\n",startpos,endpos,max_mismatches);
    printf("gbuffer: ");
    for (i = 0; i < pos5; i++) {
      printf(" ");
    }
    printf("%.*s\n",seqlength,gbuffer);
    printf("query:   %s\n",queryuc_ptr);
  });

  i = seqlength;

  nodes[i].colorscore[A] = 0;
  nodes[i].colorscore[C] = 0;
  nodes[i].colorscore[G] = 0;
  nodes[i].colorscore[T] = 0;

  nodes[i].totalscore[A] = nodes[i].ntscore[A] = 0;
  nodes[i].totalscore[C] = nodes[i].ntscore[C] = 0;
  nodes[i].totalscore[G] = nodes[i].ntscore[G] = 0;
  nodes[i].totalscore[T] = nodes[i].ntscore[T] = 0;

  bestscore = 0;
  for (i = seqlength - 1, j = pos3-1; i >= 0 && bestscore <= max_mismatches; i--, j--) {
    g = gbuffer[i];
    d = queryuc_ptr[j];

    debug2(printf("q[%d]:%c g[%d]:%c",j,d,i,g));


    /* A */
    bestscore = nodes[i+1].totalscore[A] + (d == '0' ? 0 : COLOR_MISS);
    colorscore = nodes[i+1].colorscore[A] + (d == '0' ? 0 : COLOR_MISS);
    ntscore = nodes[i+1].ntscore[A];
    prev = A;

    if ((score = nodes[i+1].totalscore[C] + (d == '1' ? 0 : COLOR_MISS)) < bestscore) {
      bestscore = score;
      colorscore = nodes[i+1].colorscore[C] + (d == '1' ? 0 : COLOR_MISS);
      ntscore = nodes[i+1].ntscore[C];
      prev = C;
    } else if (score == bestscore && nodes[i+1].ntscore[C] < ntscore) {
      bestscore = score;
      colorscore = nodes[i+1].colorscore[C] + (d == '1' ? 0 : COLOR_MISS);
      ntscore = nodes[i+1].ntscore[C];
      prev = C;
    }

    if ((score = nodes[i+1].totalscore[G] + (d == '2' ? 0 : COLOR_MISS)) < bestscore) {
      bestscore = score;
      colorscore = nodes[i+1].colorscore[G] + (d == '2' ? 0 : COLOR_MISS);
      ntscore = nodes[i+1].ntscore[G];
      prev = G;
    } else if (score == bestscore && nodes[i+1].ntscore[G] < ntscore) {
      bestscore = score;
      colorscore = nodes[i+1].colorscore[G] + (d == '2' ? 0 : COLOR_MISS);
      ntscore = nodes[i+1].ntscore[G];
      prev = G;
    }

    if ((score = nodes[i+1].totalscore[T] + (d == '3' ? 0 : COLOR_MISS)) < bestscore) {
      bestscore = score;
      colorscore = nodes[i+1].colorscore[T] + (d == '3' ? 0 : COLOR_MISS);
      ntscore = nodes[i+1].ntscore[T];
      prev = T;
    } else if (score == bestscore && nodes[i+1].ntscore[T] < ntscore) {
      bestscore = score;
      colorscore = nodes[i+1].colorscore[T] + (d == '3' ? 0 : COLOR_MISS);
      ntscore = nodes[i+1].ntscore[T];
      prev = T;
    }

    nodes[i].totalscore[A] = bestscore + (g == 'A' ? 0 : NT_MISS);
    nodes[i].colorscore[A] = colorscore;
    nodes[i].ntscore[A] = ntscore + (g == 'A' ? 0 : NT_MISS);
    nodes[i].prev[A] = prev;

    debug2(printf("\t| A: %d=%d+%d prev:%c",
		 nodes[i].totalscore[A],nodes[i].colorscore[A],nodes[i].ntscore[A],GENOME_CHARS[prev]));


    /* C */
    bestscore = nodes[i+1].totalscore[A] + (d == '1' ? 0 : COLOR_MISS);
    colorscore = nodes[i+1].colorscore[A] + (d == '1' ? 0 : COLOR_MISS);
    ntscore = nodes[i+1].ntscore[A];
    prev = A;

    if ((score = nodes[i+1].totalscore[C] + (d == '0' ? 0 : COLOR_MISS)) < bestscore) {
      bestscore = score;
      colorscore = nodes[i+1].colorscore[C] + (d == '0' ? 0 : COLOR_MISS);
      ntscore = nodes[i+1].ntscore[C];
      prev = C;
    } else if (score == bestscore && nodes[i+1].ntscore[C] < ntscore) {
      bestscore = score;
      colorscore = nodes[i+1].colorscore[C] + (d == '0' ? 0 : COLOR_MISS);
      ntscore = nodes[i+1].ntscore[C];
      prev = C;
    }

    if ((score = nodes[i+1].totalscore[G] + (d == '3' ? 0 : COLOR_MISS)) < bestscore) {
      bestscore = score;
      colorscore = nodes[i+1].colorscore[G] + (d == '3' ? 0 : COLOR_MISS);
      ntscore = nodes[i+1].ntscore[G];
      prev = G;
    } else if (score == bestscore && nodes[i+1].ntscore[G] < ntscore) {
      bestscore = score;
      colorscore = nodes[i+1].colorscore[G] + (d == '3' ? 0 : COLOR_MISS);
      ntscore = nodes[i+1].ntscore[G];
      prev = G;
    }

    if ((score = nodes[i+1].totalscore[T] + (d == '2' ? 0 : COLOR_MISS)) < bestscore) {
      bestscore = score;
      colorscore = nodes[i+1].colorscore[T] + (d == '2' ? 0 : COLOR_MISS);
      ntscore = nodes[i+1].ntscore[T];
      prev = T;
    } else if (score == bestscore && nodes[i+1].ntscore[T] < ntscore) {
      bestscore = score;
      colorscore = nodes[i+1].colorscore[T] + (d == '2' ? 0 : COLOR_MISS);
      ntscore = nodes[i+1].ntscore[T];
      prev = T;
    }

    nodes[i].totalscore[C] = bestscore + (g == 'C' ? 0 : NT_MISS);
    nodes[i].colorscore[C] = colorscore;
    nodes[i].ntscore[C] = ntscore + (g == 'C' ? 0 : NT_MISS);
    nodes[i].prev[C] = prev;

    debug2(printf("\t| C: %d=%d+%d prev:%c",
		 nodes[i].totalscore[C],nodes[i].colorscore[C],nodes[i].ntscore[C],GENOME_CHARS[prev]));

    /* G */
    bestscore = nodes[i+1].totalscore[A] + (d == '2' ? 0 : COLOR_MISS);
    colorscore = nodes[i+1].colorscore[A] + (d == '2' ? 0 : COLOR_MISS);
    ntscore = nodes[i+1].ntscore[A];
    prev = A;

    if ((score = nodes[i+1].totalscore[C] + (d == '3' ? 0 : COLOR_MISS)) < bestscore) {
      bestscore = score;
      colorscore = nodes[i+1].colorscore[C] + (d == '3' ? 0 : COLOR_MISS);
      ntscore = nodes[i+1].ntscore[C];
      prev = C;
    } else if (score == bestscore && nodes[i+1].ntscore[C] < ntscore) {
      bestscore = score;
      colorscore = nodes[i+1].colorscore[C] + (d == '3' ? 0 : COLOR_MISS);
      ntscore = nodes[i+1].ntscore[C];
      prev = C;
    }

    if ((score = nodes[i+1].totalscore[G] + (d == '0' ? 0 : COLOR_MISS)) < bestscore) {
      bestscore = score;
      colorscore = nodes[i+1].colorscore[G] + (d == '0' ? 0 : COLOR_MISS);
      ntscore = nodes[i+1].ntscore[G];
      prev = G;
    } else if (score == bestscore && nodes[i+1].ntscore[G] < ntscore) {
      bestscore = score;
      colorscore = nodes[i+1].colorscore[G] + (d == '0' ? 0 : COLOR_MISS);
      ntscore = nodes[i+1].ntscore[G];
      prev = G;
    }

    if ((score = nodes[i+1].totalscore[T] + (d == '1' ? 0 : COLOR_MISS)) < bestscore) {
      bestscore = score;
      colorscore = nodes[i+1].colorscore[T] + (d == '1' ? 0 : COLOR_MISS);
      ntscore = nodes[i+1].ntscore[T];
      prev = T;
    } else if (score == bestscore && nodes[i+1].ntscore[T] < ntscore) {
      bestscore = score;
      colorscore = nodes[i+1].colorscore[T] + (d == '1' ? 0 : COLOR_MISS);
      ntscore = nodes[i+1].ntscore[T];
      prev = T;
    }

    nodes[i].totalscore[G] = bestscore + (g == 'G' ? 0 : NT_MISS);
    nodes[i].colorscore[G] = colorscore;
    nodes[i].ntscore[G] = ntscore + (g == 'G' ? 0 : NT_MISS);
    nodes[i].prev[G] = prev;

    debug2(printf("\t| G: %d=%d+%d prev:%c",
		 nodes[i].totalscore[G],nodes[i].colorscore[G],nodes[i].ntscore[G],GENOME_CHARS[prev]));

    /* T */
    bestscore = nodes[i+1].totalscore[A] + (d == '3' ? 0 : COLOR_MISS);
    colorscore = nodes[i+1].colorscore[A] + (d == '3' ? 0 : COLOR_MISS);
    ntscore = nodes[i+1].ntscore[A];
    prev = A;

    if ((score = nodes[i+1].totalscore[C] + (d == '2' ? 0 : COLOR_MISS)) < bestscore) {
      bestscore = score;
      colorscore = nodes[i+1].colorscore[C] + (d == '2' ? 0 : COLOR_MISS);
      ntscore = nodes[i+1].ntscore[C];
      prev = C;
    } else if (score == bestscore && nodes[i+1].ntscore[C] < ntscore) {
      bestscore = score;
      colorscore = nodes[i+1].colorscore[C] + (d == '2' ? 0 : COLOR_MISS);
      ntscore = nodes[i+1].ntscore[C];
      prev = C;
    }

    if ((score = nodes[i+1].totalscore[G] + (d == '1' ? 0 : COLOR_MISS)) < bestscore) {
      bestscore = score;
      colorscore = nodes[i+1].colorscore[G] + (d == '1' ? 0 : COLOR_MISS);
      ntscore = nodes[i+1].ntscore[G];
      prev = G;
    } else if (score == bestscore && nodes[i+1].ntscore[G] < ntscore) {
      bestscore = score;
      colorscore = nodes[i+1].colorscore[G] + (d == '1' ? 0 : COLOR_MISS);
      ntscore = nodes[i+1].ntscore[G];
      prev = G;
    }

    if ((score = nodes[i+1].totalscore[T] + (d == '0' ? 0 : COLOR_MISS)) < bestscore) {
      bestscore = score;
      colorscore = nodes[i+1].colorscore[T] + (d == '0' ? 0 : COLOR_MISS);
      ntscore = nodes[i+1].ntscore[T];
      prev = T;
    } else if (score == bestscore && nodes[i+1].ntscore[T] < ntscore) {
      bestscore = score;
      colorscore = nodes[i+1].colorscore[T] + (d == '0' ? 0 : COLOR_MISS);
      ntscore = nodes[i+1].ntscore[T];
      prev = T;
    }

    nodes[i].totalscore[T] = bestscore + (g == 'T' ? 0 : NT_MISS);
    nodes[i].colorscore[T] = colorscore;
    nodes[i].ntscore[T] = ntscore + (g == 'T' ? 0 : NT_MISS);
    nodes[i].prev[T] = prev;

    debug2(printf("\t| T: %d=%d+%d prev:%c",
		 nodes[i].totalscore[T],nodes[i].colorscore[T],nodes[i].ntscore[T],GENOME_CHARS[prev]));

    /* Find overall best at position i */
    bestscore = nodes[i].totalscore[A];
    ntscore = nodes[i].ntscore[A];
    curr = A;
    if ((score = nodes[i].totalscore[C]) < bestscore) {
      bestscore = score;
      ntscore = nodes[i].ntscore[C];
      curr = C;
    } else if (score == bestscore && nodes[i].ntscore[C] < ntscore) {
      bestscore = score;
      ntscore = nodes[i].ntscore[C];
      curr = C;
    }
    if ((score = nodes[i].totalscore[G]) < bestscore) {
      bestscore = score;
      ntscore = nodes[i].ntscore[G];
      curr = G;
    } else if (score == bestscore && nodes[i].ntscore[G] < ntscore) {
      bestscore = score;
      ntscore = nodes[i].ntscore[G];
      curr = G;
    }
    if ((score = nodes[i].totalscore[T]) < bestscore) {
      bestscore = score;
      ntscore = nodes[i].ntscore[T];
      curr = T;
    } else if (score == bestscore && nodes[i].ntscore[T] < ntscore) {
      bestscore = score;
      ntscore = nodes[i].ntscore[T];
      curr = T;
    }

    debug2(printf("\t| %c: %d=x+%d\n",
		 GENOME_CHARS[curr],bestscore,ntscore));
  }

  debug2(printf("\n"));


  k = nmismatches = prevscore = bestscore;
  debug2(printf("%d mismatches:\n",k));
  for (i = i + 1, j = j + 1; i < seqlength; i++, j++) {
    g = gbuffer[i];
    d = queryuc_ptr[j];

    debug2(printf("q[%d]:%c g[%d]:%c",j,d,i,g));

    prev = nodes[i].prev[curr];
    debug2(printf("\t| %c: %d=%d+%d prev:%c",
		 GENOME_CHARS[curr],nodes[i].totalscore[curr],nodes[i].colorscore[curr],nodes[i].ntscore[curr],
		 GENOME_CHARS[prev]));

    score = nodes[i].totalscore[curr];
    if (score < prevscore) {
      mismatch_positions[--k] = j;
      debug2(printf("  Mismatch: Putting %d in index %d",j,k));
    }
#if 0    
    debug2(printf("  Putting %d in colordiffs for %d",colordiffs[k],k));
#endif
    debug2(printf("\n"));

    curr = prev;
    prevscore = score;
  }

  if (k == 1) {
    debug1(printf("Putting %d in index 0\n",pos3));
    mismatch_positions[0] = pos3;
  }

  return nmismatches;
}


