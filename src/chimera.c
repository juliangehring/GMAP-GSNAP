static char rcsid[] = "$Id: chimera.c,v 1.3 2005/02/08 00:02:27 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "chimera.h"
#include <stdio.h>
#include <stdlib.h>
#include "mem.h"

/*
#define DEBUG 1
*/

void
Chimera_bestpath (int *bestfrom, int *bestto, int *bestpos,
		  int *fromscore_3, int *toscore_5,
		  Stage3_T *stage3array, int npaths, int querylength) {
  int **matrix, *from, *to, *bestscoreatpos, i, j, pos, score, 
    bestscore;

  matrix = (int **) CALLOC(npaths,sizeof(int *));
  for (i = 0; i < npaths; i++) {
    matrix[i] = (int *) CALLOC(querylength,sizeof(int));
    Stage3_pathscores(matrix[i],stage3array[i],querylength);
  }

#ifdef DEBUG  
  for (pos = 0; pos < querylength; pos++) {
    printf("%d:",pos);
    for (i = 0; i < npaths; i++) {
      printf("\t%d",matrix[i][pos]);
    }
    printf("\n");
  }
#endif

  from = (int *) CALLOC(querylength,sizeof(int));
  to = (int *) CALLOC(querylength,sizeof(int));
  bestscoreatpos = (int *) CALLOC(querylength,sizeof(int));
  for (pos = 0; pos < querylength; pos++) {
    bestscoreatpos[pos] = -100000;
  }

  for (pos = 0; pos < querylength; pos++) {
    for (i = 0; i < npaths; i++) {
      for (j = 0; j < npaths; j++) {
	score = matrix[j][querylength-1] - matrix[j][pos] + matrix[i][pos] /* - 0 */;
	if (score > bestscoreatpos[pos]) {
	  bestscoreatpos[pos] = score;
	  from[pos] = i;
	  to[pos] = j;
	}
      }
    }
  }

  bestscore = -100000;
  for (pos = 0; pos < querylength; pos++) {
    if (bestscoreatpos[pos] > bestscore) {
      bestscore = bestscoreatpos[pos];
      *bestpos = pos;
      *bestfrom = from[pos];
      *bestto = to[pos];
    }
  }

  *fromscore_3 = matrix[*bestfrom][querylength-1] - matrix[*bestfrom][*bestpos];
  *toscore_5 = matrix[*bestto][*bestpos];

#ifdef DEBUG
  printf("From path %d to path %d at pos %d.  Score = %d.  Fromscore_3 = %d, Toscore_5 = %d\n",
	 *bestfrom,*bestto,*bestpos,bestscore,*fromscore_3,*toscore_5);
#endif


  FREE(bestscoreatpos);
  FREE(to);
  FREE(from);

  for (i = 0; i < npaths; i++) {
    FREE(matrix[i]);
  }
  FREE(matrix);

  return;
}
