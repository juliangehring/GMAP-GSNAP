static char rcsid[] = "$Id: chimera.c,v 1.7 2005/05/06 17:03:03 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "chimera.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>		/* For sqrt */
#include "mem.h"
#include "nmath.h"


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Chimera detection */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif



#define NPSEUDO 10.0

double
Chimera_detect (int *breakpoint, int *margin, Stage3_T stage3, Sequence_T queryseq) {
  double pvalue = 1.0;
  int *matchscores;
  int querylength, leftmargin, rightmargin;

  /* x signifies nmatches, y signifies nmismatches, x + y = n */
  int start, end, pos, x = 0, y, n, x_left, y_left, x_right, y_right, n_left, n_right;
  double theta, x_pseudo, theta_left, theta_right, rss, rss_left, rss_right, rss_sep;
  double fscore, tstat, min_rss_sep, best_pos = -1;

  /*
  start = Sequence_trim_start(queryseq);
  end = Sequence_trim_end(queryseq);
  */

  start = Stage3_querystart(stage3);
  end = Stage3_queryend(stage3);
  querylength = Sequence_fulllength(queryseq);

  matchscores = Stage3_matchscores(stage3,querylength);

  x = 0;
  for (pos = start; pos < end; pos++) {
    x += matchscores[pos];
  }
  n = end - start;
  y = n - x;

  /* when rss_sep == rss, fscore == 0 */
  min_rss_sep = rss = (double) x * (double) y/(double) n;
  if (rss == 0.0) {
    *breakpoint = -1;
    FREE(matchscores);
    return 1.0;
  }

  theta = (double) x/(double) n;
  x_pseudo = NPSEUDO * theta;
  debug1(printf("%d %d %d %f\n",x,y,n,theta));
  
  x_left = y_left = n_left = 0;
  x_right = x;
  y_right = y;
  n_right = n;

  debug1(printf("%s %s %s %s %s %s %s %s %s %s %s %s %s\n",
		"pos","match","x.left","y.left","n.left","x.right","y.right","n.right",
		"theta.left","theta.right","rss.left","rss.right","fscore"));

  for (pos = start; pos < end-1; pos++) {
    if (matchscores[pos] == 1) {
      x_left++;
      x_right--;
    } else if (matchscores[pos] == 0) {
      y_left++;
      y_right--;
    } else {
      abort();
    }
    n_left++;
    n_right--;
    
    theta_left = ((double) x_left + x_pseudo)/((double) n_left + NPSEUDO);
    theta_right = ((double) x_right + x_pseudo)/((double) n_right + NPSEUDO);
    rss_left = x_left*(1.0-theta_left)*(1.0-theta_left) + y_left*theta_left*theta_left;
    rss_right = x_right*(1.0-theta_right)*(1.0-theta_right) + y_right*theta_right*theta_right;
    rss_sep = rss_left + rss_right;

    if (rss_sep == 0) {
      debug1(printf("%d %d %d %d %d %d %d %d %f %f %f %f NA\n",
		    pos,matchscores[pos],x_left,y_left,n_left,x_right,y_right,n_right,
		    theta_left,theta_right,rss_left,rss_right));
    } else {
      debug1(      
	     fscore = ((double) (n - 2))*(rss - rss_sep)/rss_sep;
	     printf("%d %d %d %d %d %d %d %d %f %f %f %f %f\n",
		    pos,matchscores[pos],x_left,y_left,n_left,x_right,y_right,n_right,
		    theta_left,theta_right,rss_left,rss_right,fscore);
	     );

      /* fscore = (n-2)*(rss - rss_sep)/rss_sep = (n-2)*(rss/rss_sep -
	 1) is maximized when rss_sep is minimized */

      if (rss_sep < min_rss_sep) {
	min_rss_sep = rss_sep;
	best_pos = pos;
      }
    }
  }

  *breakpoint = best_pos;
  *margin = 0;

  fscore = ((double) (n - 2))*(rss - min_rss_sep)/min_rss_sep;
  tstat = sqrt(fscore);
  pvalue = 2*Nmath_pnormc(tstat); /* Two-sided */
  if (pvalue > 1.0) {
    pvalue = 1.0;
  }
  if ((leftmargin = *breakpoint - Sequence_trim_start(queryseq)) < 0) {
    leftmargin = 0;
  }
  if ((rightmargin = Sequence_trim_end(queryseq) - *breakpoint) < 0) {
    rightmargin = 0;
  }

  /* Return smaller margin */
  if (leftmargin < rightmargin) {
    *margin = leftmargin;
  } else {
    *margin = rightmargin;
  }

  FREE(matchscores);

  debug1(printf("at %d (margin = %d), tstat = %f, pvalue = %g\n",*breakpoint,*margin,tstat,pvalue));
  return pvalue;
}



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

  debug(
	for (pos = 0; pos < querylength; pos++) {
	  printf("%d:",pos);
	  for (i = 0; i < npaths; i++) {
	    printf("\t%d",matrix[i][pos]);
	  }
	  printf("\n");
	}
	);

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

  debug(printf("From path %d to path %d at pos %d.  Score = %d.  Fromscore_3 = %d, Toscore_5 = %d\n",
	       *bestfrom,*bestto,*bestpos,bestscore,*fromscore_3,*toscore_5));


  FREE(bestscoreatpos);
  FREE(to);
  FREE(from);

  for (i = 0; i < npaths; i++) {
    FREE(matrix[i]);
  }
  FREE(matrix);

  return;
}
