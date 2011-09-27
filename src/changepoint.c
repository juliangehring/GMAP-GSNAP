static char rcsid[] = "$Id: changepoint.c,v 1.3 2006/05/23 06:42:59 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "changepoint.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>		/* For sqrt */

#define NPSEUDO 12.0
#define THETARATIO 0.67

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


/* Goes from 3' to 5', looking for the last sharp decrease in matches */
int
Changepoint_start_change (int start, int end, int *matchscores) {
  int edge = start;

  /* x signifies nmatches, y signifies nmismatches, x + y = n */
  int pos, x = 0, y, n, x_past, y_past, x_future, y_future, n_past, n_future;
  double theta, x_pseudo, theta_past, theta_future, rss, rss_past, rss_future, rss_sep;
  double fscore, min_rss_sep, best_pos = -1, best_theta_past, best_theta_future;

  x = 0;
  for (pos = start; pos < end; pos++) {
    x += matchscores[pos];
  }
  n = end - start;
  y = n - x;

  /* when rss_sep == rss, fscore == 0 */
  min_rss_sep = rss = (double) x * (double) y/(double) n;
  if (rss == 0.0) {
    return start;
  }

  theta = (double) x/(double) n;
  x_pseudo = NPSEUDO * theta;
  debug(printf("%d %d %d %f\n",x,y,n,theta));
  
  x_past = y_past = n_past = 0;
  x_future = x;
  y_future = y;
  n_future = n;

  debug(printf("%s %s %s %s %s %s %s %s %s %s %s %s %s\n",
		"pos","match","x.past","y.past","n.past","x.future","y.future","n.future",
		"theta.past","theta.future","rss.past","rss.future","fscore"));

  for (pos = end-1; pos > start; --pos) {
    if (matchscores[pos] == 1) {
      x_past++;
      x_future--;
    } else if (matchscores[pos] == 0) {
      y_past++;
      y_future--;
    } else {
      abort();
    }
    n_past++;
    n_future--;

    theta_past = ((double) x_past + x_pseudo)/((double) n_past + NPSEUDO);
    theta_future = ((double) x_future + x_pseudo)/((double) n_future + NPSEUDO);
    rss_past = x_past*(1.0-theta_past)*(1.0-theta_past) + y_past*theta_past*theta_past;
    rss_future = x_future*(1.0-theta_future)*(1.0-theta_future) + y_future*theta_future*theta_future;
    rss_sep = rss_past + rss_future;

    if (rss_sep == 0) {
      debug(printf("%d %d %d %d %d %d %d %d %f %f %f %f NA\n",
		    pos,matchscores[pos],x_past,y_past,n_past,x_future,y_future,n_future,
		    theta_past,theta_future,rss_past,rss_future));

#if 0
    } else if (theta_future < theta_past) {
      debug(printf("Changepoint_start_change aborting early with edge=%d\n",edge));
      return edge;
#endif

    } else {
      debug(      
	     fscore = ((double) (n - 2))*(rss - rss_sep)/rss_sep;
	     printf("%d %d %d %d %d %d %d %d %f %f %f %f %f\n",
		    pos,matchscores[pos],x_past,y_past,n_past,x_future,y_future,n_future,
		    theta_past,theta_future,rss_past,rss_future,fscore);
	     );

      /* fscore = (n-2)*(rss - rss_sep)/rss_sep = (n-2)*(rss/rss_sep -
	 1) is maximized when rss_sep is minimized */

      if (theta_future/theta_past < THETARATIO) {
	debug(printf("thetadiff = %f, thetaratio = %f\n",theta_past-theta_future,theta_future/theta_past));      
	if (rss_sep < min_rss_sep) {
	  min_rss_sep = rss_sep;
	  edge = pos;
	  debug(printf("Set start_change to %d\n",pos));      
	}
      }
    }
  }

  debug(printf("Changepoint_start_change returning %d\n",edge));
  return edge;
}


/* Goes from 5' to 3', looking for the last sharp decrease in matches */
int
Changepoint_end_change (int start, int end, int *matchscores) {
  int edge = end;

  /* x signifies nmatches, y signifies nmismatches, x + y = n */
  int pos, x = 0, y, n, x_past, y_past, x_future, y_future, n_past, n_future;
  double theta, x_pseudo, theta_past, theta_future, rss, rss_past, rss_future, rss_sep;
  double fscore, min_rss_sep, best_pos = -1, best_theta_past, best_theta_future;

  x = 0;
  for (pos = start; pos < end; pos++) {
    x += matchscores[pos];
  }
  n = end - start;
  y = n - x;

  /* when rss_sep == rss, fscore == 0 */
  min_rss_sep = rss = (double) x * (double) y/(double) n;
  if (rss == 0.0) {
    return end;
  }

  theta = (double) x/(double) n;
  x_pseudo = NPSEUDO * theta;
  debug(printf("%d %d %d %f\n",x,y,n,theta));
  
  x_past = y_past = n_past = 0;
  x_future = x;
  y_future = y;
  n_future = n;

  debug(printf("%s %s %s %s %s %s %s %s %s %s %s %s %s\n",
		"pos","match","x.past","y.past","n.past","x.future","y.future","n.future",
		"theta.past","theta.future","rss.past","rss.future","fscore"));

  for (pos = start+1; pos < end; pos++) {
    if (matchscores[pos] == 1) {
      x_past++;
      x_future--;
    } else if (matchscores[pos] == 0) {
      y_past++;
      y_future--;
    } else {
      abort();
    }
    n_past++;
    n_future--;

    theta_past = ((double) x_past + x_pseudo)/((double) n_past + NPSEUDO);
    theta_future = ((double) x_future + x_pseudo)/((double) n_future + NPSEUDO);
    rss_past = x_past*(1.0-theta_past)*(1.0-theta_past) + y_past*theta_past*theta_past;
    rss_future = x_future*(1.0-theta_future)*(1.0-theta_future) + y_future*theta_future*theta_future;
    rss_sep = rss_past + rss_future;

    if (rss_sep == 0) {
      debug(printf("%d %d %d %d %d %d %d %d %f %f %f %f NA\n",
		    pos,matchscores[pos],x_past,y_past,n_past,x_future,y_future,n_future,
		    theta_past,theta_future,rss_past,rss_future));

#if 0
    } else if (theta_past < theta_future) {
      debug(printf("Changepoint_end_change aborting early with edge=%d\n",edge));
      return edge;
#endif

    } else {
      debug(      
	     fscore = ((double) (n - 2))*(rss - rss_sep)/rss_sep;
	     printf("%d %d %d %d %d %d %d %d %f %f %f %f %f\n",
		    pos,matchscores[pos],x_past,y_past,n_past,x_future,y_future,n_future,
		    theta_past,theta_future,rss_past,rss_future,fscore);
	     );

      /* fscore = (n-2)*(rss - rss_sep)/rss_sep = (n-2)*(rss/rss_sep -
	 1) is maximized when rss_sep is minimized */

      if (theta_future/theta_past < THETARATIO) {
	debug(printf("thetadiff = %f, thetaratio = %f\n",theta_past-theta_future,theta_future/theta_past));      
	if (rss_sep < min_rss_sep) {
	  min_rss_sep = rss_sep;
	  edge = pos;
	  debug(printf("Set end_change to %d\n",pos));      
	}
      }
    }
  }

  debug(printf("Changepoint_end_change returning %d\n",edge));
  return edge;
}




#if 0
int
Changepoint_middle_change (int *newstart, int *newend, Stage3_T stage3, int queryntlength, double fthreshold) {
  int breakpoint;
  double pvalue = 1.0;
  int *matchscores;

  /* x signifies nmatches, y signifies nmismatches, x + y = n */
  int start, end, pos, x = 0, y, n, x_left, y_left, x_right, y_right, n_left, n_right;
  double theta, x_pseudo, theta_left, theta_right, rss, rss_left, rss_right, rss_sep;
  double fscore, min_rss_sep, best_pos = -1, best_theta_left, best_theta_right;

  /*
  start = Sequence_trim_start(queryseq);
  end = Sequence_trim_end(queryseq);
  */

  start = Stage3_querystart(stage3);
  end = Stage3_queryend(stage3);

  matchscores = Stage3_matchscores(stage3,queryntlength);

  x = 0;
  for (pos = start; pos < end; pos++) {
    x += matchscores[pos];
  }
  n = end - start;
  y = n - x;

  /* when rss_sep == rss, fscore == 0 */
  min_rss_sep = rss = (double) x * (double) y/(double) n;
  if (rss == 0.0) {
    FREE(matchscores);
    return 0;
  }

  theta = (double) x/(double) n;
  x_pseudo = NPSEUDO * theta;
  debug(printf("%d %d %d %f\n",x,y,n,theta));
  
  x_left = y_left = n_left = 0;
  x_right = x;
  y_right = y;
  n_right = n;

  debug(printf("%s %s %s %s %s %s %s %s %s %s %s %s %s\n",
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
      debug(printf("%d %d %d %d %d %d %d %d %f %f %f %f NA\n",
		    pos,matchscores[pos],x_left,y_left,n_left,x_right,y_right,n_right,
		    theta_left,theta_right,rss_left,rss_right));
    } else {
      debug(      
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
	best_theta_left = theta_left;
	best_theta_right = theta_right;
      }
    }
  }
  FREE(matchscores);

  fscore = ((double) (n - 2))*(rss - min_rss_sep)/min_rss_sep;
  if (fscore < fthreshold) {
    return 0;
  } else {
    breakpoint = best_pos;
    debug(printf("at %d, fscore = %f\n",breakpoint,fscore));
    if (best_theta_left < best_theta_right) {
      /* trim left */
      *newstart = breakpoint;
      *newend = end;
      return breakpoint - start;
    } else {
      /* trim right */
      *newstart = start;
      *newend = breakpoint;
      return end - breakpoint;
    }
  }
}

#endif

