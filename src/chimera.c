static char rcsid[] = "$Id: chimera.c,v 1.13 2005/06/14 17:06:40 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "chimera.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>		/* For sqrt */
#include "mem.h"
#include "nmath.h"
#include "maxent.h"


/* Chimera assembly */
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

/* Exon-exon boundary detection */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif



#define T Chimera_T
struct T {
  int nonchimera_matches;
  int nonchimera_mismatches;
  int nonchimera_indels;

  int chimerapos;
  int equivpos;

  int cdna_direction;
  int exonexonpos;
  double donor_prob;
  double acceptor_prob;
};


int 
Chimera_pos (T this) {
  return this->chimerapos;
}

int
Chimera_equivpos (T this) {
  return this->equivpos;
}

int
Chimera_cdna_direction (T this) {
  return this->cdna_direction;
}

double
Chimera_donor_prob (T this) {
  if (this->exonexonpos < 0) {
    return 0.0;
  } else {
    return this->donor_prob;
  }
}

double
Chimera_acceptor_prob (T this) {
  if (this->exonexonpos < 0) {
    return 0.0;
  } else {
    return this->acceptor_prob;
  }
}


T
Chimera_new (Stage3_T nonchimericbest, int chimerapos, int chimeraequivpos) {
  T new = (T) MALLOC(sizeof(*new));

  new->nonchimera_matches = Stage3_matches(nonchimericbest);
  new->nonchimera_mismatches = Stage3_mismatches(nonchimericbest);
  new->nonchimera_indels = Stage3_indels(nonchimericbest);
  new->chimerapos = chimerapos;
  new->equivpos = chimeraequivpos;

  return new;
}

void
Chimera_free (T *old) {
  FREE(*old);
  return;
}


void
Chimera_print (T this) {
  if (this->exonexonpos >= 0) {
    printf(" *** Possible chimera with exon-exon boundary");
    if (this->cdna_direction > 0) {
      printf(" (sense)");
    } else if (this->cdna_direction < 0) {
      printf(" (antisense)");
    }
    printf(" at %d (donor_prob = %.3f, acceptor_prob = %.3f)",
	   this->exonexonpos+1,this->donor_prob,this->acceptor_prob);
  } else if (this->equivpos == this->chimerapos) {
    printf(" *** Possible chimera with breakpoint at %d",this->chimerapos+1);
  } else {
    printf(" *** Possible chimera with breakpoint at %d..%d",this->chimerapos+1,this->equivpos+1);
  }

  printf(" -- alternative has %d matches, %d mismatches, %d indels ***",
	 this->nonchimera_matches,this->nonchimera_mismatches,this->nonchimera_indels);
  return;
}


#define NPSEUDO 10.0

int
Chimera_detect (int *margin, Stage3_T stage3, Sequence_T queryseq, double fthreshold) {
  int breakpoint;
  double pvalue = 1.0;
  int *matchscores;
  int querylength, leftmargin, rightmargin;

  /* x signifies nmatches, y signifies nmismatches, x + y = n */
  int start, end, pos, x = 0, y, n, x_left, y_left, x_right, y_right, n_left, n_right;
  double theta, x_pseudo, theta_left, theta_right, rss, rss_left, rss_right, rss_sep;
  double fscore, min_rss_sep, best_pos = -1;

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
    *margin = 0;
    FREE(matchscores);
    return -1;
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
  FREE(matchscores);

  fscore = ((double) (n - 2))*(rss - min_rss_sep)/min_rss_sep;
  if (fscore < fthreshold) {
    *margin = 0;
    return -1;
  } else {
    breakpoint = best_pos;
    if ((leftmargin = breakpoint - Sequence_trim_start(queryseq)) < 0) {
      leftmargin = 0;
    }
    if ((rightmargin = Sequence_trim_end(queryseq) - breakpoint) < 0) {
      rightmargin = 0;
    }

    /* Return smaller margin */
    if (leftmargin < rightmargin) {
      *margin = leftmargin;
    } else {
      *margin = rightmargin;
    }
    debug1(printf("at %d (margin = %d), fscore = %f\n",breakpoint,*margin,fscore));
    return breakpoint;
  }
}



void
Chimera_bestpath (int *chimerapos, int *chimeraequivpos, int *bestfrom, int *bestto, 
		  Stage3_T *stage3array, int npaths, int querylength) {
  int **matrix, *from, *to, *bestscoreatpos, i, j, pos, score, 
    bestscore;
  int fromscore_3, toscore_5;

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
      *chimerapos = *chimeraequivpos = pos;
      *bestfrom = from[pos];
      *bestto = to[pos];
    } else if (bestscoreatpos[pos] == bestscore) {
      *chimeraequivpos = pos;
    }
  }


  debug(
	fromscore_3 = matrix[*bestfrom][querylength-1] - matrix[*bestfrom][*chimerapos];
	toscore_5 = matrix[*bestto][*chimerapos];
	printf("From path %d to path %d at pos %d.  Score = %d.  Fromscore_3 = %d, Toscore_5 = %d\n",
	       *bestfrom,*bestto,*bestpos,bestscore,*fromscore_3,*toscore_5);
	);

  FREE(bestscoreatpos);
  FREE(to);
  FREE(from);

  for (i = 0; i < npaths; i++) {
    FREE(matrix[i]);
  }
  FREE(matrix);

  return;
}


static double
find_exonexon_fwd (int *exonexonpos, double *donor_prob, double *acceptor_prob,
		   Stage3_T left_part, Stage3_T right_part, Genome_T genome,
		   int breakpoint_start, int breakpoint_end) {
  Sequence_T donor_genomicseg, acceptor_genomicseg;
  char gbuffer1[512], gbuffer2[512], gbuffer3[512], gbuffer4[512], *donor_ptr, *acceptor_ptr;
  int i, j;
  Genomicpos_T left, right, donor_length, acceptor_length;
  bool revcomp;
  double donor_prob_1, acceptor_prob_1, bestproduct = 0.0, product;

  *exonexonpos = -1;

  donor_length = breakpoint_end - breakpoint_start + DONOR_MODEL_LEFT_MARGIN + DONOR_MODEL_RIGHT_MARGIN;
  left = Stage3_genomicpos(left_part,breakpoint_start,/*headp*/false);
  if (Stage3_watsonp(left_part) == true) {
    left -= DONOR_MODEL_LEFT_MARGIN;
    revcomp = false;
  } else {
    left += DONOR_MODEL_LEFT_MARGIN;
    left -= donor_length;
    revcomp = true;
  }
  donor_genomicseg = Genome_get_segment(genome,left,donor_length+1,revcomp,
					gbuffer1,gbuffer2,/*gbufferlen*/512);
  donor_ptr = Sequence_fullpointer(donor_genomicseg);
  debug2(
	 for (i = 0; i < ACCEPTOR_MODEL_LEFT_MARGIN - DONOR_MODEL_LEFT_MARGIN; i++) {
	   printf(" ");
	 }
	 Sequence_print(donor_genomicseg,/*uppercasep*/true,/*wraplength*/50,/*trimmedp*/false);
	 );

  acceptor_length = breakpoint_end - breakpoint_start + ACCEPTOR_MODEL_LEFT_MARGIN + ACCEPTOR_MODEL_RIGHT_MARGIN;
  left = Stage3_genomicpos(right_part,breakpoint_end+1,/*headp*/true);
  if (Stage3_watsonp(right_part) == true) {
    left += ACCEPTOR_MODEL_RIGHT_MARGIN;
    left -= acceptor_length;
    revcomp = false;
  } else {
    left -= ACCEPTOR_MODEL_RIGHT_MARGIN;
    revcomp = true;
  }
  acceptor_genomicseg = Genome_get_segment(genome,left,acceptor_length+1,revcomp,
					   gbuffer3,gbuffer4,/*gbufferlen*/512);
  acceptor_ptr = Sequence_fullpointer(acceptor_genomicseg);
  debug2(
	 for (i = 0; i < DONOR_MODEL_LEFT_MARGIN - ACCEPTOR_MODEL_LEFT_MARGIN; i++) {
	   printf(" ");
	 }
	 Sequence_print(acceptor_genomicseg,/*uppercasep*/true,/*wraplength*/50,/*trimmedp*/false);
	 );

  *donor_prob = 0.0;
  *acceptor_prob = 0.0;
  for (i = DONOR_MODEL_LEFT_MARGIN, j = ACCEPTOR_MODEL_LEFT_MARGIN; 
       i <= donor_length - DONOR_MODEL_RIGHT_MARGIN && 
	 j <= acceptor_length - ACCEPTOR_MODEL_RIGHT_MARGIN;
       i++, j++) {
    debug2(printf("%c%c %c%c\n",donor_ptr[i+1],donor_ptr[i+2],acceptor_ptr[j-2],acceptor_ptr[j-1]));
    if (donor_ptr[i+1] == 'G' && donor_ptr[i+2] == 'T' &&
	acceptor_ptr[j-2] == 'A' && acceptor_ptr[j-1] == 'G') {
      donor_prob_1 = Maxent_donor_prob(&(donor_ptr[i+1-DONOR_MODEL_LEFT_MARGIN]));
      acceptor_prob_1 = Maxent_acceptor_prob(&(acceptor_ptr[j-ACCEPTOR_MODEL_LEFT_MARGIN]));
      if ((product = donor_prob_1*acceptor_prob_1) > bestproduct) {
	bestproduct = product;
	*donor_prob = donor_prob_1;
	*acceptor_prob = acceptor_prob_1;
	*exonexonpos = breakpoint_start - DONOR_MODEL_LEFT_MARGIN + i;
	debug2(printf("%f %f %d\n",donor_prob_1,acceptor_prob_1,*exonexonpos));
      }
    }
  }

  Sequence_free(&acceptor_genomicseg);
  Sequence_free(&donor_genomicseg);

  return bestproduct;
}

static double
find_exonexon_rev (int *exonexonpos, double *donor_prob, double *acceptor_prob,
		   Stage3_T left_part, Stage3_T right_part, Genome_T genome,
		   int breakpoint_start, int breakpoint_end) {
  Sequence_T donor_genomicseg, acceptor_genomicseg;
  char gbuffer1[512], gbuffer2[512], gbuffer3[512], gbuffer4[512], *donor_ptr, *acceptor_ptr;
  int i, j;
  Genomicpos_T left, right, donor_length, acceptor_length;
  bool revcomp;
  double donor_prob_1, acceptor_prob_1, bestproduct = 0.0, product;

  *exonexonpos = -1;

  donor_length = breakpoint_end - breakpoint_start + DONOR_MODEL_LEFT_MARGIN + DONOR_MODEL_RIGHT_MARGIN;
  left = Stage3_genomicpos(right_part,breakpoint_end+1,/*headp*/true);
  if (Stage3_watsonp(right_part) == true) {
    left += DONOR_MODEL_LEFT_MARGIN;
    left -= donor_length;
    revcomp = true;
  } else {
    left -= DONOR_MODEL_LEFT_MARGIN;
    revcomp = false;
  }
  donor_genomicseg = Genome_get_segment(genome,left,donor_length+1,revcomp,
					gbuffer1,gbuffer2,/*gbufferlen*/512);
  donor_ptr = Sequence_fullpointer(donor_genomicseg);
  debug2(
	 for (i = 0; i < ACCEPTOR_MODEL_LEFT_MARGIN - DONOR_MODEL_LEFT_MARGIN; i++) {
	   printf(" ");
	 }
	 Sequence_print(donor_genomicseg,/*uppercasep*/true,/*wraplength*/50,/*trimmedp*/false);
	 );

  acceptor_length = breakpoint_end - breakpoint_start + ACCEPTOR_MODEL_LEFT_MARGIN + ACCEPTOR_MODEL_RIGHT_MARGIN;
  left = Stage3_genomicpos(left_part,breakpoint_start,/*headp*/false);
  if (Stage3_watsonp(left_part) == true) {
    left -= ACCEPTOR_MODEL_RIGHT_MARGIN;
    revcomp = true;
  } else {
    left += ACCEPTOR_MODEL_RIGHT_MARGIN;
    left -= acceptor_length;
    revcomp = false;
  }
  acceptor_genomicseg = Genome_get_segment(genome,left,acceptor_length+1,revcomp,
					   gbuffer3,gbuffer4,/*gbufferlen*/512);
  acceptor_ptr = Sequence_fullpointer(acceptor_genomicseg);
  debug2(
	 for (i = 0; i < DONOR_MODEL_LEFT_MARGIN - ACCEPTOR_MODEL_LEFT_MARGIN; i++) {
	   printf(" ");
	 }
	 Sequence_print(acceptor_genomicseg,/*uppercasep*/true,/*wraplength*/50,/*trimmedp*/false);
	 );

  *donor_prob = 0.0;
  *acceptor_prob = 0.0;
  for (i = DONOR_MODEL_LEFT_MARGIN, j = ACCEPTOR_MODEL_LEFT_MARGIN; 
       i <= donor_length - DONOR_MODEL_RIGHT_MARGIN && 
	 j <= acceptor_length - ACCEPTOR_MODEL_RIGHT_MARGIN;
       i++, j++) {
    debug2(printf("%c%c %c%c\n",donor_ptr[i+1],donor_ptr[i+2],acceptor_ptr[j-2],acceptor_ptr[j-1]));
    if (donor_ptr[i+1] == 'G' && donor_ptr[i+2] == 'T' &&
	acceptor_ptr[j-2] == 'A' && acceptor_ptr[j-1] == 'G') {
      donor_prob_1 = Maxent_donor_prob(&(donor_ptr[i+1-DONOR_MODEL_LEFT_MARGIN]));
      acceptor_prob_1 = Maxent_acceptor_prob(&(acceptor_ptr[j-ACCEPTOR_MODEL_LEFT_MARGIN]));
      if ((product = donor_prob_1*acceptor_prob_1) > bestproduct) {
	bestproduct = product;
	*donor_prob = donor_prob_1;
	*acceptor_prob = acceptor_prob_1;
	*exonexonpos = breakpoint_end + DONOR_MODEL_LEFT_MARGIN - i;
	debug2(printf("%f %f %d\n",donor_prob_1,acceptor_prob_1,*exonexonpos));
      }
    }
  }

  Sequence_free(&acceptor_genomicseg);
  Sequence_free(&donor_genomicseg);

  return bestproduct;
}


void
Chimera_find_exonexon (T this, Stage3_T left_part, Stage3_T right_part,
		       Genome_T genome) {
  int breakpoint_start, breakpoint_end, exonexonpos_fwd, exonexonpos_rev;
  double bestproduct_fwd, bestproduct_rev, donor_prob_fwd, donor_prob_rev, acceptor_prob_fwd, acceptor_prob_rev;
  int left_cdna_direction, right_cdna_direction, try_direction;

  breakpoint_start = this->chimerapos;
  breakpoint_end = this->equivpos;

  left_cdna_direction = Stage3_cdna_direction(left_part);
  right_cdna_direction = Stage3_cdna_direction(right_part);

  if (left_cdna_direction == 0 && right_cdna_direction == 0) {
    try_direction = 0;
  } else if (left_cdna_direction >= 0 && right_cdna_direction >= 0) {
    try_direction = +1;
  } else if (left_cdna_direction <= 0 && right_cdna_direction <= 0) {
    try_direction = -1;
  } else {
    try_direction = 0;
  }

  if (try_direction == +1) {
    this->cdna_direction = +1;
    find_exonexon_fwd(&this->exonexonpos,&this->donor_prob,&this->acceptor_prob,left_part,right_part,genome,
		      breakpoint_start,breakpoint_end);
  } else if (try_direction == -1) {
    this->cdna_direction = -1;
    find_exonexon_rev(&this->exonexonpos,&this->donor_prob,&this->acceptor_prob,left_part,right_part,genome,
		      breakpoint_start,breakpoint_end);
  } else {
    bestproduct_fwd = find_exonexon_fwd(&exonexonpos_fwd,&donor_prob_fwd,&acceptor_prob_fwd,
					left_part,right_part,genome,breakpoint_start,breakpoint_end);
    bestproduct_rev = find_exonexon_rev(&exonexonpos_rev,&donor_prob_rev,&acceptor_prob_rev,
					left_part,right_part,genome,breakpoint_start,breakpoint_end);
    if (bestproduct_fwd == 0.0 && bestproduct_rev == 0.0) {
      this->cdna_direction = 0;
      this->exonexonpos = -1;
      this->donor_prob = 0.0;
      this->acceptor_prob = 0.0;
    } else if (bestproduct_fwd >= bestproduct_rev) {
      this->cdna_direction = +1;
      this->exonexonpos = exonexonpos_fwd;
      this->donor_prob = donor_prob_fwd;
      this->acceptor_prob = acceptor_prob_fwd;
    } else {
      this->cdna_direction = -1;
      this->exonexonpos = exonexonpos_rev;
      this->donor_prob = donor_prob_rev;
      this->acceptor_prob = acceptor_prob_rev;
    }
  }

  if (this->exonexonpos >= 0) {
    this->chimerapos = this->equivpos = this->exonexonpos;
  }

  return;
}
