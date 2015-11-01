static char rcsid[] = "$Id: splice.c 101271 2013-07-12 02:44:39Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "splice.h"

#include <stdio.h>
#include "sense.h"
#include "genome_hr.h"
#include "genome_sites.h"
#include "substring.h"
#include "maxent.h"
#include "maxent_hr.h"
#include "stage3hr.h"


#define LOWPROB_SUPPORT 20


/* Splice_solve_single */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

/* Splice_solve_double */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif


static bool novelsplicingp = true;


/* Do not compare against true or false */
/* Loosest criterion */
static int
sufficient_splice_prob_local (int support, int nmismatches, double spliceprob) {
  support -= 3*nmismatches;
  if (support < 14) {
    return (spliceprob > 0.95);
  } else if (support < 20) {
    return (spliceprob > 0.90);
  } else if (support < 26) {
    return (spliceprob > 0.85);
  } else {
    return (spliceprob > 0.70);
  }
}



/* Note: knowni holds joffset + j + 1, so 0 represents no known site
   and values greater than 0 represent a known site.  Need to subtract
   1 to obtain joffset + j. */

List_T
Splice_solve_single (int *found_score, int *nhits, List_T hits, List_T *lowprob,

		     bool *segmenti_usedp, bool *segmentj_usedp,
		     Univcoord_T segmenti_left, Univcoord_T segmentj_left,
		     Chrnum_T segmenti_chrnum, Univcoord_T segmenti_chroffset,
		     Univcoord_T segmenti_chrhigh, Chrpos_T segmenti_chrlength,
		     Chrnum_T segmentj_chrnum, Univcoord_T segmentj_chroffset,
		     Univcoord_T segmentj_chrhigh, Chrpos_T segmentj_chrlength,
		     
		     int querylength, Compress_T query_compress,
		     int *segmenti_donor_knownpos, int *segmentj_acceptor_knownpos,
		     int *segmentj_antidonor_knownpos, int *segmenti_antiacceptor_knownpos,
		     int *segmenti_donor_knowni, int *segmentj_acceptor_knowni,
		     int *segmentj_antidonor_knowni, int *segmenti_antiacceptor_knowni,
		     int segmenti_donor_nknown, int segmentj_acceptor_nknown,
		     int segmentj_antidonor_nknown, int segmenti_antiacceptor_nknown,
		     int splicing_penalty, int max_mismatches_allowed,
		     bool first_read_p, bool plusp, int genestrand, bool subs_or_indels_p,
		     bool sarrayp) {
  Substring_T donor, acceptor;
  int best_splice_pos, splice_pos_start, splice_pos_end, splice_pos, i, j;
  int donor_positions_alloc[MAX_READLENGTH+1], acceptor_positions_alloc[MAX_READLENGTH+1];
  int donor_knowni_alloc[MAX_READLENGTH+1], acceptor_knowni_alloc[MAX_READLENGTH+1];

  int best_nmismatches, nmismatches;
  int best_segmenti_nmismatches, best_segmentj_nmismatches, segmenti_nmismatches, segmentj_nmismatches;
  int donor_support, acceptor_support;
  int best_donor_knowni, best_acceptor_knowni;
  double best_prob, best_donor_prob, best_acceptor_prob, probi, probj;
  bool sufficient1p, sufficient2p, orig_plusp, sensep;
  int sensedir;

  int donori_nsites, acceptorj_nsites, antiacceptori_nsites, antidonorj_nsites;
  int *donori_positions, *acceptorj_positions, *antiacceptori_positions, *antidonorj_positions;
  int *donori_knowni, *acceptorj_knowni, *antiacceptori_knowni, *antidonorj_knowni;


  debug1(printf("Splice_solve_single: Getting genome at lefti %u and leftj %u (diff: %d)\n",
		 segmenti_left,segmentj_left,segmentj_left-segmenti_left));

#if 0
  int sum, lefti, righti;
  splice_pos_start = querylength;
  splice_pos_end = 0;
  for (sum = 0; sum <= max_mismatches_allowed; sum++) {
    for (lefti = 0; lefti <= sum && lefti < nmismatches_left; lefti++) {
      if ((righti = sum - lefti) < nmismatches_right &&
	  mismatch_positions_left[lefti] > mismatch_positions_right[righti]) {
	debug1(printf("At %d+%d mismatches, splice_pos using right: %d\n",lefti,righti,mismatch_positions_right[righti]+1));
	debug1(printf("At %d+%d mismatches, splice_pos using left: %d\n",lefti,righti,mismatch_positions_left[lefti]));
	if (mismatch_positions_right[righti] + 1 < splice_pos_start) {
	  splice_pos_start = mismatch_positions_right[righti] + 1;	/* This is leftmost position in righti+1 .. lefti */
	}
	if (mismatch_positions_left[lefti] > splice_pos_end) {
	  splice_pos_end = mismatch_positions_left[lefti];	/* This is rightmost position in righti+1 .. lefti */
	}
      }
    }
  }

  /* Exclude ends */
  if (splice_pos_start < min_localsplicing_end_matches) {
    splice_pos_start = min_localsplicing_end_matches;
  }
  if (splice_pos_end > querylength - min_localsplicing_end_matches) {
    splice_pos_end = querylength - min_localsplicing_end_matches;
  }
#else
  /* splice_pos_start = min_localsplicing_end_matches; */
  /* splice_pos_end = querylength - min_localsplicing_end_matches; */
  splice_pos_start = 2;
  splice_pos_end = querylength - 2;
#endif


  if (splice_pos_start <= splice_pos_end) {
    /* Originally from plus strand.  No complement.  */
    /* Sense (End 1 to End 2) or Antisense (End 5 to End 6) */
    if (novelsplicingp && segmenti_left + splice_pos_start >= DONOR_MODEL_LEFT_MARGIN) {
      donori_nsites = Genome_donor_positions(donor_positions_alloc,donor_knowni_alloc,
					     segmenti_donor_knownpos,segmenti_donor_knowni,
					     segmenti_left,splice_pos_start,splice_pos_end);
      donori_positions = donor_positions_alloc;
      donori_knowni = donor_knowni_alloc;
    } else {
      donori_nsites = segmenti_donor_nknown;
      donori_positions = segmenti_donor_knownpos;
      donori_knowni = segmenti_donor_knowni;
    }

#ifdef DEBUG1
    printf("Found %d donori sites:",donori_nsites);
    for (i = 0; i < donori_nsites; i++) {
      printf(" %d",donori_positions[i]);
      if (donori_knowni[i] >= 0) {
	printf(" (%d)",donori_knowni[i]);
      }
    }
    printf("\n");
#endif

    if (novelsplicingp && segmentj_left + splice_pos_start >= ACCEPTOR_MODEL_LEFT_MARGIN) {
      acceptorj_nsites = Genome_acceptor_positions(acceptor_positions_alloc,acceptor_knowni_alloc,
						   segmentj_acceptor_knownpos,segmentj_acceptor_knowni,
						   segmentj_left,splice_pos_start,splice_pos_end);
      acceptorj_positions = acceptor_positions_alloc;
      acceptorj_knowni = acceptor_knowni_alloc;
    } else {
      acceptorj_nsites = segmentj_acceptor_nknown;
      acceptorj_positions = segmentj_acceptor_knownpos;
      acceptorj_knowni = segmentj_acceptor_knowni;
    }

#ifdef DEBUG1
    printf("Found %d acceptorj sites:",acceptorj_nsites);
    for (i = 0; i < acceptorj_nsites; i++) {
      printf(" %d",acceptorj_positions[i]);
      if (acceptorj_knowni[i] >= 0) {
	printf(" (%d)",acceptorj_knowni[i]);
      }
    }
    printf("\n");
#endif

    best_nmismatches = max_mismatches_allowed;
    best_prob = 0.0;
    orig_plusp = true;

    i = j = 0;
    while (i < donori_nsites && j < acceptorj_nsites) {
      if ((splice_pos = donori_positions[i]) < acceptorj_positions[j]) {
	i++;
      } else if (splice_pos > acceptorj_positions[j]) {
	j++;
      } else {
	segmenti_nmismatches = Genome_count_mismatches_substring(query_compress,/*left*/segmenti_left,/*pos5*/0,/*pos3*/splice_pos,
								 plusp,genestrand);
	segmentj_nmismatches = Genome_count_mismatches_substring(query_compress,/*left*/segmentj_left,/*pos5*/splice_pos,/*pos3*/querylength,
								 plusp,genestrand);
	if ((nmismatches = segmenti_nmismatches + segmentj_nmismatches) <= best_nmismatches) {
	  if (donori_knowni[i] >= 0) {
	    probi = 1.0;
	  } else {
	    probi = Maxent_hr_donor_prob(segmenti_left + splice_pos,segmenti_chroffset);
	  }

	  if (acceptorj_knowni[j] >= 0) {
	    probj = 1.0;
	  } else {
	    probj = Maxent_hr_acceptor_prob(segmentj_left + splice_pos,segmentj_chroffset);
	  }

	  debug1(
		  if (plusp == true) {
		    printf("plus sense splice_pos  %d, i.donor %f, j.acceptor %f\n",splice_pos,probi,probj);
		  } else {
		    printf("minus antisense splice_pos  %d, i.donor %f, j.acceptor %f\n",splice_pos,probi,probj);
		  });

	  if (nmismatches < best_nmismatches ||
	      (nmismatches == best_nmismatches && probi + probj > best_prob)) {
	    /* Success */
	    best_nmismatches = nmismatches;
	    best_prob = probi + probj;

	    best_donor_knowni = donori_knowni[i];
	    best_acceptor_knowni = acceptorj_knowni[j];
	    best_donor_prob = probi;
	    best_acceptor_prob = probj;
	    best_splice_pos = splice_pos;
	    best_segmenti_nmismatches = segmenti_nmismatches;
	    best_segmentj_nmismatches = segmentj_nmismatches;
	  }
	}
	i++;
	j++;
      }
    }


    /* Originally from minus strand.  Complement. */
    /* Antisense (End 7 to End 8) or Sense (End 3 to End 4) */
    if (novelsplicingp && segmenti_left + splice_pos_start >= ACCEPTOR_MODEL_RIGHT_MARGIN) {
      antiacceptori_nsites = Genome_antiacceptor_positions(acceptor_positions_alloc,acceptor_knowni_alloc,
							   segmenti_antiacceptor_knownpos,segmenti_antiacceptor_knowni,
							   segmenti_left,splice_pos_start,splice_pos_end);
      antiacceptori_positions = acceptor_positions_alloc;
      antiacceptori_knowni = acceptor_knowni_alloc;
    } else {
      antiacceptori_nsites = segmenti_antiacceptor_nknown;
      antiacceptori_positions = segmenti_antiacceptor_knownpos;
      antiacceptori_knowni = segmenti_antiacceptor_knowni;
    }

#ifdef DEBUG1
    printf("Found %d antiacceptori sites:",antiacceptori_nsites);
    for (i = 0; i < antiacceptori_nsites; i++) {
      printf(" %d",antiacceptori_positions[i]);
      if (antiacceptori_knowni[i] >= 0) {
	printf(" (%d)",antiacceptori_knowni[i]);
      }
    }
    printf("\n");
#endif

    if (novelsplicingp && segmentj_left + splice_pos_start >= DONOR_MODEL_RIGHT_MARGIN) {
      antidonorj_nsites = Genome_antidonor_positions(donor_positions_alloc,donor_knowni_alloc,
						     segmentj_antidonor_knownpos,segmentj_antidonor_knowni,
						     segmentj_left,splice_pos_start,splice_pos_end);
      antidonorj_positions = donor_positions_alloc;
      antidonorj_knowni = donor_knowni_alloc;
    } else {
      antidonorj_nsites = segmentj_antidonor_nknown;
      antidonorj_positions = segmentj_antidonor_knownpos;
      antidonorj_knowni = segmentj_antidonor_knowni;
    }

#ifdef DEBUG1
    printf("Found %d antidonorj sites:",antidonorj_nsites);
    for (i = 0; i < antidonorj_nsites; i++) {
      printf(" %d",antidonorj_positions[i]);
      if (antidonorj_knowni[i] >= 0) {
	printf(" (%d)",antidonorj_knowni[i]);
      }
    }
    printf("\n");
#endif

    i = j = 0;
    while (i < antiacceptori_nsites && j < antidonorj_nsites) {
      if ((splice_pos = antiacceptori_positions[i]) < antidonorj_positions[j]) {
	i++;
      } else if (splice_pos > antidonorj_positions[j]) {
	j++;
      } else {
	segmenti_nmismatches = Genome_count_mismatches_substring(query_compress,/*left*/segmenti_left,/*pos5*/0,/*pos3*/splice_pos,
								 plusp,genestrand);
	segmentj_nmismatches = Genome_count_mismatches_substring(query_compress,/*left*/segmentj_left,/*pos5*/splice_pos,/*pos3*/querylength,
								 plusp,genestrand);
	if ((nmismatches = segmenti_nmismatches + segmentj_nmismatches) <= best_nmismatches) {
	  if (antiacceptori_knowni[i] >= 0) {
	    probi = 1.0;
	  } else {
	    probi = Maxent_hr_antiacceptor_prob(segmenti_left + splice_pos,segmenti_chroffset);
	  }

	  if (antidonorj_knowni[j] >= 0) {
	    probj = 1.0;
	  } else {
	    probj = Maxent_hr_antidonor_prob(segmentj_left + splice_pos,segmentj_chroffset);
	  }

	  debug1(
		  if (plusp == true) {
		    printf("plus antisense splice_pos  %d, j.donor %f, i.acceptor %f\n",splice_pos,probj,probi);
		  } else {
		    printf("minus sense splice_pos  %d, j.donor %f, i.acceptor %f\n",splice_pos,probj,probi);
		  });
	  
	  if (nmismatches < best_nmismatches ||
	      (nmismatches == best_nmismatches && probi + probj > best_prob)) {
	    /* Success */
	    best_nmismatches = nmismatches;
	    best_prob = probi + probj;

	    best_donor_knowni = antidonorj_knowni[j];
	    best_acceptor_knowni = antiacceptori_knowni[i];
	    best_donor_prob = probj;
	    best_acceptor_prob = probi;
	    best_splice_pos = splice_pos;
	    best_segmentj_nmismatches = segmentj_nmismatches;
	    best_segmenti_nmismatches = segmenti_nmismatches;
	    orig_plusp = false;
	  }
	}
	i++;
	j++;
      }
    }

    if (best_prob > 0.0) {
      debug1(printf("best_prob = %f at splice_pos %d (%d,%d)\n",
		     best_prob,best_splice_pos,best_donor_knowni,best_acceptor_knowni));
      if (orig_plusp == true) {
	/* Originally from plus strand.  No complement. */
	sensep = (plusp == true) ? true : false;
	sensedir = (plusp == true) ? SENSE_FORWARD : SENSE_ANTI;

	donor = Substring_new_donor(best_donor_knowni,/*joffset*/0,best_splice_pos,best_segmenti_nmismatches,
				    best_donor_prob,/*left*/segmenti_left,query_compress,
				    querylength,plusp,genestrand,sensep,
				    segmenti_chrnum,segmenti_chroffset,segmenti_chrhigh,segmenti_chrlength);

	acceptor = Substring_new_acceptor(best_acceptor_knowni,/*joffset*/0,best_splice_pos,best_segmentj_nmismatches,
					  best_acceptor_prob,/*left*/segmentj_left,query_compress,
					  querylength,plusp,genestrand,sensep,
					  segmentj_chrnum,segmentj_chroffset,segmentj_chrhigh,segmentj_chrlength);

	if (donor == NULL || acceptor == NULL) {
	  if (donor != NULL) Substring_free(&donor);
	  if (acceptor != NULL) Substring_free(&acceptor);
	} else {
	  *nhits += 1;
	  debug1(printf("Splice_solve_single success\n"));
	  *segmenti_usedp = *segmentj_usedp = true;

	  donor_support = best_splice_pos;
	  acceptor_support = querylength - best_splice_pos;
	  sufficient1p = sufficient_splice_prob_local(donor_support,best_segmenti_nmismatches,best_donor_prob);
	  sufficient2p = sufficient_splice_prob_local(acceptor_support,best_segmentj_nmismatches,best_acceptor_prob);

	  if (sufficient1p && sufficient2p) {
	    return List_push(hits,(void *) Stage3end_new_splice(&(*found_score),best_segmenti_nmismatches,best_segmentj_nmismatches,
								donor,acceptor,/*distance*/segmentj_left - segmenti_left,
								/*shortdistancep*/true,splicing_penalty,querylength,
								/*amb_nmatches*/0,/*ambi_left*/NULL,/*ambi_right*/NULL,
								/*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
								/*copy_donor_p*/false,/*copy_acceptor_p*/false,first_read_p,sensedir,
								sarrayp));
	  } else if (subs_or_indels_p == true) {
	    if (donor != NULL) Substring_free(&donor);
	    if (acceptor != NULL) Substring_free(&acceptor);
	    return hits;
	  } else if (donor_support < LOWPROB_SUPPORT || acceptor_support < LOWPROB_SUPPORT) {
	    if (donor != NULL) Substring_free(&donor);
	    if (acceptor != NULL) Substring_free(&acceptor);
	    return hits;
	  } else if (sufficient1p || sufficient2p) {
	    *lowprob = List_push(*lowprob,
				 (void *) Stage3end_new_splice(&(*found_score),best_segmenti_nmismatches,best_segmentj_nmismatches,
							       donor,acceptor,/*distance*/segmentj_left - segmenti_left,
							       /*shortdistancep*/true,splicing_penalty,querylength,
							       /*amb_nmatches*/0,/*ambi_left*/NULL,/*ambi_right*/NULL,
							       /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
							       /*copy_donor_p*/false,/*copy_acceptor_p*/false,first_read_p,sensedir,
							       sarrayp));
	    return hits;
	  } else {
	    if (donor != NULL) Substring_free(&donor);
	    if (acceptor != NULL) Substring_free(&acceptor);
	  }
	}

      } else {
	/* Originally from minus strand.  Complement. */
	sensep = (plusp == true) ? false : true;
	sensedir = (plusp == true) ? SENSE_ANTI : SENSE_FORWARD;

	donor = Substring_new_donor(best_donor_knowni,/*joffset*/0,best_splice_pos,best_segmentj_nmismatches,
				    best_donor_prob,/*left*/segmentj_left,query_compress,
				    querylength,plusp,genestrand,sensep,
				    segmentj_chrnum,segmentj_chroffset,segmentj_chrhigh,segmentj_chrlength);

	acceptor = Substring_new_acceptor(best_acceptor_knowni,/*joffset*/0,best_splice_pos,best_segmenti_nmismatches,
					  best_acceptor_prob,/*left*/segmenti_left,query_compress,
					  querylength,plusp,genestrand,sensep,
					  segmenti_chrnum,segmenti_chroffset,segmenti_chrhigh,segmenti_chrlength);

	if (donor == NULL || acceptor == NULL) {
	  if (donor != NULL) Substring_free(&donor);
	  if (acceptor != NULL) Substring_free(&acceptor);
	} else {
	  *nhits += 1;
	  debug1(printf("Splice_solve_single success\n"));
	  *segmenti_usedp = *segmentj_usedp = true;

	  acceptor_support = best_splice_pos;
	  donor_support = querylength - best_splice_pos;
	  sufficient1p = sufficient_splice_prob_local(acceptor_support,best_segmenti_nmismatches,best_acceptor_prob);
	  sufficient2p = sufficient_splice_prob_local(donor_support,best_segmentj_nmismatches,best_donor_prob);
	  if (sufficient1p && sufficient2p) {
	    return List_push(hits,(void *) Stage3end_new_splice(&(*found_score),best_segmentj_nmismatches,best_segmenti_nmismatches,
								donor,acceptor,/*distance*/segmentj_left - segmenti_left,
								/*shortdistancep*/true,splicing_penalty,querylength,
								/*amb_nmatches*/0,/*ambi_left*/NULL,/*ambi_right*/NULL,
								/*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
								/*copy_donor_p*/false,/*copy_acceptor_p*/false,first_read_p,sensedir,
								sarrayp));
	  } else if (subs_or_indels_p == true) {
	    if (donor != NULL) Substring_free(&donor);
	    if (acceptor != NULL) Substring_free(&acceptor);
	    return hits;
	  } else if (donor_support < LOWPROB_SUPPORT || acceptor_support < LOWPROB_SUPPORT) {
	    if (donor != NULL) Substring_free(&donor);
	    if (acceptor != NULL) Substring_free(&acceptor);
	    return hits;
	  } else if (sufficient1p || sufficient2p) {
	    *lowprob = List_push(*lowprob,
				 (void *) Stage3end_new_splice(&(*found_score),best_segmentj_nmismatches,best_segmenti_nmismatches,
							       donor,acceptor,/*distance*/segmentj_left - segmenti_left,
							       /*shortdistancep*/true,splicing_penalty,querylength,
							       /*amb_nmatches*/0,/*ambi_left*/NULL,/*ambi_right*/NULL,
							       /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
							       /*copy_donor_p*/false,/*copy_acceptor_p*/false,first_read_p,sensedir,
							       sarrayp));
	    return hits;
	  } else {
	    if (donor != NULL) Substring_free(&donor);
	    if (acceptor != NULL) Substring_free(&acceptor);
	    return hits;
	  }
	}
      }
    }
  }

  debug1(printf("Splice_solve_single fail\n"));
  return hits;
}


List_T
Splice_solve_double (int *found_score, int *nhits, List_T hits, List_T *lowprob,

		     bool *segmenti_usedp, bool *segmentm_usedp, bool *segmentj_usedp,
		     Univcoord_T segmenti_left, Univcoord_T segmentm_left, Univcoord_T segmentj_left,
		     Chrnum_T segmenti_chrnum, Univcoord_T segmenti_chroffset,
		     Univcoord_T segmenti_chrhigh, Chrpos_T segmenti_chrlength,
		     Chrnum_T segmentm_chrnum, Univcoord_T segmentm_chroffset,
		     Univcoord_T segmentm_chrhigh, Chrpos_T segmentm_chrlength,
		     Chrnum_T segmentj_chrnum, Univcoord_T segmentj_chroffset,
		     Univcoord_T segmentj_chrhigh, Chrpos_T segmentj_chrlength,

		     int querylength, Compress_T query_compress,
		     int *segmenti_donor_knownpos, int *segmentm_acceptor_knownpos, int *segmentm_donor_knownpos, int *segmentj_acceptor_knownpos,
		     int *segmentj_antidonor_knownpos, int *segmentm_antiacceptor_knownpos, int *segmentm_antidonor_knownpos, int *segmenti_antiacceptor_knownpos,
		     int *segmenti_donor_knowni, int *segmentm_acceptor_knowni, int *segmentm_donor_knowni, int *segmentj_acceptor_knowni,
		     int *segmentj_antidonor_knowni, int *segmentm_antiacceptor_knowni, int *segmentm_antidonor_knowni, int *segmenti_antiacceptor_knowni,
		     int segmenti_donor_nknown, int segmentm_acceptor_nknown, int segmentm_donor_nknown, int segmentj_acceptor_nknown,
		     int segmentj_antidonor_nknown, int segmentm_antiacceptor_nknown, int segmentm_antidonor_nknown, int segmenti_antiacceptor_nknown,
		     int splicing_penalty, int max_mismatches_allowed, bool plusp, int genestrand, bool subs_or_indels_p, bool sarrayp) {
  Substring_T donor, shortexon, acceptor;
  int best_splice_pos_1, best_splice_pos_2, splice_pos_start, splice_pos_end, splice_pos_1, splice_pos_2;
  int i, a, b, j;
  int donor1_positions_alloc[MAX_READLENGTH+1], acceptor1_positions_alloc[MAX_READLENGTH+1],
    donor2_positions_alloc[MAX_READLENGTH+1], acceptor2_positions_alloc[MAX_READLENGTH+1];
  int donor1_knowni_alloc[MAX_READLENGTH+1], acceptor1_knowni_alloc[MAX_READLENGTH+1],
    donor2_knowni_alloc[MAX_READLENGTH+1], acceptor2_knowni_alloc[MAX_READLENGTH+1];

  int best_nmismatches, nmismatches;
  int best_segmenti_nmismatches, best_segmentm_nmismatches, best_segmentj_nmismatches,
    segmenti_nmismatches, segmentm_nmismatches, segmentj_nmismatches;
  int donor_support, acceptor_support, middle_support;
  int best_donor1_knowni, best_acceptor1_knowni, best_donor2_knowni, best_acceptor2_knowni;
  double best_prob, best_donor1_prob, best_acceptor1_prob, best_donor2_prob, best_acceptor2_prob,
    probi, proba, probb, probj;
  bool sufficient1p, sufficient2p, sufficient3p, sufficient4p, orig_plusp, sensep, matchp;
  int sensedir;

  int donori_nsites, acceptora_nsites, donorb_nsites, acceptorj_nsites,
    antiacceptori_nsites, antidonora_nsites, antiacceptorb_nsites, antidonorj_nsites;
  int *donori_positions, *acceptora_positions, *donorb_positions, *acceptorj_positions,
    *antiacceptori_positions, *antidonora_positions, *antiacceptorb_positions, *antidonorj_positions;
  int *donori_knowni, *acceptora_knowni, *donorb_knowni, *acceptorj_knowni,
    *antiacceptori_knowni, *antidonora_knowni, *antiacceptorb_knowni, *antidonorj_knowni;


  debug2(printf("Splice_solve_double: Getting genome at lefti %u, leftm %u, and leftj %u\n",
		 segmenti_left,segmentm_left,segmentj_left));

  splice_pos_start = 2;
  splice_pos_end = querylength - 2;

  if (splice_pos_start <= splice_pos_end) {
    /* Originally from plus strand.  No complement. */
    /* Sense (End 1 to End 2) or Antisense (End 5 to End 6) */

    /* Segment i */
    if (novelsplicingp && segmenti_left + splice_pos_start >= DONOR_MODEL_LEFT_MARGIN) {
      donori_nsites = Genome_donor_positions(donor1_positions_alloc,donor1_knowni_alloc,
					     segmenti_donor_knownpos,segmenti_donor_knowni,
					     segmenti_left,splice_pos_start,splice_pos_end);
      donori_positions = donor1_positions_alloc;
      donori_knowni = donor1_knowni_alloc;
    } else {
      donori_nsites = segmenti_donor_nknown;
      donori_positions = segmenti_donor_knownpos;
      donori_knowni = segmenti_donor_knowni;
    }

#ifdef DEBUG2
    printf("Found %d donori sites:",donori_nsites);
    for (i = 0; i < donori_nsites; i++) {
      printf(" %d",donori_positions[i]);
      if (donori_knowni[i] >= 0) {
	printf(" (%d)",donori_knowni[i]);
      }
    }
    printf("\n");
#endif

    /* Segment m1 */
    if (novelsplicingp && segmentm_left + splice_pos_start >= ACCEPTOR_MODEL_LEFT_MARGIN) {
      acceptora_nsites = Genome_acceptor_positions(acceptor1_positions_alloc,acceptor1_knowni_alloc,
						   segmentm_acceptor_knownpos,segmentm_acceptor_knowni,
						   segmentm_left,splice_pos_start,splice_pos_end);
      acceptora_positions = acceptor1_positions_alloc;
      acceptora_knowni = acceptor1_knowni_alloc;
    } else {
      acceptora_nsites = segmentm_acceptor_nknown;
      acceptora_positions = segmentm_acceptor_knownpos;
      acceptora_knowni = segmentm_acceptor_knowni;
    }

#ifdef DEBUG2
    printf("Found %d acceptora sites:",acceptora_nsites);
    for (i = 0; i < acceptora_nsites; i++) {
      printf(" %d",acceptora_positions[i]);
      if (acceptora_knowni[i] >= 0) {
	printf(" (%d)",acceptora_knowni[i]);
      }
    }
    printf("\n");
#endif

    /* Segment m2 */
    if (novelsplicingp && segmentm_left + splice_pos_start >= DONOR_MODEL_LEFT_MARGIN) {
      donorb_nsites = Genome_donor_positions(donor2_positions_alloc,donor2_knowni_alloc,
					     segmentm_donor_knownpos,segmentm_donor_knowni,
					     segmentm_left,splice_pos_start,splice_pos_end);
      donorb_positions = donor2_positions_alloc;
      donorb_knowni = donor2_knowni_alloc;
    } else {
      donorb_nsites = segmentm_donor_nknown;
      donorb_positions = segmentm_donor_knownpos;
      donorb_knowni = segmentm_donor_knowni;
    }

#ifdef DEBUG2
    printf("Found %d donorb sites:",donorb_nsites);
    for (i = 0; i < donorb_nsites; i++) {
      printf(" %d",donorb_positions[i]);
      if (donorb_knowni[i] >= 0) {
	printf(" (%d)",donorb_knowni[i]);
      }
    }
    printf("\n");
#endif

    /* Segment j */
    if (novelsplicingp && segmentj_left + splice_pos_start >= ACCEPTOR_MODEL_LEFT_MARGIN) {
      acceptorj_nsites = Genome_acceptor_positions(acceptor2_positions_alloc,acceptor2_knowni_alloc,
						   segmentj_acceptor_knownpos,segmentj_acceptor_knowni,
						   segmentj_left,splice_pos_start,splice_pos_end);
      acceptorj_positions = acceptor2_positions_alloc;
      acceptorj_knowni = acceptor2_knowni_alloc;
    } else {
      acceptorj_nsites = segmentj_acceptor_nknown;
      acceptorj_positions = segmentj_acceptor_knownpos;
      acceptorj_knowni = segmentj_acceptor_knowni;
    }

#ifdef DEBUG2
    printf("Found %d acceptorj sites:",acceptorj_nsites);
    for (i = 0; i < acceptorj_nsites; i++) {
      printf(" %d",acceptorj_positions[i]);
      if (acceptorj_knowni[i] >= 0) {
	printf(" (%d)",acceptorj_knowni[i]);
      }
    }
    printf("\n");
#endif

    best_nmismatches = max_mismatches_allowed;
    best_prob = 0.0;
    orig_plusp = true;

    i = a = b = j = 0;
    while (i < donori_nsites && a < acceptora_nsites) {
      if ((splice_pos_1 = donori_positions[i]) < acceptora_positions[a]) {
	i++;
      } else if (splice_pos_1 > acceptora_positions[a]) {
	a++;
      } else {
	while (b < donorb_nsites && donorb_positions[b] <= splice_pos_1) {
	  b++;
	}
	while (j < acceptorj_nsites && acceptorj_positions[j] <= splice_pos_1) {
	  j++;
	}
	matchp = false;
	while (b < donorb_nsites && j < acceptorj_nsites && matchp == false) {
	  if ((splice_pos_2 = donorb_positions[b]) < acceptorj_positions[j]) {
	    b++;
	  } else if (splice_pos_2 > acceptorj_positions[j]) {
	    j++;
	  } else {
	    segmenti_nmismatches = Genome_count_mismatches_substring(query_compress,/*left*/segmenti_left,/*pos5*/0,/*pos3*/splice_pos_1,
								     plusp,genestrand);
	    segmentm_nmismatches = Genome_count_mismatches_substring(query_compress,/*left*/segmentm_left,/*pos5*/splice_pos_1,/*pos3*/splice_pos_2,
								     plusp,genestrand);
	    segmentj_nmismatches = Genome_count_mismatches_substring(query_compress,/*left*/segmentj_left,/*pos5*/splice_pos_2,/*pos3*/querylength,
								     plusp,genestrand);
	    if ((nmismatches = segmenti_nmismatches + segmentm_nmismatches + segmentj_nmismatches) <= best_nmismatches) {
	      if (donori_knowni[i] >= 0) {
		probi = 1.0;
	      } else {
		probi = Maxent_hr_donor_prob(segmenti_left + splice_pos_1,segmenti_chroffset);
	      }

	      if (acceptora_knowni[a] >= 0) {
		proba = 1.0;
	      } else {
		proba = Maxent_hr_acceptor_prob(segmentm_left + splice_pos_1,segmentm_chroffset);
	      }

	      if (donorb_knowni[b] >= 0) {
		probb = 1.0;
	      } else {
		probb = Maxent_hr_donor_prob(segmentm_left + splice_pos_2,segmentm_chroffset);
	      }
	      
	      if (acceptorj_knowni[j] >= 0) {
		probj = 1.0;
	      } else {
		probj = Maxent_hr_acceptor_prob(segmentj_left + splice_pos_2,segmentj_chroffset);
	      }

	      debug2(
		      if (plusp == true) {
			printf("plus sense splice_pos  %d, %d, i.donor %f, m.acceptor %f, m.donor %f, j.acceptor %f\n",
			       splice_pos_1,splice_pos_2,probi,proba,probb,probj);
		      } else {
			printf("minus antisense splice_pos  %d %d, i.donor %f, m.acceptor %f, m.donor %f, j.acceptor %f\n",
			       splice_pos_1,splice_pos_2,probi,proba,probb,probj);
		      });

	      if (nmismatches < best_nmismatches ||
		  (nmismatches == best_nmismatches && probi + proba + probb + probj > best_prob)) {
		/* Success */
		best_nmismatches = nmismatches;
		best_prob = probi + proba + probb + probj;

		best_donor1_knowni = donori_knowni[i];
		best_acceptor1_knowni = acceptora_knowni[a];
		best_donor2_knowni = donorb_knowni[b];
		best_acceptor2_knowni = acceptorj_knowni[j];
		best_donor1_prob = probi;
		best_acceptor1_prob = proba;
		best_donor2_prob = probb;
		best_acceptor2_prob = probj;
		best_splice_pos_1 = splice_pos_1;
		best_splice_pos_2 = splice_pos_2;
		best_segmenti_nmismatches = segmenti_nmismatches;
		best_segmentm_nmismatches = segmentm_nmismatches;
		best_segmentj_nmismatches = segmentj_nmismatches;
	      }
	    }
	    /* b++; j++; Don't advance b or j, so next i/a can match */
	    matchp = true;
	  }
	}
	i++;
	a++;
      }
    }


    /* Originally from minus strand.  Complement. */
    /* Antisense (End 7 to End 8) or Sense (End 3 to End 4) */

    /* Segment i */
    if (novelsplicingp && segmenti_left + splice_pos_start >= ACCEPTOR_MODEL_RIGHT_MARGIN) {
      antiacceptori_nsites = Genome_antiacceptor_positions(acceptor1_positions_alloc,acceptor1_knowni_alloc,
							   segmenti_antiacceptor_knownpos,segmenti_antiacceptor_knowni,
							   segmenti_left,splice_pos_start,splice_pos_end);
      antiacceptori_positions = acceptor1_positions_alloc;
      antiacceptori_knowni = acceptor1_knowni_alloc;
    } else {
      antiacceptori_nsites = segmenti_antiacceptor_nknown;
      antiacceptori_positions = segmenti_antiacceptor_knownpos;
      antiacceptori_knowni = segmenti_antiacceptor_knowni;
    }

#ifdef DEBUG2
    printf("Found %d antiacceptori sites:",antiacceptori_nsites);
    for (i = 0; i < antiacceptori_nsites; i++) {
      printf(" %d",antiacceptori_positions[i]);
      if (antiacceptori_knowni[i] >= 0) {
	printf(" (%d)",antiacceptori_knowni[i]);
      }
    }
    printf("\n");
#endif

    /* Segment m1 */
    if (novelsplicingp && segmentm_left + splice_pos_start >= DONOR_MODEL_RIGHT_MARGIN) {
      antidonora_nsites = Genome_antidonor_positions(donor1_positions_alloc,donor1_knowni_alloc,
						     segmentm_antidonor_knownpos,segmentm_antidonor_knowni,
						     segmentm_left,splice_pos_start,splice_pos_end);
      antidonora_positions = donor1_positions_alloc;
      antidonora_knowni = donor1_knowni_alloc;
    } else {
      antidonora_nsites = segmentm_antidonor_nknown;
      antidonora_positions = segmentm_antidonor_knownpos;
      antidonora_knowni = segmentm_antidonor_knowni;
    }

#ifdef DEBUG2
    printf("Found %d antidonora sites:",antidonora_nsites);
    for (i = 0; i < antidonora_nsites; i++) {
      printf(" %d",antidonora_positions[i]);
      if (antidonora_knowni[i] >= 0) {
	printf(" (%d)",antidonora_knowni[i]);
      }
    }
    printf("\n");
#endif

    /* Segment m2 */
    if (novelsplicingp && segmentm_left + splice_pos_start >= ACCEPTOR_MODEL_RIGHT_MARGIN) {
      antiacceptorb_nsites = Genome_antiacceptor_positions(acceptor2_positions_alloc,acceptor2_knowni_alloc,
							   segmentm_antiacceptor_knownpos,segmentm_antiacceptor_knowni,
							   segmentm_left,splice_pos_start,splice_pos_end);
      antiacceptorb_positions = acceptor2_positions_alloc;
      antiacceptorb_knowni = acceptor2_knowni_alloc;
    } else {
      antiacceptorb_nsites = segmentm_antiacceptor_nknown;
      antiacceptorb_positions = segmentm_antiacceptor_knownpos;
      antiacceptorb_knowni = segmentm_antiacceptor_knowni;
    }

#ifdef DEBUG2
    printf("Found %d antiacceptorb sites:",antiacceptorb_nsites);
    for (i = 0; i < antiacceptorb_nsites; i++) {
      printf(" %d",antiacceptorb_positions[i]);
      if (antiacceptorb_knowni[i] >= 0) {
	printf(" (%d)",antiacceptorb_knowni[i]);
      }
    }
    printf("\n");
#endif

    /* Segment j */
    if (novelsplicingp && segmentj_left + splice_pos_start >= DONOR_MODEL_RIGHT_MARGIN) {
      antidonorj_nsites = Genome_antidonor_positions(donor2_positions_alloc,donor2_knowni_alloc,
						     segmentj_antidonor_knownpos,segmentj_antidonor_knowni,
						     segmentj_left,splice_pos_start,splice_pos_end);
      antidonorj_positions = donor2_positions_alloc;
      antidonorj_knowni = donor2_knowni_alloc;
    } else {
      antidonorj_nsites = segmentj_antidonor_nknown;
      antidonorj_positions = segmentj_antidonor_knownpos;
      antidonorj_knowni = segmentj_antidonor_knowni;
    }

#ifdef DEBUG2
    printf("Found %d antidonorj sites:",antidonorj_nsites);
    for (i = 0; i < antidonorj_nsites; i++) {
      printf(" %d",antidonorj_positions[i]);
      if (antidonorj_knowni[i] >= 0) {
	printf(" (%d)",antidonorj_knowni[i]);
      }
    }
    printf("\n");
#endif


    i = a = b = j = 0;
    while (i < antiacceptori_nsites && a < antidonora_nsites) {
      if ((splice_pos_1 = antiacceptori_positions[i]) < antidonora_positions[a]) {
	i++;
      } else if (splice_pos_1 > antidonora_positions[a]) {
	a++;
      } else {
	while (b < antiacceptorb_nsites && antiacceptorb_positions[b] <= splice_pos_1) {
	  b++;
	}
	while (j < antidonorj_nsites && antidonorj_positions[j] <= splice_pos_1) {
	  j++;
	}
	matchp = false;
	while (b < antiacceptorb_nsites && j < antidonorj_nsites && matchp == false) {
	  if ((splice_pos_2 = antiacceptorb_positions[b]) < antidonorj_positions[j]) {
	    b++;
	  } else if (splice_pos_2 > antidonorj_positions[j]) {
	    j++;
	  } else {
	    segmenti_nmismatches = Genome_count_mismatches_substring(query_compress,/*left*/segmenti_left,/*pos5*/0,/*pos3*/splice_pos_1,
								     plusp,genestrand);
	    segmentm_nmismatches = Genome_count_mismatches_substring(query_compress,/*left*/segmentm_left,/*pos5*/splice_pos_1,/*pos3*/splice_pos_2,
								     plusp,genestrand);
	    segmentj_nmismatches = Genome_count_mismatches_substring(query_compress,/*left*/segmentj_left,/*pos5*/splice_pos_2,/*pos3*/querylength,
								     plusp,genestrand);
	    
	    if ((nmismatches = segmenti_nmismatches + segmentm_nmismatches + segmentj_nmismatches) <= best_nmismatches) {
	      if (antiacceptori_knowni[i] >= 0) {
		probi = 1.0;
	      } else {
		probi = Maxent_hr_antiacceptor_prob(segmenti_left + splice_pos_1,segmenti_chroffset);
	      }
	    
	      if (antidonora_knowni[a] >= 0) {
		proba = 1.0;
	      } else {
		proba = Maxent_hr_antidonor_prob(segmentm_left + splice_pos_1,segmentm_chroffset);
	      }

	      if (antiacceptorb_knowni[b] >= 0) {
		probb = 1.0;
	      } else {
		probb = Maxent_hr_antiacceptor_prob(segmentm_left + splice_pos_2,segmentm_chroffset);
	      }

	      if (antidonorj_knowni[j] >= 0) {
		probj = 1.0;
	      } else {
		probj = Maxent_hr_antidonor_prob(segmentj_left + splice_pos_2,segmentj_chroffset);
	      }

	      debug2(
		      if (plusp == true) {
			printf("plus antisense splice_pos  %d, %d, i.antiacceptor %f, m.antidonor %f, m.antiacceptor %f, j.antidonor %f\n",
			       splice_pos_1,splice_pos_2,probi,proba,probb,probj);
		      } else {
			printf("minus sense splice_pos  %d, %d, i.antiacceptor %f, m.antidonor %f, m.antiacceptor %f, j.antidonor %f\n",
			       splice_pos_1,splice_pos_2,probi,proba,probb,probj);
		      });

	      if (nmismatches < best_nmismatches ||
		  (nmismatches == best_nmismatches && probi + proba + probb + probj > best_prob)) {
		/* Success */
		best_nmismatches = nmismatches;
		best_prob = probi + proba + probb + probj;

		best_acceptor1_knowni = antiacceptori_knowni[i];
		best_donor1_knowni = antidonora_knowni[a];
		best_acceptor2_knowni = antiacceptorb_knowni[b];
		best_donor2_knowni = antidonorj_knowni[j];
		best_acceptor1_prob = probi;
		best_donor1_prob = proba;
		best_acceptor2_prob = probb;
		best_donor2_prob = probj;
		best_splice_pos_1 = splice_pos_1;
		best_splice_pos_2 = splice_pos_2;
		best_segmenti_nmismatches = segmenti_nmismatches;
		best_segmentm_nmismatches = segmentm_nmismatches;
		best_segmentj_nmismatches = segmentj_nmismatches;
		orig_plusp = false;
	      }
	    }
	    /* b++; j++; Don't advance b or j, so next i/a can match */
	    matchp = true;
	  }
	}
	i++;
	a++;
      }
    }


    if (best_prob > 0.0) {
      debug2(printf("best_prob = %f at splice_pos %d and %d\n",best_prob,best_splice_pos_1,best_splice_pos_2));
      if (orig_plusp == true) {
	/* Originally from plus strand.  No complement. */
	sensep = (plusp == true) ? true : false;
	sensedir = (plusp == true) ? SENSE_FORWARD : SENSE_ANTI;

	donor = Substring_new_donor(best_donor1_knowni,/*joffset*/0,best_splice_pos_1,best_segmenti_nmismatches,
				    best_donor1_prob,/*left*/segmenti_left,query_compress,
				    querylength,plusp,genestrand,sensep,
				    segmenti_chrnum,segmenti_chroffset,segmenti_chrhigh,segmenti_chrlength);

	shortexon = Substring_new_shortexon(best_acceptor1_knowni,best_donor2_knowni,/*joffset*/0,
					    /*acceptor_pos*/best_splice_pos_1,/*donor_pos*/best_splice_pos_2,best_segmentm_nmismatches,
					    /*acceptor_prob*/best_acceptor1_prob,/*donor_prob*/best_donor2_prob,
					    /*left*/segmentm_left,query_compress,
					    querylength,plusp,genestrand,sensep,/*acceptor_ambp*/false,/*donor_ambp*/false,
					    segmentm_chrnum,segmentm_chroffset,segmentm_chrhigh,segmentm_chrlength);

	acceptor = Substring_new_acceptor(best_acceptor2_knowni,/*joffset*/0,best_splice_pos_2,best_segmentj_nmismatches,
					  best_acceptor2_prob,/*left*/segmentj_left,query_compress,
					  querylength,plusp,genestrand,sensep,
					  segmentj_chrnum,segmentj_chroffset,segmentj_chrhigh,segmentj_chrlength);

	if (donor == NULL || shortexon == NULL || acceptor == NULL) {
	  if (donor != NULL) Substring_free(&donor);
	  if (shortexon != NULL) Substring_free(&shortexon);
	  if (acceptor != NULL) Substring_free(&acceptor);
	} else {
	  *nhits += 1;
	  *segmenti_usedp = *segmentm_usedp = *segmentj_usedp = true;

	  donor_support = best_splice_pos_1;
	  middle_support = best_splice_pos_2 - best_splice_pos_1;
	  acceptor_support = querylength - best_splice_pos_2;
	  sufficient1p = sufficient_splice_prob_local(donor_support,best_segmenti_nmismatches,best_donor1_prob);
	  sufficient2p = sufficient_splice_prob_local(middle_support,best_segmentm_nmismatches,best_acceptor1_prob);
	  sufficient3p = sufficient_splice_prob_local(middle_support,best_segmentm_nmismatches,best_donor2_prob);
	  sufficient4p = sufficient_splice_prob_local(acceptor_support,best_segmentj_nmismatches,best_acceptor2_prob);
	  if (sufficient1p && sufficient2p && sufficient3p && sufficient4p) {
	    hits = List_push(hits,(void *) Stage3end_new_shortexon(&(*found_score),donor,acceptor,shortexon,
								   /*acceptor_distance*/segmentm_left - segmenti_left,
								   /*donor_distance*/segmentj_left - segmentm_left,
								   /*amb_nmatches_donor*/0,/*amb_nmatches_acceptor*/0,
								   /*ambi_left*/NULL,/*ambi_right*/NULL,
								   /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
								   /*copy_donor_p*/false,/*copy_acceptor_p*/false,/*copy_shortexon_p*/false,
								   splicing_penalty,querylength,sensedir,sarrayp));
	  } else if (subs_or_indels_p == true) {
	    /* Don't alter hits */
	    if (donor != NULL) Substring_free(&donor);
	    if (shortexon != NULL) Substring_free(&shortexon);
	    if (acceptor != NULL) Substring_free(&acceptor);
	  } else if (donor_support < LOWPROB_SUPPORT || acceptor_support < LOWPROB_SUPPORT) {
	    if (donor != NULL) Substring_free(&donor);
	    if (shortexon != NULL) Substring_free(&shortexon);
	    if (acceptor != NULL) Substring_free(&acceptor);
	  } else if ((sufficient1p || sufficient2p) && (sufficient3p || sufficient4p)) {
	    *lowprob = List_push(*lowprob,
				 (void *) Stage3end_new_shortexon(&(*found_score),donor,acceptor,shortexon,
								  /*acceptor_distance*/segmentm_left - segmenti_left,
								  /*donor_distance*/segmentj_left - segmentm_left,
								  /*amb_nmatches_donor*/0,/*amb_nmatches_acceptor*/0,
								  /*ambi_left*/NULL,/*ambi_right*/NULL,
								  /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
								  /*copy_donor_p*/false,/*copy_acceptor_p*/false,/*copy_shortexon_p*/false,
								  splicing_penalty,querylength,sensedir,sarrayp));
	  } else {
	    if (donor != NULL) Substring_free(&donor);
	    if (shortexon != NULL) Substring_free(&shortexon);
	    if (acceptor != NULL) Substring_free(&acceptor);
	  }
	}

      } else {
	/* Originally from minus strand.  Complement. */
	sensep = (plusp == true) ? false : true;
	sensedir = (plusp == true) ? SENSE_ANTI : SENSE_FORWARD;

	donor = Substring_new_donor(best_donor2_knowni,/*joffset*/0,best_splice_pos_2,best_segmentj_nmismatches,
				    best_donor2_prob,/*left*/segmentj_left,query_compress,
				    querylength,plusp,genestrand,sensep,
				    segmentj_chrnum,segmentj_chroffset,segmentj_chrhigh,segmentj_chrlength);

	shortexon = Substring_new_shortexon(best_acceptor2_knowni,best_donor1_knowni,/*joffset*/0,
					    /*acceptor_pos*/best_splice_pos_2,/*donor_pos*/best_splice_pos_1,best_segmentm_nmismatches,
					    /*acceptor_prob*/best_acceptor2_prob,/*donor_prob*/best_donor1_prob,
					    /*left*/segmentm_left,query_compress,querylength,
					    plusp,genestrand,sensep,/*acceptor_ambp*/false,/*donor_ambp*/false,
					    segmentm_chrnum,segmentm_chroffset,segmentm_chrhigh,segmentm_chrlength);

	acceptor = Substring_new_acceptor(best_acceptor1_knowni,/*joffset*/0,best_splice_pos_1,best_segmenti_nmismatches,
					  best_acceptor1_prob,/*left*/segmenti_left,query_compress,
					  querylength,plusp,genestrand,sensep,
					  segmenti_chrnum,segmenti_chroffset,segmenti_chrhigh,segmenti_chrlength);

	if (donor == NULL || shortexon == NULL || acceptor == NULL) {
	  if (donor != NULL) Substring_free(&donor);
	  if (shortexon != NULL) Substring_free(&shortexon);
	  if (acceptor != NULL) Substring_free(&acceptor);
	} else {
	  *nhits += 1;
	  *segmenti_usedp = *segmentm_usedp = *segmentj_usedp = true;

	  acceptor_support = best_splice_pos_1;
	  middle_support = best_splice_pos_2 - best_splice_pos_1;
	  donor_support = querylength - best_splice_pos_2;
	  sufficient1p = sufficient_splice_prob_local(acceptor_support,best_segmenti_nmismatches,best_acceptor1_prob);
	  sufficient2p = sufficient_splice_prob_local(middle_support,best_segmentm_nmismatches,best_donor1_prob);
	  sufficient3p = sufficient_splice_prob_local(middle_support,best_segmentm_nmismatches,best_acceptor2_prob);
	  sufficient4p = sufficient_splice_prob_local(donor_support,best_segmentj_nmismatches,best_donor2_prob);
	  if (sufficient1p && sufficient2p && sufficient3p && sufficient4p) {
	    hits = List_push(hits,(void *) Stage3end_new_shortexon(&(*found_score),donor,acceptor,shortexon,
								   /*acceptor_distance*/segmentj_left - segmentm_left,
								   /*donor_distance*/segmentm_left - segmenti_left,
								   /*amb_nmatches_donor*/0,/*amb_nmatches_acceptor*/0,
								   /*ambi_left*/NULL,/*ambi_right*/NULL,
								   /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
								   /*copy_donor_p*/false,/*copy_acceptor_p*/false,/*copy_shortexon_p*/false,
								   splicing_penalty,querylength,sensedir,sarrayp));
	  } else if (subs_or_indels_p == true) {
	    /* Don't alter hits */
	    if (donor != NULL) Substring_free(&donor);
	    if (shortexon != NULL) Substring_free(&shortexon);
	    if (acceptor != NULL) Substring_free(&acceptor);
	  } else if (donor_support < LOWPROB_SUPPORT || acceptor_support < LOWPROB_SUPPORT) {
	    if (donor != NULL) Substring_free(&donor);
	    if (shortexon != NULL) Substring_free(&shortexon);
	    if (acceptor != NULL) Substring_free(&acceptor);
	  } else if ((sufficient1p || sufficient2p) && (sufficient3p || sufficient4p)) {
	    *lowprob = List_push(*lowprob,
				 (void *) Stage3end_new_shortexon(&(*found_score),donor,acceptor,shortexon,
								  /*acceptor_distance*/segmentj_left - segmentm_left,
								  /*donor_distance*/segmentm_left - segmenti_left,
								  /*amb_nmatches_donor*/0,/*amb_nmatches_acceptor*/0,
								  /*ambi_left*/NULL,/*ambi_right*/NULL,
								  /*amb_nmismatches_left*/NULL,/*amb_nmismatches_right*/NULL,
								  /*copy_donor_p*/false,/*copy_acceptor_p*/false,/*copy_shortexon_p*/false,
								  splicing_penalty,querylength,sensedir,sarrayp));
	  } else {
	    if (donor != NULL) Substring_free(&donor);
	    if (shortexon != NULL) Substring_free(&shortexon);
	    if (acceptor != NULL) Substring_free(&acceptor);
	  }
	}
      }
    }
  }

  return hits;
}


