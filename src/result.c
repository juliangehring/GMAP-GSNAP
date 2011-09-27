static char rcsid[] = "$Id: result.c,v 1.53 2006/10/09 16:58:10 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "result.h"
#include <stdlib.h>
#include "mem.h"


#define T Result_T
struct T {
  int id;
  int worker_id;

  Chimera_T chimera;		/* NULL indicates not a chimera */
  Stage3_T *array;
  int npaths;
  Failure_T failuretype;

  /* Diagnostic info */
#ifndef PMAP
  double query_oligodepth;
  int query_badoligos;
  int query_repoligos;
  int query_trimoligos;
  int query_trim_start;
  int query_trim_end;
#endif

  double stage1_runtime;
  int stage1_trial;
  int stage1_firstpair_searched_5;
  int stage1_firstpair_searched_3;
  int stage1_stutter_searched_5;
  int stage1_stutter_searched_3;
  int stage1_stutter_matchpairs;
  int stage1_stutter_matches_5;
  int stage1_stutter_matches_3;
  int stage1_dangling_5;
  int stage1_dangling_3;
  int stage1_dangling_matchpairs;
  int stage1_sampling_rounds;
  int stage1_sampling_nskip;
  int stage1_nentries_total;
  int stage1_nentries_used;
};


int
Result_id (T this) {
  return this->id;
}


int
Result_worker_id (T this) {
  return this->worker_id;
}


Chimera_T
Result_chimera (T this) {
  return this->chimera;
}


Stage3_T *
Result_array (int *npaths, T this) {
  *npaths = this->npaths;
  return this->array;
}


Failure_T
Result_failuretype (T this) {
  return this->failuretype;
}


T
Result_new (int id, int worker_id, Chimera_T chimera,
	    Stage3_T *array, int npaths, Failure_T failuretype,
#ifndef PMAP
	    double oligodepth, int badoligos, int repoligos, int trimoligos, int trim_start, int trim_end,
#endif
	    Stage1_T stage1) {
  T new = (T) MALLOC(sizeof(*new));

  new->id = id;
  new->worker_id = worker_id;
  new->chimera = chimera;
  new->array = array;
  new->npaths = npaths;
  new->failuretype = failuretype;

#ifndef PMAP
  new->query_oligodepth = oligodepth;
  new->query_badoligos = badoligos;
  new->query_repoligos = repoligos;
  new->query_trimoligos = trimoligos;
  new->query_trim_start = trim_start;
  new->query_trim_end = trim_end;;
#endif

  if (stage1 == NULL) {
    new->stage1_trial = -1;
  } else {
    new->stage1_runtime = Stage1_runtime(stage1);
    new->stage1_trial = Stage1_trial(stage1);
    new->stage1_firstpair_searched_5 = Stage1_firstpair_searched_5(stage1);
    new->stage1_firstpair_searched_3 = Stage1_firstpair_searched_3(stage1);
    new->stage1_stutter_searched_5 = Stage1_stutter_searched_5(stage1);
    new->stage1_stutter_searched_3 = Stage1_stutter_searched_3(stage1);
    new->stage1_stutter_matchpairs = Stage1_stutter_matchpairs(stage1);
    new->stage1_stutter_matches_5 = Stage1_stutter_matches_5(stage1);
    new->stage1_stutter_matches_3 = Stage1_stutter_matches_3(stage1);
    new->stage1_dangling_5 = Stage1_dangling_5(stage1);
    new->stage1_dangling_3 = Stage1_dangling_3(stage1);
#ifndef PMAP
    new->stage1_dangling_matchpairs = Stage1_dangling_matchpairs(stage1);
#endif
    new->stage1_sampling_rounds = Stage1_sampling_rounds(stage1);
    new->stage1_sampling_nskip = Stage1_sampling_nskip(stage1);
    new->stage1_nentries_total = Stage1_nentries_total(stage1);
    new->stage1_nentries_used = Stage1_nentries_used(stage1);
  }

  return new;
}

void
Result_free (T *old) {
#ifdef BETATEST
  Chimera_T chimera;
#endif
  Stage3_T stage3;
  int i;

  if (*old) {
#ifdef BETATEST    
    if ((chimera = (*old)->chimera) != NULL) {
      Chimera_free(&chimera);
    }
#endif
    if ((*old)->array) {
      for (i = 0; i < (*old)->npaths; i++) {
	stage3 = (*old)->array[i];
	Stage3_free(&stage3);
      }
      FREE((*old)->array);
    }

    FREE(*old);
  }
  return;
}


void
Result_print_diagnostics (T this) {
#ifndef PMAP
  printf("Query check oligodepth: %f\n",this->query_oligodepth);
  printf("Query check badoligos: %d/%d\n",this->query_badoligos,this->query_trimoligos);
  printf("Query check repoligos: %d/%d\n",this->query_repoligos,this->query_trimoligos);
  printf("Query check trim bounds: %d..%d\n",this->query_trim_start,this->query_trim_end);
#endif

  if (this->stage1_trial >= 0) {
    printf("Stage 1 runtime: %.3f sec\n",this->stage1_runtime);
    printf("Stage 1 trial: %d\n",this->stage1_trial);
    printf("Stage 1 firstpair_searched_5: %d\n",this->stage1_firstpair_searched_5);
    printf("Stage 1 firstpair_searched_3: %d\n",this->stage1_firstpair_searched_3);
    printf("Stage 1 stutter_searched_5: %d\n",this->stage1_stutter_searched_5);
    printf("Stage 1 stutter_searched_3: %d\n",this->stage1_stutter_searched_3);
    printf("Stage 1 stutter_matchpairs: %d\n",this->stage1_stutter_matchpairs);
    printf("Stage 1 stutter_matches_5: %d (%d dangling)\n",this->stage1_stutter_matches_5,this->stage1_dangling_5);
    printf("Stage 1 stutter_matches_3: %d (%d dangling)\n",this->stage1_stutter_matches_3,this->stage1_dangling_3);
#ifndef PMAP
    printf("Stage 1 dangling_matchpairs: %d\n",this->stage1_dangling_matchpairs);
#endif
    printf("Stage 1 sampling rounds: %d",this->stage1_sampling_rounds);
    if (this->stage1_sampling_rounds == 0) {
      printf("\n");
    } else {
      printf(" (nskip = %d)\n",this->stage1_sampling_nskip);
    }
    printf("Stage 1 nentries_total: %d\n",this->stage1_nentries_total);
    printf("Stage 1 nentries_used: %d\n",this->stage1_nentries_used);
  }

  return;
}
