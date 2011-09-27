/* $Id: chimera.h 46206 2011-08-31 22:13:13Z twu $ */
#ifndef CHIMERA_INCLUDED
#define CHIMERA_INCLUDED

typedef struct Chimera_T *Chimera_T;

#include <stdio.h>
#include "bool.h"
#include "genome.h"
#include "stage3.h"

#define T Chimera_T

extern int 
Chimera_pos (T this);
extern int
Chimera_equivpos (T this);
extern int
Chimera_cdna_direction (T this);
extern void
Chimera_print_sam_tag (FILE *fp, T this);
extern double
Chimera_donor_prob (T this);
extern double
Chimera_acceptor_prob (T this);

extern T
Chimera_new (int chimerapos, int chimeraequivpos, int exonexonpos, int cdna_direction,
	     char donor1, char donor2, char acceptor2, char acceptor1,
	     double donor_prob, double acceptor_prob);
extern void
Chimera_free (T *old);
extern void
Chimera_print (FILE *fp, T this);

extern int
Chimera_alignment_break (int *newstart, int *newend, Stage3_T stage3, int queryntlength, double fthreshold);
extern void
Chimera_bestpath (int *five_score, int *three_score, int *chimerapos, int *chimeraequivpos, int *bestfrom, int *bestto, 
		  Stage3_T *stage3array_sub1, int npaths_sub1, Stage3_T *stage3array_sub2, int npaths_sub2, 
		  int queryntlength);
#if 0
extern void
Chimera_find_exonexon_old (T this, Stage3_T left_part, Stage3_T right_part,
			   Genome_T genome, IIT_T chromosome_iit);
#endif

extern void
Chimera_find_exonexon (int *exonexonpos, int *cdna_direction,
		       char *donor1, char *donor2, char *acceptor2, char *acceptor1,
		       double *donor_prob, double *acceptor_prob,
		       Stage3_T left_part, Stage3_T right_part, Genome_T genome, IIT_T chromosome_iit,
		       int breakpoint_start, int breakpoint_end);
#undef T
#endif
