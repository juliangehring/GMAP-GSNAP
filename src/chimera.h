/* $Id: chimera.h,v 1.13 2006/03/02 22:52:26 twu Exp $ */
#ifndef CHIMERA_INCLUDED
#define CHIMERA_INCLUDED
#include "bool.h"
#include "reader.h"		/* For cDNAEnd_T */
#include "sequence.h"
#include "genome.h"
#include "stage3.h"

#define T Chimera_T
typedef struct T *T;

extern int 
Chimera_pos (T this);
extern int
Chimera_equivpos (T this);
extern int
Chimera_cdna_direction (T this);
extern double
Chimera_donor_prob (T this);
extern double
Chimera_acceptor_prob (T this);

extern T
Chimera_new (int chimerapos, int chimeraequivpos, int exonexonpos, int cdna_direction,
	     double donor_prob, double acceptor_prob);
extern void
Chimera_free (T *old);
extern void
Chimera_print (T this);

extern int
Chimera_alignment_break (int *newstart, int *newend, Stage3_T stage3, int queryntlength, double fthreshold);
extern void
Chimera_bestpath (int *five_score, int *three_score, int *chimerapos, int *chimeraequivpos, int *bestfrom, int *bestto, 
		  Stage3_T *stage3array_sub1, int npaths_sub1, Stage3_T *stage3array_sub2, int npaths_sub2, 
		  int queryntlength);
extern void
Chimera_find_exonexon (T this, Stage3_T left_part, Stage3_T right_part,
		       Genome_T genome);

extern bool
Chimera_exonexon_p (int *exonexonpos, int *cdna_direction, double *donor_prob, double *acceptor_prob,
		    Stage3_T left_part, Stage3_T right_part, Genome_T genome, int querylength);
#undef T
#endif
