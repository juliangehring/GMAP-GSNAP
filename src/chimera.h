/* $Id: chimera.h,v 1.7 2005/06/03 20:13:32 twu Exp $ */
#ifndef CHIMERA_INCLUDED
#define CHIMERA_INCLUDED
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
Chimera_new (Stage3_T nonchimericbest, int chimerapos, int chimeraequivpos);
extern void
Chimera_free (T *old);
extern void
Chimera_print (T this);

extern int
Chimera_detect (int *margin, Stage3_T stage3, Sequence_T queryseq, double fthreshold);
extern void
Chimera_bestpath (int *chimerapos, int *chimeraequivpos, int *bestfrom, int *bestto, 
		  Stage3_T *stage3array, int npaths, int querylength);
extern void
Chimera_find_exonexon (T this, Stage3_T left_part, Stage3_T right_part,
		       Genome_T genome);

#undef T
#endif
