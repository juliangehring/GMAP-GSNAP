#ifndef MAXENT_HR_INCLUDED
#define MAXENT_HR_INCLUDED
#include "types.h"
#include "genomicpos.h"

extern void
Maxent_hr_setup (UINT4 *ref_blocks_in);

extern double
Maxent_hr_donor_prob (Genomicpos_T splice_pos);

extern double
Maxent_hr_acceptor_prob (Genomicpos_T splice_pos);

extern double
Maxent_hr_antidonor_prob (Genomicpos_T splice_pos);

extern double
Maxent_hr_antiacceptor_prob (Genomicpos_T splice_pos);

#endif

