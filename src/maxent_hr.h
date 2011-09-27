#ifndef MAXENT_HR_INCLUDED
#define MAXENT_HR_INCLUDED
#include "types.h"
#include "genomicpos.h"

extern double
Maxent_hr_donor_prob (Genomicpos_T splice_pos, UINT4 *blocks);

extern double
Maxent_hr_acceptor_prob (Genomicpos_T splice_pos, UINT4 *blocks);

extern double
Maxent_hr_antidonor_prob (Genomicpos_T splice_pos, UINT4 *blocks);

extern double
Maxent_hr_antiacceptor_prob (Genomicpos_T splice_pos, UINT4 *blocks);

#endif

