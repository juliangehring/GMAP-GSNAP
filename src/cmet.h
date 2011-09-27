#ifndef CMET_INCLUDED
#define CMET_INCLUDED

#include "indexdb.h"		/* For Storedoligomer_T */

extern Storedoligomer_T
Cmet_reduce_ct (Storedoligomer_T oligo);
extern Storedoligomer_T
Cmet_reduce_ga (Storedoligomer_T oligo);
extern Storedoligomer_T
Cmet_mark_a (Storedoligomer_T high, Storedoligomer_T low);
extern Storedoligomer_T
Cmet_mark_c (Storedoligomer_T high, Storedoligomer_T low);
extern Storedoligomer_T
Cmet_mark_g (Storedoligomer_T high, Storedoligomer_T low);
extern Storedoligomer_T
Cmet_mark_t (Storedoligomer_T high, Storedoligomer_T low);

#endif

