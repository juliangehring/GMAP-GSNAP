/* $Id: intron.h,v 1.5 2005/02/07 23:56:56 twu Exp $ */
#ifndef INTRON_INCLUDED
#define INTRON_INCLUDED
#include "bool.h"

extern bool
Intron_canonical_fwd_p (char donor1, char donor2, char acceptor2, char acceptor1);
extern bool
Intron_canonical_rev_p (char donor1, char donor2, char acceptor2, char acceptor1);
extern bool
Intron_gcag_fwd_p (char donor1, char donor2, char acceptor2, char acceptor1);
extern bool
Intron_gcag_rev_p (char donor1, char donor2, char acceptor2, char acceptor1);
extern bool
Intron_atac_fwd_p (char donor1, char donor2, char acceptor2, char acceptor1);
extern bool
Intron_atac_rev_p (char donor1, char donor2, char acceptor2, char acceptor1);

#endif
