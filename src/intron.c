static char rcsid[] = "$Id: intron.c,v 1.8 2005/02/15 01:58:50 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "intron.h"


bool
Intron_canonical_fwd_p (char donor1, char donor2, char acceptor2, char acceptor1) {
  if (donor1 == 'G' && donor2 == 'T' &&
      acceptor2 == 'A' && acceptor1 == 'G') {
    return true;
  } else {
    return false;
  }
}

bool
Intron_canonical_rev_p (char donor1, char donor2, char acceptor2, char acceptor1) {
  if (donor1 == 'C' && donor2 == 'T' &&
      acceptor2 == 'A' && acceptor1 == 'C') {
    return true;
  } else {
    return false;
  }
}

bool
Intron_gcag_fwd_p (char donor1, char donor2, char acceptor2, char acceptor1) {
  if (donor1 == 'G' && donor2 == 'C' &&
      acceptor2 == 'A' && acceptor1 == 'G') {
    return true;
  } else {
    return false;
  }
}

bool
Intron_atac_fwd_p (char donor1, char donor2, char acceptor2, char acceptor1) {
  if (donor1 == 'A' && donor2 == 'T' &&
      acceptor2 == 'A' && acceptor1 == 'C') {
    return true;
  } else {
    return false;
  }
}

bool
Intron_gcag_rev_p (char donor1, char donor2, char acceptor2, char acceptor1) {
  if (donor1 == 'C' && donor2 == 'T' &&
      acceptor2 == 'G' && acceptor1 == 'C') {
    return true;
  } else {
    return false;
  }
}

bool
Intron_atac_rev_p (char donor1, char donor2, char acceptor2, char acceptor1) {
  if (donor1 == 'G' && donor2 == 'T' &&
      acceptor2 == 'A' && acceptor1 == 'T') {
    return true;
  } else {
    return false;
  }
}

