static char rcsid[] = "$Id: intron.c,v 1.10 2006/02/23 22:06:37 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "intron.h"


int
Intron_type (char left1, char left2, char right2, char right1, int cdna_direction) {
  int introntype, leftdi, rightdi;

  if (left1 == 'G' && left2 == 'T') {
    leftdi = LEFT_GT;
  } else if (left1 == 'G' && left2 == 'C') {
    leftdi = LEFT_GC;
  } else if (left1 == 'A' && left2 == 'T') {
    leftdi = LEFT_AT;
#ifndef PMAP
  } else if (left1 == 'C' && left2 == 'T') {
    leftdi = LEFT_CT;
#endif
  } else {
    return NONINTRON;
  }

  if (right2 == 'A' && right1 == 'G') {
    rightdi = RIGHT_AG;
  } else if (right2 == 'A' && right1 == 'C') {
    rightdi = RIGHT_AC;
#ifndef PMAP
  } else if (right2 == 'G' && right1 == 'C') {
    rightdi = RIGHT_GC;
  } else if (right2 == 'A' && right1 == 'T') {
    rightdi = RIGHT_AT;
#endif
  } else {
    return NONINTRON;
  }

  if ((introntype = leftdi & rightdi) == 0x00) {
    return NONINTRON;
  } else if (cdna_direction > 0) {
    if (introntype < 0x08) {
      return NONINTRON;
    } else {
      return introntype;
    }
  } else if (cdna_direction < 0) {
    if (introntype > 0x04) {
      return NONINTRON;
    } else {
      return introntype;
    }
  } else {
    abort();
  }
}


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

