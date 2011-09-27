static char rcsid[] = "$Id: assert.c,v 1.6 2005-02-07 23:56:55 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "assert.h"

const Except_T Assert_Failed = { "Assertion Failed" };

/*
void
(assert) (int e) {
  assert(e);
}
*/

