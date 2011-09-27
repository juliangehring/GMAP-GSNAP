static char rcsid[] = "$Id: assert.c 27450 2010-08-05 19:02:48Z twu $";
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

