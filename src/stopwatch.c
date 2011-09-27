static char rcsid[] = "$Id: stopwatch.c,v 1.9 2005/07/08 14:43:19 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "stopwatch.h"

#ifdef HAVE_UNISTD_H
#include <unistd.h>		/* For sysconf */
#endif
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>		/* For sys/times.h under AT&T System V Interface */
#endif
#include <sys/times.h>

static struct tms start;
static struct tms stop;
static clock_t start_elapsed;
static clock_t stop_elapsed;

void
Stopwatch_start (void) {
  start_elapsed = times(&start);
}

/* Returns user seconds elapsed since start of process. */
/* struct tms is defined in <sys/times.h> (see man times); tms values are
   in "clock ticks per second", CLK_TCK */
double 
Stopwatch_stop (void) {
  long clk_tck = sysconf(_SC_CLK_TCK);

  stop_elapsed = times(&stop);

  /* user time is in stop.tms_utime */
  return (double) (stop_elapsed - start_elapsed)/(double) clk_tck;
}

