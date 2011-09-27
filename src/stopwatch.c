static char rcsid[] = "$Id: stopwatch.c,v 1.10 2006/10/04 19:24:58 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "stopwatch.h"
#include "mem.h"

#ifdef HAVE_UNISTD_H
#include <unistd.h>		/* For sysconf */
#endif
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>		/* For sys/times.h under AT&T System V Interface */
#endif
#include <sys/times.h>

#define T Stopwatch_T
struct T {
  struct tms start;
  struct tms stop;
  clock_t start_elapsed;
  clock_t stop_elapsed;
};

T
Stopwatch_new () {
  T new = (T) MALLOC(sizeof(*new));
  return new;
}

void
Stopwatch_free (T *old) {
  if (*old) {
    FREE(*old);
  }
  return;
}

void
Stopwatch_start (T this) {
  if (this != NULL) {
    this->start_elapsed = times(&this->start);
  }
  return;
}

/* Returns user seconds elapsed since start of process. */
/* struct tms is defined in <sys/times.h> (see man times); tms values are
   in "clock ticks per second", CLK_TCK */
double 
Stopwatch_stop (T this) {
  long clk_tck = sysconf(_SC_CLK_TCK);

  if (this == NULL) {
    return 0.0;
  } else {
    this->stop_elapsed = times(&this->stop);

    /* user time is in stop.tms_utime */
    return (double) (this->stop_elapsed - this->start_elapsed)/(double) clk_tck;
  }
}

