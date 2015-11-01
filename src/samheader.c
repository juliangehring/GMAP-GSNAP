static char rcsid[] = "$Id: samheader.c 112663 2013-10-25 16:57:13Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "samheader.h"


void
SAM_header_print_HD (FILE *fp, int nworkers, bool orderedp) {
  fprintf(fp,"@HD");
  fprintf(fp,"\tVN:1.0");	/* or 1.4 */
  if (nworkers > 1 && orderedp == false) {
    fprintf(fp,"\tSO:unsorted");
  } else {
    /* Picard does not recognize type unknown */
    /* fprintf(fp,"\tSO:unknown"); */
    fprintf(fp,"\tSO:unsorted");
  }
  fprintf(fp,"\n");
  return;
}


void
SAM_header_print_PG (FILE *fp, int argc, char **argv, int optind) {
  char **argstart;
  int c;

  fprintf(fp,"@PG");
#ifdef GSNAP
  fprintf(fp,"\tID:GSNAP");
  fprintf(fp,"\tPN:gsnap");
#elif defined(PMAP)
  fprintf(fp,"\tID:PMAP");
  fprintf(fp,"\tPN:pmap");
#else
  fprintf(fp,"\tID:GMAP");
  fprintf(fp,"\tPN:gmap");
#endif
  fprintf(fp,"\tVN:%s",PACKAGE_VERSION);

  fprintf(fp,"\tCL:");
  argstart = &(argv[-optind]);
  fprintf(fp,"%s",argstart[0]);
  for (c = 1; c < argc + optind; c++) {
    fprintf(fp," %s",argstart[c]);
  }
  fprintf(fp,"\n");
  return;
}

