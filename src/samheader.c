static char rcsid[] = "$Id: samheader.c 155503 2014-12-16 22:22:55Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "samheader.h"


#define CHUNK 1024

void
SAM_header_change_HD_tosorted_stdout (FILE *fp, int headerlen) {
  char buffer[CHUNK], c, c0, c1, c2;

  /* @HD */
  while (headerlen > 0 && (c = fgetc(fp)) != '\t') {
    putchar(c);
    headerlen--;
  }
  if (headerlen > 0) {
    putchar('\t');
    headerlen--;
  }

  /* VN */
  while (headerlen > 0 && (c = fgetc(fp)) != '\t') {
    putchar(c);
    headerlen--;
  }
  if (headerlen > 0) {
    putchar('\t');
    headerlen--;
  }

  if (headerlen > 3) {
    /* SO: */
    c0 = fgetc(fp);
    c1 = fgetc(fp);
    c2 = fgetc(fp);
    printf("%c%c%c",c0,c1,c2);
    headerlen -= 3;

    if (c0 == 'S' && c1 == 'O' && c2 == ':') {
      printf("coordinate\n");
      while (headerlen > 0 && fgetc(fp) != '\n') {
	/* Skip given SO value */
	headerlen--;
      }
      headerlen--;
    }
  }

  while (headerlen > CHUNK) {
    fread(buffer,sizeof(char),CHUNK,fp);
    fwrite(buffer,sizeof(char),CHUNK,stdout);
    headerlen -= CHUNK;
  }
  if (headerlen > 0) {
    fread(buffer,sizeof(char),headerlen,fp);
    fwrite(buffer,sizeof(char),headerlen,stdout);
  }

  return;
}


void
SAM_header_change_HD_tosorted_split (FILE *fp, int headerlen, FILE **outputs, int noutputs) {
  char buffer[CHUNK], c, c0, c1, c2;
  int i;

  /* @HD */
  while (headerlen > 0 && (c = fgetc(fp)) != '\t') {
    for (i = 1; i <= noutputs; i++) {
      putc(c,outputs[i]);
    }
    headerlen--;
  }
  if (headerlen > 0) {
    for (i = 1; i <= noutputs; i++) {
      putc('\t',outputs[i]);
    }
    headerlen--;
  }

  /* VN */
  while (headerlen > 0 && (c = fgetc(fp)) != '\t') {
    for (i = 1; i <= noutputs; i++) {
      putc(c,outputs[i]);
    }
    headerlen--;
  }
  if (headerlen > 0) {
    for (i = 1; i <= noutputs; i++) {
      putc('\t',outputs[i]);
    }
    headerlen--;
  }

  if (headerlen > 3) {
    /* SO: */
    c0 = fgetc(fp);
    c1 = fgetc(fp);
    c2 = fgetc(fp);
    for (i = 1; i <= noutputs; i++) {
      fprintf(outputs[i],"%c%c%c",c0,c1,c2);
    }
    headerlen -= 3;

    if (c0 == 'S' && c1 == 'O' && c2 == ':') {
      for (i = 1; i <= noutputs; i++) {
	fprintf(outputs[i],"coordinate\n");
      }
      while (headerlen > 0 && fgetc(fp) != '\n') {
	/* Skip given SO value */
	headerlen--;
      }
      headerlen--;
    }
  }

  while (headerlen > CHUNK) {
    fread(buffer,sizeof(char),CHUNK,fp);
    for (i = 1; i <= noutputs; i++) {
      fwrite(buffer,sizeof(char),CHUNK,outputs[i]);
    }
    headerlen -= CHUNK;
  }
  if (headerlen > 0) {
    fread(buffer,sizeof(char),headerlen,fp);
    for (i = 1; i <= noutputs; i++) {
      fwrite(buffer,sizeof(char),headerlen,outputs[i]);
    }
  }

  return;
}


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

#if 0
  /* Algorithm types.  Now using XG instead. */
  fprintf(fp,"@PG");
#ifdef GSNAP
  fprintf(fp,"\tID:A");
  fprintf(fp,"\tPN:gsnap-suffix-array");
  fprintf(fp,"\n");

  fprintf(fp,"\tID:M");
  fprintf(fp,"\tPN:gsnap-gmap-method");
  fprintf(fp,"\n");

  fprintf(fp,"\tID:O");
  fprintf(fp,"\tPN:gsnap-overlap-merge");
  fprintf(fp,"\n");
#endif

#endif

  return;
}

int
SAM_header_length (int *lastchar, FILE *fp) {
  int headerlen = 0;
  int c;

  while (!feof(fp) && (c = getc(fp)) == '@') {
    headerlen++;
    while (!feof(fp) && (c = getc(fp)) != '\n') {
      headerlen++;
    }
    headerlen++;
  }
  /* headerlen++; -- Don't count in header, but as part of first SAM line */

  *lastchar = c;
  return headerlen;
}

