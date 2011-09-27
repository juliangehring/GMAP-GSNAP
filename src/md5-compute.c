static char rcsid[] = "$Id: md5-compute.c,v 1.8 2005/07/08 07:58:33 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include "bool.h"
#include "mem.h"
#include "sequence.h"
#include "md5.h"


static bool onesequencep = false;

#ifdef __STRICT_ANSI__
int getopt (int argc, char *const argv[], const char *optstring);
#endif

int
main (int argc, char *argv[]) {
  Sequence_T sequence;
  unsigned char *digest;
  int nextchar;

  int opt;
  extern int optind;
  extern char *optarg;

  while ((opt = getopt(argc,argv,"1")) != -1) {
    switch (opt) {
    case '1': onesequencep = true; break;
    }
  }

  argc -= optind;
  argv += optind;

  if (onesequencep == true) {
    sequence = Sequence_read(&nextchar,stdin);
    digest = MD5_compute((unsigned char *) Sequence_fullpointer(sequence),Sequence_fulllength(sequence));
    MD5_print(digest);
    printf("\n");
  
    FREE(digest);
    Sequence_free(&sequence);

  } else {
    while ((sequence = Sequence_read(&nextchar,stdin)) != NULL) {
      digest = MD5_compute((unsigned char *) Sequence_fullpointer(sequence),Sequence_fulllength(sequence));
      printf("%s\t",Sequence_accession(sequence));
      printf("%d\t",Sequence_fulllength(sequence));
      MD5_print(digest);
      printf("\n");
  
      FREE(digest);
      Sequence_free(&sequence);
    }
  }

  return 0;
}

