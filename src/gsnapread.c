static char rcsid[] = "$Id: gsnapread.c 40271 2011-05-28 02:29:18Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>		/* For strcpy */
#include <strings.h>		/* For rindex */
#include <ctype.h>

#include "gsnapread.h"
#include "except.h"
#include "mem.h"
#include "bool.h"


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif



bool
Gsnapread_concordantp (char *header) {
  char *p;

  p = &(header[1]);      

  /* Advance past first tab */
  while (!isspace(*p)) p++;
  if (*p == '\t') p++;
  
  /* Advance past first space */
  while (*p != '\0' && isdigit(*p)) p++;

  if (*p == '\t') {
    return false;
  } else if (*p == ' ') {
    p++;
    if (!strncmp(p,"concordant",strlen("concordant"))) {
      return true;
    } else {
      return false;
    }
  } else {
    fprintf(stderr,"Cannot parse %s",header);
    abort();
  }
}


char *
Gsnapread_quality_string (char *header) {
  char *quality_string;
  char *p, *q, *r;
  int readlength, quality_string_length;

  p = &(header[1]);      

  /* Advance past first tab */
  q = p;
  while (!isspace(*q)) q++;
  readlength = (q - p)/sizeof(char);
  if (*q == '\t') q++;

  p = q;
  /* Advance past second tab */
  while (*p != '\0' && *p != '\t') p++;
  if (*p == '\t') p++;
  
  r = q = p;			/* After second tab */
  
  /* See if there is a third tab */
  while (*r != '\0' && *r != '\t') r++;
  if (*r == '\0') {
    /* No third tab */
    return (char *) NULL;
  } else {
    while (!isspace(*q)) q++;
    quality_string_length = (q - p)/sizeof(char);
    quality_string = (char *) CALLOC(readlength+1,sizeof(char));
    if (quality_string_length == readlength) {
      strncpy(quality_string,p,readlength);
      return quality_string;
    } else {
      fprintf(stderr,"Warning quality string length not equal to read length at %s\n",p);
      return (char *) NULL;
    }
  }

}



int
Gsnapread_nhits (char *header) {
  int nhits;
  char *p;

  p = &(header[1]);      

  /* Advance past first tab */
  while (!isspace(*p)) p++;
  if (*p == '\t') p++;
  
  if (sscanf(p,"%d",&nhits) != 1) {
    fprintf(stderr,"Can't parse nhits in %s\n",header);
    abort();
  } else {
    return nhits;
  }
}

/* This gets query length, not length after chopping adapters */
int
Gsnapread_readlength_preadapter (char *header) {
  int readlength = 0;
  char *p;

  p = &(header[1]);      

  while (!isspace(*p)) {
    p++;
    readlength++;
  }

  return readlength;
}

int
Gsnapread_readlength (char *line) {
  int readlength = 0;
  char *p;

  p = &(line[1]);

  while (!isspace(*p) && *p != '*') {
    p++;
    readlength++;
  }

  return readlength;
}

char *
Gsnapread_accession (char *header) {
  char *accession;
  char *p, *q, *r;
  int len;

  p = &(header[1]);      

  /* Advance past first tab */
  while (!isspace(*p)) p++;
  if (*p == '\t') p++;

  /* Advance past second tab */
  while (*p != '\0' && *p != '\t') p++;
  if (*p == '\t') p++;
  
  q = p;
  while (*q != '\0' && !isspace(*q)) q++;

  /* See if there is a third tab */
  r = q;
  while (*r != '\0' && *r != '\t') r++;
  if (*r == '\0') {
    /* No third tab */
    len = q - p;
    accession = (char *) CALLOC(len+1,sizeof(char));
    strncpy(accession,p,len);
  } else {
    r++;			/* Advance past third tab */
    q = r;
    while (*q != '\0' && !isspace(*q)) q++;
    len = q - r;
    accession = (char *) CALLOC(len+1,sizeof(char));
    strncpy(accession,r,len);
  }

  return accession;
}


static void
get_splicesite_info (bool *sensep, double *prob, char *p) {
  
  while (*p != ':') p++;
  p++;

  if (sscanf(p,"%lf",&(*prob)) < 1) {
    fprintf(stderr,"Can't parse probability part at %s\n",p);
    abort();
  }

  if ((p = strstr(p,"dir:")) == NULL) {
    fprintf(stderr,"Can't parse dir: part of type %s\n",p);
    abort();
  } else {
    while (*p != ':') p++;
    p++;
  }
  if (!strncmp(p,"sense",strlen("sense"))) {
    *sensep = true;
  } else if (!strncmp(p,"antisense",strlen("antisense"))) {
    *sensep = false;
  } else {
    fprintf(stderr,"Don't recognize part after dir: %s\n",p);
    abort();
  }

  return;
}



void
Gsnapread_parse_line (char *line, int *query5, int *query3, char *firstend, char *secondend,
		      int *support, int *nmismatches, char *strand, char **chr,
		      Genomicpos_T *firstpos, Genomicpos_T *secondpos,
		      Genomicpos_T *donorpos, Genomicpos_T *acceptorpos, char *truestrand, bool *sensep) {
  char *p, *q;
  int len, i;
  double donorprob, acceptorprob;


  debug(printf("Entering parse_line with %s\n",line));

  /* Get query coordinates */
  /* p = AGCACTCTGTTCTCAAAACctgggggaatggagaaggcgcatacaccttagagactgcagatgcagagcaggaca	1..19	-17:16344406..16344388	acceptor:1.00,dir:antisense,sub:0,score:2	dist:820 */
  p = &(line[1]);      
  while (*p != '\0' && *p != '\t') p++;
  if (*p != '\0' && *p == '\t') p++;
  if (*p == '\0') {
    fprintf(stderr,"Can't parse query coordinates part of %s\n",line);
    abort();
  }

  /* p = 1..19	-17:16344406..16344388	acceptor:1.00,dir:antisense,sub:0,score:2	dist:820 */
  if (sscanf(p,"%d",&(*query5)) != 1) {
    fprintf(stderr,"Can't parse first query coordinate in %s\n",line);
    abort();
  }
  while (*p != '\0' && isdigit(*p)) p++;
  while (*p != '\0' && !isdigit(*p)) p++;

  /* p = 19	-17:16344406..16344388	dist:820 */
  if (sscanf(p,"%d",&(*query3)) != 1) {
    fprintf(stderr,"Can't parse second query coordinate in %s\n",line);
    abort();
  }

  *support = (*query3) - (*query5) + 1;

  /* Count mismatches */
  *nmismatches = 0;
  q = &(line[1]);
  for (i = 1; i < *query5; i++) {
    /* printf("Skipping %c\n",*q); */
    q++;
  }
  for ( ; i <= *query3; i++) {
    /* printf("i = %d, Evaluating %c\n",i,*q); */
    if (islower(*q)) {
      (*nmismatches)++;
    }
    q++;
  }
  /* printf("Done\n\n"); */

  debug(printf("Got query %d to %d, with length %d and %d mismatches\n",*query5,*query3,*support,*nmismatches));


  /* Advance to strand */
  while (*p != '\0' && *p != '\t') p++;
  if (*p != '\0' && *p == '\t') p++;

  /* p = -17:16344406..16344388	acceptor:1.00,dir:antisense,sub:0,score:2	dist:820 */
  if (*p == '\0') {
    fprintf(stderr,"Can't parse strand part of %s\n",line);
    abort();
  } else {
    *strand = *p++;
  }

  /* Get chr part */
  q = p;
  len = 0;
  while (*q != '\0' && *q != ':') {
    q++;
    len++;
  }
  if (*q == '\0') {
    fprintf(stderr,"Can't parse chr part of %s\n",line);
    abort();
  } else {
    *chr = (char *) CALLOC(len+1,sizeof(char));
    strncpy(*chr,p,len);
    (*chr)[len] = '\0';
  }

  p = ++q;
  /* p = 16344406..16344388	dist:820 */
  if (sscanf(p,"%u",&(*firstpos)) != 1) {
    fprintf(stderr,"Can't parse first chrpos in %s\n",line);
    abort();
  }

  /* Advance past first chrpos */
  while (*p != '\0' && isdigit(*p)) p++;
  while (*p != '\0' && !isdigit(*p)) p++;
  /* p = 16344388	dist:820 */
  if (sscanf(p,"%u",&(*secondpos)) != 1) {
    fprintf(stderr,"Can't parse second chrpos in %s\n",line);
    abort();
  }


  /* Advance to type */
  while (*p != '\0' && !isspace(*p)) p++;
  while (*p != '\0' && isspace(*p)) p++;
  /* p = acceptor:1.00,dir:antisense,sub:0,score:2	dist:820 */
  if (*p == '\0') {
    fprintf(stderr,"Can't parse type part of %s\n",line);
    abort();
  }

  if (!strncmp(p,"start",strlen("start"))) {
    *firstend = 'S';
#if 0
  } else if (!strncmp(p,"end",strlen("end"))) {
    *firstend = 'E';
#endif
  } else if (!strncmp(p,"amb",strlen("amb"))) {
    *firstend = 'A';
  } else if (!strncmp(p,"term",strlen("term"))) {
    *firstend = 'T';
  } else if (!strncmp(p,"ins",strlen("ins"))) {
    *firstend = 'I';
  } else if (!strncmp(p,"del",strlen("del"))) {
    *firstend = 'D';
  } else if (!strncmp(p,"donor",strlen("donor"))) {
    *firstend = '5';
  } else if (!strncmp(p,"acceptor",strlen("acceptor"))) {
    *firstend = '3';
  } else {
    fprintf(stderr,"Can't parse type %s\n",p);
    abort();
  }


  *donorpos = *acceptorpos = 0U;

  /* Check start */
  /* p = acceptor:1.00,dir:antisense,sub:0,score:2 */
  if (!strncmp(p,"donor",strlen("donor"))) {
    /* Donor */
    *donorpos = *firstpos;
    get_splicesite_info(&(*sensep),&donorprob,p);
  } else if (!strncmp(p,"acceptor",strlen("acceptor"))) {
    *acceptorpos = *firstpos;
    get_splicesite_info(&(*sensep),&acceptorprob,p);
  }


  /* Check end */
  while (*p != '\0' && (*p != '.' || *(p+1) != '.')) p++;
  p++; p++;			/* Skip ".." */

  if (!strncmp(p,"end",strlen("end"))) {
    *secondend = 'E';
#if 0
  } else if (!strncmp(p,"start",strlen("start"))) {
    *secondend = 'S';
#endif
  } else if (!strncmp(p,"amb",strlen("amb"))) {
    *secondend = 'A';
  } else if (!strncmp(p,"term",strlen("term"))) {
    *secondend = 'T';
  } else if (!strncmp(p,"ins",strlen("ins"))) {
    *secondend = 'I';
  } else if (!strncmp(p,"del",strlen("del"))) {
    *secondend = 'D';
  } else if (!strncmp(p,"donor",strlen("donor"))) {
    *secondend = '5';
  } else if (!strncmp(p,"acceptor",strlen("acceptor"))) {
    *secondend = '3';
  } else {
    fprintf(stderr,"Can't parse type %s\n",p);
    abort();
  }


  if (!strncmp(p,"donor",strlen("donor"))) {
    *donorpos = *secondpos;
    get_splicesite_info(&(*sensep),&donorprob,p);
  } else if (!strncmp(p,"acceptor",strlen("acceptor"))) {
    *acceptorpos = *secondpos;
    get_splicesite_info(&(*sensep),&acceptorprob,p);
  }


  if (*donorpos == 0U && *acceptorpos == 0U) {
    *truestrand = ' ';
  } else if (*sensep == true) {
    *truestrand = *strand;
  } else if (*strand == '+') {
    *truestrand = '-';
  } else if (*strand == '-') {
    *truestrand = '+';
  } else {
    abort();
  }

  return;
}
