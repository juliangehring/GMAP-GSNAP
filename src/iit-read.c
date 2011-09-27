static char rcsid[] = "$Id: iit-read.c,v 1.95 2007/09/26 16:41:54 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "iit-read.h"
#include "iitdef.h"

#ifdef WORDS_BIGENDIAN
#include "bigendian.h"
#else
#include "littleendian.h"
#endif

#include <stdlib.h>		/* For qsort */
#include <string.h>		/* For memset */
#include <strings.h>
#include <ctype.h>		/* For isspace */
#ifdef HAVE_UNISTD_H
#include <unistd.h>		/* For mmap on Linux */
#endif
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>		/* For open, fstat, and mmap */
#endif
/* Not sure why this was included
#include <sys/param.h>
*/
#ifdef HAVE_FCNTL_H
#include <fcntl.h>		/* For open */
#endif
#ifdef HAVE_SYS_STAT_H
#include <sys/stat.h>		/* For open and fstat */
#endif
#include <sys/mman.h>		/* For mmap and madvise */
#include <errno.h>		/* For perror */
#include "assert.h"
#include "except.h"
#include "mem.h"
#include "access.h"
#include "fopen.h"
#include "chrom.h"
#include "uintlist.h"
#include "intlist.h"

/* Integer interval tree. */

/*
 * n intervals;
 *   specified by their indices e[1..n]
 *   and endpoint-access function:
 *                low  (e[i])
 *                high (e[i])
 *        is_contained (x, e[i])
 *   eg:
 *        interval e[i]          ... "[" low (e[i]) "," high (e[i]) ")"
 *        is_contained (x, e[i]) ... (    (low (e[i]) <= x
 *                                    and (x < high (e[i]))
 */

/*--------------------------------------------------------------------------*/


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


#define T IIT_T

char *
IIT_name (T this) {
  return this->name;
}

int
IIT_version (T this) {
  return this->version;
}

int
IIT_nintervals (T this) {
  return this->nintervals;
}


int
IIT_ntypes (T this) {
  return this->ntypes;
}

int
IIT_nfields (T this) {
  return this->nfields;
}


/* Note: strings are not allocated */
char **
IIT_types (int *ntypes, T this, bool alphabetizep) {
  char **types, *typestring;
  int type, t;
  Chrom_T *chroms;

  *ntypes = this->ntypes;
  types = (char **) CALLOC(this->ntypes,sizeof(char *));

  if (alphabetizep == false) {
    for (t = 0; t < this->ntypes; t++) {
      types[t] = IIT_typestring(this,t);
    }

  } else {
    chroms = (Chrom_T *) CALLOC(this->ntypes,sizeof(Chrom_T));

    for (type = 0; type < this->ntypes; type++) {
      typestring = IIT_typestring(this,type);
      chroms[type] = Chrom_from_string(typestring);
    }
    qsort(chroms,this->ntypes,sizeof(Chrom_T),Chrom_compare);

    for (t = 0; t < this->ntypes; t++) {
      types[t] = Chrom_to_string(chroms[t]);
    }
  }

  return types;
}

char **
IIT_fields (int *nfields, T this) {
  char **fields, *fieldstring;
  int field, t;

  *nfields = this->nfields;
  fields = (char **) CALLOC(this->nfields,sizeof(char *));

  for (t = 0; t < this->nfields; t++) {
    fields[t] = IIT_fieldstring(this,t);
  }

  return fields;
}


unsigned int
IIT_length (T this, int index) {
  Interval_T interval;

  interval = &(this->intervals[index-1]);
  return Interval_length(interval);
}


unsigned int
IIT_totallength (T this) {
  unsigned int max = 0U;
  Interval_T interval;
  int i;

  for (i = 0; i < this->nintervals; i++) {
    interval = &(this->intervals[i]);
    if (Interval_high(interval) > max) {
      max = Interval_high(interval);
    }
  }
  /* Convert from zero-based coordinate */
  return max+1U;
}


Interval_T
IIT_interval (T this, int index) {
  assert(index <= this->nintervals);
  return &(this->intervals[index-1]); /* Convert to 0-based */
}

/* The iit file has a '\0' after each string, so functions know where
   it ends */
char *
IIT_typestring (T this, int type) {
  unsigned int start;

  start = this->typepointers[type];
  return &(this->typestrings[start]);
}

int
IIT_typeint (T this, char *typestring) {
  int i = 0;
  unsigned int start;

  while (i < this->ntypes) {
    start = this->typepointers[i];
    if (!strcmp(typestring,&(this->typestrings[start]))) {
      return i;
    }
    i++;
  }

  return -1;
}

char *
IIT_fieldstring (T this, int fieldint) {
  unsigned int start;

  start = this->fieldpointers[fieldint];
  return &(this->fieldstrings[start]);
}

int
IIT_fieldint (T this, char *fieldstring) {
  int i = 0;
  unsigned int start;

  while (i < this->nfields) {
    start = this->fieldpointers[i];
    if (!strcmp(fieldstring,&(this->fieldstrings[start]))) {
      return i;
    }
    i++;
  }

  return -1;
}

/* The iit file has a '\0' after each string, so functions know where
   it ends */
char *
IIT_label (T this, int index) {
  int recno;
  unsigned int start;

  recno = index - 1; /* Convert to 0-based */
  start = this->labelpointers[recno];
  return &(this->labels[start]);
}

static void
file_move_absolute (T this, unsigned int start) {
  off_t offset = this->offset + start*((off_t) sizeof(char));

  if (lseek(this->fd,offset,SEEK_SET) < 0) {
    perror("Error in gmap, file_move_absolute");
    exit(9);
  }
  return;
}

static char *
IIT_label_to_string (T this, int index, bool *allocp) {
  char *label;
  int recno;
  unsigned int start, end;

  recno = index - 1; /* Convert to 0-based */
  start = this->labelpointers[recno];
  if (this->access == FILEIO) {
    end = this->labelpointers[recno+1];
    label = (char *) CALLOC(end-start,sizeof(char));
#ifdef HAVE_PTHREAD
    pthread_mutex_lock(&this->read_mutex);
#endif
    file_move_absolute(this,start);
    read(this->fd,label,end-start);
#ifdef HAVE_PTHREAD
    pthread_mutex_unlock(&this->read_mutex);
#endif
    *allocp = true;
    return label;
  } else {
    *allocp = false;
    return &(this->labels[start]);
  }
}


/* The iit file has a '\0' after each string, so functions know where
   it ends */
char *
IIT_annotation (T this, int index, bool *allocp) {
  char *annotation;
  int recno;
  unsigned int start, end;

  recno = index - 1; /* Convert to 0-based */
  start = this->annotpointers[recno];
  if (this->access == FILEIO) {
    end = this->annotpointers[recno+1];
    annotation = (char *) CALLOC(end-start,sizeof(char));
#ifdef HAVE_PTHREAD
    pthread_mutex_lock(&this->read_mutex);
#endif
    file_move_absolute(this,start);
    read(this->fd,annotation,end-start);
#ifdef HAVE_PTHREAD
    pthread_mutex_unlock(&this->read_mutex);
#endif
    *allocp = true;
    return annotation;
  } else {
    *allocp = false;
    return &(this->annotations[start]);
  }
}

/* The iit file has a '\0' after each string, so functions know where
   it ends */
char
IIT_annotation_firstchar (T this, int index) {
  int recno;
  unsigned int start;
  char buffer[1];

  recno = index - 1; /* Convert to 0-based */
  start = this->annotpointers[recno];
  if (this->access == FILEIO) {
#ifdef HAVE_PTHREAD
    pthread_mutex_lock(&this->read_mutex);
#endif
    file_move_absolute(this,start);
    read(this->fd,buffer,1);
#ifdef HAVE_PTHREAD
    pthread_mutex_unlock(&this->read_mutex);
#endif
    return buffer[0];
  } else {
    return this->annotations[start];
  }
}

unsigned int
IIT_annotation_strlen (T this, int index) {
  int recno;
  unsigned int start, end;

  recno = index - 1; /* Convert to 0-based */
  start = this->annotpointers[recno];
  end = this->annotpointers[recno+1];

  /*
  if (strlen(&(this->annotations[start])) != (end - start - 1)) {
    printf("Problem with %s: %d != %u\n",
    &(this->labels[this->labelpointers[recno]]),strlen(&(this->annotations[start])),end-start-1);
    abort();
  } else {
    printf("Okay %s: %d == %u\n",
    &(this->labels[this->labelpointers[recno]]),strlen(&(this->annotations[start])),end-start-1);
  }
  */

  return (end - start - 1);	/* Subtract terminal '\0' */
}

/* Always allocated */
char *
IIT_fieldvalue (T this, int index, int fieldint) {
  char *fieldvalue, *annotation, *p, *q;
  int recno, fieldno = 0, fieldlen;
  unsigned int start, end;
  bool allocp;

  recno = index - 1; /* Convert to 0-based */
  start = this->annotpointers[recno];
  if (this->access == FILEIO) {
    end = this->annotpointers[recno+1];
    annotation = (char *) CALLOC(end-start,sizeof(char));
#ifdef HAVE_PTHREAD
    pthread_mutex_lock(&this->read_mutex);
#endif
    file_move_absolute(this,start);
    read(this->fd,annotation,end-start);
#ifdef HAVE_PTHREAD
    pthread_mutex_unlock(&this->read_mutex);
#endif
    allocp = true;
  } else {
    annotation = &(this->annotations[start]);
    allocp = false;
  }

  p = annotation;
  while (*p != '\0' && fieldno < fieldint) {
    if (*p == '\n') {
      fieldno++;
    }
    p++;
  }

  if (*p == '\0') {
    fieldvalue = (char *) CALLOC(1,sizeof(char));
    fieldvalue[0] = '\0';
  } else {
    q = p;
    while (*q != '\0' && *q != '\n') {
      q++;
    }
    fieldlen = (q - p)/sizeof(char);
    fieldvalue = (char *) CALLOC(fieldlen+1,sizeof(char));
    strncpy(fieldvalue,p,fieldlen);
  }

  if (allocp == true) {
    FREE(annotation);
  }

  return fieldvalue;
}


void
IIT_dump_typestrings (FILE *fp, T this) {
  int type;
  unsigned int start;

  for (type = 0; type < this->ntypes; type++) {
    start = this->typepointers[type];
    fprintf(fp,"%d\t%s\n",type,&(this->typestrings[start]));
  }
  return;
}

void
IIT_dump_fieldstrings (FILE *fp, T this) {
  int field;
  unsigned int start;

  for (field = 0; field < this->nfields; field++) {
    start = this->fieldpointers[field];
    fprintf(fp,"%d\t%s\n",field,&(this->fieldstrings[start]));
  }
  return;
}

void
IIT_dump_labels (FILE *fp, T this) {
  int i;
  unsigned int start;
  char *label;

  for (i = 0; i < this->nintervals; i++) {
    start = this->labelpointers[i];
    label = &(this->labels[start]);
    fprintf(fp,"%s ",label);
  }
  fprintf(fp,"\n");
  return;
}


void
IIT_dump (T this, bool annotationonlyp) {
  int i;
  Interval_T interval;
  char *annotation;
  bool allocp;

  for (i = 0; i < this->nintervals; i++) {
    if (annotationonlyp == false) {
      printf(">%s",IIT_label(this,i+1));

      interval = &(this->intervals[i]);
      if (Interval_sign(interval) < 0) {
	printf(" %u %u",Interval_high(interval),Interval_low(interval));
      } else {
	printf(" %u %u",Interval_low(interval),Interval_high(interval));
      }
      if (Interval_type(interval) > 0) {
	printf(" %s",IIT_typestring(this,Interval_type(interval)));
      }
      printf("\n");
    }

    annotation = IIT_annotation(this,i+1,&allocp);
    if (strlen(annotation) == 0) {
      /* Don't print anything */
    } else if (annotation[strlen(annotation)-1] == '\n') {
      printf("%s",annotation);
    } else {
      printf("%s\n",annotation);
    }
    if (allocp == true) {
      FREE(annotation);
    }
  }
  return;
}

void
IIT_dump_formatted (T this, bool directionalp) {
  int i;
  Interval_T interval;
  unsigned int start, startpos, endpos;
  char *label, firstchar;

  for (i = 0; i < this->nintervals; i++) {
    interval = &(this->intervals[i]);
    start = this->labelpointers[i];
    label = &(this->labels[start]);
    printf("%s\t",label);
    startpos = Interval_low(interval);
    endpos = startpos + Interval_length(interval) - 1U;

    if (directionalp == false) {
      printf("%u..%u\t",startpos+1U,endpos+1U);
    } else if (this->version <= 1) {
      firstchar = IIT_annotation_firstchar(this,i+1);
      if (firstchar == '-') {
	printf("%u..%u\t",endpos+1U,startpos+1U);
      } else {
	printf("%u..%u\t",startpos+1U,endpos+1U);
      }
    } else {
      if (Interval_sign(interval) < 0) {
	printf("%u..%u\t",endpos+1U,startpos+1U);
      } else {
	printf("%u..%u\t",startpos+1U,endpos+1U);
      }
    }

    printf("%u",Interval_length(interval));
    if (Interval_type(interval) > 0) {
      printf("\t%s",IIT_typestring(this,Interval_type(interval)));
    }
    printf("\n");
  }
  return;
}

static int
uint_cmp (const void *x, const void *y) {
  unsigned int a = * (unsigned int *) x;
  unsigned int b = * (unsigned int *) y;

  if (a < b) {
    return -1;
  } else if (a > b) {
    return +1;
  } else {
    return 0;
  }
}


unsigned int *
IIT_transitions (int **signs, int *nedges, T this) { 
  unsigned int *edges, *starts, *ends;
  int nintervals, i, j, k;
  Interval_T interval;
  Uintlist_T startlist = NULL, endlist = NULL;

  for (i = 0; i < this->nintervals; i++) {
    interval = &(this->intervals[i]);
    startlist = Uintlist_push(startlist,Interval_low(interval));
    endlist = Uintlist_push(endlist,Interval_high(interval));
  }

  if (Uintlist_length(startlist) == 0) {
    edges = (unsigned int *) NULL;
    *signs = (int *) NULL;
    *nedges = 0;
  } else {
    starts = Uintlist_to_array(&nintervals,startlist);
    ends = Uintlist_to_array(&nintervals,endlist);
    qsort(starts,nintervals,sizeof(unsigned int),uint_cmp);
    qsort(ends,nintervals,sizeof(unsigned int),uint_cmp);

    *nedges = nintervals+nintervals;
    *signs = (int *) CALLOC(*nedges,sizeof(int));
    edges = (unsigned int *) CALLOC(*nedges,sizeof(unsigned int));
    i = j = k = 0;
    while (i < nintervals && j < nintervals) {
      if (starts[i] <= ends[j]) {
	(*signs)[k] = +1;
	edges[k++] = starts[i++];
      } else {
	(*signs)[k] = -1;
	edges[k++] = ends[j++];
      }
    }
    while (i < nintervals) {
      (*signs)[k] = +1;
      edges[k++] = starts[i++];
    }
    while (j < nintervals) {
      (*signs)[k] = -1;
      edges[k++] = ends[j++];
    }
  }

  Uintlist_free(&endlist);
  Uintlist_free(&startlist);

  return edges;
}

static int *
sort_matches (T this, int *matches, int nmatches, bool alphabetizep) {
  int *sorted;
  int type, index, i, j, k = 0, t;
  List_T *intervallists;
  Interval_T *intervals, interval;
  int *matches1, nmatches1, nintervals;
  char *typestring;
  Chrom_T *chroms;

  if (nmatches == 0) {
    return (int *) NULL;
  } else {
    sorted = (int *) CALLOC(nmatches,sizeof(int));
  }
  
  intervallists = (List_T *) CALLOC(this->ntypes,sizeof(List_T));
  for (i = 0; i < nmatches; i++) {
    index = matches[i];
    interval = &(this->intervals[index-1]);
    type = Interval_type(interval);
    intervallists[type] = List_push(intervallists[type],(void *) interval);
  }

  if (alphabetizep == true) {
    chroms = (Chrom_T *) CALLOC(this->ntypes,sizeof(Chrom_T));

    for (type = 0; type < this->ntypes; type++) {
      typestring = IIT_typestring(this,type);
      chroms[type] = Chrom_from_string(typestring);
    }
    qsort(chroms,this->ntypes,sizeof(Chrom_T),Chrom_compare);
  }

  for (t = 0; t < this->ntypes; t++) {
    if (alphabetizep == false) {
      type = t;
      typestring = IIT_typestring(this,type);
    } else {
      typestring = Chrom_to_string(chroms[t]);
      type = IIT_typeint(this,typestring);
    }

    if ((nintervals = List_length(intervallists[type])) > 0) {
      intervals = (Interval_T *) List_to_array(intervallists[type],/*end*/NULL);
      qsort(intervals,nintervals,sizeof(Interval_T),Interval_cmp);

      i = 0;
      while (i < nintervals) {
	interval = intervals[i];
	matches1 = IIT_get_exact_multiple(&nmatches1,this,Interval_low(interval),Interval_high(interval),type);
	if (matches1 != NULL) {
	  for (j = 0; j < nmatches1; j++) {
	    sorted[k++] = matches1[j];
	  }
	  i += nmatches1;
	  FREE(matches1);
	}
      }

      FREE(intervals);
      List_free(&(intervallists[type]));
    }

  }

  if (alphabetizep == true) {
    for (t = 0; t < this->ntypes; t++) {
      Chrom_free(&(chroms[t]));
    }
    FREE(chroms);
  }

  FREE(intervallists);
  return sorted;
}



void
IIT_dump_counts (T this, bool alphabetizep) { 
  int type, index, i, j, k, t;
  Interval_T interval;
  Uintlist_T *startlists, *endlists;
  int *matches, nmatches, nintervals;
  unsigned int *starts, *ends, edge;
  char *typestring;
  Chrom_T *chroms;

  startlists = (Uintlist_T *) CALLOC(this->ntypes,sizeof(Uintlist_T));
  endlists = (Uintlist_T *) CALLOC(this->ntypes,sizeof(Uintlist_T));
  for (i = 0; i < this->nintervals; i++) {
    interval = &(this->intervals[i]);
    type = Interval_type(interval);
    startlists[type] = Uintlist_push(startlists[type],Interval_low(interval));
    endlists[type] = Uintlist_push(endlists[type],Interval_high(interval));
  }

  if (alphabetizep == true) {
    chroms = (Chrom_T *) CALLOC(this->ntypes,sizeof(Chrom_T));

    for (type = 0; type < this->ntypes; type++) {
      typestring = IIT_typestring(this,type);
      chroms[type] = Chrom_from_string(typestring);
    }
    qsort(chroms,this->ntypes,sizeof(Chrom_T),Chrom_compare);
  }

  for (t = 0; t < this->ntypes; t++) {
    if (alphabetizep == false) {
      type = t;
      typestring = IIT_typestring(this,type);
    } else {
      typestring = Chrom_to_string(chroms[t]);
      type = IIT_typeint(this,typestring);
    }

    if (Uintlist_length(startlists[type]) > 0) {
      starts = Uintlist_to_array(&nintervals,startlists[type]);
      ends = Uintlist_to_array(&nintervals,endlists[type]);
      qsort(starts,nintervals,sizeof(unsigned int),uint_cmp);
      qsort(ends,nintervals,sizeof(unsigned int),uint_cmp);

      i = j = 0;
      while (i < nintervals || j < nintervals) {
	if (i >= nintervals && j >= nintervals) {
	  /* done */
	  matches = (int *) NULL;
	} else if (i >= nintervals) {
	  /* work on remaining ends */
	  edge = ends[j++];
	  matches = IIT_get_typed(&nmatches,this,edge,edge,type,/*sortp*/false);
	  printf("%s\t%u\tend\t%d",typestring,edge,nmatches);
	  while (j < nintervals && ends[j] == edge) {
	    j++;
	  }
	} else if (j >= nintervals) {
	  /* work on remaining starts */
	  edge = starts[i++];
	  matches = IIT_get_typed(&nmatches,this,edge,edge,type,/*sortp*/false);
	  printf("%s\t%u\tstart\t%d",typestring,edge,nmatches);
	  while (i < nintervals && starts[i] == edge) {
	    i++;
	  }
	} else if (starts[i] <= ends[j]) {
	  edge = starts[i++];
	  matches = IIT_get_typed(&nmatches,this,edge,edge,type,/*sortp*/false);
	  printf("%s\t%u\tstart\t%d",typestring,edge,nmatches);
	  while (i < nintervals && starts[i] == edge) {
	    i++;
	  }
	} else {
	  edge = ends[j++];
	  matches = IIT_get_typed(&nmatches,this,edge,edge,type,/*sortp*/false);
	  printf("%s\t%u\tend\t%d",typestring,edge,nmatches);
	  while (j < nintervals && ends[j] == edge) {
	    j++;
	  }
	}

	if (matches != NULL) {
	  index = matches[0];
	  printf("\t%s",IIT_label(this,index));
	  for (k = 1; k < nmatches; k++) {
	    index = matches[k];
	    printf(",%s",IIT_label(this,index));
	  }
	  printf("\n");
	  FREE(matches);
	}
      }

      Uintlist_free(&(endlists[type]));
      Uintlist_free(&(startlists[type]));
      FREE(ends);
      FREE(starts);
    }

  }

  if (alphabetizep == true) {
    for (t = 0; t < this->ntypes; t++) {
      Chrom_free(&(chroms[t]));
    }
    FREE(chroms);
  }

  FREE(endlists);
  FREE(startlists);

  return;
}

/************************************************************************
 * For file format, see iit-write.c
 ************************************************************************/

void
IIT_free (T *old) {
  if (*old != NULL) {
    if ((*old)->name != NULL) {
      FREE((*old)->name);
    }

    if ((*old)->access == MMAPPED) {
      munmap((void *) (*old)->finfo,(*old)->flength);
      close((*old)->fd);
    } else if ((*old)->access == FILEIO) {
      close((*old)->fd);
    } else if ((*old)->access == ALLOCATED) {
      /* Nothing to close.  IIT must have been created by IIT_new. */
    } else {
      abort();
    }

    if ((*old)->access != ALLOCATED) {
      FREE((*old)->annotpointers);
      FREE((*old)->labels);
      FREE((*old)->labelpointers);
      FREE((*old)->labelorder);
      if ((*old)->fieldstrings != NULL) {
	FREE((*old)->fieldstrings);
      }
      FREE((*old)->fieldpointers);
      FREE((*old)->typestrings);
      FREE((*old)->typepointers);
    }

    FREE((*old)->intervals);
    FREE((*old)->nodes);
    FREE((*old)->omegas);
    FREE((*old)->sigmas);
    if ((*old)->alphas != NULL) {
      FREE((*old)->betas);
      FREE((*old)->alphas);
    }

    FREE(*old);

  }
  return;
}

/* This procedure leaks memory, but used only for debugging */
void
IIT_debug (char *filename) {
  T new;
  FILE *fp;
  off_t offset = 0, filesize, last_offset, stringlen;
  int i, j;
  size_t items_read;

  if ((fp = FOPEN_READ_BINARY(filename)) == NULL) {
    fprintf(stderr,"Can't open file %s\n",filename);
    return;
  } else {
    new = (T) MALLOC(sizeof(*new));
    new->version = 1;
  }

  filesize = Access_filesize(filename);

  if (FREAD_INT(&new->nintervals,fp) < 1) {
    return;
  }
  if (new->nintervals == 0) {
    /* New format */
    FREAD_INT(&new->version,fp);
    if (new->version > IIT_LATEST_VERSION) {
      fprintf(stderr,"This file is version %d, but this software can only read up to version %d\n",
	      new->version,IIT_LATEST_VERSION);
      return;
    }
    if (FREAD_INT(&new->nintervals,fp) < 1) {
      return;
    } else if ((offset += sizeof(int)+sizeof(int)) > filesize) {
      return;
    }
  }
  printf("version: %d\n",new->version);

  if (new->nintervals < 0) {
    return;
  } else if ((offset += sizeof(int)) > filesize) {
    return;
  } else {
    printf("nintervals: %d\n",new->nintervals);
  }

  if (FREAD_INT(&new->ntypes,fp) < 1) {
    return;
  } else if (new->ntypes < 0) {
    return;
  } else if ((offset += sizeof(int)) > filesize) {
    return;
  } else {
    printf("ntypes: %d\n",new->ntypes);
  }

  if (new->version < 2) {
    new->nfields = 0;
  } else {
    if (FREAD_INT(&new->nfields,fp) < 1) {
      return;
    } else if (new->nfields < 0) {
      return;
    } else if ((offset += sizeof(int)) > filesize) {
      return;
    } else {
      printf("nfields: %d\n",new->nfields);
    }
  }

  if (FREAD_INT(&new->nnodes,fp) < 1) {
    return;
  } else if (new->nnodes < 0) {
    return;
  } else if ((offset += sizeof(int)) > filesize) {
    return;
  } else {
    printf("nnodes: %d\n",new->nnodes);
  }

  if (new->version >= 2) {
    if ((offset += sizeof(int)*(new->nintervals+1)) > filesize) {
      return;
    } else {
      new->alphas = (int *) CALLOC(new->nintervals+1,sizeof(int));
      if ((items_read = FREAD_INTS(new->alphas,new->nintervals+1,fp)) != new->nintervals + 1) {
	return;
      }
    }

    if ((offset += sizeof(int)*(new->nintervals+1)) > filesize) {
      return;
    } else {
      new->betas = (int *) CALLOC(new->nintervals+1,sizeof(int));
      if ((items_read = FREAD_INTS(new->betas,new->nintervals+1,fp)) != new->nintervals + 1) {
	return;
      }
    }
  }

  if ((offset += sizeof(int)*(new->nintervals+1)) > filesize) {
    return;
  } else {
    new->sigmas = (int *) CALLOC(new->nintervals+1,sizeof(int));
    if ((items_read = FREAD_INTS(new->sigmas,new->nintervals+1,fp)) != new->nintervals + 1) {
      return;
    }
    printf("sigmas:");
    for (i = 0; i < new->nintervals+1; i++) {
      printf(" %d",new->sigmas[i]);
    }
    printf("\n");
  }

  if ((offset += sizeof(int)*(new->nintervals+1)) > filesize) {
    return;
  } else {
    new->omegas = (int *) CALLOC(new->nintervals+1,sizeof(int));
    if ((items_read = FREAD_INTS(new->omegas,new->nintervals+1,fp)) != new->nintervals + 1) {
      return;
    }
    printf("omegas:");
    for (i = 0; i < new->nintervals+1; i++) {
      printf(" %d",new->omegas[i]);
    }
    printf("\n");
  }

  new->nodes = (struct FNode_T *) CALLOC(new->nnodes,sizeof(struct FNode_T));
#ifdef WORDS_BIGENDIAN
  for (i = 0; i < new->nnodes; i++) {
    Bigendian_fread_uint(&(new->nodes[i].value),fp);
    Bigendian_fread_int(&(new->nodes[i].a),fp);
    Bigendian_fread_int(&(new->nodes[i].b),fp);
    Bigendian_fread_int(&(new->nodes[i].leftindex),fp);
    Bigendian_fread_int(&(new->nodes[i].rightindex),fp);
  }
   offset += (sizeof(unsigned int)+sizeof(int)+sizeof(int)+sizeof(int)+sizeof(int))*new->nnodes;
#else
  if (sizeof(struct FNode_T) == sizeof(unsigned int)+sizeof(int)+sizeof(int)+sizeof(int)+sizeof(int)) {
    offset += sizeof(struct FNode_T)*fread(new->nodes,sizeof(struct FNode_T),new->nnodes,fp);
  } else {
    for (i = 0; i < new->nnodes; i++) {
      fread(&(new->nodes[i].value),sizeof(unsigned int),1,fp);
      fread(&(new->nodes[i].a),sizeof(int),1,fp);
      fread(&(new->nodes[i].b),sizeof(int),1,fp);
      fread(&(new->nodes[i].leftindex),sizeof(int),1,fp);
      fread(&(new->nodes[i].rightindex),sizeof(int),1,fp);
    }
    offset += (sizeof(unsigned int)+sizeof(int)+sizeof(int)+sizeof(int)+sizeof(int))*new->nnodes;
  }
#endif
  printf("nodes:\n");
  for (i = 0; i < new->nnodes; i++) {
    printf(" value:%u a:%d b:%d lindex:%d rindex:%d\n",
	   new->nodes[i].value,new->nodes[i].a,new->nodes[i].b,
	   new->nodes[i].leftindex,new->nodes[i].rightindex);
  }
  printf("\n");

  new->intervals = (struct Interval_T *) CALLOC(new->nintervals,sizeof(struct Interval_T));
#ifdef WORDS_BIGENDIAN
  for (i = 0; i < new->nintervals; i++) {
    Bigendian_fread_uint(&(new->intervals[i].low),fp);
    Bigendian_fread_uint(&(new->intervals[i].high),fp);
    if (new->version >= 2) {
      Bigendian_fread_int(&(new->intervals[i].sign),fp);
    } else {
      new->intervals[i].sign = +1;
      /* First char of annotation occasionally used in version 1 to indicate reverse sign */
    }
    Bigendian_fread_int(&(new->intervals[i].type),fp);
  }
  if (new->version >= 2) {
    offset += (sizeof(unsigned int)+sizeof(unsigned int)+sizeof(int)+sizeof(int))*new->nintervals;
  } else {
    offset += (sizeof(unsigned int)+sizeof(unsigned int)+sizeof(int))*new->nintervals;
  }
#else
  if (new->version >= 2 && sizeof(struct Interval_T) == sizeof(unsigned int)+sizeof(unsigned int)+sizeof(int)+sizeof(int)) {
    offset += sizeof(struct Interval_T)*fread(new->intervals,sizeof(struct Interval_T),new->nintervals,fp);
  } else if (new->version <= 1 && sizeof(struct Interval_T) == sizeof(unsigned int)+sizeof(unsigned int)+sizeof(int)) {
    offset += sizeof(struct Interval_T)*fread(new->intervals,sizeof(struct Interval_T),new->nintervals,fp);
  } else {
    for (i = 0; i < new->nintervals; i++) {
      fread(&(new->intervals[i].low),sizeof(unsigned int),1,fp);
      fread(&(new->intervals[i].high),sizeof(unsigned int),1,fp);
      if (new->version >= 2) {
	fread(&(new->intervals[i].sign),sizeof(int),1,fp);
      } else {
	new->intervals[i].sign = +1;
	/* First char of annotation occasionally used in version 1 to indicate reverse sign */
      }
      fread(&(new->intervals[i].type),sizeof(int),1,fp);
    }
    if (new->version >= 2) {
      offset += (sizeof(unsigned int)+sizeof(unsigned int)+sizeof(int)+sizeof(int))*new->nintervals;
    } else {
      offset += (sizeof(unsigned int)+sizeof(unsigned int)+sizeof(int))*new->nintervals;
    }
  }
#endif
  printf("intervals:\n");
  for (i = 0; i < new->nintervals; i++) {
    printf(" low:%u high:%u",new->intervals[i].low,new->intervals[i].high);
    if (new->version >= 2) {
      if (new->intervals[i].sign > 0) {
	printf(" sign:+");
      } else if (new->intervals[i].sign < 0) {
	printf(" sign:-");
      } else {
	printf(" sign:0");
      }
    }
    printf(" type:%d\n",new->intervals[i].type);
  }
  printf("\n");


  if (new->version < 2) {
#if 0
    IIT_compute_flanking(new);
#endif
  } else {

    /* For debugging of flanking */
    printf("alpha order:\n");
    for (i = 0; i < new->nintervals+1; i++) {
      printf(" %d",j = new->alphas[i]);
#if 0
      if (j > 0) {
	printf(" [low:%u high:%u]",new->intervals[j-1].low,new->intervals[j-1].high);
      }
      printf("\n");
#endif
    }
    printf("\n");

    printf("beta order:\n");
    for (i = 0; i < new->nintervals+1; i++) {
      printf(" %d",j = new->betas[i]);
#if 0
      if (j > 0) {
	printf(" [low:%u high:%u]",new->intervals[j-1].low,new->intervals[j-1].high);
      }
      printf("\n");
#endif
    }
    printf("\n");
    /* End of debugging */
  }

  new->typepointers = (unsigned int *) CALLOC(new->ntypes+1,sizeof(unsigned int));
  offset += sizeof(int)*FREAD_UINTS(new->typepointers,new->ntypes+1,fp);
  printf("typepointers:");
  for (i = 0; i < new->ntypes+1; i++) {
    printf(" %u",new->typepointers[i]);
  }
  printf("\n");

  stringlen = new->typepointers[new->ntypes];
  printf("stringlen of typestrings = %d\n",stringlen);
  new->typestrings = (char *) CALLOC(stringlen,sizeof(char));
  offset += sizeof(char)*FREAD_CHARS(new->typestrings,stringlen,fp);
  printf("typestrings:\n");
  for (i = 0; i < stringlen; i++) {
    printf("%c",new->typestrings[i]);
  }
  printf("\n");


  new->fieldpointers = (unsigned int *) CALLOC(new->nfields+1,sizeof(unsigned int));
  if (new->version < 2) {
    new->fieldpointers[0] = '\0';
  } else {
    offset += sizeof(int)*FREAD_UINTS(new->fieldpointers,new->nfields+1,fp);
  }
  printf("fieldpointers:");
  for (i = 0; i < new->nfields+1; i++) {
    printf(" %u",new->fieldpointers[i]);
  }
  printf("\n");

  stringlen = new->fieldpointers[new->nfields];
  printf("stringlen of fieldstrings = %d\n",stringlen);
  if (stringlen == 0) {
    new->fieldstrings = (char *) NULL;
  } else {
    new->fieldstrings = (char *) CALLOC(stringlen,sizeof(char));
    offset += sizeof(char)*FREAD_CHARS(new->fieldstrings,stringlen,fp);
  }
  printf("fieldstrings:\n");
  for (i = 0; i < stringlen; i++) {
    printf("%c",new->fieldstrings[i]);
  }
  printf("\n");


  new->labelorder = (int *) CALLOC(new->nintervals,sizeof(int));
  offset += sizeof(int)*FREAD_INTS(new->labelorder,new->nintervals,fp);
  printf("labelorder:");
  for (i = 0; i < new->nintervals; i++) {
    printf(" %d",new->labelorder[i]);
  }
  printf("\n");

  new->labelpointers = (unsigned int *) CALLOC(new->nintervals+1,sizeof(unsigned int));
  offset += sizeof(int)*FREAD_UINTS(new->labelpointers,new->nintervals+1,fp);
  printf("labelpointers:");
  for (i = 0; i < new->nintervals+1; i++) {
    printf(" %u",new->labelpointers[i]);
  }
  printf("\n");
  
  stringlen = new->labelpointers[new->nintervals];
  printf("total stringlen of labels = %d\n",stringlen);
  new->labels = (char *) CALLOC(stringlen,sizeof(char));
  offset += sizeof(char)*FREAD_CHARS(new->labels,stringlen,fp);
  printf("labels:\n");
  for (i = 0; i < stringlen; i++) {
    printf("%c",new->labels[i]);
  }
  printf("\n");

  new->annotpointers = (unsigned int *) CALLOC(new->nintervals+1,sizeof(unsigned int));
  offset += sizeof(int)*FREAD_UINTS(new->annotpointers,new->nintervals+1,fp);
  printf("annotpointers:");
  for (i = 0; i < new->nintervals+1; i++) {
    printf(" %d",new->annotpointers[i]);
  }
  printf("\n");

  stringlen = new->annotpointers[new->nintervals];
  printf("total stringlen of annotations = %d\n",stringlen);
  last_offset = offset + sizeof(char)*stringlen;

  printf("filesize: %u bytes, final offset: %u bytes\n",filesize,last_offset);

  fclose(fp);
  return;
}


T
IIT_read (char *filename, char *name, bool readonlyp) {
  T new;
  FILE *fp;
  off_t offset = 0, filesize, last_offset, stringlen;
  int i, prot, oflag, mapflags;
  size_t items_read;

  if ((fp = FOPEN_READ_BINARY(filename)) == NULL) {
    fprintf(stderr,"Unable to read IIT file %s\n",filename);
    return NULL;
  } else {
    new = (T) MALLOC(sizeof(*new));
    new->version = 1;
  }

  filesize = Access_filesize(filename);

  if (name == NULL) {
    new->name = NULL;
  } else {
    new->name = (char *) CALLOC(strlen(name)+1,sizeof(char));
    strcpy(new->name,name);
  }

  if (FREAD_INT(&new->nintervals,fp) < 1) {
    fprintf(stderr,"IIT file %s appears to be empty\n",filename);
    return NULL;
  }
  if (new->nintervals == 0) {
    /* New format */
    FREAD_INT(&new->version,fp);
    if (new->version > IIT_LATEST_VERSION) {
      fprintf(stderr,"This file is version %d, but this software can only read up to version %d\n",
	      new->version,IIT_LATEST_VERSION);
      return NULL;
    }
    if (FREAD_INT(&new->nintervals,fp) < 1) {
      fprintf(stderr,"IIT file %s appears to be truncated\n",filename);
      return NULL;
    } else if ((offset += sizeof(int)+sizeof(int)) > filesize) {
      fprintf(stderr,"IIT file %s appears to have an offset that is too large\n",filename);
      return NULL;
    }
  }
  if (new->nintervals < 0) {
    fprintf(stderr,"IIT file %s appears to have a negative number of intervals\n",filename);
    return NULL;
  } else if ((offset += sizeof(int)) > filesize) {
    fprintf(stderr,"IIT file %s appears to have an offset that is too large\n",filename);
    return NULL;
  }

  if (FREAD_INT(&new->ntypes,fp) < 1) {
    fprintf(stderr,"IIT file %s appears to be truncated\n",filename);
    return NULL;
  } else if (new->ntypes < 0) {
    fprintf(stderr,"IIT file %s appears to have a negative number of types\n",filename);
    return NULL;
  } else if ((offset += sizeof(int)) > filesize) {
    fprintf(stderr,"IIT file %s appears to have an offset that is too large\n",filename);
    return NULL;
  }

  if (new->version < 2) {
    new->nfields = 0;
  } else {
    if (FREAD_INT(&new->nfields,fp) < 1) {
      fprintf(stderr,"IIT file %s appears to be truncated\n",filename);
      return NULL;
    } else if (new->nfields < 0) {
      fprintf(stderr,"IIT file %s appears to have a negative number of fields\n",filename);
      return NULL;
    } else if ((offset += sizeof(int)) > filesize) {
      fprintf(stderr,"IIT file %s appears to have an offset that is too large\n",filename);
      return NULL;
    }
  }

  if (FREAD_INT(&new->nnodes,fp) < 1) {
    fprintf(stderr,"IIT file %s appears to be truncated\n",filename);
    return NULL;
  } else if (new->nnodes < 0) {
    fprintf(stderr,"IIT file %s appears to have a negative number of nodes\n",filename);
    return NULL;
  } else if ((offset += sizeof(int)) > filesize) {
    fprintf(stderr,"IIT file %s appears to have an offset that is too large\n",filename);
    return NULL;
  }
    
  if (new->version >= 2) {
    if ((offset += sizeof(int)*(new->nintervals+1)) > filesize) {
      fprintf(stderr,"IIT file %s appears to have an offset that is too large\n",filename);
      return NULL;
    } else {
      new->alphas = (int *) CALLOC(new->nintervals+1,sizeof(int));
      if ((items_read = FREAD_INTS(new->alphas,new->nintervals+1,fp)) != new->nintervals + 1) {
	fprintf(stderr,"IIT file %s appears to be truncated\n",filename);
	return NULL;
      }
    }

    if ((offset += sizeof(int)*(new->nintervals+1)) > filesize) {
      fprintf(stderr,"IIT file %s appears to have an offset that is too large\n",filename);
      return NULL;
    } else {
      new->betas = (int *) CALLOC(new->nintervals+1,sizeof(int));
      if ((items_read = FREAD_INTS(new->betas,new->nintervals+1,fp)) != new->nintervals + 1) {
	fprintf(stderr,"IIT file %s appears to be truncated\n",filename);
	return NULL;
      }
    }
  }

  if ((offset += sizeof(int)*(new->nintervals+1)) > filesize) {
    fprintf(stderr,"IIT file %s appears to have an offset that is too large\n",filename);
    return NULL;
  } else {
    new->sigmas = (int *) CALLOC(new->nintervals+1,sizeof(int));
    if ((items_read = FREAD_INTS(new->sigmas,new->nintervals+1,fp)) != new->nintervals + 1) {
      fprintf(stderr,"IIT file %s appears to be truncated\n",filename);
      return NULL;
    }
  }

  if ((offset += sizeof(int)*(new->nintervals+1)) > filesize) {
    fprintf(stderr,"IIT file %s appears to have an offset that is too large\n",filename);
    return NULL;
  } else {
    new->omegas = (int *) CALLOC(new->nintervals+1,sizeof(int));
    if ((items_read = FREAD_INTS(new->omegas,new->nintervals+1,fp)) != new->nintervals + 1) {
      fprintf(stderr,"IIT file %s appears to be truncated\n",filename);
      return NULL;
    }
  }

  new->nodes = (struct FNode_T *) CALLOC(new->nnodes,sizeof(struct FNode_T));
#ifdef WORDS_BIGENDIAN
  for (i = 0; i < new->nnodes; i++) {
    Bigendian_fread_uint(&(new->nodes[i].value),fp);
    Bigendian_fread_int(&(new->nodes[i].a),fp);
    Bigendian_fread_int(&(new->nodes[i].b),fp);
    Bigendian_fread_int(&(new->nodes[i].leftindex),fp);
    Bigendian_fread_int(&(new->nodes[i].rightindex),fp);
  }
  offset += (sizeof(unsigned int)+sizeof(int)+sizeof(int)+sizeof(int)+sizeof(int))*new->nnodes;
#else
  if (sizeof(struct FNode_T) == sizeof(unsigned int)+sizeof(int)+sizeof(int)+sizeof(int)+sizeof(int)) {
    offset += sizeof(struct FNode_T)*fread(new->nodes,sizeof(struct FNode_T),new->nnodes,fp);
  } else {
    for (i = 0; i < new->nnodes; i++) {
      fread(&(new->nodes[i].value),sizeof(unsigned int),1,fp);
      fread(&(new->nodes[i].a),sizeof(int),1,fp);
      fread(&(new->nodes[i].b),sizeof(int),1,fp);
      fread(&(new->nodes[i].leftindex),sizeof(int),1,fp);
      fread(&(new->nodes[i].rightindex),sizeof(int),1,fp);
    }
    offset += (sizeof(unsigned int)+sizeof(int)+sizeof(int)+sizeof(int)+sizeof(int))*new->nnodes;
  }
#endif
  if (offset > filesize) {
    fprintf(stderr,"IIT file %s appears to have an offset that is too large\n",filename);
    return NULL;
  }

  new->intervals = (struct Interval_T *) CALLOC(new->nintervals,sizeof(struct Interval_T));
#ifdef WORDS_BIGENDIAN
  for (i = 0; i < new->nintervals; i++) {
    Bigendian_fread_uint(&(new->intervals[i].low),fp);
    Bigendian_fread_uint(&(new->intervals[i].high),fp);
    if (new->version >= 2) {
      Bigendian_fread_int(&(new->intervals[i].sign),fp);
    } else {
      new->intervals[i].sign = +1;
    }
    Bigendian_fread_int(&(new->intervals[i].type),fp);
  }
  if (new->version >= 2) {
    offset += (sizeof(unsigned int)+sizeof(unsigned int)+sizeof(int)+sizeof(int))*new->nintervals;
  } else {
    offset += (sizeof(unsigned int)+sizeof(unsigned int)+sizeof(int))*new->nintervals;
  }
#else
  if (new->version >= 2 && sizeof(struct Interval_T) == sizeof(unsigned int)+sizeof(unsigned int)+sizeof(int)+sizeof(int)) {
    offset += sizeof(struct Interval_T)*fread(new->intervals,sizeof(struct Interval_T),new->nintervals,fp);
  } else if (new->version <= 1 && sizeof(struct Interval_T) == sizeof(unsigned int)+sizeof(unsigned int)+sizeof(int)) {
    offset += sizeof(struct Interval_T)*fread(new->intervals,sizeof(struct Interval_T),new->nintervals,fp);
  } else {
    for (i = 0; i < new->nintervals; i++) {
      fread(&(new->intervals[i].low),sizeof(unsigned int),1,fp);
      fread(&(new->intervals[i].high),sizeof(unsigned int),1,fp);
      if (new->version >= 2) {
	fread(&(new->intervals[i].sign),sizeof(int),1,fp);
      } else {
	new->intervals[i].sign = +1;
      }
      fread(&(new->intervals[i].type),sizeof(int),1,fp);
    }
    if (new->version >= 2) {
      offset += (sizeof(unsigned int)+sizeof(unsigned int)+sizeof(int)+sizeof(int))*new->nintervals;
    } else {
      offset += (sizeof(unsigned int)+sizeof(unsigned int)+sizeof(int))*new->nintervals;
    }
  }
#endif
  if (offset > filesize) {
    fprintf(stderr,"IIT file %s appears to have an offset that is too large\n",filename);
    return NULL;
  }

  if (new->version < 2) {
#if 0
    /* Computing only if needed */
    compute_flanking(new);
#else
    new->alphas = new->betas = (int *) NULL;
#endif
  }

  new->typepointers = (unsigned int *) CALLOC(new->ntypes+1,sizeof(unsigned int));
  offset += sizeof(int)*FREAD_UINTS(new->typepointers,new->ntypes+1,fp);

  stringlen = new->typepointers[new->ntypes];
  new->typestrings = (char *) CALLOC(stringlen,sizeof(char));
  offset += sizeof(char)*FREAD_CHARS(new->typestrings,stringlen,fp);


  new->fieldpointers = (unsigned int *) CALLOC(new->nfields+1,sizeof(unsigned int));
  if (new->version < 2) {
    new->fieldpointers[0] = '\0';
  } else {
    offset += sizeof(int)*FREAD_UINTS(new->fieldpointers,new->nfields+1,fp);
  }
  stringlen = new->fieldpointers[new->nfields];
  if (stringlen == 0) {
    new->fieldstrings = (char *) NULL;
  } else {
    new->fieldstrings = (char *) CALLOC(stringlen,sizeof(char));
    offset += sizeof(char)*FREAD_CHARS(new->fieldstrings,stringlen,fp);
  }


  new->labelorder = (int *) CALLOC(new->nintervals,sizeof(int));
  offset += sizeof(int)*FREAD_INTS(new->labelorder,new->nintervals,fp);

  new->labelpointers = (unsigned int *) CALLOC(new->nintervals+1,sizeof(unsigned int));
  offset += sizeof(int)*FREAD_UINTS(new->labelpointers,new->nintervals+1,fp);
  
  stringlen = new->labelpointers[new->nintervals];
  new->labels = (char *) CALLOC(stringlen,sizeof(char));
  offset += sizeof(char)*FREAD_CHARS(new->labels,stringlen,fp);

  new->annotpointers = (unsigned int *) CALLOC(new->nintervals+1,sizeof(unsigned int));
  offset += sizeof(int)*FREAD_UINTS(new->annotpointers,new->nintervals+1,fp);

  stringlen = new->annotpointers[new->nintervals];
  last_offset = offset + sizeof(char)*stringlen;

  /* For some reason, the check 'last_offset != filesize' fails */
  if (last_offset != filesize) {
    fprintf(stderr,"Problem with last_offset (%u) not equal to filesize = (%u)\n",last_offset,filesize);
    return NULL;
  }

  fclose(fp);


#ifndef HAVE_MMAP
  new->annotations = (char *) NULL;
  new->fd = Access_fileio(filename);
  new->access = FILEIO;
  new->offset = offset;
#else  
  if (readonlyp == true) {
    new->finfo = (char *) Access_mmap(&new->fd,&new->flength,filename,sizeof(char),/*randomp*/true);
  } else {
    new->finfo = (char *) Access_mmap_rw(&new->fd,&new->flength,filename,sizeof(char),/*randomp*/true);
  }
  if (new->finfo == NULL) {
    new->annotations = (char *) NULL;
    /* fd already assigned */
    new->access = FILEIO;
    new->offset = offset;
  } else {
    new->annotations = (char *) &(new->finfo[offset]);
    new->access = MMAPPED;
    new->offset = 0U;		/* Not used */
  }
#endif

#ifdef HAVE_PTHREAD
  if (new->access == FILEIO) {
    pthread_mutex_init(&new->read_mutex,NULL);
  }
#endif

  return new;
}


/************************************************************************/

static void 
fnode_query_aux (int *min, int *max, T this, int nodeindex, unsigned int x) {
  int lambda;
  FNode_T node;

  if (nodeindex == -1) {
    return;
  }

  node = &(this->nodes[nodeindex]);
  if (x == node->value) {
    debug(printf("%uD:\n",node->value));
    if (node->a < *min) {
      *min = node->a;
    }
    if (node->b > *max) {
      *max = node->b;
    }
    return;
  } else if (x < node->value) {
    fnode_query_aux(&(*min),&(*max),this,node->leftindex,x);
    debug(printf("%uL:\n",node->value));
    if (node->a < *min) {
      *min = node->a;
    }
    for (lambda = node->a; lambda <= node->b; lambda++) {
      debug(printf("Looking at lambda %d, segment %d\n",
		   lambda,this->sigmas[lambda]));
      if (Interval_is_contained(x,this->intervals,this->sigmas[lambda]) == true) {
	if (lambda > *max) {
	  *max = lambda;
	}
      } else {
	return;
      }
    }
    return;
  } else { 
    /* (node->value < x) */
    fnode_query_aux(&(*min),&(*max),this,node->rightindex,x);
    debug(printf("%uR:\n", node->value));
    if (node->b > *max) {
      *max = node->b;
    }
    for (lambda = node->b; lambda >= node->a; lambda--) {
      debug(printf("Looking at lambda %d, segment %d\n",
		   lambda,this->omegas[lambda]));
      if (Interval_is_contained(x,this->intervals,this->omegas[lambda]) == true) {
	if (lambda < *min) {
	  *min = lambda;
	}
      } else {
	return;
      }
    }
    return;
  }
}

/************************************************************************/

int *
IIT_find (int *nmatches, T this, char *label) {
  int *matches = NULL, j;
  int low, middle, high, recno;
  bool foundp = false;
  int cmp;

  low = 0;
  high = this->nintervals;
  *nmatches = 0;

  debug(
	for (middle = low; middle < high; middle++) {
	  printf("%d:%d:%s, ",middle,this->labelorder[middle],
		 &(this->labels[this->labelpointers[this->labelorder[middle]]]));
	}
	printf("\n");
	);

  while (!foundp && low < high) {
    middle = (low+high)/2;
    debug(printf("low %d:%d:%s. middle %d:%d:%s high %d:%d:%s\n",
		 low,this->labelorder[low],
		 &(this->labels[this->labelpointers[this->labelorder[low]]]),
		 middle,this->labelorder[middle],
		 &(this->labels[this->labelpointers[this->labelorder[middle]]]),
		 high,this->labelorder[high],
		 &(this->labels[this->labelpointers[this->labelorder[high]]])));

    cmp = strcmp(label,&(this->labels[this->labelpointers[this->labelorder[middle]]]));
    if (cmp < 0) {
      high = middle;
    } else if (cmp > 0) {
      low = middle + 1;
    } else {
      foundp = true;
    }
  }

  if (foundp == true) {
    low = middle;
    while (low-1 >= 0 && 
	   !strcmp(label,&(this->labels[this->labelpointers[this->labelorder[low-1]]]))) {
      low--;
    }

    high = middle;
    while (high+1 < this->nintervals && 
	   !strcmp(label,&(this->labels[this->labelpointers[this->labelorder[high+1]]]))) {
      high++;
    }


    *nmatches = high - low + 1;
    if (*nmatches > 0) {
      matches = (int *) CALLOC(*nmatches,sizeof(int));
      j = 0;
      for (recno = low; recno <= high; recno++) {
	debug(printf("Pushing %d:%d\n",recno,this->labelorder[recno]));
	matches[j++] = this->labelorder[recno]+1;
      }
    }
  }

  return matches;
}

/* Slow.  Used before binary search method above. */
int
IIT_find_linear (T this, char *label) {
  int i;
  char *p;

  for (i = 0; i < this->nintervals; i++) {
    p = &(this->labels[this->labelpointers[i]]);
    while (isspace((int) *p)) {
      p++;
    }
    if (!strcmp(label,p)) {
      return i + 1;
    }
  }
  return -1;
}

int *
IIT_find_multiple (int *nmatches, T this, char **labels, int nlabels) {
  int *matches, *matches1, *sorted1;
  Intlist_T matchlist = NULL;
  int nmatches1, i, j;

  for (i = 0; i < nlabels; i++) {
    matches1 = IIT_find(&nmatches1,this,labels[i]);
    sorted1 = sort_matches(this,matches1,nmatches1,/*alphabetizep*/true);
    for (j = 0; j < nmatches1; j++) {
      matchlist = Intlist_push(matchlist,sorted1[j]);
    }
    FREE(sorted1);
    FREE(matches1);
  }
  matchlist = Intlist_reverse(matchlist);
  matches = Intlist_to_array(&(*nmatches),matchlist);
  Intlist_free(&matchlist);
  return matches;
}

int
IIT_find_one (T this, char *label) {
  int index;
  int *matches, nmatches;

  matches = IIT_find(&nmatches,this,label);
  if (nmatches == 0) {
    /*
    fprintf(stderr,"Expected one match for %s, but got 0\n",
	    label);
    */
    index = -1;
  } else {
    if (nmatches > 1) {
      fprintf(stderr,"Expected one match for %s, but got %d\n",
	      label,nmatches);
    }
    index = matches[0];
    FREE(matches);
  }

  return index;
}    


/************************************************************************/


static int
int_compare (const void *a, const void *b) {
  int x = * (int *) a;
  int y = * (int *) b;

  if (x < y) {
    return -1;
  } else if (y < x) {
    return 1;
  } else {
    return 0;
  }
}

/* This procedure not called currently */
int *
IIT_get_all (int *nmatches, T this, bool sortp) {
  int *sorted, *matches;
  int lambda, i, j = 0;

  matches = (int *) CALLOC(this->nintervals,sizeof(int));
  *nmatches = this->nintervals;
  if (sortp == true) {
    if (this->alphas == NULL) {
#if 0
      IIT_compute_flanking(this);
#else
      fprintf(stderr,"Flanking hits not supported on version %d of iit files.  Please use iit_update to update your file\n",
	      this->version);
      exit(9);
#endif
    }
    for (i = 1; i <= this->nintervals; i++) {
      matches[j++] = this->alphas[i];
    }
  } else {
    for (i = 1; i <= this->nintervals; i++) {
      matches[j++] = i;
    }
  }

  if (sortp == false) {
    return matches;
  } else if (this->nfields == 1) {
    return matches;
  } else {
    /* Fixes the tags */
    sorted = sort_matches(this,matches,*nmatches,/*alphabetizep*/true);
    FREE(matches);
    return sorted;
  }
}

/* This procedure not called currently */
int *
IIT_get_all_typed (int *nmatches, T this, int *types, int ntypes, bool sortp) {
  int *sorted, *matches;
  Intlist_T matchlist = NULL;
  Interval_T interval;
  int i, k;

  for (i = 1; i <= this->nintervals; i++) {
    interval = &(this->intervals[i-1]);
    k = 0;
    while (k < ntypes && Interval_type(interval) != types[k]) {
      k++;
    }
    if (k < ntypes) {
      matchlist = Intlist_push(matchlist,i);
    }
  }

  if (sortp == false) {
    matchlist = Intlist_reverse(matchlist);
    matches = Intlist_to_array(&(*nmatches),matchlist);
    Intlist_free(&matchlist);
    return matches;
  } else {
    matches = Intlist_to_array(&(*nmatches),matchlist);
    Intlist_free(&matchlist);
    sorted = sort_matches(this,matches,*nmatches,/*alphabetizep*/true);
    FREE(matches);
    return sorted;
  }
}


int *
IIT_get (int *nmatches, T this, unsigned int x, unsigned int y, bool sortp) {
  int *sorted, *matches = NULL, *uniq, neval, nuniq, i, j;
  int lambda, prev;
  int min1 = this->nintervals+1, max1 = 0, min2 = this->nintervals+1, max2 = 0;

  debug(printf("Entering IIT_get with query %u %u\n",x,y));
  fnode_query_aux(&min1,&max1,this,0,x);
  fnode_query_aux(&min2,&max2,this,0,y);
  debug(printf("min1=%d max1=%d  min2=%d max2=%d\n",min1,max1,min2,max2));

  *nmatches = 0;
  if (max2 >= min1) {
    neval = (max2 - min1 + 1) + (max2 - min1 + 1);
    matches = (int *) CALLOC(neval,sizeof(int));
    uniq = (int *) CALLOC(neval,sizeof(int));

    i = 0;
    for (lambda = min1; lambda <= max2; lambda++) {
      matches[i++] = this->sigmas[lambda];
      matches[i++] = this->omegas[lambda];
    }

    /* Eliminate duplicates */
    qsort(matches,neval,sizeof(int),int_compare);
    nuniq = 0;
    prev = 0;
    debug(printf("unique segments in lambda %d to %d:",min1,max2));
    for (i = 0; i < neval; i++) {
      if (matches[i] != prev) {
	debug(printf(" %d",matches[i]));
	uniq[nuniq++] = matches[i];
	prev = matches[i];
      }
    }
    debug(printf("\n"));

    for (i = 0; i < nuniq; i++) {
      if (Interval_overlap_p(x,y,this->intervals,uniq[i]) == true) {
	matches[(*nmatches)++] = uniq[i];
	debug(printf("Pushing overlapping segment %d (%u..%u)\n",uniq[i],
		     Interval_low(&(this->intervals[uniq[i]-1])),
		     Interval_high(&(this->intervals[uniq[i]-1]))));
      } else {
	debug(printf("Not pushing non-overlapping segment %d (%u..%u)\n",uniq[i],
		     Interval_low(&(this->intervals[uniq[i]-1])),
		     Interval_high(&(this->intervals[uniq[i]-1]))));
      }
    }

    FREE(uniq);
  }

  if (sortp == false) {
    return matches;
  } else {
    sorted = sort_matches(this,matches,*nmatches,/*alphabetizep*/true);
    FREE(matches);
    return sorted;
  }
}

static int
coord_search_low (T this, unsigned int x) {
  int low, middle, high;
  bool foundp = false;
  unsigned int middlevalue;

  low = 0;
  high = this->nintervals;

  while (!foundp && low < high) {
    middle = (low+high)/2;
    middlevalue = Interval_low(&(this->intervals[this->alphas[middle]-1]));
    if (x < middlevalue) {
      high = middle;
    } else if (x > middlevalue) {
      low = middle + 1;
    } else {
      foundp = true;
    }
  }

  if (foundp == true) {
    return middle;
  } else {
    return low;
  }
}

static int
coord_search_high (T this, unsigned int x) {
  int low, middle, high;
  bool foundp = false;
  unsigned int middlevalue;

  low = 0;
  high = this->nintervals;

  while (!foundp && low < high) {
    middle = (low+high)/2;
    middlevalue = Interval_high(&(this->intervals[this->betas[middle]-1]));
    if (x < middlevalue) {
      high = middle;
    } else if (x > middlevalue) {
      low = middle + 1;
    } else {
      foundp = true;
    }
  }

  if (foundp == true) {
    return middle;
  } else {
    return high;
  }
}


void
IIT_get_flanking (int **leftflanks, int *nleftflanks, int **rightflanks, int *nrightflanks,
		  T this, unsigned int x, unsigned int y, int nflanking, int sign) {
  int lambda;
  Interval_T interval;
  bool stopp;

  debug(printf("Entering IIT_get_flanking with query %u %u, nflanking = %d, sign %d\n",x,y,nflanking,sign));

  if (this->alphas == NULL) {
#if 0
    compute_flanking(this);
#else
    fprintf(stderr,"Flanking hits not supported on version %d of iit files.  Please use iit_update to update your file\n",
	    this->version);
    exit(9);
#endif
  }

  /* Look at alphas for right flank */
  lambda = coord_search_low(this,y);

  *rightflanks = (int *) CALLOC(nflanking,sizeof(int));
  *nrightflanks = 0;
  stopp = false;
  while (lambda <= this->nintervals && stopp == false) {
    interval = &(this->intervals[this->alphas[lambda]-1]);
    if (Interval_low(interval) <= y) {
      lambda++;
    } else if (sign != 0 && Interval_sign(interval) != sign) {
      lambda++;
    } else {
      (*rightflanks)[(*nrightflanks)++] = this->alphas[lambda];
      if (*nrightflanks < nflanking) {
	lambda++;
      } else {
	stopp = true;
      }
    }
  }

  /* Look at betas for left flank */
  lambda = coord_search_high(this,x);

  *leftflanks = (int *) CALLOC(nflanking,sizeof(int));
  *nleftflanks = 0;
  stopp = false;
  while (lambda >= 1 && stopp == false) {
    interval = &(this->intervals[this->betas[lambda]-1]);
    if (Interval_high(interval) >= x) {
      lambda--;
    } else if (sign != 0 && Interval_sign(interval) != sign) {
      lambda--;
    } else {
      (*leftflanks)[(*nleftflanks)++] = this->betas[lambda];
      if (*nleftflanks < nflanking) {
	lambda--;
      } else {
	stopp = true;
      }
    }
  }

  return;
}

void
IIT_get_flanking_typed (int **leftflanks, int *nleftflanks, int **rightflanks, int *nrightflanks,
			T this, unsigned int x, unsigned int y, int nflanking, int type) {
  int lambda;
  Interval_T interval;
  bool stopp;

  debug(printf("Entering IIT_get_flanking_typed with query %u %u\n",x,y));

  if (this->alphas == NULL) {
#if 0
    IIT_compute_flanking(this);
#else
    fprintf(stderr,"Flanking hits not supported on version %d of iit files.  Please use iit_update to update your file\n",
	    this->version);
    exit(9);
#endif
  }

  /* Look at alphas for right flank */
  lambda = coord_search_low(this,y);

  *rightflanks = (int *) CALLOC(nflanking,sizeof(int));
  *nrightflanks = 0;
  stopp = false;
  while (lambda <= this->nintervals && stopp == false) {
    interval = &(this->intervals[this->alphas[lambda]-1]);
    if (Interval_low(interval) <= y) {
      lambda++;
    } else if (Interval_type(interval) != type) {
      lambda++;
    } else {
      (*rightflanks)[(*nrightflanks)++] = this->alphas[lambda];
      if (*nrightflanks < nflanking) {
	lambda++;
      } else {
	stopp = true;
      }
    }
  }

  /* Look at betas for left flank */
  lambda = coord_search_high(this,x);

  *leftflanks = (int *) CALLOC(nflanking,sizeof(int));
  *nleftflanks = 0;
  stopp = false;
  while (lambda >= 1 && stopp == false) {
    interval = &(this->intervals[this->betas[lambda]-1]);
    if (Interval_high(interval) >= x) {
      lambda--;
    } else if (Interval_type(interval) != type) {
      lambda--;
    } else {
      (*leftflanks)[(*nleftflanks)++] = this->betas[lambda];
      if (*nleftflanks < nflanking) {
	lambda--;
      } else {
	stopp = true;
      }
    }
  }

  return;
}

void
IIT_get_flanking_multiple_typed (int **leftflanks, int *nleftflanks, int **rightflanks, int *nrightflanks,
				 T this, unsigned int x, unsigned int y, int nflanking, int *types, int ntypes) {
  int k;
  int lambda;
  Interval_T interval;
  bool stopp;

  debug(printf("Entering IIT_get_flanking_multiple_typed with query %u %u\n",x,y));

  if (this->alphas == NULL) {
#if 0
    IIT_compute_flanking(this);
#else
    fprintf(stderr,"Flanking hits not supported on version %d of iit files.  Please use iit_update to update your file\n",
	    this->version);
    exit(9);
#endif
  }

  /* Look at alphas for right flank */
  lambda = coord_search_low(this,y);

  *rightflanks = (int *) CALLOC(nflanking,sizeof(int));
  *nrightflanks = 0;
  stopp = false;
  while (lambda <= this->nintervals && stopp == false) {
    interval = &(this->intervals[this->alphas[lambda]-1]);
    if (Interval_low(interval) <= y) {
      lambda++;
    } else {
      k = 0;
      while (k < ntypes && Interval_type(interval) != types[k]) {
	k++;
      }
      if (k >= ntypes) {
	lambda++;
      } else {
	(*rightflanks)[(*nrightflanks)++] = this->alphas[lambda];
	if (*nrightflanks < nflanking) {
	  lambda++;
	} else {
	  stopp = true;
	}
      }
    }
  }


  /* Look at betas for left flank */
  lambda = coord_search_high(this,x);

  *leftflanks = (int *) CALLOC(nflanking,sizeof(int));
  *nleftflanks = 0;
  stopp = false;
  while (lambda >= 1 && stopp == false) {
    interval = &(this->intervals[this->betas[lambda]-1]);
    if (Interval_high(interval) >= x) {
      lambda--;
    } else {
      k = 0;
      while (k < ntypes && Interval_type(interval) != types[k]) {
	k++;
      }
      if (k >= ntypes) {
	lambda--;
      } else {
	(*leftflanks)[(*nleftflanks)++] = this->betas[lambda];
	if (*nleftflanks < nflanking) {
	  lambda--;
	} else {
	  stopp = true;
	}
      }
    }
  }

  return;
}


static const Except_T iit_error = { "IIT problem" };

int
IIT_get_one (T this, unsigned int x, unsigned int y) {
  int lambda;
  int min1 = this->nintervals+1, max1 = 0, min2 = this->nintervals+1, max2 = 0;

  debug(printf("Entering IIT_get_one with query %u %u\n",x,y));
  fnode_query_aux(&min1,&max1,this,0,x);
  fnode_query_aux(&min2,&max2,this,0,y);
  debug(printf("min1=%d max1=%d  min2=%d max2=%d\n",min1,max1,min2,max2));

  if (max2 >= min1) {
    for (lambda = min1; lambda <= max2; lambda++) {
      if (Interval_overlap_p(x,y,this->intervals,this->sigmas[lambda]) == true) {
	return this->sigmas[lambda];
      }
    }
    for (lambda = min1; lambda <= max2; lambda++) {
      if (Interval_overlap_p(x,y,this->intervals,this->omegas[lambda]) == true) {
	return this->omegas[lambda];
      }
    }
  }

  fprintf(stderr,"Expected one match for %u--%u, but got none\n",x,y);

  return -1;
}

/* Generally called where intervals don't overlap, like chromosomes,
   and where x == y. */
/*
int
IIT_get_one_safe (T this, unsigned int x, unsigned int y) {
  int index;
  int *matches, nmatches;

  matches = IIT_get(&nmatches,this,x,y,sortp);
  if (nmatches != 1) {
    fprintf(stderr,"Expected one match for %u--%u, but got %d\n",
	    x,y,nmatches);
    abort();
  }
  index = matches[0];
  FREE(matches);
  return index;
}
*/

int *
IIT_get_typed (int *ntypematches, T this, unsigned int x, unsigned int y, int type, bool sortp) {
  int *sorted;
  int index;
  int *typematches = NULL, *matches, nmatches, i, j;
  Interval_T interval;

  *ntypematches = 0;
  matches = IIT_get(&nmatches,this,x,y,/*sortp*/false);
  for (i = 0; i < nmatches; i++) {
    index = matches[i];
    interval = &(this->intervals[index-1]);
    if (Interval_type(interval) == type) {
      (*ntypematches)++;
    }
  }

  if (*ntypematches > 0) {
    typematches = (int *) CALLOC(*ntypematches,sizeof(int));
    j = 0;
    for (i = 0; i < nmatches; i++) {
      index = matches[i];
      interval = &(this->intervals[index-1]);
      if (Interval_type(interval) == type) {
	typematches[j++] = index;
      }
    }
  }
  
  if (matches != NULL) {
    FREE(matches);
  }

  if (sortp == false) {
    return typematches;
  } else {
    sorted = sort_matches(this,typematches,*ntypematches,/*alphabetizep*/false);
    FREE(typematches);
    return sorted;
  }
}

int *
IIT_get_multiple_typed (int *ntypematches, T this, unsigned int x, unsigned int y, 
			int *types, int ntypes, bool sortp) {
  int *sorted;
  int index;
  int *typematches = NULL, *matches, nmatches, i, j, k;
  Interval_T interval;

  *ntypematches = 0;
  matches = IIT_get(&nmatches,this,x,y,/*sortp*/false);
  for (i = 0; i < nmatches; i++) {
    index = matches[i];
    interval = &(this->intervals[index-1]);
    k = 0;
    while (k < ntypes && Interval_type(interval) != types[k]) {
      k++;
    }
    if (k < ntypes) {
      (*ntypematches)++;
    }
  }

  if (*ntypematches > 0) {
    typematches = (int *) CALLOC(*ntypematches,sizeof(int));
    j = 0;
    for (i = 0; i < nmatches; i++) {
      index = matches[i];
      interval = &(this->intervals[index-1]);
      k = 0;
      while (k < ntypes && Interval_type(interval) != types[k]) {
	k++;
      }
      if (k < ntypes) {
	typematches[j++] = index;
      }
    }
  }
  
  if (matches != NULL) {
    FREE(matches);
  }

  if (sortp == false) {
    return typematches;
  } else {
    sorted = sort_matches(this,typematches,*ntypematches,/*alphabetizep*/true);
    FREE(typematches);
    return sorted;
  }
}

int
IIT_get_exact (T this, unsigned int x, unsigned int y, int type) {
  int index;
  int *matches, nmatches, i;
  Interval_T interval;

  matches = IIT_get(&nmatches,this,x,y,/*sortp*/false);
  for (i = 0; i < nmatches; i++) {
    index = matches[i];
    interval = &(this->intervals[index-1]);
    if (Interval_low(interval) == x && Interval_high(interval) == y &&
	Interval_type(interval) == type) {
      FREE(matches);
      return index;
    }
  }

  fprintf(stderr,"IIT_get_exact failed for %u %u %d\n",x,y,type);
  exit(9);

  FREE(matches);
  return -1;
}

int *
IIT_get_exact_multiple (int *nexactmatches, T this, unsigned int x, unsigned int y, int type) {
  int *exactmatches;
  int index;
  int *matches, nmatches, i, j;
  Interval_T interval;

  *nexactmatches = 0;
  matches = IIT_get(&nmatches,this,x,y,/*sortp*/false);
  for (i = 0; i < nmatches; i++) {
    index = matches[i];
    interval = &(this->intervals[index-1]);
    if (Interval_low(interval) == x && Interval_high(interval) == y &&
	Interval_type(interval) == type) {
      (*nexactmatches)++;
    }
  }

  if (nexactmatches == 0) {
    fprintf(stderr,"IIT_get_exact_matches failed for %u %u %d\n",x,y,type);
    exit(9);
  } else {
    exactmatches = (int *) CALLOC(*nexactmatches,sizeof(int));
    j = 0;
    for (i = 0; i < nmatches; i++) {
      index = matches[i];
      interval = &(this->intervals[index-1]);
      if (Interval_low(interval) == x && Interval_high(interval) == y &&
	  Interval_type(interval) == type) {
	exactmatches[j++] = index;
      }
    }
    FREE(matches);
    return exactmatches;
  }
}


/* Assume 0-based index */
static void
print_record (T this, int recno, bool map_bothstrands_p, T chromosome_iit, int level) {
  unsigned int start, end, chrpos1, chrpos2;
  char *string, *chrstring;
  int typeint;
  Interval_T interval;
  bool allocp;

  start = this->annotpointers[recno];
  end = this->annotpointers[recno+1];
  
  if (end <= start + 1U) {
    /* No annotation; use label */
    string = IIT_label_to_string(this,recno+1,&allocp);
  } else {
    string = IIT_annotation(this,recno+1,&allocp);
  }

  printf("    %s",this->name);
  interval = &(this->intervals[recno]);
  chrstring = Chrom_string_from_position(&chrpos1,Interval_low(interval),
					 chromosome_iit);
  chrstring = Chrom_string_from_position(&chrpos2,Interval_high(interval),
					 chromosome_iit);
  if (level >= 0) {
    printf("\t%d",level);
  }
  printf("\t%s:%u..%u",chrstring,chrpos1,chrpos2);

  if (map_bothstrands_p == true) {
    if ((typeint = Interval_type(interval)) <= 0) {
      printf("\t\t%s",string);
    } else {
      printf("\t%s\t%s",IIT_typestring(this,typeint),string);
    }
  } else {
    printf("\t%s",string);
  }
  if (string[strlen(string)-1] != '\n') {
    printf("\n");
  }

  if (allocp == true) {
    FREE(string);
  }

  return;
}


void
IIT_print (T this, int *matches, int nmatches, bool map_bothstrands_p,
	   T chromosome_iit, int *levels, bool reversep) {
  int recno, i;

  if (levels == NULL) {
    if (reversep == true) {
      for (i = nmatches-1; i >= 0; i--) {
	recno = matches[i] - 1;	/* Convert to 0-based */
	print_record(this,recno,map_bothstrands_p,chromosome_iit,/*level*/-1);
      }
    } else {
      for (i = 0; i < nmatches; i++) {
	recno = matches[i] - 1;	/* Convert to 0-based */
	print_record(this,recno,map_bothstrands_p,chromosome_iit,/*level*/-1);
      }
    }
  } else {
    if (reversep == true) {
      for (i = nmatches-1; i >= 0; i--) {
	recno = matches[i] - 1;	/* Convert to 0-based */
	print_record(this,recno,map_bothstrands_p,chromosome_iit,levels[i]);
      }
    } else {
      for (i = 0; i < nmatches; i++) {
	recno = matches[i] - 1;	/* Convert to 0-based */
	print_record(this,recno,map_bothstrands_p,chromosome_iit,levels[i]);
      }
    }
  }

  return;
}


/************************************************************************/

/* Retrieves intervals from an IIT where type > 0.  Used by gmapindex to 
   construct altstrain_iit.  Here, the iit is a contig_iit.  */
List_T
IIT_intervallist_typed (List_T *labellist, Uintlist_T *seglength_list, T this) {
  List_T intervallist = NULL;
  Interval_T interval;
  char *annotation, firstchar;
  bool allocp;
  int i;
  unsigned int seglength;

  *labellist = NULL;
  *seglength_list = NULL;
  for (i = 0; i < this->nintervals; i++) {
    interval = &(this->intervals[i]);
    if (Interval_type(interval) > 0) {
      intervallist = List_push(intervallist,Interval_copy(interval));
      *labellist = List_push(*labellist,IIT_label(this,i+1));

      if (this->version <= 1) {
	/* Annotation may be negative to indicate contig is reverse complement */
	annotation = IIT_annotation(this,i+1,&allocp);
	firstchar = annotation[0];
	if (firstchar == '-') {
	  seglength = (unsigned int) strtoul(&(annotation[1]),NULL,10);
	} else {
	  seglength = (unsigned int) strtoul(annotation,NULL,10);
	  *seglength_list = Uintlist_push(*seglength_list,seglength);
	  if (allocp == true) {
	    FREE(annotation);
	  }
	}
      } else {
	seglength = (unsigned int) strtoul(annotation,NULL,10);
	*seglength_list = Uintlist_push(*seglength_list,seglength);
      }
    }
  }
  *labellist = List_reverse(*labellist);
  *seglength_list = Uintlist_reverse(*seglength_list);
  return List_reverse(intervallist);
}


List_T
IIT_typelist (T this) {
  List_T typelist = NULL;
  int i;
  char *typestring, *copy;

  for (i = 0; i < this->ntypes; i++) {
    typestring = IIT_typestring(this,i);
    copy = (char *) CALLOC(strlen(typestring)+1,sizeof(char));
    strcpy(copy,typestring);
    typelist = List_push(typelist,copy);
  }
  return List_reverse(typelist);
}
					  

