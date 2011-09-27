static char rcsid[] = "$Id: iit-read.c,v 1.80 2006/11/13 04:05:29 twu Exp $";
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
#include "fopen.h"
#include "chrom.h"
#include "uintlist.h"

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
IIT_nintervals (T this) {
  return this->nintervals;
}


int
IIT_ntypes (T this) {
  return this->ntypes;
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
IIT_dump (T this, bool annotationonlyp) {
  int i;
  Interval_T interval;
  char *annotation;
  bool allocp;

  for (i = 0; i < this->nintervals; i++) {
    if (annotationonlyp == false) {
      printf(">%s",IIT_label(this,i+1));

      interval = &(this->intervals[i]);
      printf(" %u %u",Interval_low(interval),Interval_high(interval));
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
    } else {
      firstchar = IIT_annotation_firstchar(this,i+1);
      if (firstchar == '-') {
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


void
IIT_dump_counts (T this, bool alphabetizep) { 
  int type, index, i, t, k;
  Interval_T interval;
  Uintlist_T *edgelists;
  int *matches, nmatches, nedges;
  unsigned int *edges, edge, lastedge;
  char *typestring;
  Chrom_T *chroms;

  edgelists = (Uintlist_T *) CALLOC(this->ntypes,sizeof(Uintlist_T));
  for (i = 0; i < this->nintervals; i++) {
    interval = &(this->intervals[i]);
    type = Interval_type(interval);
    edgelists[type] = Uintlist_push(edgelists[type],Interval_low(interval));
    edgelists[type] = Uintlist_push(edgelists[type],Interval_high(interval));
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

    if (Uintlist_length(edgelists[type]) > 0) {
      edges = Uintlist_to_array(&nedges,edgelists[type]);
      qsort(edges,nedges,sizeof(unsigned int),uint_cmp);

      edge = edges[0];
      matches = IIT_get_typed(&nmatches,this,edge,edge,type);
      printf("%s\t%u\t%d",typestring,edge,nmatches);

      index = matches[0];
      printf("\t%s",IIT_label(this,index));
      for (k = 1; k < nmatches; k++) {
	index = matches[k];
	printf(",%s",IIT_label(this,index));
      }
      printf("\n");

      lastedge = edge;

      for (i = 1; i < nedges; i++) {
	if ((edge = edges[i]) != lastedge) {
	  matches = IIT_get_typed(&nmatches,this,edge,edge,type);
	  printf("%s\t%u\t%d",typestring,edge,nmatches);
	  
	  index = matches[0];
	  printf("\t%s",IIT_label(this,index));
	  for (k = 1; k < nmatches; k++) {
	    index = matches[k];
	    printf(",%s",IIT_label(this,index));
	  }
	  printf("\n");
	  
	  lastedge = edge;
	}
      }

      Uintlist_free(&(edgelists[type]));
      FREE(edges);
    }

  }

  if (alphabetizep == true) {
    for (t = 0; t < this->ntypes; t++) {
      Chrom_free(&(chroms[t]));
    }
    FREE(chroms);
  }

  FREE(edgelists);

  return;
}

/************************************************************************
 * File format:
 *   nintervals: sizeof(int)
 *   ntypes: sizeof(int)
 *   nnodes: sizeof(int) 
 *   sigmas: (nintervals+1)*sizeof(int)
 *   omegas: (nintervals+1)*sizeof(int)
 *   nodes: nnodes*sizeof(struct FNode_T)
 *   intervals: nintervals*sizeof(struct Interval_T)
 *
 *   typepointers: (ntypes+1)*sizeof(unsigned int)
 *   types: ntypes*(variable length strings, including '\0')
 *
 *   labelorder: nintervals*sizeof(int);
 *   labelpointers: (nintervals+1)*sizeof(unsigned int)
 *   labels: nintervals*(variable length strings, including '\0')
 *
 *   annotpointers: (nintervals+1)*sizeof(unsigned int)
 *   annotations: nintervals*(variable length strings, including '\0')
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
    } else {
      abort();
    }

    FREE((*old)->annotpointers);
    FREE((*old)->labels);
    FREE((*old)->labelpointers);
    FREE((*old)->labelorder);
    FREE((*old)->typestrings);
    FREE((*old)->typepointers);
    FREE((*old)->intervals);
    FREE((*old)->nodes);
    FREE((*old)->omegas);
    FREE((*old)->sigmas);

    FREE(*old);

  }
  return;
}

/* This procedure leaks memory, but used only for debugging */
void
IIT_debug (char *filename) {
  T new;
  FILE *fp;
  off_t offset = 0;
  int i, stringlen;

  if ((fp = FOPEN_READ_BINARY(filename)) == NULL) {
    fprintf(stderr,"Can't open file %s\n",filename);
    return;
  } else {
    new = (T) MALLOC(sizeof(*new));
  }

  offset += sizeof(int)*FREAD_INT(&new->nintervals,fp);
  printf("nintervals: %d\n",new->nintervals);
  offset += sizeof(int)*FREAD_INT(&new->ntypes,fp);
  printf("ntypes: %d\n",new->ntypes);
  offset += sizeof(int)*FREAD_INT(&new->nnodes,fp);
  printf("nnodes: %d\n",new->nnodes);

  new->sigmas = (int *) CALLOC(new->nintervals+1,sizeof(int));
  offset += sizeof(int)*FREAD_INTS(new->sigmas,new->nintervals+1,fp);
  printf("sigmas:");
  for (i = 0; i < new->nintervals+1; i++) {
    printf(" %d",new->sigmas[i]);
  }
  printf("\n");

  new->omegas = (int *) CALLOC(new->nintervals+1,sizeof(int));
  offset += sizeof(int)*FREAD_INTS(new->omegas,new->nintervals+1,fp);
  printf("omegas:");
  for (i = 0; i < new->nintervals+1; i++) {
    printf(" %d",new->omegas[i]);
  }
  printf("\n");

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
    Bigendian_fread_int(&(new->intervals[i].type),fp);
  }
  offset += (sizeof(unsigned int)+sizeof(unsigned int)+sizeof(int))*new->nintervals;
#else
  if (sizeof(struct Interval_T) == sizeof(unsigned int)+sizeof(unsigned int)+sizeof(int)) {
    offset += sizeof(struct Interval_T)*fread(new->intervals,sizeof(struct Interval_T),new->nintervals,fp);
  } else {
    for (i = 0; i < new->nintervals; i++) {
      fread(&(new->intervals[i].low),sizeof(unsigned int),1,fp);
      fread(&(new->intervals[i].high),sizeof(unsigned int),1,fp);
      fread(&(new->intervals[i].type),sizeof(int),1,fp);
    }
    offset += (sizeof(unsigned int)+sizeof(unsigned int)+sizeof(int))*new->nintervals;
  }
#endif
  printf("intervals:\n");
  for (i = 0; i < new->nintervals; i++) {
    printf(" low:%u high:%u type:%d\n",
	   new->intervals[i].low,new->intervals[i].high,new->intervals[i].type);
  }
  printf("\n");

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
  printf("stringlen of labelpointers = %d\n",stringlen);
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

  fclose(fp);
  return;
}


T
IIT_read (char *filename, char *name, bool readonlyp) {
  T new;
  FILE *fp;
  off_t offset = 0;
  int i, stringlen, prot, oflag, mapflags;

  if ((fp = FOPEN_READ_BINARY(filename)) == NULL) {
    return NULL;
  } else {
    new = (T) MALLOC(sizeof(*new));
  }

  if (name != NULL) {
    new->name = (char *) CALLOC(strlen(name)+1,sizeof(char));
    strcpy(new->name,name);
  } else {
    new->name = NULL;
  }

  offset += sizeof(int)*FREAD_INT(&new->nintervals,fp);
  offset += sizeof(int)*FREAD_INT(&new->ntypes,fp);
  offset += sizeof(int)*FREAD_INT(&new->nnodes,fp);

  new->sigmas = (int *) CALLOC(new->nintervals+1,sizeof(int));
  offset += sizeof(int)*FREAD_INTS(new->sigmas,new->nintervals+1,fp);

  new->omegas = (int *) CALLOC(new->nintervals+1,sizeof(int));
  offset += sizeof(int)*FREAD_INTS(new->omegas,new->nintervals+1,fp);

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

  new->intervals = (struct Interval_T *) CALLOC(new->nintervals,sizeof(struct Interval_T));
#ifdef WORDS_BIGENDIAN
  for (i = 0; i < new->nintervals; i++) {
    Bigendian_fread_uint(&(new->intervals[i].low),fp);
    Bigendian_fread_uint(&(new->intervals[i].high),fp);
    Bigendian_fread_int(&(new->intervals[i].type),fp);
  }
  offset += (sizeof(unsigned int)+sizeof(unsigned int)+sizeof(int))*new->nintervals;
#else
  if (sizeof(struct Interval_T) == sizeof(unsigned int)+sizeof(unsigned int)+sizeof(int)) {
    offset += sizeof(struct Interval_T)*fread(new->intervals,sizeof(struct Interval_T),new->nintervals,fp);
  } else {
    for (i = 0; i < new->nintervals; i++) {
      fread(&(new->intervals[i].low),sizeof(unsigned int),1,fp);
      fread(&(new->intervals[i].high),sizeof(unsigned int),1,fp);
      fread(&(new->intervals[i].type),sizeof(int),1,fp);
    }
    offset += (sizeof(unsigned int)+sizeof(unsigned int)+sizeof(int))*new->nintervals;
  }
#endif

  new->typepointers = (unsigned int *) CALLOC(new->ntypes+1,sizeof(unsigned int));
  offset += sizeof(int)*FREAD_UINTS(new->typepointers,new->ntypes+1,fp);

  stringlen = new->typepointers[new->ntypes];
  new->typestrings = (char *) CALLOC(stringlen,sizeof(char));
  offset += sizeof(char)*FREAD_CHARS(new->typestrings,stringlen,fp);

  new->labelorder = (int *) CALLOC(new->nintervals,sizeof(int));
  offset += sizeof(int)*FREAD_INTS(new->labelorder,new->nintervals,fp);

  new->labelpointers = (unsigned int *) CALLOC(new->nintervals+1,sizeof(unsigned int));
  offset += sizeof(int)*FREAD_UINTS(new->labelpointers,new->nintervals+1,fp);
  
  stringlen = new->labelpointers[new->nintervals];
  new->labels = (char *) CALLOC(stringlen,sizeof(char));
  offset += sizeof(char)*FREAD_CHARS(new->labels,stringlen,fp);

  new->annotpointers = (unsigned int *) CALLOC(new->nintervals+1,sizeof(unsigned int));
  offset += sizeof(int)*FREAD_UINTS(new->annotpointers,new->nintervals+1,fp);

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

int *
IIT_get (int *nmatches, T this, unsigned int x, unsigned int y) {
  int *matches = NULL, *uniq, neval, nuniq, i, j;
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
	debug(printf("Pushing overlapping segment %d\n",uniq[i]));
      } else {
	debug(printf("Not pushing non-overlapping segment %d\n",uniq[i]));
      }
    }

    FREE(uniq);
  }

  return matches;
}

void
IIT_get_flanking (int **leftflanks, int *nleftflanks, int **rightflanks, int *nrightflanks,
		  T this, unsigned int x, unsigned int y, int nflanking) {
  int lambda;
  int min1 = this->nintervals+1, max1 = 0, min2 = this->nintervals+1, max2 = 0;

  debug(printf("Entering IIT_get_flanking with query %u %u, nflanking = %d\n",x,y,nflanking));
  fnode_query_aux(&min1,&max1,this,0,x);
  fnode_query_aux(&min2,&max2,this,0,y);
  debug(printf("min1=%d max1=%d  min2=%d max2=%d\n",min1,max1,min2,max2));

  /* Look at sigmas for right flank */
  *rightflanks = (int *) CALLOC(nflanking,sizeof(int));
  lambda = min1;
  while (lambda <= this->nintervals && Interval_low(&(this->intervals[this->sigmas[lambda]-1])) <= y) {
    lambda++;
  }

  *nrightflanks = 0;
  while (*nrightflanks < nflanking && lambda <= this->nintervals) {
    (*rightflanks)[(*nrightflanks)++] = this->sigmas[lambda];
    debug(printf("Pushed lambda %d, segment %d as right flank number %d\n",lambda,this->sigmas[lambda],*nrightflanks));
    lambda++;
  }

  /* Look at omegas for left flank */
  *leftflanks = (int *) CALLOC(nflanking,sizeof(int));
  lambda = max2;
  while (lambda >= 1 && Interval_high(&(this->intervals[this->omegas[lambda]-1])) >= x) {
    lambda--;
  }
  
  *nleftflanks = 0;
  while (*nleftflanks < nflanking && lambda >= 1) {
    (*leftflanks)[(*nleftflanks)++] = this->omegas[lambda];
    debug(printf("Pushed lambda %d, segment %d as left flank number %d\n",lambda,this->omegas[lambda],*nleftflanks));
    lambda--;
  }

  return;
}

void
IIT_get_flanking_typed (int **leftflanks, int *nleftflanks, int **rightflanks, int *nrightflanks,
			T this, unsigned int x, unsigned int y, int nflanking, int type) {
  int lambda;
  int min1 = this->nintervals+1, max1 = 0, min2 = this->nintervals+1, max2 = 0;
  Interval_T interval;
  bool stopp;

  debug(printf("Entering IIT_get with query %u %u\n",x,y));
  fnode_query_aux(&min1,&max1,this,0,x);
  fnode_query_aux(&min2,&max2,this,0,y);
  debug(printf("min1=%d max1=%d  min2=%d max2=%d\n",min1,max1,min2,max2));

  /* Look at sigmas for right flank */
  lambda = min1;
  *rightflanks = (int *) CALLOC(nflanking,sizeof(int));
  *nrightflanks = 0;
  stopp = false;
  while (lambda <= this->nintervals && stopp == false) {
    interval = &(this->intervals[this->sigmas[lambda]-1]);
    if (Interval_low(interval) <= y) {
      lambda++;
    } else if (Interval_type(interval) != type) {
      lambda++;
    } else {
      (*rightflanks)[(*nrightflanks)++] = this->sigmas[lambda];
      if (*nrightflanks < nflanking) {
	lambda++;
      } else {
	stopp = true;
      }
    }
  }

  /* Look at omegas for left flank */
  lambda = max2;
  *leftflanks = (int *) CALLOC(nflanking,sizeof(int));
  *nleftflanks = 0;
  stopp = false;
  while (lambda >= 1 && stopp == false) {
    interval = &(this->intervals[this->omegas[lambda]-1]);
    if (Interval_high(interval) >= x) {
      lambda--;
    } else if (Interval_type(interval) != type) {
      lambda--;
    } else {
      (*leftflanks)[(*nleftflanks)++] = this->omegas[lambda];
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
  int min1 = this->nintervals+1, max1 = 0, min2 = this->nintervals+1, max2 = 0;
  Interval_T interval;
  bool stopp;

  debug(printf("Entering IIT_get with query %u %u\n",x,y));
  fnode_query_aux(&min1,&max1,this,0,x);
  fnode_query_aux(&min2,&max2,this,0,y);
  debug(printf("min1=%d max1=%d  min2=%d max2=%d\n",min1,max1,min2,max2));

  /* Look at sigmas for right flank */
  lambda = min1;
  *rightflanks = (int *) CALLOC(nflanking,sizeof(int));
  *nrightflanks = 0;
  stopp = false;
  while (lambda <= this->nintervals && stopp == false) {
    interval = &(this->intervals[this->sigmas[lambda]-1]);
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
	(*rightflanks)[(*nrightflanks)++] = this->sigmas[lambda];
	if (*nrightflanks < nflanking) {
	  lambda++;
	} else {
	  stopp = true;
	}
      }
    }
  }

  /* Look at omegas for left flank */
  lambda = max2;
  *leftflanks = (int *) CALLOC(nflanking,sizeof(int));
  *nleftflanks = 0;
  stopp = false;
  while (lambda >= 1 && stopp == false) {
    interval = &(this->intervals[this->omegas[lambda]-1]);
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
	stopp = true;
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

  matches = IIT_get(&nmatches,this,x,y);
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
IIT_get_typed (int *ntypematches, T this, unsigned int x, unsigned int y, int type) {
  int index;
  int *typematches = NULL, *matches, nmatches, i, j;
  Interval_T interval;

  *ntypematches = 0;
  matches = IIT_get(&nmatches,this,x,y);
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

  return typematches;
}

int *
IIT_get_multiple_typed (int *ntypematches, T this, unsigned int x, unsigned int y, 
			int *types, int ntypes) {
  int index;
  int *typematches = NULL, *matches, nmatches, i, j, k;
  Interval_T interval;

  *ntypematches = 0;
  matches = IIT_get(&nmatches,this,x,y);
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

  return typematches;
}

int
IIT_get_exact (T this, unsigned int x, unsigned int y, int type) {
  int index;
  int *matches, nmatches, i;
  Interval_T interval;

  matches = IIT_get(&nmatches,this,x,y);
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

      /* Annotation may be negative to indicate contig is reverse complement */
      annotation = IIT_annotation(this,i+1,&allocp);
      firstchar = annotation[0];
      if (firstchar == '-') {
	seglength = (unsigned int) strtoul(&(annotation[1]),NULL,10);
      } else {
	seglength = (unsigned int) strtoul(annotation,NULL,10);
      }
      *seglength_list = Uintlist_push(*seglength_list,seglength);
      if (allocp == true) {
	FREE(annotation);
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
					  

