static char rcsid[] = "$Id: iit-read.c,v 1.58 2005/03/11 17:55:12 twu Exp $";
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
#include <sys/mman.h>		/* For mmap */
#include "assert.h"
#include "mem.h"

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

/* The iit file has a '\0' after each string, so functions know where
   it ends */
char *
IIT_annotation (T this, int index) {
  int recno;
  unsigned int start;

  recno = index - 1; /* Convert to 0-based */
  start = this->annotpointers[recno];
  return &(this->annotations[start]);
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
IIT_dump (T this) {
  int i;
  Interval_T interval;
  int start;
  char *typestring, *annotation;

  for (i = 0; i < this->nintervals; i++) {
    printf(">%s",IIT_label(this,i+1));

    interval = &(this->intervals[i]);
    printf(" %u %u",Interval_low(interval),Interval_high(interval));
    if (Interval_type(interval) > 0) {
      printf(" %s",IIT_typestring(this,Interval_type(interval)));
    }
    printf("\n");

    annotation = IIT_annotation(this,i+1);
    if (strlen(annotation) == 0) {
      /* Don't print anything */
    } else if (annotation[strlen(annotation)-1] == '\n') {
      printf("%s",annotation);
    } else {
      printf("%s\n",annotation);
    }
  }
  return;
}

void
IIT_dump_formatted (T this, bool directionalp) {
  int i;
  Interval_T interval;
  unsigned int start, startpos, endpos;
  char *label, *annotation;

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
      annotation = IIT_annotation(this,i+1);
      if (annotation[0] == '-') {
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
IIT_free_mmapped (T *old) {
  if (*old) {
    if ((*old)->name != NULL) {
      FREE((*old)->name);
    }

    munmap((void *) (*old)->finfo,(*old)->flength);
    close((*old)->fd);

    /* Annotations is mmapped, so don't free. */
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
}

/* This procedure leaks memory, but used only for debugging */
void
IIT_debug (char *filename) {
  T new;
  FILE *fp;
  int offset = 0, i, stringlen, prot, oflag, mapflags;

  if ((fp = fopen(filename,"r")) == NULL) {
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
  for (i = 0; i < new->nnodes; i++) {
    printf(" low:%u high:%u type:%d\n",
	   new->intervals[i].low,new->intervals[i].high,new->intervals[i].type);
  }
  printf("\n");

  new->typepointers = (unsigned int *) CALLOC(new->ntypes+1,sizeof(unsigned int));
  offset += sizeof(int)*FREAD_UINTS(new->typepointers,new->ntypes+1,fp);
  printf("typepointers:");
  for (i = 0; i < new->ntypes+1; i++) {
    printf(" %d",new->typepointers[i]);
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
    printf(" %d",new->labelpointers[i]);
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
  struct stat sb;
  int offset = 0, i, stringlen, prot, oflag, mapflags;

  if ((fp = fopen(filename,"r")) == NULL) {
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

  if (readonlyp == true) {
    oflag = O_RDONLY;
    prot = PROT_READ;
    mapflags = 0
#ifdef HAVE_MMAP_MAP_PRIVATE
			     |MAP_PRIVATE
#endif
#ifdef HAVE_MMAP_MAP_FILE
			     |MAP_FILE
#endif
      ;
  } else {
    oflag = O_RDWR;
    prot = PROT_READ | PROT_WRITE;
    mapflags = 0
#ifdef HAVE_MMAP_MAP_SHARED
			     |MAP_SHARED
#endif
#ifdef HAVE_MMAP_MAP_FILE
			     |MAP_FILE
#endif
      ;
  }

  if ((new->fd = open(filename,oflag,0764)) < 0) {
    fprintf(stderr,"Error: can't open file %s\n",filename);
    return NULL;
  }
  fstat(new->fd,&sb);
  new->flength = sb.st_size;

  new->finfo = (char *) mmap(NULL,sb.st_size,prot,mapflags,
			     new->fd,0);
  if (new->finfo == MAP_FAILED) {
    fprintf(stderr,"Error: mmap failed for file %s\n",
	    filename);
    return NULL;
  }
  madvise((caddr_t) new->finfo,new->flength,MADV_RANDOM);
  new->annotations = (char *) &(new->finfo[offset]);

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
      if (Interval_is_contained(this->intervals,x,this->sigmas[lambda]) == true) {
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
      if (Interval_is_contained(this->intervals,x,this->omegas[lambda]) == true) {
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

Intlist_T
IIT_find (T this, char *label) {
  Intlist_T matches = NULL;
  int low, middle, high, recno;
  bool foundp = false;
  int cmp;

  low = 0;
  high = this->nintervals;

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

    for (recno = low; recno <= high; recno++) {
      debug(printf("Pushing %d:%d\n",recno,this->labelorder[recno]));
      matches = Intlist_push(matches,this->labelorder[recno]+1);
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
    while (isspace(*p)) {
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
  Intlist_T matches;

  matches = IIT_find(this,label);
  if (Intlist_length(matches) == 0) {
    /*
    fprintf(stderr,"Expected one match for %s, but got 0\n",
	    label);
    */
    index = -1;
  } else {
    if (Intlist_length(matches) > 1) {
      fprintf(stderr,"Expected one match for %s, but got %d\n",
	      label,Intlist_length(matches));
    }
    index = Intlist_head(matches);
  }
  Intlist_free(&matches);
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

Intlist_T
IIT_get (T this, unsigned int x, unsigned int y) {
  Intlist_T matches = NULL, p;
  int *array, nmatches, i;
  int lambda, prev;
  int min1 = this->nintervals+1, max1 = 0, min2 = this->nintervals+1, max2 = 0;

  debug(printf("%u %u\n",x,y));
  fnode_query_aux(&min1,&max1,this,0,x);
  fnode_query_aux(&min2,&max2,this,0,y);
  debug(fprintf(stderr,"%d %d  %d %d\n",min1,max1,min2,max2));
  for (lambda = min1; lambda <= max2; lambda++) {
    if (Interval_overlap_p(this->intervals,x,y,this->sigmas[lambda]) == true) {
      matches = Intlist_push(matches,this->sigmas[lambda]);
      debug(printf("Pushing %d\n",this->sigmas[lambda]));
    }
  }
  for (lambda = min1; lambda <= max2; lambda++) {
    if (Interval_overlap_p(this->intervals,x,y,this->omegas[lambda]) == true) {
      matches = Intlist_push(matches,this->omegas[lambda]);
      debug(printf("Pushing %d\n",this->omegas[lambda]));
    }
  }

  /* Eliminate duplicates */
  if (matches != NULL) {
    array = Intlist_to_array(&nmatches,matches);
    Intlist_free(&matches);
    matches = (Intlist_T) NULL;

    qsort(array,nmatches,sizeof(int),int_compare);
      
    prev = array[0];
    for (i = 1; i < nmatches; i++) {
      if (array[i] == prev) {
	array[i] = 0;		/* recnos are 1-based, so 0 means dup */
      } else {
	prev = array[i];
      }
    }
    
    for (i = 0; i < nmatches; i++) {
      if (array[i] != 0) {
	matches = Intlist_push(matches,array[i]);
      }
    }
    FREE(array);
  }

  return Intlist_reverse(matches);
}

/* Generally called where intervals don't overlap, like chromosomes,
   and where x == y. */
int
IIT_get_one (T this, unsigned int x, unsigned int y) {
  int index;
  Intlist_T matches;

  matches = IIT_get(this,x,y);
  if (Intlist_length(matches) != 1) {
    fprintf(stderr,"Expected one match for %u--%u, but got %d\n",
	    x,y,Intlist_length(matches));
    abort();
  }
  index = Intlist_head(matches);
  Intlist_free(&matches);
  return index;
}

Intlist_T
IIT_get_typed (T this, unsigned int x, unsigned int y, int type) {
  int index;
  Intlist_T typematches = NULL, matches, p;
  Interval_T interval;

  matches = IIT_get(this,x,y);
  for (p = matches; p != NULL; p = Intlist_next(p)) {
    index = Intlist_head(p);
    interval = &(this->intervals[index-1]);
    if (Interval_type(interval) == type) {
      typematches = Intlist_push(typematches,index);
    }
  }
  Intlist_free(&matches);
  return Intlist_reverse(typematches);
}

int
IIT_get_exact (T this, unsigned int x, unsigned int y, int type) {
  int index;
  Intlist_T matches, p;
  Interval_T interval;

  matches = IIT_get(this,x,y);
  for (p = matches; p != NULL; p = Intlist_next(p)) {
    index = Intlist_head(p);
    interval = &(this->intervals[index-1]);
    if (Interval_low(interval) == x && Interval_high(interval) == y &&
	Interval_type(interval) == type) {
      Intlist_free(&matches);
      return index;
    }
  }

  fprintf(stderr,"IIT_get_exact failed for %u %u %d\n",x,y,type);
  exit(9);
  Intlist_free(&matches);
  return -1;
}


/* Assume 0-based index */
static void
print_record (T this, int recno, bool map_bothstrands_p) {
  unsigned int start, end;
  char *string;
  int typeint;
  Interval_T interval;

  start = this->annotpointers[recno];
  end = this->annotpointers[recno+1];
  
  string = (char *) CALLOC(end - start + 1,sizeof(char));
  strncpy(string,&(this->annotations[start]),end-start);

  if (map_bothstrands_p == true) {
    interval = &(this->intervals[recno]);
    if ((typeint = Interval_type(interval)) <= 0) {
      printf("    %s\t\t%s",this->name,string);
    } else {
      printf("    %s\t%s\t%s",this->name,IIT_typestring(this,typeint),string);
    }
  } else {
    printf("    %s\t%s",this->name,string);
  }
  if (string[strlen(string)-1] != '\n') {
    printf("\n");
  }

  FREE(string);

  return;
}


void
IIT_print (T this, Intlist_T matches, bool map_bothstrands_p) {
  Intlist_T p;
  int recno;

  for (p = matches; p != NULL; p = Intlist_next(p)) {
    recno = Intlist_head(p) - 1; /* Convert to 0-based */
    print_record(this,recno,map_bothstrands_p);
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
  char *annotation;
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
      annotation = IIT_annotation(this,i+1);
      if (annotation[0] == '-') {
	annotation += 1;
      }
      seglength = (unsigned int) strtoul(annotation,NULL,10);
      *seglength_list = Uintlist_push(*seglength_list,seglength);
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
					  

