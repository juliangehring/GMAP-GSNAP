static char rcsid[] = "$Id: iit-write.c,v 1.19 2005/05/09 22:33:12 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifndef HAVE_MEMCPY
# define memcpy(d,s,n) bcopy((s),(d),(n))
#endif
#ifndef HAVE_MEMMOVE
# define memmove(d,s,n) bcopy((s),(d),(n))
#endif

#include "iit-write.h"
#include "iitdef.h"

#ifdef WORDS_BIGENDIAN
#include "bigendian.h"
#else
#include "littleendian.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>		/* For strlen and strcmp */
#include "assert.h"
#include "bool.h"
#include "mem.h"
#include "interval.h"

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


typedef struct Node_T *Node_T;
struct Node_T {
  int index;
  unsigned int value;
  int a;
  int b;
  Node_T left;
  Node_T right;
};


static Node_T
Node_new () {
  Node_T new = (Node_T) MALLOC(sizeof(*new));

  new->index = -1;
  new->left = new->right = NULL;
  return new;
}

static void 
Node_gc (Node_T *old) {
  if (*old) {
    Node_gc(&((*old)->left));
    Node_gc(&((*old)->right));
    FREE(*old);
  }
  return;
}


static FNode_T
FNode_new (unsigned int value, int a, int b, int leftindex, int rightindex) {
  FNode_T new = (FNode_T) MALLOC(sizeof(*new));

  new->value = value;
  new->a = a;
  new->b = b;
  new->leftindex = leftindex;
  new->rightindex = rightindex;
  return new;
}

static void
FNode_free (FNode_T *old) {
  FREE(*old);
  return;
}


#define T IIT_T

/* Note: This procedure differs from the one in iit-read.c */
void
IIT_free (T *old) {
  if (*old) {
    FREE((*old)->sigmas);
    FREE((*old)->omegas);
    FREE((*old)->intervals);
    FREE(*old);
  }
}

/************************************************************************/

static int 
is_right_split (int i, int j, int r, unsigned int value, int *sigmas, 
		struct Interval_T *intervals) {
  int iota, lambda;

  iota = 0;
  for (lambda = i; lambda <= j; lambda++) {
    if (Interval_is_contained(value,intervals,sigmas[lambda]) == true) {
      iota = lambda;
    }
  }
  return (iota == r);
}  

static bool
is_sorted (int array[], int i, int j, unsigned int (*endpoint)(), struct Interval_T *intervals) {
  int lambda;

  for (lambda = i; lambda <= j - 1; lambda++) {
    if (endpoint(intervals,array[lambda]) > endpoint(intervals,array[lambda + 1])) {
      return false;
    }
  }
  return true;
}

static bool
is_empty (int array[], int i, int j) {
  int lambda;

  for (lambda = i; lambda <= j; lambda++) {
    if (array[lambda] != 0) {
      return false;
    }
  }
  return true;
}

static bool
is_valid_input (int i, int j, int *sigmas, int *omegas, struct Interval_T *intervals) {
  int lambda, iota;

  assert(is_sorted (sigmas, i, j, Interval_array_low, intervals));
  assert(is_empty (omegas, i, j));
  for (lambda = i; lambda <= j; lambda++) {
    iota = sigmas[lambda];
    assert(Interval_is_contained(Interval_array_low(intervals,iota),intervals,iota)
	   || Interval_is_contained(Interval_array_high(intervals,iota),intervals,iota));
	   
  }
  return true;
}

/************************************************************************/

static bool
Node_is_valid_output (Node_T node, int i, int j, int *sigmas, int *omegas, 
		      struct Interval_T *intervals) {
  int lambda;

  assert((i <= node->a) && (node->b <= j));
  assert(is_sorted (sigmas, node->a, node->b, Interval_array_low, intervals));
  assert(is_sorted (omegas, node->a, node->b, Interval_array_high, intervals));
  assert(is_right_split (i, j, node->b, node->value, sigmas, intervals));

  for (lambda = node->a; lambda <= node->b; lambda++) {
    assert(Interval_is_contained(node->value,intervals,sigmas[lambda])
	   && Interval_is_contained(node->value,intervals,omegas[lambda]));
  }

  for (lambda = i; lambda <= node->a - 1; lambda++) {
    assert(Interval_is_contained(node->value,intervals,sigmas[lambda]) == false);
  }
  for (lambda = node->b + 1; lambda <= j; lambda++) {
    assert(Interval_is_contained(node->value,intervals,sigmas[lambda]) == false);
  }

  return true;
}



/* Refinement of node_make().
   Select proper value and split off "right of" triangles. */
static void 
node_select (int *index, unsigned int *value, int i, int j, 
	     int *sigmas, int *omegas, struct Interval_T *intervals) {
  int r = j - (j - i) / 3;
  unsigned int k = Interval_array_low(intervals,sigmas[r]);

  while ((r < j) && (Interval_array_low(intervals,sigmas[r + 1]) == k)) {
    r ++;
  }

  if (Interval_is_contained(k,intervals,sigmas[r]) == false) {
    /* adjust r to the left for "open" intervals */
    while ((r > i) && (Interval_is_contained(k,intervals,sigmas[r-1]) == false)) {
      r --;
      printf(" basic_iit: (-)\n");
    }
    if (Interval_is_contained(k,intervals,sigmas[r]) == false) {
      r --;
      printf(" basic_iit: [-]\n");
      assert(r == i - 1);
      printf(" basic_iit WARNING: empty NODE!?!\n");
    }
  }
  assert(is_right_split (i, j, r, k, sigmas, intervals));
  *index = r;
  *value = k;
  return;
}


/* Makes node out of sigmas[i..j], and recurses. */
static Node_T
Node_make (int *nnodes, int i, int j, int *sigmas, int *omegas, struct Interval_T *intervals) {
  Node_T node;
  int lambda, iota;
  int q, r;

  assert(is_valid_input (i, j, sigmas, omegas, intervals));
  if (i > j) {
    return (Node_T) NULL;
  } else {
    node = Node_new();
    *nnodes += 1;

    /* select value & get "right of" intervals sigma[r+1..j] */
    node_select (&r, &node->value, i, j, sigmas, omegas, intervals);

    /* mark "contains" intervals in sigma[i..r] to omega[i<q+1..r] */
    q = r;
    for (lambda = r; lambda >= i; lambda--) {
      if (Interval_is_contained(node->value,intervals,sigmas[lambda]) == true) {
	omegas[q] = sigmas[lambda];
	sigmas[lambda] = 0;
	q --;
      }
    }

    /* move remaining "left of" intervals from sigma[i..r] to sigma[i..q] */
    iota = i;
    for (lambda = i; lambda <= r; lambda++) {
      if (sigmas[lambda] != 0) {
	sigmas[iota] = sigmas[lambda];
	iota ++;
      }
    }
    assert(iota == q + 1);
    
    /* copy omega[q+1..r] back to sigma[q+1..r] & sort omega[q+1..r] */
    for (lambda = q+1; lambda <= r; lambda++) {
      sigmas[lambda] = omegas[lambda];
    }
    Interval_qsort_by_omega(omegas,q+1,r,intervals);
    node->a = q + 1;
    node->b = r;

#if 0
    fprintf(stderr," NODE=%u [%d..%d], left: %d, cont: %d, right: %d\n",
	    node->value, i, j, q - i + 1, r - q, j - r);
#endif

    assert(Node_is_valid_output (node, i, j, sigmas, omegas, intervals));

    /* recurse */
    node->left  = Node_make(&(*nnodes),i,q,sigmas,omegas,intervals);
    node->right = Node_make(&(*nnodes),r+1,j,sigmas,omegas,intervals);
    
    return node;
  }
}



static void 
Node_index (Node_T node, int *index) {
  if (node != NULL) {
    node->index = (*index)++;
    Node_index(node->left,&(*index));
    Node_index(node->right,&(*index));
  }
  return;
}


static T
IIT_build (Node_T *root, List_T intervallist) {
  T new = (T) MALLOC(sizeof(*new));
  int index = 0;		/* Must be initialized to 0 */
  int i;
  List_T p;

  new->nintervals = List_length(intervallist);
  new->intervals = (struct Interval_T *) CALLOC(new->nintervals,sizeof(struct Interval_T));
  for (p = intervallist, i = 0; i < new->nintervals; i++, p = List_next(p)) {
    memcpy(&(new->intervals[i]),List_head(p),sizeof(struct Interval_T));
  }

  new->nnodes = 0;
  new->sigmas = (int *) CALLOC(new->nintervals+1,sizeof(int));
  for (i = 1; i <= new->nintervals; i++) {
    new->sigmas[i] = i;
  }

  /* Sort sigmas with respect to Interval_array_low */
  Interval_qsort_by_sigma(new->sigmas,1,new->nintervals,new->intervals);

  new->omegas = (int *) CALLOC(new->nintervals+1,sizeof(int));

  /* make first node, and recurse... */
  *root = Node_make(&new->nnodes,1,new->nintervals,new->sigmas,new->omegas,new->intervals);
  Node_index(*root,&index);

  return new;
}


/************************************************************************
 *   Output procedures
 ************************************************************************/

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
 *   labelorder: nintervals*sizeof(int)
 *   labelpointers: (nintervals+1)*sizeof(unsigned int)
 *   labels: nintervals*(variable length strings, including '\0')
 *
 *   annotpointers: (nintervals+1)*sizeof(unsigned int)
 *   annotations: nintervals*(variable length strings, including '\0')
 ************************************************************************/


/* For making labelorder */
struct Sortitem_T {
  int recno;
  char *label;
};

static int
Sortitem_cmp (const void *x, const void *y) {
  struct Sortitem_T *a = (struct Sortitem_T *) x;
  struct Sortitem_T *b = (struct Sortitem_T *) y;
  
  return strcmp(a->label,b->label);
}

static int *
get_labelorder (List_T labellist, int nintervals) {
  int *labelorder, recno, i;
  struct Sortitem_T *sortitems;
  List_T p;
  char *label;

  labelorder = (int *) CALLOC(nintervals,sizeof(int));
  sortitems = (struct Sortitem_T *) CALLOC(nintervals,sizeof(struct Sortitem_T));
  recno = 0;
  for (p = labellist; p != NULL; p = List_next(p)) {
    sortitems[recno].recno = recno;
    sortitems[recno].label = (char *) List_head(p);
    recno++;
  }

  qsort(sortitems,nintervals,sizeof(struct Sortitem_T),Sortitem_cmp);
  for (i = 0; i < nintervals; i++) {
    labelorder[i] = sortitems[i].recno;
  }
  FREE(sortitems);
  return labelorder;
}


/* Prints in DFS order */
static void
Node_fwrite (FILE *fp, Node_T node) {
  FNode_T fnode;
  int leftindex, rightindex;

  if (node != NULL) {
    if (node->left == NULL) {
      leftindex = -1;
    } else {
      leftindex = node->left->index;
    }

    if (node->right == NULL) {
      rightindex = -1;
    } else {
      rightindex = node->right->index;
    }

    fnode = FNode_new(node->value,node->a,node->b,leftindex,rightindex);
    /*
    debug(fprintf(stderr,"Writing node %d %d %d %d %d\n",
		  fnode->value,fnode->a,fnode->b,fnode->leftindex,fnode->rightindex));
    */

    FWRITE_UINT(fnode->value,fp);
    FWRITE_INT(fnode->a,fp);
    FWRITE_INT(fnode->b,fp);
    FWRITE_INT(fnode->leftindex,fp);
    FWRITE_INT(fnode->rightindex,fp);

    FNode_free(&fnode);

    Node_fwrite(fp,node->left);
    Node_fwrite(fp,node->right);
  }
  return;
}

static void
IIT_output (FILE *fp, T this, Node_T root, List_T typelist, List_T labellist, List_T annotlist,
	    Uintlist_T annot_strlen_list) {
  List_T p;
  Uintlist_T q;
  unsigned int pointer = 0U, count;
  char *type, *label, *annot, X[1], endofstring[1];
  int i;
  Interval_T interval;
  int *labelorder;

  this->ntypes = List_length(typelist);
  FWRITE_INT(this->nintervals,fp);
  FWRITE_INT(this->ntypes,fp);
  FWRITE_INT(this->nnodes,fp);
  FWRITE_INTS(this->sigmas,this->nintervals + 1,fp);
  FWRITE_INTS(this->omegas,this->nintervals + 1,fp);
  Node_fwrite(fp,root);

  for (i = 0; i < this->nintervals; i++) {
    FWRITE_UINT(this->intervals[i].low,fp);
    FWRITE_UINT(this->intervals[i].high,fp);
    FWRITE_INT(this->intervals[i].type,fp);
  }

  /* Write type pointers */
  pointer = 0U;
  FWRITE_UINT(pointer,fp);
  for (p = typelist; p != NULL; p = List_next(p)) {
    type = (char *) List_head(p);
    pointer += (unsigned int) strlen(type)+1U;	/* Add count for '\0' */
    FWRITE_UINT(pointer,fp);
  }

  /* Write types */
  for (p = typelist; p != NULL; p = List_next(p)) {
    type = (char *) List_head(p);
    FWRITE_CHARS(type,strlen(type)+1,fp); /* Write '\0' */
  }      


  /* Write labelorder */
  labelorder = get_labelorder(labellist,this->nintervals);
  FWRITE_INTS(labelorder,this->nintervals,fp);
  FREE(labelorder);

  /* Write label pointers */
  pointer = 0;
  FWRITE_INT(pointer,fp);
  for (p = labellist; p != NULL; p = List_next(p)) {
    label = (char *) List_head(p);
    pointer += (unsigned int) strlen(label)+1;	/* Add count for '\0' */
    FWRITE_INT(pointer,fp);
  }

  /* Write labels */
  for (p = labellist; p != NULL; p = List_next(p)) {
    label = (char *) List_head(p);
    FWRITE_CHARS(label,strlen(label)+1,fp); /* Write '\0' */
  }


  if (annotlist == NULL) {
    /* Special hack for altstrain IIT.  Sequence filled in later by IIT_write_annotation */
    pointer = 0U;
    FWRITE_UINT(pointer,fp);
    for (q = annot_strlen_list; q != NULL; q = Uintlist_next(q)) {
      pointer += Uintlist_head(q)+1U; /* Add 1 for '\0' */
      FWRITE_UINT(pointer,fp);
    }

    /* Fill file with X's.  This takes care of case where strain
       sequence is shorter than mapped coordinates */
    X[0] = 'x';
    endofstring[0] = '\0';
    for (q = annot_strlen_list; q != NULL; q = Uintlist_next(q)) {
      for (count = 0; count < Uintlist_head(q); count++) {
	FWRITE_CHARS(X,1,fp);
      }
      FWRITE_CHARS(endofstring,1,fp);
    }

  } else {
    /* Write annot pointers */
    pointer = 0;
    FWRITE_INT(pointer,fp);
    for (p = annotlist; p != NULL; p = List_next(p)) {
      annot = (char *) List_head(p);
      pointer += (unsigned int) strlen(annot)+1;	/* Add count for '\0' */
      FWRITE_INT(pointer,fp);
    }

    /* Write annotations */
    for (p = annotlist; p != NULL; p = List_next(p)) {
      annot = (char *) List_head(p);
      FWRITE_CHARS(annot,strlen(annot)+1,fp); /* Write '\0' */
    }
  }
  return;
}


/* Remember to List_reverse intervallist, typelist, and annotlist if desired */
/* If annotlist is NULL, X's are written */
void
IIT_write (char *iitfile, List_T intervallist, List_T typelist, List_T labellist, List_T annotlist,
	   Uintlist_T annot_strlen_list) {
  Node_T root;
  T iit;
  FILE *fp;

  if ((fp = fopen(iitfile,"wb")) == NULL) {
    fprintf(stderr,"Error: can't open file %s\n",iitfile);
    exit(9);
  } else {
    iit = IIT_build(&root,intervallist);
    IIT_output(fp,iit,root,typelist,labellist,annotlist,annot_strlen_list);
    Node_gc(&root);
    IIT_free(&iit);
    fclose(fp);
    return;
  }
}

void
IIT_backfill_sequence (T this, int index, int offset, char *Buffer) {
  int recno;
  char *ptr;

  recno = index - 1;
  ptr = &(this->annotations[this->annotpointers[recno] + offset]);
  if (strncpy(ptr,Buffer,strlen(Buffer)) == NULL) {
    fprintf(stderr,"Unable to write %s into iit file\n",Buffer);
    exit(9);
  }
  return;
}

