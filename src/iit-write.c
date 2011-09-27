static char rcsid[] = "$Id: iit-write.c,v 1.37 2009/08/29 00:36:12 twu Exp $";
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
#include "fopen.h"
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


#if 0
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
#endif


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
    node_select(&r,&node->value,i,j,sigmas,omegas,intervals);

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


static int
IIT_count_nnodes (List_T intervallist) {
  int nnodes;
  Node_T root;
  int nintervals, i;
  List_T p;
  struct Interval_T *intervals;
  int *sigmas;
  int *omegas;

  if ((nintervals = List_length(intervallist)) == 0) {
    return 0;
  } else {
    intervals = (struct Interval_T *) CALLOC(nintervals,sizeof(struct Interval_T));
    for (p = intervallist, i = 0; i < nintervals; i++, p = List_next(p)) {
      memcpy(&(intervals[i]),List_head(p),sizeof(struct Interval_T));
    }

    /* IIT ordering of intervals */
    nnodes = 0;
    sigmas = (int *) CALLOC(nintervals+1,sizeof(int));
    for (i = 1; i <= nintervals; i++) {
      sigmas[i] = i;
    }

    /* Sort sigmas with respect to Interval_array_low */
    Interval_qsort_by_sigma(sigmas,1,nintervals,intervals);
    
    omegas = (int *) CALLOC(nintervals+1,sizeof(int));

    /* make first node, and recurse... */
    root = Node_make(&nnodes,1,nintervals,sigmas,omegas,intervals);

    Node_gc(&root);
    FREE(omegas);
    FREE(sigmas);
    FREE(intervals);

    return nnodes;
  }
}

static void
IIT_build_one_div (Node_T *root, struct Interval_T **intervals, int **alphas, int **betas, int **sigmas, int **omegas,
		   int *nnodes, List_T intervallist, int nintervals) {
  int index = 0;		/* Must be initialized to 0 */
  int i;
  List_T p;

  if (nintervals > 0) {
    *intervals = (struct Interval_T *) CALLOC(nintervals,sizeof(struct Interval_T));
    for (p = intervallist, i = 0; i < nintervals; i++, p = List_next(p)) {
      memcpy(&((*intervals)[i]),List_head(p),sizeof(struct Interval_T));
    }
  }

  /* Strict ordering of intervals */
  
  *alphas = (int *) CALLOC(nintervals+1,sizeof(int));
  *betas = (int *) CALLOC(nintervals+1,sizeof(int));
  for (i = 1; i <= nintervals; i++) {
    (*alphas)[i] = (*betas)[i] = i;
  }
  Interval_qsort_by_sigma(*alphas,1,nintervals,*intervals);
  Interval_qsort_by_omega(*betas,1,nintervals,*intervals);


  /* IIT ordering of intervals */
  *nnodes = 0;
  *sigmas = (int *) CALLOC(nintervals+1,sizeof(int));
  for (i = 1; i <= nintervals; i++) {
    (*sigmas)[i] = i;
  }

  /* Sort sigmas with respect to Interval_array_low */
  Interval_qsort_by_sigma(*sigmas,1,nintervals,*intervals);

  *omegas = (int *) CALLOC(nintervals+1,sizeof(int));

  /* make first node, and recurse... */
  *root = Node_make(&(*nnodes),1,nintervals,*sigmas,*omegas,*intervals);
  Node_index(*root,&index);

  return;
}

/************************************************************************
 *   Output procedures
 ************************************************************************/

/************************************************************************
 * File format:
 *   0: sizeof(int)  [in version >= 2]
 *   version: sizeof(int)  [in version >= 2]
 *
 *   total_nintervals: sizeof(int)  [so maximum is 2 billion]
 *   ntypes: sizeof(int)
 *   nfields: sizeof(int)  [in version >= 2]
 *
 *   ndivs: sizeof(int) [in version >= 3.  In version < 3, set ndivs = 1.]
 *   nintervals: ndivs*sizeof(int) [in version >= 3]
 *   cum_nintervals: (ndivs+1)*sizeof(int) [in version >= 3]
 *   nnodes: ndivs*sizeof(int)
 *   cum_nnodes: (ndivs+1)*sizeof(int) [in version >= 3]
 *
 *   divsort: sizeof(int) [in version >= 3: 0 = none, 1 = alpha, 2 = chrom]
 *   divpointers: (ndivs+1)*sizeof(unsigned int)  [in version >= 3]
 *   divs: ndivs*(variable length strings, including '\0')  [in version >= 3]
 *
 *   for each div (recall that in version < 3, ndivs = 1):
 *     alphas: (nintervals[div]+1)*sizeof(int)  [in version >= 2]
 *     betas: (nintervals[div]+1)*sizeof(int)  [in version >= 2]
 *     sigmas: (nintervals[div]+1)*sizeof(int)
 *     omegas: (nintervals[div]+1)*sizeof(int)
 *     nodes: nnodes[div]*sizeof(struct FNode_T)
 *
 *   intervals: total_nintervals[div]*sizeof(struct Interval_T)  [3 ints in version 1; 4 ints in version >= 2]
 *
 *   typepointers: (ntypes+1)*sizeof(unsigned int)
 *   types: ntypes*(variable length strings, including '\0')
 *
 *   fieldpointers: (nfields+1)*sizeof(unsigned int)  [in version >= 2]
 *   fields: nfields*(variable length strings, including '\0')  [in version >= 2]
 *
 *   labelorder: total_nintervals*sizeof(int);
 *   labelpointers: (total_nintervals+1)*sizeof(unsigned int)  [changes to long unsigned int in v4]
 *   labels: total_nintervals*(variable length strings, including '\0')
 *
 *   annotpointers: (total_nintervals+1)*sizeof(unsigned int)  [changes to long unsigned int in v4]
 *   annotations: total_nintervals*(variable length strings, including '\0')
 *
 *   Note: labels and annotations are organized by division, so we can
 *   find the label/annot by looking up index = cum_nintervals[div]+recno.
 *   This is because alphas/betas/sigmas/omegas run from 1 to nintervals[div]
 ************************************************************************/


/* For making labelorder */
struct Sortitem_T {
  int divno;
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
get_labelorder (List_T divlist, Table_T labeltable, int *cum_nintervals, int total_nintervals) {
  int *labelorder, divno, recno, i, k = 0;
  struct Sortitem_T *sortitems;
  char *divstring;
  List_T labellist, d, p;

  labelorder = (int *) CALLOC(total_nintervals,sizeof(int));
  sortitems = (struct Sortitem_T *) CALLOC(total_nintervals,sizeof(struct Sortitem_T));
  divno = 0;

  for (d = divlist; d != NULL; d = List_next(d)) {
    divstring = (char *) List_head(d);
    labellist = (List_T) Table_get(labeltable,(void *) divstring);
    recno = 0;
    for (p = labellist; p != NULL; p = List_next(p)) {
      sortitems[k].divno = divno;
      sortitems[k].recno = recno;
      sortitems[k].label = (char *) List_head(p);
      k++;
      recno++;
    }
    divno++;
  }

  qsort(sortitems,total_nintervals,sizeof(struct Sortitem_T),Sortitem_cmp);
  for (i = 0; i < total_nintervals; i++) {
    labelorder[i] = cum_nintervals[sortitems[i].divno] + sortitems[i].recno;
  }
  FREE(sortitems);
  return labelorder;
}


/* Prints in DFS order */
static void
Node_fwrite (FILE *fp, Node_T node) {
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

    debug(printf("Writing node %u %d %d\n",node->value,node->a,node->b));
    FWRITE_UINT(node->value,fp);
    FWRITE_INT(node->a,fp);
    FWRITE_INT(node->b,fp);
    FWRITE_INT(leftindex,fp);
    FWRITE_INT(rightindex,fp);

    Node_fwrite(fp,node->left);
    Node_fwrite(fp,node->right);
  }
  return;
}

/* Stores in DFS order */
static void
Node_store (int *fnodei, struct FNode_T *fnodes, Node_T node) {
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

    fnodes[*fnodei].value = node->value;
    fnodes[*fnodei].a = node->a;
    fnodes[*fnodei].b = node->b;
    fnodes[*fnodei].leftindex = leftindex;
    fnodes[*fnodei].rightindex = rightindex;
    (*fnodei)++;

    Node_store(&(*fnodei),fnodes,node->left);
    Node_store(&(*fnodei),fnodes,node->right);
  }
  return;
}


static void
IIT_output_div_header (FILE *fp, int total_nintervals, int ntypes, int nfields, int ndivs, 
		       int *nintervals, int *nnodes, int *cum_nintervals, int *cum_nnodes, 
		       List_T divlist, Sorttype_T divsort, int version) {
  int new_format_indicator = 0;
  int divno;
  unsigned int pointer = 0U;
  List_T p;
  char *divstring;

  if (version >= 2) {
    FWRITE_INT(new_format_indicator,fp); /* Indicates new format, since nintervals > 0 */
    FWRITE_INT(version,fp);
  }

  FWRITE_INT(total_nintervals,fp);
  FWRITE_INT(ntypes,fp);
  if (version >= 2) {
    FWRITE_INT(nfields,fp);
  }

  if (version >= 3) {
    FWRITE_INT(ndivs,fp);
    for (divno = 0; divno < ndivs; divno++) {
      FWRITE_INT(nintervals[divno],fp);
    }
    for (divno = 0; divno <= ndivs; divno++) {
      FWRITE_INT(cum_nintervals[divno],fp);
    }
  }
  for (divno = 0; divno < ndivs; divno++) {
    FWRITE_INT(nnodes[divno],fp);
  }
  if (version >= 3) {
    for (divno = 0; divno <= ndivs; divno++) {
      FWRITE_INT(cum_nnodes[divno],fp);
    }
  }

  if (version >= 3) {
    FWRITE_INT(divsort,fp);
  }

  if (version >= 3) {
    /* Write div pointers */
    pointer = 0U;
    FWRITE_UINT(pointer,fp);
    for (p = divlist; p != NULL; p = List_next(p)) {
      divstring = (char *) List_head(p);
      pointer += (unsigned int) strlen(divstring)+1U;	/* Add count for '\0' */
      FWRITE_UINT(pointer,fp);
    }

    /* Write divs */
    for (p = divlist; p != NULL; p = List_next(p)) {
      divstring = (char *) List_head(p);
      FWRITE_CHARS(divstring,strlen(divstring)+1,fp); /* Write '\0' */
    }
  }

  return;
}


static void
IIT_output_one_div (FILE *fp, Node_T root, struct Interval_T *intervals, int *alphas, int *betas, int *sigmas, int *omegas,
		    int nintervals, int version) {
#ifdef DEBUG
  int i;
#endif

  if (version >= 2) {
    FWRITE_INTS(alphas,nintervals + 1,fp);
    FWRITE_INTS(betas,nintervals + 1,fp);

    debug(
	  printf("alphas:");
	  for (i = 0; i < nintervals + 1; i++) {
	    printf(" %d",alphas[i]);
	  }
	  printf("\n");

	  printf("betas:");
	  for (i = 0; i < nintervals + 1; i++) {
	    printf(" %d",betas[i]);
	  }
	  printf("\n");
	  );

  }
  FWRITE_INTS(sigmas,nintervals + 1,fp);
  FWRITE_INTS(omegas,nintervals + 1,fp);

  debug(
	printf("sigmas:");
	for (i = 0; i < nintervals + 1; i++) {
	  printf(" %d",sigmas[i]);
	}
	printf("\n");

	printf("omegas:");
	for (i = 0; i < nintervals + 1; i++) {
	  printf(" %d",omegas[i]);
	}
	printf("\n");
	);

  debug(printf("Starting to write nodes\n"));
  Node_fwrite(fp,root);
  debug(printf("\n"));

  return;
}

static void
IIT_output_div_footer (FILE *fp, int ndivs, List_T divlist, List_T typelist, List_T fieldlist,
		       Table_T intervaltable, Table_T labeltable, Table_T annottable, 
		       int *cum_nintervals, int total_nintervals, int version) {
  List_T intervallist, labellist, annotlist, d, p;
  Interval_T interval;
#ifdef HAVE_64_BIT
  UINT8 pointer8 = 0LU;
#endif
  unsigned int pointer = 0U;
  int *labelorder;
  char *divstring, *typestring, *fieldstring, *label, *annot;

#ifndef HAVE_64_BIT
  if (version >= 4) {
    fprintf(stderr,"Machine is not 64-bit, so cannot write version 4 files.  Try writing version 3 instead.\n");
    exit(9);
  }
#endif

  for (d = divlist; d != NULL; d = List_next(d)) {
    divstring = (char *) List_head(d);
    intervallist = (List_T) Table_get(intervaltable,(void *) divstring);
    for (p = intervallist; p != NULL; p = List_next(p)) {
      interval = (Interval_T) List_head(p);
      FWRITE_UINT(interval->low,fp);
      FWRITE_UINT(interval->high,fp);
      if (version >= 2) {
	FWRITE_INT(interval->sign,fp);
      }
      FWRITE_INT(interval->type,fp);
    }
  }

  /* Write type pointers */
  pointer = 0U;
  FWRITE_UINT(pointer,fp);
  for (p = typelist; p != NULL; p = List_next(p)) {
    typestring = (char *) List_head(p);
    if (typestring == NULL) {
      abort();			/* Even type 0 has an empty string */
    } else {
      pointer += (unsigned int) strlen(typestring)+1U;	/* Add count for '\0' */
      FWRITE_UINT(pointer,fp);
    }
  }

  /* Write types */
  for (p = typelist; p != NULL; p = List_next(p)) {
    typestring = (char *) List_head(p);
    FWRITE_CHARS(typestring,strlen(typestring)+1,fp); /* Write '\0' */
  }      

  if (version >= 2) {
    /* Write field pointers */
    pointer = 0U;
    FWRITE_UINT(pointer,fp);
    for (p = fieldlist; p != NULL; p = List_next(p)) {
      fieldstring = (char *) List_head(p);
      pointer += (unsigned int) strlen(fieldstring)+1U;	/* Add count for '\0' */
      FWRITE_UINT(pointer,fp);
    }

    /* Write fields */
    for (p = fieldlist; p != NULL; p = List_next(p)) {
      fieldstring = (char *) List_head(p);
      FWRITE_CHARS(fieldstring,strlen(fieldstring)+1,fp); /* Write '\0' */
    }
  }

  /* Write labelorder */
  labelorder = get_labelorder(divlist,labeltable,cum_nintervals,total_nintervals);
  FWRITE_INTS(labelorder,total_nintervals,fp);
  FREE(labelorder);
  FREE(cum_nintervals);

  /* Write label pointers */
#ifdef HAVE_64_BIT
  if (version >= 4) {
    pointer8 = 0LU;
    FWRITE_UINT8(pointer8,fp);
  } else {
    pointer = 0U;
    FWRITE_UINT(pointer,fp);
  }
#else
  pointer = 0U;
  FWRITE_UINT(pointer,fp);
#endif
  for (d = divlist; d != NULL; d = List_next(d)) {
    divstring = (char *) List_head(d);
    labellist = (List_T) Table_get(labeltable,(void *) divstring);
    for (p = labellist; p != NULL; p = List_next(p)) {
      label = (char *) List_head(p);
#ifdef HAVE_64_BIT      
      if (version >= 4) {
	pointer8 += (UINT8) strlen(label)+1LU;	/* Add count for '\0' */
	FWRITE_UINT8(pointer8,fp);
      } else {
	pointer += (unsigned int) strlen(label)+1U;	/* Add count for '\0' */
	FWRITE_UINT(pointer,fp);
      }
#else
      pointer += (unsigned int) strlen(label)+1U;	/* Add count for '\0' */
      FWRITE_UINT(pointer,fp);
#endif
    }
  }

  /* Write labels */
  for (d = divlist; d != NULL; d = List_next(d)) {
    divstring = (char *) List_head(d);
    labellist = (List_T) Table_get(labeltable,(void *) divstring);
    for (p = labellist; p != NULL; p = List_next(p)) {
      label = (char *) List_head(p);
      FWRITE_CHARS(label,strlen(label)+1,fp); /* Write '\0' */
    }
  }

  /* Write annot pointers */
#ifdef HAVE_64_BIT
  if (version >= 4) {
    pointer8 = 0LU;
    FWRITE_UINT8(pointer8,fp);
  } else {
    pointer = 0U;
    FWRITE_UINT(pointer,fp);
  }
#else
  pointer = 0U;
  FWRITE_UINT(pointer,fp);
#endif
  for (d = divlist; d != NULL; d = List_next(d)) {
    divstring = (char *) List_head(d);
    annotlist = (List_T) Table_get(annottable,(void *) divstring);
    for (p = annotlist; p != NULL; p = List_next(p)) {
      annot = (char *) List_head(p);
#ifdef HAVE_64_BIT
      if (version >= 4) {
	pointer8 += (UINT8) strlen(annot)+1LU;	/* Add count for '\0' */
	FWRITE_UINT8(pointer8,fp);
      } else {
	pointer += (unsigned int) strlen(annot)+1;	/* Add count for '\0' */
	FWRITE_UINT(pointer,fp);
      }
#else
      pointer += (unsigned int) strlen(annot)+1;	/* Add count for '\0' */
      FWRITE_UINT(pointer,fp);
#endif
    }
  }

  /* Write annotations */
  for (d = divlist; d != NULL; d = List_next(d)) {
    divstring = (char *) List_head(d);
    annotlist = (List_T) Table_get(annottable,(void *) divstring);
    for (p = annotlist; p != NULL; p = List_next(p)) {
      annot = (char *) List_head(p);
      FWRITE_CHARS(annot,strlen(annot)+1,fp); /* Write '\0' */
    }
  }

  return;
}

#if 0
static void
compute_flanking (T this) {
  int i;

  this->alphas = (int *) CALLOC(this->nintervals+1,sizeof(int));
  this->betas = (int *) CALLOC(this->nintervals+1,sizeof(int));
  for (i = 1; i <= this->nintervals; i++) {
    this->alphas[i] = this->betas[i] = i;
  }
  Interval_qsort_by_sigma(this->alphas,1,this->nintervals,this->intervals);
  Interval_qsort_by_omega(this->betas,1,this->nintervals,this->intervals);
  return;
}

void
IIT_output_direct (char *iitfile, T this, int version) {
  FILE *fp;
  off_t stringlen;
  int new_format_indicator = 0, i;
  FNode_T node;

  if ((fp = FOPEN_WRITE_BINARY(iitfile)) == NULL) {
    fprintf(stderr,"Error: can't open file %s\n",iitfile);
    exit(9);
  }

  if (version >= 2) {
    FWRITE_INT(new_format_indicator,fp); /* Indicates new format, since nintervals > 0 */
    FWRITE_INT(version,fp);
  }
  FWRITE_INT(this->nintervals,fp);
  FWRITE_INT(this->ntypes,fp);
  if (version < 2) {
    this->nfields = 0;
  } else {
    FWRITE_INT(this->nfields,fp);
  }
  FWRITE_INT(this->nnodes,fp);
  if (version >= 2) {
    if (this->alphas == NULL) {
      compute_flanking(this);
    }
    FWRITE_INTS(this->alphas,this->nintervals + 1,fp);
    FWRITE_INTS(this->betas,this->nintervals + 1,fp);
  }
  FWRITE_INTS(this->sigmas,this->nintervals + 1,fp);
  FWRITE_INTS(this->omegas,this->nintervals + 1,fp);

  /* Write nodes directly */
  for (i = 0; i < this->nnodes; i++) {
    node = &(this->nodes[i]);
    FWRITE_UINT(node->value,fp);
    FWRITE_INT(node->a,fp);
    FWRITE_INT(node->b,fp);
    FWRITE_INT(node->leftindex,fp);
    FWRITE_INT(node->rightindex,fp);
  }

  for (i = 0; i < this->nintervals; i++) {
    FWRITE_UINT(this->intervals[i].low,fp);
    FWRITE_UINT(this->intervals[i].high,fp);
    if (version >= 2) {
      FWRITE_INT(this->intervals[i].sign,fp);
    }
    FWRITE_INT(this->intervals[i].type,fp);
  }

  /* Write types directly */
  for (i = 0; i < this->ntypes+1; i++) {
    FWRITE_UINT(this->typepointers[i],fp);
  }
  if ((stringlen = this->typepointers[this->ntypes]) == 0) {
    fprintf(stderr,"Error in writing types: type stringlen is 0.\n");
    exit(9);
  } else {
    FWRITE_CHARS(this->typestrings,stringlen,fp);
  }

  if (version >= 2) {
    /* Write fields directly */
    for (i = 0; i < this->nfields+1; i++) {
      FWRITE_UINT(this->fieldpointers[i],fp);
    }
    stringlen = this->fieldpointers[this->nfields];
    if (stringlen > 0) {
      FWRITE_CHARS(this->fieldstrings,stringlen,fp);
    }
  }

  /* Write labelorder */
  FWRITE_INTS(this->labelorder,this->nintervals,fp);

  /* Write labels directly */
#ifdef HAVE_64_BIT
  if (version >= 4) {
    FWRITE_UINT8S(this->labelpointers,this->nintervals+1,fp);
  } else {
    FWRITE_UINTS(this->labelpointers,this->nintervals+1,fp);
  }
#else
  FWRITE_UINTS(this->labelpointers,this->nintervals+1,fp);
#endif
  if ((stringlen = this->labelpointers[this->nintervals]) == 0) {
    fprintf(stderr,"Error in writing labels: label stringlen is 0.\n");
    exit(9);
  } else {
    FWRITE_CHARS(this->labels,stringlen,fp);
  }

  /* Write annotations directly */
#ifdef HAVE_64_BIT
  if (version >= 4) {
    FWRITE_UINT8S(this->annotpointers,this->nintervals+1,fp);
  } else {
    FWRITE_UINTS(this->annotpointers,this->nintervals+1,fp);
  }
#else
  FWRITE_UINTS(this->annotpointers,this->nintervals+1,fp);
#endif
  if ((stringlen = this->annotpointers[this->nintervals]) == 0) {
    fprintf(stderr,"Error in writing annotations: annotation stringlen is 0.\n");
    exit(9);
  } else {
    FWRITE_CHARS(this->annotations,stringlen,fp);
  }

  fclose(fp);

  return;
}
#endif


/* If annotlist is NULL, X's are written */
void
IIT_write (char *iitfile, List_T divlist, List_T typelist, List_T fieldlist, Table_T intervaltable,
	   Table_T labeltable, Table_T annottable, Sorttype_T divsort, int version) {
  Node_T root;
  FILE *fp;
  List_T intervallist, d;
  char *divstring;
  int ndivs, total_nintervals, *nintervals, *nnodes, nnodes_one_div, divno;
  int *cum_nintervals, *cum_nnodes;
  struct Interval_T *intervals;
  int *alphas, *betas, *sigmas, *omegas;

  if ((fp = FOPEN_WRITE_BINARY(iitfile)) == NULL) {
    fprintf(stderr,"Error: can't open file %s\n",iitfile);
    exit(9);
  } else {
    ndivs = List_length(divlist);
    nintervals = (int *) CALLOC(ndivs,sizeof(int));
    nnodes = (int *) CALLOC(ndivs,sizeof(int));
    for (d = divlist, divno = 0; d != NULL; d = List_next(d), divno++) {
      divstring = (char *) List_head(d);
      intervallist = (List_T) Table_get(intervaltable,(void *) divstring);
      nintervals[divno] = List_length(intervallist);
      nnodes[divno] = IIT_count_nnodes(intervallist);
    }

    cum_nintervals = (int *) CALLOC(ndivs+1,sizeof(int));
    cum_nintervals[0] = 0;
    for (divno = 1; divno <= ndivs; divno++) {
      cum_nintervals[divno] = cum_nintervals[divno-1] + nintervals[divno-1];
    }
    total_nintervals = cum_nintervals[ndivs];
    if (total_nintervals == 0) {
      fprintf(stderr,"Error: No intervals were seen in input file\n");
      exit(9);
    }

    cum_nnodes = (int *) CALLOC(ndivs+1,sizeof(int));
    cum_nnodes[0] = 0;
    for (divno = 1; divno <= ndivs; divno++) {
      cum_nnodes[divno] = cum_nnodes[divno-1] + nnodes[divno-1];
    }

    fprintf(stderr,"Writing IIT file header information...");
    IIT_output_div_header(fp,total_nintervals,List_length(typelist),List_length(fieldlist),
			  ndivs,nintervals,nnodes,cum_nintervals,cum_nnodes,divlist,divsort,version);
    fprintf(stderr,"done\n");

    for (d = divlist, divno = 0; d != NULL; d = List_next(d), divno++) {
      divstring = (char *) List_head(d);
      intervallist = (List_T) Table_get(intervaltable,(void *) divstring);

      if (divstring[0] == '\0') {
	fprintf(stderr,"Processing null division/chromosome...sorting...");
      } else {
	fprintf(stderr,"Processing division/chromosome %s...sorting...",divstring);
      }
      IIT_build_one_div(&root,&intervals,&alphas,&betas,&sigmas,&omegas,&nnodes_one_div,intervallist,nintervals[divno]);

      fprintf(stderr,"writing...");
      IIT_output_one_div(fp,root,intervals,alphas,betas,sigmas,omegas,nintervals[divno],version);
      fprintf(stderr,"done\n");

      Node_gc(&root);
      FREE(omegas);
      FREE(sigmas);
      FREE(betas);
      FREE(alphas);
      FREE(intervals);
    }

    fprintf(stderr,"Writing IIT file footer information...");
    IIT_output_div_footer(fp,ndivs,divlist,typelist,fieldlist,intervaltable,
			  labeltable,annottable,cum_nintervals,total_nintervals,version);
    fprintf(stderr,"done\n");

    fclose(fp);
    return;
  }
}


#if 0
/* Problematic because this file no longer defines T */
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
#endif

