static char rcsid[] = "$Id: mutation.c,v 1.11 2005/07/12 16:35:03 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "mutation.h"
#include <stdio.h>
#include <stdlib.h>
#include "mem.h"


#define T Mutation_T
struct T {
  Muttype_T muttype;

  int refquerypos;
  union {
    int aapos;
    int deletionstart;
  } a;
  union {
    char aa_g;
    char aa_g1;
  } c1;
  union {
    char aa_e;
    char aa_g2;
  } c2;
  union {
    int ninsertions;
    int deletionend;
  } z;
};

Muttype_T
Mutation_type (T this) {
  return this->muttype;
}

int
Mutation_refquerypos (T this) {
  return this->refquerypos;
}

int
Mutation_aapos (T this) {
  return this->a.aapos;
}

int
Mutation_naa (T this) {
  switch (this->muttype) {
  case SUBSTITUTION: return 1;
  case INSERTION: return this->z.ninsertions;
  case DELETION: return this->z.deletionend - this->a.deletionstart + 1;
  case INVALID: return 0;
  }
  return 0;
}

bool
Mutation_at_ends_p (T segmental, T single) {
  /*
    printf("Checking ends of ");
    Mutation_print(segmental);
    printf(" and ");
    Mutation_print(single);
    printf("\n");
  */

  if (segmental->muttype == INSERTION) {
    if (single->a.aapos == segmental->a.aapos) {
      return true;
    } else if (single->a.aapos == segmental->a.aapos - 1) {
      return true;
    } else if (single->a.aapos == segmental->a.aapos + segmental->z.ninsertions - 1) {
      return true;
    } else if (single->a.aapos == segmental->a.aapos + segmental->z.ninsertions - 1 + 1) {
      return true;
    } else {
      return false;
    }
  } else if (segmental->muttype == DELETION) {
    if (single->a.aapos == segmental->a.deletionstart) {
      return true;
    } else if (single->a.aapos == segmental->a.deletionstart - 1) {
      return true;
    } else if (single->a.aapos == segmental->z.deletionend) {
      return true;
    } else if (single->a.aapos == segmental->z.deletionend + 1) {
      return true;
    } else {
      return false;
    }
  } else {
    return false;
  }
}


void
Mutation_print (T this) {

  switch (this->muttype) {
  case SUBSTITUTION:
    printf("%c%d%c",this->c1.aa_g,this->a.aapos,this->c2.aa_e);
    break;
  case INSERTION:
    if (this->z.ninsertions == 1) {
      printf("ins%d%c",this->a.aapos,this->c2.aa_e);
    } else {
      printf("ins%d+%daa",this->a.aapos,this->z.ninsertions);    }
    break;
  case DELETION:
    if (this->a.deletionstart == this->z.deletionend) {
      printf("del%c%d",this->c1.aa_g1,this->a.deletionstart);
    } else {
      printf("del%c%d-%c%d",this->c1.aa_g1,this->a.deletionstart,this->c2.aa_g2,this->z.deletionend);
    }
    break;
  case INVALID:
    break;
  }

  return;
}

int
Mutation_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  int aapos1, aapos2;

  if (x->muttype == DELETION) {
    aapos1 = x->a.deletionstart;
  } else {
    aapos1 = x->a.aapos;
  }

  if (y->muttype == DELETION) {
    aapos2 = y->a.deletionstart;
  } else {
    aapos2 = y->a.aapos;
  }

  if (aapos1 < aapos2) {
    return -1;
  } else if (aapos1 > aapos2) {
    return 1;
  } else {
    return 0;
  }
}

T
Mutation_substitution_new (int refquerypos, int aapos, char aa_g, char aa_e) {
  T new = (T) MALLOC(sizeof(*new));
  
  new->muttype = SUBSTITUTION;
  new->refquerypos = refquerypos;
  new->a.aapos = aapos;
  new->c1.aa_g = aa_g;
  new->c2.aa_e = aa_e;
  
  return new;
}

T
Mutation_insertion_new (int refquerypos, int aapos, char aa_e, int ninsertions) {
  T new = (T) MALLOC(sizeof(*new));
  
  new->muttype = INSERTION;
  new->refquerypos = refquerypos;
  new->a.aapos = aapos;
  new->c2.aa_e = aa_e;
  new->z.ninsertions = ninsertions;
  
  return new;
}

T
Mutation_deletion_new (int refquerypos, int deletionstart, char aa_g1, int deletionend, char aa_g2) {
  T new = (T) MALLOC(sizeof(*new));
  
  new->muttype = DELETION;
  new->refquerypos = refquerypos;
  new->a.deletionstart = deletionstart;
  new->z.deletionend = deletionend;
  new->c1.aa_g1 = aa_g1;
  new->c2.aa_g2 = aa_g2;
  
  return new;
}


void
Mutation_free (T *old) {
  FREE(*old);
  return;
}

T
Mutation_merge_insertions (T mutation1, T mutation2) {
  /* Have to check muttype of mutation1 and mutation2 to exclude INVALIDs */
  if (mutation1->muttype == INSERTION && mutation2->muttype == INSERTION &&
      mutation1->a.aapos == mutation2->a.aapos) {
    return Mutation_insertion_new(mutation1->refquerypos,mutation1->a.aapos,'*',
				  mutation1->z.ninsertions + mutation2->z.ninsertions);
  } else {
    return NULL;
  }
}

T
Mutation_merge_deletions (T mutation1, T mutation2) {
/* Have to check muttype of mutation1 and mutation2 to exclude INVALIDs */
  if (mutation1->muttype == DELETION && mutation2->muttype == DELETION) {
    if (mutation1->z.deletionend + 1 == mutation2->a.deletionstart) {
      return Mutation_deletion_new(mutation1->refquerypos,mutation1->a.deletionstart,mutation1->c1.aa_g1,
				   mutation2->z.deletionend,mutation2->c2.aa_g2);
      
    } else if (mutation2->z.deletionend + 1 == mutation1->a.deletionstart) {
      return Mutation_deletion_new(mutation2->refquerypos,mutation2->a.deletionstart,mutation2->c1.aa_g1,
				   mutation1->z.deletionend,mutation1->c2.aa_g2);
      
    } else {
      return NULL;
    }
  } else {
    return NULL;
  }
}


void
Mutation_invalidate (T this) {
  this->muttype = INVALID;
  return;
}

bool
Mutation_validp (T this) {
  if (this->muttype == INVALID) {
    return false;
  } else {
    return true;
  }
}

