static char rcsid[] = "$Id: mutation.c,v 1.8 2005/02/07 23:56:56 twu Exp $";
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
  };
  union {
    char aa_g;
    char aa_g1;
  };
  union {
    char aa_e;
    char aa_g2;
  };
  union {
    int ninsertions;
    int deletionend;
  };
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
  return this->aapos;
}

int
Mutation_naa (T this) {
  switch (this->muttype) {
  case SUBSTITUTION: return 1;
  case INSERTION: return this->ninsertions;
  case DELETION: return this->deletionend - this->deletionstart + 1;
  case INVALID: return 0;
  }
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
    if (single->aapos == segmental->aapos) {
      return true;
    } else if (single->aapos == segmental->aapos - 1) {
      return true;
    } else if (single->aapos == segmental->aapos + segmental->ninsertions - 1) {
    } else if (single->aapos == segmental->aapos + segmental->ninsertions - 1 + 1) {
      return true;
    } else {
      return false;
    }
  } else if (segmental->muttype == DELETION) {
    if (single->aapos == segmental->deletionstart) {
      return true;
    } else if (single->aapos == segmental->deletionstart - 1) {
      return true;
    } else if (single->aapos == segmental->deletionend) {
      return true;
    } else if (single->aapos == segmental->deletionend + 1) {
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
  char c;

  switch (this->muttype) {
  case SUBSTITUTION:
    printf("%c%d%c",this->aa_g,this->aapos,this->aa_e);
    break;
  case INSERTION:
    if (this->ninsertions == 1) {
      printf("ins%d%c",this->aapos,this->aa_e);
    } else {
      printf("ins%d+%daa",this->aapos,this->ninsertions);    }
    break;
  case DELETION:
    if (this->deletionstart == this->deletionend) {
      printf("del%c%d",this->aa_g1,this->deletionstart);
    } else {
      printf("del%c%d-%c%d",this->aa_g1,this->deletionstart,this->aa_g2,this->deletionend);
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

  if (x->aapos < y->aapos) {
    return -1;
  } else if (x->aapos > y->aapos) {
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
  new->aapos = aapos;
  new->aa_g = aa_g;
  new->aa_e = aa_e;
  
  return new;
}

T
Mutation_insertion_new (int refquerypos, int aapos, char aa_e, int ninsertions) {
  T new = (T) MALLOC(sizeof(*new));
  
  new->muttype = INSERTION;
  new->refquerypos = refquerypos;
  new->aapos = aapos;
  new->aa_e = aa_e;
  new->ninsertions = ninsertions;
  
  return new;
}

T
Mutation_deletion_new (int refquerypos, int deletionstart, char aa_g1, int deletionend, char aa_g2) {
  T new = (T) MALLOC(sizeof(*new));
  
  new->muttype = DELETION;
  new->refquerypos = refquerypos;
  new->deletionstart = deletionstart;
  new->deletionend = deletionend;
  new->aa_g1 = aa_g1;
  new->aa_g2 = aa_g2;
  
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
      mutation1->aapos == mutation2->aapos) {
    return Mutation_insertion_new(mutation1->refquerypos,mutation1->aapos,'*',
				  mutation1->ninsertions + mutation2->ninsertions);
  } else {
    return NULL;
  }
}

T
Mutation_merge_deletions (T mutation1, T mutation2) {
/* Have to check muttype of mutation1 and mutation2 to exclude INVALIDs */
  if (mutation1->muttype == DELETION && mutation2->muttype == DELETION) {
    if (mutation1->deletionend + 1 == mutation2->deletionstart) {
      return Mutation_deletion_new(mutation1->refquerypos,mutation1->deletionstart,mutation1->aa_g1,
				   mutation2->deletionend,mutation2->aa_g2);
      
    } else if (mutation2->deletionend + 1 == mutation1->deletionstart) {
      return Mutation_deletion_new(mutation2->refquerypos,mutation2->deletionstart,mutation2->aa_g1,
				   mutation1->deletionend,mutation1->aa_g2);
      
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

