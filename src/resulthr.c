static char rcsid[] = "$Id: resulthr.c 37238 2011-03-27 00:24:03Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "resulthr.h"
#include <stdlib.h>
#include "mem.h"
#include "stage3hr.h"


/* Assignment of resulttype */
#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif



#define T Result_T
struct T {
  Resulttype_T resulttype;
  bool translocationp;
  int id;
  void **array;
  int npaths;
  void **array2;
  int npaths2;
};


Resulttype_T
Result_resulttype (T this) {
  return this->resulttype;
}

bool
Result_translocationp (T this) {
  return this->translocationp;
}


char *
Resulttype_string (Resulttype_T resulttype) {
  switch (resulttype) {
  case SINGLEEND_NOMAPPING: return "singleend_nomapping";
  case PAIREDEND_NOMAPPING: return "pairedend_nomapping";
  case SINGLEEND_UNIQ: return "singleend_uniq";
  case SINGLEEND_MULT: return "singleend_mult";
  case PAIRED_UNIQ: return "paired_uniq";
  case PAIRED_MULT: return "paired_mult";
  case CONCORDANT_UNIQ: return "concordant_uniq";
  case CONCORDANT_MULT: return "concordant_mult";
  case HALFMAPPING_UNIQ: return "halfmapping_uniq";
  case HALFMAPPING_MULT: return "halfmapping_mult";
  case UNPAIRED_UNIQ: return "unpaired_uniq";
  case UNPAIRED_MULT: return "unpaired_mult";
  default: 
    fprintf(stderr,"Unknown resulttype %d\n",resulttype);
    abort();
  }
  return "";
}


int
Result_id (T this) {
  return this->id;
}


void **
Result_array (int *npaths, T this) {
  *npaths = this->npaths;
  return this->array;
}


void **
Result_array2 (int *npaths, T this) {
  *npaths = this->npaths2;
  return this->array2;
}


T
Result_single_read_new (int id, void **resultarray, int npaths) {
  T new = (T) MALLOC(sizeof(*new));

  if (npaths == 0) {
    debug(printf("npaths == 0 => SINGLEEND_NOMAPPING\n"));
    new->translocationp = false;
    new->resulttype = SINGLEEND_NOMAPPING;
  } else {
    if (npaths == 1) {
      debug(printf("npaths == 1 => SINGLEEND_UNIQ\n"));
      new->resulttype = SINGLEEND_UNIQ;
    } else {
      debug(printf("npaths == %d => SINGLEEND_MULT\n",npaths));
      new->resulttype = SINGLEEND_MULT;
    }
    if (Stage3_chrnum((Stage3_T) resultarray[0]) == 0) {
      debug(printf("chrnum of first result == 0 => translocationp is true\n"));
      new->translocationp = true;
    } else {
      debug(printf("chrnum of first result == %d => translocationp is false\n",
		   Stage3_chrnum((Stage3_T) resultarray[0])));
      new->translocationp = false;
    }
  }

  new->id = id;
  new->array = resultarray;
  new->npaths = npaths;

  return new;
}

T
Result_paired_read_new (int id, void **resultarray, int npaths, bool concordantp) {
  T new = (T) MALLOC(sizeof(*new));
  int i;

  if (npaths == 0) {
    abort();

  } else if (concordantp == false) {
    if (npaths == 1) {
      debug(printf("concordantp false and npaths == 1 => PAIRED_UNIQ\n"));
      new->resulttype = PAIRED_UNIQ;
    } else {
      debug(printf("concordantp false and npaths == %d => PAIRED_MULT\n",npaths));
      new->resulttype = PAIRED_MULT;
    }

  } else {
    if (npaths == 1) {
      debug(printf("concordantp true and npaths == 1 => CONCORDANT_UNIQ\n"));
      new->resulttype = CONCORDANT_UNIQ;
    } else {
      debug(printf("concordantp true and npaths == %d => CONCORDANT_MULT\n",npaths));
      new->resulttype = CONCORDANT_MULT;
    }
  }

  if (Stage3_chrnum(Stage3pair_hit5((Stage3pair_T) resultarray[0])) == 0) {
    debug(printf("chrnum of first hit == 0 (npaths %d) => translocationp is true\n",npaths));
    debug(
	  for (i = 0; i < npaths; i++) {
	    printf(" %d-%d",
		   Stage3_chrnum(Stage3pair_hit5((Stage3pair_T) resultarray[i])),
		   Stage3_chrnum(Stage3pair_hit3((Stage3pair_T) resultarray[i])));
	  }
	  printf("\n");
	  );

    new->translocationp = true;
    
  } else if (Stage3_chrnum(Stage3pair_hit3((Stage3pair_T) resultarray[0])) == 0) {
    debug(printf("chrnum of second hit == 0 (npaths %d) => translocationp is true\n",npaths));
    debug(
	  for (i = 0; i < npaths; i++) {
	    printf(" %d-%d",
		   Stage3_chrnum(Stage3pair_hit5((Stage3pair_T) resultarray[i])),
		   Stage3_chrnum(Stage3pair_hit3((Stage3pair_T) resultarray[i])));
	  }
	  printf("\n");
	  );

    new->translocationp = true;

  } else {
    new->translocationp = false;

  }

  new->id = id;
  new->array = resultarray;
  new->npaths = npaths;

  return new;
}

T
Result_paired_as_singles_new (int id, void **hits5, int npaths5, void **hits3, int npaths3) {
  T new = (T) MALLOC(sizeof(*new));

  if (npaths5 == 0 && npaths3 == 0) {
    new->resulttype = PAIREDEND_NOMAPPING;
    new->translocationp = false;
  } else if (npaths5 == 0 && npaths3 == 1) {
    new->resulttype = HALFMAPPING_UNIQ;
  } else if (npaths5 == 1 && npaths3 == 0) {
    new->resulttype = HALFMAPPING_UNIQ;
  } else if (npaths5 == 0 || npaths3 == 0) {
    new->resulttype = HALFMAPPING_MULT;
  } else if (npaths5 == 1 && npaths3 == 1) {
    new->resulttype = UNPAIRED_UNIQ;
  } else {
    new->resulttype = UNPAIRED_MULT;
  }

  if (npaths5 > 0 && Stage3_chrnum((Stage3_T) hits5[0]) == 0) {
    new->translocationp = true;
  } else if (npaths3 > 0 && Stage3_chrnum((Stage3_T) hits3[0]) == 0) {
    new->translocationp = true;
  } else {
    new->translocationp = false;
  }

  new->id = id;
  new->array = hits5;
  new->npaths = npaths5;
  new->array2 = hits3;
  new->npaths2 = npaths3;

  return new;
}

void
Result_free (T *old) {
  int i;
  Stage3_T stage3;
  Stage3pair_T stage3pair;

  switch ((*old)->resulttype) {
  case SINGLEEND_NOMAPPING: case PAIREDEND_NOMAPPING: break;

  case SINGLEEND_UNIQ: case SINGLEEND_MULT:
    for (i = 0; i < (*old)->npaths; i++) {
      stage3 = (*old)->array[i];
      Stage3_free(&stage3);
    }
    FREE((*old)->array);
    break;

  case PAIRED_UNIQ: case PAIRED_MULT: case CONCORDANT_UNIQ: case CONCORDANT_MULT:
    for (i = 0; i < (*old)->npaths; i++) {
      stage3pair = (*old)->array[i];
      Stage3pair_free(&stage3pair);
    }
    FREE((*old)->array);
    break;

  default:
    for (i = 0; i < (*old)->npaths2; i++) {
      stage3 = (*old)->array2[i];
      Stage3_free(&stage3);
    }
    FREE((*old)->array2);

    for (i = 0; i < (*old)->npaths; i++) {
      stage3 = (*old)->array[i];
      Stage3_free(&stage3);
    }
    FREE((*old)->array);
  }

  FREE(*old);
  return;
}


