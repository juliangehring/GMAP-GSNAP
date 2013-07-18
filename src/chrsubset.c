static char rcsid[] = "$Id: chrsubset.c 99737 2013-06-27 19:33:03Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "chrsubset.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>		/* For strlen */
#include <ctype.h>		/* For isspace */
#include "mem.h"
#include "fopen.h"


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


#define T Chrsubset_T
struct T {
  char *name;
  int nincluded;
  bool *includep;
  int *newindices;
  int *oldindices;
};


T
Chrsubset_make (Chrnum_T chrnum, Univ_IIT_T chromosome_iit) {
  T new = (T) MALLOC(sizeof(*new));
  int nchromosomes, i;

  new->name = NULL;
  new->nincluded = 1;

  nchromosomes = Univ_IIT_total_nintervals(chromosome_iit);
  new->includep = (bool *) CALLOC(nchromosomes,sizeof(bool));
  new->includep[chrnum] = true;

  new->newindices = (int *) CALLOC(nchromosomes,sizeof(int));
  for (i = 0; i < nchromosomes; i++) {
    new->newindices[i] = -1;
  }
  new->newindices[chrnum] = 1;

  new->oldindices = (int *) CALLOC(1,sizeof(int));
  new->oldindices[0] = chrnum + 1;
  
  return new;
}


void
Chrsubset_print (T this) {
  if (this != NULL && this->name != NULL) {
    printf("  [chrsubset: %s]",this->name);
  }
  return;
}

void
Chrsubset_print_chromosomes (T this, Univ_IIT_T chromosome_iit) {
  int i;
  bool firstp = true, allocp;
  char *label;

  for (i = 0; i < Univ_IIT_total_nintervals(chromosome_iit); i++) {
    if (this == NULL || this->includep[i] == true) {
      label = Univ_IIT_label(chromosome_iit,i+1,&allocp);
      if (firstp == true) {
	printf("%s",label);
	firstp = false;
      } else {
	printf(",%s",label);
      }
      if (allocp == true) {
	FREE(label);
      }
    }
  }

  return;
}

char *
Chrsubset_name (T this) {
  return this->name;
}

int
Chrsubset_nincluded (T this, Univ_IIT_T chromosome_iit) {
  if (this == NULL) {
    return Univ_IIT_total_nintervals(chromosome_iit);
  } else {
    return this->nincluded;
  }
}


/* Assumes that this != NULL.  Should call this routine with "if (!chrsubset or Chrsubset_includep()..." */
bool
Chrsubset_includep (T this, Univcoord_T position, Univ_IIT_T chromosome_iit) {
  int index;

  index = Univ_IIT_get_one(chromosome_iit,position,position);
  return this->includep[index-1];
}


/* index here is 1-based */
int
Chrsubset_newindex (T this, int index) {
  if (this == NULL) {
    return index;
  } else {
    return this->newindices[index-1];
  }
}

/* index here is 0-based */
int
Chrsubset_oldindex (T this, int index) {
  if (this == NULL) {
    return index;
  } else {
    return this->oldindices[index-1];
  }
}

#if 0
unsigned int *
Chrsubset_transitions (int **signs, int *nedges, T this, Univ_IIT_T chromosome_iit) {
  if (this == NULL) {
    return Univ_IIT_transitions(&(*signs),&(*nedges),chromosome_iit);
  } else {
    return Univ_IIT_transitions_subset(&(*signs),&(*nedges),chromosome_iit,this->oldindices,this->nincluded);
  }
}
#endif


void
Chrsubset_free (T *old) {
  if (*old != NULL) {
    if ((*old)->name != NULL) {
      FREE((*old)->name);
    }
    if ((*old)->includep != NULL) {
      FREE((*old)->includep);
    }
    if ((*old)->newindices != NULL) {
      FREE((*old)->newindices);
    }
    if ((*old)->oldindices != NULL) {
      FREE((*old)->oldindices);
    }
    FREE(*old);
  }
  return;
}

#define BUFSIZE 2048

static char *
read_header (FILE *fp, char *filename) {
  char *subsetname, Buffer[BUFSIZE], *p, *start, *end;

  Buffer[0] = '\0';

  while (Buffer[0] == '\0' || isspace((int) Buffer[0])) {
    /* Read past all empty lines */
    if (fgets(Buffer,BUFSIZE,fp) == NULL) {
      return NULL;
    }
  }
  if (Buffer[0] != '>') {
    fprintf(stderr,"In chromosome subset file %s, a '>' line was expected\n",filename);
    exit(9);
  } else {
    p = &(Buffer[1]);
    if (*p == '\0' || isspace((int) *p)) {
      fprintf(stderr,"The '>' line in chromosome subset file %s is invalid\n",filename);
      exit(9);
    } else {
      start = p;
      while (*p != '\0' && !isspace((int) *p)) {
	p++;
      }
      end = p;

      subsetname = (char *) CALLOC((size_t) (end - start +1),sizeof(char));
      strncpy(subsetname,start,(size_t) (end - start));
      debug(printf("Read header %s\n",subsetname));
      return subsetname;
    }
  }
}

static void
skip_list (FILE *fp, char *filename, char *subsetname) {
  char Buffer[BUFSIZE];

  if (fgets(Buffer,BUFSIZE,fp) == NULL) {
    fprintf(stderr,"In %s, expected a line after >%s\n",filename,subsetname);
    exit(9);
  }
  return;
}


/*
static bool *
include_all (Univ_IIT_T chromosome_iit) {
  bool *includep;
  int nchromosomes, i;

  nchromosomes = Univ_IIT_total_nintervals(chromosome_iit);
  includep = (bool *) CALLOC(nchromosomes,sizeof(bool));
  for (i = 0; i < nchromosomes; i++) {
    includep[i] = true;
  }
  return includep;
}
*/

static char *
get_next_token (char **string) {
  char *token, *start, *end, *p;

  p = *string;
  while (*p != '\0' && (isspace((int) *p) || *p == ',')) {
    p++;
  }
  if (*p == '\0') {
    *string = p;
    return NULL;
  } else {
    start = p;
  }

  while (*p != '\0' && !isspace((int) *p) && *p != ',') {
    p++;
  }
  end = p;

  token = (char *) CALLOC((size_t) (end - start + 1),sizeof(char));
  strncpy(token,start,(size_t) (end - start));
  *string = p;
  return token;
}



static bool *
process_inclusions (char *string, Univ_IIT_T chromosome_iit) {
  bool *includep;
  int nchromosomes, index, i;
  char *chrstring;

  /* Default is to exclude */
  nchromosomes = Univ_IIT_total_nintervals(chromosome_iit);
  includep = (bool *) CALLOC(nchromosomes,sizeof(bool));
  for (i = 0; i < nchromosomes; i++) {
    includep[i] = false;
  }

  while ((chrstring = get_next_token(&string)) != NULL) {
    debug(printf("Token is %s\n",chrstring));
    index = Univ_IIT_find_one(chromosome_iit,chrstring);
    includep[index-1] = true;
    FREE(chrstring);
  }
  return includep;
}

static bool *
process_exclusions (char *string, Univ_IIT_T chromosome_iit) {
  bool *includep;
  int nchromosomes, index, i;
  char *chrstring;

  /* Default is to include */
  nchromosomes = Univ_IIT_total_nintervals(chromosome_iit);
  includep = (bool *) CALLOC(nchromosomes,sizeof(bool));
  for (i = 0; i < nchromosomes; i++) {
    includep[i] = true;
  }

  while ((chrstring = get_next_token(&string)) != NULL) {
    debug(printf("Token is %s\n",chrstring));
    index = Univ_IIT_find_one(chromosome_iit,chrstring);
    includep[index-1] = false;
    FREE(chrstring);
  }
  return includep;
}


static bool *
process_list (FILE *fp, char *filename, char *subsetname, Univ_IIT_T chromosome_iit) {
  char Buffer[BUFSIZE], *p;

  if (fgets(Buffer,BUFSIZE,fp) == NULL) {
    fprintf(stderr,"In %s, expected a line after >%s\n",filename,subsetname);
    exit(9);
  } else {
    p = Buffer;
    while (*p != '\0' && isspace((int) *p)) {
      p++;
    }
    if (*p == '\0') {
      /* Blank line.  Interpret to mean inclusion of everything */
      return (bool *) NULL;

    } else if (*p == '+') {
      /* Skip blanks */
      p++;
      return process_inclusions(p,chromosome_iit);

    } else if (*p == '-') {
      /* Skip blanks */
      p++;
      while (*p != '\0' && !isspace((int) *p)) {
	p++;
      }
      return process_exclusions(p,chromosome_iit);

    } else {
      fprintf(stderr,"In %s, don't know how to interpret this line:\n%s\n",
	      filename,Buffer);
      exit(9);
    }
  }
}

static int
compute_nincluded (bool *includep, int nchromosomes) {
  int nincluded = 0, i;

  for (i = 0; i < nchromosomes; i++) {
    if (includep[i] == true) {
      nincluded++;
    }
  }
  return nincluded;
}

static int *
compute_new_indices (bool *includep, int nchromosomes) {
  int *newindices;
  int i, newi = 1;		/* Because index is 1-based */

  newindices = (int *) CALLOC(nchromosomes,sizeof(int));
  for (i = 0; i < nchromosomes; i++) {
    if (includep[i] == true) {
      newindices[i] = newi++;
    } else {
      newindices[i] = -1;	/* Actually, this could be 0 */
    }
  }

  return newindices;
}

static int *
compute_old_indices (bool *includep, int nincluded, int nchromosomes) {
  int *oldindices;
  int i, j = 0;

  oldindices = (int *) CALLOC(nincluded,sizeof(int));
  for (i = 0; i < nchromosomes; i++) {
    if (includep[i] == true) {
      oldindices[j++] = i + 1;
    }
  }

  return oldindices;
}


T
Chrsubset_new_single (Chrnum_T chrnum, Univ_IIT_T chromosome_iit) {
  T new = (T) MALLOC(sizeof(*new));

  new->name = Chrnum_to_string(chrnum,chromosome_iit);
  new->includep = (bool *) CALLOC(Univ_IIT_total_nintervals(chromosome_iit),sizeof(bool));
  new->includep[chrnum-1] = true;
  new->nincluded = 1;
  new->newindices = compute_new_indices(new->includep,Univ_IIT_total_nintervals(chromosome_iit));
  new->oldindices = compute_old_indices(new->includep,new->nincluded,Univ_IIT_total_nintervals(chromosome_iit));
  return new;
}


T
Chrsubset_read (char *user_chrsubsetfile, char *genomesubdir, char *fileroot, 
		char *user_chrsubsetname, Univ_IIT_T chromosome_iit) {
  T new = NULL;
  FILE *fp;
  char *filename, *subsetname;
  bool *includep;
#ifdef DEBUG
  int i;
#endif

  if (user_chrsubsetfile != NULL) {
    filename = (char *) CALLOC(strlen(user_chrsubsetfile)+1,sizeof(char));
    strcpy(filename,user_chrsubsetfile);
    fp = FOPEN_READ_TEXT(filename);
    if (fp == NULL) {
      fprintf(stderr,"The provided file chromosome subset file %s could not be read\n",
	      filename);
      exit(9);
    }

  } else {
    filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+strlen(fileroot)+
			       strlen(".chrsubset")+1,sizeof(char));
    sprintf(filename,"%s/%s.chrsubset",genomesubdir,fileroot);
    fp = FOPEN_READ_TEXT(filename);
  }

  if (fp == NULL) {
    debug(printf("Standard file doesn't exist.  Including all chromosomes\n"));
    new = (T) NULL;

  } else if (user_chrsubsetname == NULL) {
    if ((subsetname = read_header(fp,filename)) == NULL) {
      fprintf(stderr,"The chromosome subset file %s is empty\n",filename);
      exit(9);
    }

    if ((includep = process_list(fp,filename,subsetname,chromosome_iit)) == NULL) {
      new = (T) NULL;
      FREE(subsetname);
    } else {
      new = (T) MALLOC(sizeof(*new));
      new->name = subsetname;
      new->includep = includep;
      new->nincluded = compute_nincluded(includep,Univ_IIT_total_nintervals(chromosome_iit));
      new->newindices = compute_new_indices(includep,Univ_IIT_total_nintervals(chromosome_iit));
      new->oldindices = compute_old_indices(includep,new->nincluded,Univ_IIT_total_nintervals(chromosome_iit));
    }

    fclose(fp);
    debug(printf("User didn't specify a subset.  Using first list: %s\n",subsetname));
    
  } else {
    debug(printf("User specified subset %s\n",user_chrsubsetname));
    while ((subsetname = read_header(fp,filename)) != NULL &&
	   strcmp(subsetname,user_chrsubsetname)) {
      debug(printf("Skipping %s\n",subsetname));
      skip_list(fp,filename,subsetname);
      FREE(subsetname);
    }
    if (subsetname == NULL) {
      fprintf(stderr,"Unable to find subset %s in chromosome subset file %s\n",
	      user_chrsubsetname,filename);
      exit(9);
    } else if ((includep = process_list(fp,filename,subsetname,chromosome_iit)) == NULL) {
      new = (T) NULL;
      FREE(subsetname);
    } else {
      new = (T) MALLOC(sizeof(*new));
      new->name = subsetname;
      new->includep = includep;
      new->nincluded = compute_nincluded(includep,Univ_IIT_total_nintervals(chromosome_iit));
      new->newindices = compute_new_indices(includep,Univ_IIT_total_nintervals(chromosome_iit));
      new->oldindices = compute_old_indices(includep,new->nincluded,Univ_IIT_total_nintervals(chromosome_iit));
    }
    fclose(fp);
  }
  FREE(filename);

  debug(
	if (new != NULL) {
	  for (i = 0; i < Univ_IIT_total_nintervals(chromosome_iit); i++) {
	    printf(" %d: %d\n",i,new->includep[i]);
	  }
	}
	);

  return new;
}


