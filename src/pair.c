static char rcsid[] = "$Id: pair.c 53584 2011-12-02 18:51:17Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifndef HAVE_MEMCPY
# define memcpy(d,s,n) bcopy((s),(d),(n))
#endif
#ifndef HAVE_MEMMOVE
# define memmove(d,s,n) bcopy((s),(d),(n))
#endif

#include "pair.h"
#include "pairdef.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>		/* For memcpy */
#include <math.h>		/* For rint(), abs() */
#include <ctype.h>		/* For toupper */

#include "assert.h"
#include "except.h"
#include "mem.h"
#include "comp.h"
#include "complement.h"
#include "intron.h"
#include "listdef.h"
#include "intlist.h"
#include "separator.h"
#include "scores.h"
#include "segmentpos.h"
#include "translation.h"
#include "maxent.h"
#include "maxent_hr.h"
#include "mapq.h"

#ifndef PMAP
#include "substring.h"		/* For Endtype_T */
#include "stage3hr.h"
#include "samflags.h"
#endif

#ifdef GSNAP
#include "samprint.h"
#endif


#define ONEBASEDP 1		/* 1-based coordinates.  Also defined in segmentpos.c */

#define MIN_INTRONLEN 20	/* For deciding between N and D in cigar string */



/* Check for ANSI mode, which does not include rint */
#ifdef __STRICT_ANSI__
#define rint(x) floor(0.5+(x))
#endif

#define DEFAULT_MARGIN 14

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Phase information */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

/* PSL indels */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

/* Pair_fracidentity_max */
#ifdef DEBUG3
#define debug3(x) x
#else
#define debug3(x)
#endif

/* compute_md_string */
#ifdef DEBUG4
#define debug4(x) x
#else
#define debug4(x)
#endif

/* trimming */
#ifdef DEBUG8
#define debug8(x) x
#else
#define debug8(x)
#endif


/* binary search */
#ifdef DEBUG10
#define debug10(x) x
#else
#define debug10(x)
#endif



#define TRIM_MATCH_SCORE 1
static int trim_mismatch_score;
static int trim_indel_score;


void
Pair_setup (int trim_mismatch_score_in, int trim_indel_score_in) {
  trim_mismatch_score = trim_mismatch_score_in;
  trim_indel_score = trim_indel_score_in;
  return;
}



#define T Pair_T

int
Pair_querypos (T this) {
  return this->querypos;
}

int
Pair_genomepos (T this) {
  return this->genomepos;
}

char
Pair_cdna (T this) {
  return this->cdna;
}

char
Pair_comp (T this) {
  return this->comp;
}

char
Pair_genome (T this) {
  return this->genome;
}

bool
Pair_gapp (T this) {
  return this->gapp;
}

bool
Pair_shortexonp (T this) {
  return this->shortexonp;
}

void
Pair_set_shortexonp (T this) {
  this->shortexonp = true;
  return;
}


void
Pair_set_genomepos (struct Pair_T *pairarray, int npairs, Genomicpos_T chrpos,
		    Genomicpos_T genomiclength, bool watsonp) {
  int i;

  if (watsonp == true) {
    for (i = 0; i < npairs; i++) {
      pairarray[i].genomepos = chrpos + pairarray[i].genomepos;
    }
  } else {
    for (i = 0; i < npairs; i++) {
      pairarray[i].genomepos = chrpos + (genomiclength - 1) - pairarray[i].genomepos;
    }
  }
  return;
}


List_T
Pair_protect (List_T pairs) {
  List_T p;
  T pair;

  for (p = pairs; p != NULL; p = p->rest) {
    pair = (T) p->first;
    pair->protectedp = true;
  }

  return pairs;
}


/* Pairs now made in pairpool.c */
/*
T
Pair_new (int querypos, int genomepos, char cdna, char comp, char genome) {
  T new = (T) MALLOC(sizeof(*new));
  
  new->querypos = querypos;
  new->genomepos = genomepos;
  new->cdna = cdna;
  new->comp = comp;
  new->genome = genome;
  switch (comp) {
  case '>': case '<': case '=': case '~': case '.': case ')': case '(': case ']': case '[': case '#':
    new->gapp = true; break;
  default: new->gapp = false;
  }

  return new;
}
*/

/* Pairs now made in pairpool.c */
/*
void
Pair_free (T *old) {
  if (*old) {
    FREE(*old);
  }
  return;
}
*/


/* Print routines */

static char *RULER = "    .    :    .    :    .    :    .    :    .    :";
static void
print_top_ruler (FILE *fp, int n, int npairs, int margin, int wraplength) {
  fprintf(fp,"%*d ",margin,n);
  if (n + wraplength < npairs) {
    fprintf(fp,"%s\n",RULER);
  } else {
    fprintf(fp,"%.*s\n",npairs-n,RULER);
  }
  return;
}

/*
static void
print_bottom_ruler (int n, int npairs, int margin, int wraplength) {
  printf("%*s ",margin,"");
  if (n + wraplength < npairs) {
    printf("%s\n",RULER);
  } else {
    printf("%.*s\n",npairs-n,RULER);
  }
  return;
}
*/


static void
print_cdna_sequence (FILE *fp, struct T *ptr, int n, int npairs, int margin, int wraplength) {
  struct T *this;
  int i;

  this = ptr;
  fprintf(fp,"%*u ",margin,this->querypos + ONEBASEDP);
  for (i = 0; n < npairs && i < wraplength; n++, i++) {
    this = ptr++;
    putc(this->cdna,fp);
  }
  putc('\n',fp);
  return;
}

static int
find_aapos_in_line (struct T *ptr, int n, int npairs, int wraplength, 
		    bool genomep) {
  struct T *this, *last;

  if (npairs - n < wraplength) {
    last = &ptr[npairs - n - 1];
  } else {
    last = &ptr[wraplength - 1];
  }
  this = ptr;
  while (this <= last && (genomep ? this->aa_g : this->aa_e) == ' ') {
    this++;
  }

  if (this > last) {
    /* No aa found */
    return -1;
  } else {
    return this->aapos;
  }
}


static void
print_peptide (FILE *fp, struct T *ptr, int n, int npairs, int margin,
	       int wraplength, bool genomep) {
  struct T *this;
  int aapos, i;

  if ((aapos = find_aapos_in_line(ptr,n,npairs,wraplength,genomep)) < 0) {
    fprintf(fp,"%*s ",margin,"");
  } else {
    /* 4 is length of "aa.c" and "aa.g" */
    if (genomep == true) {
      fprintf(fp,"aa.g%*d ",margin-4,aapos);
    } else {
      fprintf(fp,"aa.c%*d ",margin-4,aapos);
    }
  }

  if (genomep == true) {
    for (i = 0; n < npairs && i < wraplength; n++, i++) {
      this = ptr++;
      putc(this->aa_g,fp);
    }
  } else {
    for (i = 0; n < npairs && i < wraplength; n++, i++) {
      this = ptr++;
      putc(this->aa_e,fp);
    }
  }

  putc('\n',fp);
  return;
}

static void
print_alignment (FILE *fp, struct T *ptr, int n, int npairs, bool diagnosticp, 
		 int margin, int wraplength) {
  struct T *this;
  int i;

  fprintf(fp,"%*s ",margin,"");
  for (i = 0; n < npairs && i < wraplength; n++, i++) {
    this = ptr++;
    
    if (diagnosticp == true) {
      /* Subtract 1 because dynprogindices start at +1 and -1 */
      if (this->comp == DYNPROG_MATCH_COMP) {
	if (this->dynprogindex > 0) {
	  fprintf(fp,"%c",(this->dynprogindex-1)%26+'a');
	} else if (this->dynprogindex < 0) {
	  fprintf(fp,"%c",(-this->dynprogindex-1)%26+'A');
	} else {
	  putc(DYNPROG_MATCH_COMP,fp);
	}
      } else if (this->shortexonp == true) {
	putc(DIAGNOSTIC_SHORTEXON_COMP,fp);
      } else {
	putc(this->comp,fp);
      }

    } else if (this->comp == DYNPROG_MATCH_COMP) {
      putc(MATCH_COMP,fp);
    } else if (this->comp == AMBIGUOUS_COMP) {
#ifdef PMAP
      putc(AMBIGUOUS_COMP,fp);
#else
      putc(MISMATCH_COMP,fp);
#endif
    } else if (this->comp == SHORTGAP_COMP) {
      putc(INDEL_COMP,fp);
    } else if (this->comp == EXTRAEXON_COMP) {
      putc(INTRONGAP_COMP,fp);
    } else {
      putc(this->comp,fp);
    }
  }

  putc('\n',fp);
  return;
}


static void
print_genomic_sequence (FILE *fp, struct T *ptr, int n, int npairs,
			char *chrstring, Genomicpos_T chroffset,
			int margin, int wraplength) {
  struct T *this;
  int i;
  char Buffer[100];

  this = ptr;
  if (chrstring == NULL) {
    sprintf(Buffer,"%u",chroffset+this->genomepos + ONEBASEDP);
  } else {
    sprintf(Buffer,"%s:%u",chrstring,this->genomepos + ONEBASEDP);
  }
  fprintf(fp,"%*s ",margin,Buffer);
  for (i = 0; n < npairs && i < wraplength; n++, i++) {
    this = ptr++;
    if (this->comp == EXTRAEXON_COMP) {
      putc(INTRONGAP_CHAR,fp);
    } else {
      putc(this->genome,fp);
    }
  }
  putc('\n',fp);
  return;
}

static int
compute_margin (struct T *start, struct T *end, char *chrstring,
		Genomicpos_T chroffset) {
  int margin;
  char Buffer[100];

  if (chrstring == NULL) {
    sprintf(Buffer,"%u",chroffset + start->genomepos + ONEBASEDP);
  } else {
    sprintf(Buffer,"%s:%u",chrstring,start->genomepos + ONEBASEDP);
  }
  margin = (int) strlen(Buffer) + 1;

  if (chrstring == NULL) {
    sprintf(Buffer,"%u",chroffset + end->genomepos + ONEBASEDP);
  } else {
    sprintf(Buffer,"%s:%u",chrstring,end->genomepos + ONEBASEDP);
  }
  if ((int) strlen(Buffer) + 1 > margin) {
    margin = (int) strlen(Buffer) + 1;
  }

  if (margin < DEFAULT_MARGIN) {
    margin = DEFAULT_MARGIN;
  }

  return margin;
}


/*
static char
intron_symbol_rev (char c) {
  switch (c) {
  case '>': return '<';
  case ')': return '(';
  case ']': return '[';
  case '<': return '>';
  case '(': return ')';
  case '[': return ']';
  default: return c;
  }
}
*/

static char complCode[128] = COMPLEMENT_LC;

static struct T *
invert_path (struct T *old, int npairs) {
  struct T *new;
  int i, j;

  new = (struct T *) MALLOC(npairs*sizeof(struct T));
  for (i = 0, j = npairs-1; i < npairs; i++, j--) {
    memcpy(&(new[j]),&(old[i]),sizeof(struct T));
    new[j].comp = complCode[(int) old[i].comp];
  }
  return new;
}

static struct T *
invert_and_revcomp_path (struct T *old, int npairs) {
  struct T *new;
  int i, j;

  new = (struct T *) MALLOC(npairs*sizeof(struct T));
  for (i = 0, j = npairs-1; i < npairs; i++, j--) {
    memcpy(&(new[j]),&(old[i]),sizeof(struct T));
    new[j].cdna = complCode[(int) old[i].cdna];
    new[j].genome = complCode[(int) old[i].genome];
    new[j].comp = complCode[(int) old[i].comp];
  }
  return new;
}


static struct T *
invert_and_revcomp_path_and_coords (struct T *old, int npairs, int querylength) {
  struct T *new;
  int i, j;

  new = (struct T *) MALLOC(npairs*sizeof(struct T));
  for (i = 0, j = npairs-1; i < npairs; i++, j--) {
    memcpy(&(new[j]),&(old[i]),sizeof(struct T));
    new[j].querypos = (querylength - 1) - old[i].querypos;
    new[j].cdna = complCode[(int) old[i].cdna];
    new[j].genome = complCode[(int) old[i].genome];
    new[j].comp = complCode[(int) old[i].comp];
  }
  return new;
}


static void
add_intronlengths (struct T *pairs, int npairs) {
  struct T *prev, *this = NULL, *ptr;
  int space, margin, i, j, k, gapstart;
  char intronstring[20], cdnabreak[20], genomicbreak[20], comp;
  int last_querypos = -1;
  Genomicpos_T last_genomepos = -1U;

  i = 0;
  while (i < npairs) {
    prev = this;
    this = &(pairs[i++]);

    if (this->extraexonp == true) {
      /* Don't add any lengths */
    } else if (this->gapp) {
      comp = this->comp;
      gapstart = i-1;
      space = 0;
      while (this->gapp) {
	this = &(pairs[i++]);
	space++;
      }

      if (comp == DUALBREAK_COMP || comp == EXTRAEXON_COMP) {
	sprintf(cdnabreak,"%d",abs(this->querypos - last_querypos)-1);
	sprintf(genomicbreak,"%d",abs(this->genomepos - last_genomepos)-1);

	margin = (space - strlen(cdnabreak))/2;
	j = gapstart;
	while (margin > 0) {
	  ptr = &(pairs[j++]);
	  margin--;
	}
	for (k = 0; k < (int) strlen(cdnabreak); k++) {
	  ptr = &(pairs[j++]);
	  ptr->cdna = cdnabreak[k];
	}

	margin = (space - strlen(genomicbreak))/2;
	j = gapstart;
	while (margin > 0) {
	  ptr = &(pairs[j++]);
	  margin--;
	}
	for (k = 0; k < (int) strlen(genomicbreak); k++) {
	  ptr = &(pairs[j++]);
	  ptr->genome = genomicbreak[k];
	}

      } else {			/* Intron */
	sprintf(intronstring,"%d",abs(this->genomepos - last_genomepos)-1);
	margin = (space - strlen(intronstring))/2;
	j = gapstart;
	while (margin > 0) {
	  ptr = &(pairs[j++]);
	  margin--;
	}
	for (k = 0; k < (int) strlen(intronstring); k++) {
	  ptr = &(pairs[j++]);
	  ptr->cdna = intronstring[k];
	}
      }
    }

    if (this->cdna != ' ') {
      last_querypos = this->querypos;
    }
    if (this->genome != ' ') {
      last_genomepos = this->genomepos;
    }
  }	
  return;
}


/* Needed to recompute translation_length in parts of chimeras */
int
Pair_translation_length (struct T *pairs, int npairs) {
  int translation_length = 0;
  int i;

  for (i = 0; i < npairs; i++) {
    if (pairs[i].aa_e == ' ') {
    } else if (pairs[i].aa_e == '*') {
    } else {
      translation_length++;
    }
  }
  return translation_length;
}


void
Pair_print_continuous (FILE *fp, struct T *pairs, int npairs, bool watsonp,
		       bool diagnosticp, bool genomefirstp, int invertmode,
		       bool nointronlenp) {
  T this;
  struct T *save = NULL, *ptr;
  int n = 0;

  if (watsonp == true) {
    ptr = pairs;
  } else if (invertmode == 0) {
    ptr = pairs;
  } else if (invertmode == 1) {
    save = ptr = invert_path(pairs,npairs);
  } else if (invertmode == 2) {
    save = ptr = invert_and_revcomp_path(pairs,npairs);
  } else {
    fprintf(stderr,"Don't recognize invert mode %d\n",invertmode);
    exit(9);
  }
  if (nointronlenp == false) {
    add_intronlengths(ptr,npairs);
  }

  if (genomefirstp == true) {
    ptr = pairs;
    for (n = 0; n < npairs; n++) {
      this = ptr++;
      putc(this->genome,fp);
    }
    putc('\n',fp);

    ptr = pairs;
    for (n = 0; n < npairs; n++) {
      this = ptr++;
      if (this->comp == MATCH_COMP) {
	putc(MATCH_COMP,fp);
      } else if (diagnosticp == false && this->comp == DYNPROG_MATCH_COMP) {
	putc(MATCH_COMP,fp);
      } else if (diagnosticp == false && this->comp == AMBIGUOUS_COMP) {
#ifdef PMAP
	putc(AMBIGUOUS_COMP,fp);
#else
	putc(MISMATCH_COMP,fp);
#endif
      } else {
	putc(this->comp,fp);
      }
    }
    putc('\n',fp);

    ptr = pairs;
    for (n = 0; n < npairs; n++) {
      this = ptr++;
      putc(this->cdna,fp);
    }
    putc('\n',fp);

  } else {
    ptr = pairs;
    for (n = 0; n < npairs; n++) {
      this = ptr++;
      putc(this->cdna,fp);
    }
    putc('\n',fp);

    ptr = pairs;
    for (n = 0; n < npairs; n++) {
      this = ptr++;
      if (this->comp == MATCH_COMP) {
	putc(MATCH_COMP,fp);
      } else if (diagnosticp == false && this->comp == DYNPROG_MATCH_COMP) {
	putc(MATCH_COMP,fp);
      } else if (diagnosticp == false && this->comp == AMBIGUOUS_COMP) {
#ifdef PMAP
	putc(AMBIGUOUS_COMP,fp);
#else
	putc(MISMATCH_COMP,fp);
#endif
      } else {
	putc(this->comp,fp);
      }
    }
    putc('\n',fp);

    ptr = pairs;
    for (n = 0; n < npairs; n++) {
      this = ptr++;
      putc(this->genome,fp);
    }
    putc('\n',fp);
  }

  if (save != NULL) {
    FREE(save);
  }
  return;
}  



void
Pair_print_continuous_byexon (FILE *fp, struct T *pairs, int npairs, bool watsonp, bool diagnosticp, int invertmode) {
  T this;
  struct T *save = NULL, *ptr;
  int i = 0, j;

  if (watsonp == true) {
    ptr = pairs;
  } else if (invertmode == 0) {
    ptr = pairs;
  } else if (invertmode == 1) {
    save = ptr = invert_path(pairs,npairs);
  } else if (invertmode == 2) {
    save = ptr = invert_and_revcomp_path(pairs,npairs);
  } else {
    fprintf(stderr,"Don't recognize invert mode %d\n",invertmode);
    exit(9);
  }

  ptr = pairs;
  while (i < npairs) {
    j = i;
    this = ptr;

    while (j < npairs && this->gapp == false) {
      putc(this->genome,fp);
      this++;
      j++;
    }
    putc('\n',fp);

    j = i;
    this = ptr;
    while (j < npairs && this->gapp == false) {
      if (this->comp == MATCH_COMP) {
	putc(MATCH_COMP,fp);
      } else if (diagnosticp == false && this->comp == DYNPROG_MATCH_COMP) {
	putc(MATCH_COMP,fp);
      } else if (diagnosticp == false && this->comp == AMBIGUOUS_COMP) {
#ifdef PMAP
	putc(AMBIGUOUS_COMP,fp);
#else
	putc(MISMATCH_COMP,fp);
#endif
      } else {
	putc(this->comp,fp);
      }
      this++;
      j++;
    }
    putc('\n',fp);

    j = i;
    this = ptr;
    while (j < npairs && this->gapp == false) {
      putc(this->cdna,fp);
      this++;
      j++;
    }
    fprintf(fp,"\n\n");

    i = j;
    while (i < npairs && this->gapp == true) {
      this++;
      i++;
    }
    ptr = this;
  }

  if (save != NULL) {
    FREE(save);
  }
  return;
}  


void
Pair_print_alignment (FILE *fp, struct T *pairs, int npairs, Chrnum_T chrnum,
		      Genomicpos_T chroffset, IIT_T chromosome_iit, bool watsonp,
		      bool diagnosticp, int invertmode, bool nointronlenp,
		      int wraplength) {
  struct T *save = NULL, *ptr;
  int n = 0, i;
  char *chrstring = NULL;
  int margin;

  if (watsonp == true) {
    ptr = pairs;

  } else if (invertmode == 0) {
    /* Given cDNA sequence, use minus genome strand */
    ptr = pairs;

  } else if (invertmode == 1) {
    /* Invert cDNA sequence, use minus genome strand */
    save = ptr = invert_path(pairs,npairs);

  } else if (invertmode == 2) {
    /* Invert cDNA sequence, use plus genome strand */
    save = ptr = invert_and_revcomp_path(pairs,npairs);

  } else {
    fprintf(stderr,"Don't recognize invert mode %d\n",invertmode);
    exit(9);
  }
  
  if (nointronlenp == false) {
    add_intronlengths(ptr,npairs);
  }
  if (chrnum != 0) {
    if (invertmode == 2) {
      chrstring = Chrnum_to_string(chrnum,chromosome_iit);
    } else {
      chrstring = Chrnum_to_string_signed(chrnum,chromosome_iit,watsonp);
    }
  }

  margin = compute_margin(&(pairs[0]),&(pairs[npairs-1]),chrstring,chroffset);

  while (n < npairs) {
    print_top_ruler(fp,n,npairs,margin,wraplength);
    print_peptide(fp,ptr,n,npairs,margin,wraplength,/*genomep*/true);
    print_genomic_sequence(fp,ptr,n,npairs,chrstring,
			   chroffset,margin,wraplength);
    print_alignment(fp,ptr,n,npairs,diagnosticp,margin,wraplength);
    print_cdna_sequence(fp,ptr,n,npairs,margin,wraplength);
    print_peptide(fp,ptr,n,npairs,margin,wraplength,/*genomep*/false);
    putc('\n',fp);
    for (i = 0; n < npairs && i < wraplength; n++, i++) {
      ptr++;
    }
  }
  if (chrstring != NULL) {
    FREE(chrstring);
  }
  if (save != NULL) {
    FREE(save);
  }
  return;
}  

void
Pair_print_pathsummary (FILE *fp, int pathnum, T start, T end, Chrnum_T chrnum,
			Genomicpos_T chroffset, IIT_T chromosome_iit, bool referencealignp, 
			IIT_T altstrain_iit, char *strain, IIT_T contig_iit, char *dbversion,
			int querylength_given, int skiplength, int trim_start, int trim_end,
			int nexons, int matches, int unknowns, int mismatches, 
			int qopens, int qindels, int topens, int tindels, int goodness,
			bool watsonp, int cdna_direction,
			int translation_start, int translation_end, int translation_length,
			int relaastart, int relaaend, bool maponlyp,
			bool diagnosticp, int stage1_genomicstart, int stage1_genomiclength,
			int stage2_source, int stage2_indexsize, double stage3_defectrate) {
  int querypos1, querypos2, den;
  double fracidentity, coverage, trimmed_coverage;
  Genomicpos_T position1, position2, chrpos1, chrpos2;
  char *refstrain, *comma1, *comma2, *chr;

  querypos1 = start->querypos;
  querypos2 = end->querypos;

  fprintf(fp,"  Path %d: ",pathnum);
  fprintf(fp,"query %d%s%d (%d bp) => ",
	 querypos1 + ONEBASEDP,SEPARATOR,querypos2 + ONEBASEDP,querypos2-querypos1+1);

  chrpos1 = start->genomepos;
  chrpos2 = end->genomepos;

  comma1 = Genomicpos_commafmt(chrpos1 + ONEBASEDP);
  comma2 = Genomicpos_commafmt(chrpos2 + ONEBASEDP);
  if (chrnum == 0) {
    if (watsonp) {
      fprintf(fp,"genome %s%s%s (%d bp)\n",
	     comma1,SEPARATOR,comma2,chrpos2-chrpos1+1);
    } else {
      fprintf(fp,"genome %s%s%s (%d bp)\n",
	     comma1,SEPARATOR,comma2,chrpos2-chrpos1-1);
    }
  } else {
    chr = Chrnum_to_string(chrnum,chromosome_iit);
    if (watsonp) {
      fprintf(fp,"genome %s:%s%s%s (%d bp)\n",chr,comma1,SEPARATOR,comma2,chrpos2-chrpos1+1);
    } else {
      fprintf(fp,"genome %s:%s%s%s (%d bp)\n",chr,comma1,SEPARATOR,comma2,chrpos2-chrpos1-1);
    }
    FREE(chr);
  }
  FREE(comma2);
  FREE(comma1);

  if (maponlyp == false) {

    if (diagnosticp == true) {
      comma1 = Genomicpos_commafmt(stage1_genomicstart + ONEBASEDP);
      comma2 = Genomicpos_commafmt(stage1_genomiclength + ONEBASEDP);
      fprintf(fp,"    Stage 1 genomicstart: %s\n",comma1);
      fprintf(fp,"    Stage 1 genomiclength: %s\n",comma2);
      FREE(comma2);
      FREE(comma1);

      /* fprintf(fp,"    Stage 2 diag runtime: %.3f sec\n",stage2_diag_runtime); */
      /* fprintf(fp,"    Stage 2 align runtime: %.3f sec\n",stage2_align_runtime); */
      fprintf(fp,"    Stage 2 source: %d\n",stage2_source);
      fprintf(fp,"    Stage 2 indexsize: %d\n",stage2_indexsize);
      /* fprintf(fp,"    Stage 3 runtime: %.3f sec\n",stage3_runtime); */
      fprintf(fp,"    Stage 3 defectrate: %f\n",stage3_defectrate);
      fprintf(fp,"    Goodness: %d\n",goodness);
    }

    fprintf(fp,"    cDNA direction: ");
    if (cdna_direction > 0) {
      fprintf(fp,"sense\n");
    } else if (cdna_direction < 0) {
      fprintf(fp,"antisense\n");
    } else {
      fprintf(fp,"indeterminate\n");
    }
  }

  if (altstrain_iit != NULL) {
    if (strain == NULL) {
      refstrain = IIT_typestring(altstrain_iit,/*straintype*/0);
      if (refstrain[0] == '\0') {
	/* Backward compatibility with old altstrain_iit */
	fprintf(fp,"    Strain: reference\n");
      } else {
	fprintf(fp,"    Strain: %s (reference)\n",refstrain);
      }
    } else {
      fprintf(fp,"    Strain: %s\n",strain);
    }
  }

  position1 = chroffset + chrpos1;
  position2 = chroffset + chrpos2;
  comma1 = Genomicpos_commafmt(position1 + ONEBASEDP);
  comma2 = Genomicpos_commafmt(position2 + ONEBASEDP);
  if (dbversion == NULL) {
    fprintf(fp,"    Genomic pos: %s%s%s",comma1,SEPARATOR,comma2);
  } else {
    fprintf(fp,"    Genomic pos: %s:%s%s%s",dbversion,comma1,SEPARATOR,comma2);
  }
  if (chrpos1 <= chrpos2) {
    fprintf(fp," (+ strand)\n");
  } else {
    fprintf(fp," (- strand)\n");
  }
  FREE(comma2);
  FREE(comma1);

  if (contig_iit != NULL) {
    if (position1 <= position2) {
      Segmentpos_print_accessions(fp,contig_iit,position1,position2,referencealignp,strain);
    } else {
      Segmentpos_print_accessions(fp,contig_iit,position2,position1,referencealignp,strain);
    }
  }
    
  if (maponlyp == false) {
    fprintf(fp,"    Number of exons: %d\n",nexons);

#ifdef PMAP
    coverage = (double) (querypos2 - querypos1 + 1)/(double) (3*(querylength_given + skiplength));
    /* coverage = (double) (matches + mismatches + qindels)/(double) (3*(querylength_given + skiplength)); */

    /* Can have coverage greater than given querylength because of added '*' at end */
    if (coverage > 1.0) {
      coverage = 1.0;
    }
#else
    /* coverage = (double) (matches + mismatches + qindels)/(double) (querylength_given + skiplength); */
    coverage = (double) (querypos2 - querypos1 + 1)/(double) (querylength_given + skiplength);
#endif
    fprintf(fp,"    Coverage: %.1f",((double) rint(1000.0*coverage))/10.0);
#ifdef PMAP
    fprintf(fp," (query length: %d aa)\n",querylength_given);
#else
    fprintf(fp," (query length: %d bp)\n",querylength_given);
    if (querypos2 + 1 > trim_end) {
      trim_end = querypos2 + 1;
    }
    if (querypos1 < trim_start) {
      trim_start = querypos1;
    }

    trimmed_coverage = (double) (querypos2 - querypos1 + 1)/(double) (trim_end - trim_start + skiplength);
    fprintf(fp,"    Trimmed coverage: %.1f",((double) rint(1000.0*trimmed_coverage))/10.0);
    fprintf(fp," (trimmed length: %d bp, trimmed region: %d..%d)",
	   trim_end-trim_start,trim_start+ONEBASEDP,trim_end-1+ONEBASEDP);
    putc('\n',fp);
#endif

    if ((den = matches + mismatches + qindels + tindels) == 0) {
      fracidentity = 1.0;
    } else {
      fracidentity = (double) matches/(double) den;
    }

    /* The definition of indels here should be consistent with Stage3_indels */
    fprintf(fp,"    Percent identity: %.1f (%d matches, %d mismatches, %d indels, %d unknowns)\n",
	   ((double) rint(1000.0*fracidentity))/10.0,matches,mismatches,qindels+tindels,unknowns);
    if (qindels + tindels > 0) {
      fprintf(fp,"    Non-intron gaps: %d openings, %d bases in cdna; %d openings, %d bases in genome\n",
	     qopens,qindels,topens,tindels);
    } 

#ifndef PMAP
    if (translation_length > 0) {
      if (cdna_direction >= 0) {
	fprintf(fp,"    Translation: %d..%d (%d aa)\n",
	       translation_start+ONEBASEDP,translation_end+ONEBASEDP,translation_length);
      } else {
	fprintf(fp,"    Translation: %d..%d (%d aa)\n",
	       translation_end+ONEBASEDP,translation_start+ONEBASEDP,translation_length);
      }
    } else if (relaastart > 0) {
      if (relaastart < relaaend) {
	fprintf(fp,"    Protein coords: %d..%d\n",relaastart,relaaend);
      } else {
	fprintf(fp,"    Protein coords: %d..%d\n",relaaend,relaastart);
      }
    }
#endif

    /* fprintf(fp,"    Defect rate (percent): %.1f\n",defect_rate*100.0); */

    /* putc('\n',fp); -- Done by caller */
  }

  return;
}


void
Pair_print_coordinates (FILE *fp, struct T *pairs, int npairs, Chrnum_T chrnum,
			Genomicpos_T chroffset, IIT_T chromosome_iit,
			bool watsonp, int invertmode) {
  T this;
  struct T *save = NULL;
  int i;
  char *chrstring = NULL;

  if (watsonp == true) {
    /* ptr = pairs; */

  } else if (invertmode == 0) {
    /* Given cDNA sequence, use minus genome strand */
    /* ptr = pairs; */

  } else if (invertmode == 1) {
    /* Invert cDNA sequence, use minus genome strand */
    save = invert_path(pairs,npairs);

  } else if (invertmode == 2) {
    /* Invert cDNA sequence, use plus genome strand */
    save = invert_and_revcomp_path(pairs,npairs);

  } else {
    fprintf(stderr,"Don't recognize invert mode %d\n",invertmode);
    exit(9);
  }
  
  if (chrnum != 0) {
    if (invertmode == 2) {
      chrstring = Chrnum_to_string(chrnum,chromosome_iit);
    } else {
      chrstring = Chrnum_to_string_signed(chrnum,chromosome_iit,watsonp);
    }
  }

  for (i = 0; i < npairs; i++) {
    this = pairs++;
    if (this->gapp == false) {
#ifdef DEBUG1
      fprintf(fp,"%d %d %c\t",this->aapos,this->aaphase_e,this->aa_e);
#else
      if (this->aaphase_e != 0) {
	fprintf(fp,"%d\t",this->aapos);
      } else {
	fprintf(fp,"%d %c\t",this->aapos,this->aa_e);
      }
#endif
      fprintf(fp,"%d %c\t",this->querypos + ONEBASEDP,this->cdna);
      if (chrstring == NULL) {
	fprintf(fp,"%u %u %c",this->genomepos + ONEBASEDP,
		chroffset + this->genomepos + ONEBASEDP,
		this->genome);
      } else {
	fprintf(fp,"%s:%u %u %c",chrstring,
		this->genomepos + ONEBASEDP,
		chroffset + this->genomepos + ONEBASEDP,
		this->genome);
      }
#ifdef DEBUG1
      fprintf(fp,"\t%d %c",this->aaphase_g,this->aa_g);
#else
      if (this->aaphase_g != 0) {
	fprintf(fp,"\t");
      } else {
	fprintf(fp,"\t%c",this->aa_g);
      }
#endif
      putc('\n',fp);
    }
  }

  if (chrstring != NULL) {
    FREE(chrstring);
  }
  if (save != NULL) {
    FREE(save);
  }
  return;
}  


void
Pair_dump_one (T this, bool zerobasedp) {
  printf("%p ",this);
  if (this->gapp == true && this->extraexonp == false) {
    printf("*** Gap: queryjump = %d, genomejump = %d, type: ",this->queryjump,this->genomejump);
    switch (this->comp) {
    case FWD_CANONICAL_INTRON_COMP: printf("> GT-AG"); break;
    case FWD_GCAG_INTRON_COMP: printf(") GC-AG"); break;
    case FWD_ATAC_INTRON_COMP: printf("] AT-AC"); break;
    case REV_ATAC_INTRON_COMP: printf("[ AT-AC"); break;
    case REV_GCAG_INTRON_COMP: printf("( GC-AG"); break;
    case REV_CANONICAL_INTRON_COMP: printf("< GT-AG"); break;
    case SHORTGAP_COMP: printf("~ shortgap"); break;
    case NONINTRON_COMP: printf("= nonintron"); break;
    default: printf("? unknown"); break;
    }

    printf(" donor:%f acceptor:%f",this->donor_prob,this->acceptor_prob);
    printf(" ***");

  } else {
    printf("%d %d %c ",
	   this->querypos + !zerobasedp,this->genomepos + !zerobasedp,this->cdna);

    /* Subtract 1 because dynprogindices start at +1 and -1 */
    if (this->dynprogindex > 0) {
      printf("%c%c",this->comp,(this->dynprogindex-1)%26+'a');
    } else if (this->dynprogindex < 0) {
      printf("%c%c",this->comp,(-this->dynprogindex-1)%26+'A');
    } else {
      putchar(this->comp);
    }
    printf(" %c",this->genome);
  }

  if (this->shortexonp == true) {
    printf(" shortexon");
  }

  if (this->state == BAD) {
    printf(" bad");
  }

  return;
}


/* Useful for debugging */
void
Pair_dump_list (List_T pairs, bool zerobasedp) {
  T this;
  List_T p;

  for (p = pairs; p != NULL; p = List_next(p)) {
    this = List_head(p);
    Pair_dump_one(this,zerobasedp);
    printf("\n");
  }
  printf("***End of list***\n");
  return;
}  

void
Pair_dump_array (struct T *pairs, int npairs, bool zerobasedp) {
  struct T *this;
  int i;

  for (i = 0; i < npairs; i++) {
    this = pairs++;
    printf("%d: %d %d %d %c ",
	   i,this->querypos + !zerobasedp,this->genomepos + !zerobasedp,this->aapos,
	   this->cdna);

    /* Subtract 1 because dynprogindices start at +1 and -1 */
    if (this->dynprogindex > 0) {
      printf("%c%c",this->comp,(this->dynprogindex-1)%26+'a');
    } else if (this->dynprogindex < 0) {
      printf("%c%c",this->comp,(-this->dynprogindex-1)%26+'A');
    } else {
      putchar(this->comp);
    }
    printf(" %c",this->genome);

    if (this->aaphase_g == 0 || this->aaphase_e == 0) {
      printf(" => %c %c",this->aa_g,this->aa_e);
    }
    printf("\n");
  }
  return;
}  


Genomicpos_T
Pair_genomicpos (struct T *pairs, int npairs, int querypos, bool headp) {
  struct T *this;
  int i;

  if (headp == true) {
    for (i = 0; i < npairs; i++) {
      this = pairs++;
      if (this->querypos == querypos) {
	return this->genomepos;
      } else if (this->querypos > querypos) {
	return 0;
      }
    }
  } else {
    pairs += npairs;
    for (i = npairs-1; i >= 0; --i) {
      this = --pairs;
      if (this->querypos == querypos) {
	return this->genomepos;
      } else if (this->querypos < querypos) {
	return 0;
      }
    }
  }
  return 0;
}

int
Pair_codon_changepos (struct T *pairs, int npairs, int aapos, int cdna_direction) {
  struct T *this, *start, *end;
  int changepos = 0, i, ngenome = 0, ncdna = 0;

  i = 0;
  this = pairs;
  while (i < npairs && this->aapos != aapos) {
    this++;
    i++;
  }
  start = this;

  while (i < npairs && (ngenome < 3 || ncdna < 3)) {
    if (this->gapp == false) {
      if (this->genome != ' ') {
	ngenome++;
      }
      if (this->cdna != ' ') {
	ncdna++;
      }
    }
    this++;
    i++;
  }
  end = --this;

  if (cdna_direction < 0) {
    for (this = end; this >= start; --this) {
      if (this->gapp == true) {
      } else if (this->genome == ' ') {
      } else if (this->cdna == ' ') {
      } else if (this->genome != this->cdna) {
	return changepos;
      } else {
	changepos++;
      }
    }
  } else {
    for (this = start; this <= end; this++) {
      if (this->gapp == true) {
      } else if (this->genome == ' ') {
      } else if (this->cdna == ' ') {
      } else if (this->genome != this->cdna) {
	return changepos;
      } else {
	changepos++;
      }
    }
  }

  return changepos;
}  


void
Pair_check_list (List_T pairs) {
  T this;
  List_T p;
  int prev_querypos, prev_genomepos;

  if (pairs == NULL) {
    return;
  } else {
    this = List_head(pairs);
    prev_querypos = this->querypos;
    prev_genomepos = this->genomepos;

    for (p = List_next(pairs); p != NULL; p = List_next(p)) {
      this = List_head(p);
      if (this->gapp == false) {
	if (this->querypos > prev_querypos) {
	  printf("Problem at querypos %d\n",this->querypos);
	}
#if 0
	/* No longer a valid check after genomepos converted to chrpos */
	if (this->genomepos > prev_genomepos) {
	  printf("Problem at genomepos %d\n",this->genomepos);
	}
#endif
	prev_querypos = this->querypos;
	prev_genomepos = this->genomepos;
      }
    }
  }
  return;
}  


bool
Pair_check_array (struct T *pairs, int npairs) {
  bool result = false;
  struct T *this;
  int prev_querypos, prev_genomepos;
  int i;

  if (npairs == 0) {
    return false;
  } else {
    this = pairs++;
    prev_querypos = this->querypos;
    prev_genomepos = this->genomepos;

    for (i = 1; i < npairs; i++) {
      this = pairs++;
      if (this->querypos < prev_querypos) {
	fprintf(stderr,"Problem at querypos %d\n",this->querypos);
	result = true;
      } else if (this->querypos - prev_querypos > 1) {
	/* Could be the result of a dual break */
	fprintf(stderr,"Jump at querypos %d\n",this->querypos);
	result = false;
      }
#if 0
      /* No longer a valid check after genomepos converted to chrpos */
      if (this->genomepos < prev_genomepos) {
	fprintf(stderr,"Problem at genomepos %d\n",this->genomepos);
	result = true;
      }
#endif
      prev_querypos = this->querypos;
      prev_genomepos = this->genomepos;
    }
  }
  return result;
}  


#if 0
static void
make_complement_buffered (char *complement, char *sequence, unsigned int length) {
  int i, j;

  /* complement = (char *) CALLOC(length+1,sizeof(char)); */
  for (i = length-1, j = 0; i >= 0; i--, j++) {
    complement[(int) j] = complCode[(int) sequence[i]];
  }
  complement[length] = '\0';
  return;
}
#endif

static void
make_complement_inplace (char *sequence, unsigned int length) {
  char temp;
  unsigned int i, j;

  for (i = 0, j = length-1; i < length/2; i++, j--) {
    temp = complCode[(int) sequence[i]];
    sequence[i] = complCode[(int) sequence[j]];
    sequence[j] = temp;
  }
  if (i == j) {
    sequence[i] = complCode[(int) sequence[i]];
  }

  return;
}


static double
donor_score (Genomicpos_T genomicpos, Genomicpos_T chroffset, bool revcomp, Genome_T genome, IIT_T chromosome_iit) {
  Genomicpos_T left;
  Chrnum_T chrnum;
  int nunknowns;
  char gbuffer[MAXENT_MAXLENGTH];
  UINT4 *genome_blocks;

  if (revcomp == false) {
    if ((genome_blocks = Genome_blocks(genome)) != NULL) {
      /* Add 1 to get from exon end to intron start */
      return Maxent_hr_donor_prob(genomicpos + 1U,chroffset);
    } else {
      left = genomicpos + 1 - DONOR_MODEL_LEFT_MARGIN; /* Add 1 to get from exon end to intron start */
      Genome_fill_buffer(&chrnum,&nunknowns,genome,left,DONOR_MODEL_LEFT_MARGIN+DONOR_MODEL_RIGHT_MARGIN+1,gbuffer,chromosome_iit);
#if 0
      printf("\n");
      printf("%s donor truestrand:+ left:%u\n",gbuffer,left);
      printf("%*s^^\n",DONOR_MODEL_LEFT_MARGIN,"");
#endif
      return Maxent_donor_prob(gbuffer);
    }

  } else {
    if ((genome_blocks = Genome_blocks(genome)) != NULL) {
      return Maxent_hr_antidonor_prob(genomicpos,chroffset);
    } else {
      left = genomicpos - DONOR_MODEL_RIGHT_MARGIN - 1;
      Genome_fill_buffer(&chrnum,&nunknowns,genome,left,DONOR_MODEL_LEFT_MARGIN+DONOR_MODEL_RIGHT_MARGIN+1,gbuffer,chromosome_iit);
      make_complement_inplace(gbuffer,DONOR_MODEL_LEFT_MARGIN+DONOR_MODEL_RIGHT_MARGIN+1);
#if 0
      printf("\n");
      printf("%s donor truestrand:- left:%u\n",gbuffer,left);
      printf("%*s^^\n",DONOR_MODEL_LEFT_MARGIN,"");
#endif
      return Maxent_donor_prob(gbuffer);
    }
  }
}


static double
acceptor_score (Genomicpos_T genomicpos, Genomicpos_T chroffset, bool revcomp, Genome_T genome, IIT_T chromosome_iit) {
  Genomicpos_T left;
  Chrnum_T chrnum;
  int nunknowns;
  char gbuffer[MAXENT_MAXLENGTH];
  UINT4 *genome_blocks;

  if (revcomp == false) {
    /* sense on plus strand, or antisense on minus strand */
    if ((genome_blocks = Genome_blocks(genome)) != NULL) {
      return Maxent_hr_acceptor_prob(genomicpos,chroffset);
    } else {
      left = genomicpos - ACCEPTOR_MODEL_LEFT_MARGIN;
      Genome_fill_buffer(&chrnum,&nunknowns,genome,left,ACCEPTOR_MODEL_LEFT_MARGIN+ACCEPTOR_MODEL_RIGHT_MARGIN+1,gbuffer,chromosome_iit);
#if 0
      printf("\n");
      printf("%s acceptor truestrand:+ left:%u\n",gbuffer,left);
      printf("%*s^^\n",ACCEPTOR_MODEL_LEFT_MARGIN-2,"");
#endif
      return Maxent_acceptor_prob(gbuffer);
    }

  } else {
    if ((genome_blocks = Genome_blocks(genome)) != NULL) {
      /* Add 1 to get from exon end to intron start */
      return Maxent_hr_antiacceptor_prob(genomicpos + 1U,chroffset);
    } else {
      left = genomicpos - ACCEPTOR_MODEL_RIGHT_MARGIN;
      Genome_fill_buffer(&chrnum,&nunknowns,genome,left,ACCEPTOR_MODEL_LEFT_MARGIN+ACCEPTOR_MODEL_RIGHT_MARGIN+1,gbuffer,chromosome_iit);
      make_complement_inplace(gbuffer,ACCEPTOR_MODEL_LEFT_MARGIN+ACCEPTOR_MODEL_RIGHT_MARGIN+1);
#if 0
      printf("\n");
      printf("%s acceptor truestrand:- left:%u\n",gbuffer,left);
      printf("%*s^^\n",ACCEPTOR_MODEL_LEFT_MARGIN-2,"");
#endif
      return Maxent_acceptor_prob(gbuffer);
    }
  }
}



static bool
unknown_base (char c) {
  switch (c) {
  case 'A': case 'C': case 'G': case 'T': case 'U':
  case 'a': case 'c': case 'g': case 't': case 'u': return false;
  default: return true;
  }
}
    
void
Pair_print_exonsummary (FILE *fp, struct T *pairs, int npairs, Chrnum_T chrnum,
			Genomicpos_T chroffset, Genome_T genome, IIT_T chromosome_iit,
			bool watsonp, int cdna_direction, bool genomefirstp, int invertmode) {
  bool in_exon = false;
  struct T *save = NULL, *ptr, *prev, *this = NULL;
  int exon_querystart = -1, exon_genomestart = -1, exon_queryend, exon_genomeend, intron_start, intron_end;
  int num = 0, den = 0, i;
  char *chrstring = NULL;
  int last_querypos = -1;
  Genomicpos_T last_genomepos = -1U;


  if (watsonp == true) {
    ptr = pairs;
  } else if (invertmode == 0) {
    ptr = pairs;
  } else if (invertmode == 1) {
    save = ptr = invert_path(pairs,npairs);
  } else if (invertmode == 2) {
    save = ptr = invert_and_revcomp_path(pairs,npairs);
  } else {
    fprintf(stderr,"Don't recognize invert mode %d\n",invertmode);
    exit(9);
  }
  
  if (chrnum != 0) {
    if (invertmode == 2) {
      chrstring = Chrnum_to_string_signed(chrnum,chromosome_iit,/*watsonp*/true);
    } else {
      chrstring = Chrnum_to_string_signed(chrnum,chromosome_iit,watsonp);
    }
  }

  debug(Pair_dump_array(pairs,npairs,/*zerobasedp*/true));

  for (i = 0; i < npairs; i++) {
    prev = this;
    this = ptr++;

    if (this->gapp) {
      if (in_exon == true) {
	exon_queryend = last_querypos + ONEBASEDP;
	exon_genomeend = last_genomepos + ONEBASEDP;
	if (watsonp) {
	  intron_start = exon_genomeend + 1;
	} else {
	  intron_start = exon_genomeend - 1;
	}
	if (genomefirstp == true) {
	  fprintf(fp,"    ");
	  if (chrnum == 0) {
	    fprintf(fp,"%u-%u",chroffset+exon_genomestart,chroffset+exon_genomeend);
	  } else {
	    fprintf(fp,"%s:%d-%d",chrstring,exon_genomestart,exon_genomeend);
	  }
	  fprintf(fp,"  (%d-%d)",exon_querystart,exon_queryend);
	} else {
	  fprintf(fp,"    %d-%d",exon_querystart,exon_queryend);
	  fprintf(fp,"  ");
	  if (chrnum == 0) {
	    fprintf(fp,"(%u-%u)",chroffset+exon_genomestart,chroffset+exon_genomeend);
	  } else {
	    fprintf(fp,"(%s:%d-%d)",chrstring,exon_genomestart,exon_genomeend);
	  }
	}
	if (den == 0) {
	  fprintf(fp,"   %d%%",100);
	} else {
	  fprintf(fp,"   %d%%",(int) floor(100.0*(double) num/(double) den));
	}
	if (this->comp == FWD_CANONICAL_INTRON_COMP) {
	  fprintf(fp," ->");
	  /* sensep = true; */
	} else if (this->comp == REV_CANONICAL_INTRON_COMP) {
	  fprintf(fp," <-");
	  /* sensep = false; */
	} else if (this->comp == FWD_GCAG_INTRON_COMP) {
	  fprintf(fp," -)");
	  /* sensep = true; */
	} else if (this->comp == REV_GCAG_INTRON_COMP) {
	  fprintf(fp," (-");
	  /* sensep = false; */
	} else if (this->comp == FWD_ATAC_INTRON_COMP) {
	  fprintf(fp," -]");
	  /* sensep = true; */
	} else if (this->comp == REV_ATAC_INTRON_COMP) {
	  fprintf(fp," [-");
	  /* sensep = false; */
	} else if (this->comp == NONINTRON_COMP) {
	  fprintf(fp," ==");
	  /* sensep = true; */
	} else {
	  fprintf(fp," ##");
	  /* sensep = true; */
	}
	in_exon = false;
      }
    } else if (this->comp == INTRONGAP_COMP) {
      /* Do nothing */
    } else {
      /* Remaining possibilities are MATCH_COMP, DYNPROG_MATCH_COMP, AMBIGUOUS_COMP, INDEL_COMP, 
	 SHORTGAP_COMP, or MISMATCH_COMP */
      if (in_exon == false) {
	exon_querystart = this->querypos + ONEBASEDP;
	exon_genomestart = this->genomepos + ONEBASEDP;
	if (watsonp) {
	  intron_end = exon_genomestart - 1;
	} else {
	  intron_end = exon_genomestart + 1;
	}
	if (i > 0) {
	  if (intron_end > intron_start) {
	    fprintf(fp,"   ...%d...",intron_end - intron_start + 1);
	  } else {
	    fprintf(fp,"   ...%d...",intron_start - intron_end + 1);
	  }

	  if (exon_querystart > exon_queryend + 1) {
	    fprintf(fp,"   ***query_skip:%d***",exon_querystart-(exon_queryend+1));
	  }

	  if (genome != NULL) {
	    if (cdna_direction >= 0) {
	      fprintf(fp,"  %.3f, %.3f",
		      donor_score(chroffset+exon_genomeend-1,chroffset,!watsonp,genome,chromosome_iit),
		      acceptor_score(chroffset+exon_genomestart-1,chroffset,!watsonp,genome,chromosome_iit));
	    } else {
	      fprintf(fp,"  %.3f, %.3f",
		      acceptor_score(chroffset+exon_genomeend-1,chroffset,watsonp,genome,chromosome_iit),
		      donor_score(chroffset+exon_genomestart-1,chroffset,watsonp,genome,chromosome_iit));
	    }
	  }

	  putc('\n',fp);
	}
	num = den = 0;
	in_exon = true;
      }
      if (this->comp == INDEL_COMP || this->comp == SHORTGAP_COMP) {
	/* Previously not counted in numerator or denominator */
	den++;
#ifndef PMAP
      } else if (unknown_base(this->cdna) || unknown_base(this->genome)) {
	/* Comp must be a space */
	/* Don't count in numerator or denominator */
#endif
      } else {
	den++;
	if (this->comp == MATCH_COMP || this->comp == DYNPROG_MATCH_COMP) {
	  num++;
	} else if (this->comp == AMBIGUOUS_COMP) {
#ifdef PMAP
	  num++;
#else
	  den--;
#endif
	}
      }
    }

    if (this->cdna != ' ') {
      last_querypos = this->querypos;
    }
    if (this->genome != ' ') {
      last_genomepos = this->genomepos;
    }
  }

  prev = this;
  exon_queryend = last_querypos + ONEBASEDP;
  exon_genomeend = last_genomepos + ONEBASEDP;
  if (genomefirstp == true) {
    fprintf(fp,"    ");
    if (chrnum == 0) {
      fprintf(fp,"%u-%u",chroffset+exon_genomestart,chroffset+exon_genomeend);
    } else {
      fprintf(fp,"%s:%d-%d",chrstring,exon_genomestart,exon_genomeend);
    }
    fprintf(fp,"  (%d-%d)",exon_querystart,exon_queryend);
  } else {
    fprintf(fp,"    %d-%d",exon_querystart,exon_queryend);
    fprintf(fp,"  ");
    if (chrnum == 0) {
      fprintf(fp,"(%u-%u)",chroffset+exon_genomestart,chroffset+exon_genomeend);
    } else {
      fprintf(fp,"(%s:%d-%d)",chrstring,exon_genomestart,exon_genomeend);
    }
  }
  if (den == 0) {
    fprintf(fp,"   %d%%",100);
  } else {
    fprintf(fp,"   %d%%",(int) floor(100.0*(double) num/(double) den));
  }
  fprintf(fp,"\n\n");

  if (chrstring != NULL) {
    FREE(chrstring);
  }
  if (save != NULL) {
    FREE(save);
  }

  return;
}

/* Tokens used by compressed and gff3 formats */

static void
print_tokens_compressed (FILE *fp, List_T tokens) {
  List_T p;
  int tokencount = 1;
  char *token, *lasttoken = NULL;
  
  for (p = tokens; p != NULL; p = List_next(p)) {
    token = (char *) List_head(p);
    if (lasttoken == NULL) {
      fprintf(fp,"\t%s",token);
      lasttoken = token;
    } else if (!strcmp(token,lasttoken)) {
      tokencount++;
    } else {
      if (tokencount > 1) {
	fprintf(fp,"!%d",tokencount);
      }
      fprintf(fp," %s",token);
      lasttoken = token;
      tokencount = 1;
    }
  }
  if (tokencount > 1) {
    fprintf(fp,"!%d",tokencount);
  }

  for (p = tokens; p != NULL; p = List_next(p)) {
    token = (char *) List_head(p);
    FREE(token);
  }

  return;
}

static void
print_tokens_gff3 (FILE *fp, List_T tokens) {
  List_T p;
  char *token;
  
  if (tokens != NULL) {
    p = tokens;
    token = (char *) List_head(p);
    fprintf(fp,"%s",token);

    for (p = List_next(p); p != NULL; p = List_next(p)) {
      token = (char *) List_head(p);
      fprintf(fp," %s",token);
    }
  }

  for (p = tokens; p != NULL; p = List_next(p)) {
    token = (char *) List_head(p);
    FREE(token);
  }

  return;
}

static List_T
push_token (List_T tokens, char *token) {
  char *copy;

  copy = (char *) CALLOC(strlen(token)+1,sizeof(char));
  strcpy(copy,token);
  return List_push(tokens,(void *) copy);
}


/* Definition of GFF3 format is at http://song.sourceforge.net/gff3.shtml */

static void
print_gff3_gene (FILE *fp, int pathnum, char *sourcename, char *accession, char *chrstring, Genomicpos_T start_genomepos, 
		 Genomicpos_T end_genomepos, bool watsonp, int cdna_direction) {

  fprintf(fp,"%s\t",chrstring);	/* 1: seqid */
  fprintf(fp,"%s\t",sourcename);	/* 2: source */
  fprintf(fp,"gene\t");		/* 3: type */

  if (start_genomepos < end_genomepos) {
    fprintf(fp,"%u\t%u\t",start_genomepos,end_genomepos); /* 4,5: start, end */
  } else {
    fprintf(fp,"%u\t%u\t",end_genomepos,start_genomepos); /* 4,5: start, end */
  }

  fprintf(fp,".\t");		/* 6: score */

  if (watsonp == true) {
    if (cdna_direction >= 0) {
      fprintf(fp,"+\t");
    } else {
      fprintf(fp,"-\t");
    }
  } else {
    if (cdna_direction >= 0) {
      fprintf(fp,"-\t");		/* 7: strand */
    } else {
      fprintf(fp,"+\t");
    }
  }

  fprintf(fp,".\t");		/* 8: phase */

  /* 9: features */
  fprintf(fp,"ID=%s.path%d;Name=%s\n",accession,pathnum,accession);

  return;
}

static void
print_gff3_mrna (FILE *fp, int pathnum, char *sourcename, char *accession, char *chrstring, Genomicpos_T start_genomepos, 
		 Genomicpos_T end_genomepos, int querylength_given, int skiplength,
		 int matches, int mismatches, int qindels, int tindels, 
		 bool watsonp, int cdna_direction) {
  int den;
  double coverage, fracidentity;

  fprintf(fp,"%s\t",chrstring);	/* 1: seqid */
  fprintf(fp,"%s\t",sourcename);	/* 2: source */
  fprintf(fp,"mRNA\t");		/* 3: type */
  if (start_genomepos < end_genomepos) {
    fprintf(fp,"%u\t%u\t",start_genomepos,end_genomepos); /* 4,5: start, end */
  } else {
    fprintf(fp,"%u\t%u\t",end_genomepos,start_genomepos); /* 4,5: start, end */
  }

  fprintf(fp,".\t");		/* 6: score */

  if (watsonp == true) {
    if (cdna_direction >= 0) {
      fprintf(fp,"+\t");
    } else {
      fprintf(fp,"-\t");
    }
  } else {
    if (cdna_direction >= 0) {
      fprintf(fp,"-\t");		/* 7: strand */
    } else {
      fprintf(fp,"+\t");
    }
  }

  fprintf(fp,".\t");		/* 8: phase */

  /* 9: features */
  fprintf(fp,"ID=%s.mrna%d;Name=%s;Parent=%s.path%d;",
	 accession,pathnum,accession,accession,pathnum);

#ifdef PMAP
    coverage = (double) (matches + mismatches + qindels)/(double) (3*(querylength_given + skiplength));
#else
    coverage = (double) (matches + mismatches + qindels)/(double) (querylength_given + skiplength);
#endif
  fprintf(fp,"Coverage=%.1f;",((double) rint(1000.0*coverage))/10.0);

  if ((den = matches + mismatches + qindels + tindels) == 0) {
    fracidentity = 1.0;
  } else {
    fracidentity = (double) matches/(double) den;
  }
  fprintf(fp,"Identity=%.1f",((double) rint(1000.0*fracidentity))/10.0);

  putc('\n',fp);

  return;
}


static void
print_gff3_exon (FILE *fp, int exonno, int pathnum, char *sourcename, char *accession, char *chrstring,
		 Genomicpos_T exon_genomestart, Genomicpos_T exon_genomeend,
		 int exon_querystart, int exon_queryend, bool watsonp, int cdna_direction,
		 int pctidentity) {

  fprintf(fp,"%s\t",chrstring);	/* 1: seqid */
  fprintf(fp,"%s\t",sourcename);	/* 2: source */
  fprintf(fp,"exon\t");		/* 3: type */
  if (exon_genomestart < exon_genomeend) {
    fprintf(fp,"%u\t%u\t",exon_genomestart,exon_genomeend); /* 4,5: start, end */
  } else {
    fprintf(fp,"%u\t%u\t",exon_genomeend,exon_genomestart); /* 4,5: start, end */
  }
  fprintf(fp,"%d\t",pctidentity);	/* 6: score */

  if (watsonp == true) {
    if (cdna_direction >= 0) {
      fprintf(fp,"+\t");
    } else {
      fprintf(fp,"-\t");
    }
  } else {
    if (cdna_direction >= 0) {
      fprintf(fp,"-\t");		/* 7: strand */
    } else {
      fprintf(fp,"+\t");
    }
  }

  fprintf(fp,".\t");		/* 8: phase */

  /* 9: features */
  fprintf(fp,"ID=%s.mrna%d.exon%d;",accession,pathnum,exonno);
  fprintf(fp,"Name=%s;",accession);
  fprintf(fp,"Parent=%s.mrna%d;",accession,pathnum);
  if (cdna_direction >= 0) {
    fprintf(fp,"Target=%s %d %d +\n",accession,exon_querystart,exon_queryend);
  } else {
    fprintf(fp,"Target=%s %d %d -\n",accession,exon_queryend,exon_querystart);
  }

  return;
}

static void
print_gff3_cds (FILE *fp, int cdsno, int pathnum, char *sourcename, char *accession, char *chrstring,
		Genomicpos_T cds_genomestart, Genomicpos_T cds_genomeend,
		int cds_querystart, int cds_queryend, bool watsonp, int cdna_direction,
		int pctidentity, int cds_phase) {

  fprintf(fp,"%s\t",chrstring);	/* 1: seqid */
  fprintf(fp,"%s\t",sourcename);	/* 2: source */
  fprintf(fp,"CDS\t");		/* 3: type */
  if (cds_genomestart < cds_genomeend) {
    fprintf(fp,"%u\t%u\t",cds_genomestart,cds_genomeend); /* 4,5: start, end */
  } else {
    fprintf(fp,"%u\t%u\t",cds_genomeend,cds_genomestart); /* 4,5: start, end */
  }
  fprintf(fp,"%d\t",pctidentity);	/* 6: score */

  if (watsonp == true) {
    if (cdna_direction >= 0) {
      fprintf(fp,"+\t");
    } else {
      fprintf(fp,"-\t");
    }
  } else {
    if (cdna_direction >= 0) {
      fprintf(fp,"-\t");		/* 7: strand */
    } else {
      fprintf(fp,"+\t");
    }
  }

  fprintf(fp,"%d\t",cds_phase);	/* 8: phase */

  /* 9: features */
  fprintf(fp,"ID=%s.mrna%d.cds%d;",accession,pathnum,cdsno);
  fprintf(fp,"Name=%s;",accession);
  fprintf(fp,"Parent=%s.mrna%d;",accession,pathnum);
  if (cdna_direction >= 0) {
    fprintf(fp,"Target=%s %d %d +\n",accession,cds_querystart,cds_queryend);
  } else {
    fprintf(fp,"Target=%s %d %d -\n",accession,cds_queryend,cds_querystart);
  }

  return;
}


static void
print_gff3_cdna_match (FILE *fp, int pathnum, char *sourcename, char *accession, char *chrstring,
		       Genomicpos_T exon_genomestart, Genomicpos_T exon_genomeend,
		       int exon_querystart, int exon_queryend, bool watsonp,
		       int pctidentity, List_T tokens) {
  
  fprintf(fp,"%s\t",chrstring);	/* 1: seqid */
  fprintf(fp,"%s\t",sourcename);	/* 2: source */
  fprintf(fp,"cDNA_match\t");		/* 3: type */
  if (exon_genomestart < exon_genomeend) {
    fprintf(fp,"%u\t%u\t",exon_genomestart,exon_genomeend); /* 4,5: start, end */
  } else {
    fprintf(fp,"%u\t%u\t",exon_genomeend,exon_genomestart); /* 4,5: start, end */
  }
  fprintf(fp,"%d\t",pctidentity);	/* 6: score */

  /* 7: strand */
  if (watsonp == true) {
    fprintf(fp,"+\t");
  } else {
    fprintf(fp,"-\t");
  }

  fprintf(fp,".\t");		/* 8: phase */

  /* 9: features */
  fprintf(fp,"ID=%s.path%d;",accession,pathnum);
  fprintf(fp,"Name=%s;",accession);
  fprintf(fp,"Target=%s %d %d;Gap=",accession,exon_querystart,exon_queryend);
  print_tokens_gff3(fp,tokens);
  putc('\n',fp);

  return;
}


static char
strand_char (int strand) {
  switch (strand) {
    case  1: return '+';
    case -1: return '-';
    case  0: return '?';
    default: return '.';
  }
}


static void
print_gff3_est_match (FILE *fp, int pathnum, char *sourcename, char *accession, char *chrstring,
		      Genomicpos_T exon_genomestart, Genomicpos_T exon_genomeend,
		      int exon_querystart, int exon_queryend,
		      int querylength_given, int skiplength, int matches, int mismatches, int qindels, int tindels,
		      bool watsonp, int cdna_direction, int pctidentity, List_T tokens) {
  int feature_strand, target_strand;
  double coverage, fracidentity;
  int den;

  fprintf(fp,"%s\t",chrstring);	/* 1: seqid */
  fprintf(fp,"%s\t",sourcename);	/* 2: source */
  fprintf(fp,"EST_match\t");	/* 3: type */
  if (exon_genomestart < exon_genomeend) {
    fprintf(fp,"%u\t%u\t",exon_genomestart,exon_genomeend); /* 4,5: start, end */
  } else {
    fprintf(fp,"%u\t%u\t",exon_genomeend,exon_genomestart); /* 4,5: start, end */
  }
  fprintf(fp,"%d\t",pctidentity);	/* 6: score */

  /* 7: strand */
  feature_strand = watsonp ? cdna_direction : -cdna_direction;
  fprintf(fp,"%c\t",strand_char(feature_strand));

  fprintf(fp,".\t");		/* 8: phase */

  /* 9: features */
  fprintf(fp,"ID=%s.path%d;",accession,pathnum);
  fprintf(fp,"Name=%s;",accession);
  target_strand = cdna_direction != 0 ? cdna_direction : (watsonp ? 1 : -1);
  fprintf(fp,"Target=%s %d %d %c;Gap=",accession,exon_querystart,exon_queryend,
      strand_char(target_strand));
  print_tokens_gff3(fp,tokens);

#ifdef PMAP
  coverage = (double) (matches + mismatches + qindels)/(double) (3*(querylength_given + skiplength));
#else
  coverage = (double) (matches + mismatches + qindels)/(double) (querylength_given + skiplength);
#endif
  fprintf(fp,";Coverage=%.1f",((double) rint(1000.0*coverage))/10.0);

  if ((den = matches + mismatches + qindels + tindels) == 0) {
    fracidentity = 1.0;
  } else {
    fracidentity = (double) matches/(double) den;
  }
  fprintf(fp,";Identity=%.1f",((double) rint(1000.0*fracidentity))/10.0);

  putc('\n',fp);
}


static void
print_gff3_exons_forward (FILE *fp, struct T *pairs, int npairs, int pathnum, char *sourcename, char *accession, char *chrstring,
			  int querylength_given, int skiplength, int matches, int mismatches,
			  int qindels, int tindels, bool watsonp, int cdna_direction,
			  bool gff_introns_p, bool gff_gene_format_p, bool gff_estmatch_format_p) {
  bool in_exon = false;
  struct T *ptr, *prev, *this = NULL;
  int exon_querystart = -1, exon_genomestart = -1, exon_queryend, exon_genomeend;
  int pctidentity, num = 0, den = 0, exonno = 0, i;
  int Mlength = 0, Ilength = 0, Dlength = 0;
  List_T tokens = NULL;
  char token[10];
  int intron_start, intron_end;
#if 0
  int intronno = 0;
#endif
  int estmatch_querystart, estmatch_queryend, estmatch_genomestart, estmatch_genomeend;
  int last_querypos = -1;
  Genomicpos_T last_genomepos = -1U;

  ptr = pairs;
  for (i = 0; i < npairs; i++) {
    prev = this;
    this = ptr++;

    if (this->gapp) {
      if (in_exon == true) {
	exon_queryend = last_querypos + 1;
	exon_genomeend = last_genomepos + 1;
	if (watsonp) {
	  intron_start = exon_genomeend + 1;
	} else {
	  intron_start = exon_genomeend - 1;
	}
	
	if (den == 0) {
	  pctidentity = 100;
	} else {
	  pctidentity = (int) floor(100.0*(double) num/(double) den);
	}
	
	if (gff_gene_format_p == true) {
	  print_gff3_exon(fp,++exonno,pathnum,sourcename,accession,chrstring,exon_genomestart,exon_genomeend,
			  exon_querystart,exon_queryend,watsonp,cdna_direction,pctidentity);
	} else {
	  if (Mlength > 0) {
	    sprintf(token,"M%d",Mlength);
	    tokens = push_token(tokens,token);
	  } else if (Ilength > 0) {
	    sprintf(token,"I%d",Ilength);
	    tokens = push_token(tokens,token);
	  } else if (Dlength > 0) {
	    sprintf(token,"D%d",Dlength);
	    tokens = push_token(tokens,token);
	  }
	  if (gff_estmatch_format_p == false) {
	    tokens = List_reverse(tokens);
	    /* ++exonno; */
	    print_gff3_cdna_match(fp,pathnum,sourcename,accession,chrstring,exon_genomestart,exon_genomeend,
				  exon_querystart,exon_queryend,watsonp,pctidentity,tokens);
	    List_free(&tokens);
	  }
	}

	Mlength = Ilength = Dlength = 0;
	in_exon = false;
      }
    } else if (this->comp == INTRONGAP_COMP) {
      /* Do nothing */
    } else {
      /* Remaining possibilities are MATCH_COMP, DYNPROG_MATCH_COMP, AMBIGUOUS_COMP, INDEL_COMP, 
	 SHORTGAP_COMP, or MISMATCH_COMP */
      if (in_exon == false) {
	exon_querystart = this->querypos + 1;
	exon_genomestart = this->genomepos + 1;
	if (watsonp) {
	  intron_end = exon_genomestart - 1;
	} else {
	  intron_end = exon_genomestart + 1;
	}

	if (gff_estmatch_format_p == true && i > 0) {
	  sprintf(token,"N%u",abs(intron_end - intron_start) + 1);
	  tokens = push_token(tokens,token);
	} else if (gff_introns_p == true) {
	  if (i > 0) {
#if 0
	    printf_gff3_intron(++intronno,pathnum,sourcename,accession,chrstring,?,?,intron_start,intron_end,watsonp);
#endif
	  }
	  putc('\n',fp);
	}

	num = den = 0;
	in_exon = true;
      }
      if (this->comp == INDEL_COMP || this->comp == SHORTGAP_COMP) {
	/* Gap in upper or lower sequence */
	if (gff_gene_format_p == true) {
	  /* Don't deal with tokens */
	} else if (this->genome == ' ') {
	  if (Mlength > 0) {
	    sprintf(token,"M%d",Mlength);
	    tokens = push_token(tokens,token);
	    Mlength = 0;
	  } else if (Dlength > 0) {
	    /* unlikely */
	    sprintf(token,"D%d",Dlength);
	    tokens = push_token(tokens,token);
	    Dlength = 0;
	  }
	  Ilength++;
	} else if (this->cdna == ' ') {
	  if (Mlength > 0) {
	    sprintf(token,"M%d",Mlength);
	    tokens = push_token(tokens,token);
	    Mlength = 0;
	  } else if (Ilength > 0) {
	    sprintf(token,"I%d",Ilength);
	    tokens = push_token(tokens,token);
	    Ilength = 0;
	  }
	  Dlength++;
	} else {
	  fprintf(stderr,"Error at %c%c%c\n",this->genome,this->comp,this->cdna);
	  exit(9);
	}

	/* Previously not counted in numerator or denominator */
	den++;

      } else {
	/* Count in token even if unknown base */

	if (gff_gene_format_p == true) {
	  /* Don't deal with tokens */
	} else if (Ilength > 0) {
	  sprintf(token,"I%d",Ilength);
	  tokens = push_token(tokens,token);
	  Ilength = 0;
	} else if (Dlength > 0) {
	  sprintf(token,"D%d",Dlength);
	  tokens = push_token(tokens,token);
	  Dlength = 0;
	}
	Mlength++;

#ifdef PMAP	
	den++;
	if (this->comp == MATCH_COMP || this->comp == DYNPROG_MATCH_COMP) {
	  num++;
	} else if (this->comp == AMBIGUOUS_COMP) {
	  num++;
	}
#else
	if (unknown_base(this->cdna) || unknown_base(this->genome)) {
	  /* Comp must be a space */
	  /* Don't count in numerator or denominator */
	} else {
	  den++;
	  if (this->comp == MATCH_COMP || this->comp == DYNPROG_MATCH_COMP) {
	    num++;
	  } else if (this->comp == AMBIGUOUS_COMP) {
	    den--;
	  }
	}
#endif

      }
    }

    if (this->cdna != ' ') {
      last_querypos = this->querypos;
    }
    if (this->genome != ' ') {
      last_genomepos = this->genomepos;
    }
  }

  prev = this;
  exon_queryend = last_querypos + 1;
  exon_genomeend = last_genomepos + 1;

  if (den == 0) {
    pctidentity = 100;
  } else {
    pctidentity = (int) floor(100.0*(double) num/(double) den);
  }
	
  if (gff_gene_format_p == true) {
    print_gff3_exon(fp,++exonno,pathnum,sourcename,accession,chrstring,exon_genomestart,exon_genomeend,
		    exon_querystart,exon_queryend,watsonp,cdna_direction,pctidentity);
  } else {
    if (Mlength > 0) {
      sprintf(token,"M%d",Mlength);
      tokens = push_token(tokens,token);
    } else if (Ilength > 0) {
      sprintf(token,"I%d",Ilength);
      tokens = push_token(tokens,token);
    } else if (Dlength > 0) {
      sprintf(token,"D%d",Dlength);
      tokens = push_token(tokens,token);
    }
    if (gff_estmatch_format_p == true) {
      estmatch_querystart = pairs->querypos + 1;
      estmatch_queryend = exon_queryend;
      estmatch_genomestart = pairs->genomepos + 1;
      estmatch_genomeend = exon_genomeend;
      if (watsonp) {
	tokens = List_reverse(tokens);
      }
      print_gff3_est_match(fp,pathnum,sourcename,accession,chrstring,
			   estmatch_genomestart,estmatch_genomeend,
			   estmatch_querystart,estmatch_queryend,
			   querylength_given,skiplength,matches,mismatches,qindels,tindels,
			   watsonp,cdna_direction,pctidentity,tokens);
    } else {
      tokens = List_reverse(tokens);
      /* ++exonno; */
      print_gff3_cdna_match(fp,pathnum,sourcename,accession,chrstring,exon_genomestart,exon_genomeend,
			    exon_querystart,exon_queryend,watsonp,pctidentity,tokens);
    }
    List_free(&tokens);
  }

  return;
}

static void
print_gff3_exons_backward (FILE *fp, struct T *pairs, int npairs, int pathnum, char *sourcename, char *accession, char *chrstring,
			   bool watsonp, int cdna_direction, bool gff_introns_p) {
  bool in_exon = false;
  struct T *ptr, *prev, *this = NULL;
  int exon_querystart = -1, exon_genomestart = -1, exon_queryend, exon_genomeend;
  int pctidentity, num = 0, den = 0, exonno = 0, i;
#if 0
  int intron_start, intron_end, intronno = 0;
#endif
  int last_querypos = -1;
  Genomicpos_T last_genomepos = -1U;

  ptr = &(pairs[npairs-1]);
  for (i = npairs-1; i >= 0; i--) {
    prev = this;
    this = ptr--;

    if (this->gapp) {
      if (in_exon == true) {
	exon_queryend = last_querypos + 1;
	exon_genomeend = last_genomepos + 1;
	
	if (den == 0) {
	  pctidentity = 100;
	} else {
	  pctidentity = (int) floor(100.0*(double) num/(double) den);
	}
	
	print_gff3_exon(fp,++exonno,pathnum,sourcename,accession,chrstring,exon_genomestart,exon_genomeend,
			exon_querystart,exon_queryend,watsonp,cdna_direction,pctidentity);

	in_exon = false;
      }
    } else if (this->comp == INTRONGAP_COMP) {
      /* Do nothing */
    } else {
      /* Remaining possibilities are MATCH_COMP, DYNPROG_MATCH_COMP, AMBIGUOUS_COMP, INDEL_COMP, 
	 SHORTGAP_COMP, or MISMATCH_COMP */
      if (in_exon == false) {
	exon_querystart = this->querypos + 1;
	exon_genomestart = this->genomepos + 1;

	if (gff_introns_p == true) {
	  if (i > 0) {
#if 0
	    printf_gff3_intron(++intronno,pathnum,sourcename,accession,chrstring,?,?,intron_start,intron_end,watsonp);
#endif
	  }
	  putc('\n',fp);
	}

	num = den = 0;
	in_exon = true;
      }
      if (this->comp == INDEL_COMP || this->comp == SHORTGAP_COMP) {
	/* Previously not counted in numerator or denominator */
	den++;

#ifndef PMAP
      } else if (unknown_base(this->cdna) || unknown_base(this->genome)) {
	/* Comp must be a space */
	/* Don't count in numerator or denominator */
#endif
      } else {
	den++;
	if (this->comp == MATCH_COMP || this->comp == DYNPROG_MATCH_COMP) {
	  num++;
	} else if (this->comp == AMBIGUOUS_COMP) {
#ifdef PMAP
	  num++;
#else
	  den--;
#endif
	}
      }
    }
    if (this->cdna != ' ') {
      last_querypos = this->querypos;
    }
    if (this->genome != ' ') {
      last_genomepos = this->genomepos;
    }
  }

  prev = this;
  exon_queryend = last_querypos + 1;
  exon_genomeend = last_genomepos + 1;

  if (den == 0) {
    pctidentity = 100;
  } else {
    pctidentity = (int) floor(100.0*(double) num/(double) den);
  }
	
  print_gff3_exon(fp,++exonno,pathnum,sourcename,accession,chrstring,exon_genomestart,exon_genomeend,
		  exon_querystart,exon_queryend,watsonp,cdna_direction,pctidentity);
  return;
}


static void
print_gff3_cdss_forward (FILE *fp, struct T *pairs, int npairs, int pathnum, char *sourcename, char *accession, char *chrstring,
			 bool watsonp, int cdna_direction) {
  bool in_cds = false;
  struct T *ptr, *prev, *this = NULL;
  int exon_querystart = -1, exon_genomestart = -1, exon_queryend, exon_genomeend, exon_phase;
  int pctidentity, num = 0, den = 0, cdsno = 0;
#if 0
  int intron_start, intron_end;
#endif
  int last_querypos = -1;
  Genomicpos_T last_genomepos = -1U;

  ptr = pairs;
  while (ptr < &(pairs[npairs])) {
    prev = this;
    this = ptr++;

    if (in_cds == true) {
      if (this->aaphase_e == -1) { /* was aaphase_g */
	/* End of cds */
	exon_queryend = last_querypos + 1;
	exon_genomeend = last_genomepos + 1;
	  
	if (den == 0) {
	  pctidentity = 100;
	} else {
	  pctidentity = (int) floor(100.0*(double) num/(double) den);
	}
	
	print_gff3_cds(fp,++cdsno,pathnum,sourcename,accession,chrstring,exon_genomestart,exon_genomeend,
		       exon_querystart,exon_queryend,watsonp,cdna_direction,pctidentity,exon_phase);

	in_cds = false;

      } else {
	/* Continuation of cds */
	if (this->comp == INDEL_COMP || this->comp == SHORTGAP_COMP) {
	  /* Previously not counted in numerator or denominator */
	  den++;
	  
#ifndef PMAP
	} else if (unknown_base(this->cdna) || unknown_base(this->genome)) {
	  /* Comp must be a space */
	  /* Don't count in numerator or denominator */
#endif
	} else {
	  den++;
	  if (this->comp == MATCH_COMP || this->comp == DYNPROG_MATCH_COMP) {
	    num++;
	  } else if (this->comp == AMBIGUOUS_COMP) {
#ifdef PMAP
	    num++;
#else
	    den--;
#endif
	  }
	}
      }

    } else {
      if (this->aaphase_e == -1) {
	/* Continuation of non-cds */
      } else {
	/* Start of cds */
	exon_querystart = this->querypos + 1;
	exon_phase = this->aaphase_e; /* ? was aaphase_g */
	exon_genomestart = this->genomepos + 1;

	num = den = 0;
	in_cds = true;
      }
    }
    if (this->cdna != ' ') {
      last_querypos = this->querypos;
    }
    if (this->genome != ' ') {
      last_genomepos = this->genomepos;
    }
  }

  if (in_cds == true) {
    exon_queryend = last_querypos + 1;
    exon_genomeend = last_genomepos + 1;

    if (den == 0) {
      pctidentity = 100;
    } else {
      pctidentity = (int) floor(100.0*(double) num/(double) den);
    }
	
    print_gff3_cds(fp,++cdsno,pathnum,sourcename,accession,chrstring,exon_genomestart,exon_genomeend,
		   exon_querystart,exon_queryend,watsonp,cdna_direction,pctidentity,exon_phase);
  }

  return;
}

static void
print_gff3_cdss_backward (FILE *fp, struct T *pairs, int npairs, int pathnum, char *sourcename, char *accession, char *chrstring,
			  bool watsonp, int cdna_direction) {
  bool in_cds = false;
  struct T *ptr, *prev, *this = NULL;
  int exon_querystart = -1, exon_genomestart = -1, exon_queryend, exon_genomeend, exon_phase;
  int pctidentity, num = 0, den = 0, cdsno = 0;
#if 0
  int intron_start, intron_end;
#endif
  int last_querypos = -1;
  Genomicpos_T last_genomepos = -1U;


  ptr = &(pairs[npairs-1]);
  while (ptr >= &(pairs[0])) {
    prev = this;
    this = ptr--;

    if (in_cds == true) {
      if (this->aaphase_e == -1) { /* was aaphase_g */
	/* End of cds */
	exon_queryend = last_querypos + 1;
	exon_genomeend = last_genomepos + 1;
	  
	if (den == 0) {
	  pctidentity = 100;
	} else {
	  pctidentity = (int) floor(100.0*(double) num/(double) den);
	}
	
	print_gff3_cds(fp,++cdsno,pathnum,sourcename,accession,chrstring,exon_genomestart,exon_genomeend,
		       exon_querystart,exon_queryend,watsonp,cdna_direction,pctidentity,exon_phase);

	in_cds = false;

      } else {
	/* Continuation of cds */
	if (this->comp == INDEL_COMP || this->comp == SHORTGAP_COMP) {
	  /* Previously not counted in numerator or denominator */
	  den++;

#ifndef PMAP
	} else if (unknown_base(this->cdna) || unknown_base(this->genome)) {
	  /* Comp must be a space */
	  /* Don't count in numerator or denominator */
#endif
	} else {
	  den++;
	  if (this->comp == MATCH_COMP || this->comp == DYNPROG_MATCH_COMP) {
	    num++;
	  } else if (this->comp == AMBIGUOUS_COMP) {
#ifdef PMAP
	    num++;
#else
	    den--;
#endif
	  }
	}
      }

    } else {
      if (this->aaphase_e == -1) { /* was aaphase_g */
	/* Continuation of non-cds */
      } else {
	/* Start of cds */
	exon_querystart = this->querypos + 1;
	exon_phase = this->aaphase_e; /* ? was aaphase_g */
	exon_genomestart = this->genomepos + 1;

	num = den = 0;
	in_cds = true;
      }
    }

    if (this->cdna != ' ') {
      last_querypos = this->querypos;
    }
    if (this->genome != ' ') {
      last_genomepos = this->genomepos;
    }
  }

  if (in_cds == true) {
    exon_queryend = last_querypos + 1;
    exon_genomeend = last_genomepos + 1;

    if (den == 0) {
      pctidentity = 100;
    } else {
      pctidentity = (int) floor(100.0*(double) num/(double) den);
    }

    print_gff3_cds(fp,++cdsno,pathnum,sourcename,accession,chrstring,exon_genomestart,exon_genomeend,
		   exon_querystart,exon_queryend,watsonp,cdna_direction,pctidentity,exon_phase);
  }

  return;
}


void
Pair_print_gff3 (FILE *fp, struct T *pairs, int npairs, int pathnum, char *accession, 
		 T start, T end, Chrnum_T chrnum, IIT_T chromosome_iit, Sequence_T usersegment,
		 int translation_end,
		 int querylength_given, int skiplength, int matches, int mismatches, 
		 int qindels, int tindels, bool watsonp, int cdna_direction,
		 bool gff_gene_format_p, bool gff_estmatch_format_p, char *sourcename) {
  char *chrstring = NULL;
  Genomicpos_T chrpos1, chrpos2;

  if (chrnum == 0) {
    chrstring = Sequence_accession(usersegment);
  } else {
    chrstring = Chrnum_to_string(chrnum,chromosome_iit);
  }

  if (sourcename == NULL) {
    sourcename = "NA";
  }

  if (gff_gene_format_p == true) {
    chrpos1 = start->genomepos;
    chrpos2 = end->genomepos;

    print_gff3_gene(fp,pathnum,sourcename,accession,chrstring,chrpos1+1,chrpos2+1,watsonp,cdna_direction);
    print_gff3_mrna(fp,pathnum,sourcename,accession,chrstring,chrpos1+1,chrpos2+1,
		    querylength_given,skiplength,matches,mismatches,qindels,tindels,
		    watsonp,cdna_direction);

    if (cdna_direction >= 0) {
      print_gff3_exons_forward(fp,pairs,npairs,pathnum,sourcename,accession,chrstring,
			       querylength_given,skiplength,matches,mismatches,qindels,tindels,
			       watsonp,cdna_direction,/*gff_introns_p*/false,/*gff_gene_format_p*/true,
			       /*gff_estmatch_format_p*/false);
      if (translation_end > 0) {
	print_gff3_cdss_forward(fp,pairs,npairs,pathnum,sourcename,accession,chrstring,watsonp,
				cdna_direction);
      }
    } else {
      print_gff3_exons_backward(fp,pairs,npairs,pathnum,sourcename,accession,chrstring,watsonp,
				cdna_direction,/*gff_introns_p*/false);
      if (translation_end > 0) {
	print_gff3_cdss_backward(fp,pairs,npairs,pathnum,sourcename,accession,chrstring,watsonp,
				 cdna_direction);
      }
    }

    fprintf(fp,"###\n");		/* Terminates gene format */
  } else {
    print_gff3_exons_forward(fp,pairs,npairs,pathnum,sourcename,accession,chrstring,
			     querylength_given,skiplength,matches,mismatches,qindels,tindels,
			     watsonp,cdna_direction,/*gff_introns_p*/false,/*gff_gene_format_p*/false,
			     gff_estmatch_format_p);
  }

  if (chrnum != 0) {
    FREE(chrstring);
  }

  return;
}


#ifndef PMAP
/************************************************************************
 *   GSNAP
 ************************************************************************/

/* Based on procedure in substring.c */
static void
print_splicesite_labels (FILE *fp, Chrnum_T chrnum, Genomicpos_T splicesitepos,
			 char *tag, IIT_T splicesites_iit, int *splicesites_divint_crosstable,
			 int typeint) {
  int *splicesites, nsplicesites, i;
  char *label;
  bool allocp;

  if (splicesites_iit != NULL) {
    /* Note: IIT_get_typed_signed_with_divno does not work here */
    splicesites = IIT_get_exact_multiple_with_divno(&nsplicesites,splicesites_iit,
						    splicesites_divint_crosstable[chrnum],
						    splicesitepos,splicesitepos+1U,typeint);
    if (nsplicesites == 0) {
#if 0
      /* GMAP can have novel splice sites in the middle */
      fprintf(stderr,"Supposed to have a splicesite at chrnum %d, %u..%u, type %d\n",
	      chrnum,splicesitepos,splicesitepos+1U,typeint);
#endif
    } else {
      fprintf(fp,",%s:",tag);
      label = IIT_label(splicesites_iit,splicesites[0],&allocp);
      fprintf(fp,"%s",label);
      if (allocp) FREE(label);

      for (i = 1; i < nsplicesites; i++) {
	label = IIT_label(splicesites_iit,splicesites[i],&allocp);
	fprintf(fp,"|%s",label);
	if (allocp) FREE(label);
      }
      FREE(splicesites);
    }
  }

  return;
}

static void
print_endtypes (FILE *fp,
		Endtype_T endtype1, int ntrim1, int nindels1, Genomicpos_T prev_splice_dist,
		Endtype_T endtype2, int ntrim2, int nindels2, Genomicpos_T splice_dist,
		int nmatches, int nmismatches, Chrnum_T chrnum, Genomicpos_T chroffset,
		Genomicpos_T exon_genomestart, Genomicpos_T exon_genomeend,
		bool watsonp, int cdna_direction,
		IIT_T splicesites_iit, int *splicesites_divint_crosstable,
		int donor_typeint, int acceptor_typeint) {
  double prob;
  Genomicpos_T prev_splicesitepos, splicesitepos;
  int typeint1, typeint2;

  if (endtype1 == END) {
    fprintf(fp,"start:%d",ntrim1);
  } else if (endtype1 == INS) {
    fprintf(fp,"ins:%d",nindels1);
  } else if (endtype1 == DEL) {
    fprintf(fp,"del:%d",nindels1);
  } else if (endtype1 == DON || endtype1 == AMB_DON) {
    typeint1 = donor_typeint;
    if (watsonp == true) {
      prev_splicesitepos = exon_genomestart-1U;
      if (cdna_direction > 0) {
	prob = Maxent_hr_donor_prob(chroffset+exon_genomestart,chroffset);
      } else {
	prob = Maxent_hr_antidonor_prob(chroffset+exon_genomestart-1U,chroffset);
      }
    } else {
      prev_splicesitepos = exon_genomestart;
      if (cdna_direction > 0) {
	prob = Maxent_hr_antidonor_prob(chroffset+exon_genomestart-1U,chroffset);
      } else {
	prob = Maxent_hr_donor_prob(chroffset+exon_genomestart,chroffset);
      }
    }
    fprintf(fp,"donor:%.2f",prob);
  } else if (endtype1 == ACC || endtype1 == AMB_ACC) {
    typeint1 = acceptor_typeint;
    if (watsonp == true) {
      prev_splicesitepos = exon_genomestart-1U;
      if (cdna_direction > 0) {
	prob = Maxent_hr_acceptor_prob(chroffset+exon_genomestart-1U,chroffset);
      } else {
	prob = Maxent_hr_antiacceptor_prob(chroffset+exon_genomestart,chroffset);
      }
    } else {
      prev_splicesitepos = exon_genomestart;
      if (cdna_direction > 0) {
	prob = Maxent_hr_antiacceptor_prob(chroffset+exon_genomestart,chroffset);
      } else {
	prob = Maxent_hr_acceptor_prob(chroffset+exon_genomestart-1U,chroffset);
      }
    }
    fprintf(fp,"acceptor:%.2f",prob);
  } else {
    fprintf(fp,"unknown");
  }

  fprintf(fp,"..");

  if (endtype2 == END) {
    fprintf(fp,"end:%d",ntrim2);
  } else if (endtype2 == INS) {
    fprintf(fp,"ins:%d",nindels2);
  } else if (endtype2 == DEL) {
    fprintf(fp,"del:%d",nindels2);
  } else if (endtype2 == DON || endtype2 == AMB_DON) {
    typeint2 = donor_typeint;
    if (watsonp == true) {
      splicesitepos = exon_genomeend;
      if (cdna_direction > 0) {
	prob = Maxent_hr_donor_prob(chroffset+exon_genomeend,chroffset);
      } else {
	prob = Maxent_hr_antidonor_prob(chroffset+exon_genomeend-1U,chroffset);
      }
    } else {
      splicesitepos = exon_genomeend-1U;
      if (cdna_direction > 0) {
	prob = Maxent_hr_antidonor_prob(chroffset+exon_genomeend-1U,chroffset);
      } else {
	prob = Maxent_hr_donor_prob(chroffset+exon_genomeend,chroffset);
      }
    }
    fprintf(fp,"donor:%.2f",prob);
  } else if (endtype2 == ACC || endtype2 == AMB_ACC) {
    typeint2 = acceptor_typeint;
    if (watsonp == true) {
      splicesitepos = exon_genomeend;
      if (cdna_direction > 0) {
	prob = Maxent_hr_acceptor_prob(chroffset+exon_genomeend-1U,chroffset);
      } else {
	prob = Maxent_hr_antiacceptor_prob(chroffset+exon_genomeend,chroffset);
      }
    } else {
      splicesitepos = exon_genomeend-1U;
      if (cdna_direction > 0) {
	prob = Maxent_hr_antiacceptor_prob(chroffset+exon_genomeend,chroffset);
      } else {
	prob = Maxent_hr_acceptor_prob(chroffset+exon_genomeend-1U,chroffset);
      }
    }
    fprintf(fp,"acceptor:%.2f",prob);
  } else {
    fprintf(fp,"unknown");
  }

  fprintf(fp,",matches:%d,sub:%d",nmatches,nmismatches);

  if (prev_splice_dist != 0 && splice_dist != 0) {
    /* Double introns */
    if (cdna_direction > 0) {
      fprintf(fp,",dir:sense,splice_type:consistent");
    } else {
      fprintf(fp,",dir:antisense,splice_type:consistent");
    }
    fprintf(fp,",splice_dist_1:%u,splice_dist_2:%u",prev_splice_dist,splice_dist);
    print_splicesite_labels(fp,chrnum,prev_splicesitepos,"label_1",splicesites_iit,
			    splicesites_divint_crosstable,typeint1);
    print_splicesite_labels(fp,chrnum,splicesitepos,"label_2",splicesites_iit,
			    splicesites_divint_crosstable,typeint2);

  } else if (prev_splice_dist != 0) {
    /* Prev intron */
    if (cdna_direction > 0) {
      fprintf(fp,",dir:sense,splice_type:consistent");
    } else {
      fprintf(fp,",dir:antisense,splice_type:consistent");
    }
    fprintf(fp,",splice_dist_1:%u",prev_splice_dist);
    print_splicesite_labels(fp,chrnum,prev_splicesitepos,"label_1",splicesites_iit,
			    splicesites_divint_crosstable,typeint1);
    if (endtype2 == AMB_DON || endtype2 == AMB_ACC) {
      print_splicesite_labels(fp,chrnum,splicesitepos,"label_2",splicesites_iit,
			      splicesites_divint_crosstable,typeint2);
    }

  } else if (splice_dist != 0) {
    /* Next intron */
    if (cdna_direction > 0) {
      fprintf(fp,",dir:sense,splice_type:consistent");
    } else {
      fprintf(fp,",dir:antisense,splice_type:consistent");
    }
    fprintf(fp,",splice_dist_2:%u",splice_dist);
    if (endtype1 == AMB_DON || endtype1 == AMB_ACC) {
      print_splicesite_labels(fp,chrnum,splicesitepos,"label_1",splicesites_iit,
			      splicesites_divint_crosstable,typeint1);
    }
    print_splicesite_labels(fp,chrnum,splicesitepos,"label_2",splicesites_iit,
			    splicesites_divint_crosstable,typeint2);
  } else {
    /* No introns */
    if (endtype1 == AMB_DON || endtype1 == AMB_ACC) {
      print_splicesite_labels(fp,chrnum,splicesitepos,"label_1",splicesites_iit,
			      splicesites_divint_crosstable,typeint1);
    }
    if (endtype2 == AMB_DON || endtype2 == AMB_ACC) {
      print_splicesite_labels(fp,chrnum,splicesitepos,"label_2",splicesites_iit,
			      splicesites_divint_crosstable,typeint2);
    }
  }

  return;
}




/* Based on print_pair_info in stage3hr.c */
static void
print_pair_info (FILE *fp, int insertlength, int pairscore, Pairtype_T pairtype) {
  fprintf(fp,"pair_score:%d",pairscore);
  fprintf(fp,",insert_length:%d",insertlength);

  switch (pairtype) {
  case CONCORDANT: break;
  case PAIRED_SCRAMBLE: fprintf(fp,",pairtype:scramble"); break;
  case PAIRED_INVERSION: fprintf(fp,",pairtype:inversion"); break;
  case PAIRED_TOOLONG: fprintf(fp,",pairtype:toolong"); break;
  case TRANSLOCATION: break;
  case PAIRED_UNSPECIFIED: abort();
  case UNPAIRED: abort();
  }

  return;
}



void
Pair_print_gsnap (FILE *fp, struct T *pairs_querydir, int npairs, int nsegments, bool invertedp,
		  Endtype_T start_endtype, Endtype_T end_endtype,
		  Chrnum_T chrnum, Genomicpos_T chroffset, Genomicpos_T chrhigh,
		  int querylength, bool watsonp, int cdna_direction, int score,
		  int insertlength, int pairscore, int mapq_score,
		  IIT_T chromosome_iit, IIT_T splicesites_iit,
		  int *splicesites_divint_crosstable, int donor_typeint, int acceptor_typeint) {
  bool in_exon = true;
  struct T *pairs, *ptr, *ptr0, *prev, *this = NULL;
  int exon_querystart = -1, exon_genomestart = -1, exon_queryend, exon_genomeend;
  int querypos, nmismatches, nmatches, ntrim_start, ntrim_end, indel_pos, nindels, prev_nindels, i;
  int last_querypos = -1;
  Genomicpos_T last_genomepos = -1U, splice_dist = 0U, prev_splice_dist, splicesitepos, prev_splicesitepos, pos;
  char *chr, strand, c;
  Endtype_T endtype, prev_endtype;
  bool allocp, firstp = true;

  if (invertedp == true) {
    pairs = invert_and_revcomp_path_and_coords(pairs_querydir,npairs,querylength);
    watsonp = !watsonp;
    cdna_direction = -cdna_direction;
  } else {
    pairs = pairs_querydir;
  }


  chr = IIT_label(chromosome_iit,chrnum,&allocp);
  if (watsonp == true) {
    strand = '+';
  } else {
    strand = '-';
  }

  fprintf(fp," ");		/* Beginning of GSNAP line */

  ptr = pairs;
  exon_querystart = ptr->querypos + 1;
  exon_genomestart = ptr->genomepos + 1;
  endtype = start_endtype;
  nmismatches = nmatches = 0;

  ntrim_start = ptr->querypos;

  /* Print leading trimmed characters */
  if (watsonp == true) {
    if (ntrim_start >= (exon_genomestart - 1)) {
      for (querypos = 0; querypos < ntrim_start - exon_genomestart + 1; querypos++) {
	fprintf(fp,"*");
      }
      pos = chroffset;
      for ( ; querypos < ntrim_start; querypos++) {
	c = Genome_get_char_blocks(pos++);
	fprintf(fp,"%c",tolower(c));
      }

    } else {
      pos = chroffset + (exon_genomestart - 1) - ntrim_start;
      for (querypos = 0; querypos < ntrim_start; querypos++) {
	c = Genome_get_char_blocks(pos++);
	fprintf(fp,"%c",tolower(c));
      }
    }

  } else {
    if ((pos = chroffset + (exon_genomestart - 1) + ntrim_start) >= chrhigh) {
      assert(ntrim_start - (int) (chrhigh - chroffset - exon_genomestart + 1) < querylength);
      for (querypos = 0; querypos <= ntrim_start - (int) (chrhigh - chroffset - exon_genomestart + 1); querypos++) {
	fprintf(fp,"*");
      }
      pos = chrhigh - 1;
      for ( ; querypos < ntrim_start; querypos++) {
	c = Genome_get_char_blocks(pos--);
	fprintf(fp,"%c",tolower(complCode[(int) c]));
      }

    } else {
      for (querypos = 0; querypos < ntrim_start; querypos++) {
	c = Genome_get_char_blocks(pos--);
	fprintf(fp,"%c",tolower(complCode[(int) c]));
      }
    }
  }

  i = 0;
  while (i < npairs) {
    prev = this;
    this = ptr++;
    i++;

    if (this->gapp) {
      if (in_exon == true) {
	/* SPLICE START */
	ptr0 = ptr;
	while (ptr0->gapp) {
	  ptr0++;
	}
	exon_queryend = last_querypos + 1;
	exon_genomeend = last_genomepos + 1;
	prev_endtype = endtype;

	if (prev_endtype == INS || prev_endtype == DEL) {
	  prev_nindels = nindels;
	  prev_splice_dist = 0U;
	} else if (prev_endtype == DON) {
	  prev_splice_dist = splice_dist;
	  prev_endtype = ACC;
	} else if (prev_endtype == ACC) {
	  prev_splice_dist = splice_dist;
	  prev_endtype = DON;
	} else {
	  prev_splice_dist = 0U;
	}

	if (cdna_direction > 0) {
	  endtype = DON;
	} else {
	  endtype = ACC;
	}

	if (watsonp == true) {
	  splice_dist = ptr0->genomepos - last_genomepos - 1;
	  prev_splicesitepos = exon_genomestart-1U;
	  splicesitepos = exon_genomeend;
	} else {
	  splice_dist = last_genomepos - ptr0->genomepos - 1;
	  prev_splicesitepos = exon_genomestart;
	  splicesitepos = exon_genomeend-1U;
	}


	fprintf(fp,"%c",tolower(ptr[-1].genome)); /* dinucleotide */
	fprintf(fp,"%c",tolower(ptr[0].genome));
	for (querypos = exon_queryend+2; querypos < querylength; querypos++) {
	  fprintf(fp,"-");
	}

	fprintf(fp,"\t%d..%d",exon_querystart,exon_queryend);
	fprintf(fp,"\t%c%s:%u..%u",strand,chr,exon_genomestart,exon_genomeend);
	fprintf(fp,"\t");
	print_endtypes(fp,prev_endtype,ntrim_start,prev_nindels,prev_splice_dist,
		       endtype,ntrim_end,/*nindels*/0,splice_dist,
		       nmatches,nmismatches,chrnum,chroffset,
		       exon_genomestart,exon_genomeend,watsonp,cdna_direction,
		       splicesites_iit,splicesites_divint_crosstable,
		       donor_typeint,acceptor_typeint);

	if (firstp == true) {
	  fprintf(fp,"\tsegs:%d,align_score:%d,mapq:%d",nsegments,score,mapq_score);
	  fprintf(fp,",method:gmap");
	  fprintf(fp,"\t");
	  print_pair_info(fp,insertlength,pairscore,/*pairtype*/CONCORDANT);
	  firstp = false;
	}

	nmismatches = nmatches = 0;
	fprintf(fp,"\n");

	in_exon = false;
      }
    } else if (this->comp == INTRONGAP_COMP) {
      /* May want to print dinucleotides */

    } else {
      /* Remaining possibilities are MATCH_COMP, DYNPROG_MATCH_COMP, AMBIGUOUS_COMP, INDEL_COMP, 
	 SHORTGAP_COMP, or MISMATCH_COMP */
      if (in_exon == false) {
	/* SPLICE CONTINUATION */
	exon_querystart = this->querypos + 1;
	exon_genomestart = this->genomepos + 1;

	fprintf(fp,",");
	for (querypos = 0; querypos < this->querypos - 2; querypos++) {
	  fprintf(fp,"-");
	}
	fprintf(fp,"%c",tolower(ptr[-3].genome)); /* dinucleotide */
	fprintf(fp,"%c",tolower(ptr[-2].genome));

	in_exon = true;
      }
      if (this->comp == INDEL_COMP || this->comp == SHORTGAP_COMP) {
	if (this->genome == ' ') {
	  /* INSERTION */
	  exon_queryend = last_querypos + 1;
	  exon_genomeend = last_genomepos + 1;
	  prev_endtype = endtype;
	  endtype = INS;

	  if (prev_endtype == INS || prev_endtype == DEL) {
	    prev_nindels = nindels;
	    prev_splice_dist = 0U;
	  } else if (prev_endtype == DON) {
	    prev_splice_dist = splice_dist;
	    prev_endtype = ACC;
	  } else if (prev_endtype == ACC) {
	    prev_splice_dist = splice_dist;
	    prev_endtype = DON;
	  } else {
	    prev_splice_dist = 0U;
	  }

	  indel_pos = this->querypos;
	  nindels = 0;
	  while (i < npairs && this->gapp == false && this->genome == ' ') {
	    nindels++;
	    this = ptr++;
	    i++;
	  }
	  ptr--;
	  i--;

	  /* Finish rest of this line */
	  for (querypos = exon_queryend; querypos < querylength; querypos++) {
	    fprintf(fp,"-");
	  }
	  fprintf(fp,"\t%d..%d",exon_querystart,exon_queryend);
	  fprintf(fp,"\t%c%s:%u..%u",strand,chr,exon_genomestart,exon_genomeend);
	  fprintf(fp,"\t");
	  print_endtypes(fp,prev_endtype,ntrim_start,prev_nindels,prev_splice_dist,
			 endtype,ntrim_end,nindels,/*splice_dist*/0U,
			 nmatches,nmismatches,chrnum,chroffset,
			 exon_genomestart,exon_genomeend,watsonp,cdna_direction,
			 splicesites_iit,splicesites_divint_crosstable,
			 donor_typeint,acceptor_typeint);

	  if (firstp == true) {
	    fprintf(fp,"\tsegs:%d,align_score:%d,mapq:%d",nsegments,score,mapq_score);
	    fprintf(fp,",method:gmap");
	    fprintf(fp,"\t");
	    print_pair_info(fp,insertlength,pairscore,/*pairtype*/CONCORDANT);
	    firstp = false;
	  }

	  fprintf(fp,"\n,");

	  this = ptr;
	  exon_querystart = this->querypos + 1;
	  exon_genomestart = this->genomepos + 1;
	  nmismatches = nmatches = 0;

	  /* Start of next line */
	  for (querypos = 1; querypos < exon_querystart; querypos++) {
	    fprintf(fp,"-");
	  }

	} else if (this->cdna == ' ') {
	  /* DELETION */
	  exon_queryend = last_querypos + 1;
	  exon_genomeend = last_genomepos + 1;
	  prev_endtype = endtype;
	  endtype = DEL;

	  if (prev_endtype == INS || prev_endtype == DEL) {
	    prev_nindels = nindels;
	    prev_splice_dist = 0U;
	  } else if (prev_endtype == DON) {
	    prev_endtype = ACC;
	    prev_splice_dist = splice_dist;
	  } else if (prev_endtype == ACC) {
	    prev_endtype = DON;
	    prev_splice_dist = splice_dist;
	  } else {
	    prev_splice_dist = 0U;
	  }

	  indel_pos = this->querypos;
	  nindels = 0;
	  while (i < npairs && this->gapp == false && this->cdna == ' ') {
	    fprintf(fp,"%c",tolower(this->genome));
	    nindels++;
	    this = ptr++;
	    i++;
	  }
	  ptr--;
	  i--;

	  /* Finish rest of this line */
	  for (querypos = exon_queryend + nindels; querypos < querylength; querypos++) {
	    fprintf(fp,"-");
	  }
	  fprintf(fp,"\t%d..%d",exon_querystart,exon_queryend);
	  fprintf(fp,"\t%c%s:%u..%u",strand,chr,exon_genomestart,exon_genomeend);
	  fprintf(fp,"\t");
	  print_endtypes(fp,prev_endtype,ntrim_start,prev_nindels,prev_splice_dist,
			 endtype,ntrim_end,nindels,/*splice_dist*/0U,
			 nmatches,nmismatches,chrnum,chroffset,
			 exon_genomestart,exon_genomeend,watsonp,cdna_direction,
			 splicesites_iit,splicesites_divint_crosstable,
			 donor_typeint,acceptor_typeint);

	  if (firstp == true) {
	    fprintf(fp,"\tsegs:%d,align_score:%d,mapq:%d",nsegments,score,mapq_score);
	    fprintf(fp,",method:gmap");
	    fprintf(fp,"\t");
	    print_pair_info(fp,insertlength,pairscore,/*pairtype*/CONCORDANT);
	    firstp = false;
	  }

	  fprintf(fp,"\n,");

	  this = ptr;
	  exon_querystart = this->querypos + 1;
	  exon_genomestart = this->genomepos + 1;
	  nmismatches = nmatches = 0;

	  /* Start of next line */
	  for (querypos = 1; querypos < exon_querystart; querypos++) {
	    fprintf(fp,"-");
	  }

	} else {
	  fprintf(stderr,"Error at %c%c%c\n",this->genome,this->comp,this->cdna);
	  exit(9);
	}

      } else {
	c = this->genome;
	if (this->genome == this->cdna) {
	  fprintf(fp,"%c",c);
	  nmatches++;
	} else {
	  fprintf(fp,"%c",tolower(c));
	  nmismatches++;
	}
      }
    }

    if (this->cdna != ' ') {
      last_querypos = this->querypos;
    }
    if (this->genome != ' ') {
      last_genomepos = this->genomepos;
    }
  }

  prev = this;
  exon_queryend = last_querypos + 1;
  exon_genomeend = last_genomepos + 1;
  prev_endtype = endtype;
  endtype = end_endtype;

  if (prev_endtype == INS || prev_endtype == DEL) {
    prev_nindels = nindels;
    prev_splice_dist = 0U;
  } else if (prev_endtype == DON) {
    prev_endtype = ACC;
    prev_splice_dist = splice_dist;
  } else if (prev_endtype == ACC) {
    prev_endtype = DON;
    prev_splice_dist = splice_dist;
  } else {
    prev_splice_dist = 0U;
  }

  ntrim_end = querylength - exon_queryend;

  /* Print trailing trimmed characters */
  if (watsonp == true) {
    pos = chroffset + (exon_genomeend - 1);
    if (pos + ntrim_end >= chrhigh) {
      assert((int) (chrhigh - chroffset - exon_genomeend) < querylength);
      for (i = 0; i < (int) (chrhigh - chroffset - exon_genomeend); i++) {
	c = Genome_get_char_blocks(++pos);
	fprintf(fp,"%c",tolower(c));
      }
      for ( ; i < ntrim_end; i++) {
	fprintf(fp,"*");
      }

    } else {
      for (i = 0; i < ntrim_end; i++) {
	c = Genome_get_char_blocks(++pos);
	fprintf(fp,"%c",tolower(c));
      }
    }

  } else {
    pos = chroffset + (exon_genomeend - 1);
    if (ntrim_end >= (exon_genomeend - 1)) {
      for (i = 0; i < exon_genomeend - 1; i++) {
	c = Genome_get_char_blocks(--pos);
	fprintf(fp,"%c",tolower(complCode[(int) c]));
      }
      for ( ; i < ntrim_end; i++) {
	fprintf(fp,"*");
      }

    } else {
      for (i = 0; i < ntrim_end; i++) {
	c = Genome_get_char_blocks(--pos);
	fprintf(fp,"%c",tolower(complCode[(int) c]));
      }
    }
  }

  fprintf(fp,"\t%d..%d",exon_querystart,exon_queryend);
  fprintf(fp,"\t%c%s:%u..%u",strand,chr,exon_genomestart,exon_genomeend);
  fprintf(fp,"\t");
  print_endtypes(fp,prev_endtype,ntrim_start,prev_nindels,prev_splice_dist,
		 /*endtype2*/end_endtype,ntrim_end,/*nindels*/0,/*splice_dist*/0U,
		 nmatches,nmismatches,chrnum,chroffset,
		 exon_genomestart,exon_genomeend,watsonp,cdna_direction,
		 splicesites_iit,splicesites_divint_crosstable,
		 donor_typeint,acceptor_typeint);

  if (firstp == true) {
    fprintf(fp,"\tsegs:%d,align_score:%d,mapq:%d",nsegments,score,mapq_score);
    fprintf(fp,",method:gmap");
    fprintf(fp,"\t");
    print_pair_info(fp,insertlength,pairscore,/*pairtype*/CONCORDANT);
    firstp = false;
  }

  fprintf(fp,"\n");

  if (allocp) {
    FREE(chr);
  }

  if (invertedp == true) {
    FREE(pairs);
  }

  return;
}


int
Pair_gsnap_nsegments (int *total_nmismatches, int *total_nindels, int *nintrons,
		      int *nindelbreaks, struct T *pairs, int npairs) {
  int nsegments = 0;
  bool in_exon = true;
  struct T *ptr, *ptr0, *this = NULL;
  int i;

  ptr = pairs;
  *total_nmismatches = 0;
  *total_nindels = 0;
  *nintrons = 0;
  *nindelbreaks = 0;

  i = 0;
  while (i < npairs) {
    this = ptr++;
    i++;

    if (this->gapp) {
      if (in_exon == true) {
	/* SPLICE START */
	ptr0 = ptr;
	while (ptr0->gapp) {
	  ptr0++;
	}

	(*nintrons) += 1;
	nsegments++;

	in_exon = false;
      }
    } else if (this->comp == INTRONGAP_COMP) {
      /* May want to print dinucleotides */

    } else {
      /* Remaining possibilities are MATCH_COMP, DYNPROG_MATCH_COMP, AMBIGUOUS_COMP, INDEL_COMP, 
	 SHORTGAP_COMP, or MISMATCH_COMP */
      if (in_exon == false) {
	/* SPLICE CONTINUATION */
	in_exon = true;
      }
      if (this->comp == INDEL_COMP || this->comp == SHORTGAP_COMP) {
	if (this->genome == ' ') {
	  /* INSERTION */
	  while (i < npairs && this->genome == ' ') {
	    (*total_nindels) += 1;
	    this = ptr++;
	    i++;
	  }
	  ptr--;
	  i--;

	  (*nindelbreaks) += 1;
	  nsegments++;

	} else if (this->cdna == ' ') {
	  /* DELETION */
	  while (i < npairs && this->cdna == ' ') {
	    (*total_nindels) -= 1;
	    this = ptr++;
	    i++;
	  }
	  ptr--;
	  i--;

	  (*nindelbreaks) += 1;
	  nsegments++;

	} else {
	  fprintf(stderr,"Error at %c%c%c\n",this->genome,this->comp,this->cdna);
	  exit(9);
	}

      } else {
	if (this->genome != this->cdna) {
	  (*total_nmismatches) += 1;
	}
      }
    }
  }

  nsegments++;

  return nsegments;
}



/************************************************************************
 *   SAM
 ************************************************************************/

/* Derived from print_tokens_gff3 */
static void
print_tokens_sam (FILE *fp, List_T tokens) {
  List_T p;
  char *token;
  
  for (p = tokens; p != NULL; p = List_next(p)) {
    token = (char *) List_head(p);
    fprintf(fp,"%s",token);
    FREE(token);
  }

  return;
}



/* Only for GMAP program */
static unsigned int
compute_sam_flag_nomate (int pathnum, int npaths, bool firstp, bool watsonp, bool sam_paired_p) {
  unsigned int flag = 0U;

  if (sam_paired_p == true) {
    flag |= PAIRED_READ;
    if (firstp == true) {
      flag |= FIRST_READ_P;
    } else {
      flag |= SECOND_READ_P;
    }
  }

  if (npaths == 0) {
    flag |= QUERY_UNMAPPED;
  } else if (watsonp == false) {
    flag |= QUERY_MINUSP;
  }

#if 0
  /* Will let external program decide what is primary */
  if (pathnum > 1) {
    flag |= NOT_PRIMARY;
  }
#endif

  return flag;
}


/* Modeled after Shortread_print_chopped */
static void
print_chopped (FILE *fp, char *contents, int querylength,
	       int hardclip_low, int hardclip_high) {
  int i;

  for (i = hardclip_low; i < querylength - hardclip_high; i++) {
    putc(contents[i],fp);
  }
  return;
}

static void
print_chopped_revcomp (FILE *fp, char *contents, int querylength,
		       int hardclip_low, int hardclip_high) {
  int i;

  for (i = querylength - 1 - hardclip_high; i >= hardclip_low; --i) {
    putc(complCode[(int) contents[i]],fp);
  }
  return;
}


/* Modeled after Shortread_print_quality */
static void
print_quality (FILE *fp, char *quality, int querylength,
	       int hardclip_low, int hardclip_high, int shift) {
  int i;
  int c;

  if (quality == NULL) {
    putc('*',fp);
  } else {
    for (i = hardclip_low; i < querylength - hardclip_high; i++) {
      if ((c = quality[i] + shift) <= 32) {
	fprintf(stderr,"Warning: With a quality-print-shift of %d, QC score %c becomes non-printable.  May need to specify --quality-protocol or --quality-print-shift\n",
		shift,quality[i]);
	abort();
      } else {
	putc(c,fp);
      }
    }
  }
  return;
}


static void
print_quality_revcomp (FILE *fp, char *quality, int querylength,
		       int hardclip_low, int hardclip_high, int shift) {
  int i;
  int c;

  if (quality == NULL) {
    putc('*',fp);
  } else {
    for (i = querylength - 1 - hardclip_high; i >= hardclip_low; --i) {
      if ((c = quality[i] + shift) <= 32) {
	fprintf(stderr,"Warning: With a quality-print-shift of %d, QC score %c becomes non-printable.  May need to specify --quality-protocol or --quality-print-shift\n",
		shift,quality[i]);
	abort();
      } else {
	putc(c,fp);
      }
    }
  }

  return;
}



/* Derived from print_gff3_cdna_match */
/* Assumes pairarray has been hard clipped already */
static void
print_sam_line (FILE *fp, bool firstp, char *accession, char *chrstring,
		bool watsonp, int cdna_direction, List_T cigar_tokens, List_T md_tokens,
		int nmismatches, int nindels, bool intronp, char *queryseq_ptr, char *quality_string,
		int hardclip_low, int hardclip_high, int querylength, Chimera_T chimera, int quality_shift,
		int pathnum, int npaths, int absmq_score, int second_absmq, unsigned int flag,
#ifdef GSNAP
		int pair_mapq_score, int end_mapq_score,
		Genomicpos_T chrpos, Genomicpos_T mate_chrpos, int pairedlength,
#else
		int mapq_score, struct T *pairarray, int npairs,
#endif
		char *sam_read_group_id) {

#ifndef GSNAP
  T pair;
  Genomicpos_T chrpos;
#endif


  /* 1. QNAME or Accession */
  fprintf(fp,"%s\t",accession);

  /* 2. Flags */
  fprintf(fp,"%u\t",flag);

  /* 3. RNAME or Chrstring */
  /* 4. POS or Chrlow */
  /* Taken from GMAP part of SAM_chromosomal_pos */
#ifndef GSNAP
  if (watsonp == true) {
    pair = &(pairarray[0]);
    chrpos = pair->genomepos + 1U;
  } else {
    pair = &(pairarray[npairs-1]);
    chrpos = pair->genomepos + 1U;
  }
#endif
  fprintf(fp,"%s\t%u\t",chrstring,chrpos /*+ 1U*/);


  /* 5. MAPQ or Mapping quality */
#ifdef GSNAP
  fprintf(fp,"%d\t",pair_mapq_score);
#else
  fprintf(fp,"%d\t",mapq_score);
#endif

  /* 6. CIGAR */
  print_tokens_sam(fp,cigar_tokens);

  /* 7. MRNM: Mate chr */
  /* 8. MPOS: Mate chrpos */
#ifdef GSNAP
  if (mate_chrpos == 0U) {
    fprintf(fp,"\t*\t0");
  } else {
    fprintf(fp,"\t=\t%u",mate_chrpos /* +1U*/);
  }
#else
  fprintf(fp,"\t*\t0");
#endif

  /* 9. ISIZE: Insert size */
#ifdef GSNAP
  if (mate_chrpos == 0) {
    fprintf(fp,"\t0");
  } else if (chrpos < mate_chrpos) {
    fprintf(fp,"\t%d",pairedlength);
  } else if (chrpos > mate_chrpos) {
    fprintf(fp,"\t%d",-pairedlength);
  } else if (firstp == true) {
    fprintf(fp,"\t%d",pairedlength);
  } else {
    fprintf(fp,"\t%d",-pairedlength);
  }
#else
  fprintf(fp,"\t0");
#endif

  /* 10. SEQ: queryseq and 11. QUAL: quality_scores */
  fprintf(fp,"\t");
  if (watsonp == true) {
    print_chopped(fp,queryseq_ptr,querylength,hardclip_low,hardclip_high);
    fprintf(fp,"\t");
    print_quality(fp,quality_string,querylength,hardclip_low,hardclip_high,
		  quality_shift);
  } else {
    print_chopped_revcomp(fp,queryseq_ptr,querylength,hardclip_low,hardclip_high);
    fprintf(fp,"\t");
    print_quality_revcomp(fp,quality_string,querylength,hardclip_low,hardclip_high,
			  quality_shift);
  }

  /* 12. TAGS: RG */
  if (sam_read_group_id != NULL) {
    fprintf(fp,"\t");
    fprintf(fp,"RG:Z:%s",sam_read_group_id);
  }

  /* 12. TAGS: MD string */
  fprintf(fp,"\t");
  fprintf(fp,"MD:Z:");
  print_tokens_sam(fp,md_tokens);

  /* 12. TAGS: NH */
  fprintf(fp,"\t");
  fprintf(fp,"NH:i:%d",npaths);

  /* 12. TAGS: HI */
  fprintf(fp,"\t");
  fprintf(fp,"HI:i:%d",pathnum);

  /* 12. TAGS: NM */
  fprintf(fp,"\t");
  fprintf(fp,"NM:i:%d",nmismatches + nindels);

  /* 12. TAGS: SM */
  fprintf(fp,"\t");
#ifdef GSNAP
  fprintf(fp,"SM:i:%d",end_mapq_score);
#else
  fprintf(fp,"SM:i:%d",40);
#endif

  /* 12. TAGS: XQ */
  fprintf(fp,"\t");
  fprintf(fp,"XQ:i:%d",absmq_score);

  /* 12. TAGS: X2 */
  fprintf(fp,"\t");
  fprintf(fp,"X2:i:%d",second_absmq);

  /* 12. TAGS: XS */
  if (intronp == true) {
    fprintf(fp,"\t");
    if (cdna_direction > 0) {
      /* sense */
      if (watsonp == true) {
	fprintf(fp,"XS:A:+");
      } else {
	fprintf(fp,"XS:A:-");
      }

    } else if (cdna_direction < 0) {
      /* antisense */
      if (watsonp == true) {
	fprintf(fp,"XS:A:-");
      } else {
	fprintf(fp,"XS:A:+");
      }

    } else {
      /* non-canonical */
      fprintf(fp,"XS:A:?");
    }
  }

  /* 12. TAGS: XT */
  if (chimera != NULL) {
    fprintf(fp,"\t");
    fprintf(fp,"XT:Z:");
    Chimera_print_sam_tag(fp,chimera);
  }

  putc('\n',fp);

  return;
}


static struct T *
hardclip_pairs (int *clipped_npairs, int *hardclip_low, int *hardclip_high,
		struct T *pairs, int npairs, int querylength) {
  struct T *clipped_pairs, *ptr;
  int i, starti;

  debug10(printf("Starting with %d pairs\n",npairs));

  i = 0;
  ptr = pairs;
  while (i < npairs && ptr->querypos < *hardclip_low) {
    i++;
    ptr++;
  }
  while (i < npairs && (ptr->gapp == true || ptr->cdna == ' ' || ptr->genome == ' ')) {
    i++;
    ptr++;
  }

  if (*hardclip_low > 0) {
    *hardclip_low = ptr->querypos;
  }

  starti = i;
  debug10(printf("starti is %d\n",starti));

  clipped_pairs = ptr;

  while (i < npairs && ptr->querypos < querylength - *hardclip_high) {
    i++;
    ptr++;
  }

  i--;
  ptr--;
  while (i >= starti && (ptr->gapp == true || ptr->cdna == ' ' || ptr->genome == ' ')) {
    i--;
    ptr--;
  }
  
  if (*hardclip_high > 0) {
    *hardclip_high = querylength - 1 - ptr->querypos;
  }
  *clipped_npairs = i - starti + 1;

  debug10(printf("Ending with %d pairs\n",*clipped_npairs));

  return clipped_pairs;
}



static List_T
compute_cigar (bool *intronp, int *hardclip_low, int *hardclip_high, struct T *pairs, int npairs, int querylength_given,
	       bool watsonp, int cdna_direction, int chimera_part) {
  List_T tokens = NULL, unique, p;
  char token[10], *curr_token, *last_token;
  int type, last_type;
  int Mlength = 0, Ilength = 0, Dlength = 0, length;
  bool in_exon = false, deletionp, duplicatep;
  struct T *ptr, *prev, *this = NULL;
  int exon_querystart = -1, exon_genomestart = -1, exon_queryend = -1, exon_genomeend;
  int intron_start, intron_end;
  int genome_gap, query_gap;
  int last_querypos = -1;
  Genomicpos_T last_genomepos = -1U;
  int i;

  /* *chimera_hardclip_low = *chimera_hardclip_high = 0; */
  *intronp = false;

  ptr = pairs;

  if (chimera_part == +1) {
    if (ptr->querypos > *hardclip_low) {
      if (ptr->querypos > 0) {
	/* Clip to beginning */
	*hardclip_low = ptr->querypos;
	sprintf(token,"%dH",*hardclip_low);
	tokens = push_token(tokens,token);
      }
    } else {
      if (*hardclip_low > 0) {
	/* Clip to hard clip boundary */
	sprintf(token,"%dH",*hardclip_low);
	tokens = push_token(tokens,token);
      }
    }
  } else {
    if (*hardclip_low > 0) {
      sprintf(token,"%dH",*hardclip_low);
      tokens = push_token(tokens,token);
    }
    if (ptr->querypos > (*hardclip_low)) {
      sprintf(token,"%dS",ptr->querypos - (*hardclip_low));
      tokens = push_token(tokens,token);
    }
  }

  this = (T) NULL;
  for (i = 0; i < npairs; i++) {
    prev = this;
    this = ptr++;

#if 0
    print_tokens_sam(stdout,tokens);
    printf("querypos %d, %c %c, exon %u..%u, intron %u..%u\n",
	   this->querypos,this->cdna,this->genome,exon_genomestart,exon_genomeend,
	   intron_start,intron_end);
#endif

    if (this->gapp) {
      if (in_exon == true) {
	exon_queryend = last_querypos + 1;
	exon_genomeend = last_genomepos + 1;
	if (watsonp) {
	  intron_start = exon_genomeend + 1;
	} else {
	  intron_start = exon_genomeend - 1;
	}
	
	if (Mlength > 0) {
	  sprintf(token,"%dM",Mlength);
	  tokens = push_token(tokens,token);
	} else if (Ilength > 0) {
	  sprintf(token,"%dI",Ilength);
	  tokens = push_token(tokens,token);
	} else if (Dlength > 0) {
	  sprintf(token,"%dD",Dlength);
	  tokens = push_token(tokens,token);
	}

	Mlength = Ilength = Dlength = 0;

	in_exon = false;
      }

    } else if (this->comp == INTRONGAP_COMP) {
      /* Do nothing */

    } else {
      /* Remaining possibilities are MATCH_COMP, DYNPROG_MATCH_COMP, AMBIGUOUS_COMP, INDEL_COMP, 
	 SHORTGAP_COMP, or MISMATCH_COMP */
      if (in_exon == false) {
	exon_querystart = this->querypos + 1;
	exon_genomestart = this->genomepos + 1;
	if (watsonp) {
	  intron_end = exon_genomestart - 1;
	} else {
	  intron_end = exon_genomestart + 1;
	}

	if (prev != NULL) {
	  /* Gap */
	  genome_gap = abs(intron_end - intron_start) + 1;

	  deletionp = false;
#ifdef CONVERT_INTRONS_TO_DELETIONS
	  if (cdna_direction > 0) {
	    if (prev->comp == FWD_CANONICAL_INTRON_COMP ||
		prev->comp == FWD_GCAG_INTRON_COMP ||
		prev->comp == FWD_ATAC_INTRON_COMP) {
	      sprintf(token,"%uN",genome_gap);
	      *intronp = true;
	    } else if (cigar_noncanonical_splices_p == true && genome_gap >= MIN_INTRONLEN) {
	      sprintf(token,"%uN",genome_gap);
	      *intronp = true;
	    } else {
	      sprintf(token,"%uD",genome_gap);
	      deletionp = true;
	    }
	  } else if (cdna_direction < 0) {
	    if (prev->comp == REV_CANONICAL_INTRON_COMP ||
		prev->comp == REV_GCAG_INTRON_COMP ||
		prev->comp == REV_ATAC_INTRON_COMP) {
	      sprintf(token,"%uN",genome_gap);
	      *intronp = true;
	    } else if (cigar_noncanonical_splices_p == true && genome_gap >= MIN_INTRONLEN) {
	      sprintf(token,"%uN",genome_gap);
	      *intronp = true;
	    } else {
	      sprintf(token,"%uD",genome_gap);
	      deletionp = true;
	    }
	  } else if (cigar_noncanonical_splices_p == true && genome_gap >= MIN_INTRONLEN){
	    sprintf(token,"%uN",genome_gap);
	    *intronp = true;
	  } else {
	    sprintf(token,"%uD",genome_gap);
	    deletionp = true;
	  }
#else
	  sprintf(token,"%uN",genome_gap);
	  *intronp = true;
#endif
	  tokens = push_token(tokens,token);

	  /* Check for dual gap.  Doesn't work for hard clipping. */
	  assert(exon_queryend >= 0);

	  query_gap = this->querypos - exon_queryend;
	  assert(query_gap >= 0);
	  if (query_gap > 0) {
	    if (deletionp == true) {
	      /* Put zero matches between deletion and insertion, since some programs will complain */
	      sprintf(token,"0M");
	      tokens = push_token(tokens,token);
	    }

	    sprintf(token,"%uI",query_gap);
	    tokens = push_token(tokens,token);
	  }
	}

	in_exon = true;
      }

      if (this->comp == INDEL_COMP || this->comp == SHORTGAP_COMP) {
	/* Gap in upper or lower sequence */
	if (this->genome == ' ') {
	  /* Insertion relative to genome */
	  if (Mlength > 0) {
	    sprintf(token,"%dM",Mlength);
	    tokens = push_token(tokens,token);
	    Mlength = 0;
	  } else if (Dlength > 0) {
	    /* unlikely */
	    sprintf(token,"%dD",Dlength);
	    tokens = push_token(tokens,token);
	    Dlength = 0;
	  }
	  Ilength++;
	} else if (this->cdna == ' ') {
	  /* Deletion relative to genome */
	  if (Mlength > 0) {
	    sprintf(token,"%dM",Mlength);
	    tokens = push_token(tokens,token);
	    Mlength = 0;
	  } else if (Ilength > 0) {
	    sprintf(token,"%dI",Ilength);
	    tokens = push_token(tokens,token);
	    Ilength = 0;
	  }
	  Dlength++;
	} else {
	  fprintf(stderr,"Error at %c%c%c\n",this->genome,this->comp,this->cdna);
	  exit(9);
	}

      } else {
	/* Count even if unknown base */

	if (Ilength > 0) {
	  sprintf(token,"%dI",Ilength);
	  tokens = push_token(tokens,token);
	  Ilength = 0;
	} else if (Dlength > 0) {
	  sprintf(token,"%dD",Dlength);
	  tokens = push_token(tokens,token);
	  Dlength = 0;
	}
	Mlength++;
      }
    }

    if (this != NULL) {
      if (this->cdna != ' ') {
	last_querypos = this->querypos;
      }
      if (this->genome != ' ') {
	last_genomepos = this->genomepos;
      }
    }
  }

  /* prev = this; */
  exon_queryend = last_querypos + 1;
  exon_genomeend = last_genomepos + 1;

  if (Mlength > 0) {
    sprintf(token,"%dM",Mlength);
    tokens = push_token(tokens,token);
  } else if (Ilength > 0) {
    sprintf(token,"%dI",Ilength);
    tokens = push_token(tokens,token);
  } else if (Dlength > 0) {
    sprintf(token,"%dD",Dlength);
    tokens = push_token(tokens,token);
  }


  /* Terminal clipping */
  if (chimera_part == -1) {
    if (last_querypos < querylength_given - 1 - (*hardclip_high)) {
      if (last_querypos < querylength_given - 1) {
	/* Clip to end */
	*hardclip_high = querylength_given - 1 - last_querypos;
	sprintf(token,"%dH",*hardclip_high);
	tokens = push_token(tokens,token);
      }
    } else {
      if (*hardclip_high > 0) {
	/* Clip to hard clip boundary */
	sprintf(token,"%dH",*hardclip_high);
	tokens = push_token(tokens,token);
      }
    }
  } else {
    if (last_querypos < querylength_given - 1 - (*hardclip_high)) {
      sprintf(token,"%dS",querylength_given - 1 - (*hardclip_high) - last_querypos);
      tokens = push_token(tokens,token);
    }
    if (*hardclip_high > 0) {
      sprintf(token,"%dH",*hardclip_high);
      tokens = push_token(tokens,token);
    }
  }


  /* Merge duplicates */
  last_type = ' ';
  duplicatep = false;
  length = 0;
  unique = (List_T) NULL;

  for (p = tokens; p != NULL; p = List_next(p)) {
    curr_token = (char *) List_head(p);
    type = curr_token[strlen(curr_token)-1];
    if (type == last_type) {
      length += atoi(last_token);
      FREE(last_token);
      duplicatep = true;
    } else {
      if (last_type == ' ') {
	/* Skip */
      } else if (duplicatep == false) {
	unique = List_push(unique,(void *) last_token);
      } else {
	length += atoi(last_token);
	FREE(last_token);
	sprintf(token,"%d%c",length,last_type);
	unique = push_token(unique,token);
      }
      last_type = type;
      duplicatep = false;
      length = 0;
    }
    last_token = curr_token;
  }
  if (last_type == ' ') {
    /* Skip */
  } else if (duplicatep == false) {
    unique = List_push(unique,(void *) last_token);
  } else {
    length += atoi(last_token);
    FREE(last_token);
    sprintf(token,"%d%c",length,last_type);
    unique = push_token(unique,token);
  }
  List_free(&tokens);


  /* Insert "0M" between adjacent I and D operations */
  last_type = ' ';
  tokens = (List_T) NULL;
  for (p = unique; p != NULL; p = List_next(p)) {
    curr_token = (char *) List_head(p);
    type = curr_token[strlen(curr_token)-1];
    if (last_type == 'I' && type == 'D') {
      tokens = push_token(tokens,"0M");
    } else if (last_type == 'D' && type == 'I') {
      tokens = push_token(tokens,"0M");
    }
    tokens = List_push(tokens,(void *) curr_token);
    last_type = type;
  }


  /* Return result */
  if (watsonp) {
    /* Put tokens in forward order */
    return List_reverse(tokens);
  } else {
    /* Keep tokens in reverse order */
    return tokens;
  }
}


typedef enum {IN_MATCHES, IN_MISMATCHES, IN_DELETION} MD_state_T;


#if 0
static List_T
compute_md_string_old (int *nmismatches, struct T *pairs, int npairs, bool watsonp) {
  List_T tokens = NULL;
  char token[10], *first_token;
  int nmatches = 0;
  struct T *ptr, *prev, *this = NULL;
  MD_state_T state = IN_MISMATCHES;
  int i;

  ptr = pairs;
  *nmismatches = 0;

  /* Ignore initial soft clipping */

  if (watsonp == true) {
    for (i = 0; i < npairs; i++) {
      prev = this;
      this = ptr++;

      if (this->comp == MATCH_COMP || this->comp == DYNPROG_MATCH_COMP || this->comp == AMBIGUOUS_COMP) {
	nmatches++;
	state = IN_MATCHES;

      } else if (this->comp == MISMATCH_COMP) {
	*nmismatches += 1;
	if (state == IN_MATCHES) {
	  if (nmatches > 0) {
	    sprintf(token,"%d",nmatches);
	    tokens = push_token(tokens,token);
	    nmatches = 0;
	  }

	} else if (state == IN_DELETION) {
	  tokens = push_token(tokens,"0");
	}
	state = IN_MISMATCHES;

	sprintf(token,"%c",watsonp ? this->genome : complCode[(int) this->genome]);
	tokens = push_token(tokens,token);

      } else if (this->comp == INDEL_COMP || this->comp == SHORTGAP_COMP) {
	if (this->genome == ' ') {
#if 0
	  /* Insertion relative to genome.  Ignored in MD string (but not in cigar). */
	  nmatches++;
	  state = IN_MATCHES;
#endif

	} else if (this->cdna == ' ') {
	  /* Deletion relative to genome */
	  if (state == IN_MATCHES) {
	    if (nmatches > 0) {
	      sprintf(token,"%d",nmatches);
	      tokens = push_token(tokens,token);
	      nmatches = 0;
	    }
	    tokens = push_token(tokens,"^");

	  } else if (state == IN_MISMATCHES) {
	    tokens = push_token(tokens,"^");

	  }
	  state = IN_DELETION;

	  sprintf(token,"%c",watsonp ? this->genome : complCode[(int) this->genome]);
	  tokens = push_token(tokens,token);
	}

      } else {
	/* Ignore */
      }
    }

    /* Ignore terminal soft clipping */

    if (nmatches > 0) {
      sprintf(token,"%d",nmatches);
      tokens = push_token(tokens,token);
    }

    /* Put tokens in forward order */
    tokens = List_reverse(tokens);

  } else {

    for (i = 0; i < npairs; i++) {
      prev = this;
      this = ptr++;

      if (this->comp == MATCH_COMP || this->comp == DYNPROG_MATCH_COMP || this->comp == AMBIGUOUS_COMP) {
	if (state == IN_DELETION) {
	  tokens = push_token(tokens,"^");
	}
	nmatches++;
	state = IN_MATCHES;

      } else if (this->comp == MISMATCH_COMP) {
	*nmismatches += 1;
	if (state == IN_MATCHES) {
	  if (nmatches > 0) {
	    sprintf(token,"%d",nmatches);
	    tokens = push_token(tokens,token);
	    nmatches = 0;
	  }

	} else if (state == IN_DELETION) {
	  tokens = push_token(tokens,"^");
	}
	state = IN_MISMATCHES;

	sprintf(token,"%c",watsonp ? this->genome : complCode[(int) this->genome]);
	tokens = push_token(tokens,token);

      } else if (this->comp == INDEL_COMP || this->comp == SHORTGAP_COMP) {
	if (this->genome == ' ') {
#if 0
	  /* Insertion relative to genome.  Ignored in MD string, but not in cigar string. */
	  if (state == IN_DELETION) {
	    tokens = push_token(tokens,"^");
	  }
	  nmatches++;
	  state = IN_MATCHES;
#endif

	} else if (this->cdna == ' ') {
	  /* Deletion relative to genome */
	  if (state == IN_MATCHES) {
	    if (nmatches > 0) {
	      sprintf(token,"%d",nmatches);
	      tokens = push_token(tokens,token);
	      nmatches = 0;
	    }

	  } else if (state == IN_MISMATCHES) {
	    tokens = push_token(tokens,"0");

	  }
	  state = IN_DELETION;

	  sprintf(token,"%c",watsonp ? this->genome : complCode[(int) this->genome]);
	  tokens = push_token(tokens,token);
	}

      } else {
	/* Ignore */
      }
    }

    /* Ignore terminal soft clipping */

    if (nmatches > 0) {
      sprintf(token,"%d",nmatches);
      tokens = push_token(tokens,token);
    }
    
    /* Keep tokens in reverse order */
  }


  /* Insert initial 0 token if necessary */
  if (tokens != NULL) {
    first_token = (char *) List_head(tokens);
    if (!isdigit(first_token[0])) {
      tokens = push_token(tokens,"0");
    }
  }
  
  return tokens;
}
#endif


static List_T
compute_md_string (int *nmismatches, int *nindels, struct T *pairs, int npairs, bool watsonp,
		   List_T cigar_tokens) {
  List_T md_tokens = NULL, p;
  char *cigar_token, token[10], *first_token, type;
  Pair_T this;
  int nmatches = 0, length;
  MD_state_T state = IN_MISMATCHES;
  int i, k = 0;

  *nmismatches = *nindels = 0;

  debug4(Pair_dump_array(pairs,npairs,true));
  debug4(printf("watsonp %d\n",watsonp));

  if (watsonp == true) {
    for (p = cigar_tokens; p != NULL; p = List_next(p)) {
      cigar_token = (char *) List_head(p);
      debug4(printf("token is %s\n",cigar_token));
      type = cigar_token[strlen(cigar_token)-1];
      length = atoi(cigar_token);
    
      if (type == 'H') {
	/* k += length; */

      } else if (type == 'S') {
	/* k += length; */

      } else if (type == 'M') {
	for (i = 0; i < length; i++, k++) {
	  this = &(pairs[k]);
	  debug4(printf("M %d/%d comp %c\n",i,length,this->comp));
	  if (this->comp == MATCH_COMP || this->comp == DYNPROG_MATCH_COMP || this->comp == AMBIGUOUS_COMP) {
	    nmatches++;
	    state = IN_MATCHES;

	  } else if (this->comp == MISMATCH_COMP) {
	    *nmismatches += 1;
	    if (state == IN_MATCHES) {
	      sprintf(token,"%d",nmatches);
	      md_tokens = push_token(md_tokens,token);
	      nmatches = 0;
	    } else if (state == IN_DELETION) {
	      md_tokens = push_token(md_tokens,"0");
	    }
	    state = IN_MISMATCHES;

	    sprintf(token,"%c",this->genome);
	    md_tokens = push_token(md_tokens,token);

	  } else {
	    fprintf(stderr,"Unexpected comp '%c'\n",this->comp);
	    exit(9);
	  }
	}

      } else if (type == 'I') {
	while (k < npairs && pairs[k].comp == INDEL_COMP && pairs[k].genome == ' ') {
	  *nindels += 1;
	  k++;
	}

      } else if (type == 'N') {
	while (k < npairs && pairs[k].gapp == true) {
	  k++;
	}

      } else if (type == 'D') {
	if (state == IN_MATCHES) {
	  if (nmatches > 0) {
	    sprintf(token,"%d",nmatches);
	    md_tokens = push_token(md_tokens,token);
	    nmatches = 0;
	  }
	}

	if (state != IN_DELETION) {
	  md_tokens = push_token(md_tokens,"^");
	}
	for (i = 0; i < length; i++, k++) {
	  this = &(pairs[k]);
	  sprintf(token,"%c",this->genome);
	  md_tokens = push_token(md_tokens,token);
	  *nindels += 1;
	}

	state = IN_DELETION;

      } else {
	fprintf(stderr,"Don't recognize type %c\n",type);
	abort();
      }
    }

    if (nmatches > 0) {
      sprintf(token,"%d",nmatches);
      md_tokens = push_token(md_tokens,token);
    }

    md_tokens = List_reverse(md_tokens);

  } else {
    cigar_tokens = List_reverse(cigar_tokens);
    for (p = cigar_tokens; p != NULL; p = List_next(p)) {
      cigar_token = (char *) List_head(p);
      debug4(printf("token is %s\n",cigar_token));
      type = cigar_token[strlen(cigar_token)-1];
      length = atoi(cigar_token);
    
      if (type == 'H') {
	/* k += length; */

      } else if (type == 'S') {
	/* k += length; */

      } else if (type == 'M') {
	if (state == IN_DELETION) {
	  md_tokens = push_token(md_tokens,"^");
	}

	for (i = 0; i < length; i++, k++) {
	  this = &(pairs[k]);
	  debug4(printf("M %d/%d comp %c\n",i,length,this->comp));
	  if (this->comp == MATCH_COMP || this->comp == DYNPROG_MATCH_COMP || this->comp == AMBIGUOUS_COMP) {
	    nmatches++;
	    state = IN_MATCHES;

	  } else if (this->comp == MISMATCH_COMP) {
	    *nmismatches += 1;
	    if (state == IN_MATCHES) {
	      sprintf(token,"%d",nmatches);
	      md_tokens = push_token(md_tokens,token);
	      nmatches = 0;
	    }
	    state = IN_MISMATCHES;

	    sprintf(token,"%c",complCode[(int) this->genome]);
	    md_tokens = push_token(md_tokens,token);

	  } else {
	    fprintf(stderr,"Unexpected comp '%c'\n",this->comp);
	    abort();
	  }
	}

      } else if (type == 'I') {
	while (k < npairs && pairs[k].comp == INDEL_COMP && pairs[k].genome == ' ') {
	  *nindels += 1;
	  k++;
	}

      } else if (type == 'N') {
	while (k < npairs && pairs[k].gapp == true) {
	  k++;
	}

      } else if (type == 'D') {
	if (state == IN_MATCHES) {
	  if (nmatches > 0) {
	    sprintf(token,"%d",nmatches);
	    md_tokens = push_token(md_tokens,token);
	    nmatches = 0;
	  }
	} else if (state == IN_MISMATCHES) {
	  md_tokens = push_token(md_tokens,"0");
	}

	for (i = 0; i < length; i++, k++) {
	  this = &(pairs[k]);
	  sprintf(token,"%c",complCode[(int) this->genome]);
	  md_tokens = push_token(md_tokens,token);
	  *nindels += 1;
	}
	state = IN_DELETION;

      } else {
	fprintf(stderr,"Don't recognize type %c\n",type);
	abort();
      }
    }

    if (nmatches > 0) {
      sprintf(token,"%d",nmatches);
      md_tokens = push_token(md_tokens,token);
    }

    /* Restore cigar_tokens */
    cigar_tokens = List_reverse(cigar_tokens);
  }

  assert(k == npairs);

  /* Insert initial 0 token if necessary */
  if (md_tokens != NULL) {
    first_token = (char *) List_head(md_tokens);
    if (!isdigit(first_token[0])) {
      md_tokens = push_token(md_tokens,"0");
    }
  }

  return md_tokens;
}




void
Pair_print_sam (FILE *fp, struct T *pairs, int npairs,
		char *accession, Chrnum_T chrnum, IIT_T chromosome_iit, Sequence_T usersegment,
		char *queryseq_ptr, char *quality_string,
		int hardclip5, int hardclip3, int querylength_given,
		bool watsonp, int cdna_direction, int chimera_part, Chimera_T chimera,
		int quality_shift, bool firstp, int pathnum, int npaths, int absmq_score, int second_absmq,
#ifdef GSNAP
		unsigned int flag, int pair_mapq_score, int end_mapq_score,
		Genomicpos_T chrpos, Genomicpos_T mate_chrpos, int pairedlength,
#else
		int mapq_score, bool sam_paired_p,
#endif
		char *sam_read_group_id) {

  char *chrstring = NULL;
#ifndef GSNAP
  unsigned int flag;
#endif

  List_T cigar_tokens = NULL, md_tokens = NULL;
  int nmismatches, nindels;
  bool intronp;
  int hardclip_low, hardclip_high;
  struct T *clipped_pairs;
  int clipped_npairs;


  if (chrnum == 0) {
    chrstring = Sequence_accession(usersegment);
  } else {
    chrstring = Chrnum_to_string(chrnum,chromosome_iit);
  }

#ifndef GSNAP
  flag = compute_sam_flag_nomate(pathnum,npaths,firstp,watsonp,sam_paired_p);
#endif


  /* Incoming hardclip5 and hardclip3 are due to overlaps, not chimera */
  if (firstp == true) {
    hardclip_low = 0;
    hardclip_high = hardclip5;
  } else {
    hardclip_low = hardclip3;
    hardclip_high = 0;
  }


  clipped_pairs = hardclip_pairs(&clipped_npairs,&hardclip_low,&hardclip_high,
				 pairs,npairs,querylength_given);

  /* Cigar updates hardclip5 and hardclip3 for chimeras */
  cigar_tokens = compute_cigar(&intronp,&hardclip_low,&hardclip_high,clipped_pairs,clipped_npairs,querylength_given,
			       watsonp,cdna_direction,chimera_part);

  md_tokens = compute_md_string(&nmismatches,&nindels,clipped_pairs,clipped_npairs,watsonp,cigar_tokens);

  print_sam_line(fp,firstp,accession,chrstring,
		 watsonp,cdna_direction,cigar_tokens,md_tokens,nmismatches,nindels,
		 intronp,queryseq_ptr,quality_string,hardclip_low,hardclip_high,
		 querylength_given,chimera,quality_shift,pathnum,npaths,
		 absmq_score,second_absmq,flag,
#ifdef GSNAP
		 pair_mapq_score,end_mapq_score,chrpos,mate_chrpos,pairedlength,
#else
		 mapq_score,clipped_pairs,clipped_npairs,
#endif
		 sam_read_group_id);

  /* Print procedures free the character strings */
  List_free(&md_tokens);
  List_free(&cigar_tokens);


  if (chrnum != 0) {
    FREE(chrstring);
  }

  return;
}


void
Pair_print_sam_nomapping (FILE *fp, char *accession, char *queryseq_ptr,
			  char *quality_string, int querylength, int quality_shift,
			  bool firstp, bool sam_paired_p, char *sam_read_group_id) {
  unsigned int flag;

#ifdef GSNAP
  fprintf(stderr,"Unexpected call to Pair_print_sam_nomapping in GSNAP\n");
  abort();
#endif

  /* 1. QNAME */
  fprintf(fp,"%s",accession);
  
  /* 2. FLAG */
  flag = compute_sam_flag_nomate(/*pathnum*/0,/*npaths*/0,firstp,/*watsonp*/true,sam_paired_p);
  fprintf(fp,"\t%u",flag);

  /* 3. RNAME: chr */
  fprintf(fp,"\t*");

  /* 4. POS: chrpos */
  fprintf(fp,"\t0");

  /* 5. MAPQ: Mapping quality */
  /* Picard says MAPQ should be 0 for an unmapped read */
  fprintf(fp,"\t0");

  /* 6. CIGAR */
  fprintf(fp,"\t*");

  /* 7. MRNM: Mate chr */
  /* 8. MPOS: Mate chrpos */
  /* 9. ISIZE: Insert size */
  fprintf(fp,"\t*\t0\t0\t");

  /* 10. SEQ: queryseq and 11. QUAL: quality scores */
  print_chopped(fp,queryseq_ptr,querylength,/*hardclip_low*/0,/*hardclip_high*/0);
  fprintf(fp,"\t");
  print_quality(fp,quality_string,querylength,/*hardclip_low*/0,/*hardclip_high*/0,
		quality_shift);

  /* 12. TAGS: RG */
  if (sam_read_group_id != NULL) {
    fprintf(fp,"\t");
    fprintf(fp,"RG:Z:%s",sam_read_group_id);
  }
  
  putc('\n',fp);

  return;
}


#endif



Uintlist_T
Pair_exonbounds (struct T *pairs, int npairs, Genomicpos_T chroffset) {
  Uintlist_T exonbounds = NULL;
  struct T *ptr, *prev, *this = NULL;
  bool in_exon = false;
  int i;
  Genomicpos_T last_genomepos = -1U;

  ptr = pairs;
  for (i = 0; i < npairs; i++) {
    prev = this;
    this = ptr++;

    if (this->gapp) {
      if (in_exon == true) {
	/* exon genomeend */
	exonbounds = Uintlist_push(exonbounds,chroffset + last_genomepos);
	in_exon = false;
      }
    } else if (this->comp == INTRONGAP_COMP) {
      /* Do nothing */
    } else {
      /* Remaining possibilities are MATCH_COMP, DYNPROG_MATCH_COMP, AMBIGUOUS_COMP, INDEL_COMP, 
	 SHORTGAP_COMP, or MISMATCH_COMP */
      if (in_exon == false) {
	/* exon genomestart */
	exonbounds = Uintlist_push(exonbounds,chroffset + this->genomepos);
	in_exon = true;
      }
    }
    if (this->genome != ' ') {
      last_genomepos = this->genomepos;
    }
  }

  prev = this;
  exonbounds = Uintlist_push(exonbounds,chroffset + last_genomepos);

  return Uintlist_reverse(exonbounds);
}


static int
count_psl_blocks_nt (Intlist_T *blockSizes, Intlist_T *qStarts, Uintlist_T *tStarts, struct T *pairs_directional,
		     int npairs, int querylength, bool watsonp) {
  int nblocks = 0, i;
  int block_querystart, block_queryend;
  struct T *ptr = pairs_directional, *prev, *this = NULL;
  bool in_block = false;
  int last_querypos = -1;
  Genomicpos_T last_genomepos = -1U;

  for (i = 0; i < npairs; i++) {
    prev = this;
    this = ptr++;

    if (this->gapp) {
      if (in_block == true) {
	nblocks++;
	block_queryend = last_querypos;
	debug2(fprintf(fp,"Block size: %d\n",abs(block_queryend-block_querystart)+1));
	*blockSizes = Intlist_push(*blockSizes,abs(block_queryend-block_querystart)+1);
	in_block = false;
      }
    } else if (this->comp == INTRONGAP_COMP) {
      /* Do nothing */

    } else if (this->cdna == ' ' || this->genome == ' ') {
      if (in_block == true) {
	nblocks++;
	block_queryend = last_querypos;
	debug2(fprintf(fp,"Block size: %d\n",abs(block_queryend-block_querystart)+1));
	*blockSizes = Intlist_push(*blockSizes,abs(block_queryend-block_querystart)+1);
	in_block = false;
      }

    } else {
      /* Remaining possibilities are MATCH_COMP, DYNPROG_MATCH_COMP, AMBIGUOUS_COMP, INDEL_COMP, 
         or SHORTGAP_COMP */
      if (in_block == false) {
	block_querystart = this->querypos;
	if (watsonp == true) {
	  debug2(fprintf(fp,"Pushing qstart: %d\n",block_querystart));
	  *qStarts = Intlist_push(*qStarts,block_querystart);
	} else {
	  debug2(fprintf(fp,"Pushing qstart: %d\n",querylength-block_querystart-1));
	  *qStarts = Intlist_push(*qStarts,querylength-block_querystart-1);
	}
	*tStarts = Uintlist_push(*tStarts,this->genomepos);
	in_block = true;
      }
    }

    if (this->cdna != ' ') {
      last_querypos = this->querypos;
    }
    if (this->genome != ' ') {
      last_genomepos = this->genomepos;
    }
  }

  if (in_block == true) {
    prev = this;
    nblocks++;
    block_queryend = last_querypos;
    debug2(fprintf(fp,"Block size: %d\n",abs(block_queryend-block_querystart)+1));
    *blockSizes = Intlist_push(*blockSizes,abs(block_queryend-block_querystart)+1);
  }

  *blockSizes = Intlist_reverse(*blockSizes);
  *qStarts = Intlist_reverse(*qStarts);
  *tStarts = Uintlist_reverse(*tStarts);

  return nblocks;
}


static int
count_psl_blocks_pro (Intlist_T *blockSizes, Intlist_T *qStarts, Uintlist_T *tStarts, struct T *pairs_directional,
		      int npairs, bool watsonp, Genomicpos_T chrlength) {
  int nblocks = 0, i;
  int naminoacids = 0;
  int block_querystart;
  struct T *ptr = pairs_directional, *this = NULL;
  bool in_block = false;
#ifdef NOGAPSINBLOCK
  struct T *prev;
#endif

  for (i = 0; i < npairs; i++) {
#ifdef NOGAPSINBLOCK
    prev = this;
#endif
    this = ptr++;

    if (this->gapp) {
      if (in_block == true) {
	nblocks++;
	*blockSizes = Intlist_push(*blockSizes,naminoacids);
	in_block = false;
	naminoacids = 0;
      }
    } else if (this->comp == INTRONGAP_COMP) {
      /* Do nothing */

#ifdef NOGAPSINBLOCK
    } else if (this->cdna == ' ' || this->genome == ' ') {
      if (in_block == true) {
	nblocks++;
	block_queryend = last_querypos;
	*blockSizes = Intlist_push(*blockSizes,block_queryend/3-(block_querystart+2)/3+1);
	in_block = false;
      }
#endif

    } else {
      /* Remaining possibilities are MATCH_COMP, DYNPROG_MATCH_COMP, AMBIGUOUS_COMP, INDEL_COMP, 
         or SHORTGAP_COMP */
      if (this->aa_e != ' ') {
	naminoacids++;
      }
      if (in_block == false) {
	block_querystart = this->querypos;
	*qStarts = Intlist_push(*qStarts,(block_querystart+2)/3);
	if (watsonp == true) {
	  *tStarts = Uintlist_push(*tStarts,this->genomepos);
	} else {
#if 0
	  /* Should be this */
	  *tStarts = Uintlist_push(*tStarts,this->genomepos);
#else
	  /* But is actually this */
	  *tStarts = Uintlist_push(*tStarts,chrlength - this->genomepos - 1);
#endif
	}
	in_block = true;
      }
    }
  }

  if (in_block == true) {
#ifdef NOGAPSINBLOCK
    prev = this;
#endif
    nblocks++;
    *blockSizes = Intlist_push(*blockSizes,naminoacids);
  }

  *blockSizes = Intlist_reverse(*blockSizes);
  *qStarts = Intlist_reverse(*qStarts);
  *tStarts = Uintlist_reverse(*tStarts);

  return nblocks;
}


static void
compute_gap_lengths_int (int *nbreaks, int *length, Intlist_T blockSizes, Intlist_T Starts, int nblocks) {
  int i;
  int start, end;
  /* Intlist_T p = blockSizes, q = Starts; */

  debug2(fprintf(fp,"Entered compute_gap_lengths_int with nblocks = %d, and Starts having length %d\n",
		nblocks,Intlist_length(Starts)));
  *nbreaks = *length = 0;
  for (i = 0; i < nblocks - 1; i++) {
    if (i > 0) {
      start = Intlist_head(Starts);
      if (start - end > 0) {
	*nbreaks += 1;
	*length += (start - end);
      }
      debug2(fprintf(fp,"%d - %d = %d, gap = %d\n",start,end,start-end,*length));
    }
    end = Intlist_head(Starts) + Intlist_head(blockSizes);
    blockSizes = Intlist_next(blockSizes);
    Starts = Intlist_next(Starts);
  }

  if (i > 0) {
    start = Intlist_head(Starts);
    if (start - end > 0) {
      *nbreaks += 1;
      *length += (start - end);
    }
    debug2(fprintf(fp,"%d - %d = %d, gap = %d\n",start,end,start-end,*length));
  }

  return;
}

static void
compute_gap_lengths_uint (int *nbreaks, int *length, Intlist_T blockSizes, Uintlist_T Starts, int nblocks) {
  int i;
  int start, end;
  /*
  Intlist_T p = blockSizes;
  Uintlist_T q = Starts;
  */

  *nbreaks = *length = 0;
  for (i = 0; i < nblocks - 1; i++) {
    if (i > 0) {
      start = Uintlist_head(Starts);
      if (start - end > 0) {
	*nbreaks += 1;
	*length += (start - end);
      }
      debug2(fprintf(fp,"%d - %d = %d, gap = %d\n",start,end,start-end,*length));
    }
    end = Uintlist_head(Starts) + Intlist_head(blockSizes);
    blockSizes = Intlist_next(blockSizes);
    Starts = Uintlist_next(Starts);
  }

  if (i > 0) {
    start = Uintlist_head(Starts);
    if (start - end > 0) {
      *nbreaks += 1;
      *length += (start - end);
    }
    debug2(fprintf(fp,"%d - %d = %d, gap = %d\n",start,end,start-end,*length));
  }

  return;
}



static void
count_matches_pro (int *matches, int *mismatches, int *unknowns, 
		   struct T *pairs, int npairs) {
  struct T *this = NULL;
  int i;

  i = 0;
  while (i < npairs) {
    /* prev = this; */
    this = &(pairs[i++]);

    if (this->gapp == false) {
      if (this->aa_g != ' ' && this->aa_e != ' ') {
	if (this->aa_g == this->aa_e) {
	  *matches += 1;
	} else if (this->aa_e == 'X') {
	  *unknowns += 1;
	} else {
	  *mismatches += 1;
	}
      }
    }
  }	
  return;
}



void
Pair_print_pslformat_nt (FILE *fp, struct T *pairs, int npairs, T start, T end,
			 Sequence_T queryseq, Chrnum_T chrnum,
			 IIT_T chromosome_iit, Sequence_T usersegment,
			 int matches, int unknowns, int mismatches, 
			 bool watsonp) {
  Genomicpos_T chrpos1, chrpos2;
  struct T *pairs_directional = NULL;
  Intlist_T blockSizes = NULL, qStarts = NULL, p;
  Uintlist_T tStarts = NULL, q;
  int nblocks;
  int qnbreaks, qlength, tnbreaks, tlength, querylength;
  char *chr;

#ifdef PMAP
    querylength = 3*Sequence_fulllength(queryseq);
#else
    querylength = Sequence_fulllength(queryseq);
#endif

  if (watsonp == true) {
    pairs_directional = pairs;
  } else {
    pairs_directional = invert_and_revcomp_path(pairs,npairs);
  }

  nblocks = count_psl_blocks_nt(&blockSizes,&qStarts,&tStarts,pairs_directional,npairs,
				querylength,watsonp);
  compute_gap_lengths_int(&qnbreaks,&qlength,blockSizes,qStarts,nblocks);
  compute_gap_lengths_uint(&tnbreaks,&tlength,blockSizes,tStarts,nblocks);

  fprintf(fp,"%d\t%d\t%d\t%d\t",matches,mismatches,/*repeatmatches*/0,unknowns);
  fprintf(fp,"%d\t%d\t%d\t%d\t",qnbreaks,qlength,tnbreaks,tlength);

  if (watsonp == true) {
    fprintf(fp,"+");
  } else {
    fprintf(fp,"-");
  }
  fprintf(fp,"\t%s\t%d",Sequence_accession(queryseq),Sequence_fulllength_given(queryseq));

  fprintf(fp,"\t%d\t%d",start->querypos,end->querypos+1);

  /* T name and T size */
  if (chrnum == 0) {
    fprintf(fp,"\t%s\t%u",Sequence_accession(usersegment),Sequence_fulllength(usersegment));
  } else {
    chr = Chrnum_to_string(chrnum,chromosome_iit);
    fprintf(fp,"\t%s\t%u",chr,Chrnum_length(chrnum,chromosome_iit));
    FREE(chr);
  }

  /* T start and T end */
  chrpos1 = start->genomepos;
  chrpos2 = end->genomepos;
  if (watsonp) {
    fprintf(fp,"\t%u\t%u",chrpos1,chrpos2+1U);
  } else {
    fprintf(fp,"\t%u\t%u",chrpos2,chrpos1+1U);
  }

  fprintf(fp,"\t%d",nblocks);

  fprintf(fp,"\t");
  for (p = blockSizes; p != NULL; p = Intlist_next(p)) {
    fprintf(fp,"%d,",Intlist_head(p));
  }

  fprintf(fp,"\t");
  for (p = qStarts; p != NULL; p = Intlist_next(p)) {
    fprintf(fp,"%d,",Intlist_head(p));
  }

  fprintf(fp,"\t");
  for (q = tStarts; q != NULL; q = Uintlist_next(q)) {
    fprintf(fp,"%u,",Uintlist_head(q));
  }

  Intlist_free(&blockSizes);
  Intlist_free(&qStarts);
  Uintlist_free(&tStarts);

  if (watsonp == false) {
    FREE(pairs_directional);
  }

  putc('\n',fp);
  return;
}

void
Pair_print_pslformat_pro (FILE *fp, struct T *pairs, int npairs, T start, T end,
			  Sequence_T queryseq, Chrnum_T chrnum,
			  IIT_T chromosome_iit, Sequence_T usersegment,
			  bool watsonp, int cdna_direction) {
  Genomicpos_T chrpos1, chrpos2, chrlength;
  Intlist_T blockSizes = NULL, qStarts = NULL, p;
  Uintlist_T tStarts = NULL, q;
  int nblocks, matches = 0, mismatches = 0, unknowns = 0;
  int qnbreaks, qlength, tnbreaks, tlength;
  char *chr;

  chrlength = Chrnum_length(chrnum,chromosome_iit);
  nblocks = count_psl_blocks_pro(&blockSizes,&qStarts,&tStarts,pairs,npairs,
				 watsonp,chrlength);
  compute_gap_lengths_int(&qnbreaks,&qlength,blockSizes,qStarts,nblocks);
  compute_gap_lengths_uint(&tnbreaks,&tlength,blockSizes,tStarts,nblocks);

  count_matches_pro(&matches,&mismatches,&unknowns,pairs,npairs);

  fprintf(fp,"%d\t%d\t%d\t%d\t",matches,mismatches,/*repeatmatches*/0,unknowns);
  fprintf(fp,"%d\t%d\t%d\t%d\t",qnbreaks,qlength,tnbreaks,tlength);

  if (cdna_direction >= 0) {
    fprintf(fp,"+");
  } else {
    fprintf(fp,"-");
  }

  if (watsonp == true) {
    fprintf(fp,"+");
  } else {
    fprintf(fp,"-");
  }
  fprintf(fp,"\t%s\t%d",Sequence_accession(queryseq),Sequence_fulllength_given(queryseq));

  fprintf(fp,"\t%d\t%d",(start->querypos+2)/3,end->querypos/3+1);

  /* T name and T size */
  if (chrnum == 0) {
    fprintf(fp,"\t%s\t%u",Sequence_accession(usersegment),Sequence_fulllength(usersegment));
  } else {
    chr = Chrnum_to_string(chrnum,chromosome_iit);
    fprintf(fp,"\tchr%s\t%u",chr,Chrnum_length(chrnum,chromosome_iit));
    FREE(chr);
  }

  /* T start and T end */
  chrpos1 = start->genomepos;
  chrpos2 = end->genomepos;
  if (watsonp) {
    fprintf(fp,"\t%u\t%u",chrpos1,chrpos2+1U);
  } else {
    fprintf(fp,"\t%u\t%u",chrpos2,chrpos1+1U);
  }

  nblocks = count_psl_blocks_pro(&blockSizes,&qStarts,&tStarts,pairs,npairs,
				 watsonp,chrlength);
  fprintf(fp,"\t%d",nblocks);
  fprintf(fp,"\t");

  for (p = blockSizes; p != NULL; p = Intlist_next(p)) {
    fprintf(fp,"%d,",Intlist_head(p));
  }

  fprintf(fp,"\t");
  for (p = qStarts; p != NULL; p = Intlist_next(p)) {
    fprintf(fp,"%d,",Intlist_head(p));
  }

  fprintf(fp,"\t");

  for (q = tStarts; q != NULL; q = Uintlist_next(q)) {
    fprintf(fp,"%u,",Uintlist_head(q));
  }

  Intlist_free(&blockSizes);
  Intlist_free(&qStarts);
  Uintlist_free(&tStarts);

  putc('\n',fp);
  return;
}

void
Pair_print_exons (FILE *fp, struct T *pairs, int npairs, int wraplength, int ngap, bool cdnap) {
  bool in_exon = false;
  struct T *ptr, *this = NULL;
  int i, exonno = 0, column = 0;

  ptr = pairs;
  for (i = 0; i < npairs; i++) {
    this = ptr++;

    if (this->gapp) {
      if (in_exon == true) {
	if (column != 0) {
	  putc('\n',fp);
	  column = 0;
	}
	fprintf(fp,"</exon>\n");
	in_exon = false;
	if (ngap > 0) {
	  fprintf(fp,"<intron %d>\n",exonno);
	  putc(this->genome,fp);
	  column = 1;
	}
      } else {
	if (ngap > 0) {
	  putc(this->genome,fp);
	  if (++column % wraplength == 0) {
	    putc('\n',fp);
	    column = 0;
	  }
	}
      }
    } else if (this->comp == INTRONGAP_COMP) {
      /* Do nothing */
    } else {
      /* Remaining possibilities are MATCH_COMP, DYNPROG_MATCH_COMP, AMBIGUOUS_COMP, INDEL_COMP,
	 SHORTGAP_COMP, or MISMATCH_COMP */
      if (in_exon == false) {
	if (ngap > 0) {
	  if (exonno > 0) {
	    if (column != 0) {
	      putc('\n',fp);
	      column = 0;
	    }
	    fprintf(fp,"</intron>\n");
	  }
	}
	fprintf(fp,"<exon %d",++exonno);
	if (cdnap == true) {
	  if (this->aaphase_e >= 0) {
	    fprintf(fp,", phase %d",this->aaphase_e);
	  }
	} else {
	  if (this->aaphase_g >= 0) {
	    fprintf(fp,", phase %d",this->aaphase_g);
	  }
	}
	fprintf(fp,">\n");
	in_exon = true;
      }
      if (cdnap == true) {
	if (this->cdna != ' ') {
	  putc(this->cdna,fp);
	  if (++column % wraplength == 0) {
	    putc('\n',fp);
	    column = 0;
	  }
	}
      } else {
	if (this->genome != ' ') {
	  putc(this->genome,fp);
	  if (++column % wraplength == 0) {
	    putc('\n',fp);
	    column = 0;
	  }
	}
      }
    }
  }
  if (column != 0) {
    putc('\n',fp);
  }
  fprintf(fp,"</exon>\n");

  return;
}


int
Pair_nmatches (List_T pairs) {
  int nmatches = 0;
  bool in_intron = false, indelp = false;
  List_T p;
  T this;

  for (p = pairs; p != NULL; p = p->rest) {
    this = p->first;
    if (this->gapp) {
      if (this->donor_prob < 0.90 && this->acceptor_prob < 0.90) {
	nmatches -= 2;
      }
      if (!in_intron) {
	in_intron = true;
      }
    } else {
      if (in_intron) {
	in_intron = false;
      }
      if (this->comp == INDEL_COMP || this->comp == SHORTGAP_COMP) {
	indelp = true;
#ifndef PMAP
      } else if (unknown_base(this->cdna) || unknown_base(this->genome) || this->comp == AMBIGUOUS_COMP) {
	/* (*unknowns)++; */
#endif
      } else if (this->comp == MATCH_COMP || this->comp == DYNPROG_MATCH_COMP || this->comp == AMBIGUOUS_COMP) {
	nmatches++;
      } else if (this->comp == MISMATCH_COMP) {
	/* (*mismatches)++; */
      } else {
	fprintf(stderr,"Can't parse comp %c, gapp %d\n",this->comp,this->gapp);
	abort();
      }
    }
  }

  if (indelp == true) {
    nmatches -= 2;
  }

  return nmatches;
}


void
Pair_fracidentity_simple (int *matches, int *unknowns, int *mismatches, List_T pairs) {
  bool in_intron = false;
  List_T p;
  T this;

  *matches = *unknowns = *mismatches = 0;
  for (p = pairs; p != NULL; p = p->rest) {
    this = p->first;
    if (this->gapp) {
      if (!in_intron) {
	in_intron = true;
      }
    } else {
      if (in_intron) {
	in_intron = false;
      }
      if (this->comp == INDEL_COMP || this->comp == SHORTGAP_COMP) {
#ifndef PMAP
      } else if (unknown_base(this->cdna) || unknown_base(this->genome) || this->comp == AMBIGUOUS_COMP) {
	(*unknowns)++;
#endif
      } else if (this->comp == MATCH_COMP || this->comp == DYNPROG_MATCH_COMP || this->comp == AMBIGUOUS_COMP) {
	(*matches)++;
      } else if (this->comp == MISMATCH_COMP) {
	(*mismatches)++;
      } else {
	fprintf(stderr,"Can't parse comp %c, gapp %d\n",this->comp,this->gapp);
	abort();
      }
    }
  }

  return;
}


void
Pair_fracidentity (int *matches, int *unknowns, int *mismatches, int *qopens, int *qindels, 
		   int *topens, int *tindels, int *ncanonical, int *nsemicanonical, int *nnoncanonical,
		   List_T pairs, int cdna_direction) {
  bool in_intron = false;
  List_T p;
  T this, prev = NULL;

  *matches = *unknowns = *mismatches = *qopens = *qindels = *topens = *tindels = 
    *ncanonical = *nsemicanonical = *nnoncanonical = 0;

  for (p = pairs; p != NULL; p = p->rest) {
    this = p->first;
    if (this->gapp) {
      if (!in_intron) {
	if (cdna_direction > 0) {
	  if (this->comp == FWD_CANONICAL_INTRON_COMP) {
	    (*ncanonical)++;
	  } else if (this->comp == FWD_GCAG_INTRON_COMP || this->comp == FWD_ATAC_INTRON_COMP) {
	    (*nsemicanonical)++;
	  } else if (this->comp == NONINTRON_COMP) {
	    (*nnoncanonical)++;
	  }

	} else if (cdna_direction < 0) {
	  if (this->comp == REV_CANONICAL_INTRON_COMP) {
	    (*ncanonical)++;
	  } else if (this->comp == REV_GCAG_INTRON_COMP || this->comp == REV_ATAC_INTRON_COMP) {
	    (*nsemicanonical)++;
	  } else if (this->comp == NONINTRON_COMP) {
	    (*nnoncanonical)++;
	  }
	}
	in_intron = true;
      }
    } else {
      if (in_intron) {
	in_intron = false;
      }
      if (this->comp == INDEL_COMP || this->comp == SHORTGAP_COMP) {
	if (this->cdna == ' ') {
	  (*tindels)++;		/* If genome has extra char, count it as a genome skip */
	  if (prev && prev->cdna != ' ') {
	    (*topens)++;
	  }
	} else if (this->genome == ' ') {
	  (*qindels)++;
	  if (prev && prev->genome != ' ') {
	    (*qopens)++;
	  }
	} else {
	  fprintf(stderr,"Can't parse comp %c, cdna %c, genome %c\n",
		  this->comp,this->cdna,this->genome);
	  abort();
	}
#ifndef PMAP
      } else if (unknown_base(this->cdna) || unknown_base(this->genome) || this->comp == AMBIGUOUS_COMP) {
	(*unknowns)++;
#endif
      } else if (this->comp == MATCH_COMP || this->comp == DYNPROG_MATCH_COMP || this->comp == AMBIGUOUS_COMP) {
	(*matches)++;
      } else if (this->comp == MISMATCH_COMP) {
	(*mismatches)++;
      } else {
	fprintf(stderr,"Can't parse comp %c, gapp %d\n",this->comp,this->gapp);
	abort();
      }
    }
    prev = this;
  }

  return;
}


#if 0
/* Called on first and last exons during distal/medial calculation */
/* Procedure seems to give random results */
int
Pair_fracidentity_changepoint (List_T pairs, int cdna_direction) {
  int changepoint = 0, maxscore = 0, score = 0;
  int i = 0;

  bool in_intron = false;
  List_T p;
  T this, prev = NULL;

  for (p = pairs; p != NULL; p = p->rest) {
    i++;
    this = p->first;
    debug3(fprintf(fp,"%d: ",i));
    debug3(Pair_dump_one(this,/*zerobasedp*/false));
    if (this->gapp) {
      if (!in_intron) {
#if 0
	/* Don't expect an intron */
	if (cdna_direction > 0) {
	  if (this->comp == FWD_CANONICAL_INTRON_COMP) {
	    (*ncanonical)++;
	  } else if (this->comp == FWD_GCAG_INTRON_COMP || this->comp == FWD_ATAC_INTRON_COMP) {
	    (*nsemicanonical)++;
	  } else if (this->comp == NONINTRON_COMP) {
	    (*nnoncanonical)++;
	  }

	} else if (cdna_direction < 0) {
	  if (this->comp == REV_CANONICAL_INTRON_COMP) {
	    (*ncanonical)++;
	  } else if (this->comp == REV_GCAG_INTRON_COMP || this->comp == REV_ATAC_INTRON_COMP) {
	    (*nsemicanonical)++;
	  } else if (this->comp == NONINTRON_COMP) {
	    (*nnoncanonical)++;
	  }
	}
#endif
	in_intron = true;
      }
    } else {
      if (in_intron) {
	in_intron = false;
      }
      if (this->comp == INDEL_COMP || this->comp == SHORTGAP_COMP) {
	if (this->cdna == ' ') {
	  score += TINDEL;
	  if (prev && prev->cdna != ' ') {
	    score += TOPEN;
	  }
	} else if (this->genome == ' ') {
	  score += QINDEL;
	  if (prev && prev->genome != ' ') {
	    score += QOPEN;
	  }
	} else {
	  fprintf(stderr,"Can't parse comp %c, cdna %c, genome %c\n",
		  this->comp,this->cdna,this->genome);
	  abort();
	}
#ifndef PMAP
      } else if (unknown_base(this->cdna) || unknown_base(this->genome) || this->comp == AMBIGUOUS_COMP) {
	/* (*unknowns)++; */
#endif
      } else if (this->comp == MATCH_COMP || this->comp == DYNPROG_MATCH_COMP || this->comp == AMBIGUOUS_COMP) {
#if 0
	score += (MATCH + MATCH); /* Give more weight to matches to allow for poor quality at ends */
#else
	score += MATCH;
#endif
	if (score > maxscore) {
	  maxscore = score;
	  changepoint = i;
	  debug3(fprintf(fp," => maxscore %d",maxscore));
	}
      } else if (this->comp == MISMATCH_COMP) {
	score += MISMATCH;
      } else {
	fprintf(stderr,"Can't parse comp %c, gapp %d\n",this->comp,this->gapp);
	abort();
      }
    }
    debug3(fprintf(fp,"\n"));
    prev = this;
  }

  return changepoint;
}
#endif


int
Pair_fracidentity_score (List_T pairs, int cdna_direction) {
  int score = 0;
  int i = 0;

  bool in_intron = false;
  List_T p;
  T this, prev = NULL;

  for (p = pairs; p != NULL; p = p->rest) {
    i++;
    this = p->first;
    debug3(fprintf(fp,"%d: ",i));
    debug3(Pair_dump_one(this,/*zerobasedp*/false));
    if (this->gapp) {
      if (!in_intron) {
	in_intron = true;
      }
    } else {
      if (in_intron) {
	in_intron = false;
      }
      if (this->comp == INDEL_COMP || this->comp == SHORTGAP_COMP) {
	if (this->cdna == ' ') {
	  score += TINDEL;
	  if (prev && prev->cdna != ' ') {
	    score += TOPEN;
	  }
	} else if (this->genome == ' ') {
	  score += QINDEL;
	  if (prev && prev->genome != ' ') {
	    score += QOPEN;
	  }
	} else {
	  fprintf(stderr,"Can't parse comp %c, cdna %c, genome %c\n",
		  this->comp,this->cdna,this->genome);
	  abort();
	}
#ifndef PMAP
      } else if (unknown_base(this->cdna) || unknown_base(this->genome) || this->comp == AMBIGUOUS_COMP) {
	/* (*unknowns)++; */
#endif
      } else if (this->comp == MATCH_COMP || this->comp == DYNPROG_MATCH_COMP || this->comp == AMBIGUOUS_COMP) {
	score += MATCH;
      } else if (this->comp == MISMATCH_COMP) {
	score += MISMATCH;
      } else {
	fprintf(stderr,"Can't parse comp %c, gapp %d\n",this->comp,this->gapp);
	abort();
      }
    }
    debug3(fprintf(fp,"\n"));
    prev = this;
  }

  return score;
}


double
Pair_frac_error (List_T pairs, int cdna_direction) {
  int matches, unknowns, mismatches, qopens, qindels,
    topens, tindels, ncanonical, nsemicanonical, nnoncanonical;
  int den;

  Pair_fracidentity(&matches,&unknowns,&mismatches,&qopens,&qindels, 
		    &topens,&tindels,&ncanonical,&nsemicanonical,&nnoncanonical,
		    pairs,cdna_direction);

  if ((den = matches + mismatches + qindels + tindels) == 0) {
    return 1.0;
  } else {
    return (double) (mismatches + qindels + tindels)/(double) den;
  }
}

void
Pair_fracidentity_bounded (int *matches, int *unknowns, int *mismatches, 
			   int *qopens, int *qindels, int *topens, int *tindels,
			   int *ncanonical, int *nsemicanonical, int *nnoncanonical,
			   struct T *ptr, int npairs, 
			   int cdna_direction, int minpos, int maxpos) {
  bool in_intron = false;
  T this, prev = NULL;
  int i;

  *matches = *unknowns = *mismatches = *qopens = *qindels = *topens = *tindels = 
    *ncanonical = *nsemicanonical = *nnoncanonical = 0;
  
  for (i = 0; i < npairs; i++) {
    this = ptr++;
    if (this->gapp) {
      if (!in_intron) {
	if (this->querypos >= minpos && this->querypos <= maxpos) {
	  if (this->comp == FWD_CANONICAL_INTRON_COMP) {
	    (*ncanonical)++;
	  } else if (this->comp == FWD_GCAG_INTRON_COMP || this->comp == FWD_ATAC_INTRON_COMP) {
	    (*nsemicanonical)++;
	  } else if (this->comp == NONINTRON_COMP) {
	    (*nnoncanonical)++;
	  }
	} else if (cdna_direction < 0) {
	  if (this->comp == REV_CANONICAL_INTRON_COMP) {
	    (*ncanonical)++;
	  } else if (this->comp == REV_GCAG_INTRON_COMP || this->comp == REV_ATAC_INTRON_COMP) {
	    (*nsemicanonical)++;
	  } else if (this->comp == NONINTRON_COMP) {
	    (*nnoncanonical)++;
	  }
	}
	in_intron = true;
      }
    } else {
      if (in_intron) {
	in_intron = false;
      }
      if (this->querypos >= minpos && this->querypos <= maxpos) {
	if (this->comp == INDEL_COMP || this->comp == SHORTGAP_COMP) {
	  if (this->cdna == ' ') {
	    (*tindels)++;		/* If genome has extra char, count it as a genome skip */
	    if (prev && prev->cdna != ' ') {
	      (*topens)++;
	    }
	  } else if (this->genome == ' ') {
	    (*qindels)++;
	    if (prev && prev->genome != ' ') {
	      (*qopens)++;
	    }
	  } else {
	    fprintf(stderr,"Can't parse comp %c, cdna %c, genome %c\n",
		    this->comp,this->cdna,this->genome);
	    abort();
	  }
#ifndef PMAP
	} else if (unknown_base(this->cdna) || unknown_base(this->genome) || this->comp == AMBIGUOUS_COMP) {
	  (*unknowns)++;
#endif
	} else if (this->comp == MATCH_COMP || this->comp == DYNPROG_MATCH_COMP || this->comp == AMBIGUOUS_COMP) {
	  (*matches)++;
	} else if (this->comp == MISMATCH_COMP) {
	  (*mismatches)++;
	} else {
	  fprintf(stderr,"Can't parse comp %c, gapp %d\n",this->comp,this->gapp);
	  abort();
	}
      }
    }
    prev = this;
  }
  return;
}

static const Except_T Array_bounds_error = { "Exceeded array bounds" };


int *
Pair_matchscores (struct T *ptr, int npairs, int querylength) {
  int *matchscores;
  T this;
  int querypos;
  int i;

  matchscores = (int *) CALLOC(querylength,sizeof(int));

  for (i = 0; i < npairs; i++) {
    this = ptr++;
    querypos = this->querypos;

    if (this->gapp) {
      matchscores[querypos] = 0;	/* Count as mismatch; make evidence support the gap */
    } else if (this->comp == MISMATCH_COMP) {
      matchscores[querypos] = 0; /* For mismatch */
    } else if (this->comp == INDEL_COMP) {
      matchscores[querypos] = -1;	/* Ignore indels */
    } else {
      matchscores[querypos] = 1; /* For match */
    }
  }

  return matchscores;
}


int *
Pair_matchscores_list (int *nmatches, int *ntotal, int *length, List_T pairs) {
  int *matchscores;
  T this;
  List_T p;
  int i = 0;

  matchscores = (int *) CALLOC(List_length(pairs),sizeof(int));
  *nmatches = *ntotal = *length = 0;

  for (p = pairs; p != NULL; p = p->rest) {
    this = p->first;
    if (this->gapp) {
      matchscores[i++] = 0;	/* Count as mismatch; make evidence support the gap */
      (*ntotal) += 1;
    } else if (this->comp == MISMATCH_COMP) {
      matchscores[i++] = 0; /* For mismatch */
      (*ntotal) += 1;
#ifndef PMAP
    } else if (this->comp == AMBIGUOUS_COMP) {
      matchscores[i++] = 0; /* For cases involving 'N' */
      (*ntotal) += 1;
#endif
    } else if (this->comp == INDEL_COMP) {
      matchscores[i++] = -1;	/* Ignore indels */
    } else {
      matchscores[i++] = 1; /* For match */
      (*nmatches) += 1;
      (*ntotal) += 1;
    }
    (*length) += 1;
  }

  return matchscores;
}


void
Pair_pathscores (bool *gapp, int *pathscores, struct T *ptr, int npairs, 
		 int cdna_direction, int querylength, cDNAEnd_T cdnaend) {
  int querypos, querystart, queryend;
  bool in_intron = false;
  T this, prev = NULL;
  int i;

  /* Determine these before ptr changes */
  this = &(ptr[0]);
  querystart = this->querypos;
  this = &(ptr[npairs-1]);
  queryend = this->querypos;

  for (i = 0; i < npairs; i++) {
    this = ptr++;

    querypos = this->querypos;
    if (querypos >= querylength) {
      fprintf(stderr,"Pair_pathscores: querypos %d >= querylength %d\n",querypos,querylength);
      Pair_dump_array(ptr,npairs,/*zerobasedp*/true);
      fflush(stdout);
      abort();
      RAISE(Array_bounds_error);
    }

    if (this->gapp) {
      gapp[querypos] = true;
      if (in_intron == false) {
	/* Adds only a single reward/penalty per intron */
	if (cdna_direction > 0) {
	  if (this->comp == FWD_CANONICAL_INTRON_COMP) {
	    pathscores[querypos] = CANONICAL_POINTS;
	  } else if (this->comp == FWD_GCAG_INTRON_COMP || this->comp == FWD_ATAC_INTRON_COMP) {
	    pathscores[querypos] = SEMICANONICAL_POINTS;
	  } else {
	    pathscores[querypos] = NONCANONICAL_POINTS; /* noncanonical */
	  }
	} else if (cdna_direction < 0) {
	  if (this->comp == REV_CANONICAL_INTRON_COMP) {
	    pathscores[querypos] = CANONICAL_POINTS;
	  } else if (this->comp == REV_GCAG_INTRON_COMP || this->comp == REV_ATAC_INTRON_COMP) {
	    pathscores[querypos] = SEMICANONICAL_POINTS;
	  } else {
	    pathscores[querypos] = NONCANONICAL_POINTS; /* noncanonical */
	  }
	} else {
	  pathscores[querypos] = NONCANONICAL_POINTS; /* indeterminate */
	}
	in_intron = true;
      }

    } else {
      if (in_intron) {
	in_intron = false;
      }
      if (this->comp == INDEL_COMP || this->comp == SHORTGAP_COMP) {
	if (this->cdna == ' ') {
	  pathscores[querypos] = TINDEL;
	  if (prev && prev->cdna != ' ') {
	    pathscores[querypos] = TOPEN;
	  }
	} else if (this->genome == ' ') {
	  pathscores[querypos] = QINDEL;
	  if (prev && prev->genome != ' ') {
	    pathscores[querypos] = QOPEN;
	  }
	} else {
	  fprintf(stderr,"Can't parse comp %c, cdna %c, genome %c\n",
		  this->comp,this->cdna,this->genome);
	  abort();
	}
#ifndef PMAP
      } else if (unknown_base(this->cdna) || unknown_base(this->genome) || this->comp == AMBIGUOUS_COMP) {
	/* (*unknowns)++; */
#endif
      } else if (this->comp == MATCH_COMP || this->comp == DYNPROG_MATCH_COMP || this->comp == AMBIGUOUS_COMP) {
	pathscores[querypos] = +1; /* For match */
      } else if (this->comp == MISMATCH_COMP) {
	pathscores[querypos] = MISMATCH;
      } else {
	fprintf(stderr,"Can't parse comp %c, gapp %d\n",this->comp,this->gapp);
	abort();
      }
    }
    prev = this;
  }

  /* Gets querystart to queryend inclusive */
  if (querystart == 0) {
    for (i = 1; i <= queryend; i++) {
      pathscores[i] += pathscores[i-1];
    }
  } else {
    for (i = querystart; i <= queryend; i++) {
      pathscores[i] += pathscores[i-1];
    }
  }

  if (cdnaend == FIVE) {
    for (i = queryend + 1; i < querylength; i++) {
      pathscores[i] = pathscores[i-1] + QINDEL;
    }
  } else if (cdnaend == THREE) {
    for (i = querystart - 1; i >= 0; --i) {
      pathscores[i] = pathscores[i+1] - QINDEL;
    }
    for (i = queryend + 1; i < querylength; i++) {
      pathscores[i] = pathscores[i-1];
    }
  }

  return;
}


int
Pair_nexons_approx (List_T pairs) {
  int nexons = 0;
  bool in_exon = false;
  T this;
  List_T p;
  
  for (p = pairs; p != NULL; p = List_next(p)) {
    this = List_head(p);
    if (this->gapp) {
      if (in_exon) {
	in_exon = false;
      }
    } else {
      if (!in_exon) {
	nexons++;
	in_exon = true;
      }
    }
  }

  return nexons;
}


int
Pair_nexons (struct T *pairs, int npairs) {
  int nexons = 0;
  struct T *ptr, *this = NULL;
  bool in_exon = false;
  int i;
  
  ptr = pairs;
  for (i = 0; i < npairs; i++) {
    this = ptr++;
    if (this->gapp) {
      if (in_exon) {
	in_exon = false;
      }
    } else if (this->comp == INTRONGAP_COMP) {
      /* Do nothing */
    } else {
      if (!in_exon) {
	nexons++;
	in_exon = true;
      }
    }
  }

  return nexons;
}

int
Pair_cdna_direction (List_T pairs) {
  int cdna_direction = 0;
  bool in_intron = false;
  T this;
  List_T p;
  
  for (p = pairs; p != NULL; p = List_next(p)) {
    this = List_head(p);
    if (this->gapp) {
      if (!in_intron) {
	if (this->comp == FWD_CANONICAL_INTRON_COMP) {
	  cdna_direction += 1;
	} else if (this->comp == REV_CANONICAL_INTRON_COMP) {
	  cdna_direction -= 1;
	}
	in_intron = true;
      }
    } else {
      if (in_intron) {
	in_intron = false;
      }
    }
  }

  return cdna_direction;
}

bool
Pair_consistentp (int *ncanonical, struct T *pairs, int npairs, int cdna_direction) {
  bool in_intron = false;
  struct T *this;
  int i;

  *ncanonical = 0;
  for (i = 0; i < npairs; i++) {
    this = pairs++;
    if (this->gapp) {
      if (!in_intron) {
	if (cdna_direction > 0) {
	  if (this->comp == REV_CANONICAL_INTRON_COMP || 
	      this->comp == REV_GCAG_INTRON_COMP || 
	      this->comp == REV_ATAC_INTRON_COMP) {
	    return false;
	  } else if (this->comp == FWD_CANONICAL_INTRON_COMP) {
	    (*ncanonical)++;
	  }
	} else if (cdna_direction < 0) {
	  if (this->comp == FWD_CANONICAL_INTRON_COMP || 
	      this->comp == FWD_GCAG_INTRON_COMP || 
	      this->comp == FWD_ATAC_INTRON_COMP) {
	    return false;
	  } else if (this->comp == REV_CANONICAL_INTRON_COMP) {
	    (*ncanonical)++;
	  }
	} else if (cdna_direction == 0) {
	  /* Set cdna_direction for next time */
	  if (this->comp == FWD_CANONICAL_INTRON_COMP || 
	      this->comp == FWD_GCAG_INTRON_COMP || 
	      this->comp == FWD_ATAC_INTRON_COMP) {
	    cdna_direction = +1;
	  } else if (this->comp == REV_CANONICAL_INTRON_COMP || 
		     this->comp == REV_GCAG_INTRON_COMP || 
		     this->comp == REV_ATAC_INTRON_COMP) {
	    cdna_direction = -1;
	  } 
	}
	in_intron = true;
      }
    } else {
      if (in_intron) {
	in_intron = false;
      }
    }
  }

  return true;
}


static void
invert_intron (char *donor, char *acceptor) {
  char temp;

  temp = donor[0];
  donor[0] = complCode[(int) acceptor[1]];
  acceptor[1] = complCode[(int) temp];

  temp = donor[1];
  donor[1] = complCode[(int) acceptor[0]];
  acceptor[0] = complCode[(int) temp];
  
  return;
}


void
Pair_print_protein_genomic (FILE *fp, struct T *ptr, int npairs, int wraplength, bool forwardp) {
  struct T *this;
  int xpos = 0, i;

  if (forwardp == true) {
    for (i = 0; i < npairs; i++) {
      this = ptr++;
      if (this->aa_g != ' ') {
	if (xpos == wraplength) {
	  putc('\n',fp);
	  xpos = 0;
	}
#ifdef PMAP
	putc(this->aa_g,fp);
	xpos++;
#else
	if (this->aa_g != '*') {
	  putc(this->aa_g,fp);
	  xpos++;
	}
#endif
      }
    }
    putc('\n',fp);

  } else {
    for (i = npairs-1; i >= 0; i--) {
      this = ptr--;
      if (this->aa_g != ' ') {
	if (xpos == wraplength) {
	  putc('\n',fp);
	  xpos = 0;
	}
#ifdef PMAP
	abort();
	putc(this->aa_g,fp);
	xpos++;
#else
	if (this->aa_g != '*') {
	  putc(this->aa_g,fp);
	  xpos++;
	}
#endif
      }
    }
    putc('\n',fp);

  }

  return;
}

#ifdef PMAP
void
Pair_print_nucleotide_cdna (FILE *fp, struct T *ptr, int npairs, int wraplength) {
  struct T *this;
  int xpos = 0, i;

  for (i = 0; i < npairs; i++) {
    this = ptr++;
    if (this->cdna != ' ') {
      if (xpos == wraplength) {
	putc('\n',fp);
	xpos = 0;
      }
      putc(this->cdna,fp);
      xpos++;
    }
  }
  putc('\n',fp);
  return;
}
#else
void
Pair_print_protein_cdna (FILE *fp, struct T *ptr, int npairs, int wraplength, bool forwardp) {
  struct T *this;
  int xpos = 0, i;

  if (forwardp == true) {
    for (i = 0; i < npairs; i++) {
      this = ptr++;
      if (this->aa_e != ' ') {
	if (xpos == wraplength) {
	  putc('\n',fp);
	  xpos = 0;
	}
	if (this->aa_e != '*') {
	  putc(this->aa_e,fp);
	  xpos++;
	}
      }
    }
    putc('\n',fp);

  } else {
    for (i = npairs-1; i >= 0; i--) {
      this = ptr--;
      if (this->aa_e != ' ') {
	if (xpos == wraplength) {
	  putc('\n',fp);
	  xpos = 0;
	}
	if (this->aa_e != '*') {
	  putc(this->aa_e,fp);
	  xpos++;
	}
      }
    }
    putc('\n',fp);
  }

  return;
}
#endif


void
Pair_print_compressed (FILE *fp, int pathnum, int npaths, T start, T end, Sequence_T queryseq, char *dbversion,
		       Sequence_T usersegment, int nexons, double fracidentity,
		       struct T *pairs, int npairs, Chrnum_T chrnum,
		       Genomicpos_T chroffset, IIT_T chromosome_iit, int querylength_given,
		       int skiplength, int trim_start, int trim_end, bool checksump,
		       int chimerapos, int chimeraequivpos, double donor_prob, double acceptor_prob,
		       int chimera_cdna_direction, char *strain, bool watsonp, int cdna_direction) {
  Genomicpos_T chrpos1, chrpos2, position1, position2;

  bool in_exon = false;
  List_T tokens = NULL;
  struct T *ptr = pairs, *prev, *this = NULL;
  int querypos1, querypos2;
  int exon_querystart = -1, exon_genomestart = -1, exon_queryend, exon_genomeend,
    intron_start, intron_end;
  int num = 0, den = 0, runlength = 0, i;
  int print_dinucleotide_p;
  char token[10], donor[3], acceptor[3], *chr;
  double coverage;
  /* double trimmed_coverage; */
  int last_querypos = -1;
  Genomicpos_T last_genomepos = -1U;

  donor[0] = donor[1] = donor[2] = '\0';
  acceptor[0] = acceptor[1] = acceptor[2] = '\0';

  querypos1 = start->querypos;
  querypos2 = end->querypos;

  fprintf(fp,">%s ",Sequence_accession(queryseq));
  if (dbversion != NULL) {
    fprintf(fp,"%s ",dbversion);
  } else if (usersegment != NULL && Sequence_accession(usersegment) != NULL) {
    fprintf(fp,"%s ",Sequence_accession(usersegment));
  } else {
    fprintf(fp,"user-provided ");
  }
#ifdef PMAP
  fprintf(fp,"%d/%d %d %d",pathnum,npaths,(querylength_given+skiplength)*3,nexons);
  coverage = (double) (querypos2 - querypos1 + 1)/(double) ((querylength_given+skiplength)*3);
  fprintf(fp," %.1f",((double) rint(1000.0*coverage)));
#else
  coverage = (double) (querypos2 - querypos1 + 1)/(double) (querylength_given+skiplength);
  if (end->querypos + 1 > trim_end) {
    trim_end = end->querypos + 1;
  }
  if (start->querypos < trim_start) {
    trim_start = start->querypos;
  }
  /*
  trimmed_coverage = (double) (end->querypos - start->querypos + 1)/(double) (trim_end - trim_start + skiplength);
  fprintf(fp,">%s %s %d/%d %d(%d) %d",
	 Sequence_accession(queryseq),dbversion,pathnum,npaths,
	 querylength_given+skiplength,trim_end-trim_start,nexons);
  fprintf(fp," %.1f(%.1f)",((double) rint(1000.0*coverage))/10.0,((double) rint(1000.0*trimmed_coverage))/10.0);
  */
  fprintf(fp,"%d/%d %d %d",pathnum,npaths,querylength_given+skiplength,nexons);
  fprintf(fp," %.1f",((double) rint(1000.0*coverage))/10.0);
#endif
  fprintf(fp," %.1f",((double) rint(1000.0*fracidentity))/10.0);

  start = &(pairs[0]);
  end = &(pairs[npairs-1]);
  fprintf(fp," %d%s%d",start->querypos + ONEBASEDP,"..",end->querypos + ONEBASEDP);

  chrpos1 = start->genomepos;
  chrpos2 = end->genomepos;
  position1 = chroffset + chrpos1;
  position2 = chroffset + chrpos2;
  fprintf(fp," %u%s%u",position1 + ONEBASEDP,"..",position2 + ONEBASEDP);

  if (chrnum == 0) {
    fprintf(fp," %u%s%u",chrpos1 + ONEBASEDP,"..",chrpos2 + ONEBASEDP);
  } else {
    chr = Chrnum_to_string(chrnum,chromosome_iit);
    fprintf(fp," %s:%u%s%u",chr,chrpos1 + ONEBASEDP,"..",chrpos2 + ONEBASEDP);
    FREE(chr);
  }

  if (chrpos1 <= chrpos2) {
    fprintf(fp," +");
  } else {
    fprintf(fp," -");
  }

  if (cdna_direction > 0) {
    fprintf(fp," dir:sense");
  } else if (cdna_direction < 0) {
    fprintf(fp," dir:antisense");
  } else {
    fprintf(fp," dir:indet");
  }

  if (checksump == true) {
    fprintf(fp," md5:");
    Sequence_print_digest(fp,queryseq);
  }

  if (chimerapos >= 0) {
    if (chimeraequivpos == chimerapos) {
      if (donor_prob > 0.0 && acceptor_prob > 0.0) {
	if (chimera_cdna_direction >= 0) {
	  fprintf(fp," chimera:%d(>)/%.3f/%.3f",chimerapos + ONEBASEDP,donor_prob,acceptor_prob);
	} else {
	  fprintf(fp," chimera:%d(<)/%.3f/%.3f",chimerapos + ONEBASEDP,donor_prob,acceptor_prob);
	}
      } else {
	fprintf(fp," chimera:%d",chimerapos + ONEBASEDP);
      }
    } else {
      fprintf(fp," chimera:%d..%d",chimerapos + ONEBASEDP,chimeraequivpos + ONEBASEDP);
    }
  }

  if (strain != NULL) {
    fprintf(fp," strain:%s",strain);
  }

  putc('\n',fp);

  for (i = 0; i < npairs; i++) {
    prev = this;
    this = ptr++;

    if (this->gapp) {
      if (in_exon == true) {
	/* Beginning of gap */
	exon_queryend = last_querypos + ONEBASEDP;
	exon_genomeend = last_genomepos + ONEBASEDP;
	if (watsonp) {
	  intron_start = exon_genomeend + 1;
	} else {
	  intron_start = exon_genomeend - 1;
	}

	fprintf(fp,"\t%u %u",exon_genomestart,exon_genomeend);
	fprintf(fp," %d %d",exon_querystart,exon_queryend);
	if (den == 0) {
	  fprintf(fp," 100");
	} else {
	  fprintf(fp," %d",(int) floor(100.0*(double) num/(double) den));
	}
	print_dinucleotide_p = 1;
	if (this->comp == FWD_CANONICAL_INTRON_COMP) {
	  sprintf(token,"%d>",runlength);
	} else if (this->comp == REV_CANONICAL_INTRON_COMP) {
	  sprintf(token,"%d<",runlength);
	  print_dinucleotide_p = -1;
	} else if (this->comp == NONINTRON_COMP) {
	  sprintf(token,"%d=",runlength);
	} else if (this->comp == FWD_GCAG_INTRON_COMP) {
	  sprintf(token,"%d)",runlength);
	} else if (this->comp == REV_GCAG_INTRON_COMP) {
	  sprintf(token,"%d(",runlength);
	  print_dinucleotide_p = -1;
	} else if (this->comp == FWD_ATAC_INTRON_COMP) {
	  sprintf(token,"%d]",runlength);
	} else if (this->comp == REV_ATAC_INTRON_COMP) {
	  sprintf(token,"%d[",runlength);
	  print_dinucleotide_p = -1;
	} else if (this->comp == DUALBREAK_COMP) {
	  sprintf(token,"%d#",runlength);
	  print_dinucleotide_p = 0;
	} else if (this->comp == EXTRAEXON_COMP) {
	  sprintf(token,"%d#",runlength);
	  print_dinucleotide_p = 0;
	} else {
	  fprintf(stderr,"Can't parse comp '%c' in compression for %s\n",
		  this->comp,Sequence_accession(queryseq));
	  abort();
	}
	tokens = push_token(tokens,token);
	tokens = List_reverse(tokens);
	print_tokens_compressed(fp,tokens);
	List_free(&tokens);
	fprintf(fp,"\t%d",exon_queryend - exon_querystart + 1);

	runlength = 0;
	donor[0] = this->genome;
	donor[1] = '\0';
	in_exon = false;
      } else if (donor[1] == '\0') {
	donor[1] = this->genome;
      } else {
	acceptor[0] = acceptor[1];
	acceptor[1] = this->genome;
      }
    } else if (this->comp == INTRONGAP_COMP) {
      /* Do nothing */
    } else {
      /* Remaining possibilities are MATCH_COMP, DYNPROG_MATCH_COMP, AMBIGUOUS_COMP, INDEL_COMP, 
	 SHORTGAP_COMP, or MISMATCH_COMP */
      if (in_exon == false) {
	exon_querystart = this->querypos + ONEBASEDP;
	exon_genomestart = this->genomepos + ONEBASEDP;
	if (watsonp) {
	  intron_end = exon_genomestart - 1;
	} else {
	  intron_end = exon_genomestart + 1;
	}
	if (i > 0) {
	  if (intron_end > intron_start) {
	    fprintf(fp,"\t%d",intron_end - intron_start + 1);
	  } else {
	    fprintf(fp,"\t%d",intron_start - intron_end + 1);
	  }
	  if (print_dinucleotide_p == -1) {
	    invert_intron(donor,acceptor);
	  }
	  if (print_dinucleotide_p != 0) {
	    if ((donor[0] == 'G' || donor[0] == 'g') &&
		(donor[1] == 'T' || donor[1] == 't') &&
		(acceptor[0] == 'A' || acceptor[0] == 'a') &&
		(acceptor[1] == 'G' || acceptor[1] == 'g')) {
	      /* Do nothing */
	    } else {
	      fprintf(fp,"\t%c%c-%c%c",toupper(donor[0]),toupper(donor[1]),toupper(acceptor[0]),toupper(acceptor[1]));
	    }
	  }
#if 0
	  if (exon_querystart > exon_queryend + 1) {
	    fprintf(fp,"***");
	  }
#endif
	  putc('\n',fp);
	}

	num = den = 0;
	in_exon = true;
      }
      if (this->comp == INDEL_COMP || this->comp == SHORTGAP_COMP) {
	/* Gap in upper or lower sequence */
	if (this->genome == ' ') {
	  sprintf(token,"%d^%c",runlength,this->cdna);
	} else if (this->cdna == ' ') {
	  sprintf(token,"%dv",runlength);
	} else {
	  fprintf(stderr,"Error at %c%c%c\n",this->genome,this->comp,this->cdna);
	  exit(9);
	}
	tokens = push_token(tokens,token);
	runlength = 0;
	/* Don't increment den */

      } else if (this->comp == MISMATCH_COMP) {
	sprintf(token,"%dx%c",runlength,this->cdna);
	tokens = push_token(tokens,token);
	runlength = 0;
	den++;

#ifndef PMAP
      } else if (this->comp == AMBIGUOUS_COMP) {
	sprintf(token,"%d:%c",runlength,this->cdna);
	tokens = push_token(tokens,token);
	runlength = 0;
	den++;
	num++;
#endif

      } else {
	runlength++;
	den++;
	if (this->comp == MATCH_COMP || this->comp == DYNPROG_MATCH_COMP) {
	  /* AMBIGUOUS_COMP handled above */
	  num++;
	}
      }
    }

    if (this->cdna != ' ') {
      last_querypos = this->querypos;
    }
    if (this->genome != ' ') {
      last_genomepos = this->genomepos;
    }
  }
  
  prev = this;
  exon_queryend = last_querypos + ONEBASEDP;
  exon_genomeend = last_genomepos + ONEBASEDP;
  
  fprintf(fp,"\t%d %d",exon_genomestart,exon_genomeend);
  fprintf(fp," %d %d",exon_querystart,exon_queryend);
  if (den == 0) {
    fprintf(fp," 100");
  } else {
    fprintf(fp," %d",(int) floor(100.0*(double) num/(double) den));
  }

  sprintf(token,"%d*",runlength);
  tokens = push_token(tokens,token);
  tokens = List_reverse(tokens);
  print_tokens_compressed(fp,tokens);
  List_free(&tokens);

  fprintf(fp,"\t%d",exon_queryend - exon_querystart + 1);
  putc('\n',fp);

  return;
}


void
Pair_print_iit_map (FILE *fp, Sequence_T queryseq, char *accession,
		    T start, T end, Chrnum_T chrnum, IIT_T chromosome_iit) {
  char *chrstring = NULL;
  Genomicpos_T chrpos1, chrpos2;

  if (chrnum == 0) {
    chrstring = "";
  } else {
    chrstring = Chrnum_to_string(chrnum,chromosome_iit);
  }

  /* Made identical to code for Pair_print_iit_exon_map */
  chrpos1 = start->genomepos + ONEBASEDP;
  chrpos2 = end->genomepos + ONEBASEDP;
  fprintf(fp,">%s %s:%u..%u\n",accession,chrstring,chrpos1,chrpos2);
  Sequence_print_header(fp,queryseq,/*checksump*/false);

  if (chrnum != 0) {
    FREE(chrstring);
  }

  return;
}


void
Pair_print_iit_exon_map (FILE *fp, struct T *pairs, int npairs, Sequence_T queryseq, char *accession,
			 T start, T end, Chrnum_T chrnum, IIT_T chromosome_iit) {
  int i;
  bool in_exon = false;
  struct T *ptr = pairs, *prev, *this = NULL;
  int exon_genomestart = -1, exon_genomeend;
  char *chrstring = NULL;
  Genomicpos_T chrpos1, chrpos2;
  Genomicpos_T last_genomepos = -1U;

  if (chrnum == 0) {
    chrstring = "";
  } else {
    chrstring = Chrnum_to_string(chrnum,chromosome_iit);
  }

  chrpos1 = start->genomepos + ONEBASEDP;
  chrpos2 = end->genomepos + ONEBASEDP;
  fprintf(fp,">%s %s:%u..%u\n",accession,chrstring,chrpos1,chrpos2);
  Sequence_print_header(fp,queryseq,/*checksump*/false);

  for (i = 0; i < npairs; i++) {
    prev = this;
    this = ptr++;

    if (this->gapp) {
      if (in_exon == true) {
	/* Beginning of gap */
	exon_genomeend = last_genomepos + ONEBASEDP;
	fprintf(fp,"%u %u\n",exon_genomestart,exon_genomeend);
	in_exon = false;
      }
    } else if (this->comp == INTRONGAP_COMP) {
      /* Do nothing */
    } else {
      /* Remaining possibilities are MATCH_COMP, DYNPROG_MATCH_COMP, AMBIGUOUS_COMP, INDEL_COMP, 
	 SHORTGAP_COMP, or MISMATCH_COMP */
      if (in_exon == false) {
	exon_genomestart = this->genomepos + ONEBASEDP;
	in_exon = true;
      }
    }
    if (this->genome != ' ') {
      last_genomepos = this->genomepos;
    }
  }
  
  prev = this;
  exon_genomeend = last_genomepos + ONEBASEDP;
  
  fprintf(fp,"%u %u\n",exon_genomestart,exon_genomeend);

  if (chrnum != 0) {
    FREE(chrstring);
  }

  return;
}


void
Pair_print_splicesites (FILE *fp, struct T *pairs, int npairs, char *accession,
			int nexons, Chrnum_T chrnum, IIT_T chromosome_iit, bool watsonp) {
  int exoni = 0, i;
  bool in_exon = false;
  struct T *ptr = pairs, *prev, *this = NULL;
  int exon_genomestart = -1, exon_genomeend;
  char *chrstring = NULL;
  Genomicpos_T last_genomepos = -1U, intron_length;

  if (chrnum == 0) {
    chrstring = "";
  } else {
    chrstring = Chrnum_to_string(chrnum,chromosome_iit);
  }

  for (i = 0; i < npairs; i++) {
    prev = this;
    this = ptr++;

    if (this->gapp) {
      if (in_exon == true) {
	/* Beginning of gap */
	exon_genomeend = last_genomepos + ONEBASEDP;
	if (watsonp) {
	  fprintf(fp,">%s.exon%d/%d %s:%u..%u donor",accession,exoni,nexons,chrstring,exon_genomeend,exon_genomeend+1U);
	} else {
	  fprintf(fp,">%s.exon%d/%d %s:%u..%u donor",accession,exoni,nexons,chrstring,exon_genomeend,exon_genomeend-1U);
	}
	in_exon = false;
      }
    } else if (this->comp == INTRONGAP_COMP) {
      /* Do nothing */
    } else {
      /* Remaining possibilities are MATCH_COMP, DYNPROG_MATCH_COMP, AMBIGUOUS_COMP, INDEL_COMP, 
	 SHORTGAP_COMP, or MISMATCH_COMP */
      if (in_exon == false) {
	exoni++;
	if (exoni > 1) {
	  exon_genomestart = this->genomepos + ONEBASEDP;
	  if (watsonp) {
	    intron_length = exon_genomestart - exon_genomeend - 1U;
	    fprintf(fp," %u\n",intron_length); /* For previous donor */
	    fprintf(fp,">%s.exon%d/%d %s:%u..%u acceptor",accession,exoni,nexons,chrstring,exon_genomestart-1U,exon_genomestart);
	    fprintf(fp," %u\n",intron_length);
	  } else {
	    intron_length = exon_genomeend - exon_genomestart - 1U;
	    fprintf(fp," %u\n",intron_length); /* For previous donor */
	    fprintf(fp,">%s.exon%d/%d %s:%u..%u acceptor",accession,exoni,nexons,chrstring,exon_genomestart+1U,exon_genomestart);
	    fprintf(fp," %u\n",intron_length);
	  }
	}

	in_exon = true;
      }
    }
    if (this->genome != ' ') {
      last_genomepos = this->genomepos;
    }
  }
  
  if (chrnum != 0) {
    FREE(chrstring);
  }

  return;
}


void
Pair_print_introns (FILE *fp, struct T *pairs, int npairs, char *accession,
		    int nexons, Chrnum_T chrnum, IIT_T chromosome_iit) {
  int exoni = 0, i;
  bool in_exon = false;
  struct T *ptr = pairs, *prev, *this = NULL;
  int exon_genomestart = -1, exon_genomeend;
  char *chrstring = NULL;
  Genomicpos_T last_genomepos = -1U;

  if (chrnum == 0) {
    chrstring = "";
  } else {
    chrstring = Chrnum_to_string(chrnum,chromosome_iit);
  }

  for (i = 0; i < npairs; i++) {
    prev = this;
    this = ptr++;

    if (this->gapp) {
      if (in_exon == true) {
	/* Beginning of gap */
	exon_genomeend = last_genomepos + ONEBASEDP;
	in_exon = false;
      }
    } else if (this->comp == INTRONGAP_COMP) {
      /* Do nothing */
    } else {
      /* Remaining possibilities are MATCH_COMP, DYNPROG_MATCH_COMP, AMBIGUOUS_COMP, INDEL_COMP, 
	 SHORTGAP_COMP, or MISMATCH_COMP */
      if (in_exon == false) {
	exoni++;
	if (exoni > 1) {
	  exon_genomestart = this->genomepos + ONEBASEDP;
	  fprintf(fp,">%s.intron%d/%d %s:%u..%u\n",accession,exoni-1,nexons-1,chrstring,exon_genomeend,exon_genomestart);
	}

	in_exon = true;
      }
    }
    if (this->genome != ' ') {
      last_genomepos = this->genomepos;
    }
  }
  
  if (chrnum != 0) {
    FREE(chrstring);
  }

  return;
}



/* goal_start < goal_end */
Genomicpos_T
Pair_binary_search_ascending (int *querypos, int lowi, int highi, struct T *pairarray,
			      Genomicpos_T goal_start, Genomicpos_T goal_end) {
  int middlei;

  debug10(printf("entered binary search_ascending with lowi=%d, highi=%d, goal=%u..%u\n",
		 lowi,highi,goal_start,goal_end));

  while (lowi < highi) {
    middlei = (lowi+highi)/2;
    while (middlei < highi && pairarray[middlei].cdna == ' ') {
      /* Go forward past pairs corresponding to gaps */
      middlei++;
    }
    if (middlei >= highi) {
      middlei = (lowi+highi)/2;
      while (middlei >= lowi && pairarray[middlei].cdna == ' ') {
	/* Go backward past pairs corresponding to gaps */
	middlei--;
      }
      if (middlei < lowi) {
	debug10(printf("all intermediate pairs are gaps\n",lowi));
#if 0
	*querypos = pairarray[lowi].querypos;
	return pairarray[lowi].genomepos;
#else
	return 0U;
#endif
      }
    }

    debug10(printf("  binary: %d:%u %d:%u %d:%u   vs. %u..%u\n",
		   lowi,pairarray[lowi].genomepos,middlei,pairarray[middlei].genomepos,
		   highi,pairarray[highi].genomepos,goal_start,goal_end));
    if (goal_end < pairarray[middlei].genomepos) {
      highi = middlei;
    } else if (goal_start > pairarray[middlei].genomepos) {
      lowi = middlei + 1;
    } else {
      debug10(printf("binary search returns %d\n",middlei));
      *querypos = pairarray[middlei].querypos;
      return pairarray[middlei].genomepos;
    }
  }

  debug10(printf("binary search returns %d\n",highi));
  return 0U;
}

/* goal_start > goal_end */
Genomicpos_T
Pair_binary_search_descending (int *querypos, int lowi, int highi, struct T *pairarray,
			       Genomicpos_T goal_start, Genomicpos_T goal_end) {
  int middlei;

  debug10(printf("entered binary search_descending with lowi=%d, highi=%d, goal=%u..%u\n",
		 lowi,highi,goal_start,goal_end));

  while (lowi < highi) {
    middlei = (lowi+highi)/2;
    while (middlei < highi && pairarray[middlei].cdna == ' ') {
      /* Go forward past pairs corresponding to gaps */
      middlei++;
    }
    if (middlei >= highi) {
      middlei = (lowi+highi)/2;
      while (middlei >= lowi && pairarray[middlei].cdna == ' ') {
	/* Go backward past pairs corresponding to gaps */
	middlei--;
      }
      if (middlei < lowi) {
	debug10(printf("all intermediate pairs are gaps\n",lowi));
#if 0
	*querypos = pairarray[lowi].querypos;
	return pairarray[lowi].genomepos;
#else
	return 0U;
#endif
      }
    }

    debug10(printf("  binary: %d:%u %d:%u %d:%u   vs. %u..%u\n",
		   lowi,pairarray[lowi].genomepos,middlei,pairarray[middlei].genomepos,
		   highi,pairarray[highi].genomepos,goal_start,goal_end));
    if (goal_end > pairarray[middlei].genomepos) {
      highi = middlei;
    } else if (goal_start < pairarray[middlei].genomepos) {
      lowi = middlei + 1;
    } else {
      debug10(printf("binary search returns %d\n",middlei));
      *querypos = pairarray[middlei].querypos;
      return pairarray[middlei].genomepos;
    }
  }

  debug10(printf("binary search returns %d\n",highi));
  return 0U;
}


#ifndef PMAP

Genomicpos_T
Pair_genomicpos_low (int *hardclip_low, int *hardclip_high, struct T *pairarray, int npairs,
		     int querylength, bool watsonp) {
  struct T *clipped_pairs;
  int clipped_npairs;
  T pair;

  clipped_pairs = hardclip_pairs(&clipped_npairs,&(*hardclip_low),&(*hardclip_high),
				 pairarray,npairs,querylength);
  if (watsonp == true) {
    pair = &(clipped_pairs[0]);
    return pair->genomepos + 1U;
  } else {
    pair = &(clipped_pairs[clipped_npairs-1]);
    return pair->genomepos + 1U;
  }
}

#endif


Genomicpos_T
Pairarray_genomicbound_from_start (struct T *pairarray, int npairs, int overlap) {
  int i;
  struct T pair;

  i = 0;
  pair = pairarray[i];
  while (i < npairs && overlap > 0) {
    pair = pairarray[i];
    if (pair.cdna != ' ') {
      overlap--;
    }
    i++;
  }

  return pair.genomepos;
}

Genomicpos_T
Pairarray_genomicbound_from_end (struct T *pairarray, int npairs, int overlap) {
  int i;
  struct T pair;

  i = npairs-1;
  pair = pairarray[i];
  while (i >= 0 && overlap > 0) {
    pair = pairarray[i];
    if (pair.cdna != ' ') {
      overlap--;
    }
    i--;
  }

  return pair.genomepos;
}




List_T
Pair_trim_ends (List_T pairs) {
  List_T trimmed = NULL;
  int trim_right = 0, trim_left = -1; /* Needs to be -1 to avoid trimming when pairs is NULL */
  int bestscore, score;
  int pairi;
  List_T p, pairptr;
  T this;
  int i;
  bool in_indelp;

  debug8(printf("Entered trim_ends\n"));

  /* Find trim_right */
  bestscore = 0;
  score = 0;
  in_indelp = false;
  for (p = pairs, pairi = 0; p != NULL; p = p->rest, pairi++) {
    this = p->first;

    if (this->comp == INDEL_COMP || this->comp == SHORTGAP_COMP) {
      if (in_indelp == false) {
	score += trim_indel_score;
	if (score < 0) {
	  score = 0;
	}
	in_indelp = true;
      }

    } else {
      in_indelp = false;
      if (this->gapp) {
	/* Don't count */
	
      } else if (this->comp == INTRONGAP_COMP) {
	/* Do nothing */

      } else if (this->comp == MATCH_COMP || this->comp == DYNPROG_MATCH_COMP || this->comp == AMBIGUOUS_COMP) {
	score += TRIM_MATCH_SCORE;
	if (score >= bestscore) { /* Want >= and not >, so extend to ends */
	  bestscore = score;
	  trim_right = pairi;
	}

      } else if (this->comp == MISMATCH_COMP) {
	score += trim_mismatch_score;
	if (score < 0) {
	  score = 0;
	} else if (score >= bestscore) { /* Want >= and not >, so extend to ends */
	  bestscore = score;
	  trim_right = pairi;
	}
	
      } else {
	fprintf(stderr,"Can't parse comp %c, gapp %d\n",this->comp,this->gapp);
	abort();
      }
    }

    debug8(printf("pairi %d, querypos %d, genomepos %d, comp %c: Trim right score %d, trim_right %d, protectedp %d\n",
		  pairi,this->querypos,this->genomepos,this->comp,score,trim_right,this->protectedp));
  }

  if (this->protectedp == true) {
    debug8(printf("Protected against trim_right\n"));
    trim_right = 0;
  } else {
    trim_right = pairi - 1 - trim_right;
    debug8(printf("Final: Trim right pairi %d, score %d, trim_right %d\n",pairi,score,trim_right));
  }
  debug8(printf("\n"));


  /* Find trim_left */
  pairs = List_reverse(pairs);
  bestscore = 0;
  score = 0;
  in_indelp = false;
  for (p = pairs, pairi = 0; p != NULL; p = p->rest, pairi++) {
    this = p->first;

    if (this->comp == INDEL_COMP || this->comp == SHORTGAP_COMP) {
      if (in_indelp == false) {
	score += trim_indel_score;
	if (score < 0) {
	  score = 0;
	}
	in_indelp = true;
      }

    } else {
      in_indelp = false;

      if (this->gapp) {
	/* Don't count */
	
      } else if (this->comp == INTRONGAP_COMP) {
	/* Do nothing */

      } else if (this->comp == MATCH_COMP || this->comp == DYNPROG_MATCH_COMP || this->comp == AMBIGUOUS_COMP) {
	score += TRIM_MATCH_SCORE;
	if (score >= bestscore) { /* Want >= and not >, so extend to ends */
	  bestscore = score;
	  trim_left = pairi;
	}
	
      } else if (this->comp == MISMATCH_COMP) {
	score += trim_mismatch_score;
	if (score < 0) {
	  score = 0;
	} else if (score >= bestscore) { /* Want >= and not >, so extend to ends */
	  bestscore = score;
	  trim_left = pairi;
	}

      } else {
	fprintf(stderr,"Can't parse comp %c, gapp %d\n",this->comp,this->gapp);
	abort();
      }
    }

    debug8(printf("pairi %d, querypos %d, genomepos %d, comp %c: Trim left score %d, trim_left %d, protectedp %d\n",
		  pairi,this->querypos,this->genomepos,this->comp,score,trim_left,this->protectedp));
  }

  if (this->protectedp == true) {
    debug8(printf("Protected against trim_left\n"));
    trim_left = pairi - 1;
  } else {
    debug8(printf("Final: Trim left pairi %d, score %d, trim_left %d\n",pairi,score,trim_left));
  }
  debug8(printf("\n"));


  /* trim */
  i = 0;
  while (i < trim_right) {
    pairs = Pairpool_pop(pairs,&this);
    i++;
  }

  while (i <= trim_left) {
    pairptr = pairs;
    pairs = Pairpool_pop(pairs,&this);
#ifdef WASTE
    path = Pairpool_push_existing(path,pairpool,pair);
#else
    trimmed = List_push_existing(trimmed,pairptr);
#endif
    i++;
  }

  debug8(Pair_dump_list(trimmed,/*zerobasedp*/true));

  return trimmed;
}


#ifdef GSNAP
static int quality_score_adj = 33; /* Default is Sanger */

/* Taken from mapq.c */

static double
mismatch_logprob[MAX_QUALITY_SCORE+1] =
  /* log(1/3*10^(-Q/10)) */
  {-1.098612,
   -1.328871, -1.559129, -1.789388, -2.019646, -2.249905,
   -2.480163, -2.710422, -2.940680, -3.170939, -3.401197,
   -3.631456, -3.861714, -4.091973, -4.322231, -4.552490,
   -4.782748, -5.013007, -5.243265, -5.473524, -5.703782,
   -5.934041, -6.164299, -6.394558, -6.624817, -6.855075,
   -7.085334, -7.315592, -7.545851, -7.776109, -8.006368,
   -8.236626, -8.466885, -8.697143, -8.927402, -9.157660,
   -9.387919, -9.618177, -9.848436, -10.078694, -10.308953};



double
Pair_compute_mapq (struct T *pairarray, int npairs, char *quality_string) {
  double loglik = 0.0;
  int Q;
  T pair;
  int querypos, i;

  for (i = 0; i < npairs; i++) {
    pair = &(pairarray[i]);
    if (pair->comp == MISMATCH_COMP) {
      /* printf("Got a mismatch at querypos %d, cdna %c\n",pair->querypos,pair->cdna); */
      querypos = pair->querypos;
      Q = (quality_string == NULL) ? MAX_QUALITY_SCORE : quality_string[querypos] - quality_score_adj;
      if (Q < 0) {
	fprintf(stderr,"Warning: quality score %c (ASCII %d) - %d (quality-zero-score) = %d, which is less than 0.  May need to specify --quality-protocol or --quality-zero-score\n",
		quality_string[querypos],(int) quality_string[querypos],quality_score_adj, Q);
	/* fprintf(stderr,"Position %d in %d-%d of %s\n",querypos,querystart,queryend,quality_string); */
	Q = 0;
      } else if (Q > MAX_QUALITY_SCORE_INPUT) {
	fprintf(stderr,"Warning: quality score %c (ASCII %d) - %d (quality-zero-score) = %d, which exceeds %d.  May need to specify --quality-protocol or --quality-zero-score\n",
		quality_string[querypos],(int) quality_string[querypos],quality_score_adj,Q,MAX_QUALITY_SCORE_INPUT);
	/* fprintf(stderr,"Position %d in %d-%d of %s\n",querypos,querystart,queryend,quality_string); */
	Q = MAX_QUALITY_SCORE;
      } else if (Q > MAX_QUALITY_SCORE) {
	Q = MAX_QUALITY_SCORE;
      }

      loglik += mismatch_logprob[Q];
    }
  }

  return loglik;
}



Overlap_T
Pair_gene_overlap (struct T *ptr /* = pairarray */, int npairs, IIT_T genes_iit, int divno,
		   bool favor_multiexon_p) {
  struct T *this = NULL;
  bool in_exon = false;
  int i;
  Genomicpos_T low, high, exon_genomestart, exon_genomeend, last_genomepos = -1U;
  bool foundp = false;
  Overlap_T overlap;

  for (i = 0; i < npairs; i++) {
    this = ptr++;

    if (this->gapp) {
      if (in_exon == true) {
	exon_genomeend = last_genomepos;
	if (exon_genomestart <= exon_genomeend) {
	  low = exon_genomestart;
	  high = exon_genomeend;
	} else {
	  low = exon_genomeend;
	  high = exon_genomestart;
	}

	overlap = IIT_gene_overlap(genes_iit,divno,low,high,favor_multiexon_p);
	if (overlap == NO_KNOWN_GENE) {
	  /* Keep searching */
	} else if (overlap == KNOWN_GENE) {
	  if (favor_multiexon_p == true) {
	    /* Keep looking for a multiexon gene */
	    foundp = true;
	  } else {
	    return KNOWN_GENE;
	  }
	} else {
	  return KNOWN_GENE_MULTIEXON;
	}
	  
	in_exon = false;
      }

    } else if (this->comp == INTRONGAP_COMP) {
      /* Do nothing */

    } else {
      /* Remaining possibilities are MATCH_COMP, DYNPROG_MATCH_COMP, AMBIGUOUS_COMP, INDEL_COMP, 
	 SHORTGAP_COMP, or MISMATCH_COMP */
      if (in_exon == false) {
	exon_genomestart = this->genomepos;
	in_exon = true;
      }
    }
    if (this->genome != ' ') {
      last_genomepos = this->genomepos;
    }
  }

  exon_genomeend = last_genomepos;
  if (exon_genomestart <= exon_genomeend) {
    low = exon_genomestart;
    high = exon_genomeend;
  } else {
    low = exon_genomeend;
    high = exon_genomestart;
  }

  overlap = IIT_gene_overlap(genes_iit,divno,low,high,favor_multiexon_p);
  if (overlap == NO_KNOWN_GENE) {
    /* Keep searching */
  } else if (overlap == KNOWN_GENE) {
    if (favor_multiexon_p == true) {
      /* Keep looking for a multiexon gene */
      foundp = true;
    } else {
      return KNOWN_GENE;
    }
  } else {
    return KNOWN_GENE_MULTIEXON;
  }

  if (foundp == true) {
    return KNOWN_GENE;
  } else {
    return NO_KNOWN_GENE;
  }
}


void
Pair_init (int quality_score_adj_in) {
  quality_score_adj = quality_score_adj_in;
  return;
}

#endif


