static char rcsid[] = "$Id: pair.c,v 1.221 2007/09/11 20:45:01 twu Exp $";
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
print_top_ruler (int n, int npairs, int margin, int wraplength) {
  printf("%*d ",margin,n);
  if (n + wraplength < npairs) {
    printf("%s\n",RULER);
  } else {
    printf("%.*s\n",npairs-n,RULER);
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
print_cdna_sequence (struct T *ptr, int n, int npairs, bool zerobasedp,
		     int margin, int wraplength) {
  struct T *this;
  int i;

  this = ptr;
  printf("%*u ",margin,this->querypos + !zerobasedp);
  for (i = 0; n < npairs && i < wraplength; n++, i++) {
    this = ptr++;
    putchar(this->cdna);
  }
  printf("\n");
  return;
}

static int
find_aapos_in_line (int *aapos, struct T *ptr, int n, int npairs, int wraplength, 
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
print_peptide (struct T *ptr, int n, int npairs, int margin,
	       int wraplength, bool genomep) {
  struct T *this, *last = NULL;
  int aapos, i;

  if ((aapos = find_aapos_in_line(&aapos,ptr,n,npairs,wraplength,genomep)) < 0) {
    printf("%*s ",margin,"");
  } else {
    /* 4 is length of "aa.c" and "aa.g" */
    if (genomep == true) {
      printf("aa.g%*d ",margin-4,aapos);
    } else {
      printf("aa.c%*d ",margin-4,aapos);
    }
  }

  if (genomep == true) {
    for (i = 0; n < npairs && i < wraplength; n++, i++) {
      this = ptr++;
      putchar(this->aa_g);
    }
  } else {
    for (i = 0; n < npairs && i < wraplength; n++, i++) {
      this = ptr++;
      putchar(this->aa_e);
    }
  }

  printf("\n");
  return;
}

static void
print_alignment (struct T *ptr, int n, int npairs, bool diagnosticp, 
		 int margin, int wraplength) {
  struct T *this;
  int i;

  printf("%*s ",margin,"");
  for (i = 0; n < npairs && i < wraplength; n++, i++) {
    this = ptr++;
    
    if (diagnosticp == true) {
      /* Subtract 1 because dynprogindices start at +1 and -1 */
      if (this->comp == DYNPROG_MATCH_COMP) {
	if (this->dynprogindex > 0) {
	  printf("%c",(this->dynprogindex-1)%26+'a');
	} else if (this->dynprogindex < 0) {
	  printf("%c",(-this->dynprogindex-1)%26+'A');
	} else {
	  putchar(DYNPROG_MATCH_COMP);
	}
      } else if (this->shortexonp == true) {
	putchar(DIAGNOSTIC_SHORTEXON_COMP);
      } else {
	putchar(this->comp);
      }

    } else if (this->comp == DYNPROG_MATCH_COMP) {
      putchar(MATCH_COMP);
#ifndef PMAP
    } else if (this->comp == AMBIGUOUS_COMP) {
      putchar(MATCH_COMP);
#endif
    } else if (this->comp == SHORTGAP_COMP) {
      putchar(INDEL_COMP);
    } else if (this->comp == EXTRAEXON_COMP) {
      putchar(INTRONGAP_COMP);
    } else {
      putchar(this->comp);
    }
  }

  printf("\n");
  return;
}


static void
print_genomic_sequence (struct T *ptr, int n, int npairs, bool zerobasedp,
			char *chrstring, Genomicpos_T chrpos, Genomicpos_T chroffset,
			Genomicpos_T genomiclength, bool watsonp, bool universalp,
			int margin, int wraplength) {
  struct T *this;
  int i;
  char Buffer[100];

  this = ptr;
  if (chrstring == NULL || universalp) {
    if (watsonp) {
      sprintf(Buffer,"%u",chroffset+chrpos+this->genomepos + !zerobasedp);
    } else {
      sprintf(Buffer,"%u",chroffset+chrpos + (genomiclength - 1) - this->genomepos + !zerobasedp);
    }
  } else {
    if (watsonp) {
      sprintf(Buffer,"%s:%u",chrstring,chrpos+this->genomepos + !zerobasedp);
    } else {
      sprintf(Buffer,"%s:%u",chrstring,chrpos + (genomiclength - 1) - this->genomepos + !zerobasedp);
    }
  }
  printf("%*s ",margin,Buffer);
  for (i = 0; n < npairs && i < wraplength; n++, i++) {
    this = ptr++;
    if (this->comp == EXTRAEXON_COMP) {
      putchar(INTRONGAP_CHAR);
    } else {
      putchar(this->genome);
    }
  }
  printf("\n");
  return;
}

static int
compute_margin (struct T *start, struct T *end, bool zerobasedp, 
		char *chrstring, Genomicpos_T chrpos, Genomicpos_T chroffset,
		Genomicpos_T genomiclength, bool watsonp, bool universalp) {
  int margin;
  char Buffer[100];

  if (chrstring == NULL || universalp) {
    if (watsonp) {
      sprintf(Buffer,"%u",chroffset+chrpos+start->genomepos + !zerobasedp);
    } else {
      sprintf(Buffer,"%u",chroffset+chrpos + (genomiclength - 1) - start->genomepos + !zerobasedp);
    }
  } else {
    if (watsonp) {
      sprintf(Buffer,"%s:%u",chrstring,chrpos+start->genomepos + !zerobasedp);
    } else {
      sprintf(Buffer,"%s:%u",chrstring,chrpos + (genomiclength - 1) - start->genomepos + !zerobasedp);
    }
  }
  margin = strlen(Buffer) + 1;

  if (chrstring == NULL || universalp) {
    if (watsonp) {
      sprintf(Buffer,"%u",chroffset+chrpos+end->genomepos + !zerobasedp);
    } else {
      sprintf(Buffer,"%u",chroffset+chrpos + (genomiclength - 1) - end->genomepos + !zerobasedp);
    }
  } else {
    if (watsonp) {
      sprintf(Buffer,"%s:%u",chrstring,chrpos+end->genomepos + !zerobasedp);
    } else {
      sprintf(Buffer,"%s:%u",chrstring,chrpos + (genomiclength - 1) - end->genomepos + !zerobasedp);
    }
  }
  if (strlen(Buffer) + 1 > margin) {
    margin = strlen(Buffer) + 1;
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
    new[j].comp = complCode[old[i].comp];
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
    new[j].cdna = complCode[old[i].cdna];
    new[j].genome = complCode[old[i].genome];
    new[j].comp = complCode[old[i].comp];
  }
  return new;
}


static void
add_intronlengths (struct T *pairs, int npairs) {
  struct T *prev, *this = NULL, *ptr;
  int space, margin, i, j, k, gapstart, genomepos;
  char intronstring[20], cdnabreak[20], genomicbreak[20], comp;

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

      if (comp == DUALBREAK_COMP) {
	sprintf(cdnabreak,"%d",abs(this->querypos - prev->querypos)-1);
	sprintf(genomicbreak,"%d",abs(this->genomepos - prev->genomepos)-1);

	margin = (space - strlen(cdnabreak))/2;
	j = gapstart;
	while (margin > 0) {
	  ptr = &(pairs[j++]);
	  margin--;
	}
	for (k = 0; k < strlen(cdnabreak); k++) {
	  ptr = &(pairs[j++]);
	  ptr->cdna = cdnabreak[k];
	}

	margin = (space - strlen(genomicbreak))/2;
	j = gapstart;
	while (margin > 0) {
	  ptr = &(pairs[j++]);
	  margin--;
	}
	for (k = 0; k < strlen(genomicbreak); k++) {
	  ptr = &(pairs[j++]);
	  ptr->genome = genomicbreak[k];
	}

      } else {			/* Intron */
	sprintf(intronstring,"%d",abs(this->genomepos - prev->genomepos)-1);
	margin = (space - strlen(intronstring))/2;
	j = gapstart;
	while (margin > 0) {
	  ptr = &(pairs[j++]);
	  margin--;
	}
	for (k = 0; k < strlen(intronstring); k++) {
	  ptr = &(pairs[j++]);
	  ptr->cdna = intronstring[k];
	}
      }
    }
  }	
  return;
}


void
Pair_print_continuous (struct T *pairs, int npairs, Chrnum_T chrnum, Genomicpos_T chrpos,
		       Genomicpos_T chroffset, Genomicpos_T genomiclength,
		       bool watsonp, int cdna_direction, bool universalp, bool zerobasedp, bool diagnosticp, 
		       bool genomefirstp, int invertmode, bool nointronlenp) {
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
      putchar(this->genome);
    }
    printf("\n");

    ptr = pairs;
    for (n = 0; n < npairs; n++) {
      this = ptr++;
      if (this->comp == MATCH_COMP) {
	putchar(MATCH_COMP);
      } else if (diagnosticp == false && this->comp == DYNPROG_MATCH_COMP) {
	putchar(MATCH_COMP);
#ifdef PMAP
      } else if (diagnosticp == false && this->comp == AMBIGUOUS_COMP) {
	putchar(MATCH_COMP);
#endif
      } else {
	putchar(this->comp);
      }
    }
    printf("\n");

    ptr = pairs;
    for (n = 0; n < npairs; n++) {
      this = ptr++;
      putchar(this->cdna);
    }
    printf("\n");

  } else {
    ptr = pairs;
    for (n = 0; n < npairs; n++) {
      this = ptr++;
      putchar(this->cdna);
    }
    printf("\n");

    ptr = pairs;
    for (n = 0; n < npairs; n++) {
      this = ptr++;
      if (this->comp == MATCH_COMP) {
	putchar(MATCH_COMP);
      } else if (diagnosticp == false && this->comp == DYNPROG_MATCH_COMP) {
	putchar(MATCH_COMP);
#ifdef PMAP
      } else if (diagnosticp == false && this->comp == AMBIGUOUS_COMP) {
	putchar(MATCH_COMP);
#endif
      } else {
	putchar(this->comp);
      }
    }
    printf("\n");

    ptr = pairs;
    for (n = 0; n < npairs; n++) {
      this = ptr++;
      putchar(this->genome);
    }
    printf("\n");
  }

  if (save != NULL) {
    FREE(save);
  }
  return;
}  



void
Pair_print_alignment (struct T *pairs, int npairs, Chrnum_T chrnum, Genomicpos_T chrpos,
		      Genomicpos_T chroffset, IIT_T chromosome_iit, Genomicpos_T genomiclength, 
		      bool watsonp, int cdna_direction, bool universalp, bool zerobasedp, bool diagnosticp, 
		      bool genomicprop, int invertmode, bool nointronlenp, int wraplength) {
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
      chrstring = Chrnum_to_string(chrnum,chromosome_iit,/*allocp*/true);
    } else {
      chrstring = Chrnum_to_string_signed(chrnum,chromosome_iit,watsonp);
    }
  }

  margin = compute_margin(&(pairs[0]),&(pairs[npairs-1]),zerobasedp,chrstring,
			  chrpos,chroffset,genomiclength,watsonp,universalp);

  while (n < npairs) {
    print_top_ruler(n,npairs,margin,wraplength);
    print_peptide(ptr,n,npairs,margin,wraplength,/*genomep*/true);
    print_genomic_sequence(ptr,n,npairs,zerobasedp,chrstring,
			   chrpos,chroffset,genomiclength,
			   watsonp,universalp,margin,wraplength);
    print_alignment(ptr,n,npairs,diagnosticp,margin,wraplength);
    print_cdna_sequence(ptr,n,npairs,zerobasedp,margin,wraplength);
    print_peptide(ptr,n,npairs,margin,wraplength,/*genomep*/false);
    printf("\n");
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
Pair_print_pathsummary (int pathnum, T start, T end, Chrnum_T chrnum, Genomicpos_T chrpos,
			Genomicpos_T chroffset, IIT_T chromosome_iit, bool referencealignp, 
			IIT_T altstrain_iit, char *strain, IIT_T contig_iit, char *dbversion,
			int querylength_given, int skiplength, int trim_start, int trim_end, Genomicpos_T genomiclength,
			int nexons, int matches, int unknowns, int mismatches, 
			int qopens, int qindels, int topens, int tindels, int goodness,
			bool watsonp, int cdna_direction, double defect_rate, 
			int translation_start, int translation_end, int translation_length,
			int relaastart, int relaaend, bool zerobasedp, bool maponlyp,
			bool diagnosticp, int stage1_genomicstart, int stage1_genomiclength,
			double stage2_runtime, int stage2_indexsize, double stage3_runtime,
			double stage3_defectrate) {
  int querypos1, querypos2, den;
  double fracidentity, coverage, trimmed_coverage;
  Genomicpos_T position1, position2, chrpos1, chrpos2;
  char *refstrain, *comma1, *comma2;

  querypos1 = start->querypos;
  querypos2 = end->querypos;

  printf("  Path %d: ",pathnum);
  printf("query %d%s%d (%d bp) => ",
	 querypos1 + !zerobasedp,SEPARATOR,querypos2 + !zerobasedp,querypos2-querypos1+1);

  if (watsonp) {
    chrpos1 = chrpos + start->genomepos;
    chrpos2 = chrpos + end->genomepos;
  } else {
    chrpos1 = chrpos + (genomiclength - 1) - start->genomepos;
    chrpos2 = chrpos + (genomiclength - 1) - end->genomepos;
  }

  comma1 = Genomicpos_commafmt(chrpos1 + !zerobasedp);
  comma2 = Genomicpos_commafmt(chrpos2 + !zerobasedp);
  if (chrnum == 0) {
    if (watsonp) {
      printf("genomic %s%s%s (%d bp)\n",
	     comma1,SEPARATOR,comma2,chrpos2-chrpos1+1);
    } else {
      printf("genomic %s%s%s (%d bp)\n",
	     comma1,SEPARATOR,comma2,chrpos2-chrpos1-1);
    }
  } else {
    if (watsonp) {
      printf("chr %s:%s%s%s (%d bp)\n",
	     Chrnum_to_string(chrnum,chromosome_iit,/*allocp*/false),
	     comma1,SEPARATOR,comma2,chrpos2-chrpos1+1);
    } else {
      printf("chr %s:%s%s%s (%d bp)\n",
	     Chrnum_to_string(chrnum,chromosome_iit,/*allocp*/false),
	     comma1,SEPARATOR,comma2,chrpos2-chrpos1-1);
    }
  }
  FREE(comma2);
  FREE(comma1);

  if (maponlyp == false) {

    if (diagnosticp == true) {
      comma1 = Genomicpos_commafmt(stage1_genomicstart + !zerobasedp);
      comma2 = Genomicpos_commafmt(stage1_genomiclength + !zerobasedp);
      printf("    Stage 1 genomicstart: %s\n",comma1);
      printf("    Stage 1 genomiclength: %s\n",comma2);
      FREE(comma2);
      FREE(comma1);

      printf("    Stage 2 runtime: %.3f sec\n",stage2_runtime);
      printf("    Stage 2 indexsize: %d\n",stage2_indexsize);
      printf("    Stage 3 runtime: %.3f sec\n",stage3_runtime);
      printf("    Stage 3 defectrate: %f\n",stage3_defectrate);
    }

    printf("    cDNA direction: ");
    if (cdna_direction > 0) {
      printf("sense\n");
    } else if (cdna_direction < 0) {
      printf("antisense\n");
    } else {
      printf("indeterminate\n");
    }
  }

  if (altstrain_iit != NULL) {
    if (strain == NULL) {
      refstrain = IIT_typestring(altstrain_iit,/*straintype*/0);
      if (refstrain[0] == '\0') {
	/* Backward compatibility with old altstrain_iit */
	printf("    Strain: reference\n");
      } else {
	printf("    Strain: %s (reference)\n",refstrain);
      }
    } else {
      printf("    Strain: %s\n",strain);
    }
  }

  position1 = chroffset + chrpos1;
  position2 = chroffset + chrpos2;
  comma1 = Genomicpos_commafmt(position1 + !zerobasedp);
  comma2 = Genomicpos_commafmt(position2 + !zerobasedp);
  if (dbversion == NULL) {
    printf("    Genomic pos: %s%s%s",comma1,SEPARATOR,comma2);
  } else {
    printf("    Genomic pos: %s:%s%s%s",dbversion,comma1,SEPARATOR,comma2);
  }
  if (chrpos1 <= chrpos2) {
    printf(" (+ strand)\n");
  } else {
    printf(" (- strand)\n");
  }
  FREE(comma2);
  FREE(comma1);

  if (contig_iit != NULL) {
    if (position1 <= position2) {
      Segmentpos_print_accessions(contig_iit,position1,position2,referencealignp,strain,zerobasedp);
    } else {
      Segmentpos_print_accessions(contig_iit,position2,position1,referencealignp,strain,zerobasedp);
    }
  }
    
  if (maponlyp == false) {
    printf("    Number of exons: %d\n",nexons);

#ifdef PMAP
    coverage = (double) (matches + mismatches + qindels)/(double) (3*(querylength_given + skiplength));
    /* Can have coverage greater than given querylength because of added '*' at end */
    if (coverage > 1.0) {
      coverage = 1.0;
    }
#else
    coverage = (double) (matches + mismatches + qindels)/(double) (querylength_given + skiplength);
#endif
    printf("    Coverage: %.1f",((double) rint(1000.0*coverage))/10.0);
#ifdef PMAP
    printf(" (query length: %d aa)\n",querylength_given);
#else
    printf(" (query length: %d bp)\n",querylength_given);
    if (querypos2 > trim_end) {
      trim_end = querypos2;
    }
    if (querypos1 < trim_start) {
      trim_start = querypos1;
    }

    trimmed_coverage = (double) (querypos2 - querypos1 + 1)/(double) (trim_end - trim_start + 1 + skiplength);
    printf("    Trimmed coverage: %.1f",((double) rint(1000.0*trimmed_coverage))/10.0);
    printf(" (trimmed length: %d bp, trimmed region: %d..%d)",
	   trim_end-trim_start+1,trim_start+!zerobasedp,trim_end+!zerobasedp);
    printf("\n");
#endif

    if ((den = matches + mismatches + qindels + tindels) == 0) {
      fracidentity = 1.0;
    } else {
      fracidentity = (double) matches/(double) den;
    }

    debug(printf("     Goodness: %d\n",goodness));

    /* The definition of indels here should be consistent with Stage3_indels */
    printf("    Percent identity: %.1f (%d matches, %d mismatches, %d indels, %d unknowns)\n",
	   ((double) rint(1000.0*fracidentity))/10.0,matches,mismatches,qindels+tindels,unknowns);
    if (qindels + tindels > 0) {
      printf("    Non-intron gaps: %d openings, %d bases in cdna; %d openings, %d bases in genome\n",
	     qopens,qindels,topens,tindels);
    } 

#ifndef PMAP
    if (translation_length > 0) {
      if (cdna_direction >= 0) {
	printf("    Translation: %d..%d (%d aa)\n",
	       translation_start+!zerobasedp,translation_end+!zerobasedp,translation_length);
      } else {
	printf("    Translation: %d..%d (%d aa)\n",
	       translation_end+!zerobasedp,translation_start+!zerobasedp,translation_length);
      }
    } else if (relaastart > 0) {
      if (relaastart < relaaend) {
	printf("    Protein coords: %d..%d\n",relaastart,relaaend);
      } else {
	printf("    Protein coords: %d..%d\n",relaaend,relaastart);
      }
    }
#endif

    /* printf("    Defect rate (percent): %.1f\n",defect_rate*100.0); */

    /* printf("\n"); -- Done by caller */
  }

  return;
}


void
Pair_print_coordinates (struct T *pairs, int npairs, Chrnum_T chrnum, Genomicpos_T chrpos,
			Genomicpos_T chroffset, IIT_T chromosome_iit, Genomicpos_T genomiclength, 
			bool watsonp, bool zerobasedp, int invertmode) {
  T this;
  struct T *save = NULL, *ptr;
  int i;
  char *chrstring = NULL;

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
  
  if (chrnum != 0) {
    if (invertmode == 2) {
      chrstring = Chrnum_to_string(chrnum,chromosome_iit,/*allocp*/true);
    } else {
      chrstring = Chrnum_to_string_signed(chrnum,chromosome_iit,watsonp);
    }
  }

  if (watsonp) {
    for (i = 0; i < npairs; i++) {
      this = pairs++;
      if (this->gapp == false) {
#ifdef DEBUG1
	printf("%d %d %c\t",this->aapos,this->aaphase_e,this->aa_e);
#else
	if (this->aaphase_e != 0) {
	  printf("%d\t",this->aapos);
	} else {
	  printf("%d %c\t",this->aapos,this->aa_e);
	}
#endif
	printf("%d %c\t",this->querypos + !zerobasedp,this->cdna);
	if (chrstring == NULL) {
	  printf("%u %u %c",chrpos+this->genomepos + !zerobasedp,
		 chroffset+chrpos+this->genomepos + !zerobasedp,
		 this->genome);
	} else {
	  printf("%s:%u %u %c",chrstring,
		 chrpos+this->genomepos + !zerobasedp,
		 chroffset+chrpos+this->genomepos + !zerobasedp,
		 this->genome);
	}
#ifdef DEBUG1
	printf("\t%d %c",this->aaphase_g,this->aa_g);
#else
	if (this->aaphase_g != 0) {
	  printf("\t");
	} else {
	  printf("\t%c",this->aa_g);
	}
#endif
	printf("\n");
      }
    }
  } else {
    for (i = 0; i < npairs; i++) {
      this = pairs++;
      if (this->gapp == false) {
#ifdef DEBUG1
	printf("%d %d %c\t",this->aapos,this->aaphase_e,this->aa_e);
#else
	if (this->aaphase_e != 0) {
	  printf("%d\t",this->aapos);
	} else {
	  printf("%d %c\t",this->aapos,this->aa_e);
	}
#endif
	printf("%d %c\t",this->querypos + !zerobasedp,this->cdna);
	if (chrstring == NULL) {
	  printf("%u %u %c",chrpos+(genomiclength-1)-this->genomepos + !zerobasedp,
		 chroffset+chrpos+(genomiclength-1)-this->genomepos + !zerobasedp,
		 this->genome);
	} else {
	  printf("%s:%u %u %c",chrstring,
		 chrpos+(genomiclength-1)-this->genomepos + !zerobasedp,
		 chroffset+chrpos+(genomiclength-1)-this->genomepos + !zerobasedp,
		 this->genome);
	}
#ifdef DEBUG1
	printf("\t%d %c",this->aaphase_g,this->aa_g);
#else
	if (this->aaphase_g != 0) {
	  printf("\t");
	} else {
	  printf("\t%c",this->aa_g);
	}
#endif
	printf("\n");
      }
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
    case REV_GCAG_INTRON_COMP: printf("[ GC-AG"); break;
    case REV_CANONICAL_INTRON_COMP: printf("< GT-AG"); break;
    case SHORTGAP_COMP: printf("~ shortgap"); break;
    case NONINTRON_COMP: printf("= nonintron"); break;
    default: printf("? unknown"); break;
    }
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
      printf("|%c",(this->dynprogindex-1)%26+'a');
    } else if (this->dynprogindex < 0) {
      printf("|%c",(-this->dynprogindex-1)%26+'A');
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
	if (this->genomepos > prev_genomepos) {
	  printf("Problem at genomepos %d\n",this->genomepos);
	}
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
      if (this->genomepos < prev_genomepos) {
	fprintf(stderr,"Problem at genomepos %d\n",this->genomepos);
	result = true;
      }
      prev_querypos = this->querypos;
      prev_genomepos = this->genomepos;
    }
  }
  return result;
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
Pair_print_exonsummary (struct T *pairs, int npairs, Chrnum_T chrnum, Genomicpos_T chrpos,
			Genomicpos_T chroffset, IIT_T chromosome_iit, Genomicpos_T genomiclength,
			bool watsonp, bool universalp, bool zerobasedp, bool genomefirstp, 
			int invertmode) {
  bool in_exon = false;
  struct T *save = NULL, *ptr, *prev, *this = NULL;
  int exon_querystart = -1, exon_genomestart = -1, exon_queryend, exon_genomeend, intron_start, intron_end;
  int num = 0, den = 0, i;
  char *chrstring = NULL;

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

  debug(Pair_dump_array(pairs,npairs,zerobasedp));

  for (i = 0; i < npairs; i++) {
    prev = this;
    this = ptr++;

    if (this->gapp) {
      if (in_exon == true) {
	exon_queryend = prev->querypos + !zerobasedp;
	if (watsonp) {
	  exon_genomeend = chrpos + prev->genomepos + !zerobasedp;
	  intron_start = exon_genomeend + 1;
	} else {
	  exon_genomeend = chrpos + (genomiclength - 1) - prev->genomepos + !zerobasedp;
	  intron_start = exon_genomeend - 1;
	}
	if (genomefirstp == true) {
	  printf("    ");
	  if (chrnum == 0 || universalp) {
	    printf("%u-%u",chroffset+exon_genomestart,chroffset+exon_genomeend);
	  } else {
	    printf("%s:%d-%d",chrstring,exon_genomestart,exon_genomeend);
	  }
	  printf("  (%d-%d)",exon_querystart,exon_queryend);
	} else {
	  printf("    %d-%d",exon_querystart,exon_queryend);
	  printf("  ");
	  if (chrnum == 0 || universalp) {
	    printf("(%u-%u)",chroffset+exon_genomestart,chroffset+exon_genomeend);
	  } else {
	    printf("(%s:%d-%d)",chrstring,exon_genomestart,exon_genomeend);
	  }
	}
	if (den == 0) {
	  printf("   %d%%",100);
	} else {
	  printf("   %d%%",(int) floor(100.0*(double) num/(double) den));
	}
	if (this->comp == FWD_CANONICAL_INTRON_COMP) {
	  printf(" ->");
	} else if (this->comp == REV_CANONICAL_INTRON_COMP) {
	  printf(" <-");
	} else if (this->comp == FWD_GCAG_INTRON_COMP) {
	  printf(" -)");
	} else if (this->comp == REV_GCAG_INTRON_COMP) {
	  printf(" (-");
	} else if (this->comp == FWD_ATAC_INTRON_COMP) {
	  printf(" -]");
	} else if (this->comp == REV_ATAC_INTRON_COMP) {
	  printf(" [-");
	} else if (this->comp == NONINTRON_COMP) {
	  printf(" ==");
	} else {
	  printf(" ##");
	}
	in_exon = false;
      }
    } else if (this->comp == INTRONGAP_COMP) {
      /* Do nothing */
    } else {
      /* Remaining possibilities are MATCH_COMP, DYNPROG_MATCH_COMP, AMBIGUOUS_COMP, INDEL_COMP, 
	 SHORTGAP_COMP, or MISMATCH_COMP */
      if (in_exon == false) {
	exon_querystart = this->querypos + !zerobasedp;
	if (watsonp) {
	  exon_genomestart = chrpos + this->genomepos + !zerobasedp;
	  intron_end = exon_genomestart - 1;
	} else {
	  exon_genomestart = chrpos + (genomiclength - 1) - this->genomepos + !zerobasedp;
	  intron_end = exon_genomestart + 1;
	}
	if (i > 0) {
	  if (intron_end > intron_start) {
	    printf("   ...%d...",intron_end - intron_start + 1);
	  } else {
	    printf("   ...%d...",intron_start - intron_end + 1);
	  }

	  if (exon_querystart > exon_queryend + 1) {
	    printf("   ***query_skip:%d***",exon_querystart-(exon_queryend+1));
	  }

	  printf("\n");
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
	if (this->comp == MATCH_COMP || this->comp == DYNPROG_MATCH_COMP || this->comp == AMBIGUOUS_COMP) {
	  num++;
	}
      }
    }
  }

  prev = this;
  exon_queryend = prev->querypos + !zerobasedp;
  if (watsonp) {
    exon_genomeend = chrpos + prev->genomepos + !zerobasedp;
  } else {
    exon_genomeend = chrpos + (genomiclength - 1) - prev->genomepos + !zerobasedp;
  }
  if (genomefirstp == true) {
    printf("    ");
    if (chrnum == 0 || universalp) {
      printf("%u-%u",chroffset+exon_genomestart,chroffset+exon_genomeend);
    } else {
      printf("%s:%d-%d",chrstring,exon_genomestart,exon_genomeend);
    }
    printf("  (%d-%d)",exon_querystart,exon_queryend);
  } else {
    printf("    %d-%d",exon_querystart,exon_queryend);
    printf("  ");
    if (chrnum == 0 || universalp) {
      printf("(%u-%u)",chroffset+exon_genomestart,chroffset+exon_genomeend);
    } else {
      printf("(%s:%d-%d)",chrstring,exon_genomestart,exon_genomeend);
    }
  }
  if (den == 0) {
    printf("   %d%%",100);
  } else {
    printf("   %d%%",(int) floor(100.0*(double) num/(double) den));
  }
  printf("\n\n");

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
print_tokens_compressed (List_T tokens) {
  List_T p;
  int tokencount = 1;
  char *token, *lasttoken = NULL;
  
  for (p = tokens; p != NULL; p = List_next(p)) {
    token = (char *) List_head(p);
    if (lasttoken == NULL) {
      printf("\t%s",token);
      lasttoken = token;
    } else if (!strcmp(token,lasttoken)) {
      tokencount++;
    } else {
      if (tokencount > 1) {
	printf("!%d",tokencount);
      }
      printf(" %s",token);
      lasttoken = token;
      tokencount = 1;
    }
  }
  if (tokencount > 1) {
    printf("!%d",tokencount);
  }

  for (p = tokens; p != NULL; p = List_next(p)) {
    token = (char *) List_head(p);
    FREE(token);
  }

  return;
}

static void
print_tokens_gff3 (List_T tokens) {
  List_T p;
  char *token;
  
  if (tokens != NULL) {
    p = tokens;
    token = (char *) List_head(p);
    printf("%s",token);

    for (p = List_next(p); p != NULL; p = List_next(p)) {
      token = (char *) List_head(p);
      printf(" %s",token);
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
print_gff3_gene (int pathnum, char *dbversion, char *accession, char *chrstring, Genomicpos_T start_genomepos, 
		 Genomicpos_T end_genomepos, bool watsonp, int cdna_direction) {

  printf("%s\t",chrstring);	/* 1: seqid */
  printf("%s\t",dbversion);	/* 2: source */
  printf("gene\t");		/* 3: type */

  if (start_genomepos < end_genomepos) {
    printf("%u\t%u\t",start_genomepos,end_genomepos); /* 4,5: start, end */
  } else {
    printf("%u\t%u\t",end_genomepos,start_genomepos); /* 4,5: start, end */
  }

  printf(".\t");		/* 6: score */

  if (watsonp == true) {
    if (cdna_direction >= 0) {
      printf("+\t");
    } else {
      printf("-\t");
    }
  } else {
    if (cdna_direction >= 0) {
      printf("-\t");		/* 7: strand */
    } else {
      printf("+\t");
    }
  }

  printf(".\t");		/* 8: phase */

  /* 9: features */
  printf("ID=%s.path%d;Name=%s\n",accession,pathnum,accession);

  return;
}

static void
print_gff3_mrna (int pathnum, char *dbversion, char *accession, char *chrstring, Genomicpos_T start_genomepos, 
		 Genomicpos_T end_genomepos, int querylength_given, int skiplength,
		 int matches, int mismatches, int qindels, int tindels, 
		 bool watsonp, int cdna_direction) {
  int den;
  double coverage, fracidentity;

  printf("%s\t",chrstring);	/* 1: seqid */
  printf("%s\t",dbversion);	/* 2: source */
  printf("mRNA\t");		/* 3: type */
  if (start_genomepos < end_genomepos) {
    printf("%u\t%u\t",start_genomepos,end_genomepos); /* 4,5: start, end */
  } else {
    printf("%u\t%u\t",end_genomepos,start_genomepos); /* 4,5: start, end */
  }

  printf(".\t");		/* 6: score */

  if (watsonp == true) {
    if (cdna_direction >= 0) {
      printf("+\t");
    } else {
      printf("-\t");
    }
  } else {
    if (cdna_direction >= 0) {
      printf("-\t");		/* 7: strand */
    } else {
      printf("+\t");
    }
  }

  printf(".\t");		/* 8: phase */

  /* 9: features */
  printf("ID=%s.mrna%d;Name=%s;Parent=%s.path%d;",
	 accession,pathnum,accession,accession,pathnum);

#ifdef PMAP
    coverage = (double) (matches + mismatches + qindels)/(double) (3*(querylength_given + skiplength));
#else
    coverage = (double) (matches + mismatches + qindels)/(double) (querylength_given + skiplength);
#endif
  printf("Coverage=%.1f;",((double) rint(1000.0*coverage))/10.0);

  if ((den = matches + mismatches + qindels + tindels) == 0) {
    fracidentity = 1.0;
  } else {
    fracidentity = (double) matches/(double) den;
  }
  printf("Identity=%.1f",((double) rint(1000.0*fracidentity))/10.0);

  printf("\n");

  return;
}


static void
print_gff3_exon (int exonno, int pathnum, char *dbversion, char *accession, char *chrstring,
		 Genomicpos_T exon_genomestart, Genomicpos_T exon_genomeend,
		 int exon_querystart, int exon_queryend, bool watsonp, int cdna_direction,
		 int pctidentity) {

  printf("%s\t",chrstring);	/* 1: seqid */
  printf("%s\t",dbversion);	/* 2: source */
  printf("exon\t");		/* 3: type */
  if (exon_genomestart < exon_genomeend) {
    printf("%u\t%u\t",exon_genomestart,exon_genomeend); /* 4,5: start, end */
  } else {
    printf("%u\t%u\t",exon_genomeend,exon_genomestart); /* 4,5: start, end */
  }
  printf("%d\t",pctidentity);	/* 6: score */

  if (watsonp == true) {
    if (cdna_direction >= 0) {
      printf("+\t");
    } else {
      printf("-\t");
    }
  } else {
    if (cdna_direction >= 0) {
      printf("-\t");		/* 7: strand */
    } else {
      printf("+\t");
    }
  }

  printf(".\t");		/* 8: phase */

  /* 9: features */
  printf("ID=%s.mrna%d.exon%d;",accession,pathnum,exonno);
  printf("Name=%s;",accession);
  printf("Parent=%s.mrna%d;",accession,pathnum);
  if (cdna_direction >= 0) {
    printf("Target=%s %d %d +\n",accession,exon_querystart,exon_queryend);
  } else {
    printf("Target=%s %d %d -\n",accession,exon_queryend,exon_querystart);
  }

  return;
}

static void
print_gff3_cds (int cdsno, int pathnum, char *dbversion, char *accession, char *chrstring,
		Genomicpos_T cds_genomestart, Genomicpos_T cds_genomeend,
		int cds_querystart, int cds_queryend, bool watsonp, int cdna_direction,
		int pctidentity, int cds_phase) {

  printf("%s\t",chrstring);	/* 1: seqid */
  printf("%s\t",dbversion);	/* 2: source */
  printf("CDS\t");		/* 3: type */
  if (cds_genomestart < cds_genomeend) {
    printf("%u\t%u\t",cds_genomestart,cds_genomeend); /* 4,5: start, end */
  } else {
    printf("%u\t%u\t",cds_genomeend,cds_genomestart); /* 4,5: start, end */
  }
  printf("%d\t",pctidentity);	/* 6: score */

  if (watsonp == true) {
    if (cdna_direction >= 0) {
      printf("+\t");
    } else {
      printf("-\t");
    }
  } else {
    if (cdna_direction >= 0) {
      printf("-\t");		/* 7: strand */
    } else {
      printf("+\t");
    }
  }

  printf("%d\t",cds_phase);	/* 8: phase */

  /* 9: features */
  printf("ID=%s.mrna%d.cds%d;",accession,pathnum,cdsno);
  printf("Name=%s;",accession);
  printf("Parent=%s.mrna%d;",accession,pathnum);
  if (cdna_direction >= 0) {
    printf("Target=%s %d %d +\n",accession,cds_querystart,cds_queryend);
  } else {
    printf("Target=%s %d %d -\n",accession,cds_queryend,cds_querystart);
  }

  return;
}


static void
print_gff3_cdna_match (int cdsno, int pathnum, char *dbversion, char *accession, char *chrstring,
		       Genomicpos_T exon_genomestart, Genomicpos_T exon_genomeend,
		       int exon_querystart, int exon_queryend, bool watsonp, int cdna_direction,
		       int pctidentity, List_T tokens) {
  
  printf("%s\t",chrstring);	/* 1: seqid */
  printf("%s\t",dbversion);	/* 2: source */
  printf("cDNA_match\t");		/* 3: type */
  if (exon_genomestart < exon_genomeend) {
    printf("%u\t%u\t",exon_genomestart,exon_genomeend); /* 4,5: start, end */
  } else {
    printf("%u\t%u\t",exon_genomeend,exon_genomestart); /* 4,5: start, end */
  }
  printf("%d\t",pctidentity);	/* 6: score */

  /* 7: strand */
  if (watsonp == true) {
    printf("+\t");
  } else {
    printf("-\t");
  }

  printf(".\t");		/* 8: phase */

  /* 9: features */
  printf("ID=%s.path%d;",accession,pathnum);
  printf("Name=%s;",accession);
  printf("Target=%s %d %d;Gap=",accession,exon_querystart,exon_queryend);
  print_tokens_gff3(tokens);
  printf("\n");

  return;
}


static void
print_gff3_exons_forward (struct T *pairs, int npairs, int pathnum, char *dbversion, char *accession, char *chrstring,
			  Genomicpos_T chrpos, Genomicpos_T genomiclength, bool watsonp, int cdna_direction,
			  bool gff_introns_p, bool gff_gene_format_p) {
  bool in_exon = false;
  struct T *ptr, *prev, *this = NULL;
  int exon_querystart = -1, exon_genomestart = -1, exon_queryend, exon_genomeend, intron_start, intron_end;
  int pctidentity, num = 0, den = 0, exonno = 0, intronno = 0, i;
  int Mlength = 0, Ilength = 0, Dlength = 0;
  List_T tokens = NULL;
  char token[10];

  ptr = pairs;
  for (i = 0; i < npairs; i++) {
    prev = this;
    this = ptr++;

    if (this->gapp) {
      if (in_exon == true) {
	exon_queryend = prev->querypos + 1;
	if (watsonp) {
	  exon_genomeend = chrpos + prev->genomepos + 1;
	  intron_start = exon_genomeend + 1;
	} else {
	  exon_genomeend = chrpos + (genomiclength - 1) - prev->genomepos + 1;
	  intron_start = exon_genomeend - 1;
	}
	
	if (den == 0) {
	  pctidentity = 100;
	} else {
	  pctidentity = (int) floor(100.0*(double) num/(double) den);
	}
	
	if (gff_gene_format_p == true) {
	  print_gff3_exon(++exonno,pathnum,dbversion,accession,chrstring,exon_genomestart,exon_genomeend,
			  exon_querystart,exon_queryend,watsonp,cdna_direction,pctidentity);
	} else {
	  if (Mlength > 0) {
	    sprintf(token,"M%d",Mlength);
	    tokens = push_token(tokens,token);
	  } else if (Ilength > 0) {
	    sprintf(token,"I%d",Ilength);
	    tokens = push_token(tokens,token);
	  } else if (Dlength > 0) {
	    sprintf(token,"D%d",Ilength);
	    tokens = push_token(tokens,token);
	  }
	  tokens = List_reverse(tokens);
	  print_gff3_cdna_match(++exonno,pathnum,dbversion,accession,chrstring,exon_genomestart,exon_genomeend,
				exon_querystart,exon_queryend,watsonp,cdna_direction,pctidentity,tokens);
	  List_free(&tokens);
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
	if (watsonp) {
	  exon_genomestart = chrpos + this->genomepos + 1;
	  intron_end = exon_genomestart - 1;
	} else {
	  exon_genomestart = chrpos + (genomiclength - 1) - this->genomepos + 1;
	  intron_end = exon_genomestart + 1;
	}

	if (gff_introns_p == true) {
	  if (i > 0) {
#if 0
	    printf_gff3_intron(++intronno,pathnum,dbversion,accession,chrstring,?,?,intron_start,intron_end,watsonp);
#endif
	  }
	  printf("\n");
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

#ifndef PMAP
      } else if (unknown_base(this->cdna) || unknown_base(this->genome)) {
	/* Comp must be a space */
	/* Don't count in numerator or denominator */
#endif
      } else {

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

	den++;
	if (this->comp == MATCH_COMP || this->comp == DYNPROG_MATCH_COMP || this->comp == AMBIGUOUS_COMP) {
	  num++;
	}
      }
    }
  }

  prev = this;
  exon_queryend = prev->querypos + 1;
  if (watsonp) {
    exon_genomeend = chrpos + prev->genomepos + 1;
  } else {
    exon_genomeend = chrpos + (genomiclength - 1) - prev->genomepos + 1;
  }

  if (den == 0) {
    pctidentity = 100;
  } else {
    pctidentity = (int) floor(100.0*(double) num/(double) den);
  }
	
  if (gff_gene_format_p == true) {
    print_gff3_exon(++exonno,pathnum,dbversion,accession,chrstring,exon_genomestart,exon_genomeend,
		    exon_querystart,exon_queryend,watsonp,cdna_direction,pctidentity);
  } else {
    if (Mlength > 0) {
      sprintf(token,"M%d",Mlength);
      tokens = push_token(tokens,token);
    } else if (Ilength > 0) {
      sprintf(token,"I%d",Ilength);
      tokens = push_token(tokens,token);
    } else if (Dlength > 0) {
      sprintf(token,"D%d",Ilength);
      tokens = push_token(tokens,token);
    }
    tokens = List_reverse(tokens);
    print_gff3_cdna_match(++exonno,pathnum,dbversion,accession,chrstring,exon_genomestart,exon_genomeend,
			  exon_querystart,exon_queryend,watsonp,cdna_direction,pctidentity,tokens);
    List_free(&tokens);
  }

  return;
}

static void
print_gff3_exons_backward (struct T *pairs, int npairs, int pathnum, char *dbversion, char *accession, char *chrstring,
			   Genomicpos_T chrpos, Genomicpos_T genomiclength, bool watsonp, int cdna_direction,
			   bool gff_introns_p) {
  bool in_exon = false;
  struct T *ptr, *prev, *this = NULL;
  int exon_querystart = -1, exon_genomestart = -1, exon_queryend, exon_genomeend, intron_start, intron_end;
  int pctidentity, num = 0, den = 0, exonno = 0, intronno = 0, i;

  ptr = &(pairs[npairs-1]);
  for (i = npairs-1; i >= 0; i--) {
    prev = this;
    this = ptr--;

    if (this->gapp) {
      if (in_exon == true) {
	exon_queryend = prev->querypos + 1;
	if (watsonp) {
	  exon_genomeend = chrpos + prev->genomepos + 1;
	  intron_start = exon_genomeend + 1;
	} else {
	  exon_genomeend = chrpos + (genomiclength - 1) - prev->genomepos + 1;
	  intron_start = exon_genomeend - 1;
	}
	
	if (den == 0) {
	  pctidentity = 100;
	} else {
	  pctidentity = (int) floor(100.0*(double) num/(double) den);
	}
	
	print_gff3_exon(++exonno,pathnum,dbversion,accession,chrstring,exon_genomestart,exon_genomeend,
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
	if (watsonp) {
	  exon_genomestart = chrpos + this->genomepos + 1;
	  intron_end = exon_genomestart - 1;
	} else {
	  exon_genomestart = chrpos + (genomiclength - 1) - this->genomepos + 1;
	  intron_end = exon_genomestart + 1;
	}

	if (gff_introns_p == true) {
	  if (i > 0) {
#if 0
	    printf_gff3_intron(++intronno,pathnum,dbversion,accession,chrstring,?,?,intron_start,intron_end,watsonp);
#endif
	  }
	  printf("\n");
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
	if (this->comp == MATCH_COMP || this->comp == DYNPROG_MATCH_COMP || this->comp == AMBIGUOUS_COMP) {
	  num++;
	}
      }
    }
  }

  prev = this;
  exon_queryend = prev->querypos + 1;
  if (watsonp) {
    exon_genomeend = chrpos + prev->genomepos + 1;
  } else {
    exon_genomeend = chrpos + (genomiclength - 1) - prev->genomepos + 1;
  }

  if (den == 0) {
    pctidentity = 100;
  } else {
    pctidentity = (int) floor(100.0*(double) num/(double) den);
  }
	
  print_gff3_exon(++exonno,pathnum,dbversion,accession,chrstring,exon_genomestart,exon_genomeend,
		  exon_querystart,exon_queryend,watsonp,cdna_direction,pctidentity);
  return;
}


static void
print_gff3_cdss_forward (struct T *pairs, int npairs, int pathnum, char *dbversion, char *accession, char *chrstring,
			 Genomicpos_T chrpos, Genomicpos_T genomiclength, bool watsonp, int cdna_direction,
			 int translation_start, int translation_end) {
  bool in_cds = false;
  struct T *ptr, *prev, *this = NULL;
  int exon_querystart = -1, exon_genomestart = -1, exon_queryend, exon_genomeend, intron_start, intron_end, exon_phase;
  int pctidentity, num = 0, den = 0, cdsno = 0;

  ptr = pairs;
  while (ptr < &(pairs[npairs])) {
    prev = this;
    this = ptr++;

    if (in_cds == true) {
      if (this->aaphase_g == -1) {
	/* End of cds */
	exon_queryend = prev->querypos + 1;
	if (watsonp) {
	  exon_genomeend = chrpos + prev->genomepos + 1;
	  intron_start = exon_genomeend + 1;
	} else {
	  exon_genomeend = chrpos + (genomiclength - 1) - prev->genomepos + 1;
	  intron_start = exon_genomeend - 1;
	}
	  
	if (den == 0) {
	  pctidentity = 100;
	} else {
	  pctidentity = (int) floor(100.0*(double) num/(double) den);
	}
	
	print_gff3_cds(++cdsno,pathnum,dbversion,accession,chrstring,exon_genomestart,exon_genomeend,
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
	  if (this->comp == MATCH_COMP || this->comp == DYNPROG_MATCH_COMP || this->comp == AMBIGUOUS_COMP) {
	    num++;
	  }
	}
      }

    } else {
      if (this->aaphase_g == -1) {
	/* Continuation of non-cds */
      } else {
	/* Start of cds */
	exon_querystart = this->querypos + 1;
	exon_phase = this->aaphase_g; /* ? was aaphase_e */
	if (watsonp) {
	  exon_genomestart = chrpos + this->genomepos + 1;
	  intron_end = exon_genomestart - 1;
	} else {
	  exon_genomestart = chrpos + (genomiclength - 1) - this->genomepos + 1;
	  intron_end = exon_genomestart + 1;
	}

	num = den = 0;
	in_cds = true;
      }
    }
  }

  if (in_cds == true) {
    prev = this;
    exon_queryend = prev->querypos + 1;
    if (watsonp) {
      exon_genomeend = chrpos + prev->genomepos + 1;
    } else {
      exon_genomeend = chrpos + (genomiclength - 1) - prev->genomepos + 1;
    }

    if (den == 0) {
      pctidentity = 100;
    } else {
      pctidentity = (int) floor(100.0*(double) num/(double) den);
    }
	
    print_gff3_cds(++cdsno,pathnum,dbversion,accession,chrstring,exon_genomestart,exon_genomeend,
		   exon_querystart,exon_queryend,watsonp,cdna_direction,pctidentity,exon_phase);
  }

  return;
}

static void
print_gff3_cdss_backward (struct T *pairs, int npairs, int pathnum, char *dbversion, char *accession, char *chrstring,
			  Genomicpos_T chrpos, Genomicpos_T genomiclength, bool watsonp, int cdna_direction,
			  int translation_start, int translation_end) {
  bool in_cds = false;
  struct T *ptr, *prev, *this = NULL;
  int exon_querystart = -1, exon_genomestart = -1, exon_queryend, exon_genomeend, intron_start, intron_end, exon_phase;
  int pctidentity, num = 0, den = 0, cdsno = 0;

  /* printf("translation_start = %d, translation_end = %d\n",translation_start,translation_end); */

  ptr = &(pairs[npairs-1]);
  while (ptr >= &(pairs[0])) {
    prev = this;
    this = ptr--;

    if (in_cds == true) {
      if (this->aaphase_g == -1) {
	/* End of cds */
	exon_queryend = prev->querypos + 1;
	if (watsonp) {
	  exon_genomeend = chrpos + prev->genomepos + 1;
	  intron_start = exon_genomeend + 1;
	} else {
	  exon_genomeend = chrpos + (genomiclength - 1) - prev->genomepos + 1;
	  intron_start = exon_genomeend - 1;
	}
	  
	if (den == 0) {
	  pctidentity = 100;
	} else {
	  pctidentity = (int) floor(100.0*(double) num/(double) den);
	}
	
	print_gff3_cds(++cdsno,pathnum,dbversion,accession,chrstring,exon_genomestart,exon_genomeend,
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
	  if (this->comp == MATCH_COMP || this->comp == DYNPROG_MATCH_COMP || this->comp == AMBIGUOUS_COMP) {
	    num++;
	  }
	}
      }

    } else {
      if (this->aaphase_g == -1) {
	/* Continuation of non-cds */
      } else {
	/* Start of cds */
	exon_querystart = this->querypos + 1;
	exon_phase = this->aaphase_g; /* ? was aaphase_e */
	if (watsonp) {
	  exon_genomestart = chrpos + this->genomepos + 1;
	  intron_end = exon_genomestart - 1;
	} else {
	  exon_genomestart = chrpos + (genomiclength - 1) - this->genomepos + 1;
	  intron_end = exon_genomestart + 1;
	}

	num = den = 0;
	in_cds = true;
      }
    }
  }

  if (in_cds == true) {
    prev = this;
    exon_queryend = prev->querypos + 1;
    if (watsonp) {
      exon_genomeend = chrpos + prev->genomepos + 1;
    } else {
      exon_genomeend = chrpos + (genomiclength - 1) - prev->genomepos + 1;
    }

    if (den == 0) {
      pctidentity = 100;
    } else {
      pctidentity = (int) floor(100.0*(double) num/(double) den);
    }
	
    print_gff3_cds(++cdsno,pathnum,dbversion,accession,chrstring,exon_genomestart,exon_genomeend,
		   exon_querystart,exon_queryend,watsonp,cdna_direction,pctidentity,exon_phase);
  }

  return;
}


void
Pair_print_gff3 (struct T *pairs, int npairs, int pathnum, char *dbversion, char *accession, 
		 T start, T end, Chrnum_T chrnum, Genomicpos_T chrpos,
		 IIT_T chromosome_iit, Genomicpos_T genomiclength,
		 int translation_start, int translation_end, 
		 int querylength_given, int skiplength, int matches, int unknowns, int mismatches, 
		 int qopens, int qindels, int topens, int tindels, 
		 bool watsonp, int cdna_direction,
		 bool gff_gene_format_p, char *user_genomicseg) {
  char *chrstring = NULL, *source;
  Genomicpos_T chrpos1, chrpos2;
  struct T *ptr;

  if (chrnum == 0) {
    chrstring = "NA";
  } else {
    chrstring = Chrnum_to_string(chrnum,chromosome_iit,/*allocp*/false);
  }

  if (dbversion != NULL) {
    source = dbversion;
  } else if (user_genomicseg != NULL) {
    source = user_genomicseg;
  } else {
    source = "NA";
  }

  if (gff_gene_format_p == true) {
    if (watsonp) {
      chrpos1 = chrpos + start->genomepos;
      chrpos2 = chrpos + end->genomepos;
    } else {
      chrpos1 = chrpos + (genomiclength - 1) - start->genomepos;
      chrpos2 = chrpos + (genomiclength - 1) - end->genomepos;
    }
    print_gff3_gene(pathnum,source,accession,chrstring,chrpos1+1,chrpos2+1,watsonp,cdna_direction);
    print_gff3_mrna(pathnum,source,accession,chrstring,chrpos1+1,chrpos2+1,
		    querylength_given,skiplength,matches,mismatches,qindels,tindels,
		    watsonp,cdna_direction);

    if (cdna_direction >= 0) {
      print_gff3_exons_forward(pairs,npairs,pathnum,source,accession,chrstring,chrpos,genomiclength,watsonp,
			       cdna_direction,/*gff_introns_p*/false,/*gff_gene_format_p*/true);
      if (translation_end > 0) {
	print_gff3_cdss_forward(pairs,npairs,pathnum,source,accession,chrstring,chrpos,genomiclength,watsonp,
				cdna_direction,translation_start,translation_end);
      }
    } else {
      print_gff3_exons_backward(pairs,npairs,pathnum,source,accession,chrstring,chrpos,genomiclength,watsonp,
				cdna_direction,/*gff_introns_p*/false);
      if (translation_end > 0) {
	print_gff3_cdss_backward(pairs,npairs,pathnum,source,accession,chrstring,chrpos,genomiclength,watsonp,
				 cdna_direction,translation_start,translation_end);
      }
    }

    printf("###\n");		/* Terminates gene format */
  } else {
    print_gff3_exons_forward(pairs,npairs,pathnum,source,accession,chrstring,chrpos,genomiclength,watsonp,
			     cdna_direction,/*gff_introns_p*/false,/*gff_gene_format_p*/false);
  }

  return;
}


Uintlist_T
Pair_exonbounds (struct T *pairs, int npairs, Genomicpos_T chrpos,
		 Genomicpos_T chroffset, Genomicpos_T genomiclength, bool watsonp) {
  Uintlist_T exonbounds = NULL;
  struct T *ptr, *prev, *this = NULL;
  bool in_exon = false;
  int i;

  ptr = pairs;
  for (i = 0; i < npairs; i++) {
    prev = this;
    this = ptr++;

    if (this->gapp) {
      if (in_exon == true) {
	/* exon genomeend */
	if (watsonp) {
	  exonbounds = Uintlist_push(exonbounds,chroffset + chrpos + prev->genomepos);
	} else {
	  exonbounds = Uintlist_push(exonbounds,chroffset + chrpos + (genomiclength - 1) - prev->genomepos);
	}
	in_exon = false;
      }
    } else if (this->comp == INTRONGAP_COMP) {
      /* Do nothing */
    } else {
      /* Remaining possibilities are MATCH_COMP, DYNPROG_MATCH_COMP, AMBIGUOUS_COMP, INDEL_COMP, 
	 SHORTGAP_COMP, or MISMATCH_COMP */
      if (in_exon == false) {
	/* exon genomestart */
	if (watsonp) {
	  exonbounds = Uintlist_push(exonbounds,chroffset + chrpos + this->genomepos);
	} else {
	  exonbounds = Uintlist_push(exonbounds,chroffset + chrpos + (genomiclength - 1) - this->genomepos);
	}
	in_exon = true;
      }
    }
  }

  prev = this;
  if (watsonp) {
    exonbounds = Uintlist_push(exonbounds,chroffset + chrpos + prev->genomepos);
  } else {
    exonbounds = Uintlist_push(exonbounds,chroffset + chrpos + (genomiclength - 1) - prev->genomepos);
  }

  return Uintlist_reverse(exonbounds);
}


static int
count_psl_blocks_nt (Intlist_T *blockSizes, Intlist_T *qStarts, Uintlist_T *tStarts, struct T *pairs_directional,
		     int npairs, int querylength, Genomicpos_T chrpos, Genomicpos_T genomiclength, bool watsonp) {
  int nblocks = 0, i;
  int block_querystart, block_queryend;
  Genomicpos_T block_genomestart;
  struct T *ptr = pairs_directional, *prev, *this = NULL;
  bool in_block = false;

  for (i = 0; i < npairs; i++) {
    prev = this;
    this = ptr++;

    if (this->gapp) {
      if (in_block == true) {
	nblocks++;
	block_queryend = prev->querypos;
	*blockSizes = Intlist_push(*blockSizes,abs(block_queryend-block_querystart)+1);
	in_block = false;
      }
    } else if (this->comp == INTRONGAP_COMP) {
      /* Do nothing */

    } else if (this->cdna == ' ' || this->genome == ' ') {
      if (in_block == true) {
	nblocks++;
	block_queryend = prev->querypos;
	*blockSizes = Intlist_push(*blockSizes,abs(block_queryend-block_querystart)+1);
	in_block = false;
      }

    } else {
      /* Remaining possibilities are MATCH_COMP, DYNPROG_MATCH_COMP, AMBIGUOUS_COMP, INDEL_COMP, 
         or SHORTGAP_COMP */
      if (in_block == false) {
	block_querystart = this->querypos;
	if (watsonp == true) {
	  *qStarts = Intlist_push(*qStarts,block_querystart);
	  *tStarts = Uintlist_push(*tStarts,chrpos + this->genomepos);
	} else {
	  *qStarts = Intlist_push(*qStarts,querylength-block_querystart-1);
	  *tStarts = Uintlist_push(*tStarts,chrpos + (genomiclength - 1) - this->genomepos);
	}
	in_block = true;
      }
    }
  }

  prev = this;
  nblocks++;
  block_queryend = prev->querypos;
  *blockSizes = Intlist_push(*blockSizes,abs(block_queryend-block_querystart)+1);

  *blockSizes = Intlist_reverse(*blockSizes);
  *qStarts = Intlist_reverse(*qStarts);
  *tStarts = Uintlist_reverse(*tStarts);

  return nblocks;
}

static int
count_psl_blocks_pro (Intlist_T *blockSizes, Intlist_T *qStarts, Uintlist_T *tStarts, struct T *pairs_directional,
		      int npairs, Genomicpos_T chrpos, Genomicpos_T genomiclength, bool watsonp,
		      Genomicpos_T chrlength) {
  int nblocks = 0, i;
  int naminoacids = 0;
  int block_querystart;
  Genomicpos_T block_genomestart;
  struct T *ptr = pairs_directional, *prev, *this = NULL;
  bool in_block = false;

  for (i = 0; i < npairs; i++) {
    prev = this;
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
	block_queryend = prev->querypos;
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
	  *tStarts = Uintlist_push(*tStarts,chrpos + this->genomepos);
	} else {
#if 0
	  /* Should be this */
	  *tStarts = Uintlist_push(*tStarts,chrpos + (genomiclength - 1) - this->genomepos);
#else
	  /* But is actually this */
	  *tStarts = Uintlist_push(*tStarts,chrlength - (chrpos + (genomiclength - 1) - this->genomepos) - 1);
#endif
	}
	in_block = true;
      }
    }
  }

  prev = this;
  nblocks++;
  *blockSizes = Intlist_push(*blockSizes,naminoacids);

  *blockSizes = Intlist_reverse(*blockSizes);
  *qStarts = Intlist_reverse(*qStarts);
  *tStarts = Uintlist_reverse(*tStarts);

  return nblocks;
}


static void
compute_gap_lengths_int (int *nbreaks, int *length, Intlist_T blockSizes, Intlist_T Starts, int nblocks) {
  int i;
  int start, end;
  Intlist_T p = blockSizes, q = Starts;

  *nbreaks = *length = 0;
  for (i = 0; i < nblocks - 1; i++) {
    if (i > 0) {
      start = Intlist_head(Starts);
      if (start - end > 0) {
	*nbreaks += 1;
	*length += (start - end);
      }
      debug2(printf("%d - %d = %d, gap = %d\n",start,end,start-end,*length));
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
    debug2(printf("%d - %d = %d, gap = %d\n",start,end,start-end,*length));
  }

  return;
}

static void
compute_gap_lengths_uint (int *nbreaks, int *length, Intlist_T blockSizes, Uintlist_T Starts, int nblocks) {
  int i;
  int start, end;
  Intlist_T p = blockSizes;
  Uintlist_T q = Starts;

  *nbreaks = *length = 0;
  for (i = 0; i < nblocks - 1; i++) {
    if (i > 0) {
      start = Uintlist_head(Starts);
      if (start - end > 0) {
	*nbreaks += 1;
	*length += (start - end);
      }
      debug2(printf("%d - %d = %d, gap = %d\n",start,end,start-end,*length));
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
    debug2(printf("%d - %d = %d, gap = %d\n",start,end,start-end,*length));
  }

  return;
}



static void
count_matches_pro (int *matches, int *mismatches, int *unknowns, 
		   struct T *pairs, int npairs) {
  struct T *prev, *this = NULL;
  int i;

  i = 0;
  while (i < npairs) {
    prev = this;
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
Pair_print_pslformat_nt (struct T *pairs, int npairs, T start, T end,
			 Sequence_T queryseq, Chrnum_T chrnum, Genomicpos_T chrpos,
			 IIT_T chromosome_iit, Genomicpos_T genomiclength,
			 int nexons, int matches, int unknowns, int mismatches, 
			 int qopens, int qindels, int topens, int tindels,
			 bool watsonp, int cdna_direction) {
  Genomicpos_T chrpos1, chrpos2;
  struct T *pairs_directional = NULL;
  Intlist_T blockSizes = NULL, qStarts = NULL, p;
  Uintlist_T tStarts = NULL, q;
  int nblocks;
  int qnbreaks, qlength, tnbreaks, tlength, querylength;

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
				querylength,chrpos,genomiclength,watsonp);
  compute_gap_lengths_int(&qnbreaks,&qlength,blockSizes,qStarts,nblocks);
  compute_gap_lengths_uint(&tnbreaks,&tlength,blockSizes,tStarts,nblocks);

  printf("%d\t%d\t%d\t%d\t",matches,mismatches,/*repeatmatches*/0,unknowns);
  printf("%d\t%d\t%d\t%d\t",qnbreaks,qlength,tnbreaks,tlength);

  if (watsonp == true) {
    printf("+");
  } else {
    printf("-");
  }
  printf("\t%s\t%d",Sequence_accession(queryseq),Sequence_fulllength_given(queryseq));

  printf("\t%d\t%d",start->querypos,end->querypos+1);

  /* T name and T size */
  printf("\t%s\t%u",
	 Chrnum_to_string(chrnum,chromosome_iit,/*allocp*/false),
	 Chrnum_length(chrnum,chromosome_iit));

  /* T start and T end */
  if (watsonp) {
    chrpos1 = chrpos + start->genomepos;
    chrpos2 = chrpos + end->genomepos;
    printf("\t%u\t%u",chrpos1,chrpos2+1U);
  } else {
    chrpos1 = chrpos + (genomiclength - 1) - start->genomepos;
    chrpos2 = chrpos + (genomiclength - 1) - end->genomepos;
    printf("\t%u\t%u",chrpos2,chrpos1+1U);
  }

  printf("\t%d",nblocks);

  printf("\t");
  for (p = blockSizes; p != NULL; p = Intlist_next(p)) {
    printf("%d,",Intlist_head(p));
  }

  printf("\t");
  for (p = qStarts; p != NULL; p = Intlist_next(p)) {
    printf("%d,",Intlist_head(p));
  }

  printf("\t");
  for (q = tStarts; q != NULL; q = Uintlist_next(q)) {
    printf("%u,",Uintlist_head(q));
  }

  Intlist_free(&blockSizes);
  Intlist_free(&qStarts);
  Uintlist_free(&tStarts);

  if (watsonp == false) {
    FREE(pairs_directional);
  }

  printf("\n");
  return;
}

void
Pair_print_pslformat_pro (struct T *pairs, int npairs, T start, T end,
			  Sequence_T queryseq, Chrnum_T chrnum, Genomicpos_T chrpos,
			  IIT_T chromosome_iit, Genomicpos_T genomiclength, int nexons, 
			  int qopens, int qindels, int topens, int tindels,
			  bool watsonp, int cdna_direction) {
  Genomicpos_T chrpos1, chrpos2, chrlength;
  Intlist_T blockSizes = NULL, qStarts = NULL, p;
  Uintlist_T tStarts = NULL, q;
  int nblocks, matches = 0, mismatches = 0, unknowns = 0;
  int qnbreaks, qlength, tnbreaks, tlength;

  chrlength = Chrnum_length(chrnum,chromosome_iit);
  nblocks = count_psl_blocks_pro(&blockSizes,&qStarts,&tStarts,pairs,npairs,
				 chrpos,genomiclength,watsonp,chrlength);
  compute_gap_lengths_int(&qnbreaks,&qlength,blockSizes,qStarts,nblocks);
  compute_gap_lengths_uint(&tnbreaks,&tlength,blockSizes,tStarts,nblocks);

  count_matches_pro(&matches,&mismatches,&unknowns,pairs,npairs);

  printf("%d\t%d\t%d\t%d\t",matches,mismatches,/*repeatmatches*/0,unknowns);
  printf("%d\t%d\t%d\t%d\t",qnbreaks,qlength,tnbreaks,tlength);

  if (cdna_direction >= 0) {
    printf("+");
  } else {
    printf("-");
  }

  if (watsonp == true) {
    printf("+");
  } else {
    printf("-");
  }
  printf("\t%s\t%d",Sequence_accession(queryseq),Sequence_fulllength_given(queryseq));

  printf("\t%d\t%d",(start->querypos+2)/3,end->querypos/3+1);

  /* T name and T size */
  printf("\tchr%s\t%u",
	 Chrnum_to_string(chrnum,chromosome_iit,/*allocp*/false),
	 Chrnum_length(chrnum,chromosome_iit));

  /* T start and T end */
  if (watsonp) {
    chrpos1 = chrpos + start->genomepos;
    chrpos2 = chrpos + end->genomepos;
    printf("\t%u\t%u",chrpos1,chrpos2+1U);
  } else {
    chrpos1 = chrpos + (genomiclength - 1) - start->genomepos;
    chrpos2 = chrpos + (genomiclength - 1) - end->genomepos;
    printf("\t%u\t%u",chrpos2,chrpos1+1U);
  }

  nblocks = count_psl_blocks_pro(&blockSizes,&qStarts,&tStarts,pairs,npairs,
				 chrpos,genomiclength,watsonp,chrlength);
  printf("\t%d",nblocks);
  printf("\t");

  for (p = blockSizes; p != NULL; p = Intlist_next(p)) {
    printf("%d,",Intlist_head(p));
  }

  printf("\t");
  for (p = qStarts; p != NULL; p = Intlist_next(p)) {
    printf("%d,",Intlist_head(p));
  }

  printf("\t");

  for (q = tStarts; q != NULL; q = Uintlist_next(q)) {
    printf("%u,",Uintlist_head(q));
  }

  Intlist_free(&blockSizes);
  Intlist_free(&qStarts);
  Uintlist_free(&tStarts);

  printf("\n");
  return;
}

void
Pair_print_exons (struct T *pairs, int npairs, int wraplength, int ngap, bool cdnap) {
  bool in_exon = false;
  struct T *ptr, *this = NULL;
  int i, exonno = 0, column = 0;

  ptr = pairs;
  for (i = 0; i < npairs; i++) {
    this = ptr++;

    if (this->gapp) {
      if (in_exon == true) {
	if (column != 0) {
	  printf("\n");
	  column = 0;
	}
	printf("</exon>\n");
	in_exon = false;
	if (ngap > 0) {
	  printf("<intron %d>\n",exonno);
	  putchar(this->genome);
	  column = 1;
	}
      } else {
	if (ngap > 0) {
	  putchar(this->genome);
	  if (++column % wraplength == 0) {
	    printf("\n");
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
	      printf("\n");
	      column = 0;
	    }
	    printf("</intron>\n");
	  }
	}
	printf("<exon %d>\n",++exonno);
	in_exon = true;
      }
      if (cdnap == true) {
	if (this->cdna != ' ') {
	  putchar(this->cdna);
	  if (++column % wraplength == 0) {
	    printf("\n");
	    column = 0;
	  }
	}
      } else {
	if (this->genome != ' ') {
	  putchar(this->genome);
	  if (++column % wraplength == 0) {
	    printf("\n");
	    column = 0;
	  }
	}
      }
    }
  }
  if (column != 0) {
    printf("\n");
  }
  printf("</exon>\n");

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
Pair_matchscores (struct T *ptr, int npairs, 
		  int cdna_direction, int querylength) {
  int *matchscores;
  int querypos;
  T this;
  int i;

  matchscores = (int *) CALLOC(querylength,sizeof(int));
  /* Make default to be a match, so any sequence oddity shows up as a mismatch */
  for (i = 0; i < querylength; i++) {
    matchscores[i] = 1;
  }

  for (i = 0; i < npairs; i++) {
    this = ptr++;

    querypos = this->querypos;
    if (querypos >= querylength) {
      fprintf(stderr,"Pair_matchscores: querypos %d >= querylength %d\n",querypos,querylength);
      RAISE(Array_bounds_error);
    }

    if (!this->gapp) {
      if (this->comp != MATCH_COMP && this->comp != DYNPROG_MATCH_COMP) {
	matchscores[querypos] = 0; /* For mismatch */
      }
    }
  }

  return matchscores;
}


void
Pair_matchscores_list (int *matchscores, List_T pairs, int npairs, int ngaps, int matchespergap, int cdna_direction) {
  T this;
  List_T p;
  int i = 0, j;

  for (p = pairs; p != NULL; p = p->rest) {
    this = p->first;
    if (!this->gapp) {
      if (this->comp == MISMATCH_COMP || this->comp == INDEL_COMP) {
	matchscores[i++] = 0; /* For mismatch */
      } else {
	matchscores[i++] = 1; /* For match */
      }
    } else if (cdna_direction > 0) {
      if (this->comp == FWD_CANONICAL_INTRON_COMP || this->comp == FWD_GCAG_INTRON_COMP || this->comp == FWD_ATAC_INTRON_COMP) {
	for (j = 0; j < matchespergap; j++) {
	  matchscores[i++] = 1;
	}
      } else {
	for (j = 0; j < matchespergap; j++) {
	  matchscores[i++] = 0;
	}
      }
    } else if (cdna_direction < 0) {
      if (this->comp == REV_CANONICAL_INTRON_COMP || this->comp == REV_GCAG_INTRON_COMP || this->comp == REV_ATAC_INTRON_COMP) {
	for (j = 0; j < matchespergap; j++) {
	  matchscores[i++] = 1;
	}
      } else {
	for (j = 0; j < matchespergap; j++) {
	  matchscores[i++] = 0;
	}
      }
    } else {
      abort();
    }
  }

  if (i != npairs + ngaps*(matchespergap-1)) {
    abort();
  }

  return;
}


void
Pair_pathscores (int *pathscores, struct T *ptr, int npairs, 
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
      if (!in_intron) {
	if (cdna_direction > 0) {
	  if (this->comp == FWD_CANONICAL_INTRON_COMP) {
	    pathscores[querypos] = CANONICAL_POINTS;
	  } else if (this->comp == FWD_GCAG_INTRON_COMP || this->comp == FWD_ATAC_INTRON_COMP) {
	    pathscores[querypos] = SEMICANONICAL_POINTS;
	  } else if (this->comp == NONINTRON_COMP) {
	    pathscores[querypos] = 0; /* noncanonical */
	  }
	} else if (cdna_direction < 0) {
	  if (this->comp == REV_CANONICAL_INTRON_COMP) {
	    pathscores[querypos] = CANONICAL_POINTS;
	  } else if (this->comp == REV_GCAG_INTRON_COMP || this->comp == REV_ATAC_INTRON_COMP) {
	    pathscores[querypos] = SEMICANONICAL_POINTS;
	  } else if (this->comp == NONINTRON_COMP) {
	    pathscores[querypos] = 0; /* noncanonical */
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
      } else if (unknown_base(this->cdna) || unknown_base(this->genome)) {
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
  donor[0] = complCode[acceptor[1]];
  acceptor[1] = complCode[temp];

  temp = donor[1];
  donor[1] = complCode[acceptor[0]];
  acceptor[0] = complCode[temp];
  
  return;
}


void
Pair_print_protein_genomic (struct T *ptr, int npairs, int wraplength) {
  struct T *this;
  int xpos = 0, i;

  for (i = 0; i < npairs; i++) {
    this = ptr++;
    if (this->aa_g != ' ') {
      if (xpos == wraplength) {
	printf("\n");
	xpos = 0;
      }
#ifdef PMAP
      putchar(this->aa_g);
      xpos++;
#else
      if (this->aa_g != '*') {
	putchar(this->aa_g);
	xpos++;
      }
#endif
    }
  }
  printf("\n");
  return;
}

#ifdef PMAP
void
Pair_print_nucleotide_cdna (struct T *ptr, int npairs, int wraplength) {
  struct T *this;
  int xpos = 0, i;

  for (i = 0; i < npairs; i++) {
    this = ptr++;
    if (this->cdna != ' ') {
      if (xpos == wraplength) {
	printf("\n");
	xpos = 0;
      }
      putchar(this->cdna);
      xpos++;
    }
  }
  printf("\n");
  return;
}
#else
void
Pair_print_protein_cdna (struct T *ptr, int npairs, int wraplength) {
  struct T *this;
  int xpos = 0, i;

  for (i = 0; i < npairs; i++) {
    this = ptr++;
    if (this->aa_e != ' ') {
      if (xpos == wraplength) {
	printf("\n");
	xpos = 0;
      }
      if (this->aa_e != '*') {
	putchar(this->aa_e);
	xpos++;
      }
    }
  }
  printf("\n");
  return;
}
#endif


void
Pair_print_compressed (int pathnum, int npaths, T start, T end, Sequence_T queryseq, char *dbversion, 
		       int nexons, double fracidentity,
		       struct T *pairs, int npairs, Chrnum_T chrnum, Genomicpos_T chrpos,
		       Genomicpos_T chroffset, IIT_T chromosome_iit, 
		       int querylength_given, int skiplength, int trim_start, int trim_end,
		       Genomicpos_T genomiclength, bool checksump,
		       int chimerapos, int chimeraequivpos, double donor_prob, double acceptor_prob,
		       int chimera_cdna_direction, char *strain, bool watsonp, int cdna_direction,
		       bool zerobasedp, int worker_id) {
  Genomicpos_T chrpos1, chrpos2, position1, position2;

  bool in_exon = false;
  List_T tokens = NULL;
  struct T *ptr = pairs, *prev, *this = NULL;
  int querypos1, querypos2;
  int exon_querystart = -1, exon_genomestart = -1, exon_queryend, exon_genomeend,
    intron_start, intron_end;
  int num = 0, den = 0, runlength = 0, i;
  int print_dinucleotide_p;
  char token[10], donor[3], acceptor[3];
  double coverage, trimmed_coverage;

  donor[0] = donor[1] = donor[2] = '\0';
  acceptor[0] = acceptor[1] = acceptor[2] = '\0';

  querypos1 = start->querypos;
  querypos2 = end->querypos;

#ifdef PMAP
  printf(">%s %s %d/%d %d %d",
	 Sequence_accession(queryseq),dbversion,pathnum,npaths,
	 (querylength_given+skiplength)*3,nexons);
  coverage = (double) (querypos2 - querypos1 + 1)/(double) ((querylength_given+skiplength)*3);
  printf(" %.1f",((double) rint(1000.0*coverage)));
#else
  coverage = (double) (querypos2 - querypos1 + 1)/(double) (querylength_given+skiplength);
  if (end->querypos > trim_end) {
    trim_end = end->querypos;
  }
  if (start->querypos < trim_start) {
    trim_start = start->querypos;
  }
  trimmed_coverage = (double) (end->querypos - start->querypos + 1)/(double) (trim_end - trim_start + 1 + skiplength);
  printf(">%s %s %d/%d %d(%d) %d",
	 Sequence_accession(queryseq),dbversion,pathnum,npaths,
	 querylength_given+skiplength,trim_end-trim_start+1,nexons);
  printf(" %.1f(%.1f)",((double) rint(1000.0*coverage))/10.0,((double) rint(1000.0*trimmed_coverage))/10.0);
#endif
  printf(" %.1f",((double) rint(1000.0*fracidentity))/10.0);

  start = &(pairs[0]);
  end = &(pairs[npairs-1]);
  printf(" %d%s%d",start->querypos + !zerobasedp,"..",end->querypos + !zerobasedp);

  if (watsonp) {
    chrpos1 = chrpos + start->genomepos;
    chrpos2 = chrpos + end->genomepos;
  } else {
    chrpos1 = chrpos + (genomiclength - 1) - start->genomepos;
    chrpos2 = chrpos + (genomiclength - 1) - end->genomepos;
  }
  position1 = chroffset + chrpos1;
  position2 = chroffset + chrpos2;
  printf(" %u%s%u",position1 + !zerobasedp,"..",position2 + !zerobasedp);

  if (chrnum == 0) {
    printf(" %u%s%u",chrpos1 + !zerobasedp,"..",chrpos2 + !zerobasedp);
  } else {
    printf(" %s:%u%s%u",
	   Chrnum_to_string(chrnum,chromosome_iit,/*allocp*/false),
	   chrpos1 + !zerobasedp,"..",chrpos2 + !zerobasedp);
  }

  if (chrpos1 <= chrpos2) {
    printf(" +");
  } else {
    printf(" -");
  }

  if (cdna_direction > 0) {
    printf(" dir:sense");
  } else if (cdna_direction < 0) {
    printf(" dir:antisense");
  } else {
    printf(" dir:indet");
  }

  if (checksump == true) {
    printf(" md5:");
    Sequence_print_digest(queryseq);
  }

  if (chimerapos >= 0) {
    if (chimeraequivpos == chimerapos) {
      if (donor_prob > 0.0 && acceptor_prob > 0.0) {
	if (chimera_cdna_direction >= 0) {
	  printf(" chimera:%d(>)/%.3f/%.3f",chimerapos + !zerobasedp,donor_prob,acceptor_prob);
	} else {
	  printf(" chimera:%d(<)/%.3f/%.3f",chimerapos + !zerobasedp,donor_prob,acceptor_prob);
	}
      } else {
	printf(" chimera:%d",chimerapos + !zerobasedp);
      }
    } else {
      printf(" chimera:%d..%d",chimerapos + !zerobasedp,chimeraequivpos + !zerobasedp);
    }
  }

  if (strain != NULL) {
    printf(" strain:%s",strain);
  }

#if 0
  printf(" thread:%d",worker_id);
#endif

  printf("\n");

  for (i = 0; i < npairs; i++) {
    prev = this;
    this = ptr++;

    if (this->gapp) {
      if (in_exon == true) {
	/* Beginning of gap */
	exon_queryend = prev->querypos + !zerobasedp;
	if (watsonp) {
	  exon_genomeend = chrpos + prev->genomepos + !zerobasedp;
	  intron_start = exon_genomeend + 1;
	} else {
	  exon_genomeend = chrpos + (genomiclength - 1) - prev->genomepos + !zerobasedp;
	  intron_start = exon_genomeend - 1;
	}

	printf("\t%u %u",exon_genomestart,exon_genomeend);
	printf(" %d %d",exon_querystart,exon_queryend);
	if (den == 0) {
	  printf(" 100");
	} else {
	  printf(" %d",(int) floor(100.0*(double) num/(double) den));
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
	} else {
	  fprintf(stderr,"Can't parse %c in compression for %s\n",
		  this->comp,Sequence_accession(queryseq));
	  abort();
	}
	tokens = push_token(tokens,token);
	tokens = List_reverse(tokens);
	print_tokens_compressed(tokens);
	List_free(&tokens);
	printf("\t%d",exon_queryend - exon_querystart + 1);

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
	exon_querystart = this->querypos + !zerobasedp;
	if (watsonp) {
	  exon_genomestart = chrpos + this->genomepos + !zerobasedp;
	  intron_end = exon_genomestart - 1;
	} else {
	  exon_genomestart = chrpos + (genomiclength - 1) - this->genomepos + !zerobasedp;
	  intron_end = exon_genomestart + 1;
	}
	if (i > 0) {
	  if (intron_end > intron_start) {
	    printf("\t%d",intron_end - intron_start + 1);
	  } else {
	    printf("\t%d",intron_start - intron_end + 1);
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
	      printf("\t%c%c-%c%c",toupper(donor[0]),toupper(donor[1]),toupper(acceptor[0]),toupper(acceptor[1]));
	    }
	  }
#if 0
	  if (exon_querystart > exon_queryend + 1) {
	    printf("***");
	  }
#endif
	  printf("\n");
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
  }
  
  prev = this;
  exon_queryend = prev->querypos + !zerobasedp;
  if (watsonp) {
    exon_genomeend = chrpos + prev->genomepos + !zerobasedp;
  } else {
    exon_genomeend = chrpos + (genomiclength - 1) - prev->genomepos + !zerobasedp;
  }
  
  printf("\t%d %d",exon_genomestart,exon_genomeend);
  printf(" %d %d",exon_querystart,exon_queryend);
  if (den == 0) {
    printf(" 100");
  } else {
    printf(" %d",(int) floor(100.0*(double) num/(double) den));
  }

  sprintf(token,"%d*",runlength);
  tokens = push_token(tokens,token);
  tokens = List_reverse(tokens);
  print_tokens_compressed(tokens);
  List_free(&tokens);

  printf("\t%d",exon_queryend - exon_querystart + 1);
  printf("\n");

  return;
}


void
Pair_print_iit_map (Sequence_T queryseq, int pathnum, char *accession,
		    T start, T end, Chrnum_T chrnum, Genomicpos_T chrpos,
		    IIT_T chromosome_iit, Genomicpos_T genomiclength, bool watsonp) {
  char *chrstring = NULL;
  Genomicpos_T chrpos1, chrpos2;

  if (chrnum == 0) {
    chrstring = "";
  } else {
    chrstring = Chrnum_to_string(chrnum,chromosome_iit,/*allocp*/false);
  }

  if (watsonp) {
    chrpos1 = chrpos + start->genomepos;
    chrpos2 = chrpos + end->genomepos;
    printf(">%s %u %u %s+\n",accession,chrpos1,chrpos2,chrstring);
  } else {
    chrpos1 = chrpos + (genomiclength - 1) - start->genomepos;
    chrpos2 = chrpos + (genomiclength - 1) - end->genomepos;
    printf(">%s %u %u %s-\n",accession,chrpos1,chrpos2,chrstring);
  }
  Sequence_print_header(queryseq,/*checksump*/false);
  return;
}


