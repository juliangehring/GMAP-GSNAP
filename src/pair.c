static char rcsid[] = "$Id: pair.c,v 1.177 2006/04/07 01:17:41 twu Exp $";
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
Pair_flip (T this) {
  char c;
  int k;

  c = this->cdna;
  this->cdna = this->genome;
  this->genome = c;

  k = this->querypos;
  this->querypos = this->genomepos;
  this->genomepos = k;

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
#ifdef PMAP
    /* Make nucleotide consistent */
    if (this->comp == AMBIGUOUS_COMP) {
      putchar(this->genome);
    } else {
      putchar(this->cdna);
    }
#else
    putchar(this->cdna);
#endif
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
    if (diagnosticp == true && this->shortexonp == true) {
      putchar(DIAGNOSTIC_SHORTEXON_COMP);
    } else if (this->comp == MATCH_COMP) {
      putchar(MATCH_COMP);
    } else if (diagnosticp == false && this->comp == DYNPROG_MATCH_COMP) {
      putchar(MATCH_COMP);
    } else if (diagnosticp == false && this->comp == AMBIGUOUS_COMP) {
      putchar(MATCH_COMP);
    } else if (this->comp == SHORTGAP_COMP) {
      putchar(INDEL_COMP);
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
    putchar(this->genome);
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

static struct T *
invert_path (struct T *old, int npairs) {
  struct T *new;
  char complCode[128] = COMPLEMENT;
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
  char complCode[128] = COMPLEMENT;
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


static void
add_intronlengths (struct T *pairs, int npairs, int ngap) {
  struct T *prev, *this = NULL, *ptr;
  int space, margin, i, j, k, gapstart, genomepos;
  char intronstring[20], cdnabreak[20], genomicbreak[20], comp;

  i = 0;
  while (i < npairs) {
    prev = this;
    this = &(pairs[i++]);

    if (this->gapp) {
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

#if 0
	/* Fix genomic positions on left end of gap */
	genomepos = prev->genomepos + 1;
	for (j = gapstart; j < gapstart + ngap; j++) {
	  ptr = &(pairs[j]);
	  ptr->genomepos = genomepos++;
	}

	/* Fix genomic positions on right end of gap */
	genomepos = this->genomepos - 1;
	for (j = gapstart + space - 1; j >= gapstart + space - ngap; --j) {
	  ptr = &(pairs[j]);
	  ptr->genomepos = genomepos--;
	}

	/* Fix genomic positions in middle of gap */
	for ( ; j >= gapstart + space - ngap - 3; --j) {
	  ptr = &(pairs[j]);
	  ptr->genomepos = genomepos;
	}
#endif

      }
    }
  }	
  return;
}


void
Pair_print_continuous (struct T *pairs, int npairs, Chrnum_T chrnum, Genomicpos_T chrpos,
		       Genomicpos_T chroffset, Genomicpos_T genomiclength,
		       bool watsonp, int cdna_direction, bool universalp, bool zerobasedp, bool diagnosticp, 
		       bool genomefirstp, int invertmode, bool nointronlenp, int ngap) {
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
    add_intronlengths(ptr,npairs,ngap);
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
		      bool genomicprop, int invertmode, bool nointronlenp, int wraplength, int ngap) {
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
    add_intronlengths(ptr,npairs,ngap);
  }
  if (chrnum != 0) {
    if (invertmode == 2) {
      chrstring = Chrnum_to_string(chrnum,chromosome_iit);
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
Pair_debug_alignment (List_T list, int ngap) {
  struct T *pairs;
  int npairs;

  npairs = List_length(list);
  pairs = Pair_block_copy(list,&npairs,ngap);
  Pair_print_alignment(pairs,npairs,/*chrnum*/0,/*chrpos*/0U,/*chroffset*/0U,/*chromosome_iit*/NULL,
		       /*genomiclength*/0,/*watsonp*/true,/*cdna_direction*/0,/*universalp*/false,
		       /*zerobasedp*/false,/*diagnosticp*/true,/*genomicprop*/true,
		       /*invertmode*/0,/*nointronlenp*/false,/*wraplength*/50,ngap);
  FREE(pairs);
  return;
}

void
Pair_print_pathsummary (int pathnum, T start, T end, Chrnum_T chrnum, Genomicpos_T chrpos,
			Genomicpos_T chroffset, IIT_T chromosome_iit, bool referencealignp, 
			IIT_T altstrain_iit, char *strain, 
			IIT_T contig_iit, char *dbversion, Genomicpos_T genomiclength,
			int nexons, double coverage, int matches, int unknowns, int mismatches, 
			int qopens, int qindels, int topens, int tindels, int goodness,
			bool watsonp, int cdna_direction, double defect_rate, 
			int translation_start, int translation_end, int translation_length,
			int relaastart, int relaaend, bool zerobasedp, bool maponlyp) {
  int querypos1, querypos2, chrpos1, chrpos2, den;
  double fracidentity;
  Genomicpos_T position1, position2;
  char *chrstring = NULL, *refstrain, *comma1, *comma2;

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
    chrstring = Chrnum_to_string(chrnum,chromosome_iit);
    if (watsonp) {
      printf("chr %s:%s%s%s (%d bp)\n",
	     chrstring,comma1,SEPARATOR,comma2,chrpos2-chrpos1+1);
    } else {
      printf("chr %s:%s%s%s (%d bp)\n",
	     chrstring,comma1,SEPARATOR,comma2,chrpos2-chrpos1-1);
    }
  }
  FREE(comma2);
  FREE(comma1);

  if (maponlyp == false) {
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
    printf("    Coverage: %.1f",((double) rint(1000.0*coverage))/10.0);
    printf("\n");

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
      printf("    Translation: %d..%d (%d aa)\n",
	     translation_start+!zerobasedp,translation_end,translation_length);
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

  if (chrstring != NULL) {
    FREE(chrstring);
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
      chrstring = Chrnum_to_string(chrnum,chromosome_iit);
    } else {
      chrstring = Chrnum_to_string_signed(chrnum,chromosome_iit,watsonp);
    }
  }

  if (watsonp) {
    for (i = 0; i < npairs; i++) {
      this = pairs++;
      if (this->gapp == false) {
	if (this->aamarker_e == false) {
	  printf("%d\t",this->aapos);
	} else {
	  printf("%d %c\t",this->aapos,this->aa_e);
	}
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
	if (this->aamarker_g == false) {
	  printf("\t");
	} else {
	  printf("\t%c",this->aa_g);
	}
	printf("\n");
      }
    }
  } else {
    for (i = 0; i < npairs; i++) {
      this = pairs++;
      if (this->gapp == false) {
	if (this->aamarker_e == false) {
	  printf("%d\t",this->aapos);
	} else {
	  printf("%d %c\t",this->aapos,this->aa_e);
	}
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
	if (this->aamarker_g == false) {
	  printf("\t");
	} else {
	  printf("\t%c",this->aa_g);
	}
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
  printf("%d %d %c %c %c\n",
	 this->querypos + !zerobasedp,this->genomepos + !zerobasedp,
	 this->cdna,this->comp,this->genome);
}


/* Useful for debugging */
void
Pair_dump_list (List_T pairs, bool zerobasedp) {
  T this;
  List_T p;

  for (p = pairs; p != NULL; p = List_next(p)) {
    this = List_head(p);
    printf("%d %d %c %c %c\n",
	   this->querypos + !zerobasedp,this->genomepos + !zerobasedp,
	   this->cdna,this->comp,this->genome);
  }
  return;
}  

void
Pair_dump_array (struct T *pairs, int npairs, bool zerobasedp) {
  struct T *this;
  int i;

  for (i = 0; i < npairs; i++) {
    this = pairs++;
    printf("%d: %d %d %d %c %c %c",
	   i,this->querypos + !zerobasedp,this->genomepos + !zerobasedp,this->aapos,
	   this->cdna,this->comp,this->genome);
    if (this->aamarker_g == true || this->aamarker_e == true) {
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


void
Pair_check_array (struct T *pairs, int npairs) {
  struct T *this;
  int prev_querypos, prev_genomepos;
  int i;

  if (npairs == 0) {
    return;
  } else {
    this = pairs++;
    prev_querypos = this->querypos;
    prev_genomepos = this->genomepos;

    for (i = 1; i < npairs; i++) {
      this = pairs++;
      if (this->querypos < prev_querypos) {
	fprintf(stderr,"Problem at querypos %d\n",this->querypos);
      }
      if (this->genomepos < prev_genomepos) {
	fprintf(stderr,"Problem at genomepos %d\n",this->genomepos);
      }
      prev_querypos = this->querypos;
      prev_genomepos = this->genomepos;
    }
  }
  return;
}  


static bool
unknown_base (char c) {
  switch (c) {
  case 'A': case 'C': case 'G': case 'T': 
  case 'a': case 'c': case 'g': case 't': return false;
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
      chrstring = Chrnum_to_string_signed(chrnum,chromosome_iit,true);
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
	  printf("\n");
	}
	num = den = 0;
	in_exon = true;
      }
      if (this->comp == INDEL_COMP || this->comp == SHORTGAP_COMP) {
	/* Don't count in numerator or denominator */
#ifndef PMAP
      } else if (unknown_base(this->cdna) || unknown_base(this->genome)) {
	/* Comp must be a space */
	/* Don't count in numerator or denominator */
#endif
      } else {
	den++;
	if (this->comp == MATCH_COMP || this->comp == DYNPROG_MATCH_COMP
#ifdef PMAP
	    || this->comp == AMBIGUOUS_COMP
#endif
	    ) {
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
		      int npairs, int querylength, Genomicpos_T chrpos, Genomicpos_T genomiclength, bool watsonp) {
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
	  *tStarts = Uintlist_push(*tStarts,chrpos + (genomiclength - 1) - this->genomepos);
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
compute_gap_lengths (int *qnbreaks, int *qlength, int *tnbreaks, int *tlength, struct T *pairs, int npairs) {
  struct T *prev, *this = NULL;
  int i;
  char comp;

  *qnbreaks = *qlength = *tnbreaks = *tlength = 0;

  i = 0;
  while (i < npairs) {
    prev = this;
    this = &(pairs[i++]);

    if (this->gapp) {
      comp = this->comp;
      while (i < npairs && this->gapp) {
	this = &(pairs[i++]);
      }

      if (comp == DUALBREAK_COMP) {
	*qnbreaks += 1;
	*tnbreaks += 1;
	*qlength += abs(this->querypos - prev->querypos)-1;
	*tlength += abs(this->genomepos - prev->genomepos)-1;

      } else {			/* Intron */
	*tnbreaks += 1;
	*tlength += abs(this->genomepos - prev->genomepos)-1;
      }
    }
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
  char *chrstring;
  struct T *pairs_directional = NULL;
  Intlist_T blockSizes = NULL, qStarts = NULL, p;
  Uintlist_T tStarts = NULL, q;
  int nblocks;
  int qnbreaks, qlength, tnbreaks, tlength, querylength;

  compute_gap_lengths(&qnbreaks,&qlength,&tnbreaks,&tlength,pairs,npairs);

  printf("%d\t%d\t%d\t%d\t",matches,mismatches,/*repeatmatches*/0,unknowns);
  printf("%d\t%d\t%d\t%d\t",qopens+qnbreaks,qindels+qlength,topens+tnbreaks,tindels+tlength);

  if (watsonp == true) {
    printf("+");
  } else {
    printf("-");
  }
  printf("\t%s\t%d",Sequence_accession(queryseq),Sequence_fulllength(queryseq));

  printf("\t%d\t%d",start->querypos,end->querypos+1);

  /* T name and T size */
  chrstring = Chrnum_to_string(chrnum,chromosome_iit);
  printf("\t%s\t%u",chrstring,Chrnum_length(chrnum,chromosome_iit));
  FREE(chrstring);

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

  if (watsonp == true) {
    pairs_directional = pairs;
  } else {
    pairs_directional = invert_and_revcomp_path(pairs,npairs);
#ifdef PMAP
    querylength = 3*Sequence_fulllength(queryseq);
#else
    querylength = Sequence_fulllength(queryseq);
#endif
  }

  nblocks = count_psl_blocks_nt(&blockSizes,&qStarts,&tStarts,pairs_directional,npairs,
				querylength,chrpos,genomiclength,watsonp);
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
  Genomicpos_T chrpos1, chrpos2;
  char *chrstring;
  Intlist_T blockSizes = NULL, qStarts = NULL, p;
  Uintlist_T tStarts = NULL, q;
  int nblocks, matches = 0, mismatches = 0, unknowns = 0;
  int qnbreaks, qlength, tnbreaks, tlength, querylength;

  compute_gap_lengths(&qnbreaks,&qlength,&tnbreaks,&tlength,pairs,npairs);
  count_matches_pro(&matches,&mismatches,&unknowns,pairs,npairs);

  printf("%d\t%d\t%d\t%d\t",matches,mismatches,/*repeatmatches*/0,unknowns);
  printf("%d\t%d\t%d\t%d\t",qopens+qnbreaks,qindels+qlength,topens+tnbreaks,tindels+tlength);

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
  printf("\t%s\t%d",Sequence_accession(queryseq),Sequence_fulllength(queryseq));

  printf("\t%d\t%d",(start->querypos+2)/3,end->querypos/3+1);

  /* T name and T size */
  chrstring = Chrnum_to_string(chrnum,chromosome_iit);
  printf("\tchr%s\t%u",chrstring,Chrnum_length(chrnum,chromosome_iit));
  FREE(chrstring);

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

  querylength = 3*Sequence_fulllength(queryseq);

  nblocks = count_psl_blocks_pro(&blockSizes,&qStarts,&tStarts,pairs,npairs,
				 querylength,chrpos,genomiclength,watsonp);
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
      if (this->comp == INDEL_COMP) {
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
      } else if (this->comp == SHORTGAP_COMP) {
	/* Do nothing */
#ifndef PMAP
      } else if (unknown_base(this->cdna) || unknown_base(this->genome) || this->comp == AMBIGUOUS_COMP) {
	(*unknowns)++;
#endif
      } else if (this->comp == MATCH_COMP || this->comp == DYNPROG_MATCH_COMP
#ifdef PMAP
		 || this->comp == AMBIGUOUS_COMP
#endif
		 ) {
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
	if (this->comp == INDEL_COMP) {
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
	} else if (this->comp == SHORTGAP_COMP) {
	  /* Do nothing */
#ifndef PMAP
	} else if (unknown_base(this->cdna) || unknown_base(this->genome) || this->comp == AMBIGUOUS_COMP) {
	  (*unknowns)++;
#endif
	} else if (this->comp == MATCH_COMP || this->comp == DYNPROG_MATCH_COMP
#ifdef PMAP
		   || this->comp == AMBIGUOUS_COMP
#endif
		   ) {
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
      if (this->comp != MATCH_COMP && this->comp != DYNPROG_MATCH_COMP
#ifdef PMAP
	  && this->comp != AMBIGUOUS_COMP
#endif
	  ) {
	matchscores[querypos] = 0; /* For mismatch */
      }
    }
  }

  return matchscores;
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
      if (this->comp == INDEL_COMP) {
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
      } else if (this->comp == SHORTGAP_COMP) {
	/* Do nothing */
#ifndef PMAP
      } else if (unknown_base(this->cdna) || unknown_base(this->genome) || this->comp == AMBIGUOUS_COMP) {
	/* (*unknowns)++; */
#endif
      } else if (this->comp == MATCH_COMP || this->comp == DYNPROG_MATCH_COMP
		 || this->comp == AMBIGUOUS_COMP
		 ) {
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
Pair_nexons (List_T pairs) {
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


/* TODO: Generalize to use alt canonical */
static void
merge_one_gap (struct T *pairs, int gapcount, int ngap) {
  int j, k;
  char newmatch;

  j = -gapcount;
  if (Intron_canonical_fwd_p(pairs[j].genome,pairs[j+1].genome,pairs[-2].genome,pairs[-1].genome)) {
    newmatch = FWD_CANONICAL_INTRON_COMP;
  } else if (Intron_canonical_rev_p(pairs[j].genome,pairs[j+1].genome,pairs[-2].genome,pairs[-1].genome)) {
    newmatch = REV_CANONICAL_INTRON_COMP;
  } else if (Intron_gcag_fwd_p(pairs[j].genome,pairs[j+1].genome,pairs[-2].genome,pairs[-1].genome)) {
    newmatch = FWD_GCAG_INTRON_COMP;
  } else if (Intron_gcag_rev_p(pairs[j].genome,pairs[j+1].genome,pairs[-2].genome,pairs[-1].genome)) {
    newmatch = REV_GCAG_INTRON_COMP;
  } else if (Intron_atac_fwd_p(pairs[j].genome,pairs[j+1].genome,pairs[-2].genome,pairs[-1].genome)) {
    newmatch = FWD_ATAC_INTRON_COMP;
  } else if (Intron_atac_rev_p(pairs[j].genome,pairs[j+1].genome,pairs[-2].genome,pairs[-1].genome)) {
    newmatch = REV_ATAC_INTRON_COMP;
  } else {
    newmatch = NONINTRON_COMP;
  }

  for (k = 0; k < ngap; k++, j++) {
    pairs[j].comp = newmatch;
  }
  j += 3;			/* For ... */
  for ( ; j < -ngap; j++) {
    pairs[j].comp = 'D';
  }
  for ( ; j < 0; j++) {
    pairs[j].comp = newmatch;
  }
  return;
}


static int
merge_gaps (struct T *pairs, int npairs, int ngap) {
  struct T *this;
  int i;
  int nmerged = 0, gapcount = 0;
  int gaplength = ngap + 3 + ngap;

  for (i = 0; i < npairs; i++) {
    this = pairs++;
    if (this->gapp) {
      gapcount++;
    } else {
      if (gapcount > gaplength) {
	merge_one_gap(pairs-1,gapcount,ngap);
	nmerged += gapcount - gaplength;
      }
      gapcount = 0;
    }
  }
  if (gapcount > gaplength) {
    merge_one_gap(pairs-1,gapcount,ngap);
    nmerged += gapcount - gaplength;
  }

  return nmerged;
}


static int
count_tail_gap (List_T list) {
  int gaplength = 0;
  T this;
  List_T p;

  for (p = list; p != NULL; p = p->rest) {
    this = p->first;
    if (this->gapp) {
      gaplength++;
    } else {
      gaplength = 0;
    }
  }
  return gaplength;
}


/* Given a list of pairs, returns a block of pairs.  Removes all
   occurences of comp == 'D' introduced by merge_gaps */
struct T *
Pair_block_copy (List_T list, int *npairs, int ngap) {
  struct T *chunk, *chunk2;
  T newpair, oldpair;
  List_T p;
  int nmerged, tailgap, i, j;

  if ((tailgap = count_tail_gap(list)) > 0) {
    *npairs -= tailgap;
  }

  newpair = chunk = (struct T *) MALLOC((*npairs)*sizeof(struct T));
  for (p = list, i = 0; p != NULL && i < *npairs; p = p->rest, i++) {
    oldpair = (T) p->first;
    /*
    newpair->querypos = oldpair->querypos;
    newpair->genomepos = oldpair->genomepos;
    newpair->cdna = oldpair->cdna;
    newpair->comp = oldpair->comp;
    newpair->genome = oldpair->genome;
    newpair->gapp = oldpair->gapp;
    newpair->shortexonp = oldpair->shortexonp;
    */
    memcpy(newpair,oldpair,sizeof(struct T));
    newpair++;
  }

  /* If ngap < 0, it means we are sure there are no gaps to merge, such as
     when we are copying a Stage3_T object. */
  if (ngap < 0 || (nmerged = merge_gaps(chunk,*npairs,ngap)) == 0) {
    return chunk;
  } else {
    chunk2 = (struct T *) MALLOC((*npairs-nmerged)*sizeof(struct T));
    j = 0;
    for (i = 0; i < *npairs; i++) {
      if (chunk[i].comp != 'D') {
	/*
	chunk2[j].querypos = chunk[i].querypos;
	chunk2[j].genomepos = chunk[i].genomepos;
	chunk2[j].cdna = chunk[i].cdna;
	chunk2[j].comp = chunk[i].comp;
	chunk2[j].genome = chunk[i].genome;
	chunk2[j].gapp = chunk[i].gapp;
	chunk2[j].shortexonp = chunk[i].shortexonp;
	*/
	memcpy(&(chunk2[j]),&(chunk[i]),sizeof(struct T));
	j++;
      }
    }
    FREE(chunk);
    *npairs -= nmerged;
    return chunk2;
  }
}


static void
print_tokens (List_T tokens) {
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

static List_T
push_token (List_T tokens, char *token) {
  char *copy;

  copy = (char *) CALLOC(strlen(token)+1,sizeof(char));
  strcpy(copy,token);
  return List_push(tokens,(void *) copy);
}

static void
invert_intron (char *donor, char *acceptor) {
  char temp, complCode[128] = COMPLEMENT;

  temp = donor[0];
  donor[0] = complCode[(int) acceptor[1]];
  acceptor[1] = complCode[(int) temp];

  temp = donor[1];
  donor[1] = complCode[(int) acceptor[0]];
  acceptor[0] = complCode[(int) temp];
  
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
Pair_print_compressed (Sequence_T queryseq, char *version, int pathnum, int npaths,
		       int nexons, double coverage, double fracidentity,
		       struct T *pairs, int npairs, Chrnum_T chrnum, Genomicpos_T chrpos,
		       Genomicpos_T chroffset, IIT_T chromosome_iit, 
		       Genomicpos_T genomiclength, bool checksump,
		       int chimerapos, int chimeraequivpos, double donor_prob, double acceptor_prob,
		       int chimera_cdna_direction, char *strain, bool watsonp, int cdna_direction,
		       bool zerobasedp, int worker_id) {
  Genomicpos_T chrpos1, chrpos2, position1, position2;
  char *chrstring = NULL;

  bool in_exon = false;
  List_T tokens = NULL;
  struct T *ptr = pairs, *prev, *this = NULL;
  T start, end;
  int exon_querystart = -1, exon_genomestart = -1, exon_queryend, exon_genomeend,
    intron_start, intron_end;
  int num = 0, den = 0, runlength = 0, i;
  int print_dinucleotide_p;
  char token[10], donor[3], acceptor[3];

  donor[0] = donor[1] = donor[2] = '\0';
  acceptor[0] = acceptor[1] = acceptor[2] = '\0';

#ifdef PMAP
  printf(">%s %s %d/%d %d %d",
	 Sequence_accession(queryseq),version,pathnum,npaths,
	 (Sequence_fulllength(queryseq)+Sequence_skiplength(queryseq))*3,nexons);
#else
  printf(">%s %s %d/%d %d %d",
	 Sequence_accession(queryseq),version,pathnum,npaths,
	 Sequence_fulllength(queryseq)+Sequence_skiplength(queryseq),nexons);
#endif
  printf(" %.1f",((double) rint(1000.0*coverage))/10.0);
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
    chrstring = Chrnum_to_string(chrnum,chromosome_iit);
    printf(" %s:%u%s%u",chrstring,chrpos1 + !zerobasedp,"..",chrpos2 + !zerobasedp);
    FREE(chrstring);
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
	print_tokens(tokens);
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
	    if (!strcmp(donor,"GT") && !strcmp(acceptor,"AG")) {
	      /* Do nothing */
	    } else {
	      printf("\t%s-%s",donor,acceptor);
	    }
	  }
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
#endif

      } else {
	runlength++;
	den++;
	if (this->comp == MATCH_COMP || this->comp == DYNPROG_MATCH_COMP
#ifdef PMAP
	    || this->comp == AMBIGUOUS_COMP
#endif
	    ) {
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
  print_tokens(tokens);
  List_free(&tokens);

  printf("\t%d",exon_queryend - exon_querystart + 1);
  printf("\n");

  return;
}


