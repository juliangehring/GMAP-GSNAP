static char rcsid[] = "$Id: pair.c,v 1.122 2005/03/09 19:23:05 twu Exp $";
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
#include <math.h>		/* For rint() */
#include "mem.h"
#include "complement.h"
#include "intron.h"
#include "listdef.h"
#include "intlist.h"
#include "separator.h"
#include "scores.h"
#include "translation.h"


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

void
Pair_free (T *old) {
  if (*old) {
    FREE(*old);
  }
  return;
}



/* Print routines */

static char *
truncate_left (char *name, int length) {
  char *ptr;

  if (strlen(name) < length) {
    return name;
  } else {
    ptr = name;
    while (*ptr != '\0') {
      ++ptr;
    }
    while (length-- > 0) {
      --ptr;
    }
  }
  return ptr;
}


static char *RULER = "    .    :    .    :    .    :    .    :    .    :";
static void
print_top_ruler (int n, int npairs, int wraplength) {
  printf("%14d ",n);
  if (n + wraplength < npairs) {
    printf("%s\n",RULER);
  } else {
    printf("%.*s\n",npairs-n,RULER);
  }
  return;
}

static void
print_bottom_ruler (int n, int npairs, int wraplength) {
  printf("%14s ","");
  if (n + wraplength < npairs) {
    printf("%s\n",RULER);
  } else {
    printf("%.*s\n",npairs-n,RULER);
  }
  return;
}


static void
print_cdna_sequence (struct T *ptr, int n, int npairs, bool zerobasedp,
		     int wraplength) {
  struct T *this;
  int i;

  this = ptr;
  printf("%14u ",this->querypos + !zerobasedp);
  for (i = 0; n < npairs && i < wraplength; n++, i++) {
    this = ptr++;
    printf("%c",this->cdna);
  }
  printf("\n");
  return;
}

static void
print_peptide (struct T *ptr, int n, int npairs, int wraplength, bool genomep) {
  struct T *this, *last = NULL;
  int aapos, i;
  bool proteinp = false;

  if (ptr->aapos > 0) {
    proteinp = true;
  } else {
    /* Check if protein starts in middle of this line */
    if (npairs - n < wraplength) {
      last = &ptr[npairs - n - 1];
    } else {
      last = &ptr[wraplength - 1];
    }
    if (last->aapos != 0) {
      proteinp = true;
    }
  }

  if (proteinp == false) {
    printf("%14s ","");
  } else {
    /* Look for first aa to determine aapos */
    this = ptr;
    if (last == NULL) {
      if (npairs - n < wraplength) {
	last = &ptr[npairs - n - 1];
      } else {
	last = &ptr[wraplength - 1];
      }
    }
    while (this <= last && (genomep ? this->aa_g : this->aa_e) == ' ') {
      this++;
    }

    if (this > last) {
      aapos = ptr->aapos;
    } else {
      aapos = this->aapos;
    }
    if (genomep == true) {
      printf("aa.g%10d ",aapos);
    } else {
      printf("aa.c%10d ",aapos);
    }
  }

  for (i = 0; n < npairs && i < wraplength; n++, i++) {
    this = ptr++;
    if (genomep == true) {
      printf("%c",this->aa_g);
    } else {
      printf("%c",this->aa_e);
    }
  }
  printf("\n");
  return;
}

static void
print_alignment (struct T *ptr, int n, int npairs, bool diagnosticp, int wraplength) {
  struct T *this;
  int i;

  printf("%14s ","");
  for (i = 0; n < npairs && i < wraplength; n++, i++) {
    this = ptr++;
    if (diagnosticp == true && this->shortexonp == true) {
      printf("s");
    } else if (diagnosticp == false && this->comp == '*') {
      printf("|");
    } else if (this->comp == '~') {
      printf("-");
    } else {
      printf("%c",this->comp);
    }
  }
  printf("\n");
  return;
}


static void
print_genomic_sequence (struct T *ptr, int n, int npairs, bool zerobasedp,
			char *chrstring, Genomicpos_T chrpos, Genomicpos_T chroffset,
			Genomicpos_T genomiclength, bool watsonp, bool universalp,
			int wraplength) {
  struct T *this;
  int i;
  char Buffer[14];

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
  printf("%14s ",Buffer);
  for (i = 0; n < npairs && i < wraplength; n++, i++) {
    this = ptr++;
    printf("%c",this->genome);
  }
  printf("\n");
  return;
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
    new[j].comp = complCode[old[i].comp];
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
    new[j].cdna = complCode[old[i].cdna];
    new[j].genome = complCode[old[i].genome];
    new[j].comp = complCode[old[i].comp];
  }
  return new;
}

static char *
reverse_string (char *old, int nchars) {
  char *new;
  int i, j;

  new = (char *) CALLOC(nchars+1,sizeof(char));
  for (i = 0, j = nchars-1; i < nchars; i++, j--) {
    new[j] = old[i];
  }
  return new;
}


static void
add_intronlengths (struct T *pairs, int npairs) {
  struct T *prev, *this = NULL, *ptr;
  int space, margin, i, j, k, gapstart;
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

      if (comp == '#') {	/* Dual break */
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
  int n = 0, i;

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
      printf("%c",this->genome);
    }
    printf("\n");

    ptr = pairs;
    for (n = 0; n < npairs; n++) {
      this = ptr++;
      if (diagnosticp == false && this->comp == '*') {
	printf("|");
      } else {
	printf("%c",this->comp);
      }
    }
    printf("\n");

    ptr = pairs;
    for (n = 0; n < npairs; n++) {
      this = ptr++;
      printf("%c",this->cdna);
    }
    printf("\n");

  } else {
    ptr = pairs;
    for (n = 0; n < npairs; n++) {
      this = ptr++;
      printf("%c",this->cdna);
    }
    printf("\n");

    ptr = pairs;
    for (n = 0; n < npairs; n++) {
      this = ptr++;
      if (diagnosticp == false && this->comp == '*') {
	printf("|");
      } else {
	printf("%c",this->comp);
      }
    }
    printf("\n");

    ptr = pairs;
    for (n = 0; n < npairs; n++) {
      this = ptr++;
      printf("%c",this->genome);
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

  while (n < npairs) {
    print_top_ruler(n,npairs,wraplength);
    if (genomicprop == true) {
      print_peptide(ptr,n,npairs,wraplength,/*genomep*/true);
      print_genomic_sequence(ptr,n,npairs,zerobasedp,chrstring,chrpos,chroffset,genomiclength,
			     watsonp,universalp,wraplength);
      print_alignment(ptr,n,npairs,diagnosticp,wraplength);
      print_cdna_sequence(ptr,n,npairs,zerobasedp,wraplength);
      print_peptide(ptr,n,npairs,wraplength,/*genomep*/false);
    } else {
      print_peptide(ptr,n,npairs,wraplength,/*genomep*/true);
      print_genomic_sequence(ptr,n,npairs,zerobasedp,chrstring,chrpos,chroffset,genomiclength,
			     watsonp,universalp,wraplength);
      print_alignment(ptr,n,npairs,diagnosticp,wraplength);
      print_cdna_sequence(ptr,n,npairs,zerobasedp,wraplength);
      print_peptide(ptr,n,npairs,wraplength,/*genomep*/false);
    }
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
		       /*invertmode*/0,/*nointronlenp*/false,/*wraplength*/50);
  FREE(pairs);
  return;
}

void
Pair_print_pathsummary (int pathnum, T start, T end, Chrnum_T chrnum, Genomicpos_T chrpos,
			Genomicpos_T chroffset, IIT_T chromosome_iit, bool referencealignp, char *strain, 
			IIT_T contig_iit, char *dbversion, Genomicpos_T genomiclength,
			int nexons, double coverage, int coverage_correction, int ntrimmed, 
			int matches, int unknowns, int mismatches, 
			int qopens, int qindels, int topens, int tindels, int goodness,
			bool watsonp, int cdna_direction, double defect_rate, 
			int translation_start, int translation_end, int translation_length,
			int relaastart, int relaaend, bool zerobasedp) {
  int querypos1, querypos2, chrpos1, chrpos2, den;
  double fracidentity;
  Genomicpos_T position1, position2;
  char *chrstring = NULL, *comma1, *comma2;

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

  printf("    cDNA direction: ");
  if (cdna_direction > 0) {
    printf("sense\n");
  } else if (cdna_direction < 0) {
    printf("antisense\n");
  } else {
    printf("indeterminate\n");
  }

  if (strain != NULL) {
    if (strain[0] == '\0') {
      /* Backward compatibility with old altstrain_iit */
      printf("    Strain: reference\n",strain);
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
    
  printf("    Number of exons: %d\n",nexons);
  printf("    Coverage: %.1f",((double) rint(1000.0*coverage))/10.0);
  if (ntrimmed > 0 && coverage_correction > 0) {
    printf(" (%d trimmed, %d genomic gap)",ntrimmed,coverage_correction);
  } else if (ntrimmed > 0) {
    printf(" (%d trimmed)",ntrimmed);
  } else if (coverage_correction > 0) {
    printf(" (%d genomic gap)",coverage_correction);
  }
  printf("\n");

  if ((den = matches + mismatches + qindels + tindels) == 0) {
    fracidentity = 1.0;
  } else {
    fracidentity = (double) matches/(double) den;
  }

  debug(printf("     Goodness: %d\n",goodness));

  printf("    Percent identity: %.1f (%d matches, %d mismatches, %d indels, %d unknowns)\n",
	 ((double) rint(1000.0*fracidentity))/10.0,matches,mismatches,qindels+tindels,unknowns);
  if (qindels + tindels > 0) {
    printf("    Non-intron gaps: %d openings, %d bases in cdna; %d openings, %d bases in genome\n",
	   qopens,qindels,topens,tindels);
  } 

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

  /* printf("    Defect rate (percent): %.1f\n",defect_rate*100.0); */

  /* printf("\n"); -- Done by caller */

  if (chrstring != NULL) {
    FREE(chrstring);
  }

  return;
}


static int
find_genomicpos (Genomicpos_T position, struct T *ptr, int npairs, Genomicpos_T chrpos,
		 Genomicpos_T chroffset, Genomicpos_T genomiclength, bool watsonp, bool zerobasedp) {
  struct T *this;
  Genomicpos_T thisposition;
  int i;

  this = ptr;
  for (i = 0; i < npairs; i++) {
    if (watsonp) {
      thisposition = chroffset+chrpos+this->genomepos + !zerobasedp;
    } else {
      thisposition = chroffset+chrpos + (genomiclength - 1) - this->genomepos + !zerobasedp;
    }
    if (thisposition == position) {
      return i;
    }
    this++;
  }

  return -1;
}

static int
find_chrpos (Genomicpos_T position, struct T *pairs, int npairs, 
	     int translation_start, int translation_end, Genomicpos_T chrpos,
	     Genomicpos_T genomiclength, bool watsonp, bool zerobasedp) {
  struct T *this, *start, *end;
  Genomicpos_T thisposition, startposition, endposition;
  int i;

  start = &(pairs[0]);
  end = &(pairs[npairs-1]);

  /* Check endpoints */
  if (watsonp) {
    startposition = chrpos+start->genomepos + !zerobasedp;
    endposition = chrpos+end->genomepos + !zerobasedp;
    if (position < startposition || position > endposition) {
      return -9;
    }
  } else {
    startposition = chrpos + (genomiclength - 1) - start->genomepos + !zerobasedp;
    endposition = chrpos + (genomiclength - 1) - end->genomepos + !zerobasedp;
    if (position < endposition || position > startposition) {
      return -9;
    }
  }

  /* Check 5' UTR */
  this = &(pairs[translation_start]);
  if (watsonp) {
    thisposition = chrpos+this->genomepos + !zerobasedp;
    if (position < thisposition) {
      return -5;
    }
  } else {
    thisposition = chrpos + (genomiclength - 1) - this->genomepos + !zerobasedp;
    if (position > thisposition) {
      return -5;
    }
  }

  /* Check 3' UTR */
  this = &(pairs[translation_end]);
  if (watsonp) {
    thisposition = chrpos+this->genomepos + !zerobasedp;
    if (position > thisposition) {
      return -3;
    }
  } else {
    thisposition = chrpos + (genomiclength - 1) - this->genomepos + !zerobasedp;
    if (position < thisposition) {
      return -3;
    }
  }

  /* Check all */
  this = pairs;
  for (i = 0; i < npairs; i++) {
    if (watsonp) {
      thisposition = chrpos+this->genomepos + !zerobasedp;
    } else {
      thisposition = chrpos + (genomiclength - 1) - this->genomepos + !zerobasedp;
    }
    if (thisposition == position) {
      if (this->gapp == true) {
	return -1;
      } else {
	return i;
      }
    }
    this++;
  }

  return -1;
}

#if 0
extern void
Pair_print_mutation (struct T *pairs, int npairs, int cdna_direction,
		     int translation_frame, int translation_start, int translation_end,
		     Genomicpos_T chrpos, Genomicpos_T genomiclength, bool watsonp, 
		     bool zerobasedp, Genomicpos_T mutposition, char *change) {
  int mutindex, codon_index;
  char newnt, codon_orig, codon_mut;

  mutindex = find_chrpos(mutposition,pairs,npairs,translation_start,translation_end,
			 chrpos,genomiclength,watsonp,zerobasedp);
  printf("%u%s\t",mutposition,change);
  if (mutindex == -9) {
    printf("Out of bounds");
  } else if (mutindex == -5) {
    printf("5' UTR");
  } else if (mutindex == -3) {
    printf("3' UTR");
  } else if (mutindex == -1) {
    printf("Intron");
  } else {
    if (change[0] == 'i' || change[0] == 'd' || change[0] == '^' || change[0] == 'v') {
      if (cdna_direction < 0) {
	Translation_mutation_effect(&codon_index,&codon_orig,&codon_mut,
				    pairs,npairs,translation_frame,translation_start,
				    translation_end,mutindex,/*newnt*/'X',
				    /*backwardsp*/true,/*revcompp*/true);
      } else {
	Translation_mutation_effect(&codon_index,&codon_orig,&codon_mut,
				    pairs,npairs,translation_frame,translation_start,
				    translation_end,mutindex,/*newnt*/'X',
				    /*backwardsp*/false,/*revcompp*/false);
      }
      if (change[0] == 'i' || change[0] == '^') {
	printf("Frameshift insertion at %c%d",codon_orig,codon_index);
      } else {
	printf("Frameshift deletion at %c%d",codon_orig,codon_index);
      }
    } else if (change[0] == 'n') {
      printf("Coding region, ambiguous");
    } else if (change[0] == 'x') {
      newnt = change[1];
      if (newnt != 'A' && newnt != 'C' && newnt != 'G' && newnt != 'T') {
	printf("Coding region, ambiguous");
      } else {
	if (cdna_direction < 0) {
	  Translation_mutation_effect(&codon_index,&codon_orig,&codon_mut,
				      pairs,npairs,translation_frame,translation_start,
				      translation_end,mutindex,newnt,
				      /*backwardsp*/true,/*revcompp*/true);
	} else {
	  Translation_mutation_effect(&codon_index,&codon_orig,&codon_mut,
				      pairs,npairs,translation_frame,translation_start,
				      translation_end,mutindex,newnt,
				      /*backwardsp*/false,/*revcompp*/false);
	}
	if (codon_orig == codon_mut) {
	  printf("Same: %c%d",codon_orig,codon_index);
	} else {
	  printf("Change: %c%d%c",codon_orig,codon_index,codon_mut);
	}
      }
    } else {
      printf("Coding region, but can't parse %s",change);
    }
  }
  printf("\n");
  return;
}
#endif


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
    if (this->aamarker == true) {
      printf(" => %c %c",this->aa_g,this->aa_e);
    }
    printf("\n");
  }
  return;
}  


void
Pair_dump_aapos (struct T *pairs, int npairs, int aapos, int cdna_direction) {
  struct T *this, *start, *end;
  int i, ngenome = 0, ncdna = 0;
  char complCode[128] = COMPLEMENT;

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
	printf("%c",this->comp);
      } else {
	printf("%c",complCode[this->genome]);
      }
    }
    printf(">");
    for (this = end; this >= start; --this) {
      if (this->gapp == true) {
      } else if (this->cdna == ' ') {
	printf("%c",this->comp);
      } else {
	printf("%c",complCode[this->cdna]);
      }
    }

  } else {

    for (this = start; this <= end; this++) {
      if (this->gapp == true) {
      } else if (this->genome == ' ') {
	printf("%c",this->comp);
      } else {
	printf("%c",this->genome);
      }
    }
    printf(">");
    for (this = start; this <= end; this++) {
      if (this->gapp == true) {
      } else if (this->cdna == ' ') {
	printf("%c",this->comp);
      } else {
	printf("%c",this->cdna);
      }
    }
  
  }
  return;
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
  int i;


  if (pairs == NULL) {
    return;
  } else {
    this = List_head(pairs);
    prev_querypos = this->querypos;
    prev_genomepos = this->genomepos;

    for (p = List_next(pairs); p != NULL; p = List_next(p)) {
      this = List_head(p);
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
  switch (toupper(c)) {
  case 'A': case 'C': case 'G': case 'T': return false;
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
  int exon_querystart = -1, exon_genomestart = -1, exon_queryend, exon_genomeend;
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
	if (this->comp == '>') {
	  printf(" ->\n");
	} else if (this->comp == '<') {
	  printf(" <-\n");
	} else if (this->comp == ')') {
	  printf(" -)\n");
	} else if (this->comp == '(') {
	  printf(" (-\n");
	} else if (this->comp == ']') {
	  printf(" -]\n");
	} else if (this->comp == '[') {
	  printf(" [-\n");
	} else if (this->comp == '=') {
	  printf(" ==\n");
	} else {
	  printf(" ##\n");
	}
	in_exon = false;
      }
    } else if (this->comp == '.') {
      /* Do nothing */
    } else {
      /* Remaining possibilities are '|', '*', '-', '~', or ' ' */
      if (in_exon == false) {
	exon_querystart = this->querypos + !zerobasedp;
	if (watsonp) {
	  exon_genomestart = chrpos + this->genomepos + !zerobasedp;
	} else {
	  exon_genomestart = chrpos + (genomiclength - 1) - this->genomepos + !zerobasedp;
	}
	num = den = 0;
	in_exon = true;
      }
      if (this->comp == '-' || this->comp == '~') {
	/* Don't count in numerator or denominator */
      } else if (unknown_base(this->cdna) || unknown_base(this->genome)) {
	/* Comp must be a space */
	/* Don't count in numerator or denominator */
      } else {
	den++;
	if (this->comp == '|' || this->comp == '*') {
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

void
Pair_print_cdna_exons (struct T *pairs, int npairs, int wraplength) {
  bool in_exon = false;
  struct T *ptr, *prev, *this = NULL;
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
      }
    } else if (this->comp == '.') {
      /* Do nothing */
    } else {
      /* Remaining possibilities are '|', '*', '-', '~', or ' ' */
      if (in_exon == false) {
	printf("<exon %d>\n",++exonno);
	in_exon = true;
      }
      printf("%c",this->cdna);
      if (++column % wraplength == 0) {
	printf("\n");
	column = 0;
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
	  if (this->comp == '>') {
	    (*ncanonical)++;
	  } else if (this->comp == ')' || this->comp == ']') {
	    (*nsemicanonical)++;
	  } else if (this->comp == '=') {
	    (*nnoncanonical)++;
	  }

	} else if (cdna_direction < 0) {
	  if (this->comp == '<') {
	    (*ncanonical)++;
	  } else if (this->comp == '(' || this->comp == '[') {
	    (*nsemicanonical)++;
	  } else if (this->comp == '=') {
	    (*nnoncanonical)++;
	  }
	}
	in_intron = true;
      }
    } else {
      if (in_intron) {
	in_intron = false;
      }
      if (this->comp == '-') {
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
      } else if (this->comp == '~') {
	/* Do nothing */
      } else if (unknown_base(this->cdna) || unknown_base(this->genome)) {
	(*unknowns)++;
      } else if (this->comp == '|' || this->comp == '*') {
	(*matches)++;
      } else if (this->comp == ' ') {
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
  List_T p;
  T this, prev = NULL;
  int i;

  *matches = *unknowns = *mismatches = *qopens = *qindels = *topens = *tindels = 
    *ncanonical = *nsemicanonical = *nnoncanonical = 0;
  
  for (i = 0, this = ptr; i < npairs; i++, this = ptr++) {
    if (this->gapp) {
      if (!in_intron) {
	if (this->querypos >= minpos && this->querypos <= maxpos) {
	  if (this->comp == '>') {
	    (*ncanonical)++;
	  } else if (this->comp == ')' || this->comp == ']') {
	    (*nsemicanonical)++;
	  } else if (this->comp == '=') {
	    (*nnoncanonical)++;
	  }
	} else if (cdna_direction < 0) {
	  if (this->comp == '<') {
	    (*ncanonical)++;
	  } else if (this->comp == '(' || this->comp == '[') {
	    (*nsemicanonical)++;
	  } else if (this->comp == '=') {
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
	if (this->comp == '-') {
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
	} else if (this->comp == '~') {
	  /* Do nothing */
	} else if (unknown_base(this->cdna) || unknown_base(this->genome)) {
	  (*unknowns)++;
	} else if (this->comp == '|' || this->comp == '*') {
	  (*matches)++;
	} else if (this->comp == ' ') {
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


void
Pair_pathscores (int *pathscores, struct T *ptr, int npairs, 
		 int cdna_direction, int querylength) {
  int querypos;
  bool in_intron = false;
  List_T p;
  T this, prev = NULL;
  int i;

  for (i = 0, this = ptr; i < npairs; i++, this = ptr++) {
    querypos = this->querypos;
    if (this->gapp) {
      if (!in_intron) {
	if (cdna_direction > 0) {
	  if (this->comp == '>') {
	    pathscores[querypos] = CANONICAL_POINTS;
	  } else if (this->comp == ')' || this->comp == ']') {
	    pathscores[querypos] = SEMICANONICAL_POINTS;
	  } else if (this->comp == '=') {
	    pathscores[querypos] = 0; /* noncanonical */
	  }
	} else if (cdna_direction < 0) {
	  if (this->comp == '<') {
	    pathscores[querypos] = CANONICAL_POINTS;
	  } else if (this->comp == '(' || this->comp == '[') {
	    pathscores[querypos] = SEMICANONICAL_POINTS;
	  } else if (this->comp == '=') {
	    pathscores[querypos] = 0; /* noncanonical */
	  }
	}
	in_intron = true;
      }
    } else {
      if (in_intron) {
	in_intron = false;
      }
      if (this->comp == '-') {
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
      } else if (this->comp == '~') {
	/* Do nothing */
      } else if (unknown_base(this->cdna) || unknown_base(this->genome)) {
	/* (*unknowns)++; */
      } else if (this->comp == '|' || this->comp == '*') {
	pathscores[querypos] = +1; /* For match */
      } else if (this->comp == ' ') {
	pathscores[querypos] = MISMATCH;
      } else {
	fprintf(stderr,"Can't parse comp %c, gapp %d\n",this->comp,this->gapp);
	abort();
      }
    }
    prev = this;
  }

  for (i = 1; i < querylength; i++) {
    pathscores[i] += pathscores[i-1];
  }

  return;
}


int
Pair_nexons (struct T *pairs, int npairs) {
  int nexons = 0;
  bool in_exon = false;
  struct T *this;
  int num = 0, den = 0, i;

  for (i = 0; i < npairs; i++) {
    this = pairs++;
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
	if (this->comp == '>') {
	  cdna_direction += 1;
	} else if (this->comp == '<') {
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
	  if (this->comp == '<' || this->comp == '(' || this->comp == '[') {
	    return false;
	  } else if (this->comp == '>') {
	    (*ncanonical)++;
	  }
	} else if (cdna_direction < 0) {
	  if (this->comp == '>' || this->comp == ')' || this->comp == ']') {
	    return false;
	  } else if (this->comp == '<') {
	    (*ncanonical)++;
	  }
	} else if (cdna_direction == 0) {
	  /* Set cdna_direction for next time */
	  if (this->comp == '>' || this->comp == ')' || this->comp == ']') {
	    cdna_direction = +1;
	  } else if (this->comp == '<' || this->comp == '(' || this->comp == '[') {
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
  struct T *ptr;
  int j, k;
  char newmatch;

  j = -gapcount;
  if (Intron_canonical_fwd_p(pairs[j].genome,pairs[j+1].genome,pairs[-2].genome,pairs[-1].genome)) {
    newmatch = '>';
  } else if (Intron_canonical_rev_p(pairs[j].genome,pairs[j+1].genome,pairs[-2].genome,pairs[-1].genome)) {
    newmatch = '<';
  } else if (Intron_gcag_fwd_p(pairs[j].genome,pairs[j+1].genome,pairs[-2].genome,pairs[-1].genome)) {
    newmatch = ')';
  } else if (Intron_gcag_rev_p(pairs[j].genome,pairs[j+1].genome,pairs[-2].genome,pairs[-1].genome)) {
    newmatch = '(';
  } else if (Intron_atac_fwd_p(pairs[j].genome,pairs[j+1].genome,pairs[-2].genome,pairs[-1].genome)) {
    newmatch = ']';
  } else if (Intron_atac_rev_p(pairs[j].genome,pairs[j+1].genome,pairs[-2].genome,pairs[-1].genome)) {
    newmatch = '[';
  } else {
    newmatch = '=';
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


int
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
  struct T *chunk, *chunk2, *ptr;
  T newpair, oldpair;
  List_T p;
  int ncopy, nmerged, tailgap, i, j;

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

  if ((nmerged = merge_gaps(chunk,*npairs,ngap)) == 0) {
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
      if (this->aa_g != '*') {
	printf("%c",this->aa_g);
	xpos++;
      }
    }
  }
  printf("\n");
  return;
}

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
	printf("%c",this->aa_e);
	xpos++;
      }
    }
  }
  printf("\n");
  return;
}


void
Pair_print_compressed (Sequence_T queryseq, char *version, int pathnum, int npaths,
		       int nexons, double coverage, double fracidentity,
		       struct T *pairs, int npairs, Chrnum_T chrnum, Genomicpos_T chrpos,
		       Genomicpos_T chroffset, IIT_T chromosome_iit, 
		       Genomicpos_T genomiclength, bool checksump,
		       bool chimerap, char *strain, bool watsonp, bool zerobasedp) {
  Genomicpos_T chrpos1, chrpos2, position1, position2;
  char *chrstring = NULL;

  bool in_exon = false;
  List_T tokens = NULL;
  struct T *ptr = pairs, *prev, *this = NULL;
  T start, end;
  int exon_querystart = -1, exon_genomestart = -1, exon_queryend, exon_genomeend,
    intron_start, intron_end;
  int num = 0, den = 0, runlength = 0, tokencount = 1, i;
  int print_dinucleotide_p;
  char token[10], donor[3], acceptor[3];
  unsigned char *digest;

  donor[0] = donor[1] = donor[2] = '\0';
  acceptor[0] = acceptor[1] = acceptor[2] = '\0';

  printf(">%s %s %d/%d %d %d",
	 Sequence_accession(queryseq),version,pathnum,npaths,
	 Sequence_length_full(queryseq),nexons);
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

  if (checksump == true) {
    printf(" md5:");
    Sequence_print_digest(queryseq);
  }

  if (chimerap == true) {
    printf(" chimera:T");
  }

  if (strain != NULL) {
    printf(" strain:%s",strain);
  }

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
	if (this->comp == '>') {
	  sprintf(token,"%d>",runlength);
	} else if (this->comp == '<') {
	  sprintf(token,"%d<",runlength);
	  print_dinucleotide_p = -1;
	} else if (this->comp == '=') {
	  sprintf(token,"%d=",runlength);
	} else if (this->comp == ')') {
	  sprintf(token,"%d)",runlength);
	} else if (this->comp == '(') {
	  sprintf(token,"%d(",runlength);
	  print_dinucleotide_p = -1;
	} else if (this->comp == ']') {
	  sprintf(token,"%d]",runlength);
	} else if (this->comp == '[') {
	  sprintf(token,"%d[",runlength);
	  print_dinucleotide_p = -1;
	} else if (this->comp == '#') {
	  sprintf(token,"%d#",runlength);
	  print_dinucleotide_p = 0;
	} else {
	  fprintf(stderr,"Can't parse %c in compression\n",this->comp);
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
    } else if (this->comp == '.') {
      /* Do nothing */
    } else {
      /* Remaining possibilities are '|', '*', '-', '~', or ' ' */
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
      if (this->comp == '-' || this->comp == '~') {
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
      } else if (this->comp == ' ') {
	/* Mismatch, possibly due to unknown in either sequence */
	if (this->genome == 'N' && this->cdna == 'N') {
	  sprintf(token,"%d?",runlength);
	} else if (this->genome == 'N') {
	  sprintf(token,"%dN%c",runlength,this->cdna);
	} else if (this->cdna == 'N') {
	  sprintf(token,"%dn",runlength);
	} else {
	  sprintf(token,"%dx%c",runlength,this->cdna);
	}
	tokens = push_token(tokens,token);
	runlength = 0;
	den++;

      } else {
	runlength++;
	den++;
	if (this->comp == '|' || this->comp == '*') {
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
  
  printf("\t%u %u",exon_genomestart,exon_genomeend);
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



