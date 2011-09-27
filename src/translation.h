/* $Id: translation.h,v 1.23 2007-04-23 16:10:23 twu Exp $ */
#ifndef TRANSLATION_INCLUDED
#define TRANSLATION_INCLUDED
#include "bool.h"
#include "pair.h"

#define T Translation_T
typedef struct T *T;

extern char
Translation_get_codon (char a, char b, char c);

#ifdef PMAP
extern void
Translation_via_cdna (int *translation_leftpos, int *translation_rightpos, int *translation_length,
		      int *relaastart, int *relaaend,
		      struct Pair_T *pairs, int npairs, char *queryaaseq_ptr, bool strictp);
#else
extern void
Translation_via_genomic (int *translation_leftpos, int *translation_rightpos, int *translation_length,
			 int *relaastart, int *relaaend,
			 struct Pair_T *pairs, int npairs, bool backwardp, bool revcompp, bool fulllengthp,
			 bool strictp);
#endif

extern void
Translation_via_reference (int *relaastart, int *relaaend,
			   struct Pair_T *pairs, int npairs, bool watsonp, bool backwardp, bool revcompp,
			   struct Pair_T *refpairs, int nrefpairs, bool refwatsonp, int genomiclength,
			   bool fixshiftp);

extern void
Translation_compare (struct Pair_T *pairs, int npairs, struct Pair_T *refpairs, int nrefpairs,
		     int cdna_direction, int relaastart, int relaaend, int maxmutations);

#undef T
#endif




 
