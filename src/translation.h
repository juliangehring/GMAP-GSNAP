/* $Id: translation.h,v 1.14 2005/02/15 01:58:51 twu Exp $ */
#ifndef TRANSLATION_INCLUDED
#define TRANSLATION_INCLUDED
#include "bool.h"
#include "pair.h"

#define T Translation_T
typedef struct T *T;

extern void
Translation_via_genomic (int *translation_leftpos, int *translation_rightpos, int *translation_length,
			 int *relaastart, int *relaaend,
			 struct Pair_T *pairs, int npairs, bool backwardp, bool revcompp, bool fulllengthp);

extern void
Translation_via_reference (int *relaastart, int *relaaend,
			   struct Pair_T *pairs, int npairs, bool watsonp, bool backwardp, bool revcompp,
			   struct Pair_T *refpairs, int nrefpairs, bool refwatsonp, int genomiclength,
			   bool fixshiftp);

extern void
Translation_compare (struct Pair_T *pairs, int npairs, struct Pair_T *refpairs, int nrefpairs,
		     int cdna_direction, int relaastart, int relaaend);

#undef T
#endif




 
