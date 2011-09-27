/* $Id: iit-read.h,v 1.47 2007/09/19 17:48:31 twu Exp $ */
#ifndef IIT_READ_INCLUDED
#define IIT_READ_INCLUDED
#include <stdio.h>
#include "bool.h"
#include "uintlist.h"
#include "list.h"
#include "interval.h"

#define T IIT_T
#ifndef IIT_TYPEDEF
#define IIT_TYPEDEF
typedef struct T *T;
#endif

extern char *
IIT_name (T this);
extern int
IIT_nintervals (T this);
extern int
IIT_ntypes (T this);
extern int
IIT_nfields (T this);
extern char **
IIT_types (int *ntypes, T this, bool alphabetizep);
extern char **
IIT_fields (int *nfields, T this);

extern unsigned int
IIT_length (T this, int index);
extern unsigned int
IIT_totallength (T this);
extern Interval_T
IIT_interval (T this, int index);
extern char *
IIT_typestring (T this, int type);
extern int
IIT_typeint (T this, char *typestring);
extern char *
IIT_fieldstring (T this, int fieldint);
extern char *
IIT_label (T this, int index);
extern char *
IIT_annotation (T this, int index, bool *allocp);
extern char
IIT_annotation_firstchar (T this, int index);
extern unsigned int
IIT_annotation_strlen (T this, int index);
extern char *
IIT_fieldvalue (T this, int index, int fieldint);
extern int
IIT_fieldint (T this, char *fieldstring);

extern void
IIT_debug (char *filename);
extern void
IIT_dump_typestrings (FILE *fp, T this);
extern void
IIT_dump_fieldstrings (FILE *fp, T this);
extern void
IIT_dump_labels (FILE *fp, T this);
extern void
IIT_dump (T this, bool annotationonlyp);
extern void
IIT_dump_formatted (T this, bool directionalp);
extern unsigned int *
IIT_transitions (int **signs, int *nedges, T this);
extern void
IIT_dump_counts (T this, bool alphabetizep);

extern void
IIT_free (T *old);
extern T
IIT_read (char *filename, char *name, bool readonlyp);

extern int *
IIT_find (int *nmatches, T this, char *label);
extern int
IIT_find_linear (T this, char *label);
extern int *
IIT_find_multiple (int *nmatches, T this, char **labels, int nlabels);
extern int
IIT_find_one (T this, char *label);

extern int *
IIT_get_all (int *nmatches, T this, bool sortp);
extern int *
IIT_get (int *nmatches, T this, unsigned int x, unsigned int y, bool sortp);
extern int *
IIT_get_all_typed (int *nmatches, T this, int *types, int ntypes, bool sortp);
extern void
IIT_get_flanking (int **leftflanks, int *nleftflanks, int **rightflanks, int *nrightflanks,
		  T this, unsigned int x, unsigned int y, int nflanking, int sign);
extern void
IIT_get_flanking_typed (int **leftflanks, int *nleftflanks, int **rightflanks, int *nrightflanks,
			T this, unsigned int x, unsigned int y, int nflanking, int type);
extern void
IIT_get_flanking_multiple_typed (int **leftflanks, int *nleftflanks, int **rightflanks, int *nrightflanks,
				 T this, unsigned int x, unsigned int y, int nflanking, int *types, int ntypes);
extern int
IIT_get_one (T this, unsigned int x, unsigned int y);
extern int *
IIT_get_typed (int *ntypematches, T this, unsigned int x, unsigned int y, int type, bool sortp);
extern int *
IIT_get_multiple_typed (int *ntypematches, T this, unsigned int x, unsigned int y, 
			int *types, int ntypes, bool sortp);
extern int
IIT_get_exact (T this, unsigned int x, unsigned int y, int type);
extern int *
IIT_get_exact_multiple (int *nmatches, T this, unsigned int x, unsigned int y, int type);
extern void
IIT_print (T this, int *matches, int nmatches, bool map_bothstrands_p,
	   T chromosome_iit, int *levels, bool reversep);

extern List_T
IIT_intervallist_typed (List_T *labellist, Uintlist_T *seglength_list, T this);
extern List_T
IIT_typelist (T this);


#undef T
#endif
