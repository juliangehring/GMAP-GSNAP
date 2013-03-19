/* $Id: iit-read.h 80072 2012-11-28 07:44:07Z twu $ */
#ifndef IIT_READ_INCLUDED
#define IIT_READ_INCLUDED
#include <stdio.h>
#include "bool.h"
#include "uintlist.h"
#include "list.h"
#include "interval.h"

typedef enum {READ_ALL, READ_ONE, READ_NONE} Divread_T;
/* READ_NONE is useful if we want to obtain an interval by name,
   rather than by coordinate */

typedef enum {NO_KNOWN_GENE, KNOWN_GENE, KNOWN_GENE_MULTIEXON} Overlap_T;


#define T IIT_T
#ifndef IIT_TYPEDEF
#define IIT_TYPEDEF
typedef struct T *T;
#endif

extern char *
IIT_name (T this);
extern int
IIT_version (T this);
extern int
IIT_total_nintervals (T this);
extern int
IIT_nintervals (T this, int divno);
extern int
IIT_ntypes (T this);
extern int
IIT_nfields (T this);

extern unsigned int
IIT_length (T this, int index);
extern unsigned int
IIT_divlength (T this, char *divstring);
extern unsigned int
IIT_totallength (T this);
extern unsigned int
IIT_genomelength (T chromosome_iit, bool with_circular_alias_p);
extern bool *
IIT_circularp (T chromosome_iit);
extern Interval_T
IIT_interval (T this, int index);
extern unsigned int
IIT_interval_low (T this, int index);
extern unsigned int
IIT_interval_high (T this, int index);
extern unsigned int
IIT_interval_length (T this, int index);
extern int
IIT_interval_type (T this, int index);
extern unsigned int
IIT_next_chrbound (T this, int index, int circular_typeint);
extern int
IIT_interval_sign (T this, int index);
extern void
IIT_interval_bounds (unsigned int *low, unsigned int *high, unsigned int *length, T this,
		     int index, int circular_typeint);
extern int
IIT_index (T this, int divno, int i);

extern int
IIT_ndivs (T this);
extern char *
IIT_divstring (T this, int divno);
extern int
IIT_divint (T this, char *divstring);
extern int *
IIT_divint_crosstable (T chromosome_iit, T iit);
extern char *
IIT_divstring_from_index (T this, int index);
extern char *
IIT_typestring (T this, int type);
extern int
IIT_typeint (T this, char *typestring);
extern char *
IIT_fieldstring (T this, int fieldint);
extern char *
IIT_label (T this, int index, bool *allocp);
extern char *
IIT_annotation (char **restofheader, T this, int index, bool *alloc_header_p);
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
IIT_dump_divstrings (FILE *fp, T this);
extern void
IIT_dump_typestrings (FILE *fp, T this);
extern void
IIT_dump_fieldstrings (FILE *fp, T this);
extern void
IIT_dump_labels (FILE *fp, T this);
extern void
IIT_dump (T this, bool annotationonlyp, bool sortp);
extern void
IIT_dump_simple (T this);
extern void
IIT_dump_sam (FILE *fp, T this, char *sam_read_group_id, char *sam_read_group_name,
	      char *sam_read_group_library, char *sam_read_group_platform);
extern void
IIT_dump_version1 (T this, T chromosome_iit, bool directionalp);
extern void
IIT_dump_formatted (T this, bool directionalp);
extern unsigned int *
IIT_transitions (int **signs, int *nedges, T this);
extern unsigned int *
IIT_transitions_subset (int **signs, int *nedges, T this, int *indices, int nindices);
extern void
IIT_dump_counts (T this, bool alphabetizep);

extern void
IIT_free (T *old);
extern int
IIT_read_divint (char *filename, char *divstring, bool add_iit_p);
extern T
/* add_iit_p means to add the ".iit" suffix to the filename */
IIT_read (char *filename, char *name, bool readonlyp, Divread_T divread, char *divstring,
	  bool add_iit_p, bool labels_read_p);

extern int *
IIT_find (int *nmatches, T this, char *label);
extern int
IIT_find_linear (T this, char *label);
extern int
IIT_find_one (T this, char *label);

extern unsigned int *
IIT_get_highs_for_low (int *nuniq, T this, int divno, unsigned int x);
extern unsigned int *
IIT_get_lows_for_high (int *nuniq, T this, int divno, unsigned int x);
extern bool
IIT_low_exists_signed_p (T this, int divno, unsigned int x, int sign);
extern bool
IIT_high_exists_signed_p (T this, int divno, unsigned int x, int sign);


extern int *
IIT_get (int *nmatches, T this, char *divstring, unsigned int x, unsigned int y, bool sortp);
extern int *
IIT_get_signed (int *nmatches, T this, char *divstring, unsigned int x, unsigned int y, int sign, bool sortp);
extern bool
IIT_exists_with_divno (T this, int divno, unsigned int x, unsigned int y);
extern bool
IIT_exists_with_divno_signed (T this, int divno, unsigned int x, unsigned int y, int sign);
extern bool
IIT_exists_with_divno_typed_signed (T this, int divno, unsigned int x, unsigned int y, int type, int sign);
extern int *
IIT_get_with_divno (int *nmatches, T this, int divno, unsigned int x, unsigned int y, bool sortp);
extern int *
IIT_get_signed_with_divno (int *nmatches, T this, int divno, unsigned int x, unsigned int y, bool sortp,
			   int sign);
extern void
IIT_get_flanking (int **leftflanks, int *nleftflanks, int **rightflanks, int *nrightflanks,
		  T this, char *divstring, unsigned int x, unsigned int y, int nflanking, int sign);
extern void
IIT_get_flanking_with_divno (int **leftflanks, int *nleftflanks, int **rightflanks, int *nrightflanks,
			     T this, int divno, unsigned int x, unsigned int y, int nflanking, int sign);
extern void
IIT_get_flanking_typed (int **leftflanks, int *nleftflanks, int **rightflanks, int *nrightflanks,
			T this, char *divstring, unsigned int x, unsigned int y, int nflanking, int type,
			int sign);
extern void
IIT_get_flanking_multiple_typed (int **leftflanks, int *nleftflanks, int **rightflanks, int *nrightflanks,
				 T this, char *divstring, unsigned int x, unsigned int y, int nflanking, int *types, int ntypes);
extern int
IIT_get_one (T this, char *divstring, unsigned int x, unsigned int y);
extern int *
IIT_get_typed (int *ntypematches, T this, char *divstring, unsigned int x, unsigned int y, int type, bool sortp);
extern int *
IIT_get_typed_with_divno (int *ntypematches, T this, int divno, unsigned int x, unsigned int y, int type, bool sortp);
extern int *
IIT_get_typed_signed (int *ntypematches, T this, char *divstring, unsigned int x, unsigned int y,
		      int type, int sign, bool sortp);
extern int *
IIT_get_typed_signed_with_divno (int *ntypematches, T this, int divno, unsigned int x, unsigned int y, 
				 int type, int sign, bool sortp);
extern int *
IIT_get_multiple_typed (int *ntypematches, T this, char *divstring, unsigned int x, unsigned int y, 
			int *types, int ntypes, bool sortp);
extern int
IIT_get_exact (T this, char *divstring, unsigned int x, unsigned int y, int type);
extern bool
IIT_exact_p (T this, char *divstring, unsigned int x, unsigned int y, int type);
extern int *
IIT_get_exact_multiple (int *nmatches, T this, char *divstring, unsigned int x, unsigned int y, int type);
extern int *
IIT_get_exact_multiple_with_divno (int *nmatches, T this, int divno, unsigned int x, unsigned int y, int type);

extern List_T
IIT_intervallist_typed (List_T *labellist, Uintlist_T *seglength_list, T this);
extern List_T
IIT_typelist (T this);

extern char *
IIT_string_from_position (unsigned int *chrpos, unsigned int position, 
			  T chromosome_iit);
extern void
IIT_print_header (FILE *fp, T this, int *matches, int nmatches, bool map_bothstrands_p,
		  char *chr, bool reversep, bool relativep, unsigned int left, bool print_comment_p);

extern Overlap_T
IIT_gene_overlap (T map_iit, int divno, unsigned int x, unsigned int y, bool favor_multiexon_p);

#undef T
#endif
