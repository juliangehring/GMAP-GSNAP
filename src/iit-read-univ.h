/* $Id: iit-read-univ.h 99737 2013-06-27 19:33:03Z twu $ */
#ifndef IIT_READ_UNIV_INCLUDED
#define IIT_READ_UNIV_INCLUDED
#include <stdio.h>
#include "bool.h"
#include "uintlist.h"
#include "list.h"
#include "univinterval.h"
#include "types.h"
#include "iit-read.h"		/* For IIT_divint */

#define T Univ_IIT_T
typedef struct T *T;

extern bool
Univ_IIT_coord_values_8p (T this);
extern int
Univ_IIT_total_nintervals (T this);
extern int
Univ_IIT_ntypes (T this);
extern Univcoord_T
Univ_IIT_length (T this, int index);
extern Univcoord_T
Univ_IIT_genomelength (T chromosome_iit, bool with_circular_alias_p);
extern bool *
Univ_IIT_circularp (T chromosome_iit);
extern Univinterval_T
Univ_IIT_interval (T this, int index);
extern Univcoord_T
Univ_IIT_interval_low (T this, int index);
extern Univcoord_T
Univ_IIT_interval_high (T this, int index);
extern Univcoord_T
Univ_IIT_interval_length (T this, int index);
extern int
Univ_IIT_interval_type (T this, int index);
extern Univcoord_T
Univ_IIT_next_chrbound (T this, int index, int circular_typeint);
extern void
Univ_IIT_interval_bounds (Univcoord_T *low, Univcoord_T *high, Chrpos_T *length, T this,
			  int index, int circular_typeint);
extern void
Univ_IIT_intervals_setup (Univcoord_T **chroffsets, Univcoord_T **chrhighs, Chrpos_T **chrlengths,
			  T this, int nchromosomes, int circular_typeint);
extern int *
Univ_IIT_divint_crosstable (T chromosome_iit, IIT_T iit);
extern char *
Univ_IIT_typestring (T this, int type);
extern int
Univ_IIT_typeint (T this, char *typestring);
extern char *
Univ_IIT_label (T this, int index, bool *allocp);
extern char *
Univ_IIT_annotation (char **restofheader, T this, int index, bool *alloc_header_p);
extern void
Univ_IIT_dump_typestrings (FILE *fp, T this);
extern void
Univ_IIT_dump_labels (FILE *fp, T this);
extern void
Univ_IIT_dump (T this);
extern void
Univ_IIT_dump_table (T this, bool zerobasedp);
extern void
Univ_IIT_dump_fai (T this);
extern void
Univ_IIT_dump_sam (FILE *fp, T this, char *sam_read_group_id, char *sam_read_group_name,
		   char *sam_read_group_library, char *sam_read_group_platform);
extern void
Univ_IIT_dump_labels (FILE *fp, T this);
extern char
Univ_IIT_annotation_firstchar (T this, int index);
extern void
Univ_IIT_dump_contigs (T this, T chromosome_iit, bool directionalp);
extern void
Univ_IIT_free (T *old);
extern T
Univ_IIT_read (char *filename, bool readonlyp, bool add_iit_p);
extern void
Univ_IIT_debug (char *filename);
extern int *
Univ_IIT_find (int *nmatches, T this, char *label);
extern int
Univ_IIT_find_linear (T this, char *label);
extern int
Univ_IIT_find_one (T this, char *label);
extern int *
Univ_IIT_get (int *nmatches, T this, Univcoord_T x, Univcoord_T y);
extern int
Univ_IIT_get_one (T this, Univcoord_T x, Univcoord_T y);
extern char *
Univ_IIT_string_from_position (Chrpos_T *chrpos, Univcoord_T position, T chromosome_iit);

#undef T
#endif

