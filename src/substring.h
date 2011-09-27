/* $Id: substring.h,v 1.4 2010/03/05 05:51:57 twu Exp $ */
#ifndef SUBSTRING_INCLUDED
#define SUBSTRING_INCLUDED

#include "genomicpos.h"
#include "chrnum.h"
#include "sequence.h"
#include "genome.h"
#include "iit-read.h"
#include "bool.h"

/* Should arrange in order of goodness, best to worst */
typedef enum {EXACT, SUB, INS, DEL, SPLICE} Hittype_T;

typedef struct Substring_T *Substring_T;

extern void
Substring_print_nsnpdiffs (bool labelsp);
extern void
Substring_print_ncolordiffs ();


extern Substring_T
Substring_new (int nmismatches, int ncolordiffs, Chrnum_T chrnum, Genomicpos_T chroffset,
	       Genomicpos_T genomicstart, Genomicpos_T genomicend,
	       int querystart, int queryend, int querylength,
	       Genomicpos_T alignstart, Genomicpos_T alignend, int genomiclength,
	       int extraleft, int extraright, char *genomicseg, char *query,
	       bool plusp, bool trim_ends_p, bool dibasep, bool cmetp);

extern void
Substring_free (Substring_T *old);

extern bool
Substring_plusp (Substring_T this);
extern int
Substring_nmismatches (Substring_T this);
extern int
Substring_ncolordiffs (Substring_T this);
extern int
Substring_trim_left (Substring_T this);
extern int
Substring_trim_right (Substring_T this);
extern int
Substring_match_length (Substring_T this);

extern Chrnum_T
Substring_chrnum (Substring_T this);
extern Genomicpos_T
Substring_chroffset (Substring_T this);
extern Genomicpos_T
Substring_genomicstart (Substring_T this);
extern Genomicpos_T
Substring_genomicend (Substring_T this);
extern Genomicpos_T
Substring_genomiclength (Substring_T this);

extern double
Substring_chimera_prob (Substring_T this);
extern int
Substring_chimera_pos (Substring_T this);
extern bool
Substring_chimera_knownp (Substring_T this);
extern bool
Substring_chimera_sensep (Substring_T this);


extern Substring_T
Substring_copy (Substring_T old);

extern Substring_T
Substring_new_donor (int donor_pos, int donor_nmismatches, int donor_ncolordiffs,
		     double donor_prob, Genomicpos_T left, int querylength,
		     bool plusp, bool sensep, char *genomicseg, char *query, Chrnum_T chrnum,
		     Genomicpos_T chroffset, bool knownp, bool trim_ends_p, bool dibasep, bool cmetp);
extern Substring_T
Substring_new_acceptor (int acceptor_pos, int acceptor_nmismatches, int acceptor_ncolordiffs,
			double acceptor_prob, Genomicpos_T left, int querylength,
			bool plusp, bool sensep, char *genomicseg, char *query, Chrnum_T chrnum,
			Genomicpos_T chroffset, bool knownp, bool trim_ends_p, bool dibasep, bool cmetp);
extern List_T
Substring_sort_chimera_halves (List_T hitlist, bool ascendingp);


extern void
Substring_print_single (Substring_T substring, Hittype_T hittype, Sequence_T queryseq,
			char *chr, int querylength,
			bool invertp, IIT_T snps_iit, int *snps_divint_crosstable);
extern void
Substring_print_insertion_1 (Substring_T substring1, Substring_T substring2, int nindels, char *chr,
			     bool invertp, IIT_T snps_iit, int *snps_divint_crosstable);
extern void
Substring_print_insertion_2 (Substring_T substring1, Substring_T substring2, int nindels, char *chr,
			     bool invertp, IIT_T snps_iit, int *snps_divint_crosstable);
extern void
Substring_print_deletion_1 (Substring_T substring1, Substring_T substring2, int nindels, 
			    char *deletion, char *chr,
			    bool invertp, IIT_T snps_iit, int *snps_divint_crosstable);
extern void
Substring_print_deletion_2 (Substring_T substring1, Substring_T substring2, int nindels, 
			    char *deletion, char *chr,
			    bool invertp, IIT_T snps_iit, int *snps_divint_crosstable);
extern void
Substring_print_donor (Substring_T donor, bool sensep, bool invertp,
		       IIT_T chromosome_iit, IIT_T snps_iit, int *snps_divint_crosstable, IIT_T splicesites_iit,
		       int donor_typeint, int *splicesites_divint_crosstable, bool endp,
		       Substring_T acceptor, Genomicpos_T chimera_distance);
extern void 
Substring_print_acceptor (Substring_T acceptor, bool sensep, bool invertp,
			  IIT_T chromosome_iit, IIT_T snps_iit, int *snps_divint_crosstable, IIT_T splicesites_iit,
			  int acceptor_typeint, int *splicesites_divint_crosstable, bool endp,
			  Substring_T donor, Genomicpos_T chimera_distance);

extern void
Substring_assign_donor_prob (Substring_T donor, Genome_T genome, IIT_T chromosome_iit);
extern void
Substring_assign_acceptor_prob (Substring_T acceptor, Genome_T genome, IIT_T chromosome_iit);


extern int
Substring_geneprob_eval_single (Substring_T substring, IIT_T geneprob_iit, IIT_T chromosome_iit);
extern int
Substring_geneprob_eval_insertion (Substring_T substring1, Substring_T substring2, int indel_pos, int nindels,
				   IIT_T geneprob_iit, IIT_T chromosome_iit, Sequence_T queryseq);
extern int
Substring_geneprob_eval_deletion (Substring_T substring1, Substring_T substring2, int indel_pos, int nindels,
				  IIT_T geneprob_iit, IIT_T chromosome_iit, Sequence_T queryseq);
extern int
Substring_geneprob_eval_splice (Substring_T donor, Substring_T acceptor, IIT_T geneprob_iit, IIT_T chromosome_iit, int querylength,
				bool sensep);


#endif


