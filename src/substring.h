/* $Id: substring.h,v 1.11 2010-07-27 00:02:38 twu Exp $ */
#ifndef SUBSTRING_INCLUDED
#define SUBSTRING_INCLUDED

#include "genomicpos.h"
#include "chrnum.h"
#include "sequence.h"
#include "genome.h"
#include "iit-read.h"
#include "bool.h"

/* Should arrange in order of goodness, best to worst */
typedef enum {EXACT, SUB, INS, DEL, SPLICE, TERMINAL} Hittype_T;

#define T Substring_T
typedef struct T *T;

extern void
Substring_print_nsnpdiffs (bool labelsp);
extern void
Substring_print_ncolordiffs ();


extern T
Substring_new (int nmismatches, int ncolordiffs, Chrnum_T chrnum, Genomicpos_T chroffset,
	       Genomicpos_T genomicstart, Genomicpos_T genomicend,
	       int querystart, int queryend, int querylength,
	       Genomicpos_T alignstart, Genomicpos_T alignend, int genomiclength,
	       int extraleft, int extraright, char *genomicseg, char *query,
	       bool plusp, bool trim_left_p, bool trim_right_p, int trim_maxlength,
	       int minlength, bool dibasep, bool cmetp);

extern void
Substring_free (T *old);

extern bool
Substring_plusp (T this);
extern int
Substring_nmismatches (T this);
extern int
Substring_nmatches (T this);
extern int
Substring_ncolordiffs (T this);
extern int
Substring_querystart (T this);
extern int
Substring_queryend (T this);
extern int
Substring_querylength (T this);
extern int
Substring_match_length (T this);

extern Chrnum_T
Substring_chrnum (T this);
extern Genomicpos_T
Substring_chroffset (T this);
extern Genomicpos_T
Substring_alignstart_trim (T this);
extern Genomicpos_T
Substring_alignend_trim (T this);
extern Genomicpos_T
Substring_genomicstart (T this);
extern Genomicpos_T
Substring_genomicend (T this);
extern Genomicpos_T
Substring_genomiclength (T this);

extern double
Substring_chimera_prob (T this);
extern int
Substring_chimera_pos (T this);
extern bool
Substring_chimera_knownp (T this);
extern bool
Substring_chimera_sensep (T this);


extern T
Substring_copy (T old);

extern T
Substring_new_donor (int donor_pos, int donor_nmismatches, int donor_ncolordiffs,
		     double donor_prob, Genomicpos_T left, int querylength,
		     bool plusp, bool sensep, char *genomicseg, char *query, Chrnum_T chrnum,
		     Genomicpos_T chroffset, bool knownp, int trim_maxlength, bool dibasep, bool cmetp);
extern T
Substring_new_acceptor (int acceptor_pos, int acceptor_nmismatches, int acceptor_ncolordiffs,
			double acceptor_prob, Genomicpos_T left, int querylength,
			bool plusp, bool sensep, char *genomicseg, char *query, Chrnum_T chrnum,
			Genomicpos_T chroffset, bool knownp, int trim_maxlength, bool dibasep, bool cmetp);

extern List_T
Substring_sort_chimera_halves (List_T hitlist, bool ascendingp);


extern void
Substring_print_single (T substring, Hittype_T hittype, Sequence_T queryseq,
			char *chr, int querylength,
			bool invertp, IIT_T snps_iit, int *snps_divint_crosstable);
extern void
Substring_print_insertion_1 (T substring1, T substring2, int nindels, 
			     Sequence_T queryseq, char *chr,
			     bool invertp, IIT_T snps_iit, int *snps_divint_crosstable);
extern void
Substring_print_insertion_2 (T substring1, T substring2, int nindels,
			     Sequence_T queryseq, char *chr,
			     bool invertp, IIT_T snps_iit, int *snps_divint_crosstable);
extern void
Substring_print_deletion_1 (T substring1, T substring2, int nindels, 
			    char *deletion, Sequence_T queryseq, char *chr,
			    bool invertp, IIT_T snps_iit, int *snps_divint_crosstable);
extern void
Substring_print_deletion_2 (T substring1, T substring2, int nindels, 
			    char *deletion, Sequence_T queryseq, char *chr,
			    bool invertp, IIT_T snps_iit, int *snps_divint_crosstable);
extern void
Substring_print_donor (T donor, bool sensep, bool invertp, Sequence_T queryseq,
		       IIT_T chromosome_iit, IIT_T snps_iit, int *snps_divint_crosstable, IIT_T splicesites_iit,
		       int donor_typeint, int *splicesites_divint_crosstable, bool endp,
		       T acceptor, Genomicpos_T chimera_distance);
extern void 
Substring_print_acceptor (T acceptor, bool sensep, bool invertp, Sequence_T queryseq,
			  IIT_T chromosome_iit, IIT_T snps_iit, int *snps_divint_crosstable, IIT_T splicesites_iit,
			  int acceptor_typeint, int *splicesites_divint_crosstable, bool endp,
			  T donor, Genomicpos_T chimera_distance);
extern void 
Substring_print_terminal (T terminal, bool invertp, Sequence_T queryseq, 
			  char *chr, IIT_T snps_iit, int *snps_divint_crosstable);

extern void
Substring_assign_donor_prob (T donor, Genome_T genome, IIT_T chromosome_iit);
extern void
Substring_assign_acceptor_prob (T acceptor, Genome_T genome, IIT_T chromosome_iit);


extern int
Substring_geneprob_eval_single (T substring, IIT_T geneprob_iit, IIT_T chromosome_iit);
extern int
Substring_geneprob_eval_insertion (T substring1, T substring2, int indel_pos, int nindels,
				   IIT_T geneprob_iit, IIT_T chromosome_iit, Sequence_T queryseq);
extern int
Substring_geneprob_eval_deletion (T substring1, T substring2, int indel_pos, int nindels,
				  IIT_T geneprob_iit, IIT_T chromosome_iit, Sequence_T queryseq);
extern int
Substring_geneprob_eval_splice (T donor, T  acceptor, IIT_T geneprob_iit, IIT_T chromosome_iit, int querylength,
				bool sensep);

#undef T
#endif


