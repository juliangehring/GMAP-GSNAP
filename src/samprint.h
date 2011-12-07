/* $Id: samprint.h 53168 2011-11-28 02:55:16Z twu $ */
#ifndef SAMPRINT_INCLUDED
#define SAMPRINT_INCLUDED

#include <stdio.h>
#include "stage3hr.h"
#include "iit-read.h"
#include "shortread.h"
#include "resulthr.h"
#include "genomicpos.h"
#include "substring.h"
#include "bool.h"

extern void
SAM_setup (bool quiet_if_excessive_p_in, int maxpaths_in);

extern Genomicpos_T
SAM_compute_chrpos (int *hardclip_low, int *hardclip_high, Stage3end_T this,
		    Substring_T substring_low, int querylength);

extern unsigned int
SAM_compute_flag (bool plusp, Stage3end_T mate, Resulttype_T resulttype,
		  bool first_read_p, int pathnum, int npaths, int npaths_mate,
		  bool invertp, bool invert_mate_p);

extern void
SAM_print_nomapping (FILE *fp, Shortread_T queryseq, Stage3end_T mate, char *acc,
		     IIT_T chromosome_iit, Resulttype_T resulttype, bool first_read_p,
		     int npaths_mate, Genomicpos_T mate_chrpos,
		     int quality_shift, char *sam_read_group_id, bool invertp, bool invert_mate_p);

extern void
SAM_print (FILE *fp, Stage3end_T this, Stage3end_T mate, char *acc, int pathnum, int npaths,
	   int absmq_score, int second_absmq, int mapq_score, IIT_T chromosome_iit, Shortread_T queryseq,
	   Shortread_T queryseq2, int pairedlength, Genomicpos_T chrpos, Genomicpos_T mate_chrpos,
	   int hardclip_low, int hardclip_high, Resulttype_T resulttype, bool first_read_p,
	   int npaths_mate, int quality_shift, char *sam_read_group_id, bool invertp, bool invert_mate_p,
	   bool merge_samechr_p);

extern void
SAM_print_paired (Result_T result, Resulttype_T resulttype,
		  IIT_T chromosome_iit, Shortread_T queryseq1, Shortread_T queryseq2,
		  bool invert_first_p, bool invert_second_p,
		  bool nofailsp, bool failsonlyp, bool fails_as_input_p,
		  bool fastq_format_p, bool clip_overlap_p, bool merge_samechr_p,
		  int quality_shift, char *sam_read_group_id,
		  FILE *fp_nomapping_1, FILE *fp_nomapping_2,
		  FILE *fp_unpaired_uniq, FILE *fp_unpaired_transloc, FILE *fp_unpaired_mult,
		  FILE *fp_halfmapping_uniq, FILE *fp_halfmapping_transloc, FILE *fp_halfmapping_mult,
		  FILE *fp_paired_uniq_inv, FILE *fp_paired_uniq_scr, FILE *fp_paired_long,
		  FILE *fp_paired_mult, FILE *fp_concordant_uniq, FILE *fp_concordant_transloc,
		  FILE *fp_concordant_mult);

#endif

