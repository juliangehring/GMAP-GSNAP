/* $Id: samprint.h 37254 2011-03-28 16:34:08Z twu $ */
#ifndef SAMPRINT_INCLUDED
#define SAMPRINT_INCLUDED

#include <stdio.h>
#include "stage3hr.h"
#include "iit-read.h"
#include "shortread.h"
#include "resulthr.h"
#include "bool.h"

extern unsigned int
SAM_compute_flag (Substring_T substring, Stage3_T mate, Resulttype_T resulttype,
		  bool first_read_p, int pathnum, int npaths, int npaths_mate,
		  bool invertp, bool invert_mate_p);

extern void
SAM_print_nomapping (FILE *fp, Shortread_T queryseq, Stage3_T mate, char *acc,
		     IIT_T chromosome_iit, Resulttype_T resulttype, bool translocationp, bool first_read_p,
		     int npaths_mate, Shortread_T queryseq_mate, int quality_shift,
		     char *sam_read_group_id, bool invertp, bool invert_mate_p);

extern void
SAM_print (FILE *fp, Stage3_T this, Stage3_T mate, char *acc, int pathnum, int npaths,
	   int score, IIT_T chromosome_iit, Shortread_T queryseq,
	   Shortread_T queryseq2, int pairedlength, Resulttype_T resulttype, bool translocationp,
	   bool first_read_p, int npaths_mate, int quality_shift,
	   char *sam_read_group_id, bool invertp, bool invert_mate_p);

extern void
SAM_print_paired (Result_T result, Resulttype_T resulttype, bool translocationp,
		  IIT_T chromosome_iit, Shortread_T queryseq1, Shortread_T queryseq2,
		  int maxpaths, bool quiet_if_excessive_p,
		  bool invert_first_p, bool invert_second_p,
		  bool nofailsp, bool failsonlyp, bool fails_as_input_p,
		  bool fastq_format_p, int quality_shift, char *sam_read_group_id,
		  FILE *fp_nomapping_1, FILE *fp_nomapping_2,
		  FILE *fp_unpaired_uniq, FILE *fp_unpaired_mult,
		  FILE *fp_halfmapping_uniq, FILE *fp_halfmapping_mult,
		  FILE *fp_paired_uniq_inv, FILE *fp_paired_uniq_scr, FILE *fp_paired_long,
		  FILE *fp_paired_mult, FILE *fp_concordant_uniq, FILE *fp_concordant_mult);

#endif

