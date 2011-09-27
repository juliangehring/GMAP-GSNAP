/* $Id: sam.h,v 1.3 2010/03/08 20:22:41 twu Exp $ */
#ifndef SAM_INCLUDED
#define SAM_INCLUDED

#include "stage3hr.h"
#include "iit-read.h"
#include "genome.h"
#include "sequence.h"
#include "resulthr.h"
#include "bool.h"

extern void
SAM_print_nomapping (Sequence_T queryseq, Stage3_T mate, char *acc,
		     IIT_T chromosome_iit, Resulttype_T resulttype, bool first_read_p,
		     int npaths_mate, Sequence_T queryseq_mate);

extern void
SAM_print (Stage3_T this, Stage3_T mate, char *acc, int pathnum, int npaths,
	   int score, Genome_T genome, IIT_T chromosome_iit, Sequence_T queryseq,
	   Sequence_T queryseq2, IIT_T snps_iit, int *snps_divint_crosstable,
	   IIT_T splicesites_iit, int donor_typeint, int acceptor_typeint,
	   int *splicesites_divint_crosstable, int pairedlength, Resulttype_T resulttype,
	   bool first_read_p, int npaths_mate);

extern void
SAM_print_paired (Result_T result, Genome_T genome, IIT_T chromosome_iit, Sequence_T queryseq1, Sequence_T queryseq2,
		  IIT_T snps_iit, int *snps_divint_crosstable,
		  IIT_T splicesites_iit, int donor_typeint, int acceptor_typeint,
		  int *splicesites_divint_crosstable, int maxpaths, bool quiet_if_excessive_p,
		  bool circularp);


#endif

