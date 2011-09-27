/* $Id: sam.h,v 1.7 2010-07-22 00:11:06 twu Exp $ */
#ifndef SAM_INCLUDED
#define SAM_INCLUDED

#include "stage3hr.h"
#include "iit-read.h"
#include "sequence.h"
#include "resulthr.h"
#include "bool.h"

#define PAIRED_READ        0x0001
#define PAIRED_MAPPING     0x0002
#define QUERY_UNMAPPED     0x0004
#define MATE_UNMAPPED      0x0008
#define QUERY_MINUSP       0x0010
#define MATE_MINUSP        0x0020
#define FIRST_READ_P       0x0040
#define SECOND_READ_P      0x0080
#define NOT_PRIMARY        0x0100
#define BAD_READ_QUALITY   0x0200
#define DUPLICATE_READ     0x0400


extern void
SAM_print_nomapping (Sequence_T queryseq, Stage3_T mate, char *acc,
		     IIT_T chromosome_iit, Resulttype_T resulttype, bool first_read_p,
		     int npaths_mate, Sequence_T queryseq_mate, int quality_shift);

extern void
SAM_print (Stage3_T this, Stage3_T mate, char *acc, int pathnum, int npaths,
	   int score, IIT_T chromosome_iit, Sequence_T queryseq,
	   Sequence_T queryseq2, IIT_T snps_iit, int *snps_divint_crosstable,
	   IIT_T splicesites_iit, int donor_typeint, int acceptor_typeint,
	   int *splicesites_divint_crosstable, int pairedlength, Resulttype_T resulttype,
	   bool first_read_p, int npaths_mate, int quality_shift);

extern void
SAM_print_paired (Result_T result, IIT_T chromosome_iit, Sequence_T queryseq1, Sequence_T queryseq2,
		  IIT_T snps_iit, int *snps_divint_crosstable,
		  IIT_T splicesites_iit, int donor_typeint, int acceptor_typeint,
		  int *splicesites_divint_crosstable, int maxpaths, bool quiet_if_excessive_p,
		  int quality_shift, bool circularp);


#endif

