/* $Id: shortread.h 36020 2011-03-03 17:06:19Z twu $ */
#ifndef SHORTREAD_INCLUDED
#define SHORTREAD_INCLUDED
#include <stdio.h>
#include "bool.h"

#ifdef HAVE_ZLIB
#include <zlib.h>
#endif


#define T Shortread_T
typedef struct T *T;

extern char *
Shortread_accession (T this);

extern int
Shortread_input_init (FILE *fp);

#ifdef HAVE_ZLIB
extern int
Shortread_input_init_gzip (gzFile fp);
#endif

extern char *
Shortread_fullpointer (T this);
extern char *
Shortread_trimpointer (T this);

extern char *
Shortread_fullpointer_uc (T this);

extern int
Shortread_barcode_length (T this);
extern char *
Shortread_barcode (T this);

extern int
Shortread_choplength (T this);
extern char *
Shortread_quality_string (T this);

extern int
Shortread_fulllength (T this);

extern void
Shortread_free (T *old);

extern void
Shortread_chop_primers (T queryseq1, T queryseq2);

extern T
Shortread_new (char *acc, char *restofheader,
	       char *sequence, int sequence_length, char *quality, int quality_length,
	       int barcode_length, bool invertp, bool copy_acc_p);

extern T
Shortread_read_fasta_shortreads (int *nextchar, T *queryseq2, FILE **input, char ***files, int *nfiles,
				 int barcode_length, bool invert_first_p, bool invert_second_p,
				 bool pc_linefeeds_p);
extern T
Shortread_read_fastq_shortreads (int *nextchar, T *queryseq2, FILE *input1, FILE *input2,
				 int barcode_length, bool invert_first_p, bool invert_second_p,
				 bool pc_linefeeds_p);

#ifdef HAVE_ZLIB
extern T
Shortread_read_fasta_shortreads_gzip (int *nextchar, T *queryseq2, gzFile input, char ***files, int *nfiles,
				      int barcode_length, bool invert_first_p, bool invert_second_p,
				      bool pc_linefeeds_p);
extern T
Shortread_read_fastq_shortreads_gzip (int *nextchar, T *queryseq2, gzFile gzipped, gzFile gzipped2,
				      int barcode_length, bool invert_first_p, bool invert_second_p,
				      bool pc_linefeeds_p);
#endif


extern void
Shortread_print_header (FILE *fp, T this);

extern void
Shortread_print_query_singleend_fasta (FILE *fp, T queryseq);
extern void
Shortread_print_query_singleend_fastq (FILE *fp, T queryseq);
extern void
Shortread_print_query_pairedend_fasta (FILE *fp, T queryseq1, T queryseq2,
				      bool invert_first_p, bool invert_second_p);
extern void
Shortread_print_query_pairedend_fastq (FILE *fp1, FILE *fp2, T queryseq1, T queryseq2,
				      bool invert_first_p, bool invert_second_p);

extern void
Shortread_print_oneline (FILE *fp, T this);
extern void
Shortread_print_oneline_revcomp (FILE *fp, T this);

extern void
Shortread_print_chopped (FILE *fp, T this, int hardclip_low, int hardclip_high);
extern void
Shortread_print_chopped_revcomp (FILE *fp, T this, int hardclip_low, int hardclip_high);

extern void
Shortread_print_chop_symbols (FILE *fp, T this);
extern void
Shortread_print_quality (FILE *fp, T this, int hardclip_low, int hardclip_high,
			int shift, bool show_chopped_p);
extern void
Shortread_print_quality_revcomp (FILE *fp, T this, int hardclip_low, int hardclip_high,
				int shift, bool show_chopped_p);
extern void
Shortread_print_oneline_uc (FILE *fp, T this);
extern void
Shortread_print_oneline_revcomp_uc (FILE *fp, T this);

#undef T
#endif
