/* $Id: shortread.h 101822 2013-07-17 18:43:45Z twu $ */
#ifndef SHORTREAD_INCLUDED
#define SHORTREAD_INCLUDED
#include <stdio.h>
#include "bool.h"

#ifdef HAVE_ZLIB
#include <zlib.h>
#endif

#ifdef HAVE_BZLIB
#include "bzip2.h"
#endif


#define T Shortread_T
typedef struct T *T;

extern void
Shortread_setup (int acc_fieldi_start_in, int acc_fieldi_end_in,
		 bool force_singled_end_p_in, bool filter_chastity_p_in,
		 bool allow_paired_end_mismatch_p_in);

extern char *
Shortread_accession (T this);
extern char *
Shortread_header (T this);
extern bool
Shortread_filterp (T this);
extern bool
Shortread_invertedp (T this);

extern int
Shortread_input_init (FILE *fp);

#ifdef HAVE_ZLIB
extern int
Shortread_input_init_gzip (gzFile fp);
#endif

#ifdef HAVE_BZLIB
extern int
Shortread_input_init_bzip2 (Bzip2_T fp);
#endif

extern char *
Shortread_fullpointer (T this);
extern char *
Shortread_trimpointer (T this);

extern char *
Shortread_fullpointer_uc (T this);
extern char *
Shortread_contents_uc (T this);

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

extern bool
Shortread_chop_primers (T queryseq1, T queryseq2);
extern bool
Shortread_find_primers (T queryseq1, T queryseq2);
extern int
Shortread_max_overlap (T queryseq1, T queryseq2);
extern int
Shortread_find_overlap (T queryseq1, T queryseq2);

extern T
Shortread_new (char *acc, char *restofheader, bool filterp,
	       char *sequence, int sequence_length, char *quality, int quality_length,
	       int barcode_length, bool invertp, bool copy_acc_p, bool skipp);

extern T
Shortread_read_fasta_shortreads (int *nextchar, T *queryseq2, FILE **input1, FILE **input2,
				 char ***files, int *nfiles, bool skipp,
				 int barcode_length, bool invert_first_p, bool invert_second_p);
extern T
Shortread_read_fastq_shortreads (int *nextchar, T *queryseq2, FILE **input1, FILE **input2,
				 char ***files, int *nfiles, bool skipp,
				 int barcode_length, bool invert_first_p, bool invert_second_p);

#ifdef HAVE_ZLIB
extern T
Shortread_read_fasta_shortreads_gzip (int *nextchar, T *queryseq2, gzFile *input1, gzFile *input2,
				      char ***files, int *nfiles, bool skipp,
				      int barcode_length, bool invert_first_p, bool invert_second_p);
extern T
Shortread_read_fastq_shortreads_gzip (int *nextchar, T *queryseq2, gzFile *input1, gzFile *input2,
				      char ***files, int *nfiles, bool skipp,
				      int barcode_length, bool invert_first_p, bool invert_second_p);
#endif

#ifdef HAVE_BZLIB
extern T
Shortread_read_fasta_shortreads_bzip2 (int *nextchar, T *queryseq2, Bzip2_T *input1, Bzip2_T *input2,
				       char ***files, int *nfiles, bool skipp,
				       int barcode_length, bool invert_first_p, bool invert_second_p);
extern T
Shortread_read_fastq_shortreads_bzip2 (int *nextchar, T *queryseq2, Bzip2_T *input1, Bzip2_T *input2,
				       char ***files, int *nfiles, bool skipp,
				       int barcode_length, bool invert_first_p, bool invert_second_p);
#endif


extern void
Shortread_print_header (FILE *fp, T queryseq1, T queryseq2);

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
Shortread_print_barcode (FILE *fp, T this);
extern void
Shortread_print_chop (FILE *fp, T this, bool invertp);
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
