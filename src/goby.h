/* $Id: goby.h 42526 2011-07-08 22:08:21Z twu $ */
#ifndef GOBY_INCLUDED
#define GOBY_INCLUDED

#include <stdio.h>
#include "bool.h"
#include "iit-read.h"
#include "shortread.h"
#include "stage3hr.h"

typedef struct Gobyreader_T *Gobyreader_T;
typedef struct Gobywriter_T *Gobywriter_T;

extern void
Goby_setup (bool show_refdiff_p_in);
extern void
Goby_shutdown ();
extern void
Goby_reader_free (Gobyreader_T *old);
extern Gobyreader_T
Goby_reader_new (char **files, int nfiles, unsigned long window_start, unsigned long window_end, bool complement_reads_p);
extern Shortread_T
Goby_read (Shortread_T *queryseq2, Gobyreader_T reader, int barcode_length,
	   bool invert_first_p, bool invert_second_p);
extern void
Goby_reader_finish (Gobyreader_T reader);

extern void
Goby_writer_free (Gobywriter_T *old);
extern Gobywriter_T
Goby_writer_new (char *output_root, char *aligner_name, char *aligner_version);
extern void
Goby_writer_add_chromosomes (Gobywriter_T writer, IIT_T chromosome_iit);
extern void
Goby_observe_aligned(Gobywriter_T writer);
extern void
Goby_print_tmh (Gobywriter_T writer, Stage3end_T stage3, Shortread_T queryseq1, int npaths);
void
Goby_print_single (Gobywriter_T writer, Stage3end_T this, int score,
		   IIT_T chromosome_iit, Shortread_T queryseq,
		   bool invertp, Stage3end_T hit5, Stage3end_T hit3,
		   int insertlength, int pairscore, Pairtype_T pairtype,
		   int mapq_score);
extern void
Goby_print_paired (Gobywriter_T writer, Result_T result, Resulttype_T resulttype,
		   IIT_T chromosome_iit, Shortread_T queryseq1, Shortread_T queryseq2,
		   int maxpaths, bool quiet_if_excessive_p,
		   bool invert_first_p, bool invert_second_p,
		   bool nofailsp, bool failsonlyp, bool fails_as_input_p,
		   bool fastq_format_p, int quality_shift, char *sam_read_group_id);
extern void
Goby_writer_finish (Gobywriter_T writer, Gobyreader_T reader);

#endif

