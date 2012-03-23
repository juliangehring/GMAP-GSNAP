/* $Id: goby.h 58459 2012-02-24 21:02:35Z twu $ */
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

extern Gobyreader_T
Goby_reader_new (char **files, int nfiles, unsigned long window_start, unsigned long window_end, bool complement_reads_p);
extern void
Goby_reader_finish (Gobyreader_T reader);
extern void
Goby_reader_free (Gobyreader_T *old);
extern Shortread_T
Goby_read (Shortread_T *queryseq2, Gobyreader_T reader, int barcode_length,
	   bool invert_first_p, bool invert_second_p);

extern Gobywriter_T
Goby_writer_new (char *output_root, char *aligner_name, char *aligner_version);
extern void
Goby_writer_finish (Gobywriter_T writer, Gobyreader_T reader);
extern void
Goby_writer_free (Gobywriter_T *old);
extern void
Goby_writer_add_chromosomes (Gobywriter_T writer, IIT_T chromosome_iit);
extern void
Goby_file_handles (FILE **fp_capture, FILE**fp_ignore, Gobywriter_T writer);
extern void
Goby_start_capture (Gobywriter_T writer);
extern void
Goby_finish_capture (Gobywriter_T writer);
extern void
Goby_print_tmh (Gobywriter_T writer, Stage3end_T stage3, Shortread_T queryseq1, int npaths);
extern void
Goby_print_pair_tmh (Gobywriter_T writer, Resulttype_T resulttype, Stage3pair_T stage3pair, Shortread_T queryseq1, int npaths);

#endif

