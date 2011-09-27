/* $Id: outbuffer.h 36834 2011-03-19 01:59:02Z twu $ */
#ifndef OUTBUFFER_INCLUDED
#define OUTBUFFER_INCLUDED

#include "types.h"
#include "bool.h"
#include "sequence.h"
#include "iit-read.h"

#include "request.h"

#ifdef GSNAP
#include "goby.h"
#include "resulthr.h"

#else
#include "stage3.h"		/* Has Printtype_T */
#include "result.h"
#include "chrsubset.h"
#include "genome.h"

#endif


#define T Outbuffer_T
typedef struct T *T;

#ifdef GSNAP

extern T
Outbuffer_new (int nread, char *sevenway_root, UINT4 *genome_blocks, IIT_T chromosome_iit,
	       bool output_sam_p, bool sam_headers_p, char *sam_read_group_id, char *sam_read_group_name,
	       Gobywriter_T gobywriter, bool nofailsp, bool failsonlyp, bool fails_as_input_p,
	       bool fastq_format_p, int maxpaths, bool quiet_if_excessive_p, int quality_shift,
	       bool invert_first_p, bool invert_second_p, Genomicpos_T pairmax);

#else

extern T
Outbuffer_new (int nread, char *sevenway_root, char *user_genomicseg, Sequence_T usersegment,
	       char *dbversion, Genome_T genome, IIT_T chromosome_iit,
	       Chrsubset_T chrsubset, IIT_T contig_iit, IIT_T altstrain_iit, IIT_T map_iit,
	       int *map_divint_crosstable, Printtype_T printtype, bool checksump, int chimera_margin,
#ifndef PMAP
	       bool sam_headers_p, int quality_shift, bool sam_paired_p, bool cigar_noncanonical_splices_p,
	       char *sam_read_group_id, char *sam_read_group_name,
#endif
	       bool nofailsp, bool failsonlyp, bool fails_as_input_p, int maxpaths, bool quiet_if_excessive_p,
	       bool map_exons_p, bool map_bothstrands_p, bool print_comment_p, int nflanking,
	       int proteinmode, int invertmode, bool nointronlenp, int wraplength, int ngap, int cds_startpos,
	       bool fulllengthp, bool truncatep, bool strictp, bool diagnosticp, bool maponlyp,
	       bool stage1debug, bool diag_debug, bool debug_graphic_p,
	       int argc, char **argv, int optind);

#endif

extern void
Outbuffer_free (T *old);

extern int
Outbuffer_nread (T this);

extern void
Outbuffer_add_nread (T this, int nread);

extern void
Outbuffer_put_result (T this, Result_T result, Request_T request);

extern void
Outbuffer_print_result (T this, Result_T result, Request_T request);

extern void *
Outbuffer_thread_anyorder (void *data);

extern void *
Outbuffer_thread_ordered (void *data);

#undef T
#endif

