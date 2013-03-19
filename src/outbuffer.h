/* $Id: outbuffer.h 87096 2013-02-22 21:04:02Z twu $ */
#ifndef OUTBUFFER_INCLUDED
#define OUTBUFFER_INCLUDED

#include "types.h"
#include "bool.h"
#include "genomicpos.h"
#include "sequence.h"
#include "iit-read.h"

#include "request.h"
#include "mem.h"		/* To get MEMUSAGE */

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
Outbuffer_new (unsigned int output_buffer_size, unsigned int nread, char *sevenway_root, bool appendp, IIT_T chromosome_iit,
	       bool timingp, bool output_sam_p, bool sam_headers_p, char *sam_read_group_id, char *sam_read_group_name,
	       char *sam_read_group_library, char *sam_read_group_platform,
	       Gobywriter_T gobywriter, bool nofailsp, bool failsonlyp, bool fails_as_input_p,
	       bool fastq_format_p, bool clip_overlap_p, bool merge_samechr_p,
	       int maxpaths_report, bool quiet_if_excessive_p, int quality_shift,
	       bool invert_first_p, bool invert_second_p, Genomicpos_T pairmax);

#else

extern T
Outbuffer_new (unsigned int output_buffer_size, unsigned int nread, char *sevenway_root, bool appendp,
	       bool chimeras_allowed_p, char *user_genomicseg, Sequence_T usersegment,
	       char *dbversion, Genome_T genome, IIT_T chromosome_iit,
	       Chrsubset_T chrsubset, IIT_T contig_iit, IIT_T altstrain_iit, IIT_T map_iit,
	       int *map_divint_crosstable, Printtype_T printtype, bool checksump, int chimera_margin,
#ifndef PMAP
	       bool sam_headers_p, int quality_shift, bool sam_paired_p,
	       char *sam_read_group_id, char *sam_read_group_name,
	       char *sam_read_group_library, char *sam_read_group_platform,
#endif
	       bool nofailsp, bool failsonlyp, bool fails_as_input_p, int maxpaths_report, bool quiet_if_excessive_p,
	       bool map_exons_p, bool map_bothstrands_p, bool print_comment_p, int nflanking,
	       int proteinmode, int invertmode, bool nointronlenp, int wraplength,
	       int ngap, int cds_startpos,
	       bool fulllengthp, bool truncatep, bool strictp, bool diagnosticp, bool maponlyp,
	       bool stage1debug, bool diag_debug, bool debug_graphic_p,
	       int argc, char **argv, int optind);

#endif

extern void
Outbuffer_free (T *old);

extern unsigned int
Outbuffer_nread (T this);

extern void
Outbuffer_add_nread (T this, unsigned int nread);

extern void
Outbuffer_put_result (T this, Result_T result, Request_T request);

#ifdef GSNAP
extern void
Outbuffer_print_result (T this, Result_T result, Request_T request
#ifdef MEMUSAGE
			, unsigned int noutput
#endif
			);
#else
extern void
Outbuffer_print_result (T this, Result_T result, Request_T request, Sequence_T headerseq
#ifdef MEMUSAGE
			, unsigned int noutput
#endif
			);
#endif

extern void *
Outbuffer_thread_anyorder (void *data);

extern void *
Outbuffer_thread_ordered (void *data);

#undef T
#endif

