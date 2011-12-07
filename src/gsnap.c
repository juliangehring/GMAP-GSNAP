static char rcsid[] = "$Id: gsnap.c 53585 2011-12-02 18:55:24Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>		/* Needed to define pthread_t on Solaris */
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>		/* For strcpy */
#include <strings.h>		/* For rindex */
#include <ctype.h>
#include <math.h>		/* For rint */

#ifdef HAVE_PTHREAD
#include <pthread.h>
#endif

#ifdef HAVE_ZLIB
#include <zlib.h>
#define GZBUFFER_SIZE 131072
#endif

#include <signal.h>

#include "assert.h"
#include "except.h"
#include "mem.h"
#include "bool.h"
#include "types.h"
#include "fopen.h"

#include "mode.h"
#include "sequence.h"
#include "shortread.h"		/* For Shortread_setup */
#include "stopwatch.h"
#include "genome.h"
#include "genome_hr.h"		/* For Genome_hr_setup */
#include "maxent_hr.h"		/* For Maxent_hr_setup */
#include "indexdb_hr.h"
#include "mapq.h"
#include "substring.h"
#include "stage3hr.h"
#include "goby.h"
#include "spanningelt.h"
#include "splicetrie_build.h"
#include "oligo.h"		/* For Oligo_setup */
#include "oligoindex_hr.h"	/* For Oligoindex_hr_setup */
#include "stage2.h"		/* For Stage2_setup */
#include "stage1hr.h"
#include "indexdb.h"
#include "resulthr.h"
#include "request.h"
#include "intlist.h"
#include "list.h"
#include "listdef.h"
#include "iit-read.h"
#include "datadir.h"
#include "inbuffer.h"
#include "outbuffer.h"
#include "samprint.h"		/* For SAM_setup */

#include "stage3.h"		/* To get EXTRAQUERYGAP */

#include "getopt.h"


#define MIN_INDEXDB_SIZE_THRESHOLD 100


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/************************************************************************
 *   GMAP parameters
 ************************************************************************/

static int gmap_mode = GMAP_PAIRSEARCH | GMAP_TERMINAL | GMAP_IMPROVEMENT;
static int nullgap = 600;
static int maxpeelback = 11;
static int maxpeelback_distalmedial = 24;
static int extramaterial_end = 10;
static int extramaterial_paired = 8;
static int max_gmap_pairsearch = 10;
static int max_gmap_terminal = 5;
static int max_gmap_improvement = 5;

static double microexon_spliceprob = 0.95;
static int suboptimal_score_start = -1; /* Determined by simulations to have minimal effect */
static int suboptimal_score_end = 3; /* Determined by simulations to have diminishing returns above 3 */

static int trigger_score_for_gmap = 5;
/* static int trigger_score_for_terminals = 5; -- obsolete */

static int min_intronlength = 9;
static int max_deletionlength = 50;


/************************************************************************
 *   Global parameters
 ************************************************************************/

static IIT_T chromosome_iit = NULL;
static int nchromosomes = 0;
static Indexdb_T indexdb = NULL;
static Indexdb_T indexdb2 = NULL; /* For cmet or atoi */
static Genome_T genome = NULL;
static Genome_T genomealt = NULL;
static UINT4 *genome_blocks = NULL;

static bool fastq_format_p = false;
static bool creads_format_p = false;
static Stopwatch_T stopwatch = NULL;

/************************************************************************
 *   Program options
 ************************************************************************/

/* Input options */
static char *user_genomedir = NULL;
static char *dbroot = NULL;
static char *dbversion = NULL;
static int part_modulus = 0;
static int part_interval = 1;
static int barcode_length = 0;
static bool invert_first_p = false;
static bool invert_second_p = true;
static int acc_fieldi_start = 0;
static int acc_fieldi_end = 0;
static bool filter_chastity_p = false;
static bool filter_if_both_p = false;
static bool gunzip_p = false;

/* Compute options */
static bool chop_primers_p = false;
static bool query_unk_mismatch_p = false;
static bool genome_unk_mismatch_p = true;
static bool novelsplicingp = false;

static int trim_mismatch_score = -3;
static int trim_indel_score = -4;


static Access_mode_T offsetscomp_access = USE_ALLOCATE;
static bool expand_offsets_p = false;

#ifdef HAVE_MMAP
static Access_mode_T positions_access = USE_MMAP_PRELOAD;
static Access_mode_T genome_access = USE_MMAP_PRELOAD;
#else
static Access_mode_T positions_access = USE_ALLOCATE;
static Access_mode_T genome_access = USE_ALLOCATE;
#endif

static int pairmax;
static int pairmax_dna = 1000;
static int pairmax_rna = 200000;

static Genomicpos_T expected_pairlength = 200;
static Genomicpos_T pairlength_deviation = 25;

#ifdef HAVE_PTHREAD
static pthread_t output_thread_id, *worker_thread_ids;
static pthread_key_t global_request_key;
static int nworkers = 1;	/* (int) sysconf(_SC_NPROCESSORS_ONLN) */
#endif

/* static Masktype_T masktype = MASK_REPETITIVE; */
static int subopt_levels = 0;

/* If negative, then hasn't been specified by user.  If between 0 and
   1, then treated as a fraction of the querylength.  Else, treated as
   an integer */
static double user_maxlevel_float = -1.0;

static int terminal_threshold = 2;

/* Really have only one indel penalty */
static int indel_penalty_middle = 2;
static int indel_penalty_end = 2;

static bool allow_end_indels_p = true;
static int max_middle_insertions = 9;
static int max_middle_deletions = 30;
static int max_end_insertions = 3;
static int max_end_deletions = 6;
static int min_indel_end_matches = 4;
static Genomicpos_T shortsplicedist = 200000;
static Genomicpos_T shortsplicedist_known;
static int localsplicing_penalty = 0;
static int distantsplicing_penalty = 3;
static int min_distantsplicing_end_matches = 16;
static double min_distantsplicing_identity = 0.95;
static int min_shortend = 2;
/* static bool find_novel_doublesplices_p = true; */
static int antistranded_penalty = 0; /* Most RNA-Seq is non-stranded */

static int index1part;
static int required_index1part = 0;
static int indexdb_size_threshold;


/* Genes IIT */
static char *genes_file = (char *) NULL;
static IIT_T genes_iit = NULL;
static int *genes_divint_crosstable = NULL;
static bool favor_multiexon_p = false;


/* Splicing IIT */
static bool knownsplicingp = false;
static bool distances_observed_p = false;
static char *user_splicingdir = (char *) NULL;
static char *splicing_file = (char *) NULL;
static IIT_T splicing_iit = NULL;
static bool amb_closest_p = false;
static bool amb_clip_p = true;

static int donor_typeint = -1;		/* for splicing_iit */
static int acceptor_typeint = -1;	/* for splicing_iit */

static int *splicing_divint_crosstable = NULL;
static UINT4 *splicecomp = NULL;
static Genomicpos_T *splicesites = NULL;
static Splicetype_T *splicetypes = NULL;
static Genomicpos_T *splicedists = NULL; /* maximum observed splice distance for given splice site */
static List_T *splicestrings = NULL;
static UINT4 *splicefrags_ref = NULL;
static UINT4 *splicefrags_alt = NULL;
static int nsplicesites = 0;


/* Splicing via splicesites */
static int *nsplicepartners_skip = NULL;
static int *nsplicepartners_obs = NULL;
static int *nsplicepartners_max = NULL;

static bool splicetrie_precompute_p = true;
static unsigned int *trieoffsets_obs = NULL;
static unsigned int *triecontents_obs = NULL;
static unsigned int *trieoffsets_max = NULL;
static unsigned int *triecontents_max = NULL;


/* Dibase and CMET */
static bool dibasep = false;
static char *user_cmetdir = NULL;
static char *user_atoidir = NULL;
static Mode_T mode = STANDARD;

/* SNPs IIT */
static char *user_snpsdir = NULL;
static char *snps_root = (char *) NULL;
static IIT_T snps_iit = NULL;
static int *snps_divint_crosstable = NULL;


/* Tally IIT */
static char *user_tallydir = NULL;
static char *tally_root = (char *) NULL;
static IIT_T tally_iit = NULL;
static int *tally_divint_crosstable = NULL;

static char *user_runlengthdir = NULL;
static char *runlength_root = (char *) NULL;
static IIT_T runlength_iit = NULL;
static int *runlength_divint_crosstable = NULL;


/* Output options */
static unsigned int output_buffer_size = 1000;
static bool output_sam_p = false;
static bool output_goby_p = false;

/* For Illumina, subtract 64.  For Sanger, subtract 33.  For Goby, subtract 0. */
/* static int quality_score_adj = 64;  -- Stored in mapq.c */

static bool user_quality_score_adj = false;
static bool user_quality_shift = false;
static int quality_shift = 0;   /* For printing, may want -31 */

static bool exception_raise_p = true;
static bool quiet_if_excessive_p = false;
static int maxpaths = 100;
static bool orderedp = false;
static bool failsonlyp = false;
static bool nofailsp = false;
static bool fails_as_input_p = false;

static bool print_ncolordiffs_p = false;
static bool print_nsnpdiffs_p = false;
static bool print_snplabels_p = false;

static bool show_refdiff_p = false;
static bool clip_overlap_p = false;
static bool merge_samechr_p = false;

/* SAM */
static int sam_headers_batch = -1;
static bool sam_headers_p = true;
static char *sam_read_group_id = NULL;
static char *sam_read_group_name = NULL;
static char *sam_read_group_library = NULL;
static char *sam_read_group_platform = NULL;


/* Goby */
static char *goby_output_root = NULL;
static unsigned long creads_window_start = 0;
static unsigned long creads_window_end = 0;
static bool creads_complement_p = false;
static Gobyreader_T gobyreader = NULL;
static Gobywriter_T gobywriter = NULL;

/* Input/output */
static char *sevenway_root = NULL;
static Outbuffer_T outbuffer;
static Inbuffer_T inbuffer;
static unsigned int inbuffer_nspaces = 1000;
static unsigned int inbuffer_maxchars = -1U; /* Currently not used by Inbuffer_T */
static bool timingp = false;
static bool unloadp = false;


/* Alignment options */
static bool uncompressedp = false;

/* getopt used alphabetically: AaBDdEeGgiJjKklMmNnOoQqstVvwYyZz7 */

static struct option long_options[] = {
  /* Input options */
  {"dir", required_argument, 0, 'D'},	/* user_genomedir */
  {"db", required_argument, 0, 'd'}, /* dbroot */
  {"kmer", required_argument, 0, 'k'}, /* required_index1part, index1part */
  {"genomefull", no_argument, 0, 'G'}, /* uncompressedp */
  {"part", required_argument, 0, 'q'}, /* part_modulus, part_interval */
  {"orientation", required_argument, 0, 'o'}, /* invert_first_p, invert_second_p */
  {"input-buffer-size", required_argument, 0, 0}, /* inbuffer_nspaces */
  {"barcode-length", required_argument, 0, 0},	  /* barcode_length */
  {"fastq-id-start", required_argument, 0, 0},	  /* acc_fieldi_start */
  {"fastq-id-end", required_argument, 0, 0},	  /* acc_fieldi_end */
  {"filter-chastity", required_argument, 0, 0},	/* filter_chastity_p, filter_if_both_p */

#ifdef HAVE_ZLIB
  {"gunzip", no_argument, 0, 0}, /* gunzip_p */
#endif

  /* Compute options */
#ifdef HAVE_MMAP
  {"batch", required_argument, 0, 'B'}, /* offsetscomp_access, positions_access, genome_access */
#endif
  {"pairmax-dna", required_argument, 0, 0}, /* pairmax_dna */
  {"pairmax-rna", required_argument, 0, 0}, /* pairmax_rna */
  {"pairexpect", required_argument, 0, 0}, /* expected_pairlength */
  {"pairdev", required_argument, 0, 0}, /* pairlength_deviation */
#ifdef HAVE_PTHREAD
  {"nthreads", required_argument, 0, 't'}, /* nworkers */
#endif
  {"adapter-strip", required_argument, 0, 'a'},	/* chop_primers_p */

  {"query-unk-mismatch", required_argument, 0, 0}, /* query_unk_mismatch_p */
  {"genome-unk-mismatch", required_argument, 0, 0}, /* genome_unk_mismatch_p */

  {"trim-mismatch-score", required_argument, 0, 0}, /* trim_mismatch_score */
  {"trim-indel-score", required_argument, 0, 0}, /* trim_indel_score */
  {"novelsplicing", required_argument, 0, 'N'}, /* novelsplicingp */

  {"max-mismatches", required_argument, 0, 'm'}, /* user_maxlevel_float */
  {"terminal-threshold", required_argument, 0, 0}, /* terminal_threshold */

#if 0
  {"indel-penalty-middle", required_argument, 0, 'i'}, /* indel_penalty_middle */
  {"indel-penalty-end", required_argument, 0, 'I'}, /* indel_penalty_end */
#else
  {"indel-penalty", required_argument, 0, 'i'}, /* indel_penalty_middle, indel_penalty_end */
#endif

  {"indel-endlength", required_argument, 0, 0}, /* min_indel_end_matches, allow_end_indels_p */

  {"max-middle-insertions", required_argument, 0, 'y'}, /* max_middle_insertions */
  {"max-middle-deletions", required_argument, 0, 'z'}, /* max_middle_deletions */
  {"max-end-insertions", required_argument, 0, 'Y'}, /* max_end_insertions */
  {"max-end-deletions", required_argument, 0, 'Z'}, /* max_end_deletions */
  {"suboptimal-levels", required_argument, 0, 'M'}, /* subopt_levels */

  {"localsplicedist", required_argument, 0, 'w'}, /* shortsplicedist */
  {"splicingdir", required_argument, 0, 0},	  /* user_splicingdir */
  {"use-splicing", required_argument, 0, 's'}, /* splicing_iit, knownsplicingp */
  {"ambig-splice-noclip", no_argument, 0, 0},  /* amb_clip_p */
  {"genes", required_argument, 0, 'g'}, /* genes_iit */
  {"favor-multiexon", no_argument, 0, 0}, /* favor_multiexon_p */

  {"local-splice-penalty", required_argument, 0, 'e'}, /* localsplicing_penalty */
  {"distant-splice-penalty", required_argument, 0, 'E'}, /* distantsplicing_penalty */
  {"distant-splice-endlength", required_argument, 0, 'K'}, /* min_distantsplicing_end_matches */
  {"shortend-splice-endlength", required_argument, 0, 'l'}, /* min_shortend */
  {"antistranded-penalty", required_argument, 0, 0},	    /* antistranded_penalty */
  {"merge-distant-samechr", no_argument, 0, 0},		    /* merge_samechr_p */

  {"cmetdir", required_argument, 0, 0}, /* user_cmetdir */
  {"atoidir", required_argument, 0, 0}, /* user_atoidir */
  {"mode", required_argument, 0, 0}, /* mode */

  {"snpsdir", required_argument, 0, 'V'},   /* user_snpsdir */
  {"use-snps", required_argument, 0, 'v'}, /* snps_root */

  {"tallydir", required_argument, 0, 0},   /* user_tallydir */
  {"use-tally", required_argument, 0, 0}, /* tally_root */

  {"runlengthdir", required_argument, 0, 0},   /* user_runlengthdir */
  {"use-runlength", required_argument, 0, 0}, /* runlength_root */

  {"gmap-mode", required_argument, 0, 0}, /* gmap_mode */
  {"trigger-score-for-gmap", required_argument, 0, 0}, /* trigger_score_for_gmap */
  {"max-gmap-pairsearch", required_argument, 0, 0}, /* max_gmap_pairsearch */
  {"max-gmap-terminal", required_argument, 0, 0}, /* max_gmap_terminal */
  {"max-gmap-improvement", required_argument, 0, 0}, /* max_gmap_improvement */
  {"microexon-spliceprob", required_argument, 0, 0}, /* microexon_spliceprob */
  {"stage2-start", required_argument, 0, 0},	     /* suboptimal_score_start */
  {"stage2-end", required_argument, 0, 0},	     /* suboptimal_score_end */

  /* Output options */
  {"output-buffer-size", required_argument, 0, 0}, /* output_buffer_size */
  {"format", required_argument, 0, 'A'}, /* output_sam_p, output_goby_p */

  {"quality-protocol", required_argument, 0, 0}, /* quality_score_adj, quality_shift */
  {"quality-zero-score", required_argument, 0, 'J'}, /* quality_score_adj */
  {"quality-print-shift", required_argument, 0, 'j'}, /* quality_shift */
  {"sam-headers-batch", required_argument, 0, 0},	/* sam_headers_batch */
  {"no-sam-headers", no_argument, 0, 0},	/* sam_headers_p */
  {"read-group-id", required_argument, 0, 0},	/* sam_read_group_id */
  {"read-group-name", required_argument, 0, 0},	/* sam_read_group_name */
  {"read-group-library", required_argument, 0, 0},	/* sam_read_group_library */
  {"read-group-platform", required_argument, 0, 0},	/* sam_read_group_platform */

  {"noexceptions", no_argument, 0, '0'}, /* exception_raise_p */
  {"npaths", required_argument, 0, 'n'}, /* maxpaths */
  {"quiet-if-excessive", no_argument, 0, 'Q'}, /* quiet_if_excessive_p */
  {"ordered", no_argument, 0, 'O'}, /* orderedp */
  {"clip-overlap", no_argument, 0, 0},	     /* clip_overlap_p */
  {"show-refdiff", no_argument, 0, 0},	       /* show_refdiff_p */
  {"print-snps", no_argument, 0, 0}, /* print_snplabels_p */
  {"failsonly", no_argument, 0, 0}, /* failsonlyp */
  {"nofails", no_argument, 0, 0}, /* nofailsp */
  {"fails-as-input", no_argument, 0, 0}, /* fails_as_input_p */
  {"split-output", required_argument, 0, 0}, /* sevenway_root */

#ifdef HAVE_GOBY
  /* Goby-specific options */
  {"goby-output", required_argument, 0, 0}, /* goby_output_root */
  {"creads-window-start", required_argument, 0, 0}, /* creads_window_start */
  {"creads-window-end", required_argument, 0, 0}, /* creads_window_end */
  {"creads-complement", no_argument, 0, 0}, /* creads_complement_p */
#endif

  /* Diagnostic options */
  {"time", no_argument, 0, 0},	/* timingp */
  {"unload", no_argument, 0, 0},	/* unloadp */

  /* Help options */
  {"version", no_argument, 0, 0}, /* print_program_version */
  {"help", no_argument, 0, 0}, /* print_program_usage */
  {0, 0, 0, 0}
};


static void
print_program_version () {
  fprintf(stdout,"\n");
  fprintf(stdout,"GSNAP: Genomic Short Nucleotide Alignment Program\n");
  fprintf(stdout,"Part of GMAP package, version %s\n",PACKAGE_VERSION);
  fprintf(stdout,"Build target: %s\n",TARGET);
  fprintf(stdout,"Features: ");
#ifdef HAVE_PTHREAD
  fprintf(stdout,"pthreads enabled, ");
#else
  fprintf(stdout,"no pthreads, ");
#endif
#ifdef HAVE_ZLIB
  fprintf(stdout,"zlib available, ");
#else
  fprintf(stdout,"no zlib, ");
#endif
#ifdef HAVE_MMAP
  fprintf(stdout,"mmap available, ");
#else
  fprintf(stdout,"no mmap, ");
#endif
#ifdef WORDS_BIGENDIAN
  fprintf(stdout,"bigendian, ");
#else
  fprintf(stdout,"littleendian, ");
#endif
#ifdef HAVE_SIGACTION
  fprintf(stdout,"sigaction available, ");
#else
  fprintf(stdout,"no sigaction, ");
#endif
#ifdef HAVE_64_BIT
  fprintf(stdout,"64 bits available");
#else
  fprintf(stdout,"64 bits not available");
#endif
  fprintf(stdout,"\n");

  fprintf(stdout,"Builtin functions:");
#ifdef HAVE_BUILTIN_CLZ
  fprintf(stdout," clz");
#endif
#ifdef HAVE_BUILTIN_CTZ
  fprintf(stdout," ctz");
#endif
#ifdef HAVE_BUILTIN_POPCOUNT
  fprintf(stdout," popcount");
#endif
  fprintf(stdout,"\n");

  fprintf(stdout,"Sizes: off_t (%lu), size_t (%lu), unsigned int (%lu), long int (%lu)\n",
	  sizeof(off_t),sizeof(size_t),sizeof(unsigned int),sizeof(long int));
  fprintf(stdout,"Default gmap directory: %s\n",GMAPDB);
  fprintf(stdout,"Maximum read length: %d\n",MAX_READLENGTH);
  fprintf(stdout,"Thomas D. Wu, Genentech, Inc.\n");
  fprintf(stdout,"Contact: twu@gene.com\n");
  fprintf(stdout,"\n");
  return;
}

/* This flag is not well-supported, and therefore hidden, but
   kept for backwards compatibility */
/*  -R, --rel=STRING               Release\n\ */


static void
print_program_usage ();


/************************************************************************/


static Result_T
process_request (Request_T request, Floors_T *floors_array,
		 Oligoindex_T *oligoindices_major, int noligoindices_major,
		 Oligoindex_T *oligoindices_minor, int noligoindices_minor,
		 Pairpool_T pairpool, Diagpool_T diagpool,
		 Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
		 Stopwatch_T worker_stopwatch) {
  int jobid;
  Shortread_T queryseq1, queryseq2;
  Stage3end_T *stage3array, *stage3array5, *stage3array3;
  Stage3pair_T *stage3pairarray;

  int npaths, npaths5, npaths3, i;
  int second_absmq, second_absmq5, second_absmq3;
  Pairtype_T final_pairtype;
  double worker_runtime;

  jobid = Request_id(request);
  queryseq1 = Request_queryseq1(request);
  queryseq2 = Request_queryseq2(request);

  /* printf("%s\n",Shortread_accession(queryseq1)); */

  if (worker_stopwatch != NULL) {
    Stopwatch_start(worker_stopwatch);
  }

  if (queryseq2 == NULL) {
    stage3array = Stage1_single_read(&npaths,&second_absmq,queryseq1,indexdb,indexdb2,indexdb_size_threshold,
				     genome,floors_array,maxpaths,user_maxlevel_float,subopt_levels,
				     indel_penalty_middle,indel_penalty_end,
				     max_middle_insertions,max_middle_deletions,
				     allow_end_indels_p,max_end_insertions,max_end_deletions,min_indel_end_matches,
				     shortsplicedist,localsplicing_penalty,distantsplicing_penalty,
				     min_distantsplicing_end_matches,min_distantsplicing_identity,min_shortend,
				     oligoindices_major,noligoindices_major,
				     oligoindices_minor,noligoindices_minor,pairpool,diagpool,
				     dynprogL,dynprogM,dynprogR,
				     /*keep_floors_p*/true);

    worker_runtime = worker_stopwatch == NULL ? 0.00 : Stopwatch_stop(worker_stopwatch);
    return Result_single_read_new(jobid,(void **) stage3array,npaths,second_absmq,worker_runtime);

  } else if ((stage3pairarray = Stage1_paired_read(&npaths,&second_absmq,&final_pairtype,
						   &stage3array5,&npaths5,&second_absmq5,
						   &stage3array3,&npaths3,&second_absmq3,
						   queryseq1,queryseq2,indexdb,indexdb2,indexdb_size_threshold,
						   genome,floors_array,/*maxpaths*/maxpaths,user_maxlevel_float,subopt_levels,
						   indel_penalty_middle,indel_penalty_end,
						   max_middle_insertions,max_middle_deletions,
						   allow_end_indels_p,max_end_insertions,max_end_deletions,min_indel_end_matches,
						   shortsplicedist,localsplicing_penalty,distantsplicing_penalty,
						   min_distantsplicing_end_matches,min_distantsplicing_identity,min_shortend,
						   oligoindices_major,noligoindices_major,
						   oligoindices_minor,noligoindices_minor,pairpool,diagpool,
						   dynprogL,dynprogM,dynprogR,
						   pairmax,/*keep_floors_p*/true)) != NULL) {
    /* Paired or concordant hits found */
    worker_runtime = worker_stopwatch == NULL ? 0.00 : Stopwatch_stop(worker_stopwatch);
    return Result_paired_read_new(jobid,(void **) stage3pairarray,npaths,second_absmq,final_pairtype,worker_runtime);

  } else if (chop_primers_p == false || Shortread_chop_primers(queryseq1,queryseq2) == false) {
    /* No paired or concordant hits found, and no adapters found */
    /* Report ends as unpaired */
    worker_runtime = worker_stopwatch == NULL ? 0.00 : Stopwatch_stop(worker_stopwatch);
    return Result_paired_as_singles_new(jobid,(void **) stage3array5,npaths5,second_absmq5,
					(void **) stage3array3,npaths3,second_absmq3,worker_runtime);

  } else {
    /* Try with potential primers chopped.  queryseq1 and queryseq2 altered by Shortread_chop_primers. */
    for (i = 0; i < npaths5; i++) {
      Stage3end_free(&(stage3array5[i]));
    }
    FREE_OUT(stage3array5);

    for (i = 0; i < npaths3; i++) {
      Stage3end_free(&(stage3array3[i]));
    }
    FREE_OUT(stage3array3);

    if ((stage3pairarray = Stage1_paired_read(&npaths,&second_absmq,&final_pairtype,
					      &stage3array5,&npaths5,&second_absmq5,
					      &stage3array3,&npaths3,&second_absmq3,
					      queryseq1,queryseq2,indexdb,indexdb2,indexdb_size_threshold,
					      genome,floors_array,/*maxpaths*/maxpaths,user_maxlevel_float,subopt_levels,
					      indel_penalty_middle,indel_penalty_end,
					      max_middle_insertions,max_middle_deletions,
					      allow_end_indels_p,max_end_insertions,max_end_deletions,min_indel_end_matches,
					      shortsplicedist,localsplicing_penalty,distantsplicing_penalty,
					      min_distantsplicing_end_matches,min_distantsplicing_identity,min_shortend,
					      oligoindices_major,noligoindices_major,
					      oligoindices_minor,noligoindices_minor,pairpool,diagpool,
					      dynprogL,dynprogM,dynprogR,
					      pairmax,/*keep_floors_p*/false)) != NULL) {
      /* Paired or concordant hits found, after chopping adapters */
      worker_runtime = worker_stopwatch == NULL ? 0.00 : Stopwatch_stop(worker_stopwatch);
      return Result_paired_read_new(jobid,(void **) stage3pairarray,npaths,second_absmq,final_pairtype,worker_runtime);

    } else {
      /* No paired or concordant hits found, after chopping adapters */
      worker_runtime = worker_stopwatch == NULL ? 0.00 : Stopwatch_stop(worker_stopwatch);
      return Result_paired_as_singles_new(jobid,(void **) stage3array5,npaths5,second_absmq5,
					  (void **) stage3array3,npaths3,second_absmq3,worker_runtime);
    }
  }
}



#ifdef HAVE_SIGACTION
static const Except_T sigfpe_error = {"SIGFPE--arithmetic exception"};
static const Except_T sigsegv_error = {"SIGSEGV--segmentation violation"};
static const Except_T sigtrap_error = {"SIGTRAP--hardware fault"};
static const Except_T misc_signal_error = {"Miscellaneous signal"};

#if 0
static void
signal_handler_old (int sig) {
  if (sig == SIGUSR1) {
#ifdef HAVE_PTHREAD
    pthread_exit(NULL);
#else
    exit(9);
#endif
  } else if (sig == SIGFPE) {
    Except_raise(&sigfpe_error,__FILE__,__LINE__);
  } else if (sig == SIGSEGV) {
    Except_raise(&sigsegv_error,__FILE__,__LINE__);
  } else if (sig == SIGTRAP) {
    Except_raise(&sigtrap_error,__FILE__,__LINE__);
  } else {
    fprintf(stderr,"Signal %d\n",sig);
    Except_raise(&misc_signal_error,__FILE__,__LINE__);
  }
  return;
}
#endif

static void
signal_handler (int sig) {
  Request_T request;
  Shortread_T queryseq1, queryseq2;

  if (sig == SIGFPE) {
    fprintf(stderr,"Signal received: Floating point error\n");
  } else if (sig == SIGSEGV) {
    fprintf(stderr,"Signal received: Segmentation fault\n");
  } else if (sig == SIGTRAP) {
    fprintf(stderr,"Signal received: Trap\n");
  } else {
    fprintf(stderr,"Signal received: %d\n",sig);
  }


#ifdef HAVE_PTHREAD
  request = (Request_T) pthread_getspecific(global_request_key);
  if (request == NULL) {
    fprintf(stderr,"Unable to retrieve request for thread\n");
  } else {
    queryseq1 = Request_queryseq1(request);
    queryseq2 = Request_queryseq2(request);
    if (queryseq1 == NULL) {
      fprintf(stderr,"Unable to retrieve queryseq for request\n");
    } else {
	fprintf(stderr,"Problem sequence: ");
	fprintf(stderr,"%s (%d bp)\n",Shortread_accession(queryseq1),Shortread_fulllength(queryseq1));
	if (queryseq2 == NULL) {
	  Shortread_print_query_singleend_fasta(stderr,queryseq1);
	} else {
	  Shortread_print_query_pairedend_fasta(stderr,queryseq1,queryseq2,
						invert_first_p,invert_second_p);
      }
    }
  }
#endif

  abort();
  return;
}

#endif


#define POOL_FREE_INTERVAL 200

static void
single_thread () {
  Floors_T *floors_array;
  Request_T request;
  Result_T result;
  Shortread_T queryseq1;
  int i;
  int noutput = 0;
  Stopwatch_T worker_stopwatch;

  /* For GMAP */
  Oligoindex_T *oligoindices_major, *oligoindices_minor;
  int noligoindices_major, noligoindices_minor;
  Dynprog_T dynprogL, dynprogM, dynprogR;
  Pairpool_T pairpool;
  Diagpool_T diagpool;
  int jobid = 0;

#ifdef MEMUSAGE
  long int memusage_constant = 0;
#endif

  oligoindices_major = Oligoindex_new_major(&noligoindices_major);
  oligoindices_minor = Oligoindex_new_minor(&noligoindices_minor);
  dynprogL = Dynprog_new(nullgap,EXTRAQUERYGAP,maxpeelback,extramaterial_end,extramaterial_paired);
  dynprogM = Dynprog_new(nullgap,EXTRAQUERYGAP,maxpeelback,extramaterial_end,extramaterial_paired);
  dynprogR = Dynprog_new(nullgap,EXTRAQUERYGAP,maxpeelback,extramaterial_end,extramaterial_paired);
  pairpool = Pairpool_new();
  diagpool = Diagpool_new();
  worker_stopwatch = (timingp == true) ? Stopwatch_new() : (Stopwatch_T) NULL;

  floors_array = (Floors_T *) CALLOC(MAX_READLENGTH+1,sizeof(Floors_T));
  /* Except_stack_create(); -- requires pthreads */

#ifdef MEMUSAGE
  memusage_constant += Mem_usage_report();
  Mem_usage_reset(0);
#endif

  while ((request = Inbuffer_get_request(inbuffer)) != NULL) {
    debug(printf("single_thread got request %d\n",Request_id(request)));

    TRY
      result = process_request(request,floors_array,oligoindices_major,noligoindices_major,
			       oligoindices_minor,noligoindices_minor,
			       pairpool,diagpool,dynprogL,dynprogM,dynprogR,worker_stopwatch);
    ELSE
      queryseq1 = Request_queryseq1(request);
      if (queryseq1 == NULL) {
	fprintf(stderr,"NULL");
      } else if (Shortread_accession(queryseq1) == NULL) {
	fprintf(stderr,"unnamed (%d bp)",Shortread_fulllength(queryseq1));
      } else {
	fprintf(stderr,"Problem sequence: ");
	fprintf(stderr,"%s (%d bp)",Shortread_accession(queryseq1),Shortread_fulllength(queryseq1));
      }
      fprintf(stderr,"\n");
      if (Request_queryseq2(request) == NULL) {
	Shortread_print_query_singleend_fasta(stderr,queryseq1);
      } else {
	Shortread_print_query_pairedend_fasta(stderr,queryseq1,Request_queryseq2(request),
					      invert_first_p,invert_second_p);
      }
      fprintf(stderr,"\n");
      fprintf(stderr,"To obtain a core dump, re-run program on problem sequence with the -0 [zero] flag\n");

      fprintf(stderr,"Exiting...\n");
      exit(9);
    RERAISE;
    END_TRY;

#ifdef MEMUSAGE
    Outbuffer_print_result(outbuffer,result,request,noutput+1);
    printf("***%s\n",Shortread_accession(Request_queryseq1(request)));
#else
    Outbuffer_print_result(outbuffer,result,request);
#endif

    Result_free(&result);

    /* Allocated by fill_buffer in Inbuffer_get_request */
    Request_free(&request);
    noutput++;

    if (jobid % POOL_FREE_INTERVAL == 0) {
      Pairpool_free_memory(pairpool);
      Diagpool_free_memory(diagpool);
    }

#ifdef MEMUSAGE
    printf("Memusage of single thread: %ld\n",Mem_usage_report());
    printf("Memusage of OUT: %ld\n",Mem_usage_out_report());
    assert(Mem_usage_report() == 0);
    assert(Mem_usage_out_report() == 0);
#endif

  }

  /* Except_stack_destroy(); -- requires pthreads */

#ifdef MEMUSAGE
  Mem_usage_add(memusage_constant);
#endif

  for (i = 0; i <= MAX_READLENGTH; i++) {
    if (floors_array[i] != NULL) {
      Floors_free_keep(&(floors_array[i]));
    }
  }
  FREE(floors_array);

  if (worker_stopwatch != NULL) {
    Stopwatch_free(&worker_stopwatch);
  }
  Diagpool_free(&diagpool);
  Pairpool_free(&pairpool);
  Dynprog_free(&dynprogR);
  Dynprog_free(&dynprogM);
  Dynprog_free(&dynprogL);
  Oligoindex_free_array(&oligoindices_minor,noligoindices_minor);
  Oligoindex_free_array(&oligoindices_major,noligoindices_major);

  return;
}



static void *
worker_thread (void *data) {
  long int worker_id = (long int) data;
  Floors_T *floors_array;
  Request_T request;
  Result_T result;
  Shortread_T queryseq1;
  int i;
  Stopwatch_T worker_stopwatch;

  /* For GMAP */
  Oligoindex_T *oligoindices_major, *oligoindices_minor;
  int noligoindices_major, noligoindices_minor;
  Dynprog_T dynprogL, dynprogM, dynprogR;
  Pairpool_T pairpool;
  Diagpool_T diagpool;
  int worker_jobid = 0;

#ifdef MEMUSAGE
  long int memusage_constant = 0, memusage;
  char threadname[12];
  sprintf(threadname,"thread-%ld",worker_id);
  Mem_usage_set_threadname(threadname);
#endif

  /* Thread-specific data and storage */
  oligoindices_major = Oligoindex_new_major(&noligoindices_major);
  oligoindices_minor = Oligoindex_new_minor(&noligoindices_minor);
  dynprogL = Dynprog_new(nullgap,EXTRAQUERYGAP,maxpeelback,extramaterial_end,extramaterial_paired);
  dynprogM = Dynprog_new(nullgap,EXTRAQUERYGAP,maxpeelback,extramaterial_end,extramaterial_paired);
  dynprogR = Dynprog_new(nullgap,EXTRAQUERYGAP,maxpeelback,extramaterial_end,extramaterial_paired);
  pairpool = Pairpool_new();
  diagpool = Diagpool_new();
  worker_stopwatch = (timingp == true) ? Stopwatch_new() : (Stopwatch_T) NULL;

  floors_array = (Floors_T *) CALLOC(MAX_READLENGTH+1,sizeof(Floors_T));
  Except_stack_create();

#ifdef MEMUSAGE
  memusage_constant += Mem_usage_report();
  Mem_usage_reset(0);
#endif

  while ((request = Inbuffer_get_request(inbuffer)) != NULL) {
    debug(printf("worker_thread %ld got request %d\n",worker_id,Request_id(request)));
    pthread_setspecific(global_request_key,(void *) request);
    if (worker_jobid % POOL_FREE_INTERVAL == 0) {
      Pairpool_free_memory(pairpool);
      Diagpool_free_memory(diagpool);
    }

#ifdef MEMUSAGE
    memusage = Mem_usage_report();
    /* printf("Memusage of worker thread %ld: %ld\n",worker_id,memusage); */
    if (memusage != 0) {
      fprintf(stderr,"Memusage of worker thread %ld: %ld\n",worker_id,memusage);
      fflush(stdout);
      exit(9);
    }
#endif

    TRY
      result = process_request(request,floors_array,oligoindices_major,noligoindices_major,
			       oligoindices_minor,noligoindices_minor,
			       pairpool,diagpool,dynprogL,dynprogM,dynprogR,worker_stopwatch);
    ELSE
      queryseq1 = Request_queryseq1(request);
      if (queryseq1 == NULL) {
	fprintf(stderr,"NULL");
      } else if (Shortread_accession(queryseq1) == NULL) {
	fprintf(stderr,"unnamed (%d bp)",Shortread_fulllength(queryseq1));
      } else {
	fprintf(stderr,"Problem sequence: ");
	fprintf(stderr,"%s (%d bp)",Shortread_accession(queryseq1),Shortread_fulllength(queryseq1));
      }
      fprintf(stderr,"\n");
      if (Request_queryseq2(request) == NULL) {
	Shortread_print_query_singleend_fasta(stderr,queryseq1);
      } else {
	Shortread_print_query_pairedend_fasta(stderr,queryseq1,Request_queryseq2(request),
					      invert_first_p,invert_second_p);
      }
      fprintf(stderr,"\n");
      fprintf(stderr,"To obtain a core dump, re-run program on problem sequence with the -0 [zero] flag\n");

      fprintf(stderr,"Exiting...\n");
      exit(9);
    RERAISE;
    END_TRY;

    debug(printf("worker_thread putting result %d\n",Result_id(result)));

    Outbuffer_put_result(outbuffer,result,request);

    /* Don't free result or request; done by outbuffer thread */
  }

#ifdef MEMUSAGE
  Mem_usage_add(memusage_constant);
#endif

  Except_stack_destroy();

  for (i = 0; i <= MAX_READLENGTH; i++) {
    if (floors_array[i] != NULL) {
      Floors_free_keep(&(floors_array[i]));
    }
  }
  FREE(floors_array);

  if (worker_stopwatch != NULL) {
    Stopwatch_free(&worker_stopwatch);
  }
  Diagpool_free(&diagpool);
  Pairpool_free(&pairpool);
  Dynprog_free(&dynprogR);
  Dynprog_free(&dynprogM);
  Dynprog_free(&dynprogL);
  Oligoindex_free_array(&oligoindices_minor,noligoindices_minor);
  Oligoindex_free_array(&oligoindices_major,noligoindices_major);

  return (void *) NULL;
}


static void
parse_part (int *part_modulus, int *part_interval, char *string) {
  char *p = string;

  if (sscanf(p,"%d",&(*part_modulus)) < 1) {
    fprintf(stderr,"Cannot parse first integer from %s\n",string);
    exit(9);
  }

  while (*p != '\0' && isdigit(*p)) {
    p++;
  }
  while (*p != '\0' && !isdigit(*p)) {
    p++;
  }
  if (sscanf(p,"%d",&(*part_interval)) < 1) {
    fprintf(stderr,"Cannot parse first integer from %s\n",string);
    exit(9);
  }
  if ((*part_modulus) >= (*part_interval)) {
    fprintf(stderr,"In %s, batch number %d must be less than the number of batches %d\n",
	    string,*part_modulus,*part_interval);
    exit(9);
  }
  if (*part_interval == 0) {
    fprintf(stderr,"Bad batch specification %s.  Batch interval cannot be 0.\n",string);
    exit(9);
  }

  return;
}


static int
add_gmap_mode (char *string) {
  if (!strcmp(string,"none")) {
    gmap_mode = 0;
    return 0;
  } else {
    if (!strcmp(string,"improve")) {
      gmap_mode |= GMAP_IMPROVEMENT;
    } else if (!strcmp(string,"terminal")) {
      gmap_mode |= GMAP_TERMINAL;
    } else if (!strcmp(string,"pairsearch")) {
      gmap_mode |= GMAP_PAIRSEARCH;
    } else {
      fprintf(stderr,"Don't recognize gmap-mode type %s\n",string);
      fprintf(stderr,"Allowed values are: none, improve, terminal, pairsearch\n");
      exit(9);
    }
    return 1;
  }
}


char *
check_valid_int (char *string) {
  char *p = string;

  if (*p == '+' || *p == '-') {
    p++;
  }

  if (!isdigit(*p)) {
    fprintf(stderr,"value %s is not a valid int\n",string);
    exit(9);
    return NULL;
  }
  while (*p != '\0' && isdigit(*p)) {
    p++;
  }

  if (*p == 'e') {
    p++;
    if (*p == '+') {
      p++;
    }
    if (!isdigit(*p)) {
      return false;
    }
    while (*p != '\0' && isdigit(*p)) {
      p++;
    }
  }

  if (*p == '\0') {
    return string;
  } else {
    fprintf(stderr,"value %s is not a valid int\n",string);
    exit(9);
    return NULL;
  }
}


char *
check_valid_float (char *string) {
  char *p = string;

  if (*p == '+' || *p == '-') {
    p++;
  }

  while (*p != '\0' && isdigit(*p)) {
    p++;
  }
  if (*p == '\0') {
    return string;
  }

  if (*p == '.') {
    p++;
  }

  if (!isdigit(*p)) {
    fprintf(stderr,"value %s is not a valid float\n",string);
    exit(9);
    return NULL;
  }
  while (*p != '\0' && isdigit(*p)) {
    p++;
  }

  if (*p == 'e') {
    p++;
    if (*p == '+' || *p == '-') {
      p++;
    }
    if (!isdigit(*p)) {
      fprintf(stderr,"value %s is not a valid float\n",string);
      exit(9);
      return NULL;
    }
    while (*p != '\0' && isdigit(*p)) {
      p++;
    }
  }

  if (*p == '\0') {
    return string;
  } else {
    fprintf(stderr,"value %s is not a valid float\n",string);
    exit(9);
    return NULL;
  }
}



int
main (int argc, char *argv[]) {
  char *genomesubdir = NULL, *snpsdir = NULL, *modedir = NULL, *mapdir = NULL, *iitfile = NULL, *fileroot = NULL;
  FILE *input = NULL, *input2 = NULL;
#ifdef HAVE_ZLIB
  gzFile gzipped = NULL, gzipped2 = NULL;
#endif

  bool multiple_sequences_p = false;
  char **files;
  int nfiles, nextchar = '\0';
  long int worker_id;

  unsigned int nread;
  double runtime;

#ifdef HAVE_PTHREAD
  int ret;
  pthread_attr_t thread_attr_join;
#ifdef WORKER_DETACH
  pthread_attr_t thread_attr_detach;
#endif
#endif

  int opt, c;
  extern int optind;
  extern char *optarg;
  int long_option_index = 0;
  const char *long_name;
  char **argstart;
  char *string;

#ifdef HAVE_SIGACTION
  struct sigaction signal_action;
#endif

#ifdef MEMUSAGE
  Mem_usage_init();
  Mem_usage_set_threadname("main");
#endif


  fprintf(stderr,"GSNAP version %s called with args:",PACKAGE_VERSION);
  argstart = &(argv[-optind]);
  for (c = 1; c < argc + optind; c++) {
    fprintf(stderr," %s",argstart[c]);
  }
  fprintf(stderr,"\n");


  while ((opt = getopt_long(argc,argv,
			    "D:d:k:Gq:o:a:N:M:m:i:y:Y:z:Z:w:E:e:J:K:l:g:s:V:v:B:t:A:j:0n:QO",
			    long_options, &long_option_index)) != -1) {
    switch (opt) {
    case 0:
      long_name = long_options[long_option_index].name;
      if (!strcmp(long_name,"version")) {
	print_program_version();
	exit(0);
      } else if (!strcmp(long_name,"help")) {
	print_program_usage();
	exit(0);

      } else if (!strcmp(long_name,"time")) {
	timingp = true;

      } else if (!strcmp(long_name,"unload")) {
	unloadp = true;

      } else if (!strcmp(long_name,"mode")) {
	if (!strcmp(optarg,"standard")) {
	  mode = STANDARD;
	} else if (!strcmp(optarg,"cmet-stranded")) {
	  mode = CMET_STRANDED;
	} else if (!strcmp(optarg,"cmet-nonstranded")) {
	  mode = CMET_NONSTRANDED;
	} else if (!strcmp(optarg,"atoi-stranded")) {
	  mode = ATOI_STRANDED;
	} else if (!strcmp(optarg,"atoi-nonstranded")) {
	  mode = ATOI_NONSTRANDED;
	} else {
	  fprintf(stderr,"--mode must be standard, cmet-stranded, cmet-nonstranded, atoi-stranded, or atoi\n");
	  exit(9);
	}

      } else if (!strcmp(long_name,"cmetdir")) {
	user_cmetdir = optarg;
      } else if (!strcmp(long_name,"atoidir")) {
	user_atoidir = optarg;

      } else if (!strcmp(long_name,"splicingdir")) {
	user_splicingdir = optarg;
      } else if (!strcmp(long_name,"ambig-splice-noclip")) {
	amb_clip_p = false;

      } else if (!strcmp(long_name,"tallydir")) {
	user_tallydir = optarg;
      } else if (!strcmp(long_name,"use-tally")) {
	tally_root = optarg;

      } else if (!strcmp(long_name,"runlengthdir")) {
	user_runlengthdir = optarg;
      } else if (!strcmp(long_name,"use-runlength")) {
	runlength_root = optarg;

      } else if (!strcmp(long_name,"gmap-mode")) {
	gmap_mode = 0;
	string = strtok(optarg,",");
	if (add_gmap_mode(string) != 0) {
	  while ((string = strtok(NULL,",")) != NULL && add_gmap_mode(string) != 0) {
	  }
	}

      } else if (!strcmp(long_name,"trigger-score-for-gmap")) {
	trigger_score_for_gmap = atoi(check_valid_int(optarg));

      } else if (!strcmp(long_name,"max-gmap-pairsearch")) {
	max_gmap_pairsearch = atoi(check_valid_int(optarg));
      } else if (!strcmp(long_name,"max-gmap-terminal")) {
	max_gmap_terminal = atoi(check_valid_int(optarg));
      } else if (!strcmp(long_name,"max-gmap-improvement")) {
	max_gmap_improvement = atoi(check_valid_int(optarg));
      } else if (!strcmp(long_name,"microexon-spliceprob")) {
	microexon_spliceprob = atof(check_valid_float(optarg));
      } else if (!strcmp(long_name,"stage2-start")) {
	suboptimal_score_start = atoi(check_valid_int(optarg));
      } else if (!strcmp(long_name,"stage2-end")) {
	suboptimal_score_end = atoi(check_valid_int(optarg));

      } else if (!strcmp(long_name,"input-buffer-size")) {
	inbuffer_nspaces = atoi(check_valid_int(optarg));
      } else if (!strcmp(long_name,"output-buffer-size")) {
	output_buffer_size = atoi(check_valid_int(optarg));
      } else if (!strcmp(long_name,"barcode-length")) {
	barcode_length = atoi(check_valid_int(optarg));
      } else if (!strcmp(long_name,"fastq-id-start")) {
	acc_fieldi_start = atoi(check_valid_int(optarg)) - 1;
	if (acc_fieldi_start < 0) {
	  fprintf(stderr,"Value for fastq-id-start must be 1 or greater\n");
	  exit(9);
	}
      } else if (!strcmp(long_name,"fastq-id-end")) {
	acc_fieldi_end = atoi(check_valid_int(optarg)) - 1;
	if (acc_fieldi_end < 0) {
	  fprintf(stderr,"Value for fastq-id-end must be 1 or greater\n");
	  exit(9);
	}
      } else if (!strcmp(long_name,"filter-chastity")) {
	if (!strcmp(optarg,"off")) {
	  filter_chastity_p = false;
	  filter_if_both_p = false;
	} else if (!strcmp(optarg,"either")) {
	  filter_chastity_p = true;
	  filter_if_both_p = false;
	} else if (!strcmp(optarg,"both")) {
	  filter_chastity_p = true;
	  filter_if_both_p = true;
	} else {
	  fprintf(stderr,"--filter-chastity values allowed: off, either, both\n");
	  exit(9);
	}

#ifdef HAVE_ZLIB
      } else if (!strcmp(long_name,"gunzip")) {
	gunzip_p = true;
#endif
      } else if (!strcmp(long_name,"split-output")) {
	sevenway_root = optarg;
      } else if (!strcmp(long_name,"fails-as-input")) {
	fails_as_input_p = true;
      } else if (!strcmp(long_name,"pairmax-dna")) {
	pairmax_dna = atoi(check_valid_int(optarg));
      } else if (!strcmp(long_name,"pairmax-rna")) {
	pairmax_rna = atoi(check_valid_int(optarg));
      } else if (!strcmp(long_name,"pairexpect")) {
	expected_pairlength = atoi(check_valid_int(optarg));
      } else if (!strcmp(long_name,"pairdev")) {
	pairlength_deviation = atoi(check_valid_int(optarg));
      } else if (!strcmp(long_name,"indel-endlength")) {
	min_indel_end_matches = atoi(check_valid_int(optarg));
	if (min_indel_end_matches > 14) {
	  allow_end_indels_p = false;
	}

      } else if (!strcmp(long_name,"terminal-threshold")) {
	terminal_threshold = atoi(check_valid_int(optarg));

      } else if (!strcmp(long_name,"antistranded-penalty")) {
	antistranded_penalty = atoi(check_valid_int(optarg));

      } else if (!strcmp(long_name,"favor-multiexon")) {
	favor_multiexon_p = true;

      } else if (!strcmp(long_name,"merge-distant-samechr")) {
	merge_samechr_p = true;

      } else if (!strcmp(long_name,"query-unk-mismatch")) {
	if (!strcmp(optarg,"1")) {
	  query_unk_mismatch_p = true;
	} else if (!strcmp(optarg,"0")) {
	  query_unk_mismatch_p = false;
	} else {
	  fprintf(stderr,"--query-unk-mismatch flag must be 0 or 1\n");
	  exit(9);
	}
      } else if (!strcmp(long_name,"genome-unk-mismatch")) {
	if (!strcmp(optarg,"1")) {
	  genome_unk_mismatch_p = true;
	} else if (!strcmp(optarg,"0")) {
	  genome_unk_mismatch_p = false;
	} else {
	  fprintf(stderr,"--genome-unk-mismatch flag must be 0 or 1\n");
	  exit(9);
	}

      } else if (!strcmp(long_name,"trim-mismatch-score")) {
	trim_mismatch_score = atoi(check_valid_int(optarg));
      } else if (!strcmp(long_name,"trim-indel-score")) {
	trim_indel_score = atoi(check_valid_int(optarg));

      } else if (!strcmp(long_name,"show-refdiff")) {
	show_refdiff_p = true;
      } else if (!strcmp(long_name,"clip-overlap")) {
	clip_overlap_p = true;
      } else if (!strcmp(long_name,"no-sam-headers")) {
	sam_headers_p = false;
      } else if (!strcmp(long_name,"sam-headers-batch")) {
	sam_headers_batch = atoi(check_valid_int(optarg));
      } else if (!strcmp(long_name,"quality-protocol")) {
	if (user_quality_score_adj == true) {
	  fprintf(stderr,"Cannot specify both -J (--quality-zero-score) and --quality-protocol\n");
	  exit(9);
	} else if (user_quality_shift == true) {
	  fprintf(stderr,"Cannot specify both -j (--quality-print-shift) and --quality-protocol\n");
	  exit(9);
	} else if (!strcmp(optarg,"illumina")) {
	  MAPQ_init(/*quality_score_adj*/64);
	  Pair_init(/*quality_score_adj*/64);
	  user_quality_score_adj = true;
	  quality_shift = -31;
	  user_quality_shift = true;
	} else if (!strcmp(optarg,"sanger")) {
	  MAPQ_init(/*quality_score_adj*/33);
	  Pair_init(/*quality_score_adj*/33);
	  user_quality_score_adj = true;
	  quality_shift = 0;
	  user_quality_shift = true;
	} else {
	  fprintf(stderr,"The only values allowed for --quality-protocol are illumina or sanger\n");
	  exit(9);
	}

      } else if (!strcmp(long_name,"read-group-id")) {
	sam_read_group_id = optarg;
      } else if (!strcmp(long_name,"read-group-name")) {
	sam_read_group_name = optarg;
      } else if (!strcmp(long_name,"read-group-library")) {
	sam_read_group_library = optarg;
      } else if (!strcmp(long_name,"read-group-platform")) {
	sam_read_group_platform = optarg;
      } else if (!strcmp(long_name,"goby-output")) {
	goby_output_root = optarg;
      } else if (!strcmp(long_name,"distant-splice-identity")) {
	min_distantsplicing_identity = atof(check_valid_float(optarg));
      } else if (!strcmp(long_name,"print-snps")) {
	print_snplabels_p = true;
      } else if (!strcmp(long_name,"failsonly")) {
	if (nofailsp == true) {
	  fprintf(stderr,"Cannot specify both --nofails and --failsonly\n");
	  exit(9);
	} else {
	  failsonlyp = true;
	}
      } else if (!strcmp(long_name,"nofails")) {
	if (failsonlyp == true) {
	  fprintf(stderr,"Cannot specify both --nofails and --failsonly\n");
	  exit(9);
	} else {
	  nofailsp = true;
	}
      } else if (!strcmp(long_name,"creads-window-start")) {
	creads_window_start = strtoul(check_valid_int(optarg),NULL,10);
      } else if (!strcmp(long_name,"creads-window-end")) {
	creads_window_end = strtoul(check_valid_int(optarg),NULL,10);
      } else if (!strcmp(long_name,"creads-complement")) {
	creads_complement_p = true;
      } else {
	/* Shouldn't reach here */
	fprintf(stderr,"Don't recognize option %s.  For usage, run 'gsnap --help'",long_name);
	exit(9);
      }
      break;

    case 'D': user_genomedir = optarg; break;
    case 'd': 
      dbroot = (char *) CALLOC(strlen(optarg)+1,sizeof(char));
      strcpy(dbroot,optarg);
      break;
    case 'k':
      required_index1part = atoi(check_valid_int(optarg));
      if (required_index1part >= 12 && required_index1part <= 15) {
	/* Okay */
      } else {
	fprintf(stderr,"The only values allowed for -k are 12, 13, 14, or 15\n");
	exit(9);
      }
      break;
    case 'G': uncompressedp = true; break;

    case 'q': parse_part(&part_modulus,&part_interval,optarg); break;
    case 'o': 
      if (!strcmp(optarg,"FR")) {
	invert_first_p = false;
	invert_second_p = true;
      } else if (!strcmp(optarg,"RF")) {
	invert_first_p = true;
	invert_second_p = false;
      } else if (!strcmp(optarg,"FF")) {
	invert_first_p = invert_second_p = false;
      } else {
	fprintf(stderr,"Currently allowed values for orientation (-o): FR (fwd-rev), RF (rev-fwd) or FF (fwd-fwd)\n");
	exit(9);
      }
      break;

    case 'a': 
      if (!strcmp(optarg,"paired")) {
	chop_primers_p = true;
      } else if (!strcmp(optarg,"off")) {
	chop_primers_p = true;
      } else {
	fprintf(stderr,"Currently allowed values for adapter stripping (-a): off, paired\n");
	exit(9);
      }
      break;

    case 'N':
      if (!strcmp(optarg,"1")) {
	novelsplicingp = true;
      } else if (!strcmp(optarg,"0")) {
	novelsplicingp = false;
      } else {
	fprintf(stderr,"Novel splicing (-N flag) must be 0 or 1\n");
	exit(9);
      }
      break;

#if 0
    case 'R': 
      if (!strcmp(optarg,"0")) {
	masktype = MASK_NONE;
      } else if (!strcmp(optarg,"1")) {
	masktype = MASK_FREQUENT;
      } else if (!strcmp(optarg,"2")) {
	masktype = MASK_REPETITIVE;
      } else if (!strcmp(optarg,"3")) {
	masktype = MASK_GREEDY_FREQUENT;
      } else if (!strcmp(optarg,"4")) {
	masktype = MASK_GREEDY_REPETITIVE;
      } else {
	fprintf(stderr,"Masking mode %s not recognized.\n",optarg);
	fprintf(stderr,"Mode 0 means no masking, mode 1 masks frequent oligomers;\n");
	fprintf(stderr,"  mode 2 masks frequent and repetitive oligomers;\n");
	fprintf(stderr,"  mode 3 does greedy masking of frequent oligomers,\n");
	fprintf(stderr,"    then no masking if necessary;\n");
	fprintf(stderr,"  mode 4 does greedy masking of frequent and repetitive oligomers,\n");
	fprintf(stderr,"    then no masking if necessary.\n");
	exit(9);
      }
      break;
#endif

    case 'M': subopt_levels = atoi(check_valid_int(optarg)); break;
    case 'm':
      user_maxlevel_float = atof(check_valid_float(optarg));
      if (user_maxlevel_float > 1.0 && user_maxlevel_float != rint(user_maxlevel_float)) {
	fprintf(stderr,"Cannot specify fractional value %f for --max-mismatches except between 0.0 and 1.0\n",user_maxlevel_float);
	exit(9);
      }
      break;

    case 'i': indel_penalty_middle = indel_penalty_end = atoi(check_valid_int(optarg)); break;

    case 'y': max_middle_insertions = atoi(check_valid_int(optarg)); break;
    case 'Y': max_end_insertions = atoi(check_valid_int(optarg)); break;
    case 'z': max_middle_deletions = atoi(check_valid_int(optarg)); break;
    case 'Z': max_end_deletions = atoi(check_valid_int(optarg)); break;

    case 'w': shortsplicedist = strtoul(optarg,NULL,10); break;

    case 'E': distantsplicing_penalty = atoi(check_valid_int(optarg)); break;
    case 'e': localsplicing_penalty = atoi(check_valid_int(optarg)); break;
    case 'K': min_distantsplicing_end_matches = atoi(check_valid_int(optarg)); break;
    case 'l': min_shortend = atoi(check_valid_int(optarg)); break;

    case 'g': genes_file = optarg; break;

    case 's': splicing_file = optarg; knownsplicingp = true; break;

    case 'V': user_snpsdir = optarg; break;
    case 'v': snps_root = optarg; break;

    case 'B':
      if (!strcmp(optarg,"0")) {
	offsetscomp_access = USE_ALLOCATE; /* was batch_offsets_p = true */
	positions_access = USE_MMAP_ONLY; /* was batch_positions_p = false */
	genome_access = USE_MMAP_ONLY; /* was batch_genome_p = false */
      } else if (!strcmp(optarg,"1")) {
	offsetscomp_access = USE_ALLOCATE; /* was batch_offsets_p = true */
	positions_access = USE_MMAP_PRELOAD; /* was batch_positions_p = true */
	genome_access = USE_MMAP_ONLY; /* was batch_genome_p = false */
      } else if (!strcmp(optarg,"2")) {
	offsetscomp_access = USE_ALLOCATE; /* was batch_offsets_p = true */
	positions_access = USE_MMAP_PRELOAD; /* was batch_positions_p = true */
	genome_access = USE_MMAP_PRELOAD; /* was batch_genome_p = true */
      } else if (!strcmp(optarg,"3")) {
	offsetscomp_access = USE_ALLOCATE;
	positions_access = USE_ALLOCATE;
	genome_access = USE_MMAP_PRELOAD; /* was batch_genome_p = true */
      } else if (!strcmp(optarg,"4")) {
	offsetscomp_access = USE_ALLOCATE;
	positions_access = USE_ALLOCATE;
	genome_access = USE_ALLOCATE;
      } else if (!strcmp(optarg,"5")) {
	expand_offsets_p = true;
	offsetscomp_access = USE_ALLOCATE; /* Doesn't matter */
	positions_access = USE_ALLOCATE;
	genome_access = USE_ALLOCATE;
      } else {
	fprintf(stderr,"Batch mode %s not recognized.  Only allow 0-5.  Run 'gsnap --help' for more information.\n",optarg);
	exit(9);
      }
      break;

#ifdef HAVE_PTHREAD
    case 't': nworkers = atoi(check_valid_int(optarg)); break;
#endif

    case 'A':
      if (!strcmp(optarg,"sam")) {
	output_sam_p = true;
      } else if (!strcmp(optarg,"goby")) {
	output_goby_p = true;
      } else {
	fprintf(stderr,"Output format %s not recognized\n",optarg);
	exit(9);
      }
      break;

    case 'j':
      if (user_quality_shift == true) {
	fprintf(stderr,"Cannot specify both -j (--quality-print-shift) and --quality-protocol\n");
	exit(9);
      } else {
	quality_shift = atoi(check_valid_int(optarg));
	user_quality_shift = true;
      }
      break;

    case 'J':
      if (user_quality_score_adj == true) {
	fprintf(stderr,"Cannot specify both -J (--quality-zero-score) and --quality-protocol\n");
	exit(9);
      } else {
	MAPQ_init(/*quality_score_adj*/atoi(check_valid_int(optarg)));
	Pair_init(/*quality_score_adj*/atoi(check_valid_int(optarg)));
	user_quality_score_adj = true;
      }
      break;

    case '0': exception_raise_p = false; break; /* Allows signals to pass through */
    case 'n': maxpaths = atoi(check_valid_int(optarg)); break;
    case 'Q': quiet_if_excessive_p = true; break;

    case 'O': orderedp = true; break;

    case '?': fprintf(stderr,"For usage, run 'gsnap --help'\n"); exit(9);
    default: exit(9);
    }
  }
  argc -= optind;
  argv += optind;


  if (exception_raise_p == false) {
    fprintf(stderr,"Allowing signals and exceptions to pass through\n");
    Except_inactivate();
  } else {
#ifdef HAVE_SIGACTION
    signal_action.sa_handler = signal_handler;
    signal_action.sa_flags = 0;
    sigfillset(&signal_action.sa_mask);

    sigaction(SIGFPE,&signal_action,NULL);
    sigaction(SIGSEGV,&signal_action,NULL);
    sigaction(SIGTRAP,&signal_action,NULL);
    sigaction(SIGUSR1,&signal_action,NULL);
#endif
  }


  if (dbroot == NULL) {
    fprintf(stderr,"Need to specify the -d flag.  For usage, run 'gsnap --help'\n");
    /* print_program_usage(); */
    exit(9);
  }

  if (acc_fieldi_end < acc_fieldi_start) {
    fprintf(stderr,"--fastq-id-end must be equal to or greater than --fastq-id-start\n");
    exit(9);
  } else {
    Shortread_setup(acc_fieldi_start,acc_fieldi_end,filter_chastity_p);
  }

  if (novelsplicingp == true && knownsplicingp == true) {
    fprintf(stderr,"Novel splicing (-N) and known splicing (-s) both turned on => assume reads are RNA-Seq\n");
    pairmax = pairmax_rna;
    shortsplicedist_known = shortsplicedist;

  } else if (knownsplicingp == true) {
    fprintf(stderr,"Known splicing (-s) turned on => assume reads are RNA-Seq\n");
    pairmax = pairmax_rna;
    shortsplicedist_known = shortsplicedist;

  } else if (novelsplicingp == true) {
    fprintf(stderr,"Novel splicing (-N) turned on => assume reads are RNA-Seq\n");
    pairmax = pairmax_rna;
    shortsplicedist_known = 0;

  } else {
    /* Appears to be DNA-Seq */
    fprintf(stderr,"Neither novel splicing (-N) nor known splicing (-s) turned on => assume reads are DNA-Seq (genomic)\n");
    pairmax = pairmax_dna;
    shortsplicedist = shortsplicedist_known = 0U;
  }

  if (distantsplicing_penalty < localsplicing_penalty) {
    fprintf(stderr,"The distant splicing penalty %d cannot be less than local splicing penalty %d\n",
	    distantsplicing_penalty,localsplicing_penalty);
    exit(9);
  }

  if (fails_as_input_p == true && (sevenway_root == NULL && failsonlyp == false)) {
    fprintf(stderr,"The --fails-as-input option makes sense only with the --split-output or --failsonly option.  Turning it off.\n");
    fails_as_input_p = false;
  }

  if (sam_headers_batch >= 0) {
    if (part_modulus == sam_headers_batch) {
      sam_headers_p = true;
    } else {
      sam_headers_p = false;
    }
  }

  if (sam_read_group_id == NULL && sam_read_group_name != NULL) {
    sam_read_group_id = sam_read_group_name;
  } else if (sam_read_group_id != NULL && sam_read_group_name == NULL) {
    sam_read_group_name = sam_read_group_id;
  }

  if (chop_primers_p == true) {
    if (invert_first_p == false && invert_second_p == true) {
      /* orientation FR */
    } else {
      fprintf(stderr,"Adapter stripping not currently implemented for given orientation\n");
      exit(9);
    }
  }


  /* Open input stream and peek at first char */
  if (argc == 0) {
    input = stdin;
    files = (char **) NULL;
    nfiles = 0;
    nextchar = Shortread_input_init(input);
    if (nextchar == 0xFF) {
      fprintf(stderr,"Input appears to be a compact-reads file, which is not allowed as stdin.\n");
      exit(9);
    }
  } else {
    files = argv;
    nfiles = argc;

    if (gunzip_p == true) {
#ifdef HAVE_ZLIB
      if ((gzipped = gzopen(files[0],"rb")) == NULL) {
	fprintf(stderr,"Cannot open gzipped file %s\n",files[0]);
	exit(9);
      } else {
#ifdef HAVE_ZLIB_GZBUFFER
	gzbuffer(gzipped,GZBUFFER_SIZE);
#endif
	nextchar = Shortread_input_init_gzip(gzipped);
      }
#endif

    } else {
      if ((input = FOPEN_READ_TEXT(files[0])) == NULL) {
	fprintf(stderr,"Cannot open file %s\n",files[0]);
	exit(9);
      } else {
	nextchar = Shortread_input_init(input);
	if (nextchar == 0xFF) {
	  fclose(input);
	  input = (FILE *) NULL;
	  gobyreader = Goby_reader_new(files,nfiles,creads_window_start,creads_window_end,creads_complement_p);
	  creads_format_p = true;
	}
      }
    }

    files++;
    nfiles--;
  }

  /* Interpret first char to determine input type */
  if (nextchar == EOF) {
    fprintf(stderr,"Input is empty\n");
    exit(9);

#ifdef HAVE_GOBY
  } else if (creads_format_p == true) {
    if (user_quality_score_adj == false) {
      /* Use Goby default of 0, keeping Phred scores.  It is not
	 recommended that you override this value with -J x when
	 reading from Goby compact reads files. */
      MAPQ_init(/*quality_score_adj*/0);
      Pair_init(/*quality_score_adj*/0);
    }
    if (user_quality_shift == false) {
      /* By default, when outputting a non-Goby compact alignment
	 format (gsnap, sam), this will output Sanger quality scores,
	 equivalent to "-j 33".  If you prefer to output Illumina
	 quality scores, use "-j 64".  Goby compact alignment output
	 always uses Phred scores, ignoring this quality_shift value. */
      quality_shift = 33;
    }
#endif

  } else if (nextchar == '@') {
    /* Looks like a FASTQ file */
    if (nfiles == 0) {
#ifdef HAVE_ZLIB
      gzipped2 = (gzFile) NULL;
#endif
      input2 = (FILE *) NULL;
    } else {
      if (gunzip_p == true) {
#ifdef HAVE_ZLIB
	if ((gzipped2 = gzopen(files[0],"rb")) == NULL) {
	  fprintf(stderr,"Cannot open gzipped file %s\n",files[0]);
	  exit(9);
	} else {
#ifdef HAVE_ZLIB_GZBUFFER
	  gzbuffer(gzipped2,GZBUFFER_SIZE);
#endif
	  /* nextchar2 = */ Shortread_input_init_gzip(gzipped2);
	}
#endif
      } else {
	if ((input2 = FOPEN_READ_TEXT(files[0])) == NULL) {
	  fprintf(stderr,"Cannot open file %s\n",files[0]);
	  exit(9);
	} else {
	  /* nextchar2 = */ Shortread_input_init(input2);
	}
      }
      files++;
      nfiles--;
    }
    fastq_format_p = true;

  } else if (nextchar == '>') {
    /* Looks like a FASTA file */

  } else {
    fprintf(stderr,"First char is %c.  Expecting either '>' for FASTA or '@' for FASTQ format.\n",nextchar);
    exit(9);
  }


  /* Read in first batch of sequences */
  inbuffer = Inbuffer_new(nextchar,input,input2,
#ifdef HAVE_ZLIB
			  gzipped,gzipped2,
#endif
#ifdef HAVE_GOBY
			  gobyreader,
#endif
			  files,nfiles,fastq_format_p,creads_format_p,
			  barcode_length,invert_first_p,invert_second_p,chop_primers_p,
			  inbuffer_nspaces,inbuffer_maxchars,part_interval,part_modulus,
			  filter_if_both_p);

  nread = Inbuffer_fill_init(inbuffer);

  if (nread > 1) {
    multiple_sequences_p = true;
    if (offsetscomp_access != USE_ALLOCATE || genome_access != USE_ALLOCATE) {
      fprintf(stderr,"Note: >1 sequence detected, so index files are being memory mapped.\n");
      fprintf(stderr,"  GSNAP can run slowly at first while the computer starts to accumulate\n");
      fprintf(stderr,"  pages from the hard disk into its cache.  To copy index files into RAM\n");
      fprintf(stderr,"  instead of memory mapping, use -B 3, -B 4, or -B 5, if you have enough RAM.\n");
#ifdef HAVE_PTHREAD
      fprintf(stderr,"  For more speed, also try multiple threads (-t <int>), if you have multiple processors or cores.");
#endif
      fprintf(stderr,"\n");
    }
  } else {
    /* multiple_sequences_p = false; */
    /* fprintf(stderr,"Note: only 1 sequence detected.  Ignoring batch (-B) command\n"); */
    expand_offsets_p = false;
    offsetscomp_access = USE_MMAP_ONLY;
    positions_access = USE_MMAP_ONLY;
    genome_access = USE_MMAP_ONLY;
  }


  /* Prepare genomic data */

  genomesubdir = Datadir_find_genomesubdir(&fileroot,&dbversion,user_genomedir,dbroot);

  iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			    strlen(fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
  sprintf(iitfile,"%s/%s.chromosome.iit",genomesubdir,fileroot);
  if ((chromosome_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
				 /*divstring*/NULL,/*add_iit_p*/false,/*labels_read_p*/true)) == NULL) {
    fprintf(stderr,"IIT file %s is not valid\n",iitfile);
    exit(9);
  } else {
    nchromosomes = IIT_total_nintervals(chromosome_iit);
  }
  FREE(iitfile);


  if (snps_root == NULL) {
    genome = Genome_new(genomesubdir,fileroot,/*snps_root*/NULL,uncompressedp,genome_access);
    if (dibasep == true) {
      fprintf(stderr,"No longer supporting 2-base encoding\n");
      exit(9);
      if ((indexdb = Indexdb_new_genome(&index1part,genomesubdir,fileroot,/*idx_filesuffix*/"dibase",/*snps_root*/NULL,
					required_index1part,/*required_interval*/0,
					expand_offsets_p,offsetscomp_access,positions_access)) == NULL) {
	fprintf(stderr,"Cannot find offsets file %s.%s*offsets, needed for GSNAP color mode\n",fileroot,"dibase");
	exit(9);
      }
      indexdb2 = indexdb;
      print_ncolordiffs_p = true;

    } else if (mode == CMET_STRANDED || mode == CMET_NONSTRANDED) {
      if (user_cmetdir == NULL) {
	modedir = genomesubdir;
      } else {
	modedir = user_cmetdir;
      }

      if ((indexdb = Indexdb_new_genome(&index1part,modedir,fileroot,/*idx_filesuffix*/"metct",/*snps_root*/NULL,
					required_index1part,/*required_interval*/3,
					expand_offsets_p,offsetscomp_access,positions_access)) == NULL) {
	fprintf(stderr,"Cannot find metct index file.  Need to run cmetindex first\n");
	exit(9);
      }

      if ((indexdb2 = Indexdb_new_genome(&index1part,modedir,fileroot,/*idx_filesuffix*/"metga",/*snps_root*/NULL,
					 required_index1part,/*required_interval*/3,
					 expand_offsets_p,offsetscomp_access,positions_access)) == NULL) {
	fprintf(stderr,"Cannot find metga index file.  Need to run cmetindex first\n");
	exit(9);
      }

    } else if (mode == ATOI_STRANDED || mode == ATOI_NONSTRANDED) {
      if (user_atoidir == NULL) {
	modedir = genomesubdir;
      } else {
	modedir = user_atoidir;
      }

      if ((indexdb = Indexdb_new_genome(&index1part,modedir,fileroot,/*idx_filesuffix*/"a2iag",/*snps_root*/NULL,
					required_index1part,/*required_interval*/3,
					expand_offsets_p,offsetscomp_access,positions_access)) == NULL) {
	fprintf(stderr,"Cannot find a2iag index file.  Need to run atoiindex first\n");
	exit(9);
      }

      if ((indexdb2 = Indexdb_new_genome(&index1part,modedir,fileroot,/*idx_filesuffix*/"a2itc",/*snps_root*/NULL,
					 required_index1part,/*required_interval*/3,
					 expand_offsets_p,offsetscomp_access,positions_access)) == NULL) {
	fprintf(stderr,"Cannot find a2itc index file.  Need to run atoiindex first\n");
	exit(9);
      }


    } else {
      /* Standard behavior */
      if ((indexdb = Indexdb_new_genome(&index1part,genomesubdir,fileroot,IDX_FILESUFFIX,/*snps_root*/NULL,
					required_index1part,/*required_interval*/0,
					expand_offsets_p,offsetscomp_access,positions_access)) == NULL) {
	fprintf(stderr,"Cannot find offsets file %s.%s*offsets, needed for GSNAP\n",fileroot,IDX_FILESUFFIX);
	exit(9);
      }
      indexdb2 = indexdb;
    }

  } else {
    if (user_snpsdir == NULL) {
      snpsdir = genomesubdir;
      mapdir = Datadir_find_mapdir(/*user_mapdir*/NULL,genomesubdir,dbroot);
    } else {
      snpsdir = user_snpsdir;
      mapdir = user_snpsdir;
    }

    /* SNPs */
    genome = Genome_new(genomesubdir,fileroot,/*snps_root*/NULL,uncompressedp,genome_access);
    genomealt = Genome_new(snpsdir,fileroot,snps_root,uncompressedp,genome_access);

    if (dibasep == true) {
      fprintf(stderr,"Currently cannot combine SNPs with 2-base encoding\n");
      exit(9);
      print_ncolordiffs_p = true;

    } else if (mode == CMET_STRANDED || mode == CMET_NONSTRANDED) {
      if (user_cmetdir == NULL) {
	modedir = snpsdir;
      } else {
	modedir = user_cmetdir;
      }

      if ((indexdb = Indexdb_new_genome(&index1part,modedir,fileroot,/*idx_filesuffix*/"metct",snps_root,
					required_index1part,/*required_interval*/3,
					expand_offsets_p,offsetscomp_access,positions_access)) == NULL) {
	fprintf(stderr,"Cannot find metct index file.  Need to run cmetindex first\n");
	exit(9);
      }
      if ((indexdb2 = Indexdb_new_genome(&index1part,modedir,fileroot,/*idx_filesuffix*/"metga",snps_root,
					 required_index1part,/*required_interval*/3,
					 expand_offsets_p,offsetscomp_access,positions_access)) == NULL) {
	fprintf(stderr,"Cannot find metga index file.  Need to run cmetindex first\n");
	exit(9);
      }

    } else if (mode == ATOI_STRANDED || mode == ATOI_NONSTRANDED) {
      if (user_atoidir == NULL) {
	modedir = snpsdir;
      } else {
	modedir = user_atoidir;
      }

      if ((indexdb = Indexdb_new_genome(&index1part,modedir,fileroot,/*idx_filesuffix*/"a2iag",snps_root,
					required_index1part,/*required_interval*/3,
					expand_offsets_p,offsetscomp_access,positions_access)) == NULL) {
	fprintf(stderr,"Cannot find a2iag index file.  Need to run atoiindex first\n");
	exit(9);
      }
      if ((indexdb2 = Indexdb_new_genome(&index1part,modedir,fileroot,/*idx_filesuffix*/"a2itc",snps_root,
					 required_index1part,/*required_interval*/3,
					 expand_offsets_p,offsetscomp_access,positions_access)) == NULL) {
	fprintf(stderr,"Cannot find a2itc index file.  Need to run atoiindex first\n");
	exit(9);
      }

    } else {
      indexdb = Indexdb_new_genome(&index1part,snpsdir,fileroot,/*idx_filesuffix*/"ref",snps_root,
				   required_index1part,/*required_interval*/3,
				   expand_offsets_p,offsetscomp_access,positions_access);
      if (indexdb == NULL) {
	fprintf(stderr,"Cannot find snps index file for %s in directory %s\n",snps_root,snpsdir);
	exit(9);
      }
      indexdb2 = indexdb;
    }

    iitfile = (char *) CALLOC(strlen(mapdir)+strlen("/")+
			      strlen(snps_root)+1,sizeof(char));
    sprintf(iitfile,"%s/%s",mapdir,snps_root);
    fprintf(stderr,"Reading SNPs file %s/%s...",mapdir,snps_root);
    if ((snps_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
			     /*divstring*/NULL,/*add_iit_p*/true,/*labels_read_p*/true)) == NULL) {
      fprintf(stderr,"SNPs file %s.iit not found in %s.\n",snps_root,mapdir);
      if (user_snpsdir == NULL) {
	fprintf(stderr,"Available files:\n");
	Datadir_list_directory(stderr,genomesubdir);
	fprintf(stderr,"Either install file %s or specify a full directory path\n",snps_root);
	fprintf(stderr,"using the -D flag to gsnap.\n");
	exit(9);
      }
    }

    print_nsnpdiffs_p = true;
    snps_divint_crosstable = IIT_divint_crosstable(chromosome_iit,snps_iit);

    fprintf(stderr,"done\n");
    FREE(iitfile);
    if (user_snpsdir == NULL) {
      FREE(mapdir);
    }
  }
  genome_blocks = Genome_blocks(genome);

  if (min_distantsplicing_end_matches < index1part) {
    fprintf(stderr,"Minimum value for distant-splice-endlength is the value for -k (kmer size) %d\n",index1part);
    exit(9);
  }

  Dynprog_init(nullgap,EXTRAQUERYGAP,maxpeelback,extramaterial_end,extramaterial_paired);
  Compoundpos_init_positions_free(Indexdb_positions_fileio_p(indexdb));
  Spanningelt_init_positions_free(Indexdb_positions_fileio_p(indexdb));
  Stage1_init_positions_free(Indexdb_positions_fileio_p(indexdb));

  indexdb_size_threshold = (int) (10*Indexdb_mean_size(indexdb,mode,index1part));
  debug(printf("Size threshold is %d\n",indexdb_size_threshold));
  if (indexdb_size_threshold < MIN_INDEXDB_SIZE_THRESHOLD) {
    indexdb_size_threshold = MIN_INDEXDB_SIZE_THRESHOLD;
  }

  if (genes_file != NULL) {
    if ((genes_iit = IIT_read(genes_file,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
			      /*divstring*/NULL,/*add_iit_p*/true,/*labels_read_p*/true)) != NULL) {
      fprintf(stderr,"Reading genes file %s locally...",genes_file);
    } else {
      mapdir = Datadir_find_mapdir(/*user_mapdir*/NULL,genomesubdir,dbroot);
      iitfile = (char *) CALLOC(strlen(mapdir)+strlen("/")+
				strlen(genes_file)+1,sizeof(char));
      sprintf(iitfile,"%s/%s",mapdir,genes_file);
      if ((genes_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
				/*divstring*/NULL,/*add_iit_p*/true,/*labels_read_p*/true)) != NULL) {
	fprintf(stderr,"Reading genes file %s...",iitfile);
	FREE(iitfile);
	FREE(mapdir);
      } else {
	fprintf(stderr,"Genes file %s.iit not found locally or in %s.  Available files:\n",genes_file,mapdir);
	Datadir_list_directory(stderr,mapdir);
	fprintf(stderr,"Either install file %s or specify a full directory path\n",genes_file);
	exit(9);
      }
    }
    genes_divint_crosstable = IIT_divint_crosstable(chromosome_iit,genes_iit);
  }


  if (splicing_file != NULL) {
    if (user_splicingdir == NULL) {
      if ((splicing_iit = IIT_read(splicing_file,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
				   /*divstring*/NULL,/*add_iit_p*/true,/*labels_read_p*/true)) != NULL) {
	fprintf(stderr,"Reading splicing file %s locally...",splicing_file);
      }
    } else {
      iitfile = (char *) CALLOC(strlen(user_splicingdir)+strlen("/")+strlen(splicing_file)+1,sizeof(char));
      sprintf(iitfile,"%s/%s",user_splicingdir,splicing_file);
      if ((splicing_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
				   /*divstring*/NULL,/*add_iit_p*/true,/*labels_read_p*/true)) != NULL) {
	fprintf(stderr,"Reading splicing file %s...",iitfile);
	FREE(iitfile);
      }
    }

    if (splicing_iit == NULL) {
      mapdir = Datadir_find_mapdir(/*user_mapdir*/NULL,genomesubdir,dbroot);
      iitfile = (char *) CALLOC(strlen(mapdir)+strlen("/")+
				strlen(splicing_file)+1,sizeof(char));
      sprintf(iitfile,"%s/%s",mapdir,splicing_file);
      if ((splicing_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
				      /*divstring*/NULL,/*add_iit_p*/true,/*labels_read_p*/true)) != NULL) {
	fprintf(stderr,"Reading splicing file %s...",iitfile);
	FREE(iitfile);
	FREE(mapdir);
      } else {
	fprintf(stderr,"Splicing file %s.iit not found locally or in %s.  Available files:\n",splicing_file,mapdir);
	Datadir_list_directory(stderr,mapdir);
	fprintf(stderr,"Either install file %s or specify a full directory path\n",splicing_file);
	exit(9);
      }
    }

    splicing_divint_crosstable = IIT_divint_crosstable(chromosome_iit,splicing_iit);
    if ((donor_typeint = IIT_typeint(splicing_iit,"donor")) >= 0 && 
	(acceptor_typeint = IIT_typeint(splicing_iit,"acceptor")) >= 0) {
      fprintf(stderr,"found donor and acceptor tags, so treating as splicesites file\n");
      splicesites = Splicetrie_retrieve_via_splicesites(&distances_observed_p,&splicecomp,&splicetypes,&splicedists,
							&splicestrings,&splicefrags_ref,&splicefrags_alt,
							&nsplicesites,splicing_iit,splicing_divint_crosstable,
							donor_typeint,acceptor_typeint,chromosome_iit,
							genome,genomealt,shortsplicedist);
      if (nsplicesites == 0) {
	fprintf(stderr,"\nWarning: No splicesites observed for genome %s.  Are you sure this splicesite file was built for this genome?  Please compare chromosomes below:\n",
		dbroot);
	fprintf(stderr,"Chromosomes in the genome: ");
	IIT_dump_labels(stderr,chromosome_iit);
	fprintf(stderr,"Chromosomes in the splicesites IIT file: ");
	IIT_dump_divstrings(stderr,splicing_iit);
	exit(9);

      } else {
	Splicetrie_npartners(&nsplicepartners_skip,&nsplicepartners_obs,&nsplicepartners_max,splicesites,splicetypes,splicedists,
			     splicestrings,nsplicesites,chromosome_iit,shortsplicedist,distances_observed_p);
#if 0
	if (multiple_sequences_p == true && splicetrie_precompute_p == true) {
#endif
	  Splicetrie_build_via_splicesites(&triecontents_obs,&trieoffsets_obs,&triecontents_max,&trieoffsets_max,
					   nsplicepartners_skip,nsplicepartners_obs,nsplicepartners_max,splicetypes,
					   splicestrings,nsplicesites);
	  FREE(nsplicepartners_max);
	  FREE(nsplicepartners_obs);
	  FREE(nsplicepartners_skip);
	  Splicestring_gc(splicestrings,nsplicesites);
#if 0
	}
#endif
      }

    } else {
      fprintf(stderr,"no donor or acceptor tags found, so treating as introns file\n");
      splicesites = Splicetrie_retrieve_via_introns(&splicecomp,&splicetypes,&splicedists,
						    &splicestrings,&splicefrags_ref,&splicefrags_alt,
						    &nsplicesites,splicing_iit,splicing_divint_crosstable,
						    chromosome_iit,genome,genomealt);
      if (nsplicesites == 0) {
	fprintf(stderr,"\nWarning: No splicesites observed for genome %s.  Are you sure this splicesite file was built for this genome?  Please compare chromosomes below:\n",
		dbroot);
	fprintf(stderr,"Chromosomes in the genome: ");
	IIT_dump_labels(stderr,chromosome_iit);
	fprintf(stderr,"Chromosomes in the splicesites IIT file: ");
	IIT_dump_divstrings(stderr,splicing_iit);
	exit(9);
      } else {
#if 0
	if (multiple_sequences_p == true && splicetrie_precompute_p == true) {
#endif
	  Splicetrie_build_via_introns(&triecontents_obs,&trieoffsets_obs,splicesites,splicetypes,
				       splicestrings,nsplicesites,chromosome_iit,splicing_iit,splicing_divint_crosstable);
	  triecontents_max = (unsigned int *) NULL;
	  trieoffsets_max =  (unsigned int *) NULL;
	  Splicestring_gc(splicestrings,nsplicesites);
#if 0
	}
#endif
      }

    }

    /* For benchmarking purposes.  Can spend time/memory to load
       splicesites, but then not use them. */
    if (unloadp == true) {
      fprintf(stderr,"unloading...");

      if (nsplicesites > 0) {
	if (splicetrie_precompute_p == true) {
	  FREE(triecontents_max);
	  FREE(trieoffsets_max);
	  FREE(triecontents_obs);
	  FREE(trieoffsets_obs);
	} else {
	  FREE(nsplicepartners_max);
	  FREE(nsplicepartners_obs);
	  FREE(nsplicepartners_skip);
	  Splicestring_gc(splicestrings,nsplicesites);
	}
	FREE(splicefrags_ref);
	FREE(splicedists);
	FREE(splicetypes);
	FREE(splicesites);
	FREE(splicecomp);
	nsplicesites = 0;
      }

      FREE(splicing_divint_crosstable);
      IIT_free(&splicing_iit);
      splicing_iit = NULL;
      knownsplicingp = false;
      splicing_file = (char *) NULL;
    }

    fprintf(stderr,"done\n");
  }


  if (tally_root != NULL) {
    if (user_tallydir == NULL) {
      if ((tally_iit = IIT_read(tally_root,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
				/*divstring*/NULL,/*add_iit_p*/true,/*labels_read_p*/true)) != NULL) {
	fprintf(stderr,"Reading tally file %s.iit locally...",tally_root);
      }
    } else {
      iitfile = (char *) CALLOC(strlen(user_tallydir)+strlen("/")+strlen(tally_root)+1,sizeof(char));
      sprintf(iitfile,"%s/%s",user_tallydir,tally_root);
      if ((tally_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
				/*divstring*/NULL,/*add_iit_p*/true,/*labels_read_p*/true)) != NULL) {
	fprintf(stderr,"Reading tally file %s...",iitfile);
	FREE(iitfile);
      }
    }

    if (tally_iit == NULL) {
      iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+strlen(tally_root)+1,sizeof(char));
      sprintf(iitfile,"%s/%s",genomesubdir,tally_root);
      if ((tally_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
				/*divstring*/NULL,/*add_iit_p*/true,/*labels_read_p*/true)) != NULL) {
	fprintf(stderr,"Reading tally file %s...",iitfile);
	FREE(iitfile);
      } else {
	fprintf(stderr,"Tally file %s.iit not found locally",tally_root);
	if (user_tallydir != NULL) {
	  fprintf(stderr," or in %s",user_tallydir);
	}
	fprintf(stderr," or in %s\n",genomesubdir);
	exit(9);
      }
    }

    tally_divint_crosstable = IIT_divint_crosstable(chromosome_iit,tally_iit);
    fprintf(stderr,"done\n");
  }


  if (runlength_root != NULL) {
    if (user_runlengthdir == NULL) {
      if ((runlength_iit = IIT_read(runlength_root,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
				/*divstring*/NULL,/*add_iit_p*/true,/*labels_read_p*/true)) != NULL) {
	fprintf(stderr,"Reading runlength file %s.iit locally...",runlength_root);
      }
    } else {
      iitfile = (char *) CALLOC(strlen(user_runlengthdir)+strlen("/")+strlen(runlength_root)+1,sizeof(char));
      sprintf(iitfile,"%s/%s",user_runlengthdir,runlength_root);
      if ((runlength_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
				/*divstring*/NULL,/*add_iit_p*/true,/*labels_read_p*/true)) != NULL) {
	fprintf(stderr,"Reading runlength file %s...",iitfile);
	FREE(iitfile);
      }
    }

    if (runlength_iit == NULL) {
      iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+strlen(runlength_root)+1,sizeof(char));
      sprintf(iitfile,"%s/%s",genomesubdir,runlength_root);
      if ((runlength_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
				/*divstring*/NULL,/*add_iit_p*/true,/*labels_read_p*/true)) != NULL) {
	fprintf(stderr,"Reading runlength file %s...",iitfile);
	FREE(iitfile);
      } else {
	fprintf(stderr,"Runlength file %s.iit not found locally",runlength_root);
	if (user_runlengthdir != NULL) {
	  fprintf(stderr," or in %s",user_runlengthdir);
	}
	fprintf(stderr," or in %s\n",genomesubdir);
	exit(9);
      }
    }

    runlength_divint_crosstable = IIT_divint_crosstable(chromosome_iit,runlength_iit);
    fprintf(stderr,"done\n");
  }

  FREE(genomesubdir);
  FREE(fileroot);
  FREE(dbroot);


  Genome_setup(genome);
  Genome_hr_setup(Genome_blocks(genome),/*snp_blocks*/genomealt ? Genome_blocks(genomealt) : NULL,
		  query_unk_mismatch_p,genome_unk_mismatch_p,mode);
  Maxent_hr_setup(Genome_blocks(genome));
  Indexdb_setup(index1part);
  Indexdb_hr_setup(index1part);
  Oligo_setup(index1part);
  Splicetrie_setup(splicecomp,splicesites,splicefrags_ref,splicefrags_alt,
		   trieoffsets_obs,triecontents_obs,trieoffsets_max,triecontents_max,
		   /*snpp*/snps_iit ? true : false,amb_closest_p,amb_clip_p,min_shortend);
  Stage1hr_setup(index1part,chromosome_iit,nchromosomes,
		 genomealt,mode,terminal_threshold,
		 splicesites,splicetypes,splicedists,nsplicesites,
		 novelsplicingp,knownsplicingp,shortsplicedist_known,
		 nullgap,maxpeelback,maxpeelback_distalmedial,
		 extramaterial_end,extramaterial_paired,gmap_mode,
		 trigger_score_for_gmap,max_gmap_pairsearch,
		 max_gmap_terminal,max_gmap_improvement,antistranded_penalty);
  Substring_setup(print_nsnpdiffs_p,print_snplabels_p,
		  show_refdiff_p,snps_iit,snps_divint_crosstable,
		  genes_iit,genes_divint_crosstable,
		  splicing_iit,splicing_divint_crosstable,
		  donor_typeint,acceptor_typeint,trim_mismatch_score,
		  output_sam_p,mode);
  Dynprog_setup(splicing_iit,splicing_divint_crosstable,donor_typeint,acceptor_typeint,
		splicesites,splicetypes,splicedists,nsplicesites,
		trieoffsets_obs,triecontents_obs,trieoffsets_max,triecontents_max);
  Oligoindex_hr_setup(Genome_blocks(genome));
  Stage2_setup(/*splicingp*/novelsplicingp == true || knownsplicingp == true,
	       suboptimal_score_start,suboptimal_score_end);
  Pair_setup(trim_mismatch_score,trim_indel_score);
  Stage3_setup(/*splicingp*/novelsplicingp == true || knownsplicingp == true,
	       splicing_iit,splicing_divint_crosstable,donor_typeint,acceptor_typeint,
	       splicesites,min_intronlength,max_deletionlength,
	       expected_pairlength,pairlength_deviation);
  Stage3hr_setup(invert_first_p,invert_second_p,genes_iit,genes_divint_crosstable,
		 tally_iit,tally_divint_crosstable,runlength_iit,runlength_divint_crosstable,
		 distances_observed_p,pairmax,expected_pairlength,pairlength_deviation,
		 antistranded_penalty,favor_multiexon_p);
  SAM_setup(quiet_if_excessive_p,maxpaths);
  Goby_setup(show_refdiff_p);


  /* Setup outbuffer */
  if (output_goby_p == true) {
    if (goby_output_root == NULL) {
      fprintf(stderr,"--goby-output must be specified for Goby output.  For usage, run 'gsnap --help'\n");
      /* print_program_usage(); */
      exit(9);
    } else if (creads_format_p == false) {
      fprintf(stderr,"Currently can only write Goby if you read from compact reads files\n");
      exit(9);
    } else {
      gobywriter = Goby_writer_new(goby_output_root,"gsnap",PACKAGE_VERSION);
      Goby_writer_add_chromosomes(gobywriter,chromosome_iit);
    }
  }

  outbuffer = Outbuffer_new(output_buffer_size,nread,sevenway_root,
#ifdef USE_OLD_MAXENT
			    genome,
#endif
			    chromosome_iit,timingp,
			    output_sam_p,sam_headers_p,sam_read_group_id,sam_read_group_name,
			    sam_read_group_library,sam_read_group_platform,
			    gobywriter,nofailsp,failsonlyp,fails_as_input_p,
			    fastq_format_p,clip_overlap_p,merge_samechr_p,
			    maxpaths,quiet_if_excessive_p,quality_shift,
			    invert_first_p,invert_second_p,pairmax);

  Inbuffer_set_outbuffer(inbuffer,outbuffer);

  fprintf(stderr,"Starting alignment\n");
  stopwatch = Stopwatch_new();
  Stopwatch_start(stopwatch);

#ifndef HAVE_PTHREAD
  single_thread();
#else
  if (nworkers == 0) {
    single_thread();

  } else if (multiple_sequences_p == false) {
    single_thread();

  } else {
#ifdef WORKER_DETACH
    pthread_attr_init(&thread_attr_detach);
    if ((ret = pthread_attr_setdetachstate(&thread_attr_detach,PTHREAD_CREATE_DETACHED)) != 0) {
      fprintf(stderr,"ERROR: pthread_attr_setdetachstate returned %d\n",ret);
      exit(1);
    }
#endif
    pthread_attr_init(&thread_attr_join);
    if ((ret = pthread_attr_setdetachstate(&thread_attr_join,PTHREAD_CREATE_JOINABLE)) != 0) {
      fprintf(stderr,"ERROR: pthread_attr_setdetachstate returned %d\n",ret);
      exit(1);
    }
    
    worker_thread_ids = (pthread_t *) CALLOC(nworkers,sizeof(pthread_t));

    Except_init_pthread();
    pthread_key_create(&global_request_key,NULL);

    if (orderedp == true) {
      pthread_create(&output_thread_id,&thread_attr_join,Outbuffer_thread_ordered,
		     (void *) outbuffer);
    } else {
      pthread_create(&output_thread_id,&thread_attr_join,Outbuffer_thread_anyorder,
		     (void *) outbuffer);
    }
    for (worker_id = 0; worker_id < nworkers; worker_id++) {
#ifdef WORKER_DETACH
      pthread_create(&(worker_thread_ids[worker_id]),&thread_attr_detach,worker_thread,(void *) worker_id);
#else
      /* Need to have worker threads finish before we call Inbuffer_free() */
      pthread_create(&(worker_thread_ids[worker_id]),&thread_attr_join,worker_thread,(void *) worker_id);
#endif
    }
    
    pthread_join(output_thread_id,NULL);
    for (worker_id = 0; worker_id < nworkers; worker_id++) {
      pthread_join(worker_thread_ids[worker_id],NULL);
    }

    pthread_key_delete(global_request_key);
    /* Do not delete global_except_key, because worker threads might still need it */
    /* Except_term_pthread(); */

    FREE(worker_thread_ids);
  }
#endif /* HAVE_PTHREAD */

  runtime = Stopwatch_stop(stopwatch);
  Stopwatch_free(&stopwatch);

  nread = Outbuffer_nread(outbuffer);
  fprintf(stderr,"Processed %u queries in %.2f seconds (%.2f queries/sec)\n",
	  nread,runtime,(double) nread/runtime);

  Outbuffer_free(&outbuffer);
  Inbuffer_free(&inbuffer);	/* Also closes inputs, except for Goby */

  if (output_goby_p == true) {
    Goby_writer_finish(gobywriter,gobyreader);
    Goby_writer_free(&gobywriter);
  }
  if (creads_format_p == true) {
    Goby_reader_finish(gobyreader);
    Goby_reader_free(&gobyreader);
  }

#ifdef HAVE_GOBY
  /* Always call this, even if not using goby */
  Goby_shutdown();
#endif

  Dynprog_term();

  if (indexdb2 != indexdb) {
    Indexdb_free(&indexdb2);
  }
  if (indexdb != NULL) {
    Indexdb_free(&indexdb);
  }
  if (dbversion != NULL) {
    FREE(dbversion);
  }
  if (genomealt != NULL) {
    Genome_free(&genomealt);
    FREE(splicefrags_alt);	/* If genomealt == NULL, then splicefrags_alt == splicefrags_ref */
  }
  if (genome != NULL) {
    Genome_free(&genome);
  }

  if (nsplicesites > 0) {
    if (splicetrie_precompute_p == true) {
      FREE(triecontents_max);
      FREE(trieoffsets_max);
      FREE(triecontents_obs);
      FREE(trieoffsets_obs);
    } else {
      FREE(nsplicepartners_max);
      FREE(nsplicepartners_obs);
      FREE(nsplicepartners_skip);
      Splicestring_gc(splicestrings,nsplicesites);
    }
    FREE(splicefrags_ref);
    FREE(splicedists);
    FREE(splicetypes);
    FREE(splicecomp);
    FREE(splicesites);
  }

  if (runlength_iit != NULL) {
    FREE(runlength_divint_crosstable);
    IIT_free(&runlength_iit);
  }

  if (tally_iit != NULL) {
    FREE(tally_divint_crosstable);
    IIT_free(&tally_iit);
  }

  if (splicing_iit != NULL) {
    FREE(splicing_divint_crosstable);
    IIT_free(&splicing_iit);
  }

  if (genes_iit != NULL) {
    FREE(genes_divint_crosstable);
    IIT_free(&genes_iit);
  }

  if (snps_iit != NULL) {
    FREE(snps_divint_crosstable);
    IIT_free(&snps_iit);
  }

  if (chromosome_iit != NULL) {
    IIT_free(&chromosome_iit);
  }

  return 0;
}


static void
print_program_usage () {
  fprintf(stdout,"\
Usage: gsnap [OPTIONS...] <FASTA file>, or\n\
       cat <FASTA file> | gmap [OPTIONS...]\n\
\n\
");

  /* Input options */
  fprintf(stdout,"Input options (must include -d)\n");
  fprintf(stdout,"\
  -D, --dir=directory            Genome directory\n\
  -d, --db=STRING                Genome database\n\
  -k, --kmer=INT                 kmer size to use in genome database\n\
                                   (allowed values: 12-15).  If not specified, the program will\n\
                                   find the highest available kmer size in the genome database\n\
  -q, --part=INT/INT             Process only the i-th out of every n sequences\n\
                                   e.g., 0/100 or 99/100 (useful for distributing jobs\n\
                                   to a computer farm).\n\
  --input-buffer-size=INT        Size of input buffer (program reads this many sequences\n\
                                   at a time for efficiency) (default 1000)\n\
  --barcode-length=INT           Amount of barcode to remove from start of read\n\
                                   (default 0)\n\
  -o, --orientation=STRING       Orientation of paired-end reads\n\
                                   Allowed values: FR (fwd-rev, or typical Illumina; default),\n\
                                   RF (rev-fwd, for circularized inserts), or FF (fwd-fwd, same strand)\n\
  --fastq-id-start=INT           Starting position of identifier in FASTQ header, space-delimited (>= 1)\n\
  --fastq-id-end=INT             Ending position of identifier in FASTQ header, space-delimited (>= 1)\n\
                                 Examples:\n\
                                   @HWUSI-EAS100R:6:73:941:1973#0/1\n\
                                      start=1, end=1 (default) => identifier is HWUSI-EAS100R:6:73:941:1973#0\n\
                                   @SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=36\n\
                                      start=1, end=1  => identifier is SRR001666.1\n\
                                      start=2, end=2  => identifier is 071112_SLXA-EAS1_s_7:5:1:817:345\n\
                                      start=1, end=2  => identifier is SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345\n\
  --filter-chastity=STRING       Skips reads marked by the Illumina chastity program.  Expecting a string\n\
                                   after the accession having a 'Y' after the first colon, like this:\n\
                                         @accession 1:Y:0:CTTGTA\n\
                                   where the 'Y' signifies filtering by chastity.\n\
                                   Values: off (default), either, both.  For 'either', a 'Y' on either end\n\
                                   of a paired-end read will be filtered.  For 'both', a 'Y' is required\n\
                                   on both ends of a paired-end read (or on the only end of a single-end read).\n\
");
#ifdef HAVE_ZLIB
  fprintf(stdout,"\
  --gunzip                       Uncompress gzipped input files\n\
");
#endif
  fprintf(stdout,"\n");

  /* Computation options */
  fprintf(stdout,"Computation options\n");
  fprintf(stdout,"\
\n\
  Note: GSNAP has an ultrafast algorithm for calculating mismatches up to and including\n\
((readlength+2)/kmer - 2) (\"ultrafast mismatches\").  The program will run fastest if\n\
max-mismatches (plus suboptimal-levels) is within that value.\n\
Also, indels, especially end indels, take longer to compute, although the algorithm\n\
is still designed to be fast.\n\
\n\
");
#ifdef HAVE_MMAP
  fprintf(stdout,"\
  -B, --batch=INT                Batch mode (default = 2)\n\
                                 Mode     Offsets       Positions       Genome\n\
                                   0      allocate      mmap            mmap\n\
                                   1      allocate      mmap & preload  mmap\n\
                      (default)    2      allocate      mmap & preload  mmap & preload\n\
                                   3      allocate      allocate        mmap & preload\n\
                                   4      allocate      allocate        allocate\n\
                                   5      expand        allocate        allocate\n\
                           Note: For a single sequence, all data structures use mmap\n\
                           If mmap not available and allocate not chosen, then will use fileio (very slow)\n\
");

#endif
  fprintf(stdout,"\
  -m, --max-mismatches=FLOAT     Maximum number of mismatches allowed (if not specified, then\n\
                                   defaults to the ultrafast level of ((readlength+2)/12 - 2))\n\
                                   If specified between 0.0 and 1.0, then treated as a fraction\n\
                                   of each read length.  Otherwise, treated as an integral number\n\
                                   of mismatches (including indel and splicing penalties)\n\
                                   For RNA-Seq, you may need to increase this value slightly\n\
                                   to align reads extending past the ends of an exon.\n\
  --query-unk-mismatch=INT       Whether to count unknown (N) characters in the query as a mismatch\n\
                                   (0=no (default), 1=yes)\n\
  --genome-unk-mismatch=INT      Whether to count unknown (N) characters in the genome as a mismatch\n\
                                   (0=no, 1=yes (default))\n\
");

#if 0
  fprintf(stdout,"\
  -i, --indel-penalty-middle=INT Penalty for an indel in middle of read (default 1).\n\
                                   Counts against mismatches allowed.  To find indels, make\n\
                                   indel-penalty less than or equal to max-mismatches\n\
  -I, --indel-penalty-end=INT    Penalty for an indel at end of read (default 2).\n\
                                   Counts against mismatches allowed.  To find indels, make\n\
                                   indel-penalty less than or equal to max-mismatches\n\
");
#else
  fprintf(stdout,"\
  --terminal-threshold=INT       Threshold for searching for a terminal alignment (from one end of the\n\
                                   read to the best possible position at the other end) (default 2).\n\
                                   For example, if this value is 2, then if GSNAP finds an exact or\n\
                                   1-mismatch alignment, it will not try to find a terminal alignment.\n\
                                   Note that this default value may not be low enough if you want to\n\
                                   obtain terminal alignments for very short reads, although such reads\n\
                                   probably don't have enough specificity for terminal alignments anyway.\n\
                                   To turn off terminal alignments, set this to a high value, greater\n\
                                   than the value for --max-mismatches.\n\
  -i, --indel-penalty=INT        Penalty for an indel (default 2).\n\
                                   Counts against mismatches allowed.  To find indels, make\n\
                                   indel-penalty less than or equal to max-mismatches.\n\
                                   A value < 2 can lead to false positives at read ends\n\
");
#endif

#if 0
  /* No longer used */
  fprintf(stdout,"\
  -R, --masking=INT              Masking of frequent/repetitive oligomers to avoid spending time\n\
                                   on non-unique or repetitive reads\n\
                                   0 = no masking (will try to find non-unique or repetitive matches)\n\
                                   1 = mask frequent oligomers\n\
                                   2 = mask frequent and repetitive oligomers (fastest) (default)\n\
                                   3 = greedy frequent: mask frequent oligomers first, then\n\
                                       try no masking if alignments not found\n\
                                   4 = greedy repetitive: mask frequent and repetitive oligomers first, then\n\
                                       try no masking if alignments not found\n\
");
#endif

  fprintf(stdout,"\
  --indel-endlength=INT          Minimum length at end required for indel alignments (default 4)\n\
  -y, --max-middle-insertions=INT  Maximum number of middle insertions allowed (default 9)\n\
  -z, --max-middle-deletions=INT Maximum number of middle deletions allowed (default 30)\n\
  -Y, --max-end-insertions=INT   Maximum number of end insertions allowed (default 3)\n\
  -Z, --max-end-deletions=INT    Maximum number of end deletions allowed (default 6)\n\
  -M, --suboptimal-levels=INT    Report suboptimal hits beyond best hit (default 0)\n\
                                   All hits with best score plus suboptimal-levels are reported\n\
  -a, --adapter-strip=STRING     Method for removing adapters from reads.  Currently allowed values: off, paired.\n\
                                   Default is \"paired\", which removes adapters from paired-end reads if a\n\
                                   concordant or paired alignment cannot be found from the original read.\n\
                                   To turn off, use the value \"off\".\n\
  --trim-mismatch-score=INT      Score to use for mismatches when trimming at ends (default is -3;\n\
                                   to turn off trimming, specify 0).  Warning: turning trimming off\n\
                                   will give false positive mismatches at the ends of reads\n\
  --trim-indel-score=INT         Score to use for indels when trimming at ends (default is -4;\n\
                                   to turn off trimming, specify 0).  Warning: turning trimming off\n\
                                   will give false positive indels at the ends of reads\n\
  -V, --snpsdir=STRING           Directory for SNPs index files (created using snpindex) (default is\n\
                                   location of genome index files specified using -D and -d)\n \
  -v, --use-snps=STRING          Use database containing known SNPs (in <STRING>.iit, built\n\
                                   previously using snpindex) for tolerance to SNPs\n\
  -cmetdir=STRING                Directory for methylcytosine index files (created using cmetindex)\n\
                                   (default is location of genome index files specified using -D, -V, and -d)\n\
  -atoidir=STRING                Directory for A-to-I RNA editing index files (created using atoiindex)\n\
                                   (default is location of genome index files specified using -D, -V, and -d)\n\
  --mode=STRING                  Alignment mode: standard (default), cmet-stranded, cmet-nonstranded,\n\
                                    atoi-stranded, or atoi-nonstranded.  Non-standard modes requires you\n\
                                    to have previously run the cmetindex or atoiindex programs on the genome\n\
  --tallydir=STRING              Directory for tally IIT file to resolve concordant multiple results (default is\n\
                                   location of genome index files specified using -D and -d).  Note: can\n\
                                   just give full path name to --use-tally instead.\n\
  --use-tally=STRING             Use this tally IIT file to resolve concordant multiple results\n\
  --runlengthdir=STRING          Directory for runlength IIT file to resolve concordant multiple results (default is\n\
                                   location of genome index files specified using -D and -d).  Note: can\n\
                                   just give full path name to --use-runlength instead.\n\
  --use-runlength=STRING         Use this runlength IIT file to resolve concordant multiple results\n\
");


#if 0
  fprintf(stdout,"\
  -2, --dibase                   Input is 2-base encoded (e.g., SOLiD), with database built\n\
                                   previously using dibaseindex)\n\
");
#endif

#ifdef HAVE_PTHREAD
  fprintf(stdout,"\
  -t, --nthreads=INT             Number of worker threads\n\
");
#endif

  fprintf(stdout,"\n");
  fprintf(stdout,"Options for GMAP alignment within GSNAP\n");
  fprintf(stdout,"\
  --gmap-mode=STRING             Cases to use GMAP for complex alignments containing multiple splices or indels\n\
                                 Allowed values: none, pairsearch, terminal, improve (or multiple,\n\
                                    separated by commas).  Default: pairsearch,terminal,improve\n\
  --trigger-score-for-gmap=INT   Try GMAP pairsearch on nearby genomic regions if best score (the total\n\
                                   of both ends if paired-end) exceeds this value (default 5)\n \
  --max-gmap-pairsearch=INT      Perform GMAP pairsearch on nearby genomic regions up to this many\n\
                                   many candidate ends (default 3).  Requires pairsearch in --gmap-mode\n\
  --max-gmap-terminal=INT        Perform GMAP terminal on nearby genomic regions up to this many\n\
                                   candidate ends (default 3).  Requires terminal in --gmap-mode\n\
  --max-gmap-improvement=INT     Perform GMAP improvement on nearby genomic regions up to this many\n\
                                   candidate ends (default 3).  Requires improve in --gmap-mode\n\
  --microexon-spliceprob=FLOAT   Allow microexons only if one of the splice site probabilities is\n\
                                   greater than this value (default 0.90)\n\
");
  fprintf(stdout,"\n");


#if 0
  fprintf(stdout,"Genes options for RNA-Seq\n");
  fprintf(stdout,"\
  -g, --genes=STRING             Look for known genes in <STRING>.iit, to be used for resolving\n\
                                   multiple mapping reads.  See README instructions for the correct\n\
                                   formatting of a genes IIT file.\n\
  --favor-multiexon              In resolving multiple mapping reads, overlaps with known\n\
                                   multi-exon genes are favored over those with known single-exon\n\
                                   genes.  This favors spliced genes over psuedogenes.\n\
");
  fprintf(stdout,"\n");
#endif


  /* Splicing options */
  fprintf(stdout,"Splicing options for RNA-Seq\n");
  fprintf(stdout,"\
  -N, --novelsplicing=INT              Look for novel splicing (0=no (default), 1=yes)\n\
  --splicingdir=STRING                 Directory for splicing involving known sites or known introns,\n\
                                         as specified by the -s or --use-splicing flag (default is\n\
                                         directory computed from -D and -d flags).  Note: can\n\
                                         just give full pathname to the -s flag instead.\n\
  -s, --use-splicing=STRING            Look for splicing involving known sites or known introns\n\
                                         (in <STRING>.iit), at short or long distances\n\
                                         See README instructions for the distinction between known sites\n\
                                         and known introns\n\
  --ambig-splice-noclip                For ambiguous known splicing at ends of the read, do not clip at the\n\
                                         splice site, but extend instead into the intron.  This flag makes\n\
                                         sense only if you provide the --use-splicing flag, and you are trying\n\
                                         to eliminate all soft clipping with --trim-mismatch-score=0\n\
  -w, --localsplicedist=INT            Definition of local novel splicing event (default 200000)\n\
  -e, --local-splice-penalty=INT       Penalty for a local splice (default 0).  Counts against mismatches allowed\n\
  -E, --distant-splice-penalty=INT     Penalty for a distant splice (default 3).  A distant splice is one where\n\
                                         the intron length exceeds the value of -w, or --localsplicedist, or is an\n\
                                         inversion, scramble, or translocation between two different chromosomes\n\
                                         Counts against mismatches allowed\n\
  -K, --distant-splice-endlength=INT   Minimum length at end required for distant spliced alignments (default 16, min\n\
                                         allowed is the value of -k, or kmer size)\n\
  -l, --shortend-splice-endlength=INT  Minimum length at end required for short-end spliced alignments (default 2,\n\
                                         but unless known splice sites are provided with the -s flag, GSNAP may still\n\
                                         need the end length to be the value of -k, or kmer size to find a given splice\n\
  --distant-splice-identity=FLOAT      Minimum identity at end required for distant spliced alignments (default 0.95)\n\
  --antistranded-penalty=INT           Penalty for antistranded splicing when using stranded RNA-Seq protocols.\n\
                                         A positive value, such as 1, expects antisense on the first read\n\
                                         and sense on the second read.  Default is 0, which treats sense and antisense\n\
                                         equally well\n\
  --merge-distant-samechr              Report distant splices on the same chromosome as a single splice, if possible.\n\
                                         Will produce a single SAM line instead of two SAM lines, which is also done\n\
                                         for translocations, inversions, and scramble events\n\
");
  fprintf(stdout,"\n");


  /* Paired-end options */
  fprintf(stdout,"Options for paired-end reads\n");
  fprintf(stdout,"\
  --pairmax-dna=INT              Max total genomic length for DNA-Seq paired reads, or other reads\n\
                                   without splicing (default 1000).  Used if -N or -s is not specified.\n\
  --pairmax-rna=INT              Max total genomic length for RNA-Seq paired reads, or other reads\n\
                                   that could have a splice (default 200000).  Used if -N or -s is specified.\n\
                                   Should probably match the value for -w, --localsplicedist.\n\
  --pairexpect=INT               Expected paired-end length, used for calling splices in medial part of\n\
                                   paired-end reads (default 200)\n\
  --pairdev=INT                  Allowable deviation from expected paired-end length, used for\n\
                                   calling splices in medial part of paired-end reads (default 25)\n\
");
  fprintf(stdout,"\n");


  /* Quality score options */
  fprintf(stdout,"Options for quality scores\n");
  fprintf(stdout,"\
  --quality-protocol=STRING      Protocol for input quality scores.  Allowed values:\n\
                                   illumina (ASCII 64-126) (equivalent to -J 64 -j -31)\n\
                                   sanger   (ASCII 33-126) (equivalent to -J 33 -j 0)\n\
                                 Default is sanger (no quality print shift)\n\
                                 SAM output files should have quality scores in sanger protocol\n\
\n\
                                 Or you can customize this behavior with these flags:\n\
  -J, --quality-zero-score=INT   FASTQ quality scores are zero at this ASCII value\n\
                                   (default is 33 for sanger protocol; for Illumina, select 64)\n\
  -j, --quality-print-shift=INT  Shift FASTQ quality scores by this amount in output\n\
                                   (default is 0 for sanger protocol; to change Illumina input\n\
                                   to Sanger output, select -31)\n\
");

  /* Output options */
  fprintf(stdout,"Output options\n");
  fprintf(stdout,"\
  -n, --npaths=INT               Maximum number of paths to print (default 100).\n\
  -Q, --quiet-if-excessive       If more than maximum number of paths are found,\n\
                                   then nothing is printed.\n\
  -O, --ordered                  Print output in same order as input (relevant\n\
                                   only if there is more than one worker thread)\n\
  --show-refdiff                 For GSNAP output in SNP-tolerant alignment, shows all differences\n\
                                   relative to the reference genome as lower case (otherwise, it shows\n\
                                   all differences relative to both the reference and alternate genome)\n\
  --clip-overlap                 For paired-end reads whose alignments overlap, clip the overlapping region.\n\
  --print-snps                   Print detailed information about SNPs in reads (works only if -v also selected)\n\
                                   (not fully implemented yet)\n\
  --failsonly                    Print only failed alignments, those with no results\n\
  --nofails                      Exclude printing of failed alignments\n\
  --fails-as-input=STRING        Print completely failed alignments as input FASTA or FASTQ format\n\
                                   Allowed values: yes, no\n\
");

#ifdef HAVE_GOBY
  fprintf(stdout,"\
  -A, --format=STRING            Another format type, other than default.\n\
                                   Currently implemented: sam, goby\n\
");
#else
  fprintf(stdout,"\
  -A, --format=STRING            Another format type, other than default.\n\
                                   Currently implemented: sam\n\
                                   Also allowed, but not installed at compile-time: goby\n\
                                   (To install, need to re-compile with appropriate options)\n\
");
#endif

  fprintf(stdout,"\
  --split-output=STRING          Basename for multiple-file output, separately for nomapping,\n\
                                   halfmapping_uniq, halfmapping_mult, unpaired_uniq, unpaired_mult,\n\
                                   paired_uniq, paired_mult, concordant_uniq, and concordant_mult results (up to 9 files,\n\
                                   or 10 if --fails-as-input is selected, or 3 for single-end reads)\n\
  --output-buffer-size=INT       Buffer size, in queries, for output thread (default 1000).  When the number\n\
                                   of results to be printed exceeds this size, the worker threads are halted\n\
                                   until the backlog is cleared\n\
");
  fprintf(stdout,"\n");

  /* SAM options */
  fprintf(stdout,"Options for SAM output\n");
  fprintf(stdout,"\
  --no-sam-headers               Do not print headers beginning with '@'\n\
  --sam-headers-batch=INT        Print headers only for this batch, as specified by -q\n\
  --read-group-id=STRING         Value to put into read-group id (RG-ID) field\n\
  --read-group-name=STRING       Value to put into read-group name (RG-SM) field\n\
  --read-group-library=STRING    Value to put into read-group library (RG-LB) field\n\
  --read-group-platform=STRING   Value to put into read-group library (RG-PL) field\n\
");
  fprintf(stdout,"\n");


#ifdef HAVE_GOBY
  /* Goby options */
  fprintf(stdout,"Options for Goby library\n");
  fprintf(stdout,"\
  --goby-output=STRING           Basename for Goby output files\n\
  --creads-window-start=INT      Compact reads window start (default: 0=start of file)\n\
  --creads-window-end=INT        Compact reads window end (default: 0=end of file)\n\
  --creads-complement            Complement read sequences (without reversing)\n\
");
  fprintf(stdout,"\n");
#endif

  /* Help options */
  fprintf(stdout,"Help options\n");
  fprintf(stdout,"\
  --version                      Show version\n\
  --help                         Show this help message\n\
");
  return;
}

