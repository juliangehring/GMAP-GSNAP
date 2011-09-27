static char rcsid[] = "$Id: gsnap.c 37254 2011-03-28 16:34:08Z twu $";
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

#include "except.h"
#include "mem.h"
#include "bool.h"
#include "types.h"
#include "fopen.h"

#include "sequence.h"
#include "stopwatch.h"
#include "genome.h"
#include "indexdb_hr.h"
#include "mapq.h"
#include "substring.h"
#include "stage3hr.h"
#include "goby.h"
#include "spanningelt.h"
#include "splicetrie.h"
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
#include "getopt.h"


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/************************************************************************
 *   Global parameters
 ************************************************************************/

static IIT_T chromosome_iit = NULL;
static Indexdb_T indexdb = NULL;
static Indexdb_T indexdb2 = NULL; /* For cmet */
static Genome_T genome = NULL;
static Genome_T genomealt = NULL;
static UINT4 *genome_blocks = NULL;

static bool fastq_format_p = false;
static bool creads_format_p = false;
static bool pc_linefeeds_p = false;
static Stopwatch_T stopwatch = NULL;

/************************************************************************
 *   Program options
 ************************************************************************/

/* Input options */
static char *user_genomedir = NULL;
static char *user_snpsdir = NULL;
static char *dbroot = NULL;
static char *dbversion = NULL;
static int part_modulus = 0;
static int part_interval = 1;
static int barcode_length = 0;
static bool invert_first_p = false;
static bool invert_second_p = true;
static bool gunzip_p = false;

/* Compute options */
static bool chop_primers_p = false;
static bool novelsplicingp = false;
static int trim_mismatch_score = -3;
static Genomicpos_T expected_pairlength = 200;
static Genomicpos_T pairlength_deviation = 50;

static Access_mode_T offsets_access = USE_ALLOCATE;

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

#ifdef HAVE_PTHREAD
static pthread_t output_thread_id, *worker_thread_ids;
static int nworkers = 1;	/* (int) sysconf(_SC_NPROCESSORS_ONLN) */
#endif

static Masktype_T masktype = MASK_REPETITIVE;
static int subopt_levels = 0;

/* If negative, then hasn't been specified by user.  If between 0 and
   1, then treated as a fraction of the querylength.  Else, treated as
   an integer */
static double usermax_level_float = -1.0;

static int terminal_penalty = 1;
static int max_terminal_length = 16; /* Should be 14 or more because need to have 12-mer q3 */

static int indel_penalty = 1;
static bool allow_end_indels_p = true;
static int max_middle_insertions = 9;
static int max_middle_deletions = 30;
static int max_end_insertions = 3;
static int max_end_deletions = 6;
static int min_indel_end_matches = 3;
static Genomicpos_T shortsplicedist = 200000;
static int localsplicing_penalty = 0;
static int distantsplicing_penalty = 3;
static int min_localsplicing_end_matches = 14;
static int min_distantsplicing_end_matches = 16;
static double min_distantsplicing_identity = 0.95;
static int min_shortend = 2;
static bool find_novel_doublesplices_p = false;
static int indexdb_size_threshold;


/* Splicesites IIT */
static bool knownsplicingp = false;
static bool distances_observed_p = false;
static char *splicesites_file = (char *) NULL;
static IIT_T splicesites_iit = NULL;

static int donor_typeint = -1;		/* for splicesites_iit */
static int acceptor_typeint = -1;	/* for splicesites_iit */

static int *splicesites_divint_crosstable = NULL;
static Genomicpos_T *splicesites = NULL;
static Splicetype_T *splicetypes = NULL;
static Genomicpos_T *splicedists = NULL; /* maximum observed splice distance for given splice site */
static int *nsplicepartners_skip = NULL;
static int *nsplicepartners_obs = NULL;
static int *nsplicepartners_max = NULL;
static List_T *splicestrings = NULL;
static UINT4 *splicefrags_ref = NULL;
static UINT4 *splicefrags_alt = NULL;
static int nsplicesites = 0;

static bool splicetrie_precompute_p = true;
static unsigned int *trieoffsets_obs = NULL;
static unsigned int *triecontents_obs = NULL;
static unsigned int *trieoffsets_max = NULL;
static unsigned int *triecontents_max = NULL;


/* SNPs IIT */
static bool dibasep = false;
static char *user_cmetdir = NULL;
static bool cmetp = false;
static char *snps_root = (char *) NULL;
static IIT_T snps_iit = NULL;
static int *snps_divint_crosstable = NULL;


/* Output options */
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

/* SAM */
static int sam_headers_batch = -1;
static bool sam_headers_p = true;
static char *sam_read_group_id = NULL;
static char *sam_read_group_name = NULL;

/* Goby */
static char *goby_output_root = NULL;
static unsigned long creads_window_start = 0;
static unsigned long creads_window_end = 0;
static Gobyreader_T gobyreader = NULL;
static Gobywriter_T gobywriter = NULL;

/* Input/output */
static char *sevenway_root = NULL;
static Outbuffer_T outbuffer;
static Inbuffer_T inbuffer;
static int inbuffer_nspaces = 1000;
static unsigned int inbuffer_maxchars = -1U; /* Currently not used by Inbuffer_T */



/* Alignment options */
static bool uncompressedp = false;

/* getopt used alphabetically: AaBCcDdEeGIiJjKklMmNnOopQqRSstVvwYyZz27 */

static struct option long_options[] = {
  /* Input options */
  {"dir", required_argument, 0, 'D'},	/* user_genomedir */
  {"db", required_argument, 0, 'd'}, /* dbroot */
  {"genomefull", no_argument, 0, 'G'}, /* uncompressedp */
  {"part", required_argument, 0, 'q'}, /* part_modulus, part_interval */
  {"orientation", required_argument, 0, 'o'}, /* invert_first_p, invert_second_p */
  {"input-buffer-size", required_argument, 0, 0}, /* inbuffer_nspaces */
  {"barcode-length", required_argument, 0, 0},	  /* barcode_length */
  {"pc-linefeeds", no_argument, 0, 0},		  /* pc_linefeeds_p */

#ifdef HAVE_ZLIB
  {"gunzip", no_argument, 0, 0}, /* gunzip_p */
#endif

  /* Compute options */
#ifdef HAVE_MMAP
  {"batch", required_argument, 0, 'B'}, /* offsets_access, positions_access, genome_access */
#endif
  {"pairmax-dna", required_argument, 0, 0}, /* pairmax_dna */
  {"pairmax-rna", required_argument, 0, 0}, /* pairmax_rna */
  {"pairexpect", required_argument, 0, 'p'}, /* expected_pairlength */
  {"pairdev", required_argument, 0, 0}, /* pairlength_deviation */
#ifdef HAVE_PTHREAD
  {"nthreads", required_argument, 0, 't'}, /* nworkers */
#endif
  {"adapter-strip", required_argument, 0, 'a'},	/* chop_primers_p */
  {"trim-mismatch-score", required_argument, 0, 0}, /* trim_mismatch_score */
  {"novelsplicing", required_argument, 0, 'N'}, /* novelsplicingp */
  {"novel-doublesplices", no_argument, 0, 0},	/* find_novel_doublesplices_p */

  {"max-mismatches", required_argument, 0, 'm'}, /* usermax_level_float */
  {"terminal-penalty", required_argument, 0, 0}, /* terminal_penalty */

  {"indel-penalty", required_argument, 0, 'i'}, /* indel_penalty */

  {"indel-endlength", required_argument, 0, 'I'}, /* min_indel_end_matches, allow_end_indels_p */

  {"max-middle-insertions", required_argument, 0, 'y'}, /* max_middle_insertions */
  {"max-middle-deletions", required_argument, 0, 'z'}, /* max_middle_deletions */
  {"max-end-insertions", required_argument, 0, 'Y'}, /* max_end_insertions */
  {"max-end-deletions", required_argument, 0, 'Z'}, /* max_end_deletions */
  {"suboptimal-levels", required_argument, 0, 'M'}, /* subopt_levels */
  {"masking", required_argument, 0, 'R'}, /* masktype */

  {"localsplicedist", required_argument, 0, 'w'}, /* shortsplicedist */
  {"splicesites", required_argument, 0, 's'}, /* splicesites_iit, knownsplicingp */
  {"splicetrie-precompute", required_argument, 0, 'S'}, /* splicetrie_precompute_p */

  {"local-splice-penalty", required_argument, 0, 'e'}, /* localsplicing_penalty */
  {"distant-splice-penalty", required_argument, 0, 'E'}, /* distantsplicing_penalty */
  {"local-splice-endlength", required_argument, 0, 'k'}, /* min_localsplicing_end_matches */
  {"distant-splice-endlength", required_argument, 0, 'K'}, /* min_distantsplicing_end_matches */
  {"shortend-splice-endlength", required_argument, 0, 'l'}, /* min_shortend */

  {"cmetdir", required_argument, 0, 'C'}, /* user_cmetdir */
  {"use-cmet", no_argument, 0, 'c'}, /* cmetp */

  {"dibase", no_argument, 0, '2'}, /* dibasep */
  {"snpsdir", required_argument, 0, 'V'},   /* user_snpsdir */
  {"use-snps", required_argument, 0, 'v'}, /* snps_root */

  /* Output options */
  {"format", required_argument, 0, 'A'}, /* output_sam_p, output_goby_p */

  {"quality-protocol", required_argument, 0, 0}, /* quality_score_adj, quality_shift */
  {"quality-zero-score", required_argument, 0, 'J'}, /* quality_score_adj */
  {"quality-print-shift", required_argument, 0, 'j'}, /* quality_shift */
  {"sam-headers-batch", required_argument, 0, 0},	/* sam_headers_batch */
  {"no-sam-headers", no_argument, 0, 0},	/* sam_headers_p */
  {"read-group-id", required_argument, 0, 0},	/* sam_read_group_id */
  {"read-group-name", required_argument, 0, 0},	/* sam_read_group_name */

  {"noexceptions", no_argument, 0, '0'}, /* exception_raise_p */
  {"npaths", required_argument, 0, 'n'}, /* maxpaths */
  {"quiet-if-excessive", no_argument, 0, 'Q'}, /* quiet_if_excessive_p */
  {"ordered", no_argument, 0, 'O'}, /* orderedp */
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
#endif

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
process_request (Request_T request, Floors_T *floors_array) {
  int jobid;
  Shortread_T queryseq1, queryseq2;
  Stage3_T *stage3array, *stage3array5, *stage3array3;
  Stage3pair_T *stage3pairarray;

  int npaths, npaths5, npaths3;
  bool concordantp;

  jobid = Request_id(request);
  queryseq1 = Request_queryseq1(request);
  queryseq2 = Request_queryseq2(request);

  /* printf("%s\n",Shortread_accession(queryseq1)); */

  if (queryseq2 == NULL) {
    stage3array = Stage1_single_read(&npaths,queryseq1,indexdb,indexdb2,indexdb_size_threshold,chromosome_iit,
				     genome,genomealt,floors_array,knownsplicingp,novelsplicingp,/*canonicalp*/true,
				     maxpaths,/*maxchimerapaths*/maxpaths,
				     usermax_level_float,subopt_levels,masktype,
				     terminal_penalty,max_terminal_length,indel_penalty,
				     max_middle_insertions,max_middle_deletions,
				     allow_end_indels_p,max_end_insertions,max_end_deletions,min_indel_end_matches,
				     shortsplicedist,
				     localsplicing_penalty,distantsplicing_penalty,min_localsplicing_end_matches,
				     min_distantsplicing_end_matches,min_distantsplicing_identity,min_shortend,
				     find_novel_doublesplices_p,splicesites,splicetypes,splicedists,
				     splicefrags_ref,splicefrags_alt,nsplicesites,nsplicepartners_skip,
				     trieoffsets_obs,triecontents_obs,nsplicepartners_obs,
				     trieoffsets_max,triecontents_max,nsplicepartners_max,
				     splicestrings,dibasep,cmetp);
    return Result_single_read_new(jobid,(void **) stage3array,npaths);

  } else if ((stage3pairarray = Stage1_paired_read(&npaths,&concordantp,&stage3array5,&npaths5,&stage3array3,&npaths3,
						   queryseq1,queryseq2,indexdb,indexdb2,indexdb_size_threshold,
						   chromosome_iit,genome,genomealt,floors_array,
						   knownsplicingp,novelsplicingp,/*canonicalp*/true,
						   /*maxpaths*/maxpaths,/*maxchimerapaths*/maxpaths,
						   usermax_level_float,subopt_levels,masktype,
						   terminal_penalty,max_terminal_length,indel_penalty,
						   max_middle_insertions,max_middle_deletions,
						   allow_end_indels_p,max_end_insertions,max_end_deletions,min_indel_end_matches,
						   shortsplicedist,
						   localsplicing_penalty,distantsplicing_penalty,min_localsplicing_end_matches,
						   min_distantsplicing_end_matches,min_distantsplicing_identity,min_shortend,
						   find_novel_doublesplices_p,splicesites,splicetypes,splicedists,
						   splicefrags_ref,splicefrags_alt,nsplicesites,nsplicepartners_skip,
						   trieoffsets_obs,triecontents_obs,nsplicepartners_obs,
						   trieoffsets_max,triecontents_max,nsplicepartners_max,
						   splicestrings,pairmax,expected_pairlength,pairlength_deviation,
						   dibasep,cmetp)) == NULL) {
    /* No paired or concordant hits found */
    return Result_paired_as_singles_new(jobid,(void **) stage3array5,npaths5,(void **) stage3array3,npaths3);

  } else {
    /* Paired or concordant hits found */
    return Result_paired_read_new(jobid,(void **) stage3pairarray,npaths,concordantp);
  }

}



#ifdef HAVE_SIGACTION
static const Except_T sigfpe_error = {"SIGFPE--arithmetic exception"};
static const Except_T sigsegv_error = {"SIGSEGV--segmentation violation"};
static const Except_T sigtrap_error = {"SIGTRAP--hardware fault"};
static const Except_T misc_signal_error = {"Miscellaneous signal"};

static void
signal_handler (int sig) {
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
single_thread () {
  Floors_T *floors_array;
  Request_T request;
  Result_T result;
  Shortread_T queryseq1;
  int i;

  floors_array = (Floors_T *) CALLOC(MAX_QUERYLENGTH+1,sizeof(Floors_T));
  /* Except_stack_create(); -- requires pthreads */

  while ((request = Inbuffer_get_request(inbuffer)) != NULL) {
    debug(printf("single_thread got request %d\n",Request_id(request)));
#ifdef LEAKCHECK
    Mem_leak_check_start(__FILE__,__LINE__);
#endif

    TRY
      result = process_request(request,floors_array);
    ELSE
      queryseq1 = Request_queryseq1(request);
      if (queryseq1 == NULL) {
	fprintf(stderr,"NULL");
      } else if (Shortread_accession(queryseq1) == NULL) {
	fprintf(stderr,"unnamed (%d bp)",Shortread_fulllength(queryseq1));
      } else {
	fprintf(stderr,"%s (%d bp)",Shortread_accession(queryseq1),Shortread_fulllength(queryseq1));
      }
      fprintf(stderr,"\n");
      Shortread_print_oneline(stderr,queryseq1);
      fprintf(stderr,"\n");
      fprintf(stderr,"To obtain a core dump, re-run program on problem sequence with the -0 [zero] flag\n");

      fprintf(stderr,"Exiting...\n");
      exit(9);
    RERAISE;
    END_TRY;

    Outbuffer_print_result(outbuffer,result,request);
    Result_free(&result);
#ifdef LEAKCHECK
    Mem_leak_check_end(__FILE__,__LINE__);
#endif

    /* Allocated by fill_buffer in Inbuffer_get_request */
    Request_free(&request);
  }

  /* Except_stack_destroy(); -- requires pthreads */

  for (i = 0; i <= MAX_QUERYLENGTH; i++) {
    if (floors_array[i] != NULL) {
      Floors_free(&(floors_array[i]));
    }
  }
  FREE(floors_array);

  return;
}



static void *
worker_thread (void *data) {
  Floors_T *floors_array;
  Request_T request;
  Result_T result;
  Shortread_T queryseq1;
  int i;

  /* Thread-specific data and storage */
  floors_array = (Floors_T *) CALLOC(MAX_QUERYLENGTH+1,sizeof(Floors_T));
  Except_stack_create();

  while ((request = Inbuffer_get_request(inbuffer)) != NULL) {
    debug(printf("worker_thread got request %d\n",Request_id(request)));

    TRY
      result = process_request(request,floors_array);
    ELSE
      queryseq1 = Request_queryseq1(request);
      if (queryseq1 == NULL) {
	fprintf(stderr,"NULL");
      } else if (Shortread_accession(queryseq1) == NULL) {
	fprintf(stderr,"unnamed (%d bp)",Shortread_fulllength(queryseq1));
      } else {
	fprintf(stderr,"%s (%d bp)",Shortread_accession(queryseq1),Shortread_fulllength(queryseq1));
      }
      fprintf(stderr,"\n");
      Shortread_print_oneline(stderr,queryseq1);
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

  Except_stack_destroy();

  for (i = 0; i <= MAX_QUERYLENGTH; i++) {
    if (floors_array[i] != NULL) {
      Floors_free(&(floors_array[i]));
    }
  }
  FREE(floors_array);

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


int
main (int argc, char *argv[]) {
  char *genomesubdir = NULL, *snpsdir = NULL, *cmetdir = NULL, *mapdir = NULL, *iitfile = NULL, *fileroot = NULL;
  FILE *input = NULL, *input2 = NULL;
#ifdef HAVE_ZLIB
  gzFile gzipped = NULL, gzipped2 = NULL;
#endif

  bool multiple_sequences_p = false;
  char **files;
  int nfiles, nextchar = '\0', i;

  int nread;
  double runtime;

#ifdef HAVE_PTHREAD
  int ret;
  pthread_attr_t thread_attr_detach, thread_attr_join;
#endif

  int opt, c;
  extern int optind;
  extern char *optarg;
  int long_option_index = 0;
  const char *long_name;
  char **argstart;

#ifdef HAVE_SIGACTION
  struct sigaction signal_action;
#endif

  while ((opt = getopt_long(argc,argv,
			    "D:d:Gq:o:a:N:R:M:m:i:I:y:Y:z:Z:w:E:e:J:K:k:l:s:S:2V:v:B:p:t:A:j:0n:QC:cO",
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
      } else if (!strcmp(long_name,"input-buffer-size")) {
	inbuffer_nspaces = atoi(optarg);
      } else if (!strcmp(long_name,"barcode-length")) {
	barcode_length = atoi(optarg);
      } else if (!strcmp(long_name,"pc-linefeeds")) {
	pc_linefeeds_p = true;
#ifdef HAVE_ZLIB
      } else if (!strcmp(long_name,"gunzip")) {
	gunzip_p = true;
#endif
      } else if (!strcmp(long_name,"split-output")) {
	sevenway_root = optarg;
      } else if (!strcmp(long_name,"fails-as-input")) {
	fails_as_input_p = true;
      } else if (!strcmp(long_name,"pairmax-dna")) {
	pairmax_dna = atoi(optarg);
      } else if (!strcmp(long_name,"pairmax-rna")) {
	pairmax_rna = atoi(optarg);
      } else if (!strcmp(long_name,"pairdev")) {
	pairlength_deviation = atoi(optarg);
      } else if (!strcmp(long_name,"novel-doublesplices")) {
	find_novel_doublesplices_p = true;
      } else if (!strcmp(long_name,"terminal-penalty")) {
	terminal_penalty = atoi(optarg);
      } else if (!strcmp(long_name,"trim-mismatch-score")) {
	trim_mismatch_score = atoi(optarg);
      } else if (!strcmp(long_name,"show-refdiff")) {
	show_refdiff_p = true;
      } else if (!strcmp(long_name,"no-sam-headers")) {
	sam_headers_p = false;
      } else if (!strcmp(long_name,"sam-headers-batch")) {
	sam_headers_batch = atoi(optarg);
      } else if (!strcmp(long_name,"quality-protocol")) {
	if (user_quality_score_adj == true) {
	  fprintf(stderr,"Cannot specify both -J (--quality-zero-score) and --quality-protocol\n");
	  exit(9);
	} else if (user_quality_shift == true) {
	  fprintf(stderr,"Cannot specify both -j (--quality-print-shift) and --quality-protocol\n");
	  exit(9);
	} else if (!strcmp(optarg,"illumina")) {
	  MAPQ_init(/*quality_score_adj*/64);
	  user_quality_score_adj = true;
	  quality_shift = -31;
	  user_quality_shift = true;
	} else if (!strcmp(optarg,"sanger")) {
	  MAPQ_init(/*quality_score_adj*/33);
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
      } else if (!strcmp(long_name,"goby-output")) {
	goby_output_root = optarg;
      } else if (!strcmp(long_name,"distant-splice-identity")) {
	min_distantsplicing_identity = atof(optarg);
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
	creads_window_start = strtoul(optarg,NULL,10);
      } else if (!strcmp(long_name,"creads-window-end")) {
	creads_window_end = strtoul(optarg,NULL,10);
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
      } else {
	fprintf(stderr,"Currently allowed values for adapter stripping (-a): paired\n");
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

    case 'M': subopt_levels = atoi(optarg); break;
    case 'm':
      usermax_level_float = atof(optarg);
      if (usermax_level_float > 1.0 && usermax_level_float != rint(usermax_level_float)) {
	fprintf(stderr,"Cannot specify fractional value %f for --max-mismatches except between 0.0 and 1.0\n",usermax_level_float);
	exit(9);
      }
      break;

    case 'i': indel_penalty = atoi(optarg); break;
    case 'I':
      min_indel_end_matches = atoi(optarg);
      if (min_indel_end_matches > 12) {
	allow_end_indels_p = false;
      }
      break;

    case 'y': max_middle_insertions = atoi(optarg); break;
    case 'Y': max_end_insertions = atoi(optarg); break;
    case 'z': max_middle_deletions = atoi(optarg); break;
    case 'Z': max_end_deletions = atoi(optarg); break;

    case 'w': shortsplicedist = strtoul(optarg,NULL,10); break;

    case 'E': distantsplicing_penalty = atoi(optarg); break;
    case 'e': localsplicing_penalty = atoi(optarg); break;
    case 'K': 
      min_distantsplicing_end_matches = atoi(optarg); 
      if (min_distantsplicing_end_matches < 14) {
	fprintf(stderr,"Minimum value for distant-splice-endlength is 14\n");
	exit(9);
      }
      break;
    case 'k': 
      min_localsplicing_end_matches = atoi(optarg); 
      if (min_localsplicing_end_matches < 14) {
	fprintf(stderr,"Minimum value for local-splice-endlength is 14\n");
	exit(9);
      }
      break;
    case 'l':  min_shortend = atoi(optarg); break;

    case 's': splicesites_file = optarg; knownsplicingp = true; break;
    case 'S': 
      if (!strcmp(optarg,"0")) {
	splicetrie_precompute_p = false;
      } else if (!strcmp(optarg,"1")) {
	splicetrie_precompute_p = true;
      } else {
	fprintf(stderr,"SNP printing mode %s not recognized.\n",optarg);
	exit(9);
      }
      break;

    case '2': dibasep = true; break;
    case 'C': user_cmetdir = optarg; break;
    case 'c': cmetp = true; break;
    case 'V': user_snpsdir = optarg; break;
    case 'v': snps_root = optarg; break;

    case 'B':
      if (!strcmp(optarg,"0")) {
	offsets_access = USE_ALLOCATE; /* was batch_offsets_p = true */
	positions_access = USE_MMAP_ONLY; /* was batch_positions_p = false */
	genome_access = USE_MMAP_ONLY; /* was batch_genome_p = false */
      } else if (!strcmp(optarg,"1")) {
	offsets_access = USE_ALLOCATE; /* was batch_offsets_p = true */
	positions_access = USE_MMAP_PRELOAD; /* was batch_positions_p = true */
	genome_access = USE_MMAP_ONLY; /* was batch_genome_p = false */
      } else if (!strcmp(optarg,"2")) {
	offsets_access = USE_ALLOCATE; /* was batch_offsets_p = true */
	positions_access = USE_MMAP_PRELOAD; /* was batch_positions_p = true */
	genome_access = USE_MMAP_PRELOAD; /* was batch_genome_p = true */
      } else if (!strcmp(optarg,"3")) {
	offsets_access = USE_ALLOCATE;
	positions_access = USE_ALLOCATE;
	genome_access = USE_MMAP_PRELOAD; /* was batch_genome_p = true */
      } else if (!strcmp(optarg,"4")) {
	offsets_access = USE_ALLOCATE;
	positions_access = USE_ALLOCATE;
	genome_access = USE_ALLOCATE;
      } else {
	fprintf(stderr,"Batch mode %s not recognized.  Only allow 0-4.  Run 'gsnap --help' for more information.\n",optarg);
	exit(9);
      }
      break;
    case 'p': expected_pairlength = atoi(optarg); break;
#ifdef HAVE_PTHREAD
    case 't': nworkers = atoi(optarg); break;
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
	quality_shift = atoi(optarg);
	user_quality_shift = true;
      }
      break;

    case 'J':
      if (user_quality_score_adj == true) {
	fprintf(stderr,"Cannot specify both -J (--quality-zero-score) and --quality-protocol\n");
	exit(9);
      } else {
	MAPQ_init(/*quality_score_adj*/atoi(optarg));
	user_quality_score_adj = true;
      }
      break;

    case '0': exception_raise_p = false; break; /* Allows signals to pass through */
    case 'n': maxpaths = atoi(optarg); break;
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

  if (novelsplicingp == true || knownsplicingp == true) {
    /* Appears to be RNA-Seq */
    pairmax = pairmax_rna;
  } else {
    pairmax = pairmax_dna;
  }

  if (fails_as_input_p == true && (sevenway_root == NULL && failsonlyp == false)) {
    fprintf(stderr,"The --fails-as-input option makes sense only with the --7-way-output or --failsonly option.  Turning it off.\n");
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


#ifdef MEMUSAGE
  Mem_usage_init();
  nworkers = 0;
  fprintf(stderr,"For memusage, setting to 0 threads\n");
#endif

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
	  gobyreader = Goby_reader_new(files,nfiles,creads_window_start,creads_window_end);
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
      /* Use Goby default of 0 */
      MAPQ_init(/*quality_score_adj*/0);
    }
    if (user_quality_shift == false) {
      quality_shift = 64;
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


  fprintf(stderr,"GSNAP version %s called with args:",PACKAGE_VERSION);
  argstart = &(argv[-optind]);
  for (c = 1; c < argc + optind; c++) {
    fprintf(stderr," %s",argstart[c]);
  }
  fprintf(stderr,"\n");


  /* Read in first batch of sequences */
  inbuffer = Inbuffer_new(nextchar,input,input2,
#ifdef HAVE_ZLIB
			  gzipped,gzipped2,
#endif
#ifdef HAVE_GOBY
			  gobyreader,
#endif
			  files,nfiles,fastq_format_p,creads_format_p,pc_linefeeds_p,
			  barcode_length,invert_first_p,invert_second_p,chop_primers_p,
			  inbuffer_nspaces,inbuffer_maxchars,part_interval,part_modulus);

  nread = Inbuffer_fill_init(inbuffer);

  if (nread > 1) {
    multiple_sequences_p = true;
    if (offsets_access != USE_ALLOCATE || genome_access != USE_ALLOCATE) {
      fprintf(stderr,"Note: >1 sequence detected, so index files are being memory mapped.\n");
      fprintf(stderr,"  GSNAP can run slowly at first while the computer starts to accumulate\n");
      fprintf(stderr,"  pages from the hard disk into its cache.  To copy index files into RAM\n");
      fprintf(stderr,"  instead of memory mapping, use -B 3 or -B 4, if you have enough RAM.\n");
#ifdef HAVE_PTHREAD
      fprintf(stderr,"  For more speed, also try multiple threads (-t <int>), if you have multiple processors or cores.");
#endif
      fprintf(stderr,"\n");
    }
  } else {
    /* multiple_sequences_p = false; */
    /* fprintf(stderr,"Note: only 1 sequence detected.  Ignoring batch (-B) command\n"); */
    offsets_access = USE_MMAP_ONLY;
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
  }
  FREE(iitfile);


  if (snps_root == NULL) {
    genome = Genome_new(genomesubdir,fileroot,/*snps_root*/NULL,uncompressedp,genome_access);
    if (dibasep == true) {
      if ((indexdb = Indexdb_new_genome(genomesubdir,fileroot,/*idx_filesuffix*/"dibase",/*snps_root*/NULL,
					/*required_interval*/0,offsets_access,positions_access)) == NULL) {
	fprintf(stderr,"Cannot find offsets file %s.%s*offsets, needed for GSNAP color mode\n",fileroot,"dibase");
	exit(9);
      }
      indexdb2 = indexdb;
      print_ncolordiffs_p = true;

    } else if (cmetp == true) {
      if (user_cmetdir == NULL) {
	cmetdir = genomesubdir;
      } else {
	cmetdir = user_cmetdir;
      }

      indexdb = Indexdb_new_genome(cmetdir,fileroot,/*idx_filesuffix*/"metct",/*snps_root*/NULL,
				   /*required_interval*/3,offsets_access,positions_access);
      indexdb2 = Indexdb_new_genome(cmetdir,fileroot,/*idx_filesuffix*/"metga",/*snps_root*/NULL,
				    /*required_interval*/3,offsets_access,positions_access);

    } else {
      /* Normal behavior */
      if ((indexdb = Indexdb_new_genome(genomesubdir,fileroot,IDX_FILESUFFIX,/*snps_root*/NULL,
					/*required_interval*/0,offsets_access,positions_access)) == NULL) {
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

    } else if (cmetp == true) {
      if (user_cmetdir == NULL) {
	cmetdir = snpsdir;
      } else {
	cmetdir = user_cmetdir;
      }

      indexdb = Indexdb_new_genome(cmetdir,fileroot,/*idx_filesuffix*/"metct",snps_root,
				   /*required_interval*/3,offsets_access,positions_access);
      indexdb2 = Indexdb_new_genome(cmetdir,fileroot,/*idx_filesuffix*/"metga",snps_root,
				    /*required_interval*/3,offsets_access,positions_access);
      if (indexdb == NULL) {
	fprintf(stderr,"Cannot find metct index file.  Need to run cmetindex first\n");
	exit(9);
      }
      if (indexdb2 == NULL) {
	fprintf(stderr,"Cannot find metga index file.  Need to run cmetindex first\n");
	exit(9);
      }

    } else {
      indexdb = Indexdb_new_genome(snpsdir,fileroot,/*idx_filesuffix*/"ref",snps_root,
				   /*required_interval*/3,offsets_access,positions_access);
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

  Compoundpos_init_positions_free(Indexdb_positions_fileio_p(indexdb));
  Spanningelt_init_positions_free(Indexdb_positions_fileio_p(indexdb));
  Stage1_init_positions_free(Indexdb_positions_fileio_p(indexdb));
  indexdb_size_threshold = (int) (10*Indexdb_mean_size(indexdb,cmetp));
  debug(printf("Size threshold is %d\n",indexdb_size_threshold));

  if (splicesites_file != NULL) {
    if ((splicesites_iit = IIT_read(splicesites_file,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
				    /*divstring*/NULL,/*add_iit_p*/false,/*labels_read_p*/true)) != NULL) {
      fprintf(stderr,"Reading splicesite file %s locally...",splicesites_file);
    } else {
      mapdir = Datadir_find_mapdir(/*user_mapdir*/NULL,genomesubdir,dbroot);
      iitfile = (char *) CALLOC(strlen(mapdir)+strlen("/")+
				strlen(splicesites_file)+1,sizeof(char));
      sprintf(iitfile,"%s/%s",mapdir,splicesites_file);
      if ((splicesites_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
				      /*divstring*/NULL,/*add_iit_p*/true,/*labels_read_p*/true)) != NULL) {
	fprintf(stderr,"Reading splicesite file %s...",iitfile);
	FREE(iitfile);
	FREE(mapdir);
      } else {
	fprintf(stderr,"Splicesite file %s.iit not found locally or in %s.  Available files:\n",splicesites_file,mapdir);
	Datadir_list_directory(stderr,mapdir);
	fprintf(stderr,"Either install file %s or specify a full directory path\n",splicesites_file);
	exit(9);
      }
    }

    if ((donor_typeint = IIT_typeint(splicesites_iit,"donor")) < 0) {
      fprintf(stderr,"Splicesite file %s.iit does not have tag 'donor'.\n",splicesites_file);
      exit(9);
    }
    if ((acceptor_typeint = IIT_typeint(splicesites_iit,"acceptor")) < 0) {
      fprintf(stderr,"Splicesite file %s.iit does not have tag 'acceptor'.\n",splicesites_file);
      exit(9);
    }

    splicesites_divint_crosstable = IIT_divint_crosstable(chromosome_iit,splicesites_iit);
    splicesites = Splicetrie_retrieve_splicesites(&distances_observed_p,&splicetypes,&splicedists,
						  &splicestrings,&splicefrags_ref,&splicefrags_alt,
						  &nsplicesites,splicesites_iit,splicesites_divint_crosstable,
						  donor_typeint,acceptor_typeint,chromosome_iit,
						  genome,genomealt,shortsplicedist);
    if (nsplicesites == 0) {
      fprintf(stderr,"\nWarning: No splicesites observed for genome %s.  Are you sure this splicesite file was built for this genome?...",
	      dbroot);
    } else {
      Splicetrie_npartners(&nsplicepartners_skip,&nsplicepartners_obs,&nsplicepartners_max,splicesites,splicetypes,splicedists,
			   splicestrings,nsplicesites,chromosome_iit,shortsplicedist,distances_observed_p);
      if (multiple_sequences_p == true && splicetrie_precompute_p == true) {
	Splicetrie_build(&triecontents_obs,&trieoffsets_obs,&triecontents_max,&trieoffsets_max,
			 nsplicepartners_skip,nsplicepartners_obs,nsplicepartners_max,splicetypes,
			 splicestrings,nsplicesites);
	FREE(nsplicepartners_max);
	FREE(nsplicepartners_obs);
	FREE(nsplicepartners_skip);
	Splicestring_gc(splicestrings,nsplicesites);
      }
    }

    fprintf(stderr,"done\n");
  }

  FREE(genomesubdir);
  FREE(fileroot);
  FREE(dbroot);


  Substring_setup(print_ncolordiffs_p,print_nsnpdiffs_p,print_snplabels_p,
		  show_refdiff_p,snps_iit,snps_divint_crosstable,
		  splicesites_iit,splicesites_divint_crosstable,
		  donor_typeint,acceptor_typeint,trim_mismatch_score,
		  output_sam_p);
  Stage3hr_setup(invert_first_p,invert_second_p);


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

  outbuffer = Outbuffer_new(nread,sevenway_root,
#ifdef USE_OLD_MAXENT
			    genome,
#else
			    genome_blocks,
#endif
			    chromosome_iit,
			    output_sam_p,sam_headers_p,sam_read_group_id,sam_read_group_name,
			    gobywriter,nofailsp,failsonlyp,fails_as_input_p,
			    fastq_format_p,maxpaths,quiet_if_excessive_p,quality_shift,
			    invert_first_p,invert_second_p,pairmax);

  Inbuffer_set_outbuffer(inbuffer,outbuffer);

#ifdef LEAKCHECK
  nworkers = 0;
#endif

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
    pthread_attr_init(&thread_attr_detach);
    if ((ret = pthread_attr_setdetachstate(&thread_attr_detach,PTHREAD_CREATE_DETACHED)) != 0) {
      fprintf(stderr,"ERROR: pthread_attr_setdetachstate returned %d\n",ret);
      exit(1);
    }
    pthread_attr_init(&thread_attr_join);
    if ((ret = pthread_attr_setdetachstate(&thread_attr_join,PTHREAD_CREATE_JOINABLE)) != 0) {
      fprintf(stderr,"ERROR: pthread_attr_setdetachstate returned %d\n",ret);
      exit(1);
    }
    
    worker_thread_ids = (pthread_t *) CALLOC(nworkers,sizeof(pthread_t));

    Except_init_pthread();

    if (orderedp == true) {
      pthread_create(&output_thread_id,&thread_attr_join,Outbuffer_thread_ordered,
		     (void *) outbuffer);
    } else {
      pthread_create(&output_thread_id,&thread_attr_join,Outbuffer_thread_anyorder,
		     (void *) outbuffer);
    }
    for (i = 0; i < nworkers; i++) {
      pthread_create(&(worker_thread_ids[i]),&thread_attr_detach,worker_thread,
		     (void *) NULL);
    }
    
    pthread_join(output_thread_id,NULL);

    /* Do not delete global_except_key, because worker threads might still need it */
    /* Except_term_pthread(); */

    FREE(worker_thread_ids);
  }
#endif /* HAVE_PTHREAD */

  runtime = Stopwatch_stop(stopwatch);
  Stopwatch_free(&stopwatch);

  nread = Outbuffer_nread(outbuffer);
  fprintf(stderr,"Processed %d queries in %.2f seconds (%.2f queries/sec)\n",
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
    FREE(splicesites);
  }

  if (splicesites_iit != NULL) {
    FREE(splicesites_divint_crosstable);
    IIT_free(&splicesites_iit);
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
  -q, --part=INT/INT             Process only the i-th out of every n sequences\n\
                                   e.g., 0/100 or 99/100 (useful for distributing jobs\n\
                                   to a computer farm).\n\
  --input-buffer-size=INT        Size of input buffer (program reads this many sequences\n\
                                   at a time for efficiency) (default 1000)\n\
  --barcode-length=INT           Amount of barcode to remove from start of read\n\
                                   (default 0)\n\
  --pc-linefeeds                 Strip PC line feeds (ASCII 13) from input\n\
  -o, --orientation=STRING       Orientation of paired-end reads\n\
                                   Allowed values: FR (fwd-rev, or typical Illumina; default),\n\
                                   FR (rev-fwd, for circularized inserts), or FF (fwd-fwd, same strand)\n\
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
((readlength+2)/12 - 2) (\"ultrafast mismatches\").  The program will run fastest if\n\
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
                         Note: For a single sequence, all data structures use mmap\n\
                         If mmap not available and allocate not chosen, then will use fileio (slow)\n\
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
  --terminal-penalty=INT         Penalty for a terminal alignment (alignment from one end of the read\n\
                                   to the best possible position at the other end) (default 1)\n\
  -i, --indel-penalty=INT        Penalty for an indel (default 1).\n\
                                   Counts against mismatches allowed.  To find indels, make\n\
                                   indel-penalty less than or equal to max-mismatches\n\
                                   For 2-base reads, need to set indel-penalty somewhat high\n\
  -I, --indel-endlength=INT      Minimum length at end required for indel alignments (default 3)\n\
  -y, --max-middle-insertions=INT  Maximum number of middle insertions allowed (default 9)\n\
  -z, --max-middle-deletions=INT Maximum number of middle deletions allowed (default 30)\n\
  -Y, --max-end-insertions=INT   Maximum number of end insertions allowed (default 3)\n\
  -Z, --max-end-deletions=INT    Maximum number of end deletions allowed (default 6)\n\
  -M, --suboptimal-levels=INT    Report suboptimal hits beyond best hit (default 0)\n\
                                   All hits with best score plus suboptimal-levels are reported\n\
  -R, --masking=INT              Masking of frequent/repetitive oligomers to avoid spending time\n\
                                   on non-unique or repetitive reads\n\
                                   0 = no masking (will try to find non-unique or repetitive matches)\n\
                                   1 = mask frequent oligomers\n\
                                   2 = mask frequent and repetitive oligomers (fastest) (default)\n\
                                   3 = greedy frequent: mask frequent oligomers first, then\n\
                                       try no masking if alignments not found\n\
                                   4 = greedy repetitive: mask frequent and repetitive oligomers first, then\n\
                                       try no masking if alignments not found\n\
  -a, --adapter-strip=STRING     Method for removing adapters from reads.  Currently allowed values: paired\n\
  --trim-mismatch-score=INT      Score to use for mismatches when trimming at ends (default is -3;\n\
                                   to turn off trimming, specify 0)\n\
  -V, --snpsdir=STRING           Directory for SNPs index files (created using snpindex) (default is\n\
                                   location of genome index files specified using -D and -d)\n \
  -v, --use-snps=STRING          Use database containing known SNPs (in <STRING>.iit, built\n\
                                   previously using snpindex) for tolerance to SNPs\n\
  -C  --cmetdir=STRING           Directory for methylcytosine index files (created using cmetindex)\n\
                                   (default is location of genome index files specified using -D, -V, and -d)\n\
  -c, --use-cmet                 Use database for methylcytosine experiments, built\n\
                                   previously using cmetindex)\n\
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


  /* Splicing options */
  fprintf(stdout,"Splicing options for RNA-Seq\n");
  fprintf(stdout,"\
  -s, --splicesites=STRING             Look for splicing involving known splice sites\n\
                                         (in <STRING>.iit), at short or long distances\n\
  -S, --splicetrie-precompute=INT      Pre-compute splicetrie for all known splice sites\n\
                                         (0=no, 1=yes (default)).  Requires --splicesites flag\n\
                                          and multiple sequence input.\n\
  -N, --novelsplicing=INT              Look for novel splicing, not in known splice sites (if -s provided)\n\
  --novel-doublesplices                Allow GSNAP to look for two splices in a single-end involving novel\n\
                                          splice sites (default is not to allow this).  Caution: this option\n\
                                          can slow down the program considerably.  A better way to detect\n\
                                          double splices is with known splice sites, using the\n\
                                          --splicesites option.\n\
  -w, --localsplicedist=INT            Definition of local novel splicing event (default 200000)\n\
  -e, --local-splice-penalty=INT       Penalty for a local splice (default 0).\n\
                                       Counts against mismatches allowed\n\
  -E, --distant-splice-penalty=INT     Penalty for a distant splice (default 3).\n\
                                       Counts against mismatches allowed\n\
  -k, --local-splice-endlength=INT     Minimum length at end required for local spliced alignments (default 15, min is 14)\n\
  -K, --distant-splice-endlength=INT   Minimum length at end required for distant spliced alignments (default 16, min is 14)\n\
  -l, --shortend-splice-endlength=INT  Minimum length at end required for short-end spliced alignments (default 2)\n\
  --distant-splice-identity=FLOAT      Minimum identity at end required for distant spliced alignments (default 0.95)\n\
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
  -p, --pairexpect=INT           Expected paired-end length (default 200)\n\
  --pairdev=INT                  Allowable deviation from expected paired-end length, used for\n\
                                    discriminating between alternative alignments (default 50)\n\
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
");
  fprintf(stdout,"\n");

  /* SAM options */
  fprintf(stdout,"Options for SAM output\n");
  fprintf(stdout,"\
  --no-sam-headers               Do not print headers beginning with '@'\n\
  --sam-headers-batch=INT        Print headers only for this batch, as specified by -q\n\
  --read-group-id=STRING         Value to put into read-group id (RG-ID) field\n\
  --read-group-name=STRING       Value to put into read-group name (RG-SM) field\n\
");
  fprintf(stdout,"\n");


#ifdef HAVE_GOBY
  /* Goby options */
  fprintf(stdout,"Options for Goby library\n");
  fprintf(stdout,"\
  --goby-output=STRING           Basename for Goby output files\n\
  --creads-window-start=INT      Compact reads window start (default: 0=start of file)\n\
  --creads-window-end=INT        Compact reads window end (default: 0=end of file)\n\
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

