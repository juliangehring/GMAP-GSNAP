static char rcsid[] = "$Id: gsnap.c,v 1.82 2010/03/10 01:32:50 twu Exp $";
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

#include <signal.h>

#include "except.h"
#include "mem.h"
#include "bool.h"
#include "fopen.h"

#include "sequence.h"
#include "stopwatch.h"
#include "genome.h"
#include "stage3hr.h"
#include "sam.h"
#include "stage1hr.h"
#include "indexdb.h"
#include "resulthr.h"
#include "request.h"
#include "reqpost.h"
#include "blackboard.h"
#include "intlist.h"
#include "list.h"
#include "listdef.h"
#include "iit-read.h"
#include "datadir.h"
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


/************************************************************************
 *   Program options
 ************************************************************************/

/* Input options */
static char *user_genomedir = NULL;
static char *dbroot = NULL;
static char *dbversion = NULL;
static int part_modulus = 0;
static int part_interval = 1;
static bool circularp = false;

/* Compute options */
static bool novelsplicingp = false;
static bool trim_ends_p = false;
static int pairlength = 200;
static bool batch_offsets_p = true;
static bool batch_positions_p = true;
static bool batch_genome_p = false;
static int pairmax = 1000;
#ifdef HAVE_PTHREAD
static pthread_t input_thread_id, output_thread_id, *worker_thread_ids;
static int nworkers = 1;	/* (int) sysconf(_SC_NPROCESSORS_ONLN) */
#endif

static Masktype_T masktype = MASK_REPETITIVE;
static int subopt_levels = 0;

/* If negative, then hasn't been specified by user.  If between 0 and
   1, then treated as a fraction of the querylength.  Else, treated as
   an integer */
static double usermax_level_float = -1.0;

static int indel_penalty = 1000;
static bool allow_end_indels_p = true;
static int max_middle_insertions = 9;
static int max_middle_deletions = 30;
static int max_end_insertions = 3;
static int max_end_deletions = 6;
static int min_indel_end_matches = 3;
static Genomicpos_T shortsplicedist = 200000;
static int localsplicing_penalty = 2;
static int distantsplicing_penalty = 3;
static int min_localsplicing_end_matches = 15;
static int min_distantsplicing_end_matches = 16;
static double min_distantsplicing_identity = 0.95;
static int indexdb_size_threshold;


/* Splicesites IIT */
static bool knownsplicingp = false;
static char *splicesites_file = (char *) NULL;
static IIT_T splicesites_iit = NULL;
static int donor_typeint;		/* for splicesites_iit */
static int acceptor_typeint;	/* for splicesites_iit */
static int *splicesites_divint_crosstable = NULL;
static Genomicpos_T *splicesites = NULL;
static Splicetype_T *splicetypes = NULL;
static int nsplicesites = 0;

/* SNPs IIT */
static bool dibasep = false;
static bool cmetp = false;
static char *snps_root = (char *) NULL;
static IIT_T snps_iit = NULL;
static int *snps_divint_crosstable = NULL;

/* Geneprob IIT */
static char *geneprob_root = (char *) NULL;
static IIT_T geneprob_iit = NULL; /* Needs to be of cumulative IIT
				     form >(count) (genomicpos),
				     usually from pairing */


/* Output options */
static bool output_sam_p = false;
static bool exception_raise_p = true;
static bool quiet_if_excessive_p = false;
static int maxpaths = 100;
static bool orderedp = false;
static bool print_snplabels_p = false;
static bool failsonlyp = false;
static bool nofailsp = false;

/* Alignment options */
static bool uncompressedp = false;
static bool invertp = false;

/* getopt used alphabetically: ABCcDdEeFfGgIiJKkLlMmNnOPpQqRSsTtVvwYyZz2 */

static struct option long_options[] = {
  /* Input options */
  {"dir", required_argument, 0, 'D'},	/* user_genomedir */
  {"db", required_argument, 0, 'd'}, /* dbroot */
  {"genomefull", no_argument, 0, 'G'}, /* uncompressedp */
  {"part", required_argument, 0, 'q'}, /* part_modulus, part_interval */
  {"circular-input", no_argument, 0, 'c'}, /* circularp */

  /* Compute options */
#ifdef HAVE_MMAP
  {"batch", required_argument, 0, 'B'}, /* batch_offsets_p, batch_positions_p, batch_genome_p */
#endif
  {"pairmax", required_argument, 0, 'P'}, /* pairmax */
  {"pairlength", required_argument, 0, 'p'}, /* pairlength */
#ifdef HAVE_PTHREAD
  {"nthreads", required_argument, 0, 't'}, /* nworkers */
#endif
  {"trim", required_argument, 0, 'T'}, /* trim_ends_p */
  {"novelsplicing", required_argument, 0, 'N'}, /* novelsplicingp */

  {"max-mismatches", required_argument, 0, 'm'}, /* usermax_level_float */
  {"indel-penalty", required_argument, 0, 'i'}, /* indel_penalty */

  {"indel-endlength", required_argument, 0, 'I'}, /* min_indel_end_matches, allow_end_indels_p */

  {"max-middle-insertions", required_argument, 0, 'y'}, /* max_middle_insertions */
  {"max-middle-deletions", required_argument, 0, 'z'}, /* max_middle_deletions */
  {"max-end-insertions", required_argument, 0, 'Y'}, /* max_end_insertions */
  {"max-end-deletions", required_argument, 0, 'Z'}, /* max_end_deletions */
  {"suboptimal-levels", required_argument, 0, 'M'}, /* subopt_levels */
  {"masking", required_argument, 0, 'R'}, /* masktype */

  {"localsplicedist", required_argument, 0, 'w'}, /* shortsplicedist */
  {"splicesites", required_argument, 0, 's'}, /* splicesites_iit */
  {"local-splice-penalty", required_argument, 0, 'e'}, /* localsplicing_penalty */
  {"distant-splice-penalty", required_argument, 0, 'E'}, /* distantsplicing_penalty */
  {"local-splice-endlength", required_argument, 0, 'k'}, /* min_localsplicing_end_matches */
  {"distant-splice-endlength", required_argument, 0, 'K'}, /* min_distantsplicing_end_matches */
  {"distant-splice-identity", required_argument, 0, 'J'}, /* min_distantsplicing_identity */

  {"cmet", no_argument, 0, 'C'}, /* cmetp */
  {"dibase", no_argument, 0, '2'}, /* dibasep */
  {"usesnps", required_argument, 0, 'V'}, /* snps_root */
  {"geneprob", required_argument, 0, 'g'}, /* geneprob_iit */

  /* Output options */
  {"format", required_argument, 0, 'A'}, /* output_sam_p */
  {"noexceptions", no_argument, 0, '0'}, /* exception_raise_p */
  {"npaths", required_argument, 0, 'n'}, /* maxpaths */
  {"quiet-if-excessive", no_argument, 0, 'Q'}, /* quiet_if_excessive_p */
  {"ordered", no_argument, 0, 'O'}, /* orderedp */
  {"print-snps", required_argument, 0, 'S'}, /* print_snplabels_p */
  {"failsonly", no_argument, 0, 'F'}, /* failsonlyp */
  {"nofails", no_argument, 0, 'f'}, /* nofailsp */

  /* Help options */
  {"version", no_argument, 0, 'v'}, /* print_program_version */
  {"help", no_argument, 0, '?'}, /* print_program_usage */
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

static int requestid = 0;

#ifdef HAVE_PTHREAD
static void *
input_thread (void *data) {
  Blackboard_T blackboard = (Blackboard_T) data;
  FILE *input = Blackboard_input(blackboard);
  char **files = Blackboard_files(blackboard);
  int nfiles = Blackboard_nfiles(blackboard);
  Sequence_T queryseq1, queryseq2;
  Request_T request;
  int inputid = 1;		/* Initial queryseq, inputid 0, already handled by main() */
  int nextchar = Blackboard_nextchar(blackboard);

  debug(printf("input_thread: Starting\n"));
  while ((queryseq1 = Sequence_read_multifile_shortreads(&nextchar,&queryseq2,&input,&files,&nfiles,
							 circularp)) != NULL) {
    if (inputid % part_interval != part_modulus) {
      Sequence_free(&queryseq1);
      if (queryseq2 != NULL) {
	Sequence_free(&queryseq2);
      }
    } else {
      debug(printf("input_thread: Putting request id %d\n",requestid));
      request = Request_new(requestid++,queryseq1,queryseq2);
      Blackboard_put_request(blackboard,request);
    }
    inputid++;
  }

  Blackboard_set_inputdone(blackboard);
  debug(printf("input_thread: Ending\n"));
  return (void *) NULL;
}
#endif


static void
print_result_sam (Result_T result, Request_T request) {
  Sequence_T queryseq1;
  Stage3_T *stage3array;
  int npaths, pathnum;

  if (Result_resulttype(result) == SINGLEEND_READ) {
    stage3array = (Stage3_T *) Result_array(&npaths,result);
    queryseq1 = Request_queryseq1(request);
    if (npaths == 0 || (quiet_if_excessive_p && npaths > maxpaths)) {
      SAM_print_nomapping(queryseq1,/*mate*/NULL,/*acc*/Sequence_accession(queryseq1),
			  chromosome_iit,/*resulttype*/SINGLEEND_READ,
			  /*first_read_p*/true,/*nhits_mate*/0,/*queryseq_mate*/NULL);
    } else {
      qsort(stage3array,npaths,sizeof(Stage3_T),Stage3_output_cmp);
      for (pathnum = 1; pathnum <= npaths && pathnum <= maxpaths; pathnum++) {
	SAM_print(stage3array[pathnum-1],/*mate*/NULL,/*acc*/Sequence_accession(queryseq1),
		  pathnum,npaths,Stage3_score(stage3array[pathnum-1]),
		  genome,chromosome_iit,queryseq1,
		  /*queryseq2*/NULL,snps_iit,snps_divint_crosstable,
		  splicesites_iit,donor_typeint,acceptor_typeint,
		  splicesites_divint_crosstable,/*pairedlength*/0,/*resulttype*/SINGLEEND_READ,
		  /*first_read_p*/true,/*npaths_mate*/0);
      }
    }

  } else {
    /* PAIRED_CONCORDANT, PAIRED_SAMECHR, PAIRED_AS_SINGLES, or PAIRED_AS_SINGLES_UNIQUE */
    SAM_print_paired(result,genome,chromosome_iit,
		     Request_queryseq1(request),Request_queryseq2(request),
		     snps_iit,snps_divint_crosstable,
		     splicesites_iit,donor_typeint,acceptor_typeint,
		     splicesites_divint_crosstable,maxpaths,
		     quiet_if_excessive_p,circularp);
  }

  return;
}



static void
print_result_gsnap (Result_T result, Request_T request) {
  Sequence_T queryseq1;
  Stage3_T *stage3array;
  int npaths, pathnum;

  if (Result_resulttype(result) == SINGLEEND_READ) {
    stage3array = (Stage3_T *) Result_array(&npaths,result);
    if (npaths == 0 || (quiet_if_excessive_p && npaths > maxpaths)) {
      if (nofailsp == true) {
	/* Skip */
      } else if (failsonlyp == true) {
	queryseq1 = Request_queryseq1(request);
	printf(">");
	Sequence_print_header(queryseq1,/*checksump*/false);
	/* printf("\n"); -- included in header */
	Sequence_print_oneline(stdout,queryseq1);
	printf("\n");
      } else {
	queryseq1 = Request_queryseq1(request);
	printf(">");
	Sequence_print_oneline(stdout,queryseq1);
	printf("\t%d\t",npaths);
	Sequence_print_header(queryseq1,/*checksump*/false);
	/* printf("\n"); -- included in header */
	printf("\n");
      }
    } else {
      /* npaths > 0 */
      if (failsonlyp == true) {
	/* Skip */
      } else {
	queryseq1 = Request_queryseq1(request);
	printf(">");
	Sequence_print_oneline(stdout,queryseq1);
	printf("\t%d\t",npaths);
	Sequence_print_header(queryseq1,/*checksump*/false);
	/* printf("\n"); -- included in header */
	
	qsort(stage3array,npaths,sizeof(Stage3_T),Stage3_output_cmp);
	
	for (pathnum = 1; pathnum <= npaths && pathnum <= maxpaths; pathnum++) {
	  Stage3_print(stage3array[pathnum-1],Stage3_score(stage3array[pathnum-1]),
		       Stage3_geneprob(stage3array[pathnum-1]),
		       genome,chromosome_iit,queryseq1,
		       snps_iit,snps_divint_crosstable,
		       splicesites_iit,donor_typeint,acceptor_typeint,
		       splicesites_divint_crosstable,/*invertp*/false,
		       /*hit5*/(Stage3_T) NULL,/*hit3*/(Stage3_T) NULL,/*pairtype*/UNPAIRED,
		       /*pairlength*/0,/*pairscore*/0);
	}
	printf("\n");
      }
    }

  } else {
    /* PAIRED_CONCORDANT, PAIRED_SAMECHR, PAIRED_AS_SINGLES, or PAIRED_AS_SINGLES_UNIQUE */
    Stage3_print_paired(result,genome,chromosome_iit,
			Request_queryseq1(request),Request_queryseq2(request),
			snps_iit,snps_divint_crosstable,
			splicesites_iit,donor_typeint,acceptor_typeint,
			splicesites_divint_crosstable,maxpaths,
			quiet_if_excessive_p,circularp,invertp);
  }

  return;
}


static void
print_result (Result_T result, Request_T request) {
  if (output_sam_p == true) {
    print_result_sam(result,request);
  } else {
    print_result_gsnap(result,request);
  }
  return;
}



#ifdef HAVE_PTHREAD
static void *
output_thread_anyorder (void *data) {
  Blackboard_T blackboard = (Blackboard_T) data;
  Result_T result;
  Request_T request;
  int outputid = 0;

  debug(printf("output_thread: Starting\n"));
  while ((result = Blackboard_get_result(&request,blackboard)) != NULL) {
    print_result(result,request);
    Result_free(&result);
    Request_free(&request);
    outputid++;
  }

  debug(printf("output_thread: Ending\n"));
  fprintf(stderr,"Processed %d queries\n",outputid);
  return (void *) NULL;
}

static List_T
queue_insert_result (List_T list, Result_T result) {
  List_T *metap;
  int id = Result_id(result);

  metap = &list;
  while (*metap != NULL && id > Result_id((*metap)->first)) {
    metap = &(*metap)->rest;
  }
  *metap = List_push(*metap,result);

  return list;
}

static List_T
queue_insert_request (List_T list, Request_T request) {
  List_T *metap;
  int id = Request_id(request);

  metap = &list;
  while (*metap != NULL && id > Request_id((*metap)->first)) {
    metap = &(*metap)->rest;
  }
  *metap = List_push(*metap,request);

  return list;
}

static void *
output_thread_ordered (void *data) {
  Blackboard_T blackboard = (Blackboard_T) data;
  Result_T result;
  Request_T request;
  void *item;
  List_T request_queue = NULL, result_queue = NULL;
  int outputid = 0;

  debug(printf("output_thread: Starting\n"));
  while ((result = Blackboard_get_result(&request,blackboard)) != NULL) {
    if (Result_id(result) != outputid) {
      result_queue = queue_insert_result(result_queue,result);
      request_queue = queue_insert_request(request_queue,request);
    } else {
      debug(printf("a resultid %d, outputid %d\n",Result_id(result),outputid));
      print_result(result,request);
      Result_free(&result);
      Request_free(&request);
      outputid++;

      while (result_queue != NULL && Result_id(List_head(result_queue)) == outputid) {
	result_queue = List_pop(result_queue,&item);
	result = (Result_T) item;
	request_queue = List_pop(request_queue,&item);
	request = (Request_T) item;

	debug(printf("b resultid %d, outputid %d\n",Result_id(result),outputid));
	print_result(result,request);
	Result_free(&result);
	Request_free(&request);
	outputid++;
      }
    }
  }

  while (result_queue != NULL) {
    result_queue = List_pop(result_queue,&item);
    result = (Result_T) item;
    request_queue = List_pop(request_queue,&item);
    request = (Request_T) item;

    debug(printf("c resultid %d, outputid %d\n",Result_id(result),outputid));
    print_result(result,request);
    Result_free(&result);
    Request_free(&request);
  }

  debug(printf("output_thread: Ending\n"));
  fprintf(stderr,"Processed %d queries\n",outputid);
  return (void *) NULL;
}

#endif
    

static void
handle_request (Request_T request, int worker_id, Floors_T *floors_array,
		Reqpost_T reqpost, Stopwatch_T stopwatch) {
  int jobid;
  Result_T result;
  Sequence_T queryseq1, queryseq2;
  Stage3_T *stage3array, *stage3array2;
  Stage3pair_T *stage3pairarray;

  List_T hits, singlehits5, singlehits3, hitpairs;
  int npaths, npaths5, npaths3;

  jobid = Request_id(request);
  debug(printf("worker_thread %d: Got request id\n",worker_id,jobid));
  queryseq1 = Request_queryseq1(request);
  queryseq2 = Request_queryseq2(request);

  if (queryseq2 == NULL) {
    hits = Stage1_single_read(queryseq1,indexdb,indexdb2,indexdb_size_threshold,geneprob_iit,chromosome_iit,
			      genome,genomealt,floors_array,knownsplicingp,novelsplicingp,/*canonicalp*/true,
			      trim_ends_p,maxpaths,/*maxchimerapaths*/maxpaths,
			      usermax_level_float,subopt_levels,masktype,indel_penalty,
			      max_middle_insertions,max_middle_deletions,
			      allow_end_indels_p,max_end_insertions,max_end_deletions,min_indel_end_matches,
			      shortsplicedist,
			      localsplicing_penalty,distantsplicing_penalty,min_localsplicing_end_matches,
			      min_distantsplicing_end_matches,min_distantsplicing_identity,
			      splicesites,splicetypes,nsplicesites,dibasep,cmetp);
    npaths = List_length(hits);
    stage3array = (Stage3_T *) List_to_array(hits,NULL);
    List_free(&hits);
      
    result = Result_single_read_new(jobid,worker_id,(void **) stage3array,npaths);

    if (reqpost == NULL) {
      print_result(result,request);
      Result_free(&result);
    } else {
#ifdef HAVE_PTHREAD
      Reqpost_put_result(reqpost,result);
#endif
    }

  } else if ((hitpairs = Stage1_paired_read(&singlehits5,&singlehits3,
					    queryseq1,queryseq2,indexdb,indexdb2,indexdb_size_threshold,
					    geneprob_iit,chromosome_iit,genome,genomealt,floors_array,
					    knownsplicingp,novelsplicingp,/*canonicalp*/true,
					    trim_ends_p,maxpaths,/*maxchimerapaths*/maxpaths,
					    usermax_level_float,subopt_levels,masktype,indel_penalty,
					    max_middle_insertions,max_middle_deletions,
					    allow_end_indels_p,max_end_insertions,max_end_deletions,min_indel_end_matches,
					    shortsplicedist,
					    localsplicing_penalty,distantsplicing_penalty,min_localsplicing_end_matches,
					    min_distantsplicing_end_matches,min_distantsplicing_identity,
					    splicesites,splicetypes,nsplicesites,dibasep,cmetp,
					    pairmax,pairlength)) == NULL) {
    /* One read failed to match, or possible translocation */
    npaths5 = List_length(singlehits5);
    stage3array = (Stage3_T *) List_to_array(singlehits5,NULL);
    List_free(&singlehits5);
    
    npaths3 = List_length(singlehits3);
    stage3array2 = (Stage3_T *) List_to_array(singlehits3,NULL);
    List_free(&singlehits3);

    result = Result_paired_as_singles_new(jobid,worker_id,(void **) stage3array,npaths5,(void **) stage3array2,npaths3);
    if (reqpost == NULL) {
      print_result(result,request);
      Result_free(&result);
    } else {
#ifdef HAVE_PTHREAD
      Reqpost_put_result(reqpost,result);
#endif
    }

  } else {
    /* Paired hits found */
    npaths = List_length(hitpairs);
    stage3pairarray = (Stage3pair_T *) List_to_array(hitpairs,NULL);
    List_free(&hitpairs);

    result = Result_paired_read_new(jobid,worker_id,(void **) stage3pairarray,npaths);
    if (reqpost == NULL) {
      print_result(result,request);
      Result_free(&result);
    } else {
#ifdef HAVE_PTHREAD
      Reqpost_put_result(reqpost,result);
#endif
    }

  }

  return;
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
    Except_raise(&sigsegv_error,__FILE__,__LINE__);
  } else {
    fprintf(stderr,"Signal %d\n",sig);
    Except_raise(&misc_signal_error,__FILE__,__LINE__);
  }
  return;
}
#endif


static void
single_thread (FILE *input, char **files, int nfiles, int nextchar, 
	       Sequence_T queryseq1, Sequence_T queryseq2) {
  Floors_T *floors_array;
  Stopwatch_T stopwatch;
  Request_T request;
  int jobid = 0, i;

  floors_array = (Floors_T *) CALLOC(MAX_QUERYLENGTH+1,sizeof(Floors_T));
  stopwatch = (Stopwatch_T) NULL;

  while (jobid == 0 || (queryseq1 = Sequence_read_multifile_shortreads(&nextchar,&queryseq2,&input,&files,&nfiles,
								       circularp)) != NULL) {
    request = Request_new(jobid++,queryseq1,queryseq2);
TRY
    handle_request(request,/*worker_id*/0,floors_array,(Reqpost_T) NULL,stopwatch);
ELSE
    if (Sequence_accession(queryseq1) == NULL) {
      fprintf(stderr,"Problem with unnamed sequence (%d bp):\n",Sequence_fulllength_given(queryseq1));
    } else {
      fprintf(stderr,"Problem with sequence %s (%d bp):\n",
	      Sequence_accession(queryseq1),Sequence_fulllength_given(queryseq1));
    }
    Sequence_print_oneline(stderr,queryseq1);
    fprintf(stderr,"\n");

    fprintf(stderr,"To obtain a core dump, re-run program on problem sequence with the -0 [zero] flag\n");
    fprintf(stderr,"Exiting...\n");
    exit(9);
  RERAISE;
END_TRY;
    Request_free(&request);
   /* Don't free queryseq; done by Request_free */
  }

  Stopwatch_free(&stopwatch);
  for (i = 0; i <= MAX_QUERYLENGTH; i++) {
    if (floors_array[i] != NULL) {
      Floors_free(&(floors_array[i]));
    }
  }
  FREE(floors_array);

  return;
}


#ifdef HAVE_PTHREAD
Blackboard_T blackboard_global;	/* Needed only for exception handling */

static void *
worker_thread (void *data) {
  Floors_T *floors_array;
  Stopwatch_T stopwatch;
  Request_T request;
  Sequence_T queryseq1;

  Reqpost_T reqpost = (Reqpost_T) data;
  int worker_id = Reqpost_id(reqpost), id, i;

  /* Thread-specific data and storage */
  floors_array = (Floors_T *) CALLOC(MAX_QUERYLENGTH+1,sizeof(Floors_T));
  stopwatch = (Stopwatch_T) NULL;
  Except_stack_create();

  while ((request = Reqpost_get_request(reqpost))) {
TRY
    handle_request(request,worker_id,floors_array,reqpost,stopwatch);
ELSE
    /* Kill output thread first, so requests (and queryseqs) can't disappear */
    fprintf(stderr,"Killing output thread...\n");
    pthread_kill(output_thread_id,SIGUSR1);

    fprintf(stderr,"Thread assignments: ");
    for (id = 0; id < nworkers; id++) {
      if (id > 0) {
	fprintf(stderr,", ");
      }
      if (id == worker_id) {
        fprintf(stderr,"**");
      }

      /* Get requests for all threads */
      request = Reqpost_get_request(Blackboard_get_reqpost(blackboard_global,id));
      if (request == NULL) {
        fprintf(stderr,"NULL");
      } else {
  	queryseq1 = Request_queryseq1(request);
	if (queryseq1 == NULL) {
	  fprintf(stderr,"NULL");
	} else if (Sequence_accession(queryseq1) == NULL) {
	  fprintf(stderr,"unnamed (%d bp)",Sequence_fulllength_given(queryseq1));
	} else {
	  fprintf(stderr,"%s (%d bp)",Sequence_accession(queryseq1),Sequence_fulllength_given(queryseq1));
	}
        if (id == worker_id) {
	  fprintf(stderr,"\n");
	  Sequence_print_oneline(stderr,queryseq1);
	  fprintf(stderr,"\n");
        }
      }
    }
    fprintf(stderr,"\n");
    fprintf(stderr,"** indicates thread receiving the exception\n");
    fprintf(stderr,"To obtain a core dump, re-run program on problem sequence with the -0 [zero] flag\n");

#if 0
    for (id = 0; id < nworkers; id++) {
      if (id != worker_id) {
	pthread_kill(worker_thread_ids[id],SIGUSR1);
      }
    }
    pthread_kill(input_thread_id,SIGUSR1);
#endif

    fprintf(stderr,"Exiting...\n");
    exit(9);

  RERAISE;
END_TRY;
    /* Don't free request; done by output thread */
    /* Don't free queryseq; done by Request_free */
  }

  Except_stack_destroy();

  Stopwatch_free(&stopwatch);
  for (i = 0; i <= MAX_QUERYLENGTH; i++) {
    if (floors_array[i] != NULL) {
      Floors_free(&(floors_array[i]));
    }
  }
  FREE(floors_array);

  Reqpost_free(&reqpost);

  return (void *) NULL;
}
#endif


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
  Sequence_T queryseq1, queryseq2;
  char *genomesubdir = NULL, *mapdir = NULL, *iitfile = NULL, *fileroot = NULL;
  FILE *input;

  bool multiple_sequences_p = false;
  char **files;
  int nfiles, nextchar = '\0';

#ifdef HAVE_PTHREAD
  int ret, i;
  pthread_attr_t thread_attr_detach, thread_attr_join;
  Blackboard_T blackboard;
  Request_T request;
#endif

  int opt;
  extern int optind;
  extern char *optarg;
  int long_option_index = 0;

#ifdef HAVE_SIGACTION
  struct sigaction signal_action;
#endif

  while ((opt = getopt_long(argc,argv,
			    "q:D:d:GT:N:R:M:m:i:I:y:Y:z:Z:w:E:e:J:K:k:s:2V:g:B:P:p:t:A:0n:QCOsSFfcv?",
			    long_options, &long_option_index)) != -1) {
    switch (opt) {
    case 'q': parse_part(&part_modulus,&part_interval,optarg); break;
    case 'D': user_genomedir = optarg; break;
    case 'd': 
      dbroot = (char *) CALLOC(strlen(optarg)+1,sizeof(char));
      strcpy(dbroot,optarg);
      break;
    case 'G': uncompressedp = true; break;

    case 'T':
      if (!strcmp(optarg,"1")) {
	trim_ends_p = true;
      } else if (!strcmp(optarg,"0")) {
	trim_ends_p = false;
      } else {
	fprintf(stderr,"Trim (-T flag) must be 0 or 1\n");
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
    case 'J': min_distantsplicing_identity = atof(optarg); break;
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

    case 's': splicesites_file = optarg; knownsplicingp = true; break;
    case '2': dibasep = true; break;
    case 'C': cmetp = true; break;
    case 'V': snps_root = optarg; break;
    case 'g': geneprob_root = optarg; break;

    case 'B':
      if (!strcmp(optarg,"0")) {
	batch_offsets_p = true; batch_positions_p = false; batch_genome_p = false;
      } else if (!strcmp(optarg,"1")) {
	batch_offsets_p = true; batch_positions_p = true; batch_genome_p = false;
      } else if (!strcmp(optarg,"2")) {
	batch_offsets_p = true; batch_positions_p = true; batch_genome_p = true;
      } else {
	fprintf(stderr,"Batch mode %s not recognized.\n",optarg);
	fprintf(stderr,"Mode 0 means no pre-loading, mode 1 preloads indices only; mode 2 preloads both indices and genome\n");
	exit(9);
      }
      break;
    case 'P': pairmax = atoi(optarg); break;
    case 'p': pairlength = atoi(optarg); break;
#ifdef HAVE_PTHREAD
    case 't': nworkers = atoi(optarg); break;
#endif

    case 'A':
      if (!strcmp(optarg,"sam")) {
	output_sam_p = true;
      } else {
	fprintf(stderr,"Output format %s not recognized\n",optarg);
      }
      break;

    case '0': exception_raise_p = false; break; /* Allows signals to pass through */
    case 'n': maxpaths = atoi(optarg); break;
    case 'Q': quiet_if_excessive_p = true; break;

    case 'O': orderedp = true; break;
    case 'S': 
      if (!strcmp(optarg,"0")) {
	print_snplabels_p = false;
      } else if (!strcmp(optarg,"1")) {
	print_snplabels_p = true;
      } else {
	fprintf(stderr,"SNP printing mode %s not recognized.\n",optarg);
	exit(9);
      }
      break;

    case 'F': failsonlyp = true; break;
    case 'f': nofailsp = true; break;
    case 'c': circularp = true; break;

    case 'v': print_program_version(); exit(0);
    case '?': print_program_usage(); exit(0);
    default: exit(9);
    }
  }
  argc -= optind;
  argv += optind;
      
  if (dbroot == NULL) {
    fprintf(stderr,"Need to specify the -d flag\n");
    print_program_usage();
    exit(9);
  }

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

  /* Open query stream */
  if (argc == 0) {
    input = stdin;
    files = (char **) NULL;
    nfiles = 0;
  } else {
    input = NULL;
    files = argv;
    nfiles = argc;
  }

  /* Read first query sequence */
  if ((queryseq1 = Sequence_read_multifile_shortreads(&nextchar,&queryseq2,&input,&files,&nfiles,
						      circularp)) == NULL) {
    fprintf(stderr,"No input detected\n");
    exit(9);
  } else if (nextchar == '>' || nfiles >= 1) {
    multiple_sequences_p = true;
#ifdef HAVE_MMAP
    if (batch_offsets_p == false && batch_positions_p == false) {
      fprintf(stderr,"Note: >1 sequence detected.  To speed up very large runs, try batch mode (-B 1 or -B 2)\n");
      fprintf(stderr,"  if you have sufficient RAM.");
#ifdef HAVE_PTHREAD
      fprintf(stderr,"  For more speed, also try multiple threads (-t <int>),\n");
      fprintf(stderr,"  if you have multiple processors.");
#endif
      fprintf(stderr,"\n");
    }
#endif
    batch_offsets_p = true;	/* Results in pre-reading of offsets file */
  } else {
    /* multiple_sequences_p = false; */
    if (batch_offsets_p == true || batch_positions_p == true || batch_genome_p == true) {
      /* fprintf(stderr,"Note: only 1 sequence detected.  Ignoring batch (-B) command\n"); */
      batch_offsets_p = false;
      batch_positions_p = false;
      batch_genome_p = false;
    }
  }

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


  if (geneprob_root != NULL) {
    iitfile = (char *) CALLOC(strlen(geneprob_root)+1,sizeof(char));
    sprintf(iitfile,"%s",geneprob_root);
    fprintf(stderr,"Reading geneprob file %s...",geneprob_root);
    if ((geneprob_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
				 /*divstring*/NULL,/*add_iit_p*/true,/*labels_read_p*/true)) == NULL) {
      fprintf(stderr,"Geneprob file %s.iit not found\n",geneprob_root);
      exit(9);
    }

    fprintf(stderr,"done\n");
    FREE(iitfile);
  }

  if (snps_root == NULL) {
    genome = Genome_new(genomesubdir,fileroot,/*snps_root*/NULL,uncompressedp,batch_genome_p);
    if (dibasep == true) {
      if ((indexdb = Indexdb_new_genome(genomesubdir,fileroot,/*idx_filesuffix*/"dibase",/*snps_root*/NULL,
					/*required_interval*/0,batch_offsets_p,batch_positions_p)) == NULL) {
	fprintf(stderr,"Cannot find offsets file %s.%s*offsets, needed for GSNAP color mode\n",fileroot,"dibase");
	exit(9);
      }
      indexdb2 = indexdb;
      Stage3hr_print_ncolordiffs();

    } else if (cmetp == true) {
      indexdb = Indexdb_new_genome(genomesubdir,fileroot,/*idx_filesuffix*/"metct",/*snps_root*/NULL,
				   /*required_interval*/3,batch_offsets_p,batch_positions_p);
      indexdb2 = Indexdb_new_genome(genomesubdir,fileroot,/*idx_filesuffix*/"metga",/*snps_root*/NULL,
				    /*required_interval*/3,batch_offsets_p,batch_positions_p);

    } else {
      /* Normal behavior */
      if ((indexdb = Indexdb_new_genome(genomesubdir,fileroot,IDX_FILESUFFIX,/*snps_root*/NULL,
					/*required_interval*/0,batch_offsets_p,batch_positions_p)) == NULL) {
	fprintf(stderr,"Cannot find offsets file %s.%s*offsets, needed for GSNAP\n",fileroot,IDX_FILESUFFIX);
	exit(9);
      }
      indexdb2 = indexdb;
    }

  } else {
    /* SNPs */
    genome = Genome_new(genomesubdir,fileroot,/*snps_root*/NULL,uncompressedp,batch_genome_p);
    genomealt = Genome_new(genomesubdir,fileroot,snps_root,uncompressedp,batch_genome_p);

    if (dibasep == true) {
      fprintf(stderr,"Currently cannot combine SNPs with 2-base encoding\n");
      exit(9);
      Stage3hr_print_ncolordiffs();

    } else if (cmetp == true) {
      indexdb = Indexdb_new_genome(genomesubdir,fileroot,/*idx_filesuffix*/"metct",snps_root,
				   /*required_interval*/3,batch_offsets_p,batch_positions_p);
      indexdb2 = Indexdb_new_genome(genomesubdir,fileroot,/*idx_filesuffix*/"metga",snps_root,
				    /*required_interval*/3,batch_offsets_p,batch_positions_p);
    } else {
      indexdb = Indexdb_new_genome(genomesubdir,fileroot,/*idx_filesuffix*/"ref",snps_root,
				   /*required_interval*/3,batch_offsets_p,batch_positions_p);
      indexdb2 = indexdb;
    }
    Stage3hr_print_nsnpdiffs(print_snplabels_p);

    mapdir = Datadir_find_mapdir(/*user_mapdir*/NULL,genomesubdir,dbroot);

    iitfile = (char *) CALLOC(strlen(mapdir)+strlen("/")+
			      strlen(snps_root)+1,sizeof(char));
    sprintf(iitfile,"%s/%s",mapdir,snps_root);
    fprintf(stderr,"Reading SNPs file %s...",snps_root);
    if ((snps_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
			     /*divstring*/NULL,/*add_iit_p*/true,/*labels_read_p*/true)) == NULL) {
      fprintf(stderr,"SNPs file %s.iit not found in %s.  Available files:\n",snps_root,mapdir);
      Datadir_list_directory(stderr,mapdir);
      fprintf(stderr,"Either install file %s or specify a full directory path\n",snps_root);
      fprintf(stderr,"using the -D flag to gsnap.\n");
      exit(9);
    }

    snps_divint_crosstable = IIT_divint_crosstable(chromosome_iit,snps_iit);

    fprintf(stderr,"done\n");
    FREE(iitfile);
    FREE(mapdir);
  }

  indexdb_size_threshold = (int) (10*Indexdb_mean_size(indexdb,cmetp));
  debug(printf("Size threshold is %d\n",indexdb_size_threshold));

  if (splicesites_file != NULL) {
    mapdir = Datadir_find_mapdir(/*user_mapdir*/NULL,genomesubdir,dbroot);
    iitfile = (char *) CALLOC(strlen(mapdir)+strlen("/")+
			      strlen(splicesites_file)+1,sizeof(char));
    sprintf(iitfile,"%s/%s",mapdir,splicesites_file);
    fprintf(stderr,"Reading splicesite file %s...",splicesites_file);
    if ((splicesites_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
				    /*divstring*/NULL,/*add_iit_p*/true,/*labels_read_p*/true)) == NULL) {
      fprintf(stderr,"Splicesite file %s.iit not found in %s.  Available files:\n",splicesites_file,mapdir);
      Datadir_list_directory(stderr,mapdir);
      fprintf(stderr,"Either install file %s or specify a full directory path\n",splicesites_file);
      fprintf(stderr,"using the -D flag to gsnap.\n");
      exit(9);
    }
    if ((donor_typeint = IIT_typeint(splicesites_iit,"donor")) < 0) {
      fprintf(stderr,"Splicesite file %s.iit does not have tag 'donor'\n",splicesites_file);
      exit(9);
    }
    if ((acceptor_typeint = IIT_typeint(splicesites_iit,"acceptor")) < 0) {
      fprintf(stderr,"Splicesite file %s.iit does not have tag 'acceptor'\n",splicesites_file);
      exit(9);
    }

    splicesites_divint_crosstable = IIT_divint_crosstable(chromosome_iit,splicesites_iit);
    splicesites = Stage1_retrieve_splicesites(&splicetypes,&nsplicesites,splicesites_iit,
					      splicesites_divint_crosstable,donor_typeint,acceptor_typeint,
					      chromosome_iit);

    fprintf(stderr,"done\n");
    FREE(iitfile);
    FREE(mapdir);
  }

  FREE(genomesubdir);
  FREE(fileroot);
  FREE(dbroot);

  if (output_sam_p == true) {
    IIT_dump_sam(chromosome_iit);
  }

#ifndef HAVE_PTHREAD
  single_thread(input,files,nfiles,nextchar,queryseq1,queryseq2);
#else
  if (multiple_sequences_p == false) {
    single_thread(input,files,nfiles,nextchar,queryseq1,queryseq2);
  } else if (nworkers == 0) {
    single_thread(input,files,nfiles,nextchar,queryseq1,queryseq2);
  } else {
    /* Make blackboard and threads */
    blackboard = Blackboard_new(input,files,nfiles,nextchar,/*usersegment*/NULL,nworkers);
    blackboard_global = blackboard; /* Needed only for exception handling */
    debug(printf("input_thread: Putting request id 0\n"));
    if (part_modulus != 0) {
      Sequence_free(&queryseq1);
      if (queryseq2 != NULL) {
	Sequence_free(&queryseq2);
      }
    } else {
      request = Request_new(requestid++,queryseq1,queryseq2);
      Blackboard_put_request(blackboard,request);
    }
    
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
      pthread_create(&output_thread_id,&thread_attr_join,output_thread_ordered,
		     (void *) blackboard);
    } else {
      pthread_create(&output_thread_id,&thread_attr_join,output_thread_anyorder,
		     (void *) blackboard);
    }
    for (i = 0; i < nworkers; i++) {
      pthread_create(&(worker_thread_ids[i]),&thread_attr_detach,worker_thread,
		     (void *) Blackboard_get_reqpost(blackboard,i));
    }
    pthread_create(&input_thread_id,&thread_attr_detach,input_thread,(void *) blackboard);
    
    pthread_join(output_thread_id,NULL);

    /* Do not delete global_except_key, because worker threads might still need it */
    /* Except_term_pthread(); */

    FREE(worker_thread_ids);
    Blackboard_free(&blackboard);
  }
#endif /* HAVE_PTHREAD */

  if (argc > 0) {
    fclose(input);
  }

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
  }
  if (genome != NULL) {
    Genome_free(&genome);
  }

  if (nsplicesites > 0) {
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
Input options (must include -d)\n\
  -D, --dir=directory            Genome directory\n\
  -d, --db=STRING                Genome database\n\
  -q, --part=INT/INT             Process only the i-th out of every n sequences\n\
                                   e.g., 0/100 or 99/100\n\
  -c, --circular-input           Circular-end data (paired reads are on same strand)\n\
\n\
");

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
  -B, --batch=INT                Batch mode (0 = no pre-loading, 1 = pre-load only indices;\n\
                                   2 (default) = pre-load both indices and genome)\n\
");
#endif
    fprintf(stdout,"\
  -m, --max-mismatches=FLOAT     Maximum number of mismatches allowed (if not specified, then\n\
                                   defaults to the ultrafast level of ((readlength+2)/12 - 2))\n\
                                   If specified between 0.0 and 1.0, then treated as a fraction\n\
                                   of each read length.  Otherwise, treated as an integral number\n\
                                   of mismatches (including indel and splicing penalties)\n\
  -i, --indel-penalty=INT        Penalty for an indel (default 1000, essentially turning it off).\n\
                                   Counts against mismatches allowed.  To find indels, make\n\
                                   indel-penalty less than or equal to max-mismatches\n\
                                   For 2-base reads, need to set indel-penalty somewhat high\n\
  -I, --indel-endlength=INT      Minimum length at end required for indel alignments (default 3)\n\
  -y, --max-middle-insertions=INT  Maximum number of middle insertions allowed (default 9)\n\
  -z, --max-middle-deletions=INT Maximum number of middle deletions allowed (default 30)\n\
  -Y, --max-end-insertions=INT   Maximum number of end insertions allowed (default 3)\n\
  -Y, --max-end-deletions=INT    Maximum number of end deletions allowed (default 6)\n\
  -M, --suboptimal-score=INT     Report suboptimal hits beyond best hit (default 0)\n\
                                   All hits with best score plus suboptimal-score are reported\n\
  -R, --masking=INT              Masking of frequent/repetitive oligomers to avoid spending time\n\
                                   on non-unique or repetitive reads\n\
                                   0 = no masking (will try to find non-unique or repetitive matches)\n\
                                   1 = mask frequent oligomers\n\
                                   2 = mask frequent and repetitive oligomers (fastest) (default)\n\
                                   3 = greedy frequent: mask frequent oligomers first, then\n\
                                       try no masking if alignments not found\n\
                                   4 = greedy repetitive: mask frequent and repetitive oligomers first, then\n\
                                       try no masking if alignments not found\n\
  -T, --trim=INT                 Trim mismatches at ends (0 = no (default), 1 = yes)\n\
  -2, --dibase                   Input is 2-base encoded (e.g., SOLiD), with database built\n\
                                   previously using dibaseindex)\n\
  -C, --cmet                     Use database for methylcytosine experiments, built\n\
                                   previously using cmetindex)\n\
  -V, --usesnps=STRING           Use database containing known SNPs (in <STRING>.iit, built\n\
                                   previously using snpindex) for tolerance to SNPs\n\
  -g, --geneprob=STRING          Use IIT file containing geneprob (in <STRING>.iit, of cumulative\n\
                                   format  >(count) (genomicpos)  to resolve ties\n\
");
#ifdef HAVE_PTHREAD
    fprintf(stdout,"\
  -t, --nthreads=INT             Number of worker threads\n\
");
#endif
    fprintf(stdout,"\n");


    fprintf(stdout,"Splicing options for RNA-Seq\n");
    fprintf(stdout,"\
  -s, --splicesites=STRING       Look for splicing involving known splice sites\n\
                                   (in <STRING>.iit), at short or long distances\n\
  -N, --novelsplicing=INT        Look for novel splicing, not in known splice sites (if -s provided)\n\
                                  within shortsplicedist (-w flag) or with novelspliceprob (-x flag)\n\
  -w, --localsplicedist=INT      Definition of local novel splicing event (default 200000)\n\
  -e, --local-splice-penalty=INT       Penalty for a local splice (default 2).\n\
                                       Counts against mismatches allowed\n\
  -E, --distant-splice-penalty=INT     Penalty for a distant splice (default 3).\n\
                                       Counts against mismatches allowed\n\
  -k, --local-splice-endlength=INT     Minimum length at end required for local spliced alignments (default 15, min is 14)\n\
  -K, --distant-splice-endlength=INT   Minimum length at end required for distant spliced alignments (default 16, min is 14)\n\
  -J, --distant-splice-identity=FLOAT  Minimum identity at end required for distant spliced alignments (default 0.95)\n\
");
    fprintf(stdout,"\n");


    fprintf(stdout,"Options for paired-end reads\n");
    fprintf(stdout,"\
  -P, --pairmax=INT              Max total genomic length for paired reads\n\
                                   (default 1000).  Should increase for RNA-Seq reads.\n\
  -p, --pairlength=INT           Expected paired-end length (default 200)\n\
");
    fprintf(stdout,"\n");


    fprintf(stdout,"\
Output options\n\
  -n, --npaths=INT               Maximum number of paths to print (default 100).\n\
  -Q, --quiet-if-excessive       If more than maximum number of paths are found,\n\
                                    then nothing is printed.\n\
  -O, --ordered                  Print output in same order as input (relevant\n\
                                   only if there is more than one worker thread)\n\
  -S, --print-snps=INT           Print detailed information about SNPs in reads (works only if -V also selected)\n\
                                   (0=no (default), 1=positions and labels)\n\
  -F, --failsonly                Print only failed alignments, those with no results\n\
  -f, --nofails                  Exclude printing of failed alignments\n\
  -A, --format=STRING            Another format type, other than default.\n\
                                   Currently implemented: sam\n\
");
    fprintf(stdout,"\n");

    fprintf(stdout,"\
Help options\n\
  -v, --version                  Show version\n\
  -?, --help                     Show this help message\n\
");
  return;
}

