static char rcsid[] = "$Id: gmap.c 53583 2011-12-02 18:23:41Z twu $";
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

#ifdef HAVE_PTHREAD
#include <pthread.h>
#endif

#include <signal.h>

#include "except.h"
#include "mem.h"
#include "bool.h"
#include "fopen.h"
#include "access.h"

#include "sequence.h"
#include "oligoindex.h"
#include "match.h"
#include "matchpool.h"
#include "pairpool.h"
#include "diagpool.h"
#include "stopwatch.h"
#include "genome.h"
#include "genome_hr.h"		/* For Genome_hr_setup */
#include "maxent_hr.h"		/* For Maxent_hr_setup */
#include "stage1.h"
#include "gregion.h"
#include "oligoindex_hr.h"	/* For Oligoindex_hr_setup */
#include "stage2.h"
#include "splicetrie.h"
#include "dynprog.h"
#include "stage3.h"
#include "chimera.h"
#ifdef PMAP
#include "backtranslation.h"
#endif
#include "indexdb.h"
#include "result.h"
#include "request.h"
#include "intlist.h"
#include "list.h"
#include "listdef.h"
#include "iit-read.h"
#include "chrsubset.h"
#include "datadir.h"
#include "inbuffer.h"
#include "outbuffer.h"
#include "getopt.h"


#define POSSIBLE_OLIGOS 65536	/* 4^8 */
#define MAX_OLIGODEPTH 3.0
#define MAX_BADOLIGOS 0.30	/* Setting to 1.0 effectively turns this check off */
#define MAX_REPOLIGOS 0.40	/* Setting to 1.0 effectively turns this check off */

#define CHIMERA_IDENTITY 0.98
#define CHIMERA_PVALUE 0.01
#define CHIMERA_FVALUE 6.634897	/* qnorm(CHIMERA_PVALUE/2)^2 */
#define CHIMERA_HANDICAP 10	/* points, for minor alignment differences based on different defect rates */
#ifdef PMAP
#define CHIMERA_SLOP 2
#else
#define CHIMERA_SLOP 6
#endif

#define MIN_MATCHES 20


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Chimera detection */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif


/************************************************************************
 *   Global variables
 ************************************************************************/

static IIT_T chromosome_iit = NULL;
static int nchrs;
static Genomicpos_T genome_totallength = 0;
static Chrsubset_T chrsubset = NULL;
static IIT_T contig_iit = NULL;
static Genome_T genome = NULL;
static Genome_T genomealt = NULL;
static UINT4 *genome_blocks = NULL;
static UINT4 *snp_blocks = NULL;

#ifdef PMAP
static Indexdb_T indexdb_fwd = NULL;
static Indexdb_T indexdb_rev = NULL;
static int index1part_aa = 6;
#else
static Indexdb_T indexdb = NULL;
static int index1part;
#endif
static int required_index1part = 0;
static int indexdb_size_threshold = 0;

static IIT_T altstrain_iit = NULL;

static char *snps_root = (char *) NULL;
static IIT_T map_iit = NULL;
static int *map_divint_crosstable = NULL;

#ifdef PMAP
#if 0
static int minindexsize = 3;	/* In stage 2; in aa */
static int maxindexsize = 6;	/* In stage 2; in aa */
#endif
static int maxpeelback = 12;	/* Needs to be at least indexsize
				   because stage 2 jumps by indexsize.
				   Also should exceed length of
				   repeated nucleotides (e.g., a
				   string of consecutive T's) */
#else
/* Making minindexsize too small can lead to spurious exons in stage 2 */
/* FOOBAR */
#if 0
static int minindexsize = 8;	/* In stage 2; in nt.  Used if sampling required in stage 1. */
static int maxindexsize = 8;	/* In stage 2; in nt */
#endif
static int maxpeelback = 11;	/* Needs to be at least indexsize
				   because stage 2 jumps by indexsize.
				   Also should exceed length of
				   repeated nucleotides (e.g., a
				   string of consecutive T's) */
#endif
static int maxpeelback_distalmedial = 100; /* Needs to be longer to fix bad end exons */

/* static int stuttercycles = 2; */
static int stutterhits = 3;
static int sufflookback = 60;
static int nsufflookback = 5;

#if 0
static int maxoligohits = 400; /* Must be smaller than ALLOC in oligoindex.c */
#endif
static int nullgap = 600;
static int extramaterial_end = 10;
static int extramaterial_paired = 8; /* Should be at least indexsize in nt */
static int extraband_single = 3; /* This is in addition to length2 -
				    length1.  If onesidegap is true in
				    dynprog.c, then this is equivalent
				    to extraband_single of 0.  Needs
				    to be > 0 to handle default
				    close_indels_mode. */
static int extraband_end = 3; /* Was 6.  Shouldn't differ from 0, since onesidegapp is true?
				 This is only on both sides of main diagonal */
static int extraband_paired = 7; /* This is in addition to length2 - length1 */
static int minendexon = 9;

static Stopwatch_T stopwatch = NULL;


/************************************************************************
 *   Program options
 ************************************************************************/

/* Input options */
static char *user_genomedir = NULL;
static char *dbroot = NULL;
static char *dbversion = NULL;
static char *user_genomicseg = NULL;
static bool user_pairalign_p = false;
static char *user_cmdline = NULL;
static Sequence_T usersegment = NULL;
static char *user_chrsubsetfile = NULL;
static int part_modulus = 0;
static int part_interval = 1;

/* Compute options */
static int min_matches;
static Access_mode_T offsetscomp_access = USE_ALLOCATE;
static bool expand_offsets_p = false;

#ifdef HAVE_MMAP
static Access_mode_T positions_access = USE_MMAP_PRELOAD;
static Access_mode_T genome_access = USE_MMAP_PRELOAD;
#else
static Access_mode_T positions_access = USE_ALLOCATE;
static Access_mode_T genome_access = USE_ALLOCATE;
#endif

static int min_intronlength = 9;
static int max_deletionlength = 50;
static int maxtotallen_bound = 2400000;
static int maxintronlen_bound = 1000000;
static int chimera_margin = 40;	/* Useful for finding readthroughs */
static bool maponlyp = false;
#ifdef PMAP
static bool userstage1p = false; /* Apply stage 1 for user-provided genomic segments.  Must be false. */
#else
static bool userstage1p = false; /* Apply stage 1 for user-provided genomic segments */
#endif
static int index1interval = 3; /* Stage 1 interval if user provides a genomic segment */
static char *referencefile = NULL;

#if 0
#ifndef PMAP
static bool literalrefp = false;
#endif
#endif

static bool altstrainp = false;
#ifdef HAVE_PTHREAD
static pthread_t output_thread_id, *worker_thread_ids;
static int nworkers = 1;	/* (int) sysconf(_SC_NPROCESSORS_ONLN) */
#endif
#ifndef PMAP
static bool prune_poor_p = false;
static bool prune_repetitive_p = false;
#endif
static int canonical_mode = 1;
static bool use_shifted_canonical_p = false; /* Use this for cross-species */
static char *user_chrsubsetname = NULL;
static int close_indels_mode = +1;
static double microexon_spliceprob = 0.95;
static int suboptimal_score_start = -1; /* Determined by simulations to have minimal effect */
static int suboptimal_score_end = 3; /* Determined by simulations to have diminishing returns above 3 */

static int trim_mismatch_score = -3;
static int trim_indel_score = -4;


/* Output options */
static unsigned int output_buffer_size = 1000;
static Printtype_T printtype = SIMPLE;
static bool exception_raise_p = true;
static bool debug_graphic_p = false;
static bool stage1debug = false;
static bool diag_debug = false;
static Stage3debug_T stage3debug = NO_STAGE3DEBUG;
static bool diagnosticp = false;
static bool checkp = false;
static int maxpaths = 5;	/* 0 means 1 if nonchimeric, 2 if chimeric */
static bool quiet_if_excessive_p = false;
static int suboptimal_score = 1000000;


/* SAM */
#ifndef PMAP
static bool sam_paired_p = false;
static bool user_quality_shift = false;
static int quality_shift = 0;
static bool sam_headers_p = true;
static char *sam_read_group_id = NULL;
static char *sam_read_group_name = NULL;
static char *sam_read_group_library = NULL;
static char *sam_read_group_platform = NULL;
#endif

static bool orderedp = false;
static bool failsonlyp = false;
static bool nofailsp = false;
static bool fails_as_input_p = false;
static bool checksump = false;
static int chimera_overlap = 0;

/* Map file options */
static char *user_mapdir = NULL;
static char *map_iitfile = NULL;
static bool map_exons_p = false;
static bool map_bothstrands_p = false;
static bool print_comment_p = false;
static int nflanking = 0;

/* Alignment options */
static bool fulllengthp = false;
static int cds_startpos = -1;
static bool truncatep = false;
static int sense_try = 0;		/* both */
static int sense_filter = 0;		/* both */
static bool strictp = true;
static int proteinmode = 1;
static bool uncompressedp = false;
static bool nointronlenp = false;
static int invertmode = 0;
static int ngap = 3;
static int wraplength = 50;


/* Splicing IIT */
static bool novelsplicingp = true; /* Can be disabled with --nosplicing flag */
static bool knownsplicingp = false;
static bool distances_observed_p = false;
static Genomicpos_T shortsplicedist = 200000;
static int min_extra_end;		      /* If knownsplicing, then equals shortsplicedist */
static char *user_splicingdir = (char *) NULL;
static char *splicing_file = (char *) NULL;
static IIT_T splicing_iit = NULL;
static bool amb_closest_p = false;

static int donor_typeint = -1;		/* for splicing_iit */
static int acceptor_typeint = -1;	/* for splicing_iit */

static int *splicing_divint_crosstable = NULL;
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


/* Input/output */
static char *sevenway_root = NULL;
static Inbuffer_T inbuffer = NULL;
static Outbuffer_T outbuffer = NULL;
static unsigned int inbuffer_nspaces = 1000;
static unsigned int inbuffer_maxchars = -1U; /* Currently not used by Inbuffer_T */


#ifdef PMAP
/* Used alphabetically: 01235789ABbCcDdEefGgHIiKkLlMmNnOoPQRSstuVvwXxYZ */
#else
/* Used alphabetically: 01235789AaBbCcDdEeFfGgHIijKkLlMmNnOoPpQRSsTtuVvwXxYZ */
#endif

static struct option long_options[] = {
  /* Input options */
  {"dir", required_argument, 0, 'D'},	/* user_genomedir */
  {"db", required_argument, 0, 'd'}, /* dbroot */
  {"kmer", required_argument, 0, 'k'}, /* required_index1part, index1part */
  {"genomefull", no_argument, 0, 'G'}, /* uncompressedp */
  {"gseg", required_argument, 0, 'g'}, /* user_genomicseg */
  {"pairalign", no_argument, 0, '2'}, /* user_pairalign_p */
  {"cmdline", required_argument, 0, 0}, /* user_cmdline */
  {"part", required_argument, 0, 'q'}, /* part_modulus, part_interval */
  {"input-buffer-size", required_argument, 0, 0}, /* inbuffer_nspaces */

  /* Compute options */
#ifdef HAVE_MMAP
  {"batch", required_argument, 0, 'B'}, /* offsetscomp_access, positions_access, genome_access */
#endif
  {"min-intronlength", required_argument, 0, 0}, /* min_intronlength */
  {"intronlength", required_argument, 0, 'K'}, /* maxintronlen_bound */
  {"totallength", required_argument, 0, 'L'}, /* maxtotallen_bound */
  {"chimera-margin", required_argument, 0, 'x'}, /* chimera_margin */
#if 0
  {"reference", required_argument, 0, 'w'}, /* referencefile */
#else
  {"localsplicedist", required_argument, 0, 'w'}, /* shortsplicedist */
#endif
#ifdef HAVE_PTHREAD
  {"nthreads", required_argument, 0, 't'}, /* nworkers */
#endif
  {"splicingdir", required_argument, 0, 0}, /* user_splicingdir */
  {"nosplicing", no_argument, 0, 0},	    /* novelsplicingp */
  {"use-splicing", required_argument, 0, 's'}, /* splicing_iit, knownsplicingp (was previously altstrainp) */
  {"chrsubsetfile", required_argument, 0, 'C'}, /* user_chrsubsetfile */
  {"chrsubset", required_argument, 0, 'c'}, /* user_chrsubsetname */
  {"trimendexons", required_argument, 0, 'H'}, /* minendexon */
  {"canonical-mode", required_argument, 0, 0}, /* canonical_mode */
  {"cross-species", no_argument, 0, 0}, /* use_shifted_canonical_p */
#ifndef PMAP
  {"prunelevel", required_argument, 0, 'p'}, /* prune_poor_p, prune_repetitive_p */
#endif
  {"allow-close-indels", required_argument, 0, 0}, /* close_indels_mode, extraband_single */
  {"microexon-spliceprob", required_argument, 0, 0}, /* microexon_spliceprob */
  {"stage2-start", required_argument, 0, 0},	     /* suboptimal_score_start */
  {"stage2-end", required_argument, 0, 0},	     /* suboptimal_score_end */

  /* Output options */
  {"output-buffer-size", required_argument, 0, 0}, /* output_buffer_size */
  {"summary", no_argument, 0, 'S'}, /* printtype */
  {"align", no_argument, 0, 'A'}, /* printtype */
  {"continuous", no_argument, 0, '3'}, /* printtype */
  {"continuous-by-exon", no_argument, 0, '4'}, /* printtype */
  {"noexceptions", no_argument, 0, '0'}, /* exception_raise_p */
  {"graphic", no_argument, 0, '6'}, /* debug_graphic_p */
  {"stage3debug", required_argument, 0, '8'}, /* stage3debug, diagnosticp */
  {"check", no_argument, 0, '9'}, /* checkp */
  {"npaths", required_argument, 0, 'n'}, /* maxpaths */
  {"quiet-if-excessive", no_argument, 0, 0}, /* quiet_if_excessive_p */
  {"format", required_argument, 0, 'f'}, /* printtype */
  {"failsonly", no_argument, 0, 0}, /* failsonlyp */
  {"nofails", no_argument, 0, 0}, /* nofailsp */
  {"fails-as-input", no_argument, 0, 0}, /* fails_as_input_p */
  {"split-output", required_argument, 0, 0}, /* sevenway_root */
  {"suboptimal-score", required_argument, 0, 0}, /* suboptimal_score */

#ifndef PMAP
  {"quality-protocol", required_argument, 0, 0}, /* quality_shift */
  {"quality-print-shift", required_argument, 0, 'j'}, /* quality_shift */
  {"no-sam-headers", no_argument, 0, 0},	/* sam_headers_p */
  {"read-group-id", required_argument, 0, 0},	/* sam_read_group_id */
  {"read-group-name", required_argument, 0, 0},	/* sam_read_group_name */
  {"read-group-library", required_argument, 0, 0}, /* sam_read_group_library */
  {"read-group-platform", required_argument, 0, 0}, /* sam_read_group_platform */
#endif

  {"compress", no_argument, 0, 'Z'}, /* printtype */
  {"ordered", no_argument, 0, 'O'}, /* orderedp */
  {"md5", no_argument, 0, '5'}, /* checksump */
  {"chimera-overlap", required_argument, 0, 'o'}, /* chimera_overlap */
  {"use-snps", required_argument, 0, 'V'}, /* snps_root */

  /* Map file options */
  {"mapdir", required_argument, 0, 'M'}, /* user_mapdir */
  {"map", required_argument, 0, 'm'},	/* map_iitfile */
  {"mapexons", no_argument, 0, 'e'}, /* map_exons_p */
  {"mapboth", no_argument, 0, 'b'}, /* map_bothstrands_p */
  {"nflanking", required_argument, 0, 'u'}, /* nflanking */
  {"print-comment", no_argument, 0, 0},	    /* print_comment_p */

  /* Alignment options */
  {"exons", required_argument, 0, 'E'}, /* printtype */
#ifdef PMAP
  {"protein_gen", no_argument, 0, 'P'}, /* printtype */
  {"nucleotide", no_argument, 0, 'Q'}, /* printtype */
#else
  {"protein_dna", no_argument, 0, 'P'}, /* printtype */
  {"protein_gen", no_argument, 0, 'Q'}, /* printtype */
  {"fulllength", no_argument, 0, 'F'}, /* fulllengthp */
  {"cdsstart", required_argument, 0, 'a'}, /* cds_startpos */
  {"truncate", no_argument, 0, 'T'}, /* truncatep */
  {"direction", required_argument, 0, 'z'}, /* sense_try, sense_filter */
#endif
  {"tolerant", no_argument, 0, 'Y'}, /* strictp */
  {"nolengths", no_argument, 0, 'N'},	/* nointronlenp */
  {"invertmode", required_argument, 0, 'I'}, /* invertmode */
  {"introngap", required_argument, 0, 'i'}, /* ngap */
  {"wraplength", required_argument, 0, 'l'}, /* wraplength */
  
  /* Help options */
  {"version", no_argument, 0, 0}, /* print_program_version */
  {"help", no_argument, 0, 0}, /* print_program_usage */
  {0, 0, 0, 0}
};


static void
print_program_version () {
  char *genomedir;

  fprintf(stdout,"\n");
#ifdef PMAP
  fprintf(stdout,"PMAP: Protein Mapping and Alignment Program\n");
#else
  fprintf(stdout,"GMAP: Genomic Mapping and Alignment Program\n");
#endif
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
  fprintf(stdout,"sigaction available\n");
#else
  fprintf(stdout,"no sigaction\n");
#endif
#ifdef PMAP
  fprintf(stdout,"Stage 1 index size: %d aa\n",index1part_aa);
#endif
  fprintf(stdout,"Sizes: off_t (%lu), size_t (%lu), unsigned int (%lu), long int (%lu)\n",
	  sizeof(off_t),sizeof(size_t),sizeof(unsigned int),sizeof(long int));
  fprintf(stdout,"Default gmap directory (compiled): %s\n",GMAPDB);
  genomedir = Datadir_find_genomedir(/*user_genomedir*/NULL);
  fprintf(stdout,"Default gmap directory (environment): %s\n",genomedir);
  FREE(genomedir);
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


static Stage3_T *
stage3array_from_list (int *npaths, int *second_absmq, List_T stage3list, bool chimerap, bool remove_overlaps_p) {
  Stage3_T *array1, *array0, x, y;
  bool *eliminate;
  int norig, i, j;
  int threshold_score;


  Stage3_recompute_goodness(stage3list);

  if ((norig = List_length(stage3list)) == 0) {
    *second_absmq = 0;
    return (Stage3_T *) NULL;

  } else if (chimerap == true) {
    array0 = (Stage3_T *) List_to_array(stage3list,NULL);
    List_free(&stage3list);
    if (norig <= 2) {
      *second_absmq = 0;
    } else {
      qsort(&(array0[2]),norig-2,sizeof(Stage3_T),Stage3_cmp);
      *second_absmq = Stage3_absmq_score(array0[2]);
    }
    *npaths = norig;
    return array0;

  } else if (remove_overlaps_p == false) {
    array0 = (Stage3_T *) List_to_array(stage3list,NULL);
    List_free(&stage3list);
    qsort(array0,norig,sizeof(Stage3_T),Stage3_cmp);

    threshold_score = Stage3_goodness(array0[0]) - suboptimal_score;
    i = 1;
    while (i < norig && Stage3_goodness(array0[i]) >= threshold_score) {
      i++;
    }
    *npaths = i;

    if (*npaths < 2) {
      *second_absmq = 0;
    } else {
      *second_absmq = Stage3_absmq_score(array0[1]);
    }

    return array0;

  } else {
    eliminate = (bool *) CALLOC(norig,sizeof(bool));

    /* Initial sort to remove subsumed alignments */
    array0 = (Stage3_T *) List_to_array(stage3list,NULL);
    List_free(&stage3list);
    qsort(array0,norig,sizeof(Stage3_T),Stage3_cmp);
    for (i = 0; i < norig; i++) {
      x = array0[i];
      for (j = i+1; j < norig; j++) {
	y = array0[j];
	if (Stage3_overlap(x,y)) {
	  eliminate[j] = true;
	}
      }
    }

    *npaths = 0;
    for (i = 0; i < norig; i++) {
      if (eliminate[i] == false) {
	(*npaths)++;
      }
    }

    array1 = (Stage3_T *) CALLOC(*npaths,sizeof(Stage3_T));
    j = 0;
    for (i = 0; i < norig; i++) {
      x = array0[i];
      if (eliminate[i] == true) {
	Stage3_free(&x,/*free_pairarray_p*/true);
      } else {
	array1[j++] = x;
      }
    }
    FREE(array0);
    FREE(eliminate);

    threshold_score = Stage3_goodness(array1[0]) - suboptimal_score;
    i = 1;
    while (i < *npaths && Stage3_goodness(array1[i]) >= threshold_score) {
      i++;
    }
    *npaths = i;

    if (*npaths < 2) {
      *second_absmq = 0;
    } else {
      *second_absmq = Stage3_absmq_score(array1[1]);
    }
    return array1;
  }
}


static List_T
update_stage3list (List_T stage3list, bool lowidentityp, Sequence_T queryseq,
#ifdef PMAP
		   Sequence_T queryntseq,
#endif
		   Sequence_T queryuc, Sequence_T genomicseg, 
		   Genomicpos_T genomicstart, Genomicpos_T genomiclength,
		   Oligoindex_T *oligoindices_major, int noligoindices_major,
		   Oligoindex_T *oligoindices_minor, int noligoindices_minor,
		   Pairpool_T pairpool, Diagpool_T diagpool, int straintype, char *strain,
		   Chrnum_T chrnum, Genomicpos_T chroffset, Genomicpos_T chrpos,
		   bool watsonp, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
		   Stopwatch_T worker_stopwatch) {
  bool do_final_p;
  int stage2_source, stage2_indexsize;

  Sequence_T genomicuc;
  Genomicpos_T genomicend;
  List_T all_paths, path, p;
  Stage3_T stage3;

  struct Pair_T *pairarray;
  List_T pairs;
  int npairs, cdna_direction, matches, unknowns, mismatches, qopens, qindels, topens, tindels,
    ncanonical, nsemicanonical, nnoncanonical;
  int nmatches_pretrim, nmatches_posttrim;
  int sensedir;
  int ambig_end_length_5, ambig_end_length_3;
  Splicetype_T ambig_splicetype_5, ambig_splicetype_3;
  double defect_rate;
  double stage3_runtime;

  if (user_genomicseg == NULL && uncompressedp == false && straintype == 0) {
    genomicuc = Sequence_alias(genomicseg);
  } else {
    genomicuc = Sequence_uppercase(genomicseg);
  }

  debug(printf("Beginning Stage2_compute with genomiclength %d\n",Sequence_fulllength(genomicseg)));
  if (genome == NULL) {
    genomicend = 0U;
  } else {
    genomicend = genomicstart + Sequence_fulllength(genomicseg);
  }

  all_paths = Stage2_compute(&stage2_source,&stage2_indexsize,
			     Sequence_trimpointer(queryseq),Sequence_trimpointer(queryuc),
			     Sequence_trimlength(queryseq),/*query_offset*/0,

			     Sequence_fullpointer(genomicseg),Sequence_fullpointer(genomicuc),
			     genomicstart,genomicend,/*mappingstart*/genomicstart,/*mappingend*/genomicend,
			     /*plusp*/watsonp,/*genomiclength*/Sequence_fulllength(genomicseg),
			     /*genomic_offset*/0,

			     oligoindices_major,noligoindices_major,/*proceed_pctcoverage*/0.5,
			     pairpool,diagpool,sufflookback,nsufflookback,maxintronlen_bound,
			     /*localp*/true,/*skip_repetitive_p*/true,use_shifted_canonical_p,
			     /*favor_right_p*/false,/*just_one_p*/false,debug_graphic_p,
			     diagnosticp,worker_stopwatch,diag_debug);

  for (p = all_paths; p != NULL; p = List_next(p)) {
    path = (List_T) List_head(p);
    if (diag_debug == true) {
      stage3list = path;		/* really diagonals */

    } else if (path != NULL) {
      debug(printf("Beginning Stage3_compute\n"));

      if (canonical_mode == 0) {
	do_final_p = false;
      } else if (canonical_mode == 1) {
	do_final_p = true;
      } else if (lowidentityp == false) {
	do_final_p = false;
      } else {
	do_final_p = true;
      }

      Stopwatch_start(worker_stopwatch);
      pairarray = Stage3_compute(&pairs,&npairs,&cdna_direction,&sensedir,&matches,
				 &nmatches_pretrim,&nmatches_posttrim,
				 &ambig_end_length_5,&ambig_end_length_3,
				 &ambig_splicetype_5,&ambig_splicetype_3,
				 &unknowns,&mismatches,&qopens,&qindels,&topens,&tindels,
				 &ncanonical,&nsemicanonical,&nnoncanonical,&defect_rate,
				 path,genomiclength,
#ifdef PMAP
				 /*queryaaseq_ptr*/Sequence_fullpointer(queryseq),
				 /*queryseq_ptr*/Sequence_fullpointer(queryntseq),
				 /*queryuc_ptr*/Sequence_fullpointer(queryntseq),
				 /*querylength*/Sequence_fulllength(queryntseq),
				 /*skiplength*/Sequence_skiplength(queryntseq),
				 /*query_subseq_offset*/Sequence_subseq_offset(queryntseq),
#else
				 /*queryseq_ptr*/Sequence_fullpointer(queryseq),
				 /*queryuc_ptr*/Sequence_fullpointer(queryuc),
				 /*querylength*/Sequence_fulllength(queryseq),
				 /*skiplength*/Sequence_skiplength(queryseq),
				 /*query_subseq_offset*/Sequence_subseq_offset(queryseq),
#endif
				 /*genomicseg_ptr*/Sequence_fullpointer(genomicseg),
				 /*genomicuc_ptr*/Sequence_fullpointer(genomicuc),
				 chrnum,chroffset,chrpos,
				 /*knownsplice_limit_low*/0U,/*knownsplice_limit_high*/-1U,
				 genome,/*usersegment_p*/usersegment ? true : false,
				 watsonp,/*jump_late_p*/watsonp ? false : true,
				 maxpeelback,maxpeelback_distalmedial,nullgap,
				 extramaterial_end,extramaterial_paired,
				 extraband_single,extraband_end,extraband_paired,
				 minendexon,pairpool,dynprogL,dynprogM,dynprogR,ngap,
				 stage3debug,diagnosticp,checkp,do_final_p,sense_try,sense_filter,
				 oligoindices_minor,noligoindices_minor,diagpool,
				 sufflookback,nsufflookback,maxintronlen_bound,close_indels_mode,
				 /*paired_favor_mode*/0,/*zero_offset*/0);
      stage3_runtime = Stopwatch_stop(worker_stopwatch);
      if (pairarray == NULL) {
	/* Skip */
      } else if (matches < min_matches) {
	FREE_OUT(pairarray);
      } else if ((stage3 = Stage3_new(pairarray,pairs,npairs,cdna_direction,genomicstart,genomiclength,
				      stage2_source,stage2_indexsize,matches,unknowns,mismatches,
				      qopens,qindels,topens,tindels,ncanonical,nsemicanonical,nnoncanonical,
				      defect_rate,chrnum,chroffset,chrpos,watsonp,
				      /*skiplength*/Sequence_skiplength(queryseq),
				      /*trimlength*/Sequence_trimlength(queryseq),
				      stage3_runtime,straintype,strain,altstrain_iit)) != NULL) {
	stage3list = List_push(stage3list,(void *) stage3);
      }
    }
  }

  List_free(&all_paths);

  Sequence_free(&genomicuc);

  return stage3list;
}

#if 0
static List_T
update_stage3list_maponlyp (List_T stage3list, Gregion_T gregion, Sequence_T queryseq, 
#ifdef PMAP
			    Sequence_T queryntseq,
#endif
			    Sequence_T queryuc, Pairpool_T pairpool, int straintype, char *strain, Genome_T genome,
			    Chrnum_T chrnum,  Genomicpos_T chroffset, Genomicpos_T chrpos, Genomicpos_T chrlength,
			    bool watsonp, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR) {
  Stage3_T stage3;

  if ((stage3 = Stage3_direct(gregion,
#ifdef PMAP
			      queryseq,queryntseq,queryntseq,
#else
			      queryseq,queryuc,
#endif
			      pairpool,genome,chrnum,chroffset,chrpos,watsonp,ngap,
			      dynprogL,dynprogR,extramaterial_end,extraband_end)) != NULL) {
    stage3list = List_push(stage3list,stage3);
  }

  return stage3list;
}
#endif


#if 0
/* This code is duplicated in get-genome.c */
static int
index_compare (const void *a, const void *b) {
  int index1 = * (int *) a;
  int index2 = * (int *) b;
  int type1, type2;
  Genomicpos_T pos1, pos2;

  type1 = Interval_type(IIT_interval(altstrain_iit,index1));
  type2 = Interval_type(IIT_interval(altstrain_iit,index2));
  
  if (type1 < type2) {
    return -1;
  } else if (type1 > type2) {
    return +1;
  } else {
    /* Store in descending genomic position, so right shifting works
       in Genome_patch_strain */
    pos1 = Interval_low(IIT_interval(altstrain_iit,index1));
    pos2 = Interval_low(IIT_interval(altstrain_iit,index2));

    if (pos1 > pos2) {
      return -1;
    } else if (pos1 < pos2) {
      return +1;
    } else {
      return 0;
    }
  }
}
#endif


static Stage3_T *
stage3_from_usersegment (int *npaths, int *second_absmq, bool lowidentityp, Sequence_T queryseq,
			 Sequence_T queryuc, Sequence_T usersegment,
			 Oligoindex_T *oligoindices_major, int noligoindices_major,
			 Oligoindex_T *oligoindices_minor, int noligoindices_minor,
			 Pairpool_T pairpool, Diagpool_T diagpool,
			 Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
			 Stopwatch_T worker_stopwatch) {
  List_T stage3list;
  Genomicpos_T chroffset, chrpos, chrlength;
  Sequence_T revcomp;
  Chrnum_T chrnum = 0;

#ifdef PMAP
  Sequence_T queryntseq;
  queryntseq = Sequence_convert_to_nucleotides(queryseq);
#endif
		    
  chroffset = chrpos = 0U;
  chrlength = Sequence_fulllength(usersegment);

  stage3list = update_stage3list(/*stage3list*/NULL,lowidentityp,queryseq,
#ifdef PMAP
				 queryntseq,
#endif
				 queryuc,usersegment,/*genomicstart*/0U,
				 /*genomiclength*/Sequence_fulllength(usersegment),
				 oligoindices_major,noligoindices_major,oligoindices_minor,noligoindices_minor,
				 pairpool,diagpool,/*straintype*/0,/*strain*/NULL,chrnum,chroffset,chrpos,
				 /*watsonp*/true,dynprogL,dynprogM,dynprogR,worker_stopwatch);

  revcomp = Sequence_revcomp(usersegment);

  stage3list = update_stage3list(stage3list,lowidentityp,queryseq,
#ifdef PMAP
				 queryntseq,
#endif
				 queryuc,revcomp,/*genomicstart*/0U,
				 /*genomiclength*/Sequence_fulllength(usersegment),
				 oligoindices_major,noligoindices_major,oligoindices_minor,noligoindices_minor,
				 pairpool,diagpool,/*straintype*/0,/*strain*/NULL,chrnum,chroffset,chrpos,
				 /*watsonp*/false,dynprogL,dynprogM,dynprogR,worker_stopwatch);

  Sequence_free(&revcomp);

#ifdef PMAP
  Sequence_free(&queryntseq);
#endif

  if (stage3list == NULL) {
    *npaths = 0;
    return NULL;
  } else {
    return stage3array_from_list(&(*npaths),&(*second_absmq),stage3list,/*chimerap*/false,/*remove_overlaps_p*/true);
  }
}


static List_T
stage3list_sort (List_T stage3list) {
  List_T sorted = NULL;
  Stage3_T *array;
  int n, i;

  if ((n = List_length(stage3list)) == 0) {
    return (List_T) NULL;
  } else if (n == 1) {
    return stage3list;
  } else {
    array = (Stage3_T *) List_to_array(stage3list,NULL);
    List_free(&stage3list);
    qsort(array,n,sizeof(Stage3_T),Stage3_cmp);
    for (i = n-1; i >= 0; i--) {
      sorted = List_push(sorted,(void *) array[i]);
    }
    FREE(array);

    return sorted;
  }
}


static List_T
stage3_from_gregions (List_T stage3list, List_T gregions, bool lowidentityp, Sequence_T queryseq,
		      Sequence_T queryuc, Sequence_T usersegment, 
		      Oligoindex_T *oligoindices_major, int noligoindices_major,
		      Oligoindex_T *oligoindices_minor, int noligoindices_minor,
		      Pairpool_T pairpool, Diagpool_T diagpool,
		      Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
		      Stopwatch_T worker_stopwatch) {
  Gregion_T gregion, *array;
  char *strain;
  Sequence_T genomicseg, genomicuc;
  int ngregions, ncovered, max_ncovered, stage2_source;
  int i;
#if 0
  int *indexarray, nindices, straintype, j;
#endif
  void *item;

#ifdef PMAP
  Sequence_T queryntseq;
  queryntseq = Sequence_convert_to_nucleotides(queryseq);
#endif
		    
  if (usersegment == NULL && (ngregions = List_length(gregions)) > 0) {
    array = (Gregion_T *) List_to_array(gregions,NULL);
    List_free(&gregions);

    for (i = 0; i < ngregions; i++) {
      gregion = array[i];
      genomicseg = Genome_get_segment(genome,Gregion_genomicstart(gregion),Gregion_genomiclength(gregion),
				      /*chromosome_iit*/NULL,Gregion_revcompp(gregion));
      genomicuc = Sequence_uppercase(genomicseg);
      ncovered = Stage2_scan(&stage2_source,Sequence_trimpointer(queryuc),Sequence_trimlength(queryseq),
			     Sequence_fullpointer(genomicuc),Sequence_fulllength(genomicseg),
			     oligoindices_major,noligoindices_major,diagpool,
			     debug_graphic_p,diagnosticp);
      Gregion_set_ncovered(gregion,ncovered,stage2_source);
      if (diagnosticp == true) {
	fprintf(stderr,"Scanned %d ncovered\n",ncovered);
      }
      Sequence_free(&genomicuc);
      Sequence_free(&genomicseg);
    }
    qsort(array,ngregions,sizeof(Gregion_T),Gregion_cmp);
    max_ncovered = Gregion_ncovered(array[0]);

    gregions = (List_T) NULL;
    i = 0;
    while (i < ngregions && Gregion_ncovered(array[i]) > 0.80*max_ncovered) {
      gregions = List_push(gregions,(void *) array[i]);
      if (diagnosticp == true) {
	fprintf(stderr,"Keeping %d ncovered relative to %d\n",Gregion_ncovered(array[i]),max_ncovered);
      }
      i++;
    }
    while (i < ngregions) {
      Gregion_free(&(array[i]));
      i++;
    }
    FREE(array);
  }

  while (gregions != NULL) {
    gregions = List_pop(gregions,&item);
    gregion = (Gregion_T) item;

    /* if (Match_usep(match) == true) { */
    if (1) {
      if (usersegment != NULL) {
	/* chrlength = Sequence_fulllength(usersegment); */
	strain = NULL;
	genomicseg = Sequence_substring(usersegment,Gregion_genomicstart(gregion),Gregion_genomiclength(gregion),
					Gregion_revcompp(gregion));
	stage3list = update_stage3list(stage3list,lowidentityp,queryseq,
#ifdef PMAP
				       queryntseq,
#endif
				       queryuc,genomicseg,Gregion_genomicstart(gregion),
				       Gregion_genomiclength(gregion),oligoindices_major,noligoindices_major,
				       oligoindices_minor,noligoindices_minor,pairpool,diagpool,
				       /*straintype*/0,/*strain*/NULL,Gregion_chrnum(gregion),
				       Gregion_chroffset(gregion),Gregion_chrpos(gregion),Gregion_plusp(gregion),
				       dynprogL,dynprogM,dynprogR,worker_stopwatch);
	Sequence_free(&genomicseg);

      } else if (maponlyp == true) {
	fprintf(stderr,"maponlyp mode not currently supported\n");
	exit(9);
#if 0
	stage3list = update_stage3list_maponlyp(stage3list,gregion,queryseq,
#ifdef PMAP
						queryntseq,
#endif
						queryuc,pairpool,/*straintype*/0,/*strain*/NULL,genome,
						Gregion_chrnum(gregion),Gregion_chroffset(gregion),
						Gregion_chrpos(gregion),Gregion_chrlength(gregion),Gregion_plusp(gregion),
						dynprogL,dynprogM,dynprogR);
#endif

      } else {
#if 0
	if (diagnosticp == true) {
	  printf("Got sequence at %u with length %u, revcomp %d\n",
		 Gregion_genomicstart(gregion),Gregion_genomiclength(gregion),Gregion_revcompp(gregion));
	}
#endif
	if (genomealt != NULL) {
	  genomicseg = Genome_get_segment_alt(genomealt,Gregion_genomicstart(gregion),Gregion_genomiclength(gregion),
					      /*chromosome_iit*/NULL,Gregion_revcompp(gregion));
	} else {
	  genomicseg = Genome_get_segment(genome,Gregion_genomicstart(gregion),Gregion_genomiclength(gregion),
					  /*chromosome_iit*/NULL,Gregion_revcompp(gregion));
	}
	stage3list = update_stage3list(stage3list,lowidentityp,queryseq,
#ifdef PMAP
				       queryntseq,
#endif
				       queryuc,genomicseg,Gregion_genomicstart(gregion),
				       Gregion_genomiclength(gregion),oligoindices_major,noligoindices_major,
				       oligoindices_minor,noligoindices_minor,pairpool,diagpool,
				       /*straintype*/0,/*strain*/NULL,Gregion_chrnum(gregion),
				       Gregion_chroffset(gregion),Gregion_chrpos(gregion),Gregion_plusp(gregion),
				       dynprogL,dynprogM,dynprogR,worker_stopwatch);
	Sequence_free(&genomicseg);

#if 0
	/* We rely upon the fact that gbuffer1 still holds the genomic segment.  This code is duplicated in get-genome.c */
	if (altstrain_iit != NULL) {
	  indexarray = IIT_get(&nindices,altstrain_iit,/*divstring*/NULL,Gregion_genomicstart(gregion)+1U,
			       Gregion_genomicstart(gregion)+Gregion_genomiclength(gregion)-1,/*sortp*/false);
	  if (nindices > 0) {
	    /* Sort according to type and genome position */
	    qsort(indexarray,nindices,sizeof(int),index_compare);
	    j = 0;
	    while (j < nindices) {
	      i = j++;
	      straintype = Interval_type(IIT_interval(altstrain_iit,indexarray[i]));
	      strain = IIT_typestring(altstrain_iit,straintype);
	      while (j < nindices && Interval_type(IIT_interval(altstrain_iit,indexarray[j])) == straintype) {
		j++;
	      }
	      /* Patch from i to j */
	      genomicseg = Genome_patch_strain(&(indexarray[i]),j-i,altstrain_iit,
					       Gregion_genomicstart(gregion),Gregion_genomiclength(gregion),
					       Gregion_revcompp(gregion),
					       Gbuffer_chars1(gbuffer),Gbuffer_chars2(gbuffer),Gbuffer_chars3(gbuffer),
					       Gbuffer_gbufferlen(gbuffer));
	      stage3list = update_stage3list(stage3list,lowidentityp,queryseq,
#ifdef PMAP
					     queryntseq,
#endif					     
					     queryuc,genomicseg,Gregion_genomicstart(gregion),
					     Gregion_genomiclength(gregion),oligoindices_major,noligoindices_major,
					     oligoindices_minor,noligoindices_minor,pairpool,diagpool,
					     straintype,strain,Gregion_chrnum(gregion),Gregion_chroffset(gregion),
					     Gregion_chrpos(gregion),Gregion_plusp(gregion),
					     dynprogL,dynprogM,dynprogR,worker_stopwatch);
	      Sequence_free(&genomicseg);
	    }
	    FREE(indexarray);
	  }
	}
#endif

      }
    }
    Gregion_free(&gregion);
  }
	
#ifdef PMAP
  Sequence_free(&queryntseq);
#endif

  

  return stage3list;		/* if diag_debug == true, really diagonals */
}


#if 0
static bool
chimeric_join5_p (Stage3_T from, int effective_start, int effective_end) {
  debug2(printf("From %d..%d -> effective %d..%d ?\n",
		Stage3_querystart(from),Stage3_queryend(from),effective_start,effective_end));
  debug2(printf("Checking that effective_end %d > endfrom %d\n",
		effective_end,Stage3_queryend(from)));
  debug2(printf("Checking that effective_start %d > startfrom %d\n",
		effective_start,Stage3_querystart(from)));
  debug2(printf("Checking that endfrom %d - effective_start %d < %d\n",
		Stage3_queryend(from),effective_start,CHIMERA_SLOP));
  debug2(printf("Checking that effective_start %d - endfrom %d < %d\n",
		effective_start,Stage3_queryend(from),CHIMERA_SLOP));

  if (effective_end > Stage3_queryend(from) &&
      effective_start > Stage3_querystart(from) &&
      Stage3_queryend(from) - effective_start < CHIMERA_SLOP &&
      effective_start - Stage3_queryend(from) < CHIMERA_SLOP) {
    debug2(printf("returning true\n\n"));
    return true;
  } else {
    debug2(printf("returning false\n\n"));
    return false;
  }
}

static bool
chimeric_join3_p (Stage3_T to, int effective_start, int effective_end) {
  debug2(printf("Effective %d..%d -> to %d..%d ?\n",
		effective_start,effective_end,Stage3_querystart(to),Stage3_queryend(to)));
  debug2(printf("Checking that endto %d > effective_end %d\n",
		Stage3_queryend(to),effective_end));
  debug2(printf("Checking that startto %d > effective_start %d\n",
		Stage3_querystart(to),effective_start));
  debug2(printf("Checking that effective_end %d - startto %d < %d\n",
		effective_end,Stage3_querystart(to),CHIMERA_SLOP));
  debug2(printf("Checking that startto %d - effective_end %d < %d\n",
		Stage3_querystart(to),effective_end,CHIMERA_SLOP));

  if (Stage3_queryend(to) > effective_end &&
      Stage3_querystart(to) > effective_start &&
      effective_end - Stage3_querystart(to) < CHIMERA_SLOP &&
      Stage3_querystart(to) - effective_end < CHIMERA_SLOP) {
    debug2(printf("returning true\n\n"));
    return true;
  } else {
    debug2(printf("returning false\n\n"));
    return false;
  }
}
#endif

static bool
chimeric_join_p (Stage3_T from, Stage3_T to) {
  debug2(printf("from %d..%d (%u..%u) -> to %d..%d (%u..%u) => ",
		Stage3_querystart(from),Stage3_queryend(from),
		Stage3_genomicstart(from),Stage3_genomicend(from),
		Stage3_querystart(to),Stage3_queryend(to),
		Stage3_genomicstart(to),Stage3_genomicend(to)));

  if (Stage3_queryend(from) - Stage3_querystart(to) < CHIMERA_SLOP &&
      Stage3_querystart(to) - Stage3_queryend(from) < CHIMERA_SLOP) {
    debug2(printf("true\n"));
    return true;
  } else {
    debug2(printf("false\n"));
    return false;
  }
}

/*
static bool
chimera_exists_p (Stage3_T *stage3array, int npaths, Genome_T genome, int querylength) {
  int i;
  int exonexonpos, cdna_direction;
  double donor_prob, acceptor_prob;

  for (i = 1; i < npaths; i++) {
    if (chimeric_join_p(stage3array[0],stage3array[i]) == true) {
      debug2(printf("Found join for 0 to %d\n",i));
      if (Chimera_exonexon_p(&exonexonpos,&cdna_direction,&donor_prob,&acceptor_prob,
			     stage3array[0],stage3array[i],genome,querylength) == true) {
	debug2(printf("exonexonpos is %d\n",exonexonpos));
	return true;
      }
    } else if (chimeric_join_p(stage3array[i],stage3array[0]) == true) {
      debug2(printf("Found join for %d to 0\n",i));
      if (Chimera_exonexon_p(&exonexonpos,&cdna_direction,&donor_prob,&acceptor_prob,
			     stage3array[i],stage3array[0],genome,querylength) == true) {
	debug2(printf("exonexonpos is %d\n",exonexonpos));
	return true;
      }
    }
  }
  return false;
}
*/


/* Returns nonjoinable */
static List_T
chimera_separate_paths (Stage3_T **stage3array_sub1, int *npaths_sub1, 
			Stage3_T **stage3array_sub2, int *npaths_sub2,
			List_T stage3list) {
  List_T nonjoinable = NULL, p, q;
  Stage3_T stage3_1, stage3_2, stage3;
  bool *joinable_left, *joinable_right;
  int npaths, i, j, k;

  debug2(printf("chimera_separate_paths called with list length %d\n",List_length(stage3list)));

  if (stage3list == NULL) {
    *stage3array_sub1 = (Stage3_T *) NULL;
    *npaths_sub1 = 0;
    *stage3array_sub2 = (Stage3_T *) NULL;
    *npaths_sub2 = 0;
    return (List_T) NULL;
  }

  npaths = List_length(stage3list);
  joinable_left = (bool *) CALLOC(npaths,sizeof(bool));
  joinable_right = (bool *) CALLOC(npaths,sizeof(bool));

  for (p = stage3list, i = 0; p != NULL; p = List_next(p), i++) {
    stage3_1 = (Stage3_T) List_head(p);

    for (q = List_next(p), j = i+1; q != NULL; q = List_next(q), j++) {
      stage3_2 = (Stage3_T) List_head(q);
      
      if (chimeric_join_p(stage3_1,stage3_2) == true) {
	debug2(printf("Found join from %d to %d\n",i,j));
	joinable_left[i] = true;
	joinable_right[j] = true;
      } else if (chimeric_join_p(stage3_2,stage3_1) == true) {
	debug2(printf("Found join from %d to %d\n",j,i));
	joinable_right[i] = true;
	joinable_left[j] = true;
      }
    }
  }

  *npaths_sub1 = *npaths_sub2 = 0;
  for (i = 0; i < npaths; i++) {
    if (joinable_left[i] == true) {
      (*npaths_sub1)++;
    }
    if (joinable_right[i] == true) {
      (*npaths_sub2)++;
    }
  }

  if (*npaths_sub1 == 0 || *npaths_sub2 == 0) {
    *stage3array_sub1 = (Stage3_T *) NULL;
    *npaths_sub1 = 0;
    *stage3array_sub2 = (Stage3_T *) NULL;
    *npaths_sub2 = 0;
  } else {
    *stage3array_sub1 = (Stage3_T *) CALLOC(*npaths_sub1,sizeof(Stage3_T));
    *stage3array_sub2 = (Stage3_T *) CALLOC(*npaths_sub2,sizeof(Stage3_T));
    j = k = 0;
    for (p = stage3list, i = 0; p != NULL; p = List_next(p), i++) {
      stage3 = (Stage3_T) List_head(p);
      if (joinable_left[i] == false && joinable_right[i] == false) {
	nonjoinable = List_push(nonjoinable,stage3);
      } else {
	/* Note: it is possible that the same stage3 object gets put into both lists */
	if (joinable_left[i] == true) {
	  (*stage3array_sub1)[j++] = stage3;
	}
	if (joinable_right[i] == true) {
	  (*stage3array_sub2)[k++] = stage3;
	}
      }
    }
  }

  FREE(joinable_right);
  FREE(joinable_left);

  return nonjoinable;
}

#if 0
#endif


/* Returns a list with only one Stage3_T object */
static List_T
merge_left_and_right_readthrough (Stage3_T *stage3array_sub1, int npaths_sub1, int bestfrom,
				  Stage3_T *stage3array_sub2, int npaths_sub2, int bestto,
				  char comp, Genomicpos_T genomegap, List_T nonjoinable, int chimerapos, int chimeraequivpos,
				  int queryntlength, Pairpool_T pairpool, Genome_T genome,
				  int ngap) {
  List_T newstage3list, p;
  Stage3_T best0, best1, *array, last, freed0 = NULL, freed1 = NULL;
  int i, k;

  best0 = stage3array_sub1[bestfrom];
  best1 = stage3array_sub2[bestto];
  debug2(printf("bestfrom %d: %p, bestto %d: %p\n",bestfrom,best0,bestto,best1));

  if (Stage3_cdna_direction(best0) != Stage3_cdna_direction(best1) &&
      Stage3_cdna_direction(best0) != 0 &&
      Stage3_cdna_direction(best1) != 0) {
    debug2(printf("cdna_directions are not compatible\n"));
    if (Stage3_npairs(best0) > Stage3_npairs(best1)) {
      newstage3list = (List_T) NULL;
      newstage3list = List_push(newstage3list,(void *) best0);
      freed1 = best1;
      Stage3_free(&best1,/*free_pairarray_p*/true);
    } else {
      newstage3list = (List_T) NULL;
      newstage3list = List_push(newstage3list,(void *) best1);
      freed0 = best0;
      Stage3_free(&best0,/*free_pairarray_p*/true);
    }
  } else {
    Stage3_merge_readthrough(best0,best1,comp,genomegap,/*minpos1*/0,/*maxpos1*/chimeraequivpos+chimera_overlap,
			     /*minpos2*/chimerapos+1-chimera_overlap,/*maxpos2*/queryntlength,
			     pairpool,genome,ngap);
    debug2(printf("Rearranging paths\n"));
    newstage3list = (List_T) NULL;
    newstage3list = List_push(newstage3list,(void *) best0);
    freed1 = best1;
    Stage3_free(&best1,/*free_pairarray_p*/false);
    debug2(printf("Pushing stage3 %p\n",best0));
  }

  if (npaths_sub1 + npaths_sub2 > 2) {
    /* Push rest of results, taking care not to have duplicates */

    array = (Stage3_T *) CALLOC(npaths_sub1 + npaths_sub2 - 2,sizeof(Stage3_T));
    k = 0;
    for (i = 0; i < npaths_sub1; i++) {
      if (i != bestfrom) {
	debug2(printf("array %d is now sub1 %d: %p\n",k,i,stage3array_sub1[i]));
	array[k++] = stage3array_sub1[i];
      }
    }
    for (i = 0; i < npaths_sub2; i++) {
      if (i != bestto) {
	debug2(printf("array %d is now sub2 %d: %p\n",k,i,stage3array_sub2[i]));
	array[k++] = stage3array_sub2[i];
      }
    }
    qsort(array,npaths_sub1+npaths_sub2-2,sizeof(Stage3_T),Stage3_identity_cmp);

    last = (Stage3_T) NULL;
    for (i = 0; i < npaths_sub1+npaths_sub2-2; i++) {
      if (array[i] == last) {
	/* Skip */
	debug2(printf("array %d: Skipping stage3 %p, because just pushed, so duplicate\n",i,array[i]));
      } else if (array[i] == best0 || array[i] == best1) {
	/* Skip */
	debug2(printf("array %d: Skipping stage3 %p, because in chimera\n",i,array[i]));
      } else if (array[i] == freed0 || array[i] == freed1) {
	/* Skip */
	debug2(printf("array %d: Skipping stage3 %p, because already freed\n",i,array[i]));
      } else {
	debug2(printf("array %d: Pushing stage3 %p\n",i,array[i]));
	newstage3list = List_push(newstage3list,(void *) array[i]);
	last = array[i];
      }
    }

    FREE(array);
  }

  for (p = nonjoinable; p != NULL; p = List_next(p)) {
    debug2(printf("Pushing nonjoinable stage3 %p\n",List_head(p)));
    newstage3list = List_push(newstage3list,(void *) List_head(p));
  }

  return List_reverse(newstage3list);
}


/* Returns a list with only two Stage3_T objects */
static List_T
merge_left_and_right_transloc (Stage3_T *stage3array_sub1, int npaths_sub1, int bestfrom,
			       Stage3_T *stage3array_sub2, int npaths_sub2, int bestto,
			       List_T nonjoinable, int chimerapos, int chimeraequivpos,
			       int queryntlength) {
  List_T newstage3list, p;
  Stage3_T best0, best1, *array, last;
  int i, k;

  best0 = stage3array_sub1[bestfrom];
  best1 = stage3array_sub2[bestto];

  Stage3_merge_chimera(best0,best1,/*minpos1*/0,/*maxpos1*/chimeraequivpos+chimera_overlap,
		       /*minpos2*/chimerapos+1-chimera_overlap,/*maxpos2*/queryntlength);

  debug2(printf("Rearranging paths\n"));
  newstage3list = (List_T) NULL;
  newstage3list = List_push(newstage3list,(void *) best0);
  newstage3list = List_push(newstage3list,(void *) best1);
  debug2(printf("Pushing stage3 %p\n",best0));
  debug2(printf("Pushing stage3 %p\n",best1));

  if (npaths_sub1 + npaths_sub2 > 2) {
    /* Push rest of results, taking care not to have duplicates */

    array = (Stage3_T *) CALLOC(npaths_sub1 + npaths_sub2 - 2,sizeof(Stage3_T));
    k = 0;
    for (i = 0; i < npaths_sub1; i++) {
      if (i != bestfrom) {
	array[k++] = stage3array_sub1[i];
      }
    }
    for (i = 0; i < npaths_sub2; i++) {
      if (i != bestto) {
	array[k++] = stage3array_sub2[i];
      }
    }
    qsort(array,npaths_sub1+npaths_sub2-2,sizeof(Stage3_T),Stage3_identity_cmp);

    last = (Stage3_T) NULL;
    for (i = 0; i < npaths_sub1+npaths_sub2-2; i++) {
      if (array[i] == last) {
	/* Skip */
	debug2(printf("Skipping stage3 %p, because just pushed\n",array[i]));
      } else if (array[i] == best0 || array[i] == best1) {
	/* Skip */
	debug2(printf("Skipping stage3 %p, because in chimera\n",array[i]));
      } else {
	debug2(printf("Pushing stage3 %p\n",array[i]));
	newstage3list = List_push(newstage3list,(void *) array[i]);
	last = array[i];
      }
    }

    FREE(array);
  }

  for (p = nonjoinable; p != NULL; p = List_next(p)) {
    debug2(printf("Pushing nonjoinable stage3 %p\n",List_head(p)));
    newstage3list = List_push(newstage3list,(void *) List_head(p));
  }

  return List_reverse(newstage3list);
}


static List_T
check_for_chimera (Chimera_T *chimera, List_T stage3list, int effective_start, int effective_end,
		   Sequence_T queryseq, Sequence_T queryuc, Sequence_T usersegment, 
		   Oligoindex_T *oligoindices_major, int noligoindices_major,
		   Oligoindex_T *oligoindices_minor, int noligoindices_minor,
		   Matchpool_T matchpool, Pairpool_T pairpool, 
		   Diagpool_T diagpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR) {
  List_T gregions = NULL, nonjoinable = NULL;
  Stage3_T *stage3array_sub1 = NULL, *stage3array_sub2 = NULL;
  Sequence_T querysubseq = NULL, querysubuc = NULL;
  Diagnostic_T diagnostic;
  int chimerapos, chimeraequivpos;
  int bestfrom, bestto;
  int five_margin, three_margin, five_score = 0, three_score = 0;
  int npaths_sub1 = 0, npaths_sub2 = 0;
  int queryntlength;
  bool lowidentityp;

#ifdef DEBUG2
  List_T p;
  Stage3_T stage3;
#endif

  int exonexonpos, cdna_direction;
  char donor1, donor2, acceptor2, acceptor1;
  double donor_prob, acceptor_prob;
  char comp;
  Genomicpos_T genomegap;

  five_margin = effective_start - Sequence_trim_start(queryseq);
  three_margin = Sequence_trim_end(queryseq) - effective_end;
  debug2(printf("Margins are %d on the 5' end and %d on the 3' end\n",
		five_margin,three_margin));

#ifdef DEBUG2
  for (p = stage3list; p != NULL; p = List_next(p)) {
    stage3 = (Stage3_T) List_head(p);
    Pair_dump_array(Stage3_pairarray(stage3),Stage3_npairs(stage3),/*zerobasedp*/true);
    printf("\n");
  }
#endif

  queryntlength = Sequence_ntlength(queryseq);
  nonjoinable = chimera_separate_paths(&stage3array_sub1,&npaths_sub1,&stage3array_sub2,&npaths_sub2,
				       stage3list);

  if (npaths_sub1 == 0 && npaths_sub2 == 0) {
    /* Need to compute on margin explicitly */
    if (five_margin < chimera_margin && three_margin < chimera_margin) {
      debug2(printf("Insufficient margins\n"));
    } else if (five_margin > three_margin) {
      if ((querysubseq = Sequence_subsequence(queryseq,0,effective_start+CHIMERA_SLOP)) != NULL) {
	if ((querysubuc = Sequence_subsequence(queryuc,0,effective_start+CHIMERA_SLOP)) != NULL) {
	  debug2(printf("5 margin > 3 margin.  "));
	  debug2(printf("Beginning Stage1_compute on 5' margin (%d..%d)\n",0,effective_start+CHIMERA_SLOP));
	  debug2(Sequence_print(stdout,querysubseq,/*uppercasep*/true,wraplength,/*trimmedp*/true));

	  diagnostic = Diagnostic_new();
	  gregions = Stage1_compute(&lowidentityp,querysubuc,
#ifdef PMAP
				    indexdb_fwd,indexdb_rev,
#else
				    indexdb,
#endif
				    indexdb_size_threshold,chromosome_iit,chrsubset,matchpool,
				    maxintronlen_bound,maxtotallen_bound,min_extra_end,
				    stutterhits,diagnostic,/*worker_stopwatch*/NULL);
	  Diagnostic_free(&diagnostic);
	  debug2(printf("Performing Stage 3 starting with list length %d\n",List_length(stage3list)));
	  stage3list = stage3_from_gregions(stage3list,gregions,lowidentityp,querysubseq,
					    querysubuc,usersegment,oligoindices_major,noligoindices_major,
					    oligoindices_minor,noligoindices_minor,pairpool,diagpool,
					    dynprogL,dynprogM,dynprogR,/*worker_stopwatch*/NULL);

	  /* Above function frees gregions */
	  Sequence_free(&querysubuc);

	  List_free(&nonjoinable);
	  nonjoinable = chimera_separate_paths(&stage3array_sub1,&npaths_sub1,&stage3array_sub2,&npaths_sub2,
					       stage3list);
	  debug2(printf("Got %d and %d paths\n",npaths_sub1,npaths_sub2));
	}
	Sequence_free(&querysubseq);
      }
	
    } else {
      if ((querysubseq = Sequence_subsequence(queryseq,effective_end-CHIMERA_SLOP,queryntlength)) != NULL) {
	if ((querysubuc = Sequence_subsequence(queryuc,effective_end-CHIMERA_SLOP,queryntlength)) != NULL) {
	  debug2(printf("5 margin <= 3 margin.  "));
	  debug2(printf("Beginning Stage1_compute on 3' margin (%d..%d)\n",effective_end-CHIMERA_SLOP,queryntlength));
	  debug2(Sequence_print(stdout,querysubseq,/*uppercasep*/true,wraplength,/*trimmedp*/true));

	  diagnostic = Diagnostic_new();
	  gregions = Stage1_compute(&lowidentityp,querysubuc,
#ifdef PMAP
				    indexdb_fwd,indexdb_rev,
#else
				    indexdb,
#endif
				    indexdb_size_threshold,chromosome_iit,chrsubset,matchpool,
				    maxintronlen_bound,maxtotallen_bound,min_extra_end,
				    stutterhits,diagnostic,/*worker_stopwatch*/NULL);
	  Diagnostic_free(&diagnostic);
	  debug2(printf("Performing Stage 3 with list length %d\n",List_length(stage3list)));
	  stage3list = stage3_from_gregions(stage3list,gregions,lowidentityp,querysubseq,
					    querysubuc,usersegment,oligoindices_major,noligoindices_major,
					    oligoindices_minor,noligoindices_minor,pairpool,diagpool,
					    dynprogL,dynprogM,dynprogR,/*worker_stopwatch*/NULL);

	  /* Above function frees gregions */
	  Sequence_free(&querysubuc);

	  List_free(&nonjoinable);
	  nonjoinable = chimera_separate_paths(&stage3array_sub1,&npaths_sub1,&stage3array_sub2,&npaths_sub2,
					       stage3list);
	  debug2(printf("Got %d and %d paths\n",npaths_sub1,npaths_sub2));
	}
	Sequence_free(&querysubseq);
      }
    }
  }

  if (npaths_sub1 == 0 || npaths_sub2 == 0) {
    *chimera = (Chimera_T) NULL;
  } else {

    Chimera_bestpath(&five_score,&three_score,&chimerapos,&chimeraequivpos,&bestfrom,&bestto,
		     stage3array_sub1,npaths_sub1,stage3array_sub2,npaths_sub2,queryntlength);

    /*
    chimeric_goodness = Stage3_chimeric_goodness(&matches0,&matches1,
						 stage3array_sub1[bestfrom],stage3array_sub2[bestto],
						 chimerapos);
    */
    debug2(printf("Chimera_bestpath returns boundary at %d..%d (switch can occur at %d..%d)\n",
		  chimerapos,chimeraequivpos,chimerapos-1,chimeraequivpos));

    Chimera_find_exonexon(&exonexonpos,&cdna_direction,
			  &donor1,&donor2,&acceptor2,&acceptor1,&donor_prob,&acceptor_prob,
			  stage3array_sub1[bestfrom],stage3array_sub2[bestto],
			  genome,chromosome_iit,chimerapos-1,chimeraequivpos);
    debug2(printf("Exon-exon boundary found at %d\n",exonexonpos));
    chimerapos = chimeraequivpos = exonexonpos;


    /* Check to see if we can merge chimeric parts */
    if (Stage3_mergeable(&comp,&genomegap,stage3array_sub1[bestfrom],stage3array_sub2[bestto],
			 exonexonpos,queryntlength,cdna_direction,donor_prob,acceptor_prob) == true) {
      debug2(printf("Mergeable! -- Merging left and right as a readthrough\n"));
      List_free(&stage3list);
      stage3list = merge_left_and_right_readthrough(stage3array_sub1,npaths_sub1,bestfrom,
						    stage3array_sub2,npaths_sub2,bestto,
						    comp,genomegap,nonjoinable,chimerapos,chimeraequivpos,
						    queryntlength,pairpool,genome,ngap);
#if 0
      best1 = stage3array[1];
      Stage3_free(&best1,/*free_pairarray_p*/true);
      for (j = 2; j < *npaths; j++) {
	stage3array[j-1] = stage3array[j];
      }
      stage3array[(*npaths)-1] = (Stage3_T) NULL;
      *npaths -= 1;
      Chimera_free(&(*chimera));
#endif
    } else if (Stage3_test_bounds(stage3array_sub1[bestfrom],0,chimeraequivpos+chimera_overlap) == true &&
	       Stage3_test_bounds(stage3array_sub2[bestto],chimerapos+1-chimera_overlap,queryntlength) == true) {
      *chimera = Chimera_new(chimerapos,chimeraequivpos,exonexonpos,cdna_direction,
			     donor1,donor2,acceptor2,acceptor1,donor_prob,acceptor_prob);
      debug2(printf("Not mergeable -- Merging left and right as a transloc\n"));
      List_free(&stage3list);
      stage3list = merge_left_and_right_transloc(stage3array_sub1,npaths_sub1,bestfrom,
						 stage3array_sub2,npaths_sub2,bestto,
						 nonjoinable,chimerapos,chimeraequivpos,
						 queryntlength);
    }
    FREE(stage3array_sub2);
    FREE(stage3array_sub1);
  }

  List_free(&nonjoinable);

  debug2(printf("check_for_chimera returning list of length %d\n",List_length(stage3list)));
  return stage3list;
}



static List_T
apply_stage3 (Chimera_T *chimera, List_T gregions, bool lowidentityp, Sequence_T queryseq, Sequence_T queryuc,
	      Sequence_T usersegment, Oligoindex_T *oligoindices_major, int noligoindices_major,
	      Oligoindex_T *oligoindices_minor, int noligoindices_minor,
	      Matchpool_T matchpool, Pairpool_T pairpool, 
	      Diagpool_T diagpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
	      Stopwatch_T worker_stopwatch) {
  List_T stage3list;
  Stage3_T nonchimericbest;
  bool testchimerap = false;
  int effective_start, effective_end;

  
  *chimera = NULL;

  stage3list = stage3_from_gregions(/*stage3list*/(List_T) NULL,gregions,lowidentityp,queryseq,queryuc,
				    usersegment,oligoindices_major,noligoindices_major,
				    oligoindices_minor,noligoindices_minor,pairpool,diagpool,
				    dynprogL,dynprogM,dynprogR,worker_stopwatch);
  Stage3_recompute_goodness(stage3list);
  stage3list = stage3list_sort(stage3list);

  debug2(printf("Initial search gives stage3list of length %d\n",List_length(stage3list)));

  if (diag_debug == true) {
    return stage3list;		/* really diagonals */
  }

  if (stage3list != NULL) {
    nonchimericbest = (Stage3_T) List_head(stage3list);

    if (chimera_margin < 0) {
      debug2(printf("turned off\n"));
      testchimerap = false;

    } else if (Stage3_domain(nonchimericbest) < chimera_margin) {
      debug2(printf("Existing alignment is too short, so won't look for chimera\n"));
      testchimerap = false;

#if 0
    } else if (Stage3_fracidentity(nonchimericbest) < CHIMERA_IDENTITY &&
	       Chimera_alignment_break(&effective_start,&effective_end,nonchimericbest,Sequence_ntlength(queryseq),CHIMERA_FVALUE) >= chimera_margin
	       ) {
      debug2(printf("Break in alignment quality at %d..%d detected, so will look for chimera\n",
		    effective_start,effective_end));
      testchimerap = true;
#endif

    } else if (Stage3_largemargin(&effective_start,&effective_end,nonchimericbest,Sequence_ntlength(queryseq)) >= chimera_margin) {
      debug2(printf("Large margin at %d..%d detected, so will look for chimera\n",
		    effective_start,effective_end));
      testchimerap = true;

    } else {
      debug2(printf("Good alignment already with identity %f, so won't look for chimera\n",
		    Stage3_fracidentity(nonchimericbest)));
      testchimerap = false;
    }

    if (testchimerap == true) {
      debug2(printf("Checking for chimera, starting with list length %d\n",List_length(stage3list)));
      stage3list = check_for_chimera(&(*chimera),stage3list,effective_start,effective_end,
				     queryseq,queryuc,usersegment,
				     oligoindices_major,noligoindices_major,oligoindices_minor,noligoindices_minor,
				     matchpool,pairpool,diagpool,dynprogL,dynprogM,dynprogR);
    }
  }

  debug2(printf("apply_stage3 returning list of length %d\n",List_length(stage3list)));
  return stage3list;
}


static Result_T
process_request (Request_T request, Matchpool_T matchpool, Pairpool_T pairpool, Diagpool_T diagpool, 
		 Oligoindex_T *oligoindices_major, int noligoindices_major,
		 Oligoindex_T *oligoindices_minor, int noligoindices_minor,
		 Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
		 Stopwatch_T worker_stopwatch) {
  Result_T result;
  int jobid;
  Diagnostic_T diagnostic;
  Sequence_T queryseq, queryuc;
  Chimera_T chimera = NULL;
  bool lowidentityp;
#ifndef PMAP
  bool repetitivep = false, poorp = false;
#endif

  List_T gregions, stage3list;
  Stage3_T *stage3array;
  int npaths, second_absmq;

  jobid = Request_id(request);
  queryseq = Request_queryseq(request);
  Matchpool_reset(matchpool);
  Pairpool_reset(pairpool);
  Diagpool_reset(diagpool);

  if (Sequence_fulllength_given(queryseq) <= 0) {
    result = Result_new(jobid,(Chimera_T) NULL,(Stage3_T *) NULL,/*npaths*/0,/*second_absmq*/0,
			/*diagnostic*/NULL,EMPTY_SEQUENCE);
      
  } else if (Sequence_fulllength_given(queryseq) < 
#ifdef PMAP
	     index1part_aa
#else
	     index1part
#endif
	     ) {
    result = Result_new(jobid,(Chimera_T) NULL,(Stage3_T *) NULL,/*npaths*/0,/*second_absmq*/0,
			/*diagnostic*/NULL,SHORT_SEQUENCE);

  } else {			/* Sequence_fulllength_given(queryseq) > 0 */
    queryuc = Sequence_uppercase(queryseq);
    diagnostic = Diagnostic_new();
    if (maponlyp == true) {
      diagnostic->query_trim_start = 0;
      diagnostic->query_trim_end = Sequence_fulllength(queryseq);
    } else {
#ifdef PMAP
      Oligoindex_set_inquery(&diagnostic->query_badoligos,&diagnostic->query_repoligos,
			     &diagnostic->query_trimoligos,&diagnostic->query_trim_start,
			     &diagnostic->query_trim_end,oligoindices_major[0],Sequence_fullpointer(queryuc),
			     Sequence_fulllength(queryuc),/*trimp*/false);
#else
      diagnostic->query_oligodepth = 
	Oligoindex_set_inquery(&diagnostic->query_badoligos,&diagnostic->query_repoligos,
			       &diagnostic->query_trimoligos,&diagnostic->query_trim_start,
			       &diagnostic->query_trim_end,oligoindices_major[0],Sequence_fullpointer(queryuc),
			       Sequence_fulllength(queryuc),/*trimp*/true);

      if (diagnostic->query_trimoligos == 0) {
	poorp = true;
      } else if (((double) diagnostic->query_badoligos/(double) diagnostic->query_trimoligos > MAX_BADOLIGOS) ||
		 (diagnostic->query_trim_end - diagnostic->query_trim_start < 80 && diagnostic->query_badoligos > 0)) {
	poorp = true;
      }
#if 0
      if (diagnostic->query_trimoligos == 0) {
	repetitivep = false;
      } else if (diagnostic->query_oligodepth > MAX_OLIGODEPTH || 
		 (double) diagnostic->query_repoligos/(double) diagnostic->query_trimoligos > MAX_REPOLIGOS) {
	repetitivep = true;
      }
#endif
      repetitivep = false;
#endif
    }

#ifndef PMAP
    if (poorp == true && prune_poor_p == true) {
      result = Result_new(jobid,(Chimera_T) NULL,(Stage3_T *) NULL,/*npaths*/0,/*second_absmq*/0,
			  diagnostic,POOR_SEQUENCE);

    } else if (repetitivep == true && prune_repetitive_p == true) {
      result = Result_new(jobid,(Chimera_T) NULL,(Stage3_T *) NULL,/*npaths*/0,/*second_absmq*/0,
			  diagnostic,REPETITIVE);

    }
#endif

    if (usersegment != NULL && userstage1p == false) {
#ifndef PMAP
#if 0
      /* Don't do Sequence_trim, because it affects sequences like NM_018406 */
      Sequence_trim(queryseq,diagnostic->query_trim_start,diagnostic->query_trim_end);
      Sequence_trim(queryuc,diagnostic->query_trim_start,diagnostic->query_trim_end);
#endif
#endif
      stage3array = stage3_from_usersegment(&npaths,&second_absmq,/*lowidentityp*/false,queryseq,queryuc,usersegment,
					    oligoindices_major,noligoindices_major,oligoindices_minor,noligoindices_minor,
					    pairpool,diagpool,dynprogL,dynprogM,dynprogR,worker_stopwatch);
      result = Result_new(jobid,(Chimera_T) NULL,stage3array,npaths,second_absmq,diagnostic,NO_FAILURE);

    } else {		/* Not user segment and not maponly */
#ifndef PMAP
#if 0
      /* Don't do Sequence_trim, because it affects sequences like NM_018406 */
      Sequence_trim(queryseq,diagnostic->query_trim_start,diagnostic->query_trim_end);
      Sequence_trim(queryuc,diagnostic->query_trim_start,diagnostic->query_trim_end);
#endif
#endif

      gregions = Stage1_compute(&lowidentityp,queryuc,
#ifdef PMAP
				indexdb_fwd,indexdb_rev,
#else
				indexdb,
#endif
				indexdb_size_threshold,chromosome_iit,chrsubset,matchpool,
				maxintronlen_bound,maxtotallen_bound,min_extra_end,
				stutterhits,diagnostic,worker_stopwatch);
      if (stage1debug == true) {
	result = Result_new_stage1debug(jobid,gregions,diagnostic,NO_FAILURE);
      } else {
	stage3list = apply_stage3(&chimera,gregions,lowidentityp,queryseq,queryuc,usersegment,
				  oligoindices_major,noligoindices_major,oligoindices_minor,noligoindices_minor,
				  matchpool,pairpool,diagpool,dynprogL,dynprogM,dynprogR,worker_stopwatch);
	if (diag_debug == true) {
	  result = Result_new_diag_debug(jobid,/*diagonals*/stage3list,diagnostic,NO_FAILURE);
	} else if (stage3list == NULL) {
	  result = Result_new(jobid,chimera,/*stage3array*/NULL,/*npaths*/0,/*second_absmq*/0,diagnostic,NO_FAILURE);
	} else if (chimera == NULL) {
	  stage3array = stage3array_from_list(&npaths,&second_absmq,stage3list,/*chimerap*/false,/*remove_overlaps_p*/true);
	  result = Result_new(jobid,/*chimera*/NULL,stage3array,npaths,second_absmq,diagnostic,NO_FAILURE);
	} else {
	  stage3array = stage3array_from_list(&npaths,&second_absmq,stage3list,/*chimerap*/true,/*remove_overlaps_p*/false);
	  result = Result_new(jobid,chimera,stage3array,npaths,second_absmq,diagnostic,NO_FAILURE);
	}
      }

      Oligoindex_clear_inquery(oligoindices_major[0]);

    } /* Matches not user segment and not maponly */

    Sequence_free(&queryuc);
  } /* Matches sequence length > 0 */

  return result;
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

#define POOL_FREE_INTERVAL 200

static void
single_thread () {
  Oligoindex_T *oligoindices_major, *oligoindices_minor;
  int noligoindices_major, noligoindices_minor;
  Dynprog_T dynprogL, dynprogM, dynprogR;
  Matchpool_T matchpool;
  Pairpool_T pairpool;
  Diagpool_T diagpool;
  Stopwatch_T worker_stopwatch;
  Request_T request;
  Result_T result;
  Sequence_T queryseq;
  int noutput = 0;
  int jobid = 0;

  oligoindices_major = Oligoindex_new_major(&noligoindices_major);
  oligoindices_minor = Oligoindex_new_minor(&noligoindices_minor);
  dynprogL = Dynprog_new(nullgap,EXTRAQUERYGAP,maxpeelback,extramaterial_end,extramaterial_paired);
  dynprogM = Dynprog_new(nullgap,EXTRAQUERYGAP,maxpeelback,extramaterial_end,extramaterial_paired);
  dynprogR = Dynprog_new(nullgap,EXTRAQUERYGAP,maxpeelback,extramaterial_end,extramaterial_paired);
  matchpool = Matchpool_new();
  pairpool = Pairpool_new();
  diagpool = Diagpool_new();
  worker_stopwatch = diagnosticp == true ? Stopwatch_new() : (Stopwatch_T) NULL;

  while ((request = Inbuffer_get_request(inbuffer)) != NULL) {
    if (jobid % POOL_FREE_INTERVAL == 0) {
      Pairpool_free_memory(pairpool);
      Diagpool_free_memory(diagpool);
      Matchpool_free_memory(matchpool);
    }
    TRY
      result = process_request(request,matchpool,pairpool,diagpool,
  			       oligoindices_major,noligoindices_major,
			       oligoindices_minor,noligoindices_minor,
			       dynprogL,dynprogM,dynprogR,worker_stopwatch);
    ELSE
      queryseq = Request_queryseq(request);
      if (Sequence_accession(queryseq) == NULL) {
        fprintf(stderr,"Problem with unnamed sequence (%d bp)\n",Sequence_fulllength_given(queryseq));
      } else {
        fprintf(stderr,"Problem with sequence %s (%d bp)\n",
  	      Sequence_accession(queryseq),Sequence_fulllength_given(queryseq));
      }
      fprintf(stderr,"To obtain a core dump, re-run program on problem sequence with the -0 [zero] flag\n");
      fprintf(stderr,"Exiting...\n");
      exit(9);
    RERAISE;
    END_TRY;

#ifdef MEMUSAGE
    Outbuffer_print_result(outbuffer,result,request,noutput+1);
#else
    Outbuffer_print_result(outbuffer,result,request);
#endif
    Result_free(&result);
    Request_free(&request);
    noutput++;

  }

  Stopwatch_free(&worker_stopwatch);
  Diagpool_free(&diagpool);
  Pairpool_free(&pairpool);
  Matchpool_free(&matchpool);
  Dynprog_free(&dynprogR);
  Dynprog_free(&dynprogM);
  Dynprog_free(&dynprogL);
  Oligoindex_free_array(&oligoindices_minor,noligoindices_minor);
  Oligoindex_free_array(&oligoindices_major,noligoindices_major);

  return;
}


static void *
worker_thread (void *data) {
  Oligoindex_T *oligoindices_major, *oligoindices_minor;
  int noligoindices_major, noligoindices_minor;
  Dynprog_T dynprogL, dynprogM, dynprogR;
  Matchpool_T matchpool;
  Pairpool_T pairpool;
  Diagpool_T diagpool;
  Stopwatch_T worker_stopwatch;
  Request_T request;
  Result_T result;
  Sequence_T queryseq;
  int jobid = 0;

  /* Thread-specific data and storage */
  oligoindices_major = Oligoindex_new_major(&noligoindices_major);
  oligoindices_minor = Oligoindex_new_minor(&noligoindices_minor);
  dynprogL = Dynprog_new(nullgap,EXTRAQUERYGAP,maxpeelback,extramaterial_end,extramaterial_paired);
  dynprogM = Dynprog_new(nullgap,EXTRAQUERYGAP,maxpeelback,extramaterial_end,extramaterial_paired);
  dynprogR = Dynprog_new(nullgap,EXTRAQUERYGAP,maxpeelback,extramaterial_end,extramaterial_paired);
  matchpool = Matchpool_new();
  pairpool = Pairpool_new();
  diagpool = Diagpool_new();
  worker_stopwatch = diagnosticp == true ? Stopwatch_new() : (Stopwatch_T) NULL;

  Except_stack_create();

  while ((request = Inbuffer_get_request(inbuffer)) != NULL) {
    if (jobid % POOL_FREE_INTERVAL == 0) {
      Pairpool_free_memory(pairpool);
      Diagpool_free_memory(diagpool);
      Matchpool_free_memory(matchpool);
    }

    TRY
      result = process_request(request,matchpool,pairpool,diagpool,
			       oligoindices_major,noligoindices_major,
			       oligoindices_minor,noligoindices_minor,
			       dynprogL,dynprogM,dynprogR,worker_stopwatch);
    ELSE
      queryseq = Request_queryseq(request);
      if (queryseq == NULL) {
	fprintf(stderr,"NULL");
      } else if (Sequence_accession(queryseq) == NULL) {
	fprintf(stderr,"unnamed (%d bp)",Sequence_fulllength_given(queryseq));
      } else {
	fprintf(stderr,"%s (%d bp)",Sequence_accession(queryseq),Sequence_fulllength_given(queryseq));
      }
      fprintf(stderr,"\n");
      fprintf(stderr,"To obtain a core dump, re-run program on problem sequence with the -0 [zero] flag\n");

      fprintf(stderr,"Exiting...\n");
      exit(9);

    RERAISE;
    END_TRY;

    Outbuffer_put_result(outbuffer,result,request);
    /* Don't free result or request; done by outbuffer thread */
  }

  Except_stack_destroy();

  Stopwatch_free(&worker_stopwatch);
  Diagpool_free(&diagpool);
  Pairpool_free(&pairpool);
  Matchpool_free(&matchpool);
  Dynprog_free(&dynprogR);
  Dynprog_free(&dynprogM);
  Dynprog_free(&dynprogL);
  Oligoindex_free_array(&oligoindices_minor,noligoindices_minor);
  Oligoindex_free_array(&oligoindices_major,noligoindices_major);

  return (void *) NULL;
}


#if 0

static void
align_relative (FILE *input, char **files, int nfiles, int nextchar,
		Sequence_T queryseq, Sequence_T referenceseq) {
  Oligoindex_T *oligoindices_major, *oligoindices_minor;
  int noligoindices_major, noligoindices_minor;
  Diagnostic_T diagnostic;
  bool lowidentityp;
#ifndef PMAP
  bool poorp, repetitivep;
#endif
  Dynprog_T dynprogL, dynprogM, dynprogR;
  Matchpool_T matchpool;
  Pairpool_T pairpool;
  Diagpool_T diagpool;
  Stopwatch_T stopwatch;

  Genomicpos_T genomicstart, genomiclength;
  Sequence_T genomicseg, queryuc, referenceuc;
  int jobid = 0;

  Chimera_T chimera = NULL;
  List_T gregions, stage3list;
  Stage3_T *stage3array, stage3, stage3ref;
  int npaths, i;

  oligoindices_major = Oligoindex_new_major(&noligoindices_major);
  oligoindices_minor = Oligoindex_new_minor(&noligoindices_minor);
  dynprogL = Dynprog_new(nullgap,EXTRAQUERYGAP,maxpeelback,extramaterial_end,extramaterial_paired);
  dynprogM = Dynprog_new(nullgap,EXTRAQUERYGAP,maxpeelback,extramaterial_end,extramaterial_paired);
  dynprogR = Dynprog_new(nullgap,EXTRAQUERYGAP,maxpeelback,extramaterial_end,extramaterial_paired);
  matchpool = Matchpool_new();
  pairpool = Pairpool_new();
  diagpool = Diagpool_new();
  stopwatch = diagnosticp == true ? Stopwatch_new() : (Stopwatch_T) NULL;

  Matchpool_reset(matchpool);
  Pairpool_reset(pairpool);
  Diagpool_reset(diagpool);

  referenceuc = Sequence_uppercase(referenceseq);

  /* Do not trim the mutation refseq */
  diagnostic = Diagnostic_new();
  Oligoindex_set_inquery(&diagnostic->query_badoligos,&diagnostic->query_repoligos,&diagnostic->query_trimoligos,
			 &diagnostic->query_trim_start,&diagnostic->query_trim_end,oligoindices_major[0],
			 Sequence_fullpointer(referenceuc),Sequence_fulllength(referenceuc),/*trimp*/false);
#ifndef PMAP
#if 0
  /* Don't do Sequence_trim, because it affects sequences like NM_018406 */
  Sequence_trim(referenceseq,diagnostic->query_trim_start,diagnostic->query_trim_end);
#endif
#endif
  gregions = Stage1_compute(&lowidentityp,referenceuc,
#ifdef PMAP
			    indexdb_fwd,indexdb_rev,
#else
			    indexdb,
#endif
			    indexdb_size_threshold,chromosome_iit,chrsubset,matchpool,
			    maxintronlen_bound,maxtotallen_bound,min_extra_end,
			    stutterhits,diagnostic,/*stopwatch*/NULL);
  stage3list = apply_stage3(&chimera,gregions,lowidentityp,referenceseq,referenceuc,/*usersegment*/NULL,
			    oligoindices_major,noligoindices_major,oligoindices_minor,noligoindices_minor,
			    matchpool,pairpool,diagpool,dynprogL,dynprogM,dynprogR,stopwatch);
  if (stage3list == NULL) {
    npaths = 0;
    stage3array = (Stage3_T *) NULL;
  } else {
    stage3array = stage3array_from_list(&npaths,stage3list,/*chimerap*/false,/*remove_overlaps_p*/true);
  }

  Diagnostic_free(&diagnostic);

  /* chimera should be NULL */
  for (i = 1; i < npaths; i++) {
    stage3 = stage3array[i];
    Stage3_free(&stage3,/*free_pairarray_p*/true);
  }
  if (npaths > 0) {
    stage3ref = stage3array[0];
#ifdef PMAP
    Stage3_translate_cdna(stage3ref,queryseq,strictp);
    Stage3_backtranslate_cdna(stage3ref,/*diagnosticp*/false);
#else
    Stage3_translate_genomic(stage3ref,/*fulllengthp*/true,/*cds_startpos*/-1,
			     Sequence_fulllength_given(queryseq),/*truncatep*/false,strictp);
#endif
    FREE(stage3array);

    Stage3_genomicbounds(&genomicstart,&genomiclength,stage3ref);
    if (genomealt != NULL) {
      genomicseg = Genome_get_segment(genomealt,genomicstart,genomiclength,chromosome_iit,/*revcomp*/false);
    } else {
      genomicseg = Genome_get_segment(genome,genomicstart,genomiclength,chromosome_iit,/*revcomp*/false);
    }

    while (jobid == 0 || (queryseq = Sequence_read_multifile(&nextchar,&input,&files,&nfiles,maponlyp)) != NULL) {
      Matchpool_reset(matchpool);
      Pairpool_reset(pairpool);
      Diagpool_reset(diagpool);

      fprintf(fp,">");
      Sequence_print_header(stdout,queryseq,checksump);
      diagnostic = Diagnostic_new();
      if (Sequence_fulllength_given(queryseq) <= 0) {
	print_npaths(fp,0,diagnostic,/*usersegment*/NULL,chrsubset,/*chimera*/NULL,EMPTY_SEQUENCE);

      } else if (Sequence_fulllength_given(queryseq) <
#ifdef PMAP
		 index1part_aa
#else
		 index1part
#endif
		 ) {
	print_npaths(fp,0,diagnostic,/*usersegment*/NULL,chrsubset,/*chimera*/NULL,SHORT_SEQUENCE);

      } else {

	queryuc = Sequence_uppercase(queryseq);
#ifdef PMAP
	Oligoindex_set_inquery(&diagnostic->query_badoligos,&diagnostic->query_repoligos,
			       &diagnostic->query_trimoligos,&diagnostic->query_trim_start,
			       &diagnostic->query_trim_end,oligoindices_major[0],Sequence_fullpointer(queryuc),
			       Sequence_fulllength(queryuc),/*trimp*/false);
#else
	diagnostic->query_oligodepth = 
	  Oligoindex_set_inquery(&diagnostic->query_badoligos,&diagnostic->query_repoligos,
				 &diagnostic->query_trimoligos,&diagnostic->query_trim_start,
				 &diagnostic->query_trim_end,oligoindices_major[0],Sequence_fullpointer(queryuc),
				 Sequence_fulllength(queryuc),/*trimp*/true);

	if (diagnostic->query_trimoligos == 0) {
	  poorp = true;
	} else if (((double) diagnostic->query_badoligos/(double) diagnostic->query_trimoligos > MAX_BADOLIGOS) ||
		   (diagnostic->query_trim_end - diagnostic->query_trim_start < 80 && diagnostic->query_badoligos > 0)) {
	  poorp = true;
	} else {
	  poorp = false;
	}
#if 0
	if (diagnostic->query_trimoligos == 0) {
	  repetitivep = false;
	} else if (diagnostic->query_oligodepth > MAX_OLIGODEPTH ||
		   (double) diagnostic->query_repoligos/(double) diagnostic->query_trimoligos > MAX_REPOLIGOS) {
	  repetitivep = true;
	} else {
	  repetitivep = false;
	}
#endif
	repetitivep = false;

	if (poorp == true && prune_poor_p == true) {
	  print_npaths(fp,0,diagnostic,/*usersegment*/NULL,chrsubset,/*chimera*/NULL,POOR_SEQUENCE);
	} else if (repetitivep == true && prune_repetitive_p == true) {
	  print_npaths(fp,0,diagnostic,/*usersegment*/NULL,chrsubset,/*chimera*/NULL,REPETITIVE);
	} else {
#endif /* PMAP */
	  stage3array = stage3_from_usersegment(&npaths,lowidentityp,queryseq,queryuc,genomicseg,
						oligoindices_major,noligoindices_major,oligoindices_minor,noligoindices_minor,
						pairpool,diagpool,dynprogL,dynprogM,dynprogR,stopwatch);
	  
	  if (npaths == 0) {
	    print_npaths(fp,0,diagnostic,/*usersegment*/NULL,chrsubset,/*chimera*/NULL,NO_FAILURE);
	  } else if (printtype == COORDS) {
	    Stage3_print_coordinates(stage3array[0],queryseq,chromosome_iit,
				     invertmode,Sequence_fulllength_given(queryseq));

	  } else {
	    /* Usual output */
	    print_npaths(fp,1,diagnostic,/*usersegment*/NULL,chrsubset,/*chimera*/NULL,NO_FAILURE);
#ifndef PMAP
	    Stage3_translate_cdna_via_reference(stage3array[0],stage3ref,literalrefp);
#endif
	    Stage3_fix_cdna_direction(stage3array[0],stage3ref);
	    Stage3_print_mutations(stage3array[0],stage3ref,chromosome_iit,queryseq,
				   dbversion,printtype,diagnosticp,proteinmode,
				   invertmode,nointronlenp,wraplength,/*maxmutations*/1000000);
	    for (i = 0; i < npaths; i++) {
	      stage3 = stage3array[i];
	      Stage3_free(&stage3,/*free_pairarray_p*/true);
	    }
	    FREE(stage3array);

	  }

#ifndef PMAP
	}
#endif

	Oligoindex_clear_inquery(oligoindices_major[0]);

	Sequence_free(&queryuc);
      }
      Sequence_free(&queryseq);
      jobid++;
    }
    Sequence_free(&genomicseg);
    Stage3_free(&stage3ref);
  }

  Stopwatch_free(&stopwatch);
  Diagpool_free(&diagpool);
  Pairpool_free(&pairpool);
  Dynprog_free(&dynprogR);
  Dynprog_free(&dynprogM);
  Dynprog_free(&dynprogL);
  Oligoindex_free_array(&oligoindices_minor,noligoindices_minor);
  Oligoindex_free_array(&oligoindices_major,noligoindices_major);

  return;
}

#endif

void
check_map_iit (IIT_T map_iit, IIT_T chromosome_iit) {
  char *typestring, *lookup, *p;
  int type, destranded_len;
  bool errorp = false;

  for (type = 1; type < IIT_ntypes(map_iit); type++) {
    lookup = typestring = IIT_typestring(map_iit,type);
    if ((p = rindex(typestring,'+')) != NULL) {
      destranded_len = (p - typestring)/sizeof(char);
      lookup = (char *) CALLOC(destranded_len+1,sizeof(char));
      strncpy(lookup,typestring,destranded_len);
    } else if ((p = rindex(typestring,'-')) != NULL) {
      destranded_len = (p - typestring)/sizeof(char);
      lookup = (char *) CALLOC(destranded_len+1,sizeof(char));
      strncpy(lookup,typestring,destranded_len);
    }

    if (IIT_find_one(chromosome_iit,lookup) < 0) {
      if (p != NULL) {
	fprintf(stderr,"Warning: In %s, type %s (without the %s) does not correspond to a known chromosome in %s.\n",
		map_iitfile,typestring,p,dbversion);
      } else {
	fprintf(stderr,"Warning: In %s, type %s does not correspond to a known chromosome in %s.\n",
		map_iitfile,typestring,dbversion);
      }
      errorp = true;
    }

    if (p != NULL) {
      FREE(lookup);
    }
  }
  if (errorp == true) {
    fprintf(stderr,"Known chromosomes: ");
    IIT_dump_labels(stderr,chromosome_iit);
  }
  return;
}


void
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
  Sequence_T referenceseq = NULL;
  char *genomesubdir = NULL, *mapdir = NULL, *iitfile = NULL, *fileroot = NULL, *p;
  FILE *input = NULL;

  int user_ngap = -1;
  bool showcontigp = true, multiple_sequences_p = false;
  char **files;
  int nfiles;
  unsigned int nread;
  double runtime;

#ifdef HAVE_PTHREAD
  int ret, i;
  pthread_attr_t thread_attr_join;
#ifdef WORKER_DETACH
  pthread_attr_t thread_attr_detach;
#endif
#endif

  int opt, len, c;
  extern int optind;
  extern char *optarg;
  int long_option_index = 0;
  const char *long_name;
  char **argstart;

#ifdef HAVE_SIGACTION
  struct sigaction signal_action;
#endif


  fprintf(stderr,"GMAP version %s called with args:",PACKAGE_VERSION);
  argstart = &(argv[-optind]);
  for (c = 1; c < argc + optind; c++) {
    fprintf(stderr," %s",argstart[c]);
  }
  fprintf(stderr,"\n");


  while ((opt = getopt_long(argc,argv,
#ifdef PMAP
			    "q:D:d:k:Gg:2C:B:K:w:L:x:1t:s:c:H:SA03468:9n:f:ZO5o:V:M:m:ebu:E:PQYNI:i:l:",
#else
			    "q:D:d:k:Gg:2C:B:K:w:L:x:1t:s:c:H:p:SA03468:9n:f:ZO5o:V:M:m:ebu:E:PQFa:Tz:j:YNI:i:l:",
#endif
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

      } else if (!strcmp(long_name,"cmdline")) {
	user_cmdline = optarg; break;

      } else if (!strcmp(long_name,"suboptimal-score")) {
	suboptimal_score = atoi(check_valid_int(optarg));

      } else if (!strcmp(long_name,"splicingdir")) {
	user_splicingdir = optarg;

      } else if (!strcmp(long_name,"nosplicing")) {
	novelsplicingp = false;

      } else if (!strcmp(long_name,"min-intronlength")) {
	min_intronlength = atoi(check_valid_int(optarg));

      } else if (!strcmp(long_name,"allow-close-indels")) {
	if (!strcmp(optarg,"0")) {
	  /* Disallow */
	  close_indels_mode = -1;
	  extraband_single = 0;
	} else if (!strcmp(optarg,"1")) {
	  /* Always allow */
	  close_indels_mode = +1;
	  extraband_single = 3;
	} else if (!strcmp(optarg,"2")) {
	  /* Allow for high-quality alignments */
	  close_indels_mode = 0;
	  extraband_single = 3;
	} else {
	  fprintf(stderr,"allow-close-indels argument %s not recognized.  Only allow 0, 1, or 2.  Run 'gsnap --help' for more information.\n",optarg);
	  exit(9);
	}
      } else if (!strcmp(long_name,"microexon-spliceprob")) {
	microexon_spliceprob = atof(check_valid_float(optarg));
      } else if (!strcmp(long_name,"stage2-start")) {
	suboptimal_score_start = atoi(check_valid_int(optarg));
      } else if (!strcmp(long_name,"stage2-end")) {
	suboptimal_score_end = atoi(check_valid_int(optarg));

      } else if (!strcmp(long_name,"canonical-mode")) {
	if (!strcmp(optarg,"0")) {
	  canonical_mode = 0;
	} else if (!strcmp(optarg,"1")) {
	  canonical_mode = 1;
	} else if (!strcmp(optarg,"2")) {
	  canonical_mode = 2;
	} else {
	  fprintf(stderr,"Canonical level %s not recognized.\n",optarg);
	  fprintf(stderr,"0=low reward for canonical introns, 1=high reward for canonical introns (default)\n");
	  fprintf(stderr,"2=low reward for high-identity seqs, high reward otherwise\n");
	  exit(9);
	}

      } else if (!strcmp(long_name,"cross-species")) {
	use_shifted_canonical_p = true;

      } else if (!strcmp(long_name,"input-buffer-size")) {
	inbuffer_nspaces = atoi(check_valid_int(optarg));
      } else if (!strcmp(long_name,"output-buffer-size")) {
	output_buffer_size = atoi(check_valid_int(optarg));
      } else if (!strcmp(long_name,"print-comment")) {
	print_comment_p = true;
      } else if (!strcmp(long_name,"failsonly")) {
	if (nofailsp == true) {
	  fprintf(stderr,"Cannot specify both --nofails and --failsonly\n");
	  exit(9);
	} else {
	  failsonlyp = true;
	}
      } else if (!strcmp(long_name,"quiet-if-excessive")) {
	quiet_if_excessive_p = true;
      } else if (!strcmp(long_name,"nofails")) {
	if (failsonlyp == true) {
	  fprintf(stderr,"Cannot specify both --nofails and --failsonly\n");
	  exit(9);
	} else {
	  nofailsp = true;
	}
      } else if (!strcmp(long_name,"split-output")) {
	sevenway_root = optarg; break;
#ifndef PMAP
      } else if (!strcmp(long_name,"no-sam-headers")) {
	sam_headers_p = false;
      } else if (!strcmp(long_name,"quality-protocol")) {
	if (user_quality_shift == true) {
	  fprintf(stderr,"Cannot specify both -j (--quality-print-shift) and --quality-protocol\n");
	  exit(9);
	} else if (!strcmp(optarg,"illumina")) {
	  quality_shift = -31;
	  user_quality_shift = true;
	} else if (!strcmp(optarg,"sanger")) {
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
#endif
      } else {
	/* Shouldn't reach here */
	fprintf(stderr,"Don't recognize option %s.  For usage, run 'gsnap --help'",long_name);
	exit(9);
      }
      break;

    case 'q': parse_part(&part_modulus,&part_interval,optarg); break;
    case 'D': user_genomedir = optarg; break;
    case 'd': 
      dbroot = (char *) CALLOC(strlen(optarg)+1,sizeof(char));
      strcpy(dbroot,optarg);
      break;
#ifdef PMAP
    case 'k': index1part_aa = atoi(check_valid_int(optarg));
      if (index1part_aa == 5) {
	/* Okay */
      } else if (index1part_aa == 6) {
	/* Okay */
      } else if (index1part_aa == 7) {
	/* Okay */
      } else if (index1part_aa == 8) {
	/* Okay */
      } else {
	fprintf(stderr,"The only values allowed for -k are 5, 6, 7, or 8\n");
	exit(9);
      }
      break;
#else
    case 'k':
      required_index1part = atoi(check_valid_int(optarg));
      if (required_index1part >= 12 && required_index1part <= 15) {
	/* Okay */
      } else {
	fprintf(stderr,"The only values allowed for -k are 12, 13, 14, or 15\n");
	exit(9);
      }
      break;
#endif
    case 'G': uncompressedp = true; break;
    case 'g': user_genomicseg = optarg; break;
    case '2': user_pairalign_p = true; break;
    case 'C': user_chrsubsetfile = optarg; break;

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

    case 'K': maxintronlen_bound = atoi(check_valid_int(optarg)); break;
    case 'w': shortsplicedist = strtoul(check_valid_int(optarg),NULL,10); break;

    case 'L': maxtotallen_bound = atoi(check_valid_int(optarg)); break;
    case 'x': 
#ifdef PMAP
      chimera_margin = atoi(check_valid_int(optarg))/3; 
#else
      chimera_margin = atoi(check_valid_int(optarg)); 
#endif
      if (chimera_margin < CHIMERA_SLOP) {
	chimera_margin = CHIMERA_SLOP;
      }
      break;
    case '1': maponlyp = true; break;
      /* case 'w': referencefile = optarg; break; */
#ifdef HAVE_PTHREAD
    case 't': nworkers = atoi(check_valid_int(optarg)); break;
#endif
    case 's': splicing_file = optarg; knownsplicingp = true; break;
    case 'c': user_chrsubsetname = optarg; break;
    case 'H': minendexon = atoi(check_valid_int(optarg)); break;

#ifndef PMAP
    case 'p': switch (atoi(check_valid_int(optarg))) {
      case 0: prune_poor_p = false, prune_repetitive_p = false; break;
      case 1: prune_poor_p = true; prune_repetitive_p = false; break;
      case 2: prune_poor_p = false; prune_repetitive_p = true; break;
      case 3: prune_poor_p = true; prune_repetitive_p = true; break;
      default: fprintf(stderr,"Prune level %s not recognized.\n",optarg);
	fprintf(stderr,"0=no pruning, 1=poor seqs, 2=repetitive seqs, 3=both poor and repetitive seqs (default)\n");
	exit(9);
      }
      break;
#endif

    case 'S': printtype = SUMMARY; break;
    case 'A': printtype = ALIGNMENT; break;
    case '0': exception_raise_p = false; break; /* Allows signals to pass through */
    case '3': printtype = CONTINUOUS; break;
    case '4': printtype = CONTINUOUS_BY_EXON; break;
    case '6': debug_graphic_p = true; diagnosticp = false; break;
    case '8': /* diagnosticp = true; */
      if (!strcmp(optarg,"stage1")) {
	stage1debug = true;
      } else if (!strcmp(optarg,"diag")) {
	diag_debug = true;
      } else if (!strcmp(optarg,"stage2")) {
	stage3debug = POST_STAGE2;
      } else if (!strcmp(optarg,"smoothing")) {
	stage3debug = POST_SMOOTHING;
      } else if (!strcmp(optarg,"singles")) {
	stage3debug = POST_SINGLES;
      } else if (!strcmp(optarg,"introns")) {
	stage3debug = POST_INTRONS;
      } else if (!strcmp(optarg,"hmm")) {
	stage3debug = POST_HMM;
      } else if (!strcmp(optarg,"dualbreaks")) {
	stage3debug = POST_DUAL_BREAKS;
      } else if (!strcmp(optarg,"cycles")) {
	stage3debug = POST_CYCLES;
      } else if (!strcmp(optarg,"canonical")) {
	stage3debug = POST_CANONICAL;
      } else if (!strcmp(optarg,"changepoint")) {
	stage3debug = POST_CHANGEPOINT;
      } else if (!strcmp(optarg,"distalmedial")) {
	stage3debug = POST_DISTAL_MEDIAL;
      } else {
	fprintf(stderr,"Allowed arguments for -8 flag are stage2, smoothing, singles, introns, hmm, dualbreaks, cycles, canonical, changepoint, distalmedial\n");
	exit(9);
      }
      break;
    case '9': checkp = true; diagnosticp = true; break;
    case 'n': maxpaths = atoi(check_valid_int(optarg)); break;
    case 'f':
      if (!strcmp(optarg,"1") || !strcmp(optarg,"psl_nt")) {
	printtype = PSL_NT;
#ifdef PMAP
      } else if (!strcmp(optarg,"0") || !strcmp(optarg,"psl_pro")) {
	printtype = PSL_PRO;
#else
      } else if (!strcmp(optarg,"psl")) {
	printtype = PSL_NT;
      } else if (!strcmp(optarg,"6") || !strcmp(optarg,"splicesites")) {
	printtype = SPLICESITES;
      } else if (!strcmp(optarg,"introns")) {
	printtype = INTRONS;
      } else if (!strcmp(optarg,"samse")) {
	printtype = SAM;
	sam_paired_p = false;
      } else if (!strcmp(optarg,"sampe")) {
	printtype = SAM;
	sam_paired_p = true;
#endif
      } else if (!strcmp(optarg,"2") || !strcmp(optarg,"gff3_gene")) {
	printtype = GFF3_GENE;
      } else if (!strcmp(optarg,"3") || !strcmp(optarg,"gff3_match_cdna")) {
	printtype = GFF3_MATCH_CDNA;
      } else if (!strcmp(optarg,"4") || !strcmp(optarg,"gff3_match_est")) {
	printtype = GFF3_MATCH_EST;
      } else if (!strcmp(optarg,"7") || !strcmp(optarg,"map_exons")) {
	printtype = MAP_EXONS;
      } else if (!strcmp(optarg,"8") || !strcmp(optarg,"map_genes")) {
	printtype = MAP_GENES;
      } else if (!strcmp(optarg,"9") || !strcmp(optarg,"coords")) {
	printtype = COORDS;
      } else {
	fprintf(stderr,"Output format %s not recognized\n",optarg);
	exit(9);
      }
      break;
    case 'Z': printtype = COMPRESSED; break;
    case 'O': orderedp = true; break;
    case '5': checksump = true; break;
    case 'o': chimera_overlap = atoi(check_valid_int(optarg)); break;
    case 'V': snps_root = optarg; break;

    case 'M': user_mapdir = optarg; break;
    case 'm': 
      map_iitfile = (char *) CALLOC(strlen(optarg)+1,sizeof(char));
      strcpy(map_iitfile,optarg);
      if ((len = strlen(map_iitfile)) > 4 && strcmp(&(map_iitfile[len-4]),".iit") == 0) {
	map_iitfile[len-4] = '\0';
      }
      break;

    case 'e': map_exons_p = true; break;
    case 'b': map_bothstrands_p = true; break;
    case 'u': nflanking = atoi(check_valid_int(optarg)); break;

    case 'E': 
      if (!strcmp(optarg,"cdna")) {
	printtype = EXONS_CDNA;
      } else if (!strncmp(optarg,"genomic",strlen("genomic"))) {
	printtype = EXONS_GENOMIC;
      } else {
	fprintf(stderr,"Argument to -E flag must be either \"cdna\" or \"genomic\"\n");
	exit(9);
      }
      break;

#ifdef PMAP
    case 'P': printtype = PROTEIN_GENOMIC; break;
    case 'Q': printtype = CDNA; break; 
#else
    case 'P': printtype = CDNA; break;
    case 'Q': printtype = PROTEIN_GENOMIC; break;
    case 'F': fulllengthp = true; break;
    case 'a': cds_startpos = atoi(check_valid_int(optarg)); break;
    case 'T': truncatep = true; fulllengthp = true; break;
    case 'z':
      if (!strcmp(optarg,"sense_force")) {
	sense_try = +1;
	sense_filter = 0;
      } else if (!strcmp(optarg,"antisense_force")) {
	sense_try = -1;
	sense_filter = 0;
      } else if (!strcmp(optarg,"sense_filter")) {
	sense_try = 0;
	sense_filter = +1;
      } else if (!strcmp(optarg,"antisense_filter")) {
	sense_try = 0;
	sense_filter = -1;
      } else if (!strcmp(optarg,"auto")) {
	sense_try = 0;
	sense_filter = 0;
      } else {
	fprintf(stderr,"direction %s not recognized.  Must be sense_force, antisense_force, sense_filter, antisense_filter, or auto\n",optarg);
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

#endif
    case 'Y': strictp = false; break;
    case 'N': nointronlenp = true; break;
    case 'I': invertmode = atoi(check_valid_int(optarg)); break;
    case 'i': user_ngap = atoi(check_valid_int(optarg)); break;
    case 'l': wraplength = atoi(check_valid_int(optarg)); break;

    case '?': fprintf(stderr,"For usage, run 'gmap --help'\n"); exit(9);
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


  if (printtype == SPLICESITES || printtype == INTRONS) {
    if (maxpaths > 1 || (sense_try != +1 && sense_filter != +1)) {
      fprintf(stderr,"For splicesites or introns output, you should probably add flags '-n 1' and either '-z sense_force' or '-z sense_filter'.\n");
    }
  }

  if (user_ngap >= 0) {
    ngap = user_ngap;
  } else if (printtype == EXONS_CDNA || printtype == EXONS_GENOMIC) {
    /* If user didn't specify, then set to zero */
    ngap = 0;
  };

  if (maxintronlen_bound > maxtotallen_bound) {
    maxintronlen_bound = maxtotallen_bound;
  }

#ifdef HAVE_PTHREAD
#ifdef USE_DIAGPOOL
  if (diag_debug == true && nworkers > 0) {
    fprintf(stderr,"For diag output, must specify 0 threads\n");
    exit(9);
  }
#endif
#endif


  /* Handle "?" command-line queries */

  if (user_cmdline != NULL) {
    nchrs = 1;
  } else if (user_pairalign_p == true) {
    nchrs = 1;
  } else if (user_genomicseg != NULL) {
    /* Ignore -D and -d flags */
    nchrs = 1;
  } else if (dbroot == NULL) {
    fprintf(stderr,"Need to specify the -d, -g, -2, or --cmdline flag\n");
    print_program_usage();
    exit(9);
  } else if (!strcmp(dbroot,"?")) {
    Datadir_avail_gmap_databases(stdout,user_genomedir);
    exit(0);
  } else {
    genomesubdir = Datadir_find_genomesubdir(&fileroot,&dbversion,user_genomedir,dbroot);

    iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			      strlen(fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
    sprintf(iitfile,"%s/%s.chromosome.iit",genomesubdir,fileroot);
    if ((chromosome_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
				   /*divstring*/NULL,/*add_iit_p*/false,/*labels_read_p*/true)) == NULL) {
      fprintf(stderr,"IIT file %s is not valid\n",iitfile);
      exit(9);
    } else {
      nchrs = IIT_total_nintervals(chromosome_iit);
    }
    FREE(iitfile);

    if (map_iitfile == NULL) {
      /* Skip */
    } else if (!strcmp(map_iitfile,"?")) {
      Datadir_avail_maps(stdout,user_mapdir,genomesubdir,fileroot);
      exit(0);
    } else {
      mapdir = Datadir_find_mapdir(user_mapdir,genomesubdir,fileroot);
      iitfile = (char *) CALLOC(strlen(mapdir)+strlen("/")+
				strlen(map_iitfile)+strlen(".iit")+1,sizeof(char));
      sprintf(iitfile,"%s/%s.iit",mapdir,map_iitfile);
      if ((map_iit = IIT_read(iitfile,/*name*/map_iitfile,/*readonlyp*/true,/*divread*/READ_ALL,
			      /*divstring*/NULL,/*add_iit_p*/true,/*labels_read_p*/true)) == NULL) {
	fprintf(stderr,"Map file %s.iit not found in %s.  Available files:\n",map_iitfile,mapdir);
	Datadir_list_directory(stderr,mapdir);
	fprintf(stderr,"Either install file %s or specify a full directory path\n",map_iitfile);
	fprintf(stderr,"using the -M flag to gmap.\n");
	exit(9);
      } else {
	map_divint_crosstable = IIT_divint_crosstable(chromosome_iit,map_iit);
      }

      check_map_iit(map_iit,chromosome_iit);

      FREE(iitfile);
      FREE(mapdir);
      FREE(map_iitfile);
    }

    if (splicing_file != NULL) {
      if (user_splicingdir == NULL) {
	if ((splicing_iit = IIT_read(splicing_file,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
				     /*divstring*/NULL,/*add_iit_p*/false,/*labels_read_p*/true)) != NULL) {
	  fprintf(stderr,"Reading splicing file %s locally...",splicing_file);
	} else {
	  iitfile = (char *) CALLOC(strlen(user_splicingdir)+strlen("/")+strlen(splicing_file)+1,sizeof(char));
	  sprintf(iitfile,"%s/%s",user_splicingdir,splicing_file);
	  if ((splicing_iit = IIT_read(splicing_file,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
				       /*divstring*/NULL,/*add_iit_p*/false,/*labels_read_p*/true)) != NULL) {
	    fprintf(stderr,"Reading splicing file %s locally...",splicing_file);
	    FREE(iitfile);
	  }
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
    }
  }


  /* Read user segment before rest of sequences, because of shared usage of sequence.c */
  if (user_cmdline != NULL) {
    p = user_cmdline;
    while (*p != '\0' && *p != ',') {
      p++;
    }
    if (*p == '\0') {
      fprintf(stderr,"--cmdline requires two strings separated by a comma");
      exit(9);
    } else {
      usersegment = Sequence_genomic_new(user_cmdline,(int) (p - user_cmdline),/*copyp*/true);
      if ((min_matches = Sequence_fulllength(usersegment)/2) > MIN_MATCHES) {
	min_matches = MIN_MATCHES;
      }
      p++;
    }

  } else if (user_pairalign_p == true) {
    /* Unfortunately, this procedure reads header of queryseq */
    usersegment = Sequence_read_unlimited(stdin);
    if ((min_matches = Sequence_fulllength(usersegment)/2) > MIN_MATCHES) {
      min_matches = MIN_MATCHES;
    }

  } else if (user_genomicseg != NULL) {
    if ((input = FOPEN_READ_TEXT(user_genomicseg)) == NULL) {
      fprintf(stderr,"Can't open file %s\n",user_genomicseg);
      exit(9);
    }
    if ((usersegment = Sequence_read_unlimited(input)) == NULL) {
      fprintf(stderr,"File %s is empty\n",user_genomicseg);
      exit(9);
    }
    if ((min_matches = Sequence_fulllength(usersegment)/2) > MIN_MATCHES) {
      min_matches = MIN_MATCHES;
    }
    fclose(input);

  } else {
    min_matches = MIN_MATCHES;
  }

  /* Read referencefile before rest of sequences, because of shared usage of sequence.c */
  if (referencefile != NULL) {
    if ((input = FOPEN_READ_TEXT(referencefile)) == NULL) {
      fprintf(stderr,"Can't open file %s\n",referencefile);
      exit(9);
    }
    if ((referenceseq = Sequence_read_unlimited(input)) == NULL) {
      fprintf(stderr,"File %s is empty\n",referencefile);
      exit(9);
    }
    fclose(input);
  }


#ifdef MEMUSAGE
  Mem_usage_init();
  nworkers = 0;
  fprintf(stderr,"For memusage, setting to 0 threads\n");
#endif

  if (user_cmdline != NULL) {
    inbuffer = Inbuffer_cmdline(p,strlen(p));
    nread = 1;
  } else {
    /* Open input stream and peek at first char */
    if (user_pairalign_p == true) {
      input = stdin;
      files = (char **) NULL;
      nfiles = 0;
    } else if (argc == 0) {
      input = stdin;
      files = (char **) NULL;
      nfiles = 0;
    } else {
      input = NULL;
      files = argv;
      nfiles = argc;
    }

    /* Read in first batch of sequences */
    inbuffer = Inbuffer_new(/*nextchar*/'\0',input,files,nfiles,maponlyp,
			    inbuffer_nspaces,inbuffer_maxchars,part_interval,part_modulus,
			    /*filter_if_both_p*/false);
    nread = Inbuffer_fill_init(inbuffer);
  }

  if (nread > 1) {
    multiple_sequences_p = true;
#ifdef HAVE_MMAP
    if (offsetscomp_access != USE_ALLOCATE || genome_access != USE_ALLOCATE) {
      fprintf(stderr,"Note: >1 sequence detected, so index files are being memory mapped.\n");
      fprintf(stderr,"  GMAP can run slowly at first while the computer starts to accumulate\n");
      fprintf(stderr,"  pages from the hard disk into its cache.  To copy index files into RAM\n");
      fprintf(stderr,"  instead of memory mapping, use -B 3, -B 4, or -B 5, if you have enough RAM.\n");
#ifdef HAVE_PTHREAD
      fprintf(stderr,"  For more speed, also try multiple threads (-t <int>), if you have multiple processors or cores.");
#endif
      fprintf(stderr,"\n");
    }
#endif
  } else {
    /* multiple_sequences_p = false; */
    /* fprintf(stderr,"Note: only 1 sequence detected.  Ignoring batch (-B) command\n"); */
    expand_offsets_p = false;
    offsetscomp_access = USE_MMAP_ONLY;
    positions_access = USE_MMAP_ONLY;
    genome_access = USE_MMAP_ONLY;
  }


  /* Prepare genomic data */

  /* Complement_init(); */
  Dynprog_init(nullgap,EXTRAQUERYGAP,maxpeelback,extramaterial_end,extramaterial_paired);
#ifdef PMAP
  Backtranslation_init();
#endif

  if (usersegment != NULL) {
    /* Map against user-provided genomic segment */
    showcontigp = false;
    /* maxpaths = 1; -- no; could have different paths against the user segment. */

    genome = (Genome_T) NULL;
    genomealt = (Genome_T) NULL;
    dbversion = (char *) NULL;
    if (userstage1p == true) {
#ifdef PMAP
      indexdb_fwd = Indexdb_new_segment(Sequence_fullpointer(usersegment),index1part_aa,/*watsonp*/true,index1interval);
      indexdb_rev = Indexdb_new_segment(Sequence_fullpointer(usersegment),index1part_aa,/*watsonp*/false,index1interval);
#else
      indexdb = Indexdb_new_segment(Sequence_fullpointer(usersegment),index1part,index1interval);
#endif
    }
    genome_totallength = Sequence_fulllength(usersegment);
    if (genome_totallength > 1000000) {
      fprintf(stderr,"Genomic sequence is unusually long (%d bp).  GMAP handles genomes better when\n",
	      genome_totallength);
      fprintf(stderr,"  they are converted into gmap databases first using gmap_setup, and then accessed\n");
      fprintf(stderr,"  with the -d flag.\n");
    }

  } else {
    /* Map against genome */
    genome_totallength = IIT_totallength(chromosome_iit);

    chrsubset = Chrsubset_read(user_chrsubsetfile,genomesubdir,fileroot,user_chrsubsetname,
			       chromosome_iit);

    if (showcontigp == true) {
      iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
				strlen(fileroot)+strlen(".contig.iit")+1,sizeof(char));
      sprintf(iitfile,"%s/%s.contig.iit",genomesubdir,fileroot);
      if ((contig_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
				 /*divstring*/NULL,/*add_iit_p*/false,/*labels_read_p*/true)) == NULL) {
	fprintf(stderr,"IIT file %s is not valid\n",iitfile);
	exit(9);
      }
      FREE(iitfile);
    }
  
#ifdef PMAP
    indexdb_fwd = Indexdb_new_genome(&index1part_aa,genomesubdir,fileroot,FWD_FILESUFFIX,/*snps_root*/NULL,
				     required_index1part,/*required_interval*/0,
				     expand_offsets_p,offsetscomp_access,positions_access);
    indexdb_rev = Indexdb_new_genome(&index1part_aa,genomesubdir,fileroot,REV_FILESUFFIX,/*snps_root*/NULL,
				     required_index1part,/*required_interval*/0,
				     expand_offsets_p,offsetscomp_access,positions_access);

#if 0
    indexdb_size_threshold = (int) (10*Indexdb_mean_size(indexdb,/*cmetp*/false,index1part_aa));
#endif
    if (indexdb_fwd == NULL || indexdb_rev == NULL) {
      fprintf(stderr,"Cannot find offsets file %s.%s*offsets or %s.%s*offsets.\n",
	      fileroot,FWD_FILESUFFIX,fileroot,REV_FILESUFFIX);
      fprintf(stderr,"You may need to run 'pmapindex -d %s' to build the indices needed for PMAP.\n",
	      fileroot);
      exit(9);
    }      
#else
    if ((indexdb = Indexdb_new_genome(&index1part,genomesubdir,fileroot,IDX_FILESUFFIX,/*snps_root*/NULL,
				      required_index1part,/*required_interval*/0,
				      expand_offsets_p,offsetscomp_access,positions_access)) == NULL) {
      /* Try older version */
      if ((indexdb = Indexdb_new_genome(&index1part,genomesubdir,fileroot,"id",/*snps_root*/NULL,
					/*required_index1part*/12,/*required_interval*/0,
					expand_offsets_p,offsetscomp_access,positions_access)) == NULL) {
	fprintf(stderr,"Cannot find offsets file %s.%s*offsets or %s.%s*offsets\n",fileroot,IDX_FILESUFFIX,fileroot,"id");
	exit(9);
      }
    }
    indexdb_size_threshold = (int) (10*Indexdb_mean_size(indexdb,/*cmetp*/false,index1part));
    debug(printf("Size threshold is %d\n",indexdb_size_threshold));
#endif
    genome = Genome_new(genomesubdir,fileroot,/*snps_root*/NULL,uncompressedp,genome_access);
    genome_blocks = Genome_blocks(genome);
    if (snps_root == NULL) {
      snp_blocks = (UINT4 *) NULL;
    } else {
      genomealt = Genome_new(genomesubdir,fileroot,snps_root,uncompressedp,genome_access);
      snp_blocks = Genome_blocks(genomealt);
    }

    if (altstrainp == true) {
      if (usersegment != NULL) {
	fprintf(stderr,"Ignoring -s flag when user segment (-g flag) is provided\n");
      } else {
	iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
				  strlen(fileroot)+strlen(".altstrain.iit")+1,sizeof(char));
	sprintf(iitfile,"%s/%s.altstrain.iit",genomesubdir,fileroot);
	if ((altstrain_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
				      /*divstring*/NULL,/*add_iit_p*/false,/*labels_read_p*/false)) == NULL) {
	  fprintf(stderr,"IIT file %s is not valid\n",iitfile);
	  exit(9);
	}
	FREE(iitfile);
      }
    }

    FREE(genomesubdir);
    FREE(fileroot);
    FREE(dbroot);
  }

  if (splicing_file != NULL && genome != NULL) {
    if (Genome_blocks(genome) == NULL) {
      fprintf(stderr,"known splicing can be used only with compressed genome\n");
    } else {
      /* TODO: Handle case for observed distances */
      /* min_extra_end no longer used by gregion.c */
      min_extra_end = shortsplicedist;

      splicing_divint_crosstable = IIT_divint_crosstable(chromosome_iit,splicing_iit);
      if ((donor_typeint = IIT_typeint(splicing_iit,"donor")) >= 0 && 
	  (acceptor_typeint = IIT_typeint(splicing_iit,"acceptor")) >= 0) {
	fprintf(stderr,"found donor and acceptor tags, so treating as splicesites file\n");
	splicesites = Splicetrie_retrieve_via_splicesites(&distances_observed_p,&splicetypes,&splicedists,
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
	splicesites = Splicetrie_retrieve_via_introns(&splicetypes,&splicedists,
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
    }

    fprintf(stderr,"done\n");
  }


  if (genome != NULL) {
    Genome_setup(genome);
    Genome_hr_setup(Genome_blocks(genome),/*snp_blocks*/genomealt ? Genome_blocks(genomealt) : NULL,
		    /*query_unk_mismatch_p*/false,/*genome_unk_mismatch_p*/true,/*mode*/STANDARD);
    Maxent_hr_setup(Genome_blocks(genome));
#ifndef PMAP
    Oligoindex_hr_setup(Genome_blocks(genome));
#endif
  }
#ifdef PMAP
  Indexdb_setup(index1part_aa);
  Stage1_setup(index1part_aa);
#else
  Indexdb_setup(index1part);
  Stage1_setup(index1part);
#endif
  Stage2_setup(/*splicingp*/novelsplicingp == true || knownsplicingp == true,
	       suboptimal_score_start,suboptimal_score_end);
  Dynprog_setup(splicing_iit,splicing_divint_crosstable,donor_typeint,acceptor_typeint,
		splicesites,splicetypes,splicedists,nsplicesites,
		trieoffsets_obs,triecontents_obs,trieoffsets_max,triecontents_max);
  Pair_setup(trim_mismatch_score,trim_indel_score);
  Stage3_setup(/*splicingp*/novelsplicingp == true || knownsplicingp == true,
	       splicing_iit,splicing_divint_crosstable,donor_typeint,acceptor_typeint,
	       splicesites,min_intronlength,max_deletionlength,
	       /*expected_pairlength*/0,/*pairlength_deviation*/0);
  Splicetrie_setup(splicesites,splicefrags_ref,splicefrags_alt,
		   trieoffsets_obs,triecontents_obs,trieoffsets_max,triecontents_max,
		   /*snpp*/false,amb_closest_p,/*amb_clip_p*/true,/*min_shortend*/2);

  /* Setup outbuffer */
#ifndef PMAP
  if (printtype == SAM) {
    if (sam_read_group_id == NULL && sam_read_group_name != NULL) {
      sam_read_group_id = sam_read_group_name;
    } else if (sam_read_group_id != NULL && sam_read_group_name == NULL) {
      sam_read_group_name = sam_read_group_id;
    }
  }
#endif

  outbuffer = Outbuffer_new(output_buffer_size,nread,sevenway_root,
			    /*chimerap*/chimera_margin > 0 ? true : false,
			    user_genomicseg,usersegment,dbversion,genome,chromosome_iit,
			    chrsubset,contig_iit,altstrain_iit,map_iit,
			    map_divint_crosstable,printtype,checksump,chimera_margin,
#ifndef PMAP
			    sam_headers_p,quality_shift,sam_paired_p,
			    sam_read_group_id,sam_read_group_name,
			    sam_read_group_library,sam_read_group_platform,
#endif
			    nofailsp,failsonlyp,fails_as_input_p,maxpaths,quiet_if_excessive_p,
			    map_exons_p,map_bothstrands_p,print_comment_p,nflanking,
			    proteinmode,invertmode,nointronlenp,wraplength,ngap,cds_startpos,
			    fulllengthp,truncatep,strictp,diagnosticp,maponlyp,
			    stage1debug,diag_debug,debug_graphic_p,argc,argv,optind);

  Inbuffer_set_outbuffer(inbuffer,outbuffer);


  stopwatch = Stopwatch_new();
  Stopwatch_start(stopwatch);

  if (referenceseq != NULL) {
    fprintf(stderr,"Relative alignment currently not implemented\n");
    exit(9);
    chimera_margin = -1;
    /* align_relative(input,files,nfiles,nextchar,queryseq,referenceseq); */
    Sequence_free(&referenceseq);

  } else {
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
	fprintf(stderr,"ERROR: pthread_attr_setdetachstate %d\n",ret);
	exit(1);
      }
#endif
      pthread_attr_init(&thread_attr_join);
      if ((ret = pthread_attr_setdetachstate(&thread_attr_join,PTHREAD_CREATE_JOINABLE)) != 0) {
	fprintf(stderr,"ERROR: pthread_attr_setdetachstate %d\n",ret);
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
#ifdef WORKER_DETACH
	pthread_create(&(worker_thread_ids[i]),&thread_attr_detach,worker_thread,(void *) NULL);
#else
	/* Need to have worker threads finish before we call Inbuffer_free() */
	pthread_create(&(worker_thread_ids[i]),&thread_attr_join,worker_thread,(void *) NULL);
#endif
      }
    
      pthread_join(output_thread_id,NULL);
      for (i = 0; i < nworkers; i++) {
	pthread_join(worker_thread_ids[i],NULL);
      }

      /* Do not delete global_except_key, because worker threads might still need it */
      /* Except_term_pthread(); */

      FREE(worker_thread_ids);

    }

#endif /* HAVE_PTHREAD */
  }

  runtime = Stopwatch_stop(stopwatch);
  Stopwatch_free(&stopwatch);

  nread = Outbuffer_nread(outbuffer);
  fprintf(stderr,"Processed %u queries in %.2f seconds (%.2f queries/sec)\n",
	  nread,runtime,(double) nread/runtime);

  Outbuffer_free(&outbuffer);
  Inbuffer_free(&inbuffer);	/* Also closes inputs */

  if (usersegment != NULL) {
    Sequence_free(&usersegment);
  }

#ifdef PMAP
  Backtranslation_term();
#endif
  Dynprog_term();


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

  if (splicing_iit != NULL) {
    FREE(splicing_divint_crosstable);
    IIT_free(&splicing_iit);
  }


#ifdef PMAP
 if (indexdb_rev != NULL) {
    Indexdb_free(&indexdb_rev);
  }
  if (indexdb_fwd != NULL) {
    Indexdb_free(&indexdb_fwd);
  }
#else
  if (indexdb != NULL) {
    Indexdb_free(&indexdb);
  }
#endif
  if (dbversion != NULL) {
    FREE(dbversion);
  }
  if (altstrain_iit != NULL) {
    IIT_free(&altstrain_iit);
  }
  if (genomealt != NULL) {
    Genome_free(&genomealt);
  }
  if (genome != NULL) {
    Genome_free(&genome);
  }
  if (map_iit != NULL) {
    IIT_free(&map_iit);
  }
  if (contig_iit != NULL) {
    IIT_free(&contig_iit);
  }
  if (chromosome_iit != NULL) {
    IIT_free(&chromosome_iit);
  }
  if (chrsubset != NULL) {
    Chrsubset_free(&chrsubset);
  }

  return 0;
}


static void
print_program_usage () {
#ifdef PMAP
    fprintf(stdout,"\
Usage: pmap [OPTIONS...] <FASTA files...>, or\n\
       cat <FASTA files...> | pmap [OPTIONS...]\n\
");
#else
    fprintf(stdout,"\
Usage: gmap [OPTIONS...] <FASTA files...>, or\n\
       cat <FASTA files...> | gmap [OPTIONS...]\n\
");
#endif

    fprintf(stdout,"\n\
Input options (must include -d or -g)\n\
  -D, --dir=directory            Genome directory\n\
  -d, --db=STRING                Genome database.  If argument is '?' (with\n\
                                   the quotes), this command lists available databases.\n\
  -k, --kmer=INT                 kmer size to use in genome database\n\
                                   (allowed values: 12-15).  If not specified, the program will\n\
                                   find the highest available kmer size in the genome database\n\
  -G, --genomefull               Use full genome (all ASCII chars allowed;\n\
                                   built explicitly during setup), not\n\
                                   compressed version\n\
  -g, --gseg=filename            User-supplied genomic segment\n\
  -2, --pairalign                Align two sequences in FASTA format via stdin, first one being\n\
                                   genomic and second one being cDNA\n\
  --cmdline=STRING,STRING        Align these two sequences provided on the command line,\n\
                                   first one being genomic and second one being cDNA\n\
  -q, --part=INT/INT             Process only the i-th out of every n sequences\n\
                                   e.g., 0/100 or 99/100 (useful for distributing jobs\n\
                                   to a computer farm).\n\
  --input-buffer-size=INT        Size of input buffer (program reads this many sequences\n\
                                   at a time for efficiency) (default 1000)\n\
\n\
");

    fprintf(stdout,"Computation options\n");
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
  --nosplicing                   Turns off splicing (useful for aligning genomic sequences\n\
                                   onto a genome)\n\
  --min-intronlength=INT         Min length for one internal intron (default 9).  Below this size,\n\
                                   a genomic gap will be considered a deletion rather than an intron.\n\
  -K, --intronlength=INT         Max length for one internal intron (default 1000000)\n\
  -w, --localsplicedist=INT      Max length for known splice sites at ends of sequence\n\
                                   (default 200000)\n\
  -L, --totallength=INT          Max total intron length (default 2400000)\n\
  -x, --chimera-margin=INT       Amount of unaligned sequence that triggers\n\
                                   search for the remaining sequence (default 40).\n\
                                   Enables alignment of chimeric reads, and may help\n\
                                   with some non-chimeric reads.  To turn off, set to\n\
                                   a large value (greater than the query length).\n\
");

#if 0
    fprintf(stdout,"\
  -w, --reference=filename       Reference cDNA sequence for relative alignment\n\
");
#endif

#ifdef HAVE_PTHREAD
    fprintf(stdout,"\
  -t, --nthreads=INT             Number of worker threads\n\
");
#endif
    fprintf(stdout,"\
  -C, --chrsubsetfile=filename   User-supplied chromosome subset file\n\
  -c, --chrsubset=string         Chromosome subset to search\n\
  -z, --direction=STRING         cDNA direction (sense_force, antisense_force,\n\
                                   sense_filter, antisense_filter,or auto (default))\n\
  -H, --trimendexons=INT         Trim end exons with fewer than given number of matches\n\
                                   (in nt, default 12)\n\
  --cross-species                For cross-species alignments, use a more sensitive search\n\
                                   for canonical splicing\n\
  --canonical-mode=INT           Reward for canonical and semi-canonical introns\n\
                                   0=low reward, 1=high reward (default), 2=low reward for\n\
                                   high-identity sequences and high reward otherwise\n\
  --allow-close-indels=INT       Allow an insertion and deletion close to each other\n\
                                   (0=no, 1=yes (default), 2=only for high-quality alignments)\n\
  --microexon-spliceprob=FLOAT   Allow microexons only if one of the splice site probabilities is\n\
                                   greater than this value (default 0.90)\n\
");

#if 0
    /* Causes seg faults, so do not advertise */
    fprintf(stdout,"\
  -s, --splicing=STRING          Look for splicing involving known sites\n\
                                   (in <STRING>.iit)\n\
");
#endif

#ifndef PMAP
    fprintf(stdout,"\
  -p, --prunelevel               Pruning level: 0=no pruning (default), 1=poor seqs,\n\
                                   2=repetitive seqs, 3=poor and repetitive\n\
");
#endif

    fprintf(stdout,"\n");
    fprintf(stdout,"\
Output types\n\
  -S, --summary                  Show summary of alignments only\n\
  -A, --align                    Show alignments\n\
  -3, --continuous               Show alignment in three continuous lines\n\
  -4, --continuous-by-exon       Show alignment in three lines per exon\n\
  -Z, --compress                 Print output in compressed format\n\
  -E, --exons=STRING             Print exons (\"cdna\" or \"genomic\")\n\
");

#ifdef PMAP    
    fprintf(stdout,"\
  -P, --protein_gen              Print protein sequence (genomic)\n\
  -Q, --nucleotide               Print inferred nucleotide sequence from protein\n\
");
#else
    fprintf(stdout,"\
  -P, --protein_dna              Print protein sequence (cDNA)\n\
  -Q, --protein_gen              Print protein sequence (genomic)\n\
");
#endif

#ifdef PMAP
    fprintf(stdout,"\
  -f, --format=INT               Other format for output (also note the -A and -S options\n\
                                   and other options listed under Output types):\n\
                                   psl_pro (or 0) = PSL format in protein coords,\n\
                                   psl_nt (or 1) = PSL format in nucleotide coords,\n\
                                   gff3_gene (or 2) = GFF3 gene format,\n\
                                   gff3_match_cdna (or 3) = GFF3 cDNA_match format,\n\
                                   gff3_match_est (or 4) = GFF3 EST_match format,\n\
                                   map_exons (or 7) = IIT FASTA exon map format,\n\
                                   map_genes (or 8) = IIT FASTA map format,\n\
                                   coords (or 9) = coords in table format\n\
");
#else
    fprintf(stdout,"\
  -f, --format=INT               Other format for output (also note the -A and -S options\n\
                                   and other options listed under Output types):\n\
                                   psl (or 1) = PSL (BLAT) format,\n\
                                   gff3_gene (or 2) = GFF3 gene format,\n\
                                   gff3_match_cdna (or 3) = GFF3 cDNA_match format,\n\
                                   gff3_match_est (or 4) = GFF3 EST_match format,\n\
                                   splicesites (or 6) = splicesites output (for GSNAP splicing file),\n\
                                   introns = introns output (for GSNAP splicing file),\n\
                                   map_exons (or 7) = IIT FASTA exon map format,\n\
                                   map_genes (or 8) = IIT FASTA map format,\n\
                                   coords (or 9) = coords in table format,\n\
                                   sampe = SAM format (setting paired_read bit in flag),\n\
                                   samse = SAM format (without setting paired_read bit)\n\
");
#endif
    fprintf(stdout,"\n");

    fprintf(stdout,"\
Output options\n\
  -n, --npaths=INT               Maximum number of paths to show.  If set to 0,\n \
                                 prints two paths if chimera detected, else one.\n\
  --quiet-if-excessive           If more than maximum number of paths are found,\n\
                                   then nothing is printed.\n\
  --suboptimal-score=INT         Report only paths whose score is within this value of the\n\
                                   best path.  By default, if this option is not provided,\n\
                                   the program prints all paths found.\n\
  -O, --ordered                  Print output in same order as input (relevant\n\
                                   only if there is more than one worker thread)\n\
  -5, --md5                      Print MD5 checksum for each query sequence\n\
  -o, --chimera-overlap          Overlap to show, if any, at chimera breakpoint\n\
  --failsonly                    Print only failed alignments, those with no results\n\
  --nofails                      Exclude printing of failed alignments\n\
  --fails-as-input=STRING        Print completely failed alignments as input FASTA or FASTQ format\n\
                                   Allowed values: yes, no\n\
  -V, --usesnps=STRING           Use database containing known SNPs (in <STRING>.iit, built\n\
                                   previously using snpindex) for reporting output\n\
");

  fprintf(stdout,"\
  --split-output=STRING          Basename for multiple-file output, separately for nomapping,\n\
                                   uniq, mult, (and chimera, if --chimera-margin is selected)\n\
  --output-buffer-size=INT       Buffer size, in queries, for output thread (default 1000).  When the number\n\
                                   of results to be printed exceeds this size, the worker threads are halted\n\
                                   until the backlog is cleared\n\
");


#ifdef PMAP    
    fprintf(stdout,"\
  -Y, --tolerant                 Translates genome with corrections for frameshifts\n\
");
#else
    fprintf(stdout,"\
  -F, --fulllength               Assume full-length protein, starting with Met\n\
  -a, --cdsstart=INT             Translate codons from given nucleotide (1-based)\n\
  -T, --truncate                 Truncate alignment around full-length protein, Met to Stop\n\
                                 Implies -F flag.\n\
  -Y, --tolerant                 Translates cDNA with corrections for frameshifts\n\
");
#endif

    fprintf(stdout,"\n");

#ifndef PMAP
  fprintf(stdout,"Options for SAM output\n");
  fprintf(stdout,"\
  --no-sam-headers               Do not print headers beginning with '@'\n\
  --read-group-id=STRING         Value to put into read-group id (RG-ID) field\n\
  --read-group-name=STRING       Value to put into read-group name (RG-SM) field\n\
  --read-group-library=STRING    Value to put into read-group library (RG-LB) field\n\
  --read-group-platform=STRING   Value to put into read-group library (RG-PL) field\n\
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
                                 Or you can specify the print shift with this flag:\n\
  -j, --quality-print-shift=INT  Shift FASTQ quality scores by this amount in output\n\
                                   (default is 0 for sanger protocol; to change Illumina input\n\
                                   to Sanger output, select -31)\n\
");
#endif

    fprintf(stdout,"\
External map file options\n\
  -M, --mapdir=directory         Map directory\n\
  -m, --map=iitfile              Map file.  If argument is '?' (with the quotes),\n\
                                   this lists available map files.\n\
  -e, --mapexons                 Map each exon separately\n\
  -b, --mapboth                  Report hits from both strands of genome\n\
  -u, --flanking=INT             Show flanking hits (default 0)\n\
  --print-comment                Show comment line for each hit\n\
");
    fprintf(stdout,"\n");

    fprintf(stdout,"\
Alignment output options\n\
  -N, --nolengths                No intron lengths in alignment\n\
  -I, --invertmode=INT           Mode for alignments to genomic (-) strand:\n\
                                   0=Don't invert the cDNA (default)\n\
                                   1=Invert cDNA and print genomic (-) strand\n\
                                   2=Invert cDNA and print genomic (+) strand\n\
  -i, --introngap=INT            Nucleotides to show on each end of intron (default=3)\n\
  -l, --wraplength=INT           Wrap length for alignment (default=50)\n\
\n\
Help options\n\
  --version                      Show version\n\
  --help                         Show this help message\n\
");

    return;
}
