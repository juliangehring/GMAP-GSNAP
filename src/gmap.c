static char rcsid[] = "$Id: gmap.c,v 1.439 2010/02/03 02:13:41 twu Exp $";
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

#include "sequence.h"
#include "oligoindex.h"
#include "match.h"
#include "matchpool.h"
#include "pairpool.h"
#include "intpool.h"
#include "diagpool.h"
#include "stopwatch.h"
#include "genome.h"
#include "gbuffer.h"
#include "stage1.h"
#include "gregion.h"
#include "stage2.h"
#include "dynprog.h"
#include "stage3.h"
#ifdef BETATEST
#include "chimera.h"
#endif
#ifdef PMAP
#include "backtranslation.h"
#endif
#include "indexdb.h"
#include "result.h"
#include "request.h"
#include "reqpost.h"
#include "blackboard.h"
#include "intlist.h"
#include "list.h"
#include "listdef.h"
#include "iit-read.h"
#include "chrsubset.h"
#include "datadir.h"
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
#define CHIMERA_SLOP 7
#else
#define CHIMERA_SLOP 25
#endif


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

/* Chimera parts */
#ifdef DEBUG3
#define debug3(x) x
#else
#define debug3(x)
#endif


/************************************************************************
 *   Global variables
 ************************************************************************/

static IIT_T chromosome_iit = NULL;
static Genomicpos_T genome_totallength = 0;
static Chrsubset_T chrsubset = NULL;
static IIT_T contig_iit = NULL;
static Genome_T genome = NULL;
static Genome_T genomealt = NULL;

#ifdef PMAP
static Indexdb_T indexdb_fwd = NULL;
static Indexdb_T indexdb_rev = NULL;
#else
static Indexdb_T indexdb = NULL;
#endif
static int indexdb_size_threshold = 0;

static IIT_T altstrain_iit = NULL;

static char *snps_root = (char *) NULL;
static IIT_T map_iit = NULL;
static int *map_divint_crosstable = NULL;

static bool map_iit_universal_p = false;
static int map_iit_forward_type = 0;
static int map_iit_reverse_type = 0;

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

static int stuttercycles = 2;
static int stutterhits = 3;
static int sufflookback = 60;
static int nsufflookback = 5;

#if 0
static int maxoligohits = 400; /* Must be smaller than ALLOC in oligoindex.c */
#endif
static int nullgap = 600;
static int extramaterial_end = 10;
static int extramaterial_paired = 8; /* Should be at least indexsize in nt */
static int extraband_single = 0; /* This is in addition to length2 -
				    length1.  If onesidegap is true in
				    dynprog.c, then this is equivalent
				    to extraband_single of 0 */
static int extraband_end = 6; /* This is only on both sides of main diagonal */
static int extraband_paired = 7; /* This is in addition to length2 - length1 */
static int minendexon = 12;


/************************************************************************
 *   Program options
 ************************************************************************/

/* Input options */
static char *user_genomedir = NULL;
static char *dbroot = NULL;
static char *dbversion = NULL;
static char *user_genomicseg = NULL;
static char *user_chrsubsetfile = NULL;
static int batch_modulus = 0;
static int batch_interval = 1;

/* Compute options */
static bool batch_offsets_p = true;
static bool batch_positions_p = true;
static bool batch_genome_p = false;
static int maxtotallen_bound = 2400000;
static int maxintronlen_bound = 1000000;
static int chimera_margin = -1;
static bool maponlyp = false;
#ifdef PMAP
static bool userstage1p = false; /* Apply stage 1 for user-provided genomic segments.  Must be false. */
#else
static bool userstage1p = false; /* Apply stage 1 for user-provided genomic segments */
#endif
static int index1interval = 3; /* Stage 1 interval if user provides a genomic segment */
static char *referencefile = NULL;
#ifndef PMAP
static bool literalrefp = false;
#endif
static bool altstrainp = false;
#ifdef HAVE_PTHREAD
static pthread_t input_thread_id, output_thread_id, *worker_thread_ids;
static int nworkers = 1;	/* (int) sysconf(_SC_NPROCESSORS_ONLN) */
#endif
#ifndef PMAP
static bool prune_poor_p = false;
static bool prune_repetitive_p = false;
#endif
static bool end_microexons_p = false;
static int canonical_mode = 1;
static char *user_chrsubsetname = NULL;


/* Output options */
static bool summaryonlyp = false;
static bool showalignp = false;
static bool exception_raise_p = true;
static bool continuousp = false;
static bool continuous_by_exon_p = false;
static bool debug_graphic_p = false;
static Stage3debug_T stage3debug = NO_STAGE3DEBUG;
static bool diagnosticp = false;
static bool checkp = false;
static int maxpaths = 5;	/* 0 means 1 if nonchimeric, 2 if chimeric */
static bool psloutput_nt_p = false;
#ifdef PMAP
static bool psloutput_pro_p = false;
#endif
static bool gffoutput_p = false;
static bool gff_gene_format_p = false;
static bool mapoutput_p = false;
static bool exon_mapoutput_p = false;
static bool splicesites_output_p = false;
static bool print_coordinates_p = false;
static bool compressoutp = false;
static bool orderedp = false;
static bool checksump = false;
static int chimera_overlap = 0;

/* Map file options */
static char *user_mapdir = NULL;
static char *map_iitfile = NULL;
static bool map_exons_p = false;
static bool map_bothstrands_p = false;
static int nflanking = 0;

/* Alignment options */
static bool cdna_exons_p = false;
static bool genomic_exons_p = false;
static bool print_cdna_p = false;
static bool protein_genomic_p = false;
static bool fulllengthp = false;
static bool truncatep = false;
static int sense = 0;		/* both */
static bool strictp = true;
static int proteinmode = 1;
static bool uncompressedp = false;
static bool nointronlenp = false;
static int invertmode = 0;
static int ngap = 3;
static int wraplength = 50;

#ifdef PMAP
/* Used alphabetically: 0135789ABbCcDdEefGgHIiKLlMmNnOoPQRSstuVvwXxYZ */
#else
/* Used alphabetically: 0135789ABbCcDdEeFfGgHIiKLlMmNnOoPpQRSsTtuVvwXxYZ */
#endif

static struct option long_options[] = {
  /* Input options */
  {"dir", required_argument, 0, 'D'},	/* user_genomedir */
  {"db", required_argument, 0, 'd'}, /* dbroot */
  {"genomefull", no_argument, 0, 'G'}, /* uncompressedp */
  {"gseg", required_argument, 0, 'g'}, /* user_genomicseg */
  {"jobdiv", required_argument, 0, 'q'}, /* batch_modulus, batch_interval */

  /* Compute options */
#ifdef HAVE_MMAP
  {"batch", required_argument, 0, 'B'}, /* batch_offsets_p, batch_positions_p, batch_genome_p */
#endif
  {"intronlength", required_argument, 0, 'K'}, /* maxintronlen_bound */
  {"totallength", required_argument, 0, 'L'}, /* maxtotallen_bound */
  {"chimera_margin", required_argument, 0, 'x'}, /* chimera_margin */
  {"reference", required_argument, 0, 'w'}, /* referencefile */
#ifdef HAVE_PTHREAD
  {"nthreads", required_argument, 0, 't'}, /* nworkers */
#endif
  {"altstrain", no_argument, 0, 's'},	/* altstrainp */
  {"chrsubsetfile", required_argument, 0, 'C'}, /* user_chrsubsetfile */
  {"chrsubset", required_argument, 0, 'c'}, /* user_chrsubsetname */
  {"trimendexons", required_argument, 0, 'H'}, /* minendexon */
  {"canonical", required_argument, 0, 'X'}, /* canonical_mode */
#ifndef PMAP
  {"prunelevel", required_argument, 0, 'p'}, /* prune_poor_p, prune_repetitive_p */
#endif

  /* Output options */
  {"summary", no_argument, 0, 'S'}, /* summaryonlyp */
  {"align", no_argument, 0, 'A'}, /* showalignp */
  {"continuous", no_argument, 0, '3'}, /* continuousp */
  {"noexceptions", no_argument, 0, '0'}, /* exception_raise_p */
  {"graphic", no_argument, 0, '6'}, /* debug_graphic_p */
  {"stage3debug", required_argument, 0, '8'}, /* stage3debug, diagnosticp */
  {"check", no_argument, 0, '9'}, /* checkp */
  {"npaths", required_argument, 0, 'n'}, /* maxpaths */
  {"format", required_argument, 0, 'f'}, /* psloutput_nt_p, psloutput_pro_p, gffoutput_p, gff_gene_format_p, mapoutput_p, splicesites_output_p, exon_mapoutput_p, print_coordinates_p */
  {"compress", no_argument, 0, 'Z'}, /* compressoutp */
  {"ordered", no_argument, 0, 'O'}, /* orderedp */
  {"md5", no_argument, 0, '5'}, /* checksump */
  {"chimera_overlap", required_argument, 0, 'o'}, /* chimera_overlap */
  {"usesnps", required_argument, 0, 'V'}, /* snp_root */

  /* Map file options */
  {"mapdir", required_argument, 0, 'M'}, /* user_mapdir */
  {"map", required_argument, 0, 'm'},	/* map_iitfile */
  {"mapexons", no_argument, 0, 'e'}, /* map_exons_p */
  {"mapboth", no_argument, 0, 'b'}, /* map_bothstrands_p */
  {"nflanking", required_argument, 0, 'u'}, /* nflanking */

  /* Alignment options */
  {"exons", required_argument, 0, 'E'}, /* cdna_exons_p, genomic_exons_p */
#ifdef PMAP
  {"protein_gen", no_argument, 0, 'P'}, /* protein_genomic_p */
  {"nucleotide", no_argument, 0, 'Q'}, /* print_cdna_p */
#else
  {"protein_dna", no_argument, 0, 'P'}, /* print_cdna_p */
  {"protein_gen", no_argument, 0, 'Q'}, /* protein_genomic_p */
  {"fulllength", no_argument, 0, 'F'}, /* fulllengthp */
  {"truncate", no_argument, 0, 'T'}, /* truncatep */
  {"direction", required_argument, 0, 'z'}, /* sense */
#endif
  {"tolerant", no_argument, 0, 'Y'}, /* strictp */
  {"nolengths", no_argument, 0, 'N'},	/* nointronlenp */
  {"invertmode", required_argument, 0, 'I'}, /* invertmode */
  {"introngap", required_argument, 0, 'i'}, /* ngap */
  {"wraplength", required_argument, 0, 'l'}, /* wraplength */
  
  /* Help options */
  {"version", no_argument, 0, 'v'}, /* print_program_version */
  {"help", no_argument, 0, '?'}, /* print_program_usage */
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
  fprintf(stdout,"Stage 1 index size: %d aa\n",INDEX1PART_AA);
#endif
  fprintf(stdout,"Sizes: off_t (%lu), size_t (%lu), unsigned int (%lu), long int (%lu)\n",
	  sizeof(off_t),sizeof(size_t),sizeof(unsigned int),sizeof(long int));
  fprintf(stdout,"Default gmap directory (compiled): %s\n",GMAPDB);
  genomedir = Datadir_find_genomedir(/*user_genomedir*/NULL);
  fprintf(stdout,"Default gmap directory (environment): %s\n",genomedir);
  FREE(genomedir);
#ifdef PMAP
  fprintf(stdout,"Thomas D. Wu, Genentech, Inc.\n");
#else
  fprintf(stdout,"Thomas D. Wu and Colin K. Watanabe, Genentech, Inc.\n");
#endif
  fprintf(stdout,"Contact: twu@gene.com\n");
  fprintf(stdout,"\n");
  return;
}

/* This flag is not well-supported, and therefore hidden, but
   kept for backwards compatibility */
/*  -R, --rel=STRING               Release\n\ */

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
  -G, --genomefull               Use full genome (all ASCII chars allowed;\n\
                                   built explicitly during setup), not\n\
                                   compressed version\n\
  -g, --gseg=filename            User-suppled genomic segment\n\
  -q, --jobdiv=INT/INT           Process only i out of every n sequences\n\
                                   e.g., 0/100 or 99/100\n\
\n\
");

    fprintf(stdout,"Computation options\n");
#ifdef HAVE_MMAP
    fprintf(stdout,"\
  -B, --batch=INT                Batch mode (0 = no pre-loading, 1 = pre-load only indices;\n\
                                   2 (default) = pre-load both indices and genome)\n\
");
#endif
    fprintf(stdout,"\
  -K, --intronlength=INT         Max length for one intron (default 1000000)\n\
  -L, --totallength=INT          Max total intron length (default 2400000)\n\
  -x, --chimera_margin=INT       Amount of unaligned sequence that triggers\n\
                                   search for a chimera (default off)\n\
  -w, --reference=filename       Reference cDNA sequence for relative alignment\n\
");
#ifdef HAVE_PTHREAD
    fprintf(stdout,"\
  -t, --nthreads=INT             Number of worker threads\n\
");
#endif
    fprintf(stdout,"\
  -s, --altstrain                Search alternate strains in addition\n\
  -C, --chrsubsetfile=filename   User-suppled chromosome subset file\n\
  -c, --chrsubset=string         Chromosome subset to search\n\
  -z, --direction=STRING         cDNA direction (sense, antisense, or auto (default))\n\
  -H, --trimendexons=INT         Trim end exons with fewer than given number of matches\n\
                                   (in nt, default 12)\n\
  -X, --canonical=INT            Reward for canonical and semi-canonical introns\n\
                                   0=low reward, 1=high reward (default), 2=low reward for\n\
                                   high-identity sequences and high reward otherwise\n\
");
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
  -4, --alignedexons             Show alignment in three lines per exon\n\
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
  -f, --format=INT               Format for output (0 = PSL format in protein coords,\n\
                                   1 = PSL format in nucleotide coords, 2 = GFF3 gene format,\n\
                                   3 = GFF3 match format, 7 = IIT FASTA exon map format,\n\
                                   8 = IIT FASTA map format, 9 = coords in table format)\n\
");
#else
    fprintf(stdout,"\
  -f, --format=INT               Format for output (1 = PSL (BLAT) format, 2 = GFF3 gene format,\n\
                                   3 = GFF3 match format, 6 = splicesites output (for GSNAP),\n\
                                   7 = IIT FASTA exon map format, 8 = IIT FASTA map format,\n\
                                   9 = coords in table format)\n\
");
#endif
    fprintf(stdout,"\n");

    fprintf(stdout,"\
Output options\n\
  -n, --npaths=INT               Maximum number of paths to show.  If set to 0,\n \
                                 prints two paths if chimera detected, else one.\n\
  -O, --ordered                  Print output in same order as input (relevant\n\
                                   only if there is more than one worker thread)\n\
  -5, --md5                      Print MD5 checksum for each query sequence\n\
  -o, --chimera_overlap          Overlap to show, if any, at chimera breakpoint\n\
  -V, --usesnps=STRING           Use database containing known SNPs (in <STRING>.iit, built\n\
                                   previously using snpindex) for reporting output\n\
");

#ifdef PMAP    
    fprintf(stdout,"\
  -Y, --tolerant                 Translates genome with corrections for frameshifts\n\
");
#else
    fprintf(stdout,"\
  -F, --fulllength               Assume full-length protein, starting with Met\n\
  -T, --truncate                 Truncate alignment around full-length protein, Met to Stop\n\
                                 Implies -F flag.\n\
  -Y, --tolerant                 Translates cDNA with corrections for frameshifts\n\
");
#endif

    fprintf(stdout,"\n");

    fprintf(stdout,"\
External map file options\n\
  -M, --mapdir=directory         Map directory\n\
  -m, --map=iitfile              Map file.  If argument is '?' (with the quotes),\n\
                                   this lists available map files.\n\
  -e, --mapexons                 Map each exon separately\n\
  -b, --mapboth                  Report hits from both strands of genome\n\
  -u, --flanking=INT             Show flanking hits (default 0)\n\
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
  -v, --version                  Show version\n\
  -?, --help                     Show this help message\n\
");

    return;
}

/************************************************************************/


static void
print_npaths (int npaths, Diagnostic_T diagnostic, Sequence_T usersegment, Chrsubset_T chrsubset, 
#ifdef BETATEST
	      Chimera_T chimera, 
#endif
	      Failure_T failuretype) {

  if (diagnosticp == true) {
    Diagnostic_print(diagnostic);
  }

  if (npaths == 0) {
    printf("Paths (0):");
  } else {
    printf("Paths (%d):",npaths);
  }
  Chrsubset_print(chrsubset);
  if (failuretype == NO_FAILURE) {
#ifdef BETATEST
    if (chimera != NULL) {
      Chimera_print(chimera);
    }
#endif
  } else if (failuretype == EMPTY_SEQUENCE) {
    printf(" *** Empty sequence ***");
  } else if (failuretype == SHORT_SEQUENCE) {
#ifdef PMAP
    printf(" *** Short sequence < %d aa ***",INDEX1PART_AA);
#else
    printf(" *** Short sequence < %d bp ***",INDEX1PART);
#endif
  } else if (failuretype == POOR_SEQUENCE) {
    printf(" *** Poor sequence (use -p flag to change pruning behavior) ***");
  } else if (failuretype == REPETITIVE) {
    printf(" *** Repetitive sequence (use -p flag to change pruning behavior) ***");
  }
  printf("\n");
  if (npaths == 0) {
    printf("\n");
  }
  return;
}


/* Called by output thread so queryaauc not available */
static void
print_path_summaries (Stage3_T *stage3array, int npaths, Sequence_T usersegment,
		      Sequence_T queryseq, Diagnostic_T diagnostic, Failure_T failuretype, char *dbversion,
		      Chimera_T chimera, bool zerobasedp) {
  int pathnum;
#ifdef PMAP
  bool fulllengthp = false, truncatep = false;
#endif

  if (npaths == 0) {
#ifdef BETATEST
    print_npaths(0,diagnostic,usersegment,chrsubset,/*chimera*/NULL,failuretype);
#else
    print_npaths(0,diagnostic,usersegment,chrsubset,failuretype);
#endif
  } else {
#ifdef BETATEST
    print_npaths(npaths,diagnostic,usersegment,chrsubset,chimera,NO_FAILURE);
#else
    print_npaths(npaths,diagnostic,usersegment,chrsubset,NO_FAILURE);
#endif
    if (maxpaths == 0) {
      /* Special mode */
      if (npaths > 0) {
	Stage3_print_pathsummary(stage3array[0],/*pathnum*/1,chromosome_iit,contig_iit,
				 altstrain_iit,queryseq,fulllengthp,truncatep,strictp,
				 dbversion,/*maxmutations*/1000000,zerobasedp,diagnosticp,maponlyp);
	if (chimera != NULL && npaths > 1) {
	  Stage3_print_pathsummary(stage3array[1],/*pathnum*/2,chromosome_iit,contig_iit,
				   altstrain_iit,queryseq,fulllengthp,truncatep,strictp,
				   dbversion,/*maxmutations*/1000000,zerobasedp,diagnosticp,maponlyp);
	}
      }
    } else {
      for (pathnum = 1; pathnum <= npaths && pathnum <= maxpaths; pathnum++) {
	Stage3_print_pathsummary(stage3array[pathnum-1],pathnum,chromosome_iit,contig_iit,
				 altstrain_iit,queryseq,fulllengthp,truncatep,strictp,
				 dbversion,/*maxmutations*/1000000,zerobasedp,diagnosticp,maponlyp);
      }
    }
  }

  return;
}


static int requestid = 0;

#ifdef HAVE_PTHREAD
static void *
input_thread (void *data) {
  Blackboard_T blackboard = (Blackboard_T) data;
  FILE *input = Blackboard_input(blackboard);
  char **files = Blackboard_files(blackboard);
  int nfiles = Blackboard_nfiles(blackboard);
  Sequence_T queryseq, usersegment;
  Request_T request;
  int inputid = 1;		/* Initial queryseq, inputid 0, already handled by main() */
  int nextchar = Blackboard_nextchar(blackboard);

  debug(printf("input_thread: Starting\n"));
  usersegment = Blackboard_usersegment(blackboard);
  while ((queryseq = Sequence_read_multifile(&nextchar,&input,&files,&nfiles,maponlyp)) != NULL) {
    if (inputid % batch_interval != batch_modulus) {
      Sequence_free(&queryseq);
    } else {
      debug(printf("input_thread: Putting request id %d\n",inputid));
      request = Request_new(requestid++,queryseq,usersegment);
      Blackboard_put_request(blackboard,request);
    }
    inputid++;
    debug(printf("advanced input id to %d\n",inputid));
  }

  Blackboard_set_inputdone(blackboard);
  debug(printf("input_thread: Ending\n"));
  return (void *) NULL;
}
#endif


static void
print_result (Result_T result, Request_T request) {
  Sequence_T queryseq, usersegment;
  Diagnostic_T diagnostic;
  Stage3_T *stage3array, stage3;
  int npaths, pathnum;
  Chimera_T chimera = NULL;
  int chimerapos, chimeraequivpos, chimera_cdna_direction, worker_id;
  double donor_prob, acceptor_prob;

  queryseq = Request_queryseq(request);
  usersegment = Request_usersegment(request);
  worker_id = Result_worker_id(result);

  stage3array = Result_array(&npaths,result);
  if (debug_graphic_p == true) {
    /* No output */
  } else if (compressoutp == false && cdna_exons_p == false && genomic_exons_p == false &&
#ifdef PMAP
      psloutput_pro_p == false &&
#endif
	     print_cdna_p == false && print_coordinates_p == false && 
	     protein_genomic_p == false && psloutput_nt_p == false && 
	     gffoutput_p == false && mapoutput_p == false && exon_mapoutput_p == false && 
	     splicesites_output_p == false) {
    printf(">");
    Sequence_print_header(queryseq,checksump);

  } else if (npaths == 0) {
    if (Result_failuretype(result) == POOR_SEQUENCE) {
      fprintf(stderr,"Accession %s skipped (poor sequence).  Use -p flag to change pruning behavior\n",Sequence_accession(queryseq));
    } else if (Result_failuretype(result) == REPETITIVE) {
      fprintf(stderr,"Accession %s skipped (repetitive sequence).  Use -p flag to change pruning behavior\n",Sequence_accession(queryseq));
    } else {
      fprintf(stderr,"No paths found for %s\n",Sequence_accession(queryseq));
    }
  }

  if (debug_graphic_p == true) {
    /* No output */
  } else if (compressoutp == true) {
    if ((chimera = Result_chimera(result)) == NULL) {
      chimerapos = chimeraequivpos = -1;
      chimera_cdna_direction = 0;
      donor_prob = acceptor_prob = 0.0;
    } else {
#ifdef BETATEST
      chimerapos = Chimera_pos(chimera);
      chimeraequivpos = Chimera_equivpos(chimera);
      donor_prob = Chimera_donor_prob(chimera);
      acceptor_prob = Chimera_acceptor_prob(chimera);
      chimera_cdna_direction = Chimera_cdna_direction(chimera);
#endif
    }

    if (maxpaths == 0) {
      /* Special mode */
      if (npaths > 0) {
	Stage3_print_compressed(stage3array[0],queryseq,chromosome_iit,
				dbversion,/*pathnum*/1,npaths,checksump,chimerapos,chimeraequivpos,
				donor_prob,acceptor_prob,chimera_cdna_direction,/*zerobasedp*/false,
				truncatep,strictp,worker_id);
	if (chimera != NULL && npaths > 1) {

	  Stage3_print_compressed(stage3array[1],queryseq,chromosome_iit,
				  dbversion,/*pathnum*/2,npaths,checksump,chimerapos,chimeraequivpos,
				  donor_prob,acceptor_prob,chimera_cdna_direction,/*zerobasedp*/false,
				  truncatep,strictp,worker_id);
	}
      }
    } else {
      for (pathnum = 1; pathnum <= npaths && pathnum <= maxpaths; pathnum++) {
	Stage3_print_compressed(stage3array[pathnum-1],queryseq,chromosome_iit,
				dbversion,pathnum,npaths,checksump,chimerapos,chimeraequivpos,
				donor_prob,acceptor_prob,chimera_cdna_direction,/*zerobasedp*/false,
				truncatep,strictp,worker_id);
      }
    }

  } else if (continuousp == true) {
    if (npaths == 0) {
      printf("\n\n\n");
    } else {
      stage3 = stage3array[0];
      Stage3_print_alignment(stage3,queryseq,genome,chromosome_iit,summaryonlyp,/*universalp*/false,
			     /*zerobasedp*/false,continuousp,continuous_by_exon_p,diagnosticp,strictp,/*flipgenomep*/true,
			     proteinmode,invertmode,nointronlenp,wraplength,maponlyp);
    }

  } else if (cdna_exons_p == true) {
    if (npaths > 0) {
      printf(">");
      Sequence_print_header(queryseq,checksump);
      Stage3_print_cdna_exons(stage3array[0],wraplength,ngap);
    }

  } else if (genomic_exons_p == true) {
    if (npaths > 0) {
      printf(">");
      Sequence_print_header(queryseq,checksump);
      Stage3_print_genomic_exons(stage3array[0],wraplength,ngap);
    }

  } else if (print_cdna_p == true) {
    if (maxpaths == 0) {
      /* Special mode */
      if (npaths > 0) {
	printf(">");
	Sequence_print_header(queryseq,checksump);
	Stage3_print_cdna(stage3array[0],queryseq,fulllengthp,truncatep,strictp,wraplength);
	if (Result_chimera(result) != NULL && npaths > 1) {
	  printf(">");
	  Sequence_print_header(queryseq,checksump);
	  Stage3_print_cdna(stage3array[1],queryseq,fulllengthp,truncatep,strictp,wraplength);
	}
      }
    } else {
      for (pathnum = 0; pathnum < npaths && pathnum < maxpaths; pathnum++) {
	printf(">");
	Sequence_print_header(queryseq,checksump);
	Stage3_print_cdna(stage3array[pathnum],queryseq,fulllengthp,truncatep,strictp,wraplength);
      }
    }

  } else if (protein_genomic_p == true) {
    if (maxpaths == 0) {
      /* Special mode */
      if (npaths > 0) {
	printf(">");
	Sequence_print_header(queryseq,checksump);
	Stage3_print_protein_genomic(stage3array[0],queryseq,fulllengthp,truncatep,strictp,wraplength);
	if (Result_chimera(result) != NULL && npaths > 1) {
	  printf(">");
	  Sequence_print_header(queryseq,checksump);
	  Stage3_print_protein_genomic(stage3array[1],queryseq,fulllengthp,truncatep,strictp,wraplength);
	}
      }
    } else {
      for (pathnum = 0; pathnum < npaths && pathnum < maxpaths; pathnum++) {
	printf(">");
	Sequence_print_header(queryseq,checksump);
	Stage3_print_protein_genomic(stage3array[pathnum],queryseq,fulllengthp,truncatep,strictp,wraplength);
      }
    }

  } else if (psloutput_nt_p == true) {
    if (maxpaths == 0) {
      /* Special mode */
      if (npaths > 0) {
	Stage3_print_pslformat_nt(stage3array[0],/*pathnum*/1,chromosome_iit,usersegment,queryseq);
	if (Result_chimera(result) != NULL && npaths > 1) {
	  Stage3_print_pslformat_nt(stage3array[1],/*pathnum*/2,chromosome_iit,usersegment,queryseq);
	}
      }
    } else {
      for (pathnum = 0; pathnum < npaths && pathnum < maxpaths; pathnum++) {
	Stage3_print_pslformat_nt(stage3array[pathnum],/*pathnum*/pathnum+1,chromosome_iit,usersegment,queryseq);
      }
    }

#ifdef PMAP
  } else if (psloutput_pro_p == true) {
    if (maxpaths == 0) {
      /* Special mode */
      if (npaths > 0) {
	Stage3_print_pslformat_pro(stage3array[0],/*pathnum*/1,chromosome_iit,usersegment,queryseq,strictp);
	if (Result_chimera(result) != NULL && npaths > 1) {
	  Stage3_print_pslformat_pro(stage3array[1],/*pathnum*/2,chromosome_iit,usersegment,queryseq,strictp);
	}
      }
    } else {
      for (pathnum = 0; pathnum < npaths && pathnum < maxpaths; pathnum++) {
	Stage3_print_pslformat_pro(stage3array[pathnum],/*pathnum*/pathnum+1,chromosome_iit,usersegment,queryseq,strictp);
      }
    }
#endif

  } else if (gffoutput_p == true) {
    if (maxpaths == 0) {
      /* Special mode */
      if (npaths > 0) {
	Stage3_print_gff3(stage3array[0],/*pathnum*/1,chromosome_iit,queryseq,dbversion,
			  diagnosticp,fulllengthp,truncatep,strictp,gff_gene_format_p,user_genomicseg);
	if (Result_chimera(result) != NULL && npaths > 1) {
	  Stage3_print_gff3(stage3array[1],/*pathnum*/2,chromosome_iit,queryseq,dbversion,
			    diagnosticp,fulllengthp,truncatep,strictp,gff_gene_format_p,user_genomicseg);
	}
      }
    } else {
      for (pathnum = 0; pathnum < npaths && pathnum < maxpaths; pathnum++) {
	Stage3_print_gff3(stage3array[pathnum],/*pathnum*/pathnum+1,chromosome_iit,queryseq,
			  dbversion,diagnosticp,fulllengthp,truncatep,strictp,gff_gene_format_p,user_genomicseg);
      }
    }

  } else if (mapoutput_p == true) {
    if (maxpaths == 0) {
      /* Special mode */
      if (npaths > 0) {
	Stage3_print_iit_map(stage3array[0],/*pathnum*/1,chromosome_iit,queryseq);
	if (Result_chimera(result) != NULL && npaths > 1) {
	  Stage3_print_iit_map(stage3array[1],/*pathnum*/2,chromosome_iit,queryseq);
	}
      }
    } else {
      for (pathnum = 0; pathnum < npaths && pathnum < maxpaths; pathnum++) {
	Stage3_print_iit_map(stage3array[pathnum],/*pathnum*/pathnum+1,chromosome_iit,queryseq);
      }
    }

  } else if (exon_mapoutput_p == true) {
    if (maxpaths == 0) {
      /* Special mode */
      if (npaths > 0) {
	Stage3_print_iit_exon_map(stage3array[0],/*pathnum*/1,chromosome_iit,queryseq);
	if (Result_chimera(result) != NULL && npaths > 1) {
	  Stage3_print_iit_exon_map(stage3array[1],/*pathnum*/2,chromosome_iit,queryseq);
	}
      }
    } else {
      for (pathnum = 0; pathnum < npaths && pathnum < maxpaths; pathnum++) {
	Stage3_print_iit_exon_map(stage3array[pathnum],/*pathnum*/pathnum+1,chromosome_iit,queryseq);
      }
    }

  } else if (splicesites_output_p == true) {
    /* Print only best path */
    if (npaths > 0) {
      Stage3_print_splicesites(stage3array[0],chromosome_iit,queryseq);
    }

  } else if (print_coordinates_p == true) {
    if (npaths > 0) {
      printf(">");
      Sequence_print_header(queryseq,checksump);
      Stage3_print_coordinates(stage3array[0],queryseq,chromosome_iit,
			       /*zerobasedp*/false,invertmode,fulllengthp,truncatep,strictp,maponlyp);
    }

  } else {
    /* Usual output */
    diagnostic = Result_diagnostic(result);
    if (npaths == 0) {
#ifdef BETATEST
      print_npaths(0,diagnostic,usersegment,chrsubset,/*chimera*/NULL,Result_failuretype(result));
#else
      print_npaths(0,diagnostic,usersegment,chrsubset,Result_failuretype(result));
#endif
    } else {
      chimera = Result_chimera(result);
      print_path_summaries(stage3array,npaths,usersegment,queryseq,diagnostic,Result_failuretype(result),
			   dbversion,chimera,/*zerobasedp*/false);
    
      if (summaryonlyp == true || showalignp == true) {
	printf("Alignments:\n");
	if (maxpaths == 0) {
	  /* Special mode */
	  if (npaths > 0) {
	    printf("  Alignment for path 1:\n\n");
	    Stage3_print_alignment(stage3array[0],queryseq,genome,chromosome_iit,summaryonlyp,/*universalp*/false,
				   /*zerobasedp*/false,continuousp,continuous_by_exon_p,diagnosticp,strictp,
				   /*flipgenomep*/true,proteinmode,invertmode,nointronlenp,wraplength,maponlyp);
	    if (chimera != NULL && npaths > 1) {
	      printf("  Alignment for path 2:\n\n");
	      Stage3_print_alignment(stage3array[1],queryseq,genome,chromosome_iit,summaryonlyp,/*universalp*/false,
				     /*zerobasedp*/false,continuousp,continuous_by_exon_p,diagnosticp,strictp,
				     /*flipgenomep*/true,proteinmode,invertmode,nointronlenp,wraplength,maponlyp);
	    }
	  }
	} else {
	  for (pathnum = 1; pathnum <= npaths && pathnum <= maxpaths; pathnum++) {
	    printf("  Alignment for path %d:\n\n",pathnum);
	    Stage3_print_alignment(stage3array[pathnum-1],queryseq,genome,chromosome_iit,summaryonlyp,/*universalp*/false,
				   /*zerobasedp*/false,continuousp,continuous_by_exon_p,diagnosticp,strictp,
				   /*flipgenomep*/true,proteinmode,invertmode,nointronlenp,wraplength,maponlyp);
	  }
	}
      }

      if (map_iit != NULL) {
	printf("Maps:\n");
	if (maxpaths == 0) {
	  /* Special mode */
	  if (npaths > 0) {
	    Stage3_print_map(stage3array[0],map_iit,map_divint_crosstable,map_iit_universal_p,
			     map_iit_forward_type,map_iit_reverse_type,
			     chromosome_iit,/*pathnum*/1,map_exons_p,map_bothstrands_p,
			     nflanking);
	    if (chimera != NULL && npaths > 1) {
	      Stage3_print_map(stage3array[1],map_iit,map_divint_crosstable,map_iit_universal_p,
			       map_iit_forward_type,map_iit_reverse_type,
			       chromosome_iit,/*pathnum*/2,map_exons_p,map_bothstrands_p,
			       nflanking);
	    }
	  }
	} else {
	  for (pathnum = 1; pathnum <= npaths && pathnum <= maxpaths; pathnum++) {
	    Stage3_print_map(stage3array[pathnum-1],map_iit,map_divint_crosstable,map_iit_universal_p,
			     map_iit_forward_type,map_iit_reverse_type,
			     chromosome_iit,pathnum,map_exons_p,map_bothstrands_p,
			     nflanking);
	  }
	}
      }

    }
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
  List_T request_queue = NULL, result_queue = NULL;
  void *item;
  int outputid = 0;

  debug(printf("output_thread: Starting\n"));
  while ((result = Blackboard_get_result(&request,blackboard)) != NULL) {
    if (Result_id(result) != outputid) {
      result_queue = queue_insert_result(result_queue,result);
      request_queue = queue_insert_request(request_queue,request);
    } else {
      print_result(result,request);
      Result_free(&result);
      Request_free(&request);
      outputid++;

      while (result_queue != NULL && Result_id(List_head(result_queue)) == outputid) {
	result_queue = List_pop(result_queue,&item);
	result = (Result_T) item;
	request_queue = List_pop(request_queue,&item);
	request = (Request_T) item;

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

    print_result(result,request);
    Result_free(&result);
    Request_free(&request);
  }

  debug(printf("output_thread: Ending\n"));
  fprintf(stderr,"Processed %d queries\n",outputid);
  return (void *) NULL;
}

#endif
    
#if 0
static void
output_one (Result_T result, Request_T request) {
  print_result(result,request);
  Result_free(&result);
  return;
}
#endif


static Stage3_T *
stage3array_from_list (int *npaths, List_T stage3list) {
  Stage3_T *array1, *array0, x, y;
  bool *eliminate;
  int norig, i, j;

  Stage3_recompute_goodness(stage3list);

  norig = List_length(stage3list);
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
      Stage3_free(&x);
    } else {
      array1[j++] = x;
    }
  }
  FREE(array0);
  FREE(eliminate);

  return array1;
}


static List_T
update_stage3list (List_T stage3list, bool lowidentityp, Sequence_T queryseq, Sequence_T queryuc, Sequence_T genomicseg, 
		   Genomicpos_T genomicstart, Genomicpos_T genomiclength,
		   Oligoindex_T *oligoindices_major, int noligoindices_major,
		   Oligoindex_T *oligoindices_minor, int noligoindices_minor,
		   Pairpool_T pairpool, Intpool_T intpool, Diagpool_T diagpool,
		   int straintype, char *strain, Genome_T genome,
		   Chrnum_T chrnum, Genomicpos_T chroffset, Genomicpos_T chrpos, Genomicpos_T chrlength,
		   bool watsonp, Gbuffer_T gbuffer, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
		   Stopwatch_T stopwatch) {
  bool do_final_p;

#ifdef PMAP
  Sequence_T queryntseq;
#endif
  Sequence_T genomicuc;
  List_T path;
  Stage3_T stage3;

#ifdef PMAP
  queryntseq = Sequence_convert_to_nucleotides(queryseq);
#endif

  if (user_genomicseg == NULL && uncompressedp == false && straintype == 0) {
    genomicuc = Sequence_alias(genomicseg);
  } else {
    genomicuc = Sequence_uppercase(genomicseg);
  }

  debug(printf("Beginning Stage2_compute\n"));
  path = Stage2_compute(Sequence_trimpointer(queryseq),Sequence_trimpointer(queryuc),
			Sequence_trimlength(queryseq),/*query_offset*/0,
			Sequence_fullpointer(genomicseg),Sequence_fullpointer(genomicuc),
			Sequence_fulllength(genomicseg),/*genomic_offset*/0,
			oligoindices_major,noligoindices_major,pairpool,diagpool,
			sufflookback,nsufflookback,maxintronlen_bound,
			/*localp*/true,/*skip_repetitive_p*/true,debug_graphic_p,diagnosticp,stopwatch);

  if (path != NULL) {
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

    if ((stage3 = Stage3_compute(path,genomicstart,genomiclength,
#ifdef PMAP
				 queryseq,queryntseq,queryntseq,
#else
				 queryseq,queryuc,
#endif
				 genomicseg,genomicuc,straintype,strain,genome,genomealt,chrnum,chroffset,
				 chrpos,chrlength,watsonp,maxpeelback,nullgap,extramaterial_end,extramaterial_paired,
				 extraband_single,extraband_end,extraband_paired,end_microexons_p,
				 minendexon,pairpool,gbuffer,dynprogL,dynprogM,dynprogR,altstrain_iit,ngap,
				 stage3debug,diagnosticp,checkp,stopwatch,do_final_p,sense,
				 oligoindices_minor,noligoindices_minor,diagpool,
				 sufflookback,nsufflookback,maxintronlen_bound)) != NULL) {


      stage3list = List_push(stage3list,stage3);
    }
  }

  Sequence_free(&genomicuc);
#ifdef PMAP
  Sequence_free(&queryntseq);
#endif

  return stage3list;
}

static List_T
update_stage3list_maponlyp (List_T stage3list, Gregion_T gregion, Sequence_T queryseq, Sequence_T queryuc,
			    Pairpool_T pairpool, int straintype, char *strain, Genome_T genome,
			    Chrnum_T chrnum,  Genomicpos_T chroffset, Genomicpos_T chrpos, Genomicpos_T chrlength,
			    bool watsonp, Gbuffer_T gbuffer, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR) {
  Stage3_T stage3;
#ifdef PMAP
  Sequence_T queryntseq;
#endif

#ifdef PMAP
  queryntseq = Sequence_convert_to_nucleotides(queryseq);
#endif

  if ((stage3 = Stage3_direct(gregion,
#ifdef PMAP
			      queryseq,queryntseq,queryntseq,
#else
			      queryseq,queryuc,
#endif
			      pairpool,genome,chrnum,chroffset,chrpos,watsonp,ngap,
			      gbuffer,dynprogL,dynprogR,extramaterial_end,extraband_end)) != NULL) {
    stage3list = List_push(stage3list,stage3);
  }

  return stage3list;
}

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


static Stage3_T *
stage3_from_usersegment (int *npaths, bool lowidentityp, Sequence_T queryseq, Sequence_T queryuc, Sequence_T usersegment,
			 Oligoindex_T *oligoindices_major, int noligoindices_major,
			 Oligoindex_T *oligoindices_minor, int noligoindices_minor,
			 Pairpool_T pairpool, Intpool_T intpool, Diagpool_T diagpool,
			 Gbuffer_T gbuffer, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
			 Stopwatch_T stopwatch) {
  List_T stage3list;
  Genomicpos_T chroffset, chrpos, chrlength;
  Sequence_T revcomp;
  Chrnum_T chrnum = 0;
		    
  chroffset = chrpos = 0U;
  chrlength = Sequence_fulllength(usersegment);

  stage3list = update_stage3list(/*stage3list*/NULL,lowidentityp,queryseq,queryuc,usersegment,
				 /*genomicstart*/0U,/*genomiclength*/Sequence_fulllength(usersegment),
				 oligoindices_major,noligoindices_major,oligoindices_minor,noligoindices_minor,
				 pairpool,intpool,diagpool,
				 /*straintype*/0,/*strain*/NULL,/*genome*/NULL,chrnum,chroffset,chrpos,chrlength,
				 /*watsonp*/true,gbuffer,dynprogL,dynprogM,dynprogR,stopwatch);

  revcomp = Sequence_revcomp(usersegment);

  stage3list = update_stage3list(stage3list,lowidentityp,queryseq,queryuc,revcomp,
				 /*genomicstart*/0U,/*genomiclength*/Sequence_fulllength(usersegment),
				 oligoindices_major,noligoindices_major,oligoindices_minor,noligoindices_minor,
				 pairpool,intpool,diagpool,
				 /*straintype*/0,/*strain*/NULL,/*genome*/NULL,chrnum,chroffset,chrpos,chrlength,
				 /*watsonp*/false,gbuffer,dynprogL,dynprogM,dynprogR,stopwatch);

  Sequence_free(&revcomp);

  if (stage3list == NULL) {
    *npaths = 0;
    return NULL;
  } else {
    return stage3array_from_list(&(*npaths),stage3list);
  }
}


static Stage3_T *
stage3_from_gregions (int *npaths, Stage3_T *oldstage3array, List_T gregions, bool lowidentityp, Sequence_T queryseq,
		      Sequence_T queryuc, Sequence_T usersegment, 
		      Oligoindex_T *oligoindices_major, int noligoindices_major,
		      Oligoindex_T *oligoindices_minor, int noligoindices_minor,
		      Pairpool_T pairpool, Intpool_T intpool, Diagpool_T diagpool,
		      Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
		      Gbuffer_T gbuffer, Stopwatch_T stopwatch) {
  List_T stage3list = NULL;
  Gregion_T gregion;
  char *strain;
  Sequence_T genomicseg;
  int *indexarray, nindices, straintype, i, j;
  void *item;
		    
  for (i = 0; i < *npaths; i++) {
    stage3list = List_push(stage3list,(void *) oldstage3array[i]);
  }
  if (oldstage3array != NULL) {
    FREE(oldstage3array);
  }

  while (gregions != NULL) {
    gregions = List_pop(gregions,&item);
    gregion = (Gregion_T) item;

    /* if (Match_usep(match) == true) { */
    if (1) {
      if (usersegment != NULL) {
	/* chrlength = Sequence_fulllength(usersegment); */
	strain = NULL;
	Gbuffer_alloc_contents(gbuffer,Gregion_genomiclength(gregion));
	genomicseg = Sequence_substring(usersegment,Gregion_genomicstart(gregion),Gregion_genomiclength(gregion),
					Gregion_revcompp(gregion),
					Gbuffer_chars1(gbuffer),Gbuffer_chars2(gbuffer),Gbuffer_gbufferlen(gbuffer));
	stage3list = update_stage3list(stage3list,lowidentityp,queryseq,queryuc,genomicseg,Gregion_genomicstart(gregion),
				       Gregion_genomiclength(gregion),oligoindices_major,noligoindices_major,
				       oligoindices_minor,noligoindices_minor,pairpool,intpool,diagpool,
				       /*straintype*/0,/*strain*/NULL,/*genome*/NULL,
				       Gregion_chrnum(gregion),Gregion_chroffset(gregion),
				       Gregion_chrpos(gregion),Gregion_chrlength(gregion),Gregion_plusp(gregion),
				       gbuffer,dynprogL,dynprogM,dynprogR,stopwatch);
	Sequence_free(&genomicseg);
	Gbuffer_free_contents(gbuffer);

      } else if (maponlyp == true) {
	Gbuffer_alloc_contents(gbuffer,Gregion_genomiclength(gregion));
	stage3list = update_stage3list_maponlyp(stage3list,gregion,queryseq,queryuc,pairpool,
						/*straintype*/0,/*strain*/NULL,genome,
						Gregion_chrnum(gregion),Gregion_chroffset(gregion),
						Gregion_chrpos(gregion),Gregion_chrlength(gregion),Gregion_plusp(gregion),
						gbuffer,dynprogL,dynprogM,dynprogR);
	Gbuffer_free_contents(gbuffer);

      } else {
	Gbuffer_alloc_contents(gbuffer,Gregion_genomiclength(gregion));
	if (diagnosticp == true) {
	  printf("Got sequence at %u with length %u, revcomp %d\n",
		 Gregion_genomicstart(gregion),Gregion_genomiclength(gregion),Gregion_revcompp(gregion));
	}
	if (genomealt != NULL) {
	  genomicseg = Genome_get_segment_alt(genomealt,Gregion_genomicstart(gregion),Gregion_genomiclength(gregion),
					      /*chromosome_iit*/NULL,Gregion_revcompp(gregion),
					      Gbuffer_chars1(gbuffer),Gbuffer_chars2(gbuffer),Gbuffer_gbufferlen(gbuffer));
	} else {
	  genomicseg = Genome_get_segment(genome,Gregion_genomicstart(gregion),Gregion_genomiclength(gregion),
					  /*chromosome_iit*/NULL,Gregion_revcompp(gregion),
					  Gbuffer_chars1(gbuffer),Gbuffer_chars2(gbuffer),Gbuffer_gbufferlen(gbuffer));
	}
	stage3list = update_stage3list(stage3list,lowidentityp,queryseq,queryuc,genomicseg,Gregion_genomicstart(gregion),
				       Gregion_genomiclength(gregion),oligoindices_major,noligoindices_major,
				       oligoindices_minor,noligoindices_minor,pairpool,intpool,diagpool,
				       /*straintype*/0,/*strain*/NULL,genome,
				       Gregion_chrnum(gregion),Gregion_chroffset(gregion),
				       Gregion_chrpos(gregion),Gregion_chrlength(gregion),Gregion_plusp(gregion),
				       gbuffer,dynprogL,dynprogM,dynprogR,stopwatch);
	Sequence_free(&genomicseg);

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
	      stage3list = update_stage3list(stage3list,lowidentityp,queryseq,queryuc,genomicseg,Gregion_genomicstart(gregion),
					     Gregion_genomiclength(gregion),oligoindices_major,noligoindices_major,
					     oligoindices_minor,noligoindices_minor,pairpool,intpool,diagpool,
					     straintype,strain,genome,
					     Gregion_chrnum(gregion),Gregion_chroffset(gregion),
					     Gregion_chrpos(gregion),Gregion_chrlength(gregion),Gregion_plusp(gregion),
					     gbuffer,dynprogL,dynprogM,dynprogR,stopwatch);
	      Sequence_free(&genomicseg);
	    }
	    FREE(indexarray);
	  }
	}

	Gbuffer_free_contents(gbuffer);
      }
    }
    Gregion_free(&gregion);
  }
	
  if (stage3list == NULL) {
    *npaths = 0;
    return NULL;
  } else {
    return stage3array_from_list(&(*npaths),stage3list);
  }
}


#ifdef BETATEST
static bool
chimeric_join5_p (Stage3_T from, int effective_start, int effective_end) {
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

static void
chimera_separate_paths (Stage3_T **stage3array_sub1, int *npaths_sub1, 
			Stage3_T **stage3array_sub2, int *npaths_sub2,
			Stage3_T *stage3array, int npaths, int five_margin, int three_margin,
			int effective_start, int effective_end) {
  int i, njoinable_left = 0, njoinable_right = 0, j = 0, k = 0;

  for (i = 1; i < npaths; i++) {
    if (chimeric_join3_p(stage3array[i],effective_start,effective_end) == true) {
      debug2(printf("Found join for 0 to %d\n",i));
      njoinable_right++;
    }
  }
  for (i = 1; i < npaths; i++) {
    if (chimeric_join5_p(stage3array[i],effective_start,effective_end) == true) {
      debug2(printf("Found join for %d to 0\n",i));
      njoinable_left++;
    }
  }

  if (njoinable_left == 0 && njoinable_right == 0) {
    *stage3array_sub1 = (Stage3_T *) NULL;
    *npaths_sub1 = 0;
    *stage3array_sub2 = (Stage3_T *) NULL;
    *npaths_sub2 = 0;

  } else if (njoinable_left == njoinable_right) {
    if (five_margin >= three_margin) {
      /* Use njoinable_left */
      njoinable_right = 0;
    } else {
      /* Use njoinable_right */
      njoinable_left = 0;

    }
  }
  debug2(printf("Separation process yields %d paths for sub1 and %d paths for sub2\n",njoinable_right,njoinable_left));

  if (njoinable_right > njoinable_left) {
    *npaths_sub2 = njoinable_right;
    *stage3array_sub2 = (Stage3_T *) CALLOC(njoinable_right,sizeof(Stage3_T));
    *npaths_sub1 = npaths - njoinable_right;
    *stage3array_sub1 = (Stage3_T *) CALLOC(npaths - njoinable_right,sizeof(Stage3_T));
    (*stage3array_sub1)[j++] = stage3array[0];
    for (i = 1; i < npaths; i++) {
      if (chimeric_join3_p(stage3array[i],effective_start,effective_end) == true) {
	(*stage3array_sub2)[k++] = stage3array[i];
      } else {
	(*stage3array_sub1)[j++] = stage3array[i];
      }
    }

  } else if (njoinable_left > njoinable_right) {
    *npaths_sub1 = njoinable_left;
    *stage3array_sub1 = (Stage3_T *) CALLOC(njoinable_left,sizeof(Stage3_T));
    *npaths_sub2 = npaths - njoinable_left;
    *stage3array_sub2 = (Stage3_T *) CALLOC(npaths - njoinable_left,sizeof(Stage3_T));
    (*stage3array_sub2)[k++] = stage3array[0];
    for (i = 1; i < npaths; i++) {
      if (chimeric_join5_p(stage3array[i],effective_start,effective_end) == true) {
	(*stage3array_sub1)[j++] = stage3array[i];
      } else {
	(*stage3array_sub2)[k++] = stage3array[i];
      }
    }
  }

  return;
}

#if 0
      /* Check to see if we can merge chimeric parts */
      if (Stage3_mergeable(stage3array[0],stage3array[1],Chimera_cdna_direction(*chimera),
			   Chimera_donor_prob(*chimera),Chimera_acceptor_prob(*chimera),
			   pairpool,genome,ngap) == true) {
	best1 = stage3array[1];
	Stage3_free(&best1);
	for (j = 2; j < *npaths; j++) {
	  stage3array[j-1] = stage3array[j];
	}
	stage3array[(*npaths)-1] = (Stage3_T) NULL;
	*npaths -= 1;
	Chimera_free(&(*chimera));
      }
#endif


static Stage3_T *
merge_left_and_right (int *npaths, Stage3_T *stage3array_sub1, int npaths_sub1, int bestfrom,
		      Stage3_T *stage3array_sub2, int npaths_sub2, int bestto, int chimerapos,
		      int chimeraequivpos, int queryntlength, bool boundp) {
  Stage3_T *newstage3array, best0, best1;
  int i, k = 2;

  if (boundp == false) {
    best0 = stage3array_sub1[bestfrom];
    best1 = stage3array_sub2[bestto];
  } else {
    best0 = Stage3_apply_bounds(stage3array_sub1[bestfrom],0,chimeraequivpos+chimera_overlap,/*revertp*/false);
    best1 = Stage3_apply_bounds(stage3array_sub2[bestto],chimerapos+1-chimera_overlap,queryntlength,/*revertp*/false);
  }

  if (best0 == NULL || best1 == NULL) {
    debug2(printf("Applying bounds kills one of the paths\n"));
    abort();
    return (Stage3_T *) NULL;
  } else {
    debug2(printf("Rearranging paths\n"));

    *npaths = npaths_sub1 + npaths_sub2;
    newstage3array = (Stage3_T *) CALLOC(*npaths,sizeof(Stage3_T));
    newstage3array[0] = best0;
    newstage3array[1] = best1;

    for (i = 0; i < npaths_sub1; i++) {
      if (i != bestfrom) {
	newstage3array[k++] = stage3array_sub1[i];
      }
    }

    for (i = 0; i < npaths_sub2; i++) {
      if (i != bestto) {
	newstage3array[k++] = stage3array_sub2[i];
      }
    }
    return newstage3array;
  }
}


static Stage3_T *
check_for_chimera (Chimera_T *chimera, int *npaths, int effective_start, int effective_end, Stage3_T nonchimericbest, 
		   Stage3_T *stage3array, Sequence_T queryseq, Sequence_T queryuc, Sequence_T usersegment, 
		   int maxpeelback, int sufflookback, int nsufflookback, int nullgap,
		   Oligoindex_T *oligoindices_major, int noligoindices_major,
		   Oligoindex_T *oligoindices_minor, int noligoindices_minor,
		   Matchpool_T matchpool, Pairpool_T pairpool, Intpool_T intpool,
		   Diagpool_T diagpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR, Gbuffer_T gbuffer) {
  List_T gregions = NULL;
  Stage3_T *stage3array_sub1 = NULL, *stage3array_sub2 = NULL;
  Sequence_T querysubseq = NULL, querysubuc = NULL;
  Diagnostic_T diagnostic;
  int chimerapos, chimeraequivpos;
  int bestfrom, bestto;
  int five_margin, three_margin, five_score = 0, three_score = 0;
  int npaths_sub1 = 0, npaths_sub2 = 0;
  int  queryntlength;
  bool lowidentityp;

  int exonexonpos, cdna_direction;
  double donor_prob, acceptor_prob;

  five_margin = effective_start - Sequence_trim_start(queryseq);
  three_margin = Sequence_trim_end(queryseq) - effective_end;
  debug2(printf("Margins are %d on the 5' end and %d on the 3' end\n",
		five_margin,three_margin));

  queryntlength = Sequence_ntlength(queryseq);
  chimera_separate_paths(&stage3array_sub1,&npaths_sub1,&stage3array_sub2,&npaths_sub2,
			 stage3array,*npaths,five_margin,three_margin,effective_start,effective_end);

  if (npaths_sub1 == 0 && npaths_sub2 == 0) {
    /* Need to compute on margin explicitly */
    if (five_margin < chimera_margin && three_margin < chimera_margin) {
      debug2(printf("Insufficient margins\n"));
    } else if (five_margin > three_margin) {
      if ((querysubseq = Sequence_subsequence(queryseq,0,effective_start+CHIMERA_SLOP)) != NULL) {
	if ((querysubuc = Sequence_subsequence(queryuc,0,effective_start+CHIMERA_SLOP)) != NULL) {
	  debug2(printf("5 margin > 3 margin.  "));
	  debug2(printf("Beginning Stage1_compute on 5' margin (%d..%d)\n",0,effective_start+CHIMERA_SLOP));
	  debug2(Sequence_print(querysubseq,/*uppercasep*/true,wraplength,/*trimmedp*/true));

	  diagnostic = Diagnostic_new();
	  gregions = Stage1_compute(&lowidentityp,querysubuc,
#ifdef PMAP
				    indexdb_fwd,indexdb_rev,
#else
				    indexdb,
#endif
				    indexdb_size_threshold,chromosome_iit,chrsubset,matchpool,
				    maxintronlen_bound,maxtotallen_bound,stuttercycles,
				    stutterhits,/*subsequencep*/true,
				    genome,genome_totallength,diagnostic,/*stopwatch*/NULL);
	  Diagnostic_free(&diagnostic);
	  debug2(printf("Performing Stage 3\n"));
	  stage3array = stage3_from_gregions(&(*npaths),stage3array,gregions,lowidentityp,querysubseq,
					     querysubuc,usersegment,oligoindices_major,noligoindices_major,
					     oligoindices_minor,noligoindices_minor,pairpool,intpool,diagpool,
					     dynprogL,dynprogM,dynprogR,gbuffer,/*stopwatch*/NULL);
	  /* Above function frees gregions */
	  Sequence_free(&querysubuc);

	  chimera_separate_paths(&stage3array_sub1,&npaths_sub1,&stage3array_sub2,&npaths_sub2,
				 stage3array,*npaths,five_margin,three_margin,effective_start,effective_end);
	}
	Sequence_free(&querysubseq);
      }
	
    } else {
      if ((querysubseq = Sequence_subsequence(queryseq,effective_end-CHIMERA_SLOP,queryntlength)) != NULL) {
	if ((querysubuc = Sequence_subsequence(queryuc,effective_end-CHIMERA_SLOP,queryntlength)) != NULL) {
	  debug2(printf("5 margin <= 3 margin.  "));
	  debug2(printf("Beginning Stage1_compute on 3' margin (%d..%d)\n",effective_end-CHIMERA_SLOP,queryntlength));
	  debug2(Sequence_print(querysubseq,/*uppercasep*/true,wraplength,/*trimmedp*/true));

	  diagnostic = Diagnostic_new();
	  gregions = Stage1_compute(&lowidentityp,querysubuc,
#ifdef PMAP
				    indexdb_fwd,indexdb_rev,
#else
				    indexdb,
#endif
				    indexdb_size_threshold,chromosome_iit,chrsubset,matchpool,
				    maxintronlen_bound,maxtotallen_bound,stuttercycles,
				    stutterhits,/*subsequencep*/true,
				    genome,genome_totallength,diagnostic,/*stopwatch*/NULL);
	  Diagnostic_free(&diagnostic);
	  debug2(printf("Performing Stage 3\n"));
	  stage3array = stage3_from_gregions(&(*npaths),stage3array,gregions,lowidentityp,querysubseq,
					     querysubuc,usersegment,oligoindices_major,noligoindices_major,
					     oligoindices_minor,noligoindices_minor,pairpool,intpool,diagpool,
					     dynprogL,dynprogM,dynprogR,gbuffer,/*stopwatch*/NULL);
	  /* Above function frees gregions */
	  Sequence_free(&querysubuc);

	  chimera_separate_paths(&stage3array_sub1,&npaths_sub1,&stage3array_sub2,&npaths_sub2,
				 stage3array,*npaths,five_margin,three_margin,effective_start,effective_end);
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
						 chimerapos,queryntlength);
    */
    debug2(printf("Chimera_bestpath returns boundary at %d..%d\n",chimerapos,chimeraequivpos));

    debug3(printf("Merging left and right early\n"));
    debug3(return merge_left_and_right(&(*npaths),stage3array_sub1,npaths_sub1,bestfrom,
				       stage3array_sub2,npaths_sub2,bestto,chimerapos,chimeraequivpos,
				       queryntlength,/*boundp*/false));
      

    if (Chimera_exonexon_p(&exonexonpos,&cdna_direction,&donor_prob,&acceptor_prob,
			   stage3array_sub1[bestfrom],stage3array_sub2[bestto],
			   genome,chromosome_iit,Sequence_fulllength(queryseq)) == true) {
      debug2(printf("Exon-exon boundary found at %d\n",exonexonpos));
      chimerapos = chimeraequivpos = exonexonpos;
    } else {
      debug2(printf("No exon-exon boundary found\n"));
      exonexonpos = -1;
      cdna_direction = 0;
      donor_prob = 0.0;
      acceptor_prob = 0.0;
    }

    if (Stage3_test_bounds(stage3array_sub1[bestfrom],0,chimeraequivpos+chimera_overlap) == true &&
	Stage3_test_bounds(stage3array_sub2[bestto],chimerapos+1-chimera_overlap,queryntlength) == true) {
      *chimera = Chimera_new(chimerapos,chimeraequivpos,exonexonpos,cdna_direction,
			     donor_prob,acceptor_prob);
      FREE(stage3array);
      debug2(printf("Merging left and right\n"));
      stage3array = merge_left_and_right(&(*npaths),stage3array_sub1,npaths_sub1,bestfrom,
					 stage3array_sub2,npaths_sub2,bestto,chimerapos,chimeraequivpos,
					 queryntlength,/*boundp*/true);
    }
    FREE(stage3array_sub2);
    FREE(stage3array_sub1);
  }

  return stage3array;
}


#endif

static Stage3_T *
apply_stage3 (Chimera_T *chimera, int *npaths, List_T gregions, bool lowidentityp, Sequence_T queryseq, Sequence_T queryuc,
	      Sequence_T usersegment, Oligoindex_T *oligoindices_major, int noligoindices_major,
	      Oligoindex_T *oligoindices_minor, int noligoindices_minor,
	      Matchpool_T matchpool, Pairpool_T pairpool, Intpool_T intpool,
	      Diagpool_T diagpool, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR, Gbuffer_T gbuffer,
	      Stopwatch_T stopwatch) {
  Stage3_T *stage3array = NULL, nonchimericbest;
  bool testchimerap = false;
#ifdef BETATEST
  int effective_start, effective_end;
#endif
  
  *npaths = 0;
  *chimera = NULL;


  stage3array = stage3_from_gregions(&(*npaths),stage3array,gregions,lowidentityp,queryseq,queryuc,
				     usersegment,oligoindices_major,noligoindices_major,
				     oligoindices_minor,noligoindices_minor,pairpool,intpool,diagpool,
				     dynprogL,dynprogM,dynprogR,gbuffer,stopwatch);

  if (stage3array != NULL) {
    nonchimericbest = stage3array[0]; 

    if (chimera_margin < 0) {
      debug2(printf("turned off\n"));
      testchimerap = false;

#ifndef BETATEST
    } else {
      fprintf(stderr,"Sorry, chimera mode is undergoing further development\n");
      exit(9);
#else
    } else if (Stage3_domain(nonchimericbest) < chimera_margin) {
      debug2(printf("Existing alignment is too short, so won't look for chimera\n"));
      testchimerap = false;

    } else if (Stage3_fracidentity(nonchimericbest) < CHIMERA_IDENTITY &&
	       Chimera_alignment_break(&effective_start,&effective_end,nonchimericbest,Sequence_ntlength(queryseq),CHIMERA_FVALUE) >= chimera_margin
	       ) {
      debug2(printf("Break in alignment quality at %d..%d detected, so will look for chimera\n",
		    effective_start,effective_end));
      testchimerap = true;

    } else if (Stage3_largemargin(&effective_start,&effective_end,nonchimericbest,Sequence_ntlength(queryseq)) >= chimera_margin) {
      debug2(printf("Large margin at %d..%d detected, so will look for chimera\n",
		    effective_start,effective_end));
      testchimerap = true;

    } else {
      debug2(printf("Good alignment already with identity %f, so won't look for chimera\n",
		    Stage3_fracidentity(nonchimericbest)));
      testchimerap = false;
#endif

    }

    if (testchimerap == true) {
#ifdef BETATEST
      debug2(printf("Checking for chimera\n"));
      stage3array = check_for_chimera(&(*chimera),&(*npaths),effective_start,effective_end,nonchimericbest,stage3array,
				      queryseq,queryuc,usersegment,
				      maxpeelback,sufflookback,nsufflookback,nullgap,
				      oligoindices_major,noligoindices_major,oligoindices_minor,noligoindices_minor,
				      matchpool,pairpool,intpool,diagpool,
				      dynprogL,dynprogM,dynprogR,gbuffer);
#endif
    }
  }

  return stage3array;
}


static void
handle_request (Request_T request, int worker_id, Matchpool_T matchpool, Pairpool_T pairpool, 
		Intpool_T intpool, Diagpool_T diagpool, 
		Oligoindex_T *oligoindices_major, int noligoindices_major,
		Oligoindex_T *oligoindices_minor, int noligoindices_minor,
		Gbuffer_T gbuffer, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR, Reqpost_T reqpost,
		Stopwatch_T stopwatch) {
  int jobid;
  Result_T result;
  Diagnostic_T diagnostic;
  Sequence_T queryseq, queryuc, usersegment;
  Chimera_T chimera = NULL;
  bool lowidentityp;
#ifndef PMAP
  bool repetitivep = false, poorp = false;
#endif

  List_T gregions;
  Stage3_T *stage3array;
  int npaths;

  jobid = Request_id(request);
  usersegment = Request_usersegment(request);
  debug(printf("worker_thread %d: Got request id\n",worker_id,jobid));
  queryseq = Request_queryseq(request);
  Matchpool_reset(matchpool);
  Pairpool_reset(pairpool);

  if (Sequence_fulllength_given(queryseq) <= 0) {
    result = Result_new(jobid,worker_id,(Chimera_T) NULL,(Stage3_T *) NULL,0,/*diagnostic*/NULL,
			EMPTY_SEQUENCE);
    if (reqpost == NULL) {
      print_result(result,request);
      Result_free(&result);
    } else {
#ifdef HAVE_PTHREAD
      Reqpost_put_result(reqpost,result);
#endif
    }
      
  } else if (Sequence_fulllength_given(queryseq) < 
#ifdef PMAP
	     INDEX1PART_AA
#else
	     INDEX1PART
#endif
	     ) {
    result = Result_new(jobid,worker_id,(Chimera_T) NULL,(Stage3_T *) NULL,0,/*diagnostic*/NULL,
			SHORT_SEQUENCE);
    if (reqpost == NULL) {
      print_result(result,request);
      Result_free(&result);
    } else {
#ifdef HAVE_PTHREAD
      Reqpost_put_result(reqpost,result);
#endif
    }

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
      result = Result_new(jobid,worker_id,(Chimera_T) NULL,(Stage3_T *) NULL,0,diagnostic,POOR_SEQUENCE);
      if (reqpost == NULL) {
	print_result(result,request);
	Result_free(&result);
      } else {
#ifdef HAVE_PTHREAD
	debug(printf("worker_thread %d: Posting result\n",worker_id));
	Reqpost_put_result(reqpost,result);
#endif
      }
    } else if (repetitivep == true && prune_repetitive_p == true) {
      result = Result_new(jobid,worker_id,(Chimera_T) NULL,(Stage3_T *) NULL,0,diagnostic,REPETITIVE);
      if (reqpost == NULL) {
	print_result(result,request);
	Result_free(&result);
      } else {
#ifdef HAVE_PTHREAD
	debug(printf("worker_thread %d: Posting result\n",worker_id));
	Reqpost_put_result(reqpost,result);
#endif
      }

    } else 
#endif	/* PMAP */
      if (usersegment != NULL && userstage1p == false) {
#ifndef PMAP
#if 0
	/* Don't do Sequence_trim, because it affects sequences like NM_018406 */
      Sequence_trim(queryseq,diagnostic->query_trim_start,diagnostic->query_trim_end);
      Sequence_trim(queryuc,diagnostic->query_trim_start,diagnostic->query_trim_end);
#endif
#endif
      stage3array = stage3_from_usersegment(&npaths,/*lowidentityp*/false,queryseq,queryuc,usersegment,
					    oligoindices_major,noligoindices_major,oligoindices_minor,noligoindices_minor,
					    pairpool,intpool,diagpool,gbuffer,
					    dynprogL,dynprogM,dynprogR,stopwatch);
      result = Result_new(jobid,worker_id,(Chimera_T) NULL,stage3array,npaths,diagnostic,NO_FAILURE);
      if (reqpost == NULL) {
	print_result(result,request);
	Result_free(&result);
      } else {
#ifdef HAVE_PTHREAD
	debug(printf("worker_thread %d: Posting result\n",worker_id));
	Reqpost_put_result(reqpost,result);
#endif	/* HAVE_PTHREAD */
      }

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
				maxintronlen_bound,maxtotallen_bound,stuttercycles,stutterhits,
				/*subsequencep*/false,genome,genome_totallength,diagnostic,stopwatch);

      stage3array = apply_stage3(&chimera,&npaths,gregions,lowidentityp,queryseq,queryuc,usersegment,
				 oligoindices_major,noligoindices_major,oligoindices_minor,noligoindices_minor,
				 matchpool,pairpool,intpool,diagpool,
				 dynprogL,dynprogM,dynprogR,gbuffer,stopwatch);

      result = Result_new(jobid,worker_id,chimera,stage3array,npaths,diagnostic,NO_FAILURE);
      if (reqpost == NULL) {
	print_result(result,request);
	Result_free(&result);
      } else {
#ifdef HAVE_PTHREAD
	debug(printf("worker_thread %d: Posting result\n",worker_id));
	Reqpost_put_result(reqpost,result);
#endif
      }

      Oligoindex_clear_inquery(oligoindices_major[0]);

    } /* Matches not user segment and not maponly */

    Sequence_free(&queryuc);
  } /* Matches sequence length > 0 */

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
	       Sequence_T queryseq, Sequence_T usersegment) {
  Oligoindex_T *oligoindices_major, *oligoindices_minor;
  int noligoindices_major, noligoindices_minor;
  Dynprog_T dynprogL, dynprogM, dynprogR;
  Gbuffer_T gbuffer;		/* Individual space for reading genome sequence,
				   revcomp, and strain calculations */
  Matchpool_T matchpool;
  Pairpool_T pairpool;
  Intpool_T intpool;
  Diagpool_T diagpool;
  Stopwatch_T stopwatch;
  Request_T request;
  int jobid = 0;

  oligoindices_major = Oligoindex_new_major(&noligoindices_major);
  oligoindices_minor = Oligoindex_new_minor(&noligoindices_minor);
  dynprogL = Dynprog_new(nullgap,EXTRAQUERYGAP,maxpeelback,extramaterial_end,extramaterial_paired);
  dynprogM = Dynprog_new(nullgap,EXTRAQUERYGAP,maxpeelback,extramaterial_end,extramaterial_paired);
  dynprogR = Dynprog_new(nullgap,EXTRAQUERYGAP,maxpeelback,extramaterial_end,extramaterial_paired);
  gbuffer = Gbuffer_new();
  matchpool = Matchpool_new();
  pairpool = Pairpool_new();
  intpool = Intpool_new();
  diagpool = Diagpool_new();
  stopwatch = diagnosticp == true ? Stopwatch_new() : (Stopwatch_T) NULL;

  while (jobid == 0 || (queryseq = Sequence_read_multifile(&nextchar,&input,&files,&nfiles,maponlyp)) != NULL) {
    request = Request_new(jobid++,queryseq,usersegment);
TRY
  handle_request(request,/*worker_id*/0,matchpool,pairpool,intpool,diagpool,
		 oligoindices_major,noligoindices_major,oligoindices_minor,noligoindices_minor,
		 gbuffer,dynprogL,dynprogM,dynprogR,(Reqpost_T) NULL,stopwatch);
ELSE
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
    Request_free(&request);
   /* Don't free queryseq; done by Request_free */
  }

  Stopwatch_free(&stopwatch);
  Diagpool_free(&diagpool);
  Intpool_free(&intpool);
  Pairpool_free(&pairpool);
  Matchpool_free(&matchpool);
  Gbuffer_free(&gbuffer);
  Dynprog_free(&dynprogR);
  Dynprog_free(&dynprogM);
  Dynprog_free(&dynprogL);
  Oligoindex_free_array(&oligoindices_minor,noligoindices_minor);
  Oligoindex_free_array(&oligoindices_major,noligoindices_major);

  return;
}


#ifdef HAVE_PTHREAD
Blackboard_T blackboard_global;	/* Needed only for exception handling */

static void *
worker_thread (void *data) {
  Oligoindex_T *oligoindices_major, *oligoindices_minor;
  int noligoindices_major, noligoindices_minor;
  Dynprog_T dynprogL, dynprogM, dynprogR;
  Gbuffer_T gbuffer;		/* Individual space for reading genome sequence,
				   revcomp, and strain calculations */
  Matchpool_T matchpool;
  Pairpool_T pairpool;
  Intpool_T intpool;
  Diagpool_T diagpool;
  Stopwatch_T stopwatch;
  Request_T request;
  Sequence_T queryseq;

  Reqpost_T reqpost = (Reqpost_T) data;
  int worker_id = Reqpost_id(reqpost), id;

  /* Thread-specific data and storage */
  oligoindices_major = Oligoindex_new_major(&noligoindices_major);
  oligoindices_minor = Oligoindex_new_minor(&noligoindices_minor);
  dynprogL = Dynprog_new(nullgap,EXTRAQUERYGAP,maxpeelback,extramaterial_end,extramaterial_paired);
  dynprogM = Dynprog_new(nullgap,EXTRAQUERYGAP,maxpeelback,extramaterial_end,extramaterial_paired);
  dynprogR = Dynprog_new(nullgap,EXTRAQUERYGAP,maxpeelback,extramaterial_end,extramaterial_paired);
  gbuffer = Gbuffer_new();
  matchpool = Matchpool_new();
  pairpool = Pairpool_new();
  intpool = Intpool_new();
  diagpool = Diagpool_new();
  stopwatch = diagnosticp == true ? Stopwatch_new() : (Stopwatch_T) NULL;
  Except_stack_create();

  while ((request = Reqpost_get_request(reqpost))) {
TRY
    handle_request(request,worker_id,matchpool,pairpool,intpool,diagpool,
		   oligoindices_major,noligoindices_major,oligoindices_minor,noligoindices_minor,
		   gbuffer,dynprogL,dynprogM,dynprogR,reqpost,stopwatch);
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
	queryseq = Request_queryseq(request);
	if (queryseq == NULL) {
	  fprintf(stderr,"NULL");
	} else if (Sequence_accession(queryseq) == NULL) {
	  fprintf(stderr,"unnamed (%d bp)",Sequence_fulllength_given(queryseq));
	} else {
	  fprintf(stderr,"%s (%d bp)",Sequence_accession(queryseq),Sequence_fulllength_given(queryseq));
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
  Diagpool_free(&diagpool);
  Intpool_free(&intpool);
  Pairpool_free(&pairpool);
  Matchpool_free(&matchpool);
  Gbuffer_free(&gbuffer);
  Dynprog_free(&dynprogR);
  Dynprog_free(&dynprogM);
  Dynprog_free(&dynprogL);
  Oligoindex_free_array(&oligoindices_minor,noligoindices_minor);
  Oligoindex_free_array(&oligoindices_major,noligoindices_major);
  Reqpost_free(&reqpost);

  return (void *) NULL;
}
#endif


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
  Gbuffer_T gbuffer;		/* Individual space for reading genome sequence,
				   revcomp, and strain calculations */
  Matchpool_T matchpool;
  Pairpool_T pairpool;
  Intpool_T intpool;
  Diagpool_T diagpool;
  Stopwatch_T stopwatch;

  Genomicpos_T genomicstart, genomiclength;
  Sequence_T genomicseg, queryuc, referenceuc;
  int jobid = 0;

  Chimera_T chimera = NULL;
  List_T gregions;
  Stage3_T *stage3array, stage3, stage3ref;
  int npaths, i;

  oligoindices_major = Oligoindex_new_major(&noligoindices_major);
  oligoindices_minor = Oligoindex_new_minor(&noligoindices_minor);
  dynprogL = Dynprog_new(nullgap,EXTRAQUERYGAP,maxpeelback,extramaterial_end,extramaterial_paired);
  dynprogM = Dynprog_new(nullgap,EXTRAQUERYGAP,maxpeelback,extramaterial_end,extramaterial_paired);
  dynprogR = Dynprog_new(nullgap,EXTRAQUERYGAP,maxpeelback,extramaterial_end,extramaterial_paired);
  gbuffer = Gbuffer_new();
  matchpool = Matchpool_new();
  pairpool = Pairpool_new();
  intpool = Intpool_new();
  diagpool = Diagpool_new();
  stopwatch = diagnosticp == true ? Stopwatch_new() : (Stopwatch_T) NULL;

  Matchpool_reset(matchpool);
  Pairpool_reset(pairpool);

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
			    maxintronlen_bound,maxtotallen_bound,stuttercycles,stutterhits,/*subsequencep*/false,
			    genome,genome_totallength,diagnostic,/*stopwatch*/NULL);
  stage3array = apply_stage3(&chimera,&npaths,gregions,lowidentityp,referenceseq,referenceuc,/*usersegment*/NULL,
			     oligoindices_major,noligoindices_major,oligoindices_minor,noligoindices_minor,
			     matchpool,pairpool,intpool,diagpool,dynprogL,dynprogM,dynprogR,gbuffer,stopwatch);
  Diagnostic_free(&diagnostic);

  /* chimera should be NULL */
  for (i = 1; i < npaths; i++) {
    stage3 = stage3array[i];
    Stage3_free(&stage3);
  }
  if (npaths > 0) {
    stage3ref = stage3array[0];
#ifdef PMAP
    stage3ref = Stage3_translate_cdna(stage3ref,queryseq,strictp);
    stage3ref = Stage3_backtranslate_cdna(stage3ref,/*diagnosticp*/false);
#else
    stage3ref = Stage3_translate_genomic(stage3ref,/*fulllengthp*/true,/*truncatep*/false,strictp);
#endif
    FREE(stage3array);

    Stage3_genomicbounds(&genomicstart,&genomiclength,stage3ref);
    Gbuffer_alloc_contents(gbuffer,genomiclength);
    if (genomealt != NULL) {
      genomicseg = Genome_get_segment(genomealt,genomicstart,genomiclength,chromosome_iit,/*revcomp*/false,
				      Gbuffer_chars1(gbuffer),Gbuffer_chars2(gbuffer),Gbuffer_gbufferlen(gbuffer));
    } else {
      genomicseg = Genome_get_segment(genome,genomicstart,genomiclength,chromosome_iit,/*revcomp*/false,
				      Gbuffer_chars1(gbuffer),Gbuffer_chars2(gbuffer),Gbuffer_gbufferlen(gbuffer));
    }

    while (jobid == 0 || (queryseq = Sequence_read_multifile(&nextchar,&input,&files,&nfiles,maponlyp)) != NULL) {
      Matchpool_reset(matchpool);
      Pairpool_reset(pairpool);

      printf(">");
      Sequence_print_header(queryseq,checksump);
      diagnostic = Diagnostic_new();
      if (Sequence_fulllength_given(queryseq) <= 0) {
#ifdef BETATEST
	print_npaths(0,diagnostic,/*usersegment*/NULL,chrsubset,/*chimera*/NULL,EMPTY_SEQUENCE);
#else
	print_npaths(0,diagnostic,/*usersegment*/NULL,chrsubset,EMPTY_SEQUENCE);
#endif

      } else if (Sequence_fulllength_given(queryseq) <
#ifdef PMAP
		 INDEX1PART_AA
#else
		 INDEX1PART
#endif
		 ) {
#ifdef BETATEST
	print_npaths(0,diagnostic,/*usersegment*/NULL,chrsubset,/*chimera*/NULL,SHORT_SEQUENCE);
#else
	print_npaths(0,diagnostic,/*usersegment*/NULL,chrsubset,SHORT_SEQUENCE);
#endif

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
#ifdef BETATEST
	  print_npaths(0,diagnostic,/*usersegment*/NULL,chrsubset,/*chimera*/NULL,POOR_SEQUENCE);
#else
	  print_npaths(0,diagnostic,/*usersegment*/NULL,chrsubset,POOR_SEQUENCE);
#endif
	} else if (repetitivep == true && prune_repetitive_p == true) {
#ifdef BETATEST
	  print_npaths(0,diagnostic,/*usersegment*/NULL,chrsubset,/*chimera*/NULL,REPETITIVE);
#else
	  print_npaths(0,diagnostic,/*usersegment*/NULL,chrsubset,REPETITIVE);
#endif
	} else {
#endif /* PMAP */
	  stage3array = stage3_from_usersegment(&npaths,lowidentityp,queryseq,queryuc,genomicseg,
						oligoindices_major,noligoindices_major,oligoindices_minor,noligoindices_minor,
						pairpool,intpool,diagpool,gbuffer,dynprogL,dynprogM,dynprogR,stopwatch);
	  
	  if (npaths == 0) {
#ifdef BETATEST
	    print_npaths(0,diagnostic,/*usersegment*/NULL,chrsubset,/*chimera*/NULL,NO_FAILURE);
#else
	    print_npaths(0,diagnostic,/*usersegment*/NULL,chrsubset,NO_FAILURE);
#endif
	  } else if (print_coordinates_p == true) {
	    Stage3_print_coordinates(stage3array[0],queryseq,chromosome_iit,
				     /*zerobasedp*/false,invertmode,fulllengthp,truncatep,strictp,maponlyp);

	  } else {
	    /* Usual output */
#ifdef BETATEST
	    print_npaths(1,diagnostic,/*usersegment*/NULL,chrsubset,/*chimera*/NULL,NO_FAILURE);
#else
	    print_npaths(1,diagnostic,/*usersegment*/NULL,chrsubset,NO_FAILURE);
#endif
#ifndef PMAP
	    Stage3_translate_cdna_via_reference(stage3array[0],stage3ref,literalrefp);
#endif
	    Stage3_fix_cdna_direction(stage3array[0],stage3ref);
	    Stage3_print_mutations(stage3array[0],stage3ref,chromosome_iit,queryseq,
				   dbversion,showalignp,/*zerobasedp*/false,continuousp,diagnosticp,proteinmode,
				   invertmode,nointronlenp,wraplength,/*maxmutations*/1000000);
	    for (i = 0; i < npaths; i++) {
	      stage3 = stage3array[i];
	      Stage3_free(&stage3);
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
  Intpool_free(&intpool);
  Pairpool_free(&pairpool);
  Gbuffer_free(&gbuffer);
  Dynprog_free(&dynprogR);
  Dynprog_free(&dynprogM);
  Dynprog_free(&dynprogL);
  Oligoindex_free_array(&oligoindices_minor,noligoindices_minor);
  Oligoindex_free_array(&oligoindices_major,noligoindices_major);

  return;
}


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
parse_batch (int *batch_modulus, int *batch_interval, char *string) {
  char *p = string;

  if (sscanf(p,"%d",&(*batch_modulus)) < 1) {
    fprintf(stderr,"Cannot parse first integer from %s\n",string);
    exit(9);
  }

  while (*p != '\0' && isdigit(*p)) {
    p++;
  }
  while (*p != '\0' && !isdigit(*p)) {
    p++;
  }
  if (sscanf(p,"%d",&(*batch_interval)) < 1) {
    fprintf(stderr,"Cannot parse first integer from %s\n",string);
    exit(9);
  }
  if ((*batch_modulus) >= (*batch_interval)) {
    fprintf(stderr,"In %s, batch number %d must be less than the number of batches %d\n",
	    string,*batch_modulus,*batch_interval);
    exit(9);
  }
  if (*batch_interval == 0) {
    fprintf(stderr,"Bad batch specification %s.  Batch interval cannot be 0.\n",string);
    exit(9);
  }

  return;
}

int
main (int argc, char *argv[]) {
  Sequence_T queryseq, usersegment = NULL, referenceseq = NULL;
  char *genomesubdir = NULL, *mapdir = NULL, *iitfile = NULL, *fileroot = NULL, **argstart;
  FILE *input;

  int user_ngap = -1;
  bool showcontigp = true, multiple_sequences_p = false;
  char **files;
  int nfiles, nextchar = '\0';

#ifdef HAVE_PTHREAD
  int ret, i;
  pthread_attr_t thread_attr_detach, thread_attr_join;
  Blackboard_T blackboard;
  Request_T request;
#endif

  int opt, len, c;
  extern int optind;
  extern char *optarg;
  int long_option_index = 0;

#ifdef HAVE_SIGACTION
  struct sigaction signal_action;
#endif

  while ((opt = getopt_long(argc,argv,
#ifdef PMAP
			    "q:D:d:Gg:C:B:K:L:x:1w:t:sc:H:X:SA03468:9n:f:ZO5o:V:M:m:ebu:E:PQYNI:i:l:v?",
#else
			    "q:D:d:Gg:C:B:K:L:x:1w:t:sc:H:X:p:SA03468:9n:f:ZO5o:V:M:m:ebu:E:PQFTz:YNI:i:l:v?",
#endif
			    long_options, &long_option_index)) != -1) {
    switch (opt) {
    case 'q': parse_batch(&batch_modulus,&batch_interval,optarg); break;
    case 'D': user_genomedir = optarg; break;
    case 'd': 
      dbroot = (char *) CALLOC(strlen(optarg)+1,sizeof(char));
      strcpy(dbroot,optarg);
      break;
    case 'G': uncompressedp = true; break;
    case 'g': user_genomicseg = optarg; break;
    case 'C': user_chrsubsetfile = optarg; break;

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
    case 'K': maxintronlen_bound = atoi(optarg); break;
    case 'L': maxtotallen_bound = atoi(optarg); break;
    case 'x': 
#ifdef PMAP
      chimera_margin = atoi(optarg)/3; 
#else
      chimera_margin = atoi(optarg); 
#endif
      if (chimera_margin < CHIMERA_SLOP) {
	chimera_margin = CHIMERA_SLOP;
      }
      break;
    case '1': maponlyp = true; break;
    case 'w': referencefile = optarg; break;
#ifdef HAVE_PTHREAD
    case 't': nworkers = atoi(optarg); break;
#endif
    case 's': altstrainp = true; break;
    case 'c': user_chrsubsetname = optarg; break;
    case 'H': minendexon = atoi(optarg); break;
    case 'X': 
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

#ifndef PMAP
    case 'p': switch (atoi(optarg)) {
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

    case 'S': summaryonlyp = true; break;
    case 'A': showalignp = true; break;
    case '0': exception_raise_p = false; break; /* Allows signals to pass through */
    case '3': continuousp = true; break;
    case '4': continuous_by_exon_p = true; showalignp = true; break;
    case '6': debug_graphic_p = true; diagnosticp = false; break;
    case '8': diagnosticp = true;
      if (!strcmp(optarg,"stage2")) {
	stage3debug = POST_STAGE2;
      } else if (!strcmp(optarg,"smoothing")) {
	stage3debug = POST_SMOOTHING;
      } else if (!strcmp(optarg,"singles")) {
	stage3debug = POST_SINGLES;
      } else if (!strcmp(optarg,"introns")) {
	stage3debug = POST_INTRONS;
      } else if (!strcmp(optarg,"canonical")) {
	stage3debug = POST_CANONICAL;
      } else if (!strcmp(optarg,"changepoint")) {
	stage3debug = POST_CHANGEPOINT;
      } else if (!strcmp(optarg,"distalmedial")) {
	stage3debug = POST_DISTAL_MEDIAL;
      } else if (!strcmp(optarg,"trimmiddle")) {
	stage3debug = POST_TRIM_MIDDLE;
      } else {
	fprintf(stderr,"Allowed arguments for -8 flag are stage2, smoothing, singles, introns, canonical, changepoint, distalmedial, trimmiddle\n");
	exit(9);
      }
      break;
    case '9': checkp = true; diagnosticp = true; break;
    case 'n': maxpaths = atoi(optarg); break;
    case 'f':
#ifdef PMAP
      if (!strcmp(optarg,"0")) {
	psloutput_pro_p = true;
      } else if (!strcmp(optarg,"1")) {
	psloutput_nt_p = true;
      } else if (!strcmp(optarg,"2")) {
	gffoutput_p = true; gff_gene_format_p = true;
      } else if (!strcmp(optarg,"3")) {
	gffoutput_p = true; gff_gene_format_p = false;
      } else if (!strcmp(optarg,"6")) {
	splicesites_output_p = true;
      } else if (!strcmp(optarg,"7")) {
	exon_mapoutput_p = true;
      } else if (!strcmp(optarg,"8")) {
	mapoutput_p = true;
      } else if (!strcmp(optarg,"9")) {
	print_coordinates_p = true;
      } else {
	fprintf(stderr,"Output format %s not recognized\n",optarg);
	exit(9);
      }
#else
      if (!strcmp(optarg,"1")) {
	psloutput_nt_p = true;
      } else if (!strcmp(optarg,"2")) {
	gffoutput_p = true; gff_gene_format_p = true;
      } else if (!strcmp(optarg,"3")) {
	gffoutput_p = true; gff_gene_format_p = false;
      } else if (!strcmp(optarg,"6")) {
	splicesites_output_p = true;
      } else if (!strcmp(optarg,"7")) {
	exon_mapoutput_p = true;
      } else if (!strcmp(optarg,"8")) {
	mapoutput_p = true;
      } else if (!strcmp(optarg,"9")) {
	print_coordinates_p = true;
      } else {
	fprintf(stderr,"Output format %s not recognized\n",optarg);
	exit(9);
      }
#endif
      break;
    case 'Z': compressoutp = true; break;
    case 'O': orderedp = true; break;
    case '5': checksump = true; break;
    case 'o': chimera_overlap = atoi(optarg); break;
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
    case 'u': nflanking = atoi(optarg); break;

    case 'E': 
      if (!strcmp(optarg,"cdna")) {
	cdna_exons_p = true;
      } else if (!strncmp(optarg,"gen",3)) {
	genomic_exons_p = true;
      } else {
	fprintf(stderr,"Argument to -E flag must be either \"cdna\" or \"genomic\"\n");
	exit(9);
      }
      break;

#ifdef PMAP
    case 'P': protein_genomic_p = true; break;
    case 'Q': print_cdna_p = true; break; 
#else
    case 'P': print_cdna_p = true; break;
    case 'Q': protein_genomic_p = true; break;
    case 'F': fulllengthp = true; break;
    case 'T': truncatep = true; fulllengthp = true; break;
    case 'z': 
      if (!strcmp(optarg,"sense")) {
	sense = +1;
      } else if (!strcmp(optarg,"antisense")) {
	sense = -1;
      } else if (!strcmp(optarg,"auto")) {
	sense = 0;
      } else {
	fprintf(stderr,"direction %s not recognized.  Must be sense, antisense, or auto\n",optarg);
	exit(9);
      }

#endif
    case 'Y': strictp = false; break;
    case 'N': nointronlenp = true; break;
    case 'I': invertmode = atoi(optarg); break;
    case 'i': user_ngap = atoi(optarg); break;
    case 'l': wraplength = atoi(optarg); break;

    case 'v': print_program_version(); exit(0);
    case '?': print_program_usage(); exit(0);
    default: exit(9);
    }
  }
  argc -= optind;
  argv += optind;
      
  if (user_genomicseg != NULL) {
    /* Ignore -D and -d flags */
  } else if (dbroot == NULL) {
    fprintf(stderr,"Need to specify the -d or the -g flag\n");
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

      if (IIT_ntypes(map_iit) == 3 &&
	  ((map_iit_forward_type = IIT_typeint(map_iit,"FWD")) >= 0 ||
	   (map_iit_forward_type = IIT_typeint(map_iit,"+")) >= 0) &&
	  ((map_iit_reverse_type = IIT_typeint(map_iit,"REV")) >= 0 ||
	   (map_iit_reverse_type = IIT_typeint(map_iit,"-")) >= 0)) {
	map_iit_universal_p = true;
      } else {
	map_iit_universal_p = false;
	check_map_iit(map_iit,chromosome_iit);
      }
      FREE(iitfile);
      FREE(mapdir);
      FREE(map_iitfile);
    }
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

  if (user_ngap >= 0) {
    ngap = user_ngap;
  } else if (cdna_exons_p == true || genomic_exons_p == true) {
    /* If user didn't specify, then set to zero */
    ngap = 0;
  };

  if (maxintronlen_bound > maxtotallen_bound) {
    maxintronlen_bound = maxtotallen_bound;
  }

  /* Read user segment before rest of sequences, because of shared usage of sequence.c */
  if (user_genomicseg != NULL) {
    if ((input = FOPEN_READ_TEXT(user_genomicseg)) == NULL) {
      fprintf(stderr,"Can't open file %s\n",user_genomicseg);
      exit(9);
    }
    if ((usersegment = Sequence_read_unlimited(input)) == NULL) {
      fprintf(stderr,"File %s is empty\n",user_genomicseg);
      exit(9);
    }
    fclose(input);
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
  if ((queryseq = Sequence_read_multifile(&nextchar,&input,&files,&nfiles,maponlyp)) == NULL) {
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
#ifndef PMAP
    batch_offsets_p = true;	/* Results in pre-reading of offsets file */
#endif
  } else {
    /* multiple_sequences_p = false; */
    if (batch_offsets_p == true || batch_positions_p == true || batch_genome_p == true) {
      /* fprintf(stderr,"Note: only 1 sequence detected.  Ignoring batch (-B) command\n"); */
      batch_offsets_p = false;
      batch_positions_p = false;
      batch_genome_p = false;
    }
  }

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
      indexdb_fwd = Indexdb_new_segment(Sequence_fullpointer(usersegment),/*watsonp*/true,index1interval);
      indexdb_rev = Indexdb_new_segment(Sequence_fullpointer(usersegment),/*watsonp*/false,index1interval);
#else
      indexdb = Indexdb_new_segment(Sequence_fullpointer(usersegment),index1interval);
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
    indexdb_fwd = Indexdb_new_genome(genomesubdir,fileroot,FWD_FILESUFFIX,/*snps_root*/NULL,
				     /*required_interval*/0,batch_offsets_p,batch_positions_p);
    indexdb_rev = Indexdb_new_genome(genomesubdir,fileroot,REV_FILESUFFIX,/*snps_root*/NULL,
				     /*required_interval*/0,batch_offsets_p,batch_positions_p);

    if (indexdb_fwd == NULL || indexdb_rev == NULL) {
      fprintf(stderr,"Cannot find offsets file %s.%s*offsets or %s.%s*offsets.\n",
	      fileroot,FWD_FILESUFFIX,fileroot,REV_FILESUFFIX);
      fprintf(stderr,"You may need to run 'pmapindex -d %s' to build the indices needed for PMAP.\n",
	      fileroot);
      exit(9);
    }      
#else
    if ((indexdb = Indexdb_new_genome(genomesubdir,fileroot,IDX_FILESUFFIX,/*snps_root*/NULL,
				      /*required_interval*/0,batch_offsets_p,batch_positions_p)) == NULL) {
      /* Try older version */
      if ((indexdb = Indexdb_new_genome(genomesubdir,fileroot,"id",/*snps_root*/NULL,
					/*required_interval*/0,batch_offsets_p,batch_positions_p)) == NULL) {
	fprintf(stderr,"Cannot find offsets file %s.%s*offsets or %s.%s*offsets\n",fileroot,IDX_FILESUFFIX,fileroot,"id");
	exit(9);
      }
    }
    indexdb_size_threshold = (int) (10*Indexdb_mean_size(indexdb,/*cmetp*/false));
    debug(printf("Size threshold is %d\n",indexdb_size_threshold));
#endif
    genome = Genome_new(genomesubdir,fileroot,/*snps_root*/NULL,uncompressedp,batch_genome_p);
    if (snps_root != NULL) {
      genomealt = Genome_new(genomesubdir,fileroot,snps_root,uncompressedp,batch_genome_p);
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

  /* Output headers */
  if (gffoutput_p == true) {
    printf("##gff-version   3\n");
    printf("# Generated by GMAP version %s using call: ",PACKAGE_VERSION);
    argstart = &(argv[-optind]);
    for (c = 0; c < argc + optind; c++) {
      printf(" %s",argstart[c]);
    }
    printf("\n");
  }

  if (referenceseq != NULL) {
    chimera_margin = -1;
    align_relative(input,files,nfiles,nextchar,queryseq,referenceseq);
    Sequence_free(&referenceseq);

  } else {
#ifndef HAVE_PTHREAD
    single_thread(input,files,nfiles,nextchar,queryseq,usersegment);
#else
    if (multiple_sequences_p == false) {
      single_thread(input,files,nfiles,nextchar,queryseq,usersegment);
    } else if (nworkers == 0) {
      single_thread(input,files,nfiles,nextchar,queryseq,usersegment);
    } else {
      /* Make blackboard and threads */
      blackboard = Blackboard_new(input,files,nfiles,nextchar,usersegment,nworkers);
      blackboard_global = blackboard; /* Needed only for exception handling */
      if (batch_modulus != 0) {
	Sequence_free(&queryseq);
      } else {
	debug(printf("input_thread: Putting request id 0\n"));
	request = Request_new(requestid++,queryseq,usersegment);
	Blackboard_put_request(blackboard,request);
      }

      pthread_attr_init(&thread_attr_detach);
      if ((ret = pthread_attr_setdetachstate(&thread_attr_detach,PTHREAD_CREATE_DETACHED)) != 0) {
	fprintf(stderr,"ERROR: pthread_attr_setdetachstate %d\n",ret);
	exit(1);
      }
      pthread_attr_init(&thread_attr_join);
      if ((ret = pthread_attr_setdetachstate(&thread_attr_join,PTHREAD_CREATE_JOINABLE)) != 0) {
	fprintf(stderr,"ERROR: pthread_attr_setdetachstate %d\n",ret);
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

      if (!Blackboard_donep(blackboard)) {
	fprintf(stderr,"Something appears to be wrong with your pthread library; not all input was processed.\n");
	exit(9);
      }

      /* Do not delete global_except_key, because worker threads might still need it */
      /* Except_term_pthread(); */

      FREE(worker_thread_ids);
      Blackboard_free(&blackboard);
    }

#endif /* HAVE_PTHREAD */
  }

  if (usersegment != NULL) {
    Sequence_free(&usersegment);
  }

#ifdef PMAP
  Backtranslation_term();
#endif
  Dynprog_term();

  if (argc > 0) {
    fclose(input);
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
