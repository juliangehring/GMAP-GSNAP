static char rcsid[] = "$Id: gmap.c,v 1.321 2006/04/07 01:18:37 twu Exp $";
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

#include "except.h"
#include "mem.h"
#include "bool.h"
#include "fopen.h"

#include "sequence.h"
#include "oligoindex.h"
#include "match.h"
#include "matchpool.h"
#include "matchpair.h"
#include "pairpool.h"
#include "genome.h"
#include "stage1.h"
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
#include "params.h"
#include "blackboard.h"
#include "intlist.h"
#include "list.h"
#include "listdef.h"
#include "iit-read.h"
#include "chrsubset.h"
#include "datadir.h"
#include "getopt.h"

#define POSSIBLE_OLIGOS 65536	/* 4^8 */
#define MAX_OLIGO_DEPTH 2.0
#define MAX_BADOLIGOS 1.0	/* Setting to 1.0 effectively turns this check off */

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

/* Sequence pruning */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
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
 *   Program options
 ************************************************************************/

/* Input options */
static char *user_genomedir = NULL;
static char *dbroot = NULL;
static char *releasestring = NULL;
static char *dbversion = NULL;
static char *user_genomicseg = NULL;
static char *user_chrsubsetfile = NULL;


/* Compute options */
static bool batch_offsets_p = false;
static bool batch_positions_p = false;
static bool batch_genome_p = false;
static int maxtotallen_bound = 1200000;
static int maxintronlen_bound = 1200000;
static int chimera_margin = -1;
static bool maponlyp = false;
#ifdef PMAP
static bool userstage1p = true; /* Apply stage 1 for user-provided genomic segments */
#else
static bool userstage1p = false; /* Apply stage 1 for user-provided genomic segments */
#endif
static int index1interval = INDEX1INTERVAL; /* Stage 1 interval if user provides a genomic segment */
static bool crossspeciesp = false;
static char *mutationref = NULL;
static bool literalrefp = false;
static bool altstrainp = false;
#ifdef HAVE_PTHREAD
static int nworkers = 1;	/* (int) sysconf(_SC_NPROCESSORS_ONLN) */
#endif
static bool extend_mismatch_p = false;
static bool end_microexons_p = false;
static char *user_chrsubsetname = NULL;


/* Output options */
static bool summaryonlyp = false;
static bool showalignp = false;
static bool continuousp = false;
static bool diagnosticp = false;
static int maxpaths = 5;	/* 0 means 1 if nonchimeric, 2 if chimeric */
static bool psloutput_nt_p = false;
#ifdef PMAP
static bool psloutput_pro_p = false;
#endif
static bool print_coordinates_p = false;
static bool compressoutp = false;
static bool orderedp = false;
static bool checksump = false;
static char *user_mapdir = NULL;
static char *map_iitfile = NULL;
static bool map_exons_p = false;
static bool map_bothstrands_p = false;
static int chimera_overlap = 0;
static int maxmutations = 20;

/* Alignment options */
static bool cdna_exons_p = false;
static bool genomic_exons_p = false;
static bool print_cdna_p = false;
static bool protein_genomic_p = false;
#ifndef PMAP
static bool fulllengthp = false;
static bool truncatep = false;
#endif
static int proteinmode = 1;
static bool uncompressedp = false;
static bool nointronlenp = false;
static int invertmode = 0;
static int ngap = 3;
static int wraplength = 50;

#ifdef PMAP
/* Used alphabetically: 1359ABbCcDdEefGgIiKLlMmNnOoPQRSstUVvwXxZ */
#else
/* Used alphabetically: 1359ABbCcDdEeFfGgIiKLlMmNnOoPQRSsTtUVvwXxZ */
#endif

static struct option long_options[] = {
  /* Input options */
  {"dir", required_argument, 0, 'D'},	/* user_genomedir */
  {"db", required_argument, 0, 'd'}, /* dbroot */
  {"genomefull", no_argument, 0, 'G'}, /* uncompressedp */
  {"rel", required_argument, 0, 'R'}, /* releasestring */
  {"gseg", required_argument, 0, 'g'}, /* user_genomicseg */

  /* Compute options */
#ifdef HAVE_MMAP
  {"batch", required_argument, 0, 'B'}, /* batch_offsets_p, batch_positions_p, batch_genome_p */
#endif
  {"intronlength", required_argument, 0, 'K'}, /* maxintronlen_bound */
  {"totallength", required_argument, 0, 'L'}, /* maxtotallen_bound */
  {"chimera_margin", required_argument, 0, 'x'}, /* chimera_margin */
  {"maponly", no_argument, 0, '1'}, /* maponlyp */
  {"cross", no_argument, 0, 'X'}, /* crossspeciesp */
  {"mutationref", required_argument, 0, 'w'},	/* mutationref */
#ifdef HAVE_PTHREAD
  {"nthreads", required_argument, 0, 't'}, /* nworkers */
#endif
  {"altstrain", no_argument, 0, 's'},	/* altstrainp */
  {"chrsubsetfile", required_argument, 0, 'C'}, /* user_chrsubsetfile */
  {"chrsubset", required_argument, 0, 'c'}, /* user_chrsubsetname */
  {"endmicro", no_argument, 0, 'U'}, /* end_microexons_p */

  /* Output options */
  {"summary", no_argument, 0, 'S'}, /* summaryonlyp */
  {"align", no_argument, 0, 'A'}, /* showalignp */
  {"continuous", no_argument, 0, '3'}, /* continuousp */
  {"diagnostic", no_argument, 0, '9'}, /* diagnosticp */
  {"npaths", required_argument, 0, 'n'}, /* maxpaths */
  {"format", required_argument, 0, 'f'}, /* psloutput_nt_p, psloutput_pro_p, print_coordinates_p */
  {"compress", no_argument, 0, 'Z'}, /* compressoutp */
  {"ordered", no_argument, 0, 'O'}, /* orderedp */
  {"md5", no_argument, 0, '5'}, /* checksump */
  {"mapdir", required_argument, 0, 'M'}, /* user_mapdir */
  {"map", required_argument, 0, 'm'},	/* map_iitfile */
  {"mapexons", no_argument, 0, 'e'}, /* map_exons_p */
  {"mapboth", no_argument, 0, 'b'}, /* map_bothstrands_p */
  {"chimera_overlap", required_argument, 0, 'o'}, /* chimera_overlap */
  {"variations", required_argument, 0, 'v'}, /* maxmutations */

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
#endif
  {"nolengths", no_argument, 0, 'N'},	/* nointronlenp */
  {"invertmode", required_argument, 0, 'I'}, /* invertmode */
  {"introngap", required_argument, 0, 'i'}, /* ngap */
  {"wraplength", required_argument, 0, 'l'}, /* wraplength */
  
  /* Help options */
  {"version", no_argument, 0, 'V'}, /* print_program_version */
  {"help", no_argument, 0, '?'}, /* print_program_usage */
  {0, 0, 0, 0}
};


static void
print_program_version () {
  fprintf(stdout,"\n");
#ifdef PMAP
  fprintf(stdout,"PMAP: Protein Mapping and Alignment Program\n");
#else
  fprintf(stdout,"GMAP: Genomic Mapping and Alignment Program\n");
#endif
  fprintf(stdout,"Part of GMAP package, version %s\n",PACKAGE_VERSION);
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
  fprintf(stdout,"bigendian\n");
#else
  fprintf(stdout,"littleendian\n");
#endif
  fprintf(stdout,"Sizes: off_t (%d), size_t (%d), unsigned int (%d), long int (%d)\n",
	  sizeof(off_t),sizeof(size_t),sizeof(unsigned int),sizeof(long int));
  fprintf(stdout,"Default gmap directory: %s\n",GMAPDB);
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
    fprintf(stdout,"\
Usage: gmap [OPTIONS...] <FASTA file>, or\n\
       cat <FASTA file> | gmap [OPTIONS...]\n\
\n\
Input options (must include -d or -g)\n\
  -D, --dir=directory            Genome directory\n\
  -d, --db=STRING                Genome database\n\
  -G, --genomefull               Use full genome (all ASCII chars allowed;\n\
                                   built explicitly during setup), not\n\
                                   compressed version\n\
  -g, --gseg=filename            User-suppled genomic segment\n\
\n\
");
    fprintf(stdout,"Computation options\n");
#ifdef HAVE_MMAP
    fprintf(stdout,"\
  -B, --batch=INT                Batch mode (1 = pre-load only indices;\n\
                                   2 = pre-load both indices and genome)\n\
");
#endif
    fprintf(stdout,"\
  -K, --intronlength=INT         Max length for one intron (default 1200000)\n\
  -L, --totallength=INT          Max total intron length (default 1200000)\n\
  -x, --chimera_margin=INT       Amount of unaligned sequence that triggers\n\
                                   search for a chimera (default off)\n\
  -1, --maponly                  Perform only genomic mapping and alignment at ends\n\
                                   (skip internal alignment)\n\
  -X, --cross                    Cross-species mode (changes parameters)\n\
  -w, --mutationref=filename     Mutation reference\n\
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
  -U, --endmicro                 Allow end microexons with arbitrarily long introns\n\
\n\
");
    fprintf(stdout,"\
Output options\n\
  -S, --summary                  Show summary of alignments only\n\
  -A, --align                    Show alignments\n\
  -3, --continuous               Show alignment in three continuous lines\n\
  -9, --diagnostic               Show diagnostics (alignment with *'s for\n\
                                   stage 3)\n\
  -n, --npaths=INT               Maximum number of paths to show.  If set to 0,\n \
                                 prints two paths if chimera detected, else one.\n\
  -Z, --compress                 Print output in compressed format\n\
  -O, --ordered                  Print output in same order as input (relevant\n\
                                   only if there is more than worker thread)\n\
  -5, --md5                      Print MD5 checksum for each query sequence\n\
  -M, --mapdir=directory         Map directory\n\
  -m, --map=iitfile              Map file\n\
  -e, --mapexons                 Map each exon separately\n\
  -b, --mapboth                  Map both strands of genome\n\
  -o, --chimera_overlap          Overlap to show, if any, at chimera breakpoint\n\
  -v, --variations=INT           Maximum number of variations to show (default 20)\n\
");
#ifdef PMAP
    fprintf(stdout,"\
  -f, --format=INT               Format for output (0 = PSL format in protein coords,\n\
                                   1 = PSL format in nucleotide coords, 9 = coords\n\
                                   table format)\n\
");
#else
    fprintf(stdout,"\
  -f, --format=INT               Format for output (1 = PSL (BLAT) format, 9 = coords\n\
                                   in table format)\n\
");
#endif
    fprintf(stdout,"\n");
    fprintf(stdout,"\
Alignment output options\n\
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
  -F, --fulllength               Assume full-length protein, starting with Met\n\
  -T, --truncate                 Truncate full-length protein, from Met to Stop\n\
                                 Implies -F flag.\n\
");
#endif
    fprintf(stdout,"\
  -N, --nolengths                No intron lengths in alignment\n\
  -I, --invertmode=INT           Mode for alignments to (-) strand: 1=Invert\n\
                                 and print sense (-) strand; 2=Invert and\n\
                                 print antisense (+) strand.\n\
  -i, --introngap=INT            Nucleotides to show on each end of intron (default=3)\n\
  -l, --wraplength=INT           Wrap length for alignment (default=50)\n\
\n\
Help options\n\
  -V, --version                  Show version\n\
  -?, --help                     Show this help message\n\
");
  return;
}

/************************************************************************/



static void
print_npaths (int npaths, Sequence_T usersegment, Chrsubset_T chrsubset, 
#ifdef BETATEST
	      Chimera_T chimera, 
#endif
	      Failure_T failuretype) {
  if (npaths == 0) {
    printf("Paths (0):");
  } else if (usersegment != NULL) {
    /* npaths will actually be 2 here, one for each direction of usersegment */
    printf("Paths (1):");
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
  } else if (failuretype == REPETITIVE) {
    printf(" *** Repetitive sequence ***");
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
		      Sequence_T queryaaseq, Failure_T failuretype,
		      IIT_T chromosome_iit, Chrsubset_T chrsubset, IIT_T contig_iit, IIT_T altstrain_iit,
		      char *dbversion, int maxpaths, Chimera_T chimera, int maxmutations,
		      bool zerobasedp) {
  int pathnum;
#ifdef PMAP
  bool fulllengthp = false, truncatep = false;
#endif

  if (npaths == 0) {
#ifdef BETATEST
    print_npaths(0,usersegment,chrsubset,/*chimera*/NULL,failuretype);
#else
    print_npaths(0,usersegment,chrsubset,failuretype);
#endif
  } else {
#ifdef BETATEST
    print_npaths(npaths,usersegment,chrsubset,chimera,NO_FAILURE);
#else
    print_npaths(npaths,usersegment,chrsubset,NO_FAILURE);
#endif
    if (maxpaths == 0) {
      /* Special mode */
      if (npaths > 0) {
	Stage3_print_pathsummary(stage3array[0],/*pathnum*/1,chromosome_iit,contig_iit,
				 altstrain_iit,queryaaseq,fulllengthp,truncatep,
				 dbversion,maxmutations,zerobasedp,diagnosticp,maponlyp);
	if (chimera != NULL && npaths > 1) {
	  Stage3_print_pathsummary(stage3array[1],/*pathnum*/2,chromosome_iit,contig_iit,
				   altstrain_iit,queryaaseq,fulllengthp,truncatep,
				   dbversion,maxmutations,zerobasedp,diagnosticp,maponlyp);
	}
      }
    } else {
      for (pathnum = 1; pathnum <= npaths && pathnum <= maxpaths; pathnum++) {
	Stage3_print_pathsummary(stage3array[pathnum-1],pathnum,chromosome_iit,contig_iit,
				 altstrain_iit,queryaaseq,fulllengthp,truncatep,
				 dbversion,maxmutations,zerobasedp,diagnosticp,maponlyp);
      }
    }
  }

  return;
}


#ifdef HAVE_PTHREAD
static void *
input_thread (void *data) {
  Blackboard_T blackboard = (Blackboard_T) data;
  FILE *input = Blackboard_input(blackboard);
  Sequence_T queryseq, usersegment;
  Request_T request;
  int inputid = 1;		/* Initial queryseq already handled by main() */
  int nextchar = Blackboard_nextchar(blackboard);

  debug(printf("input_thread: Starting\n"));
  usersegment = Blackboard_usersegment(blackboard);
  while ((queryseq = Sequence_read(&nextchar,input,maponlyp)) != NULL) {
    debug(printf("input_thread: Putting request id %d\n",inputid));
    request = Request_new(inputid++,queryseq,usersegment);
    Blackboard_put_request(blackboard,request);
  }

  Blackboard_set_inputdone(blackboard);
  debug(printf("input_thread: Ending\n"));
  return (void *) NULL;
}
#endif


static void
print_result (Result_T result, Request_T request, Params_T params) {
  Sequence_T queryseq, usersegment;
  List_T matchlist;
  Stage1_T stage1;
  Stage3_T *stage3array, stage3;
  int npaths, pathnum;
  Chimera_T chimera = NULL;
  int chimerapos, chimeraequivpos, chimera_cdna_direction, worker_id;
  double donor_prob, acceptor_prob;
#ifdef PMAP
  bool fulllengthp = false, truncatep = false;
#endif

  queryseq = Request_queryseq(request);
  usersegment = Request_usersegment(request);
  worker_id = Result_worker_id(result);

  if (compressoutp == false && cdna_exons_p == false && genomic_exons_p == false &&
#ifdef PMAP
      psloutput_pro_p == false &&
#endif
      print_cdna_p == false && print_coordinates_p == false && 
      protein_genomic_p == false && psloutput_nt_p == false) {
    Sequence_print_header(queryseq,checksump);
  }

  if (compressoutp == true) {
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

    stage3array = Result_array(&npaths,result);
    if (npaths == 0) {
      fprintf(stderr,"No paths found for %s\n",Sequence_accession(queryseq));
    } else if (maxpaths == 0) {
      /* Special mode */
      if (npaths > 0) {
	Stage3_print_compressed(stage3array[0],queryseq,Params_chromosome_iit(params),
				dbversion,/*pathnum*/1,npaths,checksump,chimerapos,chimeraequivpos,
				donor_prob,acceptor_prob,chimera_cdna_direction,/*zerobasedp*/false,
				truncatep,worker_id);
	if (chimera != NULL && npaths > 1) {
	  Stage3_print_compressed(stage3array[1],queryseq,Params_chromosome_iit(params),
				  dbversion,/*pathnum*/2,npaths,checksump,chimerapos,chimeraequivpos,
				  donor_prob,acceptor_prob,chimera_cdna_direction,/*zerobasedp*/false,
				  truncatep,worker_id);
	}
      }
    } else {
      for (pathnum = 1; pathnum <= npaths && pathnum <= maxpaths; pathnum++) {
	Stage3_print_compressed(stage3array[pathnum-1],queryseq,Params_chromosome_iit(params),
				dbversion,pathnum,npaths,checksump,chimerapos,chimeraequivpos,
				donor_prob,acceptor_prob,chimera_cdna_direction,/*zerobasedp*/false,
				truncatep,worker_id);
      }
    }

  } else if (continuousp == true) {
    stage3array = Result_array(&npaths,result);
    if (npaths == 0) {
      printf("\n\n\n");
    } else {
      stage3 = stage3array[0];
      Stage3_print_alignment(stage3,queryseq,Params_chromosome_iit(params),summaryonlyp,/*universalp*/false,
			     /*zerobasedp*/false,continuousp,diagnosticp,/*flipgenomep*/true,
			     proteinmode,invertmode,nointronlenp,wraplength,ngap,maponlyp);
    }

  } else if (cdna_exons_p == true) {
    stage3array = Result_array(&npaths,result);
    Sequence_print_header(queryseq,checksump);
    Stage3_print_cdna_exons(stage3array[0],wraplength,ngap);

  } else if (genomic_exons_p == true) {
    stage3array = Result_array(&npaths,result);
    Sequence_print_header(queryseq,checksump);
    Stage3_print_genomic_exons(stage3array[0],wraplength,ngap);

  } else if (print_cdna_p == true) {
    stage3array = Result_array(&npaths,result);
    if (maxpaths == 0) {
      /* Special mode */
      if (npaths > 0) {
	Sequence_print_header(queryseq,checksump);
	Stage3_print_cdna(stage3array[0],queryseq,fulllengthp,truncatep,wraplength);
	if (Result_chimera(result) != NULL && npaths > 1) {
	  Sequence_print_header(queryseq,checksump);
	  Stage3_print_cdna(stage3array[1],queryseq,fulllengthp,truncatep,wraplength);
	}
      }
    } else {
      for (pathnum = 0; pathnum < npaths && pathnum < maxpaths; pathnum++) {
	Sequence_print_header(queryseq,checksump);
	Stage3_print_cdna(stage3array[pathnum],queryseq,fulllengthp,truncatep,wraplength);
      }
    }

  } else if (protein_genomic_p == true) {
    stage3array = Result_array(&npaths,result);
    if (maxpaths == 0) {
      /* Special mode */
      if (npaths > 0) {
	Sequence_print_header(queryseq,checksump);
	Stage3_print_protein_genomic(stage3array[0],queryseq,fulllengthp,truncatep,wraplength);
	if (Result_chimera(result) != NULL && npaths > 1) {
	  Sequence_print_header(queryseq,checksump);
	  Stage3_print_protein_genomic(stage3array[1],queryseq,fulllengthp,truncatep,wraplength);
	}
      }
    } else {
      for (pathnum = 0; pathnum < npaths && pathnum < maxpaths; pathnum++) {
	Sequence_print_header(queryseq,checksump);
	Stage3_print_protein_genomic(stage3array[pathnum],queryseq,fulllengthp,truncatep,wraplength);
      }
    }

  } else if (psloutput_nt_p == true) {
    stage3array = Result_array(&npaths,result);
    if (maxpaths == 0) {
      /* Special mode */
      if (npaths > 0) {
	Stage3_print_pslformat_nt(stage3array[0],/*pathnum*/1,Params_chromosome_iit(params),queryseq);
	if (Result_chimera(result) != NULL && npaths > 1) {
	  Stage3_print_pslformat_nt(stage3array[1],/*pathnum*/2,Params_chromosome_iit(params),queryseq);
	}
      }
    } else {
      for (pathnum = 0; pathnum < npaths && pathnum < maxpaths; pathnum++) {
	Stage3_print_pslformat_nt(stage3array[pathnum],/*pathnum*/pathnum+1,Params_chromosome_iit(params),queryseq);
      }
    }

#ifdef PMAP
  } else if (psloutput_pro_p == true) {
    stage3array = Result_array(&npaths,result);
    if (maxpaths == 0) {
      /* Special mode */
      if (npaths > 0) {
	Stage3_print_pslformat_pro(stage3array[0],/*pathnum*/1,Params_chromosome_iit(params),queryseq);
	if (Result_chimera(result) != NULL && npaths > 1) {
	  Stage3_print_pslformat_pro(stage3array[1],/*pathnum*/2,Params_chromosome_iit(params),queryseq);
	}
      }
    } else {
      for (pathnum = 0; pathnum < npaths && pathnum < maxpaths; pathnum++) {
	Stage3_print_pslformat_pro(stage3array[pathnum],/*pathnum*/pathnum+1,Params_chromosome_iit(params),queryseq);
      }
    }
#endif
  } else if (print_coordinates_p == true) {
    stage3array = Result_array(&npaths,result);
    Sequence_print_header(queryseq,checksump);
    Stage3_print_coordinates(stage3array[0],queryseq,Params_chromosome_iit(params),
			     /*zerobasedp*/false,invertmode,fulllengthp,truncatep,maponlyp);

  } else {
    /* Usual output */
    stage3array = Result_array(&npaths,result);
    if (npaths == 0) {
#ifdef BETATEST
      print_npaths(0,usersegment,Params_chrsubset(params),/*chimera*/NULL,Result_failuretype(result));
#else
      print_npaths(0,usersegment,Params_chrsubset(params),Result_failuretype(result));
#endif
    } else {
      chimera = Result_chimera(result);
      print_path_summaries(stage3array,npaths,usersegment,queryseq,Result_failuretype(result),
			   Params_chromosome_iit(params),Params_chrsubset(params),Params_contig_iit(params),
			   Params_altstrain_iit(params),dbversion,maxpaths,chimera,
			   Params_maxmutations(params),/*zerobasedp*/false);
    
      if (summaryonlyp == true || showalignp == true) {
	printf("Alignments:\n");
	if (maxpaths == 0) {
	  /* Special mode */
	  if (npaths > 0) {
	    printf("  Alignment for path 1:\n\n");
	    Stage3_print_alignment(stage3array[0],queryseq,Params_chromosome_iit(params),summaryonlyp,/*universalp*/false,
				   /*zerobasedp*/false,continuousp,diagnosticp,
				   /*flipgenomep*/true,proteinmode,invertmode,nointronlenp,wraplength,ngap,maponlyp);
	    if (chimera != NULL && npaths > 1) {
	      printf("  Alignment for path 2:\n\n");
	      Stage3_print_alignment(stage3array[1],queryseq,Params_chromosome_iit(params),summaryonlyp,/*universalp*/false,
				     /*zerobasedp*/false,continuousp,diagnosticp,
				     /*flipgenomep*/true,proteinmode,invertmode,nointronlenp,wraplength,ngap,maponlyp);
	    }
	  }
	} else {
	  for (pathnum = 1; pathnum <= npaths && pathnum <= maxpaths; pathnum++) {
	    printf("  Alignment for path %d:\n\n",pathnum);
	    Stage3_print_alignment(stage3array[pathnum-1],queryseq,Params_chromosome_iit(params),summaryonlyp,/*universalp*/false,
				   /*zerobasedp*/false,continuousp,diagnosticp,
				   /*flipgenomep*/true,proteinmode,invertmode,nointronlenp,wraplength,ngap,maponlyp);
	  }
	}
      }

      if (Params_map_iit(params) != NULL) {
	printf("Maps:\n");
	if (maxpaths == 0) {
	  /* Special mode */
	  if (npaths > 0) {
	    Stage3_print_map(stage3array[0],Params_map_iit(params),Params_chromosome_iit(params),
			     /*pathnum*/1,map_exons_p,map_bothstrands_p);
	    if (chimera != NULL && npaths > 1) {
	      Stage3_print_map(stage3array[1],Params_map_iit(params),Params_chromosome_iit(params),
			       /*pathnum*/2,map_exons_p,map_bothstrands_p);
	    }
	  }
	} else {
	  for (pathnum = 1; pathnum <= npaths && pathnum <= maxpaths; pathnum++) {
	    Stage3_print_map(stage3array[pathnum-1],Params_map_iit(params),Params_chromosome_iit(params),
			     pathnum,map_exons_p,map_bothstrands_p);
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
  Params_T params = Blackboard_params(blackboard);
  Result_T result;
  Request_T request;

  debug(printf("output_thread: Starting\n"));
  while ((result = Blackboard_get_result(&request,blackboard)) != NULL) {
    print_result(result,request,params);
    Result_free(&result);
    Request_free(&request);
  }

  debug(printf("output_thread: Ending\n"));
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
  Params_T params = Blackboard_params(blackboard);
  Result_T result;
  Request_T request;
  List_T request_queue = NULL, result_queue = NULL;
  int outputid = 0;

  debug(printf("output_thread: Starting\n"));
  while ((result = Blackboard_get_result(&request,blackboard)) != NULL) {
    if (Result_id(result) != outputid) {
      result_queue = queue_insert_result(result_queue,result);
      request_queue = queue_insert_request(request_queue,request);
    } else {
      print_result(result,request,params);
      Result_free(&result);
      Request_free(&request);
      outputid++;

      while (result_queue != NULL && Result_id(List_head(result_queue)) == outputid) {
	result_queue = List_pop(result_queue,(void **) &result);
	request_queue = List_pop(request_queue,(void **) &request);
	print_result(result,request,params);
	Result_free(&result);
	Request_free(&request);
	outputid++;
      }
    }
  }

  while (result_queue != NULL) {
    result_queue = List_pop(result_queue,(void **) &result);
    request_queue = List_pop(request_queue,(void **) &request);
    print_result(result,request,params);
    Result_free(&result);
    Request_free(&request);
  }

  debug(printf("output_thread: Ending\n"));
  return (void *) NULL;
}

#endif
    
static void
output_one (Result_T result, Request_T request, Params_T params) {
  print_result(result,request,params);
  Result_free(&result);
  return;
}


static Stage3_T *
stage3array_from_list (int *npaths, List_T stage3list) {
  Stage3_T *array1, *array0, x, y;
  bool *eliminate;
  int norig, i, j;

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
update_stage3list (int *maxconsecutive, List_T stage3list, Sequence_T queryseq, Sequence_T queryuc, Sequence_T genomicseg, 
		   Matchpair_T matchpair, Oligoindex_T oligoindex, int indexsize, Pairpool_T pairpool, 
		   int sufflookback, int nsufflookback, int badoligos, int straintype, char *strain,
		   Chrnum_T chrnum,  Genomicpos_T chroffset, Genomicpos_T chrpos, bool watsonp, int maxpeelback, 
		   int nullgap, int extramaterial_end, int extramaterial_paired, 
		   int extraband_single, int extraband_end, int extraband_paired, int ngap, 
		   Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR, IIT_T altstrain_iit) {
#ifdef PMAP
  Sequence_T queryntseq;
#endif
  Sequence_T genomicuc;
  Stage2_T stage2;
  Stage3_T stage3;
  Matchpairend_T matchpairend;

  if (matchpair == NULL) {
    matchpairend = MIXED;
  } else {
    matchpairend = Matchpair_matchpairend(matchpair);
  }

#ifdef PMAP
  queryntseq = Sequence_convert_to_nucleotides(queryseq);
#endif

  if (user_genomicseg == NULL && uncompressedp == false && straintype == 0) {
    genomicuc = Sequence_alias(genomicseg);
  } else {
    genomicuc = Sequence_uppercase(genomicseg);
  }

  debug(printf("Beginning Stage2_compute\n"));
  stage2 = Stage2_compute(&(*maxconsecutive),queryseq,queryuc,
#ifdef PMAP
			  queryntseq,
#endif
			  genomicseg,genomicuc,oligoindex,indexsize,
			  pairpool,sufflookback,nsufflookback,maxintronlen_bound,
			  badoligos,crossspeciesp);

  if (stage2 != NULL) {
    debug(printf("Beginning Stage3_compute\n"));
    if ((stage3 = Stage3_compute(stage2,
#ifdef PMAP
				 queryseq,queryntseq,queryntseq,
#else
				 queryseq,queryuc,
#endif
				 genomicseg,genomicuc,matchpairend,straintype,strain,chrnum,chroffset,
				 chrpos,watsonp,maxpeelback,nullgap,extramaterial_end,extramaterial_paired,
				 extraband_single,extraband_end,extraband_paired,
				 extend_mismatch_p,end_microexons_p,pairpool,dynprogL,dynprogM,dynprogR,
				 altstrain_iit,indexsize,ngap)) == NULL) {
    } else {
      stage3list = List_push(stage3list,stage3);
    }
    Stage2_free(&stage2);
  }

  Sequence_free(&genomicuc);
#ifdef PMAP
  Sequence_free(&queryntseq);
#endif

  return stage3list;
}

static List_T
update_stage3list_maponlyp (int *maxconsecutive, List_T stage3list, Sequence_T queryseq, Sequence_T queryuc,
			    char *gbuffer1, char *gbuffer2, int gbufferlen,
			    Matchpair_T matchpair, Oligoindex_T oligoindex, int indexsize, Pairpool_T pairpool, 
			    Genome_T genome, int sufflookback, int nsufflookback, int straintype, char *strain,
			    Chrnum_T chrnum,  Genomicpos_T chroffset, Genomicpos_T chrpos, bool watsonp, int maxpeelback, 
			    int nullgap, int extramaterial_end, int extramaterial_paired, 
			    int extraband_single, int extraband_end, int extraband_paired, int ngap, 
			    Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR, IIT_T altstrain_iit) {
  Stage3_T stage3;
  Matchpairend_T matchpairend;
#ifdef PMAP
  Sequence_T queryntseq;
#endif

#ifdef PMAP
  queryntseq = Sequence_convert_to_nucleotides(queryseq);
#endif

  if (matchpair == NULL) {
    matchpairend = MIXED;
  } else {
    matchpairend = Matchpair_matchpairend(matchpair);
  }

  stage3 = Stage3_direct(matchpair,
#ifdef PMAP
			 queryseq,queryntseq,queryntseq,
#else
			 queryseq,queryuc,
#endif
			 pairpool,genome,matchpairend,
			 chrnum,chroffset,chrpos,watsonp,ngap,gbuffer1,gbuffer2,gbufferlen,dynprogL,dynprogR,
			 extramaterial_end,extraband_end);
  stage3list = List_push(stage3list,stage3);

  return stage3list;
}


/* This code is duplicated in get-genome.c */
static IIT_T global_altstrain_iit = NULL;

static int
index_compare (const void *a, const void *b) {
  int index1 = * (int *) a;
  int index2 = * (int *) b;
  int type1, type2;
  Genomicpos_T pos1, pos2;

  type1 = Interval_type(IIT_interval(global_altstrain_iit,index1));
  type2 = Interval_type(IIT_interval(global_altstrain_iit,index2));
  
  if (type1 < type2) {
    return -1;
  } else if (type1 > type2) {
    return +1;
  } else {
    /* Store in descending genomic position, so right shifting works
       in Genome_patch_strain */
    pos1 = Interval_low(IIT_interval(global_altstrain_iit,index1));
    pos2 = Interval_low(IIT_interval(global_altstrain_iit,index2));

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
stage3_from_usersegment (int *npaths, Sequence_T queryseq, Sequence_T queryuc,
			 Sequence_T usersegment, Params_T params, 
			 int maxpeelback, int sufflookback, int nsufflookback, 
			 int badoligos, int nullgap, 
			 Oligoindex_T oligoindex, Pairpool_T pairpool, 
			 Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR) {
  List_T stage3list;
  Genomicpos_T chroffset, chrpos;
  Sequence_T revcomp;
  int indexsize = Params_indexsize(params);
  int extramaterial_end = Params_extramaterial_end(params);
  int extramaterial_paired = Params_extramaterial_paired(params);
  int extraband_single = Params_extraband_single(params);
  int extraband_end = Params_extraband_end(params);
  int extraband_paired = Params_extraband_paired(params);
  int maxconsecutive;
  Chrnum_T chrnum = 0;
  IIT_T altstrain_iit = NULL;
		    
  chroffset = chrpos = 0U;

  maxconsecutive = 0;
  stage3list = update_stage3list(&maxconsecutive,/*stage3list*/NULL,queryseq,queryuc,
				 usersegment,(Matchpair_T) NULL,
				 oligoindex,indexsize,pairpool,sufflookback,nsufflookback,badoligos,
				 /*straintype*/0,/*strain*/NULL,chrnum,chroffset,chrpos,
				 /*watsonp*/true,maxpeelback,nullgap,extramaterial_end,extramaterial_paired,
				 extraband_single,extraband_end,extraband_paired,
				 ngap,dynprogL,dynprogM,dynprogR,altstrain_iit);

  revcomp = Sequence_revcomp(usersegment);

  maxconsecutive = 0;
  stage3list = update_stage3list(&maxconsecutive,stage3list,queryseq,queryuc,
				 revcomp,(Matchpair_T) NULL,
				 oligoindex,indexsize,pairpool,sufflookback,nsufflookback,badoligos,
				 /*straintype*/0,/*strain*/NULL,chrnum,chroffset,chrpos,
				 /*watsonp*/false,maxpeelback,nullgap,extramaterial_end,extramaterial_paired,
				 extraband_single,extraband_end,extraband_paired,
				 ngap,dynprogL,dynprogM,dynprogR,altstrain_iit);

  if (stage3list == NULL) {
    *npaths = 0;
    return NULL;
  } else {
    return stage3array_from_list(&(*npaths),stage3list);
  }
}


static Stage3_T *
stage3_from_matchlist (int *npaths, Stage3_T *oldstage3array, Stage1_T stage1,
		       List_T matchlist, Sequence_T queryseq, Sequence_T queryuc,
		       Sequence_T usersegment, Params_T params, 
		       int maxpeelback, int sufflookback, int nsufflookback, 
		       int badoligos, int nullgap, 
		       Oligoindex_T oligoindex, Pairpool_T pairpool, 
		       Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
		       char *gbuffer1, char *gbuffer2, char *gbuffer3, int gbufferlen) {
  List_T stage3list = NULL, p;
  Matchpair_T matchpair;
  Match_T start;
  Chrnum_T chrnum;
  char *strain;
  Genomicpos_T chroffset, chrlength, chrpos, genomicstart, genomiclength;
  bool watsonp;
  Sequence_T genomicseg;
  int indexsize = Params_indexsize(params);
  int extramaterial_end = Params_extramaterial_end(params);
  int extramaterial_paired = Params_extramaterial_paired(params);
  int extraband_single = Params_extraband_single(params);
  int extraband_end = Params_extraband_end(params);
  int extraband_paired = Params_extraband_paired(params);
  int maxconsecutive = 0;
  IIT_T altstrain_iit = Params_altstrain_iit(params);
  IIT_T chromosome_iit = Params_chromosome_iit(params);
  int *indexarray, nindices, choplen, straintype, i, j;
		    
  for (i = 0; i < *npaths; i++) {
    stage3list = List_push(stage3list,(void *) oldstage3array[i]);
  }
  if (oldstage3array != NULL) {
    FREE(oldstage3array);
  }

  for (p = matchlist; p != NULL; p = List_next(p)) {
    matchpair = List_head(p);
    if (Matchpair_usep(matchpair) == true) {
      start = Matchpair_bound5(matchpair);
      chrnum = Match_chrnum(start);
      if (usersegment != NULL) {
	chrlength = Sequence_fulllength(usersegment);
	chroffset = 0U;
	Matchpair_get_coords(&chrpos,&genomicstart,&genomiclength,&watsonp,matchpair,stage1,
			     chrlength,Params_maxextension(params),maponlyp);
	strain = NULL;
	genomicseg = Sequence_substring(usersegment,genomicstart,genomiclength,!watsonp,
					gbuffer1,gbuffer2,gbufferlen);
	stage3list = update_stage3list(&maxconsecutive,stage3list,queryseq,queryuc,genomicseg,matchpair,
				       oligoindex,indexsize,pairpool,sufflookback,nsufflookback,badoligos,
				       /*straintype*/0,/*strain*/NULL,chrnum,chroffset,chrpos,
				       watsonp,maxpeelback,nullgap,extramaterial_end,extramaterial_paired,
				       extraband_single,extraband_end,extraband_paired,
				       ngap,dynprogL,dynprogM,dynprogR,altstrain_iit);
	Sequence_free(&genomicseg);

      } else {
	chrlength = Chrnum_length(chrnum,chromosome_iit);
	chroffset = Chrnum_offset(chrnum,chromosome_iit);
	Matchpair_get_coords(&chrpos,&genomicstart,&genomiclength,&watsonp,matchpair,stage1,
			     chrlength,Params_maxextension(params),maponlyp);

	if (altstrain_iit == NULL) {
	  strain = NULL;
	} else {
	  strain = Params_refstrain(params);
	}

	if (maponlyp == true) {
	  stage3list = update_stage3list_maponlyp(&maxconsecutive,stage3list,queryseq,queryuc,
						  gbuffer1,gbuffer2,gbufferlen,matchpair,oligoindex,indexsize,pairpool,
						  Params_genome(params),sufflookback,nsufflookback,
						  /*straintype*/0,/*strain*/NULL,chrnum,chroffset,chrpos,
						  watsonp,maxpeelback,nullgap,extramaterial_end,extramaterial_paired,
						  extraband_single,extraband_end,extraband_paired,
						  ngap,dynprogL,dynprogM,dynprogR,altstrain_iit);
	} else {
	  genomicseg = Genome_get_segment(Params_genome(params),genomicstart,genomiclength,!watsonp,
					  gbuffer1,gbuffer2,gbufferlen);
	  stage3list = update_stage3list(&maxconsecutive,stage3list,queryseq,queryuc,genomicseg,matchpair,
					 oligoindex,indexsize,pairpool,sufflookback,nsufflookback,badoligos,
					 /*straintype*/0,/*strain*/NULL,chrnum,chroffset,chrpos,
					 watsonp,maxpeelback,nullgap,extramaterial_end,extramaterial_paired,
					 extraband_single,extraband_end,extraband_paired,
					 ngap,dynprogL,dynprogM,dynprogR,altstrain_iit);
	  Sequence_free(&genomicseg);
	}
      }

      /* We rely upon the fact that gbuffer1 still holds the genomic segment.  This code is duplicated in get-genome.c */
      if (maponlyp == true) {
	/* Don't do strain calculations */
      } else if (altstrain_iit != NULL) {
	indexarray = IIT_get(&nindices,altstrain_iit,genomicstart+1,genomicstart+genomiclength-1);
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
					     genomicstart,genomiclength,!watsonp,
					     gbuffer1,gbuffer2,gbuffer3,gbufferlen);
	    stage3list = update_stage3list(&maxconsecutive,stage3list,queryseq,queryuc,genomicseg,matchpair,
					   oligoindex,indexsize,pairpool,sufflookback,nsufflookback,badoligos,
					   straintype,strain,chrnum,chroffset,chrpos,
					   watsonp,maxpeelback,nullgap,extramaterial_end,extramaterial_paired,
					   extraband_single,extraband_end,extraband_paired,
					   ngap,dynprogL,dynprogM,dynprogR,altstrain_iit);
	    Sequence_free(&genomicseg);
	  }
	  FREE(indexarray);
	}
      }
    }
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
			   pairpool,Params_genome(params),ngap) == true) {
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
		   Params_T params, int maxpeelback, int sufflookback, int nsufflookback, int badoligos, int nullgap,
		   Oligoindex_T oligoindex, Matchpool_T matchpool, Pairpool_T pairpool,
		   Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
		   char *gbuffer1, char *gbuffer2, char *gbuffer3, int gbufferlen) {
  Stage1_T stage1sub;
  List_T matchlist = NULL;
  Stage3_T *stage3array_sub1 = NULL, *stage3array_sub2 = NULL, stage3;
  Sequence_T querysubseq = NULL, querysubuc = NULL;
  int chimerapos, chimeraequivpos;
  int bestfrom, bestto;
  int five_margin, three_margin, five_score = 0, three_score = 0;
  int npaths_sub1 = 0, npaths_sub2 = 0;
  int i, queryntlength;

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

	  stage1sub = Stage1_compute(querysubuc,
#ifdef PMAP
				     Params_indexdb_fwd(params),Params_indexdb_rev(params),
#else
				     Params_indexdb(params),
#endif
				     Params_chromosome_iit(params),Params_chrsubset(params),
				     matchpool,maxtotallen_bound,Params_stuttercycles(params),
				     Params_stutterhits(params),crossspeciesp,/*subsequencep*/true,
				     Params_genome(params));
	  debug2(printf("Getting Stage 1 matchlist\n"));
	  matchlist = Stage1_matchlist(stage1sub);
	  debug2(printf("Performing Stage 3\n"));
	  stage3array = stage3_from_matchlist(&(*npaths),stage3array,stage1sub,matchlist,querysubseq,
					      querysubuc,usersegment,params,maxpeelback,sufflookback,nsufflookback,
					      badoligos,nullgap,
					      oligoindex,pairpool,dynprogL,dynprogM,dynprogR,gbuffer1,gbuffer2,
					      gbuffer3,gbufferlen);
	  Stage1_free(&stage1sub,/*completep*/true);
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

	  stage1sub = Stage1_compute(querysubuc,
#ifdef PMAP
				     Params_indexdb_fwd(params),Params_indexdb_rev(params),
#else
				     Params_indexdb(params),
#endif
				     Params_chromosome_iit(params),Params_chrsubset(params),
				     matchpool,maxtotallen_bound,Params_stuttercycles(params),
				     Params_stutterhits(params),crossspeciesp,/*subsequencep*/true,
				     Params_genome(params));
	  debug2(printf("Getting Stage 1 matchlist\n"));
	  matchlist = Stage1_matchlist(stage1sub);
	  debug2(printf("Performing Stage 3\n"));
	  stage3array = stage3_from_matchlist(&(*npaths),stage3array,stage1sub,matchlist,querysubseq,
					      querysubuc,usersegment,params,maxpeelback,sufflookback,nsufflookback,
					      badoligos,nullgap,
					      oligoindex,pairpool,dynprogL,dynprogM,dynprogR,gbuffer1,gbuffer2,
					      gbuffer3,gbufferlen);
	  Stage1_free(&stage1sub,/*completep*/true);
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
			   Params_genome(params),Sequence_fulllength(queryseq)) == true) {
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
apply_stage3 (Chimera_T *chimera, int *npaths, Stage1_T stage1, Sequence_T queryseq, Sequence_T queryuc,
	      Sequence_T usersegment, Params_T params, int maxpeelback, int sufflookback, int nsufflookback,
	      int badoligos, int nullgap, Oligoindex_T oligoindex, Matchpool_T matchpool, Pairpool_T pairpool, 
	      Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
	      char *gbuffer1, char *gbuffer2, char *gbuffer3, int gbufferlen) {
  Stage3_T *stage3array = NULL, nonchimericbest;
  List_T matchlist = NULL;
  bool testchimerap = false;
#ifdef BETATEST
  int effective_start, effective_end;
  Genome_T genome = Params_genome(params);
#endif
  
  *npaths = 0;
  *chimera = NULL;

  matchlist = Stage1_matchlist(stage1);
  stage3array = stage3_from_matchlist(&(*npaths),stage3array,stage1,matchlist,queryseq,queryuc,
				      usersegment,params,maxpeelback,sufflookback,nsufflookback,
				      badoligos,nullgap,
				      oligoindex,pairpool,dynprogL,dynprogM,dynprogR,gbuffer1,gbuffer2,
				      gbuffer3,gbufferlen);

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
				      queryseq,queryuc,usersegment,params,
				      maxpeelback,sufflookback,nsufflookback,badoligos,nullgap,
				      oligoindex,matchpool,pairpool,dynprogL,dynprogM,dynprogR,
				      gbuffer1,gbuffer2,gbuffer3,gbufferlen);
#endif
    }
  }

  return stage3array;
}


static void
handle_request (Request_T request, int worker_id, Matchpool_T matchpool, Pairpool_T pairpool, 
		Oligoindex_T oligoindex, Params_T params,
		char *gbuffer1, char *gbuffer2, char *gbuffer3, int gbufferlen, 
		Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR, Reqpost_T reqpost) {
  int jobid;
  Sequence_T queryseq, queryuc, usersegment;
  int badoligos = 0, trimoligos, trim_start, trim_end;
  double oligodepth;
  Chimera_T chimera = NULL;
  List_T matchlist;
  bool badsequencep = false;

  Stage1_T stage1;
  Stage3_T *stage3array;
  int npaths;

  int stuttercycles = Params_stuttercycles(params);
  int stutterhits = Params_stutterhits(params);
  int maxpeelback = Params_maxpeelback(params);
  int sufflookback = Params_sufflookback(params);
  int nsufflookback = Params_nsufflookback(params);
  int nullgap = Params_nullgap(params);

  jobid = Request_id(request);
  usersegment = Request_usersegment(request);
  debug(printf("worker_thread %d: Got request id\n",worker_id,jobid));
  queryseq = Request_queryseq(request);
  Matchpool_reset(matchpool);
  Pairpool_reset(pairpool);

  if (Sequence_fulllength(queryseq) <= 0) {
    if (reqpost == NULL) {
      output_one(Result_new(jobid,worker_id,(Chimera_T) NULL,(Stage3_T *) NULL,0,EMPTY_SEQUENCE),request,params);
    } else {
#ifdef HAVE_PTHREAD
      Reqpost_put_result(reqpost,Result_new(jobid,worker_id,(Chimera_T) NULL,(Stage3_T *) NULL,0,EMPTY_SEQUENCE));
#endif
    }
      
  } else {			/* Sequence_fulllength(queryseq) > 0 */
    queryuc = Sequence_uppercase(queryseq);
    if (maponlyp == false) {
#ifdef PMAP
      Oligoindex_set_inquery(&badoligos,&trimoligos,&trim_start,&trim_end,oligoindex,queryuc,/*trimp*/false);
#else
      oligodepth = Oligoindex_set_inquery(&badoligos,&trimoligos,&trim_start,&trim_end,oligoindex,queryuc,/*trimp*/true);
      debug1(printf("oligo depth: %f\n",oligodepth));
      debug1(printf("badoligos: %d/%d\n",badoligos,trimoligos));
      debug1(printf("trimstart: %d\n",trim_start));
      debug1(printf("trimend: %d\n",trim_end));
      if (crossspeciesp == true) {
	badsequencep = false;
      } else if ((double) badoligos/(double) trimoligos > MAX_BADOLIGOS) {
	badsequencep = true;
      }
#endif
    }

    if (badsequencep == true) {
      if (reqpost == NULL) {
	output_one(Result_new(jobid,worker_id,(Chimera_T) NULL,(Stage3_T *) NULL,0,REPETITIVE),request,params);
      } else {
#ifdef HAVE_PTHREAD
	Reqpost_put_result(reqpost,Result_new(jobid,worker_id,(Chimera_T) NULL,(Stage3_T *) NULL,0,REPETITIVE));
#endif
      }
    } else if (usersegment != NULL && userstage1p == false) {
      stage3array = stage3_from_usersegment(&npaths,queryseq,queryuc,usersegment,params,
					    maxpeelback,sufflookback,nsufflookback,
					    badoligos,nullgap,
					    oligoindex,pairpool,dynprogL,dynprogM,dynprogR);
      if (reqpost == NULL) {
	output_one(Result_new(jobid,worker_id,(Chimera_T) NULL,stage3array,npaths,NO_FAILURE),request,params);
      } else {
#ifdef HAVE_PTHREAD
	Reqpost_put_result(reqpost,Result_new(jobid,worker_id,(Chimera_T) NULL,stage3array,npaths,NO_FAILURE));
#endif
      }

    } else {		/* Not user segment and not maponly */
#ifndef PMAP
      if (maponlyp == false) {
	Sequence_trim(queryseq,trim_start,trim_end);
	Sequence_trim(queryuc,trim_start,trim_end);
      }
#endif
      stage1 = Stage1_compute(queryuc,
#ifdef PMAP
			      Params_indexdb_fwd(params),Params_indexdb_rev(params),
#else
			      Params_indexdb(params),
#endif
			      Params_chromosome_iit(params),Params_chrsubset(params),matchpool,
			      maxtotallen_bound,stuttercycles,stutterhits,crossspeciesp,/*subsequencep*/false,
			      Params_genome(params));
      stage3array = apply_stage3(&chimera,&npaths,stage1,queryseq,queryuc,usersegment,params,
				 maxpeelback,sufflookback,nsufflookback,badoligos,nullgap,
				 oligoindex,matchpool,pairpool,dynprogL,dynprogM,dynprogR,
				 gbuffer1,gbuffer2,gbuffer3,gbufferlen);
      Stage1_free(&stage1,/*completep*/true);
      if (reqpost == NULL) {
	output_one(Result_new(jobid,worker_id,chimera,stage3array,npaths,NO_FAILURE),request,params);
      } else {
#ifdef HAVE_PTHREAD
	Reqpost_put_result(reqpost,Result_new(jobid,worker_id,chimera,stage3array,npaths,NO_FAILURE));
#endif
      }
    } /* Matches not user segment and not maponly */

    Sequence_free(&queryuc);
  } /* Matches sequence length > 0 */

  return;
}


static void
single_thread (FILE *input, int nextchar, Sequence_T queryseq, Params_T params, Sequence_T usersegment) {
  Oligoindex_T oligoindex;
  Dynprog_T dynprogL, dynprogM, dynprogR;
  char *gbuffer1 = NULL, *gbuffer2 = NULL, *gbuffer3 = NULL; /* Individual space for reading genome sequence,
								revcomp, and strain calculations */
  int gbufferlen = 0;
  Matchpool_T matchpool;
  Pairpool_T pairpool;
  Request_T request;
  int jobid = 0;

  int nullgap = Params_nullgap(params);
  int maxpeelback = Params_maxpeelback(params);
  int extramaterial_end = Params_extramaterial_end(params);
  int extramaterial_paired = Params_extramaterial_paired(params);

  oligoindex = Oligoindex_new();
  dynprogL = Dynprog_new(nullgap,EXTRAQUERYGAP,maxpeelback,extramaterial_end,extramaterial_paired);
  dynprogM = Dynprog_new(nullgap,EXTRAQUERYGAP,maxpeelback,extramaterial_end,extramaterial_paired);
  dynprogR = Dynprog_new(nullgap,EXTRAQUERYGAP,maxpeelback,extramaterial_end,extramaterial_paired);
  gbufferlen = MAXSEQLEN + maxtotallen_bound + 2*Params_maxextension(params) + 
    1000;			/* Add extra 1000 just to be sure */
  gbuffer1 = (char *) CALLOC(gbufferlen+1,sizeof(char));
  gbuffer2 = (char *) CALLOC(gbufferlen+1,sizeof(char));
  gbuffer3 = (char *) CALLOC(gbufferlen+1,sizeof(char));
  matchpool = Matchpool_new();
  pairpool = Pairpool_new();

  while (jobid == 0 || (queryseq = Sequence_read(&nextchar,input,maponlyp)) != NULL) {
    request = Request_new(jobid++,queryseq,usersegment);
TRY
  handle_request(request,/*worker_id*/0,matchpool,pairpool,oligoindex,params,gbuffer1,gbuffer2,gbuffer3,gbufferlen,
		   dynprogL,dynprogM,dynprogR,(Reqpost_T) NULL);
ELSE
    if (Sequence_accession(queryseq) == NULL) {
      fprintf(stderr,"Problem with unnamed sequence (%d bp)\n",Sequence_fulllength(queryseq));
    } else {
      fprintf(stderr,"Problem with sequence %s (%d bp)\n",
	      Sequence_accession(queryseq),Sequence_fulllength(queryseq));
    }
  RERAISE;
END_TRY;
    Request_free(&request);
   /* Don't free queryseq; done by Request_free */
  }

  Pairpool_free(&pairpool);
  Matchpool_free(&matchpool);
  FREE(gbuffer3);
  FREE(gbuffer2);
  FREE(gbuffer1);
  Dynprog_free(&dynprogR);
  Dynprog_free(&dynprogM);
  Dynprog_free(&dynprogL);
  Oligoindex_free(&oligoindex);

  return;
}


#ifdef HAVE_PTHREAD
static void *
worker_thread (void *data) {
  Oligoindex_T oligoindex;
  Dynprog_T dynprogL, dynprogM, dynprogR;
  char *gbuffer1 = NULL, *gbuffer2 = NULL, *gbuffer3 = NULL; /* Individual space for reading genome sequence,
								revcomp, and strain calculations */
  int gbufferlen = 0;
  Matchpool_T matchpool;
  Pairpool_T pairpool;
  Request_T request;
  Sequence_T queryseq;

  Reqpost_T reqpost = (Reqpost_T) data;
  Params_T params = Reqpost_params(reqpost);
  int worker_id = Reqpost_id(reqpost);

  int nullgap = Params_nullgap(params);
  int maxpeelback = Params_maxpeelback(params);
  int extramaterial_end = Params_extramaterial_end(params);
  int extramaterial_paired = Params_extramaterial_paired(params);

  /* Thread-specific data and storage */
  oligoindex = Oligoindex_new();
  dynprogL = Dynprog_new(nullgap,EXTRAQUERYGAP,maxpeelback,extramaterial_end,extramaterial_paired);
  dynprogM = Dynprog_new(nullgap,EXTRAQUERYGAP,maxpeelback,extramaterial_end,extramaterial_paired);
  dynprogR = Dynprog_new(nullgap,EXTRAQUERYGAP,maxpeelback,extramaterial_end,extramaterial_paired);
  gbufferlen = MAXSEQLEN + maxtotallen_bound + 2*Params_maxextension(params) + 
    1000;			/* Add extra 1000 just to be sure */
  gbuffer1 = (char *) CALLOC(gbufferlen+1,sizeof(char));
  gbuffer2 = (char *) CALLOC(gbufferlen+1,sizeof(char));
  gbuffer3 = (char *) CALLOC(gbufferlen+1,sizeof(char));
  matchpool = Matchpool_new();
  pairpool = Pairpool_new();

  while ((request = Reqpost_get_request(reqpost))) {
TRY
  handle_request(request,worker_id,matchpool,pairpool,oligoindex,params,gbuffer1,gbuffer2,gbuffer3,gbufferlen,
	  	   dynprogL,dynprogM,dynprogR,reqpost);
ELSE
    queryseq = Request_queryseq(request);
    if (Sequence_accession(queryseq) == NULL) {
      fprintf(stderr,"Problem with unnamed sequence (%d bp)\n",Sequence_fulllength(queryseq));
    } else {
      fprintf(stderr,"Problem with sequence %s (%d bp)\n",
	      Sequence_accession(queryseq),Sequence_fulllength(queryseq));
    }
  RERAISE;
END_TRY;
   /* Don't free request; done by output thread */
   /* Don't free queryseq; done by Request_free */
  }

  Pairpool_free(&pairpool);
  Matchpool_free(&matchpool);
  FREE(gbuffer3);
  FREE(gbuffer2);
  FREE(gbuffer1);
  Dynprog_free(&dynprogR);
  Dynprog_free(&dynprogM);
  Dynprog_free(&dynprogL);
  Oligoindex_free(&oligoindex);

  return (void *) NULL;
}
#endif


static void
align_mutations (FILE *input, int nextchar, Sequence_T queryseq, Params_T params, Sequence_T mutationrefseq) {
  Oligoindex_T oligoindex;
  int badoligos, trimoligos;
  Dynprog_T dynprogL, dynprogM, dynprogR;
  char *gbuffer1 = NULL, *gbuffer2 = NULL, *gbuffer3 = NULL; /* Individual space for reading genome sequence,
								revcomp, and strain calculations */
  int trim_start, trim_end;
  int gbufferlen = 0,
    stuttercycles, stutterhits, 
    maxpeelback, sufflookback, nsufflookback, nullgap;
  int extramaterial_end, extramaterial_paired;
  Matchpool_T matchpool;
  Pairpool_T pairpool;
  Genomicpos_T genomicstart, genomiclength;
  Sequence_T genomicseg, queryuc, mutationrefuc;
  int jobid = 0;

  Chimera_T chimera = NULL;
  Stage1_T stage1;
  Stage3_T *stage3array, stage3, stage3ref;
  int npaths, i;

  stuttercycles = Params_stuttercycles(params);
  stutterhits = Params_stutterhits(params);
  maxpeelback = Params_maxpeelback(params);
  sufflookback = Params_sufflookback(params);
  nsufflookback = Params_nsufflookback(params);
  nullgap = Params_nullgap(params);
  extramaterial_end = Params_extramaterial_end(params);
  extramaterial_paired = Params_extramaterial_paired(params);

  oligoindex = Oligoindex_new();
  dynprogL = Dynprog_new(nullgap,EXTRAQUERYGAP,maxpeelback,extramaterial_end,extramaterial_paired);
  dynprogM = Dynprog_new(nullgap,EXTRAQUERYGAP,maxpeelback,extramaterial_end,extramaterial_paired);
  dynprogR = Dynprog_new(nullgap,EXTRAQUERYGAP,maxpeelback,extramaterial_end,extramaterial_paired);
  gbufferlen = MAXSEQLEN + maxtotallen_bound + 2*Params_maxextension(params) + 
    1000;			/* Add extra 1000 just to be sure */
  gbuffer1 = (char *) CALLOC(gbufferlen+1,sizeof(char));
  gbuffer2 = (char *) CALLOC(gbufferlen+1,sizeof(char));
  gbuffer3 = (char *) CALLOC(gbufferlen+1,sizeof(char));
  matchpool = Matchpool_new();
  pairpool = Pairpool_new();

  Matchpool_reset(matchpool);
  Pairpool_reset(pairpool);

  mutationrefuc = Sequence_uppercase(mutationrefseq);

  /* Do not trim the mutation refseq */
  Oligoindex_set_inquery(&badoligos,&trimoligos,&trim_start,&trim_end,oligoindex,mutationrefuc,/*trimp*/false);
#ifndef PMAP
  Sequence_trim(mutationrefseq,trim_start,trim_end);
#endif
  stage1 = Stage1_compute(mutationrefuc,
#ifdef PMAP
			  Params_indexdb_fwd(params),Params_indexdb_rev(params),
#else
			  Params_indexdb(params),
#endif
			  Params_chromosome_iit(params),Params_chrsubset(params),matchpool,
			  maxtotallen_bound,stuttercycles,stutterhits,crossspeciesp,/*subsequencep*/false,
			  Params_genome(params));
  stage3array = apply_stage3(&chimera,&npaths,stage1,mutationrefseq,mutationrefuc,/*usersegment*/NULL,params,
			     maxpeelback,sufflookback,nsufflookback,badoligos,nullgap,
			     oligoindex,matchpool,pairpool,dynprogL,dynprogM,dynprogR,
			     gbuffer1,gbuffer2,gbuffer3,gbufferlen);
  /* chimera should be NULL */
  for (i = 1; i < npaths; i++) {
    stage3 = stage3array[i];
    Stage3_free(&stage3);
  }
  if (npaths > 0) {
    stage3ref = stage3array[0];
#ifdef PMAP
    stage3ref = Stage3_translate_cdna(stage3ref,queryseq);
    stage3ref = Stage3_backtranslate_cdna(stage3ref,/*diagnosticp*/false);
#else
    stage3ref = Stage3_translate_genomic(stage3ref,/*fulllengthp*/true,/*truncatep*/false);
#endif
    FREE(stage3array);

    Stage3_genomicbounds(&genomicstart,&genomiclength,stage3ref);
    genomicseg = Genome_get_segment(Params_genome(params),genomicstart,genomiclength,/*revcomp*/false,
				    gbuffer1,gbuffer2,gbufferlen);

    while (jobid == 0 || (queryseq = Sequence_read(&nextchar,input,maponlyp)) != NULL) {
      Matchpool_reset(matchpool);
      Pairpool_reset(pairpool);

      Sequence_print_header(queryseq,checksump);
      if (Sequence_fulllength(queryseq) <= 0) {
#ifdef BETATEST
	print_npaths(0,/*usersegment*/NULL,Params_chrsubset(params),/*chimera*/NULL,EMPTY_SEQUENCE);
#else
	print_npaths(0,/*usersegment*/NULL,Params_chrsubset(params),EMPTY_SEQUENCE);
#endif

      } else {
	queryuc = Sequence_uppercase(queryseq);

#ifdef PMAP
	Oligoindex_set_inquery(&badoligos,&trimoligos,&trim_start,&trim_end,oligoindex,queryuc,/*trimp*/false);
#else
	Oligoindex_set_inquery(&badoligos,&trimoligos,&trim_start,&trim_end,oligoindex,queryuc,/*trimp*/true);
	if ((double) badoligos/(double) trimoligos > MAX_BADOLIGOS) {
#ifdef BETATEST
	  print_npaths(0,/*usersegment*/NULL,Params_chrsubset(params),/*chimera*/NULL,REPETITIVE);
#else
	  print_npaths(0,/*usersegment*/NULL,Params_chrsubset(params),REPETITIVE);
#endif
	} else {
#endif
	  stage3array = stage3_from_usersegment(&npaths,queryseq,queryuc,genomicseg,params,
						maxpeelback,sufflookback,nsufflookback,
						badoligos,nullgap,
						oligoindex,pairpool,dynprogL,dynprogM,dynprogR);
	  if (npaths == 0) {
#ifdef BETATEST
	    print_npaths(0,/*usersegment*/NULL,Params_chrsubset(params),/*chimera*/NULL,NO_FAILURE);
#else
	    print_npaths(0,/*usersegment*/NULL,Params_chrsubset(params),NO_FAILURE);
#endif
	  } else {
#ifdef BETATEST
	    print_npaths(1,/*usersegment*/NULL,Params_chrsubset(params),/*chimera*/NULL,NO_FAILURE);
#else
	    print_npaths(1,/*usersegment*/NULL,Params_chrsubset(params),NO_FAILURE);
#endif
#ifndef PMAP
	    Stage3_translate_cdna_via_reference(stage3array[0],stage3ref,literalrefp);
#endif
	    Stage3_fix_cdna_direction(stage3array[0],stage3ref);
	    Stage3_print_mutations(stage3array[0],stage3ref,Params_chromosome_iit(params),
				   dbversion,showalignp,/*zerobasedp*/false,continuousp,diagnosticp,proteinmode,
				   invertmode,nointronlenp,wraplength,ngap,Params_maxmutations(params));
	    for (i = 0; i < npaths; i++) {
	      stage3 = stage3array[i];
	      Stage3_free(&stage3);
	    }
	    FREE(stage3array);

	  }
#ifndef PMAP
	}
#endif

	Sequence_free(&queryuc);
      }
      Sequence_free(&queryseq);
      jobid++;
    }
    Sequence_free(&genomicseg);
    Stage3_free(&stage3ref);
  }

  Stage1_free(&stage1,/*completep*/true);

  Pairpool_free(&pairpool);
  FREE(gbuffer3);
  FREE(gbuffer2);
  FREE(gbuffer1);
  Dynprog_free(&dynprogR);
  Dynprog_free(&dynprogM);
  Dynprog_free(&dynprogL);
  Oligoindex_free(&oligoindex);

  return;
}


/*
static void
finish_control_c (int signal) {
  exit(0);
}
*/



int
main (int argc, char *argv[]) {
  IIT_T chromosome_iit = NULL, contig_iit = NULL, altstrain_iit = NULL, map_iit = NULL;
  Chrsubset_T chrsubset = NULL;
#ifdef PMAP
  Indexdb_T indexdb_fwd = NULL, indexdb_rev = NULL;
#else
  Indexdb_T indexdb = NULL;
#endif
  Genome_T genome = NULL;
  Sequence_T queryseq, usersegment = NULL, mutationrefseq = NULL;
  char *genomesubdir = NULL, *mapdir = NULL, *iitfile = NULL, *olddbroot, *fileroot = NULL;
  FILE *input;
  Params_T params;

  int 
    stuttercycles = 2,		/* In stage 1 */
    stutterhits = 3,		/* In stage 1 */
    maxextension = 30000,	/* In stage 1 */
    indexsize,			/* In stage 2 */
    maxpeelback,
    sufflookback = 60,
    nsufflookback = 4,
    nullgap,		/* In stage 3 */
    extramaterial_end,
    extramaterial_paired,
    extraband_single,
    extraband_end,
    extraband_paired,
    user_ngap = -1;
  bool showcontigp = true, multiple_sequences_p = false;
  int nextchar;

#ifdef HAVE_PTHREAD
  int ret, i;
  pthread_t input_thread_id, output_thread_id, worker_thread_id;
  pthread_attr_t thread_attr_detach, thread_attr_join;
  Blackboard_T blackboard;
  Request_T request;
#endif

  int opt;
  extern int optind;
  extern char *optarg;
  int long_option_index = 0;

  while ((opt = getopt_long(argc,argv,
#ifdef PMAP
			    "D:d:R:Gg:C:B:K:L:x:1Xw:t:sc:USA39n:f:ZO5M:m:ebo:v:E:PQNI:i:l:V?",
#else
			    "D:d:R:Gg:C:B:K:L:x:1Xw:t:sc:USA39n:f:ZO5M:m:ebo:v:E:PQFTNI:i:l:V?",
#endif
			    long_options, &long_option_index)) != -1) {
    switch (opt) {
    case 'D': user_genomedir = optarg; break;
    case 'd': 
      dbroot = (char *) CALLOC(strlen(optarg)+1,sizeof(char));
      strcpy(dbroot,optarg);
      break;
    case 'R': releasestring = optarg; break;
    case 'G': uncompressedp = true; break;
    case 'g': user_genomicseg = optarg; break;
    case 'C': user_chrsubsetfile = optarg; break;

    case 'B': switch (atoi(optarg)) {
      case 1: batch_offsets_p = true; batch_positions_p = true; batch_genome_p = false; break;
      case 2: batch_offsets_p = true; batch_positions_p = true; batch_genome_p = true; break;
      default: fprintf(stderr,"Batch mode %s not recognized.\n",optarg);
	fprintf(stderr,"Mode 1 preloads indices only; mode 2 preloads both indices and genome\n");
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
    case 'X': crossspeciesp = true; break;
    case 'w': mutationref = optarg; break;
#ifdef HAVE_PTHREAD
    case 't': nworkers = atoi(optarg); break;
#endif
    case 's': altstrainp = true; break;
    case 'c': user_chrsubsetname = optarg; break;
    case 'U': end_microexons_p = true; break;

    case 'S': summaryonlyp = true; break;
    case 'A': showalignp = true; break;
    case '3': continuousp = true; break;
    case '9': diagnosticp = true; break;
    case 'n': maxpaths = atoi(optarg); break;
    case 'f': switch (atoi(optarg)) {
#ifdef PMAP
      case 0: psloutput_pro_p = true; break;
      case 1: psloutput_nt_p = true; break;
      case 9: print_coordinates_p = true; break;
#else
      case 1: psloutput_nt_p = true; break;
      case 9: print_coordinates_p = true; break;
#endif
      default: fprintf(stderr,"Output format %s not recognized\n",optarg); exit(9);
      }
      break;
    case 'Z': compressoutp = true; break;
    case 'O': orderedp = true; break;
    case '5': checksump = true; break;
    case 'M': user_mapdir = optarg; break;
    case 'm': map_iitfile = optarg; break;
    case 'e': map_exons_p = true; break;
    case 'b': map_bothstrands_p = true; break;
    case 'o': chimera_overlap = atoi(optarg); break;
    case 'v': maxmutations=atoi(optarg); break;

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
#endif
    case 'N': nointronlenp = true; break;
    case 'I': invertmode = atoi(optarg); break;
    case 'i': user_ngap = atoi(optarg); break;
    case 'l': wraplength = atoi(optarg); break;

    case 'V': print_program_version(); exit(0);
    case '?': print_program_usage(); exit(0);
    default: exit(9);
    }
  }
  argc -= optind;
  argv += optind;
      
  if (dbroot == NULL && user_genomicseg == NULL) {
    fprintf(stderr,"Need to specify the -d or the -g flag\n");
    print_program_usage();
    exit(9);
  }

  if (user_ngap >= 0) {
    ngap = user_ngap;
  } else if (cdna_exons_p == true || genomic_exons_p == true) {
    /* If user didn't specify, then set to zero */
    ngap = 0;
  };

  if (crossspeciesp == false) {
#ifdef PMAP
    indexsize = 12;		/* Needs to be a multiple of 3 */
    maxpeelback = 12;		/* Needs to be at least indexsize
				   because stage 2 jumps by indexsize.
				   Also should exceed length of
				   repeated nucleotides (e.g., a
				   string of consecutive T's) */
#else
    indexsize = 8;
    maxpeelback = 11;		/* Needs to be at least indexsize
				   because stage 2 jumps by indexsize.
				   Also should exceed length of
				   repeated nucleotides (e.g., a
				   string of consecutive T's) */
#endif

    if (sufflookback == 0) {
      sufflookback = 60;
      nsufflookback = 5;
    }

    nullgap = 600;
    extramaterial_end = 10;
    extramaterial_paired = 8;	/* Should be at least indexsize */
    extraband_single = 7;	/* Allows program to see indel of 2 aa */
    extraband_end = 7;
    extraband_paired = 7;	/* Allows program to see indel of 2 aa */
  } else {
#ifdef PMAP
    indexsize = 6;		/* Needs to be a multiple of 3 */
    maxpeelback = 6;		/* Needs to be at least indexsize
				   because stage 2 jumps by indexsize.
				   Also should exceed length of
				   repeated nucleotides (e.g., a
				   string of consecutive T's) */
#else
    indexsize = 6;
    maxpeelback = 6;		/* Needs to be at least indexsize because stage 2 jumps by indexsize */
#endif
    if (sufflookback == 0) {
      sufflookback = 120;
      nsufflookback = 4;
    }
    nullgap = 1200;
    extramaterial_end = 20;
    extramaterial_paired = 16;
    extraband_single = 6;
    extraband_end = 7;
    extraband_paired = 3;
  }

  if (maxintronlen_bound > maxtotallen_bound) {
    maxintronlen_bound = maxtotallen_bound;
  }

  /*
  maxlookback = calculate_maxlookback(failureprob,pctidentity);
  fprintf(stderr,"Max lookback = %d + %d for failure prob = %.3e and pctidentity = %.3f\n",
	  maxlookback,lookback_safety,failureprob,pctidentity);
  maxlookback += lookback_safety;
  */


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

  /* Read mutationref before rest of sequences, because of shared usage of sequence.c */
  if (mutationref != NULL) {
    if ((input = FOPEN_READ_TEXT(mutationref)) == NULL) {
      fprintf(stderr,"Can't open file %s\n",mutationref);
      exit(9);
    }
    if ((mutationrefseq = Sequence_read_unlimited(input)) == NULL) {
      fprintf(stderr,"File %s is empty\n",mutationref);
      exit(9);
    }
    fclose(input);
  }


  /* Open query stream */
  if (argc == 0) {
    input = stdin;
  } else if ((input = FOPEN_READ_TEXT(argv[0])) == NULL) {
    fprintf(stderr,"Can't open file %s\n",argv[0]);
    exit(9);
  }

  /* Read first query sequence */
  if ((queryseq = Sequence_read(&nextchar,input,maponlyp)) == NULL) {
    fprintf(stderr,"No input detected\n");
    exit(9);
  } else if (nextchar == '>') {
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
  }

  /* Complement_init(); */
  Oligoindex_init(indexsize);
  Dynprog_init(nullgap,EXTRAQUERYGAP,maxpeelback,extramaterial_end,extramaterial_paired);
#ifdef PMAP
  Backtranslation_init();
#endif

  if (usersegment != NULL) {
    /* Map against user-provided genomic segment */
    showcontigp = false;
    maxpaths = 1;

    genome = (Genome_T) NULL;
    dbversion = (char *) NULL;
    if (userstage1p == true) {
#ifdef PMAP
      indexdb_fwd = Indexdb_new_segment(Sequence_fullpointer(usersegment),/*watsonp*/true,index1interval);
      indexdb_rev = Indexdb_new_segment(Sequence_fullpointer(usersegment),/*watsonp*/false,index1interval);
#else
      indexdb = Indexdb_new_segment(Sequence_fullpointer(usersegment),index1interval);
#endif
    }

  } else {
    /* Map against genome */
    if (releasestring != NULL) {
      olddbroot = dbroot;
      dbroot = (char *) CALLOC(strlen(dbroot)+strlen("_R")+strlen(releasestring)+1,sizeof(char));
      sprintf(dbroot,"%s_R%s",olddbroot,releasestring);
      FREE(olddbroot);
    }
    genomesubdir = Datadir_find_genomesubdir(&fileroot,&dbversion,user_genomedir,dbroot);

    iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			      strlen(fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
    sprintf(iitfile,"%s/%s.chromosome.iit",genomesubdir,fileroot);
    chromosome_iit = IIT_read(iitfile,NULL,true);
    FREE(iitfile);

    chrsubset = Chrsubset_read(user_chrsubsetfile,genomesubdir,fileroot,user_chrsubsetname,
			       chromosome_iit);

    if (showcontigp == true) {
      iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
				strlen(fileroot)+strlen(".contig.iit")+1,sizeof(char));
      sprintf(iitfile,"%s/%s.contig.iit",genomesubdir,fileroot);
      contig_iit = IIT_read(iitfile,NULL,true);
      FREE(iitfile);
    }
  
    if (map_iitfile != NULL) {
      mapdir = Datadir_find_mapdir(user_mapdir,genomesubdir,fileroot);
      iitfile = (char *) CALLOC(strlen(mapdir)+strlen("/")+
				strlen(map_iitfile)+strlen(".iit")+1,sizeof(char));
      sprintf(iitfile,"%s/%s.iit",mapdir,map_iitfile);
      if ((map_iit = IIT_read(iitfile,map_iitfile,true)) == NULL) {
	fprintf(stderr,"Map file %s.iit not found in %s.  Available files:\n",map_iitfile,mapdir);
	Datadir_list_directory(stderr,mapdir);
	fprintf(stderr,"Either install file %s or specify a full directory path\n",map_iitfile);
	fprintf(stderr,"using the -M flag to gmap.\n");
	exit(9);
      }
      FREE(iitfile);
      FREE(mapdir);
    }

#ifdef PMAP
    indexdb_fwd = Indexdb_new_genome(genomesubdir,fileroot,/*watsonp*/true,
				     batch_offsets_p,batch_positions_p);
    indexdb_rev = Indexdb_new_genome(genomesubdir,fileroot,/*watsonp*/false,
				     batch_offsets_p,batch_positions_p);
#else
    indexdb = Indexdb_new_genome(genomesubdir,fileroot,batch_offsets_p,batch_positions_p);
#endif
    genome = Genome_new(genomesubdir,fileroot,uncompressedp,batch_genome_p);

    if (altstrainp == true) {
      iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
				strlen(fileroot)+strlen(".altstrain.iit")+1,sizeof(char));
      sprintf(iitfile,"%s/%s.altstrain.iit",genomesubdir,fileroot);
      global_altstrain_iit = altstrain_iit = IIT_read(iitfile,NULL,true);
      if (altstrain_iit == NULL) {
	fprintf(stderr,"Unable to find file %s.  Include strain information during gmap_setup process, or run gmap without the -s flag.\n",iitfile);
	exit(9);
      }
      FREE(iitfile);
    }

    FREE(genomesubdir);
    FREE(fileroot);
    FREE(dbroot);
  }

  params = Params_new(genome,altstrain_iit,
#ifdef PMAP
		      indexdb_fwd,indexdb_rev,
#else
		      indexdb,
#endif
		      chromosome_iit,chrsubset,
		      contig_iit,map_iit,maxextension,stuttercycles,stutterhits,
		      indexsize,maxpeelback,sufflookback,nsufflookback,nullgap,
		      extramaterial_end,extramaterial_paired,
		      extraband_single,extraband_end,extraband_paired,
		      maxmutations);

  if (mutationrefseq != NULL) {
    chimera_margin = -1;
    align_mutations(input,nextchar,queryseq,params,mutationrefseq);
    Sequence_free(&mutationrefseq);

  } else {
#ifndef HAVE_PTHREAD
    single_thread(input,nextchar,queryseq,params,usersegment);
#else
    if (multiple_sequences_p == false) {
      single_thread(input,nextchar,queryseq,params,usersegment);
    } else {
      /* Make blackboard and threads */
      blackboard = Blackboard_new(input,nextchar,usersegment,nworkers,params);
      debug(printf("input_thread: Putting request id 0\n"));
      request = Request_new(0,queryseq,usersegment);
      Blackboard_put_request(blackboard,request);

      pthread_attr_init(&thread_attr_detach);
      if ((ret = pthread_attr_setdetachstate(&thread_attr_detach,PTHREAD_CREATE_DETACHED)) != 0) {
	fprintf(stderr,"ERROR: pthread_attr_setdetachstate\n");
	exit(1);
      }
      pthread_attr_init(&thread_attr_join);
      if ((ret = pthread_attr_setdetachstate(&thread_attr_join,PTHREAD_CREATE_JOINABLE)) != 0) {
	fprintf(stderr,"ERROR: pthread_attr_setdetachstate\n");
	exit(1);
      }
    
      for (i = 0; i < nworkers; i++) {
	pthread_create(&worker_thread_id,&thread_attr_detach,worker_thread,
		       (void *) Blackboard_get_reqpost(blackboard,i));
      }
      if (orderedp == true) {
	pthread_create(&output_thread_id,&thread_attr_join,output_thread_ordered,
		       (void *) blackboard);
      } else {
	pthread_create(&output_thread_id,&thread_attr_join,output_thread_anyorder,
		       (void *) blackboard);
      }
      pthread_create(&input_thread_id,&thread_attr_detach,input_thread,(void *) blackboard);
    
      pthread_join(output_thread_id,NULL);
      Blackboard_free(&blackboard);
    }

#endif /* HAVE_PTHREAD */
  }

  Params_free(&params);
  
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
