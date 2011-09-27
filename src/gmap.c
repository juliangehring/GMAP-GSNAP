static char rcsid[] = "$Id: gmap.c,v 1.255 2005/06/23 22:49:12 twu Exp $";
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
#include "sequence.h"
#include "oligoindex.h"
#include "match.h"
#include "matchpair.h"
#include "genome.h"
#include "stage1.h"
#include "stage2.h"
#include "dynprog.h"
#include "stage3.h"
#ifdef BETATEST
#include "chimera.h"
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

#define OLIGO_DEPTH_SEQLENGTH 10000 /* Substantial length relative to 4^8 = 65536 */
#define MAX_OLIGO_DEPTH 2.0
#define MAX_BADOLIGOS 1.0	/* Setting to 1.0 effectively turns this check off */
#define MIN_STAGE1_SUPPORT 0.10
#define MAX_STAGE1_STRETCH 2000.0

#define CHIMERA_IDENTITY 0.98
#define CHIMERA_PVALUE 0.01
#define CHIMERA_FVALUE 6.634897	/* qnorm(CHIMERA_PVALUE/2)^2 */
#define CHIMERA_HANDICAP 10	/* points, for minor alignment differences based on different defect rates */

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
static bool batchp = false;
static int maxintronlen_bound = 1200000;
static int chimera_margin = -1;
static bool maponlyp = false;
static bool userstage1p = true;	/* Apply stage 1 for user-provided genomic segments */
static int index1interval = 6;	/* Stage 1 interval if user provides a genomic segment */
static bool crossspeciesp = false;
static char *mutationref = NULL;
static bool literalrefp = false;
static bool altstrainp = false;
static bool nworkers = 1;	/* (int) sysconf(_SC_NPROCESSORS_ONLN) */
static bool extend_mismatch_p = false;
static bool end_microexons_p = false;
static char *user_chrsubsetname = NULL;


/* Output options */
static bool summaryonlyp = false;
static bool showalignp = false;
static bool continuousp = false;
static bool diagnosticp = false;
static int maxpaths = 5;	/* 0 means 1 if nonchimeric, 2 if chimeric */
static bool compressoutp = false;
static bool orderedp = false;
static bool checksump = false;
static char *user_mapdir = NULL;
static char *map_iitfile = NULL;
static bool map_bothstrands_p = false;
static int chimera_overlap = 0;

/* Alignment options */
static bool cdna_exons_p = false;
static bool protein_cdna_p = false;
static bool protein_genomic_p = false;
static bool fulllengthp = false;
static bool truncatep = false;
static int proteinmode = 1;
static bool uncompressedp = false;
static bool nointronlenp = false;
static int invertmode = 0;
static int ngap = 3;
static int wraplength = 50;

/* Used alphabetically: 1359ABbCcDdEeFgIiLlMmNnOoPQqRSsTtUVwXxZ */

static struct option long_options[] = {
  /* Input options */
  {"dir", required_argument, 0, 'D'},	/* user_genomedir */
  {"db", required_argument, 0, 'd'}, /* dbroot */
  {"rel", required_argument, 0, 'R'}, /* releasestring */
  {"gseg", required_argument, 0, 'g'}, /* user_genomicseg */

  /* Compute options */
  {"batch", no_argument, 0, 'B'}, /* batchp */
  {"length", required_argument, 0, 'L'}, /* maxintronlen_bound */
  {"chimera_margin", required_argument, 0, 'x'}, /* chimera_margin */
  {"maponly", no_argument, 0, '1'}, /* maponlyp */
  {"indexinterval", required_argument, 0, 'q'},	/* index1interval */
  {"cross", no_argument, 0, 'X'}, /* crossspeciesp */
  {"mutationref", required_argument, 0, 'w'},	/* mutationref */
  {"nthreads", required_argument, 0, 't'}, /* nworkers */
  {"altstrain", no_argument, 0, 's'},	/* altstrainp */
  {"extend", no_argument, 0, 'e'}, /* extend_mismatch_p */
  {"chrsubsetfile", required_argument, 0, 'C'}, /* user_chrsubsetfile */
  {"chrsubset", required_argument, 0, 'c'}, /* user_chrsubsetname */
  {"endmicro", no_argument, 0, 'U'}, /* end_microexons_p */

  /* Output options */
  {"summary", no_argument, 0, 'S'}, /* summaryonlyp */
  {"align", no_argument, 0, 'A'}, /* showalignp */
  {"continuous", no_argument, 0, '3'}, /* continuousp */
  {"diagnostic", no_argument, 0, '9'}, /* diagnosticp */
  {"npaths", required_argument, 0, 'n'}, /* maxpaths */
  {"compress", no_argument, 0, 'Z'}, /* compressoutp */
  {"ordered", no_argument, 0, 'O'}, /* orderedp */
  {"md5", no_argument, 0, '5'}, /* checksump */
  {"mapdir", required_argument, 0, 'M'}, /* user_mapdir */
  {"map", required_argument, 0, 'm'},	/* map_iitfile */
  {"mapboth", no_argument, 0, 'b'}, /* map_bothstrands_p */
  {"chimera_overlap", required_argument, 0, 'o'}, /* chimera_overlap */

  /* Alignment options */
  {"exons", no_argument, 0, 'E'}, /* cdna_exons_p */
  {"protein_dna", no_argument, 0, 'P'}, /* protein_cdna_p */
  {"protein_gen", no_argument, 0, 'Q'}, /* protein_genomic_p */
  {"fulllength", no_argument, 0, 'F'}, /* fulllengthp */
  {"truncate", no_argument, 0, 'T'}, /* truncatep */
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
  fprintf(stdout,"GMAP: Genomic Mapping and Alignment Program\n");
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
  fprintf(stdout,"Thomas D. Wu and Colin K. Watanabe, Genentech, Inc.\n");
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
  -g, --gseg=filename            User-suppled genomic segment\n\
\n\
Computation options\n\
  -B, --batch                    Batch mode (pre-loads genomic index files)\n\
  -L, --length                   Max total intron length (default 1200000)\n\
  -x, --chimera_margin=INT       Amount of unaligned sequence that triggers\n\
                                 search for a chimera (default off)\n\
  -1, --maponly                  Report stage 1 only (genomic mapping only)\n\
  -X, --cross                    Cross-species mode (changes parameters)\n\
  -w, --mutationref=filename     Mutation reference\n\
  -t, --nthreads=INT             Number of worker threads\n\
  -s, --altstrain                Search alternate strains in addition\n\
  -e, --extend                   Extend alignment to show mismatch at end\n\
  -C, --chrsubsetfile=filename   User-suppled chromosome subset file\n\
  -c, --chrsubset=string         Chromosome subset to search\n\
  -U, --endmicro                 Allow end microexons with arbitrarily long introns\n\
\n\
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
  -b, --mapboth                  Map both strands of genome\n\
  -o, --chimera_overlap          Overlap to show, if any, at chimera breakpoint\n\
\n\
Alignment output options\n\
  -E, --exons                    Print cDNA exons\n\
  -P, --protein_dna              Print protein sequence (cDNA)\n\
  -Q, --protein_gen              Print protein sequence (genomic)\n\
  -F, --fulllength               Assume full-length protein, starting with Met\n\
  -T, --truncate                 Truncate full-length protein, from Met to Stop\n\
                                 Implies -F flag.\n\
  -N, --nolengths                No intron lengths in alignment\n\
  -I, --invertmode=INT           Mode for alignments to (-) strand: 1=Invert\n\
                                 and print sense (-) strand; 2=Invert and\n\
                                 print antisense (+) strand.\n\
  -i, --introngap=INT            Nucleotides to show on each end of intron\n\
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
print_path_summaries (Stage3_T *stage3array, int npaths, Failure_T failuretype,
		      IIT_T chromosome_iit, IIT_T contig_iit, 
		      char *dbversion, int maxpaths, Chimera_T chimera, bool zerobasedp) {
  int pathnum;
  Stage3_T new, old;

  if (npaths == 0) {
    printf("Paths (0):");
    if (failuretype == EMPTY_SEQUENCE) {
      printf(" *** Empty sequence ***");
    } else if (failuretype == REPETITIVE) {
      printf(" *** Repetitive sequence ***");
    }
    printf("\n\n");
  } else {
    printf("Paths (%d):",npaths);
#ifdef BETATEST
    if (chimera != NULL) {
      Chimera_print(chimera);
    }
#endif
    printf("\n");
    if (maxpaths == 0) {
      /* Special mode */
      if (npaths > 0) {
	stage3array[0] = Stage3_translate_genomic(stage3array[0],fulllengthp,truncatep);
	Stage3_print_pathsummary(stage3array[0],/*pathnum*/1,chromosome_iit,contig_iit,
				 dbversion,zerobasedp);
	if (chimera != NULL && npaths > 1) {
	  stage3array[1] = Stage3_translate_genomic(stage3array[1],fulllengthp,truncatep);
	  Stage3_print_pathsummary(stage3array[1],/*pathnum*/2,chromosome_iit,contig_iit,
				   dbversion,zerobasedp);
	}
      }
    } else {
      for (pathnum = 1; pathnum <= npaths && pathnum <= maxpaths; pathnum++) {
	stage3array[pathnum-1] = Stage3_translate_genomic(stage3array[pathnum-1],fulllengthp,truncatep);
	Stage3_print_pathsummary(stage3array[pathnum-1],pathnum,chromosome_iit,contig_iit,
				 dbversion,zerobasedp);
      }
    }
  }

  return;
}


#ifdef HAVE_PTHREAD
static void *
input_thread (void *data) {
  Blackboard_T blackboard = (Blackboard_T) data;
  Params_T params = Blackboard_params(blackboard);
  FILE *input = Blackboard_input(blackboard);
  Sequence_T queryseq, usersegment;
  Request_T request;
  int inputid = 0;
  int nextchar;

  debug(printf("input_thread: Starting\n"));
  usersegment = Blackboard_usersegment(blackboard);
  while ((queryseq = Sequence_read(&nextchar,input)) != NULL) {
    if (inputid == 0 && nextchar == '>' && batchp == false) {
      fprintf(stderr,"Note: batch mode (-B) is strongly recommended for large runs.\n");
    }
    debug(printf("input_thread: Putting request id %d\n",inputid));
    request = Request_new(inputid++,queryseq,usersegment);
    Blackboard_put_request(blackboard,request);
  }

  Blackboard_set_inputdone(blackboard);
  debug(printf("input_thread: Ending\n"));
  return NULL;
}
#endif

void
print_maponly (List_T matchlist, int npaths, IIT_T chromosome_iit, IIT_T contig_iit,
	       char *dbversion, int maxpaths, bool zerobasedp) {
  Matchpair_T matchpair;
  List_T p;
  int pathnum = 1;

  printf("Paths (%d):\n",npaths);
  if (maxpaths == 0) {
    /* Chimeras not computed in this case, so treat as 1 */
    maxpaths = 1;
  }
  for (p = matchlist; p != NULL && pathnum <= maxpaths; p = List_next(p)) {
    matchpair = List_head(p);
    Match_print_two(pathnum++,Matchpair_bound5(matchpair),Matchpair_bound3(matchpair),
		    chromosome_iit,contig_iit,dbversion,zerobasedp);
  }
  return;
}

static void
print_result (Result_T result, Request_T request, Params_T params) {
  Sequence_T queryseq;
  List_T matchlist;
  Stage1_T stage1;
  Stage3_T *stage3array, stage3;
  int npaths, pathnum;
  Chimera_T chimera;
  int chimerapos, chimeraequivpos, chimera_cdna_direction;
  double donor_prob, acceptor_prob;

  queryseq = Request_queryseq(request);
  if (compressoutp == false && cdna_exons_p == false &&
      protein_cdna_p == false && protein_genomic_p == false) {
    Sequence_print_header(queryseq,checksump);
  }

  if (maponlyp == true) {
    if ((stage1 = Result_stage1(result)) == NULL) {
      printf("Paths (0):"); 
      if (Result_failuretype(result) == EMPTY_SEQUENCE) {
	printf(" *** Empty sequence ***");
      } else if (Result_failuretype(result) == REPETITIVE) {
	printf(" *** Repetitive sequence ***");
      }
     printf("\n\n");
    } else {
      matchlist = Stage1_matchlist(stage1,Params_indexdb(params),Params_chromosome_iit(params),
				   Params_chrsubset(params),crossspeciesp);
      npaths = List_length(matchlist);
      if (npaths == 0) {
	printf("Paths (0):\n\n");
      } else {
	print_maponly(matchlist,npaths,Params_chromosome_iit(params),Params_contig_iit(params),
		      dbversion,maxpaths,/*zerobasedp*/false);
      }
    }
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

    stage3array = Result_array(&npaths,result);
    if (maxpaths == 0) {
      /* Special mode */
      if (npaths > 0) {
	if (truncatep == true) {
	  stage3array[0] = Stage3_truncate_fulllength(stage3array[0],/*translatep*/true);
	}
	Stage3_print_compressed(stage3array[0],queryseq,Params_chromosome_iit(params),
				dbversion,/*pathnum*/1,npaths,checksump,chimerapos,chimeraequivpos,
				donor_prob,acceptor_prob,chimera_cdna_direction,/*zerobasedp*/false);
	if (chimera != NULL && npaths > 1) {
	  if (truncatep == true) {
	    stage3array[1] = Stage3_truncate_fulllength(stage3array[1],/*translatep*/true);
	  }
	  Stage3_print_compressed(stage3array[1],queryseq,Params_chromosome_iit(params),
				  dbversion,/*pathnum*/2,npaths,checksump,chimerapos,chimeraequivpos,
				  donor_prob,acceptor_prob,chimera_cdna_direction,/*zerobasedp*/false);
	}
      }
    } else {
      for (pathnum = 1; pathnum <= npaths && pathnum <= maxpaths; pathnum++) {
	if (truncatep == true) {
	  stage3array[pathnum-1] = Stage3_truncate_fulllength(stage3array[pathnum-1],/*translatep*/true);
	}
	Stage3_print_compressed(stage3array[pathnum-1],queryseq,Params_chromosome_iit(params),
				dbversion,pathnum,npaths,checksump,chimerapos,chimeraequivpos,
				donor_prob,acceptor_prob,chimera_cdna_direction,/*zerobasedp*/false);
      }
    }

  } else if (continuousp == true) {
    stage3array = Result_array(&npaths,result);
    if (npaths == 0) {
      printf("\n\n\n");
    } else {
      stage3 = stage3array[0];
      Stage3_print_alignment(stage3,Params_chromosome_iit(params),summaryonlyp,/*universalp*/false,
			     /*zerobasedp*/false,continuousp,diagnosticp,/*flipgenomep*/true,
			     proteinmode,invertmode,nointronlenp,wraplength);
    }

  } else if (cdna_exons_p == true) {
    stage3array = Result_array(&npaths,result);
    Sequence_print_header(queryseq,checksump);
    Stage3_print_cdna_exons(stage3array[0],wraplength);

  } else if (protein_cdna_p == true) {
    stage3array = Result_array(&npaths,result);
    if (maxpaths == 0) {
      /* Special mode */
      if (npaths > 0) {
	Sequence_print_header(queryseq,checksump);
	stage3array[0] = Stage3_translate_genomic(stage3array[0],fulllengthp,truncatep);
	Stage3_print_protein_cdna(stage3array[0],wraplength);
	if (Result_chimera(result) != NULL && npaths > 1) {
	  Sequence_print_header(queryseq,checksump);
	  stage3array[1] = Stage3_translate_genomic(stage3array[1],fulllengthp,truncatep);
	  Stage3_print_protein_cdna(stage3array[1],wraplength);
	}
      }
    } else {
      for (pathnum = 0; pathnum < npaths && pathnum < maxpaths; pathnum++) {
	Sequence_print_header(queryseq,checksump);
	stage3array[pathnum] = Stage3_translate_genomic(stage3array[pathnum],fulllengthp,truncatep);
	Stage3_print_protein_cdna(stage3array[pathnum],wraplength);
      }
    }

  } else if (protein_genomic_p == true) {
    stage3array = Result_array(&npaths,result);
    if (maxpaths == 0) {
      /* Special mode */
      if (npaths > 0) {
	Sequence_print_header(queryseq,checksump);
	stage3array[0] = Stage3_translate_genomic(stage3array[0],fulllengthp,truncatep);
	Stage3_print_protein_genomic(stage3array[0],wraplength);
	if (Result_chimera(result) != NULL && npaths > 1) {
	  Sequence_print_header(queryseq,checksump);
	  stage3array[1] = Stage3_translate_genomic(stage3array[1],fulllengthp,truncatep);
	  Stage3_print_protein_genomic(stage3array[1],wraplength);
	}
      }
    } else {
      for (pathnum = 0; pathnum < npaths && pathnum < maxpaths; pathnum++) {
	Sequence_print_header(queryseq,checksump);
	stage3array[pathnum] = Stage3_translate_genomic(stage3array[pathnum],fulllengthp,truncatep);
	Stage3_print_protein_genomic(stage3array[pathnum],wraplength);
      }
    }

  } else {
    /* Usual output */
    stage3array = Result_array(&npaths,result);
    if (npaths == 0) {
      printf("Paths (0):");
      if (Result_failuretype(result) == EMPTY_SEQUENCE) {
	printf(" *** Empty sequence ***");
      } else if (Result_failuretype(result) == REPETITIVE) {
	printf(" *** Repetitive sequence ***");
      }
      printf("\n\n");
    } else {
      chimera = Result_chimera(result);
      print_path_summaries(stage3array,npaths,Result_failuretype(result),
			   Params_chromosome_iit(params),Params_contig_iit(params),
			   dbversion,maxpaths,chimera,/*zerobasedp*/false);
    
      if (Params_map_iit(params) != NULL) {
	printf("Maps:\n");
	if (maxpaths == 0) {
	  /* Special mode */
	  if (npaths > 0) {
	    Stage3_print_map(stage3array[0],Params_map_iit(params),/*pathnum*/1,map_bothstrands_p);
	    if (chimera != NULL && npaths > 1) {
	      Stage3_print_map(stage3array[1],Params_map_iit(params),/*pathnum*/2,map_bothstrands_p);
	    }
	  }
	} else {
	  for (pathnum = 1; pathnum <= npaths && pathnum <= maxpaths; pathnum++) {
	    Stage3_print_map(stage3array[pathnum-1],Params_map_iit(params),pathnum,map_bothstrands_p);
	  }
	}
      }

      if (summaryonlyp == true || showalignp == true) {
	printf("Alignments:\n");
	if (maxpaths == 0) {
	  /* Special mode */
	  if (npaths > 0) {
	    printf("  Alignment for path 1:\n\n");
	    Stage3_print_alignment(stage3array[0],Params_chromosome_iit(params),summaryonlyp,/*universalp*/false,
				   /*zerobasedp*/false,continuousp,diagnosticp,
				   /*flipgenomep*/true,proteinmode,invertmode,nointronlenp,wraplength);
	    if (chimera != NULL && npaths > 1) {
	      printf("  Alignment for path 2:\n\n");
	      Stage3_print_alignment(stage3array[1],Params_chromosome_iit(params),summaryonlyp,/*universalp*/false,
				     /*zerobasedp*/false,continuousp,diagnosticp,
				     /*flipgenomep*/true,proteinmode,invertmode,nointronlenp,wraplength);
	    }
	  }
	} else {
	  for (pathnum = 1; pathnum <= npaths && pathnum <= maxpaths; pathnum++) {
	    printf("  Alignment for path %d:\n\n",pathnum);
	    Stage3_print_alignment(stage3array[pathnum-1],Params_chromosome_iit(params),summaryonlyp,/*universalp*/false,
				   /*zerobasedp*/false,continuousp,diagnosticp,
				   /*flipgenomep*/true,proteinmode,invertmode,nointronlenp,wraplength);
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
  return NULL;
}
#endif

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

#ifdef HAVE_PTHREAD
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
  return NULL;
}
#endif
    
static void
output_one (Result_T result, Request_T request, Params_T params) {
  print_result(result,request,params);
  Result_free(&result);
  Request_free(&request);
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
update_stage3list (List_T stage3list, Sequence_T queryseq, Sequence_T genomicseg, 
		   Matchpair_T matchpair, int straintype, char *strain, 
		   Oligoindex_T oligoindex, int indexsize, Pairpool_T pairpool, 
		   int sufflookback, int nsufflookback, int badoligos, Chrnum_T chrnum,
		   Genomicpos_T chroffset, Genomicpos_T chrpos, bool watsonp, int maxpeelback, int nullgap,
		   int extramaterial_end, int extramaterial_paired, int extraband_single, int extraband_end, 
		   int extraband_paired, int ngap, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
		   IIT_T altstrain_iit) {
  Stage2_T stage2;
  Stage3_T stage3;
  Matchpairend_T matchpairend;

  if (matchpair == NULL) {
    matchpairend = MIXED;
  } else {
    debug1(printf("support: %d/%d\n",Matchpair_support(matchpair),Sequence_trimlength(queryseq)));
    debug1(printf("stretch: %f\n",Matchpair_stretch(matchpair)));
    if ((double) Matchpair_support(matchpair)/(double) Sequence_trimlength(queryseq) < MIN_STAGE1_SUPPORT) {
      /* Insufficient coverage of the query sequence */
      return stage3list;
    } else if (Matchpair_stretch(matchpair) > MAX_STAGE1_STRETCH) {
      /* Genomic region is too large relative to matching cDNA region */
      return stage3list;
    } else {
      matchpairend = Matchpair_matchpairend(matchpair);
    }
  }

  if ((stage2 = Stage2_compute(queryseq,genomicseg,oligoindex,indexsize,
			       pairpool,sufflookback,nsufflookback,badoligos,crossspeciesp)) != NULL) {
    if ((stage3 = Stage3_compute(stage2,queryseq,genomicseg,matchpairend,straintype,strain,chrnum,chroffset,
				 chrpos,watsonp,maxpeelback,nullgap,extramaterial_end,extramaterial_paired,
				 extraband_single,extraband_end,extraband_paired,
				 ngap,extend_mismatch_p,end_microexons_p,pairpool,dynprogL,dynprogM,dynprogR,
				 altstrain_iit)) != NULL) {
      stage3list = List_push(stage3list,stage3);
    }
    Stage2_free(&stage2);
  }
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
stage3_from_usersegment (int *npaths, Sequence_T queryseq, Sequence_T usersegment, Params_T params, 
			 int maxpeelback, int sufflookback, int nsufflookback, 
			 int badoligos, int nullgap, 
			 Oligoindex_T oligoindex, Pairpool_T pairpool, 
			 Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
			 char *gbuffer1, char *gbuffer2, char *gbuffer3, int gbufferlen) {
  List_T stage3list = NULL, p;
  Intlist_T indices;
  Genomicpos_T chroffset, chrlength, chrpos, genomicstart, genomiclength;
  bool watsonp;
  Sequence_T revcomp;
  int indexsize = Params_indexsize(params);
  int extramaterial_end = Params_extramaterial_end(params);
  int extramaterial_paired = Params_extramaterial_paired(params);
  int extraband_single = Params_extraband_single(params);
  int extraband_end = Params_extraband_end(params);
  int extraband_paired = Params_extraband_paired(params);
  Interval_T interval;
  int index, *indexarray, nindices, type, i, j;
  char *strain = NULL;
  Chrnum_T chrnum = 0;
  IIT_T altstrain_iit = NULL;
		    
  chroffset = chrpos = 0U;
  genomicstart = 0U;
  genomiclength = chrlength = Sequence_fulllength(usersegment);
  watsonp = true;
  strain = NULL;

  stage3list = update_stage3list(stage3list,queryseq,usersegment,(Matchpair_T) NULL,
				 0,NULL,oligoindex,indexsize,pairpool,
				 sufflookback,nsufflookback,badoligos,chrnum,chroffset,chrpos,
				 watsonp,maxpeelback,nullgap,extramaterial_end,extramaterial_paired,
				 extraband_single,extraband_end,extraband_paired,
				 ngap,dynprogL,dynprogM,dynprogR,altstrain_iit);
	
  revcomp = Sequence_revcomp(usersegment);
  watsonp = false;
  stage3list = update_stage3list(stage3list,queryseq,revcomp,(Matchpair_T) NULL,
				 0,NULL,oligoindex,indexsize,pairpool,
				 sufflookback,nsufflookback,badoligos,chrnum,chroffset,chrpos,
				 watsonp,maxpeelback,nullgap,extramaterial_end,extramaterial_paired,
				 extraband_single,extraband_end,extraband_paired,
				 ngap,dynprogL,dynprogM,dynprogR,altstrain_iit);
  Sequence_free(&revcomp);

  if (stage3list == NULL) {
    *npaths = 0;
    return NULL;
  } else {
    return stage3array_from_list(&(*npaths),stage3list);
  }
}


static Stage3_T *
stage3_from_matchlist (int *npaths, Stage3_T *oldstage3array, Stage1_T stage1,
		       List_T matchlist, Sequence_T queryseq, Sequence_T usersegment, Params_T params, 
		       int maxpeelback, int sufflookback, int nsufflookback, 
		       int badoligos, int nullgap, 
		       Oligoindex_T oligoindex, Pairpool_T pairpool, 
		       Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
		       char *gbuffer1, char *gbuffer2, char *gbuffer3, int gbufferlen) {
  List_T stage3list = NULL, p;
  Matchpair_T matchpair;
  Match_T start, end;
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
  IIT_T altstrain_iit = Params_altstrain_iit(params);
  IIT_T chromosome_iit = Params_chromosome_iit(params);
  Interval_T interval;
  int index, *indexarray, nindices, type, i, j;
		    
  for (i = 0; i < *npaths; i++) {
    stage3list = List_push(stage3list,(void *) oldstage3array[i]);
  }
  if (oldstage3array != NULL) {
    FREE(oldstage3array);
  }

  for (p = matchlist; p != NULL; p = List_next(p)) {
    matchpair = List_head(p);
    start = Matchpair_bound5(matchpair);
    chrnum = Match_chrnum(start);
    if (usersegment != NULL) {
      chrlength = Sequence_fulllength(usersegment);
      chroffset = 0U;
      Matchpair_get_coords(&chrpos,&genomicstart,&genomiclength,&watsonp,matchpair,stage1,queryseq,
			   chrlength,Params_maxextension(params));
      genomicseg = Sequence_substring(usersegment,genomicstart,genomiclength,!watsonp,
				      gbuffer1,gbuffer2,gbufferlen);
      strain = NULL;

    } else {
      index = IIT_get_one(chromosome_iit,Match_position(start),Match_position(start));
      interval = IIT_interval(chromosome_iit,index);
      chrlength = Interval_length(interval);
      chroffset = Interval_low(interval);
      Matchpair_get_coords(&chrpos,&genomicstart,&genomiclength,&watsonp,matchpair,stage1,queryseq,
			   chrlength,Params_maxextension(params));
      genomicseg = Genome_get_segment(Params_genome(params),genomicstart,genomiclength,!watsonp,
				      gbuffer1,gbuffer2,gbufferlen);
      if (altstrain_iit == NULL) {
	strain = NULL;
      } else {
	strain = Params_refstrain(params);
      }
    }

    stage3list = update_stage3list(stage3list,queryseq,genomicseg,matchpair,
				   0,strain,oligoindex,indexsize,pairpool,
				   sufflookback,nsufflookback,badoligos,chrnum,chroffset,chrpos,
				   watsonp,maxpeelback,nullgap,extramaterial_end,extramaterial_paired,
				   extraband_single,extraband_end,extraband_paired,
				   ngap,dynprogL,dynprogM,dynprogR,altstrain_iit);
    Sequence_free(&genomicseg);

    /* We rely upon the fact that gbuffer1 still holds the genomic segment.  This code is duplicated in get-genome.c */
    if (altstrain_iit != NULL) {
      indexarray = IIT_get(&nindices,altstrain_iit,genomicstart+1,genomicstart+genomiclength-1);
      if (nindices > 0) {
	/* Sort according to type and genome position */
	qsort(indexarray,nindices,sizeof(int),index_compare);
	j = 0;
	while (j < nindices) {
	  i = j++;
	  type = Interval_type(IIT_interval(altstrain_iit,indexarray[i]));
	  strain = IIT_typestring(altstrain_iit,type);
	  while (j < nindices && Interval_type(IIT_interval(altstrain_iit,indexarray[j])) == type) {
	    j++;
	  }
	  /* Patch from i to j */
	  genomicseg = Genome_patch_strain(&(indexarray[i]),j-i,altstrain_iit,
					   genomicstart,genomiclength,!watsonp,
					   gbuffer1,gbuffer2,gbuffer3,gbufferlen);
	  stage3list = update_stage3list(stage3list,queryseq,genomicseg,matchpair,
					 type,strain,oligoindex,indexsize,pairpool,
					 sufflookback,nsufflookback,badoligos,chrnum,chroffset,chrpos,
					 watsonp,maxpeelback,nullgap,extramaterial_end,extramaterial_paired,
					 extraband_single,extraband_end,extraband_paired,
					 ngap,dynprogL,dynprogM,dynprogR,altstrain_iit);
	  Sequence_free(&genomicseg);
	}
	FREE(indexarray);
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
static Stage3_T *
check_for_chimera (Chimera_T *chimera, int *npaths, int breakpoint, Stage3_T nonchimericbest, 
		   Stage3_T *stage3array, Sequence_T queryseq, Sequence_T usersegment, 
		   Params_T params, int maxpeelback, int sufflookback, int nsufflookback, int badoligos, int nullgap,
		   Oligoindex_T oligoindex, Pairpool_T pairpool,
		   Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
		   char *gbuffer1, char *gbuffer2, char *gbuffer3, int gbufferlen) {
  int substart, subend;
  Stage1_T stage1sub;
  List_T matchlist = NULL;
  Stage3_T *stage3array_sub, stage3, best0, best1;
  Sequence_T querysubseq1 = NULL, querysubseq2 = NULL;
  int chimerapos, chimeraequivpos;
  int npaths_sub, bestfrom, bestto, matches0, matches1;
  int newj, oldj, j, querylength;
  double donor_prob, acceptor_prob;

  substart = Sequence_trim_start(queryseq);
  subend = Sequence_trim_end(queryseq);

  /* printf("%d %d %d\n",substart,breakpoint,subend); */

  /* Need to check both margins */
  if (breakpoint - substart < chimera_margin ||
      subend - breakpoint < chimera_margin ||
      (querysubseq1 = Sequence_subsequence(queryseq,substart,breakpoint)) == NULL ||
      (querysubseq2 = Sequence_subsequence(queryseq,breakpoint,subend)) == NULL) {
    debug2(printf("Insufficient margins\n"));
    *chimera = (Chimera_T) NULL;
    if (querysubseq1 != NULL) {
      Sequence_free(&querysubseq1);
    }
    /* The following should not be necessary; if we got here, then querysubseq2 is NULL */
    if (querysubseq2 != NULL) {
      Sequence_free(&querysubseq2);
    }

  } else {
    debug2(Sequence_print(querysubseq1,/*uppercasep*/true,wraplength,/*trimmedp*/true));
    debug2(Sequence_print(querysubseq2,/*uppercasep*/true,wraplength,/*trimmedp*/true));

    *chimera = (Chimera_T) NULL;
    npaths_sub = 0;

    stage1sub = Stage1_compute(querysubseq1,Params_indexdb(params),Params_chromosome_iit(params),
			       Params_chrsubset(params),maxintronlen_bound,Params_stuttercycles(params),
			       Params_stutterhits(params),crossspeciesp);
    matchlist = Stage1_matchlist(stage1sub,Params_indexdb(params),Params_chromosome_iit(params),
				 Params_chrsubset(params),crossspeciesp);
    stage3array_sub = stage3_from_matchlist(&npaths_sub,(Stage3_T *) NULL,stage1sub,matchlist,querysubseq1,
					    usersegment,params,maxpeelback,sufflookback,nsufflookback,
					    badoligos,nullgap,
					    oligoindex,pairpool,dynprogL,dynprogM,dynprogR,gbuffer1,gbuffer2,
					    gbuffer3,gbufferlen);
    Stage1_free(&stage1sub);
    Sequence_free(&querysubseq1);

    stage1sub = Stage1_compute(querysubseq2,Params_indexdb(params),Params_chromosome_iit(params),
			       Params_chrsubset(params),maxintronlen_bound,Params_stuttercycles(params),
			       Params_stutterhits(params),crossspeciesp);
    matchlist = Stage1_matchlist(stage1sub,Params_indexdb(params),Params_chromosome_iit(params),
				 Params_chrsubset(params),crossspeciesp);
    stage3array_sub = stage3_from_matchlist(&npaths_sub,stage3array_sub,stage1sub,matchlist,querysubseq2,
					    usersegment,params,maxpeelback,sufflookback,nsufflookback,
					    badoligos,nullgap,
					    oligoindex,pairpool,dynprogL,dynprogM,dynprogR,gbuffer1,gbuffer2,
					    gbuffer3,gbufferlen);
    Stage1_free(&stage1sub);
    Sequence_free(&querysubseq2);

    if (npaths_sub == 0) {
      debug2(printf("npaths_sub == 0\n"));
      *chimera = (Chimera_T) NULL;

    } else if (npaths_sub == 1) {
      debug2(printf("npaths_sub == 1\n"));
      *chimera = (Chimera_T) NULL;
      stage3 = stage3array_sub[0];
      Stage3_free(&stage3);
      FREE(stage3array_sub);

    } else {
      querylength = Sequence_fulllength(queryseq);
      Chimera_bestpath(&chimerapos,&chimeraequivpos,&bestfrom,&bestto,stage3array_sub,npaths_sub,querylength);
      
      if (bestfrom == bestto || 
	  Stage3_chimeric_goodness(&matches0,&matches1,stage3array_sub[bestfrom],stage3array_sub[bestto],
				   chimerapos,querylength) <= Stage3_goodness(nonchimericbest) + CHIMERA_HANDICAP ||
	  matches0 < chimera_margin ||
	  matches1 < chimera_margin) {
	debug2(printf("bestfrom = %d, bestto = %d\n",bestfrom,bestto));
	*chimera = (Chimera_T) NULL;
	for (j = 0; j < npaths_sub; j++) {
	  stage3 = stage3array_sub[j];
	  Stage3_free(&stage3);
	}
	FREE(stage3array_sub);

      } else {

	*chimera = Chimera_new(nonchimericbest,chimerapos,chimeraequivpos);
	    
	best0 = stage3array_sub[bestfrom];
	best1 = stage3array_sub[bestto];

	Chimera_find_exonexon(*chimera,best0,best1,Params_genome(params));
	chimerapos = Chimera_pos(*chimera);	
	chimeraequivpos = Chimera_equivpos(*chimera);

	for (oldj = npaths_sub-1, newj = npaths_sub-1; oldj >= 0; oldj--) {
	  if (oldj != bestfrom && oldj != bestto) {
	    /* fprintf(stderr,"Assigning old %d to new %d\n",oldj,newj); */
	    stage3array_sub[newj] = stage3array_sub[oldj];
	    newj--;
	  }
	}
	if ((stage3array_sub[0] = Stage3_copy_bounded(best0,0,chimeraequivpos+chimera_overlap)) == NULL) {
	  stage3array_sub[0] = best0;
	} else {
	  Stage3_free(&best0);
	}
	if ((stage3array_sub[1] = Stage3_copy_bounded(best1,chimerapos+1-chimera_overlap,querylength)) == NULL) {
	  stage3array_sub[1] = best1;
	} else {
	  Stage3_free(&best1);
	}
	
	for (j = 0; j < *npaths; j++) {
	  stage3 = stage3array[j];
	  Stage3_free(&stage3);
	}
	FREE(stage3array);
	
	stage3array = stage3array_sub;
	*npaths = npaths_sub;
      }
    }
  }
  return stage3array;
}
#endif


static Stage3_T *
apply_stage3 (Chimera_T *chimera, int *npaths, Stage1_T stage1, Sequence_T queryseq, Sequence_T usersegment, 
	      Params_T params, int maxpeelback, int sufflookback, int nsufflookback,
	      int badoligos, int nullgap, Oligoindex_T oligoindex, Pairpool_T pairpool, 
	      Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
	      char *gbuffer1, char *gbuffer2, char *gbuffer3, int gbufferlen) {
  Stage3_T *stage3array = NULL, nonchimericbest;
  List_T matchlist = NULL;
  bool testchimerap = false;
  int breakpoint = -1, margin;
  
  *npaths = 0;

  matchlist = Stage1_matchlist(stage1,Params_indexdb(params),Params_chromosome_iit(params),
			       Params_chrsubset(params),crossspeciesp);
  stage3array = stage3_from_matchlist(&(*npaths),stage3array,stage1,matchlist,queryseq,
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
      debug2(printf("existing alignment is too short\n"));
      testchimerap = false;

    } else if (Stage3_fracidentity(nonchimericbest) < CHIMERA_IDENTITY &&
	       (breakpoint = Chimera_detect(&margin,nonchimericbest,queryseq,CHIMERA_FVALUE)) >= 0 &&
	       margin >= chimera_margin) {
      debug2(printf("break in alignment quality at %d and margin %d\n",breakpoint,margin));
      testchimerap = true;

    } else if (Stage3_largemargin(&breakpoint,nonchimericbest,queryseq) >= chimera_margin) {
      debug2(printf("large margin of %d: breakpoint = %d\n",
		    Stage3_largemargin(&breakpoint,nonchimericbest,queryseq),breakpoint));
      testchimerap = true;

    } else {
      debug2(printf("good alignment with identity %f and margin %d\n",
		    Stage3_fracidentity(nonchimericbest),
		    Stage3_largemargin(&breakpoint,nonchimericbest,queryseq)));
      testchimerap = false;
#endif

    }

    if (testchimerap == false) {
      *chimera = (Chimera_T) NULL;
    } else {
#ifdef BETATEST
      stage3array = check_for_chimera(&(*chimera),&(*npaths),breakpoint,nonchimericbest,stage3array,
				      queryseq,usersegment,params,
				      maxpeelback,sufflookback,nsufflookback,badoligos,nullgap,
				      oligoindex,pairpool,dynprogL,dynprogM,dynprogR,
				      gbuffer1,gbuffer2,gbuffer3,gbufferlen);
#endif
    }

  }

  return stage3array;
}




#ifdef HAVE_PTHREAD
static void *
worker_thread (void *data)
#else
static void
single_thread (FILE *input, Params_T params, Sequence_T usersegment)
#endif
{
  Oligoindex_T oligoindex;
  int badoligos;
  Dynprog_T dynprogL, dynprogM, dynprogR;
  int breakpoint;
  char *gbuffer1 = NULL, *gbuffer2 = NULL, *gbuffer3 = NULL; /* Individual space for reading genome sequence,
								revcomp, and strain calculations */
  int nextchar;
  int trim_start, trim_end;
  int gbufferlen = 0,
    stuttercycles, stutterhits, 
    maxpeelback, sufflookback, nsufflookback, nullgap;
  int extramaterial_end, extramaterial_paired;
  List_T matchlist;
  Pairpool_T pairpool;
  Request_T request;
  Sequence_T queryseq;
  Chimera_T chimera;
  Result_T result;
  int jobid = 0;

  Stage1_T stage1;
  Stage3_T *stage3array;
  int npaths;

#ifdef HAVE_PTHREAD
  Reqpost_T reqpost = (Reqpost_T) data;
  Params_T params = Reqpost_params(reqpost);
  Sequence_T usersegment;
#endif

  debug(printf("worker_thread: Starting\n"));
  
  /* Thread-specific data and storage */
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
  gbufferlen = MAXSEQLEN + maxintronlen_bound + 2*Params_maxextension(params) + 
    1000;			/* Add extra 1000 just to be sure */
  gbuffer1 = (char *) CALLOC(gbufferlen+1,sizeof(char));
  gbuffer2 = (char *) CALLOC(gbufferlen+1,sizeof(char));
  gbuffer3 = (char *) CALLOC(gbufferlen+1,sizeof(char));
  pairpool = Pairpool_new();

#ifdef HAVE_PTHREAD
  while ((request = Reqpost_get_request(reqpost))) {
#else
  while ((queryseq = Sequence_read(&nextchar,input)) != NULL) {
    if (jobid == 0 && nextchar == '>' && batchp == false) {
      fprintf(stderr,"Note: batch mode (-B) is strongly recommended for large runs.\n");
    }
#endif

#ifdef HAVE_PTHREAD
#else
    request = Request_new(jobid++,queryseq,usersegment);
#endif
    jobid = Request_id(request);
    usersegment = Request_usersegment(request);
    debug(printf("worker_thread: Got request id\n",jobid));
    queryseq = Request_queryseq(request);
    Pairpool_reset(pairpool);

    debug1(printf("oligo depth: %f\n",Oligoindex_set_inquery(&badoligos,&trim_start,&trim_end,oligoindex,queryseq)));
    debug1(printf("badoligos: %d/%d\n",badoligos,Sequence_fulllength(queryseq)));
    debug1(printf("trimstart: %d\n",trim_start));
    debug1(printf("trimend: %d\n",trim_end));

    if (Sequence_fulllength(queryseq) <= 0) {
#ifdef HAVE_PTHREAD
      Reqpost_put_result(reqpost,Result_new(jobid,(Chimera_T) NULL,(Stage1_T) NULL,(Stage3_T *) NULL,0,EMPTY_SEQUENCE));
#else
      output_one(Result_new(jobid,(Chimera_T) NULL,(Stage1_T) NULL,(Stage3_T *) NULL,0,EMPTY_SEQUENCE),request,params);
#endif
      
    } else if (!maponlyp && ((Oligoindex_set_inquery(&badoligos,&trim_start,&trim_end,oligoindex,queryseq) > MAX_OLIGO_DEPTH &&
			      crossspeciesp == false &&
			      Sequence_fulllength(queryseq) < OLIGO_DEPTH_SEQLENGTH) ||
			     (double) badoligos/(double) Sequence_fulllength(queryseq) > MAX_BADOLIGOS)) {
      /* Sequence is repetitive or poor */
#ifdef HAVE_PTHREAD
      Reqpost_put_result(reqpost,Result_new(jobid,(Chimera_T) NULL,(Stage1_T) NULL,(Stage3_T *) NULL,0,REPETITIVE));
#else
      output_one(Result_new(jobid,(Chimera_T) NULL,(Stage1_T) NULL,(Stage3_T *) NULL,0,REPETITIVE),request,params);
#endif

    } else {

TRY

      if (usersegment != NULL && userstage1p == false) {
	stage3array = stage3_from_usersegment(&npaths,queryseq,usersegment,params,
					      maxpeelback,sufflookback,nsufflookback,
					      badoligos,nullgap,
					      oligoindex,pairpool,dynprogL,dynprogM,dynprogR,
					      gbuffer1,gbuffer2,gbuffer3,gbufferlen);
#ifdef HAVE_PTHREAD
	Reqpost_put_result(reqpost,Result_new(jobid,(Chimera_T) NULL,(Stage1_T) NULL,stage3array,npaths,NO_FAILURE));
#else
	output_one(Result_new(jobid,(Chimera_T) NULL,(Stage1_T) NULL,stage3array,npaths,NO_FAILURE),request,params);
#endif

      } else if (maponlyp == true) {
	/* Sequence_trim(queryseq,trim_start,trim_end); */
	stage1 = Stage1_compute(queryseq,Params_indexdb(params),Params_chromosome_iit(params),
				Params_chrsubset(params),maxintronlen_bound,stuttercycles,
				stutterhits,crossspeciesp);
	matchlist = Stage1_matchlist(stage1,Params_indexdb(params),Params_chromosome_iit(params),
				     Params_chrsubset(params),crossspeciesp);
	npaths = List_length(matchlist);
#ifdef HAVE_PTHREAD
	Reqpost_put_result(reqpost,Result_new(jobid,(Chimera_T) NULL,stage1,(Stage3_T *) NULL,npaths,NO_FAILURE));
#else
	output_one(Result_new(jobid,(Chimera_T) NULL,stage1,(Stage3_T *) NULL,npaths,NO_FAILURE),request,params);
#endif

      } else {	
	Sequence_trim(queryseq,trim_start,trim_end);
	stage1 = Stage1_compute(queryseq,Params_indexdb(params),Params_chromosome_iit(params),
				Params_chrsubset(params),maxintronlen_bound,stuttercycles,
				stutterhits,crossspeciesp);
	stage3array = apply_stage3(&chimera,&npaths,stage1,queryseq,usersegment,params,
				   maxpeelback,sufflookback,nsufflookback,badoligos,nullgap,
				   oligoindex,pairpool,dynprogL,dynprogM,dynprogR,gbuffer1,gbuffer2,
				   gbuffer3,gbufferlen);
	Stage1_free(&stage1);
#ifdef HAVE_PTHREAD
	Reqpost_put_result(reqpost,Result_new(jobid,chimera,(Stage1_T) NULL,stage3array,npaths,NO_FAILURE));
#else
	output_one(Result_new(jobid,chimera,(Stage1_T) NULL,stage3array,npaths,NO_FAILURE),request,params);
#endif
      }

ELSE
  if (Sequence_accession(queryseq) == NULL) {
    fprintf(stderr,"Problem with unnamed sequence (%d bp)\n",Sequence_fulllength(queryseq));
  } else {
    fprintf(stderr,"Problem with sequence %s (%d bp)\n",
	    Sequence_accession(queryseq),Sequence_fulllength(queryseq));
  }
  RERAISE;
END_TRY;

    }
  }

  Pairpool_free(&pairpool);
  FREE(gbuffer3);
  FREE(gbuffer2);
  FREE(gbuffer1);
  Dynprog_free(&dynprogR);
  Dynprog_free(&dynprogM);
  Dynprog_free(&dynprogL);
  Oligoindex_free(&oligoindex);
#ifdef HAVE_PTHREAD
  return NULL;
#else
  return;
#endif
}


static void
align_mutations (FILE *input, Sequence_T mutationrefseq, Params_T params) {
  Oligoindex_T oligoindex;
  int badoligos;
  Dynprog_T dynprogL, dynprogM, dynprogR;
  char *gbuffer1 = NULL, *gbuffer2 = NULL, *gbuffer3 = NULL; /* Individual space for reading genome sequence,
								revcomp, and strain calculations */
  int nextchar;
  int trim_start, trim_end;
  int gbufferlen = 0,
    stuttercycles, stutterhits, 
    maxpeelback, sufflookback, nsufflookback, nullgap;
  int extramaterial_end, extramaterial_paired;
  Pairpool_T pairpool;
  Genomicpos_T genomicstart, genomiclength;
  Sequence_T genomicseg, queryseq;
  int jobid = 0;

  Chimera_T chimera;
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
  gbufferlen = MAXSEQLEN + maxintronlen_bound + 2*Params_maxextension(params) + 
    1000;			/* Add extra 1000 just to be sure */
  gbuffer1 = (char *) CALLOC(gbufferlen+1,sizeof(char));
  gbuffer2 = (char *) CALLOC(gbufferlen+1,sizeof(char));
  gbuffer3 = (char *) CALLOC(gbufferlen+1,sizeof(char));
  pairpool = Pairpool_new();

  Pairpool_reset(pairpool);

  Oligoindex_set_inquery(&badoligos,&trim_start,&trim_end,oligoindex,mutationrefseq);
  Sequence_trim(mutationrefseq,trim_start,trim_end);
  stage1 = Stage1_compute(mutationrefseq,Params_indexdb(params),Params_chromosome_iit(params),
			  Params_chrsubset(params),maxintronlen_bound,stuttercycles,stutterhits,
			  crossspeciesp);
  stage3array = apply_stage3(&chimera,&npaths,stage1,mutationrefseq,/*usersegment*/NULL,params,
			     maxpeelback,sufflookback,nsufflookback,badoligos,nullgap,
			     oligoindex,pairpool,dynprogL,dynprogM,dynprogR,gbuffer1,gbuffer2,
			     gbuffer3,gbufferlen);
  /* chimera should be NULL */
  for (i = 1; i < npaths; i++) {
    stage3 = stage3array[i];
    Stage3_free(&stage3);
  }
  if (npaths > 0) {
    stage3ref = stage3array[0];
    stage3ref = Stage3_translate_genomic(stage3ref,fulllengthp,/*truncatep*/false);
    FREE(stage3array);

    Stage3_genomicbounds(&genomicstart,&genomiclength,stage3ref);
    genomicseg = Genome_get_segment(Params_genome(params),genomicstart,genomiclength,/*revcomp*/false,
				    gbuffer1,gbuffer2,gbufferlen);

    while ((queryseq = Sequence_read(&nextchar,input)) != NULL) {
      if (jobid == 0 && nextchar == '>' && batchp == false) {
        fprintf(stderr,"Note: batch mode (-B) is strongly recommended for large runs.\n");
      }

      Pairpool_reset(pairpool);

      Sequence_print_header(queryseq,checksump);
      if (Sequence_fulllength(queryseq) <= 0) {
	printf("Paths (0): *** Empty sequence ***\n");

      } else if ((Oligoindex_set_inquery(&badoligos,&trim_start,&trim_end,oligoindex,queryseq) > MAX_OLIGO_DEPTH &&
		  crossspeciesp == false &&
		  Sequence_fulllength(queryseq) < OLIGO_DEPTH_SEQLENGTH) ||
		 (double) badoligos/(double) Sequence_fulllength(queryseq) > MAX_BADOLIGOS) {
	printf("Paths (0): *** Repetitive sequence ***\n");

      } else {
	stage3array = stage3_from_usersegment(&npaths,queryseq,genomicseg,params,
					      maxpeelback,sufflookback,nsufflookback,
					      badoligos,nullgap,
					      oligoindex,pairpool,dynprogL,dynprogM,dynprogR,
					      gbuffer1,gbuffer2,gbuffer3,gbufferlen);
	if (npaths == 0) {
	  printf("Paths (0):\n");
	} else {
	  printf("Paths (1):\n");
	  Stage3_translate_cdna(stage3array[0],stage3ref,literalrefp);
	  Stage3_fix_cdna_direction(stage3array[0],stage3ref);
	  Stage3_print_mutations(stage3array[0],stage3ref,Params_chromosome_iit(params),
				 dbversion,showalignp,/*zerobasedp*/false,continuousp,diagnosticp,proteinmode,
				 invertmode,nointronlenp,wraplength);
	  for (i = 0; i < npaths; i++) {
	    stage3 = stage3array[i];
	    Stage3_free(&stage3);
	  }
	  FREE(stage3array);
	}
      }
      Sequence_free(&queryseq);
      jobid++;
    }
    Sequence_free(&genomicseg);
    Stage3_free(&stage3ref);
  }

  Stage1_free(&stage1);

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

/*
static int
calculate_maxlookback (double failureprob, double pctidentity) {
  return rint(log(failureprob)/log(1.0 - pow(pctidentity,8.0)) + INDEXSIZE - 1);
}
*/



int
main (int argc, char *argv[]) {
  IIT_T chromosome_iit = NULL, contig_iit = NULL, altstrain_iit = NULL, map_iit = NULL;
  Chrsubset_T chrsubset = NULL;
  Indexdb_T indexdb = NULL;
  Genome_T genome = NULL;
  Sequence_T usersegment = NULL, mutationrefseq = NULL;
  char *genomesubdir = NULL, *mapdir = NULL, *iitfile = NULL, *olddbroot, *fileroot = NULL;
  int ret, i;
  FILE *input, *mutation_fp = NULL;
  Params_T params;
  Blackboard_T blackboard;

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
    extraband_paired;
  /* int lookback_safety = 2*INDEXSIZE; */
  bool showcontigp = true;
  double failureprob = 0.001, pctidentity = 0.95;

#ifdef HAVE_PTHREAD
  pthread_t input_thread_id, output_thread_id, worker_thread_id;
  pthread_attr_t thread_attr_detach, thread_attr_join;
#endif

  int opt;
  extern int optind;
  extern char *optarg;
  int long_option_index = 0;

  while ((opt = getopt_long(argc,argv,"D:d:R:g:C:BL:x:1q:Xw:t:sec:USA39n:ZO5M:m:bo:EPQFTNI:i:l:V?",
			    long_options, &long_option_index)) != -1) {
    switch (opt) {
    case 'D': user_genomedir = optarg; break;
    case 'd': 
      dbroot = (char *) CALLOC(strlen(optarg)+1,sizeof(char));
      strcpy(dbroot,optarg);
      break;
    case 'R': releasestring = optarg; break;
    case 'g': user_genomicseg = optarg; break;
    case 'C': user_chrsubsetfile = optarg; break;

    case 'B': batchp = true; break;
    case 'L': maxintronlen_bound = atoi(optarg); break;
    case 'x': chimera_margin = atoi(optarg); break;
    case '1': maponlyp = true; break;
    case 'q': 
      index1interval = atoi(optarg);
      if (INDEX1PART % index1interval != 0) {
	fprintf(stderr,"Selected indexing interval %d is not evenly divisible into the size of the oligomer %d\n",
		index1interval,INDEX1PART);
	exit(9);
      }
      break;
    case 'X': crossspeciesp = true; break;
    case 'w': mutationref = optarg; break;
    case 't': nworkers = atoi(optarg); break;
    case 's': altstrainp = true; break;
    case 'e': extend_mismatch_p = true; break;
    case 'c': user_chrsubsetname = optarg; break;
    case 'U': end_microexons_p = true; break;

    case 'S': summaryonlyp = true; break;
    case 'A': showalignp = true; break;
    case '3': continuousp = true; break;
    case '9': diagnosticp = true; break;
    case 'n': maxpaths = atoi(optarg); break;
    case 'Z': compressoutp = true; break;
    case 'O': orderedp = true; break;
    case '5': checksump = true; break;
    case 'M': user_mapdir = optarg; break;
    case 'm': map_iitfile = optarg; break;
    case 'b': map_bothstrands_p = true; break;
    case 'o': chimera_overlap = atoi(optarg); break;

    case 'E': cdna_exons_p = true; break;
    case 'P': protein_cdna_p = true; break;
    case 'Q': protein_genomic_p = true; break;
    case 'F': fulllengthp = true; break;
    case 'T': truncatep = true; fulllengthp = true; break;
    case 'N': nointronlenp = true; break;
    case 'I': invertmode = atoi(optarg); break;
    case 'i': ngap = atoi(optarg); break;
    case 'l': wraplength = atoi(optarg); break;

    case 'V': print_program_version(); exit(0);
    case '?': print_program_usage(); exit(0);
    default: exit(9);
    }
  }
  argc -= optind;
  argv += optind;
      
  if (crossspeciesp == false) {
    indexsize = 8;
    maxpeelback = 11;		/* Needs to be at least indexsize
				   because stage 2 jumps by indexsize.
				   Also should exceed length of
				   repeated nucleotides (e.g., a
				   string of consecutive T's) */

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
    indexsize = 6;
    maxpeelback = 6;		/* Needs to be at least indexsize because stage 2 jumps by indexsize */
    if (sufflookback == 0) {
      sufflookback = 120;
    }
    nullgap = 1200;
    extramaterial_end = 20;
    extramaterial_paired = 16;
    extraband_single = 6;
    extraband_end = 7;
    extraband_paired = 3;
  }

  /*
  maxlookback = calculate_maxlookback(failureprob,pctidentity);
  fprintf(stderr,"Max lookback = %d + %d for failure prob = %.3e and pctidentity = %.3f\n",
	  maxlookback,lookback_safety,failureprob,pctidentity);
  maxlookback += lookback_safety;
  */

  /* Complement_init(); */
  Oligoindex_init(indexsize);
  Dynprog_init(nullgap,EXTRAQUERYGAP,maxpeelback,extramaterial_end,extramaterial_paired);
  if (user_genomicseg != NULL) {
    /* Map against user-provided genomic segment */
    showcontigp = false;

    if ((input = fopen(user_genomicseg,"r")) == NULL) {
      fprintf(stderr,"Can't open file %s\n",user_genomicseg);
      exit(9);
    }
    if ((usersegment = Sequence_read_unlimited(input)) == NULL) {
      fprintf(stderr,"File %s is empty\n",user_genomicseg);
      exit(9);
    }
    fclose(input);

    genome = (Genome_T) NULL;
    dbversion = (char *) NULL;
    if (userstage1p == true) {
      indexdb = Indexdb_new_segment(Sequence_fullpointer(usersegment),index1interval);
    }

    if (argc == 0) {
      input = stdin;
    } else if ((input = fopen(argv[0],"r")) == NULL) {
      fprintf(stderr,"Can't open file %s\n",argv[0]);
      exit(9);
    }

  } else if (dbroot == NULL) {
    fprintf(stderr,"Need to specify the -d or the -g flag\n");
    print_program_usage();
    exit(9);

  } else {
    /* Map against genome */
    if (releasestring != NULL) {
      olddbroot = dbroot;
      dbroot = (char *) CALLOC(strlen(dbroot)+strlen("_R")+strlen(releasestring)+1,sizeof(char));
      sprintf(dbroot,"%s_R%s",olddbroot,releasestring);
      FREE(olddbroot);
    }
    genomesubdir = Datadir_find_genomesubdir(&fileroot,&dbversion,user_genomedir,dbroot);

    if (argc == 0) {
      input = stdin;
    } else if ((input = fopen(argv[0],"r")) == NULL) {
      fprintf(stderr,"Can't open file %s\n",argv[0]);
      exit(9);
    }
    
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
	fprintf(stderr,"Map file %s.iit not found in %s\n",map_iitfile,mapdir);
	fprintf(stderr,"Either install that file or specify a full directory path\n");
	fprintf(stderr,"using the -M flag to gmap.\n");
	exit(9);
      }
      FREE(iitfile);
      FREE(mapdir);
    }

    indexdb = Indexdb_new_genome(genomesubdir,fileroot,batchp);
    genome = Genome_new(genomesubdir,fileroot,uncompressedp,batchp);

    if (altstrainp == true) {
      iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
				strlen(fileroot)+strlen(".altstrain.iit")+1,sizeof(char));
      sprintf(iitfile,"%s/%s.altstrain.iit",genomesubdir,fileroot);
      global_altstrain_iit = altstrain_iit = IIT_read(iitfile,NULL,true);
      FREE(iitfile);
    }

    FREE(genomesubdir);
    FREE(fileroot);
    FREE(dbroot);
  }

  params = Params_new(genome,altstrain_iit,indexdb,chromosome_iit,chrsubset,
		      contig_iit,map_iit,maxextension,stuttercycles,stutterhits,
		      indexsize,maxpeelback,sufflookback,nsufflookback,nullgap,
		      extramaterial_end,extramaterial_paired,
		      extraband_single,extraband_end,extraband_paired);

  if (mutationref != NULL) {
    /* Mutation analysis */
    if ((mutation_fp = fopen(mutationref,"r")) == NULL) {
      fprintf(stderr,"Can't open file %s\n",mutationref);
      exit(9);
    } else {
      mutationrefseq = Sequence_read_unlimited(mutation_fp);
      chimera_margin = -1;
      align_mutations(input,mutationrefseq,params);
      fclose(mutation_fp);
      Sequence_free(&mutationrefseq);
      Params_free(&params);
    }

  } else {

    /* Make blackboard and threads */
#ifdef HAVE_PTHREAD
    blackboard = Blackboard_new(input,usersegment,nworkers,params);
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
    Params_free(&params);
    Blackboard_free(&blackboard);

#else  /* HAVE_PTHREAD */

    single_thread(input,params,usersegment);
    Params_free(&params);

#endif /* HAVE_PTHREAD */
  }
  
  Dynprog_term();

  if (argc > 0) {
    fclose(input);
  }
  if (indexdb != NULL) {
    Indexdb_free(&indexdb);
  }
  if (dbversion != NULL) {
    FREE(dbversion);
  }
  if (altstrain_iit != NULL) {
    IIT_free_mmapped(&altstrain_iit);
  }
  if (genome != NULL) {
    Genome_free(&genome);
  }
  if (map_iit != NULL) {
    IIT_free_mmapped(&map_iit);
  }
  if (contig_iit != NULL) {
    IIT_free_mmapped(&contig_iit);
  }
  if (chromosome_iit != NULL) {
    IIT_free_mmapped(&chromosome_iit);
  }
  if (chrsubset != NULL) {
    Chrsubset_free(&chrsubset);
  }

  return 0;
}
