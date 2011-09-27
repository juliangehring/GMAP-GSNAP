static char rcsid[] = "$Id: gmap.c,v 1.235 2005/03/11 17:56:43 twu Exp $";
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

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
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
static double chimera_threshold = 0.00;
static bool maponlyp = false;
static bool userstage1p = true;	/* Apply stage 1 for user-provided genomic segments */
static bool crossspeciesp = false;
static char *mutationref = NULL;
static bool literalrefp = false;
static bool aitstrainp = false;
static bool nworkers = 1;	/* (int) sysconf(_SC_NPROCESSORS_ONLN) */
static bool extend_mismatch_p = false;
static char *user_chrsubsetname = NULL;


/* Output options */
static bool summaryonlyp = false;
static bool showalignp = false;
static bool continuousp = false;
static bool diagnosticp = false;
static int maxpaths = 5;
static bool compressoutp = false;
static bool orderedp = false;
static bool checksump = false;
static char *user_mapdir = NULL;
static char *map_iitfile = NULL;
static bool map_bothstrands_p = false;


/* Alignment options */
static bool cdna_exons_p = false;
static bool protein_cdna_p = false;
static bool protein_genomic_p = false;
static bool fulllengthp = false;
static int proteinmode = 1;
static bool uncompressedp = false;
static bool nointronlenp = false;
static int invertmode = 0;
static int ngap = 3;
static int wraplength = 50;


static struct option long_options[] = {
  /* Input options */
  {"dir", required_argument, 0, 'D'},	/* user_genomedir */
  {"db", required_argument, 0, 'd'}, /* dbroot */
  {"rel", required_argument, 0, 'R'}, /* releasestring */
  {"gseg", required_argument, 0, 'g'}, /* user_genomicseg */

  /* Compute options */
  {"batch", no_argument, 0, 'B'}, /* batchp */
  {"length", required_argument, 0, 'L'}, /* maxintronlen_bound */
  {"chimera", required_argument, 0, 'x'}, /* chimera_threshold */
  {"maponly", no_argument, 0, '1'}, /* maponlyp */
  {"cross", no_argument, 0, 'X'}, /* crossspeciesp */
  {"mutationref", required_argument, 0, 'w'},	/* mutationref */
  {"nthreads", required_argument, 0, 't'}, /* nworkers */
  {"altstrain", no_argument, 0, 's'},	/* aitstrainp */
  {"extend", no_argument, 0, 'e'}, /* extend_mismatch_p */
  {"chrsubsetfile", required_argument, 0, 'C'}, /* user_chrsubsetfile */
  {"chrsubset", required_argument, 0, 'c'}, /* user_chrsubsetname */

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

  /* Alignment options */
  {"exons", no_argument, 0, 'E'}, /* cdna_exons_p */
  {"protein_dna", no_argument, 0, 'P'}, /* protein_cdna_p */
  {"protein_gen", no_argument, 0, 'Q'}, /* protein_genomic_p */
  {"fulllength", no_argument, 0, 'F'}, /* fulllengthp */
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
  -x, --chimera=DOUBLE           Threshold for chimera search (default 0.00,\n\
                                 turned off).  Turn on with 0.6 or higher.\n\
  -1, --maponly                  Report stage 1 only (genomic mapping only)\n\
  -X, --cross                    Cross-species mode (changes parameters)\n\
  -w, --mutationref=filename     Mutation reference\n\
  -t, --nthreads=INT             Number of worker threads\n\
  -s, --altstrain                Search alternate strains in addition\n\
  -e, --extend                   Extend alignment to show mismatch at end\n\
  -C, --chrsubsetfile=filename   User-suppled chromosome subset file\n\
  -c, --chrsubset=string         Chromosome subset to search\n\
\n\
Output options\n\
  -S, --summary                  Show summary of alignments only\n\
  -A, --align                    Show alignments\n\
  -3, --continuous               Show alignment in three continuous lines\n\
  -9, --diagnostic               Show diagnostics (alignment with *'s for\n\
                                 stage 3)\n\
  -n, --npaths=INT               Maximum number of paths to show\n\
  -Z, --compress                 Print output in compressed format\n\
  -O, --ordered                  Print output in same order as input (relevant\n\
                                 only if there is more than worker thread)\n\
  -5, --md5                      Print MD5 checksum for each query sequence\n\
  -M, --mapdir=directory         Map directory\n\
  -m, --map=iitfile              Map file\n\
  -b, --mapboth                  Map both strands of genome\n\
\n\
Alignment output options\n\
  -E, --exons                    Print cDNA exons\n\
  -P, --protein_dna              Print protein sequence (cDNA)\n\
  -Q, --protein_gen              Print protein sequence (genomic)\n\
  -F, --fulllength               Assume full-length protein, starting with Met\n\
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


void
print_path_summaries (Stage3_T *stage3array, int npaths, IIT_T chromosome_iit, IIT_T contig_iit, 
		      char *dbversion, int maxpaths, bool chimerap, bool zerobasedp,
		      int ntrimmed) {
  Stage3_T stage3;
  int pathnum;

  if (npaths == 0) {
    printf("Paths (0):\n\n");
  } else {
    printf("Paths (%d):",npaths);
    if (chimerap == true) {
      printf(" *** Possible chimera ***");
    }
    printf("\n");
    for (pathnum = 1; pathnum <= npaths && pathnum <= maxpaths; pathnum++) {
      stage3 = stage3array[pathnum-1];
      Stage3_print_pathsummary(stage3,pathnum,chromosome_iit,contig_iit,
			       dbversion,zerobasedp,ntrimmed,fulllengthp);
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

  debug(printf("input_thread: Starting\n"));
  usersegment = Blackboard_usersegment(blackboard);
  while ((queryseq = Sequence_read(input,/*polya_trim*/true)) != NULL) {
    debug(printf("input_thread: Putting request id %d\n",inputid));
    request = Request_new(inputid++,queryseq,usersegment);
    Blackboard_put_request(blackboard,request);
    if (inputid == 5 && batchp == false) {
      fprintf(stderr,"Note: batch mode (-B) is strongly recommended for large runs.\n");
    }
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
  bool chimerap;

  queryseq = Request_queryseq(request);
  if (compressoutp == false && cdna_exons_p == false &&
      protein_cdna_p == false && protein_genomic_p == false) {
    Sequence_print_header(queryseq,checksump);
  }

  if (maponlyp == true) {
    stage1 = Result_stage1(result);
    matchlist = Stage1_matchlist(stage1,Params_indexdb(params),Params_chromosome_iit(params),Params_chrsubset(params));
    npaths = List_length(matchlist);
    if (npaths == 0) {
      printf("Paths (0):\n\n");
    } else {
      print_maponly(matchlist,npaths,Params_chromosome_iit(params),Params_contig_iit(params),
		    dbversion,maxpaths,/*zerobasedp*/false);
    }
  } else if (compressoutp == true) {
    stage3array = Result_array(&npaths,result);
    if (Result_chimerapos(result) >= 0) {
      chimerap = true;
    } else {
      chimerap = false;
    }
    for (pathnum = 0; pathnum < npaths && pathnum < maxpaths; pathnum++) {
      Stage3_print_compressed(stage3array[pathnum],queryseq,Params_chromosome_iit(params),
			      dbversion,pathnum+1,npaths,checksump,chimerap,
			      /*zerobasedp*/false);
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
    for (pathnum = 0; pathnum < npaths && pathnum < maxpaths; pathnum++) {
      Sequence_print_header(queryseq,checksump);
      Stage3_print_protein_cdna(stage3array[pathnum],wraplength,fulllengthp);
    }

  } else if (protein_genomic_p == true) {
    stage3array = Result_array(&npaths,result);
    for (pathnum = 0; pathnum < npaths && pathnum < maxpaths; pathnum++) {
      Sequence_print_header(queryseq,checksump);
      Stage3_print_protein_genomic(stage3array[pathnum],wraplength,fulllengthp);
    }

  } else {
    /* Usual output */
    stage3array = Result_array(&npaths,result);
    if (npaths == 0) {
      printf("Paths (0):\n\n");
    } else {
      if (Result_chimerapos(result) >= 0) {
	chimerap = true;
      } else {
	chimerap = false;
      }
      print_path_summaries(stage3array,npaths,Params_chromosome_iit(params),Params_contig_iit(params),
			   dbversion,maxpaths,chimerap,
			   /*zerobasedp*/false,Sequence_ntrimmed(queryseq));
    
      if (Params_map_iit(params) != NULL) {
	printf("Maps:\n");
	for (pathnum = 0; pathnum < npaths && pathnum < maxpaths; pathnum++) {
	  stage3 = stage3array[pathnum];
	  Stage3_print_map(stage3,Params_map_iit(params),pathnum+1,map_bothstrands_p);
	}
      }

      if (summaryonlyp == true || showalignp == true) {
	printf("Alignments:\n");
	for (pathnum = 0; pathnum < npaths && pathnum < maxpaths; pathnum++) {
	  printf("  Alignment for path %d:\n\n",pathnum+1);
	  stage3 = stage3array[pathnum];
	  Stage3_print_alignment(stage3,Params_chromosome_iit(params),summaryonlyp,/*universalp*/false,
				 /*zerobasedp*/false,continuousp,diagnosticp,
				 /*flipgenomep*/true,proteinmode,invertmode,nointronlenp,wraplength);
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
stage3array_from_list (int *npaths, List_T stage3list, int querylength) {
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
		   Matchpairend_T matchpairend, int straintype, char *strain, 
		   Oligoindex_T oligoindex, int indexsize, Pairpool_T pairpool, 
		   int sufflookback, int nsufflookback, Chrnum_T chrnum,
		   Genomicpos_T chroffset, Genomicpos_T chrpos, bool watsonp, int maxpeelback, int nullgap,
		   int extramaterial_end, int extramaterial_paired, int extraband_single, int extraband_end, 
		   int extraband_paired, int ngap, Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
		   IIT_T altstrain_iit) {
  Stage2_T stage2;
  Stage3_T stage3;

  if ((stage2 = Stage2_compute(queryseq,genomicseg,oligoindex,indexsize,
			       pairpool,sufflookback,nsufflookback)) != NULL) {
    if ((stage3 = Stage3_compute(stage2,queryseq,genomicseg,matchpairend,straintype,strain,chrnum,chroffset,
				 chrpos,watsonp,maxpeelback,nullgap,extramaterial_end,extramaterial_paired,
				 extraband_single,extraband_end,extraband_paired,
				 ngap,extend_mismatch_p,pairpool,dynprogL,dynprogM,dynprogR,
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
			 int maxpeelback, int sufflookback, int nsufflookback, int nullgap, 
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
  genomiclength = chrlength = Sequence_length(usersegment);
  watsonp = true;
  strain = NULL;

  stage3list = update_stage3list(stage3list,queryseq,usersegment,MIXED,
				 0,NULL,oligoindex,indexsize,pairpool,
				 sufflookback,nsufflookback,chrnum,chroffset,chrpos,
				 watsonp,maxpeelback,nullgap,extramaterial_end,extramaterial_paired,
				 extraband_single,extraband_end,extraband_paired,
				 ngap,dynprogL,dynprogM,dynprogR,altstrain_iit);
	
  revcomp = Sequence_revcomp(usersegment);
  watsonp = false;
  stage3list = update_stage3list(stage3list,queryseq,revcomp,MIXED,
				 0,NULL,oligoindex,indexsize,pairpool,
				 sufflookback,nsufflookback,chrnum,chroffset,chrpos,
				 watsonp,maxpeelback,nullgap,extramaterial_end,extramaterial_paired,
				 extraband_single,extraband_end,extraband_paired,
				 ngap,dynprogL,dynprogM,dynprogR,altstrain_iit);
  Sequence_free(&revcomp);

  if (stage3list == NULL) {
    *npaths = 0;
    return NULL;
  } else {
    return stage3array_from_list(&(*npaths),stage3list,Sequence_length(queryseq));
  }
}


static Stage3_T *
stage3_from_matchlist (int *npaths, Stage3_T *oldstage3array, Stage1_T stage1,
		       List_T matchlist, Sequence_T queryseq, Sequence_T usersegment, Params_T params, 
		       int maxpeelback, int sufflookback, int nsufflookback, int nullgap, 
		       Oligoindex_T oligoindex, Pairpool_T pairpool, 
		       Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
		       char *gbuffer1, char *gbuffer2, char *gbuffer3, int gbufferlen) {
  List_T stage3list = NULL, p;
  Intlist_T indices;
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
      chrlength = Sequence_length(usersegment);
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

    stage3list = update_stage3list(stage3list,queryseq,genomicseg,Matchpair_matchpairend(matchpair),
				   0,strain,oligoindex,indexsize,pairpool,
				   sufflookback,nsufflookback,chrnum,chroffset,chrpos,
				   watsonp,maxpeelback,nullgap,extramaterial_end,extramaterial_paired,
				   extraband_single,extraband_end,extraband_paired,
				   ngap,dynprogL,dynprogM,dynprogR,altstrain_iit);
    Sequence_free(&genomicseg);

    /* We rely upon the fact that gbuffer1 still holds the genomic segment.  This code is duplicated in get-genome.c */
    if (altstrain_iit != NULL) {
      indices = IIT_get(altstrain_iit,genomicstart+1,genomicstart+genomiclength-1);
      indexarray = Intlist_to_array(&nindices,indices);
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
	  stage3list = update_stage3list(stage3list,queryseq,genomicseg,Matchpair_matchpairend(matchpair),
					 type,strain,oligoindex,indexsize,pairpool,
					 sufflookback,nsufflookback,chrnum,chroffset,chrpos,
					 watsonp,maxpeelback,nullgap,extramaterial_end,extramaterial_paired,
					 extraband_single,extraband_end,extraband_paired,
					 ngap,dynprogL,dynprogM,dynprogR,altstrain_iit);
	  Sequence_free(&genomicseg);
	}
	FREE(indexarray);
      }
      Intlist_free(&indices);
    }
  }
	
  if (stage3list == NULL) {
    *npaths = 0;
    return NULL;
  } else {
    return stage3array_from_list(&(*npaths),stage3list,Sequence_length(queryseq));
  }
}


static Stage3_T *
apply_stage3 (int *chimerapos, int *npaths, Stage1_T stage1, Sequence_T queryseq, Sequence_T usersegment, 
	      Params_T params, int maxpeelback, int sufflookback, int nsufflookback,
	      int nullgap, Oligoindex_T oligoindex, Pairpool_T pairpool, 
	      Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
	      char *gbuffer1, char *gbuffer2, char *gbuffer3, int gbufferlen) {
  Stage3_T *stage3array = NULL, *stage3subarray, best0, best1, tmp1, tmp2;
  Stage1_T stage1sub;
  Matchpair_T match;
  List_T matchlist = NULL;
  Sequence_T querysubseq;
  Genomicpos_T genomicstart, genomiclength;
  int substart, subend;
  int newj, oldj, querylength;
  int stage1size;
  int bestfrom, bestto, bestpos, fromscore_3, toscore_5;

  *npaths = 0;

  matchlist = Stage1_matchlist(stage1,Params_indexdb(params),Params_chromosome_iit(params),Params_chrsubset(params));
  stage3array = stage3_from_matchlist(&(*npaths),stage3array,stage1,matchlist,queryseq,
				      usersegment,params,maxpeelback,sufflookback,nsufflookback,nullgap,
				      oligoindex,pairpool,dynprogL,dynprogM,dynprogR,gbuffer1,gbuffer2,
				      gbuffer3,gbufferlen);

  *chimerapos = -1;
  if (stage3array != NULL && Stage3_coverage(stage3array[0]) < chimera_threshold) {
    /* Possible chimera */
    if (Stage3_five_end(&substart,&subend,stage3array[0],queryseq) == true) {
      querysubseq = Sequence_subsequence(queryseq,substart,subend,false);
    } else {
      querysubseq = Sequence_subsequence(queryseq,substart,subend,true);
    }
    if (querysubseq != NULL) {
      stage1sub = Stage1_compute(querysubseq,Params_indexdb(params),Params_chromosome_iit(params),
				 Params_chrsubset(params),maxintronlen_bound,Params_stuttercycles(params),
				 Params_stutterhits(params),crossspeciesp);
      matchlist = Stage1_matchlist(stage1sub,Params_indexdb(params),Params_chromosome_iit(params),Params_chrsubset(params));
      stage3array = stage3_from_matchlist(&(*npaths),stage3array,stage1sub,matchlist,queryseq,
					  usersegment,params,maxpeelback,sufflookback,nsufflookback,nullgap,
					  oligoindex,pairpool,dynprogL,dynprogM,dynprogR,gbuffer1,gbuffer2,
					  gbuffer3,gbufferlen);

      if (*npaths > 1) {
	querylength = Sequence_length(queryseq);
	Chimera_bestpath(&bestfrom,&bestto,&bestpos,&fromscore_3,&toscore_5,stage3array,*npaths,querylength);

	if ((double) bestpos/(double) querylength < chimera_threshold ||
	    (double) (querylength - bestpos)/(double) querylength < chimera_threshold) {
	  *chimerapos = bestpos;
	  best0 = stage3array[bestfrom];
	  best1 = stage3array[bestto];
	  for (oldj = *npaths-1, newj = *npaths-1; oldj >= 0; oldj--) {
	    if (oldj != bestfrom && oldj != bestto) {
	      /* fprintf(stderr,"Assigning old %d to new %d\n",oldj,newj); */
	      stage3array[newj] = stage3array[oldj];
	      newj--;
	    }
	  }
	  stage3array[1] = best1;
	  stage3array[0] = best0;
	}
      }
      
      Stage1_free(&stage1sub);
      Sequence_free(&querysubseq);
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
  Dynprog_T dynprogL, dynprogM, dynprogR;
  int chimerapos;
  char *gbuffer1 = NULL, *gbuffer2 = NULL, *gbuffer3 = NULL; /* Individual space for reading genome sequence,
								revcomp, and strain calculations */
  int gbufferlen = 0, stage1size,
    stuttercycles, stutterhits, 
    maxpeelback, sufflookback, nsufflookback, nullgap;
  int extramaterial_end, extramaterial_paired;
  List_T matchlist;
  Pairpool_T pairpool;
  Request_T request;
  Sequence_T queryseq;
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
  while ((request = Reqpost_get_request(reqpost)))
#else
  while ((queryseq = Sequence_read(input,/*polya_trim*/true)) != NULL)
#endif
  {

#ifdef HAVE_PTHREAD
#else
    request = Request_new(jobid++,queryseq,usersegment);
#endif
    jobid = Request_id(request);
    usersegment = Request_usersegment(request);
    debug(printf("worker_thread: Got request id\n",jobid));
    queryseq = Request_queryseq(request);
    Pairpool_reset(pairpool);

    if (Sequence_length(queryseq) <= 0) {
#ifdef HAVE_PTHREAD
      Reqpost_put_result(reqpost,Result_new(jobid,-1,NULL,NULL,0));
#else
      output_one(Result_new(jobid,-1,NULL,NULL,0),request,params);
#endif
      
    } else {
      if (usersegment != NULL && userstage1p == false) {
	stage3array = stage3_from_usersegment(&npaths,queryseq,usersegment,params,
					      maxpeelback,sufflookback,nsufflookback,nullgap,
					      oligoindex,pairpool,dynprogL,dynprogM,dynprogR,
					      gbuffer1,gbuffer2,gbuffer3,gbufferlen);
#ifdef HAVE_PTHREAD
	Reqpost_put_result(reqpost,Result_new(jobid,-1,NULL,stage3array,npaths));
#else
	output_one(Result_new(jobid,-1,NULL,stage3array,npaths),request,params);
#endif

      } else if (maponlyp == true) {
	stage1 = Stage1_compute(queryseq,Params_indexdb(params),Params_chromosome_iit(params),
				Params_chrsubset(params),maxintronlen_bound,stuttercycles,
				stutterhits,crossspeciesp);
	matchlist = Stage1_matchlist(stage1,Params_indexdb(params),Params_chromosome_iit(params),
				     Params_chrsubset(params));
	npaths = List_length(matchlist);
#ifdef HAVE_PTHREAD
	Reqpost_put_result(reqpost,Result_new(jobid,-1,stage1,NULL,npaths));
#else
	output_one(Result_new(jobid,-1,stage1,NULL,npaths),request,params);
#endif

      } else {	
	stage1 = Stage1_compute(queryseq,Params_indexdb(params),Params_chromosome_iit(params),
				Params_chrsubset(params),maxintronlen_bound,stuttercycles,
				stutterhits,crossspeciesp);
	stage3array = apply_stage3(&chimerapos,&npaths,stage1,queryseq,usersegment,params,
				   maxpeelback,sufflookback,nsufflookback,nullgap,
				   oligoindex,pairpool,dynprogL,dynprogM,dynprogR,gbuffer1,gbuffer2,
				   gbuffer3,gbufferlen);
	Stage1_free(&stage1);
#ifdef HAVE_PTHREAD
	Reqpost_put_result(reqpost,Result_new(jobid,chimerapos,NULL,stage3array,npaths));
#else
	output_one(Result_new(jobid,chimerapos,NULL,stage3array,npaths),request,params);
#endif
      }
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
  Dynprog_T dynprogL, dynprogM, dynprogR;
  int chimerapos;
  char *gbuffer1 = NULL, *gbuffer2 = NULL, *gbuffer3 = NULL; /* Individual space for reading genome sequence,
								revcomp, and strain calculations */
  int gbufferlen = 0, stage1size, 
    stuttercycles, stutterhits, 
    maxpeelback, sufflookback, nsufflookback, nullgap;
  int extramaterial_end, extramaterial_paired;
  Pairpool_T pairpool;
  Genomicpos_T genomicstart, genomiclength;
  Sequence_T genomicseg, queryseq;
  int jobid = 0;

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

  stage1 = Stage1_compute(mutationrefseq,Params_indexdb(params),Params_chromosome_iit(params),
			  Params_chrsubset(params),maxintronlen_bound,stuttercycles,stutterhits,
			  crossspeciesp);
  stage3array = apply_stage3(&chimerapos,&npaths,stage1,mutationrefseq,/*usersegment*/NULL,params,
			     maxpeelback,sufflookback,nsufflookback,nullgap,
			     oligoindex,pairpool,dynprogL,dynprogM,dynprogR,gbuffer1,gbuffer2,
			     gbuffer3,gbufferlen);
  for (i = 1; i < npaths; i++) {
    stage3 = stage3array[i];
    Stage3_free(&stage3);
  }
  if (npaths > 0) {
    stage3ref = stage3array[0];
    Stage3_translate_genomic(stage3ref,fulllengthp);
    FREE(stage3array);

    Stage3_genomicbounds(&genomicstart,&genomiclength,stage3ref);
    genomicseg = Genome_get_segment(Params_genome(params),genomicstart,genomiclength,/*revcomp*/false,
				    gbuffer1,gbuffer2,gbufferlen);

    while ((queryseq = Sequence_read(input,/*polya_trim*/true)) != NULL) {
      Pairpool_reset(pairpool);

      Sequence_print_header(queryseq,checksump);
      if (Sequence_length(queryseq) <= 0) {
	printf("Paths (0):\n");
      } else {
	stage3array = stage3_from_usersegment(&npaths,queryseq,genomicseg,params,
					      maxpeelback,sufflookback,nsufflookback,nullgap,
					      oligoindex,pairpool,dynprogL,dynprogM,dynprogR,
					      gbuffer1,gbuffer2,gbuffer3,gbufferlen);
	if (npaths == 0) {
	  printf("Paths (0):\n");
	} else {
	  printf("Paths (1):\n");
	  Stage3_translate_cdna(stage3array[0],stage3ref,literalrefp);
	  Stage3_fix_cdna_direction(stage3array[0],stage3ref);
	  Stage3_print_mutations(stage3array[0],stage3ref,Params_chromosome_iit(params),
				 dbversion,Sequence_ntrimmed(queryseq),
				 showalignp,/*zerobasedp*/false,continuousp,diagnosticp,proteinmode,
				 invertmode,nointronlenp,wraplength);
	  for (i = 0; i < npaths; i++) {
	    stage3 = stage3array[i];
	    Stage3_free(&stage3);
	  }
	  FREE(stage3array);
	}
      }
      Sequence_free(&queryseq);
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

  while ((opt = getopt_long(argc,argv,"D:d:R:g:C:BL:x:1Xw:t:sec:SA39n:ZO5M:m:bEPQFNI:i:l:V?",
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
    case 'x': chimera_threshold = atof(optarg); break;
    case '1': maponlyp = true; break;
    case 'X': crossspeciesp = true; break;
    case 'w': mutationref = optarg; break;
    case 't': nworkers = atoi(optarg); break;
    case 's': aitstrainp = true; break;
    case 'e': extend_mismatch_p = true; break;
    case 'c': user_chrsubsetname = optarg; break;

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

    case 'E': cdna_exons_p = true; break;
    case 'P': protein_cdna_p = true; break;
    case 'Q': protein_genomic_p = true; break;
    case 'F': fulllengthp = true; break;
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
      indexdb = Indexdb_new_segment(Sequence_pointer(usersegment));
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

    if (aitstrainp == true) {
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
  Chrsubset_free(&chrsubset);
  IIT_free_mmapped(&chromosome_iit);

  return 0;
}
