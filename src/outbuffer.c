static char rcsid[] = "$Id: outbuffer.c 109763 2013-10-02 17:12:58Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "outbuffer.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef HAVE_PTHREAD
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>		/* Needed to define pthread_t on Solaris */
#endif
#include <pthread.h>
#endif

#include "assert.h"
#include "bool.h"
#include "mem.h"
#include "samheader.h"
#include "samflags.h"		/* For output types */

#ifdef GSNAP
#include "shortread.h"
#include "samprint.h"
#include "stage3hr.h"
#endif


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif


typedef struct RRlist_T *RRlist_T;
struct RRlist_T {
  int id;
  Result_T result;
  Request_T request;
  RRlist_T next;
};


#ifdef DEBUG1
static void
RRlist_dump (RRlist_T head, RRlist_T tail) {
  RRlist_T this;

  printf("head %p\n",head);
  for (this = head; this != NULL; this = head->next) {
    printf("%p: next %p\n",this,this->next);
  }
  printf("tail %p\n",tail);
  printf("\n");
  return;
}
#endif


/* Returns new tail */
static RRlist_T
RRlist_push (RRlist_T *head, RRlist_T tail, Request_T request, Result_T result) {
  RRlist_T new;

  new = (RRlist_T) MALLOC_OUT(sizeof(*new)); /* Called by worker thread */
  new->request = request;
  new->result = result;
  new->next = (RRlist_T) NULL;
  
  if (*head == NULL) {		/* Equivalent to tail == NULL, but using *head avoids having to set tail in RRlist_pop */
    *head = new;
  } else {
    tail->next = new;
  }

  return new;
}


/* Returns new head */
static RRlist_T
RRlist_pop (RRlist_T head, Request_T *request, Result_T *result) {
  RRlist_T newhead;

  *request = head->request;
  *result = head->result;

  newhead = head->next;

  FREE_OUT(head);		/* Called by outbuffer thread */
  return newhead;
}


static RRlist_T
RRlist_insert (RRlist_T list, int id, Request_T request, Result_T result) {
  RRlist_T *p;
  RRlist_T new;
  
  p = &list;
  while (*p != NULL && id > (*p)->id) {
    p = &(*p)->next;
  }

  new = (RRlist_T) MALLOC_OUT(sizeof(*new));
  new->id = id;
  new->request = request;
  new->result = result;
  
  new->next = *p;
  *p = new;
  return list;
}

/* Returns new head */
static RRlist_T
RRlist_pop_id (RRlist_T head, int *id, Request_T *request, Result_T *result) {
  RRlist_T newhead;

  *id = head->id;
  *request = head->request;
  *result = head->result;

  newhead = head->next;

  FREE_OUT(head);		/* Called by outbuffer thread */
  return newhead;
}




#define T Outbuffer_T
struct T {

#ifndef GSNAP
  Genome_T genome;
#endif

  Univ_IIT_T chromosome_iit;

  char *sevenway_root;
  bool appendp;

#ifdef GSNAP
  bool sam_headers_p;
  char *sam_read_group_id;
  char *sam_read_group_name;
  char *sam_read_group_library;
  char *sam_read_group_platform;
  int quality_shift;
  int nworkers;
  bool orderedp;
  int argc;
  char **argv;
  int optind;
#elif defined PMAP

#else
  bool sam_headers_p;
  bool sam_paired_p;
  char *sam_read_group_id;
  char *sam_read_group_name;
  char *sam_read_group_library;
  char *sam_read_group_platform;
  int quality_shift;
  int nworkers;
  bool orderedp;
  int argc;
  char **argv;
  int optind;
#endif

#ifdef GSNAP

  FILE *fp_nomapping_1;		/* N1 */
  FILE *fp_nomapping_2;		/* N2 */
  FILE *fp_halfmapping_uniq;	/* HU */
  FILE *fp_halfmapping_circular; /* HC */
  FILE *fp_halfmapping_transloc; /* HT */
  FILE *fp_halfmapping_mult;	 /* HM */
  FILE *fp_unpaired_uniq;	 /* UU */
  FILE *fp_unpaired_circular;	 /* UC */
  FILE *fp_unpaired_transloc;	 /* UT */
  FILE *fp_unpaired_mult;	 /* UM */
  FILE *fp_paired_uniq_circular; /* PC */
  FILE *fp_paired_uniq_inv;	 /* PI */
  FILE *fp_paired_uniq_scr;	 /* PS */
  FILE *fp_paired_uniq_long;	 /* PL */
  FILE *fp_paired_mult;		 /* PM */
  FILE *fp_concordant_uniq;	 /* CU */
  FILE *fp_concordant_circular;	 /* CC */
  FILE *fp_concordant_transloc;	 /* CT */
  FILE *fp_concordant_mult;	 /* CM */

  bool timingp;
  bool output_sam_p;
  Gobywriter_T gobywriter;

  bool fastq_format_p;
  bool clip_overlap_p;
  bool merge_samechr_p;

  bool invert_first_p;
  bool invert_second_p;
  Chrpos_T pairmax;

#else

  FILE *fp_nomapping;		/* N1 */
  FILE *fp_uniq;		/* UU */
  FILE *fp_circular;		/* UC */
  FILE *fp_transloc;		/* UT */
  FILE *fp_mult;		/* UM */

  bool chimeras_allowed_p;

  char *user_genomicseg;
  Sequence_T usersegment;

  char *dbversion;
  Chrsubset_T chrsubset;
  Univ_IIT_T contig_iit;
  IIT_T altstrain_iit;
  IIT_T map_iit;
  int *map_divint_crosstable;

  Printtype_T printtype;
  bool checksump;
  int chimera_margin;

  bool map_exons_p;
  bool map_bothstrands_p;
  bool print_comment_p;
  int nflanking;

  int proteinmode;
  int invertmode;
  bool nointronlenp;

  int wraplength;
  int ngap;
  int cds_startpos;

  bool fulllengthp;
  bool truncatep;
  bool strictp;
  bool diagnosticp;
  bool maponlyp;

  bool stage1debug;
  bool diag_debug;
  bool debug_graphic_p;

#endif

  int maxpaths_report;
  bool nofailsp;
  bool failsonlyp;
  bool fails_as_input_p;
  bool quiet_if_excessive_p;

#ifdef HAVE_PTHREAD
  pthread_mutex_t lock;
#endif

  unsigned int output_buffer_size;
  unsigned int nread;
  unsigned int ntotal;
  unsigned int nprocessed;

  RRlist_T head;
  RRlist_T tail;
  
#ifdef HAVE_PTHREAD
  pthread_cond_t result_avail_p;
#endif
};


/************************************************************************
 *   File routines
 ************************************************************************/

#ifdef GSNAP

static void
sevenway_open_single (T this) {
  char *filename;
  char *write_mode;

  if (this->appendp == true) {
    write_mode = "a";
  } else {
    write_mode = "w";
  }

  if (this->fails_as_input_p == true) {
    filename = (char *) CALLOC(strlen(this->sevenway_root)+strlen(".nomapping.fq")+1,sizeof(char));
    sprintf(filename,"%s.nomapping.fq",this->sevenway_root);
    if ((this->fp_nomapping_1 = fopen(filename,write_mode)) == NULL) {
      fprintf(stderr,"Cannot open file %s for writing\n",filename);
      exit(9);
    }
    FREE(filename);

  } else {
    filename = (char *) CALLOC(strlen(this->sevenway_root)+strlen(".nomapping")+1,sizeof(char));
    sprintf(filename,"%s.nomapping",this->sevenway_root);
    if ((this->fp_nomapping_1 = fopen(filename,write_mode)) == NULL) {
      fprintf(stderr,"Cannot open file %s for writing\n",filename);
      exit(9);
    }
    FREE(filename);

    if (this->output_sam_p == true && this->sam_headers_p == true) {
      SAM_header_print_HD(this->fp_nomapping_1,this->nworkers,this->orderedp);
      SAM_header_print_PG(this->fp_nomapping_1,this->argc,this->argv,this->optind);
      Univ_IIT_dump_sam(this->fp_nomapping_1,this->chromosome_iit,this->sam_read_group_id,this->sam_read_group_name,
			this->sam_read_group_library,this->sam_read_group_platform);
    }
  }

  filename = (char *) CALLOC(strlen(this->sevenway_root)+strlen(".unpaired_uniq")+1,sizeof(char));
  sprintf(filename,"%s.unpaired_uniq",this->sevenway_root);
  if ((this->fp_unpaired_uniq = fopen(filename,write_mode)) == NULL) {
    fprintf(stderr,"Cannot open file %s for writing\n",filename);
    exit(9);
  }
  FREE(filename);

  filename = (char *) CALLOC(strlen(this->sevenway_root)+strlen(".unpaired_circular")+1,sizeof(char));
  sprintf(filename,"%s.unpaired_circular",this->sevenway_root);
  if ((this->fp_unpaired_circular = fopen(filename,write_mode)) == NULL) {
    fprintf(stderr,"Cannot open file %s for writing\n",filename);
    exit(9);
  }
  FREE(filename);

  filename = (char *) CALLOC(strlen(this->sevenway_root)+strlen(".unpaired_transloc")+1,sizeof(char));
  sprintf(filename,"%s.unpaired_transloc",this->sevenway_root);
  if ((this->fp_unpaired_transloc = fopen(filename,write_mode)) == NULL) {
    fprintf(stderr,"Cannot open file %s for writing\n",filename);
    exit(9);
  }
  FREE(filename);

  filename = (char *) CALLOC(strlen(this->sevenway_root)+strlen(".unpaired_mult")+1,sizeof(char));
  sprintf(filename,"%s.unpaired_mult",this->sevenway_root);
  if ((this->fp_unpaired_mult = fopen(filename,write_mode)) == NULL) {
    fprintf(stderr,"Cannot open file %s for writing\n",filename);
    exit(9);
  }
  FREE(filename);

  if (this->output_sam_p == true && this->sam_headers_p == true) {
    SAM_header_print_HD(this->fp_unpaired_uniq,this->nworkers,this->orderedp);
    SAM_header_print_PG(this->fp_unpaired_uniq,this->argc,this->argv,this->optind);
    Univ_IIT_dump_sam(this->fp_unpaired_uniq,this->chromosome_iit,this->sam_read_group_id,this->sam_read_group_name,
		      this->sam_read_group_library,this->sam_read_group_platform);
    SAM_header_print_HD(this->fp_unpaired_circular,this->nworkers,this->orderedp);
    SAM_header_print_PG(this->fp_unpaired_circular,this->argc,this->argv,this->optind);
    Univ_IIT_dump_sam(this->fp_unpaired_circular,this->chromosome_iit,this->sam_read_group_id,this->sam_read_group_name,
		      this->sam_read_group_library,this->sam_read_group_platform);
    SAM_header_print_HD(this->fp_unpaired_transloc,this->nworkers,this->orderedp);
    SAM_header_print_PG(this->fp_unpaired_transloc,this->argc,this->argv,this->optind);
    Univ_IIT_dump_sam(this->fp_unpaired_transloc,this->chromosome_iit,this->sam_read_group_id,this->sam_read_group_name,
		      this->sam_read_group_library,this->sam_read_group_platform);
    SAM_header_print_HD(this->fp_unpaired_mult,this->nworkers,this->orderedp);
    SAM_header_print_PG(this->fp_unpaired_mult,this->argc,this->argv,this->optind);
    Univ_IIT_dump_sam(this->fp_unpaired_mult,this->chromosome_iit,this->sam_read_group_id,this->sam_read_group_name,
		      this->sam_read_group_library,this->sam_read_group_platform);
  }

  return;
}


static void
sevenway_open_paired (T this) {
  char *filename;
  char *write_mode;

  if (this->appendp == true) {
    write_mode = "a";
  } else {
    write_mode = "w";
  }

  if (this->fails_as_input_p == true) {
    if (this->fp_nomapping_1 == NULL) {
      filename = (char *) CALLOC(strlen(this->sevenway_root)+strlen(".nomapping.1.fq")+1,sizeof(char));
      sprintf(filename,"%s.nomapping.1.fq",this->sevenway_root);
      if ((this->fp_nomapping_1 = fopen(filename,write_mode)) == NULL) {
	fprintf(stderr,"Cannot open file %s for writing\n",filename);
	exit(9);
      }
      FREE(filename);
    }

    filename = (char *) CALLOC(strlen(this->sevenway_root)+strlen(".nomapping.2.fq")+1,sizeof(char));
    sprintf(filename,"%s.nomapping.2.fq",this->sevenway_root);
    if ((this->fp_nomapping_2 = fopen(filename,write_mode)) == NULL) {
      fprintf(stderr,"Cannot open file %s for writing\n",filename);
      exit(9);
    }
    FREE(filename);

  } else {
    if (this->fp_nomapping_1 == NULL) {
      filename = (char *) CALLOC(strlen(this->sevenway_root)+strlen(".nomapping")+1,sizeof(char));
      sprintf(filename,"%s.nomapping",this->sevenway_root);
      if ((this->fp_nomapping_1 = fopen(filename,write_mode)) == NULL) {
	fprintf(stderr,"Cannot open file %s for writing\n",filename);
	exit(9);
      }
      FREE(filename);

      if (this->output_sam_p == true && this->sam_headers_p == true) {
	SAM_header_print_HD(this->fp_nomapping_1,this->nworkers,this->orderedp);
	SAM_header_print_PG(this->fp_nomapping_1,this->argc,this->argv,this->optind);
	Univ_IIT_dump_sam(this->fp_nomapping_1,this->chromosome_iit,this->sam_read_group_id,this->sam_read_group_name,
			  this->sam_read_group_library,this->sam_read_group_platform);
      }
    }
  }

  filename = (char *) CALLOC(strlen(this->sevenway_root)+strlen(".halfmapping_uniq")+1,sizeof(char));
  sprintf(filename,"%s.halfmapping_uniq",this->sevenway_root);
  if ((this->fp_halfmapping_uniq = fopen(filename,write_mode)) == NULL) {
    fprintf(stderr,"Cannot open file %s for writing\n",filename);
    exit(9);
  }
  FREE(filename);

  filename = (char *) CALLOC(strlen(this->sevenway_root)+strlen(".halfmapping_circular")+1,sizeof(char));
  sprintf(filename,"%s.halfmapping_circular",this->sevenway_root);
  if ((this->fp_halfmapping_circular = fopen(filename,write_mode)) == NULL) {
    fprintf(stderr,"Cannot open file %s for writing\n",filename);
    exit(9);
  }
  FREE(filename);

  filename = (char *) CALLOC(strlen(this->sevenway_root)+strlen(".halfmapping_transloc")+1,sizeof(char));
  sprintf(filename,"%s.halfmapping_transloc",this->sevenway_root);
  if ((this->fp_halfmapping_transloc = fopen(filename,write_mode)) == NULL) {
    fprintf(stderr,"Cannot open file %s for writing\n",filename);
    exit(9);
  }
  FREE(filename);

  filename = (char *) CALLOC(strlen(this->sevenway_root)+strlen(".halfmapping_mult")+1,sizeof(char));
  sprintf(filename,"%s.halfmapping_mult",this->sevenway_root);
  if ((this->fp_halfmapping_mult = fopen(filename,write_mode)) == NULL) {
    fprintf(stderr,"Cannot open file %s for writing\n",filename);
    exit(9);
  }
  FREE(filename);

  filename = (char *) CALLOC(strlen(this->sevenway_root)+strlen(".paired_uniq_circular")+1,sizeof(char));
  sprintf(filename,"%s.paired_uniq_circular",this->sevenway_root);
  if ((this->fp_paired_uniq_circular = fopen(filename,write_mode)) == NULL) {
    fprintf(stderr,"Cannot open file %s for writing\n",filename);
    exit(9);
  }
  FREE(filename);

  filename = (char *) CALLOC(strlen(this->sevenway_root)+strlen(".paired_uniq_inv")+1,sizeof(char));
  sprintf(filename,"%s.paired_uniq_inv",this->sevenway_root);
  if ((this->fp_paired_uniq_inv = fopen(filename,write_mode)) == NULL) {
    fprintf(stderr,"Cannot open file %s for writing\n",filename);
    exit(9);
  }
  FREE(filename);

  filename = (char *) CALLOC(strlen(this->sevenway_root)+strlen(".paired_uniq_scr")+1,sizeof(char));
  sprintf(filename,"%s.paired_uniq_scr",this->sevenway_root);
  if ((this->fp_paired_uniq_scr = fopen(filename,write_mode)) == NULL) {
    fprintf(stderr,"Cannot open file %s for writing\n",filename);
    exit(9);
  }
  FREE(filename);

  filename = (char *) CALLOC(strlen(this->sevenway_root)+strlen(".paired_uniq_long")+1,sizeof(char));
  sprintf(filename,"%s.paired_uniq_long",this->sevenway_root);
  if ((this->fp_paired_uniq_long = fopen(filename,write_mode)) == NULL) {
    fprintf(stderr,"Cannot open file %s for writing\n",filename);
    exit(9);
  }
  FREE(filename);

  filename = (char *) CALLOC(strlen(this->sevenway_root)+strlen(".paired_mult")+1,sizeof(char));
  sprintf(filename,"%s.paired_mult",this->sevenway_root);
  if ((this->fp_paired_mult = fopen(filename,write_mode)) == NULL) {
    fprintf(stderr,"Cannot open file %s for writing\n",filename);
    exit(9);
  }
  FREE(filename);

  filename = (char *) CALLOC(strlen(this->sevenway_root)+strlen(".concordant_uniq")+1,sizeof(char));
  sprintf(filename,"%s.concordant_uniq",this->sevenway_root);
  if ((this->fp_concordant_uniq = fopen(filename,write_mode)) == NULL) {
    fprintf(stderr,"Cannot open file %s for writing\n",filename);
    exit(9);
  }
  FREE(filename);

  filename = (char *) CALLOC(strlen(this->sevenway_root)+strlen(".concordant_circular")+1,sizeof(char));
  sprintf(filename,"%s.concordant_circular",this->sevenway_root);
  if ((this->fp_concordant_circular = fopen(filename,write_mode)) == NULL) {
    fprintf(stderr,"Cannot open file %s for writing\n",filename);
    exit(9);
  }
  FREE(filename);

  filename = (char *) CALLOC(strlen(this->sevenway_root)+strlen(".concordant_transloc")+1,sizeof(char));
  sprintf(filename,"%s.concordant_transloc",this->sevenway_root);
  if ((this->fp_concordant_transloc = fopen(filename,write_mode)) == NULL) {
    fprintf(stderr,"Cannot open file %s for writing\n",filename);
    exit(9);
  }
  FREE(filename);

  filename = (char *) CALLOC(strlen(this->sevenway_root)+strlen(".concordant_mult")+1,sizeof(char));
  sprintf(filename,"%s.concordant_mult",this->sevenway_root);
  if ((this->fp_concordant_mult = fopen(filename,write_mode)) == NULL) {
    fprintf(stderr,"Cannot open file %s for writing\n",filename);
    exit(9);
  }
  FREE(filename);

  if (this->output_sam_p == true && this->sam_headers_p == true) {
    SAM_header_print_HD(this->fp_halfmapping_uniq,this->nworkers,this->orderedp);
    SAM_header_print_PG(this->fp_halfmapping_uniq,this->argc,this->argv,this->optind);
    Univ_IIT_dump_sam(this->fp_halfmapping_uniq,this->chromosome_iit,this->sam_read_group_id,this->sam_read_group_name,
		      this->sam_read_group_library,this->sam_read_group_platform);
    SAM_header_print_HD(this->fp_halfmapping_circular,this->nworkers,this->orderedp);
    SAM_header_print_PG(this->fp_halfmapping_circular,this->argc,this->argv,this->optind);
    Univ_IIT_dump_sam(this->fp_halfmapping_circular,this->chromosome_iit,this->sam_read_group_id,this->sam_read_group_name,
		      this->sam_read_group_library,this->sam_read_group_platform);
    SAM_header_print_HD(this->fp_halfmapping_transloc,this->nworkers,this->orderedp);
    SAM_header_print_PG(this->fp_halfmapping_transloc,this->argc,this->argv,this->optind);
    Univ_IIT_dump_sam(this->fp_halfmapping_transloc,this->chromosome_iit,this->sam_read_group_id,this->sam_read_group_name,
		      this->sam_read_group_library,this->sam_read_group_platform);
    SAM_header_print_HD(this->fp_halfmapping_mult,this->nworkers,this->orderedp);
    SAM_header_print_PG(this->fp_halfmapping_mult,this->argc,this->argv,this->optind);
    Univ_IIT_dump_sam(this->fp_halfmapping_mult,this->chromosome_iit,this->sam_read_group_id,this->sam_read_group_name,
		      this->sam_read_group_library,this->sam_read_group_platform);
    SAM_header_print_HD(this->fp_paired_uniq_circular,this->nworkers,this->orderedp);
    SAM_header_print_PG(this->fp_paired_uniq_circular,this->argc,this->argv,this->optind);
    Univ_IIT_dump_sam(this->fp_paired_uniq_circular,this->chromosome_iit,this->sam_read_group_id,this->sam_read_group_name,
		      this->sam_read_group_library,this->sam_read_group_platform);
    SAM_header_print_HD(this->fp_paired_uniq_inv,this->nworkers,this->orderedp);
    SAM_header_print_PG(this->fp_paired_uniq_inv,this->argc,this->argv,this->optind);
    Univ_IIT_dump_sam(this->fp_paired_uniq_inv,this->chromosome_iit,this->sam_read_group_id,this->sam_read_group_name,
		      this->sam_read_group_library,this->sam_read_group_platform);
    SAM_header_print_HD(this->fp_paired_uniq_scr,this->nworkers,this->orderedp);
    SAM_header_print_PG(this->fp_paired_uniq_scr,this->argc,this->argv,this->optind);
    Univ_IIT_dump_sam(this->fp_paired_uniq_scr,this->chromosome_iit,this->sam_read_group_id,this->sam_read_group_name,
		      this->sam_read_group_library,this->sam_read_group_platform);
    SAM_header_print_HD(this->fp_paired_uniq_long,this->nworkers,this->orderedp);
    SAM_header_print_PG(this->fp_paired_uniq_long,this->argc,this->argv,this->optind);
    Univ_IIT_dump_sam(this->fp_paired_uniq_long,this->chromosome_iit,this->sam_read_group_id,this->sam_read_group_name,
		      this->sam_read_group_library,this->sam_read_group_platform);
    SAM_header_print_HD(this->fp_paired_mult,this->nworkers,this->orderedp);
    SAM_header_print_PG(this->fp_paired_mult,this->argc,this->argv,this->optind);
    Univ_IIT_dump_sam(this->fp_paired_mult,this->chromosome_iit,this->sam_read_group_id,this->sam_read_group_name,
		      this->sam_read_group_library,this->sam_read_group_platform);
    SAM_header_print_HD(this->fp_concordant_uniq,this->nworkers,this->orderedp);
    SAM_header_print_PG(this->fp_concordant_uniq,this->argc,this->argv,this->optind);
    Univ_IIT_dump_sam(this->fp_concordant_uniq,this->chromosome_iit,this->sam_read_group_id,this->sam_read_group_name,
		      this->sam_read_group_library,this->sam_read_group_platform);
    SAM_header_print_HD(this->fp_concordant_circular,this->nworkers,this->orderedp);
    SAM_header_print_PG(this->fp_concordant_circular,this->argc,this->argv,this->optind);
    Univ_IIT_dump_sam(this->fp_concordant_circular,this->chromosome_iit,this->sam_read_group_id,this->sam_read_group_name,
		      this->sam_read_group_library,this->sam_read_group_platform);
    SAM_header_print_HD(this->fp_concordant_transloc,this->nworkers,this->orderedp);
    SAM_header_print_PG(this->fp_concordant_transloc,this->argc,this->argv,this->optind);
    Univ_IIT_dump_sam(this->fp_concordant_transloc,this->chromosome_iit,this->sam_read_group_id,this->sam_read_group_name,
		      this->sam_read_group_library,this->sam_read_group_platform);
    SAM_header_print_HD(this->fp_concordant_mult,this->nworkers,this->orderedp);
    SAM_header_print_PG(this->fp_concordant_mult,this->argc,this->argv,this->optind);
    Univ_IIT_dump_sam(this->fp_concordant_mult,this->chromosome_iit,this->sam_read_group_id,this->sam_read_group_name,
		      this->sam_read_group_library,this->sam_read_group_platform);
  }

  return;
}

static void
sevenway_close (T this) {
  fclose(this->fp_unpaired_uniq);
  fclose(this->fp_unpaired_circular);
  fclose(this->fp_unpaired_transloc);
  fclose(this->fp_unpaired_mult);
  if (this->fp_nomapping_2 != NULL) {
    fclose(this->fp_nomapping_2);
  }
  if (this->fp_nomapping_1 != NULL) {
    fclose(this->fp_nomapping_1);
  }
  if (this->fp_halfmapping_uniq != NULL) {
    fclose(this->fp_halfmapping_uniq);
    fclose(this->fp_halfmapping_circular);
    fclose(this->fp_halfmapping_transloc);
    fclose(this->fp_halfmapping_mult);
    fclose(this->fp_paired_uniq_long);
    fclose(this->fp_paired_uniq_scr);
    fclose(this->fp_paired_uniq_inv);
    fclose(this->fp_paired_uniq_circular);
    fclose(this->fp_paired_mult);
    fclose(this->fp_concordant_uniq);
    fclose(this->fp_concordant_circular);
    fclose(this->fp_concordant_transloc);
    fclose(this->fp_concordant_mult);
  }

  return;
}

#else

/* GMAP version */

static void
print_gff_header (FILE *fp, int argc, char **argv, int optind) {
  char **argstart;
  int c;

  fprintf(fp,"##gff-version   3\n");
  fprintf(fp,"# Generated by GMAP version %s using call: ",PACKAGE_VERSION);
  argstart = &(argv[-optind]);
  for (c = 0; c < argc + optind; c++) {
    fprintf(fp," %s",argstart[c]);
  }
  fprintf(fp,"\n");
  return;
}


/* Taken from Univ_IIT_dump_sam */
static void
dump_sam_usersegment (FILE *fp, Sequence_T usersegment,
		      char *sam_read_group_id, char *sam_read_group_name,
		      char *sam_read_group_library, char *sam_read_group_platform) {

  fprintf(fp,"@SQ\tSN:%s",Sequence_accession(usersegment));
  fprintf(fp,"\tLN:%u",Sequence_fulllength(usersegment));
  fprintf(fp,"\n");

  if (sam_read_group_id != NULL) {
    fprintf(fp,"@RG\tID:%s",sam_read_group_id);
    if (sam_read_group_platform != NULL) {
      fprintf(fp,"\tPL:%s",sam_read_group_platform);
    }
    if (sam_read_group_library != NULL) {
      fprintf(fp,"\tLB:%s",sam_read_group_library);
    }
    fprintf(fp,"\tSM:%s",sam_read_group_name);
    fprintf(fp,"\n");
  }

  return;
}



static void
sevenway_open (T this, int argc, char **argv, int optind) {
  char *filename;
  char *write_mode;

  if (this->appendp == true) {
    write_mode = "a";
  } else {
    write_mode = "w";
  }

  if (this->fails_as_input_p == true) {
    filename = (char *) CALLOC(strlen(this->sevenway_root)+strlen(".nomapping.fa")+1,sizeof(char));
    sprintf(filename,"%s.nomapping.fa",this->sevenway_root);
    if ((this->fp_nomapping = fopen(filename,write_mode)) == NULL) {
      fprintf(stderr,"Cannot open file %s for writing\n",filename);
      exit(9);
    }
    FREE(filename);

  } else {
    filename = (char *) CALLOC(strlen(this->sevenway_root)+strlen(".nomapping")+1,sizeof(char));
    sprintf(filename,"%s.nomapping",this->sevenway_root);
    if ((this->fp_nomapping = fopen(filename,write_mode)) == NULL) {
      fprintf(stderr,"Cannot open file %s for writing\n",filename);
      exit(9);
    }
    FREE(filename);

    if (this->printtype == GFF3_GENE || this->printtype == GFF3_MATCH_CDNA || this->printtype == GFF3_MATCH_EST) {
      print_gff_header(this->fp_nomapping,argc,argv,optind);
#ifndef PMAP
    } else if (this->printtype == SAM && this->sam_headers_p == true) {
      if (this->usersegment != NULL) {
	dump_sam_usersegment(this->fp_nomapping,this->usersegment,this->sam_read_group_id,this->sam_read_group_name,
			     this->sam_read_group_library,this->sam_read_group_platform);
      } else {
	Univ_IIT_dump_sam(this->fp_nomapping,this->chromosome_iit,this->sam_read_group_id,this->sam_read_group_name,
			  this->sam_read_group_library,this->sam_read_group_platform);
      }
#endif
    }
  }

  filename = (char *) CALLOC(strlen(this->sevenway_root)+strlen(".uniq")+1,sizeof(char));
  sprintf(filename,"%s.uniq",this->sevenway_root);
  if ((this->fp_uniq = fopen(filename,write_mode)) == NULL) {
    fprintf(stderr,"Cannot open file %s for writing\n",filename);
    exit(9);
  }
  FREE(filename);

  filename = (char *) CALLOC(strlen(this->sevenway_root)+strlen(".circular")+1,sizeof(char));
  sprintf(filename,"%s.circular",this->sevenway_root);
  if ((this->fp_circular = fopen(filename,write_mode)) == NULL) {
    fprintf(stderr,"Cannot open file %s for writing\n",filename);
    exit(9);
  }
  FREE(filename);

  if (this->chimeras_allowed_p == true) {
    filename = (char *) CALLOC(strlen(this->sevenway_root)+strlen(".transloc")+1,sizeof(char));
    sprintf(filename,"%s.transloc",this->sevenway_root);
    if ((this->fp_transloc = fopen(filename,write_mode)) == NULL) {
      fprintf(stderr,"Cannot open file %s for writing\n",filename);
      exit(9);
    }
    FREE(filename);
  }

  filename = (char *) CALLOC(strlen(this->sevenway_root)+strlen(".mult")+1,sizeof(char));
  sprintf(filename,"%s.mult",this->sevenway_root);
  if ((this->fp_mult = fopen(filename,write_mode)) == NULL) {
    fprintf(stderr,"Cannot open file %s for writing\n",filename);
    exit(9);
  }
  FREE(filename);

  if (this->printtype == GFF3_GENE || this->printtype == GFF3_MATCH_CDNA || this->printtype == GFF3_MATCH_EST) {
    print_gff_header(this->fp_uniq,argc,argv,optind);
    print_gff_header(this->fp_circular,argc,argv,optind);
    print_gff_header(this->fp_mult,argc,argv,optind);

#ifndef PMAP
  } else if (this->printtype == SAM && this->sam_headers_p == true) {
    if (this->usersegment != NULL) {
      dump_sam_usersegment(this->fp_uniq,this->usersegment,
			   this->sam_read_group_id,this->sam_read_group_name,
			   this->sam_read_group_library,this->sam_read_group_platform);
      dump_sam_usersegment(this->fp_circular,this->usersegment,
			   this->sam_read_group_id,this->sam_read_group_name,
			   this->sam_read_group_library,this->sam_read_group_platform);
      dump_sam_usersegment(this->fp_mult,this->usersegment,
			   this->sam_read_group_id,this->sam_read_group_name,
			   this->sam_read_group_library,this->sam_read_group_platform);
    } else {
      Univ_IIT_dump_sam(this->fp_uniq,this->chromosome_iit,this->sam_read_group_id,this->sam_read_group_name,
			this->sam_read_group_library,this->sam_read_group_platform);
      Univ_IIT_dump_sam(this->fp_circular,this->chromosome_iit,this->sam_read_group_id,this->sam_read_group_name,
			this->sam_read_group_library,this->sam_read_group_platform);
      Univ_IIT_dump_sam(this->fp_mult,this->chromosome_iit,this->sam_read_group_id,this->sam_read_group_name,
			this->sam_read_group_library,this->sam_read_group_platform);
    }
#endif
  }

  return;
}

static void
sevenway_close (T this) {
  fclose(this->fp_mult);
  fclose(this->fp_circular);
  fclose(this->fp_uniq);
  if (this->chimeras_allowed_p == true) {
    fclose(this->fp_transloc);
  }
  fclose(this->fp_nomapping);
  return;
}

#endif



#ifdef GSNAP

T
Outbuffer_new (unsigned int output_buffer_size, unsigned int nread, char *sevenway_root, bool appendp, Univ_IIT_T chromosome_iit,
	       bool timingp, bool output_sam_p, bool sam_headers_p, char *sam_read_group_id, char *sam_read_group_name,
	       char *sam_read_group_library, char *sam_read_group_platform,
	       int nworkers, bool orderedp,
	       Gobywriter_T gobywriter, bool nofailsp, bool failsonlyp, bool fails_as_input_p,
	       bool fastq_format_p, bool clip_overlap_p, bool merge_samechr_p,
	       int maxpaths_report, bool quiet_if_excessive_p, int quality_shift,
	       bool invert_first_p, bool invert_second_p, Chrpos_T pairmax,
	       int argc, char **argv, int optind) {
  T new = (T) MALLOC(sizeof(*new));
  FILE *fp_capture = NULL, *fp_ignore = NULL;

  new->chromosome_iit = chromosome_iit;

  new->fp_nomapping_1 = NULL;
  new->fp_nomapping_2 = NULL;
  new->fp_halfmapping_uniq = NULL;
  new->fp_halfmapping_circular = NULL;
  new->fp_halfmapping_transloc = NULL;
  new->fp_halfmapping_mult = NULL;
  new->fp_unpaired_uniq = NULL;
  new->fp_unpaired_circular = NULL;
  new->fp_unpaired_transloc = NULL;
  new->fp_unpaired_mult = NULL;
  new->fp_paired_uniq_circular = NULL;
  new->fp_paired_uniq_inv = NULL;
  new->fp_paired_uniq_scr = NULL;
  new->fp_paired_uniq_long = NULL;
  new->fp_paired_mult = NULL;
  new->fp_concordant_uniq = NULL;
  new->fp_concordant_circular = NULL;
  new->fp_concordant_transloc = NULL;
  new->fp_concordant_mult = NULL;
  
  new->sevenway_root = sevenway_root;
  new->appendp = appendp;

  new->timingp = timingp;
  new->output_sam_p = output_sam_p;
  new->sam_headers_p = sam_headers_p;
  new->sam_read_group_id = sam_read_group_id;
  new->sam_read_group_name = sam_read_group_name;
  new->sam_read_group_library = sam_read_group_library;
  new->sam_read_group_platform = sam_read_group_platform;
  new->nworkers = nworkers;
  new->orderedp = orderedp;
  new->argc = argc;
  new->argv = argv;
  new->optind = optind;

  new->gobywriter = gobywriter;

  new->nofailsp = nofailsp;
  new->failsonlyp = failsonlyp;
  new->fails_as_input_p = fails_as_input_p;
  new->fastq_format_p = fastq_format_p;
  new->clip_overlap_p = clip_overlap_p;
  new->merge_samechr_p = merge_samechr_p;

  new->maxpaths_report = maxpaths_report;
  new->quiet_if_excessive_p = quiet_if_excessive_p;

  new->quality_shift = quality_shift;
  new->invert_first_p = invert_first_p;
  new->invert_second_p = invert_second_p;
  new->pairmax = pairmax;

#ifdef HAVE_PTHREAD
  pthread_mutex_init(&new->lock,NULL);
#endif

  new->output_buffer_size = output_buffer_size;
  new->nread = nread;
  new->ntotal = (unsigned int) -1U; /* Set to infinity until all reads are input */
  new->nprocessed = 0;

  new->head = (RRlist_T) NULL;
  new->tail = (RRlist_T) NULL;

#ifdef HAVE_PTHREAD
  pthread_cond_init(&new->result_avail_p,NULL);
#endif

  /* Initialize output streams */
  if (new->gobywriter != NULL) {
    Goby_file_handles(&fp_capture,&fp_ignore,new->gobywriter);
    new->fp_nomapping_1 = fp_ignore;
    new->fp_nomapping_2 = fp_ignore;
    new->fp_halfmapping_uniq = fp_capture;
    new->fp_halfmapping_circular = fp_capture;
    new->fp_halfmapping_transloc = fp_capture;
    new->fp_halfmapping_mult = fp_capture;
    new->fp_unpaired_uniq = fp_capture;
    new->fp_unpaired_circular = fp_capture;
    new->fp_unpaired_transloc = fp_capture;
    new->fp_unpaired_mult = fp_capture;
    new->fp_paired_uniq_circular = fp_capture;
    new->fp_paired_uniq_inv = fp_capture;
    new->fp_paired_uniq_scr = fp_capture;
    new->fp_paired_uniq_long = fp_capture;
    new->fp_paired_mult = fp_capture;
    new->fp_concordant_uniq = fp_capture;
    new->fp_concordant_circular = fp_capture;
    new->fp_concordant_transloc = fp_capture;
    new->fp_concordant_mult = fp_capture;

  } else if (sevenway_root != NULL) {
    sevenway_open_single(new);

  } else {
    new->fp_nomapping_1 = stdout;
    new->fp_nomapping_2 = stdout;
    new->fp_halfmapping_uniq = stdout;
    new->fp_halfmapping_circular = stdout;
    new->fp_halfmapping_transloc = stdout;
    new->fp_halfmapping_mult = stdout;
    new->fp_unpaired_uniq = stdout;
    new->fp_unpaired_circular = stdout;
    new->fp_unpaired_transloc = stdout;
    new->fp_unpaired_mult = stdout;
    new->fp_paired_uniq_circular = stdout;
    new->fp_paired_uniq_inv = stdout;
    new->fp_paired_uniq_scr = stdout;
    new->fp_paired_uniq_long = stdout;
    new->fp_paired_mult = stdout;
    new->fp_concordant_uniq = stdout;
    new->fp_concordant_circular = stdout;
    new->fp_concordant_transloc = stdout;
    new->fp_concordant_mult = stdout;

    if (output_sam_p == true && sam_headers_p == true) {
      if (fails_as_input_p == true) {
	/* Don't print chromosomes */
      } else {
	SAM_header_print_HD(stdout,nworkers,orderedp);
	SAM_header_print_PG(stdout,argc,argv,optind);
	Univ_IIT_dump_sam(stdout,chromosome_iit,sam_read_group_id,sam_read_group_name,
			  sam_read_group_library,sam_read_group_platform);
      }
    }
  }

  return new;
}

#else

T
Outbuffer_new (unsigned int output_buffer_size, unsigned int nread, char *sevenway_root, bool appendp,
	       bool chimeras_allowed_p, char *user_genomicseg, Sequence_T usersegment,
	       char *dbversion, Genome_T genome, Univ_IIT_T chromosome_iit,
	       Chrsubset_T chrsubset, Univ_IIT_T contig_iit, IIT_T altstrain_iit, IIT_T map_iit,
	       int *map_divint_crosstable, Printtype_T printtype, bool checksump, int chimera_margin,
#ifndef PMAP
	       bool sam_headers_p, int quality_shift, bool sam_paired_p,
	       char *sam_read_group_id, char *sam_read_group_name,
	       char *sam_read_group_library, char *sam_read_group_platform,
	       int nworkers, bool orderedp,
#endif
	       bool nofailsp, bool failsonlyp, bool fails_as_input_p, int maxpaths_report, bool quiet_if_excessive_p,
	       bool map_exons_p, bool map_bothstrands_p, bool print_comment_p, int nflanking,
	       int proteinmode, int invertmode, bool nointronlenp, int wraplength,
	       int ngap, int cds_startpos,
	       bool fulllengthp, bool truncatep, bool strictp, bool diagnosticp, bool maponlyp,
	       bool stage1debug, bool diag_debug, bool debug_graphic_p,
	       int argc, char **argv, int optind) {

  T new = (T) MALLOC(sizeof(*new));

  new->chimeras_allowed_p = chimeras_allowed_p;

  new->user_genomicseg = user_genomicseg;
  new->usersegment = usersegment;

  new->dbversion = dbversion;
  new->genome = genome;
  new->chromosome_iit = chromosome_iit;
  new->chrsubset = chrsubset;
  new->contig_iit = contig_iit;
  new->altstrain_iit = altstrain_iit;
  new->map_iit = map_iit;
  new->map_divint_crosstable = map_divint_crosstable;

  new->printtype = printtype;
  new->checksump = checksump;
  new->chimera_margin = chimera_margin;

  new->sevenway_root = sevenway_root;
  new->appendp = appendp;

  new->fp_nomapping = NULL;
  new->fp_uniq = NULL;
  new->fp_circular = NULL;
  new->fp_transloc = NULL;
  new->fp_mult = NULL;
  
#ifndef PMAP
  new->sam_headers_p = sam_headers_p;
  new->quality_shift = quality_shift;
  new->sam_paired_p = sam_paired_p;
  new->sam_read_group_id = sam_read_group_id;
  new->sam_read_group_name = sam_read_group_name;
  new->sam_read_group_library = sam_read_group_library;
  new->sam_read_group_platform = sam_read_group_platform;
  new->nworkers = nworkers;
  new->orderedp = orderedp;
  new->argc = argc;
  new->argv = argv;
  new->optind = optind;
#endif

  new->nofailsp = nofailsp;
  new->failsonlyp = failsonlyp;
  new->fails_as_input_p = fails_as_input_p;
  new->maxpaths_report = maxpaths_report;
  new->quiet_if_excessive_p = quiet_if_excessive_p;

  new->map_exons_p = map_exons_p;
  new->map_bothstrands_p = map_bothstrands_p;
  new->print_comment_p = print_comment_p;

  new->nflanking = nflanking;
  new->proteinmode = proteinmode;
  new->invertmode = invertmode;
  new->nointronlenp = nointronlenp;

  new->wraplength = wraplength;
  new->ngap = ngap;
  new->cds_startpos = cds_startpos;

  new->fulllengthp = fulllengthp;
  new->truncatep = truncatep;
  new->strictp = strictp;
  new->diagnosticp = diagnosticp;
  new->maponlyp = maponlyp;

  new->stage1debug = stage1debug;
  new->diag_debug = diag_debug;
  new->debug_graphic_p = debug_graphic_p;

#ifdef HAVE_PTHREAD
  pthread_mutex_init(&new->lock,NULL);
#endif

  new->output_buffer_size = output_buffer_size;
  new->nread = nread;
  new->ntotal = (unsigned int) -1U; /* Set to infinity until all reads are input */
  new->nprocessed = 0;

  new->head = (RRlist_T) NULL;
  new->tail = (RRlist_T) NULL;

#ifdef HAVE_PTHREAD
  pthread_cond_init(&new->result_avail_p,NULL);
#endif

  /* Initialize output streams */
  if (sevenway_root != NULL) {
    sevenway_open(new,argc,argv,optind);

  } else {
    new->fp_nomapping = stdout;
    new->fp_uniq = stdout;
    new->fp_circular = stdout;
    new->fp_transloc = stdout;
    new->fp_mult = stdout;

#ifndef PMAP
    if (printtype == SAM && sam_headers_p == true) {
      if (fails_as_input_p == true) {
	/* Don't print chromosomes */
      } else if (usersegment != NULL) {
	dump_sam_usersegment(stdout,usersegment,sam_read_group_id,sam_read_group_name,
			     sam_read_group_library,sam_read_group_platform);
      } else {
	SAM_header_print_HD(stdout,nworkers,orderedp);
	SAM_header_print_PG(stdout,argc,argv,optind);
	Univ_IIT_dump_sam(stdout,chromosome_iit,sam_read_group_id,sam_read_group_name,
			  sam_read_group_library,sam_read_group_platform);
      }
    }
#endif

  }

  return new;
}

#endif

void
Outbuffer_free (T *old) {
  if (*old) {
    if ((*old)->sevenway_root != NULL) {
      sevenway_close(*old);
    }

#ifdef HAVE_PTHREAD
    pthread_cond_destroy(&(*old)->result_avail_p);
    pthread_mutex_destroy(&(*old)->lock);
#endif

    FREE(*old);
  }
  return;
}



unsigned int
Outbuffer_nread (T this) {
  return this->nread;
}



void
Outbuffer_add_nread (T this, unsigned int nread) {

#ifdef HAVE_PTHREAD
  pthread_mutex_lock(&this->lock);
#endif

  if (nread == 0) {
    /* Finished reading, so able to determine total reads in input */
    this->ntotal = this->nread;
    debug(fprintf(stderr,"__Outbuffer_add_nread added 0 reads, so setting ntotal to be %u\n",this->ntotal));

#ifdef HAVE_PTHREAD
    pthread_cond_signal(&this->result_avail_p);
#endif

  } else {
    this->nread += nread;
    debug(fprintf(stderr,"__Outbuffer_add_nread added %d read, now %d\n",nread,this->nread));
  }

#ifdef HAVE_PTHREAD
  pthread_mutex_unlock(&this->lock);
#endif

  return;
}


void
Outbuffer_put_result (T this, Result_T result, Request_T request) {

#ifdef HAVE_PTHREAD
  pthread_mutex_lock(&this->lock);
#endif

  this->tail = RRlist_push(&this->head,this->tail,request,result);
  debug1(RRlist_dump(this->head,this->tail));
  this->nprocessed += 1;

#ifdef HAVE_PTHREAD
  pthread_cond_signal(&this->result_avail_p);
  pthread_mutex_unlock(&this->lock);
#endif

  return;
}



#ifdef GSNAP

/************************************************************************
 *   Print routines and threads for GSNAP
 ************************************************************************/

static void
print_query_singleend (T this, FILE *fp, Request_T request) {
  Shortread_T queryseq1;

  queryseq1 = Request_queryseq1(request);

  if (this->fastq_format_p == true) {
    Shortread_print_query_singleend_fastq(fp,queryseq1);
  } else {
    Shortread_print_query_singleend_fasta(fp,queryseq1);
  }

  return;
}

static void
print_header_singleend (T this, FILE *fp, Request_T request, bool translocationp, int npaths) {
  Shortread_T queryseq1;

  queryseq1 = Request_queryseq1(request);

  fprintf(fp,">");
  Shortread_print_oneline(fp,queryseq1);
  fprintf(fp,"\t%d",npaths);
  if (translocationp == true) {
    fprintf(fp," (transloc)");
  }

  /* No sequence inversion on single-end reads */
  if (Shortread_quality_string(queryseq1) != NULL) {
    fprintf(fp,"\t");
    Shortread_print_quality(fp,queryseq1,/*hardclip_low*/0,/*hardclip_high*/0,
			    this->quality_shift,/*show_chopped_p*/true);
  }

  fprintf(fp,"\t");
  Shortread_print_header(fp,queryseq1,/*queryseq2*/NULL);
  /* fprintf(fp,"\n"); -- included in header */

  return;
}


static void
print_result_sam (T this, Result_T result, Request_T request) {
  Resulttype_T resulttype;
  Shortread_T queryseq1;
  Stage3end_T *stage3array, stage3;
  Chrpos_T chrpos;
  int ignore = 0;
  int npaths, pathnum, first_absmq, second_absmq;
  FILE *fp;
  char *abbrev;

  resulttype = Result_resulttype(result);

  if (resulttype == SINGLEEND_NOMAPPING) {
    if (this->nofailsp == true) {
      /* Skip */
    } else {
      if (this->fails_as_input_p == true) {
	print_query_singleend(this,this->fp_nomapping_1,request);
      } else {
	queryseq1 = Request_queryseq1(request);
	SAM_print_nomapping(this->fp_nomapping_1,ABBREV_NOMAPPING_1,
			    queryseq1,/*mate*/NULL,/*acc1*/Shortread_accession(queryseq1),
			    /*acc2*/NULL,this->chromosome_iit,resulttype,
			    /*first_read_p*/true,/*nhits_mate*/0,/*mate_chrpos*/0U,
			    this->quality_shift,this->sam_read_group_id,this->invert_first_p,this->invert_second_p);
      }
    }

  } else if (resulttype == SINGLEEND_UNIQ) {
    stage3array = (Stage3end_T *) Result_array(&npaths,&first_absmq,&second_absmq,result);

    if (this->failsonlyp == true) {
      /* Skip */
    } else {
      queryseq1 = Request_queryseq1(request);
      /* Stage3end_eval_and_sort(stage3array,npaths,this->maxpaths_report,queryseq1); */

      stage3 = stage3array[0];
      chrpos = SAM_compute_chrpos(/*hardclip_low*/&ignore,/*hardclip_high*/&ignore,
				  /*clipdir*/0,/*hardclip5*/0,/*hardclip3*/0,/*first_read_p*/true,
				  stage3,Stage3end_substring_low(stage3),Shortread_fulllength(queryseq1));
      if (Stage3end_circularpos(stage3) > 0) {
	fp = this->fp_unpaired_circular;
	abbrev = ABBREV_UNPAIRED_CIRCULAR;
      } else {
	fp = this->fp_unpaired_uniq;
	abbrev = ABBREV_UNPAIRED_UNIQ;
      }
      SAM_print(fp,abbrev,stage3,/*mate*/NULL,/*acc1*/Shortread_accession(queryseq1),/*acc2*/NULL,
		/*pathnum*/1,npaths,Stage3end_absmq_score(stage3array[0]),first_absmq,second_absmq,
		Stage3end_mapq_score(stage3array[0]),
		this->chromosome_iit,queryseq1,/*queryseq2*/NULL,
		/*pairedlength*/0,chrpos,/*mate_chrpos*/0U,
		/*clipdir*/0,/*hardclip_low*/0,/*hardclip_high*/0,resulttype,
		/*first_read_p*/true,/*npaths_mate*/0,this->quality_shift,
		this->sam_read_group_id,this->invert_first_p,this->invert_second_p,
		this->merge_samechr_p);
    }

  } else if (resulttype == SINGLEEND_TRANSLOC) {
    stage3array = (Stage3end_T *) Result_array(&npaths,&first_absmq,&second_absmq,result);

    if (this->failsonlyp == true) {
      /* Skip */

    } else if (this->quiet_if_excessive_p && npaths > this->maxpaths_report) {
      queryseq1 = Request_queryseq1(request);
      /* Stage3end_eval_and_sort(stage3array,npaths,this->maxpaths_report,queryseq1); */
      SAM_print_nomapping(this->fp_unpaired_transloc,ABBREV_UNPAIRED_TRANSLOC,
			  queryseq1,/*mate*/NULL,/*acc1*/Shortread_accession(queryseq1),
			  /*acc2*/NULL,this->chromosome_iit,resulttype,
			  /*first_read_p*/true,/*nhits_mate*/0,/*mate_chrpos*/0U,
			  this->quality_shift,this->sam_read_group_id,this->invert_first_p,this->invert_second_p);

    } else {
      queryseq1 = Request_queryseq1(request);
      /* Stage3end_eval_and_sort(stage3array,npaths,this->maxpaths_report,queryseq1); */
      for (pathnum = 1; pathnum <= npaths && pathnum <= this->maxpaths_report; pathnum++) {

	stage3 = stage3array[pathnum-1];
	chrpos = SAM_compute_chrpos(/*hardclip_low*/&ignore,/*hardclip_high*/&ignore,
				    /*clipdir*/0,/*hardclip5*/0,/*hardclip3*/0,/*first_read_p*/true,
				    stage3,Stage3end_substring_low(stage3),Shortread_fulllength(queryseq1));
	SAM_print(this->fp_unpaired_transloc,ABBREV_UNPAIRED_TRANSLOC,
		  stage3,/*mate*/NULL,/*acc1*/Shortread_accession(queryseq1),
		  /*acc2*/NULL,pathnum,npaths,
		  Stage3end_absmq_score(stage3array[pathnum-1]),first_absmq,second_absmq,
		  Stage3end_mapq_score(stage3array[pathnum-1]),
		  this->chromosome_iit,queryseq1,/*queryseq2*/NULL,
		  /*pairedlength*/0,chrpos,/*mate_chrpos*/0U,
		  /*clipdir*/0,/*hardclip_low*/0,/*hardclip_high*/0,resulttype,
		  /*first_read_p*/true,/*npaths_mate*/0,this->quality_shift,
		  this->sam_read_group_id,this->invert_first_p,this->invert_second_p,
		  this->merge_samechr_p);
      }
    }

  } else if (resulttype == SINGLEEND_MULT) {
    stage3array = (Stage3end_T *) Result_array(&npaths,&first_absmq,&second_absmq,result);

    if (this->failsonlyp == true) {
      /* Skip */

    } else if (this->quiet_if_excessive_p && npaths > this->maxpaths_report) {
      queryseq1 = Request_queryseq1(request);
      /* Stage3end_eval_and_sort(stage3array,npaths,this->maxpaths_report,queryseq1); */
      SAM_print_nomapping(this->fp_unpaired_mult,ABBREV_UNPAIRED_MULT,
			  queryseq1,/*mate*/NULL,/*acc1*/Shortread_accession(queryseq1),
			  /*acc2*/NULL,this->chromosome_iit,resulttype,
			  /*first_read_p*/true,/*nhits_mate*/0,/*mate_chrpos*/0U,
			  this->quality_shift,this->sam_read_group_id,this->invert_first_p,this->invert_second_p);

    } else {
      queryseq1 = Request_queryseq1(request);
      /* Stage3end_eval_and_sort(stage3array,npaths,this->maxpaths_report,queryseq1); */
      for (pathnum = 1; pathnum <= npaths && pathnum <= this->maxpaths_report; pathnum++) {

	stage3 = stage3array[pathnum-1];
	chrpos = SAM_compute_chrpos(/*hardclip_low*/&ignore,/*hardclip_high*/&ignore,
				    /*clipdir*/0,/*hardclip5*/0,/*hardclip3*/0,/*first_read_p*/true,
				    stage3,Stage3end_substring_low(stage3),Shortread_fulllength(queryseq1));
	SAM_print(this->fp_unpaired_mult,ABBREV_UNPAIRED_MULT,
		  stage3,/*mate*/NULL,/*acc1*/Shortread_accession(queryseq1),
		  /*acc2*/NULL,pathnum,npaths,
		  Stage3end_absmq_score(stage3array[pathnum-1]),first_absmq,second_absmq,
		  Stage3end_mapq_score(stage3array[pathnum-1]),
		  this->chromosome_iit,queryseq1,/*queryseq2*/NULL,
		  /*pairedlength*/0,chrpos,/*mate_chrpos*/0U,
		  /*clipdir*/0,/*hardclip_low*/0,/*hardclip_high*/0,resulttype,
		  /*first_read_p*/true,/*npaths_mate*/0,this->quality_shift,
		  this->sam_read_group_id,this->invert_first_p,this->invert_second_p,
		  this->merge_samechr_p);
      }
    }

  } else {
    if (this->fp_concordant_uniq == NULL) {
      sevenway_open_paired(this);
    }
    SAM_print_paired(result,resulttype,this->chromosome_iit,
		     Request_queryseq1(request),Request_queryseq2(request),
		     this->invert_first_p,this->invert_second_p,
		     this->nofailsp,this->failsonlyp,this->fails_as_input_p,this->fastq_format_p,
		     this->clip_overlap_p,this->merge_samechr_p,this->quality_shift,this->sam_read_group_id,
		     this->fp_nomapping_1,this->fp_nomapping_2,
		     this->fp_unpaired_uniq,this->fp_unpaired_circular,
		     this->fp_unpaired_transloc,this->fp_unpaired_mult,
		     this->fp_halfmapping_uniq,this->fp_halfmapping_circular,
		     this->fp_halfmapping_transloc,this->fp_halfmapping_mult,
		     this->fp_paired_uniq_circular,this->fp_paired_uniq_inv,this->fp_paired_uniq_scr,
		     this->fp_paired_uniq_long,this->fp_paired_mult,
		     this->fp_concordant_uniq,this->fp_concordant_circular,
		     this->fp_concordant_transloc,this->fp_concordant_mult);
  }

  return;
}


static void
print_result_gsnap (T this, Result_T result, Request_T request) {
  Resulttype_T resulttype;
  Shortread_T queryseq1;
  Stage3end_T *stage3array, stage3;
  int npaths, pathnum, first_absmq, second_absmq;
  FILE *fp;

  resulttype = Result_resulttype(result);

  if (resulttype == SINGLEEND_NOMAPPING) {
    if (this->nofailsp == true) {
      /* Skip */
    } else {
      if (this->fails_as_input_p == true) {
	print_query_singleend(this,this->fp_nomapping_1,request);
      } else {
	print_header_singleend(this,this->fp_nomapping_1,request,/*translocationp*/false,/*npaths*/0);
	fprintf(this->fp_nomapping_1,"\n");
      }
    }

  } else if (resulttype == SINGLEEND_UNIQ) {
    stage3array = (Stage3end_T *) Result_array(&npaths,&first_absmq,&second_absmq,result);

    if (this->failsonlyp == true) {
      /* Skip */
    } else {
      stage3 = stage3array[0];
      if (Stage3end_circularpos(stage3) > 0) {
	fp = this->fp_unpaired_circular;
      } else {
	fp = this->fp_unpaired_uniq;
      }

      print_header_singleend(this,fp,request,/*translocationp*/false,/*npaths*/1);

      queryseq1 = Request_queryseq1(request);
#if 0
      Stage3end_eval_and_sort(stage3array,/*npaths*/1,this->maxpaths_report,queryseq1);
#endif
      Stage3end_print(fp,stage3,Stage3end_score(stage3),
		      this->chromosome_iit,queryseq1,this->invert_first_p,
		      /*hit5*/(Stage3end_T) NULL,/*hit3*/(Stage3end_T) NULL,
		      /*pairlength*/0,/*pairscore*/0,/*pairtype*/UNPAIRED,
		      Stage3end_mapq_score(stage3));
      fprintf(fp,"\n");
    }

  } else if (resulttype == SINGLEEND_TRANSLOC) {
    stage3array = (Stage3end_T *) Result_array(&npaths,&first_absmq,&second_absmq,result);

    if (this->failsonlyp == true) {
      /* Skip */

    } else if (this->quiet_if_excessive_p && npaths > this->maxpaths_report) {
      print_header_singleend(this,this->fp_unpaired_transloc,request,/*translocationp*/true,npaths);
      fprintf(this->fp_unpaired_transloc,"\n");

    } else {
      print_header_singleend(this,this->fp_unpaired_transloc,request,/*translocationp*/true,npaths);

      queryseq1 = Request_queryseq1(request);
#if 0
      Stage3end_eval_and_sort(stage3array,npaths,this->maxpaths_report,queryseq1);
#endif
      for (pathnum = 1; pathnum <= npaths && pathnum <= this->maxpaths_report; pathnum++) {
	stage3 = stage3array[pathnum-1];
	Stage3end_print(this->fp_unpaired_transloc,stage3,Stage3end_score(stage3),
			this->chromosome_iit,queryseq1,this->invert_first_p,
			/*hit5*/(Stage3end_T) NULL,/*hit3*/(Stage3end_T) NULL,
			/*pairlength*/0,/*pairscore*/0,/*pairtype*/UNPAIRED,
			Stage3end_mapq_score(stage3));
      }
      fprintf(this->fp_unpaired_transloc,"\n");
    }

  } else if (resulttype == SINGLEEND_MULT) {
    stage3array = (Stage3end_T *) Result_array(&npaths,&first_absmq,&second_absmq,result);

    if (this->failsonlyp == true) {
      /* Skip */

    } else if (this->quiet_if_excessive_p && npaths > this->maxpaths_report) {
      print_header_singleend(this,this->fp_unpaired_mult,request,/*translocationp*/false,npaths);
      fprintf(this->fp_unpaired_mult,"\n");

    } else {
      print_header_singleend(this,this->fp_unpaired_mult,request,/*translocationp*/false,npaths);

      queryseq1 = Request_queryseq1(request);
#if 0
      Stage3end_eval_and_sort(stage3array,npaths,this->maxpaths_report,queryseq1);
#endif
      for (pathnum = 1; pathnum <= npaths && pathnum <= this->maxpaths_report; pathnum++) {
	stage3 = stage3array[pathnum-1];
	Stage3end_print(this->fp_unpaired_mult,stage3,Stage3end_score(stage3),
			this->chromosome_iit,queryseq1,this->invert_first_p,
			/*hit5*/(Stage3end_T) NULL,/*hit3*/(Stage3end_T) NULL,
			/*pairlength*/0,/*pairscore*/0,/*pairtype*/UNPAIRED,
			Stage3end_mapq_score(stage3));
      }
      fprintf(this->fp_unpaired_mult,"\n");
    }

  } else {
    if (this->fp_concordant_uniq == NULL) {
      sevenway_open_paired(this);
    }
    Stage3pair_print(result,resulttype,this->chromosome_iit,
		     Request_queryseq1(request),Request_queryseq2(request),
		     this->maxpaths_report,this->quiet_if_excessive_p,
#if 0
		     this->invert_first_p,this->invert_second_p,
#endif
		     this->nofailsp,this->failsonlyp,this->fails_as_input_p,
		     this->fastq_format_p,this->quality_shift,
		     this->fp_nomapping_1,this->fp_nomapping_2,
		     this->fp_unpaired_uniq,this->fp_unpaired_circular,
		     this->fp_unpaired_transloc,this->fp_unpaired_mult,
		     this->fp_halfmapping_uniq,this->fp_halfmapping_circular,
		     this->fp_halfmapping_transloc,this->fp_halfmapping_mult,
		     this->fp_paired_uniq_circular,this->fp_paired_uniq_inv,this->fp_paired_uniq_scr,
		     this->fp_paired_uniq_long,this->fp_paired_mult,
		     this->fp_concordant_uniq,this->fp_concordant_circular,
		     this->fp_concordant_transloc,this->fp_concordant_mult);
  }

  return;
}


static void
print_result_goby (T this, Result_T result, Request_T request) {
  Resulttype_T resulttype;
  Shortread_T queryseq1;
  Stage3end_T *stage3array1, *stage3array2;
  Stage3pair_T *stage3pairarray;
  int npaths1 = 0, npaths2 = 0, first_absmq, second_absmq;
  bool output_alignment = true;

  resulttype = Result_resulttype(result);
  queryseq1 = Request_queryseq1(request);
  switch (resulttype) {
    /* Determine if we are in a TMH situation or some other condition where we */
    /* don't want to output the alignment. */
    case SINGLEEND_NOMAPPING:
    case PAIREDEND_NOMAPPING:
      /* Goby does nothing with no-mapping results. */
      output_alignment = false;
      break;
    case SINGLEEND_MULT:
      /* Check single end Too Many Hits (TMH) */
      stage3array1 = (Stage3end_T *) Result_array(&npaths1,&first_absmq,&second_absmq,result);
      if (npaths1 > this->maxpaths_report) {
        Goby_print_tmh(this->gobywriter,stage3array1[0],queryseq1,npaths1);
        output_alignment = false;
      }
      break;
    case SINGLEEND_UNIQ:
    case SINGLEEND_TRANSLOC:
    case CONCORDANT_UNIQ:
    case CONCORDANT_TRANSLOC:
    case UNPAIRED_UNIQ:
    case UNPAIRED_TRANSLOC:
    case PAIRED_UNIQ:
    case HALFMAPPING_UNIQ:
      /* output alignment but no need to check TMH. */
      break;
    case CONCORDANT_MULT:
    case PAIRED_MULT:
      stage3pairarray = (Stage3pair_T *) Result_array(&npaths1,&first_absmq,&second_absmq,result);
      if (npaths1 > this->maxpaths_report) {
        Goby_print_pair_tmh(this->gobywriter,resulttype,stage3pairarray[0],queryseq1,npaths1);
        output_alignment = false;
      }
      break;
    case UNPAIRED_MULT:
    case HALFMAPPING_TRANSLOC:
    case HALFMAPPING_MULT:
      stage3array1 = (Stage3end_T *) Result_array(&npaths1,&first_absmq,&second_absmq,result);
      stage3array2 = (Stage3end_T *) Result_array2(&npaths2,&first_absmq,&second_absmq,result);
      if (npaths1 >= this->maxpaths_report) {
        Goby_print_tmh(this->gobywriter,stage3array1[0],queryseq1,npaths1);
      }
      if (npaths2 >= this->maxpaths_report) {
        Goby_print_tmh(this->gobywriter,stage3array2[0],queryseq1,npaths2);
      }
      if (npaths1 >= this->maxpaths_report && npaths2 >= this->maxpaths_report) {
        output_alignment = false;
      }
      break;
  }

  if (output_alignment) {
    Goby_start_capture(this->gobywriter);
    print_result_gsnap(this,result,request);
    Goby_finish_capture(this->gobywriter);
  }
  return;
}


void
Outbuffer_print_result (T this, Result_T result, Request_T request
#ifdef MEMUSAGE
			, unsigned int noutput
#endif
			) {
  Shortread_T queryseq1;

  if (this->timingp == true) {
    queryseq1 = Request_queryseq1(request);
    printf("%s\t%.6f\n",Shortread_accession(queryseq1),Result_worker_runtime(result));
  } else if (this->output_sam_p == true) {
    print_result_sam(this,result,request);
  } else if (this->gobywriter != NULL) {
    print_result_goby(this,result,request);
  } else {
    print_result_gsnap(this,result,request);
  }

#ifdef MEMUSAGE
  printf("Memusage of IN: %ld.  Memusage of OUT: %ld.  Entries in outbuffer: %d = %d processed - %u output\n",
	 Mem_usage_in_report(),Mem_usage_out_report(),this->nprocessed - noutput,this->nprocessed,noutput);
#endif

  return;
}

#else

/************************************************************************
 *   Print routines and threads for GMAP
 ************************************************************************/

static void
print_npaths (T this, FILE *fp, int npaths, Diagnostic_T diagnostic,
	      Chrsubset_T chrsubset, bool mergedp, Chimera_T chimera, Failure_T failuretype) {

  if (this->diagnosticp == true) {
    Diagnostic_print(diagnostic);
  }

  if (npaths == 0) {
    fprintf(fp,"Paths (0):");
  } else if (mergedp == true) {
    fprintf(fp,"Paths (1):");
  } else {
    fprintf(fp,"Paths (%d):",npaths);
  }
  Chrsubset_print(chrsubset);
  if (failuretype == NO_FAILURE) {
    if (chimera != NULL) {
      Chimera_print(fp,chimera);
    }
  } else if (failuretype == EMPTY_SEQUENCE) {
    fprintf(fp," *** Empty sequence ***");
  } else if (failuretype == SHORT_SEQUENCE) {
    fprintf(fp," *** Short sequence < index oligo size ***");
  } else if (failuretype == POOR_SEQUENCE) {
    fprintf(fp," *** Poor sequence (use -p flag to change pruning behavior) ***");
  } else if (failuretype == REPETITIVE) {
    fprintf(fp," *** Repetitive sequence (use -p flag to change pruning behavior) ***");
  }
  fprintf(fp,"\n");
  if (npaths == 0) {
    fprintf(fp,"\n");
  }
  return;
}


void
Outbuffer_print_result (T this, Result_T result, Request_T request, Sequence_T headerseq
#ifdef MEMUSAGE
			, unsigned int noutput
#endif
			) {
  FILE *fp;
  char *abbrev;
  Sequence_T queryseq;
  Diagnostic_T diagnostic;
  Stage3_T *stage3array;
  int npaths, pathnum, effective_maxpaths, first_absmq, second_absmq;
  Chimera_T chimera = NULL;
  int chimerapos, chimeraequivpos, chimera_cdna_direction;
  int querylength;
  double donor_prob, acceptor_prob;
  List_T p;
  Gregion_T gregion;
  bool printp, mergedp = false;
#ifdef MEMUSAGE
  char *comma1, *comma2;
#endif

  queryseq = Request_queryseq(request);

  if (this->stage1debug == true) {
    putc('>',stdout);
    Sequence_print_header(stdout,headerseq,this->checksump);

    for (p = Result_gregionlist(result); p != NULL; p = List_next(p)) {
      gregion = (Gregion_T) List_head(p);
      Gregion_print(gregion);
    }
    return;

  } else if (this->diag_debug == true) {
    putc('>',stdout);
    Sequence_print_header(stdout,headerseq,this->checksump);

    Diag_print_segments(Result_diagonals(result),/*queryseq_ptr*/NULL,/*genomicseg_ptr*/NULL);
    return;
  }

  stage3array = Result_array(&npaths,&first_absmq,&second_absmq,result);
  querylength = Sequence_fulllength_given(queryseq);

  chimerapos = chimeraequivpos = -1;
  chimera_cdna_direction = 0;
  donor_prob = acceptor_prob = 0.0;

  /* Translation */
  if (npaths == 0) {
    effective_maxpaths = 0;
    fp = this->fp_nomapping;
    abbrev = ABBREV_NOMAPPING_1;

    if (this->nofailsp == true) {
      printp = false;
      if (this->fails_as_input_p == true) {
	putc('>',fp);
	Sequence_print_header(fp,headerseq,this->checksump);
	Sequence_print(fp,queryseq,/*uppercasep*/false,this->wraplength,/*trimmedp*/false);
      }
    } else {
      printp = true;
    }

    if (Result_failuretype(result) == POOR_SEQUENCE) {
      fprintf(stderr,"Accession %s skipped (poor sequence).  Use -p flag to change pruning behavior\n",Sequence_accession(headerseq));
    } else if (Result_failuretype(result) == REPETITIVE) {
      fprintf(stderr,"Accession %s skipped (repetitive sequence).  Use -p flag to change pruning behavior\n",Sequence_accession(headerseq));
    } else {
      fprintf(stderr,"No paths found for %s\n",Sequence_accession(headerseq));
    }

  } else if ((mergedp = Result_mergedp(result)) == true) {
    if (Stage3_circularpos(stage3array[0]) > 0) {
      fp = this->fp_circular;
      abbrev = ABBREV_UNPAIRED_CIRCULAR;
    } else {
      fp = this->fp_uniq;
      abbrev = ABBREV_UNPAIRED_UNIQ;
    }
    effective_maxpaths = 1;

    if (this->failsonlyp == true) {
      printp = false;
    } else {
      printp = true;

      for (pathnum = 1; pathnum <= /*effective_maxpaths*/1; pathnum++) {
	Stage3_translate(stage3array[pathnum-1],
#ifdef PMAP
			 queryseq,this->diagnosticp,
#endif
			 querylength,this->fulllengthp,
			 this->cds_startpos,this->truncatep,this->strictp,
			 this->maponlyp);
      }
    }

  } else if ((chimera = Result_chimera(result)) != NULL) {
    if (this->chimeras_allowed_p == true) {
      effective_maxpaths = 2;
    } else {
      effective_maxpaths = 0;
    }
    fp = this->fp_transloc;
    abbrev = ABBREV_UNPAIRED_TRANSLOC;

    if (this->failsonlyp == true) {
      printp = false;

#if 0
    } else if (this->quiet_if_excessive_p && npaths > this->maxpaths_report) {
      /* Counting a chimera as a single path */
      printp = true;
#endif

    } else {
      printp = true;

      chimerapos = Chimera_pos(chimera);
      chimeraequivpos = Chimera_equivpos(chimera);
      donor_prob = Chimera_donor_prob(chimera);
      acceptor_prob = Chimera_acceptor_prob(chimera);
      chimera_cdna_direction = Chimera_cdna_direction(chimera);

      Stage3_translate_chimera(stage3array[0],stage3array[1],
#ifdef PMAP
			       queryseq,this->diagnosticp,
#endif
			       querylength,this->fulllengthp,
			       this->cds_startpos,this->truncatep,this->strictp,
			       this->maponlyp);
    }

  } else if (this->maxpaths_report == 0) {
    effective_maxpaths = 1;
    if (npaths > 1) {
      fp = this->fp_mult;
      abbrev = ABBREV_UNPAIRED_MULT;
    } else if (Stage3_circularpos(stage3array[0]) > 0) {
      fp = this->fp_circular;
      abbrev = ABBREV_UNPAIRED_CIRCULAR;
    } else {
      fp = this->fp_uniq;
      abbrev = ABBREV_UNPAIRED_UNIQ;
    }

    if (this->failsonlyp == true) {
      printp = false;
    } else if (this->quiet_if_excessive_p && npaths > this->maxpaths_report) {
      printp = false;
    } else {
      printp = true;

      Stage3_translate(stage3array[0],
#ifdef PMAP
		       queryseq,this->diagnosticp,
#endif
		       querylength,this->fulllengthp,
		       this->cds_startpos,this->truncatep,this->strictp,
		       this->maponlyp);
    }

  } else {
    if (npaths > 1) {
      fp = this->fp_mult;
      abbrev = ABBREV_UNPAIRED_MULT;
    } else if (Stage3_circularpos(stage3array[0]) > 0) {
      fp = this->fp_circular;
      abbrev = ABBREV_UNPAIRED_CIRCULAR;
    } else {
      fp = this->fp_uniq;
      abbrev = ABBREV_UNPAIRED_UNIQ;
    }

    if (npaths < this->maxpaths_report) {
      effective_maxpaths = npaths;
    } else {
      effective_maxpaths = this->maxpaths_report;
    }

    if (this->failsonlyp == true) {
      printp = false;
    } else if (this->quiet_if_excessive_p && npaths > this->maxpaths_report) {
      printp = false;
    } else {
      printp = true;

      for (pathnum = 1; pathnum <= effective_maxpaths; pathnum++) {
	Stage3_translate(stage3array[pathnum-1],
#ifdef PMAP
			 queryseq,this->diagnosticp,
#endif
			 querylength,this->fulllengthp,
			 this->cds_startpos,this->truncatep,this->strictp,
			 this->maponlyp);
      }
    }
  }

  /* Printing */
  if (this->debug_graphic_p == true) {
    printf("q()\n");

  } else if (printp == false) {
    /* No output */

  } else if (this->printtype == SIMPLE || this->printtype == SUMMARY || this->printtype == ALIGNMENT) {
    /* Print header, even if no alignment is found */
    putc('>',fp);
    Sequence_print_header(fp,headerseq,this->checksump);

    diagnostic = Result_diagnostic(result);
    if (npaths == 0) {
      print_npaths(this,fp,0,diagnostic,this->chrsubset,/*mergedp*/false,/*chimera*/NULL,Result_failuretype(result));
    } else {
      print_npaths(this,fp,npaths,diagnostic,this->chrsubset,mergedp,chimera,NO_FAILURE);
      for (pathnum = 1; pathnum <= effective_maxpaths; pathnum++) {
	Stage3_print_pathsummary(fp,stage3array[pathnum-1],pathnum,
				 this->chromosome_iit,this->contig_iit,
				 this->altstrain_iit,queryseq,
				 this->dbversion,/*maxmutations*/1000000,
				 this->diagnosticp,this->maponlyp);
      }
    }

    if (this->printtype != SIMPLE) {
      fprintf(fp,"Alignments:\n");
      for (pathnum = 1; pathnum <= effective_maxpaths; pathnum++) {
	fprintf(fp,"  Alignment for path %d:\n\n",pathnum);
	Stage3_print_alignment(fp,stage3array[pathnum-1],
			       this->genome,this->chromosome_iit,this->printtype,
			       /*continuousp*/false,/*continuous_by_exon_p*/false,
			       this->diagnosticp,/*flipgenomep*/true,
			       this->invertmode,this->nointronlenp,
			       this->wraplength);
      }
    }

    if (this->map_iit != NULL) {
      fprintf(fp,"Maps:\n");
      for (pathnum = 1; pathnum <= effective_maxpaths; pathnum++) {
	Stage3_print_map(fp,stage3array[pathnum-1],this->map_iit,this->map_divint_crosstable,
			 this->chromosome_iit,pathnum,this->map_exons_p,this->map_bothstrands_p,
			 this->nflanking,this->print_comment_p);
      }
    }

  } else if (this->printtype == COMPRESSED) {
    for (pathnum = 1; pathnum <= effective_maxpaths; pathnum++) {
      Stage3_print_compressed(fp,stage3array[pathnum-1],queryseq,this->chromosome_iit,
			      this->dbversion,this->usersegment,pathnum,npaths,
			      this->checksump,chimerapos,chimeraequivpos,
			      donor_prob,acceptor_prob,chimera_cdna_direction);
    }

  } else if (this->printtype == CONTINUOUS) {
    putc('>',fp);
    Sequence_print_header(fp,headerseq,this->checksump);
    if (npaths == 0) {
      fprintf(fp,"\n\n\n");
    } else {
      Stage3_print_alignment(fp,stage3array[0],this->genome,this->chromosome_iit,this->printtype,
			     /*continuousp*/true,/*continuous_by_exon_p*/false,
			     this->diagnosticp,/*flipgenomep*/true,
			     this->invertmode,this->nointronlenp,
			     this->wraplength);
    }

  } else if (this->printtype == CONTINUOUS_BY_EXON) {
    diagnostic = Result_diagnostic(result);

    putc('>',fp);
    Sequence_print_header(fp,headerseq,this->checksump);
    print_npaths(this,fp,npaths,diagnostic,this->chrsubset,mergedp,chimera,NO_FAILURE);
    if (npaths == 0) {
      fprintf(fp,"\n\n\n");
    } else {
      Stage3_print_pathsummary(fp,stage3array[0],/*pathnum*/1,
			       this->chromosome_iit,this->contig_iit,
			       this->altstrain_iit,queryseq,
			       this->dbversion,/*maxmutations*/1000000,
			       this->diagnosticp,this->maponlyp);
      fprintf(fp,"Alignments:\n");
      fprintf(fp,"  Alignment for path %d:\n\n",/*pathnum*/1);
      Stage3_print_alignment(fp,stage3array[0],this->genome,this->chromosome_iit,this->printtype,
			     /*continuousp*/false,/*continuous_by_exon_p*/true,
			     this->diagnosticp,/*flipgenomep*/true,
			     this->invertmode,this->nointronlenp,
			     this->wraplength);
    }

  } else if (this->printtype == EXONS_CDNA) {
    putc('>',fp);
    Sequence_print_header(fp,headerseq,this->checksump);
    for (pathnum = 1; pathnum <= effective_maxpaths; pathnum++) {
      fprintf(fp,"<path %d>\n",pathnum);
      Pair_print_exons(fp,Stage3_pairarray(stage3array[0]),Stage3_npairs(stage3array[0]),
		       this->wraplength,this->ngap,/*cdna*/true);
      fprintf(fp,"</path>\n");
    }

  } else if (this->printtype == EXONS_GENOMIC) {
    putc('>',fp);
    Sequence_print_header(fp,headerseq,this->checksump);
    for (pathnum = 1; pathnum <= effective_maxpaths; pathnum++) {
      fprintf(fp,"<path %d>\n",pathnum);
      Pair_print_exons(fp,Stage3_pairarray(stage3array[0]),Stage3_npairs(stage3array[0]),
		       this->wraplength,this->ngap,/*cdna*/false);
      fprintf(fp,"</path>\n");
    }

  } else if (this->printtype == CDNA) {
    for (pathnum = 1; pathnum <= effective_maxpaths; pathnum++) {
      putc('>',fp);
      Sequence_print_header(fp,headerseq,this->checksump);
      Stage3_print_cdna(fp,stage3array[pathnum-1],this->wraplength);
    }

  } else if (this->printtype == PROTEIN_GENOMIC) {
    for (pathnum = 1; pathnum <= effective_maxpaths; pathnum++) {
      putc('>',fp);
      Sequence_print_header(fp,headerseq,this->checksump);
      Stage3_print_protein_genomic(fp,stage3array[pathnum-1],this->wraplength);
    }

  } else if (this->printtype == PSL_NT) {
    for (pathnum = 1; pathnum <= effective_maxpaths; pathnum++) {
      Stage3_print_pslformat_nt(fp,stage3array[pathnum-1],
				this->chromosome_iit,this->usersegment,queryseq);
    }

#ifdef PMAP
  } else if (this->printtype == PSL_PRO) {
    for (pathnum = 1; pathnum <= effective_maxpaths; pathnum++) {
      Stage3_print_pslformat_pro(fp,stage3array[pathnum-1],
				 this->chromosome_iit,this->usersegment,queryseq,this->strictp);
    }
#endif

  } else if (this->printtype == GFF3_GENE || this->printtype == GFF3_MATCH_CDNA ||
	     this->printtype == GFF3_MATCH_EST) {
    for (pathnum = 1; pathnum <= effective_maxpaths; pathnum++) {
      Stage3_print_gff3(fp,stage3array[pathnum-1],pathnum,
			this->chromosome_iit,this->usersegment,queryseq,querylength,this->printtype,
			/*sourcename*/this->usersegment ? this->user_genomicseg : this->dbversion);
    }

#ifndef PMAP
  } else if (this->printtype == SAM) {
    if (npaths == 0) {
      Pair_print_sam_nomapping(fp,abbrev,/*acc1*/Sequence_accession(headerseq),/*acc2*/NULL,
			       Sequence_fullpointer(queryseq),Sequence_quality_string(queryseq),
			       Sequence_fulllength(queryseq),this->quality_shift,
			       Sequence_firstp(queryseq),this->sam_paired_p,this->sam_read_group_id);

    } else if (this->quiet_if_excessive_p && npaths > this->maxpaths_report) {
      Pair_print_sam_nomapping(fp,abbrev,/*acc1*/Sequence_accession(headerseq),/*acc2*/NULL,
			       Sequence_fullpointer(queryseq),Sequence_quality_string(queryseq),
			       Sequence_fulllength(queryseq),this->quality_shift,
			       Sequence_firstp(queryseq),this->sam_paired_p,this->sam_read_group_id);

    } else if (mergedp == true) {
      Stage3_print_sam(fp,abbrev,stage3array[0],/*pathnum*/1,/*npaths*/1,
		       Stage3_absmq_score(stage3array[0]),first_absmq,second_absmq,
		       Stage3_mapq_score(stage3array[0]),
		       this->chromosome_iit,this->usersegment,queryseq,
		       /*chimera_part*/0,/*chimera*/NULL,this->quality_shift,this->sam_paired_p,
		       this->sam_read_group_id);

    } else if (chimera != NULL) {
      Stage3_print_sam(fp,abbrev,stage3array[0],/*pathnum*/1,npaths,
		       Stage3_absmq_score(stage3array[0]),first_absmq,second_absmq,
		       Stage3_mapq_score(stage3array[0]),
		       this->chromosome_iit,this->usersegment,queryseq,
		       /*chimera_part*/-1,chimera,this->quality_shift,this->sam_paired_p,
		       this->sam_read_group_id);
      Stage3_print_sam(fp,abbrev,stage3array[1],/*pathnum*/1,npaths,
		       Stage3_absmq_score(stage3array[0]),first_absmq,second_absmq,
		       Stage3_mapq_score(stage3array[0]),
		       this->chromosome_iit,this->usersegment,queryseq,
		       /*chimera_part*/+1,chimera,this->quality_shift,this->sam_paired_p,
		       this->sam_read_group_id);

    } else {
      for (pathnum = 1; pathnum <= effective_maxpaths; pathnum++) {
	Stage3_print_sam(fp,abbrev,stage3array[pathnum-1],pathnum,npaths,
			 Stage3_absmq_score(stage3array[pathnum-1]),first_absmq,second_absmq,
			 Stage3_mapq_score(stage3array[pathnum-1]),
			 this->chromosome_iit,this->usersegment,queryseq,
			 /*chimera_part*/0,/*chimera*/NULL,this->quality_shift,this->sam_paired_p,
			 this->sam_read_group_id);
      }
    }
#endif

  } else if (this->printtype == COORDS) {
    for (pathnum = 1; pathnum <= effective_maxpaths; pathnum++) {
      fprintf(fp,">");
      Sequence_print_header(fp,headerseq,this->checksump);
      Stage3_print_coordinates(fp,stage3array[pathnum-1],this->chromosome_iit,this->invertmode);
    }

  } else if (this->printtype == SPLICESITES) {
    /* Print only best path */
    if (npaths > 0) {
      Stage3_print_splicesites(fp,stage3array[0],this->chromosome_iit,queryseq);
    }

  } else if (this->printtype == INTRONS) {
    /* Print only best path */
    if (npaths > 0) {
      Stage3_print_introns(fp,stage3array[0],this->chromosome_iit,queryseq);
    }

  } else if (this->printtype == MAP_RANGES) {
    for (pathnum = 1; pathnum <= effective_maxpaths; pathnum++) {
      Stage3_print_iit_map(fp,stage3array[pathnum-1],this->chromosome_iit,queryseq);
    }
      
  } else if (this->printtype == MAP_EXONS) {
    for (pathnum = 1; pathnum <= effective_maxpaths; pathnum++) {
      Stage3_print_iit_exon_map(fp,stage3array[pathnum-1],this->chromosome_iit,queryseq);
    }

  } else {
    fprintf(stderr,"Unexpected printtype %d\n",this->printtype);
    abort();

  }

#ifdef MEMUSAGE
  comma1 = Genomicpos_commafmt(Mem_usage_report());
  comma2 = Genomicpos_commafmt(Mem_max_usage_report());
  printf("Memusage: %s.  Peak: %s.\n",comma1,comma2);
  FREE(comma2);
  FREE(comma1);
#endif

  return;
}

#endif


void *
Outbuffer_thread_anyorder (void *data) {
  T this = (T) data;
  unsigned int output_buffer_size = this->output_buffer_size;
  unsigned int noutput = 0;
  Result_T result;
  Request_T request;
  
#ifdef MEMUSAGE
  Mem_usage_set_threadname("outbuffer");
#endif

  while (noutput < this->ntotal) {

#ifdef HAVE_PTHREAD
    pthread_mutex_lock(&this->lock);
    while (this->head == NULL && noutput < this->ntotal) {
      debug(fprintf(stderr,"__outbuffer_thread_anyorder waiting for result_avail_p\n"));
      pthread_cond_wait(&this->result_avail_p,&this->lock);
    }
    debug(fprintf(stderr,"__outbuffer_thread_anyorder woke up\n"));
#endif

    if (this->head == NULL) {
      /* False wake up */
#ifdef HAVE_PTHREAD
      pthread_mutex_unlock(&this->lock);
#endif

    } else {
      this->head = RRlist_pop(this->head,&request,&result);
      debug1(RRlist_dump(this->head,this->tail));

#ifdef HAVE_PTHREAD
      /* Let worker threads put results while we print */
      pthread_mutex_unlock(&this->lock);
#endif
#ifdef MEMUSAGE
      Outbuffer_print_result(this,result,request,
#ifndef GSNAP
			     Request_queryseq(request),
#endif
			     noutput+1);
#else
      Outbuffer_print_result(this,result,request
#ifndef GSNAP
			     ,Request_queryseq(request)
#endif
			     );
#endif
      Result_free(&result);
      Request_free(&request);
      noutput++;

      if (this->head && this->nprocessed - noutput > output_buffer_size) {
	/* Clear out backlog */
#ifdef HAVE_PTHREAD
	pthread_mutex_lock(&this->lock);
#endif
	while (this->head && this->nprocessed - noutput > output_buffer_size) {
	  this->head = RRlist_pop(this->head,&request,&result);
	  debug1(RRlist_dump(this->head,this->tail));

#ifdef MEMUSAGE
	  Outbuffer_print_result(this,result,request,
#ifndef GSNAP
				 Request_queryseq(request),
#endif
				 noutput+1);
#else
	  Outbuffer_print_result(this,result,request
#ifndef GSNAP
				 ,Request_queryseq(request)
#endif
				 );
#endif
	  Result_free(&result);
	  Request_free(&request);
	  noutput++;
	}

#ifdef HAVE_PTHREAD
	pthread_mutex_unlock(&this->lock);
#endif
      }

    }
  }

  assert(this->head == NULL);

  return (void *) NULL;
}



void *
Outbuffer_thread_ordered (void *data) {
  T this = (T) data;
  unsigned int output_buffer_size = this->output_buffer_size;
  unsigned int noutput = 0, nqueued = 0;
  Result_T result;
  Request_T request;
  RRlist_T queue = NULL;
  int id;

#ifdef MEMUSAGE
  Mem_usage_set_threadname("outbuffer");
#endif

  while (noutput < this->ntotal) {
#ifdef HAVE_PTHREAD
    pthread_mutex_lock(&this->lock);
    while (this->head == NULL && noutput < this->ntotal) {
      pthread_cond_wait(&this->result_avail_p,&this->lock);
    }
    debug(fprintf(stderr,"__outbuffer_thread_ordered woke up\n"));
#endif

    if (this->head == NULL) {
#ifdef HAVE_PTHREAD
      /* False wake up */
      pthread_mutex_unlock(&this->lock);
#endif

    } else {
      this->head = RRlist_pop(this->head,&request,&result);

#ifdef HAVE_PTHREAD
      /* Allow workers access to the queue */
      pthread_mutex_unlock(&this->lock);
#endif
      if ((id = Result_id(result)) != (int) noutput) {
	/* Store in queue */
	queue = RRlist_insert(queue,id,request,result);
	nqueued++;
      } else {
#ifdef MEMUSAGE
	Outbuffer_print_result(this,result,request,
#ifndef GSNAP
			       Request_queryseq(request),
#endif
			       noutput+1);
#else
	Outbuffer_print_result(this,result,request
#ifndef GSNAP
			       ,Request_queryseq(request)
#endif
			       );
#endif
	Result_free(&result);
	Request_free(&request);
	noutput++;
	
	/* Print out rest of stored queue */
	while (queue != NULL && queue->id == (int) noutput) {
	  queue = RRlist_pop_id(queue,&id,&request,&result);
	  nqueued--;
#ifdef MEMUSAGE
	  Outbuffer_print_result(this,result,request,
#ifndef GSNAP
				 Request_queryseq(request),
#endif
				 noutput+1);
#else
	  Outbuffer_print_result(this,result,request
#ifndef GSNAP
				 ,Request_queryseq(request)
#endif
				 );
#endif
	  Result_free(&result);
	  Request_free(&request);
	  noutput++;
	}
      }

      if (this->head && this->nprocessed - nqueued - noutput > output_buffer_size) {
	/* Clear out backlog */
#ifdef HAVE_PTHREAD
	pthread_mutex_lock(&this->lock);
#endif

	while (this->head && this->nprocessed - nqueued - noutput > output_buffer_size) {
	  this->head = RRlist_pop(this->head,&request,&result);
	  if ((id = Result_id(result)) != (int) noutput) {
	    /* Store in queue */
	    queue = RRlist_insert(queue,id,request,result);
	    nqueued++;
	  } else {
#ifdef MEMUSAGE
	    Outbuffer_print_result(this,result,request,
#ifndef GSNAP
				   Request_queryseq(request),
#endif				   
				   noutput+1);
#else
	    Outbuffer_print_result(this,result,request
#ifndef GSNAP
				   ,Request_queryseq(request)
#endif
				   );
#endif
	    Result_free(&result);
	    Request_free(&request);
	    noutput++;
	
	    /* Print out rest of stored queue */
	    while (queue != NULL && queue->id == (int) noutput) {
	      queue = RRlist_pop_id(queue,&id,&request,&result);
	      nqueued--;
#ifdef MEMUSAGE
	      Outbuffer_print_result(this,result,request,
#ifndef GSNAP
				     Request_queryseq(request),
#endif
				     noutput+1);
#else
	      Outbuffer_print_result(this,result,request
#ifndef GSNAP
				     ,Request_queryseq(request)
#endif
				     );
#endif
	      Result_free(&result);
	      Request_free(&request);
	      noutput++;
	    }
	  }
	}

#ifdef HAVE_PTHREAD
	pthread_mutex_unlock(&this->lock);
#endif
      }

    }
  }

  assert(queue == NULL);

  return (void *) NULL;
}


