static char rcsid[] = "$Id: outbuffer.c 53170 2011-11-28 02:56:47Z twu $";
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

  IIT_T chromosome_iit;

  char *sevenway_root;

#ifdef GSNAP
  bool sam_headers_p;
  char *sam_read_group_id;
  char *sam_read_group_name;
  char *sam_read_group_library;
  char *sam_read_group_platform;
  int quality_shift;
#elif defined PMAP

#else
  bool sam_headers_p;
  bool sam_paired_p;
  char *sam_read_group_id;
  char *sam_read_group_name;
  char *sam_read_group_library;
  char *sam_read_group_platform;
  int quality_shift;
#endif

#ifdef GSNAP

  FILE *fp_nomapping_1;
  FILE *fp_nomapping_2;
  FILE *fp_halfmapping_uniq;
  FILE *fp_halfmapping_transloc;
  FILE *fp_halfmapping_mult;
  FILE *fp_unpaired_uniq;
  FILE *fp_unpaired_transloc;
  FILE *fp_unpaired_mult;
  FILE *fp_paired_uniq_inv;
  FILE *fp_paired_uniq_scr;
  FILE *fp_paired_uniq_long;
  FILE *fp_paired_mult;
  FILE *fp_concordant_uniq;
  FILE *fp_concordant_transloc;
  FILE *fp_concordant_mult;

  bool timingp;
  bool output_sam_p;
  Gobywriter_T gobywriter;

  bool fastq_format_p;
  bool clip_overlap_p;
  bool merge_samechr_p;

  bool invert_first_p;
  bool invert_second_p;
  Genomicpos_T pairmax;

#else

  FILE *fp_nomapping;
  FILE *fp_uniq;
  FILE *fp_transloc;
  FILE *fp_mult;

  bool chimerap;

  char *user_genomicseg;
  Sequence_T usersegment;

  char *dbversion;
  Chrsubset_T chrsubset;
  IIT_T contig_iit;
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

  int maxpaths;
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

  if (this->fails_as_input_p == true) {
    filename = (char *) CALLOC(strlen(this->sevenway_root)+strlen(".nomapping.fq")+1,sizeof(char));
    sprintf(filename,"%s.nomapping.fq",this->sevenway_root);
    if ((this->fp_nomapping_1 = fopen(filename,"w")) == NULL) {
      fprintf(stderr,"Cannot open file %s for writing\n",filename);
      exit(9);
    }
    FREE(filename);

  } else {
    filename = (char *) CALLOC(strlen(this->sevenway_root)+strlen(".nomapping")+1,sizeof(char));
    sprintf(filename,"%s.nomapping",this->sevenway_root);
    if ((this->fp_nomapping_1 = fopen(filename,"w")) == NULL) {
      fprintf(stderr,"Cannot open file %s for writing\n",filename);
      exit(9);
    }
    FREE(filename);

    if (this->output_sam_p == true && this->sam_headers_p == true) {
      IIT_dump_sam(this->fp_nomapping_1,this->chromosome_iit,this->sam_read_group_id,this->sam_read_group_name,
		   this->sam_read_group_library,this->sam_read_group_platform);
    }
  }

  filename = (char *) CALLOC(strlen(this->sevenway_root)+strlen(".unpaired_uniq")+1,sizeof(char));
  sprintf(filename,"%s.unpaired_uniq",this->sevenway_root);
  if ((this->fp_unpaired_uniq = fopen(filename,"w")) == NULL) {
    fprintf(stderr,"Cannot open file %s for writing\n",filename);
    exit(9);
  }
  FREE(filename);

  filename = (char *) CALLOC(strlen(this->sevenway_root)+strlen(".unpaired_transloc")+1,sizeof(char));
  sprintf(filename,"%s.unpaired_transloc",this->sevenway_root);
  if ((this->fp_unpaired_transloc = fopen(filename,"w")) == NULL) {
    fprintf(stderr,"Cannot open file %s for writing\n",filename);
    exit(9);
  }
  FREE(filename);

  filename = (char *) CALLOC(strlen(this->sevenway_root)+strlen(".unpaired_mult")+1,sizeof(char));
  sprintf(filename,"%s.unpaired_mult",this->sevenway_root);
  if ((this->fp_unpaired_mult = fopen(filename,"w")) == NULL) {
    fprintf(stderr,"Cannot open file %s for writing\n",filename);
    exit(9);
  }
  FREE(filename);

  if (this->output_sam_p == true && this->sam_headers_p == true) {
    IIT_dump_sam(this->fp_unpaired_uniq,this->chromosome_iit,this->sam_read_group_id,this->sam_read_group_name,
		 this->sam_read_group_library,this->sam_read_group_platform);
    IIT_dump_sam(this->fp_unpaired_transloc,this->chromosome_iit,this->sam_read_group_id,this->sam_read_group_name,
		 this->sam_read_group_library,this->sam_read_group_platform);
    IIT_dump_sam(this->fp_unpaired_mult,this->chromosome_iit,this->sam_read_group_id,this->sam_read_group_name,
		 this->sam_read_group_library,this->sam_read_group_platform);
  }

  return;
}


static void
sevenway_open_paired (T this) {
  char *filename;

  if (this->fails_as_input_p == true) {
    filename = (char *) CALLOC(strlen(this->sevenway_root)+strlen(".nomapping.1.fq")+1,sizeof(char));
    sprintf(filename,"%s.nomapping.1.fq",this->sevenway_root);
    if ((this->fp_nomapping_1 = fopen(filename,"w")) == NULL) {
      fprintf(stderr,"Cannot open file %s for writing\n",filename);
      exit(9);
    }
    FREE(filename);

    filename = (char *) CALLOC(strlen(this->sevenway_root)+strlen(".nomapping.2.fq")+1,sizeof(char));
    sprintf(filename,"%s.nomapping.2.fq",this->sevenway_root);
    if ((this->fp_nomapping_2 = fopen(filename,"w")) == NULL) {
      fprintf(stderr,"Cannot open file %s for writing\n",filename);
      exit(9);
    }
    FREE(filename);

  } else {
    filename = (char *) CALLOC(strlen(this->sevenway_root)+strlen(".nomapping")+1,sizeof(char));
    sprintf(filename,"%s.nomapping",this->sevenway_root);
    if ((this->fp_nomapping_1 = fopen(filename,"w")) == NULL) {
      fprintf(stderr,"Cannot open file %s for writing\n",filename);
      exit(9);
    }
    FREE(filename);

    if (this->output_sam_p == true && this->sam_headers_p == true) {
      IIT_dump_sam(this->fp_nomapping_1,this->chromosome_iit,this->sam_read_group_id,this->sam_read_group_name,
		   this->sam_read_group_library,this->sam_read_group_platform);
    }
  }

  filename = (char *) CALLOC(strlen(this->sevenway_root)+strlen(".halfmapping_uniq")+1,sizeof(char));
  sprintf(filename,"%s.halfmapping_uniq",this->sevenway_root);
  if ((this->fp_halfmapping_uniq = fopen(filename,"w")) == NULL) {
    fprintf(stderr,"Cannot open file %s for writing\n",filename);
    exit(9);
  }
  FREE(filename);

  filename = (char *) CALLOC(strlen(this->sevenway_root)+strlen(".halfmapping_transloc")+1,sizeof(char));
  sprintf(filename,"%s.halfmapping_transloc",this->sevenway_root);
  if ((this->fp_halfmapping_transloc = fopen(filename,"w")) == NULL) {
    fprintf(stderr,"Cannot open file %s for writing\n",filename);
    exit(9);
  }
  FREE(filename);

  filename = (char *) CALLOC(strlen(this->sevenway_root)+strlen(".halfmapping_mult")+1,sizeof(char));
  sprintf(filename,"%s.halfmapping_mult",this->sevenway_root);
  if ((this->fp_halfmapping_mult = fopen(filename,"w")) == NULL) {
    fprintf(stderr,"Cannot open file %s for writing\n",filename);
    exit(9);
  }
  FREE(filename);

  filename = (char *) CALLOC(strlen(this->sevenway_root)+strlen(".paired_uniq_inv")+1,sizeof(char));
  sprintf(filename,"%s.paired_uniq_inv",this->sevenway_root);
  if ((this->fp_paired_uniq_inv = fopen(filename,"w")) == NULL) {
    fprintf(stderr,"Cannot open file %s for writing\n",filename);
    exit(9);
  }
  FREE(filename);

  filename = (char *) CALLOC(strlen(this->sevenway_root)+strlen(".paired_uniq_scr")+1,sizeof(char));
  sprintf(filename,"%s.paired_uniq_scr",this->sevenway_root);
  if ((this->fp_paired_uniq_scr = fopen(filename,"w")) == NULL) {
    fprintf(stderr,"Cannot open file %s for writing\n",filename);
    exit(9);
  }
  FREE(filename);

  filename = (char *) CALLOC(strlen(this->sevenway_root)+strlen(".paired_uniq_long")+1,sizeof(char));
  sprintf(filename,"%s.paired_uniq_long",this->sevenway_root);
  if ((this->fp_paired_uniq_long = fopen(filename,"w")) == NULL) {
    fprintf(stderr,"Cannot open file %s for writing\n",filename);
    exit(9);
  }
  FREE(filename);

  filename = (char *) CALLOC(strlen(this->sevenway_root)+strlen(".paired_mult")+1,sizeof(char));
  sprintf(filename,"%s.paired_mult",this->sevenway_root);
  if ((this->fp_paired_mult = fopen(filename,"w")) == NULL) {
    fprintf(stderr,"Cannot open file %s for writing\n",filename);
    exit(9);
  }
  FREE(filename);

  filename = (char *) CALLOC(strlen(this->sevenway_root)+strlen(".concordant_uniq")+1,sizeof(char));
  sprintf(filename,"%s.concordant_uniq",this->sevenway_root);
  if ((this->fp_concordant_uniq = fopen(filename,"w")) == NULL) {
    fprintf(stderr,"Cannot open file %s for writing\n",filename);
    exit(9);
  }
  FREE(filename);

  filename = (char *) CALLOC(strlen(this->sevenway_root)+strlen(".concordant_transloc")+1,sizeof(char));
  sprintf(filename,"%s.concordant_transloc",this->sevenway_root);
  if ((this->fp_concordant_transloc = fopen(filename,"w")) == NULL) {
    fprintf(stderr,"Cannot open file %s for writing\n",filename);
    exit(9);
  }
  FREE(filename);

  filename = (char *) CALLOC(strlen(this->sevenway_root)+strlen(".concordant_mult")+1,sizeof(char));
  sprintf(filename,"%s.concordant_mult",this->sevenway_root);
  if ((this->fp_concordant_mult = fopen(filename,"w")) == NULL) {
    fprintf(stderr,"Cannot open file %s for writing\n",filename);
    exit(9);
  }
  FREE(filename);

  if (this->output_sam_p == true && this->sam_headers_p == true) {
    IIT_dump_sam(this->fp_halfmapping_uniq,this->chromosome_iit,this->sam_read_group_id,this->sam_read_group_name,
		 this->sam_read_group_library,this->sam_read_group_platform);
    IIT_dump_sam(this->fp_halfmapping_transloc,this->chromosome_iit,this->sam_read_group_id,this->sam_read_group_name,
		 this->sam_read_group_library,this->sam_read_group_platform);
    IIT_dump_sam(this->fp_halfmapping_mult,this->chromosome_iit,this->sam_read_group_id,this->sam_read_group_name,
		 this->sam_read_group_library,this->sam_read_group_platform);
    IIT_dump_sam(this->fp_paired_uniq_inv,this->chromosome_iit,this->sam_read_group_id,this->sam_read_group_name,
		 this->sam_read_group_library,this->sam_read_group_platform);
    IIT_dump_sam(this->fp_paired_uniq_scr,this->chromosome_iit,this->sam_read_group_id,this->sam_read_group_name,
		 this->sam_read_group_library,this->sam_read_group_platform);
    IIT_dump_sam(this->fp_paired_uniq_long,this->chromosome_iit,this->sam_read_group_id,this->sam_read_group_name,
		 this->sam_read_group_library,this->sam_read_group_platform);
    IIT_dump_sam(this->fp_paired_mult,this->chromosome_iit,this->sam_read_group_id,this->sam_read_group_name,
		 this->sam_read_group_library,this->sam_read_group_platform);
    IIT_dump_sam(this->fp_concordant_uniq,this->chromosome_iit,this->sam_read_group_id,this->sam_read_group_name,
		 this->sam_read_group_library,this->sam_read_group_platform);
    IIT_dump_sam(this->fp_concordant_transloc,this->chromosome_iit,this->sam_read_group_id,this->sam_read_group_name,
		 this->sam_read_group_library,this->sam_read_group_platform);
    IIT_dump_sam(this->fp_concordant_mult,this->chromosome_iit,this->sam_read_group_id,this->sam_read_group_name,
		 this->sam_read_group_library,this->sam_read_group_platform);
  }

  return;
}

static void
sevenway_close (T this) {
  fclose(this->fp_unpaired_uniq);
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
    fclose(this->fp_halfmapping_transloc);
    fclose(this->fp_halfmapping_mult);
    fclose(this->fp_paired_uniq_long);
    fclose(this->fp_paired_uniq_scr);
    fclose(this->fp_paired_uniq_inv);
    fclose(this->fp_paired_mult);
    fclose(this->fp_concordant_uniq);
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


/* Taken from IIT_dump_sam */
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

  if (this->fails_as_input_p == true) {
    filename = (char *) CALLOC(strlen(this->sevenway_root)+strlen(".nomapping.fa")+1,sizeof(char));
    sprintf(filename,"%s.nomapping.fa",this->sevenway_root);
    if ((this->fp_nomapping = fopen(filename,"w")) == NULL) {
      fprintf(stderr,"Cannot open file %s for writing\n",filename);
      exit(9);
    }
    FREE(filename);

  } else {
    filename = (char *) CALLOC(strlen(this->sevenway_root)+strlen(".nomapping")+1,sizeof(char));
    sprintf(filename,"%s.nomapping",this->sevenway_root);
    if ((this->fp_nomapping = fopen(filename,"w")) == NULL) {
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
	IIT_dump_sam(this->fp_nomapping,this->chromosome_iit,this->sam_read_group_id,this->sam_read_group_name,
		     this->sam_read_group_library,this->sam_read_group_platform);
      }
#endif
    }
  }

  filename = (char *) CALLOC(strlen(this->sevenway_root)+strlen(".uniq")+1,sizeof(char));
  sprintf(filename,"%s.uniq",this->sevenway_root);
  if ((this->fp_uniq = fopen(filename,"w")) == NULL) {
    fprintf(stderr,"Cannot open file %s for writing\n",filename);
    exit(9);
  }
  FREE(filename);

  if (this->chimerap == true) {
    filename = (char *) CALLOC(strlen(this->sevenway_root)+strlen(".transloc")+1,sizeof(char));
    sprintf(filename,"%s.transloc",this->sevenway_root);
    if ((this->fp_transloc = fopen(filename,"w")) == NULL) {
      fprintf(stderr,"Cannot open file %s for writing\n",filename);
      exit(9);
    }
    FREE(filename);
  }

  filename = (char *) CALLOC(strlen(this->sevenway_root)+strlen(".mult")+1,sizeof(char));
  sprintf(filename,"%s.mult",this->sevenway_root);
  if ((this->fp_mult = fopen(filename,"w")) == NULL) {
    fprintf(stderr,"Cannot open file %s for writing\n",filename);
    exit(9);
  }
  FREE(filename);

  if (this->printtype == GFF3_GENE || this->printtype == GFF3_MATCH_CDNA || this->printtype == GFF3_MATCH_EST) {
    print_gff_header(this->fp_uniq,argc,argv,optind);
    print_gff_header(this->fp_mult,argc,argv,optind);

#ifndef PMAP
  } else if (this->printtype == SAM && this->sam_headers_p == true) {
    if (this->usersegment != NULL) {
      dump_sam_usersegment(this->fp_uniq,this->usersegment,
			       this->sam_read_group_id,this->sam_read_group_name,
			       this->sam_read_group_library,this->sam_read_group_platform);
      dump_sam_usersegment(this->fp_mult,this->usersegment,
			       this->sam_read_group_id,this->sam_read_group_name,
			       this->sam_read_group_library,this->sam_read_group_platform);
    } else {
      IIT_dump_sam(this->fp_uniq,this->chromosome_iit,this->sam_read_group_id,this->sam_read_group_name,
		   this->sam_read_group_library,this->sam_read_group_platform);
      IIT_dump_sam(this->fp_mult,this->chromosome_iit,this->sam_read_group_id,this->sam_read_group_name,
		   this->sam_read_group_library,this->sam_read_group_platform);
    }
#endif
  }

  return;
}

static void
sevenway_close (T this) {
  fclose(this->fp_mult);
  fclose(this->fp_uniq);
  if (this->chimerap == true) {
    fclose(this->fp_transloc);
  }
  fclose(this->fp_nomapping);
  return;
}

#endif



#ifdef GSNAP

T
Outbuffer_new (unsigned int output_buffer_size, unsigned int nread, char *sevenway_root, IIT_T chromosome_iit,
	       bool timingp, bool output_sam_p, bool sam_headers_p, char *sam_read_group_id, char *sam_read_group_name,
	       char *sam_read_group_library, char *sam_read_group_platform,
	       Gobywriter_T gobywriter, bool nofailsp, bool failsonlyp, bool fails_as_input_p,
	       bool fastq_format_p, bool clip_overlap_p, bool merge_samechr_p,
	       int maxpaths, bool quiet_if_excessive_p, int quality_shift,
	       bool invert_first_p, bool invert_second_p, Genomicpos_T pairmax) {
  T new = (T) MALLOC(sizeof(*new));

  new->chromosome_iit = chromosome_iit;

  new->fp_nomapping_1 = NULL;
  new->fp_nomapping_2 = NULL;
  new->fp_halfmapping_uniq = NULL;
  new->fp_halfmapping_transloc = NULL;
  new->fp_halfmapping_mult = NULL;
  new->fp_unpaired_uniq = NULL;
  new->fp_unpaired_transloc = NULL;
  new->fp_unpaired_mult = NULL;
  new->fp_paired_uniq_inv = NULL;
  new->fp_paired_uniq_scr = NULL;
  new->fp_paired_uniq_long = NULL;
  new->fp_paired_mult = NULL;
  new->fp_concordant_uniq = NULL;
  new->fp_concordant_transloc = NULL;
  new->fp_concordant_mult = NULL;
  
  new->sevenway_root = sevenway_root;

  new->timingp = timingp;
  new->output_sam_p = output_sam_p;
  new->sam_headers_p = sam_headers_p;
  new->sam_read_group_id = sam_read_group_id;
  new->sam_read_group_name = sam_read_group_name;
  new->sam_read_group_library = sam_read_group_library;
  new->sam_read_group_platform = sam_read_group_platform;

  new->gobywriter = gobywriter;

  new->nofailsp = nofailsp;
  new->failsonlyp = failsonlyp;
  new->fails_as_input_p = fails_as_input_p;
  new->fastq_format_p = fastq_format_p;
  new->clip_overlap_p = clip_overlap_p;
  new->merge_samechr_p = merge_samechr_p;

  new->maxpaths = maxpaths;
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
  if (sevenway_root != NULL) {
    sevenway_open_single(new);

  } else {
    new->fp_nomapping_1 = stdout;
    new->fp_nomapping_2 = stdout;
    new->fp_halfmapping_uniq = stdout;
    new->fp_halfmapping_transloc = stdout;
    new->fp_halfmapping_mult = stdout;
    new->fp_unpaired_uniq = stdout;
    new->fp_unpaired_transloc = stdout;
    new->fp_unpaired_mult = stdout;
    new->fp_paired_uniq_inv = stdout;
    new->fp_paired_uniq_scr = stdout;
    new->fp_paired_uniq_long = stdout;
    new->fp_paired_mult = stdout;
    new->fp_concordant_uniq = stdout;
    new->fp_concordant_transloc = stdout;
    new->fp_concordant_mult = stdout;

    if (output_sam_p == true && sam_headers_p == true) {
      if (fails_as_input_p == true) {
	/* Don't print chromosomes */
      } else {
	IIT_dump_sam(stdout,chromosome_iit,sam_read_group_id,sam_read_group_name,
		     sam_read_group_library,sam_read_group_platform);
      }
    }
  }

  return new;
}

#else

T
Outbuffer_new (unsigned int output_buffer_size, unsigned int nread, char *sevenway_root,
	       bool chimerap, char *user_genomicseg, Sequence_T usersegment,
	       char *dbversion, Genome_T genome, IIT_T chromosome_iit,
	       Chrsubset_T chrsubset, IIT_T contig_iit, IIT_T altstrain_iit, IIT_T map_iit,
	       int *map_divint_crosstable, Printtype_T printtype, bool checksump, int chimera_margin,
#ifndef PMAP
	       bool sam_headers_p, int quality_shift, bool sam_paired_p,
	       char *sam_read_group_id, char *sam_read_group_name,
	       char *sam_read_group_library, char *sam_read_group_platform,
#endif
	       bool nofailsp, bool failsonlyp, bool fails_as_input_p, int maxpaths, bool quiet_if_excessive_p,
	       bool map_exons_p, bool map_bothstrands_p, bool print_comment_p, int nflanking,
	       int proteinmode, int invertmode, bool nointronlenp, int wraplength, int ngap, int cds_startpos,
	       bool fulllengthp, bool truncatep, bool strictp, bool diagnosticp, bool maponlyp,
	       bool stage1debug, bool diag_debug, bool debug_graphic_p,
	       int argc, char **argv, int optind) {

  T new = (T) MALLOC(sizeof(*new));

  new->chimerap = chimerap;

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
  new->fp_nomapping = NULL;
  new->fp_uniq = NULL;
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
#endif

  new->nofailsp = nofailsp;
  new->failsonlyp = failsonlyp;
  new->fails_as_input_p = fails_as_input_p;
  new->maxpaths = maxpaths;
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
	IIT_dump_sam(stdout,chromosome_iit,sam_read_group_id,sam_read_group_name,
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
  Shortread_print_header(fp,queryseq1);
  /* fprintf(fp,"\n"); -- included in header */

  return;
}


static void
print_result_sam (T this, Result_T result, Request_T request) {
  Resulttype_T resulttype;
  Shortread_T queryseq1;
  Stage3end_T *stage3array, stage3;
  Genomicpos_T chrpos;
  int ignore = 0;
  int npaths, pathnum, second_absmq;
  FILE *fp;

  resulttype = Result_resulttype(result);

  if (resulttype == SINGLEEND_NOMAPPING) {
    if (this->nofailsp == true) {
      /* Skip */
    } else {
      if (this->fails_as_input_p == true) {
	print_query_singleend(this,this->fp_nomapping_1,request);
      } else {
	queryseq1 = Request_queryseq1(request);
	SAM_print_nomapping(this->fp_nomapping_1,queryseq1,/*mate*/NULL,/*acc*/Shortread_accession(queryseq1),
			    this->chromosome_iit,resulttype,
			    /*first_read_p*/true,/*nhits_mate*/0,/*mate_chrpos*/0U,
			    this->quality_shift,this->sam_read_group_id,this->invert_first_p,this->invert_second_p);
      }
    }

  } else if (resulttype == SINGLEEND_UNIQ || resulttype == SINGLEEND_TRANSLOC) {
    if (resulttype == SINGLEEND_UNIQ) {
      fp = this->fp_unpaired_uniq;
    } else {
      fp = this->fp_unpaired_transloc;
    }
    stage3array = (Stage3end_T *) Result_array(&npaths,&second_absmq,result);

    if (this->failsonlyp == true) {
      /* Skip */
    } else {
      queryseq1 = Request_queryseq1(request);
      /* Stage3end_eval_and_sort(stage3array,npaths,this->maxpaths,queryseq1); */

      stage3 = stage3array[0];
      chrpos = SAM_compute_chrpos(/*hardclip_low*/&ignore,/*hardclip_high*/&ignore,stage3,
				  Stage3end_substring_low(stage3),Shortread_fulllength(queryseq1));
      SAM_print(fp,stage3,/*mate*/NULL,/*acc*/Shortread_accession(queryseq1),/*pathnum*/1,npaths,
		Stage3end_absmq_score(stage3array[0]),second_absmq,
		Stage3end_mapq_score(stage3array[0]),
		this->chromosome_iit,queryseq1,/*queryseq2*/NULL,
		/*pairedlength*/0,chrpos,/*mate_chrpos*/0U,
		/*hardclip_low*/0,/*hardclip_high*/0,resulttype,
		/*first_read_p*/true,/*npaths_mate*/0,this->quality_shift,
		this->sam_read_group_id,this->invert_first_p,this->invert_second_p,
		this->merge_samechr_p);
    }

  } else if (resulttype == SINGLEEND_MULT) {
    stage3array = (Stage3end_T *) Result_array(&npaths,&second_absmq,result);

    if (this->failsonlyp == true) {
      /* Skip */

    } else if (this->quiet_if_excessive_p && npaths > this->maxpaths) {
      queryseq1 = Request_queryseq1(request);
      /* Stage3end_eval_and_sort(stage3array,npaths,this->maxpaths,queryseq1); */
      SAM_print_nomapping(this->fp_unpaired_mult,queryseq1,/*mate*/NULL,/*acc*/Shortread_accession(queryseq1),
			  this->chromosome_iit,resulttype,
			  /*first_read_p*/true,/*nhits_mate*/0,/*mate_chrpos*/0U,
			  this->quality_shift,this->sam_read_group_id,this->invert_first_p,this->invert_second_p);

    } else {
      queryseq1 = Request_queryseq1(request);
      /* Stage3end_eval_and_sort(stage3array,npaths,this->maxpaths,queryseq1); */
      for (pathnum = 1; pathnum <= npaths && pathnum <= this->maxpaths; pathnum++) {

	stage3 = stage3array[pathnum-1];
	chrpos = SAM_compute_chrpos(/*hardclip_low*/&ignore,/*hardclip_high*/&ignore,stage3,
				    Stage3end_substring_low(stage3),Shortread_fulllength(queryseq1));
	SAM_print(this->fp_unpaired_mult,stage3,/*mate*/NULL,/*acc*/Shortread_accession(queryseq1),
		  pathnum,npaths,
		  Stage3end_absmq_score(stage3array[pathnum-1]),second_absmq,
		  Stage3end_mapq_score(stage3array[pathnum-1]),
		  this->chromosome_iit,queryseq1,/*queryseq2*/NULL,
		  /*pairedlength*/0,chrpos,/*mate_chrpos*/0U,
		  /*hardclip_low*/0,/*hardclip_high*/0,resulttype,
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
		     this->fp_unpaired_uniq,this->fp_unpaired_transloc,this->fp_unpaired_mult,
		     this->fp_halfmapping_uniq,this->fp_halfmapping_transloc,this->fp_halfmapping_mult,
		     this->fp_paired_uniq_inv,this->fp_paired_uniq_scr,
		     this->fp_paired_uniq_long,this->fp_paired_mult,
		     this->fp_concordant_uniq,this->fp_concordant_transloc,this->fp_concordant_mult);
  }

  return;
}


static void
print_result_goby (T this, Result_T result, Request_T request) {

  Resulttype_T resulttype;
  Shortread_T queryseq1;
  Stage3end_T *stage3array, stage3;
  int npaths, pathnum, second_absmq;

  resulttype = Result_resulttype(result);

  if (resulttype == SINGLEEND_NOMAPPING || resulttype == PAIREDEND_NOMAPPING) {
    /* No output in Goby */

  } else if (resulttype == SINGLEEND_UNIQ || resulttype == SINGLEEND_TRANSLOC) {
    if (this->failsonlyp == true) {
      /* Skip */

    } else {
      stage3array = (Stage3end_T *) Result_array(&npaths,&second_absmq,result);
      stage3 = stage3array[0];
      queryseq1 = Request_queryseq1(request);
#if 0
      Stage3end_eval_and_sort(stage3array,/*npaths*/1,this->maxpaths,queryseq1);
#endif

      Goby_observe_aligned(this->gobywriter);
      Goby_print_single(this->gobywriter,stage3,Stage3end_score(stage3),
			this->chromosome_iit,queryseq1,this->invert_first_p,
			/*hit5*/(Stage3end_T) NULL,/*hit3*/(Stage3end_T) NULL,
			/*pairlength*/0,/*pairscore*/0,/*pairtype*/UNPAIRED,
			Stage3end_mapq_score(stage3));
    }

  } else if (resulttype == SINGLEEND_MULT) {
    stage3array = (Stage3end_T *) Result_array(&npaths,&second_absmq,result);
    queryseq1 = Request_queryseq1(request);

    if (this->failsonlyp == true) {
      /* Skip */
      
    } else if (npaths > this->maxpaths) {
      Goby_observe_aligned(this->gobywriter);
      Goby_print_tmh(this->gobywriter,stage3array[0],queryseq1,npaths);

    } else {
#if 0
      Stage3end_eval_and_sort(stage3array,npaths,this->maxpaths,queryseq1);
#endif
      Goby_observe_aligned(this->gobywriter);
      for (pathnum = 1; pathnum <= npaths; pathnum++) {
	stage3 = stage3array[pathnum-1];
	Goby_print_single(this->gobywriter,stage3,Stage3end_score(stage3),
			  this->chromosome_iit,queryseq1,this->invert_first_p,
			  /*hit5*/(Stage3end_T) NULL,/*hit3*/(Stage3end_T) NULL,
			  /*pairlength*/0,/*pairscore*/0,/*pairtype*/UNPAIRED,
			  Stage3end_mapq_score(stage3));
      }
    }

  } else {
    Goby_observe_aligned(this->gobywriter);
    Goby_print_paired(this->gobywriter,result,resulttype,this->chromosome_iit,
		      Request_queryseq1(request),Request_queryseq2(request),
                      this->maxpaths,this->quiet_if_excessive_p,
                      this->invert_first_p,this->invert_second_p,
                      this->nofailsp,this->failsonlyp,this->fails_as_input_p,
                      this->fastq_format_p,this->quality_shift,
		      this->sam_read_group_id);
  }

  return;
}


static void
print_result_gsnap (T this, Result_T result, Request_T request) {
  Resulttype_T resulttype;
  bool translocationp;
  Shortread_T queryseq1;
  Stage3end_T *stage3array, stage3;
  int npaths, pathnum, second_absmq;
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

  } else if (resulttype == SINGLEEND_UNIQ || resulttype == SINGLEEND_TRANSLOC) {
    if (resulttype == SINGLEEND_UNIQ) {
      translocationp = false;
      fp = this->fp_unpaired_uniq;
    } else {
      translocationp = true;
      fp = this->fp_unpaired_transloc;
    }

    stage3array = (Stage3end_T *) Result_array(&npaths,&second_absmq,result);

    if (this->failsonlyp == true) {
      /* Skip */
    } else {
      print_header_singleend(this,this->fp_unpaired_uniq,request,translocationp,/*npaths*/1);

      queryseq1 = Request_queryseq1(request);
#if 0
      Stage3end_eval_and_sort(stage3array,/*npaths*/1,this->maxpaths,queryseq1);
#endif
      stage3 = stage3array[0];
      Stage3end_print(fp,stage3,Stage3end_score(stage3),
		      this->chromosome_iit,queryseq1,this->invert_first_p,
		      /*hit5*/(Stage3end_T) NULL,/*hit3*/(Stage3end_T) NULL,
		      /*pairlength*/0,/*pairscore*/0,/*pairtype*/UNPAIRED,
		      Stage3end_mapq_score(stage3));
      fprintf(this->fp_unpaired_uniq,"\n");
    }

  } else if (resulttype == SINGLEEND_MULT) {
    stage3array = (Stage3end_T *) Result_array(&npaths,&second_absmq,result);

    if (this->failsonlyp == true) {
      /* Skip */

    } else if (this->quiet_if_excessive_p && npaths > this->maxpaths) {
      print_header_singleend(this,this->fp_unpaired_mult,request,/*translocationp*/false,npaths);
      fprintf(this->fp_unpaired_mult,"\n");

    } else {
      print_header_singleend(this,this->fp_unpaired_mult,request,/*translocationp*/false,npaths);

      queryseq1 = Request_queryseq1(request);
#if 0
      Stage3end_eval_and_sort(stage3array,npaths,this->maxpaths,queryseq1);
#endif
      for (pathnum = 1; pathnum <= npaths && pathnum <= this->maxpaths; pathnum++) {
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
		     this->maxpaths,this->quiet_if_excessive_p,
#if 0
		     this->invert_first_p,this->invert_second_p,
#endif
		     this->nofailsp,this->failsonlyp,this->fails_as_input_p,
		     this->fastq_format_p,this->quality_shift,
		     this->fp_nomapping_1,this->fp_nomapping_2,
		     this->fp_unpaired_uniq,this->fp_unpaired_transloc,this->fp_unpaired_mult,
		     this->fp_halfmapping_uniq,this->fp_halfmapping_transloc,this->fp_halfmapping_mult,
		     this->fp_paired_uniq_inv,this->fp_paired_uniq_scr,
		     this->fp_paired_uniq_long,this->fp_paired_mult,
		     this->fp_concordant_uniq,this->fp_concordant_transloc,this->fp_concordant_mult);
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
	      Chrsubset_T chrsubset, Chimera_T chimera, Failure_T failuretype) {

  if (this->diagnosticp == true) {
    Diagnostic_print(diagnostic);
  }

  if (npaths == 0) {
    fprintf(fp,"Paths (0):");
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
Outbuffer_print_result (T this, Result_T result, Request_T request
#ifdef MEMUSAGE
			, unsigned int noutput
#endif
			) {
  FILE *fp;
  Sequence_T queryseq;
  Diagnostic_T diagnostic;
  Stage3_T *stage3array;
  int npaths, pathnum, effective_maxpaths, second_absmq;
  Chimera_T chimera = NULL;
  int chimerapos, chimeraequivpos, chimera_cdna_direction;
  int querylength;
  double donor_prob, acceptor_prob;
  List_T p;
  Gregion_T gregion;
  bool printp;

  queryseq = Request_queryseq(request);

  if (this->stage1debug == true) {
    putc('>',stdout);
    Sequence_print_header(stdout,queryseq,this->checksump);

    for (p = Result_gregionlist(result); p != NULL; p = List_next(p)) {
      gregion = (Gregion_T) List_head(p);
      Gregion_print(gregion);
    }
    return;

  } else if (this->diag_debug == true) {
    putc('>',stdout);
    Sequence_print_header(stdout,queryseq,this->checksump);

    Diag_print_segments(Result_diagonals(result),/*queryseq_ptr*/NULL,/*genomicseg_ptr*/NULL);
    return;
  }

  stage3array = Result_array(&npaths,&second_absmq,result);
  querylength = Sequence_fulllength_given(queryseq);

  chimerapos = chimeraequivpos = -1;
  chimera_cdna_direction = 0;
  donor_prob = acceptor_prob = 0.0;

  /* Translation */
  if (npaths == 0) {
    effective_maxpaths = 0;
    fp = this->fp_nomapping;

    if (this->nofailsp == true) {
      printp = false;
      if (this->fails_as_input_p == true) {
	putc('>',fp);
	Sequence_print_header(fp,queryseq,this->checksump);
	Sequence_print(fp,queryseq,/*uppercasep*/false,this->wraplength,/*trimmedp*/false);
      }
    } else {
      printp = true;
    }

    if (Result_failuretype(result) == POOR_SEQUENCE) {
      fprintf(stderr,"Accession %s skipped (poor sequence).  Use -p flag to change pruning behavior\n",Sequence_accession(queryseq));
    } else if (Result_failuretype(result) == REPETITIVE) {
      fprintf(stderr,"Accession %s skipped (repetitive sequence).  Use -p flag to change pruning behavior\n",Sequence_accession(queryseq));
    } else {
      fprintf(stderr,"No paths found for %s\n",Sequence_accession(queryseq));
    }

  } else if ((chimera = Result_chimera(result)) != NULL) {
    effective_maxpaths = 2;
    fp = this->fp_transloc;

    if (this->failsonlyp == true) {
      printp = false;

#if 0
    } else if (this->quiet_if_excessive_p && npaths > this->maxpaths) {
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
			       queryseq,querylength,this->fulllengthp,
			       this->cds_startpos,this->truncatep,this->strictp,
			       this->diagnosticp,this->maponlyp);
    }

  } else if (this->maxpaths == 0) {
    effective_maxpaths = 1;
    if (npaths == 1) {
      fp = this->fp_uniq;
    } else {
      fp = this->fp_mult;
    }

    if (this->failsonlyp == true) {
      printp = false;
    } else if (this->quiet_if_excessive_p && npaths > this->maxpaths) {
      printp = false;
    } else {
      printp = true;

      Stage3_translate(stage3array[0],queryseq,querylength,this->fulllengthp,
		       this->cds_startpos,this->truncatep,this->strictp,
		       this->diagnosticp,this->maponlyp);
    }

  } else {
    if (npaths == 1) {
      fp = this->fp_uniq;
    } else {
      fp = this->fp_mult;
    }

    if (npaths < this->maxpaths) {
      effective_maxpaths = npaths;
    } else {
      effective_maxpaths = this->maxpaths;
    }

    if (this->failsonlyp == true) {
      printp = false;
    } else if (this->quiet_if_excessive_p && npaths > this->maxpaths) {
      printp = false;
    } else {
      printp = true;

      for (pathnum = 1; pathnum <= effective_maxpaths; pathnum++) {
	Stage3_translate(stage3array[pathnum-1],queryseq,querylength,this->fulllengthp,
			 this->cds_startpos,this->truncatep,this->strictp,
			 this->diagnosticp,this->maponlyp);
      }
    }
  }

  /* Printing */
  if (this->debug_graphic_p == true) {
    /* No output */

  } else if (printp == false) {
    /* No output */

  } else if (this->printtype == SIMPLE || this->printtype == SUMMARY || this->printtype == ALIGNMENT) {
    /* Print header, even if no alignment is found */
    putc('>',fp);
    Sequence_print_header(fp,queryseq,this->checksump);

    diagnostic = Result_diagnostic(result);
    if (npaths == 0) {
      print_npaths(this,fp,0,diagnostic,this->chrsubset,/*chimera*/NULL,Result_failuretype(result));
    } else {
      print_npaths(this,fp,npaths,diagnostic,this->chrsubset,chimera,NO_FAILURE);
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
	Stage3_print_alignment(fp,stage3array[pathnum-1],queryseq,
			       this->genome,this->chromosome_iit,this->printtype,
			       /*continuousp*/false,/*continuous_by_exon_p*/false,
			       this->diagnosticp,this->strictp,/*flipgenomep*/true,
			       this->invertmode,this->nointronlenp,
			       this->wraplength,this->maponlyp);
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
			      donor_prob,acceptor_prob,chimera_cdna_direction,
			      this->truncatep,this->strictp);
    }

  } else if (this->printtype == CONTINUOUS) {
    putc('>',fp);
    Sequence_print_header(fp,queryseq,this->checksump);
    if (npaths == 0) {
      fprintf(fp,"\n\n\n");
    } else {
      Stage3_print_alignment(fp,stage3array[0],queryseq,this->genome,this->chromosome_iit,this->printtype,
			     /*continuousp*/true,/*continuous_by_exon_p*/false,
			     this->diagnosticp,this->strictp,/*flipgenomep*/true,
			     this->invertmode,this->nointronlenp,
			     this->wraplength,this->maponlyp);
    }

  } else if (this->printtype == CONTINUOUS_BY_EXON) {
    putc('>',fp);
    Sequence_print_header(fp,queryseq,this->checksump);
    print_npaths(this,fp,npaths,diagnostic,this->chrsubset,chimera,NO_FAILURE);
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
      Stage3_print_alignment(fp,stage3array[0],queryseq,this->genome,this->chromosome_iit,this->printtype,
			     /*continuousp*/false,/*continuous_by_exon_p*/true,
			     this->diagnosticp,this->strictp,/*flipgenomep*/true,
			     this->invertmode,this->nointronlenp,
			     this->wraplength,this->maponlyp);
    }

  } else if (this->printtype == EXONS_CDNA) {
    putc('>',fp);
    Sequence_print_header(fp,queryseq,this->checksump);
    for (pathnum = 1; pathnum <= effective_maxpaths; pathnum++) {
      fprintf(fp,"<path %d>\n",pathnum);
      Pair_print_exons(fp,Stage3_pairarray(stage3array[0]),Stage3_npairs(stage3array[0]),
		       this->wraplength,this->ngap,/*cdna*/true);
      fprintf(fp,"</path>\n");
    }

  } else if (this->printtype == EXONS_GENOMIC) {
    putc('>',fp);
    Sequence_print_header(fp,queryseq,this->checksump);
    for (pathnum = 1; pathnum <= effective_maxpaths; pathnum++) {
      fprintf(fp,"<path %d>\n",pathnum);
      Pair_print_exons(fp,Stage3_pairarray(stage3array[0]),Stage3_npairs(stage3array[0]),
		       this->wraplength,this->ngap,/*cdna*/false);
      fprintf(fp,"</path>\n");
    }

  } else if (this->printtype == CDNA) {
    for (pathnum = 1; pathnum <= effective_maxpaths; pathnum++) {
      putc('>',fp);
      Sequence_print_header(fp,queryseq,this->checksump);
      Stage3_print_cdna(fp,stage3array[pathnum-1],queryseq,this->wraplength);
    }

  } else if (this->printtype == PROTEIN_GENOMIC) {
    for (pathnum = 1; pathnum <= effective_maxpaths; pathnum++) {
      putc('>',fp);
      Sequence_print_header(fp,queryseq,this->checksump);
      Stage3_print_protein_genomic(fp,stage3array[pathnum-1],queryseq,this->wraplength);
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
      Pair_print_sam_nomapping(fp,Sequence_accession(queryseq),
			       Sequence_fullpointer(queryseq),Sequence_quality_string(queryseq),
			       Sequence_fulllength(queryseq),this->quality_shift,
			       Sequence_firstp(queryseq),this->sam_paired_p,this->sam_read_group_id);

    } else if (this->quiet_if_excessive_p && npaths > this->maxpaths) {
      Pair_print_sam_nomapping(fp,Sequence_accession(queryseq),
			       Sequence_fullpointer(queryseq),Sequence_quality_string(queryseq),
			       Sequence_fulllength(queryseq),this->quality_shift,
			       Sequence_firstp(queryseq),this->sam_paired_p,this->sam_read_group_id);

    } else if (chimera != NULL) {
      Stage3_print_sam(fp,stage3array[0],/*pathnum*/1,npaths,
		       Stage3_absmq_score(stage3array[0]),second_absmq,
		       Stage3_mapq_score(stage3array[0]),
		       this->chromosome_iit,this->usersegment,queryseq,
		       /*chimera_part*/-1,chimera,this->quality_shift,this->sam_paired_p,
		       this->sam_read_group_id);
      Stage3_print_sam(fp,stage3array[1],/*pathnum*/1,npaths,
		       Stage3_absmq_score(stage3array[0]),second_absmq,
		       Stage3_mapq_score(stage3array[0]),
		       this->chromosome_iit,this->usersegment,queryseq,
		       /*chimera_part*/+1,chimera,this->quality_shift,this->sam_paired_p,
		       this->sam_read_group_id);

    } else {
      for (pathnum = 1; pathnum <= effective_maxpaths; pathnum++) {
	Stage3_print_sam(fp,stage3array[pathnum-1],pathnum,npaths,
			 Stage3_absmq_score(stage3array[pathnum-1]),second_absmq,
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
      Sequence_print_header(fp,queryseq,this->checksump);
      Stage3_print_coordinates(fp,stage3array[pathnum-1],queryseq,
			       this->chromosome_iit,this->invertmode);
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

  } else if (this->printtype == MAP_GENES) {
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
  printf("Memusage: %ld\n",Mem_usage_report());
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
      Outbuffer_print_result(this,result,request,noutput+1);
#else
      Outbuffer_print_result(this,result,request);
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
	  Outbuffer_print_result(this,result,request,noutput+1);
#else
	  Outbuffer_print_result(this,result,request);
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
	Outbuffer_print_result(this,result,request,noutput+1);
#else
	Outbuffer_print_result(this,result,request);
#endif
	Result_free(&result);
	Request_free(&request);
	noutput++;
	
	/* Print out rest of stored queue */
	while (queue != NULL && queue->id == (int) noutput) {
	  queue = RRlist_pop_id(queue,&id,&request,&result);
	  nqueued--;
#ifdef MEMUSAGE
	  Outbuffer_print_result(this,result,request,noutput+1);
#else
	  Outbuffer_print_result(this,result,request);
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
	    Outbuffer_print_result(this,result,request,noutput+1);
#else
	    Outbuffer_print_result(this,result,request);
#endif
	    Result_free(&result);
	    Request_free(&request);
	    noutput++;
	
	    /* Print out rest of stored queue */
	    while (queue != NULL && queue->id == (int) noutput) {
	      queue = RRlist_pop_id(queue,&id,&request,&result);
	      nqueued--;
#ifdef MEMUSAGE
	      Outbuffer_print_result(this,result,request,noutput+1);
#else
	      Outbuffer_print_result(this,result,request);
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


