static char rcsid[] = "$Id: outbuffer.c 37254 2011-03-28 16:34:08Z twu $";
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
#else
#include "indexdb.h"		/* For INDEX1PART */
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
    printf("%p: prev %p, next %p\n",this,this->prev,this->next);
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

  new = (RRlist_T) MALLOC(sizeof(*new));
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

  FREE(head);
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

  new = (RRlist_T) MALLOC(sizeof(*new));
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

  FREE(head);
  return newhead;
}




#define T Outbuffer_T
struct T {

#ifdef GSNAP
  UINT4 *genome_blocks;
#else
  Genome_T genome;
#endif

  IIT_T chromosome_iit;

  char *sevenway_root;

#ifdef GSNAP
  bool sam_headers_p;
  char *sam_read_group_id;
  char *sam_read_group_name;
  int quality_shift;
#elif defined PMAP

#else
  bool sam_headers_p;
  bool sam_paired_p;
  bool cigar_noncanonical_splices_p;
  char *sam_read_group_id;
  char *sam_read_group_name;
  int quality_shift;
#endif

#ifdef GSNAP

  FILE *fp_nomapping_1;
  FILE *fp_nomapping_2;
  FILE *fp_halfmapping_uniq;
  FILE *fp_halfmapping_mult;
  FILE *fp_unpaired_uniq;
  FILE *fp_unpaired_mult;
  FILE *fp_paired_uniq_inv;
  FILE *fp_paired_uniq_scr;
  FILE *fp_paired_uniq_long;
  FILE *fp_paired_mult;
  FILE *fp_concordant_uniq;
  FILE *fp_concordant_mult;

  bool output_sam_p;
  Gobywriter_T gobywriter;

  bool fastq_format_p;


  bool invert_first_p;
  bool invert_second_p;
  Genomicpos_T pairmax;

#else

  FILE *fp_nomapping;
  FILE *fp_uniq;
  FILE *fp_mult;

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

  unsigned int nread;
  unsigned int ndone;
  unsigned int noutput;

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
      IIT_dump_sam(this->fp_nomapping_1,this->chromosome_iit,this->sam_read_group_id,this->sam_read_group_name);
    }
  }

  filename = (char *) CALLOC(strlen(this->sevenway_root)+strlen(".unpaired_uniq")+1,sizeof(char));
  sprintf(filename,"%s.unpaired_uniq",this->sevenway_root);
  if ((this->fp_unpaired_uniq = fopen(filename,"w")) == NULL) {
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
    IIT_dump_sam(this->fp_unpaired_uniq,this->chromosome_iit,this->sam_read_group_id,this->sam_read_group_name);
    IIT_dump_sam(this->fp_unpaired_mult,this->chromosome_iit,this->sam_read_group_id,this->sam_read_group_name);
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
      IIT_dump_sam(this->fp_nomapping_1,this->chromosome_iit,this->sam_read_group_id,this->sam_read_group_name);
    }
  }

  filename = (char *) CALLOC(strlen(this->sevenway_root)+strlen(".halfmapping_uniq")+1,sizeof(char));
  sprintf(filename,"%s.halfmapping_uniq",this->sevenway_root);
  if ((this->fp_halfmapping_uniq = fopen(filename,"w")) == NULL) {
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

  filename = (char *) CALLOC(strlen(this->sevenway_root)+strlen(".concordant_mult")+1,sizeof(char));
  sprintf(filename,"%s.concordant_mult",this->sevenway_root);
  if ((this->fp_concordant_mult = fopen(filename,"w")) == NULL) {
    fprintf(stderr,"Cannot open file %s for writing\n",filename);
    exit(9);
  }
  FREE(filename);

  if (this->output_sam_p == true && this->sam_headers_p == true) {
    IIT_dump_sam(this->fp_halfmapping_uniq,this->chromosome_iit,this->sam_read_group_id,this->sam_read_group_name);
    IIT_dump_sam(this->fp_halfmapping_mult,this->chromosome_iit,this->sam_read_group_id,this->sam_read_group_name);
    IIT_dump_sam(this->fp_paired_uniq_inv,this->chromosome_iit,this->sam_read_group_id,this->sam_read_group_name);
    IIT_dump_sam(this->fp_paired_uniq_scr,this->chromosome_iit,this->sam_read_group_id,this->sam_read_group_name);
    IIT_dump_sam(this->fp_paired_uniq_long,this->chromosome_iit,this->sam_read_group_id,this->sam_read_group_name);
    IIT_dump_sam(this->fp_paired_mult,this->chromosome_iit,this->sam_read_group_id,this->sam_read_group_name);
    IIT_dump_sam(this->fp_concordant_uniq,this->chromosome_iit,this->sam_read_group_id,this->sam_read_group_name);
    IIT_dump_sam(this->fp_concordant_mult,this->chromosome_iit,this->sam_read_group_id,this->sam_read_group_name);
  }

  return;
}

static void
sevenway_close (T this) {
  fclose(this->fp_unpaired_uniq);
  fclose(this->fp_unpaired_mult);
  if (this->fp_nomapping_2 != NULL) {
    fclose(this->fp_nomapping_2);
  }
  if (this->fp_nomapping_1 != NULL) {
    fclose(this->fp_nomapping_1);
  }
  if (this->fp_halfmapping_uniq != NULL) {
    fclose(this->fp_halfmapping_uniq);
    fclose(this->fp_halfmapping_mult);
    fclose(this->fp_paired_uniq_long);
    fclose(this->fp_paired_uniq_scr);
    fclose(this->fp_paired_uniq_inv);
    fclose(this->fp_paired_mult);
    fclose(this->fp_concordant_uniq);
    fclose(this->fp_concordant_mult);
  }

  return;
}

#else

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
      IIT_dump_sam(this->fp_nomapping,this->chromosome_iit,this->sam_read_group_id,this->sam_read_group_name);
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
    IIT_dump_sam(this->fp_uniq,this->chromosome_iit,this->sam_read_group_id,this->sam_read_group_name);
    IIT_dump_sam(this->fp_mult,this->chromosome_iit,this->sam_read_group_id,this->sam_read_group_name);
#endif
  }

  return;
}

static void
sevenway_close (T this) {
  fclose(this->fp_mult);
  fclose(this->fp_uniq);
  fclose(this->fp_nomapping);
  return;
}

#endif



#ifdef GSNAP

T
Outbuffer_new (int nread, char *sevenway_root, UINT4 *genome_blocks, IIT_T chromosome_iit,
	       bool output_sam_p, bool sam_headers_p, char *sam_read_group_id, char *sam_read_group_name,
	       Gobywriter_T gobywriter, bool nofailsp, bool failsonlyp, bool fails_as_input_p,
	       bool fastq_format_p, int maxpaths, bool quiet_if_excessive_p, int quality_shift,
	       bool invert_first_p, bool invert_second_p, Genomicpos_T pairmax) {
  T new = (T) MALLOC(sizeof(*new));

  new->genome_blocks = genome_blocks;
  new->chromosome_iit = chromosome_iit;

  new->fp_nomapping_1 = NULL;
  new->fp_nomapping_2 = NULL;
  new->fp_halfmapping_uniq = NULL;
  new->fp_halfmapping_mult = NULL;
  new->fp_unpaired_uniq = NULL;
  new->fp_unpaired_mult = NULL;
  new->fp_paired_uniq_inv = NULL;
  new->fp_paired_uniq_scr = NULL;
  new->fp_paired_uniq_long = NULL;
  new->fp_paired_mult = NULL;
  new->fp_concordant_uniq = NULL;
  new->fp_concordant_mult = NULL;
  
  new->sevenway_root = sevenway_root;

  new->output_sam_p = output_sam_p;
  new->sam_headers_p = sam_headers_p;
  new->sam_read_group_id = sam_read_group_id;
  new->sam_read_group_name = sam_read_group_name;

  new->gobywriter = gobywriter;

  new->nofailsp = nofailsp;
  new->failsonlyp = failsonlyp;
  new->fails_as_input_p = fails_as_input_p;
  new->fastq_format_p = fastq_format_p;
  new->maxpaths = maxpaths;
  new->quiet_if_excessive_p = quiet_if_excessive_p;

  new->quality_shift = quality_shift;
  new->invert_first_p = invert_first_p;
  new->invert_second_p = invert_second_p;
  new->pairmax = pairmax;

#ifdef HAVE_PTHREAD
  pthread_mutex_init(&new->lock,NULL);
#endif

  new->nread = nread;
  new->ndone = (unsigned int) -1U;
  new->noutput = 0;

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
    new->fp_halfmapping_mult = stdout;
    new->fp_unpaired_uniq = stdout;
    new->fp_unpaired_mult = stdout;
    new->fp_paired_uniq_inv = stdout;
    new->fp_paired_uniq_scr = stdout;
    new->fp_paired_uniq_long = stdout;
    new->fp_paired_mult = stdout;
    new->fp_concordant_uniq = stdout;
    new->fp_concordant_mult = stdout;

    if (output_sam_p == true && sam_headers_p == true) {
      if (fails_as_input_p == true) {
	/* Don't print chromosomes */
      } else {
	IIT_dump_sam(stdout,chromosome_iit,sam_read_group_id,sam_read_group_name);
      }
    }
  }

  return new;
}

#else

T
Outbuffer_new (int nread, char *sevenway_root, char *user_genomicseg, Sequence_T usersegment,
	       char *dbversion, Genome_T genome, IIT_T chromosome_iit,
	       Chrsubset_T chrsubset, IIT_T contig_iit, IIT_T altstrain_iit, IIT_T map_iit,
	       int *map_divint_crosstable, Printtype_T printtype, bool checksump, int chimera_margin,
#ifndef PMAP
	       bool sam_headers_p, int quality_shift, bool sam_paired_p, bool cigar_noncanonical_splices_p,
	       char *sam_read_group_id, char *sam_read_group_name,
#endif
	       bool nofailsp, bool failsonlyp, bool fails_as_input_p, int maxpaths, bool quiet_if_excessive_p,
	       bool map_exons_p, bool map_bothstrands_p, bool print_comment_p, int nflanking,
	       int proteinmode, int invertmode, bool nointronlenp, int wraplength, int ngap, int cds_startpos,
	       bool fulllengthp, bool truncatep, bool strictp, bool diagnosticp, bool maponlyp,
	       bool stage1debug, bool diag_debug, bool debug_graphic_p,
	       int argc, char **argv, int optind) {

  T new = (T) MALLOC(sizeof(*new));

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
  new->fp_mult = NULL;
  
#ifndef PMAP
  new->sam_headers_p = sam_headers_p;
  new->quality_shift = quality_shift;
  new->sam_paired_p = sam_paired_p;
  new->cigar_noncanonical_splices_p = cigar_noncanonical_splices_p;
  new->sam_read_group_id = sam_read_group_id;
  new->sam_read_group_name = sam_read_group_name;
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

  new->nread = nread;
  new->ndone = (unsigned int) -1U;
  new->noutput = 0;

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
    new->fp_mult = stdout;

#ifndef PMAP
    if (printtype == SAM && sam_headers_p == true) {
      if (fails_as_input_p == true) {
	/* Don't print chromosomes */
      } else {
	IIT_dump_sam(stdout,chromosome_iit,sam_read_group_id,sam_read_group_name);
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

#ifdef LEAKCHECK
    Mem_leak_check_deactivate();
#endif
    FREE(*old);
#ifdef LEAKCHECK
    Mem_leak_check_activate();
#endif
  }
  return;
}



int
Outbuffer_nread (T this) {
  int nread;

#ifdef HAVE_PTHREAD
  pthread_mutex_lock(&this->lock);
#endif

  nread = this->nread;

#ifdef HAVE_PTHREAD
  pthread_mutex_unlock(&this->lock);
#endif
  
  return nread;
}



void
Outbuffer_add_nread (T this, int nread) {

#ifdef HAVE_PTHREAD
  pthread_mutex_lock(&this->lock);
#endif

  if (nread == 0) {
    this->ndone = this->nread;
    debug(fprintf(stderr,"__Outbuffer_add_nread added 0 reads, so setting ndone to be %u\n",this->ndone));

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
    Shortread_print_quality(fp,queryseq1,/*hardclip_low*/0,/*hardclip_high*/0,this->quality_shift,
			   /*show_chopped_p*/true);
  }

  fprintf(fp,"\t");
  Shortread_print_header(fp,queryseq1);
  /* fprintf(fp,"\n"); -- included in header */

  return;
}


static void
print_result_sam (T this, Result_T result, Request_T request) {
  Resulttype_T resulttype;
  bool translocationp;
  Shortread_T queryseq1;
  Stage3_T *stage3array;
  int npaths, pathnum;

  resulttype = Result_resulttype(result);
  translocationp = Result_translocationp(result);

  if (resulttype == SINGLEEND_NOMAPPING) {
    if (this->nofailsp == true) {
      /* Skip */
    } else {
      if (this->fails_as_input_p == true) {
	print_query_singleend(this,this->fp_nomapping_1,request);
      } else {
	queryseq1 = Request_queryseq1(request);
	SAM_print_nomapping(this->fp_nomapping_1,queryseq1,/*mate*/NULL,/*acc*/Shortread_accession(queryseq1),
			    this->chromosome_iit,resulttype,translocationp,
			    /*first_read_p*/true,/*nhits_mate*/0,/*queryseq_mate*/NULL,
			    this->quality_shift,this->sam_read_group_id,this->invert_first_p,this->invert_second_p);
      }
    }

  } else if (resulttype == SINGLEEND_UNIQ) {
    stage3array = (Stage3_T *) Result_array(&npaths,result);

    if (this->failsonlyp == true) {
      /* Skip */
    } else {
      queryseq1 = Request_queryseq1(request);
      /* Stage3_eval_and_sort(stage3array,npaths,this->maxpaths,queryseq1); */
      SAM_print(this->fp_unpaired_uniq,stage3array[0],/*mate*/NULL,/*acc*/Shortread_accession(queryseq1),
		/*pathnum*/1,npaths,Stage3_mapq_score(stage3array[0]),
		this->chromosome_iit,queryseq1,/*queryseq2*/NULL,
		/*pairedlength*/0,resulttype,translocationp,
		/*first_read_p*/true,/*npaths_mate*/0,this->quality_shift,this->sam_read_group_id,
		this->invert_first_p,this->invert_second_p);
    }

  } else if (resulttype == SINGLEEND_MULT) {
    stage3array = (Stage3_T *) Result_array(&npaths,result);

    if (this->failsonlyp == true) {
      /* Skip */

    } else if (this->quiet_if_excessive_p && npaths > this->maxpaths) {
      queryseq1 = Request_queryseq1(request);
      /* Stage3_eval_and_sort(stage3array,npaths,this->maxpaths,queryseq1); */
      SAM_print_nomapping(this->fp_unpaired_mult,queryseq1,/*mate*/NULL,/*acc*/Shortread_accession(queryseq1),
			  this->chromosome_iit,resulttype,translocationp,
			  /*first_read_p*/true,/*nhits_mate*/0,/*queryseq_mate*/NULL,
			  this->quality_shift,this->sam_read_group_id,this->invert_first_p,this->invert_second_p);

    } else {
      queryseq1 = Request_queryseq1(request);
      /* Stage3_eval_and_sort(stage3array,npaths,this->maxpaths,queryseq1); */
      for (pathnum = 1; pathnum <= npaths && pathnum <= this->maxpaths; pathnum++) {
	SAM_print(this->fp_unpaired_mult,stage3array[pathnum-1],/*mate*/NULL,/*acc*/Shortread_accession(queryseq1),
		  pathnum,npaths,Stage3_mapq_score(stage3array[pathnum-1]),
		  this->chromosome_iit,queryseq1,/*queryseq2*/NULL,
		  /*pairedlength*/0,resulttype,translocationp,
		  /*first_read_p*/true,/*npaths_mate*/0,this->quality_shift,this->sam_read_group_id,
		  this->invert_first_p,this->invert_second_p);
      }
    }

  } else {
    if (this->fp_concordant_uniq == NULL) {
      sevenway_open_paired(this);
    }
    SAM_print_paired(result,resulttype,translocationp,this->chromosome_iit,
		     Request_queryseq1(request),Request_queryseq2(request),
		     this->maxpaths,this->quiet_if_excessive_p,this->invert_first_p,this->invert_second_p,
		     this->nofailsp,this->failsonlyp,this->fails_as_input_p,this->fastq_format_p,this->quality_shift,
		     this->sam_read_group_id,
		     this->fp_nomapping_1,this->fp_nomapping_2,
		     this->fp_unpaired_uniq,this->fp_unpaired_mult,
		     this->fp_halfmapping_uniq,this->fp_halfmapping_mult,
		     this->fp_paired_uniq_inv,this->fp_paired_uniq_scr,
		     this->fp_paired_uniq_long,this->fp_paired_mult,
		     this->fp_concordant_uniq,this->fp_concordant_mult);
  }

  return;
}


static void
print_result_goby (T this, Result_T result, Request_T request) {
  Resulttype_T resulttype;
  Shortread_T queryseq1;
  Stage3_T *stage3array;
  int npaths;

  resulttype = Result_resulttype(result);

  if (resulttype == SINGLEEND_NOMAPPING || resulttype == PAIREDEND_NOMAPPING) {
    /* No output in Goby */

  } else if (resulttype == SINGLEEND_UNIQ || resulttype == SINGLEEND_MULT) {
    stage3array = (Stage3_T *) Result_array(&npaths,result);
    queryseq1 = Request_queryseq1(request);

    if (this->failsonlyp == true) {
      /* Skip */

    } else {
#if 0
      if (this->quiet_if_excessive_p && npaths > this->maxpaths) {
	/* Don't sort */
      } else {
	Stage3_eval_and_sort(stage3array,npaths,this->maxpaths,queryseq1);
      }
#endif
      Goby_print_single(this->gobywriter,this->chromosome_iit,stage3array,queryseq1,npaths,this->maxpaths,this->quiet_if_excessive_p);
    }

  } else {
    Goby_print_paired(this->gobywriter,result,resulttype,this->chromosome_iit,
		      Request_queryseq1(request),Request_queryseq2(request),
		      this->maxpaths,this->quiet_if_excessive_p,this->invert_first_p,this->invert_second_p,
		      this->nofailsp,this->failsonlyp,this->fails_as_input_p,this->fastq_format_p,this->quality_shift,
		      this->sam_read_group_id);
  }

  return;
}



static void
print_result_gsnap (T this, Result_T result, Request_T request) {
  Resulttype_T resulttype;
  bool translocationp;
  Shortread_T queryseq1;
  Stage3_T *stage3array, stage3;
  int npaths, pathnum;

  resulttype = Result_resulttype(result);
  translocationp = Result_translocationp(result);

  if (resulttype == SINGLEEND_NOMAPPING) {
    if (this->nofailsp == true) {
      /* Skip */
    } else {
      if (this->fails_as_input_p == true) {
	print_query_singleend(this,this->fp_nomapping_1,request);
      } else {
	print_header_singleend(this,this->fp_nomapping_1,request,translocationp,/*npaths*/0);
	fprintf(this->fp_nomapping_1,"\n");
      }
    }

  } else if (resulttype == SINGLEEND_UNIQ) {
    stage3array = (Stage3_T *) Result_array(&npaths,result);

    if (this->failsonlyp == true) {
      /* Skip */
    } else {
      print_header_singleend(this,this->fp_unpaired_uniq,request,translocationp,/*npaths*/1);

      queryseq1 = Request_queryseq1(request);
#if 0
      Stage3_eval_and_sort(stage3array,/*npaths*/1,this->maxpaths,queryseq1);
#endif
      stage3 = stage3array[0];
      Stage3_print(this->fp_unpaired_uniq,stage3,Stage3_score(stage3),this->genome_blocks,
		   this->chromosome_iit,queryseq1,this->invert_first_p,
		   /*hit5*/(Stage3_T) NULL,/*hit3*/(Stage3_T) NULL,
		   /*pairlength*/0,/*pairscore*/0,/*pairtype*/UNPAIRED,
		   Stage3_mapq_score(stage3));
      fprintf(this->fp_unpaired_uniq,"\n");
    }

  } else if (resulttype == SINGLEEND_MULT) {
    stage3array = (Stage3_T *) Result_array(&npaths,result);

    if (this->failsonlyp == true) {
      /* Skip */

    } else if (this->quiet_if_excessive_p && npaths > this->maxpaths) {
      print_header_singleend(this,this->fp_unpaired_mult,request,translocationp,npaths);
      fprintf(this->fp_unpaired_mult,"\n");

    } else {
      print_header_singleend(this,this->fp_unpaired_mult,request,translocationp,npaths);

      queryseq1 = Request_queryseq1(request);
#if 0
      Stage3_eval_and_sort(stage3array,npaths,this->maxpaths,queryseq1);
#endif
      for (pathnum = 1; pathnum <= npaths && pathnum <= this->maxpaths; pathnum++) {
	stage3 = stage3array[pathnum-1];
	Stage3_print(this->fp_unpaired_mult,stage3,Stage3_score(stage3),this->genome_blocks,
		     this->chromosome_iit,queryseq1,this->invert_first_p,
		     /*hit5*/(Stage3_T) NULL,/*hit3*/(Stage3_T) NULL,
		     /*pairlength*/0,/*pairscore*/0,/*pairtype*/UNPAIRED,
		     Stage3_mapq_score(stage3));
      }
      fprintf(this->fp_unpaired_mult,"\n");
    }

  } else {
    if (this->fp_concordant_uniq == NULL) {
      sevenway_open_paired(this);
    }
    Stage3_print_paired(result,resulttype,translocationp,this->genome_blocks,this->chromosome_iit,
			Request_queryseq1(request),Request_queryseq2(request),this->pairmax,
			this->maxpaths,this->quiet_if_excessive_p,
#if 0
			this->invert_first_p,this->invert_second_p,
#endif
			this->nofailsp,this->failsonlyp,this->fails_as_input_p,
			this->fastq_format_p,this->quality_shift,
			this->fp_nomapping_1,this->fp_nomapping_2,
			this->fp_unpaired_uniq,this->fp_unpaired_mult,
			this->fp_halfmapping_uniq,this->fp_halfmapping_mult,
			this->fp_paired_uniq_inv,this->fp_paired_uniq_scr,
			this->fp_paired_uniq_long,this->fp_paired_mult,
			this->fp_concordant_uniq,this->fp_concordant_mult);
  }

  return;
}



void
Outbuffer_print_result (T this, Result_T result, Request_T request) {

#ifdef MEMUSAGE
  printf("Memusage: %lu\n",Mem_usage_report());
#endif

  if (this->output_sam_p == true) {
    print_result_sam(this,result,request);
  } else if (this->gobywriter != NULL) {
    print_result_goby(this,result,request);
  } else {
    print_result_gsnap(this,result,request);
  }
  return;
}

#else

/************************************************************************
 *   Print routines and threads for GMAP
 ************************************************************************/

static void
print_npaths (T this, FILE *fp, int npaths, Diagnostic_T diagnostic, Sequence_T usersegment,
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
#ifdef PMAP
    fprintf(fp," *** Short sequence < %d aa ***",INDEX1PART_AA);
#else
    fprintf(fp," *** Short sequence < %d bp ***",INDEX1PART);
#endif
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
Outbuffer_print_result (T this, Result_T result, Request_T request) {
  FILE *fp;
  Sequence_T queryseq;
  Diagnostic_T diagnostic;
  Stage3_T *stage3array;
  int npaths, pathnum, effective_maxpaths;
  Chimera_T chimera = NULL;
  int chimerapos, chimeraequivpos, chimera_cdna_direction;
  int querylength;
  double donor_prob, acceptor_prob;
  List_T p;
  Gregion_T gregion;
  bool printp;

  queryseq = Request_queryseq(request);

#ifdef MEMUSAGE
  printf("Memusage: %lu\n",Mem_usage_report());
#endif

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

  stage3array = Result_array(&npaths,result);
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
    fp = this->fp_uniq;

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
      print_npaths(this,fp,0,diagnostic,this->usersegment,this->chrsubset,
		   /*chimera*/NULL,Result_failuretype(result));
    } else {
      print_npaths(this,fp,npaths,diagnostic,this->usersegment,this->chrsubset,
		   chimera,NO_FAILURE);
      for (pathnum = 1; pathnum <= effective_maxpaths; pathnum++) {
	Stage3_print_pathsummary(fp,stage3array[pathnum-1],pathnum,
				 this->chromosome_iit,this->contig_iit,
				 this->altstrain_iit,queryseq,querylength,
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
			       this->proteinmode,this->invertmode,this->nointronlenp,
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
			      this->dbversion,pathnum,npaths,this->checksump,chimerapos,chimeraequivpos,
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
			     this->proteinmode,this->invertmode,this->nointronlenp,
			     this->wraplength,this->maponlyp);
    }

  } else if (this->printtype == CONTINUOUS_BY_EXON) {
    putc('>',fp);
    Sequence_print_header(fp,queryseq,this->checksump);
    if (npaths == 0) {
      fprintf(fp,"\n\n\n");
    } else {
      Stage3_print_alignment(fp,stage3array[0],queryseq,this->genome,this->chromosome_iit,this->printtype,
			     /*continuousp*/false,/*continuous_by_exon_p*/true,
			     this->diagnosticp,this->strictp,/*flipgenomep*/true,
			     this->proteinmode,this->invertmode,this->nointronlenp,
			     this->wraplength,this->maponlyp);
    }

  } else if (this->printtype == EXONS_CDNA) {
    if (npaths > 0) {
      putc('>',fp);
      Sequence_print_header(fp,queryseq,this->checksump);
      Pair_print_exons(fp,Stage3_pairarray(stage3array[0]),Stage3_npairs(stage3array[0]),
		       this->wraplength,this->ngap,/*cdna*/true);
    }

  } else if (this->printtype == EXONS_GENOMIC) {
    if (npaths > 0) {
      putc('>',fp);
      Sequence_print_header(fp,queryseq,this->checksump);
      Pair_print_exons(fp,Stage3_pairarray(stage3array[0]),Stage3_npairs(stage3array[0]),
		       this->wraplength,this->ngap,/*cdna*/false);
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
      Stage3_print_pslformat_nt(fp,stage3array[pathnum-1],pathnum,
				this->chromosome_iit,this->usersegment,queryseq);
    }

#ifdef PMAP
  } else if (this->printtype == PSL_PRO) {
    for (pathnum = 1; pathnum <= effective_maxpaths; pathnum++) {
      Stage3_print_pslformat_pro(fp,stage3array[pathnum-1],pathnum,
				 this->chromosome_iit,this->usersegment,queryseq,this->strictp);
    }
#endif

  } else if (this->printtype == GFF3_GENE || this->printtype == GFF3_MATCH_CDNA ||
	     this->printtype == GFF3_MATCH_EST) {
    for (pathnum = 1; pathnum <= effective_maxpaths; pathnum++) {
      Stage3_print_gff3(fp,stage3array[pathnum-1],pathnum,
			this->chromosome_iit,queryseq,querylength,this->printtype,
			/*sourcename*/this->usersegment ? this->user_genomicseg : this->dbversion);
    }

#ifndef PMAP
  } else if (this->printtype == SAM) {
    if (npaths == 0) {
      Pair_print_sam_nomapping(fp,/*translocationp*/false,Sequence_accession(queryseq),queryseq,
			       this->quality_shift,this->sam_paired_p,this->sam_read_group_id);

    } else if (chimera != NULL) {
      Stage3_print_sam(fp,stage3array[0],/*pathnum*/1,npaths,/*translocationp*/true,
		       this->chromosome_iit,queryseq,/*chimera_part*/-1,
		       this->quality_shift,this->sam_paired_p,
		       this->cigar_noncanonical_splices_p,this->sam_read_group_id);
      Stage3_print_sam(fp,stage3array[1],/*pathnum*/2,npaths,/*translocationp*/true,
		       this->chromosome_iit,queryseq,/*chimera_part*/+1,
		       this->quality_shift,this->sam_paired_p,
		       this->cigar_noncanonical_splices_p,this->sam_read_group_id);

    } else {
      for (pathnum = 1; pathnum <= effective_maxpaths; pathnum++) {
	Stage3_print_sam(fp,stage3array[pathnum-1],pathnum,npaths,/*translocationp*/false,
			 this->chromosome_iit,queryseq,/*chimera_part*/0,
			 this->quality_shift,this->sam_paired_p,
			 this->cigar_noncanonical_splices_p,this->sam_read_group_id);
      }
    }
#endif

  } else if (this->printtype == COORDS) {
    for (pathnum = 1; pathnum <= effective_maxpaths; pathnum++) {
      fprintf(fp,">");
      Sequence_print_header(fp,queryseq,this->checksump);
      Stage3_print_coordinates(fp,stage3array[pathnum-1],queryseq,
			       this->chromosome_iit,this->invertmode,querylength);
    }

  } else if (this->printtype == SPLICESITES) {
    /* Print only best path */
    if (npaths > 0) {
      Stage3_print_splicesites(fp,stage3array[0],this->chromosome_iit,queryseq);
    }

  } else if (this->printtype == MAP_GENES) {
    for (pathnum = 1; pathnum <= effective_maxpaths; pathnum++) {
      Stage3_print_iit_map(fp,stage3array[pathnum-1],pathnum,this->chromosome_iit,queryseq);
    }

  } else if (this->printtype == MAP_EXONS) {
    for (pathnum = 1; pathnum <= effective_maxpaths; pathnum++) {
      Stage3_print_iit_exon_map(fp,stage3array[pathnum-1],pathnum,this->chromosome_iit,queryseq);
    }

  } else {
    fprintf(stderr,"Unexpected printtype %d\n",this->printtype);
    abort();

  }

  return;
}

#endif


void *
Outbuffer_thread_anyorder (void *data) {
  T this = (T) data;
  Result_T result;
  Request_T request;
  
#ifdef HAVE_PTHREAD
  pthread_mutex_lock(&this->lock);
#endif

  while (this->noutput < this->ndone) {

#ifdef HAVE_PTHREAD
    while (this->head == NULL && this->noutput < this->ndone) {
      debug(fprintf(stderr,"__outbuffer_thread_anyorder waiting for result_avail_p\n"));
      pthread_cond_wait(&this->result_avail_p,&this->lock);
    }
    debug(fprintf(stderr,"__outbuffer_thread_anyorder woke up\n"));
#endif

    if (this->head) {
      /* Just take one output at a time to allow workers access to the queue */
      this->head = RRlist_pop(this->head,&request,&result);
      debug1(RRlist_dump(this->head,this->tail));

#ifdef HAVE_PTHREAD
      pthread_mutex_unlock(&this->lock);
#endif
      Outbuffer_print_result(this,result,request);
      Result_free(&result);
      Request_free(&request);
#ifdef HAVE_PTHREAD
      pthread_mutex_lock(&this->lock);
#endif

      this->noutput++;
    }
  }

#ifdef HAVE_PTHREAD
  pthread_mutex_unlock(&this->lock);
#endif

  assert(this->head == NULL);

  return (void *) NULL;
}



void *
Outbuffer_thread_ordered (void *data) {
  T this = (T) data;
  Result_T result;
  Request_T request;
  RRlist_T queue = NULL;
  int id;

#ifdef HAVE_PTHREAD
  pthread_mutex_lock(&this->lock);
#endif

  while (this->noutput < this->ndone) {
#ifdef HAVE_PTHREAD
    while (this->head == NULL && this->noutput < this->ndone) {
      pthread_cond_wait(&this->result_avail_p,&this->lock);
    }
    debug(fprintf(stderr,"__outbuffer_thread_ordered woke up\n"));
#endif

    if (this->head) {
      /* Just take one output at a time to allow workers access to the queue */
      this->head = RRlist_pop(this->head,&request,&result);

#ifdef HAVE_PTHREAD
      pthread_mutex_unlock(&this->lock);
#endif
    
      if ((id = Result_id(result)) != this->noutput) {
	/* Store in queue */
	queue = RRlist_insert(queue,id,request,result);
      } else {
	Outbuffer_print_result(this,result,request);
	Result_free(&result);
	Request_free(&request);
	this->noutput++;
	
	/* Print out rest of stored queue */
	while (queue != NULL && queue->id == this->noutput) {
	  queue = RRlist_pop_id(queue,&id,&request,&result);
	  Outbuffer_print_result(this,result,request);
	  Result_free(&result);
	  Request_free(&request);
	  this->noutput++;
	}
      }

#ifdef HAVE_PTHREAD
      pthread_mutex_lock(&this->lock);
#endif
    }
  }

#ifdef HAVE_PTHREAD
  pthread_mutex_unlock(&this->lock);
#endif

  assert(queue == NULL);

  return (void *) NULL;
}


