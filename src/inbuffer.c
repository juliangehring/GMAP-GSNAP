static char rcsid[] = "$Id: inbuffer.c 49442 2011-10-08 01:32:38Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "inbuffer.h"
#include <stdio.h>
#include <stdlib.h>

#ifdef HAVE_PTHREAD
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>		/* Needed to define pthread_t on Solaris */
#endif
#include <pthread.h>
#endif

#include "mem.h"

#ifdef GSNAP
#include "shortread.h"
#endif


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


#define T Inbuffer_T

struct T {
  Outbuffer_T outbuffer;
  bool filter_if_both_p;

  FILE *input;
#ifdef GSNAP
  FILE *input2;
#endif

#ifdef HAVE_ZLIB
  gzFile gzipped;
  gzFile gzipped2;
#else
  void *gzipped;
  void *gzipped2;
#endif

  char **files;
  int nfiles;
  int nextchar;

#ifdef GSNAP
  Gobyreader_T gobyreader;
  bool fastq_format_p;
  bool creads_format_p;
  int barcode_length;
  bool invert_first_p;
  bool invert_second_p;
  bool chop_primers_p;
#else
  bool maponlyp;
#endif

  int part_modulus;
  int part_interval;

  unsigned int nspaces;
  unsigned int maxchars;

#ifdef HAVE_PTHREAD
  pthread_mutex_t lock;
#endif

  Request_T *buffer;
  int ptr;
  int nleft;
  int inputid;
  int requestid;
};


#ifndef GSNAP
T
Inbuffer_cmdline (char *contents, int length) {
  T new = (T) MALLOC(sizeof(*new));

  new->input = (FILE *) NULL;
  new->files = (char **) NULL;
  new->nfiles = 0;
  new->nextchar = '\0';
  new->maponlyp = false;

  new->part_modulus = 0;
  new->part_interval = 1;

  new->nspaces = 0;
  /* new->maxchars = maxchars; */

  new->buffer = (Request_T *) CALLOC(1,sizeof(Request_T));

  new->ptr = 0;
  new->nleft = 1;
  new->inputid = 0;
  new->requestid = 0;

  new->buffer[0] = Request_new(new->requestid++,Sequence_genomic_new(contents,length,/*copyp*/true));

#ifdef HAVE_PTHREAD
  pthread_mutex_init(&new->lock,NULL);
#endif

  return new;
}
#endif



T
Inbuffer_new (int nextchar, FILE *input,
#ifdef GSNAP
	      FILE *input2,
#ifdef HAVE_ZLIB
	      gzFile gzipped, gzFile gzipped2,
#endif
#ifdef HAVE_GOBY
	      Gobyreader_T gobyreader,
#endif
#endif
	      char **files, int nfiles,
#ifdef GSNAP
	      bool fastq_format_p, bool creads_format_p,
	      int barcode_length, bool invert_first_p, bool invert_second_p,
	      bool chop_primers_p,
#else
	      bool maponlyp,
#endif
	      unsigned int nspaces, unsigned int maxchars, int part_interval, int part_modulus,
	      bool filter_if_both_p) {

  T new = (T) MALLOC(sizeof(*new));

  new->filter_if_both_p = filter_if_both_p;
  new->input = input;
#ifdef GSNAP
  new->input2 = input2;
#ifdef HAVE_ZLIB
  new->gzipped = gzipped;
  new->gzipped2 = gzipped2;
#else
  new->gzipped = (void *) NULL;
  new->gzipped2 = (void *) NULL;
#endif
#ifdef HAVE_GOBY
  new->gobyreader = gobyreader;
#endif
#endif

  new->files = files;
  new->nfiles = nfiles;
  new->nextchar = nextchar;

#ifdef GSNAP
  new->fastq_format_p = fastq_format_p;
  new->creads_format_p = creads_format_p;
  new->barcode_length = barcode_length;
  new->invert_first_p = invert_first_p;
  new->invert_second_p = invert_second_p;
  new->chop_primers_p = chop_primers_p;
#if 0
  if (chop_primers_p == true) {
    Shortread_dynprog_init(MAX_READLENGTH);
  }
#endif
#else
  new->maponlyp = maponlyp;
#endif

  new->part_modulus = part_modulus;
  new->part_interval = part_interval;

  new->nspaces = nspaces;
  new->maxchars = maxchars;

#ifdef HAVE_PTHREAD
  pthread_mutex_init(&new->lock,NULL);
#endif

  new->buffer = (Request_T *) CALLOC(nspaces,sizeof(Request_T));
  new->ptr = 0;
  new->nleft = 0;
  new->inputid = 0;
  new->requestid = 0;

  return new;
}

void
Inbuffer_set_outbuffer (T this, Outbuffer_T outbuffer) {
  this->outbuffer = outbuffer;
  return;
}


void
Inbuffer_free (T *old) {
  if (*old) {
    /* No need to close input, since done by Shortread and Sequence read procedures */

    FREE((*old)->buffer);
    
#ifdef HAVE_PTHREAD
    pthread_mutex_destroy(&(*old)->lock);
#endif

    FREE(*old);
  }
  return;
}


#ifdef GSNAP

/* Returns number of requests read */
static unsigned int
fill_buffer (T this) {
  unsigned int nread = 0;
  unsigned int nchars = 0U;
  Shortread_T queryseq1, queryseq2;

  if (this->fastq_format_p == true) {
    if (this->gzipped != NULL) {
#ifdef HAVE_ZLIB    
      /* FASTQ input, gzipped */
      while (nread < this->nspaces &&
#if 0
	     nchars < this->maxchars &&
#endif
	     (queryseq1 = Shortread_read_fastq_shortreads_gzip(&this->nextchar,&queryseq2,&this->gzipped,&this->gzipped2,
							       &this->files,&this->nfiles,
							       this->barcode_length,this->invert_first_p,this->invert_second_p)) != NULL) {
	if (this->inputid % this->part_interval != this->part_modulus) {
	  Shortread_free(&queryseq1);
	  if (queryseq2 != NULL) {
	    Shortread_free(&queryseq2);
	  }

	} else if (this->filter_if_both_p == true &&
		   Shortread_filterp(queryseq1) == true && (queryseq2 == NULL || Shortread_filterp(queryseq2) == true)) {
	  Shortread_free(&queryseq1);
	  if (queryseq2 != NULL) {
	    Shortread_free(&queryseq2);
	  }

	} else if (this->filter_if_both_p == false &&
		   (Shortread_filterp(queryseq1) == true || (queryseq2 != NULL && Shortread_filterp(queryseq2) == true))) {
	  Shortread_free(&queryseq1);
	  if (queryseq2 != NULL) {
	    Shortread_free(&queryseq2);
	  }

	} else {
	  this->buffer[nread++] = Request_new(this->requestid++,queryseq1,queryseq2);
	  nchars += Shortread_fulllength(queryseq1);
	  if (queryseq2 != NULL) {
	    nchars += Shortread_fulllength(queryseq2);
	  }
	}
	this->inputid++;
      }
#endif

    } else {
      /* FASTQ input, text */
      while (nread < this->nspaces &&
#if 0
	     nchars < this->maxchars &&
#endif
	     (queryseq1 = Shortread_read_fastq_shortreads(&this->nextchar,&queryseq2,&this->input,&this->input2,
							  &this->files,&this->nfiles,
							  this->barcode_length,this->invert_first_p,this->invert_second_p)) != NULL) {
	if (this->inputid % this->part_interval != this->part_modulus) {
	  Shortread_free(&queryseq1);
	  if (queryseq2 != NULL) {
	    Shortread_free(&queryseq2);
	  }

	} else if (this->filter_if_both_p == true &&
		   Shortread_filterp(queryseq1) == true && (queryseq2 == NULL || Shortread_filterp(queryseq2) == true)) {
	  Shortread_free(&queryseq1);
	  if (queryseq2 != NULL) {
	    Shortread_free(&queryseq2);
	  }

	} else if (this->filter_if_both_p == false &&
		   (Shortread_filterp(queryseq1) == true || (queryseq2 != NULL && Shortread_filterp(queryseq2) == true))) {
	  Shortread_free(&queryseq1);
	  if (queryseq2 != NULL) {
	    Shortread_free(&queryseq2);
	  }

	} else {
	  this->buffer[nread++] = Request_new(this->requestid++,queryseq1,queryseq2);
	  nchars += Shortread_fulllength(queryseq1);
	  if (queryseq2 != NULL) {
	    nchars += Shortread_fulllength(queryseq2);
	  }
	}
	this->inputid++;
      }
    }

  } else if (this->creads_format_p == true) {
#ifdef HAVE_GOBY
    /* GOBY input */
    while (nread < this->nspaces &&
#if 0
	   nchars < this->maxchars &&
#endif
	   (queryseq1 = Goby_read(&queryseq2,this->gobyreader,this->barcode_length,
				  this->invert_first_p,this->invert_second_p)) != NULL) {
      if (this->inputid % this->part_interval != this->part_modulus) {
	Shortread_free(&queryseq1);
	if (queryseq2 != NULL) {
	  Shortread_free(&queryseq2);
	}

      } else if (this->filter_if_both_p == true &&
		 Shortread_filterp(queryseq1) == true && (queryseq2 == NULL || Shortread_filterp(queryseq2) == true)) {
	Shortread_free(&queryseq1);
	if (queryseq2 != NULL) {
	  Shortread_free(&queryseq2);
	}

      } else if (this->filter_if_both_p == false &&
		 (Shortread_filterp(queryseq1) == true || (queryseq2 != NULL && Shortread_filterp(queryseq2) == true))) {
	Shortread_free(&queryseq1);
	if (queryseq2 != NULL) {
	  Shortread_free(&queryseq2);
	}

      } else {
	this->buffer[nread++] = Request_new(this->requestid++,queryseq1,queryseq2);
	nchars += Shortread_fulllength(queryseq1);
	if (queryseq2 != NULL) {
	  nchars += Shortread_fulllength(queryseq2);
	}
      }
      this->inputid++;
    }
#endif

  } else {
    if (this->gzipped != NULL) {
#ifdef HAVE_ZLIB    
      /* FASTA input, gzipped */
      while (nread < this->nspaces &&
#if 0
	     nchars < this->maxchars &&
#endif
	     (queryseq1 = Shortread_read_fasta_shortreads_gzip(&this->nextchar,&queryseq2,&this->gzipped,&this->files,&this->nfiles,
							       this->barcode_length,this->invert_first_p,this->invert_second_p)) != NULL) {
	if (this->inputid % this->part_interval != this->part_modulus) {
	  Shortread_free(&queryseq1);
	  if (queryseq2 != NULL) {
	    Shortread_free(&queryseq2);
	  }

	} else if (this->filter_if_both_p == true &&
		   Shortread_filterp(queryseq1) == true && (queryseq2 == NULL || Shortread_filterp(queryseq2) == true)) {
	  Shortread_free(&queryseq1);
	  if (queryseq2 != NULL) {
	    Shortread_free(&queryseq2);
	  }

	} else if (this->filter_if_both_p == false &&
		   (Shortread_filterp(queryseq1) == true || (queryseq2 != NULL && Shortread_filterp(queryseq2) == true))) {
	  Shortread_free(&queryseq1);
	  if (queryseq2 != NULL) {
	    Shortread_free(&queryseq2);
	  }

	} else {
	  this->buffer[nread++] = Request_new(this->requestid++,queryseq1,queryseq2);
	  nchars += Shortread_fulllength(queryseq1);
	  if (queryseq2 != NULL) {
	    nchars += Shortread_fulllength(queryseq2);
	  }
	}
	this->inputid++;
      }
#endif

    } else {
      /* FASTA input, text */
      while (nread < this->nspaces &&
#if 0
	     nchars < this->maxchars &&
#endif
	     (queryseq1 = Shortread_read_fasta_shortreads(&this->nextchar,&queryseq2,&this->input,&this->files,&this->nfiles,
							  this->barcode_length,this->invert_first_p,this->invert_second_p)) != NULL) {
	if (this->inputid % this->part_interval != this->part_modulus) {
	  Shortread_free(&queryseq1);
	  if (queryseq2 != NULL) {
	    Shortread_free(&queryseq2);
	  }

	} else if (this->filter_if_both_p == true &&
		   Shortread_filterp(queryseq1) == true && (queryseq2 == NULL || Shortread_filterp(queryseq2) == true)) {
	  Shortread_free(&queryseq1);
	  if (queryseq2 != NULL) {
	    Shortread_free(&queryseq2);
	  }

	} else if (this->filter_if_both_p == false &&
		   (Shortread_filterp(queryseq1) == true || (queryseq2 != NULL && Shortread_filterp(queryseq2) == true))) {
	  Shortread_free(&queryseq1);
	  if (queryseq2 != NULL) {
	    Shortread_free(&queryseq2);
	  }

	} else {
	  debug(printf("inbuffer creating request %d\n",this->requestid));
	  this->buffer[nread++] = Request_new(this->requestid++,queryseq1,queryseq2);
	  nchars += Shortread_fulllength(queryseq1);
	  if (queryseq2 != NULL) {
	    nchars += Shortread_fulllength(queryseq2);
	  }
	}
	this->inputid++;
      }
    }
  }

  this->nleft = nread;
  this->ptr = 0;

  return nread;
}

#else
	 
/* Returns number of requests read */
static unsigned int
fill_buffer (T this) {
  unsigned int nread = 0;
  unsigned int nchars = 0U;
  Sequence_T queryseq;

  while (nread < this->nspaces &&
#if 0
	 nchars < this->maxchars &&
#endif
	 (queryseq = Sequence_read_multifile(&this->nextchar,&this->input,
					     &this->files,&this->nfiles,this->maponlyp)) != NULL) {
    if (this->inputid % this->part_interval != this->part_modulus) {
      Sequence_free(&queryseq);
    } else {
      debug(printf("inbuffer creating request %d\n",this->requestid));
      this->buffer[nread++] = Request_new(this->requestid++,queryseq);
      nchars += Sequence_fulllength(queryseq);
    }
    this->inputid++;
  }

  this->nleft = nread;
  this->ptr = 0;

  return nread;
}

#endif


/* No need to lock, since only main thread calls */
/* Returns nread to give to Outbuffer_new */
unsigned int
Inbuffer_fill_init (T this) {
  unsigned int nread;

  debug(printf("inbuffer filling initially\n"));
  nread = fill_buffer(this);
  debug(printf("inbuffer read %d sequences\n",nread));

  return nread;
}
  


Request_T
Inbuffer_get_request (T this) {
  Request_T request;
  unsigned int nread;

#ifdef HAVE_PTHREAD
  pthread_mutex_lock(&this->lock);
#endif
  
  if (this->nleft > 0) {
    request = this->buffer[this->ptr++];
    this->nleft -= 1;

  } else {
    debug(printf("inbuffer filling\n"));
    nread = fill_buffer(this);
    Outbuffer_add_nread(this->outbuffer,nread);
    debug(printf("inbuffer read %d sequences\n",nread));
    
    if (nread == 0) {
      /* Still empty */
      request = NULL;
    } else {
      request = this->buffer[this->ptr++];
      this->nleft -= 1;
    }
  }

#ifdef HAVE_PTHREAD
  pthread_mutex_unlock(&this->lock);
#endif

  return request;
}


