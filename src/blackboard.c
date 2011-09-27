static char rcsid[] = "$Id: blackboard.c,v 1.18 2007/08/28 23:24:29 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_PTHREAD
#include "blackboard.h"
#include <stdlib.h>
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>		/* Needed to define pthread_t on Solaris */
#endif
#include <pthread.h>
#include "bool.h"
#include "mem.h"
#include "reqpost.h"

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

#define T Blackboard_T
struct T {
  FILE *input;
  char **files;
  int nfiles;
  int nextchar;
  Sequence_T usersegment;
  Params_T params;

  int nworkers;
  Reqpost_T *reqposts;

  pthread_mutex_t lock;

  int inputi;			/* written and read by input thread */
  int outputi;			/* written and read by output thread */

  int ninputs;			/* written by input thread, read by output thread */
  int noutputs;			/* written and read by output thread */
  bool inputdone;		/* written by input thread, read by output thread */

  bool *input_ready;		/* written and read by input and output threads */
  bool *output_ready;		/* written and read by worker and output threads */
  int n_input_ready;		/* written and read by input and output threads */
  int n_output_ready;		/* written and read by worker and output threads */
  pthread_cond_t input_ready_p;
  pthread_cond_t output_ready_p;
};

T
Blackboard_new (FILE *input, char **files, int nfiles, int nextchar, Sequence_T usersegment, int nworkers, Params_T params) {
  T new = (T) MALLOC(sizeof(*new));
  int i;

  new->input = input;
  new->files = files;
  new->nfiles = nfiles;
  new->nextchar = nextchar;
  new->usersegment = usersegment;
  new->params = params;
  new->nworkers = nworkers;
  new->reqposts = (Reqpost_T *) CALLOC(nworkers,sizeof(Reqpost_T));
  for (i = 0; i < nworkers; i++) {
    new->reqposts[i] = Reqpost_new(new,i,params);
  }

  pthread_mutex_init(&new->lock,NULL);

  new->inputi = 0;
  new->outputi = 0;

  new->ninputs = 0;
  new->noutputs = 0;
  new->inputdone = false;

  new->input_ready = (bool *) CALLOC(nworkers,sizeof(bool));
  new->output_ready = (bool *) CALLOC(nworkers,sizeof(bool));
  for (i = 0; i < nworkers; i++) {
    new->input_ready[i] = true;
    new->output_ready[i] = false;
  }
  new->n_input_ready = nworkers;
  new->n_output_ready = 0;
  pthread_cond_init(&new->input_ready_p,NULL);
  pthread_cond_init(&new->output_ready_p,NULL);

  return new;
}


void
Blackboard_free (T *old) {
  int i;

  if (*old) {
    pthread_cond_destroy(&(*old)->input_ready_p);
    pthread_cond_destroy(&(*old)->output_ready_p);
    FREE((*old)->input_ready);
    FREE((*old)->output_ready);
    
    pthread_mutex_destroy(&(*old)->lock);
    /* Let each thread destroy its own reqpost */
    /*
    for (i = 0; i < (*old)->nworkers; i++) {
      Reqpost_free(&((*old)->reqposts[i]));
    }
    */
    FREE((*old)->reqposts);
    FREE(*old);
  }
  return;
}


FILE *
Blackboard_input (T this) {
  return this->input;
}

char **
Blackboard_files (T this) {
  return this->files;
}

int
Blackboard_nfiles (T this) {
  return this->nfiles;
}

int
Blackboard_nextchar (T this) {
  return this->nextchar;
}

Sequence_T
Blackboard_usersegment (T this) {
  return this->usersegment;
}

Params_T
Blackboard_params (T this) {
  return this->params;
}

Reqpost_T
Blackboard_get_reqpost (T this, int i) {
  return this->reqposts[i];
}


/* Executed by input thread */
void
Blackboard_put_request (T this, Request_T request) {
  debug(printf("Blackboard: got input request id %d\n",Request_id(request)));

  pthread_mutex_lock(&this->lock);
  
  while (this->n_input_ready == 0) {
    debug(printf("Blackboard: cond wait for input_ready_p\n"));
    pthread_cond_wait(&this->input_ready_p,&this->lock);
  }
  debug(printf("Blackboard: cond okay for input_ready_p\n")); /*  */
    
  while (this->input_ready[this->inputi] == false) {
    if (++this->inputi >= this->nworkers) {
      this->inputi = 0;
    }
  }
  debug(printf("Blackboard: putting request into reqpost %d\n",this->inputi));
  Reqpost_put_request(this->reqposts[this->inputi],request);
  this->input_ready[this->inputi] = false;
  this->n_input_ready -= 1;
  this->ninputs++;

  if (++this->inputi >= this->nworkers) {
    this->inputi = 0;
  }
  pthread_mutex_unlock(&this->lock);

  return;
}


/* Executed by input thread */
void
Blackboard_set_inputdone (T this) {
  pthread_mutex_lock(&this->lock);
  debug(printf("Blackboard: Setting inputdone to be true\n"));
  this->inputdone = true;
  pthread_mutex_unlock(&this->lock);
  return;
}


/* Executed by worker thread */
void
Blackboard_notify_output_ready (T this, int i) {
  pthread_mutex_lock(&this->lock);
  
  this->output_ready[i] = true;
  if (this->n_output_ready++ == 0) {
    debug(printf("Blackboard: Signaling that output is ready\n"));
    pthread_cond_signal(&this->output_ready_p);
  }
  pthread_mutex_unlock(&this->lock);
  return;
}


/* Executed by output thread */
Result_T
Blackboard_get_result (Request_T *request, T this) {
  Result_T result;

  pthread_mutex_lock(&this->lock);
  if (this->inputdone == true &&
      this->noutputs == this->ninputs) {
    pthread_mutex_unlock(&this->lock);
    return NULL;

  } else {
    while (this->n_output_ready == 0) {
      debug(printf("Blackboard: Cond wait for output_ready_p\n"));
      pthread_cond_wait(&this->output_ready_p,&this->lock);
    }
    debug(printf("Blackboard: Cond okay for output_ready_p\n"));
      
    while (this->output_ready[this->outputi] == false) {
      if (++this->outputi >= this->nworkers) {
	this->outputi = 0;
      }
    }
    result = Reqpost_get_result(&(*request),this->reqposts[this->outputi]);
    this->output_ready[this->outputi] = false;
    this->n_output_ready -= 1;
    
    this->input_ready[this->outputi] = true;
    if (this->n_input_ready++ == 0) {
      debug(printf("Blackboard: Signaling that input is ready\n"));
      pthread_cond_signal(&this->input_ready_p);
    }

    this->noutputs++;
    if (++this->outputi >= this->nworkers) {
      this->outputi = 0;
    }

    pthread_mutex_unlock(&this->lock);

    debug(printf("Blackboard: Returning result id %d\n",Result_id(result)));
    return result;
  }
}

#endif /* HAVE_PTHREAD */
