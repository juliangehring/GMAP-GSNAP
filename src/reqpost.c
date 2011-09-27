static char rcsid[] = "$Id: reqpost.c,v 1.20 2008-02-28 18:12:06 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_PTHREAD
#include "reqpost.h"
#include <stdlib.h>
#ifdef HAVE_SYS_TYPES_H
# include <sys/types.h>		/* Needed to define pthread_t on Solaris */
#endif
#include <pthread.h>
#include "bool.h"
#include "mem.h"

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


#define T Reqpost_T
struct T {
  int id;
  pthread_mutex_t lock;

  Request_T request;
  Result_T result;
  bool resultp;

  pthread_cond_t request_ready_p;

  Blackboard_T blackboard;
};

T
Reqpost_new (Blackboard_T blackboard, int id) {
  T new = (T) MALLOC(sizeof(*new));

  new->id = id;
  pthread_mutex_init(&new->lock,NULL);

  new->request = NULL;
  new->result = NULL;
  new->resultp = false;

  pthread_cond_init(&new->request_ready_p,NULL);

  new->blackboard = blackboard;
  return new;
}

void
Reqpost_free (T *old) {
  pthread_cond_destroy(&(*old)->request_ready_p);
  pthread_mutex_destroy(&(*old)->lock);
  FREE(*old);
  return;
}

int
Reqpost_id (T this) {
  return this->id;
}


/* Called by input thread */
void
Reqpost_put_request (T this, Request_T request) {
  pthread_mutex_lock(&this->lock);

  this->request = request;
  pthread_cond_signal(&this->request_ready_p);

  pthread_mutex_unlock(&this->lock);
  return;
}

/* Called by worker thread */
Request_T
Reqpost_get_request (T this) {
  Request_T request;

  pthread_mutex_lock(&this->lock);
  while (this->request == NULL || this->resultp == true) {
    debug(printf("Reqpost %d: Cond wait for request_ready_p\n",this->id));
    pthread_cond_wait(&this->request_ready_p,&this->lock);
  }
  debug(printf("Reqpost %d: Cond okay for request_ready_p\n",this->id));
  request = this->request;
  pthread_mutex_unlock(&this->lock);
  return request;
}


/* Called by worker thread */
void
Reqpost_put_result (T this, Result_T result) {
  pthread_mutex_lock(&this->lock);
  this->result = result;
  this->resultp = true;
  pthread_mutex_unlock(&this->lock);
  Blackboard_notify_output_ready(this->blackboard,this->id);
  return;
}

/* Called by output thread */
Result_T
Reqpost_get_result (Request_T *request, T this) {
  Result_T result;

  pthread_mutex_lock(&this->lock);
  result = this->result;
  *request = this->request;

  this->request = NULL;
  this->result = NULL;
  this->resultp = false;
  pthread_mutex_unlock(&this->lock);

  return result;
}

#endif /* HAVE_PTHREAD */
