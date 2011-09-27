static char rcsid[] = "$Id: substring.c,v 1.6 2010/03/09 20:25:23 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "substring.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "mem.h"
#include "stage1hr.h"		/* For MAX_QUERYLENGTH */
#include "maxent.h"
#include "listdef.h"
#include "list.h"
#include "complement.h"


#define MATCH_SCORE 1
#define MISMATCH_SCORE -1

#define TRANSLOCATION_TEXT "splice_translocation"
#define INVERSION_TEXT "splice_inversion"
#define SCRAMBLE_TEXT "splice_scramble"


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* mark_mismatches */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

/* Substring_new */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

/* splice site probs */
#ifdef DEBUG4
#define debug4(x) x
#else
#define debug4(x)
#endif


/* trimming */
#ifdef DEBUG8
#define debug8(x) x
#else
#define debug8(x)
#endif



static bool print_nsnpdiffs_p = false;
static bool print_ncolordiffs_p = false;
static bool print_snplabels_p = false;

void
Substring_print_nsnpdiffs (bool labelsp) {
  print_nsnpdiffs_p = true;
  print_snplabels_p = labelsp;
  return;
}

void
Substring_print_ncolordiffs () {
  print_ncolordiffs_p = true;
  return;
}


static char complCode[128] = COMPLEMENT_LC;

static char *
make_complement_buffered (char *complement, char *sequence, unsigned int length) {
  int i, j;

  /* complement = (char *) CALLOC(length+1,sizeof(char)); */
  for (i = length-1, j = 0; i >= 0; i--, j++) {
    complement[j] = complCode[(int) sequence[i]];
  }
  complement[length] = '\0';
  return complement;
}


struct Substring_T {
  int nmismatches;
  int nsnpdiffs;
  int ncolordiffs;
  int trim_left;
  int trim_right;
  Chrnum_T chrnum;
  Genomicpos_T chroffset;

  Genomicpos_T genomicstart;	/* For region corresponding to querylength */
  Genomicpos_T genomicend;

  int querystart;		/* For part that aligns to genome */
  int queryend;
  int querylength;

  Genomicpos_T alignstart;	/* For part that aligns to genome, including part that is trimmed */
  Genomicpos_T alignend;

  Genomicpos_T alignstart_trim;	/* For part that aligns to genome, excluding part that is trimmed */
  Genomicpos_T alignend_trim;

  int genomiclength;
  bool plusp;
  char *genomicdir;		/* In same direction as query */

  /* for splices */
  bool chimera_sensep;
  bool chimera_knownp;
  int chimera_pos;
  double chimera_prob;
  Genomicpos_T chimera_modelpos;
};


static void
fill_w_dashes (char *string, int start, int end) {
  int i;

  for (i = start; i < end; i++) {
    string[i] = '-';
  }
  return;
}

static int
mark_mismatches (char *genome, char *query, int start, int end) {
  int nmismatches = 0, i;

  debug1(printf("query:  %s\n",query));
  debug1(printf("genome: %s\n",genome));
  debug1(printf("mark:   "));

  debug1(
	 for (i = 0; i < start; i++) {
	   printf(" ");
	 }
	 );

  for (i = start; i < end; i++) {
    if (query[i] != genome[i]) {
      genome[i] = (char) tolower(genome[i]);
      nmismatches++;
      debug1(printf("x"));
    } else {
      debug1(printf("*"));
    }
  }
  debug1(printf("\n"));
  debug1(printf("%d mismatches\n\n",nmismatches));

  return nmismatches;
}


static bool
dibase_mismatch_p (char leftnt, char leftcolor, char nt, char rightcolor, char rightnt) {
  char left_predict, right_predict;

  switch (toupper(leftnt)) {
  case 'A':
    switch (leftcolor) {
    case '0': left_predict = 'A'; break;
    case '1': left_predict = 'C'; break;
    case '2': left_predict = 'G'; break;
    case '3': left_predict = 'T'; break;
    }
    break;
  case 'C':
    switch (leftcolor) {
    case '0': left_predict = 'C'; break;
    case '1': left_predict = 'A'; break;
    case '2': left_predict = 'T'; break;
    case '3': left_predict = 'G'; break;
    }
    break;
  case 'G':
    switch (leftcolor) {
    case '0': left_predict = 'G'; break;
    case '1': left_predict = 'T'; break;
    case '2': left_predict = 'A'; break;
    case '3': left_predict = 'C'; break;
    }
    break;
  case 'T':
    switch (leftcolor) {
    case '0': left_predict = 'T'; break;
    case '1': left_predict = 'G'; break;
    case '2': left_predict = 'C'; break;
    case '3': left_predict = 'A'; break;
    }
    break;
  }

  switch (toupper(rightnt)) {
  case 'A':
    switch (rightcolor) {
    case '0': right_predict = 'A'; break;
    case '1': right_predict = 'C'; break;
    case '2': right_predict = 'G'; break;
    case '3': right_predict = 'T'; break;
    }
    break;
  case 'C':
    switch (rightcolor) {
    case '0': right_predict = 'C'; break;
    case '1': right_predict = 'A'; break;
    case '2': right_predict = 'T'; break;
    case '3': right_predict = 'G'; break;
    }
    break;
  case 'G':
    switch (rightcolor) {
    case '0': right_predict = 'G'; break;
    case '1': right_predict = 'T'; break;
    case '2': right_predict = 'A'; break;
    case '3': right_predict = 'C'; break;
    }
    break;
  case 'T':
    switch (rightcolor) {
    case '0': right_predict = 'T'; break;
    case '1': right_predict = 'G'; break;
    case '2': right_predict = 'C'; break;
    case '3': right_predict = 'A'; break;
    }
    break;
  }

  if (left_predict != right_predict) {
    /* printf("left_predict %c != right_predict %c\n",left_predict,right_predict); */
    return false;
  } else if (left_predict == toupper(nt)) {
    return false;
  } else {
    return true;
  }
}


/* query is numbers */
static int
mark_mismatches_dibase (char *genome, char *query, int start, int end, bool plusp) {
  int nmismatches = 0, i;

  debug1(printf("query:  %s\n",query));
  debug1(printf("genome: %s\n",genome));
  debug1(printf("mark:   "));

  if (plusp) {
    /* querypos start */
    debug1(printf("*"));
    if (dibase_mismatch_p(/*leftnt*/genome[start],/*leftcolor*/'0',genome[start],query[start],genome[start+1])) {
      genome[start] = (char) tolower(genome[start]);
      nmismatches++;
    }

    for (i = start+1; i < end-1; i++) {
      debug1(printf("*"));

      if (dibase_mismatch_p(genome[i-1],query[i-1],genome[i],query[i],genome[i+1])) {
	/* printf("Got a dibase mismatch at %d with %c %c %c %c %c\n",
	   i,genome[i-1],query[i-1],genome[i],query[i],genome[i+1]); */
	genome[i] = (char) tolower(genome[i]);
	nmismatches++;
      }
    }

    /* querypos (end-1) */
    if (dibase_mismatch_p(genome[end-2],query[end-2],genome[end-1],
			  /*rightcolor*/'0',/*rightnt*/genome[end-1])) {
      genome[end-1] = (char) tolower(genome[end-1]);
      nmismatches++;
    }

  } else {
    /* querypos start */
    debug1(printf("*"));
    if (dibase_mismatch_p(/*leftnt*/genome[start],/*leftcolor*/'0',genome[start],query[start+1],genome[start+1])) {
      genome[start] = (char) tolower(genome[start]);
      nmismatches++;
    }

    for (i = start+1; i < end-1; i++) {
      debug1(printf("*"));

      if (dibase_mismatch_p(genome[i-1],query[i],genome[i],query[i+1],genome[i+1])) {
	/* printf("Got a dibase mismatch at %d with %c %c %c %c %c\n",
	   i,genome[i-1],query[i],genome[i],query[i+1],genome[i+1]); */
	genome[i] = (char) tolower(genome[i]);
	nmismatches++;
      }
    }

    /* querypos (end-1) */
    if (dibase_mismatch_p(genome[end-2],query[end-1],genome[end-1],
			  /*rightcolor*/'0',/*rightnt*/genome[end-1])) {
      genome[end-1] = (char) tolower(genome[end-1]);
      nmismatches++;
    }
  }

  debug1(printf("\n"));
  debug1(printf("%d mismatches\n\n",nmismatches));

  return nmismatches;
}



static int
mark_mismatches_met (int *nmetdiffs, char *gbuffer, char *query, int start, int end, bool plusp) {
  int nmismatches = 0, i;
  int nunknowns;
  
  debug1(printf("query:  %s\n",query));
  debug1(printf("genome: %s\n",gbuffer));
  debug1(printf("count:  "));

  nunknowns = 0;
  *nmetdiffs = 0;

  if (plusp) {
    for (i = start; i < end; i++) {
      if (gbuffer[i] == 'C' && query[i] == 'T') {
	debug1(printf("."));
	gbuffer[i] = '.';
	(*nmetdiffs)++;
      } else if (query[i] != gbuffer[i]) {
	debug1(printf("x"));
	if (gbuffer[i] == OUTOFBOUNDS) {
	  abort();
	  nunknowns++;
	} else {
	  nmismatches++;
	}
	gbuffer[i] = (char) tolower(gbuffer[i]);
      } else {
	debug1(printf("*"));
      }
    }

  } else {
    for (i = start; i < end; i++) {
      if (gbuffer[i] == 'G' && query[i] == 'A') {
	debug1(printf("."));
	gbuffer[i] = '.';
	(*nmetdiffs)++;
      } else if (query[i] != gbuffer[i]) {
	debug1(printf("x"));
	if (gbuffer[i] == OUTOFBOUNDS) {
	  abort();
	  nunknowns++;
	} else {
	  nmismatches++;
	}
	gbuffer[i] = (char) tolower(gbuffer[i]);
      } else {
	debug1(printf("*"));
      }
    }
  }

  debug1(printf("\n"));

  return nmismatches;
}


/* querystart (typically 0) and query3 (typically querylength) are exclusive */
/* sequences may have had lower case characters marked */
static int
trim_left_end (char *genomicdir, char *query, int querystart, int queryend) {
  int bestscore, score;
  int trim5, alignlength, pos;
  char *p, *q;

  alignlength = queryend - querystart;

  p = &(query[querystart]);
  q = &(genomicdir[querystart]);

  bestscore = 0;
  score = 0;
  trim5 = 0;
  for (pos = alignlength-1; pos >= 0; pos--) {
    if (toupper(p[pos]) == toupper(q[pos])) {
      score += MATCH_SCORE;
    } else {
      score += MISMATCH_SCORE;
    }
    if (score > bestscore) {
      bestscore = score;
      trim5 = pos;
    }
    debug8(printf("Trim left pos %d, score %d, trim5 %d\n",pos,score,trim5));
  }

  debug8({
      printf("Trim left called with querystart %d and queryend %d\n",querystart,queryend);
      printf("At query ->: %.*s\n",alignlength,&(query[querystart]));
      printf("At genome->: %.*s\n",alignlength,&(genomicdir[querystart]));
      printf("trim %02d  ->: ",trim5);
      for (pos = 0; pos < trim5; pos++) {
	printf(" ");
      }
      for ( ; pos < alignlength; pos++) {
	printf("*");
      }
      printf("\n");
    });

  return trim5;
}



/* querystart (typically 0) and queryend (typically querylength) are exclusive */
/* sequences may have had lower case characters marked */
static int
trim_right_end (char *genomicdir, char *query, int querystart, int queryend) {
  int bestscore, score;
  int trim3, alignlength, pos;
  char *p, *q;

  alignlength = queryend - querystart;

  p = &(query[querystart]);
  q = &(genomicdir[querystart]);

  bestscore = 0;
  score = 0;
  trim3 = 0;
  for (pos = 0; pos < alignlength; pos++) {
    if (toupper(p[pos]) == toupper(q[pos])) {
      score += MATCH_SCORE;
    } else {
      score += MISMATCH_SCORE;
    }
    if (score > bestscore) {
      bestscore = score;
      trim3 = alignlength - pos - 1;
    }
    debug8(printf("Trim right pos %d, score %d, trim3 %d\n",pos,score,trim3));
  }

  debug8({
      printf("Trim right called with querystart %d and queryend %d\n",querystart,queryend);
      printf("At query ->: %.*s\n",alignlength,&(query[querystart]));
      printf("At genome->: %.*s\n",alignlength,&(genomicdir[querystart]));
      printf("trim %02d  ->: ",trim3);
      for (pos = 0; pos < alignlength - trim3; pos++) {
	printf("*");
      }
      for ( ; pos < alignlength; pos++) {
	printf(" ");
      }
      printf("\n");
    });

  return trim3;
}




Substring_T
Substring_new (int nmismatches, int ncolordiffs, Chrnum_T chrnum, Genomicpos_T chroffset,
	       Genomicpos_T genomicstart, Genomicpos_T genomicend,
	       int querystart, int queryend, int querylength,
	       Genomicpos_T alignstart, Genomicpos_T alignend, int genomiclength,
	       int extraleft, int extraright, char *genomicseg, char *query,
	       bool plusp, bool trim_ends_p, bool dibasep, bool cmetp) {
  Substring_T new = (Substring_T) MALLOC(sizeof(*new));
  char *genomicsegrc, *genomicdir;
  int alignoffset;
  int nmismatches_all, nmetdiffs, i, j, k;

  new->nmismatches = nmismatches;
  new->ncolordiffs = ncolordiffs;
  new->chrnum = chrnum;
  new->chroffset = chroffset;
  new->genomicstart = genomicstart;
  new->genomicend = genomicend;

  new->querystart = querystart;
  new->queryend = queryend;
  new->querylength = querylength;

  new->alignstart = new->alignstart_trim = alignstart;
  new->alignend = new->alignend_trim = alignend;

  new->genomiclength = genomiclength;
  new->plusp = plusp;

  new->trim_left = 0;
  new->trim_right = 0;
  if (genomicseg == NULL) {
    new->genomicdir = (char *) NULL;
    new->nsnpdiffs = 0;
  } else {
    if (plusp == true) {
      genomicdir = genomicseg;
      alignoffset = alignstart - genomicstart;
    } else {
      genomicsegrc = (char *) CALLOC(genomiclength+1,sizeof(char));
      genomicdir = make_complement_buffered(genomicsegrc,genomicseg,genomiclength);
      alignoffset = genomicstart - alignstart;
    }

    debug2(printf("querylength is %d, alignstart is %u, alignend is %u, genomicstart is %u, genomicend is %u, alignoffset is %d, alignlength is %d-%d=%d\n",
		  querylength,alignstart,alignend,genomicstart,genomicend,alignoffset,queryend,querystart,queryend-querystart));

    debug2(printf("q:  %s\n",query));

    new->genomicdir = (char *) CALLOC(querylength+1,sizeof(char));
    new->genomicdir[querylength] = '\0';

    fill_w_dashes(new->genomicdir,0,querystart);
    debug2(printf("g1: %s\n",new->genomicdir));
    strncpy(&(new->genomicdir[querystart]),&(genomicdir[alignoffset]),queryend-querystart);
    debug2(printf("g1: %s\n",new->genomicdir));
    fill_w_dashes(new->genomicdir,queryend,querylength);
    debug2(printf("g1: %s\n",new->genomicdir));

    for (k = 0, i = querystart-1, j = alignoffset-1; k < extraleft; k++, i--, j--) {
      new->genomicdir[i] = (char) tolower(genomicdir[j]);
    }
    for (k = 0, i = queryend, j = alignoffset+queryend-querystart; k < extraright; k++, i++, j++) {
      new->genomicdir[i] = (char) tolower(genomicdir[j]);
    }
    debug2(printf("g1: %s\n",new->genomicdir));

    if (plusp == false) {
      FREE(genomicsegrc);
    }

    if (dibasep) {
      nmismatches_all = mark_mismatches_dibase(new->genomicdir,query,querystart,queryend,plusp);
      new->nsnpdiffs = nmismatches_all - nmismatches;
    } else if (cmetp) {
      nmismatches_all = mark_mismatches_met(&nmetdiffs,new->genomicdir,query,querystart,queryend,/*plusp*/true);
      new->nsnpdiffs = nmetdiffs;
    } else {
      nmismatches_all = mark_mismatches(new->genomicdir,query,querystart,queryend);
      new->nsnpdiffs = nmismatches_all - nmismatches;
    }

    if (trim_ends_p == true && querystart == 0) {
      new->trim_left = trim_left_end(new->genomicdir,query,querystart,queryend);
      new->querystart = new->trim_left;
      if (plusp == true) {
	new->alignstart_trim += new->trim_left;
      } else {
	new->alignstart_trim -= new->trim_left;
      }

    }
    if (trim_ends_p == true && queryend == querylength) {
      new->trim_right = trim_right_end(new->genomicdir,query,querystart,queryend);
      new->queryend = querylength - new->trim_right;
      if (plusp == true) {
	new->alignend_trim -= new->trim_right;
      } else {
	new->alignend_trim += new->trim_right;
      }
    }

  }

  return new;
}


void
Substring_free (Substring_T *old) {
  if ((*old)->genomicdir != NULL) {
    FREE((*old)->genomicdir);
  }
  FREE(*old);
  return;
}

bool
Substring_plusp (Substring_T this) {
  return this->plusp;
}

int
Substring_nmismatches (Substring_T this) {
  return this->nmismatches;
}

int
Substring_ncolordiffs (Substring_T this) {
  return this->ncolordiffs;
}

int
Substring_trim_left (Substring_T this) {
  return this->trim_left;
}

int
Substring_trim_right (Substring_T this) {
  return this->trim_right;
}

int
Substring_match_length (Substring_T this) {
  return this->queryend - this->querystart; /* Values are not exclusive */
}


Chrnum_T
Substring_chrnum (Substring_T this) {
  return this->chrnum;
}

Genomicpos_T
Substring_chroffset (Substring_T this) {
  return this->chroffset;
}

Genomicpos_T
Substring_genomicstart (Substring_T this) {
  return this->genomicstart;
}

Genomicpos_T
Substring_genomicend (Substring_T this) {
  return this->genomicend;
}

Genomicpos_T
Substring_genomiclength (Substring_T this) {
  return this->genomiclength;
}

double
Substring_chimera_prob (Substring_T this) {
  return this->chimera_prob;
}

int
Substring_chimera_pos (Substring_T this) {
  return this->chimera_pos;
}

bool
Substring_chimera_knownp (Substring_T this) {
  return this->chimera_knownp;
}

bool
Substring_chimera_sensep (Substring_T this) {
  return this->chimera_sensep;
}




Substring_T
Substring_copy (Substring_T old) {
  Substring_T new;

  if (old == NULL) {
    return NULL;
  } else {
    new = (Substring_T) MALLOC(sizeof(*new));

    new->nmismatches = old->nmismatches;
    new->ncolordiffs = old->ncolordiffs;
    new->chrnum = old->chrnum;
    new->chroffset = old->chroffset;
    new->genomicstart = old->genomicstart;
    new->genomicend = old->genomicend;

    new->querystart = old->querystart;
    new->queryend = old->queryend;
    new->querylength = old->querylength;

    new->alignstart = old->alignstart;
    new->alignend = old->alignend;
    new->alignstart_trim = old->alignstart_trim;
    new->alignend_trim = old->alignend_trim;
    new->genomiclength = old->genomiclength;
    new->plusp = old->plusp;

    if (old->genomicdir == NULL) {
      new->genomicdir = (char *) NULL;
    } else {
      new->genomicdir = (char *) CALLOC(strlen(old->genomicdir)+1,sizeof(char));
      strcpy(new->genomicdir,old->genomicdir);
    }

    new->nsnpdiffs = old->nsnpdiffs;

    new->trim_left = old->trim_left;
    new->trim_right = old->trim_right;

    new->chimera_sensep = old->chimera_sensep;
    new->chimera_knownp = old->chimera_knownp;
    new->chimera_modelpos = old->chimera_modelpos;
    new->chimera_pos = old->chimera_pos;
    new->chimera_prob = old->chimera_prob;

    return new;
  }
}



Substring_T
Substring_new_donor (int donor_pos, int donor_nmismatches, int donor_ncolordiffs,
		     double donor_prob, Genomicpos_T left, int querylength,
		     bool plusp, bool sensep, char *genomicseg, char *query, Chrnum_T chrnum,
		     Genomicpos_T chroffset, bool knownp, bool trim_ends_p, bool dibasep, bool cmetp) {
  Substring_T new;
  int querystart, queryend, extraleft, extraright;
  Genomicpos_T genomicstart, genomicend, alignstart, alignend;

  if (dibasep) {
  } else {
    donor_ncolordiffs = 0;
  }

  if (plusp == true) {
    genomicstart = left;
    genomicend = left + querylength;
    if (sensep == true) {
      querystart = 0;
      queryend = donor_pos;
      extraleft = 0;
      extraright = 2;
      alignstart = genomicstart;
      alignend = genomicstart + donor_pos;
    } else {
      extraleft = 2;
      extraright = 0;
      querystart = donor_pos;
      queryend = querylength;
      alignstart = genomicstart + donor_pos;
      alignend = genomicend;
    }

  } else {
    genomicstart = left + querylength;
    genomicend = left;
    if (sensep == true) {
      querystart = 0;
      queryend = donor_pos;
      extraleft = 0;
      extraright = 2;
      alignstart = genomicstart;
      alignend = genomicstart - donor_pos;
    } else {
      extraleft = 2;
      extraright = 0;
      querystart = donor_pos;
      queryend = querylength;
      alignstart = genomicstart - donor_pos;
      alignend = genomicend;
    }
  }

  new = Substring_new(donor_nmismatches,donor_ncolordiffs,chrnum,chroffset,
		      genomicstart,genomicend,querystart,queryend,querylength,
		      alignstart,alignend,/*genomiclength*/querylength,
		      extraleft,extraright,genomicseg,query,plusp,trim_ends_p,dibasep,cmetp);

  if (plusp == true) {
    new->chimera_modelpos = left + donor_pos;
  } else {
    new->chimera_modelpos = left + (querylength - donor_pos);
  }
  new->chimera_sensep = sensep;
  new->chimera_knownp = knownp;
  new->chimera_pos = donor_pos;
  new->chimera_prob = donor_prob;

  return new;
}


Substring_T
Substring_new_acceptor (int acceptor_pos, int acceptor_nmismatches, int acceptor_ncolordiffs,
			double acceptor_prob, Genomicpos_T left, int querylength,
			bool plusp, bool sensep, char *genomicseg, char *query, Chrnum_T chrnum,
			Genomicpos_T chroffset, bool knownp, bool trim_ends_p, bool dibasep, bool cmetp) {
  Substring_T new;
  int querystart, queryend, extraleft, extraright;
  Genomicpos_T genomicstart, genomicend, alignstart, alignend;

  if (dibasep) {
  } else {
    acceptor_ncolordiffs = 0;
  }

  if (plusp == true) {
    genomicstart = left;
    genomicend = left + querylength;
    if (sensep == true) {
      extraleft = 2;
      extraright = 0;
      querystart = acceptor_pos;
      queryend = querylength;
      alignstart = genomicstart + acceptor_pos;
      alignend = genomicend;
    } else {
      querystart = 0;
      queryend = acceptor_pos;
      extraleft = 0;
      extraright = 2;
      alignstart = genomicstart;
      alignend = genomicstart + acceptor_pos;
    }

  } else {
    genomicstart = left + querylength;
    genomicend = left;
    if (sensep == true) {
      extraleft = 2;
      extraright = 0;
      querystart = acceptor_pos;
      queryend = querylength;
      alignstart = genomicstart - acceptor_pos;
      alignend = genomicend;
    } else {
      querystart = 0;
      queryend = acceptor_pos;
      extraleft = 0;
      extraright = 2;
      alignstart = genomicstart;
      alignend = genomicstart - acceptor_pos;
    }
  }

  new = Substring_new(acceptor_nmismatches,acceptor_ncolordiffs,chrnum,chroffset,
		      genomicstart,genomicend,querystart,queryend,querylength,
		      alignstart,alignend,/*genomiclength*/querylength,
		      extraleft,extraright,genomicseg,query,plusp,trim_ends_p,dibasep,cmetp);

  if (plusp == true) {
    new->chimera_modelpos = left + acceptor_pos;
  } else {
    new->chimera_modelpos = left + (querylength - acceptor_pos);
  }
  new->chimera_sensep = sensep;
  new->chimera_knownp = knownp;
  new->chimera_pos = acceptor_pos;
  new->chimera_prob = acceptor_prob;

  return new;
}


void
Substring_assign_donor_prob (Substring_T donor, Genome_T genome, IIT_T chromosome_iit) {
  Chrnum_T chrnum_ignore;
  int nunknowns;
  char gbuffer1[MAX_QUERYLENGTH+1], gbuffer2[MAX_QUERYLENGTH+1];

  if (donor == NULL) {
    return;
  } else if (donor->chimera_knownp == false) {
    /* Prob already assigned */
    return;
  } else if (donor->plusp == donor->chimera_sensep) {
    Genome_fill_buffer(&chrnum_ignore,&nunknowns,genome,donor->chimera_modelpos-DONOR_MODEL_LEFT_MARGIN,
		       DONOR_MODEL_LEFT_MARGIN+DONOR_MODEL_RIGHT_MARGIN+1,gbuffer1,chromosome_iit);
    donor->chimera_prob = Maxent_donor_prob(gbuffer1);
    debug4(printf("\ndonor, plusp %d, sensep %d, at %u: %f\n",
		  donor->plusp,donor->chimera_sensep,donor->chimera_modelpos,donor->chimera_prob));
    debug4(printf("%s\n",gbuffer1));
    debug4(printf("%*s^^\n",DONOR_MODEL_LEFT_MARGIN,""));
  } else {
    Genome_fill_buffer(&chrnum_ignore,&nunknowns,genome,donor->chimera_modelpos-DONOR_MODEL_RIGHT_MARGIN-1,
		       DONOR_MODEL_LEFT_MARGIN+DONOR_MODEL_RIGHT_MARGIN+1,gbuffer1,chromosome_iit);
    make_complement_buffered(gbuffer2,gbuffer1,DONOR_MODEL_LEFT_MARGIN+DONOR_MODEL_RIGHT_MARGIN+1);
    donor->chimera_prob = Maxent_donor_prob(gbuffer2);
    debug4(printf("\ndonor, plusp %d, sensep %d, at %u: %f\n",
		  donor->plusp,donor->chimera_sensep,donor->chimera_modelpos,donor->chimera_prob));
    debug4(printf("%s\n",gbuffer2));
    debug4(printf("%*s^^\n",DONOR_MODEL_LEFT_MARGIN,""));
  }

  return;
}

void
Substring_assign_acceptor_prob (Substring_T acceptor, Genome_T genome, IIT_T chromosome_iit) {
  Chrnum_T chrnum_ignore;
  int nunknowns;
  char gbuffer1[MAX_QUERYLENGTH+1], gbuffer2[MAX_QUERYLENGTH+1];

  if (acceptor == NULL) {
    return;
  } else if (acceptor->chimera_knownp == false) {
    /* Prob already assigned */
    return;
  } else if (acceptor->plusp == acceptor->chimera_sensep) {
    Genome_fill_buffer(&chrnum_ignore,&nunknowns,genome,acceptor->chimera_modelpos-ACCEPTOR_MODEL_LEFT_MARGIN,
		       ACCEPTOR_MODEL_LEFT_MARGIN+ACCEPTOR_MODEL_RIGHT_MARGIN+1,gbuffer1,chromosome_iit);
    acceptor->chimera_prob = Maxent_acceptor_prob(gbuffer1);
    debug4(printf("\nacceptor, plusp %d, sensep %d, at %u: %f\n",
		  acceptor->plusp,acceptor->chimera_sensep,acceptor->chimera_modelpos,acceptor->chimera_prob));
    debug4(printf("%s\n",gbuffer1));
    debug4(printf("%*s^^\n",ACCEPTOR_MODEL_LEFT_MARGIN-2,""));
  } else {
    Genome_fill_buffer(&chrnum_ignore,&nunknowns,genome,acceptor->chimera_modelpos-ACCEPTOR_MODEL_RIGHT_MARGIN-1,
		       ACCEPTOR_MODEL_LEFT_MARGIN+ACCEPTOR_MODEL_RIGHT_MARGIN+1,gbuffer1,chromosome_iit);
    make_complement_buffered(gbuffer2,gbuffer1,ACCEPTOR_MODEL_LEFT_MARGIN+ACCEPTOR_MODEL_RIGHT_MARGIN+1);
    acceptor->chimera_prob = Maxent_acceptor_prob(gbuffer2);
    debug4(printf("\nacceptor, plusp %d, sensep %d, at %u: %f\n",
		  acceptor->plusp,acceptor->chimera_sensep,acceptor->chimera_modelpos,acceptor->chimera_prob));
    debug4(printf("%s\n",gbuffer2));
    debug4(printf("%*s^^\n",ACCEPTOR_MODEL_LEFT_MARGIN-2,""));
  }

  return;
}



static int
ascending_pos_cmp (const void *a, const void *b) {
  Substring_T x = * (Substring_T *) a;
  Substring_T y = * (Substring_T *) b;

  if (x->chimera_pos < y->chimera_pos) {
    return -1;
  } else if (x->chimera_pos > y->chimera_pos) {
    return +1;
  } else if (x->genomicstart < y->genomicstart) {
    return -1;
  } else if (x->genomicstart > y->genomicstart) {
    return +1;
  } else if (x->chimera_knownp == true && y->chimera_knownp == false) {
    return -1;
  } else if (y->chimera_knownp == true && x->chimera_knownp == false) {
    return +1;
  } else {
    return 0;
  }
}

static int
descending_pos_cmp (const void *a, const void *b) {
  Substring_T x = * (Substring_T *) a;
  Substring_T y = * (Substring_T *) b;

  if (x->chimera_pos < y->chimera_pos) {
    return -1;
  } else if (x->chimera_pos > y->chimera_pos) {
    return +1;
  } else if (x->genomicstart > y->genomicstart) {
    return -1;
  } else if (x->genomicstart < y->genomicstart) {
    return +1;
  } else if (x->chimera_knownp == true && y->chimera_knownp == false) {
    return -1;
  } else if (y->chimera_knownp == true && x->chimera_knownp == false) {
    return +1;
  } else {
    return 0;
  }
}

List_T
Substring_sort_chimera_halves (List_T hitlist, bool ascendingp) {
  List_T sorted = NULL, p;
  Substring_T x, *hits;
  int n, i, j;
  bool *eliminate;

  n = List_length(hitlist);
  debug(printf("Checking %d spliceends for duplicates...",n));
  if (n == 0) {
    debug(printf("\n"));
    return NULL;
  }

  hits = (Substring_T *) CALLOC(n,sizeof(Substring_T));
  for (p = hitlist, i = 0; p != NULL; p = p->rest) {
    hits[i++] = (Substring_T) p->first;
  }
  List_free(&hitlist);

  if (ascendingp == true) {
    qsort(hits,n,sizeof(Substring_T),ascending_pos_cmp);
  } else {
    qsort(hits,n,sizeof(Substring_T),descending_pos_cmp);
  }

  /* Check for duplicates */
  eliminate = (bool *) CALLOC(n,sizeof(bool));
  for (i = 0; i < n; i++) {
    x = hits[i];
    j = i+1;
    while (j < n && hits[j]->chimera_pos == x->chimera_pos && hits[j]->genomicstart == x->genomicstart) {
      eliminate[j] = true;
      j++;
    }
  }

  debug(j = 0);
  for (i = n-1; i >= 0; i--) {
    x = hits[i];
    if (eliminate[i] == false) {
      sorted = List_push(sorted,x);
    } else {
      Substring_free(&x);
      debug(j++);
    }
  }
  debug(printf("%d eliminated\n",j));

  FREE(hits);
  FREE(eliminate);

  return sorted;
}


static void
print_snp_labels (Substring_T this, IIT_T snps_iit, int *snps_divint_crosstable) {
  int *snps, nsnps, querypos, c, i;
  char *label, *alleles;
  bool allocp, printp = false;
  Interval_T interval;
  Genomicpos_T position;

  if (this->plusp == true) {
    snps = IIT_get_with_divno(&nsnps,snps_iit,
			      snps_divint_crosstable[this->chrnum],
			      this->alignstart - this->chroffset + 1,this->alignend - this->chroffset - 1,
			      /*sortp*/false);
  } else {
    snps = IIT_get_with_divno(&nsnps,snps_iit,
			      snps_divint_crosstable[this->chrnum],
			      this->alignend - this->chroffset + 1,this->alignstart - this->chroffset - 1,
			      /*sortp*/false);
  }

  printf(",snps:");

  if (this->plusp) {

#if 0
    for (i = 0; i < nsnps; i++) {
      interval = IIT_interval(snps_iit,snps[i]);
      position = Interval_low(interval);
      querypos = position - (this->genomicstart-this->chroffset) - 1;
      printf("%d ",querypos);
    }
    printf("\n");
#endif

    for (i = 0; i < nsnps; i++) {
      label = IIT_label(snps_iit,snps[i],&allocp);
    
      interval = IIT_interval(snps_iit,snps[i]);
      position = Interval_low(interval);
      querypos = position - (this->genomicstart-this->chroffset) - 1;
      if (querypos < 0 || querypos >= this->genomiclength) {
	abort();
      }
      c = this->genomicdir[querypos];

      /* printf("\n%d%c\n",querypos,c); */
      if (islower(c)) {
	c = toupper(c);
	alleles = IIT_typestring(snps_iit,Interval_type(interval));
	
	if (c == alleles[0] || c == alleles[1]) {
	  if (printp) {
	    printf("|");
	  }
	  printf("%d@",querypos+1);
	  printf("%s",label);
	  printp = true;
	}
	
      }
      if (allocp) FREE(label);
    }

  } else {
    for (i = nsnps-1; i >= 0; i--) {
      label = IIT_label(snps_iit,snps[i],&allocp);
    
      interval = IIT_interval(snps_iit,snps[i]);
      position = Interval_low(interval);
      querypos = (this->genomicstart-this->chroffset) - position;
      if (querypos < 0 || querypos >= this->genomiclength) {
	abort();
      }
      c = complCode[(int) this->genomicdir[querypos]];

      /* printf("\n%d%c\n",querypos,c); */
      if (islower(c)) {
	c = toupper(c);
	alleles = IIT_typestring(snps_iit,Interval_type(interval));
	
	if (c == alleles[0] || c == alleles[1]) {
	  if (printp) {
	    printf("|");
	  }
	  printf("%d@",querypos+1);
	  printf("%s",label);
	  printp = true;
	}
	
      }
      if (allocp) FREE(label);
    }
  }

  FREE(snps);

  return;
}


static void
print_splicesite_labels (Substring_T this, IIT_T splicesites_iit, int *splicesites_divint_crosstable,
			 int typeint) {
  Genomicpos_T splicesitepos;
  int *splicesites, nsplicesites, i;
  char *label;
  bool allocp;

  if (this->plusp == true) {
    splicesitepos = this->genomicstart - this->chroffset + this->chimera_pos;
  } else {
    splicesitepos = this->genomicstart - this->chroffset - this->chimera_pos;
  }
  if (this->chimera_knownp == true) {
    splicesites = IIT_get_exact_multiple_with_divno(&nsplicesites,splicesites_iit,
						    splicesites_divint_crosstable[this->chrnum],
						    splicesitepos,splicesitepos+1U,typeint);
  } else {
#if 0
    splicesites = IIT_get_with_divno(&nsplicesites,splicesites_iit,
				     splicesites_divint_crosstable[this->chrnum],
				     splicesitepos,splicesitepos+1U,typeint);
#else
    nsplicesites = 0;
#endif
  }

  if (nsplicesites > 0) {
    printf(",label:");
    label = IIT_label(splicesites_iit,splicesites[0],&allocp);
    printf("%s",label);
    if (allocp) FREE(label);

    for (i = 1; i < nsplicesites; i++) {
      label = IIT_label(splicesites_iit,splicesites[i],&allocp);
      printf("|%s",label);
      if (allocp) FREE(label);
    }
    FREE(splicesites);
  }

  return;
}


static void
print_forward (char *string, int n) {
  int i;
  
  for (i = 0; i < n; i++) {
    printf("%c",string[i]);
  }
  return;
}




static void
print_lc (char *string, int n) {
  int i;
  
  for (i = 0; i < n; i++) {
    printf("%c",(char) tolower(string[i]));
  }
  return;
}



static void
print_revcomp (char *nt, int len) {
  int i;

  for (i = len-1; i >= 0; --i) {
    printf("%c",complCode[(int) nt[i]]);
  }
  return;
}

static void
print_revcomp_lc (char *nt, int len) {
  int i;

  for (i = len-1; i >= 0; --i) {
    printf("%c",(char) tolower(complCode[(int) nt[i]]));
  }
  return;
}

static void
print_genomic (Substring_T substring, char *deletion, int deletionlength, bool invertp) {

  if (invertp == false) {
    print_forward(substring->genomicdir,substring->queryend);
    if (deletion != NULL) {
      print_lc(deletion,deletionlength);
    }
    print_forward(&(substring->genomicdir[substring->queryend]),substring->querylength - substring->queryend);
    printf("\t");
    printf("%d..%d",1 + substring->querystart,substring->queryend);
  } else {
    print_revcomp(&(substring->genomicdir[substring->querystart]),substring->querylength - substring->querystart);
    if (deletion != NULL) {
      print_revcomp_lc(deletion,deletionlength);
    }
    print_revcomp(substring->genomicdir,substring->querystart);
    printf("\t");
    printf("%d..%d",1 + substring->querylength - substring->queryend,substring->querylength - substring->querystart);
  }
  return;
}


static void
print_coordinates (Substring_T substring, char *chr, bool invertp) {

  if (substring->plusp == true) {
    if (invertp == false) {
      printf("+%s:%u..%u",chr,substring->alignstart_trim - substring->chroffset + 1U,
	     substring->alignend_trim - substring->chroffset);
    } else {
      printf("-%s:%u..%u",chr,substring->alignend_trim - substring->chroffset,
	     substring->alignstart_trim - substring->chroffset + 1U);
    }
  } else {
    if (invertp == false) {
      printf("-%s:%u..%u",chr,substring->alignstart_trim - substring->chroffset,
	     substring->alignend_trim - substring->chroffset + 1U);
    } else {
      printf("+%s:%u..%u",chr,substring->alignend_trim - substring->chroffset + 1U,
	     substring->alignstart_trim - substring->chroffset);
    }
  }

  return;
}


void
Substring_print_single (Substring_T substring, Hittype_T hittype, Sequence_T queryseq,
			char *chr, int querylength,
			bool invertp, IIT_T snps_iit, int *snps_divint_crosstable) {

  if (substring->genomicdir == NULL) {
    /* Exact match */
    if (invertp == false) {
      Sequence_print_oneline_uc(stdout,queryseq);
    } else {
      Sequence_print_oneline_revcomp_uc(stdout,queryseq);
    }
  } else {
    if (invertp == false) {
      printf("%s",substring->genomicdir);
    } else {
      print_revcomp(substring->genomicdir,substring->genomiclength);
    }
  }

  printf("\t");
  if (invertp == false) {
    printf("%d..%d",1 + substring->trim_left,querylength - substring->trim_right);
  } else {
    printf("%d..%d",1 + substring->trim_right,querylength - substring->trim_left);
  }

  printf("\t");
#if 0
  if (substring->plusp == true) {
    /* Add 1 to report in 1-based coordinates */
    chrpos = substring->genomicstart - substring->chroffset + 1U;
    if (invertp == false) {
      printf("+%s:%u..%u",chr,
	     chrpos + substring->trim_left,
	     chrpos + substring->genomiclength-1 - substring->trim_right);
    } else {
      printf("-%s:%u..%u",chr,
	     chrpos + substring->genomiclength-1 - substring->trim_right,
	     chrpos + substring->trim_left);
    }
  } else {
    /* Sometimes genomicstart goes into next chromosome, so subtract 1 */
    chrpos = substring->genomicstart - substring->chroffset;
    if (invertp == false) {
      printf("-%s:%u..%u",chr,
	     chrpos - substring->trim_left,
	     chrpos - substring->genomiclength+1 + substring->trim_right);
    } else {
      printf("+%s:%u..%u",chr,
	     chrpos - substring->genomiclength+1 + substring->trim_right,
	     chrpos - substring->trim_left);
    }
  }
#else
  print_coordinates(substring,chr,invertp);
#endif


  printf("\t");
  if (hittype == EXACT) {
    /* Previously wrote exact, but now exact is found by looking for start:0..end:0 */
    printf("start:0..end:0,sub:0");
  } else if (hittype == SUB) {
    if (substring->nmismatches == 0) {
    /* Previously wrote exact, but now exact is found by looking for start:0..end:0 */
      printf("start:0..end:0,sub:0");
    } else if (invertp == false) {
      printf("start:%d..end:%d,sub:%d",
	     substring->trim_left,substring->trim_right,substring->nmismatches);
    } else {
      printf("start:%d..end:%d,sub:%d",
	     substring->trim_right,substring->trim_left,substring->nmismatches);
    }
  }
  if (print_ncolordiffs_p) {
    printf("+%d=%d",substring->ncolordiffs,substring->nmismatches+substring->ncolordiffs);
  } else if (print_nsnpdiffs_p) {
    printf("+%d=%d",substring->nsnpdiffs,substring->nmismatches+substring->nsnpdiffs);
    if (print_snplabels_p && substring->nsnpdiffs > 0) {
      print_snp_labels(substring,snps_iit,snps_divint_crosstable);
    }
  }

  return;
}


void
Substring_print_insertion_1 (Substring_T substring1, Substring_T substring2, int nindels, char *chr,
			     bool invertp, IIT_T snps_iit, int *snps_divint_crosstable) {
  Substring_T substring;

#if 0
  /* First part */
  if (invertp == false) {
    substring = substring1;
    /* printf("%.*s",this->indel_pos,substring->genomic); */
    print_forward(substring->genomicdir,indel_pos);
    for (i = indel_pos; i < querylength; i++) {
      printf("-");
    }
    printf("\t");
    printf("%d..%d",1 + substring->trim_left,indel_pos);
  } else {
    substring = substring2;
    print_revcomp(&(substring->genomicdir[indel_pos+nindels]),querylength-indel_pos-nindels);
    for (i = querylength-indel_pos-nindels; i < querylength; i++) {
      printf("-");
    }
    printf("\t");
    printf("%d..%d",1 + substring->trim_right,querylength-indel_pos-nindels);
  }
#else
  if (invertp == false) {
    substring = substring1;
    print_genomic(substring1,/*deletion*/NULL,/*deletionlength*/0,/*invertp*/false);
  } else {
    substring = substring2;
    print_genomic(substring2,/*deletion*/NULL,/*deletionlength*/0,/*invertp*/true);
  }
#endif

  printf("\t");
#if 0
  if (substring->plusp == true) {
    /* Add 1 to report in 1-based coordinates */
    chrpos = substring->genomicstart - substring->chroffset + 1U;
    if (invertp == false) {
      printf("+%s:%u..%u",chr,
	     chrpos + substring->trim_left,
	     chrpos+indel_pos-1U);
    } else {
      printf("-%s:%u..%u",chr,
	     chrpos+querylength-nindels-1U - substring->trim_right,
	     chrpos+indel_pos);
    }
  } else {
    /* Sometimes genomicstart goes into next chromosome, so subtract 1 */
    chrpos = substring->genomicstart - substring->chroffset;
    if (invertp == false) {
      printf("-%s:%u..%u",chr,
	     chrpos - substring->trim_left,
	     chrpos-indel_pos+1U);
    } else {
      printf("+%s:%u..%u",chr,
	     chrpos-querylength+nindels+1U + substring->trim_right,
	     chrpos-indel_pos);
    }
  }
#else
  print_coordinates(substring,chr,invertp);
#endif


  printf("\t");
  if (invertp == false) {
    printf("start:%d..ins:%d,sub:%d",substring->trim_left,nindels,substring->nmismatches);
  } else {
    printf("start:%d..ins:%d,sub:%d",substring->trim_right,nindels,substring->nmismatches);
  }
  if (print_ncolordiffs_p) {
    printf("+%d=%d",substring->ncolordiffs,substring->nmismatches+substring->ncolordiffs);
  } else if (print_nsnpdiffs_p) {
    printf("+%d=%d",substring->nsnpdiffs,substring->nmismatches+substring->nsnpdiffs);
    if (print_snplabels_p && substring->nsnpdiffs > 0) {
      print_snp_labels(substring,snps_iit,snps_divint_crosstable);
    }
  }


  return;
}

void
Substring_print_insertion_2 (Substring_T substring1, Substring_T substring2, int nindels, char *chr,
			     bool invertp, IIT_T snps_iit, int *snps_divint_crosstable) {
  Substring_T substring;

#if 0
  /* Second part */
  if (invertp == false) {
    substring = substring2;
    for (i = 0; i < indel_pos+nindels; i++) {
      printf("-");
    }
    printf("%s",&(substring->genomicdir[indel_pos+nindels]));
    printf("\t");
    printf("%d..%d",1+indel_pos+nindels,querylength - substring->trim_right);
  } else {
    substring = substring1;
    for (i = 0; i < querylength-indel_pos; i++) {
      printf("-");
    }
    print_revcomp(substring->genomicdir,indel_pos);
    printf("\t");
    printf("%d..%d",1+querylength-indel_pos,querylength - substring->trim_left);
  }
#else
  if (invertp == false) {
    substring = substring2;
    print_genomic(substring2,/*deletion*/NULL,/*deletionlength*/0,/*invertp*/false);
  } else {
    substring = substring1;
    print_genomic(substring1,/*deletion*/NULL,/*deletionlength*/0,/*invertp*/true);
  }
#endif


  printf("\t");
#if 0
  if (substring->plusp == true) {
    if (invertp == false) {
      chrpos = substring->genomicstart - substring->chroffset + 1U;
      printf("+%s:%u..%u",chr,
	     chrpos+indel_pos,
	     chrpos+querylength-nindels-1U - substring->trim_right);
    } else {
      printf("-%s:%u..%u",chr,
	     chrpos+indel_pos-1U,
	     chrpos + substring->trim_left);
    }
  } else {
    chrpos = substring->genomicstart - substring->chroffset + 1U;
    if (invertp == false) {
      printf("-%s:%u..%u",chr,
	     chrpos-indel_pos,
	     chrpos-querylength+nindels+1U + substring->trim_right);
    } else {
      printf("+%s:%u..%u",chr,
	     chrpos-indel_pos+1U,
	     chrpos - substring->trim_left);
    }
  }
#else
  print_coordinates(substring,chr,invertp);
#endif


  printf("\t");
  if (invertp == false) {
    printf("ins:%d..end:%d,sub:%d",nindels,substring->trim_right,substring->nmismatches);
  } else {
    printf("ins:%d..end:%d,sub:%d",nindels,substring->trim_left,substring->nmismatches);
  }
  if (print_ncolordiffs_p) {
    printf("+%d=%d",substring->ncolordiffs,substring->nmismatches+substring->ncolordiffs);
  } else if (print_nsnpdiffs_p) {
    printf("+%d=%d",substring->nsnpdiffs,substring->nmismatches+substring->nsnpdiffs);
    if (print_snplabels_p && substring->nsnpdiffs > 0) {
      print_snp_labels(substring,snps_iit,snps_divint_crosstable);
    }
  }

  return;
}


void
Substring_print_deletion_1 (Substring_T substring1, Substring_T substring2, int nindels, 
			    char *deletion, char *chr, 
			    bool invertp, IIT_T snps_iit, int *snps_divint_crosstable) {
  Substring_T substring;

#if 0
  /* First part */
  if (invertp == false) {
    substring = substring1;
    /* printf("%.*s",indel_pos,substring->genomicdir); */
    print_forward(substring->genomicdir,indel_pos);
    print_lc(deletion,nindels);
    for (i = indel_pos+nindels; i < querylength; i++) {
      printf("-");
    }
    printf("\t");
    printf("%d..%d",1 + substring->trim_left,indel_pos);
  } else {
    substring = substring2;
    print_revcomp(&(substring->genomicdir[indel_pos]),querylength-indel_pos);
    print_revcomp_lc(deletion,nindels);
    for (i = querylength-indel_pos+nindels; i < querylength; i++) {
      printf("-");
    }
    printf("\t");
    printf("%d..%d",1 + substring->trim_right,querylength-indel_pos);
  }
#else
  if (invertp == false) {
    substring = substring1;
    print_genomic(substring1,deletion,nindels,/*invertp*/false);
  } else {
    substring = substring2;
    print_genomic(substring2,deletion,nindels,/*invertp*/true);
  }
#endif

  printf("\t");
#if 0
  if (substring->plusp == true) {
    /* Add 1 to report in 1-based coordinates */
    chrpos = substring->genomicstart - substring->chroffset + 1U;
    if (invertp == false) {
      printf("+%s:%u..%u",chr,
	     chrpos + substring->trim_left,
	     chrpos+indel_pos-1U);
    } else {
      printf("-%s:%u..%u",chr,
	     chrpos+querylength+nindels-1U - substring->trim_right,
	     chrpos+indel_pos+nindels);
    }
  } else {
    /* Sometimes genomicstart goes into next chromosome, so subtract 1 */
    chrpos = substring->genomicstart - substring->chroffset;
    if (invertp == false) {
      printf("-%s:%u..%u",chr,
	     chrpos - substring->trim_left,
	     chrpos-indel_pos+1U);
    } else {
      printf("+%s:%u..%u",chr,
	     chrpos-querylength-nindels+1U + substring->trim_right,
	     chrpos-indel_pos-nindels);
    }
  }
#else
  print_coordinates(substring,chr,invertp);
#endif


  printf("\t");
  if (invertp == false) {
    printf("start:%d..del:%d,sub:%d",substring->trim_left,nindels,substring->nmismatches);
  } else {
    printf("start:%d..del:%d,sub:%d",substring->trim_right,nindels,substring->nmismatches);
  }
  if (print_ncolordiffs_p) {
    printf("+%d=%d",substring->ncolordiffs,substring->nmismatches+substring->ncolordiffs);
  } else if (print_nsnpdiffs_p) {
    printf("+%d=%d",substring->nsnpdiffs,substring->nmismatches+substring->nsnpdiffs);
    if (print_snplabels_p && substring->nsnpdiffs > 0) {
      print_snp_labels(substring,snps_iit,snps_divint_crosstable);
    }
  }

  return;
}


void
Substring_print_deletion_2 (Substring_T substring1, Substring_T substring2, int nindels, 
			    char *deletion, char *chr, 
			    bool invertp, IIT_T snps_iit, int *snps_divint_crosstable) {
  Substring_T substring;

#if 0
  /* Second part */
  if (invertp == false) {
    substring = substring2;
    for (i = 0; i < indel_pos; i++) {
      printf("-");
    }
    printf("%s",&(substring->genomicdir[indel_pos]));
    printf("\t");
    printf("%d..%d",1+indel_pos,querylength - substring->trim_right);
  } else {
    substring = substring1;
    for (i = 0; i < querylength-indel_pos; i++) {
      printf("-");
    }
    print_revcomp(substring->genomicdir,indel_pos);
    printf("\t");
    printf("%d..%d",1+querylength-indel_pos,querylength - substring->trim_left);
  }
#else
  if (invertp == false) {
    substring = substring2;
    print_genomic(substring2,/*deletion*/NULL,/*deletionlength*/0,/*invertp*/false);
  } else {
    substring = substring1;
    print_genomic(substring1,/*deletion*/NULL,/*deletionlength*/0,/*invertp*/true);
  }
#endif


  printf("\t");
#if 0
  if (substring->plusp == true) {
    chrpos = substring->genomicstart - substring->chroffset + 1U;
    if (invertp == false) {
      printf("+%s:%u..%u",chr,
	     chrpos+nindels+indel_pos,
	     chrpos+nindels+querylength-1U - substring->trim_right);
    } else {
      printf("-%s:%u..%u",chr,
	     chrpos+indel_pos-1U,
	     chrpos + substring->trim_left);
    }
  } else {
    chrpos = substring->genomicstart - substring->chroffset;
    if (invertp == false) {
      printf("-%s:%u..%u",chr,
	     chrpos-nindels-indel_pos,
	     chrpos-nindels-querylength+1U + substring->trim_right);
    } else {
      printf("+%s:%u..%u",chr,
	     chrpos-indel_pos+1U,
	     chrpos - substring->trim_left);
    }
  }
#else
  print_coordinates(substring,chr,invertp);
#endif

  printf("\t");
  if (invertp == false) {
    printf("del:%d..end:%d,sub:%d",nindels,substring->trim_right,substring->nmismatches);
  } else {
    printf("del:%d..end:%d,sub:%d",nindels,substring->trim_left,substring->nmismatches);
  }
  if (print_ncolordiffs_p) {
    printf("+%d=%d",substring->ncolordiffs,substring->nmismatches+substring->ncolordiffs);
  } else if (print_nsnpdiffs_p) {
    printf("+%d=%d",substring->nsnpdiffs,substring->nmismatches+substring->nsnpdiffs);
    if (print_snplabels_p && substring->nsnpdiffs > 0) {
      print_snp_labels(substring,snps_iit,snps_divint_crosstable);
    }
  }
  
  return;
}


static void
print_splice_distance (Substring_T donor, Substring_T acceptor, Genomicpos_T distance, bool sensep) {

  if (donor == NULL || acceptor == NULL) {
    /* Don't print anything */
  } else if (distance == 0U) {
    printf(",%s",TRANSLOCATION_TEXT);
  } else {
    printf(",splice_dist:%u",distance);
    if (donor->plusp != acceptor->plusp) {
      printf(",%s",INVERSION_TEXT);
    } else if (donor->plusp == true) {
      if (sensep == true) {
	if (acceptor->genomicstart < donor->genomicstart) {
	  printf(",%s",SCRAMBLE_TEXT);
	}
      } else {
	if (donor->genomicstart < acceptor->genomicstart) {
	  printf(",%s",SCRAMBLE_TEXT);
	}
      }
    } else {
      if (sensep == true) {
	if (donor->genomicstart < acceptor->genomicstart) {
	  printf(",%s",SCRAMBLE_TEXT);
	}
      } else {
	if (acceptor->genomicstart < donor->genomicstart) {
	  printf(",%s",SCRAMBLE_TEXT);
	}
      }
    }
  }

  return;
}



void
Substring_print_donor (Substring_T donor, bool sensep, bool invertp,
		       IIT_T chromosome_iit, IIT_T snps_iit, int *snps_divint_crosstable, IIT_T splicesites_iit,
		       int donor_typeint, int *splicesites_divint_crosstable, bool endp,
		       Substring_T acceptor, Genomicpos_T chimera_distance) {
  char *chr;
  bool allocp;

#if 0
  /* printf("donor chimera_pos is %d\n",donor->chimera_pos); */
  if (sensep == true && invertp == false) {
    /* printf("%.*s",donor->chimera_pos,donor->genomic); */
    print_forward(donor->genomicdir,donor->chimera_pos);
    print_lc(&(donor->genomicdir[donor->chimera_pos]),querylength - donor->chimera_pos);
    printf("\t");
    printf("%d..%d",1 + donor->trim_left,donor->chimera_pos);
  } else if (sensep == true && invertp == true) {
    print_revcomp_lc(&(donor->genomicdir[donor->chimera_pos]),querylength - donor->chimera_pos);
    print_revcomp(donor->genomicdir,donor->chimera_pos);
    printf("\t");
    printf("%d..%d",querylength-donor->chimera_pos+1,querylength - donor->trim_left);
  } else if (sensep == false && invertp == false) {
    print_lc(donor->genomicdir,donor->chimera_pos);
    print_forward(&(donor->genomicdir[donor->chimera_pos]),querylength-donor->chimera_pos);
    printf("\t");
    printf("%d..%d",donor->chimera_pos+1,querylength - donor->trim_right);
  } else if (sensep == false && invertp == true) {
    print_revcomp(&(donor->genomicdir[donor->chimera_pos]),querylength-donor->chimera_pos);
    print_revcomp_lc(donor->genomicdir,donor->chimera_pos);
    printf("\t");
    printf("%d..%d",1 + donor->trim_right,querylength-donor->chimera_pos);
  }
#else
  print_genomic(donor,/*deletion*/NULL,/*deletionlength*/0,invertp);
#endif

  printf("\t");
#if 0
  if (sensep == true) {
    if (donor->plusp == true) {
      chr = IIT_label(chromosome_iit,donor->chrnum,&allocp);
      /* Add 1 to report in 1-based coordinates */
      chrpos = donor->genomicstart - donor->chroffset + 1U;
      if (invertp == false) {
	printf("+%s:%u..%u",chr,
	       chrpos + donor->trim_left,
	       chrpos+donor->chimera_pos-1);
      } else {
	printf("-%s:%u..%u",chr,
	       chrpos+donor->chimera_pos-1,
	       chrpos + donor->trim_left);
      }
    } else {
      /* Sometimes genomicstart goes into next chromosome, so subtract 1 */
      chr = IIT_label(chromosome_iit,donor->chrnum,&allocp);
      chrpos = donor->genomicstart - donor->chroffset;
      if (invertp == false) {
	printf("-%s:%u..%u",chr,
	       chrpos - donor->trim_left,
	       chrpos-donor->chimera_pos+1);
      } else {
	printf("+%s:%u..%u",chr,
	       chrpos-donor->chimera_pos+1,
	       chrpos - donor->trim_left);
      }
    }

  } else if (sensep == false) {
    if (donor->plusp == true) {
      chr = IIT_label(chromosome_iit,donor->chrnum,&allocp);
      /* Add 1 to report in 1-based coordinates */
      chrpos = donor->genomicstart - donor->chroffset + 1U;
      if (invertp == false) {
	printf("+%s:%u..%u",chr,
	       chrpos+donor->chimera_pos,
	       chrpos+querylength-1 - donor->trim_right);
      } else {
	printf("-%s:%u..%u",chr,
	       chrpos+querylength-1 - donor->trim_right,
	       chrpos+donor->chimera_pos);
      }
    } else {
      /* Sometimes genomicstart goes into next chromosome, so subtract 1 */
      chr = IIT_label(chromosome_iit,donor->chrnum,&allocp);
      chrpos = donor->genomicstart - donor->chroffset;
      if (invertp == false) {
	printf("-%s:%u..%u",chr,
	       chrpos-donor->chimera_pos,
	       chrpos-querylength+1 + donor->trim_right);
      } else {
	printf("+%s:%u..%u",chr,
	       chrpos-querylength+1 + donor->trim_right,
	       chrpos-donor->chimera_pos);
      }
    }
  }
#else
  chr = IIT_label(chromosome_iit,donor->chrnum,&allocp);
  print_coordinates(donor,chr,invertp);
#endif

  /* printf("donor chimera_pos is %d\n",donor->chimera_pos); */
  if (sensep == true && invertp == false) {
    printf("\t");
    printf("start:%d..%s:%.2f,dir:sense",
	   donor->trim_left,endp ? "end_donor" : "donor",donor->chimera_prob);
    print_splice_distance(donor,acceptor,chimera_distance,sensep);
  } else if (sensep == true && invertp == true) {
    printf("\t");
    printf("%s:%.2f..end:%d,dir:antisense",
	   endp ? "start_donor" : "donor",donor->chimera_prob,donor->trim_left);
    print_splice_distance(donor,acceptor,chimera_distance,sensep);
  } else if (sensep == false && invertp == false) {
    printf("\t");
    printf("%s:%.2f..end:%d,dir:antisense",
	   endp ? "start_donor" : "donor",donor->chimera_prob,donor->trim_right);
    print_splice_distance(donor,acceptor,chimera_distance,sensep);
  } else if (sensep == false && invertp == true) {
    printf("\t");
    printf("start:%d..%s:%.2f,dir:sense",
	   donor->trim_right,endp ? "end_donor" : "donor",donor->chimera_prob);
    print_splice_distance(donor,acceptor,chimera_distance,sensep);
  }

  printf(",sub:%d",donor->nmismatches);

  if (print_ncolordiffs_p) {
    printf("+%d=%d",donor->ncolordiffs,donor->nmismatches+donor->ncolordiffs);
  } else if (print_nsnpdiffs_p) {
    printf("+%d=%d",donor->nsnpdiffs,donor->nmismatches+donor->nsnpdiffs);
    if (print_snplabels_p && donor->nsnpdiffs > 0) {
      print_snp_labels(donor,snps_iit,snps_divint_crosstable);
    }
  }

  if (donor->chimera_knownp && splicesites_iit) {
    print_splicesite_labels(donor,splicesites_iit,splicesites_divint_crosstable,donor_typeint);
  }

  if (allocp == true) {
    FREE(chr);
  }

  return;
}

void 
Substring_print_acceptor (Substring_T acceptor, bool sensep, bool invertp,
			  IIT_T chromosome_iit, IIT_T snps_iit, int *snps_divint_crosstable, IIT_T splicesites_iit,
			  int acceptor_typeint, int *splicesites_divint_crosstable, bool endp,
			  Substring_T donor, Genomicpos_T chimera_distance) {
  char *chr;
  bool allocp;

#if 0
  /* printf("acceptor chimera_pos is %d\n",acceptor->chimera_pos); */
  if (sensep == true && invertp == false) {
    print_lc(acceptor->genomicdir,acceptor->chimera_pos);
    print_forward(&(acceptor->genomicdir[acceptor->chimera_pos]),querylength-acceptor->chimera_pos);
    printf("\t");
    printf("%d..%d",acceptor->chimera_pos+1,querylength - acceptor->trim_right);
  } else if (sensep == true && invertp == true) {
    print_revcomp(&(acceptor->genomicdir[acceptor->chimera_pos]),querylength-acceptor->chimera_pos);
    print_revcomp_lc(acceptor->genomicdir,acceptor->chimera_pos);
    printf("\t");
    printf("%d..%d",1 + acceptor->trim_right,querylength-acceptor->chimera_pos);
  } else if (sensep == false && invertp == false) {
    /* printf("%.*s",acceptor->chimera_pos,acceptor->genomicdir); */
    print_forward(acceptor->genomicdir,acceptor->chimera_pos);
    print_lc(&(acceptor->genomicdir[acceptor->chimera_pos]),querylength - acceptor->chimera_pos);
    printf("\t");
    printf("%d..%d",1 + acceptor->trim_left,acceptor->chimera_pos);
  } else if (sensep == false && invertp == true) {
    print_revcomp_lc(&(acceptor->genomicdir[acceptor->chimera_pos]),querylength - acceptor->chimera_pos);
    print_revcomp(acceptor->genomicdir,acceptor->chimera_pos);
    printf("\t");
    printf("%d..%d",querylength-acceptor->chimera_pos+1,querylength - acceptor->trim_left);
  }
#else
  print_genomic(acceptor,/*deletion*/NULL,/*deletionlength*/0,invertp);
#endif

  printf("\t");
#if 0
  if (sensep == true) {
    if (acceptor->plusp == true) {
      chr = IIT_label(chromosome_iit,acceptor->chrnum,&allocp);
      /* Add 1 to report in 1-based coordinates */
      chrpos = acceptor->genomicstart - acceptor->chroffset + 1U;
      if (invertp == false) {
	printf("+%s:%u..%u",chr,
	       chrpos+acceptor->chimera_pos,
	       chrpos+querylength-1 - acceptor->trim_right);
      } else {
	printf("-%s:%u..%u",chr,
	       chrpos+querylength-1 - acceptor->trim_right,
	       chrpos+acceptor->chimera_pos);
      }
    } else {
      /* Sometimes genomicstart goes into next chromosome, so subtract 1 */
      chr = IIT_label(chromosome_iit,acceptor->chrnum,&allocp);
      chrpos = acceptor->genomicstart - acceptor->chroffset;
      if (invertp == false) {
	printf("-%s:%u..%u",chr,
	       chrpos-acceptor->chimera_pos,
	       chrpos-querylength+1 + acceptor->trim_right);
      } else {
	printf("+%s:%u..%u",chr,
	       chrpos-querylength+1 + acceptor->trim_right,
	       chrpos-acceptor->chimera_pos);
      }
    }

  } else if (sensep == false) {
    if (acceptor->plusp == true) {
      chr = IIT_label(chromosome_iit,acceptor->chrnum,&allocp);
      /* Add 1 to report in 1-based coordinates */
      chrpos = acceptor->genomicstart - acceptor->chroffset + 1U;
      if (invertp == false) {
	printf("+%s:%u..%u",chr,
	       chrpos + acceptor->trim_left,
	       chrpos+acceptor->chimera_pos-1);
      } else {
	printf("-%s:%u..%u",chr,
	       chrpos+acceptor->chimera_pos-1,
	       chrpos + acceptor->trim_left);
      }
    } else {
      /* Sometimes genomicstart goes into next chromosome, so subtract 1 */
      chr = IIT_label(chromosome_iit,acceptor->chrnum,&allocp);
      chrpos = acceptor->genomicstart - acceptor->chroffset;
      if (invertp == false) {
	printf("-%s:%u..%u",chr,
	       chrpos - acceptor->trim_left,
	       chrpos-acceptor->chimera_pos+1);
      } else {
	printf("+%s:%u..%u",chr,
	       chrpos-acceptor->chimera_pos+1,
	       chrpos - acceptor->trim_left);
      }
    }
  }
#else
  chr = IIT_label(chromosome_iit,acceptor->chrnum,&allocp);
  print_coordinates(acceptor,chr,invertp);
#endif

  /* printf("acceptor chimera_pos is %d\n",acceptor->chimera_pos); */
  if (sensep == true && invertp == false) {
    printf("\t");
    printf("%s:%.2f..end:%d,dir:sense",
	   endp ? "start_acceptor" : "acceptor",acceptor->chimera_prob,acceptor->trim_right);
    print_splice_distance(donor,acceptor,chimera_distance,sensep);
  } else if (sensep == true && invertp == true) {
    printf("\t");
    printf("start:%d..%s:%.2f,dir:antisense",
	   acceptor->trim_right,endp ? "end_acceptor" : "acceptor",acceptor->chimera_prob);
    print_splice_distance(donor,acceptor,chimera_distance,sensep);
  } else if (sensep == false && invertp == false) {
    printf("\t");
    printf("start:%d..%s:%.2f,dir:antisense",
	   acceptor->trim_left,endp ? "end_acceptor" : "acceptor",acceptor->chimera_prob);
    print_splice_distance(donor,acceptor,chimera_distance,sensep);
  } else if (sensep == false && invertp == true) {
    printf("\t");
    printf("%s:%.2f..end:%d,dir:sense",
	   endp ? "start_acceptor" : "acceptor",acceptor->chimera_prob,acceptor->trim_left);
    print_splice_distance(donor,acceptor,chimera_distance,sensep);
  }

  printf(",sub:%d",acceptor->nmismatches);

  if (print_ncolordiffs_p) {
    printf("+%d=%d",acceptor->ncolordiffs,acceptor->nmismatches+acceptor->ncolordiffs);
  } else if (print_nsnpdiffs_p) {
    printf("+%d=%d",acceptor->nsnpdiffs,acceptor->nmismatches+acceptor->nsnpdiffs);
    if (print_snplabels_p && acceptor->nsnpdiffs > 0) {
      print_snp_labels(acceptor,snps_iit,snps_divint_crosstable);
    }
  }

  if (acceptor->chimera_knownp && splicesites_iit) {
    print_splicesite_labels(acceptor,splicesites_iit,splicesites_divint_crosstable,acceptor_typeint);
  }


  if (allocp == true) {
    FREE(chr);
  }

  return;
}




int
Substring_geneprob_eval_single (Substring_T substring, IIT_T geneprob_iit, IIT_T chromosome_iit) {
  int geneprob = 0, value;
  char *acc;
  int *matches;
  int nmatches, i;

  Genomicpos_T chrpos, pos;
  char *chr;
  bool allocp, alloc2p;

  if (substring->plusp == true) {
    chr = IIT_label(chromosome_iit,substring->chrnum,&allocp);
    /* Add 1 to report in 1-based coordinates */
    chrpos = substring->genomicstart - substring->chroffset + 1U;
    for (pos = chrpos + substring->trim_left; pos <= chrpos + substring->genomiclength-1 - substring->trim_right; pos++) {
      matches = IIT_get(&nmatches,geneprob_iit,/*divstring*/chr,/*coordstart*/pos,/*coordend*/pos,/*sortp*/true);
      if (nmatches > 1) {
	fprintf(stderr,"Multiple intervals for %s:%u\n",chr,pos);
	exit(9);
      } else {
	for (i = 0; i < nmatches; i++) {
	  /* interval = IIT_interval(geneprob_iit,matches[i]); */
	  acc = IIT_label(geneprob_iit,matches[i],&alloc2p);
	  value = atoi(acc);
	  geneprob += value;
	  if (alloc2p == true) {
	    FREE(acc);
	  }
	}
      }
      FREE(matches);
    }

  } else {
    /* Sometimes genomicstart goes into next chromosome, so subtract 1 */
    chr = IIT_label(chromosome_iit,substring->chrnum,&allocp);
    chrpos = substring->genomicstart - substring->chroffset;
    for (pos = chrpos - substring->trim_left; pos >= chrpos - substring->genomiclength+1 + substring->trim_right; pos--) {
      matches = IIT_get(&nmatches,geneprob_iit,/*divstring*/chr,/*coordstart*/pos,/*coordend*/pos,/*sortp*/true);
      if (nmatches > 1) {
	fprintf(stderr,"Multiple intervals for %s:%u\n",chr,pos);
	exit(9);
      } else {
	for (i = 0; i < nmatches; i++) {
	  /* interval = IIT_interval(geneprob_iit,matches[i]); */
	  acc = IIT_label(geneprob_iit,matches[i],&alloc2p);
	  value = atoi(acc);
	  geneprob += value;
	  if (alloc2p == true) {
	    FREE(acc);
	  }
	}
      }
      FREE(matches);
    }
  }

  if (allocp == true) {
    FREE(chr);
  }

  return geneprob;
}




int
Substring_geneprob_eval_insertion (Substring_T substring1, Substring_T substring2, int indel_pos, int nindels,
				   IIT_T geneprob_iit, IIT_T chromosome_iit, Sequence_T queryseq) {
  int geneprob = 0, value;
  char *acc;
  int *matches;
  int nmatches, i;

  Genomicpos_T chrpos, pos;
  Substring_T substring;
  char *chr;
  int querylength;
  bool allocp, alloc2p;

  querylength = Sequence_fulllength(queryseq);

  /* First part */
  substring = substring1;
  if (substring->plusp == true) {
    chr = IIT_label(chromosome_iit,substring->chrnum,&allocp);
    /* Add 1 to report in 1-based coordinates */
    chrpos = substring->genomicstart - substring->chroffset + 1U;
    for (pos = chrpos; pos <= chrpos + indel_pos - 1U; pos++) {
      matches = IIT_get(&nmatches,geneprob_iit,/*divstring*/chr,/*coordstart*/pos,/*coordend*/pos,/*sortp*/true);
      if (nmatches > 1) {
	fprintf(stderr,"Multiple intervals for %s:%u\n",chr,pos);
	exit(9);
      } else {
	for (i = 0; i < nmatches; i++) {
	  /* interval = IIT_interval(geneprob_iit,matches[i]); */
	  acc = IIT_label(geneprob_iit,matches[i],&alloc2p);
	  value = atoi(acc);
	  geneprob += value;
	  if (alloc2p == true) {
	    FREE(acc);
	  }
	}
      }
      FREE(matches);
    }

  } else {
    /* Sometimes genomicstart goes into next chromosome, so subtract 1 */
    chr = IIT_label(chromosome_iit,substring->chrnum,&allocp);
    chrpos = substring->genomicstart - substring->chroffset;
    for (pos = chrpos; pos >= chrpos - indel_pos + 1U; pos--) {
      matches = IIT_get(&nmatches,geneprob_iit,/*divstring*/chr,/*coordstart*/pos,/*coordend*/pos,/*sortp*/true);
      if (nmatches > 1) {
	fprintf(stderr,"Multiple intervals for %s:%u\n",chr,pos);
	exit(9);
      } else {
	for (i = 0; i < nmatches; i++) {
	  /* interval = IIT_interval(geneprob_iit,matches[i]); */
	  acc = IIT_label(geneprob_iit,matches[i],&alloc2p);
	  value = atoi(acc);
	  geneprob += value;
	  if (alloc2p == true) {
	    FREE(acc);
	  }
	}
      }
      FREE(matches);
    }
  }

  substring = substring2;
  if (substring->plusp == true) {
    for (pos = chrpos + indel_pos; pos <= chrpos + querylength - nindels - 1U; pos++) {
      matches = IIT_get(&nmatches,geneprob_iit,/*divstring*/chr,/*coordstart*/pos,/*coordend*/pos,/*sortp*/true);
      if (nmatches > 1) {
	fprintf(stderr,"Multiple intervals for %s:%u\n",chr,pos);
	exit(9);
      } else {
	for (i = 0; i < nmatches; i++) {
	  /* interval = IIT_interval(geneprob_iit,matches[i]); */
	  acc = IIT_label(geneprob_iit,matches[i],&alloc2p);
	  value = atoi(acc);
	  geneprob += value;
	  if (alloc2p == true) {
	    FREE(acc);
	  }
	}
      }
      FREE(matches);
    }
  } else {
    for (pos = chrpos - indel_pos; pos >= chrpos - querylength + nindels + 1U; pos--) {
      matches = IIT_get(&nmatches,geneprob_iit,/*divstring*/chr,/*coordstart*/pos,/*coordend*/pos,/*sortp*/true);
      if (nmatches > 1) {
	fprintf(stderr,"Multiple intervals for %s:%u\n",chr,pos);
	exit(9);
      } else {
	for (i = 0; i < nmatches; i++) {
	  /* interval = IIT_interval(geneprob_iit,matches[i]); */
	  acc = IIT_label(geneprob_iit,matches[i],&alloc2p);
	  value = atoi(acc);
	  geneprob += value;
	  if (alloc2p == true) {
	    FREE(acc);
	  }
	}
      }
      FREE(matches);
    }
  }

  if (allocp == true) {
    FREE(chr);
  }

  return geneprob;
}


int
Substring_geneprob_eval_deletion (Substring_T substring1, Substring_T substring2, int indel_pos, int nindels,
				  IIT_T geneprob_iit, IIT_T chromosome_iit, Sequence_T queryseq) {
  int geneprob = 0, value;
  char *acc;
  int *matches;
  int nmatches, i;

  Genomicpos_T chrpos, pos;
  Substring_T substring;
  char *chr;
  bool allocp, alloc2p;
  int querylength;

  querylength = Sequence_fulllength(queryseq);

  /* First part */
  substring = substring1;

  if (substring->plusp == true) {
    chr = IIT_label(chromosome_iit,substring->chrnum,&allocp);
    /* Add 1 to report in 1-based coordinates */
    chrpos = substring->genomicstart - substring->chroffset + 1U;
    for (pos = chrpos; pos <= chrpos + indel_pos - 1U; pos++) {
      matches = IIT_get(&nmatches,geneprob_iit,/*divstring*/chr,/*coordstart*/pos,/*coordend*/pos,/*sortp*/true);
      if (nmatches > 1) {
	fprintf(stderr,"Multiple intervals for %s:%u\n",chr,pos);
	exit(9);
      } else {
	for (i = 0; i < nmatches; i++) {
	  /* interval = IIT_interval(geneprob_iit,matches[i]); */
	  acc = IIT_label(geneprob_iit,matches[i],&alloc2p);
	  value = atoi(acc);
	  geneprob += value;
	  if (alloc2p == true) {
	    FREE(acc);
	  }
	}
      }
      FREE(matches);
    }
  } else {
    /* Sometimes genomicstart goes into next chromosome, so subtract 1 */
    chr = IIT_label(chromosome_iit,substring->chrnum,&allocp);
    chrpos = substring->genomicstart - substring->chroffset;
    for (pos = chrpos; pos >= chrpos - indel_pos + 1U; pos--) {
      matches = IIT_get(&nmatches,geneprob_iit,/*divstring*/chr,/*coordstart*/pos,/*coordend*/pos,/*sortp*/true);
      if (nmatches > 1) {
	fprintf(stderr,"Multiple intervals for %s:%u\n",chr,pos);
	exit(9);
      } else {
	for (i = 0; i < nmatches; i++) {
	  /* interval = IIT_interval(geneprob_iit,matches[i]); */
	  acc = IIT_label(geneprob_iit,matches[i],&alloc2p);
	  value = atoi(acc);
	  geneprob += value;
	  if (alloc2p == true) {
	    FREE(acc);
	  }
	}
      }
      FREE(matches);
    }
  }

  /* Second part */
  substring = substring2;

  if (substring->plusp == true) {
    for (pos = chrpos + nindels + indel_pos; pos <= chrpos + nindels + querylength - 1U; pos++) {
      matches = IIT_get(&nmatches,geneprob_iit,/*divstring*/chr,/*coordstart*/pos,/*coordend*/pos,/*sortp*/true);
      if (nmatches > 1) {
	fprintf(stderr,"Multiple intervals for %s:%u\n",chr,pos);
	exit(9);
      } else {
	for (i = 0; i < nmatches; i++) {
	  /* interval = IIT_interval(geneprob_iit,matches[i]); */
	  acc = IIT_label(geneprob_iit,matches[i],&alloc2p);
	  value = atoi(acc);
	  geneprob += value;
	  if (alloc2p == true) {
	    FREE(acc);
	  }
	}
      }
      FREE(matches);
    }
  } else {
    for (pos = chrpos - nindels - indel_pos; pos >= chrpos - nindels - querylength + 1U; pos--) {
      matches = IIT_get(&nmatches,geneprob_iit,/*divstring*/chr,/*coordstart*/pos,/*coordend*/pos,/*sortp*/true);
      if (nmatches > 1) {
	fprintf(stderr,"Multiple intervals for %s:%u\n",chr,pos);
	exit(9);
      } else {
	for (i = 0; i < nmatches; i++) {
	  /* interval = IIT_interval(geneprob_iit,matches[i]); */
	  acc = IIT_label(geneprob_iit,matches[i],&alloc2p);
	  value = atoi(acc);
	  geneprob += value;
	  if (alloc2p == true) {
	    FREE(acc);
	  }
	}
      }
      FREE(matches);
    }
  }

  if (allocp == true) {
    FREE(chr);
  }

  return geneprob;
}


static int
Substring_geneprob_eval_donor (Substring_T donor, IIT_T geneprob_iit, IIT_T chromosome_iit,
			       bool sensep, int querylength) {
  int geneprob = 0, value;
  char *acc;
  int *matches;
  int nmatches, i;

  Genomicpos_T chrpos, pos;
  char *chr;
  bool allocp, alloc2p;

  if (sensep == true) {
    if (donor->plusp == true) {
      chr = IIT_label(chromosome_iit,donor->chrnum,&allocp);
      /* Add 1 to report in 1-based coordinates */
      chrpos = donor->genomicstart - donor->chroffset + 1U;
      for (pos = chrpos; pos <= chrpos + donor->chimera_pos - 1; pos++) {
	matches = IIT_get(&nmatches,geneprob_iit,/*divstring*/chr,/*coordstart*/pos,/*coordend*/pos,/*sortp*/true);
	if (nmatches > 1) {
	  fprintf(stderr,"Multiple intervals for %s:%u\n",chr,pos);
	  exit(9);
	} else {
	  for (i = 0; i < nmatches; i++) {
	    /* interval = IIT_interval(geneprob_iit,matches[i]); */
	    acc = IIT_label(geneprob_iit,matches[i],&alloc2p);
	    value = atoi(acc);
	    geneprob += value;
	    if (alloc2p == true) {
	      FREE(acc);
	    }
	  }
	}
	FREE(matches);
      }

    } else {
      /* Sometimes genomicstart goes into next chromosome, so subtract 1 */
      chr = IIT_label(chromosome_iit,donor->chrnum,&allocp);
      chrpos = donor->genomicstart - donor->chroffset;
      for (pos = chrpos; pos >= chrpos - donor->chimera_pos + 1; pos--) {
	matches = IIT_get(&nmatches,geneprob_iit,/*divstring*/chr,/*coordstart*/pos,/*coordend*/pos,/*sortp*/true);
	if (nmatches > 1) {
	  fprintf(stderr,"Multiple intervals for %s:%u\n",chr,pos);
	  exit(9);
	} else {
	  for (i = 0; i < nmatches; i++) {
	    /* interval = IIT_interval(geneprob_iit,matches[i]); */
	    acc = IIT_label(geneprob_iit,matches[i],&alloc2p);
	    value = atoi(acc);
	    geneprob += value;
	    if (alloc2p == true) {
	      FREE(acc);
	    }
	  }
	}
	FREE(matches);
      }
    }

  } else if (sensep == false) {
    if (donor->plusp == true) {
      chr = IIT_label(chromosome_iit,donor->chrnum,&allocp);
      /* Add 1 to report in 1-based coordinates */
      chrpos = donor->genomicstart - donor->chroffset + 1U;
      for (pos = chrpos + donor->chimera_pos; pos <= chrpos + querylength - 1; pos++) {
	matches = IIT_get(&nmatches,geneprob_iit,/*divstring*/chr,/*coordstart*/pos,/*coordend*/pos,/*sortp*/true);
	if (nmatches > 1) {
	  fprintf(stderr,"Multiple intervals for %s:%u\n",chr,pos);
	  exit(9);
	} else {
	  for (i = 0; i < nmatches; i++) {
	    /* interval = IIT_interval(geneprob_iit,matches[i]); */
	    acc = IIT_label(geneprob_iit,matches[i],&alloc2p);
	    value = atoi(acc);
	    geneprob += value;
	    if (alloc2p == true) {
	      FREE(acc);
	    }
	  }
	}
	FREE(matches);
      }
    } else {
      /* Sometimes genomicstart goes into next chromosome, so subtract 1 */
      chr = IIT_label(chromosome_iit,donor->chrnum,&allocp);
      chrpos = donor->genomicstart - donor->chroffset;
      for (pos = chrpos - donor->chimera_pos; pos >= chrpos - querylength - 1; pos--) {
	matches = IIT_get(&nmatches,geneprob_iit,/*divstring*/chr,/*coordstart*/pos,/*coordend*/pos,/*sortp*/true);
	if (nmatches > 1) {
	  fprintf(stderr,"Multiple intervals for %s:%u\n",chr,pos);
	  exit(9);
	} else {
	  for (i = 0; i < nmatches; i++) {
	    /* interval = IIT_interval(geneprob_iit,matches[i]); */
	    acc = IIT_label(geneprob_iit,matches[i],&alloc2p);
	    value = atoi(acc);
	    geneprob += value;
	    if (alloc2p == true) {
	      FREE(acc);
	    }
	  }
	}
	FREE(matches);
      }
    }
  }

  if (allocp == true) {
    FREE(chr);
  }

  return geneprob;
}


static int
Substring_geneprob_eval_acceptor (Substring_T acceptor, IIT_T geneprob_iit, IIT_T chromosome_iit,
				  bool sensep, int querylength) {
  int geneprob = 0, value;
  char *acc;
  int *matches;
  int nmatches, i;

  Genomicpos_T chrpos, pos;
  char *chr;
  bool allocp, alloc2p;

  if (sensep == true) {
    if (acceptor->plusp == true) {
      chr = IIT_label(chromosome_iit,acceptor->chrnum,&allocp);
      /* Add 1 to report in 1-based coordinates */
      chrpos = acceptor->genomicstart - acceptor->chroffset + 1U;
      for (pos = chrpos + acceptor->chimera_pos; pos <= chrpos + querylength - 1; pos++) {
	matches = IIT_get(&nmatches,geneprob_iit,/*divstring*/chr,/*coordstart*/pos,/*coordend*/pos,/*sortp*/true);
	if (nmatches > 1) {
	  fprintf(stderr,"Multiple intervals for %s:%u\n",chr,pos);
	  exit(9);
	} else {
	  for (i = 0; i < nmatches; i++) {
	    /* interval = IIT_interval(geneprob_iit,matches[i]); */
	    acc = IIT_label(geneprob_iit,matches[i],&alloc2p);
	    value = atoi(acc);
	    geneprob += value;
	    if (alloc2p == true) {
	      FREE(acc);
	    }
	  }
	}
	FREE(matches);
      }
    } else {
      /* Sometimes genomicstart goes into next chromosome, so subtract 1 */
      chr = IIT_label(chromosome_iit,acceptor->chrnum,&allocp);
      chrpos = acceptor->genomicstart - acceptor->chroffset;
      for (pos = chrpos - acceptor->chimera_pos; pos >= chrpos - querylength + 1; pos--) {
	matches = IIT_get(&nmatches,geneprob_iit,/*divstring*/chr,/*coordstart*/pos,/*coordend*/pos,/*sortp*/true);
	if (nmatches > 1) {
	  fprintf(stderr,"Multiple intervals for %s:%u\n",chr,pos);
	  exit(9);
	} else {
	  for (i = 0; i < nmatches; i++) {
	    /* interval = IIT_interval(geneprob_iit,matches[i]); */
	    acc = IIT_label(geneprob_iit,matches[i],&alloc2p);
	    value = atoi(acc);
	    geneprob += value;
	    if (alloc2p == true) {
	      FREE(acc);
	    }
	  }
	}
	FREE(matches);
      }
    }

  } else if (sensep == false) {
    if (acceptor->plusp == true) {
      chr = IIT_label(chromosome_iit,acceptor->chrnum,&allocp);
      /* Add 1 to report in 1-based coordinates */
      chrpos = acceptor->genomicstart - acceptor->chroffset + 1U;
      for (pos = chrpos; pos <= chrpos + acceptor->chimera_pos - 1; pos++) {
	matches = IIT_get(&nmatches,geneprob_iit,/*divstring*/chr,/*coordstart*/pos,/*coordend*/pos,/*sortp*/true);
	if (nmatches > 1) {
	  fprintf(stderr,"Multiple intervals for %s:%u\n",chr,pos);
	  exit(9);
	} else {
	  for (i = 0; i < nmatches; i++) {
	    /* interval = IIT_interval(geneprob_iit,matches[i]); */
	    acc = IIT_label(geneprob_iit,matches[i],&alloc2p);
	    value = atoi(acc);
	    geneprob += value;
	    if (alloc2p == true) {
	      FREE(acc);
	    }
	  }
	}
	FREE(matches);
      }
    } else {
      /* Sometimes genomicstart goes into next chromosome, so subtract 1 */
      chr = IIT_label(chromosome_iit,acceptor->chrnum,&allocp);
      chrpos = acceptor->genomicstart - acceptor->chroffset;
      for (pos = chrpos; pos >= chrpos - acceptor->chimera_pos + 1; pos--) {
	matches = IIT_get(&nmatches,geneprob_iit,/*divstring*/chr,/*coordstart*/pos,/*coordend*/pos,/*sortp*/true);
	if (nmatches > 1) {
	  fprintf(stderr,"Multiple intervals for %s:%u\n",chr,pos);
	  exit(9);
	} else {
	  for (i = 0; i < nmatches; i++) {
	    /* interval = IIT_interval(geneprob_iit,matches[i]); */
	    acc = IIT_label(geneprob_iit,matches[i],&alloc2p);
	    value = atoi(acc);
	    geneprob += value;
	    if (alloc2p == true) {
	      FREE(acc);
	    }
	  }
	}
	FREE(matches);
      }
    }
  }

  if (allocp == true) {
    FREE(chr);
  }

  return geneprob;
}


int
Substring_geneprob_eval_splice (Substring_T donor, Substring_T acceptor, IIT_T geneprob_iit, IIT_T chromosome_iit, int querylength,
				bool sensep) {
  if (acceptor == NULL) {
    /* Single sequence */
    return Substring_geneprob_eval_donor(donor,geneprob_iit,chromosome_iit,sensep,querylength);

  } else if (donor == NULL) {
    /* Single sequence */
    return Substring_geneprob_eval_acceptor(acceptor,geneprob_iit,chromosome_iit,sensep,querylength);

  } else {
    return Substring_geneprob_eval_donor(donor,geneprob_iit,chromosome_iit,sensep,querylength) +
      Substring_geneprob_eval_acceptor(acceptor,geneprob_iit,chromosome_iit,sensep,querylength);
  }
}

