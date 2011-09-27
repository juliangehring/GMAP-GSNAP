static char rcsid[] = "$Id: stage1.c,v 1.108 2005/10/01 15:30:17 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "stage1.h"
#include <stdio.h>
#include <stdlib.h>
#include "mem.h"
#include "reader.h"
#include "block.h"
#include "listdef.h"
#include "matchpair.h"
#include "chrsubset.h"


#ifdef PMAP
#define INDEX1PART INDEX1PART_AA
#endif

#define MININTRONLEN 20
#define VERYSHORTSEQLEN 40
#define SHORTSEQLEN 80

#define EXPANSION 1000		/* (int) Relative to trimlength */

#define MAX_FILL_IN 200
#define MAX_DANGLING_PCT 0.33

#define EXTRA_LONGEND  30000	/* Should exceed INDEX1PART */
#define EXTRA_SHORTEND 10000	/* Should exceed INDEX1PART */
#define SIZEBOUND 0.7
#define MAXMATCHPAIRS_PREUNIQUE 1000
#define MAXMATCHPAIRS_POSTUNIQUE 80

/* Debugging of scanning for 24-mers */
#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Debugging of clusterlist.  May want to turn on DEBUG1 or DEBUG2 in matchpair.c. */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

/* Debugging of 12-mer extension */ 
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

/* Detailed view of scanning for 24-mers */
#ifdef DEBUG3
#define debug3(x) x
#else
#define debug3(x)
#endif

/* Triplet matching */
#ifdef DEBUG4
#define debug4(x) x
#else
#define debug4(x)
#endif

/* connectable_p */
#ifdef DEBUG5
#define debug5(x) x
#else
#define debug5(x)
#endif



#define T Stage1_T
struct T {
  int maxintronlen;
  int querylength;
  int maxentries;

  int stage1size;
  int interval;
  List_T matchpairlist;
  List_T matchbestlist;
  Reader_T reader;
  Block_T block5;
  Block_T block3;
  List_T matches5;
  List_T matches3;
  Genomicpos_T **plus_positions;
  Genomicpos_T **minus_positions;
  int *plus_npositions;
  int *minus_npositions;
  bool *plus_matchedp;		/* For identify_matches */
  bool *minus_matchedp;
  bool *processedp;		/* For Block_process_oligo */
  int querystart;
  int queryend;
};


static T
Stage1_new (Sequence_T queryuc, int maxintronlen, int stage1size, int maxentries, 
	    int reader_overlap) {
  T new = (T) MALLOC(sizeof(*new));
  
  new->maxintronlen = maxintronlen;
  new->maxentries = maxentries;

  new->stage1size = stage1size;
  new->interval = stage1size - INDEX1PART;

  new->querylength = Sequence_fulllength(queryuc);
  new->querystart = Sequence_trim_start(queryuc);
  new->queryend = Sequence_trim_end(queryuc);

  new->matchpairlist = NULL;
  new->matchbestlist = NULL;
  new->reader = Reader_new(Sequence_fullpointer(queryuc),new->querystart,new->queryend,reader_overlap);
  debug(Sequence_print(queryuc,/*uppercasep*/true,/*wraplength*/50,/*trimmedp*/true));
  new->block5 = Block_new(FIVE,new->reader);
  new->block3 = Block_new(THREE,new->reader);
  new->matches5 = NULL;
  new->matches3 = NULL;

  new->plus_positions = (Genomicpos_T **) CALLOC(new->querylength,sizeof(Genomicpos_T *));
  new->minus_positions = (Genomicpos_T **) CALLOC(new->querylength,sizeof(Genomicpos_T *));
  new->plus_npositions = (int *) CALLOC(new->querylength,sizeof(int));
  new->minus_npositions = (int *) CALLOC(new->querylength,sizeof(int));
  new->plus_matchedp = (bool *) CALLOC(new->querylength,sizeof(bool));
  new->minus_matchedp = (bool *) CALLOC(new->querylength,sizeof(bool));
  new->processedp = (bool *) CALLOC(new->querylength,sizeof(bool));

  return new;
}

void
Stage1_free (T *old) {
  List_T p;
  Matchpair_T matchpair;
  Match_T match;
  int i;

  for (i = 0; i < (*old)->querylength; i++) {
    if ((*old)->plus_positions[i] != NULL) {
      FREE((*old)->plus_positions[i]);
    }
    if ((*old)->minus_positions[i] != NULL) {
      FREE((*old)->minus_positions[i]);
    }
  }    

  FREE((*old)->plus_positions);
  FREE((*old)->minus_positions);
  FREE((*old)->plus_npositions);
  FREE((*old)->minus_npositions);
  FREE((*old)->processedp);
  FREE((*old)->plus_matchedp);
  FREE((*old)->minus_matchedp);

  for (p = (*old)->matchbestlist; p != NULL; p = p->rest) {
    matchpair = p->first;
    Matchpair_free(&matchpair);
  }
  List_free(&(*old)->matchbestlist);

  for (p = (*old)->matchpairlist; p != NULL; p = p->rest) {
    matchpair = p->first;
    Matchpair_free(&matchpair);
  }
  List_free(&(*old)->matchpairlist);

  for (p = (*old)->matches5; p != NULL; p = p->rest) {
    match = p->first;
    Match_free(&match);
  }
  List_free(&(*old)->matches5);

  for (p = (*old)->matches3; p != NULL; p = p->rest) {
    match = p->first;
    Match_free(&match);
  }
  List_free(&(*old)->matches3);

  Block_free(&(*old)->block3);
  Block_free(&(*old)->block5);
  Reader_free(&(*old)->reader);

  FREE(*old);
  return;
}



static bool
connectable_p (Match_T match5, Match_T match3, int maxintronlen, int interval) {
  Genomicpos_T position5, position3;
  int querypos5, querypos3, exonlen;
  bool forwardp5, forwardp3;

  debug5(printf("Comparing #%d:%u at querypos %d with #%d:%u at querypos %d\n",
		Match_chrnum(match5),Match_chrpos(match5),Match_querypos(match5),
		Match_chrnum(match3),Match_chrpos(match3),Match_querypos(match3)));

  if (Match_chrnum(match5) != Match_chrnum(match3)) {
    debug5(printf("No, different chromosomes\n\n"));
    return false;
  } else {
    querypos5 = Match_querypos(match5);
    querypos3 = Match_querypos(match3);
    exonlen = querypos3 - querypos5;
    if (exonlen < interval) {
      debug5(printf("No, exonlen == 0\n\n"));
      return false;
    } else {
      position5 = Match_position(match5);
      position3 = Match_position(match3);
      /* intronlen = abs(position3 - position5) - exonlen; -- shouldn't subtract unsigned ints */
      if (position3 > position5) {
	/* intronlen = position3 - position5 - exonlen; -- Don't subtract into a signed int */
	/* The check below is equivalent to intronlen > maxintronlen */
	if (position3 > (Genomicpos_T) maxintronlen + position5 + (Genomicpos_T) exonlen) {
	  debug5(printf("No, intron too long (%u > %u + %u + %u)\n\n",
			position3,maxintronlen,position5,exonlen));
	  return false;
	}
      } else {
	/* intronlen = position5 - position3 - exonlen; -- Don't subtract into a signed int */
	/* The check below is equivalent to intronlen > maxintronlen */
	if (position5 > (Genomicpos_T) maxintronlen + position3 + (Genomicpos_T) exonlen) {
	  debug5(printf("No, intron too long (%u > %u + %u + %u)\n\n",
			position5,maxintronlen,position3,exonlen));
	  return false;
	}
      }
      forwardp5 = Match_forwardp(match5);
      forwardp3 = Match_forwardp(match3);

      if (forwardp5 != forwardp3) {
	debug5(printf("No, forwardp different\n\n"));
	return false;
      } else if (forwardp5 == true && position3 < position5) {
	debug5(printf("No, forwardp true, but positions wrong\n\n"));
	return false;
      } else if (forwardp5 == false && position5 < position3) {
	debug5(printf("No, forwardp false, but positions wrong\n\n"));
	return false;
      } else {
	debug5(printf("Yes\n\n"));
	return true;
      }
    }
  }
}


/* Updates a list of Stage1_T objects */
static List_T
pair_up (bool *foundpairp, List_T matchpairlist, List_T newmatches5, List_T newmatches3, 
	 List_T matches5, List_T matches3, int maxintronlen, int interval) {
  List_T p, q, s, newmatchpairs = NULL;
  Match_T match5, match3;
  Matchpair_T matchpair;
  
  /* Do N vs N */
  for (q = newmatches5; q != NULL; q = q->rest) {
    match5 = q->first;
    for (s = newmatches3; s != NULL; s = s->rest) {
      match3 = s->first;
      if (connectable_p(match5,match3,maxintronlen,interval)) {
	newmatchpairs = List_push(newmatchpairs,Matchpair_new(match5,match3,interval+INDEX1PART,0,MIXED));
	Match_set_pairedp(match5);
	Match_set_pairedp(match3);
      }
    }
  }

  /* Do N vs (N-1..1) */
  for (q = newmatches5; q != NULL; q = q->rest) {
    match5 = q->first;
    for (s = matches3; s != NULL; s = s->rest) {
      match3 = s->first;
      if (connectable_p(match5,match3,maxintronlen,interval)) {
	newmatchpairs = List_push(newmatchpairs,Matchpair_new(match5,match3,interval+INDEX1PART,0,MIXED));
	Match_set_pairedp(match5);
	Match_set_pairedp(match3);
      }
    }
  }

  /* Do (N-1..1) vs N */
  for (q = matches5; q != NULL; q = q->rest) {
    match5 = q->first;
    for (s = newmatches3; s != NULL; s = s->rest) {
      match3 = s->first;
      if (connectable_p(match5,match3,maxintronlen,interval)) {
	newmatchpairs = List_push(newmatchpairs,Matchpair_new(match5,match3,interval+INDEX1PART,0,MIXED));
	Match_set_pairedp(match5);
	Match_set_pairedp(match3);
      }
    }
  }

  if (newmatchpairs == NULL) {
    *foundpairp = false;

  } else if (List_length(newmatchpairs) > MAXMATCHPAIRS_PREUNIQUE) {
    debug(printf("Too many matching pairs before unique (%d > %d)\n",
		 List_length(newmatchpairs),MAXMATCHPAIRS_PREUNIQUE));
    *foundpairp = false;
    for (p = newmatchpairs; p != NULL; p = p->rest) {
      matchpair = (Matchpair_T) List_head(p);
      Matchpair_free(&matchpair);
    }
    List_free(&newmatchpairs);

  } else {
    newmatchpairs = Matchpair_filter_unique(newmatchpairs);
    if (List_length(newmatchpairs) > MAXMATCHPAIRS_POSTUNIQUE) {
      debug(printf("Too many matching pairs (%d > %d)\n",
		   List_length(newmatchpairs),MAXMATCHPAIRS_POSTUNIQUE));
      
      *foundpairp = false;
      for (p = newmatchpairs; p != NULL; p = p->rest) {
	matchpair = (Matchpair_T) List_head(p);
	Matchpair_free(&matchpair);
      }
      List_free(&newmatchpairs);
      
    } else {
      *foundpairp = true;
      for (p = newmatchpairs; p != NULL; p = p->rest) {
	debug(
	      matchpair = (Matchpair_T) List_head(p);
	      match5 = Matchpair_bound5(matchpair);
	      match3 = Matchpair_bound3(matchpair);
	      printf("Found matching pair: #%d:%u (%d)-#%d:%u (%d) for %d..%d\n",
		     Match_chrnum(match5),Match_chrpos(match5),Match_forwardp(match5),
		     Match_chrnum(match3),Match_chrpos(match3),Match_forwardp(match3),
		     Match_querypos(match5),Match_querypos(match3));
	      );
	matchpairlist = List_push(matchpairlist,(Matchpair_T) List_head(p));
      }
      List_free(&newmatchpairs);
    }
  }

  return matchpairlist;
}


static List_T
identify_singles (bool *newp, bool *overflowp, List_T matches, int querypos, 
		  Genomicpos_T *positions, int npositions, 
		  int stage1size, IIT_T chromosome_iit, 
		  Chrsubset_T chrsubset, bool forwardp, bool fivep, int maxentries) {
  List_T newmatches = NULL, p;
  Match_T match;
  Genomicpos_T position;
  int i = 0, nentries = 0;
  bool donep = false;

  if (npositions == 0) {
    *newp = *overflowp = false;
    return matches;
  } else {
    position = positions[0];
  }    

  while (!donep) {
    if (Chrsubset_includep(chrsubset,position,chromosome_iit) == true) {
      if (++nentries > maxentries) {
	donep = true;
      } else if (forwardp == true) {
	newmatches = List_push(newmatches,(void *) Match_new(querypos,true,fivep,position,chromosome_iit));
      } else {
	/* Addition of stage1size needed to make Genome_get_segment
	   symmetric on revcomp of query.  Since we are 1 past the
	   last character, must subtract 1 when looking up in IIT */
	newmatches = List_push(newmatches,(void *) Match_new(querypos,false,fivep,position+stage1size,chromosome_iit));
      }
    }

    /* Advance */
    if (++i >= npositions) {
      donep = true;
    } else {
      position = positions[i];
    }
  }

  if (nentries > maxentries) {
    /* Too many entries.  Not unique enough in genome */
    *newp = false;
    *overflowp = true;
    debug(printf("  Singles overflow at %d\n",querypos));
    for (p = newmatches; p != NULL; p = p->rest) {
      match = p->first;
      Match_free(&match);
    }
    List_free(&newmatches);
  } else {
    *newp = newmatches ? true : false;
    *overflowp = false;
    for (p = newmatches; p != NULL; p = p->rest) {
      debug(
	    printf("  Match at %d: #%d:%u (%d)\n",
		   querypos,Match_chrnum(List_head(p)),
		   Match_chrpos(List_head(p)),Match_forwardp(List_head(p)));
	    );
      matches = List_push(matches,p->first);
    }
    List_free(&newmatches);
  }

  return matches;
}


/* positions0 should be left of positions1, and expecteddist should be
 * positive.  This procedure assumes that the entries in position0 and
 * positions1 are in order.  They are assumed not to have duplicates, which
 * are removed by Indexdb_write_positions and Indexdb_read. */
static List_T
identify_doubles (bool *newp, bool *overflowp, List_T matches, int querypos0, 
		  Genomicpos_T *positions0, int npositions0, 
		  int querypos1, Genomicpos_T *positions1, int npositions1,
		  int stage1size, int expecteddist, IIT_T chromosome_iit, 
		  Chrsubset_T chrsubset, bool forwardp, bool fivep, int maxentries) {
  List_T newmatches = NULL, p;
  Match_T match;
  Genomicpos_T position0, expected1, position1;
  int i = 0, j = 0, nentries = 0;
  bool donep = false;

  if (npositions0 == 0) {
    *newp = *overflowp = false;
    return matches;
  } else {
    position0 = positions0[0];
    expected1 = position0 + expecteddist;
    if (npositions1 == 0) {
      *newp = *overflowp = false;
      return matches;
    } else {
      position1 = positions1[0];
    }
  }

  while (!donep) {
    /* debug(printf("  %d:%u %d:%u\n",i,position0,j,position1)); */
    /*
    debug(
	  if (abs(position1-position0) < 100) {
	    printf("Close: %u %u\n",position0,position1);
	  }
	  );
    */

    if (expected1 < position1) {
      /* Advance position0 */
      if (++i >= npositions0) {
	donep = true;
      } else {
	position0 = positions0[i];
	expected1 = position0 + expecteddist;
      }
    } else if (expected1 > position1) {
      /* Advance position1 */
      if (++j >= npositions1) {
	donep = true;
      } else {
	position1 = positions1[j];
      }

    } else {
      if (Chrsubset_includep(chrsubset,position0,chromosome_iit) == true) {
	if (++nentries > maxentries) {
	  donep = true;
	} else if (forwardp == true) {
	  newmatches = List_push(newmatches,(void *) Match_new(querypos0,true,fivep,position0,chromosome_iit));
	} else {
	  /* Addition of stage1size needed to make Genome_get_segment
	     symmetric on revcomp of query.  Since we are 1 past the
	     last character, must subtract 1 when looking up in IIT */
	  newmatches = List_push(newmatches,(void *) Match_new(querypos1,false,fivep,position0+stage1size,chromosome_iit));
	}
      }

      /* Advance both */
      if (++i >= npositions0) {
	donep = true;
      } else {
	position0 = positions0[i];
	expected1 = position0 + expecteddist;
	if (++j >= npositions1) {
	  donep = true;
	} else {
	  position1 = positions1[j];
	}
      }
    }
  }

  if (nentries > maxentries) {
    /* Too many entries.  Not unique enough in genome */
    *newp = false;
    *overflowp = true;
    debug(printf("  Doubles overflow at %d\n",querypos1));
    for (p = newmatches; p != NULL; p = p->rest) {
      match = p->first;
      Match_free(&match);
    }
    List_free(&newmatches);
  } else {
    *newp = newmatches ? true : false;
    *overflowp = false;
    for (p = newmatches; p != NULL; p = p->rest) {
      debug(
	    printf("  Match at %d: #%d:%u (%d)\n",
		   querypos1,Match_chrnum(List_head(p)),
		   Match_chrpos(List_head(p)),Match_forwardp(List_head(p)));
	    );
      matches = List_push(matches,p->first);
    }
    List_free(&newmatches);
  }

  return matches;
}


/* Involves search of three lists.  This procedure assumes that the
 * entries in positions0, positions1, and positions2 are in order.  They are assumed
 * not to have duplicates, which are removed by
 * Indexdb_write_positions and Indexdb_read. */
static List_T
identify_triples (bool *newp, bool *overflowp, List_T matches, int querypos0, 
		  Genomicpos_T *positions0, int npositions0, 
		  Genomicpos_T *positions1, int npositions1,
		  int querypos2, Genomicpos_T *positions2, int npositions2,
		  int stage1size, int expecteddist1, int expecteddist2, IIT_T chromosome_iit, 
		  Chrsubset_T chrsubset, bool forwardp, bool fivep, int maxentries) {
  List_T newmatches = NULL, p;
  Match_T match;
  Genomicpos_T position0, expected1, position1, expected2, position2;
  int i = 0, j = 0, k = 0, nentries = 0;
  int low2, middle2, high2;
  bool donep = false, foundp;

  if (npositions0 == 0) {
    *newp = *overflowp = false;
    return matches;
  } else {
    position0 = positions0[0];
    expected1 = position0 + expecteddist1;
    if (npositions1 == 0) {
      *newp = *overflowp = false;
      return matches;
    } else {
      position1 = positions1[0];
      expected2 = position1 + expecteddist2;
      if (npositions2 == 0) {
	*newp = *overflowp = false;
	return matches;
      } else {
	position2 = positions2[0];
      }
    }
  }
      
  while (!donep) {
    debug4(printf("  %d:%u %d:%u %d:%u\n",
		  i,position0,j,position1,k,position2));
    if (expected1 < position1) {
      /* Advance position0 */
      if (++i >= npositions0) {
	donep = true;
      } else {
	position0 = positions0[i];
	expected1 = position0 + expecteddist1;
      }
    } else if (expected1 > position1) {
      /* Advance position1 */
      if (++j >= npositions1) {
	donep = true;
      } else {
	position1 = positions1[j];
	expected2 = position1 + expecteddist2;
      }
    } else if (expected2 < position2) {
      /* Advance position1 */
      if (++j >= npositions1) {
	donep = true;
      } else {
	position1 = positions1[j];
	expected2 = position1 + expecteddist2;
      }
    } else if (expected2 > position2) {
      /* Binary search */
      low2 = k+1;
      high2 = npositions2;
      foundp = false;
      debug4(
	     for (middle2 = low2; middle2 < high2; middle2++) {
	       printf("  --%d:%u %d:%u %d:%u\n",
		      i,position0,j,position1,middle2,
		      positions2[middle2]);
	     }
	     );

      while (!foundp && low2 < high2) {
	middle2 = (low2+high2)/2;
	debug4(printf("  **%d:%u %d:%u %d:%u   vs. %u\n",
		      low2,positions2[low2],
		      middle2,positions2[middle2],
		      high2,positions2[high2],
		      expected2));
	if (expected2 < positions2[middle2]) {
	  high2 = middle2;
	} else if (expected2 > positions2[middle2]) {
	  low2 = middle2 + 1;
	} else {
	  foundp = true;
	}
      }
      if (foundp == true) {
	k = middle2;
	position2 = positions2[k];
      } else if ((k = high2) >= npositions2) {
	donep = true;
      } else {
	position2 = positions2[k];
      }

    } else {
      debug4(printf("Success at position %u\n",positions0[i]));
      if (Chrsubset_includep(chrsubset,position0,chromosome_iit) == true) {
	if (++nentries > maxentries) {
	  donep = true;
	} else if (forwardp == true) {
	  newmatches = List_push(newmatches,(void *) Match_new(querypos0,true,fivep,position0,chromosome_iit));
	} else {
	  newmatches = List_push(newmatches,(void *) Match_new(querypos2,false,fivep,position0+stage1size,chromosome_iit));
	}
      }

      /* Advance all */
      if (++i >= npositions0) {
	donep = true;
      } else {
	position0 = positions0[i];
	expected1 = position0 + expecteddist1;
	if (++j >= npositions1) {
	  donep = true;
	} else {
	  position1 = positions1[j];
	  expected2 = position1 + expecteddist2;
	  if (++k >= npositions2) {
	    donep = true;
	  } else {
	    position2 = positions2[k];
	  }
	}
      }
    }
  }

  if (nentries > maxentries) {
    *newp = false;
    *overflowp = true;
    debug(printf("  Triples overflow at %d\n",querypos2));
    for (p = newmatches; p != NULL; p = p->rest) {
      match = p->first;
      Match_free(&match);
    }
    List_free(&newmatches);
  } else {
    *newp = newmatches ? true : false;
    *overflowp = false;
    for (p = newmatches; p != NULL; p = p->rest) {
      matches = List_push(matches,p->first);
    }
    List_free(&newmatches);
  }

  return matches;
}


static List_T
identify_matches (bool *newp, bool *overflowp, List_T matches, int querypos,
		  int stage1size, int interval, Genomicpos_T **plus_positions, int *plus_npositions,
		  Genomicpos_T **minus_positions, int *minus_npositions, 
		  IIT_T chromosome_iit, Chrsubset_T chrsubset, bool forwardp, bool fivep, int maxentries) {
  int prevpos, middlepos;

  if (fivep == true) {
    prevpos = querypos - interval;
    if (forwardp == true) {
      /* fivep == true and forwardp == true */
      if (stage1size == INDEX1PART) {
	return identify_singles(&(*newp),&(*overflowp),matches,querypos,plus_positions[querypos],plus_npositions[querypos],
				stage1size,chromosome_iit,chrsubset,
				/*forwardp*/true,/*fivep*/true,maxentries);
      } else if (stage1size <= 2*INDEX1PART) {
	return identify_doubles(&(*newp),&(*overflowp),matches,prevpos,plus_positions[prevpos],plus_npositions[prevpos],
				querypos,plus_positions[querypos],plus_npositions[querypos],
				stage1size,interval,chromosome_iit,chrsubset,
				/*forwardp*/true,/*fivep*/true,maxentries);
      } else if (stage1size == 3*INDEX1PART) {
	middlepos = querypos - INDEX1PART;
	return identify_triples(&(*newp),&(*overflowp),matches,prevpos,plus_positions[prevpos],plus_npositions[prevpos],
				plus_positions[middlepos],plus_npositions[middlepos],
				querypos,plus_positions[querypos],plus_npositions[querypos],
				stage1size,INDEX1PART,INDEX1PART,chromosome_iit,chrsubset,
				/*forwardp*/true,/*fivep*/true,maxentries);
      } else {
	abort();
      }

    } else {
      /* fivep == true and forwardp == false */
      if (stage1size == INDEX1PART) {
	return identify_singles(&(*newp),&(*overflowp),matches,querypos,minus_positions[querypos],minus_npositions[querypos],
				stage1size,chromosome_iit,chrsubset,
				/*forwardp*/false,/*fivep*/true,maxentries);

      } else if (stage1size <= 2*INDEX1PART) {
	return identify_doubles(&(*newp),&(*overflowp),matches,querypos,minus_positions[querypos],minus_npositions[querypos],
				prevpos,minus_positions[prevpos],minus_npositions[prevpos],
				stage1size,interval,chromosome_iit,chrsubset,
				/*forwardp*/false,/*fivep*/true,maxentries);
      } else if (stage1size == 3*INDEX1PART) {
	middlepos = querypos - INDEX1PART;
	return identify_triples(&(*newp),&(*overflowp),matches,querypos,minus_positions[querypos],minus_npositions[querypos],
				minus_positions[middlepos],minus_npositions[middlepos],
				prevpos,minus_positions[prevpos],minus_npositions[prevpos],
				stage1size,INDEX1PART,INDEX1PART,chromosome_iit,chrsubset,
				/*forwardp*/false,/*fivep*/true,maxentries);

      } else {
	abort();
      }

    }

  } else {
    prevpos = querypos + interval;
    if (forwardp == true) {
      /* fivep == false and forwardp == true */
      if (stage1size == INDEX1PART) {
	return identify_singles(&(*newp),&(*overflowp),matches,querypos,plus_positions[querypos],plus_npositions[querypos],
				stage1size,chromosome_iit,chrsubset,
				/*forwardp*/true,/*fivep*/false,maxentries);

      } else if (stage1size <= 2*INDEX1PART) {
	return identify_doubles(&(*newp),&(*overflowp),matches,querypos,plus_positions[querypos],plus_npositions[querypos],
				prevpos,plus_positions[prevpos],plus_npositions[prevpos],
				stage1size,interval,chromosome_iit,chrsubset,
				/*forwardp*/true,/*fivep*/false,maxentries);
      } else if (stage1size == 3*INDEX1PART) {
	middlepos = querypos + INDEX1PART;
	return identify_triples(&(*newp),&(*overflowp),matches,querypos,plus_positions[querypos],plus_npositions[querypos],
				plus_positions[middlepos],plus_npositions[middlepos],
				prevpos,plus_positions[prevpos],plus_npositions[prevpos],
				stage1size,INDEX1PART,INDEX1PART,chromosome_iit,chrsubset,
				/*forwardp*/true,/*fivep*/false,maxentries);
      } else {
	abort();
      }

    } else {
      /* fivep == false and forwardp == false */
      if (stage1size == INDEX1PART) {
	return identify_singles(&(*newp),&(*overflowp),matches,prevpos,minus_positions[prevpos],minus_npositions[prevpos],
				stage1size,chromosome_iit,chrsubset,
				/*forwardp*/false,/*fivep*/false,maxentries);
      } else if (stage1size <= 2*INDEX1PART) {
	return identify_doubles(&(*newp),&(*overflowp),matches,prevpos,minus_positions[prevpos],minus_npositions[prevpos],
				querypos,minus_positions[querypos],minus_npositions[querypos],
				stage1size,interval,chromosome_iit,chrsubset,
				/*forwardp*/false,/*fivep*/false,maxentries);
      } else if (stage1size == 3*INDEX1PART) {
	middlepos = querypos + INDEX1PART;
	return identify_triples(&(*newp),&(*overflowp),matches,prevpos,minus_positions[prevpos],minus_npositions[prevpos],
				minus_positions[middlepos],minus_npositions[middlepos],
				querypos,minus_positions[querypos],minus_npositions[querypos],
				stage1size,INDEX1PART,INDEX1PART,chromosome_iit,chrsubset,
				/*forwardp*/false,/*fivep*/false,maxentries);
	
      } else {
	abort();
      }
    }
  }
}


/* Need to change these procedures from n*m to n+m */
static List_T
find_5prime_matches (bool *newp, List_T matches5, T this, Genomicpos_T **plus_positions, int *plus_npositions,
		     Genomicpos_T **minus_positions, int *minus_npositions,
		     IIT_T chromosome_iit, Chrsubset_T chrsubset, int querystart, 
		     int interval, int maxentries, bool pairedp) {
  int prevpos;			/* Note: negative values must be allowed */
  bool newplusp = false, newminusp = false, overflowp;

  if ((prevpos = querystart - interval) <= 0) {
    debug3(printf("Quitting because %d - %d <= 0\n",querystart,interval));
    *newp = false;
  } else {
    debug3(printf("Identifying 5' matches for %d-mer at %d and %d-mer at %d.  maxentries=%d.\n",
		  INDEX1PART,prevpos,INDEX1PART,querystart,maxentries));

    debug3(printf("  Previous status of plus_matchedp at %d is %d\n",querystart,this->plus_matchedp[querystart]));
    if (pairedp == true || this->plus_matchedp[querystart] == false) {
      matches5 = identify_matches(&newplusp,&overflowp,matches5,querystart,this->stage1size,interval,
				  plus_positions,plus_npositions,minus_positions,minus_npositions,
				  chromosome_iit,chrsubset,true,true,maxentries);
      if (overflowp == false) {
	debug3(printf("  No overflow so setting plus_matchedp at %d to true\n",querystart));
	this->plus_matchedp[querystart] = true;
      }
    }
    debug3(printf("  Previous status of minus_matchedp at %d is %d\n",prevpos,this->minus_matchedp[prevpos]));
    if (pairedp == true || this->minus_matchedp[prevpos] == false) {
      matches5 = identify_matches(&newminusp,&overflowp,matches5,querystart,this->stage1size,interval,
				  plus_positions,plus_npositions,minus_positions,minus_npositions,
				  chromosome_iit,chrsubset,false,true,maxentries);
      if (overflowp == false) {
	debug3(printf("  No overflow so setting minus_matchedp at %d to true\n",prevpos));
	this->minus_matchedp[prevpos] = true;
      }
    }
    *newp = newplusp || newminusp;
  }

  return matches5;
}

static List_T
find_3prime_matches (bool *newp, List_T matches3, T this, Genomicpos_T **plus_positions, int *plus_npositions,
		     Genomicpos_T **minus_positions, int *minus_npositions,
		     IIT_T chromosome_iit, Chrsubset_T chrsubset, int queryend, int querylength, 
		     int interval, int maxentries, bool pairedp) {
  int prevpos;			/* Note: negative values must be allowed */
  bool newplusp = false, newminusp = false, overflowp;

  if ((prevpos = queryend + interval) >= querylength) {
    *newp = false;
  } else {
    debug3(printf("Identifying 3' matches for %d-mer at %d and %d-mer at %d.  maxentries=%d\n",
		  INDEX1PART,queryend,INDEX1PART,prevpos,maxentries));

    /* If maxentries > 0, then we are working on pairs.  Otherwise, check for a successful previous check */
    debug3(printf("  Previous status of plus_matchedp at %d is %d\n",prevpos,this->plus_matchedp[prevpos]));
    if (pairedp == true || this->plus_matchedp[prevpos] == false) {

      matches3 = identify_matches(&newplusp,&overflowp,matches3,queryend,this->stage1size,interval,
				  plus_positions,plus_npositions,minus_positions,minus_npositions,
				  chromosome_iit,chrsubset,true,false,maxentries);
      if (overflowp == false) {
	debug3(printf("  No overflow so setting plus_matchedp at %d to true\n",prevpos));
	this->plus_matchedp[prevpos] = true;
      }
    }

    debug3(printf("  Previous status of minus_matchedp at %d is %d\n",queryend,this->minus_matchedp[queryend]));
    if (pairedp == true || this->minus_matchedp[queryend] == false) {
      matches3 = identify_matches(&newminusp,&overflowp,matches3,queryend,this->stage1size,interval,
				  plus_positions,plus_npositions,minus_positions,minus_npositions,
				  chromosome_iit,chrsubset,false,false,maxentries);
      if (overflowp == false) {
	debug3(printf("  No overflow so setting minus_matchedp at %d to true\n",queryend));
	this->minus_matchedp[queryend] = true;
      }
    }

    *newp = newplusp || newminusp;
  }

  return matches3;
}


/*
static bool
check_fraction_paired (List_T matches5, List_T matches3) {
  List_T p;
  int npaired5 = 0, nsingle5 = 0, npaired3 = 0, nsingle3 = 0;
  double fracpaired5, fracpaired3;

  for (p = matches5; p != NULL; p = p->rest) {
    if (Match_pairedp((Match_T) List_head(p)) == true) {
      npaired5++;
      } else {
	nsingle5++;
      }
  }
  for (p = matches3; p != NULL; p = p->rest) {
    if (Match_pairedp((Match_T) List_head(p)) == true) {
      npaired3++;
    } else {
      nsingle3++;
    }
  }

  debug(printf("npaired5: %d, nsingle5: %d, npaired3: %d, nsingle3: %d\n",
	       npaired5,nsingle5,npaired3,nsingle3);
	);
  if (npaired5+nsingle5 == 0) {
    fracpaired5 = (double) (npaired5+1)/(double) (npaired5+nsingle5+1);
  } else {
    fracpaired5 = (double) npaired5/(double) (npaired5+nsingle5);
  }
  if (npaired3+nsingle3 == 0) {
    fracpaired3 = (double) npaired3/(double) (npaired3+nsingle3+1);
  } else {
    fracpaired3 = (double) npaired3/(double) (npaired3+nsingle3);
  }
  if (fracpaired5 > 0.75 || fracpaired3 > 0.75) {
    return true;
  } else {
    return false;
  }
}
*/

static void
stutter (T this, 
#ifdef PMAP
	 Indexdb_T indexdb_fwd,
	 Indexdb_T indexdb_rev,
#else
	 Indexdb_T indexdb, 
#endif
	 IIT_T chromosome_iit, Chrsubset_T chrsubset, int maxintronlen, 
	 int stuttercycles, int stutterhits, int stage1size, int maxentries) {
  List_T newmatches5 = NULL, newmatches3 = NULL;
  int stutterdist5 = 0, stutterdist3 = 0, maxbases, start5, start3, n5hits = 0, n3hits = 0;
  bool newp, foundpairp;
  Match_T match;

  start5 = Block_querypos(this->block5);
  start3 = Block_querypos(this->block3);
  maxbases = stuttercycles*stage1size;
  if (maxbases > (start3 - start5)/2) {
    maxbases = (start3 - start5)/2;
  }
  debug(printf("*** Beginning stutter.  stuttercycles = %d, maxbases = %d ***\n",stuttercycles,maxbases));

  while (Block_next(this->block5) == true &&
	 (stutterdist5 < maxbases || n5hits < stutterhits)) {
    this->querystart = Block_querypos(this->block5);
    /* For stutter, shouldn't need to check for prior this->processedp */
#ifdef PMAP
    Block_process_oligo(&(this->plus_positions[this->querystart]),&(this->plus_npositions[this->querystart]),
			&(this->minus_positions[this->querystart]),&(this->minus_npositions[this->querystart]),
			this->block5,indexdb_fwd,indexdb_rev);
#else
    Block_process_oligo(&(this->plus_positions[this->querystart]),&(this->plus_npositions[this->querystart]),
			&(this->minus_positions[this->querystart]),&(this->minus_npositions[this->querystart]),
			this->block5,indexdb);
#endif
    this->processedp[this->querystart] = true;
    newmatches5 = find_5prime_matches(&newp,newmatches5,this,this->plus_positions,this->plus_npositions,
				      this->minus_positions,this->minus_npositions,
				      chromosome_iit,chrsubset,this->querystart,
				      this->interval,maxentries,true);

    stutterdist5 = Block_querypos(this->block5) - start5;
    if (newp) {
      n5hits++;
      debug(printf("n5hits = %d, stutterdist5 = %d\n",n5hits,stutterdist5));
    }
  }

  while (Block_next(this->block3) == true && 
	 (stutterdist3 < maxbases || n3hits < stutterhits)) {
    this->queryend = Block_querypos(this->block3);
    /* For stutter, shouldn't need to check for prior this->processedp */
#ifdef PMAP
    Block_process_oligo(&(this->plus_positions[this->queryend]),&(this->plus_npositions[this->queryend]),
			&(this->minus_positions[this->queryend]),&(this->minus_npositions[this->queryend]),
			this->block3,indexdb_fwd,indexdb_rev);
#else
    Block_process_oligo(&(this->plus_positions[this->queryend]),&(this->plus_npositions[this->queryend]),
			&(this->minus_positions[this->queryend]),&(this->minus_npositions[this->queryend]),
			this->block3,indexdb);
#endif
    this->processedp[this->queryend] = true;
    newmatches3 = find_3prime_matches(&newp,newmatches3,this,this->plus_positions,this->plus_npositions,
				      this->minus_positions,this->minus_npositions,
				      chromosome_iit,chrsubset,this->queryend,
				      this->querylength,this->interval,maxentries,true);

    stutterdist3 = start3 - Block_querypos(this->block3);
    if (newp) {
      n3hits++;
      debug(printf("n3hits = %d, stutterdist3 = %d\n",n3hits,stutterdist3));
    }
  }

  debug(printf("*** Ending stutter ***\n"));

  this->matchpairlist = pair_up(&foundpairp,this->matchpairlist,newmatches5,newmatches3,
				this->matches5,this->matches3,maxintronlen,this->interval);

  while (newmatches5 != NULL) {
    newmatches5 = List_pop(newmatches5,(void **) &match);
    this->matches5 = List_push(this->matches5,(void *) match);
  }
  while (newmatches3 != NULL) {
    newmatches3 = List_pop(newmatches3,(void **) &match);
    this->matches3 = List_push(this->matches3,(void *) match);
  }

  return;
}


/* Tries to find matches on the 5' end to unpaired hits from the 3' end */
static void
fill_in_5 (T this, List_T dangling3,
#ifdef PMAP
	   Indexdb_T indexdb_fwd,
	   Indexdb_T indexdb_rev,
#else
	   Indexdb_T indexdb, 
#endif
	   IIT_T chromosome_iit, Chrsubset_T chrsubset, int maxintronlen, int maxentries) {
  List_T newmatches5 = NULL;
  int fillindist5 = 0, maxbases, start5;
  bool newp, foundpairp = false;
  Match_T match;

  start5 = Block_querypos(this->block5);
  maxbases = MAX_FILL_IN;
  if (maxbases > this->querylength/2 - start5) {
    maxbases = this->querylength/2 - start5;
  }
  debug(printf("*** Beginning fill_in_5.  maxbases = %d ***\n",maxbases));

  while (Block_next(this->block5) == true && 
	 fillindist5 < maxbases && foundpairp == false) {
    this->querystart = Block_querypos(this->block5);
#ifdef PMAP
    Block_process_oligo(&(this->plus_positions[this->querystart]),&(this->plus_npositions[this->querystart]),
			&(this->minus_positions[this->querystart]),&(this->minus_npositions[this->querystart]),
			this->block5,indexdb_fwd,indexdb_rev);
#else
    Block_process_oligo(&(this->plus_positions[this->querystart]),&(this->plus_npositions[this->querystart]),
			&(this->minus_positions[this->querystart]),&(this->minus_npositions[this->querystart]),
			this->block5,indexdb);
#endif
    this->processedp[this->querystart] = true;
    newmatches5 = find_5prime_matches(&newp,newmatches5,this,this->plus_positions,this->plus_npositions,
				      this->minus_positions,this->minus_npositions,
				      chromosome_iit,chrsubset,this->querystart,
				      this->interval,maxentries,true);
    fillindist5 = Block_querypos(this->block5) - start5;

    if (newp) {
      this->matchpairlist = pair_up(&foundpairp,this->matchpairlist,newmatches5,(List_T) NULL,
				    (List_T) NULL,dangling3,maxintronlen,this->interval);
      debug(printf("   Foundpairp = %d\n",foundpairp));
    }
  }

  /* Mark newmatches5 as being pairedp if they match non-dangling matches3 */
  this->matchpairlist = pair_up(&foundpairp,this->matchpairlist,newmatches5,(List_T) NULL,
				(List_T) NULL,this->matches3,maxintronlen,this->interval);

  while (newmatches5 != NULL) {
    newmatches5 = List_pop(newmatches5,(void **) &match);
    this->matches5 = List_push(this->matches5,(void *) match);
  }

  debug(printf("*** Ending fill_in_5 ***\n"));
  return;
}


/* Tries to find matches on the 5' end to unpaired hits from the 3' end */
static void
fill_in_3 (T this, List_T dangling5, 
#ifdef PMAP
	   Indexdb_T indexdb_fwd,
	   Indexdb_T indexdb_rev,
#else
	   Indexdb_T indexdb,
#endif
	   IIT_T chromosome_iit, Chrsubset_T chrsubset, int maxintronlen, int maxentries) {
  List_T newmatches3 = NULL;
  int fillindist3 = 0, maxbases, start3;
  bool newp, foundpairp = false;
  Match_T match;

  start3 = Block_querypos(this->block3);
  maxbases = MAX_FILL_IN;
  if (maxbases > start3 - this->querylength/2) {
    maxbases = start3 - this->querylength/2;
  }
  debug(printf("*** Beginning fill_in_3.  maxbases = %d ***\n",maxbases));

  while (Block_next(this->block3) == true && 
	 fillindist3 < maxbases && foundpairp == false) {
    this->queryend = Block_querypos(this->block3);
    /* For stutter, shouldn't need to check for prior this->processedp */
#ifdef PMAP
    Block_process_oligo(&(this->plus_positions[this->queryend]),&(this->plus_npositions[this->queryend]),
			&(this->minus_positions[this->queryend]),&(this->minus_npositions[this->queryend]),
			this->block3,indexdb_fwd,indexdb_rev);
#else
    Block_process_oligo(&(this->plus_positions[this->queryend]),&(this->plus_npositions[this->queryend]),
			&(this->minus_positions[this->queryend]),&(this->minus_npositions[this->queryend]),
			this->block3,indexdb);
#endif
    this->processedp[this->queryend] = true;
    newmatches3 = find_3prime_matches(&newp,newmatches3,this,this->plus_positions,this->plus_npositions,
				      this->minus_positions,this->minus_npositions,
				      chromosome_iit,chrsubset,this->queryend,
				      this->querylength,this->interval,maxentries,true);
    fillindist3 = start3 - Block_querypos(this->block3);

    if (newp) {
      this->matchpairlist = pair_up(&foundpairp,this->matchpairlist,(List_T) NULL,newmatches3,
				    dangling5,(List_T) NULL,maxintronlen,this->interval);
      debug(printf("   Foundpairp = %d\n",foundpairp));
    }
  }

  /* Mark newmatches3 as being pairedp if they match non-dangling matches5 */
  this->matchpairlist = pair_up(&foundpairp,this->matchpairlist,(List_T) NULL,newmatches3,
				this->matches5,(List_T) NULL,maxintronlen,this->interval);

  while (newmatches3 != NULL) {
    newmatches3 = List_pop(newmatches3,(void **) &match);
    this->matches3 = List_push(this->matches3,(void *) match);
  }

  debug(printf("*** Ending fill_in_3 ***\n"));
  return;
}



static void
sample (T this, int nskip, int maxentries2, 
#ifdef PMAP
	Indexdb_T indexdb_fwd,
	Indexdb_T indexdb_rev,
#else
	Indexdb_T indexdb,
#endif
	IIT_T chromosome_iit, Chrsubset_T chrsubset) {
  List_T newmatches5 = NULL, newmatches3 = NULL;
  int n5hits = 0, n3hits = 0;
  bool donep = false, newp;
  Match_T match;

  Block_reset_ends(this->block5);
  Block_reset_ends(this->block3);

  while (!donep) {
    if (n5hits <= n3hits) {
      if (Block_next(this->block5) == false) {
	donep = true;
      } else {
	this->querystart = Block_querypos(this->block5);
	if (this->processedp[this->querystart] == false) {
#ifdef PMAP
	  Block_process_oligo(&(this->plus_positions[this->querystart]),&(this->plus_npositions[this->querystart]),
			      &(this->minus_positions[this->querystart]),&(this->minus_npositions[this->querystart]),
			      this->block5,indexdb_fwd,indexdb_rev);
#else
	  Block_process_oligo(&(this->plus_positions[this->querystart]),&(this->plus_npositions[this->querystart]),
			      &(this->minus_positions[this->querystart]),&(this->minus_npositions[this->querystart]),
			      this->block5,indexdb);
#endif
	  this->processedp[this->querystart] = true;
	}
	newmatches5 = find_5prime_matches(&newp,newmatches5,this,this->plus_positions,this->plus_npositions,
					  this->minus_positions,this->minus_npositions,
					  chromosome_iit,chrsubset,this->querystart,this->interval,
					  maxentries2,false);
	if (newp) {
	  n5hits++;
	  debug3(printf("    n5hits: %d, n3hits: %d\n",n5hits,n3hits));
	  while (newmatches5 != NULL) {
	    newmatches5 = List_pop(newmatches5,(void **) &match);
	    this->matches5 = List_push(this->matches5,(void *) match);
	  }
	  Block_skip(this->block5,nskip);
	}
      }
    } else {
      if (Block_next(this->block3) == false) {
	donep = true;
      } else {
	this->queryend = Block_querypos(this->block3);
	if (this->processedp[this->queryend] == false) {
#ifdef PMAP
	  Block_process_oligo(&(this->plus_positions[this->queryend]),&(this->plus_npositions[this->queryend]),
			      &(this->minus_positions[this->queryend]),&(this->minus_npositions[this->queryend]),
			      this->block3,indexdb_fwd,indexdb_rev);
#else
	  Block_process_oligo(&(this->plus_positions[this->queryend]),&(this->plus_npositions[this->queryend]),
			      &(this->minus_positions[this->queryend]),&(this->minus_npositions[this->queryend]),
			      this->block3,indexdb);
#endif
	  this->processedp[this->queryend] = true;
	}
	newmatches3 = find_3prime_matches(&newp,newmatches3,this,this->plus_positions,this->plus_npositions,
					  this->minus_positions,this->minus_npositions,
					  chromosome_iit,chrsubset,this->queryend,
					  this->querylength,this->interval,maxentries2,false);
	if (newp) {
	  n3hits++;
	  debug3(printf("    n5hits: %d, n3hits: %d\n",n5hits,n3hits));
	  while (newmatches3 != NULL) {
	    newmatches3 = List_pop(newmatches3,(void **) &match);
	    this->matches3 = List_push(this->matches3,(void *) match);
	  }
	  Block_skip(this->block3,nskip);
	}
      }
    }
  }

  return;
}


void
Stage1_find_extensions (Genomicpos_T *extension5, Genomicpos_T *extension3, T this, 
			Match_T match5, Match_T match3, int maxextension) {
  int querystart, queryend, newstart, newend, i, j;
  Genomicpos_T genomicstart, genomicend, position;
  int delta;
  
  querystart = Match_querypos(match5);
  queryend = Match_querypos(match3);
  genomicstart = Match_position(match5);
  genomicend = Match_position(match3);

  *extension5 = newstart = querystart;
  debug2(printf("Extending %d - %d\n",querystart - 1,0));
  for (i = querystart - 1; i >= 0; i--) {
    if (Match_forwardp(match5) == true) {
      for (j = 0; j < this->plus_npositions[i]; j++) {
	position = this->plus_positions[i][j];
	if (position < genomicstart && (delta = genomicstart - position) < maxextension && delta > 0) {
	  debug2(printf("%d +%u %d\n",i,position,delta));
	  if (delta + i >= *extension5) {
	    *extension5 = delta + i;
	    newstart = i;
	  }
	}
      }
    } else {
      for (j = 0; j < this->minus_npositions[i]; j++) {
	position = this->minus_positions[i][j] + INDEX1PART;
	if (position > genomicstart && (delta = position - genomicstart) < maxextension && delta > 0) {
	  debug2(printf("%d -%u %d\n",i,position,delta));
	  if (delta + i >= *extension5) {
	    *extension5 = delta + i;
	    newstart = i;
	  }
	}
      }
    }
  }
  if (newstart > INDEX1PART) {
    *extension5 += EXTRA_LONGEND;
  } else {
    *extension5 += EXTRA_SHORTEND;
  }
  debug2(printf("Extension5 = %d\n",*extension5));

  *extension3 = newend = this->querylength - queryend - INDEX1PART;
  debug2(printf("Extending %d - %d\n",queryend+1,this->querylength));
  for (i = queryend + 1; i < this->querylength; i++) {
    if (Match_forwardp(match3) == true) {
      for (j = 0; j < this->plus_npositions[i]; j++) {
	position = this->plus_positions[i][j];
	if (position > genomicend && (delta = position - genomicend) < maxextension && delta > 0) {
	  debug2(printf("%d +%u %d\n",i,position,delta));
	  if (delta + this->querylength - i - INDEX1PART >= *extension3) {
	    *extension3 = delta + this->querylength - i - INDEX1PART;
	    newend = this->querylength - i - INDEX1PART;
	  }
	}
      }
    } else {
      for (j = 0; j < this->minus_npositions[i]; j++) {
	position = this->minus_positions[i][j] + INDEX1PART;
	if (position < genomicend && (delta = genomicend - position) < maxextension && delta > 0) {
	  debug2(printf("%d -%u %d\n",i,position,delta));
	  if (delta + this->querylength - i - INDEX1PART >= *extension3) {
	    *extension3 = delta + this->querylength - i - INDEX1PART;
	    newend = this->querylength - i - INDEX1PART;
	  }
	}
      }
    }
  }

  if (newend > INDEX1PART) {
    *extension3 += EXTRA_LONGEND;
  } else {
    *extension3 += EXTRA_SHORTEND;
  }
  debug2(printf("Extension3 = %d\n",*extension3));

  return;
}



static bool
find_first_pair (T this, 
#ifdef PMAP
		 Indexdb_T indexdb_fwd,
		 Indexdb_T indexdb_rev,
#else
		 Indexdb_T indexdb,
#endif
		 IIT_T chromosome_iit, Chrsubset_T chrsubset, int maxintronlen, int maxentries) {
  List_T newmatches5 = NULL, newmatches3 = NULL;
  Match_T match;
  bool donep = false, newp, foundpairp = false;
  double n5hits = 0.0, n3hits = 0.0;

  debug(printf("*** Starting stage1 ***\n"));
  debug(printf("Reader querystart is %d and queryend is %d\n",
	       Reader_querystart(this->reader),Reader_queryend(this->reader)));
  while (!donep && foundpairp == false) {
    if (n5hits <= n3hits) {
      if (Block_next(this->block5) == false) {
	donep = true;
      } else {
	this->querystart = Block_querypos(this->block5);
	/* For initial scan, shouldn't need to check for prior this->processedp */
#ifdef PMAP
	Block_process_oligo(&(this->plus_positions[this->querystart]),&(this->plus_npositions[this->querystart]),
			    &(this->minus_positions[this->querystart]),&(this->minus_npositions[this->querystart]),
			    this->block5,indexdb_fwd,indexdb_rev);
#else
	Block_process_oligo(&(this->plus_positions[this->querystart]),&(this->plus_npositions[this->querystart]),
			    &(this->minus_positions[this->querystart]),&(this->minus_npositions[this->querystart]),
			    this->block5,indexdb);
#endif
	this->processedp[this->querystart] = true;
	newmatches5 = find_5prime_matches(&newp,NULL,this,this->plus_positions,this->plus_npositions,
					  this->minus_positions,this->minus_npositions,chromosome_iit,
					  chrsubset,this->querystart,this->interval,maxentries,true);
	if (newp) {
	  n5hits += 1.0/(double) (1 + List_length(newmatches5));
	  debug(printf("    n5hits: %.1f, n3hits: %.1f\n",n5hits,n3hits));
	  this->matchpairlist = pair_up(&foundpairp,this->matchpairlist,newmatches5,NULL,
					this->matches5,this->matches3,maxintronlen,this->interval);
	  while (newmatches5 != NULL) {
	    newmatches5 = List_pop(newmatches5,(void **) &match);
	    this->matches5 = List_push(this->matches5,(void *) match);
	  }
	}
      }

    } else {
      if (Block_next(this->block3) == false) {
	donep = true;
      } else {
	this->queryend = Block_querypos(this->block3);
	/* For initial scan, shouldn't need to check for prior this->processedp */
#ifdef PMAP
	Block_process_oligo(&(this->plus_positions[this->queryend]),&(this->plus_npositions[this->queryend]),
			    &(this->minus_positions[this->queryend]),&(this->minus_npositions[this->queryend]),
			    this->block3,indexdb_fwd,indexdb_rev);
#else
	Block_process_oligo(&(this->plus_positions[this->queryend]),&(this->plus_npositions[this->queryend]),
			    &(this->minus_positions[this->queryend]),&(this->minus_npositions[this->queryend]),
			    this->block3,indexdb);
#endif
	this->processedp[this->queryend] = true;
	newmatches3 = find_3prime_matches(&newp,NULL,this,this->plus_positions,this->plus_npositions,
					  this->minus_positions,this->minus_npositions,
					  chromosome_iit,chrsubset,this->queryend,this->querylength,
					  this->interval,maxentries,true);
	if (newp) {
	  n3hits += 1.0/(double) (1 + List_length(newmatches3));
	  debug(printf("    n5hits: %.1f, n3hits: %.1f\n",n5hits,n3hits));
	  this->matchpairlist = pair_up(&foundpairp,this->matchpairlist,NULL,newmatches3,
					this->matches5,this->matches3,maxintronlen,this->interval);
	  while (newmatches3 != NULL) {
	    newmatches3 = List_pop(newmatches3,(void **) &match);
	    this->matches3 = List_push(this->matches3,(void *) match);
	  }
	}
      }
    }
  }

  return foundpairp;
}


/* A dangling match is one that has not been paired up with a match
   from the other end of the cDNA sequence.  A significant fraction of
   dangling matches may indicate the need to continue finding more
   matches from the other end.  Conversely, the absence of dangling
   matches on both ends indicates that no other possibilities exist in
   the genome. */

#ifdef DEBUG
static int
count_dangling (List_T matches) {
  int ndangling = 0;
  Match_T match;
  List_T p;

  for (p = matches; p != NULL; p = p->rest) {
    match = (Match_T) p->first;
    if (Match_pairedp(match) == false) {
      ndangling++;
    }
  }
  return ndangling;
}
#endif

static double
dangling_pct (List_T matches) {
  int ndangling = 0, n = 0;
  Match_T match;
  List_T p;

  for (p = matches; p != NULL; p = p->rest) {
    match = (Match_T) p->first;
    if (Match_pairedp(match) == false) {
      ndangling++;
    }
    n++;
  }
  if (n == 0) {
    return 0.0;
  } else {
    return (double) ndangling/(double) n;
  }
}

static List_T
get_dangling (List_T matches) {
  List_T dangling = NULL, p;
  Match_T match;

  for (p = matches; p != NULL; p = p->rest) {
    match = (Match_T) p->first;
    if (Match_pairedp(match) == false) {
      dangling = List_push(dangling,match);
    }
  }
  return dangling;
}



List_T
Stage1_matchlist (T this,
#ifdef PMAP
		  Indexdb_T indexdb_fwd,
		  Indexdb_T indexdb_rev,
#else
		  Indexdb_T indexdb,
#endif
		  IIT_T chromosome_iit, Chrsubset_T chrsubset) {
  List_T clusterlist = NULL;
  int querystart, queryend;
  int nsamples, prevnskip, nskip;

  debug(printf("\nMatchbestlist requested\n"));

  querystart = 0;
  queryend = this->querylength-1;

  if (this->matchpairlist != NULL && 
      dangling_pct(this->matches5) <= MAX_DANGLING_PCT && 
      dangling_pct(this->matches3) <= MAX_DANGLING_PCT) {
    debug(printf("Dangling5 = %f, Dangling3 = %f\n",
		 dangling_pct(this->matches5),dangling_pct(this->matches3)));
    debug(printf("Returning matchpairlist\n"));
    return this->matchpairlist;

  } else if (this->querylength <= SHORTSEQLEN) {
    /* Sample exhaustively on short sequences */
    nskip = 0;
    debug1(printf("***Beginning terminal sampling.  nskip = %d\n",nskip));
#ifdef PMAP
    sample(this,nskip,2*this->maxentries,indexdb_fwd,indexdb_rev,chromosome_iit,chrsubset);
#else
    sample(this,nskip,2*this->maxentries,indexdb,chromosome_iit,chrsubset);
#endif
    clusterlist = Matchpair_find_clusters(this->matches5,this->matches3,this->interval+INDEX1PART,
					  this->maxintronlen,1,SIZEBOUND,MOVING_THRESHOLD);
    debug1(printf("***Ending terminal sampling.  Got %d clusters\n",List_length(clusterlist)));

  } else {
    nsamples = 2;			/* Should get ends and middle */
    prevnskip = nskip = (queryend - querystart - nsamples*(this->interval+INDEX1PART))/nsamples;
    if (nskip < 1) {
      prevnskip = nskip = 0;
    }

    while (nskip > 0 && clusterlist == NULL) {
      debug1(printf("***Beginning sampling.  nsamples = %d, nskip = %d\n",nsamples,nskip));
#ifdef PMAP
      sample(this,nskip,2*this->maxentries,indexdb_fwd,indexdb_rev,chromosome_iit,chrsubset);
#else
      sample(this,nskip,2*this->maxentries,indexdb,chromosome_iit,chrsubset);
#endif
      clusterlist = Matchpair_find_clusters(this->matches5,this->matches3,this->interval+INDEX1PART,
					    this->maxintronlen,2,SIZEBOUND,BY_CANDIDATES);
      debug1(printf("***Ending sampling.  Got %d clusters\n",List_length(clusterlist)));
      
      nsamples = 2*nsamples + 1;
      prevnskip = nskip;
      nskip = (queryend - querystart - nsamples*(this->interval+INDEX1PART))/nsamples;
      if (nskip < 0) {
	nskip = 0;
      }
    }

    if (nskip == 0 && clusterlist == NULL) {
      nskip = prevnskip;
      debug1(printf("***Redoing sampling.  nsamples = %d, nskip = %d\n",nsamples,nskip));
      clusterlist = Matchpair_find_clusters(this->matches5,this->matches3,this->interval+INDEX1PART,
					    this->maxintronlen,2,SIZEBOUND,MOVING_THRESHOLD);
      debug1(printf("***Ending sampling.  Got %d clusters\n",List_length(clusterlist)));

#if 0
      debug1(printf("***Beginning terminal sampling.  nsamples = %d, nskip = %d\n",nsamples,nskip));
#ifdef PMAP
      sample(this,nskip,2*this->maxentries,indexdb_fwd,indexdb_rev,chromosome_iit,chrsubset);
#else
      sample(this,nskip,2*this->maxentries,indexdb,chromosome_iit,chrsubset);
#endif
      clusterlist = Matchpair_find_clusters(this->matches5,this->matches3,this->interval+INDEX1PART,
					    this->maxintronlen,1,SIZEBOUND,MOVING_THRESHOLD);
      debug1(printf("***Ending terminal sampling.  Got %d clusters\n",List_length(clusterlist)));
#endif
    }
  }

  this->matchbestlist = clusterlist;
  return this->matchbestlist;
}


T
Stage1_compute (Sequence_T queryuc, 
#ifdef PMAP
		Indexdb_T indexdb_fwd,
		Indexdb_T indexdb_rev,
#else
		Indexdb_T indexdb,
#endif
		IIT_T chromosome_iit, Chrsubset_T chrsubset, int maxintronlen_bound, int stuttercycles, 
		int stutterhits, bool crossspeciesp) {
  T this;
  int trimlength, stage1size, maxentries, maxintronlen, reader_overlap;
  double dangling5_pct, dangling3_pct;
  List_T dangling5, dangling3;
  bool foundpairp = false;

  trimlength = Sequence_trimlength(queryuc);
#ifdef PMAP
  stage1size = INDEX1PART;	/* 7-mer */
  maxentries = 100;
  maxintronlen = 20 + EXPANSION*(trimlength*3 - SHORTSEQLEN);
  if (trimlength*3 <= VERYSHORTSEQLEN) {
    reader_overlap = INDEX1PART;
  } else {
    reader_overlap = 0;
  }
#else
  if (trimlength <= VERYSHORTSEQLEN) {
    stage1size = 12;
    maxentries = 250;
    maxintronlen = MININTRONLEN; /* Essentially, we don't want to find any introns */
    reader_overlap = INDEX1PART;
  } else if (trimlength <= SHORTSEQLEN) {
    stage1size = 18;
    maxentries = 100;
    maxintronlen = MININTRONLEN; /* Essentially, we don't want to find any introns */
    reader_overlap = 0;
  } else if (crossspeciesp == true) {
    stage1size = 18;
    maxentries = 250;
    maxintronlen = 20 + EXPANSION*(trimlength - SHORTSEQLEN);
    reader_overlap = 0;
  } else {
    stage1size = 24;
    maxentries = 40;
    maxintronlen = 20 + EXPANSION*(trimlength - SHORTSEQLEN);
    reader_overlap = 0;
  }
#endif

  if (maxintronlen > maxintronlen_bound) {
    maxintronlen = maxintronlen_bound;
  } else if (maxintronlen < MININTRONLEN) {
    maxintronlen = MININTRONLEN;
  }

  this = Stage1_new(queryuc,maxintronlen,stage1size,maxentries,reader_overlap);
#ifdef PMAP
  foundpairp = find_first_pair(this,indexdb_fwd,indexdb_rev,chromosome_iit,chrsubset,maxintronlen,maxentries);
#else
  foundpairp = find_first_pair(this,indexdb,chromosome_iit,chrsubset,maxintronlen,maxentries);
#endif

  if (foundpairp == false) {
    debug(printf("*** No pair found ***\n"));
    this->matchpairlist = NULL;
  } else {
    debug(printf("*** Found pair ***\n"));
#ifdef PMAP
    stutter(this,indexdb_fwd,indexdb_rev,chromosome_iit,chrsubset,maxintronlen,
	    stuttercycles,stutterhits,stage1size,maxentries);
#else
    stutter(this,indexdb,chromosome_iit,chrsubset,maxintronlen,
	    stuttercycles,stutterhits,stage1size,maxentries);
#endif

    dangling5_pct = dangling_pct(this->matches5);
    dangling3_pct = dangling_pct(this->matches3);
    debug(printf("Dangling on 5' end: %d/%d\n",
		 count_dangling(this->matches5),List_length(this->matches5)));
    debug(printf("Dangling on 3' end: %d/%d\n",
		 count_dangling(this->matches3),List_length(this->matches3)));

    if (dangling5_pct > MAX_DANGLING_PCT) {
      dangling5 = get_dangling(this->matches5);
#ifdef PMAP
      fill_in_3(this,dangling5,indexdb_fwd,indexdb_rev,chromosome_iit,chrsubset,maxintronlen,maxentries);
#else
      fill_in_3(this,dangling5,indexdb,chromosome_iit,chrsubset,maxintronlen,maxentries);
#endif
      List_free(&dangling5);
    }
    if (dangling3_pct > MAX_DANGLING_PCT) {
      dangling3 = get_dangling(this->matches3);
#ifdef PMAP
      fill_in_5(this,dangling3,indexdb_fwd,indexdb_rev,chromosome_iit,chrsubset,maxintronlen,maxentries);
#else
      fill_in_5(this,dangling3,indexdb,chromosome_iit,chrsubset,maxintronlen,maxentries);
#endif
      List_free(&dangling3);
    }

    this->matchpairlist = Matchpair_filter_unique(this->matchpairlist);
  }
  debug(printf("*** Returning from Stage1_compute ***\n\n"));

  return this;
}

