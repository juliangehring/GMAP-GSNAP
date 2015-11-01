static char rcsid[] = "$Id: samprint.c 112423 2013-10-23 22:18:52Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "samprint.h"
#include "samflags.h"
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include "mem.h"
#include "complement.h"
#include "stage3hr.h"
#include "assert.h"


#define SANGER_ILLUMINA_DIFF 31

/* BAM appears to truncate the H information on the ends of a cigar */
/* Also, this provides the information needed for getting term information */


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


/* compute_cigar */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif


/* print_md_string */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif


/* overlap */
#ifdef DEBUG3
#define debug3(x) x
#else
#define debug3(x)
#endif


static bool quiet_if_excessive_p;
static int maxpaths_report;
static bool sam_multiple_primaries_p;
static bool force_xs_direction_p;
static bool md_lowercase_variant_p;
static IIT_T snps_iit;

void
SAM_setup (bool quiet_if_excessive_p_in, int maxpaths_report_in, bool sam_multiple_primaries_p_in,
	   bool force_xs_direction_p_in, bool md_lowercase_variant_p_in, IIT_T snps_iit_in) {
  quiet_if_excessive_p = quiet_if_excessive_p_in;
  maxpaths_report = maxpaths_report_in;
  sam_multiple_primaries_p = sam_multiple_primaries_p_in;
  force_xs_direction_p = force_xs_direction_p_in;
  md_lowercase_variant_p = md_lowercase_variant_p_in;
  snps_iit = snps_iit_in;
  return;
}


unsigned int
SAM_compute_flag (bool plusp, Stage3end_T mate, Resulttype_T resulttype,
		  bool first_read_p, int pathnum, int npaths, int npaths_mate,
		  int absmq_score, int first_absmq, bool invertp, bool invert_mate_p) {
  unsigned int flag = 0U;

  debug(printf("Resulttype: %s\n",Resulttype_string(resulttype)));

  if (npaths == 0) {
    debug(printf("npaths = 0, so QUERY_UNMAPPED %d\n",QUERY_UNMAPPED));
    flag |= QUERY_UNMAPPED;
  } else if (plusp == invertp) {
    debug(printf("plusp %d and invertp %d, so QUERY_MINUSP %d\n",
		 plusp,invertp,QUERY_MINUSP));
    flag |= QUERY_MINUSP;
  }

  if (resulttype == SINGLEEND_NOMAPPING || resulttype == SINGLEEND_UNIQ || resulttype == SINGLEEND_TRANSLOC || resulttype == SINGLEEND_MULT) {
    /* No first or second read or mate */
  } else {
    debug(printf("PAIRED_READ %d\n",PAIRED_READ));
    flag |= PAIRED_READ;
    if (first_read_p == true) {
      debug(printf("FIRST_READ %d\n",FIRST_READ_P));
      flag |= FIRST_READ_P;
    } else {
      debug(printf("SECOND_READ %d\n",SECOND_READ_P));
      flag |= SECOND_READ_P;
    }
    if (npaths_mate == 0) {
      debug(printf("MATE_UNMAPPED %d\n",MATE_UNMAPPED));
      flag |= MATE_UNMAPPED;

    } else if (quiet_if_excessive_p && npaths_mate > maxpaths_report) {
      debug(printf("MATE_UNMAPPED %d\n",MATE_UNMAPPED));
      flag |= MATE_UNMAPPED;

    } else if (mate == NULL) {
      /* Unpaired; no mate.  Not clear if should be MATE_UNMAPPED. */

    } else if (npaths == 0) {
      /* Need to check npaths == 0 in case clipping of overlaps results in a nomapping */
      if (Stage3end_plusp(mate) == invert_mate_p) {
	debug(printf("MATE_MINUSP %d\n",MATE_MINUSP));
	flag |= MATE_MINUSP;
      }

    } else if (resulttype == CONCORDANT_UNIQ || resulttype == CONCORDANT_TRANSLOC || resulttype == CONCORDANT_MULT) {
      /* Can distinguish concordant mappings by presence of insert length */
      debug(printf("PAIRED_MAPPING %d\n",PAIRED_MAPPING));
      flag |= PAIRED_MAPPING;
      if (0 && Stage3end_chrnum(mate) == 0) {
	/* Splice without a direction.  But want the effective plusp anyway. */

      } else if (Stage3end_plusp(mate) == invert_mate_p) {
	debug(printf("MATE_MINUSP %d\n",MATE_MINUSP));
	flag |= MATE_MINUSP;
      }

    } else if (resulttype == PAIRED_UNIQ || resulttype == PAIRED_MULT) {
      /* Note: We are counting PAIRED_UNIQ and PAIRED_MULT as "paired" mappings.
	 However, we are no longer counting UNPAIRED_UNIQ as a "paired" mapping. */
      debug(printf("PAIRED_MAPPING %d\n",PAIRED_MAPPING));
      flag |= PAIRED_MAPPING;
      if (0 && Stage3end_chrnum(mate) == 0) {
	/* Splice without a direction.  But want the effective plusp anyway. */

      } else if (Stage3end_plusp(mate) == invert_mate_p) {
	debug(printf("MATE_MINUSP %d\n",MATE_MINUSP));
	flag |= MATE_MINUSP;
      }

    } else {
      if (0 && Stage3end_chrnum(mate) == 0) {
	/* Splice without a direction.  But want the effective plusp anyway. */

      } else if (Stage3end_plusp(mate) == invert_mate_p) {
	debug(printf("MATE_MINUSP %d\n",MATE_MINUSP));
	flag |= MATE_MINUSP;
      }
    }
  }

  if (pathnum > 1) {
    if (sam_multiple_primaries_p == false) {
      debug(printf("NOT_PRIMARY %d\n",NOT_PRIMARY));
      flag |= NOT_PRIMARY;
    } else if (absmq_score != first_absmq) {
      debug(printf("NOT_PRIMARY %d\n",NOT_PRIMARY));
      flag |= NOT_PRIMARY;
    } else {
      /* Just as good as first alignment, so don't mark as secondary */
    }
  }

  return flag;
}


Chrpos_T
SAM_compute_chrpos (int *hardclip_low, int *hardclip_high,
		    int clipdir, int hardclip5, int hardclip3, bool firstp,
		    Stage3end_T this, Substring_T substring_low, int querylength) {
  Chrpos_T chrpos;
  Substring_T substring;
  int querystart, queryend;
  bool plusp;

  if (this == NULL) {
    return 0U;

  } else if (Stage3end_hittype(this) == GMAP) {
    chrpos = Pair_genomicpos_low(clipdir,hardclip5,hardclip3,
				 Stage3end_pairarray(this),Stage3end_npairs(this),
				 querylength,/*watsonp*/Stage3end_plusp(this),firstp);

  } else {
    if (substring_low != NULL) {
      plusp = Substring_plusp(substring_low);
    } else {
      plusp = Stage3end_plusp(this);
    }

    if (plusp == true) {
      /* Add 1 to report in 1-based coordinates */
      if (*hardclip_low > Substring_querystart(substring_low)) {
	substring = (Substring_T) NULL;
	while (*hardclip_low < querylength && substring == NULL) {
	  substring = Stage3end_substring_containing(this,*hardclip_low);
	  (*hardclip_low)++;
	}
	(*hardclip_low)--;

	if (substring == NULL) {
	  return 0U;
	}

	chrpos = Substring_alignstart_trim(substring) - Substring_chroffset(substring) + 1U;
	querystart = Substring_querystart(substring);
	queryend = Substring_queryend(Stage3end_substring_high(this));

#if 0
	fprintf(stderr,"case 1, hardclip_low %d, hardclip_high %d, querystart %d, queryend %d\n",
		*hardclip_low,*hardclip_high,querystart,queryend);
#endif

	if (*hardclip_low >= queryend || querylength - *hardclip_high < querystart) {
	  /* fprintf(stderr,"Returning 0\n"); */
	  return 0U;
	}
#if 0
	fprintf(stderr,"Adding to chrpos %u: hardclip_low %d - querystart %d\n",
		chrpos,*hardclip_low,querystart);
#endif

	chrpos += (*hardclip_low) - querystart;

      } else {
	querystart = Substring_querystart(substring_low);
	queryend = Substring_queryend(Stage3end_substring_high(this));

#if 0
	fprintf(stderr,"case 2, hardclip_low %d, hardclip_high %d, querystart %d, queryend %d\n",
		*hardclip_low,*hardclip_high,querystart,queryend);
#endif

	if (*hardclip_low >= queryend || querylength - *hardclip_high < querystart) {
	  /* fprintf(stderr,"Returning 0\n"); */
	  return 0U;
	}

	chrpos = Substring_alignstart_trim(substring_low) - Substring_chroffset(substring_low) + 1U;
      }

    } else {
      /* Add 1 to report in 1-based coordinates */
      if (querylength - *hardclip_low < Substring_queryend(substring_low)) {
	substring = (Substring_T) NULL;
	while (*hardclip_low < querylength && substring == NULL) {
	  substring = Stage3end_substring_containing(this,querylength - 1 - (*hardclip_low));
	  (*hardclip_low)++;
	}
	(*hardclip_low)--;

	if (substring == NULL) {
	  return 0U;
	}

	chrpos = Substring_alignend_trim(substring) - Substring_chroffset(substring) + 1U;
	queryend = Substring_queryend(substring);
	querystart = Substring_querystart(Stage3end_substring_high(this));

#if 0
	fprintf(stderr,"case 3, hardclip_low %d, hardclip_high %d, querystart %d, queryend %d\n",
		*hardclip_low,*hardclip_high,querystart,queryend);
#endif

	if (*hardclip_high >= queryend || querylength - *hardclip_low < querystart) {
	  /* fprintf(stderr,"Returning 0\n"); */
	  return 0U;
	}

#if 0
	fprintf(stderr,"Adding to chrpos %u: queryend %d - (querylength %d - hardclip_high %d)\n",
		chrpos,queryend,querylength,*hardclip_high);
#endif

	chrpos += queryend - (querylength - (*hardclip_low));
      } else {
	querystart = Substring_querystart(Stage3end_substring_high(this));
	queryend = Substring_queryend(substring_low);

#if 0
	fprintf(stderr,"case 4, hardclip_low %d, hardclip_high %d, querystart %d, queryend %d\n",
		*hardclip_low,*hardclip_high,querystart,queryend);
#endif

	if (*hardclip_high >= queryend || querylength - *hardclip_low < querystart) {
	  /* fprintf(stderr,"Returning 0\n"); */
	  return 0U;
	}

	chrpos = Substring_alignend_trim(substring_low) - Substring_chroffset(substring_low) + 1U;
      }
    }
  }
    
  return chrpos;
}


static void
print_chromosomal_pos (FILE *fp, Chrnum_T chrnum, Chrpos_T chrpos,
		       Univ_IIT_T chromosome_iit) {
  bool allocp;
  char *chr;

#if 0
  if (chrpos == 0U) {
    /* No mapping */
    fprintf(fp,"\t*\t0");
    return;
  }
#endif

  if (chrnum == 0) {
    /* Interchromosomal splice */
    fprintf(stderr,"Trying to print interchrosomal splice in one line\n");
    abort();

  } else {
    chr = Univ_IIT_label(chromosome_iit,chrnum,&allocp);
    fprintf(fp,"\t%s",chr);
    if (allocp == true) {
      FREE(chr);
    }

    /* chrpos already in 1-based coordinates */
    fprintf(fp,"\t%u",chrpos /*+1U*/);
    return;
  }
}

static void
print_mate_chromosomal_pos (FILE *fp, Chrnum_T mate_chrnum, Chrnum_T mate_effective_chrnum,
			    Chrpos_T mate_chrpos, Chrnum_T anchor_chrnum, Chrpos_T anchor_chrpos,
			    Univ_IIT_T chromosome_iit) {
  bool allocp;
  char *chr;

  if (mate_chrpos == 0U) {
    /* No mapping */
    fprintf(fp,"\t*\t0");
    return;

  } else if (mate_chrnum == 0) {
    /* Interchromosomal splice.  Choose effective chrnum. */
    if (anchor_chrpos > 0U && anchor_chrnum > 0 && mate_effective_chrnum == anchor_chrnum) {
      fprintf(fp,"\t=");
    } else {
      chr = Univ_IIT_label(chromosome_iit,mate_effective_chrnum,&allocp);
      fprintf(fp,"\t%s",chr);
      if (allocp == true) {
	FREE(chr);
      }
    }
    
    /* chrpos already in 1-based coordinates */
    fprintf(fp,"\t%u",mate_chrpos /*+1U*/);
    return;

  } else {
    if (anchor_chrpos > 0U && anchor_chrnum > 0 && mate_chrnum == anchor_chrnum) {
      fprintf(fp,"\t=");
    } else {
      chr = Univ_IIT_label(chromosome_iit,mate_chrnum,&allocp);
      fprintf(fp,"\t%s",chr);
      if (allocp == true) {
	FREE(chr);
      }
    }
    
    /* chrpos already in 1-based coordinates */
    fprintf(fp,"\t%u",mate_chrpos /*+1U*/);
    return;
  }
}





static char complCode[128] = COMPLEMENT_LC;

static void
make_complement_buffered (char *complement, char *sequence, unsigned int length) {
  int i, j;

  /* complement = (char *) CALLOC(length+1,sizeof(char)); */
  for (i = length-1, j = 0; i >= 0; i--, j++) {
    complement[j] = complCode[(int) sequence[i]];
  }
  complement[length] = '\0';
  return;
}



void
SAM_print_nomapping (FILE *fp, char *abbrev, Shortread_T queryseq, Stage3end_T mate, char *acc1, char *acc2,
		     Univ_IIT_T chromosome_iit, Resulttype_T resulttype, bool first_read_p,
		     int npaths_mate, Chrpos_T mate_chrpos, int quality_shift,
		     char *sam_read_group_id, bool invertp, bool invert_mate_p) {
  unsigned int flag;


  /* 1. QNAME */
  if (acc2 == NULL) {
    fprintf(fp,"%s",acc1);
  } else {
    fprintf(fp,"%s,%s",acc1,acc2);
  }
  
  /* 2. FLAG */
  flag = SAM_compute_flag(/*plusp (NA)*/true,mate,resulttype,first_read_p,
			  /*pathnum*/0,/*npaths*/0,npaths_mate,
			  /*absmq_score*/0,/*first_absmq*/0,invertp,invert_mate_p);
  fprintf(fp,"\t%u",flag);

  /* 3. RNAME: chr */
  fprintf(fp,"\t*");

  /* 4. POS: chrpos */
  fprintf(fp,"\t0");

  /* 5. MAPQ: Mapping quality */
  /* Picard says MAPQ should be 0 for an unmapped read */
  fprintf(fp,"\t0");

  /* 6. CIGAR */
  fprintf(fp,"\t*");

  /* 7. MRNM: Mate chr */
  /* 8. MPOS: Mate chrpos */
  print_mate_chromosomal_pos(fp,Stage3end_chrnum(mate),Stage3end_effective_chrnum(mate),mate_chrpos,
			     /*anchor_chrnum*/0,/*anchor_chrpos*/0U,chromosome_iit);


  /* 9. ISIZE: Insert size */
  fprintf(fp,"\t0");

  /* 10. SEQ: queryseq and 11. QUAL: quality scores */
  /* Since there is no mapping, we print the original query sequence. */
  fprintf(fp,"\t");
  if (invertp == false) {
    Shortread_print_chopped(fp,queryseq,/*hardclip_low*/0,/*hardclip_high*/0);
    fprintf(fp,"\t");
    Shortread_print_quality(fp,queryseq,/*hardclip_low*/0,/*hardclip_high*/0,
			    quality_shift,/*show_chopped_p*/false);
  } else {
    Shortread_print_chopped_revcomp(fp,queryseq,/*hardclip_low*/0,/*hardclip_high*/0);
    fprintf(fp,"\t");
    Shortread_print_quality_revcomp(fp,queryseq,/*hardclip_low*/0,/*hardclip_high*/0,
				    quality_shift,/*show_chopped_p*/false);
  }

  /* 12. TAGS: RG */
  if (sam_read_group_id != NULL) {
    fprintf(fp,"\t");
    fprintf(fp,"RG:Z:%s",sam_read_group_id);
  }
  
  /* 12. TAGS: XB */
  Shortread_print_barcode(fp,queryseq);

  /* 12. TAGS: XP */
  Shortread_print_chop(fp,queryseq,invertp);

  /* 12. TAGS: XO */
  fprintf(fp,"\t");
  fprintf(fp,"XO:Z:%s",abbrev);

  fprintf(fp,"\n");
  return;
}


/* Derived from print_tokens_gff3 */
static void
print_tokens_sam (FILE *fp, List_T tokens) {
  List_T p;
  char *token;
  
  for (p = tokens; p != NULL; p = List_next(p)) {
    token = (char *) List_head(p);
    fprintf(fp,"%s",token);
    FREE(token);
  }

  return;
}

static List_T
push_token (List_T tokens, char *token) {
  char *copy;

  copy = (char *) CALLOC(strlen(token)+1,sizeof(char));
  strcpy(copy,token);
  return List_push(tokens,(void *) copy);
}


#if 0
/* Currently used for insertions and deletions */
static List_T
compute_cigar_old (List_T tokens, char type, int stringlength, int querypos, int querylength,
		   int hardclip_low, int hardclip_high, bool plusp, bool firstp, bool lastp) {
  char token[10];
  
  debug1(printf("\nEntering compute_cigar with type %c, stringlength %d, querypos %d, querylength %d, hardclip_low %d, hardclip_high %d, plusp %d\n",
		type,stringlength,querypos,querylength,hardclip_low,hardclip_high,plusp));

  if (firstp == true) {
    debug1(printf("firstp is true\n"));
    if (plusp == true) {
      if (hardclip_low > 0) {
	sprintf(token,"%dH",hardclip_low);
	debug1(printf("Pushing token %s\n",token));
	tokens = push_token(tokens,token);
      }
      if (querypos > hardclip_low) {
	sprintf(token,"%dS",querypos - hardclip_low);
	debug1(printf("Pushing token %s\n",token));
	tokens = push_token(tokens,token);
      }
    } else {
      if (hardclip_high > 0) {
	sprintf(token,"%dH",hardclip_high);
	debug1(printf("Pushing token %s\n",token));
	tokens = push_token(tokens,token);
      }
      if (querypos < querylength - hardclip_high) {
	sprintf(token,"%dS",querypos - hardclip_high);
	debug1(printf("Pushing token %s\n",token));
	tokens = push_token(tokens,token);
      }
    }
  }

  if (type == 'D' || type == 'N') {
    if (querypos < hardclip_low || querypos >= querylength - hardclip_high) {
      stringlength = 0;
    }

  } else if (plusp == true) {
    debug1(printf("Comparing querypos %d..%d against %d..%d\n",
		  querypos,querypos + stringlength,hardclip_low,querylength - hardclip_high));
    if (/* querypos < hardclip_low && */querypos + stringlength < hardclip_low) {
      /* Print nothing */
      stringlength = 0;
      debug1(printf("Case 1: stringlength 0\n"));
    } else if (querypos < hardclip_low) {
      if (querypos + stringlength < querylength - hardclip_high) {
	/* Print part after hardclip_low */
	stringlength = (querypos + stringlength) - hardclip_low;
	debug1(printf("Case 2: stringlength %d\n",stringlength));
      } else {
	/* Print part between hardclip_low and hardclip_high */
	stringlength = (querylength - hardclip_high) - hardclip_low;
	debug1(printf("Case 3: stringlength %d\n",stringlength));
      }
    } else if (querypos < querylength - hardclip_high) {
      if (querypos + stringlength >= querylength - hardclip_high) {
	/* Print up to hardclip_high */
	stringlength = (querylength - hardclip_high) - querypos;
	debug1(printf("Case 4: stringlength %d\n",stringlength));
      } else {
	/* Print full stringlength */
	debug1(printf("Case 5: stringlength %d\n",stringlength));
      }
    } else {
      /* Print nothing */
      stringlength = 0;
      debug1(printf("Case 6: stringlength 0\n"));
    }

  } else {
    debug1(printf("Comparing querypos %d..%d against %d..%d\n",
		  querypos,querypos - stringlength,hardclip_low,querylength - hardclip_high));
    if (/* querypos >= querylength - hardclip_high && */ querypos - stringlength >= querylength - hardclip_high) {
      /* Print nothing */
      stringlength = 0;
      debug1(printf("Case 1: stringlength 0\n"));
    } else if (querypos >= querylength - hardclip_high) {
      if (querypos - stringlength >= hardclip_low) {
	/* Print part after hardclip_high */
	stringlength = (querylength - hardclip_high) - (querypos - stringlength);
	debug1(printf("Case 2: stringlength %d\n",stringlength));
      } else {
	/* Print part between hardclip_low and hardclip_high */
	stringlength = (querylength - hardclip_high) - hardclip_low;
	debug1(printf("Case 3: stringlength %d\n",stringlength));
      }
    } else if (querypos >= hardclip_low) {
      if (querypos - stringlength < hardclip_low) {
	/* Print up to hardclip_low */
	stringlength = querypos - hardclip_low;
	debug1(printf("Case 4: stringlength %d\n",stringlength));
      } else {
	/* Print full stringlength */
	debug1(printf("Case 5: stringlength %d\n",stringlength));
      }
    } else {
      /* Print nothing */
      stringlength = 0;
      debug1(printf("Case 5: stringlength 0\n"));
    }
  }

  if (stringlength > 0) {
    sprintf(token,"%d%c",stringlength,type);
    debug1(printf("Pushing token %s\n",token));
    tokens = push_token(tokens,token);
  }

  if (lastp == true) {
    debug1(printf("lastp is true\n"));
    if (plusp == true) {
      querypos += stringlength;
      if (querypos < querylength - 1 - hardclip_high) {
	sprintf(token,"%dS",querylength - 1 - hardclip_high - querypos);
	debug1(printf("Pushing token %s\n",token));
	tokens = push_token(tokens,token);
      }
      if (hardclip_high > 0) {
	sprintf(token,"%dH",hardclip_high);
	debug1(printf("Pushing token %s\n",token));
	tokens = push_token(tokens,token);
      }
    } else {
      querypos -= stringlength;
      if (querypos > hardclip_low) {
	sprintf(token,"%dS",hardclip_low - querypos);
	debug1(printf("Pushing token %s\n",token));
	tokens = push_token(tokens,token);
      }
      if (hardclip_low > 0) {
	sprintf(token,"%dH",hardclip_low);
	debug1(printf("Pushing token %s\n",token));
	tokens = push_token(tokens,token);
      }
    }
  }

  return tokens;
}
#endif


/* Currently used for insertions and deletions */
static List_T
compute_cigar (List_T tokens, char type, int stringlength, int querypos, int querylength,
	       int hardclip_low, int hardclip_high, bool plusp, int lastp) {
  int matchlength = 0;
  int startpos, endpos;
  int cliplength = 0;
  char token[10];
  
  if (plusp == true) {
    debug1(printf("\nEntering compute_cigar with type %c, stringlength %d, querypos %d, querylength %d, hardclip_low %d, hardclip_high %d, plus\n",
		  type,stringlength,querypos,querylength,hardclip_low,hardclip_high));
    if (hardclip_low > querypos) { /* > not >= */
      startpos = hardclip_low;
      cliplength = hardclip_low;
    } else {
      startpos = querypos;
    }

    if (querylength - hardclip_high < querypos + stringlength) {
      endpos = querylength - hardclip_high;
      debug1(printf("  endpos %d = querylength %d - hardclip_high %d\n",endpos,querylength,hardclip_high));
    } else {
      endpos = querypos + stringlength;
      debug1(printf("  endpos %d = querypos %d + stringlength %d\n",endpos,querypos,stringlength));
    }

    debug1(printf("  new startpos %d, endpos %d, cliplength %d\n",startpos,endpos,cliplength));

    if (endpos >= startpos) {
      if (cliplength > 0) {
	debug1(printf("  Pushing initial %dH\n",cliplength));
	sprintf(token,"%dH",cliplength);
	debug1(printf("Pushing token %s\n",token));
	tokens = push_token(tokens,token);
      }
      matchlength = endpos - startpos;
      if (matchlength > 0) {
	debug1(printf("  Pushing %d%c\n",matchlength,type));
	sprintf(token,"%d%c",matchlength,type);
	debug1(printf("Pushing token %s\n",token));
	tokens = push_token(tokens,token);
      }
    }


    if (lastp == true) {
      /* cliplength = querypos + stringlength - endpos; */
      cliplength = querylength - endpos;
      if (cliplength > 0) {
	debug1(printf("  Pushing final %dH\n",cliplength));
	sprintf(token,"%dH",cliplength);
	debug1(printf("Pushing token %s\n",token));
	tokens = push_token(tokens,token);
      }
    }

  } else {
    debug1(printf("\nEntering compute_cigar with type %c, stringlength %d, querypos %d, querylength %d, hardclip_low %d, hardclip_high %d, minus\n",
		  type,stringlength,querypos,querylength,hardclip_low,hardclip_high));

    if (querylength - hardclip_low < querypos) {
      startpos = querylength - hardclip_low;
      cliplength = hardclip_low;
    } else {
      startpos = querypos;
    }

    if (hardclip_high >= querypos - stringlength) {
      endpos = hardclip_high;
      debug1(printf("  endpos %d = hardclip_high %d\n",endpos,hardclip_high));
    } else {
      endpos = querypos - stringlength;
      debug1(printf("  endpos %d = querypos %d - stringlength %d\n",endpos,querypos,stringlength));
    }

    debug1(printf("  new startpos %d, endpos %d, cliplength %d\n",startpos,endpos,cliplength));

    if (endpos <= startpos) {
      if (cliplength > 0) {
	debug1(printf("  Pushing initial %dH\n",cliplength));
	sprintf(token,"%dH",cliplength);
	debug1(printf("Pushing token %s\n",token));
	tokens = push_token(tokens,token);
      }
      matchlength = startpos - endpos;
      if (matchlength > 0) {
	debug1(printf("  Pushing %d%c\n",matchlength,type));
	sprintf(token,"%d%c",matchlength,type);
	debug1(printf("Pushing token %s\n",token));
	tokens = push_token(tokens,token);
      }
    }


    if (lastp == true) {
      cliplength = endpos;
      if (cliplength > 0) {
	debug1(printf("  Pushing final %dH\n",cliplength));
	sprintf(token,"%dH",cliplength);
	debug1(printf("Pushing token %s\n",token));
	tokens = push_token(tokens,token);
      }
    }
  }

  return tokens;
}


static int
print_cigar (FILE *fp, char type, int stringlength, int querypos, int querylength,
	     int hardclip_low, int hardclip_high, bool plusp, int lastp) {
  int matchlength = 0;
  int startpos, endpos;
  int cliplength = 0;
  
  if (plusp == true) {
    debug1(printf("\nEntering print_cigar with type %c, stringlength %d, querypos %d, querylength %d, hardclip_low %d, hardclip_high %d, plus\n",
		  type,stringlength,querypos,querylength,hardclip_low,hardclip_high));
    if (hardclip_low > querypos) { /* > not >= */
      startpos = hardclip_low;
      cliplength = hardclip_low;
    } else {
      startpos = querypos;
    }

    if (querylength - hardclip_high < querypos + stringlength) {
      endpos = querylength - hardclip_high;
      debug1(printf("  endpos %d = querylength %d - hardclip_high %d\n",endpos,querylength,hardclip_high));
    } else {
      endpos = querypos + stringlength;
      debug1(printf("  endpos %d = querypos %d + stringlength %d\n",endpos,querypos,stringlength));
    }

    debug1(printf("  new startpos %d, endpos %d, cliplength %d\n",startpos,endpos,cliplength));

    if (endpos >= startpos) {
      if (cliplength > 0) {
	debug1(printf("  Pushing initial %dH\n",cliplength));
	fprintf(fp,"%dH",cliplength);
      }
      matchlength = endpos - startpos;
      if (matchlength > 0) {
	debug1(printf("  Pushing %d%c\n",matchlength,type));
	fprintf(fp,"%d%c",matchlength,type);
      }
    }


    if (lastp == true) {
      /* cliplength = querypos + stringlength - endpos; */
      cliplength = querylength - endpos;
      if (cliplength > 0) {
	debug1(printf("  Pushing final %dH\n",cliplength));
	fprintf(fp,"%dH",cliplength);
      }
    }

  } else {
    debug1(printf("\nEntering print_cigar with type %c, stringlength %d, querypos %d, querylength %d, hardclip_low %d, hardclip_high %d, minus\n",
		  type,stringlength,querypos,querylength,hardclip_low,hardclip_high));

    if (querylength - hardclip_low < querypos) {
      startpos = querylength - hardclip_low;
      cliplength = hardclip_low;
    } else {
      startpos = querypos;
    }

    if (hardclip_high >= querypos - stringlength) {
      endpos = hardclip_high;
      debug1(printf("  endpos %d = hardclip_high %d\n",endpos,hardclip_high));
    } else {
      endpos = querypos - stringlength;
      debug1(printf("  endpos %d = querypos %d - stringlength %d\n",endpos,querypos,stringlength));
    }

    debug1(printf("  new startpos %d, endpos %d, cliplength %d\n",startpos,endpos,cliplength));

    if (endpos <= startpos) {
      if (cliplength > 0) {
	debug1(printf("  Pushing initial %dH\n",cliplength));
	fprintf(fp,"%dH",cliplength);
      }
      matchlength = startpos - endpos;
      if (matchlength > 0) {
	debug1(printf("  Pushing %d%c\n",matchlength,type));
	fprintf(fp,"%d%c",matchlength,type);
      }
    }


    if (lastp == true) {
      cliplength = endpos;
      if (cliplength > 0) {
	debug1(printf("  Pushing final %dH\n",cliplength));
	fprintf(fp,"%dH",cliplength);
      }
    }
  }

  return endpos;
}


static int
print_md_string (bool *printp, int *nmismatches_refdiff, int *nmismatches_bothdiff,
		 FILE *fp, int matchlength, char *genomicfwd_refdiff, char *genomicfwd_bothdiff,
		 int stringlength, int querypos, int querylength,
		 int hardclip_low, int hardclip_high, bool plusp, bool lastp) {
  int starti, endi, i;
  int local_nmismatches = 0;
  bool hardclip_end_p = false;

  if (plusp == true) {
    debug2(printf("\nEntering md_string with matchlength %d, querypos %d, querylength %d, hardclip_low %d, hardclip_high %d, plus: %s ref, %s both\n",
		  matchlength,querypos,querylength,hardclip_low,hardclip_high,genomicfwd_refdiff,genomicfwd_bothdiff));
    if (hardclip_low == 0) {
      starti = 0;
      hardclip_end_p = true;
    } else if (hardclip_low > querypos) {
      /* startpos = hardclip_low; */
      starti = hardclip_low - querypos;
      hardclip_end_p = true;
      debug2(printf("  Setting starti %d = hardclip_low %d - querypos %d\n",
		    starti,hardclip_low,querypos));
    } else {
      /* startpos = querypos; */
      starti = 0;
    }

    if (querylength - hardclip_high < querypos + stringlength) {
      /* endpos = querylength - hardclip_high; */
      endi = (querylength - hardclip_high) - querypos;
      debug2(printf("  Setting endi %d = (querylength %d - hardclip_high %d) - querypos %d\n",
		    endi,querylength,hardclip_high,querypos));
    } else {
      /* endpos = querypos + stringlength; */
      endi = stringlength;
    }

    debug2(printf("  Counting matches from %d to %d\n",starti,endi));

    if (genomicfwd_refdiff == NULL) {
      if (endi > starti) {
	matchlength += (endi - starti);
      }

    } else if (md_lowercase_variant_p == false) {
      for (i = starti; i < endi; i++) {
	if (isupper(genomicfwd_refdiff[i])) {
	  matchlength++;

	} else {
	  /* A true mismatch against both variants */
	  if (matchlength > 0 || hardclip_end_p == true) {
	    fprintf(fp,"%d",matchlength);
	    *printp = true;
	    hardclip_end_p = false;
	  }
	  fprintf(fp,"%c",toupper(genomicfwd_refdiff[i]));
	  *printp = true;
	  local_nmismatches += 1;
	  matchlength = 0;
	}
      }
      *nmismatches_refdiff += local_nmismatches;

    } else {
      for (i = starti; i < endi; i++) {
	if (isupper(genomicfwd_refdiff[i])) {
	  matchlength++;

	} else if (isupper(genomicfwd_bothdiff[i])) {
	  /* A mismatch against the reference only => alternate variant */
	  if (matchlength > 0 || hardclip_end_p == true) {
	    fprintf(fp,"%d",matchlength);
	    *printp = true;
	    hardclip_end_p = false;
	  }
	  fprintf(fp,"%c",genomicfwd_refdiff[i]); /* Leave as lower case */
	  *printp = true;
	  local_nmismatches += 1;
	  matchlength = 0;

	} else {
	  /* A true mismatch against both variants */
	  if (matchlength > 0 || hardclip_end_p == true) {
	    fprintf(fp,"%d",matchlength);
	    *printp = true;
	    hardclip_end_p = false;
	  }
	  fprintf(fp,"%c",toupper(genomicfwd_refdiff[i]));
	  *printp = true;
	  local_nmismatches += 1;
	  matchlength = 0;
	}
      }
      *nmismatches_refdiff += local_nmismatches;
    }

  } else {
    debug2(printf("\nEntering md_string with matchlength %d, querypos %d, querylength %d, hardclip_low %d, hardclip_high %d, minus: %s ref, %s both\n",
		  matchlength,querypos,querylength,hardclip_low,hardclip_high,genomicfwd_refdiff,genomicfwd_bothdiff));
    querypos = querylength - querypos - stringlength;
    debug2(printf("  Revising querypos to be %d\n",querypos));

    if (hardclip_low == 0) {
      starti = 0;
      hardclip_end_p = true;
    } else if (hardclip_low > querypos) {
      /* startpos = hardclip_low; */
      starti = hardclip_low - querypos;
      hardclip_end_p = true;
      debug2(printf("  Setting starti %d = hardclip_low %d - querypos %d\n",
		    starti,hardclip_low,querypos));
    } else {
      /* startpos = querypos; */
      starti = 0;
    }

    if (querylength - hardclip_high < querypos + stringlength) {
      /* endpos = querylength - hardclip_high; */
      endi = (querylength - hardclip_high) - querypos;
      debug2(printf("  Setting endi %d = (querylength %d - hardclip_high %d) - querypos %d\n",
		    endi,querylength,hardclip_high,querypos));
    } else {
      /* endpos = querypos + stringlength; */
      endi = stringlength;
    }

    debug2(printf("  Counting matches from %d to %d\n",starti,endi));

    if (genomicfwd_refdiff == NULL) {
      if (endi > starti) {
	matchlength += (endi - starti);
      }

    } else if (md_lowercase_variant_p == false) {
      for (i = starti; i < endi; i++) {
	if (isupper(genomicfwd_refdiff[i])) {
	  matchlength++;

	} else {
	  if (matchlength > 0 || hardclip_end_p == true) {
	    fprintf(fp,"%d",matchlength);
	    *printp = true;
	    hardclip_end_p = false;
	  }
	  fprintf(fp,"%c",toupper(genomicfwd_refdiff[i]));
	  *printp = true;
	  local_nmismatches += 1;
	  matchlength = 0;
	}
      }
      *nmismatches_refdiff += local_nmismatches;

    } else {
      for (i = starti; i < endi; i++) {
	if (isupper(genomicfwd_refdiff[i])) {
	  matchlength++;

	} else if (isupper(genomicfwd_bothdiff[i])) {
	  /* A mismatch against the reference only => alternate variant */
	  if (matchlength > 0 || hardclip_end_p == true) {
	    fprintf(fp,"%d",matchlength);
	    *printp = true;
	    hardclip_end_p = false;
	  }
	  fprintf(fp,"%c",genomicfwd_refdiff[i]); /* Leave as lower case */
	  *printp = true;
	  local_nmismatches += 1;
	  matchlength = 0;

	} else {
	  /* A true mismatch against both variants */
	  if (matchlength > 0 || hardclip_end_p == true) {
	    fprintf(fp,"%d",matchlength);
	    *printp = true;
	    hardclip_end_p = false;
	  }
	  fprintf(fp,"%c",toupper(genomicfwd_refdiff[i]));
	  *printp = true;
	  local_nmismatches += 1;
	  matchlength = 0;
	}
      }
      *nmismatches_refdiff += local_nmismatches;
    }
  }

  /* Update nmismatches_bothdiff */
  if (genomicfwd_bothdiff == NULL) {
    /* No change to nmismatches_bothdiff */
  } else if (genomicfwd_bothdiff == genomicfwd_refdiff) {
    *nmismatches_bothdiff += local_nmismatches;
  } else {
    for (i = starti; i < endi; i++) {
      if (!isupper(genomicfwd_bothdiff[i])) {
	*nmismatches_bothdiff += 1;
      }
    }
  }

  debug2(printf("  Ending with matchlength %d\n",matchlength));

  if (lastp == false) {
    return matchlength;
  } else if (matchlength > 0) {
    fprintf(fp,"%d",matchlength);
    *printp = true;
    return 0;
  } else {
    return 0;
  }
}


static void
print_single (FILE *fp, char *abbrev, Stage3end_T this, Stage3end_T mate,
	      char *acc1, char *acc2, int pathnum, int npaths,
	      int absmq_score, int first_absmq, int second_absmq, int mapq_score,
	      Univ_IIT_T chromosome_iit, Shortread_T queryseq, int pairedlength,
	      Chrpos_T chrpos, Chrpos_T mate_chrpos, int clipdir, int hardclip5, int hardclip3,
	      Resulttype_T resulttype, bool first_read_p,
	      int npaths_mate, int quality_shift,
	      char *sam_read_group_id, bool invertp, bool invert_mate_p, bool circularp) {
  unsigned int flag = 0U;
  Substring_T substring;
  int nmismatches_refdiff = 0, nmismatches_bothdiff = 0, querylength, substring_start, substring_length;
  int hardclip_low, hardclip_high;
  /* int mate_hardclip_low, mate_hardclip_high; */
  char *genomicfwd_refdiff, *genomicfwd_bothdiff, *genomicdir_refdiff, *genomicdir_bothdiff;
  bool plusp, printp;


  querylength = Shortread_fulllength(queryseq);
  plusp = Stage3end_plusp(this);
  substring = Stage3end_substring1(this);

  if (circularp == true) {
    /* clipdir should be +1 */
    if (1 || plusp == true) {
      hardclip_low = hardclip5;
      hardclip_high = hardclip3;
      /* mate_hardclip_low = 0; */
      /* mate_hardclip_high = 0; */
    } else {
      hardclip_low = hardclip3;
      hardclip_high = hardclip5;
      /* mate_hardclip_low = 0; */
      /* mate_hardclip_high = 0; */
    }

  } else {
    if (first_read_p == true) {
      if (clipdir >= 0) {
	hardclip_low = 0;
	hardclip_high = hardclip5;
      } else {
	hardclip_low = hardclip5;
	hardclip_high = 0;
      }
      /* mate_hardclip_low = hardclip3; */
      /* mate_hardclip_high = 0; */
      /* fprintf(stderr,"first read: hardclip_low = %d, hardclip_high = %d\n",hardclip_low,hardclip_high); */
    } else {
      if (clipdir >= 0) {
	hardclip_low = hardclip3;
	hardclip_high = 0;
      } else {
	hardclip_low = 0;
	hardclip_high = hardclip3;
      }
      /* mate_hardclip_low = 0; */
      /* mate_hardclip_high = hardclip5; */
      /* fprintf(stderr,"second read: hardclip_low = %d, hardclip_high = %d\n",hardclip_low,hardclip_high); */
    }
  }
  /* printf("clipdir is %d, hardclip_low %d, hardclip_high %d\n",clipdir,hardclip_low,hardclip_high); */


  /* 1. QNAME */
  if (acc2 == NULL) {
    fprintf(fp,"%s",acc1);
  } else {
    fprintf(fp,"%s,%s",acc1,acc2);
  }

  /* 2. FLAG */
  flag = SAM_compute_flag(Stage3end_plusp(this),mate,resulttype,first_read_p,
			  pathnum,npaths,npaths_mate,absmq_score,first_absmq,
			  invertp,invert_mate_p);
  fprintf(fp,"\t%u",flag);

  /* 3. RNAME: chr */
  /* 4. POS: chrpos */
  print_chromosomal_pos(fp,Stage3end_chrnum(this),chrpos,chromosome_iit);


  /* 5. MAPQ: Mapping quality */
  fprintf(fp,"\t%d",mapq_score);

  /* 6. CIGAR */
  fprintf(fp,"\t");

  if (plusp == true) {
    print_cigar(fp,/*type*/'S',Substring_querystart(substring),
		/*querypos*/0,querylength,hardclip_low,hardclip_high,
		/*plusp*/true,/*lastp*/false);
    print_cigar(fp,/*type*/'M',Substring_match_length(substring),
		/*querypos*/Substring_querystart(substring),querylength,
		hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);
    print_cigar(fp,/*type*/'S',querylength - Substring_queryend(substring),
		/*querypos*/Substring_queryend(substring),querylength,
		hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);

  } else {
    print_cigar(fp,/*type*/'S',querylength - Substring_queryend(substring),
		/*querypos*/querylength,querylength,hardclip_low,hardclip_high,
		/*plusp*/false,/*lastp*/false);
    print_cigar(fp,/*type*/'M',Substring_match_length(substring),
		/*querypos*/Substring_queryend(substring),querylength,
		hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
    print_cigar(fp,/*type*/'S',Substring_querystart(substring),
		/*querypos*/Substring_querystart(substring),querylength,
		hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
  }


  /* 7. MRNM: Mate chr */
  /* 8. MPOS: Mate chrpos */
  print_mate_chromosomal_pos(fp,Stage3end_chrnum(mate),Stage3end_effective_chrnum(mate),mate_chrpos,
			     Stage3end_chrnum(this),chrpos,chromosome_iit);


  /* 9. ISIZE: Insert size */
  if (resulttype == CONCORDANT_UNIQ || resulttype == CONCORDANT_TRANSLOC || resulttype == CONCORDANT_MULT) {
    if (plusp == invertp) {
      fprintf(fp,"\t%d",-pairedlength);
    } else {
      fprintf(fp,"\t%d",pairedlength);
    }
  } else if (mate_chrpos == 0) {
    fprintf(fp,"\t%d",pairedlength);
  } else if (chrpos < mate_chrpos) {
    fprintf(fp,"\t%d",pairedlength);
  } else if (chrpos > mate_chrpos) {
    fprintf(fp,"\t%d",-pairedlength);
  } else if (first_read_p == true) {
    fprintf(fp,"\t%d",pairedlength);
  } else {
    fprintf(fp,"\t%d",-pairedlength);
  }


  /* 10. SEQ: queryseq and 11. QUAL: quality scores */
  /* Queryseq has already been inverted, so just measure plusp relative to its current state */
  fprintf(fp,"\t");
  if (plusp == true) {
    Shortread_print_chopped(fp,queryseq,hardclip_low,hardclip_high);
    fprintf(fp,"\t");
    Shortread_print_quality(fp,queryseq,hardclip_low,hardclip_high,
			    quality_shift,/*show_chopped_p*/false);
  } else {
    Shortread_print_chopped_revcomp(fp,queryseq,hardclip_low,hardclip_high);
    fprintf(fp,"\t");
    Shortread_print_quality_revcomp(fp,queryseq,hardclip_low,hardclip_high,
				    quality_shift,/*show_chopped_p*/false);
  }

  /* 12. TAGS: RG */
  if (sam_read_group_id != NULL) {
    fprintf(fp,"\t");
    fprintf(fp,"RG:Z:%s",sam_read_group_id);
  }

  /* 12. TAGS: XB */
  Shortread_print_barcode(fp,queryseq);

  /* 12. TAGS: XP */
  Shortread_print_chop(fp,queryseq,invertp);

  /* 12. TAGS: MD */
  fprintf(fp,"\t");
  fprintf(fp,"MD:Z:");  
  printp = false;

  substring_start = Substring_querystart(substring);
  substring_length = Substring_match_length(substring);

  if ((genomicdir_bothdiff = Substring_genomic_bothdiff(substring)) == NULL) {
    if (plusp == true) {
      print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
		      fp,/*matchlength*/0,/*genomicfwd_refdiff*/NULL,/*genomicfwd_bothdiff*/NULL,
		      substring_length,/*querypos*/substring_start,querylength,
		      hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
    } else {
      print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
		      fp,/*matchlength*/0,/*genomicfwd_refdiff*/NULL,/*genomicfwd_bothdiff*/NULL,
		      substring_length,/*querypos*/substring_start,querylength,
		      hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
    }

  } else if (plusp == true) {
    genomicdir_refdiff = Substring_genomic_refdiff(substring);
    print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
		    fp,/*matchlength*/0,&(genomicdir_refdiff[substring_start]),&(genomicdir_bothdiff[substring_start]),
		    substring_length,/*querypos*/substring_start,querylength,
		    hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);

  } else if ((genomicdir_refdiff = Substring_genomic_refdiff(substring)) == genomicdir_bothdiff) {
    genomicfwd_refdiff = (char *) CALLOC(querylength+1,sizeof(char));
    make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring_start]),substring_length);
    print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
		    fp,/*matchlength*/0,genomicfwd_refdiff,/*genomicfwd_bothdiff*/genomicfwd_refdiff,
		    substring_length,/*querypos*/substring_start,querylength,
		    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
    FREE(genomicfwd_refdiff);

  } else {
    genomicfwd_refdiff = (char *) CALLOC(querylength+1,sizeof(char));
    make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring_start]),substring_length);
    genomicfwd_bothdiff = (char *) CALLOC(querylength+1,sizeof(char));
    make_complement_buffered(genomicfwd_bothdiff,&(genomicdir_bothdiff[substring_start]),substring_length);
    print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
		    fp,/*matchlength*/0,genomicfwd_refdiff,genomicfwd_bothdiff,
		    substring_length,/*querypos*/substring_start,querylength,
		    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
    FREE(genomicfwd_bothdiff);
    FREE(genomicfwd_refdiff);
  }

  if (printp == false) {
    fprintf(fp,"0");
  }
  

  /* 12. TAGS: NH */
  fprintf(fp,"\t");
  fprintf(fp,"NH:i:%d",npaths);

  /* 12. TAGS: HI */
  fprintf(fp,"\t");
  fprintf(fp,"HI:i:%d",pathnum);

  /* 12. TAGS: NM */
  fprintf(fp,"\t");
  /* fprintf(fp,"NM:i:%d",Stage3end_nmismatches_refdiff(this)); */
  fprintf(fp,"NM:i:%d",nmismatches_refdiff);

  if (snps_iit) {
    /* 12. TAGS: XW and XV */
    fprintf(fp,"\t");
    fprintf(fp,"XW:i:%d",nmismatches_bothdiff);
    fprintf(fp,"\t");
    fprintf(fp,"XV:i:%d",nmismatches_refdiff - nmismatches_bothdiff);
  }

  /* 12. TAGS: SM */
  fprintf(fp,"\t");
  fprintf(fp,"SM:i:%d",mapq_score);

  /* 12. TAGS: XQ */
  fprintf(fp,"\t");
  fprintf(fp,"XQ:i:%d",absmq_score);

  /* 12. TAGS: X2 */
  fprintf(fp,"\t");
  fprintf(fp,"X2:i:%d",second_absmq);

  /* 12. TAGS: XO */
  fprintf(fp,"\t");
  fprintf(fp,"XO:Z:%s",abbrev);

  /* 12. TAGS: PG */
  if (Stage3end_sarrayp(this) == true) {
    fprintf(fp,"\t");
    fprintf(fp,"PG:Z:A");
  } else if (Stage3end_hittype(this) == TERMINAL) {
    fprintf(fp,"\t");
    fprintf(fp,"PG:Z:T");
  }

  fprintf(fp,"\n");
  return;
}


static void
print_insertion (FILE *fp, char *abbrev, Stage3end_T this, Stage3end_T mate,
		 char *acc1, char *acc2, int pathnum, int npaths,
		 int absmq_score, int first_absmq, int second_absmq, int mapq_score,
		 Univ_IIT_T chromosome_iit, Shortread_T queryseq, int pairedlength,
		 Chrpos_T chrpos, Chrpos_T mate_chrpos, int clipdir, int hardclip5, int hardclip3,
		 Resulttype_T resulttype, bool first_read_p, int npaths_mate,
		 int quality_shift, char *sam_read_group_id, bool invertp, bool invert_mate_p,
		 bool circularp) {
  unsigned int flag = 0U;
  Substring_T substring1, substring2;
  int nmismatches_refdiff = 0, nmismatches_bothdiff = 0, querylength;
  char *genomicfwd_refdiff, *genomicfwd_bothdiff, *genomicdir_refdiff, *genomicdir_bothdiff;
  int substring1_start, substring2_start, substring1_length, substring2_length, matchlength, nindels;
  int querypos;
  int hardclip_low, hardclip_high;
  /* int mate_hardclip_low, mate_hardclip_high; */
  bool plusp, printp;
  List_T cigar_tokens = NULL;

  querylength = Shortread_fulllength(queryseq);
  plusp = Stage3end_plusp(this);

  substring1 = Stage3end_substring1(this);
  substring2 = Stage3end_substring2(this);

  if (circularp == true) {
    /* clipdir should be +1 */
    if (1 || plusp == true) {
      hardclip_low = hardclip5;
      hardclip_high = hardclip3;
      /* mate_hardclip_low = 0; */
      /* mate_hardclip_high = 0; */
    } else {
      hardclip_low = hardclip3;
      hardclip_high = hardclip5;
      /* mate_hardclip_low = 0; */
      /* mate_hardclip_high = 0; */
    }
  } else {
    if (first_read_p == true) {
      if (clipdir >= 0) {
	hardclip_low = 0;
	hardclip_high = hardclip5;
      } else {
	hardclip_low = hardclip5;
	hardclip_high = 0;
      }
      /* mate_hardclip_low = hardclip3; */
      /* mate_hardclip_high = 0; */
      /* fprintf(stderr,"first read: hardclip_low = %d, hardclip_high = %d\n",hardclip_low,hardclip_high); */
    } else {
      if (clipdir >= 0) {
	hardclip_low = hardclip3;
	hardclip_high = 0;
      } else {
	hardclip_low = 0;
	hardclip_high = hardclip3;
      }
      /* mate_hardclip_low = 0; */
      /* mate_hardclip_high = hardclip5; */
      /* fprintf(stderr,"second read: hardclip_low = %d, hardclip_high = %d\n",hardclip_low,hardclip_high); */
    }
  }

  nindels = Stage3end_nindels(this);

  /* 1. QNAME */
  if (acc2 == NULL) {
    fprintf(fp,"%s",acc1);
  } else {
    fprintf(fp,"%s,%s",acc1,acc2);
  }

  /* 2. FLAG */
  flag = SAM_compute_flag(plusp,mate,resulttype,first_read_p,
			  pathnum,npaths,npaths_mate,absmq_score,first_absmq,
			  invertp,invert_mate_p);
  fprintf(fp,"\t%u",flag);

  /* 3. RNAME: chr */
  /* 4. POS: chrpos */
  print_chromosomal_pos(fp,Stage3end_chrnum(this),chrpos,chromosome_iit);


  /* 5. MAPQ: Mapping quality */
  fprintf(fp,"\t%d",mapq_score);

  /* 6. CIGAR */
  fprintf(fp,"\t");

  if (plusp == true) {
    cigar_tokens = compute_cigar(cigar_tokens,/*type*/'S',Substring_querystart(substring1),
				 /*querypos*/0,querylength,hardclip_low,hardclip_high,
				 /*plusp*/true,/*lastp*/false);
    cigar_tokens = compute_cigar(cigar_tokens,/*type*/'M',Substring_match_length(substring1),
				 /*querypos*/Substring_querystart(substring1),querylength,
				 hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);
    cigar_tokens = compute_cigar(cigar_tokens,/*type*/'I',nindels,
				 /*querypos*/Substring_queryend(substring1),querylength,
				 hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);
    cigar_tokens = compute_cigar(cigar_tokens,/*type*/'M',Substring_match_length(substring2),
				 /*querypos*/Substring_querystart(substring2),querylength,
				 hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);
    cigar_tokens = compute_cigar(cigar_tokens,/*type*/'S',querylength - Substring_queryend(substring2),
				 /*querypos*/Substring_queryend(substring2),querylength,
				 hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);

  } else {
    cigar_tokens = compute_cigar(cigar_tokens,/*type*/'S',querylength - Substring_queryend(substring2),
				 /*querypos*/querylength,querylength,hardclip_low,hardclip_high,
				 /*plusp*/false,/*lastp*/false);
    cigar_tokens = compute_cigar(cigar_tokens,/*type*/'M',Substring_match_length(substring2),
				 /*querypos*/Substring_queryend(substring2),querylength,
				 hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
    cigar_tokens = compute_cigar(cigar_tokens,/*type*/'I',nindels,
				 /*querypos*/Substring_querystart(substring2),querylength,
				 hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
    cigar_tokens = compute_cigar(cigar_tokens,/*type*/'M',Substring_match_length(substring1),
				 /*querypos*/Substring_queryend(substring1),querylength,
				 hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
    cigar_tokens = compute_cigar(cigar_tokens,/*type*/'S',Substring_querystart(substring1),
				 /*querypos*/Substring_querystart(substring1),querylength,
				 hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
  }
  cigar_tokens = Pair_clean_cigar(cigar_tokens,/*watsonp*/true);
  print_tokens_sam(fp,cigar_tokens);
  List_free(&cigar_tokens);


  /* 7. MRNM: Mate chr */
  /* 8. MPOS: Mate chrpos */
  print_mate_chromosomal_pos(fp,Stage3end_chrnum(mate),Stage3end_effective_chrnum(mate),mate_chrpos,
			     Stage3end_chrnum(this),chrpos,chromosome_iit);


  /* 9. ISIZE: Insert size */
  if (resulttype == CONCORDANT_UNIQ || resulttype == CONCORDANT_TRANSLOC || resulttype == CONCORDANT_MULT) {
    if (plusp == invertp) {
      fprintf(fp,"\t%d",-pairedlength);
    } else {
      fprintf(fp,"\t%d",pairedlength);
    }
  } else if (mate_chrpos == 0) {
    fprintf(fp,"\t%d",pairedlength);
  } else if (chrpos < mate_chrpos) {
    fprintf(fp,"\t%d",pairedlength);
  } else if (chrpos > mate_chrpos) {
    fprintf(fp,"\t%d",-pairedlength);
  } else if (first_read_p == true) {
    fprintf(fp,"\t%d",pairedlength);
  } else {
    fprintf(fp,"\t%d",-pairedlength);
  }

  /* 10. SEQ: queryseq and 11. QUAL: quality scores */
  /* Queryseq has already been inverted, so just measure plusp relative to its current state */
  fprintf(fp,"\t");
  if (plusp == true) {
    Shortread_print_chopped(fp,queryseq,hardclip_low,hardclip_high);
    fprintf(fp,"\t");
    Shortread_print_quality(fp,queryseq,hardclip_low,hardclip_high,
			    quality_shift,/*show_chopped_p*/false);
  } else {
    Shortread_print_chopped_revcomp(fp,queryseq,hardclip_low,hardclip_high);
    fprintf(fp,"\t");
    Shortread_print_quality_revcomp(fp,queryseq,hardclip_low,hardclip_high,
				    quality_shift,/*show_chopped_p*/false);
  }

  /* 12. TAGS: RG */
  if (sam_read_group_id != NULL) {
    fprintf(fp,"\t");
    fprintf(fp,"RG:Z:%s",sam_read_group_id);
  }

  /* 12. TAGS: XB */
  Shortread_print_barcode(fp,queryseq);

  /* 12. TAGS: XP */
  Shortread_print_chop(fp,queryseq,invertp);

  /* 12. TAGS: MD */
  fprintf(fp,"\t");
  fprintf(fp,"MD:Z:");
  printp = false;

  substring1_start = Substring_querystart(substring1);
  substring1_length = Substring_match_length(substring1);
  substring2_start = Substring_querystart(substring2);
  substring2_length = Substring_match_length(substring2);

  if (plusp == true) {
    genomicfwd_refdiff = Substring_genomic_refdiff(substring1);
    genomicfwd_bothdiff = Substring_genomic_bothdiff(substring1);
    matchlength = print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,fp,/*matchlength*/0,
				  &(genomicfwd_refdiff[substring1_start]),&(genomicfwd_bothdiff[substring1_start]),
				  substring1_length,/*querypos*/substring1_start,querylength,
				  hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);

#if 0
    /* If MD string is supposed to include insertion, then uncomment this */
    matchlength += nindels;
#endif

    genomicfwd_refdiff = Substring_genomic_refdiff(substring2);
    genomicfwd_bothdiff = Substring_genomic_bothdiff(substring2);
    print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,fp,matchlength,
		    &(genomicfwd_refdiff[substring2_start]),&(genomicfwd_bothdiff[substring2_start]),
		    substring2_length,/*querypos*/substring2_start,querylength,
		    hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
  } else {
    genomicdir_refdiff = Substring_genomic_refdiff(substring2);
    genomicdir_bothdiff = Substring_genomic_bothdiff(substring2);
    if (genomicdir_bothdiff == genomicdir_refdiff) {
      genomicfwd_refdiff = (char *) CALLOC(substring2_length+1,sizeof(char));
      make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring2_start]),substring2_length);
      matchlength = print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
				    fp,/*matchlength*/0,genomicfwd_refdiff,/*genomicfwd_bothdiff*/genomicfwd_refdiff,
				    substring2_length,/*querypos*/substring2_start,querylength,
				    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
      FREE(genomicfwd_refdiff);
    } else {
      genomicfwd_refdiff = (char *) CALLOC(substring2_length+1,sizeof(char));
      make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring2_start]),substring2_length);
      genomicfwd_bothdiff = (char *) CALLOC(substring2_length+1,sizeof(char));
      make_complement_buffered(genomicfwd_bothdiff,&(genomicdir_bothdiff[substring2_start]),substring2_length);
      matchlength = print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
				    fp,/*matchlength*/0,genomicfwd_refdiff,genomicfwd_bothdiff,
				    substring2_length,/*querypos*/substring2_start,querylength,
				    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
      FREE(genomicfwd_bothdiff);
      FREE(genomicfwd_refdiff);
    }

#if 0
    /* If MD string is supposed to include insertion, then uncomment this */
    matchlength += nindels;
#endif

    genomicdir_refdiff = Substring_genomic_refdiff(substring1);
    genomicdir_bothdiff = Substring_genomic_bothdiff(substring1);
    if (genomicdir_bothdiff == genomicdir_refdiff) {
      genomicfwd_refdiff = (char *) CALLOC(substring1_length+1,sizeof(char));
      make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring1_start]),substring1_length);
      print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
		      fp,matchlength,genomicfwd_refdiff,/*genomicfwd_bothdiff*/genomicfwd_refdiff,
		      substring1_length,/*querypos*/substring1_start,querylength,
		      hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
      FREE(genomicfwd_refdiff);
    } else {
      genomicfwd_refdiff = (char *) CALLOC(substring1_length+1,sizeof(char));
      make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring1_start]),substring1_length);
      genomicfwd_bothdiff = (char *) CALLOC(substring1_length+1,sizeof(char));
      make_complement_buffered(genomicfwd_bothdiff,&(genomicdir_bothdiff[substring1_start]),substring1_length);
      print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
		      fp,matchlength,genomicfwd_refdiff,genomicfwd_bothdiff,
		      substring1_length,/*querypos*/substring1_start,querylength,
		      hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
      FREE(genomicfwd_bothdiff);
      FREE(genomicfwd_refdiff);
    }

  }
  if (printp == false) {
    fprintf(fp,"0");
  }


  /* 12. TAGS: NH */
  fprintf(fp,"\t");
  fprintf(fp,"NH:i:%d",npaths);

  /* 12. TAGS: HI */
  fprintf(fp,"\t");
  fprintf(fp,"HI:i:%d",pathnum);

  /* 12. TAGS: NM */
  fprintf(fp,"\t");
  /* fprintf(fp,"NM:i:%d",Stage3end_nmismatches_refdiff(this)); */
  fprintf(fp,"NM:i:%d",nmismatches_refdiff + nindels);

  if (snps_iit) {
    /* 12. TAGS: XW and XV */
    fprintf(fp,"\t");
    fprintf(fp,"XW:i:%d",nmismatches_bothdiff);
    fprintf(fp,"\t");
    fprintf(fp,"XV:i:%d",nmismatches_refdiff - nmismatches_bothdiff);
  }

  /* 12. TAGS: SM */
  fprintf(fp,"\t");
  fprintf(fp,"SM:i:%d",mapq_score);

  /* 12. TAGS: XQ */
  fprintf(fp,"\t");
  fprintf(fp,"XQ:i:%d",absmq_score);

  /* 12. TAGS: X2 */
  fprintf(fp,"\t");
  fprintf(fp,"X2:i:%d",second_absmq);

  /* 12. TAGS: XO */
  fprintf(fp,"\t");
  fprintf(fp,"XO:Z:%s",abbrev);

  /* 12. TAGS: PG */
  if (Stage3end_sarrayp(this) == true) {
    fprintf(fp,"\t");
    fprintf(fp,"PG:Z:A");
  }

  fprintf(fp,"\n");
  return;
}


static void
print_deletion (FILE *fp, char *abbrev, Stage3end_T this, Stage3end_T mate,
		char *acc1, char *acc2, int pathnum, int npaths,
		int absmq_score, int first_absmq, int second_absmq, int mapq_score,
		Univ_IIT_T chromosome_iit, Shortread_T queryseq, int pairedlength,
		Chrpos_T chrpos, Chrpos_T mate_chrpos, int clipdir, int hardclip5, int hardclip3,
		Resulttype_T resulttype, bool first_read_p, int npaths_mate,
		int quality_shift, char *sam_read_group_id, bool invertp, bool invert_mate_p,
		bool circularp) {
  unsigned int flag = 0U;
  Substring_T substring1, substring2;
  int nmismatches_refdiff = 0, nmismatches_bothdiff = 0, querylength;
  char *genomicfwd_refdiff, *genomicfwd_bothdiff, *genomicfwd_deletion,
    *genomicdir_refdiff, *genomicdir_bothdiff;
  int substring1_start, substring2_start, substring1_length, substring2_length, nindels;
  int hardclip_low, hardclip_high;
  /* int mate_hardclip_low, mate_hardclip_high; */
  bool plusp, printp;

  querylength = Shortread_fulllength(queryseq);
  plusp = Stage3end_plusp(this);

  substring1 = Stage3end_substring1(this);
  substring2 = Stage3end_substring2(this);

  if (circularp == true) {
    /* clipdir should be +1 */
    if (1 || plusp == true) {
      hardclip_low = hardclip5;
      hardclip_high = hardclip3;
      /* mate_hardclip_low = 0; */
      /* mate_hardclip_high = 0; */
    } else {
      hardclip_low = hardclip3;
      hardclip_high = hardclip5;
      /* mate_hardclip_low = 0; */
      /* mate_hardclip_high = 0; */
    }
  } else {
    if (first_read_p == true) {
      if (clipdir >= 0) {
	hardclip_low = 0;
	hardclip_high = hardclip5;
      } else {
	hardclip_low = hardclip5;
	hardclip_high = 0;
      }
      /* mate_hardclip_low = hardclip3; */
      /* mate_hardclip_high = 0; */
      /* fprintf(stderr,"first read: hardclip_low = %d, hardclip_high = %d\n",hardclip_low,hardclip_high); */
    } else {
      if (clipdir >= 0) {
	hardclip_low = hardclip3;
	hardclip_high = 0;
      } else {
	hardclip_low = 0;
	hardclip_high = hardclip3;
      }
      /* mate_hardclip_low = 0; */
      /* mate_hardclip_high = hardclip5; */
      /* fprintf(stderr,"second read: hardclip_low = %d, hardclip_high = %d\n",hardclip_low,hardclip_high); */
    }
  }

#if 0
  /* These cases are checked below */
  if (hardclip_low >= Substring_querystart(substring2)) {
    nindels = 0;
  } else if (querylength - hardclip_high <= Substring_queryend(substring1)) {
    nindels = 0;
  } else {
    nindels = Stage3end_nindels(this); /* nindels is positive */
  }
#else
  nindels = Stage3end_nindels(this); /* nindels is positive */
#endif


  /* 1. QNAME */
  if (acc2 == NULL) {
    fprintf(fp,"%s",acc1);
  } else {
    fprintf(fp,"%s,%s",acc1,acc2);
  }

  /* 2. FLAG */
  flag = SAM_compute_flag(plusp,mate,resulttype,first_read_p,
			  pathnum,npaths,npaths_mate,absmq_score,first_absmq,
			  invertp,invert_mate_p);
  fprintf(fp,"\t%u",flag);

  /* 3. RNAME: chr */
  /* 4. POS: chrpos */
  print_chromosomal_pos(fp,Stage3end_chrnum(this),chrpos,chromosome_iit);


  /* 5. MAPQ: Mapping quality */
  fprintf(fp,"\t%d",mapq_score);

  /* 6. CIGAR */
  fprintf(fp,"\t");

  substring1_start = Substring_querystart(substring1);
  substring1_length = Substring_match_length(substring1);
  substring2_start = Substring_querystart(substring2);
  substring2_length = Substring_match_length(substring2);

  if (plusp == true) {
    print_cigar(fp,/*type*/'S',Substring_querystart(substring1),
		/*querypos*/0,querylength,hardclip_low,hardclip_high,
		/*plusp*/true,/*lastp*/false);
    if (/*nindels > 0 &&*/ hardclip_low < substring1_start + substring1_length && hardclip_high < querylength - substring2_start) {
      print_cigar(fp,/*type*/'M',Substring_match_length(substring1),
		  /*querypos*/Substring_querystart(substring1),querylength,
		  hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);
      fprintf(fp,"%dD",nindels);
      print_cigar(fp,/*type*/'M',Substring_match_length(substring2),
		  /*querypos*/Substring_querystart(substring2),querylength,
		  hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);
    } else {
      print_cigar(fp,/*type*/'M',Substring_match_length(substring1) + Substring_match_length(substring2),
		  /*querypos*/Substring_querystart(substring1),querylength,
		  hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);
    }
    print_cigar(fp,/*type*/'S',querylength - Substring_queryend(substring2),
		/*querypos*/Substring_queryend(substring2),querylength,
		hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);

  } else {
    print_cigar(fp,/*type*/'S',querylength - Substring_queryend(substring2),
		/*querypos*/querylength,querylength,
		hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
    if (/*nindels > 0 &&*/ hardclip_low < querylength - substring2_start && hardclip_high < substring1_start + substring1_length) {
      print_cigar(fp,/*type*/'M',Substring_match_length(substring2),
		  /*querypos*/Substring_queryend(substring2),querylength,
		  hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
      fprintf(fp,"%dD",nindels);
      print_cigar(fp,/*type*/'M',Substring_match_length(substring1),
		  /*querypos*/Substring_querystart(substring2),querylength,
		  hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
    } else {
      print_cigar(fp,/*type*/'M',Substring_match_length(substring2) + Substring_match_length(substring1),
		  /*querypos*/Substring_queryend(substring2),querylength,
		  hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
    }
    print_cigar(fp,/*type*/'S',Substring_querystart(substring1),
		/*querypos*/Substring_querystart(substring1),querylength,
		hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
  }

  /* 7. MRNM: Mate chr */
  /* 8. MPOS: Mate chrpos */
  print_mate_chromosomal_pos(fp,Stage3end_chrnum(mate),Stage3end_effective_chrnum(mate),mate_chrpos,
			     Stage3end_chrnum(this),chrpos,chromosome_iit);


  /* 9. ISIZE: Insert size */
  if (resulttype == CONCORDANT_UNIQ || resulttype == CONCORDANT_TRANSLOC || resulttype == CONCORDANT_MULT) {
    if (plusp == invertp) {
      fprintf(fp,"\t%d",-pairedlength);
    } else {
      fprintf(fp,"\t%d",pairedlength);
    }
  } else if (mate_chrpos == 0) {
    fprintf(fp,"\t%d",pairedlength);
  } else if (chrpos < mate_chrpos) {
    fprintf(fp,"\t%d",pairedlength);
  } else if (chrpos > mate_chrpos) {
    fprintf(fp,"\t%d",-pairedlength);
  } else if (first_read_p == true) {
    fprintf(fp,"\t%d",pairedlength);
  } else {
    fprintf(fp,"\t%d",-pairedlength);
  }


  /* 10. SEQ: queryseq and 11. QUAL: quality scores */
  /* Queryseq has already been inverted, so just measure plusp relative to its current state */
  fprintf(fp,"\t");
  if (plusp == true) {
    Shortread_print_chopped(fp,queryseq,hardclip_low,hardclip_high);
    fprintf(fp,"\t");
    Shortread_print_quality(fp,queryseq,hardclip_low,hardclip_high,
			    quality_shift,/*show_chopped_p*/false);
  } else {
    Shortread_print_chopped_revcomp(fp,queryseq,hardclip_low,hardclip_high);
    fprintf(fp,"\t");
    Shortread_print_quality_revcomp(fp,queryseq,hardclip_low,hardclip_high,
				    quality_shift,/*show_chopped_p*/false);
  }

  /* 12. TAGS: RG */
  if (sam_read_group_id != NULL) {
    fprintf(fp,"\t");
    fprintf(fp,"RG:Z:%s",sam_read_group_id);
  }

  /* 12. TAGS: XB */
  Shortread_print_barcode(fp,queryseq);

  /* 12. TAGS: XP */
  Shortread_print_chop(fp,queryseq,invertp);

  /* 12. TAGS: MD */
  fprintf(fp,"\t");
  fprintf(fp,"MD:Z:");
  printp = false;

  if (plusp == true) {
    genomicfwd_refdiff = Substring_genomic_refdiff(substring1);
    genomicfwd_bothdiff = Substring_genomic_bothdiff(substring1);
    print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,fp,/*matchlength*/0,
		    &(genomicfwd_refdiff[substring1_start]),&(genomicfwd_bothdiff[substring1_start]),
		    substring1_length,/*querypos*/substring1_start,querylength,
		    hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);

    debug1(printf("\nhardclip_low %d, hardclip_high %d\n",hardclip_low,hardclip_high));
    debug1(printf("substring1_length %d, substring2_length %d\n",substring1_length,substring2_length));
    debug1(printf("substring1 %d..%d, substring2 %d..%d\n",
		  Substring_querystart(substring1),Substring_queryend(substring1),
		  Substring_querystart(substring2),Substring_queryend(substring2)));
    debug1(printf("trim1: %d..%d, trim2 %d..%d\n",
		  Substring_trim_left(substring1),Substring_trim_right(substring1),
		  Substring_trim_left(substring2),Substring_trim_right(substring2)));
    if (hardclip_low < substring1_start + substring1_length && hardclip_high < querylength - substring2_start) {
      fprintf(fp,"^%s",Stage3end_deletion_string(this));
    }

    genomicfwd_refdiff = Substring_genomic_refdiff(substring2);
    genomicfwd_bothdiff = Substring_genomic_bothdiff(substring2);
    print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,fp,/*matchlength*/0,
		    &(genomicfwd_refdiff[substring2_start]),&(genomicfwd_bothdiff[substring2_start]),
		    substring2_length,/*querypos*/substring2_start,querylength,
		    hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);

  } else {
    genomicdir_refdiff = Substring_genomic_refdiff(substring2);
    genomicdir_bothdiff = Substring_genomic_bothdiff(substring2);
    if (genomicdir_bothdiff == genomicdir_refdiff) {
      genomicfwd_refdiff = (char *) CALLOC(substring2_length+1,sizeof(char));
      make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring2_start]),substring2_length);
      print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
		      fp,/*matchlength*/0,genomicfwd_refdiff,/*genomicfwd_bothdiff*/genomicfwd_refdiff,
		      substring2_length,/*querypos*/substring2_start,querylength,
		      hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
      FREE(genomicfwd_refdiff);
    } else {
      genomicfwd_refdiff = (char *) CALLOC(substring2_length+1,sizeof(char));
      make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring2_start]),substring2_length);
      genomicfwd_bothdiff = (char *) CALLOC(substring2_length+1,sizeof(char));
      make_complement_buffered(genomicfwd_bothdiff,&(genomicdir_bothdiff[substring2_start]),substring2_length);
      print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
		      fp,/*matchlength*/0,genomicfwd_refdiff,genomicfwd_bothdiff,
		      substring2_length,/*querypos*/substring2_start,querylength,
		      hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
      FREE(genomicfwd_bothdiff);
      FREE(genomicfwd_refdiff);
    }


    debug1(printf("\nhardclip_low %d, hardclip_high %d\n",hardclip_low,hardclip_high));
    debug1(printf("substring2_length %d, substring1_length %d\n",substring2_length,substring1_length));
    debug1(printf("substring1 %d..%d, substring2 %d..%d\n",
		  Substring_querystart(substring1),Substring_queryend(substring1),
		  Substring_querystart(substring2),Substring_queryend(substring2)));
    debug1(printf("trim1: %d..%d, trim2 %d..%d\n",
		  Substring_trim_left(substring1),Substring_trim_right(substring1),
		  Substring_trim_left(substring2),Substring_trim_right(substring2)));
    if (hardclip_low < querylength - substring2_start && hardclip_high < substring1_start + substring1_length) {
      /* Deletion string: Potential problem if followed by a mismatch, but can be resolved by looking at CIGAR string */
      genomicfwd_deletion = (char *) CALLOC(nindels+1,sizeof(char));
      make_complement_buffered(genomicfwd_deletion,Stage3end_deletion_string(this),nindels);
      fprintf(fp,"^%s",genomicfwd_deletion);
      FREE(genomicfwd_deletion);
    }

    genomicdir_refdiff = Substring_genomic_refdiff(substring1);
    genomicdir_bothdiff = Substring_genomic_bothdiff(substring1);
    if (genomicdir_bothdiff == genomicdir_refdiff) {
      genomicfwd_refdiff = (char *) CALLOC(substring1_length+1,sizeof(char));
      make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring1_start]),substring1_length);
      print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
		      fp,/*matchlength*/0,genomicfwd_refdiff,/*genomicfwd_bothdiff*/genomicfwd_refdiff,
		      substring1_length,/*querypos*/substring1_start,querylength,
		      hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
      FREE(genomicfwd_refdiff);
    } else {
      genomicfwd_refdiff = (char *) CALLOC(substring1_length+1,sizeof(char));
      make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring1_start]),substring1_length);
      genomicfwd_bothdiff = (char *) CALLOC(substring1_length+1,sizeof(char));
      make_complement_buffered(genomicfwd_bothdiff,&(genomicdir_bothdiff[substring1_start]),substring1_length);
      print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
		      fp,/*matchlength*/0,genomicfwd_refdiff,genomicfwd_bothdiff,
		      substring1_length,/*querypos*/substring1_start,querylength,
		      hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
      FREE(genomicfwd_bothdiff);
      FREE(genomicfwd_refdiff);
    }

  }
  if (printp == false) {
    fprintf(fp,"0");
  }


  /* 12. TAGS: NH */
  fprintf(fp,"\t");
  fprintf(fp,"NH:i:%d",npaths);

  /* 12. TAGS: HI */
  fprintf(fp,"\t");
  fprintf(fp,"HI:i:%d",pathnum);

  /* 12. TAGS: NM */
  fprintf(fp,"\t");
  /* fprintf(fp,"NM:i:%d",Stage3end_nmismatches_refdiff(this)); */
  fprintf(fp,"NM:i:%d",nmismatches_refdiff + nindels);

  if (snps_iit) {
    /* 12. TAGS: XW and XV */
    fprintf(fp,"\t");
    fprintf(fp,"XW:i:%d",nmismatches_bothdiff);
    fprintf(fp,"\t");
    fprintf(fp,"XV:i:%d",nmismatches_refdiff - nmismatches_bothdiff);
  }

  /* 12. TAGS: SM */
  fprintf(fp,"\t");
  fprintf(fp,"SM:i:%d",mapq_score);

  /* 12. TAGS: XQ */
  fprintf(fp,"\t");
  fprintf(fp,"XQ:i:%d",absmq_score);

  /* 12. TAGS: X2 */
  fprintf(fp,"\t");
  fprintf(fp,"X2:i:%d",second_absmq);

  /* 12. TAGS: XO */
  fprintf(fp,"\t");
  fprintf(fp,"XO:Z:%s",abbrev);

  /* 12. TAGS: PG */
  if (Stage3end_sarrayp(this) == true) {
    fprintf(fp,"\t");
    fprintf(fp,"PG:Z:A");
  }

  fprintf(fp,"\n");
  return;
}


static void
halfdonor_dinucleotide (char *donor1, char *donor2, Substring_T donor) {
  bool sensep;
  char *genomic;
  int substring_start, substring_length;

  /* sensedir for chimera must be SENSE_FORWARD or SENSE_ANTI, not SENSE_NULL */
  sensep = Substring_chimera_sensep(donor);

  substring_start = Substring_querystart(donor);
  genomic = Substring_genomic_refdiff(donor);

  if (sensep == true) {
    substring_length = Substring_match_length(donor);
    *donor1 = toupper(genomic[substring_start+substring_length]);
    *donor2 = toupper(genomic[substring_start+substring_length+1]);

  } else {			/* sensep == false */
    *donor2 = toupper(complCode[(int) genomic[substring_start-2]]);
    *donor1 = toupper(complCode[(int) genomic[substring_start-1]]);
  }

  return;
}

static void
halfacceptor_dinucleotide (char *acceptor2, char *acceptor1, Substring_T acceptor) {
  bool sensep;
  char *genomic;
  int substring_start, substring_length;

  /* sensedir for chimera must be SENSE_FORWARD or SENSE_ANTI, not SENSE_NULL */
  sensep = Substring_chimera_sensep(acceptor);

  substring_start = Substring_querystart(acceptor);
  genomic = Substring_genomic_refdiff(acceptor);

  if (sensep == true) {
    *acceptor2 = toupper(genomic[substring_start-2]);
    *acceptor1 = toupper(genomic[substring_start-1]);

  } else {			/* sensep == false */
    substring_length = Substring_match_length(acceptor);
    *acceptor1 = toupper(complCode[(int) genomic[substring_start+substring_length]]);
    *acceptor2 = toupper(complCode[(int) genomic[substring_start+substring_length+1]]);
  }

  return;
}



static void
print_halfdonor (FILE *fp, char *abbrev, Substring_T donor, Stage3end_T this, Stage3end_T mate,
		 char *acc1, char *acc2, int pathnum, int npaths, int absmq_score, int first_absmq, int second_absmq, int mapq_score,
		 Univ_IIT_T chromosome_iit, Shortread_T queryseq, int pairedlength,
		 Chrpos_T concordant_chrpos, Chrpos_T chrpos, Chrpos_T mate_chrpos,
		 int clipdir, int hardclip5, int hardclip3, Resulttype_T resulttype, bool first_read_p, int npaths_mate,
		 int quality_shift, char *sam_read_group_id, bool invertp, bool invert_mate_p,
		 bool use_hardclip_p, bool print_xt_p, char donor1, char donor2, char acceptor2, char acceptor1,
		 double donor_prob, double acceptor_prob, bool circularp) {
  unsigned int flag = 0U;
  int nmismatches_refdiff = 0, nmismatches_bothdiff = 0, querylength;
  int sensedir;
  bool sensep;
  char *genomicfwd_refdiff, *genomicfwd_bothdiff, *genomicdir_refdiff, *genomicdir_bothdiff;
  int substring_start, substring_length;
  int hardclip_low, hardclip_high;
  /* int mate_hardclip_low, mate_hardclip_high; */
  int transloc_hardclip_low, transloc_hardclip_high;
  bool plusp, printp;

  querylength = Shortread_fulllength(queryseq);
  plusp = Substring_plusp(donor);

  if (circularp == true) {
    /* clipdir should be +1 */
    if (1 || plusp == true) {
      hardclip_low = hardclip5;
      hardclip_high = hardclip3;
      /* mate_hardclip_low = 0; */
      /* mate_hardclip_high = 0; */
    } else {
      hardclip_low = hardclip3;
      hardclip_high = hardclip5;
      /* mate_hardclip_low = 0; */
      /* mate_hardclip_high = 0; */
    }
  } else {
    if (first_read_p == true) {
      if (clipdir >= 0) {
	hardclip_low = 0;
	hardclip_high = hardclip5;
      } else {
	hardclip_low = hardclip5;
	hardclip_high = 0;
      }
      /* mate_hardclip_low = hardclip3; */
      /* mate_hardclip_high = 0; */
      /* fprintf(stderr,"first read: hardclip_low = %d, hardclip_high = %d\n",hardclip_low,hardclip_high); */
    } else {
      if (clipdir >= 0) {
	hardclip_low = hardclip3;
	hardclip_high = 0;
      } else {
	hardclip_low = 0;
	hardclip_high = hardclip3;
      }
      /* mate_hardclip_low = 0; */
      /* mate_hardclip_high = hardclip5; */
      /* fprintf(stderr,"second read: hardclip_low = %d, hardclip_high = %d\n",hardclip_low,hardclip_high); */
    }
  }


  /* 1. QNAME */
  if (acc2 == NULL) {
    fprintf(fp,"%s",acc1);
  } else {
    fprintf(fp,"%s,%s",acc1,acc2);
  }

  /* 2. FLAG */
  flag = SAM_compute_flag(plusp,mate,resulttype,first_read_p,
			  pathnum,npaths,npaths_mate,absmq_score,first_absmq,
			  invertp,invert_mate_p);
  fprintf(fp,"\t%u",flag);

  /* 3. RNAME: chr */
  /* 4. POS: chrpos */
  print_chromosomal_pos(fp,Substring_chrnum(donor),chrpos,chromosome_iit);
  

  /* 5. MAPQ: Mapping quality */
  fprintf(fp,"\t%d",mapq_score);

  /* 6. CIGAR */
  fprintf(fp,"\t");
  /* sensedir for chimera must be SENSE_FORWARD or SENSE_ANTI, not SENSE_NULL */
  sensedir = Substring_chimera_sensedir(donor);
  sensep = Substring_chimera_sensep(donor);

  if (use_hardclip_p == true) {
    if (sensep == plusp) {
      transloc_hardclip_low = 0;
      if (plusp == true) {
	/* sensep true */
	transloc_hardclip_high = querylength - Substring_queryend(donor);

      } else {
	/* sensep false */
	transloc_hardclip_high = Substring_querystart(donor);
      }

    } else { /* sensep != Substring_plusp(donor) */
      transloc_hardclip_high = 0;
      if (plusp == true) {
	transloc_hardclip_low = Substring_querystart(donor);

      } else {
	transloc_hardclip_low = querylength - Substring_queryend(donor);
      }
    }

    if (transloc_hardclip_low > hardclip_low) {
      hardclip_low = transloc_hardclip_low;
    }
    if (transloc_hardclip_high > hardclip_high) {
      hardclip_high = transloc_hardclip_high;
    }
  }


  if (sensep == plusp) {
    if (plusp == true) {
      /* sensep true */
      assert(Substring_chimera_pos(donor) == Substring_queryend(donor));
      print_cigar(fp,/*type*/'S',Substring_querystart(donor),
		  /*querypos*/0,querylength,hardclip_low,hardclip_high,
		  /*plusp*/true,/*lastp*/false);
      print_cigar(fp,/*type*/'M',Substring_match_length(donor),
		  /*querypos*/Substring_querystart(donor),querylength,
		  hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);
      print_cigar(fp,/*type*/use_hardclip_p ? 'H' : 'S',querylength - Substring_queryend(donor),
		  /*querypos*/Substring_queryend(donor),querylength,
		  hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);

    } else {
      /* sensep false */
      assert(Substring_chimera_pos(donor) == Substring_querystart(donor));
      print_cigar(fp,/*type*/'S',querylength - Substring_queryend(donor),
		  /*querypos*/querylength,querylength,
		  hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
      print_cigar(fp,/*type*/'M',Substring_match_length(donor),
		  /*querypos*/Substring_queryend(donor),querylength,
		  hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
      print_cigar(fp,/*type*/use_hardclip_p ? 'H' : 'S',Substring_querystart(donor),
		  /*querypos*/Substring_querystart(donor),querylength,
		  hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
    }

  } else { /* sensep != Substring_plusp(donor) */
    if (plusp == true) {
      assert(Substring_chimera_pos(donor) == Substring_querystart(donor));
      print_cigar(fp,/*type*/use_hardclip_p ? 'H' : 'S',Substring_querystart(donor),
		  /*querypos*/0,querylength,hardclip_low,hardclip_high,
		  /*plusp*/true,/*lastp*/false);
      print_cigar(fp,/*type*/'M',Substring_match_length(donor),
		  /*querypos*/Substring_querystart(donor),querylength,
		  hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);
      print_cigar(fp,/*type*/'S',querylength - Substring_queryend(donor),
		  /*querypos*/Substring_queryend(donor),querylength,
		  hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);

    } else {
      assert(Substring_chimera_pos(donor) == Substring_queryend(donor));
      print_cigar(fp,/*type*/use_hardclip_p ? 'H' : 'S',querylength - Substring_queryend(donor),
		  /*querypos*/querylength,querylength,
		  hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
      print_cigar(fp,/*type*/'M',Substring_match_length(donor),
		  /*querypos*/Substring_queryend(donor),querylength,
		  hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
      print_cigar(fp,/*type*/'S',Substring_querystart(donor),
		  /*querypos*/Substring_querystart(donor),querylength,hardclip_low,hardclip_high,
		  /*plusp*/false,/*lastp*/true);
    }
  }

  /* 7. MRNM: Mate chr */
  /* 8. MPOS: Mate chrpos */
  /* For anchor_chrnum, previously used Stage3end_chrnum(this), but this is 0 */
  print_mate_chromosomal_pos(fp,Stage3end_chrnum(mate),Stage3end_effective_chrnum(mate),mate_chrpos,
			     /*anchor_chrnum*/Substring_chrnum(donor),chrpos,chromosome_iit);


  /* 9. ISIZE: Insert size */
  if (resulttype == CONCORDANT_UNIQ || resulttype == CONCORDANT_TRANSLOC || resulttype == CONCORDANT_MULT) {
    if (plusp == invertp) {
      fprintf(fp,"\t%d",-pairedlength);
    } else {
      fprintf(fp,"\t%d",pairedlength);
    }
  } else if (mate_chrpos == 0) {
    fprintf(fp,"\t%d",pairedlength);
  } else if (concordant_chrpos < mate_chrpos) {
    fprintf(fp,"\t%d",pairedlength);
  } else if (concordant_chrpos > mate_chrpos) {
    fprintf(fp,"\t%d",-pairedlength);
  } else if (first_read_p == true) {
    fprintf(fp,"\t%d",pairedlength);
  } else {
    fprintf(fp,"\t%d",-pairedlength);
  }


  /* 10. SEQ: queryseq and 11. QUAL: quality scores */
  /* Queryseq has already been inverted, so just measure plusp relative to its current state */
  fprintf(fp,"\t");
  if (plusp == true) {
    Shortread_print_chopped(fp,queryseq,hardclip_low,hardclip_high);
    fprintf(fp,"\t");
    Shortread_print_quality(fp,queryseq,hardclip_low,hardclip_high,
			    quality_shift,/*show_chopped_p*/false);
  } else {
    Shortread_print_chopped_revcomp(fp,queryseq,hardclip_low,hardclip_high);
    fprintf(fp,"\t");
    Shortread_print_quality_revcomp(fp,queryseq,hardclip_low,hardclip_high,
				    quality_shift,/*show_chopped_p*/false);
  }


  /* 12. TAGS: RG */
  if (sam_read_group_id != NULL) {
    fprintf(fp,"\t");
    fprintf(fp,"RG:Z:%s",sam_read_group_id);
  }

  /* 12. TAGS: XB */
  Shortread_print_barcode(fp,queryseq);

  /* 12. TAGS: XP */
  Shortread_print_chop(fp,queryseq,invertp);

  /* 12. TAGS: MD */
  fprintf(fp,"\t");
  fprintf(fp,"MD:Z:");
  printp = false;

  substring_start = Substring_querystart(donor);
  substring_length = Substring_match_length(donor);

  if (use_hardclip_p == false) {
    genomicdir_refdiff = Substring_genomic_refdiff(donor);
    genomicdir_bothdiff = Substring_genomic_bothdiff(donor);
    if (plusp == true) {
      print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,fp,/*matchlength*/0,
		      &(genomicdir_refdiff[substring_start]),&(genomicdir_bothdiff[substring_start]),
		      substring_length,/*querypos*/substring_start,querylength,
		      hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
    } else if (genomicdir_bothdiff == genomicdir_refdiff) {
      genomicfwd_refdiff = (char *) CALLOC(querylength+1,sizeof(char));
      make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring_start]),substring_length);
      print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
		      fp,/*matchlength*/0,genomicfwd_refdiff,/*genomicfwd_bothdiff*/genomicfwd_refdiff,
		      substring_length,/*querypos*/substring_start,querylength,
		      hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
      FREE(genomicfwd_refdiff);
    } else {
      genomicfwd_refdiff = (char *) CALLOC(querylength+1,sizeof(char));
      make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring_start]),substring_length);
      genomicfwd_bothdiff = (char *) CALLOC(querylength+1,sizeof(char));
      make_complement_buffered(genomicfwd_bothdiff,&(genomicdir_bothdiff[substring_start]),substring_length);
      print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
		      fp,/*matchlength*/0,genomicfwd_refdiff,genomicfwd_bothdiff,
		      substring_length,/*querypos*/substring_start,querylength,
		      hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
      FREE(genomicfwd_bothdiff);
      FREE(genomicfwd_refdiff);
    }

  } else if (sensep == true) {
    if (plusp == true) {
      genomicfwd_refdiff = Substring_genomic_refdiff(donor);
      genomicfwd_bothdiff = Substring_genomic_bothdiff(donor);
      print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,fp,/*matchlength*/0,
		      &(genomicfwd_refdiff[substring_start]),&(genomicfwd_bothdiff[substring_start]),
		      substring_length,/*querypos*/substring_start,querylength,
		      hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
    } else {
      genomicdir_refdiff = Substring_genomic_refdiff(donor);
      genomicdir_bothdiff = Substring_genomic_bothdiff(donor);
      if (genomicdir_bothdiff == genomicdir_refdiff) {
	genomicfwd_refdiff = (char *) CALLOC(substring_length+1,sizeof(char));
	make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring_start]),substring_length);
	print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
			fp,/*matchlength*/0,genomicfwd_refdiff,/*genomicfwd_bothdiff*/genomicfwd_refdiff,
			substring_length,/*querypos*/substring_start,querylength,
			hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
	FREE(genomicfwd_refdiff);
      } else {
	genomicfwd_refdiff = (char *) CALLOC(substring_length+1,sizeof(char));
	make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring_start]),substring_length);
	genomicfwd_bothdiff = (char *) CALLOC(substring_length+1,sizeof(char));
	make_complement_buffered(genomicfwd_bothdiff,&(genomicdir_bothdiff[substring_start]),substring_length);
	print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
			fp,/*matchlength*/0,genomicfwd_refdiff,genomicfwd_bothdiff,
			substring_length,/*querypos*/substring_start,querylength,
			hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
	FREE(genomicfwd_bothdiff);
	FREE(genomicfwd_refdiff);
      }
    }

  } else {			/* sensep == false */
    if (plusp == true) {
      genomicfwd_refdiff = Substring_genomic_refdiff(donor);
      genomicfwd_bothdiff = Substring_genomic_bothdiff(donor);
      print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,fp,/*matchlength*/0,
		      &(genomicfwd_refdiff[substring_start]),&(genomicfwd_bothdiff[substring_start]),
		      substring_length,/*querypos*/substring_start,querylength,
		      hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
    } else {
      genomicdir_refdiff = Substring_genomic_refdiff(donor);
      genomicdir_bothdiff = Substring_genomic_refdiff(donor);
      if (genomicdir_bothdiff == genomicdir_refdiff) {
	genomicfwd_refdiff = (char *) CALLOC(substring_length+1,sizeof(char));
	make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring_start]),substring_length);
	print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
			fp,/*matchlength*/0,genomicfwd_refdiff,/*genomicfwd_bothdiff*/genomicfwd_refdiff,
			substring_length,/*querypos*/substring_start,querylength,
			hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
	FREE(genomicfwd_refdiff);
      } else {
	genomicfwd_refdiff = (char *) CALLOC(substring_length+1,sizeof(char));
	make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring_start]),substring_length);
	genomicfwd_bothdiff = (char *) CALLOC(substring_length+1,sizeof(char));
	make_complement_buffered(genomicfwd_bothdiff,&(genomicdir_bothdiff[substring_start]),substring_length);
	print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
			fp,/*matchlength*/0,genomicfwd_refdiff,genomicfwd_bothdiff,
			substring_length,/*querypos*/substring_start,querylength,
			hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
	FREE(genomicfwd_bothdiff);
	FREE(genomicfwd_refdiff);
      }
    }
  }
  if (printp == false) {
    fprintf(fp,"0");
  }


  /* 12. TAGS: NH */
  fprintf(fp,"\t");
  fprintf(fp,"NH:i:%d",npaths);

  /* 12. TAGS: HI */
  fprintf(fp,"\t");
  fprintf(fp,"HI:i:%d",pathnum);

  /* 12. TAGS: NM */
  fprintf(fp,"\t");
  /* fprintf(fp,"NM:i:%d",Substring_nmismatches_refdiff(donor)); */
  fprintf(fp,"NM:i:%d",nmismatches_refdiff);
  
  if (snps_iit) {
    /* 12. TAGS: XW and XV */
    fprintf(fp,"\t");
    fprintf(fp,"XW:i:%d",nmismatches_bothdiff);
    fprintf(fp,"\t");
    fprintf(fp,"XV:i:%d",nmismatches_refdiff - nmismatches_bothdiff);
  }

  /* 12. TAGS: SM */
  fprintf(fp,"\t");
  fprintf(fp,"SM:i:%d",mapq_score);

  /* 12. TAGS: XQ */
  fprintf(fp,"\t");
  fprintf(fp,"XQ:i:%d",absmq_score);

  /* 12. TAGS: X2 */
  fprintf(fp,"\t");
  fprintf(fp,"X2:i:%d",second_absmq);

  /* 12. TAGS: XO */
  fprintf(fp,"\t");
  fprintf(fp,"XO:Z:%s",abbrev);

  /* 12. TAGS: XS */
  fprintf(fp,"\t");
  /* sensedir for chimera must be SENSE_FORWARD or SENSE_ANTI, not SENSE_NULL */
  if (sensedir == SENSE_FORWARD) {
    if (plusp == true) {
      fprintf(fp,"XS:A:+");
    } else {
      fprintf(fp,"XS:A:-");
    }
  } else if (sensedir == SENSE_ANTI) {
    if (plusp == true) {
      fprintf(fp,"XS:A:-");
    } else {
      fprintf(fp,"XS:A:+");
    }
  } else if (force_xs_direction_p == true) {
    fprintf(fp,"XS:A:+");
  } else {
    fprintf(fp,"XS:A:?");
  }

  /* 12. TAGS: XT */
  if (print_xt_p == true) {
    fprintf(fp,"\t");
    fprintf(fp,"XT:Z:%c%c-%c%c,%.2f,%.2f",donor1,donor2,acceptor2,acceptor1,donor_prob,acceptor_prob);
  }

  /* 12. TAGS: PG */
  if (Stage3end_sarrayp(this) == true) {
    fprintf(fp,"\t");
    fprintf(fp,"PG:Z:A");
  }

  fprintf(fp,"\n");
  return;
}


static void
print_halfacceptor (FILE *fp, char *abbrev, Substring_T acceptor, Stage3end_T this, Stage3end_T mate,
		    char *acc1, char *acc2, int pathnum, int npaths, int absmq_score, int first_absmq, int second_absmq, int mapq_score,
		    Univ_IIT_T chromosome_iit, Shortread_T queryseq, int pairedlength,
		    Chrpos_T concordant_chrpos, Chrpos_T chrpos, Chrpos_T mate_chrpos,
		    int clipdir, int hardclip5, int hardclip3, Resulttype_T resulttype, bool first_read_p, int npaths_mate,
		    int quality_shift, char *sam_read_group_id, bool invertp, bool invert_mate_p,
		    bool use_hardclip_p, bool print_xt_p, char donor1, char donor2, char acceptor2, char acceptor1,
		    double donor_prob, double acceptor_prob, bool circularp) {
  unsigned int flag = 0U;
  int nmismatches_refdiff = 0, nmismatches_bothdiff = 0, querylength;
  int sensedir;
  bool sensep;
  char *genomicfwd_refdiff, *genomicfwd_bothdiff, *genomicdir_refdiff, *genomicdir_bothdiff;
  int substring_start, substring_length;
  int hardclip_low, hardclip_high;
  /* int mate_hardclip_low, mate_hardclip_high; */
  int transloc_hardclip_low, transloc_hardclip_high;
  bool plusp, printp;

  querylength = Shortread_fulllength(queryseq);
  plusp = Substring_plusp(acceptor);

  if (circularp == true) {
    /* clipdir should be +1 */
    if (1 || plusp == true) {
      hardclip_low = hardclip5;
      hardclip_high = hardclip3;
      /* mate_hardclip_low = 0; */
      /* mate_hardclip_high = 0; */
    } else {
      hardclip_low = hardclip3;
      hardclip_high = hardclip5;
      /* mate_hardclip_low = 0; */
      /* mate_hardclip_high = 0; */
    }
  } else {
    if (first_read_p == true) {
      if (clipdir >= 0) {
	hardclip_low = 0;
	hardclip_high = hardclip5;
      } else {
	hardclip_low = hardclip5;
	hardclip_high = 0;
      }
      /* mate_hardclip_low = hardclip3; */
      /* mate_hardclip_high = 0; */
      /* fprintf(stderr,"first read: hardclip_low = %d, hardclip_high = %d\n",hardclip_low,hardclip_high); */
    } else {
      if (clipdir >= 0) {
	hardclip_low = hardclip3;
	hardclip_high = 0;
      } else {
	hardclip_low = 0;
	hardclip_high = hardclip3;
      }
      /* mate_hardclip_low = 0; */
      /* mate_hardclip_high = hardclip5; */
      /* fprintf(stderr,"second read: hardclip_low = %d, hardclip_high = %d\n",hardclip_low,hardclip_high); */
    }
  }


  /* 1. QNAME */
  if (acc2 == NULL) {
    fprintf(fp,"%s",acc1);
  } else {
    fprintf(fp,"%s,%s",acc1,acc2);
  }

  /* 2. FLAG */
  flag = SAM_compute_flag(plusp,mate,resulttype,first_read_p,
			  pathnum,npaths,npaths_mate,absmq_score,first_absmq,
			  invertp,invert_mate_p);
  fprintf(fp,"\t%u",flag);

  /* 3. RNAME: chr */
  /* 4. POS: chrpos */
  print_chromosomal_pos(fp,Substring_chrnum(acceptor),chrpos,chromosome_iit);


  /* 5. MAPQ: Mapping quality */
  fprintf(fp,"\t%d",mapq_score);

  /* 6. CIGAR */
  fprintf(fp,"\t");
  /* sensedir for chimera must be SENSE_FORWARD or SENSE_ANTI, not SENSE_NULL */
  sensedir = Substring_chimera_sensedir(acceptor);
  sensep = Substring_chimera_sensep(acceptor);

  if (use_hardclip_p == true) {
    if (sensep != plusp) {
      transloc_hardclip_low = 0;
      if (plusp == true) {
	/* sensep false */
	transloc_hardclip_high = querylength - Substring_queryend(acceptor);

      } else {
	/* sensep true */
	transloc_hardclip_high = Substring_querystart(acceptor);
      }

    } else { /*  sensep == Substring_plusp(acceptor) */
      transloc_hardclip_high = 0;
      if (plusp == true) {
	transloc_hardclip_low = Substring_querystart(acceptor);

      } else {
	transloc_hardclip_low = querylength - Substring_queryend(acceptor);
      }
    }

    if (transloc_hardclip_low > hardclip_low) {
      hardclip_low = transloc_hardclip_low;
    }
    if (transloc_hardclip_high > hardclip_high) {
      hardclip_high = transloc_hardclip_high;
    }
  }

  if (sensep != plusp) {
    if (plusp == true) {
      /* sensep false */
      assert(Substring_chimera_pos(acceptor) == Substring_queryend(acceptor));
      print_cigar(fp,/*type*/'S',Substring_querystart(acceptor),
		  /*querypos*/0,querylength,hardclip_low,hardclip_high,
		  /*plusp*/true,/*lastp*/false);
      print_cigar(fp,/*type*/'M',Substring_match_length(acceptor),
		  /*querypos*/Substring_querystart(acceptor),querylength,
		  hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);
      print_cigar(fp,/*type*/use_hardclip_p ? 'H' : 'S',querylength - Substring_queryend(acceptor),
		  /*querypos*/Substring_queryend(acceptor),querylength,
		  hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);

    } else {
      /* sensep true */
      assert(Substring_chimera_pos(acceptor) == Substring_querystart(acceptor));
      print_cigar(fp,/*type*/'S',querylength - Substring_queryend(acceptor),
		  /*querypos*/querylength,querylength,
		  hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
      print_cigar(fp,/*type*/'M',Substring_match_length(acceptor),
		  /*querypos*/Substring_queryend(acceptor),querylength,
		  hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
      print_cigar(fp,/*type*/use_hardclip_p ? 'H' : 'S',Substring_querystart(acceptor),
		  /*querypos*/Substring_querystart(acceptor),querylength,hardclip_low,hardclip_high,
		  /*plusp*/false,/*lastp*/true);
    }

  } else { /*  sensep == Substring_plusp(acceptor) */
    if (plusp == true) {
      assert(Substring_chimera_pos(acceptor) == Substring_querystart(acceptor));
      print_cigar(fp,/*type*/use_hardclip_p ? 'H' : 'S',Substring_querystart(acceptor),
		  /*querypos*/0,querylength,hardclip_low,hardclip_high,
		  /*plusp*/true,/*lastp*/false);
      print_cigar(fp,/*type*/'M',Substring_match_length(acceptor),
		  /*querypos*/Substring_querystart(acceptor),querylength,
		  hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);
      print_cigar(fp,/*type*/'S',querylength - Substring_queryend(acceptor),
		  /*querypos*/Substring_queryend(acceptor),querylength,
		  hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);

    } else {
      assert(Substring_chimera_pos(acceptor) == Substring_queryend(acceptor));
      print_cigar(fp,/*type*/use_hardclip_p ? 'H' : 'S',querylength - Substring_queryend(acceptor),
		  /*querypos*/querylength,querylength,
		  hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
      print_cigar(fp,/*type*/'M',Substring_match_length(acceptor),
		  /*querypos*/Substring_queryend(acceptor),querylength,
		  hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
      print_cigar(fp,/*type*/'S',Substring_querystart(acceptor),
		  /*querypos*/Substring_querystart(acceptor),querylength,hardclip_low,hardclip_high,
		  /*plusp*/false,/*lastp*/true);
    }
  }


  /* 7. MRNM: Mate chr */
  /* 8. MPOS: Mate chrpos */
  /* For anchor_chrnum, previously used Stage3end_chrnum(this), but this is 0 */
  print_mate_chromosomal_pos(fp,Stage3end_chrnum(mate),Stage3end_effective_chrnum(mate),mate_chrpos,
			     /*anchor_chrnum*/Substring_chrnum(acceptor),chrpos,chromosome_iit);


  /* 9. ISIZE: Insert size */
  if (resulttype == CONCORDANT_UNIQ || resulttype == CONCORDANT_TRANSLOC || resulttype == CONCORDANT_MULT) {
    if (plusp == invertp) {
      fprintf(fp,"\t%d",-pairedlength);
    } else {
      fprintf(fp,"\t%d",pairedlength);
    }
  } else if (mate_chrpos == 0) {
    fprintf(fp,"\t%d",pairedlength);
  } else if (concordant_chrpos < mate_chrpos) {
    fprintf(fp,"\t%d",pairedlength);
  } else if (concordant_chrpos > mate_chrpos) {
    fprintf(fp,"\t%d",-pairedlength);
  } else if (first_read_p == true) {
    fprintf(fp,"\t%d",pairedlength);
  } else {
    fprintf(fp,"\t%d",-pairedlength);
  }


  /* 10. SEQ: queryseq and 11. QUAL: quality scores */
  /* Queryseq has already been inverted, so just measure plusp relative to its current state */
  fprintf(fp,"\t");
  if (plusp == true) {
    Shortread_print_chopped(fp,queryseq,hardclip_low,hardclip_high);
    fprintf(fp,"\t");
    Shortread_print_quality(fp,queryseq,hardclip_low,hardclip_high,
			    quality_shift,/*show_chopped_p*/false);
  } else {
    Shortread_print_chopped_revcomp(fp,queryseq,hardclip_low,hardclip_high);
    fprintf(fp,"\t");
    Shortread_print_quality_revcomp(fp,queryseq,hardclip_low,hardclip_high,
				    quality_shift,/*show_chopped_p*/false);
  }


  /* 12. TAGS: RG */
  if (sam_read_group_id != NULL) {
    fprintf(fp,"\t");
    fprintf(fp,"RG:Z:%s",sam_read_group_id);
  }

  /* 12. TAGS: XB */
  Shortread_print_barcode(fp,queryseq);

  /* 12. TAGS: XP */
  Shortread_print_chop(fp,queryseq,invertp);

  /* 12. TAGS: MD */
  fprintf(fp,"\t");
  fprintf(fp,"MD:Z:");
  printp = false;

  substring_start = Substring_querystart(acceptor);
  substring_length = Substring_match_length(acceptor);

  if (use_hardclip_p == false) {
    genomicdir_refdiff = Substring_genomic_refdiff(acceptor);
    genomicdir_bothdiff = Substring_genomic_bothdiff(acceptor);
    if (plusp == true) {
      print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,fp,/*matchlength*/0,
		      &(genomicdir_refdiff[substring_start]),&(genomicdir_bothdiff[substring_start]),
		      substring_length,/*querypos*/substring_start,querylength,
		      hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
    } else if (genomicdir_bothdiff == genomicdir_refdiff) {
      genomicfwd_refdiff = (char *) CALLOC(querylength+1,sizeof(char));
      make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring_start]),substring_length);
      print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
		      fp,/*matchlength*/0,genomicfwd_refdiff,/*genomicfwd_bothdiff*/genomicfwd_refdiff,
		      substring_length,/*querypos*/substring_start,querylength,
		      hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
      FREE(genomicfwd_refdiff);
    } else {
      genomicfwd_refdiff = (char *) CALLOC(querylength+1,sizeof(char));
      make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring_start]),substring_length);
      genomicfwd_bothdiff = (char *) CALLOC(querylength+1,sizeof(char));
      make_complement_buffered(genomicfwd_bothdiff,&(genomicdir_bothdiff[substring_start]),substring_length);
      print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
		      fp,/*matchlength*/0,genomicfwd_refdiff,genomicfwd_bothdiff,
		      substring_length,/*querypos*/substring_start,querylength,
		      hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
      FREE(genomicfwd_bothdiff);
      FREE(genomicfwd_refdiff);
    }

  } else if (sensep == false) {
    if (plusp == true) {
      genomicfwd_refdiff = Substring_genomic_refdiff(acceptor);
      genomicfwd_bothdiff = Substring_genomic_bothdiff(acceptor);
      print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,fp,/*matchlength*/0,
		      &(genomicfwd_refdiff[substring_start]),&(genomicfwd_bothdiff[substring_start]),
		      substring_length,/*querypos*/substring_start,querylength,
		      hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
    } else {
      genomicdir_refdiff = Substring_genomic_refdiff(acceptor);
      genomicdir_bothdiff = Substring_genomic_bothdiff(acceptor);
      if (genomicdir_bothdiff == genomicdir_refdiff) {
	genomicfwd_refdiff = (char *) CALLOC(substring_length+1,sizeof(char));
	make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring_start]),substring_length);
	print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
			fp,/*matchlength*/0,genomicfwd_refdiff,/*genomicfwd_bothdiff*/genomicfwd_refdiff,
			substring_length,/*querypos*/substring_start,querylength,
			hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
	FREE(genomicfwd_refdiff);
      } else {
	genomicfwd_refdiff = (char *) CALLOC(substring_length+1,sizeof(char));
	make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring_start]),substring_length);
	genomicfwd_bothdiff = (char *) CALLOC(substring_length+1,sizeof(char));
	make_complement_buffered(genomicfwd_bothdiff,&(genomicdir_bothdiff[substring_start]),substring_length);
	print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
			fp,/*matchlength*/0,genomicfwd_refdiff,genomicfwd_bothdiff,
			substring_length,/*querypos*/substring_start,querylength,
			hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
	FREE(genomicfwd_bothdiff);
	FREE(genomicfwd_refdiff);
      }

    }

  } else {			/* sensep true */
    if (plusp == true) {
      genomicfwd_refdiff = Substring_genomic_refdiff(acceptor);
      genomicfwd_bothdiff = Substring_genomic_bothdiff(acceptor);
      print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,fp,/*matchlength*/0,
		      &(genomicfwd_refdiff[substring_start]),&(genomicfwd_bothdiff[substring_start]),
		      substring_length,/*querypos*/substring_start,querylength,
		      hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
    } else {
      genomicdir_refdiff = Substring_genomic_refdiff(acceptor);
      genomicdir_bothdiff = Substring_genomic_bothdiff(acceptor);
      if (genomicdir_bothdiff == genomicdir_refdiff) {
	genomicfwd_refdiff = (char *) CALLOC(substring_length+1,sizeof(char));
	make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring_start]),substring_length);
	print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
			fp,/*matchlength*/0,genomicfwd_refdiff,/*genomicfwd_bothdiff*/genomicfwd_refdiff,
			substring_length,/*querypos*/substring_start,querylength,
			hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
	FREE(genomicfwd_refdiff);
      } else {
	genomicfwd_refdiff = (char *) CALLOC(substring_length+1,sizeof(char));
	make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring_start]),substring_length);
	genomicfwd_bothdiff = (char *) CALLOC(substring_length+1,sizeof(char));
	make_complement_buffered(genomicfwd_bothdiff,&(genomicdir_bothdiff[substring_start]),substring_length);
	print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
			fp,/*matchlength*/0,genomicfwd_refdiff,genomicfwd_bothdiff,
			substring_length,/*querypos*/substring_start,querylength,
			hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
	FREE(genomicfwd_bothdiff);
	FREE(genomicfwd_refdiff);
      }
    }
  }
  if (printp == false) {
    fprintf(fp,"0");
  }


  /* 12. TAGS: NH */
  fprintf(fp,"\t");
  fprintf(fp,"NH:i:%d",npaths);

  /* 12. TAGS: HI */
  fprintf(fp,"\t");
  fprintf(fp,"HI:i:%d",pathnum);

  /* 12. TAGS: NM */
  fprintf(fp,"\t");
  /* fprintf(fp,"NM:i:%d",Substring_nmismatches_refdiff(acceptor)); */
  fprintf(fp,"NM:i:%d",nmismatches_refdiff);
  
  if (snps_iit) {
    /* 12. TAGS: XW and XV */
    fprintf(fp,"\t");
    fprintf(fp,"XW:i:%d",nmismatches_bothdiff);
    fprintf(fp,"\t");
    fprintf(fp,"XV:i:%d",nmismatches_refdiff - nmismatches_bothdiff);
  }

  /* 12. TAGS: SM */
  fprintf(fp,"\t");
  fprintf(fp,"SM:i:%d",mapq_score);

  /* 12. TAGS: XQ */
  fprintf(fp,"\t");
  fprintf(fp,"XQ:i:%d",absmq_score);

  /* 12. TAGS: X2 */
  fprintf(fp,"\t");
  fprintf(fp,"X2:i:%d",second_absmq);

  /* 12. TAGS: XO */
  fprintf(fp,"\t");
  fprintf(fp,"XO:Z:%s",abbrev);

  /* 12. TAGS: XS */
  fprintf(fp,"\t");
  /* sensedir for chimera must be SENSE_FORWARD or SENSE_ANTI, not SENSE_NULL */
  if (sensedir == SENSE_FORWARD) {
    if (plusp == true) {
      fprintf(fp,"XS:A:+");
    } else {
      fprintf(fp,"XS:A:-");
    }
  } else if (sensedir == SENSE_ANTI) {
    if (plusp == true) {
      fprintf(fp,"XS:A:-");
    } else {
      fprintf(fp,"XS:A:+");
    }
  } else if (force_xs_direction_p == true) {
    fprintf(fp,"XS:A:+");
  } else {
    fprintf(fp,"XS:A:?");
  }

  /* 12. TAGS: XT */
  if (print_xt_p == true) {
    fprintf(fp,"\t");
    fprintf(fp,"XT:Z:%c%c-%c%c,%.2f,%.2f",donor1,donor2,acceptor2,acceptor1,donor_prob,acceptor_prob);
  }

  /* 12. TAGS: PG */
  if (Stage3end_sarrayp(this) == true) {
    fprintf(fp,"\t");
    fprintf(fp,"PG:Z:A");
  }

  fprintf(fp,"\n");
  return;
}



static void
print_localsplice (FILE *fp, char *abbrev, Stage3end_T this, Stage3end_T mate,
		   char *acc1, char *acc2, int pathnum, int npaths,
		   int absmq_score, int first_absmq, int second_absmq, int mapq_score,
		   Univ_IIT_T chromosome_iit, Shortread_T queryseq, int pairedlength,
		   Chrpos_T chrpos, Chrpos_T mate_chrpos, int clipdir, int hardclip5, int hardclip3,
		   Resulttype_T resulttype, bool first_read_p, int npaths_mate,
		   int quality_shift, char *sam_read_group_id, bool invertp, bool invert_mate_p,
		   bool circularp) {
  unsigned int flag = 0U;
  Substring_T substring1, substring2;
  int nmismatches_refdiff = 0, nmismatches_bothdiff = 0, querylength;
  int sensedir;
  bool sensep;
  char *genomicfwd_refdiff, *genomicfwd_bothdiff, *genomicdir_refdiff, *genomicdir_bothdiff;
  int substring1_start, substring2_start, substring1_length, substring2_length, matchlength;
  int hardclip_low, hardclip_high;
  /* int mate_hardclip_low, mate_hardclip_high; */
  bool plusp, printp;

  querylength = Shortread_fulllength(queryseq);
  plusp = Stage3end_plusp(this);

  if ((sensedir = Stage3end_sensedir(this)) == SENSE_NULL) {
    sensedir = Stage3end_sensedir(mate);
  }
  sensep = (sensedir == SENSE_FORWARD);

  if (circularp == true) {
    /* clipdir should be +1 */
    if (1 || plusp == true) {
      hardclip_low = hardclip5;
      hardclip_high = hardclip3;
      /* mate_hardclip_low = 0; */
      /* mate_hardclip_high = 0; */
    } else {
      hardclip_low = hardclip3;
      hardclip_high = hardclip5;
      /* mate_hardclip_low = 0; */
      /* mate_hardclip_high = 0; */
    }
  } else {
    if (first_read_p == true) {
      if (clipdir >= 0) {
	hardclip_low = 0;
	hardclip_high = hardclip5;
      } else {
	hardclip_low = hardclip5;
	hardclip_high = 0;
      }
      /* mate_hardclip_low = hardclip3; */
      /* mate_hardclip_high = 0; */
      /* fprintf(stderr,"first read: hardclip_low = %d, hardclip_high = %d\n",hardclip_low,hardclip_high); */
    } else {
      if (clipdir >= 0) {
	hardclip_low = hardclip3;
	hardclip_high = 0;
      } else {
	hardclip_low = 0;
	hardclip_high = hardclip3;
      }
      /* mate_hardclip_low = 0; */
      /* mate_hardclip_high = hardclip5; */
      /* fprintf(stderr,"second read: hardclip_low = %d, hardclip_high = %d\n",hardclip_low,hardclip_high); */
    }
  }


  /* 1. QNAME */
  if (acc2 == NULL) {
    fprintf(fp,"%s",acc1);
  } else {
    fprintf(fp,"%s,%s",acc1,acc2);
  }

  /* 2. FLAG */
  flag = SAM_compute_flag(plusp,mate,resulttype,first_read_p,
			  pathnum,npaths,npaths_mate,absmq_score,first_absmq,
			  invertp,invert_mate_p);
  fprintf(fp,"\t%u",flag);

  /* 3. RNAME: chr */
  /* 4. POS: chrpos */
  print_chromosomal_pos(fp,Stage3end_chrnum(this),chrpos,chromosome_iit);


  /* 5. MAPQ: Mapping quality */
  fprintf(fp,"\t%d",mapq_score);

  /* 6. CIGAR */
  fprintf(fp,"\t");
  if (sensep == plusp) {
    substring1 = /* donor */ Stage3end_substring_donor(this);
    substring2 = /* acceptor */ Stage3end_substring_acceptor(this);
  } else {
    substring1 = /* acceptor */ Stage3end_substring_acceptor(this);
    substring2 = /* donor */ Stage3end_substring_donor(this);
  }

  if (plusp == true) {
    print_cigar(fp,/*type*/'S',Substring_querystart(substring1),
		/*querypos*/0,querylength,hardclip_low,hardclip_high,
		/*plusp*/true,/*lastp*/false);
    print_cigar(fp,/*type*/'M',Substring_match_length(substring1),
		/*querypos*/Substring_querystart(substring1),querylength,
		hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);
    if (hardclip_low < Substring_queryend(substring1) &&
	querylength - hardclip_high > Substring_querystart(substring2)) {
      debug1(printf("\ncase 1: hardclip_low %d < queryend(substring1) %d && querylength %d - hardclip_high %d > querystart(substring2) %d\n",
		    hardclip_low,Substring_queryend(substring1),querylength,hardclip_high,Substring_querystart(substring2)));
      fprintf(fp,"%uN",Stage3end_distance(this));
    }
    print_cigar(fp,/*type*/'M',Substring_match_length(substring2),
		/*querypos*/Substring_querystart(substring2),querylength,
		hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);
    print_cigar(fp,/*type*/'S',querylength - Substring_queryend(substring2),
		/*querypos*/Substring_queryend(substring2),querylength,
		hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);

  } else {
    print_cigar(fp,/*type*/'S',querylength - Substring_queryend(substring1),
		/*querypos*/querylength,querylength,
		hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
    print_cigar(fp,/*type*/'M',Substring_match_length(substring1),
		/*querypos*/Substring_queryend(substring1),querylength,
		hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
    if (querylength - hardclip_low > Substring_queryend(substring2) &&
	hardclip_high < Substring_querystart(substring1)) {
      debug1(printf("\ncase 2: querylength %d - hardclip_low %d > queryend(substring2) %d && hardclip_high %d < querystart(substring1) %d\n",
		    querylength,hardclip_low,Substring_queryend(substring2),hardclip_high,Substring_querystart(substring1)));
      fprintf(fp,"%uN",Stage3end_distance(this));
    }
    print_cigar(fp,/*type*/'M',Substring_match_length(substring2),
		/*querypos*/Substring_querystart(substring1),querylength,
		hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
    print_cigar(fp,/*type*/'S',Substring_querystart(substring2),
		/*querypos*/Substring_querystart(substring2),querylength,hardclip_low,hardclip_high,
		/*plusp*/false,/*lastp*/true);
  }


  /* 7. MRNM: Mate chr */
  /* 8. MPOS: Mate chrpos */
  print_mate_chromosomal_pos(fp,Stage3end_chrnum(mate),Stage3end_effective_chrnum(mate),mate_chrpos,
			     Stage3end_chrnum(this),chrpos,chromosome_iit);


  /* 9. ISIZE: Insert size */
  if (resulttype == CONCORDANT_UNIQ || resulttype == CONCORDANT_TRANSLOC || resulttype == CONCORDANT_MULT) {
    if (plusp == invertp) {
      fprintf(fp,"\t%d",-pairedlength);
    } else {
      fprintf(fp,"\t%d",pairedlength);
    }
  } else if (mate_chrpos == 0) {
    fprintf(fp,"\t%d",pairedlength);
  } else if (chrpos < mate_chrpos) {
    fprintf(fp,"\t%d",pairedlength);
  } else if (chrpos > mate_chrpos) {
    fprintf(fp,"\t%d",-pairedlength);
  } else if (first_read_p == true) {
    fprintf(fp,"\t%d",pairedlength);
  } else {
    fprintf(fp,"\t%d",-pairedlength);
  }


  /* 10. SEQ: queryseq and 11. QUAL: quality scores */
  /* Queryseq has already been inverted, so just measure plusp relative to its current state */
  fprintf(fp,"\t");
  if (plusp == true) {
    Shortread_print_chopped(fp,queryseq,hardclip_low,hardclip_high);
    fprintf(fp,"\t");
    Shortread_print_quality(fp,queryseq,hardclip_low,hardclip_high,
			   quality_shift,/*show_chopped_p*/false);
  } else {
    Shortread_print_chopped_revcomp(fp,queryseq,hardclip_low,hardclip_high);
    fprintf(fp,"\t");
    Shortread_print_quality_revcomp(fp,queryseq,hardclip_low,hardclip_high,
				   quality_shift,/*show_chopped_p*/false);
  } 

  /* 12. TAGS: RG */
  if (sam_read_group_id != NULL) {
    fprintf(fp,"\t");
    fprintf(fp,"RG:Z:%s",sam_read_group_id);
  }

  /* 12. TAGS: XB */
  Shortread_print_barcode(fp,queryseq);

  /* 12. TAGS: XP */
  Shortread_print_chop(fp,queryseq,invertp);

  /* 12. TAGS: MD */
  fprintf(fp,"\t");
  fprintf(fp,"MD:Z:");
  printp = false;

  substring1_start = Substring_querystart(substring1);
  substring1_length = Substring_match_length(substring1);
  substring2_start = Substring_querystart(substring2);
  substring2_length = Substring_match_length(substring2);

  if (plusp == true) {
    genomicfwd_refdiff = Substring_genomic_refdiff(substring1);
    genomicfwd_bothdiff = Substring_genomic_bothdiff(substring1);
    matchlength = print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,fp,/*matchlength*/0,
				  &(genomicfwd_refdiff[substring1_start]),&(genomicfwd_bothdiff[substring1_start]),
				  substring1_length,/*querypos*/substring1_start,querylength,
				  hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);
    
#if 0
    /* Intron: Gets skipped in MD string */
    fprintf(fp,"^");
#endif

    genomicfwd_refdiff = Substring_genomic_refdiff(substring2);
    genomicfwd_bothdiff = Substring_genomic_bothdiff(substring2);
    print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,fp,matchlength,
		    &(genomicfwd_refdiff[substring2_start]),&(genomicfwd_bothdiff[substring2_start]),
		    substring2_length,/*querypos*/substring2_start,querylength,
		    hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);

  } else {
    genomicdir_refdiff = Substring_genomic_refdiff(substring1);
    genomicdir_bothdiff = Substring_genomic_bothdiff(substring1);
    if (genomicdir_bothdiff == genomicdir_refdiff) {
      genomicfwd_refdiff = (char *) CALLOC(substring1_length+1,sizeof(char));
      make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring1_start]),substring1_length);
      matchlength = print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
				    fp,/*matchlength*/0,genomicfwd_refdiff,/*genomicfwd_bothdiff*/genomicfwd_refdiff,
				    substring1_length,/*querypos*/substring1_start,querylength,
				    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
      FREE(genomicfwd_refdiff);
    } else {
      genomicfwd_refdiff = (char *) CALLOC(substring1_length+1,sizeof(char));
      make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring1_start]),substring1_length);
      genomicfwd_bothdiff = (char *) CALLOC(substring1_length+1,sizeof(char));
      make_complement_buffered(genomicfwd_bothdiff,&(genomicdir_bothdiff[substring1_start]),substring1_length);
      matchlength = print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
				    fp,/*matchlength*/0,genomicfwd_refdiff,genomicfwd_bothdiff,
				    substring1_length,/*querypos*/substring1_start,querylength,
				    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
      FREE(genomicfwd_bothdiff);
      FREE(genomicfwd_refdiff);
    }

#if 0
    /* Intron: Gets skipped in MD string */
    fprintf(fp,"^");
#endif

    genomicdir_refdiff = Substring_genomic_refdiff(substring2);
    genomicdir_bothdiff = Substring_genomic_bothdiff(substring2);
    if (genomicdir_bothdiff == genomicdir_refdiff) {
      genomicfwd_refdiff = (char *) CALLOC(substring2_length+1,sizeof(char));
      make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring2_start]),substring2_length);
      print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
		      fp,matchlength,genomicfwd_refdiff,/*genomicfwd_bothdiff*/genomicfwd_refdiff,
		      substring2_length,/*querypos*/substring2_start,querylength,
		      hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
      FREE(genomicfwd_refdiff);
    } else {
      genomicfwd_refdiff = (char *) CALLOC(substring2_length+1,sizeof(char));
      make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring2_start]),substring2_length);
      genomicfwd_bothdiff = (char *) CALLOC(substring2_length+1,sizeof(char));
      make_complement_buffered(genomicfwd_bothdiff,&(genomicdir_bothdiff[substring2_start]),substring2_length);
      print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
		      fp,matchlength,genomicfwd_refdiff,genomicfwd_bothdiff,
		      substring2_length,/*querypos*/substring2_start,querylength,
		      hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
      FREE(genomicfwd_bothdiff);
      FREE(genomicfwd_refdiff);
    }
  }
  if (printp == false) {
    fprintf(fp,"0");
  }


  /* 12. TAGS: NH */
  fprintf(fp,"\t");
  fprintf(fp,"NH:i:%d",npaths);

  /* 12. TAGS: HI */
  fprintf(fp,"\t");
  fprintf(fp,"HI:i:%d",pathnum);

  /* 12. TAGS: NM */
  fprintf(fp,"\t");
  /* fprintf(fp,"NM:i:%d",Stage3end_nmismatches_refdiff(this)); */
  fprintf(fp,"NM:i:%d",nmismatches_refdiff);
  
  if (snps_iit) {
    /* 12. TAGS: XW and XV */
    fprintf(fp,"\t");
    fprintf(fp,"XW:i:%d",nmismatches_bothdiff);
    fprintf(fp,"\t");
    fprintf(fp,"XV:i:%d",nmismatches_refdiff - nmismatches_bothdiff);
  }

  /* 12. TAGS: SM */
  fprintf(fp,"\t");
  fprintf(fp,"SM:i:%d",mapq_score);

  /* 12. TAGS: XQ */
  fprintf(fp,"\t");
  fprintf(fp,"XQ:i:%d",absmq_score);

  /* 12. TAGS: X2 */
  fprintf(fp,"\t");
  fprintf(fp,"X2:i:%d",second_absmq);

  /* 12. TAGS: XO */
  fprintf(fp,"\t");
  fprintf(fp,"XO:Z:%s",abbrev);

  /* 12. TAGS: XS */
  fprintf(fp,"\t");
  if (sensedir == SENSE_FORWARD) {
    if (plusp == true) {
      fprintf(fp,"XS:A:+");
    } else {
      fprintf(fp,"XS:A:-");
    }
  } else if (sensedir == SENSE_ANTI) {
    if (plusp == true) {
      fprintf(fp,"XS:A:-");
    } else {
      fprintf(fp,"XS:A:+");
    }
  } else if (force_xs_direction_p == true) {
    fprintf(fp,"XS:A:+");
  } else {
    fprintf(fp,"XS:A:?");
  }

  /* 12. TAGS: PG */
  if (Stage3end_sarrayp(this) == true) {
    fprintf(fp,"\t");
    fprintf(fp,"PG:Z:A");
  }

  fprintf(fp,"\n");
  return;
}


static void
print_shortexon (FILE *fp, char *abbrev, Stage3end_T shortexon, Stage3end_T mate,
		 char *acc1, char *acc2, int pathnum, int npaths,
		 int absmq_score, int first_absmq, int second_absmq, int mapq_score,
		 Univ_IIT_T chromosome_iit, Shortread_T queryseq, int pairedlength,
		 Chrpos_T chrpos, Chrpos_T mate_chrpos, int clipdir, int hardclip5, int hardclip3,
		 Resulttype_T resulttype, bool first_read_p, int npaths_mate,
		 int quality_shift, char *sam_read_group_id, bool invertp, bool invert_mate_p,
		 bool circularp) {
  unsigned int flag = 0U;
  /* substring1 is low coordinate on genome, substring2 is high */
  Substring_T substring1, substring2, substringM;
  Chrpos_T distance1, distance2;
  int nmismatches_refdiff = 0, nmismatches_bothdiff = 0, querylength;
  int sensedir;
  bool sensep;
  char *genomicfwd_refdiff, *genomicfwd_bothdiff, *genomicdir_refdiff, *genomicdir_bothdiff;
  int substring1_start, substring2_start, substringM_start,
    substring1_length, substring2_length, substringM_length, matchlength;
  int hardclip_low, hardclip_high;
  /* int mate_hardclip_low, mate_hardclip_high; */
  bool plusp, printp;

  querylength = Shortread_fulllength(queryseq);
  plusp = Stage3end_plusp(shortexon);

  if ((sensedir = Stage3end_sensedir(shortexon)) == SENSE_NULL) {
    sensedir = Stage3end_sensedir(mate);
  }
  sensep = (sensedir == SENSE_FORWARD);

  if (circularp == true) {
    /* clipdir should be +1 */
    if (1 || plusp == true) {
      hardclip_low = hardclip5;
      hardclip_high = hardclip3;
      /* mate_hardclip_low = 0; */
      /* mate_hardclip_high = 0; */
    } else {
      hardclip_low = hardclip3;
      hardclip_high = hardclip5;
      /* mate_hardclip_low = 0; */
      /* mate_hardclip_high = 0; */
    }
  } else {
    if (first_read_p == true) {
      if (clipdir >= 0) {
	hardclip_low = 0;
	hardclip_high = hardclip5;
      } else {
	hardclip_low = hardclip5;
	hardclip_high = 0;
      }
      /* mate_hardclip_low = hardclip3; */
      /* mate_hardclip_high = 0; */
      /* fprintf(stderr,"first read: hardclip_low = %d, hardclip_high = %d\n",hardclip_low,hardclip_high); */
    } else {
      if (clipdir >= 0) {
	hardclip_low = hardclip3;
	hardclip_high = 0;
      } else {
	hardclip_low = 0;
	hardclip_high = hardclip3;
      }
      /* mate_hardclip_low = 0; */
      /* mate_hardclip_high = hardclip5; */
      /* fprintf(stderr,"second read: hardclip_low = %d, hardclip_high = %d\n",hardclip_low,hardclip_high); */
    }
  }


  /* 1. QNAME */
  if (acc2 == NULL) {
    fprintf(fp,"%s",acc1);
  } else {
    fprintf(fp,"%s,%s",acc1,acc2);
  }

  /* 2. FLAG */
  flag = SAM_compute_flag(plusp,mate,resulttype,first_read_p,
			  pathnum,npaths,npaths_mate,absmq_score,first_absmq,
			  invertp,invert_mate_p);
  fprintf(fp,"\t%u",flag);

  /* 3. RNAME: chr */
  /* 4. POS: chrpos */
  print_chromosomal_pos(fp,Stage3end_chrnum(shortexon),chrpos,chromosome_iit);


  /* 5. MAPQ: Mapping quality */
  fprintf(fp,"\t%d",mapq_score);

  /* 6. CIGAR */
  fprintf(fp,"\t");
  substringM = Stage3end_substring1(shortexon);

  if (sensep == plusp) {
    substring1 = /* donor */ Stage3end_substringD(shortexon);
    distance1 = Stage3end_shortexon_acceptor_distance(shortexon);
    distance2 = Stage3end_shortexon_donor_distance(shortexon);
    substring2 = /* acceptor */ Stage3end_substringA(shortexon);
  } else {
    substring1 = /* acceptor */ Stage3end_substringA(shortexon);
    distance1 = Stage3end_shortexon_donor_distance(shortexon);
    distance2 = Stage3end_shortexon_acceptor_distance(shortexon);
    substring2 = /* donor */ Stage3end_substringD(shortexon);
  }

  if (substring1 == NULL) {
    if (plusp == true) {
      print_cigar(fp,/*type*/'S',Substring_querystart(substringM),
		  /*querypos*/0,querylength,hardclip_low,hardclip_high,
		  /*plusp*/true,/*lastp*/false);
    } else {
      print_cigar(fp,/*type*/'S',querylength - Substring_queryend(substringM),
		  /*querypos*/querylength,querylength,
		  hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
    }

  } else if (plusp == true) {
    print_cigar(fp,/*type*/'S',Substring_querystart(substring1),
		/*querypos*/0,querylength,hardclip_low,hardclip_high,
		/*plusp*/true,/*lastp*/false);
    print_cigar(fp,/*type*/'M',Substring_match_length(substring1),
		/*querypos*/Substring_querystart(substring1),querylength,
		hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);
    if (hardclip_low < Substring_queryend(substring1) &&
	querylength - hardclip_high > Substring_querystart(substringM)) {
      debug1(printf("\ncase 3: hardclip_low %d < queryend(substring1) %d && querylength %d - hardclip_high %d > querystart(substringM) %d\n",
		    hardclip_low,Substring_queryend(substring1),querylength,hardclip_high,Substring_querystart(substringM)));
      fprintf(fp,"%uN",distance1);
    }

  } else {
    print_cigar(fp,/*type*/'S',querylength - Substring_queryend(substring1),
		/*querypos*/querylength,querylength,
		hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
    print_cigar(fp,/*type*/'M',Substring_match_length(substring1),
		/*querypos*/Substring_queryend(substring1),querylength,
		hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
    if (querylength - hardclip_low > Substring_queryend(substringM) &&
	hardclip_high < Substring_querystart(substring1)) {
      debug1(printf("\ncase 4: querylength %d - hardclip_low %d > queryend(substringM) %d && hardclip_high %d < querystart(substring1) %d\n",
		    querylength,hardclip_low,Substring_queryend(substringM),hardclip_high,Substring_querystart(substring1)));
      fprintf(fp,"%uN",distance1);
    }
  }

  if (plusp == true) {
    print_cigar(fp,/*type*/'M',Substring_match_length(substringM),
		/*querypos*/Substring_querystart(substringM),querylength,
		hardclip_low,hardclip_high,plusp,/*lastp*/false);
  } else {
    print_cigar(fp,/*type*/'M',Substring_match_length(substringM),
		/*querypos*/Substring_queryend(substringM),querylength,
		hardclip_low,hardclip_high,plusp,/*lastp*/false);
  }

  if (substring2 == NULL) {
    if (plusp == true) {
      print_cigar(fp,/*type*/'S',querylength - Substring_queryend(substringM),
		  /*querypos*/Substring_queryend(substringM),querylength,
		  hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
    } else {
      print_cigar(fp,/*type*/'S',Substring_querystart(substringM),
		  /*querypos*/Substring_querystart(substringM),querylength,
		  hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
    }

  } else if (plusp == true) {
    if (hardclip_low < Substring_queryend(substringM) &&
	querylength - hardclip_high > Substring_querystart(substring2)) {
      debug1(printf("\ncase 5: hardclip_low %d < queryend(substringM) %d && querylength %d - hardclip_high %d > querystart(substring2) %d\n",
		    hardclip_low,Substring_queryend(substringM),querylength,hardclip_high,Substring_querystart(substring2)));
      fprintf(fp,"%uN",distance2);
    }
    print_cigar(fp,/*type*/'M',Substring_match_length(substring2),
		/*querypos*/Substring_querystart(substring2),querylength,
		hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);
    print_cigar(fp,/*type*/'S',querylength - Substring_queryend(substring2),
		/*querypos*/Substring_queryend(substring2),querylength,
		hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);

  } else {
    if (querylength - hardclip_low > Substring_queryend(substring2) &&
	hardclip_high < Substring_querystart(substringM)) {
      debug1(printf("\ncase 6: querylength %d - hardclip_low %d > queryend(substring2) %d && hardclip_high %d < querystart(substringM) %d\n",
		    querylength,hardclip_low,Substring_queryend(substring2),querylength,Substring_querystart(substringM)));
      fprintf(fp,"%uN",distance2);
    }
    print_cigar(fp,/*type*/'M',Substring_match_length(substring2),
		/*querypos*/Substring_queryend(substring2),querylength,
		hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
    print_cigar(fp,/*type*/'S',Substring_querystart(substring2),
		/*querypos*/Substring_querystart(substring2),querylength,hardclip_low,hardclip_high,
		/*plusp*/false,/*lastp*/true);
  }


  /* 7. MRNM: Mate chr */
  /* 8. MPOS: Mate chrpos */
  print_mate_chromosomal_pos(fp,Stage3end_chrnum(mate),Stage3end_effective_chrnum(mate),mate_chrpos,
			     Stage3end_chrnum(shortexon),chrpos,chromosome_iit);


  /* 9. ISIZE: Insert size */
  if (resulttype == CONCORDANT_UNIQ || resulttype == CONCORDANT_TRANSLOC || resulttype == CONCORDANT_MULT) {
    if (plusp == invertp) {
      fprintf(fp,"\t%d",-pairedlength);
    } else {
      fprintf(fp,"\t%d",pairedlength);
    }
  } else if (mate_chrpos == 0) {
    fprintf(fp,"\t%d",pairedlength);
  } else if (chrpos < mate_chrpos) {
    fprintf(fp,"\t%d",pairedlength);
  } else if (chrpos > mate_chrpos) {
    fprintf(fp,"\t%d",-pairedlength);
  } else if (first_read_p == true) {
    fprintf(fp,"\t%d",pairedlength);
  } else {
    fprintf(fp,"\t%d",-pairedlength);
  }


  /* 10. SEQ: queryseq and 11. QUAL: quality scores */
  /* Queryseq has already been inverted, so just measure plusp relative to its current state */
  fprintf(fp,"\t");
  if (plusp == true) {
    Shortread_print_chopped(fp,queryseq,hardclip_low,hardclip_high);
    fprintf(fp,"\t");
    Shortread_print_quality(fp,queryseq,hardclip_low,hardclip_high,
			   quality_shift,/*show_chopped_p*/false);
  } else {
    Shortread_print_chopped_revcomp(fp,queryseq,hardclip_low,hardclip_high);
    fprintf(fp,"\t");
    Shortread_print_quality_revcomp(fp,queryseq,hardclip_low,hardclip_high,
				   quality_shift,/*show_chopped_p*/false);
  } 

  /* 12. TAGS: RG */
  if (sam_read_group_id != NULL) {
    fprintf(fp,"\t");
    fprintf(fp,"RG:Z:%s",sam_read_group_id);
  }

  /* 12. TAGS: XB */
  Shortread_print_barcode(fp,queryseq);

  /* 12. TAGS: XP */
  Shortread_print_chop(fp,queryseq,invertp);

  /* 12. TAGS: MD */
  fprintf(fp,"\t");
  fprintf(fp,"MD:Z:");
  printp = false;

  substringM_start = Substring_querystart(substringM);
  substringM_length = Substring_match_length(substringM);

  if (substring1 == NULL) {
    substring1_start = 0;
    substring1_length = 0;
  } else {
    substring1_start = Substring_querystart(substring1);
    substring1_length = Substring_match_length(substring1);
  }
  if (substring2 == NULL) {
    substring2_start = 0;
    substring2_length = 0;
  } else {
    substring2_start = Substring_querystart(substring2);
    substring2_length = Substring_match_length(substring2);
  }

  if (plusp == true) {

    if (substring1 == NULL) {
      matchlength = 0;
    } else {
      genomicfwd_refdiff = Substring_genomic_refdiff(substring1);
      genomicfwd_bothdiff = Substring_genomic_bothdiff(substring1);
      matchlength = print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,fp,/*matchlength*/0,
				    &(genomicfwd_refdiff[substring1_start]),&(genomicfwd_bothdiff[substring1_start]),
				    substring1_length,/*querypos*/substring1_start,querylength,
				    hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);
    }

#if 0
    /* Intron 1: Gets skipped in MD string */
    fprintf(fp,"^");
#endif

    genomicfwd_refdiff = Substring_genomic_refdiff(substringM);
    genomicfwd_bothdiff = Substring_genomic_bothdiff(substringM);
    matchlength = print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,fp,matchlength,
				  &(genomicfwd_refdiff[substringM_start]),&(genomicfwd_bothdiff[substringM_start]),
				  substringM_length,/*querypos*/substringM_start,querylength,
				  hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);

#if 0
    /* Intron 2: Gets skipped in MD string */
    fprintf(fp,"^");
#endif

    if (substring2 == NULL) {
      /* Equivalent: if (matchlength > 0) fprintf(fp,"%d",matchlength); */
      print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
		      fp,matchlength,/*genomicfwd_refdiff*/NULL,/*genomicfwd_bothdiff*/NULL,
		      /*substring2_length*/0,/*querypos*/0,querylength,
		      hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
    } else {
      genomicfwd_refdiff = Substring_genomic_refdiff(substring2);
      genomicfwd_bothdiff = Substring_genomic_bothdiff(substring2);
      print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,fp,matchlength,
		      &(genomicfwd_refdiff[substring2_start]),&(genomicfwd_bothdiff[substring2_start]),
		      substring2_length,/*querypos*/substring2_start,querylength,
		      hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
    }

  } else {

    if (substring1 == NULL) {
      matchlength = 0;
    } else {
      genomicdir_refdiff = Substring_genomic_refdiff(substring1);
      genomicdir_bothdiff = Substring_genomic_bothdiff(substring1);
      if (genomicdir_bothdiff == genomicdir_refdiff) {
	genomicfwd_refdiff = (char *) CALLOC(substring1_length+1,sizeof(char));
	make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring1_start]),substring1_length);
	matchlength = print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
				      fp,/*matchlength*/0,genomicfwd_refdiff,/*genomicfwd_bothdiff*/genomicfwd_refdiff,
				      substring1_length,/*querypos*/substring1_start,querylength,
				      hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
	FREE(genomicfwd_refdiff);
      } else {
	genomicfwd_refdiff = (char *) CALLOC(substring1_length+1,sizeof(char));
	make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring1_start]),substring1_length);
	genomicfwd_bothdiff = (char *) CALLOC(substring1_length+1,sizeof(char));
	make_complement_buffered(genomicfwd_bothdiff,&(genomicdir_bothdiff[substring1_start]),substring1_length);
	matchlength = print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
				      fp,/*matchlength*/0,genomicfwd_refdiff,genomicfwd_bothdiff,
				      substring1_length,/*querypos*/substring1_start,querylength,
				      hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
	FREE(genomicfwd_bothdiff);
	FREE(genomicfwd_refdiff);
      }
    }

#if 0
    /* Intron 1: Gets skipped in MD string */
    fprintf(fp,"^");
#endif

    genomicdir_refdiff = Substring_genomic_refdiff(substringM);
    genomicdir_bothdiff = Substring_genomic_bothdiff(substringM);
    if (genomicdir_bothdiff == genomicdir_refdiff) {
      genomicfwd_refdiff = (char *) CALLOC(substringM_length+1,sizeof(char));
      make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substringM_start]),substringM_length);
      matchlength = print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
				    fp,matchlength,genomicfwd_refdiff,/*genomicfwd_bothdiff*/genomicfwd_refdiff,
				    substringM_length,/*querypos*/substringM_start,querylength,
				    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
      FREE(genomicfwd_refdiff);
    } else {
      genomicfwd_refdiff = (char *) CALLOC(substringM_length+1,sizeof(char));
      make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substringM_start]),substringM_length);
      genomicfwd_bothdiff = (char *) CALLOC(substringM_length+1,sizeof(char));
      make_complement_buffered(genomicfwd_bothdiff,&(genomicdir_bothdiff[substringM_start]),substringM_length);
      matchlength = print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
				    fp,matchlength,genomicfwd_refdiff,genomicfwd_bothdiff,
				    substringM_length,/*querypos*/substringM_start,querylength,
				    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
      FREE(genomicfwd_bothdiff);
      FREE(genomicfwd_refdiff);
    }

#if 0
    /* Intron 2: Not sure how to handle in MD string */
    fprintf(fp,"^");
#endif

    if (substring2 == NULL) {
      /* Equivalent: if (matchlength > 0) fprintf(fp,"%d",matchlength); */
      print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
		      fp,matchlength,/*genomicfwd_refdiff*/NULL,/*genomicfwd_bothdiff*/NULL,
		      /*substring2_length*/0,/*querypos*/0,querylength,
		      hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
    } else {
      genomicdir_refdiff = Substring_genomic_refdiff(substring2);
      genomicdir_bothdiff = Substring_genomic_bothdiff(substring2);
      if (genomicdir_bothdiff == genomicdir_refdiff) {
	genomicfwd_refdiff = (char *) CALLOC(substring2_length+1,sizeof(char));
	make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring2_start]),substring2_length);
	print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
			fp,matchlength,genomicfwd_refdiff,/*genomicfwd_bothdiff*/genomicfwd_refdiff,
			substring2_length,/*querypos*/substring2_start,querylength,
			hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
	FREE(genomicfwd_refdiff);
      } else {
	genomicfwd_refdiff = (char *) CALLOC(substring2_length+1,sizeof(char));
	make_complement_buffered(genomicfwd_refdiff,&(genomicdir_refdiff[substring2_start]),substring2_length);
	genomicfwd_bothdiff = (char *) CALLOC(substring2_length+1,sizeof(char));
	make_complement_buffered(genomicfwd_bothdiff,&(genomicdir_bothdiff[substring2_start]),substring2_length);
	print_md_string(&printp,&nmismatches_refdiff,&nmismatches_bothdiff,
			fp,matchlength,genomicfwd_refdiff,genomicfwd_bothdiff,
			substring2_length,/*querypos*/substring2_start,querylength,
			hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
	FREE(genomicfwd_bothdiff);
	FREE(genomicfwd_refdiff);
      }
    }
  }
  if (printp == false) {
    fprintf(fp,"0");
  }


  /* 12. TAGS: NH */
  fprintf(fp,"\t");
  fprintf(fp,"NH:i:%d",npaths);

  /* 12. TAGS: HI */
  fprintf(fp,"\t");
  fprintf(fp,"HI:i:%d",pathnum);

  /* 12. TAGS: NM */
  fprintf(fp,"\t");
  /* fprintf(fp,"NM:i:%d",Stage3end_nmismatches_refdiff(shortexon)); */
  fprintf(fp,"NM:i:%d",nmismatches_refdiff);
  
  if (snps_iit) {
    /* 12. TAGS: XW and XV */
    fprintf(fp,"\t");
    fprintf(fp,"XW:i:%d",nmismatches_bothdiff);
    fprintf(fp,"\t");
    fprintf(fp,"XV:i:%d",nmismatches_refdiff - nmismatches_bothdiff);
  }

  /* 12. TAGS: SM */
  fprintf(fp,"\t");
  fprintf(fp,"SM:i:%d",mapq_score);

  /* 12. TAGS: XQ */
  fprintf(fp,"\t");
  fprintf(fp,"XQ:i:%d",absmq_score);

  /* 12. TAGS: X2 */
  fprintf(fp,"\t");
  fprintf(fp,"X2:i:%d",second_absmq);

  /* 12. TAGS: XO */
  fprintf(fp,"\t");
  fprintf(fp,"XO:Z:%s",abbrev);

  /* 12. TAGS: XS */
  fprintf(fp,"\t");
  if (sensedir == SENSE_FORWARD) {
    if (plusp == true) {
      fprintf(fp,"XS:A:+");
    } else {
      fprintf(fp,"XS:A:-");
    }
  } else if (sensedir == SENSE_ANTI) {
    if (plusp == true) {
      fprintf(fp,"XS:A:-");
    } else {
      fprintf(fp,"XS:A:+");
    }
  } else if (force_xs_direction_p == true) {
    fprintf(fp,"XS:A:+");
  } else {
    fprintf(fp,"XS:A:?");
  }

  /* 12. TAGS: PG */
  if (Stage3end_sarrayp(shortexon) == true) {
    fprintf(fp,"\t");
    fprintf(fp,"PG:Z:A");
  }

  fprintf(fp,"\n");
  return;
}



/* Distant splicing, including scramble, inversion, translocation */
static void
print_exon_exon (FILE *fp, char *abbrev, Stage3end_T this, Stage3end_T mate,
		 char *acc1, char *acc2, int pathnum, int npaths,
		 int absmq_score, int first_absmq, int second_absmq, int mapq_score,
		 Univ_IIT_T chromosome_iit, Shortread_T queryseq, int pairedlength,
		 Chrpos_T mate_chrpos, int clipdir, int hardclip5, int hardclip3,
		 Resulttype_T resulttype, bool first_read_p, int npaths_mate,
		 int quality_shift, char *sam_read_group_id, bool invertp, bool invert_mate_p) {
  Chrpos_T donor_chrpos, acceptor_chrpos, concordant_chrpos;
  Substring_T donor, acceptor;
  int hardclip_low, hardclip_high;
  char donor1, donor2, acceptor2, acceptor1;
  double donor_prob, acceptor_prob;
  int circularpos, querylength;

  donor = Stage3end_substring_donor(this);
  acceptor = Stage3end_substring_acceptor(this);

#if 0
  if (first_read_p == true) {
    hardclip_low = 0;
    hardclip_high = hardclip5;
  } else {
    hardclip_low = hardclip3;
    hardclip_high = 0;
  }
#else
  /* Shouldn't have any overlap on a distant splice */
  hardclip_low = hardclip_high = 0;
#endif

  querylength = Shortread_fulllength(queryseq);
  donor_chrpos = SAM_compute_chrpos(&hardclip_low,&hardclip_high,
				    clipdir,hardclip5,hardclip3,first_read_p,
				    this,/*substring_low*/donor,querylength);
  acceptor_chrpos = SAM_compute_chrpos(&hardclip_low,&hardclip_high,
				       clipdir,hardclip5,hardclip3,first_read_p,
				       this,/*substring_low*/acceptor,querylength);
  if (Stage3end_substring_low(this) == donor) {
    concordant_chrpos = donor_chrpos;
  } else if (Stage3end_substring_low(this) == acceptor) {
    concordant_chrpos = acceptor_chrpos;
  } else {
    fprintf(stderr,"Stage3end_substring_low %p is neither donor %p or acceptor %p\n",
	    Stage3end_substring_low(this),donor,acceptor);
    concordant_chrpos = 0U;
  }

  halfdonor_dinucleotide(&donor1,&donor2,donor);
  halfacceptor_dinucleotide(&acceptor2,&acceptor1,acceptor);
  donor_prob = Substring_chimera_prob(donor);
  acceptor_prob = Substring_chimera_prob(acceptor);

  if (Stage3end_sensedir(this) == SENSE_FORWARD) {
    /* NEEDS WORK: Need to decide whether to split halfdonor or halfacceptor */
    /* Not sure if circular chromosomes should participate in distant splicing anyway */
    if (0 && (circularpos = Stage3end_circularpos(this)) > 0) {
      print_halfdonor(fp,abbrev,donor,this,mate,acc1,acc2,pathnum,npaths,
		      absmq_score,first_absmq,second_absmq,mapq_score,
		      chromosome_iit,queryseq,pairedlength,
		      concordant_chrpos,donor_chrpos,mate_chrpos,
		      /*clipdir*/+1,/*hardclip5*/0,/*hardclip3*/querylength-circularpos,resulttype,first_read_p,
		      npaths_mate,quality_shift,sam_read_group_id,
		      invertp,invert_mate_p,/*use_hardclip_p*/true,/*print_xt_p*/true,
		      donor1,donor2,acceptor2,acceptor1,donor_prob,acceptor_prob,
		      /*circularp*/true);
      print_halfdonor(fp,abbrev,donor,this,mate,acc1,acc2,pathnum,npaths,
		      absmq_score,first_absmq,second_absmq,mapq_score,
		      chromosome_iit,queryseq,pairedlength,
		      /*concordant_chrpos*/1,/*donor_chrpos*/1,mate_chrpos,
		      /*clipdir*/+1,/*hardclip5*/circularpos,/*hardclip3*/0,resulttype,first_read_p,
		      npaths_mate,quality_shift,sam_read_group_id,
		      invertp,invert_mate_p,/*use_hardclip_p*/true,/*print_xt_p*/true,
		      donor1,donor2,acceptor2,acceptor1,donor_prob,acceptor_prob,
		      /*circularp*/true);
    } else {
      print_halfdonor(fp,abbrev,donor,this,mate,acc1,acc2,pathnum,npaths,
		      absmq_score,first_absmq,second_absmq,mapq_score,
		      chromosome_iit,queryseq,pairedlength,
		      concordant_chrpos,donor_chrpos,mate_chrpos,
		      clipdir,hardclip5,hardclip3,resulttype,first_read_p,
		      npaths_mate,quality_shift,sam_read_group_id,
		      invertp,invert_mate_p,/*use_hardclip_p*/true,/*print_xt_p*/true,
		      donor1,donor2,acceptor2,acceptor1,donor_prob,acceptor_prob,
		      /*circularp*/false);
    }

    if (0 && (circularpos = Stage3end_circularpos(this)) > 0) {
      print_halfacceptor(fp,abbrev,acceptor,this,mate,acc1,acc2,pathnum,npaths,
			 absmq_score,first_absmq,second_absmq,mapq_score,
			 chromosome_iit,queryseq,pairedlength,
			 concordant_chrpos,acceptor_chrpos,mate_chrpos,
			 /*clipdir*/+1,/*hardclip5*/0,/*hardclip3*/querylength-circularpos,
			 resulttype,first_read_p,npaths_mate,quality_shift,sam_read_group_id,
			 invertp,invert_mate_p,/*use_hardclip_p*/true,/*print_xt_p*/true,
			 donor1,donor2,acceptor2,acceptor1,donor_prob,acceptor_prob,
			 /*circularp*/true);
      print_halfacceptor(fp,abbrev,acceptor,this,mate,acc1,acc2,pathnum,npaths,
			 absmq_score,first_absmq,second_absmq,mapq_score,
			 chromosome_iit,queryseq,pairedlength,
			 /*concordant_chrpos*/1,/*acceptor_chrpos*/1,mate_chrpos,
			 /*clipdir*/+1,/*hardclip5*/circularpos,/*hardclip3*/0,
			 resulttype,first_read_p,npaths_mate,quality_shift,sam_read_group_id,
			 invertp,invert_mate_p,/*use_hardclip_p*/true,/*print_xt_p*/true,
			 donor1,donor2,acceptor2,acceptor1,donor_prob,acceptor_prob,
			 /*circularp*/true);
    } else {
      print_halfacceptor(fp,abbrev,acceptor,this,mate,acc1,acc2,pathnum,npaths,
			 absmq_score,first_absmq,second_absmq,mapq_score,
			 chromosome_iit,queryseq,pairedlength,
			 concordant_chrpos,acceptor_chrpos,mate_chrpos,
			 clipdir,hardclip5,hardclip3,resulttype,first_read_p,
			 npaths_mate,quality_shift,sam_read_group_id,
			 invertp,invert_mate_p,/*use_hardclip_p*/true,/*print_xt_p*/true,
			 donor1,donor2,acceptor2,acceptor1,donor_prob,acceptor_prob,
			 /*circularp*/false);
    }

  } else if (Stage3end_sensedir(this) == SENSE_ANTI) {
    if (0 && (circularpos = Stage3end_circularpos(this)) > 0) {
      print_halfacceptor(fp,abbrev,acceptor,this,mate,acc1,acc2,pathnum,npaths,
			 absmq_score,first_absmq,second_absmq,mapq_score,
			 chromosome_iit,queryseq,pairedlength,
			 concordant_chrpos,acceptor_chrpos,mate_chrpos,
			 /*clipdir*/+1,/*hardclip5*/0,/*hardclip3*/querylength-circularpos,
			 resulttype,first_read_p,npaths_mate,quality_shift,sam_read_group_id,
			 invertp,invert_mate_p,/*use_hardclip_p*/true,/*print_xt_p*/true,
			 donor1,donor2,acceptor2,acceptor1,donor_prob,acceptor_prob,
			 /*circularp*/true);
      print_halfacceptor(fp,abbrev,acceptor,this,mate,acc1,acc2,pathnum,npaths,
			 absmq_score,first_absmq,second_absmq,mapq_score,
			 chromosome_iit,queryseq,pairedlength,
			 /*concordant_chrpos*/1,/*acceptor_chrpos*/1,mate_chrpos,
			 /*clipdir*/+1,/*hardclip5*/circularpos,/*hardclip3*/0,
			 resulttype,first_read_p,npaths_mate,quality_shift,sam_read_group_id,
			 invertp,invert_mate_p,/*use_hardclip_p*/true,/*print_xt_p*/true,
			 donor1,donor2,acceptor2,acceptor1,donor_prob,acceptor_prob,
			 /*circularp*/true);
    } else {
      print_halfacceptor(fp,abbrev,acceptor,this,mate,acc1,acc2,pathnum,npaths,
			 absmq_score,first_absmq,second_absmq,mapq_score,
			 chromosome_iit,queryseq,pairedlength,
			 concordant_chrpos,acceptor_chrpos,mate_chrpos,
			 clipdir,hardclip5,hardclip3,resulttype,first_read_p,
			 npaths_mate,quality_shift,sam_read_group_id,
			 invertp,invert_mate_p,/*use_hardclip_p*/true,/*print_xt_p*/true,
			 donor1,donor2,acceptor2,acceptor1,donor_prob,acceptor_prob,
			 /*circularp*/false);
    }

    if (0 && (circularpos = Stage3end_circularpos(this)) > 0) {
      print_halfdonor(fp,abbrev,donor,this,mate,acc1,acc2,pathnum,npaths,
		      absmq_score,first_absmq,second_absmq,mapq_score,
		      chromosome_iit,queryseq,pairedlength,
		      concordant_chrpos,donor_chrpos,mate_chrpos,
		      /*clipdir*/+1,/*hardclip5*/0,/*hardclip3*/querylength-circularpos,
		      resulttype,first_read_p,npaths_mate,quality_shift,sam_read_group_id,
		      invertp,invert_mate_p,/*use_hardclip_p*/true,/*print_xt_p*/true,
		      donor1,donor2,acceptor2,acceptor1,donor_prob,acceptor_prob,
		      /*circularp*/true);
      print_halfdonor(fp,abbrev,donor,this,mate,acc1,acc2,pathnum,npaths,
		      absmq_score,first_absmq,second_absmq,mapq_score,
		      chromosome_iit,queryseq,pairedlength,
		      /*concordant_chrpos*/1,/*donor_chrpos*/1,mate_chrpos,
		      /*clipdir*/+1,/*hardclip5*/circularpos,/*hardclip3*/0,
		      resulttype,first_read_p,npaths_mate,quality_shift,sam_read_group_id,
		      invertp,invert_mate_p,/*use_hardclip_p*/true,/*print_xt_p*/true,
		      donor1,donor2,acceptor2,acceptor1,donor_prob,acceptor_prob,
		      /*circularp*/true);
    } else {
      print_halfdonor(fp,abbrev,donor,this,mate,acc1,acc2,pathnum,npaths,
		      absmq_score,first_absmq,second_absmq,mapq_score,
		      chromosome_iit,queryseq,pairedlength,
		      concordant_chrpos,donor_chrpos,mate_chrpos,
		      clipdir,hardclip5,hardclip3,resulttype,first_read_p,
		      npaths_mate,quality_shift,sam_read_group_id,
		      invertp,invert_mate_p,/*use_hardclip_p*/true,/*print_xt_p*/true,
		      donor1,donor2,acceptor2,acceptor1,donor_prob,acceptor_prob,
		      /*circularp*/false);
    }

  } else {
    abort();
  }

  return;
}



void
SAM_print (FILE *fp, char *abbrev, Stage3end_T this, Stage3end_T mate,
	   char *acc1, char *acc2, int pathnum, int npaths,
	   int absmq_score, int first_absmq, int second_absmq, int mapq_score, Univ_IIT_T chromosome_iit, Shortread_T queryseq,
	   Shortread_T queryseq_mate, int pairedlength, Chrpos_T chrpos, Chrpos_T mate_chrpos,
	   int clipdir, int hardclip5, int hardclip3, Resulttype_T resulttype, bool first_read_p,
	   int npaths_mate, int quality_shift, char *sam_read_group_id, bool invertp, bool invert_mate_p,
	   bool merge_samechr_p) {
  Hittype_T hittype;
  Substring_T donor, acceptor;
  bool sensep, normalp;
  unsigned int flag;
  int ignore = 0;
  int circularpos, querylength;


  hittype = Stage3end_hittype(this);
  /* printf("hittype %s, chrpos %u\n",Stage3end_hittype_string(this),chrpos); */
  if (npaths == 0) {		/* was chrpos == 0, but we can actually align to chrpos 0 */
    SAM_print_nomapping(fp,abbrev,queryseq,mate,acc1,acc2,chromosome_iit,resulttype,first_read_p,
			npaths_mate,mate_chrpos,quality_shift,
			sam_read_group_id,invertp,invert_mate_p);

  } else if (hittype == EXACT || hittype == SUB || hittype == TERMINAL) {
    if ((circularpos = Stage3end_circularpos(this)) > 0) {
      querylength = Shortread_fulllength(queryseq);
#ifdef CHECK_ASSERTIONS
      if (Stage3end_plusp(this) == true) {
	assert(chrpos-Stage3end_trim_left(this)+circularpos-Stage3end_chrlength(this) == 1);
      } else {
	assert(chrpos-Stage3end_trim_right(this)+circularpos-Stage3end_chrlength(this) == 1);
      }
#endif
      print_single(fp,abbrev,this,mate,acc1,acc2,pathnum,npaths,
		   absmq_score,first_absmq,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
		   chrpos,mate_chrpos,/*clipdir*/+1,/*hardclip5*/0,/*hardclip3*/querylength-circularpos,resulttype,first_read_p,
		   npaths_mate,quality_shift,sam_read_group_id,
		   invertp,invert_mate_p,/*circularp*/true);
      print_single(fp,abbrev,this,mate,acc1,acc2,pathnum,npaths,
		   absmq_score,first_absmq,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
		   /*chrpos*/1,mate_chrpos,/*clipdir*/+1,/*hardclip5*/circularpos,/*hardclip3*/0,resulttype,first_read_p,
		   npaths_mate,quality_shift,sam_read_group_id,
		   invertp,invert_mate_p,/*circularp*/true);
    } else {
      print_single(fp,abbrev,this,mate,acc1,acc2,pathnum,npaths,
		   absmq_score,first_absmq,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
		   chrpos,mate_chrpos,clipdir,hardclip5,hardclip3,resulttype,first_read_p,
		   npaths_mate,quality_shift,sam_read_group_id,
		   invertp,invert_mate_p,/*circularp*/false);
    }

  } else if (hittype == INSERTION) {
    if ((circularpos = Stage3end_circularpos(this)) > 0) {
      querylength = Shortread_fulllength(queryseq);
      print_insertion(fp,abbrev,this,mate,acc1,acc2,pathnum,npaths,
		      absmq_score,first_absmq,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
		      chrpos,mate_chrpos,/*clipdir*/+1,/*hardclip5*/0,/*hardclip3*/querylength-circularpos,
		      resulttype,first_read_p,npaths_mate,quality_shift,sam_read_group_id,
		      invertp,invert_mate_p,/*circularp*/true);
      print_insertion(fp,abbrev,this,mate,acc1,acc2,pathnum,npaths,
		      absmq_score,first_absmq,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
		      /*chrpos*/1,mate_chrpos,/*clipdir*/+1,/*hardclip5*/circularpos,/*hardclip3*/0,
		      resulttype,first_read_p,npaths_mate,quality_shift,sam_read_group_id,
		      invertp,invert_mate_p,/*circularp*/true);
    } else {
      print_insertion(fp,abbrev,this,mate,acc1,acc2,pathnum,npaths,
		      absmq_score,first_absmq,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
		      chrpos,mate_chrpos,clipdir,hardclip5,hardclip3,resulttype,first_read_p,
		      npaths_mate,quality_shift,sam_read_group_id,
		      invertp,invert_mate_p,/*circularp*/false);
    }

  } else if (hittype == DELETION) {
    if ((circularpos = Stage3end_circularpos(this)) > 0) {
      querylength = Shortread_fulllength(queryseq);
      print_deletion(fp,abbrev,this,mate,acc1,acc2,pathnum,npaths,
		     absmq_score,first_absmq,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
		     chrpos,mate_chrpos,/*clipdir*/+1,/*hardclip5*/0,/*hardclip3*/querylength-circularpos,
		     resulttype,first_read_p,npaths_mate,quality_shift,sam_read_group_id,
		     invertp,invert_mate_p,/*circularp*/true);
      print_deletion(fp,abbrev,this,mate,acc1,acc2,pathnum,npaths,
		     absmq_score,first_absmq,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
		     /*chrpos*/1,mate_chrpos,/*clipdir*/+1,/*hardclip5*/circularpos,/*hardclip3*/0,
		     resulttype,first_read_p,npaths_mate,quality_shift,sam_read_group_id,
		     invertp,invert_mate_p,/*circularp*/true);
    } else {
      print_deletion(fp,abbrev,this,mate,acc1,acc2,pathnum,npaths,
		     absmq_score,first_absmq,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
		     chrpos,mate_chrpos,clipdir,hardclip5,hardclip3,resulttype,first_read_p,
		     npaths_mate,quality_shift,sam_read_group_id,
		     invertp,invert_mate_p,/*circularp*/false);
    }

  } else if (hittype == HALFSPLICE_DONOR) {
    if ((circularpos = Stage3end_circularpos(this)) > 0) {
      querylength = Shortread_fulllength(queryseq);
      print_halfdonor(fp,abbrev,Stage3end_substring_donor(this),this,mate,acc1,acc2,pathnum,npaths,
		      absmq_score,first_absmq,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
		      /*concordant_chrpos*/chrpos,chrpos,mate_chrpos,
		      /*clipdir*/+1,/*hardclip5*/0,/*hardclip3*/querylength-circularpos,
		      resulttype,first_read_p,npaths_mate,quality_shift,sam_read_group_id,
		      invertp,invert_mate_p,/*use_hardclip_p*/false,/*print_xt_p*/false,
		      /*donor1*/'X',/*donor2*/'X',/*acceptor2*/'X',/*acceptor1*/'X',
		      /*donor_prob*/0.0,/*acceptor_prob*/0.0,/*circularp*/true);
      print_halfdonor(fp,abbrev,Stage3end_substring_donor(this),this,mate,acc1,acc2,pathnum,npaths,
		      absmq_score,first_absmq,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
		      /*concordant_chrpos*/1,/*chrpos*/1,mate_chrpos,
		      /*clipdir*/+1,/*hardclip5*/circularpos,/*hardclip3*/0,
		      resulttype,first_read_p,npaths_mate,quality_shift,sam_read_group_id,
		      invertp,invert_mate_p,/*use_hardclip_p*/false,/*print_xt_p*/false,
		      /*donor1*/'X',/*donor2*/'X',/*acceptor2*/'X',/*acceptor1*/'X',
		      /*donor_prob*/0.0,/*acceptor_prob*/0.0,/*circularp*/true);
    } else {
      print_halfdonor(fp,abbrev,Stage3end_substring_donor(this),this,mate,acc1,acc2,pathnum,npaths,
		      absmq_score,first_absmq,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
		      /*concordant_chrpos*/chrpos,chrpos,mate_chrpos,
		      clipdir,hardclip5,hardclip3,resulttype,first_read_p,
		      npaths_mate,quality_shift,sam_read_group_id,
		      invertp,invert_mate_p,/*use_hardclip_p*/false,/*print_xt_p*/false,
		      /*donor1*/'X',/*donor2*/'X',/*acceptor2*/'X',/*acceptor1*/'X',
		      /*donor_prob*/0.0,/*acceptor_prob*/0.0,/*circularp*/false);
    }

  } else if (hittype == HALFSPLICE_ACCEPTOR) {
    if ((circularpos = Stage3end_circularpos(this)) > 0) {
      querylength = Shortread_fulllength(queryseq);
      print_halfacceptor(fp,abbrev,Stage3end_substring_acceptor(this),this,mate,acc1,acc2,pathnum,npaths,
			 absmq_score,first_absmq,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
			 /*concordant_chrpos*/chrpos,chrpos,mate_chrpos,
			 /*clipdir*/+1,/*hardclip5*/0,/*hardclip3*/querylength-circularpos,
			 resulttype,first_read_p,npaths_mate,quality_shift,sam_read_group_id,
			 invertp,invert_mate_p,/*use_hardclip_p*/false,/*print_xt_p*/false,
			 /*donor1*/'X',/*donor2*/'X',/*acceptor2*/'X',/*acceptor1*/'X',
			 /*donor_prob*/0.0,/*acceptor_prob*/0.0,/*circularp*/true);
      print_halfacceptor(fp,abbrev,Stage3end_substring_acceptor(this),this,mate,acc1,acc2,pathnum,npaths,
			 absmq_score,first_absmq,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
			 /*concordant_chrpos*/1,/*chrpos*/1,mate_chrpos,
			 /*clipdir*/+1,/*hardclip5*/circularpos,/*hardclip3*/0,
			 resulttype,first_read_p,npaths_mate,quality_shift,sam_read_group_id,
			 invertp,invert_mate_p,/*use_hardclip_p*/false,/*print_xt_p*/false,
			 /*donor1*/'X',/*donor2*/'X',/*acceptor2*/'X',/*acceptor1*/'X',
			 /*donor_prob*/0.0,/*acceptor_prob*/0.0,/*circularp*/true);
    } else {
      print_halfacceptor(fp,abbrev,Stage3end_substring_acceptor(this),this,mate,acc1,acc2,pathnum,npaths,
			 absmq_score,first_absmq,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
			 /*concordant_chrpos*/chrpos,chrpos,mate_chrpos,
			 clipdir,hardclip5,hardclip3,resulttype,first_read_p,
			 npaths_mate,quality_shift,sam_read_group_id,
			 invertp,invert_mate_p,/*use_hardclip_p*/false,/*print_xt_p*/false,
			 /*donor1*/'X',/*donor2*/'X',/*acceptor2*/'X',/*acceptor1*/'X',
			 /*donor_prob*/0.0,/*acceptor_prob*/0.0,/*circularp*/false);
    }

  } else if (hittype == SPLICE || hittype == SAMECHR_SPLICE || hittype == TRANSLOC_SPLICE) {
    /* Follows print_splice_distance() in substring.c */
    donor = Stage3end_substring_donor(this);
    acceptor = Stage3end_substring_acceptor(this);

    if (donor == NULL || acceptor == NULL) {
      abort();
    } else if (hittype == TRANSLOC_SPLICE || (hittype == SAMECHR_SPLICE && merge_samechr_p == false)) {
      /* Stage3end_chrnum(this) == 0 || Stage3end_distance(this) == 0U */
      /* distant splice */
      print_exon_exon(fp,abbrev,this,mate,acc1,acc2,pathnum,npaths,
		      absmq_score,first_absmq,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
		      mate_chrpos,clipdir,hardclip5,hardclip3,resulttype,first_read_p,
		      npaths_mate,quality_shift,sam_read_group_id,
		      invertp,invert_mate_p);
    } else {
      normalp = true;
      sensep = (Stage3end_sensedir(this) == SENSE_FORWARD);

      if (Substring_plusp(donor) != Substring_plusp(acceptor)) {
	/* inversion */
	normalp = false;
      } else if (Substring_plusp(donor) == true) {
	if (sensep == true) {
	  if (Substring_genomicstart(acceptor) < Substring_genomicstart(donor)) {
	    /* scramble */
	    normalp = false;
	  }
	} else {
	  if (Substring_genomicstart(donor) < Substring_genomicstart(acceptor)) {
	    /* scramble */
	    normalp = false;
	  }
	}
      } else {
	if (sensep == true) {
	  if (Substring_genomicstart(donor) < Substring_genomicstart(acceptor)) {
	    /* scramble */
	    normalp = false;
	  }
	} else {
	  if (Substring_genomicstart(acceptor) < Substring_genomicstart(donor)) {
	    /* scramble */
	    normalp = false;
	  }
	}
      }
      if (normalp == true) {
	if ((circularpos = Stage3end_circularpos(this)) > 0) {
	  querylength = Shortread_fulllength(queryseq);
	  print_localsplice(fp,abbrev,this,mate,acc1,acc2,pathnum,npaths,
			    absmq_score,first_absmq,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
			    chrpos,mate_chrpos,/*clipdir*/+1,/*hardclip5*/0,/*hardclip3*/querylength-circularpos,
			    resulttype,first_read_p,npaths_mate,quality_shift,sam_read_group_id,
			    invertp,invert_mate_p,/*circularp*/true);
	  print_localsplice(fp,abbrev,this,mate,acc1,acc2,pathnum,npaths,
			    absmq_score,first_absmq,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
			    /*chrpos*/1,mate_chrpos,/*clipdir*/+1,/*hardclip5*/circularpos,/*hardclip3*/0,
			    resulttype,first_read_p,npaths_mate,quality_shift,sam_read_group_id,
			    invertp,invert_mate_p,/*circularp*/true);
	} else {
	  print_localsplice(fp,abbrev,this,mate,acc1,acc2,pathnum,npaths,
			    absmq_score,first_absmq,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
			    chrpos,mate_chrpos,clipdir,hardclip5,hardclip3,
			    resulttype,first_read_p,npaths_mate,quality_shift,sam_read_group_id,
			    invertp,invert_mate_p,/*circularp*/false);
	}
      } else {
	print_exon_exon(fp,abbrev,this,mate,acc1,acc2,pathnum,npaths,
			absmq_score,first_absmq,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
			mate_chrpos,clipdir,hardclip5,hardclip3,
			resulttype,first_read_p,npaths_mate,quality_shift,sam_read_group_id,
			invertp,invert_mate_p);
      }
    }
      
  } else if (hittype == ONE_THIRD_SHORTEXON || hittype == TWO_THIRDS_SHORTEXON || hittype == SHORTEXON) {
    if ((circularpos = Stage3end_circularpos(this)) > 0) {
      querylength = Shortread_fulllength(queryseq);
      print_shortexon(fp,abbrev,this,mate,acc1,acc2,pathnum,npaths,
		      absmq_score,first_absmq,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
		      chrpos,mate_chrpos,/*clipdir*/+1,/*hardclip5*/0,/*hardclip3*/querylength-circularpos,
		      resulttype,first_read_p,npaths_mate,quality_shift,sam_read_group_id,
		      invertp,invert_mate_p,/*circularp*/true);
      print_shortexon(fp,abbrev,this,mate,acc1,acc2,pathnum,npaths,
		      absmq_score,first_absmq,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
		      /*chrpos*/1,mate_chrpos,/*clipdir*/+1,/*hardclip5*/circularpos,/*hardclip3*/0,
		      resulttype,first_read_p,npaths_mate,quality_shift,sam_read_group_id,
		      invertp,invert_mate_p,/*circularp*/true);
    } else {
      print_shortexon(fp,abbrev,this,mate,acc1,acc2,pathnum,npaths,
		      absmq_score,first_absmq,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
		      chrpos,mate_chrpos,clipdir,hardclip5,hardclip3,
		      resulttype,first_read_p,npaths_mate,quality_shift,sam_read_group_id,
		      invertp,invert_mate_p,/*circularp*/false);
    }

  } else if (hittype == GMAP) {
    /* Note: sam_paired_p must be true because we are calling GMAP only on halfmapping uniq */

    if (mate == NULL) {
      chrpos = SAM_compute_chrpos(/*hardclip_low*/&ignore,/*hardclip_high*/&hardclip5,
				  clipdir,hardclip5,hardclip3,first_read_p,
				  this,/*substring_low*/NULL,Shortread_fulllength(queryseq));
      mate_chrpos = 0U;
      hardclip3 = 0;

    } else if (first_read_p == true) {
      if (clipdir >= 0) {
	chrpos = SAM_compute_chrpos(/*hardclip_low*/&ignore,/*hardclip_high*/&hardclip5,
				    clipdir,hardclip5,hardclip3,/*first_read_p*/true,
				    this,/*substring_low*/NULL,Shortread_fulllength(queryseq));
	mate_chrpos = SAM_compute_chrpos(/*hardclip_low*/&hardclip3,/*hardclip_high*/&ignore,
					 clipdir,hardclip5,hardclip3,/*first_read_p*/false,
					 mate,Stage3end_substring_low(mate),
					 Shortread_fulllength(queryseq_mate));
      } else {
	chrpos = SAM_compute_chrpos(/*hardclip_low*/&hardclip5,/*hardclip_high*/&ignore,
				    clipdir,hardclip5,hardclip3,/*first_read_p*/true,
				    this,/*substring_low*/NULL,Shortread_fulllength(queryseq));
	mate_chrpos = SAM_compute_chrpos(/*hardclip_low*/&ignore,/*hardclip_high*/&hardclip3,
					 clipdir,hardclip5,hardclip3,/*first_read_p*/false,
					 mate,Stage3end_substring_low(mate),
					 Shortread_fulllength(queryseq_mate));
      }
    } else {
      if (clipdir >= 0) {
	chrpos = SAM_compute_chrpos(/*hardclip_low*/&hardclip3,/*hardclip_high*/&ignore,
				    clipdir,hardclip5,hardclip3,/*first_read_p*/false,
				    this,/*substring_low*/NULL,Shortread_fulllength(queryseq));
	mate_chrpos = SAM_compute_chrpos(/*hardclip_low*/&ignore,/*hardclip_high*/&hardclip5,
					 clipdir,hardclip5,hardclip3,/*first_read_p*/true,
					 mate,Stage3end_substring_low(mate),
					 Shortread_fulllength(queryseq_mate));
      } else {
	chrpos = SAM_compute_chrpos(/*hardclip_low*/&ignore,/*hardclip_high*/&hardclip3,
				    clipdir,hardclip5,hardclip3,/*first_read_p*/false,
				    this,/*substring_low*/NULL,Shortread_fulllength(queryseq));
	mate_chrpos = SAM_compute_chrpos(/*hardclip_low*/&hardclip5,/*hardclip_high*/&ignore,
					 clipdir,hardclip5,hardclip3,/*first_read_p*/true,
					 mate,Stage3end_substring_low(mate),
					 Shortread_fulllength(queryseq_mate));
      }
    }
    assert(ignore == 0);

    flag = SAM_compute_flag(Stage3end_plusp(this),mate,resulttype,first_read_p,
			    pathnum,npaths,npaths_mate,absmq_score,first_absmq,
			    invertp,invert_mate_p);

    if ((circularpos = Stage3end_circularpos(this)) > 0) {
      querylength = Shortread_fulllength(queryseq);
      Pair_print_sam(fp,abbrev,Stage3end_pairarray(this),Stage3end_npairs(this),
		     acc1,acc2,Stage3end_chrnum(this),chromosome_iit,
		     /*usersegment*/(Sequence_T) NULL,
		     Shortread_fullpointer(queryseq),Shortread_quality_string(queryseq),
		     /*clipdir*/+1,/*hardclip5*/0,/*hardclip3*/querylength-circularpos,Shortread_fulllength(queryseq),
		     Stage3end_plusp(this),Stage3end_cdna_direction(this),
		     /*chimera_part*/0,/*chimera*/NULL,quality_shift,first_read_p,
		     pathnum,npaths,absmq_score,first_absmq,second_absmq,chrpos,
		     resulttype,flag,/*pair_mapq_score*/mapq_score,/*end_mapq_score*/mapq_score,
		     Stage3end_chrnum(mate),Stage3end_effective_chrnum(mate),mate_chrpos,
		     /*mate_cdna_direction*/Stage3end_cdna_direction(mate),
		     pairedlength,sam_read_group_id,invertp,/*circularp*/true);
      Pair_print_sam(fp,abbrev,Stage3end_pairarray(this),Stage3end_npairs(this),
		     acc1,acc2,Stage3end_chrnum(this),chromosome_iit,
		     /*usersegment*/(Sequence_T) NULL,
		     Shortread_fullpointer(queryseq),Shortread_quality_string(queryseq),
		     /*clipdir*/+1,/*hardclip5*/circularpos,/*hardclip3*/0,Shortread_fulllength(queryseq),
		     Stage3end_plusp(this),Stage3end_cdna_direction(this),
		     /*chimera_part*/0,/*chimera*/NULL,quality_shift,first_read_p,
		     pathnum,npaths,absmq_score,first_absmq,second_absmq,/*chrpos*/1,
		     resulttype,flag,/*pair_mapq_score*/mapq_score,/*end_mapq_score*/mapq_score,
		     Stage3end_chrnum(mate),Stage3end_effective_chrnum(mate),mate_chrpos,
		     /*mate_cdna_direction*/Stage3end_cdna_direction(mate),
		     pairedlength,sam_read_group_id,invertp,/*circularp*/true);
    } else {
      Pair_print_sam(fp,abbrev,Stage3end_pairarray(this),Stage3end_npairs(this),
		     acc1,acc2,Stage3end_chrnum(this),chromosome_iit,
		     /*usersegment*/(Sequence_T) NULL,
		     Shortread_fullpointer(queryseq),Shortread_quality_string(queryseq),
		     clipdir,hardclip5,hardclip3,Shortread_fulllength(queryseq),
		     Stage3end_plusp(this),Stage3end_cdna_direction(this),
		     /*chimera_part*/0,/*chimera*/NULL,quality_shift,first_read_p,
		     pathnum,npaths,absmq_score,first_absmq,second_absmq,chrpos,
		     resulttype,flag,/*pair_mapq_score*/mapq_score,/*end_mapq_score*/mapq_score,
		     Stage3end_chrnum(mate),Stage3end_effective_chrnum(mate),mate_chrpos,
		     /*mate_cdna_direction*/Stage3end_cdna_direction(mate),
		     pairedlength,sam_read_group_id,invertp,/*circularp*/false);
    }
  } else {
    abort();
  }

  return;
}



void
SAM_print_paired (Result_T result, Resulttype_T resulttype,
		  Univ_IIT_T chromosome_iit, Shortread_T queryseq1, Shortread_T queryseq2,
		  bool invert_first_p, bool invert_second_p,
		  bool nofailsp, bool failsonlyp, bool fails_as_input_p,
		  bool fastq_format_p, bool clip_overlap_p, bool merge_samechr_p,
		  int quality_shift, char *sam_read_group_id,
		  FILE *fp_nomapping_1, FILE *fp_nomapping_2,
		  FILE *fp_unpaired_uniq, FILE *fp_unpaired_circular,
		  FILE *fp_unpaired_transloc, FILE *fp_unpaired_mult,
		  FILE *fp_halfmapping_uniq, FILE *fp_halfmapping_circular,
		  FILE *fp_halfmapping_transloc, FILE *fp_halfmapping_mult,
		  FILE *fp_paired_uniq_circular, FILE *fp_paired_uniq_inv, FILE *fp_paired_uniq_scr,
		  FILE *fp_paired_uniq_long, FILE *fp_paired_mult,
		  FILE *fp_concordant_uniq, FILE *fp_concordant_circular,
		  FILE *fp_concordant_transloc, FILE *fp_concordant_mult) {
  Stage3pair_T *stage3pairarray, stage3pair;
  Stage3end_T *stage3array1, *stage3array2, stage3, mate, hit5, hit3;
  Chrpos_T chrpos5, chrpos3;
  int npaths, npaths1, npaths2, pathnum;
  int first_absmq, second_absmq, first_absmq1, second_absmq1, first_absmq2, second_absmq2;
  int hardclip5 = 0, hardclip3 = 0, ignore = 0, clipdir;
  char *acc1, *acc2;
  Pairtype_T pairtype;
  FILE *fp;
  char *abbrev;

  acc1 = Shortread_accession(queryseq1);
  acc2 = Shortread_accession(queryseq2); /* NULL, unless --allow-pe-name-mismatch is specified */

  if (resulttype == PAIREDEND_NOMAPPING) {
    if (nofailsp == true) {
      /* No output */
      return;
      
    } else if (fails_as_input_p == true) {
      if (fastq_format_p == true) {
	Shortread_print_query_pairedend_fastq(fp_nomapping_1,fp_nomapping_2,queryseq1,queryseq2,
					     invert_first_p,invert_second_p);
      } else {
	Shortread_print_query_pairedend_fasta(fp_nomapping_1,queryseq1,queryseq2,
					     invert_first_p,invert_second_p);
      }

    } else {
      SAM_print_nomapping(fp_nomapping_1,ABBREV_NOMAPPING_1,queryseq1,/*mate*/(Stage3end_T) NULL,
			  acc1,acc2,chromosome_iit,resulttype,
			  /*first_read_p*/true,/*npaths_mate*/0,
			  /*mate_chrpos*/0U,quality_shift,
			  sam_read_group_id,invert_first_p,invert_second_p);
      SAM_print_nomapping(fp_nomapping_1,ABBREV_NOMAPPING_1,queryseq2,/*mate*/(Stage3end_T) NULL,
			  acc1,acc2,chromosome_iit,resulttype,
			  /*first_read_p*/false,/*npaths_mate*/0,
			  /*mate_chrpos*/0U,quality_shift,
			  sam_read_group_id,invert_second_p,invert_first_p);
    }

  } else {
    if (failsonlyp == true) {
      /* Unwanted success: skip */

    } else if (resulttype == CONCORDANT_UNIQ) {
      stage3pairarray = (Stage3pair_T *) Result_array(&npaths,&first_absmq,&second_absmq,result);
      /* Stage3pair_eval(stage3pairarray,npaths,maxpaths_report,queryseq1,queryseq2); */

      stage3pair = stage3pairarray[0];
      hardclip5 = hardclip3 = 0;
      if (clip_overlap_p == false) {
	clipdir = 0;
      } else if (Stage3pair_circularp(stage3pair) == true) {
	/* Don't resolve overlaps on a circular alignment */
	clipdir = 0;
      } else {
	clipdir = Stage3pair_overlap(&hardclip5,&hardclip3,stage3pair);
	debug3(printf("clipdir %d with hardclip5 = %d, hardclip3 = %d\n",clipdir,hardclip5,hardclip3));
      }

      hit5 = Stage3pair_hit5(stage3pair);
      hit3 = Stage3pair_hit3(stage3pair);
      if (clipdir >= 0) {
	debug3(printf("clipping %d on hit5 high and %d on hit3 low\n",hardclip5,hardclip3));
	chrpos5 = SAM_compute_chrpos(/*hardclip_low*/&ignore,/*hardclip_high*/&hardclip5,
				     clipdir,hardclip5,hardclip3,/*first_read_p*/true,hit5,
				     Stage3end_substring_low(hit5),Shortread_fulllength(queryseq1));
	chrpos3 = SAM_compute_chrpos(/*hardclip_low*/&hardclip3,/*hardclip_high*/&ignore,
				     clipdir,hardclip5,hardclip3,/*first_read_p*/false,hit3,
				     Stage3end_substring_low(hit3),Shortread_fulllength(queryseq2));
      } else {
	debug3(printf("clipping %d on hit5 low and %d on hit3 high\n",hardclip5,hardclip3));
	chrpos5 = SAM_compute_chrpos(/*hardclip_low*/&hardclip5,/*hardclip_high*/&ignore,
				     clipdir,hardclip5,hardclip3,/*first_read_p*/true,hit5,
				     Stage3end_substring_low(hit5),Shortread_fulllength(queryseq1));
	chrpos3 = SAM_compute_chrpos(/*hardclip_low*/&ignore,/*hardclip_high*/&hardclip3,
				     clipdir,hardclip5,hardclip3,/*first_read_p*/false,hit3,
				     Stage3end_substring_low(hit3),Shortread_fulllength(queryseq2));
      }

      /* Note: it is possible that hard clipping can changed a concordant_uniq to an unpaired_uniq,
	 or a concordant_transloc to an unpaired_transloc */

      if (Stage3pair_circularp(stage3pair) == true) {
	fp = fp_concordant_circular;
	abbrev = ABBREV_CONCORDANT_CIRCULAR;
      } else {
	fp = fp_concordant_uniq;
	abbrev = ABBREV_CONCORDANT_UNIQ;
      }

      /* print first end */
      SAM_print(fp,abbrev,hit5,/*mate*/hit3,acc1,acc2,/*pathnum*/1,/*npaths*/1,
		Stage3pair_absmq_score(stage3pair),first_absmq,/*second_absmq*/0,
		Stage3pair_mapq_score(stage3pair),chromosome_iit,
		/*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
		Stage3pair_pairlength(stage3pair),/*chrpos*/chrpos5,/*mate_chrpos*/chrpos3,
		clipdir,hardclip5,hardclip3,resulttype,/*first_read_p*/true,/*npaths_mate*/npaths,
		quality_shift,sam_read_group_id,invert_first_p,invert_second_p,
		merge_samechr_p);

      /* print second end */
      SAM_print(fp,abbrev,hit3,/*mate*/hit5,acc1,acc2,/*pathnum*/1,/*npaths*/1,
		Stage3pair_absmq_score(stage3pair),first_absmq,/*second_absmq*/0,
		Stage3pair_mapq_score(stage3pair),chromosome_iit,
		/*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
		Stage3pair_pairlength(stage3pair),/*chrpos*/chrpos3,/*mate_chrpos*/chrpos5,
		clipdir,hardclip5,hardclip3,resulttype,/*first_read_p*/false,/*npaths_mate*/npaths,
		quality_shift,sam_read_group_id,invert_second_p,invert_first_p,
		merge_samechr_p);

    } else if (resulttype == CONCORDANT_TRANSLOC) {
      stage3pairarray = (Stage3pair_T *) Result_array(&npaths,&first_absmq,&second_absmq,result);

      if (quiet_if_excessive_p && npaths > maxpaths_report) {
	/* Print as nomapping, but send to fp_concordant_transloc */
	SAM_print_nomapping(fp_concordant_transloc,ABBREV_CONCORDANT_TRANSLOC,
			    queryseq1,/*mate*/(Stage3end_T) NULL,
			    acc1,acc2,chromosome_iit,resulttype,
			    /*first_read_p*/true,/*npaths_mate*/npaths,
			    /*mate_chrpos*/0U,quality_shift,
			    sam_read_group_id,invert_first_p,invert_second_p);
	SAM_print_nomapping(fp_concordant_transloc,ABBREV_CONCORDANT_TRANSLOC,
			    queryseq2,/*mate*/(Stage3end_T) NULL,
			    acc1,acc2,chromosome_iit,resulttype,
			    /*first_read_p*/false,/*npaths_mate*/npaths,
			    /*mate_chrpos*/0U,quality_shift,
			    sam_read_group_id,invert_second_p,invert_first_p);

      } else {
	/* Stage3pair_eval(stage3pairarray,npaths,maxpaths_report,queryseq1,queryseq2); */
	for (pathnum = 1; pathnum <= npaths && pathnum <= maxpaths_report; pathnum++) {

	  stage3pair = stage3pairarray[pathnum-1];
	  hardclip5 = hardclip3 = 0;
	  if (clip_overlap_p == false) {
	    clipdir = 0;
	  } else if (Stage3pair_circularp(stage3pair) == true) {
	    /* Don't resolve overlaps on a circular alignment */
	    clipdir = 0;
	  } else {
	    clipdir = Stage3pair_overlap(&hardclip5,&hardclip3,stage3pair);
	    debug3(printf("clipdir %d with hardclip5 = %d, hardclip3 = %d\n",clipdir,hardclip5,hardclip3));
	  }

	  hit5 = Stage3pair_hit5(stage3pair);
	  hit3 = Stage3pair_hit3(stage3pair);
	  if (clipdir >= 0) {
	    chrpos5 = SAM_compute_chrpos(/*hardclip_low*/&ignore,/*hardclip_high*/&hardclip5,
					 clipdir,hardclip5,hardclip3,/*first_read_p*/true,hit5,
					 Stage3end_substring_low(hit5),Shortread_fulllength(queryseq1));
	    chrpos3 = SAM_compute_chrpos(/*hardclip_low*/&hardclip3,/*hardclip_high*/&ignore,
					 clipdir,hardclip5,hardclip3,/*first_read_p*/false,hit3,
					 Stage3end_substring_low(hit3),Shortread_fulllength(queryseq2));
	  } else {
	    chrpos5 = SAM_compute_chrpos(/*hardclip_low*/&hardclip5,/*hardclip_high*/&ignore,
					 clipdir,hardclip5,hardclip3,/*first_read_p*/true,hit5,
					 Stage3end_substring_low(hit5),Shortread_fulllength(queryseq1));
	    chrpos3 = SAM_compute_chrpos(/*hardclip_low*/&ignore,/*hardclip_high*/&hardclip3,
					 clipdir,hardclip5,hardclip3,/*first_read_p*/false,hit3,
					 Stage3end_substring_low(hit3),Shortread_fulllength(queryseq2));
	  }

	  /* print first end */
	  SAM_print(fp_concordant_transloc,ABBREV_CONCORDANT_TRANSLOC,
		    hit5,/*mate*/hit3,acc1,acc2,pathnum,npaths,
		    Stage3pair_absmq_score(stage3pair),first_absmq,second_absmq,
		    Stage3pair_mapq_score(stage3pair),chromosome_iit,
		    /*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
		    Stage3pair_pairlength(stage3pair),/*chrpos*/chrpos5,/*mate_chrpos*/chrpos3,
		    clipdir,hardclip5,hardclip3,resulttype,/*first_read_p*/true,/*npaths_mate*/npaths,
		    quality_shift,sam_read_group_id,invert_first_p,invert_second_p,
		    merge_samechr_p);

	  /* print second end */
	  SAM_print(fp_concordant_transloc,ABBREV_CONCORDANT_TRANSLOC,
		    hit3,/*mate*/hit5,acc1,acc2,pathnum,npaths,
		    Stage3pair_absmq_score(stage3pair),first_absmq,second_absmq,
		    Stage3pair_mapq_score(stage3pair),chromosome_iit,
		    /*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
		    Stage3pair_pairlength(stage3pair),/*chrpos*/chrpos3,/*mate_chrpos*/chrpos5,
		    clipdir,hardclip5,hardclip3,resulttype,/*first_read_p*/false,/*npaths_mate*/npaths,
		    quality_shift,sam_read_group_id,invert_second_p,invert_first_p,
		    merge_samechr_p);
	}
      }
    
    } else if (resulttype == CONCORDANT_MULT) {
      stage3pairarray = (Stage3pair_T *) Result_array(&npaths,&first_absmq,&second_absmq,result);

      if (quiet_if_excessive_p && npaths > maxpaths_report) {
	/* Print as nomapping, but send to fp_concordant_mult */
	SAM_print_nomapping(fp_concordant_mult,ABBREV_CONCORDANT_MULT,
			    queryseq1,/*mate*/(Stage3end_T) NULL,
			    acc1,acc2,chromosome_iit,resulttype,
			    /*first_read_p*/true,/*npaths_mate*/npaths,
			    /*mate_chrpos*/0U,quality_shift,
			    sam_read_group_id,invert_first_p,invert_second_p);
	SAM_print_nomapping(fp_concordant_mult,ABBREV_CONCORDANT_MULT,
			    queryseq2,/*mate*/(Stage3end_T) NULL,
			    acc1,acc2,chromosome_iit,resulttype,
			    /*first_read_p*/false,/*npaths_mate*/npaths,
			    /*mate_chrpos*/0U,quality_shift,
			    sam_read_group_id,invert_second_p,invert_first_p);

      } else {
	/* Stage3pair_eval(stage3pairarray,npaths,maxpaths_report,queryseq1,queryseq2); */
	for (pathnum = 1; pathnum <= npaths && pathnum <= maxpaths_report; pathnum++) {

	  stage3pair = stage3pairarray[pathnum-1];
	  hardclip5 = hardclip3 = 0;
	  if (clip_overlap_p == false) {
	    clipdir = 0;
	  } else if (Stage3pair_circularp(stage3pair) == true) {
	    /* Don't resolve overlaps on a circular alignment */
	    clipdir = 0;
	  } else {
	    clipdir = Stage3pair_overlap(&hardclip5,&hardclip3,stage3pair);
	    debug3(printf("clipdir %d with hardclip5 = %d, hardclip3 = %d\n",clipdir,hardclip5,hardclip3));
	  }

	  hit5 = Stage3pair_hit5(stage3pair);
	  hit3 = Stage3pair_hit3(stage3pair);
	  if (clipdir >= 0) {
	    chrpos5 = SAM_compute_chrpos(/*hardclip_low*/&ignore,/*hardclip_high*/&hardclip5,
					 clipdir,hardclip5,hardclip3,/*first_read_p*/true,hit5,
					 Stage3end_substring_low(hit5),Shortread_fulllength(queryseq1));
	    chrpos3 = SAM_compute_chrpos(/*hardclip_low*/&hardclip3,/*hardclip_high*/&ignore,
					 clipdir,hardclip5,hardclip3,/*first_read_p*/false,hit3,
					 Stage3end_substring_low(hit3),Shortread_fulllength(queryseq2));
	  } else {
	    chrpos5 = SAM_compute_chrpos(/*hardclip_low*/&hardclip5,/*hardclip_high*/&ignore,
					 clipdir,hardclip5,hardclip3,/*first_read_p*/true,hit5,
					 Stage3end_substring_low(hit5),Shortread_fulllength(queryseq1));
	    chrpos3 = SAM_compute_chrpos(/*hardclip_low*/&ignore,/*hardclip_high*/&hardclip3,
					 clipdir,hardclip5,hardclip3,/*first_read_p*/false,hit3,
					 Stage3end_substring_low(hit3),Shortread_fulllength(queryseq2));
	  }

	  /* print first end */
	  SAM_print(fp_concordant_mult,ABBREV_CONCORDANT_MULT,
		    hit5,/*mate*/hit3,acc1,acc2,pathnum,npaths,
		    Stage3pair_absmq_score(stage3pair),first_absmq,second_absmq,
		    Stage3pair_mapq_score(stage3pair),chromosome_iit,
		    /*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
		    Stage3pair_pairlength(stage3pair),/*chrpos*/chrpos5,/*mate_chrpos*/chrpos3,
		    clipdir,hardclip5,hardclip3,resulttype,/*first_read_p*/true,/*npaths_mate*/npaths,
		    quality_shift,sam_read_group_id,invert_first_p,invert_second_p,
		    merge_samechr_p);

	  /* print second end */
	  SAM_print(fp_concordant_mult,ABBREV_CONCORDANT_MULT,
		    hit3,/*mate*/hit5,acc1,acc2,pathnum,npaths,
		    Stage3pair_absmq_score(stage3pair),first_absmq,second_absmq,
		    Stage3pair_mapq_score(stage3pair),chromosome_iit,
		    /*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
		    Stage3pair_pairlength(stage3pair),/*chrpos*/chrpos3,/*mate_chrpos*/chrpos5,
		    clipdir,hardclip5,hardclip3,resulttype,/*first_read_p*/false,/*npaths_mate*/npaths,
		    quality_shift,sam_read_group_id,invert_second_p,invert_first_p,
		    merge_samechr_p);
	}
      }

    } else if (resulttype == PAIRED_UNIQ) {
      stage3pairarray = (Stage3pair_T *) Result_array(&npaths,&first_absmq,&second_absmq,result);
      /* Stage3pair_eval(stage3pairarray,npaths,maxpaths_report,queryseq1,queryseq2); */

      stage3pair = stage3pairarray[0];
      if (Stage3pair_circularp(stage3pair) == true) {
	fp = fp_paired_uniq_circular;
	abbrev = ABBREV_PAIRED_UNIQ_CIRCULAR;
      } else if ((pairtype = Stage3pair_pairtype(stage3pair)) == PAIRED_INVERSION) {
	fp = fp_paired_uniq_inv;
	abbrev = ABBREV_PAIRED_UNIQ_INV;
      } else if (pairtype == PAIRED_SCRAMBLE) {
	fp = fp_paired_uniq_scr;
	abbrev = ABBREV_PAIRED_UNIQ_SCR;
      } else if (pairtype == PAIRED_TOOLONG) {
	fp = fp_paired_uniq_long;
	abbrev = ABBREV_PAIRED_UNIQ_LONG;
      } else {
	fprintf(stderr,"Unexpected pairtype %d\n",pairtype);
	abort();
      }

      hardclip5 = hardclip3 = 0;

      hit5 = Stage3pair_hit5(stage3pair);
      hit3 = Stage3pair_hit3(stage3pair);
      chrpos5 = SAM_compute_chrpos(/*hardclip_low*/&ignore,/*hardclip_high*/&hardclip5,
				   /*clipdir*/0,hardclip5,hardclip3,/*first_read_p*/true,hit5,
				   Stage3end_substring_low(hit5),Shortread_fulllength(queryseq1));
      chrpos3 = SAM_compute_chrpos(/*hardclip_low*/&hardclip3,/*hardclip_high*/&ignore,
				   /*clipdir*/0,hardclip5,hardclip3,/*first_read_p*/false,hit3,
				   Stage3end_substring_low(hit3),Shortread_fulllength(queryseq2));

      /* print first end */
      SAM_print(fp,abbrev,hit5,/*mate*/hit3,acc1,acc2,/*pathnum*/1,/*npaths*/1,
		Stage3pair_absmq_score(stage3pair),first_absmq,/*second_absmq*/0,
		Stage3pair_mapq_score(stage3pair),chromosome_iit,
		/*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
		Stage3pair_pairlength(stage3pair),/*chrpos*/chrpos5,/*mate_chrpos*/chrpos3,
		/*clipdir*/0,hardclip5,hardclip3,resulttype,/*first_read_p*/true,/*npaths_mate*/npaths,
		quality_shift,sam_read_group_id,invert_first_p,invert_second_p,
		merge_samechr_p);

      /* print second end */
      SAM_print(fp,abbrev,hit3,/*mate*/hit5,acc1,acc2,/*pathnum*/1,/*npaths*/1,
		Stage3pair_absmq_score(stage3pair),first_absmq,/*second_absmq*/0,
		Stage3pair_mapq_score(stage3pair),chromosome_iit,
		/*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
		Stage3pair_pairlength(stage3pair),/*chrpos*/chrpos3,/*mate_chrpos*/chrpos5,
		/*clipdir*/0,hardclip5,hardclip3,resulttype,/*first_read_p*/false,/*npaths_mate*/npaths,
		quality_shift,sam_read_group_id,invert_second_p,invert_first_p,
		merge_samechr_p);

    } else if (resulttype == PAIRED_MULT) {
      stage3pairarray = (Stage3pair_T *) Result_array(&npaths,&first_absmq,&second_absmq,result);

      if (quiet_if_excessive_p && npaths > maxpaths_report) {
	/* Print as nomapping, but send to fp_concordant_mult */
	SAM_print_nomapping(fp_paired_mult,ABBREV_PAIRED_MULT,
			    queryseq1,/*mate*/(Stage3end_T) NULL,
			    acc1,acc2,chromosome_iit,resulttype,
			    /*first_read_p*/true,/*npaths_mate*/npaths,
			    /*mate_chrpos*/0U,quality_shift,
			    sam_read_group_id,invert_first_p,invert_second_p);
	SAM_print_nomapping(fp_paired_mult,ABBREV_PAIRED_MULT,
			    queryseq2,/*mate*/(Stage3end_T) NULL,
			    acc1,acc2,chromosome_iit,resulttype,
			    /*first_read_p*/false,/*npaths_mate*/npaths,
			    /*mate_chrpos*/0U,quality_shift,
			    sam_read_group_id,invert_second_p,invert_first_p);

      } else {
	/* Stage3pair_eval(stage3pairarray,npaths,maxpaths_report,queryseq1,queryseq2); */
	for (pathnum = 1; pathnum <= npaths && pathnum <= maxpaths_report; pathnum++) {

	  stage3pair = stage3pairarray[pathnum-1];
	  hardclip5 = hardclip3 = 0;

	  hit5 = Stage3pair_hit5(stage3pair);
	  hit3 = Stage3pair_hit3(stage3pair);
	  chrpos5 = SAM_compute_chrpos(/*hardclip_low*/&ignore,/*hardclip_high*/&hardclip5,
				       /*clipdir*/0,hardclip5,hardclip3,/*first_read_p*/true,hit5,
				       Stage3end_substring_low(hit5),Shortread_fulllength(queryseq1));
	  chrpos3 = SAM_compute_chrpos(/*hardclip_low*/&hardclip3,/*hardclip_high*/&ignore,
				       /*clipdir*/0,hardclip5,hardclip3,/*first_read_p*/false,hit3,
				       Stage3end_substring_low(hit3),Shortread_fulllength(queryseq2));

	  /* print first end */
	  SAM_print(fp_paired_mult,ABBREV_PAIRED_MULT,
		    hit5,/*mate*/hit3,acc1,acc2,pathnum,npaths,
		    Stage3pair_absmq_score(stage3pair),first_absmq,second_absmq,
		    Stage3pair_mapq_score(stage3pair),chromosome_iit,
		    /*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
		    Stage3pair_pairlength(stage3pair),/*chrpos*/chrpos5,/*mate_chrpos*/chrpos3,
		    /*clipdir*/0,hardclip5,hardclip3,resulttype,/*first_read_p*/true,/*npaths_mate*/npaths,
		    quality_shift,sam_read_group_id,invert_first_p,invert_second_p,
		    merge_samechr_p);

	  /* print second end */
	  SAM_print(fp_paired_mult,ABBREV_PAIRED_MULT,
		    hit3,/*mate*/hit5,acc1,acc2,pathnum,npaths,
		    Stage3pair_absmq_score(stage3pair),first_absmq,second_absmq,
		    Stage3pair_mapq_score(stage3pair),chromosome_iit,
		    /*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
		    Stage3pair_pairlength(stage3pair),/*chrpos*/chrpos3,/*mate_chrpos*/chrpos5,
		    /*clipdir*/0,hardclip5,hardclip3,resulttype,/*first_read_p*/false,/*npaths_mate*/npaths,
		    quality_shift,sam_read_group_id,invert_second_p,invert_first_p,
		    merge_samechr_p);
	}
      }

    } else if (resulttype == UNPAIRED_UNIQ) {
      /* Even though they are not related, we should print mate information in this situation */
      stage3array1 = (Stage3end_T *) Result_array(&npaths1,&first_absmq1,&second_absmq1,result);
      stage3array2 = (Stage3end_T *) Result_array2(&npaths2,&first_absmq2,&second_absmq2,result);

      hardclip5 = hardclip3 = 0;

      hit5 = stage3array1[0];
      hit3 = stage3array2[0];
      chrpos5 = SAM_compute_chrpos(/*hardclip_low*/&ignore,/*hardclip_high*/&hardclip5,
				   /*clipdir*/0,hardclip5,hardclip3,/*first_read_p*/true,hit5,
				   Stage3end_substring_low(hit5),Shortread_fulllength(queryseq1));
      chrpos3 = SAM_compute_chrpos(/*hardclip_low*/&hardclip3,/*hardclip_high*/&ignore,
				   /*clipdir*/0,hardclip5,hardclip3,/*first_read_p*/false,hit3,
				   Stage3end_substring_low(hit3),Shortread_fulllength(queryseq2));

      if (Stage3end_circularpos(hit5) > 0 || Stage3end_circularpos(hit3) > 0) {
	fp = fp_unpaired_circular;
	abbrev = ABBREV_UNPAIRED_CIRCULAR;
      } else {
	fp = fp_unpaired_uniq;
	abbrev = ABBREV_UNPAIRED_UNIQ;
      }

      /* print first end */
      /* Stage3end_eval_and_sort(stage3array1,npaths1,maxpaths_report,queryseq1); */
      SAM_print(fp,abbrev,hit5,/*mate*/hit3,acc1,acc2,/*pathnum*/1,/*npaths*/1,
		Stage3end_absmq_score(stage3array1[0]),first_absmq1,/*second_absmq*/0,
		Stage3end_mapq_score(stage3array1[0]),chromosome_iit,
		/*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
		/*pairedlength*/0U,/*chrpos*/chrpos5,/*mate_chrpos*/chrpos3,
		/*clipdir*/0,/*hardclip5*/0,/*hardclip3*/0,resulttype,/*first_read_p*/true,
		/*npaths_mate*/1,quality_shift,sam_read_group_id,invert_first_p,invert_second_p,
		merge_samechr_p);

      /* print second end */
      /* Stage3end_eval_and_sort(stage3array2,npaths2,maxpaths_report,queryseq2); */
      SAM_print(fp,abbrev,hit3,/*mate*/hit5,acc1,acc2,/*pathnum*/1,/*npaths*/1,
		Stage3end_absmq_score(stage3array2[0]),first_absmq2,/*second_absmq*/0,
		Stage3end_mapq_score(stage3array2[0]),chromosome_iit,
		/*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
		/*pairedlength*/0U,/*chrpos*/chrpos3,/*mate_chrpos*/chrpos5,
		/*clipdir*/0,/*hardclip5*/0,/*hardclip3*/0,resulttype,/*first_read_p*/false,
		/*npaths_mate*/1,quality_shift,sam_read_group_id,invert_second_p,invert_first_p,
		merge_samechr_p);

    } else if (resulttype == UNPAIRED_MULT || resulttype == UNPAIRED_TRANSLOC) {
      if (resulttype == UNPAIRED_MULT) {
	fp = fp_unpaired_mult;
	abbrev = ABBREV_UNPAIRED_MULT;
      } else {
	fp = fp_unpaired_transloc;
	abbrev = ABBREV_UNPAIRED_TRANSLOC;
      }

      stage3array1 = (Stage3end_T *) Result_array(&npaths1,&first_absmq1,&second_absmq1,result);
      stage3array2 = (Stage3end_T *) Result_array2(&npaths2,&first_absmq2,&second_absmq2,result);

#if 0
      /* Do eval and sorting first */
      if (npaths1 == 1) {
	/* Stage3end_eval_and_sort(stage3array1,npaths1,maxpaths_report,queryseq1); */
      } else if (quiet_if_excessive_p && npaths1 > maxpaths_report) {
	/* Don't sort */
      } else {
	/* Stage3end_eval_and_sort(stage3array1,npaths1,maxpaths_report,queryseq1); */
      }
#endif

#if 0
      if (npaths2 == 1) {
	/* Stage3end_eval_and_sort(stage3array2,npaths2,maxpaths_report,queryseq2); */
      } else if (quiet_if_excessive_p && npaths2 > maxpaths_report) {
	/* Don't sort */
      } else {
	/* Stage3end_eval_and_sort(stage3array2,npaths2,maxpaths_report,queryseq2); */
      }
#endif

      /* print first end results */
      if (npaths2 == 0) {
	mate = (Stage3end_T) NULL;
	chrpos3 = 0U;
      } else if (quiet_if_excessive_p && npaths2 > maxpaths_report) {
	mate = (Stage3end_T) NULL;
	chrpos3 = 0U;
      } else {
	mate = stage3array2[0];
	hardclip3 = 0;
	chrpos3 = SAM_compute_chrpos(/*hardclip_low*/&hardclip3,/*hardclip_high*/&ignore,
				     /*clipdir*/0,hardclip5,hardclip3,/*first_read_p*/false,mate,
				     Stage3end_substring_low(mate),Shortread_fulllength(queryseq2));
      }

      if (npaths1 == 1) {
	stage3 = stage3array1[0];
	hardclip5 = 0;
	chrpos5 = SAM_compute_chrpos(/*hardclip_low*/&ignore,/*hardclip_high*/&hardclip5,
				     /*clipdir*/0,hardclip5,hardclip3,/*first_read_p*/true,stage3,
				     Stage3end_substring_low(stage3),Shortread_fulllength(queryseq1));

	SAM_print(fp,abbrev,stage3,mate,acc1,acc2,/*pathnum*/1,npaths1,
		  Stage3end_absmq_score(stage3),first_absmq1,second_absmq1,
		  Stage3end_mapq_score(stage3),chromosome_iit,
		  /*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
		  /*pairedlength*/0U,/*chrpos*/chrpos5,/*mate_chrpos*/chrpos3,
		  /*clipdir*/0,/*hardclip5*/0,/*hardclip3*/0,resulttype,/*first_read_p*/true,
		  /*npaths_mate*/npaths2,quality_shift,sam_read_group_id,
		  invert_first_p,invert_second_p,merge_samechr_p);

      } else if (quiet_if_excessive_p && npaths1 > maxpaths_report) {
	SAM_print_nomapping(fp,abbrev,queryseq1,mate,acc1,acc2,chromosome_iit,
			    resulttype,/*first_read_p*/true,/*npaths_mate*/npaths2,
			    /*mate_chrpos*/chrpos3,
			    quality_shift,sam_read_group_id,invert_first_p,invert_second_p);

      } else {
	for (pathnum = 1; pathnum <= npaths1 && pathnum <= maxpaths_report; pathnum++) {
	  stage3 = stage3array1[pathnum-1];
	  hardclip5 = 0;
	  chrpos5 = SAM_compute_chrpos(/*hardclip_low*/&ignore,/*hardclip_high*/&hardclip5,
				       /*clipdir*/0,hardclip5,hardclip3,/*first_read_p*/true,stage3,
				       Stage3end_substring_low(stage3),Shortread_fulllength(queryseq1));

	  SAM_print(fp,abbrev,stage3,mate,acc1,acc2,pathnum,npaths1,
		    Stage3end_absmq_score(stage3),first_absmq1,second_absmq1,
		    Stage3end_mapq_score(stage3),chromosome_iit,
		    /*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
		    /*pairedlength*/0U,/*chrpos*/chrpos5,/*mate_chrpos*/chrpos3,
		    /*clipdir*/0,/*hardclip5*/0,/*hardclip3*/0,resulttype,/*first_read_p*/true,
		    /*npaths_mate*/npaths2,quality_shift,sam_read_group_id,
		    invert_first_p,invert_second_p,merge_samechr_p);
	}
      }
			  
      /* print second end results */
      if (npaths1 == 0) {
	mate = (Stage3end_T) NULL;
	chrpos5 = 0U;
      } else if (quiet_if_excessive_p && npaths1 > maxpaths_report) {
	mate = (Stage3end_T) NULL;
	chrpos5 = 0U;
      } else {
	mate = stage3array1[0];
	hardclip5 = 0;
	chrpos5 = SAM_compute_chrpos(/*hardclip_low*/&ignore,/*hardclip_high*/&hardclip5,
				     /*clipdir*/0,hardclip5,hardclip3,/*first_read_p*/true,mate,
				     Stage3end_substring_low(mate),Shortread_fulllength(queryseq1));
      }

      if (npaths2 == 1) {
	stage3 = stage3array2[0];
	hardclip3 = 0;
	chrpos3 = SAM_compute_chrpos(/*hardclip_low*/&hardclip3,/*hardclip_high*/&ignore,
				     /*clipdir*/0,hardclip5,hardclip3,/*first_read_p*/false,stage3,
				     Stage3end_substring_low(stage3),Shortread_fulllength(queryseq2));

	SAM_print(fp,abbrev,stage3,mate,acc1,acc2,/*pathnum*/1,npaths2,
		  Stage3end_absmq_score(stage3),first_absmq2,second_absmq2,
		  Stage3end_mapq_score(stage3),chromosome_iit,
		  /*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
		  /*pairedlength*/0U,/*chrpos*/chrpos3,/*mate_chrpos*/chrpos5,
		  /*clipdir*/0,/*hardclip5*/0,/*hardclip3*/0,resulttype,/*first_read_p*/false,
		  /*npaths_mate*/npaths1,quality_shift,sam_read_group_id,
		  invert_second_p,invert_first_p,merge_samechr_p);

      } else if (quiet_if_excessive_p && npaths2 > maxpaths_report) {
	SAM_print_nomapping(fp,abbrev,queryseq2,mate,acc1,acc2,chromosome_iit,
			    resulttype,/*first_read_p*/false,/*npaths_mate*/npaths1,
			    /*mate_chrpos*/chrpos5,
			    quality_shift,sam_read_group_id,invert_second_p,invert_first_p);

      } else {
	for (pathnum = 1; pathnum <= npaths2 && pathnum <= maxpaths_report; pathnum++) {
	  stage3 = stage3array2[pathnum-1];
	  hardclip3 = 0;
	  chrpos3 = SAM_compute_chrpos(/*hardclip_low*/&hardclip3,/*hardclip_high*/&ignore,
				       /*clipdir*/0,hardclip5,hardclip3,/*first_read_p*/false,stage3,
				       Stage3end_substring_low(stage3),Shortread_fulllength(queryseq2));

	  SAM_print(fp,abbrev,stage3,mate,acc1,acc2,pathnum,npaths2,
		    Stage3end_absmq_score(stage3),first_absmq2,second_absmq2,
		    Stage3end_mapq_score(stage3),chromosome_iit,
		    /*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
		    /*pairedlength*/0U,/*chrpos*/chrpos3,/*mate_chrpos*/chrpos5,
		    /*clipdir*/0,/*hardclip5*/0,/*hardclip3*/0,resulttype,/*first_read_p*/false,
		    /*npaths_mate*/npaths1,quality_shift,sam_read_group_id,
		    invert_second_p,invert_first_p,merge_samechr_p);
	}
      }

    } else {
      stage3array1 = (Stage3end_T *) Result_array(&npaths1,&first_absmq1,&second_absmq1,result);
      stage3array2 = (Stage3end_T *) Result_array2(&npaths2,&first_absmq2,&second_absmq2,result);

      if (resulttype == HALFMAPPING_UNIQ) {
	if (npaths1 == 1 && Stage3end_circularpos(stage3array1[0]) > 0) {
	  fp = fp_halfmapping_circular;
	  abbrev = ABBREV_HALFMAPPING_CIRCULAR;
	} else if (npaths2 == 1 && Stage3end_circularpos(stage3array2[0]) > 0) {
	  fp = fp_halfmapping_circular;
	  abbrev = ABBREV_HALFMAPPING_CIRCULAR;
	} else {
	  fp = fp_halfmapping_uniq;
	  abbrev = ABBREV_HALFMAPPING_UNIQ;
	}
      } else if (resulttype == HALFMAPPING_TRANSLOC) {
	fp = fp_halfmapping_transloc;
	abbrev = ABBREV_HALFMAPPING_TRANSLOC;
      } else if (resulttype == HALFMAPPING_MULT) {
	fp = fp_halfmapping_mult;
	abbrev = ABBREV_HALFMAPPING_MULT;
      } else {
	abort();
      }

#if 0
      /* Do eval and sorting first */
      if (npaths1 == 0) {
	/* Nothing to sort */
      } else if (npaths1 == 1) {
	/* Stage3end_eval_and_sort(stage3array1,npaths1,maxpaths_report,queryseq1); */
      } else if (quiet_if_excessive_p && npaths1 > maxpaths_report) {
	/* Don't sort */
      } else {
	/* Stage3end_eval_and_sort(stage3array1,npaths1,maxpaths_report,queryseq1); */
      }
#endif

#if 0
      if (npaths2 == 0) {
	/* Nothing to sort */
      } else if (npaths2 == 1) {
	/* Stage3end_eval_and_sort(stage3array2,npaths2,maxpaths_report,queryseq2); */
      } else if (quiet_if_excessive_p && npaths2 > maxpaths_report) {
	/* Don't sort */
      } else {
	/* Stage3end_eval_and_sort(stage3array2,npaths2,maxpaths_report,queryseq2); */
      }
#endif


      /* print first end results */
      if (npaths2 == 0) {
	mate = (Stage3end_T) NULL;
	chrpos3 = 0U;
      } else if (quiet_if_excessive_p && npaths2 > maxpaths_report) {
	mate = (Stage3end_T) NULL;
	chrpos3 = 0U;
      } else {
	mate = stage3array2[0];
	hardclip3 = 0;
	chrpos3 = SAM_compute_chrpos(/*hardclip_low*/&hardclip3,/*hardclip_high*/&ignore,
				     /*clipdir*/0,hardclip5,hardclip3,/*first_read_p*/false,mate,
				     Stage3end_substring_low(mate),Shortread_fulllength(queryseq2));
      }

      if (npaths1 == 0) {
	/* mate should be non-NULL here */
	SAM_print_nomapping(fp,abbrev,queryseq1,mate,acc1,acc2,chromosome_iit,resulttype,
			    /*first_read_p*/true,/*npaths_mate*/npaths2,
			    /*mate_chrpos*/chrpos3,
			    quality_shift,sam_read_group_id,invert_first_p,invert_second_p);

      } else if (npaths1 == 1) {
	/* mate should be NULL here */

	stage3 = stage3array1[0];
	hardclip5 = 0;
	chrpos5 = SAM_compute_chrpos(/*hardclip_low*/&ignore,/*hardclip_high*/&hardclip5,
				     /*clipdir*/0,hardclip5,hardclip3,/*first_read_p*/true,stage3,
				     Stage3end_substring_low(stage3),Shortread_fulllength(queryseq1));

	SAM_print(fp,abbrev,stage3,mate,acc1,acc2,/*pathnum*/1,npaths1,
		  Stage3end_absmq_score(stage3),first_absmq1,/*second_absmq1*/0,
		  Stage3end_mapq_score(stage3),chromosome_iit,
		  /*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
		  /*pairedlength*/0U,/*chrpos*/chrpos5,/*mate_chrpos*/chrpos3,
		  /*clipdir*/0,/*hardclip5*/0,/*hardclip3*/0,resulttype,/*first_read_p*/true,
		  /*npaths_mate*/npaths2,quality_shift,sam_read_group_id,
		  invert_first_p,invert_second_p,merge_samechr_p);

      } else if (quiet_if_excessive_p && npaths1 > maxpaths_report) {
	/* mate should be NULL here */
	SAM_print_nomapping(fp,abbrev,queryseq1,mate,acc1,acc2,chromosome_iit,resulttype,
			    /*first_read_p*/true,/*npaths_mate*/npaths2,
			    /*mate_chrpos*/chrpos3,
			    quality_shift,sam_read_group_id,invert_first_p,invert_second_p);

      } else {
	/* mate should be NULL here */
	for (pathnum = 1; pathnum <= npaths1 && pathnum <= maxpaths_report; pathnum++) {
	  stage3 = stage3array1[pathnum-1];
	  hardclip5 = 0;
	  chrpos5 = SAM_compute_chrpos(/*hardclip_low*/&ignore,/*hardclip_high*/&hardclip5,
				       /*clipdir*/0,hardclip5,hardclip3,/*first_read_p*/true,stage3,
				       Stage3end_substring_low(stage3),Shortread_fulllength(queryseq1));

	  SAM_print(fp,abbrev,stage3,mate,acc1,acc2,pathnum,npaths1,
		    Stage3end_absmq_score(stage3),first_absmq1,second_absmq1,
		    Stage3end_mapq_score(stage3),chromosome_iit,
		    /*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
		    /*pairedlength*/0U,/*chrpos*/chrpos5,/*mate_chrpos*/chrpos3,
		    /*clipdir*/0,/*hardclip5*/0,/*hardclip3*/0,resulttype,/*first_read_p*/true,
		    /*npaths_mate*/npaths2,quality_shift,sam_read_group_id,
		    invert_first_p,invert_second_p,merge_samechr_p);
	}
      }
			  
      /* print second end results */
      if (npaths1 == 0) {
	mate = (Stage3end_T) NULL;
	chrpos5 = 0U;
      } else if (quiet_if_excessive_p && npaths1 > maxpaths_report) {
	mate = (Stage3end_T) NULL;
	chrpos5 = 0U;
      } else {
	mate = stage3array1[0];
	hardclip5 = 0;
	chrpos5 = SAM_compute_chrpos(/*hardclip_low*/&ignore,/*hardclip_high*/&hardclip5,
				     /*clipdir*/0,hardclip5,hardclip3,/*first_read_p*/true,mate,
				     Stage3end_substring_low(mate),Shortread_fulllength(queryseq1));
      }

      if (npaths2 == 0) {
	/* mate should be non-NULL here */
	SAM_print_nomapping(fp,abbrev,queryseq2,mate,acc1,acc2,chromosome_iit,resulttype,
			    /*first_read_p*/false,/*npaths_mate*/npaths1,
			    /*mate_chrpos*/chrpos5,
			    quality_shift,sam_read_group_id,invert_second_p,invert_first_p);

      } else if (npaths2 == 1) {
	/* mate should be NULL here */

	stage3 = stage3array2[0];
	hardclip3 = 0;
	chrpos3 = SAM_compute_chrpos(/*hardclip_low*/&hardclip3,/*hardclip_high*/&ignore,
				     /*clipdir*/0,hardclip5,hardclip3,/*first_read_p*/false,stage3,
				     Stage3end_substring_low(stage3),Shortread_fulllength(queryseq2));

	SAM_print(fp,abbrev,stage3,mate,acc1,acc2,/*pathnum*/1,npaths2,
		  Stage3end_absmq_score(stage3),first_absmq2,/*second_absmq2*/0,
		  Stage3end_mapq_score(stage3),chromosome_iit,
		  /*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
		  /*pairedlength*/0U,/*chrpos*/chrpos3,/*mate_chrpos*/chrpos5,
		  /*clipdir*/0,/*hardclip5*/0,/*hardclip3*/0,resulttype,/*first_read_p*/false,
		  /*npaths_mate*/npaths1,quality_shift,sam_read_group_id,
		  invert_second_p,invert_first_p,merge_samechr_p);

      } else if (quiet_if_excessive_p && npaths2 > maxpaths_report) {
	/* mate should be NULL here */
	SAM_print_nomapping(fp,abbrev,queryseq2,mate,acc1,acc2,chromosome_iit,resulttype,
			    /*first_read_p*/false,/*npaths_mate*/npaths1,
			    /*mate_chrpos*/chrpos5,
			    quality_shift,sam_read_group_id,invert_second_p,invert_first_p);

      } else {
	/* mate should be NULL here */
	for (pathnum = 1; pathnum <= npaths2 && pathnum <= maxpaths_report; pathnum++) {
	  stage3 = stage3array2[pathnum-1];
	  hardclip3 = 0;
	  chrpos3 = SAM_compute_chrpos(/*hardclip_low*/&hardclip3,/*hardclip_high*/&ignore,
				       /*clipdir*/0,hardclip5,hardclip3,/*first_read_p*/false,stage3,
				       Stage3end_substring_low(stage3),Shortread_fulllength(queryseq2));

	  SAM_print(fp,abbrev,stage3,mate,acc1,acc2,pathnum,npaths2,
		    Stage3end_absmq_score(stage3),first_absmq2,second_absmq2,
		    Stage3end_mapq_score(stage3),chromosome_iit,
		    /*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
		    /*pairedlength*/0U,/*chrpos*/chrpos3,/*mate_chrpos*/chrpos5,
		    /*clipdir*/0,/*hardclip5*/0,/*hardclip3*/0,resulttype,/*first_read_p*/false,
		    /*npaths_mate*/npaths1,quality_shift,sam_read_group_id,
		    invert_second_p,invert_first_p,merge_samechr_p);
	}
      }

    }
  }

  return;
}



