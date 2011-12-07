static char rcsid[] = "$Id: samprint.c 53332 2011-11-29 22:20:52Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "samprint.h"
#include "samflags.h"
#include <stdlib.h>
#include <ctype.h>

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


/* print_cigar */
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
static int maxpaths;

void
SAM_setup (bool quiet_if_excessive_p_in, int maxpaths_in) {
  quiet_if_excessive_p = quiet_if_excessive_p_in;
  maxpaths = maxpaths_in;
  return;
}


unsigned int
SAM_compute_flag (bool plusp, Stage3end_T mate, Resulttype_T resulttype,
		  bool first_read_p, int pathnum, int npaths, int npaths_mate,
		  bool invertp, bool invert_mate_p) {
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

    } else if (quiet_if_excessive_p && npaths_mate > maxpaths) {
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
      if (Stage3end_chrnum(mate) == 0) {
	/* Splice without a direction */

      } else if (Stage3end_plusp(mate) == invert_mate_p) {
	debug(printf("MATE_MINUSP %d\n",MATE_MINUSP));
	flag |= MATE_MINUSP;
      }

    } else if (resulttype == PAIRED_UNIQ || resulttype == PAIRED_MULT) {
      /* Note: We are counting PAIRED_UNIQ and PAIRED_MULT as "paired" mappings.
	 However, we are no longer counting UNPAIRED_UNIQ as a "paired" mapping. */
      debug(printf("PAIRED_MAPPING %d\n",PAIRED_MAPPING));
      flag |= PAIRED_MAPPING;
      if (Stage3end_chrnum(mate) == 0) {
	/* Splice without a direction */

      } else if (Stage3end_plusp(mate) == invert_mate_p) {
	debug(printf("MATE_MINUSP %d\n",MATE_MINUSP));
	flag |= MATE_MINUSP;
      }

    } else {
      if (Stage3end_chrnum(mate) == 0) {
	/* Splice without a direction */
      } else if (Stage3end_plusp(mate) == invert_mate_p) {
	debug(printf("MATE_MINUSP %d\n",MATE_MINUSP));
	flag |= MATE_MINUSP;
      }
    }
  }

  if (pathnum > 1) {
    debug(printf("NOT_PRIMARY %d\n",NOT_PRIMARY));
    flag |= NOT_PRIMARY;
  }

  return flag;
}


Genomicpos_T
SAM_compute_chrpos (int *hardclip_low, int *hardclip_high, Stage3end_T this,
		    Substring_T substring_low, int querylength) {
  Genomicpos_T chrpos;
  Substring_T substring;
  int querystart, queryend;
  bool plusp;

  if (this == NULL) {
    return 0U;

  } else if (Stage3end_hittype(this) == GMAP) {
    chrpos = Pair_genomicpos_low(&(*hardclip_low),&(*hardclip_high),
				 Stage3end_pairarray(this),Stage3end_npairs(this),
				 querylength,/*watsonp*/Stage3end_plusp(this));

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
      if (querylength - *hardclip_high < Substring_queryend(substring_low)) {
	substring = (Substring_T) NULL;
	while (*hardclip_high < querylength && substring == NULL) {
	  substring = Stage3end_substring_containing(this,querylength - 1 - (*hardclip_high));
	  (*hardclip_high)++;
	}
	(*hardclip_high)--;

	chrpos = Substring_alignend_trim(substring) - Substring_chroffset(substring) + 1U;
	queryend = Substring_queryend(substring);
	querystart = Substring_querystart(Stage3end_substring_high(this));

#if 0
	fprintf(stderr,"case 3, hardclip_low %d, hardclip_high %d, querystart %d, queryend %d\n",
		*hardclip_low,*hardclip_high,querystart,queryend);
#endif

	if (*hardclip_low >= queryend || querylength - *hardclip_high < querystart) {
	  /* fprintf(stderr,"Returning 0\n"); */
	  return 0U;
	}

#if 0
	fprintf(stderr,"Adding to chrpos %u: queryend %d - (querylength %d - hardclip_high %d)\n",
		chrpos,queryend,querylength,*hardclip_high);
#endif

	chrpos += queryend - (querylength - (*hardclip_high));
      } else {
	querystart = Substring_querystart(Stage3end_substring_high(this));
	queryend = Substring_queryend(substring_low);

#if 0
	fprintf(stderr,"case 4, hardclip_low %d, hardclip_high %d, querystart %d, queryend %d\n",
		*hardclip_low,*hardclip_high,querystart,queryend);
#endif

	if (*hardclip_low >= queryend || querylength - *hardclip_high < querystart) {
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
print_chromosomal_pos (FILE *fp, Chrnum_T chrnum, Genomicpos_T chrpos,
		       IIT_T chromosome_iit) {
  bool allocp;
  char *chr;

  if (chrpos == 0U) {
    /* No mapping */
    fprintf(fp,"\t*\t0");
    return;

  } else if (chrnum == 0) {
    /* Interchromosomal splice */
    fprintf(stderr,"Trying to print interchrosomal splice in one line\n");
    abort();

  } else {
    chr = IIT_label(chromosome_iit,chrnum,&allocp);
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
			    Genomicpos_T mate_chrpos, Chrnum_T anchor_chrnum, Genomicpos_T anchor_chrpos,
			    IIT_T chromosome_iit) {
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
      chr = IIT_label(chromosome_iit,mate_effective_chrnum,&allocp);
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
      chr = IIT_label(chromosome_iit,mate_chrnum,&allocp);
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
SAM_print_nomapping (FILE *fp, Shortread_T queryseq, Stage3end_T mate, char *acc,
		     IIT_T chromosome_iit, Resulttype_T resulttype, bool first_read_p,
		     int npaths_mate, Genomicpos_T mate_chrpos,
		     int quality_shift, char *sam_read_group_id, bool invertp, bool invert_mate_p) {
  unsigned int flag;


  /* 1. QNAME */
  fprintf(fp,"%s",acc);
  
  /* 2. FLAG */
  flag = SAM_compute_flag(/*plusp (NA)*/true,mate,resulttype,first_read_p,
			  /*pathnum*/0,/*npaths*/0,npaths_mate,invertp,invert_mate_p);
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

  fprintf(fp,"\n");
  return;
}


static void
print_cigar (FILE *fp, char type, int stringlength, int querypos, int querylength,
	     int hardclip_low, int hardclip_high, bool plusp, int lastp) {
  int matchlength = 0;
  int startpos, endpos;
  int cliplength = 0;
  
  if (plusp == true) {
    debug1(printf("\nEntering print_cigar with stringlength %d, querypos %d, querylength %d, hardclip_low %d, hardclip_high %d, plus\n",
		  stringlength,querypos,querylength,hardclip_low,hardclip_high));
    if (hardclip_low >= querypos) {
      startpos = hardclip_low;
      cliplength = hardclip_low;
    } else {
      startpos = querypos;
    }

    if (querylength - hardclip_high < querypos + stringlength) {
      endpos = querylength - hardclip_high;
    } else {
      endpos = querypos + stringlength;
    }

  } else {
    debug1(printf("\nEntering print_cigar with stringlength %d, querypos %d, querylength %d, hardclip_low %d, hardclip_high %d, minus\n",
		  stringlength,querypos,querylength,hardclip_low,hardclip_high));

    querypos = querylength - querypos - stringlength;
    debug1(printf("  Revising querypos to be %d\n",querypos));

    if (hardclip_high >= querypos) {
      startpos = hardclip_high;
      cliplength = hardclip_high;
    } else {
      startpos = querypos;
    }

    if (querylength - hardclip_low < querypos + stringlength) {
      endpos = querylength - hardclip_low;
    } else {
      endpos = querypos + stringlength;
    }
  }

  debug1(printf("  startpos %d, endpos %d, cliplength %d\n",startpos,endpos,cliplength));

  if (endpos > startpos) {
    matchlength = endpos - startpos;
    if (matchlength > 0) {
      if (cliplength > 0) {
	fprintf(fp,"%dH",cliplength);
      }
      fprintf(fp,"%d%c",matchlength,type);
    }
  }


  if (lastp == true) {
    cliplength = querypos + stringlength - endpos;
    if (cliplength > 0) {
      fprintf(fp,"%dH",cliplength);
    }
  }

  return;
}


static int
print_md_string (int *nmismatches, FILE *fp, int matchlength, char *genomicfwd, int stringlength,
		 int querypos, int querylength, int hardclip_low, int hardclip_high,
		 bool plusp, bool lastp) {
  int starti, endi, i;
  bool hardclip_end_p = false;

  if (plusp == true) {
    debug2(printf("\nEntering md_string with matchlength %d, querypos %d, querylength %d, hardclip_low %d, hardclip_high %d, plus: %s\n",
		  matchlength,querypos,querylength,hardclip_low,hardclip_high,genomicfwd));
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

    if (genomicfwd == NULL) {
      if (endi > starti) {
	matchlength += (endi - starti);
      }
    } else {
      for (i = starti; i < endi; i++) {
	if (isupper(genomicfwd[i])) {
	  matchlength++;
	} else {
	  if (matchlength > 0 || hardclip_end_p == true) {
	    fprintf(fp,"%d",matchlength);
	    hardclip_end_p = false;
	  }
	  fprintf(fp,"%c",toupper(genomicfwd[i]));
	  (*nmismatches) += 1;
	  matchlength = 0;
	}
      }
    }

  } else {
    debug2(printf("\nEntering md_string with matchlength %d, querypos %d, querylength %d, hardclip_low %d, hardclip_high %d, minus: %s\n",
		  matchlength,querypos,querylength,hardclip_low,hardclip_high,genomicfwd));
    querypos = querylength - querypos - stringlength;
    debug2(printf("  Revising querypos to be %d\n",querypos));

    if (hardclip_high == 0) {
      starti = 0;
      hardclip_end_p = true;
    } else if (hardclip_high > querypos) {
      /* startpos = hardclip_high; */
      starti = hardclip_high - querypos;
      hardclip_end_p = true;
      debug2(printf("  Setting starti %d = hardclip_high %d - querypos %d\n",
		    starti,hardclip_high,querypos));
    } else {
      /* startpos = querypos; */
      starti = 0;
    }

    if (querylength - hardclip_low < querypos + stringlength) {
      /* endpos = querylength - hardclip_low; */
      endi = (querylength - hardclip_low) - querypos;
      debug2(printf("  Setting endi %d = (querylength %d - hardclip_low %d) - querypos %d\n",
		    endi,querylength,hardclip_low,querypos));
    } else {
      /* endpos = querypos + stringlength; */
      endi = stringlength;
    }

    debug2(printf("  Counting matches from %d to %d\n",starti,endi));

    if (genomicfwd == NULL) {
      if (endi > starti) {
	matchlength += (endi - starti);
      }
    } else {
      for (i = starti; i < endi; i++) {
	if (isupper(genomicfwd[i])) {
	  matchlength++;
	} else {
	  if (matchlength > 0 || hardclip_end_p == true) {
	    fprintf(fp,"%d",matchlength);
	    hardclip_end_p = false;
	  }
	  fprintf(fp,"%c",toupper(genomicfwd[i]));
	  (*nmismatches) += 1;
	  matchlength = 0;
	}
      }
    }
  }

  debug2(printf("  Ending with matchlength %d\n",matchlength));

  if (lastp == false) {
    return matchlength;
  } else if (matchlength > 0) {
    fprintf(fp,"%d",matchlength);
    return 0;
  } else {
    return 0;
  }
}


static void
print_single (FILE *fp, Stage3end_T this, Stage3end_T mate, char *acc, int pathnum, int npaths,
	      int absmq_score, int second_absmq, int mapq_score, IIT_T chromosome_iit, Shortread_T queryseq, int pairedlength,
	      Genomicpos_T chrpos, Genomicpos_T mate_chrpos, int hardclip5, int hardclip3,
	      Resulttype_T resulttype, bool first_read_p,
	      int npaths_mate, int quality_shift,
	      char *sam_read_group_id, bool invertp, bool invert_mate_p) {
  unsigned int flag = 0U;
  Substring_T substring;
  int nmismatches = 0, querylength, substring_start, substring_length;
  int hardclip_low, hardclip_high, mate_hardclip_low, mate_hardclip_high;
  char *genomicfwd, *genomicdir;
  bool plusp;


  querylength = Shortread_fulllength(queryseq);
  plusp = Stage3end_plusp(this);
  substring = Stage3end_substring1(this);

  if (first_read_p == true) {
    hardclip_low = 0;
    hardclip_high = hardclip5;
    mate_hardclip_low = hardclip3;
    mate_hardclip_high = 0;
    /* fprintf(stderr,"first read: hardclip_low = %d, hardclip_high = %d\n",hardclip_low,hardclip_high); */
  } else {
    hardclip_low = hardclip3;
    hardclip_high = 0;
    mate_hardclip_low = 0;
    mate_hardclip_high = hardclip5;
    /* fprintf(stderr,"second read: hardclip_low = %d, hardclip_high = %d\n",hardclip_low,hardclip_high); */
  }


  /* 1. QNAME */
  fprintf(fp,"%s",acc);

  /* 2. FLAG */
  flag = SAM_compute_flag(Stage3end_plusp(this),mate,resulttype,first_read_p,
			  pathnum,npaths,npaths_mate,invertp,invert_mate_p);
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
		/*querypos*/Substring_queryend(substring),querylength,hardclip_low,hardclip_high,
		/*plusp*/false,/*lastp*/false);
    print_cigar(fp,/*type*/'M',Substring_match_length(substring),
		/*querypos*/Substring_querystart(substring),querylength,
		hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
    print_cigar(fp,/*type*/'S',Substring_querystart(substring),
		/*querypos*/0,querylength,
		hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
  }


  /* 7. MRNM: Mate chr */
  /* 8. MPOS: Mate chrpos */
  print_mate_chromosomal_pos(fp,Stage3end_chrnum(mate),Stage3end_effective_chrnum(mate),mate_chrpos,
			     Stage3end_chrnum(this),chrpos,chromosome_iit);


  /* 9. ISIZE: Insert size */
#ifdef PAIRED_ZERO_ISIZE
  if (resulttype == CONCORDANT_UNIQ || resulttype == CONCORDANT_TRANSLOC || resulttype == CONCORDANT_MULT) {
    if (plusp == invertp) {
      fprintf(fp,"\t%d",-pairedlength);
    } else {
      fprintf(fp,"\t%d",pairedlength);
    }
  } else {
    fprintf(fp,"\t0");
  }
#else
  if (mate_chrpos == 0) {
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
#endif


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
  substring_start = Substring_querystart(substring);
  substring_length = Substring_match_length(substring);

  fprintf(fp,"MD:Z:");  
  if ((genomicdir = Substring_genomic_refdiff(substring)) == NULL) {
    if (plusp == true) {
      print_md_string(&nmismatches,fp,/*matchlength*/0,/*genomicfwd*/NULL,substring_length,
		      /*querypos*/substring_start,querylength,hardclip_low,hardclip_high,
		      /*plusp*/true,/*lastp*/true);
    } else {
      print_md_string(&nmismatches,fp,/*matchlength*/0,/*genomicfwd*/NULL,substring_length,
		      /*querypos*/substring_start,querylength,
		      hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
    }

  } else if (plusp == true) {
    print_md_string(&nmismatches,fp,/*matchlength*/0,&(genomicdir[substring_start]),substring_length,
		    /*querypos*/substring_start,querylength,hardclip_low,hardclip_high,
		    /*plusp*/true,/*lastp*/true);

  } else {
    genomicfwd = (char *) CALLOC(querylength+1,sizeof(char));
    make_complement_buffered(genomicfwd,&(genomicdir[substring_start]),substring_length);
    print_md_string(&nmismatches,fp,/*matchlength*/0,genomicfwd,substring_length,
		    /*querypos*/substring_start,querylength,
		    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
    FREE(genomicfwd);
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
  fprintf(fp,"NM:i:%d",nmismatches);

  /* 12. TAGS: SM */
  fprintf(fp,"\t");
  fprintf(fp,"SM:i:%d",mapq_score);

  /* 12. TAGS: XQ */
  fprintf(fp,"\t");
  fprintf(fp,"XQ:i:%d",absmq_score);

  /* 12. TAGS: X2 */
  fprintf(fp,"\t");
  fprintf(fp,"X2:i:%d",second_absmq);

  fprintf(fp,"\n");
  return;
}


static void
print_insertion (FILE *fp, Stage3end_T this, Stage3end_T mate, char *acc, int pathnum, int npaths,
		 int absmq_score, int second_absmq, int mapq_score, IIT_T chromosome_iit, Shortread_T queryseq, int pairedlength,
		 Genomicpos_T chrpos, Genomicpos_T mate_chrpos, int hardclip5, int hardclip3,
		 Resulttype_T resulttype, bool first_read_p, int npaths_mate,
		 int quality_shift, char *sam_read_group_id, bool invertp, bool invert_mate_p) {
  unsigned int flag = 0U;
  Substring_T substring1, substring2;
  int nmismatches = 0, querylength;
  char *genomicfwd, *genomicdir;
  int substring1_start, substring2_start, substring1_length, substring2_length, matchlength, nindels;
  int hardclip_low, hardclip_high, mate_hardclip_low, mate_hardclip_high;
  bool plusp;

  querylength = Shortread_fulllength(queryseq);
  plusp = Stage3end_plusp(this);

  substring1 = Stage3end_substring1(this);
  substring2 = Stage3end_substring2(this);

  if (first_read_p == true) {
    hardclip_low = 0;
    hardclip_high = hardclip5;
    mate_hardclip_low = hardclip3;
    mate_hardclip_high = 0;
    /* fprintf(stderr,"first read: hardclip_low = %d, hardclip_high = %d\n",hardclip_low,hardclip_high); */
  } else {
    hardclip_low = hardclip3;
    hardclip_high = 0;
    mate_hardclip_low = 0;
    mate_hardclip_high = hardclip5;
    /* fprintf(stderr,"second read: hardclip_low = %d, hardclip_high = %d\n",hardclip_low,hardclip_high); */
  }


  /* 1. QNAME */
  fprintf(fp,"%s",acc);

  /* 2. FLAG */
  flag = SAM_compute_flag(plusp,mate,resulttype,first_read_p,
			  pathnum,npaths,npaths_mate,invertp,invert_mate_p);
  fprintf(fp,"\t%u",flag);

  /* 3. RNAME: chr */
  /* 4. POS: chrpos */
  print_chromosomal_pos(fp,Stage3end_chrnum(this),chrpos,chromosome_iit);


  /* 5. MAPQ: Mapping quality */
  fprintf(fp,"\t%d",mapq_score);

  /* 6. CIGAR */
  fprintf(fp,"\t");

  nindels = Stage3end_nindels(this);
  if (plusp == true) {
    print_cigar(fp,/*type*/'S',Substring_querystart(substring1),
		/*querypos*/0,querylength,hardclip_low,hardclip_high,
		/*plusp*/true,/*lastp*/false);
    print_cigar(fp,/*type*/'M',Substring_match_length(substring1),
		/*querypos*/Substring_querystart(substring1),querylength,
		hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);
    print_cigar(fp,/*type*/'I',nindels,
		/*querypos*/Substring_queryend(substring1),querylength,
		hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);
    print_cigar(fp,/*type*/'M',Substring_match_length(substring2),
		/*querypos*/Substring_querystart(substring2),querylength,
		hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);
    print_cigar(fp,/*type*/'S',querylength - Substring_queryend(substring2),
		/*querypos*/Substring_queryend(substring2),querylength,
		hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);

  } else {
    print_cigar(fp,/*type*/'S',querylength - Substring_queryend(substring2),
		/*querypos*/Substring_queryend(substring2),querylength,hardclip_low,hardclip_high,
		/*plusp*/false,/*lastp*/false);
    print_cigar(fp,/*type*/'M',Substring_match_length(substring2),
		/*querypos*/Substring_querystart(substring2),querylength,
		hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
    print_cigar(fp,/*type*/'I',nindels,
		/*querypos*/Substring_queryend(substring1),querylength,
		hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);
    print_cigar(fp,/*type*/'M',Substring_match_length(substring1),
		/*querypos*/Substring_querystart(substring1),querylength,
		hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
    print_cigar(fp,/*type*/'S',Substring_querystart(substring1),
		/*querypos*/0,querylength,
		hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
  }

  /* 7. MRNM: Mate chr */
  /* 8. MPOS: Mate chrpos */
  print_mate_chromosomal_pos(fp,Stage3end_chrnum(mate),Stage3end_effective_chrnum(mate),mate_chrpos,
			     Stage3end_chrnum(this),chrpos,chromosome_iit);


  /* 9. ISIZE: Insert size */
#ifdef PAIRED_ZERO_ISIZE
  if (resulttype == CONCORDANT_UNIQ || resulttype == CONCORDANT_TRANSLOC || resulttype == CONCORDANT_MULT) {
    if (plusp == invertp) {
      fprintf(fp,"\t%d",-pairedlength);
    } else {
      fprintf(fp,"\t%d",pairedlength);
    }
  } else {
    fprintf(fp,"\t0");
  }
#else
  if (mate_chrpos == 0) {
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
#endif

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

  substring1_start = Substring_querystart(substring1);
  substring1_length = Substring_match_length(substring1);
  substring2_start = Substring_querystart(substring2);
  substring2_length = Substring_match_length(substring2);

  if (plusp == true) {
    genomicfwd = Substring_genomic_refdiff(substring1);
    matchlength = print_md_string(&nmismatches,fp,/*matchlength*/0,&(genomicfwd[substring1_start]),substring1_length,
				  /*querypos*/substring1_start,querylength,hardclip_low,hardclip_high,
				  /*plusp*/true,/*lastp*/false);

#if 0
    /* If MD string is supposed to include insertion, then uncomment this */
    matchlength += nindels;
#endif

    genomicfwd = Substring_genomic_refdiff(substring2);
    print_md_string(&nmismatches,fp,matchlength,&(genomicfwd[substring2_start]),substring2_length,
		    /*querypos*/substring2_start,querylength,hardclip_low,hardclip_high,
		    /*plusp*/true,/*lastp*/true);
  } else {
    genomicfwd = (char *) CALLOC(substring2_length+1,sizeof(char));
    genomicdir = Substring_genomic_refdiff(substring2);
    make_complement_buffered(genomicfwd,&(genomicdir[substring2_start]),substring2_length);
    matchlength = print_md_string(&nmismatches,fp,/*matchlength*/0,genomicfwd,substring2_length,
				  /*querypos*/substring2_start,querylength,
				  hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
    FREE(genomicfwd);

#if 0
    /* If MD string is supposed to include insertion, then uncomment this */
    matchlength += nindels;
#endif

    genomicfwd = (char *) CALLOC(substring1_length+1,sizeof(char));
    genomicdir = Substring_genomic_refdiff(substring1);
    make_complement_buffered(genomicfwd,&(genomicdir[substring1_start]),substring1_length);
    print_md_string(&nmismatches,fp,matchlength,genomicfwd,substring1_length,
		    /*querypos*/substring1_start,querylength,
		    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
    FREE(genomicfwd);
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
  fprintf(fp,"NM:i:%d",nmismatches + nindels);

  /* 12. TAGS: SM */
  fprintf(fp,"\t");
  fprintf(fp,"SM:i:%d",mapq_score);

  /* 12. TAGS: XQ */
  fprintf(fp,"\t");
  fprintf(fp,"XQ:i:%d",absmq_score);

  /* 12. TAGS: X2 */
  fprintf(fp,"\t");
  fprintf(fp,"X2:i:%d",second_absmq);


  fprintf(fp,"\n");
  return;
}


static void
print_deletion (FILE *fp, Stage3end_T this, Stage3end_T mate, char *acc, int pathnum, int npaths,
		int absmq_score, int second_absmq, int mapq_score, IIT_T chromosome_iit, Shortread_T queryseq, int pairedlength,
		Genomicpos_T chrpos, Genomicpos_T mate_chrpos, int hardclip5, int hardclip3,
		Resulttype_T resulttype, bool first_read_p, int npaths_mate,
		int quality_shift, char *sam_read_group_id, bool invertp, bool invert_mate_p) {
  unsigned int flag = 0U;
  Substring_T substring1, substring2;
  int nmismatches = 0, querylength;
  char *genomicfwd, *genomicdir;
  int substring1_start, substring2_start, substring1_length, substring2_length, nindels;
  int hardclip_low, hardclip_high, mate_hardclip_low, mate_hardclip_high;
  bool plusp;

  querylength = Shortread_fulllength(queryseq);
  plusp = Stage3end_plusp(this);

  substring1 = Stage3end_substring1(this);
  substring2 = Stage3end_substring2(this);

  if (first_read_p == true) {
    hardclip_low = 0;
    hardclip_high = hardclip5;
    mate_hardclip_low = hardclip3;
    mate_hardclip_high = 0;
    /* fprintf(stderr,"first read: hardclip_low = %d, hardclip_high = %d\n",hardclip_low,hardclip_high); */
  } else {
    hardclip_low = hardclip3;
    hardclip_high = 0;
    mate_hardclip_low = 0;
    mate_hardclip_high = hardclip5;
    /* fprintf(stderr,"second read: hardclip_low = %d, hardclip_high = %d\n",hardclip_low,hardclip_high); */
  }


  /* 1. QNAME */
  fprintf(fp,"%s",acc);

  /* 2. FLAG */
  flag = SAM_compute_flag(plusp,mate,resulttype,first_read_p,
			  pathnum,npaths,npaths_mate,invertp,invert_mate_p);
  fprintf(fp,"\t%u",flag);

  /* 3. RNAME: chr */
  /* 4. POS: chrpos */
  print_chromosomal_pos(fp,Stage3end_chrnum(this),chrpos,chromosome_iit);


  /* 5. MAPQ: Mapping quality */
  fprintf(fp,"\t%d",mapq_score);

  /* 6. CIGAR */
  fprintf(fp,"\t");

  nindels = Stage3end_nindels(this); /* nindels is positive */
  if (plusp == true) {
    print_cigar(fp,/*type*/'S',Substring_querystart(substring1),
		/*querypos*/0,querylength,hardclip_low,hardclip_high,
		/*plusp*/true,/*lastp*/false);
    print_cigar(fp,/*type*/'M',Substring_match_length(substring1),
		/*querypos*/Substring_querystart(substring1),querylength,
		hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);
    if (hardclip_low < Substring_querystart(substring2)) {
      fprintf(fp,"%dD",nindels);
    }
    print_cigar(fp,/*type*/'M',Substring_match_length(substring2),
		/*querypos*/Substring_querystart(substring2),querylength,
		hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);
    print_cigar(fp,/*type*/'S',querylength - Substring_queryend(substring2),
		/*querypos*/Substring_queryend(substring2),querylength,
		hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);

  } else {
    print_cigar(fp,/*type*/'S',querylength - Substring_queryend(substring2),
		/*querypos*/Substring_queryend(substring2),querylength,
		hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
    print_cigar(fp,/*type*/'M',Substring_match_length(substring2),
		/*querypos*/Substring_querystart(substring2),querylength,
		hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
    if (querylength - hardclip_high > Substring_queryend(substring1)) {
      fprintf(fp,"%dD",nindels);
    }
    print_cigar(fp,/*type*/'M',Substring_match_length(substring1),
		/*querypos*/Substring_querystart(substring1),querylength,
		hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
    print_cigar(fp,/*type*/'S',Substring_querystart(substring1),
		/*querypos*/0,querylength,
		hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
  }

  /* 7. MRNM: Mate chr */
  /* 8. MPOS: Mate chrpos */
  print_mate_chromosomal_pos(fp,Stage3end_chrnum(mate),Stage3end_effective_chrnum(mate),mate_chrpos,
			     Stage3end_chrnum(this),chrpos,chromosome_iit);


  /* 9. ISIZE: Insert size */
#ifdef PAIRED_ZERO_ISIZE
  if (resulttype == CONCORDANT_UNIQ || resulttype == CONCORDANT_TRANSLOC || resulttype == CONCORDANT_MULT) {
    if (plusp == invertp) {
      fprintf(fp,"\t%d",-pairedlength);
    } else {
      fprintf(fp,"\t%d",pairedlength);
    }
  } else {
    fprintf(fp,"\t0");
  }
#else
  if (mate_chrpos == 0) {
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
#endif


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

  substring1_start = Substring_querystart(substring1);
  substring1_length = Substring_match_length(substring1);
  substring2_start = Substring_querystart(substring2);
  substring2_length = Substring_match_length(substring2);

  if (plusp == true) {
    genomicfwd = Substring_genomic_refdiff(substring1);
    print_md_string(&nmismatches,fp,/*matchlength*/0,&(genomicfwd[substring1_start]),substring1_length,
		    /*querypos*/substring1_start,querylength,hardclip_low,hardclip_high,
		    /*plusp*/true,/*lastp*/true);

    if (hardclip_low < substring2_start) {
      /* Deletion string: Potential problem if followed by a mismatch, but can be resolved by looking at CIGAR string */
      fprintf(fp,"^%s",Stage3end_deletion_string(this));
    }

    genomicfwd = Substring_genomic_refdiff(substring2);
    print_md_string(&nmismatches,fp,/*matchlength*/0,&(genomicfwd[substring2_start]),substring2_length,
		    /*querypos*/substring2_start,querylength,hardclip_low,hardclip_high,
		    /*plusp*/true,/*lastp*/true);
  } else {
    genomicfwd = (char *) CALLOC(substring2_length+1,sizeof(char));
    genomicdir = Substring_genomic_refdiff(substring2);
    make_complement_buffered(genomicfwd,&(genomicdir[substring2_start]),substring2_length);
    print_md_string(&nmismatches,fp,/*matchlength*/0,genomicfwd,substring2_length,
		    /*querypos*/substring2_start,querylength,
		    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
    FREE(genomicfwd);

    if (querylength - hardclip_high > substring1_start + substring1_length) {
      /* Deletion string: Potential problem if followed by a mismatch, but can be resolved by looking at CIGAR string */
      genomicfwd = (char *) CALLOC(nindels+1,sizeof(char));
      make_complement_buffered(genomicfwd,Stage3end_deletion_string(this),nindels);
      fprintf(fp,"^%s",genomicfwd);
      FREE(genomicfwd);
    }

    genomicfwd = (char *) CALLOC(substring1_length+1,sizeof(char));
    genomicdir = Substring_genomic_refdiff(substring1);
    make_complement_buffered(genomicfwd,&(genomicdir[substring1_start]),substring1_length);
    print_md_string(&nmismatches,fp,/*matchlength*/0,genomicfwd,substring1_length,
		    /*querypos*/substring1_start,querylength,
		    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
    FREE(genomicfwd);
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
  fprintf(fp,"NM:i:%d",nmismatches + nindels);

  /* 12. TAGS: SM */
  fprintf(fp,"\t");
  fprintf(fp,"SM:i:%d",mapq_score);

  /* 12. TAGS: XQ */
  fprintf(fp,"\t");
  fprintf(fp,"XQ:i:%d",absmq_score);

  /* 12. TAGS: X2 */
  fprintf(fp,"\t");
  fprintf(fp,"X2:i:%d",second_absmq);

  fprintf(fp,"\n");
  return;
}


static void
halfdonor_dinucleotide (char *donor1, char *donor2, Substring_T donor) {
  bool sensep;
  char *genomic;
  int substring_start, substring_length;

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
print_halfdonor (FILE *fp, Substring_T donor, Stage3end_T this, Stage3end_T mate,
		 char *acc, int pathnum, int npaths, int absmq_score, int second_absmq, int mapq_score,
		 IIT_T chromosome_iit, Shortread_T queryseq, int pairedlength,
		 Genomicpos_T concordant_chrpos, Genomicpos_T chrpos, Genomicpos_T mate_chrpos,
		 int hardclip5, int hardclip3, Resulttype_T resulttype, bool first_read_p, int npaths_mate,
		 int quality_shift, char *sam_read_group_id, bool invertp, bool invert_mate_p,
		 bool use_hardclip_p, bool print_xt_p, char donor1, char donor2, char acceptor2, char acceptor1,
		 double donor_prob, double acceptor_prob) {
  unsigned int flag = 0U;
  int nmismatches = 0, querylength;
  bool sensep;
  char *genomicfwd, *genomicdir;
  int substring_start, substring_length;
  int hardclip_low, hardclip_high, mate_hardclip_low, mate_hardclip_high;
  int transloc_hardclip_low, transloc_hardclip_high;
  bool plusp;

  querylength = Shortread_fulllength(queryseq);
  plusp = Substring_plusp(donor);

  if (first_read_p == true) {
    hardclip_low = 0;
    hardclip_high = hardclip5;
    mate_hardclip_low = hardclip3;
    mate_hardclip_high = 0;
    /* fprintf(stderr,"first read: hardclip_low = %d, hardclip_high = %d\n",hardclip_low,hardclip_high); */
  } else {
    hardclip_low = hardclip3;
    hardclip_high = 0;
    mate_hardclip_low = 0;
    mate_hardclip_high = hardclip5;
    /* fprintf(stderr,"second read: hardclip_low = %d, hardclip_high = %d\n",hardclip_low,hardclip_high); */
  }


  /* 1. QNAME */
  fprintf(fp,"%s",acc);

  /* 2. FLAG */
  flag = SAM_compute_flag(plusp,mate,resulttype,first_read_p,
			  pathnum,npaths,npaths_mate,invertp,invert_mate_p);
  fprintf(fp,"\t%u",flag);

  /* 3. RNAME: chr */
  /* 4. POS: chrpos */
  print_chromosomal_pos(fp,Substring_chrnum(donor),chrpos,chromosome_iit);
  

  /* 5. MAPQ: Mapping quality */
  fprintf(fp,"\t%d",mapq_score);

  /* 6. CIGAR */
  fprintf(fp,"\t");
  sensep = Substring_chimera_sensep(donor);

  if (use_hardclip_p == true) {
    if (sensep == plusp) {
      if (plusp == true) {
	/* sensep true */
	transloc_hardclip_low = 0;
	transloc_hardclip_high = querylength - Substring_queryend(donor);

      } else {
	/* sensep false */
	transloc_hardclip_low = Substring_querystart(donor);
	transloc_hardclip_high = 0;
      }

    } else { /* sensep != Substring_plusp(donor) */
      if (plusp == true) {
	transloc_hardclip_low = Substring_querystart(donor);
	transloc_hardclip_high = 0;

      } else {
	transloc_hardclip_low = 0;
	transloc_hardclip_high = querylength - Substring_queryend(donor);
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
		  /*querypos*/Substring_queryend(donor),querylength,
		  hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
      print_cigar(fp,/*type*/'M',Substring_match_length(donor),
		  /*querypos*/Substring_querystart(donor),querylength,
		  hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
      print_cigar(fp,/*type*/use_hardclip_p ? 'H' : 'S',Substring_querystart(donor),
		  /*querypos*/0,querylength,
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
		  /*querypos*/Substring_queryend(donor),querylength,
		  hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
      print_cigar(fp,/*type*/'M',Substring_match_length(donor),
		  /*querypos*/Substring_querystart(donor),querylength,
		  hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
      print_cigar(fp,/*type*/'S',Substring_querystart(donor),
		  /*querypos*/0,querylength,hardclip_low,hardclip_high,
		  /*plusp*/false,/*lastp*/true);
    }
  }

  /* 7. MRNM: Mate chr */
  /* 8. MPOS: Mate chrpos */
  print_mate_chromosomal_pos(fp,Stage3end_chrnum(mate),Stage3end_effective_chrnum(mate),mate_chrpos,
			     Stage3end_chrnum(this),chrpos,chromosome_iit);


  /* 9. ISIZE: Insert size */
#ifdef PAIRED_ZERO_ISIZE
  if (resulttype == CONCORDANT_UNIQ || resulttype == CONCORDANT_TRANSLOC || resulttype == CONCORDANT_MULT) {
    if (plusp == invertp) {
      fprintf(fp,"\t%d",-pairedlength);
    } else {
      fprintf(fp,"\t%d",pairedlength);
    }
  } else {
    fprintf(fp,"\t0");
  }
#else
  if (mate_chrpos == 0) {
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
#endif


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

  substring_start = Substring_querystart(donor);
  substring_length = Substring_match_length(donor);

  if (use_hardclip_p == false) {
    genomicdir = Substring_genomic_refdiff(donor);
    if (plusp == true) {
      print_md_string(&nmismatches,fp,/*matchlength*/0,&(genomicdir[substring_start]),substring_length,
		      /*querypos*/substring_start,querylength,hardclip_low,hardclip_high,
		      /*plusp*/true,/*lastp*/true);
    } else {
      genomicfwd = (char *) CALLOC(querylength+1,sizeof(char));
      make_complement_buffered(genomicfwd,&(genomicdir[substring_start]),substring_length);
      print_md_string(&nmismatches,fp,/*matchlength*/0,genomicfwd,substring_length,
		      /*querypos*/substring_start,querylength,
		      hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
      FREE(genomicfwd);
    }

  } else if (sensep == true) {
    if (plusp == true) {
      genomicfwd = Substring_genomic_refdiff(donor);
      print_md_string(&nmismatches,fp,/*matchlength*/0,&(genomicfwd[substring_start]),substring_length,
		      /*querypos*/substring_start,querylength,hardclip_low,hardclip_high,
		      /*plusp*/true,/*lastp*/true);
    } else {
      genomicfwd = (char *) CALLOC(substring_length+1,sizeof(char));
      genomicdir = Substring_genomic_refdiff(donor);
      make_complement_buffered(genomicfwd,&(genomicdir[substring_start]),substring_length);
      print_md_string(&nmismatches,fp,/*matchlength*/0,genomicfwd,substring_length,
		      /*querypos*/substring_start,querylength,
		      hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
      FREE(genomicfwd);
    }

  } else {			/* sensep == false */
    if (plusp == true) {
      genomicfwd = Substring_genomic_refdiff(donor);
      print_md_string(&nmismatches,fp,/*matchlength*/0,&(genomicfwd[substring_start]),substring_length,
		      /*querypos*/substring_start,querylength,hardclip_low,hardclip_high,
		      /*plusp*/true,/*lastp*/true);
    } else {
      genomicfwd = (char *) CALLOC(substring_length+1,sizeof(char));
      genomicdir = Substring_genomic_refdiff(donor);
      make_complement_buffered(genomicfwd,&(genomicdir[substring_start]),substring_length);
      print_md_string(&nmismatches,fp,/*matchlength*/0,genomicfwd,substring_length,
		      /*querypos*/substring_start,querylength,
		      hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
      FREE(genomicfwd);
    }
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
  fprintf(fp,"NM:i:%d",nmismatches);
  
  /* 12. TAGS: SM */
  fprintf(fp,"\t");
  fprintf(fp,"SM:i:%d",mapq_score);

  /* 12. TAGS: XQ */
  fprintf(fp,"\t");
  fprintf(fp,"XQ:i:%d",absmq_score);

  /* 12. TAGS: X2 */
  fprintf(fp,"\t");
  fprintf(fp,"X2:i:%d",second_absmq);

  /* 12. TAGS: XS */
  fprintf(fp,"\t");
  if (sensep == plusp) {
    fprintf(fp,"XS:A:+");
  } else {
    fprintf(fp,"XS:A:-");
  }

  /* 12. TAGS: XT */
  if (print_xt_p == true) {
    fprintf(fp,"\t");
    fprintf(fp,"XT:Z:%c%c-%c%c,%.2f,%.2f",donor1,donor2,acceptor2,acceptor1,donor_prob,acceptor_prob);
  }

  fprintf(fp,"\n");
  return;
}


static void
print_halfacceptor (FILE *fp, Substring_T acceptor, Stage3end_T this, Stage3end_T mate,
		    char *acc, int pathnum, int npaths, int absmq_score, int second_absmq, int mapq_score,
		    IIT_T chromosome_iit, Shortread_T queryseq, int pairedlength,
		    Genomicpos_T concordant_chrpos, Genomicpos_T chrpos, Genomicpos_T mate_chrpos,
		    int hardclip5, int hardclip3, Resulttype_T resulttype, bool first_read_p, int npaths_mate,
		    int quality_shift, char *sam_read_group_id, bool invertp, bool invert_mate_p,
		    bool use_hardclip_p, bool print_xt_p, char donor1, char donor2, char acceptor2, char acceptor1,
		    double donor_prob, double acceptor_prob) {
  unsigned int flag = 0U;
  int nmismatches = 0, querylength;
  bool sensep;
  char *genomicfwd, *genomicdir;
  int substring_start, substring_length;
  int hardclip_low, hardclip_high, mate_hardclip_low, mate_hardclip_high;
  int transloc_hardclip_low, transloc_hardclip_high;
  bool plusp;

  querylength = Shortread_fulllength(queryseq);
  plusp = Substring_plusp(acceptor);

  if (first_read_p == true) {
    hardclip_low = 0;
    hardclip_high = hardclip5;
    mate_hardclip_low = hardclip3;
    mate_hardclip_high = 0;
    /* fprintf(stderr,"first read: hardclip_low = %d, hardclip_high = %d\n",hardclip_low,hardclip_high); */
  } else {
    hardclip_low = hardclip3;
    hardclip_high = 0;
    mate_hardclip_low = 0;
    mate_hardclip_high = hardclip5;
    /* fprintf(stderr,"second read: hardclip_low = %d, hardclip_high = %d\n",hardclip_low,hardclip_high); */
  }


  /* 1. QNAME */
  fprintf(fp,"%s",acc);

  /* 2. FLAG */
  flag = SAM_compute_flag(plusp,mate,resulttype,first_read_p,
			  pathnum,npaths,npaths_mate,invertp,invert_mate_p);
  fprintf(fp,"\t%u",flag);

  /* 3. RNAME: chr */
  /* 4. POS: chrpos */
  print_chromosomal_pos(fp,Substring_chrnum(acceptor),chrpos,chromosome_iit);


  /* 5. MAPQ: Mapping quality */
  fprintf(fp,"\t%d",mapq_score);

  /* 6. CIGAR */
  fprintf(fp,"\t");
  sensep = Substring_chimera_sensep(acceptor);

  if (use_hardclip_p == true) {
    if (sensep != plusp) {
      if (plusp == true) {
	/* sensep false */
	transloc_hardclip_low = 0;
	transloc_hardclip_high = querylength - Substring_queryend(acceptor);

      } else {
	/* sensep true */
	transloc_hardclip_low = Substring_querystart(acceptor);
	transloc_hardclip_high = 0;
      }

    } else { /*  sensep == Substring_plusp(acceptor) */
      if (plusp == true) {
	transloc_hardclip_low = Substring_querystart(acceptor);
	transloc_hardclip_high = 0;

      } else {
	transloc_hardclip_low = 0;
	transloc_hardclip_high = querylength - Substring_queryend(acceptor);
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
		  /*querypos*/Substring_queryend(acceptor),querylength,
		  hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
      print_cigar(fp,/*type*/'M',Substring_match_length(acceptor),
		  /*querypos*/Substring_querystart(acceptor),querylength,
		  hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
      print_cigar(fp,/*type*/use_hardclip_p ? 'H' : 'S',Substring_querystart(acceptor),
		  /*querypos*/0,querylength,hardclip_low,hardclip_high,
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
		  /*querypos*/Substring_queryend(acceptor),querylength,
		  hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
      print_cigar(fp,/*type*/'M',Substring_match_length(acceptor),
		  /*querypos*/Substring_querystart(acceptor),querylength,
		  hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
      print_cigar(fp,/*type*/'S',Substring_querystart(acceptor),
		  /*querypos*/0,querylength,hardclip_low,hardclip_high,
		  /*plusp*/false,/*lastp*/true);
    }
  }


  /* 7. MRNM: Mate chr */
  /* 8. MPOS: Mate chrpos */
  print_mate_chromosomal_pos(fp,Stage3end_chrnum(mate),Stage3end_effective_chrnum(mate),mate_chrpos,
			     Stage3end_chrnum(this),chrpos,chromosome_iit);


  /* 9. ISIZE: Insert size */
#ifdef PAIRED_ZERO_ISIZE
  if (resulttype == CONCORDANT_UNIQ || resulttype == CONCORDANT_TRANSLOC || resulttype == CONCORDANT_MULT) {
    if (plusp == invertp) {
      fprintf(fp,"\t%d",-pairedlength);
    } else {
      fprintf(fp,"\t%d",pairedlength);
    }
  } else {
    fprintf(fp,"\t0");
  }
#else
  if (mate_chrpos == 0) {
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
#endif


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

  substring_start = Substring_querystart(acceptor);
  substring_length = Substring_match_length(acceptor);

  if (use_hardclip_p == false) {
    genomicdir = Substring_genomic_refdiff(acceptor);
    if (plusp == true) {
      print_md_string(&nmismatches,fp,/*matchlength*/0,&(genomicdir[substring_start]),substring_length,
		      /*querypos*/substring_start,querylength,hardclip_low,hardclip_high,
		      /*plusp*/true,/*lastp*/true);
    } else {
      genomicfwd = (char *) CALLOC(querylength+1,sizeof(char));
      make_complement_buffered(genomicfwd,&(genomicdir[substring_start]),substring_length);
      print_md_string(&nmismatches,fp,/*matchlength*/0,genomicfwd,substring_length,
		      /*querypos*/substring_start,querylength,
		      hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
      FREE(genomicfwd);
    }

  } else if (sensep == false) {
    if (plusp == true) {
      genomicfwd = Substring_genomic_refdiff(acceptor);
      print_md_string(&nmismatches,fp,/*matchlength*/0,&(genomicfwd[substring_start]),substring_length,
		      /*querypos*/substring_start,querylength,hardclip_low,hardclip_high,
		      /*plusp*/true,/*lastp*/true);
    } else {
      genomicfwd = (char *) CALLOC(substring_length+1,sizeof(char));
      genomicdir = Substring_genomic_refdiff(acceptor);
      make_complement_buffered(genomicfwd,&(genomicdir[substring_start]),substring_length);
      print_md_string(&nmismatches,fp,/*matchlength*/0,genomicfwd,substring_length,
		      /*querypos*/substring_start,querylength,
		      hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
      FREE(genomicfwd);
    }

  } else {			/* sensep true */
    if (plusp == true) {
      genomicfwd = Substring_genomic_refdiff(acceptor);
      print_md_string(&nmismatches,fp,/*matchlength*/0,&(genomicfwd[substring_start]),substring_length,
		      /*querypos*/substring_start,querylength,hardclip_low,hardclip_high,
		      /*plusp*/true,/*lastp*/true);
    } else {
      genomicfwd = (char *) CALLOC(substring_length+1,sizeof(char));
      genomicdir = Substring_genomic_refdiff(acceptor);
      make_complement_buffered(genomicfwd,&(genomicdir[substring_start]),substring_length);
      print_md_string(&nmismatches,fp,/*matchlength*/0,genomicfwd,substring_length,
		      /*querypos*/substring_start,querylength,
		      hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
      FREE(genomicfwd);
    }
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
  fprintf(fp,"NM:i:%d",nmismatches);
  
  /* 12. TAGS: SM */
  fprintf(fp,"\t");
  fprintf(fp,"SM:i:%d",mapq_score);

  /* 12. TAGS: XQ */
  fprintf(fp,"\t");
  fprintf(fp,"XQ:i:%d",absmq_score);

  /* 12. TAGS: X2 */
  fprintf(fp,"\t");
  fprintf(fp,"X2:i:%d",second_absmq);

  /* 12. TAGS: XS */
  fprintf(fp,"\t");
  if (sensep == plusp) {
    fprintf(fp,"XS:A:+");
  } else {
    fprintf(fp,"XS:A:-");
  }

  /* 12. TAGS: XT */
  if (print_xt_p == true) {
    fprintf(fp,"\t");
    fprintf(fp,"XT:Z:%c%c-%c%c,%.2f,%.2f",donor1,donor2,acceptor2,acceptor1,donor_prob,acceptor_prob);
  }

  fprintf(fp,"\n");
  return;
}



static void
print_localsplice (FILE *fp, Stage3end_T this, Stage3end_T mate, char *acc, int pathnum, int npaths,
		   int absmq_score, int second_absmq, int mapq_score, IIT_T chromosome_iit, Shortread_T queryseq, int pairedlength,
		   Genomicpos_T chrpos, Genomicpos_T mate_chrpos, int hardclip5, int hardclip3,
		   Resulttype_T resulttype, bool first_read_p, int npaths_mate,
		   int quality_shift, char *sam_read_group_id, bool invertp, bool invert_mate_p) {
  unsigned int flag = 0U;
  Substring_T substring1, substring2;
  int nmismatches = 0, querylength;
  bool sensep;
  char *genomicfwd, *genomicdir;
  int substring1_start, substring2_start, substring1_length, substring2_length, matchlength;
  int hardclip_low, hardclip_high, mate_hardclip_low, mate_hardclip_high;
  bool plusp;

  querylength = Shortread_fulllength(queryseq);
  plusp = Stage3end_plusp(this);

  sensep = (Stage3end_sensedir(this) == SENSE_FORWARD);

  if (first_read_p == true) {
    hardclip_low = 0;
    hardclip_high = hardclip5;
    mate_hardclip_low = hardclip3;
    mate_hardclip_high = 0;
    /* fprintf(stderr,"first read: hardclip_low = %d, hardclip_high = %d\n",hardclip_low,hardclip_high); */
  } else {
    hardclip_low = hardclip3;
    hardclip_high = 0;
    mate_hardclip_low = 0;
    mate_hardclip_high = hardclip5;
    /* fprintf(stderr,"second read: hardclip_low = %d, hardclip_high = %d\n",hardclip_low,hardclip_high); */
  }


  /* 1. QNAME */
  fprintf(fp,"%s",acc);

  /* 2. FLAG */
  flag = SAM_compute_flag(plusp,mate,resulttype,first_read_p,
			  pathnum,npaths,npaths_mate,invertp,invert_mate_p);
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
    if (hardclip_low < Substring_queryend(substring1) && querylength - hardclip_high > Substring_querystart(substring2)) {
      debug1(printf("\nhardclip_low %d < queryend(substring1) %d && querylength %d - hardclip_high %d > querystart(substring2) %d\n",
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
		/*querypos*/Substring_queryend(substring1),querylength,
		hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
    print_cigar(fp,/*type*/'M',Substring_match_length(substring1),
		/*querypos*/Substring_querystart(substring1),querylength,
		hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
    if (hardclip_low < Substring_queryend(substring2) && querylength - hardclip_high > Substring_querystart(substring1)) {
      debug1(printf("\nhardclip_low %d < queryend(substring2) %d && querylength %d - hardclip_high %d > querystart(substring1) %d\n",
		    hardclip_low,Substring_queryend(substring2),querylength,hardclip_high,Substring_querystart(substring1)));
      fprintf(fp,"%uN",Stage3end_distance(this));
    }
    print_cigar(fp,/*type*/'M',Substring_match_length(substring2),
		/*querypos*/Substring_querystart(substring2),querylength,
		hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
    print_cigar(fp,/*type*/'S',Substring_querystart(substring2),
		/*querypos*/0,querylength,hardclip_low,hardclip_high,
		/*plusp*/false,/*lastp*/true);
  }


  /* 7. MRNM: Mate chr */
  /* 8. MPOS: Mate chrpos */
  print_mate_chromosomal_pos(fp,Stage3end_chrnum(mate),Stage3end_effective_chrnum(mate),mate_chrpos,
			     Stage3end_chrnum(this),chrpos,chromosome_iit);


  /* 9. ISIZE: Insert size */
#ifdef PAIRED_ZERO_ISIZE
  if (resulttype == CONCORDANT_UNIQ || resulttype == CONCORDANT_TRANSLOC || resulttype == CONCORDANT_MULT) {
    if (plusp == invertp) {
      fprintf(fp,"\t%d",-pairedlength);
    } else {
      fprintf(fp,"\t%d",pairedlength);
    }
  } else {
    fprintf(fp,"\t0");
  }
#else
  if (mate_chrpos == 0) {
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
#endif

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


  substring1_start = Substring_querystart(substring1);
  substring1_length = Substring_match_length(substring1);
  substring2_start = Substring_querystart(substring2);
  substring2_length = Substring_match_length(substring2);

  if (plusp == true) {
    genomicfwd = Substring_genomic_refdiff(substring1);
    matchlength = print_md_string(&nmismatches,fp,/*matchlength*/0,&(genomicfwd[substring1_start]),substring1_length,
				  /*querypos*/substring1_start,querylength,hardclip_low,hardclip_high,
				  /*plusp*/true,/*lastp*/false);
    
#if 0
    /* Intron: Gets skipped in MD string */
    fprintf(fp,"^");
#endif

    genomicfwd = Substring_genomic_refdiff(substring2);
    print_md_string(&nmismatches,fp,matchlength,&(genomicfwd[substring2_start]),substring2_length,
		    /*querypos*/substring2_start,querylength,hardclip_low,hardclip_high,
		    /*plusp*/true,/*lastp*/true);

  } else {
    genomicfwd = (char *) CALLOC(substring1_length+1,sizeof(char));
    genomicdir = Substring_genomic_refdiff(substring1);
    make_complement_buffered(genomicfwd,&(genomicdir[substring1_start]),substring1_length);
    matchlength = print_md_string(&nmismatches,fp,/*matchlength*/0,genomicfwd,substring1_length,
				  /*querypos*/substring1_start,querylength,
				  hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
    FREE(genomicfwd);

#if 0
    /* Intron: Gets skipped in MD string */
    fprintf(fp,"^");
#endif

    genomicfwd = (char *) CALLOC(substring2_length+1,sizeof(char));
    genomicdir = Substring_genomic_refdiff(substring2);
    make_complement_buffered(genomicfwd,&(genomicdir[substring2_start]),substring2_length);
    print_md_string(&nmismatches,fp,matchlength,genomicfwd,substring2_length,
		    /*querypos*/substring2_start,querylength,
		    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
    FREE(genomicfwd);
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
  fprintf(fp,"NM:i:%d",nmismatches);
  
  /* 12. TAGS: SM */
  fprintf(fp,"\t");
  fprintf(fp,"SM:i:%d",mapq_score);

  /* 12. TAGS: XQ */
  fprintf(fp,"\t");
  fprintf(fp,"XQ:i:%d",absmq_score);

  /* 12. TAGS: X2 */
  fprintf(fp,"\t");
  fprintf(fp,"X2:i:%d",second_absmq);

  /* 12. TAGS: XS */
  fprintf(fp,"\t");
  if (sensep == plusp) {
    fprintf(fp,"XS:A:+");
  } else {
    fprintf(fp,"XS:A:-");
  }

  fprintf(fp,"\n");
  return;
}


static void
print_shortexon (FILE *fp, Stage3end_T shortexon, Stage3end_T mate, char *acc, int pathnum, int npaths,
		 int absmq_score, int second_absmq, int mapq_score, IIT_T chromosome_iit, Shortread_T queryseq, int pairedlength,
		 Genomicpos_T chrpos, Genomicpos_T mate_chrpos, int hardclip5, int hardclip3,
		 Resulttype_T resulttype, bool first_read_p, int npaths_mate,
		 int quality_shift, char *sam_read_group_id, bool invertp, bool invert_mate_p) {
  unsigned int flag = 0U;
  /* substring1 is low coordinate on genome, substring2 is high */
  Substring_T substring1, substring2, substringM;
  Genomicpos_T distance1, distance2;
  int nmismatches = 0, querylength;
  bool sensep;
  char *genomicfwd, *genomicdir;
  int substring1_start, substring2_start, substringM_start,
    substring1_length, substring2_length, substringM_length, matchlength;
  int hardclip_low, hardclip_high, mate_hardclip_low, mate_hardclip_high;
  bool plusp;

  querylength = Shortread_fulllength(queryseq);
  plusp = Stage3end_plusp(shortexon);

  sensep = (Stage3end_sensedir(shortexon) == SENSE_FORWARD);

  if (first_read_p == true) {
    hardclip_low = 0;
    hardclip_high = hardclip5;
    mate_hardclip_low = hardclip3;
    mate_hardclip_high = 0; 
    /* fprintf(stderr,"first read: hardclip_low = %d, hardclip_high = %d\n",hardclip_low,hardclip_high); */
 } else {
    hardclip_low = hardclip3;
    hardclip_high = 0;
    mate_hardclip_low = 0;
    mate_hardclip_high = hardclip5;
    /* fprintf(stderr,"second read: hardclip_low = %d, hardclip_high = %d\n",hardclip_low,hardclip_high); */
  }


  /* 1. QNAME */
  fprintf(fp,"%s",acc);

  /* 2. FLAG */
  flag = SAM_compute_flag(plusp,mate,resulttype,first_read_p,
			  pathnum,npaths,npaths_mate,invertp,invert_mate_p);
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
		  /*querypos*/Substring_queryend(substringM),querylength,
		  hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
    }

  } else if (plusp == true) {
    print_cigar(fp,/*type*/'S',Substring_querystart(substring1),
		/*querypos*/0,querylength,hardclip_low,hardclip_high,
		/*plusp*/true,/*lastp*/false);
    print_cigar(fp,/*type*/'M',Substring_match_length(substring1),
		/*querypos*/Substring_querystart(substring1),querylength,
		hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/false);
    if (hardclip_low < Substring_queryend(substring1) && querylength - hardclip_high > Substring_querystart(substringM)) {
      debug1(printf("\nhardclip_low %d < queryend(substring1) %d && querylength %d - hardclip_high %d > querystart(substringM) %d\n",
		    hardclip_low,Substring_queryend(substring1),querylength,hardclip_high,Substring_querystart(substringM)));
      fprintf(fp,"%uN",distance1);
    }

  } else {
    print_cigar(fp,/*type*/'S',querylength - Substring_queryend(substring1),
		/*querypos*/Substring_queryend(substring1),querylength,
		hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
    print_cigar(fp,/*type*/'M',Substring_match_length(substring1),
		/*querypos*/Substring_querystart(substring1),querylength,
		hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
    if (hardclip_low < Substring_queryend(substringM) && querylength - hardclip_high > Substring_querystart(substring1)) {
      debug1(printf("\nhardclip_low %d < queryend(substringM) %d && querylength %d - hardclip_high %d > querystart(substring1) %d\n",
		    hardclip_low,Substring_queryend(substringM),querylength,hardclip_high,Substring_querystart(substring1)));
      fprintf(fp,"%uN",distance1);
    }
  }

  print_cigar(fp,/*type*/'M',Substring_match_length(substringM),
		/*querypos*/Substring_querystart(substringM),querylength,
		hardclip_low,hardclip_high,plusp,/*lastp*/false);

  if (substring2 == NULL) {
    if (plusp == true) {
      print_cigar(fp,/*type*/'S',querylength - Substring_queryend(substringM),
		  /*querypos*/Substring_queryend(substringM),querylength,
		  hardclip_low,hardclip_high,/*plusp*/true,/*lastp*/true);
    } else {
      print_cigar(fp,/*type*/'S',Substring_querystart(substringM),
		  /*querypos*/0,querylength,
		  hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
    }

  } else if (plusp == true) {
    if (hardclip_low < Substring_queryend(substringM) && querylength - hardclip_high > Substring_querystart(substring2)) {
      debug1(printf("\nhardclip_low %d < queryend(substringM) %d && querylength %d - hardclip_high %d > querystart(substring2) %d\n",
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
    if (hardclip_low < Substring_queryend(substring2) && querylength - hardclip_high > Substring_querystart(substringM)) {
      debug1(printf("\nhardclip_low %d < queryend(substring2) %d && querylength %d - hardclip_high %d > querystart(substringM) %d\n",
		    hardclip_low,Substring_queryend(substring2),querylength,hardclip_high,Substring_querystart(substringM)));
      fprintf(fp,"%uN",distance2);
    }
    print_cigar(fp,/*type*/'M',Substring_match_length(substring2),
		/*querypos*/Substring_querystart(substring2),querylength,
		hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
    print_cigar(fp,/*type*/'S',Substring_querystart(substring2),
		/*querypos*/0,querylength,hardclip_low,hardclip_high,
		/*plusp*/false,/*lastp*/true);
  }


  /* 7. MRNM: Mate chr */
  /* 8. MPOS: Mate chrpos */
  print_mate_chromosomal_pos(fp,Stage3end_chrnum(mate),Stage3end_effective_chrnum(mate),mate_chrpos,
			     Stage3end_chrnum(shortexon),chrpos,chromosome_iit);


  /* 9. ISIZE: Insert size */
#ifdef PAIRED_ZERO_ISIZE
  if (resulttype == CONCORDANT_UNIQ || resulttype == CONCORDANT_TRANSLOC || resulttype == CONCORDANT_MULT) {
    if (plusp == invertp) {
      fprintf(fp,"\t%d",-pairedlength);
    } else {
      fprintf(fp,"\t%d",pairedlength);
    }
  } else {
    fprintf(fp,"\t0");
  }
#else
  if (mate_chrpos == 0) {
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
#endif


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
      genomicfwd = Substring_genomic_refdiff(substring1);
      matchlength = print_md_string(&nmismatches,fp,/*matchlength*/0,&(genomicfwd[substring1_start]),substring1_length,
				    /*querypos*/substring1_start,querylength,hardclip_low,hardclip_high,
				    /*plusp*/true,/*lastp*/false);
    }

#if 0
    /* Intron 1: Gets skipped in MD string */
    fprintf(fp,"^");
#endif

    genomicfwd = Substring_genomic_refdiff(substringM);
    matchlength = print_md_string(&nmismatches,fp,matchlength,&(genomicfwd[substringM_start]),substringM_length,
				  /*querypos*/substringM_start,querylength,hardclip_low,hardclip_high,
				  /*plusp*/true,/*lastp*/false);

#if 0
    /* Intron 2: Gets skipped in MD string */
    fprintf(fp,"^");
#endif

    if (substring2 == NULL) {
      /* Equivalent: if (matchlength > 0) fprintf(fp,"%d",matchlength); */
      print_md_string(&nmismatches,fp,matchlength,/*genomicfwd*/NULL,/*substring2_length*/0,
		      /*querypos*/0,querylength,hardclip_low,hardclip_high,
		      /*plusp*/true,/*lastp*/true);
    } else {
      genomicfwd = Substring_genomic_refdiff(substring2);
      print_md_string(&nmismatches,fp,matchlength,&(genomicfwd[substring2_start]),substring2_length,
		      /*querypos*/substring2_start,querylength,hardclip_low,hardclip_high,
		      /*plusp*/true,/*lastp*/true);
    }

  } else {

    if (substring1 == NULL) {
      matchlength = 0;
    } else {
      genomicfwd = (char *) CALLOC(substring1_length+1,sizeof(char));
      genomicdir = Substring_genomic_refdiff(substring1);
      make_complement_buffered(genomicfwd,&(genomicdir[substring1_start]),substring1_length);
      matchlength = print_md_string(&nmismatches,fp,/*matchlength*/0,genomicfwd,substring1_length,
				    /*querypos*/substring1_start,querylength,
				    hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
      FREE(genomicfwd);
    }

#if 0
    /* Intron 1: Gets skipped in MD string */
    fprintf(fp,"^");
#endif

    genomicfwd = (char *) CALLOC(substringM_length+1,sizeof(char));
    genomicdir = Substring_genomic_refdiff(substringM);
    make_complement_buffered(genomicfwd,&(genomicdir[substringM_start]),substringM_length);
    matchlength = print_md_string(&nmismatches,fp,matchlength,genomicfwd,substringM_length,
				  /*querypos*/substringM_start,querylength,
				  hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/false);
    FREE(genomicfwd);

#if 0
    /* Intron 2: Not sure how to handle in MD string */
    fprintf(fp,"^");
#endif

    if (substring2 == NULL) {
      /* Equivalent: if (matchlength > 0) fprintf(fp,"%d",matchlength); */
      print_md_string(&nmismatches,fp,matchlength,/*genomicfwd*/0,/*substring2_length*/0,
		      /*querypos*/0,querylength,hardclip_low,hardclip_high,
		      /*plusp*/false,/*lastp*/true);
    } else {
      genomicfwd = (char *) CALLOC(substring2_length+1,sizeof(char));
      genomicdir = Substring_genomic_refdiff(substring2);
      make_complement_buffered(genomicfwd,&(genomicdir[substring2_start]),substring2_length);
      print_md_string(&nmismatches,fp,matchlength,genomicfwd,substring2_length,
		      /*querypos*/substring2_start,querylength,
		      hardclip_low,hardclip_high,/*plusp*/false,/*lastp*/true);
      FREE(genomicfwd);
    }
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
  fprintf(fp,"NM:i:%d",nmismatches);
  
  /* 12. TAGS: SM */
  fprintf(fp,"\t");
  fprintf(fp,"SM:i:%d",mapq_score);

  /* 12. TAGS: XQ */
  fprintf(fp,"\t");
  fprintf(fp,"XQ:i:%d",absmq_score);

  /* 12. TAGS: X2 */
  fprintf(fp,"\t");
  fprintf(fp,"X2:i:%d",second_absmq);

  /* 12. TAGS: XS */
  fprintf(fp,"\t");
  if (sensep == plusp) {
    fprintf(fp,"XS:A:+");
  } else {
    fprintf(fp,"XS:A:-");
  }

  fprintf(fp,"\n");
  return;
}



/* Distant splicing, including scramble, inversion, translocation */
static void
print_exon_exon (FILE *fp, Stage3end_T this, Stage3end_T mate, char *acc, int pathnum, int npaths,
		 int absmq_score, int second_absmq, int mapq_score, IIT_T chromosome_iit, Shortread_T queryseq, int pairedlength,
		 Genomicpos_T mate_chrpos, int hardclip5, int hardclip3,
		 Resulttype_T resulttype, bool first_read_p, int npaths_mate,
		 int quality_shift, char *sam_read_group_id, bool invertp, bool invert_mate_p) {
  Genomicpos_T donor_chrpos, acceptor_chrpos, concordant_chrpos;
  Substring_T donor, acceptor;
  int hardclip_low, hardclip_high;
  char donor1, donor2, acceptor2, acceptor1;
  double donor_prob, acceptor_prob;

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
  hardclip_low = hardclip_high = 0;
#endif

  donor_chrpos = SAM_compute_chrpos(&hardclip_low,&hardclip_high,this,/*substring_low*/donor,
				    Shortread_fulllength(queryseq));
  acceptor_chrpos = SAM_compute_chrpos(&hardclip_low,&hardclip_high,this,/*substring_low*/acceptor,
				       Shortread_fulllength(queryseq));
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
    print_halfdonor(fp,donor,this,mate,acc,pathnum,npaths,
		    absmq_score,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
		    concordant_chrpos,donor_chrpos,mate_chrpos,
		    hardclip5,hardclip3,resulttype,first_read_p,
		    npaths_mate,quality_shift,sam_read_group_id,
		    invertp,invert_mate_p,/*use_hardclip_p*/true,/*print_xt_p*/true,
		    donor1,donor2,acceptor2,acceptor1,donor_prob,acceptor_prob);

    print_halfacceptor(fp,acceptor,this,mate,acc,pathnum,npaths,
		       absmq_score,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
		       concordant_chrpos,acceptor_chrpos,mate_chrpos,
		       hardclip5,hardclip3,resulttype,first_read_p,
		       npaths_mate,quality_shift,sam_read_group_id,
		       invertp,invert_mate_p,/*use_hardclip_p*/true,/*print_xt_p*/true,
		       donor1,donor2,acceptor2,acceptor1,donor_prob,acceptor_prob);

  } else if (Stage3end_sensedir(this) == SENSE_ANTI) {
    print_halfacceptor(fp,acceptor,this,mate,acc,pathnum,npaths,
		       absmq_score,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
		       concordant_chrpos,acceptor_chrpos,mate_chrpos,
		       hardclip5,hardclip3,resulttype,first_read_p,
		       npaths_mate,quality_shift,sam_read_group_id,
		       invertp,invert_mate_p,/*use_hardclip_p*/true,/*print_xt_p*/true,
		       donor1,donor2,acceptor2,acceptor1,donor_prob,acceptor_prob);

    print_halfdonor(fp,donor,this,mate,acc,pathnum,npaths,
		    absmq_score,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
		    concordant_chrpos,donor_chrpos,mate_chrpos,
		    hardclip5,hardclip3,resulttype,first_read_p,
		    npaths_mate,quality_shift,sam_read_group_id,
		    invertp,invert_mate_p,/*use_hardclip_p*/true,/*print_xt_p*/true,
		    donor1,donor2,acceptor2,acceptor1,donor_prob,acceptor_prob);

  } else {
    abort();
  }

  return;
}



void
SAM_print (FILE *fp, Stage3end_T this, Stage3end_T mate, char *acc, int pathnum, int npaths,
	   int absmq_score, int second_absmq, int mapq_score, IIT_T chromosome_iit, Shortread_T queryseq,
	   Shortread_T queryseq_mate, int pairedlength, Genomicpos_T chrpos, Genomicpos_T mate_chrpos,
	   int hardclip5, int hardclip3, Resulttype_T resulttype, bool first_read_p,
	   int npaths_mate, int quality_shift, char *sam_read_group_id, bool invertp, bool invert_mate_p,
	   bool merge_samechr_p) {
  Hittype_T hittype;
  Substring_T donor, acceptor;
  bool sensep, normalp;
  unsigned int flag;
  int ignore = 0;


  hittype = Stage3end_hittype(this);
  if (chrpos == 0U) {
    SAM_print_nomapping(fp,queryseq,mate,acc,chromosome_iit,resulttype,first_read_p,
			npaths_mate,mate_chrpos,quality_shift,
			sam_read_group_id,invertp,invert_mate_p);

  } else if (hittype == EXACT || hittype == SUB || hittype == TERMINAL) {
    print_single(fp,this,mate,acc,pathnum,npaths,
		 absmq_score,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
		 chrpos,mate_chrpos,hardclip5,hardclip3,resulttype,first_read_p,
		 npaths_mate,quality_shift,sam_read_group_id,
		 invertp,invert_mate_p);

  } else if (hittype == INSERTION) {
    print_insertion(fp,this,mate,acc,pathnum,npaths,
		    absmq_score,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
		    chrpos,mate_chrpos,hardclip5,hardclip3,resulttype,first_read_p,
		    npaths_mate,quality_shift,sam_read_group_id,
		    invertp,invert_mate_p);

  } else if (hittype == DELETION) {
    print_deletion(fp,this,mate,acc,pathnum,npaths,
		   absmq_score,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
		   chrpos,mate_chrpos,hardclip5,hardclip3,resulttype,first_read_p,
		   npaths_mate,quality_shift,sam_read_group_id,
		   invertp,invert_mate_p);

  } else if (hittype == HALFSPLICE_DONOR) {
    print_halfdonor(fp,Stage3end_substring_donor(this),this,mate,acc,pathnum,npaths,
		    absmq_score,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
		    /*concordant_chrpos*/chrpos,chrpos,mate_chrpos,
		    hardclip5,hardclip3,resulttype,first_read_p,
		    npaths_mate,quality_shift,sam_read_group_id,
		    invertp,invert_mate_p,/*use_hardclip_p*/false,/*print_xt_p*/false,
		    /*donor1*/'X',/*donor2*/'X',/*acceptor2*/'X',/*acceptor1*/'X',
		    /*donor_prob*/0.0,/*acceptor_prob*/0.0);

  } else if (hittype == HALFSPLICE_ACCEPTOR) {
    print_halfacceptor(fp,Stage3end_substring_acceptor(this),this,mate,acc,pathnum,npaths,
		       absmq_score,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
		       /*concordant_chrpos*/chrpos,chrpos,mate_chrpos,
		       hardclip5,hardclip3,resulttype,first_read_p,
		       npaths_mate,quality_shift,sam_read_group_id,
		       invertp,invert_mate_p,/*use_hardclip_p*/false,/*print_xt_p*/false,
		       /*donor1*/'X',/*donor2*/'X',/*acceptor2*/'X',/*acceptor1*/'X',
		       /*donor_prob*/0.0,/*acceptor_prob*/0.0);

  } else if (hittype == SPLICE || hittype == SAMECHR_SPLICE || hittype == TRANSLOC_SPLICE) {
    /* Follows print_splice_distance() in substring.c */
    donor = Stage3end_substring_donor(this);
    acceptor = Stage3end_substring_acceptor(this);

    if (donor == NULL || acceptor == NULL) {
      abort();
    } else if (hittype == TRANSLOC_SPLICE || (hittype == SAMECHR_SPLICE && merge_samechr_p == false)) {
      /* Stage3end_chrnum(this) == 0 || Stage3end_distance(this) == 0U */
      /* distant splice */
      print_exon_exon(fp,this,mate,acc,pathnum,npaths,
		      absmq_score,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
		      mate_chrpos,hardclip5,hardclip3,resulttype,first_read_p,
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
	print_localsplice(fp,this,mate,acc,pathnum,npaths,
			  absmq_score,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
			  chrpos,mate_chrpos,hardclip5,hardclip3,resulttype,first_read_p,
			  npaths_mate,quality_shift,sam_read_group_id,
			  invertp,invert_mate_p);
      } else {
	print_exon_exon(fp,this,mate,acc,pathnum,npaths,
			absmq_score,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
			mate_chrpos,hardclip5,hardclip3,resulttype,first_read_p,
			npaths_mate,quality_shift,sam_read_group_id,
			invertp,invert_mate_p);
      }
    }
      
  } else if (hittype == ONE_THIRD_SHORTEXON || hittype == TWO_THIRDS_SHORTEXON || hittype == SHORTEXON) {
    print_shortexon(fp,this,mate,acc,pathnum,npaths,
		    absmq_score,second_absmq,mapq_score,chromosome_iit,queryseq,pairedlength,
		    chrpos,mate_chrpos,hardclip5,hardclip3,resulttype,first_read_p,
		    npaths_mate,quality_shift,sam_read_group_id,
		    invertp,invert_mate_p);

  } else if (hittype == GMAP) {
    /* Note: sam_paired_p must be true because we are calling GMAP only on halfmapping uniq */

    if (mate == NULL) {
      chrpos = SAM_compute_chrpos(/*hardclip_low*/&ignore,/*hardclip_high*/&hardclip5,this,
				  /*substring_low*/NULL,Shortread_fulllength(queryseq));
      mate_chrpos = 0U;
      hardclip3 = 0;

    } else if (first_read_p == true) {
      chrpos = SAM_compute_chrpos(/*hardclip_low*/&ignore,/*hardclip_high*/&hardclip5,this,
				  /*substring_low*/NULL,Shortread_fulllength(queryseq));
      mate_chrpos = SAM_compute_chrpos(/*hardclip_low*/&hardclip3,/*hardclip_high*/&ignore,mate,
				       Stage3end_substring_low(mate),
				       Shortread_fulllength(queryseq_mate));
    } else {
      chrpos = SAM_compute_chrpos(/*hardclip_low*/&hardclip3,/*hardclip_high*/&ignore,this,
				  /*substring_low*/NULL,Shortread_fulllength(queryseq));
      mate_chrpos = SAM_compute_chrpos(/*hardclip_low*/&ignore,/*hardclip_high*/&hardclip5,mate,
				       Stage3end_substring_low(mate),
				       Shortread_fulllength(queryseq_mate));
    }
    assert(ignore == 0);

    flag = SAM_compute_flag(Stage3end_plusp(this),mate,resulttype,first_read_p,
			    pathnum,npaths,npaths_mate,invertp,invert_mate_p);

    Pair_print_sam(fp,Stage3end_pairarray(this),Stage3end_npairs(this),
		   acc,Stage3end_chrnum(this),chromosome_iit,
		   /*usersegment*/(Sequence_T) NULL,
		   Shortread_fullpointer(queryseq),Shortread_quality_string(queryseq),
		   hardclip5,hardclip3,Shortread_fulllength(queryseq),
		   Stage3end_plusp(this),Stage3end_cdna_direction(this),
		   /*chimera_part*/0,/*chimera*/NULL,quality_shift,first_read_p,
		   pathnum,npaths,absmq_score,second_absmq,flag,
		   /*pair_mapq_score*/mapq_score,/*end_mapq_score*/mapq_score,
		   chrpos,mate_chrpos,pairedlength,sam_read_group_id);
  } else {
    abort();
  }

  return;
}



void
SAM_print_paired (Result_T result, Resulttype_T resulttype,
		  IIT_T chromosome_iit, Shortread_T queryseq1, Shortread_T queryseq2,
		  bool invert_first_p, bool invert_second_p,
		  bool nofailsp, bool failsonlyp, bool fails_as_input_p,
		  bool fastq_format_p, bool clip_overlap_p, bool merge_samechr_p,
		  int quality_shift, char *sam_read_group_id,
		  FILE *fp_nomapping_1, FILE *fp_nomapping_2,
		  FILE *fp_unpaired_uniq, FILE *fp_unpaired_transloc, FILE *fp_unpaired_mult,
		  FILE *fp_halfmapping_uniq, FILE *fp_halfmapping_transloc, FILE *fp_halfmapping_mult,
		  FILE *fp_paired_uniq_inv, FILE *fp_paired_uniq_scr, FILE *fp_paired_uniq_long,
		  FILE *fp_paired_mult, FILE *fp_concordant_uniq, FILE *fp_concordant_transloc,
		  FILE *fp_concordant_mult) {
  Stage3pair_T *stage3pairarray, stage3pair;
  Stage3end_T *stage3array1, *stage3array2, stage3, mate, hit5, hit3;
  Genomicpos_T chrpos5, chrpos3;
  int overlap;
  int npaths, npaths1, npaths2, pathnum;
  int second_absmq, second_absmq1, second_absmq2;
  int hardclip5 = 0, hardclip3 = 0, ignore = 0;
  char *acc;
  Pairtype_T pairtype;
  FILE *fp;

  acc = Shortread_accession(queryseq1);

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
      SAM_print_nomapping(fp_nomapping_1,queryseq1,/*mate*/(Stage3end_T) NULL,
			  acc,chromosome_iit,resulttype,
			  /*first_read_p*/true,/*npaths_mate*/0,
			  /*mate_chrpos*/0U,quality_shift,
			  sam_read_group_id,invert_first_p,invert_second_p);
      SAM_print_nomapping(fp_nomapping_1,queryseq2,/*mate*/(Stage3end_T) NULL,
			  acc,chromosome_iit,resulttype,
			  /*first_read_p*/false,/*npaths_mate*/0,
			  /*mate_chrpos*/0U,quality_shift,
			  sam_read_group_id,invert_second_p,invert_first_p);
    }

  } else {
    if (failsonlyp == true) {
      /* Unwanted success: skip */

    } else if (resulttype == CONCORDANT_UNIQ || resulttype == CONCORDANT_TRANSLOC) {
      if (resulttype == CONCORDANT_UNIQ) {
	fp = fp_concordant_uniq;
      } else {
	fp = fp_concordant_transloc;
      }

      stage3pairarray = (Stage3pair_T *) Result_array(&npaths,&second_absmq,result);
      /* Stage3pair_eval(stage3pairarray,npaths,maxpaths,queryseq1,queryseq2); */

      stage3pair = stage3pairarray[0];
      hardclip5 = hardclip3 = 0;
      if (clip_overlap_p == true && (overlap = Stage3pair_overlap(stage3pair)) > 0) {
	debug3(printf("overlap = %d\n",overlap));
	hardclip5 = overlap/2;
	hardclip3 = overlap - hardclip5;
      }

      hit5 = Stage3pair_hit5(stage3pair);
      hit3 = Stage3pair_hit3(stage3pair);
      chrpos5 = SAM_compute_chrpos(/*hardclip_low*/&ignore,/*hardclip_high*/&hardclip5,hit5,
				   Stage3end_substring_low(hit5),Shortread_fulllength(queryseq1));
      chrpos3 = SAM_compute_chrpos(/*hardclip_low*/&hardclip3,/*hardclip_high*/&ignore,hit3,
				   Stage3end_substring_low(hit3),Shortread_fulllength(queryseq2));

      /* Note: it is possible that hard clipping can changed a concordant_uniq to an unpaired_uniq,
	 or a concordant_transloc to an unpaired_transloc */

      /* print first end */
      SAM_print(fp,hit5,/*mate*/hit3,acc,/*pathnum*/1,/*npaths*/1,
		Stage3pair_absmq_score(stage3pair),/*second_absmq*/0,
		Stage3pair_mapq_score(stage3pair),chromosome_iit,
		/*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
		Stage3pair_pairlength(stage3pair),/*chrpos*/chrpos5,/*mate_chrpos*/chrpos3,
		hardclip5,hardclip3,resulttype,/*first_read_p*/true,/*npaths_mate*/npaths,
		quality_shift,sam_read_group_id,invert_first_p,invert_second_p,
		merge_samechr_p);

      /* print second end */
      SAM_print(fp,hit3,/*mate*/hit5,acc,/*pathnum*/1,/*npaths*/1,
		Stage3pair_absmq_score(stage3pair),/*second_absmq*/0,
		Stage3pair_mapq_score(stage3pair),chromosome_iit,
		/*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
		Stage3pair_pairlength(stage3pair),/*chrpos*/chrpos3,/*mate_chrpos*/chrpos5,
		hardclip5,hardclip3,resulttype,/*first_read_p*/false,/*npaths_mate*/npaths,
		quality_shift,sam_read_group_id,invert_second_p,invert_first_p,
		merge_samechr_p);
    
    } else if (resulttype == CONCORDANT_MULT) {
      stage3pairarray = (Stage3pair_T *) Result_array(&npaths,&second_absmq,result);

      if (quiet_if_excessive_p && npaths > maxpaths) {
	/* Print as nomapping, but send to fp_concordant_mult */
	SAM_print_nomapping(fp_concordant_mult,queryseq1,/*mate*/(Stage3end_T) NULL,
			    acc,chromosome_iit,resulttype,
			    /*first_read_p*/true,/*npaths_mate*/npaths,
			    /*mate_chrpos*/0U,quality_shift,
			    sam_read_group_id,invert_first_p,invert_second_p);
	SAM_print_nomapping(fp_concordant_mult,queryseq2,/*mate*/(Stage3end_T) NULL,
			    acc,chromosome_iit,resulttype,
			    /*first_read_p*/false,/*npaths_mate*/npaths,
			    /*mate_chrpos*/0U,quality_shift,
			    sam_read_group_id,invert_second_p,invert_first_p);

      } else {
	/* Stage3pair_eval(stage3pairarray,npaths,maxpaths,queryseq1,queryseq2); */
	for (pathnum = 1; pathnum <= npaths && pathnum <= maxpaths; pathnum++) {

	  stage3pair = stage3pairarray[pathnum-1];
	  hardclip5 = hardclip3 = 0;
	  if (clip_overlap_p == true && (overlap = Stage3pair_overlap(stage3pair)) > 0) {
	    debug3(printf("overlap = %d\n",overlap));
	    hardclip5 = overlap/2;
	    hardclip3 = overlap - hardclip5;
	  }

	  hit5 = Stage3pair_hit5(stage3pair);
	  hit3 = Stage3pair_hit3(stage3pair);
	  chrpos5 = SAM_compute_chrpos(/*hardclip_low*/&ignore,/*hardclip_high*/&hardclip5,hit5,
				       Stage3end_substring_low(hit5),Shortread_fulllength(queryseq1));
	  chrpos3 = SAM_compute_chrpos(/*hardclip_low*/&hardclip3,/*hardclip_high*/&ignore,hit3,
				       Stage3end_substring_low(hit3),Shortread_fulllength(queryseq2));

	  /* print first end */
	  SAM_print(fp_concordant_mult,hit5,/*mate*/hit3,acc,pathnum,npaths,
		    Stage3pair_absmq_score(stage3pair),second_absmq,
		    Stage3pair_mapq_score(stage3pair),chromosome_iit,
		    /*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
		    Stage3pair_pairlength(stage3pair),/*chrpos*/chrpos5,/*mate_chrpos*/chrpos3,
		    hardclip5,hardclip3,resulttype,/*first_read_p*/true,/*npaths_mate*/npaths,
		    quality_shift,sam_read_group_id,invert_first_p,invert_second_p,
		    merge_samechr_p);

	  /* print second end */
	  SAM_print(fp_concordant_mult,hit3,/*mate*/hit5,acc,pathnum,npaths,
		    Stage3pair_absmq_score(stage3pair),second_absmq,
		    Stage3pair_mapq_score(stage3pair),chromosome_iit,
		    /*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
		    Stage3pair_pairlength(stage3pair),/*chrpos*/chrpos3,/*mate_chrpos*/chrpos5,
		    hardclip5,hardclip3,resulttype,/*first_read_p*/false,/*npaths_mate*/npaths,
		    quality_shift,sam_read_group_id,invert_second_p,invert_first_p,
		    merge_samechr_p);
	}
      }

    } else if (resulttype == PAIRED_UNIQ) {
      stage3pairarray = (Stage3pair_T *) Result_array(&npaths,&second_absmq,result);
      /* Stage3pair_eval(stage3pairarray,npaths,maxpaths,queryseq1,queryseq2); */

      stage3pair = stage3pairarray[0];
      if ((pairtype = Stage3pair_pairtype(stage3pair)) == PAIRED_INVERSION) {
	fp = fp_paired_uniq_inv;
      } else if (pairtype == PAIRED_SCRAMBLE) {
	fp = fp_paired_uniq_scr;
      } else if (pairtype == PAIRED_TOOLONG) {
	fp = fp_paired_uniq_long;
      } else {
	fprintf(stderr,"Unexpected pairtype %d\n",pairtype);
	abort();
      }

      if (clip_overlap_p == true) {
	overlap = Stage3pair_overlap(stage3pair);
      }
      hardclip5 = hardclip3 = 0;

      hit5 = Stage3pair_hit5(stage3pair);
      hit3 = Stage3pair_hit3(stage3pair);
      chrpos5 = SAM_compute_chrpos(/*hardclip_low*/&ignore,/*hardclip_high*/&hardclip5,hit5,
				   Stage3end_substring_low(hit5),Shortread_fulllength(queryseq1));
      chrpos3 = SAM_compute_chrpos(/*hardclip_low*/&hardclip3,/*hardclip_high*/&ignore,hit3,
				   Stage3end_substring_low(hit3),Shortread_fulllength(queryseq2));

      /* print first end */
      SAM_print(fp,hit5,/*mate*/hit3,acc,/*pathnum*/1,/*npaths*/1,
		Stage3pair_absmq_score(stage3pair),/*second_absmq*/0,
		Stage3pair_mapq_score(stage3pair),chromosome_iit,
		/*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
		Stage3pair_pairlength(stage3pair),/*chrpos*/chrpos5,/*mate_chrpos*/chrpos3,
		hardclip5,hardclip3,resulttype,/*first_read_p*/true,/*npaths_mate*/npaths,
		quality_shift,sam_read_group_id,invert_first_p,invert_second_p,
		merge_samechr_p);

      /* print second end */
      SAM_print(fp,hit3,/*mate*/hit5,acc,/*pathnum*/1,/*npaths*/1,
		Stage3pair_absmq_score(stage3pair),/*second_absmq*/0,
		Stage3pair_mapq_score(stage3pair),chromosome_iit,
		/*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
		Stage3pair_pairlength(stage3pair),/*chrpos*/chrpos3,/*mate_chrpos*/chrpos5,
		hardclip5,hardclip3,resulttype,/*first_read_p*/false,/*npaths_mate*/npaths,
		quality_shift,sam_read_group_id,invert_second_p,invert_first_p,
		merge_samechr_p);

    } else if (resulttype == PAIRED_MULT) {
      stage3pairarray = (Stage3pair_T *) Result_array(&npaths,&second_absmq,result);

      if (quiet_if_excessive_p && npaths > maxpaths) {
	/* Print as nomapping, but send to fp_concordant_mult */
	SAM_print_nomapping(fp_paired_mult,queryseq1,/*mate*/(Stage3end_T) NULL,
			    acc,chromosome_iit,resulttype,
			    /*first_read_p*/true,/*npaths_mate*/npaths,
			    /*mate_chrpos*/0U,quality_shift,
			    sam_read_group_id,invert_first_p,invert_second_p);
	SAM_print_nomapping(fp_paired_mult,queryseq2,/*mate*/(Stage3end_T) NULL,
			    acc,chromosome_iit,resulttype,
			    /*first_read_p*/false,/*npaths_mate*/npaths,
			    /*mate_chrpos*/0U,quality_shift,
			    sam_read_group_id,invert_second_p,invert_first_p);

      } else {
	/* Stage3pair_eval(stage3pairarray,npaths,maxpaths,queryseq1,queryseq2); */
	for (pathnum = 1; pathnum <= npaths && pathnum <= maxpaths; pathnum++) {

	  stage3pair = stage3pairarray[pathnum-1];
	  if (clip_overlap_p == true) {
	    overlap = Stage3pair_overlap(stage3pair);
	  }
	  hardclip5 = hardclip3 = 0;

	  hit5 = Stage3pair_hit5(stage3pair);
	  hit3 = Stage3pair_hit3(stage3pair);
	  chrpos5 = SAM_compute_chrpos(/*hardclip_low*/&ignore,/*hardclip_high*/&hardclip5,hit5,
				       Stage3end_substring_low(hit5),Shortread_fulllength(queryseq1));
	  chrpos3 = SAM_compute_chrpos(/*hardclip_low*/&hardclip3,/*hardclip_high*/&ignore,hit3,
				       Stage3end_substring_low(hit3),Shortread_fulllength(queryseq2));

	  /* print first end */
	  SAM_print(fp_paired_mult,hit5,/*mate*/hit3,acc,pathnum,npaths,
		    Stage3pair_absmq_score(stage3pair),second_absmq,
		    Stage3pair_mapq_score(stage3pair),chromosome_iit,
		    /*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
		    Stage3pair_pairlength(stage3pair),/*chrpos*/chrpos5,/*mate_chrpos*/chrpos3,
		    hardclip5,hardclip3,resulttype,/*first_read_p*/true,/*npaths_mate*/npaths,
		    quality_shift,sam_read_group_id,invert_first_p,invert_second_p,
		    merge_samechr_p);

	  /* print second end */
	  SAM_print(fp_paired_mult,hit3,/*mate*/hit5,acc,pathnum,npaths,
		    Stage3pair_absmq_score(stage3pair),second_absmq,
		    Stage3pair_mapq_score(stage3pair),chromosome_iit,
		    /*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
		    Stage3pair_pairlength(stage3pair),/*chrpos*/chrpos3,/*mate_chrpos*/chrpos5,
		    hardclip5,hardclip3,resulttype,/*first_read_p*/false,/*npaths_mate*/npaths,
		    quality_shift,sam_read_group_id,invert_second_p,invert_first_p,
		    merge_samechr_p);
	}
      }

    } else if (resulttype == UNPAIRED_UNIQ || resulttype == UNPAIRED_TRANSLOC) {
      if (resulttype == UNPAIRED_UNIQ) {
	fp = fp_unpaired_uniq;
      } else {
	fp = fp_unpaired_transloc;
      }

      /* Even though they are not related, we should print mate information in this situation */
      stage3array1 = (Stage3end_T *) Result_array(&npaths1,&second_absmq1,result);
      stage3array2 = (Stage3end_T *) Result_array2(&npaths2,&second_absmq2,result);

      hardclip5 = hardclip3 = 0;

      hit5 = stage3array1[0];
      hit3 = stage3array2[0];
      chrpos5 = SAM_compute_chrpos(/*hardclip_low*/&ignore,/*hardclip_high*/&hardclip5,hit5,
				   Stage3end_substring_low(hit5),Shortread_fulllength(queryseq1));
      chrpos3 = SAM_compute_chrpos(/*hardclip_low*/&hardclip3,/*hardclip_high*/&ignore,hit3,
				   Stage3end_substring_low(hit3),Shortread_fulllength(queryseq2));

      /* print first end */
      /* Stage3end_eval_and_sort(stage3array1,npaths1,maxpaths,queryseq1); */
      SAM_print(fp,hit5,/*mate*/hit3,acc,/*pathnum*/1,/*npaths*/1,
		Stage3end_absmq_score(stage3array1[0]),/*second_absmq*/0,
		Stage3end_mapq_score(stage3array1[0]),chromosome_iit,
		/*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
		/*pairedlength*/0U,/*chrpos*/chrpos5,/*mate_chrpos*/chrpos3,
		/*hardclip5*/0,/*hardclip3*/0,resulttype,/*first_read_p*/true,
		/*npaths_mate*/1,quality_shift,sam_read_group_id,invert_first_p,invert_second_p,
		merge_samechr_p);

      /* print second end */
      /* Stage3end_eval_and_sort(stage3array2,npaths2,maxpaths,queryseq2); */
      SAM_print(fp,hit3,/*mate*/hit5,acc,/*pathnum*/1,/*npaths*/1,
		Stage3end_absmq_score(stage3array2[0]),/*second_absmq*/0,
		Stage3end_mapq_score(stage3array2[0]),chromosome_iit,
		/*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
		/*pairedlength*/0U,/*chrpos*/chrpos3,/*mate_chrpos*/chrpos5,
		/*hardclip5*/0,/*hardclip3*/0,resulttype,/*first_read_p*/false,
		/*npaths_mate*/1,quality_shift,sam_read_group_id,invert_second_p,invert_first_p,
		merge_samechr_p);

    } else if (resulttype == UNPAIRED_MULT) {
      stage3array1 = (Stage3end_T *) Result_array(&npaths1,&second_absmq1,result);
      stage3array2 = (Stage3end_T *) Result_array2(&npaths2,&second_absmq2,result);

#if 0
      /* Do eval and sorting first */
      if (npaths1 == 1) {
	/* Stage3end_eval_and_sort(stage3array1,npaths1,maxpaths,queryseq1); */
      } else if (quiet_if_excessive_p && npaths1 > maxpaths) {
	/* Don't sort */
      } else {
	/* Stage3end_eval_and_sort(stage3array1,npaths1,maxpaths,queryseq1); */
      }
#endif

#if 0
      if (npaths2 == 1) {
	/* Stage3end_eval_and_sort(stage3array2,npaths2,maxpaths,queryseq2); */
      } else if (quiet_if_excessive_p && npaths2 > maxpaths) {
	/* Don't sort */
      } else {
	/* Stage3end_eval_and_sort(stage3array2,npaths2,maxpaths,queryseq2); */
      }
#endif

      /* print first end results */
      if (npaths2 == 0) {
	mate = (Stage3end_T) NULL;
	chrpos3 = 0U;
      } else if (quiet_if_excessive_p && npaths2 > maxpaths) {
	mate = (Stage3end_T) NULL;
	chrpos3 = 0U;
      } else {
	mate = stage3array2[0];
	hardclip3 = 0;
	chrpos3 = SAM_compute_chrpos(/*hardclip_low*/&hardclip3,/*hardclip_high*/&ignore,mate,
				     Stage3end_substring_low(mate),Shortread_fulllength(queryseq2));
      }

      if (npaths1 == 1) {
	stage3 = stage3array1[0];
	hardclip5 = 0;
	chrpos5 = SAM_compute_chrpos(/*hardclip_low*/&ignore,/*hardclip_high*/&hardclip5,stage3,
				     Stage3end_substring_low(stage3),Shortread_fulllength(queryseq1));

	SAM_print(fp_unpaired_mult,stage3,mate,acc,/*pathnum*/1,npaths1,
		  Stage3end_absmq_score(stage3),second_absmq1,
		  Stage3end_mapq_score(stage3),chromosome_iit,
		  /*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
		  /*pairedlength*/0U,/*chrpos*/chrpos5,/*mate_chrpos*/chrpos3,
		  /*hardclip5*/0,/*hardclip3*/0,resulttype,/*first_read_p*/true,
		  /*npaths_mate*/npaths2,quality_shift,sam_read_group_id,
		  invert_first_p,invert_second_p,merge_samechr_p);

      } else if (quiet_if_excessive_p && npaths1 > maxpaths) {
	SAM_print_nomapping(fp_unpaired_mult,queryseq1,mate,acc,chromosome_iit,
			    resulttype,/*first_read_p*/true,/*npaths_mate*/npaths2,
			    /*mate_chrpos*/chrpos3,
			    quality_shift,sam_read_group_id,invert_first_p,invert_second_p);

      } else {
	for (pathnum = 1; pathnum <= npaths1 && pathnum <= maxpaths; pathnum++) {
	  stage3 = stage3array1[pathnum-1];
	  hardclip5 = 0;
	  chrpos5 = SAM_compute_chrpos(/*hardclip_low*/&ignore,/*hardclip_high*/&hardclip5,stage3,
				       Stage3end_substring_low(stage3),Shortread_fulllength(queryseq1));

	  SAM_print(fp_unpaired_mult,stage3,mate,acc,pathnum,npaths1,
		    Stage3end_absmq_score(stage3),second_absmq1,
		    Stage3end_mapq_score(stage3),chromosome_iit,
		    /*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
		    /*pairedlength*/0U,/*chrpos*/chrpos5,/*mate_chrpos*/chrpos3,
		    /*hardclip5*/0,/*hardclip3*/0,resulttype,/*first_read_p*/true,
		    /*npaths_mate*/npaths2,quality_shift,sam_read_group_id,
		    invert_first_p,invert_second_p,merge_samechr_p);
	}
      }
			  
      /* print second end results */
      if (npaths1 == 0) {
	mate = (Stage3end_T) NULL;
	chrpos5 = 0U;
      } else if (quiet_if_excessive_p && npaths1 > maxpaths) {
	mate = (Stage3end_T) NULL;
	chrpos5 = 0U;
      } else {
	mate = stage3array1[0];
	hardclip5 = 0;
	chrpos5 = SAM_compute_chrpos(/*hardclip_low*/&ignore,/*hardclip_high*/&hardclip5,mate,
				     Stage3end_substring_low(mate),Shortread_fulllength(queryseq1));
      }

      if (npaths2 == 1) {
	stage3 = stage3array2[0];
	hardclip3 = 0;
	chrpos3 = SAM_compute_chrpos(/*hardclip_low*/&hardclip3,/*hardclip_high*/&ignore,stage3,
				     Stage3end_substring_low(stage3),Shortread_fulllength(queryseq2));

	SAM_print(fp_unpaired_mult,stage3,mate,acc,/*pathnum*/1,npaths2,
		  Stage3end_absmq_score(stage3),second_absmq2,
		  Stage3end_mapq_score(stage3),chromosome_iit,
		  /*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
		  /*pairedlength*/0U,/*chrpos*/chrpos3,/*mate_chrpos*/chrpos5,
		  /*hardclip5*/0,/*hardclip3*/0,resulttype,/*first_read_p*/false,
		  /*npaths_mate*/npaths1,quality_shift,sam_read_group_id,
		  invert_second_p,invert_first_p,merge_samechr_p);

      } else if (quiet_if_excessive_p && npaths2 > maxpaths) {
	SAM_print_nomapping(fp_unpaired_mult,queryseq2,mate,acc,chromosome_iit,
			    resulttype,/*first_read_p*/false,/*npaths_mate*/npaths1,
			    /*mate_chrpos*/chrpos5,
			    quality_shift,sam_read_group_id,invert_second_p,invert_first_p);

      } else {
	for (pathnum = 1; pathnum <= npaths2 && pathnum <= maxpaths; pathnum++) {
	  stage3 = stage3array2[pathnum-1];
	  hardclip3 = 0;
	  chrpos3 = SAM_compute_chrpos(/*hardclip_low*/&hardclip3,/*hardclip_high*/&ignore,stage3,
				       Stage3end_substring_low(stage3),Shortread_fulllength(queryseq2));

	  SAM_print(fp_unpaired_mult,stage3,mate,acc,pathnum,npaths2,
		    Stage3end_absmq_score(stage3),second_absmq2,
		    Stage3end_mapq_score(stage3),chromosome_iit,
		    /*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
		    /*pairedlength*/0U,/*chrpos*/chrpos3,/*mate_chrpos*/chrpos5,
		    /*hardclip5*/0,/*hardclip3*/0,resulttype,/*first_read_p*/false,
		    /*npaths_mate*/npaths1,quality_shift,sam_read_group_id,
		    invert_second_p,invert_first_p,merge_samechr_p);
	}
      }

    } else {
      if (resulttype == HALFMAPPING_UNIQ) {
	fp = fp_halfmapping_uniq;
      } else if (resulttype == HALFMAPPING_TRANSLOC) {
	fp = fp_halfmapping_transloc;
      } else if (resulttype == HALFMAPPING_MULT) {
	fp = fp_halfmapping_mult;
      } else {
	abort();
      }

      stage3array1 = (Stage3end_T *) Result_array(&npaths1,&second_absmq1,result);
      stage3array2 = (Stage3end_T *) Result_array2(&npaths2,&second_absmq2,result);

#if 0
      /* Do eval and sorting first */
      if (npaths1 == 0) {
	/* Nothing to sort */
      } else if (npaths1 == 1) {
	/* Stage3end_eval_and_sort(stage3array1,npaths1,maxpaths,queryseq1); */
      } else if (quiet_if_excessive_p && npaths1 > maxpaths) {
	/* Don't sort */
      } else {
	/* Stage3end_eval_and_sort(stage3array1,npaths1,maxpaths,queryseq1); */
      }
#endif

#if 0
      if (npaths2 == 0) {
	/* Nothing to sort */
      } else if (npaths2 == 1) {
	/* Stage3end_eval_and_sort(stage3array2,npaths2,maxpaths,queryseq2); */
      } else if (quiet_if_excessive_p && npaths2 > maxpaths) {
	/* Don't sort */
      } else {
	/* Stage3end_eval_and_sort(stage3array2,npaths2,maxpaths,queryseq2); */
      }
#endif


      /* print first end results */
      if (npaths2 == 0) {
	mate = (Stage3end_T) NULL;
	chrpos3 = 0U;
      } else if (quiet_if_excessive_p && npaths2 > maxpaths) {
	mate = (Stage3end_T) NULL;
	chrpos3 = 0U;
      } else {
	mate = stage3array2[0];
	hardclip3 = 0;
	chrpos3 = SAM_compute_chrpos(/*hardclip_low*/&hardclip3,/*hardclip_high*/&ignore,mate,
				     Stage3end_substring_low(mate),Shortread_fulllength(queryseq2));
      }

      if (npaths1 == 0) {
	/* mate should be non-NULL here */
	SAM_print_nomapping(fp,queryseq1,mate,acc,chromosome_iit,resulttype,
			    /*first_read_p*/true,/*npaths_mate*/npaths2,
			    /*mate_chrpos*/chrpos3,
			    quality_shift,sam_read_group_id,invert_first_p,invert_second_p);

      } else if (npaths1 == 1) {
	/* mate should be NULL here */

	stage3 = stage3array1[0];
	hardclip5 = 0;
	chrpos5 = SAM_compute_chrpos(/*hardclip_low*/&ignore,/*hardclip_high*/&hardclip5,stage3,
				     Stage3end_substring_low(stage3),Shortread_fulllength(queryseq1));

	SAM_print(fp,stage3,mate,acc,/*pathnum*/1,npaths1,
		  Stage3end_absmq_score(stage3),second_absmq1,
		  Stage3end_mapq_score(stage3),chromosome_iit,
		  /*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
		  /*pairedlength*/0U,/*chrpos*/chrpos5,/*mate_chrpos*/chrpos3,
		  /*hardclip5*/0,/*hardclip3*/0,resulttype,/*first_read_p*/true,
		  /*npaths_mate*/npaths2,quality_shift,sam_read_group_id,
		  invert_first_p,invert_second_p,merge_samechr_p);

      } else if (quiet_if_excessive_p && npaths1 > maxpaths) {
	/* mate should be NULL here */
	SAM_print_nomapping(fp,queryseq1,mate,acc,chromosome_iit,resulttype,
			    /*first_read_p*/true,/*npaths_mate*/npaths2,
			    /*mate_chrpos*/chrpos3,
			    quality_shift,sam_read_group_id,invert_first_p,invert_second_p);

      } else {
	/* mate should be NULL here */
	for (pathnum = 1; pathnum <= npaths1 && pathnum <= maxpaths; pathnum++) {
	  stage3 = stage3array1[pathnum-1];
	  hardclip5 = 0;
	  chrpos5 = SAM_compute_chrpos(/*hardclip_low*/&ignore,/*hardclip_high*/&hardclip5,stage3,
				       Stage3end_substring_low(stage3),Shortread_fulllength(queryseq1));

	  SAM_print(fp,stage3,mate,acc,pathnum,npaths1,
		    Stage3end_absmq_score(stage3),second_absmq1,
		    Stage3end_mapq_score(stage3),chromosome_iit,
		    /*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
		    /*pairedlength*/0U,/*chrpos*/chrpos5,/*mate_chrpos*/chrpos3,
		    /*hardclip5*/0,/*hardclip3*/0,resulttype,/*first_read_p*/true,
		    /*npaths_mate*/npaths2,quality_shift,sam_read_group_id,
		    invert_first_p,invert_second_p,merge_samechr_p);
	}
      }
			  
      /* print second end results */
      if (npaths1 == 0) {
	mate = (Stage3end_T) NULL;
	chrpos5 = 0U;
      } else if (quiet_if_excessive_p && npaths1 > maxpaths) {
	mate = (Stage3end_T) NULL;
	chrpos5 = 0U;
      } else {
	mate = stage3array1[0];
	hardclip5 = 0;
	chrpos5 = SAM_compute_chrpos(/*hardclip_low*/&ignore,/*hardclip_high*/&hardclip5,mate,
				     Stage3end_substring_low(mate),Shortread_fulllength(queryseq1));
      }

      if (npaths2 == 0) {
	/* mate should be non-NULL here */
	SAM_print_nomapping(fp,queryseq2,mate,acc,chromosome_iit,resulttype,
			    /*first_read_p*/false,/*npaths_mate*/npaths1,
			    /*mate_chrpos*/chrpos5,
			    quality_shift,sam_read_group_id,invert_second_p,invert_first_p);

      } else if (npaths2 == 1) {
	/* mate should be NULL here */

	stage3 = stage3array2[0];
	hardclip3 = 0;
	chrpos3 = SAM_compute_chrpos(/*hardclip_low*/&hardclip3,/*hardclip_high*/&ignore,stage3,
				     Stage3end_substring_low(stage3),Shortread_fulllength(queryseq2));

	SAM_print(fp,stage3,mate,acc,/*pathnum*/1,npaths2,
		  Stage3end_absmq_score(stage3),second_absmq2,
		  Stage3end_mapq_score(stage3),chromosome_iit,
		  /*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
		  /*pairedlength*/0U,/*chrpos*/chrpos3,/*mate_chrpos*/chrpos5,
		  /*hardclip5*/0,/*hardclip3*/0,resulttype,/*first_read_p*/false,
		  /*npaths_mate*/npaths1,quality_shift,sam_read_group_id,
		  invert_second_p,invert_first_p,merge_samechr_p);

      } else if (quiet_if_excessive_p && npaths2 > maxpaths) {
	/* mate should be NULL here */
	SAM_print_nomapping(fp,queryseq2,mate,acc,chromosome_iit,resulttype,
			    /*first_read_p*/false,/*npaths_mate*/npaths1,
			    /*mate_chrpos*/chrpos5,
			    quality_shift,sam_read_group_id,invert_second_p,invert_first_p);

      } else {
	/* mate should be NULL here */
	for (pathnum = 1; pathnum <= npaths2 && pathnum <= maxpaths; pathnum++) {
	  stage3 = stage3array2[pathnum-1];
	  hardclip3 = 0;
	  chrpos3 = SAM_compute_chrpos(/*hardclip_low*/&hardclip3,/*hardclip_high*/&ignore,stage3,
				       Stage3end_substring_low(stage3),Shortread_fulllength(queryseq2));

	  SAM_print(fp,stage3,mate,acc,pathnum,npaths2,
		    Stage3end_absmq_score(stage3),second_absmq2,
		    Stage3end_mapq_score(stage3),chromosome_iit,
		    /*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
		    /*pairedlength*/0U,/*chrpos*/chrpos3,/*mate_chrpos*/chrpos5,
		    /*hardclip5*/0,/*hardclip3*/0,resulttype,/*first_read_p*/false,
		    /*npaths_mate*/npaths1,quality_shift,sam_read_group_id,
		    invert_second_p,invert_first_p,merge_samechr_p);
	}
      }

    }
  }

  return;
}



