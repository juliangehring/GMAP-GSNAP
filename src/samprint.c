static char rcsid[] = "$Id: samprint.c 37254 2011-03-28 16:34:08Z twu $";
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
#include "stage1hr.h"		/* for MAX_QUERYLENGTH */


#define SANGER_ILLUMINA_DIFF 31

/* BAM appears to truncate the H information on the ends of a cigar */
/* Also, this provides the information needed for getting term information */
/* #define HARDCLIP_ON_TERMINAL 1 */



#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


/* print_md_string */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif


unsigned int
SAM_compute_flag (Substring_T substring, Stage3_T mate, Resulttype_T resulttype,
		  bool first_read_p, int pathnum, int npaths, int npaths_mate,
		  bool invertp, bool invert_mate_p) {
  unsigned int flag = 0U;

  debug(printf("Resulttype: %s\n",Resulttype_string(resulttype)));

  if (npaths == 0) {
    debug(printf("npaths = 0, so QUERY_UNMAPPED %d\n",QUERY_UNMAPPED));
    flag |= QUERY_UNMAPPED;
  } else if (Substring_plusp(substring) == invertp) {
    debug(printf("plusp %d and invertp %d, so QUERY_MINUSP %d\n",
		 Substring_plusp(substring),invertp,QUERY_MINUSP));
    flag |= QUERY_MINUSP;
  }

  if (resulttype == SINGLEEND_NOMAPPING || resulttype == SINGLEEND_UNIQ || resulttype == SINGLEEND_MULT) {
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

    } else if (mate == NULL) {
      /* Unpaired; no mate.  Not clear if should be MATE_UNMAPPED. */

    } else if (resulttype == CONCORDANT_UNIQ || resulttype == CONCORDANT_MULT) {
      /* Can distinguish concordant mappings by presence of insert length */
      debug(printf("PAIRED_MAPPING %d\n",PAIRED_MAPPING));
      flag |= PAIRED_MAPPING;
      if (Stage3_chrnum(mate) == 0) {
	/* Splice without a direction */

      } else if (Stage3_plusp(mate) == invert_mate_p) {
	debug(printf("MATE_MINUSP %d\n",MATE_MINUSP));
	flag |= MATE_MINUSP;
      }

    } else if (resulttype == PAIRED_UNIQ || resulttype == PAIRED_MULT) {
      /* Note: We are counting PAIRED_UNIQ and PAIRED_MULT as "paired" mappings.
	 However, we are no longer counting UNPAIRED_UNIQ as a "paired" mapping. */
      debug(printf("PAIRED_MAPPING %d\n",PAIRED_MAPPING));
      flag |= PAIRED_MAPPING;
      if (Stage3_chrnum(mate) == 0) {
	/* Splice without a direction */

      } else if (Stage3_plusp(mate) == invert_mate_p) {
	debug(printf("MATE_MINUSP %d\n",MATE_MINUSP));
	flag |= MATE_MINUSP;
      }

    } else {
      if (Stage3_chrnum(mate) == 0) {
	/* Splice without a direction */
      } else if (Stage3_plusp(mate) == invert_mate_p) {
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


static Genomicpos_T
print_chromosomal_pos (FILE *fp, Stage3_T this, Stage3_T mate, IIT_T chromosome_iit, bool matep) {
  Substring_T substring;
  Genomicpos_T chrpos;
  Chrnum_T chrnum;
  char *chr;
  bool allocp;

  if (this == NULL) {
    fprintf(fp,"\t*\t0");
    return 0U;

  } else if (matep == true && mate != NULL && Stage3_chrnum(mate) > 0 && 
	     Stage3_chrnum(this) == Stage3_chrnum(mate)) {
    fprintf(fp,"\t=");

  } else if ((chrnum = Stage3_chrnum(this)) == 0) {
    /* Interchromosomal splice */
    if (matep == true) {
      /* No way to express two chromosomes for mate */
      fprintf(fp,"\t*\t0");
      return 0U;
    } else {
      fprintf(stderr,"Trying to print interchrosomal splice in one line\n");
      abort();
    }

  } else {
    chr = IIT_label(chromosome_iit,chrnum,&allocp);
    fprintf(fp,"\t%s",chr);
    if (allocp == true) {
      FREE(chr);
    }
  }

#if 0
  hittype = Stage3_hittype(this);
  if (hittype == SPLICE) {
    if (Stage3_distance(this) > 0U) {
      sensep = (Stage3_sensedir(this) == SENSE_FORWARD);
      if (sensep == Stage3_plusp(this)) {
	substring = /* donor */ Stage3_substring1(this);
      } else {
	substring = /* acceptor */ Stage3_substring2(this);
      }
    } else if (/*acceptor*/Stage3_substring2(this) == NULL) {
      substring = Stage3_substring1(this);
    } else if (/*donor*/Stage3_substring1(this) == NULL) {
      substring = Stage3_substring2(this);
    } else {
      /* Distant splicing.  Should not have reached here.  Should print each end separately. */
      fprintf(stderr,"Distant splicing of %u.  Should not have reached here\n",Stage3_distance(this));
      fprintf(stderr,"Genomicstart is %u, Genomicend is %u\n",Stage3_genomicstart(this),Stage3_genomicend(this));
      abort();
    }
  } else if (Stage3_substring2(this) == NULL) {
    substring = Stage3_substring1(this);
  } else {
    if (Stage3_plusp(this) == true) {
      if (Substring_alignstart_trim(Stage3_substring1(this)) < Substring_alignstart_trim(Stage3_substring2(this))) {
	substring = Stage3_substring1(this);
      } else {
	substring = Stage3_substring2(this);
      }
    } else {
      if (Substring_alignend_trim(Stage3_substring1(this)) < Substring_alignend_trim(Stage3_substring2(this))) {
	substring = Stage3_substring1(this);
      } else {
	substring = Stage3_substring2(this);
      }
    }
  }
#else
  substring = Stage3_substring_low(this);
#endif

  if (Stage3_plusp(this) == true) {
    chrpos = Substring_alignstart_trim(substring) - Substring_chroffset(substring);
  } else {
    chrpos = Substring_alignend_trim(substring) - Substring_chroffset(substring);
  }
  /* Add 1 to report in 1-based coordinates */
  fprintf(fp,"\t%u",chrpos+1U);

  return chrpos;
}


static Genomicpos_T
print_substring_pos (FILE *fp, Substring_T substring, IIT_T chromosome_iit) {
  Genomicpos_T chrpos;
  char *chr;
  bool allocp;

  chr = IIT_label(chromosome_iit,Substring_chrnum(substring),&allocp);
  fprintf(fp,"\t%s",chr);
  if (allocp == true) {
    FREE(chr);
  }

  if (Substring_plusp(substring) == true) {
    chrpos = Substring_alignstart_trim(substring) - Substring_chroffset(substring);
    /* Add 1 to report in 1-based coordinates */
    fprintf(fp,"\t%u",chrpos+1U);
  } else {
    chrpos = Substring_alignend_trim(substring) - Substring_chroffset(substring);
    /* Add 1 to report in 1-based coordinates */
    fprintf(fp,"\t%u",chrpos+1U);
  }

  return chrpos;
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
SAM_print_nomapping (FILE *fp, Shortread_T queryseq, Stage3_T mate, char *acc,
		     IIT_T chromosome_iit, Resulttype_T resulttype, bool translocationp,
		     bool first_read_p, int npaths_mate, Shortread_T queryseq_mate,
		     int quality_shift, char *sam_read_group_id, bool invertp,
		     bool invert_mate_p) {
  unsigned int flag;
  char *barcode;

  /* 1. QNAME */
  fprintf(fp,"%s",acc);
  
  /* 2. FLAG */
  flag = SAM_compute_flag(/*substring*/NULL,mate,resulttype,first_read_p,
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
  print_chromosomal_pos(fp,mate,/*other*/NULL,chromosome_iit,/*matep*/true);

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
  if ((barcode = Shortread_barcode(queryseq)) != NULL) {
    fprintf(fp,"\t");
    fprintf(fp,"XB:Z:%s",barcode);
  }

  /* 12. TAGS: XT string (translocation) */
  if (translocationp == true) {
    fprintf(fp,"\t");
    fprintf(fp,"XT:i:1");
  }

  fprintf(fp,"\n");
  return;
}


static int
print_md_string (FILE *fp, int matchlength, char *genomicfwd, int querylength,
		 int hardclip_low, int hardclip_high, bool lastp) {
  int i;

#if 0
  if (genomicfwd == NULL) {
    printf("Entering md_string with matchlength %d\n",matchlength);
  } else {
    printf("Entering md_string with matchlength %d, querylength %d, hardclip_low %d, hardclip_high %d: %s\n",
	   matchlength,querylength,hardclip_low,hardclip_high,genomicfwd);
  }
#endif

  for (i = hardclip_low; i < querylength - hardclip_high; i++) {
    if (isupper(genomicfwd[i])) {
      matchlength++;
    } else {
      if (matchlength > 0 || i == hardclip_low) {
	fprintf(fp,"%d",matchlength);
      }
      fprintf(fp,"%c",toupper(genomicfwd[i]));
      matchlength = 0;
    }
  }

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
print_single (FILE *fp, Stage3_T this, Stage3_T mate, char *acc, int pathnum, int npaths,
	      int mapq_score, IIT_T chromosome_iit, Shortread_T queryseq, int pairedlength,
	      Resulttype_T resulttype, bool translocationp, bool first_read_p,
	      int npaths_mate, Shortread_T queryseq_mate, int quality_shift,
	      char *sam_read_group_id, bool invertp, bool invert_mate_p) {
  unsigned int flag = 0U;
  Substring_T substring;
  int querylength, substring_start, substring_length;
  int hardclip_low = 0, hardclip_high = 0;
  char *genomicfwd, *genomicdir;
  char *barcode;
  Genomicpos_T chrpos, mate_chrpos;

  querylength = Shortread_fulllength(queryseq);
  substring = Stage3_substring1(this);

  /* 1. QNAME */
  fprintf(fp,"%s",acc);

  /* 2. FLAG */
  flag = SAM_compute_flag(substring,mate,resulttype,first_read_p,
			  pathnum,npaths,npaths_mate,invertp,invert_mate_p);
  fprintf(fp,"\t%u",flag);

  /* 3. RNAME: chr */
  /* 4. POS: chrpos */
  chrpos = print_chromosomal_pos(fp,this,mate,chromosome_iit,/*matep*/false);

  /* 5. MAPQ: Mapping quality */
  fprintf(fp,"\t%d",mapq_score);

  /* 6. CIGAR */
  fprintf(fp,"\t");
  if (Stage3_plusp(this) == true) {
    if (Substring_querystart(substring) > 0) {
#ifdef HARDCLIP_ON_TERMINAL
      if (Substring_left_endtype(substring) == TERM) {
	hardclip_low = Substring_querystart(substring);
	fprintf(fp,"%dH",hardclip_low);
      } else {
	fprintf(fp,"%dS",Substring_querystart(substring));
      }
#else
      fprintf(fp,"%dS",Substring_querystart(substring));
#endif
    }
  } else {
    if (Substring_queryend(substring) < querylength) {
#ifdef HARDCLIP_ON_TERMINAL
      if (Substring_right_endtype(substring) == TERM) {
	hardclip_low = querylength - Substring_queryend(substring);
	fprintf(fp,"%dH",hardclip_low);
      } else {
	fprintf(fp,"%dS",querylength - Substring_queryend(substring));
      }
#else
      fprintf(fp,"%dS",querylength - Substring_queryend(substring));
#endif
    }
  }

  fprintf(fp,"%dM",Substring_match_length(substring));

  if (Stage3_plusp(this) == true) {
    if (Substring_queryend(substring) < querylength) {
#ifdef HARDCLIP_ON_TERMINAL
      if (Substring_right_endtype(substring) == TERM) {
	hardclip_high = querylength - Substring_queryend(substring);
	fprintf(fp,"%dH",hardclip_high);
      } else {
	fprintf(fp,"%dS",querylength - Substring_queryend(substring));
      }
#else
      fprintf(fp,"%dS",querylength - Substring_queryend(substring));
#endif
    }
  } else {
    if (Substring_querystart(substring) > 0) {
#ifdef HARDCLIP_ON_TERMINAL
      if (Substring_left_endtype(substring) == TERM) {
	hardclip_high = Substring_querystart(substring);
	fprintf(fp,"%dH",hardclip_high);
      } else {
	fprintf(fp,"%dS",Substring_querystart(substring));
      }
#else
      fprintf(fp,"%dS",Substring_querystart(substring));
#endif
    }
  }


  /* 7. MRNM: Mate chr */
  /* 8. MPOS: Mate chrpos */
  mate_chrpos = print_chromosomal_pos(fp,mate,this,chromosome_iit,/*matep*/true);

  /* 9. ISIZE: Insert size */
#ifdef PAIRED_ZERO_ISIZE
  if (resulttype == CONCORDANT_UNIQ || resulttype == CONCORDANT_MULT) {
    if (Stage3_plusp(this) == invertp) {
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
  } else {
    fprintf(fp,"\t%d",-pairedlength);
  }
#endif



  /* 10. SEQ: queryseq and 11. QUAL: quality scores */
  /* Queryseq has already been inverted, so just measure plusp relative to its current state */
  fprintf(fp,"\t");
  if (Stage3_plusp(this) == true) {
    Shortread_print_chopped(fp,queryseq,hardclip_low,hardclip_high);
    fprintf(fp,"\t");
    Shortread_print_quality(fp,queryseq,hardclip_low,hardclip_high,quality_shift,/*show_chopped_p*/false);
  } else {
    Shortread_print_chopped_revcomp(fp,queryseq,hardclip_high,hardclip_low);
    fprintf(fp,"\t");
    Shortread_print_quality_revcomp(fp,queryseq,hardclip_high,hardclip_low,quality_shift,/*show_chopped_p*/false);
  }

  /* 12. TAGS: RG */
  if (sam_read_group_id != NULL) {
    fprintf(fp,"\t");
    fprintf(fp,"RG:Z:%s",sam_read_group_id);
  }

  /* 12. TAGS: XB */
  if ((barcode = Shortread_barcode(queryseq)) != NULL) {
    fprintf(fp,"\t");
    fprintf(fp,"XB:Z:%s",barcode);
  }

  /* 12. TAGS: NM */
  fprintf(fp,"\t");
  fprintf(fp,"NM:i:%d",Stage3_nmismatches_refdiff(this));

  /* 12. TAGS: MD */
  fprintf(fp,"\t");

#ifdef MD_INCLUDES_SOFTCLIP
  substring_start = 0;
  substring_length = querylength;
#else
  substring_start = Substring_querystart(substring);
  substring_length = Substring_match_length(substring);
#endif

  fprintf(fp,"MD:Z:");  
  if ((genomicdir = Substring_genomic_refdiff(substring)) == NULL) {
    fprintf(fp,"%d",querylength);

  } else if (Stage3_plusp(this) == true) {
    print_md_string(fp,/*matchlength*/0,&(genomicdir[substring_start]),substring_length,
		    hardclip_low,hardclip_high,/*lastp*/true);

  } else {
    genomicfwd = (char *) CALLOC(querylength+1,sizeof(char));
    make_complement_buffered(genomicfwd,&(genomicdir[substring_start]),substring_length);
    print_md_string(fp,/*matchlength*/0,genomicfwd,substring_length,
		    hardclip_low,hardclip_high,/*lastp*/true);
    FREE(genomicfwd);
  }
  
  /* 12. TAGS: SM */
  fprintf(fp,"\t");
  fprintf(fp,"SM:i:%d",mapq_score);

  /* 12. TAGS: XT string (translocation) */
  if (translocationp == true) {
    fprintf(fp,"\t");
    fprintf(fp,"XT:i:1");
  }

  fprintf(fp,"\n");
  return;
}


static void
print_insertion (FILE *fp, Stage3_T this, Stage3_T mate, char *acc, int pathnum, int npaths,
		 int mapq_score, IIT_T chromosome_iit, Shortread_T queryseq, int pairedlength,
		 Resulttype_T resulttype, bool translocationp, bool first_read_p,
		 int npaths_mate, Shortread_T queryseq_mate, int quality_shift,
		 char *sam_read_group_id, bool invertp, bool invert_mate_p) {
  unsigned int flag = 0U;
  Substring_T substring1, substring2;
  int querylength;
  char *genomicfwd, *genomicdir;
  int substring1_start, substring2_start, substring1_length, substring2_length, nindels, matchlength;
  char *barcode;
  Genomicpos_T chrpos, mate_chrpos;

  querylength = Shortread_fulllength(queryseq);

  substring1 = Stage3_substring1(this);
  substring2 = Stage3_substring2(this);

  /* 1. QNAME */
  fprintf(fp,"%s",acc);

  /* 2. FLAG */
  flag = SAM_compute_flag(substring1,mate,resulttype,first_read_p,
			  pathnum,npaths,npaths_mate,invertp,invert_mate_p);
  fprintf(fp,"\t%u",flag);

  /* 3. RNAME: chr */
  /* 4. POS: chrpos */
  chrpos = print_chromosomal_pos(fp,this,mate,chromosome_iit,/*matep*/false);

  /* 5. MAPQ: Mapping quality */
  fprintf(fp,"\t%d",mapq_score);

  /* 6. CIGAR */
  fprintf(fp,"\t");

  nindels = Stage3_nindels(this);
  if (Stage3_plusp(this) == true) {
    if (Substring_querystart(substring1) > 0) {
      fprintf(fp,"%dS",Substring_querystart(substring1));
    }
    fprintf(fp,"%dM",Substring_match_length(substring1));
    fprintf(fp,"%dI",nindels);
    fprintf(fp,"%dM",Substring_match_length(substring2));
    if (Substring_queryend(substring2) < querylength) {
      fprintf(fp,"%dS",querylength - Substring_queryend(substring2));
    }

  } else {
    if (Substring_queryend(substring2) < querylength) {
      fprintf(fp,"%dS",querylength - Substring_queryend(substring2));
    }
    fprintf(fp,"%dM",Substring_match_length(substring2));
    fprintf(fp,"%dI",nindels);
    fprintf(fp,"%dM",Substring_match_length(substring1));
    if (Substring_querystart(substring1) > 0) {
      fprintf(fp,"%dS",Substring_querystart(substring1));
    }
  }

  /* 7. MRNM: Mate chr */
  /* 8. MPOS: Mate chrpos */
  mate_chrpos = print_chromosomal_pos(fp,mate,this,chromosome_iit,/*matep*/true);

  /* 9. ISIZE: Insert size */
#ifdef PAIRED_ZERO_ISIZE
  if (resulttype == CONCORDANT_UNIQ || resulttype == CONCORDANT_MULT) {
    if (Stage3_plusp(this) == invertp) {
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
  } else {
    fprintf(fp,"\t%d",-pairedlength);
  }
#endif

  /* 10. SEQ: queryseq and 11. QUAL: quality scores */
  /* Queryseq has already been inverted, so just measure plusp relative to its current state */
  fprintf(fp,"\t");
  if (Stage3_plusp(this) == true) {
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
  if ((barcode = Shortread_barcode(queryseq)) != NULL) {
    fprintf(fp,"\t");
    fprintf(fp,"XB:Z:%s",barcode);
  }

  /* 12. TAGS: NM */
  fprintf(fp,"\t");
  fprintf(fp,"NM:i:%d",Stage3_nmismatches_refdiff(this));
  
  /* 12. TAGS: MD */
  fprintf(fp,"\t");
  fprintf(fp,"MD:Z:");


#ifdef MD_INCLUDES_SOFTCLIP
  substring1_start = Substring_querystart_orig(substring1);
  substring1_length = Substring_match_length_orig(substring1);
  substring2_start = Substring_querystart_orig(substring2);
  substring2_length = Substring_match_length_orig(substring2);
#else
  substring1_start = Substring_querystart(substring1);
  substring1_length = Substring_match_length(substring1);
  substring2_start = Substring_querystart(substring2);
  substring2_length = Substring_match_length(substring2);
#endif

  if (Stage3_plusp(this) == true) {
    genomicfwd = Substring_genomic_refdiff(substring1);
    matchlength = print_md_string(fp,/*matchlength*/0,&(genomicfwd[substring1_start]),substring1_length,
				  /*hardclip_low*/0,/*hardclip_high*/0,/*lastp*/false);
    genomicfwd = Substring_genomic_refdiff(substring2);
    print_md_string(fp,matchlength,&(genomicfwd[substring2_start]),substring2_length,
		    /*hardclip_low*/0,/*hardclip_high*/0,/*lastp*/true);
  } else {
    genomicfwd = (char *) CALLOC(substring2_length+1,sizeof(char));
    genomicdir = Substring_genomic_refdiff(substring2);
    make_complement_buffered(genomicfwd,&(genomicdir[substring2_start]),substring2_length);
    matchlength = print_md_string(fp,/*matchlength*/0,genomicfwd,substring2_length,
				  /*hardclip_low*/0,/*hardclip_high*/0,/*lastp*/false);
    FREE(genomicfwd);

    genomicfwd = (char *) CALLOC(substring1_length+1,sizeof(char));
    genomicdir = Substring_genomic_refdiff(substring1);
    make_complement_buffered(genomicfwd,&(genomicdir[substring1_start]),substring1_length);
    print_md_string(fp,matchlength,genomicfwd,substring1_length,
		    /*hardclip_low*/0,/*hardclip_high*/0,/*lastp*/true);
    FREE(genomicfwd);
  }

  /* 12. TAGS: SM */
  fprintf(fp,"\t");
  fprintf(fp,"SM:i:%d",mapq_score);

  /* 12. TAGS: XT string (translocation) */
  if (translocationp == true) {
    fprintf(fp,"\t");
    fprintf(fp,"XT:i:1");
  }

  fprintf(fp,"\n");
  return;
}


static void
print_deletion (FILE *fp, Stage3_T this, Stage3_T mate, char *acc, int pathnum, int npaths,
		int mapq_score, IIT_T chromosome_iit, Shortread_T queryseq, int pairedlength,
		Resulttype_T resulttype, bool translocationp, bool first_read_p,
		int npaths_mate, Shortread_T queryseq_mate, int quality_shift,
		char *sam_read_group_id, bool invertp, bool invert_mate_p) {
  unsigned int flag = 0U;
  Substring_T substring1, substring2;
  int querylength;
  char *genomicfwd, *genomicdir;
  int substring1_start, substring2_start, substring1_length, substring2_length, nindels;
  char *barcode;
  Genomicpos_T chrpos, mate_chrpos;

  querylength = Shortread_fulllength(queryseq);

  substring1 = Stage3_substring1(this);
  substring2 = Stage3_substring2(this);

  /* 1. QNAME */
  fprintf(fp,"%s",acc);

  /* 2. FLAG */
  flag = SAM_compute_flag(substring1,mate,resulttype,first_read_p,
			  pathnum,npaths,npaths_mate,invertp,invert_mate_p);
  fprintf(fp,"\t%u",flag);

  /* 3. RNAME: chr */
  /* 4. POS: chrpos */
  chrpos = print_chromosomal_pos(fp,this,mate,chromosome_iit,/*matep*/false);

  /* 5. MAPQ: Mapping quality */
  fprintf(fp,"\t%d",mapq_score);

  /* 6. CIGAR */
  fprintf(fp,"\t");

  nindels = Stage3_nindels(this); /* nindels is positive */
  if (Stage3_plusp(this) == true) {
    if (Substring_querystart(substring1) > 0) {
      fprintf(fp,"%dS",Substring_querystart(substring1));
    }
    fprintf(fp,"%dM",Substring_match_length(substring1));
    fprintf(fp,"%dD",nindels);
    fprintf(fp,"%dM",Substring_match_length(substring2));
    if (Substring_queryend(substring2) < querylength) {
      fprintf(fp,"%dS",querylength - Substring_queryend(substring2));
    }

  } else {
    if (Substring_queryend(substring2) < querylength) {
      fprintf(fp,"%dS",querylength - Substring_queryend(substring2));
    }
    fprintf(fp,"%dM",Substring_match_length(substring2));
    fprintf(fp,"%dD",nindels);
    fprintf(fp,"%dM",Substring_match_length(substring1));
    if (Substring_querystart(substring1) > 0) {
      fprintf(fp,"%dS",Substring_querystart(substring1));
    }
  }

  /* 7. MRNM: Mate chr */
  /* 8. MPOS: Mate chrpos */
  mate_chrpos = print_chromosomal_pos(fp,mate,this,chromosome_iit,/*matep*/true);

  /* 9. ISIZE: Insert size */
#ifdef PAIRED_ZERO_ISIZE
  if (resulttype == CONCORDANT_UNIQ || resulttype == CONCORDANT_MULT) {
    if (Stage3_plusp(this) == invertp) {
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
  } else {
    fprintf(fp,"\t%d",-pairedlength);
  }
#endif


  /* 10. SEQ: queryseq and 11. QUAL: quality scores */
  /* Queryseq has already been inverted, so just measure plusp relative to its current state */
  fprintf(fp,"\t");
  if (Stage3_plusp(this) == true) {
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
  if ((barcode = Shortread_barcode(queryseq)) != NULL) {
    fprintf(fp,"\t");
    fprintf(fp,"XB:Z:%s",barcode);
  }

  /* 12. TAGS: NM */
  fprintf(fp,"\t");
  fprintf(fp,"NM:i:%d",Stage3_nmismatches_refdiff(this));
  
  /* 12. TAGS: MD */
  fprintf(fp,"\t");
  fprintf(fp,"MD:Z:");

#ifdef MD_INCLUDES_SOFTCLIP
  substring1_start = Substring_querystart_orig(substring1);
  substring1_length = Substring_match_length_orig(substring1); /* match_length + querystart */
  substring2_start = Substring_querystart_orig(substring2);
  substring2_length = Substring_match_length_orig(substring2); /* querylength - queryend + match_length */
#else
  substring1_start = Substring_querystart(substring1);
  substring1_length = Substring_match_length(substring1);
  substring2_start = Substring_querystart(substring2);
  substring2_length = Substring_match_length(substring2);
#endif

  if (Stage3_plusp(this) == true) {
    genomicfwd = Substring_genomic_refdiff(substring1);
    print_md_string(fp,/*matchlength*/0,&(genomicfwd[substring1_start]),substring1_length,
		    /*hardclip_low*/0,/*hardclip_high*/0,/*lastp*/true);

    /* Deletion string: Potential problem if followed by a mismatch, but can be resolved by looking at CIGAR string */
    fprintf(fp,"^%s",Stage3_deletion_string(this));

    genomicfwd = Substring_genomic_refdiff(substring2);
    print_md_string(fp,/*matchlength*/0,&(genomicfwd[substring2_start]),substring2_length,
		    /*hardclip_low*/0,/*hardclip_high*/0,/*lastp*/true);
  } else {
    genomicfwd = (char *) CALLOC(substring2_length+1,sizeof(char));
    genomicdir = Substring_genomic_refdiff(substring2);
    make_complement_buffered(genomicfwd,&(genomicdir[substring2_start]),substring2_length);
    print_md_string(fp,/*matchlength*/0,genomicfwd,substring2_length,
		    /*hardclip_low*/0,/*hardclip_high*/0,/*lastp*/true);
    FREE(genomicfwd);

    /* Deletion string: Potential problem if followed by a mismatch, but can be resolved by looking at CIGAR string */
    genomicfwd = (char *) CALLOC(nindels+1,sizeof(char));
    make_complement_buffered(genomicfwd,Stage3_deletion_string(this),nindels);
    fprintf(fp,"^%s",genomicfwd);
    FREE(genomicfwd);

    genomicfwd = (char *) CALLOC(substring1_length+1,sizeof(char));
    genomicdir = Substring_genomic_refdiff(substring1);
    make_complement_buffered(genomicfwd,&(genomicdir[substring1_start]),substring1_length);
    print_md_string(fp,/*matchlength*/0,genomicfwd,substring1_length,
		    /*hardclip_low*/0,/*hardclip_high*/0,/*lastp*/true);
    FREE(genomicfwd);
  }

  /* 12. TAGS: SM */
  fprintf(fp,"\t");
  fprintf(fp,"SM:i:%d",mapq_score);

  /* 12. TAGS: XT string (translocation) */
  if (translocationp == true) {
    fprintf(fp,"\t");
    fprintf(fp,"XT:i:1");
  }

  fprintf(fp,"\n");
  return;
}


static void
print_halfdonor (FILE *fp, Substring_T donor, Stage3_T this, Stage3_T mate,
		 char *acc, int pathnum, int npaths,
		 int mapq_score, IIT_T chromosome_iit, Shortread_T queryseq, int pairedlength,
		 Resulttype_T resulttype, bool translocationp, bool first_read_p,
		 int npaths_mate, Shortread_T queryseq_mate, int quality_shift,
		 char *sam_read_group_id, bool invertp, bool invert_mate_p,
		 bool use_hardclip_p) {
  unsigned int flag = 0U;
  int querylength;
  int clip_low = 0, clip_high = 0;
  bool sensep;
  char *genomicfwd, *genomicdir;
  int substring_start, substring_length;
  char *barcode;
  Genomicpos_T chrpos, mate_chrpos;

  querylength = Shortread_fulllength(queryseq);

  /* 1. QNAME */
  fprintf(fp,"%s",acc);

  /* 2. FLAG */
  flag = SAM_compute_flag(donor,mate,resulttype,first_read_p,
			  pathnum,npaths,npaths_mate,invertp,invert_mate_p);
  fprintf(fp,"\t%u",flag);

  /* 3. RNAME: chr */
  /* 4. POS: chrpos */
#if 0
  print_chromosomal_pos(fp,this,mate,chromosome_iit,/*matep*/false);
#else
  /* Need to do this, because print_chromosomal_pos cannot handle translocations */
  chrpos = print_substring_pos(fp,donor,chromosome_iit);
#endif

  /* 5. MAPQ: Mapping quality */
  fprintf(fp,"\t%d",mapq_score);

  /* 6. CIGAR */
  fprintf(fp,"\t");
  sensep = Substring_chimera_sensep(donor);

  if (sensep == Substring_plusp(donor)) {
    if (Substring_plusp(donor) == true) {
      /* sensep true */
      if (Substring_querystart(donor) > 0) {
	fprintf(fp,"%dS",Substring_querystart(donor));
      }
    } else {
      /* sensep false */
      if (Substring_queryend(donor) < querylength) {
	fprintf(fp,"%dS",querylength - Substring_queryend(donor));
      }
    }
    fprintf(fp,"%dM",Substring_match_length(donor));
    if (sensep == true) {
      clip_high = querylength - Substring_chimera_pos(donor);
    } else {
      clip_high = Substring_chimera_pos(donor);
    }
    if (use_hardclip_p == true) {
      fprintf(fp,"%dH",clip_high); /* intron */
    } else {
      fprintf(fp,"%dS",clip_high); /* intron */
    }


  } else { /* sensep != Substring_plusp(donor) */
    if (sensep == true) {
      clip_low = querylength - Substring_chimera_pos(donor);
    } else {
      clip_low = Substring_chimera_pos(donor);
    }
    if (use_hardclip_p == true) {
      fprintf(fp,"%dH",clip_low); /* intron */
    } else {
      fprintf(fp,"%dS",clip_low); /* intron */
    }


    fprintf(fp,"%dM",Substring_match_length(donor));
    if (Substring_plusp(donor) == true) {
      /* sensep false */
      if (Substring_queryend(donor) < querylength) {
	fprintf(fp,"%dS",querylength - Substring_queryend(donor));
      }
    } else {
      /* sensep true */
      if (Substring_querystart(donor) > 0) {
	fprintf(fp,"%dS",Substring_querystart(donor));
      }
    }
  }

  /* 7. MRNM: Mate chr */
  /* 8. MPOS: Mate chrpos */
  mate_chrpos = print_chromosomal_pos(fp,mate,this,chromosome_iit,/*matep*/true);

  /* 9. ISIZE: Insert size */
#ifdef PAIRED_ZERO_ISIZE
  if (resulttype == CONCORDANT_UNIQ || resulttype == CONCORDANT_MULT) {
    if (Stage3_plusp(this) == invertp) {
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
  } else {
    fprintf(fp,"\t%d",-pairedlength);
  }
#endif


  /* 10. SEQ: queryseq and 11. QUAL: quality scores */
  /* Queryseq has already been inverted, so just measure plusp relative to its current state */
  fprintf(fp,"\t");
  if (use_hardclip_p == true) {
    if (Substring_plusp(donor) == true) {
      Shortread_print_chopped(fp,queryseq,clip_low,clip_high);
      fprintf(fp,"\t");
      Shortread_print_quality(fp,queryseq,clip_low,clip_high,
			      quality_shift,/*show_chopped_p*/false);
    } else {
      Shortread_print_chopped_revcomp(fp,queryseq,clip_low,clip_high);
      fprintf(fp,"\t");
      Shortread_print_quality_revcomp(fp,queryseq,clip_low,clip_high,
				      quality_shift,/*show_chopped_p*/false);
    }
  } else {
    if (Substring_plusp(donor) == true) {
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
  }


  /* 12. TAGS: RG */
  if (sam_read_group_id != NULL) {
    fprintf(fp,"\t");
    fprintf(fp,"RG:Z:%s",sam_read_group_id);
  }

  /* 12. TAGS: XB */
  if ((barcode = Shortread_barcode(queryseq)) != NULL) {
    fprintf(fp,"\t");
    fprintf(fp,"XB:Z:%s",barcode);
  }

  /* 12. TAGS: NM */
  fprintf(fp,"\t");
  fprintf(fp,"NM:i:%d",Substring_nmismatches_refdiff(donor));
  
  /* 12. TAGS: MD */
  fprintf(fp,"\t");
  fprintf(fp,"MD:Z:");

#ifdef MD_INCLUDES_SOFTCLIP
  substring_start = Substring_querystart_orig(donor);
  substring_length = Substring_match_length_orig(donor);
#else
  substring_start = Substring_querystart(donor);
  substring_length = Substring_match_length(donor);
#endif

  if (use_hardclip_p == false) {
    genomicdir = Substring_genomic_refdiff(donor);
    if (Substring_plusp(donor) == true) {
      print_md_string(fp,/*matchlength*/0,&(genomicdir[substring_start]),substring_length,
		      /*hardclip_low*/0,/*hardclip_high*/0,/*lastp*/true);
    } else {
      genomicfwd = (char *) CALLOC(querylength+1,sizeof(char));
      make_complement_buffered(genomicfwd,&(genomicdir[substring_start]),substring_length);
      print_md_string(fp,/*matchlength*/0,genomicfwd,substring_length,
		      /*hardclip_low*/0,/*hardclip_high*/0,/*lastp*/true);
      FREE(genomicfwd);
    }

  } else if (sensep == true) {
    if (Substring_plusp(donor) == true) {
      genomicfwd = Substring_genomic_refdiff(donor);
      print_md_string(fp,/*matchlength*/0,&(genomicfwd[substring_start]),substring_length,
		      /*hardclip_low*/0,/*hardclip_high*/0,/*lastp*/true);
    } else {
      genomicfwd = (char *) CALLOC(substring_length+1,sizeof(char));
      genomicdir = Substring_genomic_refdiff(donor);
      make_complement_buffered(genomicfwd,&(genomicdir[substring_start]),substring_length);
      print_md_string(fp,/*matchlength*/0,genomicfwd,substring_length,
		      /*hardclip_low*/0,/*hardclip_high*/0,/*lastp*/true);
      FREE(genomicfwd);
    }

  } else {			/* sensep == false */
    if (Substring_plusp(donor) == true) {
      genomicfwd = Substring_genomic_refdiff(donor);
      print_md_string(fp,/*matchlength*/0,&(genomicfwd[substring_start]),substring_length,
		      /*hardclip_low*/0,/*hardclip_high*/0,/*lastp*/true);
    } else {
      genomicfwd = (char *) CALLOC(substring_length+1,sizeof(char));
      genomicdir = Substring_genomic_refdiff(donor);
      make_complement_buffered(genomicfwd,&(genomicdir[substring_start]),substring_length);
      print_md_string(fp,/*matchlength*/0,genomicfwd,substring_length,
		      /*hardclip_low*/0,/*hardclip_high*/0,/*lastp*/true);
      FREE(genomicfwd);
    }
  }

  /* 12. TAGS: SM */
  fprintf(fp,"\t");
  fprintf(fp,"SM:i:%d",mapq_score);

  /* 12. TAGS: XT string (translocation) */
  if (translocationp == true) {
    fprintf(fp,"\t");
    fprintf(fp,"XT:i:1");
  }

  /* 12. TAGS: XS */
  fprintf(fp,"\t");
  if (sensep == Substring_plusp(donor)) {
    fprintf(fp,"XS:A:+");
  } else {
    fprintf(fp,"XS:A:-");
  }

  fprintf(fp,"\n");
  return;
}


static void
print_halfacceptor (FILE *fp, Substring_T acceptor, Stage3_T this, Stage3_T mate,
		    char *acc, int pathnum, int npaths,
		    int mapq_score, IIT_T chromosome_iit, Shortread_T queryseq, int pairedlength,
		    Resulttype_T resulttype, bool translocationp, bool first_read_p,
		    int npaths_mate, Shortread_T queryseq_mate, int quality_shift,
		    char *sam_read_group_id, bool invertp, bool invert_mate_p,
		    bool use_hardclip_p) {
  unsigned int flag = 0U;
  int querylength;
  int clip_low = 0, clip_high = 0;
  bool sensep;
  char *genomicfwd, *genomicdir;
  int substring_start, substring_length;
  char *barcode;
  Genomicpos_T chrpos, mate_chrpos;

  querylength = Shortread_fulllength(queryseq);

  /* 1. QNAME */
  fprintf(fp,"%s",acc);

  /* 2. FLAG */
  flag = SAM_compute_flag(acceptor,mate,resulttype,first_read_p,
			  pathnum,npaths,npaths_mate,invertp,invert_mate_p);
  fprintf(fp,"\t%u",flag);

  /* 3. RNAME: chr */
  /* 4. POS: chrpos */
#if 0
  print_chromosomal_pos(fp,this,mate,chromosome_iit,/*matep*/false);
#else
  /* Need to do this, because print_chromosomal_pos cannot handle translocations */
  chrpos = print_substring_pos(fp,acceptor,chromosome_iit);
#endif

  /* 5. MAPQ: Mapping quality */
  fprintf(fp,"\t%d",mapq_score);

  /* 6. CIGAR */
  fprintf(fp,"\t");
  sensep = Substring_chimera_sensep(acceptor);

  if (sensep != Substring_plusp(acceptor)) {
    if (Substring_plusp(acceptor) == true) {
      /* sensep false */
      if (Substring_querystart(acceptor) > 0) {
	fprintf(fp,"%dS",Substring_querystart(acceptor));
      }
    } else {
      /* sensep true */
      if (Substring_queryend(acceptor) < querylength) {
	fprintf(fp,"%dS",querylength - Substring_queryend(acceptor));
      }
    }
    fprintf(fp,"%dM",Substring_match_length(acceptor));

    if (sensep == true) {
      clip_high = Substring_chimera_pos(acceptor);
    } else {
      clip_high = querylength - Substring_chimera_pos(acceptor);
    }
    if (use_hardclip_p == true) {
      fprintf(fp,"%dH",clip_high);  /* intron */
    } else {
      fprintf(fp,"%dS",clip_high);  /* intron */
    }


  } else { /*  sensep == Substring_plusp(acceptor) */
    if (sensep == true) {
      clip_low = Substring_chimera_pos(acceptor);
    } else {
      clip_low = querylength - Substring_chimera_pos(acceptor);
    }
    if (use_hardclip_p == true) {
      fprintf(fp,"%dH",clip_low);  /* intron */
    } else {
      fprintf(fp,"%dS",clip_low);  /* intron */
    }

    fprintf(fp,"%dM",Substring_match_length(acceptor));
    if (Substring_plusp(acceptor) == true) {
      /* sensep true */
      if (Substring_queryend(acceptor) < querylength) {
	fprintf(fp,"%dS",querylength - Substring_queryend(acceptor));
      }
    } else {
      /* sensep false */
      if (Substring_querystart(acceptor) > 0) {
	fprintf(fp,"%dS",Substring_querystart(acceptor));
      }
    }
  }


  /* 7. MRNM: Mate chr */
  /* 8. MPOS: Mate chrpos */
  mate_chrpos = print_chromosomal_pos(fp,mate,this,chromosome_iit,/*matep*/true);

  /* 9. ISIZE: Insert size */
#ifdef PAIRED_ZERO_ISIZE
  if (resulttype == CONCORDANT_UNIQ || resulttype == CONCORDANT_MULT) {
    if (Stage3_plusp(this) == invertp) {
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
  } else {
    fprintf(fp,"\t%d",-pairedlength);
  }
#endif


  /* 10. SEQ: queryseq and 11. QUAL: quality scores */
  /* Queryseq has already been inverted, so just measure plusp relative to its current state */
  fprintf(fp,"\t");
  if (use_hardclip_p == true) {
    if (Substring_plusp(acceptor) == true) {
      Shortread_print_chopped(fp,queryseq,clip_low,clip_high);
      fprintf(fp,"\t");
      Shortread_print_quality(fp,queryseq,clip_low,clip_high,
			      quality_shift,/*show_chopped_p*/false);
    } else {
      Shortread_print_chopped_revcomp(fp,queryseq,clip_low,clip_high);
      fprintf(fp,"\t");
      Shortread_print_quality_revcomp(fp,queryseq,clip_low,clip_high,
				      quality_shift,/*show_chopped_p*/false);
    }
  } else {
    if (Substring_plusp(acceptor) == true) {
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
  }


  /* 12. TAGS: RG */
  if (sam_read_group_id != NULL) {
    fprintf(fp,"\t");
    fprintf(fp,"RG:Z:%s",sam_read_group_id);
  }

  /* 12. TAGS: XB */
  if ((barcode = Shortread_barcode(queryseq)) != NULL) {
    fprintf(fp,"\t");
    fprintf(fp,"XB:Z:%s",barcode);
  }

  /* 12. TAGS: NM */
  fprintf(fp,"\t");
  fprintf(fp,"NM:i:%d",Substring_nmismatches_refdiff(acceptor));
  
  /* 12. TAGS: MD */
  fprintf(fp,"\t");
  fprintf(fp,"MD:Z:");

#ifdef MD_INCLUDES_SOFTCLIP
  substring_start = Substring_querystart_orig(acceptor);
  substring_length = Substring_match_length_orig(acceptor);
#else
  substring_start = Substring_querystart(acceptor);
  substring_length = Substring_match_length(acceptor);
#endif


  if (use_hardclip_p == false) {
    genomicdir = Substring_genomic_refdiff(acceptor);
    if (Substring_plusp(acceptor) == true) {
      print_md_string(fp,/*matchlength*/0,&(genomicdir[substring_start]),substring_length,
		      /*hardclip_low*/0,/*hardclip_high*/0,/*lastp*/true);
    } else {
      genomicfwd = (char *) CALLOC(querylength+1,sizeof(char));
      make_complement_buffered(genomicfwd,&(genomicdir[substring_start]),substring_length);
      print_md_string(fp,/*matchlength*/0,genomicfwd,substring_length,
		      /*hardclip_low*/0,/*hardclip_high*/0,/*lastp*/true);
      FREE(genomicfwd);
    }

  } else if (sensep == false) {
    if (Substring_plusp(acceptor) == true) {
      genomicfwd = Substring_genomic_refdiff(acceptor);
      print_md_string(fp,/*matchlength*/0,&(genomicfwd[substring_start]),substring_length,
		      /*hardclip_low*/0,/*hardclip_high*/0,/*lastp*/true);
    } else {
      genomicfwd = (char *) CALLOC(substring_length+1,sizeof(char));
      genomicdir = Substring_genomic_refdiff(acceptor);
      make_complement_buffered(genomicfwd,&(genomicdir[substring_start]),substring_length);
      print_md_string(fp,/*matchlength*/0,genomicfwd,substring_length,
		      /*hardclip_low*/0,/*hardclip_high*/0,/*lastp*/true);
      FREE(genomicfwd);
    }

  } else {			/* sensep true */
    if (Substring_plusp(acceptor) == true) {
      genomicfwd = Substring_genomic_refdiff(acceptor);
      print_md_string(fp,/*matchlength*/0,&(genomicfwd[substring_start]),substring_length,
		      /*hardclip_low*/0,/*hardclip_high*/0,/*lastp*/true);
    } else {
      genomicfwd = (char *) CALLOC(substring_length+1,sizeof(char));
      genomicdir = Substring_genomic_refdiff(acceptor);
      make_complement_buffered(genomicfwd,&(genomicdir[substring_start]),substring_length);
      print_md_string(fp,/*matchlength*/0,genomicfwd,substring_length,
		      /*hardclip_low*/0,/*hardclip_high*/0,/*lastp*/true);
      FREE(genomicfwd);
    }
  }

  /* 12. TAGS: SM */
  fprintf(fp,"\t");
  fprintf(fp,"SM:i:%d",mapq_score);

  /* 12. TAGS: XT string (translocation) */
  if (translocationp == true) {
    fprintf(fp,"\t");
    fprintf(fp,"XT:i:1");
  }

  /* 12. TAGS: XS */
  fprintf(fp,"\t");
  if (sensep == Substring_plusp(acceptor)) {
    fprintf(fp,"XS:A:+");
  } else {
    fprintf(fp,"XS:A:-");
  }

  fprintf(fp,"\n");
  return;
}



static void
print_localsplice (FILE *fp, Stage3_T chimera, Stage3_T mate, char *acc, int pathnum, int npaths,
		   int mapq_score, IIT_T chromosome_iit, Shortread_T queryseq, int pairedlength,
		   Resulttype_T resulttype, bool translocationp, bool first_read_p,
		   int npaths_mate, Shortread_T queryseq_mate, int quality_shift,
		   char *sam_read_group_id, bool invertp, bool invert_mate_p) {
  unsigned int flag = 0U;
  Substring_T substring1, substring2;
  int querylength;
  bool sensep;
  char *genomicfwd, *genomicdir;
  int substring1_start, substring2_start, substring1_length, substring2_length, matchlength;
  char *barcode;
  Genomicpos_T chrpos, mate_chrpos;

  querylength = Shortread_fulllength(queryseq);

  sensep = (Stage3_sensedir(chimera) == SENSE_FORWARD);


  /* 1. QNAME */
  fprintf(fp,"%s",acc);

  /* 2. FLAG */
  flag = SAM_compute_flag(Stage3_substring1(chimera),mate,resulttype,first_read_p,
			  pathnum,npaths,npaths_mate,invertp,invert_mate_p);
  fprintf(fp,"\t%u",flag);

  /* 3. RNAME: chr */
  /* 4. POS: chrpos */
  chrpos = print_chromosomal_pos(fp,chimera,mate,chromosome_iit,/*matep*/false);

  /* 5. MAPQ: Mapping quality */
  fprintf(fp,"\t%d",mapq_score);

  /* 6. CIGAR */
  fprintf(fp,"\t");
  if (sensep == Stage3_plusp(chimera)) {
    substring1 = /* donor */ Stage3_substringD(chimera);
    substring2 = /* acceptor */ Stage3_substringA(chimera);
  } else {
    substring1 = /* acceptor */ Stage3_substringA(chimera);
    substring2 = /* donor */ Stage3_substringD(chimera);
  }
  if (Stage3_plusp(chimera) == true) {
    if (Substring_querystart(substring1) > 0) {
      fprintf(fp,"%dS",Substring_querystart(substring1));
    }
  } else {
    if (Substring_queryend(substring1) < querylength) {
      fprintf(fp,"%dS",querylength - Substring_queryend(substring1));
    }
  }
    
  fprintf(fp,"%dM",Substring_match_length(substring1));
  fprintf(fp,"%uN",Stage3_distance(chimera));
  fprintf(fp,"%dM",Substring_match_length(substring2));

  if (Stage3_plusp(chimera) == true) {
    if (Substring_queryend(substring2) < querylength) {
      fprintf(fp,"%dS",querylength - Substring_queryend(substring2));
    }
  } else {
    if (Substring_querystart(substring2) > 0) {
      fprintf(fp,"%dS",Substring_querystart(substring2));
    }
  }


  /* 7. MRNM: Mate chr */
  /* 8. MPOS: Mate chrpos */
  mate_chrpos = print_chromosomal_pos(fp,mate,chimera,chromosome_iit,/*matep*/true);

  /* 9. ISIZE: Insert size */
#ifdef PAIRED_ZERO_ISIZE
  if (resulttype == CONCORDANT_UNIQ || resulttype == CONCORDANT_MULT) {
    if (Stage3_plusp(chimera) == invertp) {
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
  } else {
    fprintf(fp,"\t%d",-pairedlength);
  }
#endif

  /* 10. SEQ: queryseq and 11. QUAL: quality scores */
  /* Queryseq has already been inverted, so just measure plusp relative to its current state */
  fprintf(fp,"\t");
  if (Stage3_plusp(chimera) == true) {
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
  if ((barcode = Shortread_barcode(queryseq)) != NULL) {
    fprintf(fp,"\t");
    fprintf(fp,"XB:Z:%s",barcode);
  }

  /* 12. TAGS: NM */
  fprintf(fp,"\t");
  fprintf(fp,"NM:i:%d",Stage3_nmismatches_refdiff(chimera));
  
  /* 12. TAGS: MD */
  fprintf(fp,"\t");
  fprintf(fp,"MD:Z:");


#ifdef MD_INCLUDES_SOFTCLIP
  /* These lengths include part that is trimmed */
  substring1_start = Substring_querystart_orig(substring1);
  substring1_length = Substring_match_length_orig(substring1);
  substring2_start = Substring_querystart_orig(substring2);
  substring2_length = Substring_match_length_orig(substring2);
#else
  substring1_start = Substring_querystart(substring1);
  substring1_length = Substring_match_length(substring1);
  substring2_start = Substring_querystart(substring2);
  substring2_length = Substring_match_length(substring2);
#endif

  if (Stage3_plusp(chimera) == true) {
    genomicfwd = Substring_genomic_refdiff(substring1);
    debug2(printf("genomicfwd: %s  start %d\n",genomicfwd,substring1_start));
    matchlength = print_md_string(fp,/*matchlength*/0,&(genomicfwd[substring1_start]),substring1_length,
				  /*hardclip_low*/0,/*hardclip_high*/0,/*lastp*/false);
    
#if 0
    /* Intron: Gets skipped in MD string */
    fprintf(fp,"^");
#endif

    genomicfwd = Substring_genomic_refdiff(substring2);
    debug2(printf("genomicfwd: %s  start %d\n",genomicfwd,substring2_start));
    print_md_string(fp,matchlength,&(genomicfwd[substring2_start]),substring2_length,
		    /*hardclip_low*/0,/*hardclip_high*/0,/*lastp*/true);

  } else {
    genomicfwd = (char *) CALLOC(substring1_length+1,sizeof(char));
    genomicdir = Substring_genomic_refdiff(substring1);
    debug2(printf("genomicdir: %s  start %d\n",genomicfwd,substring1_start));
    make_complement_buffered(genomicfwd,&(genomicdir[substring1_start]),substring1_length);
    matchlength = print_md_string(fp,/*matchlength*/0,genomicfwd,substring1_length,
				  /*hardclip_low*/0,/*hardclip_high*/0,/*lastp*/false);
    FREE(genomicfwd);

#if 0
    /* Intron: Gets skipped in MD string */
    fprintf(fp,"^");
#endif

    genomicfwd = (char *) CALLOC(substring2_length+1,sizeof(char));
    genomicdir = Substring_genomic_refdiff(substring2);
    debug2(printf("genomicdir: %s  start %d\n",genomicfwd,substring1_start));
    make_complement_buffered(genomicfwd,&(genomicdir[substring2_start]),substring2_length);
    print_md_string(fp,matchlength,genomicfwd,substring2_length,
		    /*hardclip_low*/0,/*hardclip_high*/0,/*lastp*/true);
    FREE(genomicfwd);
  }

  /* 12. TAGS: SM */
  fprintf(fp,"\t");
  fprintf(fp,"SM:i:%d",mapq_score);

  /* 12. TAGS: XT string (translocation) */
  if (translocationp == true) {
    fprintf(fp,"\t");
    fprintf(fp,"XT:i:1");
  }

  /* 12. TAGS: XS */
  fprintf(fp,"\t");
  if (sensep == Stage3_plusp(chimera)) {
    fprintf(fp,"XS:A:+");
  } else {
    fprintf(fp,"XS:A:-");
  }

  fprintf(fp,"\n");
  return;
}


static void
print_shortexon (FILE *fp, Stage3_T shortexon, Stage3_T mate, char *acc, int pathnum, int npaths,
		 int mapq_score, IIT_T chromosome_iit, Shortread_T queryseq, int pairedlength,
		 Resulttype_T resulttype, bool translocationp, bool first_read_p,
		 int npaths_mate, Shortread_T queryseq_mate, int quality_shift,
		 char *sam_read_group_id, bool invertp, bool invert_mate_p) {
  unsigned int flag = 0U;
  /* substring1 is low coordinate on genome, substring2 is high */
  Substring_T substring1, substring2, substringM;
  Genomicpos_T distance1, distance2;
  int querylength;
  bool sensep;
  char *genomicfwd, *genomicdir;
  int substring1_start, substring2_start, substringM_start,
    substring1_length, substring2_length, substringM_length, matchlength;
  char *barcode;
  Genomicpos_T chrpos, mate_chrpos;

  querylength = Shortread_fulllength(queryseq);

  sensep = (Stage3_sensedir(shortexon) == SENSE_FORWARD);


  /* 1. QNAME */
  fprintf(fp,"%s",acc);

  /* 2. FLAG */
  flag = SAM_compute_flag(Stage3_substring1(shortexon),mate,resulttype,first_read_p,
			  pathnum,npaths,npaths_mate,invertp,invert_mate_p);
  fprintf(fp,"\t%u",flag);

  /* 3. RNAME: chr */
  /* 4. POS: chrpos */
  chrpos = print_chromosomal_pos(fp,shortexon,mate,chromosome_iit,/*matep*/false);

  /* 5. MAPQ: Mapping quality */
  fprintf(fp,"\t%d",mapq_score);

  /* 6. CIGAR */
  fprintf(fp,"\t");
  substringM = Stage3_substring1(shortexon);

  if (sensep == Stage3_plusp(shortexon)) {
    substring1 = /* donor */ Stage3_substringD(shortexon);
    distance1 = Stage3_shortexon_acceptor_distance(shortexon);
    distance2 = Stage3_shortexon_donor_distance(shortexon);
    substring2 = /* acceptor */ Stage3_substringA(shortexon);
  } else {
    substring1 = /* acceptor */ Stage3_substringA(shortexon);
    distance1 = Stage3_shortexon_donor_distance(shortexon);
    distance2 = Stage3_shortexon_acceptor_distance(shortexon);
    substring2 = /* donor */ Stage3_substringD(shortexon);
  }

  if (substring1 == NULL) {
    if (Stage3_plusp(shortexon) == true) {
      fprintf(fp,"%dS",Substring_querystart(substringM));
    } else {
      fprintf(fp,"%dS",querylength - Substring_queryend(substringM));
    }

  } else if (Stage3_plusp(shortexon) == true) {
    if (Substring_querystart(substring1) > 0) {
      fprintf(fp,"%dS",Substring_querystart(substring1));
    }
    fprintf(fp,"%dM",Substring_match_length(substring1));
    fprintf(fp,"%uN",distance1);

  } else {
    if (Substring_queryend(substring1) < querylength) {
      fprintf(fp,"%dS",querylength - Substring_queryend(substring1));
    }
    fprintf(fp,"%dM",Substring_match_length(substring1));
    fprintf(fp,"%uN",distance1);

  }
    
  fprintf(fp,"%dM",Substring_match_length(substringM));



  if (substring2 == NULL) {
    if (Stage3_plusp(shortexon) == true) {
      fprintf(fp,"%dS",querylength - Substring_queryend(substringM));
    } else {
      fprintf(fp,"%dS",Substring_querystart(substringM));
    }

  } else if (Stage3_plusp(shortexon) == true) {
    fprintf(fp,"%uN",distance2);
    fprintf(fp,"%dM",Substring_match_length(substring2));
    if (Substring_queryend(substring2) < querylength) {
      fprintf(fp,"%dS",querylength - Substring_queryend(substring2));
    }
  } else {
    fprintf(fp,"%uN",distance2);
    fprintf(fp,"%dM",Substring_match_length(substring2));
    if (Substring_querystart(substring2) > 0) {
      fprintf(fp,"%dS",Substring_querystart(substring2));
    }
  }


  /* 7. MRNM: Mate chr */
  /* 8. MPOS: Mate chrpos */
  mate_chrpos = print_chromosomal_pos(fp,mate,shortexon,chromosome_iit,/*matep*/true);

  /* 9. ISIZE: Insert size */
#ifdef PAIRED_ZERO_ISIZE
  if (resulttype == CONCORDANT_UNIQ || resulttype == CONCORDANT_MULT) {
    if (Stage3_plusp(shortexon) == invertp) {
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
  } else {
    fprintf(fp,"\t%d",-pairedlength);
  }
#endif


  /* 10. SEQ: queryseq and 11. QUAL: quality scores */
  /* Queryseq has already been inverted, so just measure plusp relative to its current state */
  fprintf(fp,"\t");
  if (Stage3_plusp(shortexon) == true) {
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
  if ((barcode = Shortread_barcode(queryseq)) != NULL) {
    fprintf(fp,"\t");
    fprintf(fp,"XB:Z:%s",barcode);
  }

  /* 12. TAGS: NM */
  fprintf(fp,"\t");
  fprintf(fp,"NM:i:%d",Stage3_nmismatches_refdiff(shortexon));
  
  /* 12. TAGS: MD */
  fprintf(fp,"\t");
  fprintf(fp,"MD:Z:");

#ifdef MD_INCLUDES_SOFTCLIP
  substringM_start = Substring_querystart_orig(substringM);
  substringM_length = Substring_match_length_orig(substringM);
#else
  substringM_start = Substring_querystart(substringM);
  substringM_length = Substring_match_length(substringM);
#endif

  if (substring1 == NULL) {
    substring1_length = 0;
#ifdef MD_INCLUDES_SOFTCLIP
    if (Stage3_plusp(shortexon) == true) {
      substringM_start = 0;
      substringM_length += Substring_querystart(substringM);
    } else {
      substringM_length += querylength - Substring_queryend(substringM);
    }
#endif
  } else {
#ifdef MD_INCLUDES_SOFTCLIP
    substring1_start = Substring_querystart_orig(substring1);
    substring1_length = Substring_match_length_orig(substring1);
#else
    substring1_start = Substring_querystart(substring1);
    substring1_length = Substring_match_length(substring1);
#endif
  }
  if (substring2 == NULL) {
    substring2_length = 0;
#ifdef MD_INCLUDES_SOFTCLIP
    if (Stage3_plusp(shortexon) == true) {
      substringM_length += querylength - Substring_queryend(substringM);
    } else {
      substringM_start = 0;
      substringM_length += Substring_querystart(substringM);
    }
#endif
  } else {
#ifdef MD_INCLUDES_SOFTCLIP
    substring2_start = Substring_querystart_orig(substring2);
    substring2_length = Substring_match_length_orig(substring2);
#else
    substring2_start = Substring_querystart(substring2);
    substring2_length = Substring_match_length(substring2);
#endif
  }

  if (Stage3_plusp(shortexon) == true) {

    if (substring1 == NULL) {
      matchlength = 0;
    } else {
      genomicfwd = Substring_genomic_refdiff(substring1);
      debug2(printf("genomicfwd: %s  start %d\n",genomicfwd,substring1_start));
      matchlength = print_md_string(fp,/*matchlength*/0,&(genomicfwd[substring1_start]),substring1_length,
				    /*hardclip_low*/0,/*hardclip_high*/0,/*lastp*/false);
    }

#if 0
    /* Intron 1: Gets skipped in MD string */
    fprintf(fp,"^");
#endif

    genomicfwd = Substring_genomic_refdiff(substringM);
    debug2(printf("genomicfwd: %s  start %d\n",genomicfwd,substringM_start));
    matchlength = print_md_string(fp,matchlength,&(genomicfwd[substringM_start]),substringM_length,
				  /*hardclip_low*/0,/*hardclip_high*/0,/*lastp*/false);

#if 0
    /* Intron 2: Gets skipped in MD string */
    fprintf(fp,"^");
#endif

    if (substring2 == NULL) {
      print_md_string(fp,matchlength,/*genomicfwd*/NULL,/*substring2_length*/0,
		      /*hardclip_low*/0,/*hardclip_high*/0,/*lastp*/true);
    } else {
      genomicfwd = Substring_genomic_refdiff(substring2);
      debug2(printf("genomicfwd: %s  start %d\n",genomicfwd,substring2_start));
      print_md_string(fp,matchlength,&(genomicfwd[substring2_start]),substring2_length,
		      /*hardclip_low*/0,/*hardclip_high*/0,/*lastp*/true);
    }

  } else {

    if (substring1 == NULL) {
      matchlength = 0;
    } else {
      genomicfwd = (char *) CALLOC(substring1_length+1,sizeof(char));
      genomicdir = Substring_genomic_refdiff(substring1);
      debug2(printf("substring1 genomicdir: %s  start %d length %d\n",genomicdir,substring1_start,substring1_length));
      make_complement_buffered(genomicfwd,&(genomicdir[substring1_start]),substring1_length);
      matchlength = print_md_string(fp,/*matchlength*/0,genomicfwd,substring1_length,
				    /*hardclip_low*/0,/*hardclip_high*/0,/*lastp*/false);
      FREE(genomicfwd);
    }

#if 0
    /* Intron 1: Gets skipped in MD string */
    fprintf(fp,"^");
#endif

    genomicfwd = (char *) CALLOC(substringM_length+1,sizeof(char));
    genomicdir = Substring_genomic_refdiff(substringM);
    debug2(printf("substringM genomicdir: %s  start %d length %d\n",genomicdir,substringM_start,substringM_length));
    make_complement_buffered(genomicfwd,&(genomicdir[substringM_start]),substringM_length);
    matchlength = print_md_string(fp,matchlength,genomicfwd,substringM_length,
				  /*hardclip_low*/0,/*hardclip_high*/0,/*lastp*/false);
    FREE(genomicfwd);

#if 0
    /* Intron 2: Not sure how to handle in MD string */
    fprintf(fp,"^");
#endif

    if (substring2 == NULL) {
      print_md_string(fp,matchlength,/*genomicfwd*/0,/*substring2_length*/0,
		      /*hardclip_low*/0,/*hardclip_high*/0,/*lastp*/true);
    } else {
      genomicfwd = (char *) CALLOC(substring2_length+1,sizeof(char));
      genomicdir = Substring_genomic_refdiff(substring2);
      debug2(printf("substring2 genomicdir: %s  start %d length %d\n",genomicdir,substring2_start,substring2_length));
      make_complement_buffered(genomicfwd,&(genomicdir[substring2_start]),substring2_length);
      print_md_string(fp,matchlength,genomicfwd,substring2_length,
		      /*hardclip_low*/0,/*hardclip_high*/0,/*lastp*/true);
      FREE(genomicfwd);
    }
  }

  /* 12. TAGS: SM */
  fprintf(fp,"\t");
  fprintf(fp,"SM:i:%d",mapq_score);

  /* 12. TAGS: XT string (translocation) */
  if (translocationp == true) {
    fprintf(fp,"\t");
    fprintf(fp,"XT:i:1");
  }

  /* 12. TAGS: XS */
  fprintf(fp,"\t");
  if (sensep == Stage3_plusp(shortexon)) {
    fprintf(fp,"XS:A:+");
  } else {
    fprintf(fp,"XS:A:-");
  }

  fprintf(fp,"\n");
  return;
}



/* Distant splicing, including scramble, inversion, translocation */
static void
print_exon_exon (FILE *fp, Stage3_T this, Stage3_T mate, char *acc, int pathnum, int npaths,
		 int mapq_score, IIT_T chromosome_iit, Shortread_T queryseq, int pairedlength,
		 Resulttype_T resulttype, bool translocationp, bool first_read_p,
		 int npaths_mate, Shortread_T queryseq_mate, int quality_shift,
		 char *sam_read_group_id, bool invertp, bool invert_mate_p) {
  if (Stage3_sensedir(this) == SENSE_FORWARD) {
    print_halfdonor(fp,Stage3_substring1(this),this,mate,acc,pathnum,npaths,
		    mapq_score,chromosome_iit,queryseq,pairedlength,
		    resulttype,translocationp,first_read_p,
		    npaths_mate,queryseq_mate,quality_shift,sam_read_group_id,
		    invertp,invert_mate_p,/*use_hardclip_p*/true);
    print_halfacceptor(fp,Stage3_substring2(this),this,mate,acc,pathnum,npaths,
		       mapq_score,chromosome_iit,queryseq,pairedlength,
		       resulttype,translocationp,first_read_p,
		       npaths_mate,queryseq_mate,quality_shift,sam_read_group_id,
		       invertp,invert_mate_p,/*use_hardclip_p*/true);
  } else if (Stage3_sensedir(this) == SENSE_ANTI) {
    print_halfacceptor(fp,Stage3_substring2(this),this,mate,acc,pathnum,npaths,
		       mapq_score,chromosome_iit,queryseq,pairedlength,
		       resulttype,translocationp,first_read_p,
		       npaths_mate,queryseq_mate,quality_shift,sam_read_group_id,
		       invertp,invert_mate_p,/*use_hardclip_p*/true);
    print_halfdonor(fp,Stage3_substring1(this),this,mate,acc,pathnum,npaths,
		    mapq_score,chromosome_iit,queryseq,pairedlength,
		    resulttype,translocationp,first_read_p,
		    npaths_mate,queryseq_mate,quality_shift,sam_read_group_id,
		    invertp,invert_mate_p,/*use_hardclip_p*/true);
  } else {
    abort();
  }

  return;
}



void
SAM_print (FILE *fp, Stage3_T this, Stage3_T mate, char *acc, int pathnum, int npaths,
	   int mapq_score, IIT_T chromosome_iit, Shortread_T queryseq,
	   Shortread_T queryseq_mate, int pairedlength, Resulttype_T resulttype,
	   bool translocationp, bool first_read_p, int npaths_mate, int quality_shift,
	   char *sam_read_group_id, bool invertp, bool invert_mate_p) {
  Hittype_T hittype;
  Substring_T donor, acceptor;
  bool sensep, normalp;

  hittype = Stage3_hittype(this);
  if (hittype == EXACT || hittype == SUB || hittype == TERMINAL) {
    print_single(fp,this,mate,acc,pathnum,npaths,
		 mapq_score,chromosome_iit,queryseq,pairedlength,
		 resulttype,translocationp,first_read_p,
		 npaths_mate,queryseq_mate,quality_shift,sam_read_group_id,
		 invertp,invert_mate_p);

  } else if (hittype == INSERTION) {
    print_insertion(fp,this,mate,acc,pathnum,npaths,
		    mapq_score,chromosome_iit,queryseq,pairedlength,
		    resulttype,translocationp,first_read_p,
		    npaths_mate,queryseq_mate,quality_shift,sam_read_group_id,
		    invertp,invert_mate_p);

  } else if (hittype == DELETION) {
    print_deletion(fp,this,mate,acc,pathnum,npaths,
		   mapq_score,chromosome_iit,queryseq,pairedlength,
		   resulttype,translocationp,first_read_p,
		   npaths_mate,queryseq_mate,quality_shift,sam_read_group_id,
		   invertp,invert_mate_p);

  } else if (hittype == HALFSPLICE_DONOR) {
    print_halfdonor(fp,Stage3_substring1(this),this,mate,acc,pathnum,npaths,
		    mapq_score,chromosome_iit,queryseq,pairedlength,
		    resulttype,translocationp,first_read_p,
		    npaths_mate,queryseq_mate,quality_shift,sam_read_group_id,
		    invertp,invert_mate_p,/*use_hardclip_p*/false);

  } else if (hittype == HALFSPLICE_ACCEPTOR) {
    print_halfacceptor(fp,Stage3_substring1(this),this,mate,acc,pathnum,npaths,
		       mapq_score,chromosome_iit,queryseq,pairedlength,
		       resulttype,translocationp,first_read_p,
		       npaths_mate,queryseq_mate,quality_shift,sam_read_group_id,
		       invertp,invert_mate_p,/*use_hardclip_p*/false);

  } else if (hittype == SPLICE) {
    /* Follows print_splice_distance() in substring.c */
    donor = Stage3_substring1(this);
    acceptor = Stage3_substring2(this);

    if (donor == NULL || acceptor == NULL) {
      abort();
    } else if (Stage3_distance(this) == 0U) {
      /* translocation */
      print_exon_exon(fp,this,mate,acc,pathnum,npaths,
		      mapq_score,chromosome_iit,queryseq,pairedlength,
		      resulttype,translocationp,first_read_p,
		      npaths_mate,queryseq_mate,quality_shift,sam_read_group_id,
		      invertp,invert_mate_p);
    } else {
      normalp = true;
      sensep = (Stage3_sensedir(this) == SENSE_FORWARD);

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
			  mapq_score,chromosome_iit,queryseq,pairedlength,
			  resulttype,translocationp,first_read_p,
			  npaths_mate,queryseq_mate,quality_shift,sam_read_group_id,
			  invertp,invert_mate_p);
      } else {
	print_exon_exon(fp,this,mate,acc,pathnum,npaths,
			mapq_score,chromosome_iit,queryseq,pairedlength,
			resulttype,translocationp,first_read_p,
			npaths_mate,queryseq_mate,quality_shift,sam_read_group_id,
			invertp,invert_mate_p);
      }
    }

  } else if (hittype == ONE_THIRD_SHORTEXON || hittype == TWO_THIRDS_SHORTEXON || hittype == SHORTEXON) {
    print_shortexon(fp,this,mate,acc,pathnum,npaths,
		    mapq_score,chromosome_iit,queryseq,pairedlength,
		    resulttype,translocationp,first_read_p,
		    npaths_mate,queryseq_mate,quality_shift,sam_read_group_id,
		    invertp,invert_mate_p);
  } else {
    abort();
  }

  return;
}



void
SAM_print_paired (Result_T result, Resulttype_T resulttype, bool translocationp,
		  IIT_T chromosome_iit, Shortread_T queryseq1, Shortread_T queryseq2,
		  int maxpaths, bool quiet_if_excessive_p,
		  bool invert_first_p, bool invert_second_p,
		  bool nofailsp, bool failsonlyp, bool fails_as_input_p,
		  bool fastq_format_p, int quality_shift, char *sam_read_group_id,
		  FILE *fp_nomapping_1, FILE *fp_nomapping_2,
		  FILE *fp_unpaired_uniq, FILE *fp_unpaired_mult,
		  FILE *fp_halfmapping_uniq, FILE *fp_halfmapping_mult,
		  FILE *fp_paired_uniq_inv, FILE *fp_paired_uniq_scr, FILE *fp_paired_uniq_long,
		  FILE *fp_paired_mult, FILE *fp_concordant_uniq, FILE *fp_concordant_mult) {
  Stage3pair_T *stage3pairarray, stage3pair;
  Stage3_T *stage3array1, *stage3array2, stage3, mate;
  int npaths, npaths1, npaths2, pathnum;
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
      SAM_print_nomapping(fp_nomapping_1,queryseq1,/*mate*/(Stage3_T) NULL,
			  acc,chromosome_iit,resulttype,translocationp,
			  /*first_read_p*/true,/*npaths_mate*/0,
			  /*queryseq_mate*/queryseq2,quality_shift,sam_read_group_id,
			  invert_first_p,invert_second_p);
      SAM_print_nomapping(fp_nomapping_1,queryseq2,/*mate*/(Stage3_T) NULL,
			  acc,chromosome_iit,resulttype,translocationp,
			  /*first_read_p*/false,/*npaths_mate*/0,
			  /*queryseq_mate*/queryseq1,quality_shift,sam_read_group_id,
			  invert_second_p,invert_first_p);
    }

  } else {
    if (failsonlyp == true) {
      /* Unwanted success: skip */

    } else if (resulttype == CONCORDANT_UNIQ) {
      stage3pairarray = (Stage3pair_T *) Result_array(&npaths,result);
      /* Stage3pair_eval(stage3pairarray,npaths,maxpaths,queryseq1,queryseq2); */

      stage3pair = stage3pairarray[0];

      /* print first end */
      SAM_print(fp_concordant_uniq,Stage3pair_hit5(stage3pair),/*mate*/Stage3pair_hit3(stage3pair),acc,
		/*pathnum*/1,/*npaths*/1,Stage3pair_mapq_score(stage3pair),
		chromosome_iit,/*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
		Stage3pair_pairlength(stage3pair),resulttype,translocationp,
		/*first_read_p*/true,/*npaths_mate*/npaths,quality_shift,sam_read_group_id,
		invert_first_p,invert_second_p);

      /* print second end */
      SAM_print(fp_concordant_uniq,Stage3pair_hit3(stage3pair),/*mate*/Stage3pair_hit5(stage3pair),acc,
		/*pathnum*/1,/*npaths*/1,Stage3pair_mapq_score(stage3pair),
		chromosome_iit,/*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
		Stage3pair_pairlength(stage3pair),resulttype,translocationp,
		/*first_read_p*/false,/*npaths_mate*/npaths,quality_shift,sam_read_group_id,
		invert_second_p,invert_first_p);
    
    } else if (resulttype == CONCORDANT_MULT) {
      stage3pairarray = (Stage3pair_T *) Result_array(&npaths,result);

      if (quiet_if_excessive_p && npaths > maxpaths) {
	/* Print as nomapping, but send to fp_concordant_mult */
	SAM_print_nomapping(fp_concordant_mult,queryseq1,/*mate*/(Stage3_T) NULL,
			    acc,chromosome_iit,resulttype,translocationp,
			    /*first_read_p*/true,/*npaths_mate*/npaths,
			    /*queryseq_mate*/queryseq2,quality_shift,sam_read_group_id,
			    invert_first_p,invert_second_p);
	SAM_print_nomapping(fp_concordant_mult,queryseq2,/*mate*/(Stage3_T) NULL,
			    acc,chromosome_iit,resulttype,translocationp,
			    /*first_read_p*/false,/*npaths_mate*/npaths,
			    /*queryseq_mate*/queryseq1,quality_shift,sam_read_group_id,
			    invert_second_p,invert_first_p);

      } else {
	/* Stage3pair_eval(stage3pairarray,npaths,maxpaths,queryseq1,queryseq2); */
	for (pathnum = 1; pathnum <= npaths && pathnum <= maxpaths; pathnum++) {
	  stage3pair = stage3pairarray[pathnum-1];

	  /* print first end */
	  SAM_print(fp_concordant_mult,Stage3pair_hit5(stage3pair),/*mate*/Stage3pair_hit3(stage3pair),acc,
		    pathnum,npaths,Stage3pair_mapq_score(stage3pair),
		    chromosome_iit,/*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
		    Stage3pair_pairlength(stage3pair),resulttype,translocationp,
		    /*first_read_p*/true,/*npaths_mate*/npaths,quality_shift,sam_read_group_id,
		    invert_first_p,invert_second_p);

	  /* print second end */
	  SAM_print(fp_concordant_mult,Stage3pair_hit3(stage3pair),/*mate*/Stage3pair_hit5(stage3pair),acc,
		    pathnum,npaths,Stage3pair_mapq_score(stage3pair),
		    chromosome_iit,/*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
		    Stage3pair_pairlength(stage3pair),resulttype,translocationp,
		    /*first_read_p*/false,/*npaths_mate*/npaths,quality_shift,sam_read_group_id,
		    invert_second_p,invert_first_p);
	}
      }

    } else if (resulttype == PAIRED_UNIQ) {
      stage3pairarray = (Stage3pair_T *) Result_array(&npaths,result);
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

      /* print first end */
      SAM_print(fp,Stage3pair_hit5(stage3pair),/*mate*/Stage3pair_hit3(stage3pair),acc,
		/*pathnum*/1,/*npaths*/1,Stage3pair_mapq_score(stage3pair),
		chromosome_iit,/*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
		Stage3pair_pairlength(stage3pair),resulttype,translocationp,
		/*first_read_p*/true,/*npaths_mate*/npaths,quality_shift,sam_read_group_id,
		invert_first_p,invert_second_p);

      /* print second end */
      SAM_print(fp,Stage3pair_hit3(stage3pair),/*mate*/Stage3pair_hit5(stage3pair),acc,
		/*pathnum*/1,/*npaths*/1,Stage3pair_mapq_score(stage3pair),
		chromosome_iit,/*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
		Stage3pair_pairlength(stage3pair),resulttype,translocationp,
		/*first_read_p*/false,/*npaths_mate*/npaths,quality_shift,sam_read_group_id,
		invert_second_p,invert_first_p);

    } else if (resulttype == PAIRED_MULT) {
      stage3pairarray = (Stage3pair_T *) Result_array(&npaths,result);

      if (quiet_if_excessive_p && npaths > maxpaths) {
	/* Print as nomapping, but send to fp_concordant_mult */
	SAM_print_nomapping(fp_paired_mult,queryseq1,/*mate*/(Stage3_T) NULL,
			    acc,chromosome_iit,resulttype,translocationp,
			    /*first_read_p*/true,/*npaths_mate*/npaths,
			    /*queryseq_mate*/queryseq2,quality_shift,sam_read_group_id,
			    invert_first_p,invert_second_p);
	SAM_print_nomapping(fp_paired_mult,queryseq2,/*mate*/(Stage3_T) NULL,
			    acc,chromosome_iit,resulttype,translocationp,
			    /*first_read_p*/false,/*npaths_mate*/npaths,
			    /*queryseq_mate*/queryseq1,quality_shift,sam_read_group_id,
			    invert_second_p,invert_first_p);

      } else {
	/* Stage3pair_eval(stage3pairarray,npaths,maxpaths,queryseq1,queryseq2); */
	for (pathnum = 1; pathnum <= npaths && pathnum <= maxpaths; pathnum++) {
	  stage3pair = stage3pairarray[pathnum-1];

	  /* print first end */
	  SAM_print(fp_paired_mult,Stage3pair_hit5(stage3pair),/*mate*/Stage3pair_hit3(stage3pair),acc,
		    pathnum,npaths,Stage3pair_mapq_score(stage3pair),
		    chromosome_iit,/*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
		    Stage3pair_pairlength(stage3pair),resulttype,translocationp,
		    /*first_read_p*/true,/*npaths_mate*/npaths,quality_shift,sam_read_group_id,
		    invert_first_p,invert_second_p);

	  /* print second end */
	  SAM_print(fp_paired_mult,Stage3pair_hit3(stage3pair),/*mate*/Stage3pair_hit5(stage3pair),acc,
		    pathnum,npaths,Stage3pair_mapq_score(stage3pair),
		    chromosome_iit,/*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
		    Stage3pair_pairlength(stage3pair),resulttype,translocationp,
		    /*first_read_p*/false,/*npaths_mate*/npaths,quality_shift,sam_read_group_id,
		    invert_second_p,invert_first_p);
	}
      }

    } else if (resulttype == UNPAIRED_UNIQ) {
      /* Should print mate information in this situation */
      stage3array1 = (Stage3_T *) Result_array(&npaths1,result);
      stage3array2 = (Stage3_T *) Result_array2(&npaths2,result);

      /* print first end */
      /* Stage3_eval_and_sort(stage3array1,npaths1,maxpaths,queryseq1); */
      SAM_print(fp_unpaired_uniq,stage3array1[0],/*mate*/stage3array2[0],acc,
		/*pathnum*/1,/*npaths*/1,Stage3_mapq_score(stage3array1[0]),
		chromosome_iit,/*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
		/*pairedlength*/0U,resulttype,translocationp,
		/*first_read_p*/true,/*npaths_mate*/1,quality_shift,sam_read_group_id,
		invert_first_p,invert_second_p);

      /* print second end */
      /* Stage3_eval_and_sort(stage3array2,npaths2,maxpaths,queryseq2); */
      SAM_print(fp_unpaired_uniq,stage3array2[0],/*mate*/stage3array1[0],acc,
		/*pathnum*/1,/*npaths*/1,Stage3_mapq_score(stage3array2[0]),
		chromosome_iit,/*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
		/*pairedlength*/0U,resulttype,translocationp,
		/*first_read_p*/false,/*npaths_mate*/1,quality_shift,sam_read_group_id,
		invert_second_p,invert_first_p);

    } else if (resulttype == UNPAIRED_MULT) {
      stage3array1 = (Stage3_T *) Result_array(&npaths1,result);
      stage3array2 = (Stage3_T *) Result_array2(&npaths2,result);

#if 0
      /* Do eval and sorting first */
      if (npaths1 == 1) {
	/* Stage3_eval_and_sort(stage3array1,npaths1,maxpaths,queryseq1); */
      } else if (quiet_if_excessive_p && npaths1 > maxpaths) {
	/* Don't sort */
      } else {
	/* Stage3_eval_and_sort(stage3array1,npaths1,maxpaths,queryseq1); */
      }
#endif

#if 0
      if (npaths2 == 1) {
	/* Stage3_eval_and_sort(stage3array2,npaths2,maxpaths,queryseq2); */
      } else if (quiet_if_excessive_p && npaths2 > maxpaths) {
	/* Don't sort */
      } else {
	/* Stage3_eval_and_sort(stage3array2,npaths2,maxpaths,queryseq2); */
      }
#endif

      /* print first end results */
      mate = (npaths2 == 0) ? (Stage3_T) NULL : stage3array2[0];

      if (npaths1 == 1) {
	stage3 = stage3array1[0];
	SAM_print(fp_unpaired_mult,stage3,mate,acc,/*pathnum*/1,npaths1,Stage3_mapq_score(stage3),
		  chromosome_iit,/*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
		  /*pairedlength*/0U,resulttype,translocationp,/*first_read_p*/true,
		  /*npaths_mate*/npaths2,quality_shift,sam_read_group_id,
		  invert_first_p,invert_second_p);

      } else if (quiet_if_excessive_p && npaths1 > maxpaths) {
	SAM_print_nomapping(fp_unpaired_mult,queryseq1,mate,acc,chromosome_iit,
			    resulttype,translocationp,/*first_read_p*/true,
			    /*npaths_mate*/npaths2,/*queryseq_mate*/queryseq2,
			    quality_shift,sam_read_group_id,invert_first_p,invert_second_p);

      } else {
	for (pathnum = 1; pathnum <= npaths1 && pathnum <= maxpaths; pathnum++) {
	  stage3 = stage3array1[pathnum-1];
	  SAM_print(fp_unpaired_mult,stage3,mate,acc,pathnum,npaths1,Stage3_mapq_score(stage3),
		    chromosome_iit,/*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
		    /*pairedlength*/0U,resulttype,translocationp,/*first_read_p*/true,
		    /*npaths_mate*/npaths2,quality_shift,sam_read_group_id,
		    invert_first_p,invert_second_p);
	}
      }
			  
      /* print second end results */
      mate = (npaths1 == 0) ? (Stage3_T) NULL : stage3array1[0];

      if (npaths2 == 1) {
	stage3 = stage3array2[0];
	SAM_print(fp_unpaired_mult,stage3,mate,acc,/*pathnum*/1,npaths2,Stage3_mapq_score(stage3),
		  chromosome_iit,/*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
		  /*pairedlength*/0U,resulttype,translocationp,/*first_read_p*/false,
		  /*npaths_mate*/npaths1,quality_shift,sam_read_group_id,
		  invert_second_p,invert_first_p);

      } else if (quiet_if_excessive_p && npaths2 > maxpaths) {
	SAM_print_nomapping(fp_unpaired_mult,queryseq2,mate,acc,chromosome_iit,
			    resulttype,translocationp,/*first_read_p*/false,
			    /*npaths_mate*/npaths1,/*queryseq_mate*/queryseq1,
			    quality_shift,sam_read_group_id,invert_second_p,invert_first_p);

      } else {
	for (pathnum = 1; pathnum <= npaths2 && pathnum <= maxpaths; pathnum++) {
	  stage3 = stage3array2[pathnum-1];
	  SAM_print(fp_unpaired_mult,stage3,mate,acc,pathnum,npaths2,Stage3_mapq_score(stage3),
		    chromosome_iit,/*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
		    /*pairedlength*/0U,resulttype,translocationp,/*first_read_p*/false,
		    /*npaths_mate*/npaths1,quality_shift,sam_read_group_id,
		    invert_second_p,invert_first_p);
	}
      }

    } else {
      if (resulttype == HALFMAPPING_UNIQ) {
	fp = fp_halfmapping_uniq;
      } else if (resulttype == HALFMAPPING_MULT) {
	fp = fp_halfmapping_mult;
      } else {
	abort();
      }

      stage3array1 = (Stage3_T *) Result_array(&npaths1,result);
      stage3array2 = (Stage3_T *) Result_array2(&npaths2,result);

#if 0
      /* Do eval and sorting first */
      if (npaths1 == 0) {
	/* Nothing to sort */
      } else if (npaths1 == 1) {
	/* Stage3_eval_and_sort(stage3array1,npaths1,maxpaths,queryseq1); */
      } else if (quiet_if_excessive_p && npaths1 > maxpaths) {
	/* Don't sort */
      } else {
	/* Stage3_eval_and_sort(stage3array1,npaths1,maxpaths,queryseq1); */
      }
#endif

#if 0
      if (npaths2 == 0) {
	/* Nothing to sort */
      } else if (npaths2 == 1) {
	/* Stage3_eval_and_sort(stage3array2,npaths2,maxpaths,queryseq2); */
      } else if (quiet_if_excessive_p && npaths2 > maxpaths) {
	/* Don't sort */
      } else {
	/* Stage3_eval_and_sort(stage3array2,npaths2,maxpaths,queryseq2); */
      }
#endif


      /* print first end results */
      mate = (npaths2 == 0) ? (Stage3_T) NULL : stage3array2[0];

      if (npaths1 == 0) {
	/* mate should be non-NULL here */
	SAM_print_nomapping(fp,queryseq1,mate,acc,chromosome_iit,resulttype,translocationp,
			    /*first_read_p*/true,/*npaths_mate*/npaths2,
			    /*queryseq_mate*/queryseq2,quality_shift,sam_read_group_id,
			    invert_first_p,invert_second_p);

      } else if (npaths1 == 1) {
	/* mate should be NULL here */

	stage3 = stage3array1[0];
	SAM_print(fp,stage3,mate,acc,/*pathnum*/1,npaths1,Stage3_mapq_score(stage3),
		  chromosome_iit,/*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
		  /*pairedlength*/0U,resulttype,translocationp,/*first_read_p*/true,
		  /*npaths_mate*/npaths2,quality_shift,sam_read_group_id,
		  invert_first_p,invert_second_p);

      } else if (quiet_if_excessive_p && npaths1 > maxpaths) {
	/* mate should be NULL here */
	SAM_print_nomapping(fp,queryseq1,mate,acc,chromosome_iit,resulttype,translocationp,
			    /*first_read_p*/true,/*npaths_mate*/npaths2,
			    /*queryseq_mate*/queryseq2,quality_shift,sam_read_group_id,
			    invert_first_p,invert_second_p);

      } else {
	/* mate should be NULL here */
	for (pathnum = 1; pathnum <= npaths1 && pathnum <= maxpaths; pathnum++) {
	  stage3 = stage3array1[pathnum-1];
	  SAM_print(fp,stage3,mate,acc,pathnum,npaths1,Stage3_mapq_score(stage3),
		    chromosome_iit,/*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
		    /*pairedlength*/0U,resulttype,translocationp,/*first_read_p*/true,
		    /*npaths_mate*/npaths2,quality_shift,sam_read_group_id,
		    invert_first_p,invert_second_p);
	}
      }
			  
      /* print second end results */
      mate = (npaths1 == 0) ? (Stage3_T) NULL : stage3array1[0];

      if (npaths2 == 0) {
	/* mate should be non-NULL here */
	SAM_print_nomapping(fp,queryseq2,mate,acc,chromosome_iit,resulttype,translocationp,
			    /*first_read_p*/false,/*npaths_mate*/npaths1,
			    /*queryseq_mate*/queryseq1,quality_shift,sam_read_group_id,
			    invert_second_p,invert_first_p);

      } else if (npaths2 == 1) {
	/* mate should be NULL here */

	stage3 = stage3array2[0];
	SAM_print(fp,stage3,mate,acc,/*pathnum*/1,npaths2,Stage3_mapq_score(stage3),
		  chromosome_iit,/*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
		  /*pairedlength*/0U,resulttype,translocationp,/*first_read_p*/false,
		  /*npaths_mate*/npaths1,quality_shift,sam_read_group_id,
		  invert_second_p,invert_first_p);

      } else if (quiet_if_excessive_p && npaths2 > maxpaths) {
	/* mate should be NULL here */
	SAM_print_nomapping(fp,queryseq2,mate,acc,chromosome_iit,resulttype,translocationp,
			    /*first_read_p*/false,/*npaths_mate*/npaths1,
			    /*queryseq_mate*/queryseq1,quality_shift,sam_read_group_id,
			    invert_second_p,invert_first_p);

      } else {
	/* mate should be NULL here */
	for (pathnum = 1; pathnum <= npaths2 && pathnum <= maxpaths; pathnum++) {
	  stage3 = stage3array2[pathnum-1];
	  SAM_print(fp,stage3,mate,acc,pathnum,npaths2,Stage3_mapq_score(stage3),
		    chromosome_iit,/*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
		    /*pairedlength*/0U,resulttype,translocationp,/*first_read_p*/false,
		    /*npaths_mate*/npaths1,quality_shift,sam_read_group_id,
		    invert_second_p,invert_first_p);
	}
      }

    }
  }

  return;
}



