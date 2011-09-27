static char rcsid[] = "$Id: sam.c,v 1.10 2010-07-22 00:11:06 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "sam.h"
#include <stdio.h>
#include <stdlib.h>

#include "mem.h"
#include "complement.h"
#include "stage3hr.h"
#include "stage1hr.h"		/* for MAX_QUERYLENGTH */


#define SANGER_ILLUMINA_DIFF 31


static unsigned int
compute_flag (Substring_T substring, Stage3_T mate, Resulttype_T resulttype,
	      bool first_read_p, int pathnum, int npaths, int npaths_mate) {
  unsigned int flag = 0U;

  if (npaths == 0) {
    flag |= QUERY_UNMAPPED;
  } else {
    if (Substring_plusp(substring) != first_read_p) {
      flag |= QUERY_MINUSP;
    }
  }

  if (resulttype == SINGLEEND_READ) {
  } else {
    flag |= PAIRED_READ;
    if (first_read_p == true) {
      flag |= FIRST_READ_P;
    } else {
      flag |= SECOND_READ_P;
    }
    if (npaths_mate == 0) {
      flag |= MATE_UNMAPPED;

    } else if (mate == NULL) {
      /* Unpaired; no mate.  Not clear if should be MATE_UNMAPPED. */

    } else if (resulttype == PAIREDEND_CONCORDANT || resulttype == PAIREDEND_SAMECHR_SINGLE || 
	       resulttype == PAIREDEND_SAMECHR_MULTIPLE || resulttype == PAIREDEND_AS_SINGLES_UNIQUE) {
      flag |= PAIRED_MAPPING;
      if (Stage3_chrnum(mate) == 0) {
	/* Splice without a direction */
      } else if (Stage3_plusp(mate) == first_read_p) {
	flag |= MATE_MINUSP;
      }

    } else {
      if (Stage3_chrnum(mate) == 0) {
	/* Splice without a direction */
      } else if (Stage3_plusp(mate) == first_read_p) {
	flag |= MATE_MINUSP;
      }
    }
  }

  if (pathnum > 1) {
    flag |= NOT_PRIMARY;
  }

  return flag;
}


static void
print_chromosomal_pos (Stage3_T this, Stage3_T other, IIT_T chromosome_iit, bool matep) {
  Substring_T substring;
  Genomicpos_T chrpos;
  Chrnum_T chrnum;
  char *chr;
  bool allocp;

  if (this == NULL) {
    printf("\t*\t0");
    return;

  } else if (matep == true && other != NULL && Stage3_chrnum(other) > 0 && 
	     Stage3_chrnum(this) == Stage3_chrnum(other)) {
    printf("\t=");

  } else if ((chrnum = Stage3_chrnum(this)) == 0) {
    /* Interchromosomal splice */
    if (matep == true) {
      /* No way to express two chromosomes for mate */
      printf("\t*\t0");
      return;
    } else {
      fprintf(stderr,"Trying to print interchrosomal splice in one line\n");
      abort();
    }

  } else {
    chr = IIT_label(chromosome_iit,chrnum,&allocp);
    printf("\t%s",chr);
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
  printf("\t%u",chrpos+1U);

  return;
}


static void
print_substring_pos (Substring_T substring, IIT_T chromosome_iit) {
  Genomicpos_T chrpos;
  char *chr;
  bool allocp;

  chr = IIT_label(chromosome_iit,Substring_chrnum(substring),&allocp);
  printf("\t%s",chr);
  if (allocp == true) {
    FREE(chr);
  }

  if (Substring_plusp(substring) == true) {
    chrpos = Substring_alignstart_trim(substring) - Substring_chroffset(substring);
    /* Add 1 to report in 1-based coordinates */
    printf("\t%u",chrpos+1U);
  } else {
    chrpos = Substring_alignend_trim(substring) - Substring_chroffset(substring);
    /* Add 1 to report in 1-based coordinates */
    printf("\t%u",chrpos+1U);
  }

  return;
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
SAM_print_nomapping (Sequence_T queryseq, Stage3_T mate, char *acc,
		     IIT_T chromosome_iit, Resulttype_T resulttype, bool first_read_p,
		     int npaths_mate, Sequence_T queryseq_mate, int quality_shift) {
  unsigned int flag;

  /* 1. QNAME */
  printf("%s",acc);
  
  /* 2. FLAG */
  flag = compute_flag(/*substring*/NULL,/*mate*/NULL,resulttype,first_read_p,
		      /*pathnum*/0,/*npaths*/0,npaths_mate);
  printf("\t%u",flag);

  /* 3. RNAME: chr */
  printf("\t*");

  /* 4. POS: chrpos */
  printf("\t0");

  /* 5. MAPQ: Mapping quality */
  /* Picard says MAPQ should be 0 for an unmapped read */
  printf("\t0");

  /* 6. CIGAR */
  printf("\t*");

  /* 7. MRNM: Mate chr */
  /* 8. MPOS: Mate chrpos */
  print_chromosomal_pos(mate,/*other*/NULL,chromosome_iit,/*matep*/true);

  /* 9. ISIZE: Insert size */
  printf("\t0");

  /* 10. SEQ: queryseq and 11. QUAL: quality scores */
  /* Since there is no mapping, we print the original query sequence. */
  printf("\t");
  if (first_read_p == true) {
    Sequence_print_chopped(stdout,queryseq,/*hardclip_low*/0,/*hardclip_high*/0);
    printf("\t");
    Sequence_print_quality(stdout,queryseq,/*hardclip_low*/0,/*hardclip_high*/0,
			   quality_shift);
  } else {
    Sequence_print_chopped_revcomp(stdout,queryseq,/*hardclip_low*/0,/*hardclip_high*/0);
    printf("\t");
    Sequence_print_quality_revcomp(stdout,queryseq,/*hardclip_low*/0,/*hardclip_high*/0,
				   quality_shift);
  }

  /* 12. TAGS */
#if 0
  /* Do not terminate line with tab */
  printf("\t");
#endif
  
  printf("\n");
  return;
}


static void
print_single (Stage3_T this, Stage3_T mate, char *acc, int pathnum, int npaths,
	      int score, IIT_T chromosome_iit, Sequence_T queryseq, int pairedlength,
	      IIT_T snps_iit, int *snps_divint_crosstable, Resulttype_T resulttype, bool first_read_p,
	      int npaths_mate, Sequence_T queryseq_mate, int quality_shift) {
  unsigned int flag = 0U;
  Substring_T substring;
  int querylength;

  querylength = Sequence_fulllength(queryseq);
  substring = Stage3_substring1(this);

  /* 1. QNAME */
  printf("%s",acc);

  /* 2. FLAG */
  flag = compute_flag(substring,mate,resulttype,first_read_p,
		      pathnum,npaths,npaths_mate);
  printf("\t%u",flag);

  /* 3. RNAME: chr */
  /* 4. POS: chrpos */
  print_chromosomal_pos(this,mate,chromosome_iit,/*matep*/false);

  /* 5. MAPQ: Mapping quality */
  printf("\t255");

  /* 6. CIGAR */
  printf("\t");
  if (Stage3_plusp(this) == true) {
    if (Substring_querystart(substring) > 0) {
      printf("%dS",Substring_querystart(substring));
    }
  } else {
    if (Substring_queryend(substring) < querylength) {
      printf("%dS",querylength - Substring_queryend(substring));
    }
  }
  printf("%dM",Substring_match_length(substring));
  if (Stage3_plusp(this) == true) {
    if (Substring_queryend(substring) < querylength) {
      printf("%dS",querylength - Substring_queryend(substring));
    }
  } else {
    if (Substring_querystart(substring) > 0) {
      printf("%dS",Substring_querystart(substring));
    }
  }


  /* 7. MRNM: Mate chr */
  /* 8. MPOS: Mate chrpos */
  if (resulttype == PAIREDEND_CONCORDANT || 
      resulttype == PAIREDEND_SAMECHR_SINGLE ||
      resulttype == PAIREDEND_SAMECHR_MULTIPLE ||
      resulttype == PAIREDEND_AS_SINGLES_UNIQUE) {
    print_chromosomal_pos(mate,this,chromosome_iit,/*matep*/true);
  } else {
    printf("\t*\t0");
  }


  /* 9. ISIZE: Insert size */
  if (resulttype != PAIREDEND_CONCORDANT && resulttype != PAIREDEND_SAMECHR_SINGLE && 
      resulttype != PAIREDEND_SAMECHR_MULTIPLE && resulttype != PAIREDEND_AS_SINGLES_UNIQUE) {
    printf("\t0");
  } else if (Stage3_plusp(this) == first_read_p) {
    printf("\t%d",pairedlength);
  } else {
    printf("\t%d",-pairedlength);
  }
  
  /* 10. SEQ: queryseq and 11. QUAL: quality scores */
  /* Queryseq has already been inverted, so just measure plusp relative to its current state */
  printf("\t");
  if (Stage3_plusp(this) == true) {
    Sequence_print_chopped(stdout,queryseq,/*hardclip_low*/0,/*hardclip_high*/0);
    printf("\t");
    Sequence_print_quality(stdout,queryseq,/*hardclip_low*/0,/*hardclip_high*/0,
			   quality_shift);
  } else {
    Sequence_print_chopped_revcomp(stdout,queryseq,/*hardclip_low*/0,/*hardclip_high*/0);
    printf("\t");
    Sequence_print_quality_revcomp(stdout,queryseq,/*hardclip_low*/0,/*hardclip_high*/0,
				   quality_shift);
  }

  /* 12. TAGS */
  printf("\t");
  printf("NM:i:%d",Stage3_nmismatches(this));

  printf("\n");
  return;
}

static void
print_insertion (Stage3_T this, Stage3_T mate, char *acc, int pathnum, int npaths,
		 int score, IIT_T chromosome_iit, Sequence_T queryseq, int pairedlength,
		 IIT_T snps_iit, int *snps_divint_crosstable, Resulttype_T resulttype, bool first_read_p,
		 int npaths_mate, Sequence_T queryseq_mate, int quality_shift) {
  unsigned int flag = 0U;
  Substring_T substring1, substring2;
  int querylength;

  querylength = Sequence_fulllength(queryseq);

  substring1 = Stage3_substring1(this);
  substring2 = Stage3_substring2(this);

  /* 1. QNAME */
  printf("%s",acc);

  /* 2. FLAG */
  flag = compute_flag(substring1,mate,resulttype,first_read_p,
		      pathnum,npaths,npaths_mate);
  printf("\t%u",flag);

  /* 3. RNAME: chr */
  /* 4. POS: chrpos */
  print_chromosomal_pos(this,mate,chromosome_iit,/*matep*/false);

  /* 5. MAPQ: Mapping quality */
  printf("\t255");

  /* 6. CIGAR */
  printf("\t");

  if (Stage3_plusp(this) == true) {
    if (Substring_querystart(substring1) > 0) {
      printf("%dS",Substring_querystart(substring1));
    }
    printf("%dM",Substring_match_length(substring1));
    printf("%dI",Stage3_nindels(this));
    printf("%dM",Substring_match_length(substring2));
    if (Substring_queryend(substring2) < querylength) {
      printf("%dS",querylength - Substring_queryend(substring2));
    }

  } else {
    if (Substring_queryend(substring2) < querylength) {
      printf("%dS",querylength - Substring_queryend(substring2));
    }
    printf("%dM",Substring_match_length(substring2));
    printf("%dI",Stage3_nindels(this));
    printf("%dM",Substring_match_length(substring1));
    if (Substring_querystart(substring1) > 0) {
      printf("%dS",Substring_querystart(substring1));
    }
  }

  /* 7. MRNM: Mate chr */
  /* 8. MPOS: Mate chrpos */
  if (resulttype == PAIREDEND_CONCORDANT || resulttype == PAIREDEND_SAMECHR_SINGLE || 
      resulttype == PAIREDEND_SAMECHR_MULTIPLE || resulttype == PAIREDEND_AS_SINGLES_UNIQUE) {
    print_chromosomal_pos(mate,this,chromosome_iit,/*matep*/true);
  } else {
    printf("\t*\t0");
  }

  /* 9. ISIZE: Insert size */
  if (resulttype != PAIREDEND_CONCORDANT && resulttype != PAIREDEND_SAMECHR_SINGLE && 
      resulttype != PAIREDEND_SAMECHR_MULTIPLE && resulttype != PAIREDEND_AS_SINGLES_UNIQUE) {
    printf("\t0");
  } else if (Stage3_plusp(this) == first_read_p) {
    printf("\t%d",pairedlength);
  } else {
    printf("\t%d",-pairedlength);
  }
  
  /* 10. SEQ: queryseq and 11. QUAL: quality scores */
  /* Queryseq has already been inverted, so just measure plusp relative to its current state */
  printf("\t");
  if (Stage3_plusp(this) == true) {
    Sequence_print_chopped(stdout,queryseq,/*hardclip_low*/0,/*hardclip_high*/0);
    printf("\t");
    Sequence_print_quality(stdout,queryseq,/*hardclip_low*/0,/*hardclip_high*/0,
			   quality_shift);
  } else {
    Sequence_print_chopped_revcomp(stdout,queryseq,/*hardclip_low*/0,/*hardclip_high*/0);
    printf("\t");
    Sequence_print_quality_revcomp(stdout,queryseq,/*hardclip_low*/0,/*hardclip_high*/0,
				   quality_shift);
  }

  /* 12. TAGS */
  printf("\t");
  printf("NM:i:%d",Stage3_nmismatches(this));
  
  printf("\n");
  return;
}


static void
print_deletion (Stage3_T this, Stage3_T mate, char *acc, int pathnum, int npaths,
		int score, IIT_T chromosome_iit, Sequence_T queryseq, int pairedlength,
		IIT_T snps_iit, int *snps_divint_crosstable, Resulttype_T resulttype, bool first_read_p,
		int npaths_mate, Sequence_T queryseq_mate, int quality_shift) {
  unsigned int flag = 0U;
  Substring_T substring1, substring2;
  int querylength;

  querylength = Sequence_fulllength(queryseq);

  substring1 = Stage3_substring1(this);
  substring2 = Stage3_substring2(this);

  /* 1. QNAME */
  printf("%s",acc);

  /* 2. FLAG */
  flag = compute_flag(substring1,mate,resulttype,first_read_p,
		      pathnum,npaths,npaths_mate);
  printf("\t%u",flag);

  /* 3. RNAME: chr */
  /* 4. POS: chrpos */
  print_chromosomal_pos(this,mate,chromosome_iit,/*matep*/false);

  /* 5. MAPQ: Mapping quality */
  printf("\t255");

  /* 6. CIGAR */
  printf("\t");
  if (Stage3_plusp(this) == true) {
    if (Substring_querystart(substring1) > 0) {
      printf("%dS",Substring_querystart(substring1));
    }
    printf("%dM",Substring_match_length(substring1));
    printf("%dD",Stage3_nindels(this)); /* nindels is positive */
    printf("%dM",Substring_match_length(substring2));
    if (Substring_queryend(substring2) < querylength) {
      printf("%dS",querylength - Substring_queryend(substring2));
    }

  } else {
    if (Substring_queryend(substring2) < querylength) {
      printf("%dS",querylength - Substring_queryend(substring2));
    }
    printf("%dM",Substring_match_length(substring2));
    printf("%dD",Stage3_nindels(this)); /* nindels is positive */
    printf("%dM",Substring_match_length(substring1));
    if (Substring_querystart(substring1) > 0) {
      printf("%dS",Substring_querystart(substring1));
    }
  }

  /* 7. MRNM: Mate chr */
  /* 8. MPOS: Mate chrpos */
  if (resulttype == PAIREDEND_CONCORDANT || resulttype == PAIREDEND_SAMECHR_SINGLE || 
      resulttype == PAIREDEND_SAMECHR_MULTIPLE || resulttype == PAIREDEND_AS_SINGLES_UNIQUE) {
    print_chromosomal_pos(mate,this,chromosome_iit,/*matep*/true);
  } else {
    printf("\t*\t0");
  }

  /* 9. ISIZE: Insert size */
  if (resulttype != PAIREDEND_CONCORDANT && resulttype != PAIREDEND_SAMECHR_SINGLE && 
      resulttype != PAIREDEND_SAMECHR_MULTIPLE && resulttype != PAIREDEND_AS_SINGLES_UNIQUE) {
    printf("\t0");
  } else if (Stage3_plusp(this) == first_read_p) {
    printf("\t%d",pairedlength);
  } else {
    printf("\t%d",-pairedlength);
  }
  
  /* 10. SEQ: queryseq and 11. QUAL: quality scores */
  /* Queryseq has already been inverted, so just measure plusp relative to its current state */
  printf("\t");
  if (Stage3_plusp(this) == true) {
    Sequence_print_chopped(stdout,queryseq,/*hardclip_low*/0,/*hardclip_high*/0);
    printf("\t");
    Sequence_print_quality(stdout,queryseq,/*hardclip_low*/0,/*hardclip_high*/0,
			   quality_shift);
  } else {
    Sequence_print_chopped_revcomp(stdout,queryseq,/*hardclip_low*/0,/*hardclip_high*/0);
    printf("\t");
    Sequence_print_quality_revcomp(stdout,queryseq,/*hardclip_low*/0,/*hardclip_high*/0,
				   quality_shift);
  }

  /* 12. TAGS */
  printf("\t");
  printf("NM:i:%d",Stage3_nmismatches(this));
  
  printf("\n");
  return;
}


#define HARDCLIP_ON_HALFSPLICE


static void
print_halfdonor (Substring_T donor, Stage3_T this, Stage3_T mate, char *acc, int pathnum, int npaths,
		 int score, IIT_T chromosome_iit, Sequence_T queryseq, int pairedlength,
		 IIT_T snps_iit, int *snps_divint_crosstable, Resulttype_T resulttype, bool first_read_p,
		 int npaths_mate, Sequence_T queryseq_mate, int quality_shift) {
  unsigned int flag = 0U;
  int querylength;
  int hardclip_low = 0, hardclip_high = 0;
  bool sensep;

  querylength = Sequence_fulllength(queryseq);

  /* 1. QNAME */
  printf("%s",acc);

  /* 2. FLAG */
  flag = compute_flag(donor,mate,resulttype,first_read_p,
		      pathnum,npaths,npaths_mate);
  printf("\t%u",flag);

  /* 3. RNAME: chr */
  /* 4. POS: chrpos */
#if 0
  print_chromosomal_pos(this,mate,chromosome_iit,/*matep*/false);
#else
  /* Need to do this, because print_chromosomal_pos cannot handle translocations */
  print_substring_pos(donor,chromosome_iit);
#endif

  /* 5. MAPQ: Mapping quality */
  printf("\t255");

  /* 6. CIGAR */
  printf("\t");
  sensep = Substring_chimera_sensep(donor);

  if (sensep == Substring_plusp(donor)) {
    if (Substring_plusp(donor) == true) {
      if (Substring_querystart(donor) > 0) {
	printf("%dS",Substring_querystart(donor));
      }
    } else {
      if (Substring_queryend(donor) < querylength) {
	printf("%dS",querylength - Substring_queryend(donor));
      }
    }
    printf("%dM",Substring_match_length(donor));
    if (sensep == true) {
#ifdef HARDCLIP_ON_HALFSPLICE
      hardclip_high = querylength - Substring_chimera_pos(donor);
      printf("%dH",hardclip_high); /* intron */
#else
      printf("%dS",querylength - Substring_chimera_pos(donor));	      /* intron */
#endif
    } else {
#ifdef HARDCLIP_ON_HALFSPLICE
      hardclip_high = Substring_chimera_pos(donor);
      printf("%dH",hardclip_high); /* intron */
#else
      printf("%dS",Substring_chimera_pos(donor));  /* intron */
#endif
    }

  } else if (sensep != Substring_plusp(donor)) {
    if (sensep == true) {
#ifdef HARDCLIP_ON_HALFSPLICE
      hardclip_low = querylength - Substring_chimera_pos(donor);
      printf("%dH",hardclip_low); /* intron */
#else
      printf("%dS",querylength - Substring_chimera_pos(donor));	      /* intron */
#endif
    } else {
#ifdef HARDCLIP_ON_HALFSPLICE
      hardclip_low = Substring_chimera_pos(donor);
      printf("%dH",hardclip_low); /* intron */
#else
      printf("%dS",Substring_chimera_pos(donor));  /* intron */
#endif
    }

    printf("%dM",Substring_match_length(donor));
    if (Substring_plusp(donor) == true) {
      if (Substring_queryend(donor) < querylength) {
	printf("%dS",querylength - Substring_queryend(donor));
      }
    } else {
      if (Substring_querystart(donor) > 0) {
	printf("%dS",Substring_querystart(donor));
      }
    }
  }

  /* 7. MRNM: Mate chr */
  /* 8. MPOS: Mate chrpos */
  if (resulttype == PAIREDEND_CONCORDANT || resulttype == PAIREDEND_SAMECHR_SINGLE || 
      resulttype == PAIREDEND_SAMECHR_MULTIPLE || resulttype == PAIREDEND_AS_SINGLES_UNIQUE) {
    print_chromosomal_pos(mate,this,chromosome_iit,/*matep*/true);
  } else {
    printf("\t*\t0");
  }

  /* 9. ISIZE: Insert size */
  if (resulttype != PAIREDEND_CONCORDANT && resulttype != PAIREDEND_SAMECHR_SINGLE && 
      resulttype != PAIREDEND_SAMECHR_MULTIPLE && resulttype != PAIREDEND_AS_SINGLES_UNIQUE) {
    printf("\t0");
  } else if (Substring_plusp(donor) == first_read_p) {
    printf("\t%d",pairedlength);
  } else {
    printf("\t%d",-pairedlength);
  }
  
  /* 10. SEQ: queryseq and 11. QUAL: quality scores */
  /* Queryseq has already been inverted, so just measure plusp relative to its current state */
  printf("\t");
  if (Substring_plusp(donor) == true) {
    Sequence_print_chopped(stdout,queryseq,hardclip_low,hardclip_high);
    printf("\t");
    Sequence_print_quality(stdout,queryseq,hardclip_low,hardclip_high,
			   quality_shift);
  } else {
    Sequence_print_chopped_revcomp(stdout,queryseq,hardclip_low,hardclip_high);
    printf("\t");
    Sequence_print_quality_revcomp(stdout,queryseq,hardclip_low,hardclip_high,
				   quality_shift);
  }

  /* 12. TAGS */
  printf("\t");
  printf("NM:i:%d",Substring_nmismatches(donor));
  
  printf("\n");
  return;
}


static void
print_halfacceptor (Substring_T acceptor, Stage3_T this, Stage3_T mate, char *acc, int pathnum, int npaths,
		    int score, IIT_T chromosome_iit, Sequence_T queryseq, int pairedlength,
		    IIT_T snps_iit, int *snps_divint_crosstable, Resulttype_T resulttype, bool first_read_p,
		    int npaths_mate, Sequence_T queryseq_mate, int quality_shift) {
  unsigned int flag = 0U;
  int querylength;
  int hardclip_low = 0, hardclip_high = 0;
  bool sensep;

  querylength = Sequence_fulllength(queryseq);

  /* 1. QNAME */
  printf("%s",acc);

  /* 2. FLAG */
  flag = compute_flag(acceptor,mate,resulttype,first_read_p,
		      pathnum,npaths,npaths_mate);
  printf("\t%u",flag);

  /* 3. RNAME: chr */
  /* 4. POS: chrpos */
#if 0
  print_chromosomal_pos(this,mate,chromosome_iit,/*matep*/false);
#else
  /* Need to do this, because print_chromosomal_pos cannot handle translocations */
  print_substring_pos(acceptor,chromosome_iit);
#endif

  /* 5. MAPQ: Mapping quality */
  printf("\t255");

  /* 6. CIGAR */
  printf("\t");
  sensep = Substring_chimera_sensep(acceptor);

  if (sensep != Substring_plusp(acceptor)) {
    if (Substring_plusp(acceptor) == true) {
      if (Substring_querystart(acceptor) > 0) {
	printf("%dS",Substring_querystart(acceptor));
      }
    } else {
      if (Substring_queryend(acceptor) < querylength) {
	printf("%dS",querylength - Substring_queryend(acceptor));
      }
    }
    printf("%dM",Substring_match_length(acceptor));

    if (sensep == true) {
#ifdef HARDCLIP_ON_HALFSPLICE
      hardclip_high = Substring_chimera_pos(acceptor);
      printf("%dH",hardclip_high);  /* intron */
#else
      printf("%dS",Substring_chimera_pos(acceptor));  /* intron */
#endif
    } else {
#ifdef HARDCLIP_ON_HALFSPLICE
      hardclip_high = querylength - Substring_chimera_pos(acceptor);
      printf("%dH",hardclip_high);  /* intron */
#else
      printf("%dS",querylength - Substring_chimera_pos(acceptor));  /* intron */
#endif
    }

  } else if (sensep == Substring_plusp(acceptor)) {
    if (sensep == true) {
#ifdef HARDCLIP_ON_HALFSPLICE
      hardclip_low = Substring_chimera_pos(acceptor);
      printf("%dH",hardclip_low);  /* intron */
#else
      printf("%dS",Substring_chimera_pos(acceptor));  /* intron */
#endif
    } else {
#ifdef HARDCLIP_ON_HALFSPLICE
      hardclip_low = querylength - Substring_chimera_pos(acceptor);
      printf("%dH",hardclip_low);  /* intron */
#else
      printf("%dS",querylength - Substring_chimera_pos(acceptor));  /* intron */
#endif
    }

    printf("%dM",Substring_match_length(acceptor));
    if (Substring_plusp(acceptor) == true) {
      if (Substring_queryend(acceptor) < querylength) {
	printf("%dS",querylength - Substring_queryend(acceptor));
      }
    } else {
      if (Substring_querystart(acceptor) > 0) {
	printf("%dS",Substring_querystart(acceptor));
      }
    }
  }


  /* 7. MRNM: Mate chr */
  /* 8. MPOS: Mate chrpos */
  if (resulttype == PAIREDEND_CONCORDANT || resulttype == PAIREDEND_SAMECHR_SINGLE || 
      resulttype == PAIREDEND_SAMECHR_MULTIPLE || resulttype == PAIREDEND_AS_SINGLES_UNIQUE) {
    print_chromosomal_pos(mate,this,chromosome_iit,/*matep*/true);
  } else {
    printf("\t*\t0");
  }

  /* 9. ISIZE: Insert size */
  if (resulttype != PAIREDEND_CONCORDANT && resulttype != PAIREDEND_SAMECHR_SINGLE && 
      resulttype != PAIREDEND_SAMECHR_MULTIPLE && resulttype != PAIREDEND_AS_SINGLES_UNIQUE) {
    printf("\t0");
  } else if (Substring_plusp(acceptor) == first_read_p) {
    printf("\t%d",pairedlength);
  } else {
    printf("\t%d",-pairedlength);
  }
  
  /* 10. SEQ: queryseq and 11. QUAL: quality scores */
  /* Queryseq has already been inverted, so just measure plusp relative to its current state */
  printf("\t");
  if (Substring_plusp(acceptor) == true) {
    Sequence_print_chopped(stdout,queryseq,hardclip_low,hardclip_high);
    printf("\t");
    Sequence_print_quality(stdout,queryseq,hardclip_low,hardclip_high,
			   quality_shift);
  } else {
    Sequence_print_chopped_revcomp(stdout,queryseq,hardclip_low,hardclip_high);
    printf("\t");
    Sequence_print_quality_revcomp(stdout,queryseq,hardclip_low,hardclip_high,
				   quality_shift);
  }

  /* 12. TAGS */
  printf("\t");
  printf("NM:i:%d",Substring_nmismatches(acceptor));
  
  printf("\n");
  return;
}



static void
print_localsplice (Stage3_T chimera, Stage3_T mate, char *acc, int pathnum, int npaths,
		   int score, IIT_T chromosome_iit, Sequence_T queryseq, int pairedlength,
		   IIT_T snps_iit, int *snps_divint_crosstable, Resulttype_T resulttype, bool first_read_p,
		   int npaths_mate, Sequence_T queryseq_mate, int quality_shift) {
  unsigned int flag = 0U;
  Substring_T substring1, substring2;
  int querylength;
  bool sensep;

  querylength = Sequence_fulllength(queryseq);

  sensep = (Stage3_sensedir(chimera) == SENSE_FORWARD);


  /* 1. QNAME */
  printf("%s",acc);

  /* 2. FLAG */
  flag = compute_flag(Stage3_substring1(chimera),mate,resulttype,first_read_p,
		      pathnum,npaths,npaths_mate);
  printf("\t%u",flag);

  /* 3. RNAME: chr */
  /* 4. POS: chrpos */
  print_chromosomal_pos(chimera,mate,chromosome_iit,/*matep*/false);

  /* 5. MAPQ: Mapping quality */
  printf("\t255");

  /* 6. CIGAR */
  printf("\t");
  if (sensep == Stage3_plusp(chimera)) {
    substring1 = /* donor */ Stage3_substring1(chimera);
    substring2 = /* acceptor */ Stage3_substring2(chimera);
  } else {
    substring1 = /* acceptor */ Stage3_substring2(chimera);
    substring2 = /* donor */ Stage3_substring1(chimera);
  }
  if (Stage3_plusp(chimera) == true) {
    if (Substring_querystart(substring1) > 0) {
      printf("%dS",Substring_querystart(substring1));
    }
  } else {
    if (Substring_queryend(substring1) < querylength) {
      printf("%dS",querylength - Substring_queryend(substring1));
    }
  }
    
  printf("%dM",Substring_match_length(substring1));
  printf("%uN",Stage3_distance(chimera));
  printf("%dM",Substring_match_length(substring2));

  if (Stage3_plusp(chimera) == true) {
    if (Substring_queryend(substring2) < querylength) {
      printf("%dS",querylength - Substring_queryend(substring2));
    }
  } else {
    if (Substring_querystart(substring2) > 0) {
      printf("%dS",Substring_querystart(substring2));
    }
  }


  /* 7. MRNM: Mate chr */
  /* 8. MPOS: Mate chrpos */
  if (resulttype == PAIREDEND_CONCORDANT || resulttype == PAIREDEND_SAMECHR_SINGLE || 
      resulttype == PAIREDEND_SAMECHR_MULTIPLE || resulttype == PAIREDEND_AS_SINGLES_UNIQUE) {
    print_chromosomal_pos(mate,chimera,chromosome_iit,/*matep*/true);
  } else {
    printf("\t*\t0");
  }

  /* 9. ISIZE: Insert size */
  if (resulttype != PAIREDEND_CONCORDANT && resulttype != PAIREDEND_SAMECHR_SINGLE && 
      resulttype != PAIREDEND_SAMECHR_MULTIPLE && resulttype != PAIREDEND_AS_SINGLES_UNIQUE) {
    printf("\t0");
  } else if (Stage3_plusp(chimera) == first_read_p) {
    printf("\t%d",pairedlength);
  } else {
    printf("\t%d",-pairedlength);
  }
  
  /* 10. SEQ: queryseq and 11. QUAL: quality scores */
  /* Queryseq has already been inverted, so just measure plusp relative to its current state */
  printf("\t");
  if (Stage3_plusp(chimera) == true) {
    Sequence_print_chopped(stdout,queryseq,/*hardclip_low*/0,/*hardclip_high*/0);
    printf("\t");
    Sequence_print_quality(stdout,queryseq,/*hardclip_low*/0,/*hardclip_high*/0,
			   quality_shift);
  } else {
    Sequence_print_chopped_revcomp(stdout,queryseq,/*hardclip_low*/0,/*hardclip_high*/0);
    printf("\t");
    Sequence_print_quality_revcomp(stdout,queryseq,/*hardclip_low*/0,/*hardclip_high*/0,
				   quality_shift);
  } 

  /* 12. TAGS */
  printf("\t");
  printf("NM:i:%d",Stage3_nmismatches(chimera));
  
  printf("\n");
  return;
}


static void
print_terminal (Stage3_T this, Stage3_T mate, char *acc, int pathnum, int npaths,
		int score, IIT_T chromosome_iit, Sequence_T queryseq, int pairedlength,
		IIT_T snps_iit, int *snps_divint_crosstable, Resulttype_T resulttype, bool first_read_p,
		int npaths_mate, Sequence_T queryseq_mate, int quality_shift) {
  Substring_T substring;
  unsigned int flag = 0U;
  int querylength;

  querylength = Sequence_fulllength(queryseq);

  substring = Stage3_substring1(this);

  /* 1. QNAME */
  printf("%s",acc);

  /* 2. FLAG */
  flag = compute_flag(substring,mate,resulttype,first_read_p,
		      pathnum,npaths,npaths_mate);
  printf("\t%u",flag);

  /* 3. RNAME: chr */
  /* 4. POS: chrpos */
  print_chromosomal_pos(this,mate,chromosome_iit,/*matep*/false);

  /* 5. MAPQ: Mapping quality */
  printf("\t255");

  /* 6. CIGAR */
  printf("\t");
  if (Substring_plusp(substring) == true) {
    if (Substring_querystart(substring) > 0) {
      printf("%dS",Substring_querystart(substring));
    }
    printf("%dM",Substring_match_length(substring));
    if (Substring_queryend(substring) < querylength) {
      printf("%dS",querylength - Substring_queryend(substring));
    }

  } else {
    if (Substring_queryend(substring) < querylength) {
      printf("%dS",querylength - Substring_queryend(substring));
    }
    printf("%dM",Substring_match_length(substring));
    if (Substring_querystart(substring) > 0) {
      printf("%dS",Substring_querystart(substring));
    }
  }

  /* 7. MRNM: Mate chr */
  /* 8. MPOS: Mate chrpos */
  if (resulttype == PAIREDEND_CONCORDANT || resulttype == PAIREDEND_SAMECHR_SINGLE || 
      resulttype == PAIREDEND_SAMECHR_MULTIPLE || resulttype == PAIREDEND_AS_SINGLES_UNIQUE) {
    print_chromosomal_pos(mate,this,chromosome_iit,/*matep*/true);
  } else {
    printf("\t*\t0");
  }

  /* 9. ISIZE: Insert size */
  if (resulttype != PAIREDEND_CONCORDANT && resulttype != PAIREDEND_SAMECHR_SINGLE && 
      resulttype != PAIREDEND_SAMECHR_MULTIPLE && resulttype != PAIREDEND_AS_SINGLES_UNIQUE) {
    printf("\t0");
  } else if (Substring_plusp(substring) == first_read_p) {
    printf("\t%d",pairedlength);
  } else {
    printf("\t%d",-pairedlength);
  }
  
  /* 10. SEQ: queryseq and 11. QUAL: quality scores */
  /* Queryseq has already been inverted, so just measure plusp relative to its current state */
  printf("\t");
  if (Substring_plusp(substring) == true) {
    Sequence_print_chopped(stdout,queryseq,/*hardclip_low*/0,/*hardclip_high*/0);
    printf("\t");
    Sequence_print_quality(stdout,queryseq,/*hardclip_low*/0,/*hardclip_high*/0,
			   quality_shift);
  } else {
    Sequence_print_chopped_revcomp(stdout,queryseq,/*hardclip_low*/0,/*hardclip_high*/0);
    printf("\t");
    Sequence_print_quality_revcomp(stdout,queryseq,/*hardclip_low*/0,/*hardclip_high*/0,
				   quality_shift);
  }

  /* 12. TAGS */
  printf("\t");
  printf("NM:i:%d",Substring_nmismatches(substring));
  
  printf("\n");
  return;
}


/* Distant splicing, including scramble, inversion, translocation */
static void
print_exon_exon (Stage3_T this, Stage3_T mate, char *acc, int pathnum, int npaths,
		 int score, IIT_T chromosome_iit, Sequence_T queryseq, int pairedlength,
		 IIT_T snps_iit, int *snps_divint_crosstable, Resulttype_T resulttype, bool first_read_p,
		 int npaths_mate, Sequence_T queryseq_mate, int quality_shift) {
  if (Stage3_sensedir(this) == SENSE_FORWARD) {
    print_halfdonor(Stage3_substring1(this),this,mate,acc,pathnum,npaths,
		    score,chromosome_iit,queryseq,pairedlength,
		    snps_iit,snps_divint_crosstable,resulttype,first_read_p,
		    npaths_mate,queryseq_mate,quality_shift);
    print_halfacceptor(Stage3_substring2(this),this,mate,acc,pathnum,npaths,
		       score,chromosome_iit,queryseq,pairedlength,
		       snps_iit,snps_divint_crosstable,resulttype,first_read_p,
		       npaths_mate,queryseq_mate,quality_shift);
  } else if (Stage3_sensedir(this) == SENSE_ANTI) {
    print_halfacceptor(Stage3_substring2(this),this,mate,acc,pathnum,npaths,
		       score,chromosome_iit,queryseq,pairedlength,
		       snps_iit,snps_divint_crosstable,resulttype,first_read_p,
		       npaths_mate,queryseq_mate,quality_shift);
    print_halfdonor(Stage3_substring1(this),this,mate,acc,pathnum,npaths,
		    score,chromosome_iit,queryseq,pairedlength,
		    snps_iit,snps_divint_crosstable,resulttype,first_read_p,
		    npaths_mate,queryseq_mate,quality_shift);
  } else {
    abort();
  }

  return;
}


void
SAM_print (Stage3_T this, Stage3_T mate, char *acc, int pathnum, int npaths,
	   int score, IIT_T chromosome_iit, Sequence_T queryseq,
	   Sequence_T queryseq_mate, IIT_T snps_iit, int *snps_divint_crosstable,
	   IIT_T splicesites_iit, int donor_typeint, int acceptor_typeint,
	   int *splicesites_divint_crosstable, int pairedlength, Resulttype_T resulttype,
	   bool first_read_p, int npaths_mate, int quality_shift) {
  Hittype_T hittype;
  Substring_T donor, acceptor;
  bool sensep, normalp;

  hittype = Stage3_hittype(this);
  if (hittype == EXACT || hittype == SUB) {
    print_single(this,mate,acc,pathnum,npaths,
		 score,chromosome_iit,queryseq,pairedlength,
		 snps_iit,snps_divint_crosstable,resulttype,first_read_p,
		 npaths_mate,queryseq_mate,quality_shift);
  } else if (hittype == INS) {
    print_insertion(this,mate,acc,pathnum,npaths,
		    score,chromosome_iit,queryseq,pairedlength,
		    snps_iit,snps_divint_crosstable,resulttype,first_read_p,
		    npaths_mate,queryseq_mate,quality_shift);
  } else if (hittype == DEL) {
    print_deletion(this,mate,acc,pathnum,npaths,
		   score,chromosome_iit,queryseq,pairedlength,
		   snps_iit,snps_divint_crosstable,resulttype,first_read_p,
		   npaths_mate,queryseq_mate,quality_shift);

  } else if (hittype == SPLICE) {
    /* Follows print_splice_distance() in substring.c */
    donor = Stage3_substring1(this);
    acceptor = Stage3_substring2(this);

    if (donor == NULL || acceptor == NULL) {
      abort();
    } else if (Stage3_distance(this) == 0U) {
      /* translocation */
      print_exon_exon(this,mate,acc,pathnum,npaths,
		      score,chromosome_iit,queryseq,pairedlength,
		      snps_iit,snps_divint_crosstable,resulttype,first_read_p,
		      npaths_mate,queryseq_mate,quality_shift);
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
	print_localsplice(this,mate,acc,pathnum,npaths,
			  score,chromosome_iit,queryseq,pairedlength,
			  snps_iit,snps_divint_crosstable,resulttype,first_read_p,
			  npaths_mate,queryseq_mate,quality_shift);
      } else {
	print_exon_exon(this,mate,acc,pathnum,npaths,
			score,chromosome_iit,queryseq,pairedlength,
			snps_iit,snps_divint_crosstable,resulttype,first_read_p,
			npaths_mate,queryseq_mate,quality_shift);
      }
    }

  } else if (hittype == TERMINAL) {
    print_terminal(this,mate,acc,pathnum,npaths,
		   score,chromosome_iit,queryseq,pairedlength,
		   snps_iit,snps_divint_crosstable,resulttype,first_read_p,
		   npaths_mate,queryseq_mate,quality_shift);

  } else {
    abort();
  }

  return;
}



static void
print_paired_hits (Result_T result, IIT_T chromosome_iit,
		   Sequence_T queryseq1, Sequence_T queryseq2,
		   IIT_T snps_iit, int *snps_divint_crosstable,
		   IIT_T splicesites_iit, int donor_typeint, int acceptor_typeint,
		   int *splicesites_divint_crosstable, int maxpaths, bool quiet_if_excessive_p,
		   int quality_shift) {
  Stage3pair_T *stage3pairarray, stage3pair;
  Stage3_T *stage3array1, *stage3array2, stage3, mate;
  Resulttype_T resulttype;
  int npaths, npaths1, npaths2, pathnum;
  char *acc;

  acc = Sequence_accession(queryseq1);

  resulttype = Result_resulttype(result);
  if (resulttype == PAIREDEND_CONCORDANT || resulttype == PAIREDEND_SAMECHR_SINGLE || resulttype == PAIREDEND_SAMECHR_MULTIPLE) {
    /* Information is in pairs */
    stage3pairarray = (Stage3pair_T *) Result_array(&npaths,result);

    if (npaths == 0) {
      abort();

    } else if (quiet_if_excessive_p && npaths > maxpaths) {
      /* Not sure what to do here.  Print as nomapping?.  Leaving as blank for now. */

    } else {
      for (pathnum = 1; pathnum <= npaths && pathnum <= maxpaths; pathnum++) {
	stage3pair = stage3pairarray[pathnum-1];

	/* print 5', sequence is not inverted */
	SAM_print(Stage3pair_hit5(stage3pair),/*mate*/Stage3pair_hit3(stage3pair),acc,
		  pathnum,npaths,Stage3pair_hit5_score(stage3pair),
		  chromosome_iit,/*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
		  snps_iit,snps_divint_crosstable,
		  splicesites_iit,donor_typeint,acceptor_typeint,splicesites_divint_crosstable,
		  Stage3pair_pairlength(stage3pair),resulttype,
		  /*first_read_p*/true,/*npaths_mate*/npaths,quality_shift);

	/* print 3', sequence is inverted */
	SAM_print(Stage3pair_hit3(stage3pair),/*mate*/Stage3pair_hit5(stage3pair),acc,
		  pathnum,npaths,Stage3pair_hit3_score(stage3pair),
		  chromosome_iit,/*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
		  snps_iit,snps_divint_crosstable,
		  splicesites_iit,donor_typeint,acceptor_typeint,splicesites_divint_crosstable,
		  Stage3pair_pairlength(stage3pair),resulttype,
		  /*first_read_p*/false,/*npaths_mate*/npaths,quality_shift);
      }
    }

  } else if (resulttype == PAIREDEND_AS_SINGLES_UNIQUE) {
    /* BAM prints mate information in this situation */
    stage3array1 = (Stage3_T *) Result_array(&npaths1,result);
    stage3array2 = (Stage3_T *) Result_array2(&npaths2,result);

    /* print 5', sequence is not inverted */
    SAM_print(stage3array1[0],/*mate*/stage3array2[0],acc,
	      pathnum,/*npaths*/1,Stage3_score(stage3array1[0]),
	      chromosome_iit,/*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
	      snps_iit,snps_divint_crosstable,
	      splicesites_iit,donor_typeint,acceptor_typeint,splicesites_divint_crosstable,
	      /*pairedlength*/0U,/*resulttype*/PAIREDEND_AS_SINGLES_UNIQUE,
	      /*first_read_p*/true,/*npaths_mate*/1,quality_shift);

    /* print 3', sequence is inverted */
    SAM_print(stage3array2[0],/*mate*/stage3array1[0],acc,
	      pathnum,/*npaths*/1,Stage3_score(stage3array2[0]),
	      chromosome_iit,/*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
	      snps_iit,snps_divint_crosstable,
	      splicesites_iit,donor_typeint,acceptor_typeint,splicesites_divint_crosstable,
	      /*pairedlength*/0U,/*resulttype*/PAIREDEND_AS_SINGLES_UNIQUE,
	      /*first_read_p*/false,/*npaths_mate*/1,quality_shift);

  } else if (resulttype == PAIREDEND_AS_SINGLES) {
    stage3array1 = (Stage3_T *) Result_array(&npaths1,result);
    stage3array2 = (Stage3_T *) Result_array2(&npaths2,result);

    /* print 5' results, sequence is not inverted */
    if (npaths1 == 0 || (quiet_if_excessive_p && npaths1 > maxpaths)) {
      if (npaths2 == 1) {
	mate = stage3array2[0];
      } else {
	mate = (Stage3_T) NULL;
      }
      SAM_print_nomapping(queryseq1,mate,acc,chromosome_iit,resulttype,
			  /*first_read_p*/true,/*npaths_mate*/npaths2,
			  /*queryseq_mate*/queryseq2,quality_shift);
    } else {
      for (pathnum = 1; pathnum <= npaths1 && pathnum <= maxpaths; pathnum++) {
	stage3 = stage3array1[pathnum-1];
	SAM_print(stage3,/*mate*/NULL,acc,pathnum,npaths1,Stage3_score(stage3),
		  chromosome_iit,/*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
		  snps_iit,snps_divint_crosstable,
		  splicesites_iit,donor_typeint,acceptor_typeint,splicesites_divint_crosstable,
		  /*pairedlength*/0U,resulttype,/*first_read_p*/true,/*npaths_mate*/npaths2,
		  quality_shift);
      }
    }
			  
    /* print 3', sequence is inverted */
    if (npaths2 == 0 || (quiet_if_excessive_p && npaths2 > maxpaths)) {
      if (npaths1 == 1) {
	mate = stage3array1[0];
      } else {
	mate = (Stage3_T) NULL;
      }
      SAM_print_nomapping(queryseq2,mate,acc,chromosome_iit,resulttype,
			  /*first_read_p*/false,/*npaths_mate*/npaths1,
			  /*queryseq_mate*/queryseq1,quality_shift);
    } else {
      for (pathnum = 1; pathnum <= npaths2 && pathnum <= maxpaths; pathnum++) {
	stage3 = stage3array2[pathnum-1];
	SAM_print(stage3,/*mate*/NULL,acc,pathnum,npaths2,Stage3_score(stage3),
		  chromosome_iit,/*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
		  snps_iit,snps_divint_crosstable,
		  splicesites_iit,donor_typeint,acceptor_typeint,splicesites_divint_crosstable,
		  /*pairedlength*/0U,resulttype,/*first_read_p*/false,/*npaths_mate*/npaths1,
		  quality_shift);
      }
    }

  } else {
    abort();
  }

  return;
}



void
SAM_print_paired (Result_T result, IIT_T chromosome_iit, Sequence_T queryseq1, Sequence_T queryseq2,
		  IIT_T snps_iit, int *snps_divint_crosstable,
		  IIT_T splicesites_iit, int donor_typeint, int acceptor_typeint,
		  int *splicesites_divint_crosstable, int maxpaths, bool quiet_if_excessive_p,
		  int quality_shift, bool circularp) {

  if (circularp == false) {
    print_paired_hits(result,chromosome_iit,queryseq1,queryseq2,
		      snps_iit,snps_divint_crosstable,
		      splicesites_iit,donor_typeint,acceptor_typeint,splicesites_divint_crosstable,
		      maxpaths,quiet_if_excessive_p,quality_shift);
  } else {
    fprintf(stderr,"SAM output not yet implemented for circular paired-end reads\n");
    exit(9);

#if 0
    /* First sequence.  query2 has already been revcomp'd */
    print_paired_hits(result,'>',/*fivep*/false,
		      chromosome_iit,queryseq2,/*headerseq*/queryseq2,
		      snps_iit,snps_divint_crosstable,
		      splicesites_iit,donor_typeint,acceptor_typeint,splicesites_divint_crosstable,
		      maxpaths,quiet_if_excessive_p,/*invertp*/false);

    /* Second sequence */
    print_paired_hits(result,'<',/*fivep*/true,
		      chromosome_iit,queryseq1,/*headerseq*/queryseq2,
		      snps_iit,snps_divint_crosstable,
		      splicesites_iit,donor_typeint,acceptor_typeint,splicesites_divint_crosstable,
		      maxpaths,quiet_if_excessive_p,/*invertp*/true);
#endif
  }

  return;
}

