/* $Id: samflags.h 154089 2014-11-25 21:03:16Z twu $ */
#ifndef SAMFLAGS_INCLUDED
#define SAMFLAGS_INCLUDED

#define PAIRED_READ        0x0001 /* 1 */
#define PAIRED_MAPPING     0x0002 /* 2 */
#define QUERY_UNMAPPED     0x0004 /* 4 */
#define MATE_UNMAPPED      0x0008 /* 8 */
#define QUERY_MINUSP       0x0010 /* 16 */
#define MATE_MINUSP        0x0020 /* 32 */
#define FIRST_READ_P       0x0040 /* 64 */
#define SECOND_READ_P      0x0080 /* 128 */
#define NOT_PRIMARY        0x0100 /* 256 */
#define BAD_READ_QUALITY   0x0200 /* 512 */
#define DUPLICATE_READ     0x0400 /* 1024 */


/* 83 = first read, minus strand for paired */
/* 99 = first read, plus strand for paired */
/* 147 = second read, minus strand for paired */
/* 163 = second read, plus strand for paired */

/* For forcing a read to be primary */
#define SET_PRIMARY        0xFEFF /* do a logical-and (&) with this */


/* XO tag for output type */
#define ABBREV_NOMAPPING_1 "NM"
#define ABBREV_NOMAPPING_2 "NM"
#define ABBREV_HALFMAPPING_UNIQ "HU"
#define ABBREV_HALFMAPPING_CIRCULAR "HC"
#define ABBREV_HALFMAPPING_TRANSLOC "HT"
#define ABBREV_HALFMAPPING_MULT "HM"
#define ABBREV_HALFMAPPING_MULT_XS "HX"
#define ABBREV_UNPAIRED_UNIQ "UU"
#define ABBREV_UNPAIRED_CIRCULAR "UC"
#define ABBREV_UNPAIRED_TRANSLOC "UT"
#define ABBREV_UNPAIRED_MULT "UM"
#define ABBREV_UNPAIRED_MULT_XS "UX"
#define ABBREV_PAIRED_UNIQ_CIRCULAR "PC"
#define ABBREV_PAIRED_UNIQ_INV "PI"
#define ABBREV_PAIRED_UNIQ_SCR "PS"
#define ABBREV_PAIRED_UNIQ_LONG "PL"
#define ABBREV_PAIRED_MULT "PM"
#define ABBREV_PAIRED_MULT_XS "PX"
#define ABBREV_CONCORDANT_UNIQ "CU"
#define ABBREV_CONCORDANT_CIRCULAR "CC"
#define ABBREV_CONCORDANT_TRANSLOC "CT"
#define ABBREV_CONCORDANT_MULT "CM"
#define ABBREV_CONCORDANT_MULT_XS "CX"

typedef enum {OUTPUT_NONE,
	      OUTPUT_NM,	/* nomapping */
	      OUTPUT_HU,	/* halfmapping_uniq */
	      OUTPUT_HC,	/* halfmapping_circular */
	      OUTPUT_HT,	/* halfmapping_transloc */
	      OUTPUT_HM,	/* halfmapping_mult */
	      OUTPUT_HX,	/* halfmapping_mult_xs */
	      OUTPUT_UU,	/* unpaired_uniq */
	      OUTPUT_UC,	/* unpaired_circular */
	      OUTPUT_UT,	/* unpaired_transloc */
	      OUTPUT_UM,	/* unpaired_mult */
	      OUTPUT_UX,	/* unpaired_mult_xs */
	      OUTPUT_PC,	/* paired_uniq_circular */
	      OUTPUT_PI,	/* paired_uniq_inv */
	      OUTPUT_PS,	/* paired_uniq_scr */
	      OUTPUT_PL,	/* paired_uniq_long */
	      OUTPUT_PM,	/* paired_mult */
	      OUTPUT_PX,	/* paired_mult_xs */
	      OUTPUT_CU,	/* concordant_uniq */
	      OUTPUT_CC,	/* concordant_circular */
	      OUTPUT_CT,	/* concordant_transloc */
	      OUTPUT_CM,	/* concordant_mult */
	      OUTPUT_CX}	/* concordant_mult_xs */
  SAM_split_output_type;


#endif

