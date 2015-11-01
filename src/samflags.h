/* $Id: samflags.h 106634 2013-09-03 17:01:14Z twu $ */
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
#define ABBREV_NOMAPPING_1 "N1"
#define ABBREV_NOMAPPING_2 "N2"
#define ABBREV_HALFMAPPING_UNIQ "HU"
#define ABBREV_HALFMAPPING_CIRCULAR "HC"
#define ABBREV_HALFMAPPING_TRANSLOC "HT"
#define ABBREV_HALFMAPPING_MULT "HM"
#define ABBREV_UNPAIRED_UNIQ "UU"
#define ABBREV_UNPAIRED_CIRCULAR "UC"
#define ABBREV_UNPAIRED_TRANSLOC "UT"
#define ABBREV_UNPAIRED_MULT "UM"
#define ABBREV_PAIRED_UNIQ_CIRCULAR "PC"
#define ABBREV_PAIRED_UNIQ_INV "PI"
#define ABBREV_PAIRED_UNIQ_SCR "PS"
#define ABBREV_PAIRED_UNIQ_LONG "PL"
#define ABBREV_PAIRED_MULT "PM"
#define ABBREV_CONCORDANT_UNIQ "CU"
#define ABBREV_CONCORDANT_CIRCULAR "CC"
#define ABBREV_CONCORDANT_TRANSLOC "CT"
#define ABBREV_CONCORDANT_MULT "CM"


#endif

