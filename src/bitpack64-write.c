static char rcsid[] = "$Id: bitpack64-write.c 109823 2013-10-02 22:28:36Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "bitpack64-write.h"

#ifdef WORDS_BIGENDIAN
#include "bigendian.h"		/* For FWRITE_UINTS */
#else
#include "littleendian.h"	/* For FWRITE_UINTS */
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>		/* For memset */

#ifdef HAVE_SSE2
#include <emmintrin.h>
#endif

/* #define ALLOW_ODD_PACKSIZES 1 */

/* Note: For offset pointers, where we need fast cumulative sums, we
   use vertical format (where successive values are in different
   packed unsigned ints).  For lcp, we want raw values, and vertical
   format is still slightly more efficient than horizontal format. */


#ifdef HAVE_SSE2
static int
write_reg_buffered_vert (FILE *offsetscomp_fp, Positionsptr_T *offsets_buffer,
			 int offsets_buffer_size, int offsets_buffer_i, __m128i OutReg) {

#if 0
  /* Type casting method (when we passed in pointer to OutReg).  Needs a memory fence. */
  UINT4 *buffer = (UINT4 *) OutReg;
  _mm_lfence();  /* Needed to avoid storing incorrect values into offsets_buffer */
#else
  /* Storing method.  Safer.  */
  UINT4 buffer[4];
  _mm_store_si128((__m128i *) buffer,OutReg);
#endif

  /* printf("Writing %08X %08X %08X %08X\n",buffer[0],buffer[1],buffer[2],buffer[3]); */

  offsets_buffer[offsets_buffer_i++] = buffer[0];
  if (offsets_buffer_i == offsets_buffer_size) {
    FWRITE_UINTS(offsets_buffer,offsets_buffer_size,offsetscomp_fp);
    offsets_buffer_i = 0;
  }

  offsets_buffer[offsets_buffer_i++] = buffer[1];
  if (offsets_buffer_i == offsets_buffer_size) {
    FWRITE_UINTS(offsets_buffer,offsets_buffer_size,offsetscomp_fp);
    offsets_buffer_i = 0;
  }

  offsets_buffer[offsets_buffer_i++] = buffer[2];
  if (offsets_buffer_i == offsets_buffer_size) {
    FWRITE_UINTS(offsets_buffer,offsets_buffer_size,offsetscomp_fp);
    offsets_buffer_i = 0;
  }

  offsets_buffer[offsets_buffer_i++] = buffer[3];
  if (offsets_buffer_i == offsets_buffer_size) {
    FWRITE_UINTS(offsets_buffer,offsets_buffer_size,offsetscomp_fp);
    offsets_buffer_i = 0;
  }

  return offsets_buffer_i;
}
#else
static int
write_reg_buffered_vert (FILE *offsetscomp_fp, Positionsptr_T *offsets_buffer,
			 int offsets_buffer_size, int offsets_buffer_i,
			 UINT4 *horizontal, int nwritten) {
  UINT4 vertical[64];
  int nrows = nwritten/4, row, column, k;

  /* Convert to horizontal */
  for (column = 0; column < 4; column++) {
    k = column;
    for (row = 0; row < nrows; row++) {
      vertical[k] = *horizontal++;
      k += 4;
    }
  }
    
  /* Send to output buffer */
  for (k = 0; k < nwritten; k++) {
    /* printf("Writing %08X\n",vertical[k]); */
    offsets_buffer[offsets_buffer_i++] = vertical[k];
    if (offsets_buffer_i == offsets_buffer_size) {
      FWRITE_UINTS(offsets_buffer,offsets_buffer_size,offsetscomp_fp);
      offsets_buffer_i = 0;
    }
  }

  return offsets_buffer_i;
}
#endif

static int
write_reg_buffered_horiz (FILE *offsetscomp_fp, Positionsptr_T *offsets_buffer,
			  int offsets_buffer_size, int offsets_buffer_i,
			  UINT4 *values, int nwritten) {
  int k;

  /* Send to output buffer */
  for (k = 0; k < nwritten; k++) {
    /* printf("Writing %08X\n",values[k]); */
    offsets_buffer[offsets_buffer_i++] = values[k];
    if (offsets_buffer_i == offsets_buffer_size) {
      FWRITE_UINTS(offsets_buffer,offsets_buffer_size,offsetscomp_fp);
      offsets_buffer_i = 0;
    }
  }

  return offsets_buffer_i;
}




#ifdef HAVE_SSE2
static __m128i mask1, mask2, mask3, mask4, mask5, mask6, mask7, mask8,
  mask9, mask10, mask11, mask12, mask13, mask14, mask15, mask16,
  mask17, mask18, mask19, mask20, mask21, mask22, mask23, mask24,
  mask25, mask26, mask27, mask28, mask29, mask30, mask31;
#endif


void
Bitpack64_write_setup () {

#ifdef HAVE_SSE2
  mask1 = _mm_set1_epi32(1U);
  mask2 = _mm_set1_epi32(3U);
  mask3 =  _mm_set1_epi32(7U);
  mask4 =  _mm_set1_epi32(15U);
  mask5 =  _mm_set1_epi32(31U);
  mask6 =  _mm_set1_epi32(63U);
  mask7 =  _mm_set1_epi32(127U);
  mask8 =  _mm_set1_epi32(255U);
  mask9 =  _mm_set1_epi32(511U);
  mask10 =  _mm_set1_epi32(1023U);
  mask11 =  _mm_set1_epi32(2047U);
  mask12 =  _mm_set1_epi32(4095U);
  mask13 =  _mm_set1_epi32(8191U);
  mask14 =  _mm_set1_epi32(16383U);
  mask15 =  _mm_set1_epi32(32767U);
  mask16 =  _mm_set1_epi32(65535U);
  mask17 =  _mm_set1_epi32(131071U);
  mask18 =  _mm_set1_epi32(262143U);
  mask19 =  _mm_set1_epi32(524287U);
  mask20 =  _mm_set1_epi32(1048575U);
  mask21 =  _mm_set1_epi32(2097151U);
  mask22 =  _mm_set1_epi32(4194303U);
  mask23 =  _mm_set1_epi32(8388607U);
  mask24 =  _mm_set1_epi32(16777215U);
  mask25 =  _mm_set1_epi32(33554431U);
  mask26 =  _mm_set1_epi32(67108863U);
  mask27 =  _mm_set1_epi32(134217727U);
  mask28 =  _mm_set1_epi32(268435455U);
  mask29 =  _mm_set1_epi32(536870911U);
  mask30 =  _mm_set1_epi32(1073741823U);
  mask31 =  _mm_set1_epi32(2147483647U);
#endif

  return;
}

#ifdef ALLOW_ODD_PACKSIZES
/* nwritten = 1 * 4 = 4 */
static int
write_01_vert (FILE *offsetscomp_fp, Positionsptr_T *offsets_buffer, int offsets_buffer_size, int offsets_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask1);
    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask1);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 1));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask1);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 2));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask1);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 3));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask1);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 4));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask1);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 5));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask1);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 6));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask1);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 7));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask1);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask1);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 9));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask1);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 10));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask1);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 11));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask1);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 12));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask1);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 13));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask1);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 14));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask1);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 15));
    offsets_buffer_i = write_reg_buffered(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					  OutReg);

    return offsets_buffer_i;
}
#endif


#ifdef HAVE_SSE2
/* nwritten = 1 * 4 = 4 */
static int
write_02_vert (FILE *offsetscomp_fp, Positionsptr_T *offsets_buffer, int offsets_buffer_size, int offsets_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask2);
    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask2);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 2));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask2);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 4));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask2);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 6));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask2);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask2);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 10));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask2);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 12));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask2);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 14));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask2);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask2);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 18));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask2);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 20));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask2);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 22));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask2);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask2);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 26));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask2);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 28));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask2);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 30));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    return offsets_buffer_i;
}
#endif

static int
write_02_horiz (UINT4 *out, const UINT4 *in) {
  int column;

  for (column = 0; column < 4; column++) {
    *out |= (*in)   % (1U << 2 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 2 )  ) <<  2 ;
    ++in;
    *out |= ( (*in)   % (1U << 2 )  ) <<  4 ;
    ++in;
    *out |= ( (*in)   % (1U << 2 )  ) <<  6 ;
    ++in;
    *out |= ( (*in)   % (1U << 2 )  ) <<  8 ;
    ++in;
    *out |= ( (*in)   % (1U << 2 )  ) <<  10 ;
    ++in;
    *out |= ( (*in)   % (1U << 2 )  ) <<  12 ;
    ++in;
    *out |= ( (*in)   % (1U << 2 )  ) <<  14 ;
    ++in;
    *out |= ( (*in)   % (1U << 2 )  ) <<  16 ;
    ++in;
    *out |= ( (*in)   % (1U << 2 )  ) <<  18 ;
    ++in;
    *out |= ( (*in)   % (1U << 2 )  ) <<  20 ;
    ++in;
    *out |= ( (*in)   % (1U << 2 )  ) <<  22 ;
    ++in;
    *out |= ( (*in)   % (1U << 2 )  ) <<  24 ;
    ++in;
    *out |= ( (*in)   % (1U << 2 )  ) <<  26 ;
    ++in;
    *out |= ( (*in)   % (1U << 2 )  ) <<  28 ;
    ++in;
    *out |= ( (*in)   % (1U << 2 )  ) <<  30 ;
    ++out;
    ++in;
  }

  return 4;
}



#ifdef ALLOW_ODD_PACKSIZES
/* nwritten = 2 * 4 = 8 */
static int
write_03_vert (FILE *offsetscomp_fp, Positionsptr_T *offsets_buffer, int offsets_buffer_size, int offsets_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask3);
    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask3);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 3));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask3);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 6));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask3);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 9));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask3);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 12));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask3);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 15));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask3);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 18));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask3);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 21));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask3);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask3);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 27));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask3);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 30));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 3 - 1);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask3);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 1));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask3);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 4));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask3);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 7));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask3);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 10));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask3);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 13));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    return offsets_buffer_i;
}
#endif



#ifdef ALLOW_ODD_PACKSIZES
/* nwritten = 3 * 4 = 12 */
static int
write_05_vert (FILE *offsetscomp_fp, Positionsptr_T *offsets_buffer, int offsets_buffer_size, int offsets_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask5);
    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask5);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 5));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask5);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 10));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask5);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 15));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask5);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 20));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask5);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 25));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask5);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 30));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 5 - 3);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask5);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 3));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask5);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask5);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 13));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask5);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 18));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask5);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 23));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask5);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 28));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 5 - 1);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask5);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 1));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask5);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 6));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask5);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 11));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    return offsets_buffer_i;
}
#endif


#ifdef HAVE_SSE2
/* nwritten = 3 * 4 = 12 */
static int
write_06_vert (FILE *offsetscomp_fp, Positionsptr_T *offsets_buffer, int offsets_buffer_size, int offsets_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask6);
    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask6);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 6));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask6);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 12));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask6);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 18));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask6);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask6);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 30));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 6 - 4);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask6);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 4));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask6);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 10));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask6);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask6);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 22));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask6);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 28));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 6 - 2);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask6);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 2));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask6);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask6);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 14));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask6);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 20));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask6);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 26));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    return offsets_buffer_i;
}
#endif

static int
write_06_horiz (UINT4 *out, const UINT4 *in) {
  int column;

  for (column = 0; column < 4; column++) {
    *out |= (*in)   % (1U << 6 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 6 )  ) <<  6 ;
    ++in;
    *out |= ( (*in)   % (1U << 6 )  ) <<  12 ;
    ++in;
    *out |= ( (*in)   % (1U << 6 )  ) <<  18 ;
    ++in;
    *out |= ( (*in)   % (1U << 6 )  ) <<  24 ;
    ++in;
    *out |= ( (*in)   % (1U << 6 )  ) <<  30 ;
    ++out;
    *out |=  ( (*in)   % (1U << 6 ) ) >> ( 6  -  4 );
    ++in;
    *out |= ( (*in)   % (1U << 6 )  ) <<  4 ;
    ++in;
    *out |= ( (*in)   % (1U << 6 )  ) <<  10 ;
    ++in;
    *out |= ( (*in)   % (1U << 6 )  ) <<  16 ;
    ++in;
    *out |= ( (*in)   % (1U << 6 )  ) <<  22 ;
    ++in;
    *out |= ( (*in)   % (1U << 6 )  ) <<  28 ;
    ++out;
    *out |=  ( (*in)   % (1U << 6 ) ) >> ( 6  -  2 );
    ++in;
    *out |= ( (*in)   % (1U << 6 )  ) <<  2 ;
    ++in;
    *out |= ( (*in)   % (1U << 6 )  ) <<  8 ;
    ++in;
    *out |= ( (*in)   % (1U << 6 )  ) <<  14 ;
    ++in;
    *out |= ( (*in)   % (1U << 6 )  ) <<  20 ;
    ++in;
    *out |= ( (*in)   % (1U << 6 )  ) <<  26 ;
    ++out;
    ++in;
  }

  return 12;
}



#ifdef ALLOW_ODD_PACKSIZES
/* nwritten = 4 * 4 = 16 */
static int
write_07_vert (FILE *offsetscomp_fp, Positionsptr_T *offsets_buffer, int offsets_buffer_size, int offsets_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask7);
    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask7);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 7));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask7);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 14));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask7);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 21));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask7);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 28));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 7 - 3);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask7);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 3));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask7);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 10));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask7);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 17));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask7);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask7);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 31));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 7 - 6);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask7);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 6));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask7);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 13));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask7);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 20));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask7);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 27));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 7 - 2);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask7);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 2));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask7);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 9));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    return offsets_buffer_i;
}
#endif


#ifdef ALLOW_ODD_PACKSIZES
/* nwritten = 5 * 4 = 20 */
static int
write_09_vert (FILE *offsetscomp_fp, Positionsptr_T *offsets_buffer, int offsets_buffer_size, int offsets_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask9);
    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask9);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 9));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask9);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 18));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask9);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 27));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 9 - 4);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask9);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 4));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask9);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 13));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask9);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 22));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask9);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 31));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 9 - 8);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask9);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask9);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 17));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask9);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 26));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 9 - 3);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask9);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 3));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask9);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 12));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask9);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 21));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask9);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 30));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 9 - 7);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask9);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 7));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    return offsets_buffer_i;
}
#endif


#ifdef HAVE_SSE2
/* nwritten = 5 * 4 = 20 */
static int
write_10_vert (FILE *offsetscomp_fp, Positionsptr_T *offsets_buffer, int offsets_buffer_size, int offsets_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask10);
    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask10);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 10));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask10);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 20));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask10);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 30));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 10 - 8);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask10);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask10);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 18));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask10);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 28));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 10 - 6);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask10);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 6));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask10);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask10);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 26));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 10 - 4);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask10);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 4));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask10);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 14));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask10);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 10 - 2);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask10);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 2));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask10);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 12));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask10);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 22));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    return offsets_buffer_i;
}
#endif

static int
write_10_horiz (UINT4 *out, const UINT4 *in) {
  int column;

  for (column = 0; column < 4; column++) {
    *out |= (*in)   % (1U << 10 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 10 )  ) <<  10 ;
    ++in;
    *out |= ( (*in)   % (1U << 10 )  ) <<  20 ;
    ++in;
    *out |= ( (*in)   % (1U << 10 )  ) <<  30 ;
    ++out;
    *out |=  ( (*in)   % (1U << 10 ) ) >> ( 10  -  8 );
    ++in;
    *out |= ( (*in)   % (1U << 10 )  ) <<  8 ;
    ++in;
    *out |= ( (*in)   % (1U << 10 )  ) <<  18 ;
    ++in;
    *out |= ( (*in)   % (1U << 10 )  ) <<  28 ;
    ++out;
    *out |=  ( (*in)   % (1U << 10 ) ) >> ( 10  -  6 );
    ++in;
    *out |= ( (*in)   % (1U << 10 )  ) <<  6 ;
    ++in;
    *out |= ( (*in)   % (1U << 10 )  ) <<  16 ;
    ++in;
    *out |= ( (*in)   % (1U << 10 )  ) <<  26 ;
    ++out;
    *out |=  ( (*in)   % (1U << 10 ) ) >> ( 10  -  4 );
    ++in;
    *out |= ( (*in)   % (1U << 10 )  ) <<  4 ;
    ++in;
    *out |= ( (*in)   % (1U << 10 )  ) <<  14 ;
    ++in;
    *out |= ( (*in)   % (1U << 10 )  ) <<  24 ;
    ++out;
    *out |=  ( (*in)   % (1U << 10 ) ) >> ( 10  -  2 );
    ++in;
    *out |= ( (*in)   % (1U << 10 )  ) <<  2 ;
    ++in;
    *out |= ( (*in)   % (1U << 10 )  ) <<  12 ;
    ++in;
    *out |= ( (*in)   % (1U << 10 )  ) <<  22 ;
    ++out;
    ++in;
  }

  return 20;
}




#ifdef ALLOW_ODD_PACKSIZES
/* nwritten = 6 * 4 = 24 */
static int
write_11_vert (FILE *offsetscomp_fp, Positionsptr_T *offsets_buffer, int offsets_buffer_size, int offsets_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask11);
    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask11);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 11));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask11);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 22));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 11 - 1);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask11);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 1));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask11);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 12));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask11);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 23));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 11 - 2);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask11);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 2));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask11);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 13));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask11);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 11 - 3);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask11);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 3));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask11);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 14));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask11);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 25));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 11 - 4);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask11);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 4));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask11);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 15));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask11);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 26));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 11 - 5);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask11);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 5));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    return offsets_buffer_i;
}
#endif


#ifdef HAVE_SSE2
/* nwritten = 6 * 4 = 24 */
static int
write_12_vert (FILE *offsetscomp_fp, Positionsptr_T *offsets_buffer, int offsets_buffer_size, int offsets_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask12);
    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask12);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 12));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask12);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 12 - 4);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask12);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 4));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask12);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask12);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 28));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 12 - 8);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask12);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask12);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 20));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    InReg = _mm_and_si128(_mm_load_si128(++in), mask12);

    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask12);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 12));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask12);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 12 - 4);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask12);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 4));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask12);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask12);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 28));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 12 - 8);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask12);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask12);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 20));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    return offsets_buffer_i;
}
#endif

static int
write_12_horiz (UINT4 *out, const UINT4 *in) {
  int column;

  for (column = 0; column < 4; column++) {

    *out |= (*in)   % (1U << 12 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 12 )  ) <<  12 ;
    ++in;
    *out |= ( (*in)   % (1U << 12 )  ) <<  24 ;
    ++out;
    *out |=  ( (*in)   % (1U << 12 ) ) >> ( 12  -  4 );
    ++in;
    *out |= ( (*in)   % (1U << 12 )  ) <<  4 ;
    ++in;
    *out |= ( (*in)   % (1U << 12 )  ) <<  16 ;
    ++in;
    *out |= ( (*in)   % (1U << 12 )  ) <<  28 ;
    ++out;
    *out |=  ( (*in)   % (1U << 12 ) ) >> ( 12  -  8 );
    ++in;
    *out |= ( (*in)   % (1U << 12 )  ) <<  8 ;
    ++in;
    *out |= ( (*in)   % (1U << 12 )  ) <<  20 ;
    ++out;
    ++in;
    *out |= (*in)   % (1U << 12 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 12 )  ) <<  12 ;
    ++in;
    *out |= ( (*in)   % (1U << 12 )  ) <<  24 ;
    ++out;
    *out |=  ( (*in)   % (1U << 12 ) ) >> ( 12  -  4 );
    ++in;
    *out |= ( (*in)   % (1U << 12 )  ) <<  4 ;
    ++in;
    *out |= ( (*in)   % (1U << 12 )  ) <<  16 ;
    ++in;
    *out |= ( (*in)   % (1U << 12 )  ) <<  28 ;
    ++out;
    *out |=  ( (*in)   % (1U << 12 ) ) >> ( 12  -  8 );
    ++in;
    *out |= ( (*in)   % (1U << 12 )  ) <<  8 ;
    ++in;
    *out |= ( (*in)   % (1U << 12 )  ) <<  20 ;
    ++out;
    ++in;
  }

  return 24;
}




#ifdef ALLOW_ODD_PACKSIZES
/* nwritten = 7 * 4 = 28 */
static int
write_13_vert (FILE *offsetscomp_fp, Positionsptr_T *offsets_buffer, int offsets_buffer_size, int offsets_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask13);
    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask13);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 13));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask13);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 26));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 13 - 7);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask13);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 7));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask13);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 20));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 13 - 1);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask13);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 1));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask13);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 14));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask13);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 27));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 13 - 8);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask13);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask13);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 21));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 13 - 2);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask13);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 2));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask13);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 15));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask13);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 28));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 13 - 9);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask13);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 9));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask13);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 22));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 13 - 3);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask13);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 3));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    return offsets_buffer_i;
}
#endif


#ifdef HAVE_SSE2
/* nwritten = 7 * 4 = 28 */
static int
write_14_vert (FILE *offsetscomp_fp, Positionsptr_T *offsets_buffer, int offsets_buffer_size, int offsets_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask14);
    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask14);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 14));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask14);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 28));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 14 - 10);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask14);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 10));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask14);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 14 - 6);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask14);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 6));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask14);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 20));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 14 - 2);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask14);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 2));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask14);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask14);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 30));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 14 - 12);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask14);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 12));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask14);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 26));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 14 - 8);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask14);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask14);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 22));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 14 - 4);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask14);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 4));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask14);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 18));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    return offsets_buffer_i;
}
#endif


static int
write_14_horiz (UINT4 *out, const UINT4 *in) {
  int column;

  for (column = 0; column < 4; column++) {
    *out |= (*in)   % (1U << 14 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 14 )  ) <<  14 ;
    ++in;
    *out |= ( (*in)   % (1U << 14 )  ) <<  28 ;
    ++out;
    *out |=  ( (*in)   % (1U << 14 ) ) >> ( 14  -  10 );
    ++in;
    *out |= ( (*in)   % (1U << 14 )  ) <<  10 ;
    ++in;
    *out |= ( (*in)   % (1U << 14 )  ) <<  24 ;
    ++out;
    *out |=  ( (*in)   % (1U << 14 ) ) >> ( 14  -  6 );
    ++in;
    *out |= ( (*in)   % (1U << 14 )  ) <<  6 ;
    ++in;
    *out |= ( (*in)   % (1U << 14 )  ) <<  20 ;
    ++out;
    *out |=  ( (*in)   % (1U << 14 ) ) >> ( 14  -  2 );
    ++in;
    *out |= ( (*in)   % (1U << 14 )  ) <<  2 ;
    ++in;
    *out |= ( (*in)   % (1U << 14 )  ) <<  16 ;
    ++in;
    *out |= ( (*in)   % (1U << 14 )  ) <<  30 ;
    ++out;
    *out |=  ( (*in)   % (1U << 14 ) ) >> ( 14  -  12 );
    ++in;
    *out |= ( (*in)   % (1U << 14 )  ) <<  12 ;
    ++in;
    *out |= ( (*in)   % (1U << 14 )  ) <<  26 ;
    ++out;
    *out |=  ( (*in)   % (1U << 14 ) ) >> ( 14  -  8 );
    ++in;
    *out |= ( (*in)   % (1U << 14 )  ) <<  8 ;
    ++in;
    *out |= ( (*in)   % (1U << 14 )  ) <<  22 ;
    ++out;
    *out |=  ( (*in)   % (1U << 14 ) ) >> ( 14  -  4 );
    ++in;
    *out |= ( (*in)   % (1U << 14 )  ) <<  4 ;
    ++in;
    *out |= ( (*in)   % (1U << 14 )  ) <<  18 ;
    ++out;
    ++in;
  }

  return 28;
}


#ifdef ALLOW_ODD_PACKSIZES
/* nwritten = 8 * 4 = 32 */
static int
write_15_vert (FILE *offsetscomp_fp, Positionsptr_T *offsets_buffer, int offsets_buffer_size, int offsets_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask15);
    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask15);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 15));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask15);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 30));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 15 - 13);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask15);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 13));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask15);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 28));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 15 - 11);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask15);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 11));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask15);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 26));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 15 - 9);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask15);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 9));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask15);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 15 - 7);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask15);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 7));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask15);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 22));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 15 - 5);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask15);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 5));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask15);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 20));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 15 - 3);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask15);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 3));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask15);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 18));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 15 - 1);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask15);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 1));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    return offsets_buffer_i;
}
#endif



#ifdef ALLOW_ODD_PACKSIZES
/* nwritten = 9 * 4 = 36 */
static int
write_17_vert (FILE *offsetscomp_fp, Positionsptr_T *offsets_buffer, int offsets_buffer_size, int offsets_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask17);
    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask17);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 17));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 17 - 2);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask17);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 2));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask17);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 19));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 17 - 4);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask17);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 4));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask17);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 21));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 17 - 6);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask17);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 6));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask17);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 23));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 17 - 8);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask17);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask17);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 25));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 17 - 10);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask17);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 10));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask17);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 27));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 17 - 12);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask17);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 12));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask17);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 29));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 17 - 14);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask17);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 14));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask17);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 31));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);


    OutReg = _mm_srli_epi32(InReg, 17 - 16);
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    return offsets_buffer_i;
}
#endif



#ifdef HAVE_SSE2
/* nwritten = 9 * 4 = 36 */
static int
write_18_vert (FILE *offsetscomp_fp, Positionsptr_T *offsets_buffer, int offsets_buffer_size, int offsets_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask18);
    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask18);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 18));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 18 - 4);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask18);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 4));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask18);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 22));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 18 - 8);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask18);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask18);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 26));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 18 - 12);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask18);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 12));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask18);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 30));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 18 - 16);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask18);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 18 - 2);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask18);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 2));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask18);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 20));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 18 - 6);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask18);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 6));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask18);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 18 - 10);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask18);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 10));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask18);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 28));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 18 - 14);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask18);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 14));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    return offsets_buffer_i;
}
#endif

static int
write_18_horiz (UINT4 *out, const UINT4 *in) {
  int column;

  for (column = 0; column < 4; column++) {
    *out |= (*in)   % (1U << 18 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 18 )  ) <<  18 ;
    ++out;
    *out |=  ( (*in)   % (1U << 18 ) ) >> ( 18  -  4 );
    ++in;
    *out |= ( (*in)   % (1U << 18 )  ) <<  4 ;
    ++in;
    *out |= ( (*in)   % (1U << 18 )  ) <<  22 ;
    ++out;
    *out |=  ( (*in)   % (1U << 18 ) ) >> ( 18  -  8 );
    ++in;
    *out |= ( (*in)   % (1U << 18 )  ) <<  8 ;
    ++in;
    *out |= ( (*in)   % (1U << 18 )  ) <<  26 ;
    ++out;
    *out |=  ( (*in)   % (1U << 18 ) ) >> ( 18  -  12 );
    ++in;
    *out |= ( (*in)   % (1U << 18 )  ) <<  12 ;
    ++in;
    *out |= ( (*in)   % (1U << 18 )  ) <<  30 ;
    ++out;
    *out |=  ( (*in)   % (1U << 18 ) ) >> ( 18  -  16 );
    ++in;
    *out |= ( (*in)   % (1U << 18 )  ) <<  16 ;
    ++out;
    *out |=  ( (*in)   % (1U << 18 ) ) >> ( 18  -  2 );
    ++in;
    *out |= ( (*in)   % (1U << 18 )  ) <<  2 ;
    ++in;
    *out |= ( (*in)   % (1U << 18 )  ) <<  20 ;
    ++out;
    *out |=  ( (*in)   % (1U << 18 ) ) >> ( 18  -  6 );
    ++in;
    *out |= ( (*in)   % (1U << 18 )  ) <<  6 ;
    ++in;
    *out |= ( (*in)   % (1U << 18 )  ) <<  24 ;
    ++out;
    *out |=  ( (*in)   % (1U << 18 ) ) >> ( 18  -  10 );
    ++in;
    *out |= ( (*in)   % (1U << 18 )  ) <<  10 ;
    ++in;
    *out |= ( (*in)   % (1U << 18 )  ) <<  28 ;
    ++out;
    *out |=  ( (*in)   % (1U << 18 ) ) >> ( 18  -  14 );
    ++in;
    *out |= ( (*in)   % (1U << 18 )  ) <<  14 ;
    ++out;
    ++in;
  }

  return 36;
}



#ifdef ALLOW_ODD_PACKSIZES
/* nwritten = 10 * 4 = 40 */
static int
write_19_vert (FILE *offsetscomp_fp, Positionsptr_T *offsets_buffer, int offsets_buffer_size, int offsets_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask19);
    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask19);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 19));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 19 - 6);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask19);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 6));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask19);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 25));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 19 - 12);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask19);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 12));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask19);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 31));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 19 - 18);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask19);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 18));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 19 - 5);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask19);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 5));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask19);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 19 - 11);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask19);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 11));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask19);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 30));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 19 - 17);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask19);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 17));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 19 - 4);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask19);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 4));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask19);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 23));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 19 - 10);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask19);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 10));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask19);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 29));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);


    OutReg = _mm_srli_epi32(InReg, 19 - 16);
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    return offsets_buffer_i;
}
#endif


#ifdef HAVE_SSE2
/* nwritten = 10 * 4 = 40 */
static int
write_20_vert (FILE *offsetscomp_fp, Positionsptr_T *offsets_buffer, int offsets_buffer_size, int offsets_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask20);
    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask20);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 20));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 20 - 8);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask20);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask20);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 28));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 20 - 16);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask20);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 20 - 4);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask20);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 4));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask20);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 20 - 12);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask20);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 12));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    InReg = _mm_and_si128(_mm_load_si128(++in), mask20);

    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask20);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 20));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 20 - 8);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask20);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask20);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 28));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 20 - 16);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask20);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 20 - 4);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask20);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 4));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask20);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 20 - 12);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask20);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 12));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    return offsets_buffer_i;
}
#endif

static int
write_20_horiz (UINT4 *out, const UINT4 *in) {
  int column;

  for (column = 0; column < 4; column++) {
    *out |= (*in)   % (1U << 20 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 20 )  ) <<  20 ;
    ++out;
    *out |=  ( (*in)   % (1U << 20 ) ) >> ( 20  -  8 );
    ++in;
    *out |= ( (*in)   % (1U << 20 )  ) <<  8 ;
    ++in;
    *out |= ( (*in)   % (1U << 20 )  ) <<  28 ;
    ++out;
    *out |=  ( (*in)   % (1U << 20 ) ) >> ( 20  -  16 );
    ++in;
    *out |= ( (*in)   % (1U << 20 )  ) <<  16 ;
    ++out;
    *out |=  ( (*in)   % (1U << 20 ) ) >> ( 20  -  4 );
    ++in;
    *out |= ( (*in)   % (1U << 20 )  ) <<  4 ;
    ++in;
    *out |= ( (*in)   % (1U << 20 )  ) <<  24 ;
    ++out;
    *out |=  ( (*in)   % (1U << 20 ) ) >> ( 20  -  12 );
    ++in;
    *out |= ( (*in)   % (1U << 20 )  ) <<  12 ;
    ++out;
    ++in;
    *out |= (*in)   % (1U << 20 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 20 )  ) <<  20 ;
    ++out;
    *out |=  ( (*in)   % (1U << 20 ) ) >> ( 20  -  8 );
    ++in;
    *out |= ( (*in)   % (1U << 20 )  ) <<  8 ;
    ++in;
    *out |= ( (*in)   % (1U << 20 )  ) <<  28 ;
    ++out;
    *out |=  ( (*in)   % (1U << 20 ) ) >> ( 20  -  16 );
    ++in;
    *out |= ( (*in)   % (1U << 20 )  ) <<  16 ;
    ++out;
    *out |=  ( (*in)   % (1U << 20 ) ) >> ( 20  -  4 );
    ++in;
    *out |= ( (*in)   % (1U << 20 )  ) <<  4 ;
    ++in;
    *out |= ( (*in)   % (1U << 20 )  ) <<  24 ;
    ++out;
    *out |=  ( (*in)   % (1U << 20 ) ) >> ( 20  -  12 );
    ++in;
    *out |= ( (*in)   % (1U << 20 )  ) <<  12 ;
    ++out;
    ++in;
  }

  return 40;
}



#ifdef ALLOW_ODD_PACKSIZES
/* nwritten = 11 * 4 = 44 */
static int
write_21_vert (FILE *offsetscomp_fp, Positionsptr_T *offsets_buffer, int offsets_buffer_size, int offsets_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask21);
    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask21);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 21));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 21 - 10);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask21);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 10));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask21);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 31));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 21 - 20);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask21);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 20));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 21 - 9);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask21);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 9));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask21);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 30));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 21 - 19);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask21);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 19));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 21 - 8);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask21);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask21);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 29));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 21 - 18);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask21);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 18));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 21 - 7);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask21);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 7));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask21);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 28));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 21 - 17);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask21);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 17));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 21 - 6);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask21);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 6));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask21);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 27));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);


    OutReg = _mm_srli_epi32(InReg, 21 - 16);
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    return offsets_buffer_i;
}
#endif


#ifdef HAVE_SSE2
/* nwritten = 11 * 4 = 44 */
static int
write_22_vert (FILE *offsetscomp_fp, Positionsptr_T *offsets_buffer, int offsets_buffer_size, int offsets_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask22);
    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask22);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 22));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 22 - 12);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask22);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 12));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 22 - 2);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask22);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 2));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask22);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 22 - 14);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask22);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 14));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 22 - 4);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask22);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 4));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask22);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 26));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 22 - 16);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask22);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 22 - 6);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask22);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 6));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask22);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 28));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 22 - 18);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask22);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 18));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 22 - 8);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask22);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask22);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 30));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 22 - 20);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask22);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 20));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 22 - 10);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask22);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 10));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    return offsets_buffer_i;
}
#endif

static int
write_22_horiz (UINT4 *out, const UINT4 *in) {
  int column;

  for (column = 0; column < 4; column++) {
    *out |= (*in)   % (1U << 22 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 22 )  ) <<  22 ;
    ++out;
    *out |=  ( (*in)   % (1U << 22 ) ) >> ( 22  -  12 );
    ++in;
    *out |= ( (*in)   % (1U << 22 )  ) <<  12 ;
    ++out;
    *out |=  ( (*in)   % (1U << 22 ) ) >> ( 22  -  2 );
    ++in;
    *out |= ( (*in)   % (1U << 22 )  ) <<  2 ;
    ++in;
    *out |= ( (*in)   % (1U << 22 )  ) <<  24 ;
    ++out;
    *out |=  ( (*in)   % (1U << 22 ) ) >> ( 22  -  14 );
    ++in;
    *out |= ( (*in)   % (1U << 22 )  ) <<  14 ;
    ++out;
    *out |=  ( (*in)   % (1U << 22 ) ) >> ( 22  -  4 );
    ++in;
    *out |= ( (*in)   % (1U << 22 )  ) <<  4 ;
    ++in;
    *out |= ( (*in)   % (1U << 22 )  ) <<  26 ;
    ++out;
    *out |=  ( (*in)   % (1U << 22 ) ) >> ( 22  -  16 );
    ++in;
    *out |= ( (*in)   % (1U << 22 )  ) <<  16 ;
    ++out;
    *out |=  ( (*in)   % (1U << 22 ) ) >> ( 22  -  6 );
    ++in;
    *out |= ( (*in)   % (1U << 22 )  ) <<  6 ;
    ++in;
    *out |= ( (*in)   % (1U << 22 )  ) <<  28 ;
    ++out;
    *out |=  ( (*in)   % (1U << 22 ) ) >> ( 22  -  18 );
    ++in;
    *out |= ( (*in)   % (1U << 22 )  ) <<  18 ;
    ++out;
    *out |=  ( (*in)   % (1U << 22 ) ) >> ( 22  -  8 );
    ++in;
    *out |= ( (*in)   % (1U << 22 )  ) <<  8 ;
    ++in;
    *out |= ( (*in)   % (1U << 22 )  ) <<  30 ;
    ++out;
    *out |=  ( (*in)   % (1U << 22 ) ) >> ( 22  -  20 );
    ++in;
    *out |= ( (*in)   % (1U << 22 )  ) <<  20 ;
    ++out;
    *out |=  ( (*in)   % (1U << 22 ) ) >> ( 22  -  10 );
    ++in;
    *out |= ( (*in)   % (1U << 22 )  ) <<  10 ;
    ++out;
    ++in;
  }

  return 44;
}


#ifdef ALLOW_ODD_PACKSIZES
/* nwritten = 12 * 4 = 48 */
static int
write_23_vert (FILE *offsetscomp_fp, Positionsptr_T *offsets_buffer, int offsets_buffer_size, int offsets_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask23);
    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask23);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 23));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 23 - 14);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask23);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 14));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 23 - 5);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask23);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 5));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask23);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 28));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 23 - 19);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask23);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 19));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 23 - 10);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask23);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 10));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 23 - 1);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask23);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 1));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask23);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 23 - 15);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask23);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 15));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 23 - 6);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask23);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 6));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask23);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 29));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 23 - 20);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask23);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 20));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 23 - 11);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask23);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 11));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 23 - 2);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask23);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 2));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask23);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 25));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 23 - 16);
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    return offsets_buffer_i;
}
#endif


#ifdef HAVE_SSE2
/* nwritten = 12 * 4 = 48 */
static int
write_24_vert (FILE *offsetscomp_fp, Positionsptr_T *offsets_buffer, int offsets_buffer_size, int offsets_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask24);
    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask24);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 24 - 16);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask24);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 24 - 8);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask24);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    InReg = _mm_and_si128(_mm_load_si128(++in), mask24);

    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask24);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 24 - 16);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask24);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 24 - 8);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask24);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    InReg = _mm_and_si128(_mm_load_si128(++in), mask24);

    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask24);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 24 - 16);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask24);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 24 - 8);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask24);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    InReg = _mm_and_si128(_mm_load_si128(++in), mask24);

    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask24);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 24 - 16);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask24);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 24 - 8);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask24);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    return offsets_buffer_i;
}
#endif

static int
write_24_horiz (UINT4 *out, const UINT4 *in) {
  int column;

  for (column = 0; column < 4; column++) {
    *out |= (*in)   % (1U << 24 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 24 )  ) <<  24 ;
    ++out;
    *out |=  ( (*in)   % (1U << 24 ) ) >> ( 24  -  16 );
    ++in;
    *out |= ( (*in)   % (1U << 24 )  ) <<  16 ;
    ++out;
    *out |=  ( (*in)   % (1U << 24 ) ) >> ( 24  -  8 );
    ++in;
    *out |= ( (*in)   % (1U << 24 )  ) <<  8 ;
    ++out;
    ++in;
    *out |= (*in)   % (1U << 24 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 24 )  ) <<  24 ;
    ++out;
    *out |=  ( (*in)   % (1U << 24 ) ) >> ( 24  -  16 );
    ++in;
    *out |= ( (*in)   % (1U << 24 )  ) <<  16 ;
    ++out;
    *out |=  ( (*in)   % (1U << 24 ) ) >> ( 24  -  8 );
    ++in;
    *out |= ( (*in)   % (1U << 24 )  ) <<  8 ;
    ++out;
    ++in;
    *out |= (*in)   % (1U << 24 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 24 )  ) <<  24 ;
    ++out;
    *out |=  ( (*in)   % (1U << 24 ) ) >> ( 24  -  16 );
    ++in;
    *out |= ( (*in)   % (1U << 24 )  ) <<  16 ;
    ++out;
    *out |=  ( (*in)   % (1U << 24 ) ) >> ( 24  -  8 );
    ++in;
    *out |= ( (*in)   % (1U << 24 )  ) <<  8 ;
    ++out;
    ++in;
    *out |= (*in)   % (1U << 24 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 24 )  ) <<  24 ;
    ++out;
    *out |=  ( (*in)   % (1U << 24 ) ) >> ( 24  -  16 );
    ++in;
    *out |= ( (*in)   % (1U << 24 )  ) <<  16 ;
    ++out;
    *out |=  ( (*in)   % (1U << 24 ) ) >> ( 24  -  8 );
    ++in;
    *out |= ( (*in)   % (1U << 24 )  ) <<  8 ;
    ++out;
    ++in;
  }

  return 48;
}


#ifdef ALLOW_ODD_PACKSIZES
/* nwritten = 13 * 4 = 52 */
static int
write_25_vert (FILE *offsetscomp_fp, Positionsptr_T *offsets_buffer, int offsets_buffer_size, int offsets_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask25);
    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask25);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 25));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 25 - 18);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask25);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 18));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 25 - 11);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask25);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 11));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 25 - 4);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask25);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 4));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask25);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 29));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 25 - 22);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask25);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 22));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 25 - 15);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask25);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 15));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 25 - 8);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask25);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 25 - 1);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask25);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 1));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask25);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 26));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 25 - 19);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask25);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 19));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 25 - 12);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask25);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 12));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 25 - 5);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask25);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 5));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask25);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 30));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 25 - 23);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask25);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 23));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);


    OutReg = _mm_srli_epi32(InReg, 25 - 16);
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    return offsets_buffer_i;
}
#endif


#ifdef HAVE_SSE2
/* nwritten = 13 * 4 = 52 */
static int
write_26_vert (FILE *offsetscomp_fp, Positionsptr_T *offsets_buffer, int offsets_buffer_size, int offsets_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask26);
    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask26);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 26));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 26 - 20);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask26);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 20));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 26 - 14);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask26);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 14));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 26 - 8);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask26);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 26 - 2);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask26);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 2));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask26);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 28));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 26 - 22);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask26);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 22));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 26 - 16);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask26);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 26 - 10);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask26);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 10));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 26 - 4);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask26);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 4));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask26);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 30));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 26 - 24);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask26);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 26 - 18);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask26);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 18));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 26 - 12);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask26);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 12));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 26 - 6);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask26);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 6));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    return offsets_buffer_i;
}
#endif


static int
write_26_horiz (UINT4 *out, const UINT4 *in) {
  int column;

  for (column = 0; column < 4; column++) {
    *out |= (*in)   % (1U << 26 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 26 )  ) <<  26 ;
    ++out;
    *out |=  ( (*in)   % (1U << 26 ) ) >> ( 26  -  20 );
    ++in;
    *out |= ( (*in)   % (1U << 26 )  ) <<  20 ;
    ++out;
    *out |=  ( (*in)   % (1U << 26 ) ) >> ( 26  -  14 );
    ++in;
    *out |= ( (*in)   % (1U << 26 )  ) <<  14 ;
    ++out;
    *out |=  ( (*in)   % (1U << 26 ) ) >> ( 26  -  8 );
    ++in;
    *out |= ( (*in)   % (1U << 26 )  ) <<  8 ;
    ++out;
    *out |=  ( (*in)   % (1U << 26 ) ) >> ( 26  -  2 );
    ++in;
    *out |= ( (*in)   % (1U << 26 )  ) <<  2 ;
    ++in;
    *out |= ( (*in)   % (1U << 26 )  ) <<  28 ;
    ++out;
    *out |=  ( (*in)   % (1U << 26 ) ) >> ( 26  -  22 );
    ++in;
    *out |= ( (*in)   % (1U << 26 )  ) <<  22 ;
    ++out;
    *out |=  ( (*in)   % (1U << 26 ) ) >> ( 26  -  16 );
    ++in;
    *out |= ( (*in)   % (1U << 26 )  ) <<  16 ;
    ++out;
    *out |=  ( (*in)   % (1U << 26 ) ) >> ( 26  -  10 );
    ++in;
    *out |= ( (*in)   % (1U << 26 )  ) <<  10 ;
    ++out;
    *out |=  ( (*in)   % (1U << 26 ) ) >> ( 26  -  4 );
    ++in;
    *out |= ( (*in)   % (1U << 26 )  ) <<  4 ;
    ++in;
    *out |= ( (*in)   % (1U << 26 )  ) <<  30 ;
    ++out;
    *out |=  ( (*in)   % (1U << 26 ) ) >> ( 26  -  24 );
    ++in;
    *out |= ( (*in)   % (1U << 26 )  ) <<  24 ;
    ++out;
    *out |=  ( (*in)   % (1U << 26 ) ) >> ( 26  -  18 );
    ++in;
    *out |= ( (*in)   % (1U << 26 )  ) <<  18 ;
    ++out;
    *out |=  ( (*in)   % (1U << 26 ) ) >> ( 26  -  12 );
    ++in;
    *out |= ( (*in)   % (1U << 26 )  ) <<  12 ;
    ++out;
    *out |=  ( (*in)   % (1U << 26 ) ) >> ( 26  -  6 );
    ++in;
    *out |= ( (*in)   % (1U << 26 )  ) <<  6 ;
    ++out;
    ++in;
  }

  return 52;
}



#ifdef ALLOW_ODD_PACKSIZES
/* nwritten = 14 * 4 = 56 */
static int
write_27_vert (FILE *offsetscomp_fp, Positionsptr_T *offsets_buffer, int offsets_buffer_size, int offsets_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask27);
    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask27);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 27));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 27 - 22);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask27);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 22));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 27 - 17);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask27);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 17));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 27 - 12);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask27);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 12));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 27 - 7);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask27);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 7));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 27 - 2);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask27);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 2));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask27);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 29));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 27 - 24);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask27);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 27 - 19);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask27);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 19));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 27 - 14);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask27);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 14));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 27 - 9);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask27);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 9));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 27 - 4);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask27);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 4));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask27);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 31));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 27 - 26);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask27);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 26));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 27 - 21);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask27);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 21));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 27 - 16);
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    return offsets_buffer_i;
}
#endif


#ifdef HAVE_SSE2
/* nwritten = 14 * 4 = 56 */
static int
write_28_vert (FILE *offsetscomp_fp, Positionsptr_T *offsets_buffer, int offsets_buffer_size, int offsets_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask28);
    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask28);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 28));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 28 - 24);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask28);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 28 - 20);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask28);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 20));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 28 - 16);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask28);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 28 - 12);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask28);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 12));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 28 - 8);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask28);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 28 - 4);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask28);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 4));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    InReg = _mm_and_si128(_mm_load_si128(++in), mask28);

    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask28);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 28));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 28 - 24);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask28);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 28 - 20);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask28);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 20));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 28 - 16);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask28);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 28 - 12);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask28);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 12));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 28 - 8);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask28);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 28 - 4);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask28);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 4));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    return offsets_buffer_i;
}
#endif


static int
write_28_horiz (UINT4 *out, const UINT4 *in) {
  int column;

  for (column = 0; column < 4; column++) {
    *out |= (*in)   % (1U << 28 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 28 )  ) <<  28 ;
    ++out;
    *out |=  ( (*in)   % (1U << 28 ) ) >> ( 28  -  24 );
    ++in;
    *out |= ( (*in)   % (1U << 28 )  ) <<  24 ;
    ++out;
    *out |=  ( (*in)   % (1U << 28 ) ) >> ( 28  -  20 );
    ++in;
    *out |= ( (*in)   % (1U << 28 )  ) <<  20 ;
    ++out;
    *out |=  ( (*in)   % (1U << 28 ) ) >> ( 28  -  16 );
    ++in;
    *out |= ( (*in)   % (1U << 28 )  ) <<  16 ;
    ++out;
    *out |=  ( (*in)   % (1U << 28 ) ) >> ( 28  -  12 );
    ++in;
    *out |= ( (*in)   % (1U << 28 )  ) <<  12 ;
    ++out;
    *out |=  ( (*in)   % (1U << 28 ) ) >> ( 28  -  8 );
    ++in;
    *out |= ( (*in)   % (1U << 28 )  ) <<  8 ;
    ++out;
    *out |=  ( (*in)   % (1U << 28 ) ) >> ( 28  -  4 );
    ++in;
    *out |= ( (*in)   % (1U << 28 )  ) <<  4 ;
    ++out;
    ++in;
    *out |= (*in)   % (1U << 28 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 28 )  ) <<  28 ;
    ++out;
    *out |=  ( (*in)   % (1U << 28 ) ) >> ( 28  -  24 );
    ++in;
    *out |= ( (*in)   % (1U << 28 )  ) <<  24 ;
    ++out;
    *out |=  ( (*in)   % (1U << 28 ) ) >> ( 28  -  20 );
    ++in;
    *out |= ( (*in)   % (1U << 28 )  ) <<  20 ;
    ++out;
    *out |=  ( (*in)   % (1U << 28 ) ) >> ( 28  -  16 );
    ++in;
    *out |= ( (*in)   % (1U << 28 )  ) <<  16 ;
    ++out;
    *out |=  ( (*in)   % (1U << 28 ) ) >> ( 28  -  12 );
    ++in;
    *out |= ( (*in)   % (1U << 28 )  ) <<  12 ;
    ++out;
    *out |=  ( (*in)   % (1U << 28 ) ) >> ( 28  -  8 );
    ++in;
    *out |= ( (*in)   % (1U << 28 )  ) <<  8 ;
    ++out;
    *out |=  ( (*in)   % (1U << 28 ) ) >> ( 28  -  4 );
    ++in;
    *out |= ( (*in)   % (1U << 28 )  ) <<  4 ;
    ++out;
    ++in;
  }

  return 56;
}



#ifdef ALLOW_ODD_PACKSIZES
/* nwritten = 15 * 4 = 60 */
static int
write_29_vert (FILE *offsetscomp_fp, Positionsptr_T *offsets_buffer, int offsets_buffer_size, int offsets_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask29);
    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask29);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 29));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 29 - 26);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask29);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 26));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 29 - 23);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask29);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 23));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 29 - 20);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask29);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 20));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 29 - 17);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask29);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 17));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 29 - 14);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask29);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 14));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 29 - 11);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask29);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 11));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 29 - 8);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask29);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 29 - 5);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask29);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 5));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 29 - 2);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask29);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 2));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask29);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 31));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 29 - 28);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask29);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 28));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 29 - 25);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask29);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 25));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 29 - 22);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask29);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 22));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 29 - 19);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask29);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 19));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);


    OutReg = _mm_srli_epi32(InReg, 29 - 16);
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    return offsets_buffer_i;
}
#endif


#ifdef HAVE_SSE2
/* nwritten = 15 * 4 = 60 */
static int
write_30_vert (FILE *offsetscomp_fp, Positionsptr_T *offsets_buffer, int offsets_buffer_size, int offsets_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask30);
    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask30);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 30));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 30 - 28);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask30);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 28));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 30 - 26);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask30);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 26));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 30 - 24);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask30);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 30 - 22);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask30);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 22));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 30 - 20);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask30);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 20));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 30 - 18);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask30);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 18));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 30 - 16);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask30);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 30 - 14);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask30);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 14));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 30 - 12);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask30);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 12));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 30 - 10);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask30);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 10));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 30 - 8);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask30);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 30 - 6);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask30);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 6));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 30 - 4);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask30);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 4));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 30 - 2);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask30);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 2));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    return offsets_buffer_i;
}
#endif


static int
write_30_horiz (UINT4 *out, const UINT4 *in) {
  int column;

  for (column = 0; column < 4; column++) {
    *out |= (*in)   % (1U << 30 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 30 )  ) <<  30 ;
    ++out;
    *out |=  ( (*in)   % (1U << 30 ) ) >> ( 30  -  28 );
    ++in;
    *out |= ( (*in)   % (1U << 30 )  ) <<  28 ;
    ++out;
    *out |=  ( (*in)   % (1U << 30 ) ) >> ( 30  -  26 );
    ++in;
    *out |= ( (*in)   % (1U << 30 )  ) <<  26 ;
    ++out;
    *out |=  ( (*in)   % (1U << 30 ) ) >> ( 30  -  24 );
    ++in;
    *out |= ( (*in)   % (1U << 30 )  ) <<  24 ;
    ++out;
    *out |=  ( (*in)   % (1U << 30 ) ) >> ( 30  -  22 );
    ++in;
    *out |= ( (*in)   % (1U << 30 )  ) <<  22 ;
    ++out;
    *out |=  ( (*in)   % (1U << 30 ) ) >> ( 30  -  20 );
    ++in;
    *out |= ( (*in)   % (1U << 30 )  ) <<  20 ;
    ++out;
    *out |=  ( (*in)   % (1U << 30 ) ) >> ( 30  -  18 );
    ++in;
    *out |= ( (*in)   % (1U << 30 )  ) <<  18 ;
    ++out;
    *out |=  ( (*in)   % (1U << 30 ) ) >> ( 30  -  16 );
    ++in;
    *out |= ( (*in)   % (1U << 30 )  ) <<  16 ;
    ++out;
    *out |=  ( (*in)   % (1U << 30 ) ) >> ( 30  -  14 );
    ++in;
    *out |= ( (*in)   % (1U << 30 )  ) <<  14 ;
    ++out;
    *out |=  ( (*in)   % (1U << 30 ) ) >> ( 30  -  12 );
    ++in;
    *out |= ( (*in)   % (1U << 30 )  ) <<  12 ;
    ++out;
    *out |=  ( (*in)   % (1U << 30 ) ) >> ( 30  -  10 );
    ++in;
    *out |= ( (*in)   % (1U << 30 )  ) <<  10 ;
    ++out;
    *out |=  ( (*in)   % (1U << 30 ) ) >> ( 30  -  8 );
    ++in;
    *out |= ( (*in)   % (1U << 30 )  ) <<  8 ;
    ++out;
    *out |=  ( (*in)   % (1U << 30 ) ) >> ( 30  -  6 );
    ++in;
    *out |= ( (*in)   % (1U << 30 )  ) <<  6 ;
    ++out;
    *out |=  ( (*in)   % (1U << 30 ) ) >> ( 30  -  4 );
    ++in;
    *out |= ( (*in)   % (1U << 30 )  ) <<  4 ;
    ++out;
    *out |=  ( (*in)   % (1U << 30 ) ) >> ( 30  -  2 );
    ++in;
    *out |= ( (*in)   % (1U << 30 )  ) <<  2 ;
    ++out;
    ++in;
  }

  return 60;
}



#ifdef ALLOW_ODD_PACKSIZES
/* nwritten = 16 * 4 = 64 */
static int
write_31_vert (FILE *offsetscomp_fp, Positionsptr_T *offsets_buffer, int offsets_buffer_size, int offsets_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask31);
    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask31);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 31));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 31 - 30);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask31);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 30));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 31 - 29);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask31);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 29));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 31 - 28);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask31);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 28));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 31 - 27);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask31);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 27));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 31 - 26);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask31);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 26));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 31 - 25);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask31);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 25));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 31 - 24);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask31);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 31 - 23);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask31);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 23));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 31 - 22);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask31);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 22));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 31 - 21);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask31);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 21));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 31 - 20);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask31);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 20));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 31 - 19);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask31);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 19));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 31 - 18);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask31);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 18));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    OutReg = _mm_srli_epi32(InReg, 31 - 17);
    InReg = _mm_and_si128(_mm_load_si128(++in), mask31);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 17));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);


    OutReg = _mm_srli_epi32(InReg, 31 - 16);
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    return offsets_buffer_i;
}
#endif


#ifdef HAVE_SSE2
/* nwritten = 16 * 4 = 64 */
static int
write_32_vert (FILE *offsetscomp_fp, Positionsptr_T *offsets_buffer, int offsets_buffer_size, int offsets_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_load_si128(in);
    OutReg = InReg;
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    InReg = _mm_load_si128(++in);

    OutReg = InReg;
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    InReg = _mm_load_si128(++in);

    OutReg = InReg;
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    InReg = _mm_load_si128(++in);

    OutReg = InReg;
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    InReg = _mm_load_si128(++in);

    OutReg = InReg;
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    InReg = _mm_load_si128(++in);

    OutReg = InReg;
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    InReg = _mm_load_si128(++in);

    OutReg = InReg;
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    InReg = _mm_load_si128(++in);

    OutReg = InReg;
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    InReg = _mm_load_si128(++in);

    OutReg = InReg;
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    InReg = _mm_load_si128(++in);

    OutReg = InReg;
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    InReg = _mm_load_si128(++in);

    OutReg = InReg;
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    InReg = _mm_load_si128(++in);

    OutReg = InReg;
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    InReg = _mm_load_si128(++in);

    OutReg = InReg;
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    InReg = _mm_load_si128(++in);

    OutReg = InReg;
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    InReg = _mm_load_si128(++in);

    OutReg = InReg;
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    InReg = _mm_load_si128(++in);

    OutReg = InReg;
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    return offsets_buffer_i;
}
#endif


static int
write_32_horiz (UINT4 *out, const UINT4 *in) {
  int column;

  for (column = 0; column < 4; column++) {
    *out = *in;
    ++out;
    ++in;
    *out = *in;
    ++out;
    ++in;
    *out = *in;
    ++out;
    ++in;
    *out = *in;
    ++out;
    ++in;
    *out = *in;
    ++out;
    ++in;
    *out = *in;
    ++out;
    ++in;
    *out = *in;
    ++out;
    ++in;
    *out = *in;
    ++out;
    ++in;
    *out = *in;
    ++out;
    ++in;
    *out = *in;
    ++out;
    ++in;
    *out = *in;
    ++out;
    ++in;
    *out = *in;
    ++out;
    ++in;
    *out = *in;
    ++out;
    ++in;
    *out = *in;
    ++out;
    ++in;
    *out = *in;
    ++out;
    ++in;
    *out = *in;
    ++out;
    ++in;
  }

  return 64;
}



#ifdef HAVE_SSE2
/* nwritten = 2 * 4 = 8 */
static int
write_04_vert (FILE *offsetscomp_fp, Positionsptr_T *offsets_buffer, int offsets_buffer_size, int offsets_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask4);
    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask4);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 4));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask4);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask4);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 12));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask4);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask4);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 20));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask4);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask4);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 28));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    InReg = _mm_and_si128(_mm_load_si128(++in), mask4);

    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask4);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 4));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask4);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask4);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 12));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask4);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask4);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 20));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask4);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask4);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 28));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    return offsets_buffer_i;
}
#endif


static int
write_04_horiz (UINT4 *out, const UINT4 *in) {
  int column;

  for (column = 0; column < 4; column++) {
    *out |= (*in)   % (1U << 4 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 4 )  ) <<  4 ;
    ++in;
    *out |= ( (*in)   % (1U << 4 )  ) <<  8 ;
    ++in;
    *out |= ( (*in)   % (1U << 4 )  ) <<  12 ;
    ++in;
    *out |= ( (*in)   % (1U << 4 )  ) <<  16 ;
    ++in;
    *out |= ( (*in)   % (1U << 4 )  ) <<  20 ;
    ++in;
    *out |= ( (*in)   % (1U << 4 )  ) <<  24 ;
    ++in;
    *out |= ( (*in)   % (1U << 4 )  ) <<  28 ;
    ++out;
    ++in;
    *out |= (*in)   % (1U << 4 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 4 )  ) <<  4 ;
    ++in;
    *out |= ( (*in)   % (1U << 4 )  ) <<  8 ;
    ++in;
    *out |= ( (*in)   % (1U << 4 )  ) <<  12 ;
    ++in;
    *out |= ( (*in)   % (1U << 4 )  ) <<  16 ;
    ++in;
    *out |= ( (*in)   % (1U << 4 )  ) <<  20 ;
    ++in;
    *out |= ( (*in)   % (1U << 4 )  ) <<  24 ;
    ++in;
    *out |= ( (*in)   % (1U << 4 )  ) <<  28 ;
    ++out;
    ++in;
  }

  return 8;
}


#ifdef HAVE_SSE2
/* nwritten = 4 * 4 = 16 */
static int
write_08_vert (FILE *offsetscomp_fp, Positionsptr_T *offsets_buffer, int offsets_buffer_size, int offsets_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask8);
    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask8);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask8);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask8);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    InReg = _mm_and_si128(_mm_load_si128(++in), mask8);

    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask8);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask8);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask8);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    InReg = _mm_and_si128(_mm_load_si128(++in), mask8);

    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask8);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask8);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask8);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    InReg = _mm_and_si128(_mm_load_si128(++in), mask8);

    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask8);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 8));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask8);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    InReg = _mm_and_si128(_mm_load_si128(++in), mask8);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 24));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    return offsets_buffer_i;
}
#endif


static int
write_08_horiz (UINT4 *out, const UINT4 *in) {
  int column;

  for (column = 0; column < 4; column++) {
    *out |= (*in)   % (1U << 8 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 8 )  ) <<  8 ;
    ++in;
    *out |= ( (*in)   % (1U << 8 )  ) <<  16 ;
    ++in;
    *out |= ( (*in)   % (1U << 8 )  ) <<  24 ;
    ++out;
    ++in;
    *out |= (*in)   % (1U << 8 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 8 )  ) <<  8 ;
    ++in;
    *out |= ( (*in)   % (1U << 8 )  ) <<  16 ;
    ++in;
    *out |= ( (*in)   % (1U << 8 )  ) <<  24 ;
    ++out;
    ++in;
    *out |= (*in)   % (1U << 8 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 8 )  ) <<  8 ;
    ++in;
    *out |= ( (*in)   % (1U << 8 )  ) <<  16 ;
    ++in;
    *out |= ( (*in)   % (1U << 8 )  ) <<  24 ;
    ++out;
    ++in;
    *out |= (*in)   % (1U << 8 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 8 )  ) <<  8 ;
    ++in;
    *out |= ( (*in)   % (1U << 8 )  ) <<  16 ;
    ++in;
    *out |= ( (*in)   % (1U << 8 )  ) <<  24 ;
    ++out;
    ++in;
  }

  return 16;
}


#ifdef HAVE_SSE2
/* nwritten = 8 * 4 = 32 */
static int
write_16_vert (FILE *offsetscomp_fp, Positionsptr_T *offsets_buffer, int offsets_buffer_size, int offsets_buffer_i, const UINT4 *_in) {
    const __m128i *in = (const __m128i *) _in;
    __m128i OutReg;

    __m128i InReg = _mm_and_si128(_mm_load_si128(in), mask16);

    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask16);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    InReg = _mm_and_si128(_mm_load_si128(++in), mask16);

    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask16);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    InReg = _mm_and_si128(_mm_load_si128(++in), mask16);

    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask16);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    InReg = _mm_and_si128(_mm_load_si128(++in), mask16);

    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask16);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    InReg = _mm_and_si128(_mm_load_si128(++in), mask16);

    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask16);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    InReg = _mm_and_si128(_mm_load_si128(++in), mask16);

    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask16);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    InReg = _mm_and_si128(_mm_load_si128(++in), mask16);

    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask16);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    InReg = _mm_and_si128(_mm_load_si128(++in), mask16);

    OutReg = InReg;
    InReg = _mm_and_si128(_mm_load_si128(++in), mask16);

    OutReg =  _mm_or_si128(OutReg,_mm_slli_epi32(InReg, 16));
    offsets_buffer_i = write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,
					       OutReg);

    return offsets_buffer_i;
}
#endif


static int
write_16_horiz (UINT4 *out, const UINT4 *in) {
  int column;

  for (column = 0; column < 4; column++) {
    *out |= (*in)   % (1U << 16 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 16 )  ) <<  16 ;
    ++out;
    ++in;
    *out |= (*in)   % (1U << 16 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 16 )  ) <<  16 ;
    ++out;
    ++in;
    *out |= (*in)   % (1U << 16 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 16 )  ) <<  16 ;
    ++out;
    ++in;
    *out |= (*in)   % (1U << 16 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 16 )  ) <<  16 ;
    ++out;
    ++in;
    *out |= (*in)   % (1U << 16 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 16 )  ) <<  16 ;
    ++out;
    ++in;
    *out |= (*in)   % (1U << 16 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 16 )  ) <<  16 ;
    ++out;
    ++in;
    *out |= (*in)   % (1U << 16 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 16 )  ) <<  16 ;
    ++out;
    ++in;
    *out |= (*in)   % (1U << 16 ) ;
    ++in;
    *out |= ( (*in)   % (1U << 16 )  ) <<  16 ;
    ++out;
    ++in;
  }

  return 32;
}


/* Vertical format is better for cumulative sums, which require all
   values in a block to be decoded */
#ifdef HAVE_SSE2
int
Bitpack64_write_vert (FILE *offsetscomp_fp, Positionsptr_T *offsets_buffer, int offsets_buffer_size, int offsets_buffer_i,
		      const UINT4 *_in, int packsize) {

#if 0
  int i;

  printf("Entering with packsize %d\n",packsize);
  for (i = 0; i < 64; i++) {
    printf("%d ",_in[i]);
  }
  printf("\n");
#endif

  switch (packsize) {
#ifdef ALLOW_ODD_PACKSIZES
  case 0: return offsets_buffer_i;
  case 1: return write_01_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,_in);
  case 2: return write_02_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,_in);
  case 3: return write_03_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,_in);
  case 4: return write_04_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,_in);
  case 5: return write_05_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,_in);
  case 6: return write_06_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,_in);
  case 7: return write_07_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,_in);
  case 8: return write_08_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,_in);
  case 9: return write_09_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,_in);
  case 10: return write_10_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,_in);
  case 11: return write_11_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,_in);
  case 12: return write_12_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,_in);
  case 13: return write_13_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,_in);
  case 14: return write_14_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,_in);
  case 15: return write_15_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,_in);
  case 16: return write_16_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,_in);
  case 17: return write_17_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,_in);
  case 18: return write_18_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,_in);
  case 19: return write_19_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,_in);
  case 20: return write_20_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,_in);
  case 21: return write_21_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,_in);
  case 22: return write_22_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,_in);
  case 23: return write_23_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,_in);
  case 24: return write_24_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,_in);
  case 25: return write_25_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,_in);
  case 26: return write_26_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,_in);
  case 27: return write_27_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,_in);
  case 28: return write_28_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,_in);
  case 29: return write_29_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,_in);
  case 30: return write_30_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,_in);
  case 31: return write_31_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,_in);
  case 32: return write_32_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,_in);
#else
  case 0: return offsets_buffer_i;
  case 2: return write_02_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,_in);
  case 4: return write_04_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,_in);
  case 6: return write_06_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,_in);
  case 8: return write_08_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,_in);
  case 10: return write_10_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,_in);
  case 12: return write_12_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,_in);
  case 14: return write_14_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,_in);
  case 16: return write_16_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,_in);
  case 18: return write_18_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,_in);
  case 20: return write_20_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,_in);
  case 22: return write_22_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,_in);
  case 24: return write_24_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,_in);
  case 26: return write_26_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,_in);
  case 28: return write_28_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,_in);
  case 30: return write_30_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,_in);
  case 32: return write_32_vert(offsetscomp_fp,offsets_buffer,offsets_buffer_size,offsets_buffer_i,_in);
#endif
  default: fprintf(stderr,"packsize of %d not allowed\n",packsize); abort();
  }
}

#else

static void
reorder_values_vertically (Positionsptr_T *vertical, Positionsptr_T *horizontal) {
  int column, row, k = 0;
  Positionsptr_T *out;

  out = &(vertical[0]);
  for (column = 0; column < 4; column++) {
    k = column;
    for (row = 0; row < 16; row++) {
      *out++ = horizontal[k];
      k += 4;
    }
  }

#if 0
  printf("horizontal\n");
  for (k = 0; k < 64; k++) {
    if (k % 4 == 0) {
      printf("\n");
    }
    printf("%u ",horizontal[k]);
  }
  printf("\n");

  printf("vertical\n");
  for (k = 0; k < 64; k++) {
    if (k % 16 == 0) {
      printf("\n");
    }
    printf("%u ",vertical[k]);
  }
  printf("\n");
#endif

  return;
}


/* Non-SIMD code cannot write vertical format easily, so using
   horizontal code and conversions */
int
Bitpack64_write_vert (FILE *offsetscomp_fp, Positionsptr_T *offsets_buffer, int offsets_buffer_size, int offsets_buffer_i,
		      const UINT4 *horizontal, int packsize) {
  int nwritten;
  UINT4 buffer[64], vertical[64];

#if 0
  int i;

  printf("Entering with packsize %d\n",packsize);
  for (i = 0; i < 64; i++) {
    printf("%d ",_in[i]);
  }
  printf("\n");
#endif

  reorder_values_vertically(vertical,horizontal);
  memset((void *) buffer,0,64*sizeof(UINT4));

  switch (packsize) {
  case 0: return offsets_buffer_i;
  case 2: nwritten = write_02_horiz(buffer,&(vertical[0])); break;
  case 4: nwritten = write_04_horiz(buffer,&(vertical[0])); break;
  case 6: nwritten = write_06_horiz(buffer,&(vertical[0])); break;
  case 8: nwritten = write_08_horiz(buffer,&(vertical[0])); break;
  case 10: nwritten = write_10_horiz(buffer,&(vertical[0])); break;
  case 12: nwritten = write_12_horiz(buffer,&(vertical[0])); break;
  case 14: nwritten = write_14_horiz(buffer,&(vertical[0])); break;
  case 16: nwritten = write_16_horiz(buffer,&(vertical[0])); break;
  case 18: nwritten = write_18_horiz(buffer,&(vertical[0])); break;
  case 20: nwritten = write_20_horiz(buffer,&(vertical[0])); break;
  case 22: nwritten = write_22_horiz(buffer,&(vertical[0])); break;
  case 24: nwritten = write_24_horiz(buffer,&(vertical[0])); break;
  case 26: nwritten = write_26_horiz(buffer,&(vertical[0])); break;
  case 28: nwritten = write_28_horiz(buffer,&(vertical[0])); break;
  case 30: nwritten = write_30_horiz(buffer,&(vertical[0])); break;
  case 32: nwritten = write_32_horiz(buffer,&(vertical[0])); break;
  default: fprintf(stderr,"packsize of %d not allowed\n",packsize); abort();
  }

  return write_reg_buffered_vert(offsetscomp_fp,offsets_buffer,
				 offsets_buffer_size,offsets_buffer_i,
				 buffer,nwritten);
}
#endif


/* Horizontal format is slightly better for random access of individual values */
int
Bitpack64_write_horiz (FILE *offsetscomp_fp, Positionsptr_T *offsets_buffer, int offsets_buffer_size, int offsets_buffer_i,
		       const UINT4 *horizontal, int packsize) {
  int nwritten;
  UINT4 buffer[64];

#if 0
  int i;

  printf("Entering with packsize %d\n",packsize);
  for (i = 0; i < 64; i++) {
    printf("%d ",_in[i]);
  }
  printf("\n");
#endif

  memset((void *) buffer,0,64*sizeof(UINT4));

  switch (packsize) {
  case 0: return offsets_buffer_i;
  case 2: nwritten = write_02_horiz(buffer,&(horizontal[0])); break;
  case 4: nwritten = write_04_horiz(buffer,&(horizontal[0])); break;
  case 6: nwritten = write_06_horiz(buffer,&(horizontal[0])); break;
  case 8: nwritten = write_08_horiz(buffer,&(horizontal[0])); break;
  case 10: nwritten = write_10_horiz(buffer,&(horizontal[0])); break;
  case 12: nwritten = write_12_horiz(buffer,&(horizontal[0])); break;
  case 14: nwritten = write_14_horiz(buffer,&(horizontal[0])); break;
  case 16: nwritten = write_16_horiz(buffer,&(horizontal[0])); break;
  case 18: nwritten = write_18_horiz(buffer,&(horizontal[0])); break;
  case 20: nwritten = write_20_horiz(buffer,&(horizontal[0])); break;
  case 22: nwritten = write_22_horiz(buffer,&(horizontal[0])); break;
  case 24: nwritten = write_24_horiz(buffer,&(horizontal[0])); break;
  case 26: nwritten = write_26_horiz(buffer,&(horizontal[0])); break;
  case 28: nwritten = write_28_horiz(buffer,&(horizontal[0])); break;
  case 30: nwritten = write_30_horiz(buffer,&(horizontal[0])); break;
  case 32: nwritten = write_32_horiz(buffer,&(horizontal[0])); break;
  default: fprintf(stderr,"packsize of %d not allowed\n",packsize); abort();
  }

  return write_reg_buffered_horiz(offsetscomp_fp,offsets_buffer,
				  offsets_buffer_size,offsets_buffer_i,
				  buffer,nwritten);
}

