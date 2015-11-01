static char rcsid[] = "$Id: bitpack64-read.c 109826 2013-10-02 22:32:35Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "bitpack64-read.h"

#include <stdio.h>
#include <stdlib.h>

#ifdef HAVE_SSE2
#include <emmintrin.h>
#endif


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


/* #define ALLOW_ODD_PACKSIZES 1 */

#ifdef HAVE_SSE2
#ifdef DEBUG
/* For debugging */
static void
print_vector_hex (__m128i x) {
  UINT4 *s = (UINT4 *) &x;

  printf("%08X %08X %08X %08X\n",s[0],s[1],s[2],s[3]);
  return;
}

static void
print_vector (__m128i x) {
  UINT4 *s = (UINT4 *) &x;

  printf("%u %u %u %u\n",s[0],s[1],s[2],s[3]);
  return;
}
#endif
#endif


#ifdef HAVE_SSE2
#ifdef ALLOW_ODD_PACKSIZES
static __m128i mask1, mask2, mask3, mask4, mask5, mask6, mask7, mask8,
  mask9, mask10, mask11, mask12, mask13, mask14, mask15, mask16,
  mask17, mask18, mask19, mask20, mask21, mask22, mask23, mask24,
  mask25, mask26, mask27, mask28, mask29, mask30, mask31;
#else
static __m128i mask2, mask4, mask6, mask8, mask10, mask12, mask14, mask16,
  mask18, mask20, mask22, mask24, mask26, mask28, mask30;
#endif
#endif


static UINT4 *bitpackptrs;
static UINT4 *offsetscomp;
/* static Blocksize_T offsetscomp_blocksize = 64; */
#define OFFSETSCOMP_BLOCKSIZE 64

void
Bitpack64_read_setup (UINT4 *bitpackptrs_in, UINT4 *offsetscomp_in,
		      Blocksize_T offsetscomp_blocksize_in) {
  bitpackptrs = bitpackptrs_in;
  offsetscomp = offsetscomp_in;
  /* offsetscomp_blocksize = offsetscomp_blocksize_in; */

#ifdef HAVE_SSE2
#ifdef ALLOW_ODD_PACKSIZES
  mask1 = _mm_set1_epi32(1U);
  mask3 =  _mm_set1_epi32(7U);
  mask5 =  _mm_set1_epi32(31U);
  mask7 =  _mm_set1_epi32(127U);
  mask9 =  _mm_set1_epi32(511U);
  mask11 =  _mm_set1_epi32(2047U);
  mask13 =  _mm_set1_epi32(8191U);
  mask15 =  _mm_set1_epi32(32767U);
  mask17 =  _mm_set1_epi32(131071U);
  mask19 =  _mm_set1_epi32(524287U);
  mask21 =  _mm_set1_epi32(2097151U);
  mask23 =  _mm_set1_epi32(8388607U);
  mask25 =  _mm_set1_epi32(33554431U);
  mask27 =  _mm_set1_epi32(134217727U);
  mask29 =  _mm_set1_epi32(536870911U);
  mask31 =  _mm_set1_epi32(2147483647U);
#endif
  mask2 = _mm_set1_epi32(3U);
  mask4 =  _mm_set1_epi32(15U);
  mask6 =  _mm_set1_epi32(63U);
  mask8 =  _mm_set1_epi32(255U);
  mask10 =  _mm_set1_epi32(1023U);
  mask12 =  _mm_set1_epi32(4095U);
  mask14 =  _mm_set1_epi32(16383U);
  mask16 =  _mm_set1_epi32(65535U);
  mask18 =  _mm_set1_epi32(262143U);
  mask20 =  _mm_set1_epi32(1048575U);
  mask22 =  _mm_set1_epi32(4194303U);
  mask24 =  _mm_set1_epi32(16777215U);
  mask26 =  _mm_set1_epi32(67108863U);
  mask28 =  _mm_set1_epi32(268435455U);
  mask30 =  _mm_set1_epi32(1073741823U);
#endif

  return;
}


#ifdef HAVE_SSE2
static void
unpack_00 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
  __m128i total = _mm_set1_epi32(0U);

  _mm_store_si128(out++, total);
  _mm_store_si128(out++, total);
  _mm_store_si128(out++, total);
  _mm_store_si128(out++, total);
  _mm_store_si128(out++, total);
  _mm_store_si128(out++, total);
  _mm_store_si128(out++, total);
  _mm_store_si128(out++, total);

  return;
}
#else
static void
unpack_00 (UINT4* __restrict__ out, const UINT4* __restrict__ in) {
  int i;

  for (i = 0; i < 64; i++) {
    *out++ = 0;
  }

  return;
}
#endif


#if 0
static void
unpack_01 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg1 = _mm_load_si128(in);
    __m128i InReg2 = InReg1;
    __m128i OutReg1, OutReg2, OutReg3, OutReg4;

    unsigned shift = 0, i;

    for (i = 0; i < 4; i++) {
        OutReg1 = _mm_and_si128(  _mm_srli_epi32(InReg1,shift++) , mask1);
        OutReg2 = _mm_and_si128(  _mm_srli_epi32(InReg2,shift++) , mask1);
        OutReg3 = _mm_and_si128(  _mm_srli_epi32(InReg1,shift++) , mask1);
        OutReg4 = _mm_and_si128(  _mm_srli_epi32(InReg2,shift++) , mask1);
        _mm_store_si128(out++, OutReg1);
        _mm_store_si128(out++, OutReg2);
        _mm_store_si128(out++, OutReg3);
        _mm_store_si128(out++, OutReg4);
    }

    return;
}
#endif

#ifdef ALLOW_ODD_PACKSIZES
static void
unpack_01 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg, total;

    total = /* OutReg = */ _mm_and_si128( InReg , mask1);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,1) , mask1);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,2) , mask1);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,3) , mask1);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask1);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,5) , mask1);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,6) , mask1);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,7) , mask1);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask1);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,9) , mask1);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,10) , mask1);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,11) , mask1);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask1);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,13) , mask1);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,14) , mask1);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,15) ;
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}
#endif


#ifdef HAVE_SSE2
static void
unpack_02_fwd (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg, total;

    total = /* OutReg = */ _mm_and_si128( InReg , mask2);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,2) , mask2);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask2);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,6) , mask2);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask2);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,10) , mask2);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask2);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,14) , mask2);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_02_rev (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;

    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask2);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,18) , mask2);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,20) , mask2);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,22) , mask2);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,24) , mask2);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,26) , mask2);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,28) , mask2);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,30) ;
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

#else

static void
unpack_02 (UINT4* __restrict__ out, const UINT4* __restrict__ in) {
  unsigned int column;
  const UINT4 *bitpack = in;

  for (column = 0; column < 4; column++) {
    in = &(bitpack[column]);

    *out = ( (*in) >>  0  )   % (1U << 2 ) ;
    out++;
    *out = ( (*in) >>  2  )   % (1U << 2 ) ;
    out++;
    *out = ( (*in) >>  4  )   % (1U << 2 ) ;
    out++;
    *out = ( (*in) >>  6  )   % (1U << 2 ) ;
    out++;
    *out = ( (*in) >>  8  )   % (1U << 2 ) ;
    out++;
    *out = ( (*in) >>  10  )   % (1U << 2 ) ;
    out++;
    *out = ( (*in) >>  12  )   % (1U << 2 ) ;
    out++;
    *out = ( (*in) >>  14  )   % (1U << 2 ) ;
    out++;
    *out = ( (*in) >>  16  )   % (1U << 2 ) ;
    out++;
    *out = ( (*in) >>  18  )   % (1U << 2 ) ;
    out++;
    *out = ( (*in) >>  20  )   % (1U << 2 ) ;
    out++;
    *out = ( (*in) >>  22  )   % (1U << 2 ) ;
    out++;
    *out = ( (*in) >>  24  )   % (1U << 2 ) ;
    out++;
    *out = ( (*in) >>  26  )   % (1U << 2 ) ;
    out++;
    *out = ( (*in) >>  28  )   % (1U << 2 ) ;
    out++;
    *out = ( (*in) >>  30  )   % (1U << 2 ) ;
    out++;
  }

  return;
}
#endif



#ifdef ALLOW_ODD_PACKSIZES
static void
unpack_03 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg, total;

    total = /* OutReg = */ _mm_and_si128( InReg , mask3);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,3) , mask3);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,6) , mask3);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,9) , mask3);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask3);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,15) , mask3);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,18) , mask3);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,21) , mask3);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,24) , mask3);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,27) , mask3);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,30) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask1), 3-1));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,1) , mask3);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask3);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,7) , mask3);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,10) , mask3);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,13) , mask3);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}
#endif


#ifdef HAVE_SSE2
static void
unpack_04_fwd (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg, total;

    total = /* OutReg = */ _mm_and_si128( InReg , mask4);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask4);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask4);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask4);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask4);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,20) , mask4);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,24) , mask4);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_04_rev (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;

    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128( InReg , mask4);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask4);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask4);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask4);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask4);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,20) , mask4);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,24) , mask4);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

#else
static void
unpack_04 (UINT4* __restrict__ out, const UINT4* __restrict__ in) {
  UINT4 outer, inwordpointer;
  unsigned int column;
  const UINT4 *bitpack = in;

  for (column = 0; column < 4; column++) {
    in = &(bitpack[column]);

    for (outer = 0; outer < 2 ; outer++) {
      for (inwordpointer = 0; inwordpointer < 32; inwordpointer +=  4) {
	*(out++) = ( (*in) >> inwordpointer )   % (1U << 4 ) ;
      }
      in += 4;
    }
  }

  return;
}
#endif


#ifdef ALLOW_ODD_PACKSIZES
static void
unpack_05 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg, total;

    total = /* OutReg = */ _mm_and_si128( InReg , mask5);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,5) , mask5);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,10) , mask5);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,15) , mask5);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,20) , mask5);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,25) , mask5);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,30) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask3), 5-3));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,3) , mask5);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask5);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,13) , mask5);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,18) , mask5);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,23) , mask5);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask1), 5-1));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,1) , mask5);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,6) , mask5);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,11) , mask5);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}
#endif


#ifdef HAVE_SSE2
static void
unpack_06_fwd (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg, total;

    total = /* OutReg = */ _mm_and_si128( InReg , mask6);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,6) , mask6);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask6);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,18) , mask6);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,24) , mask6);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,30) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask4), 6-4));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask6);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,10) , mask6);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_06_rev (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;

    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask6);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,22) , mask6);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask2), 6-2));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,2) , mask6);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask6);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,14) , mask6);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,20) , mask6);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,26) ;
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}
#else
static void
unpack_06 (UINT4* __restrict__ out, const UINT4* __restrict__ in) {
  unsigned int column;
  const UINT4 *bitpack = in;

  for (column = 0; column < 4; column++) {
    in = &(bitpack[column]);

    *out = ( (*in) >>  0  )   % (1U << 6 ) ;
    out++;
    *out = ( (*in) >>  6  )   % (1U << 6 ) ;
    out++;
    *out = ( (*in) >>  12  )   % (1U << 6 ) ;
    out++;
    *out = ( (*in) >>  18  )   % (1U << 6 ) ;
    out++;
    *out = ( (*in) >>  24  )   % (1U << 6 ) ;
    out++;
    *out = ( (*in) >>  30  )   % (1U << 6 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 4 ))<<( 6 - 4 );
    out++;
    *out = ( (*in) >>  4  )   % (1U << 6 ) ;
    out++;
    *out = ( (*in) >>  10  )   % (1U << 6 ) ;
    out++;
    *out = ( (*in) >>  16  )   % (1U << 6 ) ;
    out++;
    *out = ( (*in) >>  22  )   % (1U << 6 ) ;
    out++;
    *out = ( (*in) >>  28  )   % (1U << 6 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 2 ))<<( 6 - 2 );
    out++;
    *out = ( (*in) >>  2  )   % (1U << 6 ) ;
    out++;
    *out = ( (*in) >>  8  )   % (1U << 6 ) ;
    out++;
    *out = ( (*in) >>  14  )   % (1U << 6 ) ;
    out++;
    *out = ( (*in) >>  20  )   % (1U << 6 ) ;
    out++;
    *out = ( (*in) >>  26  )   % (1U << 6 ) ;
    out++;
  }

  return;
}
#endif


#ifdef ALLOW_ODD_PACKSIZES
static void
unpack_07 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg, total;

    total = /* OutReg = */ _mm_and_si128( InReg , mask7);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,7) , mask7);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,14) , mask7);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,21) , mask7);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask3), 7-3));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,3) , mask7);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,10) , mask7);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,17) , mask7);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,24) , mask7);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,31) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask6), 7-6));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,6) , mask7);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,13) , mask7);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,20) , mask7);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,27) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask2), 7-2));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,2) , mask7);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,9) , mask7);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}
#endif



#ifdef HAVE_SSE2
static void
unpack_08_fwd (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg, total;

    total = /* OutReg = */ _mm_and_si128( InReg , mask8);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask8);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask8);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128( InReg , mask8);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask8);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask8);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_08_rev (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;

    /* ++in; InReg = _mm_load_si128(++in); */
    in += 2;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128( InReg , mask8);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask8);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask8);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128( InReg , mask8);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask8);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask8);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

#else
static void
unpack_08 (UINT4* __restrict__ out, const UINT4* __restrict__ in) {
  UINT4 outer, inwordpointer;
  unsigned int column;
  const UINT4 *bitpack = in;

  for (column = 0; column < 4; column++) {
    in = &(bitpack[column]);

    for (outer = 0; outer < 4; outer++) {
      for (inwordpointer = 0; inwordpointer < 32; inwordpointer += 8) {
	*(out++) = ( (*in) >> inwordpointer )   % (1U << 8 ) ;
      }
      in += 4;
    }
  }

  return;
}
#endif



#ifdef ALLOW_ODD_PACKSIZES
static void
unpack_09 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg, total;

    total = /* OutReg = */ _mm_and_si128( InReg , mask9);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,9) , mask9);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,18) , mask9);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,27) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask4), 9-4));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask9);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,13) , mask9);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,22) , mask9);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,31) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 9-8));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask9);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,17) , mask9);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,26) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask3), 9-3));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,3) , mask9);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask9);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,21) , mask9);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,30) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask7), 9-7));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,7) , mask9);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}
#endif



#ifdef HAVE_SSE2
static void
unpack_10_fwd (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg, total;

    total = /* OutReg = */ _mm_and_si128( InReg , mask10);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,10) , mask10);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,20) , mask10);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,30) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 10-8));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask10);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,18) , mask10);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask6), 10-6));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,6) , mask10);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_10_rev (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;

    /* ++in; InReg = _mm_load_si128(++in); */
    in += 2;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask10);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,26) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask4), 10-4));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask10);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,14) , mask10);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask2), 10-2));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,2) , mask10);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask10);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,22) ;
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

#else
static void
unpack_10 (UINT4* __restrict__ out, const UINT4* __restrict__ in) {
  unsigned int column;
  const UINT4 *bitpack = in;

  for (column = 0; column < 4; column++) {
    in = &(bitpack[column]);

    *out = ( (*in) >>  0  )   % (1U << 10 ) ;
    out++;
    *out = ( (*in) >>  10  )   % (1U << 10 ) ;
    out++;
    *out = ( (*in) >>  20  )   % (1U << 10 ) ;
    out++;
    *out = ( (*in) >>  30  )   % (1U << 10 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 8 ))<<( 10 - 8 );
    out++;
    *out = ( (*in) >>  8  )   % (1U << 10 ) ;
    out++;
    *out = ( (*in) >>  18  )   % (1U << 10 ) ;
    out++;
    *out = ( (*in) >>  28  )   % (1U << 10 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 6 ))<<( 10 - 6 );
    out++;
    *out = ( (*in) >>  6  )   % (1U << 10 ) ;
    out++;
    *out = ( (*in) >>  16  )   % (1U << 10 ) ;
    out++;
    *out = ( (*in) >>  26  )   % (1U << 10 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 4 ))<<( 10 - 4 );
    out++;
    *out = ( (*in) >>  4  )   % (1U << 10 ) ;
    out++;
    *out = ( (*in) >>  14  )   % (1U << 10 ) ;
    out++;
    *out = ( (*in) >>  24  )   % (1U << 10 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 2 ))<<( 10 - 2 );
    out++;
    *out = ( (*in) >>  2  )   % (1U << 10 ) ;
    out++;
    *out = ( (*in) >>  12  )   % (1U << 10 ) ;
    out++;
    *out = ( (*in) >>  22  )   % (1U << 10 ) ;
    out++;

  }
  return;
}
#endif


#ifdef ALLOW_ODD_PACKSIZES
static void
unpack_11 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg, total;

    total = /* OutReg = */ _mm_and_si128( InReg , mask11);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,11) , mask11);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,22) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask1), 11-1));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,1) , mask11);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask11);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,23) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask2), 11-2));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,2) , mask11);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,13) , mask11);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask3), 11-3));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,3) , mask11);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,14) , mask11);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,25) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask4), 11-4));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask11);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,15) , mask11);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,26) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask5), 11-5));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,5) , mask11);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}
#endif


#ifdef HAVE_SSE2
static void
unpack_12_fwd (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg, total;

    total = /* OutReg = */ _mm_and_si128( InReg , mask12);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask12);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask4), 12-4));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask12);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask12);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 12-8));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask12);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_12_rev (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;

    /* ++in; ++in; InReg = _mm_load_si128(++in); */
    in += 3;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128( InReg , mask12);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask12);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask4), 12-4));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask12);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask12);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 12-8));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask12);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

#else
static void
unpack_12 (UINT4* __restrict__ out, const UINT4* __restrict__ in) {
  unsigned int column;
  const UINT4 *bitpack = in;

  for (column = 0; column < 4; column++) {
    in = &(bitpack[column]);

    *out = ( (*in) >>  0  )   % (1U << 12 ) ;
    out++;
    *out = ( (*in) >>  12  )   % (1U << 12 ) ;
    out++;
    *out = ( (*in) >>  24  )   % (1U << 12 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 4 ))<<( 12 - 4 );
    out++;
    *out = ( (*in) >>  4  )   % (1U << 12 ) ;
    out++;
    *out = ( (*in) >>  16  )   % (1U << 12 ) ;
    out++;
    *out = ( (*in) >>  28  )   % (1U << 12 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 8 ))<<( 12 - 8 );
    out++;
    *out = ( (*in) >>  8  )   % (1U << 12 ) ;
    out++;
    *out = ( (*in) >>  20  )   % (1U << 12 ) ;
    in += 4;
    out++;
    *out = ( (*in) >>  0  )   % (1U << 12 ) ;
    out++;
    *out = ( (*in) >>  12  )   % (1U << 12 ) ;
    out++;
    *out = ( (*in) >>  24  )   % (1U << 12 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 4 ))<<( 12 - 4 );
    out++;
    *out = ( (*in) >>  4  )   % (1U << 12 ) ;
    out++;
    *out = ( (*in) >>  16  )   % (1U << 12 ) ;
    out++;
    *out = ( (*in) >>  28  )   % (1U << 12 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 8 ))<<( 12 - 8 );
    out++;
    *out = ( (*in) >>  8  )   % (1U << 12 ) ;
    out++;
    *out = ( (*in) >>  20  )   % (1U << 12 ) ;
    out++;
  }

  return;
}
#endif


#ifdef ALLOW_ODD_PACKSIZES
static void
unpack_13 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg, total;

    total = /* OutReg = */ _mm_and_si128( InReg , mask13);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,13) , mask13);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,26) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask7), 13-7));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,7) , mask13);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask1), 13-1));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,1) , mask13);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,14) , mask13);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,27) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 13-8));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask13);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,21) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask2), 13-2));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,2) , mask13);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,15) , mask13);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask9), 13-9));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,9) , mask13);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,22) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask3), 13-3));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,3) , mask13);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}
#endif


#ifdef HAVE_SSE2
static void
unpack_14_fwd (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg, total;

    total = /* OutReg = */ _mm_and_si128( InReg , mask14);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,14) , mask14);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask10), 14-10));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,10) , mask14);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask6), 14-6));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,6) , mask14);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask2), 14-2));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,2) , mask14);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_14_rev (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;

    /* ++in; ++in; InReg = _mm_load_si128(++in); */
    in += 3;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128(  _mm_srli_epi32(InReg,16) , mask14);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,30) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask12), 14-12));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask14);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,26) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 14-8));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask14);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,22) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask4), 14-4));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask14);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,18) ;
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

#else
static void
unpack_14 (UINT4* __restrict__ out, const UINT4* __restrict__ in) {
  unsigned int column;
  const UINT4 *bitpack = in;

  for (column = 0; column < 4; column++) {
    in = &(bitpack[column]);

    *out = ( (*in) >>  0  )   % (1U << 14 ) ;
    out++;
    *out = ( (*in) >>  14  )   % (1U << 14 ) ;
    out++;
    *out = ( (*in) >>  28  )   % (1U << 14 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 10 ))<<( 14 - 10 );
    out++;
    *out = ( (*in) >>  10  )   % (1U << 14 ) ;
    out++;
    *out = ( (*in) >>  24  )   % (1U << 14 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 6 ))<<( 14 - 6 );
    out++;
    *out = ( (*in) >>  6  )   % (1U << 14 ) ;
    out++;
    *out = ( (*in) >>  20  )   % (1U << 14 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 2 ))<<( 14 - 2 );
    out++;
    *out = ( (*in) >>  2  )   % (1U << 14 ) ;
    out++;
    *out = ( (*in) >>  16  )   % (1U << 14 ) ;
    out++;
    *out = ( (*in) >>  30  )   % (1U << 14 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 12 ))<<( 14 - 12 );
    out++;
    *out = ( (*in) >>  12  )   % (1U << 14 ) ;
    out++;
    *out = ( (*in) >>  26  )   % (1U << 14 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 8 ))<<( 14 - 8 );
    out++;
    *out = ( (*in) >>  8  )   % (1U << 14 ) ;
    out++;
    *out = ( (*in) >>  22  )   % (1U << 14 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 4 ))<<( 14 - 4 );
    out++;
    *out = ( (*in) >>  4  )   % (1U << 14 ) ;
    out++;
    *out = ( (*in) >>  18  )   % (1U << 14 ) ;
    out++;
  }

  return;
}
#endif



#ifdef ALLOW_ODD_PACKSIZES
static void
unpack_15 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg, total;

    total = /* OutReg = */ _mm_and_si128( InReg , mask15);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,15) , mask15);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,30) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask13), 15-13));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,13) , mask15);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask11), 15-11));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,11) , mask15);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,26) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask9), 15-9));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,9) , mask15);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask7), 15-7));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,7) , mask15);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,22) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask5), 15-5));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,5) , mask15);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask3), 15-3));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,3) , mask15);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,18) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask1), 15-1));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,1) , mask15);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}
#endif


#ifdef HAVE_SSE2
static void
unpack_16_fwd (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg, total;

    total = /* OutReg = */ _mm_and_si128( InReg , mask16);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128( InReg , mask16);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128( InReg , mask16);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128( InReg , mask16);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_16_rev (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;

    /* ++in; ++in; ++in; InReg = _mm_load_si128(++in); */
    in += 4;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128( InReg , mask16);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128( InReg , mask16);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128( InReg , mask16);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128( InReg , mask16);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

#else
static void
unpack_16 (UINT4* __restrict__ out, const UINT4* __restrict__ in) {
  UINT4 outer, inwordpointer;
  unsigned int column;
  const UINT4 *bitpack = in;

  for (column = 0; column < 4; column++) {
    in = &(bitpack[column]);

    for (outer = 0; outer < 8; outer++) {
      for(inwordpointer =  0; inwordpointer <32; inwordpointer += 16) {
	*(out++) = ( (*in) >> inwordpointer )   % (1U << 16 ) ;
      }
      in += 4;
    }
  }

  return;
}
#endif


#ifdef ALLOW_ODD_PACKSIZES
static void
unpack_17 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg, total;

    total = /* OutReg = */ _mm_and_si128( InReg , mask17);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,17) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask2), 17-2));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,2) , mask17);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,19) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask4), 17-4));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask17);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,21) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask6), 17-6));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,6) , mask17);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,23) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 17-8));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask17);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,25) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask10), 17-10));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,10) , mask17);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,27) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask12), 17-12));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask17);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,29) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask14), 17-14));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,14) , mask17);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,31) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask16), 17-16));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}
#endif


#ifdef HAVE_SSE2
static void
unpack_18_fwd (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg, total;

    total = /* OutReg = */ _mm_and_si128( InReg , mask18);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,18) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask4), 18-4));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask18);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,22) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 18-8));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask18);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,26) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask12), 18-12));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask18);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,30) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask16), 18-16));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_18_rev (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;

    /* ++in; ++in; ++in; InReg = _mm_load_si128(++in); */
    in += 4;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(total /*OutReg*/, _mm_slli_epi32(_mm_and_si128(InReg, mask2), 18-2));
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,2) , mask18);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask6), 18-6));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,6) , mask18);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask10), 18-10));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,10) , mask18);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask14), 18-14));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,14) ;
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

#else
static void
unpack_18 (UINT4* __restrict__ out, const UINT4* __restrict__ in) {
  unsigned int column;
  const UINT4 *bitpack = in;

  for (column = 0; column < 4; column++) {
    in = &(bitpack[column]);

    *out = ( (*in) >>  0  )   % (1U << 18 ) ;
    out++;
    *out = ( (*in) >>  18  )   % (1U << 18 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 4 ))<<( 18 - 4 );
    out++;
    *out = ( (*in) >>  4  )   % (1U << 18 ) ;
    out++;
    *out = ( (*in) >>  22  )   % (1U << 18 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 8 ))<<( 18 - 8 );
    out++;
    *out = ( (*in) >>  8  )   % (1U << 18 ) ;
    out++;
    *out = ( (*in) >>  26  )   % (1U << 18 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 12 ))<<( 18 - 12 );
    out++;
    *out = ( (*in) >>  12  )   % (1U << 18 ) ;
    out++;
    *out = ( (*in) >>  30  )   % (1U << 18 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 16 ))<<( 18 - 16 );
    out++;
    *out = ( (*in) >>  16  )   % (1U << 18 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 2 ))<<( 18 - 2 );
    out++;
    *out = ( (*in) >>  2  )   % (1U << 18 ) ;
    out++;
    *out = ( (*in) >>  20  )   % (1U << 18 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 6 ))<<( 18 - 6 );
    out++;
    *out = ( (*in) >>  6  )   % (1U << 18 ) ;
    out++;
    *out = ( (*in) >>  24  )   % (1U << 18 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 10 ))<<( 18 - 10 );
    out++;
    *out = ( (*in) >>  10  )   % (1U << 18 ) ;
    out++;
    *out = ( (*in) >>  28  )   % (1U << 18 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 14 ))<<( 18 - 14 );
    out++;
    *out = ( (*in) >>  14  )   % (1U << 18 ) ;
    out++;
  }

  return;
}
#endif


#ifdef ALLOW_ODD_PACKSIZES
static void
unpack_19 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg, total;

    total = /* OutReg = */ _mm_and_si128( InReg , mask19);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,19) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask6), 19-6));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,6) , mask19);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,25) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask12), 19-12));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,12) , mask19);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,31) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask18), 19-18));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,18) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask5), 19-5));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,5) , mask19);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask11), 19-11));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,11) , mask19);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,30) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask17), 19-17));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,17) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask4), 19-4));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask19);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,23) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask10), 19-10));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,10) , mask19);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,29) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask16), 19-16));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}
#endif


#ifdef HAVE_SSE2
static void
unpack_20_fwd (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg, total;

    total = /* OutReg = */ _mm_and_si128( InReg , mask20);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 20-8));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask20);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask16), 20-16));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask4), 20-4));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask20);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask12), 20-12));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,12) ;
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_20_rev (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;

    /* ++in; ++in; ++in; ++in; InReg = _mm_load_si128(++in); */
    in += 5;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128( InReg , mask20);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 20-8));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask20);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask16), 20-16));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask4), 20-4));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask20);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask12), 20-12));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,12) ;
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}


#else
static void
unpack_20 (UINT4* __restrict__ out, const UINT4* __restrict__ in) {
  unsigned int column;
  const UINT4 *bitpack = in;

  for (column = 0; column < 4; column++) {
    in = &(bitpack[column]);

    *out = ( (*in) >>  0  )   % (1U << 20 ) ;
    out++;
    *out = ( (*in) >>  20  )   % (1U << 20 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 8 ))<<( 20 - 8 );
    out++;
    *out = ( (*in) >>  8  )   % (1U << 20 ) ;
    out++;
    *out = ( (*in) >>  28  )   % (1U << 20 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 16 ))<<( 20 - 16 );
    out++;
    *out = ( (*in) >>  16  )   % (1U << 20 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 4 ))<<( 20 - 4 );
    out++;
    *out = ( (*in) >>  4  )   % (1U << 20 ) ;
    out++;
    *out = ( (*in) >>  24  )   % (1U << 20 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 12 ))<<( 20 - 12 );
    out++;
    *out = ( (*in) >>  12  )   % (1U << 20 ) ;
    in += 4;
    out++;
    *out = ( (*in) >>  0  )   % (1U << 20 ) ;
    out++;
    *out = ( (*in) >>  20  )   % (1U << 20 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 8 ))<<( 20 - 8 );
    out++;
    *out = ( (*in) >>  8  )   % (1U << 20 ) ;
    out++;
    *out = ( (*in) >>  28  )   % (1U << 20 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 16 ))<<( 20 - 16 );
    out++;
    *out = ( (*in) >>  16  )   % (1U << 20 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 4 ))<<( 20 - 4 );
    out++;
    *out = ( (*in) >>  4  )   % (1U << 20 ) ;
    out++;
    *out = ( (*in) >>  24  )   % (1U << 20 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 12 ))<<( 20 - 12 );
    out++;
    *out = ( (*in) >>  12  )   % (1U << 20 ) ;
    out++;
  }

  return;
}
#endif


#ifdef ALLOW_ODD_PACKSIZES
static void
unpack_21 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg, total;

    total = /* OutReg = */ _mm_and_si128( InReg , mask21);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,21) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask10), 21-10));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,10) , mask21);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,31) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask20), 21-20));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask9), 21-9));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,9) , mask21);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,30) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask19), 21-19));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,19) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 21-8));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask21);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,29) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask18), 21-18));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,18) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask7), 21-7));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,7) , mask21);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask17), 21-17));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,17) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask6), 21-6));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,6) , mask21);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,27) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask16), 21-16));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}
#endif


#ifdef HAVE_SSE2
static void
unpack_22_fwd (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg, total;

    total = /* OutReg = */ _mm_and_si128( InReg , mask22);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,22) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask12), 22-12));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,12) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask2), 22-2));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,2) , mask22);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask14), 22-14));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,14) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask4), 22-4));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask22);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,26) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask16), 22-16));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_22_rev (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;

    /* ++in; ++in; ++in; ++in; InReg = _mm_load_si128(++in); */
    in += 5;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */  _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(total /*OutReg*/, _mm_slli_epi32(_mm_and_si128(InReg, mask6), 22-6));
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,6) , mask22);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask18), 22-18));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,18) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 22-8));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,8) , mask22);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,30) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask20), 22-20));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask10), 22-10));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,10) ;
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

#else
static void
unpack_22 (UINT4* __restrict__ out, const UINT4* __restrict__ in) {
  unsigned int column;
  const UINT4 *bitpack = in;

  for (column = 0; column < 4; column++) {
    in = &(bitpack[column]);

    *out = ( (*in) >>  0  )   % (1U << 22 ) ;
    out++;
    *out = ( (*in) >>  22  )   % (1U << 22 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 12 ))<<( 22 - 12 );
    out++;
    *out = ( (*in) >>  12  )   % (1U << 22 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 2 ))<<( 22 - 2 );
    out++;
    *out = ( (*in) >>  2  )   % (1U << 22 ) ;
    out++;
    *out = ( (*in) >>  24  )   % (1U << 22 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 14 ))<<( 22 - 14 );
    out++;
    *out = ( (*in) >>  14  )   % (1U << 22 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 4 ))<<( 22 - 4 );
    out++;
    *out = ( (*in) >>  4  )   % (1U << 22 ) ;
    out++;
    *out = ( (*in) >>  26  )   % (1U << 22 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 16 ))<<( 22 - 16 );
    out++;
    *out = ( (*in) >>  16  )   % (1U << 22 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 6 ))<<( 22 - 6 );
    out++;
    *out = ( (*in) >>  6  )   % (1U << 22 ) ;
    out++;
    *out = ( (*in) >>  28  )   % (1U << 22 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 18 ))<<( 22 - 18 );
    out++;
    *out = ( (*in) >>  18  )   % (1U << 22 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 8 ))<<( 22 - 8 );
    out++;
    *out = ( (*in) >>  8  )   % (1U << 22 ) ;
    out++;
    *out = ( (*in) >>  30  )   % (1U << 22 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 20 ))<<( 22 - 20 );
    out++;
    *out = ( (*in) >>  20  )   % (1U << 22 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 10 ))<<( 22 - 10 );
    out++;
    *out = ( (*in) >>  10  )   % (1U << 22 ) ;
    out++;
  }

  return;
}
#endif


#ifdef ALLOW_ODD_PACKSIZES
static void
unpack_23 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg, total;

    total = /* OutReg = */ _mm_and_si128( InReg , mask23);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,23) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask14), 23-14));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,14) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask5), 23-5));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,5) , mask23);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask19), 23-19));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,19) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask10), 23-10));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,10) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask1), 23-1));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,1) , mask23);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask15), 23-15));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,15) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask6), 23-6));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,6) , mask23);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,29) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask20), 23-20));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask11), 23-11));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,11) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask2), 23-2));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,2) , mask23);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,25) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask16), 23-16));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}
#endif


#ifdef HAVE_SSE2
static void
unpack_24_fwd (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg, total;

    total = /* OutReg = */ _mm_and_si128( InReg , mask24);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask16), 24-16));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 24-8));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,8) ;
    InReg = _mm_load_si128(++in);

    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128( InReg , mask24);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask16), 24-16));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 24-8));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,8) ;
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_24_rev (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;

    /* ++in; ++in; ++in; ++in; ++in; InReg = _mm_load_si128(++in); */
    in += 6;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128( InReg , mask24);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask16), 24-16));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 24-8));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,8) ;
    InReg = _mm_load_si128(++in);

    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128( InReg , mask24);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask16), 24-16));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 24-8));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,8) ;
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

#else
static void
unpack_24 (UINT4* __restrict__ out, const UINT4* __restrict__ in) {
  unsigned int column;
  const UINT4 *bitpack = in;

  for (column = 0; column < 4; column++) {
    in = &(bitpack[column]);

    *out = ( (*in) >>  0  )   % (1U << 24 ) ;
    out++;
    *out = ( (*in) >>  24  )   % (1U << 24 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 16 ))<<( 24 - 16 );
    out++;
    *out = ( (*in) >>  16  )   % (1U << 24 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 8 ))<<( 24 - 8 );
    out++;
    *out = ( (*in) >>  8  )   % (1U << 24 ) ;
    in += 4;
    out++;
    *out = ( (*in) >>  0  )   % (1U << 24 ) ;
    out++;
    *out = ( (*in) >>  24  )   % (1U << 24 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 16 ))<<( 24 - 16 );
    out++;
    *out = ( (*in) >>  16  )   % (1U << 24 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 8 ))<<( 24 - 8 );
    out++;
    *out = ( (*in) >>  8  )   % (1U << 24 ) ;
    in += 4;
    out++;
    *out = ( (*in) >>  0  )   % (1U << 24 ) ;
    out++;
    *out = ( (*in) >>  24  )   % (1U << 24 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 16 ))<<( 24 - 16 );
    out++;
    *out = ( (*in) >>  16  )   % (1U << 24 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 8 ))<<( 24 - 8 );
    out++;
    *out = ( (*in) >>  8  )   % (1U << 24 ) ;
    in += 4;
    out++;
    *out = ( (*in) >>  0  )   % (1U << 24 ) ;
    out++;
    *out = ( (*in) >>  24  )   % (1U << 24 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 16 ))<<( 24 - 16 );
    out++;
    *out = ( (*in) >>  16  )   % (1U << 24 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 8 ))<<( 24 - 8 );
    out++;
    *out = ( (*in) >>  8  )   % (1U << 24 ) ;
    out++;
  }

  return;
}
#endif


#ifdef ALLOW_ODD_PACKSIZES
static void
unpack_25 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg, total;

    total = /* OutReg = */ _mm_and_si128( InReg , mask25);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,25) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask18), 25-18));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,18) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask11), 25-11));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,11) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask4), 25-4));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask25);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,29) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask22), 25-22));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,22) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask15), 25-15));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,15) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 25-8));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,8) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask1), 25-1));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,1) , mask25);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,26) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask19), 25-19));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,19) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask12), 25-12));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,12) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask5), 25-5));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,5) , mask25);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,30) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask23), 25-23));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,23) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask16), 25-16));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}
#endif


#ifdef HAVE_SSE2
static void
unpack_26_fwd (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg, total;

    total = /* OutReg = */ _mm_and_si128( InReg , mask26);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,26) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask20), 26-20));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask14), 26-14));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,14) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 26-8));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,8) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask2), 26-2));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,2) , mask26);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask22), 26-22));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,22) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask16), 26-16));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_26_rev (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;

    /* ++in; ++in; ++in; ++in; ++in; InReg = _mm_load_si128(++in); */
    in += 6;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(total /*OutReg*/, _mm_slli_epi32(_mm_and_si128(InReg, mask10), 26-10));
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,10) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask4), 26-4));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask26);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,30) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask24), 26-24));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask18), 26-18));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,18) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask12), 26-12));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,12) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask6), 26-6));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,6) ;
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

#else
static void
unpack_26 (UINT4* __restrict__ out, const UINT4* __restrict__ in) {
  unsigned int column;
  const UINT4 *bitpack = in;

  for (column = 0; column < 4; column++) {
    in = &(bitpack[column]);

    *out = ( (*in) >>  0  )   % (1U << 26 ) ;
    out++;
    *out = ( (*in) >>  26  )   % (1U << 26 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 20 ))<<( 26 - 20 );
    out++;
    *out = ( (*in) >>  20  )   % (1U << 26 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 14 ))<<( 26 - 14 );
    out++;
    *out = ( (*in) >>  14  )   % (1U << 26 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 8 ))<<( 26 - 8 );
    out++;
    *out = ( (*in) >>  8  )   % (1U << 26 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 2 ))<<( 26 - 2 );
    out++;
    *out = ( (*in) >>  2  )   % (1U << 26 ) ;
    out++;
    *out = ( (*in) >>  28  )   % (1U << 26 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 22 ))<<( 26 - 22 );
    out++;
    *out = ( (*in) >>  22  )   % (1U << 26 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 16 ))<<( 26 - 16 );
    out++;
    *out = ( (*in) >>  16  )   % (1U << 26 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 10 ))<<( 26 - 10 );
    out++;
    *out = ( (*in) >>  10  )   % (1U << 26 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 4 ))<<( 26 - 4 );
    out++;
    *out = ( (*in) >>  4  )   % (1U << 26 ) ;
    out++;
    *out = ( (*in) >>  30  )   % (1U << 26 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 24 ))<<( 26 - 24 );
    out++;
    *out = ( (*in) >>  24  )   % (1U << 26 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 18 ))<<( 26 - 18 );
    out++;
    *out = ( (*in) >>  18  )   % (1U << 26 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 12 ))<<( 26 - 12 );
    out++;
    *out = ( (*in) >>  12  )   % (1U << 26 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 6 ))<<( 26 - 6 );
    out++;
    *out = ( (*in) >>  6  )   % (1U << 26 ) ;
    out++;
  }
  
  return;
}
#endif


#ifdef ALLOW_ODD_PACKSIZES
static void
unpack_27 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg, total;

    total = /* OutReg = */ _mm_and_si128( InReg , mask27);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,27) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask22), 27-22));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,22) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask17), 27-17));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,17) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask12), 27-12));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,12) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask7), 27-7));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,7) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask2), 27-2));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,2) , mask27);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,29) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask24), 27-24));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask19), 27-19));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,19) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask14), 27-14));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,14) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask9), 27-9));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,9) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask4), 27-4));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,4) , mask27);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,31) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask26), 27-26));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,26) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask21), 27-21));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,21) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask16), 27-16));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}
#endif


#ifdef HAVE_SSE2
static void
unpack_28_fwd (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg, total;

    total = /* OutReg = */ _mm_and_si128( InReg , mask28);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask24), 28-24));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask20), 28-20));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask16), 28-16));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask12), 28-12));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,12) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 28-8));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,8) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask4), 28-4));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,4) ;
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_28_rev (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;

    /* ++in; ++in; ++in; ++in; ++in; ++in; InReg = _mm_load_si128(++in); */
    in += 7;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */ _mm_and_si128( InReg , mask28);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask24), 28-24));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask20), 28-20));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask16), 28-16));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask12), 28-12));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,12) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 28-8));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,8) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask4), 28-4));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,4) ;
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

#else
static void
unpack_28 (UINT4* __restrict__ out, const UINT4* __restrict__ in) {
  unsigned int column;
  const UINT4 *bitpack = in;

  for (column = 0; column < 4; column++) {
    in = &(bitpack[column]);

    *out = ( (*in) >>  0  )   % (1U << 28 ) ;
    out++;
    *out = ( (*in) >>  28  )   % (1U << 28 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 24 ))<<( 28 - 24 );
    out++;
    *out = ( (*in) >>  24  )   % (1U << 28 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 20 ))<<( 28 - 20 );
    out++;
    *out = ( (*in) >>  20  )   % (1U << 28 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 16 ))<<( 28 - 16 );
    out++;
    *out = ( (*in) >>  16  )   % (1U << 28 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 12 ))<<( 28 - 12 );
    out++;
    *out = ( (*in) >>  12  )   % (1U << 28 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 8 ))<<( 28 - 8 );
    out++;
    *out = ( (*in) >>  8  )   % (1U << 28 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 4 ))<<( 28 - 4 );
    out++;
    *out = ( (*in) >>  4  )   % (1U << 28 ) ;
    in += 4;
    out++;
    *out = ( (*in) >>  0  )   % (1U << 28 ) ;
    out++;
    *out = ( (*in) >>  28  )   % (1U << 28 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 24 ))<<( 28 - 24 );
    out++;
    *out = ( (*in) >>  24  )   % (1U << 28 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 20 ))<<( 28 - 20 );
    out++;
    *out = ( (*in) >>  20  )   % (1U << 28 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 16 ))<<( 28 - 16 );
    out++;
    *out = ( (*in) >>  16  )   % (1U << 28 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 12 ))<<( 28 - 12 );
    out++;
    *out = ( (*in) >>  12  )   % (1U << 28 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 8 ))<<( 28 - 8 );
    out++;
    *out = ( (*in) >>  8  )   % (1U << 28 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 4 ))<<( 28 - 4 );
    out++;
    *out = ( (*in) >>  4  )   % (1U << 28 ) ;
    out++;
  }

  return;
}
#endif


#ifdef ALLOW_ODD_PACKSIZES
static void
unpack_29 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg, total;

    total = /* OutReg = */ _mm_and_si128( InReg , mask29);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,29) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask26), 29-26));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,26) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask23), 29-23));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,23) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask20), 29-20));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask17), 29-17));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,17) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask14), 29-14));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,14) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask11), 29-11));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,11) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 29-8));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,8) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask5), 29-5));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,5) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask2), 29-2));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg = _mm_and_si128(  _mm_srli_epi32(InReg,2) , mask29);
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,31) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask28), 29-28));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask25), 29-25));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,25) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask22), 29-22));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,22) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask19), 29-19));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,19) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask16), 29-16));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}
#endif


#ifdef HAVE_SSE2
static void
unpack_30_fwd (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg, total;

    total = /* OutReg = */ _mm_and_si128( InReg , mask30);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,30) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask28), 30-28));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask26), 30-26));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,26) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask24), 30-24));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask22), 30-22));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,22) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask20), 30-20));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask18), 30-18));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,18) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask16), 30-16));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_30_rev (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg;
    __m128i OutReg, total;

    /* ++in; ++in; ++in; ++in; ++in; ++in; InReg = _mm_load_si128(++in); */
    in += 7;
    InReg = _mm_load_si128(in);


    total = /* OutReg = */   _mm_srli_epi32(InReg,16) ;
    InReg = _mm_load_si128(++in);

    total = /* OutReg = */ _mm_or_si128(total /*OutReg*/, _mm_slli_epi32(_mm_and_si128(InReg, mask14), 30-14));
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,14) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask12), 30-12));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,12) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask10), 30-10));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,10) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask8), 30-8));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,8) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask6), 30-6));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,6) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask4), 30-4));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,4) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask2), 30-2));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,2) ;
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}

#else
static void
unpack_30 (UINT4* __restrict__ out, const UINT4* __restrict__ in) {
  unsigned int column;
  const UINT4 *bitpack = in;

  for (column = 0; column < 4; column++) {
    in = &(bitpack[column]);

    *out = ( (*in) >>  0  )   % (1U << 30 ) ;
    out++;
    *out = ( (*in) >>  30  )   % (1U << 30 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 28 ))<<( 30 - 28 );
    out++;
    *out = ( (*in) >>  28  )   % (1U << 30 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 26 ))<<( 30 - 26 );
    out++;
    *out = ( (*in) >>  26  )   % (1U << 30 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 24 ))<<( 30 - 24 );
    out++;
    *out = ( (*in) >>  24  )   % (1U << 30 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 22 ))<<( 30 - 22 );
    out++;
    *out = ( (*in) >>  22  )   % (1U << 30 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 20 ))<<( 30 - 20 );
    out++;
    *out = ( (*in) >>  20  )   % (1U << 30 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 18 ))<<( 30 - 18 );
    out++;
    *out = ( (*in) >>  18  )   % (1U << 30 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 16 ))<<( 30 - 16 );
    out++;
    *out = ( (*in) >>  16  )   % (1U << 30 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 14 ))<<( 30 - 14 );
    out++;
    *out = ( (*in) >>  14  )   % (1U << 30 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 12 ))<<( 30 - 12 );
    out++;
    *out = ( (*in) >>  12  )   % (1U << 30 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 10 ))<<( 30 - 10 );
    out++;
    *out = ( (*in) >>  10  )   % (1U << 30 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 8 ))<<( 30 - 8 );
    out++;
    *out = ( (*in) >>  8  )   % (1U << 30 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 6 ))<<( 30 - 6 );
    out++;
    *out = ( (*in) >>  6  )   % (1U << 30 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 4 ))<<( 30 - 4 );
    out++;
    *out = ( (*in) >>  4  )   % (1U << 30 ) ;
    in += 4;
    *out |= ((*in) % (1U<< 2 ))<<( 30 - 2 );
    out++;
    *out = ( (*in) >>  2  )   % (1U << 30 ) ;
    out++;
  }

  return;
}
#endif


#ifdef ALLOW_ODD_PACKSIZES
static void
unpack_31 (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i InReg = _mm_load_si128(in);
    __m128i OutReg, total;

    total = /* OutReg = */ _mm_and_si128( InReg , mask31);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,31) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask30), 31-30));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,30) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask29), 31-29));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,29) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask28), 31-28));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,28) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask27), 31-27));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,27) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask26), 31-26));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,26) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask25), 31-25));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,25) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask24), 31-24));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,24) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask23), 31-23));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,23) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask22), 31-22));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,22) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask21), 31-21));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,21) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask20), 31-20));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,20) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask19), 31-19));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,19) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask18), 31-18));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,18) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask17), 31-17));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    OutReg =   _mm_srli_epi32(InReg,17) ;
    InReg = _mm_load_si128(++in);

    OutReg = _mm_or_si128(OutReg, _mm_slli_epi32(_mm_and_si128(InReg, mask16), 31-16));
    total = _mm_add_epi32(total, OutReg);
    _mm_store_si128(out++, total);

    return;
}
#endif



#ifdef HAVE_SSE2
static void
unpack_32_fwd (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i total;

    total = _mm_load_si128(in++);
    _mm_store_si128(out++, total);

    total = _mm_add_epi32(total, _mm_load_si128(in++));
    _mm_store_si128(out++, total);

    total = _mm_add_epi32(total, _mm_load_si128(in++));
    _mm_store_si128(out++, total);

    total = _mm_add_epi32(total, _mm_load_si128(in++));
    _mm_store_si128(out++, total);

    total = _mm_add_epi32(total, _mm_load_si128(in++));
    _mm_store_si128(out++, total);

    total = _mm_add_epi32(total, _mm_load_si128(in++));
    _mm_store_si128(out++, total);

    total = _mm_add_epi32(total, _mm_load_si128(in++));
    _mm_store_si128(out++, total);

    total = _mm_add_epi32(total, _mm_load_si128(in++));
    _mm_store_si128(out++, total);

    return;
}

static void
unpack_32_rev (__m128i* __restrict__ out, const __m128i* __restrict__ in) {
    __m128i total;

    in += 8;


    total = _mm_add_epi32(total, _mm_load_si128(in++));
    _mm_store_si128(out++, total);

    total = _mm_add_epi32(total, _mm_load_si128(in++));
    _mm_store_si128(out++, total);

    total = _mm_add_epi32(total, _mm_load_si128(in++));
    _mm_store_si128(out++, total);

    total = _mm_add_epi32(total, _mm_load_si128(in++));
    _mm_store_si128(out++, total);

    total = _mm_add_epi32(total, _mm_load_si128(in++));
    _mm_store_si128(out++, total);

    total = _mm_add_epi32(total, _mm_load_si128(in++));
    _mm_store_si128(out++, total);

    total = _mm_add_epi32(total, _mm_load_si128(in++));
    _mm_store_si128(out++, total);

    total = _mm_add_epi32(total, _mm_load_si128(in++));
    _mm_store_si128(out++, total);

    return;
}

#else
static void
unpack_32 (UINT4* __restrict__ out, const UINT4* __restrict__ in) {
  unsigned int column;
  const UINT4 *bitpack = in;

  for (column = 0; column < 4; column++) {
    in = &(bitpack[column]);

    *out = ( (*in) >>  0  )   ;
    in += 4;
    out++;
    *out = ( (*in) >>  0  )   ;
    in += 4;
    out++;
    *out = ( (*in) >>  0  )   ;
    in += 4;
    out++;
    *out = ( (*in) >>  0  )   ;
    in += 4;
    out++;
    *out = ( (*in) >>  0  )   ;
    in += 4;
    out++;
    *out = ( (*in) >>  0  )   ;
    in += 4;
    out++;
    *out = ( (*in) >>  0  )   ;
    in += 4;
    out++;
    *out = ( (*in) >>  0  )   ;
    in += 4;
    out++;
    *out = ( (*in) >>  0  )   ;
    in += 4;
    out++;
    *out = ( (*in) >>  0  )   ;
    in += 4;
    out++;
    *out = ( (*in) >>  0  )   ;
    in += 4;
    out++;
    *out = ( (*in) >>  0  )   ;
    in += 4;
    out++;
    *out = ( (*in) >>  0  )   ;
    in += 4;
    out++;
    *out = ( (*in) >>  0  )   ;
    in += 4;
    out++;
    *out = ( (*in) >>  0  )   ;
    in += 4;
    out++;
    *out = ( (*in) >>  0  )   ;
    out++;
  }

  return;
}
#endif



#ifdef HAVE_SSE2
typedef void (*Unpacker_T) (__m128i* __restrict__, const __m128i* __restrict__);
#else
typedef void (*Unpacker_T) (UINT4* __restrict__, const UINT4* __restrict__);
#endif


#ifdef ALLOW_ODD_PACKSIZES
static Unpacker_T unpacker_table[33] =
  {unpack_00,
   unpack_01, unpack_02, unpack_03, unpack_04,
   unpack_05, unpack_06, unpack_07, unpack_08,
   unpack_09, unpack_10, unpack_11, unpack_12,
   unpack_13, unpack_14, unpack_15, unpack_16,
   unpack_17, unpack_18, unpack_19, unpack_20,
   unpack_21, unpack_22, unpack_23, unpack_24,
   unpack_25, unpack_26, unpack_27, unpack_28,
   unpack_29, unpack_30, unpack_31, unpack_32};
#else
#ifdef HAVE_SSE2
static Unpacker_T unpacker_table[34] =
  {unpack_00, unpack_00,
   unpack_02_fwd, unpack_02_rev, unpack_04_fwd, unpack_04_rev,
   unpack_06_fwd, unpack_06_rev, unpack_08_fwd, unpack_08_rev,
   unpack_10_fwd, unpack_10_rev, unpack_12_fwd, unpack_12_rev,
   unpack_14_fwd, unpack_14_rev, unpack_16_fwd, unpack_16_rev,
   unpack_18_fwd, unpack_18_rev, unpack_20_fwd, unpack_20_rev,
   unpack_22_fwd, unpack_22_rev, unpack_24_fwd, unpack_24_rev,
   unpack_26_fwd, unpack_26_rev, unpack_28_fwd, unpack_28_rev,
   unpack_30_fwd, unpack_30_rev, unpack_32_fwd, unpack_32_rev};
#else
static Unpacker_T unpacker_table[33] =
  {unpack_00,
   unpack_00, unpack_02, unpack_00, unpack_04,
   unpack_00, unpack_06, unpack_00, unpack_08,
   unpack_00, unpack_10, unpack_00, unpack_12,
   unpack_00, unpack_14, unpack_00, unpack_16,
   unpack_00, unpack_18, unpack_00, unpack_20,
   unpack_00, unpack_22, unpack_00, unpack_24,
   unpack_00, unpack_26, unpack_00, unpack_28,
   unpack_00, unpack_30, unpack_00, unpack_32};
#endif
#endif


#define METAINFO_SIZE 2
#define TOTAL_ROWS 16

Positionsptr_T
Bitpack64_offsetptr (Positionsptr_T *end0, Storedoligomer_T oligo) {
  UINT4 *info, nwritten;
  Positionsptr_T offset0, offset1;
  int packsize, remainder;
#ifdef HAVE_SSE2
  __m128i diffs[9], *bitpack;
  UINT4 *_diffs;
#else
  Positionsptr_T ptr;
  int column, row, k;
  UINT4 diffs[65], *bitpack;
#endif
#ifdef DEBUG
  Positionsptr_T offsets[65];
  int i;
#endif

  info = &(bitpackptrs[oligo/OFFSETSCOMP_BLOCKSIZE * METAINFO_SIZE]);

  nwritten = info[0];
#ifdef HAVE_SSE2  
  bitpack = (__m128i *) &(offsetscomp[nwritten]);
#else
  bitpack = (UINT4 *) &(offsetscomp[nwritten]);
#endif

#ifdef ALLOW_ODD_PACKSIZES
  /* Cannot allow odd packsizes, because they destroy to one-to-one
     correspondence between nwritten and packsize */
  abort();
#else
  /* packsize 1 => 1 _mm128i => nwritten 4 UINT4s */
  /* packsize 2 => 1 _mm128i => nwritten 4 UINT4s */
  packsize = (info[2] - nwritten)/2;
#endif

  debug(printf("nwritten0: %u, offset: %u, next_nwritten: %u, packsize %d, remainder %d\n",
	       info[0],info[1],info[2],packsize,remainder));

#ifdef HAVE_SSE2
  if ((remainder = oligo % OFFSETSCOMP_BLOCKSIZE) < 32) {
    /* Unpack fwd 32 cumulative sums under SIMD */
    (unpacker_table[packsize])(&(diffs[1]),bitpack);
    offset0 = info[1];

#ifdef DEBUG
    printf("oligo: %08X, remainder %d, offset0 %u\n",oligo,remainder,offset0);
    printf("bitpack:\n");
    for (i = 0; i < packsize/4; i++) {
      print_vector_hex(bitpack[i]);
    }
    printf("\n");

    /* diffs[0] is used only to hold the result 0 when remainder is 0 */
    for (i = 1; i <= 8; i++) {
      print_vector(diffs[i]);
    }
    printf("end of diffs\n");
#endif  

    _diffs = (UINT4 *) &diffs;
    _diffs[3] = 0;

    *end0 = offset0 + _diffs[4+remainder];

    return offset0 + _diffs[3+remainder];

  } else {
    /* Unpack rev 32 cumulative sums under SIMD */
    (unpacker_table[packsize+1])(&(diffs[1]),bitpack);
    offset1 = info[METAINFO_SIZE+1];

#ifdef DEBUG
    printf("oligo: %08X, remainder %d, offset1 %u\n",oligo,remainder,offset1);
    printf("bitpack:\n");
    for (i = 0; i < packsize/4; i++) {
      print_vector_hex(bitpack[i]);
    }
    printf("\n");

    /* diffs[0] is used only to hold the result 0 when remainder is 0 */
    for (i = 1; i <= 8; i++) {
      print_vector(diffs[i]);
    }
    printf("end of diffs\n");
#endif  

    _diffs = (UINT4 *) &diffs;
    _diffs[3] = 0;

    *end0 = offset1 - _diffs[66-remainder];

    return offset1 - _diffs[67-remainder];
  }

#else

  debug(Bitpack64_block_offsets(offsets,bitpackptrs,offsetscomp,OFFSETSCOMP_BLOCKSIZE,oligo));


  /* Unpack all 64 diffs for non-SIMD */
  (unpacker_table[packsize])(&(diffs[1]),bitpack);

#ifdef DEBUG
  printf("oligo: %08X, remainder %d, offset0 %u, offset1 %u\n",
	 oligo,oligo % OFFSETSCOMP_BLOCKSIZE,info[1],info[METAINFO_SIZE+1]);
  printf("bitpack:\n");

  for (i = 1; i <= 64; i++) {
    printf("%d ",diffs[i]);
    if (i % 16 == 0) {
      printf("\n");
    } else if (i % 8 == 0) {
      printf("| ");
    }
  }
  printf("\n");
  printf("end of diffs\n");
#endif  

  if ((remainder = oligo % OFFSETSCOMP_BLOCKSIZE) == 0) {
    offset0 = info[1];
    *end0 = offset0 + diffs[1];
    return offset0;

  } else if (remainder < 32) {
    /* Compute necessary cumulative sums */
    /* Add 1 for start at diffs[1], and 1 to leave the first element intact */
    diffs[0] = 0;
    offset0 = info[1];

    column = (remainder - 1) % 4; /* Goes from 0 to 3 */
    /* row = ((remainder + 3) & ~3) / 4 - 1; */
    row = ((remainder + 3) >> 2) - 1; /* Assert remainder > 0 */
    for (k = column*TOTAL_ROWS + /*initial row*/0 + 1 + /*skip first diff*/ 1; k <= column*TOTAL_ROWS + row + 1; k++) {
      diffs[k] += diffs[k-1];
    };
    /* ptr = offset0 + diffs[column * TOTAL_ROWS + row + 1]; */
    ptr = offset0 + diffs[k-1];
    debug(printf("ptr remainder = %d => column %d, row %d => index %d, value %u\n",
		 remainder,column,row,column*TOTAL_ROWS + row + 1,diffs[column*TOTAL_ROWS + row + 1]));

    remainder++;
    column = (remainder - 1) % 4; /* Goes from 0 to 3 */
    /* row = ((remainder + 3) & ~3) / 4 - 1; */
    row = ((remainder + 3) >> 2) - 1;
    for (k = column*TOTAL_ROWS + /*initial row*/0 + 1 + /*skip first diff*/ 1; k <= column*TOTAL_ROWS + row + 1; k++) {
      diffs[k] += diffs[k-1];
    };
    /* *end0 = offset0 + diffs[column * TOTAL_ROWS + row + 1]; */
    *end0 = offset0 + diffs[k-1];
    debug(printf("end0 remainder = %d => column %d, row %d => index %d, value %u\n",
		 remainder,column,row,column*TOTAL_ROWS + row + 1,diffs[column*TOTAL_ROWS + row + 1]));

    return ptr;

  } else if (remainder == 63) {
    diffs[0] = 0;
    offset1 = info[METAINFO_SIZE+1];

    column = 0;
    row = 8;
    k = 10;
    ptr = offset1 - diffs[k-1];
    debug(printf("ptr remainder = %d => column %d, row %d => index %d, final k %d, diffs[%d] %u => offset %u\n",
		 remainder,column,row,column*TOTAL_ROWS + row + 1,k,column*TOTAL_ROWS + row + 1,diffs[column*TOTAL_ROWS + row + 1],ptr));

    *end0 = offset1;
    return ptr;

  } else {
    diffs[0] = 0;
    offset1 = info[METAINFO_SIZE+1];

    column = (63 - remainder) % 4; /* Goes from 0 to 3.  Assert remainder < 64 */
    /* row = ((63 - remainder) & ~3) / 4; */
    row = ((63 - remainder) >> 2) + 8;
    for (k = column*TOTAL_ROWS + /*initial row*/8 + 1 + /*skip first diff*/ 1; k <= column*TOTAL_ROWS + row + 1; k++) {
      diffs[k] += diffs[k-1];
    };
    /* ptr = offset0 + diffs[column * TOTAL_ROWS + row + 1]; */
    ptr = offset1 - diffs[k-1];
    debug(printf("ptr remainder = %d => column %d, row %d => index %d, final k %d, diffs[%d] %u => offset %u\n",
		 remainder,column,row,column*TOTAL_ROWS + row + 1,k,column*TOTAL_ROWS + row + 1,diffs[column*TOTAL_ROWS + row + 1],ptr));

    remainder++;
    column = (63 - remainder) % 4; /* Goes from 0 to 3 */
    /* row = ((63 - remainder) & ~3) / 4; */
    row = ((63 - remainder) >> 2) + 8;
    for (k = column*TOTAL_ROWS + /*initial row*/8 + 1 + /*skip first diff*/ 1; k <= column*TOTAL_ROWS + row + 1; k++) {
      diffs[k] += diffs[k-1];
    };
    /* ptr = offset0 + diffs[column * TOTAL_ROWS + row + 1]; */
    *end0 = offset1 - diffs[k-1];
    debug(printf("end0 remainder = %d => column %d, row %d => index %d, final k %d, diffs[%d] %u => offset %u\n",
		 remainder,column,row,column*TOTAL_ROWS + row + 1,k,column*TOTAL_ROWS + row + 1,diffs[column*TOTAL_ROWS + row + 1],*end0));

    return ptr;
  }

#endif

}



/* Needed for poly-T to avoid computing on metablock after the last
   one to find end0 */
Positionsptr_T
Bitpack64_offsetptr_only (Storedoligomer_T oligo) {
  UINT4 *info, nwritten;
  int packsize, remainder;
#ifdef HAVE_SSE2
  __m128i diffs[9], *bitpack;
  UINT4 *_diffs;
#else
  UINT4 diffs[65], *bitpack;
  int column, row, k;
#endif

  info = &(bitpackptrs[oligo/OFFSETSCOMP_BLOCKSIZE * METAINFO_SIZE]);

  if ((remainder = oligo % OFFSETSCOMP_BLOCKSIZE) == 0) {
    return /*offset0*/ info[1];
    
  } else {
    nwritten = info[0];
#ifdef HAVE_SSE2  
    bitpack = (__m128i *) &(offsetscomp[nwritten]);
#else
    bitpack = (UINT4 *) &(offsetscomp[nwritten]);
#endif

#ifdef ALLOW_ODD_PACKSIZES
    /* Cannot allow odd packsizes, because they destroy to one-to-one
       correspondence between nwritten and packsize */
    abort();
#else
    /* packsize 1 => 1 _mm128i => nwritten 4 UINT4s */
    /* packsize 2 => 1 _mm128i => nwritten 4 UINT4s */
    packsize = (info[2] - nwritten)/2;
#endif

#ifdef HAVE_SSE2

    if (remainder < 32) {
      /* Unpack fwd 32 cumulative sums under SIMD */
      (unpacker_table[packsize])(&(diffs[1]),bitpack);
      
      _diffs = (UINT4 *) &diffs;
      _diffs[3] = 0;
      
      return /*offset0*/ info[1] + _diffs[3+remainder];
      
    } else {
      /* Unpack rev 32 cumulative sums under SIMD */
      (unpacker_table[packsize+1])(&(diffs[1]),bitpack);
      
      _diffs = (UINT4 *) &diffs;
      _diffs[3] = 0;
      
      return /*offset1*/ info[METAINFO_SIZE+1] - _diffs[67-remainder];
    }

#else

    /* Unpack all 64 diffs for non-SIMD */
    (unpacker_table[packsize])(&(diffs[1]),bitpack);
    diffs[0] = 0;

    if (remainder < 32) {
      column = (remainder - 1) % 4; /* Goes from 0 to 3 */
      /* row = ((remainder + 3) & ~3) / 4 - 1; */
      row = ((remainder + 3) >> 2) - 1; /* Assert remainder > 0 */
      for (k = column*TOTAL_ROWS + /*initial row*/0 + 1 + /*skip first diff*/ 1; k <= column*TOTAL_ROWS + row + 1; k++) {
	diffs[k] += diffs[k-1];
      };
      /* ptr = offset0 + diffs[column * TOTAL_ROWS + row + 1]; */
      return /*offset0*/ info[1] + diffs[k-1];

    } else {
      column = (63 - remainder) % 4; /* Goes from 0 to 3.  Assert remainder < 64 */
      /* row = ((63 - remainder) & ~3) / 4; */
      row = ((63 - remainder) >> 2) + 8;
      for (k = column*TOTAL_ROWS + /*initial row*/8 + 1 + /*skip first diff*/ 1; k <= column*TOTAL_ROWS + row + 1; k++) {
	diffs[k] += diffs[k-1];
      };
      /* ptr = offset0 + diffs[column * TOTAL_ROWS + row + 1]; */
      return /*offset1*/ info[METAINFO_SIZE+1] - diffs[k-1];
    }

#endif
  }
}


/* Unpack all offsets.  Can treat offset0 as a special case */
void
Bitpack64_block_offsets (Positionsptr_T *offsets, UINT4 *bitpackptrs, Offsetscomp_T *offsetscomp,
			 Blocksize_T offsetscomp_blocksize, Storedoligomer_T oligo) {
  UINT4 *info, nwritten;
  Positionsptr_T offset0, offset1, temp;
  int packsize, k;
#ifdef HAVE_SSE2
  __m128i diffs[8], *bitpack;
  UINT4 *_diffs;
#else
  int column, row;
  UINT4 diffs[64], *bitpack, *vertical;
#endif
#ifdef DEBUG
  int i;
#endif


  info = &(bitpackptrs[oligo/OFFSETSCOMP_BLOCKSIZE * METAINFO_SIZE]);
  nwritten = info[0];
#ifdef HAVE_SSE2
  bitpack = (__m128i *) &(offsetscomp[nwritten]);
#else
  bitpack = (UINT4 *) &(offsetscomp[nwritten]);
#endif
  offset0 = info[1];
  offset1 = info[METAINFO_SIZE+1];
  packsize = (info[2] - nwritten)/2;


#ifdef HAVE_SSE2
#ifdef DEBUG
  printf("oligo: %08X, nwritten %u, offset0 %u, offset1 %u, packsize %d\n",
	 oligo,nwritten,offset0,offset1,packsize);
  printf("bitpack:\n");
  for (i = 0; i < packsize/2; i++) {
    print_vector_hex(bitpack[i]);
  }
  printf("\n");
#endif  

  _diffs = (UINT4 *) &(diffs[0]);

  /* Unpack fwd 32 cumulative sums under SIMD */
  (unpacker_table[packsize])(&(diffs[0]),bitpack);

#ifdef DEBUG
  for (i = 0; i < 8; i++) {
    print_vector(diffs[i]);
  }
  printf("end of fwd diffs\n");
#endif

  offsets[0] = offset0;
  for (k = 0; k < 32; k++) {
    offsets[k+1] = offset0 + _diffs[k];
  }

  /* Unpack rev 32 cumulative sums under SIMD */
  (unpacker_table[packsize+1])(&(diffs[0]),bitpack);

#ifdef DEBUG
  for (i = 0; i < 8; i++) {
    print_vector(diffs[i]);
  }
  printf("end of rev diffs\n");
#endif

  for (k = 0; k < 32; k++) {
    offsets[63-k] = offset1 - _diffs[k];
  }
  offsets[64] = offset1;


#ifdef DEBUG
  printf("%u\n",offsets[i]);
  for (i = 1; i <= 64; i += 4) {
    printf("%u %u %u %u\n",offsets[i],offsets[i+1],offsets[i+2],offsets[i+3]);
  }
  printf("end of offsets\n");
#endif


#else

  /* Unpack all 64 diffs for non-SIMD */
  (unpacker_table[packsize])(&(diffs[0]),bitpack);

#ifdef DEBUG
  for (i = 0; i < 64; i += 16) {
    printf("%u %u %u %u ",diffs[i],diffs[i+1],diffs[i+2],diffs[i+3]);
    printf("%u %u %u %u ",diffs[i+4],diffs[i+5],diffs[i+6],diffs[i+7]);
    printf("%u %u %u %u ",diffs[i+8],diffs[i+9],diffs[i+10],diffs[i+11]);
    printf("%u %u %u %u\n",diffs[i+12],diffs[i+13],diffs[i+14],diffs[i+15]);
  }
  printf("end of diffs horizontal (because non-SIMD unpackers are horizontal)\n");
#endif

  /* Convert to horizontal, shifting values right by 1 */
  vertical = &(diffs[0]);
  for (column = 0; column < 4; column++) {
    k = column;
    for (row = 0; row < 16; row++) {
      offsets[k+1] = *vertical++;
      k += 4;
    }
  }

#ifdef DEBUG
  printf("%u\n",offset0);
  for (i = 1; i <= 64; i += 4) {
    printf("%u %u %u %u\n",offsets[i],offsets[i+1],offsets[i+2],offsets[i+3]);
  }
  printf("end of diffs vertical\n");
#endif

  /* Perform cumulative sum */
  offsets[0] = offset0;
  offsets[1] += offset0;
  offsets[2] += offset0;
  offsets[3] += offset0;
  offsets[4] += offset0;
  for (k = 5; k <= 32; k++) {
    offsets[k] += offsets[k-4];
  }

  /* Skip offsets[33] through offsets[36] */
  for (k = 37; k <= 64; k++) {
    offsets[k] += offsets[k-4];
  }

  /* Now swap offsets */
  for (k = 33; k <= 48; k++) {
    temp = offsets[96-k];
    offsets[96-k] = offset1 - offsets[k];
    offsets[k] = offset1 - temp;
  }
  offsets[64] = offset1;


#ifdef DEBUG
  printf("%u\n",offsets[0]);
  for (i = 1; i <= 32; i += 4) {
    printf("%u %u %u %u\n",offsets[i],offsets[i+1],offsets[i+2],offsets[i+3]);
  }
  printf("\n");
  for (i = 33; i <= 64; i += 4) {
    printf("%u %u %u %u\n",offsets[i],offsets[i+1],offsets[i+2],offsets[i+3]);
  }
  printf("end of offsets\n");
#endif

#endif

  return;
}

