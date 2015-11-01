/* $Id: popcount.h 116713 2013-11-27 19:46:01Z twu $ */
#ifndef POPCOUNT_INCLUDED
#define POPCOUNT_INCLUDED

#ifndef HAVE_BUILTIN_CTZ
extern const int mod_37_bit_position[];
#endif

#ifndef HAVE_BUILTIN_POPCOUNT
extern const int count_bits[];
#endif

#ifndef HAVE_BUILTIN_CLZ
extern const int clz_table[];
#endif

#endif

