# ===========================================================================
#          http://www.gnu.org/software/autoconf-archive/ax_ext.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_EXT
#
# DESCRIPTION
#
#   Find supported SIMD extensions by requesting cpuid. When an SIMD
#   extension is found, the -m"simdextensionname" is added to SIMD_FLAGS if
#   compiler supports it. For example, if "sse2" is available, then "-msse2"
#   is added to SIMD_FLAGS.
#
#   This macro calls:
#
#     AC_SUBST(SIMD_FLAGS)
#
#   And defines:
#
#     HAVE_MMX / HAVE_SSE / HAVE_SSE2 / HAVE_SSE3 / HAVE_SSSE3 / HAVE_SSE4.1 / HAVE_SSE4.2 / HAVE_AVX
#
# LICENSE
#
#   Copyright (c) 2007 Christophe Tournayre <turn3r@users.sourceforge.net>
#   Copyright (c) 2013 Michael Petch <mpetch@capp-sysware.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

#serial 13

AC_DEFUN([AX_EXT],
[
  AC_REQUIRE([AC_CANONICAL_HOST])

  case $host_cpu in
    powerpc*)
      AC_CACHE_CHECK([whether altivec is enabled and supported], [ax_cv_have_altivec_ext],
          [
            if test `/usr/sbin/sysctl -a 2>/dev/null| grep -c hw.optional.altivec` != 0; then
                if test `/usr/sbin/sysctl -n hw.optional.altivec` = 1; then
                  ax_cv_have_altivec_ext=yes
                fi
            fi
          ])

          if test "$ax_cv_have_altivec_ext" = yes; then
            AC_DEFINE(HAVE_ALTIVEC,1,[Define to 1 if you support Altivec instructions])
            AX_CHECK_COMPILE_FLAG(-faltivec, SIMD_FLAGS="$SIMD_FLAGS -faltivec", [])
          fi
    ;;


    i[[3456]]86*|x86_64*|amd64*)

      AC_REQUIRE([AX_GCC_X86_CPUID])
      AC_REQUIRE([AX_GCC_X86_AVX_XGETBV])

      AX_GCC_X86_CPUID(0x00000001)
      ecx=0
      edx=0
      if test "$ax_cv_gcc_x86_cpuid_0x00000001" != "unknown";
      then
        ecx=`echo $ax_cv_gcc_x86_cpuid_0x00000001 | cut -d ":" -f 3`
        edx=`echo $ax_cv_gcc_x86_cpuid_0x00000001 | cut -d ":" -f 4`
      fi

      AC_CACHE_CHECK([whether mmx is enabled and supported], [ax_cv_have_mmx_ext],
      [
        ax_cv_have_mmx_ext=no
        if test "$((0x$edx>>23&0x01))" = 1; then
          ax_cv_have_mmx_ext=yes
        fi
      ])

      AC_CACHE_CHECK([whether sse is enabled and supported], [ax_cv_have_sse_ext],
      [
        ax_cv_have_sse_ext=no
        if test "$((0x$edx>>25&0x01))" = 1; then
          ax_cv_have_sse_ext=yes
        fi
      ])

      AC_CACHE_CHECK([whether sse2 is enabled and supported], [ax_cv_have_sse2_ext],
      [
        ax_cv_have_sse2_ext=no
        if test "$ax_cv_want_sse2_ext" = yes; then
          if test "$((0x$edx>>26&0x01))" = 1; then
            ax_cv_have_sse2_ext=yes
          fi
        fi
      ])

      AC_CACHE_CHECK([whether sse3 is enabled and supported], [ax_cv_have_sse3_ext],
      [
        ax_cv_have_sse3_ext=no
        if test "$((0x$ecx&0x01))" = 1; then
          ax_cv_have_sse3_ext=yes
        fi
      ])

      AC_CACHE_CHECK([whether ssse3 is enabled and supported], [ax_cv_have_ssse3_ext],
      [
        ax_cv_have_ssse3_ext=no
        if test "$ax_cv_want_ssse3_ext" = yes; then
          if test "$((0x$ecx>>9&0x01))" = 1; then
            ax_cv_have_ssse3_ext=yes
          fi
        fi
      ])

      AC_CACHE_CHECK([whether sse4.1 is enabled and supported], [ax_cv_have_sse41_ext],
      [
        ax_cv_have_sse41_ext=no
        if test "$ax_cv_want_sse41_ext" = yes; then
          if test "$((0x$ecx>>19&0x01))" = 1; then
            ax_cv_have_sse41_ext=yes
          fi
        fi
      ])

      AC_CACHE_CHECK([whether sse4.2 is enabled and supported], [ax_cv_have_sse42_ext],
      [
        ax_cv_have_sse42_ext=no
        if test "$((0x$ecx>>20&0x01))" = 1; then
          ax_cv_have_sse42_ext=yes
        fi
      ])

      AC_CACHE_CHECK([whether avx is enabled and supported by processor], [ax_cv_have_avx_cpu_ext],
      [
        ax_cv_have_avx_cpu_ext=no
        if test "$((0x$ecx>>28&0x01))" = 1; then
          ax_cv_have_avx_cpu_ext=yes
        fi
      ])

      if test x"$ax_cv_have_avx_cpu_ext" = x"yes"; then
        AX_GCC_X86_AVX_XGETBV(0x00000000)

        xgetbv_eax="0"
        if test x"$ax_cv_gcc_x86_avx_xgetbv_0x00000000" != x"unknown"; then
          xgetbv_eax=`echo $ax_cv_gcc_x86_avx_xgetbv_0x00000000 | cut -d ":" -f 1`
        fi

        AC_CACHE_CHECK([whether avx is supported by operating system], [ax_cv_have_avx_ext],
        [
          ax_cv_have_avx_ext=no

          if test "$((0x$ecx>>27&0x01))" = 1; then
            if test "$((0x$xgetbv_eax&0x6))" = 6; then
              ax_cv_have_avx_ext=yes
            fi
          fi
        ])
        if test x"$ax_cv_have_avx_ext" = x"no"; then
          AC_MSG_WARN([Your processor supports AVX, but your operating system doesn't])
        fi
      fi


      AC_CACHE_CHECK([whether popcnt is enabled and supported], [ax_cv_have_popcnt_ext],
      [
        ax_cv_have_popcnt_ext=no
        if test "$((0x$ecx>>23&0x01))" = 1; then
          ax_cv_have_popcnt_ext=yes
        fi
      ])


      AX_GCC_X86_CPUID(0x80000001)
      ecx=`echo $ax_cv_gcc_x86_cpuid_0x80000001 | cut -d ":" -f 3`
      edx=`echo $ax_cv_gcc_x86_cpuid_0x80000001 | cut -d ":" -f 4`

      AC_CACHE_CHECK([whether lzcnt is enabled and supported], [ax_cv_have_lzcnt_ext],
      [
        ax_cv_have_lzcnt_ext=no
        if test "$((0x$ecx>>5&0x01))" = 1; then
          ax_cv_have_lzcnt_ext=yes
        fi
      ])


      AX_GCC_X86_CPUID(0x00000007)
      ebx=`echo $ax_cv_gcc_x86_cpuid_0x00000007 | cut -d ":" -f 2`
      ecx=`echo $ax_cv_gcc_x86_cpuid_0x00000007 | cut -d ":" -f 3`
      edx=`echo $ax_cv_gcc_x86_cpuid_0x00000007 | cut -d ":" -f 4`

      AC_CACHE_CHECK([whether avx2 is enabled and supported], [ax_cv_have_avx2_ext],
      [
        ax_cv_have_avx2_ext=no
        if test "$((0x$ebx>>5&0x01))" = 1; then
          ax_cv_have_avx2_ext=yes
        fi
      ])

      AC_CACHE_CHECK([whether bmi1 is enabled and supported], [ax_cv_have_bmi1_ext],
      [
        ax_cv_have_bmi1_ext=no
        if test "$((0x$ebx>>3&0x01))" = 1; then
          ax_cv_have_bmi1_ext=yes
        fi
      ])

      AC_CACHE_CHECK([whether bmi2 is enabled and supported], [ax_cv_have_bmi2_ext],
      [
        ax_cv_have_bmi2_ext=no
        if test "$((0x$ebx>>8&0x01))" = 1; then
          ax_cv_have_bmi2_ext=yes
        fi
      ])



      if test "$ax_cv_have_mmx_ext" = yes; then
        AX_CHECK_COMPILE_FLAG(-mmmx, ax_cv_support_mmx_ext=yes, [])
        if test x"$ax_cv_support_mmx_ext" = x"yes"; then
          SIMD_FLAGS="$SIMD_FLAGS -mmmx"
          AC_DEFINE(HAVE_MMX,1,[Define to 1 if you support mmx instructions])
        else
          AC_MSG_WARN([Your processor supports mmx instructions but not your compiler.  Can you try another compiler or update yours?])
        fi
      fi

      if test "$ax_cv_have_sse_ext" = yes; then
        AX_CHECK_COMPILE_FLAG(-msse, ax_cv_support_sse_ext=yes, [])
        if test x"$ax_cv_support_sse_ext" = x"yes"; then
          SIMD_FLAGS="$SIMD_FLAGS -msse"
          AC_DEFINE(HAVE_SSE,1,[Define to 1 if you support SSE (Streaming SIMD Extensions) instructions])
        else
          AC_MSG_WARN([Your processor supports sse instructions but not your compiler.  Can you try another compiler or update yours?])
        fi
      fi

      if test "$ax_cv_have_sse2_ext" = yes; then
	AX_CHECK_COMPILE_FLAG(-msse2, ax_cv_support_sse2_ext=yes, [])
	if test x"$ax_cv_support_sse2_ext" = x"yes"; then
	  SIMD_FLAGS="$SIMD_FLAGS -msse2"
	  AC_DEFINE(HAVE_SSE2,1,[Define to 1 if you support SSE2 (Streaming SIMD Extensions 2) instructions])
	else
	  AC_MSG_WARN([Your processor supports sse2 instructions but not your compiler.  Can you try another compiler or update yours?])
	fi
      fi

      if test "$ax_cv_have_sse3_ext" = yes; then
        AX_CHECK_COMPILE_FLAG(-msse3, ax_cv_support_sse3_ext=yes, [])
        if test x"$ax_cv_support_sse3_ext" = x"yes"; then
          SIMD_FLAGS="$SIMD_FLAGS -msse3"
          AC_DEFINE(HAVE_SSE3,1,[Define to 1 if you support SSE3 (Streaming SIMD Extensions 3) instructions])
        else
          AC_MSG_WARN([Your processor supports sse3 instructions but not your compiler.  Can you try another compiler or update yours?])
        fi
      fi

      if test "$ax_cv_have_ssse3_ext" = yes; then
        AX_CHECK_COMPILE_FLAG(-mssse3, ax_cv_support_ssse3_ext=yes, [])
        if test x"$ax_cv_support_ssse3_ext" = x"yes"; then
          SIMD_FLAGS="$SIMD_FLAGS -mssse3"
          AC_DEFINE(HAVE_SSSE3,1,[Define to 1 if you support SSSE3 (Supplemental Streaming SIMD Extensions 3) instructions])
        else
          AC_MSG_WARN([Your processor supports ssse3 instructions but not your compiler.  Can you try another compiler or update yours?])
        fi
      fi

      if test "$ax_cv_have_sse41_ext" = yes; then
	AX_CHECK_COMPILE_FLAG(-msse4.1, ax_cv_support_sse41_ext=yes, [])
	if test x"$ax_cv_support_sse41_ext" = x"yes"; then
	  SIMD_FLAGS="$SIMD_FLAGS -msse4.1"
	  AC_DEFINE(HAVE_SSE4_1,1,[Define to 1 if you support SSE4.1 (Streaming SIMD Extensions 4.1) instructions])
	else
	  AC_MSG_WARN([Your processor supports sse4.1 instructions but not your compiler.  Can you try another compiler or update yours?])
	fi
      fi

      if test "$ax_cv_have_sse42_ext" = yes; then
        AX_CHECK_COMPILE_FLAG(-msse4.2, ax_cv_support_sse42_ext=yes, [])
        if test x"$ax_cv_support_sse42_ext" = x"yes"; then
          SIMD_FLAGS="$SIMD_FLAGS -msse4.2"
          AC_DEFINE(HAVE_SSE4_2,1,[Define to 1 if you support SSE4.2 (Streaming SIMD Extensions 4.2) instructions])
        else
          AC_MSG_WARN([Your processor supports sse4.2 instructions but not your compiler.  Can you try another compiler or update yours?])
        fi
      fi

      AC_MSG_CHECKING(for immintrin.h header file)
      AC_TRY_LINK([#include <immintrin.h>],
                  [ ],
                  [ax_cv_have_immintrin_h=yes])


      if test x"$a_cv_have_immintrin_h" = x"yes"; then
        AC_MSG_RESULT([yes])

        if test "$ax_cv_have_avx_ext" = yes; then
          AX_CHECK_COMPILE_FLAG(-mavx, ax_cv_support_avx_ext=yes, [])
          if test x"$ax_cv_support_avx_ext" = x"yes"; then
            SIMD_FLAGS="$SIMD_FLAGS -mavx"
            AC_DEFINE(HAVE_AVX,1,[Define to 1 if you support AVX (Advanced Vector Extensions) instructions])
          else
            AC_MSG_WARN([Your processor supports avx instructions but not your compiler.  Can you try another compiler or update yours?])
          fi
        fi

        if test "$ax_cv_have_popcnt_ext" = yes; then
          AC_DEFINE(HAVE_POPCNT,1,[Define to 1 if you support Intel intrinsic _popcnt instruction])
        fi

        if test "$ax_cv_have_lzcnt_ext" = yes; then
          AC_DEFINE(HAVE_LZCNT,1,[Define to 1 if you support Intel intrinsic _lzcnt instruction])
        fi

        if test "$ax_cv_have_bmi1_ext" = yes; then
          AC_DEFINE(HAVE_BMI1,1,[Define to 1 if you support BMI1 (Bit Manipulation Instruction set 1)])
        fi
      else
        AC_MSG_RESULT([no])
        if test "$ax_cv_have_avx_ext" = yes; then
          AC_MSG_WARN([Your processor supports AVX but your compiler cannot find immintrin.h.  Will not use AVX.])
        fi
  
        if test "$ax_cv_have_popcnt_ext" = yes; then
          AC_MSG_WARN([Your processor supports _popcnt instructions but your compiler cannot find immintrin.h.  Will try another method.])
        fi
        if test "$ax_cv_have_lzcnt_ext" = yes; then
          AC_MSG_WARN([Your processor supports _lzcnt instructions but your compiler cannot find immintrin.h.  Will try another method.])
        fi
        if test "$ax_cv_have_bmi1_ext" = yes; then
          AC_MSG_WARN([Your processor supports bmi instructions but your compiler cannot find immintrin.h.  Will try another method.])
        fi
      fi

      if test "$ax_cv_have_bmi2_ext" = yes; then
        AC_DEFINE(HAVE_BMI2,1,[Define to 1 if you support BMI2 (Bit Manipulation Instruction set 2)])
      fi

#      if test "$ax_cv_have_avx2_ext" = yes; then
#        AX_CHECK_COMPILE_FLAG(-mavx2, ax_cv_support_avx2_ext=yes, [])
#        if test x"$ax_cv_support_avx2_ext" = x"yes"; then
#          SIMD_FLAGS="$SIMD_FLAGS -mavx2"
#          AC_DEFINE(HAVE_AVX2,1,[Define to 1 if you support AVX2 (Advanced Vector Extensions 2) instructions])
#        else
#          AC_MSG_WARN([Your processor supports avx2 instructions but not your compiler.  Can you try another compiler or update yours?])
#        fi
#      fi

  ;;
  esac

  AC_SUBST(SIMD_FLAGS)
])
