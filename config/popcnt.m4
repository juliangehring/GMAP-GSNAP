
AC_DEFUN([ACX_POPCNT], [
  AC_REQUIRE([AC_CANONICAL_HOST])
  AC_LANG_SAVE
  AC_LANG_C

  POPCNT_CFLAGS=""
  acx_popcnt_ok=no

  save_CFLAGS="$CFLAGS"
  CFLAGS="$CFLAGS -mpopcnt"

  AC_MSG_CHECKING(whether -mpopcnt compiler flag works)
  AC_RUN_IFELSE(
  [AC_LANG_PROGRAM([[#include <stdio.h>
#include <stdlib.h>]],
                   [[unsigned int x = rand();
printf("%08X ",x);
#ifdef HAVE_BUILTIN_CLZ
printf("clz=%d ",__builtin_clz(x));
#endif
#ifdef HAVE_BUILTIN_CTZ
printf("ctz=%d ",__builtin_ctz(x));
#endif
#ifdef HAVE_BUILTIN_POPCOUNT
printf("popcount=%d ",__builtin_popcount(x));
#endif
]])],
  [acx_popcnt_ok=yes])

  CFLAGS="$save_CFLAGS"

  AC_MSG_RESULT($acx_popcnt_ok)
  if test "x$acx_popcnt_ok" = xyes; then
    POPCNT_CFLAGS="-mpopcnt"
  fi

  AC_SUBST(POPCNT_CFLAGS)

AC_LANG_RESTORE
])dnl

