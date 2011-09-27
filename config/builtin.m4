
AC_DEFUN([ACX_BUILTIN], [
AC_LANG_SAVE
AC_LANG(C)

AC_MSG_CHECKING(for __builtin_popcount)
AC_RUN_IFELSE(
  [AC_LANG_PROGRAM([[]],
                   [[return (__builtin_popcount(0xffffffffu) == 32) ? 0 : 9;]])],
  [AC_MSG_RESULT(yes)
   AC_DEFINE([HAVE_BUILTIN_POPCOUNT],[1],[Define to 1 if __builtin_popcount works.])],
  [AC_MSG_RESULT(no)])

AC_MSG_CHECKING(for __builtin_clz)
AC_RUN_IFELSE(
  [AC_LANG_PROGRAM([[]],
                   [[return (__builtin_clz(0x1u) == 31) ? 0 : 9;]])],
  [AC_MSG_RESULT(yes)
   AC_DEFINE([HAVE_BUILTIN_CLZ],[1],[Define to 1 if __builtin_clz works.])],
  [AC_MSG_RESULT(no)])

AC_MSG_CHECKING(for __builtin_ctz)
AC_RUN_IFELSE(
  [AC_LANG_PROGRAM([[]],
                   [[return (__builtin_ctz(0x80000000u) == 31) ? 0 : 9;]])],
  [AC_MSG_RESULT(yes)
   AC_DEFINE([HAVE_BUILTIN_CTZ],[1],[Define to 1 if __builtin_ctz works.])],
  [AC_MSG_RESULT(no)])

AC_LANG_RESTORE
])


