# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.57)

## We want to get the version number from file 'VERSION' (rather than
## specifying the version info in 'configure.ac'.  Hence, we need a bit
## of M4 magic.  Note that M4 esyscmd has a trailing newline because the
## shell command output has one, hence the patsubst() trick.
#m4_define([PKG_VERSION],
#          m4_bpatsubst(m4_esyscmd([cat VERSION]),
#                       [\([0-9.]*\)\(\w\|\W\)*],
#                       [\1]))
m4_define([PKG_VERSION],
          m4_bpatsubst(m4_esyscmd([cat VERSION]),
                       [\([0-9.\-]*\)\(\w\|\W\)*],
                       [\1]))
AC_INIT([gmap],PKG_VERSION,[Thomas Wu <twu@gene.com>])
AC_CONFIG_SRCDIR([src/gmap.c])
AC_CONFIG_HEADER([src/config.h])
AC_CONFIG_AUX_DIR([config])
# For autotest: AC_CONFIG_TESTDIR([tests],[tests:src:util])

AC_CANONICAL_SYSTEM
AC_ARG_PROGRAM

#AM_INIT_AUTOMAKE([no-dependencies])
AM_INIT_AUTOMAKE(AC_PACKAGE_NAME, AC_PACKAGE_VERSION)

AM_CONDITIONAL(FULLDIST,test "$enable_fulldist" = yes)
AC_ARG_ENABLE([fulldist],
	       AC_HELP_STRING([--enable-fulldist],
                              [For use by program maintainer]),
	      [enable_fulldist="yes"],
              [enable_fulldist="no"])

### Read defaults
if test -z "$CONFIG_SITE"; then
  CONFIG_SITE=""
fi
CONFIG_SITE="./config.site $CONFIG_SITE"
for ac_site_file in $CONFIG_SITE; do
  if test -r "$ac_site_file"; then
    echo "loading script $ac_site_file"
    . "$ac_site_file"
  fi
done

# For Perl scripts, because otherwise configure inserts 
# @bindir@ => "${exec_prefix}/bin";
AC_MSG_CHECKING(bindir)
ACX_EXPAND(BINDIR,$bindir)
AC_SUBST(BINDIR)
AC_MSG_RESULT($BINDIR)

# Genome directory
AC_MSG_CHECKING(gmapdb)
AC_ARG_WITH([gmapdb],
            AC_HELP_STRING([--with-gmapdb=DIR],
                           [default GMAP database directory]),
            [answer="$withval"],
            [answer=""])
if test x"$answer" == x; then
	ACX_EXPAND(GMAPDB,'${prefix}/share')
else
	GMAPDB=$answer
fi
AC_SUBST(GMAPDB)
AC_MSG_RESULT($GMAPDB)


# Works in conjunction with AC_PROG_LIBTOOL -- Commented out because no libraries being built
#AC_MSG_CHECKING(whether to link statically)
#AC_ARG_ENABLE([static-linking],
#               AC_HELP_STRING([--enable-static-linking],
#                              [Link binaries statically (default=no)]),
#              [answer="$enableval"],
#              [answer=""])
#case x"$answer" in
#	xyes)
#	AC_MSG_RESULT(enabled)
#	STATIC_LDFLAG="-all-static"
#	;;
#	
#	xno)
#	AC_MSG_RESULT(disabled)
#	;;
#
#	x)
#	AC_MSG_RESULT(not specified so disabled by default)
#	;;
#esac	
#AC_SUBST(STATIC_LDFLAG)


# Checks for programs.
#AC_PATH_PROG([PERL],[perl])
ACX_PATH_PERL
AC_PROG_CC
AC_PROG_LIBTOOL

# Checks for libraries.
AC_CHECK_LIB(m, rint)

#AC_CHECK_LIB(popt, poptGetContext, [answer="yes"], [answer="no"])
#if test x"$answer" == xyes; then
#	AC_DEFINE([HAVE_LIBPOPT],[1],[Define to 1 if you have the 'popt' library (-lpopt)])
#	POPT_LIBS="-lpopt"
#else
#	POPT_LIBS=""
#fi
#AC_SUBST(POPT_LIBS)

# In call to ACX_PTHREAD, don't provide ACTION-IF-FOUND; otherwise,
# HAVE_PTHREAD won't be defined
AC_MSG_CHECKING(for pthreads feature)
AC_ARG_ENABLE([pthreads],
              AC_HELP_STRING([--enable-pthreads],
                             [Enable pthreads (default=yes)]),
              [answer="$enableval"],
              [answer=""])
case x"$answer" in
	xyes)
        AC_MSG_RESULT(enabled)
	ACX_PTHREAD([],[pthread_warning=yes])
	;;
	
	xno)
	AC_MSG_RESULT(disabled)
	;;

	x)
	AC_MSG_RESULT([not specified so enabled by default])
	ACX_PTHREAD([],[pthread_warning=yes])
	;;
esac	

# Checks for header files.
AC_HEADER_STDC
AC_HEADER_DIRENT
AC_CHECK_HEADERS([fcntl.h limits.h stddef.h stdlib.h string.h strings.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_BIGENDIAN
AC_C_CONST
AC_C_VOLATILE
AC_TYPE_SIZE_T
AC_TYPE_OFF_T
AC_CHECK_SIZEOF(unsigned long)

# Checks for library functions.
#AC_FUNC_MALLOC -- Causes use of rpl_malloc
#AC_FUNC_REALLOC -- Causes use of rpl_realloc
#AC_FUNC_CLOSEDIR_VOID -- We are checking return value

AC_FUNC_MMAP
ACX_MMAP_FLAGS

AC_CHECK_FUNCS([ceil floor log madvise memcpy memmove memset munmap pow rint strtoul sysconf sysctl])

ACX_PAGESIZE

AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([src/Makefile])
AC_CONFIG_FILES([util/Makefile])
AC_CONFIG_FILES([util/gmap_compress.pl])
AC_CONFIG_FILES([util/gmap_uncompress.pl])
AC_CONFIG_FILES([util/gmap_setup.pl])
AC_CONFIG_FILES([util/md_coords.pl])
AC_CONFIG_FILES([util/fa_coords.pl])
AC_CONFIG_FILES([tests/Makefile])
AC_CONFIG_FILES([tests/align.test],[chmod +x tests/align.test])
AC_CONFIG_FILES([tests/coords1.test],[chmod +x tests/coords1.test])
AC_CONFIG_FILES([tests/setup1.test],[chmod +x tests/setup1.test])
AC_CONFIG_FILES([tests/setup2.test],[chmod +x tests/setup2.test])
AC_CONFIG_FILES([tests/map.test],[chmod +x tests/map.test])
AC_CONFIG_FILES([tests/iit_store.test],[chmod +x tests/iit_store.test])
AC_CONFIG_FILES([tests/iit_get.test],[chmod +x tests/iit_get.test])
AC_CONFIG_FILES([tests/iit_dump.test],[chmod +x tests/iit_dump.test])

AC_OUTPUT

if test x"$perl_warning" = xyes; then
  AC_MSG_WARN([

WARNING: You don't seem to have perl installed with working
versions of the IO::File and Getopt::Std modules.  If you know
where such a version of Perl exists, please specify its full
path under PERL in the config.site file and run configure again.
])
fi

if test x"$pthread_warning" = xyes; then
  AC_MSG_WARN([

WARNING: You don't seem to have headers and libraries for pthreads.
GMAP will work without these, but the program can run faster on
multiprocessor machines if POSIX threads are available.  If you wish
to enable pthreads, please specify the appropriate values for CC,
PTHREAD_CFLAGS, and PTHREAD_LIBS in the config.site file and run
configure again.  Otherwise, you may proceed with the rest of the
installation.
])

fi

