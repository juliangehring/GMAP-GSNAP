
SUBDIRS = src util tests

EXTRA_DIST = VERSION config.site

CLEANFILES = so_locations

install-data-local:
	$(mkinstalldirs) $(GMAPDB)

dist-hook:
	if test -d CVS; then \
	/usr/local/tools/bin/cvs2cl.pl; \
	fi

