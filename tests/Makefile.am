

EXTRA_DIST = defs ss.chr17test ss.her2 \
             align.test.ok coords1.test.ok \
             setup.genomecomp.ok setup.idxpositions.ok \
             map.test.ok \
             fa.iittest iittest.iit.ok iit_get.out.ok

CLEANFILES = align.test.out \
             fa_coords coords.chr17test \
             gmap_setup gmap_process gmapindex chr17test.* Makefile.chr17test \
             map.test.out \
             iittest.iit iit_get.out iit_dump.out

if MAINTAINER
TESTS = align.test coords1.test setup1.test iit.test \
        bioinfo.splicing.test
else
TESTS = align.test coords1.test setup1.test iit.test 
endif

TESTS_ENVIRONMENT = top_builddir=$(top_builddir)

testsubdir = testSubDir

distclean-local:
	rm -rf $(testsubdir)

okay:
	../src/gmap -A -g ss.chr17test ss.her2 > align.test.ok
	./fa_coords -o coords.chr17test ss.chr17test
	./gmap_setup -D . -d chr17test ss.chr17test && make -f Makefile.chr17test all
	../src/gmap -D . -d chr17test ss.her2 > map.test.ok
	mv -f coords.chr17test coords1.test.ok
	mv -f chr17test.genomecomp setup.genomecomp.ok
	mv -f chr17test.idxpositions setup.idxpositions.ok
	../src/iit_store -o iittest fa.iittest
	../src/iit_get iittest 20000 > iit_get.out.ok
	mv iittest.iit iittest.iit.ok
