

EXTRA_DIST = ss.chr17test ss.her2 \
             align.test.ok coords1.test.ok \
             setup.genomecomp.ok setup.idxpositions.ok \
             map.test.ok \
             fa.iittest iittest.iit.ok iit_get.out.ok

CLEANFILES = align.test.out \
             fa_coords coords.txt \
             gmap_setup gmap_process gmapindex chr17test.* Makefile.chr17test \
             map.test.out \
             iittest.iit iit_get.out iit_dump.out

TESTS = align.test coords1.test setup1.test map.test \
        iit_store.test iit_get.test iit_dump.test

okay:
	../src/gmap -A -g ss.chr17test ss.her2 > align.test.ok
	./fa_coords -o coords.txt ss.chr17test
	./gmap_setup -D . -d chr17test ss.chr17test && make -f Makefile.chr17test all
	../src/gmap -D . -d chr17test ss.her2 > map.test.ok
	mv -f coords.txt coords1.test.ok
	mv -f chr17test.genomecomp setup.genomecomp.ok
	mv -f chr17test.idxpositions setup.idxpositions.ok
	../src/iit_store -o iittest fa.iittest
	../src/iit_get iittest 20000 > iit_get.out.ok
	mv iittest.iit iittest.iit.ok
