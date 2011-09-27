
if FULLDIST
bin_PROGRAMS = gmap get-genome gmapindex iit_store iit_get iit_dump oligo-count md5-compute
else
bin_PROGRAMS = gmap get-genome gmapindex iit_store iit_get iit_dump 
endif



GMAP_SOURCES = bool.h types.h separator.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 intlist.c intlist.h listdef.h list.c list.h \
 littleendian.h bigendian.c bigendian.h \
 interval.c interval.h uintlist.c uintlist.h \
 iitdef.h iit-read.c iit-read.h \
 md5.c md5.h complement.h sequence.c sequence.h reader.c reader.h \
 stopwatch.c stopwatch.h \
 genomicpos.c genomicpos.h compress.c compress.h \
 genome.c genome.h indexdb.c indexdb.h oligo.c oligo.h block.c block.h \
 chrom.c chrom.h segmentpos.c segmentpos.h \
 chrnum.c chrnum.h matchdef.h match.c match.h matchpair.c matchpair.h \
 chrsubset.c chrsubset.h stage1.c stage1.h \
 oligoindex.c oligoindex.h \
 scores.h intron.c intron.h pairdef.h pair.c pair.h \
 pairpool.c pairpool.h stage2.c stage2.h \
 smooth.c smooth.h \
 boyer-moore.c boyer-moore.h dynprog.c dynprog.h \
 mutation.c mutation.h translation.c translation.h \
 chimera.c chimera.h stage3.c stage3.h \
 request.c request.h reqpost.c reqpost.h \
 result.c result.h blackboard.c blackboard.h \
 params.c params.h datadir.c datadir.h \
 getopt.c getopt1.c getopt.h gmap.c

gmap_CC = $(PTHREAD_CC)
gmap_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DGMAPDB=\"$(GMAPDB)\"
gmap_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG) $(PTHREAD_CFLAGS)
gmap_LDADD = $(PTHREAD_LIBS)
dist_gmap_SOURCES = $(GMAP_SOURCES)


GET_GENOME_SOURCES = bool.h types.h separator.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 intlist.c intlist.h listdef.h list.c list.h \
 littleendian.h bigendian.c bigendian.h \
 interval.c interval.h uintlist.c uintlist.h \
 iitdef.h iit-read.c iit-read.h \
 md5.c md5.h complement.h sequence.c sequence.h \
 stopwatch.c stopwatch.h genome.c genome.h \
 genomicpos.c genomicpos.h chrom.c chrom.h \
 datadir.c datadir.h getopt.c getopt1.c getopt.h get-genome.c

get_genome_CFLAGS = $(AM_CFLAGS) -DGMAPDB=\"$(GMAPDB)\"
get_genome_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
get_genome_LDADD = 
dist_get_genome_SOURCES = $(GET_GENOME_SOURCES)


IIT_STORE_SOURCES = bool.h types.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 intlist.c intlist.h listdef.h list.c list.h \
 littleendian.h bigendian.c bigendian.h \
 interval.c interval.h uintlist.c uintlist.h \
 iitdef.h iit-write.c iit-write.h \
 tableint.c tableint.h iit_store.c
iit_store_CFLAGS = $(AM_CFLAGS)
iit_store_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
dist_iit_store_SOURCES = $(IIT_STORE_SOURCES)


IIT_GET_SOURCES = bool.h types.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 intlist.c intlist.h listdef.h list.c list.h \
 littleendian.h bigendian.c bigendian.h \
 interval.c interval.h uintlist.c uintlist.h \
 iitdef.h iit-read.c iit-read.h \
 iit_get.c
iit_get_CFLAGS = $(AM_CFLAGS)
iit_get_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
iit_get_LDADD = 
dist_iit_get_SOURCES = $(IIT_GET_SOURCES)


IIT_DUMP_SOURCES = bool.h types.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 littleendian.h bigendian.c bigendian.h \
 intlist.c intlist.h list.c list.h \
 interval.c interval.h uintlist.c uintlist.h \
 iitdef.h iit-read.c iit-read.h \
 iit_dump.c

iit_dump_CFLAGS = $(AM_CFLAGS)
iit_dump_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
dist_iit_dump_SOURCES = $(IIT_DUMP_SOURCES)


OLIGO_COUNT_SOURCES = bool.h types.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 intlist.c intlist.h listdef.h list.c list.h \
 interval.c interval.h uintlist.c uintlist.h \
 iitdef.h iit-read.c iit-read.h \
 md5.c md5.h complement.h sequence.c sequence.h reader.c reader.h \
 genomicpos.c genomicpos.h \
 stopwatch.c stopwatch.h \
 compress.c compress.h indexdb.c indexdb.h oligo.c oligo.h \
 block.c block.h \
 oligo-count.c
oligo_count_CFLAGS = $(AM_CFLAGS)
oligo_count_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
nodist_oligo_count_SOURCES = $(OLIGO_COUNT_SOURCES)


MD5_COMPUTE_SOURCES = bool.h types.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 intlist.c intlist.h \
 md5.c md5.h complement.h sequence.c sequence.h \
 md5-compute.c
md5_compute_CFLAGS = $(AM_CFLAGS)
md5_compute_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
nodist_md5_compute_SOURCES = $(MD5_COMPUTE_SOURCES)



GMAPINDEX_SOURCES = bool.h types.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 intlist.c intlist.h listdef.h list.c list.h \
 littleendian.h bigendian.c bigendian.h \
 interval.c interval.h uintlist.c uintlist.h \
 iitdef.h iit-read.c iit-read.h iit-write.c iit-write.h \
 genomicpos.c genomicpos.h compress.c compress.h genome-write.c genome-write.h \
 stopwatch.c stopwatch.h indexdb.c indexdb.h \
 table.c table.h tableint.c tableint.h \
 chrom.c chrom.h segmentpos.c segmentpos.h \
 gmapindex.c

gmapindex_CFLAGS = $(AM_CFLAGS)
gmapindex_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
dist_gmapindex_SOURCES = $(GMAPINDEX_SOURCES)



