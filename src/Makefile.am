
if FULLDIST
bin_PROGRAMS = gmap get-genome gmapindex iit_store iit_get iit_dump oligo-count md5-compute genomeplot
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
 stage3.c stage3.h \
 request.c request.h reqpost.c reqpost.h \
 result.c result.h blackboard.c blackboard.h \
 params.c params.h datadir.c datadir.h \
 getopt.c getopt1.c getopt.h gmap.c

# In front of stage3
GMAP_SOURCES_BETA = nmath.c nmath.h maxent.c maxent.h chimera.c chimera.h 

gmap_CC = $(PTHREAD_CC)
if FULLDIST
gmap_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DGMAPDB=\"$(GMAPDB)\" -DBETATEST=1
else
gmap_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DGMAPDB=\"$(GMAPDB)\"
endif
gmap_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG) $(PTHREAD_CFLAGS)
gmap_LDADD = $(PTHREAD_LIBS)
if FULLDIST
dist_gmap_SOURCES = $(GMAP_SOURCES) $(GMAP_SOURCES_BETA)
else
dist_gmap_SOURCES = $(GMAP_SOURCES)
endif


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

get_genome_CC = $(PTHREAD_CC)
get_genome_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DGMAPDB=\"$(GMAPDB)\"
get_genome_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG) $(PTHREAD_CFLAGS)
get_genome_LDADD = $(PTHREAD_LIBS)
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
iit_store_LDADD =
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
iit_dump_LDADD =
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

oligo_count_CC = $(PTHREAD_CC)
oligo_count_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS)
oligo_count_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG) $(PTHREAD_CFLAGS)
oligo_count_LDADD = $(PTHREAD_LIBS)
if FULLDIST
dist_oligo_count_SOURCES = $(OLIGO_COUNT_SOURCES)
else
nodist_oligo_count_SOURCES = $(OLIGO_COUNT_SOURCES)
endif


MD5_COMPUTE_SOURCES = bool.h types.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 intlist.c intlist.h \
 md5.c md5.h complement.h sequence.c sequence.h \
 md5-compute.c
md5_compute_CFLAGS = $(AM_CFLAGS)
md5_compute_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
md5_compute_LDADD =
if FULLDIST
dist_md5_compute_SOURCES = $(MD5_COMPUTE_SOURCES)
else
nodist_md5_compute_SOURCES = $(MD5_COMPUTE_SOURCES)
endif


GENOMEPLOT_SOURCES = bool.h types.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 listdef.h list.c list.h interval.c interval.h \
 uintlist.c uintlist.h iitdef.h iit-read.c iit-read.h \
 chrsubset.c chrsubset.h color.c color.h \
 doublelist.c doublelist.h datum.c datum.h plotdata.c plotdata.h \
 intlist.c intlist.h nr-x.c nr-x.h chrsegment.c \
 datadir.c datadir.h \
 getopt.c getopt1.c getopt.h genomeplot.c

genomeplot_CFLAGS = $(AM_CFLAGS) -DGMAPDB=\"$(GMAPDB)\"
genomeplot_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
genomeplot_LDADD =
if FULLDIST
dist_genomeplot_SOURCES = $(GENOMEPLOT_SOURCES)
else
nodist_genomeplot_SOURCES = $(GENOMEPLOT_SOURCES)
endif


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

gmapindex_CC = $(PTHREAD_CC)
gmapindex_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS)
gmapindex_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG) $(PTHREAD_CFLAGS)
gmapindex_LDADD = $(PTHREAD_LIBS)
dist_gmapindex_SOURCES = $(GMAPINDEX_SOURCES)

