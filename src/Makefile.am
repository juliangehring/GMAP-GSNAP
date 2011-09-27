
if FULLDIST
bin_PROGRAMS = gmap pmap get-genome gmapindex pmapindex iit_store iit_get iit_dump oligo-count md5-compute genomeplot pdlimage
else
bin_PROGRAMS = gmap get-genome gmapindex iit_store iit_get iit_dump 
endif


GMAP_FILES = fopen.h bool.h types.h separator.h comp.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 intlist.c intlist.h listdef.h list.c list.h \
 littleendian.h bigendian.c bigendian.h \
 interval.c interval.h uintlist.c uintlist.h \
 stopwatch.c stopwatch.h access.c access.h iitdef.h iit-read.c iit-read.h \
 md5.c md5.h complement.h sequence.c sequence.h reader.c reader.h \
 genomicpos.c genomicpos.h compress.c compress.h \
 gbuffer.c gbuffer.h genome.c genome.h indexdb.c indexdb.h \
 oligo.c oligo.h block.c block.h \
 chrom.c chrom.h segmentpos.c segmentpos.h \
 chrnum.c chrnum.h chrsubset.c chrsubset.h matchdef.h match.c match.h matchpair.c matchpair.h \
 matchpool.c matchpool.h stage1.c stage1.h \
 oligoindex.c oligoindex.h \
 scores.h intron.c intron.h pairdef.h pair.c pair.h \
 pairpool.c pairpool.h stage2.c stage2.h \
 smooth.c smooth.h \
 boyer-moore.c boyer-moore.h dynprog.c dynprog.h \
 translation.c translation.h \
 changepoint.c changepoint.h stage3.c stage3.h \
 request.c request.h reqpost.c reqpost.h \
 result.c result.h blackboard.c blackboard.h \
 params.c params.h datadir.c datadir.h \
 getopt.c getopt1.c getopt.h gmap.c


# Note: dist_ commands get read by bootstrap, and don't follow the flags

# In front of stage3
#if FULLDIST
GMAP_FILES_BETA =
#else
#GMAP_FILES_BETA = maxent.c maxent.h chimera.c chimera.h 
#endif

gmap_CC = $(PTHREAD_CC)
if FULLDIST
gmap_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" -DBETATEST=1
else
gmap_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\"
endif
gmap_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG) $(PTHREAD_CFLAGS)
gmap_LDADD = $(PTHREAD_LIBS)
dist_gmap_SOURCES = $(GMAP_FILES) $(GMAP_FILES_BETA)


PMAP_FILES = fopen.h bool.h types.h separator.h comp.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 littleendian.h bigendian.c bigendian.h \
 interval.c interval.h uintlist.c uintlist.h \
 list.c list.h listdef.h \
 stopwatch.c stopwatch.h access.c access.h iitdef.h iit-read.c iit-read.h \
 genomicpos.c genomicpos.h compress.c compress.h \
 gbuffer.c gbuffer.h genome.c genome.h indexdb.c indexdb.h \
 complement.h intlist.c intlist.h md5.c md5.h sequence.c sequence.h \
 reader.c reader.h oligop.c oligop.h block.c block.h \
 chrom.c chrom.h segmentpos.c segmentpos.h \
 chrnum.c chrnum.h chrsubset.c chrsubset.h matchdef.h match.c match.h matchpair.c matchpair.h \
 matchpool.c matchpool.h stage1.c stage1.h \
 oligoindex.c oligoindex.h \
 intron.c intron.h scores.h pairdef.h pair.c pair.h \
 pairpool.c pairpool.h stage2.c stage2.h \
 smooth.c smooth.h dynprog.c dynprog.h boyer-moore.c boyer-moore.h \
 translation.c translation.h backtranslation.c backtranslation.h \
 changepoint.c changepoint.h stage3.c stage3.h \
 datadir.c datadir.h request.c request.h reqpost.c reqpost.h result.c result.h \
 blackboard.c blackboard.h params.c params.h \
 getopt.c getopt1.c getopt.h gmap.c

# In front of stage3
PMAP_FILES_BETA = maxent.c maxent.h chimera.c chimera.h 

pmap_CC = $(PTHREAD_CC)
pmap_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" -DBETATEST=1 -DPMAP=1
pmap_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG) $(PTHREAD_CFLAGS)
pmap_LDADD = $(PTHREAD_LIBS)
dist_pmap_SOURCES = 
nodist_pmap_SOURCES = $(PMAP_FILES) $(PMAP_FILES_BETA)


GET_GENOME_FILES = fopen.h bool.h types.h separator.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 intlist.c intlist.h listdef.h list.c list.h \
 littleendian.h bigendian.c bigendian.h \
 interval.c interval.h uintlist.c uintlist.h \
 stopwatch.c stopwatch.h access.c access.h iitdef.h iit-read.c iit-read.h \
 md5.c md5.h complement.h sequence.c sequence.h \
 genome.c genome.h \
 genomicpos.c genomicpos.h chrom.c chrom.h \
 chrnum.c chrnum.h chrsubset.c chrsubset.h \
 datadir.c datadir.h getopt.c getopt1.c getopt.h get-genome.c

get_genome_CC = $(PTHREAD_CC)
get_genome_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DGMAPDB=\"$(GMAPDB)\"
get_genome_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG) $(PTHREAD_CFLAGS)
get_genome_LDADD = $(PTHREAD_LIBS)
dist_get_genome_SOURCES = $(GET_GENOME_FILES)


IIT_STORE_FILES = fopen.h bool.h types.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 intlist.c intlist.h listdef.h list.c list.h \
 littleendian.h bigendian.c bigendian.h \
 interval.c interval.h uintlist.c uintlist.h \
 stopwatch.c stopwatch.h access.c access.h iitdef.h iit-write.c iit-write.h \
 tableint.c tableint.h \
 getopt.c getopt1.c getopt.h iit_store.c
iit_store_CFLAGS = $(AM_CFLAGS)  $(PTHREAD_CFLAGS)
iit_store_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG) $(PTHREAD_CFLAGS)
iit_store_LDADD =
dist_iit_store_SOURCES = $(IIT_STORE_FILES)


IIT_GET_FILES = fopen.h bool.h types.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 intlist.c intlist.h listdef.h list.c list.h \
 littleendian.h bigendian.c bigendian.h \
 interval.c interval.h uintlist.c uintlist.h \
 chrom.c chrom.h stopwatch.c stopwatch.h access.c access.h iitdef.h iit-read.c iit-read.h \
 getopt.c getopt1.c getopt.h iit_get.c
iit_get_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS)
iit_get_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG) $(PTHREAD_CFLAGS)
iit_get_LDADD = 
dist_iit_get_SOURCES = $(IIT_GET_FILES)


IIT_DUMP_FILES = fopen.h bool.h types.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 littleendian.h bigendian.c bigendian.h \
 intlist.c intlist.h list.c list.h \
 interval.c interval.h uintlist.c uintlist.h \
 chrom.c chrom.h stopwatch.c stopwatch.h access.c access.h iitdef.h iit-read.c iit-read.h \
 getopt.c getopt1.c getopt.h iit_dump.c

iit_dump_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS)
iit_dump_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG) $(PTHREAD_CFLAGS)
iit_dump_LDADD =
dist_iit_dump_SOURCES = $(IIT_DUMP_FILES)


OLIGO_COUNT_FILES = fopen.h bool.h types.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 intlist.c intlist.h listdef.h list.c list.h \
 interval.c interval.h uintlist.c uintlist.h \
 chrom.c chrom.h stopwatch.c stopwatch.h access.c access.h iitdef.h iit-read.c iit-read.h \
 md5.c md5.h complement.h sequence.c sequence.h reader.c reader.h \
 genomicpos.c genomicpos.h \
 indexdb.c indexdb.h \
 compress.c compress.h oligo.c oligo.h \
 block.c block.h \
 oligo-count.c

oligo_count_CC = $(PTHREAD_CC)
oligo_count_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS)
oligo_count_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG) $(PTHREAD_CFLAGS)
oligo_count_LDADD = $(PTHREAD_LIBS)
dist_oligo_count_SOURCES =
nodist_oligo_count_SOURCES = $(OLIGO_COUNT_FILES)


MD5_COMPUTE_FILES = fopen.h bool.h types.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 intlist.c intlist.h \
 md5.c md5.h complement.h sequence.c sequence.h \
 md5-compute.c
md5_compute_CFLAGS = $(AM_CFLAGS)
md5_compute_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
md5_compute_LDADD =
dist_md5_compute_SOURCES =
nodist_md5_compute_SOURCES = $(MD5_COMPUTE_FILES)


GENOMEPLOT_FILES = fopen.h bool.h types.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 listdef.h list.c list.h interval.c interval.h \
 uintlist.c uintlist.h \
 chrom.c chrom.h stopwatch.c stopwatch.h access.c access.h iitdef.h iit-read.c iit-read.h \
 chrsubset.c chrsubset.h color.c color.h \
 doublelist.c doublelist.h datum.c datum.h plotdata.c plotdata.h \
 plotgenes.c plotgenes.h \
 intlist.c intlist.h nr-x.c nr-x.h chrsegment.c chrsegment.h \
 chrnum.c chrnum.h datadir.c datadir.h \
 pdldata.c pdldata.h \
 getopt.c getopt1.c getopt.h genomeplot.c

genomeplot_CFLAGS = $(AM_CFLAGS)  $(PTHREAD_CFLAGS) -DGMAPDB=\"$(GMAPDB)\"
genomeplot_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG) $(PTHREAD_CFLAGS)
genomeplot_LDADD =
dist_genomeplot_SOURCES =
nodist_genomeplot_SOURCES = $(GENOMEPLOT_FILES)


PDLIMAGE_FILES = fopen.h bool.h types.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 listdef.h list.c list.h \
 color.c color.h \
 doublelist.c doublelist.h \
 intlist.c intlist.h \
 stopwatch.c stopwatch.h access.c access.h pdldata.c pdldata.h \
 getopt.c getopt1.c getopt.h pdlimage.c

pdlimage_CFLAGS = $(AM_CFLAGS)
pdlimage_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
pdlimage_LDADD =
dist_pdlimage_SOURCES =
nodist_pdlimage_SOURCES = $(PDLIMAGE_FILES)


GMAPINDEX_FILES = fopen.h bool.h types.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 intlist.c intlist.h listdef.h list.c list.h \
 littleendian.h bigendian.c bigendian.h \
 interval.c interval.h uintlist.c uintlist.h \
 stopwatch.c stopwatch.h access.c access.h iitdef.h iit-read.c iit-read.h iit-write.c iit-write.h \
 genomicpos.c genomicpos.h compress.c compress.h genome-write.c genome-write.h \
 indexdb.c indexdb.h \
 table.c table.h tableint.c tableint.h \
 chrom.c chrom.h segmentpos.c segmentpos.h \
 gmapindex.c

gmapindex_CC = $(PTHREAD_CC)
gmapindex_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS)
gmapindex_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG) $(PTHREAD_CFLAGS)
gmapindex_LDADD = $(PTHREAD_LIBS)
dist_gmapindex_SOURCES = $(GMAPINDEX_FILES)


PMAPINDEX_FILES = fopen.h bool.h types.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 littleendian.h bigendian.c bigendian.h \
 interval.c interval.h uintlist.c uintlist.h \
 list.c list.h listdef.h \
 chrom.c chrom.h stopwatch.c stopwatch.h access.c access.h iitdef.h iit-read.c iit-read.h iit-write.c iit-write.h \
 genomicpos.c genomicpos.h compress.c compress.h \
 indexdb.c indexdb.h \
 pmapindex.c

pmapindex_CC = $(PTHREAD_CC)
pmapindex_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DPMAP=1
pmapindex_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG) $(PTHREAD_CFLAGS)
pmapindex_LDADD = $(PTHREAD_LIBS)
dist_pmapindex_SOURCES = $(PMAPINDEX_FILES)
nodist_pmapindex_SOURCES = 
