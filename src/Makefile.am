
bin_PROGRAMS = gmap get-genome gmapindex iit_store iit_get iit_dump \
               gsnap snpindex cmetindex dibaseindex gsnap_tally 

GMAP_FILES = fopen.h bool.h types.h separator.h comp.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 intlistdef.h intlist.c intlist.h listdef.h list.c list.h \
 littleendian.c littleendian.h bigendian.c bigendian.h \
 interval.c interval.h uintlist.c uintlist.h \
 stopwatch.c stopwatch.h access.c access.h iitdef.h iit-read.c iit-read.h \
 md5.c md5.h complement.h sequence.c sequence.h reader.c reader.h \
 genomicpos.c genomicpos.h compress.c compress.h \
 gbuffer.c gbuffer.h genome.c genome.h indexdbdef.h indexdb.c indexdb.h indexdb_hr.c indexdb_hr.h \
 oligo.c oligo.h block.c block.h \
 chrom.c chrom.h segmentpos.c segmentpos.h \
 chrnum.c chrnum.h chrsubset.c chrsubset.h gregion.c gregion.h \
 matchdef.h match.c match.h matchpool.c matchpool.h \
 diagnostic.c diagnostic.h stage1.c stage1.h \
 intpool.c intpool.h diagdef.h diag.c diag.h diagpool.c diagpool.h \
 orderstat.c orderstat.h oligoindex.c oligoindex.h \
 scores.h intron.c intron.h maxent.c maxent.h pairdef.h pair.c pair.h \
 pairpool.c pairpool.h stage2.c stage2.h \
 smooth.c smooth.h \
 boyer-moore.c boyer-moore.h dynprog.c dynprog.h \
 translation.c translation.h \
 pbinom.c pbinom.h changepoint.c changepoint.h stage3.c stage3.h \
 request.c request.h reqpost.c reqpost.h \
 result.c result.h blackboard.c blackboard.h \
 datadir.c datadir.h \
 getopt.c getopt1.c getopt.h gmap.c


# Note: dist_ commands get read by bootstrap, and don't follow the flags

# In front of stage3
GMAP_FILES_BETA = chimera.c chimera.h 

gmap_CC = $(PTHREAD_CC)
gmap_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" -DBETATEST=1
gmap_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG) $(PTHREAD_CFLAGS)
gmap_LDADD = $(PTHREAD_LIBS)

dist_gmap_SOURCES = $(GMAP_FILES) $(GMAP_FILES_BETA)


GSNAP_FILES = fopen.h bool.h types.h separator.h comp.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 intlistdef.h intlist.c intlist.h listdef.h list.c list.h \
 littleendian.c littleendian.h bigendian.c bigendian.h \
 interval.c interval.h uintlist.c uintlist.h \
 stopwatch.c stopwatch.h access.c access.h iitdef.h iit-read.c iit-read.h \
 md5.c md5.h complement.h sequence.c sequence.h reader.c reader.h \
 genomicpos.c genomicpos.h compress.c compress.h \
 genome.c genome.h dibase.c dibase.h genome_hr.c genome_hr.h \
 indexdbdef.h indexdb.c indexdb.h indexdb_hr.c indexdb_hr.h \
 oligo.c oligo.h \
 chrom.c chrom.h segmentpos.c segmentpos.h \
 chrnum.c chrnum.h \
 maxent.c maxent.h sam.c sam.h substring.c substring.h stage3hr.c stage3hr.h \
 spanningelt.c spanningelt.h cmet.c cmet.h stage1hr.c stage1hr.h \
 request.c request.h reqpost.c reqpost.h \
 resulthr.c resulthr.h blackboard.c blackboard.h \
 datadir.c datadir.h \
 getopt.c getopt1.c getopt.h gsnap.c


# Note: dist_ commands get read by bootstrap, and don't follow the flags

gsnap_CC = $(PTHREAD_CC)
gsnap_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" -DGSNAP=1
gsnap_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG) $(PTHREAD_CFLAGS)
gsnap_LDADD = $(PTHREAD_LIBS)

dist_gsnap_SOURCES = $(GSNAP_FILES)


GMAPINDEX_FILES = fopen.h bool.h types.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 intlistdef.h intlist.c intlist.h listdef.h list.c list.h \
 littleendian.c littleendian.h bigendian.c bigendian.h \
 interval.c interval.h uintlist.c uintlist.h \
 stopwatch.c stopwatch.h access.c access.h iitdef.h iit-read.c iit-read.h iit-write.c iit-write.h \
 genomicpos.c genomicpos.h compress.c compress.h genome-write.c genome-write.h \
 indexdbdef.h indexdb.c indexdb.h \
 table.c table.h tableint.c tableint.h \
 chrom.c chrom.h segmentpos.c segmentpos.h \
 gmapindex.c

gmapindex_CC = $(PTHREAD_CC)
gmapindex_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS)
gmapindex_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG) $(PTHREAD_CFLAGS)
gmapindex_LDADD = $(PTHREAD_LIBS)
dist_gmapindex_SOURCES = $(GMAPINDEX_FILES)


SNPINDEX_FILES = fopen.h bool.h types.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 intlistdef.h intlist.c intlist.h listdef.h list.c list.h \
 littleendian.c littleendian.h bigendian.c bigendian.h \
 interval.c interval.h uintlist.c uintlist.h \
 stopwatch.c stopwatch.h access.c access.h iitdef.h iit-read.c iit-read.h \
 genomicpos.c genomicpos.h compress.c compress.h \
 indexdbdef.h indexdb.c indexdb.h \
 chrom.c chrom.h \
 complement.h md5.c md5.h sequence.c sequence.h genome.c genome.h \
 getopt.c getopt1.c getopt.h snpindex.c

snpindex_CC = $(PTHREAD_CC)
snpindex_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS)
snpindex_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG) $(PTHREAD_CFLAGS)
snpindex_LDADD = $(PTHREAD_LIBS)
dist_snpindex_SOURCES = $(SNPINDEX_FILES)


CMETINDEX_FILES = fopen.h bool.h types.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 littleendian.c littleendian.h bigendian.c bigendian.h \
 genomicpos.c genomicpos.h \
 stopwatch.c stopwatch.h access.c access.h \
 indexdbdef.h indexdb.h cmet.c cmet.h \
 getopt.c getopt1.c getopt.h cmetindex.c

cmetindex_CC = $(PTHREAD_CC)
cmetindex_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS)
cmetindex_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG) $(PTHREAD_CFLAGS)
cmetindex_LDADD = $(PTHREAD_LIBS)
dist_cmetindex_SOURCES = $(CMETINDEX_FILES)


DIBASEINDEX_FILES = fopen.h bool.h types.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 littleendian.c littleendian.h bigendian.c bigendian.h \
 interval.c interval.h uintlist.c uintlist.h \
 intlist.c intlist.h list.c list.h listdef.h \
 chrom.c chrom.h stopwatch.c stopwatch.h access.c access.h iitdef.h iit-read.c iit-read.h \
 table.c table.h iit-write.c iit-write.h \
 genomicpos.c genomicpos.h compress.c compress.h \
 indexdbdef.h indexdb_dibase.c indexdb_dibase.h \
 dibaseindex.c

dibaseindex_CC = $(PTHREAD_CC)
dibaseindex_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS)
dibaseindex_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG) $(PTHREAD_CFLAGS)
dibaseindex_LDADD = $(PTHREAD_LIBS)

dist_dibaseindex_SOURCES = $(DIBASEINDEX_FILES)



GET_GENOME_FILES = fopen.h bool.h types.h separator.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 intlistdef.h intlist.c intlist.h listdef.h list.c list.h \
 littleendian.c littleendian.h bigendian.c bigendian.h \
 interval.c interval.h uintlist.c uintlist.h \
 stopwatch.c stopwatch.h access.c access.h iitdef.h iit-read.c iit-read.h \
 md5.c md5.h complement.h sequence.c sequence.h \
 genome.c genome.h \
 genomicpos.c genomicpos.h chrom.c chrom.h \
 chrnum.c chrnum.h chrsubset.c chrsubset.h \
 datadir.c datadir.h parserange.c parserange.h \
 getopt.c getopt1.c getopt.h get-genome.c

get_genome_CC = $(PTHREAD_CC)
get_genome_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DGMAPDB=\"$(GMAPDB)\"
get_genome_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG) $(PTHREAD_CFLAGS)
get_genome_LDADD = $(PTHREAD_LIBS)
dist_get_genome_SOURCES = $(GET_GENOME_FILES)


IIT_STORE_FILES = fopen.h bool.h types.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 intlistdef.h intlist.c intlist.h listdef.h list.c list.h \
 littleendian.c littleendian.h bigendian.c bigendian.h \
 interval.c interval.h uintlist.c uintlist.h \
 stopwatch.c stopwatch.h access.c access.h iitdef.h iit-write.c iit-write.h \
 tableint.c tableint.h table.c table.h chrom.c chrom.h \
 getopt.c getopt1.c getopt.h iit_store.c
iit_store_CFLAGS = $(AM_CFLAGS)  $(PTHREAD_CFLAGS)
iit_store_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG) $(PTHREAD_CFLAGS)
iit_store_LDADD = $(PTHREAD_LIBS)
dist_iit_store_SOURCES = $(IIT_STORE_FILES)


IIT_GET_FILES = fopen.h bool.h types.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 intlistdef.h intlist.c intlist.h listdef.h list.c list.h \
 littleendian.c littleendian.h bigendian.c bigendian.h \
 interval.c interval.h uintlist.c uintlist.h \
 stopwatch.c stopwatch.h access.c access.h iitdef.h iit-read.c iit-read.h \
 complement.h parserange.c parserange.h getopt.c getopt1.c getopt.h iit_get.c
iit_get_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS)
iit_get_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG) $(PTHREAD_CFLAGS)
iit_get_LDADD = $(PTHREAD_LIBS)
dist_iit_get_SOURCES = $(IIT_GET_FILES)

IIT_DUMP_FILES = fopen.h bool.h types.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 littleendian.c littleendian.h bigendian.c bigendian.h \
 intlistdef.h intlist.c intlist.h list.c list.h \
 interval.c interval.h uintlist.c uintlist.h \
 stopwatch.c stopwatch.h access.c access.h iitdef.h iit-read.c iit-read.h \
 getopt.c getopt1.c getopt.h iit_dump.c

iit_dump_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS)
iit_dump_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG) $(PTHREAD_CFLAGS)
iit_dump_LDADD = $(PTHREAD_LIBS)
dist_iit_dump_SOURCES = $(IIT_DUMP_FILES)


GSNAP_TALLY_FILES = fopen.h bool.h types.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 intlistdef.h intlist.c intlist.h listdef.h list.c list.h \
 littleendian.c littleendian.h bigendian.c bigendian.h \
 interval.c interval.h uintlist.c uintlist.h \
 stopwatch.c stopwatch.h access.c access.h iitdef.h iit-read.c iit-read.h \
 complement.h chrom.c chrom.h genomicpos.c genomicpos.h \
 md5.c md5.h complement.h sequence.c sequence.h genome.c genome.h \
 datadir.h datadir.c \
 getopt.c getopt1.c getopt.h gsnap_tally.c

gsnap_tally_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\"
gsnap_tally_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG) $(PTHREAD_CFLAGS)
gsnap_tally_LDADD = $(PTHREAD_LIBS)
dist_gsnap_tally_SOURCES = $(GSNAP_TALLY_FILES)


