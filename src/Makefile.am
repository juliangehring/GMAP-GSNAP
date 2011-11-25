# Copied from ../util so "gmap_setup -B src" works during "make check"
CLEANFILES = fa_coords gmap_process


bin_PROGRAMS = gmap get-genome gmapindex iit_store iit_get iit_dump \
               gsnap snpindex cmetindex atoiindex gsnap_tally 

GMAP_FILES = fopen.h bool.h types.h separator.h comp.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 intlistdef.h intlist.c intlist.h listdef.h list.c list.h \
 littleendian.c littleendian.h bigendian.c bigendian.h \
 interval.c interval.h uintlist.c uintlist.h \
 stopwatch.c stopwatch.h access.c access.h iitdef.h iit-read.c iit-read.h \
 md5.c md5.h complement.h sequence.c sequence.h reader.c reader.h \
 genomicpos.c genomicpos.h compress.c compress.h \
 gbuffer.c gbuffer.h genome.c genome.h genome_hr.c genome_hr.h \
 indexdbdef.h indexdb.c indexdb.h indexdb_hr.c indexdb_hr.h \
 oligo.c oligo.h block.c block.h \
 chrom.c chrom.h segmentpos.c segmentpos.h \
 chrnum.c chrnum.h chrsubset.c chrsubset.h uinttable.c uinttable.h gregion.c gregion.h \
 matchdef.h match.c match.h matchpool.c matchpool.h \
 diagnostic.c diagnostic.h stage1.c stage1.h \
 diagdef.h diag.c diag.h diagpool.c diagpool.h \
 orderstat.c orderstat.h oligoindex.c oligoindex.h oligoindex_hr.c oligoindex_hr.h \
 scores.h intron.c intron.h maxent.c maxent.h maxent_hr.c maxent_hr.h samflags.h pairdef.h pair.c pair.h \
 pairpool.c pairpool.h stage2.c stage2.h \
 smooth.c smooth.h \
 splicetrie_build.c splicetrie_build.h splicetrie.c splicetrie.h \
 boyer-moore.c boyer-moore.h dynprog.c dynprog.h \
 translation.c translation.h \
 pbinom.c pbinom.h changepoint.c changepoint.h sense.h stage3.c stage3.h \
 request.c request.h result.c result.h \
 inbuffer.c inbuffer.h outbuffer.c outbuffer.h \
 chimera.c chimera.h datadir.c datadir.h \
 getopt.c getopt1.c getopt.h gmap.c


# Note: dist_ commands get read by bootstrap, and don't follow the flags

gmap_CC = $(PTHREAD_CC)
gmap_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) $(POPCNT_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\"
gmap_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG) $(PTHREAD_CFLAGS)
gmap_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS)

dist_gmap_SOURCES = $(GMAP_FILES)


GSNAP_FILES = fopen.h bool.h types.h separator.h comp.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 intlistdef.h intlist.c intlist.h listdef.h list.c list.h \
 littleendian.c littleendian.h bigendian.c bigendian.h \
 interval.c interval.h uintlist.c uintlist.h \
 stopwatch.c stopwatch.h access.c access.h iitdef.h iit-read.c iit-read.h \
 md5.c md5.h complement.h sequence.c sequence.h reader.c reader.h \
 genomicpos.c genomicpos.h compress.c compress.h \
 genome.c genome.h genome_hr.c genome_hr.h \
 indexdbdef.h indexdb.c indexdb.h indexdb_hr.c indexdb_hr.h \
 oligo.c oligo.h \
 chrom.c chrom.h segmentpos.c segmentpos.h \
 chrnum.c chrnum.h \
 maxent_hr.c maxent_hr.h samflags.h samprint.c samprint.h \
 mapq.c mapq.h shortread.c shortread.h substring.c substring.h stage3hr.c stage3hr.h \
 goby.c goby.h \
 spanningelt.c spanningelt.h cmet.c cmet.h atoi.c atoi.h \
 comp.h maxent.c maxent.h pairdef.h pair.c pair.h pairpool.c pairpool.h diag.c diag.h diagpool.c diagpool.h \
 orderstat.c orderstat.h oligoindex.c oligoindex.h oligoindex_hr.c oligoindex_hr.h stage2.c stage2.h \
 intron.c intron.h boyer-moore.c boyer-moore.h changepoint.c changepoint.h pbinom.c pbinom.h \
 dynprog.c dynprog.h gbuffer.c gbuffer.h translation.c translation.h smooth.c smooth.h \
 chimera.c chimera.h sense.h stage3.c stage3.h \
 splicetrie_build.c splicetrie_build.h splicetrie.c splicetrie.h \
 stage1hr.c stage1hr.h \
 request.c request.h resulthr.c resulthr.h \
 inbuffer.c inbuffer.h outbuffer.c outbuffer.h \
 datadir.c datadir.h mode.h \
 getopt.c getopt1.c getopt.h gsnap.c


# Note: dist_ commands get read by bootstrap, and don't follow the flags

gsnap_CC = $(PTHREAD_CC)
gsnap_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) $(POPCNT_CFLAGS) $(GOBY_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" -DMAX_READLENGTH=$(MAX_READLENGTH) -DGSNAP=1
gsnap_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG) $(PTHREAD_CFLAGS) $(GOBY_LDFLAGS)
gsnap_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(GOBY_LIBS)

dist_gsnap_SOURCES = $(GSNAP_FILES)


GMAPINDEX_FILES = fopen.h bool.h types.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 intlistdef.h intlist.c intlist.h listdef.h list.c list.h \
 littleendian.c littleendian.h bigendian.c bigendian.h \
 interval.c interval.h uintlist.c uintlist.h \
 stopwatch.c stopwatch.h access.c access.h iitdef.h iit-read.c iit-read.h iit-write.c iit-write.h \
 genomicpos.c genomicpos.h compress.c compress.h genome-write.c genome-write.h \
 genome_hr.c genome_hr.h indexdbdef.h indexdb.c indexdb.h \
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
 genome_hr.c genome_hr.h indexdbdef.h indexdb.c indexdb.h \
 chrom.c chrom.h \
 complement.h md5.c md5.h sequence.c sequence.h genome.c genome.h \
 datadir.c datadir.h \
 getopt.c getopt1.c getopt.h snpindex.c

snpindex_CC = $(PTHREAD_CC)
snpindex_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DGMAPDB=\"$(GMAPDB)\"
snpindex_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG) $(PTHREAD_CFLAGS)
snpindex_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS)
dist_snpindex_SOURCES = $(SNPINDEX_FILES)


CMETINDEX_FILES = fopen.h bool.h types.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 littleendian.c littleendian.h bigendian.c bigendian.h \
 genomicpos.c genomicpos.h \
 stopwatch.c stopwatch.h access.c access.h \
 interval.h interval.c iitdef.h iit-read.h iit-read.c complement.h compress.c compress.h \
 genome_hr.c genome_hr.h indexdbdef.h indexdb.h indexdb.c cmet.c cmet.h \
 list.c list.h listdef.h datadir.c datadir.h \
 getopt.c getopt1.c getopt.h cmetindex.c

cmetindex_CC = $(PTHREAD_CC)
cmetindex_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DGMAPDB=\"$(GMAPDB)\"
cmetindex_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG) $(PTHREAD_CFLAGS)
cmetindex_LDADD = $(PTHREAD_LIBS)
dist_cmetindex_SOURCES = $(CMETINDEX_FILES)


ATOIINDEX_FILES = fopen.h bool.h types.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 littleendian.c littleendian.h bigendian.c bigendian.h \
 genomicpos.c genomicpos.h \
 stopwatch.c stopwatch.h access.c access.h \
 interval.h interval.c iitdef.h iit-read.h iit-read.c complement.h compress.c compress.h \
 genome_hr.c genome_hr.h indexdbdef.h indexdb.h indexdb.c atoi.c atoi.h \
 list.c list.h listdef.h datadir.c datadir.h \
 getopt.c getopt1.c getopt.h atoiindex.c

atoiindex_CC = $(PTHREAD_CC)
atoiindex_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DGMAPDB=\"$(GMAPDB)\"
atoiindex_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG) $(PTHREAD_CFLAGS)
atoiindex_LDADD = $(PTHREAD_LIBS)
dist_atoiindex_SOURCES = $(ATOIINDEX_FILES)


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
get_genome_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS)
dist_get_genome_SOURCES = $(GET_GENOME_FILES)


IIT_STORE_FILES = fopen.h bool.h types.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 intlistdef.h intlist.c intlist.h listdef.h list.c list.h \
 littleendian.c littleendian.h bigendian.c bigendian.h \
 interval.c interval.h uintlist.c uintlist.h \
 stopwatch.c stopwatch.h access.c access.h iitdef.h iit-write.c iit-write.h \
 tableint.c tableint.h table.c table.h chrom.c chrom.h \
 getopt.c getopt1.c getopt.h iit_store.c

iit_store_CC = $(PTHREAD_CC)
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

iit_get_CC = $(PTHREAD_CC)
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

iit_dump_CC = $(PTHREAD_CC)
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
 gsnapread.c gsnapread.h \
 datadir.h datadir.c \
 getopt.c getopt1.c getopt.h gsnap_tally.c

gsnap_tally_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\"
gsnap_tally_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG) $(PTHREAD_CFLAGS)
gsnap_tally_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS)
dist_gsnap_tally_SOURCES = $(GSNAP_TALLY_FILES)


