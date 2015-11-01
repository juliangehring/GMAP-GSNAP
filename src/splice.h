/* $Id: splice.h 100404 2013-07-03 21:32:22Z twu $ */
#ifndef SPLICE_INCLUDED
#define SPLICE_INCLUDED
#include "bool.h"
#include "list.h"
#include "chrnum.h"
#include "genomicpos.h"
#include "compress.h"


extern List_T
Splice_solve_single (int *found_score, int *nhits, List_T hits, List_T *lowprob,

		     bool *segmenti_usedp, bool *segmentj_usedp,
		     Univcoord_T segmenti_left, Univcoord_T segmentj_left,
		     Chrnum_T segmenti_chrnum, Univcoord_T segmenti_chroffset,
		     Univcoord_T segmenti_chrhigh, Chrpos_T segmenti_chrlength,
		     Chrnum_T segmentj_chrnum, Univcoord_T segmentj_chroffset,
		     Univcoord_T segmentj_chrhigh, Chrpos_T segmentj_chrlength,
		     
		     int querylength, Compress_T query_compress,
		     int *segmenti_donor_knownpos, int *segmentj_acceptor_knownpos,
		     int *segmentj_antidonor_knownpos, int *segmenti_antiacceptor_knownpos,
		     int *segmenti_donor_knowni, int *segmentj_acceptor_knowni,
		     int *segmentj_antidonor_knowni, int *segmenti_antiacceptor_knowni,
		     int segmenti_donor_nknown, int segmentj_acceptor_nknown,
		     int segmentj_antidonor_nknown, int segmenti_antiacceptor_nknown,
		     int splicing_penalty, int max_mismatches_allowed,
		     bool first_read_p, bool plusp, int genestrand, bool subs_or_indels_p, bool sarrayp);

extern List_T
Splice_solve_double (int *found_score, int *nhits, List_T hits, List_T *lowprob,

		     bool *segmenti_usedp, bool *segmentm_usedp, bool *segmentj_usedp,
		     Univcoord_T segmenti_left, Univcoord_T segmentm_left, Univcoord_T segmentj_left,
		     Chrnum_T segmenti_chrnum, Univcoord_T segmenti_chroffset,
		     Univcoord_T segmenti_chrhigh, Chrpos_T segmenti_chrlength,
		     Chrnum_T segmentm_chrnum, Univcoord_T segmentm_chroffset,
		     Univcoord_T segmentm_chrhigh, Chrpos_T segmentm_chrlength,
		     Chrnum_T segmentj_chrnum, Univcoord_T segmentj_chroffset,
		     Univcoord_T segmentj_chrhigh, Chrpos_T segmentj_chrlength,

		     int querylength, Compress_T query_compress,
		     int *segmenti_donor_knownpos, int *segmentm_acceptor_knownpos, int *segmentm_donor_knownpos, int *segmentj_acceptor_knownpos,
		     int *segmentj_antidonor_knownpos, int *segmentm_antiacceptor_knownpos, int *segmentm_antidonor_knownpos, int *segmenti_antiacceptor_knownpos,
		     int *segmenti_donor_knowni, int *segmentm_acceptor_knowni, int *segmentm_donor_knowni, int *segmentj_acceptor_knowni,
		     int *segmentj_antidonor_knowni, int *segmentm_antiacceptor_knowni, int *segmentm_antidonor_knowni, int *segmenti_antiacceptor_knowni,
		     int segmenti_donor_nknown, int segmentm_acceptor_nknown, int segmentm_donor_nknown, int segmentj_acceptor_nknown,
		     int segmentj_antidonor_nknown, int segmentm_antiacceptor_nknown, int segmentm_antidonor_nknown, int segmenti_antiacceptor_nknown,
		     int splicing_penalty, int max_mismatches_allowed, bool plusp, int genestrand, bool subs_or_indels_p, bool sarrayp);

#endif

