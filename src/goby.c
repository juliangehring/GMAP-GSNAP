static char rcsid[] = "$Id: goby.c 63197 2012-05-03 17:41:52Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "goby.h"

#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>

#include "assert.h"
#include "mem.h"
#include "chrnum.h"
#include "substring.h"
#include "samflags.h"
#include "samprint.h"


/* #define DEBUG */
#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


static bool show_refdiff_p;


#ifdef HAVE_GOBY
#include <goby/C_Reads.h>
#include <goby/C_Alignments.h>
#include <goby/C_CompactHelpers.h>

struct Gobyreader_T {
  CReadsHelper *helper;
  bool complement_reads_p;
};

struct Gobywriter_T {
  CAlignmentsWriterHelper *helper;
};

#endif


static char complCode[129] = "???????????????????????????????? ??#$%&')(*+,-./0123456789:;>=<??TVGHEFCDIJMLKNOPQYSAABWXRZ]?[^_`tvghefcdijmlknopqysaabwxrz}|{~?";


/**
 * Duplicate a string with an optional length to copy. If the incoming
 * str is NULL, this will return null. If length is -1, the length
 * that is copied is strlen(str). The size of the returned buffer is
 * always length + 1 (to include the trailing '\0'). The caller is required
 * to FREE the string.
 * @param str the string to copy
 * @param the maximum length to copy or -1 for the whole string
 * @return the duplicate string.
 */
static char *
copy_string(char *str, int length) {
  int copy_length = length;
  char *new_str = (char *) NULL;

  if (str != NULL) {
    if (copy_length == -1) {
      copy_length = strlen(str);
    }
    new_str = (char *) CALLOC(copy_length + 1, sizeof(char));
    strncpy(new_str, str, copy_length);
    new_str[copy_length] = '\0';
  }

  return new_str;
}

void
Goby_setup (bool show_refdiff_p_in) {
  show_refdiff_p = show_refdiff_p_in;
  return;
}

void
Goby_shutdown () {
#ifdef HAVE_GOBY
  goby_shutdownProtobuf();
#endif
  return;
}


/************************************************************************
 *   Reader
 ************************************************************************/


Gobyreader_T
Goby_reader_new (char **files, int nfiles, unsigned long window_start, unsigned long window_end, bool complement_reads_p) {
#ifdef HAVE_GOBY
  Gobyreader_T new = (Gobyreader_T) MALLOC(sizeof(*new));

  new->complement_reads_p = complement_reads_p;
  fprintf(stderr,"Opening %s start=%lu, end=%lu\n",files[0], window_start, window_end);
  gobyReads_openReadsReaderWindowed(files,nfiles,/*circularp*/false,window_start,window_end,&new->helper);
  gobyReads_avoidZeroQuals(new->helper, 1);
  return new;
#else
  return NULL;
#endif
}


void
Goby_reader_finish (Gobyreader_T reader) {
#ifdef HAVE_GOBY
  gobyReads_finished(reader->helper);
#endif
  return;
}


void
Goby_reader_free (Gobyreader_T *old) {
  FREE(*old);
  return;
}


Shortread_T
Goby_read (Shortread_T *queryseq2, Gobyreader_T reader, int barcode_length,
	   bool invert_first_p, bool invert_second_p) {
#ifdef HAVE_GOBY
  unsigned long goby_read_index;
  char *acc, *read_identifier = NULL, *description = NULL;
  char *sequence1, *quality1, *sequence2, *quality2;
  int sequence1_length, quality1_length, sequence2_length, quality2_length, acc_length;
  int i;

  sequence1_length = 0;
  while (sequence1_length == 0) {
    /* Ignore empty sequences */
    if (gobyReads_hasNext(reader->helper) != 1) {
      return (Shortread_T) NULL;
    }
    goby_read_index = 
      gobyReads_nextSequencePair(reader->helper,&read_identifier,&description,
				 &sequence1,&sequence1_length,
				 &quality1,&quality1_length,
				 &sequence2,&sequence2_length,
				 &quality2,&quality2_length);
    if (sequence1_length != 0) {
      acc = (char *) CALLOC(25,sizeof(char));
      sprintf(acc, "%lu", goby_read_index);
      description = copy_string(description, -1);
    }
  }

  if (reader->complement_reads_p) {
    debug(
	  if (sequence1_length > 0) {
	    fprintf(stderr,"before complement, sequence1: %s\n", sequence1);
	  }
	  if (sequence2_length > 0) {
	    fprintf(stderr,"before complement, sequence2: %s\n", sequence2);
	  }
	  );
    for (i = 0; i < sequence1_length; i++) {
      sequence1[i] = complCode[(int) sequence1[i]];
    }
    for (i = 0; i < sequence2_length; i++) {
      sequence2[i] = complCode[(int) sequence2[i]];
    }
    debug(
	  if (sequence1_length > 0) {
	    fprintf(stderr," after complement, sequence1: %s\n", sequence1);
	  }
	  if (sequence2_length > 0) {
	    fprintf(stderr," after complement, sequence2: %s\n", sequence2);
	  }
	  );
  }

  *queryseq2 = Shortread_new(/*acc*/NULL,/*description*/NULL,/*filterp*/false,
			     sequence2,sequence2_length,quality2,quality2_length,
			     barcode_length,invert_second_p, /*copy_acc*/false);

  return Shortread_new(acc,description,/*filterp*/false,
		       sequence1,sequence1_length,quality1,quality1_length,
		       barcode_length,invert_first_p,/*copy_acc*/false);
#else
  return (Shortread_T) NULL;
#endif
}


/************************************************************************
 *   Writer
 ************************************************************************/


Gobywriter_T
Goby_writer_new (char *output_root, char *aligner_name, char *aligner_version) {
#ifdef HAVE_GOBY
  Gobywriter_T new = (Gobywriter_T) MALLOC(sizeof(*new));

  gobyAlignments_openAlignmentsWriterDefaultEntriesPerChunk(output_root,&new->helper);
  gobyAlignments_setAlignerName(new->helper,aligner_name);
  gobyAlignments_setAlignerVersion(new->helper,aligner_version);
  gobyGsnap_startAlignment(new->helper);
  gobyCapture_open(new->helper, 1);
  return new;
#else
  return NULL;
#endif
}


void
Goby_writer_finish (Gobywriter_T writer, Gobyreader_T reader) {
#ifdef HAVE_GOBY
  gobyAlignments_finished(writer->helper,reader->helper->numberOfReads);
#endif
  return;
}


void
Goby_writer_free (Gobywriter_T *old) {
#ifdef HAVE_GOBY
  gobyCapture_close((*old)->helper);
#endif
  FREE(*old);
  return;
}


void
Goby_writer_add_chromosomes (Gobywriter_T writer, IIT_T chromosome_iit) {
#ifdef HAVE_GOBY
  int nintervals, gsnap_target_index;
  char *gsnap_target_label;
  bool allocp;
  Interval_T interval;
  unsigned int length;

  nintervals = IIT_total_nintervals(chromosome_iit);
  for (gsnap_target_index = 1; gsnap_target_index <= nintervals; gsnap_target_index++) {
    gsnap_target_label = IIT_label(chromosome_iit,gsnap_target_index,&allocp);
    interval = IIT_interval(chromosome_iit,gsnap_target_index);
    length = Interval_length(interval);
    /* goby_target_index is 0-based, gsnap_target_index is 1-based. */
    gobyAlignments_addTarget(writer->helper,gsnap_target_index - 1,gsnap_target_label,length);
    debug(fprintf(stderr, "%u is %s\n", gsnap_target_index - 1, gsnap_target_label));
    if (allocp == true) {
      FREE(gsnap_target_label);
    }
  }
#endif
  return;
}


void
Goby_file_handles (FILE **fp_capture, FILE **fp_ignore, Gobywriter_T writer) {
#ifdef HAVE_GOBY
  *fp_capture = gobyCapture_fileHandle(writer->helper);
  *fp_ignore = gobyCapture_ignoredFileHandle(writer->helper);
#else
  *fp_capture = NULL;
  *fp_ignore = NULL;
#endif
  return;
}


void
Goby_observe_aligned(Gobywriter_T writer) {
#ifdef HAVE_GOBY
  writer->helper->numberOfAlignedReads++;
#endif /* HAVE_GOBY */
  return;
}


void
Goby_start_capture (Gobywriter_T writer) {
#ifdef HAVE_GOBY
  gobyCapture_startNew(writer->helper);
#endif
}


void
Goby_finish_capture (Gobywriter_T writer) {
#ifdef HAVE_GOBY
  char *capturedData;
  char *ignoredData;
  gobyCapture_flush(writer->helper);
  capturedData = gobyCapture_capturedData(writer->helper);
  if (strlen(capturedData) > 0) {
      Goby_observe_aligned(writer);
      gobyGsnap_parse(writer->helper, capturedData);
  }
#endif
  return;
}


void
Goby_print_tmh (Gobywriter_T writer, Stage3end_T stage3, Shortread_T queryseq, int npaths) {
#ifdef HAVE_GOBY
  unsigned long goby_read_index;
  UINT4 query_aligned_length;

  Goby_observe_aligned(writer);

  goby_read_index = (unsigned long) strtoul(Shortread_accession(queryseq), NULL, 10);
  query_aligned_length = Stage3end_query_alignment_length(stage3);
  gobyAlEntry_appendTooManyHits(writer->helper,goby_read_index,query_aligned_length,npaths);
#endif /* HAVE_GOBY */
  return;
}

void
Goby_print_pair_tmh (Gobywriter_T writer, Resulttype_T resulttype, Stage3pair_T stage3pair, Shortread_T queryseq, int npaths) {
#ifdef HAVE_GOBY
  unsigned long goby_read_index;
  UINT4 query_aligned_length;
    
  Goby_observe_aligned(writer);

  /* TODO: Is this correct for both cases (PAIRED_MULT, CONCORDANT_MULT)? */
  /* TODO: Q: Should we be outputting BOTH primary and mate TMH? */
  goby_read_index = (unsigned long) strtoul(Shortread_accession(queryseq), NULL, 10);
  query_aligned_length = Stage3end_query_alignment_length(Stage3pair_hit5(stage3pair));
  gobyAlEntry_appendTooManyHits(writer->helper,goby_read_index,query_aligned_length,npaths);

  goby_read_index = (unsigned long) strtoul(Shortread_accession(queryseq), NULL, 10);
  query_aligned_length = Stage3end_query_alignment_length(Stage3pair_hit3(stage3pair));
  gobyAlEntry_appendTooManyHits(writer->helper,goby_read_index,query_aligned_length,npaths);
#endif /* HAVE_GOBY */
    return;
}
