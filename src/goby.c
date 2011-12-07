static char rcsid[] = "$Id: goby.c 52336 2011-11-14 18:00:38Z twu $";
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
#include "complement.h"

/* #define DEBUG */
#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


static bool show_refdiff_p;

void
Goby_setup (bool show_refdiff_p_in) {
  show_refdiff_p = show_refdiff_p_in;
  return;
}


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


static char complCode[128] = COMPLEMENT_LC;

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
Goby_shutdown () {
#ifdef HAVE_GOBY
  goby_shutdownProtobuf();
#endif
  return;
}


/************************************************************************
 *   Reader
 ************************************************************************/

void
Goby_reader_free (Gobyreader_T *old) {
  FREE(*old);
  return;
}

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


void
Goby_reader_finish (Gobyreader_T reader) {
#ifdef HAVE_GOBY
  gobyReads_finished(reader->helper);
#endif
  return;
}


/************************************************************************
 *   Writer
 ************************************************************************/

void
Goby_writer_free (Gobywriter_T *old) {
#ifdef HAVE_GOBY
  /* gobyCapture_close((*old)->helper); */
#endif
  FREE(*old);
  return;
}

Gobywriter_T
Goby_writer_new (char *output_root, char *aligner_name, char *aligner_version) {
#ifdef HAVE_GOBY
  Gobywriter_T new = (Gobywriter_T) MALLOC(sizeof(*new));

  gobyAlignments_openAlignmentsWriterDefaultEntriesPerChunk(output_root,&new->helper);
  gobyAlignments_setAlignerName(new->helper,aligner_name);
  gobyAlignments_setAlignerVersion(new->helper,aligner_version);
  /* gobyCapture_open(new->helper, 1); */
  return new;
#else
  return NULL;
#endif
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
Goby_writer_finish (Gobywriter_T writer, Gobyreader_T reader) {
#ifdef HAVE_GOBY
  gobyAlignments_finished(writer->helper,reader->helper->numberOfReads);
#endif
  return;
}


/* This ifdef HAVE_GOBY spans MANY methods that are never called outside goby.c and are static */
#ifdef HAVE_GOBY

/**
 * Reverse and optionally complement a string of a specified length.
 * @param str the string to reverse (complement)
 * @param length the length of the string to reverse (complement)
 * @param if the characters should be complemented during the reversal
 * @return str, just for convenience
 */
static char *
reverse_complement (char *str, int length, bool complement) {
  bool is_odd;
  int midpoint, i;
  char temp;

  if (str == NULL || length == 0) {
    return (char *) NULL;
  }
  is_odd = (length % 2 == 1);
  midpoint = length / 2;
  for (i = 0; i < midpoint; i++) {
    temp = (complement ? complCode[(int) str[i]] : str[i]);
    str[i] = (complement ? complCode[(int) str[length - i - 1]] : str[length - i - 1]);
    str[length - i - 1] = temp;
  }
  if (is_odd && complement) {
    str[i] = complCode[(int) str[i]];
  }

  return str;
}

static char *
Goby_print_forward (char *string, int n) {
  return copy_string(string, n);
}

static char *
Goby_print_lc (char *string, int n) {
  char *copy = copy_string(string, n);
  int i;
  for (i = 0; i < n; i++) {
    copy[i] = tolower(copy[i]);
  }
  return copy;
}

static char *
Goby_print_revcomp (char *nt, int len) {
  char *copy = copy_string(nt, len);
  return reverse_complement(copy, len, true);
}

static char *
Goby_print_revcomp_lc (char *nt, int len) {
  int i;
  char *copy = Goby_print_revcomp (nt, len);
  for (i = 0; i < len; i++) {
    copy[i] = tolower(copy[i]);
  }
  return copy;
}

static char *
Goby_Shortread_print_oneline_uc (Shortread_T this) {
  return Goby_print_forward(Shortread_contents_uc(this), Shortread_fulllength(this));
}


static char *
Goby_Shortread_print_oneline_revcomp_uc (Shortread_T this) {
  return Goby_print_revcomp(Shortread_contents_uc(this), Shortread_fulllength(this));
}

static char *
merge_and_free_three (char *a, char *b, char *c) {
  int len_a = 0, len_b = 0, len_c = 0;

  if (a != NULL) {
    len_a = strlen(a);
  }
  if (b != NULL) {
    len_b = strlen(b);
  }
  if (c != NULL) {
    len_c = strlen(c);
  }

  int len_abc = len_a + len_b + len_c;
  char *copy = CALLOC(len_abc + 1, sizeof(char));
  copy[len_abc] = '\0';
  int pos = 0;
  int i;

  for (i = 0; i < len_a; i++) {
    copy[pos++] = a[i];
  }
  for (i = 0; i < len_b; i++) {
    copy[pos++] = b[i];
  }
  for (i = 0; i < len_c; i++) {
    copy[pos++] = c[i];
  }
  if (a != NULL) {
    FREE(a);
  }
  if (b != NULL) {
    FREE(b);
  }
  if (c != NULL) {
    FREE(c);
  }
  return copy;
}

/**
 * Even if this match is on the reverse strand, by the time the sequences in the match
 * get here they have been reverse complemented back to the forward strand.
 * output all sequence variations (insert, delete, mutations).
 */
static void
output_subs (Gobywriter_T writer, Hittype_T hittype, char *genomic, char *query, char *quality_string,
	     char *fasta_query, bool reverse_strand, int padding_left, int padding_right) {

  int nmismatches = 0, nindels = 0, i, query_i;
  char genomic_char, read_char, quality_char = '\0';
  int has_quality;
  unsigned long read_index, ref_position;
  int genomic_length = strlen(genomic);
  int fasta_length = strlen(fasta_query);
  int padded_length = padding_left + genomic_length + padding_right;
  bool too_big = false;

  /* Account for alignment padding on left and right */
  unsigned long *ref_positions = CALLOC(padded_length, sizeof(unsigned long));
  unsigned long *read_indexes = CALLOC(padded_length, sizeof(unsigned long));

  ref_position = 0;
  read_index = 0;
  for (i = 0; i < padded_length; i++) {
    if (i < padding_left) {
      /* In alignment padding ref_position doesn't increment but read position does */
      read_index++;
    } else if (i >= (genomic_length + padding_left)) {
      /* In alignment padding/clipping ref_position doesn't increment but read position does */
      read_index++;
    } else {
      genomic_char = toupper(genomic[i - padding_left]);
      if (genomic_char != '-') {
	ref_position++;
      }
      if (reverse_strand) {
	read_char = toupper(query[genomic_length - (i - padding_left) -  1]);
      } else {
	read_char = toupper(query[i - padding_left]);
      }
      if (read_char != '-') {
	read_index++;
      }
    }
    ref_positions[i] = ref_position;
    if (reverse_strand) {
      read_indexes[padded_length - i - 1] = read_index;
    } else {
      read_indexes[i] = read_index;
    }
  }

  debug(
	gobyAlignments_debugSequences(writer->helper, 0, genomic, query, padding_left, padding_right);
        fprintf(stderr, "::  pos=");
        for (i = 0; i < padded_length; i++) {
	  fprintf(stderr, "%lu", ref_positions[i] % 10);
        }
        fprintf(stderr, "\n::   ri=");
        for (i = 0; i < padded_length; i++) {
	  fprintf(stderr, "%lu", read_indexes[i] % 10);
        }
        fprintf(stderr, "\n");

        fprintf(stderr, "ref with positions\n");
        for (i = 0; i < padded_length; i++) {
	  if (i < padding_left) {
	    genomic_char = '_';
	  } else if (i >= (genomic_length + padding_left)) {
	    genomic_char = '_';
	  } else {
	    genomic_char = genomic[i - padding_left];
	  }
	  fprintf(stderr, "%03d:%c:%03lu  ", i, genomic_char, ref_positions[i]);
        }
        fprintf(stderr, "\n");

        fprintf(stderr, "read with positions\n");
        for (i = 0; i < padded_length; i++) {
	  if (i < padding_left) {
	    read_char = '_';
	  } else if (i >= (padding_left + genomic_length)) {
	    read_char = '_';
	  } else {
	    read_char = query[i - padding_left];
	  }
	  fprintf(stderr, "%03d:%c:%03lu  ", i, read_char, read_indexes[i]);
        }
        fprintf(stderr, "\n");
        fprintf(stderr, "padding_left=%d, padding_right=%d\n", padding_left, padding_right);
	);

  for (query_i = padding_left; query_i < padding_left + genomic_length; query_i++) {
    ref_position = ref_positions[query_i];
    read_index = read_indexes[query_i];
    i = query_i - padding_left;
    if (read_index > fasta_length) {
      too_big = true;
    }
    genomic_char = toupper(genomic[i]);
    read_char = toupper(query[i]);
    if (quality_string != NULL && read_char != '-') {
      quality_char = quality_string[i];
      has_quality = 1;
    } else {
      has_quality = 0;
    }
    if (genomic_char != read_char) {
      debug(
	    fprintf(stderr, "   Seqvar:\n");
	    fprintf(stderr, "      read_index:%u\n", read_index);
	    fprintf(stderr, "         ref_pos:%u\n", ref_position);
	    fprintf(stderr, "       read_char:%c\n", read_char);
	    fprintf(stderr, "    genomic_char:%c\n", genomic_char);
	    );

      gobyAlEntry_addSequenceVariation(writer->helper, read_index, ref_position, genomic_char, read_char,
				       has_quality, quality_char);
      if (genomic_char == '-' || read_char == '-') {
	nindels++;
      } else {
	nmismatches++;
      }
    }
  }

  if (too_big) {
    fprintf(stderr, " *** read_index [%lu] or ref_position [%lu] is too large! ***\n",
            read_index, ref_position);
    fprintf(stderr, ">%u\n", gobyAlEntry_getQueryIndex(writer->helper));
    if (fasta_query) {
      fprintf(stderr, "%s\n", fasta_query);
    } else {
      fprintf(stderr, "fasta_query was NULL\n");
    }
  }

  debug(
	fprintf(stderr, "     nmismatches:%d\n", nmismatches);
	fprintf(stderr, "           score:%d\n", genomic_length - nindels - nmismatches);
        );

  gobyAlEntry_setNumberOfMismatches(writer->helper,nmismatches);
  gobyAlEntry_setScoreInt(writer->helper,genomic_length - nindels - nmismatches);

  FREE(ref_positions);
  FREE(read_indexes);

  return;
}


static char *
merge_ref_substrings (Stage3end_T stage3) {
  Substring_T substring1, substring2;
  char *genomic1, *genomic2, *result;
  int i;

  substring1 = Stage3end_substring1(stage3);
  substring2 = Stage3end_substring2(stage3);

  genomic1 = Substring_genomic_refdiff(substring1);
  result = (char *) CALLOC(strlen(genomic1) + 1, sizeof(char));
  strcpy(result,genomic1);

  if (substring2 != NULL) {
    genomic2 = Substring_genomic_refdiff(substring2);
    for (i = 0; i < strlen(result); i++) {
      if (result[i] == '-') {
	result[i] = genomic2[i];
      }
    }
  }

  return result;
}

/**
 * In the case of a DELETION, insert the deletion bases into the ref
 * and '-'s into the read. genomic, query, (and qual if it exists) will
 * be replaced with new strings when deletion_size > 0.
 * @return the size of the deletion
 */
static int
compose_deletion (Stage3end_T stage3, char **genomic, char **query, char **qual) {
  int start, i,  deletion_length, new_length, query_length, cur_pos;
  char *deletion, *new_genomic, *new_query, *new_qual;

  start = Stage3end_indel_pos(stage3);
  deletion = Stage3end_deletion_string(stage3);
  deletion_length = strlen(deletion);
  if (deletion_length == 0) {
    /* ? No actual deletion */
    return 0;
  }
  query_length = strlen(*query);
  new_length = query_length + deletion_length;
  new_query = (char *) CALLOC(new_length + 1,sizeof(char));
  new_genomic = (char *) CALLOC(new_length + 1,sizeof(char));
  if (*qual != NULL) {
    new_qual = (char *) CALLOC(new_length + 1,sizeof(char));
  }

  for (cur_pos = 0, i = 0; i < start; i++, cur_pos++) {
    new_query[cur_pos] = (*query)[i];
    new_genomic[cur_pos] = (*genomic)[i];
    if (*qual != NULL) {
      new_qual[cur_pos] = (*qual)[i];
    }
  }

  for (i = 0; i < deletion_length; i++, cur_pos++) {
    new_query[cur_pos] = '-';
    new_genomic[cur_pos] = deletion[i];
    if (*qual != NULL) {
      new_qual[cur_pos] = '\0';
    }
  }

  for (i = start; i < query_length; i++, cur_pos++) {
    new_query[cur_pos] = (*query)[i];
    new_genomic[cur_pos] = (*genomic)[i];
    if (*qual != NULL) {
      new_qual[cur_pos] = (*qual)[i];
    }
  }
  new_query[cur_pos] = '\0';
  new_genomic[cur_pos] = '\0';
  if (*qual != NULL) {
    new_qual[cur_pos] = '\0';
  }

  /* ONLY free genomic, the other two came from gsnap. */
  FREE(*genomic);

  *genomic = new_genomic;
  *query = new_query;
  if (*qual != NULL) {
    *qual = new_qual;
  }
    
  return deletion_length;
}

static void
output_result (Gobywriter_T writer, Stage3end_T stage3, Shortread_T queryseq) {
  char *genomic, *query, *quality, *raw_genomic, *raw_query, *raw_quality, *fasta_query;
  int startpos = 0, length = 0, padding_left = 0, padding_right = 0;
  int temp;
  Hittype_T hittype;
  bool free_sequences, reverse_strand;

  hittype = Stage3end_hittype(stage3);
  reverse_strand = !Stage3end_plusp(stage3);
  if (hittype == EXACT) {
    startpos = 0;
    genomic = query = raw_query = Shortread_fullpointer(queryseq);
    length = strlen(genomic);
    quality = Shortread_quality_string(queryseq);
    free_sequences = false;
  } else {
    padding_left = startpos = Substring_querystart(Stage3end_substring1(stage3));
    raw_genomic = merge_ref_substrings(stage3);  /* must be free'd */
    fasta_query = raw_query = Shortread_fullpointer(queryseq);
    raw_quality = Shortread_quality_string(queryseq);
    length = Stage3end_query_alignment_length(stage3);
    padding_right = strlen(raw_genomic) - length  - startpos;
    debug(fprintf(stderr, ":: length=%d query_al_len=%d nindels=%d\n",length,Stage3end_query_alignment_length(stage3),Stage3end_nindels(stage3)));

    /* Compse the deletion into the sequences */
    debug(gobyAlignments_debugSequences(writer->helper,hittype,raw_genomic,raw_query, /*padding_left*/0, /*padding_left*/0));
    if (hittype == DELETION) {
      length += compose_deletion(stage3, &raw_genomic, &raw_query,  &raw_quality);
    }
  
    /* Remove the trimmed portion from the sequences */
    genomic = copy_string(&(raw_genomic[startpos]), length);
    FREE(raw_genomic);
    query = copy_string(&(raw_query[startpos]), length);
    if (hittype == DELETION) {
      FREE(raw_query);
    }
    if (raw_quality != NULL) {
      quality = copy_string(&(raw_quality[startpos]), length);
      if (hittype == DELETION) {
	FREE(raw_quality);
      }
    } else {
      quality = (char *) NULL;
    }
    free_sequences = true;
  }
  debug(gobyAlignments_debugSequences(writer->helper,hittype,genomic,query, padding_left, padding_right));

  if (hittype == DELETION || hittype == SUB || hittype == TERMINAL || hittype == INSERTION) {
    if (reverse_strand) {
      /* Both ref and reads are reverse complemented, put them back to the forward strand. */
      reverse_complement(query, length, true);
      if (genomic != query) {
	reverse_complement(genomic, length, true);
      }
      reverse_complement(quality, length, false);
      temp = padding_left;
      padding_left = padding_right;
      padding_right = temp;
      debug(fprintf(stderr, ":: Reverse complemented length=%d.\n", length));
      debug(gobyAlignments_debugSequences(writer->helper,hittype,genomic,query,padding_left,padding_right));
    }
    output_subs(writer,hittype,genomic,query,quality,fasta_query,reverse_strand,padding_left,padding_right);
  } else {
    /* Splicing results not supported */
  }

  if (free_sequences) {
    if (quality != NULL) {
      FREE(quality);
    }
    if (genomic != query) {
      FREE(genomic);
    }
    FREE(query);
  }

  return;
}

/**
 * Obtain the goby_chr_index from a stage3. The value will
 * be 0-based.
 */
static unsigned int
stage3_to_goby_chr_index(Stage3end_T stage3, IIT_T chromosome_iit) {
        Chrnum_T gsnap_chr_index;
        unsigned int goby_chr_index;
        char *label;
        bool allocp;

        /* Convert the gsnap chrom number to the goby chrom number */
        allocp = false;
        gsnap_chr_index = Stage3end_chrnum(stage3);
        debug(label = IIT_label(chromosome_iit,gsnap_chr_index,&allocp));
        goby_chr_index = gsnap_chr_index - 1;
        if (allocp) {
            FREE(label);
        }
        return goby_chr_index;
}

/**
 * Based on print_coordinates() in substring.c. This version differs because it always
 * returns such that start_coordinate < end_coordinate.
 */
static void
Goby_obtain_coordinates (Substring_T substring, bool invertp, Genomicpos_T *start_coordinate, Genomicpos_T *end_coordinate) {
        Genomicpos_T temp_coordinate;
        if (Substring_plusp(substring) == true) {
                if (invertp == false) {
                        *start_coordinate = Substring_alignstart_trim(substring) - Substring_chroffset(substring) + 1U;
                        *end_coordinate = Substring_alignend_trim(substring) - Substring_chroffset(substring);
                } else {
                        *start_coordinate = Substring_alignend_trim(substring) - Substring_chroffset(substring);
                        *end_coordinate = Substring_alignstart_trim(substring) - Substring_chroffset(substring) + 1U;
                }
        } else {
                if (invertp == false) {
                        *start_coordinate = Substring_alignstart_trim(substring) - Substring_chroffset(substring);
                        *end_coordinate = Substring_alignend_trim(substring) - Substring_chroffset(substring) + 1U;
                } else {
                        *start_coordinate = Substring_alignend_trim(substring) - Substring_chroffset(substring) + 1U;
                        *end_coordinate = Substring_alignstart_trim(substring) - Substring_chroffset(substring);
                }
        }

        if (*end_coordinate < *start_coordinate) {
                temp_coordinate = *start_coordinate;
                *start_coordinate = *end_coordinate;
                *end_coordinate = temp_coordinate;
        }

        return;
}

/**
 * TODO: Perhaps can be made faster by not reverse complementing and not including the lowercase and "---"'s.
 * This duplicates the functionality of print_genomic() but instead of printing it returns strings.
 */
static char *
Goby_print_genomic (Gobywriter_T writer, Substring_T substring, char *deletion, int deletionlength, bool invertp, Shortread_T queryseq) {
        char *a = NULL, *b = NULL, *c = NULL, *result;

        if (invertp == false) {
                if (Substring_genomic_bothdiff(substring) == NULL) {
                        /* Exact match */
                        result = Goby_Shortread_print_oneline_uc(queryseq);

                } else if (show_refdiff_p == true) {
                        a = Goby_print_forward(Substring_genomic_refdiff(substring),Substring_queryend(substring));
                        if (deletion != NULL) {
                                b = Goby_print_lc(deletion,deletionlength);
                        }
                        c = Goby_print_forward(&(Substring_genomic_refdiff(substring)[Substring_queryend(substring)]),Substring_querylength(substring) - Substring_queryend(substring));
                        result = merge_and_free_three(a, b, c);
                } else {
                        a = Goby_print_forward(Substring_genomic_bothdiff(substring),Substring_queryend(substring));
                        if (deletion != NULL) {
                                b = Goby_print_lc(deletion,deletionlength);
                        }
                        c = Goby_print_forward(&(Substring_genomic_bothdiff(substring)[Substring_queryend(substring)]),Substring_querylength(substring) - Substring_queryend(substring));
                        result = merge_and_free_three(a, b, c);
                }

        } else {
                if (Substring_genomic_bothdiff(substring) == NULL) {
                        /* Exact match */
                        result = Goby_Shortread_print_oneline_revcomp_uc(queryseq);

                } else if (show_refdiff_p == true) {
                        a = Goby_print_revcomp(&(Substring_genomic_refdiff(substring)[Substring_querystart(substring)]),Substring_querylength(substring) - Substring_querystart(substring));
                        if (deletion != NULL) {
                                b = Goby_print_revcomp_lc(deletion,deletionlength);
                        }
                        c = Goby_print_revcomp(Substring_genomic_refdiff(substring),Substring_querystart(substring));
                        result = merge_and_free_three(a, b, c);

                } else {
                        a = Goby_print_revcomp(&(Substring_genomic_bothdiff(substring)[Substring_querystart(substring)]),Substring_querylength(substring) - Substring_querystart(substring));
                        if (deletion != NULL) {
                                b = Goby_print_revcomp_lc(deletion,deletionlength);
                        }
                        c = Goby_print_revcomp(Substring_genomic_bothdiff(substring),Substring_querystart(substring));
                        result = merge_and_free_three(a, b, c);

                }
        }

        return result;
}

/**
 * Output a donor or an acceptor for a SPLICE, this one method
 * is combined to be able to output either.
 * This is loosely on Substring_print_donor() from stage3hr.c.
 */
static void
Goby_print_donor_or_acceptor (Gobywriter_T writer, Stage3end_T chimera,
                        Substring_T prev_dora, Substring_T current_dora, Substring_T next_dora,
                        unsigned int prev_frag_index, unsigned int current_frag_index, unsigned int next_frag_index,
                        bool invertp, Shortread_T queryseq,
                        IIT_T chromosome_iit,
                        int score, int mapq_score) {

        Genomicpos_T current_goby_position, prev_goby_position, next_goby_position, end_coordinate;
        char *genomic, *genomic_alloc;
        double chimera_prob;
        int trim_left = 0, trim_right = 0;
        int query_start, query_end, fasta_query_length, aligned_length;
        unsigned long goby_read_index;
        unsigned int goby_chr_index;
        unsigned int query_position;
        int nmismatches_bothdiff;
        unsigned int spliced_flags;
        char *read_comment;
        char *query, *quality, *fasta_query;
        char *query_alloc, *quality_alloc;
        bool reverse_strand;

        if (current_dora == NULL) {
	  /* Nothing to print */
	  return;
        }

        reverse_strand = (Stage3end_plusp(chimera) == true ? false : true);
        read_comment = Shortread_header(queryseq);
        goby_read_index = (unsigned long) strtoul(Shortread_accession(queryseq), NULL, 10);
        goby_chr_index = stage3_to_goby_chr_index(chimera, chromosome_iit);

        debug(fprintf(stderr, "  read_comment:%s\n",read_comment));

        genomic = genomic_alloc = Goby_print_genomic(writer, current_dora,/*deletion*/NULL,/*deletionlength*/0,invertp,queryseq);
        fasta_query = Shortread_fullpointer(queryseq);
        query = query_alloc = copy_string(fasta_query, -1);  /* Will become only the matching portion of the query */
        quality = quality_alloc = Shortread_quality_string(queryseq); /* Will become only the matching portion of the quality */
        if (quality_alloc != NULL) {
                quality = quality_alloc = copy_string(quality, strlen(query));
        }

        debug(
                fprintf(stderr, "   bc genomic:[%03d]%s\n", strlen(genomic), genomic);
                fprintf(stderr, "   bc   query:[%03d]%s\n", strlen(query), query);
        );

        Goby_obtain_coordinates(current_dora, invertp, &current_goby_position, &end_coordinate);
        current_goby_position--;

        if (prev_dora != NULL) {
                Goby_obtain_coordinates(prev_dora, invertp, &prev_goby_position, &end_coordinate);
                prev_goby_position--;
        }
        if (next_dora != NULL) {
                Goby_obtain_coordinates(next_dora, invertp, &next_goby_position, &end_coordinate);
                next_goby_position--;
        }

        chimera_prob = Substring_chimera_prob(current_dora);
        if (invertp == false) {
                query_start = 1 + Substring_querystart(current_dora);
                query_end = Substring_queryend(current_dora);
        } else {
                query_start = 1 + Substring_querylength(current_dora) - Substring_queryend(current_dora);
                query_end = Substring_querylength(current_dora) - Substring_querystart(current_dora);
        }

        fasta_query_length = strlen(fasta_query);
        aligned_length = query_end - query_start + 1;
        genomic[query_end] = '\0';
        genomic += (query_start - 1);

        query[query_end] = '\0';
        query += (query_start - 1);

        if (quality != NULL) {
                quality[query_end] = '\0';
                quality += (query_start - 1);
        }

        nmismatches_bothdiff = Substring_nmismatches_bothdiff(current_dora);
        spliced_flags = 1;  /* 1 for normal 2 for novel, if detected */

        debug(
                fprintf(stderr, "  invertp:%d\n",invertp);
                fprintf(stderr, "  reverse_strand:%d\n",(reverse_strand == true ? 1 : 0));
                fprintf(stderr, "  goby_read_index:%u\n",goby_read_index);
                fprintf(stderr, "  goby_chr_index:%u\n",goby_chr_index);
                fprintf(stderr, "      genomic:[%03d]%s\n", strlen(genomic), genomic);
                fprintf(stderr, "        query:[%03d]%s\n", strlen(query), query);
                fprintf(stderr, "     quaylity:[%03d]\n", quality != NULL ? strlen(quality) : 0);
                fprintf(stderr, "  fasta_query:[%03d]\n", strlen(fasta_query), fasta_query);
        );
        if (reverse_strand == true) {
	  /* With reverse strand matches */
                reverse_complement(genomic, strlen(genomic), true);
                reverse_complement(query, strlen(query), true);
                if (quality != NULL) {
                        reverse_complement(quality, strlen(query), false);
                }
                trim_left = fasta_query_length - query_end;
                trim_right = query_start - 1;
                debug(
                        fprintf(stderr, "      --Reverse Complement back to forward for reference--\n");
                        fprintf(stderr, "      genomic:[%03d]%s\n", strlen(genomic), genomic);
                        fprintf(stderr, "        query:[%03d]%s\n", strlen(query), query);
                );
        } else {
                trim_right = fasta_query_length - query_end;
                trim_left = query_start - 1;
        }
        query_position = trim_left;
        debug(
                fprintf(stderr, "  current_goby_position:%u\n", current_goby_position);
                fprintf(stderr, "  nmismatches_bothdiff:%d\n", nmismatches_bothdiff);
                fprintf(stderr, "  trim_left:%d\n", trim_left);
                fprintf(stderr, "  trim_right:%d\n", trim_right);
                fprintf(stderr, "  query_start:%d\n", query_start);
                fprintf(stderr, "  query_end:%d\n", query_end);
                fprintf(stderr, "  query_position:%d\n", query_position);
                fprintf(stderr, "  fasta_query_length:%d\n", fasta_query_length);
                fprintf(stderr, "  aligned_length:%d\n", aligned_length);
                fprintf(stderr, "  chimera_prob:%f\n", chimera_prob);
                fprintf(stderr, "  nindels:%d\n", Stage3end_nindels(chimera));
                fprintf(stderr, "  current_frag_index:%u\n",current_frag_index);
                if (prev_dora != NULL || next_dora != NULL) {
                        fprintf(stderr, "  Splice link:\n");
                        fprintf(stderr, "     spliced_flags:%u\n",spliced_flags);
                        if (prev_dora != NULL) {
                                fprintf(stderr, "     prev_frag_index:%u\n",prev_frag_index);
                                fprintf(stderr, "     prev_goby_position:%u\n",prev_goby_position);
                        }
                        if (next_dora != NULL) {
                                fprintf(stderr, "     next_frag_index:%u\n",next_frag_index);
                                fprintf(stderr, "     next_goby_position:%u\n",next_goby_position);
                        }
                }
        );

        gobyAlignments_appendEntry(writer->helper);
        gobyAlEntry_setMultiplicity(writer->helper, 1);
        gobyAlEntry_setQueryIndex(writer->helper, goby_read_index);
        gobyAlEntry_setTargetIndex(writer->helper, goby_chr_index);
        gobyAlEntry_setPosition(writer->helper, current_goby_position);
        gobyAlEntry_setMatchingReverseStrand(writer->helper, reverse_strand == true ? 1 : 0);
        gobyAlEntry_setQueryPosition(writer->helper, query_position);  /* query_start - 1 FIX THIS, the position of the start of match in query, probably has to do with trim */
        gobyAlEntry_setNumberOfIndels(writer->helper, 0);  /* There are no indels in this reporting method */
        gobyAlEntry_setQueryAlignedLength(writer->helper, aligned_length); /* As there are no indels in this reporting method, this should be accurate */
        gobyAlEntry_setTargetAlignedLength(writer->helper, aligned_length); /* As there are no indels in this reporting method, this should be accurate */
        gobyAlEntry_setQueryLength(writer->helper, fasta_query_length);  /* This is the fragment length, should be right */
        gobyAlEntry_setMappingQuality(writer->helper, mapq_score);
        gobyAlEntry_setFragmentIndex(writer->helper, current_frag_index);

        if (prev_dora != NULL || next_dora != NULL) {
                gobyAlEntry_setSplicedFlags(writer->helper, spliced_flags);
                if (prev_dora != NULL) {
                        gobyAlEntry_setSplicedBackwardFragmentIndex(writer->helper, prev_frag_index);
                        gobyAlEntry_setSplicedBackwardPosition(writer->helper, prev_goby_position);
                        gobyAlEntry_setSplicedBackwardTargetIndex(writer->helper, goby_chr_index);
                }
                if (next_dora != NULL) {
                        gobyAlEntry_setSplicedForwardFragmentIndex(writer->helper, next_frag_index);
                        gobyAlEntry_setSplicedForwardPosition(writer->helper, next_goby_position);
                        gobyAlEntry_setSplicedForwardTargetIndex(writer->helper, goby_chr_index);
                }
        }

        debug(
                fprintf(stderr, "     ---\n",current_goby_position);
                fprintf(stderr, "     gobyAlEntry_setPosition:%u\n", current_goby_position);
                fprintf(stderr, "     gobyAlEntry_setQueryPosition:%d\n", query_position);
                fprintf(stderr, "     gobyAlEntry_setQueryAlignedLength:%d\n", aligned_length);
                fprintf(stderr, "     gobyAlEntry_setTargetAlignedLength:%d\n", aligned_length);
                fprintf(stderr, "     gobyAlEntry_setQueryLength:%d\n", fasta_query_length);
        );
        output_subs(writer, Stage3end_hittype(chimera),
                        genomic, query, quality, fasta_query,
                        reverse_strand, trim_left, trim_right);
        FREE(genomic_alloc);
        FREE(query_alloc);
        if (quality_alloc != NULL) {
                FREE(quality_alloc);
        }
        return;
}

/**
 * Based on print_pair_info in stage3hr.c. Only validation details kept.
 */
static void
Goby_validate_pair_info (Stage3end_T hit5, Stage3end_T hit3, Pairtype_T pairtype) {

        assert(Stage3end_effective_chrnum(hit5) == Stage3end_effective_chrnum(hit3)); /* Same chromosomes */

#if 0
        /* Doesn't hold for paired (inversion) */
        assert(Stage3end_plusp(hit5) == Stage3end_plusp(hit3));   /* Same direction */
#endif

        switch (pairtype) {
                case CONCORDANT: break;
                case PAIRED_SCRAMBLE: break;
                case PAIRED_INVERSION: break;
                case PAIRED_TOOLONG: break;
                case TRANSLOCATION: break;
                case PAIRED_UNSPECIFIED: abort();
                case UNPAIRED: abort();
        }

        return;
}

void sort_output_blocks(int array_size, Substring_T substrings[], Genomicpos_T starts[]) {
        int i, j;
        Substring_T temp_substring;
        Genomicpos_T temp_start;

        /* Perform a simple bubble sort based on genomic position. */
        for (i = (array_size - 1); i > 0; i--) {
                for (j = 1; j <= i; j++) {
                        if (starts[j - 1] > starts[j]) {
                                temp_start = starts[j - 1];
                                starts[j - 1] = starts[j];
                                starts[j] = temp_start;
                                temp_substring = substrings[j - 1];
                                substrings[j - 1] = substrings[j];
                                substrings[j] = temp_substring;
                        }
                }
        }
}

/**
 * Based on print_splice in stage3hr.c.
 */
static void
Goby_print_splice (Gobywriter_T writer, Stage3end_T chimera, int score,
              IIT_T chromosome_iit, Shortread_T queryseq, bool invertp, Stage3end_T hit5, Stage3end_T hit3,
              int insertlength, int pairscore, Pairtype_T pairtype, int mapq_score) {

        debug(fprintf(stderr, "-------- NEW SPLICE --------\n"));
        int max_blocks = 2;
        Substring_T output_substrings[max_blocks];
        Genomicpos_T output_start_coords[max_blocks], start_coordinate, end_coordinate;
        int nblocks = 0, i;

        /* Pairing validation */
        if (hit5 != NULL && hit3 != NULL) {
                Goby_validate_pair_info(hit5,hit3,pairtype);
        }

        /* Initialize our array */
        for (i = 0; i < max_blocks; i++) {
                output_substrings[i] = NULL;
                output_start_coords[i] = 0;
        }

        if (Stage3end_hittype(chimera) == HALFSPLICE_DONOR) {
                output_substrings[nblocks] = Stage3end_substring1(chimera);
                Substring_assign_donor_prob(output_substrings[nblocks]);
                nblocks++;
        } else if (Stage3end_hittype(chimera) == HALFSPLICE_ACCEPTOR) {
                output_substrings[nblocks] = Stage3end_substring1(chimera);
                Substring_assign_acceptor_prob(output_substrings[nblocks]);
                nblocks++;
        } else {
                output_substrings[nblocks] = Stage3end_substring1(chimera);
                Substring_assign_donor_prob(output_substrings[nblocks]);
                nblocks++;
                output_substrings[nblocks] = Stage3end_substring2(chimera);
                Substring_assign_acceptor_prob(output_substrings[nblocks]);
                nblocks++;
                }
        /* Determine start coordinates for each of the blocks */
        for (i = 0; i < nblocks; i++) {
                Goby_obtain_coordinates(output_substrings[i], invertp, &start_coordinate, &end_coordinate);
                output_start_coords[i] = start_coordinate;
        }

        /* Order the pieces to output by genomic position */
        sort_output_blocks(nblocks, output_substrings, output_start_coords);

        /* Output the pieces */
        Goby_print_donor_or_acceptor(writer, chimera,
                        /*prev_dora*/NULL, /*current_dora*/output_substrings[0], /*next_dora*/output_substrings[1],
                /*previous_frag_index*/ 0, /*current_frag_index*/ 0, /*next_frag_index*/1,
                        invertp, queryseq, chromosome_iit, score, mapq_score);
        if (nblocks == 2) {
                Goby_print_donor_or_acceptor(writer, chimera,
                                /*prev_dora*/output_substrings[0], /*current_dora*/output_substrings[1], /*next_dora*/NULL,
                                /*previous_frag_index*/0, /*current_frag_index*/1, /*next_frag_index*/2,
                                invertp, queryseq, chromosome_iit, score, mapq_score);
        }

        return;
}

/**
 * Based on print_shortexon in stage3hr.c.
 */
static void
Goby_print_shortexon (Gobywriter_T writer, Stage3end_T chimera, int score,
		      IIT_T chromosome_iit, Shortread_T queryseq, bool invertp, Stage3end_T hit5, Stage3end_T hit3,
		      int insertlength, int pairscore, Pairtype_T pairtype, int mapq_score) {
        Substring_T donor=NULL, acceptor=NULL, shortexon;
        int max_blocks = 3;
        Substring_T output_substrings[max_blocks];
        Genomicpos_T output_start_coords[max_blocks], start_coordinate, end_coordinate;

        bool firstp = true;
        int nblocks = 0, i;

        debug(fprintf(stderr, "-------- NEW SHORTEXON --------\n"));

        /* Initialize our array */
        for (i = 0; i < max_blocks; i++) {
                output_substrings[i] = NULL;
                output_start_coords[i] = 0;
        }

        shortexon = Stage3end_substring1(chimera);
        Substring_assign_shortexon_prob(shortexon);
        output_substrings[0] = shortexon;
        nblocks++;
        if ((donor = Stage3end_substringD(chimera)) != NULL) {
                Substring_assign_donor_prob(donor);
                output_substrings[nblocks] = donor;
                nblocks++;
        }
        if ((acceptor = Stage3end_substringA(chimera)) != NULL) {
                Substring_assign_acceptor_prob(acceptor);
                output_substrings[nblocks] = acceptor;
                nblocks++;
        }
        /* Determine start coordinates for each of the blocks */
        for (i = 0; i < nblocks; i++) {
                Goby_obtain_coordinates(output_substrings[i], invertp, &start_coordinate, &end_coordinate);
                output_start_coords[i] = start_coordinate;
        }

        /* Order the pieces to output by genomic position */
        sort_output_blocks(nblocks, output_substrings, output_start_coords);

        /* Output the pieces */
        Goby_print_donor_or_acceptor(writer, chimera,
                        /*prev_dora*/NULL, /*current_dora*/output_substrings[0], /*next_dora*/output_substrings[1],
                        /*previous_frag_index*/0, /*current_frag_index*/0, /*next_frag_index*/1,
                        invertp, queryseq, chromosome_iit, score, mapq_score);
        if (nblocks >= 2) {
                Goby_print_donor_or_acceptor(writer, chimera,
                                /*prev_dora*/output_substrings[0], /*current_dora*/output_substrings[1], /*next_dora*/output_substrings[2],
                                /*previous_frag_index*/0, /*current_frag_index*/1, /*next_frag_index*/2,
                                invertp, queryseq, chromosome_iit, score, mapq_score);
        }
        if (nblocks == 3) {
                Goby_print_donor_or_acceptor(writer, chimera,
                                /*prev_dora*/output_substrings[1], /*current_dora*/output_substrings[2], /*next_dora*/NULL,
                                /*previous_frag_index*/1, /*current_frag_index*/2, /*next_frag_index*/0,
                                invertp, queryseq, chromosome_iit, score, mapq_score);
        }

  return;
}

/**
 * Based on Pair_print_gsnap. 
 * I am passing in stage3 and queryseq. Other parameters were
 * passed to Pair_print_gsnap() will be generated from stage3 and queryseq.
 */
static void
Goby_print_gmap_pairarray (Gobywriter_T writer, Stage3end_T stage3, Shortread_T queryseq, int insertlength, int pairscore, IIT_T chromosome_iit) {

        debug(fprintf(stderr, "-------- NEW GMAP --------\n"));

  return;
}

#endif /* HAVE_GOBY */

void
Goby_observe_aligned(Gobywriter_T writer) {
#ifdef HAVE_GOBY
        writer->helper->numberOfAlignedReads++;
#endif /* HAVE_GOBY */
        return;
}


void
Goby_print_tmh (Gobywriter_T writer, Stage3end_T stage3, Shortread_T queryseq, int npaths) {
#ifdef HAVE_GOBY
        unsigned long goby_read_index;
        UINT4 query_aligned_length;

        goby_read_index = (unsigned long) strtoul(Shortread_accession(queryseq), NULL, 10);
        query_aligned_length = Stage3end_query_alignment_length(stage3);
        gobyAlEntry_appendTooManyHits(writer->helper,goby_read_index,query_aligned_length,npaths);
#endif /* HAVE_GOBY */
        return;
}


/* Assume that stage3array has already been sorted */
void
Goby_print_single (Gobywriter_T writer, Stage3end_T stage3, int score,
		   IIT_T chromosome_iit, Shortread_T queryseq,
		   bool invertp, Stage3end_T hit5, Stage3end_T hit3,
		   int insertlength, int pairscore, Pairtype_T pairtype,
		   int mapq_score) {
#ifdef HAVE_GOBY
  Substring_T substring1;
        unsigned int goby_chr_index;
  unsigned long goby_read_index;
  Hittype_T hittype;

  hittype = Stage3end_hittype(stage3);

  goby_read_index = (unsigned long) strtoul(Shortread_accession(queryseq), NULL, 10);

        debug(fprintf(stderr,"+ [readIndex=%u] Writing hittype=%d, %s\n",
                        goby_read_index, hittype, Stage3end_hittype_string(stage3)));
  if (hittype == ONE_THIRD_SHORTEXON || hittype == TWO_THIRDS_SHORTEXON || hittype == SHORTEXON) {
    Goby_print_shortexon(writer,stage3,score,
			 chromosome_iit,queryseq,invertp,hit5,hit3,insertlength,
			 pairscore,pairtype,mapq_score);

  } else if (hittype == HALFSPLICE_DONOR || hittype == HALFSPLICE_ACCEPTOR || hittype == SPLICE) {
    Goby_print_splice(writer,stage3,score,chromosome_iit,queryseq,
		      invertp,hit5,hit3,insertlength,pairscore,
		      pairtype,mapq_score);

  } else {
    substring1 = Stage3end_substring1(stage3);
                goby_chr_index = stage3_to_goby_chr_index(stage3, chromosome_iit);
    
    gobyAlignments_appendEntry(writer->helper);
    gobyAlEntry_setMultiplicity(writer->helper,1);
    gobyAlEntry_setQueryIndex(writer->helper,goby_read_index);
    gobyAlEntry_setTargetIndex(writer->helper,goby_chr_index);
    gobyAlEntry_setPosition(writer->helper,Stage3end_chrpos_low_trim(stage3));
    gobyAlEntry_setMatchingReverseStrand(writer->helper,(Stage3end_plusp(stage3) == 0 ? 1 : 0));
    gobyAlEntry_setQueryPosition(writer->helper,Substring_querystart(substring1));
    gobyAlEntry_setScoreInt(writer->helper,Substring_querylength(substring1) - score);
    gobyAlEntry_setNumberOfMismatches(writer->helper,Stage3end_nmismatches_refdiff(stage3));
    gobyAlEntry_setNumberOfIndels(writer->helper,Stage3end_nindels(stage3));
    gobyAlEntry_setQueryAlignedLength(writer->helper,(UINT4) Stage3end_query_alignment_length(stage3));
    gobyAlEntry_setTargetAlignedLength(writer->helper,Stage3end_genomic_alignment_length(stage3));
    gobyAlEntry_setQueryLength(writer->helper,Substring_querylength(substring1));
                gobyAlEntry_setFragmentIndex(writer->helper, 0);
    gobyAlEntry_setMappingQuality(writer->helper,mapq_score);
    if (hittype != GMAP) {
      output_result(writer,stage3,queryseq);
    }
  }

#endif
  return;
}

/* Assumption: mate will be NULL if the second half of the pair doesn't align */
void
Goby_print_pair (Gobywriter_T writer, Stage3end_T this, Stage3end_T mate, char *acc, int pathnum, int npaths,
		 int mapq_score, IIT_T chromosome_iit, Shortread_T queryseq,
		 Shortread_T queryseq_mate, int pairedlength, Resulttype_T resulttype,
		 bool first_read_p, int npaths_mate, int quality_shift,
		 char *sam_read_group_id, bool invertp, bool invert_mate_p) {
#ifdef HAVE_GOBY
  Substring_T substring1;
  unsigned long goby_read_index;
  unsigned int sam_flags = 0U;
  unsigned int fragment_index = 0;
  unsigned int m_fragment_index = 0;
        unsigned int goby_chr_index;
  bool plusp;

  substring1 = Stage3end_substring1(this);
  plusp = Stage3end_plusp(this);
  sam_flags = SAM_compute_flag(plusp,mate,resulttype,first_read_p,
			       pathnum,npaths,npaths_mate,invertp,invert_mate_p);

  if (sam_flags & FIRST_READ_P) {
    fragment_index = 0;
    m_fragment_index = 1;
  } else if (sam_flags & SECOND_READ_P) {
    fragment_index = 1;
    m_fragment_index = 0;
  }

  goby_read_index = (unsigned long) strtoul(acc, NULL, 10);
  goby_chr_index = stage3_to_goby_chr_index(this, chromosome_iit);

  if (Stage3end_hittype(this) == SPLICE) {
                debug(fprintf(stderr,"Goby does not yet support hittype of SPLICE in conjunction with paired-end alignments\n"));
  } else {
    gobyAlignments_appendEntry(writer->helper);
    gobyAlEntry_setMultiplicity(writer->helper,1);
    gobyAlEntry_setQueryIndex(writer->helper,goby_read_index);
    gobyAlEntry_setTargetIndex(writer->helper,goby_chr_index);
    gobyAlEntry_setPosition(writer->helper,Stage3end_chrpos_low_trim(this));
    gobyAlEntry_setMatchingReverseStrand(writer->helper,(plusp == 0 ? 1 : 0));
    gobyAlEntry_setQueryPosition(writer->helper,Substring_querystart(substring1));
    gobyAlEntry_setScoreInt(writer->helper,Substring_querylength(substring1) - Stage3end_score(this));
    gobyAlEntry_setNumberOfMismatches(writer->helper,Stage3end_nmismatches_refdiff(this));
    gobyAlEntry_setNumberOfIndels(writer->helper,Stage3end_nindels(this));
    gobyAlEntry_setQueryAlignedLength(writer->helper,(UINT4) Stage3end_query_alignment_length(this));
    gobyAlEntry_setTargetAlignedLength(writer->helper,Stage3end_genomic_alignment_length(this));
    gobyAlEntry_setQueryLength(writer->helper,Substring_querylength(substring1));
    gobyAlEntry_setFragmentIndex(writer->helper,fragment_index);
    gobyAlEntry_setMappingQuality(writer->helper,mapq_score);

    gobyAlEntry_setPairFlags(writer->helper,sam_flags);
    if (mate) {
      gobyAlEntry_setPairFragmentIndex(writer->helper,m_fragment_index);
      gobyAlEntry_setPairTargetIndex(writer->helper,Stage3end_chrnum(mate) - 1);
      gobyAlEntry_setPairPosition(writer->helper,Stage3end_chrpos_low_trim(mate));
    }

    output_result(writer,this,queryseq);
  }
#endif

  return;
}

void
Goby_print_paired (Gobywriter_T writer, Result_T result, Resulttype_T resulttype,
		   IIT_T chromosome_iit, Shortread_T queryseq1, Shortread_T queryseq2,
		   int maxpaths, bool quiet_if_excessive_p,
		   bool invert_first_p, bool invert_second_p,
		   bool nofailsp, bool failsonlyp, bool fails_as_input_p,
		   bool fastq_format_p, int quality_shift, char *sam_read_group_id) {
#ifdef HAVE_GOBY
  Stage3pair_T *stage3pairarray, stage3pair;
  Stage3end_T *stage3array1, *stage3array2, stage3, mate;
  int npaths, npaths1, npaths2, pathnum;
  int second_absmq;
  char *acc;
  UINT4 query_aligned_length;
  unsigned long goby_read_index;
    
  acc = Shortread_accession(queryseq1);

  if (resulttype == CONCORDANT_UNIQ || resulttype == CONCORDANT_TRANSLOC) {
    stage3pairarray = (Stage3pair_T *) Result_array(&npaths,&second_absmq,result);
    /* Stage3pair_eval(stage3pairarray,npaths,maxpaths,queryseq1,queryseq2); */

    stage3pair = stage3pairarray[0];

    /* print first end */
    Goby_print_pair(writer,Stage3pair_hit5(stage3pair),/*mate*/Stage3pair_hit3(stage3pair),acc,
		    /*pathnum*/1,/*npaths*/1,Stage3pair_mapq_score(stage3pair),
		    chromosome_iit,/*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
		    Stage3pair_pairlength(stage3pair),resulttype,
		    /*first_read_p*/true,/*npaths_mate*/npaths,quality_shift,sam_read_group_id,
		    invert_first_p,invert_second_p);

    /* print second end */
    Goby_print_pair(writer,Stage3pair_hit3(stage3pair),/*mate*/Stage3pair_hit5(stage3pair),acc,
		    /*pathnum*/1,/*npaths*/1,Stage3pair_mapq_score(stage3pair),
		    chromosome_iit,/*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
		    Stage3pair_pairlength(stage3pair),resulttype,
		    /*first_read_p*/false,/*npaths_mate*/npaths,quality_shift,sam_read_group_id,
		    invert_second_p,invert_first_p);

  } else if (resulttype == CONCORDANT_MULT) {
    stage3pairarray = (Stage3pair_T *) Result_array(&npaths,&second_absmq,result);

    if (npaths > maxpaths) {
      /* No output if excessive for Gobyweb, but output TMH. */
      /* TODO: Q: Should we be outputting BOTH primary and mate TMH? */
      goby_read_index = (unsigned long) strtoul(Shortread_accession(queryseq1), NULL, 10);
      query_aligned_length = Stage3end_query_alignment_length(Stage3pair_hit5(stage3pair));
      gobyAlEntry_appendTooManyHits(writer->helper,goby_read_index,query_aligned_length,npaths);

      goby_read_index = (unsigned long) strtoul(Shortread_accession(queryseq2), NULL, 10);
      query_aligned_length = Stage3end_query_alignment_length(Stage3pair_hit3(stage3pair));
      gobyAlEntry_appendTooManyHits(writer->helper,goby_read_index,query_aligned_length,npaths);

    } else {
      /* Stage3pair_eval(stage3pairarray,npaths,maxpaths,queryseq1,queryseq2); */
      for (pathnum = 1; pathnum <= npaths && pathnum <= maxpaths; pathnum++) {
	stage3pair = stage3pairarray[pathnum-1];

	/* print first end */
	Goby_print_pair(writer,Stage3pair_hit5(stage3pair),/*mate*/Stage3pair_hit3(stage3pair),acc,
                        pathnum,npaths,Stage3pair_mapq_score(stage3pair),
                        chromosome_iit,/*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
                        Stage3pair_pairlength(stage3pair),resulttype,
                        /*first_read_p*/true,/*npaths_mate*/npaths,quality_shift,sam_read_group_id,
                        invert_first_p,invert_second_p);

	/* print second end */
	Goby_print_pair(writer,Stage3pair_hit3(stage3pair),/*mate*/Stage3pair_hit5(stage3pair),acc,
                        pathnum,npaths,Stage3pair_mapq_score(stage3pair),
                        chromosome_iit,/*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
                        Stage3pair_pairlength(stage3pair),resulttype,
                        /*first_read_p*/false,/*npaths_mate*/npaths,quality_shift,sam_read_group_id,
                        invert_second_p,invert_first_p);
      }
    }

  } else if (resulttype == UNPAIRED_UNIQ || resulttype == UNPAIRED_TRANSLOC) {
    /* Should print mate information in this situation */
    stage3array1 = (Stage3end_T *) Result_array(&npaths1,&second_absmq,result);
    stage3array2 = (Stage3end_T *) Result_array2(&npaths2,&second_absmq,result);

    /* print first end */
    /* Stage3end_eval_and_sort(stage3array1,npaths1,maxpaths,queryseq1); */
    Goby_print_pair(writer,stage3array1[0],/*mate*/stage3array2[0],acc,
		    /*pathnum*/1,/*npaths*/1,Stage3end_mapq_score(stage3array1[0]),
		    chromosome_iit,/*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
		    /*pairedlength*/0U,resulttype,
		    /*first_read_p*/true,/*npaths_mate*/1,quality_shift,sam_read_group_id,
		    invert_first_p,invert_second_p);

    /* print second end */
    /* Stage3end_eval_and_sort(stage3array2,npaths2,maxpaths,queryseq2); */
    Goby_print_pair(writer,stage3array2[0],/*mate*/stage3array1[0],acc,
		    /*pathnum*/1,/*npaths*/1,Stage3end_mapq_score(stage3array2[0]),
		    chromosome_iit,/*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
		    /*pairedlength*/0U,resulttype,
		    /*first_read_p*/false,/*npaths_mate*/1,quality_shift,sam_read_group_id,
		    invert_second_p,invert_first_p);

  } else if (resulttype == UNPAIRED_MULT) {
    stage3array1 = (Stage3end_T *) Result_array(&npaths1,&second_absmq,result);
    stage3array2 = (Stage3end_T *) Result_array2(&npaths2,&second_absmq,result);

#if 0
    /* Do eval and sorting first */
    if (npaths1 == 1) {
      /* Stage3end_eval_and_sort(stage3array1,npaths1,maxpaths,queryseq1); */
    } else if (npaths1 > maxpaths) {
      /* Don't sort */
    } else {
      /* Stage3end_eval_and_sort(stage3array1,npaths1,maxpaths,queryseq1); */
    }
#endif

#if 0
    if (npaths2 == 1) {
      /* Stage3end_eval_and_sort(stage3array2,npaths2,maxpaths,queryseq2); */
    } else if (npaths2 > maxpaths) {
      /* Don't sort */
    } else {
      /* Stage3end_eval_and_sort(stage3array2,npaths2,maxpaths,queryseq2); */
    }
#endif

    /* print first end results */
    mate = (npaths2 == 0) ? (Stage3end_T) NULL : stage3array2[0];

    if (npaths1 == 1) {
      stage3 = stage3array1[0];
      Goby_print_pair(writer,stage3,mate,acc,/*pathnum*/1,npaths1,Stage3end_mapq_score(stage3),
		      chromosome_iit,/*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
		      /*pairedlength*/0U,resulttype,/*first_read_p*/true,
                      /*npaths_mate*/npaths2,quality_shift,sam_read_group_id,
                      invert_first_p,invert_second_p);

    } else if (npaths1 > maxpaths) {
      /** No output if excessive for Gobyweb, but output TMH. */
      goby_read_index = (unsigned long) strtoul(Shortread_accession(queryseq1), NULL, 10);
      query_aligned_length = Stage3end_query_alignment_length(stage3array1[0]);
      gobyAlEntry_appendTooManyHits(writer->helper,goby_read_index,query_aligned_length,npaths1);

    } else {
      for (pathnum = 1; pathnum <= npaths1 && pathnum <= maxpaths; pathnum++) {
	stage3 = stage3array1[pathnum-1];
	Goby_print_pair(writer,stage3,mate,acc,pathnum,npaths1,Stage3end_mapq_score(stage3),
                        chromosome_iit,/*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
                        /*pairedlength*/0U,resulttype,/*first_read_p*/true,/*npaths_mate*/npaths2,
                        quality_shift,sam_read_group_id,invert_first_p,invert_second_p);
      }
    }

    /* print second end results */
    mate = (npaths1 == 0) ? (Stage3end_T) NULL : stage3array1[0];

    if (npaths2 == 1) {
      stage3 = stage3array2[0];
      Goby_print_pair(writer,stage3,mate,acc,/*pathnum*/1,npaths2,Stage3end_mapq_score(stage3),
		      chromosome_iit,/*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
		      /*pairedlength*/0U,resulttype,/*first_read_p*/false,/*npaths_mate*/npaths1,
		      quality_shift,sam_read_group_id,invert_second_p,invert_first_p);

    } else if (npaths2 > maxpaths) {
      /** No output if excessive for Gobyweb, but output TMH. */
      goby_read_index = (unsigned long) strtoul(Shortread_accession(queryseq2), NULL, 10);
      query_aligned_length = Stage3end_query_alignment_length(stage3array2[0]);
      gobyAlEntry_appendTooManyHits(writer->helper,goby_read_index,query_aligned_length,npaths2);

    } else {
      for (pathnum = 1; pathnum <= npaths2 && pathnum <= maxpaths; pathnum++) {
	stage3 = stage3array2[pathnum-1];
	Goby_print_pair(writer,stage3,mate,acc,pathnum,npaths2,Stage3end_mapq_score(stage3),
                        chromosome_iit,/*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
                        /*pairedlength*/0U,resulttype,/*first_read_p*/false,
                        /*npaths_mate*/npaths1,quality_shift,sam_read_group_id,
                        invert_second_p,invert_first_p);
      }
    }

  } else {
    if (resulttype == HALFMAPPING_UNIQ || resulttype == HALFMAPPING_TRANSLOC || resulttype == HALFMAPPING_MULT) {
      /* These are the last two resulttypes we can deal with */
    } else {
      abort();
    }

    stage3array1 = (Stage3end_T *) Result_array(&npaths1,&second_absmq,result);
    stage3array2 = (Stage3end_T *) Result_array2(&npaths2,&second_absmq,result);

#if 0
    /* Do eval and sorting first */
    if (npaths1 == 0) {
      /* Nothing to sort */
    } else if (npaths1 == 1) {
      /* Stage3end_eval_and_sort(stage3array1,npaths1,maxpaths,queryseq1); */
    } else if (npaths1 > maxpaths) {
      /* Don't sort */
    } else {
      /* Stage3end_eval_and_sort(stage3array1,npaths1,maxpaths,queryseq1); */
    }
#endif

#if 0
    if (npaths2 == 0) {
      /* Nothing to sort */
    } else if (npaths2 == 1) {
      /* Stage3end_eval_and_sort(stage3array2,npaths2,maxpaths,queryseq2); */
    } else if (npaths2 > maxpaths) {
      /* Don't sort */
    } else {
      /* Stage3end_eval_and_sort(stage3array2,npaths2,maxpaths,queryseq2); */
    }
#endif

    /* print first end results */
    mate = (npaths2 == 0) ? (Stage3end_T) NULL : stage3array2[0];

    if (npaths1 == 0) {
      /** No output for Gobyweb. */

    } else if (npaths1 == 1) {
      /* mate should be NULL here */

      stage3 = stage3array1[0];
      Goby_print_pair(writer,stage3,mate,acc,/*pathnum*/1,npaths1,Stage3end_mapq_score(stage3),
		      chromosome_iit,/*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
		      /*pairedlength*/0U,resulttype,/*first_read_p*/true,
                      /*npaths_mate*/npaths2, quality_shift,sam_read_group_id,
                      invert_first_p,invert_second_p);

    } else if (npaths1 > maxpaths) {
      /** No output if excessive for Gobyweb, but output TMH. */
      goby_read_index = (unsigned long) strtoul(Shortread_accession(queryseq1), NULL, 10);
      query_aligned_length = Stage3end_query_alignment_length(stage3array1[0]);
      gobyAlEntry_appendTooManyHits(writer->helper,goby_read_index,query_aligned_length,npaths1);

    } else {
      /* mate should be NULL here */
      for (pathnum = 1; pathnum <= npaths1 && pathnum <= maxpaths; pathnum++) {
	stage3 = stage3array1[pathnum-1];
	Goby_print_pair(writer,stage3,mate,acc,pathnum,npaths1,Stage3end_mapq_score(stage3),
                        chromosome_iit,/*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
                        /*pairedlength*/0U,resulttype,/*first_read_p*/true,
                        /*npaths_mate*/npaths2,quality_shift,sam_read_group_id,
                        invert_first_p,invert_second_p);
      }
    }

    /* print second end results */
    mate = (npaths1 == 0) ? (Stage3end_T) NULL : stage3array1[0];

    if (npaths2 == 0) {
      /* No output for Gobyweb. */

    } else if (npaths2 == 1) {
      /* mate should be NULL here */

      stage3 = stage3array2[0];
      Goby_print_pair(writer,stage3,mate,acc,/*pathnum*/1,npaths2,Stage3end_mapq_score(stage3),
		      chromosome_iit,/*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
		      /*pairedlength*/0U,resulttype,/*first_read_p*/false,
                      /*npaths_mate*/npaths1, quality_shift,sam_read_group_id,
                      invert_second_p,invert_first_p);

    } else if (npaths2 > maxpaths) {
      /* No output if excessive for Gobyweb, but output TMH. */
      goby_read_index = (unsigned long) strtoul(Shortread_accession(queryseq2), NULL, 10);
      query_aligned_length = Stage3end_query_alignment_length(stage3array2[0]);
      gobyAlEntry_appendTooManyHits(writer->helper,goby_read_index,query_aligned_length,npaths2);

    } else {
      /* mate should be NULL here */
      for (pathnum = 1; pathnum <= npaths2 && pathnum <= maxpaths; pathnum++) {
	stage3 = stage3array2[pathnum-1];
	Goby_print_pair(writer,stage3,mate,acc,pathnum,npaths2,Stage3end_mapq_score(stage3),
                        chromosome_iit,/*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
                        /*pairedlength*/0U,resulttype,/*first_read_p*/false,
                        /*npaths_mate*/npaths1, quality_shift,sam_read_group_id,
                        invert_second_p,invert_first_p);
      }
    }

  }
#endif
    
    return;
}
