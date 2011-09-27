static char rcsid[] = "$Id: goby.c 37254 2011-03-28 16:34:08Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "goby.h"

#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>

#include "mem.h"
#include "chrnum.h"
#include "substring.h"
#include "samflags.h"
#include "samprint.h"
#include "complement.h"

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

#ifdef HAVE_GOBY
#include <goby/C_Reads.h>
#include <goby/C_Alignments.h>
#include <goby/C_CompactHelpers.h>

struct Gobyreader_T {
  CReadsHelper *helper;
};

struct Gobywriter_T {
  CAlignmentsWriterHelper *helper;
};

#endif


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
Goby_reader_new (char **files, int nfiles, unsigned long window_start, unsigned long window_end) {
#ifdef HAVE_GOBY
  Gobyreader_T new = (Gobyreader_T) MALLOC(sizeof(*new));
  fprintf(stderr,"Opening %s\n",files[0]);
  gobyReads_openReadsReaderWindowed(files,nfiles,/*circularp*/false,window_start,window_end,&new->helper);
  return new;
#else
  return NULL;
#endif
}

static char *
copy_string (char *str, int length) {
  int copy_length = length;
  char *new_str = (char *) NULL;

  if (str != NULL) {
    if (copy_length == -1) {
      copy_length = strlen(str);
    }
    new_str = CALLOC(copy_length + 1, sizeof(char));
    strncpy(new_str, str, copy_length);
    new_str[copy_length] = '\0';
  }

  return new_str;
}

Shortread_T
Goby_read (Shortread_T *queryseq2, Gobyreader_T reader, int barcode_length,
	   bool invert_first_p, bool invert_second_p) {
#ifdef HAVE_GOBY
  unsigned long goby_read_index;
  char *acc, *read_identifier = NULL, *description = NULL;
  char *sequence1, *quality1, *sequence2, *quality2;
  int sequence1_length, quality1_length, sequence2_length, quality2_length, acc_length;

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

  *queryseq2 = Shortread_new(/*acc*/NULL,/*description*/NULL,
			     sequence2,sequence2_length,quality2,quality2_length,
			     barcode_length,invert_second_p, /*copy_acc*/false);

  return Shortread_new(acc,description,
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
#if 0
  gobyAlignments_closeIntermediateOutputFiles((*old)->helper);
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
#if 0
  gobyAlignments_openIntermediateOutputFiles(new->helper, 1);
#endif
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


#ifdef HAVE_GOBY
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
        fprintf(stderr, "::  pos=");
        for (i = 0; i < padded_length; i++) {
	  fprintf(stderr, "%d", ref_positions[i] % 10);
        }
        fprintf(stderr, "\n::   ri=");
        for (i = 0; i < padded_length; i++) {
	  fprintf(stderr, "%d", read_indexes[i] % 10);
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
	  fprintf(stderr, "%02d:%c:%02u  ", i, genomic_char, ref_positions[i]);
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
	  fprintf(stderr, "%02d:%c:%02u  ", i, read_char, read_indexes[i]);
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
    fprintf(stderr, " *** read_index [%u] or ref_position [%u] is too large! ***\n",
            read_index, ref_position);
    fprintf(stderr, ">%u\n", gobyAlEntry_getQueryIndex(writer->helper));
    if (fasta_query) {
      fprintf(stderr, "%s\n", fasta_query);
    } else {
      fprintf(stderr, "fasta_query was NULL\n");
    }
  }

  gobyAlEntry_setNumberOfMismatches(writer->helper,nmismatches);
  gobyAlEntry_setScoreInt(writer->helper,genomic_length - nindels - nmismatches);

  FREE(ref_positions);
  FREE(read_indexes);

  return;
}


static char *
merge_ref_substrings (Stage3_T stage3) {
  Substring_T substring1, substring2;
  char *genomic1, *genomic2, *result;
  int i;

  substring1 = Stage3_substring1(stage3);
  substring2 = Stage3_substring2(stage3);

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
compose_deletion (Stage3_T stage3, char **genomic, char **query, char **qual) {
  int start, i,  deletion_length, new_length, query_length, cur_pos;
  char *deletion, *new_genomic, *new_query, *new_qual;

  start = Stage3_indel_pos(stage3);
  deletion = Stage3_deletion_string(stage3);
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


static char
complement_char(char val) {
  return COMPLEMENT_LC[(int) val];
}


static void
reverse_complement(char *str, int length, bool complement) {
  bool is_odd;
  int midpoint, i;
  char temp;
  if (str == NULL || length == 0) {
    return;
  }
  is_odd = (length % 2 == 1);
  midpoint = length / 2;
  for (i = 0; i < midpoint; i++) {
    temp = (complement ? complement_char(str[i]) : str[i]);
    str[i] = (complement ? complement_char(str[length - i - 1]) : str[length - i - 1]);
    str[length - i - 1] = temp;
  }
  if (is_odd && complement) {
    str[i] = complement_char(str[i]);
  }

  return;
}

static void
output_result (Gobywriter_T writer, Stage3_T stage3, Shortread_T queryseq) {
  char *genomic, *query, *quality, *raw_genomic, *raw_query, *raw_quality, *fasta_query;
  int startpos = 0, length = 0, padding_left = 0, padding_right = 0, raw_length = 0;
  int num_deleted = 0, temp;;
  Hittype_T hittype;
  Substring_T substring1, substring2;
  bool free_sequences, reverse_strand;

  hittype = Stage3_hittype(stage3);
  reverse_strand = !Stage3_plusp(stage3);
  if (hittype == EXACT) {
    startpos = 0;
    genomic = query = raw_query = Shortread_fullpointer(queryseq);
    length = strlen(genomic);
    quality = Shortread_quality_string(queryseq);
    free_sequences = false;
  } else {
    padding_left = startpos = Substring_querystart(Stage3_substring1(stage3));
    raw_genomic = merge_ref_substrings(stage3);  // must be free'd
    fasta_query = raw_query = Shortread_fullpointer(queryseq);
    raw_quality = Shortread_quality_string(queryseq);
    length = Stage3_query_alignment_length(stage3);
    padding_right = strlen(raw_genomic) - length  - startpos;
    debug(fprintf(stderr, ":: length=%d query_al_len=%d nindels=%d\n",length,Stage3_query_alignment_length(stage3),Stage3_nindels(stage3)));

    /* Compse the deletion into the sequences */

    gobyAlignments_debugSequences(writer->helper,hittype,raw_genomic,raw_query, /*padding_left*/0, /*padding_left*/0);
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
  gobyAlignments_debugSequences(writer->helper,hittype,genomic,query, padding_left, padding_right);

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
      gobyAlignments_debugSequences(writer->helper,hittype,genomic,query,padding_left,padding_right);
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
      free(genomic);
    }
    FREE(query);
  }

  return;
}
#endif


FILE *
Goby_intermediateOutputFileHandle (Gobywriter_T writer) {
#ifdef HAVE_GOBY
  return gobyAlignments_intermediateOutputFileHandle(writer->helper);
#endif
  return NULL;
}

FILE *
Goby_intermediateIgnoredOutputFileHandle (Gobywriter_T writer) {
#ifdef HAVE_GOBY
  return gobyAlignments_intermediateIgnoredOutputFileHandle(writer->helper);
#endif
  return NULL;
}

/* Assume that stage3array has already been sorted */
void
Goby_print_single (Gobywriter_T writer, IIT_T chromosome_iit, Stage3_T *stage3array, Shortread_T queryseq1,
		   int npaths, int maxpaths, bool quiet_if_excessive_p) {
#ifdef HAVE_GOBY
  Stage3_T stage3;
  Substring_T substring1;
  UINT4 query_aligned_length;
  int pathnum, gsnap_chr_index, goby_chr_index;
  unsigned long goby_read_index;
  char *label;
  bool allocp;
  Hittype_T hittype;

  goby_read_index = (unsigned long) strtoul(Shortread_accession(queryseq1), NULL, 10);
  writer->helper->numberOfAlignedReads++;

  if (npaths > maxpaths) {
    query_aligned_length = Stage3_query_alignment_length(stage3array[0]);
    gobyAlEntry_appendTooManyHits(writer->helper,goby_read_index,query_aligned_length,npaths);
  } else {
    for (pathnum = 1; pathnum <= npaths && pathnum <= maxpaths; pathnum++) {
      stage3 = stage3array[pathnum-1];
      substring1 = Stage3_substring1(stage3);
    
      if ((hittype = Stage3_hittype(stage3)) == HALFSPLICE_DONOR || hittype == HALFSPLICE_ACCEPTOR || hittype == SPLICE ||
	  hittype == ONE_THIRD_SHORTEXON || hittype == TWO_THIRDS_SHORTEXON || hittype == SHORTEXON) {
	fprintf(stderr,"Goby does not yet support splicing\n");
      } else {
	gobyAlignments_appendEntry(writer->helper);
	gobyAlEntry_setMultiplicity(writer->helper,1);
	gobyAlEntry_setQueryIndex(writer->helper,goby_read_index);

	/* Convert the gsnap chrom number to the goby chrom number */
	gsnap_chr_index = Stage3_chrnum(stage3);
	gobyAlEntry_setTargetIndex(writer->helper,gsnap_chr_index - 1);
	gobyAlEntry_setPosition(writer->helper,Stage3_chrpos_low_trim(stage3));
	gobyAlEntry_setMatchingReverseStrand(writer->helper,(Stage3_plusp(stage3) == 0 ? 1 : 0));
	gobyAlEntry_setQueryPosition(writer->helper,Substring_querystart(substring1));
	gobyAlEntry_setScoreInt(writer->helper,Substring_querylength(substring1) - Stage3_score(stage3));
	gobyAlEntry_setNumberOfMismatches(writer->helper,Stage3_nmismatches_refdiff(stage3));
	gobyAlEntry_setNumberOfIndels(writer->helper,Stage3_nindels(stage3));
	gobyAlEntry_setQueryAlignedLength(writer->helper,(UINT4) Stage3_query_alignment_length(stage3));
	gobyAlEntry_setTargetAlignedLength(writer->helper,Stage3_genomic_alignment_length(stage3));
	gobyAlEntry_setQueryLength(writer->helper,Substring_querylength(substring1));
      
	output_result(writer,stage3,queryseq1);
      }
    }
  }

#endif
  return;
}

/* Assumption: mate will be NULL if the second half of the pair doesn't align */
void
Goby_print_pair (Gobywriter_T writer, Stage3_T this, Stage3_T mate, char *acc, int pathnum, int npaths,
		 int mapq_score, IIT_T chromosome_iit, Shortread_T queryseq,
		 Shortread_T queryseq_mate, int pairedlength, Resulttype_T resulttype,
		 bool first_read_p, int npaths_mate, int quality_shift,
		 char *sam_read_group_id, bool invertp, bool invert_mate_p) {
#ifdef HAVE_GOBY
  Substring_T substring1;
  UINT4 query_aligned_length;
  unsigned long goby_read_index;
  unsigned int sam_flags = 0U;
  unsigned int fragment_index = 0;
  unsigned int m_fragment_index = 0;
  int gsnap_chr_index, goby_chr_index;
  bool allocp;
  char *label;

  substring1 = Stage3_substring1(this);
  sam_flags = SAM_compute_flag(substring1,mate,resulttype,first_read_p,
			       pathnum,npaths,npaths_mate,invertp,invert_mate_p);

  if (sam_flags & FIRST_READ_P) {
    fragment_index = 0;
    m_fragment_index = 1;
  } else if (sam_flags & SECOND_READ_P) {
    fragment_index = 1;
    m_fragment_index = 0;
  }

  goby_read_index = (unsigned long) strtoul(acc, NULL, 10);
  writer->helper->numberOfAlignedReads++;


  if (Stage3_hittype(this) == SPLICE) {
    fprintf(stderr,"Goby does not yet support hittype of SPLICE\n");
  } else {
    gobyAlignments_appendEntry(writer->helper);
    gobyAlEntry_setMultiplicity(writer->helper,1);
    gobyAlEntry_setQueryIndex(writer->helper,goby_read_index);

    /* Convert the gsnap chrom number to the goby chrom number */
    gsnap_chr_index = Stage3_chrnum(this);
    label = IIT_label(chromosome_iit,gsnap_chr_index,&allocp);
    goby_chr_index = atoi(label);
    if (allocp) {
      FREE(label);
    }

    gobyAlEntry_setTargetIndex(writer->helper,goby_chr_index);
    gobyAlEntry_setPosition(writer->helper,Stage3_chrpos_low_trim(this));
    gobyAlEntry_setMatchingReverseStrand(writer->helper,(Stage3_plusp(this) == 0 ? 1 : 0));
    gobyAlEntry_setQueryPosition(writer->helper,Substring_querystart(substring1));
    gobyAlEntry_setScoreInt(writer->helper,Substring_querylength(substring1) - Stage3_score(this));
    gobyAlEntry_setNumberOfMismatches(writer->helper,Stage3_nmismatches_refdiff(this));
    gobyAlEntry_setNumberOfIndels(writer->helper,Stage3_nindels(this));
    gobyAlEntry_setQueryAlignedLength(writer->helper,(UINT4) Stage3_query_alignment_length(this));
    gobyAlEntry_setTargetAlignedLength(writer->helper,Stage3_genomic_alignment_length(this));
    gobyAlEntry_setQueryLength(writer->helper,Substring_querylength(substring1));
    gobyAlEntry_setFragmentIndex(writer->helper,fragment_index);

    gobyAlEntry_setPairFlags(writer->helper,sam_flags);
    if (mate) {
      gobyAlEntry_setPairFragmentIndex(writer->helper,m_fragment_index);
      gobyAlEntry_setPairTargetIndex(writer->helper,Stage3_chrnum(mate) - 1);
      gobyAlEntry_setPairPosition(writer->helper,Stage3_chrpos_low_trim(mate));
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
  Stage3_T *stage3array1, *stage3array2, stage3, mate;
  int npaths, npaths1, npaths2, pathnum;
  char *acc;
  UINT4 query_aligned_length;
  unsigned long goby_read_index;
    
  acc = Shortread_accession(queryseq1);

  if (resulttype == CONCORDANT_UNIQ) {
    stage3pairarray = (Stage3pair_T *) Result_array(&npaths,result);
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
    stage3pairarray = (Stage3pair_T *) Result_array(&npaths,result);

    if (npaths > maxpaths) {
      /* No output if excessive for Gobyweb, but output TMH. */
      /* TODO: Q: Should we be outputting BOTH primary and mate TMH? */
      goby_read_index = (unsigned long) strtoul(Shortread_accession(queryseq1), NULL, 10);
      query_aligned_length = Stage3_query_alignment_length(Stage3pair_hit5(stage3pair));
      gobyAlEntry_appendTooManyHits(writer->helper,goby_read_index,query_aligned_length,npaths);

      goby_read_index = (unsigned long) strtoul(Shortread_accession(queryseq2), NULL, 10);
      query_aligned_length = Stage3_query_alignment_length(Stage3pair_hit3(stage3pair));
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

  } else if (resulttype == UNPAIRED_UNIQ) {
    /* Should print mate information in this situation */
    stage3array1 = (Stage3_T *) Result_array(&npaths1,result);
    stage3array2 = (Stage3_T *) Result_array2(&npaths2,result);

    /* print first end */
    /* Stage3_eval_and_sort(stage3array1,npaths1,maxpaths,queryseq1); */
    Goby_print_pair(writer,stage3array1[0],/*mate*/stage3array2[0],acc,
		    /*pathnum*/1,/*npaths*/1,Stage3_mapq_score(stage3array1[0]),
		    chromosome_iit,/*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
		    /*pairedlength*/0U,resulttype,
		    /*first_read_p*/true,/*npaths_mate*/1,quality_shift,sam_read_group_id,
		    invert_first_p,invert_second_p);

    /* print second end */
    /* Stage3_eval_and_sort(stage3array2,npaths2,maxpaths,queryseq2); */
    Goby_print_pair(writer,stage3array2[0],/*mate*/stage3array1[0],acc,
		    /*pathnum*/1,/*npaths*/1,Stage3_mapq_score(stage3array2[0]),
		    chromosome_iit,/*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
		    /*pairedlength*/0U,resulttype,
		    /*first_read_p*/false,/*npaths_mate*/1,quality_shift,sam_read_group_id,
		    invert_second_p,invert_first_p);

  } else if (resulttype == UNPAIRED_MULT) {
    stage3array1 = (Stage3_T *) Result_array(&npaths1,result);
    stage3array2 = (Stage3_T *) Result_array2(&npaths2,result);

#if 0
    /* Do eval and sorting first */
    if (npaths1 == 1) {
      /* Stage3_eval_and_sort(stage3array1,npaths1,maxpaths,queryseq1); */
    } else if (npaths1 > maxpaths) {
      /* Don't sort */
    } else {
      /* Stage3_eval_and_sort(stage3array1,npaths1,maxpaths,queryseq1); */
    }
#endif

#if 0
    if (npaths2 == 1) {
      /* Stage3_eval_and_sort(stage3array2,npaths2,maxpaths,queryseq2); */
    } else if (npaths2 > maxpaths) {
      /* Don't sort */
    } else {
      /* Stage3_eval_and_sort(stage3array2,npaths2,maxpaths,queryseq2); */
    }
#endif

    /* print first end results */
    mate = (npaths2 == 0) ? (Stage3_T) NULL : stage3array2[0];

    if (npaths1 == 1) {
      stage3 = stage3array1[0];
      Goby_print_pair(writer,stage3,mate,acc,/*pathnum*/1,npaths1,Stage3_mapq_score(stage3),
		      chromosome_iit,/*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
		      /*pairedlength*/0U,resulttype,/*first_read_p*/true,
                      /*npaths_mate*/npaths2,quality_shift,sam_read_group_id,
                      invert_first_p,invert_second_p);

    } else if (npaths1 > maxpaths) {
      /** No output if excessive for Gobyweb, but output TMH. */
      goby_read_index = (unsigned long) strtoul(Shortread_accession(queryseq1), NULL, 10);
      query_aligned_length = Stage3_query_alignment_length(stage3array1[0]);
      gobyAlEntry_appendTooManyHits(writer->helper,goby_read_index,query_aligned_length,npaths1);

    } else {
      for (pathnum = 1; pathnum <= npaths1 && pathnum <= maxpaths; pathnum++) {
	stage3 = stage3array1[pathnum-1];
	Goby_print_pair(writer,stage3,mate,acc,pathnum,npaths1,Stage3_mapq_score(stage3),
                        chromosome_iit,/*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
                        /*pairedlength*/0U,resulttype,/*first_read_p*/true,/*npaths_mate*/npaths2,
                        quality_shift,sam_read_group_id,invert_first_p,invert_second_p);
      }
    }

    /* print second end results */
    mate = (npaths1 == 0) ? (Stage3_T) NULL : stage3array1[0];

    if (npaths2 == 1) {
      stage3 = stage3array2[0];
      Goby_print_pair(writer,stage3,mate,acc,/*pathnum*/1,npaths2,Stage3_mapq_score(stage3),
		      chromosome_iit,/*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
		      /*pairedlength*/0U,resulttype,/*first_read_p*/false,/*npaths_mate*/npaths1,
		      quality_shift,sam_read_group_id,invert_second_p,invert_first_p);

    } else if (npaths2 > maxpaths) {
      /** No output if excessive for Gobyweb, but output TMH. */
      goby_read_index = (unsigned long) strtoul(Shortread_accession(queryseq2), NULL, 10);
      query_aligned_length = Stage3_query_alignment_length(stage3array2[0]);
      gobyAlEntry_appendTooManyHits(writer->helper,goby_read_index,query_aligned_length,npaths2);

    } else {
      for (pathnum = 1; pathnum <= npaths2 && pathnum <= maxpaths; pathnum++) {
	stage3 = stage3array2[pathnum-1];
	Goby_print_pair(writer,stage3,mate,acc,pathnum,npaths2,Stage3_mapq_score(stage3),
                        chromosome_iit,/*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
                        /*pairedlength*/0U,resulttype,/*first_read_p*/false,
                        /*npaths_mate*/npaths1,quality_shift,sam_read_group_id,
                        invert_second_p,invert_first_p);
      }
    }

  } else {
    if (resulttype == HALFMAPPING_UNIQ || resulttype == HALFMAPPING_MULT) {
      /* These are the last two resulttypes we can deal with */
    } else {
      abort();
    }

    stage3array1 = (Stage3_T *) Result_array(&npaths1,result);
    stage3array2 = (Stage3_T *) Result_array2(&npaths2,result);

#if 0
    /* Do eval and sorting first */
    if (npaths1 == 0) {
      /* Nothing to sort */
    } else if (npaths1 == 1) {
      /* Stage3_eval_and_sort(stage3array1,npaths1,maxpaths,queryseq1); */
    } else if (npaths1 > maxpaths) {
      /* Don't sort */
    } else {
      /* Stage3_eval_and_sort(stage3array1,npaths1,maxpaths,queryseq1); */
    }
#endif

#if 0
    if (npaths2 == 0) {
      /* Nothing to sort */
    } else if (npaths2 == 1) {
      /* Stage3_eval_and_sort(stage3array2,npaths2,maxpaths,queryseq2); */
    } else if (npaths2 > maxpaths) {
      /* Don't sort */
    } else {
      /* Stage3_eval_and_sort(stage3array2,npaths2,maxpaths,queryseq2); */
    }
#endif

    /* print first end results */
    mate = (npaths2 == 0) ? (Stage3_T) NULL : stage3array2[0];

    if (npaths1 == 0) {
      /** No output for Gobyweb. */

    } else if (npaths1 == 1) {
      /* mate should be NULL here */

      stage3 = stage3array1[0];
      Goby_print_pair(writer,stage3,mate,acc,/*pathnum*/1,npaths1,Stage3_mapq_score(stage3),
		      chromosome_iit,/*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
		      /*pairedlength*/0U,resulttype,/*first_read_p*/true,
                      /*npaths_mate*/npaths2, quality_shift,sam_read_group_id,
                      invert_first_p,invert_second_p);

    } else if (npaths1 > maxpaths) {
      /** No output if excessive for Gobyweb, but output TMH. */
      goby_read_index = (unsigned long) strtoul(Shortread_accession(queryseq1), NULL, 10);
      query_aligned_length = Stage3_query_alignment_length(stage3array1[0]);
      gobyAlEntry_appendTooManyHits(writer->helper,goby_read_index,query_aligned_length,npaths1);

    } else {
      /* mate should be NULL here */
      for (pathnum = 1; pathnum <= npaths1 && pathnum <= maxpaths; pathnum++) {
	stage3 = stage3array1[pathnum-1];
	Goby_print_pair(writer,stage3,mate,acc,pathnum,npaths1,Stage3_mapq_score(stage3),
                        chromosome_iit,/*queryseq*/queryseq1,/*queryseq_mate*/queryseq2,
                        /*pairedlength*/0U,resulttype,/*first_read_p*/true,
                        /*npaths_mate*/npaths2,quality_shift,sam_read_group_id,
                        invert_first_p,invert_second_p);
      }
    }

    /* print second end results */
    mate = (npaths1 == 0) ? (Stage3_T) NULL : stage3array1[0];

    if (npaths2 == 0) {
      /* No output for Gobyweb. */

    } else if (npaths2 == 1) {
      /* mate should be NULL here */

      stage3 = stage3array2[0];
      Goby_print_pair(writer,stage3,mate,acc,/*pathnum*/1,npaths2,Stage3_mapq_score(stage3),
		      chromosome_iit,/*queryseq*/queryseq2,/*queryseq_mate*/queryseq1,
		      /*pairedlength*/0U,resulttype,/*first_read_p*/false,
                      /*npaths_mate*/npaths1, quality_shift,sam_read_group_id,
                      invert_second_p,invert_first_p);

    } else if (npaths2 > maxpaths) {
      /* No output if excessive for Gobyweb, but output TMH. */
      goby_read_index = (unsigned long) strtoul(Shortread_accession(queryseq2), NULL, 10);
      query_aligned_length = Stage3_query_alignment_length(stage3array2[0]);
      gobyAlEntry_appendTooManyHits(writer->helper,goby_read_index,query_aligned_length,npaths2);

    } else {
      /* mate should be NULL here */
      for (pathnum = 1; pathnum <= npaths2 && pathnum <= maxpaths; pathnum++) {
	stage3 = stage3array2[pathnum-1];
	Goby_print_pair(writer,stage3,mate,acc,pathnum,npaths2,Stage3_mapq_score(stage3),
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
