static char rcsid[] = "$Id: get-genome.c 85655 2013-02-07 17:26:37Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>		/* For rindex */
#include <ctype.h>

#include "bool.h"
#include "mem.h"
#include "access.h"
#include "genomicpos.h"
#include "sequence.h"
#include "chrom.h"
#include "genome.h"
#include "iit-read.h"
#include "chrsubset.h"
#include "datadir.h"
#include "parserange.h"
#include "separator.h"
#include "complement.h"
#include "getopt.h"


#define INFINITY -1


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


/* Program Options */
#if 0
static bool replacex = false;
#endif
static int print_snps_mode = 3;	/* 1 = alt version only, 2 = snps marked, 3 = both */
static char *snps_root = NULL;

static char *user_genomedir = NULL;
static char *user_snpsdir = NULL;

static bool uncompressedp = false;
static bool uppercasep = false;
static bool revcomp = false;
static bool coordp = false;
static char *genomesubdir = NULL;
static char *dbroot = NULL;
static char *dbversion = NULL;
static char *user_typestring = NULL;
static int wraplength = 60;
static char *header = NULL;
static bool levelsp = false;
static bool rawp = false;

static char *user_mapdir = NULL;
static char *map_iitfile = NULL;
static int nflanking = 0;
static bool exonsp = false;
static bool sequencep = false;
static bool signedp = false;

static bool sortp = false;


/* Dump options */
static bool dumpallp = false;
static bool dumpchrp = false;
static bool dumpchr_forsam_p = false;
static bool dumpsegsp = false;
static bool dumpchrsubsetsp = false;

static struct option long_options[] = {
  /* Input options */
  {"dir", required_argument, 0, 'D'}, /* user_genomedir */
  {"db", required_argument, 0, 'd'}, /* dbroot */

  /* Output options */
  {"coords", no_argument, 0, 'C'}, /* coordp */
  {"uppercase", no_argument, 0, 'U'}, /* uppercasep */
  {"wraplength", required_argument, 0, 'l'}, /* wraplength */
  {"fullgenome", no_argument, 0, 'G'}, /* uncompressedp */
  {"header", required_argument, 0, 'h'}, /* header */
  {"snpsdir", required_argument, 0, 'V'}, /* user_snpsdir */
  {"usesnps", required_argument, 0, 'v'}, /* snps_root */
  {"snpformat", required_argument, 0, 'f'}, /* print_snps */

  /* External map options */
  {"mapdir", required_argument, 0, 'M'}, /* user_mapdir */
  {"map", required_argument, 0, 'm'},	/* map_iitfile */
  {"ranks", no_argument, 0, 'k'}, /* levelsp */
  {"raw", no_argument, 0, 'r'}, /* rawp */
  {"flanking", required_argument, 0, 'u'}, /* nflanking */
  {"exons", no_argument, 0, 'E'},	   /* exonsp */
  {"sequence", no_argument, 0, 'S'},	   /* sequencep */
  {"signed", no_argument, 0, 's'},	   /* signedp */

  /* Dump options */
  {"dump", no_argument, 0, 'A'},	/* dumpallp */
  {"chromosomes", no_argument, 0, 'L'},	/* dumpchrp */
  {"forsam", no_argument, 0, 0},	/* dumpchr_forsam_p */
  {"contigs", no_argument, 0, 'I'}, /* dumpsegsp */
  {"chrsubsets", no_argument, 0, 'c'}, /* dumpchrsubsetsp */
  
  /* Help options */
  {"version", no_argument, 0, '^'}, /* print_program_version */
  {"help", no_argument, 0, '?'}, /* print_program_usage */
  {0, 0, 0, 0}
};


static void
print_program_version () {
  fprintf(stdout,"\n");
  fprintf(stdout,"get-genome: retrieval utility for genomic segments\n");
  fprintf(stdout,"Part of GMAP package, version %s\n",PACKAGE_VERSION);
  fprintf(stdout,"Thomas D. Wu and Colin K. Watanabe, Genentech, Inc.\n");
  fprintf(stdout,"Contact: twu@gene.com\n");
  fprintf(stdout,"\n");
  return;
}

static void
print_program_usage () {
  fprintf(stdout,"\
Usage: get-genome [OPTIONS...] -d genome [genome:]range, or\n\
       get-genome [OPTIONS...] -d genome chromosome:range, or\n\
       get-genome [OPTIONS...] -d genome contig[:range]\n\
where\n\
   range is startposition..endposition (endpos < startpos means - strand)\n\
         or startposition+length (+ strand)\n\
         or startposition+-length (- strand)\n\
\n\
Input options\n\
  -D, --dir=STRING        Genome directory\n\
  -d, --db=STRING         Genome database.  If argument is '?' (with\n\
                            the quotes), this command lists available databases.\n\
\n\
Output options\n\
  -2, --dibase            Use dibase version of genome\n\
  -C, --coords            Show coordinates only\n\
  -U, --uppercase         Convert sequence to uppercase\n\
  -l, --wraplength=INT    Wrap length for sequence (default=60)\n\
  -G, --fullgenome        Use full (uncompressed) version of genome\n\
  -h, --header=STRING     Desired header line\n\
  -V, --snpsdir=STRING    Directory for SNPs index files (created using snpindex)\n\
  -v, --usesnps=STRING    Use snp version (built by snpindex)\n\
  -f, --snpformat=INT     Print snp information from database built previously\n\
                            using snpindex (0=none, 1=alternate version only\n\
                            2=both versions merged (using N), 3=both versions separate (default)\n\
\n\
External map file options\n\
  -M, --mapdir=directory  Map directory\n\
  -m, --map=iitfile       Map file.  If argument is '?' (with the quotes),\n\
                            this lists available map files.\n\
  -S, --sequence          For a gene map file, prints the sequence\n\
  -E, --exons             For a gene map file, prints the sequence, one exon per line\n\
  -k, --ranks             Prints levels for non-overlapping printing of map hits\n\
  -r, --raw               Prints sequence as ASCII numeric codes\n\
  -u, --flanking=INT      Show flanking hits (default 0)\n\
  -t, --maptype=STRING    Show only intervals with given type\n\
  -s, --signed            Show only intervals with same direction as query.  If flanking hits\n\
                            are also requested, show only flanking hits downstream in direction of\n\
                            query.\n\
\n\
Dump options\n\
  -A, --dump              Dump entire genome in FASTA format\n\
  -L, --chromosomes       List all chromosomes with universal coordinates\n\
  --forsam                List all chromosomes for use in a SAM file\n\
  -I, --contigs           List all contigs with universal coordinates\n\
\n\
Help options\n\
  -^, --version           Show version\n\
  -?, --help              Show this help message\n\
");
  return;
}

/* Parsing functions */

/* Note that isnumber is a function in ctype.h on some systems */
static bool
isnumberp (unsigned int *result, char *string) {
  char *p = string;

  *result = 0U;
  while (*p != '\0') {
    if (*p == ',') {
      /* Skip commas */
    } else if (!isdigit((int) *p)) {
      return false;
    } else {
      *result = (*result) * 10 + (*p - '0');
    }
    p++;
  }
  return true;
}

/* Returns coordinates as zero-based */
static bool
isrange (unsigned int *left, unsigned int *length, bool *revcomp, char *string) {
  bool result;
  char *copy, *startstring, *endstring;
  unsigned int start, end;

  copy = (char *) CALLOC(strlen(string)+1,sizeof(char));
  strcpy(copy,string);

  if (index(copy,'.')) {
    startstring = strtok(copy,"..");
    endstring = strtok(NULL,"..");
    if (!isnumberp(&start,startstring) || !isnumberp(&end,endstring)) {
      result = false;
    } else if (start <= end) {
      *length = end - start + 1;
      *left = start - 1;
      *revcomp = false;
      debug(printf("(..) "));
      result = true;
    } else {
      *length = start - end + 1;
      *left = end - 1;
      *revcomp = true;
      debug(printf("(..) "));
      result = true;
    }

  } else if (index(copy,'+')) {
    startstring = strtok(copy,"+");
    endstring = strtok(NULL,"+");
    if (!isnumberp(&start,startstring)) {
      result = false;
    } else if (endstring[0] == '-' && isnumberp(&(*length),&(endstring[1]))) {
      *left = start - (*length);
      *revcomp = true;
      debug(printf("(-) "));
      result = true;
    } else if (!isnumberp(&(*length),endstring)) {
      result = false;
    } else {
      *left = start - 1;
      *revcomp = false;
      debug(printf("(+) "));
      result = true;
    }

  } else if (index(copy,'-')) {
    /* Old notation */
    startstring = strtok(copy,"--");
    endstring = strtok(NULL,"--");
    if (!isnumberp(&start,startstring) || !isnumberp(&end,endstring)) {
      result = false;
    } else if (start <= end) {
      *length = end - start + 1;
      *left = start - 1;
      *revcomp = false;
      debug(printf("(--) "));
      result = true;
    } else {
      *length = start - end + 1;
      *left = end - 1;
      *revcomp = true;
      debug(printf("(--) "));
      result = true;
    }

    /* Don't allow this yet ...
  } else if (index(copy,'-')) {
    startstring = strtok(copy,"-");
    endstring = strtok(NULL,"-");
    if (!isnumberp(&start,startstring) || !isnumberp(&end,endstring)) {
      result = false;
    } else if (end > start - 1) {
      result = false;
    } else {
      *left = start - 1 - end;
      *length = end;
      *revcomp = true;
      result = true;
    }
    */
    
  } else {
    result = false;
  }

  FREE(copy);
  return result;
}
  

/* Printing functions */

static void
print_two_coords (unsigned int left, unsigned int length, IIT_T chromosome_iit) {
  char *chromosome;
  unsigned int chrpos;

  printf("%u%s%u\t",left+1,SEPARATOR,left+length);
  chromosome = IIT_string_from_position(&chrpos,left,chromosome_iit);
  printf("%s:%u\t",chromosome,chrpos+1U);
  FREE(chromosome);

  chromosome = IIT_string_from_position(&chrpos,left+length-1,chromosome_iit);
  printf("%s:%u\n",chromosome,chrpos+1U);
  FREE(chromosome);
  return;
}


/* Retrieval functions */

static int
translate_chromosomepos_universal (unsigned int *genomicstart, unsigned int *genomiclength, 
				   char *chromosome, unsigned int left, unsigned int length,
				   IIT_T chromosome_iit) {
  int rc = 1, index;
  Interval_T interval;
  
  if ((index = IIT_find_linear(chromosome_iit,chromosome)) >= 0) {
    interval = IIT_interval(chromosome_iit,index);
    *genomicstart = Interval_low(interval)+left;
    if (length == 0) {
      *genomiclength = Interval_length(interval)-left;
    } else {
      *genomiclength = length;
    }
    rc = 0;
  }
  
  return rc;
}

static int
translate_chromosomepos_segment (unsigned int *segmentstart, unsigned int *segmentlength, 
				 char *chromosome, unsigned int left, unsigned int length,
				 IIT_T chromosome_iit) {
  int rc = 1, index;
  Interval_T interval;
  
  if ((index = IIT_find_linear(chromosome_iit,chromosome)) >= 0) {
    interval = IIT_interval(chromosome_iit,index);
    *segmentstart = left;
    if (length == 0) {
      *segmentlength = Interval_length(interval)-left;
    } else {
      *segmentlength = length;
    }
    rc = 0;
  }
  
  return rc;
}


static int
translate_contig (unsigned int *genomicstart, unsigned int *genomiclength,
		  char *contig, unsigned int left, unsigned int length, IIT_T contig_iit) {
  int rc = 1, index;
  Interval_T interval;
  
  if ((index = IIT_find_one(contig_iit,contig)) >= 0) {
    interval = IIT_interval(contig_iit,index);
    *genomicstart = Interval_low(interval)+left;
    if (length == 0) {
      *genomiclength = Interval_length(interval)-left;
    } else {
      *genomiclength = length;
    }
    rc = 0;
  }

  return rc;
}


#if 0
static bool
parse_query_universal (char **segment, Genomicpos_T *genomicstart, Genomicpos_T *genomiclength, bool *revcomp,
		       char *query, char *genomesubdir, char *fileroot, char *version) {
  char *coords, *filename;
  Genomicpos_T result, left, length;
  IIT_T chromosome_iit, contig_iit;
  int rc;
  
  *revcomp = false;
  if (index(query,':')) {
    /* Segment must be a genome, chromosome, or contig */
    debug(printf("Parsed query %s into ",query));
    *segment = strtok(query,":");
    coords = strtok(NULL,":");
    debug(printf("segment %s and coords %s\n",*segment,coords));

    /* Try chromosome first */
    filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+strlen(fileroot)+
			       strlen(".chromosome.iit")+1,sizeof(char));
    sprintf(filename,"%s/%s.chromosome.iit",genomesubdir,fileroot);
    chromosome_iit = IIT_read(filename,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
			      /*divstring*/NULL,/*add_iit_p*/false,/*labels_read_p*/true);
    FREE(filename);

    debug(printf("Interpreting segment %s as a chromosome\n",*segment));
    if (coords == NULL) {
      debug(printf("  entire chromosome\n"));
      rc = translate_chromosomepos_universal(&(*genomicstart),&(*genomiclength),*segment,left=0,length=0,chromosome_iit);
    } else if (isnumberp(&result,coords)) {
      debug(printf("  and coords %s as a number\n",coords));
      rc = translate_chromosomepos_universal(&(*genomicstart),&(*genomiclength),*segment,left=result-1,length=1,chromosome_iit);
    } else if (isrange(&left,&length,&(*revcomp),coords)) {
      debug(printf("  and coords %s as a range starting at %u with length %u and revcomp = %d\n",
		   coords,left,length,*revcomp));
      rc = translate_chromosomepos_universal(&(*genomicstart),&(*genomiclength),*segment,left,length,chromosome_iit);
    } else {
      debug(printf("  but coords %s is neither a number nor a range\n",coords));
      rc = -1;
    }

    IIT_free(&chromosome_iit);

    if (rc != 0) {
      /* Try contig */
      filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+strlen(fileroot)+
				 strlen(".contig.iit")+1,sizeof(char));
      sprintf(filename,"%s/%s.contig.iit",genomesubdir,fileroot);
      contig_iit = IIT_read(filename,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
			    /*divstring*/NULL,/*add_iit_p*/false,/*labels_read_p*/true);
      FREE(filename);

      debug(printf("Interpreting segment %s as a contig\n",*segment));
      if (isnumberp(&result,coords)) {
	debug(printf("  and coords %s as a number\n",coords));
	rc = translate_contig(&(*genomicstart),&(*genomiclength),*segment,left=result-1,length=1,contig_iit);
      } else if (isrange(&left,&length,&(*revcomp),coords)) {
	debug(printf("  and coords %s as a range starting at %u with length %u and revcomp = %d\n",
		     coords,left,length,*revcomp));
	rc = translate_contig(&(*genomicstart),&(*genomiclength),*segment,left,length,contig_iit);
      } else {
	debug(printf("  but coords %s is neither a number nor a range\n",coords));
	rc = -1;
      }

      IIT_free(&contig_iit);
    }

    if (rc != 0) {
      fprintf(stderr,"Can't find coordinates %s:%s\n",*segment,coords);
      return false;
    } else {
      return true;
    }

  } else {
    /* Query must be a genomic position, genomic range, or contig */
    debug(printf("Parsed query %s as atomic ",query));
    if (isnumberp(&result,query)) {
      debug(printf("number\n"));
      *genomicstart = result-1;
      *genomiclength = 1;
      return true;
    } else if (isrange(&left,&length,&(*revcomp),query)) {
      debug(printf("range\n"));
      *genomicstart = left;
      *genomiclength = length;
      return true;
    } else {

      debug(printf("contig\n"));
#if 0
      filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+strlen(fileroot)+
				 strlen(".contig.iit")+1,sizeof(char));
      sprintf(filename,"%s/%s.contig.iit",genomesubdir,fileroot);
      contig_iit = IIT_read(filename,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
			    /*divstring*/NULL,/*add_iit_p*/false,/*labels_read_p*/true);
      FREE(filename);

      rc = translate_contig(&(*genomicstart),&(*genomiclength),query,left=0,length=0,contig_iit);
      IIT_free(&contig_iit);
#endif
      return false;
    }
  }
}
#endif


static bool
parse_query (char **divstring, unsigned int *coordstart, unsigned int *coordend, bool *revcomp,
	     char *query_original, char *filename) {
  char *coords;
  unsigned int result, left, length;
  char *query, *temp;
  
  *revcomp = false;
  if (index(query_original,':')) {
    /* Query may have a div */
    query = (char *) CALLOC(strlen(query_original)+1,sizeof(char));
    strcpy(query,query_original);

    debug(printf("Parsed query %s into ",query));
    temp = strtok(query,":");

    *divstring = (char *) CALLOC(strlen(temp)+1,sizeof(char));
    strcpy(*divstring,temp);

    coords = strtok(NULL,":");
    debug(printf("divstring %s and coords %s\n",*divstring,coords));

    if (IIT_read_divint(filename,*divstring,/*add_iit_p*/true) < 0) {
      debug(printf("  but divstring not found in %s, so treat as label\n",filename));
      FREE(*divstring);
      FREE(query);
      return false;
    } else if (coords == NULL) {
      debug(printf("  entire div\n"));
      *coordstart = 0;
      *coordend = INFINITY;
      FREE(query);
      return true;
    } else if (isnumberp(&result,coords)) {
      debug(printf("  and coords %s as a number\n",coords));
      *coordstart = result;
      *coordend = result;
      FREE(query);
      return true;
    } else if (isrange(&left,&length,&(*revcomp),coords)) {
      debug(printf("  and coords %s as a range starting at %u with length %u and revcomp = %d\n",
		   coords,left,length,*revcomp));
      *coordstart = left;
      *coordend = left + length;
      FREE(query);
      return true;
    } else {
      debug(printf("  but coords %s is neither a number nor a range.  Interpret as a label.\n",coords));
      FREE(*divstring);
      FREE(query);
      return false;
    }

  } else {
    /* No div.  Query must be a number, range, or label */
    query = query_original;

    *divstring = NULL;
    debug(printf("Parsed query %s without a div ",query));
    if (isnumberp(&result,query)) {
      debug(printf("number\n"));
      *coordstart = result;
      *coordend = result;
      return true;
    } else if (isrange(&left,&length,&(*revcomp),query)) {
      debug(printf("range\n"));
      *coordstart = left;
      *coordend = left + length;
      return true;
    } else {
      debug(printf("label\n"));
      return false;
    }
  }
}


/* This code is duplicated in gmap.c */
static IIT_T global_altstrain_iit;

static int
index_compare (const void *a, const void *b) {
  int index1 = * (int *) a;
  int index2 = * (int *) b;
  int type1, type2;
  Genomicpos_T pos1, pos2;

  type1 = Interval_type(IIT_interval(global_altstrain_iit,index1));
  type2 = Interval_type(IIT_interval(global_altstrain_iit,index2));
  
  if (type1 < type2) {
    return -1;
  } else if (type1 > type2) {
    return +1;
  } else {
    /* Store in descending genomic position, so right shifting works
       in Genome_patch_strain */
    pos1 = Interval_low(IIT_interval(global_altstrain_iit,index1));
    pos2 = Interval_low(IIT_interval(global_altstrain_iit,index2));

    if (pos1 > pos2) {
      return -1;
    } else if (pos1 < pos2) {
      return +1;
    } else {
      return 0;
    }
  }
}


static void
print_sequence (Genome_T genome, Genome_T genomealt, Genomicpos_T genomicstart, Genomicpos_T genomiclength,
		IIT_T chromosome_iit, bool whole_chromosome_p) {
  char *chromosome1, *chromosome2;
  Sequence_T genomicseg, genomicseg_alt, genomicseg_snp;
  Genomicpos_T chrpos;

  /* Handle reference strain */
  if (user_typestring != NULL) {
    /* Don't print a header */
  } else if (header != NULL) {
    printf(">%s\n",header);
  } else if (revcomp == true) {
    printf(">");
    chromosome1 = IIT_string_from_position(&chrpos,genomicstart+genomiclength-1,chromosome_iit);
    if (whole_chromosome_p == true) {
      printf("%s\n",chromosome1);
    } else {
      printf("%s:%u",chromosome1,chrpos+1U);
      printf("%s",SEPARATOR);
      chromosome2 = IIT_string_from_position(&chrpos,genomicstart+1,chromosome_iit);
      if (strcmp(chromosome2,chromosome1) == 0) {
	printf("%u",chrpos);
      } else {
	printf("%s:%u",chromosome2,chrpos);
      }
      printf(" %s:%u%s%u\n",dbversion,genomicstart+genomiclength,SEPARATOR,genomicstart+1U);
      FREE(chromosome2);
    }
    FREE(chromosome1);
  } else {
    printf(">");
    chromosome1 = IIT_string_from_position(&chrpos,genomicstart,chromosome_iit);
    if (whole_chromosome_p == true) {
      printf("%s\n",chromosome1);
    } else {
      printf("%s:%u",chromosome1,chrpos+1U);
      printf("%s",SEPARATOR);
      chromosome2 = IIT_string_from_position(&chrpos,genomicstart+genomiclength-1,chromosome_iit);
      if (strcmp(chromosome2,chromosome1) == 0) {
	printf("%u",chrpos+1U);
      } else {
	printf("%s:%u",chromosome2,chrpos+1U);
      }
      printf(" %s:%u%s%u\n",dbversion,genomicstart+1U,SEPARATOR,genomicstart+genomiclength);
      FREE(chromosome2);
    }
    FREE(chromosome1);
  }

  if (snps_root == NULL || print_snps_mode == 0 || print_snps_mode == 2) {
    genomicseg = Genome_get_segment(genome,genomicstart,genomiclength,chromosome_iit,revcomp);
    if (user_typestring == NULL) {
      if (rawp == true) {
	Sequence_print_raw(genomicseg);
      } else {
	Sequence_print(stdout,genomicseg,uppercasep,wraplength,/*trimmedp*/false);
      }
      Sequence_free(&genomicseg);
    }

  } else if (print_snps_mode == 1) {
    /* Handle both reference and alternate versions */
    genomicseg = Genome_get_segment(genome,genomicstart,genomiclength,chromosome_iit,revcomp);
    genomicseg_alt = Genome_get_segment_alt(genomealt,genomicstart,genomiclength,chromosome_iit,revcomp);
    genomicseg_snp = Genome_get_segment_snp(genomealt,genomicstart,genomiclength,chromosome_iit,revcomp);
    if (user_typestring == NULL) {
      if (rawp == true) {
	Sequence_print_raw(genomicseg);
      } else {
	Sequence_print_alt(genomicseg,genomicseg_alt,genomicseg_snp,uppercasep,wraplength);
      }
      Sequence_free(&genomicseg_snp);
      Sequence_free(&genomicseg_alt);
      Sequence_free(&genomicseg);
    }

  } else {
    /* Handle both reference and alternate versions */
    genomicseg = Genome_get_segment(genome,genomicstart,genomiclength,chromosome_iit,revcomp);
    genomicseg_snp = Genome_get_segment_snp(genomealt,genomicstart,genomiclength,chromosome_iit,revcomp);
    if (user_typestring == NULL) {
      if (rawp == true) {
	Sequence_print_raw(genomicseg);
      } else {
	Sequence_print_two(genomicseg,genomicseg_snp,uppercasep,wraplength);
      }
      Sequence_free(&genomicseg_snp);
      Sequence_free(&genomicseg);
    }
  }


#if 0
  /* Handle alternate strains */
  if (altstrainp == true || user_typestring != NULL) {
    iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			      strlen(fileroot)+strlen(".altstrain.iit")+1,sizeof(char));
    sprintf(iitfile,"%s/%s.altstrain.iit",genomesubdir,fileroot);
    global_altstrain_iit = altstrain_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,
						    /*divread*/READ_ALL,/*divstring*/NULL,/*add_iit_p*/false,
						    /*labels_read_p*/false);
    FREE(iitfile);
  } else {
    global_altstrain_iit = altstrain_iit = (IIT_T) NULL;
  }

  /* We rely upon the fact that gbuffer1 still holds the genomic segment.  This code is duplicated in gmap.c. */
  if (altstrain_iit != NULL) {
    if (user_typestring == NULL) {
      /* No user-specified strain.  Get all indices. */
      indexarray = IIT_get(&nindices,altstrain_iit,/*divstring*/NULL,
			   genomicstart+1,genomicstart+genomiclength-1,/*sortp*/false);
    } else {
      user_type = IIT_typeint(altstrain_iit,user_typestring);
      if (user_type < 0) {
	/* Invalid user-specified strain.  Print nothing. */
	fprintf(stderr,"No such type as %s.  Allowed strains are:\n",user_typestring);
	IIT_dump_typestrings(stderr,altstrain_iit);
	indexarray = NULL;
	nindices = 0;
	Sequence_free(&genomicseg);
      } else {
	/* Valid user-specified strain.  Get subset of indices. */
	indexarray = IIT_get_typed(&nindices,altstrain_iit,/*divstring*/NULL,
				   genomicstart+1,genomicstart+genomiclength-1,user_type,/*sortp*/false);
	if (nindices == 0) {
	  /* Print reference strain */
	  if (rawp == true) {
	    Sequence_print_raw(genomicseg);
	  } else {
	    Sequence_print(stdout,genomicseg,uppercasep,wraplength,/*trimmedp*/false);
	  }
	  Sequence_free(&genomicseg);
	}
      }
    }

    if (nindices > 0) {
      /* Sort according to type and genome position*/
      qsort(indexarray,nindices,sizeof(int),index_compare);

      for (j = 0; j < nindices; ) {
	i = j++;
	type = Interval_type(interval = IIT_interval(altstrain_iit,indexarray[i]));
	strain = IIT_typestring(altstrain_iit,type);
	sourcelength = IIT_annotation_strlen(altstrain_iit,indexarray[i]);
	if (sourcelength > Interval_length(interval)) {
	  extra = sourcelength - Interval_length(interval);
	} else {
	  extra = 0;
	}
	while (j < nindices && Interval_type(interval = IIT_interval(altstrain_iit,indexarray[j])) == type) {
	  sourcelength = IIT_annotation_strlen(altstrain_iit,indexarray[j]);
	  if (sourcelength > Interval_length(interval)) {
	    extra += sourcelength - Interval_length(interval);
	  }
	  j++;
	}
	/* Patch from i to j */
	gbuffer3 = (char *) CALLOC(genomiclength+extra+1,sizeof(char));
	genomicseg = Genome_patch_strain(&(indexarray[i]),j-i,altstrain_iit,
					 genomicstart,genomiclength,
					 revcomp,gbuffer1,gbuffer2,gbuffer3,
					 genomiclength+extra);

	if (header != NULL) {
	  printf(">%s\n",header);
	} else if (revcomp == true) {
	  printf(">%s:%u%s%u_rc variant:%s\n",
		 dbversion,genomicstart+1,SEPARATOR,genomicstart+genomiclength,strain);
	} else {
	  printf(">%s:%u%s%u variant:%s\n",
		 dbversion,genomicstart+1,SEPARATOR,genomicstart+genomiclength,strain);
	}
	if (rawp == true) {
	  Sequence_print_raw(genomicseg);
	} else {
	  Sequence_print(stdout,genomicseg,uppercasep,wraplength,/*trimmedp*/false);
	}
	Sequence_free(&genomicseg);
	FREE(gbuffer3);
      }
      FREE(indexarray);
    }
    IIT_free(&altstrain_iit);
  }
#endif

  return;
}


#if 0
static int *
parse_and_get_matches (int *nmatches, char **divstring, unsigned int *coordstart, unsigned int *coordend,
		       int **leftflanks, int *nleftflanks, int **rightflanks, int *nrightflanks,
		       char *query, char *typestring, IIT_T *iit, char *filename) {
  int *matches;
  bool revcomp;
  int typeint;

  debug(printf("Entering get_matches with query %s.\n",query));
  if (parse_query(&(*divstring),&(*coordstart),&(*coordend),&revcomp,query,filename) == false) {
    /* Treat query as a label */
    *divstring = (char *) NULL;
    if (*iit == NULL) {
      if ((*iit = IIT_read(filename,/*name*/NULL,true,/*divread*/READ_NONE,/*divstring*/NULL,
			   /*add_iit_p*/true,/*labels_read_p*/true)) == NULL) {
	fprintf(stderr,"File %s appears to be an invalid IIT file\n",filename);
	exit(9);
      }
    }
    matches = IIT_find(&(*nmatches),*iit,query);

  } else if (typestring == NULL) {
    /* Treat query as coordinates, without a typestring */
    if (*iit == NULL) {
      if ((*iit = IIT_read(filename,/*name*/NULL,true,/*divread*/READ_ONE,*divstring,
			   /*add_iit_p*/true,/*labels_read_p*/false)) == NULL) {
	fprintf(stderr,"File %s appears to be an invalid IIT file\n",filename);
	exit(9);
      }
    }
    matches = IIT_get(&(*nmatches),*iit,*divstring,*coordstart,*coordend,sortp);
    if (nflanking > 0) {
      IIT_get_flanking(&(*leftflanks),&(*nleftflanks),&(*rightflanks),&(*nrightflanks),*iit,*divstring,
		       *coordstart,*coordend,nflanking,/*sign*/0);
    }

  } else if ((typeint = IIT_typeint(*iit,typestring)) < 0) {
    fprintf(stderr,"No such type as %s.  Ignoring the type.\n",typestring);
    /* Treat query as coordinates, without a typestring */
    if (*iit == NULL) {
      if ((*iit = IIT_read(filename,/*name*/NULL,true,/*divread*/READ_ONE,*divstring,
			   /*add_iit_p*/true,/*labels_read_p*/false)) == NULL) {
	fprintf(stderr,"File %s appears to be an invalid IIT file\n",filename);
	exit(9);
      }
    }
    matches = IIT_get(&(*nmatches),*iit,*divstring,*coordstart,*coordend,sortp);
    if (nflanking > 0) {
      IIT_get_flanking(&(*leftflanks),&(*nleftflanks),&(*rightflanks),&(*nrightflanks),*iit,*divstring,
		       *coordstart,*coordend,nflanking,/*sign*/0);
    }

  } else {
    /* Treat query as coordinates, with a typestring */
    if (*iit == NULL) {
      if ((*iit = IIT_read(filename,/*name*/NULL,true,/*divread*/READ_ONE,*divstring,
			   /*add_iit_p*/true,/*labels_read_p*/false)) == NULL) {
	fprintf(stderr,"File %s appears to be an invalid IIT file\n",filename);
	exit(9);
      }
    }
    matches = IIT_get_typed(&(*nmatches),*iit,*divstring,*coordstart,*coordend,typeint,sortp);
    if (nflanking > 0) {
      IIT_get_flanking_typed(&(*leftflanks),&(*nleftflanks),&(*rightflanks),&(*nrightflanks),*iit,*divstring,
			     *coordstart,*coordend,nflanking,typeint);
    }
  }

  return matches;
}
#endif


static int *
get_matches (int *nmatches, int *sign, char *divstring, unsigned int coordstart, unsigned int coordend, bool revcomp,
	     int **leftflanks, int *nleftflanks, int **rightflanks, int *nrightflanks,
	     char *typestring, IIT_T *iit, char *filename) {
  int *matches;
  int typeint;

  if (signedp == false) {
    *sign = 0;
  } else if (revcomp == true) {
    *sign = -1;
  } else if (coordend == coordstart) {
    *sign = 0;
  } else {
    *sign = +1;
  }

  if (typestring == NULL) {
    /* Treat query as coordinates, without a typestring */
    if (*iit == NULL) {
      /* This call can give a warning if there are no entries in the map IIT file for the given div */
      if ((*iit = IIT_read(filename,/*name*/NULL,true,/*divread*/READ_ONE,divstring,
			   /*add_iit_p*/true,/*labels_read_p*/false)) == NULL) {
	fprintf(stderr,"File %s appears to be an invalid IIT file\n",filename);
	exit(9);
      }
    }

    matches = IIT_get_signed(&(*nmatches),*iit,divstring,coordstart,coordend,*sign,sortp);
    if (nflanking > 0) {
      IIT_get_flanking(&(*leftflanks),&(*nleftflanks),&(*rightflanks),&(*nrightflanks),*iit,divstring,
		       coordstart,coordend,nflanking,*sign);
    }

  } else if ((typeint = IIT_typeint(*iit,typestring)) < 0) {
    fprintf(stderr,"No such type as %s.  Ignoring the type.\n",typestring);
    /* Treat query as coordinates, without a typestring */
    if (*iit == NULL) {
      if ((*iit = IIT_read(filename,/*name*/NULL,true,/*divread*/READ_ONE,divstring,
			   /*add_iit_p*/true,/*labels_read_p*/false)) == NULL) {
	fprintf(stderr,"File %s appears to be an invalid IIT file\n",filename);
	exit(9);
      }
    }
    matches = IIT_get_signed(&(*nmatches),*iit,divstring,coordstart,coordend,*sign,sortp);
    if (nflanking > 0) {
      IIT_get_flanking(&(*leftflanks),&(*nleftflanks),&(*rightflanks),&(*nrightflanks),*iit,divstring,
		       coordstart,coordend,nflanking,*sign);
    }

  } else {
    /* Treat query as coordinates, with a typestring */
    if (*iit == NULL) {
      if ((*iit = IIT_read(filename,/*name*/NULL,true,/*divread*/READ_ONE,divstring,
			   /*add_iit_p*/true,/*labels_read_p*/false)) == NULL) {
	fprintf(stderr,"File %s appears to be an invalid IIT file\n",filename);
	exit(9);
      }
    }
    matches = IIT_get_typed_signed(&(*nmatches),*iit,divstring,coordstart,coordend,typeint,*sign,sortp);
    if (nflanking > 0) {
      IIT_get_flanking_typed(&(*leftflanks),&(*nleftflanks),&(*rightflanks),&(*nrightflanks),*iit,divstring,
			     coordstart,coordend,nflanking,typeint,*sign);
    }
  }

  return matches;
}


static int *
get_matches_multiple_typed (int *nmatches, char **divstring, unsigned int *coordstart, unsigned int *coordend,
			    int **leftflanks, int *nleftflanks, int **rightflanks, int *nrightflanks,
			    char *query, int *types, int ntypes, IIT_T *iit, char *filename) {
  int *matches;
  bool revcomp;

  if (parse_query(&(*divstring),&(*coordstart),&(*coordend),&revcomp,query,filename) == false) {
    /* Not expecting a label */
    abort();
  }

  if ((*iit = IIT_read(filename,/*name*/NULL,true,/*divread*/READ_ONE,*divstring,/*add_iit_p*/true,
		       /*labels_read_p*/false)) == NULL) {
    fprintf(stderr,"File %s appears to be an invalid IIT file\n",filename);
    exit(9);
  }

  matches = IIT_get_multiple_typed(&(*nmatches),*iit,*divstring,*coordstart,*coordend,types,ntypes,sortp);
  if (nflanking > 0) {
    IIT_get_flanking_multiple_typed(&(*leftflanks),&(*nleftflanks),&(*rightflanks),&(*nrightflanks),*iit,*divstring,
				    *coordstart,*coordend,nflanking,types,ntypes);
  }

  return matches;
}

static char complCode[128] = COMPLEMENT_LC;

static void
make_complement_buffered (char *complement, char *sequence, Genomicpos_T length) {
  int i, j;

  for (i = length-1, j = 0; i >= 0; i--, j++) {
    complement[j] = complCode[(int) sequence[i]];
  }
  complement[length] = '\0';
  return;
}

static void
make_complement_inplace (char *sequence, Genomicpos_T length) {
  char temp;
  Genomicpos_T i, j;

  for (i = 0, j = length-1; i < length/2; i++, j--) {
    temp = complCode[(int) sequence[i]];
    sequence[i] = complCode[(int) sequence[j]];
    sequence[j] = temp;
  }
  if (i == j) {
    sequence[i] = complCode[(int) sequence[i]];
  }

  return;
}



static void
genemap_print_exons (char *annot, Genomicpos_T chroffset, Genome_T genome) {
  char *p;
  char *gbuffer;
  Genomicpos_T exonstart, exonend, exonlow, exonhigh, exonlength;

  /* Skip header */
  p = annot;
  while (*p != '\0' && *p != '\n') {
    putchar(*p);
    p++;
  }
  if (*p == '\n') p++;
  printf("\n");

  while (*p != '\0') {
    if (sscanf(p,"%u %u",&exonstart,&exonend) != 2) {
      fprintf(stderr,"Can't parse exon coordinates in %s\n",p);
      abort();
    } else {
      if (exonstart <= exonend) {
	exonlow = exonstart;
	exonhigh = exonend;
	exonlength = exonhigh - exonlow + 1;
	gbuffer = (char *) CALLOC(exonlength+1,sizeof(char));
	Genome_fill_buffer_simple(genome,/*left*/chroffset-1U+exonlow,exonlength,gbuffer);
      } else {
	exonlow = exonend;
	exonhigh = exonstart;
	exonlength = exonhigh - exonlow + 1;
	gbuffer = (char *) CALLOC(exonlength+1,sizeof(char));
	Genome_fill_buffer_simple(genome,/*left*/chroffset-1U+exonlow,exonlength,gbuffer);
	make_complement_inplace(gbuffer,exonlength);
      }
      printf("%s\n",gbuffer);
      FREE(gbuffer);
    }

    while (*p != '\0' && *p != '\n') p++;
    if (*p == '\n') p++;
  }

  return;
}


static void
genemap_print_sequence (char *annot, Genomicpos_T chroffset, Genome_T genome) {
  char *p;
  char *gbuffer;
  Genomicpos_T exonstart, exonend, exonlow, exonhigh, exonlength;
  int col, i;

  /* Skip header */
  p = annot;
  while (*p != '\0' && *p != '\n') {
    putchar(*p);
    p++;
  }
  if (*p == '\n') p++;
  printf("\n");

  col = 0;
  while (*p != '\0') {
    if (sscanf(p,"%u %u",&exonstart,&exonend) != 2) {
      fprintf(stderr,"Can't parse exon coordinates in %s\n",p);
      abort();
    } else {
      if (exonstart <= exonend) {
	exonlow = exonstart;
	exonhigh = exonend;
	exonlength = exonhigh - exonlow + 1;
	gbuffer = (char *) CALLOC(exonlength+1,sizeof(char));
	Genome_fill_buffer_simple(genome,/*left*/chroffset-1U+exonlow,exonlength,gbuffer);
      } else {
	exonlow = exonend;
	exonhigh = exonstart;
	exonlength = exonhigh - exonlow + 1;
	gbuffer = (char *) CALLOC(exonlength+1,sizeof(char));
	Genome_fill_buffer_simple(genome,/*left*/chroffset-1U+exonlow,exonlength,gbuffer);
	make_complement_inplace(gbuffer,exonlength);
      }
      
      for (i = 0; i < (int) exonlength; i++) {
	putchar(gbuffer[i]);
	if (++col >= wraplength) {
	  printf("\n");
	  col = 0;
	}
      }
      FREE(gbuffer);
    }

    while (*p != '\0' && *p != '\n') p++;
    if (*p == '\n') p++;
  }

  if (col > 0) {
    printf("\n");
  }

  return;
}


static void
print_interval (char *divstring, int index, IIT_T iit, int ndivs, IIT_T chromosome_iit,
		Genome_T genome, int fieldint) {
  Interval_T interval;
  char *label, *annotation, *restofheader;
  bool allocp;
  bool annotationonlyp = false, signed_output_p = true;
  int divno;
  Genomicpos_T chroffset;

  if (exonsp == true || sequencep == true) {
    label = IIT_label(iit,index,&allocp);
    printf(">%s ",label);
    if (allocp == true) {
      FREE(label);
    }

    annotation = IIT_annotation(&restofheader,iit,index,&allocp);
    printf("%s",restofheader);
    if (IIT_version(iit) < 5) {
      printf("\n");
    }

    if (ndivs > 1) {
      if (divstring == NULL) {
	/* For example, if interval was retrieved by label */
	divstring = IIT_divstring_from_index(iit,index);
      }
    }

    if ((divno = IIT_find_one(chromosome_iit,divstring)) < 0) {
      fprintf(stderr,"Cannot find chromosome %s in chromosome IIT file\n",divstring);
      /* exit(9); */
    } else {
      interval = IIT_interval(chromosome_iit,divno);
      chroffset = Interval_low(interval);
      if (exonsp == true) {
	genemap_print_exons(annotation,chroffset,genome);
      } else if (sequencep == true) {
	genemap_print_sequence(annotation,chroffset,genome);
      }
    }
    if (allocp == true) {
      FREE(restofheader);
    }

  } else {
    if (annotationonlyp == false) {
      label = IIT_label(iit,index,&allocp);
      printf(">%s ",label);
      if (allocp == true) {
	FREE(label);
      }
      
      if (ndivs > 1) {
	if (divstring == NULL) {
	  /* For example, if interval was retrieved by label */
	  divstring = IIT_divstring_from_index(iit,index);
	}
	printf("%s:",divstring);
      }

      interval = IIT_interval(iit,index);
      if (signed_output_p == false) {
	printf("%u..%u",Interval_low(interval),Interval_high(interval));
      } else if (Interval_sign(interval) < 0) {
	printf("%u..%u",Interval_high(interval),Interval_low(interval));
      } else {
	printf("%u..%u",Interval_low(interval),Interval_high(interval));
      }
      if (Interval_type(interval) > 0) {
	printf(" %s",IIT_typestring(iit,Interval_type(interval)));
      }
#if 0
      /* Unnecessary because of "\n" after restofheader below */
      if (IIT_version(iit) < 5) {
	printf("\n");
      }
#endif

    }

    if (fieldint < 0) {
      annotation = IIT_annotation(&restofheader,iit,index,&allocp);
      printf("%s\n",restofheader);
      printf("%s",annotation);
      if (allocp == true) {
	FREE(restofheader);
      }

    } else {
      annotation = IIT_fieldvalue(iit,index,fieldint);
      printf("%s",annotation);
      FREE(annotation);
    }
  }

  return;
}


#define BUFFERLEN 1024


int
main (int argc, char *argv[]) {
  char *snpsdir = NULL;
  char *iitfile, *chrsubsetfile;
  FILE *fp;
  Genome_T genome = NULL, genomealt = NULL;
  Genomicpos_T genomicstart, genomiclength;
  Genomicpos_T chroffset, chrlength, chrstart, chrend;
  char *mapdir = NULL, *fileroot = NULL, *p;
  char *divstring, *divstring2;
  IIT_T chromosome_iit, contig_iit, map_iit = NULL;
  Chrsubset_T chrsubset;
  char Buffer[BUFFERLEN], subsetname[BUFFERLEN], *segment;

  int fieldint = -1;
  int *matches, nmatches, ndivs, i, *leftflanks, *rightflanks, nleftflanks = 0, nrightflanks = 0;
  int sign;
  unsigned int coordstart, coordend;

  char *chr, *with_colon;
  int indx;
  bool allocp;

  int opt;
  extern int optind;
  extern char *optarg;
  int long_option_index = 0;
  const char *long_name;

  while ((opt = getopt_long(argc,argv,"D:d:CUl:Gh:V:v:f:M:m:kru:ESALIc^?",
			    long_options,&long_option_index)) != -1) {
    switch (opt) {
    case 0:
      long_name = long_options[long_option_index].name;
      if (!strcmp(long_name,"version")) {
	print_program_version();
	exit(0);
      } else if (!strcmp(long_name,"help")) {
	print_program_usage();
	exit(0);

      } else if (!strcmp(long_name,"forsam")) {
	dumpchrp = true;
	dumpchr_forsam_p = true;

      } else {
	/* Shouldn't reach here */
	fprintf(stderr,"Don't recognize option %s.  For usage, run 'get-genome --help'",long_name);
	exit(9);
      }
      break;

    case 'D': user_genomedir = optarg; break;
    case 'd': 
      dbroot = (char *) CALLOC(strlen(optarg)+1,sizeof(char));
      strcpy(dbroot,optarg);
      break;

    case 'C': coordp = true; break;
    case 'U': uppercasep = true; break;
    case 'l': wraplength = atoi(optarg); break;
    case 'G': uncompressedp = true; break;
    case 'h': header = optarg; break;
    case 'V': user_snpsdir = optarg; break;
    case 'v': snps_root = optarg; break;
    case 'f': print_snps_mode = atoi(optarg); break;

    case 'M': user_mapdir = optarg; break;
    case 'm': map_iitfile = optarg; break;
    case 'k': levelsp = true; break;
    case 'r': uncompressedp = true; rawp = true; break;
    case 'u': nflanking = atoi(optarg); break;
    case 'E': exonsp = true; break;
    case 'S': sequencep = true; break;
    case 's': signedp = true; break;

    case 'A': dumpallp = true; break;
    case 'L': dumpchrp = true; dumpchr_forsam_p = false; break;
    case 'I': dumpsegsp = true; break;
    case 'c': dumpchrsubsetsp = true; break;

    case '^': print_program_version(); exit(0);
    case '?': print_program_usage(); exit(0);
    default: exit(9);
    }
  }
  argc -= optind;
  argv += optind;
  
  if (dbroot == NULL) {
    print_program_usage();
    exit(9);
  } else if (!strcmp(dbroot,"?")) {
    Datadir_avail_gmap_databases(stdout,user_genomedir);
    exit(0);
  } else {
    genomesubdir = Datadir_find_genomesubdir(&fileroot,&dbversion,user_genomedir,dbroot);
  }

  if (user_snpsdir == NULL) {
    snpsdir = genomesubdir;
  } else {
    snpsdir = user_snpsdir;
  }

  if (dumpallp == true) {
    iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			      strlen(fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
    sprintf(iitfile,"%s/%s.chromosome.iit",genomesubdir,fileroot);
    chromosome_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
			      /*divstring*/NULL,/*add_iit_p*/false,/*labels_read_p*/true);
    FREE(iitfile);

    if (snps_root == NULL || print_snps_mode == 0) {
      genome = Genome_new(genomesubdir,fileroot,/*snps_root*/NULL,uncompressedp,
			  /*access*/USE_MMAP_ONLY);
    } else if (print_snps_mode == 2) {
      genome = Genome_new(snpsdir,fileroot,snps_root,uncompressedp,
			  /*access*/USE_MMAP_ONLY);
    } else if (print_snps_mode == 1 || print_snps_mode == 3) {
      genome = Genome_new(genomesubdir,fileroot,/*snps_root*/NULL,uncompressedp,
			  /*access*/USE_MMAP_ONLY);
      genomealt = Genome_new(snpsdir,fileroot,snps_root,uncompressedp,
			     /*access*/USE_MMAP_ONLY);
    }

    for (indx = 1; indx <= IIT_total_nintervals(chromosome_iit); indx++) {
      chr = IIT_label(chromosome_iit,indx,&allocp);
      with_colon = (char *) CALLOC(strlen(chr)+strlen(":")+1,sizeof(char));
      sprintf(with_colon,"%s:",chr);
      if (allocp == true) {
	FREE(chr);
      }
      if (Parserange_universal(&segment,&revcomp,&genomicstart,&genomiclength,&chrstart,&chrend,
			       &chroffset,&chrlength,with_colon,genomesubdir,fileroot) == true) {
	print_sequence(genome,genomealt,genomicstart,genomiclength,chromosome_iit,
		       /*whole_chromosome_p*/true);
      }
      FREE(with_colon);
    }

    if (genomealt != NULL) {
      Genome_free(&genomealt);
    }
    Genome_free(&genome);

    IIT_free(&chromosome_iit);

    return 0;

  } else if (dumpchrp == true) {
    iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			      strlen(fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
    sprintf(iitfile,"%s/%s.chromosome.iit",genomesubdir,fileroot);
    chromosome_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
			      /*divstring*/NULL,/*add_iit_p*/false,/*labels_read_p*/true);
    FREE(iitfile);

    if (dumpchr_forsam_p == true) {
      IIT_dump_sam(/*fp*/stdout,chromosome_iit,/*sam_read_group_id*/NULL,/*sam_read_group_name*/NULL,
		   /*sam_read_group_library*/NULL,/*sam_read_group_platform*/NULL);
    } else {
      IIT_dump_simple(chromosome_iit);
    }
    IIT_free(&chromosome_iit);
    return 0;

  } else if (dumpsegsp == true) {
    iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			      strlen(fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
    sprintf(iitfile,"%s/%s.chromosome.iit",genomesubdir,fileroot);
    chromosome_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
			      /*divstring*/NULL,/*add_iit_p*/false,/*labels_read_p*/true);
    FREE(iitfile);

    iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			      strlen(fileroot)+strlen(".contig.iit")+1,sizeof(char));
    sprintf(iitfile,"%s/%s.contig.iit",genomesubdir,fileroot);
    contig_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
			  /*divstring*/NULL,/*add_iit_p*/false,/*labels_read_p*/true);
    FREE(iitfile);

    IIT_dump_version1(contig_iit,chromosome_iit,/*directionalp*/true);
    IIT_free(&contig_iit);
    return 0;

  } else if (dumpchrsubsetsp == true) {
    iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			      strlen(fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
    sprintf(iitfile,"%s/%s.chromosome.iit",genomesubdir,fileroot);
    chromosome_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
			      /*divstring*/NULL,/*add_iit_p*/false,/*labels_read_p*/true);
    FREE(iitfile);

    chrsubsetfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
				    strlen(fileroot)+strlen(".chrsubset")+1,sizeof(char));
    sprintf(chrsubsetfile,"%s/%s.chrsubset",genomesubdir,fileroot);
    if ((fp = fopen(chrsubsetfile,"r")) == NULL) {
      fprintf(stderr,"Cannot open file %s\n",chrsubsetfile);
      exit(9);
    } else {
      while (fgets(Buffer,BUFFERLEN,fp) != NULL) {
	if (Buffer[0] == '>') {
	  if (Buffer[1] == '\0' || isspace(Buffer[1])) {
	    /* Skip */
	  } else {
	    if ((p = rindex(Buffer,'\n')) != NULL) {
	      *p = '\0';
	    }
	    sscanf(&(Buffer[1]),"%s",subsetname);
	    printf("%s\t",subsetname);
	    chrsubset = Chrsubset_read(chrsubsetfile,genomesubdir,fileroot,subsetname,chromosome_iit);
	    Chrsubset_print_chromosomes(chrsubset,chromosome_iit);
	    Chrsubset_free(&chrsubset);
	    printf("\n");
	  }
	}
      }
      fclose(fp);
    }

    FREE(chrsubsetfile);
    IIT_free(&chromosome_iit);

    return 0;

  }

#if 0
  if (replacex == true) {
    Genome_replace_x();
  }
#endif

  iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			    strlen(fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
  sprintf(iitfile,"%s/%s.chromosome.iit",genomesubdir,fileroot);
  chromosome_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
			    /*divstring*/NULL,/*add_iit_p*/false,/*labels_read_p*/true);
  FREE(iitfile);


  if (argc >= 1) {
    if (coordp == true) {
      debug(printf("coordp is true\n"));
      if (Parserange_universal(&segment,&revcomp,&genomicstart,&genomiclength,&chrstart,&chrend,
			       &chroffset,&chrlength,argv[0],genomesubdir,fileroot) == true) {
	debug(printf("Query %s parsed as: genomicstart = %u, genomiclength = %u, revcomp = %d\n",
		     argv[0],genomicstart,genomiclength,revcomp));
	print_two_coords(genomicstart,genomiclength,chromosome_iit);
      }
      
    } else if (map_iitfile == NULL) {
      debug(printf("No map file\n"));

      if (snps_root == NULL || print_snps_mode == 0) {
	genome = Genome_new(genomesubdir,fileroot,/*snps_root*/NULL,uncompressedp,
			    /*access*/USE_MMAP_ONLY);
      } else if (print_snps_mode == 2) {
	genome = Genome_new(snpsdir,fileroot,snps_root,uncompressedp,
			    /*access*/USE_MMAP_ONLY);
      } else if (print_snps_mode == 1 || print_snps_mode == 3) {
	genome = Genome_new(genomesubdir,fileroot,/*snps_root*/NULL,uncompressedp,
			    /*access*/USE_MMAP_ONLY);
	genomealt = Genome_new(snpsdir,fileroot,snps_root,uncompressedp,
			       /*access*/USE_MMAP_ONLY);
      }

      if (Parserange_universal(&segment,&revcomp,&genomicstart,&genomiclength,&chrstart,&chrend,
			       &chroffset,&chrlength,argv[0],genomesubdir,fileroot) == true) {
	debug(printf("Query %s parsed as: genomicstart = %u, genomiclength = %u, revcomp = %d\n",
		     argv[0],genomicstart,genomiclength,revcomp));
	print_sequence(genome,genomealt,genomicstart,genomiclength,chromosome_iit,
		       /*whole_chromosome_p*/false);
      }

    } else if (!strcmp(map_iitfile,"?")) {
      Datadir_avail_maps(stdout,user_mapdir,genomesubdir,fileroot);
      exit(0);

    } else {
      debug(printf("Map file\n"));

      mapdir = Datadir_find_mapdir(user_mapdir,genomesubdir,fileroot);
      iitfile = (char *) CALLOC(strlen(mapdir)+strlen("/")+
				strlen(map_iitfile)+strlen(".iit")+1,sizeof(char));
      sprintf(iitfile,"%s/%s.iit",mapdir,map_iitfile);
      if (Access_file_exists_p(iitfile) == false) {
	fprintf(stderr,"Map file %s.iit not found in %s.  Available files:\n",map_iitfile,mapdir);
	Datadir_list_directory(stderr,mapdir);
	fprintf(stderr,"Either install file %s.iit or specify a full directory path\n",map_iitfile);
	fprintf(stderr,"using the -M flag to gmap.\n");
	exit(9);
      }

      if (exonsp == true || sequencep == true) {
	if (snps_root == NULL || print_snps_mode == 0) {
	  genome = Genome_new(genomesubdir,fileroot,/*snps_root*/NULL,uncompressedp,
			      /*access*/USE_MMAP_ONLY);
	} else if (print_snps_mode == 2) {
	  genome = Genome_new(snpsdir,fileroot,snps_root,uncompressedp,
			      /*access*/USE_MMAP_ONLY);
	} else if (print_snps_mode == 1 || print_snps_mode == 3) {
	  genome = Genome_new(genomesubdir,fileroot,/*snps_root*/NULL,uncompressedp,
			      /*access*/USE_MMAP_ONLY);
#if 0
	  genomealt = Genome_new(snpsdir,fileroot,snps_root,uncompressedp,
				 /*access*/USE_MMAP_ONLY);
#endif
	}
      }

      if (Parserange_universal(&segment,&revcomp,&genomicstart,&genomiclength,&chrstart,&chrend,
			       &chroffset,&chrlength,argv[0],genomesubdir,fileroot) == true) {
	debug(printf("Query %s parsed as: genomicstart = %u, genomiclength = %u, revcomp = %d\n",
		     argv[0],genomicstart,genomiclength,revcomp));
	divstring = IIT_string_from_position(&coordstart,genomicstart,chromosome_iit);
	divstring2 = IIT_string_from_position(&coordend,genomicstart+genomiclength-1U,chromosome_iit);
	if (strcmp(divstring,divstring2)) {
	  fprintf(stderr,"Coordinates cross chromosomal boundary\n");
	  exit(9);
	} else {
	  coordstart += 1U;
	  coordend += 1U;
	  debug(printf("Query translated to %s:%u..%u\n",divstring,coordstart,coordend));
	}
	matches = get_matches(&nmatches,&sign,divstring,coordstart,coordend,revcomp,
			      &leftflanks,&nleftflanks,&rightflanks,&nrightflanks,
			      /*typestring*/NULL,&map_iit,/*filename*/iitfile);
	ndivs = IIT_ndivs(map_iit);
	if (nflanking > 0) {
	  if (sign != +1) {
	    for (i = nleftflanks-1; i >= 0; i--) {
	      print_interval(divstring,leftflanks[i],map_iit,ndivs,chromosome_iit,genome,fieldint);
	    }
	  }
	  printf("====================\n");
	  FREE(leftflanks);
	}

	for (i = 0; i < nmatches; i++) {
	  print_interval(divstring,matches[i],map_iit,ndivs,chromosome_iit,genome,fieldint);
	}

	if (nflanking > 0) {
	  printf("====================\n");
	  if (sign != -1) {
	    for (i = 0; i < nrightflanks; i++) {
	      print_interval(divstring,rightflanks[i],map_iit,ndivs,chromosome_iit,genome,fieldint);
	    }
	  }
	  FREE(rightflanks);
	}

      } else {
	/* Must have been a label */
#if 0
	if ((*iit = IIT_read(filename,/*name*/NULL,true,/*divread*/READ_NONE,/*divstring*/NULL,
			       /*add_iit_p*/true,/*labels_read_p*/true)) == NULL) {
	}
#endif
	if ((map_iit = IIT_read(iitfile,/*name*/NULL,true,/*divread*/READ_NONE,/*divstring*/NULL,
				/*add_iit_p*/true,/*labels_read_p*/true)) == NULL) {
	  fprintf(stderr,"Cannot open IIT file %s\n",iitfile);
	  exit(9);
	}
	matches = IIT_find(&nmatches,map_iit,argv[0]);
	debug(printf("Looking up label %s => %d matches\n",argv[0],nmatches));
	ndivs = IIT_ndivs(map_iit);

#if 0
	/* Not doing flanking */
	if (nflanking > 0) {
	  for (i = nleftflanks-1; i >= 0; i--) {
	    print_interval(/*divstring*/NULL,leftflanks[i],map_iit,ndivs,chromosome_iit,genome,fieldint);
	  }
	  printf("====================\n");
	  FREE(leftflanks);
	}
#endif

	for (i = 0; i < nmatches; i++) {
	  print_interval(/*divstring*/NULL,matches[i],map_iit,ndivs,chromosome_iit,genome,fieldint);
	}
	FREE(matches);

#if 0
	/* Not doing flanking */
	if (nflanking > 0) {
	  printf("====================\n");
	  for (i = 0; i < nrightflanks; i++) {
	    print_interval(/*divstring*/NULL,rightflanks[i],map_iit,ndivs,chromosome_iit,genome,fieldint);
	  }
	  FREE(rightflanks);
	}
#endif

      }
      FREE(iitfile);
      FREE(mapdir);
    }

  } else {
    /* Read from stdin */

    if (map_iitfile == NULL) {
      /* Skip */

    } else if (!strcmp(map_iitfile,"?")) {
      Datadir_avail_maps(stdout,user_mapdir,genomesubdir,fileroot);
      exit(0);

    } else {
      mapdir = Datadir_find_mapdir(user_mapdir,genomesubdir,fileroot);
      iitfile = (char *) CALLOC(strlen(mapdir)+strlen("/")+
				strlen(map_iitfile)+strlen(".iit")+1,sizeof(char));
      sprintf(iitfile,"%s/%s.iit",mapdir,map_iitfile);

      if ((map_iit = IIT_read(iitfile,/*name*/map_iitfile,/*readonlyp*/true,
			      /*divread*/READ_ALL,/*divstring*/NULL,/*add_iit_p*/false,
			      /*labels_read_p*/false)) == NULL) {
	fprintf(stderr,"Map file %s.iit not found in %s.  Available files:\n",map_iitfile,mapdir);
	Datadir_list_directory(stderr,mapdir);
	fprintf(stderr,"Either install file %s.iit or specify a full directory path\n",map_iitfile);
	fprintf(stderr,"using the -M flag to gmap.\n");
	exit(9);
#if 0
      } else if (map_typestring != NULL) {
	map_type = IIT_typeint(map_iit,map_typestring);
	if (map_type < 0) {
	  /* Invalid user-specified strain.  Print nothing. */
	  fprintf(stderr,"No such type as %s.  Allowed strains are:\n",map_typestring);
	  IIT_dump_typestrings(stderr,map_iit);
	  exit(9);
	}
#endif
      } else {
	ndivs = IIT_ndivs(map_iit);
      }
    }
    FREE(iitfile);

    if (snps_root == NULL || print_snps_mode == 0) {
      genome = Genome_new(genomesubdir,fileroot,/*snps_root*/NULL,uncompressedp,
			  /*access*/USE_MMAP_ONLY);
    } else if (print_snps_mode == 2) {
      genome = Genome_new(snpsdir,fileroot,snps_root,uncompressedp,
			  /*access*/USE_MMAP_ONLY);
    } else if (print_snps_mode == 1 || print_snps_mode == 3) {
      genome = Genome_new(genomesubdir,fileroot,/*snps_root*/NULL,uncompressedp,
			  /*access*/USE_MMAP_ONLY);
      genomealt = Genome_new(snpsdir,fileroot,snps_root,uncompressedp,
			     /*access*/USE_MMAP_ONLY);
    }

    iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+strlen(fileroot)+
			      strlen(".contig.iit")+1,sizeof(char));
    sprintf(iitfile,"%s/%s.contig.iit",genomesubdir,fileroot);
    contig_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
			  /*divstring*/NULL,/*add_iit_p*/false,/*labels_read_p*/true);
    FREE(iitfile);

    while (fgets(Buffer,BUFFERLEN,stdin) != NULL) {
      if ((p = rindex(Buffer,'\n')) != NULL) {
	*p = '\0';
      }
      fprintf(stdout,"# Query: %s\n",Buffer);
      if (Parserange_universal_iit(&segment,&revcomp,&genomicstart,&genomiclength,&chrstart,&chrend,
				   &chroffset,&chrlength,Buffer,chromosome_iit,contig_iit) == true) {
	debug(printf("Query %s parsed as: genomicstart = %u, genomiclength = %u, revcomp = %d\n",
		     Buffer,genomicstart,genomiclength,revcomp));
	divstring = IIT_string_from_position(&coordstart,genomicstart,chromosome_iit);
	divstring2 = IIT_string_from_position(&coordend,genomicstart+genomiclength-1U,chromosome_iit);
	if (strcmp(divstring,divstring2)) {
	  fprintf(stderr,"Coordinates cross chromosomal boundary\n");
	  exit(9);
	} else {
	  coordstart += 1U;
	  coordend += 1U;
	  debug(printf("Query translated to %s:%u..%u\n",divstring,coordstart,coordend));
	}

	if (map_iit != NULL) {
	  matches = get_matches(&nmatches,&sign,divstring,coordstart,coordend,revcomp,
				&leftflanks,&nleftflanks,&rightflanks,&nrightflanks,
				/*typestring*/NULL,&map_iit,/*filename*/NULL);
	  if (nflanking > 0) {
	    if (sign != +1) {
	      for (i = nleftflanks-1; i >= 0; i--) {
		print_interval(divstring,leftflanks[i],map_iit,ndivs,chromosome_iit,genome,fieldint);
	      }
	    }
	    printf("====================\n");
	    FREE(leftflanks);
	  }

	  for (i = 0; i < nmatches; i++) {
	    print_interval(divstring,matches[i],map_iit,ndivs,chromosome_iit,genome,fieldint);
	  }

	  if (nflanking > 0) {
	    printf("====================\n");
	    if (sign != -1) {
	      for (i = 0; i < nrightflanks; i++) {
		print_interval(divstring,rightflanks[i],map_iit,ndivs,chromosome_iit,genome,fieldint);
	      }
	    }
	    FREE(rightflanks);
	  }

	  fflush(stdout);
	  FREE(matches);
	  
	} else if (coordp == true) {
	  print_two_coords(genomicstart,genomiclength,chromosome_iit);
	  
	} else {
	  print_sequence(genome,genomealt,genomicstart,genomiclength,chromosome_iit,
			 /*whole_chromosome_p*/false);
	}

      } else {
	/* Must have been a label */

	if (map_iit != NULL) {
	  matches = IIT_find(&nmatches,map_iit,Buffer);
	  debug(printf("Looking up label %s => %d matches\n",Buffer,nmatches));

	  /* Not doing flanking */

	  for (i = 0; i < nmatches; i++) {
	    print_interval(/*divstring*/NULL,matches[i],map_iit,ndivs,chromosome_iit,genome,fieldint);
	  }

	  /* Not doing flanking */

	  fflush(stdout);
	  FREE(matches);
	}
      }

      fprintf(stdout,"# End\n");
    }

    IIT_free(&contig_iit);

    if (map_iitfile != NULL) {
      FREE(iitfile);
      FREE(mapdir);
    }

  }
    
  if (genomealt != NULL) {
    Genome_free(&genomealt);
  }
  Genome_free(&genome);

  IIT_free(&chromosome_iit);

  if (map_iitfile != NULL) {
    IIT_free(&map_iit);
  }

  FREE(dbversion);
  FREE(genomesubdir);
  FREE(fileroot);
  FREE(dbroot);

  return 0;
}

