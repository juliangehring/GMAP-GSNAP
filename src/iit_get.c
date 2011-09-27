static char rcsid[] = "$Id: iit_get.c 40271 2011-05-28 02:29:18Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>		/* For getopt */
#endif
#include <string.h>
#include <strings.h>		/* For index, rindex */
#include <ctype.h>
#include <math.h>		/* For log */
#include "bool.h"
#include "mem.h"
#include "access.h"
#include "intlist.h"
#include "interval.h"
#include "iit-read.h"
#include "complement.h"
#include "fopen.h"
#include "parserange.h"
#include "getopt.h"

#define BUFLEN 1024


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/************************************************************************
 *   Program options
 ************************************************************************/

static char *fieldstring = NULL;
static bool annotationonlyp = false;
static bool sortp = false;
static bool signedp = true;
static int nflanking = 0;
static bool centerp = false;
static int centerlength = 0;
static bool centeruc = false;
static bool runlengthp = false;
static bool tallyp = false;
static bool zeroesp = false;
static bool statsp = false;
static bool labelp = false;
static bool overall_total_p = false;

static struct option long_options[] = {
  /* Input options */
  {"field", required_argument, 0, 'f'},	/* fieldstring */
  {"annotonly", no_argument, 0, 'A'},	/* annotationonlyp */
  {"sort", no_argument, 0, 'S'},	/* sortp */
  {"unsigned", no_argument, 0, 'U'},	/* signedp */
  {"flanking", required_argument, 0, 'u'},	/* nflanking */
  {"center", required_argument, 0, 'c'}, /* centerp, centerlength */
  {"centeruc", no_argument, 0, 'H'}, /* centeruc */
  {"runlength", no_argument, 0, 'R'}, /* runlengthp */
  {"tally", no_argument, 0, 'T'}, /* tallyp */
  {"zeroes", no_argument, 0, 'Z'}, /* zeroesp */
  {"stats", no_argument, 0, 'N'}, /* statsp */
  {"label", no_argument, 0, 'L'}, /* labelp */

  /* Help options */
  {"version", no_argument, 0, 'V'}, /* print_program_version */
  {"help", no_argument, 0, '?'}, /* print_program_usage */
  {0, 0, 0, 0}
};

static void
print_program_version () {
  fprintf(stdout,"\n");
  fprintf(stdout,"iit_get: retrieval utility for Interval Index Trees\n");
  fprintf(stdout,"Part of GMAP package, version %s\n",PACKAGE_VERSION);
  fprintf(stdout,"Thomas D. Wu, Genentech, Inc.\n");
  fprintf(stdout,"Contact: twu@gene.com\n");
  fprintf(stdout,"\n");
  return;
}

static void
print_program_usage () {
  fprintf(stdout,"\
Usage: iit_get [OPTIONS...] iitfile query\n\
       iit_get [OPTIONS...] iitfile query types...\n\
       iit_get [OPTIONS...] iitfile\n\
\n\
where query is one of the following forms:\n\
\n\
   chr:start..end\n\
   chr:start\n\
   chr:\n\
   start..end\n\
   start\n\
   label\n\
\n\
Options\n\
  -f, --field=STRING      Show given field part of the annotation\n\
  -L, --label             Interpret query as a label, even if it is numeric\n\
  -A, --annotonly         Show annotation lines only (no headers)\n\
  -S, --sort              Sort results by coordinates\n\
  -U, --unsigned          Print all intervals as low..high, even those entered as reverse (high < low)\n\
  -u, --flanking=INT      Show flanking segments on left and right\n\
\n\
Options for specific IIT formats\n\
  -c, --center=INT        Align reads so given position is centered at given column\n\
  -H, --centeruc          Report only reads with upper-case letter at given position\n\
  -R, --runlength         Report runlength IIT file in tally format\n\
  -T, --tally             Report tally IIT file in tally format\n\
  -Z, --zeroes            Include zeroes in tally format\n\
  -N, --stats             Statistics (count, npositions, mean) of tally format\n\
\n\
  -V, --version           Show version\n\
  -?, --help              Show this help message\n\
\n\
\n\
The iit_get program retrieves segments from an iit file that overlap a\n\
given coordinate or pair of coordinates.  Retrieval is done in\n\
logarithmic time (for more details on IIT files, see Wu and Watanabe,\n\
Bioinformatics 21:1859-1875, 2005).  The start coordinate should be\n\
less than or equal to the end coordinate.  If they are not, the program\n\
will reverse them for you.  If only a single coordinate is provided,\n\
this is equivalent to providing the same number for the start and end\n\
coordinate.\n\
\n\
The given iit file may contain types (which can be displayed by using\n\
the -T flag of the iit_dump program).  These types may be used to\n\
filter the output of the iit file.  Multiple types may be specified,\n\
which indicates a disjunctive query, such that iit_get returns entries\n\
that match any one of the given tags.\n\
\n\
If no query is provided on the command line, the program will expect\n\
one or more queries from stdin, one per line.\n\
\n\
See also: iit_store, iit_dump\n\
");
  return;
}


#ifdef __STRICT_ANSI__
int getopt (int argc, char *const argv[], const char *optstring);
#endif

static char complCode[128] = COMPLEMENT_LC;

static void
print_forward (char *sequence, int start, int end) {
  int i;

  for (i = start; i <= end; i++) {
    printf("%c",sequence[i]);
  }
  return;
}

static void
print_complement (char *sequence, int end, int start) {
  int i;

  for (i = end; i >= start; i--) {
    printf("%c",complCode[(int) sequence[i]]);
  }
  return;
}

static void
print_spaces (int n) {
  while (--n >= 0) {
    printf(" ");
  }
  return;
}


/* Need to store just the part of the query specified (e.g., 1..10) */
static void
print_interval_centered (char *divstring, unsigned int coordstart, int index, IIT_T iit, int fieldint) {
  Interval_T interval;
  char *label, *annotation, *restofheader, centerchar;
  bool allocp;
  int annotlength, left, centerpos;

  if (fieldint < 0) {
    annotation = IIT_annotation(&restofheader,iit,index,&allocp);
    if (allocp == true) {
      FREE(restofheader);
    }
  } else {
    annotation = IIT_fieldvalue(iit,index,fieldint);
    allocp = true;
  }
  annotlength = strlen(annotation);
  if (annotation[annotlength-1] == '\n') {
    annotlength--;
  }

  interval = IIT_interval(iit,index);
  left = coordstart - Interval_low(interval); /* + length(query) - queryend */
  if (Interval_sign(interval) < 0) {
    centerpos = annotlength-left-1;
  } else {
    centerpos = left;
  }
  centerchar = annotation[centerpos];

  if (centeruc == true && islower(centerchar)) {
    if (fieldint >= 0 && allocp == true) {
      FREE(annotation);
    }
  } else {
    print_spaces(centerlength-left);
    if (Interval_sign(interval) < 0) {
      print_complement(annotation,annotlength-1,centerpos+1);
      printf("[%c]",complCode[(int) centerchar]);
      print_complement(annotation,centerpos-1,0);
    } else {
      print_forward(annotation,0,centerpos-1);
      printf("[%c]",centerchar);
      print_forward(annotation,centerpos+1,annotlength-1);
    }
    print_spaces(centerlength+left-annotlength);
    if (fieldint >= 0 && allocp == true) {
      FREE(annotation);
    }
  
    printf("\t");
    if (Interval_type(interval) > 0) {
      printf("%s\t",IIT_typestring(iit,Interval_type(interval)));
    }

    if (divstring != NULL) {
      if (Interval_sign(interval) < 0) {
	printf("-%s:",divstring);
      } else {
	printf("+%s:",divstring);
      }
    }

    if (signedp == false) {
      printf("%u..%u",Interval_low(interval),Interval_high(interval));
    } else if (Interval_sign(interval) < 0) {
      printf("%u..%u",Interval_high(interval),Interval_low(interval));
    } else {
      printf("%u..%u",Interval_low(interval),Interval_high(interval));
    }
    printf("\t");

    label = IIT_label(iit,index,&allocp);
    printf("%s",label);
    if (allocp == true) {
      FREE(label);
    }
    printf("\n");
  }

  return;
}

static char *
get_total_tally (long int *tally, char *ptr) {
  int n;
  char *end;

  if ((end = index(ptr,'\n')) == NULL) {
    fprintf(stderr,"Premature end of line %s\n",ptr);
    return 0;
  }
  /* fprintf(stderr,"Getting tally for %.*s\n",end-ptr,ptr); */

  while (ptr < end) {
    while (ptr < end && !isdigit((int) *ptr)) {
      ptr++;
    }
    if (ptr < end) {
      sscanf(ptr,"%d",&n);
      (*tally) += n;
      while (ptr < end && !isspace(*ptr)) {
	ptr++;
      }
      while (ptr < end && isspace(*ptr)) {
	ptr++;
      }
    }
  }

  return ptr;
}


static void
print_line (char *ptr) {
  while (*ptr != '\0' && *ptr != '\n') {
    printf("%c",*ptr);
    ptr++;
  }
  return;
}



/* Need to store just the part of the query specified (e.g., 1..10) */
static void
print_interval_tally (unsigned int *lastcoord, char *divstring, unsigned int coordstart, unsigned int coordend,
		      int indexi, IIT_T iit, bool zeroesp) {
  Interval_T interval;
  char *annotation, *restofheader, *ptr, *nextptr;
  bool allocp;
  long int total = 0;
  unsigned int chrpos, intervalend;

  annotation = IIT_annotation(&restofheader,iit,indexi,&allocp);

  interval = IIT_interval(iit,indexi);
  chrpos = Interval_low(interval);
  intervalend = Interval_high(interval);

  ptr = annotation;

  if (zeroesp == true) {
    while (*lastcoord < chrpos) {
      printf("%s\t%u\t%d\t\n",divstring,*lastcoord,0);
      (*lastcoord)++;
    }
  }

  while (chrpos < coordstart) {
    if ((ptr = index(ptr,'\n')) == NULL) {
      fprintf(stderr,"Premature end of tally from %u to %u\n",
	      Interval_low(interval),Interval_high(interval));
      return;
    } else {
      ptr++;
    }
    chrpos++;
  }

  while (chrpos <= intervalend && chrpos <= coordend) {
    total = 0;
    nextptr = get_total_tally(&total,ptr);
    if (total > 0 || zeroesp == true) {
      printf("%s\t%u\t%ld\t",divstring,chrpos,total);
      print_line(ptr);
      printf("\n");
    }
    ptr = nextptr;
    if ((ptr = index(ptr,'\n')) == NULL) {
      fprintf(stderr,"Premature end of tally from %u to %u\n",
	      Interval_low(interval),Interval_high(interval));
      return;
    } else {
      ptr++;
    }
    chrpos++;
  }

  *lastcoord = chrpos;

  if (allocp == true) {
    FREE(restofheader);
  }
  
  return;
}


/* Need to store just the part of the query specified (e.g., 1..10) */
static void
print_interval_runlength (unsigned int *lastcoord, char *divstring, unsigned int coordstart, unsigned int coordend,
			  int indexi, IIT_T iit, bool zeroesp) {
  Interval_T interval;
  char *label;
  bool allocp;
  unsigned int chrpos, intervalend;
  int value;

  label = IIT_label(iit,indexi,&allocp);
  value = atoi(label);

  interval = IIT_interval(iit,indexi);
  chrpos = Interval_low(interval);
  intervalend = Interval_high(interval);

  if (zeroesp == true) {
    while (*lastcoord < chrpos) {
      printf("%s\t%u\t%d\n",divstring,*lastcoord,0);
      (*lastcoord)++;
    }
  }

  while (chrpos < coordstart) {
    chrpos++;
  }

  while (chrpos <= intervalend && chrpos <= coordend) {
    printf("%s\t%u\t%d\n",divstring,chrpos,value);
    chrpos++;
  }

  *lastcoord = chrpos;

  if (allocp == true) {
    FREE(label);
  }
  
  return;
}


/************************************************************************
 *   Totals
 ************************************************************************/

/* Need to store just the part of the query specified (e.g., 1..10) */
static void
compute_totals_tally (long int *total, unsigned int *n, unsigned int coordstart, unsigned int coordend,
		      int indexi, IIT_T iit) {
  Interval_T interval;
  char *annotation, *restofheader, *ptr;
  bool allocp;
  unsigned int chrpos, intervalend;

  annotation = IIT_annotation(&restofheader,iit,indexi,&allocp);

  interval = IIT_interval(iit,indexi);
  chrpos = Interval_low(interval);
  intervalend = Interval_high(interval);

  ptr = annotation;

  while (chrpos < coordstart) {
    if ((ptr = index(ptr,'\n')) == NULL) {
      fprintf(stderr,"Premature end of tally from %u to %u\n",
	      Interval_low(interval),Interval_high(interval));
      return;
    } else {
      ptr++;
    }
    chrpos++;
  }

  while (chrpos <= intervalend && chrpos <= coordend) {
    ptr = get_total_tally(&(*total),ptr);
    *n += 1;
    ptr++;
    chrpos++;
  }

  if (allocp == true) {
    FREE(restofheader);
  }
  
  return;
}


/* Need to store just the part of the query specified (e.g., 1..10) */
static double
compute_logtotal_tally (long int *total, int *n, unsigned int coordstart, unsigned int coordend,
			int indexi, IIT_T iit) {
  double logtotal = 0.0;
  Interval_T interval;
  char *annotation, *restofheader, *ptr;
  bool allocp;
  unsigned int chrpos, intervalend;
  long int count;

  annotation = IIT_annotation(&restofheader,iit,indexi,&allocp);

  interval = IIT_interval(iit,indexi);
  chrpos = Interval_low(interval);
  intervalend = Interval_high(interval);

  ptr = annotation;

  while (chrpos < coordstart) {
    if ((ptr = index(ptr,'\n')) == NULL) {
      fprintf(stderr,"Premature end of tally from %u to %u\n",
	      Interval_low(interval),Interval_high(interval));
      return logtotal;
    } else {
      ptr++;
    }
    chrpos++;
  }

  while (chrpos <= intervalend && chrpos <= coordend) {
    count = 0;
    ptr = get_total_tally(&count,ptr);

    logtotal += log((double) count + 1.0);
    *total += count;
    *n += 1;
    ptr++;
    chrpos++;
  }

  if (allocp == true) {
    FREE(restofheader);
  }
  
  return logtotal;
}



/* coordstart used only if centerp or tallyp is true */
static void
print_interval (unsigned int *lastcoord, char *divstring, unsigned int coordstart, unsigned int coordend, 
		int index, IIT_T iit, int ndivs, int fieldint) {
  Interval_T interval;
  char *label, *annotation, *restofheader;
  bool allocp;

  if (centerp == true) {
    print_interval_centered(divstring,coordstart,index,iit,fieldint);
    return;
  } else if (tallyp == true) {
    print_interval_tally(&(*lastcoord),divstring,coordstart,coordend,index,iit,zeroesp);
    return;
  } else if (runlengthp == true) {
    print_interval_runlength(&(*lastcoord),divstring,coordstart,coordend,index,iit,zeroesp);
    return;
  }

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
    if (signedp == false) {
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

  return;
}


static int *
get_matches (int *nmatches, char **divstring, unsigned int *coordstart, unsigned int *coordend,
	     int **leftflanks, int *nleftflanks, int **rightflanks, int *nrightflanks,
	     char *query, char *typestring, IIT_T *iit, char *filename) {
  int *matches;
  bool revcomp;
  int typeint;

  debug(printf("Entering get_matches with query %s.\n",query));


  if (labelp == true || Parserange_query(&(*divstring),&(*coordstart),&(*coordend),&revcomp,query,filename) == false) {
    /* Treat query as a label */
    *divstring = (char *) NULL;
    if (*iit == NULL) {
      /* Read no divs */
      if ((*iit = IIT_read(filename,/*name*/NULL,true,/*divread*/READ_NONE,/*divstring*/NULL,
			   /*add_iit_p*/true,/*labels_read_p*/true)) == NULL) {
	if (Access_file_exists_p(filename) == false) {
	  fprintf(stderr,"Cannot read file %s\n",filename);
	} else {
	  fprintf(stderr,"File %s appears to be an invalid IIT file\n",filename);
	}
	exit(9);
      }
    }
    matches = IIT_find(&(*nmatches),*iit,query);

  } else {
    if (*iit == NULL) {
      if ((*iit = IIT_read(filename,/*name*/NULL,true,/*divread*/READ_ONE,*divstring,
			   /*add_iit_p*/true,/*labels_read_p*/false)) == NULL) {
	if (Access_file_exists_p(filename) == false) {
	  fprintf(stderr,"Cannot read file %s\n",filename);
	} else {
	  fprintf(stderr,"File %s appears to be an invalid IIT file\n",filename);
	}
	exit(9);
      }
    }

    if (typestring == NULL) {
      /* Treat query as coordinates, without a typestring */
      matches = IIT_get(&(*nmatches),*iit,*divstring,*coordstart,*coordend,sortp);
      if (nflanking > 0) {
	IIT_get_flanking(&(*leftflanks),&(*nleftflanks),&(*rightflanks),&(*nrightflanks),*iit,*divstring,
			 *coordstart,*coordend,nflanking,/*sign*/0);
      }

    } else if ((typeint = IIT_typeint(*iit,typestring)) < 0) {
      fprintf(stderr,"No such type as %s.\n",typestring);
#if 0
      /* Treat query as coordinates, without a typestring */
      matches = IIT_get(&(*nmatches),*iit,*divstring,*coordstart,*coordend,sortp);
      if (nflanking > 0) {
	IIT_get_flanking(&(*leftflanks),&(*nleftflanks),&(*rightflanks),&(*nrightflanks),*iit,*divstring,
			 *coordstart,*coordend,nflanking,/*sign*/0);
      }
#else
      matches = (int *) NULL;
      nmatches = 0;
#endif      

    } else {
      /* Treat query as coordinates, with a typestring */
      matches = IIT_get_typed(&(*nmatches),*iit,*divstring,*coordstart,*coordend,typeint,sortp);
      if (nflanking > 0) {
	IIT_get_flanking_typed(&(*leftflanks),&(*nleftflanks),&(*rightflanks),&(*nrightflanks),*iit,*divstring,
			       *coordstart,*coordend,nflanking,typeint);
      }
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

  if (labelp == true || Parserange_query(&(*divstring),&(*coordstart),&(*coordend),&revcomp,query,filename) == false) {
    /* Not expecting a label */
    abort();
  }

  if ((*iit = IIT_read(filename,/*name*/NULL,true,/*divread*/READ_ONE,*divstring,/*add_iit_p*/true,
		       /*labels_read_p*/false)) == NULL) {
    if (Access_file_exists_p(filename) == false) {
      fprintf(stderr,"Cannot read file %s\n",filename);
    } else {
      fprintf(stderr,"File %s appears to be an invalid IIT file\n",filename);
    }
    exit(9);
  }

  matches = IIT_get_multiple_typed(&(*nmatches),*iit,*divstring,*coordstart,*coordend,types,ntypes,sortp);
  if (nflanking > 0) {
    IIT_get_flanking_multiple_typed(&(*leftflanks),&(*nleftflanks),&(*rightflanks),&(*nrightflanks),*iit,*divstring,
				    *coordstart,*coordend,nflanking,types,ntypes);
  }

  return matches;
}


int 
main (int argc, char *argv[]) {
  char *filename;
  char *divstring, *lasttypestring, *ptr;
  unsigned int coordstart, coordend, lastcoord = 0U;
  char Buffer[BUFLEN], nocomment[BUFLEN], query[BUFLEN], typestring[BUFLEN];
  int typeint, *types, c;
  int fieldint = -1;
  int nargs, ntypes, ndivs;
  int *matches, nmatches, i, *leftflanks, *rightflanks, nleftflanks = 0, nrightflanks = 0;
  long int total;
  unsigned int n;
  IIT_T iit = NULL;
  bool skipp;
  
  int opt;
  extern int optind;
  extern char *optarg;
  int long_option_index = 0;

  while ((opt = getopt_long(argc,argv,"f:LASUu:c:HRTZN",long_options,&long_option_index)) != -1) {
    switch (opt) {
    case 'f': fieldstring = optarg; break;
    case 'L': labelp = true; break;
    case 'A': annotationonlyp = true; break;
    case 'S': sortp = true; break;
    case 'U': signedp = false; break;
    case 'u': nflanking = atoi(optarg); break;

    case 'c': centerp = true; centerlength = atoi(optarg); break;
    case 'H': centeruc = true; break;
    case 'R': runlengthp = true; break;
    case 'T': tallyp = true; break;
    case 'Z': zeroesp = true; break;
    case 'N': statsp = true; break;

    case 'V': print_program_version(); exit(0);
    case '?': print_program_usage(); exit(0);
    default: exit(9);
    }
  }

  argc -= (optind - 1);
  argv += (optind - 1);

  if (argc <= 1) {
    fprintf(stderr,"Need to specify an iit file.  Type \"iit_get --help\" for help.\n");
    exit(9);
  } else {
    filename = argv[1];
  }

  if (statsp == true && argc == 2) {
    /* Want total over entire IIT */
    if ((iit = IIT_read(filename,NULL,true,/*divread*/READ_ALL,/*divstring*/NULL,/*add_iit_p*/true,
			/*labels_read_p*/false)) == NULL) {
      if (Access_file_exists_p(filename) == false) {
	fprintf(stderr,"Cannot read file %s\n",filename);
      } else {
	fprintf(stderr,"File %s appears to be an invalid IIT file\n",filename);
      }
      exit(9);
    }

    total = 0;
    n = 0;
    for (i = 0; i < IIT_total_nintervals(iit); i++) {
      debug(printf("index = %d\n",matches[i]));
      compute_totals_tally(&total,&n,/*coordstart*/0,/*coordend*/-1U,i,iit);
    }
    printf("counts:%ld non-zero-positions:%u mean-over-nonzero:%.3f\n",total,n,(double) total/(double) n);
      
    IIT_free(&iit);
    return 0;

  } else if (argc == 2) {

    /* Expecting input from stdin */
    if ((iit = IIT_read(filename,NULL,true,/*divread*/READ_ALL,/*divstring*/NULL,/*add_iit_p*/true,
			/*labels_read_p*/true)) == NULL) {
      if (Access_file_exists_p(filename) == false) {
	fprintf(stderr,"Cannot read file %s\n",filename);
      } else {
	fprintf(stderr,"File %s appears to be an invalid IIT file\n",filename);
      }
      exit(9);
    }

    if (fieldstring != NULL) {
      if ((fieldint = IIT_fieldint(iit,fieldstring)) < 0) {
	fprintf(stderr,"No field %s defined in iit file.\n",fieldstring);
	exit(9);
      }
    }

    while (fgets(Buffer,BUFLEN,stdin) != NULL) {
      if ((ptr = rindex(Buffer,'\n')) != NULL) {
	*ptr = '\0';
      }
      strcpy(nocomment,Buffer);
      if ((ptr = rindex(nocomment,'#')) != NULL) {
	*ptr = '\0';
      }

      skipp = false;

      if ((nargs = sscanf(nocomment,"%s %s",query,typestring)) == 2) {
	matches = get_matches(&nmatches,&divstring,&coordstart,&coordend,
			      &leftflanks,&nleftflanks,&rightflanks,&nrightflanks,
			      query,typestring,&iit,filename);

      } else if (nargs == 1) {
	matches = get_matches(&nmatches,&divstring,&coordstart,&coordend,
			      &leftflanks,&nleftflanks,&rightflanks,&nrightflanks,
			      query,/*typestring*/NULL,&iit,filename);
      } else {
	fprintf(stderr,"Can't parse line %s.  Ignoring.\n",nocomment);
	skipp = true;
      }
	
      if (skipp == false) {
	fprintf(stdout,"# Query: %s\n",Buffer);
	ndivs = IIT_ndivs(iit);
	lastcoord = coordstart;
	for (i = 0; i < nmatches; i++) {
	  print_interval(&lastcoord,divstring,coordstart,coordend,matches[i],iit,ndivs,fieldint);
	}
	if (zeroesp == true) {
	  while (lastcoord <= coordend) {
	    printf("%s\t%u\t%d\n",divstring,lastcoord,0);
	    lastcoord++;
	  }
	}
      }
      FREE(matches);
      fprintf(stdout,"# End\n");
      fflush(stdout);
    }

  } else {
    if (argc == 3) {
      /* Try as 0:<iitfile> 1:<query> */
      matches = get_matches(&nmatches,&divstring,&coordstart,&coordend,
			    &leftflanks,&nleftflanks,&rightflanks,&nrightflanks,
			    argv[2],/*typestring*/NULL,&iit,filename);
    } else if (argc == 4) {
      /* Try as 0:<iitfile> 1:<query> 2:<type> */
      matches = get_matches(&nmatches,&divstring,&coordstart,&coordend,
			    &leftflanks,&nleftflanks,&rightflanks,&nrightflanks,
			    argv[2],argv[3],&iit,filename);
    } else {
      types = (int *) CALLOC(argc-3,sizeof(int));
      for (c = 3, ntypes = 0; c < argc; c++) {
	if ((typeint = IIT_typeint(iit,argv[c])) < 0) {
	  fprintf(stderr,"No such type as %s.  Ignoring the type.\n",argv[c]);
	} else {
	  types[ntypes++] = typeint;
	  lasttypestring = argv[c];
	}
      }
      if (ntypes == 0) {
	matches = get_matches(&nmatches,&divstring,&coordstart,&coordend,
			      &leftflanks,&nleftflanks,&rightflanks,&nrightflanks,
			      argv[2],/*typestring*/NULL,&iit,filename);
      } else if (ntypes == 1) {
	matches = get_matches(&nmatches,&divstring,&coordstart,&coordend,
			      &leftflanks,&nleftflanks,&rightflanks,&nrightflanks,
			      argv[2],lasttypestring,&iit,filename);
      } else {
	matches = get_matches_multiple_typed(&nmatches,&divstring,&coordstart,&coordend,
					     &leftflanks,&nleftflanks,&rightflanks,&nrightflanks,
					     argv[2],types,ntypes,&iit,filename);
      }
    }

#if 0
    if (centerp == true) {
      print_spaces(centerlength);
      printf("*");
      print_spaces(centerlength-1);
      printf("\n");
    }
#endif

    if (statsp == true) {
      total = 0;
      n = 0;
      for (i = 0; i < nmatches; i++) {
	debug(printf("index = %d\n",matches[i]));
	compute_totals_tally(&total,&n,coordstart,coordend,matches[i],iit);
      }
      n = coordend - coordstart + 1;
      printf("counts:%ld width:%u mean:%.3f\n",total,n,(double)total/(double) n);

#if 0
    } else if (geomeanp == true) {
      logtotal = 0.0;
      total = 0;
      n = 0;
      for (i = 0; i < nmatches; i++) {
	debug(printf("index = %d\n",matches[i]));
	logtotal = compute_logtotal_tally(&total,&n,coordstart,coordend,matches[i],iit);
      }
      printf("geomean:%f totalcounts:%ld posrange:%d\n",
	     exp(logtotal/(double) (coordend - coordstart + 1)) - 1.0,total,n);
#endif

    } else {
      ndivs = IIT_ndivs(iit);
      if (nflanking > 0) {
	for (i = nleftflanks-1; i >= 0; i--) {
	  print_interval(&lastcoord,divstring,coordstart,coordend,leftflanks[i],iit,ndivs,fieldint);
	}
	printf("====================\n");
	FREE(leftflanks);
      }

      lastcoord = coordstart;
      for (i = 0; i < nmatches; i++) {
	debug(printf("index = %d\n",matches[i]));
	print_interval(&lastcoord,divstring,coordstart,coordend,matches[i],iit,ndivs,fieldint);
      }
      
      if (nflanking > 0) {
	printf("====================\n");
	for (i = 0; i < nrightflanks; i++) {
	  print_interval(&lastcoord,divstring,coordstart,coordend,rightflanks[i],iit,ndivs,fieldint);
	}
	FREE(rightflanks);
      }
      
      FREE(matches);
    }
  }

  if (divstring != NULL) {
    FREE(divstring);
  }

  IIT_free(&iit);

  return 0;
}
