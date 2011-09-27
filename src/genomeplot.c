static char rcsid[] = "$Id: genomeplot.c,v 1.10 2005/06/23 22:47:01 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mem.h"
#include "datadir.h"
#include "iit-read.h"
#include "chrsubset.h"
#include "uintlist.h"
#include "doublelist.h"
#include "plotdata.h"
#include "chrsegment.h"
#include "getopt.h"

#define SHORTAXIS 8.5*72
#define LONGAXIS 11*72

static char *user_genomedir = NULL;
static char *dbroot = NULL;

static bool logp = false;

static char *user_chrsubsetfile = NULL;
static char *user_chrsubsetname = NULL;

static char *title = NULL;

static bool segmentp = false;
static bool outputsegmentsp = false;
static bool skipvaluesp = false;
static bool gifp = false;
static bool landscapep = false;
static bool logop = false;
static double valuefactor = 2.0;

static struct option long_options[] = {
  /* Input options */
  {"dir", required_argument, 0, 'D'},	/* user_genomedir */
  {"db", required_argument, 0, 'd'}, /* dbroot */
  {"chrsubsetfile", required_argument, 0, 'C'}, /* user_chrsubsetfile */
  {"chrsubset", required_argument, 0, 'c'}, /* user_chrsubsetname */

  /* Output options */
  {"title", required_argument, 0, 't'}, /* title */
  {"log", no_argument, 0, 'L'}, /* logp */
  {"mag", required_argument, 0, 'm'}, /* valuefactor */
  {"segment", no_argument, 0, 'S'}, /* segmentp */
  {"outputsegs", no_argument, 0, 'O'}, /* oututsegmentsp */
  {"novalues", no_argument, 0, 'v'}, /* skipvaluesp */

  /* Help options */
  {"help", no_argument, 0, '?'}, /* print_program_usage */
  {0, 0, 0, 0}
};


static void
print_program_usage () {
  fprintf(stdout,"\
Usage: genomeplot [OPTIONS...] -d genome file, or\n\
       cat file | genomeplot [OPTIONS...] -d genome\n\
\n\
Input options\n\
  -D, --dir=STRING               Data directory\n\
  -d, --db=STRING                Database (e.g., NHGD)\n\
  -C, --chrsubsetfile=filename   User-suppled chromosome subset file\n\
  -c, --chrsubset=string         Chromosome subset to search\n\
\n\
Output options\n\
  -t, --title=STRING             Title\n\
  -L, --log                      Log scale\n\
  -m, --mag=DOUBLE               Magnification factor (default 2)\n\
  -S, --segment                  Perform segmentation\n\
  -O, --outputsegs               Output segments, rather than graph\n\
  -v, --novalues                 Show segments only\n\
\n\
");
  return;
}


static void
compute_page_size (double *top, double *bottom, double *left, double *right,
		   double *annotheight, double *annotwidth, 
		   bool gifp, int paperxdim, int paperydim) {

  if (gifp == true) {
    *top = paperydim - 100.0;
    *bottom = 100.0;
    *left = 50.0;
    *right = paperxdim - 20.0;
    *annotheight = 65.0;
    *annotwidth = paperxdim - 450.0; /* Smaller space for annotation */
  
  } else {
    *top = paperydim - 100.0;	/* Leave room for annot */
    *bottom = 80.0;
    *left = 50.0;
    *right = paperxdim - 50.0;
    *annotheight = 65.0;
    *annotwidth = paperxdim - 200.0;
  }

  return;
}

static void
compute_page_size_landscape (double *top, double *bottom, double *left, double *right,
			     double *annotheight, double *annotwidth, bool gifp) {
  compute_page_size(&(*top),&(*bottom),&(*left),&(*right),&(*annotheight),&(*annotwidth),gifp,
		    /*paperxdim*/SHORTAXIS,/*paperydim*/LONGAXIS);

  return;
}

static void
compute_page_size_portrait (double *top, double *bottom, double *left, double *right,
			    double *annotheight, double *annotwidth, bool gifp) {
  compute_page_size(&(*top),&(*bottom),&(*left),&(*right),&(*annotheight),&(*annotwidth),gifp,
		    /*paperxdim*/LONGAXIS,/*paperydim*/SHORTAXIS);

  return;
}


static void
print_header () {

  printf("%%!PS-Adobe-2.0\n");
  printf("%%%%Pages: (atend)\n");
  printf("%%%%BoundingBox: 0 0 612 792\n");
  printf("%%%%Orientation: Portrait\n");
  printf("%%%%DocumentFonts: Genentech Helvetica\n");
  printf("%%%%EndComments\n");

  printf("%%%%BeginProcSet\n");
  /* Code to allow pdfmark statements */
  printf("/pdfmark where {pop} {userdict /pdfmark /cleartomark load put} ifelse\n");
  printf("/languagelevel where {pop languagelevel}{1} ifelse\n");
  printf("2 lt {\n");
  printf("  userdict (<<) cvn ([) cvn load put\n");
  printf("  userdict (>>) cvn (]) cvn load put\n");
  printf("} if\n");

  /* Landscape code for PDF: */
  /*
  printf("<</EndPage\n");
  printf("{exch pop 2 ne\n");
  printf("{[{ThisPage}<</Rotate 90>>/PUT pdfmark true}\n");
  printf("{false}ifelse\n");
  printf("}\n");
  printf(">>setpagedevice\n");
  */
  /* << /PageSize [612 792] /ImagingBBox null /Orientation 1 >> setpagedevice */

  if (logop) {
#if 0
    Logo_procset(stdout);
#endif
  }
  printf("%%%%EndProcSet\n");
  if (logop) {
#if 0
    Logo_font(stdout);
#endif
  }
  print_linebreak_alg();
  printf("%%%%EndProlog\n");
  printf("%%%%BeginSetup\n");
  printf("%%%%EndSetup\n");

  return;
}


int
main (int argc, char *argv[]) {
  FILE *fp;
  char *filename = NULL;
  char *genomesubdir = NULL, *dbversion = NULL, *fileroot = NULL;
  char *chromosome;
  Plotdata_T plotdata;
  IIT_T chromosome_iit;
  Chrsubset_T chrsubset = NULL;
  unsigned int *chrpositions;
  double *values;
  int nvalues;
  int newc, oldc;

  Uintlist_T *segment_startpositions = NULL, *segment_endpositions = NULL, p, q;
  Doublelist_T *segment_means = NULL, r;

  int nincluded;
 
  unsigned int mincoord = 0U, maxcoord = 0U;
  double top, bottom, left, right, annotheight, annotwidth;

  int opt;
  extern int optind;
  extern char *optarg;
  int long_option_index = 0;

  while ((opt = getopt_long(argc,argv,"D:d:C:c:t:Lm:SOv?",
			    long_options, &long_option_index)) != -1) {
    switch (opt) {
    case 'D': user_genomedir = optarg; break;
    case 'd': 
      dbroot = (char *) CALLOC(strlen(optarg)+1,sizeof(char));
      strcpy(dbroot,optarg);
      break;
    case 'C': user_chrsubsetfile = optarg; break;
    case 'c': user_chrsubsetname = optarg; break;
    case 't': title = optarg; break;
    case 'L': logp = true; break;
    case 'm': valuefactor = atof(optarg); break;
    case 'S': segmentp = true; break;
    case 'O': outputsegmentsp = true; break;
    case 'v': skipvaluesp = true; break;

    case '?': print_program_usage(); exit(0);
    default: exit(9);
    }
  }

  argc -= (optind - 1);
  argv += (optind - 1);

  if (dbroot == NULL) {
    print_program_usage();
    exit(9);
  } else {
    genomesubdir = Datadir_find_genomesubdir(&fileroot,&dbversion,user_genomedir,dbroot);
  }

  if (argc <= 1) {
    fp = stdin;
  } else {
    filename = argv[1];
    fp = fopen(filename,"r");
    if (!fp) {
      fprintf(stderr,"Can't open %s\n",filename);
      exit(9);
    }
  }

  /* Read genome files */
  filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+strlen(fileroot)+
			     strlen(".chromosome.iit")+1,sizeof(char));
  sprintf(filename,"%s/%s.chromosome.iit",genomesubdir,fileroot);
  if ((chromosome_iit = IIT_read(filename,NULL,true)) == NULL) {
    fprintf(stderr,"Can't read IIT file %s\n",filename);
    exit(9);
  }
  FREE(filename);
  
  chrsubset = Chrsubset_read(user_chrsubsetfile,genomesubdir,fileroot,user_chrsubsetname,
			     chromosome_iit);

  FREE(fileroot);
  FREE(dbversion);
  FREE(genomesubdir);

  /* Read data */
  plotdata = Plotdata_read(fp,chromosome_iit,chrsubset);
  fprintf(stderr,"Read %d items\n",Plotdata_ngenes(plotdata));

  nincluded = Chrsubset_nincluded(chrsubset,chromosome_iit);
  for (newc = 1; newc <= nincluded; newc++) {
    oldc = Chrsubset_oldindex(chrsubset,newc);
    if (IIT_length(chromosome_iit,oldc) > maxcoord) {
      maxcoord = IIT_length(chromosome_iit,oldc);
    }
  }
  fprintf(stderr,"Max coord is %u\n",maxcoord);

  /* Compute segments */
  segment_startpositions = (Uintlist_T *) CALLOC(nincluded,sizeof(Uintlist_T));
  segment_endpositions = (Uintlist_T *) CALLOC(nincluded,sizeof(Uintlist_T));
  segment_means = (Doublelist_T *) CALLOC(nincluded,sizeof(Doublelist_T));

  for (newc = 1; newc <= nincluded; newc++) {
    segment_startpositions[newc-1] = (Uintlist_T) NULL;
    segment_endpositions[newc-1] = (Uintlist_T) NULL;
    segment_means[newc-1] = (Doublelist_T) NULL;
    if (segmentp == true) {
      if ((nvalues = Plotdata_nvalues(plotdata,newc)) > 0) {
	values = Plotdata_values(plotdata,newc);
	chrpositions = Plotdata_chrpositions(plotdata,newc);
	oldc = Chrsubset_oldindex(chrsubset,newc);
	chromosome = IIT_label(chromosome_iit,oldc);
	fprintf(stderr,"Computing segments for chromosome %s\n",chromosome);
	Chrsegment_compute(&(segment_startpositions[newc-1]),&(segment_endpositions[newc-1]),
			   &(segment_means[newc-1]),values,chrpositions,nvalues);
	if (outputsegmentsp == true) {
	  for (p = segment_startpositions[newc-1], q = segment_endpositions[newc-1], r = segment_means[newc-1];
	       p != NULL; p = Uintlist_next(p), q = Uintlist_next(q), r = Doublelist_next(r)) {
	    printf("%s:%u..%u\t%f\n",chromosome,Uintlist_head(p),Uintlist_head(q),Doublelist_head(r));
	  }
	}
      }
    }
  }


  if (outputsegmentsp == false) {
    if (landscapep == true) {
      compute_page_size_landscape(&top,&bottom,&left,&right,&annotheight,&annotwidth,gifp);
    } else {
      compute_page_size_portrait(&top,&bottom,&left,&right,&annotheight,&annotwidth,gifp);
    }

    print_header();
    Plotdata_one_signature(plotdata,0,segment_startpositions,segment_endpositions,segment_means,
			   title,logp,chromosome_iit,chrsubset,mincoord,maxcoord,
			   nincluded,valuefactor,top,bottom,left,right,annotheight,annotwidth,
			   segmentp,skipvaluesp);

  }

  for (newc = 1; newc <= nincluded; newc++) {
    Uintlist_free(&(segment_startpositions[newc-1]));
    Uintlist_free(&(segment_endpositions[newc-1]));
    Doublelist_free(&(segment_means[newc-1]));
  }
  FREE(segment_means);
  FREE(segment_endpositions);
  FREE(segment_startpositions);

  Plotdata_free(&plotdata);
  Chrsubset_free(&chrsubset);
  IIT_free_mmapped(&chromosome_iit);

  if (filename != NULL) {
    fclose(fp);
  }

  return 0;
}

int
print_linebreak_alg () {
  printf("\
/wordbreak ( ) def\n\
/BreakIntoLines \n\
{ /proc exch def \n\
/linewidth exch def \n\
/textstring exch def \n\
/breakwidth wordbreak stringwidth pop def \n\
/curwidth 0 def \n\
/lastwordbreak 0 def \n\
/startchar 0 def \n\
/restoftext textstring def \n\
{ restoftext wordbreak search \n\
{ /nextword exch def pop \n\
/restoftext exch def\n\
/wordwidth nextword stringwidth pop def \n\
curwidth wordwidth add linewidth gt \n\
{ textstring startchar lastwordbreak startchar sub \n\
getinterval proc \n\
/startchar lastwordbreak def \n\
/curwidth wordwidth breakwidth add def \n\
} \n\
{ /curwidth curwidth wordwidth add \n\
breakwidth add def \n\
} ifelse \n\
/lastwordbreak lastwordbreak nextword length add 1 add def \n\
} \n\
{ pop exit } \n\
ifelse \n\
} loop \n\
/lastchar textstring length def \n\
textstring startchar lastchar startchar sub \n\
getinterval proc \n\
} def \n");

  return 0;
}


