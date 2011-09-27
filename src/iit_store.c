static char rcsid[] = "$Id: iit_store.c,v 1.30 2007/09/18 20:55:09 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif
#include <string.h>		/* For strlen */
#include <strings.h>		/* For rindex */
#include <ctype.h>
#include "bool.h"
#include "assert.h"
#include "mem.h"
#include "fopen.h"

#include "list.h"
#include "interval.h"
#include "tableint.h"
#include "iit-write.h"
#include "getopt.h"

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

#define LINELENGTH 8192

/************************************************************************
 *   Program options
 ************************************************************************/

static char *outputfile = NULL;
static bool gff3_format_p = false;
static char *labelid = "ID";
static bool fieldsp = false;
static char iit_version = 0;	/* Means latest version */

static struct option long_options[] = {
  /* Input options */
  {"output", required_argument, 0, 'o'}, /* outputfile */
  {"fields", no_argument, 0, 'F'}, /* fieldsp */
  {"gff", no_argument, 0, 'G'}, /* gff3_format_p */
  {"label", required_argument, 0, 'l'}, /* labelid */
  {"iitversion", required_argument, 0, 'v'}, /* iit_version */

  /* Help options */
  {"version", no_argument, 0, 'V'}, /* print_program_version */
  {"help", no_argument, 0, '?'}, /* print_program_usage */
  {0, 0, 0, 0}
};

static void
print_program_version () {
  fprintf(stdout,"\n");
  fprintf(stdout,"iit_store: indexing utility for Interval Index Trees\n");
  fprintf(stdout,"Part of GMAP package, version %s\n",PACKAGE_VERSION);
  fprintf(stdout,"Thomas D. Wu, Genentech, Inc.\n");
  fprintf(stdout,"Contact: twu@gene.com\n");
  fprintf(stdout,"\n");
  return;
}

static void
print_program_usage () {
  fprintf(stdout,"\
Usage: iit_store [OPTIONS...] -o outputfile inputfile, or\n\
       cat inputfile | iit_store [OPTIONS...] -o outputfile\n\
where\n\
   outputfile is the desired filename for the iit file\n\
       (.iit will be added as a suffix if necessary), and\n\
   inputfile is in either FASTA or GFF3 format, as described below.\n\
\n\
Options\n\
  -o, --output=STRING       Name of output iit file\n\
  -F, --fields              Annotation consists of separate fields\n\
  -G, --gff                 Parse input file in gff3 format\n\
  -l, --label=STRING        For gff input, the feature attribute to use (default is ID)\n\
  -v, --iitversion=STRING   Desired iit version for output iit\n\
                            (default = 0, which means latest version)\n\
\n\
  -V, --version             Show version\n\
  -?, --help                Show this help message\n\
\n\
\n\
Description of input format:\n\
\n\
The FASTA format for input files should be\n\
\n\
    >label start end optional_tag\n\
    optional_annotation (which may be zero, one, or multiple lines)\n\
\n\
For example, the label may be a sequence accession, with the start and end\n\
numbers representing chromosomal coordinates, and the tag indicating the\n\
chromosome and strand, such as '17'.\n\
\n\
In version 2 and later of IIT files, intervals may have directions.\n\
To indicate a forward direction, the start coordinate should be less\n\
than the end coordinate.  To indicate a reverse direction,\n\
the start coordinate should be greater than the end coordinate.\n\
If they are the same, then no direction is implied.\n\
\n\
A header would therefore look like\n\
\n\
    >NM_004448 35138441 35109780 17\n\
\n\
which indicates an interval on tag (i.e., chromosome) 17 in the reverse\n\
direction.\n\
\n\
If the -F flag is provided, IIT files may store annotation for each interval\n\
as separate fields.  The input must contain the names of the fields, one per\n\
line, before the first interval header.  Each interval then contains annotation\n\
corresponding to each field, one value per line.\n\
\n\
The GFF3 format requires the -G flag and optionally the -l flag.\n\
The iit_store program will parse the chromosome from column 1, the start\n\
coordinate from column 4, the end coordinate from column 5, the strand\n\
from column 7, an if possible, the label from column 9.  The -l flag\n\
will indicate which feature from column 9 to retrieve, such as ID, Name,\n\
or Parent.  Appropriate choice of label may be helpful later on, because\n\
the iit_get program can retrieve information by label, as well as by\n\
coordinates.\n\
\n\
Limitations: Start and end coordinates must be non-negative integers, and are\n\
limited to the domain of an unsigned int.  For machines where unsigned ints\n\
are 4 bytes, this means coordinates must be less than 2^32 = 4294967296.\n\
\n\
See also: iit_get, iit_dump\n\
");
  return;
}

/* Empties contents of lines */
static char *
concatenate_lines (List_T lines, int content_size) {
  char *string, *temp;
  List_T l;

  string = (char *) CALLOC(content_size+1,sizeof(char));
  for (l = lines; l; l = List_next(l)) {
    temp = (char *) List_head(l);
    strcat(string,temp);
    FREE(temp);
  }
  
  /* Keep last return
  if (string[content_size-1] == '\n') {
    string[content_size-1] = '\0';
  }
  */

  return string;
}


static int
string_compare (const void *x, const void *y) {
  char *a = (char *) x;
  char *b = (char *) y;

  return strcmp(a,b);
}

/* This is the X31 hash function */
static unsigned int
string_hash (const void *x) {
  unsigned int h = 0U;
  const char *p;
  
  for (p = x; *p != '\0'; p++) {
    h = (h << 5) - h + *p;
  }
  return h;
}


static List_T
scan_header (List_T typelist, Tableint_T typetable, 
	     char **label, unsigned int *start, unsigned int *end, int *type,
	     char *header) {
  char Buffer[1024], *typestring, *p, *ptr;

  /* Example: >A 1 10 FWD */
  if (sscanf(header,">%s %u %u\n",Buffer,&(*start),&(*end)) < 3) {
    fprintf(stderr,"Error parsing %s.  Expecting a FASTA type header with a label, two coordinates, and optional tag.\n",header);
    exit(9);
  } else {
    *label = (char *) CALLOC(strlen(Buffer)+1,sizeof(char));
    strcpy(*label,Buffer);

    p = header;
    while (!isspace((int) *p)) { p++; } /* First word */
    while (isspace((int) *p)) { p++; } /* First space */
    while (!isspace((int) *p)) { p++; } /* Second word */
    while (isspace((int) *p)) { p++; } /* Second space */
    while (!isspace((int) *p)) { p++; } /* Third word */
    while (*p != '\0' && isspace((int) *p)) { p++; } /* Third space */
    
    if (*p == '\0') {
      *type = 0;		/* Empty type string */
    } else {
      if ((ptr = rindex(p,'\n')) != NULL) {
	while (isspace((int) *ptr)) { ptr--; } /* Erase empty space */
	ptr++;
	*ptr = '\0';
      }
      if ((*type = Tableint_get(typetable,(void *) p)) == 0) {
	/* Store types as 1-based */
	*type = Tableint_length(typetable) + 1;
	typestring = (char *) CALLOC(strlen(p)+1,sizeof(char));
	strcpy(typestring,p);
	Tableint_put(typetable,typestring,*type);
	typelist = List_push(typelist,typestring);
	/* debug(printf("Entering new type %s.\n",typestring)); */
      }
    }
  }
  return typelist;
}

static List_T
parse_fieldlist (char *firstchar, FILE *fp) {
  List_T fieldlist = NULL;
  char Buffer[LINELENGTH], *fieldname, *p;

  while ((*firstchar = fgetc(fp)) != '>') {
    Buffer[0] = *firstchar;
    fgets(&(Buffer[1]),LINELENGTH-1,fp);
    if ((p = rindex(Buffer,'\n')) == NULL) {
      fprintf(stderr,"Line %s exceeds maximum length of %d\n",Buffer,LINELENGTH);
      exit(9);
    } else {
      *p = '\0';
    }
    fieldname = (char *) CALLOC(strlen(Buffer)+1,sizeof(char));
    strcpy(fieldname,Buffer);
    fieldlist = List_push(fieldlist,fieldname);
  }

  return List_reverse(fieldlist);
}


static void
parse_fasta (List_T *intervallist, List_T *typelist, List_T *labellist, List_T *annotlist,
	     FILE *fp, Tableint_T typetable, char firstchar) {
  char Buffer[LINELENGTH], *label, *tempstring;
  unsigned int start, end;
  List_T lines;
  int content_size, type;

  if (firstchar == '\0') {
    fgets(Buffer,LINELENGTH,fp);
  } else {
    Buffer[0] = firstchar;
    fgets(&(Buffer[1]),LINELENGTH-1,fp);
  }
  *typelist = scan_header(*typelist,typetable,&label,&start,&end,&type,Buffer);
  *labellist = List_push(NULL,label);
  lines = NULL;
  content_size = 0;

  while (fgets(Buffer,LINELENGTH,fp) != NULL) {
    if (Buffer[0] == '>') {

      *intervallist = List_push(*intervallist,(void *) Interval_new(start,end,type));
      lines = List_reverse(lines);
      *annotlist = List_push(*annotlist,concatenate_lines(lines,content_size));
      List_free(&lines);

      *typelist = scan_header(*typelist,typetable,&label,&start,&end,&type,Buffer);
      *labellist = List_push(*labellist,label);
      lines = NULL;
      content_size = 0;

    } else {
      tempstring = (char *) CALLOC(strlen(Buffer)+1,sizeof(char));
      strcpy(tempstring,Buffer);
      lines = List_push(lines,(void *) tempstring);
      content_size += strlen(Buffer);
    }
  }

  *intervallist = List_push(*intervallist,(void *) Interval_new(start,end,type));
  lines = List_reverse(lines);
  *annotlist = List_push(*annotlist,concatenate_lines(lines,content_size));
  List_free(&lines);
  
  *intervallist = List_reverse(*intervallist);
  *typelist = List_reverse(*typelist);
  *labellist = List_reverse(*labellist);
  *annotlist = List_reverse(*annotlist);

  return;
}


static int
assign_columns (char **columns, char *Buffer, int maxfields) {
  char *token;
  int nfields = 0;
  
  columns[nfields++] = token = strtok(Buffer,"\t");
  while ((token = strtok(NULL,"\t")) != NULL && nfields < maxfields) {
    columns[nfields++] = token;
  }
  return nfields;
}


#define CHRCOLUMN 0
#define STARTCOLUMN 3
#define ENDCOLUMN 4
#define STRANDCOLUMN 6
#define FEATURECOLUMN 8
#define GFF3_COLUMNS 9

/* Modifies feature */
static char *
gff3_feature_id (char *feature, char *labelstr, int labelstrlen, int lineno) {
  char *token, *value, *p;

  token = strtok(feature,";");
  if (!strncmp(token,labelstr,labelstrlen)) {
    value = &(token[labelstrlen]);
    if (value[0] != '"') {
      return value;
    } else {
      value = &(value[1]);
      /* Quotation marks */
      if ((p = rindex(value,'"')) == NULL) {
	fprintf(stderr,"Error in line %d: Saw no matching quotation in %s\n",lineno,token);
	exit(9);
      } else {
	*p = '\0';
      }
      return value;
    }
  } else {
    while ((token = strtok(NULL,";")) != NULL) {
      if (!strncmp(token,labelstr,labelstrlen)) {
	value = &(token[labelstrlen]);
	if (value[0] != '"') {
	  return value;
	} else {
	  value = &(value[1]);
	  /* Quotation marks */
	  if ((p = rindex(value,'"')) == NULL) {
	    fprintf(stderr,"Error in line %d: Saw no matching quotation in %s\n",lineno,token);
	    exit(9);
	  } else {
	    *p = '\0';
	  }
	  return value;
	}
      }
    }
    return NULL;
  }
}

static bool
empty_line_p (char *line) {
  char *p = line;

  while (*p != '\0' && isspace(*p)) {
    p++;
  }
  if (*p == '\0') {
    return true;
  } else {
    return false;
  }
}

static void
parse_gff3 (List_T *intervallist, List_T *typelist, List_T *labellist, List_T *annotlist,
	    FILE *fp, Tableint_T typetable) {
  char Buffer[LINELENGTH], Space[1000], *columns[GFF3_COLUMNS],
    *typestring, *label, *p, *chr, *line, *idptr;
  unsigned int start, end;
  int nfields, type, lineno = 0, row = 0, labelstrlen;
  char strandchar;
  char *labelstr;

  labelstr = (char *) CALLOC(strlen(labelid) + strlen("=") + 1,sizeof(char));
  sprintf(labelstr,"%s=",labelid);
  labelstrlen = strlen(labelstr);

  while (fgets(Buffer,LINELENGTH,fp) != NULL) {
    lineno++;
    if (Buffer[0] == '#') {
      /* Skip comment */
    } else if (empty_line_p(Buffer) == true) {
      /* Skip empty line */
    } else {
      line = (char *) CALLOC(strlen(Buffer)+1,sizeof(char));
      strcpy(line,Buffer);

      if ((p = rindex(Buffer,'\n')) == NULL) {
	fprintf(stderr,"Line exceeds maximum length of %d\n",LINELENGTH);
	exit(9);
      } else {
	*p = '\0';
      }
      nfields = assign_columns(columns,Buffer,GFF3_COLUMNS); /* destroys Buffer */

      if (nfields < GFF3_COLUMNS-1) {
	/* Subract 1 to allow for an empty feature column */
	fprintf(stderr,"Skipping line %d with only %d fields: %s\n",lineno,nfields,line);
	FREE(line);
      } else {
	chr = columns[CHRCOLUMN];
	sprintf(typestring,"%s",chr);
	
	if ((strandchar = columns[STRANDCOLUMN][0]) == '+') {
	  start = atof(columns[STARTCOLUMN]);
	  end = atof(columns[ENDCOLUMN]);
	} else if (strandchar == '-') {
	  start = atof(columns[ENDCOLUMN]);
	  end = atof(columns[STARTCOLUMN]);
	} else if (strandchar == '.' || strandchar == '?') {
	  start = atof(columns[STARTCOLUMN]);
	  end = atof(columns[ENDCOLUMN]);
	} else {
	  start = atof(columns[STARTCOLUMN]);
	  end = atof(columns[ENDCOLUMN]);
	}

	if ((type = Tableint_get(typetable,(void *) typestring)) != 0) {
	  FREE(typestring);
	} else {
	  type = Tableint_length(typetable) + 1;
	  Tableint_put(typetable,typestring,type);
	  *typelist = List_push(*typelist,typestring);
	}
	*intervallist = List_push(*intervallist,(void *) Interval_new(start,end,type));

	if (nfields <= FEATURECOLUMN) {
	  sprintf(Space,"gff.%d",row);
	  label = (char *) CALLOC(strlen(Space)+1,sizeof(char));
	  strcpy(label,Space);
	} else if ((idptr = gff3_feature_id(columns[FEATURECOLUMN],labelstr,labelstrlen,lineno)) == NULL) {
	  sprintf(Space,"gff.%d",row);
	  label = (char *) CALLOC(strlen(Space)+1,sizeof(char));
	  strcpy(label,Space);
	} else {
	  label = (char *) CALLOC(strlen(idptr)+1,sizeof(char));
	  strcpy(label,idptr);
	}
	*labellist = List_push(*labellist,(void *) label);
	*annotlist = List_push(*annotlist,line);

	row++;
      }
    }
  }

  *intervallist = List_reverse(*intervallist);
  *typelist = List_reverse(*typelist);
  *labellist = List_reverse(*labellist);
  *annotlist = List_reverse(*annotlist);

  FREE(labelstr);

  return;
}


#ifdef __STRICT_ANSI__
int getopt (int argc, char *const argv[], const char *optstring);
#endif

int 
main (int argc, char *argv[]) {
  char *inputfile = NULL, *iitfile, *tempstring, *typestring, *p;
  char Buffer[8192], firstchar;
  List_T lines = NULL, l, intervallist = NULL, typelist = NULL, labellist = NULL, fieldlist = NULL, annotlist = NULL;
  FILE *fp;
  Interval_T interval;
  Tableint_T typetable;

  int opt;
  extern int optind;
  extern char *optarg;
  int long_option_index = 0;

  while ((opt = getopt_long(argc,argv,"o:FGl:v:",
			    long_options,&long_option_index)) != -1) {
    switch (opt) {
    case 'o': outputfile = optarg; break;
    case 'F': fieldsp = true; break;
    case 'G': gff3_format_p = true; break;
    case 'l': labelid = optarg; break;
    case 'v': iit_version = atoi(optarg); break;

    case 'V': print_program_version(); exit(0);
    case '?': print_program_usage(); exit(0);
    default: exit(9);
    }
  }
  argc -= optind;
  argv += optind;

  if (outputfile == NULL) {
    fprintf(stderr,"Need to specify an output file with the -o flag\n");
    exit(9);
  }

  if (argc < 1) {
    fp = stdin;
  } else {
    inputfile = argv[0];
    fp = FOPEN_READ_TEXT(inputfile);
    if (!fp) {
      fprintf(stderr,"Can't open file %s\n",inputfile);
      exit(9);
    }
  }

  typetable = Tableint_new(1000,string_compare,string_hash);
  /* The zeroth type is empty */
  typestring = (char *) CALLOC(1,sizeof(char));
  typestring[0] = '\0';
  typelist = List_push(typelist,typestring);

  if (gff3_format_p == true) {
    parse_gff3(&intervallist,&typelist,&labellist,&annotlist,fp,typetable);
  } else {
    fieldlist = parse_fieldlist(&firstchar,fp);
    parse_fasta(&intervallist,&typelist,&labellist,&annotlist,fp,typetable,firstchar);
  }

  if (inputfile != NULL) {
    fclose(fp);
  }

  /* Figure out name of iit file */
  if (strlen(outputfile) < 4) {
    iitfile = (char *) CALLOC(strlen(outputfile)+strlen(".iit")+1,sizeof(char));
    sprintf(iitfile,"%s.iit",outputfile);
  } else {
    p = &(outputfile[strlen(outputfile)]);
    p -= 4;
    if (!strcmp(p,".iit")) {
      iitfile = (char *) CALLOC(strlen(outputfile)+1,sizeof(char));
      strcpy(iitfile,outputfile);
    } else {
      iitfile = (char *) CALLOC(strlen(outputfile)+strlen(".iit")+1,sizeof(char));
      sprintf(iitfile,"%s.iit",outputfile);
    }
  }

  IIT_write(iitfile,intervallist,typelist,labellist,fieldlist,annotlist,NULL,iit_version);
  FREE(iitfile);

  for (l = annotlist; l != NULL; l = List_next(l)) {
    tempstring = (char *) List_head(l);
    FREE(tempstring);
  }
  List_free(&annotlist);

  for (l = fieldlist; l != NULL; l = List_next(l)) {
    tempstring = (char *) List_head(l);
    FREE(tempstring);
  }
  List_free(&fieldlist);

  for (l = labellist; l != NULL; l = List_next(l)) {
    tempstring = (char *) List_head(l);
    FREE(tempstring);
  }
  List_free(&labellist);

  for (l = typelist; l != NULL; l = List_next(l)) {
    tempstring = (char *) List_head(l);
    FREE(tempstring);
  }
  List_free(&typelist);

  Tableint_free(&typetable);

  for (l = intervallist; l != NULL; l = List_next(l)) {
    interval = (Interval_T) List_head(l);
    Interval_free(&interval);
  }
  List_free(&intervallist);

  return 0;
}
