static char rcsid[] = "$Id: gmapindex.c,v 1.90 2005/03/09 19:25:50 twu Exp $";
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
#include <ctype.h>
#include <string.h>
#include <strings.h>		/* For rindex */
#include "bool.h"
#include "mem.h"
#include "table.h"
#include "tableint.h"
#include "genomicpos.h"
#include "compress.h"
#include "chrom.h"
#include "segmentpos.h"
#include "iit-write.h"
#include "iit-read.h"
#include "genome-write.h"
#include "indexdb.h"

#define BUFFERSIZE 8192

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


typedef enum{NONE,AUXFILES,GENOME,COMPRESS,UNCOMPRESS,OFFSETS,POSITIONS} Action_T;


/************************************************************************
 *   Creating aux file
 ************************************************************************/

/* accsegmentpos_table: char *accession -> Segmentpos_T segmentpos
   chrlength_table:     Chrom_T chrom -> Genomicpos_T chrlength
*/

static void
chrlength_update (Tableint_T chrlength_table, Chrom_T chrom, Genomicpos_T segend) {
  Genomicpos_T oldsegend;

  if ((oldsegend = (Genomicpos_T) Tableint_get(chrlength_table,chrom)) == (Genomicpos_T) 0) {
    /* Initial entry for this chromosome */
    Tableint_put(chrlength_table,chrom,segend);

  } else if (segend > oldsegend) {
    /* Revise */
    Tableint_put(chrlength_table,chrom,segend);
  }
  return;
}

static void
store_accession (Table_T accsegmentpos_table, Tableint_T chrlength_table,
		 char *accession, char *chr_string, Genomicpos_T chrpos1, 
		 Genomicpos_T chrpos2, bool revcompp, Genomicpos_T seglength, 
		 int contigtype) {
  int ret;
  Chrom_T chrom;
  Genomicpos_T segend;
  Segmentpos_T segmentpos;
  char *p;
  int subaccessioni;

  chrom = Chrom_from_string(chr_string);

  segmentpos = Segmentpos_new(chrom,chrpos1,chrpos2,revcompp,seglength,contigtype);
  Table_put(accsegmentpos_table,(void *) accession,(void *) segmentpos);

  /* Update chrlength */
  if (chrpos2 > chrpos1 + seglength) {
    chrlength_update(chrlength_table,chrom,chrpos2);
  } else {
    chrlength_update(chrlength_table,chrom,chrpos1+seglength);
  }

  return;
}


static int
count_sequence () {
  int c;
  int seglength = 0;

  while ((c = getc(stdin)) != EOF) {
    if (c == '>') {
      return seglength;
    } else if (c == '\n') {
      c = getc(stdin);
      if (c == EOF || c == '>') {
	return seglength;
      } else {
	seglength++;
      }
    } else {
      seglength++;
    }
  }
  return seglength;
}


static bool
process_sequence_aux (List_T *contigtypelist, Table_T accsegmentpos_table, Tableint_T contigtype_table, 
		      Tableint_T chrlength_table) {
  char Buffer[BUFFERSIZE], accession_p[BUFFERSIZE], *accession, 
    chrpos_string[100], *chr_string, *coords, *typestring, *p, *ptr;
  Genomicpos_T chrpos1, chrpos2, lower, upper, seglength;
  int contigtype;
  bool revcompp;

  /* Store sequence info */
  if (fgets(Buffer,BUFFERSIZE,stdin) == NULL) {
    return false;
  }

  if (sscanf(Buffer,"%s %s",accession_p,chrpos_string) < 2) {
    fprintf(stderr,"Can't parse line %s\n",Buffer);
    exit(1);
  } else {
    fprintf(stderr,"Processing contig %s at %s\n",accession_p,chrpos_string);
    if (!index(chrpos_string,':')) {
      fprintf(stderr,"Can't parse chromosomal coordinates %s\n",chrpos_string);
      exit(1);
    } else {
      chr_string = strtok(chrpos_string,":");
      coords = strtok(NULL,":");
      if (sscanf(coords,"%u..%u",&chrpos1,&chrpos2) == 2) {
	/* 1:3..5, one-based, inclusive => (2,5), zero-based, boundaries */
	if (chrpos1 <= chrpos2) {
	  chrpos1--;
	  revcompp = false;
	  lower = chrpos1;
	  upper = chrpos2;
	} else {
	  chrpos2--;
	  revcompp = true;
	  lower = chrpos2;
	  upper = chrpos1;
	}
      } else if (sscanf(coords,"%u",&chrpos1) == 1) {
	/* 1:3, one-based, inclusive => (3,3), zero-based, boundaries */
	revcompp = false;
	lower = upper = chrpos1;
      } else {
	fprintf(stderr,"Can't parse chromosomal coordinates %s\n",coords);
	exit(1);
      }
    }

    p = Buffer;
    while (*p != '\0' && !isspace(*p)) { p++; } /* Skip to first space */
    while (*p != '\0' && isspace(*p)) { p++; } /* Skip past first space */
    while (*p != '\0' && !isspace(*p)) { p++; } /* Skip to second space */
    while (*p != '\0' && isspace(*p)) { p++; } /* Skip past second space */

    if (*p == '\0') {
      contigtype = 0;		/* Empty type string */
    } else {
      if ((ptr = rindex(p,'\n')) != NULL) {
	while (isspace(*ptr)) { ptr--; } /* Erase empty space */
	ptr++;
	*ptr = '\0';
      }
      if ((contigtype = Tableint_get(contigtype_table,(void *) p)) == 0) {
	debug(printf("Storing type %s.\n",p));
	/* Store types as 1-based */
	contigtype = Tableint_length(contigtype_table) + 1;
	typestring = (char *) CALLOC(strlen(p)+1,sizeof(char));
	strcpy(typestring,p);
	Tableint_put(contigtype_table,(void *) typestring,contigtype);
	*contigtypelist = List_push(*contigtypelist,typestring);
      }
    }

    accession = (char *) CALLOC(strlen(accession_p)+1,sizeof(char));
    strcpy(accession,accession_p);
  }

  seglength = count_sequence();
  if (seglength != upper - lower) {
    fprintf(stderr,"%s has expected sequence length %u-%u=%u but actual length %u\n",
	    accession,upper,lower,upper-lower,seglength);
  }
  store_accession(accsegmentpos_table,chrlength_table,
		  accession,chr_string,lower,upper,revcompp,
		  seglength,contigtype);

  return true;
}


/************************************************************************
 *   Creating genome and related files
 ************************************************************************/

/* Modifies chrlength_table to store offsets, rather than chrlengths */
static void
write_chromosome_file (char *genomesubdir, char *fileroot, Tableint_T chrlength_table) {
  FILE *textfp;
  char *textfile, *iitfile, *chr_string, *typestring, emptystring[1];
  int n, i, j;
  Chrom_T *chroms, chrom;
  Genomicpos_T chroffset = 0, chrlength;
  List_T intervallist = NULL, chrtypelist = NULL, labellist = NULL, annotlist = NULL, p;
  Interval_T interval;

  emptystring[0] = '\0';

  /* Get chromosomes in order */
  chroms = (Chrom_T *) Tableint_keys(chrlength_table,0U);
  n = Tableint_length(chrlength_table);
  qsort(chroms,n,sizeof(Chrom_T),Chrom_compare);

  /* Write chromosome text file */
  textfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			     strlen(fileroot)+strlen(".chromosome")+1,sizeof(char));
  sprintf(textfile,"%s/%s.chromosome",genomesubdir,fileroot);
  if ((textfp = fopen(textfile,"w")) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",textfile);
    exit(9);
  }
  FREE(textfile);

  chrtypelist = List_push(chrtypelist,"");
  for (i = 0; i < n; i++) {
    chrlength = (Genomicpos_T) Tableint_get(chrlength_table,chroms[i]);
    chr_string = Chrom_to_string(chroms[i]);
    fprintf(stderr,"Chromosome %s has universal coordinates %u..%u\n",
	    chr_string,chroffset+1,chroffset+1+chrlength-1);

    fprintf(textfp,"%s\t%u..%u\t%u\n",
	    chr_string,chroffset+1,chroffset+chrlength,chrlength);
    intervallist = List_push(intervallist,(void *) Interval_new(chroffset,chroffset+chrlength-1U,0));
    labellist = List_push(labellist,(void *) chr_string);
    annotlist = List_push(annotlist,(void *) emptystring); /* No annotations */
    Tableint_put(chrlength_table,chroms[i],chroffset);
    chroffset += chrlength;
  }
  FREE(chroms);
  intervallist = List_reverse(intervallist);
  labellist = List_reverse(labellist);

  /* Write chromosome IIT file */
  iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			    strlen(fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
  sprintf(iitfile,"%s/%s.chromosome.iit",genomesubdir,fileroot);
  IIT_write(iitfile,intervallist,chrtypelist,labellist,annotlist,NULL);
  FREE(iitfile);

  List_free(&annotlist);

  for (p = labellist; p != NULL; p = List_next(p)) {
    chr_string = (char *) List_head(p);
    FREE(chr_string);
  }
  List_free(&labellist);

  /* chrtypelist has no dynamically allocated strings */
  List_free(&chrtypelist);

  for (p = intervallist; p != NULL; p = List_next(p)) {
    interval = (Interval_T) List_head(p);
    Interval_free(&interval);
  }
  List_free(&intervallist);

  return;
}

static Table_T current_accsegmentpos_table;

/* x is really a char ** */
static int
bysegmentpos_compare (const void *x, const void *y) {
  char *acc1 = * (char **) x;
  char *acc2 = * (char **) y;
  Segmentpos_T a = (Segmentpos_T) Table_get(current_accsegmentpos_table,(void *) acc1);
  Segmentpos_T b = (Segmentpos_T) Table_get(current_accsegmentpos_table,(void *) acc2);

  return Segmentpos_compare(&a,&b);
}

static void
write_contig_file (char *genomesubdir, char *fileroot, 
		   Table_T accsegmentpos_table, Tableint_T contigtype_table, List_T contigtypelist, 
		   Tableint_T chrlength_table) {
  FILE *textfp;
  char *textfile, *iitfile, *chr_string, seglength[24], *annot;
  void **keys;
  int *values;
  int naccessions, ntypes, i, j;
  char **accessions, **contigtypes, *accession;
  Segmentpos_T segmentpos;
  Chrom_T chrom;
  Genomicpos_T chroffset, universalpos1, universalpos2;
  List_T intervallist = NULL, labellist = NULL, annotlist = NULL, p;
  Interval_T interval;
  
  /* Get accessions in order */
  accessions = (char **) Table_keys(accsegmentpos_table,NULL);
  naccessions = Table_length(accsegmentpos_table);
  current_accsegmentpos_table = accsegmentpos_table;
  qsort(accessions,naccessions,sizeof(char *),bysegmentpos_compare);

  /* Get types in order */
  keys = Tableint_keys(contigtype_table,NULL);
  values = Tableint_values(contigtype_table,0);
  ntypes = Tableint_length(contigtype_table);
  contigtypes = (char **) CALLOC(ntypes+1,sizeof(char *)); /* Add 1 for type 0 */
  contigtypes[0] = "";
  for (j = 0; j < ntypes; j++) {
    contigtypes[values[j]] = keys[j];
  }
  FREE(values);
  FREE(keys);

  /* Write contig text file */
  textfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			     strlen(fileroot)+strlen(".contig")+1,sizeof(char));
  sprintf(textfile,"%s/%s.contig",genomesubdir,fileroot);
  if ((textfp = fopen(textfile,"w")) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",textfile);
    exit(9);
  }
  FREE(textfile);

  for (i = 0; i < naccessions; i++) {
    segmentpos = (Segmentpos_T) Table_get(accsegmentpos_table,(void *) accessions[i]);
    chrom = Segmentpos_chrom(segmentpos);
    chr_string = Chrom_to_string(chrom);
    chroffset = (Genomicpos_T) Tableint_get(chrlength_table,chrom);
    universalpos1 = chroffset + Segmentpos_chrpos1(segmentpos);
    universalpos2 = chroffset + Segmentpos_chrpos2(segmentpos);

    /* Print as 1-based, inclusive [a,b] */
    if (Segmentpos_revcompp(segmentpos) == true) {
      fprintf(textfp,"%s\t%u..%u\t%s:%u..%u\t%u",
	      accessions[i],universalpos2+1U,universalpos1,
	      chr_string,Segmentpos_chrpos2(segmentpos)+1U,Segmentpos_chrpos1(segmentpos),
	      Segmentpos_length(segmentpos));
    } else {
      fprintf(textfp,"%s\t%u..%u\t%s:%u..%u\t%u",
	      accessions[i],universalpos1+1U,universalpos2,
	      chr_string,Segmentpos_chrpos1(segmentpos)+1U,Segmentpos_chrpos2(segmentpos),
	      Segmentpos_length(segmentpos));
    }
    if (Segmentpos_type(segmentpos) > 0) {
      fprintf(textfp,"\t%s",contigtypes[Segmentpos_type(segmentpos)]);
    }

    fprintf(textfp,"\n");
    FREE(chr_string);

    /* Store as 0-based, inclusive [a,b] */
    intervallist = List_push(intervallist, 
			     (void *) Interval_new(universalpos1,universalpos2-1U,
						   Segmentpos_type(segmentpos)));
    labellist = List_push(labellist,(void *) accessions[i]);
    /* The negative sign in the annotation is the only indication that the contig 
       was reverse complement */
    if (Segmentpos_revcompp(segmentpos) == true) {
      sprintf(seglength,"-%u",Segmentpos_length(segmentpos));
    } else {
      sprintf(seglength,"%u",Segmentpos_length(segmentpos));
    }
    annot = (char *) CALLOC(strlen(seglength)+1,sizeof(char));
    strcpy(annot,seglength);
    annotlist = List_push(annotlist,(void *) annot);
  }

  FREE(contigtypes);
  FREE(accessions);
  intervallist = List_reverse(intervallist);
  /* contigtypelist = List_reverse(contigtypelist); -- Done by caller */ 
  labellist = List_reverse(labellist);
  annotlist = List_reverse(annotlist);

  /* Write contig IIT file */
  iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			    strlen(fileroot)+strlen(".contig.iit")+1,sizeof(char));
  sprintf(iitfile,"%s/%s.contig.iit",genomesubdir,fileroot);
  IIT_write(iitfile,intervallist,contigtypelist,labellist,annotlist,NULL);
  FREE(iitfile);

  List_free(&annotlist);

  /* Labels (accessions) are freed by accsegmentpos_table_gc */
  List_free(&labellist);

  for (p = intervallist; p != NULL; p = List_next(p)) {
    interval = (Interval_T) List_head(p);
    Interval_free(&interval);
  }
  List_free(&intervallist);

  return;
}

/************************************************************************/


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
  char *a = (char *) x;
  const char *p;
  
  for (p = a; *p != '\0'; p++) {
    h = (h << 5) - h + *p;
  }
  return h;
}


static void
stringlist_gc (List_T *list) {
  List_T p;
  char *string;

  for (p = *list; p != NULL; p = List_next(p)) {
    string = (char *) List_head(p);
    FREE(string);
  }
  List_free(&(*list));
  return;
}


static void
chrlength_table_gc (Tableint_T *chrlength_table) {
  /* Don't free chrom entries in table, because they are freed by Segmentpos_free */
  Tableint_free(&(*chrlength_table));
  return;
}

static void
accsegmentpos_table_gc (Table_T *accsegmentpos_table) {
  int n, i = 0;
  char *accession;
  Segmentpos_T segmentpos;
  void **keys, **values;

  n = Table_length(*accsegmentpos_table);
  keys = Table_keys(*accsegmentpos_table,NULL);
  values = Table_values(*accsegmentpos_table,NULL);
  for (i = 0; i < n; i++) {
    accession = (char *) keys[i];
    FREE(accession);
  }
  for (i = 0; i < n; i++) {
    segmentpos = (Segmentpos_T) values[i];
    Segmentpos_free(&segmentpos);
  }
  FREE(values);
  FREE(keys);
  Table_free(&(*accsegmentpos_table));
  return;
}

static char *
remove_slashes (char *buffer) {
  char *copy, *p;

  copy = (char *) CALLOC(strlen(buffer)+1,sizeof(char));
  strcpy(copy,buffer);
  p = copy;
  while (*p != '\0') {
    if (*p == '/') {
      *p = '_';
    }
    p++;
  }
  return copy;
}


int 
main (int argc, char *argv[]) {
  Action_T action = NONE;
  bool binaryp = false, uncompressedp = false, writefilep = false;
  int contigtype;
  char *refstrain = "reference";

  Table_T accsegmentpos_table;
  Tableint_T chrlength_table;
  Tableint_T contigtype_table;
  List_T altintervallist = NULL, contigtypelist = NULL, labellist = NULL;
  Uintlist_T seglength_list = NULL;
  IIT_T chromosome_iit, contig_iit, altstrain_iit;
  void **pairs;
  int c, i;
  unsigned int genomelength;
  char *sourcedir = ".", *destdir = ".", *fileroot = NULL, 
    *chromosomefile, *iitfile, *iittypefile, *offsetsfile, *positionsfile,
    *typestring, *copy;
  FILE *offsets_fp, *fp;

  extern int optind;
  extern char *optarg;

  while ((c = getopt(argc,argv,"F:D:d:AGgCUOoPps:W")) != -1) {
    switch (c) {
    case 'F': sourcedir = optarg; break;
    case 'D': destdir = optarg; break;
    case 'd': fileroot = optarg; break;
    case 'A': action = AUXFILES; break;
    case 'G': action = GENOME; uncompressedp = false; break;
    case 'g': action = GENOME; uncompressedp = true; break;
    case 'C': action = COMPRESS; break;
    case 'U': action = UNCOMPRESS; break;
    case 'O': action = OFFSETS; uncompressedp = false; break;
    case 'o': action = OFFSETS; uncompressedp = true; break;
    case 'P': action = POSITIONS; uncompressedp = false; break;
    case 'p': action = POSITIONS; uncompressedp = true; break;
    case 's': refstrain = optarg; break;
    case 'W': writefilep = true; break;
    }
  }
  argc -= (optind - 1);
  argv += (optind - 1);

  if (action != COMPRESS && action != UNCOMPRESS) {
    if (fileroot == NULL) {
      fprintf(stderr,"Missing name of genome database.  Must specify with -d flag.\n");
      exit(9);
    }
  }

  if (action == AUXFILES) {
    /* Usage: cat <fastafile> | gmapindex [-F <sourcedir>] [-D <destdir>] -d <dbname> -A [-s <refstrain>]
       Requires <fastafile> in appropriate format
       Writes <destdir>/<dbname>.chromosome and <destdir>/<dbname>.contig files 
       and corresponding .iit files */

    if (getc(stdin) != '>') {
      fprintf(stderr,"Expected file to start with '>'\n");
      exit(9);
    }

    /* keys are strings; values are structs. */
    accsegmentpos_table = Table_new(4000,string_compare,string_hash);
    /* keys are Chrom_Ts; values are uints. */
    chrlength_table = Tableint_new(40,Chrom_compare_table,Chrom_hash_table);
    /* keys are strings; values are ints */
    contigtype_table = Tableint_new(100,string_compare,string_hash);

    typestring = (char *) CALLOC(sizeof(refstrain)+1,sizeof(char));
    strcpy(typestring,refstrain);
    contigtypelist = List_push(NULL,typestring);
    while (process_sequence_aux(&contigtypelist,accsegmentpos_table,contigtype_table,
				chrlength_table) == true) {
    }

    contigtypelist = List_reverse(contigtypelist);
    write_chromosome_file(destdir,fileroot,chrlength_table);
    write_contig_file(destdir,fileroot,accsegmentpos_table,contigtype_table,contigtypelist,
		      chrlength_table);

    stringlist_gc(&contigtypelist);
    Tableint_free(&contigtype_table);
    chrlength_table_gc(&chrlength_table);
    accsegmentpos_table_gc(&accsegmentpos_table);

  } else if (action == GENOME) {
    /* Usage: cat <fastafile> | gmapindex [-F <sourcedir>] [-D <destdir>] -d <dbname> -G
       Requires <fastafile> in appropriate format and <sourcedir>/<dbname>.chromosome.iit 
       and <sourcedir>/<dbname>.contig.iit files.
       Creates <destdir>/<dbname>.genome and if strain information is present,
       <destdir>/<dbname>.altstrain.iit */

    chromosomefile = (char *) CALLOC(strlen(sourcedir)+strlen("/")+
				     strlen(fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
    sprintf(chromosomefile,"%s/%s.chromosome.iit",sourcedir,fileroot);
    chromosome_iit = IIT_read(chromosomefile,NULL,true);
    genomelength = IIT_totallength(chromosome_iit);
    FREE(chromosomefile);

    iitfile = (char *) CALLOC(strlen(sourcedir)+strlen("/")+
			      strlen(fileroot)+strlen(".contig.iit")+1,sizeof(char));
    sprintf(iitfile,"%s/%s.contig.iit",sourcedir,fileroot);
    contig_iit = IIT_read(iitfile,NULL,true);
    FREE(iitfile);

    if (IIT_ntypes(contig_iit) == 1) {
      Genome_write(destdir,fileroot,stdin,contig_iit,NULL,uncompressedp,writefilep,genomelength);
    } else if (IIT_ntypes(contig_iit) > 1) {
      iitfile = (char *) CALLOC(strlen(destdir)+strlen("/")+
				strlen(fileroot)+strlen(".altstrain.iit")+1,sizeof(char));
      sprintf(iitfile,"%s/%s.altstrain.iit",destdir,fileroot);
      altintervallist = IIT_intervallist_typed(&labellist,&seglength_list,contig_iit);
      contigtypelist = IIT_typelist(contig_iit);

      /* For altstrain_iit, we provide seglength now, and write actual sequence later
	 (in Genome_write) */
      fprintf(stderr,"Writing alternate strain file...\n");
      IIT_write(iitfile,altintervallist,contigtypelist,labellist,
		/*annotlist*/(List_T) NULL,seglength_list);
      fprintf(stderr,"Done writing alternate strain file\n");

      altstrain_iit = IIT_read(iitfile,NULL,false);
      Genome_write(destdir,fileroot,stdin,contig_iit,altstrain_iit,uncompressedp,writefilep,genomelength);
      FREE(iitfile);

      /* Write .altstrain.type file */
      iittypefile = (char *) CALLOC(strlen(destdir)+strlen("/")+
				    strlen(fileroot)+strlen(".altstrain.type")+1,sizeof(char));
      sprintf(iittypefile,"%s/%s.altstrain.type",destdir,fileroot);
      fp = fopen(iittypefile,"w");
      if (fp != NULL) {
	IIT_dump_typestrings(fp,altstrain_iit);
	fclose(fp);
      }
      FREE(iittypefile);

      IIT_free_mmapped(&altstrain_iit);
    }

    IIT_free_mmapped(&contig_iit);

  } else if (action == COMPRESS) {
    /* Usage: cat <genomefile> | gmapindex -C > <genomecompfile> */

    Compress_compress(stdin);

  } else if (action == UNCOMPRESS) {
    /* Usage: cat <genomecompfile> | gmapindex -U > <genomefile> */
    
    Compress_uncompress(stdin);

  } else if (action == OFFSETS) {
    /* Usage: cat <genomefile> | gmapindex [-F <sourcedir>] [-D <destdir>] -d <dbname> -O
       If alternate strains are present, requires <sourcedir>/<dbname>.altstrain.iit
       Creates <destdir>/<dbname>.idxoffsets */

    /* Reference strain */
    offsetsfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
				  strlen(".idxoffsets")+1,sizeof(char));
    sprintf(offsetsfile,"%s/%s.idxoffsets",destdir,fileroot);
    if ((offsets_fp = fopen(offsetsfile,"w")) == NULL) {
      fprintf(stderr,"Can't write to file %s\n",offsetsfile);
      exit(9);
    }

    iitfile = (char *) CALLOC(strlen(sourcedir)+strlen("/")+
			      strlen(fileroot)+strlen(".altstrain.iit")+1,sizeof(char));
    sprintf(iitfile,"%s/%s.altstrain.iit",sourcedir,fileroot);
    altstrain_iit = IIT_read(iitfile,NULL,true);

    Indexdb_write_offsets(offsets_fp,stdin,altstrain_iit,uncompressedp);

    IIT_free_mmapped(&altstrain_iit);
    FREE(iitfile);
    fclose(offsets_fp);
    FREE(offsetsfile);

  } else if (action == POSITIONS) {
    /* Usage: cat <genomefile> | gmapindex [-F <sourcedir>] [-D <destdir>] -d <dbname> -P
       Requires <sourcedir>/<dbname>.idxoffsets, and possibly <sourcedir>/<dbname>.altstrain.iit
       Creates <destdir>/<dbname>.idxpositions */

    /* Reference strain */
    offsetsfile = (char *) CALLOC(strlen(sourcedir)+strlen("/")+strlen(fileroot)+
				  strlen(".idxoffsets")+1,sizeof(char));
    sprintf(offsetsfile,"%s/%s.idxoffsets",sourcedir,fileroot);
    if ((offsets_fp = fopen(offsetsfile,"r")) == NULL) {
      fprintf(stderr,"Can't open file %s\n",offsetsfile);
      exit(9);
    }
    positionsfile = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
				    strlen(".idxpositions")+1,sizeof(char));
    sprintf(positionsfile,"%s/%s.idxpositions",destdir,fileroot);

    iitfile = (char *) CALLOC(strlen(sourcedir)+strlen("/")+
			      strlen(fileroot)+strlen(".altstrain.iit")+1,sizeof(char));
    sprintf(iitfile,"%s/%s.altstrain.iit",sourcedir,fileroot);
    altstrain_iit = IIT_read(iitfile,NULL,true);

    Indexdb_write_positions(positionsfile,offsets_fp,stdin,altstrain_iit,uncompressedp,writefilep);

    IIT_free_mmapped(&altstrain_iit);
    FREE(iitfile);
    fclose(offsets_fp);
    FREE(positionsfile);
    FREE(offsetsfile);
  }

  return 0;
}

