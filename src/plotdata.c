static char rcsid[] = "$Id: plotdata.c,v 1.7 2005/06/23 22:47:18 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "plotdata.h"

#include <stdlib.h>
#include <string.h>
#include <strings.h>		/* For index */
#include <math.h>
#include "mem.h"
#include "datum.h"
#include "color.h"


#define T Plotdata_T
struct T {
  int nsamples;
  int ngenes;
  
  int nincluded;
  int *nvalues;
  unsigned int **chrpositions;
  double **values;
};

void
Plotdata_free (T *old) {
  int c, g;
  
  for (c = 0; c < (*old)->nincluded; c++) {
    if ((*old)->nvalues[c] > 0) {
      FREE((*old)->values[c]);
      FREE((*old)->chrpositions[c]);
    }
  }
  FREE((*old)->values);
  FREE((*old)->chrpositions);
  FREE((*old)->nvalues);
  FREE(*old);
  return;
}


int
Plotdata_nsamples (T this) {
  return this->nsamples;
}

int
Plotdata_ngenes (T this) {
  return this->ngenes;
}

int
Plotdata_nvalues (T this, int c) {
  return this->nvalues[c-1];
}


unsigned int *
Plotdata_chrpositions (T this, int c) {
  return this->chrpositions[c-1];
}

double *
Plotdata_values (T this, int c) {
  return this->values[c-1];
}


double
Plotdata_get (unsigned int *chrpos, T this, int c, int g, int s) {
  double value;
  
  *chrpos = this->chrpositions[c-1][g];
  value = this->values[c-1][g];
  return value;
}


/* Note that isnumber is a function in ctype.h on some systems */
static bool
isnumberp (unsigned int *result, char *string) {
  char *p = string;

  *result = 0U;
  while (*p != '\0') {
    if (*p == ',') {
      /* Skip commas */
    } else if (!isdigit(*p)) {
      return false;
    } else {
      *result = (*result) * 10 + (*p - '0');
    }
    p++;
  }
  return true;
}


T
Plotdata_read (FILE *fp, IIT_T chromosome_iit, Chrsubset_T chrsubset) {
  T new = (T) MALLOC(sizeof(*new));
  char Buffer[1024], position[512];
  int chrindex, newchrindex, newc, nvalues, g;
  unsigned int chrpos;
  double value;
  char *chromosome, *coords;
  List_T *datalists;
  Datum_T *array, datum;

  new->nsamples = 1;
  new->ngenes = 0;
  new->nincluded = Chrsubset_nincluded(chrsubset,chromosome_iit);
  datalists = (List_T *) CALLOC(new->nincluded,sizeof(List_T));

  while (fgets(Buffer,1024,fp) != NULL) {
    if (sscanf(Buffer,"%s %lf",position,&value) != 2) {
      fprintf(stderr,"Skipping line %s",Buffer);
    } else {
      if (index(position,':') == NULL) {
	fprintf(stderr,"Can't parse position %s\n",position);
	exit(9);
      } else {
	/* Segment must be a genome, chromosome, or contig */
	chromosome = strtok(position,":");
	if ((chrindex = IIT_find_linear(chromosome_iit,chromosome)) < 0) {
	  fprintf(stderr,"Don't recognize chromosome %s\n",chromosome);
	  exit(9);
	}
	if ((newchrindex = Chrsubset_newindex(chrsubset,chrindex)) > 0) {
	  coords = strtok(NULL,":");
	  if (isnumberp(&chrpos,coords) == false) {
	    fprintf(stderr,"Can't find number in %s\n",coords);
	    exit(9);
	  }

	  datalists[newchrindex-1] = List_push(datalists[newchrindex-1],(void *) Datum_new(chrpos,value));
	}
      }
    }
  }
  
  new->chrpositions = (unsigned int **) CALLOC(new->nincluded,sizeof(unsigned int *));
  new->values = (double **) CALLOC(new->nincluded,sizeof(double *));
  new->nvalues = (int *) CALLOC(new->nincluded,sizeof(int));
  for (newc = 0; newc < new->nincluded; newc++) {
    nvalues = new->nvalues[newc] = List_length(datalists[newc]);
    if (nvalues == 0) {
      new->chrpositions[newc] = (unsigned int *) NULL;
      new->values[newc] = (double *) NULL;
    } else {
      new->ngenes += nvalues;
      array = (Datum_T *) List_to_array(datalists[newc],NULL);
      qsort(array,nvalues,sizeof(Datum_T),Datum_cmp);
      new->chrpositions[newc] = (unsigned int *) CALLOC(nvalues,sizeof(unsigned int));
      new->values[newc] = (double *) CALLOC(nvalues,sizeof(double));
      for (g = 0; g < nvalues; g++) {
	datum = array[g];
	new->chrpositions[newc][g] = Datum_chrpos(datum);
	new->values[newc][g] = Datum_value(datum);
	Datum_free(&datum);
      }
      FREE(array);
      List_free(&(datalists[newc]));
    }
  }
  FREE(datalists);

  return new;
}


/* Value (e.g., log ratio) */
static double
calc_xpos (int lineindex, double linesep, double value, double valuefactor,
	   double left) {
  /* Note: Must subtract value, so that positive numbers go upward */
  return left + linesep*(lineindex+1) - value*valuefactor;
}

/* Coordinate */
static double
calc_ypos (unsigned int coord, unsigned int maxcoord, unsigned int mincoord,
	   double top, double bottom) {
  return (top-bottom)/(maxcoord-mincoord)*(coord-mincoord) + bottom;
}


static char *
coord_to_string (unsigned int coord) {
  char *string, Buffer[100];

  if (coord > 1000000) {
    sprintf(Buffer,"%.1fM",coord/(double) 1000000.0);
  } else if (coord > 1000) {
    sprintf(Buffer,"%dK",coord/1000);
  } else {
    sprintf(Buffer,"%d",coord);
  }
  string = (char *) CALLOC(strlen(Buffer)+1,sizeof(char));
  strcpy(string,Buffer);

  return string;
}


static void
make_coordinate_axis (unsigned int maxcoord, unsigned int mincoord, 
		      int nincluded, double chromsep,
		      double top, double bottom, double left) {
  double xleft, xright, y, numticks, interval;
  unsigned int coord, startcoord;
  int i;
  char *coordstring;

  printf("gsave\n");
  printf("0.5 setgray\n");

  /* Pick interval so numticks is <= 50 */
  interval = 10;
  numticks = (maxcoord - mincoord)/interval;
  while (numticks > 50) {
    interval *= 10.0;
    numticks = (maxcoord - mincoord)/interval;
  }
  startcoord = (mincoord/interval + 1)*interval;

  xleft = left - 10.0;
  xright = calc_xpos(nincluded-1,chromsep,0.0,1.0,left);
  for (coord = startcoord, i = 0; coord < maxcoord; coord += interval, i++) {
    printf("newpath\n");
    y = calc_ypos(coord,maxcoord,mincoord,top,bottom);
    printf("%.2f %.2f moveto\n",xleft,y);
    printf("-6 0 rlineto\n");
#if 0
    if (gridp == true) {
      printf("%.2f %.2f lineto\n",xright,y);
    }
#endif
    printf("1 setlinewidth\n");
    printf("stroke\n");
    printf("/Helvetica findfont 6 scalefont setfont\n");
    printf("%.2f %.2f moveto\n",xleft-10,y-3);
    coordstring = coord_to_string(coord);
    printf("(%s) stringwidth pop neg 0 rmoveto\n",coordstring);
    printf("(%s) show\n",coordstring);
    FREE(coordstring);
  }

  printf("newpath\n");
  printf("%.2f %.2f moveto\n",xleft,top);
  printf("%.2f %.2f lineto\n",xleft,bottom);
  printf("0.5 setlinewidth\n");
  printf("stroke\n");
  printf("grestore\n");
  return;
}

static void
print_chromosomes (IIT_T chromosome_iit, Chrsubset_T chrsubset,
		   unsigned int mincoord, unsigned int maxcoord, 
		   double chromsep, double top, double bottom, double left) {
  int nincluded, newc, oldc;
  unsigned int chrlength;
  double chrend, x, y;
  Interval_T interval;

  printf("0 setgray\n");
  nincluded = Chrsubset_nincluded(chrsubset,chromosome_iit);
  for (newc = 1; newc <= nincluded; newc++) {
    oldc = Chrsubset_oldindex(chrsubset,newc);
    x = left + chromsep*(double) newc;
    y = bottom - 10;
    printf("gsave\n");
    printf("%.2f %.2f moveto\n",x,y);
    printf("90 rotate\n");
    printf("(%s) stringwidth pop neg -2 rmoveto\n",
	   IIT_label(chromosome_iit,oldc));
    printf("(%s) show\n",IIT_label(chromosome_iit,oldc));
    printf("grestore\n");

    interval = IIT_interval(chromosome_iit,oldc);
    chrlength = Interval_length(interval);
    chrend = calc_ypos(chrlength,maxcoord,mincoord,top,bottom);
    if (chrend > top) {
      chrend = top;
    }

    printf("newpath\n");
    printf("%.2f %.2f moveto\n",x,bottom);
    printf("%.2f %.2f lineto\n",x,chrend);
    printf("0.5 setlinewidth\n");
    printf("stroke\n");
  }
  return;
}


static void
print_title (char *title, double top, double left, double annotheight, double annotwidth) {

  printf("/Helvetica findfont 12 scalefont setfont\n");
  printf("/yline %.2f def\n",top+annotheight);
  printf("(");
  printf("%s",title);
  printf(") %.2f\n",annotwidth);
  printf("{%.2f yline moveto show\n",left);
  printf(" /yline yline 14 sub def}\n");
  printf("BreakIntoLines\n");
  printf("/Helvetica findfont 8 scalefont setfont\n");

  return;
}


static void
start_genome_page (T this, int pagenum, char *title,
		   IIT_T chromosome_iit, Chrsubset_T chrsubset,
		   unsigned int mincoord, unsigned int maxcoord,
		   int nincluded, double chromsep, 
		   double top, double bottom, double left, double right,
		   double annotheight, double annotwidth) {

  printf("%%%%Page: %d %d\n",pagenum,pagenum);
  printf("/Helvetica findfont 8 scalefont setfont\n");
  printf("180 rotate\n");
  printf("0 -792 translate\n");
#if 0
  if (magnification != 1.0) {
    printf("%.2f %.2f scale\n",magnification,magnification);
    /* Don't understand why, but this works for 1.5 */
    /* printf("0 -275 translate\n"); */
    /* Don't understand why, but this works for 2.0 */
    printf("0 -550 translate\n");
  }
#endif
  printf("90 rotate\n");

#if 0
  header_print(logop,landscapep,paperxdim,paperydim);
#endif

  if (title != NULL) {
    print_title(title,top,left,annotheight,annotwidth);
  }

  make_coordinate_axis(maxcoord,mincoord,nincluded,chromsep,
		       top,bottom,left);

  return;
}



static void
end_genome_page (T this, int pagenum, char *title,
		   IIT_T chromosome_iit, Chrsubset_T chrsubset,
		   unsigned int mincoord, unsigned int maxcoord,
		   int nincluded, double chromsep, 
		   double top, double bottom, double left, double right,
		   double annotheight, double annotwidth) {

  print_chromosomes(chromosome_iit,chrsubset,mincoord,maxcoord,chromsep,
		    top,bottom,left);
  printf("showpage\n");
  return;
}


void
Plotdata_one_signature (T this, int s, Uintlist_T *segment_startpositions,
			Uintlist_T *segment_endpositions, Doublelist_T *segment_means,
			char *title, bool logp,
			IIT_T chromosome_iit, Chrsubset_T chrsubset,
			unsigned int mincoord, unsigned int maxcoord,
			int nincluded, double valuefactor, 
			double top, double bottom, double left, double right, 
			double annotheight, double annotwidth, 
			bool segmentp, bool skipvaluesp) {
  double red, green, blue;
  int c, g;
  unsigned int chrpos;
  double value;
  double x, x0, x1, y, y1, y2, chromsep;
  Uintlist_T p, q;
  Doublelist_T r;
  
  chromsep = (right - left)/(double) nincluded;
  start_genome_page(this,/*pagenum*/1,title,chromosome_iit,chrsubset,mincoord,maxcoord,
		    nincluded,chromsep,top,bottom,left,right,annotheight,annotwidth);

  if (segmentp == true && skipvaluesp == false) {
    red = green = blue = 0.7;
  } else {
    Color_rgb(&red,&green,&blue,1.0,s,this->nsamples);
  }

  for (c = 1; c <= this->nincluded; c++) {
    if (skipvaluesp == false) {
      printf("0.25 setlinewidth\n");
      printf("%f %f %f setrgbcolor\n",red,green,blue);
      for (g = 0; g < this->nvalues[c-1]; g++) {
	chrpos = this->chrpositions[c-1][g];
	value = this->values[c-1][g];
	if (c == 2) {
	  printf("%f\n",value);
	}
	if (logp == true) {
	  if (value <= 0.0) {
	    value = 0.0;
	  } else {
	    value = log(value);
	  }
	}

	y = calc_ypos(chrpos,maxcoord,mincoord,top,bottom);
	x0 = calc_xpos(c-1,chromsep,0.0,valuefactor,left);
	x1 = calc_xpos(c-1,chromsep,value,valuefactor,left);

	printf("newpath\n");
	printf("%.2f %.2f moveto\n",x0,y);
	printf("%.2f %.2f rlineto\n",x1-x0,y-y);
	printf("stroke\n");
      }
    }

    if (segmentp == true) {
      printf("1.0 setlinewidth\n");
      printf("%f %f %f setrgbcolor\n",0.0,0.0,1.0);
      printf("/Helvetica findfont 7 scalefont setfont\n");
      for (p = segment_startpositions[c-1], q = segment_endpositions[c-1], r = segment_means[c-1];
	   p != NULL; p = Uintlist_next(p), q = Uintlist_next(q), r = Doublelist_next(r)) {
	y1 = calc_ypos(Uintlist_head(p),maxcoord,mincoord,top,bottom);
	y2 = calc_ypos(Uintlist_head(q),maxcoord,mincoord,top,bottom);
	x = calc_xpos(c-1,chromsep,Doublelist_head(r),valuefactor,left);

	printf("newpath\n");
	printf("%.2f %.2f moveto\n",x,y1);
	printf("%.2f %.2f rlineto\n",x-x,y2-y1);
	printf("stroke\n");

	x = calc_xpos(c-1,chromsep,0.0,valuefactor,left);
	y = (y1 + y2)/2.0;
	printf("%.2f %.2f moveto\n",x+3,y-2);
	printf("(%.1f) show\n",Doublelist_head(r));
      }
    }
  }

  end_genome_page(this,/*pagenum*/1,title,chromosome_iit,chrsubset,mincoord,maxcoord,
		  nincluded,chromsep,top,bottom,left,right,annotheight,annotwidth);

  return;
}




