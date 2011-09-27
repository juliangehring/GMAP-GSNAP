/* $Id: datadir.h,v 1.11 2005/10/14 19:03:16 twu Exp $ */
#ifndef DATADIR_INCLUDED
#define DATADIR_INCLUDED
#include <stdio.h>

extern char *
Datadir_find_genomesubdir (char **fileroot, char **dbversion,
			   char *user_genomedir, char *dbroot);
extern char *
Datadir_find_mapdir (char *user_mapdir, char *genomesubdir, char *fileroot);

extern void
Datadir_list_directory (FILE *fp, char *directory);

#endif

