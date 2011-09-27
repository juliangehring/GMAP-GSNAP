/* $Id: datadir.h,v 1.8 2005/02/15 01:53:23 twu Exp $ */
#ifndef DATADIR_INCLUDED
#define DATADIR_INCLUDED

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

extern char *
Datadir_find_genomesubdir (char **fileroot, char **dbversion,
			   char *user_genomedir, char *dbroot);
extern char *
Datadir_find_mapdir (char *user_mapdir, char *genomesubdir, char *fileroot);

#endif

