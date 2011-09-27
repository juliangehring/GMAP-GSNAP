static char rcsid[] = "$Id: datadir.c,v 1.19 2005/07/08 07:58:28 twu Exp $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "datadir.h"
#include <stdio.h>
#include <stdlib.h>		/* For getenv */
#include <string.h>
#include <strings.h>		/* For rindex */
#include <pwd.h>
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>		/* Needed for dirent.h */
#endif
#if HAVE_DIRENT_H
# include <dirent.h>
# define NAMLEN(dirent) strlen((dirent)->d_name)
#else
# define dirent direct
# define NAMLEN(dirent) (dirent)->d_namlen
# if HAVE_SYS_NDIR_H
#  include <sys/ndir.h>
# endif
# if HAVE_SYS_DIR_H
#  include <sys/dir.h>
# endif
# if HAVE_NDIR_H
#  include <ndir.h>
# endif
#endif
#include "mem.h"

/* Note: GMAPDB is defined externally by configure */
#ifndef GMAPDB
#error A default value for GMAPDB was not provided to configure.  Please do so, or edit the Makefile
#endif

static char *
read_config_file (FILE *fp, char *tag) {
  char *directory, seentag[1024], dirbuffer[1024], Buffer[1024];

  while (fgets(Buffer,1024,fp) != NULL) {
    if (sscanf(Buffer,"%s=%s",seentag,dirbuffer) > 0 && !strcmp(seentag,tag)) {
      directory = (char *) CALLOC(strlen(dirbuffer)+1,sizeof(char));
      strcpy(directory,dirbuffer);
      return directory;
    }
  }
  return NULL;
}

static FILE *
find_homedir_config () {
  FILE *fp = NULL;
  struct passwd *p;
  char *user, *configfile;

  if ((user = getenv("USER")) != NULL &&
      (p = getpwnam(user)) != NULL) {
    configfile = (char *) CALLOC(strlen(p->pw_dir)+strlen("/")+strlen(".gmaprc")+1,sizeof(char));
    sprintf(configfile,"%s/.gmaprc",p->pw_dir);
    fp = fopen(configfile,"r");
    FREE(configfile);
  }
  return fp;
}
    

static char *
find_fileroot (char *genomesubdir) {
  char *fileroot, *filename, *p;
  struct dirent *entry;
  DIR *dp;

  if ((dp = opendir(genomesubdir)) == NULL) {
    fprintf(stderr,"Unable to open directory %s\n",genomesubdir);
    exit(9);
  }
  while ((entry = readdir(dp)) != NULL) {
    filename = entry->d_name;
    if ((p = rindex(filename,'.')) != NULL) {
      if (!strcmp(p,".version")) {
	fileroot = (char *) CALLOC(p - &(filename[0]) + 1,sizeof(char));
	strncpy(fileroot,filename,p-&(filename[0]));
	if (closedir(dp) < 0) {
	  fprintf(stderr,"Unable to close directory %s\n",genomesubdir);
	}
	return fileroot;
      }
    }
  }

  if (closedir(dp) < 0) {
    fprintf(stderr,"Unable to close directory %s\n",genomesubdir);
  }
  fprintf(stderr,"Unable to find file ending with .version in directory %s\n",genomesubdir);
  exit(9);
}




static char *
get_dbversion (char *filename) {
  char *dbversion = NULL, Buffer[100], *p;
  FILE *fp;

  fp = fopen(filename,"r");
  if (!fp) {
    return NULL;
  } else if (fgets(Buffer,100,fp) == NULL) {
    fclose(fp);
    return NULL;
  } else {
    if ((p = rindex(Buffer,'\n')) != NULL) {
      *p = '\0';
    }
    fclose(fp);
  }
  
  dbversion = (char *) CALLOC(strlen(Buffer)+1,sizeof(char));
  strcpy(dbversion,Buffer);
  return dbversion;
}


/* Allocates space for genomesubdir, fileroot, and dbversion */
char *
Datadir_find_genomesubdir (char **fileroot, char **dbversion,
			   char *user_genomedir, char *dbroot) {
  FILE *fp;
  char *genomesubdir, *genomedir, *filename, *p, *dbrootdir, *newgenomedir;

  /* First get genomedir */
  if (user_genomedir != NULL) {
    genomedir = (char *) CALLOC(strlen(user_genomedir)+1,sizeof(char));
    strcpy(genomedir,user_genomedir);

  } else if (getenv("GMAPDB") != NULL) {
    /* Use genomedir provided by environment variable */
    genomedir = (char *) CALLOC(strlen(getenv("GMAPDB"))+1,sizeof(char));
    strcpy(genomedir,getenv("GMAPDB"));

  } else if ((fp = fopen("./.gmaprc","r")) != NULL) {
    genomedir = read_config_file(fp,"GMAPDB");
    fclose(fp);

  } else if ((fp = find_homedir_config()) != NULL) {
    genomedir = read_config_file(fp,"GMAPDB");
    fclose(fp);

  } else {
    genomedir = (char *) CALLOC(strlen(GMAPDB)+1,sizeof(char));
    strcpy(genomedir,GMAPDB);
  }

  /* Append directory part of dbroot */
  if ((p = rindex(dbroot,'/')) != NULL) {
    *p = '\0';
    p++;
    dbrootdir = dbroot;
    dbroot = p;
    
    newgenomedir = (char *) CALLOC(strlen(genomedir)+strlen(dbrootdir)+1,sizeof(char));
    sprintf(newgenomedir,"%s/%s",genomedir,dbrootdir);
    FREE(genomedir);
    genomedir = newgenomedir;
  }

  /* Find version file */
  filename = (char *) CALLOC(strlen(genomedir) + strlen("/") + strlen(dbroot) + 
			     strlen(".version") + 1,sizeof(char));
  sprintf(filename,"%s/%s.version",genomedir,dbroot);
  if ((fp = fopen(filename,"r")) != NULL) {
    /* Found in top-level genomedir */
    fclose(fp);
    FREE(filename);
    genomesubdir = genomedir;
    *fileroot = (char *) CALLOC(strlen(dbroot)+1,sizeof(char));
    strcpy(*fileroot,dbroot);

  } else {
    FREE(filename);
    genomesubdir = (char *) CALLOC(strlen(genomedir) + strlen("/") + strlen(dbroot) + 1,sizeof(char));
    sprintf(genomesubdir,"%s/%s",genomedir,dbroot);
    FREE(genomedir);

    if ((*fileroot = find_fileroot(genomesubdir)) != NULL) {
      /* Found in subdirectory */
    } else {
      fprintf(stderr,"Error: Can't open genome files in %s or %s.\n",genomedir,genomesubdir);
      fprintf(stderr,"       Please specify directory using -D flag, GMAPDB environment variable,\n");
      fprintf(stderr,"       or a configuration file .gmaprc with the line GMAPDB=<directory>;\n");
      fprintf(stderr,"       or have GMAP recompiled using the --with-gmapdb flag to configure.\n");
      exit(9);
    }
  }
    
  filename = (char *) CALLOC(strlen(genomesubdir) + strlen("/") + strlen(*fileroot) + strlen(".version") + 1,
			     sizeof(char));
  sprintf(filename,"%s/%s.version",genomesubdir,*fileroot);
  if ((*dbversion = get_dbversion(filename)) == NULL) {
    /* Something wrong with version file.  Use dbroot instead */
    *dbversion = (char *) CALLOC(strlen(dbroot)+1,sizeof(char));
    strcpy(*dbversion,dbroot);
  }

  FREE(filename);

  return genomesubdir;
}


char *
Datadir_find_mapdir (char *user_mapdir, char *genomesubdir, char *fileroot) {
  char *mapdir;

  if (user_mapdir != NULL) {
    mapdir = (char *) CALLOC(strlen(user_mapdir)+1,sizeof(char));
    strcpy(mapdir,user_mapdir);
  } else {
    mapdir = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+strlen(fileroot)+
			     strlen(".maps")+1,sizeof(char));
    sprintf(mapdir,"%s/%s.maps",genomesubdir,fileroot);
  }

  return mapdir;
}


