/*--------------------------------------------------------------------
  parms_ViewerCreate    : create a parms_Viewer object.
  parms_ViewerFree      : free the memory for a parms_Viewer object. 
  parms_ViewerGetFP     : retrieve the file pointer. 
  parms_ViewerGetFname  : get the file name.
  parms_ViewerStoreFP   : restore the file pointer.

  A code fragment for using viewer functions:

  parms_Viewer v;
  parms_Map    map;

  // create a parms_Map object.
  parms_MapCreateFromGlobal(&map, ...);
  // create a parms_Viewer object.
  parms_ViewerCreate(&v, "foo");
  // dump map to a file via v. i.e. foo0.dat on PE0, f001.dat on PE1 
  parms_MapView(map, v);
  ...
  // free the memory for the parms_Viewer object v.
  pams_ViewerFree(&v);

  $Id: parms_viewer.c,v 1.1.1.1 2006-11-27 22:28:01 zzli Exp $
  ------------------------------------------------------------------*/
#include "./include/parms_viewer_impl.h"

/** 
 * Create a parms_Viewer object.
 *
 * If PARMS_COUT and PARMS_CERR are input as fname, they stand for
 * standard output and standard error, respectively. Otherwise, each
 * PE create a file "fnameID.dat". ID stands for ID of PE.
 *
 * @param self  A pointer to the parms_Viewer object.
 * @param fname A file name to store data.
 * 
 * @return 0 on success.
 */
int parms_ViewerCreate(parms_Viewer *self, char *filename)
{
  parms_Viewer new_viewer;
  int pid, length;

  PARMS_NEW0((new_viewer));
  new_viewer->ref = 1;

  length = strlen(filename) + 6;
  MPI_Comm_rank(MPI_COMM_WORLD, &pid);
  length += pid / 10 + 1;
  PARMS_NEWARRAY(new_viewer->fname, length);
  new_viewer->isstd = false;
  if (!strcmp(filename, PARMS_COUT)) {
    new_viewer->isstd = true;
    new_viewer->fp = stdout;
  }
  else if (!strcmp(filename, PARMS_CERR)) {
    new_viewer->isstd = true;
    new_viewer->fp = stderr;
  }
  else {
    strcpy(new_viewer->fname, filename);
    sprintf(new_viewer->fname, "%s_%d.dat", new_viewer->fname, pid);
  }

  *self = new_viewer;
  if (filename == NULL) {
    return -1;
  }
  else {
    return 0;
  }
}

/** 
 * Free the memory of the parms_Viewer object.
 * 
 * @param self A pointer to the parms_Viewer object.
 * 
 * @return 0 on success.
 */
int parms_ViewerFree(parms_Viewer *self)
{
  assert(self && *self);
  (*self)->ref--;
  if ((*self)->ref == 0 ) {
    if ((*self)->fname != NULL) {
      if (strcmp((*self)->fname, PARMS_COUT) && strcmp((*self)->fname, 
						       PARMS_CERR)) { 
	PARMS_FREE((*self)->fname);
      }
    }
    PARMS_FREE(*self);
  }
  return 0;
}

/** 
 * Get a pointer to file pointer. 
 * 
 * @param self A parms_Viewer object.           
 * @param fp   A pointer to the file pointer.   
 * 
 * @return 0 on success.
 */
int parms_ViewerGetFP(parms_Viewer self, FILE **fp)
{
  if (!self->isstd) {
    self->fp = fopen(self->fname, "a+");
    if (self->fp == NULL) {
      fprintf(stderr, "cannot open file %s\n", self->fname);
    }
  }

  *fp = self->fp;
  return 0;
}

/** 
 * Store fp to the parms_Viewer object.
 * 
 * @param self A parms_Viewer object. 
 * @param fp   A file pointer.     
 * 
 * @return 0 on success.
 */
int parms_ViewerStoreFP(parms_Viewer self, FILE *fp)
{
  if (!self->isstd) {
    fclose(fp);
  }
  return 0;
}

/** 
 * Retrieve the file name.
 * 
 * @param self  A parms_Viewer object.  
 * @param fname The file name retrieved.
 * 
 * @return 0 on success.
 */
int parms_ViewerGetFname(parms_Viewer self, char **fname)
{
  *fname = self->fname;
  return 0;
}


