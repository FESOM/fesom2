/*--------------------------------------------------------------------
  parms_calloc : allocate an array of count entries of size bytes
                 each. 
  parms_free   : free the memory pointed by a pointer.
  parms_malloc : return a pointer to the allocated memory.
  parms_resize : rechange the size of the memory pointed to by a
                 pointer. 
  
  Those functions are NOT used directly. You should use macros
  PARMS_NEW, PARMS_NEWARRAY, PARMS_NEWARRAY0 instead.

  $Id: parms_mem.c,v 1.1.1.1 2006-11-27 22:28:01 zzli Exp $
--------------------------------------------------------------------*/

#include "parms_mem.h"

#ifdef C99

/** 
 * Return a pointer to the allocated memory.
 * 
 * @param size  The size of memory allocated.                 
 * @param line  The line number at which a error occurs.      
 * @param func  The function name within which a error occurs.
 * @param fname The file name in which a error occurs.
 * 
 * @return A pinter to the allocated memory.
 */
void *parms_malloc(long size, int line, const char *func, 
		   const char *fname)
{
  void *ptr;

  if (size <= 0) {
    fprintf(stderr, "Error: size = %ld at line %d in function %s at file %s\n", size, line, func, fname);
    exit(1);
  }
  ptr = malloc(size);
  if (ptr == NULL) {
    fprintf(stderr, "Error: size = %ld at line %d in function %s at file %s\n", size, line, func, fname);
    exit(1);
  }
  return ptr;
}

/** 
 * Allocate a memory for count etries each with size bytes.
 * 
 * @param count  The number of entries allocated.              
 * @param size 	 The size in bytes of each entry.              
 * @param line 	 The line number at which a error occurs.      
 * @param func 	 The function name within which a error occurs.
 * @param fname  The file name in which a error occurs.
 * 
 * @return A pointer to the allocated memory.
 */
void *parms_calloc(long count, long size, const int line, const char *func, 
		   const char *fname)
{
  void *ptr;

  if (count <= 0 || size <= 0) {
    fprintf(stderr, "Error: count = %ld, size = %ld at line %d in function %s at file %s\n", count, size, line, func, fname);
    exit(1);
  }

  ptr = calloc(count, size);
  if (ptr == NULL) {
    fprintf(stderr, "Error: count = %ld, size = %ld at line %d in function %s at file %s\n", count, size, line, func, fname);
    exit(1);
  }
  return ptr;
}

/** 
 * Rechange the size of the memory allocated.
 * 
 * @param ptr    The original pointer.                         
 * @param size 	 The new size of the memory.                   
 * @param line 	 The line number at which a error occurs.      
 * @param func 	 The function name within which a error occurs.
 * @param fname  The file name in which a error occurs.
 * 
 * @return A pointer to the memory of size bytes.
 */
void *parms_resize(void *ptr, long size, const int line, 
		   const char *func, const char *fname)
{
  if (ptr == NULL) {
    fprintf(stderr, "Error: ptr is a NULL pointer when reallocing memory at line %d in function %s at file %s\n", line, func, fname);
    exit(1);
  }
  if (size <= 0) {
    fprintf(stderr, "Error: size = %ld at line %d in function %s at file %s\n", size, line, func, fname);
    exit(1);
  }
  ptr = realloc(ptr, size);
  return ptr;
}

/** 
 * Free the memory pointed by ptr.
 * 
 * @param ptr   A pointer.                                    
 * @param line 	The line number at which a error occurs.      
 * @param func 	The function name within which a error occurs.
 * @param fname The file name in which a error occurs.
 */
void parms_free(void *ptr, const int line, const char *func, const char *fname)
{
  if (ptr == NULL) {
    fprintf(stderr, "Cannot free NULL at line %d in function %s at file %s\n", line, func, fname);
    exit(1);
  }
  free(ptr);
}

#else

void *parms_malloc(long size, int line, const char *fname)
{
  void *ptr;

  if (size <= 0) {
    fprintf(stderr, "Error: size = %ld at line %d in file %s\n", size, line, fname);
    exit(1);
  }
  ptr = malloc(size);
  if (ptr == NULL) {
    fprintf(stderr, "Error: size = %ld at line %d in file %s\n", size, line, fname);
    exit(1);
  }
  return ptr;
}

void *parms_calloc(long count, long size, const int line, const char *fname)
{
  void *ptr;

  if (count <= 0 || size <= 0) {
    fprintf(stderr, "Error: count = %ld, size = %ld at line %d in file %s\n", count, size, line, fname);
    exit(1);
  }

  ptr = calloc(count, size);
  if (ptr == NULL) {
    fprintf(stderr, "Error: count = %ld, size = %ld at line %d in file %s\n", count, size, line, fname);
    exit(1);
  }
  return ptr;
}

void *parms_resize(void *ptr, long size, const int line, 
		   const char *fname)
{
  if (ptr == NULL) {
    fprintf(stderr, "Error: ptr is a NULL pointer when reallocing memory at line %d in file %s\n", line, fname);
    exit(1);
  }
  if (size <= 0) {
    fprintf(stderr, "Error: size = %ld at line %d in file %s\n", size, line, fname);
    exit(1);
  }
  ptr = realloc(ptr, size);
  return ptr;
}

void parms_free(void *ptr, const int line, const char *fname)
{
  if (ptr == NULL) {
    fprintf(stderr, "Cannot free NULL at line %d in file %s\n", line, fname);
    exit(1);
  }
  free(ptr);
}

#endif 

