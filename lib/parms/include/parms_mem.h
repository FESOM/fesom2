/**
 * @file   parms_mem.h
 * @author Zhongze Li
 * @date   Tue Oct 17 11:48:51 2006
 * 
 * @brief  Macros and functions for memory allocation in pARMS.  
 * 
 * - PARMS_NEW and PARMS_NEW0 allocate memory for an object of struct
 *   type.  
 * - PARMS_NEWARRAY and PARMS_NEWARRYA0 allocate memory for an array.  
 * - PARMS_ALLOC(n) is used for allocating n bytes, which is useful for
 *   reducing the number of invoking malloc function. by combining
 *   multiple malloc calls together.  
 * - PARMS_RESIZE changes the size of the allocated memory. 
 * - PARMS_FREE frees the memory allocated.
 * 
 */

#ifndef _PARMS_MEM_H_
#define _PARMS_MEM_H_

#include <stdlib.h>
#include <string.h>
#include "parms_sys.h"

PARMS_CXX_BEGIN

#ifdef __STDC__
#if defined(__STDC_VERSION__) && __STDC_VERSION__ >=199901L
#define C99
#elif defined(__STDC_VERSION__) && __STDC_VERSION__>=199409L
#define C89
#else
#define C89
#endif
#endif

#ifdef C99
extern void *parms_malloc(long size, const int line, const char *func, 
    				  const char *fname);
extern void *parms_calloc(long count, long size, int line, 
			  const char *func, const char *fname);
extern void *parms_resize(void *ptr, long size, const int line, 
			  const char *func, const char *fname);
extern void parms_free(void *ptr, const int line, const char *func,
		       const char *fname);
#define PARMS_ALLOC(n) parms_malloc((long)(n), __LINE__, __func__, __FILE__)
#define PARMS_CALLOC(n, size) parms_calloc((long)(n), (long)(size), __LINE__, __func__, __FILE__)
#define PARMS_NEW(p) ((p) = PARMS_ALLOC(sizeof *(p)))
#define PARMS_NEW0(p) ((p) = PARMS_CALLOC(1, sizeof *(p)))
#define PARMS_NEWARRAY(p, n) ((p) = PARMS_ALLOC((n)*sizeof *(p)))
#define PARMS_NEWARRAY0(p, n) ((p) = PARMS_CALLOC((n), sizeof *(p)))
#define PARMS_RESIZE(p, size) ((p) = parms_resize((p), ((size) * (long)sizeof (*(p))), __LINE__, __func__, __FILE__))
#define PARMS_FREE(p)  (parms_free((p), __LINE__, __func__, __FILE__), (p) = NULL)


#else

extern void *parms_malloc(long size, const int line, const char *fname);
extern void *parms_calloc(long count, long size, int line, 
			  const char *fname);
extern void *parms_resize(void *ptr, long size, const int line, 
			  const char *fname);
extern void parms_free(void *ptr, const int line, const char *fname);

#define PARMS_ALLOC(n) parms_malloc((long)(n), __LINE__,  __FILE__)
#define PARMS_CALLOC(n, size) parms_calloc((long)(n), (long)(size), __LINE__,__FILE__)
#define PARMS_NEW(p) ((p) = PARMS_ALLOC(sizeof *(p)))
#define PARMS_NEW0(p) ((p) = PARMS_CALLOC(1, sizeof *(p)))
#define PARMS_NEWARRAY(p, n) ((p) = PARMS_ALLOC((n)*sizeof *(p)))
#define PARMS_NEWARRAY0(p, n) ((p) = PARMS_CALLOC((n), sizeof *(p)))
#define PARMS_RESIZE(p, size) ((p) = parms_resize((p), ((size) * (long)sizeof (*(p))), __LINE__, __FILE__))
#define PARMS_FREE(p)  (parms_free((p), __LINE__, __FILE__), (p) = NULL)


#endif 
#define PARMS_MEMCPY(dest, src, size) memcpy((dest), (src), (size)*(long)sizeof *(src))


PARMS_CXX_END

#endif 
