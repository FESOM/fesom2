#include "parms_sys.h"
#include "DDPQ/protos.h"
#if defined(C99)
#include <tgmath.h>
#else
#include <math.h>
#endif 

int qsplitC(FLOAT *a, int *ind, int n, int ncut)
{
  /*----------------------------------------------------------------------
    |     does a quick-sort split of a real array.
    |     on input a[0 : (n-1)] is a real array
    |     on output is permuted such that its elements satisfy:
    |
    |     abs(a[i]) >= abs(a[ncut-1]) for i < ncut-1 and
    |     abs(a[i]) <= abs(a[ncut-1]) for i > ncut-1
    |
    |     ind[0 : (n-1)] is an integer array permuted in the same way as a.
    |---------------------------------------------------------------------*/
  REAL abskey;
  FLOAT  tmp;
  int j, itmp, first, mid, last;
  first = 0;
  last = n-1;
  if (ncut<first || ncut > last) return 0;
  /* outer loop -- while mid != ncut */
 label1:
  mid = first;
  abskey = ABS_VALUE(a[mid]);
  for (j=first+1; j<=last; j++) {
    if (ABS_VALUE(a[j]) > abskey) {
      tmp = a[++mid];
      itmp = ind[mid];
      a[mid] = a[j];
      ind[mid] = ind[j];
      a[j]  = tmp;
      ind[j] = itmp;
    }
  }
  /* interchange */
  tmp = a[mid];
  a[mid] = a[first];
  a[first]  = tmp;
  itmp = ind[mid];
  ind[mid] = ind[first];
  ind[first] = itmp;
  /* test for while loop */
  if (mid == ncut) return 0;
  if (mid > ncut) 
    last = mid-1;
  else
    first = mid+1;
  goto label1;
}
