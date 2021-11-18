#include <stdlib.h>
#include <time.h>

/* Missing sys_timer for shared libraries */
double sys_timer() {
  clock_t t;
  t = clock();
  return ((double)t) / CLOCKS_PER_SEC;
}
