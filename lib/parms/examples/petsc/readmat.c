#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define BUFLEN 100

int read_matrix(char mname[BUFLEN])
{
  FILE *fmat;
  char buf[BUFLEN], *tmpstring;

/*------- Read a matrix from list of Matrices -------------*/
#if defined(DBL_CMPLX)
  if (NULL == (fmat = fopen("./matfileCmplx", "r"))) {
    fprintf(stderr, "cannot open file: matfile\n");
    exit(1);
  }
#else  
    if (NULL == (fmat = fopen("./matfileReal", "r"))) {
    fprintf(stderr, "cannot open file: matfile\n");
    exit(1);
  }
#endif  
    if(fgets(buf,BUFLEN,fmat) == NULL){
      fprintf(stderr, "Error reading matrix list");
      exit(1);
    }
    tmpstring = (char *) strtok(buf, " ");
    tmpstring++;  
    tmpstring[strlen(tmpstring)-1]='\0';
    strcpy(mname, --tmpstring);

    fclose(fmat);
 
  return 0;
}
  
