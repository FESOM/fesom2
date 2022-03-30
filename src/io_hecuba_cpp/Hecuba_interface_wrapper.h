// this is created by the CMAKE
#include "HecubaWrapperFortranMacros.h"

extern "C"
{
   void init_datamodel_c(int index);
   void hecuba_output_c(int ts, char* variable, char* description, char* units, int freq, char* freq_unit, int accuracy );
}

