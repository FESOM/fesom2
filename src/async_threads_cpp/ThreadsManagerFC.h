#ifndef ThreadsManager_47775EBD_FF03_4546_90E9_992F050CA628
#define ThreadsManager_47775EBD_FF03_4546_90E9_992F050CA628

#include "ThreadsManagerFCMacros.h"


extern "C"
{
   void init_ccall(const int *index);
   void begin_ccall(const int *index);
   void end_ccall(const int *index);
}

#endif
