#include "ThreadsManagerFC.h"
#include "ThreadsManager.h"

static AWI::ThreadsManager threadsManager;


void init_ccall(const int *index)
{
   threadsManager.addThread(*index);
}


void begin_ccall(const int *index)
{
   threadsManager.begin(*index);
}


void end_ccall(const int *index)
{
   threadsManager.end(*index);
}
