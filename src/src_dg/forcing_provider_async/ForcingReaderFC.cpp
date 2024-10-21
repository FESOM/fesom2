#include "ForcingReaderFC.h"
#include "ForcingReaderManager.h"

static AWI::ForcingReaderManager forcingReaderManager;


void init_ccall(const int *index)
{
   forcingReaderManager.addForcingField(*index);
}


void begin_ccall(const int *index)
{
   forcingReaderManager.begin(*index);
}


void end_ccall(const int *index)
{
   forcingReaderManager.end(*index);
}
