#include "Hecuba_interface_wrapper.h"
#include "Hecuba_interface.h"

static Hecuba::Hecuba_interface hecubaManager;

void init_datamodel_c(int index)
{
   hecubaManager.init_datamodel(index);
}

void hecuba_output_c(int ts, char* variable, char* description, char* units, int freq, char* freq_unit, int accuracy )
{
   hecubaManager.save_data(ts, variable, description, units, freq, freq_unit, accuracy );
}


