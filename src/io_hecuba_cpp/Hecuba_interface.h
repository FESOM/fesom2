#include <string>
#include <hecuba/HecubaSession.h>


namespace Hecuba
{
   class Hecuba_interface
   {
     public:
        void save_data(int ts, char* variable, char* description, char* units, int freq, char* freq_unit, int accuracy );
        void init_datamodel(int index_id);
        HecubaSession session;
     private:
        int connection_status;
        
   };
}