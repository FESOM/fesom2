#ifndef ForcingReaderManager_C6DCC9D9_E462_4BAF_8C3B_00B843F08A75
#define ForcingReaderManager_C6DCC9D9_E462_4BAF_8C3B_00B843F08A75

#include <string>
#include <map>
#include <memory>
#include <vector>
#include <thread>


namespace AWI
{
   class ForcingReader
   {
   public:
      ForcingReader(const int index_id);
      void readForcing();
   
   private:
      const int index_id;
   };

   
   class ForcingReaderManager
   {
      public:
      void addForcingField(const int index_id);
      void begin(const int index_id);
      void end(const int index_id);
      
   private:
     std::map< std::string, std::unique_ptr<ForcingReader> > forcing_readers;
     std::map<std::string, std::thread*> worker_threads;
   };
   
}

#endif
