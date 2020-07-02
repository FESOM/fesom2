#ifndef ThreadsManager_C6DCC9D9_E462_4BAF_8C3B_00B843F08A75
#define ThreadsManager_C6DCC9D9_E462_4BAF_8C3B_00B843F08A75

#include <string>
#include <map>
#include <memory>
#include <thread>


namespace AWI
{
   class FortranCallback
   {
   public:
      FortranCallback(const int index_id);
      void executeCallback();
   
   private:
      const int index_id;
   };

   
   class ThreadsManager
   {
      public:
      void addThread(const int index_id);
      void begin(const int index_id);
      void end(const int index_id);
      
   private:
     std::map< std::string, std::unique_ptr<FortranCallback> > callbacks;
     std::map<std::string, std::thread*> worker_threads;
   };
   
}

#endif
