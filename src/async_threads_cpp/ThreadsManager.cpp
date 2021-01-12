#include "ThreadsManager.h"
#include <stdexcept>
#include <string>
#include <map>
#include <memory>
#include <thread>
#include <sstream>

extern "C" {
    void async_threads_execute_fcall(const int index_id);
}


using namespace std;

namespace AWI
{
   FortranCallback::FortranCallback(const int index_id_)
   : index_id(index_id_)
   {
   }


   void FortranCallback::executeCallback()
   {
      int idx = index_id;
      async_threads_execute_fcall(idx);
   }
   
   
   void ThreadsManager::addThread(const int index_id)
   {
      string name = std::to_string(index_id); // todo: we do not seem to need a string here, use int in the map
      map<string, unique_ptr<FortranCallback> >::iterator it = callbacks.find(name);
      if(it != callbacks.end())
      {
         std::stringstream exceptionMessage;
         exceptionMessage << __FILE__ << ":" << __LINE__ <<" thread already exists: "<<name;
         throw std::runtime_error(exceptionMessage.str());
      }
      
      callbacks[name] = unique_ptr<FortranCallback>(new FortranCallback(index_id));
   }

   
   void ThreadsManager::begin(const int index_id)
   {
      string name = std::to_string(index_id);
      map<string, unique_ptr<FortranCallback> >::iterator it = callbacks.find(name);
      if(it != callbacks.end())
      {
         unique_ptr<FortranCallback> &w = it->second;         
         worker_threads[name] = new thread(&FortranCallback::executeCallback, w.get());
      }
      else
      {
         std::stringstream exceptionMessage;
         exceptionMessage << __FILE__ << ":" << __LINE__ <<" unknown thread: "<<name;
         throw std::runtime_error(exceptionMessage.str());
      }
   }
   
   
   void ThreadsManager::end(const int index_id)
   {
      string name = std::to_string(index_id);
      map<string, thread*>::iterator it = worker_threads.find(name);
      if(it != worker_threads.end())
      {
         thread *t = it->second;
         t->join();
         worker_threads.erase(it);
         delete t;
      }
      else
      {
         std::stringstream exceptionMessage;
         exceptionMessage << __FILE__ << ":" << __LINE__ <<" unknown thread: "<<name;
         throw std::runtime_error(exceptionMessage.str());
      }
   }
   
}
