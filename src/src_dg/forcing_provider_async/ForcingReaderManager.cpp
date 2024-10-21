#include "ForcingReaderManager.h"
#include <stdexcept>
#include <string>
#include <map>
#include <memory>
#include <vector>
#include <thread>
#include <sstream>
#include <iostream>

extern "C" {
    void fortran_call(const int index_id);
}


using namespace std;

namespace AWI
{
   ForcingReader::ForcingReader(const int index_id_)
   : index_id(index_id_)
   {
   }


   void ForcingReader::readForcing()
   {
      int idx = index_id;
      fortran_call(idx);
   }
   
   
   void ForcingReaderManager::addForcingField(const int index_id)
   {
      string name = std::to_string(index_id);
      map<string, unique_ptr<ForcingReader> >::iterator it = forcing_readers.find(name);
      if(it != forcing_readers.end())
      {
         std::stringstream exceptionMessage;
         exceptionMessage << __FILE__ << ":" << __LINE__ <<" forcing reader already exists: "<<name;
         throw std::runtime_error(exceptionMessage.str());
      }
      
      forcing_readers[name] = unique_ptr<ForcingReader>(new ForcingReader(index_id));
   }

   
   void ForcingReaderManager::begin(const int index_id)
   {
      string name = std::to_string(index_id);
      map<string, unique_ptr<ForcingReader> >::iterator it = forcing_readers.find(name);
      if(it != forcing_readers.end())
      {
         unique_ptr<ForcingReader> &w = it->second;         
         worker_threads[name] = new thread(&ForcingReader::readForcing, w.get());
         //w.get()->readForcing();
      }
      else
      {
         std::stringstream exceptionMessage;
         exceptionMessage << __FILE__ << ":" << __LINE__ <<" unknown forcing reader: "<<name;
         throw std::runtime_error(exceptionMessage.str());
      }
   }
   
   
   void ForcingReaderManager::end(const int index_id)
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
         exceptionMessage << __FILE__ << ":" << __LINE__ <<" unknown forcing reader: "<<name;
         throw std::runtime_error(exceptionMessage.str());
      }
   }
   
}
