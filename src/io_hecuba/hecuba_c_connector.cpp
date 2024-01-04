#include <iostream>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include <StorageNumpy.h>
#include <StorageDict.h>
#include <KeyClass.h>
#include <ValueClass.h>

#ifdef __cplusplus
extern "C" {
#endif

using StrKeyClass = KeyClass<std::string>;
using MyValueClass = ValueClass<StorageNumpy>;
//class metricsDict:public StorageDict <MyKeyClass, MyValueClass, metricsDict>{};
class metricsDict:public StorageDict <StrKeyClass, MyValueClass, metricsDict>{};

using StrValueClass = ValueClass<std::string>;// to store metadata
class MetaDictClass: public StorageDict <StrKeyClass,StrValueClass, MetaDictClass> { };

// static shared across a session	
static MetaDictClass md;
static metricsDict metrics;
static bool datamodel_loaded = false;

void initHecubaSession(const char * expname) {
    md.make_persistent(expname);
    std::string exp_vardict = "";
    exp_vardict = exp_vardict + expname + "_variables" ;//"aa"+"_variables";
    metrics.make_persistent(exp_vardict); 	//creates the tables in Cassandra. Generates a python file with the class sp
    datamodel_loaded = true;
} 

void setExpMetaData(const char * expname, const char * key, const char * value_str) {
    std::cout<< "+ setting metadata"<<std::endl;
    StrKeyClass ky = StrKeyClass(key);
    StrValueClass kv = StrValueClass(value_str);
    md[ky]=kv;
    std::cout<< "+ set "<<key<<": "<<value_str<<std::endl;
   }	

void sendMetricsToHecuba(const char * expname, const char * varname, int32_t timestep_i, int32_t chunk_id, void *data, uint32_t metadata) {
    if (not datamodel_loaded) {
      std::string expid = "";
      expid += expname;
      expid +="_variables";
      metrics.getByAlias(expid);
      datamodel_loaded = true;
    }
    StorageNumpy metrics_data(data,{metadata,1}); // instantiates a StorageNumpy with the specified info


    std::string generatedName = "";
    generatedName += expname;
    generatedName += "/";
    generatedName += varname;
    // Add the timestep_i to generatedName only if it is not -1
    if (timestep_i != -1) {
        generatedName += "/";
        generatedName += std::to_string(timestep_i);
    }

    // Add the chunk_id to generatedName only if it is not -1
    if (chunk_id != -1) {
        generatedName += "/";
        generatedName += std::to_string(chunk_id);
    }

    metrics_data.make_persistent(generatedName);  // in future releases this make_persistent will be implicit as part
					// of the assignment to a persistent object (the metricsDict)
    StrKeyClass k = StrKeyClass("/"+generatedName);
    MyValueClass v(metrics_data);

    metrics[k]=v; 				//store asynchronously and sends data through the stream
}

#ifdef __cplusplus
}
#endif
