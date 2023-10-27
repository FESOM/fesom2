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

void initHecubaSession(const char * expname) {
    md.make_persistent(expname);
    std::string exp_vardict = "";
    exp_vardict = exp_vardict + expname + "_variables" ;//"aa"+"_variables";
    metrics.make_persistent(exp_vardict); 	//creates the tables in Cassandra. Generates a python file with the class sp
} 

void setExpMetaData(std::string expname, std::string key, std::string value) {
    std::cout<< "+ here in set metadata"<<std::endl;
    StrKeyClass ky = StrKeyClass(key);
    StrValueClass kv = StrValueClass(value);
    md[ky]=kv;
    std::cout<< "+ end here in set metadata"<<std::endl;
   }	

void sendMetricsToHecuba(const char * expname, const char * varname, int32_t timestep_i, int32_t chunk_id, void *data, uint32_t metadata) {
    StorageNumpy metrics_data(data,{metadata,1}); // instantiates a StorageNumpy with the specified info

    std::string generatedName = "";
    generatedName = generatedName + expname + "/" + varname + "/" +std::to_string(timestep_i) + "/" +std::to_string(chunk_id); //random_generated_name();
    metrics_data.make_persistent(generatedName);  // in future releases this make_persistent will be implicit as part
					// of the assignment to a persistent object (the metricsDict)
    StrKeyClass k = StrKeyClass("/"+generatedName);
    MyValueClass v(metrics_data);

    metrics[k]=v; 				//store asynchronously and sends data through the stream
}

#ifdef __cplusplus
}
#endif
