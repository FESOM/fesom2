#include <hecuba/HecubaSession.h>
#include <hecuba/IStorage.h>
#include <hecuba/UUID.h>
#include <iostream>
#include "hecuba_c_connector.h"

#ifdef __cplusplus
extern "C" {
#endif

std::string UUID2str(uint64_t* c_uuid) {
    /* This MUST match with the 'cass_statement_bind_uuid' result */
    char str[37] = {};
    unsigned char* uuid = reinterpret_cast<unsigned char*>(c_uuid);
    //std::cout<< "HecubaSession: uuid2str: BEGIN "<<std::hex<<c_uuid[0]<<c_uuid[1]<<std::endl;
    sprintf(str,
        "%02x%02x%02x%02x-%02x%02x-%02x%02x-%02x%02x-%02x%02x%02x%02x%02x%02x",
        uuid[0], uuid[1], uuid[2], uuid[3],
        uuid[4], uuid[5],
        uuid[6], uuid[7],
        uuid[8], uuid[9],
        uuid[10], uuid[11], uuid[12], uuid[13], uuid[14], uuid[15]
        );
    //std::cout<< "HecubaSession: uuid2str: "<<str<<std::endl;
    return std::string(str);
}

// static shared across a session	
static HecubaSession *hsession = NULL;

void start_hecuba_session() {
    std::cout<< "+ STARTING C++ APP"<<std::endl;
    hsession = new HecubaSession();
    std::cout<< "+ Session started"<<std::endl;
}


void load_datamodel(){
    
    if (hsession == NULL) {
        start_hecuba_session();
    } else {

        std::cout<< "Session exists"<<std::endl;
    }
    
    //(*hsession).loadDataModel("model_class.yaml","model_class.py");

    // or 
    hsession->loadDataModel("model_class.yaml","model_class.py");
    std::cout<< "+ Data Model loaded"<<std::endl;

}


char * generateKey(float ctime, int chunk) {

    char * key = (char*) malloc (sizeof(float) + sizeof(int));
    float *time_key = (float*) key;
    *time_key = ctime;
    int *chunk_key = (int*) (key + sizeof(float));
    *chunk_key = chunk;
    std::cout << " generatekey sizeof(float) "<< sizeof(float) << " sizeof(int) " << sizeof(int)<< std::endl;
    return key;
}

//generate 1d array useful for scalars, not sure what is the point of shape
//in the end it is a blob of size, unless it is performance issue use 1d
//interestingly can use f-ordering by sending cols, then rows 
char * generateMetas(size_t arrsize, int dtypesize) {

    unsigned int * metas = (unsigned int *) malloc(sizeof(unsigned int) * 3);
    //using 2 because fesom scalars are 1d..
    metas[0]=2; // number of dimmensions
    metas[1]=4; //=(unsigned int) arrsize/dtypesize; //dtypesize; // number of elements in the first dimmension
    metas[2]=4;

    return (char *) metas;
}
char * generateMetas2(int array_len) {

    unsigned int * metas = (unsigned int *) malloc(sizeof(unsigned int) * 3);
    //using 2 because fesom scalars are 1d..
    metas[0]=2; // number of dimmensions
    
    // this makes 1d data into multi  cols when casted as list, like [[data[0]], [data[1]], ...]
    //metas[1]=array_len; //=(unsigned int) arrsize/dtypesize; //dtypesize; // number of elements in the first dimmension
    //metas[2]=1;
    // this makes 1d data into rows like 
    metas[1]=1;
    metas[2]=array_len; //=(unsigned int) arrsize/dtypesize; //dtypesize; // number of elements in the first dimmension

    return (char *) metas;
}

//add dtype size
void hecuba_put_array_val_C(char *key, void *valueC, long unsigned int *arr_size){
   char * numpymeta;
   load_datamodel();
   numpymeta = generateMetas(*arr_size, sizeof(float));

   IStorage *mi_sn = (*hsession).createObject("hecuba.hnumpy.StorageNumpy",key, numpymeta,valueC);
   (*mi_sn).sync();
}

//add dtype size
void hecuba_put_array_val_C2(char *varname, double ctime, int chunk, void *valueC, long unsigned int *array_size){
   char * numpymeta;
   char * key;
   char * random_key;

   load_datamodel();
   numpymeta = generateMetas2(*array_size);
   
   //if we use same key in another instance/time then the table will be replaced so
   //generate a random name,
   //optionally can just add ctime, chunk info to varname and use it
   
   //std::string random_str = (hsession->UUID2str(hsession->generateUUID()));
   // if we use above UUID then it complains error as no metas found, some how
   // (i am guessing if there is checking for  -, because using aaa-bbb as key wont work,
   // also . wont work)
   // it knows it is a UUID string so making substring using last 12 chars 
   std::string random_str = varname;
   //use substring
   random_str.append("_");
   //random_str.append((hsession->UUID2str(hsession->generateUUID())).substr(24));
   //random_str.append((UUID::UUID2str(UUID::generateUUID())).substr(24));
   random_str.append(UUID2str(UUID::generateUUID()).substr(24));
   //cast as char array
   random_key= &random_str[0];
   
   std::cout<< "+ random_key: "<<random_key<<std::endl;

   IStorage *mi_sn = (*hsession).createObject("hecuba.hnumpy.StorageNumpy",random_key, numpymeta,valueC);
    
   (*mi_sn).sync(); //delay sync 
   std::cout<< "+ done_sync: "<<random_key<<std::endl;
   // usually doesn't replace better read as get by alias
   IStorage *var = (*hsession).createObject("midict", varname);
   key = generateKey(ctime, chunk);
   var->setItem((void*)key, mi_sn);
   var->sync();
}

#ifdef __cplusplus
}
#endif

//int main(){
//  start_hecuba_session();
//  load_datamodel();
//}