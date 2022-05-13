#ifndef HECUBA_C_CONNECTOR_H 
#define HECUBA_C_CONNECTOR_H 

#ifdef __cplusplus
extern "C" {
#endif
 
void start_hecuba_session();
void load_datamodel();
void hecuba_put_array_val_C(char *key, void *valueC, long unsigned int *arr_size);
void hecuba_put_array_val_C2(char *varname, double ctime, int chunk, void *valueC, long unsigned int *arr_size);

#ifdef __cplusplus
}
#endif


#endif

