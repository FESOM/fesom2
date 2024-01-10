#ifndef HECUBA_C_CONNECTOR_H 
#define HECUBA_C_CONNECTOR_H 

#ifdef __cplusplus
extern "C" {
#endif
 
void initHecubaSession(const char * expname);
void sendMetricsToHecuba(const char * expname, const char * varname, int32_t timestep, int32_t chunk_id, void *data, int32_t metadata);
void getPruneState(const char * expname, char * state);

#ifdef __cplusplus
}
#endif


#endif

