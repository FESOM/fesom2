from hecuba import StorageObj

class myclass (StorageObj):
   '''
   @ClassField sim_id text
   @ClassField sim_info model_complex.info
   @ClassField submetrics model_complex.metrics
   '''
from hecuba import StorageObj

class info (StorageObj):
   '''
   @ClassField total_ts int
   @ClassField output_freq int
   '''
from hecuba import StorageDict

class metrics (StorageDict):
   '''
   @TypeSpec dict <<lat:float,ts:int>,mvalues:numpy.ndarray>
   '''