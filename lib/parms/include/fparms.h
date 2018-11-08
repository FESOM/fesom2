#ifdef VOID_POINTER_SIZE_4
#define parms_Map integer*4
#define parms_Vec integer*4
#define parms_Mat integer*4
#define parms_PC integer*4
#define parms_Solver integer*4
#define parms_Viewer integer*4
#define parms_Timer integer*4
//#elif defined(VOID_POINTER_SIZE_8)
#else
#define parms_Map integer*8
#define parms_Vec integer*8
#define parms_Mat integer*8
#define parms_PC integer*8
#define parms_Solver integer*8
#define parms_Viewer integer*8
#define parms_Timer integer*8
#endif
	integer :: NONINTERLACED, INTERLACED
	parameter(NONINTERLACED=1, INTERLACED=0)
	integer :: INSERT, ADD
	parameter(INSERT=0, ADD=1)
	integer :: SAME_NONZERO_STRUCTURE, DIFFERENT_NONZERO_STRUCTURE
	parameter(SAME_NONZERO_STRUCTURE=0,DIFFERENT_NONZERO_STRUCTURE=1)	
        integer :: PCBJ, PCSCHUR, PCRAS, PCSCHURRAS
	integer :: PCILU0, PCILUK, PCILUT, PCARMS
	parameter(PCBJ=0, PCSCHUR=1, PCRAS=2, PCSCHURRAS=3)
	parameter(PCILU0=0, PCILUK=1, PCILUT=2, PCARMS=3)
	  integer :: SOLFGMRES, SOLGMRES, SOLBICGS, SOLCG, SOLPBICGS, SOLPBICGS_RAS, SOLBICGS_RAS
	integer :: MAXITS, KSIZE, DTOL, NEIG
	  parameter(SOLFGMRES=0,SOLGMRES=1, SOLBICGS=2, SOLCG=3, SOLPBICGS=4, SOLPBICGS_RAS=5, SOLBICGS_RAS=6)
	parameter(MAXITS=0,KSIZE=1,DTOL=2,NEIG=3)

