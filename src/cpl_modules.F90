module cpl_config
   implicit none
   save

  ! Presetting of coupling parameters, should later be delivered by ECHAM xml - file

  ! *** enable / disable OASIS4 coupling
  real(kind=8)                  :: CplLonMin=-180.
  real(kind=8)                  :: CplLonMax=180.
  real(kind=8)                  :: CplLatMin=-90.
  real(kind=8)                  :: CplLatMax=90. 
  integer                       :: CplIDim=960!192   !320
  integer                       :: CplJDim=480!96    !160
  namelist /cpl_params/ CplLonMin,CplLonMax,CplLatMin,CplLatMax,CplIDim,CplJDim

  real(kind=8), dimension(:,:),   allocatable   :: a2o_fcorr_stat  !flux correction statistics for the output
 !real(kind=8)                                  :: time_send(2), time_recv(2)
  integer                                       :: o2a_call_count=0 
  real(kind=8), allocatable, dimension(:,:)     :: cplsnd
end module cpl_config
