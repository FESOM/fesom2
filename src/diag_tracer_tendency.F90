module mod_diag_tracer_tendency
   ! Total tracer tendency for the ltra_diag budget diagnostics [mmol/m2/s]:
   !    tra_tendency = (c*h_new at end of step - c*h_new at start of step) / dt
   ! store_tracer_old is called at the start of the step (before REcoM),
   ! tracer_tendency after oce_timestep_ale. hnode = hnode_new at the end of
   ! every step (update_thickness_ale), so both use the thickness that matches
   ! the tracer values at that time.
   USE MOD_MESH
   USE MOD_TRACER
   USE MOD_PARTIT
   implicit none
   private
   public :: store_tracer_old, tracer_tendency

contains

subroutine store_tracer_old(tracers, partit, mesh)
   type(t_tracer), intent(inout) :: tracers
   type(t_partit), intent(in)    :: partit
   type(t_mesh)  , intent(in)    :: mesh

   integer :: n, nz, nu1, nl1, tr_num

   do tr_num=1, tracers%num_tracers
      if (.not. tracers%data(tr_num)%ltra_diag) cycle
      do n = 1, partit%myDim_nod2D
         nu1 = mesh%ulevels_nod2D(n)
         nl1 = mesh%nlevels_nod2D(n)
         do nz = nu1, nl1-1
            tracers%work%tracer_old(nz,n,tr_num) = tracers%data(tr_num)%values(nz,n)*mesh%hnode_new(nz,n)
         end do
      end do
   end do
end subroutine store_tracer_old

subroutine tracer_tendency(tracers, partit, mesh)
   use g_config, only : dt
   type(t_tracer), intent(inout) :: tracers
   type(t_partit), intent(in)    :: partit
   type(t_mesh)  , intent(in)    :: mesh

   integer :: n, nz, nu1, nl1, tr_num

   do tr_num=1, tracers%num_tracers
      if (.not. tracers%data(tr_num)%ltra_diag) cycle
      do n = 1, partit%myDim_nod2D
         nu1 = mesh%ulevels_nod2D(n)
         nl1 = mesh%nlevels_nod2D(n)
         do nz = nu1, nl1-1
            tracers%work%tra_tendency(nz,n,tr_num) = &
                 (tracers%data(tr_num)%values(nz,n)*mesh%hnode_new(nz,n) &
                - tracers%work%tracer_old(nz,n,tr_num)) / dt
         end do
      end do
   end do
end subroutine tracer_tendency

end module mod_diag_tracer_tendency
