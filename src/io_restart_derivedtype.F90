module restart_derivedtype_module
    interface
        subroutine write_all_bin_restarts(ctarr, path_in, pathi_in, ice, dynamics, tracers, partit, mesh)
            use MOD_ICE
            use MOD_DYN
            use MOD_TRACER
            use MOD_PARTIT
            use MOD_MESH
            integer, dimension(3) , intent(in)              :: ctarr
            character(len=*), intent(in)                    :: path_in
            character(len=*), intent(in)                    :: pathi_in
            type(t_ice)   , intent(inout), target, optional :: ice
            type(t_dyn)   , intent(inout), target, optional :: dynamics
            type(t_tracer), intent(inout), target, optional :: tracers
            type(t_partit), intent(inout), target, optional :: partit
            type(t_mesh)  , intent(inout), target, optional :: mesh
        end subroutine
        
        subroutine read_all_bin_restarts(path_in, ice, dynamics, tracers, partit, mesh)
            use MOD_ICE
            use MOD_DYN
            use MOD_TRACER
            use MOD_PARTIT
            use MOD_MESH
            character(len=*), intent(in)                    :: path_in
            type(t_ice)   , intent(inout), target, optional :: ice
            type(t_dyn)   , intent(inout), target, optional :: dynamics
            type(t_tracer), intent(inout), target, optional :: tracers
            type(t_partit), intent(inout), target, optional :: partit
            type(t_mesh)  , intent(inout), target, optional :: mesh
        end subroutine
    end interface
end module    
!
!
!_______________________________________________________________________________
subroutine write_all_bin_restarts(ctarr, path_in, pathi_in, ice, dynamics, tracers, partit, mesh)
end subroutine
!
!
!_______________________________________________________________________________
! read derived type binary restart files, depending on input (see optional) not 
! all derived type binaries are read --> functionalitiy for dwarfs !
subroutine read_all_bin_restarts(path_in, ice, dynamics, tracers, partit, mesh)
end subroutine
  
