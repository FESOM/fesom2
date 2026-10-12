!> @author
!> Sergei Danilov
!> @brief
!> mod_mesh_utils.F90
!! Mesh setup routines used by the domain partitioner. These are stripped down
!! versions of the routines in the main setup -- do not mix them up.

module mod_mesh_utils
    USE MOD_MESH
    USE o_PARAM
    USE g_CONFIG
    USE g_rotate_grid
    USE MOD_PARTIT
    USE omp_lib
    USE iso_c_binding, only: idx_t=>C_INT32_T

    use mod_gen_comm, only: communication_nodn, communication_elemn
    use oce_local_module, only: save_dist_mesh

    implicit none

    private
    public :: read_mesh_ini, read_mesh_cavity, test_tri_ini, &
              find_edges_ini, find_levels, find_levels_cavity, &
              edge_center, elem_center, find_elem_neighbors_ini, &
              stiff_mat_ini, communication_ini, set_par_support_ini, &
              check_partitioning

contains

!=============================================================================
!> @brief
!> Reads mesh files 
subroutine read_mesh_ini(mesh)
!
IMPLICIT NONE
!
type(t_mesh), intent(inout), target :: mesh
INTEGER                             :: nq
INTEGER                             :: n1,n2,n3
INTEGER                             :: n, nz, exit_flag
REAL(kind=WP)                       :: x1, x2, gx1, gx2
INTEGER	                            :: tag
INTEGER, allocatable                :: elem_data(:)
INTEGER                             :: i_error
! ===================
! Surface mesh
! ===================
 print *, achar(27)//'[1m'  //'____________________________________________________________'//achar(27)//'[0m'
 print *, achar(27)//'[7;1m' //' -->: read elem2d.out & nod2d.out                           '//achar(27)//'[0m'
   
  open (20,file=trim(meshpath)//'nod2d.out', status='old')
  open (21,file=trim(meshpath)//'elem2d.out', status='old')
  READ(20,*) mesh%nod2D
  ALLOCATE(mesh%coord_nod2D(2,mesh%nod2D)) 
    
  do n=1, mesh%nod2D
     read(20,*) nq, x1, x2, tag
     x1=x1*rad
     x2=x2*rad
     if (force_rotation) then
        gx1=x1
        gx2=x2
        call g2r(gx1, gx2, x1, x2)
     end if      
     mesh%coord_nod2D(1,n)=x1
     mesh%coord_nod2D(2,n)=x2
  end do
  CLOSE(20)
  READ(21,*)  mesh%elem2D    
  ALLOCATE(mesh%elem2D_nodes(4,mesh%elem2D))
  ALLOCATE(elem_data(4*mesh%elem2D))
  elem_data(:)=-1
  
  ! meshes with quads have 4 columns, but TsunAWI grids may be
  ! purely triangular, with 3 columns each. Test, how many
  ! columns there are!  
  read(21,*,iostat=i_error) elem_data(1:4*mesh%elem2D)
  ! NR Do NOT use reshape() here: ifort places the array temporary it
  !    creates on the stack, which overflows for large meshes (ng5 etc.)
  !    and causes a segfault in read_mesh_ini. Copy element by element
  !    instead, which needs no temporary and works for any mesh size.
  if (i_error == 0) then      ! There is a fourth column => quad or mixed mesh (not working yet!)
     do n=1, mesh%elem2D
        mesh%elem2D_nodes(1,n) = elem_data(4*n-3)
        mesh%elem2D_nodes(2,n) = elem_data(4*n-2)
        mesh%elem2D_nodes(3,n) = elem_data(4*n-1)
        mesh%elem2D_nodes(4,n) = elem_data(4*n)
     end do
  else     ! No fourth column => triangles only
     do n=1, mesh%elem2D
        mesh%elem2D_nodes(1,n) = elem_data(3*n-2)
        mesh%elem2D_nodes(2,n) = elem_data(3*n-1)
        mesh%elem2D_nodes(3,n) = elem_data(3*n)
        mesh%elem2D_nodes(4,n) = elem_data(3*n-2)  ! 4th node = 1st node for triangles
     end do
  end if
    
  deallocate(elem_data)
  CLOSE(21)
   	   	 
write(*,*) '========================='
write(*,*) 'Mesh is read'
write(*,*) '========================='
END SUBROUTINE read_mesh_ini
!=============================================================================
!> @brief
!> Reads mesh files 
subroutine read_mesh_cavity(mesh)
    implicit none

    type(t_mesh), intent(inout), target :: mesh
    integer                             :: node, auxi
    character(len=MAX_PATH)             :: fname
    logical                             :: file_exist=.False.
    
    !___________________________________________________________________________    
    print *, achar(27)//'[1m'  //'____________________________________________________________'//achar(27)//'[0m'
    print *, achar(27)//'[7;1m' //' -->: read cavity depth                                     '//achar(27)//'[0m'
  
    !___________________________________________________________________________
    ! read depth of cavity-ocean boundary
    if (use_cavityonelem) then
        fname = trim(meshpath)//'cavity_depth@elem.out'
    else
        fname = trim(meshpath)//'cavity_depth@node.out'
    end if     
    file_exist=.False.
    inquire(file=trim(fname),exist=file_exist) 
    if (file_exist) then
        open (21,file=fname, status='old')
        if (use_cavityonelem) then
            allocate(mesh%cavity_depth(mesh%elem2d))
        else
            allocate(mesh%cavity_depth(mesh%nod2D))
        end if 
    else
        write(*,*) '____________________________________________________________________'
        write(*,*) ' ERROR: could not find cavity file:', fname
        write(*,*) '        --> stop partitioning here !'
        write(*,*) '____________________________________________________________________'    
        stop 
    end if
    
    !___________________________________________________________________________
    auxi=mesh%nod2D
    if (use_cavityonelem) auxi=mesh%elem2d
    do node=1, auxi
        read(21,*) mesh%cavity_depth(node)
    end do
    
    !___________________________________________________________________________
    close(21)
    
end subroutine read_mesh_cavity

!=======================================================================
!> @brief 
!> Check the order of nodes in triangles; correct it if necessary to make
!! it same sense (clockwise) 
SUBROUTINE test_tri_ini(mesh)
IMPLICIT NONE

real(kind=WP)                       ::  a(2), b(2), c(2),  r
integer                             ::  n, nx, elnodes(3)
type(t_mesh), intent(inout), target :: mesh
   
   DO n=1, mesh%elem2D
      elnodes=mesh%elem2D_nodes(1:3,n)
	  
          a=mesh%coord_nod2D(:,elnodes(1))
	  b=mesh%coord_nod2D(:,elnodes(2))-a
	  c=mesh%coord_nod2D(:,elnodes(3))-a
          
	  if(b(1)>cyclic_length/2.) b(1)=b(1)-cyclic_length
          if(b(1)<-cyclic_length/2.) b(1)=b(1)+cyclic_length
	  if(c(1)>cyclic_length/2.) c(1)=c(1)-cyclic_length
          if(c(1)<-cyclic_length/2.) c(1)=c(1)+cyclic_length
	  
	    
	  r=b(1)*c(2)-b(2)*c(1)
	  if (r>0) then
	  ! Vector b is to right of c
	  ! Exchange second and third nodes:
	  
	  nx=elnodes(2)
	  elnodes(2)=elnodes(3)
	  elnodes(3)=nx
	  mesh%elem2D_nodes(1:3,n)=elnodes
      end if
   END DO
END SUBROUTINE  test_tri_ini
!=========================================================================
!> @brief
!> Finds edges. Creates 3 files: edgenum.out, edges.out, edge_tri.out
SUBROUTINE find_edges_ini(mesh)
IMPLICIT NONE


integer, allocatable                  :: aux1(:), ne_num(:), ne_pos(:,:), nn_num(:), nn_pos(:,:)
integer                               :: counter, counter_in, n, k, q
integer                               :: elem, elem1, elems(2), q1, q2
integer                               :: elnodes(4), ed(2), flag, eledges(4)
integer                               :: temp(100), node
real(kind=WP)                         :: xc(2), xe(2), ax(3), amin
logical                               :: file_exist1, file_exist2, file_exist3
type(t_mesh), intent(inout), target   :: mesh
! ====================
! (a) find edges. To make the procedure fast 
! one needs neighbourhood arrays
! ====================
print *, achar(27)//'[1m'  //'____________________________________________________________'//achar(27)//'[0m'
print *, achar(27)//'[7;1m' //' -->: compute edge connectivity                             '//achar(27)//'[0m'

allocate(ne_num(mesh%nod2d), ne_pos(MAX_ADJACENT, mesh%nod2D), nn_num(mesh%nod2D))
ne_num=0
DO n=1,mesh%elem2D
    elnodes=mesh%elem2D_nodes(:,n)
    q1=4
    if(elnodes(1)==elnodes(4)) q1=3
    DO q=1,q1
       ne_num(elnodes(q))=ne_num(elnodes(q))+1
       if (ne_num(elnodes(q)) > MAX_ADJACENT ) then
          print *,'Recompile with larger value for MAX_ADJACENT.'
          stop
       else
          ne_pos(ne_num(elnodes(q)),elnodes(q))=n
       endif
    END Do
END DO                 ! neighbor elements are found 

! count neighbour nodes				     
! In quads we should count the nodes that are 
! connected by edges!
allocate(aux1(mesh%nod2D))                   
aux1=0

DO n=1, mesh%nod2D
    counter=0
    DO k=1, ne_num(n)
        elem=ne_pos(k,n)
        elnodes=mesh%elem2D_nodes(:,elem)        
        if(elnodes(1)==elnodes(4)) then
            DO q=1,3
                if(elnodes(q)==n) CYCLE 
                if(aux1(elnodes(q)).ne.1) then
                    counter=counter+1
                    aux1(elnodes(q))=1
                    temp(counter)=elnodes(q)
                end if
            END DO
        else
        ! Find the position of n in elnodes:
            if(elnodes(1)==n .or. elnodes(3)==n) then
                ed(1)=elnodes(2)
                ed(2)=elnodes(4)
            else
                ed(1)=elnodes(1)
                ed(2)=elnodes(3)
            end if
            DO q=1,2
                if(aux1(ed(q)).ne.1) then
                    counter=counter+1
                    aux1(ed(q))=1
                    temp(counter)=ed(q)
                end if
            END DO
        end if
    END DO
    nn_num(n)=counter         
    aux1(temp(1:counter))=0
END DO              
 
allocate(nn_pos(maxval(nn_num)+1,mesh%nod2D))
nn_pos = -1
aux1=0
 
DO n=1, mesh%nod2D
    counter=0
    DO k=1, ne_num(n)
        elem=ne_pos(k,n)
        elnodes=mesh%elem2D_nodes(:,elem)        
        if(elnodes(1)==elnodes(4)) then 
            DO q=1,3
                if(elnodes(q)==n) CYCLE 
                if(aux1(elnodes(q)).ne.1) then
                    counter=counter+1
                    aux1(elnodes(q))=1
                    temp(counter)=elnodes(q)
                end if
            END DO
        else
            ! Find the position of n in elnodes:
            if(elnodes(1)==n .or. elnodes(3)==n) then
                ed(1)=elnodes(2)
                ed(2)=elnodes(4)
            else
                ed(1)=elnodes(1)
                ed(2)=elnodes(3)
            end if
            DO q=1,2
                if(aux1(ed(q)).ne.1) then
                    counter=counter+1
                    aux1(ed(q))=1
                    temp(counter)=ed(q)
                end if
            END DO
        end if
    END DO
    nn_num(n)=counter+1
    aux1(temp(1:counter))=0
    nn_pos(1,n)=n 
    nn_pos(2:counter+1,n)=temp(1:counter)
END DO              
deallocate(aux1)
! neighboring nodes are found. First in the list is the node itself 

! ====================
! (b) Find edges and elements containing them.
!     Write information to auxiliary file
! ====================
! open(10, file='edges.out')

 ! Count edges: 
 ! ==================== 
 ! form edges with n by cycling over neighbor
 ! nodes (if edges are not accounted yet). 
 ! New edges are added only if neighbor>n  
 ! ====================
 counter = 0
 DO n=1,mesh%nod2D
    counter = counter + count(nn_pos(2:nn_num(n),n)>nn_pos(1,n))
 end do
 mesh%edge2D=counter

 allocate(mesh%edges   (2, mesh%edge2D))
 allocate(mesh%edge_tri(2, mesh%edge2D))
 counter_in=0 

 DO n=1,mesh%nod2D
    DO q=2,nn_num(n)
       node=nn_pos(q,n)
       if(node<n) CYCLE
       counter=counter+1   ! new edge (n,node)
       ! find elements containing n and node
       flag=0
       DO k=1, ne_num(n)
          elem=ne_pos(k,n)
          elnodes=mesh%elem2D_nodes(:,elem)
          q2=4
          if(elnodes(1)==elnodes(4)) q2=3
          DO q1=1,q2
             if (elnodes(q1)==node) then
                flag=flag+1
                elems(flag)=elem
                EXIT
             end if
          END DO
       END DO
       if(flag==2) then
          counter_in=counter_in+1
          mesh%edges(1,counter_in)=n
          mesh%edges(2,counter_in)=node
          mesh%edge_tri(:,counter_in)=elems
!          write(10,'(4I10)') n, node, elems
!       else if (flag==1) then
!          write(10,'(4I10)') n, node, elems(1), -999
!       else
!          write(*,*) 'flag'
       end if
    END DO
 END DO
 mesh%edge2D_in=counter_in

 ! Repeat to collect boundary edges:   
 counter=0
 DO n=1,mesh%nod2D
    DO q=2,nn_num(n)
       node=nn_pos(q,n)
       if(node<n) CYCLE
       counter=counter+1   ! new edge (n,node)
       ! find triangles containing n and node
       flag=0
       DO k=1, ne_num(n)
          elem=ne_pos(k,n)
          elnodes=mesh%elem2D_nodes(:,elem)
          q2=4
          if(elnodes(1)==elnodes(4)) q2=3
          DO q1=1,q2
             if (elnodes(q1)==node) then
                flag=flag+1
                elems(flag)=elem
                EXIT
             end if
          END DO
       END DO
       if(flag==1) then
          counter_in=counter_in+1
          mesh%edges(1,counter_in)=n
          mesh%edges(2,counter_in)=node
          elems(2)=-999
          mesh%edge_tri(:,counter_in)=elems
       end if
    END DO
 END DO
 !close(10)
 ! Edges from edge2D_in+1 to edge2D lie on the horizontal boundary
 ! The rest (1:edge2D_in) are internal edges

 ! ====================
 ! (d) the list of elements on both sides of edge e, edge_tri(:,e)
 ! should be ordered so that the first is a triangle which is to the left of 
 ! the edge (vector pointing from the first to the second node of the edge. 
 ! If the edge is on the boundary, there is only the first triangle.  
 ! ====================

 DO n=1, mesh%edge2D
    ed=mesh%edges(:,n)
    if (mesh%edge_tri(1,n)<=0) then
       elem=mesh%edge_tri(1,n)
       elem1=mesh%edge_tri(2,n)
       mesh%edge_tri(1,n)=elem1
       mesh%edge_tri(2,n)=elem
    endif
    call elem_center(mesh%edge_tri(1,n), xc(1), xc(2), mesh)
    xc=xc-mesh%coord_nod2D(:,ed(1))
    xe=mesh%coord_nod2D(:,ed(2))-mesh%coord_nod2D(:,ed(1))
    call trim_cyclic(xe(1))
    call trim_cyclic(xc(1))
    if(xc(1)*xe(2)-xc(2)*xe(1)>0.0_WP) then
       ! Vector drawn to the center of the first triangle is to the right
       ! of the edge vector. Triangles have to be exchanged:
       elem=mesh%edge_tri(1,n)
       elem1=mesh%edge_tri(2,n)
       if(elem1>0) then  !TODO
          mesh%edge_tri(1,n)=elem1
          mesh%edge_tri(2,n)=elem
       else
          if (elem<=0) write(*,*) '2 neighbouring elems are on the ground.'
          elem=mesh%edges(2,n)         ! change the order of nodes
          mesh%edges(2,n)=mesh%edges(1,n)
          mesh%edges(1,n)=elem
       end if
    end if
 END DO

 ! ====================
 ! (e) We need an array inverse to edge_tri listing edges
 ! of a given triangle 
 ! ====================
 allocate(mesh%elem_edges(4,mesh%elem2D))
 allocate(aux1(mesh%elem2D))
 aux1=0
 DO n=1, mesh%edge2D
    DO k=1,2
       q=mesh%edge_tri(k,n)   ! triangle number
       if (q>0) then
	  aux1(q)=aux1(q)+1
	  mesh%elem_edges(aux1(q),q)=n
       end if
    END DO
 END DO
 deallocate(aux1)
 ! We order the edges in this list so that they
 ! are listed in the same rotation sense as nodes.
 ! First is the edge formed by elnodes(1:2), and so on
 DO elem=1,mesh%elem2D
    elnodes=mesh%elem2D_nodes(:,elem)
    q1=4
    if(elnodes(1)==elnodes(4)) q1=3
    eledges=mesh%elem_edges(:,elem)
    DO q=1,q1-1
       DO k=1,q1
          if(((mesh%edges(1,eledges(k))==elnodes(q)).and. &
               (mesh%edges(2,eledges(k))==elnodes(q+1))).or. &
               ((mesh%edges(1,eledges(k))==elnodes(q+1)).and. &
               (mesh%edges(2,eledges(k))==elnodes(q)))) then
             mesh%elem_edges(q,elem)=eledges(k)
             exit
          end if
       END DO
    END DO
    DO k=1,q1
       if(((mesh%edges(1,eledges(k))==elnodes(q1)).and. &
            (mesh%edges(2,eledges(k))==elnodes(1))).or. &
            ((mesh%edges(1,eledges(k))==elnodes(1)).and. &
            (mesh%edges(2,eledges(k))==elnodes(q1)))) then
          mesh%elem_edges(q1,elem)=eledges(k)
          exit
       end if
    END DO
    if(q1==3) mesh%elem_edges(4,elem)=mesh%elem_edges(1,elem)
 END DO

 !> The edge and elem lists agree in the sense that edge 1 does not
 !! contain node 1 and so on
 ! Never overwrite edge files that already exist in the mesh directory: existing
 ! experiments (and their restarts) depend on them. A newer algorithm rewriting
 ! them in place silently invalidates every restart written before. Delete the
 ! files by hand if you really want to regenerate them.
 inquire(file=trim(meshpath)//'edgenum.out',  exist=file_exist1)
 inquire(file=trim(meshpath)//'edges.out',    exist=file_exist2)
 inquire(file=trim(meshpath)//'edge_tri.out', exist=file_exist3)
 if (file_exist1 .or. file_exist2 .or. file_exist3) then
    print *, achar(27)//'[1;33m'//' WARNING: edgenum.out/edges.out/edge_tri.out already exist in '//trim(meshpath)
    print *, '          keeping the EXISTING files untouched; the partitioning itself uses the'
    print *, '          newly computed in-memory edges. Delete the files first to regenerate them.'//achar(27)//'[0m'
 else
    open(11, file=trim(meshpath)//'edgenum.out')
    write(11,*) mesh%edge2D
    write(11,*) mesh%edge2D_in
    close(11)
    open(10, file=trim(meshpath)//'edges.out')
    open(12, file=trim(meshpath)//'edge_tri.out')
    do n=1,mesh%edge2D
       write(10,*) mesh%edges(:,n)
       write(12,*) mesh%edge_tri(:,n)
    end do
    close(10)
    close(12)
 end if
 deallocate(ne_num, ne_pos)
END SUBROUTINE find_edges_ini
!=========================================================================
!> @brief
!> Finds elemental and nodal levels.
!> Does some thresholding: if (depth>zbar(4)) x=zbar(4)
!> Fixes rough topography, by converting some oceans cells to ground cell(reflected by changing levels arrays)
!> Creates 2 files: elvls.out, nlvls.out
subroutine find_levels(mesh)
    implicit none
    INTEGER :: nodes(3), elems(3), eledges(3)
    integer :: elem1, j, n, nneighb, q, node, i, nz, auxi
    integer :: count_iter, count_neighb_open, exit_flag, fileID=111
    real(kind=WP) :: x, dmean
    logical :: file_exist
    integer :: max_iter=1000
    character(MAX_PATH) :: file_name
    type(t_mesh), intent(inout), target :: mesh


    print *, achar(27)//'[1m'  //'____________________________________________________________'//achar(27)//'[0m'
    print *, achar(27)//'[7;1m' //' -->: read bottom depth                                     '//achar(27)//'[0m'
    
    !___________________________________________________________________________
    ! allocate depth
    if (use_depthonelem) then 
        allocate(mesh%depth(mesh%elem2D))
    else    
        allocate(mesh%depth(mesh%nod2D))
    end if     
    
    !______________________________________________________________________________
    ! read depth from aux3d.out
    if (trim(use_depthfile)=='aux3d') then
        ! check if aux3d.out file does exist
        file_exist=.False.
        file_name=trim(meshpath)//'aux3d.out' 
        inquire(file=trim(file_name),exist=file_exist)
        !_______________________________________________________________________
        if (file_exist) then
            write(*," (A, A)") ' read file:',trim(file_name)
            !___________________________________________________________________
            ! load fesom2.0 aux3d.out file 
            open(fileID, file=file_name)
            
            ! read the number of levels 
            read(fileID,*) mesh%nl          
            allocate(mesh%zbar(mesh%nl))         ! their standard depths
            
            ! read full depth levels
            read(fileID,*) mesh%zbar
            if(mesh%zbar(2)>0) mesh%zbar=-mesh%zbar   ! zbar is negative 
            
            ! compute mid depth levels
            allocate(mesh%Z(mesh%nl-1))
            mesh%Z=mesh%zbar(1:mesh%nl-1)+mesh%zbar(2:mesh%nl)  ! mid-depths of cells
            mesh%Z=0.5_WP*mesh%Z
        else
            write(*,*) '____________________________________________________________________'
            write(*,*) ' ERROR: You want to use aux3d.out file to define your depth, but '
            write(*,*) '        the file seems not to exist'
            write(*,*) '        --> check in namelist.config, the flag use_depthfile must be'
            write(*,*) '            use_depthfile= "aux3d" or "depth@" '
            write(*,*) '        --> model stops here'
            write(*,*) '____________________________________________________________________'
            stop
        end if 
    !___________________________________________________________________________
    ! read depth from depth@node.out or depth@elem.out    
    elseif (trim(use_depthfile)=='depth@') then
        !_______________________________________________________________________
        ! load file depth_zlev.out --> contains number of model levels and full depth
        ! levels
        file_exist=.False.
        file_name=trim(meshpath)//'depth_zlev.out' 
        inquire(file=trim(file_name),exist=file_exist) 
        if (file_exist) then
            write(*," (A, A)") ' read file:',trim(file_name)
            !___________________________________________________________________
            ! load fesom2.0 aux3d.out file 
            open(fileID, file=file_name)
            
            ! read the number of levels 
            read(fileID,*) mesh%nl          
            allocate(mesh%zbar(mesh%nl))         ! their standard depths
            
            ! read full depth levels
            read(fileID,*) mesh%zbar
            if(mesh%zbar(2)>0) mesh%zbar=-mesh%zbar   ! zbar is negative 
            
            ! compute mid depth levels
            allocate(mesh%Z(mesh%nl-1))
            mesh%Z=mesh%zbar(1:mesh%nl-1)+mesh%zbar(2:mesh%nl)  ! mid-depths of cells
            mesh%Z=0.5_WP*mesh%Z
            
            close(fileID)
        else
            write(*,*) '____________________________________________________________________'
            write(*,*) ' ERROR: You want to use depth@elem.out or depth@node.out file, therefore'
            write(*,*) '        you also need the file depth_zlev.out which contains the model '
            write(*,*) '        number of layers and the depth of model levels. This file seems '
            write(*,*) '        not to exist'
            write(*,*) '        --> check in namelist.config, the flag use_depthfile must be'
            write(*,*) '            use_depthfile= "aux3d" or "depth@" and your meshfolder'
            write(*,*) '        --> model stops here'
            write(*,*) '____________________________________________________________________'
            stop
        endif 
        
        !_______________________________________________________________________
        ! load file depth@elem.out or depth@node.out contains topography either at 
        ! nodes or elements
        if (use_depthonelem) then 
            file_name=trim(meshpath)//'depth@elem.out' 
        else
            file_name=trim(meshpath)//'depth@node.out' 
        end if 
        inquire(file=trim(file_name),exist=file_exist) 
        if (file_exist) then 
            write(*," (A, A)") ' read file:',trim(file_name)
            open(fileID, file=file_name)
        else    
            write(*,*) '____________________________________________________________________'
            write(*,*) ' ERROR: You want to use depth@elem.out or depth@node.out file to '
            write(*,*) '        define your depth, but the file seems not to exist'
            write(*,*) '        --> check in namelist.config, the flag use_depthfile must be'
            write(*,*) '            use_depthfile= "aux3d" or "depth@" and your meshfolder '
            write(*,*) '        --> model stops here'
            write(*,*) '____________________________________________________________________'
            stop
        end if 
    end if 
    
    !___________________________________________________________________________
    ! read topography from file
    auxi = mesh%nod2d
    if (use_depthonelem) auxi = mesh%elem2d
!     write(*,*) ' use_depthonelem = ',use_depthonelem
!     write(*,*)  ' auxi =',auxi 
    DO n = 1, auxi
        read(fileID,*) x
        if (x>0) x=-x
        if (x>mesh%zbar(thers_zbar_lev)) x=mesh%zbar(thers_zbar_lev) !TODO KK thresholding for depth
        mesh%depth(n)=x
    END DO
    close(fileID)
    if(mesh%depth(2)>0) mesh%depth=-mesh%depth  ! depth is negative

    !___________________________________________________________________________
    print *, achar(27)//'[1m'  //'____________________________________________________________'//achar(27)//'[0m'
    print *, achar(27)//'[7;1m' //' -->: compute elem, vertice bottom depth index              '//achar(27)//'[0m'
     
    allocate(mesh%nlevels(mesh%elem2D))
    allocate(mesh%nlevels_nod2D(mesh%nod2D))

    !___________________________________________________________________________
    ! Compute the initial number number of elementa levels, based on the vertice
    ! depth information 
    do n=1, mesh%elem2D
        
        !_______________________________________________________________________
        if (use_depthonelem) then 
            dmean = mesh%depth(n) ! depth is already defined on elements 
        else
            nodes=mesh%elem2D_nodes(1:3,n)
            !___________________________________________________________________
            ! depth of element is  shallowest depth of sorounding vertices
            if     (trim(which_depth_n2e) .eq. 'min') then ; dmean=maxval(mesh%depth(nodes))
            ! depth of element is deepest depth of sorounding vertices    
            elseif (trim(which_depth_n2e) .eq. 'max') then ; dmean=minval(mesh%depth(nodes))
            ! DEFAULT: depth of element is  mean depth of sorounding vertices
            elseif (trim(which_depth_n2e) .eq. 'mean') then; dmean=sum(mesh%depth(nodes))/3.0
            end if 
        end if 
        
        !_______________________________________________________________________
        exit_flag=0
        do nz=1,mesh%nl-1
            if(mesh%Z(nz)<dmean) then
                exit_flag=1
                mesh%nlevels(n)=nz
                exit
            end if
        end do
        if((exit_flag==0).and.(dmean<0)) mesh%nlevels(n)=mesh%nl
        if(dmean>=0) mesh%nlevels(n)=thers_zbar_lev
        
        ! set minimum number of levels to --> thers_lev=5
        if(mesh%nlevels(n)<thers_zbar_lev) mesh%nlevels(n)=thers_zbar_lev
    end do ! --> do n=1, elem2D
    
    !___________________________________________________________________________
    ! write out vertical level indices before iterative geometric adaption to 
    ! exclude isolated cells 
    !_______________________________________________________________________
    file_name=trim(meshpath)//'elvls_raw.out'
    inquire(file=trim(file_name), exist=file_exist)
    if (file_exist) then
        print *, achar(27)//'[1;33m'//' WARNING: '//trim(file_name)//' already exists, keeping the existing file'//achar(27)//'[0m'
    else
        open(fileID, file=file_name)
          do n=1,mesh%elem2D
             write(fileID,*) mesh%nlevels(n)
          end do
        close(fileID)
    end if
    
    !___________________________________________________________________________
    ! check for isolated cells (cells with at least two boundary faces or three 
    ! boundary vertices) and eliminate them --> FESOM2.0 doesn't like these kind
    ! of cells 
    do nz=thers_zbar_lev+1,mesh%nl
        exit_flag=0
        count_iter=0
        
        !_______________________________________________________________________
        ! iteration loop within each layer
        do while((exit_flag==0).and.(count_iter<max_iter))
            exit_flag=1
            count_iter=count_iter+1
            
            !___________________________________________________________________
            ! loop over triangles
            do n=1,mesh%elem2D
                ! merge: result = merge(truesource, falsesource, mask)
                ! --> if elem2D_nodes(1,n) == elem2D_nodes(4,n): True  --> q=3 --> triangular mesh
                ! --> if elem2D_nodes(1,n) == elem2D_nodes(4,n): False --> q=4 --> quad mesh
                nneighb = merge(3,4,mesh%elem2D_nodes(1,n) == mesh%elem2D_nodes(4,n))
                !                       
                !                         +---isolated bottom cell
                !   ._______________      |          _______________________.
                !   |###|###|###|###|___  |      ___|###|###|###|###|###|###|
                !   |###|###|###|###|###| |  ___|###|###|###|###|###|###|###|
                !   |###|###|###|###|###| | |###|###|###|###|  BOTTOM   |###|
                !   |###|###|###|###|###|_v_|###|###|###|###|###|###|###|###|
                !   |###|###|###|###|###|###|###|###|###|###|###|###|###|###|
                !
                if (mesh%nlevels(n)>=nz) then
                    count_neighb_open=0
                    elems=mesh%elem_neighbors(1:3,n)
                    !___________________________________________________________
                    ! loop over neighbouring triangles
                    do i=1,nneighb
                        if (elems(i)>0) then
                            if (mesh%nlevels(elems(i))>=nz) then
                                !count neighbours
                                count_neighb_open=count_neighb_open+1
                            endif
                        endif
                    enddo
                    
                    !___________________________________________________________
                    ! check how many open faces to neighboring triangles the cell 
                    ! has, if there are less than 2 its isolated (a cell should 
                    ! have at least 2 valid neighbours)
                    if (count_neighb_open<2) then
                        ! if cell is "isolated", and the one levels shallower bottom 
                        ! cell would be shallower than the minimum vertical level 
                        ! treshhold (thers_lev). --> in this make sorrounding elements 
                        ! one level deeper to reconnect the isolated cell
                        if (nz-1<thers_zbar_lev) then 
                            do i=1,nneighb
                                if (elems(i)>0) then
                                    mesh%nlevels(elems(i)) = max(mesh%nlevels(elems(i)),nz)
                                end if     
                            end do    
                            
                        !if cell is "isolated" convert to one level shallower bottom cell
                        else
                            mesh%nlevels(n)=nz-1
                        end if  
                        !force recheck for all current ocean cells
                        exit_flag=0
                        
                    end if
                end if ! --> if (nlevels(n)>=nz) then
            end do ! --> do n=1,elem2D
        end do ! --> do while((exit_flag==0).and.(count1<1000))
        write(*,"(A, I5, A, i5, A, I3)") '  -[iter ]->: nlevel, iter/maxiter=',count_iter,'/',max_iter,', nz=',nz
    end do ! --> do nz=4,nl

    !___________________________________________________________________________
    ! vertical vertice level index of ocean bottom boundary
    write(*,"(A)"                  ) '  -[compu]->: nlevels_nod2D '
    mesh%nlevels_nod2D=0
    do n=1,mesh%elem2D
        q = merge(3,4,mesh%elem2D_nodes(1,n) == mesh%elem2D_nodes(4,n))
        do j=1,q
            node=mesh%elem2D_nodes(j,n)
            if(mesh%nlevels_nod2D(node)<mesh%nlevels(n)) then
            mesh%nlevels_nod2D(node)=mesh%nlevels(n)
            end if
        end do
    end do

    !___________________________________________________________________________
    ! write vertical level indices into file
    ! Never overwrite elvls.out/nlvls.out if they already exist: the runtime model
    ! reads them, and every restart ever written with the existing files depends on
    ! the level structure they define. Overwriting them with the output of a newer
    ! (even slightly different) level-finding algorithm makes old restarts blow up
    ! with zero-thickness layers. Delete the files by hand to regenerate them.
    inquire(file=trim(meshpath)//'elvls.out', exist=file_exist)
    if (.not. file_exist) inquire(file=trim(meshpath)//'nlvls.out', exist=file_exist)
    if (file_exist) then
        print *, achar(27)//'[1;33m'//' WARNING: elvls.out/nlvls.out already exist in '//trim(meshpath)
        print *, '          keeping the EXISTING files untouched; the partitioning itself uses the'
        print *, '          newly computed in-memory levels (only load-balance weights are affected).'
        print *, '          Delete the files first if you really want to regenerate them.'//achar(27)//'[0m'
    else
        file_name=trim(meshpath)//'elvls.out'
        open(fileID, file=file_name)
        do n=1,mesh%elem2D
           write(fileID,*) mesh%nlevels(n)
        end do
        close(fileID)

        file_name=trim(meshpath)//'nlvls.out'
        open(fileID, file=file_name)
          do n=1,mesh%nod2D
             write(fileID,*) mesh%nlevels_nod2D(n)
          end do
        close(fileID)
    end if
        
    !_______________________________________________________________________
    write(*,*) '========================='
    write(*,*) 'Mesh is read : ', 'nod2D=', mesh%nod2D,' elem2D=', mesh%elem2D, ' nl=', mesh%nl
    write(*,*) 'Min/max depth on mype: ',  -mesh%zbar(minval(mesh%nlevels)),-mesh%zbar(maxval(mesh%nlevels))
    write(*,*) '3D tracer nodes on mype ', sum(mesh%nlevels_nod2d)-(mesh%elem2D)
    write(*,*) '========================='
end subroutine find_levels

!
!
!_______________________________________________________________________________
! finds elemental and nodal levels of cavity-ocean boundary.
! Creates 2 files: cavity_elvls.out, cavity_nlvls.out
subroutine find_levels_cavity(mesh)
    implicit none
    integer        :: nodes(3), elems(3)
    integer        :: elem, node, nz, j
    integer        :: min_nlvl, max_ulvl, idx, idx2, val, val2
    integer        :: count_neighb_open, nneighb, cavity_maxlev, count_isoelem
    integer        :: exit_flag1, count_iter, max_iter=1000, exit_flag2, count_iter2, max_iter2=25
    real(kind=WP)  :: dmean
    character(MAX_PATH) :: file_name
    logical        :: file_exist
    integer, allocatable, dimension(:)   :: aux_arr, aux_idx
    integer, allocatable, dimension(:)   :: numelemtonode
    logical, allocatable, dimension(:)   :: elemreducelvl, elemfixlvl
    type(t_mesh), intent(inout), target  :: mesh
    !___________________________________________________________________________
    print *, achar(27)//'[1m'  //'____________________________________________________________'//achar(27)//'[0m'
    print *, achar(27)//'[7;1m' //' -->: compute elem, vertice cavity depth index               '//achar(27)//'[0m'
    
    !___________________________________________________________________________
    allocate(mesh%ulevels(mesh%elem2D))
    allocate(mesh%ulevels_nod2D(mesh%nod2D))
        
    !___________________________________________________________________________
    ! Compute level position of ocean-cavity boundary
    cavity_maxlev=0
    do elem=1, mesh%elem2D
        
        !_______________________________________________________________________
        if (use_cavityonelem) then 
            dmean = mesh%cavity_depth(elem)
        else
            nodes=mesh%elem2D_nodes(1:3,elem)
            !_______________________________________________________________________
            ! depth of element is  shallowest depth of sorounding vertices
            if     (trim(which_depth_n2e) .eq. 'min')  then ; dmean=maxval(mesh%cavity_depth(nodes))
            ! depth of element is deepest depth of sorounding vertices    
            elseif (trim(which_depth_n2e) .eq. 'max')  then ; dmean=minval(mesh%cavity_depth(nodes))
            ! DEFAULT: depth of element is  mean depth of sorounding vertices
            elseif (trim(which_depth_n2e) .eq. 'mean') then ; dmean=sum(mesh%cavity_depth(nodes))/3.0
            end if 
        end if 
        
        !_______________________________________________________________________
        ! vertical elem level index of cavity-ocean boundary
        mesh%ulevels(elem) = 1
        if (dmean<0.0_WP) mesh%ulevels(elem) = 2
        
        do nz=1,mesh%nlevels(elem)-1
            if (mesh%Z(nz)<dmean .or. mesh%nlevels(elem)-nz<=3) then
                mesh%ulevels(elem)=nz         !    to compute shechpetkin PGF
                exit
            end if
        end do
        cavity_maxlev = max(cavity_maxlev,mesh%ulevels(elem))
    end do
    
    !___________________________________________________________________________
    ! write out cavity mesh files for vertice and elemental position of 
    ! vertical cavity-ocean boundary before the iterative geometric adaption to 
    ! eliminate isolated cells
    ! write out elemental cavity-ocean boundary level
    file_name=trim(meshpath)//'cavity_elvls_raw.out'
    inquire(file=trim(file_name), exist=file_exist)
    if (file_exist) then
        print *, achar(27)//'[1;33m'//' WARNING: '//trim(file_name)//' already exists, keeping the existing file'//achar(27)//'[0m'
    else
        open(20, file=file_name)
        do elem=1,mesh%elem2D
            write(20,*) mesh%ulevels(elem)
        enddo
        close(20)
    end if
    
    !___________________________________________________________________________
    ! Eliminate cells that have two cavity boundary faces --> should not be 
    ! possible in FESOM2.0
    ! loop over all cavity levels
    allocate(elemreducelvl(mesh%elem2d),elemfixlvl(mesh%elem2d))
    allocate(numelemtonode(mesh%nl))
        
    !___________________________________________________________________________
    ! outer iteration loop    
    count_iter2 = 0
    exit_flag2  = 0
    elemfixlvl  = .False.
    do while((exit_flag2==0) .and. (count_iter2<max_iter2))
        count_iter2   = count_iter2+1
        elemreducelvl = .False.
        !_______________________________________________________________________
        ! Loop oper cavity occupying levels
        do nz=1,cavity_maxlev
            exit_flag1 = 0
            count_iter = 0
            !___________________________________________________________________
            ! iteration loop within each layer
            do while((exit_flag1==0).and.(count_iter<max_iter))
                exit_flag1    = 1
                count_iter    = count_iter+1
                count_isoelem = 0
                !_______________________________________________________________
                ! loop over triangles
                do elem=1,mesh%elem2D
                    ! nneighb = 3 --> tri mesh, nneighb = 4 --> quad mesh
                    nneighb = merge(3,4,mesh%elem2D_nodes(1,elem) == mesh%elem2D_nodes(4,elem))
                    elems   = mesh%elem_neighbors(1:3,elem)
                    !
                    !   .___________________________.~~~~~~~~~~~~~~~~~~~~~~~~~~
                    !   |###|###|###|###|###|###|###|
                    !   |#  CAVITY  |###| . |###|###|                    OCEAN
                    !   |###|###|###|###|/|\|###| 
                    !   |###|###|         |
                    !   |###|             +-- Not good can lead to isolated cells  
                    !
                    
                    !___________________________________________________________
                    ! (1st) Ask the Question: is nz at element elem an here an 
                    ! valid layer in the ocean
                    if ( nz >= mesh%ulevels(elem) .and. nz<mesh%nlevels(elem)) then
                        count_neighb_open=0
                        
                        !_______________________________________________________
                        ! (2nd) loop over the neighbouring elements of this valid ocean 
                        ! prism and check if their are further valid ocean element
                        ! neighbours at this level nz
                        ! loop over neighbouring triangles
                        do j = 1, nneighb
                            if (elems(j)>0) then ! if its a valid boundary triangle, 0=missing value
                                ! check for isolated cell
                                if ( mesh%ulevels(elems(j))<= nz .and. & 
                                     mesh%nlevels(elems(j))> nz ) then
                                    !count the open faces to neighboring cells 
                                    count_neighb_open=count_neighb_open+1
                                endif
                            end if 
                        end do ! --> do i = 1, nneighb
                        
                        !_______________________________________________________
                        ! (3rd) check how many open faces to neighboring triangles the cell 
                        ! has, if there are less than 2 its isolated (a cell should 
                        ! have at least 2 valid neighbours)
                        ! --> in this case shift cavity-ocean interface one level down
                        if (count_neighb_open<=1) then
                            count_isoelem = count_isoelem+1
                            
                            ! (4th.1 ) if cell elem is isolated convert it to a deeper ocean level
                            ! except when this levels would remain less than 3 valid 
                            ! bottom levels --> in case make the levels of all sorounding
                            ! triangles shallower
                            if ( (mesh%nlevels(elem)-(nz+1))>=3           .and. &
                                 (elemreducelvl(elem) .eqv. .False.) .and. &
                                 (elemfixlvl(   elem) .eqv. .False.)       &
                                ) then 
                                mesh%ulevels(elem)=nz+1
                                
                            ! (4th.2) can not increase depth anymore to eleminate 
                            ! isolated cell, otherwise lessthan 3 valid layers 
                            ! therefor reduce depth of ONE!!! of the neighbouring 
                            ! triangles. Choose triangle whos depth is already closest
                            ! to nz    
                            else 
                                !PS replace this with do j=1,3... because of 
                                !PS indice -999 conflict in elems, ulevels(-999)
                                !PS not possible 
                                !PS idx = minloc(ulevels(elems)-nz, 1, MASK=( (elems>0) .and. ((ulevels(elems)-nz)>0) ) )
                                val=mesh%nl
                                do j = 1, 3
                                    if (elems(j)>0) then ! if its a valid boundary triangle, 0=missing value
                                        if (mesh%ulevels(elems(j))-nz>0 .and. mesh%ulevels(elems(j))-nz<val) then
                                            val=mesh%ulevels(elems(j))-nz
                                            idx=j
                                        end if 
                                    end if 
                                end do ! --> do i = 1, nneighb
                                
                                mesh%ulevels(      elems(idx)) = nz
                                elemreducelvl(elems(idx)) = .True.
                            end if    
                            
                            ! force recheck for all current ocean cells
                            exit_flag1=0
                        end if ! --> if (count_neighb_open<2) then
                        
                    end if ! --> if ( nz >= ulevels(elem) .and. nz<nlevels(elem)) then
                    
                end do ! --> do elem=1,elem2D  
                
            end do ! --> do while((exit_flag==0).and.(count_iter<1000))
            write(*,"(A, I5, A, i5, A, I3)") '  -[iter ]->: ulevel, iter/maxiter=',count_iter,'/',max_iter,', nz=',nz
        end do ! --> do nz=1,cavity_maxlev 
        
        !_______________________________________________________________________
        ! compute vertical vertice level index of cavity_ocean boundary
        write(*,"(A)"                  ) '  -[compu]->: ulevels_nod2D '
        mesh%ulevels_nod2D = mesh%nl
        do elem=1,mesh%elem2D
            nneighb = merge(3,4,mesh%elem2D_nodes(1,elem) == mesh%elem2D_nodes(4,elem))
            !___________________________________________________________________
            ! loop over neighbouring triangles
            do j=1,nneighb
                node=mesh%elem2D_nodes(j,elem)
                mesh%ulevels_nod2D(node)=min(mesh%ulevels_nod2D(node),mesh%ulevels(elem))
            end do
        end do ! --> do elem=1,elem2D
        
        !_______________________________________________________________________
        ! check if all constrains for ulevel and ulevels_nod2D is fullfilled
        exit_flag2 = 1
        do elem=1,mesh%elem2D
            if (mesh%ulevels(elem)>=mesh%nlevels(elem)) then 
                write(*,*) ' -[check]->: elem cavity depth deeper or equal bottom depth, elem=',elem                
                exit_flag2 = 0
            end if 
            
            if (mesh%nlevels(elem)-mesh%ulevels(elem)<3) then 
                write(*,*) ' -[check]->: less than three valid elem ocean layers, elem=',elem    
                exit_flag2 = 0
            end if 
        end do ! --> do elem=1,elem2D
        
        do node=1,mesh%nod2D
            if (mesh%ulevels_nod2D(node)>=mesh%nlevels_nod2D(node)) then 
                write(*,*) ' -[check]->: vertice cavity depth deeper or equal bottom depth, node=', node
                exit_flag2 = 0
            end if
            
            if (mesh%nlevels_nod2D(node)-mesh%ulevels_nod2D(node)<3) then 
                write(*,*) ' -[check]->: less than three valid vertice ocean layers, node=', node
                exit_flag2 = 0
            end if
        end do ! --> do node=1,nod2D
        
        do elem=1,mesh%elem2D
            if (mesh%ulevels(elem)< maxval(mesh%ulevels_nod2D(mesh%elem2D_nodes(:,elem))) ) then 
                write(*,*) ' -[check]->: found elem cavity shallower than its valid maximum cavity vertice depths, elem=', mesh%elem2d
                exit_flag2 = 0
            end if 
        end do ! --> do elem=1,elem2D
        
        !_______________________________________________________________________
        ! compute how many triangle elements contribute to every vertice in every 
        ! layer
        !
        ! --> What can happen is that a node point in the middle of the vertical 
        !     domain can become isolated due to the cavity constrains. The model
        !     would not be able to deal with this kind of situation. So we must 
        !     prevent it by adapting ulevels!
        !                                 O  node 
        !                                _._ 
        !                              _/ | \_
        !                            _/   |   \_
        !                          _/     |     \_
        !                      elem(1)  elem(2)   elem(3)...  <--elem=nod_in_elem2D(j,node)
        !                             ._______. ulevel(elem2)=30
        !                             |_______|
        !                             |_______|
        !                             |_______|
        !                             |_______|
        !                             |_______| nlevel(elem2)=38
        !
        !                 In this possible gap node points
        !                 would have no neighboring elements
        !                 
        !   ulevel(elem1)=42 ._______.         ._______. ulevel(elem3)=42
        !                    |_______|         |_______|
        !                    |_______|         |_______|
        !                    |_______|         |_______|
        !                    |_______|         |_______|
        !   nlevel(elem1)=46 |_______|         |_______|
        !                                      |_______| nlevel(elem3)=48
        !
        ! --> Problem here is we want to keep nlevels fixed so what we can do is
        !     to set ulevels(elem1) and ulevels(elem3) towards nlevel(elem2)
        count_iter=0
        do node=1, mesh%nod2D
            !___________________________________________________________________
            ! check if there is a possible gap as described above that would
            ! allow for node points without neighbors
            min_nlvl = mesh%nl
            max_ulvl = 1
            do j=1, mesh%nod_in_elem2D_num(node)
                elem=mesh%nod_in_elem2D(j, node)
                min_nlvl = min(min_nlvl, mesh%nlevels(elem))
                max_ulvl = max(max_ulvl, mesh%ulevels(elem))
            end do    
            
            ! found a potential gap
            if (min_nlvl < max_ulvl) then 
                
                !_______________________________________________________________
                ! compute how many triangle elements contribute to vertice in 
                ! every layer check if there are layers where a node point has 
                ! only one or even zero neighboring triangles.
                numelemtonode=0
                do j=1, mesh%nod_in_elem2D_num(node)
                    elem=mesh%nod_in_elem2D(j, node)
                    do nz=mesh%ulevels(elem), mesh%nlevels(elem)-1
                        numelemtonode(nz) = numelemtonode(nz) + 1
                    end do
                end do
                
                !_______________________________________________________________
                ! check in which depth level is an isolated node 
                nzloop: do nz = mesh%ulevels_nod2D(node), mesh%nlevels_nod2D(node)-1
                    !___________________________________________________________
                    ! nodes has zero neighbouring triangles and is completely 
                    ! isolated need to adapt ulevels --> inflicts another 
                    ! outher iteration loop (exit_flag2=0). It needs at least 
                    ! two neighboring elements so the node is considered as 
                    ! connected. Search the index of the two elements where ulevels>nz
                    ! but that are closest to nz
                    if (numelemtonode(nz)==0) then 
                        count_iter = count_iter+1
                        write(*,"( A, I1, A, I7, A, I3)") '  -[check]->: node has only ', numelemtonode(nz) ,' neighb. triangle: n=', node, ', nz=',nz
                        !_______________________________________________________
                        allocate(aux_arr(mesh%nod_in_elem2D_num(node)), aux_idx(mesh%nod_in_elem2D_num(node)))
                        aux_arr(:) = mesh%ulevels(mesh%nod_in_elem2D(1:mesh%nod_in_elem2D_num(node),node))
                        aux_arr(:) = aux_arr(:) - nz
                        ! fill array with index of element
                        do j=1, mesh%nod_in_elem2D_num(node)
                            aux_idx(j) = j
                        end do
                        ! index of closest elem to nz where ulevel>nz
                        idx  = minloc(aux_arr, 1, MASK=((aux_arr>0)) )
                        ! index of second closest elem to nz where ulevel>nz
                        idx2 = minloc(aux_arr, 1, MASK=((aux_arr>0) .and. (aux_idx/=idx)) )
                        deallocate(aux_arr, aux_idx)
                        mesh%ulevels(   mesh%nod_in_elem2D(idx ,node)) = nz
                        mesh%ulevels(   mesh%nod_in_elem2D(idx2,node)) = nz
                        elemfixlvl(mesh%nod_in_elem2D(idx ,node)) = .True.
                        elemfixlvl(mesh%nod_in_elem2D(idx2,node)) = .True.
                        
                        !_______________________________________________________
                        ! inflict another outer iteration loop
                        exit_flag2 = 0
                        
                        !_______________________________________________________
                        ! if the upper most isolated layer is fixed all layers below should be fixed as well
                        ! --> exit loop do nz = ulevels_nod2D(node), nlevels_nod2D(node)-1
                        exit nzloop
                        
                    !___________________________________________________________
                    ! nodes has one neighbouring triangles it needs at least 
                    ! another neighboring elements so the node is considered as 
                    ! connected     
                    elseif (numelemtonode(nz)==1) then
                        count_iter = count_iter+1
                        write(*,"( A, I1, A, I7, A, I3)") '  -[check]->: node has only ', numelemtonode(nz) ,'neighb. triangle: n=', node, ', nz=',nz
                        !_______________________________________________________
                        allocate(aux_arr(mesh%nod_in_elem2D_num(node)), aux_idx(mesh%nod_in_elem2D_num(node)))
                        aux_arr(:) = mesh%ulevels(mesh%nod_in_elem2D(1:mesh%nod_in_elem2D_num(node),node))
                        aux_arr(:) = aux_arr(:) - nz
                        ! fill array with index of element
                        do j=1, mesh%nod_in_elem2D_num(node)
                            aux_idx(j) = j
                        end do
                        ! index of closest elem to nz where ulevel>nz
                        idx  = minloc(aux_arr, 1, MASK=((aux_arr>0)) )
                        deallocate(aux_arr, aux_idx)
                        mesh%ulevels(   mesh%nod_in_elem2D(idx,node)) = nz
                        elemfixlvl(mesh%nod_in_elem2D(idx,node)) = .True.
                        
                        !_______________________________________________________
                        ! inflict another outer iteration loop
                        exit_flag2 = 0
                        
                        !_______________________________________________________
                        ! if the upper most isolated layer is fixed all layers below should be fixed as well
                        ! --> exit loop do nz = ulevels_nod2D(node), nlevels_nod2D(node)-1
                        exit nzloop
                    
                    end if 
                end do nzloop ! --> do nz = ulevels_nod2D(node), nlevels_nod2D(node)-1
            end if ! --> if (min_nlvl < max_ulvl) then 
        end do ! --> do node=1, nod2D
        
        !_______________________________________________________________________
        ! check if cavity geometry did converge 
        if (exit_flag2 == 0) then 
            print *, achar(27)//'[33m'  //'____________________________________________________________'//achar(27)//'[0m'
            print *, ' -['//achar(27)//'[33m'//'WARN'//achar(27)//'[0m'//']->: Cavity geom. not converged yet, do further outer iteration'
            write(*,"(A, I3, A, I3)") '             iter step ', count_iter2,'/', max_iter2
            write(*,*)
        end if     
        
        !_______________________________________________________________________
    end do ! --> do while((exit_flag2==0) .and. (count_iter2<max_iter2))
    deallocate(elemreducelvl,elemfixlvl)
    deallocate(numelemtonode)
    
    !___________________________________________________________________________
    ! check if cavity geometry totaly converged or failed to converge in the later
    ! case will break up model
    if (exit_flag2 == 0) then 
        write(*,*)
        print *, achar(27)//'[31m'  //'____________________________________________________________'//achar(27)//'[0m'
        print *, achar(27)//'[7;31m'//' -[ERROR]->: Cavity geometry constrains did not converge !!! *\(>︿<)/*'//achar(27)//'[0m'
        write(*,*)
        stop
    else    
        write(*,*)
        print *, achar(27)//'[32m'  //'____________________________________________________________'//achar(27)//'[0m'
        print *, ' -['//achar(27)//'[7;32m'//' OK  '//achar(27)//'[0m'//']->: Cavity geometry constrains did converge !!! *\(^o^)/*'
    end if     
 
    !___________________________________________________________________________
    ! write out cavity mesh files for vertice and elemental position of 
    ! vertical cavity-ocean boundary
    ! Never overwrite cavity_elvls.out/cavity_nlvls.out if they already exist:
    ! same reasoning as for elvls.out/nlvls.out -- existing restarts depend on them.
    inquire(file=trim(meshpath)//'cavity_elvls.out', exist=file_exist)
    if (.not. file_exist) inquire(file=trim(meshpath)//'cavity_nlvls.out', exist=file_exist)
    if (file_exist) then
        print *, achar(27)//'[1;33m'//' WARNING: cavity_elvls.out/cavity_nlvls.out already exist in '//trim(meshpath)
        print *, '          keeping the EXISTING files untouched.'
        print *, '          Delete the files first if you really want to regenerate them.'//achar(27)//'[0m'
    else
        ! write out elemental cavity-ocean boundary level
        file_name=trim(meshpath)//'cavity_elvls.out'
        open(20, file=file_name)
        do elem=1,mesh%elem2D
            write(20,*) mesh%ulevels(elem)
        enddo
        close(20)

        ! write out vertice cavity-ocean boundary level + yes/no cavity flag
        file_name=trim(meshpath)//'cavity_nlvls.out'
        open(20, file=file_name)
        do node=1,mesh%nod2D
            write(20,*) mesh%ulevels_nod2D(node)
        enddo
        close(20)
    end if

end subroutine find_levels_cavity


!===================================================================

subroutine edge_center(n1, n2, x, y, mesh)
!
! Returns coordinates of edge center in x and y
! 
implicit none
integer       :: n1, n2   ! nodes of the edge
real(kind=WP) :: x, y, a(2), b(2)
type(t_mesh), intent(inout), target :: mesh

a=mesh%coord_nod2D(:,n1)
b=mesh%coord_nod2D(:,n2)
if(a(1)-b(1)>cyclic_length/2.0) a(1)=a(1)-cyclic_length
if(a(1)-b(1)<-cyclic_length/2.0) b(1)=b(1)-cyclic_length
x=0.5_WP*(a(1)+b(1))
y=0.5_WP*(a(2)+b(2))
end subroutine edge_center
!====================================================================
subroutine elem_center(elem, x, y, mesh)
!
! Returns coordinates of elem center in x and y
!
implicit none
integer, intent(in)  :: elem
integer              ::  elnodes(3), k    
real(kind=WP), intent(out) :: x, y
real(kind=WP)        ::  ax(3), amin
type(t_mesh), intent(in), target :: mesh

   elnodes=mesh%elem2D_nodes(1:3,elem)
   ax=mesh%coord_nod2D(1, elnodes)
   amin=minval(ax)
   DO k=1,3
   if(ax(k)-amin>cyclic_length/2.0) ax(k)=ax(k)-cyclic_length
   END DO
   x=sum(ax)/3.0_WP
   y=sum(mesh%coord_nod2D(2,elnodes))/3.0_WP
   
end subroutine elem_center
!=======================================================================
SUBROUTINE find_elem_neighbors_ini(mesh)
! For each element three its element neighbors are found
implicit none
integer    :: elem, eledges(3), elem1, j, n, elnodes(3)
type(t_mesh), intent(inout), target :: mesh

allocate(mesh%elem_neighbors(4,mesh%elem2D))
mesh%elem_neighbors=0
DO elem=1,mesh%elem2D
   
   eledges=mesh%elem_edges(1:3,elem)
   DO j=1,3
   elem1=mesh%edge_tri(1,eledges(j))
   if(elem1==elem) elem1=mesh%edge_tri(2,eledges(j))
   mesh%elem_neighbors(j,elem)=elem1
   END DO
   
END DO
 ! Among elem_neighbors there can be negative numbers. These correspond to 
 ! boundary elements for which neighbours are absent. However, an element 
 ! should have at least two valid neighbors

 ! Test that there are at least two neighbors at the surface:

DO elem=1,mesh%elem2D
   elem1=0
   DO j=1,3
   if(mesh%elem_neighbors(j,elem)>0) elem1=elem1+1
   END DO
   if (elem1<2) then
    write(*,*) 'find_elem_neighbors_ini:Insufficient number of neighbors ',elem
    write(*,*) 'find_elem_neighbors_ini:Elem neighbors ',mesh%elem_neighbors(:,elem)
    write(*,*) '____________________________________________________________________'
    write(*,*) ' ERROR: The mesh you want to partitioning contains triangles that'
    write(*,*) '        have just one neighbor, this was OK for FESOM1.4 but not'
    write(*,*) '        for FESOM2.0.                                           '
    write(*,*) '                                                                '
    write(*,*) '          #########################################             '
    write(*,*) '          ################### o ###################             '
    write(*,*) '          ################# ./|\. #################             '
    write(*,*) '              Land    ### ./|||||\. ###    Land                 '
    write(*,*) '          ############## /|||||||||\ ##############             '
    write(*,*) '          --o-----------o-----------o-----------o--             '
    write(*,*) '          ./ \.       ./ \.       ./ \.       ./ \.             '
    write(*,*) '               \.   ./     \.   ./     \.   ./                  '
    write(*,*) '                 \ /         \ /         \ /                    '
    write(*,*) '            ------o-----------o-----------o-------              '
    write(*,*) '                ./ \.       ./ \.       ./ \.                   '
    write(*,*) '                                                                '
    write(*,*) '        Take a programm of your choice (Python, Matlab ...) and '
    write(*,*) '        eliminate these triangles and the corresponding         '
    write(*,*) '        unconnected vertice and try to re-partitioning again    '
    write(*,*) '____________________________________________________________________'
    STOP
   end if
END DO    

 ! The rotation sense: corresponds to edges, and edges correspond 
 ! to nodes

 ! =============
 ! To facilitate computations the neibourhood
 ! information is assembled
 ! =============	 
 allocate(mesh%nod_in_elem2D_num(mesh%nod2D))
 mesh%nod_in_elem2D_num=0
 do n=1,mesh%elem2D
    elnodes=mesh%elem2D_nodes(1:3,n)
    mesh%nod_in_elem2D_num(elnodes)=mesh%nod_in_elem2D_num(elnodes)+1
 end do
 allocate(mesh%nod_in_elem2D(maxval(mesh%nod_in_elem2D_num),mesh%nod2D))
 mesh%nod_in_elem2D=0
 
 mesh%nod_in_elem2D_num=0
 do n=1,mesh%elem2D   
    elnodes=mesh%elem2D_nodes(1:3,n)
    mesh%nod_in_elem2D_num(elnodes)=mesh%nod_in_elem2D_num(elnodes)+1
    do j=1, 3
       mesh%nod_in_elem2D(mesh%nod_in_elem2D_num(elnodes(j)),elnodes(j))=n
    end do
 end do  
END SUBROUTINE find_elem_neighbors_ini
!===================================================================
! Stiffness matrix for the elevation
subroutine stiff_mat_ini(mesh)
  
  !
  implicit none
  integer                :: i, j, n, q, el, elem_nodes_max, nod(4)
  integer, allocatable   :: num_ne(:), ne(:,:)
  !
  type(t_mesh), intent(inout), target :: mesh

  mesh%ssh_stiff%dim = mesh%nod2D   
  allocate(mesh%ssh_stiff%rowptr(mesh%nod2D+1))

  allocate(num_ne(mesh%nod2D), ne(MAX_ADJACENT,mesh%nod2D))
  num_ne(:)           = 0
  ne(:,:)             = -1
 
  ! Check the maximum number of nodes in an element (FESOM triangular meshes = 3, Hybrid meshes = 4)
  elem_nodes_max = size(mesh%elem2D_nodes, 1)

  ! Fill node adjacency info
  ! all nodes in an element are adjacent in the sense of being halo nodes
  ! (also the opposite nodes of a quad: no edge, but the indirect connection
  ! should be taken into account by metis domain decomposition)
  do el=1,mesh%elem2D
     nod(1:elem_nodes_max) = mesh%elem2D_nodes(1:elem_nodes_max,el)  ! Fortran-numbering
     q = elem_nodes_max
     if (nod(1) == nod(elem_nodes_max)) q = q-1  ! triangle
    
     do i=2,q
        do j=1,i-1
           if (all(ne(:,nod(i)) /= nod(j))) then
              num_ne(nod(i)) = num_ne(nod(i)) + 1
              num_ne(nod(j)) = num_ne(nod(j)) + 1

              if (max(num_ne(nod(i)), num_ne(nod(j))) > MAX_ADJACENT ) then
                 print *,'Recompile with larger value for MAX_ADJACENT.'
                 stop
              else
                 ne(num_ne(nod(i)), nod(i)) = nod(j)
                 ne(num_ne(nod(j)), nod(j)) = nod(i)
              endif
           endif
        end do
     end do
  end do

! copy adjacency matrix to CSR-format
  mesh%ssh_stiff%rowptr(1) = 1
  do n=1,mesh%nod2D
     mesh%ssh_stiff%rowptr(n+1) = mesh%ssh_stiff%rowptr(n)+num_ne(n)
  end do

  allocate(mesh%ssh_stiff%colind(mesh%ssh_stiff%rowptr(mesh%nod2D+1)-1))  

  !required after the allocation, otherwise the pointer remains undefined
  do n=1,mesh%nod2D
     mesh%ssh_stiff%colind(mesh%ssh_stiff%rowptr(n):mesh%ssh_stiff%rowptr(n+1)-1) = ne(1:num_ne(n),n)
  end do

  deallocate(num_ne, ne)

end subroutine stiff_mat_ini

!===================================================================
! Setup of communication arrays
subroutine communication_ini(partit, mesh)
  implicit none

  integer        :: n
  character*10   :: npes_string
  character(MAX_PATH)  :: dist_mesh_dir
  LOGICAL        :: L_EXISTS
  type(t_mesh),   intent(inout), target :: mesh
  type(t_partit), intent(inout), target :: partit

  print *, achar(27)//'[1m'  //'____________________________________________________________'//achar(27)//'[0m'
  print *, achar(27)//'[7;1m' //' -->: compute communication arrays                          '//achar(27)//'[0m'
 
  ! Create the distributed mesh subdirectory
  write(npes_string,"(I10)") partit%npes
  dist_mesh_dir=trim(meshpath)//'dist_'//trim(ADJUSTL(npes_string))//'/'
  INQUIRE(file=trim(dist_mesh_dir), EXIST=L_EXISTS)
  if (.not. L_EXISTS) call system('mkdir '//trim(dist_mesh_dir))

#ifdef OMP_MAX_THREADS
!$OMP PARALLEL NUM_THREADS(OMP_MAX_THREADS)
  if (OMP_GET_THREAD_NUM() == 0) then
     write(*,*) 'Setting up communication arrays using ', OMP_GET_NUM_THREADS(), ' threads'
  endif
#else
!$OMP PARALLEL NUM_THREADS(1)
  write(*,*) 'Setting up communication arrays using 1 thread (serially)'
#endif
  
!$OMP DO
  do n = 0, partit%npes-1
     partit%mype = n ! mype is threadprivate and must not be iterator
     call communication_nodn(partit, mesh)
     call communication_elemn(partit, mesh)
     call save_dist_mesh(partit, mesh)         ! Write out communication file com_infoxxxxx.out
  end do
!$OMP END DO
!$OMP END PARALLEL

  deallocate(mesh%elem_neighbors)
  deallocate(mesh%elem_edges)
  deallocate(partit%part)
  write(*,*) 'Communication arrays have been set up'   
end subroutine communication_ini
!=================================================================
subroutine set_par_support_ini(partit, mesh)
  implicit none

  integer         :: n, j, k, nini, nend, ierr
  integer(idx_t)  :: np(10)
  type(t_mesh),   intent(inout), target :: mesh
  type(t_partit), intent(inout), target :: partit
  interface 
     subroutine do_partit(n,ptr,adj,wgt,np,part) bind(C)
       use iso_c_binding, only: idx_t=>C_INT32_T
       integer(idx_t), intent(in)  :: n, ptr(*), adj(*), wgt(*), np(*)
       integer(idx_t), intent(out) :: part(*)
     end subroutine do_partit
  end interface


  if (partit%mype==0) then
        print *, achar(27)//'[1m'  //'____________________________________________________________'//achar(27)//'[0m'
        print *, achar(27)//'[7;1m' //' -->: compute partitioning                                  '//achar(27)//'[0m'
  end if 
  
  ! Construct partitioning vector
  if (n_levels<1 .OR. n_levels>10) then
     print *,'Number of hierarchic partition levels is out of range [1-10]! Aborting...'
     stop
  end if

  np(:) = n_part(:)               ! Number of partitions on each hierarchy level
  if (n_part(1) == 0) then        ! Backward compatibility case: Take the number of 
     np(1) = partit%npes                 ! partitions from the number of MPI processes
     n_levels = 1
  end if
  if (n_levels < 10) then         ! 0 is an indicator of the last hierarchy level
     np(n_levels+1) = 0 
  end if

  allocate(partit%part(mesh%nod2D))
  partit%part=0

  partit%npes = PRODUCT(np(1:n_levels))
  if(partit%npes<1) then
     print *,'Total number of parallel partitions is less than one! Aborting...'
     stop
  end if
  
  write(*,*) 'Calling partit for npes=', np
  call do_partit(mesh%ssh_stiff%dim, mesh%ssh_stiff%rowptr, mesh%ssh_stiff%colind, &
       mesh%nlevels_nod2D, np, partit%part)

write(*,*) 'DONE'
write(*,*) size(partit%part)
  call check_partitioning(partit, mesh)

  write(*,*) 'Partitioning is done.'

! The stiffness matrix is no longer needed. 
  deallocate(mesh%ssh_stiff%rowptr)
  deallocate(mesh%ssh_stiff%colind)
        
  !NR No longer needed - last use was as weight for partitioning
  deallocate(mesh%nlevels_nod2D)
end subroutine set_par_support_ini
!=======================================================================
subroutine check_partitioning(partit, mesh)

  ! In general, METIS 5 has several advantages compared to METIS 4, e.g.,
  !   * neighbouring tasks get neighbouring partitions (important for multicore computers!)
  !   * lower maximum of weights per partition (better load balancing)
  !   * lower memory demand
  !
  ! BUT: there might be outliers, single nodes connected to their partition by
  !      only one edge or even completely isolated. This spoils everything :-(
  !
  ! This routine checks for isolated nodes and moves them to an adjacent partition,
  ! trying not to spoil the load balance.

  type(t_partit), intent(inout), target :: partit
  type(t_mesh),   intent(inout), target :: mesh
  integer :: i, j, k, n, n_iso, n_iter, is, ie, kmax, np, cnt
  integer :: nod_per_partition(2,0:partit%npes-1)
  integer :: max_nod_per_part(2), min_nod_per_part(2)
  integer :: average_nod_per_part(2), node_neighb_part(100)
  logical :: already_counted, found_part

  integer :: max_adjacent_nodes
  integer, allocatable :: ne_part(:), ne_part_num(:), ne_part_load(:,:)

  if (partit%mype==0) then
        print *, achar(27)//'[1m'  //'____________________________________________________________'//achar(27)//'[0m'
        print *, achar(27)//'[7;1m' //' -->: check partitioning                                    '//achar(27)//'[0m'
  end if

  ! Check load balancing
  do i=0, partit%npes-1
     nod_per_partition(1,i) = count(partit%part(:) == i)
     nod_per_partition(2,i) = sum(mesh%nlevels_nod2D, partit%part(:) == i)
  enddo

  min_nod_per_part(1) = minval( nod_per_partition(1,:))
  min_nod_per_part(2) = minval( nod_per_partition(2,:))

  max_nod_per_part(1) = maxval( nod_per_partition(1,:))
  max_nod_per_part(2) = maxval( nod_per_partition(2,:))

  average_nod_per_part(1) = mesh%nod2D / partit%npes
  average_nod_per_part(2) = sum(mesh%nlevels_nod2D(:)) / partit%npes

  ! Now check for isolated nodes (connect by one or even no edge to other
  ! nodes of its partition) and repair, if possible

  max_adjacent_nodes = maxval(mesh%ssh_stiff%rowptr(2:mesh%nod2D+1) - mesh%ssh_stiff%rowptr(1:mesh%nod2D))
  allocate(ne_part(max_adjacent_nodes), ne_part_num(max_adjacent_nodes), &
       ne_part_load(2,max_adjacent_nodes))

  print *,' '
  print *,'Check for isolated nodes ========'
  n_iso = 0
  checkloop: do n=1,mesh%nod2D
     is = mesh%ssh_stiff%rowptr(n)
     ie = mesh%ssh_stiff%rowptr(n+1) -1
     cnt = ie-is+1
     
     node_neighb_part(1:cnt) = partit%part(mesh%ssh_stiff%colind(is:ie))
     if (count(node_neighb_part(1:cnt) == partit%part(n)) <= 1) then

        n_iso = n_iso+1
        print *,'Isolated node',n, 'in partition', partit%part(n)
        print *,'Neighbouring nodes are in partitions',  node_neighb_part(1:cnt)

        ! count the adjacent nodes of the other PEs

        np=1
        ne_part(1) = node_neighb_part(1)
        ne_part_num(1) = 1
        ne_part_load(1,1) = nod_per_partition(1,ne_part(1)) + 1
        ne_part_load(2,1) = nod_per_partition(2,ne_part(1)) + mesh%nlevels_nod2D(n)

        do i=1,cnt
           if (node_neighb_part(i)==partit%part(n)) cycle
           already_counted = .false.
           do k=1,np
              if (node_neighb_part(i) == ne_part(k)) then
                 ne_part_num(k) = ne_part_num(k) + 1
                 already_counted = .true.
                 exit
              endif
           enddo
           if (.not. already_counted) then
              np = np+1
              ne_part(np) = node_neighb_part(i)
              ne_part_num(np) = 1
              ne_part_load(1,np) = nod_per_partition(1,ne_part(np)) + 1
              ne_part_load(2,np) = nod_per_partition(2,ne_part(np)) + mesh%nlevels_nod2D(n)
           endif
        enddo

        ! Now, check for two things: The load balance, and if 
        ! there is more than one node of that partition.
        ! Otherwise, it would become isolated again.

        ! Best choice would be the partition with most adjacent nodes (edgecut!)
        ! Choose, if it does not decrease the load balance. 
        !        (There might be two partitions with the same number of adjacent
        !         nodes. Don't care about this here)

        kmax = maxloc(ne_part_num(1:np),1)

        if (ne_part_num(kmax) <= 1) then
           ! No chance - this is probably a boundary
           ! node that has only two neighbors.
           cycle checkloop
        endif

        if  (ne_part_load(1,kmax) <= max_nod_per_part(1) .and. &
             ne_part_load(2,kmax) <= max_nod_per_part(2) ) then
           k = kmax
        else
           ! Don't make it too compicated. Reject partitions that have only one
           ! adjacent node. Take the next not violating the load balance.
           found_part = .false.
           do k=1,np
              if (ne_part_num(k)==1 .or. k==kmax) cycle

              if  (ne_part_load(1,k) <= max_nod_per_part(1) .and. &
                   ne_part_load(2,k) <= max_nod_per_part(2) ) then

                 found_part = .true.
                 exit
              endif
           enddo

           if (.not. found_part) then
              ! Ok, don't think to much. Simply go for minimized edge cut.
              k = kmax
           endif
        endif

        ! Adjust the load balancing

        nod_per_partition(1,ne_part(k)) = nod_per_partition(1,ne_part(k)) + 1 
        nod_per_partition(2,ne_part(k)) = nod_per_partition(2,ne_part(k)) + mesh%nlevels_nod2D(n)
        nod_per_partition(1,partit%part(n))    = nod_per_partition(1,partit%part(n)) - 1
        nod_per_partition(2,partit%part(n))    = nod_per_partition(2,partit%part(n)) - mesh%nlevels_nod2D(n)

        ! And, finally, move nod n to other partition        
        partit%part(n) = ne_part(k)
        print *,'Node',n,'is moved to part',partit%part(n)
     endif
  enddo checkloop

  deallocate(ne_part, ne_part_num, ne_part_load)

  print *,'=== LOAD BALANCING ==='
  print *,'2D nodes: min, aver, max per part',min_nod_per_part(1), &
       average_nod_per_part(1),max_nod_per_part(1)

  write(*,"('2D nodes: percent min, aver, max ',f8.3,'%, 100%, ',f8.3,'%')") &
       100.*real(min_nod_per_part(1)) / real(average_nod_per_part(1)), &
       100.*real(max_nod_per_part(1)) / real(average_nod_per_part(1))

  print *,'3D nodes: Min, aver, max per part',min_nod_per_part(2), &
       average_nod_per_part(2),max_nod_per_part(2)
  write(*,"('3D nodes: percent min, aver, max ',f8.3,'%, 100%, ',f8.3,'%')") &
       100.*real(min_nod_per_part(2)) / real(average_nod_per_part(2)), &
       100.*real(max_nod_per_part(2)) / real(average_nod_per_part(2))

end subroutine check_partitioning

end module mod_mesh_utils
