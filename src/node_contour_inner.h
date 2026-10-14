nn=0
do while (nn<mesh%nod_in_elem2D_num(n))
  if (elem==0) then
write(*,*) 'error diag:', nn, elem, el(1), el(2), bEdge_left, bEdge_right
  end if
   if (mesh%elem2D_nodes(1,elem)==n) then
edge_left=mesh%elem_edges(3,elem)
edge_right=mesh%elem_edges(2,elem)
   elseif (mesh%elem2D_nodes(2,elem)==n) then
edge_left=mesh%elem_edges(1,elem)
edge_right=mesh%elem_edges(3,elem)
   else
edge_left=mesh%elem_edges(2,elem)
edge_right=mesh%elem_edges(1,elem)
   end if
   nn=nn+1
   nedges(nn)=edge_left
   nelems(nn)=elem
   el=mesh%edge_tri(:,edge_right)
   if (el(1)==elem) then
elem=el(2)
   else
elem=el(1)
   end if
end do