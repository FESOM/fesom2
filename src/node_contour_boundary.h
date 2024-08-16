flag=1
nn=0
do while (flag==1)
   if (elem2D_nodes(1,elem)==n) then
      edge_left=elem_edges(3,elem)
      edge_right=elem_edges(2,elem)
   elseif (elem2D_nodes(2,elem)==n) then
      edge_left=elem_edges(1,elem)
      edge_right=elem_edges(3,elem)
   else
      edge_left=elem_edges(2,elem)
      edge_right=elem_edges(1,elem)
   end if
   nn=nn+1
   nedges(nn)=edge_left
   nelems(nn)=elem
   el=edge_tri(:,edge_right)
   if (el(2)>0) then
      if (el(1)==elem) then
  elem=el(2)
      else
  elem=el(1)
      end if
   else    !the last element
      nedges(nn+1)=edge_right
      flag=0
   end if
end do