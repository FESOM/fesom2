.. _spatial_discretization


Spatial discretization of equations (mostly in horizontal directions)
*********************************************************************

Finite-volume basics
====================

To obtain the finite-volume discretization, the governing equations are integrated over the selected control volumes. The flux divergence terms are then, by virtue of Gauss theorem, transformed to the net fluxes leaving the control volumes. Basic discrete fields such as horizontal velocities and scalar tracers are considered to be  mean over the respective finite volumes. We already described layered equations, i.e., the basic part of vertical discretization, where fluxes through top and bottom layer surfaces were appearing. In this section we discuss only the horizontal part.  The discrete fields are understood as

.. math::
   A_ch_c{\bf u}_c=\int_c {\bf u}dV,

and similarly for the temperature and other scalars,

.. math::
   A_{kv}h_{kv}T_{kv}=\int_{kv}TdV.

Here :math:`A_c` and :math:`A_{kv}` are the horizontal areas of cells and scalar prisms. The scalar areas vary with depth, hence the index :math:`k` in :math:`A_{kv}` in the formula above (the index :math:`k` will be suppressed in some cases). For layer :math:`k`, :math:`A_{kv}` is the area of the prism :math:`kv` including its top face. The area of bottom face is :math:`A_{(k+1)v}` and may differ from that of the top one if the bottom is encountered. To be consistent in spherical geometry, we us

.. math::
   A_{kv}=\sum_{c\in\overline C(v)}A_c/3,

where :math:`\overline C(v)` is the set of wet prisms containing :math:`v` in layer :math:`k`.

Since the horizontal velocity is at centroids, its cell-mean value :math:`{\bf u}_c` can be identified with the value of the field :math:`{\bf u}` at the centroid of cell :math:`c` with the second order of spatial accuracy. For scalar quantities a similar rule is valid only on uniform meshes, but even in this case it is violated in the vicinity of boundaries or topography. This has some implications for the accuracy of transport operators.

Horizontal operators
====================

Scalar gradient
---------------

Scalar gradient takes vertex values of a scalar field and returns the gradient at the cell center:

.. math::
   A_c(\nabla p)_c=\int_c\nabla pdS=\sum_{e\in E(c)}l_e{\bf n}_e\sum_{v\in V(e)}p_v/2,

where :math:`{\bf n}_e` is the outer normal to cell :math:`c`. Clearly :math:`l_e{\bf n}_e=-{\bf k}\times{\bf l}_e` if :math:`c` is the first (left) cell of :math:`c(e)`. This procedure introduces :math:`{\bf G}_{cv}=(G^x_{cv},G^y_{cv})` with the :math:`x`- and :math:`y` component matrices :math:`G^x_{cv}` and :math:`G^y_{cv}`. They have three non-zero entries for each cell (triangle) which are stored in the order of triangle vertices, first for :math:`x` and then for :math:`y` component in the array ``gradient_sca(1:6,1:myDim\_elem2D)`` in the code. In contrast to FESOM1.4, where similar arrays are stored for each tetrahedron (and for 4 vertices and 3 directions), here only surface cells are involved.


Vector gradient
---------------

Vector gradient takes the values of velocity components at neighbor triangles of a given triangle and returns their derivatives at this triangle. They are computed through the least squares fit based on the velocities on neighboring cells sharing edges with cell :math:`c`. Their set is :math:`N(c)`. The derivatives :math:`(\alpha_x, \alpha_y)` of the velocity component :math:`u` are found by minimizing

.. math::
   \mathcal{L}=\sum_{n(c)}(u_{c}-u_{n}+(\alpha_x, \alpha_y)\cdot{\bf r}_{cn})^2={\rm min}.

Here :math:`{\bf r}_{cn}=(x_{cn}, y_{cn})` is the vector connecting the center of :math:`c` to that of its neighbor :math:`n`. The solution of the minimization problem can be represented as two matrices :math:`g_{cn}^x` and :math:`g_{cn}^y`, acting on velocity differences :math:`u_n-u_c` and returning the derivatives. Computations for :math:`v`-component result in the same matrices. The explicit expressions for matrix entries are:

.. math::
   g_{cn}^x=(x_{cn}Y^2-y_{cn}XY)/d, {} \\ {} g_{cn}^y=(y_{cn}X^2-x_{cn}XY)/d. {}

Here :math:`d=X^2Y^2-(XY)^2`, :math:`X^2=\sum_{n\in N(c)} x_{cn}^2`,
:math:`Y^2=\sum_{n(c)} y_{cn}^2` and :math:`XY=\sum_{n\in N(c)} x_{cn}y_{cn}`. The matrix entries stored in the array ``gradient_vec(1:6,1:myDim_elem2D)`` in the order the neighbor elements are listed, first for :math:`x` and then for :math:`y` component. They are computed only once.

On the cells touching the lateral walls or bottom topography we use ghost cells (mirror reflections with respect to boundary edge). Their velocities are computed either as :math:`{\bf u}_{n}=-{\bf u}_{c}` or :math:`{\bf u}_{n}={\bf u}_{c}-2({\bf u}_{c}\cdot{\bf n}_{nc}){\bf n}_{nc}` for the no-slip or free-slip cases respectively. Here :math:`n` is the index of the ghost cell, and :math:`{\bf n}_{nc}` is the vector of unit normal to the edge between cells :math:`c` and :math:`n`. Note that filing ghost cells takes additional time, but allows using matrices :math:`g_{cn}^x` and :math:`g_{cn}^y` related to the surface cells only. Otherwise separate matrices will be needed for each layer. Note also that ghost cells are insufficient to implement the free-slip condition. In addition, the tangent component of viscous stress should be eliminated directly.


We stress that matrices :math:`g_{cn}^x` and :math:`g_{cn}^y` return derivatives of velocity components, and not the components of the tensor of velocity derivatives. The latter includes additional metric terms that need to be taken into account separately.

Flux divergence
---------------

Flux divergence takes fluxes defined at boundaries of scalar control volume (and in the simplest case just on cells) and returns their divergence on scalar control volumes:

.. math::
   A_{kv}(\nabla\cdot {\bf F})_vh_v=\sum_{e\in E(v)}\sum_{c\in C(e)}{\bf F}_ch_c\cdot {\bf n}_{ec}d_{ec},

where :math:`{\bf n}_{ec}` is the outer normal to control volume :math:`v`. Clearly, if :math:`v` is the first vertex in the set :math:`v(e)`, :math:`{\bf n}_{ec}d_{ec}=-{\bf k}\times{\bf d}_{ec}` if :math:`c` is the first in the set :math:`c(e)` (signs are changed accordingly in other cases). While these rules may sound difficult to memorize, in practice computations are done in a cycle over edges, in which case signs are obvious.

In contrast to the scalar gradient operator, the operator of divergence depends on the layer (because of bottom topography), which is one of the reasons why it is not stored in advance. Besides, except for the simplest cases such as :math:`\mathbf{F}=\mathbf{u}`, the fluxes :math:`{\bf F}` involve estimates of the scalar quantity being transported. Computing these estimates requires a cycle over edges in any case, so there would be no economy even if the matrices of the divergence operator were introduced.


Velocity curl
-------------

Velocity curl takes velocities at cells and returns the relative vorticity at vertices using the circulation theorem:

.. math::
   A_{kv}\int_v(\nabla\times{\bf u})\cdot {\bf k}dS=\sum_{e\in E(v)}\sum_{c\in C(e)}{\bf u}_c\cdot{\bf t}_{ec}d_{ec},

where :math:`{\bf t}_{ec}` is the unit vector along :math:`{\bf d}_{ec}` oriented so as to make an anticlockwise turn around vertex :math:`v`. If :math:`v` is the first in the set :math:`v(e)` and :math:`c` is the first in the set :math:`c(e)`,  :math:`{\bf t}_{ec}d_{ec}={\bf d}_{ec}`. This operator also depends on the layer and is not stored.

Mimetic properties
------------------

It can be verified that the operators introduced above are mimetic, which means that they satisfy the properties of their continuous analogs. For example, the scalar gradient and divergence are negative adjoint of each other in the sense of scalar products that provide energy norm, and the curl operator applied to the scalar gradient operator gives identically zero. The latter property implies that the curl of discrete pressure gradient is identically zero, which is a prerequisite of a PV conserving discretization.

Momentum advection
==================

FESOM2.0 has three options for momentum advection. Two of them are different implementations of the flux form and the third one relies on the vector invariant form. In spherical geometry the flux form acquires an additional term :math:`M{\bf k}\times{\bf u}` on the lhs, where :math:`M=u\tan\lambda/r_E` is the metric frequency, with :math:`\lambda` the latitude and :math:`r_E` the Earth radius. All the options are based on the understanding that the cell-vertex discretization has an excessive number of velocity degrees of freedom on triangular meshes. The implementation of momentum advection must contain certain averaging in order to suppress the appearance of grid-scale noise.

Flux form with velocities averaged to vertices
----------------------------------------------

The momentum flux is computed based on vertex velocities. We compute vertex velocities by averaging from cell to vertex locations

.. math::
   A_{kv}{\bf u}_{kv}h_{kv}=\sum_{c\in\overline C(v)}{\bf u}_{kc}h_{kc}A_c/3,

and use them to compute the divergence of horizontal momentum flux:

.. math::
   A_c(\nabla\cdot(h{\bf u u}))_c=\sum_{e\in E(c)}l_e(\sum_{v\in V(e)}{\bf n}_e\cdot{\bf u}_vh_v)(\sum_{v\in V(e)}{\bf u}_v/4).

Here :math:`{\bf n}_e` is the external normal and :math:`l_e{\bf n}_e=-{\bf k}\times{\bf l}_e` if :math:`c` is the first one in the set :math:`c(e)`. Since the horizontal velocity appears as the product with the thickness, the expressions here can be rewritten in terms of transports :math:`{\bf U}={\bf u}h^*`.

The fluxes through the top and bottom faces are computed with :math:`w_c=\sum_{v\in V(c)}w_v/3` using either the second or fourth order centered, or high-order upwind algorithms.

Flux form relying on scalar control volumes
-------------------------------------------

Instead of using vector control volumes, we assemble the flux divergence on the scalar control volumes and then average the result from the vertices to the cells. Here the same idea of averaging as in the previous case is applied to the momentum advection term instead of velocities. For the horizontal part,

.. math::
   A_{v}(\nabla\cdot(h{\bf u u}))_v=\sum_{e(v)}\sum_{c\in C(e)}{\bf u}_ch_c\cdot {\bf n}_{ec}{\bf u}_cd_{ec},

with the same rule for the normals as in the computations of the divergence operator.
The contributions from the top and bottom faces of scalar control volume are obtained by summing the contributions from the cells:

.. math::
   A_{v}(w_v {\bf u}^t)=w_v\sum_{c\in \overline C(v)}{\bf u}^t_cA_c/3

for the top surface, and similarly for the bottom one. The estimate of :math:`{\bf u}^t` can be either centered or upwind as above.

This option of momentum advection is special in the sense that the continuity is treated here in the same way as for the scalar quantities.


Vector-invariant form
=====================

The relative vorticity in the cell-vertex discretization is defined on vertices, and so should be the Coriolis parameter. We use the following representation

.. math::
   ((\omega+f){\bf k}\times{\bf u})_c=\sum_{v\in V(c)}(\omega+f)_v{\bf k}\times{\bf u}_c/3.

The representation with the thicknesses,

.. math::
   ((\omega+f){\bf k}\times{\bf u})_c=\sum_{v\in V(c)}\frac{\omega_v+f_v}{3h_v}{\bf k}\times{\bf u}_ch_c

is reserved for future.
The gradient of kinetic energy should be computed in the same way as the pressure gradient, which necessitates computations of :math:`{\bf u}^2` at vertices. This is done as

.. math::
   A_{v}{\bf u}^2_v=\sum_{c\in \overline C(v)}A_c{\bf u}^2_c/3.

The vertical part follows :eq:`eq_mom_vei`,

.. math::
   (w\partial_z{\bf u})^t_{c}=2({\bf u}_{(k-1)c}-{\bf u}_{kc})/(h_{(k-1)c}+h_{kc})\sum_{v(c)}w_{kv}/3

for the top surface and similarly for the bottom. Note that the contributions from the curl of horizontal velocity, the gradient of kinetic energy and the vertical part involve the same stencil of horizontal velocities.

The three options above behave similarly in simple tests on triangular meshes, but their effect on flow-topography interactions or eddy dynamics still needs to be studied.

Horizontal viscosity operators
==============================

Because the cell placement of velocities, there are more velocity degrees of freedom than needed for vertex scalars. This leads to spurious grid-scale oscillations, and the task of horizontal viscosity operator is to eliminate them.

The derivatives of horizontal velocity can be estimated at cell locations and then averaged to edges of triangles enabling computation of the viscous stress tensor :math:`\sigma_{ij}=\nu_hs_{ij}`, :math:`s_{ij}=(\partial_iu_j+\partial_ju_i)/2`, where the indices :math:`i,j` imply the horizontal directions, :math:`s_{ij}` is the strain rate tensor and :math:`\nu_h` is the harmonic horizontal viscosity coefficient. Its divergence will give the viscous force. This would be the standard way of introducing viscosity. It turns out that for cell velocities such a viscous force is insensitive to the difference in the nearest velocities (see :cite:`DanilovKutsenko2019`). To eliminate grid-scale fluctuations FESOM is bound to use other discretizations.

FESOM relies on a simplified stresses :math:`\sigma_{ij}=\nu_h\partial_iu_j`. As discussed by :cite:`Griffiesbook`, their divergence still ensures energy dissipation, but is nonzero for solid-body rotations if :math:`\nu_h` is variable. In spite of this drawback, using the simplified form is much more convenient for numerical reasons: since the divergence of stresses reduces to fluxes over vertical faces of triangular prisms, only contraction of stresses with normal vector appears, i.e., :math:`\nu_hn_i\partial_iu_j`, which is the derivative in the direction of :math:`{\bf n}`.

'Canonical' operators
=====================

'Canonical' harmonic viscosity operator is written as

.. math::
   (D_{u}\mathbf{u})_cA_ch_c=\int_c\nabla\cdot(\nu_h\nabla\mathbf{u})h_cdS_c=\sum_{e\in E(c)}l_eh_e\mathbf{n}_e\cdot(\nu_h\nabla\mathbf{u})_e.
   :label: eq_viscL1

The operator :math:`D_{uh}` appearing in :eq:`eq_mom_fl` will only differ by the absence of :math:`h_c` on the lhs, and will not be written separately.

Let :math:`n` be the cell sharing edge :math:`e` with cell :math:`c`. We formally write for the edge normal vector :math:`{\bf n}_e={\bf r}_{cn}/|{\bf r}_{cn}|+({\bf n}-{\bf r}_{cn}/|{\bf r}_{cn}|)`, where :math:`{\bf r}_{cn}={\bf d}_{en}-{\bf d}_{ec}` is the vector connecting the centroids of cells :math:`c` and :math:`n`. Then

.. math::
   \mathbf{n}_e\cdot(\nu_h\nabla\mathbf{u})_e=(\nu_h)_e\frac{\mathbf{u}_n-\mathbf{u}_c}{|{\bf r}_{cn}|}+({\bf n}-{\bf r}_{cn}/|{\bf r}_{cn}|)\cdot(\nu_h\nabla\mathbf{u})_e.
   :label: eq_viscL2

The first term on the rhs of :eq:`eq_viscL2` depends on the velocity difference across the edge, and this ensures that the viscous operator :eq:`eq_viscL1`-:eq:`eq_viscL2` will act to reduce the difference. The viscosity :math:`\nu_h` has to be estimated at edges.

The operations are implemented in two cycles. The first one estimates velocity derivatives at cells. The second one is over edges, and edge values of velocity derivatives are obtained by averaging between the cells sharing the edge.

The 'canonical' biharmonic viscous operator is obtained by computing first the field

.. math::
   \mathbf{B}_c=\frac{-1}{A_ch_c}\sum_{e\in E(c)}l_eh_e(\nu_h)_e\left(\frac{\mathbf{u}_n-\mathbf{u}_c}{|{\bf r}_{cn}|}+({\bf n}-{\bf r}_{cn}/|{\bf r}_{cn}|)\cdot(\nabla\mathbf{u})_e\right)

and then applying the operator :math:`D_{u}` :eq:`eq_viscL1`-:eq:`eq_viscL2` to this field under the agreement that :math:`\nu_h=(\nu_{bh})^{1/2}`, where :math:`\nu_{bh}` is positive biharmonic viscosity coefficient.

FESOM relies on an alternative version of biharmonic operator, where the viscosity coefficient is taken at :math:`c` locations and the :math:`\mathbf{B}` field becomes

.. math::
   \mathbf{B}_c=\frac{-(\nu_{bh})_c}{A_ch_c}\sum_{e\in E(c)}l_eh_e\left(\frac{\mathbf{u}_n-\mathbf{u}_c}{|{\bf r}_{cn}|}+({\bf n}-{\bf r}_{cn}/|{\bf r}_{cn}|)\cdot(\nabla\mathbf{u})_e\right).
   :label: eq_viscB1

In this case :math:`D_{u}` with :math:`\nu_h=1` is applied to :math:`\mathbf{B}`. The advantage of this form is that the combination :math:`(\nu_{bh})_c/A_c` in :eq:`eq_viscB1` has the dimension of harmonic viscosity :math:`\nu_h` and thus can be replaced by :math:`(\nu_h)_c`, avoiding the need to specify :math:`(\nu_{bh})_c`. Even if viscosity is flow-dependent, one needs a single routine to compute :math:`\nu_h`.

Selection of viscosity coefficient
==================================

:cite:`FoxKemperMenemenlis2008` review common recipes for the viscosity coefficient and provide necessary references. Below we list the options available in FESOM. They can be used with both harmonic and biharmonic operators as explained above.


- **Simple viscosity**
  A frequent practice in ocean modeling community is that :math:`\nu_h` is selected as :math:`\nu_h=Vl`, where :math:`V` is a velocity scale (commonly about 1 cm/s) and :math:`l` is the cell size. In FESOM this rule is replaced by :math:`(\nu_h)_c=VA_c^{1/2}`, and the edge values are obtained by averaging of values at neighboring cells. Note that FESOM choice  automatically provides scaling.

- **Smagorinsky viscosity**
  In this case :math:`\nu_h` is taken as

  .. math::
     (\nu_h)_c=C_{Smag}(A_c/\pi^2)(s_{ij}^2)^{1/2},

  where :math:`C_{Smag}` is a dimensionless factor about one, and :math:`A_c/\pi^2` is assumed to be an estimate of maximum wavenumber squared. In reality :math:`C_{Smag}` is a tunable parameter available through FESOM namelists. The velocity derivatives are needed at :math:`c` locations to compute :math:`s_{ij}^2`, see the section on velocity gradients.

- **Leith viscosity and its modification**
  In this case

  .. math::
     (\nu_h)_c=C_{Leith}(A_c/\pi^2)^{3/2}(|\nabla \omega|+C_{div}|\nabla\nabla\cdot\mathbf{u}|).

  Here :math:`C_{Leith}` and :math:`C_{div}` are dimensionless factors of order one. The term with :math:`C_{div}` was added by :cite:`FoxKemperMenemenlis2008`. The entire construction is called the modified Leith viscosity. Here vorticity and divergence should be computed first at vertices. Then their scalar gradients are estimated, giving values at cells.

Note that the Smagorinsky and (modified) Leith viscosities rely on additional computations which take time. Viscosities are computed at :math:`c` points and averaged to edges whenever needed.


Numerical viscosities (viscosity filters)
=========================================

The `canonical` viscosity, especially its biharmonic version, proves to be costly. Its expensive part involving the computation of velocity derivatives at cells is only needed on deformed meshes, but even there it contributes little to penalizing differences between the nearest velocities. This leads to the idea of simplified numerical operators based only on the nearest neighbors. We refer to them as viscosity filters. They resemble `canonic` operators, but deviate from them on irregular meshes.


Quasi-harmonic viscosity filter
-------------------------------

An analog of harmonic viscosity operator in FESOM is introduced as

.. math::
  A_ch_c(D_{u}\mathbf{u})_c=\sum_{n\in N(c)}({\bf u}_n-{\bf u}_c) h_{nc}\frac{l_{nc}}{|{\bf r}_{nc}|}\frac{\nu_n+\nu_c}{2}.

Here :math:`l_{nc}` is the length of the edge between cells :math:`n` and :math:`c`, :math:`h_{nc}` the layer thickness interpolated to the edge and :math:`{\bf r}_{nc}` is a vector connecting the cell centroids (we use :math:`nc` to identify edges here). If mesh cells (triangles) are equilateral, :math:`({\bf u}_n-{\bf u}_c)/|{\bf r}_{nc}|` is the velocity gradient in the direction of the normal to the edge between :math:`n` and :math:`c`. The expression above is then the sum of viscous fluxes leaving cell :math:`c`, which gives a harmonic viscosity operator :math:`D_u\mathbf{u}=\nabla\nu\nabla{\bf u}`. This interpretation fails on general meshes, but the formula above for viscous operator
provides the most efficient penalty of velocity difference between the nearest neighbors.

The appearance of :math:`A_ch_c` on the left hand side is critical and is needed for conservation on general meshes (see below). Since :math:`l_{nc}/|{\bf r}_{nc}` and the mean viscosity are related to the edge between :math:`n` and :math:`c`, they can be incorporated in generalized viscosity :math:`\nu_{nc}` associated to edges (hence :math:`\nu_{nc}=\nu_{cn}`), which gives

.. math::
  h_cA_c(D_u\mathbf{u})_c=\sum_{n\in N(c)}({\bf u}_n-{\bf u}_c)\nu_{nc}h_{nc}.
  :label: eq_viscQh


Since :math:`\sum_cA_ch_c(D_u\mathbf{u})_c=0` (the difference between :math:`n` and :math:`c` velocities appear with opposite signs in equations for :math:`n` and :math:`c`), the viscosity operator does not violate momentum conservation. Furthermore, taking :math:`\sum_c{\bf u}_cA_ch_c(V)_c`, we see that the contribution from :math:`n` and
:math:`c` comes twice, one time as :math:`h_{nc}\nu_{nc}{\bf u}_c({\bf u}_n-{\bf u}_c)` and the other time as  :math:`\nu_{nc}h_{nc}{\bf u}_n({\bf u}_c-{\bf u}_n)`, which sums to :math:`-\nu_{nc}({\bf u}_n-{\bf u}_c)^2`. This proves that area mean kinetic energy dissipation is non-positive.

Quasi-biharmonic viscosity filter
---------------------------------

It is only a bit more complicated with biharmonic implementation. First, we define :math:`b_{nc}=((\nu_{bh})_{nc}^b)^{1/2}`. We write

.. math::
   A_ch_c(D_u^b\mathbf{u})_c=-\sum_{n\in N(c)}({\bf L}_n-{\bf L}_c)b_{nc}h_{nc}^{1/2},

where

.. math::
   {\bf L}_c=\sum_{n\in N(c)}({\bf u}_n-{\bf u}_c)b_{nc}h_{nc}^{1/2}.

The momentum is conserved for the same reason as above. To see why the form above leads to kinetic energy dissipation we introduce matrix :math:`B` with the dimension of the number of cells, such that its entries are :math:`B_{nc}=b_{nc}h_{nc}^{1/2}` for those :math:`n` and :math:`c` that have common edges. It is a symmetric matrix. Then, we put the sum of row entries with the opposite sign at its diagonal. Obviously, in terms of this matrix :math:`{\bf L}_c=B_{cn}{\bf u}_n` (we assume summation over the repeating indices in matrix-vector multiplications). The energy dissipation is :math:`\sum_c{\bf u}_cA_ch_c(D_u^b\mathbf{u})_c` can be written as :math:`-{\bf u}_mB_{mc}B_{cn}{\bf u}_n=-{\bf L}^T{\bf L}\le 0`. Note that area appears here only at the last step of the procedure, which makes the two steps different.

Similarly to the `canonic` biharmonic viscosity, the split of biharmonic viscosity between the two harmonic steps can be avoided. In this case we write

.. math::
   {\bf L}_c=\sum_{n\in N(c)}({\bf u}_n-{\bf u}_c),\quad \mathbf{L}'_c=h_c(\nu_{bh})_c\mathbf{L}_c,
   :label: eq_viscQbh1

and

.. math::
   A_ch_c(D_u^b\mathbf{u})_c=-\sum_{n\in N(c)}({\bf L}'_n-{\bf L}'_c),
   :label: eq_viscQbh2

so that :math:`{\bf L}_c=B_{cn}{\bf u}_n` where :math:`B_{cn}=1` if :math:`c` and :math:`n` are neighbors, and :math:`B_{cc}=-\sum_{n\ne c}B_{nc}`. We introduce a diagonal matrix :math:`N` such that :math:`N_{cc}=(\nu_{bh})_ch_c`. The biharmonic operator can be written then as :math:`A_ch_c(D_u^b\mathbf{u})_c=-B_{cn}N_{nm}B_{mj}{\bf u}_j`, and making a dot product with :math:`{\bf u}_c` we will have a nonpositive energy dissipation. This alternative form is default in FESOM.

In all viscosity filter cases the main assembly is in the cycle(s) over edges.

Viscosities in viscosity filters
--------------------------------

Although viscosity filters can be using the standard viscosities specified above, there are less expensive and more convenient options. For the quasi-harmonic viscosity filter :eq:`eq_viscQh` the edge harmonic viscosity :math:`(\nu_h)_{nc}` can be estimated  in the same cycle where velocity differences are assembled as

.. math::
   (\nu_h)_{nc}=C_h|{\bf u}_n-{\bf u}_c|l_{nc},

where :math:`C_h` is a small factor determined experimentally (about 1/20). Since it depends on velocity differences, it is an analog of the Smagorinsky viscosity. For the quasi-biharmonic viscosity filter :eq:`eq_viscQbh1`,:eq:`eq_viscQbh2`, once :math:`\mathbf{L}_c` is computed, the cell biharmonic viscosity can be estimated as

.. math::
   (\nu_{bh})_c=C_{bh}A_c^{3/2}|\mathbf{L}_c|\quad ((\nu_{h})_c=C_{bh}A_c^{1/2}|\mathbf{L}_c|).

Since :math:`\mathbf{L}_c` contains double differences (on uniform meshes), it is an analog of the Leith viscosity. The advantage of these simplified expressions is that they use the already computed terms, leading to a more economical implementation.

Transport of scalar quantities. Horizontal part
===============================================

Horizontal and vertical fluxes are taken into account together, without operator splitting commonly used on structured meshes. However, the mesh is unstructured in the horizontal direction, the computation of the horizontal and vertical fluxes is different. They are therefore presented separately.

FESOM will provide the following horizontal advection schemes:

- A third-fourth order scheme based on gradient estimate (GE34),
- A third-fourth order scheme based on quadratic polynomial reconstruction (QR34)
- A compact scheme (C34)

GE34 is available at present, and the two other will be made available in 2020 (their description will be added).

All of them can be run with flux corrected transport (FCT) limiter. The FCT operates on full three-dimensional fluxes.

GE34
----

Consider edge :math:`e`. Let :math:`V(e)=(v_1,v_2)` and :math:`C(e)=(c_1, c_2)`. The advective flux of scalar quantity :math:`T` through the face of scalar volume associated to this edge is

.. math::
   F_e=T_e(-h_{c_1}{\bf d}_{ec_1}\times{\bf u}_{c_1}+h_{c_2}{\bf d}_{ec_2}\times{\bf u}_{c_2})\cdot{\bf k}=T_eQ_e.

The quantity :math:`Q_e` is the volume flux associated with edge :math:`e` which leaves the control volume :math:`v_1` and enters the control volume :math:`v_2`. To compute the mid-edge tracer estimate :math:`T_e`, for each edge :math:`e` the indices of the cells up or down this edge in the edge direction :math:`{\bf l}_e` are stored. Two estimates

.. math::
   T_e^+=T_{v_1}+(1/2){\bf l}_e(\nabla T)_e^+, \quad (\nabla T)_e^+=(2/3)(\nabla T)^c+(1/3)(\nabla T)^u,

and

.. math::
   T_e^-=T_{v_2}-(1/2){\bf l}_e(\nabla T)_e^-, \quad (\nabla T)_e^-=(2/3)(\nabla T)^c+(1/3)(\nabla T)^d

are computed. In these expressions, :math:`(\nabla T)^c{\bf l}_e=T_{v2}-T_{v1}` (:math:`c` for `centered`), while :math:`u` and :math:`d` imply the  gradients on up- and down-edge cells. They are computed and stored in a separate cycle by applying the scalar gradient operator.

The estimate

.. math::
   2T_eQ_e=(Q_e+|Q_e|)T_e^++(Q_e-|Q_e|)T_e^-

provides the standard third-order upwind method, and the estimate

.. math::
   2T_e=T_e^++T_e^-

provides the fourth-order centered method. The combination

.. math::
   2Q_eT_e=(T_e^++T_e^-)Q_e+(1-\gamma)(T_e^+-T_e^-)|Q_e|


takes the fourth-order part with the weight :math:`\gamma` and the third order part, with :math:`1-\gamma`. :math:`\gamma=0.75` will reducing the upwind dissipation by a factor of 4 compared to the third order method. The question whether this reduction is needed depends on applications.

The high order of the scheme above is only achieved on uniform meshes. Since :math:`T_e` is computed through linear reconstruction, the second order is warranted on general meshes.

The scheme requires preliminary computation of scalar gradients on cells before the main cycle over edges. An extended halo exchange is needed to make these gradients available during flux assembly.

Edges touching the topography lack either :math:`u` or :math:`d` cells. In this case the simplest choice is either to use the central estimate or the estimate based on the mean vertex gradient :math:`A_v(\nabla T)_v=\sum_{c\in \overline C(v)}A_c(\nabla T)_c/3`. The associated logistics is expensive and increases the CPU cost of the scheme.

Vertical advection of scalars
=============================

The approach advocated by the ROMS community is to avoid dissipation in vertical advection, i. e. use high-order centered schemes. Although vertical grids are not uniform, this does not affect the convergence rate very much if level spacing is smoothly varying. For discussion, see Treguier et al. (1998). What is affected, of course, is the magnitude of residual errors. Many high-order algorithms essentially rely on mesh uniformity (some error cancellation). Nevertheless, they are applied on nonuniform meshes. We describe the selection of schemes available now or soon in FESOM. One-dimensional approach is used.

Centered second order
---------------------

Let :math:`T` be the scalar to be advected. A bar will be used to denote the interface value: :math:`\overline{T}_{kv}` is the value reconstructed to the level :math:`k` separating layer :math:`k` above from layer :math:`k+1` below. The 'vertical' flux :math:`F` through the level :math:`k` for the centered scheme is

.. math::
   F_{kv}=w_{kv}A_{kv}(T_{(k-1)v}+T_{kv})/2\quad (\overline{T}_{kv}=(T_{(k-1)v}+T_{kv})/2).

Vertical GE34
-------------

We write

.. math::
   \overline{T}^+_{kv}=T_{kv}+(2G^c+G^u)h_{kv}/6, \quad w_{kv}>0,

.. math::
   \overline{T}^-_{kv}=T_{(k-1)v}-(2G^c+G^u)h_{(k-1)v}/6, \quad w_{kv}<0,

where :math:`G^c=(T_{(k-1)v}-T_{kv})/(Z_{(k-1)v}-Z_{kv})` is the central gradient estimate and :math:`G^u=(T_{kv}-T_{(k+1)v})/(Z_{kv}-Z_{(k+1)v})` for positive and :math:`G^u=(T_{(k-2)v}-T_{(k-1)v})/(Z_{(k-2)v}-Z_{(k-1)v})` for negative :math:`w_{kv}`. Note that our estimates of gradients are based on values that are mean over control volume. So the estimates themselves are not very accurate. It is the combination (of central and upwind) values that is accurate.

Using

.. math::
   2w_{kv}\overline{T}_{kv}=w_{kv}(\overline{T}^+_{kv}+\overline{T}^-_{kv})+(1-\gamma)|w_{kv}|(\overline{T}^+_{kv}-\overline{T}^-_{kv})

will give a fourth-order scheme on a uniform mesh if :math:`\gamma=1`. A blended third-fourth order scheme follows for :math:`0\le\gamma<1`.

Compact scheme (also the Parabolic Spline Method
================================================

We need scalar values at interfaces. An elegant way to find them is to use splines, requiring continuity of reconstruction and first derivatives at level locations. The result is

.. math::
   \overline{T}_{k+1}\frac{1}{h_k}+2\overline{T}_{k}\left(\frac{1}{h_k}+\frac{1}{h_{k-1}}\right)+\overline{T}_{k-1}\frac{1}{h_{k-1}}=3\left(T_k\frac{1}{h_k}+T_{k-1}\frac{1}{h_{k-1}}\right).

The boundary conditions are those of natural spline, i. e.,

.. math::
   2\overline{T}_{1}+\overline{T}_{2}=3T_1,\quad 2\overline{T}_{N+1}+\overline{T}_{N}=3T_N.

This method requires three-diagonal solve, which takes the same time as two vertical loops. The name `compact` reflects the fact that the equation above involves stencil of minimum size. It becomes the PSM method if used with semi-Lagrangian time stepping, as in PPM.

The result is more accurate than PPM (see further). It is of the fourth order as PPM on uniform grid, but has a smaller residual term. Those who learned piecewise linear finite elements may see some analogies in the reconstruction procedure. ROMS uses this method for vertical advection of both tracers and momentum.

Piecewise Parabolic Method
--------------------------

To be written


FCT
---

The FCT limiter in FESOM2 uses the first-order upwind method as the low-order monotonic method and a combination of methods above as the high-order one. The low-order solution and the antidiffusive fluxes (the difference between the high-order and low-order fluxes) are assembled in one pass (in a cycle over edges for the horizontal part and over vertices for the vertical part). We experimented with separate pre-limiting of horizontal and vertical antidiffusive fluxes and found that commonly this leads to an increased dissipation, for the horizontal admissible bounds are in many cases too tight. For this reason, the computation of admissible bounds and limiting is three-dimensional. As a result, it will not necessarily fully eliminate non-monotonic behavior in the horizontal direction. The basic difference from the FCT algorithm used in FESOM1.4 is the construction of low-order solution. In FESOM1.4 the low-order solution is obtained by adding an artificial diffusion to the high-order right hand side. Using the FCT roughly doubles the cost of transport algorithm, but makes the code more stabe in practice.

Vertical velocity splitting
---------------------------

As demonstrated in :cite:`Lemarie2015`, the strongest practical Courant number limitation is imposed by vertical advection in isolated patches adjacent to the coast. The code numerical efficiency can be improved if measures are taken to stabilize it with respect to sporadic events with large vertical velocities. Unstructured meshes may even be more vulnerable to such events because mesh irregularity can easily provoke a noisy pattern in :math:`w` just on its own. FESOM offers the approach proposed by :cite:`Shchepetkin2015` according to which the vertical transport velocity is split into two contributions :math:`w=w_{ex}+w_{im}` where the first one is determined by the maximum admissible Courant number, and the second one takes the rest. The advection with :math:`w_{ex}` is done explicitly using schemes mentioned above. The advection with :math:`w_{im}` is implicit. It uses the first-order upwind (backward Euler in time). This method leads to an operator that is  diagonally dominant. The implicit advective terms are added to the implicit vertical mixing terms and the resulting three-diagonal system of equations is solved with the standard sweep algorithm. Because of this, only very small additional costs incur if this algorithm is used. Although the first order upwind scheme is dissipative, it is applied only in critical cases to excessively large velocities.

Operator splitting
------------------

FESOM2 does not use operator splitting at present and takes the horizontal and vertical fluxes in a single step. However, from the viewpoint of increasing admissible time steps it is worthwhile to provide the implementation of advection in which tracers are updated separately for horizontal and vertical contributions. As is well known, the sequence horizontal-vertical should alternate with vertical-horizontal in this case. This work is planned, and this section will be updated in due course.

GM and isoneutral operators
===========================

The eddy-induced transport
--------------------------

FESOM2 follows  the algorithm proposed by :cite:`Ferrari2010` to implement the Gent-McWilliams (GM) parameterization :cite:`GentMcWilliams1990`,:cite:`Gent1995`. FESOM1.4 operates with skewsion (see :cite:`Griffiesbook` for mathematical detail). While working with skewsion is convenient in FESOM1.4 due to its variational formulation, it is less straightforward in FESOM2. Besides, the algorithm by :cite:`Ferrari2010` provides an explicit expression for the eddy bolus velocity streamfunction.

The bolus velocity :math:`{\bf v}^*=({\bf u}^*,w^*)` is expressed in terms of eddy-induced streamfunction :math:`\boldsymbol{\Psi}`,

.. math::
   {\bf v}^*=\nabla_3\times\boldsymbol{\Psi}, \quad \boldsymbol{\Psi}=\boldsymbol{\gamma}\times{\bf k},

where :math:`\boldsymbol{\gamma}` is a two-dimensional vector. In agreement with :cite:`Ferrari2010`, it is computed by solving

.. math::
   (c^2\partial_{zz}-N^2)\boldsymbol{\gamma}=(g/\rho_0)\kappa\nabla_z\sigma
   :label: eq_gm

with boundary conditions :math:`\boldsymbol{\gamma}=0` at the surface and ocean bottom. In this expression, :math:`c` is the speed of the first baroclinic mode, :math:`\sigma` the isoneutral density, :math:`\kappa` the thickness diffusivity, :math:`N` the Brunt–Väisälä frequency, and the index :math:`z` means that the gradient is computed for fixed :math:`z` (it differs from the gradient along layers, :math:`\nabla_z\sigma=\nabla\sigma-\partial_z\sigma\nabla Z`). In terms of the vector :math:`\boldsymbol{\gamma}` the components of eddy-induced velocity are computed as

.. math::
   {\bf u}^*=\partial_z\boldsymbol{\gamma}, \quad w^*=-\nabla\cdot\boldsymbol{\gamma}.

It is easy to see that solving :eq:`eq_gm` plays a role of tapering, for the solution is a smooth function satisfying boundary conditions.
The residual velocity :math:`{\bf u}_r={\bf u}+{\bf u}^*`, :math:`w_r=w+w^*` which is the sum of the eddy-induced velocity and the mean velocity :math:`({\bf u},w)` is consistent with :math:`\overline h` because the vertically integrated divergence of :math:`{\bf u}^*` is zero. The inclusion of eddy-induced velocity implies that the thickness and tracer equations are now written for the residual velocity :math:`{\bf u}_r`.

Although the natural placement for :math:`\boldsymbol{\gamma}` is at the cell centroids, it is moved to the mesh vertices in order to reduce the amount of computations. The vertical location is at full levels (layer interfaces). The horizontal bolus velocities are then computed at cell centroids as

.. math::
   {\bf u}^*_{c}=(1/3) \partial_z \sum_{v(c)}\boldsymbol{\gamma}_{v}.

The vertical bolus velocity :math:`w^*` is then found together with :math:`w` at the end of the ALE step and the full residual velocity is used to advect tracers.

We compute the speed :math:`c` in the WKB approximation as

.. math::
   c=\frac{1}{\pi}\int_{-H}^0Ndz.

Among other factors, the magnitude of the thickness diffusivity :math:`\kappa` depends on the resolution :math:`r` and the local Rossby radius :math:`L_R=c/f`:

.. math::
   \kappa=\kappa_0 f_{\kappa}(r/L_R),

where :math:`f_{\kappa}` is a cut-off function that tends to 0 if :math:`r/L_R<1` and to 1 otherwise. The resolution is defined as a square root of the area of the scalar control volume. On general meshes it may exhibit substantial local variations, so smoothing over the neighbor vertices is done.

Isoneutral diffusion
--------------------

Assuming that the slope of isopycnals is small, the diffusivity tensor can be written as

.. math::
   {\bf K}=
   \begin{pmatrix} K_i & 0 &s_xK_i \\
   0 & K_i & s_yK_i\\
   s_xK_i & s_yK_i & s^2K_i+K_d
   \end{pmatrix}
   :label: eq_kiso

Here :math:`K_i` and :math:`K_d` are the isoneutral and diapycnal diffusivities, and :math:`{\bf s}` is the isoneutral slope vector. Its derivatives are computed along layers,

.. math::
   {\bf s}=(s_x,s_y)=-\nabla\sigma/\partial_z\sigma.

If layer interfaces deviate substantially from geopotential surfaces, for example, if layers follow the bottom topography, the slope vector can be substantially larger than typically found on :math:`z`-coordinate meshes. Mixed derivatives in :math:`\nabla_3 h {\bf K}\nabla_3` operator in this case can limit  the time step :cite:`Lemarie2012a`. To maintain stability, the term :math:`h\partial_z(s^2K_i+K_d )\partial_z` is treated implicitly, as suggested by :cite:`Lemarie2012a`. Appendix :math:`app:isoneutral` shows the details of the numerical discretization of isoneutral diffusion.

Equation of state
-----------------

FESOM still works with potential temperature. The conservative temperature and TEOS10 will be made available soon. The estimates of density by the equation of state are made columnwise. To facilitate these estimates, for each column the arrays are computed of quantities appearing with different powers in :math:`z`. Then they are combined to estimate the in-situ density and pressure as well as to compute the Brunt–Väisälä frequency in the same routine.