.. _vertical_discretization

Vertical discretization: Layer thicknesses and layer equations
**************************************************************

FESOM2 uses Arbitrary Lagrangian Eulerian (ALE) vertical coordinate. This implies that level surfaces are allowed to move. ALE vertical coordinate on its own is only the framework enabling moving level surfaces. The way how they are moving depends on one's particular goal and may require additional algorithmic steps. Two limiting cases obviously include the case of fixed :math:`z`-levels and the case when levels are isopycnal surfaces. At present only vertical coordinates that slightly deviate from :math:`z`-surfaces are supported in FESOM, but many other options will follow.

The implementation of ALE vertical coordinate in FESOM2 basically follows  :cite:`Ringler2013`. An alternative approach, used by MOM6, see, e.g., :cite:`Adcroft_Hallberg_2006` :cite:`Adcroft2019`, is in the exploratory phase.

The essential step toward the ALE vertical coordinate lies in confining equations of section :ref:`sec_cequations` to model layers.

- Introduce layer thicknesses :math:`h_k=h_k(x,y,t)`, where :math:`k=1:K` is the layer index and :math:`K` the total number of layers. They are functions of the horizontal coordinates and time in a general case. Each layer consists of prisms defined by the surface mesh but partly masked by bottom topography.

- Layers communicate via the transport velocities :math:`w_{kv}` through the top and bottom boundaries of the prisms. The transport velocities are the differences between the physical velocities in the direction normal to the layer interfaces and the velocities due to the motion of the interfaces. These velocities are defined at the interfaces (the yellow points in :numref:`vertical`). For layer :math:`k` the top interface has index :math:`k` and the bottom one is :math:`k+1`. Note that :math:`w_{kv}` coincides with the vertical velocity only if the level surfaces are flat.

- All other quantities - horizontal velocities :math:`{\bf u}`, temperature :math:`T`, salinity :math:`S` and pressure :math:`p` are defined at mid-layers. Their depths will be denoted as :math:`Z_k`, and the notation :math:`z_k` is kept for the depths of mesh levels (the layer interfaces). They are functions of horizontal coordinates and time in a general case.

The equations of motion, continuity and tracer balance are integrated vertically over the layers. We will use :math:`T` as a representative of an arbitrary tracer.


- The continuity equation becomes the equation on layer thicknesses

.. math::
   \partial_t h_k+\nabla\cdot({\bf u}h)_k+(w^{t}-w^b)_k+W\delta_{k1}=0,
   :label: eq_thickness


- and the tracer equation becomes

.. math::
   \partial_t(hT)_k+\nabla\cdot({\bf u}hT)_k+(w^{t}T^t-w^bT^b)_k+WT_W\delta_{k1}=\nabla\cdot h_k{\bf K}\nabla T_k.
   :label: eq_tracer


Here, :math:`W` is the water flux leaving the ocean at the surface, it contributes to the first layer only (hence the delta-function); :math:`T_W` is the property transported with the surface water flux and the indices :math:`t` and :math:`b` imply the top and the bottom of the layer.

The right hand side of :eq:`eq_tracer` contains the 3 by 3 diffusivity tensor :math:`{\bf K}`. We still use :math:`\nabla` in :eq:`eq_tracer` for the 3D divergence (the outer :math:`\nabla`) for brevity, but assume the discrete form :math:`\nabla_h(...)+((...)^t-(...)^b)/h_k`, where :math:`(...)` are the placeholders for the horizontal and vertical components of 3D vector it acts on. A correct discretization of the diffusivity term is cumbersome and will be explained below.

- Vertical sum of :eq:`eq_thickness` over layers with account that :math:`w^t=0` at the free surface and :math:`w_b=0` at the bottom gives the 'layered' version of the elevation equation

   .. math::
      \partial_t\eta+\nabla_h\cdot\sum_kh_k{\bf u}_k+W=0.
      :label: eq_eta

- The layer-integrated momentum equation in the flux form is

      .. math::
         \partial_t(h{\bf u})+\nabla_h\cdot(h{\bf u u})+w^t{\bf u}^t-w^b{\bf u}^b+
         f{\bf k}\times{\bf u}h +h(\nabla_h p+g\rho\nabla Z)/\rho_0=  \nonumber \\ D_{uh}{\bf u}+(\nu_v\partial_z{\bf u})^t-(\nu_v\partial_z{\bf u})^b,
         :label: eq_mom_fl

  with :math:`D_{uh}{\bf u}` the horizontal viscosity operator for the flux form (to be specified later), :math:`\nu_v` the vertical viscosity coefficient, :math:`f` the Coriolis parameter and :math:`{\bf k}` a unit vertical vector. We ignore the momentum source due to the added water :math:`W` at the surface. Note that it could be more natural to formulate the solution procedure in terms of the horizontal layer transport velocities :math:`{\bf U}=h{\bf u}` in this case, but the present implementation in FESOM deals with separate :math:`h` and :math:`\mathbf{u}`.

- The pressure field is expressed as

      .. math::
         p=g\rho_0\eta+P, \quad P_{1}=p_a+g\rho_1h_1/2, \quad P_k=P_{k-1}+g(\rho_{k-1}h_{k-1}+ \rho_kh_k)/2.
         :label: eq_pressure

  with :math:`p_a` the atmospheric pressure, :math:`\rho` the deviation of density from its reference value :math:`\rho_0`, and :math:`P` is the variable hydrostatic pressure due to :math:`\rho`. The pressure gradient in continuous equations :eq:`eq_cmom` has to be computed at constant :math:`z`. The model levels deviate from surfaces :math:`z=\rm{const}`. The term :math:`g\rho\nabla Z`, appearing together with the horizontal pressure gradient in :eq:`eq_mom_fl` compensates for the deviation. The quantity :math:`Z` appearing in this term is the :math:`z`-coordinate of the midplane of the layer with the thickness :math:`h`.

.. note::
   Although :math:`\nabla p+g\rho\nabla Z` gives a formally correct estimate of pressure gradient at constant :math:`z`, the errors of discretization of the two terms in this expression become an issue if level surfaces deviate from :math:`z`-surfaces. They are known as pressure gradient errors and require special care.  FESOM will propose a selection of known algorithms, including the finite-volume algorithms of pressure gradient force that follows :cite:`Engwirda2017` but is adapted to the triangular prisms of FESOM mesh.

- Instead of using the flux form of momentum equation :eq:`eq_mom_fl` representing momentum balance in the layer one can work without layer integration. Of particular interest is the vector-invariant form written as

  .. math::
     \partial_t{\bf u}+\frac{\omega+f}{h}{\bf k}\times{\bf u}h+((w\partial_z{\bf u})^t+(w\partial_z{\bf u})^b)/2 +\nabla (p/\rho_0+{\bf u}^2/2)+g\rho\nabla Z/\rho_0= \nonumber \\ D_u{\bf u}+((\nu_v\partial_z{\bf u})^t-(\nu_v\partial_z{\bf u})^b)/h.
     :label: eq_mom_vei

  Here, the identity

  .. math::
     {\bf u}\cdot\nabla{\bf u}=\omega{\bf k}\times{\bf u}+\nabla({\bf u}^2/2),\quad \omega={\bf k}\cdot(\nabla\times{\bf u})

  was used.

- The second term on the lhs of :eq:`eq_mom_vei` includes division and multiplication with the layer thickness, and in doing so, it introduces the layer potential vorticity (PV), :math:`q=(\omega+f)/h` and its transport :math:`{\bf u}h`. The layer thickness formally drops out from the equation :eq:`eq_mom_vei` which is still continuous in the horizontal direction. However, in the discrete case, the location of vorticity points (vertices) and velocity points is different. By keeping separate :math:`h` the equation will then operate on the same horizontal transports as the thickness equations. This is the prerequisite for developing discretizations that conserve potential vorticity.

- One more form is possible where the vector-invariant representation is not used

      .. math::
         \partial_t({\bf u})+\nabla\cdot({\bf u u})+(w^t{\bf u}^t-w^b{\bf u}^b)/h+
         f{\bf k}\times{\bf u} +(\nabla p+g\rho\nabla Z)/\rho_0=  \nonumber \\ D_{u}{\bf u}+(A_v\partial_z{\bf u})^t-(A_v\partial_z{\bf u})^b/h.
         :label: eq_mom_f2

The default version in FESOM2 is :eq:`eq_mom_fl`. Although the versions are derived from the same continuous equations, they are not equivalent in the discrete case.

