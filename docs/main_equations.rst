.. _main_equations:

Main equations
**************

.. _sec_cequations:

Main equations
==============

The code solves the standard set of equations derived under the standard set of approximations (Boussinesq, hydrostatic, and traditional approximations).
These equations include the momentum equation for horizontal velocity

.. math::
   \partial_t\mathbf{u}+f\mathbf{e}_z\times\mathbf{u}+(\mathbf{u}\cdot\nabla_h+w\partial_z)\mathbf{u}+\nabla_h p/\rho_0=D_h\mathbf{u}+\partial_z\nu_v\partial_z\mathbf{u},
   :label: eq_cmom

the hydrostatic equation

.. math::
   \partial_zp=-g\rho,
   :label: eq_chydrost

the Boussinesq form of the continuity equation

.. math::
   \partial_zw=-\nabla\cdot\mathbf{u},
   :label: eq_ccont

and the equations for potential temperature (FESOM is still using potential temperature, which will be replaced by conservative temperature in the nearest future) and salinity

.. math::
   \partial_t T+\nabla\cdot (\mathbf{u}T)+\partial_z(wT)=\nabla\cdot\mathbf{K}\nabla T,
   :label: eq_cT}


.. math::
   \partial_t S+\nabla\cdot (\mathbf{u}S)+\partial_z(wS)=\nabla\cdot\mathbf{K}\nabla S.
   :label: eq_cS}

In these equations :math:`\mathbf{u}=(u,v)` is the horizontal velocity, :math:`f` the Coriolis parameter, :math:`\mathbf{e}_z` the unit vertical vector, :math:`\rho_0` the reference density, :math:`p` the pressure, :math:`D_h` the horizontal viscosity operators to be specified further, :math:`\nu_v` the vertical viscosity coefficient, :math:`g` the gravitational acceleration, :math:`T, S`, the potential temperature and salinity and :math:`\mathbf{K}` the diffusivity tensor directing mixing in deep ocean to be isoneutral. The operator :math:`\nabla` is two-dimensional, :math:`\nabla_h=(\partial_x,\partial_y)`, and :math:`\nabla=(\nabla, \partial_z)`. The equations above have to be complemented by the equation of state connecting density with the temperature, salinity and pressure. In the Boussinesq approximation, the pressure featuring in the equation of state is the fluid depth up to a factor :math:`g\rho_0`, so we formally write

.. math::
   \rho=\rho(T,S,z).

They also need appropriate initial and boundary conditions. The walls and bottom of the ocean basin are traditionally considered as isolated and impermeable, implying no flux boundary conditions. Flux conditions are imposed on the surface. Bottom acts as a momentum sink through the drag force applied there.

Ocean free surface denoted :math:`\eta` varies in space and time. An equation governing it is obtained by integrating :eq:`eq_ccont` vertically from the bottom at :math:`z=-H(x,y)` to the top at :math:`z=\eta(x,y,t)`

.. math::
   \partial_t\eta+\nabla_h\int^{\eta}_{-H}\mathbf{u}dz=-W,
   :label: eq_ceta

where :math:`W` is water flux leaving the ocean through the surface (specified as a part of boundary conditions). We stress that the last equation is not an independent one, but the consequence of the equations written before.

Transition from continuous to discrete equation includes the steps of spatial and temporal discretization. The spatial discretization is very different for the vertical and horizontal direction and is treated separately. We begin with vertical discretization, followed by temporal discretization and then by the horizontal discretization.

