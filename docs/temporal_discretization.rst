.. _temporal_discretization:

Temporal discretization: Asynchronous time stepping
***************************************************

FESOM2 uses asynchronous time stepping which means that velocity and scalar fields are shifted by a half time step. It is assumed that velocity is known at integer time steps and scalars at half-integer time steps. The time index is denoted by :math:`n`. The advantage  of this simple arrangement is that the discrete velocity is centered in time for the time step of scalar quantities, and pressure gradient is centered for the velocity time step. All other terms in equations require some interpolations to warrant second-order accuracy in time, and routinely the second-order Adams-Bashforth estimate is used in FESOM2. Other variants of time stepping are explored, and will be added in future.

Thickness and tracers
=====================

We advance thickness and scalar quantities as

.. math::
   h^{n+1/2}-h^{n-1/2}=-\tau[\nabla\cdot({\bf u}^nh^{*})+w^t-w^b+ W^{n-1/2}\delta_{k1}]
   :label: eq_ht

and

.. math::
   h^{n+1/2}T^{n+1/2}-h^{n-1/2}T^{n-1/2}=-\tau[\nabla\cdot({\bf u}^nh^{*}T^n)+w^tT^t-w^bT^b+W^{n-1/2}T_W\delta_{k1}]+ D_T.
   :label: eq_tracert

Here :math:`\tau` is the time step and :math:`D_T` stands for the terms related to diffusion. The time index on :math:`w` is not written, for :math:`w` is diagnosed by :math:`{\bf u}` and :math:`h`.

.. note::
   Note that :math:`\mathbf{u}` and :math:`h` enter similarly in the equations for thickness and tracer. This warrants consistency: if :math:`T=\rm{const}`, the tracer equation reduces to the thickness equation.

Since the horizontal velocity is centered in time, these equations will be of the second order for the advective terms if :math:`h^*` is centered (:math:`h^*=h^n`) or if any other high-order estimate is used. The treatment of :math:`h^*` depends on options. At present FESOM allows only the simplest choices when :math:`h` vary only slightly with respect to unperturbed thickness and :math:`\partial_th` is prescribed by the evolution of :math:`\eta`. We take :math:`h^*=h^{n-1/2}` in this case because this choice can easily be made consistent with elevation, see further.

Although this formally reduces the time order to the first, the elevation is usually computed with the accuracy shifted to the first-order in large-scale ocean models, including this one.

.. note::
   Other options, including those allowing :math:`h` to follow isopycnal dynamic will be gradually added to FESOM. They will differ in the way how :math:`h^*` is computed.


Elevation
=========

We introduce

.. math::
   \overline{h}=\sum_k h_k-H,

where :math:`H` is the unperturbed ocean thickness. :math:`\overline h` would be identical to the elevation :math:`\eta` in the continuous world, but not in the discrete formulation here.
The difference between these two quantities is that the elevation is defined at integer time levels. More importantly, it has to be computed so that the fast external mode is filtered. FESOM2 uses implicit method for that. The point is to make :math:`\overline h` and :math:`\eta` fully consistent.

In order to filter the external mode :math:`\eta` is advanced implicitly in time. For :math:`h^*=h^{n-1/2}` we write for the elevation

.. math::
   \eta^{n+1}-\eta^n=-\tau(\alpha(\nabla\cdot\sum_k{h}_k^{n+1/2}{\bf u}_k^{n+1}+W^{n+1/2})+(1-\alpha)(\nabla\cdot\sum_k{h}_k^{n-1/2}{\bf u}_k^{n}+W^{n-1/2})).
   :label: eq_etat

Here :math:`\alpha` is the implicitness parameter (:math:`0.5\le\alpha\le1`) in the continuity equation. Note that the velocities at different time steps are taken with their respective thicknesses in the same way as they appear in the thickness equation :eq:`eq_ht`. The approach below is inspired by :cite:`Campin2004`. The equation for thicknesses can be vertically integrated giving, under the condition that the surface value of :math:`w^t` vanishes,

.. math::
   \overline{h}^{n+1/2}-\overline{h}^{n-1/2}=-\tau\nabla\cdot\sum_k{h}_k^{n-1/2}{\bf u}_k^n-\tau W^{n-1/2}.
   :label: eq_hbar

Expressing the rhs in the formula for :math:`\eta` :eq:`eq_etat` through the difference in surface displacements :math:`\overline{h}` from the last formula we see that :math:`\eta` and :math:`\overline{h}` can be made consistent if we require

.. math::
   \eta^n=\alpha \overline{h}^{n+1/2}+(1-\alpha)\overline{h}^{n-1/2}.
   :label: eq_etan


To eliminate the possibility for :math:`\eta` and :math:`\overline{h}` to diverge, we compute :math:`\eta^n` from the last formula, then estimate :math:`\eta^{n+1}` by solving dynamical equations (equation :eq:`eq_etat` requires :math:`\mathbf{u}^{n+1}`, so it is solved simultaneously with the momentum equation), and use it only to compute :math:`{\bf u}^{n+1}`. On the new time step a 'copy' of :math:`\eta^{n+1}` will be created from the respective fields :math:`\overline{h}`.
We commonly select :math:`\alpha=1/2`, in this case :math:`\eta^n` is just the interpolation between the two adjacent values of :math:`\overline{h}`. Note that :eq:`eq_etan` will be valid for any estimate :math:`h^*` as far as it is used consistently in the product with the horizontal velocity.

.. note:;
   The implicit way of solving for :math:`\eta` means that FESOM uses iterative solver at this step. Such solvers are thought to be a major bottleneck in massively parallel applications. This is the reason why many groups abandon solvers and go for subcycling+filtering algorithms for the treatment of external (approximately barotropic) dynamics. Such an option will be added to FESOM. Its elementary variant is described in Appendix. We work on more advanced variants.

Momentum equation
=================

The momentum equation has to be solved together with the elevation equation :eq:`eq_etat`, which is done with a predictor-corrector method. The method is largely an adaptation of pressure correction algorithm from computational fluid dynamics.

Assuming the forms :eq:`eq_mom_vei` or :eq:`eq_mom_f2` we write (using :math:`\partial_z` for brevity instead of the top-bottom differences)

.. math::
   {\bf u}^{n+1}-{\bf u}^{n}=\tau({\bf R}^{n+1/2}_u+\partial_z\nu_v\partial_z{\bf u}^{n+1}-g\nabla(\theta\eta^{n+1}+(1-\theta)\eta^n)).

Here :math:`\theta` is the implicitness parameter for the elevation, :math:`{\bf R}^{n+1/2}_u` includes all the terms except for vertical viscosity and the contribution from the elevation which are treated implicitly. To compute :math:`{\bf R}^{n+1/2}_u`, we use the second-order Adams-Bashforth (AB) method for the terms related to the momentum advection and Coriolis acceleration. The AB estimate of quantity :math:`q` is

.. math::
   q^{AB}=(3/2+\epsilon)q^n-(1/2+\epsilon)q^{n-1}.

Here :math:`\epsilon` is a small parameter (:math:`\le0.1`) needed to ensure stability in the case of advection operators. The contribution of pressure :math:`P` does not need the AB interpolation (because it is centered). The horizontal viscosity is estimated on the level :math:`n` because this term is commonly selected from numerical, not physical reasons.

- We write the predictor equation

   .. math::
      {\bf u}^{*}-{\bf u}^{n}-\tau\partial_z\nu_v\partial_z({\bf u}^{*}-{\bf u}^n)=\tau({\bf R}^{n+1/2}_u+\partial_z\nu_v\partial_z{\bf u}^{n}-g\nabla\eta^n).
      :label: eq_predict

  The operator on the lhs connects three vertical levels, leading to three-diagonal linear problem for :math:`\Delta {\bf u}=\mathbf{u}_k^*-\mathbf{u}_k^n` for each vertical column. Solving it we find the predicted velocity update :math:`\Delta {\bf u}`. (The vertical viscosity contribution on the rhs is added during the assembly of the operator on the lhs.)


- The corrector step is written as

  .. math::
     {\bf u}^{n+1}-{\bf u}^{*}=-g\tau\theta\nabla(\eta^{n+1}-\eta^n).
     :label: eq_cor

- Expressing the new velocity from the last equation and substituting it into the equation for the elevation :eq:`eq_etat`, we find

  .. math::
     \frac{1}{\tau}(\eta^{n+1}-\eta^n)-\alpha\theta g\tau\nabla\cdot(\overline{h}^{n+1/2}+H)\nabla(\eta^{n+1}-\eta^n)dz= \nonumber \\
     -\alpha(\nabla\cdot\sum_k{h}_k^{n+1/2}({\bf u}^n+\Delta{\bf u})_k+W^{n+1/2})-(1-\alpha)(\nabla\cdot\sum{h}_k^{n-1/2}{\bf u}_k^n+W^{n-1/2}).
     :label: eq_etaa

  Here, the operator part depends on :math:`h^{n+1/2}`, which is the current value of thickness. The last term on the rhs is taken from the thickness computations on the previous time step.

The overall solution strategy is as follows.

- Compute :math:`\eta^n` from :eq:`eq_etan`. Once it is known, compute :math:`\Delta {\bf u}` from :eq:`eq_predict`.

- Update the matrix of the operator on the lhs of :eq:`eq_etaa`. Solve :eq:`eq_etaa` for :math:`\eta^{n+1}-\eta^n` using an iterative solver and estimate the new horizontal velocity from :eq:`eq_cor`.

- Compute :math:`\overline{h}^{n+3/2}` from :eq:`eq_hbar`.

- Determine layer thicknesses and :math:`w` according to the options described below.

- Advance the tracers. The implementation of implicit vertical diffusion will be detailed below.

Momentum equation in form :eq:`eq_mom_fl`
=========================================

Here an additional difficulty is the presence of :math:`h` in the time derivative and on the rhs. The rule is that :math:`\mathbf{u}` should appear with the same :math:`h^*` as in the thickness or tracer equation. We used thus far the choice :math:`h^*\mathbf{u}^n=h^{n-1/2}\mathbf{u}^n` in these equations, which implies that the time derivative will be

.. math::
   \partial_t(\mathbf{u}h)=( h^{n+1/2}\mathbf{u}^{n+1}-h^{n-1/2}\mathbf{u}^n)/\tau,

:math:`h^{n-1/2}` will be used on the rhs with pressure gradient term, and the predictor equation will be written for :math:`h^{n-1/2}\Delta\mathbf{u}`. In this case :math:`h^{n-1/2}` can be factored out of the lhs, which will make predictor solution similar. The corrector step will be modified to

.. math::
   h^{n+1/2}{\bf u}^{n+1}-h^{n-1/2}{\bf u}^{*}=-gh^{n-1/2}\tau\theta\nabla(\eta^{n+1}-\eta^n).
   :label: eq_corf1


It will lead to the replacement of :math:`h^{n+1/2}` in the lhs of :eq:`eq_etaa` by :math:`h^{n-1/2}`. We stress once again that the expressions in this section are for the particular choice of :math:`h^*`.

Current options for the vertical coordinate
===========================================

The following options for the vertical coordinate are available at present:

- Linear free surface: If we keep the layer thicknesses fixed, the time derivative in :eq:`eq_ht` drops out, and it becomes the standard equation to compute :math:`w`, starting from the bottom and continuing to the top,


  .. math::
     w^t-w^b+\nabla\cdot({h\bf u})=0.


 If this option is applied also to the first layer, the freshwater flux cannot be taken into account in the thickness equation. Its contribution to the salinity equation is added through a virtual salinity flux. In this option, :math:`w` at the (fixed) ocean surface differs from zero, and so do the tracer fluxes. They do not necessarily integrate to zero over the ocean surface which is why tracer conservation is violated.

- Full (nonlinear) free surface: We adjust the thickness of the upper layer, while the thicknesses of all other layers are kept fixed, :math:`\partial_th_k=0` for :math:`k>1`. The thickness equations are used to compute :math:`w` on levels :math:`k=2:K_v` starting from the bottom. The change in the thickness of the first layer :math:`h^{n+3/2}_1-h^{n+1/2}_1` is given by :eq:`eq_hbar` written for the respective time interval. In this case there is no transport through the upper moving surface (the transport velocity :math:`w_1` is identically zero). This option requires minimum adjustment with respect to the linear free surface. However, the matrix of the operator in :eq:`eq_etaa` needs to be re-assembled on each time step.

- We can distribute the total change in height :math:`\partial_t\overline h` between several or all eligible layers. Due to our implementation, at  *each* scalar horizontal location they can only be the layers that do not touch the bottom topography. If all eligible layers are involved we estimate


  .. math::
     \partial_t h_k=(h_k^0/\tilde H)\partial_t\overline h,


  where :math:`h_k^0` are the unperturbed layer thicknesses and :math:`\tilde H` is their sum for all eligible layers. The thickness of the layers adjacent to the topography is kept fixed. The equation on thickness, written for each layer, is used to compute transport velocities :math:`w` starting from zero bottom value. This variant gives the so-called :math:`z^*`-coordinate.

- Additional options will be gradually added. Layer thicknesses can vary in many ways provided that their tendencies sum to :math:`\partial_t\overline h` over the layers. In particular, requiring that transport velocities :math:`w` are zero, isopycnal layers can be introduced. The levels can move with high-pass vertical velocities, leading to the so called :math:`\tilde z` coordinate, see :cite:`Leclair2011`, :cite:`Petersen2015` or follow density gradients as in :cite:`Hofmeister2010`. The unperturbed layer thicknesses need not follow the geopotential surfaces and can be terrain following for example.

- The ALE vertical coordinate is only a framework where many options are in principle possible. Additional measures may be required in each particular case, such as computations of pressure gradients with reduced errors. Updated transport algorithms may be needed (in the spirit of :cite:`Lemarie2012b` to minimize spurious numerical mixing in terrain-following layers. These generalizations are among the topics of ongoing work.


Implicit vertical diffusion
===========================

We return to the tracer equation :eq:`eq_tracert`. The vertical diffusion in this equation may present a CFL limitation and is treated implicitly.

Because of varying layer thicknesses, the implementation of implicit vertical diffusion needs slight adjustment compared to the common case of fixed layers. We write, considering time levels :math:`n-1/2` and :math:`n+1/2`,

.. math::
   h^{n+1/2}T^{n+1/2}-h^{n-1/2}T^{n-1/2}=\tau(R_T^{n}+(K_{33}\partial_zT^{n+1/2})^t-(K_{33}\partial_zT^{n+1/2})^b)


and split it into

.. math::
   h^{n+1/2}T^{*}-h^{n-1/2}T^{n-1/2}=\tau R_T^{n}

and

.. math::
   h^{n+1/2}(T^{n+1/2}-T^{*})=\tau(K_{33}\partial_z(T^{n+1/2}-T^*)+K_{33}\partial_zT^*)|^t_b.

Here :math:`R_T` contains all advection terms and the terms due to the diffusion tensor except for the diagonal term with :math:`K_{33}`. The preliminary computation of :math:`T^*` is necessary to guarantee that a uniform tracer distribution stays uniform (some significant digits will be lost otherwise).
