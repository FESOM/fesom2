.. _time_stepping_transport:

Time stepping for the transport :math:`\mathbf{U}=h\mathbf{u}` instead of velocity
**********************************************************************************

For the momentum equation in the form :eq:`eq_mom_fl` an alternative variant of time stepping is possible when the quantity :math:`\mathbf{U}=h\mathbf{u}` is advanced instead of  velocity :math:`\mathbf{u}`. This will simultaneously imply that :math:`h^*=h^{n}`, making the thickness and transport equations centered with respect to time step. The thickness appearing with the pressure gradient should be them :math:`h^{n+1/2}`, which provides a centered estimate. The advection and Coriolis terms are computed through AB2 (or AB3) time stepping, and if needed, the Coriolis term can be time stepped semiimplicitly.

The time stepping algorithm can be formulated as follows

.. math::
   {\bf U}^{n+1}-{\bf U}^{n}=\tau({\bf R}_{U}^{n+1/2}-gh^{n+1/2}\nabla(\theta\eta^{n+1}+(1-\theta)\eta^n)+(\nu_v\partial_z{\bf u}^{n+1})^t-(\nu_v\partial_z{\bf u}^{n+1})^b)

with

.. math::
   {\bf R}_{U}^{n+1/2}=({\bf R}_{U}^*)^{AB}-h^{n+1/2}(\nabla p_h+g\rho\nabla Z)/\rho_0,

and

.. math::
   {\bf R}_{U}^*=-\nabla\cdot({\bf U}^n{\bf u}^n)-(w^t{\bf u}^t-w^b{\bf u}^b)^n-f{\bf k}\times{\bf U}^n.

The last expression combines the terms that need the AB method for stability and the second order. We use :math:`h^{n+1/2}` to compute :math:`Z` and follow the same rule as :eq:`eq_etan` to compute :math:`\eta^n`. The steps are:

- Do the predictor step and compute :math:`\Delta \tilde{\bf U}=\tau{\bf R}_U^{n+1/2}-\tau gh^{n+1/2}\nabla\eta^n`.

- Update for implicit viscosity.

  .. math::
     \partial_t\Delta{\bf U}-(\nu_v\partial_z(\Delta{\bf U}/h^{n+1/2}))|^t_b=\Delta\tilde{\bf U}+(\nu_v\partial_z({\bf U}^n/h^{n+1/2}))|^t_b.

- Solve for new elevation. We write first

  .. math::
     \overline{\bf U}=\sum_k{\bf U},

  and similarly for other quantities, getting

  .. math::
     \overline{\bf U}^{n+1}-\overline{\bf U}^n=\overline{\Delta{\bf U}}-g\tau(H+\overline h^{n+1/2})\theta\nabla(\eta^{n+1}-\eta^n)
     :label: eq_baru

and

  .. math::
     \eta^{n+1}-\eta^n=-\tau\nabla\cdot(\alpha\overline{\bf U}^{n+1}+(1-\alpha)\overline{\bf U}^{n})-\tau(\alpha W^{n+1/2}+(1-\alpha)W^{n-1/2}).
     :label: eq_etaU

Eliminating :math:`\overline{\bf U}^{n+1}` between these two equations, one gets the equation on elevation increment :math:`\Delta\eta=\eta^{n+1}-\eta^n`

  .. math::
     \Delta\eta-g\tau^2\theta\alpha\nabla\cdot((H+\overline h^{n+1/2})\nabla\Delta\eta)=-\tau\nabla\cdot(\alpha\overline{\Delta{\bf U}}+\overline{\bf U}^n)-\tau(\alpha W^{n+1/2}+(1-\alpha)W^{n-1/2})

  In reality, everything remains similar to the vector-invariant case, and the matrix to be inverted is the same.

- Correct the transport velocities as

  .. math::
     {\bf U}^{n+1}-{\bf U}^n={\Delta{\bf U}}-g\tau h^{n+1/2}\theta\nabla\Delta\eta.
     :label: eq_corrU

- Proceed with ALE and determine :math:`w^{n+1}`, :math:`h^{n+3/2}`, :math:`T^{n+3/2}`.

- The new velocities are estimated as

  .. math::
     {\bf u}^{n+1}={\bf U}^{n+1}/h^{n+1}.

  Here :math:`h^{n+1}` can be computed either in the agreement with the ALE procedure (:math:`\eta^{n+1}` is already known) or interpolating between :math:`n+1/2` and :math:`n+3/2` time levels.


This alternative form of time stepping is more elegant. The horizontal velocity appears in most places in the layer equations as the product with respective thickness, and the alternative form takes this naturally into account. It will be added in due time together with the development of ALE options.