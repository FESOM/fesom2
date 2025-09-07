.. _subcycling_instead_solver:

Subcycling instead of solver
****************************

Semi-implicit treatment of external mode has a drawback of suboptimal parallel scalability in the limit of very small partitions. An alternative approach is to use split explicit time stepping when the external mode (elevation and vertically integrated velocity) are time stepped with a small step (subcycling), and then filtered to remove fast contributions. This option will be added in the future, and at present optimal algorithms are explored. The description of this section gives one possibility. Flux form of momentum advection is used. We take

.. math::
   \eta^n=(\overline h^{n-1/2}+\overline h^{n+1/2})/2,

since it provides the second-order accurate estimate.

An easiest approach is to run subcycles between time levels :math:`n` and :math:`n+2`, with subsequent averaging to level :math:`n+1`.

The contribution from the elevation :math:`\eta^n` is kept in the predicting :math:`\Delta \tilde{\bf U}` because it also incorporates the implicit solve for vertical viscosity. Then the compensation term with :math:`\eta^n` appears in :eq:`eq_barus` below. This can be avoided if implicit vertical viscosity substep is moved to the end of velocity step.

Instead of :eq:`eq_baru` and :eq:`eq_etaU` we introduce subcycles indexed with :math:`j`, :math:`j=0:2J`, with :math:`\eta^{n+j/J}` shortcut to :math:`\eta^j` and same for :math:`\overline{\bf U}` in several formulas below. The simplest form of subcycling looks like

.. math::
   \eta^{j+1}-\eta^j=-(\nabla\cdot\overline{\bf U}^{j}+W^j)\tau/J.
   :label: eq_etas

.. math::
   \overline{\bf U}^{j+1}-\overline{\bf U}^j=\overline{\Delta{\bf U}}/J-g(\tau/J)(H+\overline h^{n+1/2})\nabla(\eta^{j+1}-\eta^n).
   :label: eq_barus

This is a forward--backward scheme.

Other forms of subcycling can be used to increase stability and reduce the number of subcycles :math:`2J+1`. Many of them are discussed by :cite:`Shchepetkin2005`. In particular, an AB3-AM4 scheme (see also :cite:`Lemarie2015` is demonstrated to provide good accuracy and stability.

On completing sybcycles one is at time level :math:`n+2`. In order to eliminate possible high frequencies, averaging is done to time level :math:`n+1`:

.. math::
   \overline{\bf U}^{n+1}=(2J+1)^{-1}\sum_j\overline{\bf U}^j,\quad \eta^{n+1}=(2J+1)^{-1}\sum_j\eta^j.

The common further action is to use :math:`\overline{\bf U}^{n+1}` for the barotropic transport combined with the baroclinic transport diagnosed from :math:`{\bf U}^{n+1}`. We introduce first the new baroclinic transport by writing

.. math::
   {\bf U}^*_k={\bf U}^n_k+\Delta{\bf U}_k,

.. math::
   \tilde{\bf U}^{n+1}_k={\bf U}^*_k
   -\overline{{\bf U}}^*\frac{h^{n+1}_k}{H+\eta^{n+1}}.

It is then updated to the full transport velocity by

.. math::
   {\bf U}^{n+1}_k=\tilde{\bf U}^{n+1}_k+\overline{{\bf U}}^{n+1}\frac{h^{n+1}_k}{H+\eta^{n+1}}.

Here :math:`h_k^{n+1}` is an estimate of layer thickness at time step :math:`n+1`.

A recent suggestion is to replace the time stepping in :eq:`eq_etas`-:eq:`eq_barus` by a dissipative one modifying :eq:`eq_barus` as

.. math::
   \overline{\bf U}^{j+1}-\overline{\bf U}^j=\overline{\Delta{\bf U}}/J-g(\tau/J)(H+\overline h^{n+1/2})\nabla((1+\lambda)\eta^{j+1}-\lambda \eta^{j}-\eta^n).
   :label: eq_barusm

The parameter :math:`0\le \lambda<1` controls the dissipation which alone can be sufficient to remove the high-frequency component in :math:`\overline{\bf U}` and :math:`\eta`. It remains to be seen whether this is sufficient to fully eliminate averaging and shorten integration just to :math:`n+1` instead of :math:`n+2`.