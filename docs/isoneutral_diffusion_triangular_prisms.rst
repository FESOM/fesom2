.. _isoneutral_diffusion_triangular_prisms:

Isoneutral diffusion on triangular prisms
*****************************************

A rigorous implementation of isoneutral diffusion operator is a difficult task because one has to ensure that it provides decay of discrete variance and that it is exactly zero when applied to the isoneutral density.

The horizontal and vertical components of fluxes related to (\ref{eq:kiso}) are

.. math::
   {\bf F}_h(T)=-K_i(\nabla_h T+{\bf s}\partial_z T),


.. math::
   F_z(T)=-K_i({\bf s}\nabla_h T+s^2\partial_zT)-K_d\partial_zT.

The terms including :math:`K_i` are referred to as the isoneutral flux, the
remaining term with :math:`K_d` is the dianeutral flux.
To complete the description, the slope :math:`{\bf s}` has to be expressed in terms of thermal expansion and saline contraction coefficients :math:`\alpha`
and :math:`\beta`,

.. math::
   {\bf s}=-\frac{-\alpha\nabla_h T+\beta \nabla_h S}{-\alpha\partial_z
   T+\beta\partial_z S}.

(Do not mix this :math:`\alpha` with the other one used in the time stepping.)

The implementation difficulty stems from the fact that the tracers together with :math:`\alpha` and :math:`\beta` are located at mid-layers, the vertical derivatives are located at the level surfaces, and the horizontal derivatives are at mid-layers, but at cells instead of vertices. The estimate of slope at a single point is impossible without extra interpolation, which will break full consistency. The solution involves 'triads' (see, e.g., :cite:`Griffiesbook` and :cite:`Lemarie2012a` and variational formulation. Note, however, that the implicit time stepping of the contribution with :math:`s^2K_i` in the vertical flux, needed for stability reasons :cite:`Lemarie2012a`, will introduce some errors even in this case.

First, we split each triangular prism of our mesh into subvolumes characterized by unique values of the expansion/contraction coefficients, vertical gradients and horizontal gradients, to form the triads. We obtain 6 subprisms per prism, formed by sections along midplane and by vertical planes passing through centroids and mid-edges.

Next, the dissipation functional is written. We will use different, but equivalent formulation, which would follow if tracer equations were written in a weak form. Consider the bilinear form

.. math::
   6\mathcal{F}(\tilde T,T)=-\sum_{k,c}\sum_{p=1}^{p=6}A_ch_{kc}(\nabla \tilde{T}
   {\bf K}\nabla T)_{kcp}.

Here the first summation is over mesh prisms (cells and layers), and the second one, over the subprisms :math:`p`. The volume of each subprism is 1/6 of the volume of the full prism (hence the factor 6 on the lhs). Clearly, :math:`2\mathcal{F}(T,T)` corresponds to total variance dissipation. If :math:`T` is the isoneutral density and its gradients are expressed in terms of :math:`\alpha` and :math:`\beta` as for the slope above, :math:`\mathcal{F}` vanishes.

The last step is to compute the contribution to the rhs of scalar equation from the diffusion term

.. math::
   (R_T)_{kv}=(1/A_{kv})\partial\mathcal{F}/\partial \tilde T_{kv}.

Since we deal with layer-integrated equations, the division is over the area of scalar cell :math:`v` instead of division by volume. Writing down the expression for :math:`R_T` is a rather tedious task. The result can be reformulated in terms of the discrete divergence of discrete flux. Indeed, :math:`(R_T)_{kv}A_{kv}` is the volume-integrated rhs, i. e., the sum of fluxes through the faces.

Note that since :math:`\mathcal{F}` is a bilinear form, the definition of the rhs is always globally consistent. Indeed, the total variance
dissipation is :math:`\sum_{k,v}T_{kv}(R_T)_{kv}A_{kv}
=2\sum_{k,v}T_{kv}\partial\mathcal{F}/\partial \tilde T_{kv}=2\mathcal{F}(T,T)`.

In summary, the variational formulation originally proposed for quadrilaterals can easily be extended to triangular meshes. All symmetry properties will be granted if computations
are local on subprisms.

Substituting :math:`{\bf K}` in the form :math:`\mathcal{F}` we get

.. math::
   \mathcal{F} =\sum_{k,c}\sum_p[-K_i\nabla_h \tilde T\cdot\nabla_h T-K_i\nabla_h \tilde
   T\cdot{\bf s}\partial_zT-K_i\partial_z\tilde T{\bf s}\cdot \nabla_h T-(K_d +s^2K_i)\partial_z \tilde
   T\partial_zT]_{kcp}(A_ch_{kc}/6).

The first term does not involve the slope and will not be considered.

Let us start from the third term and compute its contribution to :math:`\partial
\mathcal{F}/\partial \tilde T_{kv}`. The vertical derivative at level :math:`k` (the top surface of layer :math:`k`) is

.. math::
   (\partial_zT)_{kv} = \frac{T_{(k-1)v}-T_{kv}}{Z_{(k-1)v}-Z_{kv}},

and :math:`\nabla_h T` is defined on cell :math:`c`

.. math::
   (\nabla_h T)_{kc} = \sum_{v(c)}{\bf G}_{cv}T_{kv},

Hence it follows for the contribution from layer :math:`k` and element :math:`c`

.. math::
   \frac{\partial\mathcal{F}}{\partial \tilde T_{kv}} :\quad
   \frac{1}{6}A_ch_{kc}\left[\frac{-1}{Z_{k-1}-Z_k}(-K_i{\bf s})^t_{kcv}(\nabla_h T)_{kc} +
   \frac{1}{Z_k-Z_{k+1}} (-K_i{\bf s})^b_{kcv}\cdot(\nabla_h T)_{kc}\right],

.. math::
   \frac{\partial\mathcal{F}}{\partial \tilde T_{(k-1)v}} :\quad
   \frac{1}{6}A_ch_{kc}
   \frac{1}{Z_{k-1}-Z_k}(-K_i{\bf s})^t_{kcv}\cdot(\nabla_h T)_{kc},


.. math::
   \frac{\partial\mathcal{F}}{\partial \tilde T_{(k+1)v}} :\quad
   \frac{1}{6}A_ch_{kc}\frac{-1}{Z_{k}-Z_{k+1}}(-K_i{\bf s})^b_{kcv}\cdot(\nabla_g T)_{kc}.


In the
expressions above, indices :math:`k` and :math:`c` identify the triangular prism, and the index of vertex :math:`v` together with the upper index :math:`t` or :math:`b` identify the subprism (related to :math:`v` and either top or bottom of the full prism). The expression :math:`(K_i{\bf s})^t_{kcv}` means that :math:`K_i` is estimated on level :math:`k` and vertex :math:`v`, and the slope involves the triplet with :math:`\alpha,\beta` at :math:`kv`, the vertical derivatives at :math:`kv` and the horizontal derivatives at :math:`kc`. For :math:`(K_i{\bf s})^b_{kcv}`, the pairs of indices are :math:`(k+1)v,\, kv,\,(k+1)v` and :math:`kc` respectively.

Now, we combine the contributions from the column associated with cell :math:`c`
that enter the rhs of equation on :math:`T_{kv}` (they come from prisms :math:`(k-1)c`, :math:`kc` and :math:`(k + 1)c`)

.. math::
   \frac{\partial\mathcal{F}}{\partial\tilde T_{kv}}:\quad \frac{A_c}{6}\left[
   \frac{h_{kc}}{Z_{k-1}-Z_k}(K_i{\bf s}\cdot\nabla_h T)^t_{kcv}+
   \frac{h_{(k-1)c}}{Z_{k-1}-Z_k}(K_i{\bf s}\cdot\nabla_h T)^b_{(k-1)cv}\right.

.. math::
   \left.-\frac{h_{kc}}{Z_{k}-Z_{k+1}}(K_i{\bf s}\cdot\nabla_h T)^b_{kcv}-
   \frac{h_{(k+1)c}}{Z_{k}-Z_{k+1}}(K_i{\bf s}\cdot\nabla T)^t_{(k+1)cv}\right].

We easily recognize here the fluxes through the upper and lower surfaces of scalar prism :math:`kv` coming from the part shared with prism :math:`kc`. They are thickness-weighed over the cells on both sides. Indeed, :math:`2(Z_{k-1}-Z_k) = h_{kc}+h_{(k-1)c}` for the top surface and similarly for the bottom.

We continue with the
contribution from :math:`-s^2K_i\partial_z \tilde T\partial_zT`.
The contribution to equation at (:math:`kv`) from prisms :math:`(k-1)c`, :math:`kc` and :math:`(k+1)c` may come from the following terms in :math:`\mathcal{F}`

.. math::
   \frac{A_c}{6}\left[(-s^2K_i)^t_{kcv}\frac{\tilde T_{(k-1)v}- \tilde
   T_{kv}}{Z_{k-1}-Z_k}\frac{T_{(k-1)v}-T_{kv}}{Z_{k-1}-Z_k}h_{kc}+\right.

.. math::
   (-s^2K_i)^b_{kcv}\frac{\tilde T_{kv}- \tilde
   T_{(k+1)v}}{Z_{k}-Z_{k+1}}\frac{T_{kv}-T_{(k+1)v}}{Z_{k}-Z_{k+1}}h_{kc}+

.. math::
   (-s^2K_i)^b_{(k-1)cv}\frac{\tilde T_{(k-1)v}-\tilde
   T_{kv}}{Z_{k-1}-Z_k}\frac{T_{(k-1)v}-T_{kv}}{Z_{k-1}-Z_k}h_{(k-1)c}+

.. math::
   \left.(-s^2K_i)^t_{(k+1)cv}\frac{\tilde T_{kv}- \tilde
   T_{(k+1)v}}{Z_{k}-Z_{k+1}}\frac{T_{kv}-T_{(k+1)v}}{Z_{k}-Z_{k+1}}h_{(k+1)c}\right].


Now, performing differentiation with respect to :math:`T_{kv}`, we find

.. math::
   \frac{\partial\mathcal{F}}{\partial \tilde T_{kv}} = \frac{A_c}{6} \left[
   \left( \frac{h_{kc}}{Z_{k-1}-Z_k}
   (s^2K_i ))^t_{kcv}+ \frac{h_{(k-1)c}}{Z_{k-1}-Z_k} (s^2K_i
   ))^b_{(k-1)cv}\right)\frac{ T_{k-1}- T_k}{Z_{k-1}-Z_k}\right.

.. math::
   +\left.\left(-\frac{h_{kc}}{Z_k-Z_{k+1}}
   (s^2K_i))^b_{kcv}-\frac{h_{(k+1)c}}{Z_k-Z_{k+1}}(s^2K_i))^t_{(k+1)cv}\right)\frac{T_k-T_{k+1}}{Z_k-Z_{k+1}}\right].

The result is the standard scheme for the vertical diffusion, but the
estimates of :math:`s^2K_i` are thickness-weighted over contributing layers. The fluxes
through the top and bottom surfaces can conveniently be assembled in a cycle over cells and layers.

We return to the horizontal part in the expression for :math:`\mathcal{F}`. Layer :math:`k` and cell :math:`c` contribute to :math:`\mathcal{F}` as

.. math::
   \frac{A_c}{6}h_{kc}(\sum_{v(c)}{\bf G}_{cv}\tilde
   T_{kv})\cdot\left[\sum_{v(c)}\frac{T_{(k-1)v}-
   T_{kv}}{Z_{k-1}-Z_k}(-K_i{\bf s})^t_{kcv}+\right.

.. math::
   \left.\sum_{v(c)}\frac{T_{kv}- T_{(k+1)v}}{Z_k-Z_{k+1}}(-K_i
   {\bf s})^b_{kcv}\right].

For the contribution into equation :math:`kv` from :math:`\partial \mathcal{F}/\partial\tilde T_{kv}` it is straightforward to prove that it corresponds to the flux of the quantity in the square brackets through the segments bounding the control volume around :math:`v` inside triangle :math:`c`. Indeed, for geometrical reasons :math:`{\bf G}_{cv}` is :math:`{\bf n}_{cv}/h_{cv}` with :math:`{\bf n}_{cv}` the normal to the edge of :math:`c` opposing vertex :math:`v` directed from this vertex (outer for :math:`c`) and :math:`h_{cv}` the height in :math:`c` drawn from :math:`v`. This implies that :math:`A_c{\bf G}_{cv}={\bf n}_{cv}l_{cv}/2`, where :math:`l_{cv}` is the length of the opposing edge. Obviously, for the two segments bounding the control volume :math:`v` inside cell :math:`c` the sum of normal vectors multiplied with the lengths of segments is :math:`{\bf n}_{cv}l_{cv}/2`. Thus, we arrive at flux representation.

Although computations as written are possible, FESOM at present follows a simplified scheme which deals with the slope vector averaged over the prism (instead of considering 6 different slope vectors). The motivation for this step is purely numerical -- it is more computationally efficient and more stable. The associated dianeutral mixing is the subject of study. The implementation of full scheme is delayed.
