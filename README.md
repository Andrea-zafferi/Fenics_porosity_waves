# Fenics_porosity_waves
# Collection of Fenics simulations for Lagrangian/Eulerian hyperelastic solids and ultimately for a two-phase model describing the travelling porosity wave phenomenon. 

# Extension to concentration
## Lagrangian formulation

We consider here an extension of the previous model/code to include species diffusion effects. Thus we extend the space state by introducing the *Lagrangian* concentration variable $\bar{c}:\bar{\Omega}\to\mathbb{R}^i$, where, for the moment, $n=1$. We denote with $\bar{q}=(\bar{u},\bar{c})$ the Lagrangian state vector and with $\bar{\mathcal{Q}}$ the state space. The free energy of the system is updated with the inclusion of  *chemical* and a *mixed* energy density:
$$\begin{equation}
\bar{\mathcal{H}}(\bar{q}):= \int_{\bar{\Omega}} \bar{H}_{\rm el}(\bar{F}) + \bar{H}_{\rm ch}(\bar{J},\bar{c}) + \bar{H}_{\rm mix}(\bar{u},\bar{c})\,{\rm d}\bar{x}\,.
\end{equation}$$
The elatic energy denisty $\bar{H}_{\rm el}$ is chosen as before, while for the chemical and mix energy densities we have
$$
\begin{align}
\bar{H}_{\rm ch}(\bar{c},\bar{J})&:= \bar{c}(\log \frac{\bar{c}}{c^*\bar{J}} - 1)\,,
\\
\bar{H}_{\rm mix}(\bar{u}, \bar{c})&:= \bar{f}(\bar{\varrho}+\gamma\bar{c})(\bar{u}_d + \bar{x}_d)\,,
\end{align}
$$
Besides the Darcy potential $\bar{\mathcal{R}}_{\rm Darcy}(\bar{u})$ presented previously we consider the diffusive $H^{-1}$-norm dissipation:
$$\begin{equation}
\bar{\mathcal{R}}_{\rm diff}(\bar{q},\bar{\eta}_{\bar{c}}):=\int_{\bar{\Omega}}\tfrac12 \bar{F}^{-\top}\bar{\nabla}\bar{\eta}_{\bar{c}}\cdot \bar{D}(\bar{c})\bar{F}^{-\top}\bar{\nabla}\bar{\eta}_{\bar{c}}\bar{J}\,{\rm d}\bar{x}\,,
\end{equation}$$
where $\bar{\eta}_{\bar{c}}$ satisfies the equation
$$
\begin{equation}
-\bar{{\rm div}}({\rm Cof}\bar{F}^{\top}\bar{D}(\bar{c})\bar{F}^{-\top}\bar{\nabla}\bar{\eta}_{\bar{c}}) = \dot{\bar{c}}\,,
\end{equation}
$$
with $\bar{D}(\bar{c})$ being the diffusive matrix.
Finally, we formulate the system highlighting the saddle point structure of the equations for all $(\bar{v}_{\bar{u}}, \bar{v}_{\bar{c}}, \bar{\xi}_{\bar{c}})\in\bar{\mathcal{Q}}$:
$$
\begin{align}
-\left\langle\partial_{\dot{\bar{u}}}\bar{\mathcal{R}}_{\rm Darcy}(\bar{u}, \dot{\bar{u}}), \bar{v}_{\bar{u}}\right\rangle &= \left\langle{\rm D}_{\bar{u}}\bar{\mathcal{H}}(\bar{u},\bar{c}),\bar{v}_{\bar{u}}\right\rangle\,,
\\
-\left\langle \bar{\xi}_{\bar{c}}, \partial_{\dot{\bar{c}}}\bar{\mathcal{R}}_{\rm diff}^*(\bar{q}, \bar{\eta}_{\dot{\bar{c}}})\right\rangle -\left\langle\bar{\xi}_{\bar{c}}, \dot{\bar{c}}\right\rangle &= 0\,,
\\
\left\langle\bar{\eta}_{\dot{\bar{c}}}, \bar{v}_{\bar{c}}\right\rangle &=\left\langle{\rm D}_{\bar{c}}\bar{\mathcal{H}}(\bar{u}, \bar{c}), \bar{v}_{\bar{c}}\right\rangle\,.
\end{align}
$$

## Eulerian formulation

As done in the previous example the domain $\bar{\Omega}$ is transformed through the flow map $\chi$ into the actual domain $\Omega$. The new state vector is $q=(u_{\alpha}, c)^{\top}$ where $u_\alpha$ was defined in the previous section and $c \circ \chi = \frac{\bar{c}}{\bar{J}}$. The free energy of the system is transformed accordingly
\begin{equation}
  \mathcal{H}(q):=\int_{\Omega} H_{\rm el}(u_{\alpha}) + {H}_{\rm ch}({c}) + {H}_{\rm mix}(u_{\alpha},{c})\,{\rm d}{x}\,,
\end{equation}
where the elastic energy $H_{\rm el}$ as well as the Darcy dissipation are defined as in the previous section and
$$
\begin{align}
{H}_{\rm ch}(c)&:= {c}(\log \frac{c}{c^*} - 1)\,,
\\
{H}_{\rm mix}(u_{\alpha}, c)&:= f_0(\varrho+\gamma c)x_d\,.
\end{align}
$$
The diffusive dissipation reads now
$$\begin{equation}
{\mathcal{R}}_{\rm diff}({q},{\eta}_{{c}}):=\int_{{\Omega}}\tfrac12 {\nabla}{\eta}_{{c}}\cdot {D}({c}){\nabla}{\eta}_{{c}}\,{\rm d}{x}\,,
\end{equation}$$
where ${\eta}_{{c}}$ satisfies the equation
$$
\begin{equation}
-{{\rm div}}({D}({c}){\nabla}{\eta}_{{c}}) = \dot{{c}}\,,
\end{equation}
$$
The PDEs system finally reads:
$$
\begin{align}
-\left\langle\partial_{\dot{u}}{\mathcal{R}}_{\rm Darcy}({u}, \dot{{u}}), {v}_{{u}}\right\rangle + \left\langle \eta_{\dot{c}},\mathrm{div}(c F \dot{v}_{u}) \right\rangle&= \left\langle{\rm D}_{{u}}{\mathcal{H}}({u},{c}),{v}_{{u}}\right\rangle\,,
\\
\left\langle {\xi}_{{c}}, \partial_{\dot{{c}}}{\mathcal{R}}_{\rm diff}^*({q}, {\eta}_{\dot{{c}}})\right\rangle +\left\langle{\xi}_{{c}}, \dot{{c}}\right\rangle + \left\langle{\xi}_{{c}}, \mathrm{div}(c F \dot{u}_{\alpha}) \right\rangle &= 0\,,
\\
\left\langle{\eta}_{\dot{{c}}}, {v}_{{c}}\right\rangle &=\left\langle{\rm D}_{{c}}{\mathcal{H}}({u}, {c}), {v}_{{c}}\right\rangle\,.
\end{align}
$$
where, in the second equation, we have artificially introduced the convective derivative for $c$.

## TODO

1.   Check that the Lagrangian concentration is exponentially decaying in the steady-state case, i.e., ${\rm D}\bar{\mathcal{H}}=0$ and compare with numerical results in the case of (almost) inelastic material ($\lambda, \mu >>1$) this should go in the paper. Is it possible also for the displacement?
4.  Error analysis: plot errors (H^1 norms) for different mesh- and timestep sizes. Compare with analytical solution of 1.
6. Modify the system to get porosity waves? (See pic) and try to search in the literature for anisotropic mobility in porous media, possibly strain dependent.
