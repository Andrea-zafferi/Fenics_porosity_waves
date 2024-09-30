# Nonlinear Eulerian elasticity using the inverse map

**Description:** This is a collection of nonlinear PDEs for Lagrangian/Eulerian hyperelastic solids and ultimately for a two-phase model describing the travelling porosity wave phenomenon. 

**Authors:** Andrea Zafferi, Dirk Peschka

## Lagrangian formulation

We consider here an extension of the previous model/code to include species diffusion effects. Thus we extend the space state by introducing the *Lagrangian* concentration variable $\bar{c}:\bar{\Omega}\to\mathbb{R}^i$, where, for the moment, $n=1$. We denote with $\bar{q}=(\bar{u},\bar{c})$ the Lagrangian state vector and with $\bar{\mathcal{Q}}$ the state space. The free energy of the system is updated with the inclusion of  *chemical* and a *mixed* energy density:

```math
\bar{\mathcal{H}}(\bar{q}):= \int_{\bar{\Omega}} \bar{H}_{\rm el}(\bar{F}) + \bar{H}_{\rm ch}(\bar{J},\bar{c}) + \bar{H}_{\rm mix}(\bar{u},\bar{c})\,{\rm d}\bar{x}\,.
```

