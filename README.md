# Nonlinear Eulerian elasticity using the inverse map

**Description:** This is a collection of nonlinear PDEs for Lagrangian/Eulerian hyperelastic solids and ultimately for a two-phase model describing the travelling porosity wave phenomenon. 

**Authors:** Andrea Zafferi, Dirk Peschka

## Setting

Nonlinear elasicity is based on a flow map $\chi:\bar{\Omega}\to\mathbb{R}^d$, where $\bar{\Omega}\subset\mathbb{R}^d$ is a reference domain and 

```math
\Omega = \chi(\bar{\Omega})=\{x=\chi(\bar{x}):\bar{x}\in\bar{\Omega}\},
```

denotes the deformed domain. For simplicity we will assume $\Omega=\bar{\Omega}$ throughout this work. Then, for a nonlinear hyperelastic problem with deformation gradient

```math
\boldsymbol{F}=\nabla\chi
```

the observed deformation emerges from the minimization of an energy

```math
\mathscr{F}(q) = \int_{\bar{\Omega}} W_\mathrm{elast}(x,\boldsymbol{F})\,\mathrm{d}\bar{x}.
```

over $q=\chi$ from a suitable function space.
