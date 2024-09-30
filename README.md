# Nonlinear Eulerian elasticity using the inverse map

**Description:** This is a collection of nonlinear PDEs for Lagrangian/Eulerian hyperelastic solids and ultimately for a two-phase model describing the travelling porosity wave phenomenon. 

**Authors:** Andrea Zafferi, Dirk Peschka

## Stationary Lagrangian setting

Nonlinear elasicity is based on a flow map $\bar{\chi}:\bar{\Omega}\to\mathbb{R}^d$, where $\bar{\Omega}\subset\mathbb{R}^d$ is a reference domain and 

```math
\Omega = \bar{\chi}(\bar{\Omega})=\{x=\bar{\chi}(\bar{x}):\bar{x}\in\bar{\Omega}\},
```

denotes the deformed domain. For simplicity we will assume $\Omega=\bar{\Omega}$ throughout this work. Then, for a nonlinear hyperelastic problem with deformation gradient

```math
\bar{\boldsymbol{F}}=\bar{\nabla}\bar{\chi},
```

the observed deformation emerges from the minimization of an energy

```math
\bar{\mathscr{F}}(\bar{q}) = \int_{\bar{\Omega}} W_\mathrm{elast}(\bar{\boldsymbol{F}})\,\mathrm{d}\bar{x},
```

over $\bar{q}=\bar{\chi}$ from a suitable function space.

## Stationary Eulerian setting

Assuming $\bar{\chi}$ in invertible, then for any given function $\bar{f}$ defined on $\bar{\Omega}$ we can introduce a corresponding function $f$ defined on $\Omega$ by

```math
f(x)=\bar{f}(\bar{\chi}^{-1}(x)).
```

If we want to rewrite the energy minimization above in terms of Eulerian functions, then the flow map itself is not useful since $\chi=\bar{\chi}(\bar{\chi}^{-1}(x))=x$ is not a viable unknown function defined on $\Omega$. Instead we can use the return map (inverse map)

```math
\alpha(x) = \bar{\chi}^{-1}(x),
```
which allows us to express

```math
\boldsymbol{F} = (\nabla\alpha)^{-1}(x),
```

and using $q=\alpha$ formulate a corresponding energy minimization of 

```math
\mathscr{F}(q) = \int_{\Omega} \frac{1}{J}W_\mathrm{elast}(\boldsymbol{F})\,\mathrm{d}x,
```

where $J=\det(\boldsymbol{F})$.
