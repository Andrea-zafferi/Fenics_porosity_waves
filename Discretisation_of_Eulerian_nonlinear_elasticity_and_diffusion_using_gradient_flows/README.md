# Discretisation of Eulerian nonlinear elasticity and diffusion using gradient flows

**Authors:** Andrea Zafferi, Dirk Peschka

### Overview

This is a collection of nonlinear PDEs for Lagrangian/Eulerian hyperelastic solids and ultimately for a two-phase model describing the travelling porosity wave phenomenon. It covers the implementation shown in 

This directory contains FEniCS-based simulations for elastic-diffusive material with porosity waves used in the paper

A. Zafferi & D. Peschka. *Discretisation of Eulerian nonlinear elasticity and diffusion using gradient flows*, WIAS Preprint ??, (2026) DOI []()


**Abstract:** In this study, we introduce a general energy-based modelling approach for viscous poroelastic materials that feature diffusive
transport in both Lagrangian and Eulerian frames. Our research produces refined weak formulations by using
the reference map concept within the Eulerian configuration. We propose and implement a novel structure-preserving
discretisation strategy, utilising mixed finite element methods. This paper highlights the spatial and temporal numerical
convergence of our methods through a comparative analysis of Lagrangian and Eulerian schemes, thereby proving the
robustness and usability of our approach. Furthermore, in the context of Eulerian multiphase flow, specifically of the
quasi-static Euler-Euler type, our study demonstrates the existence of solitary fluid waves within poroviscoelastic media.
This energy-based approach forms a basis for a deeper understanding of thermodynamical modelling and corresponding
discretisation schemes for coupled poroelasticity, flow, and diffusion.

**Note:** Slight code conversion and documentations were done with AI assistance using GitHub Copilot.

The directory is divided into subdirectiories, each one containing the codes that generated the figures of the corresponding section. Within each subdirectory, there is a *README.md* file explaining the requirements and the use of the contained codes.
