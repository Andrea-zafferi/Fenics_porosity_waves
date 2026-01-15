#!/usr/bin/env python3
"""
Eulerian model for elastic-diffusive material - Porosity Wave Simulation

This script implements a FEniCS-based simulation for elastic-diffusive material
with porosity waves. It loads parameters from JSON files and runs the simulation
until completion.

Usage:
    python porosity_wave.py <parameter_file>
    
Example:
    python porosity_wave.py parameter_sets/pars_diffusive.json
    python porosity_wave.py parameter_sets/pars_flowing.json
"""

import sys
import os
import logging
import pprint
from fenics import *
import numpy as np
from tqdm.auto import tqdm
from utils import *

def energy(u, c, pars):
    """
    Compute the Lagrangian energy functional.
    
    Args:
        u: displacement field
        c: concentration field
        pars: parameter dictionary
        
    Returns:
        Total energy functional
    """
    x = SpatialCoordinate(mesh)
    I = Identity(dim)
    F = inv(I - grad(u))
    J = det(F)
    C = F.T*F / J # isochoric right Cauchy-Green tensor
    
    E     = 0.5*pars['μ']*tr(C - I) 
    Egr   = c * inner( Constant((0,pars['g0'])) , x )
    p     = (c - pars['c0']) - (J - 1)
    Emix  = 0.5*p**2

    Etot = E + Egr + Emix + 0.5*pars['ε']*inner(grad(c),grad(c))
    return Etot

def a(F, c, eta_c, deta_c, pars):
    """
    Bilinear form for diffusion.
    """
    return 0.5 * pars['mob'] * (c**2) * inner(grad(eta_c), grad(deta_c))

def b(F, c, dot_u, dot_c, eta_c, pars):
    """
    Bilinear form for coupling.
    """
    return inner(dot_c, eta_c) - inner(grad(eta_c),c*F*dot_u)

def s(F, c, dot_u, du, pars):
    """
    Bilinear form for viscosity.
    """
    return 0.5 * pars['visc'] * inner(sym(grad(F*dot_u)), sym(grad(F*du)))

def evolve(old_q, tau, pars):
    """
    Solve one time step of the Lagrangian problem.
    
    Args:
        old_q: solution from previous time step
        tau: time step size
        pars: parameter dictionary
        
    Returns:
        Updated solution
    """
    q = Function(Q)
    dq = TestFunction(Q)

    u, c, eta_c = split(q)
    old_u, old_c, old_eta_c = split(old_q)
    du, dc, deta_c = TestFunctions(Q)

    I = Identity(dim)
    old_F = inv(I - grad(old_u))
    old_J = det(old_F)

    H = energy(u, c, pars) * dx
    dot_u = (u - old_u) / tau
    dot_c = (c - old_c) / tau

    Res =   b(old_F, old_c, du, dc, eta_c, pars) * dx - derivative(H, q, dq)
    Res += -s(old_F, old_c, dot_u, du, pars) * dx
    Res +=  b(old_F, old_c, dot_u, dot_c, deta_c, pars) * dx + a(old_F,old_c, eta_c, deta_c, pars) * dx

    q.assign(old_q)
    solve(Res == 0, q, bc, solver_parameters={"newton_solver": {"linear_solver": "mumps",
                                                      "absolute_tolerance": 1e-6,
                                                      "maximum_iterations": 30}})
    return q

def get_dissipation(q, old_q, tau, pars):
    """
    Compute dissipation rates.
    
    Args:
        q: current solution
        old_q: previous solution
        tau: time step size
        pars: parameter dictionary
        
    Returns:
        Tuple of (viscous dissipation, mobility dissipation)
    """
    u, c, eta_c = split(q)
    old_u, old_c, old_eta_c = split(old_q)

    I = Identity(dim)
    old_F = inv(I - grad(old_u))
    old_J = det(old_F)

    dot_u = (u - old_u) / tau
    dot_c = (c - old_c) / tau

    dissi =  assemble(s(old_F, old_c, dot_u, dot_u, pars) * dx)
    mobbi =  assemble(a(old_F, old_c, eta_c, eta_c, pars) * dx)

    return dissi, mobbi

def solve_initial_data(cinit, pars):
    """
    Create proper initial data for the displacement field.
    
    Args:
        cinit: initial concentration field
        pars: parameter dictionary
        
    Returns:
        Initial displacement field
    """
    Qvec = VectorElement('P', triangle, FEu_deg)
    Vtemp = FunctionSpace(mesh, Qvec)

    utemp, dutemp = Function(Vtemp), TestFunction(Vtemp)

    I = Identity(dim)
    Ftemp = inv(I - grad(utemp))

    Etemp = energy(utemp, cinit, pars) * dx
    Etemp += inner(grad(utemp), grad(utemp)) * dx

    bc1 = DirichletBC(Vtemp, Constant((0, 0)), 'on_boundary && (near(x[1],0,1e-3)||near(x[1],2,1e-3))')
    bc2 = DirichletBC(Vtemp.sub(0), Constant(0), 'on_boundary && (near(x[0],0,1e-3) || near(x[0],1,1e-3))')
    bc = [bc1, bc2]

    Restemp = derivative(Etemp, utemp, dutemp)
    solve(Restemp == 0, utemp, bc, solver_parameters={"newton_solver": {"linear_solver": "mumps",
                                                      "absolute_tolerance": 1e-6,
                                                      "maximum_iterations": 30}})
    return utemp

def main():
    """
    Main simulation function.
    """
    global mesh, Q, bc, dim, FEu_deg, FE_deg, tau
    
    # Parse command line arguments
    if len(sys.argv) != 2:
        print("Usage: python porosity_wave.py <parameter_file>")
        print("Example: python porosity_wave.py parameter_sets/pars_diffusive.json")
        sys.exit(1)
    
    parameter_file = sys.argv[1]
    
    # Check if parameter file exists
    if not os.path.exists(parameter_file):
        print(f"Error: Parameter file '{parameter_file}' not found.")
        sys.exit(1)
    
    # Load parameters
    pars = read_dictionary(parameter_file)
    print(f"Loaded parameters from: {parameter_file}")
    pprint.pprint(pars)
    
    # Set up output directory
    fout = pars['fout']
    if not os.path.exists(fout):
        os.makedirs(fout)
    
    # Configure logging
    logging.getLogger('FFC').setLevel(logging.WARNING)
    logging.getLogger('UFL').setLevel(logging.WARNING)
    set_log_active(False) 
    parameters["form_compiler"]["quadrature_degree"] = 5
    
    # -----------------------------
    # Simulation parameters
    # -----------------------------
    print(f"Simulation {pars['basename']} started")
    
    STRUCTURED = pars['STRUCTURED']
    n_steps    = pars['n_steps']
    T          = pars['T']
    tau        = T / n_steps      # Time-step size
    k          = pars['k']        # Mesh resolution factor
    
    # -----------------------------
    # Finite Element Setup
    # -----------------------------
    dim     = 2                   # Spatial dimension
    FEu_deg = 2                   # Degree for displacement
    FE_deg  = 1                   # Degree for concentration and forces
    
    FEu     = VectorElement('CG', triangle, FEu_deg)
    FEc     = FiniteElement('CG', triangle, FE_deg)
    FEmu    = FiniteElement('CG', triangle, FE_deg)
    
    # -----------------------------
    # Mesh generation
    # -----------------------------
    if STRUCTURED:
        mesh = RectangleMesh(Point(0, 0), Point(1, 2), int(k), int(2 * k), diagonal="left/right")
    else:
        domain = Polygon([Point(0, 0), Point(1, 0), Point(1, 2), Point(0, 2)])
        mesh = generate_mesh(domain, int(k))
    
    # Function space for displacement, concentration, and force
    Q = FunctionSpace(mesh, MixedElement([FEu, FEc, FEmu]))
    
    # -----------------------------
    # Boundary conditions
    # -----------------------------
    bc1 = DirichletBC(Q.sub(0), Constant((0, 0)), 'on_boundary && (near(x[1], 0, 1e-3) || near(x[1], 2, 1e-3))')
    bc2 = DirichletBC(Q.sub(0).sub(0), Constant(0), 'on_boundary && (near(x[0], 0, 1e-3) || near(x[0], 1, 1e-3))')
    bc  = [bc1, bc2]
    
    # -----------------------------
    # Initial concentration field
    # -----------------------------
    cexp = Expression('c0 + mag * exp(-alpha * ((pow(x[0]-0.5, 2)) + (pow(x[1]-0.5, 2))))',
                         degree=4, t=0.0,
                         c0=pars['c0'], alpha=pars['α'], mag=pars['mag'])
    cinit = project(cexp, FunctionSpace(mesh, 'CG', 1))
    
    # -----------------------------
    # Initial displacement solution
    # -----------------------------
    print('Computing initial displacement...')
    u_init = solve_initial_data(cinit, pars)
    print('Initial displacement done.')
    
    # -----------------------------
    # Initialize and assign solution
    # -----------------------------
    old_q = Function(Q)
    assign(old_q.sub(0), u_init)
    assign(old_q.sub(1), project(cinit, FunctionSpace(mesh, 'CG', 1)))
    
    # -----------------------------
    # Initial time step
    # -----------------------------
    t   = 0.0
    dt0 = 0.5e-3
    q = evolve(old_q, dt0, pars)
    old_q.assign(q)
    
    write_state(mesh, old_q, f"{pars['fout']}stateP1_{int(k)}_{n_steps}_0", 0, t)
    
    # -----------------------------
    # Time evolution loop
    # -----------------------------
    for i in tqdm(range(n_steps), desc="Solving PDE", unit="step"):
        q = evolve(old_q, tau, pars)
        old_q.assign(q)
        t += tau
        write_state(mesh, old_q, f"{pars['fout']}stateP1_{int(k)}_{n_steps}_{i+1}", i+1, t)
    
    print("Simulation finished.")

if __name__ == "__main__":
    main()
