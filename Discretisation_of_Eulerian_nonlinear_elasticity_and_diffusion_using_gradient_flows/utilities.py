import json
from fenics import *
import numpy as np

def read_dictionary(fname):
    with open(fname) as json_file:
        pars = json.load(json_file)
    return pars

def write_dictionary(fname,mydict):
    with open(fname, 'w') as fp:
        json.dump(mydict, fp,sort_keys=True, indent=4)
        
#def output_mesh(mesh,fname):
#    ff2=HDF5File(mesh.mpi_comm(),fname  + '.h5', 'w')
#    ff2.write(mesh,"mesh")   
#    ff2.close()

def output_state(mesh,q,fname,it=0,t=0):
    ff1 = open(fname + '.time','w')
    ff1.write(str(it)+"  "+str(t)+'\n')
    ff1.close()
    ff2=HDF5File(mesh.mpi_comm(),fname  + '.h5', 'w')
    ff2.write(mesh,"mesh")
    ff2.write(q,"q",t)
    ff2.close()
     
def get_space(input_mesh):
    Qvec        = VectorElement("CG", input_mesh.ufl_cell(), 1) # Function space for solid alpha
    Qscal       = FiniteElement("CG", input_mesh.ufl_cell(), 1) # Function space for Lagrange multiplier
    return FunctionSpace(input_mesh, MixedElement([Qvec,Qvec,Qscal,Qscal]))

def read_state(fname,FE):
    f = open(fname+'.time','r')
    fin = f.read().split(' ')
    
    it = int(fin[0])
    t  = float(fin[2])
    
    f.close()
    
    mesh = Mesh()
    f=HDF5File(mesh.mpi_comm(),fname+'.h5', 'r')  
    
    f.read(mesh,"mesh",False)   
    V = FunctionSpace(mesh,FE)
    q = Function(V)
    f.read(q,"q")
    f.close()
    
    return mesh,q,it,t


def read_state_mesh(fname,default_functionspace=None):

    f = open(fname+'.time','r')
    fin = f.read().split(' ')
    
    it = int(fin[0])
    t  = float(fin[2])
    
    f.close()
    
    mesh = Mesh()
    f=HDF5File(mesh.mpi_comm(),fname+'.h5', 'r')  
    
    if default_functionspace is not None:
        V = default_functionspace
    else:
        f.read(mesh,"mesh",False)   
        V = get_space(mesh)
        
    q = Function(V)
    f.read(q,"q")
    f.close()
    return mesh,q,it,t

#def read_mesh(fname):
#    mesh = Mesh()
#    f=HDF5File(mesh.mpi_comm(),fname+'.h5', 'r')
#    f.read(mesh,"mesh",False)   
#    f.close()
#    return mesh

def mesh_gen(x0,y0,k,l):
    mesh   = RectangleMesh(Point(0,0),Point(l,2*l),k,2*k,diagonal="left/right")
    return mesh

def refine_circle(x0,y0,k):
    # mesh   = RectangleMesh(Point(0,-40),Point(20,40),8,4*8,diagonal="left/right")
    mesh   = RectangleMesh(Point(0,0),Point(20,40),k,4*k,diagonal="left/right")
    fac = 4/5
    RR = 4/(fac**4)
    for i in range(4):
        cm = MeshFunction("bool", mesh, 2)
        cm.set_all(False)
    
        for cell in cells(mesh):
            p = cell.midpoint()
            if ((p[0]-x0)**2 + (p[1]-y0)**2) < RR**2:
                cm[cell] = True
        RR *= fac
        mesh = refine(mesh, cm)
    return mesh

def refine_mesh(q0,mesh0):
    V0 = FunctionSpace(mesh0,'CG',1)
    phi = project(q0.sub(2),V0)
    max_index = np.argmax(phi.vector()[:])
    x= V0.tabulate_dof_coordinates()
    coords = x[max_index]
    print(" --> refine with max at :",coords[0],coords[1])
    mesh1 = refine_circle(coords[0],coords[1]+2.0)
    V1 = get_space(mesh1)
    q1 = project(q0,V1)
    return q1,mesh1
