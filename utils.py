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
        
def write_state(mesh,q,fname,it=0,t=0):
    data = {
            'step': it,
            'time': t
    }
    write_dictionary(fname + '.json',data)
    
    f=HDF5File(mesh.mpi_comm(),fname  + '.h5', 'w')
    f.write(mesh,"mesh")
    f.write(q,"q",t)
    f.close()
     
def get_space(input_mesh,FE):
    return FunctionSpace(input_mesh,FE)

def read_state(fname,FE):

    data = read_dictionary(fname + '.json')
    it = data['step']
    t  = data['time']
       
    mesh = Mesh()

    f=HDF5File(mesh.mpi_comm(),fname+'.h5', 'r')  
    f.read(mesh,"mesh",False)   
    V = get_space(mesh,FE)
    q = Function(V)
    f.read(q,"q")
    f.close()
    return mesh,q,it,t
