#!/usr/bin/env python
from dolfin import *
from mshr import *
from nonlinear_solver import *

Nt  = 3	
T   = 6
I   = None
rho = 1.0
gamma = 0
constant = 1.e-10
boundary_value= 0.5
#mesh_filename = "P000_lh.hdf5"

#mesh = Mesh("../CerebroMesh/UIKlhpial.xml")
#f = HDF5File(mpi_comm_world(),mesh_filename, 'r')
#f.read(mesh,mesh_filename, False)



sphere = Sphere(dolfin.Point(0.00,0.00, 0.00), 10.0)
mesh = generate_mesh(sphere, 20, "cgal")

V = FunctionSpace(mesh,"CG",1)

def boundary(x,on_boundary):
		return on_boundary
BCS = [DirichletBC(V,Constant(boundary_value),boundary)]




def diffusion_constant(u_bar):
	return constant

def save(du,t,V,dt):
	File("Saves/concentration_%.5d"%t+"s.pvd") << du 	
						
nonlinear_solver(V, Nt,T, I  , rho ,gamma,alpha=diffusion_constant,method="BE" , user_action=save ,BCS=BCS)


