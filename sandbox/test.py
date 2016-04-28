from dolfin import *

#default values
No_refinements=0 
dt_val=0.1 
#hack to read input arguments 
import sys
for s in sys.argv[1:]:
  exec(s)

print "dt_val ", dt_val, " No_ref ", No_refinements


def check_or_download(filename):
  print filename 
  import os
  if os.path.isfile(filename): 
    print "skipping download"
  else:
    print "downloading ", filename
    import urllib
    urllib.urlretrieve("http://folk.uio.no/kent-and/meshes/%s" % filename, filename)
  

def bounday(x, on_bounday): 
  return on_bounday

check_or_download("pial_mesh.xdmf")
check_or_download("pial_mesh.h5")
mesh = Mesh("pial_mesh.xdmf")
V = FunctionSpace(mesh, "Lagrange", 1)
 
for i in range(No_refinements): mesh = refine(mesh) 

u = TrialFunction(V)
v = TestFunction(V)

U = Function(V)   # current U
U_ = Function(V)  # previous U 

c_val = 1.0 
Tend = 0.2 
t = 0.0
c = Constant(c_val)
dt = Constant(dt_val)

a = u*v*dx + c*dt*inner(grad(u), grad(v))*dx 
L = U_*v*dx 

bc = DirichletBC(V, Constant(1.0), bounday)

A = assemble(a)

U_file = File("U_ref%d_dt=%e.xdmf" % (No_refinements, dt_val))

while t < Tend-dt_val/2:
  t += dt_val 
  bc.apply(A)
  b = assemble(L)
  bc.apply(b)

  solve(A, U.vector(), b, "gmres", "amg")
  
  U_.assign(U)

  print "t ", t

  U_file << U 
  


