
from fenics import *

class WholeBoundary(SubDomain):
  def inside(self, x, on_boundary):
    return on_boundary #and (x[1] < 1.0 - tol)

mesh = Mesh("left_cerebro.xml.gz")
V = FunctionSpace(mesh, "Lagrange", 1)

u = TrialFunction(V)
v = TestFunction(V)

U = Function(V)
U_prev = Function(V)

D = Constant(1.0)
dt = Constant(0.001)
f = Constant(0)

a = D*inner(grad(u), grad(v))*dx 
bc = DirichletBC(V, Constant(1), WholeBoundary())
L = dt*f*v*dx +U_prev*v*dx   

A = assemble(a) 
bc.apply(A)

t = 0; T = 0.002  
U_file = File("u.pvd")
while t<=T: 
  t += dt 
  b = assemble(L)
  bc.apply(b)
  
  solve(A, U.vector(), b, "gmres", "amg")

  U_prev.vector()[:] = U.vector()[:]
  print "time ", t 

  U_file  << U 

  

  
  




