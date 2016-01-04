#!/usr/bin/env python
from dolfin import *
from mshr import *
"""
	PDE :
		rho u_t = div( alpha(u) grad(u)  ) + f 	  in	Omega x (0 ,T]
			u = I 				  on    dOmega x (t=0)
			du/dn = 0 			  on 	dOmega  

	Weak: 
		(rho u_t ,v)  =  (alpha(u)grad(u),grad(v)) +  (f,v) , 

	
"""
	


def nonlinear_solver(function_space, Nt,T, I  , rho ,gamma,alpha =None ,Dalpha=None, f=None,  method="BE" , user_action=None ,tol=10**(-13), k_max =1000, omega=1.0 ): 
	"""
	==============	 ======================================================
	Argument                Explanation
	==============   ======================================================
	function_space	 The constructed FunctionSpace provided mesh in dolfin. 
	Nt 	 	 The number of termporal discretization points.
	T		 The length of the temporal domain.
	I		 The initial condition at the initial time.
	rho		 Constant, rho in R.
	gamma		 The switch between Picard iteration (0) and Newtons method(1).
	alpha 		 The diffusion coefficient, can be nonlinear.
	Dalpha 		 The derivitve of the diffusion coefficient, used in Nwetons method.
	f 		 The source term in the PDE.
	method 		 Selcetion of the temporal discretization, 
			 Options : BE, 2BS ,CN ,CN2
	user_action      User specific task to at each time step.	
	tol 		 The convergence criteria for Newtons and Pircard.
	k_max 		 The maximum number of iterations to obtain convergence in
			 Newton or Picard.
	omega 		 The relaxation parameter
	"""
	import numpy
	set_log_level(30)
	def const(u):
		return 1.0
	
	time , dt = numpy.linspace(0,T,Nt+1,retstep=True)

	if alpha==None:
		alpha = const
	if Dalpha==None:
		gamma = 0.0
		Dalpha = const
	if f==None:
		f   = Expression("0.0",t=0,dt=dt)



	V= function_space

	u = TrialFunction(V) 
	v = TestFunction(V)

	u_bar = Function(V)
	u_1 =Function(V)

	u_2 =Function(V)
	
	u_1.interpolate(I)

	fn = Function(V)
	f1 = Function(V)


	if method=="CN" or method=="2BS" :		
		a = Constant(rho)*inner(u,v)*dx  + Constant(dt/2)*inner( alpha(u_bar)*grad(u),grad(v) )*dx  + Constant(gamma)*Constant(dt/2)*inner( Dalpha(u_bar)*grad(u_bar)*u,grad(v) )*dx
		L = Constant(rho)*inner(u_1,v)*dx -Constant(dt/2)*inner( alpha(u_1)*grad(u_1),grad(v) )*dx   +\
						 Constant(dt/2)*inner(fn,v)*dx + Constant(dt/2)*inner(f1,v)*dx   + Constant(gamma)*Constant(dt/2)*inner(u_bar*Dalpha(u_bar)*grad(u_bar),grad(v) )*dx
	elif method=="CN2":
		a = Constant(rho)*inner(u,v)*dx  + Constant(dt/2)*inner( alpha(Constant(0.5)*(u_bar+u_1) )*grad(u),grad(v) )*dx +\
									Constant(gamma)*Constant(dt/4)*inner( Dalpha( 0.5*(u_bar+u_1) )*grad(u_bar)*u,grad(v) )*dx
		L = Constant(rho)*inner(u_1,v)*dx -  Constant(dt/2)*inner( alpha( Constant(0.5)*(u_bar+u_1) )*grad(u_1),grad(v) )*dx  + Constant(dt/2)*inner((fn+f1),v)*dx +\
									Constant(gamma)*Constant(dt/4)*inner( Dalpha( 0.5*(u_bar+u_1) )*grad(u_bar)*u_bar,grad(v) )*dx

	elif method=="BE":
		a = Constant(rho)*inner(u,v)*dx  + Constant(dt)*inner( alpha(u_bar)*grad(u),grad(v) )*dx  + Constant(gamma)*Constant(dt)*inner(Dalpha(u_bar)*grad(u_bar)*u   , grad(v))*dx
		L = Constant(rho)*inner(u_1,v)*dx  + Constant(dt)*inner(fn,v)*dx  			  + Constant(gamma)*Constant(dt)*inner(u_bar*Dalpha(u_bar)*grad(u_bar),grad(v))*dx
	
	else : #BE is default
		a = Constant(rho)*inner(u,v)*dx  + Constant(dt)*inner( alpha(u_bar)*grad(u),grad(v) )*dx  + Constant(gamma)*Constant(dt)*inner(Dalpha(u_bar)*grad(u_bar)*u   , grad(v))*dx
		L = Constant(rho)*inner(u_1,v)*dx  + Constant(dt)*inner(fn,v)*dx  			  + Constant(gamma)*Constant(dt)*inner(u_bar*Dalpha(u_bar)*grad(u_bar),grad(v))*dx

	du =Function(V)
	


	for n in time[1::]:
		u_bar.vector()[:] = u_1.vector()[:] 

		f1.interpolate(f)
		f.t= n
		fn.interpolate(f)

		
		k=0
		eps= 1.0
		while eps > tol and k< k_max :
			A = assemble(a)
			b = assemble(L)
			solve( A , du.vector(), b )  
			eps =numpy.linalg.norm(u_bar.vector().array()-du.vector().array(), ord=numpy.Inf)
			u_bar.vector()[:]  = omega*du.vector()[:] +(1-omega)*u_bar.vector()[:]
			k+=1
			
			
		if user_action!=None:
			user_action(du,n,V,dt)	
		
		u_2.vector()[:] = u_1.vector()
		u_1.vector()[:] = du.vector()
		
		if method=="2BS":
			a = Constant(rho*1.5)*inner(u,v)*dx  + Constant(dt)*inner( alpha(u_bar)*grad(u),grad(v) )*dx +Constant(gamma)*Constant(dt)*inner(u*Dalpha(u_bar)*grad(u_bar),grad(v))*dx   
			L = Constant(rho*2.0)*inner(u_1,v)*dx-Constant(rho/2.0)*inner(u_2,v)*dx+Constant(dt)*inner(fn,v)*dx +Constant(gamma)*Constant(dt)*inner(Dalpha(u_bar)*grad(u_bar)*u_bar,grad(v))*dx
		
	return u_1	
		
	
			

	
