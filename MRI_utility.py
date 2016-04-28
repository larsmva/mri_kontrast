#!/usr/bin/env python
# -*- coding: utf-8 -*-
from dolfin import *

"""
DICOM header : http://dicomlookup.com/lookup.asp?sw=Ttable&q=C.8-4


Comments: 

	There exsist different cut-off concentration for different tissue/fluid.

	The cut-off concentration is defined as 
		
			(math.log( (TR*r1+TE*r2)/(TE*r2) ) -TR*R01)/(TR*r1) --> Check

	and is the value for the signal derivative of the conctration is zero.
	Thus the maximum signal intensity possible is obtained at the cut-off 
	concentration, at larger concentration the signal intensities decays.
	


"""
def mric_solver(k1,k2,G,I ,max_num=3000 , tol=1.0e-9, cut_off=True): # TO DO cut_off of I
	"""
	Description:
		    Solver for F(c) =  ..... 

		    The derivative of F(c) is modified 
	Comments:
		There are different formulas for different scan sequences. This mehtod is 
		for spin echo (SE). Type of scan sequence can be read from dicom header.
	==============	 ==============================================================
	Argument                Explanation
	==============	 ==============================================================
	k1		 Defined as TR*r1
	k2		 defind as TE*r2
	G		 defined exp( -TR*R1)
	I		 The relative intensities 

	"""
	import numpy as np
	c_max = -np.log( k2/(G*k1+G*k2)) /k1 
	c_ = c_max*((I-1.)/I )  # initial guess


	def F(c):
		return (1.-G*np.exp(-k1*c))*np.exp(-k2*c) -I*(1.-G) 
	def H(c):
		return - F(c)/(k2-(k1+k2)*G*np.exp(-k1*c))  # ~ F(c)/F'(c)
	def RI(c):
		return (1.-G*np.exp(-k1*c))*np.exp(-k2*c)/(1.-G) 

	# CUT_OFF
	if cut_off==True:
		max_I = RI(c_max) 
		I[I>max_I] = max_I

	num_iter=0
	e = 1.0
	while e > tol and num_iter < max_num :
		c = c_ - H(c_)
		num_iter+=1
		e = abs(c_ - c).max()  # ensures that all values must converge
		c_=c
		if num_iter > max_num:
			raise ValueError('Convergence error')

	return c


def compute_concentration_SE(SI0,SI1,TR,TE,r1,r2,R01,R02,cut_off=True):
		import numpy 
		import math
		"""
		Comments:
			There are different formulas for different scan sequences. This mehtod is 
			for spin echo (SE). Type of scan sequence can be read from dicom header.
		==============	 ==============================================================
		Argument                Explanation
		==============	 ==============================================================
		SI0		 Signal itensiites for MRI without contrast agent. ARRAY OR FUNCTION ?
		SI1		 Signal itensiites for MRI with contrast agent
		TR		 Repettiton time, can be read from header of dicom file. * dim= s
		TE		 Time Echo,, can be read from header of dicom file. *  dim =s
		R01 		 The inverted T1 relaxation time, without contrast agent. dim=s^-1
		Scrap(R02)       The inverted T2 relaxation time, without contrast agent. dim=s^-1
		r1		 The  relaxivity constants for T1, dim = l/(mmol*s)       
		r2 	         The  relaxivity constants for T2, dim = l/(mmol*s)
				*=  Assumed to be the same for SI1 and SI0.

		Formula:
			
			SI = Constant*(1-exp(-TR(R01+r1*C))*exp(-TE(R02+r2*C))
			with C as concenration which is >=0. For SI0 is the value 
			obtianed with C=0.
			
		
		"""
	
		R = SI1/SI0
		R[(R==float('INF')) | (R < 1)]=1 


		C = mric_solver(TR*r1,TE*r2,math.exp(-TR*R01),R)
		
		return C






def get_header_info(path2dicom):
	"""
	Reads the listed info from the dicom-file:
		TE  = echo time
		TR  = Repetition time
		magnetric strength
		scanning sequence  
	""" 
	import dicom # REQUIRE DICOM
	ds = dicom.read_file(path2dicom)
	TE = float(ds.EchoTime)/1000 		#ms->s
	TR = float(ds.RepetitionTime)/1000      #ms->s
	magnetic_strength = float(ds.MagneticFieldStrength)
	scanning_sequence = str(ds.ScanningSequence)		
	
	if scanning_sequence!='SE':
		pass	# Implement error or different methods


	r1,r2,R01,R02 = values_from_database(magnetic_strength)
	
	return TR,TE,r1,r2,R01,R02

def values_from_database(magnetic_strength,tissue="white"): 
	# r1 and r2 taken from MICAD : Molecular Imaging and Contrast Agent Database  DIMENSION ARE l/(mmol s)
	# R01 and R02 taken from : T1, T2 relaxation and magnetization transfer in tissue at 3T , Stanisz, Greg J. et al. * 1/s
	
	v = {"white":{1.5:[4.7,6.8,1.13,13.9],3.0:[ 3.6 , 6.3,0.92, 14.5]},"grey":{1.5:[4.7,6.8,0.89,10.5],3.0:[3.6,6.3,0.55,10.1] } }

	#r1     r2       R01  	R02	
	return v[tissue][magnetic_strength]



def compute_relative_intensities_SE(C,TR,TE,r1,r2,R01,R02):
		import numpy
		return compute_intensities_SE(C,TR,TE,r1,r2,R01,R02,1)/compute_intensities_SE(0,TR,TE,r1,r2,R01,R02,1)
def compute_intensities_SE(C,TR,TE,r1,r2,R01,R02,constant):
		import numpy
		return constant*(1-numpy.exp(-TR*(R01+r1*C)))*numpy.exp(-TE*(R02+r2*C))


def get_MRI_values(path2mri_file,function_space,mesh, offset=None):
		"""
		==============	 ==============================================================
		Argument                Explanation
		==============	 ==============================================================
		path2_mri_file	The path to the MRI file with extension mgz. typically orig.mgz
		function_space	The function space of the mesh.		
		mesh 	 	The mesh of the brain	 
		function	Instead of a return value, update the function inside. 	
		offset 		The translation of the origin of the mesh , need to be np.array
		"""
		import nibabel
		from nibabel.affines import apply_affine
		import numpy.linalg as npl
		import numpy as np

		img= nibabel.load(path2mri_file) 
		inv_aff = npl.inv ( img.get_header().get_vox2ras_tkr() )
		data = img.get_data()
		if offset==None :
			offset=np.array([0,0,0])

		xyz = function_space.dofmap().tabulate_all_coordinates(mesh).reshape((function_space.dim(),-1)) - offset

		i,j,k = apply_affine(inv_aff,xyz).T
	
		# 4.7 = 5 and not 4
		i= map(round,i) 
		j= map(round,j) 
		k= map(round,k) 
	
		return np.array(map(data.item,i,j,k),dtype=float)






