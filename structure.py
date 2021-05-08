import math
import cmath 
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.pyplot as plt
from Gaussian_quadrature import quadrature

class Structure ():
	def __init__ (self,filename, D , x_end, length_damage): 

		df = pd.read_excel(filename, sheet_name ='NENJNS')
		self.nj = int ( df['nj'] )
		# print ("self.nj",self.nj)
		self.ne = int ( df['ne'] )
		# print ("self.ne",self.ne)
		self.ns = int ( df['ns'] )
		# print ("self.ns",self.ns)
		self.lj = int ( df['lj'] )
		self.nf = self.nj * 3 - self.ns
		self.coor = np.array ( pd.read_excel(filename, sheet_name ='joints') )
		self.l = np.array ( [ self.coor[i] - self.coor[i-1] for i in range(1, len(self.coor)) ] )
		df =  pd.read_excel(filename, sheet_name ='Connectivity') 
		self.be = np.array (  df['be']  )
		self.en = np.array (  df['en']  )
		df = pd.read_excel(filename, sheet_name ='Elements')
		self.e = np.array (  df['E']  )
		self.a = np.array (  df['A']  )
		self.density = np.array (  df['Rho']  )
		self.ii =  np.array (  df['I']  )
		df = pd.read_excel(filename, sheet_name ='Fixed')
		self.fix = np.array (  df['Constraints']  )
		self.fix -= 1 		
		# print("fix",self.fix)
		self.added_dofs = []
		middle = self.coor[-1] / 2

		if  D != 0 :
			if ( x_end  < middle ):
				# print ('x_mid_damage + length_damage/2  < middle')
				self.nj += 2
				self.ne += 2
				self.lj += 6
				self.nf += 6
				self.coor = np.insert (self.coor, 1, x_end, axis=None)
				self.coor = np.insert (self.coor, 1, x_end - length_damage, axis=None)
				self.l = np.array ( [ self.coor[i] - self.coor[i-1] for i in range(1, len(self.coor)) ] )
				self.be = [1,2,3,4]
				self.en = [2,3,4,5]
				self.e = np.insert (self.e, 1, (1-D) * self.e[0] , axis = None)
				self.e = np.insert (self.e, 2, self.e[0] , axis = None)
				self.a = np.insert (self.a, 1, self.a[0] , axis = None)
				self.a = np.insert (self.a, 2, self.a[0] , axis = None)
				self.density = np.insert (self.density, 1, self.density[0] , axis = None)
				self.density = np.insert (self.density, 2, self.density[0] , axis = None)
				self.ii = np.insert (self.ii, 1, self.ii[0] , axis = None)
				self.ii = np.insert (self.ii, 2, self.ii[0] , axis = None)
				for i in range ( 0, np.size (self.fix) ):
					if ( self.fix [i] >= 3 ):
						self.fix [i] += 6
				self.added_dofs = [3,4,5,6,7,8]

			if x_end  == middle :
				# print('x_mid_damage + length_damage/2  == middle')
				self.nj += 1
				self.ne += 1
				self.lj += 3
				self.nf += 3
				self.coor = np.insert (self.coor, 1, x_end - length_damage, axis=None)
				self.l = np.array ( [ self.coor[i] - self.coor[i-1] for i in range(1, len(self.coor)) ] )
				self.be = [1,2,3]
				self.en = [2,3,4]
				self.e = np.insert (self.e, 1, (1-D) * self.e[0] , axis = None)
				self.a = np.insert (self.a, 1, self.a[0] , axis = None)
				self.density = np.insert (self.density, 1, self.density[0] , axis = None)
				self.ii = np.insert (self.ii, 1, self.ii[0] , axis = None)
				for i in range ( 0, np.size (self.fix) ):
					if ( self.fix [i] >= 3 ):
						self.fix [i] += 3
				self.added_dofs = [3,4,5]

			if x_end  > middle and x_end - length_damage  < middle:
				# print ('x_mid_damage + length_damage/2  > middle and x_mid_damage - length_damage/2  < middle')
				self.nj += 2
				self.ne += 2
				self.lj += 6
				self.nf += 6
				self.coor = np.insert (self.coor, 2, x_end, axis=None)
				self.coor = np.insert (self.coor, 1, x_end - length_damage, axis=None)
				self.l = np.array ( [ self.coor[i] - self.coor[i-1] for i in range(1, len(self.coor)) ] )
				self.be = [1,2,3,4]
				self.en = [2,3,4,5]
				self.e = np.insert (self.e, 1, (1-D) * self.e[0] , axis = None)
				self.e = np.insert (self.e, 2, (1-D) * self.e[0] , axis = None)
				self.a = np.insert (self.a, 1, self.a[0] , axis = None)
				self.a = np.insert (self.a, 2, self.a[0] , axis = None)
				self.density = np.insert (self.density, 1, self.density[0] , axis = None)
				self.density = np.insert (self.density, 2, self.density[0] , axis = None)
				self.ii = np.insert (self.ii, 1, self.ii[0] , axis = None)
				self.ii = np.insert (self.ii, 2, self.ii[0] , axis = None)
				for i in range ( 0, np.size (self.fix) ):
					if ( self.fix [i] >= 3 ):
						self.fix [i] += 6
				self.added_dofs = [3,4,5,9,10,11]

			if x_end  > middle and x_end - length_damage  == middle:
				# print ('x_mid_damage + length_damage/2  > middle and x_mid_damage - length_damage/2  == middle')
				self.nj += 1
				self.ne += 1
				self.nf += 3
				self.coor = np.insert (self.coor, 2, x_end, axis=None)
				self.l = np.array ( [ self.coor[i] - self.coor[i-1] for i in range(1, len(self.coor)) ] )
				self.be = [1,2,3]
				self.en = [2,3,4]
				self.e = np.insert (self.e, 1, (1-D) * self.e[0] , axis = None)
				self.a = np.insert (self.a, 1, self.a[0] , axis = None)
				self.density = np.insert (self.density, 1, self.density[0] , axis = None)
				self.ii = np.insert (self.ii, 1, self.ii[0] , axis = None)
				for i in range ( 0, np.size (self.fix) ):
					if ( self.fix [i] >= 3 ):
						self.fix [i] += 3
				self.added_dofs = [6,7,8]

			if x_end  > middle and x_end - length_damage  > middle :
				# print ('x_mid_damage + length_damage/2  < middle ')
				self.nj += 2
				self.ne += 2
				self.nf += 6
				self.coor = np.insert (self.coor, 2, x_end, axis=None)
				self.coor = np.insert (self.coor, 2, x_end - length_damage, axis=None)
				self.l = np.array ( [ self.coor[i] - self.coor[i-1] for i in range(1, len(self.coor)) ] )
				self.be = [1,2,3,4]
				self.en = [2,3,4,5]
				self.e = np.insert (self.e, 1, self.e[0] , axis = None)
				self.e = np.insert (self.e, 2, (1-D) *self.e[0] , axis = None)
				self.a = np.insert (self.a, 1, self.a[0] , axis = None)
				self.a = np.insert (self.a, 2, self.a[0] , axis = None)
				self.density = np.insert (self.density, 1, self.density[0] , axis = None)
				self.density = np.insert (self.density, 2, self.density[0] , axis = None)
				self.ii = np.insert (self.ii, 1, self.ii[0] , axis = None)
				self.ii = np.insert (self.ii, 2, self.ii[0] , axis = None)
				for i in range ( 0, np.size (self.fix) ):
					if ( self.fix [i] >= 3 ):
						self.fix [i] += 6
				self.added_dofs = [6,7,8,9,10,11]

			# print (f'number of joints is {self.nj}')
			# print (f'lj is {self.lj}')
			# print (f'number of elements is {self.ne} ')
			# print (f'nf is {self.nf} ')
			# print (f'coor is {self.coor} ')
			# print (f'lenghts are {self.l} ')
			# print (f'E are {self.e} ')
			# print (f'A are {self.a} ')
			# print (f'I are {self.ii} ')
			# print (f'fix is {self.fix}')


	def element ( self, ielem, om ) :
		self.k_spectral_element = np.zeros( (6, 6) )
		k1 = ( self.density [ ielem ]  * self.a [ ielem ] * (om**2) / ( self.e[ielem] * self.ii [ielem] ) ) ** 0.25    
		kl = k1 * self.l [ielem]  
		det = 1 - math.cos ( kl ) * math.cosh (kl)
		# if det == 0 :
		# 	print ("det = 0",det)
		# 	print ("k1=",kl)
		# 	print ("element length=", self.l[ielem] )
		# 	print ("coordinates",self.coor)
		stf = self.e[ ielem ] * self.ii [ ielem ] / ( self.l[ ielem ] ** 3 )
		alfa =  ( math.cos( kl ) * math.sinh( kl ) + math.sin( kl ) * math.cosh( kl ) ) * ( kl ** 3 ) / det
		alfabar = ( math.sin (kl) + math.sinh ( kl ) ) * ( kl ** 3 ) / det    
		beta = ( (-1) * math.cos(kl) * math.sinh(kl) + math.sin(kl) * math.cosh(kl) ) * (kl) / det
		betabar = ( (-1) * math.sin(kl) + math.sinh(kl) ) * (kl) / det
		gama = ( (-1) * math.cos(kl) + math.cosh(kl) ) * ( kl ** 2 ) / det
		gamabar = math.sin(kl) * math.sinh(kl) * (kl**2) / det

		self.k_spectral_element [1,1] = stf * alfa
		self.k_spectral_element [1,2] = stf * gamabar * self.l[ielem]
		self.k_spectral_element [1,4] = stf * (-1) * alfabar
		self.k_spectral_element [1,5] = stf * gama * self.l[ielem]
		self.k_spectral_element [2,2] = stf * beta * ( self.l[ ielem ] ** 2 )
		self.k_spectral_element [2,4] = stf * (-1) * gama * self.l[ielem] #!!!!!! It was written gammabar in doyles paper which is wrong
		self.k_spectral_element [2,5] = stf * betabar * ( self.l[ ielem ] ** 2 )
		self.k_spectral_element [4,4] = stf * alfa
		self.k_spectral_element [4,5] = stf * (-1) * gamabar * self.l[ ielem ] 
		self.k_spectral_element [5,5] = stf * beta * (self.l[ ielem ] ** 2 ) 
		
		self.k_spectral_element [2,1] = self.k_spectral_element[1,2]
		self.k_spectral_element [4,1] = self.k_spectral_element[1,4]
		self.k_spectral_element [5,1] = self.k_spectral_element[1,5]
		self.k_spectral_element [4,2] = self.k_spectral_element[2,4]
		self.k_spectral_element [5,2] = self.k_spectral_element[2,5]
		self.k_spectral_element [5,4] = self.k_spectral_element[4,5]

		#BAR ELEMENT*************************************

		ak = om * ( self.density [ielem] * self.a [ielem] / (self.e[ielem] * self.a[ielem] ) ) ** 0.5
		exp_1 = cmath.exp( (1j) * (-1) * ak * self.l [ielem] )
		exp_2 = cmath.exp( (1j) * ak * self.l[ielem] )  

		self.k_spectral_element [0,0] = self.e[ielem] * self.a[ielem] / self.l[ielem] * (ak * self.l[ielem] / math.sin (ak*self.l[ielem]) ) * math.cos ( ak * self.l[ielem] )
		self.k_spectral_element [0,3] = (-1) * self.e[ielem] * self.a[ielem] / self.l[ielem] * (ak * self.l[ielem] / math.sin (ak*self.l[ielem]) )		             
		self.k_spectral_element [3,0] = self.k_spectral_element [0,3]
		self.k_spectral_element [3,3] = self.k_spectral_element [0,0]

	def k_global (self,ielem):
		if ielem == 0:
			self.k_spectral_global  = np.zeros ( (self.nj*3,self.nj*3) )
		ndofn = 3 
		nnode = 2 
		for inode in range (1,nnode+1):
			if inode == 1 :
				nodei = self.be[ielem]
			if inode == 2 :
				nodei = self.en[ielem]
			for idofn in range (1,ndofn+1):
				nrows = ( nodei - 1 ) * ndofn + idofn - 1
				nrowe = ( inode - 1 ) * ndofn + idofn - 1	
				for jnode in range (1,nnode+1):
					if jnode == 1 :
						nodej = self.be[ielem]
					if jnode == 2 :
						nodej = self.en[ielem]
					for jdofn in range (1,ndofn+1):
						ncols = ( nodej - 1 ) * ndofn + jdofn - 1
						ncole = ( jnode - 1 ) * ndofn + jdofn - 1 
						self.k_spectral_global [nrows,ncols]= self.k_spectral_global [nrows,ncols] + self.k_spectral_element[nrowe,ncole]
		
	def constraint (self) :
		self.k_spectral_global =  np.delete (np.delete (self.k_spectral_global, self.fix, 0), self.fix, 1)

	def element_mass ( self, ielem ) :
		self.mass_element = np.zeros( (6, 6) )
		mm = self.density[ielem] * self.a[ielem] * self.l[ielem] / 420

		self.mass_element [1,1] = mm * 156
		self.mass_element [1,2] = mm * 22 * self.l[ielem]
		self.mass_element [1,4] = mm * 54
		self.mass_element [1,5] = mm * -13 * self.l[ielem]
		self.mass_element [2,2] = mm * 4 * self.l[ielem]**2
		self.mass_element [2,4] = mm * 13 * self.l[ielem]
		self.mass_element [2,5] = mm * -3 * self.l[ielem]**2
		self.mass_element [4,4] = mm * 156
		self.mass_element [4,5] = mm * -22 * self.l[ielem]
		self.mass_element [5,5] = mm * 4 * self.l[ielem]**2
		
		self.mass_element [2,1] = self.mass_element[1,2]
		self.mass_element [4,1] = self.mass_element[1,4]
		self.mass_element [5,1] = self.mass_element[1,5]
		self.mass_element [4,2] = self.mass_element[2,4]
		self.mass_element [5,2] = self.mass_element[2,5]
		self.mass_element [5,4] = self.mass_element[4,5]

		#BAR ELEMENT*************************************
		mm = (1/6) * self.density[ielem]*self.a[ielem]*self.l[ielem]

		self.mass_element [0,0] = mm * 2
		self.mass_element [0,3] = mm * 1
		self.mass_element [3,0] = mm * 1
		self.mass_element [3,3] = mm * 2

	def m_global_matrix (self,ielem):
		if ielem == 0:
			self.m_global  = np.zeros ( (self.nj*3,self.nj*3) )
		ndofn = 3 
		nnode = 2 
		for inode in range (1,nnode+1):
			if inode == 1 :
				nodei = self.be[ielem]
			if inode == 2 :
				nodei = self.en[ielem]
			for idofn in range (1,ndofn+1):
				nrows = ( nodei - 1 ) * ndofn + idofn - 1
				nrowe = ( inode - 1 ) * ndofn + idofn - 1	
				for jnode in range (1,nnode+1):
					if jnode == 1 :
						nodej = self.be[ielem]
					if jnode == 2 :
						nodej = self.en[ielem]
					for jdofn in range (1,ndofn+1):
						ncols = ( nodej - 1 ) * ndofn + jdofn - 1
						ncole = ( jnode - 1 ) * ndofn + jdofn - 1 
						self.m_global [nrows,ncols]= self.m_global [nrows,ncols] + self.mass_element[nrowe,ncole]

	def m_constraint (self) :
		self.m_global =  np.delete (np.delete (self.m_global, self.fix, 0), self.fix, 1)


	def mass_spectral_element ( self, ielem , om,type=2):#,iterations, type
		num_points=10
		iterations = 10
		# type = 1
		if type == 1 :
			#trapozoidal method
			iterations = 10
			kl = om * ( ( self.density [ielem] * self.a[ielem] ) / ( self.e[ielem] * self.a[ielem] ) ) ** 0.5
			# l_bar = kf * self.l[ielem]
			# eta = 2 * kf * ( 1 - mat.cos( l_bar ) * math.cosh (l_bar) )
			n_bar = np.ndarray ( shape = (1,2) , dtype=float)
			self.m_spectral_element = np.ndarray ( shape = (6,6) , dtype=float)
			self.m_spectral_element.fill(0)
			for i in range ( 1, iterations + 1):
				delta_x = self.l[ielem] / iterations
				x = i * delta_x
				n_bar [0,0] =  ( 1 / math.sin ( kl * self.l[ielem] ) ) * math.sin ( kl *( self.l[ielem] - x  ) )
				n_bar [0,1] =  ( 1 / math.sin ( kl  *  self.l[ielem] ) ) * math.sin ( kl * x )
				n_bar_T_n_bar =np.matrix.transpose (n_bar)  @ n_bar 
				self.m_spectral_element [0,0] += n_bar_T_n_bar [0,0] * delta_x 
				self.m_spectral_element [0,3] += n_bar_T_n_bar [0,1] * delta_x 
				self.m_spectral_element [3,0] += n_bar_T_n_bar [1,0] * delta_x 
				self.m_spectral_element [3,3] += n_bar_T_n_bar [1,1] * delta_x 
		else: 
			kl = om * ( ( self.density [ielem] * self.a[ielem] ) / ( self.e[ielem] * self.a[ielem] ) ) ** 0.5
			n_bar = np.ndarray ( shape = (1,2) , dtype=float)
			self.m_spectral_element = np.ndarray ( shape = (6,6) , dtype=float)
			self.m_spectral_element.fill(0)
			points,weights = quadrature(0,self.l[ielem],num_points)
			for i in range ( 0,num_points ):

				x = points[i]
				w = weights[i]

				n_bar [0,0] =  ( 1 / math.sin ( kl * self.l[ielem] ) ) * math.sin ( kl *( self.l[ielem] - x  ) )
				n_bar [0,1] =  ( 1 / math.sin ( kl  *  self.l[ielem] ) ) * math.sin ( kl * x )
				n_bar_T_n_bar =np.matrix.transpose (n_bar)  @ n_bar 

				self.m_spectral_element [0,0] += n_bar_T_n_bar [0,0] * w
				self.m_spectral_element [0,3] += n_bar_T_n_bar [0,1] * w
				self.m_spectral_element [3,0] += n_bar_T_n_bar [1,0] * w
				self.m_spectral_element [3,3] += n_bar_T_n_bar [1,1] * w

			#Jacobian Scaling
			self.m_spectral_element [0,0] *= self.l[ielem] / 2
			self.m_spectral_element [0,3] *= self.l[ielem] / 2
			self.m_spectral_element [3,0] *= self.l[ielem] / 2
			self.m_spectral_element [3,3] *= self.l[ielem] / 2


		self.m_spectral_element [0,0] *= self.density[ielem]  * self.a[ielem] 
		self.m_spectral_element [0,3] *= self.density[ielem]  * self.a[ielem] 
		self.m_spectral_element [3,0] *= self.density[ielem]  * self.a[ielem] 
		self.m_spectral_element [3,3] *= self.density[ielem]  * self.a[ielem] 

		#-----------------------------------------------
		if type ==1:
			n_matrix = np.zeros((4,4))
			kf = ( om**0.5 ) * ( ( self.density [ielem] * self.a[ielem] ) / ( self.e[ielem] * self.ii[ielem] ) ) ** 0.25		
			l_bar = kf * self.l[ielem]
			eta = 2 * kf * ( 1 - math.cos( l_bar ) * math.cosh (l_bar) )
			n_beam = np.ndarray ( shape = (1,4) , dtype=float)
			for i in range ( 1, iterations + 1):
				delta_x = self.l[ielem] / iterations
				x = i * delta_x
				x_bar = kf * x
				n_beam [0,0] =  (1/eta) * kf * ( math.cos(x_bar) - math.cos(l_bar-x_bar)*math.cosh(l_bar) - math.cos(l_bar)*math.cosh(l_bar-x_bar) + math.cosh(x_bar) + math.sin(l_bar-x_bar)*math.sinh(l_bar) - math.sin(l_bar)*math.sinh(l_bar-x_bar) )
				n_beam [0,1] =  (1/eta) * ( -1*math.cosh(l_bar-x_bar)*math.sin(l_bar) +  math.cosh(l_bar)*math.sin(l_bar-x_bar) + math.sin(x_bar) - math.cos(l_bar-x_bar)* math.sinh(l_bar) +math.cos(l_bar)* math.sinh(l_bar-x_bar) +  math.sinh(x_bar)  )
				n_beam [0,2] =  (1/eta) * kf * ( math.cos(l_bar-x_bar) - math.cos(x_bar)*math.cosh(l_bar) - math.cos(l_bar)*math.cosh(x_bar) + math.cosh(l_bar-x_bar) +  math.sin(x_bar)* math.sinh(l_bar) -  math.sin(l_bar)* math.sinh(x_bar)  )
				n_beam [0,3] = -1 * (1/eta) * (-1*math.cosh(x_bar)*math.sin(l_bar) + math.cosh(l_bar)* math.sin(x_bar) +  math.sin(l_bar-x_bar) - math.cos(x_bar)* math.sinh(l_bar) + math.cos(l_bar)* math.sinh(x_bar) +  math.sinh(l_bar-x_bar) )
				n_beam_t_n_beam = np.matrix.transpose (n_beam) @ n_beam 
				n_matrix [0,0] += n_beam_t_n_beam[0,0] * delta_x
				n_matrix [0,1] += n_beam_t_n_beam[0,1] * delta_x
				n_matrix [0,2] += n_beam_t_n_beam[0,2] * delta_x
				n_matrix [0,3] += n_beam_t_n_beam[0,3] * delta_x
				n_matrix [1,0] += n_beam_t_n_beam[1,0] * delta_x
				n_matrix [1,1] += n_beam_t_n_beam[1,1] * delta_x
				n_matrix [1,2] += n_beam_t_n_beam[1,2] * delta_x
				n_matrix [1,3] += n_beam_t_n_beam[1,3] * delta_x
				n_matrix [2,0] += n_beam_t_n_beam[2,0] * delta_x
				n_matrix [2,1] += n_beam_t_n_beam[2,1] * delta_x
				n_matrix [2,2] += n_beam_t_n_beam[2,2] * delta_x
				n_matrix [2,3] += n_beam_t_n_beam[2,3] * delta_x
				n_matrix [3,0] += n_beam_t_n_beam[3,0] * delta_x
				n_matrix [3,1] += n_beam_t_n_beam[3,1] * delta_x
				n_matrix [3,2] += n_beam_t_n_beam[3,2] * delta_x
				n_matrix [3,3] += n_beam_t_n_beam[3,3] * delta_x 

				self.m_spectral_element [1,1] +=  n_beam_t_n_beam[0,0] * delta_x 
				self.m_spectral_element [1,2] +=  n_beam_t_n_beam[0,1] * delta_x 
				self.m_spectral_element [1,4] +=  n_beam_t_n_beam[0,2] * delta_x 
				self.m_spectral_element [1,5] +=  n_beam_t_n_beam[0,3] * delta_x 
				
				self.m_spectral_element [2,1] +=  n_beam_t_n_beam[1,0] * delta_x 
				self.m_spectral_element [2,2] +=  n_beam_t_n_beam[1,1] * delta_x 
				self.m_spectral_element [2,4] +=  n_beam_t_n_beam[1,2] * delta_x 
				self.m_spectral_element [2,5] +=  n_beam_t_n_beam[1,3] * delta_x 
				
				self.m_spectral_element [4,1] +=  n_beam_t_n_beam[2,0] * delta_x 
				self.m_spectral_element [4,2] +=  n_beam_t_n_beam[2,1] * delta_x 
				self.m_spectral_element [4,4] +=  n_beam_t_n_beam[2,2] * delta_x 
				self.m_spectral_element [4,5] +=  n_beam_t_n_beam[2,3] * delta_x 
				
				self.m_spectral_element [5,1] +=  n_beam_t_n_beam[3,0] * delta_x 
				self.m_spectral_element [5,2] +=  n_beam_t_n_beam[3,1] * delta_x 
				self.m_spectral_element [5,4] +=  n_beam_t_n_beam[3,2] * delta_x 
				self.m_spectral_element [5,5] +=  n_beam_t_n_beam[3,3] * delta_x 

		else:
			#Gauss integration
			kf = ( om**0.5 ) * ( ( self.density [ielem] * self.a[ielem] ) / ( self.e[ielem] * self.ii[ielem] ) ) ** 0.25
			l_bar = kf * self.l[ielem]
			eta = 2 * kf * ( 1 - math.cos( l_bar ) * math.cosh (l_bar) )
			n_beam = np.ndarray ( shape = (1,4) , dtype=float)
			for i in range (0,num_points):

				x = points[i]
				x_bar = kf * x
				w = weights[i]

				n_beam [0,0] =  (1/eta) * kf * ( math.cos(x_bar) - math.cos(l_bar-x_bar)*math.cosh(l_bar) - math.cos(l_bar)*math.cosh(l_bar-x_bar) + math.cosh(x_bar) + math.sin(l_bar-x_bar)*math.sinh(l_bar) - math.sin(l_bar)*math.sinh(l_bar-x_bar) )
				n_beam [0,1] =  (1/eta) * ( -1*math.cosh(l_bar-x_bar)*math.sin(l_bar) +  math.cosh(l_bar)*math.sin(l_bar-x_bar) + math.sin(x_bar) - math.cos(l_bar-x_bar)* math.sinh(l_bar) +math.cos(l_bar)* math.sinh(l_bar-x_bar) +  math.sinh(x_bar)  )
				n_beam [0,2] =  (1/eta) * kf * ( math.cos(l_bar-x_bar) - math.cos(x_bar)*math.cosh(l_bar) - math.cos(l_bar)*math.cosh(x_bar) + math.cosh(l_bar-x_bar) +  math.sin(x_bar)* math.sinh(l_bar) -  math.sin(l_bar)* math.sinh(x_bar)  )
				n_beam [0,3] = -1 * (1/eta) * (-1*math.cosh(x_bar)*math.sin(l_bar) + math.cosh(l_bar)* math.sin(x_bar) +  math.sin(l_bar-x_bar) - math.cos(x_bar)* math.sinh(l_bar) + math.cos(l_bar)* math.sinh(x_bar) +  math.sinh(l_bar-x_bar) )

				n_beam_t_n_beam = np.matrix.transpose (n_beam) @ n_beam  

				self.m_spectral_element [1,1] +=  n_beam_t_n_beam[0,0] * w
				self.m_spectral_element [1,2] +=  n_beam_t_n_beam[0,1] * w
				self.m_spectral_element [1,4] +=  n_beam_t_n_beam[0,2] * w
				self.m_spectral_element [1,5] +=  n_beam_t_n_beam[0,3] * w
				
				self.m_spectral_element [2,1] +=  n_beam_t_n_beam[1,0] * w
				self.m_spectral_element [2,2] +=  n_beam_t_n_beam[1,1] * w
				self.m_spectral_element [2,4] +=  n_beam_t_n_beam[1,2] * w
				self.m_spectral_element [2,5] +=  n_beam_t_n_beam[1,3] * w
				
				self.m_spectral_element [4,1] +=  n_beam_t_n_beam[2,0] * w
				self.m_spectral_element [4,2] +=  n_beam_t_n_beam[2,1] * w
				self.m_spectral_element [4,4] +=  n_beam_t_n_beam[2,2] * w
				self.m_spectral_element [4,5] +=  n_beam_t_n_beam[2,3] * w
				
				self.m_spectral_element [5,1] +=  n_beam_t_n_beam[3,0] * w
				self.m_spectral_element [5,2] +=  n_beam_t_n_beam[3,1] * w
				self.m_spectral_element [5,4] +=  n_beam_t_n_beam[3,2] * w
				self.m_spectral_element [5,5] +=  n_beam_t_n_beam[3,3] * w

			#Jacobian Scaling
			self.m_spectral_element [1,1] *= self.l[ielem] / 2 
			self.m_spectral_element [1,2] *= self.l[ielem] / 2 
			self.m_spectral_element [1,4] *= self.l[ielem] / 2 
			self.m_spectral_element [1,5] *= self.l[ielem] / 2 
			self.m_spectral_element [2,1] *= self.l[ielem] / 2 
			self.m_spectral_element [2,2] *= self.l[ielem] / 2 
			self.m_spectral_element [2,4] *= self.l[ielem] / 2 
			self.m_spectral_element [2,5] *= self.l[ielem] / 2 
			self.m_spectral_element [4,1] *= self.l[ielem] / 2 
			self.m_spectral_element [4,2] *= self.l[ielem] / 2 
			self.m_spectral_element [4,4] *= self.l[ielem] / 2 
			self.m_spectral_element [4,5] *= self.l[ielem] / 2 
			self.m_spectral_element [5,1] *= self.l[ielem] / 2 
			self.m_spectral_element [5,2] *= self.l[ielem] / 2 
			self.m_spectral_element [5,4] *= self.l[ielem] / 2 
			self.m_spectral_element [5,5] *= self.l[ielem] / 2 


		self.m_spectral_element [1,1] *=  self.density[ielem] * self.a[ielem]
		self.m_spectral_element [1,2] *=  self.density[ielem] * self.a[ielem]
		self.m_spectral_element [1,4] *=  self.density[ielem] * self.a[ielem]
		self.m_spectral_element [1,5] *=  self.density[ielem] * self.a[ielem]
		
		self.m_spectral_element [2,1] *=  self.density[ielem] * self.a[ielem]
		self.m_spectral_element [2,2] *=  self.density[ielem] * self.a[ielem]
		self.m_spectral_element [2,4] *=  self.density[ielem] * self.a[ielem]
		self.m_spectral_element [2,5] *=  self.density[ielem] * self.a[ielem]

		self.m_spectral_element [4,1] *=  self.density[ielem] * self.a[ielem]
		self.m_spectral_element [4,2] *=  self.density[ielem] * self.a[ielem]
		self.m_spectral_element [4,4] *=  self.density[ielem] * self.a[ielem]
		self.m_spectral_element [4,5] *=  self.density[ielem] * self.a[ielem]

		self.m_spectral_element [5,1] *=  self.density[ielem] * self.a[ielem]
		self.m_spectral_element [5,2] *=  self.density[ielem] * self.a[ielem]
		self.m_spectral_element [5,4] *=  self.density[ielem] * self.a[ielem]
		self.m_spectral_element [5,5] *=  self.density[ielem] * self.a[ielem]
		
	def m_spectral_global_matrix (self,ielem):
		if ielem == 0:
			self.m_spectral_global  = np.zeros ( (self.nj*3,self.nj*3) )
		ndofn = 3 
		nnode = 2 
		for inode in range (1,nnode+1):
			if inode == 1 :
				nodei = self.be[ielem]
			if inode == 2 :
				nodei = self.en[ielem]
			for idofn in range (1,ndofn+1):
				nrows = ( nodei - 1 ) * ndofn + idofn - 1
				nrowe = ( inode - 1 ) * ndofn + idofn - 1	
				for jnode in range (1,nnode+1):
					if jnode == 1 :
						nodej = self.be[ielem]
					if jnode == 2 :
						nodej = self.en[ielem]
					for jdofn in range (1,ndofn+1):
						ncols = ( nodej - 1 ) * ndofn + jdofn - 1
						ncole = ( jnode - 1 ) * ndofn + jdofn - 1 
						self.m_spectral_global [nrows,ncols]= self.m_spectral_global [nrows,ncols] + self.m_spectral_element[nrowe,ncole]

	def m_spectral_constraint (self) :
		self.m_spectral_global =  np.delete (np.delete (self.m_spectral_global, self.fix, 0), self.fix, 1)


	def spectral_bar_element  ( self, ielem, om ) :

		self.k_spectral_bar_element = np.zeros( (2, 2) )

		ak = om * ( self.density [ielem] * self.a [ielem] / (self.e[ielem] * self.a[ielem] ) ) ** 0.5
		exp_1 = cmath.exp( (1j) * (-1) * ak * self.l [ielem] )
		exp_2 = cmath.exp( (1j) * ak * self.l[ielem] )  

		self.k_spectral_bar_element [0,0] = self.e[ielem] * self.a[ielem] / self.l[ielem] * (ak * self.l[ielem] / math.sin (ak*self.l[ielem]) ) * math.cos ( ak * self.l[ielem] )
		
		self.k_spectral_bar_element [0,1] = (-1) * self.e[ielem] * self.a[ielem] / self.l[ielem] * (ak * self.l[ielem] / math.sin (ak*self.l[ielem]) )		             
		
		self.k_spectral_bar_element [1,0] = self.k_spectral_bar_element [0,1]
		
		self.k_spectral_bar_element [1,1] = self.k_spectral_bar_element [0,0]

	def assemble_spectral_bar (self,ielem):
		if ielem == 0:
			self.k_spectral_bar_global  = np.zeros ( (self.nj,self.nj) )
		ndofn = 1
		nnode = 2 
		for inode in range (1,nnode+1):
			if inode == 1 :
				nodei = self.be[ielem]
			if inode == 2 :
				nodei = self.en[ielem]
			for idofn in range (1,ndofn+1):
				nrows = ( nodei - 1 ) * ndofn + idofn - 1
				nrowe = ( inode - 1 ) * ndofn + idofn - 1	
				for jnode in range (1,nnode+1):
					if jnode == 1 :
						nodej = self.be[ielem]
					if jnode == 2 :
						nodej = self.en[ielem]
					for jdofn in range (1,ndofn+1):
						ncols = ( nodej - 1 ) * ndofn + jdofn - 1
						ncole = ( jnode - 1 ) * ndofn + jdofn - 1 
						self.k_spectral_bar_global [nrows,ncols]= self.k_spectral_bar_global [nrows,ncols] + self.k_spectral_bar_element[nrowe,ncole]

	def constraint_bar (self,fix) :
		self.k_spectral_bar_global =  np.delete (np.delete (self.k_spectral_bar_global, fix, 0), fix, 1)

	def spectral_baem_element ( self, ielem, om ) :

		self.k_spectral_beam_element = np.zeros( (4, 4) )
		k1 = ( self.density [ ielem ]  * self.a [ ielem ] * (om**2) / ( self.e[ielem] * self.ii [ielem] ) ) ** 0.25    
		kl = k1 * self.l [ielem]  
		det = 1 - math.cos ( kl ) * math.cosh (kl)

		stf = self.e[ ielem ] * self.ii [ ielem ] / ( self.l[ ielem ] ** 3 )
		alfa =  ( math.cos( kl ) * math.sinh( kl ) + math.sin( kl ) * math.cosh( kl ) ) * ( kl ** 3 ) / det
		alfabar = ( math.sin (kl) + math.sinh ( kl ) ) * ( kl ** 3 ) / det    
		beta = ( (-1) * math.cos(kl) * math.sinh(kl) + math.sin(kl) * math.cosh(kl) ) * (kl) / det
		betabar = ( (-1) * math.sin(kl) + math.sinh(kl) ) * (kl) / det
		gama = ( (-1) * math.cos(kl) + math.cosh(kl) ) * ( kl ** 2 ) / det
		gamabar = math.sin(kl) * math.sinh(kl) * (kl**2) / det

		self.k_spectral_beam_element [0,0] = stf * alfa
		self.k_spectral_beam_element [0,1] = stf * gamabar * self.l[ielem]
		self.k_spectral_beam_element [0,2] = stf * (-1) * alfabar
		self.k_spectral_beam_element [0,3] = stf * gama * self.l[ielem]
		self.k_spectral_beam_element [1,1] = stf * beta * ( self.l[ ielem ] ** 2 )
		self.k_spectral_beam_element [1,2] = stf * (-1) * gama * self.l[ielem] #!!!!!! It was written gammabar in doyles paper which is wrong
		self.k_spectral_beam_element [1,3] = stf * betabar * ( self.l[ ielem ] ** 2 )
		self.k_spectral_beam_element [2,2] = stf * alfa
		self.k_spectral_beam_element [2,3] = stf * (-1) * gamabar * self.l[ ielem ] 
		self.k_spectral_beam_element [3,3] = stf * beta * (self.l[ ielem ] ** 2 ) 
		
		self.k_spectral_beam_element [1,0] = self.k_spectral_beam_element[0,1]
		self.k_spectral_beam_element [2,0] = self.k_spectral_beam_element[0,2]
		self.k_spectral_beam_element [2,1] = self.k_spectral_beam_element[1,2]
		self.k_spectral_beam_element [3,0] = self.k_spectral_beam_element[0,3]
		self.k_spectral_beam_element [3,1] = self.k_spectral_beam_element[1,3]
		self.k_spectral_beam_element [3,2] = self.k_spectral_beam_element[2,3]

	def assemble_spectral_beam (self,ielem):
		if ielem == 0:
			self.k_spectral_beam_global  = np.zeros ( (self.nj*2,self.nj*2) )
		ndofn = 2
		nnode = 2 
		for inode in range (1,nnode+1):
			if inode == 1 :
				nodei = self.be[ielem]
			if inode == 2 :
				nodei = self.en[ielem]
			for idofn in range (1,ndofn+1):
				nrows = ( nodei - 1 ) * ndofn + idofn - 1
				nrowe = ( inode - 1 ) * ndofn + idofn - 1	
				for jnode in range (1,nnode+1):
					if jnode == 1 :
						nodej = self.be[ielem]
					if jnode == 2 :
						nodej = self.en[ielem]
					for jdofn in range (1,ndofn+1):
						ncols = ( nodej - 1 ) * ndofn + jdofn - 1
						ncole = ( jnode - 1 ) * ndofn + jdofn - 1 
						self.k_spectral_beam_global [nrows,ncols]= self.k_spectral_beam_global [nrows,ncols] + self.k_spectral_beam_element[nrowe,ncole]

	def constraint_beam (self,fix) :
		self.k_spectral_beam_global =  np.delete (np.delete (self.k_spectral_beam_global, fix, 0), fix, 1)

	def element_k_fe ( self, ielem ) :
		self.k_fe_element = np.zeros( (6, 6) )
		stf = self.e[ielem] * self.ii[ielem] / ( self.l[ielem] ** 3 )

		self.k_fe_element [1,1] = stf * 12
		self.k_fe_element [1,2] = stf * 6 * self.l[ielem]
		self.k_fe_element [1,4] = stf*(-12)
		self.k_fe_element [1,5] = stf*6*self.l[ielem]

		self.k_fe_element [2,2] = stf*4*self.l[ielem]**2
		self.k_fe_element [2,4] = -1*stf*6*self.l[ielem]
		self.k_fe_element [2,5] = stf*2*self.l[ielem]**2

		self.k_fe_element [4,4] = stf*12
		self.k_fe_element [4,5] = stf*(-6)*self.l[ielem]

		self.k_fe_element [5,5] = stf*4*self.l[ielem]**2
		
		self.k_fe_element [2,1] = self.k_fe_element[1,2]
		self.k_fe_element [4,1] = self.k_fe_element[1,4]
		self.k_fe_element [5,1] = self.k_fe_element[1,5]
		self.k_fe_element [4,2] = self.k_fe_element[2,4]
		self.k_fe_element [5,2] = self.k_fe_element[2,5]
		self.k_fe_element [5,4] = self.k_fe_element[4,5]

		#BAR ELEMENT*************************************

		self.k_fe_element [0,0] =  self.e[ielem]*self.a[ielem]/self.l[ielem]
		self.k_fe_element [0,3] = -self.e[ielem]*self.a[ielem]/self.l[ielem]
		self.k_fe_element [3,0] = -self.e[ielem]*self.a[ielem]/self.l[ielem]
		self.k_fe_element [3,3] = self.e[ielem]*self.a[ielem]/self.l[ielem]

	def k_fe_global_matrix (self,ielem):
		if ielem == 0:
			self.k_fe_global  = np.zeros ( (self.nj*3,self.nj*3) )
		ndofn = 3 
		nnode = 2 
		for inode in range (1,nnode+1):
			if inode == 1 :
				nodei = self.be[ielem]
			if inode == 2 :
				nodei = self.en[ielem]
			for idofn in range (1,ndofn+1):
				nrows = ( nodei - 1 ) * ndofn + idofn - 1
				nrowe = ( inode - 1 ) * ndofn + idofn - 1	
				for jnode in range (1,nnode+1):
					if jnode == 1 :
						nodej = self.be[ielem]
					if jnode == 2 :
						nodej = self.en[ielem]
					for jdofn in range (1,ndofn+1):
						ncols = ( nodej - 1 ) * ndofn + jdofn - 1
						ncole = ( jnode - 1 ) * ndofn + jdofn - 1 
						self.k_fe_global [nrows,ncols]= self.k_fe_global [nrows,ncols] + self.k_fe_element[nrowe,ncole]

	def k_fe_constraint (self) :
		self.k_fe_global =  np.delete (np.delete (self.k_fe_global, self.fix, 0), self.fix, 1)

	def shape_functions (self,ielem,om):

		num_points = int ( self.l[ielem] )
		n_beam = np.empty (shape=(num_points,4))

		kf = ( om**0.5 ) * ( ( self.density [ielem] * self.a[ielem] ) / ( self.e[ielem] * self.ii[ielem] ) ) ** 0.25
		l_bar = kf * self.l[ielem]
		eta = 2 * kf * ( 1 - math.cos( l_bar ) * math.cosh (l_bar) )

		for point in range (0,num_points):

			x = ( point / num_points ) * self.l[ielem]
			x_bar = kf * x

			n_beam [point,0] =  (1/eta) * kf * ( math.cos(x_bar) - math.cos(l_bar-x_bar)*math.cosh(l_bar) - math.cos(l_bar)*math.cosh(l_bar-x_bar) + math.cosh(x_bar) + math.sin(l_bar-x_bar)*math.sinh(l_bar) - math.sin(l_bar)*math.sinh(l_bar-x_bar) )
			n_beam [point,1] = (1/eta) * ( -1*math.cosh(l_bar-x_bar)*math.sin(l_bar) +  math.cosh(l_bar)*math.sin(l_bar-x_bar) + math.sin(x_bar) - math.cos(l_bar-x_bar)* math.sinh(l_bar) +math.cos(l_bar)* math.sinh(l_bar-x_bar) +  math.sinh(x_bar)  )
			n_beam [point,2] = (1/eta) * kf * ( math.cos(l_bar-x_bar) - math.cos(x_bar)*math.cosh(l_bar) - math.cos(l_bar)*math.cosh(x_bar) + math.cosh(l_bar-x_bar) +  math.sin(x_bar)* math.sinh(l_bar) -  math.sin(l_bar)* math.sinh(x_bar)  )
			n_beam [point,3] = -1 * (1/eta) * (-1*math.cosh(x_bar)*math.sin(l_bar) + math.cosh(l_bar)* math.sin(x_bar) +  math.sin(l_bar-x_bar) - math.cos(x_bar)* math.sinh(l_bar) + math.cos(l_bar)* math.sinh(x_bar) +  math.sinh(l_bar-x_bar) )

		return n_beam

	def lumpedmass_element (self,ielem):

		self.lumpedmass_element_matrix = np.zeros( shape = (6,6) )
		mm = self.density[ielem] * self.a[ielem] * self.l[ielem]
		self.lumpedmass_element_matrix [0,0] = mm / 2 
		self.lumpedmass_element_matrix [1,1] = mm / 2 
		self.lumpedmass_element_matrix [3,3] = mm / 2 
		self.lumpedmass_element_matrix [4,4] = mm / 2 

	def lumped_mass_global (self,ielem):
		if ielem == 0:
				self.lumpedmass_global_matrix  = np.zeros ( (self.nj*3,self.nj*3) )
		ndofn = 3 
		nnode = 2 
		for inode in range (1,nnode+1):
			if inode == 1 :
				nodei = self.be[ielem]
			if inode == 2 :
				nodei = self.en[ielem]
			for idofn in range (1,ndofn+1):
				nrows = ( nodei - 1 ) * ndofn + idofn - 1
				nrowe = ( inode - 1 ) * ndofn + idofn - 1	
				for jnode in range (1,nnode+1):
					if jnode == 1 :
						nodej = self.be[ielem]
					if jnode == 2 :
						nodej = self.en[ielem]
					for jdofn in range (1,ndofn+1):
						ncols = ( nodej - 1 ) * ndofn + jdofn - 1
						ncole = ( jnode - 1 ) * ndofn + jdofn - 1 
						self.lumpedmass_global_matrix [nrows,ncols]= self.lumpedmass_global_matrix [nrows,ncols] + self.lumpedmass_element_matrix[nrowe,ncole]

	def lumped_mass_constraint (self) :
		self.lumpedmass_global_matrix =  np.delete (np.delete (self.lumpedmass_global_matrix, self.fix, 0), self.fix, 1)

#---------------Getting Kfe ,Mlumpedm MConsistant
# my_structure = Structure ('10Elements.xlsx',0,0,0)
# for ielem in range (0,my_structure.ne):
# 	#create Kfe
# 	my_structure.element_k_fe(ielem)
# 	my_structure.k_fe_global_matrix(ielem)
# 	#create lumped mass
# 	my_structure.lumpedmass_element(ielem)
# 	my_structure.lumped_mass_global(ielem)
# 	#create consistant mass  
# 	my_structure.element_mass(ielem)
# 	my_structure.m_global_matrix(ielem)

# #applying constraints
# my_structure.k_fe_constraint()
# my_structure.lumped_mass_constraint()
# my_structure.m_constraint()

# #output
# print (f'kfe is :{my_structure.k_fe_global}')
# print (f'm lumped is :{my_structure.lumpedmass_global_matrix}')
# print (f'm consistant is :{my_structure.m_global}')

# #store
# np.savetxt('k.txt', my_structure.k_fe_global, delimiter='\n')
# np.savetxt('m_lumped.txt', my_structure.lumpedmass_global_matrix, delimiter='\n')
# np.savetxt('m_consistant.txt',my_structure.m_global, delimiter='\n')
# k_file = open ('k.txt', 'w')
# for i in range (0,np.shape( my_structure.k_fe_global )[0] ):
# 	for j in range (0,np.shape( my_structure.k_fe_global )[0] ):



# for ielem in range (0,my_structure.ne):
# 	my_structure.element_k_fe (ielem)
# 	my_structure.k_fe_global_matrix (ielem)
# print (my_structure.k_fe_global)
# my_structure.k_fe_constraint()


# print (my_structure.k_fe_global)


# #Finite Element Mass matrix
# my_structure.m_global = 0	
# for ielem in range (0,my_structure.ne):
# 	my_structure.element_mass (ielem)
# 	my_structure.m_global_matrix (ielem)
# my_structure.m_constraint()
# print (my_structure.m_global)

# Spectral Mass matrix
# x = []
# iterations = []
# for iteration in range(1,100):
# 	om = 100
# 	my_structure.mass_spectral_element (1,om,iteration)
# 	x.append(my_structure.m_spectral_element[4,1])
# 	iterations.append(iteration) 
# 	print (f'for iteration is:{iteration}, x is :{x}')

# plt.plot(iterations,x, label = 'convergence')
# plt.legend()
# plt.show()
# Spectral Mass matrix is wrong. Use gauss integration and increase integration points until convergence
# x = []
# oms = []
# for om in range (1000,10000,1):
# 	for ielem in range (0,my_structure.ne):
# 		my_structure.mass_spectral_element (ielem,om)
# 		my_structure.m_spectral_global_matrix (ielem)
# 	my_structure.m_spectral_constraint()
# 	print (f'for om is {om}, m is {my_structure.m_spectral_global[0,0]}')
# 	x.append(my_structure.m_spectral_global[0,0])
# 	oms.append(om)
# plt.plot(oms,x, label = 'convergence')
# plt.legend()
# plt.show()


# for ielem in range (0,my_structure.ne):
# 	my_structure.mass_spectral_element (ielem,om,1)
# 	my_structure.m_spectral_global_matrix (ielem)
# my_structure.m_spectral_constraint()
# print (f'trappzoidal mglobal is {my_structure.m_spectral_global}')
# for ielem in range (0,my_structure.ne):
# 	my_structure.element_mass(ielem)
# 	my_structure.m_global_matrix(ielem)
# my_structure.m_constraint()
# print (f'static m is {my_structure.m_global}')

# temp = np.array ( my_structure.mass_spectral_element ( 1 , 784, 1) )
# print (temp)
# print (temp[0,0:4,0:4])


# ----------Shape functions
# xs=[]
# oms = []
# for om in range(1100,1200):
# 	x =  my_structure.mass_spectral_element (1,om,1)
# 	xs.append(x)
# xs = np.array(xs)
# # print (np.shape(xs))
# for j in range (0,100):
# 	print (f'om is {1100+j}')
# 	n11 = [ xs[j,i,0,0] for i in range(0,100) ]
# 	# print (f'n11 is {n11}')
# 	plt.plot(np.linspace(0,30,num=100),n11)
# 	plt.show()
# ---------------------------------------------

# rangee = np.linspace(1,10000,num=1000)
# om = 1766.4
# # for om in rangee:
# my_structure.mass_spectral_element (1,om,1)

# om = 1500

# my_structure.mass_spectral_element (1,om,1)

# om = 1900

# my_structure.mass_spectral_element (1,om,1)

# Plotting [M] by increasing num of Gauss points
# my_structure = Structure ('Healthy.xlsx',0,0,0)
# num_points = []
# m0_0 = []
# oms=[]
# om = 1764
# for num_point in range(1,16):
# 	for ielem in range (0,my_structure.ne):
# 		my_structure.mass_spectral_element (ielem,om,num_point)
# 		my_structure.m_spectral_global_matrix (ielem)
# 	my_structure.m_spectral_constraint()
# 	num_points.append(num_point)
# 	m0_0.append(my_structure.m_spectral_global[5,5])
# 	print(f'for n={num_point}, m is {my_structure.m_spectral_global}')
# plt.plot(num_points,m0_0)
# plt.xlabel('number of points')
# plt.ylabel('M[5,5]')
# plt.show()

# Plotting [M] by increasing num of Gauss points
# my_structure = Structure ('Healthy.xlsx',0,0,0)
# num_points = []
# m0_0 = []
# oms=[]
# om = 1764
# for num_point in range(1,16):
# 	for ielem in range (0,my_structure.ne):
# 		my_structure.mass_spectral_element (ielem,om,num_point)
# 		my_structure.m_spectral_global_matrix (ielem)
# 	my_structure.m_spectral_constraint()
# 	num_points.append(num_point)
# 	m0_0.append(my_structure.m_spectral_global[5,5])
# 	print(f'for n={num_point}, m is {my_structure.m_spectral_global}')
# plt.plot(num_points,m0_0)
# plt.xlabel('number of points')
# plt.ylabel('M[5,5]')
# plt.show()


#Plotting [M]Spectral by changing om and num_integration points=10
# my_structure = Structure ('Healthy.xlsx',0,0,0)
# num_point = 10
# m0_0 = []
# oms=[]

# for om in range(100,5500,1):
# 	for ielem in range (0,my_structure.ne):
# 		my_structure.mass_spectral_element (ielem,om,num_point)
# 		my_structure.m_spectral_global_matrix (ielem)
# 	my_structure.m_spectral_constraint()
# 	oms.append(om)
# 	m0_0.append(my_structure.m_spectral_global[5,5])
# 	# print(f'for om, m is {my_structure.m_spectral_global}')
# plt.plot(oms,m0_0)
# plt.xlabel('omega')
# plt.ylabel('M[5,5]')
# plt.show()

#Calculating [M] by trapozoidal rule and gauss quadrature
# my_structure = Structure ('Healthy.xlsx',0,0,0)
# om = 1766
# for ielem in range (0,my_structure.ne):
# 	my_structure.mass_spectral_element (ielem,om,type=1)
# 	my_structure.m_spectral_global_matrix (ielem)
# my_structure.m_spectral_constraint()
# print (f'for trapozoidal M is:{my_structure.m_spectral_global}')
# m1 = my_structure.m_spectral_global
# for ielem in range (0,my_structure.ne):
# 	my_structure.mass_spectral_element (ielem,om,type=2)
# 	my_structure.m_spectral_global_matrix (ielem)
# my_structure.m_spectral_constraint()
# print (f'for gauss quadrature M is:{my_structure.m_spectral_global}')
# m2 = my_structure.m_spectral_global
# print (f'ratio is {m1/m2}')

#--------------------Printinh N for beam---------------------------------
# my_structure = Structure ('Healthy.xlsx',0,0,0)
# for om in range (1750,1800,50):
# 	n_values = my_structure.shape_functions (0,om)
# 	x = np.linspace(0,30,30)
# 	# print (x)

# 	fig,axes = plt.subplots(nrows=2,ncols=2)

# 	axes[0,0].plot(x,n_values[:,0])
# 	axes[0,1].plot(x,n_values[:,1])
# 	axes[1,0].plot(x,n_values[:,2])
# 	axes[1,1].plot(x,n_values[:,3])

# 	axes[0,0].set_title(f'N1({om})')
# 	axes[0,1].set_title(f'N2({om})')
# 	axes[1,0].set_title(f'N3({om})')
# 	axes[1,1].set_title(f'N4({om})')
	
	# plt.savefig(f'C://Users//15146//Desktop//FE thesis backup//FE(Thesis)//Thesis New (11 March)//Beam analysis in Python//1Element Simply Supported//{om}.png')

	# plt.tight_layout()
	# plt.show()

#-----------------Create Lumped mass -------------------- 
# for ielem in range(0,my_structure.ne):
# 	my_structure.lumpedmass_element(ielem)
# 	my_structure.lumped_mass_global(ielem)
# my_structure.lumped_mass_constraint()
