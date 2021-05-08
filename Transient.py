from scipy import linalg
import math 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plts

def mass_matrix (my_structure) :
	my_structure.m_global = 0	
	for ielem in range (0,my_structure.ne):
			my_structure.element_mass (ielem)
			my_structure.m_global_matrix (ielem)
	my_structure.m_constraint()

def spectral_mass_matrix (my_structure,om) :
	my_structure.m_global = 0
	for ielem in range (0,my_structure.ne):
		my_structure.mass_spectral_element (ielem, om)
		my_structure.m_spectral_global_matrix (ielem)
	my_structure.m_spectral_constraint()

def lumped_mass_matrix (my_structure) :
	for ielem in range(0,my_structure.ne):
		my_structure.lumpedmass_element(ielem)
		my_structure.lumped_mass_global(ielem)
	my_structure.lumped_mass_constraint()


def alfaa ( steady , l , my_structure , eigenvalues, eigenvectors ):
	y_0 = steady.iloc [ 0 ] . to_numpy ( )
	y_1 = steady.iloc [ 1 ] . to_numpy ( )
	diff = np.zeros ( my_structure.nf ) - y_0
	diff_dot = ( y_1 - y_0 ) / l.samplingInterval
	alfa = np.ndarray ( shape = ( 2 , np.size (eigenvalues ) ) ,dtype = float )
	for i in range ( 0 , np.size (eigenvalues) ):
		om = eigenvalues[i] 
		spectral_mass_matrix ( my_structure , om )
		# mass_matrix (my_structure)
		# lumped_mass_matrix (my_structure)
		m = my_structure.m_spectral_global
		# m = my_structure.m_global
		# m = my_structure.lumpedmass_global_matrix
		# print (f'for om is{om}, m is{m}')
		phi = np.array ( eigenvectors [eigenvalues[i]] )
		phi = np.delete ( phi, my_structure.fix, axis=0 )
		phi_t = np.matrix.transpose (phi)
		alfa_0 = ( phi_t @ m @ diff ) / ( phi_t @ m @ phi )
		# print (f'for mode {i}, ( phi_t @ m @ diff ) is {( phi_t @ m @ diff )} and ( phi_t @ m @ phi ) is {( phi_t @ m @ phi )}')
		alfa_dot_0 = ( phi_t @ m @ diff_dot ) / ( phi_t @ m @ phi )
		alfa [0,i] = alfa_0
		alfa [1,i] = alfa_dot_0
	return alfa

def transient (steady, l , my_structure, eigenvalues, eigenvectors):

	alfa = alfaa ( steady , l , my_structure , eigenvalues, eigenvectors )
	# print (f'alfa is {alfa}')
	transient = pd.DataFrame()
	for i in range ( 0 , np.size ( l.time ) ) :
		t = l.time[i]
		tr = 0
		for j in range (0, np.size (eigenvalues)  ):
			alfa_0_i = alfa[0,j]
			alfa_dot_0_i = alfa[1,j]
			phi = eigenvectors[eigenvalues[j]].values
			alfa_t = alfa_0_i * math.cos( eigenvalues[j] * t ) + alfa_dot_0_i * math.sin ( eigenvalues[j] * t ) / eigenvalues[j]
			tr = tr + alfa_t * phi
		transient [i] = tr
	return transient














