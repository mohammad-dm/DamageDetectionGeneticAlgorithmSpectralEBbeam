import scipy
import scipy.linalg 
import numpy as np
import math
from scipy import linalg
import pandas as pd
from structure import Structure

#Read
#because when this module was used in genetice module, it was returnin infinity values for modes when searching so I changed not to give modes
# If you want to use the modes to get the response do these:
# 1) in  w_w_eigenvalues comment return (df)  and uncomment return (np.array (df.columns.values)) 
# 2) in w_w_eigenvectors_search in the fist line change flag='freq' to flag='modes'


def w_w_eigenvalues (my_structure,num_modes) :
	beam_eigenvalues = w_w_eigenvalues_search ( my_structure, num_modes )
	return (beam_eigenvalues) 
	#return (df) 


def w_w_eigenvectors (my_structure, eigenvalues,num_modes) :
	df = w_w_eigenvectors_search (my_structure ,eigenvalues, num_modes)
	return df

def w_w_eigenvalues_search ( my_structure , num_modes) :
	# finding fixed dofs for bar and beam 
	fixed_flexural = my_structure.fix.tolist()
	fixed_axial = []
	for dof in fixed_flexural :
		if ( dof % 3 == 0 ):
			fixed_flexural.pop(dof)
			fixed_axial.append(dof)
	for i, dof in enumerate(fixed_flexural):
		fixed_flexural [i] -= int (dof/3) + 1

	
	# # # input number of frequencies to find
	# num_freq = 4
	# eigenvalues_bar = np.empty((num_freq))
	# om_u = 0
	# step = 10
	# eps = 0.01
	# num_iteration = 0
	# for i in range (1,num_freq+1):
	# 	k=0
	# 	om = om_u
	# 	num_iteration = 0
	# 	while  num_iteration == 0 :
	# 		k=k+1
	# 		om_l = om
	# 		om_u = om_u + step
	# 		error=1
	# 		j=0
	# 		while error > eps or  j > i  :
	# 				if (num_iteration==0) :
	# 					om =om_u
	# 				else:
	# 					om = 0.5 * ( om_u + om_l )
	# 				k_matrix (my_structure, om, 'bar', fixed_axial)
	# 				j0 = j_0 (my_structure, om, 'bar')
	# 				sign = upper_triangular ( my_structure.k_spectral_bar_global )
	# 				my_structure.k_spectral_bar_global = 0
	# 				j = j0 + sign
	# 				if j >= i  :
	# 					om_u = om
	# 					num_iteration = num_iteration + 1
	# 				if  j < i :
	# 					om_l = om
	# 				error = om_u - om_l
	# 	eigenvalues_bar [ i - 1 ] = om
	# 	om_l = om_u
	
	#Redo for beam
	num_freq = num_modes
	eigenvalues_beam = np.empty((num_freq))
	om_u = 0
	step = 1000
	eps = 0.01
	num_iteration = 0
	for i in range (1,num_freq+1):
		k=0
		om = om_u
		num_iteration = 0
		while  num_iteration == 0 :
			k=k+1
			om_l = om
			om_u = om_u + step
			error=1
			j=0
			while error > eps or  j > i  :
					if (num_iteration==0) :
						om =om_u
	
					else:
						om = 0.5 * ( om_u + om_l )
					k_matrix (my_structure, om, 'beam', fixed_flexural)
					j0 = j_0 (my_structure, om, 'beam')
					sign = upper_triangular ( my_structure.k_spectral_beam_global )

					my_structure.k_spectral_beam_global = 0
					j = j0 + sign
					if j >= i  :
						om_u = om
						num_iteration = num_iteration + 1
					if  j < i :
						om_l = om
					error = om_u - om_l
		eigenvalues_beam [ i - 1 ] = om
		om_l = om_u
	# num_freq = 7
	# eigenvalues = np.sort ( np.concatenate ( (eigenvalues_bar,eigenvalues_beam)  ) )[0:num_freq]
	return eigenvalues_beam#,eigenvalues_bar



def k_matrix (my_structure, om, flag, fix):
	if flag == 'beam':
		my_structure.k_spectral_beam_global = 0
		for ielem in range (0,my_structure.ne):
			my_structure.spectral_baem_element(ielem, om)
			my_structure.assemble_spectral_beam (ielem)
		my_structure.constraint_beam (fix)

	if flag == 'bar':
		my_structure.k_spectral_bar_global = 0
		for ielem in range (0,my_structure.ne):
			my_structure.spectral_bar_element(ielem, om)
			my_structure.assemble_spectral_bar (ielem)
		my_structure.constraint_bar (fix)
def j_0 (my_structure, om, flag):
	j0 = 0
	if flag == 'bar' :
		for ielem in range(0,my_structure.ne):
			bar_wavenumber = om * ( (my_structure.density[ielem]*my_structure.a[ielem]/(my_structure.e[ielem]*my_structure.a[ielem] ) ) ** 0.5 ) * my_structure.l [ielem] 
			ja = int ( bar_wavenumber / np.pi )
			j0 += ja
		return j0

	if flag == 'beam' :
		j0 = 0
		for ielem in range(0,my_structure.ne):
			beam_wavenumber =my_structure.l[ielem] * ( my_structure.density[ielem] * my_structure.a[ielem] * (om**2) / ( my_structure.e[ielem] * my_structure.ii[ielem] ) ) ** 0.25 
			i = int ( beam_wavenumber / np.pi )
			r = 1 - math.cos (beam_wavenumber) * math.cosh (beam_wavenumber)
			if r >= 0.0 :
				sign = 1
			else:
				sign = -1
			jb = i - ( 1 - sign * (-1) ** i  ) / 2
			j0 += jb
		return j0

def upper_triangular ( kk ) :
	#make upper triangular and count
	sign = 0
	matrix = np.copy (kk)
	nf = np.shape (kk)[0]
	for k in range ( 0 , nf ):
		for j in range ( k+1 , nf ) :
			# if matrix[k][k] == 0 :
			# 		print ("matrix[k][k] = 0")
			c = matrix[j][k] / matrix[k][k]
			for p in range ( 0 , nf ):
				matrix [j][p] = matrix[j][p] - c * matrix [k][p]
	max = np.amax (matrix)
	for i in range (0,nf):
		if matrix[i][i] < 0:
				sign = sign + 1
	return sign

def w_w_eigenvectors_search (my_structure, eigenvalues, num_modes) :

	eigenvalues_beam = eigenvalues
	df_axial_modes = pd.DataFrame ()
	df_flexural_modes =pd.DataFrame( )
	#looking for indices of fixed axial dofs
	fixed_flexural = my_structure.fix.tolist()
	fixed_axial = [] 
	for dof in fixed_flexural :
		if ( dof % 3 == 0 ):
			fixed_flexural.pop(dof)
			fixed_axial.append(dof)
	for i, dof in enumerate(fixed_flexural):
		fixed_flexural [i] -= int (dof/3) + 1

	all_axials = [*range(0,my_structure.nj)]
	cnt=0
	for dof in (all_axials):
		if dof == fixed_axial[cnt] :
			all_axials.pop(dof)
	free_axials = all_axials


	all_dofs = [*range(0,my_structure.nj*3)]
	all_axials = list (np.arange (0,my_structure.nj*3,step = 3))
	all_flexural = []
	for i in range (0,my_structure.nj*3):
		if ( i % 3 != 0):
			all_flexural.append(i)

	#go throught axial eigenvalues
	# for i, om in enumerate(eigenvalues_bar):
	# 	#build k(omega_n)
	# 	my_structure.k_spectral_bar_global = 0
	# 	for ielem in range (0,my_structure.ne):
	# 		my_structure.spectral_bar_element(ielem, om)
	# 		my_structure.assemble_spectral_bar (ielem)
	# 	my_structure.constraint_bar (fixed_axial)

	# 	#calculate mode
	# 	last_dof = np.shape (my_structure.k_spectral_bar_global)[0] - 1
	# 	f = (-1) * my_structure.k_spectral_bar_global [:,-1]
	# 	f = np.delete (f,[last_dof],0)
	# 	my_structure.k_spectral_bar_global = np.delete (np.delete (my_structure.k_spectral_bar_global, [last_dof], 0), [last_dof], 1)

	# 	d = linalg.solve (my_structure.k_spectral_bar_global, f)

	# 	#add removed dof (last one)
	# 	d = np.append (d,[1])
	# 	#add fixed axial dof as 0
	# 	for j in fixed_axial:
	# 		d = np.insert (d , j, 0 )

	# 	#Adding flexural dofs as 0
	# 	for j in all_flexural:
	# 		d = np.insert (d , j, 0 )


	# 	#Store mode in database
	# 	d = pd.DataFrame(d)
	# 	df_axial_modes = pd.concat ([df_axial_modes,d], axis = 1) 
		


	#Go through flexural dofs

	for i, om in enumerate(eigenvalues_beam):
		#build k(omega_n)
		my_structure.k_spectral_beam_global = 0
		for ielem in range (0,my_structure.ne):
			my_structure.spectral_baem_element(ielem, om)
			my_structure.assemble_spectral_beam (ielem)
		my_structure.constraint_beam (fixed_flexural)

		#calculate mode
		last_dof = np.shape (my_structure.k_spectral_beam_global)[0] - 1
		f = (-1) * my_structure.k_spectral_beam_global [:,-1]
		f = np.delete (f,[last_dof],0)
		my_structure.k_spectral_beam_global = np.delete (np.delete (my_structure.k_spectral_beam_global, [last_dof], 0), [last_dof], 1)
		
		d = linalg.solve (my_structure.k_spectral_beam_global, f)

		#add removed dof (last one)
		d = np.append (d,[1])
		#add fixed axial dof as 0
		for j in fixed_flexural:
			d = np.insert (d , j, 0 )

		#Adding flexural dofs as 0
		for j in all_axials:
			d = np.insert (d , j, 0 )
		#Store mode in database
		d = pd.DataFrame(d)
		df_flexural_modes = pd.concat ([df_flexural_modes,d], axis = 1) 

	num_freqs = num_modes
	eigenvectors = df_flexural_modes
	eigenvectors.columns = eigenvalues
	return eigenvectors



# my_structure = Structure ("Healthy.xlsx" , 0, 0, 0)
# w_w_eigenvectors ( my_structure )
	
# my_structure = Structure ('Healthy.xlsx',0,0,0)
# eigenvalues = w_w_eigenvalues (my_structure)
# eigenvectors =w_w_eigenvectors (my_structure, eigenvalues)
# print (f'eigenvalues are: {eigenvalues}')
# print (f'eigenvectors are :{eigenvectors}')