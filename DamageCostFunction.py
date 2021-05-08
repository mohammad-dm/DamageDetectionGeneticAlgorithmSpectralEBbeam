import numpy as np

def mac (vector_1,vector_2): 
    x = ( vector_1.dot(vector_2) )**2
    # print (f'x is {x}')
    y = ( vector_1.dot(vector_1) )
    z = ( vector_2.dot(vector_2) )
    mac = x/(y*z)
    # print (f'mac is :{mac}')
    return mac

def costfunction ( base_freqs, damaged_freqs, base_eigenvectors , damaged_eigenvectors, damaged_added_dofs, base_added_dofs):
	cost = 0
	# print (f'damaged_eigenvectors is damaged_eigenvectors{damaged_eigenvectors}')
	# print (f'damaged_added_dofs is {damaged_added_dofs}')
	damaged_eigenvectors.drop(damaged_added_dofs, inplace = True)
	# print (f'base_eigenvectors is damaged_eigenvectors{base_eigenvectors}')
	# print (f'base_added_dofs is {base_added_dofs}')
	base_eigenvectors = base_eigenvectors.drop(base_added_dofs)
	for i in range ( 0 , np.size(base_freqs) ):
		x = ( ( damaged_freqs [i] - base_freqs [i] ) / (damaged_freqs [i]) ) ** 2
		cost += x	
		#Normalized mode shape difference Using MAC (Mode Accuracy Criterion)
		mode_damaged = np.array( base_eigenvectors [ base_freqs [i] ] )
		mode_base = np.array( damaged_eigenvectors [ damaged_freqs [i] ] )
		macc = mac (mode_damaged,mode_base)
		mode_cost = ( (1-macc**0.5)/(macc**0.5) )**2
		cost += mode_cost
	return cost