from scipy import linalg
import math 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def steady_state ( my_structure, l  ):
	pi = 3.14159265359
	load_real = l.fft_load.real
	load_imag = l.fft_load.imag
	alldofs = np.arange(start=1, stop=my_structure.nf + 1, step=1)
	real_response = pd.DataFrame ()
	imag_response = pd.DataFrame ()
	fr = np.array ( my_structure.nf )
	fi = np.array ( my_structure.nf )
	for i in range (1, int(l.n//2 + 1) + 1 ) : 
		my_structure.k_spectral_global = 0
		# om = l.samplingFrequency * l.n / i
		om = 2 * pi * l.samplingFrequency / l.n * i
		for ielem in range (0,my_structure.ne):
			my_structure.element (ielem, om)
			my_structure.k_global (ielem)
		#appying boundary conditions
		my_structure.constraint()
		fr = l.load_vecs_real [i]
		fi = l.load_vecs_imag [i]
		solver( my_structure, fr, i , om, real_response, l.freqs)
		solver( my_structure, fi, i , om, imag_response, l.freqs)
	# print("imag_response:")
	# print(imag_response)
	# print("real_response:")
	# print(real_response)
	#read a specified degree of freedom from the dataframe used in solver
	# print("nyquist")
	# print(nyquist)

	u_steady = pd.DataFrame()
	u_hat = pd.DataFrame()
	for i in range (0,my_structure.nf):
		user_dof = i
		uhat_real = real_response.iloc[user_dof].to_numpy() 
		uhat_imag = imag_response.iloc[user_dof].to_numpy() 
		uhat =uhat_real + 1j * uhat_imag
		# print("steady response")
		u = np.fft.irfft( uhat )
		# print ( u )
		u_steady [i] = u
		u_hat [i] = uhat
	#add zero to fix degrees of freedom
	# zeros = np.zeros ( shape = ( np.size(l.time) ) ) 
	# for i in range (0,np.size(my_structure.fix)):
	# 	u_steady.insert ( int (my_structure.fix[i]) , "my_structure.fix[i]", "zeros" )
	return u_steady #remove uhat later


	
	# fig, ax = plt.subplots()
	# ax.plot(l.time, u)
	# ax.set(xlabel='time (s)', ylabel='u (in)',
	# title='About as simple as it gets, folks')
	# ax.grid()
	# plt.show()
	# print( np.max(u))















	# maxx = np.max(u)
	# return maxx



	# uhat_pos_real = real_response.loc[user_dof][ 1 : nyquist ]
	# uhat_neg_real = real_response.loc[user_dof][ nyquist  : 2*nyquist - 1 ]
	# uhat_pos_imag = imag_response.loc[user_dof][ 1 : nyquist ]
	# uhat_neg_imag = -1 * imag_response.loc[user_dof][ nyquist  : 2*nyquist - 1 ]
	# # print ("uhat pos real")
	# # print (uhat_pos_real)
	# # print ("uhat neg real")
	# # print (uhat_neg_real)
	# # print ("uhat pos image")
	# # print (uhat_pos_imag)
	# # print ("uhat neg imag")
	# # print (uhat_neg_imag)
	# uhat_real = np.concatenate ( (uhat_pos_real ,uhat_neg_real), axis = 0)
	# uhat_imag = np.concatenate ( (uhat_pos_imag ,uhat_neg_imag), axis = 0)
	# uhat =  uhat_real + 1j*uhat_imag  
	# ut = np.fft.ifft(uhat).real /l.n
	# print(ut)




	# uhat = np.array( 2 * l.n )
	# u_omega = np.array (my_structure.nf)
	# cnt = 0
	# for i in range ( 0, l.n ):
	# 	u_omega = real_response[ l.freqs[i] ]
	# 	end = cnt + my_structure.nf - 1
	# 	print(cnt)
	# 	print(end)
	# 	uhat [ cnt : end ] = u_omega
	# 	cnt += my_structure.nf
	# print('u_omega')
	# print(u_omega)


	#perform ifft to get the response of that degree of freedom

def solver (my_structure, f, i, om, df, freqs):
	#solving
	d = linalg.solve(my_structure.k_spectral_global, f)
	df [i] = d
	# print("d(",om,")=",d)