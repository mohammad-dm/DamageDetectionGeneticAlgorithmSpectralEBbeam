from Input import Load
from structure import Structure 
from Solver import steady_state
import numpy as np
from Analytical import analytical
import matplotlib.pyplot as plt
from W_W_EigenValues import w_w_eigenvalues , w_w_eigenvectors
from Transient import transient
import pandas as pd
from DamageCostFunction import costfunction
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits import mplot3d


# -------------------------------------------------------Creating cost function for different modes 
# num_modes = 3
# base_model = Structure ("Healthy.xlsx",0,0,0)

# ds = np.array([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
# locations = np.linspace (6, 59, 9)
# length = 5

# eigenvalues_base_model =  w_w_eigenvalues (base_model, num_modes) 

# df = pd.DataFrame ()
# df_mode = []
# print (df)
# for mode in range(0,num_modes) :
# 	for d in ds :
# 		for location in locations :	
# 			damaged_structure = Structure ("Healthy.xlsx",d,location,length)
# 			eigenvalues_damaged_model =  w_w_eigenvalues (damaged_structure, num_modes) 
# 			cost = 0
# 			for i, omega in enumerate (eigenvalues_damaged_model) :
# 				cost = ( ( omega - eigenvalues_base_model[i] ) / eigenvalues_base_model[i] ) ** 2
# 				df_mode.append(cost)
# 			print (df_mode)
# 			df = df.append (df_mode)
# 			df_mode = []
# print (df)


# print (df_mode)






#-------------------------------------------------------------------------------
# def f(x, y):
#     return np.sin(np.sqrt(x ** 2 + y ** 2))

# x = np.linspace(-6, 6, 30)
# y = np.linspace(-6, 6, 30)

# X, Y = np.meshgrid(x, y)
# Z = f(X, Y)
# print (x)
# ax = plt.axes(projection='3d')
# ax.plot_surface(X, Y, Z, rstride=1, cstride=1,
#                 cmap='viridis', edgecolor='none')
# ax.set_title('surface');

# ax.view_init(60, 35)

# results_time = pd.DataFrame()
# results_freq = pd.DataFrame()

# base_model = Structure ("Healthy.xlsx",0,0,0)
# l = Load (1, "Healthy.xlsx", base_model  )

# results_time['Time'] = l.time

# results_time['Load'] = 100*np.sin(2*np.pi*1*l.time)

# fft_l =  np.fft.rfft (results_time['Load'].values)
# results_freq ["Frequencies"] = l.freqs[0:33]
# results_freq ['Load'] = abs(fft_l)
# print (results_freq)

# # results_freq ['Load'] = abs(fft_l)



# steady,uhat = steady_state (base_model,l)
# user_dof = 3
# u_steady = steady [user_dof-1]
# uhat = abs ( uhat [user_dof-1] )
# results_freq ["Response in Freq"] = uhat
# results_time["Response in Time"] = u_steady

# print (results_time)
# print (results_freq)

# fig,axes = plt.subplots(nrows=2,ncols=2)
# axes[0,0].plot(results_time["Time"].values,results_time["Load"].values)
# axes[0,1].plot(results_freq["Frequencies"].values,results_freq["Load"].values)
# axes[1,0].plot(results_freq["Frequencies"].values,results_freq["Response in Freq"].values)
# axes[1,1].plot(results_time["Time"],results_time["Response in Time"])

# axes[0,0].set_title(f'Load in Time Domain')
# axes[0,1].set_title(f'Fourier Transform of Load')
# axes[1,0].set_title(f'Response in Frequency Domain')
# axes[1,1].set_title(f'Response in Time Domain')

# axes[0,0].set_xlabel("t(s)")
# axes[0,1].set_xlabel("Frequency(Hz)")
# axes[1,0].set_xlabel("Frequency(Hz)")
# axes[1,1].set_xlabel("t(s)")

# axes[0,0].set_ylabel("P(t)")
# axes[0,1].set_ylabel("P(f)")
# axes[1,0].set_ylabel("D(f)")
# axes[1,1].set_ylabel("D(t)")
# plt.tight_layout()
# plt.savefig(f'C://Users//15146//Desktop//FE thesis backup//FE(Thesis)//Thesis New (11 March)//Beam analysis in Python//Doyle_Calculation.png')
# plt.show()

#---------------------Moving element------------------------------------------------------
num_modes = 3
base_model = Structure ("Healthy.xlsx",0,0,0)
ds = np.array([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
locations = np.linspace (6, 59, 54)
length = 5

eigenvalues_base_model =  w_w_eigenvalues (base_model, num_modes) 
eigenvectors_base_model = w_w_eigenvectors ( base_model, eigenvalues_base_model, num_modes)

df_mode = pd.DataFrame()

for mode in range(0,num_modes):
	for i,d in enumerate(ds):
		difs = []
		for location in locations :
			damaged_model = Structure ("Healthy.xlsx",d,location,length)
			eigenvalues_damaged_model =  w_w_eigenvalues (damaged_model, num_modes) 
			relative_dif = (eigenvalues_damaged_model[mode] - eigenvalues_base_model[mode]) / eigenvalues_base_model[mode]
			difs.append(relative_dif)
		df_mode[i] = difs
print (df_mode	)
	# fig,axes = plt.subplots(nrows=3,ncols=3)
	# fig.set_size_inches(20, 15)
	# axes[0,0].plot(locations,df_mode[0])
	# axes[0,1].plot(locations,df_mode[1])
	# axes[0,2].plot(locations,df_mode[2])
	# axes[1,0].plot(locations,df_mode[3])
	# axes[1,1].plot(locations,df_mode[4])
	# axes[1,2].plot(locations,df_mode[5])
	# axes[2,0].plot(locations,df_mode[6])
	# axes[2,1].plot(locations,df_mode[7])
	# axes[2,2].plot(locations,df_mode[8])

	# axes[0,0].set_title(f'Effect of damage on natural frequency of \n mode {mode+1} when d is:{ds[0]}',fontsize=18)
	# axes[0,1].set_title(f'Effect of damage on natural frequency of \n mode {mode+1} when d is:{ds[1]}',fontsize=18)
	# axes[0,2].set_title(f'Effect of damage on natural frequency of \n mode {mode+1} when d is:{ds[2]}',fontsize=18)
	# axes[1,0].set_title(f'Effect of damage on natural frequency of \n mode {mode+1} when d is:{ds[3]}',fontsize=18)
	# axes[1,1].set_title(f'Effect of damage on natural frequency of \n mode {mode+1} when d is:{ds[4]}',fontsize=18)
	# axes[1,2].set_title(f'Effect of damage on natural frequency of \n mode {mode+1} when d is:{ds[5]}',fontsize=18)
	# axes[2,0].set_title(f'Effect of damage on natural frequency of \n mode {mode+1} when d is:{ds[6]}',fontsize=18)
	# axes[2,1].set_title(f'Effect of damage on natural frequency of \n mode {mode+1} when d is:{ds[7]}',fontsize=18)
	# axes[2,2].set_title(f'Effect of damage on natural frequency of \n mode {mode+1} when d is:{ds[8]}',fontsize=18)

	# axes[0,0].set_ylabel('Relative redection in natural frequency (%)',fontsize=15)	
	# axes[0,1].set_ylabel('Relative redection in natural frequency (%)',fontsize=15)	
	# axes[0,2].set_ylabel('Relative redection in natural frequency (%)',fontsize=15)	
	# axes[1,0].set_ylabel('Relative redection in natural frequency (%)',fontsize=15)	
	# axes[1,1].set_ylabel('Relative redection in natural frequency (%)',fontsize=15)	
	# axes[1,2].set_ylabel('Relative redection in natural frequency (%)',fontsize=15)	
	# axes[2,0].set_ylabel('Relative redection in natural frequency (%)',fontsize=15)	
	# axes[2,1].set_ylabel('Relative redection in natural frequency (%)',fontsize=15)	
	# axes[2,2].set_ylabel('Relative redection in natural frequency (%)',fontsize=15)	

	# axes[0,0].set_xlabel('Location (in)',fontsize=15)
	# axes[0,1].set_xlabel('Location (in)',fontsize=15)
	# axes[0,2].set_xlabel('Location (in)',fontsize=15)
	# axes[1,0].set_xlabel('Location (in)',fontsize=15)
	# axes[1,1].set_xlabel('Location (in)',fontsize=15)
	# axes[1,2].set_xlabel('Location (in)',fontsize=15)
	# axes[2,0].set_xlabel('Location (in)',fontsize=15)
	# axes[2,1].set_xlabel('Location (in)',fontsize=15)
	# axes[2,2].set_xlabel('Location (in)',fontsize=15)

	# plt.tight_layout()
	# plt.savefig(f'C://Users//15146//Desktop//FE thesis backup//FE(Thesis)//Thesis New (11 March)//Beam analysis in Python//Effect of damage on natural frequencies//{mode+1}.png')
	# plt.show()


#--------------------Checking cost function 
# num_modes = 3 
# base_model = Structure ("Healthy.xlsx",0,0,0)
# damaged_model = Structure ("Healthy.xlsx",0.3,15,5)

# eigenvalues_base_model =  w_w_eigenvalues (base_model, num_modes) 
# eigenvalues_damaged_model =  w_w_eigenvalues (damaged_model, num_modes) 

# eigenvectors_base_model = w_w_eigenvectors ( base_model, eigenvalues_base_model, num_modes)
# eigenvectors_damaged_model = w_w_eigenvectors ( damaged_model, eigenvalues_damaged_model, num_modes)

# print ( f'cost is {costfunction (eigenvalues_base_model,eigenvalues_damaged_model,eigenvectors_base_model,eigenvectors_damaged_model, damaged_model.added_dofs, base_model.added_dofs )}' )
# # #------------------printing modes and frequencies
# num_modes = 3
# base_model = Structure ("Healthy.xlsx",0,0,0)
# eigenvalues = w_w_eigenvalues (base_model, num_modes) 
# eigenvectors = w_w_eigenvectors ( base_model, eigenvalues, num_modes)
# print (eigenvectors)

#-------------------------------Printing the step load-----------------------------------


# num_modes = 1
# base_model = Structure ("Healthy.xlsx",0,0,0)

# l = Load (1, "Healthy.xlsx", base_model  )
# steady = steady_state (base_model,l)

# user_dof = 3
# u_steady = steady [user_dof-1]
# print (f'u_steady is {u_steady}')
# u_analytical = analytical( base_model, l , "step",num_modes)
# print (f'u_analytical is {u_analytical}')
# eigenvalues = w_w_eigenvalues (base_model, num_modes) 
# eigenvectors = w_w_eigenvectors ( base_model, eigenvalues, num_modes)

# print ("natural frequencies (Hz)",eigenvalues / (2*3.14) )
# print ("eigenvectors",eigenvectors)

# user_dof = 5
# transientt = transient (steady, l , base_model, eigenvalues, eigenvectors)
# print (f'transient is {transientt}')
# u_transient = transientt.iloc [user_dof-1].to_numpy() 

# u_total = u_steady + u_transient

# plt.plot(l.time,u_total,label='Spectral Element')
# plt.title('Response under step Load')
# plt.legend(loc='upper right')
# plt.ylabel('Displacement of middle (in)')
# plt.xlabel('t(s)')
# plt.show()



# #------------------------------------------------Printing sin load with different frequencies
# df_analytical = pd.DataFrame()
# df_total = pd.DataFrame()
# loads = pd.DataFrame()
# num_modes = 1
# base_model = Structure ("Healthy.xlsx",0,0,0)
# for f in range(1,3):
# 	if f == 1:
# 		l = Load (f, "Healthy.xlsx", base_model, step=1  )
# 	else:
# 		l = Load (1, "Healthy.xlsx", base_model, step=0  )
# 	steady = steady_state (base_model,l)
# 	# print("steady",steady)
# 	user_dof = 3
# 	u_steady = steady [user_dof-1]
# 	# print("i",i)
# 	# print("steady",steady)
# 	# print("usteady of ", user_dof)
# 	# print(u_steady)
	
# 	# plt.plot(l.time,u_steady, label = 'Spectral Steady')
# 	# plt.plot(l.time,u_analytical, label = 'analytical')
# 	# plt.legend()
# 	# plt.show()

# 	eigenvalues = w_w_eigenvalues (base_model, num_modes) 
# 	eigenvectors = w_w_eigenvectors ( base_model, eigenvalues, num_modes)

# 	print ("natural frequencies (Hz)",eigenvalues / (2*3.14) )
# 	print ("eigenvectors",eigenvectors)

# 	user_dof = 5
# 	transientt = transient (steady, l , base_model, eigenvalues, eigenvectors)
# 	print (f'transient is {transientt}')
# 	u_transient = transientt.iloc [user_dof-1].to_numpy() 
# 	# print("l.fi",l.fi)
# 	# print("transient",transient)
# 	# print("u_transient of",u_transient)
# 	# print("u_steady",u_steady)
# 	u_total = u_steady + u_transient
# 	df_total [f] = u_total
# 	loads [f] = l.amplitude
	 
# 	if f == 1:
# 		u_analytical = analytical( base_model, l , "step",num_modes)
# 		df_analytical [f] = u_analytical
# 	else:
# 		u_analytical = analytical( base_model, l , "sin",num_modes)
# 		df_analytical [f] = u_analytical
# 	# plt.plot(l.time,u_total, label = 'Spectral Analysis')
# 	# plt.plot(l.time,u_analytical, label = 'Analytical')
# 	# plt.xlabel ('t(s)')
# 	# plt.ylabel ('Displacement of middle (in)')
# 	# plt.legend()
# 	# plt.show()
# print (df_total)
# fig,axes = plt.subplots(nrows=2,ncols=2)
# f = 0
# font_size_label = 14
# font_size_title = 16
# print (l.time)
# axes [0,1].plot (l.time,df_total[1])
# axes [0,1].plot (l.time,df_analytical[1])
# axes [0,1].set_ylabel("Displacement of the middle (in)",fontsize = 12)
# axes [0,1].set_xlabel('t(s)',fontsize = font_size_label)
# axes [0,1].set_ylim([0.0,1.7])
# axes [0,1].set_xlim([0.0,1.0])
# axes [0,1].legend(['Spectral Element','Analytical'],loc='lower right',fontsize = 10)
# axes [0,1].set_xticklabels(axes [0,1].get_xticks(),fontsize=14)
# axes [0,1].set_yticklabels(axes [0,1].get_yticks(),fontsize=14)
# axes [0,1].yaxis.set_major_formatter(FormatStrFormatter('%0.01f'))
# axes [0,1].xaxis.set_major_formatter(FormatStrFormatter('%0.01f'))


# axes [0,0].plot (l.time,loads[1])
# axes [0,0].set_ylabel('P(t) (lb)',fontsize=font_size_label)
# axes [0,0].set_xlabel('t(s)',fontsize=font_size_label)
# axes [0,0].set_xlim([0.0,1.0])
# axes [0,0].set_xticklabels(axes [0,0].get_xticks(),fontsize=14)
# axes [0,0].set_yticklabels(axes [0,0].get_yticks(),fontsize=14)
# axes [0,0].yaxis.set_major_formatter(FormatStrFormatter('%0.01f'))
# axes [0,0].xaxis.set_major_formatter(FormatStrFormatter('%0.01f'))


# axes [1,1].plot (l.time,df_total[2])
# axes [1,1].plot (l.time,df_analytical[2])
# axes [1,1].set_ylabel("Displacement of the middle (in)",fontsize = 12)
# axes [1,1].set_xlabel('t(s)',fontsize = font_size_label)
# axes [1,1].legend(['Spectral Element','Analytical'],loc='upper right',fontsize = 10)
# axes [1,1].set_xlim ( [ 0.0 ,1.0 ] )
# axes [1,1].set_xticklabels(axes [1,1].get_xticks(),fontsize=14)
# axes [1,1].set_yticklabels(axes [1,1].get_yticks(),fontsize=14)
# axes [1,1].yaxis.set_major_formatter(FormatStrFormatter('%0.01f'))
# axes [1,1].xaxis.set_major_formatter(FormatStrFormatter('%0.01f'))


# axes [1,0].plot (l.time,loads[2])
# axes [1,0].set_ylabel('P(t) (lb)',fontsize = font_size_label)
# axes [1,0].set_xlabel('t(s)',fontsize = font_size_label)
# axes [1,0].set_xlim([0.0,1.0])
# axes [1,0].set_xticklabels(axes [1,0].get_xticks(),fontsize=14)
# axes [1,0].set_yticklabels(axes [1,0].get_yticks(),fontsize=14)
# axes [1,0].yaxis.set_major_formatter(FormatStrFormatter('%0.01f'))
# axes [1,0].xaxis.set_major_formatter(FormatStrFormatter('%0.01f'))

# axes [0,1].set_title ('Response under step load',fontsize = font_size_title)
# axes [0,0].set_title ('Step load',fontsize = font_size_title)
# axes [1,1].set_title ('Response under sine load',fontsize = font_size_title)
# axes [1,0].set_title ('Sine load',fontsize = font_size_title)


# fig.set_size_inches(8	, 6)
# plt.xticks(fontsize=14)
# plt.yticks(fontsize=14)
# plt.tight_layout()
# plt.savefig(f'C://Users//15146//Desktop//FE thesis backup//FE(Thesis)//Thesis New (11 March)//Beam analysis in Python//Doyle_Calculation.png')
# plt.show()

# fft_spectral = pd.DataFrame()
# fft_analytical = pd.DataFrame()

# for i in range(1,2):

# 	fft_spectral[i] = np.fft.rfft (df_total[i])
# 	fft_analytical[i] = np.fft.rfft (df_analytical[i])

# 	fft_spectral [i] = abs(fft_spectral[i])
# 	fft_analytical [i] = abs(fft_analytical[i])

# f = np.arange (start = 0, stop = l.samplingFrequency/2 + 1, step = l.samplingFrequency/l.n )

# fig,axes = plt.subplots(nrows=2,ncols=3)
# ff = 0
# for irow in range(0,2):
# 	for icol in range(0,3):
# 		ff = icol + irow * 3 + 1
# 		axes [irow,icol].plot (f,fft_spectral[ff])
# 		axes [irow,icol].plot (f,fft_analytical[ff])
# 		axes [irow,icol].set_title(f'Loading frequency:{ff} (Hz)')
# 		axes [irow,icol].legend(['Spectral Element','Analytical'],loc='upper right')
# 		axes [irow,icol].set_xlabel('f(Hz)')
# 		axes [irow,icol].set_ylabel('Abs ')


# plt.show()
# plt.plot(f,fft_spectral, label = 'Spectral Analysis + Veletsos')
# plt.plot(f,fft_analytical , label = 'analytical')
# plt.legend()
# plt.show()





# # #------------------------------------------------------
# create a damaged model as your real structure (D is constant)
# D = 0.25
# x_end = 35 
# l_damage = 5
# num_modes = 3
# #calculate eigenvalues and store
# my_structure1 = Structure ("Healthy.xlsx" , D , x_end, l_damage)
# eigenvalues_exp = w_w_eigenvalues ( my_structure1,num_modes )
# eigenvectors_exp = w_w_eigenvectors ( my_structure1, eigenvalues_exp, num_modes)
# #create different cases with different locations but constant 
# damage_index = np.arange ( 1,80,step=10 )/100
# damaged_element_end = np.arange ( 6 , 59 ,step= 1 )
# length_damage = np.array ( [5] )


# df = pd.DataFrame ( index = 'Location Width Cost'.split() )

# test = np.zeros  (  shape= (3)  )

# cnt = 0
# for j in range ( 0 , np.size(damage_index) ):
# 	for i in range ( 0 , np.size (damaged_element_end) ) :
# 		test [0] = damaged_element_end [ i ]
# 		test [1] = damage_index [ j ]
# 		# print ("damaged_element_end[i]=",damaged_element_end[i],"length_damage[j]=",damage_index [ j ])
# 		my_structure = Structure ("Healthy.xlsx" , damage_index [ j ] , damaged_element_end [i], length_damage[0] )
# 		eigenvalues_case = w_w_eigenvalues ( my_structure,num_modes )
# 		eigenvectors_case = w_w_eigenvectors ( my_structure, eigenvalues_case, num_modes)
# 		cost = costfunction ( eigenvalues_exp , eigenvalues_case, eigenvectors_exp, eigenvectors_case,my_structure.added_dofs, my_structure1.added_dofs)
# 		test [2] = cost
# 		df [cnt] = test
# 		cnt += 1
# 		print (cnt)

# print ( "df" , df)
# # df.to_excel("output.xlsx")

# fig = plt.figure()
# ax = plt.axes(projection="3d")


# ax.plot(df.iloc[0].to_numpy() ,df.iloc[1].to_numpy() ,df.iloc[2].to_numpy() ,color='r')
# ax.set_xlabel('Location (in)')
# ax.set_ylabel('Damage Index (%)')
# ax.set_zlabel('Cost')

# plt.show()
#-----------------------------------------------------
# create a damaged model as your real structure (D is a variable)
# D = 0.3
# x_end = 50 
# l_damage = 8

# #calculate eigenvalues and store
# my_structure = Structure ("Healthy.xlsx" , D , x_end, l_damage)
# eigenvalues_exp = w_w_eigenvalues ( my_structure,5 )

# #create different cases with different locations but constant 
# damage_index = np.arange ( 1,100,step=10 )/100
# damaged_element_end = np.arange ( 10 , 60 ,step= 1 )
# length_damage = np.array ( [8] )

# df = pd.DataFrame ( index = 'Location Width Cost'.split() )
# # D = 0.3
# test = np.zeros  (  shape= (3)  )

# cnt = 0
# for j in range ( 0 , np.size(damage_index) ):
# 	for i in range ( 0 , np.size (damaged_element_end) ) :
# 		print ( i )
# 		test [0] = damaged_element_end [ i ]
# 		test [1] = damage_index [ j ]
# 		print ("damaged_element_end[i]=",damaged_element_end[i],"damage_index[j]=",damage_index[j])
# 		my_structure = Structure ("Healthy.xlsx" , damage_index[j] , damaged_element_end[i], length_damage[0])
# 		eigenvalues_case = w_w_eigenvalues ( my_structure,5 )
# 		cost = costfunction ( eigenvalues_exp , eigenvalues_case  )
# 		test [2] = cost
# 		df [cnt] = test
# 		cnt += 1
# 		print (cnt)

# print ( "df" , df)
# df.to_excel("output.xlsx")

# fig = plt.figure()
# ax = plt.axes(projection="3d")


# ax.plot(df.iloc[0].to_numpy() ,df.iloc[1].to_numpy() ,df.iloc[2].to_numpy() ,color='r')
# ax.set_xlabel('Location')
# ax.set_ylabel('damage index')
# ax.set_zlabel('Cost')

# plt.show()

















#----------------------------------------------------------------------------
# damage_index = np.array ( [  0.3 ] )
# damaged_element_end = np.array ( [ 4 , 7, 10, 13, 16, 19, 22, 25 ] )
# length_damage = np.array ( [ 3 ] )
# df = pd.DataFrame ()

# cnt = 0
# for i in range ( 0, np.size (length_damage) ) :
# 	# base model 
# 	D = 0
# 	x_end = 0
# 	l_damage = 0
# 	my_structure = Structure ("Healthy.xlsx" , D , x_end, l_damage)
# 	eigenvalues_basemodel = w_w_eigenvalues ( my_structure )
# 	props = np.array ( [ D , x_end, l_damage] )
# 	data = np.concatenate( (props, eigenvalues_basemodel), axis=0)
# 	df [cnt] = data
# 	cnt += 1
# 	#Creating DataBase for Damage Detection
# 	l_damage = length_damage [i]
# 	for j in range ( 0, np.size (damaged_element_end) ) :
# 		x_end = damaged_element_end [j]
# 		for q in range  ( 0, np.size ( damage_index ) ):
# 			D = damage_index [q]
# 			my_structure = Structure ("Healthy.xlsx" , D , x_end , l_damage )
# 			eigenvalues = w_w_eigenvalues ( my_structure )  
# 			# eigenvectors = w_w_eigenvectors ( my_structure, eigenvalues)
# 			cost = 
# 			props = np.array ( [ D , x_end, l_damage] )
# 			data = np.concatenate( (props, eigenvalues), axis=0)
# 			df [cnt] = data
# 			cnt += 1
# 			print(cnt)

# df.to_csv('Classification DataSet CSV.csv',index = None, header=True)

# print ("eigenvalues",eigenvalues )
# print ("eigenvectors",eigenvectors )


#-------------------------------------------to plot response of damaged structure

# damage_index = np.array ( [0.1 , 0.4, 0.7] )
# damaged_element_end = np.array ( [ 6 , 10, 14, 18, 22, 26 ] )
# length_damage = np.array ( [ 1 , 2, 3 , 4, 5 ] )
# df = pd.DataFrame ()

# cnt = 0

# for i in range ( 0, np.size (length_damage) ) :
# 	#base model 
# 	D = 0
# 	x_end = 0
# 	l_damage = 0
# 	my_structure = Structure ("Healthy.xlsx" , D , x_end, l_damage)
# 	eigenvalues_basemodel = w_w_eigenvalues ( my_structure )
# 	props = np.array ( [ D , x_end, l_damage] )
# 	data = np.concatenate( (props, eigenvalues_basemodel), axis=0)
# 	df [cnt] = data
# 	cnt += 1
# 	#Creating DataBase for Damage Detection
# 	l_damage = length_damage [i]
# 	for j in range ( 0, np.size (damaged_element_end) ) :
# 		x_end = damaged_element_end [j]
# 		for q in range  ( 0, np.size ( damage_index ) ):
# 			D = damage_index [q]
# 			my_structure = Structure ("Healthy.xlsx" , D , x_end , l_damage )
# 			l = Load (3, "BaseModel.xlsx" )
# 			steady = steady_state (my_structure,l)
# 			user_dof = 3
# 			u_steady = steady [user_dof-1]
# 			eigenvalues = w_w_eigenvalues ( my_structure )  - eigenvalues_basemodel
# 			eigenvectors = w_w_eigenvectors ( my_structure, eigenvalues)
# 			u_analytical = analytical( my_structure, l , "sin")
# 			user_dof = 5
# 			transientt = transient (steady, l , my_structure, eigenvalues, eigenvectors)
# 			u_transient = transientt.iloc [user_dof-1].to_numpy()
# 			u_total = u_steady + u_transient
# 			plt.plot(l.time,u_total, label = 'total')
# 			plt.plot(l.time,u_analytical, label = 'analytical')
# 			plt.legend()
# 			plt.show() 
# 			props = np.array ( [ D , x_end, l_damage] )
# 			data = np.concatenate( (props, eigenvalues), axis=0)
# 			df [cnt] = data
# 			cnt += 1
# 			print(cnt)

#-------------------------------------------------------------------------------

# ratio = np.zeros(shape = (5) )
# freqs = np.arange(start = 0, stop =5, step =1 )
# for i in freqs:
# 	print("i",i)
# 	l = Load (i, "Healthy.xlsx" )
# 	base_model = Structure ("Healthy.xlsx",0,0,0)
# 	steady = steady_state (base_model,l)
# 	# print("steady",steady)
# 	user_dof = 3
# 	u_steady = steady [user_dof-1]
# 	# print("i",i)
# 	# print("steady",steady)
# 	# print("usteady of ", user_dof)
# 	# print(u_steady)
# 	u_analytical = analytical( base_model, l , "step")

# 	# plt.plot(l.time,u_steady, label = 'Spectral Steady')
# 	# plt.plot(l.time,u_analytical, label = 'analytical')
# 	# plt.legend()
# 	# plt.show()

# 	eigenvalues = w_w_eigenvalues (base_model) 
# 	eigenvectors = w_w_eigenvectors ( base_model, eigenvalues)

# 	print ("natural frequencies (Hz)",eigenvalues / (2*3.14) )
# 	# print ("eigenvectors",eigenvectors)

# 	user_dof = 5
# 	transientt = transient (steady, l , base_model, eigenvalues, eigenvectors)
# 	u_transient = transientt.iloc [user_dof-1].to_numpy() 
# 	# print("l.fi",l.fi)
# 	# print("transient",transient)
# 	# print("u_transient of",u_transient)
# 	# print("u_steady",u_steady)
# 	u_total = u_steady + u_transient

# 	plt.plot(l.time,u_total, label = 'total')
# 	plt.plot(l.time,u_analytical, label = 'analytical')
# 	plt.legend()
# 	plt.show()


# 	fft_spectral = np.fft.rfft (u_total)
# 	fft_analytical = np.fft.rfft (u_analytical)

# 	fft_spectral = abs(fft_spectral)
# 	fft_analytical = abs(fft_analytical )

# 	f = np.arange (start = 0, stop = l.samplingFrequency/2 + 1, step = l.samplingFrequency/l.n )
# 	# print("i",i)
# 	ratio [i] = np.amax (u_total)  /  np.amax (u_analytical)
# 	# print (ratio[i])
# 	#plt.plot(f,fft_spectral, label = 'total')
# 	#plt.plot(f,fft_analytical , label = 'analytical')
# 	#plt.legend()
# 	#plt.show()

# # print("max_total",max_total)
# # print ("max_analytical ",max_analytical )


# plt.plot(freqs , ratio   , label = 'ratio')
# plt.legend()
# plt.show()

#-------------------Printing natural frequencies
# my_structure = Structure ('Healthy.xlsx',0,0,0)
# eigenvalues = w_w_eigenvalues (my_structure,4)
# print (eigenvalues) 
# eigenvectors = w_w_eigenvectors ( my_structure, eigenvalues, 4)
# print (eigenvectors)
# eigenvectors.to_excel("modes.xlsx")


















































































