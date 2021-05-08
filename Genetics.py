import numpy as np
import random
from structure import Structure 
from W_W_EigenValues import w_w_eigenvalues ,w_w_eigenvectors
from DamageCostFunction import costfunction
import matplotlib.pyplot as plt
import pandas as pd
def compare (genome1,genome2):
	n = np.shape (genome1)[0]
	compared = 0
	for i in range (0,n):
		if genome1[i] != genome2[i] :
			compared += 1
	return compared

def diversity (population):
	pop = population
	n = np.shape(pop)[0]
	diversity = np.zeros ( shape =(n,n) ) 
	for i in range (0,n):
		for j in range (n-1,i-1,-1):
			genome1 = pop [i,:]
			genome2 = pop [j,:]
			compared = compare (genome1,genome2)
			diversity[i,j] = compared 
			diversity[j,i] = compared
	d = np.sum(diversity,axis=1)
	arr1inds = d.argsort()
	d = d [arr1inds[::-1]]
	pop = pop [arr1inds[::-1]]
	return  pop[0,:] , pop[1,:] 

def genome_reader (genome ,cut, cut2 ):
	location = genome [0:cut]
	width = genome [cut:cut2]
	d_index = genome [cut2:]
	l=0
	for i in range (0,np.size(location)):
		index = np.size(location) - 1 - i
		l += 2**i * location [index]

	w=0
	for i in range (0,np.size(width)):
		index = np.size(width) - 1 - i
		w += 2** i * width [index]
	d=0
	for i in range (0,np.size(d_index)):
		index = np.size(d_index) - 1 - i
		d += 2** i * d_index [index]
	return l,w,d

def dec_to_binary ( n ) :
	arr_int = [int (char) for char in bin(n)[2:] ]
	arr = np.array (arr_int)
	return arr 

def generate_genome ( max_location, max_width, max_d) :
	#create an array containig 2 values: Location of damage in binary + width of damage in binary
	max_l_bin = np.size ( dec_to_binary ( max_location ) )
	max_w_bin = np.size ( dec_to_binary ( max_width ) )
	max_d_bin = np.size ( dec_to_binary ( max_d ) )
	genome_size = max_l_bin + max_w_bin + max_d_bin
	a = np.zeros (max_l_bin)
	c = np.zeros (max_w_bin)
	f = np.zeros (max_d_bin)

	location_dec = random.randint ( max_width + 1 , max_location )
	width_dec = random.randint ( 1 , max_width ) 

	# print (f'generate genome gives:location is {location_dec} and width is {width_dec}')
	while location_dec - width_dec <= 0: 
		print ('generate genome giving bad solutions:')
		print (f'location is {location_dec} and width is {width}')
		location_dec = random.randint ( max_width + 1 , max_location )
		width_dec = random.randint ( 1 , max_width ) 
	d_dec = random.randint(1,100)
	location = dec_to_binary ( location_dec ) 
	width = dec_to_binary ( width_dec )
	d_index = dec_to_binary (d_dec)
	a [ max_l_bin - np.size (location) : ] = location
	c [ max_w_bin - np.size ( width) : ] = width 
	f [ max_d_bin - np.size ( d_index) : ] = d_index 
	loc_width = np.concatenate ( (a, c), axis=None )
	genome = np.concatenate ( (loc_width, f), axis=None )
	cut = max_l_bin
	cut2 = max_l_bin + max_w_bin
	return genome, cut, genome_size, max_location, max_width, max_d, cut2

def generate_population ( pop_size , eigenvalue_exp, max_location, max_width, max_d ) :

	genome , cut , genome_size, max_location, max_width, max_d, cut2 = generate_genome ( max_location, max_width, max_d )
	population = np.empty ( (pop_size, genome_size ) )
	for i in range (0,pop_size) :
		genome , cut , genome_size, max_location, max_width, max_d,cut2 = generate_genome ( max_location, max_width, max_d)
		population [i,0:] = genome
	return population, cut


def fitness (genome, cut, max_location, max_width, eigenvalues_exp, cut2, max_d, eigenvectors_base_model,exp_added_dofs) :
	#check if genome is in proper bounds
	location,width,d_index = genome_reader(genome,cut,cut2)
	d_index /= 100
	# print (f'damage location is {location} and width is {width} and d is {d_index}')
	if (location - width <= 0 ):
		quit ()
	if ( location > max_location ) or ( width > max_width ) or ( d_index > max_d):
		cost = 10**7 # a very big cost so it would not evolve
	else :
		my_structure = Structure ("Healthy.xlsx" ,  d_index , location, width)
		eigenvalues_case = w_w_eigenvalues ( my_structure,3 )
		#added for New cost with modes
		eigenvectors_case = w_w_eigenvectors ( my_structure, eigenvalues_case, 3)
		#cost = costfunction ( eigenvalues_exp , eigenvalues_case  )
		cost = costfunction (eigenvalues_exp,eigenvalues_case,eigenvectors_base_model,eigenvectors_case, my_structure.added_dofs,exp_added_dofs )
	return cost



def selection_pair (population, weights) :
	weights = ( weights  / np.amax ( weights ) )
	weights = 1 - weights
	return random.choices ( 
		population = population,
		weights = weights,
		k=2
		)



def cross_over ( parents , cut ,cut2, max_location, max_width,max_d,type  ) :
	if (type == 1 ):
		genome_a = parents [0,:]
		genome_b = parents [1,:]

		if  len ( genome_a ) != len ( genome_b ) : 
			print ("Genomes a and b must be of the same siez")
		length = np.size(genome_a)
		if length > 2 :
			#Do cross over
			p = random.randint ( 1 , cut )

			genome_a_copy = np.copy (genome_a)
			genome_a = np.append ( genome_a [0:p] , genome_b [p:cut])
			genome_a = np.append ( genome_a [0:cut], genome_a_copy[cut:] )

			genome_b_copy = np.copy (genome_b)
			genome_b = np.append ( genome_b [0:p] , genome_a_copy [p:cut])
			genome_b = np.append ( genome_b [0:cut], genome_b_copy[cut:] )

			#Crossover on width
			p = random.randint (cut, cut2)

			genome_a = np.append ( genome_a [0:p] , genome_b [p:] )

			genome_b = np.append ( genome_b [0:p] , genome_a_copy [p:] )

			#Crossover on damage index
			p = random.randint (cut2, length)

			genome_a = np.append ( genome_a [0:p] , genome_b [p:] )

			genome_b = np.append ( genome_b [0:p] , genome_a_copy [p:] )

			#Check if the solutions are inside the constraints
			l_a,w_a,d_a =genome_reader(genome_a, cut, cut2)
			l_b,w_b,d_b =genome_reader(genome_b, cut, cut2)

			while (l_a - w_a <= 0 or l_b - w_b <= 0 or l_a > max_location or l_b > max_location or w_a >max_width or w_b > max_width ):

				#Do cross over
				p = random.randint ( 1 , cut )
				# print (f'mutation produced negative : l is {l}, w is {w}, d is {d}')
				genome_a_copy = np.copy (genome_a)
				genome_a = np.append ( genome_a [0:p] , genome_b [p:cut])
				genome_a = np.append ( genome_a [0:cut], genome_a_copy[cut:] )

				genome_b_copy = np.copy (genome_b)
				genome_b = np.append ( genome_b [0:p] , genome_a_copy [p:cut])
				genome_b = np.append ( genome_b [0:cut], genome_b_copy[cut:] )

				#Crossover on width
				p = random.randint (cut, cut2)

				genome_a = np.append ( genome_a [0:p] , genome_b [p:] )

				genome_b = np.append ( genome_b [0:p] , genome_a_copy [p:] )

				#Crossover on damage index
				p = random.randint (cut2, length)

				genome_a = np.append ( genome_a [0:p] , genome_b [p:] )

				genome_b = np.append ( genome_b [0:p] , genome_a_copy [p:] )

				#Check if the solutions are inside the constraints
				l_a,w_a,d_a =genome_reader(genome_a, cut, cut2)
				l_b,w_b,d_b =genome_reader(genome_b, cut, cut2)
	if (type==2):
		genome_a = parents [0,:]
		genome_b = parents [1,:]

		if  len ( genome_a ) != len ( genome_b ) : 
			print ("Genomes a and b must be of the same siez")
		length = np.size(genome_a)
		if length > 2 :
			a = random.randint(0,2)
			if a==0:
				p=cut
				print (f'swap locations')
				genome_a_copy = np.copy (genome_a)
				genome_a = np.append ( genome_b [0:p] , genome_a [p:])
				genome_b = np.append ( genome_a_copy [0:p] , genome_b [p:])
			if a==1:
				print (f'swap widths')
				genome_a_copy = np.copy (genome_a)
				genome_a = np.append( genome_a[0:cut] , genome_b[cut:cut2] )
				genome_a = np.append( genome_a , genome_a_copy[cut2:] )

				genome_b_copy = np.copy (genome_b)
				genome_b = np.append( genome_b[0:cut] , genome_a_copy [cut:cut2] )
				genome_b = np.append( genome_b , genome_b_copy [cut2:] )
			if a==2:
				print (f'swap damage indices')
				genome_a_copy = np.copy (genome_a)
				genome_a = np.append ( genome_a [0:cut2] , genome_b [cut2:])
				genome_b = np.append ( genome_b [0:cut2] , genome_a_copy [cut2:])

			#Check if the solutions are inside the constraints
			l_a,w_a,d_a =genome_reader(genome_a, cut, cut2)
			l_b,w_b,d_b =genome_reader(genome_b, cut, cut2)
			alert = 0
			while (l_a - w_a <= 0 or l_b - w_b <= 0 or l_a > max_location or l_b > max_location or w_a >max_width or w_b > max_width ):
					
				print (f'Im stucked in loop for the {alert}th time')
				print (f'l_a is{l_a}, l_b is {l_b}')
				print (f'w_a is{w_a}, w_b is {w_b}')
				print (f'd_a is{d_a}, d_b is {d_b}')
				alert+=1
				a = random.randint(0,2)
				if a==0:
					p=cut
					print (f'swap locations')
					genome_a_copy = np.copy (genome_a)
					genome_a = np.append ( genome_b [0:p] , genome_a [p:])
					genome_b = np.append ( genome_a_copy [0:p] , genome_b [p:])
				if a==1:
					print (f'swap widths')
					genome_a_copy = np.copy (genome_a)
					genome_a = np.append( genome_a[0:cut] , genome_b[cut:cut2] )
					genome_a = np.append( genome_a , genome_a_copy[cut2:] )

					genome_b_copy = np.copy (genome_b)
					genome_b = np.append( genome_b[0:cut] , genome_a_copy [cut:cut2] )
					genome_b = np.append( genome_b , genome_b_copy [cut2:] )
				if a==2:
					print (f'swap damage indices')
					genome_a_copy = np.copy (genome_a)
					genome_a = np.append ( genome_a [0:cut2] , genome_b [cut2:])
					genome_b = np.append ( genome_b [0:cut2] , genome_a_copy [cut2:])

				l_a,w_a,d_a =genome_reader(genome_a, cut, cut2)
				l_b,w_b,d_b =genome_reader(genome_b, cut, cut2)

	return genome_a, genome_b



def mutation ( genome , cut, max_location, max_width, cut2 ) :
	probability = 0.5
	#increasing the probabilty of changing values on the left more
	indices = np.arange (0, len (genome) )
	weights = indices [ : : -1 ]
	index = random.choices (
		population = indices,
		weights = weights,
		k=1
		)
	index = random.randrange ( len (genome)  )
	genome [ index ] = genome [index] if ( random.randint(1, 100) / 100 ) < probability else abs ( genome [ index ] - 1 )
	l,w,d =genome_reader(genome , cut, cut2)
	while ( l > max_location or w > max_width or l - w <= 0 or d > 100 or w <= 0 or d > 98 ):
		# print (f'mutation produced negative : l is {l}, w is {w}, d is {d}')
		index = random.randrange ( len (genome)  )
		genome [ index ] = genome [index] if ( random.randint(1, 100) / 100 ) < probability else abs ( genome [ index ] - 1 )
		l,w,d =genome_reader(genome , cut, cut2)
	# print (f'mutation correct change is :l is {l}, w is {w}, d is {d}')
	# trying multiple mutation

	index = random.randrange ( len (genome)  )
	# print (index)
	genome [ index ] = genome [index] if ( random.randint(1, 100) / 100 ) < probability else abs ( genome [ index ] - 1 )
	l,w,d =genome_reader(genome , cut, cut2)
	while ( l > max_location or w > max_width or l - w <= 0 or d > 100 or w <= 0 or d > 98 ):
		# print (f'mutation produced negative : l is {l}, w is {w}, d is {d}')
		index = random.randrange ( len (genome)  )
		genome [ index ] = genome [index] if ( random.randint(1, 100) / 100 ) < probability else abs ( genome [ index ] - 1 )
		l,w,d =genome_reader(genome , cut, cut2)
	index = random.randrange ( len (genome) )
	# print (index)
	genome [ index ] = genome [index] if ( random.randint(1, 100) / 100 ) < probability else abs ( genome [ index ] - 1 )
	l,w,d =genome_reader(genome , cut, cut2)
	while ( l > max_location or w > max_width or l - w <= 0 or d > 100 or w <= 0  or d > 98 ):
		# print (f'mutation produced negative : l is {l}, w is {w}, d is {d}')
		index = random.randrange ( len (genome)  )
		genome [ index ] = genome [index] if ( random.randint(1, 100) / 100 ) < probability else abs ( genome [ index ] - 1 )
		l,w,d =genome_reader(genome , cut, cut2)
	return genome

# #Do an experiment
print('This is the test') 
location = 35
width = 5
d = 0.25
my_structure = Structure ("Healthy.xlsx" , d , location, width)
eigenvalue_exp = w_w_eigenvalues ( my_structure, 3 )
eigenvectors_base_model = w_w_eigenvectors ( my_structure,eigenvalue_exp, 3)
# eigenvectors_base_model.to_excel("output.xlsx")
exp_added_dofs = my_structure.added_dofs
#Setting Problem constraints
generation_limit = 500
pop_size = 8
max_width = 20
max_d = 99
max_location = 59
fitness_limit = 0.0000000000001																
results = np.empty ( shape=(1,2) )
#get size of the genome and cut
genome , cut , genome_size, max_location, max_width,max_d,cut2 = generate_genome (max_location,max_width,max_d )
print (f'cut is {cut} and cut2 is {cut2}')
next_generation = np.zeros ( shape = ( pop_size , genome_size ) )

for i_generation in range ( 0, generation_limit):
	#create population
	print ("generation number= ", i_generation)
	if ( i_generation == 0 ):
		#create the first population with random values
		population , cut = generate_population (pop_size , eigenvalue_exp ,max_location,max_width,max_d)
		print (f'first population out of generate_population is {population}')
		cost = np.zeros ( shape = (pop_size) )
		#calculating cost
		for i in range (0, pop_size) :
			cost [i] = fitness (population [i,:], cut, max_location, max_width, eigenvalue_exp, cut2,max_d,eigenvectors_base_model,exp_added_dofs)
		#Sorting the populatin in ascending order based on cost of each genome
		cost = cost * -1
		arr1inds = cost.argsort()
		cost = cost [arr1inds[::-1]]
		population = population [arr1inds[::-1]]
		cost = cost * -1
	else :
		population = next_generation
		cost = np.zeros ( shape = (pop_size) )
		#calculating cost
		for i in range (0, pop_size) :
			cost [i] = fitness (population [i,:], cut, max_location, max_width, eigenvalue_exp, cut2, max_d,eigenvectors_base_model,exp_added_dofs)
		#Sorting the populatin in ascending order based on cost of each genome
		cost = cost * -1
		arr1inds = cost.argsort()
		cost = cost [arr1inds[::-1]]
		population = population [arr1inds[::-1]]
		cost = cost * -1

	#checking if the population has already reached a good value
	if  cost [0] < fitness_limit : 
		break
	# print ("best solution", population[])
	l,w,d = genome_reader (population[0,:], cut, cut2)
	print (f'Best Solution: location is {l}  and width is {w} and d is {d}%')
	# #addding just to print decimal values of population
	# for i in range (pop_size):
	# 	l,w,d = genome_reader (population[i,:], cut, cut2)
	# 	print (f'for genome number {i} location is {l}  and width is {w} and d is {d}%')
	# print ("cost", cost)

#creating the next generation
	#keeping the first best two solutions from previous generation
	next_generation [ 0, : ] = population [0,:]
	next_generation [ 1, : ] = population [1,:]
			### take the two most diverse solutions and send to the next generatio
			### So this time, do crossover and mutation on the last 4 genomes instead of the last six
	# most_diverse,second_most_diverse = diversity (population)
	# next_generation [ 2, : ] = most_diverse
	# next_generation [ 3, : ] = second_most_diverse
	# cnt=4
	cnt = 2
	for i_parents in range (0 , int ( pop_size / 2 ) - 1) :#for diversitty -2
		#picking two parents two reproduce
		parents = np.array ( selection_pair ( population, cost ) )
	#creating offsprings
		#crossover
		if i_generation <100:
			type=1
		else:
			type=2
		offspring_a , offspring_b = cross_over ( parents ,cut, cut2, max_location, max_width,max_d,type=1 )
		#mutation
		offspring_a = mutation ( offspring_a ,cut ,max_location,max_width,cut2)
		offspring_b = mutation ( offspring_b ,cut ,max_location,max_width,cut2)
	#putting offsprings into the new generation
		next_generation [ cnt, : ] = offspring_a
		cnt += 1
		next_generation [ cnt, : ] = offspring_b
		cnt += 1
	result = np.array ([ i_generation, cost[0] ])
	result = np.reshape(result,newshape=(1,2))
	# print (f'results of generation{i_generation} is {result}')
	plt.xlabel('Generation Number')
	plt.ylabel('Cost function')
	plt.scatter(i_generation, cost[0] )
	plt.pause(0.01)
	results = np.append(results,result,axis = 0)

print (f'results is {results}')
print ("Winner population", population)
l,w,d = genome_reader (population[0,:],cut, cut2)
print (f"Location of damge is {l} and width is {w} and d is {d}%")
results = np.delete(results,0,0)
plt.scatter(results[:,0],results[:,1])
plt.xlabel('Generation Number')
plt.ylabel('Cost function')
# plt.plot(results[:,0],results[:,1])
plt.show()



	