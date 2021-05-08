from Input import Load
from structure import Structure 
import matplotlib.pyplot as plt
import math
import numpy as np


def analytical ( my_structure , l , load_type, nfreq ) :
	if load_type == "sin" :
		pi = 3.141592653594
		omega_bar = 2 *pi * l.signalFrequency
		# print ("analyitical load frequency=", l.signalFrequency)
		p0 = 100 # Amplitude of load
		rho = my_structure.density[0]
		a = my_structure.a [0]
		m = rho * a
		e = my_structure.e [0]
		i = my_structure.ii [0]
		ll = 60
		delta_t = l.samplingInterval
		# print("analytical")
		# print("a=",a,"m",m,"e",e,"i",i,"ll",ll,"delta_t",delta_t,"omega_bar",omega_bar)
		t = np.arange ( start =0.0, stop = (l.n) * l.samplingInterval, step=l.samplingInterval)
		# print("time",t)
		resp = np.empty(l.n)
		q = 0
		for j in t:
			un = 0
			for n in np.arange (start=1, stop =nfreq + 1, step = 1):
				omega =( ( (n**2) * (pi**2) ) / (ll**2) ) * ((e*i/m) ** 0.5)
				# print (f'omega used is ;{omega}')
				#print ("omega",omega)
				kn = ( (n**2)*(pi**4)*e*i ) / (2*ll**3)
				v = math.sin (n*pi/2)
				w = p0 * v / kn
				x = ( 1 / ( 1 - (omega_bar/omega)**2 ) )
				y = (omega_bar/omega) * x
				un = un + v * ( w * x * math.sin ( omega_bar * j )  ) - w * y * math.sin(omega * j)
			resp[q] = un
			q += 1
		return (resp)	

	if load_type == "step" :

		pi = 3.14159265359
		p0 = 1000 # Amplitude of load
		rho = my_structure.density[0]
		a = my_structure.a [0]
		m = rho * a
		e = my_structure.e [0]
		i = my_structure.ii [0]
		ll = 60
		t = np.arange ( start =0.0, stop = (l.n) * l.samplingInterval, step=l.samplingInterval)
		resp = np.empty(l.n)
		q = 0
		for j in t:
			un = 0
			for n in np.arange (start=1, stop =nfreq + 1, step = 1):
				omega =( ( (n**2) * (pi**2) ) / (ll**2) ) * ((e*i/m) ** 0.5)
				kn = ( (n**2)*(pi**4)*e*i ) / (2*ll**3)
				v = math.sin (n*pi/2)
				w = p0 * v / kn
				un = un + w * ( 1  - math.cos (omega * j ) )
			resp[q] = un
			q += 1
		return (resp)	

