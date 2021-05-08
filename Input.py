import math 
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from structure import Structure

class Load (Structure):

	def __init__ (self,signal_freq, filename, my_structure, step):

		#input from user
		num_points = 256
		self.samplingInterval = 1/(num_points)
		self.n = num_points # Number of Data points
		self.samplingFrequency = 1 / self.samplingInterval
		self.nyquist =int (self.samplingFrequency // 2 ) 
		self.signalFrequency = signal_freq
		self.time = np.arange(0, (self.n) * self.samplingInterval , self.samplingInterval)
		if step == 0 :
			self.amplitude =    np.sin(2*np.pi*self.signalFrequency*self.time)
			self.amplitude *= 100 #Amlitude of sin load

		#step load
		if step == 1:
			self.amplitude = np.array([1000]*(num_points) )

		self.fft_load = np.fft.rfft(self.amplitude)
		self.load_real = self.fft_load.real
		self.load_imag = self.fft_load.imag

		df = pd.read_excel(filename, sheet_name ='NENJNS')
		nj = my_structure.nj
		ns = my_structure.ns
		nf = nj * 3 - ns
		self.fr = np.zeros (nj*3) #the size is 9. It should not be!!!
		self.fi = np.zeros (nj*3)
		lj = my_structure.lj
		alldofs = np.arange(start=1, stop=nf + 1, step=1)
		self.load_vecs_real = pd.DataFrame ( )
		self.load_vecs_imag = pd.DataFrame ( )
		for i in range (1, int ( self.n//2 + 1 ) + 1 ) : 
			self.fr [ lj - 1 ] = self.load_real [i - 1]
			self.fi [ lj - 1 ] = self.load_imag [i - 1]
			self.load_vecs_real [i] = self.fr
			self.load_vecs_imag [i] = self.fi
		df = pd.read_excel(filename, sheet_name ='Fixed')
		fix = my_structure.fix
		self.load_vecs_real = self.load_vecs_real.drop ( fix , axis=0)
		self.load_vecs_imag = self.load_vecs_imag.drop ( fix , axis=0)
		self.fr = np.delete ( self.fr, fix , 0)
		self.fi = np.delete ( self.fi, fix , 0)
		self.freqs = np.arange(start=1, stop=self.n + 1, step=1) * self.samplingFrequency / self.n

