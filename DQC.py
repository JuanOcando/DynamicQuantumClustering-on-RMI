#import libraries
import numpy as np 
#import scipy as sp
#import sympy as sy
import nibabel as nib
import matplotlib.pyplot as plt
import sys
from backports import *

def show_slices(slices):
	fig, axes = plt.subplots(1,len(slices))
	
	for i, slice in enumerate(slices):
		axes[i].imshow(slice.T, cmap = "gray", origin = "lower")

	plt.show()

def get_voxel_slices(data,x,y,z):
	slice_x = data[-x,:,:]
	slice_y = data[:,-y,:]
	slice_z = data[:,:,z ]
	print(slice_x.shape)
	print(slice_y.shape)
	print(slice_z.shape)
	return slice_x, slice_y, slice_z

def one_dim_wave_funct(data,sigma):

	#To do: Reduce decimals

	time =  np.linspace(0,600, num=200)
	W = np.zeros(200)
	spec = np.zeros(200)
	V = np.zeros(200)

#	for k in range(time.shape[0]):
#		print(time[k])
	

	for k in range(time.shape[0]):

		for i in range(data.shape[0]):
			
	 		for j in range(data.shape[1]):
	 		
			 	#if(data[i,j]>0):
			 	W[k] = W[k] +  np.exp(-(time[k]-data[i,j])*(time[k]-data[i,j])/(2*sigma*sigma))
			 	spec[k] = spec[k] + (time[k]-data[i,j])*(time[k]-data[i,j])*np.exp(-(time[k]-data[i,j])*(time[k]-data[i,j])/(2*sigma*sigma))


	for k in range(time.shape[0]):
		#if(data[i,j]>0):
		V[k] = (spec[k]/(W[k]*sigma*sigma) - 1)/(2*sigma**4)
		#V[k] = (spec[k]/(W[k]*sigma*sigma) - 1.0)*0.5
		#W[k] = W[k]/(sigma*2.5066)
		print(V[k])

	energy = np.amin(V)

	for k in range(time.shape[0]):
		V[k] = V[k] - energy


	return time, W, V

def minima(V):

	mini = np.zeros(1)

	for k in range(1,V.shape[0]-1):

		if V[k-1] > V[k] and V[k] < V[k+1] :
			np.append(mini, k)

	return mini







epi_img = nib.load('Brats17_CBICA_AQJ_1_t1.nii.gz')
#epi_img = nib.load('Brats17_CBICA_AQJ_1_t2.nii.gz')
#epi_img = nib.load('Brats17_CBICA_AQJ_1_t1ce.nii.gz')

epi_img_data = epi_img.get_data()

slice_0, slice_1, slice_2 = get_voxel_slices(epi_img_data,123,121,90)

t, W, V = one_dim_wave_funct(slice_0,1)

plt.subplot(211)
plt.plot(t,W)
plt.suptitle("w") 

plt.suptitle("Probabilidad vs Energia Potencial") 

plt.subplot(212)
plt.plot(t,V)
plt.suptitle("v") 

show_slices([slice_0, slice_1, slice_2])
plt.suptitle("Center slices for EPI image") 

print(minima(V))

