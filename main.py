import numpy as np 
import matplotlib.pyplot as plt 
from math import *

class Wave(object):
	def __init__(self,s0,w,k):
		self.s0 = s0
		self.w = w
		self.k = k

def get_wave_function(wave,x,y,t):
    r = np.asarray(np.sqrt(np.power(x,2) + np.power(y,2)),dtype=np.float)
    print len(r),len(r[1])
    return (wave.s0/r)*np.exp(wave.k*r*1j - wave.w*t*1j)
 
#%%
onde = Wave(1,2*pi*440,sqrt(2*pi*440/340.))
#%%
X = np.linspace(-2,2,1000)
Y = np.linspace(-2,2,1000)
XX,YY = np.meshgrid(X,Y)
#%%
result = get_wave_function(onde,XX,YY,1)
plt.imshow(result.astype("float"))