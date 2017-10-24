#coding:utf-8
import numpy as np 
import matplotlib.pyplot as plt 
from math import *

class Wave(object):
	"""Classe contenant les attributs d'une onde harmonique sphérique
	s0: Amplitude de l'onde à 1m de la source
	w: pulsation de l'onde (radian)
	k: nombre d'onde (c² = w² / k²)"""
	def __init__(self,s0,w,k):
		self.s0 = s0
		self.w = w
		self.k = k

def get_wave_function(wave,x,y,t,dx=0,dy=0,opposite_phase=False):
	"""Retourne la fonction d'onde solution du rayonnement d'une source sphérique
	wave: object de type Wave, contenant s0,w,k
	x: coordonnée x
	y: coordonnée y
	t: temps
	dx: décalage en x de la source par rapport à l'origine
	dy: décalage en y de la source par rapport à l'origine
	opposite_phase: si True, inverse la phase de la fonction d'onde"""
	#calcul de la distance à la source
	r = np.asarray(np.sqrt(np.power(x-dx,2) + np.power(y-dy,2)),dtype=np.float)
	
	#calcul de l'amplitude de l'onde en fonction de sa phase
	if opposite_phase:
		return -((wave.s0/r)*np.exp(wave.k*r*1j - wave.w*t*1j)).astype('float')
	else:
		return ((wave.s0/r)*np.exp(wave.k*r*1j - wave.w*t*1j)).astype('float')
        
def animation(wave,timeline,vmin,vmax):
	"""Crée une série d'image
	wave: fonction retournant une amplitude à un temps t
	timeline: array de valeurs t
	vmin: valeur maximale affichée
	vmax: valeur minimale affichée"""
	for i,elm in enumerate(timeline):
		plt.imsave('image_{:06d}'.format(i),wave(elm),vmin=vmin,vmax=vmax,cmap='plasma')
 
#%%
#On crée deux types d'ondes
onde1 = Wave(1,2*pi*440,sqrt(2*pi*440/340.))
onde2 = Wave(1,2*pi*440,sqrt(2*pi*440/340.))
#%%
#On crée le maillage
X = np.linspace(-10,10,1000)
Y = np.linspace(-10,10,1000)
XX,YY = np.meshgrid(X,Y)
#%%
#On calcule les amplitudes des deux poles
pole1 = get_wave_function(onde1,XX,YY,1)
pole2 = get_wave_function(onde2,XX,YY,1,dx=-.25,opposite_phase=True)
#On les additionnes
bipole = pole1 + pole2
#On les affiche
plt.imshow(bipole,vmin=-.5,vmax=.5,cmap='plasma')
plt.colorbar()
plt.show()