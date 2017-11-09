#coding:utf-8
import numpy as np 
import matplotlib.pyplot as plt 
from math import *
from mpl_toolkits.mplot3d import Axes3D



class Wave(object):
	"""Classe contenant les attributs d'une onde harmonique sphérique
	s0: Amplitude de l'onde à 1m de la source
	w: pulsation de l'onde (radian)
	k: nombre d'onde (c² = w² / k²)"""
	def __init__(self,s0,w,k):
		self.s0 = s0
		self.w = w
		self.k = k

def get_wave_function(wave,x,y,t,dx=0,dy=0,dz=0,opposite_phase=False,\
	eq_plan=lambda x,y:0):
	"""Retourne la fonction d'onde solution du rayonnement d'une source sphérique
	wave: object de type Wave, contenant s0,w,k
	x: coordonnée x
	y: coordonnée y
	t: temps
	dx: décalage en x de la source par rapport à l'origine
	dy: décalage en y de la source par rapport à l'origine
	dz: décalage en z de la source par rapport à l'origine
	opposite_phase: si True, inverse la phase de la fonction d'onde
	eq_plan: fonction qui a x,y associe une valeur z"""
	#calcul de la distance à la source
	r = np.asarray(np.sqrt(np.power(x-dx,2) + np.power(y+dy,2) +\
				np.power(eq_plan(x-dx,y+dy)-dz,2)),dtype=np.float)
	
	#calcul de l'amplitude de l'onde en fonction de sa phase
	if opposite_phase:
		return -((wave.s0/r)*np.exp(wave.k*r*1j - wave.w*t*1j)).astype('float')
	else:
		return (np.real((wave.s0/r)*np.exp(wave.k*r*1j - wave.w*t*1j))).astype('float')
        
def animation(wave,timeline,vmin,vmax):
	"""Crée une série d'image
	wave: fonction retournant une amplitude à un temps t
	timeline: array de valeurs t
	vmin: valeur maximale affichée
	vmax: valeur minimale affichée"""
	for i,elm in enumerate(timeline):
		plt.imsave('image_{:06d}'.format(i),wave(elm),vmin=vmin,vmax=vmax,cmap='Blues')
		
def representation_physique(poles,X,Y,Z,simulation):
	"""Fonction représentant dans un espace 3d les résultats de la simulation.
	poles: position des poles
	X: mesh X
	Y: mesh Y
	Z: mesh Z
	simulation: mesh 2D contenant les valeurs d'amplitude"""
	
	simulation = abs(np.real(simulation))
	simulation = (100*simulation/(np.max(simulation))).astype('int')/100.
	
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	
	for elm in poles:
		ax.scatter(elm[0],elm[1],elm[2])
		
	
	
	colors = np.empty(simulation.shape,dtype=tuple)
	for i in range(len(simulation)):
		for j in range(len(simulation[0])):
			colors[i,j] = .5*simulation[i,j],.5*simulation[i,j],simulation[i,j]
	
	

	ax.plot_surface(X,Y,Z,facecolors=colors)
	plt.title("Representation de la simulation")
	plt.show()
	
	plt.imshow(simulation)
	plt.colorbar()
	plt.title("Rayonnement de la source sur le plan")
	plt.show()	
	
 
#%%
#On crée deux types d'ondes

#%%
#On crée le maillage

#%%
def simulation_2pole_hors_phase():
	#création du maillage
	X = np.linspace(-10,10,1000)
	Y = np.linspace(-10,10,1000)
	XX,YY = np.meshgrid(X,Y)
	#création de deux ondes
	onde1 = Wave(1,2*pi*440,sqrt(2*pi*440/340.))
	onde2 = Wave(1,2*pi*440,sqrt(2*pi*440/340.))
	#On calcule les amplitudes des deux poles
	pole1 = get_wave_function(onde1,XX,YY,1)
	pole2 = get_wave_function(onde2,XX,YY,1,dx=-.25,opposite_phase=True)
	#On les additionnes
	bipole = pole1 + pole2
	#On les affiche
	plt.imshow(bipole,vmin=-.5,vmax=.5,cmap='Blues')
	plt.colorbar()
	plt.title("2 poles en opposition de phase")
	plt.show()

#%%
def simulation_5poles(f,t):
    """Simule la pulsation de 5 poles à une fréquence f au temps t"""
    X = np.linspace(-1,9,1000)
    Y = np.linspace(-7,3,1000)
    XX,YY = np.meshgrid(X,Y)
    #création d'une onde
    onde = Wave(1,2*pi*f,sqrt(2*pi*f/340.))
    #création des 5 poles
    pole1 = get_wave_function(onde,XX,YY,t,dx=0, dy=0,opposite_phase=True)
    pole2 = get_wave_function(onde,XX,YY,t,dx=.15, dy=1,opposite_phase=True)
    pole3 = get_wave_function(onde,XX,YY,t,dx=.25, dy=2,opposite_phase=True)
    pole4 = get_wave_function(onde,XX,YY,t,dx=.15, dy=3,opposite_phase=True)  
    pole5 = get_wave_function(onde,XX,YY,t,dx=0, dy=4,opposite_phase=True)
	
    source = pole1 + pole2 + pole3 + pole4 + pole5
	
    return source

def simulation_poles_n_point(f,t,poles):
    onde = Wave(1,2*pi*f,sqrt(2*pi*f/340.))
    def sim(x,y):
        return sum([get_wave_function(onde,x,y,t,dx=elm[0],dy=elm[1]) for elm in poles])
    return sim
	
def simulation_n_poles(f,array,xmin,xmax,ymin,ymax,t,resolution,\
					eq_plan=lambda x,y:0):
	"""Simule une source composée de n monopoles
	f: fréquence en Hz d'oscillation
	array: liste de tuples contenant la position dans l'espace des poles
	xmin-xmax: positionement x de la fenetre d'observation
	ymin-ymax: positionement y de la fenetre d'obeservation
	t: temps utilisé pour la simulation
	resolution: resolution de la fenetre d'observation
	eq_plan= equation du plan de la fenetre d'observation"""
	#Création d'une onde
	onde = Wave(1,2*pi*f,sqrt(2*pi*f/340.))
	#Création du maillage
	XX,YY = np.meshgrid(np.linspace(xmin,xmax,resolution),np.linspace(ymin,\
	ymax,resolution))
	#Création et somme des poles
	source = sum([get_wave_function(onde,XX,YY,t,dx=elm[0],dy=elm[1],dz=elm[2],\
				eq_plan=eq_plan) for elm in array])
	return source
	
def intensity_over_time(simulation, duree, pas):
    """Donne une approximation de l'intensitée acoustique """
    space = np.linspace(0,duree,pas)
    intensite = np.zeros_like(simulation(0))
    for elm in space:
        intensite += simulation(elm)
    return intensite/len(space)
	