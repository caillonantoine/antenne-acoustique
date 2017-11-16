#coding:utf-8
import numpy as np
import matplotlib.pyplot as plt
from math import *

class Wave(object):
    """Définit une onde, comportant s0, sa pulsation, son nombre d'onde k"""
    def __init__(self,s0,w,k):
        self.s0 = s0
        self.w = w
        self.k = k
        
class Pole(object):
    """Définit un monopole par sa position x,y 
    et sa potentielle opposition de phase"""
    def __init__(self,x,y,opp_phase,phase=0,ponderation=1):
        self.x = x
        self.y = y
        self.opp_phase = opp_phase
        self.phase  = phase
        self.ponderation = ponderation
    

def get_pressure(wave,x,y,t,dx=0,dy=0,opposite_phase=False,phase=0,ponderation=1):
    """renvoie une pression instantanée par rapport à:
    -une onde (class Wave)
    -un pole (class Pole), dx,dy,opposition de phase
    -un point d'évalution x,y,t"""
    r = np.asarray(np.sqrt(np.power(x-dx,2) + np.power(y+dy,2)),dtype=np.float)
    if opposite_phase:
        return -ponderation*(np.real((wave.s0/r)*np.exp(wave.k*r*1j -\
        (wave.w + phase)*t*1j))).astype('float')
    else:
        return ponderation*(np.real((wave.s0/r)*np.exp(wave.k*r*1j -\
        (wave.w + phase)*t*1j))).astype('float')
    
def get_intensity_around_circle(r,pressure,pas,period):
    """r distance d'évaluation de l'intensité
    pressure: fonction x,y,t -> pression instantanée
    pas: precision de la simulation, higher = better
    period : période d'oscillation de la source"""
    space = np.linspace(0,2*pi,1000)
    circlex = np.asarray([r*cos(theta) for theta in space])
    circley = np.asarray([r*sin(theta) for theta in space])
    #On somme les pressions pour tous les intervalles de temps
    intensity = sum([np.power(pressure(circlex,circley,t),2)\
    for t in np.linspace(0,period,pas)])/float(pas)
    return space,np.sqrt(intensity)
    
#%%
    
if __name__ == "__main__":
    f=500. #Définition de la fréquence
    longueur_onde = 340/f
    
    onde = Wave(1,2*pi*f,sqrt(2*pi*f/340.)) #Définition d'une classe d'onde
    
    poles = [] #Liste de poles
    
    #Ajout des monopoles
    for k in np.arange(-3,3,1):
        poles.append(Pole(0,k*10e-2,True))
    
    #On définit une fonction donnant la pression en x,y,t
    get_source_pressure = lambda x,y,t: sum([get_pressure(onde,x,y,t,p.x,p.y,\
    opposite_phase=p.opp_phase,phase=p.phase,ponderation=p.ponderation) for p in poles])
    
    zoom_factor = 1
    #On représente la source d'émissions
    x = np.linspace(-zoom_factor,zoom_factor,1000) #SCALE X
    y = np.linspace(-zoom_factor,zoom_factor,1000) #SCALE Y
    xx,yy = np.meshgrid(x,y) #SCALE XY
    
    source = get_source_pressure(xx,yy,1) #On évalue la pression instantanée pour t=1
    #PLOT DE LA SOURCE
    plt.imshow(source,cmap='Blues',vmin=-50,vmax=50,extent=[x[0],x[-1],y[0],y[-1]])
    plt.title("Representation 2D de la source pour $f={}$".format(f))
    plt.xlabel("$x$")
    plt.ylabel("$y$")
    plt.show()
    
    #PLOT DE LA DIRECTIVITE DE LA SOURCE    
    #On récupère les valeurs d'intensité
    intensity = get_intensity_around_circle(10,get_source_pressure,1000,1/f)
    plt.polar(*intensity,color='orange')
    #%%
    #PLOT DE LA DIRECTIVITE THEORIQUE D'UN DIPOLE
    space = np.linspace(0,2*pi,1000)
    directivite = abs(np.cos(space))
    plt.polar(space,directivite,color='blue')
    plt.legend(["Mesured","Theoric"])
    plt.savefig("polarplot.eps")
    plt.show()
    
    

    
    