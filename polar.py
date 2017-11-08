#coding:utf-8
import numpy as np
import matplotlib.pyplot as plt
from math import *

class Wave(object):
    def __init__(self,s0,w,k):
        self.s0 = s0
        self.w = w
        self.k = k
        
class Pole(object):
    def __init__(self,x,y,opp_phase):
        self.x = x
        self.y = y
        self.opp_phase = opp_phase
    

def get_pressure(wave,x,y,t,dx=0,dy=0,dz=0,opposite_phase=False,\
                 eq_plan=lambda x,y:0):
    
    r = np.asarray(np.sqrt(np.power(x-dx,2) + np.power(y+dy,2) + np.power(eq_plan(x-dx,y+dy)-dz,2)),dtype=np.float)
    if opposite_phase:
        return -(np.real((wave.s0/r)*np.exp(wave.k*r*1j - wave.w*t*1j))).astype('float')
    else:
        return (np.real((wave.s0/r)*np.exp(wave.k*r*1j - wave.w*t*1j))).astype('float')
    
def get_intensity_around_circle(r,pressure,pas,period):
    """r distance d'évaluation de la pression efficace
    pressure fonction x,y,t -> pression instantanée
    pas: precision de la simulation, higher = better
    period : période d'oscillation de la source"""
    space = np.linspace(0,2*pi,1000)
    circlex = np.asarray([r*cos(theta) for theta in space])
    circley = np.asarray([r*sin(theta) for theta in space])
				
				
    intensity = sum([np.power(pressure(circlex,circley,t),2) for t in np.linspace(0,period,pas)])/float(pas)
    return space,np.sqrt(intensity)
    
#%%
    
if __name__ == "__main__":
    f=10000.    
    
    onde = Wave(1,2*pi*f,sqrt(2*pi*f/340.))
    
    poles = []
    poles.append(Pole(0,1.5,True))
    poles.append(Pole(0,.5,True))
    poles.append(Pole(0,-.5,True))
    poles.append(Pole(0,-1.5,True))
    
    poles.append(Pole(-.2,1.5,False))
    poles.append(Pole(-.2,.5,False))
    poles.append(Pole(-.2,-.5,False))
    poles.append(Pole(-.2,-1.5,False))
    
    get_source_pressure = lambda x,y,t: sum([get_pressure(onde,x,y,t,p.x,p.y, opposite_phase=p.opp_phase) for p in poles])
    
    #PLOT DE LA SOURCE

    x = np.linspace(-5,5,1000)
    y = np.linspace(-5,5,1000)
    xx,yy = np.meshgrid(x,y)
    source = get_source_pressure(xx,yy,1)
    plt.imshow(source,vmin=-5,vmax=5,cmap='Blues')
    plt.title("Representation 2D de la source pour $f={}$".format(f))
    plt.show()
    
    #PLOT DE LA DIRECTIVITE DE LA SOURCE    
    
    intensity = get_intensity_around_circle(3,get_source_pressure,200,1/f)
    plt.polar(*intensity)
    plt.savefig("polar_bipole.png")
    plt.title("Directivite de la source pour $f={}$".format(f))
    plt.show()
    
    
    

    
    