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
    
def get_intensity_around_circle(r,pressure,duree,pas):
    space = np.linspace(0,2*pi,1000)
    circlex = np.asarray([r*cos(theta) for theta in space])
    circley = np.asarray([r*sin(theta) for theta in space])
    
    intensity = sum([pressure(circlex,circley,t) for t in np.linspace(0,duree,pas)])
    return space,abs(intensity)
    
    
if __name__ == "__main__":
    onde = Wave(1,2*pi*440,sqrt(2*pi*440/340.))
    
    poles = []
    poles.append(Pole(0,0,False))
    poles.append(Pole(-.1,0,True))
    
    x = np.linspace(-10,10,1000)
    y = np.linspace(-10,10,1000)
    xx,yy = np.meshgrid(x,y)
    
    get_source_pressure = lambda x,y,t: sum([get_pressure(onde,x,y,t,p.x,p.y, opposite_phase=p.opp_phase) for p in poles])
    
    intensity = get_intensity_around_circle(1,get_source_pressure,.5,100)
    plt.polar(*intensity)
    

    
    