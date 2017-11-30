#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sympy
z = sympy.Symbol("z")
poly = sympy.poly
#%%
import pickle as pk
import numpy as np
import functools

def multiplie(liste):
    y = 1
    for elm in liste:
        y *= elm
    return y

def main(k,n,d):

	def get_root(i):
		return np.exp(-np.pi*1j*(-2*k*d*i/(float(n-1))))

	expression = multiplie([z-get_root(i) for i in np.arange(1,n+1)])

	coef = poly(expression,z).coeffs()
	with open("coef.txt","wb") as output:
            pk.dump([complex(elm) for elm in coef],output,0)
            print("ok")

main(340/5000.,11,.1)
