import pickle as pk
import numpy as np
def main(z,poly):
	k = 2/float(input("lambda? "))
	n = int(input("n? "))
	d = float(input("d? "))

	def get_root(i):
		return np.exp(-np.pi*1j*(-2*k*d*i/(float(n-1))))

	expression = reduce(lambda x,y:x*y,[z-get_root(i) for i in np.arange(1,6)])

	coef = poly(expression,z).coeffs()
	with open("coef.txt","w") as output:
		 pk.dump([complex(elm) for elm in coef],output,0)
