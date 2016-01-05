#!/usr/bin/env python
# -*- coding: utf-8 -*-
import math as mt
import numpy as np
import scipy as sp
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import pprint

def f0(x):
	if(abs(x-1.) <= 2.):
		return np.exp(-x**2 + 2*x)
	return 0

def l0(t):
	return 0
	
def r0(t):
	return 0

#Constantes D y b

#Diferencias finitas para EDP parabolica. (Nota: Xmin < Xmax)
#D: Constante D que acompaña al termino Uxx(x, y)
#b: Constante b que acompaña al termino Ux(x, y)
#Xmin: Limite inferior x
#Xmax: Limite superior x
#Tmax: Maximo tiempo a evolucionar la solucion
#dx: Distancia entre los puntos discretizados para el eje x
#dt: Distancia entre los puntos discretizados para el eje t
#f: Condicion inicial U(x, 0) = f(x)
#l: Condicion inicial U(Xmin, t) = l(t)
#r: condicion inicial U(Xmax, t) = r(t)
#Return: retorna una matriz de puntos aproximados y los intervalos para los ejes X e Y

def CrankNicholsonFiniteDifference(D, b, Xmin, Xmax, Tmax, dx, dt, f, l, r):
	Nx = int(mt.ceil((Xmax-Xmin)/dx))
	Nt = int(mt.ceil(Tmax/dt))

	sigma = D*dt/dx**2
	beta = b*dt/dx

	W = np.zeros([Nt+1, Nx+1])

	x = np.linspace(Xmin, Xmax, Nx+1)
	t = np.linspace(0, Tmax, Nt+1)

	#Matriz de coeficientes izquierdo
	A = np.zeros([Nx+1, Nx+1])
	for j in xrange(0, Nx+1):
		if(j > 0 and j < Nx):
			A[j][j-1] = -sigma
			A[j][j] = 2 + 2*sigma
			A[j][j+1] = -sigma
		else:
			A[j][j] = 1
	


	#Matriz de coeficientes derecho
	B = np.zeros([Nx+1, Nx+1])
	for j in xrange(0, Nx+1):
		if(j > 0 and j < Nx):
			B[j][j-1] = sigma - beta
			B[j][j] = 2 - 2*sigma
			B[j][j+1] = sigma + beta



	#Llenado vectores iniciales
	for i in xrange(0, Nx+1):
		W[0][i] = f(i*dx)

	for n in xrange(0, Nt):
		V = np.zeros([Nx+1, 1])
		V[0] = l((n+1)*dt)
		V[Nx] = r((n+1)*dt)
		c = np.dot(B, W[n].reshape([Nx+1, 1]))
		b = np.add(c, V)
		R = np.linalg.solve(A, b)
		W[n+1] = R.reshape([1, Nx+1])
	return W, x, t
		

w0, x0, t0 = CrankNicholsonFiniteDifference(4., -1., -4., 4.,  2., 0.5, 0.05, f0, l0, r0)
X0, T0 = np.meshgrid(x0, t0)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d', title="Grafico de la superficie de aproximacion")
ax.plot_wireframe(X0, T0, w0)
plt.show()
