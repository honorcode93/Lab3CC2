#!/usr/bin/env python
# -*- coding: utf-8 -*-
import math as mt
import numpy as np
import scipy as sp
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
import pprint

def DeltaKronecker(k, m):
	if k == m:
		return 1.
	return 0.

#Retorna la siguiente matriz
#T0: Matriz en el tiempo n-1
#T1: Matriz en el tiempo n
#c: Velocidad sonido
#dt: delta t
#h: delta h
#Xmin: Minimo x
#Xmax: Maximo x
#Ymin: Minimo y
#Ymax: Maximo y
#w: Frecuencia angular
#p: Valor de p para el delta Kronecker
#q: Valor de q para el delta Kronecker
#n: numero de iteracion
#F: Function F(t)
def GetNextMatrix(T0, T1, c, dt, h, Xmin, Xmax, Ymin, Ymax, wa, wb, p, q, m, n, it, Ua, Ub):
	sigma = c*dt/h
	Nx = int(mt.ceil((Xmax-Xmin)/h))
	Ny = int(mt.ceil((Ymax-Ymin)/h))

	#Inicializar siguiente matriz
	T2 = np.zeros((Ny+1, Nx+1))
	for i in xrange(0, Nx):
		for j in xrange(0, Ny):

			#Lado izquierdo
			if i == 0:
				T2[j][0] = T1[j][0] + T1[j][1] - T0[j][1] + sigma*(T1[j][1] - T1[j][0]-(T0[j][2]-T0[j][1]))

			#Lado derecho
			elif i == Nx:
				T2[j][Nx] = T1[j][Nx] + T1[j][Nx-1] - T0[j][Nx-1] + sigma*(T1[j][Nx-1] - T1[j][Nx] - (T0[j][Nx-2] - T0[j][Nx-1]))
			#Lado inferior
			elif j == 0:
				T2[0][i] = T1[0][i] + T1[1][i] - T0[1][i] + sigma*(T1[1][i] - T1[0][i] - (T0[2][i] - T0[1][i]))
			#Lado superior
			elif j == Ny:
				T2[Ny][i] = T1[Ny][i] + T1[Ny-1][i] - T0[Ny-1][i] + sigma(T1[Ny-1][i] - T1[Ny][i] - (T0[Ny-2][i] - T0[Ny-1][i]))
			#Discretizacion propia
			else:
				T2[j][i] = sigma**2*(T1[j][i-1] + T1[j][i+1] + T1[j-1][i] + T1[j+1][i]) + (2-4*sigma**2)*T1[j][i] - T0[j][i] + Ua*wa**2*np.cos(wa*dt*it)*DeltaKronecker(i, p)*DeltaKronecker(j, q) + Ub*wb**2*np.cos(wb*dt*it)*DeltaKronecker(i, m)*DeltaKronecker(j, n)
	return T2

#Realiza la simulacion
#fa: Frecuencia de la fuente A
#fb: Frecuencia de la fuente B
#c: Velocidad del sonido
#dt: Delta t
#Xmin: Minimo x
#Xmax: Maximo x
#Ymin: Minimo y
#Ymax: Maximo y
#p: Posicion x de la fuente A
#q: Posicion y de la fuente A
#m: Posicion x de la fuente B
#n: Posicion y de la fuente B
#Ua: Presion inicial de la fuente A
#Ub: Presion inicial de la fuente B
def Simulacion(fa, fb, c, h, dt, Xmin, Xmax, Ymin, Ymax, p, q, m, n, Ua, Ub):
	Nx = int(mt.ceil((Xmax-Xmin)/h))
	Ny = int(mt.ceil((Ymax-Ymin)/h))
	T0 = np.zeros((Ny+1, Nx+1))
	T1 = np.zeros((Ny+1, Nx+1))
	#T0 = np.full((Ny+1, Nx+1), U0)
	#T1 = np.full((Ny+1, Nx+1), U0)
	sigma = c*dt/h
	
	if(sigma >= 1/mt.sqrt(2)):
		print "Metodo inestable con el valor actual de sigma: "+str(sigma)
		return
	

	#Frecuencia angular
	wa = 2*np.pi*fa
	wb = 2*np.pi*fb

	#Equiespaciado
	x = np.linspace(Xmin, Xmax, Nx+1)
	y = np.linspace(Ymin, Ymax, Ny+1)

	#Plot
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	ax.set_zlim3d(-2*10**9, 2*10**9)
	X, Y = np.meshgrid(x, y)

	max_iters = int(mt.ceil(1/dt))

	wire = ax.plot_wireframe(X, Y, T0)
	#surf = ax.plot_surface(X, Y, T0, cstride=h, rstride=h, cmap=cm.jet)

	for i in xrange(0, max_iters):
		print "t = "+str(i*dt)
		T2 = GetNextMatrix(T0, T1, c, dt, h, Xmin, Xmax, Ymin, Ymax, wa, wb, p, q, m, n, i, Ua, Ub)
		T0 = T1
		T1 = T2

		wire.remove()
		#surf.remove()
		wire = ax.plot_wireframe(X, Y, T0)
		#surf = ax.plot_surface(X, Y, T0, cstride=h, rstride=h, cmap=cm.jet)
	
		plt.pause(dt)
	plt.close()

#Simulacion(fa, fb, c, h, dt, Xmin, Xmax, Ymin, Ymax, p, q, m, n, Ua, Ub):
Simulacion(1000, 1000, 1500, 1, 0.0001, 0, 200, 0, 200, 100, 100, 100, 50, 50, 50)


