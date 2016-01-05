#!/usr/bin/env python
# -*- coding: utf-8 -*-
import math as mt
import numpy as np
import scipy as sp
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import pprint

#Funcion de prueba
def F0(x, y):
	return x
def Z0(x, y):
	return 0
def up0(x):
	return np.sin(np.pi*x)
def down0(x):
	return np.sin(np.pi*x)
def left0(y):
	return 0
def right0(y):
	return 0

#Funcion 1
def F1(x, y):
	return np.sin(np.pi*x*y)

def Z1(x, y):
	return x**2+y**2

def up1(x):
	return 0

def down1(x):
	return 0

def left1(y):
	return 0

def right1(y):
	return 0

#Funcion 2
def F2(x, y):
	return np.cos(np.pi*x*y) + np.exp(2*x*y)

def Z2(x, y):
	return 0

def up2(x):
	return 0

def down2(x):
	return 0

def left2(y):
	return 0

def right2(y):
	return 0

#Funcion 3
def F3(x, y):
	return x

def Z3(x, y):
	return 0

def up3(x):
	return np.cos(np.pi*x)

def down3(x):
	return np.sin(np.pi*x)

def left3(y):
	return 0

def right3(y):
	return 0

#Diferencias finitas para EDP eliptica. (Nota: Xmin < Xmax, Ymin < Ymax)
#F: Funcion F(x, y) que esta a la derecha de la ecuacion.
#Z: Funcion Z(x, y)  que acompana al termino u(x, y) en la ecuacion (si es de helmholtz)
#Xmin: Limite inferior x
#Xmax: Limite superior x
#Ymin: Limite inferior y
#Ymax: Limite superior y
#Nx: Puntos requeridos eje x (h = (Xmax-Xmin)/Nx)
#Ny: Puntos requeridos eje y (k = (Ymax-Ymin)/Ny)
#up: condicion inicial arriba u(x, Ymax) = up(x)
#down: condicion inicial abajo u(x, Ymin) = down(x)
#left: condicion inicial izquierda u(Xmin, y) = left(y)
#right: conidicion inicial derecha u(Xmax, y) = right(y)
#
#Return: retorna una matriz de puntos aproximados y los intervalos para los ejes X e Y
def ElipticFiniteDifference(F, Z, Xmin, Xmax, Ymin, Ymax, Nx, Ny, up, down, left, right):
	#Discretizacion en x e y, o bien, creacion de particiones.
	x = np.linspace(Xmin, Xmax, Nx+1)
	y = np.linspace(Ymin, Ymax, Ny+1)

	#Definicion de step size, o bien, distancia entre los puntos
	#de la discretizacion.
	H = x[1]-x[0]
	K = y[1]-y[0]

	#Crear matrices del sistema lineal.
	mn = (Nx+1)*(Ny+1)
	A = np.zeros([mn, mn])
	b = np.zeros([mn, 1])

	#Matriz para la obtencion de los indices
	def index(i, j, m=Nx+1):
		return j + i*m

	#Diferencias finitas, llenado de la matriz
	for i in xrange(Nx+1):
		for j in xrange(Ny+1):

			#Obtencion del indice actual.
			k = index(i, j)

			#Esta en la frontera up
			if j == Ny:
				A[k, k] = 1
				b[k] = up(x[i])

			#Esta en la frontera down
			elif j == 0:
				A[k, k] = 1
				b[k] = down(x[i])

			#Esta en la frontera left
			elif i == 0:
				A[k, k] = 1
				b[k] = left(y[j])

			#Esta en la frontera right
			elif i == Nx:
				A[k, k] = 1
				b[k] = right(y[j])
			else:
				#Resto de los puntos desconocidos
				A[k, k] = -2/H**2 -2/K**2 -Z(x[i], y[j])
				A[k, index(i+1, j)] = 1/H**2
				A[k, index(i-1, j)] = 1/H**2
				A[k, index(i, j-1)] = 1/K**2
				A[k, index(i, j+1)] = 1/K**2
				b[k] = F(x[i], y[j])
	w = np.linalg.solve(A, b)
	return w, x, y

#Se debe hacer para n = 10, 100 y 20


#Esto es un test
'''
wt, xt, yt = ElipticFiniteDifference(F0, Z0, 0., 1., 0., 1., 20, 20, up0, down0, left0, right0)
Xt, Yt = np.meshgrid(xt, yt)
Wt = wt.reshape(Xt.shape)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d', title="Funcion de prueba")
ax.plot_wireframe(Xt, Yt, Wt)
plt.show()
'''

#F1

w1a, x1a, y1a = ElipticFiniteDifference(F1, Z1, 0., 1., 0., 1., 10, 10, up1, down1, left1, right1)
w1b, x1b, y1b = ElipticFiniteDifference(F1, Z1, 0., 1., 0., 1., 100, 100, up1, down1, left1, right1)
w1c, x1c, y1c = ElipticFiniteDifference(F1, Z1, 0., 1., 0., 1., 20, 20, up1, down1, left1, right1)

#F2
w2a, x2a, y2a = ElipticFiniteDifference(F2, Z2, 0., 1., 0., 1., 10, 10, up2, down2, left2, right2)
w2b, x2b, y2b = ElipticFiniteDifference(F2, Z2, 0., 1., 0., 1., 100, 100, up2, down2, left2, right2)
w2c, x2c, y2c = ElipticFiniteDifference(F2, Z2, 0., 1., 0., 1., 20, 20, up2, down2, left2, right2)


#F3
w3a, x3a, y3a = ElipticFiniteDifference(F3, Z3, 0., 1., 0., 1., 10, 10, up3, down3, left3, right3)
w3b, x3b, y3b = ElipticFiniteDifference(F3, Z3, 0., 1., 0., 1., 100, 100, up3, down3, left3, right3)
w3c, x3c, y3c = ElipticFiniteDifference(F3, Z3, 0., 1., 0., 1., 20, 20, up3, down3, left3, right3)

#Para graficar
X1a, Y1a = np.meshgrid(x1a, y1a)
X1b, Y1b = np.meshgrid(x1b, y1b)
X1c, Y1c = np.meshgrid(x1c, y1c)
X2a, Y2a = np.meshgrid(x2a, y2a)
X2b, Y2b = np.meshgrid(x2b, y2b)
X2c, Y2c = np.meshgrid(x2c, y2c)
X3a, Y3a = np.meshgrid(x3a, y3a)
X3b, Y3b = np.meshgrid(x3b, y3b)
X3c, Y3c = np.meshgrid(x3c, y3c)


#Convertir lista de puntos a arreglo
W1a = w1a.reshape(X1a.shape)
W1b = w1b.reshape(X1b.shape)
W1c = w1c.reshape(X1c.shape)
W2a = w2a.reshape(X2a.shape)
W2b = w2b.reshape(X2b.shape)
W2c = w2c.reshape(X2c.shape)
W3a = w3a.reshape(X3a.shape)
W3b = w3b.reshape(X3b.shape)
W3c = w3c.reshape(X3c.shape)

#Grafico F1a
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d', title="Funcion 1 (h = 0.1)")
ax.plot_wireframe(X1a, Y1a, W1a)
plt.show()

#Grafico F1b
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d', title="Funcion 1 (h = 0.01)")
ax.plot_wireframe(X1b, Y1b, W1b)
plt.show()

#Grafico F1c
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d', title="Funcion 1 (h = 0.05)")
ax.plot_wireframe(X1c, Y1c, W1c)
plt.show()

#Grafico F2a
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d', title="Funcion 2 (h = 0.1)")
ax.plot_wireframe(X2a, Y2a, W2a)
plt.show()

#Grafico F2b
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d', title="Funcion 2 (h = 0.01)")
ax.plot_wireframe(X2b, Y2b, W2b)
plt.show()

#Grafico F2c
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d', title="Funcion 2 (h = 0.05)")
ax.plot_wireframe(X2c, Y2c, W2c)
plt.show()

#Grafico F3a
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d', title="Funcion 3 (h = 0.1)")
ax.plot_wireframe(X3a, Y3a, W3a)
plt.show()

#Grafico F3b
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d', title="Funcion 3 (h = 0.01)")
ax.plot_wireframe(X3b, Y3b, W3b)
plt.show()

#Grafico F3c
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d', title="Funcion 3 (h = 0.05)")
ax.plot_wireframe(X3c, Y3c, W3c)
plt.show()





