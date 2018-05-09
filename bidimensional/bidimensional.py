#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Aluno: Thales Carl Lavoratti (151000656)
Código do problema bidimensional de uma aleta
"""

import math as mt
import numpy as np
import matplotlib.pyplot as plt

##############################
#  Setting input paramethers
##############################

w = 1.0 #[m]
L = 1.0 #[m]
t = 0.1 #[m]
k = 230.0 #[W/mºC]
h = 26.3 #[W/m^2]
T0 = 100.0 #[ºC]
Tinf = 20.0 #[ºC]
yNumberOfNodes = 3
xNumberOfNodes = 4

#############################
# Mesh generation in y
#############################

yNodesPositions = []
deltaY = 0.5*t/(yNumberOfNodes - 1);
ySum = 0.0;
for i in range(yNumberOfNodes):
    yNodesPositions.append(ySum)
    ySum += deltaY
ySurfacePositions = []
ySum = 0.0
ySurfacePositions.append(ySum)
ySum = 0.5*deltaY
for i in range(yNumberOfNodes - 1):
    ySurfacePositions.append(ySum)
    ySum += deltaY
ySum -= 0.5*deltaY
ySurfacePositions.append(ySum)

#############################
# Mesh generation in x
#############################

xNodesPositions = []
deltaX = L/(float(xNumberOfNodes) - 0.5);
xSum = 0.5*deltaX;
for i in range(xNumberOfNodes):
    xNodesPositions.append(xSum)
    xSum += deltaX
xSurfacePositions = []
xSum = 0.0
for i in range(xNumberOfNodes):
    xSurfacePositions.append(xSum)
    xSum += deltaX
xSum -= 0.5*deltaX
xSurfacePositions.append(xSum)

##############################################################
# Generating the matrix of coefficients and independent vector
##############################################################
A = []
b = []
numberOfNodes = xNumberOfNodes * yNumberOfNodes

for i in range(numberOfNodes):
    A.append([])
    for j in range(numberOfNodes):
        A[i].append(0.0)
    b.append(0.0)
#bottom
for j in range(xNumberOfNodes):
    if  j == 0:
        aEast = k*w*deltaY/deltaX
        aWest = (k*w*deltaY)/(0.5*deltaX)
        aSouth = 0.0
        aNorth = k*w*deltaX/deltaY
        ap = aEast + aNorth + aWest + aSouth
        A[j][j] = ap
        A[j][j+1] = -aEast
        A[j][j+xNumberOfNodes] = -aNorth
        b[j] = aWest * T0
        
    elif j == (xNumberOfNodes - 1):
        aEast = 0.0
        aWest = k*w*deltaY/deltaX
        aSouth = 0.0
        aNorth = k*w*0.5*deltaX/deltaY
        ap = aEast + aNorth + aWest + aSouth
        A[j][j-1] = -aWest
        A[j][j] = ap
        A[j][j+xNumberOfNodes] = -aNorth
    else:
        aEast = k*w*deltaY/deltaX
        aWest = (k*w*deltaY)/(deltaX)
        aSouth = 0.0
        aNorth = k*w*deltaX/deltaY
        ap = aEast + aNorth + aWest + aSouth
        A[j][j-1] = -aWest
        A[j][j] = ap
        A[j][j+1] = -aEast        
        A[j][j+xNumberOfNodes] = -aNorth
        
        
#center
for i in range(1,yNumberOfNodes):
    for j in range(xNumberOfNodes):
        if  j == 0:
            print(i)
            print(j)
            print(" temperatura prescrita esquerda")
        elif j == (xNumberOfNodes - 1):
            print(i)
            print(j)
            print(" fluxo prescrito direita")
        else:
            print("normal")
#            aEast  =  k*deltaY/deltaX
#            aWest  =  k*deltaY/deltaX
#            aSouth =  k*deltaX/deltaY
#            aNorth =  k*deltaX/deltaY
#            ap = aEast + aWest + aSouth + aNorth
            
#top
for j in range(xNumberOfNodes):
    if  j == 0:
        print(j)
        print(" temperatura prescrita esquerda e convecção")
    elif j == (xNumberOfNodes - 1):
        print(j)
        print(" fluxo prescrito direita e convecção")
    else:
        print(j)
        print(" convecção")           