#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 23 16:01:21 2018

@author: tclavoratti
"""
import math as mt
import solverTDMA as tdma
pi = mt.pi
##############################
#  Setting input paramethers
##############################

rho = 2700.0
cp = 900.0
k = 230
h = 20000
Tinf = 80
r0 = 0.05
wallLength = 0.005
length = 1
numberOfNodes = 5
prescribedTemperature = 200
deltaT = 0.001


#############################
# Mesh generation
#############################

nodePositions = []
surfacePositions = []
accumulator = r0
deltaR = wallLength/(numberOfNodes - 1)
for i in range(numberOfNodes):
    nodePositions.append(accumulator)
    accumulator += deltaR
accumulator = r0
surfacePositions.append(r0)
accumulator += 0.5*deltaR;
for i in range(numberOfNodes - 1):
    surfacePositions.append(accumulator)
    accumulator += deltaR
accumulator -= 0.5*deltaR;
surfacePositions.append(accumulator)

################################
# Setting beginTemperatureFields
################################
oldTemperatureField = []
for i in range(numberOfNodes):
    oldTemperatureField.append(80)


########################################
# Matrix coefficients and linear Vector
########################################

A=[]
b = []
A.append([])

#First Volume
A[0].append(1)
for i in range(numberOfNodes - 1):
    A[0].append(0)
b.append(prescribedTemperature)

#Central Volumes
for i in range(1,numberOfNodes - 1):
    ap0 = rho*cp*2*pi*nodePositions[i]*length * deltaR/deltaT
    aw = k*2*pi*surfacePositions[i] * length/deltaR
    ae = k*2*pi*surfacePositions[i+1] * length/deltaR
    ap = ap0 + aw + ae
    A.append([])
    for j in range(numberOfNodes):
        A[i].append(0)
    A[i][i-1] = -aw
    A[i][i] = ap
    A[i][i+1] = -ae
    b.append(ap0 * oldTemperatureField[i])
    
#Last Volume
ap0 = rho*cp*2*pi*0.5*nodePositions[numberOfNodes-1]*length * deltaR/deltaT
aw = 2*pi*k*surfacePositions[numberOfNodes-1]*length/deltaR
hAx = 2*pi*nodePositions[numberOfNodes-1]*h
ap = aw + ap0 + hAx
A.append([])
for j in range(numberOfNodes):
    A[numberOfNodes-1].append(0)
A[numberOfNodes-1][numberOfNodes-2] = -aw
A[numberOfNodes-1][numberOfNodes-1] = ap
b.append(ap0 * oldTemperatureField[numberOfNodes-1] + hAx*Tinf)   

temperatureField = tdma.solve(A,b)