#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 27 08:56:01 2018

@author: tclavoratti
"""

import math as mt
import numpy as np
import matplotlib.pyplot as plt
import solverTDMA as tdma
pi = mt.pi
def iterationsOnTransienteImplicito(numbOfNodes,deltaTime,tolerance):
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
    numberOfNodes = numbOfNodes
    initialTemperature = 80
    prescribedTemperature = 200
    deltaT = deltaTime
    temperatureFuncionPeriod = 0.04
    toleranceBetweenPeriods = tolerance
    
    
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
    
    ##########################################################
    # Setting beginTemperatureFields and iteration paramethers
    ##########################################################
    oldTemperatureField = []
    oldTemperatureField.append(200)
    for i in range(numberOfNodes - 1):
        oldTemperatureField.append(initialTemperature)
    
    timePosition = 0
    iterationCounter = 0
    transientTemperatureField = []
    transientTemperatureField.append([])
    transientTemperatureField[iterationCounter].append(timePosition)
    transientTemperatureField[iterationCounter].extend(oldTemperatureField)
    iterationPerPeriod = int(temperatureFuncionPeriod/deltaT)
    
    difference = 1
    
    while(difference>toleranceBetweenPeriods and iterationCounter < 10000):
        timePosition += deltaT
        iterationCounter += 1
        prescribedTemperature = 150 + 50*mt.cos(50*pi*timePosition)
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
        
        ###########################
        #Solving the linear system
        ##########################    
        temperatureField = tdma.solve(A,b)
        transientTemperatureField.append([])
        transientTemperatureField[iterationCounter].append(timePosition)
        transientTemperatureField[iterationCounter].extend(temperatureField)
        oldTemperatureField = temperatureField
        
        ###########################################
        #Comparing with the results a period before
        ###########################################
        soma = 0;
        if(iterationCounter>iterationPerPeriod):        
            for i in range(len(temperatureField)):
                tp = temperatureField[i]
                tp0 = transientTemperatureField[iterationCounter - iterationPerPeriod][i+1]
                auxiliar = tp -tp0
                soma += auxiliar*auxiliar
            difference = mt.sqrt(soma/numberOfNodes)
   
    return iterationCounter

#Variando o número de volumes de controle
x1 = []
y1 = []
for i in range(6):
    numberOfNodes = int(4*mt.pow(2,i))
    x1.append(numberOfNodes)
    y1.append(iterationsOnTransienteImplicito(numberOfNodes,0.001,0.01))
plt.plot(x1,y1,'ko',)
plt.xlabel('Número de volumes de controle')
plt.ylabel('Número de iterações')
plt.grid(b=1)
plt.savefig('./results_implicit/iterations_per_cv.png')
    
#Variando a tolerância
x2 = []
y2 = []
x3 = []
y3 = []
tolerance = 0.1
for i in range(5):
    x2.append(tolerance)
    x3.append(mt.log(tolerance,10))
    y2.append(iterationsOnTransienteImplicito(4,0.001,tolerance))
    y3.append(mt.log(y2[i]))
    tolerance /= 10

plt.cla()
plt.plot(x2,y2,'ko')
plt.xlabel('Tolerância')
plt.ylabel('Número de iterações')
plt.grid(b=1)
plt.savefig('./results_implicit/iterations_per_tolerance.png')

plt.cla()
plt.plot(x3,y2,'ko')
plt.xlabel('log(tolerância)')
plt.ylabel('Número de iterações')
plt.grid(b=1)
plt.savefig('./results_implicit/iterations_per_log_tolerance.png')