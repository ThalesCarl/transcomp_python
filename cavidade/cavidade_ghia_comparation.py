#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  6 13:59:10 2018

@author: thales
"""

import numpy as np
import auxiliar as aux
import matplotlib.pyplot as plt

def runCavidade(yNodes,method):
    ##############################
    #  Setting input paramethers
    ##############################
    xNodes = yNodes
    L = 1.0 
    topWallVelocity = 1.0 #[m/s]
    mi = 0.001 #[Pa * s]
    rho = 1.0 #[kg/m^3]
    
    
    timeStep = 30 #[s]
    numberOfTimeStepsMaximum = 60
    maxIterationCounter = 50
    tolerance = 1e-6
    
    #####################################
    #Setting the initial velocity Fields
    ####################################
    initial_U_Field = np.zeros((xNodes ,yNodes+1))
    initial_V_Field = np.zeros((xNodes+1,yNodes))
    initial_P_Field = np.zeros((xNodes,yNodes))
    
    uFields = []
    uFields.append(initial_U_Field)
    vFields = []
    vFields.append(initial_V_Field)
    pFields = []
    pFields.append(initial_P_Field)
    uCounter = (xNodes+1)*yNodes
    vCounter = xNodes*(yNodes+1)
    
    ###################
    #Appling the method
    ###################
    u0 = initial_U_Field
    v0 = initial_V_Field
    counterTime = 0
    differenceWithPast = 1.0
    while(differenceWithPast  > tolerance and counterTime < numberOfTimeStepsMaximum):
        uStar = u0
        vStar = v0
        difference = 1.0
        counterIteration = 0
        while(difference > tolerance and counterIteration < maxIterationCounter):
            A, b = aux.gen(xNodes,yNodes,L,topWallVelocity,timeStep,rho,mi,u0,uStar,v0,vStar,method)
            solution = np.linalg.solve(A,b)
            uLinear,vLinear,pLinear = np.split(solution,[uCounter,uCounter+vCounter])
            u = aux.fieldBuilder(uLinear,yNodes,xNodes+1)
            v = aux.fieldBuilder(vLinear,yNodes+1,xNodes)
            pressureField = aux.fieldBuilder(pLinear,yNodes,xNodes)        
            uDifference = abs(u - uStar)
            vDifference = abs(v - vStar)
            uMaximumDiference = uDifference.max()
            vMaximumDiference = vDifference.max()
            difference = max(uMaximumDiference,vMaximumDiference)
            uStar = u
            vStar = v
            counterIteration += 1
        print("sai da iteração")
        print(counterIteration)
        uFields.append(u)
        vFields.append(v) 
        pFields.append(pressureField)
        uDifferenceWithPast = abs(u-u0)
        vDifferenceWithPast = abs(v-v0)
        uMaximumDiferenceWithPast = uDifferenceWithPast.max()
        vMaximumDiferenceWithPast = vDifferenceWithPast.max()
        differenceWithPast = max(uMaximumDiferenceWithPast,vMaximumDiferenceWithPast)    
        u0 = u
        v0 = v
        counterTime += 1
    return u, v

yGhia = [0.0,0.0547,0.0625,0.0703,0.1016,0.1719,0.2813,0.4531,0.5,0.6172,0.7344,0.8516,0.9531,0.9609,0.9688,0.9766,1.0]
xGhia = [0.0,0.0625,0.0703,0.0781,0.0938,0.1563,0.2266,0.2344,0.5,0.8047,0.8594,0.9063,0.9453,0.9531,0.9609,0.9688,1.0]
uGhia = [0.0,-0.18109,-0.20196,-0.22220,-0.29730,-0.38289,-0.27805,-0.10648,-0.06080,0.05702,0.18719,0.33304,0.46604,0.51117,0.57492,0.65928,1.0]
vGhia = [0.0,0.27485,0.29012,0.30353,0.32627,0.37095,0.33075,0.32235,0.02526,-0.31966,-0.42665,-0.51550,-0.39188,-0.33714,-0.27669,-0.21388,0.0]
aGhia = [0.0,0.27485,0.29012,0.30353,0.32627,0.37095,0.33075,0.32235,0.02526,-0.31966,-0.42665,-0.51550,-0.39188,-0.33714,-0.27669,-0.21388,0.0]

#--------------------------------------------------------------UDS -------------------

nodes = 10
u10, v10 = runCavidade(nodes,1)
x10, y10 = aux.meshGenerator(1.0,nodes,nodes)
center10 = int(0.5*nodes)
u10 = u10.transpose()

nodes = 20
u20, v20 = runCavidade(nodes,1)
x20, y20 = aux.meshGenerator(1.0,nodes,nodes)
center20 = int(0.5*nodes)
u20 = u20.transpose()

nodes = 70
u70 = np.loadtxt('./results/uFields.csv',dtype=float,delimiter=',',skiprows=0)
v70 = np.loadtxt('./results/vFields.csv',dtype=float,delimiter=',',skiprows=0)
x70, y70 = aux.meshGenerator(1.0,nodes,nodes)
center70 = int(0.5*nodes)
u70 = u70.transpose()

plt.cla()
plt.plot(uGhia,yGhia,"ok",label = 'Ghia et al.')
plt.plot(u10[center10],y10,"+k",label = '10 volumes de controle')
plt.plot(u20[center20],y20,"*k",label = '20 volumes de controle')
plt.plot(u70[center70],y70,"-k",label = '70 volumes de controle')
legend = plt.legend(loc='lower right')
plt.show()

plt.cla()
plt.plot(vGhia,xGhia,"ok",label = 'Ghia et al.')
plt.plot(v10[center10],x10,"+k",label = '10 volumes de controle')
plt.plot(v20[center20],x20,"*k",label = '20 volumes de controle')
plt.plot(v70[center70],x70,"-k",label = '70 volumes de controle')
legend = plt.legend(loc='lower left')
plt.show()

#---------------------------------------------CDS ----------------------------
nodes = 10
u10, v10 = runCavidade(nodes,0)
x10, y10 = aux.meshGenerator(1.0,nodes,nodes)
center10 = int(0.5*nodes)
u10 = u10.transpose()

nodes = 20
u20, v20 = runCavidade(nodes,0)
x20, y20 = aux.meshGenerator(1.0,nodes,nodes)
center20 = int(0.5*nodes)
u20 = u20.transpose()

plt.cla()
plt.plot(uGhia,yGhia,"ok",label = 'Ghia et al.')
plt.plot(u10[center10],y10,"+k",label = '10 volumes de controle')
plt.plot(u20[center20],y20,"*k",label = '20 volumes de controle')
legend = plt.legend(loc='lower right')
plt.show()

plt.cla()
plt.plot(vGhia,xGhia,"ok",label = 'Ghia et al.')
plt.plot(v10[center10],x10,"+k",label = '10 volumes de controle')
plt.plot(v20[center20],x20,"*k",label = '20 volumes de controle')
legend = plt.legend(loc='lower left')
plt.show()