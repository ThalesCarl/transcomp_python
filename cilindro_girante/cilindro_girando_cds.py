"""
Aluno: Thales Carl Lavoratti (151000656)
Código do problema do cilindro girante cds
"""

import math as mt
import numpy as np
import matplotlib.pyplot as plt
import auxiliar as aux

##############################
#  Setting input paramethers
##############################
yNumberOfNodes = 4
xNumberOfNodes = 4
w = 1.0 #[m]
re = 0.1 #[m]
ri = 0.04 #[m]
L =  mt.sin(0.25*mt.pi)*(re-ri)#[m]
k = 15.0#[W/mºC]
Ti = 250.0 #[ºC]
Te = 30.0 #[ºC]
omega = 2 #[rad/s]


#############################
# Mesh generation in y
#############################

yNodesPositions = []
deltaY = L/(yNumberOfNodes - 1);
ySum = mt.sin(0.25*mt.pi)*ri
for i in range(yNumberOfNodes):
    yNodesPositions.append(ySum)
    ySum += deltaY
ySurfacePositions = []
ySum = mt.sin(0.25*mt.pi)*ri + 0.5*deltaY
for i in range(yNumberOfNodes - 1):
    ySurfacePositions.append(ySum)
    ySum += deltaY
#############################
# Mesh generation in x
#############################

xNodesPositions = []
deltaX = L/(float(xNumberOfNodes) - 1);
xSum = mt.cos(0.25*mt.pi)*ri;
for i in range(xNumberOfNodes):
    xNodesPositions.append(xSum)
    xSum += deltaX
xSurfacePositions = []
xSum = mt.sin(0.25*mt.pi)*ri + 0.5*deltaX
for i in range(xNumberOfNodes - 1):
    xSurfacePositions.append(xSum)
    xSum += deltaX

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
    ap = 1.0
    A[j][j] = ap
    b[j] = aux.analyticSolution(xNodesPositions[j],yNodesPositions[0])        
        
#center
for i in range(1,yNumberOfNodes-1):
    for j in range(xNumberOfNodes):
        if  j == 0 or j == (xNumberOfNodes - 1):
            ap = 1.0
            A[i*xNumberOfNodes+j][i*xNumberOfNodes+j] = ap
            b[i*xNumberOfNodes+j] = aux.analyticSolution(xNodesPositions[j],yNodesPositions[i])
        else:
            uw = aux.uVelocity(omega,xSurfacePositions[j-1],yNodesPositions[i])
            ue = aux.uVelocity(omega,xSurfacePositions[j],yNodesPositions[i])
            vs = aux.vVelocity(omega,xNodesPositions[j],ySurfacePositions[i-1])
            vn = aux.vVelocity(omega,xNodesPositions[j],ySurfacePositions[i])
            aEast  =  k*w*deltaY/deltaX
            aWest  =  k*w*deltaY/deltaX
            aSouth =  k*w*deltaX/deltaY
            aNorth =  k*w*deltaX/deltaY
            ap = aEast + aWest + aSouth + aNorth
            A[i*xNumberOfNodes+j][j+xNumberOfNodes*(i-1)] = -aSouth
            A[i*xNumberOfNodes+j][i*xNumberOfNodes+j-1] = -aWest
            A[i*xNumberOfNodes+j][i*xNumberOfNodes+j] = ap
            A[i*xNumberOfNodes+j][i*xNumberOfNodes+j+1] = -aEast
            A[i*xNumberOfNodes+j][j+xNumberOfNodes*(i+1)] = -aNorth
#top
for j in range(xNumberOfNodes):    
    ap = 1.0          
    A[(yNumberOfNodes-1)*xNumberOfNodes+j][(yNumberOfNodes-1)*xNumberOfNodes+j] = ap
    b[(yNumberOfNodes-1)*xNumberOfNodes+j] = aux.analyticSolution(xNodesPositions[j],yNodesPositions[yNumberOfNodes - 1])

#############################
# Solving the linear system
############################
solution = np.linalg.solve(np.array(A),np.array(b))
temperatureField = []
for i in range(yNumberOfNodes):
    temperatureField.append([])
    for j in range(xNumberOfNodes):
        temperatureField[i].append(solution[i*xNumberOfNodes+j])

