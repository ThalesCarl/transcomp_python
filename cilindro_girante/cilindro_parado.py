"""
Aluno: Thales Carl Lavoratti (151000656)
Código do problema bidimensional de uma aleta
"""

import math as mt
import numpy as np
import matplotlib.pyplot as plt
import analyticSolution as soluc

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


#############################
# Mesh generation in y
#############################

yNodesPositions = []
deltaY = L/(yNumberOfNodes - 1);
ySum = mt.sin(0.25*mt.pi)*ri;
for i in range(yNumberOfNodes):
    yNodesPositions.append(ySum)
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
    b[j] = soluc.analyticSolution(xNodesPositions[j],yNodesPositions[0])
        
        
#center
for i in range(1,yNumberOfNodes-1):
    for j in range(xNumberOfNodes):
        if  j == 0:
            ap = 1.0
            A[i*xNumberOfNodes+j][i*xNumberOfNodes+j] = ap
            b[i*xNumberOfNodes+j] = soluc.analyticSolution(xNodesPositions[j],yNodesPositions[i])
        elif j == (xNumberOfNodes - 1):
            ap = 1.0
            A[i*xNumberOfNodes+j][i*xNumberOfNodes+j] = ap
            b[i*xNumberOfNodes+j] = soluc.analyticSolution(xNodesPositions[j],yNodesPositions[i])
        else:
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
##top
#for j in range(xNumberOfNodes):
#    
#    if  j == 0:
#        aEast = k*w*0.5*deltaY/deltaX
#        aWest = k*w*0.5*deltaY/(0.5*deltaX)
#        aSouth = k*w*deltaX/deltaY
#        aNorth = h*w*deltaX
#        ap = aEast + aNorth + aWest + aSouth
#        A[(yNumberOfNodes-1)*xNumberOfNodes+j][j+(yNumberOfNodes-2)*xNumberOfNodes] = -aSouth
#        A[(yNumberOfNodes-1)*xNumberOfNodes+j][(yNumberOfNodes-1)*xNumberOfNodes+j] = ap
#        A[(yNumberOfNodes-1)*xNumberOfNodes+j][(yNumberOfNodes-1)*xNumberOfNodes+j+1] = -aEast
#        b[(yNumberOfNodes-1)*xNumberOfNodes+j] = aWest * T0 + aNorth*Tinf
#    elif j == (xNumberOfNodes - 1):
#        aEast =  0.0
#        aWest =  k*w*0.5*deltaY/deltaX
#        aSouth = k*w*0.5*deltaX/deltaY
#        aNorth = h*w*0.5*deltaX
#        ap = aEast + aNorth + aWest + aSouth          
#        A[(yNumberOfNodes-1)*xNumberOfNodes+j][j+(yNumberOfNodes-2)*xNumberOfNodes] = -aSouth
#        A[(yNumberOfNodes-1)*xNumberOfNodes+j][(yNumberOfNodes-1)*xNumberOfNodes+j-1] = -aWest
#        A[(yNumberOfNodes-1)*xNumberOfNodes+j][(yNumberOfNodes-1)*xNumberOfNodes+j] = ap
#        b[(yNumberOfNodes-1)*xNumberOfNodes+j] = aNorth*Tinf        
#    else:
#        aEast =  k*w*0.5*deltaY/deltaX
#        aWest =  k*w*0.5*deltaY/deltaX
#        aSouth = k*w*deltaX/deltaY
#        aNorth = h*w*deltaX
#        ap = aEast + aNorth + aWest + aSouth          
#        A[(yNumberOfNodes-1)*xNumberOfNodes+j][j+(yNumberOfNodes-2)*xNumberOfNodes] = -aSouth
#        A[(yNumberOfNodes-1)*xNumberOfNodes+j][(yNumberOfNodes-1)*xNumberOfNodes+j-1] = -aWest
#        A[(yNumberOfNodes-1)*xNumberOfNodes+j][(yNumberOfNodes-1)*xNumberOfNodes+j] = ap
#        A[(yNumberOfNodes-1)*xNumberOfNodes+j][(yNumberOfNodes-1)*xNumberOfNodes+j+1] = -aEast
#        b[(yNumberOfNodes-1)*xNumberOfNodes+j] = aNorth*Tinf
##############################
## Solving the linear system
#############################
#solution = np.linalg.solve(np.array(A),np.array(b))
#temperatureField = []
#for i in range(yNumberOfNodes):
#    temperatureField.append([])
#    for j in range(xNumberOfNodes):
#        temperatureField[i].append(solution[i*xNumberOfNodes+j])
#    
##############################################
##Mirror of the temperature field and the mesh
##############################################
#auxTemp = []
#auxMesh = []
#for i in range(1,yNumberOfNodes):
#    auxTemp.append(temperatureField[i])
#    auxMesh.append(-yNodesPositions[i])
#auxTemp = np.flip(auxTemp,0)
#auxMesh = np.flip(auxMesh,0)  
#for i in range(len(auxTemp)):
#    temperatureField.insert(i,auxTemp[i])
#    yNodesPositions.insert(i,auxMesh[i])
#temperatureField = np.array(temperatureField)
#yNodesPositions = np.array(yNodesPositions)
#    
#    
###################################
##Plotting the solution
#########################
#xx, yy = np.meshgrid(xNodesPositions,yNodesPositions)
#plt.contourf(xx,yy,np.array(temperatureField))
#plt.colorbar(orientation="vertical")
#plt.xlabel("x[m]")
#
