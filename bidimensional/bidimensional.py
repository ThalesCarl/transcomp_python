"""
Aluno: Thales Carl Lavoratti (151000656)
Código do problema bidimensional de uma aleta
"""

import math as mt
import numpy as np
import matplotlib.pyplot as plt
def getTemperatureField(xN,yN):
    ##############################
    #  Setting input paramethers
    ##############################
    
    w = 1.0 #[m]
    L = 1.2 #[m]
    t = 0.15 #[m]
    k = 230 #[W/mºC]
    h = 27.4 #[W/m^2ºC]
    T0 = 80.0 #[ºC]
    Tinf = 10.0 #[ºC]
    yNumberOfNodes = int(0.5*(1+yN))
    xNumberOfNodes = xN
    
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
            aEast = k*w*0.5*deltaY/deltaX
            aWest = (k*w*0.5*deltaY)/(0.5*deltaX)
            aSouth = 0.0
            aNorth = k*w*deltaX/deltaY
            ap = aEast + aNorth + aWest + aSouth
            A[j][j] = ap
            A[j][j+1] = -aEast
            A[j][j+xNumberOfNodes] = -aNorth
            b[j] = aWest * T0
            
        elif j == (xNumberOfNodes - 1):
            aEast = 0.0
            aWest = k*w*0.5*deltaY/deltaX
            aSouth = 0.0
            aNorth = k*w*0.5*deltaX/deltaY
            ap = aEast + aNorth + aWest + aSouth
            A[j][j-1] = -aWest
            A[j][j] = ap
            A[j][j+xNumberOfNodes] = -aNorth
        else:
            aEast = k*w*0.5*deltaY/deltaX
            aWest = k*w*0.5*deltaY/deltaX
            aSouth = 0.0
            aNorth = k*w*deltaX/deltaY
            ap = aEast + aNorth + aWest + aSouth
            A[j][j-1] = -aWest
            A[j][j] = ap
            A[j][j+1] = -aEast        
            A[j][j+xNumberOfNodes] = -aNorth
            
            
    #center
    for i in range(1,yNumberOfNodes-1):
        for j in range(xNumberOfNodes):
            if  j == 0:
                aEast =  k*w*deltaY/deltaX
                aWest =  k*w*deltaY/(0.5*deltaX)
                aSouth = k*w*deltaX/deltaY
                aNorth = k*w*deltaX/deltaY
                ap = aEast + aNorth + aWest + aSouth
                A[i*xNumberOfNodes+j][j+xNumberOfNodes*(i-1)] = -aSouth
                A[i*xNumberOfNodes+j][i*xNumberOfNodes+j] = ap
                A[i*xNumberOfNodes+j][i*xNumberOfNodes+j+1] = -aEast
                A[i*xNumberOfNodes+j][j+xNumberOfNodes*(i+1)] = -aNorth
                b[i*xNumberOfNodes+j] = aWest * T0
            elif j == (xNumberOfNodes - 1):
                aEast  =  0.0
                aWest  =  k*w*deltaY/deltaX
                aSouth =  k*w*0.5*deltaX/deltaY
                aNorth =  k*w*0.5*deltaX/deltaY
                ap = aEast + aWest + aSouth + aNorth
                A[i*xNumberOfNodes+j][j+xNumberOfNodes*(i-1)] = -aSouth
                A[i*xNumberOfNodes+j][i*xNumberOfNodes+j-1] = -aWest
                A[i*xNumberOfNodes+j][i*xNumberOfNodes+j] = ap
                A[i*xNumberOfNodes+j][j+xNumberOfNodes*(i+1)] = -aNorth
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
    #top
    for j in range(xNumberOfNodes):
        
        if  j == 0:
            aEast = k*w*0.5*deltaY/deltaX
            aWest = k*w*0.5*deltaY/(0.5*deltaX)
            aSouth = k*w*deltaX/deltaY
            aNorth = h*w*deltaX
            ap = aEast + aNorth + aWest + aSouth
            A[(yNumberOfNodes-1)*xNumberOfNodes+j][j+(yNumberOfNodes-2)*xNumberOfNodes] = -aSouth
            A[(yNumberOfNodes-1)*xNumberOfNodes+j][(yNumberOfNodes-1)*xNumberOfNodes+j] = ap
            A[(yNumberOfNodes-1)*xNumberOfNodes+j][(yNumberOfNodes-1)*xNumberOfNodes+j+1] = -aEast
            b[(yNumberOfNodes-1)*xNumberOfNodes+j] = aWest * T0 + aNorth*Tinf
        elif j == (xNumberOfNodes - 1):
            aEast =  0.0
            aWest =  k*w*0.5*deltaY/deltaX
            aSouth = k*w*0.5*deltaX/deltaY
            aNorth = h*w*0.5*deltaX
            ap = aEast + aNorth + aWest + aSouth          
            A[(yNumberOfNodes-1)*xNumberOfNodes+j][j+(yNumberOfNodes-2)*xNumberOfNodes] = -aSouth
            A[(yNumberOfNodes-1)*xNumberOfNodes+j][(yNumberOfNodes-1)*xNumberOfNodes+j-1] = -aWest
            A[(yNumberOfNodes-1)*xNumberOfNodes+j][(yNumberOfNodes-1)*xNumberOfNodes+j] = ap
            b[(yNumberOfNodes-1)*xNumberOfNodes+j] = aNorth*Tinf        
        else:
            aEast =  k*w*0.5*deltaY/deltaX
            aWest =  k*w*0.5*deltaY/deltaX
            aSouth = k*w*deltaX/deltaY
            aNorth = h*w*deltaX
            ap = aEast + aNorth + aWest + aSouth          
            A[(yNumberOfNodes-1)*xNumberOfNodes+j][j+(yNumberOfNodes-2)*xNumberOfNodes] = -aSouth
            A[(yNumberOfNodes-1)*xNumberOfNodes+j][(yNumberOfNodes-1)*xNumberOfNodes+j-1] = -aWest
            A[(yNumberOfNodes-1)*xNumberOfNodes+j][(yNumberOfNodes-1)*xNumberOfNodes+j] = ap
            A[(yNumberOfNodes-1)*xNumberOfNodes+j][(yNumberOfNodes-1)*xNumberOfNodes+j+1] = -aEast
            b[(yNumberOfNodes-1)*xNumberOfNodes+j] = aNorth*Tinf
    #############################
    # Solving the linear system
    ############################
    solution = np.linalg.solve(np.array(A),np.array(b))
    temperatureField = []
    for i in range(yNumberOfNodes):
        temperatureField.append([])
        for j in range(xNumberOfNodes):
            temperatureField[i].append(solution[i*xNumberOfNodes+j])
        
    #############################################
    #Mirror of the temperature field and the mesh
    #############################################
    auxTemp = []
    auxMesh = []
    for i in range(1,yNumberOfNodes):
        auxTemp.append(temperatureField[i])
        auxMesh.append(-yNodesPositions[i])
    auxTemp = np.flip(auxTemp,0)
    auxMesh = np.flip(auxMesh,0)  
    for i in range(len(auxTemp)):
        temperatureField.insert(i,auxTemp[i])
        yNodesPositions.insert(i,auxMesh[i])
    temperatureField = np.array(temperatureField)
    yNodesPositions = np.array(yNodesPositions)
    return temperatureField, xNodesPositions, yNodesPositions
    
##################################
#Plotting the solution
########################
temp,x,y = getTemperatureField(4,5)
xx, yy = np.meshgrid(x,y)
plt.contourf(xx,yy,np.array(temp))
plt.colorbar(orientation="vertical")
plt.xlabel("x[m]")

#########################
##Refinando a malha
#########################
#refineFactor = 6
#xNumbNodes = 4
#yNumbNodes = 5
#oldTemp, x, y = getTemperatureField(xNumbNodes,yNumbNodes)
#diff = 100
#count = 0
#temp,x,y = getTemperatureField(refineFactor*xNumbNodes,refineFactor*yNumbNodes)
#aux = temp[0][refineFactor*xNumbNodes-1]-oldTemp[0,xNumbNodes-1]
#soma = aux*aux
#index = int(0.5*(1+yNumbNodes)-1)
#index2 = int(0.5*(1+refineFactor*yNumbNodes)-1)
#aux = temp[index2][refineFactor*xNumbNodes-1] - oldTemp[index][xNumbNodes-1]
#soma += aux*aux
#diff = mt.sqrt(soma/2)
#oldTemp = temp
#count += 1
#refineFactor += 1
#print(count)
#print(diff)

#############################
#Heat transfer
############################


    

