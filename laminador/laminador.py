"""
Aluno: Thales Carl Lavoratti (151000656)
Código do problema bidimensional de um laminador
"""


import numpy as np
import matplotlib.pyplot as plt


##############################
#  Setting input paramethers
##############################
yN = 5
xN = 5
w = 1.0 #[m]
L = 1.0 #[m]
t = 0.01 #[m]
k = 55.0#[W/mºC]
h = 1200 #[W/m^2ºC]
T0 = 1250.0 #[ºC]
Tinf = 25.0 #[ºC]
cp = 523.0 #[J/kgK]
u = 0.0 #[m/s]
rho = 7843.0 #[kg/m^3]
cds = 0#boolean que escolhe se usa método cds ou uds
yNumberOfNodes = int(0.5*(1+yN))
xNumberOfNodes = xN
gamma = k/cp
M = rho*u

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
        De = gamma*0.5*deltaY/deltaX
        Dw = gamma*0.5*deltaY/(0.5*deltaX)
        Dn = gamma*deltaX/deltaY
        Me = M*0.5*deltaY
        Mw = M*0.5*deltaY
        if(cds):
            aEast  =  De - 0.5*Me
            aWest  =  Dw + Mw
            aSouth =  0.0
            aNorth =  Dn
        else:
            aEast = De
            aWest = Dw + Mw
            aNorth = Dn
            aSouth = 0.0            
        ap = aEast + aNorth + aWest + aSouth
        A[j][j] = ap
        A[j][j+1] = -aEast
        A[j][j+xNumberOfNodes] = -aNorth
        b[j] = aWest * T0        
    elif j == (xNumberOfNodes - 1):
        Dw = gamma*0.5*deltaY/deltaX
        Dn = gamma*0.5*deltaX/deltaY
        Me = M*0.5*deltaY
        Mw = M*0.5*deltaY
        if(cds):
            aEast  =  0.0
            aWest  =  Dw + 0.5*Mw
            aSouth =  0.0
            aNorth =  Dn
        else:
            aEast = 0.0
            aWest = Dw + Mw
            aNorth = Dn
            aSouth = 0.0        
        ap = aEast + aNorth + aWest + aSouth
        A[j][j-1] = -aWest
        A[j][j] = ap
        A[j][j+1] = -aEast        
        A[j][j+xNumberOfNodes] = -aNorth
        ap = aEast + aNorth + aWest + aSouth
        A[j][j-1] = -aWest
        A[j][j] = ap
        A[j][j+xNumberOfNodes] = -aNorth
    else:
        De = gamma*0.5*deltaY/deltaX
        Dw = gamma*0.5*deltaY/deltaX
        Dn = gamma*deltaX/deltaY
        Me = M*0.5*deltaY
        Mw = M*0.5*deltaY
        if(cds):
            aEast  =  De - 0.5*Me
            aWest  =  Dw + 0.5*Mw
            aSouth =  0.0
            aNorth =  Dn
        else:
            aEast = De
            aWest = Dw + Mw
            aNorth = Dn
            aSouth = 0.0        
        ap = aEast + aNorth + aWest + aSouth
        A[j][j-1] = -aWest
        A[j][j] = ap
        A[j][j+1] = -aEast        
        A[j][j+xNumberOfNodes] = -aNorth
        
        
#center
for i in range(1,yNumberOfNodes-1):
    for j in range(xNumberOfNodes):
        if  j == 0:
            De = gamma*deltaY/deltaX
            Dw = gamma*deltaY/(0.5*deltaX)
            Dn = gamma*deltaX/deltaY
            Ds = gamma*deltaX/deltaY
            Me = M*deltaY
            Mw = M*deltaY
            if(cds):
                aEast  =  De - 0.5*Me
                aWest  =  Dw + Mw
                aSouth =  Ds
                aNorth =  Dn
            else:
                aEast = De
                aWest = Dw + Mw
                aNorth = Dn
                aSouth = Ds        
            ap = aEast + aNorth + aWest + aSouth
            A[i*xNumberOfNodes+j][j+xNumberOfNodes*(i-1)] = -aSouth
            A[i*xNumberOfNodes+j][i*xNumberOfNodes+j] = ap
            A[i*xNumberOfNodes+j][i*xNumberOfNodes+j+1] = -aEast
            A[i*xNumberOfNodes+j][j+xNumberOfNodes*(i+1)] = -aNorth
            b[i*xNumberOfNodes+j] = aWest * T0
        elif j == (xNumberOfNodes - 1):
            Dw = gamma*deltaY/deltaX
            Dn = gamma*0.5*deltaX/deltaY
            Ds = gamma*0.5*deltaX/deltaY
            Me = M*deltaY
            Mw = M*deltaY
            if(cds):
                aEast  =  0.0
                aWest  =  Dw + 0.5*Mw
                aSouth =  Ds
                aNorth =  Dn
            else:
                aEast = 0.0
                aWest = Dw + Mw
                aNorth = Dn
                aSouth = Ds            
            ap = aEast + aWest + aSouth + aNorth
            A[i*xNumberOfNodes+j][j+xNumberOfNodes*(i-1)] = -aSouth
            A[i*xNumberOfNodes+j][i*xNumberOfNodes+j-1] = -aWest
            A[i*xNumberOfNodes+j][i*xNumberOfNodes+j] = ap
            A[i*xNumberOfNodes+j][j+xNumberOfNodes*(i+1)] = -aNorth
        else:
            De = gamma*deltaY/deltaX
            Dw = gamma*deltaY/deltaX
            Dn = gamma*deltaX/deltaY
            Ds = gamma*deltaX/deltaY
            Me = M*deltaY
            Mw = M*deltaY
            if(cds):
                aEast  =  De - 0.5*Me
                aWest  =  Dw + 0.5*Mw
                aSouth =  Ds
                aNorth =  Dn
            else:
                aEast = De
                aWest = Dw + Mw
                aNorth = Dn
                aSouth = Ds            
            ap = aEast + aWest + aSouth + aNorth
            A[i*xNumberOfNodes+j][j+xNumberOfNodes*(i-1)] = -aSouth
            A[i*xNumberOfNodes+j][i*xNumberOfNodes+j-1] = -aWest
            A[i*xNumberOfNodes+j][i*xNumberOfNodes+j] = ap
            A[i*xNumberOfNodes+j][i*xNumberOfNodes+j+1] = -aEast
            A[i*xNumberOfNodes+j][j+xNumberOfNodes*(i+1)] = -aNorth
#top
for j in range(xNumberOfNodes):    
    if  j == 0:
        De = gamma*0.5*deltaY/deltaX
        Dw = gamma*0.5*deltaY/(0.5*deltaX)
        Ds = gamma*deltaX/deltaY
        Me = M*0.5*deltaY
        Mw = M*0.5*deltaY
        if(cds):
            aEast  =  De - 0.5*Me
            aWest  =  Dw + Mw
            aSouth =  Ds
            aNorth =  h*deltaX/cp
        else:
            aEast = De
            aWest = Dw + Mw
            aNorth = h*deltaX/cp
            aSouth = Ds
        ap = aEast + aNorth + aWest + aSouth
        A[(yNumberOfNodes-1)*xNumberOfNodes+j][j+(yNumberOfNodes-2)*xNumberOfNodes] = -aSouth
        A[(yNumberOfNodes-1)*xNumberOfNodes+j][(yNumberOfNodes-1)*xNumberOfNodes+j] = ap
        A[(yNumberOfNodes-1)*xNumberOfNodes+j][(yNumberOfNodes-1)*xNumberOfNodes+j+1] = -aEast
        b[(yNumberOfNodes-1)*xNumberOfNodes+j] = aWest * T0 + aNorth*Tinf
    elif j == (xNumberOfNodes - 1):
        Dw = gamma*0.5*deltaY/deltaX
        Ds = gamma*0.5*deltaX/deltaY
        Me = M*0.5*deltaY
        Mw = M*0.5*deltaY
        if(cds):
            aEast  =  0.0
            aWest  =  Dw + 0.5*Mw
            aSouth =  Ds
            aNorth =  h*0.5*deltaX/cp
        else:
            aEast = 0.0
            aWest = Dw + Mw
            aNorth = h*0.5*deltaX/cp
            aSouth = Ds
        ap = aEast + aNorth + aWest + aSouth          
        A[(yNumberOfNodes-1)*xNumberOfNodes+j][j+(yNumberOfNodes-2)*xNumberOfNodes] = -aSouth
        A[(yNumberOfNodes-1)*xNumberOfNodes+j][(yNumberOfNodes-1)*xNumberOfNodes+j-1] = -aWest
        A[(yNumberOfNodes-1)*xNumberOfNodes+j][(yNumberOfNodes-1)*xNumberOfNodes+j] = ap
        b[(yNumberOfNodes-1)*xNumberOfNodes+j] = aNorth*Tinf        
    else:
        De = gamma*0.5*deltaY/deltaX
        Dw = gamma*0.5*deltaY/deltaX
        Ds = gamma*deltaX/deltaY
        Me = M*0.5*deltaY
        Mw = M*0.5*deltaY
        if(cds):
            aEast  =  De - 0.5*Me
            aWest  =  Dw + 0.5*Mw
            aSouth =  Ds
            aNorth =  h*deltaX/cp
        else:
            aEast = De
            aWest = Dw + Mw
            aNorth = h*deltaX/cp
            aSouth = Ds
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
    
    
##################################
#Plotting the solution
########################
xx, yy = np.meshgrid(xNodesPositions,yNodesPositions)
plt.contourf(xx,yy,np.array(temperatureField))
plt.colorbar(orientation="vertical")
plt.xlabel("x[m]")

###############################
#Exporting the solution to csv
###############################

import csv

with open("./results/temperature_field_cds.csv","w") as output:
    writer = csv.writer(output,lineterminator='\n')
    writer.writerows(temperatureField)


