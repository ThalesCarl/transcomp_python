#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  4 15:34:11 2018

@author: thales
"""

import numpy as np

#def gen(xNodes,yNodes,L,topWallVelocity,deltaT,rho,mi,u0,uStar,v0,vStar,method):
uStar = [[0,2,-3,0],[0,6,-7,0],[0,10,11,0]]
vStar = [[0,0,0],[9,-2,1],[-7,4,-5],[0,0,0]]
u0 = [[0,0,0,0],[0,0,0,0],[0,1,1,0]]
v0 = [[0,0,0],[0,0,0],[0,0,0],[0,0,0]]

yNodes = 3
xNodes = 3
L = 1.0 
topWallVelocity = 1.0 #[m/s]
mi = 1.0 #[Pa * s]
rho = 1.0 #[kg/m^3]
method = 1 # 0 = CDS; 1 = UDS
deltaT = 1

uCounter = (xNodes+1)*yNodes
vCounter = xNodes*(yNodes+1)
pCounter = xNodes*yNodes
deltaX = L/xNodes
deltaY = L/yNodes
ap0 = rho*deltaX*deltaY/deltaT
"""
A matriz A  foi dividida em nove quadrantes para facilitar o código
todos os quadrantes são gerados separadamente e depois concatenados.
O vetor b foi dividido em trÊs subvetores para facilitar o código e que
em seguida foram concatenados
"""
A00 = np.zeros((uCounter,uCounter))
A01 = np.zeros((uCounter,vCounter))
A02 = np.zeros((uCounter,pCounter))
A10 = np.zeros((vCounter,uCounter))
A11 = np.zeros((vCounter,vCounter))
A12 = np.zeros((vCounter,pCounter))
A20 = np.zeros((pCounter,uCounter))
A21 = np.zeros((pCounter,vCounter))
A22 = np.zeros((pCounter,pCounter))

b0 = np.zeros(uCounter)
b1 = np.zeros(vCounter)
b2 = np.zeros(pCounter)

################################
#primeito quadrante (matriz A00)
#################################
#bottom
for j in range(xNodes+1):
    if j == 0 or j == xNodes:
        ap = 1.0
        A00[j,j] = ap
        b0[j] = 0.0
    else:
        De = mi*deltaY/deltaX
        Dw = mi*deltaY/deltaX
        Dn = mi*deltaX/deltaY
        Ds = mi*deltaX/(0.5*deltaY)
        Mw = rho*deltaY*0.5*(uStar[0][j-1]+uStar[0][j])
        Me = rho*deltaY*0.5*(uStar[0][j]+uStar[0][j+1])
        Ms = rho*deltaX*0.5*(vStar[0][j-1]+vStar[0][j]) # nulo
        Mn = rho*deltaX*0.5*(vStar[1][j-1]+vStar[1][j])
        if method == 0:
            aEast  = De - 0.5*Me
            aWest  = Dw + 0.5*Mw
            aNorth = Dn - 0.5*Mn
            aSouth = Ds + 0.5*Ms
        elif method == 1:
            aEast  = De + max(0,-Me)
            aWest  = Dw + max(0, Mw)
            aNorth = Dn + max(0,-Mn)
            aSouth = Ds + max(0, Ms)
        ap = aEast +aWest + aSouth +aNorth + ap0 + Me - Mw + Mn - Ms
        A00[j][j-1] = -aWest
        A00[j][j] = ap
        A00[j][j+1] = -aEast        
        A00[j][j+xNodes+1] = -aNorth
        b0[j] = ap0 * u0[0][j]

#center
for i in range(1,yNodes-1):
    for j in range(xNodes+1):
        if j == 0 or j == xNodes:
            ap = 1.0
            A00[i*(xNodes+1)+j,i*(xNodes+1)+j] = ap
            b0[i*(xNodes+1)+j] = 0.0 
        else:
            De = mi*deltaY/deltaX
            Dw = mi*deltaY/deltaX
            Dn = mi*deltaX/deltaY
            Ds = mi*deltaX/deltaY
            Mw = rho*deltaY*0.5*(uStar[i][j-1]+uStar[i][j])
            Me = rho*deltaY*0.5*(uStar[i][j]+uStar[i][j+1])
            Ms = rho*deltaX*0.5*(vStar[i][j-1]+vStar[i][j]) 
            Mn = rho*deltaX*0.5*(vStar[i+1][j-1]+vStar[i+1][j])
            if method == 0:
                aEast  = De - 0.5*Me
                aWest  = Dw + 0.5*Mw
                aNorth = Dn - 0.5*Mn
                aSouth = Ds + 0.5*Ms
            elif method == 1:
                aEast  = De + max(0,-Me)
                aWest  = Dw + max(0, Mw)
                aNorth = Dn + max(0,-Mn)
                aSouth = Ds + max(0, Ms)
            ap = aEast +aWest + aSouth +aNorth + ap0 + Me - Mw + Mn - Ms
            A00[i*(xNodes+1)+j,(xNodes+1)*(i-1)+j] = -aSouth
            A00[i*(xNodes+1)+j,i*(xNodes+1)+j-1] = -aWest
            A00[i*(xNodes+1)+j,i*(xNodes+1)+j] = ap
            A00[i*(xNodes+1)+j,i*(xNodes+1)+j+1] = -aEast        
            A00[i*(xNodes+1)+j,(xNodes+1)*(i+1)+j] = -aNorth
            b0[i*(xNodes+1)+j] = ap0 *u0[i][j]
            
#top
for j in range(xNodes+1):
    i=yNodes-1
    if j == 0 or j == xNodes:
        ap = 1.0
        A00[(yNodes-1)*(xNodes+1)+j][(yNodes-1)*(xNodes+1)+j] = ap
        b0[(yNodes-1)*(xNodes+1)+j] = 0.0
    else:
        De = mi*deltaY/deltaX
        Dw = mi*deltaY/deltaX
        Dn = mi*deltaX/(deltaY*0.5)
        Ds = mi*deltaX/deltaY
        Mw = rho*deltaY*0.5*(uStar[i][j-1]+uStar[i][j])
        Me = rho*deltaY*0.5*(uStar[i][j]+uStar[i][j+1])
        Ms = rho*deltaX*0.5*(vStar[i][j-1]+vStar[i][j]) 
        Mn = rho*deltaX*0.5*(vStar[i+1][j-1]+vStar[i+1][j])
        if method == 0:
            aEast  = De - 0.5*Me
            aWest  = Dw + 0.5*Mw
            aNorth = Dn - 0.5*Mn
            aSouth = Ds + 0.5*Ms
        elif method == 1:
            aEast  = De + max(0,-Me)
            aWest  = Dw + max(0, Mw)
            aNorth = Dn + max(0,-Mn)
            aSouth = Ds + max(0, Ms)
        ap = aEast +aWest + aSouth +aNorth + ap0 + Me - Mw + Mn - Ms
        A00[(yNodes-1)*(xNodes+1)+j,(yNodes-2)*(xNodes+1)+j] = -aSouth
        A00[(yNodes-1)*(xNodes+1)+j,(yNodes-1)*(xNodes+1)+j-1] = -aWest
        A00[(yNodes-1)*(xNodes+1)+j,(yNodes-1)*(xNodes+1)+j] = ap
        A00[(yNodes-1)*(xNodes+1)+j,(yNodes-1)*(xNodes+1)+j+1] = -aEast        
        b0[(yNodes-1)*(xNodes+1)+j] = ap0 *u0[i][j] + aNorth * topWallVelocity
        

#segundo quadrante (matriz A01) é uma matriz nula


################################
#terceiro quadrante (matriz A02)
###############################
countRow = 0
countCol = 0
while(countRow < uCounter and countCol<pCounter):
    countRow += 1
    for j in range(xNodes-1):
        A02[countRow,countCol]=-deltaY
        countCol += 1
        A02[countRow,countCol]=deltaX
        countRow += 1
    countCol += 1
    countRow += 1
    

#quarto quadrante (matriz A10) é uma matriz nula


#####################################
#quinto quadrante (matriz A11)
#####################################
for i in range(yNodes+1):
    for j in range(xNodes):
        if i==0 or i==yNodes:
            ap = 1.0
            A11[i*(xNodes)+j,i*(xNodes)+j] = ap
            b1[i*(xNodes)+j] = 0.0
        else:
            Mw = rho*deltaY*0.5*(uStar[i-1][j]+uStar[i][j])
            Me = rho*deltaY*0.5*(uStar[i-1][j+1]+uStar[i][j+1])
            Mn = rho*deltaX*0.5*(vStar[i][j]+vStar[i+1][j])
            Ms = rho*deltaX*0.5*(vStar[i][j]+vStar[i-1][j])
            if j == 0:
                #borda esquerda
                De = mi*deltaY/deltaX
                Dw = mi*deltaY/(0.5*deltaX)
                Dn = mi*deltaX/deltaY
                Ds = mi*deltaX/deltaY
                if method == 0:
                    aEast  = De - 0.5*Me
                    aWest  = Dw + 0.5*Mw
                    aNorth = Dn - 0.5*Mn
                    aSouth = Ds + 0.5*Ms
                elif method == 1:
                    aEast  = De + max(0,-Me)
                    aWest  = Dw + max(0, Mw)
                    aNorth = Dn + max(0,-Mn)
                    aSouth = Ds + max(0, Ms)
                ap = aEast +aWest + aSouth +aNorth + ap0 + Me - Mw + Mn - Ms
                A11[i*xNodes+j,xNodes*(i-1)+j] = -aSouth
                A11[i*xNodes+j,i*xNodes+j] = ap
                A11[i*xNodes+j,i*xNodes+j+1] = -aEast
                A11[i*xNodes+j,xNodes*(i+1)+j] = -aNorth
                b1[i*xNodes+j] = ap0 *v0[i][j]
            elif j == (xNodes-1):
                #borda direita
                De = mi*deltaY/(0.5*deltaX)
                Dw = mi*deltaY/deltaX
                Dn = mi*deltaX/deltaY
                Ds = mi*deltaX/deltaY
                if method == 0:
                    aEast  = De - 0.5*Me
                    aWest  = Dw + 0.5*Mw
                    aNorth = Dn - 0.5*Mn
                    aSouth = Ds + 0.5*Ms
                elif method == 1:
                    aEast  = De + max(0,-Me)
                    aWest  = Dw + max(0, Mw)
                    aNorth = Dn + max(0,-Mn)
                    aSouth = Ds + max(0, Ms)
                ap = aEast +aWest + aSouth +aNorth + ap0 + Me - Mw + Mn - Ms
                A11[i*xNodes+j,xNodes*(i-1)+j] = -aSouth
                A11[i*xNodes+j,i*xNodes+j-1] = -aWest
                A11[i*xNodes+j,i*xNodes+j] = ap
                A11[i*xNodes+j,xNodes*(i+1)+j] = -aNorth
                b1[i*xNodes+j] = ap0 *v0[i][j]
            else:
                De = mi*deltaY/deltaX
                Dw = mi*deltaY/deltaX
                Dn = mi*deltaX/deltaY
                Ds = mi*deltaX/deltaY
                if method == 0:
                    aEast  = De - 0.5*Me
                    aWest  = Dw + 0.5*Mw
                    aNorth = Dn - 0.5*Mn
                    aSouth = Ds + 0.5*Ms
                elif method == 1:
                    aEast  = De + max(0,-Me)
                    aWest  = Dw + max(0, Mw)
                    aNorth = Dn + max(0,-Mn)
                    aSouth = Ds + max(0, Ms)
                ap = aEast +aWest + aSouth +aNorth + ap0 + Me - Mw + Mn - Ms
                A11[i*xNodes+j,xNodes*(i-1)+j] = -aSouth
                A11[i*xNodes+j,i*xNodes+j-1] = -aWest
                A11[i*xNodes+j,i*xNodes+j] = ap
                A11[i*xNodes+j,i*xNodes+j+1] = -aEast
                A11[i*xNodes+j,xNodes*(i+1)+j] = -aNorth
                b1[i*xNodes+j] = ap0 *v0[i][j]

##############################
#Sexto quadrante (matriz A12)
##############################
for i in range(xNodes,vCounter-xNodes):
    for j in range(pCounter):
        if j==i:
            A12[i][j]=deltaX
        elif j == (i-xNodes):
            A12[i][j] = -deltaX
        else: continue

##############################
#Sétimo quadrante (matriz A20)
##############################
countRow = 0
countCol = 0
while(countRow < pCounter and countCol<uCounter):
    for i in range(xNodes):
        A20[countRow,countCol]=-deltaY
        countCol += 1
        A20[countRow,countCol]=deltaX
        countRow += 1
    countCol += 1
A20[0,0] = 0.0
A20[0,1] = 0.0
    
##############################
#Oitavo quadrante (matriz A21)
##############################
for i in range(pCounter):
    for j in range(vCounter):
        if j==i:
            A21[i][j]=-deltaX
        elif j == (i+xNodes):
            A21[i][j] = deltaX
        else: continue    
A21[0,0] = 0.0
A21[0,xNodes]= 0.0
#Nono quadrante (matriz A10) é uma matriz nula,exceto pelo primeiro elemento
A22[0,0] = 1.0
#################################
#Concatenando a matriz e o vetor
#################################
A0 = np.hstack((A00,A01,A02))
A1 = np.hstack((A10,A11,A12))
A2 = np.hstack((A20,A21,A22))
A = np.vstack((A0,A1,A2))

b = np.hstack((b0,b1,b2))

#return A, b

#def fieldBuilder(linearVector,firstDirection,secondDirection):
#    field, remainder = np.split(linearVector,[secondDirection])        
#    for i in range(firstDirection-1):
#        builtInVector,remainder= np.split(remainder,[secondDirection])
#        field = np.vstack((field,builtInVector))
#    return field