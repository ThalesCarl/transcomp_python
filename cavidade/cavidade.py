"""
Aluno: Thales Carl Lavoratti (15100656)
Código do problema da cavidade
"""

import math as mt
import numpy as np
import matplotlib.pyplot as plt
import auxiliar as aux

##############################
#  Setting input paramethers
##############################
yNodes = 3
xNodes = 3
L = 1.0 
topWallVelocity = 1.0 #[m/s]
mi = 0.001 #[Pa * s]
rho = 1.0 #[kg/m^3]
method = 1 # 0 = CDS; 1 = UDS

deltaT = 0.1 #[s]
numberOfTimeSteps = 10
maxCounter = 20
tolerance = 0.01

#####################################
#Setting the initial velocity Fields
####################################
initial_U_Field = np.zeros((xNodes + 1,yNodes))
initial_V_Field = np.zeros((xNodes,yNodes+1))

uFields = []
uFields.append(initial_U_Field)
vFields = []
vFields.append(initial_V_Field)
uCounter = (xNodes+1)*yNodes
vCounter = xNodes*(yNodes+1)

###################
#Appling the method
###################
u0 = initial_U_Field
v0 = initial_V_Field
for i in range(numberOfTimeSteps):
    uStar = u0
    vStar = v0
    difference = 1.0
    counter = 0
    while(difference > tolerance and counter < maxCounter):
        A, b = aux.gen(xNodes,yNodes,L,topWallVelocity,deltaT,rho,mi,u0,uStar,v0,vStar,method)
        #solution u and v
        solution = np.linalg.solve(A,b)
        u,v,P = np.split(solution,[uCounter+1,uCounter+vCounter+1])
        uStar = u
        vStar = v
        difference = 0.001
        ++counter
    uFields.append(u)
    vFields.append(v)        
    u0 = u
    v0 = v