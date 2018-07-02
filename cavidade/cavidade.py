"""
Aluno: Thales Carl Lavoratti (15100656)
Código do problema da cavidade
"""

import numpy as np
import matplotlib.pyplot as plt
import auxiliar as aux

##############################
#  Setting input paramethers
##############################
yNodes = 30
xNodes = 30
L = 1.0 
topWallVelocity = 1.0 #[m/s]
mi = 0.001 #[Pa * s]
rho = 1.0 #[kg/m^3]
method = 1 # 0 = CDS; 1 = UDS

deltaT = 10 #[s]
numberOfTimeSteps = 1
maxCounter = 5
tolerance = 0.01

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
for i in range(numberOfTimeSteps):
    uStar = u0
    vStar = v0
    difference = 1.0
    counter = 0
    while(counter < maxCounter):
        A, b = aux.gen(xNodes,yNodes,L,topWallVelocity,deltaT,rho,mi,u0,uStar,v0,vStar,method)
        solution = np.linalg.solve(A,b)
        uLinear,vLinear,pLinear = np.split(solution,[uCounter,uCounter+vCounter])
        u = aux.fieldBuilder(uLinear,yNodes,xNodes+1)
        v = aux.fieldBuilder(vLinear,yNodes+1,xNodes)
        pressureField = aux.fieldBuilder(pLinear,yNodes,xNodes)        
        uStar = u
        vStar = v
        counter += 1
        difference = 0.001
        print(counter)
    
    uFields.append(u)
    vFields.append(v) 
    pFields.append(pressureField)       
    u0 = u
    v0 = v