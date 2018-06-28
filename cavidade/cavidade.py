"""
Aluno: Thales Carl Lavoratti (15100656)
Código do problema da cavidade
"""

import math as mt
import numpy as np
import matplotlib.pyplot as plt

##############################
#  Setting input paramethers
##############################
yNodes = 3
xNodes = 3
L = 1.0 
w = 1.0 #[m]
topWallVelocity = 1.0 #[m/s]
mi = 0.001 #[Pa * s]
rho = 1.0 #[kg/m^3]

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

####################
##Appling the method
####################
#u0 = initial_U_Field
#v0 = initial_V_Field
#for i in range(numberOfTimeSteps):
#    uStar = u0
#    vStar = v0
#    difference = 1.0
#    counter = 0
#    while(difference > tolerance and counter < maxCounter):
#        A, b = aux.generator(xNodes,yNodes,L,w,topWallVelocity,deltaT,rho,mi,u0,uStar,v0,vStar)
#        #solution u and v
#        #bota a solução em algum tipo de matriz com indice i
#        uFields.append(u)
#        uStar = u
#        vStar = v
#        ++counter
#    u0 = u
#    v0 = v