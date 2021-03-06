"""
Aluno: Thales Carl Lavoratti (15100656)
Código do problema da cavidade
"""

import numpy as np
import auxiliar as aux

##############################
#  Setting input paramethers
##############################
yNodes = 10
xNodes = yNodes
L = 1.0 
topWallVelocity = 1.0 #[m/s]
mi = 0.001 #[Pa * s]
rho = 1.0 #[kg/m^3]
method = 1 # 0 = CDS; 1 = UDS

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
    
###############################
#Generating the mesh to plot
###############################
xNodesPositions, yNodesPositions = aux.meshGenerator(L,xNodes,yNodes)
xx, yy = np.meshgrid(xNodesPositions,yNodesPositions)

#######################################
#Plotting the field
#######################################
aux.plotTheField(u,v,xNodes,yNodes,"./results/streamline20x20_cds.png",xx,yy)
aux.plotPressureField(pressureField,xNodes,yNodes,L,"./results/presure10x10_uds.png")

