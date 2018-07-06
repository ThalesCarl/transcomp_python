#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  6 11:53:08 2018

@author: thales
"""
import numpy as np
import auxiliar as aux

u = np.loadtxt('./results/uFields.csv',dtype=float,delimiter=',',skiprows=0)
v = np.loadtxt('./results/vFields.csv',dtype=float,delimiter=',',skiprows=0)
yNodes = 70
xNodes = yNodes
L = 1.0

###############################
#Generating the mesh to plot
###############################
xNodesPositions = []
deltaX = L/(float(xNodes));
xSum = 0.5*deltaX;
for i in range(xNodes):
    xNodesPositions.append(xSum)
    xSum += deltaX
yNodesPositions = xNodesPositions
xNodesPositions = np.array(xNodesPositions)
yNodesPositions = np.array(yNodesPositions)
xx, yy = np.meshgrid(xNodesPositions,yNodesPositions)


#######################################
#Plotting the field
#######################################
aux.plotTheField(u,v,xNodes,yNodes,"./results/streamline70x70.png",xx,yy)