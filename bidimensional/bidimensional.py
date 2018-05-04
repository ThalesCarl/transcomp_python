#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Aluno: Thales Carl Lavoratti (151000656)
Código do problema bidimensional de uma aleta
"""

import math as mt
import numpy as np
import matplotlib.pyplot as plt

##############################
#  Setting input paramethers
##############################

w = 1.0 #[m]
L = 1.0 #[m]
t = 0.1 #[m]
h = 26.3 #[W/m^2]
T0 = 100.0 #[ºC]
Tinf = 20.0 #[ºC]
yNumberOfNodes = 4
xNumberOfNodes = 8

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

