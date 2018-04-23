#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 23 16:01:21 2018

@author: tclavoratti
"""

##############################
#  Setting input paramethers
##############################

rho = 2700.0
cp = 900.0
k = 230
h = 20000
r0 = 0.05
wallLength = 0.005
numberOfNodes = 5;

#############################
# Mesh generation
#############################

nodePositions = []
surfacePositions = []
accumulator = r0
delta = wallLength/(numberOfNodes - 1)
for i in range(numberOfNodes):
    nodePositions.append(accumulator)
    accumulator += delta
accumulator = r0
surfacePositions.append(r0)
accumulator += 0.5*delta;
for i in range(numberOfNodes - 1):
    surfacePositions.append(accumulator)
    accumulator += delta
accumulator -= 0.5*delta;
surfacePositions.append(accumulator)

###########################
# Matrix coefficients
###########################

A=[]
A.append([])
A[0].append(1)
for i in range(numberOfNodes - 1):
    A[0].append(0)
    
for i in range(1,numberOfNodes - 1):
    