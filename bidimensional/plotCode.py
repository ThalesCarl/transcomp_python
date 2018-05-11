#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 10 21:08:33 2018

@author: tclavoratti
"""

xx, yy = np.meshgrid(xNodesPositions,yNodesPositions)
plt.contourf(xx,yy,np.array(temperatureField))
plt.colorbar(orientation="horizontal",shrink=0.8,ticks = [100,90,80,70,60,50,40])