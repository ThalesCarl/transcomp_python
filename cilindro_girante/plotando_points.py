#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  8 09:30:51 2018

@author: thales
"""

import numpy as np
parado = 0.000479969940433

cds = np.loadtxt("./results/maximum_error_cds_points.csv",dtype=float,delimiter=',')
uds = np.loadtxt("./results/maximum_error_uds_points.csv",dtype=float,delimiter=',')
exp = np.loadtxt("./results/maximum_error_exp_points.csv",dtype=float,delimiter=',')
points = np.loadtxt("./results/maximum_error_points.csv",dtype=float,delimiter=',')

import matplotlib.pyplot as plt


plt.plot(points,cds,label = "CDS")
plt.plot(points,uds,label = "UDS")
plt.plot(points,exp,label = "exponecial")
plt.xlabel('Número de pontos')
plt.ylabel("Erro")

legend = plt.legend(loc='upper right')
plt.savefig("./results/plot_points.png")