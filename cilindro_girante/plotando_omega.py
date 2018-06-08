#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  8 09:30:51 2018

@author: thales
"""

import numpy as np
parado = 0.000479969940433

cds = np.loadtxt("./results/maximum_error_cds_w.csv",dtype=float,delimiter=',')
uds = np.loadtxt("./results/maximum_error_uds_w.csv",dtype=float,delimiter=',')
exp = np.loadtxt("./results/maximum_error_exp_w.csv",dtype=float,delimiter=',')
omega = np.loadtxt("./results/maximum_error_velocity.csv",dtype=float,delimiter=',')

import matplotlib.pyplot as plt

plt.plot(0.0,parado,'ko', label = "Parado")
plt.plot(omega,cds,label = "CDS")
plt.plot(omega,uds,label = "UDS")
plt.plot(omega,exp,label = "exponecial")
plt.xlabel('w [rad/s]')
plt.ylabel("Erro")
plt.axis([-0.1,1,0.0004,0.0022])
legend = plt.legend(loc='lower right')
#plt.savefig("./results_implicit/plot" + str (numberOfNodes) + "nodes.png")