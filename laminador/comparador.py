#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 25 13:06:29 2018

@author: thales
"""

firstFile = "./results/aleta.csv"
secondFile = "./results/temperature_field_uds.csv"

#import csv
#
#aleta = []
#
#with open(firstFile,'r') as f1:
#    reader = csv.reader(f1)
#    i = 0
#    for row in reader:
#        aleta.append([])
#        for value in row:
#            aleta[i].append(float(value))
#        i=i+1

import numpy as np
aleta = np.loadtxt(firstFile,dtype=float,delimiter=',')
laminador = np.loadtxt(secondFile,dtype=float,delimiter=',')

difference = aleta - laminador

import csv

with open("./results/difference_uds.csv","w") as output:
    writer = csv.writer(output,lineterminator='\n')
    writer.writerows(difference)
