#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 24 16:55:30 2018

@author: tclavoratti
"""
def solverTDMA(A,B):
    a = []
    b = []
    c = []
    d = []
    x = []
    
    for j in range(len(A)):
        x.append(0)
    a.append(A[0][0])
    b.append(A[0][1])
    c.append(float('nan'))
    d.append(B[0])
    for i in range(1,len(A)-1):
        a.append(A[i][i])
        b.append(A[i][i+1])
        c.append(A[i][i-1])
        d.append(B[i])
    a.append(A[len(A)-1][len(A)-1])
    b.append(float('nan'))
    c.append(A[len(A)-1][len(A)-2])
    d.append(B[len(B)-1])
    n = len(a)
    for i in range(1,n):
        aux = c[i]/a[i-1]
        a[i] = a[i]-aux * b[i-1]
        d[i] = d[i]-aux * d[i-1]
        c[i] = 0
    x[n-1] = d[n-1]/a[n-1]
    for i in range(n,1,-1):
        x[i-2] = (d[i-2]-b[i-2]*x[i-1])/a[i-2]
    print(n)
    
solverTDMA([[1,-1,0,0,0],[1,1,-1,0,0],[0,1,-1,1,0],[0,0,-1,1,1],[0,0,0,-1,2]],[0,1,2,-1,-2])