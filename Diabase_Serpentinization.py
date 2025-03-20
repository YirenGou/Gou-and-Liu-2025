# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 14:23:01 2023

@author: Yiren Gou Email: yg7km@umsystem.edu
"""

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.family'] = 'Arial'

'''Global parameters'''
# The geothermal gradient is (k in m)
gg = (619-298)/10000
# The density of the crust is 3000 kg/m3
rho = 3000
# Create the linespace for the depth
d = np.linspace(0,30000,1001)
# Suppose a strain rate (s-1)
eps = 1e-14 
# Gas constant (J/mol/K)
R = 8.314
# Gravitational acceleration (m/s2)
g = 9.8 


'''Rehology parameters'''
# activation energy (J/mol)
Ea = [220000,66000]
# stress exponent
n = [2.6,2.8]
# pre-exponential viscous factor (Pa s)
yita_0 = [1.99e17,1e23]
# activation volume
Va = [0,0]

'''Calculate the viscosity profile'''
# Generate pressure space (Pa)
P = rho*g*d
# Generate temperature space (k)
T = 298+gg*d

for j in range(len(Ea)):        
    yita = np.zeros(len(d))
    for i in range(len (yita)):
        yita[i] = eps**((1-n[j])/n[j])*yita_0[j]**(1/n[j])*np.exp((Ea[j]+P[i]*Va[j])/(n[j]*R*T[i]))
    if j == 0:
        yita_g = yita
        plt.plot(np.log10(yita),d/1000,label='Greenstone sequences',color='red', linewidth=2)
    if j == 1:
        plt.plot(np.log10(yita),d/1000,label='Serpentinite',color='green')
        yita_s = yita

'''Calculate the mixture'''
x_g =[0.8, 0.6, 0.4, 0.2] # Greenstone rocks
x_s = [0.2, 0.4, 0.6, 0.8] # Serpentinite

for i in range(len(x_g)):
    yita = np.zeros(len(d))
    for j in range(len(yita)):
        yita[j] = (x_g[i]*(yita_g[j]**(1/3))+x_s[i]*(yita_s[j]**(1/3)))**(3)
    if i == 0:            
        plt.plot(np.log10(yita),d/1000,label='20 wt% serpentinite',color='lightcoral', linewidth=2)
    if i == 1:
        plt.plot(np.log10(yita),d/1000,label='40 wt% serpentinite',color='darkorange', linewidth=2)
    if i == 2:
        plt.plot(np.log10(yita),d/1000,label='60 wt% serpentinite',color='goldenrod', linewidth=2)
    if i == 3:
        plt.plot(np.log10(yita),d/1000,label='80 wt% serpentinite',color='olive', linewidth=2)


'''Calculate the plastic part'''
#Cohesion
#C = [43.91,49.73,31.99]
C = [10,0.1]
#Frictional angle
#fi = [28.54,33.64,13]
fi = [8.65,1]

yita = np.zeros(len(d))
for i in range(len (yita)):
    yita[i] = (C[0]*10**6+P[i]*np.sin(np.deg2rad(fi[0]))*1)/2/eps
    if yita[i]<1e18:
        yita[i]=1e18
plt.plot((np.log10(yita))[:435],(d/1000)[:435],label='Plastic criteria_greenstone',linestyle='--',color='red')
     
yita = np.zeros(len(d))
for i in range(len (yita)):
    yita[i] = (C[1]*10**6+P[i]*np.sin(np.deg2rad(fi[1]))*1)/2/eps
    if yita[i]<1e18:
        yita[i]=1e18
plt.plot((np.log10(yita))[:145],(d/1000)[:145],label='Plastic criteria_serpentinite',linestyle='--',color='green')
'''
yita = np.zeros(len(d))
for i in range(len (yita)):
    yita[i] = (C[0]+P[i]*np.sin(np.deg2rad(fi[0]))*0.05)/2/eps
    if yita[i]<1e18:
        yita[i]=1e18

plt.plot(np.log10(yita),d/1000,label='Frictional criterion-wet',linestyle='--',color='grey')'''

plt.xticks([18, 20, 22, 25, 30, 35],fontsize=12)
plt.yticks(fontsize=12)
plt.gca().invert_yaxis()
plt.xlabel('Viscosity (log Pa·s)',fontsize=12)
plt.ylabel('Depth (km)',fontsize=12)
plt.legend(loc='best',fontsize=12)
plt.savefig('Mixture.png', dpi=1200)
    