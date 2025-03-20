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
Ea = [220000,485000,260000,351900,532000,470000]
# stress exponent
n = [2.6,4.7,3.4,4.97,3.5,4]
# pre-exponential viscous factor (Pa s)
yita_0 = [1.99e17,1.98e27,1.26e24,1.08e34,3.98e16,5.01e20]
# activation volume
Va = [0,0,0,0,0,0]

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
        plt.plot(np.log10(yita),d/1000,label='Greenstone sequences',color='red', linewidth=2)
    if j == 1:
        plt.plot(np.log10(yita),d/1000,label='Dry diabase',color='yellowgreen')
    if j == 2:
        plt.plot(np.log10(yita),d/1000,label='Wet diabase1',color='forestgreen')
    if j == 3:
        plt.plot(np.log10(yita),d/1000,label='Wet diabase2',color='lightgreen')        
    if j == 4:
        plt.plot(np.log10(yita),d/1000,label='Dry olivine',color='mediumorchid')
    if j == 5:
        plt.plot(np.log10(yita),d/1000,label='Wet olivine',color='blueviolet')      

yita = [1e18]*len(d)
plt.plot(np.log10(yita),d/1000,label='Serpentinite',color='black')    

'''Calculate the plastic part'''
#Cohesion
#C = [43.91,49.73,31.99]
C = [10]
#Frictional angle
#fi = [28.54,33.64,13]
fi = [8.65]

yita = np.zeros(len(d))
for i in range(len (yita)):
    yita[i] = (C[0]+P[i]*np.sin(np.deg2rad(fi[0]))*1)/2/eps
    if yita[i]<1e18:
        yita[i]=1e18
plt.plot(np.log10(yita),d/1000,label='Plastic criteria',linestyle='--',color='black')
     
'''
yita = np.zeros(len(d))
for i in range(len (yita)):
    yita[i] = (C[0]+P[i]*np.sin(np.deg2rad(fi[0]))*0.05)/2/eps
    if yita[i]<1e18:
        yita[i]=1e18

plt.plot(np.log10(yita),d/1000,label='Frictional criterion-wet',linestyle='--',color='grey')'''

#plt.text(21, 2, 'When σ>σ_yield', fontsize=10, color='black')
#plt.annotate('', xy = (21.5, 3), xytext = (26, 3), arrowprops = dict(arrowstyle = '->', color = 'black'))
plt.xticks([18, 20, 22, 25, 30, 35],fontsize=12)
plt.yticks(fontsize=12)
plt.gca().invert_yaxis()
plt.xlabel('Viscosity (log Pa·s)',fontsize=12)
plt.ylabel('Depth (km)',fontsize=12)
plt.legend(loc='best',fontsize=12)
plt.savefig('Greenstone Rheology.png', dpi=1200)
    