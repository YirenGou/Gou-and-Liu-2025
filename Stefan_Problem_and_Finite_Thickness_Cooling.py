# -*- coding: utf-8 -*-
"""
The Stefan Problem described in Turcotte and Schubert

Created on Mon Apr  3 11:04:01 2023

@author: Yiren Gou Email: yg7km@umsystem.edu
"""

import numpy as np
import math
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 12})
# Set font
font = {'family': 'Arial'}
plt.figure(figsize=(80/25.4, 80/25.4*3/4))

# Define parameters
L = 418 # latent heat of fusion kJ/kg
T0 = 0 # temperature of the upper surface K
Tm = T0+1200 # The solidification temperature for diabase K
t = np.linspace(0,280000000000,1001) # Time s
alpha = 1 # thermal diffusivity mm^2/s
c = 1040 # heat capacity J/(kg*K)
const = L*1000*np.sqrt(np.pi)/c/(Tm-T0) # value of the left side
lamda1 = np.linspace(0.2,1.6,1000)
right = []
for i in range(len(lamda1)):
    right.append(np.exp(-lamda1[i]**2)/lamda1[i]/math.erf(lamda1[i]))
for i in range(len(right)):
    if (right[i]-const)<0:
        lamda = lamda1[i]
        break   
ym=[]
for i in range(len(t)):
    ym.append((2*lamda*np.sqrt(t[i]*alpha*0.0001))/1000)


h=np.linspace(0,10,1000)
T=np.linspace(1473,1473,1000)
T[0]=273
T = T-273
plt.plot(T,h,label='time = '+'0'+' kyr',linestyle='dashed')


for i in range(100,len(t),200):
    T=[0]
    if ym[i]<10:
        h=np.linspace(0,ym[i],1000)
    else:
        break
    for j in range(1,len(h)):
        a = h[j]/ym[i]
        for k in range(1,10):
            a = a+2/np.pi/k*np.exp(-alpha*0.0001*k**2*np.pi**2*t[i]/(ym[i]*1000)**2)*np.sin(k*np.pi*h[j]/ym[i])
        T.append(T0+(Tm-T0)*a)
    T.append(Tm)
    h=np.append(h,10)
    plt.plot(T,h,label='time = '+str(round(t[i]/60/60/24/365/1000,1))+' kyr')

plt.legend(loc='lower left', fontsize=12,prop = {'size':7.4})
plt.show
ax = plt.gca()
ax.invert_yaxis()
ax.set_xlabel('Temperature (\u00B0C)',fontdict=font,fontsize=12)
ax.set_ylabel('Depth (km)',fontdict=font,fontsize=12)
plt.yticks([0, 2, 4, 6, 8, 10])
plt.xticks(fontname='Arial',fontsize=10)
plt.yticks(fontname='Arial',fontsize=10)
plt.savefig("Stefan and Plate Cooling Model.png",dpi=1200,bbox_inches = 'tight')
