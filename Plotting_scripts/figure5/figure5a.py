#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 11:39:48 2021

@author: daan
"""

import pandas as pd

import numpy as np
import matplotlib.pyplot as plt
import math

#%%
LW=5
FS=25
# %%
data1 = pd.read_csv('bifurcation.txt',sep='\s+',header=None)
data1 = pd.DataFrame(data1)

# %%
AMOC_1=data1[4]*10**-6
CO2_1=data1[6]*10**6
C1_1=data1[7]
A1_1=data1[8]
P1_1=data1[9]
C2_1=data1[10]
A2_1=data1[11]
P2_1=data1[12]
C3_1=data1[13]
A3_1=data1[14]
P3_1=data1[15]
C5_1=data1[16]
A5_1=data1[17]
P5_1=data1[18]
C6_1=data1[19]
A6_1=data1[20]
P6_1=data1[21]
C7_1=data1[22]
A7_1=data1[23]
P7_1=data1[24]
TC_1=data1[25]

#%%
for i in range(19):
   if np.size(data1.loc[data1[i+6] < 0])>0:
       print(i+6)

#%%
CO2=CO2_1
Tadjust=np.zeros((np.size(CO2),1))
Tadjust5=np.zeros((np.size(CO2),1))
for j in range (np.size(CO2)):
    Tadjust[j]=2.085*0.54*5.35*np.log((CO2[j])/145) # Adjustment temperature CO2 concentration

T7=0.33+Tadjust   
x_2=np.where(T7<-1)

if np.size(x_2)>0:
    print('2')
    print(x_2)

eps_TEMP=1.5
eps1_base=0.9
eps2_base=1.25
eps7_base=0.62
eps5_base=0.35
    
eps1=eps_TEMP*Tadjust*(-0.1)+eps1_base
eps2=eps_TEMP*Tadjust*(-0.1)+eps2_base
eps5=eps_TEMP*Tadjust*(-0.1/3.25)+eps5_base
eps7=eps_TEMP*Tadjust*(-0.1)+eps7_base
    
if np.size(np.where(eps1<0))>0:
    print('eps1')
if np.size(np.where(eps2<0))>0:
    print('eps2')
if np.size(np.where(eps5<0))>0:
    print('eps5')
if np.size(np.where(eps7<0))>0:
    print('eps7')

#%%
x=data1.loc[data1[2] == 3]
xx=np.abs(x[1])
xHB=min(xx)-1
xST=11179
xP=12862
xSW=11872

#%%
AMOC_LC=np.array([AMOC_1[xHB],15,16.5,18,19.5,21,22.5,24,25.5,27,28.5,30])
max_LC=np.array([CO2_1[xHB],237,253,266,278,288,298,307,315,323,330,337])
min_LC=np.array([CO2_1[xHB],165,158,153,149,146,145,144,143,142,142,142])

AMOC_LC2=np.array([AMOC_1[xST],7.629,7.626,7.58,7.54,7.07,6.48,6.38,6.37])
max_LC2=np.array([CO2_1[xST],188,188,192,191,186,178,176,176])
min_LC2=np.array([CO2_1[xST],172.5,172,168,168,171,175,176,176])
#%%
x_A=[15,15]
y_A=[0,1000]

figure1 = plt.figure(1,figsize=(10,6.5))
plt.plot(x_A,y_A,'k',linewidth=LW-2,zorder=1,label='_nolegend_')
plt.plot(AMOC_1[0:xHB],CO2_1[0:xHB],color='tab:red', linestyle='dashed',zorder=1,linewidth=LW)
plt.plot(AMOC_1[xHB:xST],CO2_1[xHB:xST],color='tab:blue',zorder=1,linewidth=LW)
plt.scatter(AMOC_1[xHB],CO2_1[xHB],color='black',marker='s',s=250,zorder=2)
plt.scatter(AMOC_1[xST],CO2_1[xST],color='black',marker='s',s=100,zorder=2)
plt.plot(AMOC_LC,max_LC,color='black',zorder=1,linewidth=LW)
plt.plot(AMOC_LC,min_LC,color='black',zorder=1,linewidth=LW,label='_nolegend_')
plt.plot(AMOC_LC2,max_LC2,color='black',zorder=1,linewidth=LW,label='_nolegend_')
plt.plot(AMOC_LC2,min_LC2,color='black',zorder=1,linewidth=LW,label='_nolegend_')
plt.plot(AMOC_1[xST:xSW],CO2_1[xST:xSW],color='tab:red', linestyle='dashed',zorder=1,linewidth=LW,label='_nolegend_')
plt.plot(AMOC_1[xSW:xP],CO2_1[xSW:xP],color='tab:blue',zorder=1,linewidth=LW,label='_nolegend_')

plt.grid()
plt.title('Bifurcation diagram L-HB',fontsize=FS)
plt.ylabel('CO$_2$ [ppm]',fontsize=FS)
plt.xlabel('AMOC [Sv]',fontsize=FS)
plt.xticks(fontsize=FS-3)
plt.yticks(fontsize=FS-3)
plt.legend(['Unstable FP','Stable FP','Stable LC','HB'],loc="upper left",fontsize=FS-5,ncol=1)
plt.xlim(10,30)
plt.ylim(100,340)
plt.savefig('figure5a.png', format='png', dpi=300)