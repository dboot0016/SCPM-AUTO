#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 30 17:59:21 2021

@author: daan
"""

import pandas as pd

import numpy as np
import matplotlib.pyplot as plt
import math

#%%
LW=5
FS=25
#%%
data2 = pd.read_csv('oscillation.txt',sep='\s+',header=None)
data2 = pd.DataFrame(data2)

#%%
per_x=1.83361E+11/(365*86400)

HB2=np.zeros((21,2001))

for i in range(2001):
    HB2[0:7,i]=data2.loc[3*i,:]
    HB2[7:14,i]=data2.loc[3*i+1,:]
    HB2[14:21,i]=data2.loc[3*i+2,:]

# %%
Vat=1.7604418363824646*10**20;
V1=2.71425*10**16*0.97
V2=9.047500*10**15*0.97
V3=2.442825*10**17*0.97
V4=5.699925*10**17*0.97
V5=4.523750*10**16*0.97
V6=5.428500*10**17*0.97
V7=9.047500*10**15*0.97

# %%

HB1=HB2

t11=HB1[0,:]*per_x
t_2=t11
ti_2=np.linspace(0,np.max(t_2),10000)
CO2_2=HB1[1,:]



C1_2=HB1[2,:]
C2_2=HB1[5,:]
C3_2=HB1[8,:]
C5_2=HB1[11,:]
C6_2=HB1[14,:]
C7_2=HB1[17,:]

A1_2=HB1[3,:]
A2_2=HB1[6,:]
A3_2=HB1[9,:]
A5_2=HB1[12,:]
A6_2=HB1[15,:]
A7_2=HB1[18,:]

P1_2=HB1[4,:]
P2_2=HB1[7,:]
P3_2=HB1[10,:]
P5_2=HB1[13,:]
P6_2=HB1[16,:]
P7_2=HB1[19,:]

TC_2=HB1[20,:]

# %%
CO2_2i1=np.interp(ti_2,t_2,CO2_2)

j=(np.where(CO2_2i1==min(CO2_2i1)))
jj=int(j[0][0])

CO2_2i=np.roll(CO2_2i1,-jj)


C1_2i=np.interp(ti_2,t_2,C1_2)
C2_2i=np.interp(ti_2,t_2,C2_2)
C3_2i=np.interp(ti_2,t_2,C3_2)
C5_2i=np.interp(ti_2,t_2,C5_2)
C6_2i=np.interp(ti_2,t_2,C6_2)
C7_2i=np.interp(ti_2,t_2,C7_2)

A1_2i=np.interp(ti_2,t_2,A1_2)
A2_2i=np.interp(ti_2,t_2,A2_2)
A3_2i=np.interp(ti_2,t_2,A3_2)
A5_2i=np.interp(ti_2,t_2,A5_2)
A6_2i=np.interp(ti_2,t_2,A6_2)
A7_2i=np.interp(ti_2,t_2,A7_2)

P1_2i=np.interp(ti_2,t_2,P1_2)
P2_2i=np.interp(ti_2,t_2,P2_2)
P3_2i=np.interp(ti_2,t_2,P3_2)
P5_2i=np.interp(ti_2,t_2,P5_2)
P6_2i=np.interp(ti_2,t_2,P6_2)
P7_2i=np.interp(ti_2,t_2,P7_2)

TC_2i=np.interp(ti_2,t_2,TC_2)
TA_2i=(3446475.31695093e+12-(3349332.54317694e+12-TC_2i*(10**22))*2)
# %%
TC=CO2_2i*Vat+C6_2i*V6+C2_2i*V2+C3_2i*V3+C5_2i*V5+C1_2i*V1+C7_2i*V7                        
C4_2i=1/V4*((TC_2i*10**22)-TC) 

TA=A6_2i*V6+A2_2i*V2+A3_2i*V3+A5_2i*V5+A1_2i*V1+A7_2i*V7                        
A4_2i=1/V4*((3446475.31695093e+12-(3358734.46282621e+12-TC_2i*(10**22))*2)-TA) 

TP=P6_2i*V6+P2_2i*V2+P3_2i*V3+P5_2i*V5+P1_2i*V1+P7_2i*V7                        
P4_2i=1/V4*((2.8733503182616036*10**15)-TP)

#%%
for i in range(20):
    x=np.where(HB1[i]<0)
    if np.size(x)>0:
        print(i)
#%%       
CO2=CO2_2i*10**6
Tadjust=np.zeros((np.size(CO2),1))
for j in range (np.size(CO2)):
    Tadjust[j]=3.25*0.54*5.35*np.log((CO2[j])/145) # Adjustment temperature CO2 concentration

T7=0.33+Tadjust   
x_2=np.where(T7<-1)

if np.size(x_2)>0:
    print('2')

eps_TEMP=0.9
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
figure1 = plt.figure(1,figsize=(10,6.5))

plt.plot(ti_2,CO2_2i*10**6,linewidth=LW)

plt.xlim(0,max(ti_2))

plt.grid()
plt.title('Oscillation at 15Sv',fontsize=FS)
plt.ylabel('CO$_2$ [ppm]',fontsize=FS)
plt.xlabel('Time [yr]',fontsize=FS)
plt.xticks(fontsize=FS-3)
plt.yticks(fontsize=FS-3)
plt.ylim(160,240)
plt.savefig('figure5b.png', format='png', dpi=300)