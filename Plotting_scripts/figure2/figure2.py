#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 13:21:32 2021

@author: daan
"""

import pandas as pd

import numpy as np
import matplotlib.pyplot as plt
import math
# %%
name0='P-0'
name1='P-1'
name2='L-0'
name3='L-1'

# %%

data0 = pd.read_csv('P_0.txt',sep='\s+',header=None)
data1 = pd.read_csv('P_1.txt',sep='\s+',header=None)

data0 = pd.DataFrame(data0)
data1 = pd.DataFrame(data1)

AMOC_0=data0[4]*10**-6
AMOC_1=data1[4]*10**-6

CO2_0=data0[6]*10**6
CO2_1=data1[6]*10**6

C1_0=data0[7]
A1_0=data0[8]
P1_0=data0[9]
C2_0=data0[10]
A2_0=data0[11]
P2_0=data0[12]
C3_0=data0[13]
A3_0=data0[14]
P3_0=data0[15]
C5_0=data0[16]
A5_0=data0[17]
P5_0=data0[18]
C6_0=data0[19]
A6_0=data0[20]
P6_0=data0[21]
C7_0=data0[22]
A7_0=data0[23]
P7_0=data0[24]
TC_0=data0[25]

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

for i in range(19):
   if np.size(data0.loc[data0[i+6] < 0])>0:
       print('0')
       print(i+6)
   if np.size(data1.loc[data1[i+6] < 0])>0:
       print('1')
       print(i+6)

xP0=7133       
xP1=12068     


#%%
data2 = pd.read_csv('L_0.txt',sep='\s+',header=None)
data3 = pd.read_csv('L_1.txt',sep='\s+',header=None)

data2 = pd.DataFrame(data2)
data3 = pd.DataFrame(data3)

AMOC_2=data2[4]*10**-6
AMOC_3=data3[4]*10**-6

CO2_2=data2[6]*10**6
CO2_3=data3[6]*10**6

C1_2=data2[7]
A1_2=data2[8]
P1_2=data2[9]
C2_2=data2[10]
A2_2=data2[11]
P2_2=data2[12]
C3_2=data2[13]
A3_2=data2[14]
P3_2=data2[15]
C5_2=data2[16]
A5_2=data2[17]
P5_2=data2[18]
C6_2=data2[19]
A6_2=data2[20]
P6_2=data2[21]
C7_2=data2[22]
A7_2=data2[23]
P7_2=data2[24]
TC_2=data2[25]

C1_3=data3[7]
A1_3=data3[8]
P1_3=data3[9]
C2_3=data3[10]
A2_3=data3[11]
P2_3=data3[12]
C3_3=data3[13]
A3_3=data3[14]
P3_3=data3[15]
C5_3=data3[16]
A5_3=data3[17]
P5_3=data3[18]
C6_3=data3[19]
A6_3=data3[20]
P6_3=data3[21]
C7_3=data3[22]
A7_3=data3[23]
P7_3=data3[24]
TC_3=data3[25]

for i in range(19):
   if np.size(data2.loc[data2[i+6] < 0])>0:
       print('2')
       print(i+6)
   if np.size(data3.loc[data3[i+6] < 0])>0:
       print('3')
       print(i+6)

xP2=6900       
xP3=12055     

#%%
LW=5
FS=25
#%%
figure1 = plt.figure(1,figsize=(10,6.5))

plt.plot(AMOC_0[:xP0],CO2_0[:xP0],'tab:blue',linewidth=LW)
plt.plot(AMOC_1[:xP1],CO2_1[:xP1],'tab:red',linewidth=LW)
plt.plot(AMOC_2[:xP2],CO2_2[:xP2],linestyle='--',color='tab:blue',linewidth=LW)
plt.plot(AMOC_3[:xP3],CO2_3[:xP3],linestyle='--',color='tab:red',linewidth=LW)

plt.grid()
#plt.title('AMOC vs CO$_2$',fontsize=FS)
plt.ylabel('CO$_2$ [ppm]',fontsize=FS)
plt.xlabel('AMOC [Sv]',fontsize=FS)
plt.xticks(fontsize=FS-3)
plt.yticks(fontsize=FS-3)
plt.legend([name0,name1,name2,name3],loc="lower center",fontsize=FS-6,ncol=4)
plt.xlim(0,30)
plt.ylim(0,300)
#plt.savefig('figure2.png', format='png', dpi=300)

#%%
plt.plot(AMOC_1,CO2_0[:xP0])
plt.plot(AMOC_1,CO2_2[:xP2]+50)