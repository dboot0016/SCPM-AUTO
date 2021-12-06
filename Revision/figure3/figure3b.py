#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 15:01:57 2021

@author: daan
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 14:08:22 2021

@author: daan
"""

import pandas as pd

import numpy as np
import matplotlib.pyplot as plt
import math

from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm

# %%
name0='L-CTL'
name1='L-BIO'
name2='L-TEMP'
name3='L-PV'
name4='L-BALK'
name5='L-RAIN'
name6='L-BIO-TEMP'
name7='L-BIO-PV'
name8='L-BIO-BALK'
name9='L-BIO-RAIN'
name10='L-BIO-TEMP-PV'
name11='L-BIO-TEMP(A)-EFF'
name12='L-ALL'
# %%

data0 = pd.read_csv('L_0.txt',sep='\s+',header=None)
data1 = pd.read_csv('L_1.txt',sep='\s+',header=None)

data0 = pd.DataFrame(data0)
data1 = pd.DataFrame(data1)

AMOC_0=data0[4]*10**-6
AMOC_1=data1[4]*10**-6

CO2_0=data0[6]*10**6
CO2_1=data1[6]*10**6

xP0=7133       
xP1=12068  

#%%
data2 = pd.read_csv('L_2.txt',sep='\s+',header=None)
data3 = pd.read_csv('L_4.txt',sep='\s+',header=None)
data4 = pd.read_csv('L_5.txt',sep='\s+',header=None)
data5 = pd.read_csv('L_6.txt',sep='\s+',header=None)
data6 = pd.read_csv('L_11.txt',sep='\s+',header=None)
data7 = pd.read_csv('L_13.txt',sep='\s+',header=None)
data8 = pd.read_csv('L_14.txt',sep='\s+',header=None)
data9 = pd.read_csv('L_15.txt',sep='\s+',header=None)
data10 = pd.read_csv('L_22.txt',sep='\s+',header=None)
data11 = pd.read_csv('L_26.txt',sep='\s+',header=None)
data12 = pd.read_csv('L_ALL.txt',sep='\s+',header=None)

data2 = pd.DataFrame(data2)
data3 = pd.DataFrame(data3)
data4 = pd.DataFrame(data4)
data5 = pd.DataFrame(data5)
data6 = pd.DataFrame(data6)
data7 = pd.DataFrame(data7)
data8 = pd.DataFrame(data8)
data9 = pd.DataFrame(data9)
data10 = pd.DataFrame(data10)
data11 = pd.DataFrame(data11)
data12 = pd.DataFrame(data12)

AMOC_2=data2[4]*10**-6
AMOC_3=data3[4]*10**-6
AMOC_4=data4[4]*10**-6
AMOC_5=data5[4]*10**-6
AMOC_6=data6[4]*10**-6
AMOC_7=data7[4]*10**-6
AMOC_8=data8[4]*10**-6
AMOC_9=data9[4]*10**-6
AMOC_10=data10[4]*10**-6
AMOC_11=data11[4]*10**-6
AMOC_12=data12[4]*10**-6

CO2_2=data2[6]*10**6
CO2_3=data3[6]*10**6
CO2_4=data4[6]*10**6
CO2_5=data5[6]*10**6
CO2_6=data6[6]*10**6
CO2_7=data7[6]*10**6
CO2_8=data8[6]*10**6
CO2_9=data9[6]*10**6
CO2_10=data10[6]*10**6
CO2_11=data11[6]*10**6
CO2_12=data12[6]*10**6

#%%
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

#%%
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

#%%
C1_4=data4[7]
A1_4=data4[8]
P1_4=data4[9]
C2_4=data4[10]
A2_4=data4[11]
P2_4=data4[12]
C3_4=data4[13]
A3_4=data4[14]
P3_4=data4[15]
C5_4=data4[16]
A5_4=data4[17]
P5_4=data4[18]
C6_4=data4[19]
A6_4=data4[20]
P6_4=data4[21]
C7_4=data4[22]
A7_4=data4[23]
P7_4=data4[24]
TC_4=data4[25]

#%%
C1_5=data5[7]
A1_5=data5[8]
P1_5=data5[9]
C2_5=data5[10]
A2_5=data5[11]
P2_5=data5[12]
C3_5=data5[13]
A3_5=data5[14]
P3_5=data5[15]
C5_5=data5[16]
A5_5=data5[17]
P5_5=data5[18]
C6_5=data5[19]
A6_5=data5[20]
P6_5=data5[21]
C7_5=data5[22]
A7_5=data5[23]
P7_5=data5[24]
TC_5=data5[25]

#%%
C1_6=data6[7]
A1_6=data6[8]
P1_6=data6[9]
C2_6=data6[10]
A2_6=data6[11]
P2_6=data6[12]
C3_6=data6[13]
A3_6=data6[14]
P3_6=data6[15]
C5_6=data6[16]
A5_6=data6[17]
P5_6=data6[18]
C6_6=data6[19]
A6_6=data6[20]
P6_6=data6[21]
C7_6=data6[22]
A7_6=data6[23]
P7_6=data6[24]
TC_6=data6[25]

#%%
C1_7=data7[7]
A1_7=data7[8]
P1_7=data7[9]
C2_7=data7[10]
A2_7=data7[11]
P2_7=data7[12]
C3_7=data7[13]
A3_7=data7[14]
P3_7=data7[15]
C5_7=data7[16]
A5_7=data7[17]
P5_7=data7[18]
C6_7=data7[19]
A6_7=data7[20]
P6_7=data7[21]
C7_7=data7[22]
A7_7=data7[23]
P7_7=data7[24]
TC_7=data7[25]

#%%
C1_8=data8[7]
A1_8=data8[8]
P1_8=data8[9]
C2_8=data8[10]
A2_8=data8[11]
P2_8=data8[12]
C3_8=data8[13]
A3_8=data8[14]
P3_8=data8[15]
C5_8=data8[16]
A5_8=data8[17]
P5_8=data8[18]
C6_8=data8[19]
A6_8=data8[20]
P6_8=data8[21]
C7_8=data8[22]
A7_8=data8[23]
P7_8=data8[24]
TC_8=data8[25]

#%%
C1_9=data9[7]
A1_9=data9[8]
P1_9=data9[9]
C2_9=data9[10]
A2_9=data9[11]
P2_9=data9[12]
C3_9=data9[13]
A3_9=data9[14]
P3_9=data9[15]
C5_9=data9[16]
A5_9=data9[17]
P5_9=data9[18]
C6_9=data9[19]
A6_9=data9[20]
P6_9=data9[21]
C7_9=data9[22]
A7_9=data9[23]
P7_9=data9[24]
TC_9=data9[25]

#%%
C1_10=data10[7]
A1_10=data10[8]
P1_10=data10[9]
C2_10=data10[10]
A2_10=data10[11]
P2_10=data10[12]
C3_10=data10[13]
A3_10=data10[14]
P3_10=data10[15]
C5_10=data10[16]
A5_10=data10[17]
P5_10=data10[18]
C6_10=data10[19]
A6_10=data10[20]
P6_10=data10[21]
C7_10=data10[22]
A7_10=data10[23]
P7_10=data10[24]
TC_10=data10[25]

#%%
C1_11=data11[7]
A1_11=data11[8]
P1_11=data11[9]
C2_11=data11[10]
A2_11=data11[11]
P2_11=data11[12]
C3_11=data11[13]
A3_11=data11[14]
P3_11=data11[15]
C5_11=data11[16]
A5_11=data11[17]
P5_11=data11[18]
C6_11=data11[19]
A6_11=data11[20]
P6_11=data11[21]
C7_11=data11[22]
A7_11=data11[23]
P7_11=data11[24]
TC_11=data11[25]

#%%
C1_12=data12[7]
A1_12=data12[8]
P1_12=data12[9]
C2_12=data12[10]
A2_12=data12[11]
P2_12=data12[12]
C3_12=data12[13]
A3_12=data12[14]
P3_12=data12[15]
C5_12=data12[16]
A5_12=data12[17]
P5_12=data12[18]
C6_12=data12[19]
A6_12=data12[20]
P6_12=data12[21]
C7_12=data12[22]
A7_12=data12[23]
P7_12=data12[24]
TC_12=data12[25]

#%%
for i in range(19):
   if np.size(data2.loc[data2[i+6] < 0])>0:
       print('2')
       print(i+6)
   if np.size(data3.loc[data3[i+6] < 0])>0:
       print('3')
       print(i+6)
   if np.size(data4.loc[data4[i+6] < 0])>0:
       print('4')
       print(i+6)
   if np.size(data5.loc[data5[i+6] < 0])>0:
       print('5')
       print(i+6)
   if np.size(data6.loc[data6[i+6] < 0])>0:
       print('6')
       print(i+6)
   if np.size(data7.loc[data7[i+6] < 0])>0:
       print('7')
       print(i+6)
   if np.size(data8.loc[data8[i+6] < 0])>0:
       print('8')
       print(i+6)
   if np.size(data9.loc[data9[i+6] < 0])>0:
       print('9')
       print(i+6)
   if np.size(data10.loc[data10[i+6] < 0])>0:
       print('10')
       print(i+6)
   if np.size(data11.loc[data11[i+6] < 0])>0:
       print('11')
       print(i+6)

#%%
CO2=CO2_2
Tadjust=np.zeros((np.size(CO2),1))
Tadjust5=np.zeros((np.size(CO2),1))
for j in range (np.size(CO2)):
    Tadjust[j]=3.25*0.54*5.35*np.log((CO2[j])/145) # Adjustment temperature CO2 concentration

T7=0.33+Tadjust   
x_2=np.where(T7<-1)

if np.size(x_2)>0:
    print('2')
#%%    
CO2=CO2_6
Tadjust=np.zeros((np.size(CO2),1))
Tadjust5=np.zeros((np.size(CO2),1))
for j in range (np.size(CO2)):
    Tadjust[j]=3.25*0.54*5.35*np.log((CO2[j])/145) # Adjustment temperature CO2 concentration

T7=0.33+Tadjust   
x_2=np.where(T7<-1)

if np.size(x_2)>0:
    print('6')

#%%
CO2=CO2_10
Tadjust=np.zeros((np.size(CO2),1))
Tadjust5=np.zeros((np.size(CO2),1))
for j in range (np.size(CO2)):
    Tadjust[j]=3.25*0.54*5.35*np.log((CO2[j])/145) # Adjustment temperature CO2 concentration

T7=0.33+Tadjust   
x_2=np.where(T7<-1)

if np.size(x_2)>0:
    print('10')
#%%
CO2=CO2_11
Tadjust=np.zeros((np.size(CO2),1))
Tadjust5=np.zeros((np.size(CO2),1))
for j in range (np.size(CO2)):
    Tadjust[j]=3.25*0.54*5.35*np.log((CO2[j])/145) # Adjustment temperature CO2 concentration

T7=0.33+Tadjust   
x_2=np.where(T7<-1)

if np.size(x_2)>0:
    print('11')
    
eps_TEMP=0.9
eps1_base=0.9
eps2_base=1.25
eps7_base=0.62
eps5_base=0.35
    
eps1=Tadjust*(-0.1)+eps1_base
eps2=Tadjust*(-0.1)+eps2_base
eps5=Tadjust*(-0.1)+eps5_base
eps7=Tadjust*(-0.1)+eps7_base
    
if np.size(np.where(eps1<0))>0:
    print('eps1')
if np.size(np.where(eps2<0))>0:
    print('eps2')
if np.size(np.where(eps5<0))>0:
    print('eps5')
if np.size(np.where(eps7<0))>0:
    print('eps7')

#%%
CO2=CO2_12
Tadjust=np.zeros((np.size(CO2),1))
Tadjust5=np.zeros((np.size(CO2),1))
for j in range (np.size(CO2)):
    Tadjust[j]=3.25*0.54*5.35*np.log((CO2[j])/145) # Adjustment temperature CO2 concentration

T7=0.33+Tadjust   
x_2=np.where(T7<-1)

if np.size(x_2)>0:
    print('12')
    
eps_TEMP=0.9
eps1_base=0.9
eps2_base=1.25
eps7_base=0.62
eps5_base=0.35
    
eps1=Tadjust*(-0.1)+eps1_base
eps2=Tadjust*(-0.1)+eps2_base
eps5=Tadjust*(-0.1)+eps5_base
eps7=Tadjust*(-0.1)+eps7_base
    
if np.size(np.where(eps1<0))>0:
    print('eps1')
if np.size(np.where(eps2<0))>0:
    print('eps2')
if np.size(np.where(eps5<0))>0:
    print('eps5')
if np.size(np.where(eps7<0))>0:
    print('eps7')

#%%
xP0=6900       
xP1=12055   
xP2=6900       
xP3=6900  
xP4=6900       
xP5=6900  
xP6=11787       
xP7=np.size(CO2_7)  
xP8=12062       
xP9=np.size(CO2_9)
xP10=11769#np.size(CO2_10)       
xP11=14736 
xP12=14530

#%%
LW=15
FS=25

#%%
t0=CO2_0[:xP0]
x0=AMOC_0[:xP0]
y0=np.zeros(np.size(x0))+1

t1=CO2_1[:xP1]
x1=AMOC_1[:xP1]
y1=np.zeros(np.size(x1))+1.005

t2=CO2_2[:xP2]
x2=AMOC_2[:xP2]
y2=np.zeros(np.size(x2))+1.001

t3=CO2_3[:xP3]
x3=AMOC_3[:xP3]
y3=np.zeros(np.size(x3))+1.002

t4=CO2_4[:xP4]
x4=AMOC_4[:xP4]
y4=np.zeros(np.size(x4))+1.003

t5=CO2_5[:xP5]
x5=AMOC_5[:xP5]
y5=np.zeros(np.size(x5))+1.004

t6=CO2_6[:xP6]
x6=AMOC_6[:xP6]
y6=np.zeros(np.size(x6))+1.006

t7=CO2_7[:xP7]
x7=AMOC_7[:xP7]
y7=np.zeros(np.size(x7))+1.007

t8=CO2_8[:xP8]
x8=AMOC_8[:xP8]
y8=np.zeros(np.size(x8))+1.008

t9=CO2_9[:xP9]
x9=AMOC_9[:xP9]
y9=np.zeros(np.size(x9))+1.009

t10=CO2_10[:xP10]
x10=AMOC_10[:xP10]
y10=np.zeros(np.size(x10))+1.01

t11=CO2_11[:xP11]
x11=AMOC_11[:xP11]
y11=np.zeros(np.size(x11))+1.011

t12=CO2_12[:xP12]
x12=AMOC_12[:xP12]
y12=np.zeros(np.size(x12))+1.012

points0 = np.array([x0, y0]).T.reshape(-1, 1, 2)
segments0 = np.concatenate([points0[:-1], points0[1:]], axis=1)

points1 = np.array([x1, y1]).T.reshape(-1, 1, 2)
segments1 = np.concatenate([points1[:-1], points1[1:]], axis=1)

points2 = np.array([x2, y2]).T.reshape(-1, 1, 2)
segments2 = np.concatenate([points2[:-1], points2[1:]], axis=1)

points3 = np.array([x3, y3]).T.reshape(-1, 1, 2)
segments3 = np.concatenate([points3[:-1], points3[1:]], axis=1)

points4 = np.array([x4, y4]).T.reshape(-1, 1, 2)
segments4 = np.concatenate([points4[:-1], points4[1:]], axis=1)

points5 = np.array([x5, y5]).T.reshape(-1, 1, 2)
segments5 = np.concatenate([points5[:-1], points5[1:]], axis=1)

points6 = np.array([x6, y6]).T.reshape(-1, 1, 2)
segments6 = np.concatenate([points6[:-1], points6[1:]], axis=1)

points7 = np.array([x7, y7]).T.reshape(-1, 1, 2)
segments7 = np.concatenate([points7[:-1], points7[1:]], axis=1)

points8 = np.array([x8, y8]).T.reshape(-1, 1, 2)
segments8 = np.concatenate([points8[:-1], points8[1:]], axis=1)

points9 = np.array([x9, y9]).T.reshape(-1, 1, 2)
segments9 = np.concatenate([points9[:-1], points9[1:]], axis=1)

points10 = np.array([x10, y10]).T.reshape(-1, 1, 2)
segments10 = np.concatenate([points10[:-1], points10[1:]], axis=1)

points11 = np.array([x11, y11]).T.reshape(-1, 1, 2)
segments11 = np.concatenate([points11[:-1], points11[1:]], axis=1)

points12 = np.array([x12, y12]).T.reshape(-1, 1, 2)
segments12 = np.concatenate([points12[:-1], points12[1:]], axis=1)

lc0 = LineCollection(segments0, cmap=plt.get_cmap('RdYlBu_r'),
    norm=plt.Normalize(120,200))
lc0.set_array(t0)
lc0.set_linewidth(LW)

lc1 = LineCollection(segments1, cmap=plt.get_cmap('RdYlBu_r'),
    norm=plt.Normalize(120,200))
lc1.set_array(t1)
lc1.set_linewidth(LW)

lc2 = LineCollection(segments2, cmap=plt.get_cmap('RdYlBu_r'),
    norm=plt.Normalize(120,200))
lc2.set_array(t2)
lc2.set_linewidth(LW)

lc3 = LineCollection(segments3, cmap=plt.get_cmap('RdYlBu_r'),
    norm=plt.Normalize(120,200))
lc3.set_array(t3)
lc3.set_linewidth(LW)

lc4 = LineCollection(segments4, cmap=plt.get_cmap('RdYlBu_r'),
    norm=plt.Normalize(120,200))
lc4.set_array(t4)
lc4.set_linewidth(LW)

lc5 = LineCollection(segments5, cmap=plt.get_cmap('RdYlBu_r'),
    norm=plt.Normalize(120,200))
lc5.set_array(t5)
lc5.set_linewidth(LW)

lc6 = LineCollection(segments6, cmap=plt.get_cmap('RdYlBu_r'),
    norm=plt.Normalize(120,200))
lc6.set_array(t6)
lc6.set_linewidth(LW)

lc7 = LineCollection(segments7, cmap=plt.get_cmap('RdYlBu_r'),
    norm=plt.Normalize(120,200))
lc7.set_array(t7)
lc7.set_linewidth(LW)

lc8 = LineCollection(segments8, cmap=plt.get_cmap('RdYlBu_r'),
    norm=plt.Normalize(120,200))
lc8.set_array(t8)
lc8.set_linewidth(LW)

lc9 = LineCollection(segments9, cmap=plt.get_cmap('RdYlBu_r'),
    norm=plt.Normalize(120,200))
lc9.set_array(t9)
lc9.set_linewidth(LW)

lc10 = LineCollection(segments10, cmap=plt.get_cmap('RdYlBu_r'),
    norm=plt.Normalize(120,200))
lc10.set_array(t10)
lc10.set_linewidth(LW)

lc11 = LineCollection(segments11, cmap=plt.get_cmap('RdYlBu_r'),
    norm=plt.Normalize(120,200))
lc11.set_array(t11)
lc11.set_linewidth(LW)

lc12 = LineCollection(segments12, cmap=plt.get_cmap('RdYlBu_r'),
    norm=plt.Normalize(120,200))
lc12.set_array(t12)
lc12.set_linewidth(LW)

locs=[1,1.001,1.002,1.003,1.004,1.005,1.006,1.007,1.008,1.009,1.01,1.011,1.012]

xline=np.linspace(0,30,np.size(y7))
x_A=[15,15]
y_A=[0,10]

fig, ax = plt.subplots(figsize=(10,6.5))
plt.plot(x_A,y_A,'k',linewidth=2.5, zorder=1)
im=plt.gca().add_collection(lc0)
plt.gca().add_collection(lc1)
plt.gca().add_collection(lc2)
plt.gca().add_collection(lc3)
plt.gca().add_collection(lc4)
plt.gca().add_collection(lc5)
plt.gca().add_collection(lc6)
plt.gca().add_collection(lc7)
plt.gca().add_collection(lc8)
plt.gca().add_collection(lc9)
plt.gca().add_collection(lc10)
plt.gca().add_collection(lc11)
plt.gca().add_collection(lc12)

cbar=plt.colorbar(im,ax=ax,extend='both')
cbar.set_label('CO$_2$ [ppm]',fontsize=FS)
cbar.ax.tick_params(labelsize=FS-3)
plt.xlim(0, 30)
plt.ylim(0.999, 1.013)
plt.yticks(locs,[name0,name2,name3,name4,name5,name1,name6,name7,name8,name9,name10,name11,name12],fontsize=FS-3)
plt.xlabel('AMOC [Sv]',fontsize=FS)
plt.xticks(fontsize=FS-3)
plt.title('Last Glacial Maximum configuration',fontsize=FS)
plt.grid()
plt.show()   
fig.savefig('figure3b_rev.png', format='png', dpi=300, bbox_inches = 'tight')