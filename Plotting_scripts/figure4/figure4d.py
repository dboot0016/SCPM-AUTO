#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 11:15:59 2021

@author: daan
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 10:18:14 2021

@author: daan
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 08:46:54 2021

@author: daan
"""

import pandas as pd

import numpy as np
import matplotlib.pyplot as plt
import math

from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm

#%%
name0='LC-0'
name1='LC-1'
name2='LC-6'
name3='LC-9'
name4='LC-10'
name5='LC-11'

#%%
LW=12.5
FS=25
log=1

#%%
data0 = pd.read_csv('L_A_Z_0a.txt',sep='\s+',header=None)
data1 = pd.read_csv('L_A_Z_1a.txt',sep='\s+',header=None)
data2 = pd.read_csv('L_A_Z_2a.txt',sep='\s+',header=None)
data3 = pd.read_csv('L_A_Z_3a.txt',sep='\s+',header=None)
data4 = pd.read_csv('L_A_Z_4a.txt',sep='\s+',header=None)
data5 = pd.read_csv('L_A_Z_5a.txt',sep='\s+',header=None)

data0 = pd.DataFrame(data0)
data1 = pd.DataFrame(data1)
data2 = pd.DataFrame(data2)
data3 = pd.DataFrame(data3)
data4 = pd.DataFrame(data4)
data5 = pd.DataFrame(data5)

data10 = pd.read_csv('L_A_Z_0b.txt',sep='\s+',header=None)
data11 = pd.read_csv('L_A_Z_1b.txt',sep='\s+',header=None)
data12 = pd.read_csv('L_A_Z_2b.txt',sep='\s+',header=None)
data13 = pd.read_csv('L_A_Z_3b.txt',sep='\s+',header=None)
data14 = pd.read_csv('L_A_Z_4b.txt',sep='\s+',header=None)
data15 = pd.read_csv('L_A_Z_5b.txt',sep='\s+',header=None)

data10 = pd.DataFrame(data10)
data11 = pd.DataFrame(data11)
data12 = pd.DataFrame(data12)
data13 = pd.DataFrame(data13)
data14 = pd.DataFrame(data14)
data15 = pd.DataFrame(data15)

data20 = pd.read_csv('L_A_Z_0e.txt',sep='\s+',header=None)
data21 = pd.read_csv('L_A_Z_1e.txt',sep='\s+',header=None)
data22 = pd.read_csv('L_A_Z_2e.txt',sep='\s+',header=None)
data23 = pd.read_csv('L_A_Z_3e.txt',sep='\s+',header=None)
data24 = pd.read_csv('L_A_Z_4e.txt',sep='\s+',header=None)
data25 = pd.read_csv('L_A_Z_5e.txt',sep='\s+',header=None)

data20 = pd.DataFrame(data20)
data21 = pd.DataFrame(data21)
data22 = pd.DataFrame(data22)
data23 = pd.DataFrame(data23)
data24 = pd.DataFrame(data24)
data25 = pd.DataFrame(data25)

PAR_0=10**data0[4]
PAR_1=10**data1[4]
PAR_2=10**data2[4]
PAR_3=10**data3[4]
PAR_4=10**data4[4]
PAR_5=10**data5[4]

CO2_0=data0[6]*10**6
CO2_1=data1[6]*10**6
CO2_2=data2[6]*10**6
CO2_3=data3[6]*10**6
CO2_4=data4[6]*10**6
CO2_5=data5[6]*10**6

PAR_10=10**data10[4]
PAR_11=10**data11[4]
PAR_12=10**data12[4]
PAR_13=10**data13[4]
PAR_14=10**data14[4]
PAR_15=10**data15[4]

CO2_10=data10[6]*10**6
CO2_11=data11[6]*10**6
CO2_12=data12[6]*10**6
CO2_13=data13[6]*10**6
CO2_14=data14[6]*10**6
CO2_15=data15[6]*10**6

PAR_20=10**data20[4]
PAR_21=10**data21[4]
PAR_22=10**data22[4]
PAR_23=10**data23[4]
PAR_24=10**data24[4]
PAR_25=10**data25[4]

CO2_20=data20[6]*10**6
CO2_21=data21[6]*10**6
CO2_22=data22[6]*10**6
CO2_23=data23[6]*10**6
CO2_24=data24[6]*10**6
CO2_25=data25[6]*10**6

#%%
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

#%%
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
C1_13=data13[7]
A1_13=data13[8]
P1_13=data13[9]
C2_13=data13[10]
A2_13=data13[11]
P2_13=data13[12]
C3_13=data13[13]
A3_13=data13[14]
P3_13=data13[15]
C5_13=data13[16]
A5_13=data13[17]
P5_13=data13[18]
C6_13=data13[19]
A6_13=data13[20]
P6_13=data13[21]
C7_13=data13[22]
A7_13=data13[23]
P7_13=data13[24]
TC_13=data13[25]

#%%
C1_14=data14[7]
A1_14=data14[8]
P1_14=data14[9]
C2_14=data14[10]
A2_14=data14[11]
P2_14=data14[12]
C3_14=data14[13]
A3_14=data14[14]
P3_14=data14[15]
C5_14=data14[16]
A5_14=data14[17]
P5_14=data14[18]
C6_14=data14[19]
A6_14=data14[20]
P6_14=data14[21]
C7_14=data14[22]
A7_14=data14[23]
P7_14=data14[24]
TC_14=data14[25]

#%%
C1_15=data15[7]
A1_15=data15[8]
P1_15=data15[9]
C2_15=data15[10]
A2_15=data15[11]
P2_15=data15[12]
C3_15=data15[13]
A3_15=data15[14]
P3_15=data15[15]
C5_15=data15[16]
A5_15=data15[17]
P5_15=data15[18]
C6_15=data15[19]
A6_15=data15[20]
P6_15=data15[21]
C7_15=data15[22]
A7_15=data15[23]
P7_15=data15[24]
TC_15=data15[25]

#%%
C1_20=data20[7]
A1_20=data20[8]
P1_20=data20[9]
C2_20=data20[10]
A2_20=data20[11]
P2_20=data20[12]
C3_20=data20[13]
A3_20=data20[14]
P3_20=data20[15]
C5_20=data20[16]
A5_20=data20[17]
P5_20=data20[18]
C6_20=data20[19]
A6_20=data20[20]
P6_20=data20[21]
C7_20=data20[22]
A7_20=data20[23]
P7_20=data20[24]
TC_20=data20[25]

#%%
C1_21=data21[7]
A1_21=data21[8]
P1_21=data21[9]
C2_21=data21[10]
A2_21=data21[11]
P2_21=data21[12]
C3_21=data21[13]
A3_21=data21[14]
P3_21=data21[15]
C5_21=data21[16]
A5_21=data21[17]
P5_21=data21[18]
C6_21=data21[19]
A6_21=data21[20]
P6_21=data21[21]
C7_21=data21[22]
A7_21=data21[23]
P7_21=data21[24]
TC_21=data21[25]

#%%
C1_22=data22[7]
A1_22=data22[8]
P1_22=data22[9]
C2_22=data22[10]
A2_22=data22[11]
P2_22=data22[12]
C3_22=data22[13]
A3_22=data22[14]
P3_22=data22[15]
C5_22=data22[16]
A5_22=data22[17]
P5_22=data22[18]
C6_22=data22[19]
A6_22=data22[20]
P6_22=data22[21]
C7_22=data22[22]
A7_22=data22[23]
P7_22=data22[24]
TC_22=data22[25]

#%%
C1_23=data23[7]
A1_23=data23[8]
P1_23=data23[9]
C2_23=data23[10]
A2_23=data23[11]
P2_23=data23[12]
C3_23=data23[13]
A3_23=data23[14]
P3_23=data23[15]
C5_23=data23[16]
A5_23=data23[17]
P5_23=data23[18]
C6_23=data23[19]
A6_23=data23[20]
P6_23=data23[21]
C7_23=data23[22]
A7_23=data23[23]
P7_23=data23[24]
TC_23=data23[25]

#%%
C1_24=data24[7]
A1_24=data24[8]
P1_24=data24[9]
C2_24=data24[10]
A2_24=data24[11]
P2_24=data24[12]
C3_24=data24[13]
A3_24=data24[14]
P3_24=data24[15]
C5_24=data24[16]
A5_24=data24[17]
P5_24=data24[18]
C6_24=data24[19]
A6_24=data24[20]
P6_24=data24[21]
C7_24=data24[22]
A7_24=data24[23]
P7_24=data24[24]
TC_24=data24[25]

#%%
C1_25=data25[7]
A1_25=data25[8]
P1_25=data25[9]
C2_25=data25[10]
A2_25=data25[11]
P2_25=data25[12]
C3_25=data25[13]
A3_25=data25[14]
P3_25=data25[15]
C5_25=data25[16]
A5_25=data25[17]
P5_25=data25[18]
C6_25=data25[19]
A6_25=data25[20]
P6_25=data25[21]
C7_25=data25[22]
A7_25=data25[23]
P7_25=data25[24]
TC_25=data25[25]

#%%
for i in range(19):
   if np.size(data0.loc[data0[i+6] < 0])>0:
       print('0')
       print(i+6)
   if np.size(data1.loc[data1[i+6] < 0])>0:
       print('1')
       print(i+6)
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
   if np.size(data10.loc[data10[i+6] < 0])>0:
       print('10')
       print(i+6)
   if np.size(data11.loc[data11[i+6] < 0])>0:
       print('11')
       print(i+6)
   if np.size(data12.loc[data12[i+6] < 0])>0:
       print('12')
       print(i+6)
   if np.size(data13.loc[data13[i+6] < 0])>0:
       print('13')
       print(i+6)
   if np.size(data14.loc[data14[i+6] < 0])>0:
       print('14')
       print(i+6)
   if np.size(data15.loc[data15[i+6] < 0])>0:
       print('15')
       print(i+6)
   if np.size(data20.loc[data20[i+6] < 0])>0:
       print('20')
       print(i+6)
   if np.size(data21.loc[data21[i+6] < 0])>0:
       print('21')
       print(i+6)
   if np.size(data22.loc[data22[i+6] < 0])>0:
       print('22')
       print(i+6)
   if np.size(data23.loc[data23[i+6] < 0])>0:
       print('23')
       print(i+6)
   if np.size(data24.loc[data24[i+6] < 0])>0:
       print('24')
       print(i+6)
   if np.size(data25.loc[data25[i+6] < 0])>0:
       print('25')
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
    print(x_2)

CO2=CO2_3
Tadjust=np.zeros((np.size(CO2),1))
Tadjust5=np.zeros((np.size(CO2),1))
for j in range (np.size(CO2)):
    Tadjust[j]=3.25*0.54*5.35*np.log((CO2[j])/145) # Adjustment temperature CO2 concentration

T7=0.33+Tadjust   
x_2=np.where(T7<-1)

if np.size(x_2)>0:
    print('3')
    print(x_2)
  
CO2=CO2_5
Tadjust=np.zeros((np.size(CO2),1))
Tadjust5=np.zeros((np.size(CO2),1))
for j in range (np.size(CO2)):
    Tadjust[j]=3.25*0.54*5.35*np.log((CO2[j])/145) # Adjustment temperature CO2 concentration

T7=0.33+Tadjust   
x_2=np.where(T7<-1)

if np.size(x_2)>0:
    print('5')
    print(x_2)
    
eps_TEMP=1
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
    print(x_2)

CO2=CO2_13
Tadjust=np.zeros((np.size(CO2),1))
Tadjust5=np.zeros((np.size(CO2),1))
for j in range (np.size(CO2)):
    Tadjust[j]=3.25*0.54*5.35*np.log((CO2[j])/145) # Adjustment temperature CO2 concentration

T7=0.33+Tadjust   
x_2=np.where(T7<-1)

if np.size(x_2)>0:
    print('13')
    print(x_2)
    
CO2=CO2_15
Tadjust=np.zeros((np.size(CO2),1))
Tadjust5=np.zeros((np.size(CO2),1))
for j in range (np.size(CO2)):
    Tadjust[j]=3.25*0.54*5.35*np.log((CO2[j])/145) # Adjustment temperature CO2 concentration

T7=0.33+Tadjust   
x_2=np.where(T7<-1)

if np.size(x_2)>0:
    print('15')
    print(x_2)
    
eps_TEMP=1
eps1_base=0.9
eps2_base=1.25
eps7_base=0.62
eps5_base=0.35
    
eps1=Tadjust*(-0.1)+eps1_base
eps2=Tadjust*(-0.1)+eps2_base
eps5=Tadjust*(-0.1)+eps5_base
eps7=Tadjust*(-0.1)+eps7_base
    
if np.size(np.where(eps1<0))>0:
    print('1eps1')
if np.size(np.where(eps2<0))>0:
    print('1eps2')
if np.size(np.where(eps5<0))>0:
    print('1eps5')
if np.size(np.where(eps7<0))>0:
    print('1eps7')

#%%    
CO2=CO2_22
Tadjust=np.zeros((np.size(CO2),1))
Tadjust5=np.zeros((np.size(CO2),1))
for j in range (np.size(CO2)):
    Tadjust[j]=3.25*0.54*5.35*np.log((CO2[j])/145) # Adjustment temperature CO2 concentration

T7=0.33+Tadjust   
x_2=np.where(T7<-1)

if np.size(x_2)>0:
    print('22')
    print(x_2)

CO2=CO2_23
Tadjust=np.zeros((np.size(CO2),1))
Tadjust5=np.zeros((np.size(CO2),1))
for j in range (np.size(CO2)):
    Tadjust[j]=3.25*0.54*5.35*np.log((CO2[j])/145) # Adjustment temperature CO2 concentration

T7=0.33+Tadjust   
x_2=np.where(T7<-1)

if np.size(x_2)>0:
    print('23')
    print(x_2)
   
CO2=CO2_25
Tadjust=np.zeros((np.size(CO2),1))
Tadjust5=np.zeros((np.size(CO2),1))
for j in range (np.size(CO2)):
    Tadjust[j]=3.25*0.54*5.35*np.log((CO2[j])/145) # Adjustment temperature CO2 concentration

T7=0.33+Tadjust   
x_2=np.where(T7<-1)

if np.size(x_2)>0:
    print('25')
    print(x_2)
    
eps_TEMP=1
eps1_base=0.9
eps2_base=1.25
eps7_base=0.62
eps5_base=0.35
    
eps1=Tadjust*(-0.1)+eps1_base
eps2=Tadjust*(-0.1)+eps2_base
eps5=Tadjust*(-0.1)+eps5_base
eps7=Tadjust*(-0.1)+eps7_base
    
if np.size(np.where(eps1<0))>0:
    print('2eps1')
if np.size(np.where(eps2<0))>0:
    print('2eps2')
if np.size(np.where(eps5<0))>0:
    print('2eps5')
if np.size(np.where(eps7<0))>0:
    print('2eps7')
#%%
xP0=np.size(CO2_0)       
xP1=np.size(CO2_1)  
xP2=130      
xP3=126
xP4=np.size(CO2_4)       
xP5=127 

xP001=0
xP011=0
xP021=0
xP031=0
xP041=0 
xP051=0

xP10=np.size(CO2_10)       
xP11=np.size(CO2_11)  
xP12=132       
xP13=129  
xP14=np.size(CO2_14)       
xP15=129  

xP101=1
xP111=0
xP121=0
xP131=0
xP141=0 
xP151=0

xP20=np.size(CO2_20)       
xP21=np.size(CO2_21)  
xP22=137    
xP23=135
xP24=np.size(CO2_24)       
xP25=135 

xP201=0
xP211=1
xP221=1
xP231=0
xP241=1 
xP251=0

#%%
t0=CO2_0[xP001:xP0]
x0=PAR_0[xP001:xP0]
y0=np.zeros(np.size(x0))+1

t10=CO2_10[xP101:xP10]
x10=PAR_10[xP101:xP10]
y10=np.zeros(np.size(x10))+1.001

t20=CO2_20[xP201:xP20]
x20=PAR_20[xP201:xP20]
y20=np.zeros(np.size(x20))+1.002

t1=CO2_1[xP011:xP1]
x1=PAR_1[xP011:xP1]
y1=np.zeros(np.size(x1))+1.004

t11=CO2_11[xP111:xP11]
x11=PAR_11[xP111:xP11]
y11=np.zeros(np.size(x11))+1.005

t21=CO2_21[xP211:xP21]
x21=PAR_21[xP211:xP21]
y21=np.zeros(np.size(x21))+1.006

t2=CO2_2[xP021:xP2]
x2=PAR_2[xP021:xP2]
y2=np.zeros(np.size(x2))+1.008

t12=CO2_12[xP121:xP12]
x12=PAR_12[xP121:xP12]
y12=np.zeros(np.size(x12))+1.009

t22=CO2_22[xP221:xP22]
x22=PAR_22[xP221:xP22]
y22=np.zeros(np.size(x22))+1.01

t3=CO2_3[xP031:xP3]
x3=PAR_3[xP031:xP3]
y3=np.zeros(np.size(x3))+1.012

t13=CO2_13[xP131:xP13]
x13=PAR_13[xP131:xP13]
y13=np.zeros(np.size(x13))+1.013

t23=CO2_23[xP231:xP23]
x23=PAR_23[xP231:xP23]
y23=np.zeros(np.size(x23))+1.014

t4=CO2_4[xP041:xP4]
x4=PAR_4[xP041:xP4]
y4=np.zeros(np.size(x4))+1.016

t14=CO2_14[xP141:xP14]
x14=PAR_14[xP141:xP14]
y14=np.zeros(np.size(x14))+1.017

t24=CO2_24[xP241:xP24]
x24=PAR_24[xP241:xP24]
y24=np.zeros(np.size(x24))+1.018

t5=CO2_5[xP051:xP5]
x5=PAR_5[xP051:xP5]
y5=np.zeros(np.size(x5))+1.02

t15=CO2_15[xP151:xP15]
x15=PAR_15[xP151:xP15]
y15=np.zeros(np.size(x15))+1.021

t25=CO2_25[xP251:xP25]
x25=PAR_25[xP251:xP25]
y25=np.zeros(np.size(x25))+1.022

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

points10 = np.array([x10, y10]).T.reshape(-1, 1, 2)
segments10 = np.concatenate([points10[:-1], points10[1:]], axis=1)

points11 = np.array([x11, y11]).T.reshape(-1, 1, 2)
segments11 = np.concatenate([points11[:-1], points11[1:]], axis=1)

points12 = np.array([x12, y12]).T.reshape(-1, 1, 2)
segments12 = np.concatenate([points12[:-1], points12[1:]], axis=1)

points13 = np.array([x13, y13]).T.reshape(-1, 1, 2)
segments13 = np.concatenate([points13[:-1], points13[1:]], axis=1)

points14 = np.array([x14, y14]).T.reshape(-1, 1, 2)
segments14 = np.concatenate([points14[:-1], points14[1:]], axis=1)

points15 = np.array([x15, y15]).T.reshape(-1, 1, 2)
segments15 = np.concatenate([points15[:-1], points15[1:]], axis=1)

points20 = np.array([x20, y20]).T.reshape(-1, 1, 2)
segments20 = np.concatenate([points20[:-1], points20[1:]], axis=1)

points21 = np.array([x21, y21]).T.reshape(-1, 1, 2)
segments21 = np.concatenate([points21[:-1], points21[1:]], axis=1)

points22 = np.array([x22, y22]).T.reshape(-1, 1, 2)
segments22 = np.concatenate([points22[:-1], points22[1:]], axis=1)

points23 = np.array([x23, y23]).T.reshape(-1, 1, 2)
segments23 = np.concatenate([points23[:-1], points23[1:]], axis=1)

points24 = np.array([x24, y24]).T.reshape(-1, 1, 2)
segments24 = np.concatenate([points24[:-1], points24[1:]], axis=1)

points25 = np.array([x25, y25]).T.reshape(-1, 1, 2)
segments25 = np.concatenate([points25[:-1], points25[1:]], axis=1)

lc0 = LineCollection(segments0, cmap=plt.get_cmap('RdYlBu_r'),
    norm=plt.Normalize(70,200))
lc0.set_array(t0)
lc0.set_linewidth(LW)

lc1 = LineCollection(segments1, cmap=plt.get_cmap('RdYlBu_r'),
    norm=plt.Normalize(70,200))
lc1.set_array(t1)
lc1.set_linewidth(LW)

lc2 = LineCollection(segments2, cmap=plt.get_cmap('RdYlBu_r'),
    norm=plt.Normalize(70,200))
lc2.set_array(t2)
lc2.set_linewidth(LW)

lc3 = LineCollection(segments3, cmap=plt.get_cmap('RdYlBu_r'),
    norm=plt.Normalize(70,200))
lc3.set_array(t3)
lc3.set_linewidth(LW)

lc4 = LineCollection(segments4, cmap=plt.get_cmap('RdYlBu_r'),
    norm=plt.Normalize(70,200))
lc4.set_array(t4)
lc4.set_linewidth(LW)

lc5 = LineCollection(segments5, cmap=plt.get_cmap('RdYlBu_r'),
    norm=plt.Normalize(70,200))
lc5.set_array(t5)
lc5.set_linewidth(LW)

lc10 = LineCollection(segments10, cmap=plt.get_cmap('RdYlBu_r'),
    norm=plt.Normalize(70,200))
lc10.set_array(t10)
lc10.set_linewidth(LW)

lc11 = LineCollection(segments11, cmap=plt.get_cmap('RdYlBu_r'),
    norm=plt.Normalize(70,200))
lc11.set_array(t11)
lc11.set_linewidth(LW)

lc12 = LineCollection(segments12, cmap=plt.get_cmap('RdYlBu_r'),
    norm=plt.Normalize(70,200))
lc12.set_array(t12)
lc12.set_linewidth(LW)

lc13 = LineCollection(segments13, cmap=plt.get_cmap('RdYlBu_r'),
    norm=plt.Normalize(70,200))
lc13.set_array(t13)
lc13.set_linewidth(LW)

lc14 = LineCollection(segments14, cmap=plt.get_cmap('RdYlBu_r'),
    norm=plt.Normalize(70,200))
lc14.set_array(t14)
lc14.set_linewidth(LW)

lc15 = LineCollection(segments15, cmap=plt.get_cmap('RdYlBu_r'),
    norm=plt.Normalize(70,200))
lc15.set_array(t15)
lc15.set_linewidth(LW)

lc20 = LineCollection(segments20, cmap=plt.get_cmap('RdYlBu_r'),
    norm=plt.Normalize(70,200))
lc20.set_array(t20)
lc20.set_linewidth(LW)

lc21 = LineCollection(segments21, cmap=plt.get_cmap('RdYlBu_r'),
    norm=plt.Normalize(70,200))
lc21.set_array(t21)
lc21.set_linewidth(LW)

lc22 = LineCollection(segments22, cmap=plt.get_cmap('RdYlBu_r'),
    norm=plt.Normalize(70,200))
lc22.set_array(t22)
lc22.set_linewidth(LW)

lc23 = LineCollection(segments23, cmap=plt.get_cmap('RdYlBu_r'),
    norm=plt.Normalize(70,200))
lc23.set_array(t23)
lc23.set_linewidth(LW)

lc24 = LineCollection(segments24, cmap=plt.get_cmap('RdYlBu_r'),
    norm=plt.Normalize(70,200))
lc24.set_array(t24)
lc24.set_linewidth(LW)

lc25 = LineCollection(segments25, cmap=plt.get_cmap('RdYlBu_r'),
    norm=plt.Normalize(70,200))
lc25.set_array(t25)
lc25.set_linewidth(LW)

locs=[1,1.001,1.002,1.004,1.005,1.006,1.008,1.009,1.01,1.012,1.013,1.014,1.016,1.017,1.018,1.02,1.021,1.022]

Y=np.zeros(np.size(x1))+1
xline=np.linspace(0.1,10,np.size(Y))

fig, ax = plt.subplots(figsize=(10,6.5))
im=plt.gca().add_collection(lc0)
plt.gca().add_collection(lc1)
plt.gca().add_collection(lc2)
plt.gca().add_collection(lc3)
plt.gca().add_collection(lc4)
plt.gca().add_collection(lc5)

plt.gca().add_collection(lc10)
plt.gca().add_collection(lc11)
plt.gca().add_collection(lc12)
plt.gca().add_collection(lc13)
plt.gca().add_collection(lc14)
plt.gca().add_collection(lc15)

plt.gca().add_collection(lc20)
plt.gca().add_collection(lc21)
plt.gca().add_collection(lc22)
plt.gca().add_collection(lc23)
plt.gca().add_collection(lc24)
plt.gca().add_collection(lc25)

plt.plot(xline,Y+0.003,'k',linewidth=0.75)
plt.plot(xline,Y+0.007,'k',linewidth=0.75)
plt.plot(xline,Y+0.011,'k',linewidth=0.75)
plt.plot(xline,Y+0.015,'k',linewidth=0.75)
plt.plot(xline,Y+0.019,'k',linewidth=0.75)

if log==1:
    plt.xscale('log')
    plt.text(5.5,0.99925,name0,fontsize=FS-7)
    plt.text(5.5,1.00325,name1,fontsize=FS-7)
    plt.text(5.5,1.00725,name2,fontsize=FS-7)
    plt.text(5.5,1.01125,name3,fontsize=FS-7)
    plt.text(5.5,1.01525,name4,fontsize=FS-7)
    plt.text(5.5,1.01925,name5,fontsize=FS-7)
else:
    plt.text(8,0.9995,name0)
    plt.text(8,1.0035,name1)
    plt.text(8,1.0075,name2)
    plt.text(8,1.0115,name3)
    plt.text(8,1.0155,name4)
    plt.text(8,1.0195,name5)

plt.title('LGM: Biological production',fontsize=FS)
cbar=plt.colorbar(im,ax=ax,extend='both')
cbar.set_label('CO$_2$ [ppm]',fontsize=FS)
cbar.ax.tick_params(labelsize=FS-3)
plt.xticks(fontsize=FS-3)
plt.xlim(0.1, 10)
plt.ylim(0.999, 1.023)
plt.yticks(locs,['0','1','4','0','1','4','0','1','4','0','1','4','0','1','4','0','1','4'],fontsize=FS-7)
plt.ylabel('$\lambda_A$',fontsize=FS)
plt.xlabel('Multiplier',fontsize=FS)
plt.grid(True, which="major", ls="-")
plt.grid(True, which="minor", ls="--")
plt.show()
fig.savefig('figure4d.png', format='png', dpi=300)