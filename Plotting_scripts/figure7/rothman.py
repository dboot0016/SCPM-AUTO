#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 11 13:38:30 2021

@author: daan
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 08:41:54 2021

@author: daan
"""

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
LW=3.5
FS=15
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

C1_2i1=np.interp(ti_2,t_2,C1_2)
C2_2i1=np.interp(ti_2,t_2,C2_2)
C3_2i1=np.interp(ti_2,t_2,C3_2)
C5_2i1=np.interp(ti_2,t_2,C5_2)
C6_2i1=np.interp(ti_2,t_2,C6_2)
C7_2i1=np.interp(ti_2,t_2,C7_2)

C1_2i=np.roll(C1_2i1,-jj)
C2_2i=np.roll(C2_2i1,-jj)
C3_2i=np.roll(C3_2i1,-jj)
C5_2i=np.roll(C5_2i1,-jj)
C6_2i=np.roll(C6_2i1,-jj)
C7_2i=np.roll(C7_2i1,-jj)

A1_2i1=np.interp(ti_2,t_2,A1_2)
A2_2i1=np.interp(ti_2,t_2,A2_2)
A3_2i1=np.interp(ti_2,t_2,A3_2)
A5_2i1=np.interp(ti_2,t_2,A5_2)
A6_2i1=np.interp(ti_2,t_2,A6_2)
A7_2i1=np.interp(ti_2,t_2,A7_2)

A1_2i=np.roll(A1_2i1,-jj)
A2_2i=np.roll(A2_2i1,-jj)
A3_2i=np.roll(A3_2i1,-jj)
A5_2i=np.roll(A5_2i1,-jj)
A6_2i=np.roll(A6_2i1,-jj)
A7_2i=np.roll(A7_2i1,-jj)

P1_2i1=np.interp(ti_2,t_2,P1_2)
P2_2i1=np.interp(ti_2,t_2,P2_2)
P3_2i1=np.interp(ti_2,t_2,P3_2)
P5_2i1=np.interp(ti_2,t_2,P5_2)
P6_2i1=np.interp(ti_2,t_2,P6_2)
P7_2i1=np.interp(ti_2,t_2,P7_2)

P1_2i=np.roll(P1_2i1,-jj)
P2_2i=np.roll(P2_2i1,-jj)
P3_2i=np.roll(P3_2i1,-jj)
P5_2i=np.roll(P5_2i1,-jj)
P6_2i=np.roll(P6_2i1,-jj)
P7_2i=np.roll(P7_2i1,-jj)

TC_2i1=np.interp(ti_2,t_2,TC_2)
TC_2i=np.roll(TC_2i1,-jj)
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
CO2=CO2_2i

C4=np.zeros(10000)
C1=C1_2i
C2=C2_2i
C3=C3_2i
C4=C4_2i
C5=C5_2i
C6=C6_2i
C7=C7_2i

Alk1=A1_2i
Alk2=A2_2i
Alk3=A3_2i
Alk4=A4_2i
Alk5=A5_2i
Alk6=A6_2i
Alk7=A7_2i

P1=P1_2i
P2=P2_2i
P3=P3_2i
P4=P4_2i
P5=P5_2i
P6=P6_2i
P7=P7_2i
# %%
psi2=np.zeros(10000)
RCP=np.zeros(10000)
RCP1=np.zeros(10000)
RCP2=np.zeros(10000)
RCP3=np.zeros(10000)
RCP4=np.zeros(10000)
RCP5=np.zeros(10000)
RCP6=np.zeros(10000)
RCP7=np.zeros(10000)
RPC1=np.zeros(10000)
RPC2=np.zeros(10000)
RPC3=np.zeros(10000)
RPC4=np.zeros(10000)
RPC5=np.zeros(10000)
RPC6=np.zeros(10000)
RPC7=np.zeros(10000)
Tadjust=np.zeros(10000)
T1=np.zeros(10000)
T2=np.zeros(10000)
T5=np.zeros(10000)
T7=np.zeros(10000)

K10=np.zeros(10000)
K01=np.zeros(10000)
K11=np.zeros(10000)
K21=np.zeros(10000)
K111=np.zeros(10000)
K1=np.zeros(10000)
dVC1=np.zeros(10000)
dKC1=np.zeros(10000)
Kspres1=np.zeros(10000)
Ksp1=np.zeros(10000)
eps1=np.zeros(10000)
PV1=np.zeros(10000)
Z1=np.zeros(10000)
DUM1=np.zeros(10000)
air1=np.zeros(10000)

K20=np.zeros(10000)
K02=np.zeros(10000)
K12=np.zeros(10000)
K22=np.zeros(10000)
K222=np.zeros(10000)
K2=np.zeros(10000)
dVC2=np.zeros(10000)
dKC2=np.zeros(10000)
Kspres2=np.zeros(10000)
Ksp2=np.zeros(10000)
eps2=np.zeros(10000)
PV2=np.zeros(10000)
Z2=np.zeros(10000)
DUM2=np.zeros(10000)
air2=np.zeros(10000)

K50=np.zeros(10000)
K05=np.zeros(10000)
K15=np.zeros(10000)
K25=np.zeros(10000)
K55=np.zeros(10000)
K5=np.zeros(10000)
dVC5=np.zeros(10000)
dKC5=np.zeros(10000)
Kspres5=np.zeros(10000)
Ksp5=np.zeros(10000)
eps5=np.zeros(10000)
Z5=np.zeros(10000)
DUM5=np.zeros(10000)
air5=np.zeros(10000)

K70=np.zeros(10000)
K07=np.zeros(10000)
K17=np.zeros(10000)
K27=np.zeros(10000)
K77=np.zeros(10000)
K7=np.zeros(10000)
dVC7=np.zeros(10000)
dKC7=np.zeros(10000)
Kspres7=np.zeros(10000)
Ksp7=np.zeros(10000)
eps7=np.zeros(10000)
PV7=np.zeros(10000)
Z7=np.zeros(10000)
DUM7=np.zeros(10000)
air7=np.zeros(10000)

RA1=np.zeros(10000)
DM1=np.zeros(10000)
H1=np.zeros(10000)
CO31=np.zeros(10000)
PhysC1=np.zeros(10000)
BioC1=np.zeros(10000)
CarbC1=np.zeros(10000)
PhysAlk1=np.zeros(10000)
BioAlk1=np.zeros(10000)
CarbAlk1=np.zeros(10000)
BioP1=np.zeros(10000)
PhysP1=np.zeros(10000)
RiverC=np.zeros(10000)
RiverAlk=np.zeros(10000)

RA2=np.zeros(10000)
DM2=np.zeros(10000)
H2=np.zeros(10000)
CO32=np.zeros(10000)
PhysC2=np.zeros(10000)
BioC2=np.zeros(10000)
CarbC2=np.zeros(10000)
PhysAlk2=np.zeros(10000)
BioAlk2=np.zeros(10000)
CarbAlk2=np.zeros(10000)
BioP2=np.zeros(10000)
PhysP2=np.zeros(10000)

RA3=np.zeros(10000)
DM3=np.zeros(10000)
H3=np.zeros(10000)
CO33=np.zeros(10000)
PhysC3=np.zeros(10000)
BioC3=np.zeros(10000)
CarbC3=np.zeros(10000)
PhysAlk3=np.zeros(10000)
BioAlk3=np.zeros(10000)
CarbAlk3=np.zeros(10000)
BioP3=np.zeros(10000)
PhysP3=np.zeros(10000)

RA4=np.zeros(10000)
DM4=np.zeros(10000)
H4=np.zeros(10000)
CO34=np.zeros(10000)
PhysC4=np.zeros(10000)
BioC41=np.zeros(10000)
BioC4=np.zeros(10000)
CarbC4=np.zeros(10000)
PhysAlk4=np.zeros(10000)
BioAlk4=np.zeros(10000)
CarbAlk4=np.zeros(10000)
BioP4=np.zeros(10000)
PhysP4=np.zeros(10000)

RA5=np.zeros(10000)
DM5=np.zeros(10000)
H5=np.zeros(10000)
CO35=np.zeros(10000)
PhysC5=np.zeros(10000)
BioC5=np.zeros(10000)
CarbC5=np.zeros(10000)
PhysAlk5=np.zeros(10000)
BioAlk5=np.zeros(10000)
CarbAlk5=np.zeros(10000)
BioP5=np.zeros(10000)
PhysP5=np.zeros(10000)
PV5=np.zeros(10000)

RA6=np.zeros(10000)
DM6=np.zeros(10000)
H6=np.zeros(10000)
CO36=np.zeros(10000)
PhysC6=np.zeros(10000)
BioC6=np.zeros(10000)
CarbC6=np.zeros(10000)
CarbC61=np.zeros(10000)
CarbC62=np.zeros(10000)
PhysAlk6=np.zeros(10000)
BioAlk6=np.zeros(10000)
CarbAlk6=np.zeros(10000)
BioP6=np.zeros(10000)
PhysP6=np.zeros(10000)

RA7=np.zeros(10000)
DM7=np.zeros(10000)
H7=np.zeros(10000)
CO37=np.zeros(10000)
PhysC7=np.zeros(10000)
BioC7=np.zeros(10000)
CarbC7=np.zeros(10000)
PhysAlk7=np.zeros(10000)
BioAlk7=np.zeros(10000)
CarbAlk7=np.zeros(10000)
BioP7=np.zeros(10000)
PhysP7=np.zeros(10000)

CC=np.zeros(10000)

HCO31=np.zeros(10000)
HCO32=np.zeros(10000)
HCO33=np.zeros(10000)
HCO34=np.zeros(10000)
HCO35=np.zeros(10000)
HCO36=np.zeros(10000)
HCO37=np.zeros(10000)

H2CO31=np.zeros(10000)
H2CO32=np.zeros(10000)
H2CO33=np.zeros(10000)
H2CO34=np.zeros(10000)
H2CO35=np.zeros(10000)
H2CO36=np.zeros(10000)
H2CO37=np.zeros(10000)

# %%
for j in range(10000):
    AMOC=0
    RCP_VAR=0
    PV_VAR=1
    FCA_VAR=0
    FCA_TEMP=0
    RCP_CO2=0
    eps_TEMP=1.5
    TEMP=2.085
    BIO_ALK=0
    BIO=1
            
    Vat=1.7604418363824646*(10**20)      # Volume atmosphere
    rho=1029       # Sea water density
    FCA=0.070025243216#0.07      # Rain ratio
    n=1.0           # Order dissolution kinetics
    PerC=0.12     # Percentage of C in CaCO3 or something
    DC=2.75*(10**(-13))      # Constant dissolution
    RSiP=16.      # Redfield Si to P ration
    WSC=2.378234398782344*(10**-12)     # Constant silicate weathering
    WSV=1.585489599188229*(10**-8)     # Variable silicate weathering
    WCV=6.341958396752917*(10**-8)    # Variable CaCO3 weathering
    kCa=4.398148148148148*(10**-6)     # ?
    Pin=15356.44572697869    # Phosphate input via rivers
    b=0.75       # Martin np.exponent
    d0=100      # Base depth for Martin profile
    fract=0.5    # Fraction of GOC going into box 7 (from box 4)
    yr=31536000.      # Seconds per year
    day=86400.     # Seconds per day
    CO2baseHol=145.53186722*(10**-6) # CO2 concentration Holocene
    gamma1=31*(10**6)  # Upwelling box 4 and 6
    gamma2=40*(10**6)  # Upwelling box 1 and 3
    psi1=1.8*10**7#18*(10**6)    # GOC
    psi2_base=15*10**6
    PV=3.      # Piston velocity
    Ant= 0.
    eta=0.81
    a_RCP=3/70
    b_RCP=117
    RCP_base=130
    
    RCP[j]=RCP_CO2*(a_RCP*CO2[j]*10**6+b_RCP)+(1-RCP_CO2)*RCP_base        # Redfield C to P ratio
    RCPD=331        # Redfield C to P ratio Diazotrophs (CSIRO Mk3L COAL)
    
    RCP1[j]=RCP_VAR*(0.9*(1/(6.9/1000*P1[j]*1000+6/1000))+0.1*RCPD)+(1-RCP_VAR)*RCP[j]
    RCP2[j]=RCP_VAR*(0.9*1/(6.9/1000*P2[j]*1000+6/1000)+0.1*RCPD)+(1-RCP_VAR)*RCP[j]
    RCP3[j]=RCP_VAR*(0.9*1/(6.9/1000*P3[j]*1000+6/1000)+0.1*RCPD)+(1-RCP_VAR)*RCP[j]
    RCP4[j]=RCP_VAR*(0.9*1/(6.9/1000*P4[j]*1000+6/1000)+0.1*RCPD)+(1-RCP_VAR)*RCP[j]
    RCP5[j]=RCP_VAR*(0.9*1/(6.9/1000*P5[j]*1000+6/1000)+0.1*RCPD)+(1-RCP_VAR)*RCP[j]
    RCP6[j]=RCP_VAR*(0.9*1/(6.9/1000*P6[j]*1000+6/1000)+0.1*RCPD)+(1-RCP_VAR)*RCP[j]
    RCP7[j]=RCP_VAR*(0.9*1/(6.9/1000*P7[j]*1000+6/1000)+0.1*RCPD)+(1-RCP_VAR)*RCP[j]
    
    RPC1[j]=1/RCP1[j]      # Redfield P to C ratio
    RPC2[j]=1/RCP2[j]      # Redfield P to C ratio
    RPC3[j]=1/RCP3[j]      # Redfield P to C ratio
    RPC4[j]=1/RCP4[j]      # Redfield P to C ratio
    RPC5[j]=1/RCP5[j]      # Redfield P to C ratio
    RPC6[j]=1/RCP6[j]      # Redfield P to C ratio
    RPC7[j]=1/RCP7[j]      # Redfield P to C ratio
    
    psi2[j]=15000000+0*CO2[j]#AMOC1*10**6#psi2_base*(1-AMOC*0.1*5.35*0.54*np.log((CO2[j])/CO2baseHol)) # AMOC strength
    Tadjust[j]=TEMP*0.54*5.35*np.log((CO2[j])/CO2baseHol) # Adjustment temperature CO2 concentration
    
    # For K0
    A1_C=-60.3409
    A2_C=93.4517
    A3_C=23.3585
    B1_C=0.023517
    B2_C=-0.023656
    B3_C=0.0047036
    
    # Box 1
    T1[j]=17.34+Tadjust[j]      # Temperature
    Sal1=36.25   # Salinity
    BT1=1.179e-5*Sal1 # Total Boron
    Ca1=0.01028*Sal1/35 # Calcium
    K10[j]=A1_C+A2_C*(100.0/((T1[j]+273.16)))+A3_C*np.log((T1[j]+273.16)/100.0)
    K01[j]=np.exp((K10[j]+(Sal1)*(B1_C+B2_C*((T1[j]+273.16)/100.0)+B3_C*(((T1[j]+273.16)/100.0)**2)))) # Weiss (1974) in mols kg-1 atm-1
    K11[j]=10**-((3633.86/(T1[j]+273.16))-61.2172+(9.67770*np.log(T1[j]+273.16))-(0.011555*Sal1)+(0.0001152*(Sal1**2))) #Lueker et al (2000) in mol kg-1
    K21[j]=10**-((471.78/((T1[j]+273.16)))+25.9290-(3.16967*np.log((T1[j]+273.16)))-(0.01781*Sal1)+(0.0001122*(Sal1**2)))# Lueker et al (2000) in mol kg-1
    K111[j]=-395.8293+(6537.773/(T1[j]+273.16))+71.595*np.log((T1[j]+273.16))-0.17959*(T1[j]+273.16)
    K1[j]=np.exp(K111[j]+(-1.78938+410.64/(T1[j]+273.16)+0.0065453*(T1[j]+273.16))*Sal1**0.5-0.17755*Sal1+0.0094979*Sal1**(3.0/2.0)) # Mucci (1983)
    dVC1[j]=-65.28+0.397*T1[j]-0.005155*T1[j]**2+(19.816-0.0441*T1[j]-0.00017*T1[j]**2)*(Sal1/35)**0.5
    dKC1[j]=0.01847+0.0001956*T1[j]-0.000002212*T1[j]**2+(-0.03217-0.0000711*T1[j]-0.000002212*T1[j]**2)*(Sal1/35)**0.5
    Kspres1[j]=(-dVC1[j]+0.5*dKC1[j]*5)*5/(83.144621*(T1[j]+273.16))
    Ksp1[j]=K1[j]*np.exp(Kspres1[j])
    
    Zo1=1.1      # Bionp.logical production
    S1=2.6328225e+14    # Surface
    V1=2.6328225e+16      # Volume
    df1=100.      # Floor depth
    eps1_base=0.9     # Bionp.logical efficiency
    
    # Box 2
    T2[j]=3.1+Tadjust[j]       # Temperature
    Sal2=35.27    # Salinity
    BT2=1.179e-5*Sal2 # Total Boron
    Ca2=0.01028*Sal2/35 # Calcium
    K20[j]=A1_C+A2_C*(100.0/((T2[j]+273.16)))+A3_C*np.log((T2[j]+273.16)/100.0)
    K02[j]=np.exp((K20[j]+(Sal2)*(B1_C+B2_C*((T2[j]+273.16)/100.0)+B3_C*(((T2[j]+273.16)/100.0)**2)))) # Weiss (1974) in mols kg-1 atm-1
    
    K12[j] = 10**-((3633.86/(T2[j]+273.16)) - 61.2172 + (9.67770 * np.log(T2[j]+273.16)) - (0.011555 * Sal2) + (0.0001152 * (Sal2**2)));
    K22[j]=10**-((471.78/((T2[j]+273.16)))+25.9290-(3.16967*np.log((T2[j]+273.16)))-(0.01781*Sal2)+(0.0001122*(Sal2**2)))# Lueker et al (2000) in mol kg-1
    K222[j]=-395.8293+(6537.773/(T2[j]+273.16))+71.595*np.log((T2[j]+273.16))-0.17959*(T2[j]+273.16)
    K2[j]=np.exp(K222[j]+(-1.78938+410.64/(T2[j]+273.16)+0.0065453*(T2[j]+273.16))*Sal2**0.5-0.17755*Sal2+0.0094979*Sal2**(3.0/2.0)) # Mucci (1983)
    dVC2[j]=-65.28+0.397*T2[j]-0.005155*T2[j]**2+(19.816-0.0441*T2[j]-0.00017*T2[j]**2)*(Sal2/35)**0.5
    dKC2[j]=0.01847+0.0001956*T2[j]-0.000002212*T2[j]**2+(-0.03217-0.0000711*T2[j]-0.000002212*T2[j]**2)*(Sal2/35)**0.5
    Kspres2[j]=(-dVC2[j]+0.5*dKC2[j]*0)*0/(83.144621*(T2[j]+273.16))
    Ksp2[j]=K2[j]*np.exp(Kspres2[j])
    
    Zo2=4.5      # Bionp.logical production
    S2=3.5104300e+13      # Surface
    V2=8.77607500e+15     # Volume
    df2=250.    # Floor depth
    eps2_base=1.25    # Bionp.logical efficiency
    
    # Box 3
    T3=11.28     # Temperature
    Sal3=34.91    # Salinity
    BT3=1.179e-5*Sal3 # Total Boron
    Ca3=0.01028*Sal3/35 # Calcium
    K30=A1_C+A2_C*(100.0/((T3+273.16)))+A3_C*np.log((T3+273.16)/100.0)
    K03=np.exp((K30+(Sal3)*(B1_C+B2_C*((T3+273.16)/100.0)+B3_C*(((T3+273.16)/100.0)**2)))) # Weiss (1974) in mols kg-1 atm-1
    K13=10**-((3633.86/(T3+273.16)) - 61.2172 + (9.67770 * np.log(T3+273.16)) - (0.011555 * Sal3) + (0.0001152 * (Sal3**2))) #Lueker et al (2000) in mol kg-1
    K23=10**-((471.78/((T3+273.16))) + 25.9290 - (3.16967 * np.log((T3+273.16))) - (0.01781 * Sal3) + (0.0001122 * (Sal3**2)))# Lueker et al (2000) in mol kg-1
    K33=-395.8293+(6537.773/(T3+273.16))+71.595*np.log((T3+273.16))-0.17959*(T3+273.16)
    K3=np.exp(K33+(-1.78938+410.64/(T3+273.16)+0.0065453*(T3+273.16))*Sal3**0.5-0.17755*Sal3+0.0094979*Sal3**(3.0/2.0)) # Mucci (1983)
    dVC3=-65.28+0.397*T3-0.005155*T3**2+(19.816-0.0441*T3-0.00017*T3**2)*(Sal3/35)**0.5
    dKC3=0.01847+0.0001956*T3-0.000002212*T3**2+(-0.03217-0.0000711*T3-0.000002212*T3**2)*(Sal3/35)**0.5
    Kspres3=(-dVC3+0.5*dKC3*55)*55/(83.144621*(T3+273.16))
    Ksp3=K3*np.exp(Kspres3)
    
    S3=2.6328225e+14     # Surface
    V3=2.36954025e+17      # Volume
    df3=1000.      # Floor depth
    
    # Box 4
    T4=3.24      # Temperature
    Sal4=34.76   # Salinity
    BT4=1.179e-5*Sal4 # Total Boron
    Ca4=0.01028*Sal4/35 # Calcium
    K40=A1_C+A2_C*(100.0/((T4+273.16)))+A3_C*np.log((T4+273.16)/100.0)
    K04=np.exp((K40+(Sal4)*(B1_C+B2_C*((T4+273.16)/100.0)+B3_C*(((T4+273.16)/100.0)**2)))) # Weiss (1974) in mols kg-1 atm-1
    K14=10**-((3633.86/(T4+273.16))-61.2172+(9.67770*np.log(T4+273.16))-(0.011555*Sal4)+(0.0001152*(Sal4**2))) #Lueker et al (2000) in mol kg-1
    K24=10**-((471.78/((T4+273.16)))+25.9290-(3.16967*np.log((T4+273.16)))-(0.01781*Sal4)+(0.0001122*(Sal4**2)))# Lueker et al (2000) in mol kg-1
    K44=-395.8293+(6537.773/(T4+273.16))+71.595*np.log((T4+273.16))-0.17959*(T4+273.16)
    K4=np.exp(K44+(-1.78938+410.64/(T4+273.16)+0.0065453*(T4+273.16))*Sal4**0.5-0.17755*Sal4+0.0094979*Sal4**(3.0/2.0)) # Mucci (1983)
    dVC4=-65.28+0.397*T4-0.005155*T4**2+(19.816-0.0441*T4-0.00017*T4**2)*(Sal4/35)**0.5
    dKC4=0.01847+0.0001956*T4-0.000002212*T4**2+(-0.03217-0.0000711*T4-0.000002212*T4**2)*(Sal4/35)**0.5
    Kspres4=(-dVC4+0.5*dKC4*152.632)*152.632/(83.144621*(T4+273.16))
    Ksp4=K4*np.exp(Kspres4)
      
    S4=3.3349085e+14      # Surface 3.4434785e+14
    V4=5.52892725e+17     # Volume
    df4=2500.      # Floor depth
    
    # Box 5
    T5[j]=0.93+Tadjust[j]/TEMP      # Temperature
    Sal5=35.43   # Salinity
    BT5=1.179e-5*Sal5 # Total Boron
    Ca5=0.01028*Sal5/35 # Calcium
    K50[j]=A1_C+A2_C*(100.0/((T5[j]+273.16)))+A3_C*np.log((T5[j]+273.16)/100.0)
    K05[j]=np.exp((K50[j] + (Sal5) * (B1_C + B2_C *((T5[j]+273.16)/100.0) + B3_C * (((T5[j]+273.16)/100.0)**2)))) # Weiss (1974) in mols kg-1 atm-1
    K15[j]=10**-((3633.86/(T5[j]+273.16))-61.2172+(9.67770*np.log(T5[j]+273.16))-(0.011555*Sal5)+(0.0001152*(Sal5**2))) #Lueker et al (2000) in mol kg-1
    K25[j]=10**-((471.78/((T5[j]+273.16)))+25.9290-(3.16967*np.log((T5[j]+273.16)))-(0.01781*Sal5)+(0.0001122*(Sal5**2)))# Lueker et al (2000) in mol kg-1
    K55[j]=-395.8293+(6537.773/(T5[j]+273.16))+71.595*np.log((T5[j]+273.16))-0.17959*(T5[j]+273.16)
    K5[j]=np.exp(K55[j]+(-1.78938+410.64/(T5[j]+273.16)+0.0065453*(T5[j]+273.16))*Sal5**0.5-0.17755*Sal5+0.0094979*Sal5**(3.0/2.0)) # Mucci (1983)
    dVC5[j]=-65.28+0.397*T5[j]-0.005155*T5[j]**2+(19.816-0.0441*T5[j]-0.00017*T5[j]**2)*(Sal5/35)**0.5
    dKC5[j]=0.01847+0.0001956*T5[j]-0.000002212*T5[j]**2+(-0.03217-0.0000711*T5[j]-0.000002212*T5[j]**2)*(Sal5/35)**0.5
    Kspres5[j]=(-dVC5[j]+0.5*dKC5[j]*0)*0/(83.144621*(T5[j]+273.16))
    Ksp5[j]=K5[j]*np.exp(Kspres5[j])
       
    Zo5=1.75      # Bionp.logical production
    S5=1.7552150e+13# Surface
    V5=4.38803750e+16      # Volume
    df5=2500.     # Floor depth
    eps5_base=0.35    # Bionp.logical efficiency
    
    # Box 6
    T6=1.8      # Temperature
    Sal6=34.77    # Salinity
    BT6=1.179e-5*Sal6 # Total Boron
    Ca6=0.01028*Sal6/35 # Calcium
    K60=A1_C+A2_C*(100.0/((T6+273.16)))+A3_C*np.log((T6+273.16)/100.0)
    K06=np.exp((K60+(Sal6)*(B1_C+B2_C*((T6+273.16)/100.0)+B3_C*(((T6+273.16)/100.0)**2)))) # Weiss (1974) in mols kg-1 atm-1
    K16=10**-((3633.86/(T6+273.16))-61.2172+(9.67770*np.log(T6+273.16))-(0.011555*Sal6)+(0.0001152*(Sal6**2))) #Lueker et al (2000) in mol kg-1
    K26=10**-((471.78/((T6+273.16)))+25.9290-(3.16967*np.log((T6+273.16)))-(0.01781*Sal6)+(0.0001122*(Sal6**2)))# Lueker et al (2000) in mol kg-1
    K66=-395.8293+(6537.773/(T6+273.16))+71.595*np.log((T6+273.16))-0.17959*(T6+273.16)
    K6=np.exp(K66+(-1.78938+410.64/(T6+273.16)+0.0065453*(T6+273.16))*Sal6**0.5-0.17755*Sal6+0.0094979*Sal6**(3.0/2.0)) # Mucci (1983)
    dVC6=-65.28+0.397*T6-0.005155*T6**2+(19.816-0.0441*T6-0.00017*T6**2)*(Sal6/35)**0.5
    dKC6=0.01847+0.0001956*T6-0.000002212*T6**2+(-0.03217-0.0000711*T6-0.000002212*T6**2)*(Sal6/35)**0.5
    Kspres6=(-dVC6+0.5*dKC6*276.25)*276.25/(83.144621*(T6+273.16))
    Ksp6=K6*np.exp(Kspres6)
    KspSed=K6*np.exp((-dVC6+0.5*dKC6*400)*400/(83.144621*(T6+273.16)))
      
    S6=3.5104300e+14      # Surface
    V6=5.26564500e+17      # Volume
    df6=4000.    # Floor depth
    
    # Box 7
    T7[j]=0.33+Tadjust[j]      # Temperature
    Sal7=35.17   # Salinity
    BT7=1.179e-5*Sal7 # Total Boron
    Ca7=0.01028*Sal7/35 # Calcium
    K70[j]=A1_C+A2_C*(100.0/((T7[j]+273.16)))+A3_C*np.log((T7[j]+273.16)/100.0)
    K07[j]=np.exp((K70[j]+(Sal7)*(B1_C+B2_C*((T7[j]+273.16)/100.0)+B3_C*(((T7[j]+273.16)/100.0)**2)))) # Weiss (1974) in mols kg-1 atm-1
    K17[j]=10**-((3633.86/(T7[j]+273.16))-61.2172+(9.67770*np.log(T7[j]+273.16))-(0.011555*Sal7)+(0.0001152*(Sal7**2))) #Lueker et al (2000) in mol kg-1
    K27[j]=10**-((471.78/((T7[j]+273.16)))+25.9290-(3.16967*np.log((T7[j]+273.16)))-(0.01781*Sal7)+(0.0001122*(Sal7**2)))# Lueker et al (2000) in mol kg-1   
    K77[j]=-395.8293+(6537.773/(T7[j]+273.16))+71.595*np.log((T7[j]+273.16))-0.17959*(T7[j]+273.16)
    K7[j]=np.exp(K77[j]+(-1.78938+410.64/(T7[j]+273.16)+0.0065453*(T7[j]+273.16))*Sal7**0.5-0.17755*Sal7+0.0094979*Sal7**(3.0/2.0)) # Mucci (1983)
    dVC7[j]=-65.28+0.397*T7[j]-0.005155*T7[j]**2+(19.816-0.0441*T7[j]-0.00017*T7[j]**2)*(Sal7/35)**0.5
    dKC7[j]=0.01847+0.0001956*T7[j]-0.000002212*T7[j]**2+(-0.03217-0.0000711*T7[j]-0.000002212*T7[j]**2)*(Sal7/35)**0.5
    Kspres7[j]=(-dVC7[j]+0.5*dKC7[j]*12.5)*12.5/(83.144621*(T7[j]+273.16))
    Ksp7[j]=K7[j]*np.exp(Kspres7[j])
    
    Zo7=5.325      # Bionp.logical production
    S7=3.5104300e+13  # Surface
    V7=8.77607500e+15      # Volume
    df7=250.     # Floor depth
    eps7_base=0.62    # Bionp.logical efficiency
    
    # Bionp.logical efficiency feedback
    eps1[j]=eps_TEMP*Tadjust[j]*(-0.1)+eps1_base
    eps2[j]=eps_TEMP*Tadjust[j]*(-0.1)+eps2_base
    eps5[j]=eps_TEMP*Tadjust[j]*(-0.1)/TEMP+eps5_base
    eps7[j]=eps_TEMP*Tadjust[j]*(-0.1)+eps7_base
    
       # Piston velocity
    PV1[j]=PV*PV_VAR*(((2116.8-136.25*T1[j]+4.7353*T1[j]**2-0.092307*T1[j]**3+0.0007555*T1[j]**4)/660)**-0.5)+PV*(1-PV_VAR)
    PV2[j]=PV*PV_VAR*(((2116.8-136.25*T2[j]+4.7353*T2[j]**2-0.092307*T2[j]**3+0.0007555*T2[j]**4)/660)**-0.5)+PV*(1-PV_VAR)
    PV5[j]=(PV-2)*PV_VAR*(((2116.8-136.25*T5[j]+4.7353*T5[j]**2-0.092307*T5[j]**3+0.0007555*T5[j]**4)/660)**-0.5)+(PV-2)*(1-PV_VAR)
    PV7[j]=PV*PV_VAR*(((2116.8-136.25*T7[j]+4.7353*T7[j]**2-0.092307*T7[j]**3+0.0007555*T7[j]**4)/660)**-0.5)+PV*(1-PV_VAR)

    
       # Carbonate chemistry (Williams and Follows, 2011)
       # Def of RA_i only needed to start from trivial solution 
    RA1[j] = (C1[j]/rho)/(Alk1[j]/rho)
    RA2[j] = (C2[j]/rho)/(Alk2[j]/rho)
    RA3[j] = (C3[j]/rho)/(Alk3[j]/rho)
    RA4[j] = (C4[j]/rho)/(Alk4[j]/rho)
    RA5[j] = (C5[j]/rho)/(Alk5[j]/rho)
    RA6[j] = (C6[j]/rho)/(Alk6[j]/rho)
    RA7[j] = (C7[j]/rho)/(Alk7[j]/rho)
    
    DM1[j]=(1-RA1[j])**2*K11[j]**2-4*K11[j]*K21[j]*(1-2*RA1[j])
    DM2[j]=(1-RA2[j])**2*K12[j]**2-4*K12[j]*K22[j]*(1-2*RA2[j])
    DM3[j]=(1-RA3[j])**2*K13**2-4*K13*K23*(1-2*RA3[j])
    DM4[j]=(1-RA4[j])**2*K14**2-4*K14*K24*(1-2*RA4[j])
    DM5[j]=(1-RA5[j])**2*K15[j]**2-4*K15[j]*K25[j]*(1-2*RA5[j])
    DM6[j]=(1-RA6[j])**2*K16**2-4*K16*K26*(1-2*RA6[j])
    DM7[j]=(1-RA7[j])**2*K17[j]**2-4*K17[j]*K27[j]*(1-2*RA7[j])
    
    H1[j]=0.5*((RA1[j]-1)*K11[j]+DM1[j]**0.5)
    H2[j]=0.5*((RA2[j]-1)*K12[j]+DM2[j]**0.5)
    H3[j]=0.5*((RA3[j]-1)*K13+DM3[j]**0.5)
    H4[j]=0.5*((RA4[j]-1)*K14+DM4[j]**0.5)
    H5[j]=0.5*((RA5[j]-1)*K15[j]+DM5[j]**0.5)
    H6[j]=0.5*((RA6[j]-1)*K16+DM6[j]**0.5)
    H7[j]=0.5*((RA7[j]-1)*K17[j]+DM7[j]**0.5)
    
    CO31[j]=(C1[j]/rho/K01[j])*(H1[j]**2/(H1[j]**2+H1[j]*K11[j]+K11[j]*K21[j]))*K01[j]*K11[j]/H1[j]*K21[j]/H1[j]
    CO32[j]=(C2[j]/rho/K02[j])*(H2[j]**2/(H2[j]**2+H2[j]*K12[j]+K12[j]*K22[j]))*K02[j]*K12[j]/H2[j]*K22[j]/H2[j]
    CO33[j]=(C3[j]/rho/K03)*(H3[j]**2/(H3[j]**2+H3[j]*K13+K13*K23))*K03*K13/H3[j]*K23/H3[j]
    CO34[j]=(C4[j]/rho/K04)*(H4[j]**2/(H4[j]**2+H4[j]*K14+K14*K24))*K04*K14/H4[j]*K24/H4[j]
    CO35[j]=(C5[j]/rho/K05[j])*(H5[j]**2/(H5[j]**2+H5[j]*K15[j]+K15[j]*K25[j]))*K05[j]*K15[j]/H5[j]*K25[j]/H5[j]
    CO36[j]=(C6[j]/rho/K06)*(H6[j]**2/(H6[j]**2+H6[j]*K16+K16*K26))*K06*K16/H6[j]*K26/H6[j]
    CO37[j]=(C7[j]/rho/K07[j])*(H7[j]**2/(H7[j]**2+H7[j]*K17[j]+K17[j]*K27[j]))*K07[j]*K17[j]/H7[j]*K27[j]/H7[j]
    
    HCO31[j]=CO31[j]*H1[j]/K21[j]
    HCO32[j]=CO32[j]*H2[j]/K22[j]
    HCO33[j]=CO33[j]*H3[j]/K23
    HCO34[j]=CO34[j]*H4[j]/K24
    HCO35[j]=CO35[j]*H5[j]/K25[j]
    HCO36[j]=CO36[j]*H6[j]/K26
    HCO37[j]=CO37[j]*H7[j]/K27[j]
    
    H2CO31[j]=HCO31[j]*H1[j]/K11[j]
    H2CO32[j]=HCO32[j]*H2[j]/K12[j]
    H2CO33[j]=HCO33[j]*H3[j]/K13
    H2CO34[j]=HCO34[j]*H4[j]/K14
    H2CO35[j]=HCO35[j]*H5[j]/K15[j]
    H2CO36[j]=HCO36[j]*H6[j]/K16
    H2CO37[j]=HCO37[j]*H7[j]/K17[j]
    
    # Rain ratio
    FCA1=FCA
    FCA2=FCA
    FCA5=FCA
    FCA7=FCA
    
    FCA5=FCA_TEMP*FCA5*(T5[j]+2)/7+(1-FCA_TEMP)*FCA5
    
    # New bionp.logy
    Z1[j]=(1-BIO)*Zo1+BIO*(gamma2*P3[j]*yr*1+Pin*yr)*eps1[j]*RCP1[j]*1/S1
    Z2[j]=(1-BIO)*Zo2+BIO*psi2[j]*P3[j]*yr*eps2[j]*RCP2[j]*1/S2
    Z5[j]=(1-BIO)*Zo5+BIO*fract*psi1*P7[j]*yr*eps5[j]*RCP5[j]*1/S5
    Z7[j]=(1-BIO)*Zo7+BIO*(psi2[j]+fract*psi1)*P4[j]*yr*eps7[j]*RCP7[j]*1/S7
    
    # Define the different components (physical, bionp.logical, air exchange, carbonate, and river)
    # Box 1
    PhysC1[j]=gamma2*(C3[j]-C1[j])/V1
    BioC1[j]=-(Z1[j]*S1/(V1*yr))*(df1/d0)**(-b)
    DUM1[j] = (C1[j]/rho/K01[j])*((H1[j]**2)/((H1[j]**2)+H1[j]*K11[j]+K11[j]*K21[j]))
    air1[j]=(rho*K01[j]*PV1[j]/day*S1*(CO2[j]-DUM1[j]))/V1
    CarbC1[j]=-(Z1[j]*S1/(V1*yr))*FCA1+(CO31[j]*Ca1*rho*kCa*(1-min((CO31[j]*Ca1/Ksp1[j]),n))**n*PerC+DC)
    RiverC[j]=WSC+(WSV+WCV)*CO2[j]
    
    PhysAlk1[j]=gamma2*(Alk3[j]-Alk1[j])/V1
    CarbAlk1[j]=2*CarbC1[j]
    RiverAlk[j]=2*RiverC[j]
    BioAlk1[j]=BIO_ALK*(-BioC1[j])*16/RCP1[j]
    
    PhysP1[j]=gamma2*(P3[j]-P1[j])/V1
    BioP1[j]=BioC1[j]*RPC1[j]
    RiverP=Pin/V1
    
    # Box 2
    PhysC2[j]=psi2[j]*(C3[j]-C2[j])/V2
    BioC2[j]=-(Z2[j]*S2/(V2*yr))*(df2/d0)**(-b)
    DUM2[j] = (C2[j]/rho/K02[j])*((H2[j]**2)/((H2[j]**2)+H2[j]*K12[j]+K12[j]*K22[j]))
    air2[j]=(rho*K02[j]*PV2[j]/day*S2*(CO2[j]-DUM2[j]))/V2
    CarbC2[j]=-(Z2[j]*S2/(V2*yr))*FCA2+(CO32[j]*Ca2*rho*kCa*(1-min((CO32[j]*Ca2/Ksp2[j]),n))**n*PerC+DC)
    
    PhysAlk2[j]=psi2[j]*(Alk3[j]-Alk2[j])/V2
    CarbAlk2[j]=2*CarbC2[j]
    BioAlk2[j]=BIO_ALK*(-BioC2[j])*16/RCP2[j]
    
    PhysP2[j]=psi2[j]*(P3[j]-P2[j])/V2
    BioP2[j]=BioC2[j]*RPC2[j]
    
    # Box 3
    PhysC3[j]=gamma2*(C1[j]-C3[j])/V3+psi2[j]*(C7[j]-C3[j])/V3
    BioC3[j]=Z1[j]*S1/(V3*yr)*((df1/d0)**(-b)-(df3/d0)**(-b))
    CarbC3[j]=(CO33[j]*Ca3*rho*kCa*(1-min((CO33[j]*Ca3/Ksp3),n))**n*PerC+DC)
    
    PhysAlk3[j]=gamma2*(Alk1[j]-Alk3[j])/V3+psi2[j]*(Alk7[j]-Alk3[j])/V3
    CarbAlk3[j]=2*CarbC3[j]
    BioAlk3[j]=BIO_ALK*(-BioC3[j])*16/RCP3[j]
    
    PhysP3[j]=gamma2*(P1[j]-P3[j])/V3+psi2[j]*(P7[j]-P3[j])/V3
    BioP3[j]=BioC3[j]*RPC3[j]
    
    # Box 4
    PhysC4[j]=gamma1*(C6[j]-C4[j])/V4+psi1*(C6[j]-C4[j])/V4+psi2[j]*(C2[j]-C4[j])/V4
    BioC41[j]=(Z1[j]*S1*((df3/d0)**(-b)-(df4/d0)**(-b))+Z2[j]*S2*((df2/d0)**(-b)-(df4/d0)**(-b))+Z7[j]*S7*((df7/d0)**(-b)-(df4/d0)**(-b)))
    BioC4[j]=BioC41[j]/(V4*yr)
    CarbC4[j]=(CO34[j]*Ca4*rho*kCa*(1-min((CO34[j]*Ca4/Ksp4),n))**n*PerC+DC)
    
    PhysAlk4[j]=gamma1*(Alk6[j]-Alk4[j])/V4+psi1*(Alk6[j]-Alk4[j])/V4+psi2[j]*(Alk2[j]-Alk4[j])/V4
    CarbAlk4[j]=2*CarbC4[j]
    BioAlk4[j]=BIO_ALK*(-BioC4[j])*16/RCP4[j]
    
    PhysP4[j]=gamma1*(P6[j]-P4[j])/V4+psi1*(P6[j]-P4[j])/V4+psi2[j]*(P2[j]-P4[j])/V4
    BioP4[j]=BioC4[j]*RPC4[j]
    
    # Box 5
    PhysC5[j]=(1-fract)*psi1*(C4[j]-C5[j])/V5+fract*psi1*(C7[j]-C5[j])/V5
    BioC5[j]=-(Z5[j]*S5/(V5*yr))*(df5/d0)**(-b)
    DUM5[j] =(C5[j]/rho/K05[j])*((H5[j]**2)/((H5[j]**2)+H5[j]*K15[j]+K15[j]*K25[j]))
    air5[j]=(rho*K05[j]*PV5[j]/day*S5*(CO2[j]-DUM5[j]))/V5
    CarbC5[j]=-(Z5[j]*S5/(V5*yr))*FCA5+(CO35[j]*Ca5*rho*kCa*(1-min((CO35[j]*Ca5/Ksp5[j]),n))**n*PerC+DC)
    
    PhysAlk5[j]=(1-fract)*psi1*(Alk4[j]-Alk5[j])/V5+fract*psi1*(Alk7[j]-Alk5[j])/V5
    CarbAlk5[j]=2*CarbC5[j]
    BioAlk5[j]=BIO_ALK*(-BioC5[j])*16/RCP5[j]
    
    PhysP5[j]=(1-fract)*psi1*(P4[j]-P5[j])/V5+fract*psi1*(P7[j]-P5[j])/V5
    BioP5[j]=BioC5[j]*RPC5[j]
    
    # Box 6
    PhysC6[j]=psi1*(C5[j]-C6[j])/V6+gamma1*(C4[j]-C6[j])/V6
    BioC6[j]=(Z1[j]*S1+Z2[j]*S2+Z5[j]*S5+Z7[j]*S7)/(V6*yr)*((df4/d0)**(-b))#-(4000/d0)**(-b))
    #CarbC6=(CO36*Ca6*rho*kCa*((1-min((CO36*Ca6/Ksp6),n))**n+(1-min((CO36*Ca6/KspSed),n))**n)*PerC+DC)
    CarbC61[j]=(CO36[j]*Ca6*rho*kCa*((1-min((CO36[j]*Ca6/Ksp6),n))**n))*PerC+DC;
    CarbC62[j]=(CO36[j]*Ca6*rho*kCa*((1-min((CO36[j]*Ca6/KspSed),n))**n))*PerC+DC;   
    CarbC6[j]=CarbC61[j]+CarbC62[j];
    
    PhysAlk6[j]=psi1*(Alk5[j]-Alk6[j])/V6+gamma1*(Alk4[j]-Alk6[j])/V6
    CarbAlk6[j]=2*CarbC6[j]
    BioAlk6[j]=BIO_ALK*(-BioC6[j])*16/RCP6[j]
    
    PhysP6[j]=psi1*(P5[j]-P6[j])/V6+gamma1*(P4[j]-P6[j])/V6
    BioP6[j]=BioC6[j]*RPC6[j]
    SedP6=-Pin/V6
    
    # Box 7
    PhysC7[j]=(psi2[j]+fract*psi1)*(C4[j]-C7[j])/V7
    BioC7[j]=-(Z7[j]*S7)/(V7*yr)*(df7/d0)**(-b)
    DUM7[j] = (C7[j]/rho/K07[j])*((H7[j]**2)/((H7[j]**2)+H7[j]*K17[j]+K17[j]*K27[j]))
    air7[j]=(rho*K07[j]*PV7[j]/day*S7*(CO2[j]-DUM7[j]))/V7    
    CarbC7[j]=-(Z7[j]*S7/(V7*yr))*FCA7+(CO37[j]*Ca7*rho*kCa*(1-min((CO37[j]*Ca7/Ksp7[j]),n))**n*PerC+DC)
    
    PhysAlk7[j]=(psi2[j]+fract*psi1)*(Alk4[j]-Alk7[j])/V7
    CarbAlk7[j]=CarbC7[j]*2
    BioAlk7[j]=BIO_ALK*(-BioC7[j])*16/RCP7[j]
    
    PhysP7[j]=(psi2[j]+fract*psi1)*(P4[j]-P7[j])/V7
    BioP7[j]=BioC7[j]*RPC7[j]
    
    CC[j]=(CarbC1[j])*V1+(CarbC2[j])*V2+(CarbC3[j])*V3+(CarbC4[j])*V4+(CarbC5[j])*V5+(CarbC6[j])*V6+(CarbC7[j])*V7+RiverC[j]*V1
 
#%%
Source=RiverAlk*V1
Sink=CarbAlk1*V1+CarbAlk2*V2+CarbAlk3*V3+CarbAlk4*V4+CarbAlk5*V5+CarbAlk6*V6+CarbAlk7*V7

b1=1
b2=1
b3=400/900
b4=(V2+V7)/V4
b5=500/2500
b6=0
b7=1

Vtot=V1*b1+V2*b2+V3*b3+V4*b4+V5*b5+V6*b6+V7*b7

DIC=(C1*V1*b1+C2*V2*b2+C3*V3*b3+C4*V4*b4+C5*V5*b5+C6*V6*b6+C7*V7*b7)/Vtot*1e6/rho
CO3=(CO31*V1*b1+CO32*V2*b2+CO33*V3*b3+CO34*V4*b4+CO35*V5*b5+CO36*V6*b6+CO37*V7*b7)/Vtot*1e6
HCO3=(HCO31*V1*b1+HCO32*V2*b2+HCO33*V3*b3+HCO34*V4*b4+HCO35*V5*b5+HCO36*V6*b6+HCO37*V7*b7)/Vtot*1e6
H2CO3=(H2CO31*V1*b1+H2CO32*V2*b2+H2CO33*V3*b3+H2CO34*V4*b4+H2CO35*V5*b5+H2CO36*V6*b6+H2CO37*V7*b7)/Vtot*1e6

pH=-np.log10((H1*V1*b1+H2*V2*b2+H3*V3*b3+H4*V4*b4+H5*V5*b5+H6*V6*b6+H7*V7*b7)/Vtot)
Alk=(Alk1*V1*b1+Alk2*V2*b2+Alk3*V3*b3+Alk4*V4*b4+Alk5*V5*b5+Alk6*V6*b6+Alk7*V7*b7)/Vtot*1e6/rho

BIO=(BioC1*V1*b1+BioC2*V2*b2+BioC3*V3*b3+BioC4*V4*b4+BioC5*V5*b5+BioC6*V6*b6+BioC7*V7*b7)/Vtot
CARB=(CarbC1*V1*b1+CarbC2*V2*b2+CarbC3*V3*b3+CarbC4*V4*b4+CarbC5*V5*b5+CarbC6*V6*b6+CarbC7*V7*b7)/Vtot
#%%
Y=2.85


x_tick = [0,0.57, 2*0.57, 3*0.57, 4*0.57,Y]
labels = [0,2,4,6,8,10]

min_PH=np.where(pH==min(pH))


ti_3=np.linspace(0,3*max(ti_2),3*10000)

figure1 = plt.figure(1,figsize=(10,12.5))

ax1=figure1.add_subplot(421) 
plt.plot(ti_3/max(ti_2),np.roll(np.tile(DIC,3),-min_PH[0]),linewidth=LW)
plt.xlim(0,Y)
plt.grid()
plt.xticks(x_tick, labels)
plt.ylabel('[DIC] [$\mu$mol/kg]',fontsize=FS)
plt.xlabel('Time [yr]',fontsize=FS)
plt.xticks(fontsize=FS-3)
plt.yticks(fontsize=FS-3)

ax2=figure1.add_subplot(422) 
plt.plot(ti_3/max(ti_2),np.roll(np.tile(CO3,3),-min_PH[0]),linewidth=LW)
plt.xlim(0,Y)
plt.grid()
plt.xticks(x_tick, labels)
plt.ylabel('[CO$_3^{2-}$] [$\mu$mol/kg]',fontsize=FS)
plt.xlabel('Time [yr]',fontsize=FS)
plt.xticks(fontsize=FS-3)
plt.yticks(fontsize=FS-3)

ax3=figure1.add_subplot(423) 
plt.plot(ti_3/max(ti_2),np.roll(np.tile(H2CO3,3),-min_PH[0]),linewidth=LW)
plt.xlim(0,Y)
plt.grid()
plt.xticks(x_tick, labels)
plt.ylabel('[CO$_2$] (oceanic) [$\mu$mol/kg]',fontsize=FS)
plt.xlabel('Time [yr]',fontsize=FS)
plt.xticks(fontsize=FS-3)
plt.yticks(fontsize=FS-3)

ax4=figure1.add_subplot(424) 
plt.plot(ti_3/max(ti_2),np.roll(np.tile(HCO3,3),-min_PH[0]),linewidth=LW)
plt.xlim(0,Y)
plt.grid()
plt.xticks(x_tick, labels)
plt.ylabel('[HCO$_3^-$] [$\mu$mol/kg]',fontsize=FS)
plt.xlabel('Time [yr]',fontsize=FS)
plt.xticks(fontsize=FS-3)
plt.yticks(fontsize=FS-3)

ax5=figure1.add_subplot(425) 
plt.plot(ti_3/max(ti_2),np.roll(np.tile(pH,3),-min_PH[0]),linewidth=LW)
plt.xlim(0,Y)
plt.grid()
plt.xticks(x_tick, labels)
plt.ylabel('pH',fontsize=FS)
plt.xlabel('Time [yr]',fontsize=FS)
plt.xticks(fontsize=FS-3)
plt.yticks(fontsize=FS-3)

ax6=figure1.add_subplot(426) 
plt.plot(ti_3/max(ti_2),np.roll(np.tile(Alk,3),-min_PH[0]),linewidth=LW)
plt.xlim(0,Y)
plt.grid()
plt.xticks(x_tick, labels)
plt.ylabel('[Alk] [$\mu$mol/kg]',fontsize=FS)
plt.xlabel('Time [yr]',fontsize=FS)
plt.xticks(fontsize=FS-3)
plt.yticks(fontsize=FS-3)

ax7=figure1.add_subplot(427) 
plt.plot(ti_3/max(ti_2),np.roll(np.tile(BIO,3),-min_PH[0]),linewidth=LW)
plt.xlim(0,Y)
plt.grid()
plt.xticks(x_tick, labels)
plt.ylabel('Alk source [Mmol/s]',fontsize=FS)
plt.xlabel('Time [yr]',fontsize=FS)
plt.xticks(fontsize=FS-3)
plt.yticks(fontsize=FS-3)

ax8=figure1.add_subplot(428) 
plt.plot(ti_3/max(ti_2),-np.roll(np.tile(CARB,3),-min_PH[0]),linewidth=LW)
plt.xlim(0,Y)
plt.grid()
plt.xticks(x_tick, labels)
plt.ylabel('|Alk sink| [Mmol/s]',fontsize=FS)
plt.xlabel('Time [yr]',fontsize=FS)
plt.xticks(fontsize=FS-3)
plt.yticks(fontsize=FS-3)

plt.tight_layout()
#plt.savefig('figure7b.png', format='png', dpi=300)

#%%
dc=0.5*(max(CO3)-min(CO3))
dw=0.5*(max(DIC)-min(DIC))

dC=80
dW=550

amp=2*np.sqrt(dC*dW)

riv=np.mean(RiverC)
river=riv/rho*1e6*365*86400*10000

tau=amp/(river)+1
print(tau*10000)

