!----------------------------------------------------------------------
!----------------------------------------------------------------------
!   SCPM: Simple Carbon Project model (PI-configuration - AMOC)
!----------------------------------------------------------------------
!----------------------------------------------------------------------

      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)
!     ---------- ----

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)

      DOUBLE PRECISION CO2,C1,Alk1,P1,C2,Alk2,P2,C3,Alk3,P3,C4,Alk4,P4,C5,Alk5,P5,C6,Alk6,P6,C7,Alk7,P7
      DOUBLE PRECISION Vat,rho,FCA,n,PerC,DC,RCP,RPC,RSiP,WSC,WSV,WCV,kCa,Pin,b,d0,fract,yr,day,gamma1,gamma2,psi1,psi2,PV,KspSed
      DOUBLE PRECISION T1,Sal1,BT1,Ca1,K01,K11,K21,KB1,KW1,KS1,KP11,KP21,KP31,Ksp1,S1,Z1,Zo1,V1,df1,eps1,fg1,Dm1,H1,CO31
      DOUBLE PRECISION T2,Sal2,BT2,Ca2,K02,K12,K22,KB2,KW2,KS2,KP12,KP22,KP32,Ksp2,S2,Z2,Zo2,V2,df2,eps2,fg2,Dm2,H2,CO32
      DOUBLE PRECISION T3,Sal3,BT3,Ca3,K03,K13,K23,KB3,KW3,KS3,KP13,KP23,KP33,Ksp3,S3,V3,df3,fg3,Dm3,H3,CO33
      DOUBLE PRECISION T4,Sal4,BT4,Ca4,K04,K14,K24,KB4,KW4,KS4,KP14,KP24,KP34,Ksp4,S4,V4,df4,fg4,Dm4,H4,CO34
      DOUBLE PRECISION T5,Sal5,BT5,Ca5,K05,K15,K25,KB5,KW5,KS5,KP15,KP25,KP35,Ksp5,S5,Z5,Zo5,V5,df5,eps5,fg5,Dm5,H5,CO35
      DOUBLE PRECISION T6,Sal6,BT6,Ca6,K06,K16,K26,KB6,KW6,KS6,KP16,KP26,KP36,Ksp6,S6,V6,df6,fg6,Dm6,H6,CO36
      DOUBLE PRECISION T7,Sal7,BT7,Ca7,K07,K17,K27,KB7,KW7,KS7,KP17,KP27,KP37,Ksp7,S7,Z7,Zo7,V7,df7,eps7,fg7,Dm7,H7,CO37
      DOUBLE PRECISION PhysC1,BioC1,air1,CarbC1,RiverC,PhysAlk1,CarbAlk1,RiverAlk,PhysP1,BioP1,RiverP
      DOUBLE PRECISION PhysC2,BioC2,air2,CarbC2,PhysAlk2,CarbAlk2,PhysP2,BioP2
      DOUBLE PRECISION PhysC3,BioC3,CarbC3,PhysAlk3,CarbAlk3,PhysP3,BioP3
      DOUBLE PRECISION PhysC4,BioC4,BioC41,CarbC4,PhysAlk4,CarbAlk4,PhysP4,BioP4
      DOUBLE PRECISION PhysC5,BioC5,air5,CarbC5,PhysAlk5,CarbAlk5,PhysP5,BioP5
      DOUBLE PRECISION PhysC6,BioC6,CarbC6,CarbC61,PhysAlk6,CarbAlk6,PhysP6,BioP6,SedP6
      DOUBLE PRECISION PhysC7,BioC7,air7,CarbC7,PhysAlk7,CarbAlk7,PhysP7,BioP7,TC,TC1,TP,TP1,TA,TA1
      DOUBLE PRECISION RA1,RA2,RA3,RA4,RA5,RA6,RA7,epstol,DUM1,DUM2,DUM5,DUM7,A1_C,A2_C,A3_C,B1_C,B2_C,B3_C
      DOUBLE PRECISION K1,K2,K3,K4,K5,K6,K7,dKC1,dKC2,dKC3,dKC4,dKC5,dKC6,dKC7,dVC1,dVC2,dVC3,dVC4,dVC5,dVC6,dVC7
      DOUBLE PRECISION Kspres1,Kspres2,Kspres3,Kspres4,Kspres5,Kspres6,Kspres7,K111,K222,K33,K44,K55,K66,K77
      DOUBLE PRECISION K10,K20,K30,K40,K50,K60,K70,BIO,AMOC,TEMP,psi2_base,CO2baseHol,Tadjust,Ant
      DOUBLE PRECISION RCP_VAR,PV_VAR,BIO_ALK,FCA_VAR,FCA_TEMP,RCP_CO2,eps_TEMP,PV1,PV2,PV5,PV7
      DOUBLE PRECISION RCP1,RCP2,RCP3,RCP4,RCP5,RCP6,RCP7,RCPD,FCA1,FCA2,FCA5,FCA7,eta
      DOUBLE PRECISION BioAlk1,BioAlk2,BioAlk3,BioAlk4,BioAlk5,BioAlk6,BioALk7,CarbC62
      DOUBLE PRECISION RPC1,RPC2,RPC3,RPC4,RPC5,RPC6,RPC7,a_RCP,b_RCP,ComC,RCP_base
      DOUBLE PRECISION eps1_base,CC1,eps2_base,eps5_base,ComA,eps7_base
      DOUBLE PRECISION hallo1,CC2,hallo2,TEMP5,Tadjust5

        
        ! Switches
        BIO = PAR(5)        ! lambda_bi
        AMOC = 0.           ! lambda_a
        TEMP = PAR(6)       ! lambda_t
        RCP_VAR=PAR(7)      ! not used
        PV_VAR=PAR(8)       ! lambda_p
        BIO_ALK=PAR(9)      ! lambda_ba
        FCA_VAR=PAR(10)     ! lambda_f
        FCA_TEMP=PAR(12)    ! not used
        RCP_CO2=PAR(13)     ! not used
        eps_TEMP=PAR(15)    ! lambda_epsilon
        T5_TEMP=PAR(16)     ! lambda_t5

        ! Start with defining the state variables
        CO2=U(1)         ! Atmospheric CO2

        C1=U(2)          ! Carbon box 1
        Alk1=U(3)        ! Alkalinity box 1
        P1=U(4)          ! Phosphate box 1
        
        C2=U(5)          ! Carbon box 2
        Alk2=U(6)        ! Alkalinity box 2
        P2=U(7)          ! Phosphate box 2
        
        C3=U(8)          ! Carbon box 3
        Alk3=U(9)        ! Alkalinity box 3
        P3=U(10)          ! Phosphate box 3

        ! Box 4 is eliminated using conservation laws (see below)
        !C4=U(8)          ! Carbon box 4
        !Alk4=U(9)        ! Alkalinity box 4
        !P4=U(10)          ! Phosphate box 4
        
        C5=U(11)          ! Carbon box 5
        Alk5=U(12)        ! Alkalinity box 5
        P5=U(13)          ! Phosphate box 5
        
        C6=U(14)          ! Carbon box 6
        Alk6=U(15)        ! Alkalinity box 6
        P6=U(16)          ! Phosphate box 6
    
        C7=U(17)          ! Carbon box 7
        Alk7=U(18)        ! Alkalinity box 7
        P7=U(19)          ! Phosphate box 7

        ComC=U(20)        ! Total carbon
        
! State the parameters of the model (first general, than per box)
        
        Vat=1.7604418363824646*(1d20)      ! Volume atmosphere
        rho=1029       ! Sea water density
        FCA=0.07      ! Rain ratio
        n=1.0           ! Order dissolution kinetics
        PerC=0.12     ! Percentage of C in CaCO3 or something
        DC=2.75*(1d-13)      ! Constant dissolution
        RSiP=16.      ! Redfield Si to P ration
        WSC=2.378234398782344*(1d-12)     ! Constant silicate weathering
        WSV=1.585489599188229*(1d-8)     ! Variable silicate weathering
        WCV=6.341958396752917*(1d-8)    ! Variable CaCO3 weathering
        kCa=4.398148148148148*(1d-6)     ! ?
        Pin=15356.44572697869    ! Phosphate input via rivers
        b=0.75       ! Martin exponent
        d0=100      ! Base depth for Martin profile
        fract=0.5    ! Fraction of GOC going into box 7 (from box 4)
        yr=31536000.      ! Seconds per year
        day=86400.     ! Seconds per day
        CO2baseHol=244.44631272*(1d-6) ! CO2 concentration Holocene
        gamma1=(29)*(1d6)  ! Upwelling box 4 and 6
        gamma2=40*(1d6)  ! Upwelling box 1 and 3
        psi1=29*(1d6)    ! GOC
        psi2_base=par(4)!19*(1d6)    ! AMOC base
        PV=3.      ! Piston velocity
        Ant= par(17)
        eta=0.81
        a_RCP=3/70
        b_RCP=117
        RCP_base=130

        RCP=RCP_CO2*(a_RCP*CO2*1d6+b_RCP)+(1-RCP_CO2)*RCP_base        ! Redfield C to P ratio
        RCPD=331        ! Redfield C to P ratio Diazotrophs (CSIRO Mk3L COAL)

        IF (P1>0) then
           RCP1=RCP_VAR*(0.9*(1/(6.9/1000*P1*1000+6/1000))+0.1*RCPD)+(1-RCP_VAR)*RCP
        ELSE
           RCP1=RCP
        ENDIF

        IF (P2>0) then
           RCP2=RCP_VAR*(0.9*1/(6.9/1000*P2*1000+6/1000)+0.1*RCPD)+(1-RCP_VAR)*RCP
        ELSE
           RCP2=RCP
        ENDIF

        IF (P3>0) then
           RCP3=RCP_VAR*(0.9*1/(6.9/1000*P3*1000+6/1000)+0.1*RCPD)+(1-RCP_VAR)*RCP
        ELSE
           RCP3=RCP
        ENDIF

        IF (P4>0) then
           RCP4=RCP_VAR*(0.9*1/(6.9/1000*P4*1000+6/1000)+0.1*RCPD)+(1-RCP_VAR)*RCP
        ELSE
           RCP4=RCP
        ENDIF

        IF (P5>0) then
           RCP5=RCP_VAR*(0.9*1/(6.9/1000*P5*1000+6/1000)+0.1*RCPD)+(1-RCP_VAR)*RCP
        ELSE
           RCP5=RCP
        ENDIF

        IF (P6>0) then
           RCP6=RCP_VAR*(0.9*1/(6.9/1000*P6*1000+6/1000)+0.1*RCPD)+(1-RCP_VAR)*RCP
        ELSE
           RCP6=RCP
        ENDIF

        IF (P7>0) then
           RCP7=RCP_VAR*(0.9*1/(6.9/1000*P7*1000+6/1000)+0.1*RCPD)+(1-RCP_VAR)*RCP
        ELSE
           RCP7=RCP
        ENDIF

        RPC1=1/RCP1      ! Redfield P to C ratio
        RPC2=1/RCP2      ! Redfield P to C ratio
        RPC3=1/RCP3      ! Redfield P to C ratio
        RPC4=1/RCP4      ! Redfield P to C ratio
        RPC5=1/RCP5      ! Redfield P to C ratio
        RPC6=1/RCP6      ! Redfield P to C ratio
        RPC7=1/RCP7      ! Redfield P to C ratio

        IF (CO2>0) then
           psi2=psi2_base*(1-AMOC*0.1*5.35*0.54*log((CO2)/CO2baseHol)) ! AMOC strength
           Tadjust=TEMP*0.54*5.35*log((CO2)/CO2baseHol) ! Adjustment temperature CO2 concentration
           Tadjust5=TEMP5*0.54*5.35*log((CO2)/CO2baseHol)
        ELSE
           psi2=psi2_base
           Tadjust=0
           Tadjust5=0
        ENDIF
        
        ! For K0
        A1_C=-60.3409
        A2_C=93.4517
        A3_C=23.3585
        B1_C=0.023517
        B2_C=-0.023656
        B3_C=0.0047036

        ! Box 1
        T1=23.34+Tadjust      ! Temperature
        Sal1=35.25   ! Salinity
        BT1=1.179e-5*Sal1 ! Total Boron
        Ca1=0.01028*Sal1/35 ! Calcium
        K10=A1_C+A2_C*(100.0/((T1+273.16)))+A3_C*log((T1+273.16)/100.0)
        K01=exp((K10+(Sal1)*(B1_C+B2_C*((T1+273.16)/100.0)+B3_C*(((T1+273.16)/100.0)**2)))) ! Weiss (1974) in mols kg-1 atm-1
        K11=10**-((3633.86/(T1+273.16))-61.2172+(9.67770*log(T1+273.16))-(0.011555*Sal1)+(0.0001152*(Sal1**2))) !Lueker et al (2000) in mol kg-1
        K21=10**-((471.78/((T1+273.16)))+25.9290-(3.16967*log((T1+273.16)))-(0.01781*Sal1)+(0.0001122*(Sal1**2)))! Lueker et al (2000) in mol kg-1
        K111=-395.8293+(6537.773/(T1+273.16))+71.595*log((T1+273.16))-0.17959*(T1+273.16)
        K1=exp(K111+(-1.78938+410.64/(T1+273.16)+0.0065453*(T1+273.16))*Sal1**0.5-0.17755*Sal1+0.0094979*Sal1**(3.0/2.0)) ! Mucci (1983)
        dVC1=-65.28+0.397*T1-0.005155*T1**2+(19.816-0.0441*T1-0.00017*T1**2)*(Sal1/35)**0.5
        dKC1=0.01847+0.0001956*T1-0.000002212*T1**2+(-0.03217-0.0000711*T1-0.000002212*T1**2)*(Sal1/35)**0.5
        Kspres1=(-dVC1+0.5*dKC1*5)*5/(83.144621*(T1+273.16))
        Ksp1=K1*exp(Kspres1)
        
        Zo1=1.1      ! Biological production
        S1=2.71425*(1d14)      ! Surface
        V1=2.71425*(1d16)      ! Volume
        df1=100.      ! Floor depth
        eps1_base=0.9     ! Biological efficiency

        ! Box 2
        T2=9.1+Tadjust       ! Temperature
        Sal2=34.27    ! Salinity
        BT2=1.179e-5*Sal2 ! Total Boron
        Ca2=0.01028*Sal2/35 ! Calcium
        K20=A1_C+A2_C*(100.0/((T2+273.16)))+A3_C*log((T2+273.16)/100.0)
        K02=exp((K20+(Sal2)*(B1_C+B2_C*((T2+273.16)/100.0)+B3_C*(((T2+273.16)/100.0)**2)))) ! Weiss (1974) in mols kg-1 atm-1
        
        K12 = 10**-((3633.86/(T2+273.16)) - 61.2172 + (9.67770 * log(T2+273.16)) - (0.011555 * Sal2) + (0.0001152 * (Sal2**2)));
        K22=10**-((471.78/((T2+273.16)))+25.9290-(3.16967*log((T2+273.16)))-(0.01781*Sal2)+(0.0001122*(Sal2**2)))! Lueker et al (2000) in mol kg-1
        K222=-395.8293+(6537.773/(T2+273.16))+71.595*log((T2+273.16))-0.17959*(T2+273.16)
        K2=exp(K222+(-1.78938+410.64/(T2+273.16)+0.0065453*(T2+273.16))*Sal2**0.5-0.17755*Sal2+0.0094979*Sal2**(3.0/2.0)) ! Mucci (1983)
        dVC2=-65.28+0.397*T2-0.005155*T2**2+(19.816-0.0441*T2-0.00017*T2**2)*(Sal2/35)**0.5
        dKC2=0.01847+0.0001956*T2-0.000002212*T2**2+(-0.03217-0.0000711*T2-0.000002212*T2**2)*(Sal2/35)**0.5
        Kspres2=(-dVC2+0.5*dKC2*0)*0/(83.144621*(T2+273.16))
        Ksp2=K2*exp(Kspres2)
    
        Zo2=4.5      ! Biological production
        S2=3.61900*(1d13)      ! Surface
        V2=9.047500*(1d15)      ! Volume
        df2=250.    ! Floor depth
        eps2_base=1.25    ! Biological efficiency

        ! Box 3
        T3=11.28     ! Temperature
        Sal3=34.91    ! Salinity
        BT3=1.179e-5*Sal3 ! Total Boron
        Ca3=0.01028*Sal3/35 ! Calcium
        K30=A1_C+A2_C*(100.0/((T3+273.16)))+A3_C*log((T3+273.16)/100.0)
        K03=exp((K30+(Sal3)*(B1_C+B2_C*((T3+273.16)/100.0)+B3_C*(((T3+273.16)/100.0)**2)))) ! Weiss (1974) in mols kg-1 atm-1
        K13=10**-((3633.86/(T3+273.16)) - 61.2172 + (9.67770 * log(T3+273.16)) - (0.011555 * Sal3) + (0.0001152 * (Sal3**2))) !Lueker et al (2000) in mol kg-1
        K23=10**-((471.78/((T3+273.16))) + 25.9290 - (3.16967 * log((T3+273.16))) - (0.01781 * Sal3) + (0.0001122 * (Sal3**2)))! Lueker et al (2000) in mol kg-1
        K33=-395.8293+(6537.773/(T3+273.16))+71.595*log((T3+273.16))-0.17959*(T3+273.16)
        K3=exp(K33+(-1.78938+410.64/(T3+273.16)+0.0065453*(T3+273.16))*Sal3**0.5-0.17755*Sal3+0.0094979*Sal3**(3.0/2.0)) ! Mucci (1983)
        dVC3=-65.28+0.397*T3-0.005155*T3**2+(19.816-0.0441*T3-0.00017*T3**2)*(Sal3/35)**0.5
        dKC3=0.01847+0.0001956*T3-0.000002212*T3**2+(-0.03217-0.0000711*T3-0.000002212*T3**2)*(Sal3/35)**0.5
        Kspres3=(-dVC3+0.5*dKC3*55)*55/(83.144621*(T3+273.16))
        Ksp3=K3*exp(Kspres3)
        
        S3=2.71425*(1d14)      ! Surface
        V3=2.442825*(1d17)      ! Volume
        df3=1000.      ! Floor depth

        ! Box 4
        T4=3.24      ! Temperature
        Sal4=34.76   ! Salinity
        BT4=1.179e-5*Sal4 ! Total Boron
        Ca4=0.01028*Sal4/35 ! Calcium
        K40=A1_C+A2_C*(100.0/((T4+273.16)))+A3_C*log((T4+273.16)/100.0)
        K04=exp((K40+(Sal4)*(B1_C+B2_C*((T4+273.16)/100.0)+B3_C*(((T4+273.16)/100.0)**2)))) ! Weiss (1974) in mols kg-1 atm-1
        K14=10**-((3633.86/(T4+273.16))-61.2172+(9.67770*log(T4+273.16))-(0.011555*Sal4)+(0.0001152*(Sal4**2))) !Lueker et al (2000) in mol kg-1
        K24=10**-((471.78/((T4+273.16)))+25.9290-(3.16967*log((T4+273.16)))-(0.01781*Sal4)+(0.0001122*(Sal4**2)))! Lueker et al (2000) in mol kg-1
        K44=-395.8293+(6537.773/(T4+273.16))+71.595*log((T4+273.16))-0.17959*(T4+273.16)
        K4=exp(K44+(-1.78938+410.64/(T4+273.16)+0.0065453*(T4+273.16))*Sal4**0.5-0.17755*Sal4+0.0094979*Sal4**(3.0/2.0)) ! Mucci (1983)
        dVC4=-65.28+0.397*T4-0.005155*T4**2+(19.816-0.0441*T4-0.00017*T4**2)*(Sal4/35)**0.5
        dKC4=0.01847+0.0001956*T4-0.000002212*T4**2+(-0.03217-0.0000711*T4-0.000002212*T4**2)*(Sal4/35)**0.5
        Kspres4=(-dVC4+0.5*dKC4*152.632)*152.632/(83.144621*(T4+273.16))
        Ksp4=K4*exp(Kspres4)
      
        S4=3.43805*(1d14)      ! Surface
        V4=5.699925*(1d17)     ! Volume
        df4=2500.      ! Floor depth

        ! Box 5
        T5=0.93+Tadjust5      ! Temperature
        Sal5=34.43    ! Salinity
        BT5=1.179e-5*Sal5 ! Total Boron
        Ca5=0.01028*Sal5/35 ! Calcium
        K50=A1_C+A2_C*(100.0/((T5+273.16)))+A3_C*log((T5+273.16)/100.0)
        K05=exp((K50 + (Sal5) * (B1_C + B2_C *((T5+273.16)/100.0) + B3_C * (((T5+273.16)/100.0)**2)))) ! Weiss (1974) in mols kg-1 atm-1


        K15=10**-((3633.86/(T5+273.16))-61.2172+(9.67770*log(T5+273.16))-(0.011555*Sal5)+(0.0001152*(Sal5**2))) !Lueker et al (2000) in mol kg-1
        K25=10**-((471.78/((T5+273.16)))+25.9290-(3.16967*log((T5+273.16)))-(0.01781*Sal5)+(0.0001122*(Sal5**2)))! Lueker et al (2000) in mol kg-1
        K55=-395.8293+(6537.773/(T5+273.16))+71.595*log((T5+273.16))-0.17959*(T5+273.16)
        K5=exp(K55+(-1.78938+410.64/(T5+273.16)+0.0065453*(T5+273.16))*Sal5**0.5-0.17755*Sal5+0.0094979*Sal5**(3.0/2.0)) ! Mucci (1983)
        dVC5=-65.28+0.397*T5-0.005155*T5**2+(19.816-0.0441*T5-0.00017*T5**2)*(Sal5/35)**0.5
        dKC5=0.01847+0.0001956*T5-0.000002212*T5**2+(-0.03217-0.0000711*T5-0.000002212*T5**2)*(Sal5/35)**0.5
        Kspres5=(-dVC5+0.5*dKC5*0)*0/(83.144621*(T5+273.16))
        Ksp5=K5*exp(Kspres5)
       
        Zo5=1.75      ! Biological production
        S5=1.80950*(1d13)! Surface
        V5=4.523750*(1d16)      ! Volume
        df5=2500.     ! Floor depth
        eps5_base=0.35    ! Biological efficiency

        ! Box 6
        T6=1.8      ! Temperature
        Sal6=34.77    ! Salinity
        BT6=1.179e-5*Sal6 ! Total Boron
        Ca6=0.01028*Sal6/35 ! Calcium
        K60=A1_C+A2_C*(100.0/((T6+273.16)))+A3_C*log((T6+273.16)/100.0)
        K06=exp((K60+(Sal6)*(B1_C+B2_C*((T6+273.16)/100.0)+B3_C*(((T6+273.16)/100.0)**2)))) ! Weiss (1974) in mols kg-1 atm-1
        K16=10**-((3633.86/(T6+273.16))-61.2172+(9.67770*log(T6+273.16))-(0.011555*Sal6)+(0.0001152*(Sal6**2))) !Lueker et al (2000) in mol kg-1
        K26=10**-((471.78/((T6+273.16)))+25.9290-(3.16967*log((T6+273.16)))-(0.01781*Sal6)+(0.0001122*(Sal6**2)))! Lueker et al (2000) in mol kg-1
        K66=-395.8293+(6537.773/(T6+273.16))+71.595*log((T6+273.16))-0.17959*(T6+273.16)
        K6=exp(K66+(-1.78938+410.64/(T6+273.16)+0.0065453*(T6+273.16))*Sal6**0.5-0.17755*Sal6+0.0094979*Sal6**(3.0/2.0)) ! Mucci (1983)
        dVC6=-65.28+0.397*T6-0.005155*T6**2+(19.816-0.0441*T6-0.00017*T6**2)*(Sal6/35)**0.5
        dKC6=0.01847+0.0001956*T6-0.000002212*T6**2+(-0.03217-0.0000711*T6-0.000002212*T6**2)*(Sal6/35)**0.5
        Kspres6=(-dVC6+0.5*dKC6*276.25)*276.25/(83.144621*(T6+273.16))
        Ksp6=K6*exp(Kspres6)
        KspSed=K6*exp((-dVC6+0.5*dKC6*400)*400/(83.144621*(T6+273.16)))
      
        S6=3.61900*(1d14)      ! Surface
        V6=5.428500*(1d17)      ! Volume
        df6=4000.    ! Floor depth

        ! Box 7
        T7=5.83+Tadjust      ! Temperature
        Sal7=34.17    ! Salinity
        BT7=1.179e-5*Sal7 ! Total Boron
        Ca7=0.01028*Sal7/35 ! Calcium
        K70=A1_C+A2_C*(100.0/((T7+273.16)))+A3_C*log((T7+273.16)/100.0)
        K07=exp((K70+(Sal7)*(B1_C+B2_C*((T7+273.16)/100.0)+B3_C*(((T7+273.16)/100.0)**2)))) ! Weiss (1974) in mols kg-1 atm-1
        K17=10**-((3633.86/(T7+273.16))-61.2172+(9.67770*log(T7+273.16))-(0.011555*Sal7)+(0.0001152*(Sal7**2))) !Lueker et al (2000) in mol kg-1
        K27=10**-((471.78/((T7+273.16)))+25.9290-(3.16967*log((T7+273.16)))-(0.01781*Sal7)+(0.0001122*(Sal7**2)))! Lueker et al (2000) in mol kg-1

        K77=-395.8293+(6537.773/(T7+273.16))+71.595*log((T7+273.16))-0.17959*(T7+273.16)
        K7=exp(K77+(-1.78938+410.64/(T7+273.16)+0.0065453*(T7+273.16))*Sal7**0.5-0.17755*Sal7+0.0094979*Sal7**(3.0/2.0)) ! Mucci (1983)
        dVC7=-65.28+0.397*T7-0.005155*T7**2+(19.816-0.0441*T7-0.00017*T7**2)*(Sal7/35)**0.5
        dKC7=0.01847+0.0001956*T7-0.000002212*T7**2+(-0.03217-0.0000711*T7-0.000002212*T7**2)*(Sal7/35)**0.5
        Kspres7=(-dVC7+0.5*dKC7*12.5)*12.5/(83.144621*(T7+273.16))
        Ksp7=K7*exp(Kspres7)
        
        Zo7=5.325      ! Biological production
        S7=3.61900*(1d13)    ! Surface
        V7=9.047500*(1d15)      ! Volume
        df7=250.     ! Floor depth
        eps7_base=0.62    ! Biological efficiency

        ! Conservation laws to eliminate variable
        ! Carbon
        TC=CO2*Vat+C1*V1+C2*V2+C3*V3+C5*V5+C6*V6+C7*V7   ! Sum of carbon without eliminated variable
        TC1=ComC*(1d22)*par(1)+Ant                      !  Total carbon in system (PgC)
        C4=1/V4*((TC1)-TC)             ! Eliminated variable
        
        IF (P1<0) then
            P1=0
        ENDIF
        ! Phosphate
        TP=P1*V1+P2*V2+P3*V3+P5*V5+P6*V6+P7*V7      ! Sum of phosphate without eliminated variable
        TP1=2.9622168229501007*(1d15)*par(1)           ! Total phosphate in system
        P4=(1/V4)*(TP1-TP)                              ! Eliminated variable
   
        !Alkalinity
        TA=Alk1*V1+Alk2*V2+Alk3*V3+Alk5*V5+Alk6*V6+Alk7*V7 ! Sum of alkalinity without eliminated variable
        TA1=3486246.73610888e+12-(3427697.69670697E+12-ComC*(1d22))*2*par(1)         ! Total alkalinity in system
        Alk4=(1/V4)*(TA1-TA)                                ! Eliminated variable

        ! Biological efficiency feedback
        eps1=eps_TEMP*Tadjust*(-0.1)+eps1_base
        eps2=eps_TEMP*Tadjust*(-0.1)+eps2_base
        eps5=eps_TEMP*Tadjust5*(-0.1)+eps5_base
        eps7=eps_TEMP*Tadjust*(-0.1)+eps7_base

       ! Piston velocity
        PV1=PV*PV_VAR*(((2116.8-136.25*T1+4.7353*T1**2-0.092307*T1**3+0.0007555*T1**4)/660)**-0.5)+PV*(1-PV_VAR)
        PV2=PV*PV_VAR*(((2116.8-136.25*T2+4.7353*T2**2-0.092307*T2**3+0.0007555*T2**4)/660)**-0.5)+PV*(1-PV_VAR)
        PV5=PV*PV_VAR*(((2116.8-136.25*T5+4.7353*T5**2-0.092307*T5**3+0.0007555*T5**4)/660)**-0.5)+PV*(1-PV_VAR)
        PV7=PV*PV_VAR*(((2116.8-136.25*T7+4.7353*T7**2-0.092307*T7**3+0.0007555*T7**4)/660)**-0.5)+PV*(1-PV_VAR)
        
       ! Carbonate chemistry (Williams and Follows, 2011)
       ! Def of RA_i only needed to start from trivial solution
        epstol = 1.0d-06
        IF (abs(Alk1).gt.epstol) then
           RA1 = (C1/rho)/(Alk1/rho)
        ELSE
           RA1 = 0.0
        ENDIF
        IF (abs(Alk2).gt.epstol) then
           RA2 = (C2/rho)/(Alk2/rho)
        ELSE
           RA2 = 0.0
        ENDIF
        IF (abs(Alk3).gt.epstol) then
           RA3 = (C3/rho)/(Alk3/rho)
        ELSE
           RA3 = 0.0
        ENDIF
        IF (abs(Alk4).gt.epstol) then
           RA4 = (C4/rho)/(Alk4/rho)
        ELSE
           RA4 = 0.0
        ENDIF
        IF (abs(Alk5).gt.epstol) then
           RA5 = (C5/rho)/(Alk5/rho)
        ELSE
           RA5 = 0.0
        ENDIF
        IF (abs(Alk6).gt.epstol) then
           RA6 = (C6/rho)/(Alk6/rho)
        ELSE
           RA6 = 0.0
        ENDIF
        IF (abs(Alk7).gt.epstol) then
           RA7 = (C7/rho)/(Alk7/rho)
        ELSE
           RA7 = 0.0
        ENDIF
        
        DM1=(1-RA1)**2*K11**2-4*K11*K21*(1-2*RA1)
        DM2=(1-RA2)**2*K12**2-4*K12*K22*(1-2*RA2)
        DM3=(1-RA3)**2*K13**2-4*K13*K23*(1-2*RA3)
        DM4=(1-RA4)**2*K14**2-4*K14*K24*(1-2*RA4)
        DM5=(1-RA5)**2*K15**2-4*K15*K25*(1-2*RA5)
        DM6=(1-RA6)**2*K16**2-4*K16*K26*(1-2*RA6)
        DM7=(1-RA7)**2*K17**2-4*K17*K27*(1-2*RA7)

        H1=0.5*((RA1-1)*K11+DM1**0.5)
        H2=0.5*((RA2-1)*K12+DM2**0.5)
        H3=0.5*((RA3-1)*K13+DM3**0.5)
        H4=0.5*((RA4-1)*K14+DM4**0.5)
        H5=0.5*((RA5-1)*K15+DM5**0.5)
        H6=0.5*((RA6-1)*K16+DM6**0.5)
        H7=0.5*((RA7-1)*K17+DM7**0.5)

        CO31=(C1/rho/K01)*(H1**2/(H1**2+H1*K11+K11*K21))*K01*K11/H1*K21/H1
        CO32=(C2/rho/K02)*(H2**2/(H2**2+H2*K12+K12*K22))*K02*K12/H2*K22/H2
        CO33=(C3/rho/K03)*(H3**2/(H3**2+H3*K13+K13*K23))*K03*K13/H3*K23/H3
        CO34=(C4/rho/K04)*(H4**2/(H4**2+H4*K14+K14*K24))*K04*K14/H4*K24/H4
        CO35=(C5/rho/K05)*(H5**2/(H5**2+H5*K15+K15*K25))*K05*K15/H5*K25/H5
        CO36=(C6/rho/K06)*(H6**2/(H6**2+H6*K16+K16*K26))*K06*K16/H6*K26/H6
        CO37=(C7/rho/K07)*(H7**2/(H7**2+H7*K17+K17*K27))*K07*K17/H7*K27/H7

        ! Rain ratio

        IF (Ca1*CO31/Ksp1-1>0) then
           FCA1=FCA_VAR*0.022*(Ca1*CO31/Ksp1-1)**eta+(1-FCA_VAR)*FCA
        ELSE
           FCA1=FCA
        ENDIF

        IF (Ca2*CO32/Ksp2-1>0) then
           FCA2=FCA_VAR*0.022*(Ca2*CO32/Ksp2-1)**eta+(1-FCA_VAR)*FCA
        ELSE
           FCA2=FCA
        ENDIF

        IF (Ca5*CO35/Ksp5-1>0) then
           FCA5=FCA_VAR*0.022*(Ca5*CO35/Ksp5-1)**eta+(1-FCA_VAR)*FCA
        ELSE
           FCA5=FCA
        ENDIF

        IF (Ca7*CO37/Ksp7-1>0) then
           FCA7=FCA_VAR*0.022*(Ca7*CO37/Ksp7-1)**eta+(1-FCA_VAR)*FCA
        ELSE
           FCA7=FCA
        ENDIF

        FCA5=FCA_TEMP*FCA5*(T5+2)/7+(1-FCA_TEMP)*FCA5
        
        ! New biology
        Z1=(1-BIO)*Zo1+BIO*(gamma2*P3*yr*1+Pin*yr)*eps1*RCP1*1/S1
        Z2=(1-BIO)*Zo2+BIO*psi2*P3*yr*eps2*RCP2*1/S2
        Z5=(1-BIO)*Zo5+BIO*fract*psi1*P7*yr*eps5*RCP5*1/S5
        Z7=(1-BIO)*Zo7+BIO*(psi2+fract*psi1)*P4*yr*eps7*RCP7*1/S7

        ! Define the different components (physical, biological, air exchange, carbonate, and river)
        ! Box 1
        PhysC1=gamma2*(C3-C1)/V1
        BioC1=-(Z1*S1/(V1*yr))*(df1/d0)**(-b)
        DUM1 = par(3)*(C1/rho/K01)*((H1**2)/((H1**2)+H1*K11+K11*K21))+(1-par(3))*C1
        air1=(rho*K01*PV1/day*S1*(CO2-DUM1))/V1
        CarbC1=-(Z1*S1/(V1*yr))*FCA1+(CO31*Ca1*rho*kCa*(1-min((CO31*Ca1/Ksp1),n))**n*PerC+DC)
        RiverC=WSC+(WSV+WCV)*CO2

        PhysAlk1=gamma2*(Alk3-Alk1)/V1
        CarbAlk1=2*CarbC1
        RiverAlk=2*RiverC
        BioAlk1=BIO_ALK*(-BioC1)*16/RCP1

        PhysP1=gamma2*(P3-P1)/V1
        BioP1=BioC1*RPC1
        RiverP=Pin/V1

        ! Box 2
        PhysC2=psi2*(C3-C2)/V2
        BioC2=-(Z2*S2/(V2*yr))*(df2/d0)**(-b)
        DUM2 = par(3)*(C2/rho/K02)*((H2**2)/((H2**2)+H2*K12+K12*K22))+ (1-par(3))*C2
        air2=(rho*K02*PV2/day*S2*(CO2-DUM2))/V2
        CarbC2=-(Z2*S2/(V2*yr))*FCA2+(CO32*Ca2*rho*kCa*(1-min((CO32*Ca2/Ksp2),n))**n*PerC+DC)

        PhysAlk2=psi2*(Alk3-Alk2)/V2
        CarbAlk2=2*CarbC2
        BioAlk2=BIO_ALK*(-BioC2)*16/RCP2

        PhysP2=psi2*(P3-P2)/V2
        BioP2=BioC2*RPC2

        ! Box 3
        PhysC3=gamma2*(C1-C3)/V3+psi2*(C7-C3)/V3
        BioC3=Z1*S1/(V3*yr)*((df1/d0)**(-b)-(df3/d0)**(-b))
        CarbC3=(CO33*Ca3*rho*kCa*(1-min((CO33*Ca3/Ksp3),n))**n*PerC+DC)

        PhysAlk3=gamma2*(Alk1-Alk3)/V3+psi2*(Alk7-Alk3)/V3
        CarbAlk3=2*CarbC3
        BioAlk3=BIO_ALK*(-BioC3)*16/RCP3

        PhysP3=gamma2*(P1-P3)/V3+psi2*(P7-P3)/V3
        BioP3=BioC3*RPC3

        ! Box 4
        PhysC4=gamma1*(C6-C4)/V4+psi1*(C6-C4)/V4+psi2*(C2-C4)/V4
        BioC41=(Z1*S1*((df3/d0)**(-b)-(df4/d0)**(-b))+Z2*S2*((df2/d0)**(-b)-(df4/d0)**(-b))+Z7*S7*((df7/d0)**(-b)-(df4/d0)**(-b)))
        BioC4=BioC41/(V4*yr)
        CarbC4=(CO34*Ca4*rho*kCa*(1-min((CO34*Ca4/Ksp4),n))**n*PerC+DC)

        PhysAlk4=gamma1*(Alk6-Alk4)/V4+psi1*(Alk6-Alk4)/V4+psi2*(Alk2-Alk4)/V4
        CarbAlk4=2*CarbC4
        BioAlk4=BIO_ALK*(-BioC4)*16/RCP4

        PhysP4=gamma1*(P6-P4)/V4+psi1*(P6-P4)/V4+psi2*(P2-P4)/V4
        BioP4=BioC4*RPC4

        ! Box 5
        PhysC5=(1-fract)*psi1*(C4-C5)/V5+fract*psi1*(C7-C5)/V5
        BioC5=-(Z5*S5/(V5*yr))*(df5/d0)**(-b)
        DUM5 = par(3)*(C5/rho/K05)*((H5**2)/((H5**2)+H5*K15+K15*K25))+(1-par(3))*C5
        air5=(rho*K05*PV5/day*S5*(CO2-DUM5))/V5
        CarbC5=-(Z5*S5/(V5*yr))*FCA5+(CO35*Ca5*rho*kCa*(1-min((CO35*Ca5/Ksp5),n))**n*PerC+DC)

        PhysAlk5=(1-fract)*psi1*(Alk4-Alk5)/V5+fract*psi1*(Alk7-Alk5)/V5
        CarbAlk5=2*CarbC5
        BioAlk5=BIO_ALK*(-BioC5)*16/RCP5

        PhysP5=(1-fract)*psi1*(P4-P5)/V5+fract*psi1*(P7-P5)/V5
        BioP5=BioC5*RPC5

        ! Box 6
        PhysC6=psi1*(C5-C6)/V6+gamma1*(C4-C6)/V6
        BioC6=(Z1*S1+Z2*S2+Z5*S5+Z7*S7)/(V6*yr)*((df4/d0)**(-b))!-(4000/d0)**(-b))
        
        CarbC61=(CO36*Ca6*rho*kCa*((1-min((CO36*Ca6/Ksp6),n))**n))*PerC+DC;
        CarbC62=(CO36*Ca6*rho*kCa*((1-min((CO36*Ca6/KspSed),n))**n))*PerC+DC;
        CarbC6=CarbC61+CarbC62;

        PhysAlk6=psi1*(Alk5-Alk6)/V6+gamma1*(Alk4-Alk6)/V6
        CarbAlk6=2*CarbC6
        BioAlk6=BIO_ALK*(-BioC6)*16/RCP6

        PhysP6=psi1*(P5-P6)/V6+gamma1*(P4-P6)/V6
        BioP6=BioC6*RPC6
        SedP6=-Pin/V6

        ! Box 7
        PhysC7=(psi2+fract*psi1)*(C4-C7)/V7
        BioC7=-(Z7*S7)/(V7*yr)*(df7/d0)**(-b)
        DUM7 = par(3)*(C7/rho/K07)*((H7**2)/((H7**2)+H7*K17+K17*K27))+(1-par(3))*C7
        air7=(rho*K07*PV7/day*S7*(CO2-DUM7))/V7
        CarbC7=-(Z7*S7/(V7*yr))*FCA7+(CO37*Ca7*rho*kCa*(1-min((CO37*Ca7/Ksp7),n))**n*PerC+DC)

        PhysAlk7=(psi2+fract*psi1)*(Alk4-Alk7)/V7
        CarbAlk7=CarbC7*2
        BioAlk7=BIO_ALK*(-BioC7)*16/RCP7

        PhysP7=(psi2+fract*psi1)*(P4-P7)/V7
        BioP7=BioC7*RPC7

        CC1=(CarbC1)*V1+(CarbC2)*V2+(CarbC3)*V3+(CarbC4)*V4
        CC2=(CarbC5)*V5+(CarbC6)*V6+(CarbC7)*V7+RiverC*V1
        
        ! State the ODEs
        F(1)=(-(air1*V1+air2*V2+air5*V5+air7*V7)/Vat)     ! CO2

        F(2)=(PhysC1+air1+par(2)*(BioC1+CarbC1+RiverC))      ! C1
        F(3)=(PhysAlk1+par(2)*(CarbAlk1+RiverAlk+BioAlk1))           ! Alk1
        F(4)=(PhysP1+par(2)*(BioP1+RiverP))      ! P1

        F(5)=(PhysC2+air2+par(2)*(BioC2+CarbC2))               ! C2
        F(6)=(PhysAlk2+par(2)*(CarbAlk2+BioAlk2))                   ! Alk2
        F(7)=(PhysP2+par(2)*BioP2)                        ! P2
    
        
        F(8)=(PhysC3+par(2)*(BioC3+CarbC3))                    ! C3
        F(9)=(PhysAlk3+par(2)*(CarbAlk3+BioAlk3))                    ! Alk3
        F(10)=(PhysP3+par(2)*BioP3)                         ! P3

		! Box 4 is eliminated using conservation laws (see above)
        !F(8)=(PhysC4+par(2)*(BioC4+CarbC4))                   ! C4
        !F(9)=(PhysAlk4+par(2)*(CarbAlk4+BioAlk4))                    ! Alk4
        !F(10)=(PhysP4+par(2)*BioP4)                       ! P4

        F(11)=(PhysC5+air5+par(2)*(BioC5+CarbC5))              ! C5
        F(12)=(PhysAlk5+par(2)*(CarbAlk5+BioAlk5))                   ! Alk5
        F(13)=(PhysP5+par(2)*BioP5)                        ! P5

        F(14)=PhysC6+par(2)*(BioC6+CarbC6)                   ! C6
        F(15)=PhysAlk6+par(2)*(CarbAlk6+BioAlk6)                   ! Alk6
        F(16)=PhysP6+par(2)*(BioP6+SedP6)                  ! P6
        
        F(17)=(PhysC7+air7+par(2)*(BioC7+CarbC7))              ! C7
        F(18)=(PhysAlk7+par(2)*(CarbAlk7+BioAlk7))                   ! Alk7
        F(19)=(PhysP7+par(2)*BioP7)                        ! P7

        F(20)=(CC1+CC2)/(1d22)                              ! Total carbon

      END SUBROUTINE FUNC
!----------------------------------------------------------------------

      SUBROUTINE STPNT(NDIM,U,PAR,T)
!     ---------- -----

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      INTEGER I
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T
      

! Parameters
      PAR(1) = 1.0  ! Homotopy for total quantities
      PAR(2) = 1.0  ! Homotopy for non-physical processes
      PAR(3) = 1.0  ! Homotopy carbon chemistry
      PAR(4) = 19*(1d6)
      PAR(5) = 0
      PAR(6) = 0
      PAR(7) = 0
      PAR(8) = 0
      PAR(9) = 0
      PAR(10) = 0
      PAR(12) = 0
      PAR(13) = 0
      PAR(15) = 0
      PAR(16) = 0
      
! Special solution from transient code
! NF_PI_v11: initial conditions
       U(1)=244.44631272/(1d6)

       U(2)=2052564.19581549/(1d6)
       U(3)=2380069.11789182/(1d6)
       U(4)=265.47289588/(1d6)
       
       U(5)=2133025.57003298/(1d6)
       U(6)=2345516.44628329/(1d6)
       U(7)=650.65687622/(1d6)
       
       U(8)=2239245.96504157/(1d6)
       U(9)=2383305.79466709/(1d6)
       U(10)=1702.2373908/(1d6)
       
       !U(8)=2380702.42810271/(1d6)
       !U(9)=2408437.85299461/(1d6)
       !U(10)=2281.38038779/(1d6)
       
       U(11)=2257315.75411338/(1d6)
       U(12)=2391753.61173295/(1d6)
       U(13)=1904.67641726/(1d6)
       
       U(14)=2364319.45614937/(1d6)
       U(15)=2423644.93665859/(1d6)
       U(16)=2086.22950376/(1d6)
       
       U(17)=2219044.2894214/(1d6)
       U(18)=2383048.51550966/(1d6)
       U(19)=1575.61814193/(1d6)

       U(20)=3427697.69670697/(1d10)
       !TA=3486246.73610888
!     Trivial Solution
      !DO I = 1, 19
       ! U(I) = 0.0
      !ENDDO
      END SUBROUTINE STPNT
!----------------------------------------------------------------------

      SUBROUTINE BCND
      END SUBROUTINE BCND

    SUBROUTINE ICND(NDIM,PAR,ICP,NINT,U,UOLD,UDOT,UPOLD,FI,IJAC,DINT)
      !---------- ----

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NDIM, ICP(*), NINT, IJAC
    DOUBLE PRECISION, INTENT(IN) :: PAR(*)
    DOUBLE PRECISION, INTENT(IN) :: U(NDIM), UOLD(NDIM), UDOT(NDIM), UPOLD(NDIM)
    DOUBLE PRECISION, INTENT(OUT) :: FI(NINT)
    DOUBLE PRECISION, INTENT(INOUT) :: DINT(NINT,*)
    

    END SUBROUTINE ICND

      SUBROUTINE FOPT
      END SUBROUTINE FOPT

      SUBROUTINE PVLS
      END SUBROUTINE PVLS
