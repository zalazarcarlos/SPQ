#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import matplotlib.pyplot as plt
from matplotlib import style
style.use('fivethirtyeight')
import numpy as np
import pandas as pd
from scipy.optimize import fsolve
#from IPython.display impor display
from SPQentrada import *

#curva de la bomba
def bomba1(H):
    a=6.5035119
    b=-2.657738e-018
    c=0.48508288
    d=-0.00080261045
    g=2.3254058
    y = a + b * np.exp(H/c) + d * np.exp(H/g)
    return y


#Datos iniciales
p1 = 1
p2 = 1.5

vp_ini = 0.5
Kv_bar = 5.0

d_in = 1.0 #pulgada
L1 = 25.0
L2 = 25.0
L3 = 25.0
fd = 0.03

rho = 1000  #kg/m3
p_vap_mmHG = 100
p_vap_pasc = p_vap_mmHG * 133.3

dt = 2


#Datos Escalon P1
t0_e1= 500
A_e1 = 0

#Datos Escalon P2
t0_e2= 500
A_e2 = 0

#Datos Escalon valvula
t0_ev= 500
A_ev = 0.25

#Datos Rampa P1
t0_r1 = 500
pend1 = 0
dt1 = 500

#Datos Rampa P2
t0_r2 = 500
pend2 = 0.00
dt2 = 500

#Datos Rampa valvula
t0_rv = 500
pendv = 0
dtv = 2000

#Datos ExpDecreciente P1
t0_exd1 = 500
tau1 = 100
A_exd1 = 0

#Datos ExpDecreciente P2
t0_exd2 = 500
tau2 = 100
A_exd2 = 0

#Datos ExpDecreciente valvula
t0_exdv = 500
tauv = 100
A_exdv = 0

#límites "y". Cambio de escala
escalapresion = (0,10)
escalaF = (1,7)
escalaVP = (0,1)



def ejercicio2():
    
    #Autocalculados
    Kv_SI = Kv_bar * 0.000000878
    G = rho / 1000
    d_m = d_in * 0.0254
    area = np.pi * d_m**2 / 4
    
    T = []
    P1 = []
    P2 = []
    VP = []
    H = []
    F = []
    ANPA = []
    PA_bar = []
    PB_bar =[]
    PC_bar =[]
    PD_bar =[]
    V = []

    for i in range(0,1001,1):
        t = i * dt
        T = np.append(T, t)
        p1t = p1 \
        + ESCALON(t0_e1, A_e1, t)\
        + RAMPA(t0_r1, pend1, dt1, t)\
        + ExpDecr(t0_exd1, tau1, A_exd1, t)
        P1 = np.append(P1, p1t)
    
        p2t = p2 \
        + ESCALON(t0_e2, A_e2, t)\
        + RAMPA(t0_r2, pend2, dt2, t)\
        + ExpDecr(t0_exd2, tau2, A_exd2, t)
        P2 = np.append(P2, p2t)
        
        vpt = vp_ini \
        + ESCALON(t0_ev, A_ev, t)\
        + RAMPA(t0_rv, pendv, dtv, t)\
        + ExpDecr(t0_exdv, tauv, A_exdv, t)
        VP = np.append(VP, vpt)
        
    
        p1_pascal = p1t * 100000
        p2_pascal = p2t * 100000
        
      
        def FOsolv (x):
            
            v = x[0]
            pA_pascal = x[1]
            pB_pascal = x[2]
            pC_pascal = x[3]
            pD_pascal = x[4]
            
            #caida presión primer tramo
            F01 = (p1_pascal - pA_pascal) / (rho*9.8) \
                  - fd * L1/d_m * v**2/(2*9.8)
            #bomba      
            h = (pB_pascal - pA_pascal)/(9.8 * rho)
            
            f = bomba1(h)
            
            F02 = f / (area * 3600) - v # m/s
            
            #caida presión segundo tramo
            F03 = (pB_pascal - pC_pascal) / (rho*9.8) \
                  - fd * L2/d_m * v**2/(2*9.8)
            
            #Valvula
            F04 = (Kv_SI * vpt * ((pC_pascal -  pD_pascal)/ G)**0.5)\
                  - v*area
            
            #caida presion tercer tramo
            F05 = (pD_pascal - p2_pascal) / (rho*9.8) \
                  - fd * L3/d_m * v**2/(2*9.8)
            
            return [F01, F02, F03, F04, F05] 
        
        semilla = [1.5, 60000, 200000, 190000, 170000]
    
        sol = fsolve(FOsolv, semilla)
        v = sol[0]
        V = np.append(V,v)
        
        pA_pascal = sol[1]
        PA_bar = np.append(PA_bar, pA_pascal/100000)
        
        pB_pascal = sol[2]
        PB_bar = np.append(PB_bar, pB_pascal/100000)
        
        pC_pascal = sol[3]
        PC_bar = np.append(PC_bar, pC_pascal/100000)
        
        pD_pascal = sol[4]
        PD_bar = np.append(PD_bar, pD_pascal/100000)
        
        h = (pB_pascal - pA_pascal)/(9.8 * rho)
        H = np. append(H, h)
        
        f = v*area
        
        F = np. append(F, f)
        
        
        hsuc = v**2 / (2 * 9.8) + p1_pascal / (rho * 9.8)
        
        ANPAi = hsuc - p_vap_pasc / (9.8 * rho)
        ANPA = np.append(ANPA, ANPAi)
        
        ANPA_bar = ANPA * rho * 9.8 / 100000
        P1_cavit = P1 - ANPA_bar
        
        #print(' Pa=', pA_pascal/100000, ' Pb=', pB_pascal/100000, ' H=', h)

    F_m3h= F * 3600

    #creacion de tabla df para pandas:
    valores = np.vstack((T, P1, PA_bar, H, PB_bar, PC_bar,
                         VP, PD_bar, P2, V, F_m3h))
    valores = valores.T
    
    columnas = ['tiempo_s', 'p1_bar', 'pA', 'H_m', 'pB','pC', 
                'VP', 'pD', 'p2_bar', 'vel', 'F_m3/h']
    df = pd.DataFrame(valores, columns=columnas)

    #gráficos
    fig = plt.figure(figsize=(8,5))
    
    ax1 = plt.subplot2grid((6,1),(0,0), rowspan=4)
#    axv = plt.subplot2grid((6,1),(4,0), rowspan=2)
    
    ax1.plot(T, P1, label="P1")
    ax1.plot(T, P2, label="P2")
    ax1.plot(T, P1_cavit, label="P1 cavitac")
    ax1.plot(T, ANPA, label = "ANPA")
    ax1.set_ylabel('$Presion (bar)$')
    ax1.set_xticklabels([])

    ax2 = ax1.twinx()
    ax2.plot(T, F_m3h, 'g-')
    ax2.set_ylabel('$F (m^3/h)$', color='g')
    ax2.tick_params('y', colors='g')
    
    #limites en y (Corrección de escala)
    ax1.set_ylim(escalapresion)
    ax2.set_ylim(escalaF)
    ax1.legend()
    
#    axv.plot(T, VP, 'm-')
#    axv.set_ylabel('$apert$ $valvula$')
#    axv.set_xlabel('$Tiempo (s)$')
#    axv.set_ylim(escalaVP)

    tabla = df[(df.tiempo_s == 450) | (df.tiempo_s == 600) | (df.tiempo_s == 1000)]
    print (tabla)
    
    plt.show()
    
ejercicio2()