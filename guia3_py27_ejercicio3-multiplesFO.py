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

#Datos iniciales
p1 = 2
p2 = 1
vp_ini = 0.5

d_in = 1.0 #pulgada
L1 = 25.0
L2 = 25.0
fd = 0.03

Kv_bar = 5.0
rho = 1000  #kg/m3

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
pend2 = 0
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
escalapresion = (0,3)
escalaF = (1,7)
escalaVP = (0,1)



def ejercicio1():
    
    #Autocalculados
    Kv_SI = Kv_bar * 0.000000878
    G = rho / 1000
    d_m = d_in * 0.0254
    area = np.pi * d_m**2 / 4
    
    T = []
    P1 = []
    P2 = []
    VP = []
    F = []

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
        
        # revisar caso
        def FOsolv (x):
            
            v = x[0]
            pA_pascal = x[1]
            pB_pascal = x[2]
            
            F01 = (p1_pascal - pA_pascal) / (rho*9.8) \
                  - fd * L1/d_m * v**2/(2*9.8)
            
            F02 = (pB_pascal - p2_pascal) / (rho*9.8) \
                  - fd * L2/d_m * v**2/(2*9.8)
            
            F03 = (Kv_SI * vpt * ((pA_pascal -  pB_pascal)/ G)**0.5) - v*area
            
            return [F01, F02, F03] 
        
        semilla = [1.5, 170000, 120000]
    
        sol = fsolve(FOsolv, semilla)
        v = sol[0]
        
        pA_pascal = sol[1]
        pB_pascal = sol[2]
        print (v)
        print (pA_pascal)
        print (pB_pascal)
        
        #sQr = (Kv_SI * vpt * ((pA_pascal -  pB_pascal)/ G)**0.5) - v*area
        
        f = v*area
        
        F = np. append(F, f)

    F_m3h= F * 3600

    #creacion de tabla df para pandas:
    valores = {'tiempo_s':T,
               'p1_bar':P1, 'p2_bar':P2,
               'apert_val':VP,
               'F_m3/s':F, 'F_m3/h':F_m3h}
    columnas = ['tiempo_s', 'p1_bar', 'p2_bar',
                'apert_val', 'F_m3/s', 'F_m3/h']
    df = pd.DataFrame(valores, columns=columnas)

    #gráficos
    fig = plt.figure(figsize=(8,5))
    
    ax1 = plt.subplot2grid((6,1),(0,0), rowspan=4)
    axv = plt.subplot2grid((6,1),(4,0), rowspan=2)
    
    ax1.plot(T, P1, label="P1")
    ax1.plot(T, P2, label="P2")
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
    
    axv.plot(T, VP, 'm-')
    axv.set_ylabel('$apert$ $valvula$')
    axv.set_xlabel('$Tiempo (s)$')
    axv.set_ylim(escalaVP)

    tabla = df[(df.tiempo_s == 450) | (df.tiempo_s == 600) | (df.tiempo_s == 1000)]
    print (tabla)
    
    plt.show()
    
ejercicio1()