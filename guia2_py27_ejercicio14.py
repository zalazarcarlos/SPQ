#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import matplotlib.pyplot as plt
from matplotlib import style
style.use('fivethirtyeight')
import numpy as np
import pandas as pd
from scipy.optimize import newton
#from IPython.display impor display
from SPQentrada import *

#Datos iniciales
p1 = 2
p2 = 1
d_in = 1 
L = 51
dt = 2
fd = 0.03
rho = 1000  #kg/m3

K = 0.9 #codo 90º

L_mano1 = 10  #tramo hasta el 1º manómetro
L_mano2 = 10 + 7 
L_mano3 = 10 + 14
L_mano4 = 10 + 21

#Datos Escalon P1
t0_e1= 500
A_e1 = 1

#Datos Escalon P2
t0_e2= 500
A_e2 = 0

#Datos Rampa P1
t0_r1 = 500
pend1 = 0
dt1 = 500

#Datos Rampa P2
t0_r2 = 500
pend2 = 0.001
dt2 = 500

#Datos ExpDecreciente P1
t0_exd1 = 500
tau1 = 100
A_exd1 = 0

#Datos ExpDecreciente P2
t0_exd2 = 500
tau2 = 100
A_exd2 = 0

#límites "y". Cambio de escala
escalapresion = (0,4)
escalaF = (1,7)



def ejercicio14():
    
    #conversión automática de datos
    d_m = d_in * 0.0254
    area = np.pi * d_m**2 / 4
    
    T = []
    P1 = []
    P2 = []
    P_MANO1 = []
    P_MANO2 = []
    P_MANO3 = []
    P_MANO4 = []
    V = []
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
    
        p1_pascal = p1t * 100000
        p2_pascal = p2t * 100000
        
        #En lugar de despejar v, resuelvo para F0=0
        def FOsolv (v):
            
            F0 = (p1_pascal - p2_pascal) / rho \
            - fd * L * v ** 2 * 0.5 / d_m \
            - 2 * (K * (v**2) / 2)
            return F0
        
        sol = newton(FOsolv, 1.5)
        v = sol
    
        V = np.append(V, v)
        
        
        # 1º Manómetro
        p_mano_pasca1 = rho * (p1_pascal/rho \
                               - fd * (L_mano1 / d_m) * (v**2 /2) \
                               - K * (v**2 /2))
                               
        p_mano1 = p_mano_pasca1 / 100000
        P_MANO1 = np.append(P_MANO1, p_mano1)
        
        # 2º Manómetro
        p_mano_pasca2 = rho * (p1_pascal/rho \
                               - fd * (L_mano2 / d_m) * (v**2 /2) \
                               - K * (v**2 /2))
        
        p_mano2 = p_mano_pasca2 / 100000
        P_MANO2 = np.append(P_MANO2, p_mano2)
                               
        # 3º Manómetro
        p_mano_pasca3 = rho * (p1_pascal/rho \
                               - fd * (L_mano3 / d_m) * (v**2 /2) \
                               - K * (v**2 /2))
        
        p_mano3 = p_mano_pasca3 / 100000
        P_MANO3 = np.append(P_MANO3, p_mano3)
                               
        # 4º Manómetro
        p_mano_pasca4 = rho * (p1_pascal/rho \
                               - fd * (L_mano4 / d_m) * (v**2 /2) \
                               - K * (v**2 /2))
        
        p_mano4 = p_mano_pasca4 / 100000
        P_MANO4 = np.append(P_MANO4, p_mano4)
    
        f = v * area
        F = np. append(F, f)

    F_m3h= F * 3600

    #creacion de tabla df para pandas:
    valores = {'tiempo_s':T,
               'p1_bar':P1,
               'pm1_bar':P_MANO1,
               'pm2_bar':P_MANO2,
               'pm3_bar':P_MANO3,
               'pm4_bar':P_MANO4,
               'p2_bar':P2,
               'velc_m/s':V,
               'F_m3/s':F, 'F_m3/h':F_m3h}
    columnas = ['tiempo_s','p1_bar',
                'pm1_bar','pm2_bar',
                'pm3_bar','pm4_bar',
                'p2_bar', 'velc_m/s', 'F_m3/s', 'F_m3/h']
    df = pd.DataFrame(valores, columns=columnas)

    #gráficos
    fig, ax1 = plt.subplots()

    ax1.plot(T, P1, label="P1")
    ax1.plot(T, P2, label="P2")
    ax1.set_xlabel('$Tiempo (s)$')
    ax1.set_ylabel('$Presion (bar)$')

    ax2 = ax1.twinx()
    ax2.plot(T, F_m3h, 'g-')
    ax2.set_ylabel('$F (m^3/h)$', color='g')
    ax2.tick_params('y', colors='g')

    #limites en y (Corrección de escala)
    ax1.set_ylim(escalapresion)
    ax2.set_ylim(escalaF)
    ax1.legend()

    tabla = df[(df.tiempo_s == 450) | (df.tiempo_s == 600) | (df.tiempo_s == 1000)]
    print (tabla)
    
    plt.show()
    
ejercicio14()