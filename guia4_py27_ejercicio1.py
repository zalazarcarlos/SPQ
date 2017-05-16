# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
from matplotlib import style
style.use('fivethirtyeight')
import numpy as np
import pandas as pd
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
p1 = 1.0
p2 = 2.0

p_vap_mmHG = 100
p_vap_pasc = p_vap_mmHG * 133.3


dt = 2
rho = 1000.0  #kg/m3
d_in = 1 #inch

#Datos Escalon P1
t0_e1= 500
A_e1 = 0

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
dt2 = 950

#Datos ExpDecreciente P1
t0_exd1 = 500
tau1 = 100
A_exd1 = 0

#Datos ExpDecreciente P2
t0_exd2 = 500
tau2 = 100
A_exd2 = 0

#límites "y". Cambio de escala
escalapresion = (0,10)
escalaF = (0,8)



def ejercicio1():
    
    #conversión automática de datos
    d_m = d_in * 0.0254
    area = np.pi * d_m**2 / 4
    
    
    T = []
    P1 = []
    P2 = []
    H = []
    F = []
    ANPA = []

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
        
        h = (p2_pascal - p1_pascal)/(9.8 * rho)
        #En lugar de despejar v, resuelvo para F0=0
       
    
        H = np.append(H, h)
    
        f = bomba1(h)
        F = np. append(F, f)
        
        v = f / (area * 3600)  # m/s
        
        hs = v**2 / (2 * 9.8) + p1_pascal / (rho * 9.8)
        
        ANPAi = hs - p_vap_pasc / (9.8 * rho)
        ANPA = np.append(ANPA, ANPAi)
        
        ANPA_bar = ANPA * rho * 9.8 / 100000
        P1_cavit = P1 + ANPA_bar


    #creacion de tabla df para pandas:
    valores = np.vstack((T, P1, P2, H, F))
    valores = valores.T
    columnas = ['tiempo_s', 'p1_bar', 'p2_bar', 'H_m', 'F_m3/h']
    df = pd.DataFrame(valores, columns=columnas)
   
    #gráficos
    fig, ax1 = plt.subplots()

    ax1.plot(T, P1, label="P1")
    ax1.plot(T, P2, label="P2")
    ax1.plot(T, P1_cavit, label="P1 cavitac")
    ax1.plot(T, ANPA, label = "ANPA")
    ax1.set_xlabel('$Tiempo (s)$')
    ax1.set_ylabel('$Presion (bar)$')

    ax2 = ax1.twinx()
    ax2.plot(T, F, 'g-')
    ax2.set_ylabel('$F (m^3/h)$', color='g')
    ax2.tick_params('y', colors='g')

    #limites en y (Corrección de escala)
    ax1.set_ylim(escalapresion)
    ax2.set_ylim(escalaF)
    ax1.legend()

    tabla = df[(df.tiempo_s == 450) | (df.tiempo_s == 600) | (df.tiempo_s == 1000)]
    print tabla
    
    plt.show()
    
ejercicio1()