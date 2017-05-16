#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import matplotlib.pyplot as plt
from matplotlib import style
style.use('fivethirtyeight')
import numpy as np
import pandas as pd
#from scipy.optimize import fsolve
#from IPython.display impor display
from SPQentrada import *

#Datos iniciales
f1 = 5.0 #m3/h
f2 = 5.0 #m3/h
altura0 = 4

D_m = 1.5 #metros
d_in = 1.0 #pulgada
L = 10.0
fd = 0.03

rho = 1000  #kg/m3

dt = 2

#Datos Escalon F1
t0_e1= 500
A_e1 = 0

#Datos Escalon F2
t0_e2= 500
A_e2 = 1

##Datos Escalon valvula
#t0_ev= 500
#A_ev = 0

#Datos Rampa F1
t0_r1 = 500
pend1 = 0
dt1 = 500

#Datos Rampa F2
t0_r2 = 500
pend2 = 0
dt2 = 500

##Datos Rampa valvula
#t0_rv = 500
#pendv = 0
#dtv = 2000

#Datos ExpDecreciente F1
t0_exd1 = 500
tau1 = 100
A_exd1 = 0

#Datos ExpDecreciente F2
t0_exd2 = 500
tau2 = 100
A_exd2 = 0

##Datos ExpDecreciente valvula
#t0_exdv = 500
#tauv = 100
#A_exdv = 0

#límites "y". Cambio de escala
escalah = (3,7)
escalaF = (3,7)
#escalaVP = (0,1)



def ejercicio1():
    
    #Autocalculados

    d_m = d_in * 0.0254
    area = np.pi * d_m**2 / 4
    AREA = np.pi * D_m**2 / 4

    
    N = range(0,1001,1)
    
    t = np.zeros(len(N))
    F1 = np.zeros(len(N))
    F2 = np.zeros(len(N))
    F1_m3s = np.zeros(len(N))
    F2_m3s = np.zeros(len(N))
    h = np.zeros(len(N))
    dhdt = np.zeros(len(N))
    
    
    
    h[0] = altura0
    
    dhdt[0] = 0

    for i in range(len(N)):
        t[i] = i * dt
        F1[i] = f1 \
        + ESCALON(t0_e1, A_e1, t[i])\
        + RAMPA(t0_r1, pend1, dt1, t[i])\
        + ExpDecr(t0_exd1, tau1, A_exd1, t[i])
    
        F2[i] = f2 \
        + ESCALON(t0_e2, A_e2, t[i])\
        + RAMPA(t0_r2, pend2, dt2, t[i])\
        + ExpDecr(t0_exd2, tau2, A_exd2, t[i])
        
        
    
        F1_m3s[i] = F1[i] / 3600.0
        F2_m3s[i] = F2[i] / 3600.0
        
        if i < (len(N)-1):
            
            y = h[i] + dhdt[i] * dt
                 
            h[i+1] = y
            
            dydt = (F1_m3s[i] - F2_m3s[i]) / AREA
            
            dhdt[i+1] = dydt
        
        #print"F1=",F1[i]," F2=", F2[i], " h=", h[i]
        
        

    #creacion de tabla df para pandas:
    valores = {'tiempo_s':t,
               'F1_m3/s':F1_m3s, 'F1_m3/h':F1,
               'F2_m3/s':F2_m3s, 'F2_m3/h':F2,
               'altura':h}
    columnas = ['tiempo_s', 'F1_m3/s', 'F2_m3/s',
                'F1_m3/h', 'F2_m3/h', 'altura']
    df = pd.DataFrame(valores, columns=columnas)

    #gráficos
#    plt.plot(t, F1, label='F1')
#    plt.plot(t, F2, label='F2')
#    plt.legend()
#    print F1
    fig = plt.figure(figsize=(8,5))
#    
    ax1 = plt.subplot2grid((6,1),(0,0), rowspan=4)
#    #axv = plt.subplot2grid((6,1),(4,0), rowspan=2)
#    
    ax1.plot(t, F1, label="F1")
    ax1.plot(t, F2, label="F2")
    ax1.set_ylabel('$Caudal\;(m^3/h)$')
    ax1.set_xticklabels([])

    ax2 = ax1.twinx()
    ax2.plot(t, h, 'g-')
    ax2.set_ylabel('$altura\;(m)$', color='g')
    ax2.tick_params('y', colors='g')
    
    #limites en y (Corrección de escala)
    ax1.set_ylim(escalaF)
    ax2.set_ylim(escalah)
    ax1.legend()
    
#    axv.plot(T, VP, 'm-')
#    axv.set_ylabel('$apert$ $valvula$')
#    axv.set_xlabel('$Tiempo (s)$')
#    axv.set_ylim(escalaVP)

    tabla = df[(df.tiempo_s == 450) | (df.tiempo_s == 600) | (df.tiempo_s == 1000)]
    print (tabla)
    
    plt.show()
    
ejercicio1()