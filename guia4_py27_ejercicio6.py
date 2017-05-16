import matplotlib.pyplot as plt
from matplotlib import style
style.use('fivethirtyeight')
import numpy as np
from scipy.optimize import fsolve
import pandas as pd
#from IPython.display impor display
from SPQentrada import *

#######revisar

#Datos iniciales
p1 = 2.0
p2 = 1.0
#p3 = 1.0

L1 = 5.0
L2 = 5.0
L3 = 10.0
L4 = 5.0
L5 = 5.0

d_in = 1
fd = 0.03
rho = 1000.0  #kg/m3

dt = 2

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
pend2 = 0
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
escalapresion = (0,4)
escalaF = (0,8)



def ejercicio1():
    
    #conversión automática de datos
    d_m = d_in * 0.0254
    area = np.pi * d_m**2 / 4
    
    #Guardado de variables
    T = []
    P1 = []
    P2 = []
    #P3 = []
    V1 = []
    V2 = []
    V3 = []
    F1 = []
    F2 = []
    F3 = []

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
        
#        p3t = p2 \
#        + ESCALON(t0_e2, A_e2, t)\
#        + RAMPA(t0_r2, pend2, dt2, t)\
#        + ExpDecr(t0_exd2, tau2, A_exd2, t)
#        P3 = np.append(P3, p3t)
    
        p1_pascal = p1t * 100000
        p2_pascal = p2t * 100000
#        p3_pascal = p2t * 100000
        
        fricc1 = fd * L1/d_m
        fricc2 = fd * L2/d_m
        fricc3 = fd * L3/d_m
        fricc4 = fd * L4/d_m
        fricc5 = fd * L5/d_m
            
        #deltaP_h = (p1_pascal - p2_pascal)/rho #g se simplifico y anulo
    

#%%
                    
        def FOsolv (v1):
            
            
#            v2 = np.sqrt(((p1_pascal-p2_pascal)- fricc1 * v1**2 - fricc5 * v1**2)\
#                  / (2*fricc2 + fricc4))
#            
#            v3 = np.sqrt(((2*fricc2 + fricc4) * v2**2)\
#                  / (2*fricc2 + fricc4))

            Numerador = 2 * ((p1_pascal - p2_pascal)/rho\
                             - v1**2/2 * (fricc1 + fricc5)) 
           
            v2 = np.sqrt(Numerador / (2 * fricc2 + fricc4))
            
            v3 = np.sqrt(Numerador / (2 * fricc3 + fricc4))
                    
            
            F0 = v1 - v2 - v3
            
            return F0
        
        semilla = 1.5
    
        sol = fsolve(FOsolv, semilla)
        
        v1 = sol
        
        Numerador = 2 * ((p1_pascal - p2_pascal)/rho\
                             - v1**2/2 * (fricc1 + fricc5)) 
           
        v2 = np.sqrt(Numerador / (2 * fricc2 + fricc4))
            
        v3 = np.sqrt(Numerador / (2 * fricc3 + fricc4))
        
        
    
        V1 = np. append(V1, v1)
        f1 = v1*area
        F1 = np. append(F1, f1)
        
        V2 = np. append(V2, v2)
        f2 = v2*area
        F2 = np. append(F2, f2)
        
        V3 = np. append(V3, v3)
        f3 = v3*area
        F3 = np. append(F3, f3)

    F1_m3h= F1 * 3600
    
    F2_m3h= F2 * 3600
    
    F3_m3h= F3 * 3600
       

    #creacion de tabla df para pandas:
    valores = np.vstack((T, P1, P2, V1, V2, V3, F1_m3h, F2_m3h, F3_m3h))
    valores = valores.T
    columnas = ['tiempo_s', 'p1_bar', 'p2_bar',
    'V1', 'V2', 'V3', 'F1_m3/h','F2_m3/h','F3_m3/h']
    df = pd.DataFrame(valores, columns=columnas)

    #gráficos
    fig, ax1 = plt.subplots()

    ax1.plot(T, P1, label="P1")
    ax1.plot(T, P2, label="P2")
    ax1.set_xlabel('$Tiempo (s)$')
    ax1.set_ylabel('$Presion (bar)$')

    ax2 = ax1.twinx()
    ax2.plot(T, F1_m3h, 'g-')
    ax2.plot(T, F2_m3h, 'g-')
    ax2.plot(T, F3_m3h, 'g-')
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