import numpy

def ESCALON (t0, A, t):
    if t <= t0:
        y = 0
    else:
        y = A
    return y  



def RAMPA (t0, pend, dt, t):
    if t <= t0:
        y = 0
    elif t < t0 + dt:
        y = pend * (t - t0)
    else:
        y = pend * dt
    return y


def SENO (t0, periodo, A, dt, t):
    if t <= t0:
        SENO = 0
    elif t <= t0 + dt:
        SENO = A * numpy.sin(2*numpy.pi / periodo * (t - t0))
    else:
        SENO = 0
    return SENO


def ExpDecr (t0, tau, A, t):
    """Exponencial Decreciente"""
    if t <= t0:
        Y = 0
    else:
        Y = A * (1 - numpy.exp(-(t - t0) / tau))
    return Y


def DobleRampa(A, t0_rampa, dt_r1, tf_plano, dt_r2, t):
    if t <= t0_rampa:
        DobleRampa = 0
    elif t < t0_rampa + dt_r1:
        DobleRampa = (A / dt_r1) * (t - t0_rampa)
    elif t <= tf_plano:
        DobleRampa = A
    elif t <= tf_plano + dt_r2:
        DobleRampa = A + (-1 * A / dt_r2) * (t - tf_plano)
    else:
        DobleRampa = 0
    return DobleRampa