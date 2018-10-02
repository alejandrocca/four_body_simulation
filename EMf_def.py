import numpy as np
from scipy.optimize import newton

# Eccentricity anomaly - True anomaly conversion

def f_to_E(fase, e):
    """ NOTE: if we have r, eccentric anomaly is much simpler and inexpensive to compute:
    E = np.arccos((1.-r/a)/e)
    but needs to change sign according to v*r"""
    E = np.arctan2(np.sqrt(1-e**2) * np.sin(fase),( e + np.cos(fase)))
    if(E < 0):
        E = 2*np.pi + E
    #E[E < 0] = 2*np.pi +  E[E < 0]
    return E

def E_to_f(E, e):
    fase = 2 *np.arctan2(np.sqrt(1+e) * np.sin(E/2.),( np.sqrt(1-e) * np.cos(E/2.)))
    return fase

# Eccentricity anomaly - Mean anomaly conversion

def ecc_to_mean(E, M, e):
    x = E - e * np.sin(E) - M
    return x

def ecc_to_mean_prime(E, M, e):
    x = 1. - e * np.cos(E)
    return x

def ecc_to_mean_prime2(E, M, e):
    x = e * np.sin(E)
    return x

def M_to_E(M, e):
    #if(e<0.8):
    #    E0 = M
    #else:
    #    E0 = np.pi
    E0 = M + np.sign(np.sin(M)) * 0.85 * e
    E = newton(ecc_to_mean, E0, fprime=ecc_to_mean_prime, fprime2=ecc_to_mean_prime2, args = (M, e))
    return E

def E_to_M(E, e):
    M = E - e*np.sin(E)
    return M





