import numpy as np
from amuse.datamodel import *

MIN_ERR = 1e-15

"""
Takes:
 - secondary: particle (required attributes: mass)
 - primary: particle (required attributes: mass, 3d position, 3d velocity)
 - a: semimajor axis (in units of length)
 - e: eccentricity
 - inc: inclination (radians)
 - omega: argument of pericenter (radians)
 - Omega: longitude of ascending node (radians)
 - f: true anomaly (radians)
""" 

def kepl_to_cart(secondary, primary, a, e, inc, omega, Omega, f):   
    #ERROR DEALING
    if(e == 1.):
        raise ValueError("Can't initialize a radial orbit (e=1) with orbital elements")
    elif(e < 0.):
        raise ValueError("Eccentricity must be greater than or equal to zero")
    elif(e > 1.):
        if(a.number > 0.):
            raise ValueError("Bound orbit (a > 0) must have e < 1")
    else:
        if(a.number < 0):
            raise ValueError("Unbound orbit (a < 0) must have e > 1")
    if(e*np.cos(f)< -1.):
        raise ValueError("Unbound orbit can't have f set beyond the range allowed by the asymptotes set by the parabola.")
    
    if(e < MIN_ERR):
        omega = 0.
    #if((inc < MIN_INC) or (inc > np.pi-MIN_INC)):
    #    Omega = 0.
    
    mu = constants.G*(secondary.mass+primary.mass)
    r = a*(1-e*e)/(1 + e*np.cos(f))
    v0 = (mu/a/(1.-e*e)).sqrt() # in this form it works for elliptical and hyperbolic orbits

    cO = np.cos(Omega)
    sO = np.sin(Omega)
    co = np.cos(omega)
    so = np.sin(omega)
    cf = np.cos(f)
    sf = np.sin(f)
    ci = np.cos(inc)
    si = np.sin(inc)
    
    # Murray & Dermott Eq 2.122
    secondary.x = primary.x + r*(cO*(co*cf-so*sf) - sO*(so*cf+co*sf)*ci)
    secondary.y = primary.y + r*(sO*(co*cf-so*sf) + cO*(so*cf+co*sf)*ci)
    secondary.z = primary.z + r*(so*cf+co*sf)*si

    # Murray & Dermott Eq. 2.36 after applying the 3 rotation matrices from Sec. 2.8 to the velocities in the orbital plane
    secondary.vx = primary.vx + v0*((e+cf)*(-ci*co*sO - cO*so) - sf*(co*cO - ci*so*sO))
    secondary.vy = primary.vy + v0*((e+cf)*(ci*co*cO - sO*so)  - sf*(co*sO + ci*so*cO))
    secondary.vz = primary.vz + v0*((e+cf)*co*si - sf*si*so)

    return secondary

def kepl_to_cart2(secondary, primary, a, e, inc, omega, Omega, r):   
    #ERROR DEALING
    if(e == 1.):
        raise ValueError("Can't initialize a radial orbit (e=1) with orbital elements")
    elif(e < 0.):
        raise ValueError("Eccentricity must be greater than or equal to zero")
    elif(e > 1.):
        if(a.number > 0.):
            raise ValueError("Bound orbit (a > 0) must have e < 1")
    else:
        if(a.number < 0):
            raise ValueError("Unbound orbit (a < 0) must have e > 1")
    if(e*np.cos(f)< -1.):
        raise ValueError("Unbound orbit can't have f set beyond the range allowed by the asymptotes set by the parabola.")
    
    if(e < MIN_ERR):
        omega = 0.
    #if((inc < MIN_INC) or (inc > np.pi-MIN_INC)):
    #    Omega = 0.
    
    mu = constants.G*(secondary.mass+primary.mass)
    r = a*(1-e*e)/(1 + e*np.cos(f))

    f = 
    v0 = (mu/a/(1.-e*e)).sqrt() # in this form it works for elliptical and hyperbolic orbits

    cO = np.cos(Omega)
    sO = np.sin(Omega)
    co = np.cos(omega)
    so = np.sin(omega)
    cf = np.cos(f)
    sf = np.sin(f)
    ci = np.cos(inc)
    si = np.sin(inc)
    
    # Murray & Dermott Eq 2.122
    secondary.x = primary.x + r*(cO*(co*cf-so*sf) - sO*(so*cf+co*sf)*ci)
    secondary.y = primary.y + r*(sO*(co*cf-so*sf) + cO*(so*cf+co*sf)*ci)
    secondary.z = primary.z + r*(so*cf+co*sf)*si

    # Murray & Dermott Eq. 2.36 after applying the 3 rotation matrices from Sec. 2.8 to the velocities in the orbital plane
    secondary.vx = primary.vx + v0*((e+cf)*(-ci*co*sO - cO*so) - sf*(co*cO - ci*so*sO))
    secondary.vy = primary.vy + v0*((e+cf)*(ci*co*cO - sO*so)  - sf*(co*sO + ci*so*cO))
    secondary.vz = primary.vz + v0*((e+cf)*co*si - sf*si*so)

    return secondary
