from bin_gen import *
from amuse.units import units
from amuse.units import constants
from amuse.units import nbody_system
from amuse.datamodel import Particles
from amuse.datamodel import Particle
import numpy as np

# OldSystem(Movement of Binaries in Center of Mass) Initialization

def CoM_System (K, M1, M2, mu1, mu2, a1, a2, e1, e2, d):
	GM1 = constants.G * M1
	GM2 = constants.G * M2
	E_b_1 = GM1 * mu1 / (2 * a1)
	E_b_2 = GM2 * mu2 / (2 * a2) 
	T = K * (E_b_1 + E_b_2) # get rid of K
	v2 = np.sqrt(2 * T * M1 / (M2 * (M1 + M2)))
	v1 = (M2/M1) * v2
	OldSys = Particles(2)
	OldSys.x = [-(M2 / (M1 + M2))*d, (M1 / (M1 + M2))*d]
	OldSys.y = [0,0] | units.AU
	OldSys.z = [0,0] | units.AU
	OldSys.mass = [M1, M2]
	OldSys.vx = [v1, -v2]
	OldSys.vy = [0,0] | units.AU/units.s
	OldSys.vz = [0,0] | units.AU/units.s # use impact parameter to determine velocity
	return OldSys

# Check K
def check_k(bin1, bin2, a1, a2):
	v1 = bin1.center_of_mass_velocity
	v2 = bin2.center_of_mass_velocity
	M1 = bin1.mass.sum()
	M2 = bin2.mass.sum()
	mu1 = (bin1[0].mass * bin1[1].mass)/bin1.mass.sum()
	mu2 = (bin2[0].mass * bin2[1].mass)/bin2.mass.sum()
	T1 = 0.5 * M1 * (v1 ** 2)
	T2 = 0.5 * M2 * (v2 ** 2)
	E_b_1 = constants.G * M1 * mu1 / (2 * a1)
	E_b_2 = constants.G * M2 * mu2 / (2 * a2)
	K = (T1 + T2)/(E_b_1 + E_b_2)
	return K

# 2b-CoM conversion function

def expand_body(newsys, oldcom, syst):

	for mem in newsys:
		mem.position += oldcom.position
		mem.velocity += oldcom.velocity
       
	## removing fictional particle
	## first check if m1+m2 = p.mas
	if(newsys.mass.sum() != oldcom.mass):
		raise ValueError("components mass not equal to center of mass\n{:s}!={:s}".format(newsys.mass.sum(), oldcom.mass))
	syst.remove_particle(oldcom)
	syst.add_particles(newsys)
   
	return syst

# 2*2b->4b conversion

def four_b_i(B1, B2, comSys): #comSys should consist two P's
	B1.move_to_center()
	B2.move_to_center()
	com1 = comSys[0]
	com2 = comSys[1]
	expand_body(B1, com1, comSys)
	expand_body(B2, com2, comSys)
	return comSys

# Output Function
def write_mar_input_tides(p, filename, mass_unit=units.MSun, length_unit=units.parsec, silent=False):
    """
    p is particle dataset with:
    x
    y
    z
    vx
    vy
    vz
    mass
    radius
    polytropic_index
    eps
    """
    conv=nbody_system.nbody_to_si(mass_unit, length_unit)
    velocity_unit = conv.to_si(nbody_system.length / nbody_system.time).as_unit()
   
    np.savetxt(filename, np.transpose([p.x.value_in(length_unit),
                                       p.y.value_in(length_unit),
                                       p.z.value_in(length_unit),
                                       p.vx.value_in(velocity_unit),
                                       p.vy.value_in(velocity_unit),
                                       p.vz.value_in(velocity_unit),
                                       p.mass.value_in(mass_unit),
                                       p.eps.value_in(length_unit),
                                       p.radius.value_in(length_unit),
                                       p.polytropic_index]),
               fmt='%.13g %.13g %.13g %.13g %.13g %.13g %.13g %.13g %.13g %.13g')
    if(not silent):
        print("written {:s}\n".format(filename))
        print("1 time unit = {:s}".format(conv.to_si(nbody_system.time).as_quantity_in(units.Myr)))

# Nbody units is relative to G = 1 -> convert into si.
# G = length ** 3 / mass * time ** 2 


# Schwarzschild Radius conversion

def schwarzschild(M):
	r = 2 * constants.G * M/ (constants.c ** 2)
	r.as_quantity_in(units.RSun)
	return r

# Main

sampleSize = 1

for n in range(sampleSize):
	b1 = Pset_n[0]
	b2 = Pset_n[1]
	ComSys = CoM_System(0, M_tot[0], M_tot[1], redu_m[0], redu_m[1], a, a, e, e, d)
	FBI = four_b_i(b1, b2, ComSys)
	FBI.eps = 0 | units.RSun
	FBI.radius = schwarzschild(FBI.mass)
	FBI.polytropic_index = 0 
	write_mar_input_tides(FBI, "4b/%04d.txt" % (n,))
	T_crossing = constants.G * FBI.mass.sum() ** 2.5/ (2* FBI.kinetic_energy() ** 1.5)
	print T_crossing.as_quantity_in(units.yr)
	Pset_n.remove(Pset_n[0])
	Pset_n.remove(Pset_n[1])
