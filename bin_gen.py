import numpy as np
from amuse.units import units
from amuse.datamodel import Particles
from amuse.datamodel import Particle
from kepl2cart_def import * #kepl_to_cart(secondary, primary, a, e, inc, omega, Omega, f)
from EMf_def import * #transform M to f

## Fixed constants

a = 1.0 | units.AU
e = 0.0
d = 100 * a
m1 = 30.0 | units.MSun
m2 = 30.0 | units.MSun
n = 10000 # Size of samples

## Generate n sets of i, omega, OMEGA

# - omg: degree of rotation on orbital plane

omg_n = np.random.uniform(0, 2*np.pi, n)
np.asarray(omg_n)

# - i: degree of inclination w.r.t reference plane

cos_i = np.random.uniform(-1,1,n)
np.asarray(cos_i)
i_n = np.arccos(cos_i)

# - OMG: degree of rotation w.r.t axis normal to reference plane

OMG_n = np.random.uniform(0, 2*np.pi, n)
np.asarray(OMG_n)

## Generating n phases from randomly selected Mean Anomaly
M_n = np.random.uniform(0, 2*np.pi, n)
np.asarray(M_n)
E_n = np.empty(n)
for i in range(n):
	E_n[i] = M_to_E(M_n[i], e)
f_n = E_to_f(E_n, e)

## Generate n sets of binaries
# - fixed primary particle
primary = Particle(mass = m1, position = [0,0,0] | units.AU, velocity = [0,0,0] | units.AU/units.s)
secondary = Particle(mass = m2)

# - generated secondary particle and formed set
Pset_n = []
Pcenter_n = []
M_tot = []
redu_m = []
for i in range(n):
	M_tot.append(primary.mass + secondary.mass)
	redu_m.append(primary.mass*secondary.mass/(primary.mass + secondary.mass))
	Pset = Particles(0)
	Pset.add_particle(primary)
	secondary_i = kepl_to_cart(secondary, primary, a, e, i_n[i], omg_n[i], OMG_n[i], f_n[i])
	Pset.add_particle(secondary_i)
	Pset_n.append(Pset)
	Pcenter = Particles(0)
	Pcenter.add_particle(primary)
	Pcenter.add_particle(secondary_i)
	Pcenter.move_to_center()
	Pcenter_n.append(Pcenter)

