import numpy as np
from amuse.units import units
from amuse.units import constants
from amuse.datamodel import Particles
from amuse.datamodel import Particle

# Implementation of class Binary

class Binary:
	def __init__(self, p1 = Particle(), p2 = Particle()):
		self.p1 = p1
		self.p2 = p2
		self.r = p1.position - p2.position
		self.v = p1.velocity - p2.velocity
		self.M = p1.mass + p2.mass
		self.h = np.cross(self.r,self.v)
		self.epsi = (np.dot(self.v,self.v)/2) - ((constants.G * self.M)/np.linalg.norm(self.r))
		self.a = (- constants.G * self.M / (2 * self.epsi)).as_quantity_in(units.parsec)
		self.e = np.sqrt(1 + (2 * self.epsi * np.dot(self.h, self.h))/(constants.G * self.M) ** 2)
	def if_bound(self):
		return self.a > 0 | units.AU

# Non-member function of class Binary

def CoM (binary):
	newsys = Particle()
	bina = Particles()
	bina.add_particle(binary.p1)
	bina.add_particle(binary.p2)
	newsys.position = bina.center_of_mass()
	newsys.velocity = bina.center_of_mass_velocity()
	newsys.mass = bina.mass.sum()
	return newsys
	
# Auxiliary functions of class Sample

def bi_in_tri(bin12, sin, nbi, ns):
	triple = [bin12.p1, bin12.p2, sin]
	ntri = [nbi[0], nbi[1], ns]
	b_in_t = bin12
	nb_in_t = [nbi, ns]
	for i in range(3):
		for j in range(i+1, 3):
			k = np.array(range(3))
			k = k[np.logical_and(k!=j,k!=i)][0]
			btemp = Binary(triple[i],triple[j])
			nbtemp = [[ntri[i], ntri[j]], ntri[k]]
			tritemp = Binary(CoM(btemp), triple[k])
				
			if (btemp.if_bound() and tritemp.if_bound() and btemp.a < b_in_t.a):
				b_in_t = btemp
				nb_in_t = nbtemp
	return nb_in_t

def degree_tri(bi, combi_i,nbin,n1,n2):
	t_bs1 = Binary(combi_i[0], combi_i[1])
	t_bs2 = Binary(combi_i[0], combi_i[2])
	bin_ss = Binary(combi_i[1], combi_i[2])

	ytrip1 = t_bs1.if_bound()
	ytrip2 = t_bs2.if_bound()
	ybin2 = bin_ss.if_bound()
	
	if ytrip1:
		# Triple: (nbin[0], nbin[1], n1)
		comtri = Binary(CoM(t_bs1), combi_i[2])
		yquad = comtri.if_bound()
		tri_conf = bi_in_tri(bi, combi_i[1], nbin, n1)
		if yquad:
			# ***** Status: Quadruple *****
			return [42, [tri_conf,n2], 0, 0, 0]
		else:
			# ***** Status: Triple and a Single ***** 
			return [31, tri_conf, n2, 0, 0]
	elif ytrip2:
		# Triple: (nbin[0], nbin[1], n2)
		comtri = Binary(CoM(t_bs2), combi_i[1])
		yquad = comtri.if_bound()
		tri_conf = bi_in_tri(bi, combi_i[2], nbin, n2)
		if yquad:
			# ***** Status: Quadruple *****
			return [41, [tri_conf,n1], 0, 0, 0]
		else:
			# ***** Status: Triple and a Single *****
			return [30, tri_conf, n1, 0, 0]
	elif ybin2:
		bin_bin = Binary(combi_i[0], CoM(bin_ss))
		yquad = bin_bin.if_bound()
		if yquad:
			# ***** Status: Quadruple ***** 
			return [40, [nbin, [n1,n2]], 0 ,0 ,0]
		else: 
			# ***** Status: Two Binaries ***** 
			return [21, nbin, [n1,n2], 0, 0]
	else:
		# ***** Status: Binary and two Singles *****
		return [20, nbin, n1, n2, 0]


def stability(inner, outer):
	lhs = outer.a / inner.a
	q_out = inner.M / (outer.p2).mass
	i_m = np.arccos(np.dot(outer.h,inner.h) / (np.linalg.norm(outer.h) * np.linalg.norm(inner.h)))
	rhs = (2.8 / (1-outer.e)) * ( ((1+q_out) * (1+outer.e) / (1-outer.e) ** 0.5) ** (0.4) ) * (1 - 0.3 * (i_m / np.pi))
	return (lhs<rhs) #if unstable
		
# Implementation of class Sample			

class Sample:
	def __init__(self, num, pset):
		self.num = num
		self.pset = pset
	
	def degree(self):
		deg_of_combi = []
		
		for i in range(len(self.pset)):
			
			for j in range(i+1, len(self.pset)):
			
				combi = Particles()
				
				# Grouping a binary
				Bin_ij = Binary(self.pset[i],self.pset[j])
				
				
				# Assigning the remaining single particles
				k = range(len(self.pset))
				k.remove(i)
				k.remove(j)
				
				Sin = []
				
				for l in range(len(k)):
					Sin.append(self.pset[k[l]])
				
				# Check if the binary is unbound
				ybin = Bin_ij.if_bound()
				
				if ybin: # bound -> binary = T
					combi.add_particle(CoM(Bin_ij))
					for x in range(len(Sin)):
						combi.add_particle(Sin[x])
					
					
					nbin = [i, j]
					d = degree_tri(Bin_ij, combi,nbin,k[0],k[1])
					deg_of_combi.append(d)
					
					
		if len(deg_of_combi) != 0:
			# Select Maximum Degree
			deg = deg_of_combi[0][0]
			max_deg = deg_of_combi[0]
			for n in range(len(deg_of_combi)):
				if deg < deg_of_combi[n][0]:
					deg = deg_of_combi[n][0]
					max_deg = deg_of_combi[n]
			return max_deg
			
			
		else:
			# ***** Status: Four Singles *****
			deg = [10, 0, 1, 2, 3]
			return deg
			
	def analysis(self):
		d = self.degree()
		if d[0] == 10:
			return [self.num, 's', 'NA', 'NA', 's', 'NA', 'NA', 's', 'NA', 'NA', '(%d, %d, %d, %d)' % (d[0]+1,d[1]+1,d[2]+1,d[3]+1), False, 0]
		elif d[0] == 20:
			b = Binary(self.pset[d[1][0]], self.pset[d[1][1]])
			return [self.num, 'b', (b.a).value_in(units.parsec), b.e, 's', 'NA', 'NA', 's', 'NA', 'NA', '([%d, %d], %d, %d)' % (d[1][0]+1,d[1][1]+1,d[2]+1,d[3]+1), False, 0]
		elif d[0] == 21:
			bin1 = Binary(self.pset[d[1][0]], self.pset[d[1][1]])
			bin2 = Binary(self.pset[d[2][0]], self.pset[d[2][1]])
			unbound12 = Binary(CoM(bin1), CoM(bin2))
			return [self.num, 'b', (bin1.a).value_in(units.parsec), bin1.e, 'b', (bin2.a).value_in(units.parsec), bin2.e, 'u', (unbound12.a).value_in(units.parsec), unbound12.e, '([%d, %d], [%d, %d])' % (d[1][0]+1,d[1][1]+1,d[2][0]+1,d[2][1]+1), False, 0]
		elif d[0] >= 30 and d[0] < 40:
			b_in_t = Binary(self.pset[d[1][0][0]], self.pset[d[1][0][1]])
			tri = Binary(CoM(b_in_t), self.pset[d[1][1]])
			unboundts = Binary(CoM(tri), self.pset[d[2]])
			return [self.num, 'b', (b_in_t.a).value_in(units.parsec), b_in_t.e, 't', (tri.a).value_in(units.parsec), tri.e, 'u', (unboundts.a).value_in(units.parsec), unboundts.e, '([[%d, %d], %d], %d)' % (d[1][0][0]+1, d[1][0][1]+1, d[1][1]+1, d[2]+1), stability(b_in_t, tri), 0]
		elif d[0] == 40:
			b1 = Binary(self.pset[d[1][0][0]],self.pset[d[1][0][1]]) 
			b2 = Binary(self.pset[d[1][1][0]],self.pset[d[1][1][1]])
			q = Binary(CoM(b1), CoM(b2))
			return [self.num, 'b', (b1.a).value_in(units.parsec), b1.e, 'b', (b2.a).value_in(units.parsec), b2.e, 'q', (q.a).value_in(units.parsec), q.e, '[[%d, %d], [%d, %d]]' % (d[1][0][0]+1, d[1][0][1]+1, d[1][1][0]+1, d[1][1][1]+1), (stability(b1, q) and stability(b2, q)), 0]
		elif d[0] >= 41:
			b = Binary(self.pset[d[1][0][0][0]], self.pset[d[1][0][0][1]])
			tri = Binary(CoM(b), self.pset[d[1][0][1]])
			q = Binary(CoM(tri), self.pset[d[1][1]])
			return [self.num, 'b', (b.a).value_in(units.parsec), b.e, 't', (tri.a).value_in(units.parsec), tri.e, 'q', (q.a).value_in(units.parsec), q.e, '[[[%d, %d], %d], %d]' % (d[1][0][0][0]+1, d[1][0][0][1]+1, d[1][0][1]+1, d[1][1]+1), (stability(b,tri) and stability(tri, q)), 0] 
	
	# output format:
	# sample_no	type_bin1 a1	e1	type_bin2 a2	e2	type_bin3	a3	e3	configuration	collision_flag
	# type_bin = b(single-single), t(binary-single), q(bin-bin or tri-sin), u(unbound) and s(single) 
	
			
		

