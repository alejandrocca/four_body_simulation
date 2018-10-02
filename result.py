import glob
import errno
import numpy as np
import re
import os
from amuse.units import units
from amuse.datamodel import Particles
from amuse.datamodel import Particle


### Analysis Functions

from classBinary import * 	

### Main

import csv
from io import BytesIO

path = '****_output.dat'

files = glob.glob(path)

sample_list = []


# Read in files

for name in files:
	

	num = re.findall('\d+', name)[0]
	col = '{:}_collision.dat'.format(num)
	
	if (not os.path.isfile(col)):
	
		try:
			with open(name) as f:	

				lines = '\n'.join(f.readlines()[-4:])
				data = np.genfromtxt(BytesIO(lines))
				Pset = Particles()
				for x in range(4):
					Pset.add_particle(Particle(x = data[x][0]|units.parsec, y = data[x][1]|units.parsec, z = data[x][2]|units.parsec, vx = data[x][3]|units.parsec/(14.9080299942 * units.Myr), vy = data[x][4]|units.parsec/(14.9080299942 * units.Myr), vz = data[x][5]|units.parsec/(14.9080299942 * units.Myr), mass = data[x][6]|units.MSun))
				sample_i = Sample(name, Pset)
				sample_list.append(sample_i.analysis())
		except IOError as exc:
			if exc.errno != errno.EISDIR:
				raise
			
	
	else: 
		sample_list.append([name, 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 0, 1])
		
# Write out file

with open('table.csv', 'w') as dat:
	writer = csv.writer(dat, delimiter = '\t', lineterminator = '\n')
	writer.writerows(sample_list)

	#writer = csv.writer(csvfile, delimiter='\t')
	#[writer.writerow(r) for r in sample_list]
