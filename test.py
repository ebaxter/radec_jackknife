import numpy as np
import pdb 

import radecJackknife as jk


num_rand = 10000
#Generate random catalogs with ra \in (0,90) and dec \in (0,90)
ra_random = np.random.rand(num_rand)*90.
dec_random = np.random.rand(num_rand)*90.

num_data = 100
ra_data = np.random.rand(num_rand)*90
dec_data = np.random.rand(num_rand)*90

#Desired number of jackknife regions (output number of jackknife regions might not match this exactly)
N_jack = 20
#Supply minimum ra of survey 
min_ra = -5.
#Create jackknife object
jackObject = jk.radecJackknife(min_ra)
tol = 0.01
max_iter = 100
#Create regions
jackObject.generate_regions(ra_random, dec_random, N_jack, tol, max_iter)

#Assign data points to regions
data_labels = jackObject.label_pts(ra_data, dec_data)

#Create figure to display points assigned to jackknife regions
import matplotlib.pyplot as pl
fig, ax = pl.subplots(2,1)
ax[0].plot(ra_random, dec_random, 'bs')
ax[1].scatter(ra_data, dec_data, c = data_labels)

pdb.set_trace()




