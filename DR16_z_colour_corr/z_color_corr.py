import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from matplotlib.patches import Rectangle
from matplotlib         import colors as mcolors

from astropy.io import ascii
from astropy.io import fits

import sys
import warnings

if not sys.warnoptions:
    warnings.simplefilter("ignore")


path     = 'DR16/'
filename = 'DR16Q_v4.fits'


hdul = fits.open(path+filename)  # open a FITS file
data = hdul[1].data  # assume the first extension is a table
# print(hdul[0].header)
## Putting the SDSS 5-bands into separate arrays
uband_psfmag = data['PSFMAG'][:,0]
gband_psfmag = data['PSFMAG'][:,1]
rband_psfmag = data['PSFMAG'][:,2]
iband_psfmag = data['PSFMAG'][:,3]
zband_psfmag = data['PSFMAG'][:,4]
Yband_flux = data['YFLUX']

## Selecting all non-zero optical data
optical_data = data[(uband_psfmag>0) &
                  (gband_psfmag>0) &
                  (rband_psfmag>0) &
                  (iband_psfmag>0) &
                  (zband_psfmag>0) &
                  (Yband_flux>0)]


# ## Selecting all non-zero Y-band data
# Yband_data = optical_data['YFLUX']



# convert Yband flux to AB mag
optical_data['Yflux'] = - 2.5*np.log10(1e3*optical_data['Yflux']) - 48.6
# print(optical_data['Yflux'][:20])

# plt.hist(optical_data['Z'], bins = 200)
# plt.show()


## Settting up the min/max redshift, the bin width and the number of bins
zmin = 0.00
zmax = 5.00
delta_zbin = 0.05
nbins = (zmax-zmin)/delta_zbin

## Setting up the colour arrays   
u_minus_g = np.zeros(int(nbins))
g_minus_r = np.zeros(int(nbins))
r_minus_i = np.zeros(int(nbins))
i_minus_z = np.zeros(int(nbins))
z_minus_y = np.zeros(int(nbins))

z_lo_list = []
z_hi_list = []
## Looping over each redshift bin, calculating the mean colour in each bin
## number of quasars in each bin also reported
print('z-bin, (u-g), (g-r), (r-i), (i-z), (z-y), N_Q')

for ii in range(int(nbins)):
  z_lo =  round(ii    * delta_zbin,3)
  z_hi = round((ii+1) * delta_zbin,3)
  z_lo_list.append(z_lo)
  z_hi_list.append(z_hi)
  zcut_optical_data = optical_data[(optical_data['Z'] >= z_lo) &  (optical_data['Z'] < z_hi)]

#    print(ii, z_lo, z_hi, len( zcut_optical_data))
  u_minus_g[ii]  = round(np.mean(zcut_optical_data['PSFMAG'][:,0] - zcut_optical_data['PSFMAG'][:,1]),3)
  g_minus_r[ii]  = round(np.mean(zcut_optical_data['PSFMAG'][:,1] - zcut_optical_data['PSFMAG'][:,2]),3)
  r_minus_i[ii]  = round(np.mean(zcut_optical_data['PSFMAG'][:,2] - zcut_optical_data['PSFMAG'][:,3]),3)
  i_minus_z[ii]  = round(np.mean(zcut_optical_data['PSFMAG'][:,3] - zcut_optical_data['PSFMAG'][:,4]),3)
  z_minus_y[ii]  = round(np.mean(zcut_optical_data['PSFMAG'][:,4] - zcut_optical_data['Yflux']),3)

  # print('{0:3f} {1:4f} {2:4f} {3:4f} {4:4f} {5:4f} {6:4f}'.format(((z_lo+z_hi)/2.), u_minus_g[ii], g_minus_r[ii], r_minus_i[ii], i_minus_z[ii], z_minus_y[ii], len( zcut_optical_data)))


colour_redshift_corr = pd.DataFrame(list(zip(z_lo_list, z_hi_list, u_minus_g,  g_minus_r,  r_minus_i,  i_minus_z, z_minus_y)), 
  columns = ['low_z','high_z', 'u-g', 'g-r', 'r-i', 'i-z', 'z-y'])

colour_redshift_corr.to_csv('colour_redshift_corr.csv')

  