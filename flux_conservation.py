import numpy as np
import matplotlib.pyplot as plt
from scipy.special import ellipk, ellipe
import matplotlib.ticker as ticker
from coil_functions import *
from flux_surf_volume import get_flux_vol
import scipy.interpolate as interp
import matplotlib.patches as patches
#from scipy.optimize import fmin_cobyla
plt.ion()
plt.style.use('ggplot') #'dark_background'
import sys
try:
    import nlopt
except:
    raise ValueError('nlopt not available. Use Python virtualenv!')

g = 9.80665
charging_pos=-0.15
Ips=10 # current from power supply
d_wire = 2.0e-3 # wire diameter

# L-coil
nr_l=9; nz_l=47 # fixed during construction. Dependent on bobbin size
num_turns_l = nr_l * nz_l
tz_l=100.0e-3; tr_l=22.0e-3
I_l=nr_l*nz_l*Ips #5.0e3
r_l=0.085 #0.106 #0.08# expected L-coil winding # 0.08 #former test
z_l=0.09 #0.07 #10

# F-coil
nr_f=100  #dense model, default
nz_f=100  #dense model, default
tz_f = 25.0e-3; tr_f = 15.0e-3 #fixed/measured
I_f=3895.32  #value from NIFS test
r_f=0.0445 # Major radius 52mm - 15mm/2.0  #fixed/measured
m_f=0.33 #Includes LN2 0.2845 #fixed/measured
z_f=0.0   #fixed

# C coil 
nr_c=10
nz_c=92
I_c=nr_c*nz_c*Ips #9140 AT
r_c=0.085         # attempt larger than F-coil
z_c = z_f #-0.05 # guess for optimization should go to above
tz_c=200.0e-3; tr_c=20.0e-3

#Estimate reduction in I_f from flux freezing when close to powered L-coil

# Measured value of Bz_max:
Bz_max=32.3e-3
If1 = If_from_Bz_onaxis(Bz_max,r_f,z_f,0,0)
If1 = If1*2.0  #7000
I_l = I_l * 2.0

zarr=np.linspace(0.01, 0.10, 10)
weights=[]
for zzl in zarr:

    # Compute using flux freezing:
    X,Z,Bxm,Bzm,Bs,rAm = multicoil_fields([I_l],[zzl],[r_l],[tz_l],[tr_l],[nz_l],[nr_l])
    xidx=np.argmin(np.abs(X[0,:]-z_f))
    zidx=np.argmin(np.abs(Z[:,0]-r_f+tr_f/2.0))
    rAm_0=rAm[:,xidx]

    X_f,Z_f,Bxm_f,Bzm_f,Bs_f,rAm_f = multicoil_fields([I_f],[z_f],[r_f],[tz_f],[tr_f],[nz_f],[nr_f])
    xidx_f=np.argmin(np.abs(X_f[0,:]-z_f))
    zidx_f=np.argmin(np.abs(Z_f[:,0]-r_f+tr_f/2.0))
    rAm_f=rAm_f[:,xidx_f]

    frac= rAm_0[zidx] / rAm_f[zidx_f]
    If2= frac* If1

    force_fv = coil_force_fv(r_l,r_f,If2,I_l,z_f,zzl,nr_l,nz_l)
    ww = (force_fv-m_f*g)/g
    weights.append(ww)

plt.figure(111)
plt.plot(zarr,weights,'*-')
plt.xlabel(r'$z_l$ [m]')
plt.ylabel('Weight [kg]')
