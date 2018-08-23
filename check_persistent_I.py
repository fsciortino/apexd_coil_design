import numpy as np
import matplotlib.pyplot as plt
import coil_field

from coil_functions import *
from flux_surf_volume import get_flux_vol
plt.ion()
plt.style.use('ggplot')
try:
    import nlopt
except:
    raise ValueError('nlopt not available. Use Python virtualenv!')


Bmax=43.8*1.0e-3 #43.8 milli Tesla

Bz_val = Bz(Ic,Rc,Zc,r,z)
X,Z,Bxm,Bzm,Bs,rAm = multicoil_fields([If],[zf],[rf],[tzf],[trf])
xidx=np.argmin(np.abs(X[0,:]-rf))
zidx=np.argmin(np.abs(Z[:,0]-zf))
rAtheta_obj=rAm[xidx,zidx]
    
