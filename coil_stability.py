import numpy as np
import matplotlib.pyplot as plt
from scipy.special import ellipk, ellipe
import matplotlib.ticker as ticker
from coil_functions import *
plt.ion()
plt.style.use('ggplot') #'dark_background'

# params
r_f=1.0 #0.25
z_f=0.0                                               
r_l=0.2 #0.4                                               
z_l=4.6#0.612                                               
I_f=2.0e3#250e3                                               
I_l=12.0e3#28.8e3
m_f=1.0#110.0

#r_f=0.25
#r_l=0.6
#z_f=0.0                                               
#z_l=0.312                                               
#I_f=250e3                                               
#I_l=12.8e3
#m_f=16.8

g= 9.80665
                                               
# -------------------------------
force1 = -I_f*I_l *dMdz(r_f,r_l,z_l,z_f) 
freq = 1.0/(2*np.pi) * np.sqrt(-2*np.pi*r_f *I_f* dBrdz(I_l,r_l,z_l,r_f,z_f)/m_f)
                                               
force2=- 2*np.pi*r_f*I_f*Br(I_l,r_l,z_l,r_f,z_f)

force_imbalance=force1 - m_f * g
print force1, freq, force_imbalance

# -------------------------------
# Stability plots

# stability
# fix (r,z) for SC F-coil and vary L-coil parameters
r_l_arr=np.arange(0.05, 5.0, 0.05)
z_l_arr=np.arange(0.05, 5.0, 0.05)
dBrdz_arr=np.zeros((r_l_arr.shape[0],z_l_arr.shape[0]))
Br_arr=np.zeros((r_l_arr.shape[0],z_l_arr.shape[0]))
dBzdr_arr=np.zeros((r_l_arr.shape[0],z_l_arr.shape[0]))
tilt_var_arr=np.zeros((r_l_arr.shape[0],z_l_arr.shape[0]))
for i in range(r_l_arr.shape[0]):
    for j in range(z_l_arr.shape[0]):
        r=r_l_arr[i]; z=z_l_arr[j]
        Br_arr[i,j]=Br(I_l,r,z,r_f,0)
        dBrdz_arr[i,j]=dBrdz(I_l,r,z,r_f,0)
        dBzdr_arr[i,j]=dBzdr(I_l,r,z,r_f,0)
        tilt_var_arr[i,j]=Bz(I_l,r,z,r_f,0)+r_f * dBrdz_arr[i,j]

        
fig=plt.figure(figsize=(8,12))
ax1=plt.subplot(3,1,1)
X,Y=np.meshgrid(r_l_arr,z_l_arr)
levs=np.asarray([round(i,15) for i in np.logspace(np.log10(1e-5), np.log10(1e-2),10)])
levs=np.concatenate((levs,-levs))
levs=[float('{:.2}'.format(i)) for i in levs]
linest=['solid' if x>0 else 'dashed' for x in levs]
cols=['b' if x>0 else 'r' for x in levs]
fmt = '%r'
plt.contour(Y,X,Br_arr,levels=np.linspace(np.abs(max(Br_arr.flatten())),np.abs(min(Br_arr.flatten())),10),linewidth=2,linestyles='dashed',colors='g')
CS=plt.contour(Y,X,dBrdz_arr, levels=levs,linestyles=linest, linewidths=2, colors=cols)
CS.clabel(levs,fmt=fmt)
ax1.set_title('Vertical stability')
ax1.set_ylabel(r'$z_L$ [$r_F$ units]')

ax2=plt.subplot(3,1,2)
#levs=np.asarray([round(i,15) for i in np.logspace(np.log10(1e-5), np.log10(1e-2),10)])
levs=np.asarray([-1e-5,-0.005,-2e-3,1e-3,5e-4,-1e-4,-5e-5, -2e-5, 0.005, 0.0002,0.0001, 5e-5, 2e-5,1e-5])
levs=np.sort(levs)
#levs=np.concatenate((levs,-levs))
levs=[float('{:.2}'.format(i)) for i in levs]
linest=['dashed' if x>0 else 'solid' for x in levs]
cols=['r' if x>0 else 'b' for x in levs]
fmt = '%r' 
CS=plt.contour(Y,X,dBzdr_arr, levels=levs,linestyles=linest, linewidths=2, colors=cols)
CS.clabel(levs,fmt=fmt)
ax2.set_title('Slide stability')
ax2.set_ylabel(r'$z_L$ [$r_F$ units]')

ax3=plt.subplot(3,1,3)
levs=np.asarray([round(i,15) for i in np.logspace(np.log10(1e-5), np.log10(1e-2),10)])
levs=np.concatenate((levs,-levs))
levs=[float('{:.2}'.format(i)) for i in levs]
linest=['solid' if x>0 else 'dashed' for x in levs]
cols=['b' if x>0 else 'r' for x in levs]
fmt = '%r' 
CS=plt.contour(Y,X,tilt_var_arr, levels=levs,linestyles=linest, linewidths=2, colors=cols)
CS.clabel(levs,fmt=fmt)
ax3.set_title('Tilt stability')
ax3.set_ylabel(r'$z_L$ [$r_F$ units]') 
ax3.set_xlabel(r'$r_L$ [$r_F$ units]')

# Total stability map
Bs=np.zeros((r_l_arr.shape[0],z_l_arr.shape[0]))
vert_term=np.zeros((r_l_arr.shape[0],z_l_arr.shape[0]))
slide_term=np.zeros((r_l_arr.shape[0],z_l_arr.shape[0]))
tilt_term=np.zeros((r_l_arr.shape[0],z_l_arr.shape[0]))
selector=np.zeros((r_l_arr.shape[0],z_l_arr.shape[0]))
for i in range(r_l_arr.shape[0]):
    for j in range(z_l_arr.shape[0]):
        r=r_l_arr[i]; z=z_l_arr[j]
        Bs[i,j]=np.sqrt(Br(I_l,r,z,r_f,0)**2+Bz(I_l,r,z,r_f,0)**2) 
        vert_term[i,j] = np.abs(r_f/Br(I_l,r,z,r_f,0)) * dBrdz_arr[i,j]
        slide_term[i,j] = np.abs(r_f/Bz(I_l,r,z,r_f,0)) * dBzdr_arr[i,j]
        #tilt_term[i,j] = np.abs(1.0/Bs[i,j])*tilt_var_arr[i,j]
        tilt_term[i,j] = np.abs(1.0/Bz(I_l,r,z,r_f,0))*tilt_var_arr[i,j]
        if slide_term[i,j]<0 and tilt_term[i,j]>0:
            selector[i,j] = 1.0
        else:
            selector[i,j] = np.nan
alpha=1.0
beta=1.0
gamma=1.0
vert_term_norm = vert_term/max(abs(vert_term.flatten()))
slide_term_norm = slide_term/max(abs(slide_term.flatten()))
tilt_term_norm = tilt_term/max(abs(tilt_term.flatten()))
tot= (alpha * vert_term_norm + beta * slide_term_norm + gamma * tilt_term_norm)*selector

fig2=plt.figure()
axx1=plt.subplot(1,1,1)
CS=plt.contour(Y,X,tot, levels=np.linspace(np.nanmin(tot),np.nanmax(tot),300))
axx1.set_title('Experimental parameter space')
plt.contour(Y,X, tot, levels=[0.0], colors=['black'], linewidths=[5.0] )
axx1.set_ylabel(r'$z_L$ [$r_F$ units]') 
axx1.set_xlabel(r'$r_L$ [$r_F$ units]')


# normalized stability domains
fig=plt.figure(figsize=(8,12))
ax1=plt.subplot(3,1,1)
plt.contour(Y,X,vert_term_norm, levels=np.linspace(-0.75, 0.75, 500))
#plt.colorbar()
ax1.set_title('Vertical stability')
plt.contour(Y,X,vert_term_norm, levels=[0.0], colors=['black'], linewidths=[4.0])
ax1.set_ylabel(r'$z_L$ [$r_F$ units]') 

ax2=plt.subplot(3,1,2)
CS=plt.contour(Y,X,slide_term_norm, levels=np.linspace(-0.01, 0.01, 500))
#plt.colorbar()
ax2.set_title('Slide stability')
plt.contour(Y,X,slide_term_norm, levels=[0.0], colors=['black'], linewidths=[4.0])
ax2.set_ylabel(r'$z_L$ [$r_F$ units]') 

ax3=plt.subplot(3,1,3)
CS=plt.contour(Y,X,tilt_term_norm,levels=np.linspace(-0.05, 0.05, 500))
#plt.colorbar()
ax3.set_title('Tilt stability')
plt.contour(Y,X,tilt_term_norm, levels=[0.0], colors=['black'], linewidths=[4.0])
ax3.set_ylabel(r'$z_L$ [$r_F$ units]') 
ax3.set_xlabel(r'$r_L$ [$r_F$ units]') 


######
# Most recent parameter space visualization
fig=plt.figure(figsize=(8,12))
ax1=plt.subplot(3,1,1)
X,Y=np.meshgrid(r_l_arr,z_l_arr)
CS=plt.contourf(Y,X,np.log10(-dBrdz_arr),100)
plt.colorbar(label='log-amplitude')
ax1.set_title('Vertically unstable')
ax1.set_ylabel(r'$z_L$ [$r_F$ units]')

ax2=plt.subplot(3,1,2)
CS=plt.contourf(Y,X,np.log10(dBzdr_arr),100)
plt.colorbar(label='log-amplitude')
ax2.set_title('Slide unstable')
ax2.set_ylabel(r'$z_L$ [$r_F$ units]')

ax3=plt.subplot(3,1,3)
CS=plt.contourf(Y,X,np.log10(-tilt_var_arr),100)
cbar=plt.colorbar(label='log-amplitude')
ax3.set_title('Tilt unstable')
ax3.set_ylabel(r'$z_L$ [$r_F$ units]') 
ax3.set_xlabel(r'$r_L$ [$r_F$ units]') 
