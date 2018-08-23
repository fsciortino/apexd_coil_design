''' Script to assess levitation procedure that begins with both C- and L-coils exerting magnetic forces on the F-coil from below and above.

Note that once the F-coil reaches z~0 with small velocity and acceleration, the PID controller circuit may take control over and stabilize vertical motion. At that stage, the C-coil may be turned off and used as a safety mechanism in case of L-coil or PID circuit failure. 

F.Sciortino, 7/6/18
'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import ellipk, ellipe
import matplotlib.ticker as ticker
from coil_functions import *
from flux_surf_volume import get_flux_vol
import matplotlib as mpl
mpl.rcParams['xtick.labelsize']=12
mpl.rcParams['ytick.labelsize']=12

plt.ion()
plt.style.use('ggplot')

g = 9.80665
charging_pos=-0.15
starting_pos=-0.15
option=1

# Values from test optimization
I_l=3640.625
I_f=3900.79779375
r_l= 0.100108424677 # explore twice as large as F-coil
z_l= 0.100765730835 
r_f=0.0445 # Major radius 52mm - 15mm/2.0  #fixed/measured
m_f=0.2845 #fixed/measured
tz_f = 25.0e-3; tr_f = 15.0e-3 #fixed/measured
tz_l=20.0e-3; tr_l=10.0e-3  # chosen parameters
z_f=0.0   #fixed

# C coil parameters
I_c= - 369.05
r_c=0.065473130123         # attempt slightly larger than F-coil
z_c = -0.025 + charging_pos # position well below floating equilibrium location
tz_c=25.0e-3; tr_c=30.0e-3  # chosen parameters

# Test 2-coil lifting
force_L = coil_force(r_l,r_f,I_f,I_l,charging_pos,z_l)
force_C = coil_force(r_c,r_f,I_f,I_c,charging_pos,z_c)
force_g = - m_f * g

print force_L, force_C, force_g
force_imbalance = force_L + force_C + force_g 

mu0=4.0*np.pi*10**(-7)

########
def dBrdz_comb(I_l,I_c,r_l,r_c,z_l,z_c,r_f,z_f):
    delta=1e-6
    B1=Br(I_l,r_l,z_l,r_f,z_f+delta)+Br(I_c,r_c,z_c,r_f,z_f+delta)
    B2=Br(I_l,r_l,z_l,r_f,z_f)+Br(I_c,r_c,z_c,r_f,z_f)
    return (B1-B2)/delta

def dBzdr_comb(I_l,I_c,r_l,r_c,z_l,z_c,r_f,z_f):
    delta=1e-6
    B1=Bz(I_l,r_l,z_l,r_f+delta,z_f)+Bz(I_c,r_c,z_c,r_f+delta,z_f)
    B2=Bz(I_l,r_l,z_l,r_f,z_f)+Bz(I_c,r_c,z_c,r_f,z_f)
    return (B1-B2)/delta

########
# C and L coils working together
if option==0:
    v_init=0.0
    z_init= charging_pos
    force_L = coil_force(r_l,r_f,I_f,I_l,z_init,z_l)
    force_C = coil_force(r_c,r_f,I_f,I_c,z_init,z_c)
    force_g = - m_f * g
    force_imbalance = force_L + force_C + force_g
    a_init=force_imbalance/g

    n_steps=5000
    p=[] #np.zeros((n_steps,3))
    force_L_list=[]; force_C_list=[]; force_g_list=[]
    force_L_list.append(force_L); force_C_list.append(force_C); force_g_list.append(force_g)
    p.append((z_init,v_init,a_init))
    dt=1e-3
    time=np.arange(0,n_steps*dt,dt)
    dBrdz_arr=np.zeros(n_steps)
    dBzdr_arr=np.zeros(n_steps)
    tilt_var_arr=np.zeros(n_steps)
    dBrdz_arr[0]=dBrdz_comb(I_l,I_c,r_l,r_c,z_l,z_c,r_f,p[-1][0])
    dBzdr_arr[0]=dBzdr_comb(I_l,I_c,r_l,r_c,z_l,z_c,r_f,p[-1][0])
    tilt_var_arr[0]=Bz(I_l,r_l,z_l,r_f,p[-1][0])+Bz(I_c,r_c,z_c,r_f,p[-1][0])+r_f * dBrdz_arr[0]

    for i in range(1,n_steps):
        # fix currents and only follow motion
        force_L = coil_force(r_l,r_f,I_f,I_l,p[-1][0],z_l)
        force_C = coil_force(r_c,r_f,I_f,I_c,p[-1][0],z_c)
        force_g = - m_f * g
        force_L_list.append(force_L)
        force_C_list.append(force_C)
        force_g_list.append(force_g)
        force_imbalance = force_L + force_C + force_g
        a_new=force_imbalance/g
        v_new=p[-1][1]+dt*p[-1][2]
        z_new=p[-1][0]+dt*p[-1][1]

        p.append((z_new,v_new,a_new))

        dBrdz_arr[i]=dBrdz_comb(I_l,I_c,r_l,r_c,z_l,z_c,r_f,p[-1][0])
        dBzdr_arr[i]=dBzdr_comb(I_l,I_c,r_l,r_c,z_l,z_c,r_f,p[-1][0])

        tilt_var_arr[i]=Bz(I_l,r_l,z_l,r_f,p[-1][0])+Bz(I_c,r_c,z_c,r_f,p[-1][0])+r_f * dBrdz_arr[i]

    p=np.asarray(p)
    if p.shape[0]!=3: p=np.asarray(p).T

    plt.figure(figsize=(14,10))
    plt.suptitle(r'$I_F=%d A, I_L= %d A, I_C=%d A, r_F=%f m, r_L=%f m, r_C=%f m, z_{init}=%f m, z_L=%f m$'%(I_f,I_l,I_c,r_f,r_l,r_c,charging_pos,z_l), fontsize=15)
    ax1=plt.subplot(321)
    plt.plot(time,p[0,:],'ro')
    plt.ylabel(r'$z [m]$',fontsize=16)

    ax2=plt.subplot(323,sharex=ax1)
    plt.plot(time,p[1,:],'ro')
    plt.ylabel(r'$v [m/s]$',fontsize=16)

    ax3=plt.subplot(325,sharex=ax1)
    plt.plot(time,p[2,:],'ro')
    plt.ylabel(r'$a [m/s^2]$',fontsize=16)
    plt.xlabel(r't [s]',fontsize=16)

    axx1=plt.subplot(322)
    plt.plot(time,force_L_list,'ro')
    plt.ylabel(r'$F_L$',fontsize=16)

    axx2=plt.subplot(324,sharex=axx1)
    plt.plot(time,force_C_list,'ro')
    plt.ylabel(r'$F_C$',fontsize=16)

    axx3=plt.subplot(326,sharex=axx1)
    plt.plot(time,force_g_list,'ro')
    plt.ylabel(r'$F_g$',fontsize=16)
    plt.xlabel(r't [s]',fontsize=16)

    # stability
    plt.figure(figsize=(14,10))
    axxx1=plt.subplot(311)
    plt.plot(time,dBrdz_arr,'bo')
    plt.plot([min(time),max(time)],[0.0,0.0],'r--')
    plt.ylabel(r'$dB_r/dz [T/m]$',fontsize=16)

    axxx2=plt.subplot(312,sharex=axxx1)
    plt.plot(time,-dBzdr_arr,'bo')
    plt.plot([min(time),max(time)],[0.0,0.0],'r--')
    plt.ylabel(r'$-dB_z/dr [T/m]$',fontsize=16)

    axxx3=plt.subplot(313,sharex=ax1)
    plt.plot(time,tilt_var_arr,'bo')
    plt.plot([min(time),max(time)],[0.0,0.0],'r--')
    plt.ylabel(r'$B_{z0}+r_F \ \frac{\partial B_r}{\partial z}|_{z=0} [T]$',fontsize=16)
    plt.xlabel(r't [s]',fontsize=16)






    
# All work done by L-coil
elif option==1:
    v_init=0.0
    z_init= charging_pos
    n_steps=5000
    dt=1e-3
    time=np.arange(0,n_steps*dt,dt)
    
    Ifinal=I_l
    fact=0.9
    def Ilift(x,Ifinal,z_init):
        # modify to fix aimed position of z_F 3cm below z=0
        z_init=z_init+0.03
        global fac
        val = fact*Ifinal*abs(((z_init+x)/z_init))**3
        return val
        
    xlist=np.linspace(z_init,0.10,1000)
    Il_arr_tmp=np.asarray([Ilift(x,Ifinal,z_init) for x in xlist])
    plt.figure(); plt.plot(xlist,Il_arr_tmp)
    plt.xlabel(r'$\Delta z$ $[m]$'); plt.ylabel(r'$I_L$ $[AT]$')
    
    v_init=0.0
    z_init= charging_pos
    Iinit=Ilift(z_init,Ifinal,z_init)
    force_L = coil_force(r_l,r_f,I_f,Iinit,z_init,z_l)
    force_g = - m_f * g
    force_imbalance = force_L +  force_g
    a_init=force_imbalance/g

    p=[]
    p.append((z_init,v_init,a_init))
      
    force_L_list=[]; force_g_list=[]
    force_L_list.append(force_L); force_g_list.append(force_g)

    Il_arr=np.zeros(n_steps)
    dBrdz_arr=np.zeros(n_steps)
    dBzdr_arr=np.zeros(n_steps)
    tilt_var_arr=np.zeros(n_steps)
    dBrdz_arr[0]=dBrdz(Iinit,r_l,z_l,r_f,p[-1][0])
    dBzdr_arr[0]=dBzdr(Iinit,r_l,z_l,r_f,p[-1][0])
    tilt_var_arr[0]=Bz(Iinit,r_l,z_l,r_f,p[-1][0])+r_f * dBrdz_arr[0]

    for i in range(1,n_steps):
        # fix currents and only follow motion
        Il_new=Ilift(p[-1][0],Ifinal,z_init)
        Il_arr[i]=Il_new
        force_L = coil_force(r_l,r_f,I_f,Il_new,p[-1][0],z_l)
        force_g = - m_f * g
        force_L_list.append(force_L)
        force_g_list.append(force_g)
        force_imbalance = force_L + force_g
        a_new=force_imbalance/g

        # Euler method
        v_new=p[-1][1]+dt*p[-1][2]
        z_new=p[-1][0]+dt*p[-1][1]

        p.append((z_new,v_new,a_new))

        dBrdz_arr[i]=dBrdz(Il_new,r_l,z_l,r_f,p[-1][0])
        dBzdr_arr[i]=dBzdr(Il_new,r_l,z_l,r_f,p[-1][0])
        tilt_var_arr[i]=Bz(Il_new,r_l,z_l,r_f,p[-1][0])+r_f*dBrdz_arr[i]

    p=np.asarray(p)
    if p.shape[0]!=3: p=np.asarray(p).T

    plt.figure(figsize=(14,10))
    plt.suptitle(r'$I_F=%d A, I_L= %d A, I_C=%d A, r_F=%f m, r_L=%f m, r_C=%f m, z_{init}=%f m, z_L=%f m$'%(I_f,I_l,I_c,r_f,r_l,r_c,charging_pos,z_l), fontsize=15)
    ax1=plt.subplot(321)
    plt.plot(time,p[0,:],'ro')
    plt.ylabel(r'$z [m]$',fontsize=16)

    ax2=plt.subplot(323,sharex=ax1)
    plt.plot(time,p[1,:],'ro')
    plt.ylabel(r'$v [m/s]$',fontsize=16)

    ax3=plt.subplot(325,sharex=ax1)
    plt.plot(time,p[2,:],'ro')
    plt.ylabel(r'$a [m/s^2]$',fontsize=16)
    plt.xlabel(r't [s]',fontsize=16)

    axx1=plt.subplot(322)
    plt.plot(time,force_L_list,'ro')
    plt.ylabel(r'$F_L$',fontsize=16)

    axx2=plt.subplot(324,sharex=axx1)
    plt.plot(time,force_g_list,'ro')
    plt.ylabel(r'$F_g$',fontsize=16)
    
    axx3=plt.subplot(326,sharex=axx1)
    plt.plot(time,Il_arr,'ro')
    plt.ylabel(r'$I_L$',fontsize=16)
    plt.xlabel(r't [s]',fontsize=16)

    

    # stability
    plt.figure(figsize=(14,10))
    axxx1=plt.subplot(311)
    plt.plot(time,dBrdz_arr,'bo')
    plt.plot([min(time),max(time)],[0.0,0.0],'r--')
    plt.ylabel(r'$dB_r/dz [T/m]$',fontsize=16)

    axxx2=plt.subplot(312,sharex=axxx1)
    plt.plot(time,-dBzdr_arr,'bo')
    plt.plot([min(time),max(time)],[0.0,0.0],'r--')
    plt.ylabel(r'$-dB_z/dr [T/m]$',fontsize=16)

    axxx3=plt.subplot(313,sharex=ax1)
    plt.plot(time,tilt_var_arr,'bo')
    plt.plot([min(time),max(time)],[0.0,0.0],'r--')
    plt.ylabel(r'$B_{z0}+r_F \ \frac{\partial B_r}{\partial z}|_{z=0} [T]$',fontsize=16)
    plt.xlabel(r't [s]',fontsize=16)


