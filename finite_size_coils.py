import numpy as np
import matplotlib.pyplot as plt
from scipy.special import ellipk, ellipe
import matplotlib.ticker as ticker
from coil_functions import *
from flux_surf_volume import get_flux_vol
import scipy.interpolate as interp
import matplotlib.patches as patches
plt.ion()
plt.style.use('ggplot') #'dark_background'
import sys
try:
    import nlopt
except:
    raise ValueError('nlopt not available. Use Python virtualenv!')

g = 9.80665
charging_pos=-0.15

# Reasonable test values (also initial guesses for optimization)
I_l=5.0e3
I_f=3895.32  #value from NIFS test
r_l=0.08
z_l=0.10  
r_f=0.0445 # Major radius 52mm - 15mm/2.0  #fixed/measured
m_f=0.2845 #fixed/measured
tz_f = 25.0e-3; tr_f = 15.0e-3 #fixed/measured
tz_l=20.0e-3; tr_l=10.0e-3
z_f=0.0   #fixed

# C coil parameters
I_c=10.0e3
r_c=0.08         # attempt larger than F-coil
z_c = z_f #-0.05 # guess for optimization should go to above
tz_c=25.0e-3; tr_c=30.0e-3  # chosen parameters, only for visualization (optimization approximates as line of charge)

# vector for normalization of parameters
x_init = [I_l,I_f,r_l,z_l]
                                               
# Optimization setup
def force_imbalance(x):    
    force_mag = - 2*np.pi*r_f*x[1]*Br(x[0],x[2],x[3],r_f,0.0)
    force_grav = m_f * g
    return abs(force_mag - force_grav)

def gauss(xx,xbar,s):
    x0=xx[0]/r_f; x1=xx[1]/r_f
    return np.exp(-((x0-xbar[0])**2+(x1-xbar[1])**2)/s**2)

xbar_slide=[2.5,3.0]; s_slide=1.5
xbar_tilt=[3.0,3.0]; s_tilt=1.5
opt1_cnt=0
tot_list=[];term1_list=[];term2_list=[];term3_list=[];term4_list=[];term5_list=[]

def tot_objective(x, grad, eta, x_init):
    # objective function for optimization, with regularization by parameters in eta
    # Normalization required for effective optimization
    I_l = x_init[0]*x[0]
    I_f = x_init[1]*x[1]
    r_l = x_init[2]*x[2]
    z_l = x_init[3]*x[3]
    term1=force_imbalance([I_l,I_f,r_l,z_l])
    term2 = -eta[0]*gauss([r_l,z_l],xbar_slide,s_slide)
    term3 = -eta[1]*gauss([r_l,z_l],xbar_tilt,s_tilt)
    term4= eta[2]*I_l
    term5= eta[3]*I_f

    global opt1_cnt; opt1_cnt+=1    
    tot=term1+term2+term3+term4+term5
    
    tot_list.append(tot); term1_list.append(term1); term2_list.append(term2)
    term3_list.append(term3); term4_list.append(term4); term5_list.append(term5)
    print tot, "---", term1,term2,term3,term4,term5
    return tot


def vertical_freq(x):
    return 1.0/(2*np.pi)*np.sqrt(-2*np.pi*r_f*x[1]*dBrdz(x[0],x[2],z_f,r_f,x[3])/m_f)

  
print "Optimizing L-coil parameters"
opt = nlopt.opt(nlopt.LN_SBPLX, 4)  # LN_SBPLX
opt.set_lower_bounds([0.1,0.99,0.3,0.2])# * opt.get_dimension())
opt.set_upper_bounds([3.0,1.1,2.0,2.0])# * opt.get_dimension())
opt.set_xtol_abs(1e-5)
 
eta_1 = [1.0e-5,5.0e-4,5.0e-5,1.0e-7]
objective= lambda x,grad: tot_objective(x,grad,eta_1,x_init)

opt.set_min_objective(objective)
guess = np.asarray([1.0,1.0,1.0,1.0])
res = opt.optimize(guess)

# Save optimized parameters:
I_l_opt = res[0]*x_init[0]
I_f_opt = res[1]*x_init[1]
r_l_opt = res[2]*x_init[2]
z_l_opt = res[3]*x_init[3]

print "*******************************************"
print "(real units) I_l=",I_l_opt, ", I_f=",I_f_opt, ", r_l=",r_l_opt, ", z_l=",z_l_opt
print "(F units) r_l=",r_l_opt/r_f, ", z_l=",z_l_opt/r_f
print "Force imbalance: ", force_imbalance([I_l_opt,I_f_opt,r_l_opt,z_l_opt])
print "Vertical instability frequency: ",  vertical_freq([I_l_opt,I_f_opt,r_l_opt,z_l_opt])
print "********************************************"

'''
n_points=6
tz_l_list=np.linspace(0.001,0.05,n_points)
Bs_Fcoil=[]
fields_plot=True

for i in range(len(tz_l_list)):
    tz_l=tz_l_list[i]
    
    # multicoil_fields requires central spatial positions
    X,Z,Bxm,Bzm,Bs,rAm=multicoil_fields([I_l_opt,I_f_opt],[z_l_opt,z_f],[r_l_opt,r_f],[tz_l,tz_f],[tr_l,tr_f])
    phi= 2*np.pi*rAm;
    if fields_plot:
        # Plot L and F coil flux surfaces. These are now taken to have finite size
        fig1 = plt.figure(1, figsize=(14,8))
        axx=fig1.add_subplot(n_points/2,2,i+1, axisbg='b')
    
        cs = plt.contourf(Z,X,phi*1e3,100,cmap=plt.cm.jet, extend="both")
        cbar = plt.colorbar()
        cbar.set_label(r' $2 \pi r A_{\theta} \times 10^3 [T m^2]$',fontsize=16)
        CS = plt.contour(Z,X,phi*1.0e3)
        plt.clabel(CS, inline=1, fontsize=13, linewidth=3, colors='white')
        plt.xlabel('r [m]'); plt.ylabel('z [m]')
    
        p1=patches.Rectangle((-r_f-tr_f/2.0,z_f-tz_f/2.0),tr_f,tz_f,edgecolor='red',facecolor='white')
        p2=patches.Rectangle((r_f-tr_f/2.0,z_f-tz_f/2.0),tr_f,tz_f,edgecolor='red',facecolor='white')
        axx.add_patch(p1); axx.add_patch(p2)

        p1=patches.Rectangle((-r_l_opt-tr_l/2.0,z_l_opt-tz_l/2.0),tr_l,tz_l,edgecolor='red',facecolor='white')
        p2=patches.Rectangle((r_l_opt-tr_l/2.0,z_l_opt-tz_l/2.0),tr_l,tz_l,edgecolor='red',facecolor='white')
        axx.add_patch(p1); axx.add_patch(p2)

        # Plot total magnetic field
        fig2 = plt.figure(2, figsize=(14,8))
        axx2=fig2.add_subplot(n_points/2,2,i+1, axisbg='b')
        cs = plt.contourf(Z,X,Bs,100,cmap=plt.cm.jet, extend="both")
        cbar = plt.colorbar()
        cbar.set_label(r'$B_{tot} [T]$', fontsize=16)
        CS = plt.contour(Z,X,Bs)
        plt.clabel(CS, inline=1, fontsize=13, linewidth=3, colors='white')
        plt.xlabel('r [m]'); plt.ylabel('z [m]')
    
        # draw shapes of coils
        p1=patches.Rectangle((-r_f-tr_f/2.0,z_f-tz_f/2.0),tr_f,tz_f,edgecolor='red',facecolor='white')
        p2=patches.Rectangle((r_f-tr_f/2.0,z_f-tz_f/2.0),tr_f,tz_f,edgecolor='red',facecolor='white')
        axx2.add_patch(p1); axx2.add_patch(p2)
    
        p1=patches.Rectangle((-r_l_opt-tr_l/2.0,z_l_opt-tz_l/2.0),tr_l,tz_l,edgecolor='red',facecolor='white')
        p2=patches.Rectangle((r_l_opt-tr_l/2.0,z_l_opt-tz_l/2.0),tr_l,tz_l,edgecolor='red',facecolor='white')
        axx2.add_patch(p1); axx2.add_patch(p2)

    # Select Bs at F-coil center:
    xidx=np.argmin(np.abs(X[0,:]-r_f))
    zidx=np.argmin(np.abs(Z[:,0]-z_f))
    Bs_Fcoil.append(Bs[xidx,zidx])


fig3=plt.figure()
axx3=fig3.add_subplot(111)
plt.plot(tz_l_list, Bs_Fcoil, 'ro')
plt.xlabel('Vertical width of L-coil  [m]')
plt.ylabel(r'$B_{tot}$ at F-coil center [T]')

'''


tz_l_list=[i*1.0e-3 for i in [36,25,50,18,75,12,30]]
tr_l_list=[i*1.0e-3 for i in [25,36,18,50,12,75,30]]

n_points=len(tz_l_list)
if n_points%2 is not 0:
    n_points+=1
    
Bs_Fcoil=[]
fields_plot=True

for i in range(len(tz_l_list)):
    tz_l=tz_l_list[i]
    tr_l=tr_l_list[i]
    
    # multicoil_fields requires central spatial positions
    X,Z,Bxm,Bzm,Bs,rAm=multicoil_fields([I_l_opt,I_f_opt],[z_l_opt,z_f],[r_l_opt,r_f],[tz_l,tz_f],[tr_l,tr_f])
    phi= 2*np.pi*rAm;
    if fields_plot:
        # Plot L and F coil flux surfaces. These are now taken to have finite size
        fig1 = plt.figure(1, figsize=(14,8))
        axx=fig1.add_subplot(n_points/2,2,i+1, axisbg='b')
    
        cs = plt.contourf(Z,X,phi*1e3,100,cmap=plt.cm.jet, extend="both")
        cbar = plt.colorbar()
        cbar.set_label(r' $2 \pi r A_{\theta} \times 10^3 [T m^2]$',fontsize=16)
        CS = plt.contour(Z,X,phi*1.0e3)
        plt.clabel(CS, inline=1, fontsize=13, linewidth=3, colors='white')
        plt.xlabel('r [m]'); plt.ylabel('z [m]')
        plt.title('tz_l = %f, tr_l = %f'%(tz_l,tr_l))
        
        p1=patches.Rectangle((-r_f-tr_f/2.0,z_f-tz_f/2.0),tr_f,tz_f,edgecolor='red',facecolor='white')
        p2=patches.Rectangle((r_f-tr_f/2.0,z_f-tz_f/2.0),tr_f,tz_f,edgecolor='red',facecolor='white')
        axx.add_patch(p1); axx.add_patch(p2)

        p1=patches.Rectangle((-r_l_opt-tr_l/2.0,z_l_opt-tz_l/2.0),tr_l,tz_l,edgecolor='red',facecolor='white')
        p2=patches.Rectangle((r_l_opt-tr_l/2.0,z_l_opt-tz_l/2.0),tr_l,tz_l,edgecolor='red',facecolor='white')
        axx.add_patch(p1); axx.add_patch(p2)

        # Plot total magnetic field
        fig2 = plt.figure(2, figsize=(14,8))
        axx2=fig2.add_subplot(n_points/2,2,i+1, axisbg='b')
        cs = plt.contourf(Z,X,Bs,100,cmap=plt.cm.jet, extend="both")
        cbar = plt.colorbar()
        cbar.set_label(r'$B_{tot} [T]$', fontsize=16)
        CS = plt.contour(Z,X,Bs)
        plt.clabel(CS, inline=1, fontsize=13, linewidth=3, colors='white')
        plt.xlabel('r [m]'); plt.ylabel('z [m]')
        plt.title('tz_l = %f, tr_l = %f'%(tz_l,tr_l))
        
        # draw shapes of coils
        p1=patches.Rectangle((-r_f-tr_f/2.0,z_f-tz_f/2.0),tr_f,tz_f,edgecolor='red',facecolor='white')
        p2=patches.Rectangle((r_f-tr_f/2.0,z_f-tz_f/2.0),tr_f,tz_f,edgecolor='red',facecolor='white')
        axx2.add_patch(p1); axx2.add_patch(p2)
    
        p1=patches.Rectangle((-r_l_opt-tr_l/2.0,z_l_opt-tz_l/2.0),tr_l,tz_l,edgecolor='red',facecolor='white')
        p2=patches.Rectangle((r_l_opt-tr_l/2.0,z_l_opt-tz_l/2.0),tr_l,tz_l,edgecolor='red',facecolor='white')
        axx2.add_patch(p1); axx2.add_patch(p2)

    # Select Bs at F-coil center:
    xidx=np.argmin(np.abs(X[0,:]-r_f))
    zidx=np.argmin(np.abs(Z[:,0]-z_f))
    Bs_Fcoil.append(Bs[xidx,zidx])


fig3=plt.figure()
axx3=fig3.add_subplot(111)
plt.plot(tz_l_list, Bs_Fcoil, 'ro')
plt.xlabel('Vertical width of L-coil  [m]')
plt.ylabel(r'$B_{tot}$ at F-coil center [T]')


