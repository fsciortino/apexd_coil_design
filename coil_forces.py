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

try: fix_l=bool(sys.argv([1]))
except: fix_l=False
try: fix_c=bool(sys.argv([2]))
except: fix_c=True
try: plot_stab=bool(sys.argv([3]))
except: plot_stab=True
try: plot_bias=bool(sys.argv([4]))
except: plot_bias=False
try: plot_opt1=bool(sys.argv([5]))
except: plot_opt1=True
try: plot_Bs=bool(sys.argv([6]))
except: plot_Bs=True
try: get_contour_vol=bool(sys.argv([7]))
except: get_contour_vol=False
try: get_persistent_If = bool(sys.argv([8]))
except: get_persistent_If=False
try: get_persistent_If_2 = bool(sys.argv([9]))
except: get_persistent_If_2 = True
try: fv_calc = bool(sys.argv([10]))
except: fv_calc = False

g = 9.80665
charging_pos=-0.15
Ips=20 # current from power supply
d_wire = 2.0e-3 # wire diameter

# L-coil
nr_l=9  ; nz_l=47 # fixed during construction. Dependent on bobbin size
num_turns_l = nr_l * nz_l
tz_l=100.0e-3; tr_l=22.0e-3
I_l=nr_l*nz_l*Ips #5.0e3
r_l=0.085 #0.106 #0.08# expected L-coil winding # 0.08 #former test
z_l=0.083 #0.07 #10

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

# vector for normalization of parameters
x_init = [I_l,I_f,r_l,z_l]

#### If appropriate, include flux conservation for coil force balance calculation
flux_cons=True
if flux_cons:
    # Compute using flux freezing:
    X,Z,Bxm,Bzm,Bs,rAm = multicoil_fields([I_l],[z_l],[r_l],[tz_l],[tr_l],[nz_l],[nr_l])
    xidx=np.argmin(np.abs(X[0,:]-z_f))
    zidx=np.argmin(np.abs(Z[:,0]-r_f+tr_f/2.0))
    rAm_0=rAm[:,xidx]
    
    X_f,Z_f,Bxm_f,Bzm_f,Bs_f,rAm_f = multicoil_fields([I_f],[z_f],[r_f],[tz_f],[tr_f],[nz_f],[nr_f])
    xidx_f=np.argmin(np.abs(X_f[0,:]-z_f))
    zidx_f=np.argmin(np.abs(Z_f[:,0]-r_f+tr_f/2.0))
    rAm_f=rAm_f[:,xidx_f]
    
    frac= rAm_0[zidx] / rAm_f[zidx_f]
    I_f= frac* I_f
###
# ------------- ------------------
force1 = -I_f*I_l *dMdz(r_f,r_l,z_f,z_l) #inverted 
freq = 1.0/(2*np.pi) * np.sqrt(-2*np.pi*r_f *I_f* dBrdz(I_l,r_l,z_f,r_f,z_l)/m_f) # unclear why -ve sign needed

#force2 = -2*np.pi*r_f*I_f*Br(I_l,r_l,z_l,r_f,z_f) 
force2 = coil_force(r_l,r_f,I_f,I_l,z_f,z_l)
force2_fv = coil_force_fv(r_l,r_f,I_f,I_l,z_f,z_l,nr_l,nz_l)

force_imb=force2 - m_f * g
print "Initial vertical resonance frequency: ", freq
print "Initial forces (mag vs grav): ", force2, m_f*g
print "Initial force imbalance: ", force_imb

'''
# check if increase of current by 1 A (consider 100 turns) can counter balance this force
force3=coil_force(r_l,r_f,I_f,I_l+100,z_f,z_l)
force4=coil_force(r_l,r_f,I_f,I_l+100,z_f,z_l) 
print "Feasible force with +100A:", force3
if np.sign(force3-m_f*g)!=np.sign(force2-m_f*g) or np.sign(force4-m_f*g)!=np.sign(force2-m_f*g):
    print "Reversible imbalance!"
'''

# ------------------------------
# Stability plots (set r_f=1.0)

# fix (r,z) for SC F-coil and vary L-coil parameters
if plot_stab:
    r_l_arr=np.arange(0.05, 5.0, 0.03)
    z_l_arr=np.arange(0.05, 5.0, 0.03)
    dBrdz_arr=np.zeros((r_l_arr.shape[0],z_l_arr.shape[0]))
    Br_arr=np.zeros((r_l_arr.shape[0],z_l_arr.shape[0]))
    dBzdr_arr=np.zeros((r_l_arr.shape[0],z_l_arr.shape[0]))
    tilt_var_arr=np.zeros((r_l_arr.shape[0],z_l_arr.shape[0]))
    unstable_param=-np.ones((r_l_arr.shape[0],z_l_arr.shape[0]))
    for i in range(r_l_arr.shape[0]):
        for j in range(z_l_arr.shape[0]):
            r=r_l_arr[i]; z=z_l_arr[j]
            if fv_calc:
                Br_arr[i,j]=Br_fv(I_l,r,z,1.0,0,nr_l,nz_l)
                dBrdz_arr[i,j]=dBrdz_fv(I_l,r,z,1.0,0,nr_l,nz_l)
                dBzdr_arr[i,j]=dBzdr_fv(I_l,r,z,1.0,0,nr_l,nz_l)
                tilt_var_arr[i,j]=Bz_fv(I_l,r,z,1.0,0,nr_l,nz_l)+1.0 * dBrdz_arr[i,j]
            else:
                Br_arr[i,j]=Br(I_l,r,z,1.0,0)
                dBrdz_arr[i,j]=dBrdz(I_l,r,z,1.0,0) # two inversions cancel out
                dBzdr_arr[i,j]=dBzdr(I_l,r,z,1.0,0)
                tilt_var_arr[i,j]=Bz(I_l,r,z,1.0,0)+1.0 * dBrdz_arr[i,j]
            if dBzdr_arr[i,j]>0:
                unstable_param[i,j]=dBzdr_arr[i,j]
            if -tilt_var_arr[i,j]>0: #note -ve sign
                unstable_param[i,j]+= tilt_var_arr[i,j]
            # note that magnitudes of unstable_param are arbitrary
            
    # plot 
    fig=plt.figure(figsize=(8,10))
    ax1=plt.subplot(3,1,1)
    X,Y=np.meshgrid(r_l_arr,z_l_arr)
    CS=plt.contourf(Y,X,np.log10(-dBrdz_arr),100)
    plt.colorbar(label='log-amplitude')
    ax1.set_title('Vertically unstable')
    ax1.set_ylabel(r'z [$r_F$]')
    
    ax2=plt.subplot(3,1,2)
    CS=plt.contourf(Y,X,np.log10(dBzdr_arr),100)
    plt.colorbar(label='log-amplitude')
    ax2.set_title('Slide unstable')
    ax2.set_ylabel(r'z $[r_F]$')
    
    ax3=plt.subplot(3,1,3)
    CS=plt.contourf(Y,X,np.log10(-tilt_var_arr),100)
    cbar=plt.colorbar(label='log-amplitude')
    ax3.set_title('Tilt unstable')
    ax3.set_ylabel(r'z $[r_F]$')
    ax3.set_xlabel(r'r $[r_F]$')
    
    plt.figure()
    ax11=plt.subplot(1,1,1)
    plt.contourf(Y,X,np.log10(unstable_param),100)
    plt.colorbar(label='log-amplitude')
    ax11.set_title('Vertically unstable')

# -------------------------------
def force_imbalance_fv(x):
    I_l = x[0]
    I_f = x[1]
    r_l = x[2]
    z_l = x[3]
    nr_l= x[4] # use  number of turns rather than thickness, assuming wire with d=2mm
    nz_l=x[5]

    # Current through each turn of L-coil:
    num_turns=nz_l*nr_l
    I_l_per_turn = float(I_l) /num_turns # units: A (not AT)

    # get magnetic force from each turn
    fm = 0.0
    for iz in range(nz_l):
        # Vertical position of turn (recall that z_l and r_l are positions of coil center)
        z_turn = z_l - d_wire * float(nz_l)/2.0 + iz * d_wire
        
        for ir in range(nr_l):
            # radial position of turn
            r_turn = r_l - d_wire * float(nr_l)/2.0 + ir * d_wire

            # Compute force from turn in L-coil on F-coil
            fm += - 2*np.pi*r_f*I_f*Br(I_l_per_turn,r_turn,z_turn,r_f,z_f)

    # subtract gravitational force
    fg = m_f * g

    # compare to computation that does not resolve finite volume of L-coil:
    ###force_mag = - 2*np.pi*r_f*I_f*Br(I_l,r_l,z_l,r_f,z_f)

    ###print "With finite volume effects: ", fm
    ###print "Without: ", force_mag
    
    return abs(fm - fg)
    
    
# Optimization setup
def force_imbalance(x):    
    force_mag = - 2*np.pi*r_f*x[1]*Br(x[0],x[2],x[3],r_f,z_f)
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
    # 'grad' is needed by nlopt, but can be left empty in our case (optimization will not use derivatives of the objective function)
    I_l = x_init[0]*x[0]
    I_f = x_init[1]*x[1]
    r_l = x_init[2]*x[2]
    z_l = x_init[3]*x[3]

    ###########
    if flux_cons:
        # Obtain I_f correction, in consideration of flux freezing:
        X,Z,Bxm,Bzm,Bs,rAm = multicoil_fields([I_l],[z_l],[r_l],[tz_l],[tr_l],[nz_l],[nr_l])
        xidx=np.argmin(np.abs(X[0,:]-z_f))
        zidx=np.argmin(np.abs(Z[:,0]-r_f+tr_f/2.0))
        rAm_0=rAm[:,xidx]

        X_f,Z_f,Bxm_f,Bzm_f,Bs_f,rAm_f = multicoil_fields([I_f],[z_f],[r_f],[tz_f],[tr_f],[nz_f],[nr_f])
        xidx_f=np.argmin(np.abs(X_f[0,:]-z_f))
        zidx_f=np.argmin(np.abs(Z_f[:,0]-r_f+tr_f/2.0))
        rAm_f=rAm_f[:,xidx_f]

        frac= rAm_0[zidx] / rAm_f[zidx_f]
        I_f= frac* I_f  # corrected I_f <-- tends to be much smaller than I_f if flux from L-coil is large

    ############
    term1=force_imbalance_fv([I_l,I_f,r_l,z_l,nr_l,nz_l])
    #term1=force_imbalance_fv([I_l,I_f,r_l,z_l,nr_l,nz_l]) # number of turns not currently optimized
    term2 = -eta[0]*gauss([r_l,z_l],xbar_slide,s_slide)
    term3 = -eta[1]*gauss([r_l,z_l],xbar_tilt,s_tilt)
    term4= eta[2]*I_l
    term5= eta[3]*I_f

    global opt1_cnt; opt1_cnt+=1    
    tot=term1+term2+term3+term4+term5

    tot_list.append(tot); term1_list.append(term1); term2_list.append(term2)
    term3_list.append(term3); term4_list.append(term4); term5_list.append(term5)
    #print tot, "---", term1,term2,term3,term4,term5
    return tot


if plot_bias:
    plt.figure()
    x1=np.linspace(0.0,5.0,100); x2=np.linspace(0.0,5.0,100);
    stab_bias_slide=np.zeros((len(x1),len(x2)))
    stab_bias_tilt=np.zeros((len(x1),len(x2)))
    for i in range(len(x1)):
        for j in range(len(x2)):
            stab_bias_slide[i,j]= np.exp(-((x1[i]-xbar_slide[0])**2+(x2[j]-xbar_slide[1])**2)/s_slide**2)
            stab_bias_tilt[i,j]=np.exp(-((x1[i]-xbar_tilt[0])**2+(x2[j]-xbar_tilt[1])**2)/s_tilt**2)
    X1,X2=np.meshgrid(x1,x2)
    axx1=plt.subplot(2,1,1)
    plt.contourf(X2,X1,stab_bias_slide,100)
    cbar=plt.colorbar(label='slide bias amplitude')
    plt.ylabel(r'z [$r_f$ units]')
    axx2=plt.subplot(2,1,2)
    plt.contourf(X2,X1,stab_bias_tilt,100)
    cbar=plt.colorbar(label='tilt bias amplitude')
    plt.ylabel(r'z [r$_f$ units]')
    plt.xlabel(r'r [r$_f$ units]')
  
def vertical_freq(x):
    return 1.0/(2*np.pi)*np.sqrt(-2*np.pi*r_f*x[1]*dBrdz(x[0],x[2],z_f,r_f,x[3])/m_f)

  
if not fix_l:    
    print "Optimizing L-coil parameters"
    opt = nlopt.opt(nlopt.LN_SBPLX, 4)  
    #opt.set_lower_bounds([0.3,0.3,0.7,0.7])# * opt.get_dimension())
    #opt.set_upper_bounds([2.0,2.0,2.0,2.0])# * opt.get_dimension())
    opt.set_lower_bounds([0.1,0.99,0.99,0.2])# * opt.get_dimension())
    opt.set_upper_bounds([3.0,1.01,1.01,2.0])# * opt.get_dimension())
    opt.set_xtol_abs(1e-5)
    
    eta_1 = [1.0e-5,5.0e-4,1.0e-6,1.0e-7]
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
   # print "(initial units) I_l=",res[0], ", I_f=",res[1], ", r_l=",res[2], ", z_l=",res[3]
    print "(real units) I_l=",I_l_opt, ", I_f=",I_f_opt, ", r_l=",r_l_opt, ", z_l=",z_l_opt
    print "(F units) r_l=",r_l_opt/r_f, ", z_l=",z_l_opt/r_f
    print "Force imbalance: ", force_imbalance_fv([I_l_opt,I_f_opt,r_l_opt,z_l_opt,nr_l,nz_l])
    print "Vertical instability frequency: ",  vertical_freq([I_l_opt,I_f_opt,r_l_opt,z_l_opt])
    print "********************************************"
    
    # Plot L and F coil flux surfaces. These are now taken to have finite size
    fig1 = plt.figure()
    axx=fig1.add_subplot(1,1,1, axisbg='b')
    
    # multicoil_fields requires central spatial positions
    X,Z,Bxm,Bzm,Bs,rAm=multicoil_fields([I_l_opt,I_f_opt],[z_l_opt,z_f],[r_l_opt,r_f],[tz_l,tz_f],[tr_l,tr_f],[nz_l,nz_f],[nr_l,nr_f])
   
    cs = plt.contourf(Z,X,rAm*1e3,100,cmap=plt.cm.jet, extend="both")
    cbar = plt.colorbar()
    cbar.set_label(r' $2 \pi r A_{\theta} \times 10^3 [T m^2]$',fontsize=16)
    CS = plt.contour(Z,X,rAm*1.0e3)
    plt.clabel(CS, inline=1, fontsize=13, linewidth=3, colors='white')
    plt.xlabel('r [m]'); plt.ylabel('z [m]')
    
    p1=patches.Rectangle((-r_f-tr_f/2.0,z_f-tz_f/2.0),tr_f,tz_f,edgecolor='red',facecolor='white')
    p2=patches.Rectangle((r_f-tr_f/2.0,z_f-tz_f/2.0),tr_f,tz_f,edgecolor='red',facecolor='white')
    axx.add_patch(p1); axx.add_patch(p2)
    
    p1=patches.Rectangle((-r_l_opt-tr_l/2.0,z_l_opt-tz_l/2.0),tr_l,tz_l,edgecolor='red',facecolor='white')
    p2=patches.Rectangle((r_l_opt-tr_l/2.0,z_l_opt-tz_l/2.0),tr_l,tz_l,edgecolor='red',facecolor='white')
    axx.add_patch(p1); axx.add_patch(p2)

    # Plot total magnetic field
    if plot_Bs:
        fig1 = plt.figure()
        axx=fig1.add_subplot(1,1,1, axisbg='b')
        cs = plt.contourf(Z,X,Bs,100,cmap=plt.cm.jet, extend="both")
        cbar = plt.colorbar()
        cbar.set_label(r'$B_{tot} [T]$', fontsize=16)
        CS = plt.contour(Z,X,Bs)
        plt.clabel(CS, inline=1, fontsize=13, linewidth=3, colors='white')
        plt.xlabel('r [m]'); plt.ylabel('z [m]')
    
        # draw shapes of coils
        p1=patches.Rectangle((-r_f-tr_f/2.0,z_f-tz_f/2.0),tr_f,tz_f,edgecolor='red',facecolor='white')
        p2=patches.Rectangle((r_f-tr_f/2.0,z_f-tz_f/2.0),tr_f,tz_f,edgecolor='red',facecolor='white')
        axx.add_patch(p1); axx.add_patch(p2)
    
        p1=patches.Rectangle((-r_l_opt-tr_l/2.0,z_l_opt-tz_l/2.0),tr_l,tz_l,edgecolor='red',facecolor='white')
        p2=patches.Rectangle((r_l_opt-tr_l/2.0,z_l_opt-tz_l/2.0),tr_l,tz_l,edgecolor='red',facecolor='white')
        axx.add_patch(p1); axx.add_patch(p2)
else:
    print "Using initial parameters for I_l, I_f, r_l, z_l"
    I_l_opt = I_l
    I_f_opt = I_f
    r_l_opt = r_l
    z_l_opt = z_l
    

########## C-coil Optimization  #############
                 
if not fix_c:
    print "Optimizing C-coil parameters"

    # Ic,rc,zc
    c_bounds=[(500.0,20000.0),(r_f+tr_f/2.0+tr_c/2.0,0.3),(-0.1,0.1)]
    c_bounds=[(I_c*0.99,I_c*1.01),(r_c*0.99, r_c*1.01),(1e-5*0.99,1e-5*1.01)]

    # set to 1 any parameters in [Ic,rc,zc] that should be kept fixed to the original value
    fixed_params=[0,1,1] 
    
    # keep all parameters +ve
    I_c_norm=I_c/c_bounds[0][1]; r_c_norm=r_c/c_bounds[1][1]; z_c_norm=abs(z_c)/c_bounds[2][1]
    c_params_guess = I_c_norm,r_c_norm,z_c_norm,tz_c,tr_c,nz_c,nr_c

    nu_c=[8.0e-9]
    f_params = I_f_opt,r_f,z_f,tz_f,tr_f,nz_f,nr_f
    I_c_tmp,r_c_tmp,z_c_tmp=c_coil_opt(f_params,c_params_guess,c_bounds,nu_c, fixed_params)
    I_c_opt=I_c_tmp*c_bounds[0][1]
    r_c_opt=r_c_tmp*c_bounds[1][1]
    z_c_opt=-z_c_tmp*c_bounds[2][1]

    # # # # # Check the flux imbalance: # # # # #
    X_f,Z_f,Bxm_f,Bzm_f,Bs_f,rAm_f = multicoil_fields([I_f_opt],[z_f],[r_f],[tz_f],[tr_f],[nz_f],[nr_f])
    xidx_f=np.argmin(np.abs(X_f[0,:]-r_f))
    zidx_f=np.argmin(np.abs(Z_f[:,0]-z_f))
    rAm_xz_2=rAm_f[xidx_f,zidx_f]
    
    X,Z,Bxm,Bzm,Bs,rAm = multicoil_fields([I_c_opt],[z_c_opt],[r_c_opt],[tz_c],[tr_c],[nz_c],[nr_c])
    xidx=np.argmin(np.abs(X[0,:]-r_f))
    zidx=np.argmin(np.abs(Z[:,0]-z_f))
    rAm_xz_1=rAm[xidx,zidx]
    
    flux_imbalance = np.abs(rAm_xz_1 - rAm_xz_2)
    # # # # # # # # # # # # # # # # # # # # # # #
    
    print "************************"
    print "Optimized coil parameters:"
    print "I_l=",I_l_opt, ", I_f=",I_f_opt, ", I_c=", I_c_opt
    print "r_l=",r_l_opt, ", r_f=",r_f, ", r_c=",r_c_opt
    print "z_l=", z_l_opt,", z_f=",z_f, ", z_c=",z_c_opt
    print "Resulting force imbalance: ", force_imbalance([I_l_opt,I_f_opt,r_l_opt,z_l_opt,nz_l,nr_l])
    print "Resulting flux imbalance: ", flux_imbalance
    print "Vertical instability frequency: ",  vertical_freq([I_l_opt,I_f_opt,r_l_opt,z_l_opt])
    print "************************"
    
    fig2 = plt.figure()
    axxx=fig2.add_subplot(1,1,1, axisbg='b')
    cs = plt.contourf(Z,X+charging_pos,rAm*1.0e3,100,cmap=plt.cm.jet, extend="both")
    cbar = plt.colorbar()
    cbar.set_label(r' $2 \pi r A_{\theta} \times 10^3 [T m^2]$',fontsize=16)
    CS = plt.contour(Z,X+charging_pos,rAm*1.0e3)
    plt.clabel(CS, inline=1, fontsize=13, linewidth=3, colors='white')
    plt.xlabel('r [m]'); plt.ylabel('z [m]')

    # draw shapes of coils
    p1=patches.Rectangle((-r_f-tr_f/2.0,z_f-tz_f/2.0+charging_pos),tr_f,tz_f,edgecolor='red',facecolor='white')
    p2=patches.Rectangle((r_f-tr_f/2.0,z_f-tz_f/2.0+charging_pos),tr_f,tz_f,edgecolor='red',facecolor='white')
    axxx.add_patch(p1); axxx.add_patch(p2)
    
    p1=patches.Rectangle((-r_c_opt-tr_c/2.0,z_c_opt-tz_c/2.0+charging_pos),tr_c,tz_c,edgecolor='red',facecolor='white')
    p2=patches.Rectangle((r_c_opt-tr_c/2.0,z_c_opt-tz_c/2.0+charging_pos),tr_c,tz_c,edgecolor='red',facecolor='white')
    axxx.add_patch(p1); axxx.add_patch(p2)

    if plot_Bs:
        fig22 = plt.figure()
        axxx=fig22.add_subplot(1,1,1, axisbg='b')
        cs = plt.contourf(Z,X+charging_pos,Bs,100,cmap=plt.cm.jet, extend="both")
        cbar = plt.colorbar()
        cbar.set_label(r'$B_{tot} [T]$', fontsize=16)
        CS = plt.contour(Z,X+charging_pos,rAm*1.0e3)
        plt.clabel(CS, inline=1, fontsize=13, linewidth=3, colors='white')
        plt.xlabel('r [m]'); plt.ylabel('z [m]')

        # draw shapes of coils
        p1=patches.Rectangle((-r_f-tr_f/2.0,z_f-tz_f/2.0+charging_pos),tr_f,tz_f,edgecolor='red',facecolor='white')
        p2=patches.Rectangle((r_f-tr_f/2.0,z_f-tz_f/2.0+charging_pos),tr_f,tz_f,edgecolor='red',facecolor='white')
        axxx.add_patch(p1); axxx.add_patch(p2)
    
        p1=patches.Rectangle((-r_c_opt-tr_c/2.0,z_c_opt-tz_c/2.0+charging_pos),tr_c,tz_c,edgecolor='red',facecolor='white')
        p2=patches.Rectangle((r_c_opt-tr_c/2.0,z_c_opt-tz_c/2.0+charging_pos),tr_c,tz_c,edgecolor='red',facecolor='white')
        axxx.add_patch(p1); axxx.add_patch(p2)
else:
    I_c_opt = I_c
    r_c_opt = r_c
    z_c_opt = z_c


if get_contour_vol:
    vol_yn=str(raw_input('Do you wish to compute confinement volumes in this configuration? (y/n) \n'))
    if vol_yn=='y':
        params=[I_l_opt,I_f_opt],[z_l_opt,z_f],[r_l_opt,r_f],[tz_l,tz_f],[tr_l,tr_f],[nz_l,nz_f],[nr_l,nr_f]
        get_flux_vol(params)
    else:
        print "Skipping computation of confinement volumes"

#########

if plot_opt1:
    plt.figure()
    plt.plot(range(opt1_cnt),tot_list,'ko', label='tot')
    plt.plot(range(opt1_cnt),term1_list,'ro', label='term #1')
    plt.plot(range(opt1_cnt),term2_list,'bo', label='term #2')
    plt.plot(range(opt1_cnt),term3_list,'go', label='term #3')
    plt.plot(range(opt1_cnt),term4_list,'mo', label='term #4')
    plt.plot(range(opt1_cnt),term5_list,'co', label='term #5')
    #plt.ylim([-1.5e-5,1e-4])
    plt.legend().draggable()
    
############################################
############################################

'''Method 1: assume you don't know that flux is proportional to current. What value of If matches the 
flux produced by the C-coil at the mid-radius of the F-coil during charging? '''

get_persistent_If_opt=False
if get_persistent_If_opt:
    
    # flux from C-coil
    X,Z,Bxm,Bzm,Bs,rAm = multicoil_fields([I_c_opt],[z_c_opt],[r_c_opt],[tz_c],[tr_c],[nz_c],[nr_c])
    xidx=np.argmin(np.abs(X[0,:]-r_f))
    zidx=np.argmin(np.abs(Z[:,0]-z_f))
    rAm_0=rAm[xidx,zidx]

    def persistent_f(x,grad,rAm_0):
        If=x
        X_f,Z_f,Bxm_f,Bzm_f,Bs_f,rAm_f = multicoil_fields([If],[z_f],[r_f],[tz_f],[tr_f],[nz_f],[nr_f])
        xidx_f=np.argmin(np.abs(X_f[0,:]-r_f))
        zidx_f=np.argmin(np.abs(Z_f[:,0]-z_f))
        rAm_1=rAm_f[xidx_f,zidx_f]

        print abs(rAm_0-rAm_1)
        return abs(rAm_0-rAm_1)

    opt_If = nlopt.opt(nlopt.LN_SBPLX, 1)
    opt_If.set_lower_bounds([1.0e3,])
    opt_If.set_upper_bounds([1.0e4,])
    opt_If.set_xtol_abs(1e-3)

    object_f=lambda x,grad: persistent_f(x,grad,rAm_0)

    opt_If.set_min_objective(object_f)
    guess_f=[I_f]
    If_pers=opt_If.optimize(guess_f)

    print "Persistent I_f: ", If_pers
 

'''  Method 2: the flux produced by a coil is directly proportional to its current. Use this fact to infer
the persistent If by taking the ratio of the flux produced at the location 'freezing_pos' (also inferred within an optimization) that gives perfect flux freezing.'''

if get_persistent_If:
    # Measured value:
    Bz_max=32.3e-3
    If_meas=If_from_Bz_onaxis(Bz_max,r_f,z_f,0,0)

    def freezing_pos(x,grad):
        f_factor = x[0]
        
        # Compute using flux freezing:
        X,Z,Bxm,Bzm,Bs,rAm = multicoil_fields([I_c_opt],[z_c_opt],[r_c_opt],[tz_c],[tr_c],[nz_c],[nr_c])
        xidx=np.argmin(np.abs(X[0,:]-z_f))
        zidx=np.argmin(np.abs(Z[:,0]-r_f+f_factor*tr_f/2.0))
        rAm_0=rAm[:,xidx]

        X_f,Z_f,Bxm_f,Bzm_f,Bs_f,rAm_f = multicoil_fields([I_f],[z_f],[r_f],[tz_f],[tr_f],[nz_f],[nr_f])
        xidx_f=np.argmin(np.abs(X_f[0,:]-z_f))
        zidx_f=np.argmin(np.abs(Z_f[:,0]-r_f+f_factor*tr_f/2.0))
        rAm_f=rAm_f[:,xidx_f]

        I_f_pers = I_f * rAm_0[zidx] / rAm_f[zidx_f]

        diff=abs(If_meas-I_f_pers)
        print "Flux freezing: ", f_factor, If_meas,I_f_pers, diff
        return diff

    opt_If = nlopt.opt(nlopt.LN_SBPLX, 1)
    opt_If.set_lower_bounds([-3.0,])
    opt_If.set_upper_bounds([3.0,])
    opt_If.set_xtol_abs(1e-3)

    opt_If.set_min_objective(freezing_pos)
    guess_factor=[0.3]
    freeze_factor=opt_If.optimize(guess_factor)

    # Use this freezing factor for final value of If:
    X,Z,Bxm,Bzm,Bs,rAm = multicoil_fields([I_c_opt],[z_c_opt],[r_c_opt],[tz_c],[tr_c],[nz_c],[nr_c])
    xidx=np.argmin(np.abs(X[0,:]-z_f))
    zidx=np.argmin(np.abs(Z[:,0]-r_f+freeze_factor*tr_f/2.0))
    rAm_0=rAm[:,xidx]

    X_f,Z_f,Bxm_f,Bzm_f,Bs_f,rAm_f = multicoil_fields([I_f],[z_f],[r_f],[tz_f],[tr_f],[nz_f],[nr_f])
    xidx_f=np.argmin(np.abs(X_f[0,:]-z_f))
    zidx_f=np.argmin(np.abs(Z_f[:,0]-r_f+freeze_factor*tr_f/2.0))
    rAm_f=rAm_f[:,xidx_f]

    I_f_pers = I_f * rAm_0[zidx] / rAm_f[zidx_f]

    plt.figure()
    plt.plot(Z[:,0],rAm_0*1.0e3,'b-')
    plt.plot(Z_f[:,0],rAm_f*rAm_0[zidx]/rAm_f[zidx_f]*1.0e3,'g-') 

    plt.xlim([0,0.3])
    ax=plt.gca()
    plt.plot([r_f-tr_f/2.0,r_f-tr_f/2.0],[0.0,ax.get_ylim()[1]],color='k',linestyle='--')
    plt.plot([r_f,r_f],[0.0,ax.get_ylim()[1]],'k')
    plt.plot([r_f+tr_f/2.0,r_f+tr_f/2.0],[0.0,ax.get_ylim()[1]],color='k',linestyle='--')
    plt.xlabel('r [m]')
    plt.ylabel(r'$\Psi$ [mWb]')

    print "****************"
    print "Freeze factor: ", freeze_factor[0]
    print "Persistent I_f (1st method): ", I_f_pers
    print "Compare this persistent current with value from measurement of Bz on axis"
    print "max(B_z) = %f mT"%Bz_max
    print "If_meas = ", If_meas
    print "****************"

'''
Method 3: this is similar to Method 2, but assumes that flux freezing occurs at outer radius of the F-coil.
Note that Method 2 method instead infers the radius at which flux freezing occurs by optimizing the radial location that seems consistent with the maximum Bz measured on-axis.
'''

if get_persistent_If_2:
    # Measured value:
    Bz_max=32.3e-3
    If_meas=If_from_Bz_onaxis(Bz_max,r_f,z_f,0,0)

    # Compute using flux freezing:
    X,Z,Bxm,Bzm,Bs,rAm = multicoil_fields([I_c_opt],[z_c_opt],[r_c_opt],[tz_c],[tr_c],[nz_c],[nr_c])
    xidx=np.argmin(np.abs(X[0,:]-z_f))
    zidx=np.argmin(np.abs(Z[:,0]-r_f+tr_f/2.0))
    rAm_0=rAm[:,xidx]

    X_f,Z_f,Bxm_f,Bzm_f,Bs_f,rAm_f = multicoil_fields([I_f],[z_f],[r_f],[tz_f],[tr_f],[nz_f],[nr_f])
    xidx_f=np.argmin(np.abs(X_f[0,:]-z_f))
    zidx_f=np.argmin(np.abs(Z_f[:,0]-r_f+tr_f/2.0))
    rAm_f=rAm_f[:,xidx_f]

    I_f_pers = I_f * rAm_0[zidx] / rAm_f[zidx_f]
    
    plt.figure()
    plt.plot(Z[:,0],rAm_0*1.0e3,'b-')
    plt.plot(Z_f[:,0],rAm_f* rAm_0[zidx] / rAm_f[zidx_f]*1.0e3,'g-') #* rAm_0[zidx] / rAm_f[zidx_f]

    plt.xlim([0,0.3])
    ax=plt.gca()
    plt.plot([r_f-tr_f/2.0,r_f-tr_f/2.0],[0.0,ax.get_ylim()[1]],color='k',linestyle='--')
    plt.plot([r_f,r_f],[0.0,ax.get_ylim()[1]],'k')
    plt.plot([r_f+tr_f/2.0,r_f+tr_f/2.0],[0.0,ax.get_ylim()[1]],color='k',linestyle='--')
    plt.xlabel('r [m]')
    plt.ylabel(r'$\Psi$ [mWb]')

    print "****************"
    print "Persistent I_f (2nd method): ", I_f_pers
    print "Compare this persistent current with value from measurement of Bz on axis"
    print "max(B_z) = %f mT"%Bz_max
    print "If_meas = ", If_meas
    print "****************"


''' Estimate reduction in I_f from flux freezing when close to powered L-coil '''

# Measured value:
Bz_max=32.3e-3
If1 = If_from_Bz_onaxis(Bz_max,r_f,z_f,0,0)
#If1=20000
#z_l1=0.2

# Compute using flux freezing:
X,Z,Bxm,Bzm,Bs,rAm = multicoil_fields([I_l_opt],[z_l1],[r_l_opt],[tz_l],[tr_l],[nz_l],[nr_l])
xidx=np.argmin(np.abs(X[0,:]-z_f))
zidx=np.argmin(np.abs(Z[:,0]-r_f+tr_f/2.0))
rAm_0=rAm[:,xidx]

X_f,Z_f,Bxm_f,Bzm_f,Bs_f,rAm_f = multicoil_fields([I_f],[z_f],[r_f],[tz_f],[tr_f],[nz_f],[nr_f])
xidx_f=np.argmin(np.abs(X_f[0,:]-z_f))
zidx_f=np.argmin(np.abs(Z_f[:,0]-r_f+tr_f/2.0))
rAm_f=rAm_f[:,xidx_f]

frac= rAm_0[zidx] / rAm_f[zidx_f]
If2= frac* If1

print " ******************* "
print "If before: ", If1
print "If after: ", If2


   

    
