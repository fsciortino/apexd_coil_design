import numpy as np
import matplotlib.pyplot as plt
from scipy.special import ellipk, ellipe
import matplotlib.ticker as ticker
import scipy.interpolate as interp
import pdb
import coil_turns_fields
plt.ion()

mu0=4.0*np.pi*10**(-7)
d_wire = 2.0e-3 # wire diameter

def k(Rc,Zc,r,z):
    return np.sqrt((4.0 *Rc *r)/((r+Rc)**2+(z-Zc)**2))

def EK(Rc,Zc,r,z):
    return ellipk(k(Rc,Zc,r,z)**2)

def EE(Rc,Zc,r,z):
    return ellipe(k(Rc,Zc,r,z)**2)

def phi(Ic,Rc,Zc,r,z):
    term1= (r*mu0*Ic/(np.pi*k(Rc,Zc,r,z)))*(Rc/r)**(0.5)
    term2 = (1.0-0.5*(k(Rc,Zc,r,z)**2))*EK(Rc,Zc,r,z)-EE(Rc,Zc,r,z)
    return term1+term2

def Br(Ic,Rc,Zc,r,z):
    return mu0*Ic/(2*np.pi)*(z-Zc)/(r*((r+Rc)**2+(z-Zc)**2)**(0.5))*(-EK(Rc,Zc,r,z)+(Rc**2+r**2+(z-Zc)**2)/((Rc-r)**2+(z-Zc)**2)*EE(Rc,Zc,r, z)) 

def Br_fv(Ic,Rc,Zc,r,z,nr,nz):
    '''Calculation of Br including finite spatial distribution of coil turns.
    Rc and Zc are the mean R and Z of a coil whose field we want to compute.
    Finite volume effects are taken into consideration by considering nr and nz radial and vertical loops, respectively.'''
    num_turns=nz*nr
    Ic_per_turn = float(Ic) / num_turns # units: A (not AT)
    
    br = 0
    for iz in range(nz):
        # Vertical position of turn (recall that z_l and r_l are positions of coil center)
        z_turn = Zc - d_wire * float(nz)/2.0 + iz * d_wire
        
        for ir in range(nr):
            # radial position of turn
            r_turn = Rc - d_wire * float(nr)/2.0 + ir * d_wire
            
            br+=Br(Ic_per_turn,r_turn,z_turn,r,z)
            
    return br
            
def Bz(Ic,Rc,Zc,r,z):
    return mu0*Ic/(2.0*np.pi*((r+Rc)**2+(z-Zc)**2)**(0.5))*(EK(Rc,Zc,r,z)+(Rc**2 -r**2-(z-Zc)**2)/((Rc-r)**2+(z-Zc)**2)*EE(Rc,Zc,r,z))

def Bz_fv(Ic,Rc,Zc,r,z,nr,nz):
    '''Calculation of Bz including finite spatial distribution of coil turns  '''
    num_turns=nz*nr
    Ic_per_turn = float(Ic) / num_turns # units: A (not AT)
    
    bz = 0
    for iz in range(nz):
        # Vertical position of turn (recall that z_l and r_l are positions of coil center)
        z_turn = Zc - d_wire * float(nz)/2.0 + iz * d_wire
        
        for ir in range(nr):
            # radial position of turn
            r_turn = Rc - d_wire * float(nr)/2.0 + ir * d_wire
            
            bz+=Bz(Ic_per_turn,r_turn,z_turn,r,z)
            
    return bz

def If_from_Bz_onaxis(Bz,Rc,Zc,r,z):
    return Bz*(mu0/(2.0*np.pi*((r+Rc)**2+(z-Zc)**2)**(0.5))*(EK(Rc,Zc,r,z)+(Rc**2 -r**2-(z-Zc)**2)/((Rc-r)**2+(z-Zc)**2)*EE(Rc,Zc,r,z)))**(-1)

def M(RA,RB,z1,z):
    return mu0 * (Ra*RB)**(0.5) * ((2.0/k(Rc,Zc,r,z)-k(Rc,Zc,r,z))*EK(Rc,Zc,r,z)-2.0/k(Rc,Zc,r,z)*EK(Rc,Zc,r,z))


def dMdz(RA,RB,z1,z):
    return (-mu0*(z-z1))/(4.0*np.sqrt(RA*RB))*(-2.0*k(RA,z1,RB,z)*EK(RA,z1,RB,z)+k(RA,z1,RB,z)*(2.0-k(RA,z1,RB,z)**2)/(1.0-k(RA,z1,RB,z)**2)*EE(RA,z1,RB,z))


def dBrdz(Ic,Rc,Zc,r,z):
    delta=1e-6
    return (Br(Ic,Rc,Zc,r,z+delta)-Br(Ic,Rc,Zc,r,z))/delta

def dBrdz_fv(Ic,Rc,Zc,r,z,nr,nz):
    delta=1e-6
    return (Br_fv(Ic,Rc,Zc,r,z+delta,nr,nz)-Br_fv(Ic,Rc,Zc,r,z,nr,nz))/delta


def dBzdr(Ic,Rc,Zc,r,z):
    delta=1e-6
    return (Bz(Ic,Rc,Zc,r+delta,z)-Bz(Ic,Rc,Zc,r,z))/delta

def dBzdr_fv(Ic,Rc,Zc,r,z,nr,nz):
    delta=1e-6
    return (Bz_fv(Ic,Rc,Zc,r+delta,z,nr,nz)-Bz_fv(Ic,Rc,Zc,r,z,nr,nz))/delta

               
def coil_force(R1,R2,I1,I2,z1,z2):
    '''# 1 is coil which feels force from external field (from #2)'''
    return - 2*np.pi*R1*I1*Br(I2,R2,z2,R1,z1) #inverted

def coil_force_fv(R1,R2,I1,I2,z1,z2,nr2,nz2):
    ''' Force on current loop from coil of finite volume/number of turns.
    # 1 is coil which feels force from external field (from #2) '''
    return - 2*np.pi*R1*I1*Br_fv(I2,R2,z2,R1,z1,nr2,nz2) #inverted


def multicoil_fields(I,z,r,tz,tr,nz,nr):
    ''' Superpose magnetic fields from multiple coils. 
    Parameters:
    - I: list or array, total currents
    - z: list or array, center vertical positions
    - r: list or array, major radii 
    - tz: list or array, thicknesses in the z-direction
    - tr: list or array, thicknesses in the r-direction
    - nz: list or array, number of turns in the z-direction
    - nr: list or array, number of turns in the r-direction
    '''
   
    # [#grid points x,z; y-scale lims, x-scale lims] # Keep symmetric
    grid_defs=[50,50,-300e-3,300e-3,-280e-3,280e-3]
    #grid_defs=[200,200,-200e-3,200e-3,-300e-3,300e-3]
    I=np.asarray(I); z=np.asarray(z); r=np.asarray(r)
    tz=np.asarray(tz); tr=np.asarray(tr)
    
    # All outputs of coil_fields() are linear, so can be added
    # 1st coil
    zs1 = z[0] - tz[0]/2.0; zs2 = z[0] + tz[0]/2.0
    rs1 = r[0] - tr[0]/2.0; rs2 = r[0] + tr[0]/2.0
    
    xm,zm,Bxm,Bzm,Bs,rAm=coil_turns_fields.coil_fields(I[0],rs1,rs2,zs1,zs2,nz[0],nr[0],*grid_defs)
    Bxm=Bxm[1:np.argmax(xm),1:np.argmax(zm)]
    Bzm=Bzm[1:np.argmax(xm),1:np.argmax(zm)]
    Bs=Bs[1:np.argmax(xm),1:np.argmax(zm)]
    rAm=rAm[1:np.argmax(xm),1:np.argmax(zm)]

    # add any other coil
    if I.shape[0]>1:
        for j in range(1,I.shape[0]):
            zs1 = z[j] - tz[j]/2.0; zs2 = z[j] + tz[j]/2.0
            rs1 = r[j] - tr[j]/2.0; rs2 = r[j] + tr[j]/2.0
    
            xm,zm,Bxm2,Bzm2,Bs2,rAm2=coil_turns_fields.coil_fields(I[j],rs1,rs2,zs1,zs2,nz[j],nr[j],*grid_defs)
            Bxm2=Bxm2[1:np.argmax(xm),1:np.argmax(zm)]
            Bzm2=Bzm2[1:np.argmax(xm),1:np.argmax(zm)]
            Bs2=Bs2[1:np.argmax(xm),1:np.argmax(zm)]
            rAm2=rAm2[1:np.argmax(xm),1:np.argmax(zm)]

            # total
            Bxm=Bxm+Bxm2; Bzm=Bzm+Bzm2; Bs=Bs+Bs2; rAm=rAm+rAm2

    xm=xm[1:np.argmax(xm)]; zm=zm[1:np.argmax(zm)]
    
    # Create dense grid:
    xmd=np.linspace(min(xm),max(xm),900)
    zmd=np.linspace(min(zm),max(zm),1000)
    X,Z=np.meshgrid(xm,zm)
    Xd,Zd = np.meshgrid(xmd,zmd)

    # Interpolate fields on dense grid
    Bxmd = interp.interp2d(xm,zm,Bxm,kind='cubic')(xmd,zmd)
    Bzmd = interp.interp2d(xm,zm,Bzm,kind='cubic')(xmd,zmd)
    Bsd = interp.interp2d(xm,zm,Bs,kind='cubic')(xmd,zmd)
    rAmd = interp.interp2d(xm,zm,rAm,kind='cubic')(xmd,zmd)

    return Xd,Zd,Bxmd,Bzmd,Bsd,rAmd


def c_objective(x,grad,params_c,params_f,rAtheta,c_bounds, nu_c,fixed_params, fixed_param_values):
    """ Objective function for C-coil parameter optimization
        x: list or array -- Ic, rc and zc. Some of these parameters may be kept fixed for a specific optimization and specified in fixed_param_values. 
        grad: requirement for nlopt, leave empty, []
        params_c: parameters of C-coil that will be kept fixed -- zc,tzc,trc,nrc,nzc
        params_f: parameters of F-coil line-of-charge -- zf,rf,nrf,nzf
        rAtheta: float, value of r*A_{\theta} at major radius of F-coil (real objective)       
        nu_c: list or array of floats, regularization parameters
        fixed_params: list of bools, indicating which parameters are kept fixed
        fixed_param_values: list of floats, containing values of fixed parameters
    """
    if not fixed_params[0]: Ic_norm=x[0]
    else: Ic_norm=fixed_param_values[0]
    
    if not fixed_params[1]: rc_norm=x[1]
    else: rc_norm=fixed_param_values[1]
    
    if not fixed_params[2]: zc_norm=x[2]
    else: zc_norm=fixed_param_values[2]
    
    # Recover real units of guesses:
    Ic=Ic_norm*c_bounds[0][1]
    rc=rc_norm*c_bounds[1][1]
    
    # Recover sign of zc:
    zc= - zc_norm*c_bounds[2][1]

    # unwrap C and F-coil parameters
    tzc,trc,nzc,nrc = params_c
    zf,rf = params_f

    # get fields from C-coil parameters
    X,Z,Bxm,Bzm,Bs,rAm = multicoil_fields([Ic],[zc],[rc],[tzc],[trc],[nzc],[nrc])

    # find rAm at zf and rf position
    xidx=np.argmin(np.abs(X[0,:]-rf))
    zidx=np.argmin(np.abs(Z[:,0]-zf))
    
    rAm_xz=rAm[xidx,zidx]
    out = np.abs(rAtheta - rAm_xz)
    print out,nu_c[0]*Ic, Ic,rc,zc
    return out+nu_c[0]*Ic


def c_coil_opt(f_params,c_params_guess,c_bounds,nu_c,fixed_params):
    ''' Given current and dimension of an F-coil, calculate 
    appropriate current and dimensions of required C-coil. 
    Note that we approximate the F-coil by a line of charge here.

    Parameters:
    - f_params: list, [If,rf,zf,tz_f,tr_f,nzf,nrf]
    - c_params_guess: list, [Ic_norm,rc_norm,zc_norm,tzc,tcf,nzc,nrc]. 
    The latter are just guesses for the C-coil optimization
    - nu_c: list or array of floats, regularization parameters
    - fixed_params: lists of bools indicating which variables should be kept fixed in optimization
    '''
    try:
        import nlopt
    except:
        raise ValueError('Could not import nlopt. Try using a virtualenv')
    
    Ic_norm_in,rc_norm_in,zc_norm_in,tzc_in,trc_in,nzc,nrc=c_params_guess
    If,rf,zf,tzf,trf,nzf,nrf = f_params
        
    # Match flux in F-coil inner radius with flux from C-coil
    # this is a lower estimate of flux
    X,Z,Bxm,Bzm,Bs,rAm = multicoil_fields([If],[zf],[rf],[tzf],[trf],[nzf],[nrf])
    xidx=np.argmin(np.abs(X[0,:]-rf))
    zidx=np.argmin(np.abs(Z[:,0]-zf))
    rAtheta_obj=rAm[xidx,zidx]
 
    params_c=tzc_in,trc_in,nzc,nrc
    params_f=zf,rf

    # Use fixed_params to determine which parameters should be kept free to vary
    n_fixed=sum(fixed_params)
    print "Number of free parameters: ", 3-n_fixed
    
    n_free=3 - n_fixed
    fixed_param_values=[Ic_norm_in,rc_norm_in,zc_norm_in]
    
    opt_c = nlopt.opt(nlopt.LN_SBPLX, n_free)

    # set lower bounds only for values to be optimized
    b_tmp=[c_bounds[0][0]/c_bounds[0][1],c_bounds[1][0]/c_bounds[1][1],c_bounds[2][0]/c_bounds[2][1]]
    
    lowerbounds=[i for indx,i in enumerate(b_tmp) if not fixed_params[indx]]
        
    opt_c.set_lower_bounds(lowerbounds)
    opt_c.set_upper_bounds([1.0,]*n_free)
    opt_c.set_ftol_abs(1e-5)
    opt_c.set_xtol_abs(1e-4)
    
    C_objective= lambda x,grad: c_objective(x,grad,params_c,params_f,rAtheta_obj,c_bounds,nu_c,fixed_params,fixed_param_values)
    opt_c.set_min_objective(C_objective)

    c_guess=np.asarray([i for indx,i in enumerate([Ic_norm_in,rc_norm_in,zc_norm_in]) if not fixed_params[indx]])

    # launch optimization
    res_c = opt_c.optimize(c_guess)

    # get optimization output:
    if not fixed_params[0]: I_c_opt=res_c[0]
    else: I_c_opt=fixed_param_values[0]
    if not fixed_params[1]: r_c_opt=res_c[1]
    else: r_c_opt=fixed_param_values[1]
    if not fixed_params[2]: z_c_opt=res_c[2]
    else: z_c_opt=fixed_param_values[2]
    
    return I_c_opt,r_c_opt, z_c_opt
