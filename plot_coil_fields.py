import numpy as np
from coil_field import *
import matplotlib.pyplot as plt
plt.ion()
import matplotlib.patches as patches
import sys
import scipy.interpolate as interp

# base parameters
ci_1,rs1_1, rs2_1, zs1_1, zs2_1 = 3e4, 66e-3, 84e-3, -9e-3, 9e-3 # floating coil
ci_2,rs1_2, rs2_2, zs1_2, zs2_2 = 3.5e4, 1e-1, 1.5e-1, 1e-1, 1.1e-1 # levitation coil

ci_1,rs1_1, rs2_1, zs1_1, zs2_1 = 3e4, 66e-3, 84e-3, -9e-3, 9e-3 # floating coil
ci_2,rs1_2, rs2_2, zs1_2, zs2_2 = 0.5e4, 1e-1, 1.5e-1, 1e-1, 1.1e-1 # levitation coil

ci_1,rs1_1, rs2_1, zs1_1, zs2_1 = 0.0, 66e-3, 84e-3, -9e-3, 9e-3 # floating coil
ci_2,rs1_2, rs2_2, zs1_2, zs2_2 = 1e4, 84e-3, 100e-3,-30e-3, -18e-3 # charging coil
ww_1=rs2_1-rs1_1; hh_1=zs2_1-zs1_1
ww_2=rs2_2-rs1_2; hh_2=zs2_2-zs1_2

# grid_defs = [nx,nz,xmin,xmax,zmin,zmax]
grid_defs=[50,50,-200e-3,200e-3,-200e-3,200e-3]
try:
    option=int(sys.argv[1])
except:
    option=1


fig1=plt.figure(figsize=(12,8))

# Floating coil
xm, zm, Bxm_f, Bzm_f, Bs_f, rAm_f= coil_fields(ci_1,rs1_1,rs2_1,zs1_1,zs2_1, *grid_defs)
Bxm_f=Bxm_f[1:np.argmax(xm),1:np.argmax(zm)]
Bzm_f=Bzm_f[1:np.argmax(xm),1:np.argmax(zm)]
Bs_f=Bs_f[1:np.argmax(xm),1:np.argmax(zm)]
rAm_f=rAm_f[1:np.argmax(xm),1:np.argmax(zm)]

# Levitation coil
xm_l, zm_l, Bxm_l, Bzm_l, Bs_l, rAm_l= coil_fields(ci_2,rs1_2,rs2_2,zs1_2,zs2_2, *grid_defs)
Bxm_l=Bxm_l[1:np.argmax(xm),1:np.argmax(zm)]
Bzm_l=Bzm_l[1:np.argmax(xm),1:np.argmax(zm)]
Bs_l=Bs_l[1:np.argmax(xm),1:np.argmax(zm)]
rAm_l=rAm_l[1:np.argmax(xm),1:np.argmax(zm)]

xm=xm[1:np.argmax(xm)]; zm=zm[1:np.argmax(zm)]
xm_l=xm_l[1:np.argmax(xm_l)]; zm_l=zm_l[1:np.argmax(zm_l)] # should be redundant

# total
Bxm=Bxm_f+Bxm_l; Bzm=Bzm_f+Bzm_l; Bs=Bs_f+Bs_l; rAm=rAm_f+rAm_l
xmd=np.linspace(min(xm),max(xm),1000)
zmd=np.linspace(min(zm),max(zm),1000)
X,Z=np.meshgrid(xm,zm)
Xd,Zd = np.meshgrid(xmd,zmd)
        
axx=fig1.add_subplot(1,1,1, axisbg='b')

if option==0:
    Bsd = interp.interp2d(xm,zm,Bs,kind='cubic')(xmd,zmd)
    
    cs = plt.contourf(Zd,Xd,Bsd,100,cmap=plt.cm.jet, extend="both")
    #plt.title('I = {:04.1f} A, rs1={:04.3f} m, rs2 ={:04.3f} m, zs1 ={:04.3f} m, zs2 ={:04.3f} m'.format(ci_1,rs1,rs2,zs1,zs2))
    cbar = plt.colorbar()
    cbar.set_label(r' $B_{tot}$', fontsize=16)

    CS = plt.contour(Z,X,Bs)
    plt.clabel(CS, inline=1, fontsize=13, linewidth=3, colors='white')

    
elif option==1:
    phi= 2*np.pi*rAm;
    phid = interp.interp2d(xm,zm,phi,kind='cubic')(xmd,zmd)
                
    cs = plt.contourf(Zd,Xd,phid,100,cmap=plt.cm.jet, extend="both")
    #plt.title('I = {:04.1f} A, rs1={:04.3f} m, rs2 ={:04.3f} m, zs1 ={:04.3f} m, zs2 ={:04.3f} m'.format(ci,rs1,rs2,zs1,zs2))
    cbar = plt.colorbar()
    cbar.set_label(r' $2 \pi r A_{\theta}$',fontsize=16)
    CS = plt.contour(Z,X,phi)
    plt.clabel(CS, inline=1, fontsize=13, linewidth=3, colors='white')

    plt.xlim([min(xm),max(xm)])
    plt.ylim([min(zm),max(zm)])
    plt.show()


p1=patches.Rectangle((-rs2_1,zs1_1),ww_1,hh_1,edgecolor='red',facecolor='white')
p2=patches.Rectangle((rs1_1,zs1_1),ww_1,hh_1,edgecolor='red',facecolor='white')
axx.add_patch(p1); axx.add_patch(p2)
p1=patches.Rectangle((-rs2_2,zs1_2),ww_2,hh_2,edgecolor='red',facecolor='white')
p2=patches.Rectangle((rs1_2,zs1_2),ww_2,hh_2,edgecolor='red',facecolor='white')
axx.add_patch(p1); axx.add_patch(p2)
    
