import numpy as np
from coil_field import *
import matplotlib.pyplot as plt
plt.ion()
import matplotlib.patches as patches
import sys
import scipy.interpolate as interp

# base parameters
ci1,rs1, rs2, zs1, zs2 = 3e4, 66e-3, 84e-3, -9e-3, 9e-3
ci1_2,rs1_2, rs2_2, zs1_2, zs2_2 = 3e4, 66e-3, 84e-3, -9e-3, 9e-3
ww=rs2-rs1; hh=zs2-zs1
# grid_defs = [nx,nz,xmin,xmax,zmin,zmax]
grid_defs=[50,50,-200e-3,200e-3,-200e-3,200e-3]
try:
        option=int(sys.argv[1])
except:
        option=0

# scan over ci1
if option==0 or option ==1:

	ci1_arr = np.linspace(1e3, 1e5, 2)
	Bs_scan = []

	fig1=plt.figure(figsize=(12,8))
	    
	for k in range(len(ci1_arr)):
	    ci=ci1_arr[k]
	    xm, zm, Bxm, Bzm, Bs, rAm= coil_fields(ci,rs1,rs2,zs1,zs2, *grid_defs)

            Bxm=Bxm[1:np.argmax(xm),1:np.argmax(zm)]
            Bzm=Bzm[1:np.argmax(xm),1:np.argmax(zm)]
            Bs=Bs[1:np.argmax(xm),1:np.argmax(zm)]
            rAm=rAm[1:np.argmax(xm),1:np.argmax(zm)]
            xm=xm[1:np.argmax(xm)]; zm=zm[1:np.argmax(zm)]

            xmd=np.linspace(min(xm),max(xm),1000)
            zmd=np.linspace(min(zm),max(zm),1000)
            X,Z=np.meshgrid(xm,zm)
            Xd,Zd = np.meshgrid(xmd,zmd)
        
	    axx=fig1.add_subplot(len(ci1_arr),1,k, axisbg='b')

	    if option==0:
                Bsd = interp.interp2d(xm,zm,Bs,kind='cubic')(xmd,zmd)
                cs = plt.contourf(Zd,Xd,Bsd,100,cmap=plt.cm.jet, extend="both")
                plt.title('I = {:04.1f} A, rs1={:04.3f} m, rs2 ={:04.3f} m, zs1 ={:04.3f} m, zs2 ={:04.3f} m'.format(ci,rs1,rs2,zs1,zs2))
                cbar = plt.colorbar()
                cbar.set_label(r' $B_{tot}$', fontsize=16)

                CS = plt.contour(Z,X,Bs)
                plt.clabel(CS, inline=1, fontsize=13, linewidth=3, colors='white')
	    elif option==1:
                phi= 2*np.pi*rAm;
                phid = interp.interp2d(xm,zm,phi,kind='cubic')(xmd,zmd)
                
	    	cs = plt.contourf(Zd,Xd,phid,100,cmap=plt.cm.jet, extend="both")
                plt.title('I = {:04.1f} A, rs1={:04.3f} m, rs2 ={:04.3f} m, zs1 ={:04.3f} m, zs2 ={:04.3f} m'.format(ci,rs1,rs2,zs1,zs2))
                cbar = plt.colorbar()
                cbar.set_label(r' $2 \pi r A_{\theta}$',fontsize=16)
                CS = plt.contour(Z,X,phi)
                plt.clabel(CS, inline=1, fontsize=13, linewidth=3, colors='white')

            p1=patches.Rectangle((-rs2,zs1),ww,hh,edgecolor='red',facecolor='white')
            p2=patches.Rectangle((rs1,zs1),ww,hh,edgecolor='red',facecolor='white')
            axx.add_patch(p1); axx.add_patch(p2)
        
	    plt.xlim([min(xm),max(xm)])
	    plt.ylim([min(zm),max(zm)])
	    plt.show()
	    
            
	    Bs_scan.append(max(Bs.flatten()))

	#plt.clim([0,max(Bs.flatten())])
	#plt.figure()
	#plt.plot(ci1_arr, Bs_scan)
	#plt.show()

# make quiver plot  
if option==2:   
	ci=30000
	xm, zm, Bxm, Bzm, Bs, rAm= coil_fields(ci,rs1,rs2,zs1,zs2, *grid_defs)
        Bxm=Bxm[:np.argmax(xm),:][:,:np.argmax(zm)]
        Bzm=Bzm[:np.argmax(xm),:][:,:np.argmax(zm)]
        Bs=Bs[:np.argmax(xm),:][:,:np.argmax(zm)]
        rAm=rAm[:np.argmax(xm),:][:,:np.argmax(zm)]
        xm=xm[:np.argmax(xm)]; zm=zm[:np.argmax(xm)]
        X, Z = np.meshgrid(xm, zm)
	fig = plt.figure()
        ax = fig.add_subplot(111)
	ax.quiver(X,Z,Bxm,Bzm)

        p1=patches.Rectangle((-rs2,zs1),ww,hh)
        p2=patches.Rectangle((rs1,zs1),ww,hh)
        ax.add_patch(p1); ax.add_patch(p2)
	plt.xlim([min(xm),max(xm)])
	plt.ylim([min(zm),max(zm)])
	plt.show()    


