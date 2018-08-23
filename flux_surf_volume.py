# Function to obtain 2D area and 3D volume in APEX-D coil design within a specific flux surface indicated by the user

import matplotlib.pyplot as plt
import numpy as np
from coil_functions import *
import pdb

def area(vs):
	""" 
	Application of Green's theorem to get a lower bound on the area of a flux surface 
	"""
	a=0
	x0,y0=vs[0]
	for [x1,y1] in vs[1:]:
		dx=x1-x0
		dy=y1-y0
		a+=0.5*(y0*dx-x0*dy)
		x0=x1
		y0=y1
	return a

def get_flux_vol(params):
	''' Get 2D area and volume of flux surface '''
	[I_l_opt,I_f_opt],[z_l_opt,z_f],[r_l_opt,r_f],[tz_l,tz_f],[tr_l,tr_f],[nz_l,nz_f],[nr_l,nr_f] = params
	X,Z,Bxm,Bzm,Bs,rAm=multicoil_fields([I_l_opt,I_f_opt],[z_l_opt,z_f],[r_l_opt,r_f],[tz_l,tz_f],[tr_l,tr_f],[nz_l,nz_f],[nr_l,nr_f])
	if plt.fignum_exists(123): plt.clf()
	fig1 = plt.figure(123)
	axx=fig1.add_subplot(1,1,1, axisbg='b')
	phi= 2*np.pi*rAm;
	cs = plt.contourf(Z,X,phi*1e3,100,cmap=plt.cm.jet, extend="both")
	cbar=plt.colorbar()
	cbar.set_label(r' $2 \pi r A_{\theta} \times 10^3$',fontsize=16)
	levels=np.linspace(max(phi.flatten()*1.0e3)/5, max(phi.flatten()*1.0e3),100)
	CS = plt.contour(Z,X,phi*1.0e3,levels=levels,colors='w')
	plt.clabel(CS, inline=1, fontsize=13, linewidth=3, colors='white')
	plt.xlabel('r [m]'); plt.ylabel('z [m]')

	x_min,x_max=axx.get_xlim()
	y_min,y_max=axx.get_ylim()

	user_req=1
	while user_req==1:
		contour_lv=float(raw_input('Within which contour level do you want to compute the area and volume? \n'))
		print "Click to get 4 corners within which to compute contour properties"
		             
		coords=plt.ginput(4,show_clicks=True)
		plt.show(block=False)
		cc=np.asarray(coords)

		y_max=cc[:,1].max()
		y_min=cc[:,1].min()
		x_max=cc[:,0].max()
		x_min=cc[:,0].min()

		if plt.fignum_exists(222): plt.clf()
		fig2=plt.figure(222);
		y_max_idx=np.argmin(np.abs(X[0,:]-y_max))
		y_min_idx=np.argmin(np.abs(X[0,:]-y_min))
		x_max_idx=np.argmin(np.abs(Z[:,0]-x_max))
		x_min_idx=np.argmin(np.abs(Z[:,0]-x_min))

		contour_idx=np.argmin(np.abs(levels-contour_lv))

		plt.contourf(Z[x_min_idx:x_max_idx,y_min_idx:y_max_idx],X[x_min_idx:x_max_idx,y_min_idx:y_max_idx],phi[x_min_idx:x_max_idx,y_min_idx:y_max_idx]*1.0e3, levels=levels)
		cs = plt.contour(Z[x_min_idx:x_max_idx,y_min_idx:y_max_idx],X[x_min_idx:x_max_idx,y_min_idx:y_max_idx],phi[x_min_idx:x_max_idx,y_min_idx:y_max_idx]*1.0e3,levels=levels,colors='w')
		plt.clabel(cs, inline=1, fontsize=13, linewidth=3, colors='white')
		plt.xlabel('r [m]'); plt.ylabel('z [m]')

		contour=cs.collections[contour_idx]
		vs=contour.get_paths()[0].vertices
		cs = plt.contour(Z[x_min_idx:x_max_idx,y_min_idx:y_max_idx],X[x_min_idx:x_max_idx,y_min_idx:y_max_idx],phi[x_min_idx:x_max_idx,y_min_idx:y_max_idx]*1.0e3,levels=[levels[contour_idx],],colors='w', linewidths=4.0)

		a=area(vs)
        
                R=(np.mean(vs[:,0])**2+np.mean(vs[:,1])**2)**(0.5)
                print "R = ", R
                sub_area=tr_f*tz_f
		print "Contour of value "+str(levels[contour_idx])+" has (roughly) A="+str(round(a,4))
		print "V = "+str(round(2*np.pi*R*a,4)) + "m-3, " +str(round(2*np.pi*R*a,4)*1000) + " l"
		print "(lower estimate using Green's theorem)"
                print "Confinement volume: ", round(2*np.pi*R*a-sub_area,4)

		print "*****************"
		user_in=raw_input('Enter 1 to look at a different contour, or anything else to stop \n')
		if user_in=='1':
		    user_req=1
		else:
		    user_req=0
		print " - "*20
    
    
#for i in range(1,len(levels)):
#    contour=cs.collections[i]
#    vs=contour.get_paths()[0].vertices
#    a=area(vs)
#    print "contour of value "+str(levels[i])+" has area "+str(a)
   



    
# area calculation using Green's theorem gives lower estimates  
def test_area_calc():
    plt.figure()
    delta=0.001
    x=np.arange(-3.1,3.1,delta)
    y=np.arange(-3.1,3.1,delta)
    X,Y=np.meshgrid(x,y)
    r=np.sqrt(X**2+Y**2)

    levels=[1.0,2.0,3.0,4.0]
    levels=np.linspace(0.1,3,20)
    cs=plt.contour(X,Y,r,levels=levels)
    plt.clabel(cs,inline=1,fontsize=10)

    for i in range(len(levels)):
        contour=cs.collections[i]
        vs=contour.get_paths()[0].vertices
        a=area(vs)
        print "r= " +str(levels[i]) + ": a=" +str(a)
        print "Analytical A="+str(np.pi*levels[i]**2)
        print "r= " +str(levels[i])+"error: " + str(np.pi*levels[i]**2-a)
        print "r= " +str(levels[i])+"error fe: " + str((np.pi*levels[i]**2-a)/(np.pi*levels[i])**2)
        print "************"
        #r= 1.0: a=2.84106591043
        #r= 2.0: a=12.0016303219
        #r= 3.0: a=27.4123559723
	return
