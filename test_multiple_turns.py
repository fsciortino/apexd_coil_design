import coil_turns_fields
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interp
plt.ion()

grid_defs=[50,50,-300e-3,300e-3,-280e-3,280e-3]
nt_x=2
nt_z=2
I=1.0e3
rs1=0.04
rs2=0.045
zs1=-0.05
zs2=0.05
xm,zm,Bxm,Bzm,Bs,rAm=coil_turns_fields.coil_fields(I,rs1,rs2,zs1,zs2,nt_x,nt_z,*grid_defs)

Bxm=Bxm[1:np.argmax(xm),1:np.argmax(zm)]
Bzm=Bzm[1:np.argmax(xm),1:np.argmax(zm)]
Bs=Bs[1:np.argmax(xm),1:np.argmax(zm)]
rAm=rAm[1:np.argmax(xm),1:np.argmax(zm)]

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

fig=plt.figure()
axxx=fig.add_subplot(1,1,1, axisbg='b')
cs = plt.contourf(Zd,Xd,Bsd,100,cmap=plt.cm.jet, extend="both")
cbar = plt.colorbar()
cbar.set_label(r'$B_{tot} [T]$', fontsize=16)
CS = plt.contour(Zd,Xd,rAmd*1.0e3)
plt.clabel(CS, inline=1, fontsize=13, linewidth=3, colors='white')
plt.xlabel('r [m]'); plt.ylabel('z [m]')
        
