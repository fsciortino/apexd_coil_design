import numpy as np
from scipy.special import ellipk, ellipe
from mpmath import ellippi

import matplotlib.pyplot as plt
plt.ion()
import mpmath

mu_0=8.85418782e-12
Is=1.0e3
L=0.4
a0=0.05

def kk(rho,xi):
    return np.sqrt(4*a0*rho/((a0+rho)**2+xi**2))

def hh(rho):
    return np.sqrt(4*a0*rho/((a0+rho)**2))

def xxi(z,sign): # sign =[-1,+1]
    return z+sign*L/2.0

Bpol=lambda rho,z,k,h: ((mu_0*Is)/(4*np.pi))*(1.0/L)*np.sqrt(a0/rho)*(((k**2-2)/(k))*ellipk(k**2)+(2/k)*ellipe(k**2))

Bz=lambda rho,z,k,h: -((mu_0*Is)/(4*np.pi))*(1.0/(2*L))*np.sqrt(1.0/(a0*rho))*(xi*k*(ellipk(k**2)+(a0-rho)/(a0+rho)*ellippi(h**2,k**2)))


zz=np.linspace(-0.3,0.3,50)
xx=np.linspace(-0.25,0.25,60)

BBpol=np.zeros((len(xx),len(zz)))
BBz=np.zeros((len(xx),len(zz)))

for iz in range(len(zz)):
    z=zz[iz]
    for ix in range(len(xx)):
        x=xx[ix]
        for sign in [1,-1]:
            xi=xxi(z,sign)
            k=kk(np.abs(x),xi)
            h=hh(np.abs(x))
            BBpol[ix,iz]+=sign*Bpol(np.abs(x),z,k,h)
            BBz[ix,iz]+=sign*Bz(np.abs(x),z,k,h)

plt.figure()
X,Z=np.meshgrid(zz,xx)
plt.contourf(X,Z,BBz)
plt.colorbar()
plt.title(r'$B_{z} [T]$')
plt.xlabel('x [m]')
plt.ylabel('y [m]')

plt.figure()
plt.contourf(X,Z,BBpol)
plt.colorbar()
plt.title(r'$B_{pol} [T]$')
plt.xlabel('x [m]')
plt.ylabel('y [m]')
