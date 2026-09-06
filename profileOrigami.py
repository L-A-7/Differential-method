# Generate an Origami profile file
# W profile creation with rounded edges


import numpy as np
import matplotlib.pyplot as plt


Nx=32
alpha=44.7*np.pi/180

R=30
L=200

Nr=round(R*sin(alpha)/L*Nx)

x1=np.arange(-Nr+1,Nr-1,1)*L/Nx
pr1=R-sqrt(R*R-x1*x1) + (sin(alpha)*tan(alpha)+cos(alpha)-1)*R

x2=np.arange(Nr,Nx/2-Nr,1)*L/Nx
pr2=x2*tan(alpha)

x3=np.arange(-Nr+1,Nr-1,1)*L/Nx
pr3=tan(alpha)*L/2 - (sin(alpha)*tan(alpha)+cos(alpha)-1)*R -R + sqrt(R*R-x3*x3)

x4=np.arange(Nr,Nx/2-Nr,1)*L/Nx
pr4=tan(alpha)*L/2 - x4*tan(alpha)

x=np.concatenate((x1, x2, x3+L/2, x4+L/2)) + (Nr-1)*L/Nx
pr=np.concatenate((pr1, pr2, pr3, pr4))

plt.plot(x,pr)


