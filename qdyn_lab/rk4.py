import numpy as np
def step(psi,H,dt):
 H=np.asarray(H,complex); f=lambda y:-1j*H@y; k1=f(psi); k2=f(psi+dt*k1/2); k3=f(psi+dt*k2/2); k4=f(psi+dt*k3); return psi+dt*(k1+2*k2+2*k3+k4)/6
