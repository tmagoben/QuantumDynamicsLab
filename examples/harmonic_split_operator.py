import numpy as np
from qdyn_lab import scalar_step
from qdyn_lab.diagnostics import norm_scalar,scalar_energy
x=np.linspace(-12,12,1024,endpoint=False);m=w=1.;dt=.002;V=.5*m*w*w*x*x;x0=-2.;p0=1.4;psi=np.exp(-.5*m*w*(x-x0)**2+1j*p0*x);psi/=np.sqrt(norm_scalar(psi,x));E0=scalar_energy(psi,x,V,m)
for _ in range(1000): psi=scalar_step(psi,x,V,m,dt)
print('norm',norm_scalar(psi,x));print('energy drift',scalar_energy(psi,x,V,m)-E0)
