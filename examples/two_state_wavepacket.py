import numpy as np
from qdyn_lab import multistate_step
from qdyn_lab.diagnostics import norm_multistate
x=np.linspace(-20,20,1024,endpoint=False); mass=100.; dt=.02; V=np.zeros((len(x),2,2)); V[:,0,0]=.01*x; V[:,1,1]=-.01*x; V[:,0,1]=V[:,1,0]=.005; psi=np.zeros((len(x),2),complex); psi[:,0]=np.exp(-.2*(x+7)**2+1j*2*x); psi/=np.sqrt(norm_multistate(psi,x));
for _ in range(500): psi=multistate_step(psi,x,V,mass,dt)
print('norm',norm_multistate(psi,x));print('populations',(x[1]-x[0])*np.sum(np.abs(psi)**2,axis=0))
