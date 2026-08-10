import numpy as np
from qdyn_lab import scalar_step,multistate_step
from qdyn_lab.diagnostics import norm_scalar,norm_multistate,scalar_energy,phase_aligned_error
def initial(x):
 p=np.exp(-.5*(x+1.3)**2+.7j*x); return p/np.sqrt(norm_scalar(p,x))
def test_scalar_norm_and_energy():
 x=np.linspace(-12,12,512,endpoint=False);V=.5*x*x;psi=initial(x);n0=norm_scalar(psi,x);E0=scalar_energy(psi,x,V,1.)
 for _ in range(400): psi=scalar_step(psi,x,V,1.,.002)
 assert abs(norm_scalar(psi,x)-n0)<2e-12; assert abs(scalar_energy(psi,x,V,1.)-E0)<2e-6
def test_multistate_norm():
 x=np.linspace(-10,10,512,endpoint=False);V=np.zeros((len(x),2,2));V[:,0,0]=.01*x;V[:,1,1]=-.01*x;V[:,0,1]=V[:,1,0]=.004;psi=np.zeros((len(x),2),complex);psi[:,0]=initial(x);n0=norm_multistate(psi,x)
 for _ in range(100): psi=multistate_step(psi,x,V,50.,.01)
 assert abs(norm_multistate(psi,x)-n0)<3e-12
def propagate(dt,T):
 x=np.linspace(-10,10,512,endpoint=False);V=.5*x*x;psi=initial(x);steps=round(T/dt)
 for _ in range(steps): psi=scalar_step(psi,x,V,1.,dt)
 return x,psi
def test_second_order_convergence():
 x,r=propagate(.000625,.2);_,a=propagate(.005,.2);_,b=propagate(.0025,.2); assert phase_aligned_error(r,b,x)<phase_aligned_error(r,a,x)/3.0
