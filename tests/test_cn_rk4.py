import numpy as np
from qdyn_lab.crank_nicolson import periodic_hamiltonian,step as cn_step
from qdyn_lab.rk4 import step as rk4_step
from qdyn_lab.diagnostics import norm_scalar
def test_crank_nicolson_norm():
 x=np.linspace(-8,8,128,endpoint=False);V=.5*x*x;H=periodic_hamiltonian(x,V,1.);psi=np.exp(-x*x/2).astype(complex);psi/=np.sqrt(norm_scalar(psi,x));n0=norm_scalar(psi,x)
 for _ in range(30): psi=cn_step(psi,H,.002)
 assert abs(norm_scalar(psi,x)-n0)<1e-12
def test_rk4_converges_for_two_level_system():
 H=np.array([[.2,.03],[.03,-.1]],complex);psi=np.array([1.,0.],complex);T=.2
 def run(dt):
  p=psi.copy()
  for _ in range(round(T/dt)): p=rk4_step(p,H,dt)
  return p
 exact=np.linalg.eigh(H);e,U=exact;target=U@(np.exp(-1j*e*T)*(U.conj().T@psi)); assert np.linalg.norm(run(.01)-target)<2e-10
