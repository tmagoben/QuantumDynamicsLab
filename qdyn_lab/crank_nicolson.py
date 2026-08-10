import numpy as np
from .diagnostics import uniform_dx
def periodic_hamiltonian(x,V,mass):
 N=len(x);dx=uniform_dx(x); lap=np.zeros((N,N));
 for i in range(N): lap[i,i]=-2; lap[i,(i-1)%N]=1; lap[i,(i+1)%N]=1
 T=-(1/(2*mass*dx*dx))*lap; return T+np.diag(np.asarray(V,float))
def step(psi,H,dt):
 H=np.asarray(H,complex); I=np.eye(H.shape[0],dtype=complex); return np.linalg.solve(I+0.5j*dt*H,(I-0.5j*dt*H)@psi)
