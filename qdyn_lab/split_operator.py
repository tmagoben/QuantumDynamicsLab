import numpy as np
from .diagnostics import uniform_dx
def _kinetic_half(psi,x,mass,dt):
 dx=uniform_dx(x); k=2*np.pi*np.fft.fftfreq(len(x),d=dx); phase=np.exp(-1j*(k*k/(2*mass))*dt/2); return np.fft.ifft(phase[:,None]*np.fft.fft(psi,axis=0),axis=0) if psi.ndim==2 else np.fft.ifft(phase*np.fft.fft(psi))
def scalar_step(psi,x,V,mass,dt):
 psi=_kinetic_half(np.asarray(psi,complex),x,mass,dt); psi=np.exp(-1j*np.asarray(V)*dt)*psi; return _kinetic_half(psi,x,mass,dt)
def multistate_step(psi,x,V,mass,dt):
 psi=np.asarray(psi,complex); V=np.asarray(V,complex)
 if V.shape!=(len(x),psi.shape[1],psi.shape[1]): raise ValueError('V must be (Ngrid,Nstate,Nstate)')
 if not np.allclose(V,V.conj().transpose(0,2,1),atol=1e-12): raise ValueError('potential matrices must be Hermitian')
 psi=_kinetic_half(psi,x,mass,dt); e,U=np.linalg.eigh(V); coeff=np.einsum('xmi,xm->xi',U.conj(),psi); coeff*=np.exp(-1j*e*dt); psi=np.einsum('xmi,xi->xm',U,coeff); return _kinetic_half(psi,x,mass,dt)
