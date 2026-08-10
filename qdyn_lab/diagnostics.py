import numpy as np
def uniform_dx(x):
 dx=np.diff(x);
 if not np.allclose(dx,dx[0]): raise ValueError('grid must be uniform')
 return float(dx[0])
def norm_scalar(psi,x): return float(uniform_dx(x)*np.sum(np.abs(psi)**2))
def norm_multistate(psi,x): return float(uniform_dx(x)*np.sum(np.abs(psi)**2))
def scalar_energy(psi,x,V,mass):
 dx=uniform_dx(x); k=2*np.pi*np.fft.fftfreq(len(x),d=dx); Tpsi=np.fft.ifft((k*k/(2*mass))*np.fft.fft(psi)); return float((dx*np.vdot(psi,Tpsi+V*psi)).real)
def phase_aligned_error(reference,candidate,x):
 dx=uniform_dx(x); ov=dx*np.vdot(reference,candidate); phase=1.0 if abs(ov)==0 else ov/abs(ov); return float(np.sqrt(dx*np.sum(np.abs(candidate-phase*reference)**2)))
