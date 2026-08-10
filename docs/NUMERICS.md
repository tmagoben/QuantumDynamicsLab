# Numerical notes

The FFT split operator implements Strang splitting,

$$
 e^{-iH\Delta t} = e^{-iT\Delta t/2}\; e^{-iV\Delta t}\; e^{-iT\Delta t/2} + O(\Delta t^3),
$$

which gives second-order global convergence for the time-dependent Schrödinger equation with \(H=T+V\). Each factor is unitary for Hermitian operators \(T\) and \(V\), so norm conservation is preserved up to the time-stepping error.

Crank–Nicolson is also unitary for a Hermitian discrete Hamiltonian. RK4 is included as a general-purpose ODE integrator but is not exactly unitary, so it should be used with care for long-time quantum propagation.
