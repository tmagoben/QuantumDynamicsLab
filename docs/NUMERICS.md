# Numerical notes

The FFT split operator implements second-order Strang splitting,

$$
e^{-iH\Delta t}=
e^{-iT\Delta t/2}\,
e^{-iV\Delta t}\,
e^{-iT\Delta t/2}
+
\mathcal{O}(\Delta t^3),
$$

for $H = T + V$. The $\mathcal{O}(\Delta t^3)$ term is the local splitting error, so
the propagated state has second-order global convergence. For Hermitian $T$ and $V$,
each exponential factor is unitary; consequently, the split step preserves the norm
up to numerical roundoff even though the state still has time-discretization error.

Crank-Nicolson is also unitary for a Hermitian discrete Hamiltonian. RK4 is included
as a general-purpose ODE integrator but is not exactly unitary, so norm drift must be
monitored during long-time quantum propagation.
