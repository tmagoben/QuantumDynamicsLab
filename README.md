# QuantumDynamicsLab

Transparent reference implementations for time-dependent quantum propagation.

## Methods

- second-order Strang split operator on a periodic FFT grid;
- multistate split operator with local diagonalization of a Hermitian potential matrix;
- Crank-Nicolson on a periodic second-order finite-difference grid;
- RK4 for finite-dimensional TDSE examples, with the explicit warning that RK4 is
  not exactly unitary;
- norm, energy, and phase-aligned convergence diagnostics.

```bash
pip install -e ".[dev]"
python examples/harmonic_split_operator.py
python examples/two_state_wavepacket.py
pytest -q
```

Atomic units and $\hbar=1$ are used.
