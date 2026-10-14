"""Magnetic suspension (maglev) model, linearised around an equilibrium.

States:  x1 = ball position (m), x2 = velocity (m/s), x3 = coil current (A)
Input:   u  = coil voltage (V)
All variables below are DEVIATIONS from the equilibrium (x1e = 0.05 m).
"""
import numpy as np
from scipy.signal import cont2discrete

# Linearised model (from the B.Eng thesis, Chapter 4)
A = np.array([[0.0,    1.0,     0.0],
              [196.2, -0.01,   -3.1312],
              [0.0,    12.5284, -40.0]])
B = np.array([[0.0], [0.0], [40.0]])
C = np.array([[1.0, 0.0, 0.0]])
D = np.array([[0.0]])

TS = 0.0064  # sampling time (s)


def discretize(ts=TS):
    """Exact zero-order-hold discretisation: x[k+1] = Ad x[k] + Bd u[k]."""
    Ad, Bd, _, _, _ = cont2discrete((A, B, C, D), ts)
    return Ad, Bd


AD, BD = discretize()


def steady_state(ref):
    """State and input needed to hold the ball at position deviation `ref`.

    Row 2 of A with zero velocity:  196.2*ref - 3.1312*i = 0
    Row 3 of A with di/dt = 0:      -40*i + 40*u = 0  ->  u = i
    """
    i_ss = A[1, 0] / -A[1, 2] * ref
    return np.array([ref, 0.0, i_ss]), i_ss
