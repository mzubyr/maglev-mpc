import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
from model import A, AD
from physics import numerical_jacobians
from controllers import MPC, PID, simulate


def test_open_loop_unstable():
    assert np.linalg.eigvals(A).real.max() > 13
    assert np.abs(np.linalg.eigvals(AD)).max() > 1


def test_linearisation_matches_nonlinear_model():
    An, _ = numerical_jacobians()
    assert np.abs(An - A).max() < 0.01


def test_both_controllers_stabilise_small_offset():
    for c in (PID(umax=1.0), MPC(umax=1.0)):
        ys, _ = simulate(c, x0=(0.001, 0, 0), steps=300)
        assert abs(ys[-1]) < 1e-5
