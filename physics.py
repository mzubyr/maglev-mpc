"""Nonlinear maglev model and a numerical check of the linearisation.

The parameter values below were reverse-engineered from the linearised
matrices (they reproduce every entry). VERIFY them against the thesis
parameter table before quoting them.
"""
import numpy as np

M_BALL, K_FRIC, G = 0.1, 0.001, 9.81
L0, L1, A_GAP, R = 0.01, 0.02, 0.05, 1.0


def inductance(x1):
    return L1 + L0 / (1.0 + x1 / A_GAP)


def f(x, u):
    """Nonlinear dynamics, x = [position, velocity, current]."""
    x1, x2, x3 = x
    dx1 = x2
    dx2 = G - (K_FRIC / M_BALL) * x2 - (L0 * A_GAP * x3**2) / (2 * M_BALL * (A_GAP + x1)**2)
    dx3 = (-R * x3 + L0 * A_GAP * x2 * x3 / (A_GAP + x1)**2 + u) / inductance(x1)
    return np.array([dx1, dx2, dx3])


def equilibrium(x1e=0.05):
    x3e = np.sqrt(2 * M_BALL * G * (A_GAP + x1e)**2 / (L0 * A_GAP))
    return np.array([x1e, 0.0, x3e]), R * x3e


def numerical_jacobians(h=1e-6):
    """Central-difference Jacobians of f at the equilibrium."""
    xe, ue = equilibrium()
    A = np.zeros((3, 3))
    for j in range(3):
        d = np.zeros(3); d[j] = h
        A[:, j] = (f(xe + d, ue) - f(xe - d, ue)) / (2 * h)
    B = ((f(xe, ue + h) - f(xe, ue - h)) / (2 * h)).reshape(3, 1)
    return A, B


if __name__ == "__main__":
    from model import A as A_thesis, B as B_thesis
    An, Bn = numerical_jacobians()
    xe, ue = equilibrium()
    print("equilibrium state:", xe, " equilibrium voltage:", ue)
    print("numerical A:\n", np.round(An, 4))
    print("max |A_numerical - A_thesis| =", np.abs(An - A_thesis).max())
    print("max |B_numerical - B_thesis| =", np.abs(Bn - B_thesis).max())
