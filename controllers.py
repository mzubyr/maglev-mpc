"""PID and MPC controllers for the discretised maglev model."""
import numpy as np
from scipy.optimize import lsq_linear
from model import AD, BD, TS, steady_state

LOST = 0.05  # if |position deviation| exceeds 50 mm the ball is considered lost


class MPC:
    """Linear MPC with input bounds, solved as a bounded least-squares problem.

    cost = sum_k  q-weighted (x_k - x_ss)^2  +  r * (u - u_ss)^2
    subject to |u| <= umax. Inputs after the control horizon Nc are held.
    """

    def __init__(self, Np=29, Nc=4, q=(1e6, 1e2, 0.0), r=1e-2, umax=1.0):
        self.Np, self.Nc, self.r, self.umax = Np, Nc, r, umax
        self.F = np.vstack([np.linalg.matrix_power(AD, i) for i in range(1, Np + 1)])
        G = np.zeros((3 * Np, Nc))
        for i in range(1, Np + 1):
            for j in range(i):
                G[3*(i-1):3*i, min(j, Nc-1)] += (np.linalg.matrix_power(AD, i-1-j) @ BD).ravel()
        self.Qh = np.kron(np.eye(Np), np.diag(np.sqrt(q)))
        self.M = np.vstack([self.Qh @ G, np.sqrt(r) * np.eye(Nc)])

    def control(self, x, ref):
        xs, uss = steady_state(ref)
        Xs = np.tile(xs, self.Np)
        rhs = np.concatenate([self.Qh @ (Xs - self.F @ x),
                              np.sqrt(self.r) * uss * np.ones(self.Nc)])
        sol = lsq_linear(self.M, rhs, bounds=(-self.umax, self.umax), method="bvls")
        return sol.x[0]  # receding horizon: apply only the first input


class PID:
    """Discrete PID with output saturation and conditional-integration anti-windup."""

    def __init__(self, gains=(-152.0, -190.0, -16.0), umax=1.0):
        self.Kp, self.Ki, self.Kd = gains
        self.umax = umax
        self.I, self.e_prev = 0.0, None

    def control(self, x, ref):
        e = ref - x[0]
        d = 0.0 if self.e_prev is None else (e - self.e_prev) / TS
        self.e_prev = e
        u_raw = self.Kp * e + self.Ki * (self.I + e * TS) + self.Kd * d
        u = float(np.clip(u_raw, -self.umax, self.umax))
        if u == u_raw:          # integrate only when not saturated
            self.I += e * TS
        return u


def simulate(controller, x0=(0.0, 0.0, 0.0), ref=0.0, steps=400):
    """Closed-loop simulation. Returns position (m) and voltage arrays."""
    x = np.array(x0, dtype=float)
    ys, us = [], []
    for _ in range(steps):
        u = controller.control(x, ref)
        x = AD @ x + BD.ravel() * u
        ys.append(x[0]); us.append(u)
        if abs(x[0]) > LOST:
            break
    return np.array(ys), np.array(us)
