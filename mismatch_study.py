"""Model-mismatch study: an MPC designed on the nominal maglev model is run on plants with
perturbed parameters. Reports 1 mm recovery, step-response steady-state error and the largest
recoverable offset. Self-contained; run with: python mismatch_study.py
"""
import numpy as np
from scipy.signal import cont2discrete
from scipy.optimize import lsq_linear

A0 = np.array([[0, 1, 0], [196.2, -0.01, -3.1312], [0, 12.5284, -40]])
B0 = np.array([[0], [0], [40.0]])
C = np.array([[1, 0, 0.]]); D = np.array([[0.]])
Ts = 0.0064

def disc(A, B):
    Ad, Bd, _, _, _ = cont2discrete((A, B, C, D), Ts)
    return Ad, Bd

Ad0, Bd0 = disc(A0, B0)              # nominal model: what the controller knows

# ---- MPC built from the nominal model ----
Np, Nc, umax = 29, 4, 0.3
q = np.array([1e6, 1e2, 0.0]); r = 1e-2
F = np.vstack([np.linalg.matrix_power(Ad0, i) for i in range(1, Np + 1)])
G = np.zeros((3 * Np, Nc))
for i in range(1, Np + 1):
    for j in range(i):
        G[3*(i-1):3*i, min(j, Nc-1)] += (np.linalg.matrix_power(Ad0, i-1-j) @ Bd0).ravel()
Qh = np.kron(np.eye(Np), np.diag(np.sqrt(q)))
M = np.vstack([Qh @ G, np.sqrt(r) * np.eye(Nc)])

def mpc(x, ref):
    i_ss = A0[1, 0] / -A0[1, 2] * ref          # steady state according to the NOMINAL model
    xs = np.array([ref, 0.0, i_ss])
    rhs = np.concatenate([Qh @ (np.tile(xs, Np) - F @ x), np.sqrt(r) * i_ss * np.ones(Nc)])
    return lsq_linear(M, rhs, bounds=(-umax, umax), method="bvls").x[0]

def simulate(Ad_true, Bd_true, x0=(0, 0, 0), ref=0.0, steps=400):
    """The controller uses the nominal model; the 'real' plant is Ad_true, Bd_true."""
    x = np.array(x0, float); ys = []
    for _ in range(steps):
        x = Ad_true @ x + Bd_true.ravel() * mpc(x, ref)
        ys.append(x[0])
        if abs(x[0]) > 0.05: break                 # ball lost
    return np.array(ys)

def evaluate(name, A, B):
    Adt, Bdt = disc(A, B)
    y = simulate(Adt, Bdt, x0=(0.001, 0, 0), steps=300)
    ok = len(y) == 300 and abs(y[-1]) < 1e-5
    ys = simulate(Adt, Bdt, ref=0.001, steps=400)
    err = "lost" if len(ys) < 400 else f"{(0.001 - ys[-1]) / 0.001 * 100:.1f}"
    best = 0.0
    for x0mm in np.arange(0.5, 8.01, 0.5):
        yy = simulate(Adt, Bdt, x0=(x0mm / 1000, 0, 0), steps=300)
        if len(yy) == 300 and abs(yy[-1]) < 1e-5: best = x0mm
        else: break
    print(f"{name:26s} 1 mm recovered: {ok!s:5s} | step error (%): {err:>5s} | max offset (mm): {best}")

evaluate("no mismatch", A0, B0)
for name, a21, bgain, a23 in [("B gain -30 %", 1, 0.7, 1), ("a21 +30 %", 1.3, 1, 1),
                              ("a23 -30 %", 1, 1, 0.7), ("a23 -50 %", 1, 1, 0.5),
                              ("a23 -60 %", 1, 1, 0.4), ("B gain -60 %", 1, 0.4, 1)]:
    A = A0.copy(); B = B0.copy()
    A[1, 0] *= a21; B[2, 0] *= bgain; A[1, 2] *= a23
    evaluate(name, A, B)
