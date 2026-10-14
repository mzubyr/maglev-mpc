"""Data-driven model update: collect noisy closed-loop data from a mismatched plant, fit a corrected
discrete model by ridge regression pulled toward the nominal model, and rebuild the MPC.
Self-contained; run with: python data_driven_update.py
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

Ad0, Bd0 = disc(A0, B0)

def steady_state(Ad, Bd, ref):
    """x_ss, u_ss that hold position = ref, for the model (Ad, Bd)."""
    K = np.block([[Ad - np.eye(3), Bd], [C, np.zeros((1, 1))]])
    sol = np.linalg.solve(K, np.concatenate([np.zeros(3), [ref]]))
    return sol[:3], sol[3]

def make_mpc(Ad, Bd, Np=29, Nc=4, q=(1e6, 1e2, 0.0), r=1e-2, umax=0.3):
    F = np.vstack([np.linalg.matrix_power(Ad, i) for i in range(1, Np + 1)])
    G = np.zeros((3 * Np, Nc))
    for i in range(1, Np + 1):
        for j in range(i):
            G[3*(i-1):3*i, min(j, Nc-1)] += (np.linalg.matrix_power(Ad, i-1-j) @ Bd).ravel()
    Qh = np.kron(np.eye(Np), np.diag(np.sqrt(q)))
    M = np.vstack([Qh @ G, np.sqrt(r) * np.eye(Nc)])
    def control(x, ref):
        xs, uss = steady_state(Ad, Bd, ref)
        rhs = np.concatenate([Qh @ (np.tile(xs, Np) - F @ x), np.sqrt(r) * uss * np.ones(Nc)])
        return lsq_linear(M, rhs, bounds=(-umax, umax), method="bvls").x[0]
    return control

def simulate(ctrl, Adt, Bdt, x0=(0, 0, 0), ref=0.0, steps=400):
    x = np.array(x0, float); ys = []
    for _ in range(steps):
        x = Adt @ x + Bdt.ravel() * ctrl(x, ref)
        ys.append(x[0])
        if abs(x[0]) > 0.05: break
    return np.array(ys)

def evaluate(ctrl, Adt, Bdt):
    ys = simulate(ctrl, Adt, Bdt, ref=0.001, steps=400)
    err = None if len(ys) < 400 else (0.001 - ys[-1]) / 0.001 * 100
    best = 0.0
    for x0mm in np.arange(0.5, 8.01, 0.5):
        yy = simulate(ctrl, Adt, Bdt, x0=(x0mm / 1000, 0, 0), steps=300)
        if len(yy) == 300 and abs(yy[-1]) < 1e-5: best = x0mm
        else: break
    return err, best

def collect(ctrl, Adt, Bdt, rng, n=3000, dither=0.1, noise=(2e-6, 2e-4, 1e-4), umax=0.3):
    """Closed-loop data from the 'real' plant: nominal MPC + random voltage dither + sensor noise."""
    x = np.array([0.0005, 0, 0]); X, U, Xn = [], [], []
    for _ in range(n):
        u = np.clip(ctrl(x, 0.0) + dither * rng.standard_normal(), -umax, umax)
        xn = Adt @ x + Bdt.ravel() * u
        X.append(x + noise * rng.standard_normal(3)); U.append(u)
        Xn.append(xn + noise * rng.standard_normal(3))
        x = xn
        if abs(x[0]) > 0.02: x = np.array([0.0005, 0, 0])
    return np.array(X), np.array(U), np.array(Xn)

def identify(X, U, Xn, lam):
    """Ridge regression pulled toward the nominal model: x+ = Ad x + Bd u."""
    s = np.array([1e-3, 5e-2, 5e-3, 0.1])                 # typical sizes of x1, x2, x3, u
    Z = np.hstack([X, U[:, None]]) / s
    Th0 = np.hstack([Ad0, Bd0]).T * s[:, None]
    Th = np.linalg.solve(Z.T @ Z + lam * np.eye(4), Z.T @ Xn + lam * Th0) / s[:, None]
    return Th[:3].T, Th[3:].T

def fmt(err): return "lost" if err is None else f"{err:+.2f}"

nominal = make_mpc(Ad0, Bd0)
for name, a23, bgain in [("a23 -30 %, B gain -20 % (mild)", 0.7, 0.8),
                         ("a23 -50 %, B gain -30 % (hard)", 0.5, 0.7)]:
    A = A0.copy(); B = B0.copy(); A[1, 2] *= a23; B[2, 0] *= bgain
    Adt, Bdt = disc(A, B)
    print("\nTrue plant:", name)
    e, b = evaluate(nominal, Adt, Bdt)
    print(f"  nominal MPC         step error {fmt(e):>6s} %   max offset {b} mm")
    for seed in range(5):
        rng = np.random.default_rng(seed)
        X, U, Xn = collect(nominal, Adt, Bdt, rng)
        Adh, Bdh = identify(X, U, Xn, lam=10)
        e, b = evaluate(make_mpc(Adh, Bdh), Adt, Bdt)
        print(f"  learned MPC, run {seed}  step error {fmt(e):>6s} %   max offset {b} mm")
    e, b = evaluate(make_mpc(Adt, Bdt), Adt, Bdt)
    print(f"  oracle (true model) step error {fmt(e):>6s} %   max offset {b} mm")
