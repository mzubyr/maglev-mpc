"""Runs every experiment and writes figures + a summary to results/."""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from model import AD, TS
from controllers import MPC, PID, simulate
from metrics import step_metrics, max_recoverable_offset

OUT = "results"
os.makedirs(OUT, exist_ok=True)
lines = []


def log(s=""):
    print(s); lines.append(s)


# 1. Open loop is unstable
x = np.array([0.001, 0, 0]); pos = [x[0]]
for _ in range(78):
    x = AD @ x; pos.append(x[0])
plt.figure(); plt.plot(np.arange(79) * TS, np.array(pos) * 1000)
plt.xlabel("time (s)"); plt.ylabel("position deviation (mm)")
plt.title("Open loop: 1 mm offset grows"); plt.grid(True)
plt.savefig(f"{OUT}/open_loop.png", dpi=150); plt.close()
log(f"Open loop: 1 mm offset grows to {pos[-1]*1000:.0f} mm in 0.5 s")

# 2. Regulation from a 1 mm offset
log("\n## Regulation from a 1 mm offset (settling within 2 % of start)")
log("| Voltage limit | Controller | Settling (s) | Overshoot (%) |")
log("|---|---|---|---|")
for umax in [1.0, 0.3, 0.2]:
    for name, make in [("PID", lambda: PID(umax=umax)), ("MPC", lambda: MPC(umax=umax))]:
        ys, _ = simulate(make(), x0=(0.001, 0, 0), steps=300)
        y = ys / 0.001
        ok = len(ys) == 300 and abs(ys[-1]) < 1e-5
        if ok:
            ts = next(i for i in range(len(y)) if np.all(np.abs(y[i:]) < 0.02)) * TS
            log(f"| {umax} V | {name} | {ts:.3f} | {max(0.0, -y.min())*100:.1f} |")
        else:
            log(f"| {umax} V | {name} | not recovered | - |")

# 3. Step reference
log("\n## 1 mm step reference")
log("| Voltage limit | Controller | Rise (s) | Settling (s) | Overshoot (%) | Max voltage (V) |")
log("|---|---|---|---|---|---|")
curves = {}
for umax in [1.0, 0.3]:
    for name, make in [("PID", lambda: PID(umax=umax)), ("MPC", lambda: MPC(umax=umax))]:
        ys, us = simulate(make(), ref=0.001, steps=400)
        m = step_metrics(ys, 0.001)
        log(f"| {umax} V | {name} | {m[0]:.3f} | {m[1]:.3f} | {m[2]:.1f} | {abs(us).max():.2f} |")
        if umax == 0.3:
            curves[name] = (ys, us)
t = np.arange(400) * TS
fig, ax = plt.subplots(2, 1, sharex=True)
for name, (ys, us) in curves.items():
    ax[0].plot(t, ys * 1000, label=name); ax[1].plot(t, us, label=name)
ax[0].axhline(1.0, color="gray", ls="--", label="reference")
ax[0].set_ylabel("position (mm)"); ax[0].legend(); ax[0].grid(True)
ax[1].set_ylabel("voltage (V)"); ax[1].set_xlabel("time (s)"); ax[1].grid(True)
fig.suptitle("1 mm step, voltage limit 0.3 V")
plt.savefig(f"{OUT}/step_response.png", dpi=150); plt.close()

# 4. Recoverable offset
log("\n## Largest starting offset recovered (0.5 mm resolution)")
log("| Voltage limit | PID (mm) | MPC (mm) |")
log("|---|---|---|")
for umax in [0.5, 0.3, 0.2]:
    p = max_recoverable_offset(lambda: PID(umax=umax))
    m = max_recoverable_offset(lambda: MPC(umax=umax))
    log(f"| {umax} V | {p} | {m} |")

open(f"{OUT}/summary.md", "w").write("\n".join(lines) + "\n")
