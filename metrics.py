"""Performance measures and the recoverable-offset search."""
import numpy as np
from model import TS


def step_metrics(ys, ref):
    """Rise time (10-90 %), settling time (2 %), overshoot (%). None if no valid response."""
    y = ys / ref
    t = np.arange(len(y)) * TS
    if not (y >= 0.9).any():
        return None
    rise = t[np.argmax(y >= 0.9)] - t[np.argmax(y >= 0.1)]
    try:
        settle = next(i for i in range(len(y)) if np.all(np.abs(y[i:] - 1) < 0.02)) * TS
    except StopIteration:
        return None
    return rise, settle, max(0.0, y.max() - 1) * 100


def recovered(ys, steps):
    return len(ys) == steps and abs(ys[-1]) < 1e-5


def max_recoverable_offset(make_controller, steps=300, step_mm=0.5, limit_mm=12.0):
    """Largest starting offset (mm, in 0.5 mm steps) the controller brings back to zero."""
    from controllers import simulate
    best = 0.0
    for x0_mm in np.arange(step_mm, limit_mm + 1e-9, step_mm):
        ys, _ = simulate(make_controller(), x0=(x0_mm / 1000, 0, 0), steps=steps)
        if recovered(ys, steps):
            best = x0_mm
        else:
            break
    return best
