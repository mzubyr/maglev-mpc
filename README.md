# Maglev MPC: PID vs Model Predictive Control in Python

A Python re-implementation and extension of my B.Eng thesis
(*Design of a Model Predictive Controller for a Magnetic Suspension System*,
originally done in MATLAB/Simulink). A steel ball is held in the air by an
electromagnet. The plant is nonlinear and open-loop unstable, so it needs feedback.

**Everything here is simulation. No hardware was used.**

## What this project does

- Linearises the nonlinear maglev model around an equilibrium and checks the
  result numerically (`physics.py`).
- Discretises the model (Ts = 0.0064 s) and shows the open loop is unstable.
- Implements a **PID** controller (gains from the thesis) and a **linear MPC**
  with input (voltage) constraints.
- Compares them on regulation, step tracking and **how large a starting offset
  each can recover** under a voltage limit.
- Tests how the MPC behaves when the real plant differs from the model, and
  whether a **data-driven model update** (regularised least squares on noisy
  closed-loop data) recovers step-tracking accuracy.

## Model

States: position x1, velocity x2, coil current x3. Input: coil voltage u.
Linearised around x1 = 0.05 m (equilibrium current and voltage about 6.264 A / 6.264 V):

```
A = [[0, 1, 0], [196.2, -0.01, -3.1312], [0, 12.5284, -40]]
B = [0, 0, 40]^T      C = [1, 0, 0]
```

Eigenvalues of A are about +13.6, -14.8 and -38.8. The positive one makes the
open loop unstable (a 1 mm offset grows to roughly 470 mm in 0.5 s).

`python physics.py` compares A and B with a numerical linearisation of the
nonlinear model. The maximum difference is about 0.0009 (one entry, from rounding).

## Controllers

- **PID**: P = -152, I = -190, D = -16 (from the thesis), with output
  saturation and conditional-integration anti-windup.
- **MPC**: prediction horizon 29, control horizon 4, cost on position/velocity
  errors and on voltage, bounded input. Solved at every step as a bounded
  least-squares problem (`scipy.optimize.lsq_linear`); only the first input is applied.

## Results (reproduce with `python run_experiments.py`)

Regulation from a 1 mm offset (settling within 2 % of the start):

| Voltage limit | Controller | Settling (s) | Overshoot (%) |
|---|---|---|---|
| 1.0 V | PID | 1.402 | 19.9 |
| 1.0 V | MPC | 0.077 | 0.0 |
| 0.3 V | PID | 1.402 | 19.9 |
| 0.3 V | MPC | 0.115 | 0.0 |
| 0.2 V | PID | 1.402 | 19.9 |
| 0.2 V | MPC | 0.134 | 0.0 |

1 mm step reference:

| Voltage limit | Controller | Rise (s) | Settling (s) | Overshoot (%) | Max voltage (V) |
|---|---|---|---|---|---|
| 1.0 V | PID | 0.070 | 1.478 | 61.3 | 0.15 |
| 1.0 V | MPC | 0.045 | 0.077 | 0.0 | 1.00 |
| 0.3 V | PID | 0.070 | 1.478 | 61.3 | 0.15 |
| 0.3 V | MPC | 0.058 | 0.109 | 0.0 | 0.30 |

Largest starting offset recovered (0.5 mm resolution):

| Voltage limit | PID (mm) | MPC (mm) |
|---|---|---|
| 0.5 V | 5.0 | 5.5 |
| 0.3 V | 3.0 | 3.5 |
| 0.2 V | 2.0 | 2.0 |

Figures: `results/open_loop.png`, `results/step_response.png`.

## Model mismatch and data-driven update

`python mismatch_study.py` runs the MPC (built from the nominal model) on plants
with perturbed parameters, with a 0.3 V limit:

| True plant | Step error | Max recoverable offset |
|---|---|---|
| No mismatch | 0.0 % | 3.5 mm |
| B gain -30 % | -0.5 % | 2.0 mm |
| a21 +30 % | -1.7 % | 2.5 mm |
| a23 -30 % | -2.4 % | 2.0 mm |
| a23 -50 % | -5.8 % | 1.5 mm |
| a23 -60 % | step response lost | 1.0 mm |
| B gain -60 % | step response lost | 1.0 mm |

Feedback tolerates moderate model errors, but the steady-state error grows
(the MPC has no integral action and takes its steady-state voltage from the
nominal model) and the recoverable range shrinks.

`python data_driven_update.py` collects 3000 samples of closed-loop data from
the mismatched plant (nominal MPC plus 0.1 V random dither, simulated sensor
noise), fits a corrected discrete model by ridge regression pulled toward the
nominal model, and rebuilds the MPC. Five random seeds each:

| True plant | Nominal MPC | Learned MPC (5 runs) | Oracle (true model) |
|---|---|---|---|
| a23 -30 %, B -20 % | -2.80 % | -0.19 % to -0.57 % | 0.0 % |
| a23 -50 %, B -30 % | step response lost | -0.40 % to +0.07 % | 0.0 % |

The recoverable offset did **not** improve (1.5 mm and 1.0 mm in both cases). It is
limited by the voltage limit, not by the model. Plain (unregularised) least squares
gave a model that was *worse* than the nominal one with noisy closed-loop data,
which is why the update is pulled toward the nominal model.

### What the results show

- MPC settles much faster and without overshoot, but part of that gap comes
  from my MPC weights and the fact that I did not re-tune the PID.
- Rise time and overshoot depend strongly on the MPC voltage weight `r`
  (larger `r` gives a slower, gentler response, and a very large `r` fails to stabilise).
- With a bounded input, an unstable plant can only be recovered from a limited
  set of starting states. Lowering the voltage limit shrinks that set for *both*
  controllers, and MPC's advantage there is small (0 to 0.5 mm).

## Run

```
pip install -r requirements.txt
python physics.py
python run_experiments.py        # PID vs MPC, figures and tables in results/
python mismatch_study.py         # nominal MPC on mismatched plants
python data_driven_update.py     # ridge-regression model update from noisy data
python -m pytest tests        # optional
```

## Limitations

- Linear model valid only near the equilibrium; the ball is considered lost at 50 mm.
- No measurement noise, disturbances or model mismatch.
- The MPC penalises the voltage itself rather than its change, which is a
  simplification compared with the MATLAB toolbox used in the thesis.
- The physical parameters in `physics.py` were reverse-engineered from the
  linearised matrices. TODO: confirm against the thesis parameter table.
- PID gains are taken from the thesis and were not re-tuned.
- The data-driven update assumes the full state is measured. A real maglev
  usually measures position only and would need an observer.
- Sensor noise and mismatch are simulated, and only two parameters were perturbed.
- A simple integral term in the MPC would also remove the steady-state error
  without any learning.

## Next steps

- Constrained MPC for a 2D quadrotor (position and no-go-zone constraints).
- Step tracking with measurement noise and model mismatch.

## Author

Muhammad Ibrahim Zubair, Kano, Nigeria.
