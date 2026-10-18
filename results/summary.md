Open loop: 1 mm offset grows to 468 mm in 0.5 s

## Regulation from a 1 mm offset (settling within 2 % of start)
| Voltage limit | Controller | Settling (s) | Overshoot (%) |
|---|---|---|---|
| 1.0 V | PID | 1.402 | 19.9 |
| 1.0 V | MPC | 0.077 | 0.0 |
| 0.3 V | PID | 1.402 | 19.9 |
| 0.3 V | MPC | 0.115 | 0.0 |
| 0.2 V | PID | 1.402 | 19.9 |
| 0.2 V | MPC | 0.134 | 0.0 |

## 1 mm step reference
| Voltage limit | Controller | Rise (s) | Settling (s) | Overshoot (%) | Max voltage (V) |
|---|---|---|---|---|---|
| 1.0 V | PID | 0.070 | 1.478 | 61.3 | 0.15 |
| 1.0 V | MPC | 0.045 | 0.077 | 0.0 | 1.00 |
| 0.3 V | PID | 0.070 | 1.478 | 61.3 | 0.15 |
| 0.3 V | MPC | 0.058 | 0.109 | 0.0 | 0.30 |

## Largest starting offset recovered (0.5 mm resolution)
| Voltage limit | PID (mm) | MPC (mm) |
|---|---|---|
| 0.5 V | 5.0 | 5.5 |
| 0.3 V | 3.0 | 3.5 |
| 0.2 V | 2.0 | 2.0 |
