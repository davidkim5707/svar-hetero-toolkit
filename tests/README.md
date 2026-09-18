# Tests

```matlab
run tests/run_tests.m
```

`run_tests.m` runs four self-tests of the zero-restriction routines in `code/core/zero/`. Each test builds a synthetic model, runs one routine, and prints its checks. A failing check raises an error. The four tests together take under a minute.

| Test | What it checks |
|---|---|
| `test_build_zero_constraint` | A column in the null space of the constraint sets the targeted impulse response to zero. Two restricted shocks give a feasible build order. An empty set returns no constraint. Duplicate restrictions collapse to one row. |
| `test_draw_Q_zero_columnwise` | The drawn rotation is orthonormal, satisfies the zeros and the impact-sign screen, and equals the rotation that `gaussian_to_Q_zero` builds from the stored Gaussian matrix. One and two restricted shocks. |
| `test_ess_draw_Q_zero_columnwise` | Along many elliptical slice sweeps every rotation stays orthonormal, satisfies the zeros, and passes the impact-sign screen, and the chain moves. One and two restricted shocks. |
| `test_log_volume_element_zero` | The log volume element is zero when zeros restrict a single shock. With two restricted shocks it equals the value computed directly with the routines of Arias, Rubio-Ramirez, and Waggoner (2018). Two feasible build orders give different volume elements. |

The tests cover the zero-restriction routines only. The examples exercise the two samplers end to end.
