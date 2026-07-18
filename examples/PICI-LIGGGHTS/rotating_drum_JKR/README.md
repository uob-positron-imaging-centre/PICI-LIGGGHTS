# JKR rotating drum

This example demonstrates adhesive particles in the GranuTools-style rotating
drum. PICI-LIGGGHTS provides two compatible JKR normal-model choices.

## Legacy lookup-table model

The existing syntax remains the default compatibility path:

```text
fix resolution all property/global resolutionJKR scalar 1e-4
pair_style gran model jkr tangential history
fix walls all wall/gran model jkr tangential history mesh ...
```

It constructs and reports the original lookup table at startup. Lower
`resolutionJKR` values give finer interpolation but cover a smaller normalized
overlap range with the fixed table size.

## Analytical model

The opt-in analytical implementation removes the lookup table:

```text
pair_style gran model jkr/analytical tangential history
fix walls all wall/gran model jkr/analytical tangential history mesh ...
```

Both variants require `workOfAdhesion`; only legacy `jkr` requires
`resolutionJKR`. All other mechanics and material properties are shared.

`jkr_test.sim` continues to use legacy `jkr` so existing cases do not change.
The reproducible force validation, common-checkpoint drum comparison, plots,
and benchmarks for both styles are in `../jkr_comparison/`:

```sh
python3 validate_force_curve.py
python3 compare_and_benchmark.py
```

The comparison uses the existing drum meshes, settles one common nonadhesive
600-particle state, then runs each model for half a revolution from the same
particle IDs, positions, and velocities.

## Sources

- Eidevåg, Abrahamsson, Eng, and Rasmuson, *Modeling of dry snow adhesion
  during normal impact with surfaces*, Powder Technology 361 (2020), 1081–1092.
- Marshall, *Discrete-element modeling of particulate aerosol flows*, Journal
  of Computational Physics 228 (2009), 1541–1561.
