# Legacy/analytical JKR verification and benchmarks

This directory provides the reproducible checks used for the parallel JKR
implementations. Set `LIGGGHTS_BINARY` near the top of each Python script if
the executable is not in one of the listed build locations.

## Force-law validation

```sh
python3 validate_force_curve.py
```

The two-sphere case unloads through adhesive necking and snap-off, then reloads
through contact. Both JKR variants are compared with an independent bisection
solution of the quartic. The script runs on one and two MPI ranks, requires an
analytical normalized-force error below `1e-10`, requires rank agreement within
`1e-12`, and reports rather than suppresses legacy LUT interpolation error.

## Drum comparison and benchmarks

```sh
python3 compare_and_benchmark.py
```

The script first creates one nonadhesive, settled 600-particle checkpoint in
the existing rotating-drum mesh. Legacy and analytical JKR then start from that
same data file and run 100,000 steps at 30 rpm. Completed checkpoint and drum
trajectories are reused on a resumed run; benchmark repetitions are always
fresh.

Generated artifacts are ignored by Git:

- `validation-results/`: force errors, MPI comparison, and summary;
- `comparison-results/`: ID-aligned trajectory dumps, divergence CSV, SVG
  figures, final VTK states, and summary;
- `benchmark-results/`: fresh-process startup/RSS and frozen-contact timings.

The benchmark constants and repetition counts are intentionally kept at the
top of `compare_and_benchmark.py`, and model order alternates between measured
repetitions.
