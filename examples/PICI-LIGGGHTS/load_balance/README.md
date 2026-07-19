# Dynamic particle-count balancing

Run the example on four MPI ranks:

```sh
mpirun -np 4 liggghts -in in.load_balance
```

The particles are deliberately created in the left fifth of the box, so the
initial Cartesian decomposition is strongly imbalanced. `fix balance` shifts
the x cutting planes and the thermo output reports its post-balance imbalance,
maximum rank population, iteration count, and pre-balance imbalance.

The `custom/vtk` dump writes each particle's current zero-based owning rank as
the integer `proc` field. To select one owner only, enable the commented
`dump_modify ... thresh proc` line.
