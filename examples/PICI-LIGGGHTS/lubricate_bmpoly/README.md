# Unequal-sphere lubrication example

`in.lubricate_bmpoly` runs two unequal spheres with the `lubricate/bmpoly`
pair style from Cheal and Ness, arXiv:1801.03805v3.  It exercises normal and
tangential relative motion, particle rotation, background-fluid drag, and
the finite-volume correction.

Run it with a built executable, from this directory:

```text
/path/to/liggghts -in in.lubricate_bmpoly
```

The input uses a 1 micrometre inner gap and a 50 micrometre outer gap.  The
cutoffs are surface-to-surface gaps; the implementation converts them to
centre-distance neighbour cutoffs using the particle radii.

This pair style computes hydrodynamic lubrication forces and torques only. It
does not supply a repulsive mechanical contact law, so a production DEM case
that permits overlap should add and validate an appropriate contact model.

For the combined mechanical-contact and lubrication setup, see
`granular_plus_lubrication/in.granular_plus_lubrication`.
