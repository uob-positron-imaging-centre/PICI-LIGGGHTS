# Granular contact plus lubrication

`in.granular_plus_lubrication` demonstrates
`pair_style hybrid/overlay` with a Hertzian granular contact model and the
unequal-sphere `lubricate/bmpoly` resistance model.

The two spheres have radii 0.5 mm and 0.7 mm.  The initial centre separation
is 5 micrometres smaller than the sum of their radii, so the Hertz sub-style
has a small mechanical overlap while the lubrication sub-style evaluates its
regularised inner gap.  This is intentional for a compact overlay smoke test;
production cases should choose the mechanical contact stiffness and timestep
for their own material model.

The example uses `flagfld=0`, so only pairwise lubrication is added.  Set it to
1 when background-fluid drag, rotational drag, and stresslet contributions are
also desired.

Run it with a built executable from this directory:

```text
/path/to/liggghts -in in.granular_plus_lubrication
```

Both sub-styles are assigned to the same type pair:

```text
pair_style hybrid/overlay gran model hertz tangential history lubricate/bmpoly ...
pair_coeff * * gran
pair_coeff * * lubricate/bmpoly
```
