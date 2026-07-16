# Contact heat-generation example

The in.pp_collision input runs an isolated, damped two-particle collision.
Normal damping contributes to P_diss; the central contact hook converts it to
contactHeatSource, splits it equally, and the existing fix heat/gran transport
equation advances Temp.

Run the example from this directory:

    ../../../src/lmp_auto < in.pp_collision

Adjust the executable path for the selected build. energy.dat contains
mechanical energy, total particle thermal energy, and instantaneous generated
contact power at every timestep. final.dump exposes Temp and contactHeatSource
per particle.

## Automated validation

validate.py uses only the Python standard library. Edit the constants at the
top once for the local executable and MPI launcher, then keep the invocation
stable:

    python3 validate.py

The suite covers:

- disabled behavior against an unconfigured baseline;
- zero efficiency and an elastic collision;
- damped normal and sliding-friction particle-particle contacts;
- transient power accumulation with transport integration every five steps;
- particle-wall partitions of zero and one;
- one, two, and four MPI ranks for pair and wall contacts;
- coexistence with a persistent external heatSource;
- restart continuation;
- invalid, non-finite efficiency/partition/capacity values and missing thermal
  prerequisites.

The generated validation-results directory retains each exact input and screen
log but is ignored by Git.

## Reference energy balance

The following values were produced by the supplied suite with dt = 1e-6 s and
an MPICH build. Generated thermal energy is the timestep sum of
totalContactHeatPower multiplied by dt.

| Case | Initial mechanical | Final mechanical | Integrated generated power | Particle thermal gain | Difference | Relative error |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| damped particle-particle | 3.2724923474893676e-4 | 8.180353572234938e-5 | 2.4544548826667292e-4 | 2.4544548809046773e-4 | -1.7620518759364767e-13 | 7.179e-10 |
| sliding particle-particle | 8.181230868723419e-5 | 5.9262333796638004e-5 | 3.288824659223288e-5 | 3.2888246522588815e-5 | -6.964406536278528e-14 | 2.118e-9 |
| full particle-wall partition | 1.8980455615438333e-4 | 5.4380601176270303e-5 | 1.3245866331388048e-4 | 1.3245866336575318e-4 | 5.1872704265668035e-14 | 3.916e-10 |

The wall mechanical-energy difference also includes final rotational energy
and small integration error; the accounting criterion compares integrated
contact power directly with particle thermal-energy gain. Final pair and wall
state rows were identical across one, two, and four ranks in this validation.
