#!/usr/bin/env python3
"""Reproducible contact heat-generation validation suite.

Edit the constants below if the executable is in a different build directory.
The script needs only the Python standard library.
"""

from pathlib import Path
import os
import shutil
import subprocess


ROOT = Path(__file__).resolve().parents[3]
LIGGGHTS_BINARY = ROOT / "build" / "liggghts"
RESULTS_DIRECTORY = Path(__file__).resolve().parent / "validation-results"
RUN_MPI = True
MPI_RANKS = (1, 2, 4)
MPI_COMMAND = "mpiexec"
MPI_ENVIRONMENT = {}

DT = 1.0e-6
PP_STEPS = 6000
WALL_STEPS = 5000
ENERGY_TOLERANCE = 5.0e-9
MPI_TOLERANCE = 2.0e-12


PP_TEMPLATE = r"""
atom_style granular
atom_modify map array
boundary f f f
newton off
communicate single vel yes
units si
region domain block -0.02 0.02 -0.01 0.01 -0.01 0.01 units box
create_box 1 domain
create_atoms 1 single {x1} 0.0 0.0 units box
create_atoms 1 single {x2} 0.0 0.0 units box
set atom 1 diameter 0.01 density 2500.0 vx {vx1} vy {vy1}
set atom 2 diameter 0.01 density 2500.0 vx {vx2} vy {vy2}
neighbor 0.002 bin
neigh_modify delay 0
fix young all property/global youngsModulus peratomtype 5.0e6
fix poisson all property/global poissonsRatio peratomtype 0.30
fix restitution all property/global coefficientRestitution peratomtypepair 1 {restitution}
fix friction all property/global coefficientFriction peratomtypepair 1 {friction}
pair_style gran model hertz tangential history {model_settings}
pair_coeff * *
fix conductivity all property/global thermalConductivity peratomtype 0.1
fix capacity all property/global thermalCapacity peratomtype 100.0
fix thermal all heat/gran initial_temperature 300.0
fix integrate all nve/sphere
{external_source}
timestep 1.0e-6
compute kinetic all ke
compute temp_min all reduce min f_Temp
compute temp_max all reduce max f_Temp
compute source_sum all reduce sum f_heatSource
thermo_style custom step c_kinetic f_thermal c_temp_min c_temp_max c_source_sum {heat_columns}
thermo 1
thermo_modify lost error norm no format float %.17g
run {steps}
"""


WALL_TEMPLATE = r"""
atom_style granular
atom_modify map array
boundary f f f
newton off
communicate single vel yes
units si
region domain block -0.02 0.02 -0.01 0.01 -0.01 0.01 units box
create_box 1 domain
create_atoms 1 single -0.004 0.0 0.0 units box
set atom 1 diameter 0.01 density 2500.0 vx -0.5 vy 0.2
neighbor 0.002 bin
neigh_modify delay 0
fix young all property/global youngsModulus peratomtype 5.0e6
fix poisson all property/global poissonsRatio peratomtype 0.30
fix restitution all property/global coefficientRestitution peratomtypepair 1 0.50
fix friction all property/global coefficientFriction peratomtypepair 1 0.40
pair_style gran model hertz tangential history
pair_coeff * *
fix wall all wall/gran model hertz tangential history primitive type 1 xplane -0.01 heating_normal_hertz on heating_tangential_history on heat_generation on heat_generation_efficiency 1.0 heat_partition_wall {wall_partition} heat_generation_tracking on
fix conductivity all property/global thermalConductivity peratomtype 0.1
fix capacity all property/global thermalCapacity peratomtype 100.0
fix thermal all heat/gran initial_temperature 300.0
fix integrate all nve/sphere
timestep 1.0e-6
compute kinetic all ke
compute rotational all erotate/sphere
compute temp_min all reduce min f_Temp
compute temp_max all reduce max f_Temp
compute source_sum all reduce sum f_heatSource
thermo_style custom step c_kinetic c_rotational f_thermal c_temp_min c_temp_max c_source_sum f_totalContactHeatPower
thermo 1
thermo_modify lost error norm no format float %.17g
run {steps}
"""


def run_input(name, text, ranks=1, expect_success=True):
    case_directory = RESULTS_DIRECTORY / name
    case_directory.mkdir(parents=True, exist_ok=True)
    input_file = case_directory / "in.test"
    input_file.write_text(text)
    command = [str(LIGGGHTS_BINARY)]
    if ranks > 1:
        command = [MPI_COMMAND, "-n", str(ranks)] + command
    environment = os.environ.copy()
    environment.update(MPI_ENVIRONMENT)
    result = subprocess.run(
        command,
        input=text,
        text=True,
        cwd=case_directory,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        env=environment,
    )
    (case_directory / "screen.txt").write_text(result.stdout)
    if expect_success and result.returncode != 0:
        raise RuntimeError("{} failed; see {}".format(name, case_directory / "screen.txt"))
    if not expect_success and result.returncode == 0:
        raise RuntimeError("{} unexpectedly succeeded".format(name))
    return result.stdout, case_directory


def parse_thermo(output, columns):
    rows = []
    collecting = False
    for line in output.splitlines():
        fields = line.split()
        if fields and fields[0] == "Step":
            collecting = True
            continue
        if not collecting:
            continue
        if len(fields) != columns:
            if rows:
                break
            continue
        try:
            row = [float(value) for value in fields]
        except ValueError:
            if rows:
                break
            continue
        rows.append(row)
    if not rows:
        raise RuntimeError("No thermodynamic rows found")
    return rows


def pp_input(model_settings, restitution=0.5, friction=0.0,
             x1=-0.006, x2=0.006, vx1=0.5, vx2=-0.5,
             vy1=0.0, vy2=0.0, external_source="", heat=True,
             steps=PP_STEPS):
    heat_columns = "f_totalContactHeatPower" if heat else ""
    return PP_TEMPLATE.format(
        x1=x1, x2=x2, vx1=vx1, vx2=vx2, vy1=vy1, vy2=vy2,
        restitution=restitution, friction=friction,
        model_settings=model_settings, external_source=external_source,
        heat_columns=heat_columns, steps=steps,
    )


def assert_close(label, actual, expected, tolerance):
    scale = max(1.0, abs(actual), abs(expected))
    difference = abs(actual - expected)
    if difference > tolerance * scale:
        raise AssertionError(
            "{}: {} != {} (difference {})".format(
                label, actual, expected, difference
            )
        )


def energy_result(rows, kinetic_column, thermal_column, power_column):
    initial_mechanical = rows[0][kinetic_column]
    final_mechanical = rows[-1][kinetic_column]
    initial_thermal = rows[0][thermal_column]
    final_thermal = rows[-1][thermal_column]
    integrated_power = sum(row[power_column] for row in rows[1:]) * DT
    generated_thermal = final_thermal - initial_thermal
    difference = generated_thermal - integrated_power
    relative_error = abs(difference) / max(abs(integrated_power), 1.0e-30)
    return {
        "initial_mechanical": initial_mechanical,
        "final_mechanical": final_mechanical,
        "integrated_power": integrated_power,
        "generated_thermal": generated_thermal,
        "difference": difference,
        "relative_error": relative_error,
    }


def write_restart_reader(restart_file):
    return r"""
clear
read_restart "{restart}"
newton off
communicate single vel yes
neighbor 0.002 bin
neigh_modify delay 0
fix young all property/global youngsModulus peratomtype 5.0e6
fix poisson all property/global poissonsRatio peratomtype 0.30
fix restitution all property/global coefficientRestitution peratomtypepair 1 0.50
fix friction all property/global coefficientFriction peratomtypepair 1 0.0
pair_style gran model hertz tangential history heating_normal_hertz on heat_generation on heat_generation_efficiency 1.0 heat_partition_pp 0.5 heat_generation_tracking on
pair_coeff * *
fix conductivity all property/global thermalConductivity peratomtype 0.1
fix capacity all property/global thermalCapacity peratomtype 100.0
fix thermal all heat/gran initial_temperature 300.0
fix integrate all nve/sphere
timestep 1.0e-6
compute kinetic all ke
compute temp_min all reduce min f_Temp
compute temp_max all reduce max f_Temp
compute source_sum all reduce sum f_heatSource
thermo_style custom step c_kinetic f_thermal c_temp_min c_temp_max c_source_sum f_totalContactHeatPower
thermo 1
thermo_modify lost error norm no format float %.17g
run 3200
""".format(restart=restart_file)


def main():
    if not LIGGGHTS_BINARY.is_file():
        raise SystemExit("Set LIGGGHTS_BINARY at the top of validate.py")
    if RESULTS_DIRECTORY.exists():
        shutil.rmtree(RESULTS_DIRECTORY)
    RESULTS_DIRECTORY.mkdir()

    baseline_text = pp_input("heating_normal_hertz on", heat=False)
    disabled_text = pp_input("heating_normal_hertz on heat_generation off", heat=False)
    baseline = parse_thermo(run_input("A_baseline", baseline_text)[0], 6)
    disabled = parse_thermo(run_input("A_disabled", disabled_text)[0], 6)
    for column in (1, 2, 3, 4, 5):
        assert_close("disabled output", disabled[-1][column], baseline[-1][column], 0.0)

    zero_settings = "heating_normal_hertz on heat_generation on heat_generation_efficiency 0.0 heat_generation_tracking on"
    zero = parse_thermo(run_input("B_zero_efficiency", pp_input(zero_settings))[0], 7)
    assert_close("zero-efficiency temperature", zero[-1][2], zero[0][2], ENERGY_TOLERANCE)
    assert_close("zero-efficiency power", max(row[6] for row in zero), 0.0, 0.0)

    elastic_settings = "heating_normal_hertz on heat_generation on heat_generation_tracking on"
    elastic = parse_thermo(
        run_input("C_elastic", pp_input(elastic_settings, restitution=1.0))[0], 7
    )
    assert_close("elastic temperature", elastic[-1][2], elastic[0][2], ENERGY_TOLERANCE)
    assert_close("elastic power", max(row[6] for row in elastic), 0.0, 1.0e-15)

    damped_settings = "heating_normal_hertz on heat_generation on heat_generation_efficiency 1.0 heat_partition_pp 0.5 heat_generation_tracking on"
    damped_text = pp_input(damped_settings)
    damped = parse_thermo(run_input("D_damped_pp", damped_text)[0], 7)
    damped_energy = energy_result(damped, 1, 2, 6)
    if damped_energy["relative_error"] > ENERGY_TOLERANCE:
        raise AssertionError("particle-particle energy balance failed")
    assert_close("equal final particle temperatures", damped[-1][3], damped[-1][4], ENERGY_TOLERANCE)

    nevery_text = damped_text.replace(
        "fix thermal all heat/gran initial_temperature 300.0\n",
        "fix thermal all heat/gran initial_temperature 300.0\n"
        "fix_modify ste_heattransfer every 5\n",
    )
    nevery = parse_thermo(run_input("D_nevery", nevery_text)[0], 7)
    nevery_energy = energy_result(nevery, 1, 2, 6)
    if nevery_energy["relative_error"] > ENERGY_TOLERANCE:
        raise AssertionError("nevery contact-power accumulation failed")

    sliding_settings = "heating_normal_hertz off heating_tangential_history on heat_generation on heat_generation_efficiency 1.0 heat_partition_pp 0.5 heat_generation_tracking on"
    sliding_text = pp_input(
        sliding_settings, friction=0.5, x1=-0.0049, x2=0.0049,
        vx1=0.0, vx2=0.0, vy1=0.25, vy2=-0.25, steps=2000,
    )
    sliding = parse_thermo(run_input("E_sliding_pp", sliding_text)[0], 7)
    sliding_energy = energy_result(sliding, 1, 2, 6)
    if sliding_energy["generated_thermal"] <= 0.0:
        raise AssertionError("sliding contact did not generate heat")
    if sliding_energy["relative_error"] > ENERGY_TOLERANCE:
        raise AssertionError("sliding energy balance failed")

    wall_zero_text = WALL_TEMPLATE.format(wall_partition=0.0, steps=WALL_STEPS)
    wall_full_text = WALL_TEMPLATE.format(wall_partition=1.0, steps=WALL_STEPS)
    wall_zero = parse_thermo(run_input("F_wall_zero", wall_zero_text)[0], 8)
    wall_full = parse_thermo(run_input("F_wall_full", wall_full_text)[0], 8)
    assert_close("wall zero partition", wall_zero[-1][3], wall_zero[0][3], ENERGY_TOLERANCE)
    if wall_full[-1][3] <= wall_full[0][3]:
        raise AssertionError("wall full partition did not heat particle")
    assert_close("wall mechanical partition independence", wall_zero[-1][1], wall_full[-1][1], MPI_TOLERANCE)
    assert_close("wall rotational partition independence", wall_zero[-1][2], wall_full[-1][2], MPI_TOLERANCE)
    wall_energy = energy_result(wall_full, 1, 3, 7)
    if wall_energy["relative_error"] > ENERGY_TOLERANCE:
        raise AssertionError("particle-wall energy balance failed")

    source_text = pp_input(
        damped_settings,
        external_source="run 0\nset group all property/atom heatSource 0.25",
    )
    source = parse_thermo(run_input("H_persistent_source", source_text)[0], 7)
    assert_close("persistent heat source", source[-1][5], 0.5, 0.0)
    source_contact_energy = sum(row[6] for row in source[1:]) * DT
    expected_source_energy = 0.5 * PP_STEPS * DT
    actual_source_energy = source[-1][2] - source[0][2]
    assert_close(
        "combined external and contact heat",
        actual_source_energy,
        source_contact_energy + expected_source_energy,
        ENERGY_TOLERANCE,
    )

    writer_text = pp_input(damped_settings, steps=2800) + "\nwrite_restart restart.bin\n"
    _, writer_directory = run_input("I_restart_write", writer_text)
    restart_file = writer_directory / "restart.bin"
    restart_output, _ = run_input(
        "I_restart_read", write_restart_reader(restart_file)
    )
    restarted = parse_thermo(restart_output, 7)
    for column in (1, 2, 3, 4, 5):
        assert_close(
            "restart final output",
            restarted[-1][column], damped[-1][column], ENERGY_TOLERANCE,
        )

    invalid_settings = (
        ("invalid_efficiency_low", "heat_generation_efficiency -0.1"),
        ("invalid_efficiency_high", "heat_generation_efficiency 1.1"),
        ("invalid_partition_pp_low", "heat_partition_pp -0.1"),
        ("invalid_partition_pp_high", "heat_partition_pp 1.1"),
        ("invalid_partition_wall_low", "heat_partition_wall -0.1"),
        ("invalid_partition_wall_high", "heat_partition_wall 1.1"),
        ("invalid_efficiency_nan", "heat_generation_efficiency nan"),
        ("invalid_partition_pp_nan", "heat_partition_pp nan"),
        ("invalid_partition_wall_nan", "heat_partition_wall nan"),
    )
    for name, setting in invalid_settings:
        text = pp_input(
            "heating_normal_hertz on heat_generation on " + setting
        )
        run_input(name, text, expect_success=False)

    missing_thermal = damped_text.replace(
        "fix thermal all heat/gran initial_temperature 300.0\n", ""
    )
    missing_thermal = missing_thermal.replace(
        "compute temp_min all reduce min f_Temp\n", ""
    )
    missing_thermal = missing_thermal.replace(
        "compute temp_max all reduce max f_Temp\n", ""
    )
    missing_thermal = missing_thermal.replace(
        "compute source_sum all reduce sum f_heatSource\n", ""
    )
    missing_thermal = missing_thermal.replace(
        "thermo_style custom step c_kinetic f_thermal c_temp_min c_temp_max "
        "c_source_sum f_totalContactHeatPower",
        "thermo_style custom step c_kinetic f_totalContactHeatPower",
    )
    run_input("invalid_missing_thermal", missing_thermal, expect_success=False)

    missing_capacity = damped_text.replace(
        "fix capacity all property/global thermalCapacity peratomtype 100.0\n", ""
    )
    run_input("invalid_missing_capacity", missing_capacity, expect_success=False)
    for name, value in (
        ("invalid_zero_capacity", "0.0"),
        ("invalid_negative_capacity", "-1.0"),
        ("invalid_nan_capacity", "nan"),
    ):
        text = damped_text.replace(
            "thermalCapacity peratomtype 100.0",
            "thermalCapacity peratomtype " + value,
        )
        run_input(name, text, expect_success=False)

    mpi_results = []
    if RUN_MPI:
        for ranks in MPI_RANKS:
            pp_rows = parse_thermo(
                run_input("G_pp_{}r".format(ranks), damped_text, ranks)[0], 7
            )
            wall_rows = parse_thermo(
                run_input("G_wall_{}r".format(ranks), wall_full_text, ranks)[0], 8
            )
            mpi_results.append((ranks, pp_rows[-1], wall_rows[-1]))
        reference = mpi_results[0]
        for ranks, pp_final, wall_final in mpi_results[1:]:
            for column in (1, 2, 3, 4, 5):
                assert_close(
                    "MPI particle-particle {} ranks".format(ranks),
                    pp_final[column], reference[1][column], MPI_TOLERANCE,
                )
            for column in (1, 2, 3, 4, 5, 6):
                assert_close(
                    "MPI particle-wall {} ranks".format(ranks),
                    wall_final[column], reference[2][column], MPI_TOLERANCE,
                )

    print("Contact heat-generation validation passed")
    print("case initial_mechanical final_mechanical integrated_dissipation generated_thermal difference relative_error")
    for name, result in (
        ("damped_pp", damped_energy),
        ("sliding_pp", sliding_energy),
        ("wall_full", wall_energy),
    ):
        print(
            "{} {:.17g} {:.17g} {:.17g} {:.17g} {:.17g} {:.3e}".format(
                name,
                result["initial_mechanical"],
                result["final_mechanical"],
                result["integrated_power"],
                result["generated_thermal"],
                result["difference"],
                result["relative_error"],
            )
        )
    if not RUN_MPI:
        print("MPI cases skipped; set RUN_MPI = True to enable them")


if __name__ == "__main__":
    main()
