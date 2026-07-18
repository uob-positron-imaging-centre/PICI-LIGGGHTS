#!/usr/bin/env python3
"""Validate legacy and analytical JKR force curves with no extra packages."""

from pathlib import Path
import csv
import math
import shutil
import subprocess


ROOT = Path(__file__).resolve().parents[3]
BINARY_CANDIDATES = (
    ROOT / "build" / "liggghts",
    ROOT / "src" / "lmp_auto",
    Path("/tmp/pici_liggghts_build_mpich3/liggghts"),
    Path("/tmp/pici_liggghts_build_cxx14/liggghts"),
)
LIGGGHTS_BINARY = next((path for path in BINARY_CANDIDATES if path.is_file()), BINARY_CANDIDATES[0])
MPI_COMMAND = "/opt/homebrew/opt/mpich/bin/mpiexec"
MPI_OPTIONS = ()
MPI_RANKS = (1, 2)
RESULTS_DIRECTORY = Path(__file__).resolve().parent / "validation-results"

YOUNGS_MODULUS = 5.0e6
POISSONS_RATIO = 0.30
WORK_OF_ADHESION = 0.10
PARTICLE_RADIUS = 0.005
ANALYTICAL_TOLERANCE = 1.0e-10
MPI_TOLERANCE = 1.0e-12


def run_case(model, ranks):
    case_directory = RESULTS_DIRECTORY / "{}_rank{}".format(model, ranks)
    case_directory.mkdir(parents=True, exist_ok=True)
    input_file = Path(__file__).resolve().parent / "in.force_curve_{}".format(model)
    command = [str(LIGGGHTS_BINARY), "-in", str(input_file)]
    if ranks > 1:
        command = [MPI_COMMAND] + list(MPI_OPTIONS) + ["-n", str(ranks)] + command
    result = subprocess.run(
        command,
        cwd=case_directory,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    (case_directory / "screen.txt").write_text(result.stdout)
    dump_file = case_directory / "force_curve.dump"
    macos_mpich_finalize_error = (
        result.returncode == 143 and
        "Fatal error in internal_Finalize" in result.stdout and
        result.stdout.count("Loop time of") == 2 and
        dump_file.is_file()
    )
    if result.returncode != 0 and not macos_mpich_finalize_error:
        raise RuntimeError("{} failed; see {}".format(model, case_directory / "screen.txt"))
    snapshots = parse_dump(dump_file)
    if max(snapshots) != 6000:
        raise RuntimeError("{} did not complete 6000 steps".format(model))
    return snapshots, macos_mpich_finalize_error


def parse_dump(path):
    snapshots = {}
    lines = path.read_text().splitlines()
    index = 0
    while index < len(lines):
        if lines[index] != "ITEM: TIMESTEP":
            raise RuntimeError("Unexpected dump format at line {}".format(index + 1))
        step = int(lines[index + 1])
        atom_count = int(lines[index + 3])
        header_index = index + 8
        columns = lines[header_index].split()[2:]
        atoms = {}
        for row_index in range(header_index + 1, header_index + 1 + atom_count):
            values = lines[row_index].split()
            atom = {name: float(value) for name, value in zip(columns, values)}
            atoms[int(atom["id"])] = atom
        snapshots[step] = atoms
        index = header_index + 1 + atom_count
    return snapshots


def reference_contact_ratio(normalized_overlap):
    if normalized_overlap < -1.0:
        return 0.0, 0.0
    cbrt_six = 6.0 ** (1.0 / 3.0)
    coefficient_a = 2.0 * cbrt_six
    coefficient_b = 4.0 * cbrt_six / 3.0
    lower = (1.0 / 6.0) ** (1.0 / 3.0)

    def polynomial(value):
        return coefficient_a * value ** 4 - coefficient_b * value - normalized_overlap

    if abs(polynomial(lower)) < 1.0e-15:
        root = lower
    else:
        upper = max(1.0, (max(0.0, normalized_overlap) / coefficient_a) ** 0.25 + 1.0)
        while polynomial(upper) <= 0.0:
            upper *= 2.0
        for unused in range(200):
            midpoint = 0.5 * (lower + upper)
            if polynomial(midpoint) <= 0.0:
                lower = midpoint
            else:
                upper = midpoint
        root = 0.5 * (lower + upper)

    a_a0 = root * root
    fn_fc = 4.0 * (a_a0 ** 3 - root ** 3)
    return a_a0, fn_fc


def physical_scales():
    effective_radius = PARTICLE_RADIUS / 2.0
    effective_youngs = 1.0 / (
        2.0 * (1.0 - POISSONS_RATIO ** 2) / YOUNGS_MODULUS
    )
    a0 = (
        4.5 * WORK_OF_ADHESION * math.pi * effective_radius ** 2 /
        effective_youngs
    ) ** (1.0 / 3.0)
    delta_c = a0 ** 2 / (2.0 * 6.0 ** (1.0 / 3.0) * effective_radius)
    fc = 1.5 * math.pi * WORK_OF_ADHESION * effective_radius
    return delta_c, fc


def expected_force(step, atoms):
    delta_c, fc = physical_scales()
    separation = atoms[2]["x"] - atoms[1]["x"]
    normalized_overlap = (2.0 * PARTICLE_RADIUS - separation) / delta_c
    active = (step <= 3000 and normalized_overlap > -1.0) or (
        step > 3000 and normalized_overlap > 0.0
    )
    if not active:
        return normalized_overlap, 0.0
    unused_a, fn_fc = reference_contact_ratio(normalized_overlap)
    return normalized_overlap, fn_fc * fc


def force_errors(snapshots):
    unused_delta_c, fc = physical_scales()
    rows = []
    maximum_normalized_error = 0.0
    for step in sorted(snapshots):
        atoms = snapshots[step]
        normalized_overlap, reference = expected_force(step, atoms)
        measured = atoms[2]["fx"]
        normalized_error = abs(measured - reference) / fc
        maximum_normalized_error = max(maximum_normalized_error, normalized_error)
        rows.append((step, normalized_overlap, measured, reference, normalized_error))
    return maximum_normalized_error, rows


def maximum_snapshot_difference(first, second):
    maximum = 0.0
    if set(first) != set(second):
        raise AssertionError("MPI runs produced different timesteps")
    for step in first:
        if set(first[step]) != set(second[step]):
            raise AssertionError("MPI runs produced different particle IDs")
        for atom_id in first[step]:
            for field in ("x", "y", "z", "fx", "fy", "fz", "vx", "vy", "vz"):
                maximum = max(
                    maximum,
                    abs(first[step][atom_id][field] - second[step][atom_id][field]),
                )
    return maximum


def write_rows(path, legacy_rows, analytical_rows):
    legacy_by_step = {row[0]: row for row in legacy_rows}
    analytical_by_step = {row[0]: row for row in analytical_rows}
    with path.open("w", newline="") as output:
        writer = csv.writer(output)
        writer.writerow((
            "step", "normalized_overlap", "reference_force",
            "legacy_force", "legacy_normalized_error",
            "analytical_force", "analytical_normalized_error",
        ))
        for step in sorted(legacy_by_step):
            legacy = legacy_by_step[step]
            analytical = analytical_by_step[step]
            writer.writerow((
                step, analytical[1], analytical[3], legacy[2], legacy[4],
                analytical[2], analytical[4],
            ))


def main():
    if not LIGGGHTS_BINARY.is_file():
        raise SystemExit("Set LIGGGHTS_BINARY near the top of this script")
    if RESULTS_DIRECTORY.exists():
        shutil.rmtree(RESULTS_DIRECTORY)
    RESULTS_DIRECTORY.mkdir()

    results = {}
    for ranks in MPI_RANKS:
        for model in ("legacy", "analytical"):
            snapshots, finalize_warning = run_case(model, ranks)
            results[(model, ranks)] = snapshots
            results[(model, ranks, "finalize_warning")] = finalize_warning

    legacy_error, legacy_rows = force_errors(results[("legacy", 1)])
    analytical_error, analytical_rows = force_errors(results[("analytical", 1)])
    if analytical_error > ANALYTICAL_TOLERANCE:
        raise AssertionError(
            "Analytical normalized force error {} exceeds {}".format(
                analytical_error, ANALYTICAL_TOLERANCE
            )
        )

    mpi_differences = {}
    for model in ("legacy", "analytical"):
        difference = maximum_snapshot_difference(results[(model, 1)], results[(model, 2)])
        mpi_differences[model] = difference
        if difference > MPI_TOLERANCE:
            raise AssertionError("{} MPI difference {}".format(model, difference))

    write_rows(RESULTS_DIRECTORY / "force_errors.csv", legacy_rows, analytical_rows)
    summary = """# JKR force-curve validation

- Analytical maximum normalized force error: `{:.6e}`
- Legacy maximum normalized force error: `{:.6e}`
- Legacy 1-rank/2-rank maximum state difference: `{:.6e}`
- Analytical 1-rank/2-rank maximum state difference: `{:.6e}`
- Analytical tolerance: `{:.6e}`
- MPI tolerance: `{:.6e}`
- macOS MPICH finalization warning after completed 2-rank runs: `{}`
""".format(
        analytical_error,
        legacy_error,
        mpi_differences["legacy"],
        mpi_differences["analytical"],
        ANALYTICAL_TOLERANCE,
        MPI_TOLERANCE,
        any(results[(model, 2, "finalize_warning")] for model in ("legacy", "analytical")),
    )
    (RESULTS_DIRECTORY / "summary.md").write_text(summary)
    print(summary)


if __name__ == "__main__":
    main()
