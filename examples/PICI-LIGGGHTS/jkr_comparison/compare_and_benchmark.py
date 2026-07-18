#!/usr/bin/env python3
"""Run the common-state drum comparison and stable JKR benchmarks."""

from pathlib import Path
import csv
import math
import re
import shutil
import statistics
import subprocess
import time


CASE_DIRECTORY = Path(__file__).resolve().parent
ROOT = CASE_DIRECTORY.parents[2]
BINARY_CANDIDATES = (
    ROOT / "build" / "liggghts",
    ROOT / "src" / "lmp_auto",
    Path("/tmp/pici_liggghts_build_cxx14/liggghts"),
    Path("/tmp/pici_liggghts_build_mpich3/liggghts"),
)
LIGGGHTS_BINARY = next((path for path in BINARY_CANDIDATES if path.is_file()), BINARY_CANDIDATES[0])
RESULTS_DIRECTORY = CASE_DIRECTORY / "comparison-results"
BENCHMARK_DIRECTORY = CASE_DIRECTORY / "benchmark-results"

TIMESTEP = 1.0e-5
DRUM_STEPS = 100000
SETTLE_STEPS = 100000
BENCHMARK_STEPS = 20000
STARTUP_REPETITIONS = 5
FORCE_REPETITIONS = 5
EXPECTED_PARTICLES = 600
MODELS = ("legacy", "analytical")


def run_liggghts(input_name, output_path, measure_memory=False):
    command = [str(LIGGGHTS_BINARY), "-in", str(CASE_DIRECTORY / input_name)]
    if measure_memory:
        command = ["/usr/bin/time", "-l"] + command
    started = time.perf_counter()
    result = subprocess.run(
        command,
        cwd=CASE_DIRECTORY,
        text=True,
        errors="replace",
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    wall_time = time.perf_counter() - started
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(result.stdout)
    macos_mpich_finalize_error = (
        result.returncode == 143 and
        "Fatal error in internal_Finalize" in result.stdout and
        "Loop time of" in result.stdout
    )
    if result.returncode != 0 and not macos_mpich_finalize_error:
        raise RuntimeError("{} failed; see {}".format(input_name, output_path))
    if "nan" in result.stdout.lower():
        raise RuntimeError("{} produced NaN output".format(input_name))
    return result.stdout, wall_time


def parse_dump(path):
    snapshots = {}
    lines = path.read_text().splitlines()
    index = 0
    while index < len(lines):
        if lines[index] != "ITEM: TIMESTEP":
            raise RuntimeError("Unexpected dump format at line {} in {}".format(index + 1, path))
        step = int(lines[index + 1])
        atom_count = int(lines[index + 3])
        header_index = index + 8
        columns = lines[header_index].split()[2:]
        atoms = {}
        for row_index in range(header_index + 1, header_index + 1 + atom_count):
            values = lines[row_index].split()
            atom = {name: float(value) for name, value in zip(columns, values)}
            atom_id = int(atom["id"])
            atoms[atom_id] = atom
        snapshots[step] = atoms
        index = header_index + 1 + atom_count
    return snapshots


def particle_mass(atom):
    return atom["density"] * 4.0 * math.pi * atom["radius"] ** 3 / 3.0


def vector_distance(first, second, fields):
    return math.sqrt(sum((first[field] - second[field]) ** 2 for field in fields))


def snapshot_statistics(legacy, analytical):
    if set(legacy) != set(analytical):
        raise AssertionError("Legacy and analytical runs have different particle IDs")
    position_errors = []
    velocity_errors = []
    legacy_mass = 0.0
    analytical_mass = 0.0
    legacy_com = [0.0, 0.0, 0.0]
    analytical_com = [0.0, 0.0, 0.0]
    legacy_ke = 0.0
    analytical_ke = 0.0
    for atom_id in legacy:
        old = legacy[atom_id]
        new = analytical[atom_id]
        position_errors.append(vector_distance(old, new, ("x", "y", "z")))
        velocity_errors.append(vector_distance(old, new, ("vx", "vy", "vz")))
        old_mass = particle_mass(old)
        new_mass = particle_mass(new)
        legacy_mass += old_mass
        analytical_mass += new_mass
        for component, field in enumerate(("x", "y", "z")):
            legacy_com[component] += old_mass * old[field]
            analytical_com[component] += new_mass * new[field]
        legacy_ke += 0.5 * old_mass * sum(old[field] ** 2 for field in ("vx", "vy", "vz"))
        analytical_ke += 0.5 * new_mass * sum(new[field] ** 2 for field in ("vx", "vy", "vz"))
    legacy_com = [value / legacy_mass for value in legacy_com]
    analytical_com = [value / analytical_mass for value in analytical_com]
    return {
        "rms_position": math.sqrt(sum(value ** 2 for value in position_errors) / len(position_errors)),
        "max_position": max(position_errors),
        "rms_velocity": math.sqrt(sum(value ** 2 for value in velocity_errors) / len(velocity_errors)),
        "com_difference": math.sqrt(sum((a - b) ** 2 for a, b in zip(legacy_com, analytical_com))),
        "legacy_ke": legacy_ke,
        "analytical_ke": analytical_ke,
        "ke_difference": analytical_ke - legacy_ke,
    }


def write_metrics(path, legacy_snapshots, analytical_snapshots):
    if set(legacy_snapshots) != set(analytical_snapshots):
        raise AssertionError("Legacy and analytical dumps have different timesteps")
    rows = []
    for step in sorted(legacy_snapshots):
        if len(legacy_snapshots[step]) != EXPECTED_PARTICLES:
            raise AssertionError("Particle count changed at step {}".format(step))
        values = snapshot_statistics(legacy_snapshots[step], analytical_snapshots[step])
        values["step"] = step
        values["time"] = step * TIMESTEP
        rows.append(values)
    fields = (
        "step", "time", "rms_position", "max_position", "rms_velocity",
        "com_difference", "legacy_ke", "analytical_ke", "ke_difference",
    )
    with path.open("w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    return rows


def svg_plot(path, rows, series, title, y_label):
    width, height = 900, 520
    left, right, top, bottom = 85, 25, 45, 70
    plot_width = width - left - right
    plot_height = height - top - bottom
    xs = [row["time"] for row in rows]
    values = [row[key] for key, unused_label, unused_color in series for row in rows]
    x_max = max(xs) if xs else 1.0
    y_max = max(values) if values and max(values) > 0.0 else 1.0
    colors = {key: color for key, unused_label, color in series}
    labels = {key: label for key, label, unused_color in series}

    def sx(value):
        return left + plot_width * value / x_max

    def sy(value):
        return top + plot_height * (1.0 - value / y_max)

    content = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="{}" height="{}" viewBox="0 0 {} {}">'.format(width, height, width, height),
        '<rect width="100%" height="100%" fill="white"/>',
        '<text x="{}" y="28" text-anchor="middle" font-family="sans-serif" font-size="20">{}</text>'.format(width / 2, title),
        '<line x1="{}" y1="{}" x2="{}" y2="{}" stroke="#222"/>'.format(left, top + plot_height, left + plot_width, top + plot_height),
        '<line x1="{}" y1="{}" x2="{}" y2="{}" stroke="#222"/>'.format(left, top, left, top + plot_height),
    ]
    for tick in range(6):
        x_value = x_max * tick / 5.0
        y_value = y_max * tick / 5.0
        content.append('<text x="{:.1f}" y="{}" text-anchor="middle" font-family="sans-serif" font-size="12">{:.3g}</text>'.format(sx(x_value), top + plot_height + 22, x_value))
        content.append('<text x="{}" y="{:.1f}" text-anchor="end" dominant-baseline="middle" font-family="sans-serif" font-size="12">{:.3g}</text>'.format(left - 8, sy(y_value), y_value))
    content.append('<text x="{}" y="{}" text-anchor="middle" font-family="sans-serif">time (s)</text>'.format(left + plot_width / 2, height - 20))
    content.append('<text x="20" y="{}" text-anchor="middle" transform="rotate(-90 20 {})" font-family="sans-serif">{}</text>'.format(top + plot_height / 2, top + plot_height / 2, y_label))
    for series_index, (key, unused_label, color) in enumerate(series):
        points = " ".join("{:.2f},{:.2f}".format(sx(row["time"]), sy(row[key])) for row in rows)
        content.append('<polyline fill="none" stroke="{}" stroke-width="2" points="{}"/>'.format(color, points))
        legend_x = left + 15 + series_index * 190
        content.append('<line x1="{}" y1="{}" x2="{}" y2="{}" stroke="{}" stroke-width="3"/>'.format(legend_x, top + 12, legend_x + 25, top + 12, colors[key]))
        content.append('<text x="{}" y="{}" font-family="sans-serif" font-size="13">{}</text>'.format(legend_x + 32, top + 16, labels[key]))
    content.append("</svg>")
    path.write_text("\n".join(content))


def overlay_svg(path, legacy, analytical):
    width, height = 700, 700
    margin = 40
    values = [atom[field] for state in (legacy, analytical) for atom in state.values() for field in ("x", "z")]
    limit = max(0.05, max(abs(value) for value in values))

    def transform(x, z):
        scale = (width - 2 * margin) / (2 * limit)
        return width / 2 + x * scale, height / 2 - z * scale

    content = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="{}" height="{}">'.format(width, height),
        '<rect width="100%" height="100%" fill="white"/>',
        '<text x="350" y="24" text-anchor="middle" font-family="sans-serif" font-size="18">Final particle positions (x-z projection)</text>',
    ]
    for atom_id in sorted(legacy):
        old_x, old_z = transform(legacy[atom_id]["x"], legacy[atom_id]["z"])
        new_x, new_z = transform(analytical[atom_id]["x"], analytical[atom_id]["z"])
        content.append('<line x1="{:.2f}" y1="{:.2f}" x2="{:.2f}" y2="{:.2f}" stroke="#999" stroke-width="0.6"/>'.format(old_x, old_z, new_x, new_z))
        content.append('<circle cx="{:.2f}" cy="{:.2f}" r="1.8" fill="#2166ac"/>'.format(old_x, old_z))
        content.append('<circle cx="{:.2f}" cy="{:.2f}" r="1.8" fill="#b2182b"/>'.format(new_x, new_z))
    content.extend([
        '<circle cx="55" cy="675" r="4" fill="#2166ac"/><text x="65" y="680" font-family="sans-serif">legacy</text>',
        '<circle cx="145" cy="675" r="4" fill="#b2182b"/><text x="155" y="680" font-family="sans-serif">analytical</text>',
        "</svg>",
    ])
    path.write_text("\n".join(content))


def write_vtk(path, atoms, title):
    ordered = [atoms[atom_id] for atom_id in sorted(atoms)]
    lines = [
        "# vtk DataFile Version 3.0", title, "ASCII", "DATASET POLYDATA",
        "POINTS {} double".format(len(ordered)),
    ]
    lines.extend("{:.17g} {:.17g} {:.17g}".format(atom["x"], atom["y"], atom["z"]) for atom in ordered)
    lines.extend(("VERTICES {} {}".format(len(ordered), 2 * len(ordered)),))
    lines.extend("1 {}".format(index) for index in range(len(ordered)))
    lines.extend(("POINT_DATA {}".format(len(ordered)), "SCALARS particle_id int 1", "LOOKUP_TABLE default"))
    lines.extend(str(int(atom["id"])) for atom in ordered)
    lines.extend(("VECTORS velocity double",))
    lines.extend("{:.17g} {:.17g} {:.17g}".format(atom["vx"], atom["vy"], atom["vz"]) for atom in ordered)
    path.write_text("\n".join(lines) + "\n")


def parse_loop_time(output):
    matches = re.findall(r"Loop time of ([0-9.eE+-]+) on", output)
    if not matches:
        raise RuntimeError("Could not find LIGGGHTS loop time")
    return float(matches[-1])


def parse_contact_count(output):
    match = re.search(r"CONTACT_COUNT\s+([0-9.eE+-]+)", output)
    if not match:
        raise RuntimeError("Could not find contact count")
    return float(match.group(1))


def parse_peak_rss(output):
    match = re.search(r"([0-9]+)\s+maximum resident set size", output)
    return int(match.group(1)) if match else 0


def ordered_models(repetition):
    return MODELS if repetition % 2 == 0 else tuple(reversed(MODELS))


def run_benchmarks():
    startup_rows = []
    for repetition in range(STARTUP_REPETITIONS):
        for model in ordered_models(repetition):
            output, wall_time = run_liggghts(
                "in.startup_{}".format(model),
                BENCHMARK_DIRECTORY / "startup_{}_{}.txt".format(model, repetition + 1),
                measure_memory=True,
            )
            startup_rows.append({
                "model": model, "repetition": repetition + 1,
                "order": ordered_models(repetition).index(model) + 1,
                "wall_time": wall_time, "peak_rss_bytes": parse_peak_rss(output),
            })

    for model in MODELS:
        run_liggghts(
            "in.benchmark_{}".format(model),
            BENCHMARK_DIRECTORY / "force_warmup_{}.txt".format(model),
        )

    force_rows = []
    for repetition in range(FORCE_REPETITIONS):
        for model in ordered_models(repetition):
            output, wall_time = run_liggghts(
                "in.benchmark_{}".format(model),
                BENCHMARK_DIRECTORY / "force_{}_{}.txt".format(model, repetition + 1),
            )
            loop_time = parse_loop_time(output)
            contacts = parse_contact_count(output)
            force_rows.append({
                "model": model, "repetition": repetition + 1,
                "order": ordered_models(repetition).index(model) + 1,
                "wall_time": wall_time, "loop_time": loop_time,
                "timesteps_per_second": BENCHMARK_STEPS / loop_time,
                "contact_count": contacts,
                "contacts_per_second": BENCHMARK_STEPS * contacts / loop_time,
            })

    for filename, rows in (("startup.csv", startup_rows), ("force_loop.csv", force_rows)):
        with (BENCHMARK_DIRECTORY / filename).open("w", newline="") as output:
            writer = csv.DictWriter(output, fieldnames=rows[0].keys())
            writer.writeheader()
            writer.writerows(rows)
    return startup_rows, force_rows


def median_for(rows, model, key):
    return statistics.median(row[key] for row in rows if row["model"] == model)


def spread_for(rows, model, key):
    values = [row[key] for row in rows if row["model"] == model]
    return max(values) - min(values)


def main():
    if not LIGGGHTS_BINARY.is_file():
        raise SystemExit("Set LIGGGHTS_BINARY near the top of this script")
    RESULTS_DIRECTORY.mkdir(exist_ok=True)
    if BENCHMARK_DIRECTORY.exists():
        shutil.rmtree(BENCHMARK_DIRECTORY)
    BENCHMARK_DIRECTORY.mkdir()
    for model in MODELS:
        (RESULTS_DIRECTORY / model).mkdir(exist_ok=True)

    checkpoint = RESULTS_DIRECTORY / "common_settled.data"
    prepare_screen = RESULTS_DIRECTORY / "prepare_screen.txt"
    if checkpoint.is_file() and prepare_screen.is_file():
        preparation_output = prepare_screen.read_text(errors="replace")
    else:
        preparation_output, unused_preparation_wall = run_liggghts(
            "in.prepare_drum", prepare_screen
        )
    if not re.search(r"\n\s*100000\s+600\s", preparation_output):
        raise AssertionError("Checkpoint preparation did not finish with 600 particles")

    drum_loop_times = {}
    for model in MODELS:
        dump_path = RESULTS_DIRECTORY / model / "trajectory.dump"
        screen_path = RESULTS_DIRECTORY / model / "screen.txt"
        complete = False
        if dump_path.is_file() and screen_path.is_file():
            try:
                complete = max(parse_dump(dump_path)) == DRUM_STEPS
            except (RuntimeError, ValueError):
                complete = False
        if complete:
            drum_output = screen_path.read_text(errors="replace")
        else:
            if dump_path.exists():
                dump_path.unlink()
            drum_output, unused_wall_time = run_liggghts(
                "in.drum_{}".format(model), screen_path
            )
        drum_loop_times[model] = parse_loop_time(drum_output)

    legacy = parse_dump(RESULTS_DIRECTORY / "legacy" / "trajectory.dump")
    analytical = parse_dump(RESULTS_DIRECTORY / "analytical" / "trajectory.dump")
    rows = write_metrics(RESULTS_DIRECTORY / "divergence.csv", legacy, analytical)
    final_step = max(legacy)
    if final_step != DRUM_STEPS:
        raise AssertionError("Drum runs did not complete {} steps".format(DRUM_STEPS))
    final_legacy = legacy[final_step]
    final_analytical = analytical[final_step]
    svg_plot(
        RESULTS_DIRECTORY / "position_divergence.svg", rows,
        (("rms_position", "RMS position", "#2166ac"), ("max_position", "maximum position", "#b2182b")),
        "JKR particle-position divergence", "distance (m)",
    )
    svg_plot(
        RESULTS_DIRECTORY / "velocity_divergence.svg", rows,
        (("rms_velocity", "RMS velocity", "#4d9221"),),
        "JKR particle-velocity divergence", "velocity difference (m/s)",
    )
    svg_plot(
        RESULTS_DIRECTORY / "center_of_mass_divergence.svg", rows,
        (("com_difference", "center of mass", "#762a83"),),
        "JKR center-of-mass divergence", "distance (m)",
    )
    energy_rows = [dict(row, absolute_ke_difference=abs(row["ke_difference"])) for row in rows]
    svg_plot(
        RESULTS_DIRECTORY / "kinetic_energy_divergence.svg", energy_rows,
        (("absolute_ke_difference", "absolute KE difference", "#e08214"),),
        "JKR kinetic-energy divergence", "energy difference (J)",
    )
    overlay_svg(RESULTS_DIRECTORY / "final_position_overlay.svg", final_legacy, final_analytical)
    write_vtk(RESULTS_DIRECTORY / "legacy_final.vtk", final_legacy, "Legacy JKR final state")
    write_vtk(RESULTS_DIRECTORY / "analytical_final.vtk", final_analytical, "Analytical JKR final state")

    startup_rows, force_rows = run_benchmarks()
    startup_legacy = median_for(startup_rows, "legacy", "wall_time")
    startup_analytical = median_for(startup_rows, "analytical", "wall_time")
    force_legacy = median_for(force_rows, "legacy", "timesteps_per_second")
    force_analytical = median_for(force_rows, "analytical", "timesteps_per_second")
    contacts_legacy = median_for(force_rows, "legacy", "contacts_per_second")
    contacts_analytical = median_for(force_rows, "analytical", "contacts_per_second")
    final = rows[-1]
    summary = """# Legacy and analytical JKR comparison

## Simulation

- Common checkpoint: {} particles, non-adhesive settling for {} steps.
- Comparison: 30 rpm, `dt = 1e-5 s`, `workOfAdhesion = 0.3 J/m^2`, {} steps (half a revolution).
- All immutable particle IDs were retained; no NaNs or lost-particle solver errors occurred.
- LIGGGHTS half-revolution loop time: legacy `{:.6g} s`, analytical `{:.6g} s`.

## Final divergence

- RMS position difference: `{:.6e} m`.
- Maximum position difference: `{:.6e} m`.
- RMS velocity difference: `{:.6e} m/s`.
- Center-of-mass difference: `{:.6e} m`.
- Kinetic-energy difference (analytical - legacy): `{:.6e} J`.

Trajectory divergence is a measured outcome of the chaotic many-particle system, not a force-law failure.

## Benchmarks

- Five fresh-process `run 0` starts: legacy median `{:.6g} s` (spread `{:.3g} s`, peak RSS median `{:.3f} MiB`), analytical median `{:.6g} s` (spread `{:.3g} s`, peak RSS median `{:.3f} MiB`).
- Median legacy-minus-analytical startup difference (the observed LUT-related overhead): `{:.6g} s`.
- Frozen geometry, one warm-up plus five measured 20,000-step runs per model.
- Legacy median: `{:.6g} timestep/s` (spread `{:.3g}`), `{:.6g} contact-evaluation/s`.
- Analytical median: `{:.6g} timestep/s` (spread `{:.3g}`), `{:.6g} contact-evaluation/s`.
- Analytical/legacy steady-state speed ratio: `{:.6g}`.

Raw timings, order, peak RSS, and spreads are in `../benchmark-results/startup.csv` and `../benchmark-results/force_loop.csv`.
""".format(
        EXPECTED_PARTICLES, SETTLE_STEPS, DRUM_STEPS,
        drum_loop_times["legacy"], drum_loop_times["analytical"],
        final["rms_position"], final["max_position"], final["rms_velocity"],
        final["com_difference"], final["ke_difference"],
        startup_legacy, spread_for(startup_rows, "legacy", "wall_time"),
        median_for(startup_rows, "legacy", "peak_rss_bytes") / (1024.0 ** 2),
        startup_analytical, spread_for(startup_rows, "analytical", "wall_time"),
        median_for(startup_rows, "analytical", "peak_rss_bytes") / (1024.0 ** 2),
        startup_legacy - startup_analytical,
        force_legacy, spread_for(force_rows, "legacy", "timesteps_per_second"), contacts_legacy,
        force_analytical, spread_for(force_rows, "analytical", "timesteps_per_second"), contacts_analytical,
        force_analytical / force_legacy,
    )
    (RESULTS_DIRECTORY / "summary.md").write_text(summary)
    print(summary)


if __name__ == "__main__":
    main()
