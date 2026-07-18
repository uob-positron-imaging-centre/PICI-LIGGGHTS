#!/usr/bin/env python3
"""Run both JKR solvers and plot total and JKR-only normal forces."""

from pathlib import Path
import csv
import math
import shutil
import subprocess


CASE_DIRECTORY = Path(__file__).resolve().parent
ROOT = CASE_DIRECTORY.parents[2]
BINARY_CANDIDATES = (
    ROOT / "build" / "liggghts",
    ROOT / "src" / "lmp_auto",
    Path("/tmp/pici_liggghts_build_cxx14/liggghts"),
    Path("/tmp/pici_liggghts_build_mpich3/liggghts"),
)
LIGGGHTS_BINARY = next((path for path in BINARY_CANDIDATES if path.is_file()), BINARY_CANDIDATES[0])
RESULTS_DIRECTORY = CASE_DIRECTORY / "force-component-results"

TIMESTEP = 1.0e-6
TURNAROUND_STEP = 3000
FINAL_STEP = 6000
PARTICLE_RADIUS = 0.005
YOUNGS_MODULUS = 5.0e6
POISSONS_RATIO = 0.30
WORK_OF_ADHESION = 0.10
def run_case(model, jkr_only=False):
    source_input = CASE_DIRECTORY / "in.force_components_{}".format(model)
    input_path = source_input
    suffix = "{}_jkr".format(model) if jkr_only else model
    if jkr_only:
        generated_input = source_input.read_text()
        generated_input = generated_input.replace(
            "coefficientRestitution peratomtypepair 1 0.50",
            "coefficientRestitution peratomtypepair 1 1.0",
        )
        generated_input = generated_input.replace(
            "force-component-results/{}_state.dump".format(model),
            "force-component-results/{}_jkr_state.dump".format(model),
        )
        generated_input = generated_input.replace(
            "force-component-results/{}_contact.dump".format(model),
            "force-component-results/{}_jkr_contact.dump".format(model),
        )
        input_path = RESULTS_DIRECTORY / "generated_{}_jkr_only.in".format(model)
        input_path.write_text(generated_input)
    command = [
        str(LIGGGHTS_BINARY),
        "-in",
        str(input_path),
    ]
    result = subprocess.run(
        command,
        cwd=CASE_DIRECTORY,
        text=True,
        errors="replace",
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    (RESULTS_DIRECTORY / "{}_screen.txt".format(suffix)).write_text(result.stdout)
    if result.returncode != 0:
        raise RuntimeError("{} failed; see its screen output".format(model))
    if "nan" in result.stdout.lower():
        raise RuntimeError("{} produced NaN output".format(model))


def parse_atom_dump(path):
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


def assemble_rows(legacy, analytical, legacy_jkr, analytical_jkr):
    if not (set(legacy) == set(analytical) == set(legacy_jkr) == set(analytical_jkr)):
        raise AssertionError("Models produced different timesteps")
    if min(legacy) != 0 or max(legacy) != FINAL_STEP or len(legacy) != FINAL_STEP + 1:
        raise AssertionError("Expected every timestep from 0 through 6000")
    delta_c, unused_fc = physical_scales()
    rows = []
    for step in sorted(legacy):
        old_atoms = legacy[step]
        new_atoms = analytical[step]
        old_jkr_atoms = legacy_jkr[step]
        new_jkr_atoms = analytical_jkr[step]
        if not all(set(atoms) == {1, 2} for atoms in (old_atoms, new_atoms, old_jkr_atoms, new_jkr_atoms)):
            raise AssertionError("Particle IDs changed")
        old_overlap = 2.0 * PARTICLE_RADIUS - (old_atoms[2]["x"] - old_atoms[1]["x"])
        new_overlap = 2.0 * PARTICLE_RADIUS - (new_atoms[2]["x"] - new_atoms[1]["x"])
        old_jkr_overlap = 2.0 * PARTICLE_RADIUS - (old_jkr_atoms[2]["x"] - old_jkr_atoms[1]["x"])
        new_jkr_overlap = 2.0 * PARTICLE_RADIUS - (new_jkr_atoms[2]["x"] - new_jkr_atoms[1]["x"])
        if max(abs(old_overlap - value) for value in (new_overlap, old_jkr_overlap, new_jkr_overlap)) > 1.0e-15:
            raise AssertionError("Prescribed overlaps differ between models")
        normalized_overlap = old_overlap / delta_c
        legacy_total = old_atoms[2]["fx"]
        analytical_total = new_atoms[2]["fx"]
        legacy_jkr_force = old_jkr_atoms[2]["fx"]
        analytical_jkr_force = new_jkr_atoms[2]["fx"]
        rows.append({
            "step": step,
            "time": step * TIMESTEP,
            "phase": "unload" if step <= TURNAROUND_STEP else "reload",
            "overlap": old_overlap,
            "normalized_overlap": normalized_overlap,
            "legacy_total_force": legacy_total,
            "legacy_jkr_force": legacy_jkr_force,
            "legacy_damping_force": legacy_total - legacy_jkr_force,
            "analytical_total_force": analytical_total,
            "analytical_jkr_force": analytical_jkr_force,
            "analytical_damping_force": analytical_total - analytical_jkr_force,
            "total_force_difference": analytical_total - legacy_total,
            "jkr_force_difference": analytical_jkr_force - legacy_jkr_force,
        })
    return rows


def write_csv(path, rows):
    with path.open("w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def write_svg(path, rows, series, title, x_key, x_label, y_label, x_scale, y_scale):
    width, height = 1000, 600
    left, right, top, bottom = 105, 30, 55, 85
    plot_width = width - left - right
    plot_height = height - top - bottom
    xs = [row[x_key] * x_scale for row in rows]
    ys = [row[key] * y_scale for key, unused_label, unused_color in series for row in rows]
    x_min, x_max = min(xs), max(xs)
    y_min, y_max = min(ys), max(ys)
    x_pad = 0.04 * (x_max - x_min or 1.0)
    y_pad = 0.08 * (y_max - y_min or 1.0)
    x_min, x_max = x_min - x_pad, x_max + x_pad
    y_min, y_max = y_min - y_pad, y_max + y_pad

    def sx(value):
        return left + plot_width * (value - x_min) / (x_max - x_min)

    def sy(value):
        return top + plot_height * (y_max - value) / (y_max - y_min)

    content = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="{}" height="{}">'.format(width, height),
        '<rect width="100%" height="100%" fill="white"/>',
        '<text x="{}" y="30" text-anchor="middle" font-family="sans-serif" font-size="21">{}</text>'.format(width / 2, title),
        '<rect x="{}" y="{}" width="{}" height="{}" fill="none" stroke="#222"/>'.format(left, top, plot_width, plot_height),
    ]
    for tick in range(6):
        x_value = x_min + (x_max - x_min) * tick / 5.0
        y_value = y_min + (y_max - y_min) * tick / 5.0
        content.append('<text x="{:.2f}" y="{}" text-anchor="middle" font-family="sans-serif" font-size="12">{:.4g}</text>'.format(sx(x_value), top + plot_height + 23, x_value))
        content.append('<text x="{}" y="{:.2f}" text-anchor="end" dominant-baseline="middle" font-family="sans-serif" font-size="12">{:.4g}</text>'.format(left - 9, sy(y_value), y_value))
    if x_min <= 0.0 <= x_max:
        content.append('<line x1="{0:.2f}" y1="{1}" x2="{0:.2f}" y2="{2}" stroke="#bbb"/>'.format(sx(0.0), top, top + plot_height))
    if y_min <= 0.0 <= y_max:
        content.append('<line x1="{0}" y1="{1:.2f}" x2="{2}" y2="{1:.2f}" stroke="#bbb"/>'.format(left, sy(0.0), left + plot_width))
    content.append('<text x="{}" y="{}" text-anchor="middle" font-family="sans-serif">{}</text>'.format(left + plot_width / 2, height - 24, x_label))
    content.append('<text x="22" y="{}" text-anchor="middle" transform="rotate(-90 22 {})" font-family="sans-serif">{}</text>'.format(top + plot_height / 2, top + plot_height / 2, y_label))
    dash_patterns = ("", "8,5", "3,3", "10,4,2,4")
    for series_index, (key, label, color) in enumerate(series):
        points = " ".join(
            "{:.2f},{:.2f}".format(sx(row[x_key] * x_scale), sy(row[key] * y_scale))
            for row in rows
        )
        dash = dash_patterns[series_index % len(dash_patterns)]
        dash_attribute = ' stroke-dasharray="{}"'.format(dash) if dash else ""
        content.append('<polyline fill="none" stroke="{}" stroke-width="1.8"{} points="{}"/>'.format(color, dash_attribute, points))
        legend_x = left + 10 + (series_index % 2) * 300
        legend_y = top + 18 + (series_index // 2) * 22
        content.append('<line x1="{}" y1="{}" x2="{}" y2="{}" stroke="{}" stroke-width="2"{}/>'.format(legend_x, legend_y, legend_x + 32, legend_y, color, dash_attribute))
        content.append('<text x="{}" y="{}" font-family="sans-serif" font-size="13">{}</text>'.format(legend_x + 39, legend_y + 4, label))
    content.append("</svg>")
    path.write_text("\n".join(content))


def main():
    if not LIGGGHTS_BINARY.is_file():
        raise SystemExit("Set LIGGGHTS_BINARY near the top of this script")
    if RESULTS_DIRECTORY.exists():
        shutil.rmtree(RESULTS_DIRECTORY)
    RESULTS_DIRECTORY.mkdir()
    for model in ("legacy", "analytical"):
        run_case(model, jkr_only=False)
        run_case(model, jkr_only=True)

    legacy = parse_atom_dump(RESULTS_DIRECTORY / "legacy_state.dump")
    analytical = parse_atom_dump(RESULTS_DIRECTORY / "analytical_state.dump")
    legacy_jkr = parse_atom_dump(RESULTS_DIRECTORY / "legacy_jkr_state.dump")
    analytical_jkr = parse_atom_dump(RESULTS_DIRECTORY / "analytical_jkr_state.dump")
    rows = assemble_rows(legacy, analytical, legacy_jkr, analytical_jkr)
    write_csv(RESULTS_DIRECTORY / "force_components.csv", rows)
    write_svg(
        RESULTS_DIRECTORY / "force_components_vs_overlap.svg",
        rows,
        (
            ("legacy_total_force", "legacy total", "#2166ac"),
            ("analytical_total_force", "analytical total", "#b2182b"),
            ("legacy_jkr_force", "legacy JKR only", "#1b9e77"),
            ("analytical_jkr_force", "analytical JKR only", "#e08214"),
        ),
        "JKR lookup-table and analytical force components",
        "overlap", "overlap (micrometers)", "force (mN)", 1.0e6, 1.0e3,
    )
    write_svg(
        RESULTS_DIRECTORY / "force_components_vs_time.svg",
        rows,
        (
            ("legacy_total_force", "legacy total", "#2166ac"),
            ("analytical_total_force", "analytical total", "#b2182b"),
            ("legacy_jkr_force", "legacy JKR only", "#1b9e77"),
            ("analytical_jkr_force", "analytical JKR only", "#e08214"),
        ),
        "JKR force components through unload and reload",
        "time", "time (ms)", "force (mN)", 1.0e3, 1.0e3,
    )
    write_svg(
        RESULTS_DIRECTORY / "force_differences_vs_overlap.svg",
        rows,
        (
            ("total_force_difference", "analytical - legacy total", "#762a83"),
            ("jkr_force_difference", "analytical - legacy JKR", "#e66101"),
        ),
        "Analytical minus lookup-table force",
        "overlap", "overlap (micrometers)", "force difference (microN)", 1.0e6, 1.0e6,
    )

    maximum_total_difference = max(abs(row["total_force_difference"]) for row in rows)
    maximum_jkr_difference = max(abs(row["jkr_force_difference"]) for row in rows)
    summary = """# JKR force-component comparison

- Saved timesteps: `{}` (every step from 0 through 6000).
- Overlap range: `{:.6e} m` to `{:.6e} m`.
- Maximum absolute total-force difference: `{:.6e} N`.
- Maximum absolute JKR-only force difference: `{:.6e} N`.
- Total-force runs use restitution `0.5`; JKR-only runs use restitution `1.0`, which makes normal damping exactly zero.
- Raw `pair/gran/local` dumps retain overlap, force vector, and geometrically normal force vector for every active contact in all four runs.

The CSV contains every prescribed overlap, directly measured total force, directly measured JKR-only force, damping remainder, and analytical-minus-legacy difference.
""".format(
        len(rows),
        min(row["overlap"] for row in rows),
        max(row["overlap"] for row in rows),
        maximum_total_difference,
        maximum_jkr_difference,
    )
    (RESULTS_DIRECTORY / "summary.md").write_text(summary)
    print(summary)


if __name__ == "__main__":
    main()
