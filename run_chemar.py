import argparse
from bisect import bisect_left
import json
import math
from pathlib import Path
import subprocess


def find_executable():
    project_directory = Path(__file__).resolve().parent
    build_directory = project_directory / "build"
    for folder in (build_directory, build_directory / "Release", project_directory):
        for name in ("chemar.exe", "chemar"):
            candidate = folder / name
            if candidate.is_file():
                return candidate
    raise FileNotFoundError("Сначала соберите CHEMAR по инструкции из README.md.")


def add_species_legend(axis):
    if not axis.lines:
        return
    legend = axis.legend(
        loc="upper center", bbox_to_anchor=(0.5, -0.15), ncols=min(5, len(axis.lines))
    )
    legend_entries = {}
    for curve, legend_line, legend_label in zip(axis.lines, legend.get_lines(), legend.get_texts()):
        legend_line.set_picker(5)
        legend_label.set_picker(True)
        entry = (curve, legend_line, legend_label)
        legend_entries[legend_line] = entry
        legend_entries[legend_label] = entry

    canvas = axis.figure.canvas

    def toggle_species(event):
        entry = legend_entries.get(event.artist)
        if entry is None:
            return
        curve, legend_line, legend_label = entry
        visible = not curve.get_visible()
        curve.set_visible(visible)
        for item in (legend_line, legend_label):
            item.set_alpha(1.0 if visible else 0.3)
        canvas.draw_idle()

    canvas.mpl_connect("pick_event", toggle_species)


def plot_results(path, save_only, plot_until=None):
    import matplotlib

    if save_only:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    with path.open(encoding="utf-8") as stream:
        document = json.load(stream)
    states = document["states"]
    if not states:
        raise ValueError("В results.json нет ни одного состояния.")

    coordinate_is_distance = document["coordinate"]["name"] == "distance"
    coordinate_field = "distance_m" if coordinate_is_distance else "time_s"
    coordinates = [state[coordinate_field] for state in states]
    temperature_k = [state["temperature_k"] for state in states]
    pressure_pa = [state["pressure_pa"] for state in states]
    plot_end = coordinates[-1]
    if plot_until is not None:
        plot_end = min(plot_until, plot_end)
    plot_row_count = bisect_left(coordinates, plot_end) + 1
    if coordinate_is_distance:
        coordinate_scale, horizontal_label = 1, "Расстояние, м"
    elif plot_end and plot_end <= 1e-3:
        coordinate_scale, horizontal_label = 1e6, "Время, мкс"
    else:
        coordinate_scale, horizontal_label = 1e3, "Время, мс"
    plot_coordinates = [coordinate * coordinate_scale for coordinate in coordinates[:plot_row_count]]
    single_point_marker = "o" if plot_row_count == 1 else None

    _, temperature_axes = plt.subplots(num="Температура", figsize=(10, 4.8), layout="constrained")
    temperature_axes.plot(
        plot_coordinates, temperature_k[:plot_row_count], color="tab:red", marker=single_point_marker
    )

    _, pressure_axes = plt.subplots(num="Давление", figsize=(10, 4.8), layout="constrained")
    pressure_kpa = [value / 1000 for value in pressure_pa[:plot_row_count]]
    pressure_axes.plot(plot_coordinates, pressure_kpa, color="tab:blue", marker=single_point_marker)

    _, composition_axes = plt.subplots(
        num="Состав смеси", figsize=(10, 5.5), layout="constrained"
    )
    chemical_subscripts = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
    for species_name in document["species"]:
        molar_amounts = [state["composition_mol_kg"][species_name] for state in states]
        if any(molar_amounts):
            composition_axes.plot(
                plot_coordinates, molar_amounts[:plot_row_count],
                label=species_name.translate(chemical_subscripts), marker=single_point_marker
            )
    composition_axes.set_ylim(bottom=0)
    add_species_legend(composition_axes)

    suffix = "-detail" if plot_until is not None else ""
    plots = (
        ("temperature", temperature_axes, "Температура", "К"),
        ("pressure", pressure_axes, "Давление", "кПа"),
        ("composition", composition_axes, "Состав смеси", "Количество вещества, моль/кг")
    )
    for file_stem, axes, title, vertical_label in plots:
        axes.set_title(title, loc="left")
        axes.set_xlabel(horizontal_label)
        axes.set_ylabel(vertical_label)
        axes.grid(alpha=0.25)
        if plot_end > 0:
            axes.set_xlim(0, plot_end * coordinate_scale)
        if save_only:
            axes.figure.savefig(path.with_name(f"{file_stem}{suffix}.png"), dpi=160)

    if not save_only:
        plt.show()
    for _, axes, _, _ in plots:
        plt.close(axes.figure)


def main():
    parser = argparse.ArgumentParser(description="Запустить CHEMAR и построить графики.")
    parser.add_argument("data", type=Path, nargs="?", default=Path("data"))
    parser.add_argument("--output", type=Path, default=Path("output"))
    parser.add_argument("--steps", type=int)
    parser.add_argument("--save-only", action="store_true", help="Сохранить PNG без окон.")
    parser.add_argument("--plot-until", type=float, help="Показать начало до указанного времени в секундах или расстояния в метрах.")
    args = parser.parse_args()
    if args.plot_until is not None and not (math.isfinite(args.plot_until) and args.plot_until > 0):
        parser.error("--plot-until должен быть положительным конечным числом.")
    try:
        command = [find_executable(), "--data", args.data, "--output", args.output]
        if args.steps is not None:
            command.extend(["--steps", str(args.steps)])
        subprocess.run(command, check=True)
        plot_results(args.output / "results.json", args.save_only, args.plot_until)
    except subprocess.CalledProcessError as error:
        return error.returncode
    except (OSError, ImportError, KeyError, ValueError, json.JSONDecodeError) as error:
        parser.exit(1, f"{error}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
