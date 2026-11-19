import argparse
import csv
import math
from pathlib import Path
import matplotlib.pyplot as plt


def time_scale(time):
    end = time[-1]
    if end <= 2.0e-3:
        return [value * 1.0e6 for value in time], 'Время, мкс'
    if end <= 2.0:
        return [value * 1.0e3 for value in time], 'Время, мс'
    return time, 'Время, с'


def finish_axis(axis):
    axis.grid(True, alpha=.25)


def scalar_figure(time, values, time_label, title, ylabel, color, minimum_span=0):
    figure, axis = plt.subplots(figsize=(8, 5))
    axis.plot(time, values, color=color)
    if minimum_span > 0:
        lower = min(values)
        upper = max(values)
        center = (lower + upper) / 2
        span = max(upper - lower, minimum_span)
        axis.set_ylim(center - span / 2, center + span / 2)
    axis.ticklabel_format(axis='y', style='plain', useOffset=False)
    axis.set_xlabel(time_label)
    axis.set_ylabel(ylabel)
    axis.set_title(title)
    finish_axis(axis)
    figure.tight_layout()


def fraction_figure(time, columns, selected, colors, prefix,
                    time_label, title, ylabel):
    figure, axis = plt.subplots(figsize=(9, 6))
    for name in selected:
        axis.plot(time, columns[prefix + name], color=colors.get(name), label=name)
    axis.set_xlabel(time_label)
    axis.set_ylabel(ylabel)
    axis.set_title(title)
    axis.set_ylim(bottom=0)
    finish_axis(axis)
    axis.legend(ncol=3, frameon=False)
    figure.tight_layout()


def reaction_rate_figure(time, columns, reaction_count, suffix,
                         time_label, title, logarithmic):
    figure, axis = plt.subplots(figsize=(10, 6))
    color_map = plt.get_cmap('tab20')
    largest = 0

    for reaction in range(1, reaction_count + 1):
        values = columns[f'reaction_{reaction}_{suffix}_mol_m3_s']
        largest = max(largest, max(abs(value) for value in values))
        if logarithmic:
            values = [value if value > 0 else math.nan for value in values]
        axis.plot(time, values, color=color_map((reaction - 1) % 20), label=f'R{reaction}')

    if logarithmic:
        axis.set_yscale('log')
    else:
        axis.set_yscale('symlog', linthresh=max(largest * 1.0e-8, 1.0e-30))

    axis.set_xlabel(time_label)
    axis.set_ylabel('Скорость, моль/(м³·с)')
    axis.set_title(title)
    finish_axis(axis)
    axis.legend(ncol=2, frameon=False, fontsize=8,
                bbox_to_anchor=(1.01, 1), loc='upper left')
    figure.subplots_adjust(left=.1, right=.82, bottom=.12, top=.92)


def main():
    parser = argparse.ArgumentParser(add_help=False, description='Построение графиков по CSV-файлу CHEMAR.')
    parser.add_argument('csv_file', nargs='?', default='output/results.csv', metavar='файл_CSV',
                        help='Путь к CSV с результатами расчёта')
    parser.add_argument('--species', nargs='+', metavar='вещество', help='Имена веществ из заголовка CSV')
    parser.add_argument('-h', '--help', action='help', help='Показать справку и выйти')
    parser._positionals.title = 'Позиционные аргументы'
    parser._optionals.title = 'Параметры'
    args = parser.parse_args()
    path = Path(args.csv_file)
    if not path.is_file():
        parser.error(f'Файл не найден: {path}')
    with path.open(newline='', encoding='utf-8') as stream:
        reader = csv.DictReader(stream)
        headers = reader.fieldnames or []
        required = {'time_s', 'temperature_K', 'pressure_Pa', 'density_kg_m3'}
        if not required.issubset(headers) or len(headers) != len(set(headers)):
            parser.error('Отсутствуют обязательные столбцы или повторяются заголовки')
        species = {h[2:-7]: h for h in headers if h.startswith('c_') and h.endswith('_mol_m3')}
        if not species:
            parser.error('Столбцы концентраций не найдены')
        diagnostics = {
            'specific_enthalpy_J_kg', 'specific_internal_energy_J_kg',
            'specific_entropy_J_kg_K', 'mean_molar_mass_kg_mol',
        }
        diagnostics.update(f'mole_fraction_{name}' for name in species)
        diagnostics.update(f'mass_fraction_{name}' for name in species)
        if not diagnostics.issubset(headers):
            parser.error('В CSV нет новых термодинамических столбцов; перезапустите chemar')
        reaction_count = 0
        while f'reaction_{reaction_count + 1}_net_mol_m3_s' in headers:
            reaction_count += 1
        if reaction_count == 0:
            parser.error('В CSV отсутствуют скорости химических реакций')
        for reaction in range(1, reaction_count + 1):
            for direction in ('forward', 'reverse', 'net'):
                name = f'reaction_{reaction}_{direction}_mol_m3_s'
                if name not in headers:
                    parser.error(f'В CSV отсутствует столбец {name}')
        selected = args.species or list(species)
        if len(selected) != len(set(selected)) or any(s not in species for s in selected):
            parser.error(f'Выберите разные вещества из списка: {", ".join(species)}')
        columns = {h: [] for h in headers}
        try:
            for row in reader:
                if None in row or any(row[h] is None for h in headers):
                    raise ValueError('Неверное число столбцов')
                values = {h: float(row[h]) for h in headers}
                if not all(math.isfinite(v) for v in values.values()):
                    raise ValueError('В CSV обнаружены NaN или бесконечность')
                if columns['time_s'] and values['time_s'] <= columns['time_s'][-1]:
                    raise ValueError('Время должно строго возрастать')
                if any(values[h] < 0 for h in species.values()):
                    raise ValueError('Отрицательная концентрация')
                for h in headers:
                    columns[h].append(values[h])
        except (ValueError, TypeError) as error:
            parser.error(f'Некорректные данные CSV в строке {reader.line_num}: {error}')
    if not columns['time_s']:
        parser.error('Пустой CSV')
    time, time_label = time_scale(columns['time_s'])
    plt.style.use('seaborn-v0_8-whitegrid')
    plt.rcParams.update({
        'axes.spines.top': False,
        'axes.spines.right': False,
        'font.size': 10,
        'lines.linewidth': 2,
    })
    colors = {
        'H2': '#1f77b4', 'O2': '#d62728', 'H2O': '#2ca02c', 'N2': '#7f7f7f',
        'H': '#e377c2', 'O': '#9467bd', 'OH': '#ff7f0e',
        'HO2': '#8c564b', 'H2O2': '#17becf',
    }

    scalar_figure(time, columns['temperature_K'], time_label,
                  'Температура', 'Температура, К', '#d62728')
    scalar_figure(time, [value / 1000 for value in columns['pressure_Pa']], time_label,
                  'Давление', 'Давление, кПа', '#1f77b4')

    figure, axis = plt.subplots(figsize=(9, 6))
    for name in selected:
        axis.plot(time, columns[species[name]], color=colors.get(name), label=name)
    axis.set_xlabel(time_label)
    axis.set_ylabel('Концентрация, моль/м³')
    axis.set_title('Все вещества')
    finish_axis(axis)
    axis.legend(ncol=3, frameon=False)
    figure.tight_layout()

    column_count = min(3, len(selected))
    row_count = (len(selected) + column_count - 1) // column_count
    figure, axes = plt.subplots(row_count, column_count,
                                figsize=(4 * column_count, 2.8 * row_count),
                                sharex=True, squeeze=False)
    axis_list = list(axes.flat)
    for index, name in enumerate(selected):
        axis = axis_list[index]
        axis.plot(time, columns[species[name]], color=colors.get(name))
        axis.set_title(name, loc='left')
        axis.ticklabel_format(axis='y', style='sci', scilimits=(-3, 3), useMathText=True)
        finish_axis(axis)
        if index // column_count == row_count - 1:
            axis.set_xlabel(time_label)
        if index % column_count == 0:
            axis.set_ylabel('моль/м³')
    for index in range(len(selected), len(axis_list)):
        axis_list[index].set_visible(False)
    figure.suptitle('Каждое вещество отдельно')
    figure.tight_layout(rect=(0, 0, 1, .97))

    scalar_figure(time, [value / 1.0e6 for value in columns['specific_enthalpy_J_kg']],
                  time_label, 'Удельная энтальпия h(t)', 'h, МДж/кг', '#bc5090')
    scalar_figure(time, [value / 1.0e6 for value in columns['specific_internal_energy_J_kg']],
                  time_label, 'Удельная внутренняя энергия e(t)', 'e, МДж/кг',
                  '#58508d', .001)
    scalar_figure(time, [value / 1000 for value in columns['specific_entropy_J_kg_K']],
                  time_label, 'Удельная энтропия s(t)', 's, кДж/(кг·К)', '#ff6361')
    scalar_figure(time, [value * 1000 for value in columns['mean_molar_mass_kg_mol']],
                  time_label, 'Средняя молярная масса', 'M, г/моль', '#ffa600')

    reaction_rate_figure(time, columns, reaction_count, 'net', time_label,
                         'Результирующие скорости реакций', False)
    reaction_rate_figure(time, columns, reaction_count, 'forward', time_label,
                         'Прямые скорости реакций', True)
    reaction_rate_figure(time, columns, reaction_count, 'reverse', time_label,
                         'Обратные скорости реакций', True)

    fraction_figure(time, columns, selected, colors, 'mass_fraction_', time_label,
                    'Массовые доли', 'Массовая доля')
    fraction_figure(time, columns, selected, colors, 'mole_fraction_', time_label,
                    'Мольные доли', 'Мольная доля')
    plt.show()


if __name__ == '__main__':
    main()
