import argparse
import csv
import math
from pathlib import Path
import matplotlib.pyplot as plt


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
    time = columns['time_s']
    for column, label in [
        ('temperature_K', 'Температура, К'),
        ('pressure_Pa', 'Давление, Па'),
    ]:
        plt.figure(figsize=(7, 4))
        plt.plot(time, columns[column], label=label.split(',')[0])
        plt.xlabel('Время, с')
        plt.ylabel(label)
        plt.grid(True, alpha=.3)
        plt.legend()
        plt.ticklabel_format(axis='x', style='sci', scilimits=(-3, 3))
        plt.tight_layout()
    plt.figure(figsize=(8, 5))
    for name in selected:
        plt.plot(time, columns[species[name]], label=name)
    plt.xlabel('Время, с')
    plt.ylabel('Концентрация, моль/м³')
    plt.grid(True, alpha=.3)
    plt.legend(ncol=3)
    plt.ticklabel_format(axis='x', style='sci', scilimits=(-3, 3))
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    main()
