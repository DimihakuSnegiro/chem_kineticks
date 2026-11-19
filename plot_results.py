"""Plot an existing CHEMAR CSV. Requires only matplotlib and the standard library."""
import argparse
import csv
import math
from pathlib import Path
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('csv_file', nargs='?', default='output/results.csv')
    parser.add_argument('--species', nargs='+', help='Species names from the CSV header')
    args = parser.parse_args()
    path = Path(args.csv_file)
    if not path.is_file():
        parser.error(f'File not found: {path}')
    with path.open(newline='', encoding='utf-8') as stream:
        reader = csv.DictReader(stream)
        headers = reader.fieldnames or []
        required = {'time_s', 'temperature_K', 'pressure_Pa', 'density_kg_m3'}
        if not required.issubset(headers) or len(headers) != len(set(headers)):
            parser.error('Missing required columns or duplicate headers')
        species = {h[2:-7]: h for h in headers if h.startswith('c_') and h.endswith('_mol_m3')}
        if not species:
            parser.error('No concentration columns found')
        selected = args.species or list(species)
        if len(selected) != len(set(selected)) or any(s not in species for s in selected):
            parser.error(f'Choose distinct species from: {", ".join(species)}')
        columns = {h: [] for h in headers}
        try:
            for row in reader:
                if None in row or any(row[h] is None for h in headers):
                    raise ValueError('Wrong column count')
                values = {h: float(row[h]) for h in headers}
                if not all(math.isfinite(v) for v in values.values()):
                    raise ValueError('NaN or Inf in CSV')
                if columns['time_s'] and values['time_s'] <= columns['time_s'][-1]:
                    raise ValueError('Time must increase strictly')
                if any(values[h] < 0 for h in species.values()):
                    raise ValueError('Negative concentration')
                for h in headers:
                    columns[h].append(values[h])
        except (ValueError, TypeError) as error:
            parser.error(f'Invalid CSV at line {reader.line_num}: {error}')
    if not columns['time_s']:
        parser.error('Empty CSV')
    time = columns['time_s']
    for column, label in [
        ('temperature_K', 'Temperature, K'),
        ('pressure_Pa', 'Pressure, Pa'),
    ]:
        plt.figure(figsize=(7, 4))
        plt.plot(time, columns[column], label=label.split(',')[0])
        plt.xlabel('Time, s')
        plt.ylabel(label)
        plt.grid(True, alpha=.3)
        plt.legend()
        plt.ticklabel_format(axis='x', style='sci', scilimits=(-3, 3))
        plt.tight_layout()
    plt.figure(figsize=(8, 5))
    for name in selected:
        plt.plot(time, columns[species[name]], label=name)
    plt.xlabel('Time, s')
    plt.ylabel('Concentration, mol/m³')
    plt.grid(True, alpha=.3)
    plt.legend(ncol=3)
    plt.ticklabel_format(axis='x', style='sci', scilimits=(-3, 3))
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    main()
