"""Преобразование старых файлов CHEMAR в три понятных JSON-файла."""

import argparse
import json
from pathlib import Path


def parse_number(token):
    return float(token.replace("D", "E").replace("d", "e"))


class LegacyReader:
    def __init__(self, path):
        self.path = path
        self.lines = path.read_text(encoding="utf-8-sig").splitlines()
        self.position = 0

    def skip_line(self):
        self.position += 1

    def read_record(self, minimum):
        values = []
        while len(values) < minimum:
            if self.position >= len(self.lines):
                raise ValueError(f"конец файла {self.path}")
            line = self.lines[self.position]
            self.position += 1
            index = 0
            while index < len(line):
                if line[index] in "!<":
                    break
                if line[index].isspace() or line[index] == ",":
                    index += 1
                    continue
                if line[index] in "'\"":
                    quote = line[index]
                    index += 1
                    end = line.index(quote, index)
                    values.append(line[index:end].strip())
                    index = end + 1
                    continue
                begin = index
                while index < len(line) and not line[index].isspace() and line[index] not in ",!<":
                    index += 1
                token = line[begin:index]
                if "*" in token and token.split("*", 1)[0].isdigit():
                    count, value = token.split("*", 1)
                    values.extend([value] * int(count))
                else:
                    values.append(token)
        return values


def read_species(path):
    reader = LegacyReader(path)
    species_count, element_count = map(int, reader.read_record(2)[:2])
    # В старом формате позиции начинаются с C, H, O, N, а число в заголовке
    # показывает количество реально используемых элементов, а не длину массива.
    legacy_slot_names = ["C", "H", "O", "N"] + [
        f"element_slot_{index + 1}" for index in range(4, 32)
    ]

    species = []
    id_to_name = {}
    for _ in range(species_count):
        fields = reader.read_record(5)
        species_id = int(fields[0])
        name = fields[1]
        element_values = [parse_number(value) for value in reader.read_record(32)]
        interval_count = int(reader.read_record(1)[0])
        intervals = []
        for _ in range(interval_count):
            row = reader.read_record(14)
            converted = [0.0] * 13
            for index in range(1, 13):
                value = parse_number(row[index - 1])
                converted[index] = value * (1000.0 if index <= 4 else 1.0)
            intervals.append({
                "minimum_temperature_k": parse_number(row[12]),
                "maximum_temperature_k": parse_number(row[13]),
                "reference_enthalpy_j_mol": converted[3] - converted[4],
                "enthalpy_coefficients": [
                    converted[3] - converted[4] - 1e4 * converted[8],
                    -2e4 * converted[7],
                    1e4 * converted[6],
                    1e4 * converted[9],
                    2e4 * converted[10],
                    3e4 * converted[11],
                ],
                "entropy_coefficients": [
                    converted[5] + converted[6],
                    converted[6],
                    -converted[7],
                    2 * converted[9],
                    3 * converted[10],
                    4 * converted[11],
                ],
            })
        species.append({
            "id": species_id,
            "name": name,
            "molar_mass_kg_mol": parse_number(fields[3]) * 1e-3,
            "elements": element_values,
            "thermodynamics": intervals,
        })
        id_to_name[species_id] = name
    active_slots = [
        index for index in range(32)
        if any(entry["elements"][index] != 0 for entry in species)
    ]
    if len(active_slots) != element_count:
        raise ValueError(
            "Количество используемых элементов не совпадает с заголовком bufermm.h2"
        )
    element_names = [legacy_slot_names[index] for index in active_slots]
    for entry in species:
        counts = entry["elements"]
        entry["elements"] = {
            legacy_slot_names[index]: counts[index]
            for index in active_slots
            if counts[index] != 0
        }
    return element_names, species, id_to_name


def read_calculation(path, species):
    reader = LegacyReader(path)
    reader.skip_line()
    row = reader.read_record(9)
    reader.skip_line()
    amounts = [parse_number(value) for value in reader.read_record(len(species))]
    modes = {
        1: "fixed_density_temperature",
        2: "fixed_pressure_temperature",
        3: "fixed_density_internal_energy",
        4: "fixed_pressure_enthalpy",
        5: "post_shock_flow",
    }
    calculation = {
        "mode": modes[int(row[4])],
        "velocity_m_s": parse_number(row[0]),
        "density_kg_m3": parse_number(row[1]),
        "pressure_pa": parse_number(row[2]),
        "temperature_k": parse_number(row[3]),
        "output_intervals": int(row[6]),
        "output_interval": parse_number(row[7]),
        "use_tabulated_reverse_rates": bool(int(row[8])),
    }
    composition = {entry["name"]: amount for entry, amount in zip(species, amounts)}
    return calculation, composition


def read_reactions(path, id_to_name):
    layouts = {
        1: (3, 3), 2: (3, 3), 3: (3, 2), 4: (3, 2), 5: (3, 1),
        6: (2, 3), 7: (2, 3), 8: (2, 2), 9: (2, 1), 10: (1, 3), 11: (1, 2),
    }
    reader = LegacyReader(path)
    reaction_count = int(reader.read_record(1)[0])

    def species_name(identifier):
        if identifier == -2:
            return "M"
        return id_to_name[identifier]

    def arrhenius(fields, start):
        return {
            "prefactor": parse_number(fields[start]),
            "temperature_exponent": parse_number(fields[start + 1]),
            "activation_temperature_k": parse_number(fields[start + 2]),
        }

    reactions = []
    for _ in range(reaction_count):
        fields = reader.read_record(11)
        reaction_type = int(fields[0])
        reactant_count, product_count = layouts[reaction_type]
        reaction = {
            "reversible": int(fields[1]) == 0,
            "reactants": [species_name(int(value)) for value in fields[2:2 + reactant_count]],
            "products": [species_name(int(value)) for value in fields[5:5 + product_count]],
            "forward_rate": arrhenius(fields, 8),
        }
        if reaction["reversible"] and len(fields) >= 14:
            reaction["reverse_rate"] = arrhenius(fields, 11)
        reactions.append(reaction)
    return reactions


def main():
    parser = argparse.ArgumentParser(description="Преобразовать старые файлы CHEMAR в JSON.")
    parser.add_argument("source", type=Path, nargs="?", default=Path("data"))
    parser.add_argument("destination", type=Path, nargs="?", default=Path("data"))
    args = parser.parse_args()

    elements, species, id_to_name = read_species(args.source / "bufermm.h2")
    calculation, composition = read_calculation(args.source / "DHIM.h2", species)
    calculation_document = {
        "calculation": calculation,
        "initial_composition": composition,
    }
    species_document = {
        "elements": elements,
        "species": species,
    }
    reactions_document = {
        "reactions": read_reactions(args.source / "buferk1", id_to_name),
    }
    args.destination.mkdir(parents=True, exist_ok=True)
    documents = {
        "calculation.json": calculation_document,
        "species.json": species_document,
        "reactions.json": reactions_document,
    }
    for filename, document in documents.items():
        path = args.destination / filename
        path.write_text(
            json.dumps(document, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
        )
        print(f"Создан {path}")


if __name__ == "__main__":
    main()
