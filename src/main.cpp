#include "solver.hpp"

#include <iomanip>
#include <iostream>

int main(int argc, char** argv) {
    try {
        if (argc > 2) {
            std::cerr << "Использование: chemar [каталог_данных]\n";
            return 2;
        }

        std::filesystem::path directory = "data";

        if (argc == 2) {
            directory = argv[1];
        }

        const auto input = load(directory);
        const auto result = integrate(input);

        std::cout << std::setprecision(12)
                  << "CSV: " << input.csv.string()
                  << "\nПринято шагов: " << result.accepted
                  << ", отклонено: " << result.rejected
                  << "\nКонечное время [с]: " << result.final.time
                  << "\nКонечная температура [К]: " << result.final.temperature
                  << "\nКонечное давление [Па]: " << pressure(result.final)
                  << "\nМаксимальные относительные ошибки массы / энергии / числа атомов: " << result.max_mass_error
                  << " / " << result.max_energy_error
                  << " / " << result.max_atom_error << '\n';
    } catch (const std::exception& e) {
        std::cerr << "CHEMAR: " << e.what() << '\n';
        return 1;
    }
}
