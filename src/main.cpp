#include "solver.hpp"

#include <iomanip>
#include <iostream>

int main(int argc, char** argv) {
    try {
        if (argc > 2) {
            std::cerr << "Usage: chemar [data_directory]\n";
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
                  << "\nAccepted: " << result.accepted
                  << ", rejected: " << result.rejected
                  << "\nFinal time [s]: " << result.final.time
                  << "\nFinal T [K]: " << result.final.temperature
                  << "\nFinal p [Pa]: " << pressure(result.final)
                  << "\nMax relative mass / energy / atom errors: " << result.max_mass_error
                  << " / " << result.max_energy_error
                  << " / " << result.max_atom_error << '\n';
    } catch (const std::exception& e) {
        std::cerr << "CHEMAR: " << e.what() << '\n';
        return 1;
    }
}
