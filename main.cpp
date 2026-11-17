#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <nlohmann/json.hpp>

namespace chemar {
namespace fs = std::filesystem;
using Vector = std::vector<double>;

constexpr double gas_constant = 8.3144; // Дж/(моль·К)
constexpr double standard_pressure = 101325.0; // Па
constexpr int third_body = -1;

bool positive(double value) {
    return std::isfinite(value) && value > 0.0;
}


struct MolarThermoProperties {
    double enthalpy, heat_capacity, entropy;
};

struct EnthalpyAndHeatCapacity {
    double enthalpy;
    double heat_capacity;
};

struct ResidualAndSlope {
    double residual;
    double slope;
};

struct ThermodynamicPolynomial {
    double min_temperature, max_temperature, reference_enthalpy;
    std::array<double, 6> enthalpy_coeffs, entropy_coeffs;
    MolarThermoProperties evaluate(double temperature) const;
};

struct Species {
    int id;
    std::string name;
    double molar_mass; // Молярная масса, кг/моль.
    // Позиции соответствуют именам из массива elements в species.json.
    std::vector<double> element_counts;
    std::vector<ThermodynamicPolynomial> thermo_polynomials;

    const ThermodynamicPolynomial& polynomial_for_temperature(double temperature) const;
};

struct ArrheniusParameters {
    double prefactor, temperature_exponent, activation_temperature;

    double log_rate_constant(double temperature) const {
        if (prefactor == 0.0) {
            return -std::numeric_limits<double>::infinity();
        }
        // Ниже 300 К используем линейную экстраполяцию исходного CHEMAR.
        const double evaluation_temperature = std::max(temperature, 300.0);
        return std::log(prefactor) + temperature_exponent*std::log(evaluation_temperature)
               - activation_temperature/evaluation_temperature
               + (std::log(temperature) - std::log(evaluation_temperature));
    }
};

struct ReactionTerm {
    int species_index; // third_body для универсального третьего тела M.
    int stoichiometric_coefficient;
    int kinetic_order;
};

struct StoichiometricChange {
    std::size_t species_index;
    int coefficient;
};

struct Reaction {
    bool reversible;
    ArrheniusParameters forward_parameters;
    std::optional<ArrheniusParameters> reverse_parameters;
    std::vector<ReactionTerm> reactants;
    std::vector<ReactionTerm> products;
    std::vector<StoichiometricChange> net_stoichiometry;
    int mole_change = 0; // Изменение общего числа молей на моль реакции.
};

struct RateConstants {
    double forward;
    double reverse;
};

enum class CalculationMode {
    fixed_density_temperature = 1,
    fixed_pressure_temperature,
    fixed_density_internal_energy,
    fixed_pressure_enthalpy,
    post_shock_flow
};

bool uses_input_density(CalculationMode mode) {
    return mode == CalculationMode::fixed_density_temperature ||
           mode == CalculationMode::fixed_density_internal_energy ||
           mode == CalculationMode::post_shock_flow;
}

struct InitialConditions {
    double velocity, density, pressure, temperature;
    CalculationMode mode;
    int output_intervals;
    double output_interval; // Секунды, а для режима потока — метры.
    bool use_tabulated_reverse_rates;
    Vector relative_molar_amounts;
};

struct SimulationInput {
    std::vector<Species> species;
    InitialConditions conditions;
    std::vector<Reaction> reactions;
};

using Json = nlohmann::json;

[[noreturn]] void json_error(const fs::path& file, const std::string& field,
                             const std::string& message) {
    throw std::runtime_error(file.filename().u8string() + ": " + field + ": " + message);
}

Json load_json_file(const fs::path& file) {
    std::ifstream input(file, std::ios::binary);
    if (!input) {
        json_error(file, "$", "cannot open file");
    }
    try {
        Json document;
        input >> document;
        return document;
    } catch (const nlohmann::json::parse_error& error) {
        json_error(file, "$", std::string("invalid JSON: ") + error.what());
    }
}

const Json* require_field(const Json& object, const std::string& name,
                          const fs::path& file, const std::string& object_path) {
    if (!object.is_object()) {
        json_error(file, object_path, "expected object");
    }
    if (!object.contains(name)) {
        json_error(file, object_path + "." + name, "missing required field");
    }
    return &object.at(name);
}

const Json& require_object(const Json& value, const fs::path& file, const std::string& path) {
    if (!value.is_object()) {
        json_error(file, path, "expected object");
    }
    return value;
}

const Json& require_array(const Json& value, const fs::path& file, const std::string& path) {
    if (!value.is_array()) {
        json_error(file, path, "expected array");
    }
    return value;
}

std::string require_string(const Json& value, const fs::path& file, const std::string& path) {
    if (!value.is_string()) {
        json_error(file, path, "expected string");
    }
    return value.get<std::string>();
}

bool require_boolean(const Json& value, const fs::path& file, const std::string& path) {
    if (!value.is_boolean()) {
        json_error(file, path, "expected boolean");
    }
    return value.get<bool>();
}

double require_number(const Json& value, const fs::path& file, const std::string& path) {
    if (!value.is_number()) {
        json_error(file, path, "expected number");
    }
    double result = value.get<double>();
    if (!std::isfinite(result)) {
        json_error(file, path, "expected finite number");
    }
    return result;
}

int require_integer(const Json& value, const fs::path& file, const std::string& path) {
    if (!value.is_number_integer()) {
        json_error(file, path, "expected integer");
    }
    try {
        long long result = value.get<long long>();
        if (result < std::numeric_limits<int>::min() || result > std::numeric_limits<int>::max()) {
            json_error(file, path, "integer is outside the supported range");
        }
        return static_cast<int>(result);
    } catch (const nlohmann::json::exception&) {
        json_error(file, path, "integer is outside the supported range");
    }
}

std::array<double, 6> read_coefficients(const Json& value, const fs::path& file,
                                        const std::string& path) {
    const Json& values = require_array(value, file, path);
    if (values.size() != 6) {
        json_error(file, path, "expected exactly 6 coefficients");
    }
    std::array<double, 6> result{};
    for (std::size_t i = 0; i < result.size(); ++i) {
        result[i] = require_number(values[i], file, path + "[" + std::to_string(i) + "]");
    }
    return result;
}

CalculationMode read_calculation_mode(const std::string& name, const fs::path& file,
                                      const std::string& path) {
    if (name == "fixed_density_temperature") return CalculationMode::fixed_density_temperature;
    if (name == "fixed_pressure_temperature") return CalculationMode::fixed_pressure_temperature;
    if (name == "fixed_density_internal_energy") return CalculationMode::fixed_density_internal_energy;
    if (name == "fixed_pressure_enthalpy") return CalculationMode::fixed_pressure_enthalpy;
    if (name == "post_shock_flow") return CalculationMode::post_shock_flow;
    json_error(file, path, "unknown calculation mode '" + name + "'");
}

ArrheniusParameters read_arrhenius(const Json& value, const fs::path& file,
                                   const std::string& path) {
    const Json& object = require_object(value, file, path);
    ArrheniusParameters parameters{
        require_number(*require_field(object, "prefactor", file, path), file, path + ".prefactor"),
        require_number(*require_field(object, "temperature_exponent", file, path), file,
                       path + ".temperature_exponent"),
        require_number(*require_field(object, "activation_temperature_k", file, path), file,
                       path + ".activation_temperature_k")
    };
    if (parameters.prefactor < 0.0) {
        json_error(file, path + ".prefactor", "must be nonnegative");
    }
    return parameters;
}

std::vector<Species> read_species_json(const fs::path& file,
                                       std::vector<std::string>& element_names) {
    Json document = load_json_file(file);
    require_object(document, file, "$");

    const Json& elements = *require_field(document, "elements", file, "$");
    require_array(elements, file, "$.elements");
    if (elements.empty()) {
        json_error(file, "$.elements", "must contain at least one element name");
    }
    std::map<std::string, std::size_t> element_index;
    for (std::size_t i = 0; i < elements.size(); ++i) {
        std::string path = "$.elements[" + std::to_string(i) + "]";
        std::string name = require_string(elements[i], file, path);
        if (name.empty()) {
            json_error(file, path, "element name must not be empty");
        }
        if (!element_index.emplace(name, i).second) {
            json_error(file, path, "duplicate element name '" + name + "'");
        }
        element_names.push_back(name);
    }

    const Json& species_array = *require_field(document, "species", file, "$");
    require_array(species_array, file, "$.species");
    if (species_array.empty()) {
        json_error(file, "$.species", "must contain at least one species");
    }

    std::vector<Species> species;
    std::map<int, std::string> used_ids;
    std::map<std::string, int> used_names;
    for (std::size_t i = 0; i < species_array.size(); ++i) {
        std::string path = "$.species[" + std::to_string(i) + "]";
        const Json& object = require_object(species_array[i], file, path);
        Species entry{};
        entry.id = require_integer(*require_field(object, "id", file, path), file, path + ".id");
        entry.name = require_string(*require_field(object, "name", file, path), file, path + ".name");
        entry.molar_mass = require_number(*require_field(object, "molar_mass_kg_mol", file, path),
                                          file, path + ".molar_mass_kg_mol");
        if (entry.id <= 0) json_error(file, path + ".id", "must be positive");
        if (entry.name.empty()) json_error(file, path + ".name", "must not be empty");
        if (entry.name == "M") json_error(file, path + ".name", "'M' is reserved for a third body");
        if (!positive(entry.molar_mass)) {
            json_error(file, path + ".molar_mass_kg_mol", "must be positive and finite");
        }
        if (!used_ids.emplace(entry.id, entry.name).second) {
            json_error(file, path + ".id", "duplicate species id " + std::to_string(entry.id));
        }
        if (!used_names.emplace(entry.name, entry.id).second) {
            json_error(file, path + ".name", "duplicate species name '" + entry.name + "'");
        }

        entry.element_counts.assign(element_names.size(), 0.0);
        const Json& counts = *require_field(object, "elements", file, path);
        require_object(counts, file, path + ".elements");
        for (Json::const_iterator count = counts.begin(); count != counts.end(); ++count) {
            std::string count_path = path + ".elements." + count.key();
            std::map<std::string, std::size_t>::const_iterator position = element_index.find(count.key());
            if (position == element_index.end()) {
                json_error(file, count_path, "element is not listed in $.elements");
            }
            double number = require_number(count.value(), file, count_path);
            if (number < 0.0) json_error(file, count_path, "must be nonnegative");
            entry.element_counts[position->second] = number;
        }

        const Json& intervals = *require_field(object, "thermodynamics", file, path);
        require_array(intervals, file, path + ".thermodynamics");
        if (intervals.empty()) {
            json_error(file, path + ".thermodynamics", "must contain at least one interval");
        }
        for (std::size_t j = 0; j < intervals.size(); ++j) {
            std::string interval_path = path + ".thermodynamics[" + std::to_string(j) + "]";
            const Json& interval = require_object(intervals[j], file, interval_path);
            ThermodynamicPolynomial polynomial{
                require_number(*require_field(interval, "minimum_temperature_k", file, interval_path),
                               file, interval_path + ".minimum_temperature_k"),
                require_number(*require_field(interval, "maximum_temperature_k", file, interval_path),
                               file, interval_path + ".maximum_temperature_k"),
                require_number(*require_field(interval, "reference_enthalpy_j_mol", file, interval_path),
                               file, interval_path + ".reference_enthalpy_j_mol"),
                read_coefficients(*require_field(interval, "enthalpy_coefficients", file, interval_path),
                                  file, interval_path + ".enthalpy_coefficients"),
                read_coefficients(*require_field(interval, "entropy_coefficients", file, interval_path),
                                  file, interval_path + ".entropy_coefficients")
            };
            if (!positive(polynomial.min_temperature)) {
                json_error(file, interval_path + ".minimum_temperature_k", "must be positive");
            }
            if (!(polynomial.max_temperature > polynomial.min_temperature)) {
                json_error(file, interval_path + ".maximum_temperature_k",
                           "must be greater than minimum_temperature_k");
            }
            if (!entry.thermo_polynomials.empty()) {
                const ThermodynamicPolynomial& previous = entry.thermo_polynomials.back();
                if (polynomial.min_temperature < previous.min_temperature ||
                    polynomial.min_temperature > previous.max_temperature ||
                    polynomial.max_temperature <= previous.max_temperature) {
                    json_error(file, interval_path, "temperature intervals must be ordered without gaps");
                }
            }
            entry.thermo_polynomials.push_back(polynomial);
        }
        species.push_back(std::move(entry));
    }
    return species;
}

InitialConditions read_calculation_json(const fs::path& file,
                                        const std::vector<Species>& species) {
    Json document = load_json_file(file);
    require_object(document, file, "$");
    const Json& calculation = *require_field(document, "calculation", file, "$");
    require_object(calculation, file, "$.calculation");
    std::string mode_name = require_string(*require_field(calculation, "mode", file, "$.calculation"),
                                           file, "$.calculation.mode");
    InitialConditions conditions{
        require_number(*require_field(calculation, "velocity_m_s", file, "$.calculation"),
                       file, "$.calculation.velocity_m_s"),
        require_number(*require_field(calculation, "density_kg_m3", file, "$.calculation"),
                       file, "$.calculation.density_kg_m3"),
        require_number(*require_field(calculation, "pressure_pa", file, "$.calculation"),
                       file, "$.calculation.pressure_pa"),
        require_number(*require_field(calculation, "temperature_k", file, "$.calculation"),
                       file, "$.calculation.temperature_k"),
        read_calculation_mode(mode_name, file, "$.calculation.mode"),
        require_integer(*require_field(calculation, "output_intervals", file, "$.calculation"),
                        file, "$.calculation.output_intervals"),
        require_number(*require_field(calculation, "output_interval", file, "$.calculation"),
                       file, "$.calculation.output_interval"),
        require_boolean(*require_field(calculation, "use_tabulated_reverse_rates", file,
                                      "$.calculation"),
                        file, "$.calculation.use_tabulated_reverse_rates"),
        {}
    };
    if (!positive(conditions.temperature)) {
        json_error(file, "$.calculation.temperature_k", "must be positive");
    }
    if (!positive(conditions.output_interval)) {
        json_error(file, "$.calculation.output_interval", "must be positive");
    }
    if (conditions.output_intervals < 0) {
        json_error(file, "$.calculation.output_intervals", "must be nonnegative");
    }
    if (uses_input_density(conditions.mode) && !positive(conditions.density)) {
        json_error(file, "$.calculation.density_kg_m3", "must be positive for this mode");
    }
    if (!uses_input_density(conditions.mode) && !positive(conditions.pressure)) {
        json_error(file, "$.calculation.pressure_pa", "must be positive for this mode");
    }
    if (conditions.mode == CalculationMode::post_shock_flow && !positive(conditions.velocity)) {
        json_error(file, "$.calculation.velocity_m_s", "must be positive for post_shock_flow");
    }

    const Json& composition = *require_field(document, "initial_composition", file, "$");
    require_object(composition, file, "$.initial_composition");
    std::map<std::string, std::size_t> species_index;
    for (std::size_t i = 0; i < species.size(); ++i) species_index.emplace(species[i].name, i);
    for (Json::const_iterator item = composition.begin(); item != composition.end(); ++item) {
        if (species_index.find(item.key()) == species_index.end()) {
            json_error(file, "$.initial_composition." + item.key(), "unknown species");
        }
    }
    for (const Species& entry : species) {
        Json::const_iterator position = composition.find(entry.name);
        std::string path = "$.initial_composition." + entry.name;
        if (position == composition.end()) json_error(file, path, "missing initial amount");
        double amount = require_number(*position, file, path);
        if (amount < 0.0) json_error(file, path, "must be nonnegative");
        conditions.relative_molar_amounts.push_back(amount);
    }
    double total_amount = std::accumulate(conditions.relative_molar_amounts.begin(),
                                          conditions.relative_molar_amounts.end(), 0.0);
    if (!positive(total_amount)) {
        json_error(file, "$.initial_composition", "sum of amounts must be positive and finite");
    }
    return conditions;
}

std::vector<ReactionTerm> read_reaction_side(const Json& value, const fs::path& file,
                                             const std::string& path,
                                             const std::map<std::string, int>& species_index) {
    const Json& names = require_array(value, file, path);
    if (names.empty()) json_error(file, path, "reaction side must not be empty");

    std::vector<ReactionTerm> terms;
    for (std::size_t i = 0; i < names.size(); ++i) {
        std::string item_path = path + "[" + std::to_string(i) + "]";
        std::string name = require_string(names[i], file, item_path);
        int index = third_body;
        if (name != "M") {
            std::map<std::string, int>::const_iterator position = species_index.find(name);
            if (position == species_index.end()) json_error(file, item_path, "unknown species '" + name + "'");
            index = position->second;
        }

        bool found = false;
        for (ReactionTerm& term : terms) {
            if (term.species_index == index) {
                ++term.stoichiometric_coefficient;
                ++term.kinetic_order;
                found = true;
                break;
            }
        }
        if (!found) terms.push_back({index, 1, 1});
    }
    return terms;
}

int third_body_order(const std::vector<ReactionTerm>& terms) {
    for (const ReactionTerm& term : terms) {
        if (term.species_index == third_body) return term.kinetic_order;
    }
    return 0;
}

std::vector<Reaction> read_reactions_json(const fs::path& file,
                                          const std::vector<Species>& species,
                                          const std::vector<std::string>& element_names,
                                          bool use_tabulated_reverse_rates) {
    Json document = load_json_file(file);
    require_object(document, file, "$");
    const Json& reaction_array = *require_field(document, "reactions", file, "$");
    require_array(reaction_array, file, "$.reactions");
    std::map<std::string, int> species_index;
    for (std::size_t i = 0; i < species.size(); ++i) {
        species_index.emplace(species[i].name, static_cast<int>(i));
    }

    std::vector<Reaction> reactions;
    for (std::size_t reaction_index = 0; reaction_index < reaction_array.size(); ++reaction_index) {
        std::string path = "$.reactions[" + std::to_string(reaction_index) + "]";
        const Json& object = require_object(reaction_array[reaction_index], file, path);
        Reaction reaction{};
        reaction.reversible = require_boolean(*require_field(object, "reversible", file, path),
                                              file, path + ".reversible");
        reaction.reactants = read_reaction_side(*require_field(object, "reactants", file, path),
                                                file, path + ".reactants", species_index);
        reaction.products = read_reaction_side(*require_field(object, "products", file, path),
                                               file, path + ".products", species_index);
        reaction.forward_parameters = read_arrhenius(*require_field(object, "forward_rate", file, path),
                                                     file, path + ".forward_rate");

        Json::const_iterator reverse_rate = object.find("reverse_rate");
        if (reverse_rate != object.end()) {
            ArrheniusParameters parameters = read_arrhenius(*reverse_rate, file, path + ".reverse_rate");
            if (reaction.reversible && use_tabulated_reverse_rates) {
                reaction.reverse_parameters = parameters;
            }
        } else if (reaction.reversible && use_tabulated_reverse_rates) {
            json_error(file, path + ".reverse_rate",
                       "required when use_tabulated_reverse_rates is true");
        }

        if (third_body_order(reaction.reactants) != third_body_order(reaction.products)) {
            json_error(file, path, "third body M must have the same multiplicity on both sides");
        }

        // увеличивает стехиометрический коэффициент и кинетический порядок на единицу.
        std::vector<int> net(species.size(), 0);
        for (const ReactionTerm& term : reaction.reactants) {
            if (term.species_index != third_body) {
                net[static_cast<std::size_t>(term.species_index)] -= term.stoichiometric_coefficient;
            }
        }
        for (const ReactionTerm& term : reaction.products) {
            if (term.species_index != third_body) {
                net[static_cast<std::size_t>(term.species_index)] += term.stoichiometric_coefficient;
            }
        }

        std::vector<double> atom_balance(element_names.size(), 0.0);
        double mass_balance = 0.0;
        double mass_scale = 0.0;
        for (std::size_t i = 0; i < net.size(); ++i) {
            if (net[i] == 0) continue;
            reaction.net_stoichiometry.push_back({i, net[i]});
            reaction.mole_change += net[i];
            mass_balance += net[i] * species[i].molar_mass;
            mass_scale += std::abs(net[i]) * species[i].molar_mass;
            for (std::size_t element = 0; element < atom_balance.size(); ++element) {
                atom_balance[element] += net[i] * species[i].element_counts[element];
            }
        }
        for (std::size_t element = 0; element < atom_balance.size(); ++element) {
            if (std::abs(atom_balance[element]) > 1e-9) {
                json_error(file, path, "reaction does not conserve element '" +
                           element_names[element] + "'");
            }
        }
        if (std::abs(mass_balance) > 1e-5 * mass_scale) {
            json_error(file, path, "reaction does not conserve mass");
        }
        reactions.push_back(std::move(reaction));
    }
    return reactions;
}

SimulationInput load_simulation_input(const fs::path& directory) {
    std::vector<std::string> element_names;
    std::vector<Species> species = read_species_json(directory / "species.json", element_names);
    InitialConditions conditions = read_calculation_json(directory / "calculation.json", species);
    std::vector<Reaction> reactions = read_reactions_json(
        directory / "reactions.json", species, element_names,
        conditions.use_tabulated_reverse_rates);
    return {std::move(species), std::move(conditions), std::move(reactions)};
}


MolarThermoProperties ThermodynamicPolynomial::evaluate(double temperature) const {
    // За пределами интервала продолжаем h и s с теплоёмкостью на его границе.
    const double boundary_temperature = std::clamp(temperature, min_temperature, max_temperature);
    const double scaled_temperature = boundary_temperature * 1e-4;
    const double inverse_square = 1.0 / (scaled_temperature * scaled_temperature);
    const double enthalpy =
        (((enthalpy_coeffs[5]*scaled_temperature + enthalpy_coeffs[4])*scaled_temperature
          + enthalpy_coeffs[3])*scaled_temperature + enthalpy_coeffs[2])*scaled_temperature
        + enthalpy_coeffs[1]/scaled_temperature + enthalpy_coeffs[0];
    const double heat_capacity =
        (((4*enthalpy_coeffs[5]*scaled_temperature + 3*enthalpy_coeffs[4])*scaled_temperature
          + 2*enthalpy_coeffs[3])*scaled_temperature + enthalpy_coeffs[2]
         - enthalpy_coeffs[1]*inverse_square)*1e-4;
    const double entropy =
        ((entropy_coeffs[5]*scaled_temperature + entropy_coeffs[4])*scaled_temperature
         + entropy_coeffs[3])*scaled_temperature + entropy_coeffs[2]*inverse_square
        + entropy_coeffs[1]*std::log(scaled_temperature) + entropy_coeffs[0];
    return {enthalpy + heat_capacity*(temperature-boundary_temperature),
            heat_capacity, entropy + heat_capacity*std::log(temperature/boundary_temperature)};
}

const ThermodynamicPolynomial& Species::polynomial_for_temperature(double temperature) const {
    const ThermodynamicPolynomial* selected = &thermo_polynomials.front();
    for (const ThermodynamicPolynomial& polynomial : thermo_polynomials) {
        if (temperature > polynomial.min_temperature) {
            selected = &polynomial;
        } else {
            break;
        }
    }
    return *selected;
}


struct GasState {
    double velocity, density, pressure, temperature, sound_speed;
    // Энтальпия, внутренняя энергия и энтропия отнесены к килограмму смеси.
    double enthalpy, internal_energy, entropy;
    double molar_mass, heat_capacity_ratio;
};

struct MixtureEvaluation {
    GasState gas_state;
    Vector composition_derivative;
};

class GasKineticsModel {
public:
    GasKineticsModel(const InitialConditions& conditions, const std::vector<Species>& species,
                    const std::vector<Reaction>& reactions)
        : conditions_(conditions), species_(species), reactions_(reactions),
          initial_composition_(conditions.relative_molar_amounts) {
        // Делим исходные количества на массу смеси, чтобы получить состав в моль/кг.
        double mixture_mass = 0.0;
        for (std::size_t i = 0; i < initial_composition_.size(); ++i) {
            mixture_mass += initial_composition_[i] * species_[i].molar_mass;
        }
        if (!positive(mixture_mass)) {
            throw std::runtime_error("Invalid mass of the initial mixture");
        }
        for (double& molar_amount : initial_composition_) {
            molar_amount /= mixture_mass;
        }
        const CalculationMode mode = conditions.mode;
        const double mixture_r = mixture_gas_constant(initial_composition_);
        const bool given_density = uses_input_density(mode);
        const double pressure =
            given_density ? conditions.density*mixture_r*conditions.temperature : conditions.pressure;
        const double density = given_density ? conditions.density : pressure/(mixture_r*conditions.temperature);
        const double velocity = mode == CalculationMode::post_shock_flow ? conditions.velocity : 1.0;
        initial_gas_state_ = compute_gas_properties(
            initial_composition_, velocity, density, pressure, conditions.temperature, mixture_r);
        // В режиме потока сохраняются потоки массы и импульса и полная энтальпия.
        mass_flux_ = density*velocity;
        momentum_flux_ = pressure + mass_flux_*velocity;
        stagnation_enthalpy_ = initial_gas_state_.enthalpy + 0.5*velocity*velocity;
        if (mode == CalculationMode::post_shock_flow) {
            initial_gas_state_ = post_shock_state(mixture_r);
        }
    }

    const Vector& initial_composition() const { return initial_composition_; }
    const GasState& initial_gas_state() const { return initial_gas_state_; }

    RateConstants rate_constants(const Reaction& reaction, double temperature,
                                 const Vector& potentials) const {
        const double log_forward_constant = reaction.forward_parameters.log_rate_constant(temperature);
        const double forward_constant = std::exp(log_forward_constant);
        double log_reverse_constant = -std::numeric_limits<double>::infinity();
        if (reaction.reverse_parameters) {
            log_reverse_constant = reaction.reverse_parameters->log_rate_constant(temperature);
        } else if (reaction.reversible && reaction.forward_parameters.prefactor > 0.0) {
            // Без табличных коэффициентов обратную константу получаем из равновесия.
            double gibbs_energy_change = 0.0;
            for (const StoichiometricChange& change : reaction.net_stoichiometry) {
                gibbs_energy_change += change.coefficient * potentials[change.species_index];
            }
            const double log_reverse_to_forward =
                reaction.mole_change*std::log(gas_constant*temperature/standard_pressure)
                + gibbs_energy_change/(gas_constant*temperature);
            log_reverse_constant = log_forward_constant + log_reverse_to_forward;
        }
        const double reverse_constant = std::exp(log_reverse_constant);
        if (!std::isfinite(forward_constant) || !std::isfinite(reverse_constant)) {
            throw std::domain_error("Nonfinite reaction rate constant");
        }
        return {forward_constant, reverse_constant};
    }

    Vector chemical_potentials(double temperature) const {
        Vector potentials;
        potentials.reserve(species_.size());
        for (const Species& entry : species_) {
            const MolarThermoProperties values =
                entry.polynomial_for_temperature(temperature).evaluate(temperature);
            potentials.push_back(values.enthalpy - temperature*values.entropy);
        }
        return potentials;
    }

    MixtureEvaluation evaluate(const Vector& composition, const GasState& guess) const {
        MixtureEvaluation result{state_for_composition(composition, guess), Vector(composition.size(), 0.0)};
        const GasState& state = result.gas_state;
        Vector concentrations;
        concentrations.reserve(composition.size());
        for (double molar_amount : composition) {
            concentrations.push_back(state.density*std::max(0.0, molar_amount)); // Из моль/кг получаем моль/м³.
        }
        const double total_concentration =
            std::accumulate(concentrations.begin(), concentrations.end(), 0.0);
        const Vector potentials =
            conditions_.use_tabulated_reverse_rates ? Vector{} : chemical_potentials(state.temperature);
        for (const Reaction& reaction : reactions_) {
            RateConstants constants = rate_constants(reaction, state.temperature, potentials);
            const double forward_progress = reaction_progress(
                constants.forward, reaction.reactants, concentrations, total_concentration);
            const double reverse_progress = reaction_progress(
                constants.reverse, reaction.products, concentrations, total_concentration);
            // Деление на плотность переводит скорость в моль/(кг·с).
            // В режиме потока деление на скорость газа даёт производную состава по расстоянию.
            const double net_rate = (forward_progress - reverse_progress) /
                                    (state.density * state.velocity);
            if (!std::isfinite(net_rate)) {
                throw std::domain_error("Nonfinite reaction rate");
            }
            for (const StoichiometricChange& change : reaction.net_stoichiometry) {
                result.composition_derivative[change.species_index] += change.coefficient * net_rate;
            }
        }
        return result;
    }

private:
    const InitialConditions& conditions_;
    const std::vector<Species>& species_;
    const std::vector<Reaction>& reactions_;
    Vector initial_composition_;
    GasState initial_gas_state_{};
    double mass_flux_{}, momentum_flux_{}, stagnation_enthalpy_{};

    static double reaction_progress(double rate_constant, const std::vector<ReactionTerm>& side,
                                    const Vector& concentrations, double total_concentration) {
        for (const ReactionTerm& term : side) {
            double concentration = total_concentration;
            if (term.species_index != third_body) {
                concentration = concentrations[static_cast<std::size_t>(term.species_index)];
            }
            // Повторное умножение сохраняет исходный закон скорости и порядок округления.
            for (int order = 0; order < term.kinetic_order; ++order) {
                rate_constant *= concentration;
            }
        }
        return rate_constant;
    }


    // Восстанавливаем состояние газа из состава и величин, заданных режимом расчёта.
    GasState state_for_composition(const Vector& composition, const GasState& guess) const {
        const double mixture_r = mixture_gas_constant(composition);
        double pressure = initial_gas_state_.pressure;
        double density = initial_gas_state_.density;
        double temperature = initial_gas_state_.temperature;
        double velocity = initial_gas_state_.velocity;
        switch (conditions_.mode) {
        case CalculationMode::fixed_density_temperature:
            pressure = density*mixture_r*temperature;
            break;
        case CalculationMode::fixed_pressure_temperature:
            density = pressure/(mixture_r*temperature);
            break;
        case CalculationMode::fixed_density_internal_energy:
            temperature = solve_temperature(
                composition, guess.temperature, mixture_r, initial_gas_state_.internal_energy);
            pressure = density*mixture_r*temperature;
            break;
        case CalculationMode::fixed_pressure_enthalpy:
            temperature = solve_temperature(composition, guess.temperature, 0.0, initial_gas_state_.enthalpy);
            density = pressure/(mixture_r*temperature);
            break;
        case CalculationMode::post_shock_flow:
            velocity = solve_flow_velocity(composition, mixture_r, guess.velocity);
            pressure = momentum_flux_-mass_flux_*velocity;
            density = mass_flux_/velocity;
            temperature = pressure/(density*mixture_r);
            break;
        }
        return compute_gas_properties(composition, velocity, density, pressure, temperature, mixture_r);
    }

    static double mixture_gas_constant(const Vector& composition) {
        const double value = gas_constant * std::accumulate(composition.begin(), composition.end(), 0.0);
        if (!positive(value)) {
            throw std::domain_error("Nonpositive total molar content");
        }
        return value;
    }

    EnthalpyAndHeatCapacity enthalpy_and_heat_capacity(const Vector& composition,
                                                       double temperature) const {
        double enthalpy = 0.0;
        double heat_capacity = 0.0;
        for (std::size_t i = 0; i < composition.size(); ++i) {
            const MolarThermoProperties values =
                species_[i].polynomial_for_temperature(temperature).evaluate(temperature);
            enthalpy += composition[i]*values.enthalpy;
            heat_capacity += composition[i]*values.heat_capacity;
        }
        return {enthalpy, heat_capacity};
    }

    GasState compute_gas_properties(const Vector& composition, double velocity, double density,
                                    double pressure, double temperature, double mixture_r) const {
        if (!positive(temperature) || !positive(density) || !positive(pressure) || !positive(velocity)) {
            throw std::domain_error("Invalid temperature, density, pressure or velocity");
        }
        double enthalpy = 0.0;
        double heat_capacity = 0.0;
        double entropy = 0.0;
        for (std::size_t i = 0; i < composition.size(); ++i) {
            const MolarThermoProperties values =
                species_[i].polynomial_for_temperature(temperature).evaluate(temperature);
            enthalpy += composition[i]*values.enthalpy;
            heat_capacity += composition[i]*values.heat_capacity;
            if (composition[i] > 0.0) {
                const double log_relative_partial_pressure = std::log(pressure/standard_pressure)
                    + std::log(composition[i]) + std::log(gas_constant/mixture_r);
                entropy += composition[i]*(values.entropy - gas_constant*log_relative_partial_pressure);
            }
        }
        if (!positive(heat_capacity-mixture_r) || !std::isfinite(enthalpy) || !std::isfinite(entropy)) {
            throw std::domain_error("Invalid mixture thermodynamic properties");
        }
        const double heat_capacity_ratio = heat_capacity/(heat_capacity-mixture_r);
        const double sound_speed = std::sqrt(heat_capacity_ratio*pressure/density);
        if (!positive(sound_speed)) {
            throw std::domain_error("Invalid speed of sound");
        }
        return {velocity, density, pressure, temperature, sound_speed,
                enthalpy, enthalpy-pressure/density, entropy, gas_constant/mixture_r, heat_capacity_ratio};
    }

    GasState post_shock_state(double mixture_r) const {
        const GasState& upstream = initial_gas_state_;
        if (upstream.velocity <= upstream.sound_speed) {
            throw std::domain_error("IPR=5 requires a supersonic upstream state (U > speed of sound)");
        }
        const double limit = momentum_flux_/mass_flux_;
        const auto temperature_at = [&](double velocity) { return velocity*(limit-velocity)/mixture_r; };

        // Исключаем корень, соответствующий состоянию до скачка, и ищем состояние после него.
        double low = 0.0;
        double high = upstream.velocity;
        for (int i = 0; i < 60; ++i) {
            const double velocity = low + (high-low)/2.0;
            if (velocity == low || velocity == high) {
                break;
            }
            const double enthalpy = enthalpy_and_heat_capacity(
                initial_composition_, temperature_at(velocity)).enthalpy;
            const double residual = (enthalpy-upstream.enthalpy)/(velocity-upstream.velocity)
                                    + 0.5*(velocity+upstream.velocity);
            if (residual > 0.0) {
                low = velocity;
            } else {
                high = velocity;
            }
        }
        const double velocity = low + (high-low)/2.0;
        const GasState downstream = compute_gas_properties(
            initial_composition_, velocity, mass_flux_/velocity, momentum_flux_-mass_flux_*velocity,
            temperature_at(velocity), mixture_r);
        const double energy_error = downstream.enthalpy + 0.5*velocity*velocity - stagnation_enthalpy_;
        if (velocity >= downstream.sound_speed || velocity >= upstream.velocity ||
            std::abs(energy_error) > 2e-12*std::max(1e4, std::abs(stagnation_enthalpy_))) {
            throw std::domain_error("Could not determine the post-shock state from the input data");
        }
        return downstream;
    }

    double solve_temperature(const Vector& composition, double guess, double mixture_r, double target) const {
        double temperature = guess;
        const auto residual = [&](double trial_temperature) -> ResidualAndSlope {
            // В режиме постоянной энергии подбираем температуру из e = h - mixture_r*T.
            // При mixture_r = 0 уравнение задаёт постоянную энтальпию.
            const EnthalpyAndHeatCapacity properties =
                enthalpy_and_heat_capacity(composition, trial_temperature);
            return {std::fma(-mixture_r, trial_temperature, properties.enthalpy)-target,
                    properties.heat_capacity-mixture_r};
        };
        for (int iteration = 0; iteration < 100; ++iteration) {
            const ResidualAndSlope current = residual(temperature);
            const double error = current.residual;
            const double slope = current.slope;
            // Точное решение по температуре уменьшает численный шум в скоростях реакций.
            const double scale = std::max({1e4, std::abs(target), std::abs(temperature*slope), mixture_r*temperature});
            const double tolerance = 8*std::numeric_limits<double>::epsilon()*scale;
            if (std::abs(error) <= tolerance) {
                return temperature;
            }
            if (!positive(slope)) {
                break;
            }
            double correction = error/slope;
            bool improved = false;
            for (int trial = 0; trial < 24; ++trial) {
                const double next = temperature-correction;
                if (positive(next) && std::abs(residual(next).residual) < std::abs(error)) {
                    temperature = next;
                    improved = true;
                    break;
                }
                correction *= 0.5;
            }
            if (!improved) {
                break;
            }
        }
        throw std::domain_error("Temperature solve did not converge at fixed energy or enthalpy");
    }

    double solve_flow_velocity(const Vector& composition, double mixture_r, double guess) const {
        const double limit = momentum_flux_/mass_flux_;
        const auto residual = [&](double velocity) -> ResidualAndSlope {
            const double temperature = velocity*(limit-velocity)/mixture_r;
            const EnthalpyAndHeatCapacity properties =
                enthalpy_and_heat_capacity(composition, temperature);
            return {properties.enthalpy+0.5*velocity*velocity-stagnation_enthalpy_,
                    properties.heat_capacity*(limit-2*velocity)/mixture_r+velocity};
        };
        double velocity = guess;
        for (int iteration = 0; iteration < 100; ++iteration) {
            const ResidualAndSlope current = residual(velocity);
            const double error = current.residual;
            const double slope = current.slope;
            const double scale = std::max({1e4, std::abs(stagnation_enthalpy_),
                                           std::abs(velocity*slope), velocity*velocity});
            const double tolerance = 8*std::numeric_limits<double>::epsilon()*scale;
            if (std::abs(error) <= tolerance) {
                return velocity;
            }
            if (!std::isfinite(slope) || slope == 0.0) {
                break;
            }
            double correction = error/slope;
            bool improved = false;
            for (int trial = 0; trial < 24; ++trial) {
                const double next = velocity-correction;
                if (next > 0.0 && next < limit) {
                    const ResidualAndSlope next_residual = residual(next);
                    // Сохраняем знак производной, чтобы оставаться на той же ветви до звуковой точки.
                    if (std::signbit(slope) == std::signbit(next_residual.slope) &&
                        std::abs(next_residual.residual) < std::abs(error)) {
                        velocity = next;
                        improved = true;
                        break;
                    }
                }
                correction *= 0.5;
            }
            if (!improved) {
                break;
            }
        }
        throw std::domain_error("Flow solve did not converge on the current velocity branch");
    }
};




bool solve_linear(Vector& matrix, Vector& rhs) {
    const std::size_t n = rhs.size();
    for (std::size_t col = 0; col < n; ++col) {
        std::size_t pivot = col;
        for (std::size_t row = col+1; row < n; ++row) {
            if (std::abs(matrix[row*n+col]) > std::abs(matrix[pivot*n+col])) {
                pivot = row;
            }
        }
        if (!positive(std::abs(matrix[pivot*n+col]))) {
            return false;
        }
        if (pivot != col) {
            for (std::size_t k = col; k < n; ++k) {
                std::swap(matrix[col*n+k], matrix[pivot*n+k]);
            }
            std::swap(rhs[col], rhs[pivot]);
        }
        for (std::size_t row = col+1; row < n; ++row) {
            const double factor = matrix[row*n+col]/matrix[col*n+col];
            for (std::size_t k = col+1; k < n; ++k) {
                matrix[row*n+k] -= factor*matrix[col*n+k];
            }
            rhs[row] -= factor*rhs[col];
        }
    }
    for (std::size_t row = n; row-- > 0;) {
        for (std::size_t k = row+1; k < n; ++k) {
            rhs[row] -= matrix[row*n+k]*rhs[k];
        }
        rhs[row] /= matrix[row*n+row];
        if (!std::isfinite(rhs[row])) {
            return false;
        }
    }
    return true;
}



struct IntegrationState {
    Vector composition;
    MixtureEvaluation evaluation;
};

struct IntegrationTolerances {
    double relative = 1e-9;
    double absolute = 1e-12; // моль/кг

    double error_scale(double first_amount, double second_amount) const {
        return absolute + relative*std::max(std::abs(first_amount), std::abs(second_amount));
    }
};

class AdaptiveIntegrator {
public:
    AdaptiveIntegrator(const GasKineticsModel& model, IntegrationTolerances tolerances)
        : model_(model), tolerances_(tolerances) {}

    // Проходим один интервал DLX и возвращаем число принятых подшагов.
    int advance_output_interval(IntegrationState& current, double output_interval) {
        const double minimum_step = 16*std::numeric_limits<double>::epsilon()*output_interval;
        double remaining_interval = output_interval;
        double step =
            suggested_step_ > 0.0 ? std::min(suggested_step_, output_interval) : output_interval;
        int accepted_substeps = 0;
        std::string last_error = "composition or error tolerance";
        for (int attempt = 0; attempt < 100000; ++attempt) {
            try {
                const std::optional<IntegrationState> full_step = take_sdirk_step(current, step);
                const std::optional<IntegrationState> first_half = take_sdirk_step(current, 0.5*step);
                std::optional<IntegrationState> two_half_steps =
                    first_half ? take_sdirk_step(*first_half, 0.5*step) : std::nullopt;
                if (full_step && two_half_steps) {
                    const double error = estimate_step_error(
                        current.composition, full_step->composition, two_half_steps->composition);
                    // Локальная погрешность пропорциональна кубу шага.
                    const double factor = error == 0.0 ? 4.0
                        : std::clamp(0.9/std::cbrt(error), 0.2, 4.0);
                    if (error <= 1.0) {
                        current = std::move(*two_half_steps);
                        accepted_substeps += 2;
                        remaining_interval -= step;
                        suggested_step_ = step*factor;
                        if (remaining_interval == 0.0) {
                            return accepted_substeps;
                        }
                        step = std::min(suggested_step_, remaining_interval);
                        // Доводим последний подшаг до точки записи без отдельного шага на остаток округления.
                        if (remaining_interval-step <= minimum_step) {
                            step = remaining_interval;
                        }
                        continue;
                    }
                    step *= factor;
                    last_error = "local error above tolerance";
                } else {
                    step *= 0.5;
                    last_error = "implicit solve or nonnegative composition";
                }
            } catch (const std::domain_error& e) {
                step *= 0.5;
                last_error = e.what();
            }
            if (step <= minimum_step) {
                throw std::runtime_error("Adaptive integration failed: " + last_error);
            }
        }
        throw std::runtime_error("Adaptive integration exceeded its substep budget");
    }

private:
    const GasKineticsModel& model_;
    IntegrationTolerances tolerances_;
    double suggested_step_ = 0.0;
    static constexpr double diagonal_ = 0.2928932188134524756; // 1 - 1/sqrt(2)

    double estimate_step_error(const Vector& start, const Vector& full, const Vector& refined) const {
        // Для метода второго порядка разность, делённая на 3, оценивает ошибку двух полушагов.
        double error = 0.0;
        for (std::size_t i = 0; i < start.size(); ++i) {
            const double scale = tolerances_.error_scale(start[i], refined[i]);
            error = std::max(error, std::abs(refined[i]-full[i])/(3.0*scale));
        }
        return error;
    }

    std::optional<IntegrationState> take_sdirk_step(const IntegrationState& start, double step) const {
        const double implicit_step = diagonal_*step;
        std::optional<IntegrationState> first = solve_stage(start.composition, implicit_step, start);
        if (!first) {
            return std::nullopt;
        }
        Vector base_composition = start.composition;
        for (std::size_t i = 0; i < base_composition.size(); ++i) {
            base_composition[i] += (1.0-diagonal_)*step*first->evaluation.composition_derivative[i];
        }
        return solve_stage(base_composition, implicit_step, *first);
    }

    std::optional<IntegrationState> solve_stage(const Vector& base_composition, double step,
                                                 const IntegrationState& guess) const {
        constexpr double stage_tolerance = 0.02; // Доля допуска интегрирования на невязку стадии.
        IntegrationState current = guess;
        const std::size_t species_count = base_composition.size();
        const auto residual = [&](const IntegrationState& trial_state, Vector& values) {
            double norm = 0.0;
            for (std::size_t i = 0; i < species_count; ++i) {
                values[i] = trial_state.composition[i]-base_composition[i]
                            - step*trial_state.evaluation.composition_derivative[i];
                if (!std::isfinite(values[i])) {
                    return std::numeric_limits<double>::infinity();
                }
                const double scale =
                    tolerances_.error_scale(trial_state.composition[i], base_composition[i]);
                norm = std::max(norm, std::abs(values[i])/scale);
            }
            return norm;
        };
        Vector correction(species_count);
        for (int iteration = 0; iteration < 12; ++iteration) {
            const double norm = residual(current, correction);
            if (norm <= stage_tolerance) {
                bool clipped = false;
                for (double& molar_amount : current.composition) {
                    // Округляем к нулю только малые отрицательные значения в пределах допуска.
                    if (molar_amount < -0.002*tolerances_.absolute) {
                        return std::nullopt;
                    }
                    if (molar_amount < 0.0) {
                        molar_amount = 0.0;
                        clipped = true;
                    }
                }
                if (clipped) {
                    current.evaluation = model_.evaluate(current.composition, current.evaluation.gas_state);
                    if (residual(current, correction) > stage_tolerance) {
                        return std::nullopt;
                    }
                }
                return current;
            }
            // Оцениваем матрицу Якоби, по очереди слегка меняя количество каждого вещества.
            Vector jacobian(species_count*species_count);
            Vector shifted = current.composition;
            for (std::size_t col = 0; col < species_count; ++col) {
                shifted[col] += 1.49e-8*std::max(1.0, std::abs(current.composition[col]));
                const double offset = shifted[col] - current.composition[col];
                const MixtureEvaluation value = model_.evaluate(shifted, current.evaluation.gas_state);
                shifted[col] = current.composition[col];
                for (std::size_t row = 0; row < species_count; ++row) {
                    jacobian[row*species_count+col] = (row == col ? 1.0 : 0.0) -
                        step*(value.composition_derivative[row]-current.evaluation.composition_derivative[row])/offset;
                }
            }
            for (double& value : correction) {
                value = -value;
            }
            if (!solve_linear(jacobian, correction)) {
                return std::nullopt;
            }
            IntegrationState trial;
            trial.composition.resize(species_count);
            Vector trial_residual(species_count);
            bool improved = false;
            for (double damping = 1.0; damping >= 1.0/128.0; damping *= 0.5) {
                for (std::size_t i = 0; i < species_count; ++i) {
                    trial.composition[i] = current.composition[i]+damping*correction[i];
                }
                try {
                    trial.evaluation = model_.evaluate(trial.composition, current.evaluation.gas_state);
                    if (residual(trial, trial_residual) < norm) {
                        current = std::move(trial);
                        improved = true;
                        break;
                    }
                } catch (const std::domain_error&) {
                    // Повторяем попытку с меньшей поправкой Ньютона.
                }
            }
            if (!improved) {
                return std::nullopt;
            }
        }
        return std::nullopt;
    }
};

struct RunOptions {
    fs::path input_directory = "data";
    fs::path output_directory = "output";
    std::optional<int> output_intervals;
    IntegrationTolerances integration_tolerances;
    bool save_diagnostics = false;
    bool show_help = false;
};

const fs::path& next_option_value(const std::vector<fs::path>& arguments, std::size_t& index,
                                  const std::string& option) {
    ++index;
    if (index == arguments.size()) {
        throw std::runtime_error("Missing value for " + option);
    }
    return arguments[index];
}

int parse_option_integer(const fs::path& value, const std::string& option) {
    std::string text = value.u8string();
    std::size_t consumed = 0;
    try {
        int result = std::stoi(text, &consumed);
        if (consumed == text.size()) return result;
    } catch (const std::exception&) {
    }
    throw std::runtime_error("Invalid integer for " + option + ": " + text);
}

double parse_option_number(const fs::path& value, const std::string& option) {
    std::string text = value.u8string();
    std::size_t consumed = 0;
    try {
        double result = std::stod(text, &consumed);
        if (consumed == text.size() && std::isfinite(result)) return result;
    } catch (const std::exception&) {
    }
    throw std::runtime_error("Invalid number for " + option + ": " + text);
}

RunOptions parse_options(const std::vector<fs::path>& arguments) {
    RunOptions options;
    for (std::size_t i = 1; i < arguments.size(); ++i) {
        const std::string key = arguments[i].u8string();
        if (key == "--help" || key == "-h") {
            options.show_help = true;
            continue;
        }
        if (key == "--diagnostics") {
            options.save_diagnostics = true;
            continue;
        }
        if (key == "--data") {
            options.input_directory = next_option_value(arguments, i, key);
        } else if (key == "--output") {
            options.output_directory = next_option_value(arguments, i, key);
        } else if (key == "--steps") {
            options.output_intervals = parse_option_integer(next_option_value(arguments, i, key), key);
        } else if (key == "--rtol") {
            options.integration_tolerances.relative =
                parse_option_number(next_option_value(arguments, i, key), key);
        } else if (key == "--atol") {
            options.integration_tolerances.absolute =
                parse_option_number(next_option_value(arguments, i, key), key);
        } else {
            throw std::runtime_error("Unknown option: " + key);
        }
    }
    if ((options.output_intervals && *options.output_intervals < 0) || !positive(options.integration_tolerances.relative) ||
        !positive(options.integration_tolerances.absolute)) {
        throw std::runtime_error("Step count must be nonnegative and tolerances must be positive");
    }
    return options;
}

std::ofstream create_output(const fs::path& path) {
    std::ofstream file(path);
    if (!file) {
        throw std::runtime_error("Cannot write " + path.u8string());
    }
    file.exceptions(std::ios::badbit | std::ios::failbit);
    return file;
}

std::string csv_field(const std::string& value) {
    std::string result = "\"";
    for (char character : value) {
        if (character == '"') result += '"';
        result += character;
    }
    result += '"';
    return result;
}

void write_csv_state(std::ostream& output, double coordinate,
                     const IntegrationState& solution, int accepted_substeps) {
    const GasState& gas = solution.evaluation.gas_state;
    output << coordinate << ',' << gas.velocity << ',' << gas.pressure << ','
           << gas.temperature << ',' << gas.density << ',' << 1.0 / gas.density << ','
           << gas.sound_speed << ',' << gas.enthalpy << ',' << gas.internal_energy << ','
           << gas.entropy << ',' << gas.molar_mass << ',' << gas.heat_capacity_ratio << ','
           << accepted_substeps << ',' << solution.composition.size();
    for (double amount : solution.composition) output << ',' << amount;
    output << '\n';
}

void write_json_string(std::ostream& output, const std::string& value) {
    output << '"';
    constexpr char hex[] = "0123456789abcdef";
    for (const unsigned char ch : value) {
        switch (ch) {
        case '"': output << "\\\""; break;
        case '\\': output << "\\\\"; break;
        case '\b': output << "\\b"; break;
        case '\f': output << "\\f"; break;
        case '\n': output << "\\n"; break;
        case '\r': output << "\\r"; break;
        case '\t': output << "\\t"; break;
        default:
            if (ch < 0x20) {
                output << "\\u00" << hex[ch >> 4] << hex[ch & 0x0f];
            } else {
                output << static_cast<char>(ch);
            }
        }
    }
    output << '"';
}

const char* calculation_mode_name(CalculationMode mode) {
    switch (mode) {
    case CalculationMode::fixed_density_temperature: return "fixed_density_temperature";
    case CalculationMode::fixed_pressure_temperature: return "fixed_pressure_temperature";
    case CalculationMode::fixed_density_internal_energy: return "fixed_density_internal_energy";
    case CalculationMode::fixed_pressure_enthalpy: return "fixed_pressure_enthalpy";
    case CalculationMode::post_shock_flow: return "post_shock_flow";
    }
    throw std::logic_error("Unknown calculation mode");
}

void write_json_state(std::ostream& output, double coordinate, CalculationMode mode,
                      const IntegrationState& solution, int accepted_substeps,
                      const std::vector<Species>& species) {
    const GasState& gas = solution.evaluation.gas_state;
    output << "    {\n      \""
           << (mode == CalculationMode::post_shock_flow ? "distance_m" : "time_s")
           << "\": " << coordinate
           << ",\n      \"velocity_m_s\": " << gas.velocity
           << ",\n      \"pressure_pa\": " << gas.pressure
           << ",\n      \"temperature_k\": " << gas.temperature
           << ",\n      \"density_kg_m3\": " << gas.density
           << ",\n      \"specific_volume_m3_kg\": " << 1.0 / gas.density
           << ",\n      \"sound_speed_m_s\": " << gas.sound_speed
           << ",\n      \"specific_enthalpy_j_kg\": " << gas.enthalpy
           << ",\n      \"specific_internal_energy_j_kg\": " << gas.internal_energy
           << ",\n      \"specific_entropy_j_kg_k\": " << gas.entropy
           << ",\n      \"molar_mass_kg_mol\": " << gas.molar_mass
           << ",\n      \"heat_capacity_ratio\": " << gas.heat_capacity_ratio
           << ",\n      \"accepted_substeps\": " << accepted_substeps
           << ",\n      \"composition_mol_kg\": {";
    for (std::size_t i = 0; i < species.size(); ++i) {
        if (i != 0) output << ',';
        output << "\n        ";
        write_json_string(output, species[i].name);
        output << ": " << solution.composition[i];
    }
    output << "\n      }\n    }";
}

void write_reaction_side(std::ostream& output, const std::vector<ReactionTerm>& terms,
                         const std::vector<Species>& species) {
    bool first = true;
    for (const ReactionTerm& term : terms) {
        for (int occurrence = 0; occurrence < term.stoichiometric_coefficient; ++occurrence) {
            if (!first) output << " + ";
            first = false;
            if (term.species_index == third_body) {
                output << 'M';
            } else {
                output << species[static_cast<std::size_t>(term.species_index)].name;
            }
        }
    }
}

void write_diagnostics(const fs::path& out, const GasKineticsModel& model,
                       const std::vector<Species>& species, const std::vector<Reaction>& reactions) {
    std::ofstream thermo_report = create_output(out / "cont_td.res");
    thermo_report << std::scientific << std::setprecision(8);
    for (const Species& entry : species) {
        const double reference_enthalpy = entry.thermo_polynomials.front().reference_enthalpy;
        const MolarThermoProperties values =
            entry.polynomial_for_temperature(298.15).evaluate(298.15);
        thermo_report << "Species: " << entry.name << "\nH298=" << values.enthalpy
                      << " CP298=" << values.heat_capacity << " S298=" << values.entropy
                      << "\nT  CP  S  H-Href  -(G-Href)/T\n";
        const auto write_thermo_row = [&](double temperature) {
            const MolarThermoProperties properties =
                entry.polynomial_for_temperature(temperature).evaluate(temperature);
            const double enthalpy_difference = properties.enthalpy-reference_enthalpy;
            thermo_report << temperature << ' ' << properties.heat_capacity << ' '
                          << properties.entropy << ' ' << enthalpy_difference << ' '
                          << properties.entropy-enthalpy_difference/temperature << '\n';
        };
        for (int j = 1; j <= 101; ++j) {
            const double temperature = 100.0*j;
            if (entry.thermo_polynomials.size() == 2) {
                const double boundary = entry.thermo_polynomials[1].min_temperature;
                if (temperature-100.0 <= boundary && boundary < temperature) {
                    write_thermo_row(std::max(boundary*0.99, boundary-0.01));
                    write_thermo_row(boundary+0.01);
                }
            }
            write_thermo_row(temperature);
        }
    }
    thermo_report.close();

    std::ofstream reaction_report = create_output(out / "cont_ck.res");
    reaction_report << std::scientific << std::setprecision(8) << "kr=" << reactions.size() << '\n';
    for (std::size_t j = 0; j < reactions.size(); ++j) {
        const Reaction& reaction = reactions[j];
        reaction_report << "Reaction " << j+1 << ": ";
        write_reaction_side(reaction_report, reaction.reactants, species);
        reaction_report << (reaction.reversible ? " <=> " : " => ");
        write_reaction_side(reaction_report, reaction.products, species);
        reaction_report << "\nT  k_forward  k_reverse\n";
        for (int i = 0; i <= 100; ++i) {
            const double temperature = 2000.0+50.0*i;
            RateConstants constants =
                model.rate_constants(reaction, temperature, model.chemical_potentials(temperature));
            reaction_report << temperature << ' ' << constants.forward << ' ' << constants.reverse << '\n';
        }
    }
    reaction_report.close();
}



int run_simulation(const RunOptions& options) {
    SimulationInput input = load_simulation_input(options.input_directory);
    const std::vector<Species>& species = input.species;
    const InitialConditions& conditions = input.conditions;
    const std::vector<Reaction>& reactions = input.reactions;
    std::cout << "Loaded " << species.size() << " species and " << reactions.size()
              << " reactions from JSON.\n";
    const int requested_intervals = options.output_intervals.value_or(conditions.output_intervals);
    if (!std::isfinite(requested_intervals*conditions.output_interval)) {
        throw std::runtime_error("Requested calculation interval is too large");
    }
    const GasKineticsModel model(conditions, species, reactions);
    IntegrationState integration_state{
        model.initial_composition(), model.evaluate(model.initial_composition(), model.initial_gas_state())
    };
    AdaptiveIntegrator integrator(model, options.integration_tolerances);
    // Открываем результаты текущего прогона с перезаписью содержимого.
    fs::create_directories(options.output_directory);
    std::ofstream results = create_output(options.output_directory / "results.json");
    std::ofstream csv = create_output(options.output_directory / "results.csv");
    if (options.save_diagnostics) {
        write_diagnostics(options.output_directory, model, species, reactions);
    }

    const bool has_spatial_coordinate = conditions.mode == CalculationMode::post_shock_flow;
    results << std::setprecision(std::numeric_limits<double>::max_digits10)
            << "{\n  \"format_version\": 1,\n  \"calculation_mode\": ";
    write_json_string(results, calculation_mode_name(conditions.mode));
    results << ",\n  \"coordinate\": {\"name\": \""
            << (has_spatial_coordinate ? "distance" : "time") << "\", \"unit\": \""
            << (has_spatial_coordinate ? "m" : "s") << "\"},\n  \"species\": [";
    for (std::size_t i = 0; i < species.size(); ++i) {
        if (i != 0) results << ", ";
        write_json_string(results, species[i].name);
    }
    results << "],\n  \"states\": [\n";
    csv << (has_spatial_coordinate ? "x" : "time")
        << ",u,p,t,rho,v,a,h,e,s,vm,gamma,substeps,km";
    for (const Species& entry : species) csv << ',' << csv_field(entry.name);
    csv << '\n' << std::setprecision(std::numeric_limits<double>::max_digits10);
    // Сохраняем начальное состояние и по одной строке после каждого завершённого интервала.
    write_json_state(results, 0.0, conditions.mode, integration_state, 0, species);
    write_csv_state(csv, 0.0, integration_state, 0);
    int completed_intervals = 0;
    std::string error;
    while (completed_intervals < requested_intervals) {
        int accepted_substeps = 0;
        try {
            accepted_substeps =
                integrator.advance_output_interval(integration_state, conditions.output_interval);
        } catch (const std::exception& e) {
            error = e.what();
            break;
        }
        ++completed_intervals;
        const double coordinate = completed_intervals*conditions.output_interval;
        results << ",\n";
        write_json_state(results, coordinate, conditions.mode, integration_state,
                         accepted_substeps, species);
        write_csv_state(csv, coordinate, integration_state, accepted_substeps);
    }
    const double last_written_coordinate = completed_intervals*conditions.output_interval;
    results << "\n  ],\n  \"completed_intervals\": " << completed_intervals
            << ",\n  \"completed\": " << (error.empty() ? "true" : "false");
    if (!error.empty()) {
        results << ",\n  \"error\": ";
        write_json_string(results, error);
    }
    results << "\n}\n";
    results.close();
    csv.close();
    if (!error.empty()) {
        std::cerr << "CHEMAR stopped at output step " << static_cast<long long>(completed_intervals)+1
                  << ". Last saved " << (has_spatial_coordinate ? "distance" : "time") << ": "
                  << last_written_coordinate << (has_spatial_coordinate ? " m. " : " s. ") << error << '\n';
        return 1;
    }
    std::cout << "Completed " << completed_intervals << " output steps. Results: "
              << (options.output_directory / "results.csv").u8string() << '\n';
    return 0;
}

}

#ifdef _WIN32
int wmain(int argc, wchar_t* argv[]) {
#else
int main(int argc, char* argv[]) {
#endif
    try {
        const std::vector<chemar::fs::path> arguments(argv, argv + argc);
        const chemar::RunOptions options = chemar::parse_options(arguments);
        if (options.show_help) {
            std::cout << "CHEMAR - nonequilibrium chemical kinetics\n"
                      << "Usage: chemar [--data DIR] [--output DIR] [--steps N]\n"
                      << "              [--rtol R] [--atol A] [--diagnostics]\n"
                      << "Defaults: --data data --output output --rtol 1e-9 --atol 1e-12\n"
                      << "Input: calculation.json, species.json and reactions.json.\n"
                      << "Results: output/results.csv and output/results.json.\n";
            return 0;
        }
        return chemar::run_simulation(options);
    } catch (const std::exception& e) {
        std::cerr << "CHEMAR: " << e.what() << '\n';
        return 1;
    }
}
