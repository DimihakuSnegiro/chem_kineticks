#include "solver.hpp"

#include <yaml-cpp/yaml.h>

#include <cmath>
#include <limits>
#include <set>
#include <stdexcept>

static double number(const YAML::Node& n) {
    const double x = n.as<double>();

    if (!std::isfinite(x)) {
        throw std::runtime_error("Nonfinite YAML number");
    }

    return x;
}

static size_t positive_integer(const YAML::Node& n) {
    const auto v = n.as<long long>();

    if (v <= 0) {
        throw std::runtime_error("Expected positive integer");
    }

    return size_t(v);
}

Input load(const std::filesystem::path& directory) {
    Input in;
    const auto sp = YAML::LoadFile((directory / "species.yaml").string());

    if (sp["model"].as<std::string>() != "Glushko_Gurvich") {
        throw std::runtime_error("Unsupported thermo model");
    }

    std::map<std::string, int> index;

    for (const auto& node : sp["species"]) {
        Species s;
        s.name = node["name"].as<std::string>();

        if (s.name.empty()
            || s.name.find_first_of(",\r\n\"") != std::string::npos
            || index.count(s.name)) {
            throw std::runtime_error("Invalid or duplicate species name");
        }

        index[s.name] = int(in.mechanism.species.size());
        s.molar_mass = number(node["molar_mass_kg_mol"]);

        for (const auto& a : node["atoms"]) {
            const int count = a.second.as<int>();

            if (count <= 0) {
                throw std::runtime_error("Invalid atom count");
            }

            s.atoms.emplace(a.first.as<std::string>(), count);
        }

        for (const auto& p : node["thermo"]) {
            ThermoInterval interval{number(p["T_min"]), number(p["T_max"]), {}};

            if (p["coefficients"].size() != 12) {
                throw std::runtime_error("Expected 12 Gurvich coefficients");
            }

            for (size_t i = 0; i < 12; ++i) {
                interval.b[i] = number(p["coefficients"][i]);
            }

            s.thermo.push_back(interval);
        }

        in.mechanism.species.push_back(s);
    }

    auto find = [&](const std::string& name) {
        auto it = index.find(name);

        if (it == index.end()) {
            throw std::runtime_error("Unknown species: " + name);
        }

        return it->second;
    };

    auto terms = [&](const YAML::Node& n) {
        std::vector<Term> result;
        std::set<int> seen;

        for (const auto& item : n) {
            const int i = find(item.first.as<std::string>());
            const int v = item.second.as<int>();

            if (v <= 0 || !seen.insert(i).second) {
                throw std::runtime_error("Invalid stoichiometric term");
            }

            result.push_back({i, v});
        }

        if (result.empty()) {
            throw std::runtime_error("Empty reaction side");
        }

        return result;
    };

    auto arr = [](const YAML::Node& n) {
        Arrhenius a{
            number(n["A"]),
            number(n["n"]),
            number(n["activation_temperature_K"])
        };

        if (a.A < 0) {
            throw std::runtime_error("Negative Arrhenius prefactor");
        }

        return a;
    };

    const auto reactions = YAML::LoadFile((directory / "reactions.yaml").string());

    if (reactions["units"].as<std::string>() != "SI_mol_m_s_K") {
        throw std::runtime_error("Expected SI rate units");
    }

    for (const auto& node : reactions["reactions"]) {
        Reaction r;
        r.reactants = terms(node["reactants"]);
        r.products = terms(node["products"]);
        r.forward = arr(node["forward"]);

        if (node["reverse"]) {
            r.reverse = arr(node["reverse"]);
        }

        r.reversible = node["reversible"].as<bool>();

        if (const auto body = node["third_body"]) {
            const double base = number(body["default_efficiency"]);

            if (base < 0) {
                throw std::runtime_error("Negative collider efficiency");
            }

            r.efficiencies.assign(index.size(), base);

            for (const auto& e : body["efficiencies"]) {
                const double v = number(e.second);

                if (v < 0) {
                    throw std::runtime_error("Negative collider efficiency");
                }

                r.efficiencies[find(e.first.as<std::string>())] = v;
            }
        }

        in.mechanism.reactions.push_back(r);
    }

    prepare(in.mechanism);

    const auto init = YAML::LoadFile((directory / "initial.yaml").string());

    if (init["mode"].as<std::string>() != "constant_density_energy") {
        throw std::runtime_error("Only IPR=3 is implemented");
    }

    in.initial.temperature = number(init["temperature_K"]);
    in.initial.density = number(init["density_kg_m3"]);
    in.initial.c.assign(index.size(), 0);

    double sum = 0;
    double mean_mass = 0;
    std::set<int> composition;

    for (const auto& x : init["mole_fractions"]) {
        const int i = find(x.first.as<std::string>());
        const double v = number(x.second);

        if (v < 0 || !composition.insert(i).second) {
            throw std::runtime_error("Invalid composition");
        }

        in.initial.c[i] = v;
        sum += v;
        mean_mass += v * in.mechanism.species[i].molar_mass;
    }

    if (!(sum > 0 && mean_mass > 0 && in.initial.density > 0)) {
        throw std::runtime_error("Invalid initial state");
    }

    // Normalization cancels: c_i = rho*x_i / sum(x_k*M_k).
    for (double& c : in.initial.c) {
        c *= in.initial.density / mean_mass;
    }

    if (in.initial.temperature < in.mechanism.Tmin
        || in.initial.temperature > in.mechanism.Tmax) {
        throw std::runtime_error("Initial temperature outside common thermo range");
    }

    in.end_time = number(init["time_end_s"]);
    const auto mode = init["reverse_rates"].as<std::string>();

    if (mode != "thermodynamic" && mode != "explicit") {
        throw std::runtime_error("Invalid reverse_rates");
    }

    in.explicit_reverse = mode == "explicit";

    if (in.explicit_reverse) {
        for (const auto& r : in.mechanism.reactions) {
            if (r.reversible && !r.reverse) {
                throw std::runtime_error("Missing reverse Arrhenius parameters");
            }
        }
    }

    const auto solver = init["solver"];
    auto& s = in.solver;
    s.theta = number(solver["theta"]);
    s.initial_step = number(solver["initial_step_s"]);
    s.min_step = number(solver["min_step_s"]);
    s.max_step = number(solver["max_step_s"]);
    s.rtol = number(solver["relative_tolerance"]);
    s.atol = number(solver["absolute_tolerance_mol_m3"]);
    s.nonlinear_tolerance = number(solver["nonlinear_tolerance"]);

    const size_t iterations = positive_integer(solver["max_newton_iterations"]);

    if (iterations > 1000) {
        throw std::runtime_error("Too many Newton iterations");
    }

    s.max_newton_iterations = int(iterations);
    s.max_steps = positive_integer(solver["max_steps"]);

    if (!(s.theta >= 0.5 && s.theta <= 1
          && s.min_step > 0 && s.max_step >= s.min_step
          && s.initial_step >= s.min_step && s.initial_step <= s.max_step
          && s.atol > 0 && s.rtol > 0
          && s.nonlinear_tolerance > 0 && s.nonlinear_tolerance < 0.01
          && in.end_time > 0)) {
        throw std::runtime_error("Invalid integration settings");
    }

    in.csv = init["output"]["csv"].as<std::string>();

    if (in.csv.empty()) {
        throw std::runtime_error("Empty CSV path");
    }

    in.save_every = positive_integer(init["output"]["save_every"]);
    return in;
}
