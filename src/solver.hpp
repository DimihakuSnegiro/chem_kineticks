#pragma once

#include "chemistry.hpp"

#include <filesystem>

struct Settings {
    double theta = 0.6;
    double initial_step = 1e-8;
    double min_step = 1e-18;
    double max_step = 1e-6;
    double rtol = 1e-6;
    double atol = 1e-11;
    double nonlinear_tolerance = 1e-5;
    int max_newton_iterations = 30;
    size_t max_steps = 2000000;
};

struct Input {
    Mechanism mechanism;
    State initial;
    Settings solver;
    double end_time = 0;
    bool explicit_reverse = false;
    std::filesystem::path csv;
    size_t save_every = 10;
};

struct Statistics {
    size_t accepted = 0;
    size_t rejected = 0;
    double max_mass_error = 0;
    double max_energy_error = 0;
    double max_atom_error = 0;
    State final;
};

Input load(const std::filesystem::path& directory);

// Exposed for verification of the discrete theta equation.
State implicit_step(const Mechanism& mechanism, const State& old,
                    double dt, double energy0,
                    const Settings& settings, bool explicit_reverse);

Statistics integrate(const Input& input);
