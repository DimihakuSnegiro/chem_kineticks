#include "solver.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <stdexcept>

static std::vector<double> linear_solve(std::vector<std::vector<double>> a, std::vector<double> b) {
    const size_t n = b.size();

    for (size_t k = 0; k < n; ++k) {
        size_t pivot = k;

        for (size_t i = k + 1; i < n; ++i) {
            if (std::abs(a[i][k]) > std::abs(a[pivot][k])) {
                pivot = i;
            }
        }

        if (!std::isfinite(a[pivot][k]) || std::abs(a[pivot][k]) < 1e-30) {
            throw std::runtime_error("Вырожденная матрица метода Ньютона");
        }

        std::swap(a[k], a[pivot]);
        std::swap(b[k], b[pivot]);

        for (size_t i = k + 1; i < n; ++i) {
            const double factor = a[i][k] / a[k][k];

            for (size_t j = k + 1; j < n; ++j) {
                a[i][j] -= factor * a[k][j];
            }

            b[i] -= factor * b[k];
        }
    }

    std::vector<double> x(n);

    for (size_t k = n; k-- > 0;) {
        double v = b[k];

        for (size_t j = k + 1; j < n; ++j) {
            v -= a[k][j] * x[j];
        }

        x[k] = v / a[k][k];

        if (!std::isfinite(x[k])) {
            throw std::runtime_error("Поправка метода Ньютона не является конечным числом");
        }
    }

    return x;
}

static std::vector<double> residual(const State& v, const State& old, double dt,
                                    const Settings& settings, const std::vector<double>& f0,
                                    const std::vector<double>& f) {
    std::vector<double> r(old.c.size());

    for (size_t i = 0; i < old.c.size(); ++i) {
        r[i] = v.c[i] - old.c[i]
               - dt * ((1 - settings.theta) * f0[i] + settings.theta * f[i]);
    }

    return r;
}

static double residual_norm(const std::vector<double>& r, const State& v,
                            const State& old, const Settings& settings) {
    double value = 0;

    for (size_t i = 0; i < old.c.size(); ++i) {
        value = std::max(value, std::abs(r[i])
            / (settings.atol + settings.rtol * std::max(old.c[i], v.c[i])));
    }

    return value;
}

State implicit_step(const Mechanism& m, const State& old, double dt, double e0,
                    const Settings& settings, bool reverse) {
    State s = old;
    const auto f0 = rates(m, old, reverse).source;
    const size_t n = old.c.size();

    for (int it = 0; it < settings.max_newton_iterations; ++it) {
        s.temperature = temperature(m, s.c, e0, s.temperature);
        auto f = rates(m, s, reverse, true);
        auto r = residual(s, old, dt, settings, f0, f.source);
        const double before = residual_norm(r, s, old, settings);

        if (before <= settings.nonlinear_tolerance) {
            s.time = old.time + dt;
            return s;
        }

        for (size_t i = 0; i < n; ++i) {
            r[i] = -r[i];

            for (size_t j = 0; j < n; ++j) {
                double identity = 0.0;

                if (i == j) {
                    identity = 1.0;
                }

                f.jacobian[i][j] = identity - dt * settings.theta * f.jacobian[i][j];
            }
        }

        const auto delta = linear_solve(f.jacobian, r);
        bool improved = false;
        double damping = 1;

        for (int line = 0; line < 45; ++line, damping *= 0.5) {
            State trial = s;
            bool nonnegative = true;

            for (size_t i = 0; i < n; ++i) {
                trial.c[i] += damping * delta[i];
                nonnegative = nonnegative && trial.c[i] >= 0 && std::isfinite(trial.c[i]);
            }

            if (!nonnegative) {
                continue;
            }

            try {
                trial.temperature = temperature(m, trial.c, e0, s.temperature);
                const double after = residual_norm(
                    residual(trial, old, dt, settings, f0, rates(m, trial, reverse).source),
                    trial, old, settings);

                if (after < before * (1 - 1e-4 * damping)) {
                    s = trial;
                    improved = true;
                    break;
                }
            } catch (const std::runtime_error&) {
            }
        }

        if (!improved) {
            throw std::runtime_error("Не удалось подобрать шаг метода Ньютона");
        }
    }

    throw std::runtime_error("Метод Ньютона для химической системы не сошёлся");
}

static void csv_header(std::ofstream& out, const Mechanism& m) {
    out << "time_s,temperature_K,pressure_Pa,density_kg_m3"
        << ",specific_enthalpy_J_kg,specific_internal_energy_J_kg"
        << ",specific_entropy_J_kg_K,mean_molar_mass_kg_mol";

    for (const auto& s : m.species) {
        out << ",c_" << s.name << "_mol_m3";
    }

    for (const auto& s : m.species) {
        out << ",mole_fraction_" << s.name;
    }

    for (const auto& s : m.species) {
        out << ",mass_fraction_" << s.name;
    }

    for (size_t j = 0; j < m.reactions.size(); ++j) {
        out << ",reaction_" << j + 1 << "_forward_mol_m3_s"
            << ",reaction_" << j + 1 << "_reverse_mol_m3_s"
            << ",reaction_" << j + 1 << "_net_mol_m3_s";
    }

    out << '\n';
}

struct CsvDiagnostics {
    double specific_enthalpy = 0;
    double specific_internal_energy = 0;
    double specific_entropy = 0;
    double mean_molar_mass = 0;
    std::vector<double> mole_fractions;
    std::vector<double> mass_fractions;
    Rates reaction_rates;
};

static CsvDiagnostics csv_diagnostics(const Mechanism& m, const State& s,
                                      bool explicit_reverse) {
    const double total_concentration = std::accumulate(s.c.begin(), s.c.end(), 0.0);
    const double total_mass = mass(m, s);

    if (!(total_concentration > 0) || !(total_mass > 0)) {
        throw std::runtime_error("Некорректный состав для вычисления свойств смеси");
    }

    CsvDiagnostics result;
    result.mole_fractions.resize(s.c.size());
    result.mass_fractions.resize(s.c.size());
    result.specific_internal_energy = energy(m, s) / total_mass;
    result.mean_molar_mass = total_mass / total_concentration;

    for (size_t i = 0; i < s.c.size(); ++i) {
        const Thermo th = thermo(m.species[i], s.temperature);
        result.specific_enthalpy += s.c[i] * th.h;

        if (s.c[i] > 0) {
            const double partial_pressure = s.c[i] * R * s.temperature;
            result.specific_entropy += s.c[i]
                * (th.s - R * std::log(partial_pressure / standard_pressure));
        }

        result.mole_fractions[i] = s.c[i] / total_concentration;
        result.mass_fractions[i] = s.c[i] * m.species[i].molar_mass / total_mass;
    }

    result.specific_enthalpy /= total_mass;
    result.specific_entropy /= total_mass;
    result.reaction_rates = rates(m, s, explicit_reverse);

    const double scalar_sum = result.specific_enthalpy + result.specific_internal_energy
        + result.specific_entropy + result.mean_molar_mass;

    if (!std::isfinite(scalar_sum)) {
        throw std::runtime_error("Свойства смеси для записи в CSV содержат неконечные числа");
    }

    return result;
}

static void csv_row(std::ofstream& out, const Mechanism& m, const State& s,
                    bool explicit_reverse) {
    const double p = pressure(s);

    if (!std::isfinite(s.time + s.temperature + s.density + p)) {
        throw std::runtime_error("Состояние для записи в CSV содержит неконечные числа");
    }

    for (double c : s.c) {
        if (c < 0 || !std::isfinite(c)) {
            throw std::runtime_error("Некорректная концентрация для записи в CSV");
        }
    }

    const CsvDiagnostics diagnostics = csv_diagnostics(m, s, explicit_reverse);

    out << s.time << ',' << s.temperature << ',' << p << ',' << s.density
        << ',' << diagnostics.specific_enthalpy
        << ',' << diagnostics.specific_internal_energy
        << ',' << diagnostics.specific_entropy
        << ',' << diagnostics.mean_molar_mass;

    for (double c : s.c) {
        out << ',' << c;
    }

    for (double fraction : diagnostics.mole_fractions) {
        out << ',' << fraction;
    }

    for (double fraction : diagnostics.mass_fractions) {
        out << ',' << fraction;
    }

    for (size_t j = 0; j < diagnostics.reaction_rates.forward.size(); ++j) {
        out << ',' << diagnostics.reaction_rates.forward[j]
            << ',' << diagnostics.reaction_rates.reverse[j]
            << ',' << diagnostics.reaction_rates.forward[j]
                          - diagnostics.reaction_rates.reverse[j];
    }

    out << '\n';

    if (!out) {
        throw std::runtime_error("Ошибка записи CSV");
    }
}

static std::map<std::string, double> atom_totals(const Mechanism& m, const State& s) {
    std::map<std::string, double> totals;

    for (size_t i = 0; i < s.c.size(); ++i) {
        for (const auto& [a, n] : m.species[i].atoms) {
            totals[a] += n * s.c[i];
        }
    }

    return totals;
}

Statistics integrate(const Input& in) {
    if (in.csv.has_parent_path()) {
        std::filesystem::create_directories(in.csv.parent_path());
    }

    std::ofstream out(in.csv);

    if (!out) {
        throw std::runtime_error("Не удалось открыть CSV: " + in.csv.string());
    }

    out << std::setprecision(17);
    csv_header(out, in.mechanism);

    State state = in.initial;
    csv_row(out, in.mechanism, state, in.explicit_reverse);

    const double e0 = energy(in.mechanism, state);
    const double mass0 = mass(in.mechanism, state);
    const auto atoms0 = atom_totals(in.mechanism, state);
    const auto& cfg = in.solver;

    Statistics stats;
    double dt = cfg.initial_step;
    size_t attempts = 0;
    std::string reason;

    while (state.time < in.end_time) {
        if (++attempts > cfg.max_steps) {
            throw std::runtime_error("Достигнут предел числа попыток шага");
        }

        dt = std::min(dt, in.end_time - state.time);

        if (state.time + dt == state.time) {
            throw std::runtime_error("Шаг слишком мал для увеличения времени при текущей точности вычислений");
        }

        double error = 0;
        State accepted;

        try {
            const State full = implicit_step(in.mechanism, state, dt, e0,
                                             cfg, in.explicit_reverse);
            const State half = implicit_step(in.mechanism, state, dt / 2, e0,
                                             cfg, in.explicit_reverse);
            accepted = implicit_step(in.mechanism, half, dt / 2, e0,
                                     cfg, in.explicit_reverse);

            for (size_t i = 0; i < state.c.size(); ++i) {
                error = std::max(error, std::abs(accepted.c[i] - full.c[i])
                    / (cfg.atol + cfg.rtol * std::max(state.c[i], accepted.c[i])));
            }

            error = std::max(error, std::abs(accepted.temperature - full.temperature)
                / (cfg.rtol * accepted.temperature));

            if (error > 1) {
                reason = "локальная ошибка превышает допуск";
            }
        } catch (const std::runtime_error& e) {
            error = std::numeric_limits<double>::infinity();
            reason = e.what();
        }

        if (error > 1) {
            ++stats.rejected;
            dt *= 0.5;

            if (dt < cfg.min_step) {
                throw std::runtime_error("Достигнут минимальный шаг при t="
                    + std::to_string(state.time) + ": " + reason);
            }

            continue;
        }
        if (dt == in.end_time - state.time) {
            accepted.time = in.end_time;
        } else {
            accepted.time = state.time + dt;
        }

        state = accepted;
        ++stats.accepted;

        stats.max_mass_error = std::max(stats.max_mass_error,
            std::abs(mass(in.mechanism, state) / mass0 - 1));

        stats.max_energy_error = std::max(stats.max_energy_error,
            std::abs(energy(in.mechanism, state) - e0) / std::max(1.0, std::abs(e0)));

        for (const auto& [a, total] : atom_totals(in.mechanism, state)) {
            stats.max_atom_error = std::max(stats.max_atom_error,
                std::abs(total - atoms0.at(a)) / std::max(1.0, atoms0.at(a)));
        }

        if (stats.max_mass_error > 1e-8
            || stats.max_energy_error > 1e-9
            || stats.max_atom_error > 1e-8) {
            throw std::runtime_error("Проверка законов сохранения не пройдена");
        }

        if (stats.accepted % in.save_every == 0 || state.time == in.end_time) {
            csv_row(out, in.mechanism, state, in.explicit_reverse);
        }

        dt = std::min(cfg.max_step,
            dt * std::clamp(0.85 / std::sqrt(std::max(error, 1e-8)), 0.5, 1.5));
    }

    out.close();

    if (!out) {
        throw std::runtime_error("Ошибка закрытия или сброса буфера CSV");
    }

    stats.final = state;
    return stats;
}
