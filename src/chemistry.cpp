#include "chemistry.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

Thermo thermo(const Species& sp, double T) {
    if (!std::isfinite(T) || T < sp.thermo.front().Tmin || T > sp.thermo.back().Tmax) {
        throw std::runtime_error("Температура вне термодинамического диапазона для " + sp.name);
    }

    const ThermoInterval* interval = &sp.thermo.back();

    for (const auto& p : sp.thermo) {
        if (T <= p.Tmax) {
            interval = &p;
            break;
        }
    }

    const auto& b = interval->b;
    const double x = T * 1e-4;

    const double h = b[2] - b[3] - 1e4 * b[7] - 2e4 * b[6] / x
                     + x * (1e4 * b[5]
                            + x * (1e4 * b[8]
                                   + x * (2e4 * b[9] + x * 3e4 * b[10])));

    const double cp = b[5] + 2 * b[6] / (x * x)
                      + x * (2 * b[8] + x * (6 * b[9] + x * 12 * b[10]));

    const double s = b[4] + b[5] + b[5] * std::log(x) - b[6] / (x * x)
                     + x * (2 * b[8] + x * (3 * b[9] + x * 4 * b[10]));

    if (!(cp > R) || !std::isfinite(h + s + cp)) {
        throw std::runtime_error("Некорректные термодинамические данные для " + sp.name);
    }

    return {cp, cp - R, h, s, h - R * T};
}

void prepare(Mechanism& m) {
    if (m.species.empty() || m.reactions.empty()) {
        throw std::runtime_error("Пустой механизм реакций");
    }

    m.Tmin = 0;
    m.Tmax = std::numeric_limits<double>::max();

    for (const auto& s : m.species) {
        if (!(s.molar_mass > 0) || s.atoms.empty() || s.thermo.empty()) {
            throw std::runtime_error("Некорректные данные вещества " + s.name);
        }

        for (size_t k = 0; k < s.thermo.size(); ++k) {
            const auto& p = s.thermo[k];

            if (!(p.Tmin > 0 && p.Tmax > p.Tmin)
                || (k && p.Tmin != s.thermo[k - 1].Tmax)) {
                throw std::runtime_error("Некорректные термодинамические интервалы для " + s.name);
            }

            for (double v : p.b) {
                if (!std::isfinite(v)) {
                    throw std::runtime_error("Коэффициент не является конечным числом");
                }
            }
        }

        m.Tmin = std::max(m.Tmin, s.thermo.front().Tmin);
        m.Tmax = std::min(m.Tmax, s.thermo.back().Tmax);
    }

    if (m.Tmin >= m.Tmax) {
        throw std::runtime_error("Нет общего температурного интервала");
    }

    m.stoichiometry.clear();

    for (const auto& r : m.reactions) {
        std::map<int, int> net;

        for (const auto& t : r.reactants) {
            net[t.species] -= t.coefficient;
        }

        for (const auto& t : r.products) {
            net[t.species] += t.coefficient;
        }

        std::map<std::string, int> atoms;
        std::vector<Term> terms;
        double dm = 0;
        double scale = 0;

        for (auto [i, nu] : net) {
            if (i < 0 || size_t(i) >= m.species.size()) {
                throw std::runtime_error("Индекс вещества вне допустимого диапазона");
            }

            if (nu) {
                terms.push_back({i, nu});
            }

            for (const auto& [a, n] : m.species[i].atoms) {
                atoms[a] += nu * n;
            }

            dm += nu * m.species[i].molar_mass;
            scale += std::abs(nu) * m.species[i].molar_mass;
        }

        for (const auto& [a, n] : atoms) {
            if (n) {
                throw std::runtime_error("Реакция не сохраняет число атомов элемента " + a);
            }
        }

        if (std::abs(dm) > 1e-12 * scale) {
            throw std::runtime_error("Реакция не сохраняет массу");
        }

        m.stoichiometry.push_back(terms);
    }
}

double arrhenius(const Arrhenius& a, double T) {
    if (a.A == 0) {
        return 0;
    }

    const double t = std::max(T, 300.0);
    const double k = std::exp(std::log(a.A) + a.n * std::log(t)
                              - a.activation_temperature / t) * (T / t);

    if (!std::isfinite(k)) {
        throw std::runtime_error("Переполнение при вычислении константы Аррениуса");
    }

    return k;
}

RateConstants rate_constants(const Mechanism& m, size_t j, double T, bool explicit_reverse) {
    const auto& r = m.reactions[j];
    double dg = 0;
    double dn = 0;

    for (auto t : m.stoichiometry[j]) {
        const auto th = thermo(m.species[t.species], T);
        dg += t.coefficient * (th.h - T * th.s);
        dn += t.coefficient;
    }

    const double logK = -dg / (R * T) + dn * std::log(standard_pressure / (R * T));
    const double kf = arrhenius(r.forward, T);
    double kr = 0;

    if (r.reversible) {
        if (explicit_reverse) {
            if (!r.reverse) {
                throw std::runtime_error("Не заданы явные параметры обратной реакции");
            }

            kr = arrhenius(*r.reverse, T);
        } else if (kf > 0) {
            kr = std::exp(std::log(kf) - logK);
        }
    }

    if (!std::isfinite(kr)) {
        throw std::runtime_error("Переполнение константы скорости обратной реакции");
    }

    return {kf, kr, logK};
}

static double product(const std::vector<Term>& terms,
                      const std::vector<double>& c,
                      int derivative = -1) {
    double value = 1;

    for (auto t : terms) {
        int power = t.coefficient;

        if (t.species == derivative) {
            power -= 1;
            value *= t.coefficient;
        }

        for (int k = 0; k < power; ++k) {
            value *= c[t.species];
        }
    }

    return value;
}

static double dlogk(const Arrhenius& a, double T) {
    if (T < 300) {
        return 1 / T;
    }

    return a.n / T + a.activation_temperature / (T * T);
}

Rates rates(const Mechanism& m, const State& s, bool explicit_reverse, bool jacobian) {
    const size_t n = m.species.size();
    Rates result;
    result.source.assign(n, 0);

    if (jacobian) {
        result.jacobian.assign(n, std::vector<double>(n, 0));
    }

    std::vector<double> u(n);
    double cv = 0;

    for (size_t i = 0; i < n; ++i) {
        if (s.c[i] < 0 || !std::isfinite(s.c[i])) {
            throw std::runtime_error("Некорректная концентрация");
        }

        const auto th = thermo(m.species[i], s.temperature);
        u[i] = th.u;
        cv += s.c[i] * th.cv;
    }

    if (!(cv > 0)) {
        throw std::runtime_error("Неположительная объёмная теплоёмкость");
    }

    for (size_t j = 0; j < m.reactions.size(); ++j) {
        const auto& r = m.reactions[j];
        const auto k = rate_constants(m, j, s.temperature, explicit_reverse);
        const double pf = product(r.reactants, s.c);
        const double pr = product(r.products, s.c);
        double collider = 1;

        if (!r.efficiencies.empty()) {
            collider = std::inner_product(s.c.begin(), s.c.end(),
                                          r.efficiencies.begin(), 0.0);
        }

        const double wf = k.forward * pf * collider;
        const double wr = k.reverse * pr * collider;

        if (!std::isfinite(wf + wr)) {
            throw std::runtime_error("Переполнение скорости реакции");
        }

        result.forward.push_back(wf);
        result.reverse.push_back(wr);

        for (auto t : m.stoichiometry[j]) {
            result.source[t.species] += t.coefficient * (wf - wr);
        }

        if (!jacobian) {
            continue;
        }

        double du = 0;

        for (auto t : m.stoichiometry[j]) {
            du += t.coefficient * u[t.species];
        }
        double dr;

        if (explicit_reverse && r.reversible) {
            dr = dlogk(*r.reverse, s.temperature);
        } else {
            dr = dlogk(r.forward, s.temperature) - du / (R * s.temperature * s.temperature);
        }

        const double dwdT = wf * dlogk(r.forward, s.temperature) - wr * dr;

        for (size_t i = 0; i < n; ++i) {
            double dw = 0;

            for (auto t : r.reactants) {
                if (size_t(t.species) == i) {
                    dw += k.forward * product(r.reactants, s.c, int(i)) * collider;
                }
            }

            for (auto t : r.products) {
                if (size_t(t.species) == i) {
                    dw -= k.reverse * product(r.products, s.c, int(i)) * collider;
                }
            }

            if (!r.efficiencies.empty()) {
                dw += r.efficiencies[i] * (k.forward * pf - k.reverse * pr);
            }

            dw -= dwdT * u[i] / cv;

            for (auto t : m.stoichiometry[j]) {
                result.jacobian[t.species][i] += t.coefficient * dw;
            }
        }
    }

    return result;
}

double energy(const Mechanism& m, const State& s) {
    double value = 0;

    for (size_t i = 0; i < s.c.size(); ++i) {
        value += s.c[i] * thermo(m.species[i], s.temperature).u;
    }

    return value;
}

double mass(const Mechanism& m, const State& s) {
    double value = 0;

    for (size_t i = 0; i < s.c.size(); ++i) {
        value += s.c[i] * m.species[i].molar_mass;
    }

    return value;
}

double pressure(const State& s) {
    return R * s.temperature * std::accumulate(s.c.begin(), s.c.end(), 0.0);
}

static std::array<double, 2> energy_residual(const Mechanism& m,
                                              const std::vector<double>& c,
                                              double target, double T) {
    double e = 0;
    double cv = 0;

    for (size_t i = 0; i < c.size(); ++i) {
        const auto th = thermo(m.species[i], T);
        e += c[i] * th.u;
        cv += c[i] * th.cv;
    }

    return std::array<double, 2>{e - target, cv};
}

double temperature(const Mechanism& m, const std::vector<double>& c,
                   double target, double guess) {
    double lo = m.Tmin;
    double hi = m.Tmax;

    if (energy_residual(m, c, target, lo)[0] > 0 || energy_residual(m, c, target, hi)[0] < 0) {
        throw std::runtime_error("Решение уравнения энергии вне термодинамического диапазона");
    }

    double T = std::clamp(guess, lo, hi);

    for (int it = 0; it < 150; ++it) {
        const auto v = energy_residual(m, c, target, T);

        if (!(v[1] > 0)) {
            throw std::runtime_error("Неположительная теплоёмкость при решении уравнения энергии");
        }
        if (std::abs(v[0]) <= 1e-14 * std::max(1.0, std::abs(target))) {
            return T;
        }

        if (v[0] > 0) {
            hi = T;
        } else {
            lo = T;
        }

        const double next = T - std::clamp(v[0] / v[1], -100.0, 100.0);

        if (next > lo && next < hi) {
            T = next;
        } else {
            T = (lo + hi) / 2;
        }
    }

    throw std::runtime_error("Метод Ньютона для уравнения энергии не сошёлся");
}
