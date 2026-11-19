#pragma once

#include <array>
#include <map>
#include <optional>
#include <string>
#include <vector>

inline constexpr double R = 8.3144;
inline constexpr double standard_pressure = 101325.0;

struct Term {
    int species;
    int coefficient;
};

using Stoichiometry = std::vector<std::vector<Term>>;

struct ThermoInterval {
    double Tmin;
    double Tmax;
    std::array<double, 12> b;
};

struct Species {
    std::string name;
    double molar_mass;
    std::map<std::string, int> atoms;
    std::vector<ThermoInterval> thermo;
};
struct Arrhenius {
    double A;
    double n;
    double activation_temperature;
};

struct Reaction {
    std::vector<Term> reactants;
    std::vector<Term> products;
    Arrhenius forward;
    std::optional<Arrhenius> reverse;
    bool reversible = true;
    std::vector<double> efficiencies;
};

struct Mechanism {
    std::vector<Species> species;
    std::vector<Reaction> reactions;
    Stoichiometry stoichiometry;
    double Tmin = 0;
    double Tmax = 0;
};

struct Thermo {
    double cp;
    double cv;
    double h;
    double s;
    double u;
};

struct State {
    double time = 0;
    double temperature = 0;
    double density = 0;
    std::vector<double> c;
};

struct RateConstants {
    double forward;
    double reverse;
    double log_Kc;
};

struct Rates {
    std::vector<double> forward;
    std::vector<double> reverse;
    std::vector<double> source;
    std::vector<std::vector<double>> jacobian;
};

Thermo thermo(const Species& species, double T);
void prepare(Mechanism& mechanism);
double arrhenius(const Arrhenius& parameters, double T);

RateConstants rate_constants(const Mechanism& mechanism, size_t reaction,
                             double T, bool explicit_reverse);

Rates rates(const Mechanism& mechanism, const State& state,
            bool explicit_reverse, bool jacobian = false);

double energy(const Mechanism& mechanism, const State& state);
double mass(const Mechanism& mechanism, const State& state);
double pressure(const State& state);

double temperature(const Mechanism& mechanism, const std::vector<double>& c,
                   double target_energy, double guess);
