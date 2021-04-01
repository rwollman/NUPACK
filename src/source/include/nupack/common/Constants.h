/**
 * @brief Defines physical constants, default random number generator
 *
 * @file Constants.h
 * @author Mark Fornace
 * @date 2018-05-31
 */
#pragma once
#include "Config.h"

namespace nupack {

/// Occasionally useful variable to put in some value without recompiling the whole project
extern std::string HackHelper;

/******************************************************************************************/

constexpr real const ZeroCinK = 273.15;
constexpr real const DefaultTemperature = ZeroCinK + 37.0;
constexpr real const LogOf2 = 0.6931471805599453;
constexpr real const LogOf10 = 2.302585092994046;
constexpr real const Pi = M_PI;

//constexpr real Kb = 0.001987204118; // This is the correct value
constexpr real const Kb = 0.00198717; // This agrees with NUPACK 3
constexpr real const DefaultKT= Kb * DefaultTemperature;

/******************************************************************************************/

/// Boltzmann factor from energy and beta
inline real boltzmann_factor(real beta, real energy) {return std::exp(-beta * energy);}
/// Energy from Boltzmann factor and beta
inline real inverse_boltzmann(real beta, real factor) {return -std::log(factor)/beta;}

/******************************************************************************************/

/// A constant DNA sequence, used for testing when a random one isn't desired
extern std::string ReferenceSequence;
/// Get a constant DNA sequence, used for testing when a random one isn't desired
std::string reference_dna(std::size_t length);

/// moles per liter of water from temperature in Kelvin
real water_molarity(real t);

// Stabilization energy from salt concentration for each loop
real dna_salt_correction(real t, real na, real mg, bool long_helix=false);

/******************************************************************************************/

}
