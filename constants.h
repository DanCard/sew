#ifndef SEW_CONSTANTS_H
#define SEW_CONSTANTS_H
#include <algorithm>
#include <cassert>
#include <cmath>   // For M_PI constant and other match functions such as round.
#include <chrono>  // For logging every 500 ms
#include <cstdlib> // For rand()
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

typedef float SFloat;

namespace sew {

const SFloat kC = 299792458;  // m/s.  Speed of light
const SFloat kCoulomb = 8987551787.3681764;       // N * m^2 / C^2
                                    // https://en.wikipedia.org/wiki/Electron
const SFloat kQ = 1.602176634e-19;  // Charge of a particle in Coulombs
             // Planck's constant eV*s  https://en.wikipedia.org/wiki/Planck_constant
const SFloat kHEv = 4.135667696E-15;
// const SFloat kH = 6.62607015E-34;    // Planck's constant m^2 * kg / s
const SFloat kEMassMEv =    510998.95;  // eV / c^2  https://en.wikipedia.org/wiki/Electron
const SFloat kPMassMEv = 938272088.16;  // eV / c^2  https://en.wikipedia.org/wiki/Electronvolt#Mass
const SFloat kEMassKg = 9.1093837015e-31;   // kg
const SFloat kPMassKg = 1.67262192369e-27;  // kg
const SFloat kBohrRadius = 5.29177210903e-11;  // Meters
const SFloat kBohrRadiusProton = kBohrRadius / 2;  // value = swag / trial and error
const SFloat kLithiumAtomSize = 152e-12;  // 152 picometers
// The nuclear radius (r) can be estimated using the formula r = r0 * A^(1/3),
// where r0 is approximately 1.2 femtometers (fm) and 
// A is the mass number of the nucleus (protons + neutrons).
// https://en.wikipedia.org/wiki/Nuclear_radius

// const SFloat kBohrMagneton = kQ * kH / (4 * M_PI * kEMassKg); // https://en.wikipedia.org/wiki/Bohr_magneton
// const SFloat kProtonMagneticMoment = 1.41060679736e-26;  // J/T . https://en.wikipedia.org/wiki/Proton_magnetic_moment
const SFloat kEFrequency = kEMassMEv / kHEv;
const SFloat kPFrequency = kPMassMEv / kHEv;
const int    kMaxParticles = 16;
// For increased speed avoid lowering.  Instead increase threading.
const int    kPFrequencySubDivisions = 128;
// Ranges for dt = delta time
// Slow the simulation when there are huge forces that create huge errors.
const SFloat kShortDt = 1 / ( kPFrequency * kPFrequencySubDivisions );
// Use a long dt to make the simulation faster.
const SFloat kLongDt  = 1 / ( kEFrequency * kPFrequencySubDivisions );  // Seconds
// We change the simulation style when particle gets near the speed of light.
// Instead of using dt, we just simulate the trajectory of the particle.
// Needed because simulation creates huge errors when there are huge forces.
// https://en.wikipedia.org/wiki/Energy_drift
const SFloat kMaxSpeedElectron = kC / 4;
const SFloat kMaxSpeedProton = kMaxSpeedElectron * (kEMassMEv / kPMassMEv);  // swag / trial and error
// Can't have large dt with large forces, otherwise huge digital error.
// If electron and proton are closer than this then there is trouble due to simulation error.
// Causes short duration dt , logging to increase, and simulation speed to slow down.
const SFloat kCloseToTrouble = 2.3e-13;
const SFloat kForceTooHigh = 0.3;  // 0.3 from trial and error.
const bool   kHoldingProtonSteady = false;  // Don't move the proton(s).

} // namespace

#endif  // SEW_CONSTANTS_H
