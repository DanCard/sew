#ifndef PARTICLE_H
#define PARTICLE_H
#include <cassert>
#include <cmath>   // For M_PI constant and other match functions such as round.
#include <fstream>
#include <iostream>
#include <string>
// #include <thread>
#include <vector>

#include "atom.h"
#include "constants.h"
#include "logger.h"
#include "tee_stream.cpp"

namespace sew {

// Forward declarations
class Atom;
class Logger;


class Particle {
public:
  const int id;            // Unique id for each particle.
  const bool is_electron;  // Only used for logging.
  Atom* a_;                // The atom this particle belongs to.

  const double mass_mev;
  const double mass_kg;
  // mass = energy / c^2
  // https://en.wikipedia.org/wiki/Planck_relation
  // E = h * f
  // f = E / h
  const double frequency = mass_mev / kHEv;
  const double avg_q;         // Average charge. -e for electron. +e for proton.
  // Have the electron charge start random between 0 and +-2e.
  // Since charge is based on a sine wave and will oscillate between 0 and +-2e,
  // need a random starting point between 0 and 2 pi
  const double q_amplitude;   // Amplitude of charge.  Set to e = kQ

  // initial charge ?= (static_cast<double>(rand()) / RAND_MAX) * 2 * M_PI; // NOLINT(*-msc50-cpp)
  // See ChargeSinusoidal for how this is used:
  // return avg_q + (q_amplitude * sin((frequency * time_ * 2 * M_PI) + initial_charge));
  // Want paired particles to be at opposite ends of the sine wave.
  // In other words to be a half cycle away.
  // This causes one electron to be at -2e charge, while the other is at 0 charge.
  const double initial_charge;
  double freq_charge;

  double pos[3] = {0, 0, 0};   // Position
  double pos_magnitude;        // Set for proton.
  double pos_unit_vec[3];      // Set for proton.

  double pos_change[3];        // Change in position.
  double pos_change_magnitude;

  volatile bool dist_calcs_done[sew::kMaxParticles];
  volatile double dist_mag_all [kMaxParticles];  // Cache distance calcs to use by other threads.
  volatile double dist_mag_all2[kMaxParticles];
  // Distance from closest attracted particle.
  // Usually this is the distance from particle exerting largest force.
  // Exception is when closest has zero charge.
  double dist_closest[3];
  double dist_mag_closest;
  double dist_unit_vec[3];   // Unit vector from closest
  double distance_mag_from_origin = 0;
  Particle* par_closest;
  int    par_closest_id;

  // Limit the distance of the particle from 0 to combat energy gain.
  const double max_dist_allow;
  // Used as a hack to limit the problem of energy gain.
  // https://en.wikipedia.org/wiki/Energy_drift
  // Limit the speed of the particle to combat energy gain
  const double max_speed_allow;

  double eforce[3];           // Force caused by electric field
  double forces[3];           // Sum of all forces from all particles.
  double force_magnitude;     // Sum of all forces from all particles.
  double force_mag_closest;   // Force mag from closest attracted particle.  Only used for logging.

  double magnet_fs[3];          // Sum of all magnetic forces from all particles.
  double b_force_by_oth_vel[3]; // Force caused by others velocity.
  double b_f_intrinsic[3];      // Force caused by others intrinsic magnetic field.

  double vel[3] = {0, 0, 0};
  double vel_mag = 0;
  double vel_mag2 = 0;        // Velocity magnitude squared.
  double vel_unit_vec[3];
  double new_dt = 0;          // Global dt = smallest new dt.

  double acceleration[3];
  // Use variable dt to avoid digital simulation errors.
  // Fast fraction informs were we are in the dt range between long and short.
  double fast_fraction = 0;  // Percent we are between long dt and short dt.

  double orig_vel_dot_prod;  // Positive = moving away from origin.  Negative = moving towards origin.
  // Determine if electron is coming or going relative to closest proton.
  // For proton indicates if coming or going relative to center.
  double       dist_vel_dot_prod;
  double       dis_vel_dot_prod_old;
  bool flipped_dis_vel_dot_prod = false;  // If true then log
  bool energy_dissipated = false;      // When leaving proton, dissipate energy.
  bool energy_dissipated_prev = false;  // When different from energy_dissipated, then log.

  // Infinite force when particles are superimposed on each other.
  // To combat limitations of simulation, teleport the electron to other side of proton.
  double v_when_teleported = 0;
  double prev_fast_fraction = 0;
  int    log_count;      // Force logging around an event.
  static const int kMaxPrevLogLines = 4;
  std::vector<std::string> prev_log_lines = std::vector<std::string>(kMaxPrevLogLines);
  int    prev_log_line_index = 0;
  unsigned char color[3];

  Logger* logger;
  // std::ofstream log_file;
  // TeeLogger tee_logger("log.txt");
  TeeLogger tee_logger;
  std::ostream& tee;
  std::ofstream& log_file = tee_logger.get_file_stream();

  explicit Particle(int id, bool is_electron, Atom* a,
                    double mass_mev, double mass_kg, double avg_q, double q_amplitude,
                    double max_dist_allowed, double max_speed_allowed,
                    Logger *logger
                    );
  void logToBuffer(const std::string &s);
  void log_prev_log_lines(int max_lines_to_log = kMaxPrevLogLines);
  void ConsiderLoggingToFile(int count);
  void HandleEscape();
  void CheckForEscape();
  // Get the charge that is based on freq and thus time.
  double ChargeSinusoidal() const;

  static double* cross(const double* a, const double* b, double* c) {
    c[0] = a[1]*b[2] - b[1]*a[2];
    c[1] = a[2]*b[0] - b[2]*a[0];
    c[2] = a[0]*b[1] - b[0]*a[1];
    return c;
  }

  void InitVarsToCalcForces();

  // Calculate the magnetic force caused by:
  //  1. particle(s) moving
  //  2. Intrinsic magnetic field of the particles.
  // https://phys.libretexts.org/Bookshelves/University_Physics/Calculus-Based_Physics_(Schnick)/Volume_B%3A_Electricity_Magnetism_and_Optics/B17%3A_Magnetic_Field%3A_Causes
  // https://en.wikipedia.org/wiki/Lorentz_force
  // Calculate magnetic field at q1 from q2:
  //  1. Magnetic field of q2 due to it moving
  //  2. Magnetic field of q2 due to intrinsic magnetic field.
  // cross that with
  //  1. Magnetic field of q1 due to it moving
  //  2. Magnetic field of q1 due to intrinsic magnetic field.
  // Does the magnetic field of EM wave 1 interact with the magnetic field of EM wave 2?
  // The electric fields do, so you would think the magnetic fields should also.
  void
  MagneticForce(Particle *oth, double e_field_magnitude_oth, double q2,
                         double *dist_vector, double dist_mag_2, const double *dist_unit_vector);

  void CalcForcesFromParticle(Particle* oth /* other particle */);

  // Is there a formula equivalent to below?
  static double dummiesInverseExponential(double fast_fraction) {
    // Slow = how close to short dt.
    // Fast = how close to long dt.
    double slow_fraction = 1 - fast_fraction;

    // from 0 to 90% return 0 to 90%
    // from 90% to 95% return 90% to 99%
    // from 95% to 97% return 99% to 99.9%
    // from 97% to 99% return 99.9% to 99.99%
    // from 99% to 100% return 99.99% to 99.999%

    if (slow_fraction < 0.90) {
      return slow_fraction;
    }
    if (slow_fraction < 0.95) {
      return 0.9 + (slow_fraction - 0.90) * (0.09 / 0.05);
    }
    if (slow_fraction < 0.97) {
      return 0.99 + (slow_fraction - 0.95) * (0.009 / 0.02);
    }
    if (slow_fraction < 0.99) {
      return 0.999 + (slow_fraction - 0.97) * (0.0009 / 0.02);
    }
    if (slow_fraction < 0.999) {
      return 0.9999 + (slow_fraction - 0.99) * (0.00009 / 0.009);
    }
    return 1;
  }

  // When forces are very high we should slow down delta time (dt) to make calculations more accurate.
  void NewDt(double * fast_fraction_ptr);

  // With little or no space between particles, forces approach infinity.
  // If the simulation won't work because of large errors because of huge forces.
  void TeleportIfTooCloseToProton();

  void ApplyForcesToParticle();
};


class Electron : public Particle {
public:
  explicit Electron(int id, Atom* a, Logger *logger) : Particle(id, true, a,
                    kEMassMEv, kEMassKg, -kQ, -kQ,
                    kBohrRadius * 2, kMaxSpeedElectron, logger) {}
};

class Proton : public Particle {
public:
  // amplitude = kQ * kBohrMagneton / kProtonMagneticMoment = guess / theory.
  // Could amplitude be kQ?
  explicit Proton(int id, Atom* a, Logger *logger) : Particle(id, false, a,
                  kPMassMEv, kPMassKg, kQ, kQ,
                  kBohrRadiusProton, kMaxSpeedProton, logger) {}
};

} // namespace

#endif  // PARTICLE_H
