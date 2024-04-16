#ifndef PARTICLE_H
#define PARTICLE_H
#include <atomic>
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

  const SFloat mass_mev;
  const SFloat mass_kg;
  // mass = energy / c^2
  // https://en.wikipedia.org/wiki/Planck_relation
  // E = h * f
  // f = E / h
  const SFloat frequency = mass_mev / kHEv;
  const SFloat avg_q;         // Average charge. -e for electron. +e for proton.
  // Have the electron charge start random between 0 and +-2e.
  // Since charge is based on a sine wave and will oscillate between 0 and +-2e,
  // need a random starting point between 0 and 2 pi
  const SFloat q_amplitude;   // Amplitude of charge.  Set to e = kQ

  // initial charge ?= (static_cast<SFloat>(rand()) / RAND_MAX) * 2 * M_PI; // NOLINT(*-msc50-cpp)
  // See ChargeSinusoidal for how this is used:
  // return avg_q + (q_amplitude * sin((frequency * time_ * 2 * M_PI) + initial_charge));
  // Want paired particles to be at opposite ends of the sine wave.
  // In other words to be a half cycle away.
  // This causes one electron to be at -2e charge, while the other is at 0 charge.
  const SFloat initial_charge;
  SFloat freq_charge;

  SFloat pos[3] = {0, 0, 0};   // Position
  SFloat pos_magnitude;        // Set for proton.
  SFloat pos_unit_vec[3];      // Set for proton.

  SFloat pos_change[3];        // Change in position.
  SFloat pos_change_magnitude;

  volatile bool dist_calcs_done[sew::kMaxParticles];
  volatile SFloat dist_mag_all [kMaxParticles];  // Cache distance calcs to use by other threads.
  volatile SFloat dist_mag_all2[kMaxParticles];
  // Distance from closest attracted particle.
  // Usually this is the distance from particle exerting largest force.
  // Exception is when closest has zero charge.
  SFloat dist_closest[3];
  SFloat dist_mag_closest;
  SFloat dist_unit_vec[3];   // Unit vector from closest
  SFloat distance_mag_from_origin = 0;
  Particle* par_closest;
  int    par_closest_id;

  // Limit the distance of the particle from 0 to combat energy gain.
  const SFloat max_dist_allow;
  // Used as a hack to limit the problem of energy gain.
  // https://en.wikipedia.org/wiki/Energy_drift
  // Limit the speed of the particle to combat energy gain
  const SFloat max_speed_allow;

  SFloat eforce[3];           // Force caused by electric field
  SFloat forces[3];           // Sum of all forces from all particles.
  SFloat force_magnitude;     // Sum of all forces from all particles.
  SFloat force_mag_closest;   // Force mag from closest attracted particle.  Only used for logging.

  SFloat magnet_fs[3];          // Sum of all magnetic forces from all particles.
  SFloat b_force_by_oth_vel[3]; // Force caused by others velocity.
  SFloat b_f_intrinsic[3];      // Force caused by others intrinsic magnetic field.

  SFloat vel[3] = {0, 0, 0};
  SFloat vel_mag = 0;
  SFloat vel_mag2 = 0;        // Velocity magnitude squared.
  SFloat vel_unit_vec[3];
  SFloat new_dt = 0;          // Global dt = smallest new dt.

  SFloat acceleration[3];
  // Use variable dt to avoid digital simulation errors.
  // Fast fraction informs were we are in the dt range between long and short.
  SFloat fast_fraction = 0;  // Percent we are between long dt and short dt.

  SFloat orig_vel_dot_prod;  // Positive = moving away from origin.  Negative = moving towards origin.
  // Determine if electron is coming or going relative to closest proton.
  // For proton indicates if coming or going relative to center.
  SFloat       dist_vel_dot_prod;
  bool energy_dissipated = false;      // When leaving proton, dissipate energy.
  SFloat percent_energy_dissipated;
  int num_allowed_escapes_for_energy_capping;


  // Infinite force when particles are superimposed on each other.
  // To combat limitations of simulation, teleport the electron to other side of proton.
  SFloat v_when_teleported = 0;
  SFloat prev_fast_fraction = 0;
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
  // std::atomic<SFloat> dist_traveled_since_last_trail_update{0};
  volatile SFloat dist_traveled_since_last_trail_update;

  explicit Particle(int id, bool is_electron, Atom* a,
                    SFloat mass_mev, SFloat mass_kg, SFloat avg_q, SFloat q_amplitude,
                    SFloat max_dist_allowed, SFloat max_speed_allowed,
                    Logger *logger
                    );
  void logToBuffer(const std::string &s);
  void log_prev_log_lines(int max_lines_to_log = kMaxPrevLogLines);
  void ConsiderLoggingToFile(int count);
  void HandleEscape();
  void CheckForEscape();
  // Get the charge that is based on freq and thus time.
  SFloat ChargeSinusoidal() const;

  static SFloat* cross(const SFloat* a, const SFloat* b, SFloat* c) {
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
  MagneticForce(Particle *oth, SFloat e_field_magnitude_oth, SFloat q2,
                         SFloat *dist_vector, SFloat dist_mag_2, const SFloat *dist_unit_vector);

  void CalcForcesFromParticle(Particle* oth /* other particle */);

  // Is there a formula equivalent to below?
  static SFloat dummiesInverseExponential(SFloat fast_fraction) {
    // Slow = how close to short dt.
    // Fast = how close to long dt.
    SFloat slow_fraction = 1 - fast_fraction;

    // from 0 to 90% return 0 to 90%
    // from 90% to 95% return 90% to 99%
    // from 95% to 97% return 99% to 99.9%
    // from 97% to 99% return 99.9% to 99.99%
    // from 99% to 100% return 99.99% to 99.999%

    if (slow_fraction < 0.90f) {
      return slow_fraction;
    }
    if (slow_fraction < 0.95f) {
      return 0.9f + (slow_fraction - 0.90f) * (0.09f / 0.05f);
    }
    if (slow_fraction < 0.97f) {
      return 0.99f + (slow_fraction - 0.95f) * (0.009f / 0.02f);
    }
    if (slow_fraction < 0.99f) {
      return 0.999f + (slow_fraction - 0.97f) * (0.0009f / 0.02f);
    }
    if (slow_fraction < 0.999f) {
      return 0.9999f + (slow_fraction - 0.99f) * (0.00009f / 0.009f);
    }
    return 1;
  }

  // When forces are very high we should slow down delta time (dt) to make calculations more accurate.
  void NewDt(SFloat * fast_fraction_ptr);

  // With little or no space between particles, forces approach infinity.
  // If the simulation won't work because of large errors because of huge forces.
  void TeleportIfTooCloseToProton();
  void ApplyForces();
  void atomic_update(std::atomic<SFloat>* atomic1, SFloat magnitude);
};


class Electron : public Particle {
public:
  explicit Electron(int id, Atom *a, Logger *logger, SFloat max_allowed_dist) :
           Particle(id, true, a, kEMassMEv, kEMassKg, -kQ, -kQ,
                    max_allowed_dist, kMaxSpeedElectron, logger) {}
};

class Proton : public Particle {
public:
  // amplitude = kQ * kBohrMagneton / kProtonMagneticMoment = guess / theory.
  // Could amplitude be kQ?
  explicit Proton(int id, Atom *a, Logger *logger, SFloat max_allowed_dist) :
           Particle(id, false, a, kPMassMEv, kPMassKg, kQ, kQ,
                    max_allowed_dist, kMaxSpeedProton, logger) {}
};

} // namespace

#endif  // PARTICLE_H
