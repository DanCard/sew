// #include <algorithm>
#include <cassert>
#include <cmath>   // For M_PI constant and other match functions such as round.
#include <cstdlib> // For rand()
//#include <deque>
#include <fstream>
#include <iostream>
#include <iomanip>
// #include <Magnum/Math/Matrix3.h>
#include <sstream>
#include <string>
#include <vector>

#include "tee_stream.cpp"

#pragma clang diagnostic push
#pragma ide diagnostic ignored "UnreachableCallsOfFunction"
namespace {

const double kC = 299792458;  // m/s.  Speed of light
const double kCoulomb = 8987551787.3681764;       // N * m^2 / C^2
const double kQ = 1.602176634e-19;  // Charge of a particle
             // Planck's constant eV*s  https://en.wikipedia.org/wiki/Planck_constant
const double kHEv = 4.135667696E-15;
const double kH = 6.62607015E-34;    // Planck's constant m^2 * kg / s
const double kEMassMEv =    510998.95;  // eV / c^2  https://en.wikipedia.org/wiki/Electron
const double kPMassMEv = 938272088.16;  // eV / c^2  https://en.wikipedia.org/wiki/Electronvolt#Mass
const double kEMassKg = 9.1093837015e-31;   // kg
const double kPMassKg = 1.67262192369e-27;  // kg
const double kBohrRadius = 5.29177210903e-11;  // Meters
                                               // https://en.wikipedia.org/wiki/Bohr_magneton
const double kBohrMagneton = kQ * kH / (4 * M_PI * kEMassKg);
const double kProtonMagneticMoment = 1.41060679736e-26;  // J/T . https://en.wikipedia.org/wiki/Proton_magnetic_moment
const double kBohrRadiusProton = kBohrRadius / 20;  // value = swag / trial and error
const double kEFrequency = kEMassMEv / kHEv;
const double kPFrequency = kPMassMEv / kHEv;
const int    kMaxParticles = 4;
// Ranges for dt = delta time
// Slow the simulation when there are huge forces that create huge errors.
const double kShortDt = .1 / kPFrequency;  // Seconds
// Use a long dt to make the simulation faster.
const double kLongDt  = .1 / kEFrequency;  // Seconds
// We change the simulation style when particle gets near the speed of light.
// Instead of using dt, we just simulate the trajectory of the particle.
// Needed because simulation creates huge errors when there are huge forces.
// https://en.wikipedia.org/wiki/Energy_drift
const double kMaxSpeedElectron = kC / 4;
const double kMaxSpeedProton = kMaxSpeedElectron * (kEMassMEv / kPMassMEv);  // swag / trial and error
// Can't have large dt with large forces, otherwise huge digital error.
// If electron and proton are closer than this then there is trouble due to simulation error.
const double kCloseToTrouble = 5e-14;
const double kForceTooHigh = 0.3;  // From trial and error.


// Global variables
double time_ = 0;
  // Delta time in seconds.
double dt = kLongDt - kShortDt / 2;
int count = 0;         // Invocation count of moving particles.

class Particle {
public:
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
  const double q_amplitude;   // Amplitude of charge.  qe for electron.  e * UB / Up for proton.

  // initial charge ?= (static_cast<double>(rand()) / RAND_MAX) * 2 * M_PI; // NOLINT(*-msc50-cpp)
  // See GetSinusoidalChargeValue for how this is used:
  // return avg_q + (q_amplitude * sin((frequency * time_ * 2 * M_PI) + initial_charge));
  // Want paired particles to be at opposite ends of the sine wave.
  // In other words to be a half cycle away.
  // This causes one electron to be at -2e charge, while the other is at 0 charge.
  const double initial_charge;
  double freq_charge;

  // Used as a hack to limit the problem of energy gain.
  // https://en.wikipedia.org/wiki/Energy_drift
  // Limit the speed of the particle to combat energy gain
  const double max_speed_allow;
  // Limit the distance of the particle from 0 to combat energy gain.
  const double max_dist_allow;
  double pos[3] = {0, 0, 0};  // Position

  double distance_mag_from_origin = 0;
  double distance[kMaxParticles][3];

  double forces[3];           // Sum of all forces from all particles.
  double vel[3] = {0, 0, 0};
  double vel_mag = 0;
  double vel_mag2 = 0;        // Velocity magnitude squared.
  double new_dt = 0;          // Global dt = smallest new dt.
  double pos_change_magnitude;

  // Only used for debug.  Not used in simulation calculations.
  const int id;            // Unique id for each particle.  Just used for logging.
  const bool is_electron;  // Only used for logging.

  // std::ofstream log_file;
  // TeeLogger logger("log.txt");
  TeeLogger logger;
  std::ostream& tee;
  std::ofstream& log_file = logger.get_file_stream();
  // Infinite force when particles are superimposed on each other.
  // To combat limitations of simulation, teleport the electron to other side of proton.
  double v_when_teleported = 0;
  double prev_vel_mag = 0;
  double prev_fast_fraction = 0;
  double prev_energy = 0;
  int log_count = 4;      // Force logging around an event.
  static const int max_prev_log_lines = 4;
  std::vector<std::string> prev_log_lines = std::vector<std::string>(max_prev_log_lines);
  int prev_log_line_index = 0;
  double m_bg_f[3];         // Magnetic force1 from background magnetic field.
  double mf_q2[3];          // Force from q2.
  double acceleration[3];
  // Used for logging.  Not used for simulation calcs.
  double forces_by_each_p[kMaxParticles][3];  // Force from each particle for one iteration.
  double largest_force_mag;   // Largest force mag from a single particle.
  Particle* p_exerting_largest_force;
  int    p_exerting_largest_force_id;
  double dist_mag_from_largest;
  // Use variable dt to avoid digital simulation errors.
  // Fast fraction informs were we are in the dt range between long and short.
  double fast_fraction = 0;  // Percent we are between long dt and short dt.
  bool is_force_too_high = false;



  explicit Particle(double mass_mev, double mass_kg, double avg_q, double q_amplitude,
                    double max_speed_allowed, double max_dist_allowed,
                    int id, bool is_electron) :
          mass_mev(mass_mev), mass_kg(mass_kg), avg_q(avg_q), q_amplitude(q_amplitude),
          // Want paired particles to be at opposite ends of the sine wave.
          // In other words to be a half cycle away.
          // This causes one electron to be at -2e charge, while the other is at 0 charge.
          initial_charge(id%2 == 0 ? 0 : M_PI),
          max_speed_allow(max_speed_allowed),
          max_dist_allow(max_dist_allowed),
          id(id),
          is_electron(is_electron),
          logger(is_electron ? ("electron_" + std::to_string(id) + ".log")
                                : "proton_" + std::to_string(id) + ".log"), tee(logger.get_stream())
          {
    tee << "\t\t frequency " << frequency << "  "
              << (is_electron ? "electron" : "proton") << " mass mev " << mass_mev << std::endl;
    tee << "\t\t initial_charge " << initial_charge << std::endl;
  }

  void logToBuffer(const std::string &s) {
    prev_log_lines[prev_log_line_index] = s;
    prev_log_line_index = (prev_log_line_index + 1) % max_prev_log_lines;
  }

  void log_prev_log_lines(int max_lines_to_log = 8) {
    int index = prev_log_line_index;
    // int limit = std::min(max_lines_to_log, max_prev_log_lines);
    int limit = max_lines_to_log < max_prev_log_lines ? max_lines_to_log : max_prev_log_lines;
    for (int i = 0; i < limit; ++i) {
      if (!prev_log_lines[index].empty()) {
        tee << prev_log_lines[index] << std::endl;
      }
      index = (index + 1) % max_prev_log_lines;
    }
  }


  static bool NotTrivialParameter(const double *d3, int i) {
    if (i == 0) return d3[0] != 0;
    if (d3[0] == 0) return d3[i] != 0;
    // If within 1E4 of zero index.
    return std::abs(d3[0] / d3[i]) < 10;
  }

  static std::string Log3dArray(double *d3, const std::string &name, int width = 8) {
    std::ostringstream log_line;
    // Set precision based on width.  If width is 6, then precision is 0.
    // If width is 8, then precision is 1.  Cause 1 space for the decimal point.
    // If width is 10, then precision is 3.
    int precision = width - 7;
    if (width < 8) precision = 0;
    log_line << "  " << name << std::scientific << std::setprecision(precision);
    if (NotTrivialParameter(d3, 0)) log_line << " x " << std::setw(width) << d3[0];
    if (NotTrivialParameter(d3, 1)) log_line << " y " << std::setw(width) << d3[1];
    if (NotTrivialParameter(d3, 2)) log_line << " z " << std::setw(width) << d3[2];
    if (d3[0] == 0 && d3[1] == 0 && d3[2] == 0) {
      log_line << " 0 " << std::setw(width) << " ";
    }
    return log_line.str();
  }

  void LogStuff() {
    double pos_chng[3];     // Change in position.
    for (int i = 0; i < 3; ++i) {
      pos_chng[i] = this->vel[i] * dt;
      this->pos[i] += pos_chng[i];
    }
    pos_change_magnitude = sqrt(pow(pos_chng[0], 2) + pow(pos_chng[1], 2) + pow(pos_chng[2], 2));

    double danger_speed = max_speed_allow / 2;
    // double potential_energy = kCoulomb * freq_charge * other_charge / dist_mag;
    double kinetic_energy = 0.5 * mass_kg * vel_mag2;
    double total_energy = 0; // potential_energy + kinetic_energy;

    std::ostringstream log_line;
    log_line
      << std::setw( 7) << count << (is_electron ? " e" : " p") << id
      << std::scientific << std::setprecision(1)
      << Log3dArray(pos     , "pos")
      << " largest f from: " << (p_exerting_largest_force->is_electron ? "e" : "p")
                             << p_exerting_largest_force_id
      << " f " << largest_force_mag
      << Log3dArray(forces  , "fs"  ) << std::setprecision(1)
   // << " d mag" << std::setw(9) << dist_mag_from_largest << std::setprecision(1)
   // << Log3dArray(dist    , "dis")
   // << Log3dArray(magnet_f, "B"  )
   // << Log3dArray(m_bg_f, "Bg"   )
   // << Log3dArray(mf_q2, "Bq2"   ) << std::setprecision(1)
      << "  a "    << std::setw( 8) << acceleration[0]
      << " v mag " << vel_mag
      << Log3dArray(vel     , "v", 8)
   // << " chng"  << std::setw(10) << std::setprecision(3) << pos_change
      << Log3dArray(pos_chng, "chng") << std::setprecision(1) << std::fixed
      /*
      << " pe"    << std::setw(10) << potential_energy
      << " ke"    << std::setw(10) << kinetic_energy
      << " te"    << std::setw(10) << total_energy << std::fixed << std::setprecision(1)
      */
   // << " dt"    << std::setw( 9) << /* dt << */ " new " << new_dt
      << " fast"  << std::setw( 5) << round(fast_fraction * 10) * 10 << '%'
      << " chrg"  << std::setw( 4) << int(round((freq_charge/avg_q)*100)) << '%'
   // << " oth"   << std::setw( 4) << int(round((other_charge / oth->q_amplitude) * 100)) << '%'
   // << " inv"   << std::setw(12) << inverse_exponential
    ;
    std::string log_line_str = log_line.str();
    bool particles_are_very_close = dist_mag_from_largest < kCloseToTrouble;
    bool fast_fraction_changed_significantly = (prev_fast_fraction / fast_fraction > 1.1) &&
            (prev_fast_fraction - fast_fraction > 0.1);
    bool energy_changed = false;
    /* && is_electron && (
             std::abs(total_energy / prev_energy)  > 2  // If energy doubles then log.
             ||      (total_energy * prev_energy) <= 0  // If energy flips sign then log.
             );
             */
    const int ll = 2;  // Limit logging to once every x lines.
    if (log_count > 0 // || true
      ||  count%(ll*16) == 0
      || (count%(ll* 8) == 0 && fast_fraction  > 0.01)
      || (count%(ll* 4) == 0 && fast_fraction == 0.0)
      || (count%(ll* 2) == 0 && fast_fraction == 0.0 && particles_are_very_close)
      || (count% ll     == 0 && (is_force_too_high && vel_mag > (danger_speed * 0.95)))
      || fast_fraction_changed_significantly
      || energy_changed
    ) {
      std::string log_source = " ?";
      if (log_count > 0) {
        log_source = " L";
      } else if (fast_fraction_changed_significantly) {
        log_source = " F";
      } else if (energy_changed) {
        log_source = " E";
      } else {
        for (int i = 16; i > 0; i /= 2) {
          if (count % (ll*i) == 0) {
            log_source = " " + std::to_string(i);
            break;
          }
        }
      }
      log_line_str = log_line_str + log_source;
      tee << log_line_str << std::endl;
      prev_energy = total_energy;
      prev_fast_fraction = fast_fraction;
      log_count--;
    }
    logToBuffer(log_line_str);
    prev_vel_mag = vel_mag;
  }

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
    assert(false);
  }

  // When forces are very high we should slow down delta time (dt)
  double CalcNewDt(double * fast_fraction_ptr) {
    // If we are inside the trouble zone that slowdown should be max.
    // In other words dt should be shortest time.
    // fast_fraction < 0 when distance shorter than kCloseToTrouble.
    // fast_fraction > 1 when distance is further than kBohrRadius.
    fast_fraction = (dist_mag_from_largest - kCloseToTrouble) / (kBohrRadius - kCloseToTrouble);
    if (fast_fraction > 1) {
      fast_fraction = 1;
    } else if (fast_fraction < 0) {
      fast_fraction = 0;
    }
    *fast_fraction_ptr = fast_fraction;
    double inverse_exponential = dummiesInverseExponential(fast_fraction);
    if (inverse_exponential > 1) {
      inverse_exponential = 1;
    } else if (inverse_exponential < 0) {
      inverse_exponential = 0;
    }
    double new_dt = kLongDt + ((kShortDt - kLongDt) * inverse_exponential);
    return new_dt;
  }

  void HandleEscape() {
    log_prev_log_lines();
    tee << "\tEscaped past "
      << (is_electron ? " bohr radius of " : " proton allowed radius of ") << max_dist_allow
      << "  Distance magnitude: " << distance_mag_from_origin
      << "  V when previously teleported: " << v_when_teleported
      << "  Current velocity " << vel[0] << " " << vel[1] << " " << vel[2]
      ;
    // tee << "  pos " << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
    assert(pos[0] != 0 || pos[1] != 0 || pos[2] != 0);
    // We want the magnitude to be under the bohr radius.
    // Distance to shave off.  Subtract 1% to keep things safe.
    double dist_to_shave = ( max_dist_allow / distance_mag_from_origin ) - 0.01;
    // tee << "  dist_to_shave " << dist_to_shave;
    for (double & po : pos) {
      po = po * dist_to_shave;
    }
    // tee << "  pos " << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
    // Zero out the biggest velocity or all?
    for (double & v : vel) v = 0;
    log_count = 2;  // Force logging around this event.
    tee << std::endl;
  }

  // Get the charge that is based on freq and thus time.
  double GetSinusoidalChargeValue() const {
    // if dt > 0.1 / frequency, then charge is constant.
    // Can't simulate sinusoidal charge when dt is large.
    /*
    if (dt > (0.1/frequency)) {
      return avg_q;
    }
    */
    /*
    if (!is_electron) {
      return avg_q;
    }
    if (is_electron) {
      return avg_q;
    }
    */
    // return avg_q;
    // Need to calculate where the particle is in its frequency.
    // Theory is that charge varies between 0 and -2e each time period = 1 / frequency.
    // Frequency is kEFrequency or kPFrequency.
    // 
    // Charge = e + e * sin(2 * pi * frequency * time + initial_charge)
    return avg_q + (q_amplitude * sin((frequency * time_ * 2 * M_PI) + initial_charge));
  }

  static double* cross(const double* a, const double* b, double* c) {
    c[0] = a[1]*b[2] - b[1]*a[2];
    c[1] = a[2]*b[0] - b[2]*a[0];
    c[2] = a[0]*b[1] - b[0]*a[1];
    return c;
  }

  // Calculate the magnetic force1 caused by:
  // 1. particle(s) moving
  // 2. Intrinsic magnetic field of the particles.
  //    Not implemented yet.
  // 3. background magnetic field.  Negligible.  Could delete it.
  // https://phys.libretexts.org/Bookshelves/University_Physics/Calculus-Based_Physics_(Schnick)/Volume_B%3A_Electricity_Magnetism_and_Optics/B17%3A_Magnetic_Field%3A_Causes
  // https://en.wikipedia.org/wiki/Lorentz_force
  // Calculate magnetic field at point of q1:
  // 1. Magnetic field of q2 due to it moving
  // 2. Magnetic field of q2 due to intrinsic magnetic field.
  // cross that with
  // 1. Magnetic field of q1 due to it moving
  // 2. Magnetic field of q1 due to intrinsic magnetic field.
  void
  CalculateMagneticForce(double *force, Particle *oth, double q1, double q2,
                         double *dist_vector, double dist_mag_2) {
    // Constants for calculating magnetic force1.
    // https://academic.mu.edu/phys/matthysd/web004/l0220.htm
    // Permeability of free space.  https://en.wikipedia.org/wiki/Vacuum_permeability
    const double kPermeability  = 4 * M_PI * 1e-7;  // T * m / A
    const double kPermeabilityDividedBy4Pi = 1e-7;  // T * m / A
    const double kBackgroundMagneticField = 1e-5;   // T = Tesla
    const double kBackgroundMagneticFieldVector[3] = {0, 0,
                                                         kBackgroundMagneticField};
    // tee << "\t vel " << std::setprecision(1) << vel[0] << " " << vel[1] << " " << vel[2] << '\t';

    double v1_cross_background[3];
    double v2_cross_r[3];
    cross(     vel, kBackgroundMagneticFieldVector, v1_cross_background);
    cross(oth->vel,                    dist_vector, v2_cross_r);
    double b_field_by_q2[3];  // Magnetic field created by q2.
    for (int i = 0; i < 3; ++i) {
      m_bg_f[i] = q1 * v1_cross_background[i];
      // https://en.wikipedia.org/wiki/Biot%E2%80%93Savart_law
      // https://docs.google.com/document/d/1Tc1AZltK30YA-u-corns3hVO7Zcx22uYE1TyT4i_NNU
      b_field_by_q2[i] = kPermeabilityDividedBy4Pi * q2 * v2_cross_r[i] / dist_mag_2;
      assert(!std::isnan(b_field_by_q2[i]));
    }
    double v1_cross_field_by_q2[3];
    cross(vel, b_field_by_q2, v1_cross_field_by_q2);
    for (int i = 0; i < 3; ++i) {
      // Lorentz force1 from q2 magnetic field = q * v x B
      mf_q2[i] = q1 * v1_cross_field_by_q2[i];
      force[i] = m_bg_f[i] + mf_q2[i];
    }
    /*
    tee << " magnetic field q2 "
              << b_field_by_q2[0] << ' ' << b_field_by_q2[1] << ' ' << b_field_by_q2[2];
    tee << " magnetic force1 from q2 " << mf_q2[0] << ' ' << mf_q2[1] << ' ' << mf_q2[2];
    tee << "\t background f " << m_bg_f[0] << ' ' << m_bg_f[1] << ' ' << m_bg_f[2] << std::endl;
    */
  }


  // Calculate the force1 between two particles and update position and velocity.
  void SetPosition(Particle* oth /* other particle */) {
    distance_mag_from_origin = sqrt(pow(pos[0], 2) + pow(pos[1], 2) + pow(pos[2], 2));
    // If particle escapes, then zero out the velocity.
    // This is a hack to limit the problem of energy gain.
    // https://en.wikipedia.org/wiki/Energy_drift
    if (distance_mag_from_origin > max_dist_allow) {
      HandleEscape();
    }
    double dist[3];
    for (int i = 0; i < 3; ++i) {
      dist[i] = pos[i] - oth->pos[i];
    }
    double dist_mag2 = pow(dist[0], 2) + pow(dist[1], 2) + pow(dist[2], 2);
    double dist_mag = sqrt(dist_mag2);
    double dist_unit_vector[3];
    for (int i = 0; i < 3; ++i) {
      dist_unit_vector[i] = dist[i] / dist_mag;
    }
    // Current charge varies based on frequency and time.
    double freq_q = GetSinusoidalChargeValue();
    double oth_charge = oth->GetSinusoidalChargeValue();
    double force[3];
    // When opposite charges, force1 is negative.  When same charges, force1 is positive.
    double eforce_magnitude = kCoulomb * freq_q * oth_charge / dist_mag2;
    double magnet_f[3];     // Magnetic force1.
    CalculateMagneticForce(magnet_f, oth, freq_q, oth_charge, dist, dist_mag2);
    for (int i = 0; i < 3; ++i) {
      force[i] = (eforce_magnitude * dist_unit_vector[i]) + magnet_f[i];
      acceleration[i] = force[i] / this->mass_kg;
      this->vel[i] += acceleration[i] * dt;
    }
    vel_mag2 = pow(vel[0], 2) + pow(vel[1], 2) + pow(vel[2], 2);
    vel_mag = sqrt(vel_mag2);
    // Determine xyz component with highest velocity.
    int max_v_component = 0;
    {
      double abs_max_v = 0;
      for (int i = 0; i < 3; ++i) {
        if (std::abs(vel[i]) > abs_max_v) {
          abs_max_v = std::abs(vel[i]);
          max_v_component = i;
        }
      }
    }
    // With little or no space between particles, forces approach infinity.
    // If the simulation won't work because of large errors because of huge forces.
    // 0.2 = number pulled out of thin air.
    // Higher number = faster particle gets ejected.
    // Lower number = more likely to get teleported, loop around proton, get closer.
    bool is_force_too_high = std::abs(eforce_magnitude) > 0.3;
    bool is_speed_too_high = vel_mag > max_speed_allow;
    bool force_same_dir_as_v = (vel[max_v_component] > 0 && force[max_v_component] > 0) ||
                               (vel[max_v_component] < 0 && force[max_v_component] < 0);
    if (is_force_too_high &&
        is_speed_too_high &&
        // If the force1 is opposite the direction of the velocity, then don't worry about it.
        force_same_dir_as_v) {
      oth->log_prev_log_lines(1);
      log_prev_log_lines();
      tee << "\t\t " << (is_electron ? "electron" : "proton") << id
        << " force1 " << force[0] << " " << force[1] << " " << force[2] << std::endl;
      tee << "\t Force too high.  Dist " << dist_mag << "  Teleporting to other side of "
        << (oth->is_electron ? "electron" : "proton")
        << ".  pos " << pos[0] << " " << pos[1] << " " << pos[2]
        << "  oth pos " << oth->pos[0] << " " << oth->pos[1] << " " << oth->pos[2]
        << std::endl;
      // Teleport the particle to other side of the other particle.
      for (int i = 0; i < 3; ++i) {
        pos[i] += -2 * dist[i];
      }
      tee << "\t Teleported to " << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
      v_when_teleported = vel_mag;
      log_count = 4;               // Force extra logging after this event.
      oth->log_file << std::endl;  // Drop a hint that something important happened.
      // oth->log_count = 2;          // Force extra logging for proton.
    }

    double fast_fraction;
    new_dt = CalcNewDt(&fast_fraction);

    LogStuff();
  }

  void CalcForcesFromParticle(Particle* oth /* other particle */) {
    int oth_id = oth->id;
    distance_mag_from_origin = sqrt(pow(pos[0], 2) + pow(pos[1], 2) + pow(pos[2], 2));
    // If particle escapes, then zero out the velocity.
    // This is a hack to limit the problem of energy gain.
    // https://en.wikipedia.org/wiki/Energy_drift
    if (distance_mag_from_origin > max_dist_allow) {
      HandleEscape();
    }
    double dist[3];
    for (int i = 0; i < 3; ++i) {
      dist[i] = pos[i] - oth->pos[i];
      distance[oth_id][i] = dist[i];
    }
    double dist_mag2 = pow(dist[0], 2) + pow(dist[1], 2) + pow(dist[2], 2);
    double dist_mag  = sqrt(dist_mag2);
    double dist_unit_vector[3];
    for (int i = 0; i < 3; ++i) {
      dist_unit_vector[i] = dist[i] / dist_mag;
    }
    // Current charge varies based on frequency and time.
    freq_charge = GetSinusoidalChargeValue();
    double other_charge = oth->GetSinusoidalChargeValue();
    double force[3];
    // When opposite charges, force1 is negative.  When same charges, force1 is positive.
    double eforce_magnitude = kCoulomb * freq_charge * other_charge / dist_mag2;
    double magnet_f[3];     // Magnetic force1.
    CalculateMagneticForce(magnet_f, oth, freq_charge, other_charge, dist, dist_mag2);
    double force_magnitude;
    for (int i = 0; i < 3; ++i) {
      force[i] = (eforce_magnitude * dist_unit_vector[i]) + magnet_f[i];
      // Below are the main return values used.
      forces[i] += force[i];  // Used for logging.  Not used for simulation calcs.
      forces_by_each_p[oth_id][i] += force[i];  // Used for logging.  Not used for simulation calcs.
    }
    force_magnitude = sqrt(pow(force[0], 2) + pow(force[1], 2) + pow(force[2], 2));
    if (force_magnitude > largest_force_mag) {
      p_exerting_largest_force = oth;
      p_exerting_largest_force_id = oth_id;
      largest_force_mag = force_magnitude;
      is_force_too_high = std::abs(force_magnitude) > kForceTooHigh;
    }
  }

  void InitVarsToCalcForces() {
    largest_force_mag = 0;
    is_force_too_high = false;
    p_exerting_largest_force = nullptr;
    for (int i = 0; i < kMaxParticles; ++i) {
      for (int j = 0; j < 3; ++j) {
        forces_by_each_p[i][j] = 0;
      }
    }
    for (int i=0; i<3; ++i) {
      forces[i] = 0;
    }
  }

  void ClearLoggingForMultipleIterations() {
    return;
    for (int i = 0; i < kMaxParticles; ++i) {
      for (int j = 0; j < 3; ++j) {
      }
    }
  }

  void ApplyForcesToParticle(double * forces_ptr) {
    for (int i = 0; i < 3; ++i) {
      acceleration[i] = forces_ptr[i] / this->mass_kg;
      this->vel[i] += acceleration[i] * dt;
    }
    vel_mag2 = pow(vel[0], 2) + pow(vel[1], 2) + pow(vel[2], 2);
    vel_mag  = sqrt(vel_mag2);
    double fast_fraction;
    new_dt = CalcNewDt(&fast_fraction);
  }
};


class Electron : public Particle {
public:
  explicit Electron(int id) : Particle(kEMassMEv, kEMassKg, -kQ, -kQ,
   kMaxSpeedElectron, kBohrRadius, id, true) {}
};

class Proton : public Particle {
public:
  // amplitude = kQ * kBohrMagneton / kProtonMagneticMoment = guess / theory.
  // Could amplitude be kQ?
  explicit Proton(int id) : Particle(kPMassMEv, kPMassKg, kQ, kQ,
   kMaxSpeedProton, kBohrRadiusProton, id, false) {}
};

class __attribute__((unused)) Particles {
public:
  Particle * pars[kMaxParticles];  // par[ticle]s
  int num_particles = 0;

  explicit Particles(int numParticles) : num_particles(numParticles) {
    if (num_particles == 0) return;
    std::cout << "\t kEFrequency " << kEFrequency << "  kPFrequency " << kPFrequency << std::endl;
    std::cout << "\t\t kBohrMagneton " << kBohrMagneton
              << "  kProtonMagneticMoment " << kProtonMagneticMoment << std::endl;

    std::cout << "\t\t kShortDt " << kShortDt << "  kLongDt " << kLongDt << std::endl;
    assert(numParticles <= kMaxParticles);
    for (int i = 0; i < numParticles; ++i) {
      if (i < numParticles/2) {
        pars[i] = new Electron(i);
      } else {
        pars[i] = new Proton(i);
      }
      /* Set up:
          Want to have one electron with charge 0 all the way to the left.
          Two protons one high and one low.
          Electron in the middle with charge -2e.
      */
      if (i == 0) {
        pars[0]->pos[0] = -kBohrRadius * 0.75;
      } else if (i > 1) {
        pars[i]->pos[1] = (kBohrRadiusProton / 2) * (i%2 == 0 ? 1 : -1);
        /*
        for (double & po : pars[i]->pos) {
          po = (kBohrRadiusProton / 2) * (static_cast<double>((float) std::rand() / RAND_MAX));
        }
        */
      }
      std::cout << "\t" << (pars[i]->is_electron ? "electron" : "proton") << " " << i
                << "  pos " << pars[i]->pos[0] << " " << pars[i]->pos[1] << " " << pars[i]->pos[2]
                << std::endl;
    }
    assert(pars[1]->frequency <= kEFrequency);
  }

  void LogParticles() {
    for (int i=0; i<num_particles; ++i) {
      pars[i]->LogStuff();
      /*
      Particle * part_ptr = pars[i];
      part_ptr->log_file << "Forces from each particle for one iteration." << std::endl;
      for (int j=0; j<num_particles; ++j) {
        part_ptr->log_file << "  " << j << "  " << force_mag[i][j] << std::endl;
      }
      part_ptr->log_file << "Forces from each particle over multiple iterations." << std::endl;
      part_ptr->log_file << "Largest force from a single particle: " << largest_force_from[i] << std::endl;
       */
    }
    count++;
  }

  void CalcForcesOnParticle(Particle * part_ptr, double * forces_ptr) {
    int part_num = part_ptr->id;
    part_ptr->InitVarsToCalcForces();
    for (int i = 0; i < num_particles; ++i) {
      if (i == part_num) continue;
      part_ptr->CalcForcesFromParticle(pars[i]);
    }
  }

  void moveParticles() {
    const double min_pos_change_desired = 1e-14;
    const int kMaxTimesToGetSignificantMovement = 1024;
    double pos_change[kMaxParticles];
    std::cout << "\t\t\t\t time: " << time_ << std::endl;
    for (int j = 0; j < num_particles; ++j) {
      pos_change[j] = 0;
    }
    double forces[num_particles][3];
    int part_with_most_movement = 0;
    for (int i = 0; i < num_particles; ++i) {
      pars[i]->ClearLoggingForMultipleIterations();
    }
    for (int i = 0; i < kMaxTimesToGetSignificantMovement; ++i) {
      for (int j = 0; j < num_particles; ++j) {
        Particle * part_ptr = pars[j];
        CalcForcesOnParticle(part_ptr, forces[j]);
        part_ptr->ApplyForcesToParticle(forces[j]);
        pos_change[j] += std::abs(pars[j]->pos_change_magnitude);
      }
      LogParticles();
      // tee << " change " << pos_change << std::endl;
      time_ += dt;
      // Find the shortest dt and set the new dt to that.
      dt = pars[0]->new_dt;
      int shortest_dt = 0;
      for (int j = 1; j < num_particles; ++j) {
        if (pars[j]->new_dt < dt) {
          dt = pars[j]->new_dt;
          shortest_dt = j;
        }
      }
      std::cout << "\t shortest_dt " << shortest_dt << std::endl;
      for (int j = 0; j < num_particles; ++j) {
        if (pos_change[j] > min_pos_change_desired) {
          part_with_most_movement = j;
          goto exit_loop;
        }
      }
    }
    exit_loop:
    // if (pos_change < min_pos_change_desired) tee << std::endl;
    // tee << "  pos: " << pos << std::endl;
  }

};

} // namespace

#pragma clang diagnostic pop