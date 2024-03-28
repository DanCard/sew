// #include <algorithm>
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

#include "tee_stream.cpp"

#pragma clang diagnostic push
#pragma ide diagnostic ignored "UnreachableCallsOfFunction"
namespace {

const double kC = 299792458;  // m/s.  Speed of light
const double kCoulomb = 8987551787.3681764;       // N * m^2 / C^2
                                    // https://en.wikipedia.org/wiki/Electron
const double kQ = 1.602176634e-19;  // Charge of a particle in Coulombs
             // Planck's constant eV*s  https://en.wikipedia.org/wiki/Planck_constant
const double kHEv = 4.135667696E-15;
const double kH = 6.62607015E-34;    // Planck's constant m^2 * kg / s
const double kEMassMEv =    510998.95;  // eV / c^2  https://en.wikipedia.org/wiki/Electron
const double kPMassMEv = 938272088.16;  // eV / c^2  https://en.wikipedia.org/wiki/Electronvolt#Mass
const double kEMassKg = 9.1093837015e-31;   // kg
const double kPMassKg = 1.67262192369e-27;  // kg
const double kBohrRadius = 5.29177210903e-11;  // Meters
const double kBohrRadiusProton = kBohrRadius / 2;  // value = swag / trial and error
                                               
const double kBohrMagneton = kQ * kH / (4 * M_PI * kEMassKg); // https://en.wikipedia.org/wiki/Bohr_magneton
const double kProtonMagneticMoment = 1.41060679736e-26;  // J/T . https://en.wikipedia.org/wiki/Proton_magnetic_moment
const double kEFrequency = kEMassMEv / kHEv;
const double kPFrequency = kPMassMEv / kHEv;
const int    kMaxParticles = 8;
const int    kPFrequencySubDivisions = 16;
// Ranges for dt = delta time
// Slow the simulation when there are huge forces that create huge errors.
const double kShortDt = 1 / ( kPFrequency * kPFrequencySubDivisions );
// Use a long dt to make the simulation faster.
const double kLongDt  = 1 / ( kEFrequency * kPFrequencySubDivisions );  // Seconds
// We change the simulation style when particle gets near the speed of light.
// Instead of using dt, we just simulate the trajectory of the particle.
// Needed because simulation creates huge errors when there are huge forces.
// https://en.wikipedia.org/wiki/Energy_drift
const double kMaxSpeedElectron = kC / 4;
const double kMaxSpeedProton = kMaxSpeedElectron * (kEMassMEv / kPMassMEv);  // swag / trial and error
// Can't have large dt with large forces, otherwise huge digital error.
// If electron and proton are closer than this then there is trouble due to simulation error.
// Causes short duration dt , logging to increase, and simulation speed to slow down.
const double kCloseToTrouble = 2.3e-13;
const double kForceTooHigh = 0.3;  // 0.3 from trial and error.
const bool   kHoldingProtonSteady = false;  // Don't move the proton(s).

// Global variables
// Delta time in seconds.
double dt = kLongDt;   // Hopefully particles are far enough apart to use long dt.
double time_ = 0;
int count = 0;         // Invocation count of moving particles.

int num_particles_ = 0;
std::chrono::_V2::system_clock::time_point last_log_time;

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
  const double q_amplitude;   // Amplitude of charge.  Set to e = kQ

  // initial charge ?= (static_cast<double>(rand()) / RAND_MAX) * 2 * M_PI; // NOLINT(*-msc50-cpp)
  // See GetSinusoidalChargeValue for how this is used:
  // return avg_q + (q_amplitude * sin((frequency * time_ * 2 * M_PI) + initial_charge));
  // Want paired particles to be at opposite ends of the sine wave.
  // In other words to be a half cycle away.
  // This causes one electron to be at -2e charge, while the other is at 0 charge.
  const double initial_charge;
  double freq_charge;

  double pos[3] = {0, 0, 0};  // Position
  double pos_change_3d[3];     // Change in position.
  double pos_change_magnitude;

  volatile double dist_mag [kMaxParticles];  // Cache distance calcs to use by other threads.
  volatile double dist_mag2[kMaxParticles];
  // Distance from closest attracted particle.
  // Usually this is the distance from particle exerting largest force.
  // Exception is when closest has zero charge.
  double dist_closest[3];
  double dist_mag_closest;
  double dist_unit_vec[3];   // Unit vector from closest
  double distance_mag_from_origin = 0;
  Particle* p_closest_attracted;
  int    p_closest_attracted_id;

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

  double magnet_fs[3];  // Sum of all magnetic forces from all particles.
  double b_force_by_oth_vel[3];        // Force caused by others velocity.
  double b_force_by_oth_intrinsic[3];  // Force caused by others intrinsic magnetic field.

  double vel[3] = {0, 0, 0};
  double vel_mag = 0;
  double vel_mag2 = 0;        // Velocity magnitude squared.
  double vel_unit_vec[3];
  double new_dt = 0;          // Global dt = smallest new dt.

  // Only used for debug.  Not used in simulation calculations.
  const int id;            // Unique id for each particle.  Just used for logging.
  const bool is_electron;  // Only used for logging.

  double acceleration[3];
  // Use variable dt to avoid digital simulation errors.
  // Fast fraction informs were we are in the dt range between long and short.
  double fast_fraction = 0;  // Percent we are between long dt and short dt.

  // Used for logging to determine of electron is coming or going relative to closest proton.
  double       dis_vel_dot_prod;
  double       dis_vel_dot_prod_old;
  bool flipped_dis_vel_dot_prod = false;

  // std::ofstream log_file;
  // TeeLogger logger("log.txt");
  TeeLogger logger;
  std::ostream& tee;
  std::ofstream& log_file = logger.get_file_stream();
  // Infinite force when particles are superimposed on each other.
  // To combat limitations of simulation, teleport the electron to other side of proton.
  double v_when_teleported = 0;
  double prev_fast_fraction = 0;
  int    log_count;      // Force logging around an event.
  static const int kMaxPrevLogLines = 4;
  std::vector<std::string> prev_log_lines = std::vector<std::string>(kMaxPrevLogLines);
  int    prev_log_line_index = 0;
  unsigned char color[3];
  double total_kinetic_e_all_w = 0;

  volatile bool dist_calcs_done[kMaxParticles];

  explicit Particle(double mass_mev, double mass_kg, double avg_q, double q_amplitude,
                    double max_speed_allowed, double max_dist_allowed,
                    int id, bool is_electron) :
          mass_mev(mass_mev), mass_kg(mass_kg), avg_q(avg_q), q_amplitude(q_amplitude),
          // Want paired particles to be at opposite ends of the sine wave.
          // In other words to be a half cycle away.
          // This causes one electron to be at -2e charge, while the other is at 0 charge.
          initial_charge(id%2 == 0 ? 0 : M_PI),
          max_dist_allow(max_dist_allowed),
          max_speed_allow(max_speed_allowed),
          id(id),
          is_electron(is_electron),
          logger(is_electron ? ("e" + std::to_string(id) + ".log")
                              : "p" + std::to_string(id) + ".log"), tee(logger.get_stream())
  {
    log_count = is_electron ? 2 : 0;
    std::cout << "\t" << (is_electron ? "electron" : "proton") << " " << id
      << "\t frequency " << frequency << "  "
              << (is_electron ? "electron" : "proton") << " mass mev " << mass_mev << std::endl;
    tee << "\t\t initial_charge " << initial_charge << std::endl;
    for (double & p_energy_cycle_ : pot_energy_cycle) {
      p_energy_cycle_ = 0;
    }
  }

  void logToBuffer(const std::string &s) {
    prev_log_lines[prev_log_line_index] = s;
    prev_log_line_index = (prev_log_line_index + 1) % kMaxPrevLogLines;
  }

  void SetColorForConsole() {
    printf("\x1b[38;2;%d;%d;%dm", color[0], color[1], color[2]);
  }

  void log_prev_log_lines(int max_lines_to_log = 8) {
    SetColorForConsole();
    int index = prev_log_line_index;
    // int limit = std::min(max_lines_to_log, kMaxPrevLogLines);
    int limit = max_lines_to_log < kMaxPrevLogLines ? max_lines_to_log : kMaxPrevLogLines;
    for (int i = 0; i < limit; ++i) {
      if (!prev_log_lines[index].empty()) {
        tee << prev_log_lines[index] << " P" << std::endl;
      }
      index = (index + 1) % kMaxPrevLogLines;
    }
  }


  static bool IsSignificantParameter(const double *d3, int i) {
    if (std::abs(d3[i]) < 1e-20) return false;
    if (i == 0) return d3[0] != 0;
    if (d3[0] == 0) return d3[i] != 0;
    // If within 1E4 of zero index.
    return std::abs(d3[0] / d3[i]) < 10;
  }

  static std::string Log3dArray(double *d3, const std::string &name, int width = 6) {
    std::ostringstream log_line;
    // Set precision based on width.  If width is 6, then precision is 0.
    // If width is 8, then precision is 1.  Cause 1 space for the decimal point.
    // If width is 10, then precision is 3.
    int precision = width - 7;
    if (width < 8) precision = 0;
    log_line << "  " << name << std::scientific << std::setprecision(precision);
    if (IsSignificantParameter(d3, 0)) log_line << " x " << std::setw(width) << d3[0];
    else                               log_line << "   " << std::setw(width) << " ";
    if (IsSignificantParameter(d3, 1)) log_line << " y " << std::setw(width) << d3[1];
    else                               log_line << "   " << std::setw(width) << " ";
    if (IsSignificantParameter(d3, 2)) log_line << " z " << std::setw(width) << d3[2];
    else                               log_line << "   " << std::setw(width) << " ";
    return log_line.str();
  }


  double potential_energy_average;
  // Below block of variables only used in below method.
  static const int kManyBuffersMultiplier = 32;
  double pot_energy_cycle[kPFrequencySubDivisions];
  int    pot_energy_cycle_index = 0;
  double sum_p_energy = 0;        // Only used in below method.

  void CalcAveragePotentialEnergy() {
    sum_p_energy += pot_energy_cycle[pot_energy_cycle_index++];
    pot_energy_cycle_index = pot_energy_cycle_index % kPFrequencySubDivisions;
    // Subtract the potential energy that is going to roll off our average.
    sum_p_energy -= pot_energy_cycle[pot_energy_cycle_index];
    // -1 because of the buffer subtracted above.
    potential_energy_average = sum_p_energy / (kPFrequencySubDivisions - 1);
  }
  
  void ConsiderLoggingToFile(int num_drawing_event_already, int num_wait_for_drawing_event) {
    double danger_speed = max_speed_allow / 2;

    bool particles_are_close      = dist_mag_closest < (kCloseToTrouble * 2);
    bool particles_are_close_very = dist_mag_closest < kCloseToTrouble;
    bool fast_fraction_changed_significantly = (prev_fast_fraction / fast_fraction > 1.1) &&
            (prev_fast_fraction - fast_fraction > 0.1);
    const int ll = 1;  // Limit logging to once every x lines.
    if (log_count > 0 // || true
      ||  count%(ll*32) == 0
      || (count%(ll*16) == 0 && fast_fraction < 0.1)
      || (count%(ll* 8) == 0 && fast_fraction < 0.02)
      || (count%(ll* 4) == 0 && fast_fraction < 0.002  && particles_are_close)
      || (count%(ll* 2) == 0 && fast_fraction < 0.0002 && particles_are_close_very)
      || (count% ll     == 0 && vel_mag > danger_speed && particles_are_close_very)
      || fast_fraction_changed_significantly
    ) {
      std::string log_source = " ?";
      if (log_count > 0) {
        log_source = " L";
      } else if (fast_fraction_changed_significantly) {
        log_source = " F";
      } else {
        for (int i = 16; i > 0; i /= 2) {
          if (count % (ll*i) == 0) {
            log_source = " " + std::to_string(i);
            break;
          }
        }
      }
      std::string log_line_str = FormatLogLine(num_drawing_event_already, num_wait_for_drawing_event);
      log_file << log_line_str << log_source << std::endl;
      prev_fast_fraction = fast_fraction;
      log_count--;
      logToBuffer(log_line_str);
    }
  }

  std::string FormatLogLine(int num_drawing_event_already, int num_wait_for_drawing_event) {
    double charge_of_closest = p_closest_attracted->freq_charge;
    double total_energy = potential_energy_average + total_kinetic_e_all_w;

    if (dt < kShortDt*0.9) {
      if (dt < kShortDt / 3 ) {
        std::cout << "dt " << dt << " is too short.  " << std::endl;
      }
    }
    std::ostringstream log_line;
    log_line
      << std::setw( 7) << count <<   (is_electron ? " e" : " p") << id
      << "⋅" << (p_closest_attracted->is_electron ? "e" : "p") << p_closest_attracted_id
      << std::scientific << std::setprecision(3)
      << "  dis"  << std::setw(10) << dist_mag_closest
      << "  vel " << std::setw(10) << vel_mag
   // << Log3dArray(vel, "v")
      << std::fixed << std::setprecision(1)
      << "  d⋅v " << std::setw( 4) << dis_vel_dot_prod    // -1 = approaching, 1 = leaving
      << std::scientific << std::setprecision(2)
      << (dt < (kShortDt*0.9) ? " *" : "  ")
      << "dt"     << std::setw( 9) << dt // << " new " << new_dt
      << " fast"  << std::setw(10) << std::setprecision(3) << fast_fraction
      << " f "    << std::setw( 9) << force_mag_closest
   // << Log3dArray(forces  , " fs")
      << " Bv " << sqrt(pow(b_force_by_oth_vel[0], 2) + pow(b_force_by_oth_vel[1], 2) +
                        pow(b_force_by_oth_vel[2], 2))
      << " Bi " << std::setw( 9) << sqrt(pow(b_force_by_oth_intrinsic[0], 2) +
                                         pow(b_force_by_oth_intrinsic[1], 2) +
                                         pow(b_force_by_oth_intrinsic[2], 2))
      << Log3dArray(b_force_by_oth_intrinsic, "Bi")
      << Log3dArray(magnet_fs, "tB"  )
   // << "  B "   << sqrt(magnet_fs[0]*magnet_fs[0] + magnet_fs[1]*magnet_fs[1] + magnet_fs[2]*magnet_fs[2])
   // << Log3dArray(acceleration, "a")
   // << " chng"  << std::setw(10) << std::setprecision(3) << pos_change_magnitude
   // << Log3dArray(pos_change_3d, "chng") << std::setprecision(1)
   // << Log3dArray(pos     , "pos")
   // << " min pos change " << min_pos_change_desired
   // << round(fast_fraction * 10) * 10 << '%'
                  << std::setw( 6) << std::setprecision(1) << std::fixed 
      << " chrg"  << std::setw( 4) << int(round((freq_charge/avg_q)*100)) << '%'
      << " oth"   << std::setw( 4)
      << int(round((charge_of_closest / p_closest_attracted->q_amplitude) * 100)) << '%'
   // << " inv"   << std::setw(12) << inverse_exponential
      << std::scientific << std::setprecision(2)
      // P energy goes negative, that is why width is larger.
      << "  pe"   << std::setw(10) << potential_energy_average
   // << " "      << std::setw( 8) << p_energy_many_average
      << " ke"    << std::setw( 9) << total_kinetic_e_all_w
      << " te"    << std::setw( 9) << total_energy
      << " L "    << std::setw( 2) << num_drawing_event_already   // Late.  Drawing event already occurred.
      << " E "    << std::setw( 2) << num_wait_for_drawing_event  // Calcs were early.  Waited on drawing event.
    ;
    /*
    for (int i=0; i<num_particles_; ++i) {
      if (i == id || i == p_closest_attracted_id) continue;
      log_line << (i >= (num_particles_/2) ? " p" : "  e") << i << " " << dist_mag[i];
    }
    */
    return log_line.str();
  }


  void HandleEscape() {
    log_prev_log_lines();
    tee << "\t" << (is_electron ? "electron" : "proton")
        << " escaped past radius of " << max_dist_allow
        << "  Dist mag from origin : " << distance_mag_from_origin
        << "  Current velocity " << vel[0] << " " << vel[1] << " " << vel[2];
    if (v_when_teleported > 0) {
        tee << "  V when previously teleported: " << v_when_teleported;
    }
    // tee << "  pos " << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
    assert(pos[0] != 0 || pos[1] != 0 || pos[2] != 0);
    // We want the magnitude to be under the bohr radius.
    // Distance to shave off.  Subtract 1% to keep things safe.
    double dist_to_shave = ( max_dist_allow / distance_mag_from_origin ) - 0.05;
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

  void CheckForEscape() {
    distance_mag_from_origin = sqrt(pow(pos[0], 2) + pow(pos[1], 2) + pow(pos[2], 2));
    // If particle escapes, then zero out the velocity.
    // This is a hack to limit the problem of energy gain.
    // https://en.wikipedia.org/wiki/Energy_drift
    if (distance_mag_from_origin > max_dist_allow) {
      HandleEscape();
    }
  }


  // Get the charge that is based on freq and thus time.
  double GetSinusoidalChargeValue() const {
    // if dt > 0.25 / frequency, then charge is constant.
    // Can't simulate sinusoidal charge when dt is large.
    /*
    if (dt > (0.25/frequency)) {
      return avg_q;
    }
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

  void InitVarsToCalcForces() {
    dist_mag_closest  = 1;      // 1 represents a large number that will be overwritten
                                // quickly when checking for closest.
    p_closest_attracted = nullptr;
    for (double & force : forces) {
      force = 0;
    }
    for (double & magnet_force : magnet_fs) {
      magnet_force = 0;
    }
    for (int i=0; i<num_particles_; ++i) {
      dist_mag [i] = 0;       // Will assert if we test and find these values.
      dist_mag2[i] = 0;
    }
  }

  // Calculate the magnetic force caused by:
  // 1. particle(s) moving
  // 2. Intrinsic magnetic field of the particles.
  //    Not implemented yet.
  // https://phys.libretexts.org/Bookshelves/University_Physics/Calculus-Based_Physics_(Schnick)/Volume_B%3A_Electricity_Magnetism_and_Optics/B17%3A_Magnetic_Field%3A_Causes
  // https://en.wikipedia.org/wiki/Lorentz_force
  // Calculate magnetic field at point of q1:
  // 1. Magnetic field of q2 due to it moving
  // 2. Magnetic field of q2 due to intrinsic magnetic field.
  // cross that with
  // 1. Magnetic field of q1 due to it moving
  // 2. Magnetic field of q1 due to intrinsic magnetic field.
  void
  MagneticForce(Particle *oth, double e_field_magnitude_oth, double q2,
                         double *dist_vector, double dist_mag_2) {
    // Constants for calculating magnetic force.
    // https://academic.mu.edu/phys/matthysd/web004/l0220.htm
    // Permeability of free space.  https://en.wikipedia.org/wiki/Vacuum_permeability
    const double kPermeability  = 4 * M_PI * 1e-7;  // T * m / A
    const double kPermeabilityDividedBy4Pi = 1e-7;  // T * m / A

    double vel_oth_cross_dis[3];
    cross(oth->vel, dist_vector, vel_oth_cross_dis);
    double b_field_by_oth_vel[3];    // Magnetic field created by other particles velocity
    for (int i = 0; i < 3; ++i) {
      // https://en.wikipedia.org/wiki/Biot%E2%80%93Savart_law
      // https://docs.google.com/document/d/1Tc1AZltK30YA-u-corns3hVO7Zcx22uYE1TyT4i_NNU
      b_field_by_oth_vel[i] = kPermeabilityDividedBy4Pi * q2 * vel_oth_cross_dis[i] / dist_mag_2;
      assert(!std::isnan(b_field_by_oth_vel[i]));
    }
    double vel_cross_b_field[3];
    cross(vel, b_field_by_oth_vel, vel_cross_b_field);
    for (int i = 0; i < 3; ++i) {
      // Lorentz force from oth particle magnetic field = q * v x B
      b_force_by_oth_vel[i] = freq_charge * vel_cross_b_field[i];
      magnet_fs[i] += b_force_by_oth_vel[i];   // Only used for logging.
    }
    // b = e / c
    // http://hyperphysics.phy-astr.gsu.edu/hbase/Waves/emwv.html
    // https://www.se.edu/kfrinkle/wp-content/uploads/sites/89/2013/10/main23v130303.pdf
    // Deriving E=cB and B=µ0·ε0·c·E using Maxwell's Equations: https://www.youtube.com/watch?v=dW1q8XQcEdI
    // Magnetic field created by other.
    double b_field_caused_by_intrinsic_oth_magnitude = e_field_magnitude_oth / kC;
    for (int i = 0; i < 3; ++i) {
      // (i+1)%3 = rotate 90 degrees.
      b_force_by_oth_intrinsic[(i+1)%3] = dist_vector[i] * b_field_caused_by_intrinsic_oth_magnitude;
    }
    for (int i = 0; i < 3; ++i) {
      magnet_fs[i] += b_force_by_oth_intrinsic[i];
    }
    /*
    tee << " magnetic field q2 "
              << b_field_by_oth_vel[0] << ' ' << b_field_by_oth_vel[1] << ' ' << b_field_by_oth_vel[2];
    tee << " magnetic force from q2 " << magnet_f_oth[0] << ' ' << magnet_f_oth[1] << ' ' << magnet_f_oth[2];
    tee << "\t background f " << m_bg_f[0] << ' ' << m_bg_f[1] << ' ' << m_bg_f[2] << std::endl;
    */
  }

  void CalcForcesFromParticle(Particle* oth /* other particle */) {
    int oth_id = oth->id;
    if (oth_id == id) return;
    double dist[3];
    for (int i = 0; i < 3; ++i) {
      dist[i] = pos[i] - oth->pos[i];
    }
    double dist_magn;
    double dist_magn2;
    // If we have already done the calcs, then use the cached values in the other standing wave.
    if (oth->dist_calcs_done[oth_id]) {
      dist_magn  = oth->dist_mag [id];
      dist_magn2 = oth->dist_mag2[id];
      dist_mag [oth_id] = dist_magn;
      dist_mag2[oth_id] = dist_magn2;
      dist_calcs_done[id] = true;
    } else {
      dist_magn2 = pow(dist[0], 2) + pow(dist[1], 2) + pow(dist[2], 2);
      dist_magn  = sqrt(dist_magn2);
      dist_mag [oth_id] = dist_magn;  // Alternatively could just set dist_mag for other rather than mess with arrays,
      dist_mag2[oth_id] = dist_magn2; // but that may have concurrency issue(s).
      dist_calcs_done[oth_id] = true;
    }
    assert(dist_magn  > 0);  // Since we can't handle infinite forces, lets assume this is always true.
    assert(dist_magn2 > 0);  // Since we can't handle infinite forces, lets assume this is always true.

    double dist_unit_vector[3];
    for (int i = 0; i < 3; ++i) {
      dist_unit_vector[i] = dist[i] / dist_magn;
    }
    // Current charge varies based on frequency and time.
    double other_charge = oth->freq_charge;
    double e_field_magnitude_oth = kCoulomb * other_charge / dist_magn;
    // When opposite charges, force is negative.  When same charges, force is positive.
    double eforce_magnitude = freq_charge * e_field_magnitude_oth / dist_magn;
    for (int i = 0; i < 3; ++i) {
      eforce[i] = eforce_magnitude * dist_unit_vector[i];
    }
    // Possible optimization is to calculate magnetic force in separate thread.
    // Another possibility is to skip it or use cached values,
    // it doesn't change significantly, since it is negligible / insignificant.
    MagneticForce(oth, e_field_magnitude_oth, other_charge, dist, dist_magn2);

    for (int i = 0; i < 3; ++i) {
      forces[i] += eforce[i] + magnet_fs[i];
    }
    bool forces_attract = is_electron != oth->is_electron;
    // Determine which attractive particle we are closest to.
    // This is often the particle exerting the largest force.
    // Used for logging.
    if (dist_magn < dist_mag_closest && forces_attract) {
      p_closest_attracted    = oth;
      p_closest_attracted_id = oth_id;
      dist_mag_closest = dist_magn;
      for (int i = 0; i < 3; ++i) {
        dist_closest[i] = dist[i];
      }
      force_mag_closest = force_magnitude;
    }
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
  }

  // When forces are very high we should slow down delta time (dt) to make calculations more accurate.
  void NewDt(double * fast_fraction_ptr) {
    // Simulation will slow down alot below this distance.  Will make simulation more accurate.
    const double kSmallDtDistance = kCloseToTrouble * 2;
    // Simulation will slow down when particles are closer than this.  Will make simulation more accurate.
    const double safe_distance = kBohrRadius;
    // If we are inside the trouble zone that slowdown should be max.
    // In other words dt should be shortest time.
    // fast_fraction < 0 when distance shorter than kCloseToTrouble.
    // fast_fraction > 1 when distance is further than (kBohrRadius / 8).
    fast_fraction = (dist_mag_closest - kSmallDtDistance) / (safe_distance - kSmallDtDistance);
    double inverse_exponential;
    if (fast_fraction > 1) {
      fast_fraction = 1;
      inverse_exponential = 0;
    } else if (fast_fraction < 0) {
      fast_fraction = 0;
      inverse_exponential = 1;
    } else {
      inverse_exponential = dummiesInverseExponential(fast_fraction);
      if (inverse_exponential > 1) {
        inverse_exponential = 1;
      } else if (inverse_exponential < 0) {
        inverse_exponential = 0;
      }
    }
    *fast_fraction_ptr = fast_fraction;
    new_dt = kLongDt + ((kShortDt - kLongDt) * inverse_exponential);
    // bool heading_away = is_electron ? dis_vel_dot_prod >= 0 : dis_vel_dot_prod <= 0;
    // Trouble with energy gain.  To compensate lower simulation speed, increase accuracy
    // when leaving close proton.  Decrease dt so that we can lower speed as much as was gained coming
    // in.  Hack to limit problem of energy gain.
    if (dist_mag_closest < (kSmallDtDistance*32) && dis_vel_dot_prod >= 0 && is_electron) {
      // new_dt = new_dt / 2;
      if (dis_vel_dot_prod_old < 0) {
        flipped_dis_vel_dot_prod = true;
      }
    }
  }

  // With little or no space between particles, forces approach infinity.
  // If the simulation won't work because of large errors because of huge forces.
  void TeleportIfTooCloseToProton() {
    if (!is_electron) return;
    // Higher number = faster particle gets ejected.
    // Lower number = more likely to get teleported, loop around proton, get closer.
    bool is_force_too_high = force_magnitude > kForceTooHigh;
    if (!is_force_too_high) return;
    if (vel_mag < max_speed_allow) return;

    // Force needs to be in the same direction as velocity.
    // If not then we don't to worry about particles heading in opposite directions.
    // Need unit vector for velocity and force, then do dot product.
    Particle* close = p_closest_attracted;
    // If the vectors are pointing in the same direction (aligned), their dot product is 1.
    // If they are perpendicular, it's 0. If they are pointing in opposite directions, it's -1.
    // Since we are using distance as a proxy for force, sign is reversed.
    bool force_same_dir_as_v = dis_vel_dot_prod < -0.75;    // -0.75 = guess
    if (!force_same_dir_as_v) return;

    if (num_particles_ > 2)
    close->log_prev_log_lines(2);
           log_prev_log_lines();
    tee << "\t Force too high. " << (is_electron ? "e" : "p") << id << "  dot product " << dis_vel_dot_prod
        << "  f mag " << force_magnitude
        << " x " << forces[0] << " y " << forces[1] << " z " << forces[2] << std::endl;
    tee << "\t\t\t velocity unit vector: " << vel_unit_vec[0]
                                  << " y " << vel_unit_vec[1]
                                 << " "    << vel_unit_vec[2]
        << "\t distance unit vector: " << dist_unit_vec[0]
                              << " y " << dist_unit_vec[1]
                              << " "   << dist_unit_vec[2] << std::endl;
    tee << "  Dist " << dist_mag_closest << "  Teleporting to other side of "
        << (close->is_electron ? "electron" : "proton")
        << ".  pos " << pos[0] << " " << pos[1] << " " << pos[2]
        << "  close pos " << close->pos[0] << " " << close->pos[1] << " " << close->pos[2]
        << std::endl;
    // Teleport the particle to other side of the close by particle.
    for (int i = 0; i < 3; ++i) {
      pos[i] += -2 * dist_closest[i];
    }
    tee << "\t Teleported to " << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
    v_when_teleported = vel_mag;
    log_count = 8;               // Force extra logging after this event.
    close->log_file << "\t\t\t" << (is_electron ? "e" : "p") << id << " teleported "  << std::endl;
    if (num_particles_ > 2)
      close->log_count = 2;          // Force extra logging for other particle.
  }

  void ApplyForcesToParticle() {
    force_magnitude = sqrt(pow(forces[0], 2) + pow(forces[1], 2) + pow(forces[2], 2));
    for (int i = 0; i < 3; ++i) {
      acceleration[i] = forces[i] / this->mass_kg;
      vel[i] += acceleration[i] * dt;
    }
    vel_mag2 = pow(vel[0], 2) + pow(vel[1], 2) + pow(vel[2], 2);
    vel_mag  = sqrt(vel_mag2);
    for (int i = 0; i < 3; ++i) {
      pos_change_3d[i]  = vel[i] * dt;
      if (!kHoldingProtonSteady || is_electron) {
        pos        [i] += pos_change_3d[i];
      }
      vel_unit_vec [i]  = vel[i] / vel_mag;
      dist_unit_vec[i]  = dist_closest[i] / dist_mag_closest;
    }
    dis_vel_dot_prod_old = dis_vel_dot_prod;
    flipped_dis_vel_dot_prod = false;

    // dot product between dist of closest and velocity.
    dis_vel_dot_prod = dist_unit_vec[0] * vel_unit_vec[0] +
                       dist_unit_vec[1] * vel_unit_vec[1] +
                       dist_unit_vec[2] * vel_unit_vec[2];
    // Use dot product between velocity and distance of closest to determine
    // if coming or going relative to closest attracting particle.
    // If the vectors are pointing in the same direction (aligned), their dot product is 1.
    // If they are perpendicular, it's 0. If they are pointing in opposite directions, it's -1.
    NewDt(&fast_fraction);

    // pos_change_magnitude used to decide if we have moved enough for current frame.
    pos_change_magnitude = sqrt(pow(pos_change_3d[0], 2) + pow(pos_change_3d[1], 2) + pow(pos_change_3d[2], 2));
    // When particles are on top of each other forces approach infinity.
    // To get around that problem we teleport to other side of particle.
    // Alternative approach is to have a preprogrammed path and not bother with force calcs.
    if (is_electron) {
      TeleportIfTooCloseToProton();
    }
  }
};


class Electron : public Particle {
public:
  explicit Electron(int id) : Particle(kEMassMEv, kEMassKg, -kQ, -kQ,
   kBohrRadius * 2, kMaxSpeedElectron, id, true) {}
};

class Proton : public Particle {
public:
  // amplitude = kQ * kBohrMagneton / kProtonMagneticMoment = guess / theory.
  // Could amplitude be kQ?
  explicit Proton(int id) : Particle(kPMassMEv, kPMassKg, kQ, kQ,
   kBohrRadiusProton, kMaxSpeedProton, id, false) {}
};


class Particles {
public:
  Particle * pars[kMaxParticles];  // par[ticle]s
  int num_particles = 0;

  volatile bool screen_draw_event_occurred = true;
  // Only used for logging.
  // If we are doing lots of calcs because electron close to proton, then we will wait for
  // draw event before finishing calcs.  That way we get maximum movement per frame.
  // In other words: simulation is too damn slow, and we want speed up so use all of the
  // compute budget per frame.
  volatile int num_drawing_event_already  = 0; // Drawing event already occurred, no need to wait.
  volatile int num_wait_for_drawing_event = 0;

  explicit Particles(int numParticles) : num_particles(numParticles) {
    if (num_particles == 0) return;
    const int num_threads = std::thread::hardware_concurrency();
    std::cout << "\t\t\t num threads " << num_threads << std::endl;
    num_particles_ = num_particles;
    std::cout << "\t\t\t max speed electron " << kMaxSpeedElectron << "  kPFrequencySubDivisions " << kPFrequencySubDivisions
              << "\t kShortDt " << kShortDt << "  kLongDt " << kLongDt << std::endl;
    std::cout << "\t\t\t\t kEFrequency " << kEFrequency << "  kPFrequency " << kPFrequency << std::endl;
    // std::cout << "\t\t kBohrMagneton " << kBohrMagneton << "  kProtonMagneticMoment " << kProtonMagneticMoment << std::endl;

    assert(numParticles <= kMaxParticles);
    int divider;  // Prefer bright colors, but with many particles becomes indistinguishable.
         if (num_particles <= 2)  divider = 2;
    else if (num_particles <= 4)  divider = 3;
    else if (num_particles <= 6)  divider = 4;
    else                          divider = 5;  // With more particles don't brighten as much.
    for (int i = 0; i < numParticles; ++i) {
      Particle* p;
      if (i < numParticles / 2) {
        pars[i] = new Electron(i);
        p = pars[i];
        // Set pseudo random colors.  Electrons tend to be more red.  Protons tend to be more blue.
        // Prefer bright colors over dark colors.
        p->color[0] = 151 + (std::rand() % 105);
        p->color[1] = 100 + (std::rand() % 155);
        p->color[2] =   0 + (std::rand() % 240);
        if (i == 0) {
          p->pos[0] = - kBohrRadius;
        }
      } else {
        pars[i] = new Proton(i);
        p = pars[i];
        p->color[0] =   0 + (std::rand() % 210);
        p->color[1] =   0 + (std::rand() % 240);
        p->color[2] = 151 + (std::rand() % 105);
      }
      std::cout << "\t\t color " << int(p->color[0]) << " " << int(p->color[1])
                << " " << int(p->color[2]);
      // double max_dist = p->max_dist_allow * 0.5;
      for (int j = 0; j < 3; ++j) {
        if (num_particles > 2) {
          // Set random locations
          p->pos[j] = (std::rand() / (RAND_MAX + 1.0) - 0.9) * kBohrRadiusProton;
        }
        // Increase brightness
        int increase = p->color[j] / divider;
        if (p->color[j] + increase > 255) p->color[j]  = 255;
        else                              p->color[j] += increase;
      }
      std::cout << "\t\t color " << int(p->color[0]) << " " << int(p->color[1])
                << " " << int(p->color[2]) << std::endl;
      std::cout << "\t\t pos " << p->pos[0] << " " << p->pos[1] << " " << p->pos[2] << std::endl;
    }
  }

  int w_to_log_id = 0;  // Rotate through particles to log to screen, when we don't

  // Log a particle and misc info.
  void LogStuff(Particle* w) {
    // Log based on time interval to console.
    if (w_to_log_id == w->id || w->log_count > 0 || w->flipped_dis_vel_dot_prod) {
      auto now = std::chrono::system_clock::now(); // Get current time and see if we logged more than xxx milliseconds ago.
      auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(now - last_log_time);
      if (w->log_count > 0 || duration.count() > 800 || w->flipped_dis_vel_dot_prod) {
        w->SetColorForConsole();
        if (w->flipped_dis_vel_dot_prod) {
          w->log_prev_log_lines(1);
        }
        std::string log_line_str = w->FormatLogLine(num_drawing_event_already,
                                                    num_wait_for_drawing_event);
        num_drawing_event_already  = 0;
        num_wait_for_drawing_event = 0;
        w->tee << log_line_str << std::endl;
        last_log_time = now;
        w_to_log_id = (w_to_log_id + 1) % (num_particles_ / 2);   // Divide by 2 to skip logging protons to screen.
        w->log_count--;
        w->logToBuffer(log_line_str);
        return;
      }
    }
    // Didn't log to screen so consider logging to just file.
    w->ConsiderLoggingToFile(num_drawing_event_already, num_wait_for_drawing_event);
  }


  void CalcEnergyAndLog() {
    // Calculate potential energy.
    double potential_energy_all_w = 0;
    double total_kinetic_e_all_w = 0;
    double closest = 1;  // meters.  Just setting to a large number.

    // Potential energy = sum of all potential energies.
    // Set the potential energy for each pair of particles.
    for (int i = 0; i < num_particles; ++i) {
      Particle* w1 = pars[i];
      total_kinetic_e_all_w += 0.5 * w1->mass_kg * w1->vel_mag2;
      for (int j = i+1; j < num_particles; ++j) {
        Particle* w2 = pars[j];
        assert(w1->dist_mag[j] != 0);
        potential_energy_all_w += kCoulomb * w1->freq_charge * w2->freq_charge / w1->dist_mag[j];
      }
      // If the electron is interesting because it is close to a proton than give it preferential logging.
      if (w1->is_electron && w1->dist_mag_closest < kCloseToTrouble*4
       && w1->dist_mag_closest < closest) {
        closest = w1->dist_mag_closest;
        w_to_log_id = w1->id;
      }
    }
    for (int i = 0; i < num_particles; ++i) {
      // w = wave.  Don't use p because it is confusing with p for potential energy.
      Particle* w = pars[i];
      w->pot_energy_cycle[w->pot_energy_cycle_index] = potential_energy_all_w;
      w->pot_energy_cycle_index = (w->pot_energy_cycle_index + 1) % kPFrequencySubDivisions;
      w->total_kinetic_e_all_w = total_kinetic_e_all_w;
      w->CalcAveragePotentialEnergy();
      LogStuff(w);
    }
  }

  void AllForcesOnParticle(Particle * part_ptr) {
    int part_num = part_ptr->id;
    part_ptr->InitVarsToCalcForces();
    for (int i = 0; i < num_particles; ++i) {
      if (i == part_num) continue;
      part_ptr->CalcForcesFromParticle(pars[i]);
    }
  }

  // Called once for every screen draw.
  void moveParticles() {

    double pos_change_per_particle[kMaxParticles];
    for (int j = 0; j < num_particles; ++j) {
      pos_change_per_particle[j] = 0;
    }
    for (int i=0; i<num_particles; ++i) {
      // Only check for escape once per significant movement to save on CPU.
      pars[i]->CheckForEscape();
    }
    // Do until we get significant movement, then wait for screen draw.
    for (int iter = 0; iter<4096 && !screen_draw_event_occurred; ++iter) {
      for (int i = 0; i < num_particles; ++i) {
        Particle * wave_ptr = pars[i];
        wave_ptr->freq_charge = wave_ptr->GetSinusoidalChargeValue();
        for (int j = 0; j < num_particles; ++j) {
          wave_ptr->dist_calcs_done[j] = false;
        }
      }
      for (int j = 0; j < num_particles; ++j) {
        Particle * part_ptr = pars[j];
        AllForcesOnParticle(part_ptr);
        part_ptr->ApplyForcesToParticle();
        pos_change_per_particle[j] += std::abs(pars[j]->pos_change_magnitude);
      }

      // Find the shortest dt and set the new dt to that.
      dt = pars[0]->new_dt;
      for (int j = 1; j < num_particles; ++j) {
        dt = std::min(pars[j]->new_dt, dt);
      }
      // std::cout << "\t dt " << std::scientific << dt << std::endl;
      CalcEnergyAndLog();
      time_ += dt;
      count++;      // Num iterations of move particles.
      if (count >= 1124831) {
        if (dt > 1) std::cout << ".";
      }

      for (int j = 0; j < num_particles; ++j) {
        // Lets not move faster than it would take an electron to go from center to edge,
        // faster than two seconds.
        const double max_pos_change_desired = kBohrRadius / (60*2);  // 60 fps
        if (pos_change_per_particle[j] > max_pos_change_desired) {
          return;
        }
      }
      // If screen_draw_event_occurred then we are over our computation budget.
    }
  }
};

} // namespace

#pragma clang diagnostic pop