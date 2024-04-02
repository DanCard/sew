#include "particle.h"

#include <cassert>
#include <cmath>   // For M_PI constant and other match functions such as round.
#include <iostream>
#include <string>
// #include <thread>

#include "atom.h"
#include "logger.h"

namespace sew {

Particle::Particle(int id, bool is_electron, Atom* a,
                    double mass_mev, double mass_kg, double avg_q, double q_amplitude,
                    double max_dist_allowed, double max_speed_allowed,
                    Logger *logger
                    ) :
          id(id), is_electron(is_electron), a_(a),
          mass_mev(mass_mev), mass_kg(mass_kg), avg_q(avg_q), q_amplitude(q_amplitude),
          // Want paired particles to be at opposite ends of the sine wave.
          // In other words to be a half cycle away.
          // This causes one electron to be at -2e charge, while the other is at 0 charge.
          initial_charge(id%2 == 0 ? 0 : M_PI),
          max_dist_allow(max_dist_allowed),
          max_speed_allow(max_speed_allowed),
          logger(logger),
          tee_logger(is_electron ? "e" + std::to_string(id) + ".log"
                                 : "p" + std::to_string(id) + ".log"),
                     tee(tee_logger.get_stream())
  {
    log_count = is_electron ? 2 : 0;
    std::cout << "\t" << (is_electron ? "electron" : "proton") << " " << id
      << "\t frequency " << frequency << "  "
              << (is_electron ? "electron" : "proton") << " mass mev " << mass_mev << std::endl;
    tee << "\t\t initial_charge " << initial_charge << std::endl;
  }

void Particle::logToBuffer(const std::string &s) {
    prev_log_lines[prev_log_line_index] = s;
    prev_log_line_index = (prev_log_line_index + 1) % kMaxPrevLogLines;
  }

void Particle::log_prev_log_lines(int max_lines_to_log) {
    Logger::SetColorForConsole(color[0], color[1], color[2]);
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


void Particle::ConsiderLoggingToFile(int count) {
    double danger_speed = max_speed_allow / 2;

    bool particles_are_close      = dist_mag_closest < (kCloseToTrouble * 2);
    bool particles_are_close_very = dist_mag_closest < kCloseToTrouble;
    bool fast_fraction_changed_significantly = (prev_fast_fraction / fast_fraction > 1.1) &&
            (prev_fast_fraction - fast_fraction > 0.1);
    // Make it difficult for log files to reach gigabyte sizes on long running simulation.
    static int local_count = 0;
    static int ll = 1;  // Limit logging to once every x lines.
    if (++local_count > 20000000  && ll <= 16384) {
      local_count = 0;
      ll *= 2;
      // std::cout << "\t\t\t\t ll " << ll << std::endl;
    }
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
      std::string log_line_str = logger->FormatLogLine(this);
      log_file << log_line_str << log_source << std::endl;
      prev_fast_fraction = fast_fraction;
      log_count--;
      logToBuffer(log_line_str);
    }
  }

  void Particle::HandleEscape() {
    log_prev_log_lines(1);
    tee << "\t" << (is_electron ? "electron" : "proton")
        << " escaped past radius of " << max_dist_allow
        << "  Dist mag from origin : " << distance_mag_from_origin
        << "  Current velocity " << sqrt(vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2])
        << "  x " << vel[0] << " y " << vel[1] << " z " << vel[2]
        << "  total energy " << a_->total_energy;
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
    log_count = 1;  // Force logging around this event.
    tee << std::endl;
    if (a_->total_energy_cap == 0) {
      a_->total_energy_cap = a_->total_energy * 1.2;  // * 1.1 because protons continue gaining energy.
    }
  }


  void Particle::CheckForEscape() {
    distance_mag_from_origin = sqrt(pow(pos[0], 2) + pow(pos[1], 2) + pow(pos[2], 2));
    // If particle escapes, then zero out the velocity.
    // This is a hack to limit the problem of energy gain.
    // https://en.wikipedia.org/wiki/Energy_drift
    if (distance_mag_from_origin > max_dist_allow) {
      HandleEscape();
    }
  }


  // Get the charge that is based on freq and thus time.
  double Particle::ChargeSinusoidal() const {
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
    return avg_q + (q_amplitude * sin((frequency * a_->time_ * 2 * M_PI) + initial_charge));
  }

  void Particle::InitVarsToCalcForces() {
    dist_mag_closest  = 1;      // 1 represents a large number that will be overwritten
                                // quickly when checking for closest.
    par_closest = nullptr;
    energy_dissipated = false;
    for (double & force : forces) {
      force = 0;
    }
    for (double & magnet_force : magnet_fs) {
      magnet_force = 0;
    }
    for (int i=0; i<a_->num_particles; ++i) {
      dist_mag [i] = 0;       // Will assert if we test and find these values.
      dist_mag2[i] = 0;
    }
  }

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
  void Particle::
  MagneticForce(Particle *oth, double e_field_magnitude_oth, double q2,
                         double *dist_vector, double dist_mag_2, double *dist_unit_vector) {
    // Save some CPU processing time when we are in slow mode, by ignoring the insignificant magnetic force.
    if (a_->dt <= kShortDt) return;
    // Constants for calculating magnetic force.
    // https://academic.mu.edu/phys/matthysd/web004/l0220.htm
    // Permeability of free space.  https://en.wikipedia.org/wiki/Vacuum_permeability
    const double kPermeabilityDividedBy4Pi = 1e-7;  // T * m / A
    // Biot-Savart law.  https://en.wikipedia.org/wiki/Biot%E2%80%93Savart_law
    // B = μ0/4π * q * (cross product of v and r) / |r|^2
    double vel_oth_cross_dis[3];
    cross(oth->vel, dist_vector, vel_oth_cross_dis);
    double b_field_by_oth_vel[3];    // Magnetic field created by other particles velocity
    for (int i = 0; i < 3; ++i) {
      // https://en.wikipedia.org/wiki/Biot%E2%80%93Savart_law
      // https://docs.google.com/document/d/1Tc1AZltK30YA-u-corns3hVO7Zcx22uYE1TyT4i_NNU
      b_field_by_oth_vel[i] = kPermeabilityDividedBy4Pi * q2 * vel_oth_cross_dis[i] / dist_mag_2;
      assert(!std::isnan(b_field_by_oth_vel[i]));
    }
    // Lorentz force from oth particle magnetic field = q * v x B
    double v_cross_b_field_oth_v[3];
    cross(vel, b_field_by_oth_vel, v_cross_b_field_oth_v);
    for (int i = 0; i < 3; ++i) {
      b_force_by_oth_vel[i] = freq_charge * v_cross_b_field_oth_v[i];
      magnet_fs[i] += b_force_by_oth_vel[i];   // Only used for logging.
    }
    // Above we calculated the magnetic force due to moving charged particle.
    // Below we will calculate the magnetic force due to the other particle's intrinsic magnetic field.

    // b = e / c
    // http://hyperphysics.phy-astr.gsu.edu/hbase/Waves/emwv.html
    // https://www.se.edu/kfrinkle/wp-content/uploads/sites/89/2013/10/main23v130303.pdf
    // Deriving E=cB and B=µ0·ε0·c·E using Maxwell's Equations: https://www.youtube.com/watch?v=dW1q8XQcEdI
    // Magnetic field created by other.
    double b_field_caused_by_intrinsic_oth_magnitude = e_field_magnitude_oth / kC;
    double b_field_by_oth_intrinsic[3];
    for (int i = 0; i < 3; ++i) {
      // (i+1)%3 = rotate 90 degrees.  To do: Figure out correct way of doing this.
      b_field_by_oth_intrinsic[(i+1)%3] = dist_unit_vector[i] * b_field_caused_by_intrinsic_oth_magnitude;
    }
    double vel_cross_b_field_intrinsic[3];
    cross(vel, b_field_by_oth_intrinsic, vel_cross_b_field_intrinsic);
    for (int i = 0; i < 3; ++i) {
      // Lorentz force from oth particle magnetic field = q * v x B
      b_f_intrinsic[i] = freq_charge * vel_cross_b_field_intrinsic[i];
    }
    for (int i = 0; i < 3; ++i) {
      magnet_fs[i] += b_f_intrinsic[i];
    }
    /*
    tee << " magnetic field q2 "
              << b_field_by_oth_vel[0] << ' ' << b_field_by_oth_vel[1] << ' ' << b_field_by_oth_vel[2];
    tee << " magnetic force from q2 " << magnet_f_oth[0] << ' ' << magnet_f_oth[1] << ' ' << magnet_f_oth[2];
    tee << "\t background f " << m_bg_f[0] << ' ' << m_bg_f[1] << ' ' << m_bg_f[2] << std::endl;
    */
  }

  void Particle::CalcForcesFromParticle(Particle* oth /* other particle */) {
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
    MagneticForce(oth, e_field_magnitude_oth, other_charge, dist, dist_magn2, dist_unit_vector);

    for (int i = 0; i < 3; ++i) {
      forces[i] += eforce[i] + magnet_fs[i];
    }
    bool forces_attract = is_electron != oth->is_electron;
    // Determine which attractive particle we are closest to.
    // This is often the particle exerting the largest force.
    // Used for logging.
    if (dist_magn < dist_mag_closest && (forces_attract || !is_electron)) {
      par_closest    = oth;
      par_closest_id = oth_id;
      dist_mag_closest = dist_magn;
      for (int i = 0; i < 3; ++i) {
        dist_closest[i] = dist[i];
      }
      force_mag_closest = force_magnitude;
    }
  }


  // When forces are very high we should slow down delta time (dt) to make calculations more accurate.
  void Particle::NewDt(double * fast_fraction_ptr) {
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
    // Dissipate energy when heading away.
    if (dist_vel_dot_prod >= 0
     && orig_vel_dot_prod >= 0
     && dist_mag_closest > (kCloseToTrouble*16)) {
      if (dis_vel_dot_prod_old < 0) {
        flipped_dis_vel_dot_prod = true;  // If true then log
      } else if (a_->total_energy_cap != 0
              && a_->total_energy_cap < a_->total_energy) {
        energy_dissipated = true;
        for (int i = 0; i < 3; ++i) {
          // Will cause total energy to drop.
          vel[i] *= 0.999999;  // Combat energy gain.  Dissipate energy when heading away.
        }
      }
    }
  }

  // With little or no space between particles, forces approach infinity.
  // If the simulation won't work because of large errors because of huge forces.
  void Particle::TeleportIfTooCloseToProton() {
    if (!is_electron) return;
    // Higher number = faster particle gets ejected.
    // Lower number = more likely to get teleported, loop around proton, get closer.
    bool is_force_too_high = force_magnitude > kForceTooHigh;
    if (!is_force_too_high) return;
    if (vel_mag < max_speed_allow) return;

    // Force needs to be in the same direction as velocity.
    // If not then we don't to worry about particles heading in opposite directions.
    // Need unit vector for velocity and force, then do dot product.
    Particle* close = par_closest;
    // If the vectors are pointing in the same direction (aligned), their dot product is 1.
    // If they are perpendicular, it's 0. If they are pointing in opposite directions, it's -1.
    // Since we are using distance as a proxy for force, sign is reversed.
    bool force_same_dir_as_v = dist_vel_dot_prod < -0.75;    // -0.75 = guess
    if (!force_same_dir_as_v) return;

    close->log_prev_log_lines(1);
           log_prev_log_lines(1);
    tee << "\t Teleporting because force too high. " << (is_electron ? "e" : "p") << id
        << "  distance velocity dot product " << dist_vel_dot_prod
        << "  force magnitude " << force_magnitude
        << " x " << forces[0] << " y " << forces[1] << " z " << forces[2] << std::endl;
    tee << " velocity unit vector: x " << vel_unit_vec[0]
                              << " y " << vel_unit_vec[1]
                              << " z " << vel_unit_vec[2]
        << " distance unit vector: x " << dist_unit_vec[0]
                              << " y " << dist_unit_vec[1]
                              << " z " << dist_unit_vec[2] << std::endl;
    tee << "  Dist " << dist_mag_closest << " x " << dist_closest[0] << " y " << dist_closest[1] << " z " << dist_closest[2]
        << "  Teleporting to other side of "
        << (close->is_electron ? "electron" : "proton")
        << ".  Position " << pos[0] << " " << pos[1] << " " << pos[2]
        << "  proton pos x " << close->pos[0] << " y " << close->pos[1] << " z " << close->pos[2]
        << std::endl;
    // Teleport the particle to other side of the close by particle.
    for (int i = 0; i < 3; ++i) {
      pos[i] += -2 * dist_closest[i];
    }
    tee << "\t Teleported to " << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
    v_when_teleported = vel_mag;
    log_count = 1;               // Force extra logging after this event.
    close->log_file << "\t\t\t" << (is_electron ? "e" : "p") << id << " teleported "  << std::endl;
    // if (num_particles_ > 2)
    close->log_count = 1;          // Force extra logging for other particle.
  }

  void Particle::ApplyForcesToParticle() {
    force_magnitude = sqrt(pow(forces[0], 2) + pow(forces[1], 2) + pow(forces[2], 2));
    for (int i = 0; i < 3; ++i) {
      acceleration[i] = forces[i] / this->mass_kg;
      vel[i] += acceleration[i] * a_->dt;
    }
    vel_mag2 = pow(vel[0], 2) + pow(vel[1], 2) + pow(vel[2], 2);
    vel_mag  = sqrt(vel_mag2);
    for (int i = 0; i < 3; ++i) {
      pos_change[i]  = vel[i] * a_->dt;
      if (!kHoldingProtonSteady || is_electron) {
        pos        [i] += pos_change[i];
      }
      vel_unit_vec [i]  = vel[i] / vel_mag;
      dist_unit_vec[i]  = dist_closest[i] / dist_mag_closest;
    }
    dis_vel_dot_prod_old = dist_vel_dot_prod;
    flipped_dis_vel_dot_prod = false;

    if (is_electron) {
      // dot product between dist of closest and velocity.
      dist_vel_dot_prod = dist_unit_vec[0] * vel_unit_vec[0] +
                         dist_unit_vec[1] * vel_unit_vec[1] +
                         dist_unit_vec[2] * vel_unit_vec[2];
    }
    // Calc position unit vector
    pos_magnitude = sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
    for (int i = 0; i < 3; ++i) {
      pos_unit_vec[i] = pos[i] / pos_magnitude;
    }
    orig_vel_dot_prod = pos_unit_vec[0] * vel_unit_vec[0] +
                        pos_unit_vec[1] * vel_unit_vec[1] +
                        pos_unit_vec[2] * vel_unit_vec[2];
    // If it is a proton don't care if it is heading away or towards electrons,
    // Since don't have ejection issue with protons.
    if (!is_electron) {
      // Only consider if proton is heading away or from origin when considering to dissipate energy.
      dist_vel_dot_prod = orig_vel_dot_prod;
    }
    // Use dot product between velocity and distance of closest to determine
    // if coming or going relative to closest attracting particle.
    // If the vectors are pointing in the same direction (aligned), their dot product is 1.
    // If they are perpendicular, it's 0. If they are pointing in opposite directions, it's -1.
    NewDt(&fast_fraction);

    // pos_change_magnitude used to decide if we have moved enough for current frame.
    pos_change_magnitude = sqrt(pow(pos_change[0], 2) + pow(pos_change[1], 2) + pow(pos_change[2], 2));
    // When particles are on top of each other forces approach infinity.
    // To get around that problem we teleport to other side of particle.
    // Alternative approach is to have a preprogrammed path and not bother with force calcs.
    if (is_electron) {
      TeleportIfTooCloseToProton();
    }
  }

} // namespace