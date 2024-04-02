#include "logger.h"

#include <cmath>   // For M_PI constant and other match functions such as round.
#include <chrono>  // For logging every 500 ms
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>

#include "atom.h"
#include "particle.h"

namespace sew {

  Logger::Logger(Atom* a) : a_(a) {}

  void Logger::EnergyLoggingToggle() {
    energy_logging = !energy_logging;
  }

  void Logger::VelocityLoggingToggle() {
    velocity_logging = !velocity_logging;
  }

  // A bit of a mess because we have particle data and particles(system) data that we are logging.
  std::string Logger::FormatLogLine(Particle* w) {
    Particle* par_closest = w->par_closest;
    double charge_of_closest = par_closest->freq_charge;

    std::ostringstream log_line;
    log_line
      << std::setw(8) << a_->count << (w->is_electron ? " e" : " p") << w->id
      << "⋅" << (par_closest->is_electron ? "e" : "p") << w->par_closest_id
      << std::scientific << std::setprecision(3)
      << "  dis"  << std::setw(10) << (w->is_electron ? w->dist_mag_closest : w->pos_magnitude)
      << "  vel " << std::setw(10) << w->vel_mag
      << (velocity_logging ? Log3dArray(w->vel, "v") : "")
      << (w->energy_dissipated ? " *" : "  ") << std::fixed << std::setprecision(1)
      << "d⋅v " << std::setw( 4) << w->dist_vel_dot_prod    // -1 = approaching, 1 = leaving
      << std::scientific << std::setprecision(2)
      << "  dt"   << std::setw( 9) << a_->dt    // << " new " << new_dt
      << " fast"  << std::setw(10) << std::setprecision(3) << w->fast_fraction << std::setprecision(2)
      << " f "    << std::setw( 9) << w->force_mag_closest
   // << Log3dArray(forces  , " fs")
   // << " Bv " << sqrt(pow(b_force_by_oth_vel[0], 2) + pow(b_force_by_oth_vel[1], 2) + pow(b_force_by_oth_vel[2], 2))
   // << " Bi " << std::setw( 9) << sqrt(pow(b_f_intrinsic[0], 2) + pow(b_f_intrinsic[1], 2) + pow(b_f_intrinsic[2], 2))
   // << Log3dArray(b_f_intrinsic, "Bi")
   // << Log3dArray(magnet_fs, "tB"  )
   // << "  B "   << sqrt(magnet_fs[0]*magnet_fs[0] + magnet_fs[1]*magnet_fs[1] + magnet_fs[2]*magnet_fs[2])
   // << Log3dArray(acceleration, "a")
   // << " chng"  << std::setw(10) << std::setprecision(3) << pos_change_magnitude
   // << Log3dArray(pos_change, "chng") << std::setprecision(1)
   // << Log3dArray(pos     , "pos")
   // << " min pos change " << min_pos_change_desired
   // << round(fast_fraction * 10) * 10 << '%'
                  << std::setw( 6) << std::setprecision(1) << std::fixed 
      << " chrg"  << std::setw( 4) << int(round((w->freq_charge/w->avg_q)*100)) << '%'
      << " oth"   << std::setw( 4)
      << int(round((charge_of_closest / par_closest->q_amplitude) * 100)) << '%';
   // << " inv"   << std::setw(12) << inverse_exponential
    if (energy_logging) {
      log_line << std::scientific << std::setprecision(2)
        // P energy goes negative, that is why width is larger.
        << "  pe"   << std::setw(10) << a_->potential_energy_average
        << " ke"    << std::setw( 9) << a_->total_kinetic_energy
        << " te"    << std::setw( 9) << a_->total_energy;
    }
    log_line
      << " L "    << std::setw( 2) << a_->num_drawing_event_already   // Late.  Drawing event already occurred.
      << " E "    << std::setw( 2) << a_->num_wait_for_drawing_event  // Calcs were early.  Waited on drawing event.
    ;
    /*
    for (int i=0; i<num_particles_; ++i) {
      if (i == id || i == par_closest_id) continue;
      log_line << (i >= (num_particles_/2) ? " p" : "  e") << i << " " << dist_mag[i];
    }
    */
    return log_line.str();
  }

  // Log a particle and misc info.
  void Logger::LogStuff(Particle* w) {
    static std::chrono::_V2::system_clock::time_point last_log_time;

    // Log based on time interval to console.
    auto now = std::chrono::system_clock::now();
    bool do_log =  w->log_count > 0
     //        || (w->energy_dissipated && !w->energy_dissipated_prev)
               ||  w->flipped_dis_vel_dot_prod;
    if (do_log || (w->id == w_to_log_id &&
        std::chrono::duration_cast<std::chrono::milliseconds>(
        now - last_log_time).count() > 1000)) {
      SetColorForConsole(w->color[0], w->color[1], w->color[2]);
      if (w->flipped_dis_vel_dot_prod) {
        w->log_prev_log_lines(1);
      }
      std::string log_line_str = FormatLogLine(w);
      a_->num_drawing_event_already  = 0;
      a_->num_wait_for_drawing_event = 0;
      w->tee << log_line_str
             << (w->energy_dissipated_prev != w->energy_dissipated ? " *" : "  ")
             << std::endl;
      last_log_time = now;
      w_to_log_id = (w_to_log_id + 1) % a_->num_particles;
      w->log_count--;
      w->logToBuffer(log_line_str);
      w->energy_dissipated_prev = w->energy_dissipated;
      return;
    }
    // Didn't log to screen so consider logging to just file.
    w->ConsiderLoggingToFile(a_->count);
  }

} // namespace