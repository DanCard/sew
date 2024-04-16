#include "logger.h"

#include <cmath>   // For M_PI constant and other match functions such as round.
#include <chrono>  // For logging every 500 ms
#include <iomanip>
#include <sstream>
#include <string>

#include "atom.h"
#include "particle.h"

namespace sew {

  Logger::Logger(Atom* a) :
    a_(a), start_time_(std::chrono::_V2::system_clock::now()) {
  }

  // A bit of a mess because we have particle data and particles(system) data that we are logging.
  std::string Logger::FormatLogLine(const Particle* w, bool to_file) const {
    Particle* par_closest = w->par_closest;
    SFloat charge_of_closest = par_closest->freq_charge;

    std::ostringstream log_line;
    log_line
      << std::setw(8) << a_->count << (w->is_electron ? " e" : " p") << w->id
      << "⋅" << (par_closest->is_electron ? "e" : "p") << w->par_closest_id
      << std::scientific << std::setprecision(3)
      << "  dis"  << std::setw(10) << (w->is_electron ? w->dist_mag_closest : w->pos_magnitude);
    if (velocity_logging) {
      log_line << "  vel " << std::setw(10) << w->vel_mag;
    }
    if (velocity_component_logging) log_line << (velocity_logging ? Log3dArray(w->vel, "v") : "");
    if (percent_energy_dissipated_logging)
         log_line << std::setw(10) << std::setprecision(3) << w->percent_energy_dissipated << '%';
    else log_line << (w->energy_dissipated ? " *" : "  ");
    assert(w->energy_dissipated ? w->percent_energy_dissipated > 0 : w->percent_energy_dissipated == 0);
    if (dv_logging) {
      log_line << "d⋅v " << std::setw( 4) << std::fixed
               << w->dist_vel_dot_prod    // -1 = approaching, 1 = leaving
               << std::scientific;
    }
    if (dt_logging) {
      log_line << "  dt" << std::setw( 8) << std::setprecision(1)
               << a_->dt;    // << " new " << new_dt
    }
    if (fast_logging) {
      log_line << " fast"  << std::setw(10) << std::setprecision(3)
               << w->fast_fraction;
    }
    log_line
      << std::setprecision(2)
      << "  f"    << std::setw(9) << w->force_mag_closest
   // << Log3dArray(forces  , " fs")
   // << " Bv " << sqrt(pow(b_force_by_oth_vel[0], 2) + pow(b_force_by_oth_vel[1], 2) + pow(b_force_by_oth_vel[2], 2))
   // << " Bi " << std::setw( 9) << sqrt(pow(b_f_intrinsic[0], 2) + pow(b_f_intrinsic[1], 2) + pow(b_f_intrinsic[2], 2))
   // << Log3dArray(b_f_intrinsic, "Bi")
   // << Log3dArray(magnet_fs, "tB"  )
   // << "  B "   << sqrt(magnet_fs[0]*magnet_fs[0] + magnet_fs[1]*magnet_fs[1] + magnet_fs[2]*magnet_fs[2])
   // << Log3dArray(acceleration, "a")
   // << Log3dArray(pos_change, "chng") << std::setprecision(1)
      << (position_logging ? Log3dArray(w->pos, "pos") : "");
   // << " min pos change " << min_pos_change_desired
   // << round(fast_fraction * 10) * 10 << '%'
    if (charge_logging) {
      log_line    << std::setw( 6) << std::setprecision(1) << std::fixed
      << " chrg"  << std::setw( 4) << static_cast<int>(std::round((w->freq_charge / w->avg_q) * 100.0f)) << '%'
      << " oth"   << std::setw( 4)
      << static_cast<int>(std::round((charge_of_closest / par_closest->q_amplitude) * 100.0f)) << '%';
    }
   // << " inv"   << std::setw(12) << inverse_exponential
    if (energy_logging) {
      log_line << std::scientific << std::setprecision(2)
        // P energy goes negative, that is why width is larger.
        << "  pe" << std::setw(10) << a_->potential_energy_average
        << " ke"  << std::setw( 9) << a_->total_kinetic_energy
        << " te"  << std::setw( 9) << a_->total_energy;
    }
    if (frame_draw_statistics_logging) {
      log_line
              // Late.  Drawing event already occurred.
              << " L " << std::setw(2)
              << a_->n_times_per_screen_log_MoveParticles_not_compl_before_next_frame_draw_event
              // Early. Waited on drawing event.
              << " E " << std::setw(2)
              << a_->n_times_per_screen_log_MoveParticles_completed_before_next_frame_draw_event;
    }
    if (iterations_logging) {
      log_line << " it" << std::setw(5) << a_->iter;
    }
    if (wall_clock_time_logging) {
      auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(
              std::chrono::_V2::system_clock::now() - start_time_);
      log_line << " mi " << std::setw(6) << elapsed_time.count();
    }
    if (time_logging) {
      log_line << " t " << std::scientific << std::setprecision(3) << std::setw(9) << a_->time_;
    }
    log_line << " chng"  << std::setw(8) << std::setprecision(1) << w->pos_change_magnitude
             << " dtsltu " << std::setprecision(1) << w->dist_traveled_since_last_trail_update;
    /*
    for (int i=0; i<num_particles_; ++i) {
      if (i == id || i == par_closest_id) continue;
      log_line << (i >= (num_particles_/2) ? " p" : "  e") << i << " " << dist_mag_all[i];
    }
    */
    return log_line.str();
  }

  // Log a particle and misc info.
  void Logger::LogStuff(Particle* w) {
    static int milliseconds_to_wait_before_logging = 400;
    static std::chrono::_V2::system_clock::time_point last_log_time;
    /*
    bool dist_reset = w->prev_dist_traveled_since_last_trail_update
                        >  w->dist_traveled_since_last_trail_update;
    */

    // Log based on time interval to console.
    auto now = std::chrono::system_clock::now();
    bool do_log =  w->log_count > 0
     //        || (w->energy_dissipated && !w->energy_dissipated_prev)
     //        ||  w->flipped_dis_vel_dot_prod
               ;
    if (do_log
     // || (dist_reset && w->id == 0)
     || (w->id == w_to_log_id &&
        std::chrono::duration_cast<std::chrono::milliseconds>(
        now - last_log_time).count() > milliseconds_to_wait_before_logging)) {
      SetColorForConsole(w->color[0], w->color[1], w->color[2]);
      // std::cout << std::endl;
      std::string log_line_str = FormatLogLine(w, false);
      a_->n_times_per_screen_log_MoveParticles_not_compl_before_next_frame_draw_event = 0;
      a_->n_times_per_screen_log_MoveParticles_completed_before_next_frame_draw_event = 0;
      w->tee << log_line_str << std::endl;
      last_log_time = now;
      w_to_log_id = (w_to_log_id + 1) % a_->num_particles;
      w->log_count--;
      w->logToBuffer(log_line_str);
      if (milliseconds_to_wait_before_logging < 1200)
        milliseconds_to_wait_before_logging += 50;
      return;
    }
    // Didn't log to screen so consider logging to just file.
    w->ConsiderLoggingToFile(a_->count);
  }

  void Logger::ChargeLoggingToggle() {charge_logging = !charge_logging;  }
  void Logger::DtLoggingToggle    () {dt_logging     = !dt_logging    ;  }
  void Logger::EnergyLoggingToggle() {energy_logging = !energy_logging;  }
  void Logger::FastLoggingToggle  () {fast_logging   = !fast_logging  ;  }
  void Logger::DvModeToggle       () {dv_logging     = !dv_logging    ;  }
  void Logger::FrameDrawStatisticsLogToggle() {
    frame_draw_statistics_logging = !frame_draw_statistics_logging;
  }
  void Logger::IterationsLoggingToggle() {
    iterations_logging = !iterations_logging;
  }
  void Logger::PositionLoggingToggle() {position_logging = !position_logging; }
  void Logger::PercentEnergyDissipatedToggle() {
    percent_energy_dissipated_logging = !percent_energy_dissipated_logging;
  }
  void Logger::VelocityComponentsLogToggle() {velocity_component_logging = !velocity_component_logging;}
  void Logger::VelocityLoggingToggle() {velocity_logging = !velocity_logging; }
  void Logger::TimeLoggingToggle    () {time_logging     = !time_logging    ; }
  void Logger::WallClockLoggingToggle() {
    wall_clock_time_logging = !wall_clock_time_logging;
  }



} // namespace
