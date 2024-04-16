#include "atom.h"

#include <cstdlib> // For rand()
#include <iostream>
#include <thread>
#include <vector>

#include "logger.h"
#include "particle.h"
#include "thread_pool.h"


namespace sew {

Atom::Atom(int numParticles) :
         num_particles(numParticles), long_dt(kLongDt), short_dt(kShortDtSlow) {
    if (num_particles == 0) return;
    std::cout << "\t Electron period: " << kEPeriod << " proton period: " << kPPeriod << std::endl;
    logger = new sew::Logger(this);
    std::cout << "  kMaxPosChangeDesiredPerFrame " << kMaxPosChangeDesiredPerFrame
      << "  kBhorRadius " << kBohrRadius << "  kBohrRadiusProton " << kBohrRadiusProton
      << "  kLithiumAtomSize " << kLithiumAtomSize
      << std::endl;
    std::cout << "  max speed electron " << kMaxSpeedElectron
      << "  kPFrequencySubDivisions " << kPFrequencySubDivisions
      << "  kEFrequencySubDivisions " << kEFrequencySubDivisions
              << "\t kShortDt " << kShortDtSlow << "  kLongDt " << kLongDt << std::endl;
    std::cout << "\t\t\t\t kEFrequency " << kEFrequency << "  kPFrequency " << kPFrequency << std::endl;
    // std::cout << "\t\t kBohrMagneton " << kBohrMagneton << "  kProtonMagneticMoment " << kProtonMagneticMoment << std::endl;
    const SFloat nucleus_initial_radius = kBohrRadiusProton * (num_particles / 2);
    const SFloat electron_initial_radius = nucleus_initial_radius * 4;
    std::cout << "\t\t\t\t nucleus initial radius " << nucleus_initial_radius << "  electron initial radius " << electron_initial_radius << std::endl;
    int divider;  // Prefer bright colors, but with many particles becomes indistinguishable.
         if (num_particles <= 2)  divider = 1;
    else if (num_particles <= 4)  divider = 3;
    else if (num_particles <= 6)  divider = 4;
    else                          divider = 8;  // With more particles don't brighten as much.
    bool is_electron = true;
    // Create the subatomic particles.
    for (int i = 0; i < numParticles; ++i) {
      Particle* p;
      if (i < numParticles / 2) {
        pars[i] = new Electron(i, this, logger, numParticles > 8 ? kLithiumAtomSize : kBohrRadius * 2);
        p = pars[i];
        // Set pseudo random colors.  Electrons tend to be more red.  Protons tend to be more blue.
        // Prefer bright colors over dark colors.
        // p->color[1] =  51 + (std::rand() % 155);
        p->color[0] = std::rand() % 256;
        p->color[1] = std::rand() % 256;
        p->color[2] = std::rand() % 245;
        if (i == 0) {
          // p->vel[1] = -1e4;
          // p->pos[0] = -p->max_dist_allow * 0.95f;
          p->color[0] = 255;
          p->color[1] = 120;
          p->color[2] = 120;
        }
        for (int j = 0; j < 3; ++j)
            p->pos[j] = ((float)(std::rand() / (RAND_MAX + 1.0)) - 0.5f) * electron_initial_radius;
      } else {
        is_electron = false;
        pars[i] = new Proton(i, this, logger, nucleus_initial_radius * 2);
        p = pars[i];
        p->color[0] =   0 + (std::rand() % 210);
        p->color[1] =   0 + (std::rand() % 240);
        p->color[2] = 151 + (std::rand() % 105);
        p->vel[1] = 1e4;
        if (i == 1) {
          p->pos[0] = 0;
        }
        for (int j = 0; j < 3; ++j)
            p->pos[j] = ((float)(std::rand() / (RAND_MAX + 1.0)) - 0.5f) * nucleus_initial_radius;
      }
      std::cout << "\t\t color " << int(p->color[0]) << " " << int(p->color[1])
                << " " << int(p->color[2]);
      // SFloat max_dist = p->max_dist_allow * 0.5;
      for (int j = 0; j < 3; ++j) {
        // Increase brightness
        int increase = p->color[j] / divider;
        if (p->color[j] + increase > 255) p->color[j]  = 255;
        else                              p->color[j] += increase;
      }
      /*
      std::cout << "\t\t color " << int(p->color[0]) << " " << int(p->color[1])
                << " "           << int(p->color[2]) << std::endl;
      */
      std::cout << "\t\t pos " << p->pos[0] << " " << p->pos[1] << " " << p->pos[2] << std::endl;
    }
    for (SFloat & p_energy_cycle_ : pot_energy_cycle) {
      p_energy_cycle_ = 0;
    }
    const int num_threads = std::min((int)std::thread::hardware_concurrency(), num_particles*2);
    std::cout << "\t\t\t num threads " << num_threads << std::endl;
    thread_pool = new ThreadPool(num_threads);
}

  // Potential energy changes because of sinusoidal charge frequency.
void Atom::CalcAveragePotentialEnergy() {
    sum_p_energy += total_potential_energy;
    pot_energy_cycle[pot_energy_cycle_index] = total_potential_energy;
    pot_energy_cycle_index = (pot_energy_cycle_index + 1) % kEFrequencySubDivisions;
    // Subtract the potential energy that is rolling off our average.
    sum_p_energy -= pot_energy_cycle[pot_energy_cycle_index];
    // -1 because of the value subtracted above.
    potential_energy_average = sum_p_energy / (kEFrequencySubDivisions - 1);
  }

void Atom::CalcEnergyOfAtom() {
    // Calculate potential energy.
    total_potential_energy = 0;
    total_kinetic_energy = 0;
    SFloat closest = 1;  // meters.  Just setting to a large number.

    // Potential energy = sum of all potential energies.
    // Set the potential energy for each pair of particles.
    for (int i = 0; i < num_particles; ++i) {
      Particle *w1 = pars[i];   // Wave 1
      total_kinetic_energy += 0.5f * w1->mass_kg * w1->vel_mag2;
      for (int j = i + 1; j < num_particles; ++j) {
        Particle *w2 = pars[j];   // Wave 2
        if(w1->dist_mag_all[j] == 0) {
          std::cout << "\t\t w1->dist_mag_all[j] " << w1->dist_mag_all[j] << std::endl;
        }
        total_potential_energy += kCoulomb * w1->freq_charge * w2->freq_charge / w1->dist_mag_all[j];
      }
      // If the electron is interesting because it is close to a proton than give it preferential logging.
      if (w1->is_electron && w1->dist_mag_closest < kCloseToTrouble
          && w1->dist_mag_closest < closest) {
        closest = w1->dist_mag_closest;
        logger->w_to_log_id = w1->id;
      }
    }

    // Potential energy changes because of sinusoidal charge frequency.
    // Alternatively could just use average charge.
    CalcAveragePotentialEnergy();
    total_energy = potential_energy_average + total_kinetic_energy;
  }

void Atom::AllForcesOnParticle(Particle * part_ptr) {
    int part_num = part_ptr->id;
    for (int i = 0; i < num_particles; ++i) {
      if (i == part_num) continue;
      part_ptr->CalcForcesFromParticle(pars[i]);
    }
  }

void Atom::CalcForcesAndApply(Particle * part_ptr) {
  Atom::AllForcesOnParticle(part_ptr);
  part_ptr->ApplyForces();
}

// Twice as slow to use a thread pool.  Why?
// #define use_thread_pool

  // Called once for every screen draw.
void Atom::MoveParticles() {
    SFloat pos_change_per_particle[kMaxParticles];
    for (int j = 0; j < num_particles; ++j) {
      pos_change_per_particle[j] = 0;
    }
    for (int i=0; i<num_particles; ++i) {
      // Only check for escape once per significant movement to save on CPU.
      pars[i]->CheckForEscape();
    }
    // Do until we get significant movement, then wait for screen draw.
    // If frame_draw_event_occurred then we are over our computation budget.
    for (iter = 0; iter<(4096*2) && !frame_draw_event_occurred; ++iter) {
      // Initialize thread unsafe variables.
      for (int i = 0; i < num_particles; ++i) {
        pars[i]->InitVarsToCalcForces();
      }

#ifdef use_thread_pool
      std::future<void> future_results[kMaxParticles];
#endif
      for (int j = 0; j < num_particles; ++j) {
        Particle * part_ptr = pars[j];
#ifdef use_thread_pool
        future_results[j] = thread_pool->AddTask(&Atom::CalcForcesAndApply, this, part_ptr);
#else
        CalcForcesAndApply(part_ptr);
#endif
      }

      dt = pars[0]->new_dt;
      for (int j = 0; j < num_particles; ++j) {
#ifdef use_thread_pool
        future_results[j].get();      // Wait for task in thread pool to finish.
#endif
        Particle * part_ptr = pars[j];
        pos_change_per_particle[j] += std::abs(part_ptr->pos_change_magnitude);
        // Find the shortest dt and set the new dt to that.
        dt = std::min(part_ptr->new_dt, dt);
      }

      CalcEnergyOfAtom();

      for (int i = 0; i < num_particles; ++i) {
        logger->LogStuff(pars[i]);
      }
      time_ += dt;
      count++;      // Num iterations of move particles.

      for (int j = 0; j < num_particles; ++j) {
        if (pos_change_per_particle[j] > kMaxPosChangeDesiredPerFrame) {
          return;
        }
      }
    }
  }

  void Atom::ChargeLoggingToggle() {logger->ChargeLoggingToggle();  }
  void Atom::DtLoggingToggle    () {logger->DtLoggingToggle    ();  }
  void Atom::EnergyLoggingToggle() {logger->EnergyLoggingToggle();  }
  void Atom::FastLoggingToggle  () {logger->FastLoggingToggle  ();  }
  void Atom::DvModeToggle       () {logger->DvModeToggle       ();  }  // g key
  void Atom::IterationsLoggingToggle() {logger->IterationsLoggingToggle();  }
  void Atom::FrameDrawStatisticsLogToggle() { logger->FrameDrawStatisticsLogToggle(); }
  void Atom::PositionLoggingToggle() {logger->PositionLoggingToggle();  }
  void Atom::PercentEnergyDissipatedToggle() { logger->PercentEnergyDissipatedToggle(); }
  void Atom::VelocityComponentsLogToggle() {logger->VelocityComponentsLogToggle();}
  void Atom::VelocityLoggingToggle() {logger->VelocityLoggingToggle();}
  void Atom::FastModeToggle     () {short_dt = (kShortDtSlow == short_dt) ? kLongDtFast : kShortDtSlow;  }
  void Atom::SlowMode           () {short_dt = kShortDtSlow;  }
  void Atom::TimeLoggingToggle  () {logger->TimeLoggingToggle();  }
  void Atom::WallClockLogToggle () {logger->WallClockLoggingToggle();  }

} // namespace
