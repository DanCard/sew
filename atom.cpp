#include "atom.h"

#include <cassert>
#include <cstdlib> // For rand()
#include <iostream>
#include <thread>
#include <vector>

#include "logger.h"
#include "particle.h"


namespace sew {

Atom::Atom(int numParticles) : num_particles(numParticles) {
    if (num_particles == 0) return;
    logger = new sew::Logger(this);
    const int num_threads = std::thread::hardware_concurrency();
    std::cout << "\t\t\t num threads " << num_threads << std::endl;
    std::cout << "\t\t\t max speed electron " << kMaxSpeedElectron << "  kPFrequencySubDivisions " << kPFrequencySubDivisions
              << "\t kShortDt " << kShortDt << "  kLongDt " << kLongDt << std::endl;
    std::cout << "\t\t\t\t kEFrequency " << kEFrequency << "  kPFrequency " << kPFrequency << std::endl;
    // std::cout << "\t\t kBohrMagneton " << kBohrMagneton << "  kProtonMagneticMoment " << kProtonMagneticMoment << std::endl;

    assert(numParticles <= kMaxParticles);
    int divider;  // Prefer bright colors, but with many particles becomes indistinguishable.
         if (num_particles <= 2)  divider = 1;
    else if (num_particles <= 4)  divider = 3;
    else if (num_particles <= 6)  divider = 4;
    else                          divider = 8;  // With more particles don't brighten as much.
    for (int i = 0; i < numParticles; ++i) {
      Particle* p;
      if (i < numParticles / 2) {
        pars[i] = new Electron(i, this, logger);
        p = pars[i];
        // Set pseudo random colors.  Electrons tend to be more red.  Protons tend to be more blue.
        // Prefer bright colors over dark colors.
        // p->color[0] = 151 + (std::rand() % 105);
        // p->color[1] =  51 + (std::rand() % 155);
        // p->color[2] =   0 + (std::rand() % 240);
        p->color[0] = std::rand() % 256;
        p->color[1] = std::rand() % 256;
        p->color[2] = std::rand() % 245;
        if (i == 0) {
          p->vel[1] = -1e4;
          p->pos[0] = - kBohrRadius;
          p->color[0] = 255;
          p->color[1] = 120;
          p->color[2] = 120;
        }
      } else {
        pars[i] = new Proton(i, this, logger);
        p = pars[i];
        p->color[0] =   0 + (std::rand() % 210);
        p->color[1] =   0 + (std::rand() % 240);
        p->color[2] = 151 + (std::rand() % 105);
        p->vel[1] = 1e4;
      }
      std::cout << "\t\t color " << int(p->color[0]) << " " << int(p->color[1])
                << " " << int(p->color[2]);
      // double max_dist = p->max_dist_allow * 0.5;
      for (int j = 0; j < 3; ++j) {
        if (num_particles > 2) {
          // Set random locations
          p->pos[j] = (std::rand() / (RAND_MAX + 1.0) - 0.5) * kBohrRadiusProton;
        }
        // Increase brightness
        int increase = p->color[j] / divider;
        if (p->color[j] + increase > 255) p->color[j]  = 255;
        else                              p->color[j] += increase;
      }
      std::cout << "\t\t color " << int(p->color[0]) << " " << int(p->color[1])
                << " "           << int(p->color[2]) << std::endl;
      std::cout << "\t\t pos " << p->pos[0] << " " << p->pos[1] << " " << p->pos[2] << std::endl;
      std::cout << "\t\t this " << this << "  p->a_ " << p->a_ << std::endl;
      assert(p->a_ == this);
    }
    for (double & p_energy_cycle_ : pot_energy_cycle) {
      p_energy_cycle_ = 0;
    }
    assert(pars[0]->a_ == this);
    assert(pars[1]->a_ == this);
  }

  // Potential energy changes because of sinusoidal charge frequency.
void Atom::CalcAveragePotentialEnergy() {
    sum_p_energy += pot_energy_cycle[pot_energy_cycle_index++];
    pot_energy_cycle_index = pot_energy_cycle_index % kPFrequencySubDivisions;
    // Subtract the potential energy that is going to roll off our average.
    sum_p_energy -= pot_energy_cycle[pot_energy_cycle_index];
    // -1 because of the value subtracted above.
    potential_energy_average = sum_p_energy / (kPFrequencySubDivisions - 1);
  }

void Atom::CalcEnergyAndLog() {
    // Calculate potential energy.
    total_potential_energy = 0;
    total_kinetic_energy = 0;
    double closest = 1;  // meters.  Just setting to a large number.

    // Potential energy = sum of all potential energies.
    // Set the potential energy for each pair of particles.
    for (int i = 0; i < num_particles; ++i) {
      Particle *w1 = pars[i];   // Wave 1
      total_kinetic_energy += 0.5 * w1->mass_kg * w1->vel_mag2;
      for (int j = i + 1; j < num_particles; ++j) {
        Particle *w2 = pars[j];   // Wave 2
        assert(w1->dist_mag[j] != 0);
        total_potential_energy += kCoulomb * w1->freq_charge * w2->freq_charge / w1->dist_mag[j];
      }
      // If the electron is interesting because it is close to a proton than give it preferential logging.
      if (w1->is_electron && w1->dist_mag_closest < kCloseToTrouble
          && w1->dist_mag_closest < closest) {
        closest = w1->dist_mag_closest;
        logger->w_to_log_id = w1->id;
      }
    }
    assert(num_particles > 0);

    pot_energy_cycle[pot_energy_cycle_index] = total_potential_energy;
    pot_energy_cycle_index = (pot_energy_cycle_index + 1) % kPFrequencySubDivisions;
    // Potential energy changes because of sinusoidal charge frequency.
    // Alternatively could just use average charge.
    CalcAveragePotentialEnergy();
    total_energy = potential_energy_average + total_kinetic_energy;

    for (int i = 0; i < num_particles; ++i) {
      logger->LogStuff(pars[i]);
    }
  }

void Atom::AllForcesOnParticle(Particle * part_ptr) {
    assert(num_particles > 0);
    assert(part_ptr->a_ == this);

    int part_num = part_ptr->id;
    part_ptr->InitVarsToCalcForces();
    for (int i = 0; i < num_particles; ++i) {
      if (i == part_num) continue;
      part_ptr->CalcForcesFromParticle(pars[i]);
    }
  }

  // Called once for every screen draw.
void Atom::moveParticles() {
    assert(pars[0]->a_ == this);

    double pos_change_per_particle[kMaxParticles];
    for (int j = 0; j < num_particles; ++j) {
      pos_change_per_particle[j] = 0;
    }
    for (int i=0; i<num_particles; ++i) {
      // Only check for escape once per significant movement to save on CPU.
      pars[i]->CheckForEscape();
    }
    // Do until we get significant movement, then wait for screen draw.
    for (int iter = 0; iter<(4096 * 2) && !screen_draw_event_occurred; ++iter) {
      for (int i = 0; i < num_particles; ++i) {
        Particle * wave_ptr = pars[i];
        wave_ptr->freq_charge = wave_ptr->ChargeSinusoidal();
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

} // namespace
