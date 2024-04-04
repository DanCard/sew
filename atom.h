#ifndef ATOM_H
#define ATOM_H
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

#include "constants.h"
#include "logger.h"
#include "particle.h"
#include "thread_pool.h"

namespace sew {

class Atom {
public:
  Particle * pars[kMaxParticles];  // par[ticle]s
  int num_particles = 0;

  volatile bool frame_draw_event_occurred;
  // Only used for logging.
  // If we are doing lots of calcs because electron close to proton, then we will wait for
  // draw event before finishing calcs.  That way we get maximum movement per frame.
  // In other words: simulation is too damn slow, and we want speed up so use all of the
  // compute budget per frame.
  // Number of times per screen logging that MoveParticles() completed before
  // next frame draw event.
  volatile int n_times_per_screen_log_MoveParticles_completed_before_next_frame_draw_event = 0;
  volatile int n_times_per_screen_log_MoveParticles_not_compl_before_next_frame_draw_event = 0;

  // Delta time in seconds.
  SFloat dt = kShortDt;
  SFloat time_ = 0;
  int    count = 0;         // Invocation count of moving particles.
  SFloat total_potential_energy = 0;
  SFloat total_kinetic_energy = 0;
  SFloat total_energy = 0;    // Total energy of all particles = potential energy + kinetic energy
  // Energy of all particles when we had an escape.  Above this energy we start dissipating energy.
  SFloat total_energy_cap = 1;
  SFloat potential_energy_average;

  SFloat pot_energy_cycle[kEFrequencySubDivisions];
  int    pot_energy_cycle_index = 0;
  SFloat sum_p_energy = 0;        // Only used in below method.
  Logger* logger;
  int iter = 0;   // Number of iterations to get significant movement.
  SFloat long_dt;
  SFloat short_dt;

  explicit Atom(int numParticles);

  void MoveParticles();

  void DtLoggingToggle();
  void EnergyLoggingToggle();
  void FastLoggingToggle();
  void FrameDrawStatisticsLoggingToggle();
  void IterationsLoggingToggle();
  void PositionLoggingToggle();
  void VelocityLoggingToggle();
  void PercentEnergyDissipatedLoggingToggle();
  void FastModeToggle();

private:
  void CalcAveragePotentialEnergy();
  void CalcEnergyOfAtom();
  void AllForcesOnParticle(Particle * part_ptr);
  void CalcForcesAndApply(Particle *part_ptr);
  ThreadPool* thread_pool;
};

} // namespace

#endif  // ATOM_H
