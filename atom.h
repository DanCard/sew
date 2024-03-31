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


namespace sew {

class Atom {
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

  // Delta time in seconds.
  double dt = kShortDt;
  double time_ = 0;
  int    count = 0;         // Invocation count of moving particles.
  double total_potential_energy = 0;
  double total_kinetic_energy = 0;
  double total_energy = 0;    // Total energy of all particles = potential energy + kinetic energy
  // Energy of all particles when we had an escape.  Above this energy we start dissipating energy.
  double total_energy_escape = 0;
  double potential_energy_average;

  double pot_energy_cycle[kPFrequencySubDivisions];
  int    pot_energy_cycle_index = 0;
  double sum_p_energy = 0;        // Only used in below method.
  Logger* logger;

  Atom(int numParticles);

  void moveParticles();

private:
  void CalcAveragePotentialEnergy();
  void CalcEnergyAndLog();
  void AllForcesOnParticle(Particle * part_ptr);
};

} // namespace

#endif  // ATOM_H
