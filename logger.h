#ifndef LOGGING_H
#define LOGGING_H
// #include <algorithm>
#include <cmath>   // For M_PI constant and other match functions such as round.
#include <chrono>  // For logging every 500 ms
#include <cstdlib> // For rand()
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>

#include "atom.h"
#include "constants.h"
#include "particle.h"

namespace sew {

// Forward declarations
class Atom;
class Particle;

// Forward declarations
class Atom;

class Logger {
public:
  sew::Atom* a_;

  explicit Logger(Atom* a);

  bool energy_logging = false;
  bool velocity_logging = true;
  void EnergyLoggingToggle();
  void VelocityLoggingToggle();


  static void SetColorForConsole(unsigned char r, unsigned char g, unsigned char b) {
    printf("\x1b[38;2;%d;%d;%dm", r, g, b);
  }

  static bool IsSignificantParameter(const SFloat *d3, int i) {
    if (std::abs(d3[i]) < 1e-20) return false;
    if (i == 0) return d3[0] != 0;
    if (d3[0] == 0) return d3[i] != 0;
    // If within 1E4 of zero index.
    return std::abs(d3[0] / d3[i]) < 10;
  }

  static std::string Log3dArray(SFloat *d3, const std::string &name, int width = 6) {
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

  // A bit of a mess because we have particle data and particles(system) data that we are logging.
  std::string FormatLogLine(Particle* w, bool to_file) const;

  int w_to_log_id = 0;  // Rotate through particles to log to screen, when we don't

  // Log a particle and misc info.
  void LogStuff(Particle* w);
};

} // namespace

#endif  // LOGGING_H
