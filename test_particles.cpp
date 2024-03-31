#include <cassert>
#include "particle.h"

namespace sew {

static void testDummiesInverseExponential() {
  double result1 = Particle::dummiesInverseExponential(0.5);
  // std::cout << result1 << std::endl;
  assert(result1 == 0.5);

  // Test case 2: fast_fraction = 0.2
  double result2 = Particle::dummiesInverseExponential(0.2);
  // std::cout << result2 << std::endl;
  assert(result2 == 0.8);

  // Test case 3: fast_fraction = 0.09
  double result3 = Particle::dummiesInverseExponential(0.1);
  std::cout << result3 << std::endl;
  assert(result3 == .9);

  std::cout << "\t test 4" << std::endl;
  // Test case 4: fast_fraction = 0.04
  double result4 = Particle::dummiesInverseExponential(0.09);
  std::cout << "\t fast = 9%  inverse: " << result4 << std::endl;
  assert(result4 == .918);

  double result5 = Particle::dummiesInverseExponential(0.05);
  std::cout << "\t fast = 5%  inverse: " << result5 << std::endl;
  assert(result5 == .99);

  double result6 = Particle::dummiesInverseExponential(0.03);
  std::cout << "\t fast = 3%  slow = 97%  inverse: " << result6 << std::endl;
  assert(result6 == .999);

  double result7 = Particle::dummiesInverseExponential(0.01);
  std::cout << "\t fast = 1%  inverse: " << result7 << std::endl;
  assert(result7 == .9999);

  double result8 = Particle::dummiesInverseExponential(0.0011);
  std::cout << "\t fast = 0.11%  inverse: " << result8 << std::endl;
  assert(result8 == 0.999989);  // Is it correct?

  double result9 = Particle::dummiesInverseExponential(0.001);
  std::cout << "\t fast = 0.1%  inverse: " << result9 << std::endl;
  assert(result9 == 1.0);

  // Add more test cases as needed...
}

  
}

int main() {
  sew::testDummiesInverseExponential();
  return 0;
}
