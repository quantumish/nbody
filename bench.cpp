#include <cmath>
#include "nbody.hpp"

int main()
{
    Sim sim (200, TreePM, Leapfrog);
    sim.add_body(pow(10, 9), {0, pow(10, 5), 0}, {1, 0, 0}, {0, 0, 0});
    sim.add_body(pow(10, 15), {0, 0, 0}, {0, 0, 0}, {0, 0, 0});
    sim.update();
}
