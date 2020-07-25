#include <vector>
#include <Eigen/Dense>

class Body {
public:
  double mass;
  Eigen::Vector3d position;
  Eigen::Vector3d velocity;
  Eigen::Vector3d acceleration;
  Eigen::Vector3d net_force;
  Body(double m, Eigen::Vector3d, x, Eigen::Vector3d v, Eigen::Vector3d a);
};

Body::Body(double m, Eigen::Vector3d x, Eigen::Vector3d v, Eigen::Vector3d a)
  :mass(m), position(x), velocity(v), acceleration(a)
{
  net_force = {0,0,0};
  assert(mass >= 0);
}


class Simulation {
public:
  std::vector<Body> bodies;
  Simulation();
  void add_body(double m, Eigen::Vector3d, x, Eigen::Vector3d v, Eigen::Vector3d a);
};
  
Simulation::Simulation()
{
  __asm__("nop");
}

Simulation::add_body(double m, Eigen::Vector3d, x, Eigen::Vector3d v, Eigen::Vector3d a)
{
  bodies.emplace_back(m,x,v,a);
}
