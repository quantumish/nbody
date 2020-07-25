#include <vector>
#include <Eigen/Dense>

#ifdef PYTHON
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
namespace py = pybind11;
#endif

#define GRAV_CONST (6.674 * pow(10,-11))

class Body {
public:
  double mass;
  Eigen::Vector3d position;
  Eigen::Vector3d velocity;
  Eigen::Vector3d acceleration;
  Eigen::Vector3d net_force;
  Body(double m, Eigen::Vector3d, x, Eigen::Vector3d v, Eigen::Vector3d a);
  bool operator==(Body& other);
};

Body::Body(double m, Eigen::Vector3d x, Eigen::Vector3d v, Eigen::Vector3d a)
  :mass(m), position(x), velocity(v), acceleration(a)
{
  net_force = {0,0,0};
  assert(mass >= 0);
}

Body::operator==(Body& other)
{
  if (this == &other) return true;
  else return false;
}


class Sim {
public:
  std::vector<Body> bodies;
  double dt;
  double t = 0;
  
  Sim();
  void add_body(double m, Eigen::Vector3d, x, Eigen::Vector3d v, Eigen::Vector3d a);
  void calc_net_force(Body body);
  void update();
};
  
Sim::Sim(double delta_t)
  :dt(delta_t)
{
  __asm__("nop");
}

void Sim::add_body(double m, Eigen::Vector3d, x, Eigen::Vector3d v, Eigen::Vector3d a)
{
  bodies.emplace_back(m,x,v,a);
}

void Sim::calc_net_force(Body& body)
{
  for (Body other : bodies) {
    if (body == other) continue;
    double distance = sqrt(pow(body.position[0] - other.position[0],2)+pow(body.position[1] - other.position[1],2))+pow(body.position[2] - other.position[2],2);
    body.net_force += (GRAV_CONST * body.mass * other.mass)/(pow(distance,2))
  }
}

void Sim::update()
{
  for (Body body : bodies) {
    body.position += body.velocity * dt;
    body.velocity += (body.net_force / m) * dt;
    body.net_force = {0,0,0};
  }
  t+=dt;
}

#ifdef PYTHON
PYBIND11_MODULE(astrosim, m) {
  m.doc() = "Simulate exoplanets with C++.";
  
  py::class_<Matter>(m, "Matter")
    .def(py::init<double, double, Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d>())
    .def_readonly("mass", &Matter::mass)
    .def_readonly("position", &Matter::position)
    .def_readonly("velocity", &Matter::velocity)
    .def_readonly("acceleration", &Matter::acceleration);
  
  py::class_<Universe>(m, "Universe")
    .def(py::init<double>())
    .def("add_body", &Universe::add_matter, py::arg("m"), py::arg("x"), py::arg("v"), py::arg("a"))
    .def("update", &Universe::update)
    .def_readonly("t", &Universe::ticks)
    .def_readonly("bodies", &Universe::bodies);
}
#endif
