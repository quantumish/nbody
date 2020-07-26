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
  Body(double m, Eigen::Vector3d x, Eigen::Vector3d v, Eigen::Vector3d a);
  bool operator==(Body& other);
};

Body::Body(double m, Eigen::Vector3d x, Eigen::Vector3d v, Eigen::Vector3d a)
  :mass(m), position(x), velocity(v), acceleration(a)
{
  net_force = {0,0,0};
  assert(mass >= 0);
}

bool Body::operator==(Body& other)
{
  if (this == &other) return true;
  else return false;
}

enum ForceMethod {Direct, Tree, FMM, Mesh, P3M};
enum TimeMethod {Euler, Leapfrog, Hermite};

class Sim {
  ForceMethod force_method;
  TimeMethod time_method;
  void direct_calc(Body& body);
  void tree_calc(Body& body);
  void fmm_calc(Body& body);
  void mesh_calc(Body& body);
  void p3m_calc(Body& body);
  void calc_net_force(Body& body);
  void leapfrog_update(Body& body);
public:
  std::vector<Body> bodies;
  double dt;
  double t = 0;
  
  Sim(double delta_t, ForceMethod fm, TimeMethod tm);
  void add_body(double m, Eigen::Vector3d x, Eigen::Vector3d v, Eigen::Vector3d a);
  void add_body(Body body);
  void set_bodies(std::vector<Body> new_bodies);
  void update();
};
  
Sim::Sim(double delta_t, ForceMethod fm, TimeMethod tm)
  :dt(delta_t), force_method(fm), time_method(tm)
{
  __asm__("nop");
}

void Sim::add_body(double m, Eigen::Vector3d x, Eigen::Vector3d v, Eigen::Vector3d a)
{
  bodies.emplace_back(m,x,v,a);
}

void Sim::add_body(Body body)
{
  bodies.push_back(body);
}

void Sim::set_bodies(std::vector<Body> new_bodies)
{
  bodies = new_bodies;
}

void Sim::direct_calc(Body& body)
{
  Eigen::Vector3d sum = {0,0,0};
  for (Body& other : bodies) {
    if (body == other) continue;
    double distance = sqrt(pow(body.position[0] - other.position[0],2)+pow(body.position[1] - other.position[1],2))+pow(body.position[2] - other.position[2],2);
    double magnitude = (GRAV_CONST * body.mass * other.mass)/(pow(distance,2));
    sum += (other.position - body.position).normalized() * magnitude;
  }
  body.net_force = sum;
}

void Sim::tree_calc(Body& body)
{
}

void Sim::fmm_calc(Body& body)
{
}

void Sim::mesh_calc(Body& body)
{
}

void Sim::p3m_calc(Body& body)
{
}

void Sim::calc_net_force(Body& body)
{
  switch (force_method) {
  case Direct:
    direct_calc(body);
  case Tree:
    tree_calc(body);
  case FMM:
    fmm_calc(body);
  case Mesh:
    mesh_calc(body);
  case P3M:
    p3m_calc(body);
  }
}

void Sim::leapfrog_update(Body& body)
{
  body.position += (body.velocity * dt) + (0.5 * body.acceleration * pow(dt,2));
  Eigen::Vector3d a_0 = body.acceleration;
  calc_net_force(body);
  body.acceleration = (body.net_force / body.mass);
  body.velocity += 0.5 * (a_0 + body.acceleration) * dt;
}

void Sim::update()
{
  for (Body& body : bodies) {
    switch (time_method) {
    case Euler:
      body.position += body.velocity * dt;
      body.velocity += body.acceleration * dt;
      calc_net_force(body);
      body.acceleration = (body.net_force / body.mass);
      break;
    case Leapfrog:
      // TODO: The derivation of this scares me. It's probably a good idea to look at it more though...
      leapfrog_update(body);
    case Hermite:
      __asm__("nop");
    }
  }
  t+=dt;
}

#ifdef PYTHON
PYBIND11_MODULE(nbody, m) {
  m.doc() = "Simulate exoplanets with C++.";
  py::enum_<ForceMethod>(m, "ForceMethod")
    .value("Direct", ForceMethod::Direct)
    .value("Tree", ForceMethod::Tree)
    .value("FMM", ForceMethod::FMM)
    .value("Mesh", ForceMethod::Mesh)
    .value("P3M", ForceMethod::P3M)
    .export_values();
  py::enum_<TimeMethod>(m, "TimeMethod")
    .value("Euler", TimeMethod::Euler)
    .value("Leapfrog", TimeMethod::Leapfrog)
    .value("Hermite", TimeMethod::Hermite)
    .export_values();
  py::class_<Body>(m, "Body")
    .def(py::init<double, Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d>())
    .def_readonly("mass", &Body::mass)
    .def_readonly("position", &Body::position)
    .def_readonly("velocity", &Body::velocity)
    .def_readonly("acceleration", &Body::acceleration)
    .def_readonly("net_force", &Body::net_force);
  py::class_<Sim>(m, "Sim")
    .def(py::init<double, ForceMethod, TimeMethod>())
    .def("add_body", (void (Sim::*)(double,Eigen::Vector3d,Eigen::Vector3d,Eigen::Vector3d)) &Sim::add_body, py::arg("m"), py::arg("x"), py::arg("v"), py::arg("a"))
    .def("add_body", (void (Sim::*)(Body)) &Sim::add_body, py::arg("body"))
    .def("set_bodies", &Sim::set_bodies, py::arg("new_bodies"))
    .def("update", &Sim::update)
    .def_readonly("t", &Sim::t)
    .def_readonly("bodies", &Sim::bodies);
}
#endif
