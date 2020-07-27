#ifdef PYTHON
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
namespace py = pybind11;
#endif

#include "nbody.hpp"

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
  
Sim::Sim(double delta_t, ForceMethod fm, TimeMethod tm)
  :dt(delta_t), force_method(fm), time_method(tm), bound_max{0,0,0}, bound_min{0,0,0}
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
  get_box();
  Eigen::Vector3d sum = {0,0,0};
  for (Body& other : bodies) {
    if (body == other) continue;
    double distance = sqrt(pow(body.position[0] - other.position[0],2)+pow(body.position[1] - other.position[1],2))+pow(body.position[2] - other.position[2],2);
    double magnitude = (GRAV_CONST * body.mass * other.mass)/(pow(distance,2));
    sum += (other.position - body.position).normalized() * magnitude;
  }
  body.net_force = sum;
}

void Sim::get_box()
{
  for (Body& i : bodies) {
    for (int j = 0; j < 3; j++) {
      if (i.position[j] < bound_min[j]) bound_min[j] = i.position[j];
      if (i.position[j] > bound_max[j]) bound_max[j] = i.position[j];
    }
  }
}

void Sim::check_for_planet(struct Node node, Eigen::Vector3d corner1, Eigen::Vector3d corner2)
{
  for (Body i : bodies) {
    bool greater = true;
    bool less = true;
    for (int j = 0; j < 3; j++) {
      std::cout << i.position[j] << " vs. " << corner1[j] << " and " << corner2[j] << "\n";
      if (i.position[j] < corner1[j]) {
        std::cout << "Failed LESS!\n";
        less = false;
      }
      if (i.position[j] > corner2[j]) {
        std::cout << "Failed GREATER!\n";
        greater = false;
      }
    }
    if (less && greater) {
      std::cout << "Ladies and gentlemen... we got 'em\n";
      node.bodies.push_back(&i);
    }
  }
}

void Sim::generate_tree(struct Node head)
{
  struct Node current = head;
  Eigen::Vector3d box_max = bound_max;
  Eigen::Vector3d box_min = bound_min;
  for (int i = 0; i < 1; i++) {
    double distance = sqrt(pow(box_max[0] - box_min[0],2)+pow(box_max[1] - box_min[1],2))+pow(box_max[2] - box_min[2],2);
    std::vector<Body*> empty;
    Node temp = {empty, NULL};
    current.children[0] = &temp;
    check_for_planet(*current.children[0], box_min, (box_min.array()+(distance/2)).matrix());
    for (Body* j : current.children[0]->bodies) {
      std::cout << j->mass << " MASS? \n";
    }
  }
}

void Sim::tree_calc(Body& body)
{
  get_box();
  struct Node head = {};
  generate_tree(head);
  assert(1<0);
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
    break;
  case Tree:
    tree_calc(body);
    break;
  case FMM:
    fmm_calc(body);
    break;
  case Mesh:
    mesh_calc(body);
    break;
  case P3M:
    p3m_calc(body);
    break;
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
    .def_readonly("bound_max", &Sim::bound_max)
    .def_readonly("bound_min", &Sim::bound_min)
    .def_readonly("bodies", &Sim::bodies);
}
#endif
