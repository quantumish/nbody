#ifdef PYTHON
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
namespace py = pybind11;
#endif

#include "nbody.hpp"
#include <iostream>

#define DIST_THRESHOLD 10000

Body::Body(double m, Eigen::Vector3d x, Eigen::Vector3d v, Eigen::Vector3d a)
    :mass(m), position(x), velocity(v), acceleration(a)
{
    net_force = Eigen::Vector3d::Zero();
    assert(mass >= 0);
}

bool Body::operator==(Body& other)
{
    if (this == &other) return true;
    else return false;
}
  
Sim::Sim(double delta_t, ForceMethod fm, TimeMethod tm)
    :dt(delta_t), force_method(fm), time_method(tm), t(0)
{
    if (fm == Tree) {
        head = init_head();
        generate_tree(head);
    }
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
    Eigen::Vector3d sum = Eigen::Vector3d::Zero();
    for (Body& other : bodies) {
        if (body == other) continue;
        double distance = sqrt(pow(body.position[0] - other.position[0],2)+pow(body.position[1] - other.position[1],2))+pow(body.position[2] - other.position[2],2);
        double magnitude = (GRAV_CONST * body.mass * other.mass)/(pow(distance,2));
        sum += (other.position - body.position).normalized() * magnitude;
    }
    body.net_force = sum;
}

void Sim::initialize_children(struct Node& node)
{
    Eigen::Vector3d half ((node.max-node.min)/2);
    for (int i = 0; i < 8; i++) {
        Eigen::Vector3d child_min (min+(half[0]*((i & 1) == 1)), min+(half[1]*((i & 2) == 2)),
                                   min+(half[2]*((i & 4) == 4)));
        octree.emplace_back(child_min, child_min + half, {0, 0, 0}, 0, nullptr, nullptr);
        node.children[i] = &octree[octree.size()-1];
    }
}

void Sim::initialize_octree()
{
    octree.clear();
    for (Body& body : bodies) {
        for (int i = 0; i < 3; i++) {
            if (body.position[i] > max[i]) max[i] = body.position[i];
            if (body.position[i] < min[i]) min[i] = body.position[i];
        }       
    }
    octree.emplace_back(min, max, bodies[0].position, bodies[0].mass, &bodies[0], nullptr);
}

bool vector_within(Eigen::Vector3d v, Eigen::Vector3d a, Eigen::Vector3d b)
{
    return v[0] > a[0] && v[1] > a[1] && v[2] > a[2] && v[0] < b[0] && v[1] < b[1] && v[2] < b[2];
}

void Sim::insert_body(Body& body)
{
    std::stack<Node> nodes;
    stack.push(octree[0]);
    while (stack.size()) {
        if (!vector_within(body.position, stack.top.min, stack.top.max) {
            stack.pop();
            continue;
        }            
        if (stack.top.body == nullptr && stack.top.children == nullptr) {
            stack.top.body = body;
            break;
        }
        if (stack.top.body != nullptr && stack.top.children == nullptr) {
            // init children, move stuff
        }
        for (int i = 0; i < 8; i++) {
            stack.push(stack.top.children[i]);
        }
        stack.pop();
    }
}

void Sim::calc_net_force(Body& body)
{
    switch(fm) {
    case Direct:
        direct_calc();
        break;
    case TreePM:
        break;
    }
}

// TODO: Look into Leapfrog derivation | -m The derivation of this scares me. It's probably a good idea to look at it more though...
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
