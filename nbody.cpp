#ifdef PYTHON
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
namespace py = pybind11;
#endif

#include "nbody.hpp"
#include <iostream>

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

struct Node Sim::init_head()
{
    struct Node head = {{0,0,0}, {0,0,0}, {0,0,0}, 0, nullptr};
    head.children = (struct Node*)malloc(8*sizeof(struct Node));
    for (Body& i : bodies) {
        for (int j = 0; j < 3; j++) {
            if (i.position[j] < head.min[j]) head.min[j] = i.position[j];
            if (i.position[j] > head.max[j]) head.max[j] = i.position[j];
        }
    }
    return head;
}

// TODO: Fix variable names -m It's too confusing.
int Sim::check_bodies(Eigen::Vector3d corner1, Eigen::Vector3d corner2)
{
    int num_bodies = 0;
    for (Body i : bodies) {
        bool greater = true;
        bool less = true;
        for (int j = 0; j < 3; j++) {
            if (i.position[j] < corner1[j]) less = false;
            if (i.position[j] > corner2[j]) greater = false;
        }
        if (less && greater) num_bodies++;
    }
    return num_bodies;
}

void Sim::make_child(struct Node& head, int iter, Eigen::Vector3d min, Eigen::Vector3d max)
{
    struct Node child = {{0,0,0}, {0,0,0}, {0,0,0}, 0, nullptr};
    child.children = (struct Node*)malloc(8*sizeof(struct Node));
    for (Body i : bodies) {
        bool inside = true;
        for (int j = 0; j < 3; j++) {
            if (i.position[j] < min[j]) inside = false;
            if (i.position[j] > max[j]) inside = false;
        }
        if (inside) {
            child.mass += i.mass;
            child.center += i.mass * i.position;
        }
    }
    // TODO: Investigate whether this has any real consequences if there are no bodies
    child.center /= child.mass;
    head.children[iter] = child;
}

void Sim::generate_tree(struct Node& head)
{
    if (check_bodies(head.min, head.max) < 2) return;
    Eigen::Vector3d perturb = (head.max - head.min)/2;
    for (int i = 0; i < 2; i++) {
        Eigen::Vector3d start = {0,0,0};
        start[2] = head.min[2] + perturb[2] * i;
        for (int j = 0; j < 2; j++) {
            start[1] = head.min[1] + perturb[1] * j;
            for (int k = 0; k < 2; k++) {
                start[0] = head.min[0] + perturb[0] * k;
                int iter = (i * 4) + (j * 2) + k;
                make_child(head, iter, start, start+perturb);
                generate_tree(head.children[iter]);
            }            
        }
    }
}

void Sim::tree_calc(Body& body)
{
    struct Node head = init_head();
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
        .def_readonly("bodies", &Sim::bodies);
}
#endif
