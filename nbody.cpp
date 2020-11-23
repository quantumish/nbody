#ifdef PYTHON
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
namespace py = pybind11;
#endif

#include "nbody.hpp"
#include <iostream>
#include <stack>

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
        double distance = sqrt(pow(body.position[0] - other.position[0],2)+pow(body.position[1] - other.position[1],2)+pow(body.position[2] - other.position[2],2));
        double magnitude = (GRAV_CONST * body.mass * other.mass)/(pow(distance,2));
        sum += (other.position - body.position).normalized() * magnitude;
    }
    body.net_force = sum;
}

bool vector_within(Eigen::Vector3d v, Eigen::Vector3d a, Eigen::Vector3d b)
{
    std::cout << v << "(vec)\n\n";
    std::cout << a << "\n\n";
    std::cout << b << "\n\n";
    //std::cout << (v[0] > a[0] && v[1] > a[1] && v[2] > a[2] && v[0] < b[0] && v[1] < b[1] && v[2] < b[2]) << "\n";
    return v[0] >= a[0] && v[1] >= a[1] && v[2] >= a[2] && v[0] < b[0] && v[1] < b[1] && v[2] < b[2];
}

void Sim::initialize_children(struct Node& node)
{
    Eigen::Vector3d half = (node.max-node.min)/2;
    for (int i = 0; i < 8; i++) {
        Eigen::Vector3d child_min(node.min[0] + (half[0] * ((i & 1) == 1)),
                                  node.min[1] + (half[1] * ((i & 2) == 2)),
                                  node.min[2] + (half[2] * ((i & 4) == 4)));
        octree.push_back({child_min, child_min + half, {0, 0, 0}, 0, nullptr, nullptr});
        node.children[i] = &octree[octree.size()-1];
    }
    std::cout << "done" << "\n";
    for (int i = 0; i < 8; i++) {
        if (vector_within(node.body->position, node.children[i]->min, node.children[i]->max)) {
            std::cout << "hi2" << "\n";
            //node.children[i]->body = node.body;
            std::cout << i << "\n";
            std::cout << node.children[i] << "\n";
            node.children[i]->body = 0x0;
            std::cout << "hi3" << "\n";
            node.body = nullptr;
            break;
        }
    }
    std::cout << node.children << "\n\n";
    for (int i = 0; i < 8; i++) std::cout << node.children[i] << "\n";
    std::cout << "??!" << "\n";
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
    for (int i = 0; i < 3; i++) {
        max[i] += 1;
        min[i] -= 1;
    }
    octree.push_back({min, max, {0, 0, 0}, 0, nullptr, nullptr});
}

void Sim::insert_body(Body& body)
{
    std::cout << body.position << "pos\n";
    std::stack<Node*> stack;
    stack.push(&octree[0]);
    while (stack.size()) {
        if (!vector_within(body.position, stack.top()->min, stack.top()->max)) {
            stack.pop();
            continue;
        }
        std::cout << "Made it!" << "\n";
        std::cout << stack.top()->body << " " << stack.top()->children[0] << "\n";
        if (stack.top()->body == nullptr && stack.top()->children[0] == nullptr) {
            std::cout << "hmm" << "\n";
            stack.top()->body = &body;
            break;
        }
        std::cout << "uhoh" << "\n";
        if (stack.top()->body != nullptr && stack.top()->children[0] == nullptr) {
            std::cout << stack.top()->children << "\n\n";
            for (int i = 0; i < 8; i++) std::cout << stack.top()->children[i] << "\n";
            initialize_children(*stack.top());
        }
        //std::cout << stack.top()->children[2] << "\n\n";
        std::cout << stack.top()->children << "\n\n";
        for (int i = 0; i < 8; i++) {
            std::cout << stack.top()->children[i] << "\n";
            //stack.push(stack.top()->children[i]);
        }
        std::cout << "23" << "\n";
        stack.pop();
    }
    std::cout << "??" << "\n";
}

void Sim::calc_center_mass(Node& node)
{
    if (node.children[0] == nullptr) {
        if (node.body == nullptr) return;
        node.center = node.body->position;
        node.mass = node.body->mass;
        return;
    }
    for (int i = 0; i < 8; i++) calc_center_mass(*node.children[i]);
    double mass_sum = 0;
    Eigen::Vector3d weighted_sum = {0, 0, 0};
    for (int i = 0; i < 8; i++) {
        weighted_sum += node.children[i]->mass * node.children[i]->center;
        mass_sum += node.children[i]->mass;
    }
    node.center = weighted_sum/mass_sum;
    node.mass = mass_sum;
}

#define THETA 1
void Sim::tree_calc(Body& body)
{
    std::stack<Node*> stack;
    Eigen::Vector3d sum;
    stack.push(&octree[0]);
    double distance = sqrt(pow(body.position[0] - stack.top()->center[0],2) +
                           pow(body.position[1] - stack.top()->center[1],2) +
                           pow(body.position[2] - stack.top()->center[2],2));
    while (stack.size()) {
        double maxdim = 0;
        for (int i = 0; i < 3; i++) {
            if (stack.top()->max[i]-stack.top()->min[i] > maxdim) {
                maxdim = stack.top()->max[i]-stack.top()->min[i];
            }
        } 
        if (distance/maxdim >= THETA) {
            double magnitude = (GRAV_CONST * body.mass * stack.top()->mass)/(pow(distance,2));
            sum += (stack.top()->center - body.position).normalized() * magnitude;
            stack.pop();
            continue;
        }
        for (int i = 0; i < 8; i++) stack.push(stack.top()->children[i]);
        stack.pop();
    }
    body.net_force = sum;
}

void Sim::calc_net_force(Body& body)
{
    switch(force_method) {
    case Direct:
        direct_calc(body);
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
    if (force_method == TreePM) {
        initialize_octree();
        for (Body& body : bodies) insert_body(body);
        std::cout << "all done" << "\n";
        calc_center_mass(octree[0]);
        for (Node node : octree) {
            std::cout << "ihudwahkjadsjk" << node.min << "\n\n" << node.max << "\n\n" << node.body << "\n\n" << node.center << "\n\n" << node.mass << "\n\n\n\n";
        }
    }
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
        .value("TreePM", ForceMethod::TreePM)
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
