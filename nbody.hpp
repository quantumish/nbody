#include <vector>
#include <array>
#include <Eigen/Dense>

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

struct Node {
    Eigen::Vector3d min;
    Eigen::Vector3d max;
    Eigen::Vector3d center;
    double mass;
    Body* body;
    Node* children[8];    
    Node(Eigen::Vector3d mini, Eigen::Vector3d maxi, Eigen::Vector3d c,
           double m, Body* b, Node* child);
};

enum ForceMethod {Direct, Tree, FMM, Mesh, P3M, TreePM};
enum TimeMethod {Euler, Leapfrog, Hermite};

class Sim {
    Eigen::Vector3d min {0,0,0};
    Eigen::Vector3d max {0,0,0};
    Node* octree;
    ForceMethod force_method;
    TimeMethod time_method;

    void initialize_octree();
    void dump_tree();
    void purge_tree();
    void insert_body(Body& body);
    void initialize_children(Node& node);
    void calc_center_mass(Node& node);
    void direct_calc(Body& body);
    void tree_calc(Body& body);
    void calc_net_force(Body& body);
    void leapfrog_update(Body& body);
public:
    std::vector<Body> bodies;
    double dt;
    double t;
  
    Sim(double delta_t, ForceMethod fm, TimeMethod tm);
    void add_body(double m, Eigen::Vector3d x, Eigen::Vector3d v, Eigen::Vector3d a);
    void add_body(Body body);
    void set_bodies(std::vector<Body> new_bodies);
    void update();
};

#define GRAV_CONST (6.674 * pow(10,-11))
