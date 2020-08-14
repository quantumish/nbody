#include <vector>
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
    Node* children;
};

enum ForceMethod {Direct, Tree, FMM, Mesh, P3M};
enum TimeMethod {Euler, Leapfrog, Hermite};

class Sim {
    struct Node init_head();
    int check_bodies(Eigen::Vector3d corner1, Eigen::Vector3d corner2);
    void make_child(struct Node& head, int iter, Eigen::Vector3d min, Eigen::Vector3d max);
    void generate_tree(struct Node& head);

    struct Node head;
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
    double t;
  
    Sim(double delta_t, ForceMethod fm, TimeMethod tm);
    void add_body(double m, Eigen::Vector3d x, Eigen::Vector3d v, Eigen::Vector3d a);
    void add_body(Body body);
    void set_bodies(std::vector<Body> new_bodies);
    void update();
};

#define GRAV_CONST (6.674 * pow(10,-11))
