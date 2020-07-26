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

#define GRAV_CONST (6.674 * pow(10,-11))
