#include <vector>
#include "/Users/alexmann/Library/Mobile Documents/com~apple~CloudDocs/Documents/School/Carleton/Masters/program_working_dir/inviscid_model/spline.h"

#ifndef FILAMENT_H
#define FILAMENT_H

struct Node {
  int id;
  int born_on;
  double x;
  double y;
  double z;
  double vx;
  double vy;
  double vz;
  double ivx;
  double ivy;
  double ivz;
  double x_prev;
  double y_prev;
  double z_prev;
  double vcx;
  double vcy;
  double vcz;
  double sx;
  double sy;
  double sz;
  double sl;
  double vol;
  double r;
  Node* prev;
  Node* next;
};

Node* new_node(double x, double y, double z, double r, Node* prev, Node* next);

//sinebackground class:
class SineBackgroundGenerator {
public:
  tk::spline bottom_ufunc;
  //default(empty) constructor
  SineBackgroundGenerator() {}
  
  //constructor(implementation in filament.cpp)
  SineBackgroundGenerator(std::string background_velocity_file,
                          std::vector<double> amplitudes,
                          std::vector<double> wavelengths,
                          std::vector<double> phases,
                          std::vector<int> directions);
  
  //(implementation in filament.cpp)
  std::array<double, 3> evaluate(double x, double y, double z, int t);

  private:
    std::vector<double> amplitudes;
    std::vector<double> wavelengths;
    std::vector<double> phases;
    std::vector<int> directions;
    int n;
    std::pair<std::vector<double>, std::vector<double>> data;
};


class Filament {
  int node_count;
  double circulation;
  double core_radius;
  double max_sl;
  double min_sl;
  double max_range;
  int timestep;
  double t;
  int nodes_added;
  int nodes_deleted;

public:
  //public data member (for evaluation)
  SineBackgroundGenerator background_velocity;

  Node* head;
  Node* tail;
  //default(empty) constructor
  Filament(){}

  //constructor (implementation in filament.cpp)
  Filament(double circulation, double core_radius, double max_sl, double min_sl, double max_range,
           std::vector<std::vector<double>> coords, SineBackgroundGenerator background_velocity, bool circular);

  ~Filament(); //declare destructor (implementation in filament.cpp)
  //declare other member functions (implementation in filament.cpp)
  void step(double dt);
  void compute_segments_and_repartition_if_needed();
  void split_segment(Node* n1, Node* n2);
  void merge_segment(Node* n1, Node* n2);
  
  bool influenced_by_segment(double x, double y, double z, Node* n1);
  double perpendicular_distance_to_segment(double x, double y, double z, Node* n1);
  void compute_single_segment_geometric_property(Node* n1);
  std::string save_single_timestep(std::string fn_template);
  Node* getHead() const;
  Node* getTail() const;
  int get_node_count() const;
  double get_circulation() const;
  double get_time() const;
  int get_timestep() const;
  int get_nodes_added() const;
  double get_max_range() const;
  void setTimestep(int newTimestep);
  void setT(double newT);
  void og_step(double dt);
  void og_move_nodes_via_induction(double dt);
};

//Declare standalone functions used within Filament class (implementation in filament.cpp)
double norm(double x, double y, double z);
std::array<double, 3> cross(double x1, double y1, double z1, double x2, double y2, double z2);
std::pair<std::vector<double>, std::vector<double> > load_from_file(const std::string &filename);
std::vector<std::vector<double>> get_fil_coords(bool fil_type, double x, double y, double domain_width, double dz_i);
std::vector<std::vector<double>> get_fil_coords2(bool fil_type, double x1, double y1, double z1, double x2, double y2, double z2, double dz_i);
void calculate_induced_velocity(Node& control_node, std::vector<std::unique_ptr<Filament>>& filaments, double& ivx, double& ivy, double& ivz, int node_num);
void move_nodes_via_induction(std::vector<std::unique_ptr<Filament>>& filaments, double dt);
void step(std::vector<std::unique_ptr<Filament>>& filaments, double dt);
std::string save_multiple_fil_timestep(std::string fn_template, const std::vector<std::unique_ptr<Filament>>& filaments);
void read_vtk(std::string fn, std::vector<std::unique_ptr<Filament>>& filaments, double circulation, double core_radius, double max_sl, double min_sl, double max_range, SineBackgroundGenerator background_vel);

#endif // FILAMENT_H