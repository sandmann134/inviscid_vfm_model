#include <cmath>
#include <vtkUnstructuredGrid.h>
#include <vtkFieldData.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkSmartPointer.h>
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkFieldData.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <array>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <map>
#include <tuple>
#include "/Users/alexmann/Library/Mobile Documents/com~apple~CloudDocs/Documents/School/Carleton/Masters/program_working_dir/inviscid_model/spline.h"
#include "/Users/alexmann/Library/Mobile Documents/com~apple~CloudDocs/Documents/School/Carleton/Masters/program_working_dir/inviscid_model/filament.h"

//vector magnitude 'norm' fxn:
double norm(double x, double y, double z){
    return sqrt(x*x + y*y + z*z);
}

//cross product function for 2 3D vector/arrays (defined by individual components here) 
std::array<double, 3> cross(double x1, double y1, double z1, double x2, double y2, double z2) {
  try{
    return {
        y1 * z2 - z1 * y2,
        z1 * x2 - x1 * z2,
        x1 * y2 - y1 * x2
    };
  }catch(std::exception &e){
    std::cerr << "An error occurred in cross(): " << e.what() << std::endl;
    return {0.0, 0.0, 0.0};
  }
}

std::pair<std::vector<double>, std::vector<double> > load_from_file(const std::string &filename) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error("Failed to open file: " + filename);
  }

  std::vector<double> y, u;
  double y_val, u_val;
  while (file >> u_val >> y_val) {
    y.push_back(y_val);
    u.push_back(u_val);
  }

  return std::make_pair(y, u);

}

//new function get_fil_coords2 which is more generalized and takes start and end points as input:
 std::vector<std::vector<double>> get_fil_coords2(bool fil_type, double x1, double y1, double z1, double x2, double y2, double z2, double dz_i) {
  std::vector<std::vector<double>> coords;
  
  if(fil_type){
    //for straight line between two points:
    //first, calculate total length of line:
    double line_length = norm(x2-x1, y2-y1, z2-z1);
    //then, calculate number of nodes needed:
    int node_count = std::round(line_length/dz_i);
    //then, calculate increment in each direction:
    double x_inc = (x2-x1)/node_count;
    double y_inc = (y2-y1)/node_count;
    double z_inc = (z2-z1)/node_count;
    //then, calculate coordinates for each node:
    for(int i=0; i<node_count; i++){
      std::vector<double> temp_coords {x1+i*x_inc, y1+i*y_inc, z1+i*z_inc};
      coords.push_back(temp_coords);
    }
  }else{    //for circular filaments
    std::cout << "Error - circular filaments not yet implemented in get_fil_coords2()...try get_fil_coords()\n";
  }
  return coords;
}

std::vector<std::vector<double>> get_fil_coords(bool fil_type, double x, double y, double domain_width, double dz_i) {
  std::vector<std::vector<double>> coords;
  
  if(fil_type){
    //first, for straight, z oriented line, starting at z=0:
    for(double z=0; z<domain_width+dz_i/2; z+=dz_i){
      std::vector<double> temp_coords {x, y, z};
      coords.push_back(temp_coords);
    }
  }else{    //for circular filaments.  Here assumed to be oriented in x direction and centered at z=0;
    //initial node count based on desired initial spacing dz_i and total circumference:
      //(domain_width taken as diameter for circles):
    double node_count1 = std::round(M_PI*domain_width/dz_i);
    //std::cout << "node count: " << node_count1 << "\n";
    //increment in radians (done in two steps to ensure even node spacing):
    double rad_inc = 2*M_PI/node_count1;
    //establish node coordinates, evenly spaced around circular circumference(including repeating end point):
    for(double r=0; r<2*M_PI+rad_inc/2; r+=rad_inc){
      double y_temp = y + domain_width*sin(r)/2;
      double z_temp = domain_width*cos(r)/2;

      std::vector<double> temp_coords = {x, y_temp, z_temp};
      coords.push_back(temp_coords);
    }
  }

  return coords;
}

void calculate_induced_velocity(Node& control_node, std::vector<std::unique_ptr<Filament>>& filaments, double& ivx, double& ivy, double& ivz, int node_num) {

  //create and open log file for each individual induced velocity (used for troubleshooting):
  //std::ofstream log_file;
  //log_file.open("log_file.txt", std::ios_base::app);

  // Initialize induced velocity to zero
  ivx = ivy = ivz = 0.0;

  //domain width for periodic boundary conditions (to be coded into config file)
  double width = 0.04;
  //periodic boolean, again to be coded into config file later
  bool peri = 1;

  for (std::unique_ptr<Filament>& filament : filaments) {

    //log_file << "Calculating induced velocity on node " << node_num << " (" << control_node.x << ", " 
              //<< control_node.y << ", " << control_node.z << ") of " << filament->get_node_count() << "...\n";

    Node* vortex_node = filament->getHead();
    int loop_count = 0;
    while (vortex_node != nullptr && loop_count < filament->get_node_count()) {
      
        // Calculate distance between control node and segment
        double rx = control_node.x - vortex_node->vcx;
        double ry = control_node.y - vortex_node->vcy;
        double rz = control_node.z - vortex_node->vcz;
        double rl = norm(rx, ry, rz);


        // Calculate k value
        double k = filament->get_circulation() / (4 * M_PI * pow(rl, 3));

          // Calculate cross product of segment and distance
        std::array<double, 3> cross_result = cross(vortex_node->sx, vortex_node->sy, vortex_node->sz,
              rx, ry, rz);
        double u = cross_result[0];
        double v = cross_result[1];
        double w = cross_result[2];

      
        u *= k;
        v *= k;
        w *= k;

        //log_file << " n(" << vortex_node->vcx << ", " << vortex_node->vcy << ", " << vortex_node->vcz << "): V(" << u << ", " << v << ", " << w << ")";

      if (filament->influenced_by_segment(control_node.x, control_node.y, control_node.z, vortex_node)) {
        // Update induced velocity
        ivx += u;
        ivy += v;
        ivz += w;
        //log_file << " (influenced)" << std::endl;
      }

      //if calculating including periodic boundary conditions:
      if(peri == 1){
        //repeat same process for one domain on either side:
        double rz_p1 = rz + width;
        double rl_p1 = norm(rx, ry, rz_p1);

        if(rl_p1 < filament->get_max_range() && rl_p1 > vortex_node->r){
          k = filament->get_circulation() / (4 * M_PI * pow(rl_p1, 3));
          cross_result = cross(vortex_node->sx, vortex_node->sy, vortex_node->sz,
            rx, ry, rz_p1);
          u = cross_result[0]*k;    //times k here to save lines for readability
          v = cross_result[1]*k;
          w = cross_result[2]*k;

          ivx += u;   //one domain at a time to reuse variables
          ivy += v;
          ivz += w;

          //log_file << "   n(" << vortex_node->vcx << ", " << vortex_node->vcy << ", " << rz_p1 << "): V(" << u << ", " << v << ", " << w << ")\n";
        }
          
        
        double rz_p2 = rz - width;
        double rl_p2 = norm(rx, ry, rz_p2);

        if(rl_p2 < filament->get_max_range() && rl_p2 > vortex_node->r){
          k = filament->get_circulation() / (4 * M_PI * pow(rl_p2, 3));
          cross_result = cross(vortex_node->sx, vortex_node->sy, vortex_node->sz,
            rx, ry, rz_p2);
          u = cross_result[0]*k; 
          v = cross_result[1]*k;
          w = cross_result[2]*k;

          ivx += u;  
          ivy += v;
          ivz += w;

          //log_file << "   n(" << vortex_node->vcx << ", " << vortex_node->vcy << ", " << rz_p2 << "): V(" << u << ", " << v << ", " << w << ")\n";
        }

        //repeating domains again for 5 domain total span, as a test
        double rz_p3 = rz + 2*width;
        double rl_p3 = norm(rx, ry, rz_p3);

        if(rl_p3 < filament->get_max_range() && rl_p3 > vortex_node->r){
          k = filament->get_circulation() / (4 * M_PI * pow(rl_p3, 3));
          cross_result = cross(vortex_node->sx, vortex_node->sy, vortex_node->sz,
            rx, ry, rz_p3);
          u = cross_result[0]*k;    //*k here to save lines for readability
          v = cross_result[1]*k;
          w = cross_result[2]*k;

          ivx += u;   //one domain at a time to reuse variables
          ivy += v;
          ivz += w;
        }
          
        
        double rz_p4 = rz - 2*width;
        double rl_p4 = norm(rx, ry, rz_p4);

        if(rl_p4 < filament->get_max_range() && rl_p4 > vortex_node->r){
          k = filament->get_circulation() / (4 * M_PI * pow(rl_p4, 3));
          cross_result = cross(vortex_node->sx, vortex_node->sy, vortex_node->sz,
            rx, ry, rz_p4);
          u = cross_result[0]*k; 
          v = cross_result[1]*k;
          w = cross_result[2]*k;

          ivx += u;  
          ivy += v;
          ivz += w;
        }

      }
    
      vortex_node = vortex_node->next;
      loop_count++;
    }
  }

  //close log file
  //log_file.close();
  
}

void move_nodes_via_induction(std::vector<std::unique_ptr<Filament>>& filaments, double dt) {
    int fil_count = 1;
    for (std::unique_ptr<Filament>& filament : filaments) {
        std::cout << "    On filament " << fil_count++ << " of " << filaments.size() << "...\n";
        Node* ctrl_node = filament->getHead();
        int node_count = 0;
        while (ctrl_node != nullptr && node_count < filament->get_node_count()) {
            // Calculate background velocity(&temporarily replace total v with it)
            auto evaluated_result = filament->background_velocity.evaluate(ctrl_node->x, ctrl_node->y, ctrl_node->z, filament->get_timestep());
            ctrl_node->vx = std::get<0>(evaluated_result);
            ctrl_node->vy = std::get<1>(evaluated_result);
            ctrl_node->vz = std::get<2>(evaluated_result);

            // Calculate induced velocity
            double ivx, ivy, ivz; 
            calculate_induced_velocity(*ctrl_node, filaments, ivx, ivy, ivz, node_count);

            // Update induced velocity (ramped for ends of braid vortices)
            //double c_ramp = 0.5;    //ramp coefficient for end values
            if( (fil_count%(filaments.size()/7)) < 8 ){   //if filament belongs to KH, update as usual      ***HARD CODED # BILLOWS***
              ctrl_node->ivx = ivx;
              ctrl_node->ivy = ivy;
              ctrl_node->ivz = ivz;

            /*COMMENTED OUT UNTIL ISSUE IS FIXED (SEE DEC1/4 RESULTS)     ---***---HERE:
            }else if( (node_count/(1.0*(filament->get_node_count()))) < 0.05 ){   //if node is in first 5% of filament, ramp to 50% by tip
              ctrl_node->ivx = ivx*(c_ramp + (1-c_ramp)*node_count/(0.05*filament->get_node_count()));
              ctrl_node->ivy = ivy*(c_ramp + (1-c_ramp)*node_count/(0.05*filament->get_node_count()));
              ctrl_node->ivz = ivz*(c_ramp + (1-c_ramp)*node_count/(0.05*filament->get_node_count()));
            }else if( (node_count/(1.0*(filament->get_node_count()))) > 0.95 ){   //if node is in last 5% of filament, ramp down to c_ramp (50% rn)
              ctrl_node->ivx = ivx*(1 - (c_ramp-1)*(node_count/(filament->get_node_count())-0.95)/0.05);
              ctrl_node->ivy = ivy*(1 - (c_ramp-1)*(node_count/(filament->get_node_count())-0.95)/0.05);
              ctrl_node->ivz = ivz*(1 - (c_ramp-1)*(node_count/(filament->get_node_count())-0.95)/0.05);
            */

              //ctrl_node->ivx = ivx*((filament->get_node_count()-node_count)/(0.05*filament->get_node_count()));
              //ctrl_node->ivy = ivy*((filament->get_node_count()-node_count)/(0.05*filament->get_node_count()));
              //ctrl_node->ivz = ivz*((filament->get_node_count()-node_count)/(0.05*filament->get_node_count()));
            }else{    //if node is in middle 90% of filament
              ctrl_node->ivx = ivx;
              ctrl_node->ivy = ivy;
              ctrl_node->ivz = ivz;
            }

            // Update total velocity
            ctrl_node->vx += ctrl_node->ivx;
            ctrl_node->vy += ctrl_node->ivy;
            ctrl_node->vz += ctrl_node->ivz;

            // Update position
            ctrl_node->x += ctrl_node->vx * dt;
            ctrl_node->y += ctrl_node->vy * dt;
            ctrl_node->z += ctrl_node->vz * dt;

            ctrl_node = ctrl_node->next;
            node_count++;
        }
        //bookkeeping
        filament->setTimestep(filament->get_timestep()+1);
        filament->setT(filament->get_time()+dt);
    } 
}

void step(std::vector<std::unique_ptr<Filament>>& filaments, double dt) {
    // Compute the vortex center and lengths for the current timestep.
    // This is done before movement so the vortex segments are fixed when
    // induced velocities.
    std::cout << "Stepping on from timestep " << filaments[0]->get_timestep() << "...\n"; //first vortex (0) chosen for output as they should all be at the same timestep (so it doesn't matter)
    // Also splits a segment into two if it grows too long
    int fil_count = 0;
    for (std::unique_ptr<Filament>& filament : filaments) {
        std::cout << "  On filament " << fil_count << "...\n";
        int ns = filament->get_nodes_added();
        std::cout << "    Computing segments and repartitioning if needed...\n";
        filament->compute_segments_and_repartition_if_needed();
        //if >0 segments split, go back through to see if any need splitting again
        std::cout << "    Done computing first segments and repartitioning if needed...\n";
        while (filament->get_nodes_added() > ns){
            ns = filament->get_nodes_added();
            filament->compute_segments_and_repartition_if_needed();
        }
        std::cout << "    Done computing segments and repartitioning if needed...\n";
        fil_count++;
    }
    
    std::cout << "  Moving nodes via induction...\n";
    // Actually move the nodes via induction
    move_nodes_via_induction(filaments, dt);
    //old, to be deleted:
    std::cout << "Now on timestep " << filaments[0]->get_timestep() << "...\n";
}

/*OLD VERSION FOR VECTOR OF FILS, NOT POINTERS:
std::string save_multiple_fil_timestep(std::string fn_template, const std::vector<Filament>& filaments) {
    int timestep_id = filaments[0].get_timestep();
    std::cout << "Saving timestep " << timestep_id << "...\n";

    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();

    std::vector<double> point_velocities;
    std::vector<double> point_induced_velocities;
    std::vector<double> point_circulation;
    std::vector<double> point_core_radiuses;
    std::vector<int> point_filament_ids;

    for (int f_idx = 0; f_idx < filaments.size(); ++f_idx) {
        const Filament& f = filaments[f_idx];

        int node_count = f.get_node_count();
        vtkIdType* node_ids_vtk = new vtkIdType[node_count];

        std::cout << "node_count = " << node_count << "...\n";

        Node* node = f.getHead();
        if (node == nullptr) {
              std::cout << "getHead node is null" << std::endl;
              break;
            }

        for (int i = 0; i < node_count; i++) {

          std::cout << "i = " << i << " OK start...\n";

            node_ids_vtk[i] = points->InsertNextPoint(node->x, node->y, node->z);
            node = node->next;

            if (node == nullptr) {    //this is catching something that needs fixing, but the catch allows it to work for now...
              std::cout << "node is null" << std::endl;
              break;
            }
            std::cout << "i = " << i << " OK finish...\n";
        }

        std::cout << "OK1\n";

        cells->InsertNextCell(node_count, node_ids_vtk);
        delete[] node_ids_vtk;

        std::cout << "OK2\n";

        node = f.getHead();
        for (int i = 0; i < node_count; i++) {
            point_velocities.push_back(node->vx);
            point_velocities.push_back(node->vy);
            point_velocities.push_back(node->vz);

            point_induced_velocities.push_back(node->ivx);
            point_induced_velocities.push_back(node->ivy);
            point_induced_velocities.push_back(node->ivz);

            point_circulation.push_back(f.get_circulation());

            point_core_radiuses.push_back(node->r);

            point_filament_ids.push_back(f_idx);

            node = node->next;
        }
    }

    std::cout << "OK3\n";

    vtkSmartPointer<vtkUnstructuredGrid> ugrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    ugrid->SetPoints(points);
    ugrid->SetCells(VTK_LINE, cells);

    vtkSmartPointer<vtkIntArray> filament_id_array = vtkSmartPointer<vtkIntArray>::New();
    filament_id_array->SetName("FilamentId");
    filament_id_array->SetNumberOfComponents(1);
    filament_id_array->SetNumberOfTuples(point_filament_ids.size());
    for (int i = 0; i < point_filament_ids.size(); ++i) {
        filament_id_array->SetValue(i, point_filament_ids[i]);
    }

    vtkSmartPointer<vtkDoubleArray> velocity_array = vtkSmartPointer<vtkDoubleArray>::New();
    velocity_array->SetNumberOfComponents(3);
    velocity_array->SetNumberOfTuples(point_velocities.size() / 3);
    velocity_array->SetName("Velocity");
    velocity_array->SetArray(point_velocities.data(), point_velocities.size(), 1);

    vtkSmartPointer<vtkDoubleArray> induced_velocity_array = vtkSmartPointer<vtkDoubleArray>::New();
    induced_velocity_array->SetNumberOfComponents(3);
    induced_velocity_array->SetNumberOfTuples(point_induced_velocities.size() / 3);
    induced_velocity_array->SetName("InducedVelocity");
    induced_velocity_array->SetArray(point_induced_velocities.data(), point_induced_velocities.size(), 1);

    vtkSmartPointer<vtkDoubleArray> circulation_array = vtkSmartPointer<vtkDoubleArray>::New();
    circulation_array->SetName("Circulation");
    circulation_array->SetNumberOfComponents(1);
    circulation_array->SetNumberOfTuples(point_circulation.size());
    circulation_array->SetArray(point_circulation.data(), point_circulation.size(), 1);

    vtkSmartPointer<vtkDoubleArray> core_radius_array = vtkSmartPointer<vtkDoubleArray>::New();
    core_radius_array->SetName("CoreRadius");
    core_radius_array->SetNumberOfComponents(1);
    core_radius_array->SetNumberOfTuples(point_core_radiuses.size());
    core_radius_array->SetArray(point_core_radiuses.data(), point_core_radiuses.size(), 1);

    ugrid->GetPointData()->AddArray(filament_id_array);
    ugrid->GetPointData()->AddArray(velocity_array);
    ugrid->GetPointData()->AddArray(induced_velocity_array);
    ugrid->GetPointData()->AddArray(circulation_array);
    ugrid->GetPointData()->AddArray(core_radius_array);

    std::stringstream ss;
    ss << fn_template << "_" << std::setw(7) << std::setfill('0') << timestep_id << ".vtu";

    std::string fn = ss.str();

    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(fn.c_str());
    writer->SetInputData(ugrid);
    writer->Write();

    std::cout << "Saved to file " << fn << std::endl;

    return fn;
}
*/

//updated version for unique_ptr
std::string save_multiple_fil_timestep(std::string fn_template, const std::vector<std::unique_ptr<Filament>>& filaments) {
    //int timestep_id = filaments[0]->get_timestep();
    std::cout << "Saving timestep...\n";

    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();

    std::vector<double> point_velocities;
    std::vector<double> point_induced_velocities;
    std::vector<double> point_circulation;
    std::vector<double> point_core_radiuses;
    std::vector<int> point_filament_ids;

    for (int f_idx = 0; f_idx < filaments.size(); ++f_idx) {
        const std::unique_ptr<Filament>& f = filaments[f_idx];

        int node_count = f->get_node_count();
        vtkIdType* node_ids_vtk = new vtkIdType[node_count];

        std::cout << "Fil " << f_idx << " node_count = " << node_count << "...\n";

        Node* node = f->getHead();
        if (node == nullptr) {
              std::cout << "getHead node is null" << std::endl;
              break;
            }

        for (int i = 0; i < node_count; i++) {

            node_ids_vtk[i] = points->InsertNextPoint(node->x, node->y, node->z);
            node = node->next;

            if (node == nullptr) {    //this catches nullptr on last node of filament
              //std::cout << "node is null" << std::endl;
              break;
            }
        }

        //cell is filament in this case
        cells->InsertNextCell(node_count, node_ids_vtk);
        delete[] node_ids_vtk;

        node = f->getHead();
        for (int i = 0; i < node_count; i++) {
            point_velocities.push_back(node->vx);
            point_velocities.push_back(node->vy);
            point_velocities.push_back(node->vz);

            point_induced_velocities.push_back(node->ivx);
            point_induced_velocities.push_back(node->ivy);
            point_induced_velocities.push_back(node->ivz);

            point_circulation.push_back(f->get_circulation());

            point_core_radiuses.push_back(node->r);

            point_filament_ids.push_back(f_idx);

            node = node->next;
        }
    }

    vtkSmartPointer<vtkUnstructuredGrid> ugrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    ugrid->SetPoints(points);
    ugrid->SetCells(VTK_LINE, cells);

    vtkSmartPointer<vtkIntArray> filament_id_array = vtkSmartPointer<vtkIntArray>::New();
    filament_id_array->SetName("FilamentId");
    filament_id_array->SetNumberOfComponents(1);
    filament_id_array->SetNumberOfTuples(point_filament_ids.size());
    for (int i = 0; i < point_filament_ids.size(); ++i) {
        filament_id_array->SetValue(i, point_filament_ids[i]);
    }

    vtkSmartPointer<vtkDoubleArray> velocity_array = vtkSmartPointer<vtkDoubleArray>::New();
    velocity_array->SetNumberOfComponents(3);
    velocity_array->SetNumberOfTuples(point_velocities.size() / 3);
    velocity_array->SetName("Velocity");
    velocity_array->SetArray(point_velocities.data(), point_velocities.size(), 1);

    vtkSmartPointer<vtkDoubleArray> induced_velocity_array = vtkSmartPointer<vtkDoubleArray>::New();
    induced_velocity_array->SetNumberOfComponents(3);
    induced_velocity_array->SetNumberOfTuples(point_induced_velocities.size() / 3);
    induced_velocity_array->SetName("InducedVelocity");
    induced_velocity_array->SetArray(point_induced_velocities.data(), point_induced_velocities.size(), 1);

    vtkSmartPointer<vtkDoubleArray> circulation_array = vtkSmartPointer<vtkDoubleArray>::New();
    circulation_array->SetName("Circulation");
    circulation_array->SetNumberOfComponents(1);
    circulation_array->SetNumberOfTuples(point_circulation.size());
    circulation_array->SetArray(point_circulation.data(), point_circulation.size(), 1);

    vtkSmartPointer<vtkDoubleArray> core_radius_array = vtkSmartPointer<vtkDoubleArray>::New();
    core_radius_array->SetName("CoreRadius");
    core_radius_array->SetNumberOfComponents(1);
    core_radius_array->SetNumberOfTuples(point_core_radiuses.size());
    core_radius_array->SetArray(point_core_radiuses.data(), point_core_radiuses.size(), 1);

    ugrid->GetPointData()->AddArray(filament_id_array);
    ugrid->GetPointData()->AddArray(velocity_array);
    ugrid->GetPointData()->AddArray(induced_velocity_array);
    ugrid->GetPointData()->AddArray(circulation_array);
    ugrid->GetPointData()->AddArray(core_radius_array);

    std::stringstream ss;
    ss << fn_template << ".vtu";

    std::string fn = ss.str();

    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(fn.c_str());
    writer->SetInputData(ugrid);
    writer->Write();

    std::cout << "Saved to file " << fn << std::endl;

    return fn;
}

//Function to read vtk output from above and intialize filaments vector for simulation restarts:

/*//original try
void read_vtk(std::string fn, std::vector<std::unique_ptr<Filament>>& filaments) {
    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName(fn.c_str());
    reader->Update();

    vtkSmartPointer<vtkUnstructuredGrid> ugrid = reader->GetOutput();

    vtkSmartPointer<vtkPoints> points = ugrid->GetPoints();
    vtkSmartPointer<vtkCellArray> cells = ugrid->GetCells();

    vtkSmartPointer<vtkIntArray> filament_id_array = vtkSmartPointer<vtkIntArray>::New();
    filament_id_array = vtkIntArray::SafeDownCast(ugrid->GetPointData()->GetArray("FilamentId"));

    vtkSmartPointer<vtkDoubleArray> velocity_array = vtkSmartPointer<vtkDoubleArray>::New();
    velocity_array = vtkDoubleArray::SafeDownCast(ugrid->GetPointData()->GetArray("Velocity"));

    vtkSmartPointer<vtkDoubleArray> induced_velocity_array = vtkSmartPointer<vtkDoubleArray>::New();
    induced_velocity_array = vtkDoubleArray::SafeDownCast(ugrid->GetPointData()->GetArray("InducedVelocity"));

    vtkSmartPointer<vtkDoubleArray> circulation_array = vtkSmartPointer<vtkDoubleArray>::New();
    circulation_array = vtkDoubleArray::SafeDownCast(ugrid->GetPointData()->GetArray("Circulation"));

    vtkSmartPointer<vtkDoubleArray> core_radius_array = vtkSmartPointer<vtkDoubleArray>::New();
    core_radius_array = vtkDoubleArray::SafeDownCast(ugrid->GetPointData()->GetArray("CoreRadius"));

    int num_cells = cells->GetNumberOfCells();
    int num_points = points->GetNumberOfPoints();

    std::cout << "num_cells = " << num_cells << std::endl;
    std::cout << "num_points = " << num_points << std::endl;

    int cell_count = 0;
    int point_count = 0;

    for (int i = 0; i < num_cells; i++) {
        vtkSmartPointer<vtkIdList> cell = vtkSmartPointer<vtkIdList>::New();
        cells->GetNextCell(cell);

        int num_points_in_cell = cell->GetNumberOfIds();

        std::cout << "cell " << i << " has " << num_points_in_cell << " points" << std::endl;
        //finish above function to read in filaments from vtk file:
        std::unique_ptr<Filament> filament = std::make_unique<Filament>(num_points_in_cell);
        for (int j = 0; j < num_points_in_cell; j++) {
            int point_id = cell->GetId(j);

            double* coords = points->GetPoint(point_id);
            double x = coords[0];
            double y = coords[1];
            double z = coords[2];

            double velocity = velocity_array->GetValue(point_count);
            double induced_velocity = induced_velocity_array->GetValue(point_count);
            double circulation = circulation_array->GetValue(point_count);
            double core_radius = core_radius_array->GetValue(point_count);

            int filament_id = filament_id_array->GetValue(point_count);

            if (j == 0) {
                filament->initialize(x, y, z, velocity, induced_velocity, circulation, core_radius, filament_id, true);
            } else if (j == num_points_in_cell - 1) {
                filament->initialize(x, y, z, velocity, induced_velocity, circulation, core_radius, filament_id, false);
            } else {
                filament->add_node(x, y, z, velocity, induced_velocity, circulation, core_radius, filament_id);
            }

            point_count++;
        }

        filaments.emplace_back(std::move(filament));
        cell_count++;
    }
}
*/

void read_vtk(std::string fn, std::vector<std::unique_ptr<Filament>>& filaments, double circulation, double core_radius, double max_sl, double min_sl, double max_range, SineBackgroundGenerator background_vel) {
    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName(fn.c_str());
    reader->Update();

    vtkSmartPointer<vtkUnstructuredGrid> ugrid = reader->GetOutput();

    vtkSmartPointer<vtkPoints> points = ugrid->GetPoints();
    vtkSmartPointer<vtkCellArray> cells = ugrid->GetCells();

    vtkSmartPointer<vtkIntArray> filament_id_array = vtkSmartPointer<vtkIntArray>::New();
    filament_id_array = vtkIntArray::SafeDownCast(ugrid->GetPointData()->GetArray("FilamentId"));

    vtkSmartPointer<vtkDoubleArray> velocity_array = vtkSmartPointer<vtkDoubleArray>::New();
    velocity_array = vtkDoubleArray::SafeDownCast(ugrid->GetPointData()->GetArray("Velocity"));

    vtkSmartPointer<vtkDoubleArray> induced_velocity_array = vtkSmartPointer<vtkDoubleArray>::New();
    induced_velocity_array = vtkDoubleArray::SafeDownCast(ugrid->GetPointData()->GetArray("InducedVelocity"));

    vtkSmartPointer<vtkDoubleArray> circulation_array = vtkSmartPointer<vtkDoubleArray>::New();
    circulation_array = vtkDoubleArray::SafeDownCast(ugrid->GetPointData()->GetArray("Circulation"));

    vtkSmartPointer<vtkDoubleArray> core_radius_array = vtkSmartPointer<vtkDoubleArray>::New();
    core_radius_array = vtkDoubleArray::SafeDownCast(ugrid->GetPointData()->GetArray("CoreRadius"));

    int num_cells = cells->GetNumberOfCells();
    int num_points = points->GetNumberOfPoints();

    std::cout << "num_cells = " << num_cells << std::endl;
    std::cout << "num_points = " << num_points << std::endl;

    int point_count = 0;
    std::vector<std::vector<double>> filament_node_data;

    for (int i = 0; i < num_cells; i++) {
        vtkSmartPointer<vtkIdList> cell = vtkSmartPointer<vtkIdList>::New();
        cells->GetNextCell(cell);

        int num_points_in_cell = cell->GetNumberOfIds();

        std::vector<std::vector<double>> node_coords;

        std::cout << "cell " << i << " has " << num_points_in_cell << " points" << std::endl;

        std::vector<double> node_data;
        for (int j = 0; j < num_points_in_cell; j++) {
            int point_id = cell->GetId(j);

            double* coords = points->GetPoint(point_id);
            double x = coords[0];
            double y = coords[1];
            double z = coords[2];
            
            std::vector<double> temp_coords = {x, y, z};
            node_coords.push_back(temp_coords);
            /*//don't think I need any of this since initialization will figure it all out, 
            //just need current coords for initialization
            double velocity = velocity_array->GetValue(point_count);
            double induced_velocity = induced_velocity_array->GetValue(point_count);
            double circulation = circulation_array->GetValue(point_count);
            double core_radius = core_radius_array->GetValue(point_count);

            int filament_id = filament_id_array->GetValue(point_count);
            
            node_data.push_back(x);
            node_data.push_back(y);
            node_data.push_back(z);
            node_data.push_back(velocity);
            node_data.push_back(induced_velocity);
            node_data.push_back(circulation);
            node_data.push_back(core_radius);
            node_data.push_back(filament_id);

            */

            point_count++;
        }

        auto temp_f = std::make_unique<Filament>(circulation, core_radius, max_sl, min_sl, max_range, node_coords, background_vel);
        filaments.push_back(std::move(temp_f));
        //ideally lookup/savetimestep:
    }
}


/*another failed attempt at finishing:
    for (const auto& node_data : filament_node_data) {
    int filament_id = node_data.back();
    double circulation = node_data.at(5);
    double core_radius = node_data.at(6);
    SineBackgroundGenerator background_velocity;

    // Get all node coordinates
    std::vector<std::vector<double>> coords;
    for (int i = 0; i < node_data.size() - 3; i += 3) {
        coords.push_back({node_data.at(i), node_data.at(i + 1), node_data.at(i + 2)});
    }

    // Create new filament object
    std::unique_ptr<Filament> filament = std::make_unique<Filament>(
        circulation, core_radius, 0.0, 0.0, 0.0, coords, background_velocity, true);

    // Add filament object to the vector
    bool added_to_existing_filament = false;
    for (const auto& existing_filament : filaments) {
        if (existing_filament->get_id() == filament_id) {
            existing_filament->add_node(filament->get_nodes()[0]);
            added_to_existing_filament = true;
            break;
        }
    }
    if (!added_to_existing_filament) {
        filaments.push_back(std::move(filament));
    }
}
*/




Node* new_node(double x, double y, double z, double r, Node* prev, Node* next) {
  Node* node = new Node();
  node->born_on = 0;
  node->x = x;
  node->y = y;
  node->z = z;
  node->vx = 0.0;
  node->vy = 0.0;
  node->vz = 0.0;
  node->ivx = 0.0;
  node->ivy = 0.0;
  node->ivz = 0.0;

  // We cannot compute most of these as the filament needs to be completely
  // constructed before they can be built.
  node->x_prev = x;
  node->y_prev = y;
  node->z_prev = z;
  node->vcx = 0.0;
  node->vcy = 0.0;
  node->vcz = 0.0;
  node->sx = 0.0;
  node->sy = 0.0;
  node->sz = 0.0;
  node->sl = 0.0;
  node->vol = -1.0;
  node->r = r;
  node->prev = prev;
  node->next = next;
  return node;
}

//sinebackground member functions
SineBackgroundGenerator::SineBackgroundGenerator(std::string background_velocity_file,
                          std::vector<double> amplitudes,
                          std::vector<double> wavelengths,
                          std::vector<double> phases,
                          std::vector<int> directions)
      : amplitudes(amplitudes),
        wavelengths(wavelengths),
        phases(phases),
        directions(directions),
        n(directions.size()), 
        bottom_ufunc() {
    data = load_from_file(background_velocity_file);
    std::vector<double> y = data.first;
    std::vector<double> u = data.second;
    if (y.empty() || u.empty() || y.size() != u.size()) {
        throw std::runtime_error("Error: y and u vectors must be non-empty and the same size.");
    }
    /*for(int i=0; i<y.size(); i++){
        std::cout << y[i] << ", " << u[i] << "\n" ;
    }*/
    //tk::spline bottom_ufunc(y, u);
    bottom_ufunc.set_points(y, u);
    /*double test_y = 0.003;
    double test_val = bottom_ufunc(test_y); 
    std::cout << test_val << "\n" ;*/
}

std::array<double, 3> SineBackgroundGenerator::evaluate(double x, double y, double z, int t) {
    double u = bottom_ufunc(y) ;
    //std::cout << "y: " << y << std::endl;
    //std::cout << "u: " << u << std::endl;
    std::array<double, 3> v = {u, 0.0, 0.0};
    for (int i = 0; i < n; i++) {
      //std::cout << "amplitudes[" << i << "]: " << amplitudes[i] << std::endl;
      //std::cout << "wavelengths[" << i << "]: " << wavelengths[i] << std::endl;
      //std::cout << "phases[" << i << "]: " << phases[i] << std::endl;
      //std::cout << "directions[" << i << "]: " << directions[i] << std::endl;
      /*if(z<0.01){       //trying out a 'fade in/out' to limit end effects
        v[directions[i]] += amplitudes[i] * (z/0.01)*sin(
          2 * M_PI / wavelengths[i] * (z + phases[i])
        );
      }else if(z>0.03){
        v[directions[i]] += amplitudes[i] * ((0.04-z)/0.01)*sin(
          2 * M_PI / wavelengths[i] * (z + phases[i])
        );
      }else{
        v[directions[i]] += amplitudes[i] * sin(
          2 * M_PI / wavelengths[i] * (z + phases[i])
        );
      }*/
      v[directions[i]] += amplitudes[i] * sin(
        2 * M_PI / wavelengths[i] * (z + phases[i])
      );
    }
    //std::cout << "n: " << n << std::endl;
    return v;
}

// Filament member function implementations
Filament::Filament(double circulation, double core_radius, double max_sl, double min_sl, double max_range,
           std::vector<std::vector<double> > coords, SineBackgroundGenerator background_velocity, bool straight = 0)
    : timestep(0), t(0), nodes_added(0), nodes_deleted(0) {
    if (coords.size() < 2) {
      throw std::runtime_error("Must have at least two nodes");
    }

    // Create a file object and open the file in write mode (troubleshooting)
    //std::ofstream fout("constructor_test_output.txt");

    this->node_count = coords.size();
    //fout << "node count: " << this->node_count << "\n"; //troubleshooting
    this->circulation = circulation;
    //fout << "circ: " << this->circulation << "\n"; //troubleshooting
    this->core_radius = core_radius;
    //fout << "c_r: " << this->core_radius << "\n"; //troubleshooting
    this->max_sl = max_sl;
    this->min_sl = min_sl;
    this->max_range = max_range;
    this->background_velocity = background_velocity;

    double x = coords[0][0];
    double y = coords[0][1];
    double z = coords[0][2];

    this->head = new_node(x, y, z, core_radius, nullptr, nullptr);
    std::array<double, 3> velocity = this->background_velocity.evaluate(x, y, z, 0);
    this->head->vx = velocity[0];
    this->head->vy = velocity[1];
    this->head->vz = velocity[2];

    Node* prev_node = this->head;
    Node* current_node = nullptr;

    for (int i = 1; i < this->node_count; i++) {
      x = coords[i][0];
      y = coords[i][1];
      z = coords[i][2];

      current_node = new_node(x, y, z, core_radius, prev_node, nullptr);
      //fout << "(x, y, z): (" << current_node->x << ", " << current_node->y << ", " << current_node->z << ")\n"; //troubleshooting
      velocity = this->background_velocity.evaluate(x, y, z, 0);
      current_node->vx = velocity[0];
      current_node->vy = velocity[1];
      current_node->vz = velocity[2];

      // Hook up nodes and compute segment property
      prev_node->next = current_node;
      this->compute_single_segment_geometric_property(prev_node);
      prev_node->vol = M_PI * prev_node->r * prev_node->r * prev_node->sl;
      prev_node = current_node;
    }
    this->tail = current_node;

    //if circular, connect head and tail
    if (!straight) {
      this->tail->next = this->head;
      this->head->prev = this->tail;
      this->compute_single_segment_geometric_property(this->tail);
      this->tail->vol = M_PI * this->tail->r * this->tail->r * this->tail->sl;
    }

    //fout << "Tail: " << this->tail->x << ", " << this->tail->y << ", " << this->tail->z << "\n";
    //fout.close();
}

Filament::~Filament() {
    Node* current = head;
    while (current != nullptr) {
        Node* next = current->next;
        delete current;
        current = next;
    }
}

void Filament::compute_segments_and_repartition_if_needed() {
    Node* n1 = head;
    Node* n2 = head->next;
    while (n2 != nullptr && n2 != head) {
      compute_single_segment_geometric_property(n1);
      if (n1->sl > max_sl) {
        split_segment(n1, n2);
        nodes_added += 1;
      } else if (n1->sl < min_sl && n1 != head && n2 != tail) {
        merge_segment(n1, n2);
        nodes_deleted += 1;
      }
      n1 = n2;
      n2 = n2->next;
    }
}

void Filament::split_segment(Node* n1, Node* n2) {
    Node* nmid = new_node(n1->vcx, n1->vcy, n1->vcz, core_radius, n1, n2);
    nmid->born_on = timestep;

    n1->next = nmid;
    n2->prev = nmid;
    node_count++;

    compute_single_segment_geometric_property(n1);
    compute_single_segment_geometric_property(nmid);
}

void Filament::merge_segment(Node* n1, Node* n2) {
    n2->x = n1->vcx;
    n2->y = n1->vcy;
    n2->z = n1->vcz;
      
    n2->born_on = timestep;
      
    n2->prev = n1->prev;
    n2->prev->next = n2;
      
    node_count -= 1;
    
    delete n1;

    compute_single_segment_geometric_property(n2);
}

bool Filament::influenced_by_segment(double x, double y, double z, Node* n1) {
    double radius_of_influence, distance;

    radius_of_influence = sqrt(pow(n1->r, 2) + pow(n1->sl / 2, 2));
    distance = sqrt(pow(x - n1->vcx, 2) + pow(y - n1->vcy, 2) + pow(z - n1->vcz, 2));

    if (distance < radius_of_influence) {
      distance = perpendicular_distance_to_segment(x, y, z, n1);
    }

    return distance > n1->r && distance < max_range;
}

double Filament::perpendicular_distance_to_segment(double x, double y, double z, Node* n1) {
    double vec12x, vec12y, vec12z, vec01x, vec01y, vec01z;

    vec12x = n1->sx;
    vec12y = n1->sy;
    vec12z = n1->sz;

    vec01x = n1->x_prev - x;
    vec01y = n1->y_prev - y;
    vec01z = n1->z_prev - z;

    std::array<double, 3> cross_product = cross(vec12x, vec12y, vec12z, vec01x, vec01y, vec01z);

    return norm(cross_product[0], cross_product[1], cross_product[2]) / norm(vec12x, vec12y, vec12z);
}

void Filament::compute_single_segment_geometric_property(Node* n1) {
    Node* n2 = n1->next;
    if(n2 == nullptr){
      std::cout << "n2 is null\n";
      if(n1 == tail){
        std::cout << "n1 is tail\n";
      }
    }

    n1->vcx = (n1->x + n2->x) / 2.0;
    n1->vcy = (n1->y + n2->y) / 2.0;
    n1->vcz = (n1->z + n2->z) / 2.0;

    n1->sx = n2->x - n1->x;
    n1->sy = n2->y - n1->y;
    n1->sz = n2->z - n1->z;

    n1->sl = norm(n1->sx, n1->sy, n1->sz);
    
    n1->x_prev = n1->x;
    n1->y_prev = n1->y;
    n1->z_prev = n1->z;
}

std::string Filament::save_single_timestep(std::string fn_template) {
    std::cout << "Saving timestep " << timestep << "...\n";
    int i = 0;
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(node_count);

    vtkSmartPointer<vtkUnstructuredGrid> ugrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    ugrid->Allocate(1);

    std::vector<int> point_ids(node_count);
    std::vector<std::array<double, 3> > velocities(node_count);
    std::vector<std::array<double, 3> > induced_velocities(node_count);
    std::vector<double> core_radiuses(node_count);

    Node *current = head;
    while (current && i < node_count) {
      points->SetPoint(i, current->x, current->y, current->z);
      point_ids[i] = i;

      velocities[i][0] = current->vx;
      velocities[i][1] = current->vy;
      velocities[i][2] = current->vz;

      induced_velocities[i][0] = current->ivx;
      induced_velocities[i][1] = current->ivy;
      induced_velocities[i][2] = current->ivz;

      core_radiuses[i] = current->r;

      current = current->next;
      i++;
    }

    vtkSmartPointer<vtkDoubleArray> velocity_array = vtkSmartPointer<vtkDoubleArray>::New();
    velocity_array->SetNumberOfComponents(3);
    velocity_array->SetNumberOfTuples(node_count);
    velocity_array->SetName("Velocity");
    for (int i = 0; i < node_count; i++) {
      velocity_array->SetTuple(i, velocities[i].data());
    }

    vtkSmartPointer<vtkDoubleArray> induced_velocity_array = vtkSmartPointer<vtkDoubleArray>::New();
    induced_velocity_array->SetNumberOfComponents(3);
    induced_velocity_array->SetNumberOfTuples(node_count);
    induced_velocity_array->SetName("InducedVelocity");
    for (int i = 0; i < node_count; i++) {
      induced_velocity_array->SetTuple(i, induced_velocities[i].data());
    }

    vtkSmartPointer<vtkDoubleArray> circulation_array = vtkSmartPointer<vtkDoubleArray>::New();
    circulation_array->SetNumberOfTuples(node_count);
    circulation_array->SetName("Circulation");
    for (int i = 0; i < node_count; i++) {
      circulation_array->SetTuple1(i, circulation);
    }

    vtkSmartPointer<vtkDoubleArray> core_radius_array = vtkSmartPointer<vtkDoubleArray>::New();
    core_radius_array->SetNumberOfTuples(node_count);
    core_radius_array->SetName("CoreRadius");
    for (int i = 0; i < node_count; i++) {
      core_radius_array->SetTuple1(i, core_radiuses[i]);
    }

    //vtkSmartPointer<vtkUnstructuredGrid> ugrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    ugrid->SetPoints(points);
    vtkIdType* point_ids_vtk = new vtkIdType[node_count];
    for (int i = 0; i < node_count; i++) {
      point_ids_vtk[i] = point_ids[i];
    }

    ugrid->InsertNextCell(VTK_LINE, node_count, point_ids_vtk);
    delete [] point_ids_vtk;

    ugrid->GetPointData()->AddArray(velocity_array);
    ugrid->GetPointData()->AddArray(induced_velocity_array);
    ugrid->GetPointData()->AddArray(circulation_array);
    ugrid->GetPointData()->AddArray(core_radius_array);

    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    std::string fn = fn_template;
    fn = fn_template.substr(0, fn_template.size() - 4) + std::to_string(timestep) + ".vtu";
    writer->SetFileName(fn.c_str());
    writer->SetInputData(ugrid);
    writer->Write();

    return fn;
}

Node* Filament::getHead() const {
  return head;
}

Node* Filament::getTail() const {
  return tail;
}

int Filament::get_node_count() const {
  return node_count;
}

double Filament::get_circulation() const {
    return circulation;
}

int Filament::get_timestep() const {
    return timestep;
}

int Filament::get_nodes_added() const {
    return nodes_added;
}

void Filament::setTimestep(int newTimestep) {
  timestep = newTimestep;
}

void Filament::setT(double newT) {
  t = newT;
}

double Filament::get_time() const {
    return t;
}

double Filament::get_max_range() const {
    return max_range;
}

void Filament::og_step(double dt) {
    // Compute the vortex center and lengths for the current timestep.
    // This is done before movement so the vortex segments are fixed when
    // induced velocities.
    std::cout << "Stepping on from timestep " << timestep << "...\n";
    // Also splits a segment into two if it grows too long
    int ns = nodes_added;
    compute_segments_and_repartition_if_needed();
    //if >0 segments split, go back through to see if any need splitting again
    while (nodes_added>ns){
      ns = nodes_added;
      compute_segments_and_repartition_if_needed();
    }
    // Actually move the nodes via induction
    og_move_nodes_via_induction(dt);
    // Bookkeeping
    timestep++;
    t += dt;
    std::cout << "Now on timestep " << timestep << "...\n";
  }

void Filament::og_move_nodes_via_induction(double dt) {
    Node* ctrl_node = head;
    Node* vortex_node1;
    Node* vortex_node2;
    double rx, ry, rz, rl, k, u, v, w;

    while (ctrl_node != nullptr) {
        // Evaluate background velocity
        auto evaluated_result = background_velocity.evaluate(ctrl_node->x, ctrl_node->y, ctrl_node->z, timestep);
        ctrl_node->vx = std::get<0>(evaluated_result);
        ctrl_node->vy = std::get<1>(evaluated_result);
        ctrl_node->vz = std::get<2>(evaluated_result);

        ctrl_node->ivx = ctrl_node->ivy = ctrl_node->ivz = 0.0;

        vortex_node1 = head;
        vortex_node2 = head->next;
        while (vortex_node2 != nullptr) {
          if (influenced_by_segment(ctrl_node->x, ctrl_node->y, ctrl_node->z, vortex_node1)) {
              // Calculate distance between control node and segment
              rx = ctrl_node->x - vortex_node1->vcx;
              ry = ctrl_node->y - vortex_node1->vcy;
              rz = ctrl_node->z - vortex_node1->vcz;
              rl = norm(rx, ry, rz);

              // Calculate k value
              k = circulation / (4 * M_PI * pow(rl, 3));

              // Calculate cross product of segment and distance
              std::array<double, 3> cross_result = cross(vortex_node1->sx, vortex_node1->sy, vortex_node1->sz,
                    rx, ry, rz);
              u = cross_result[0];
              v = cross_result[1];
              w = cross_result[2];

              u *= k;
              v *= k;
              w *= k;

              // Update induced velocity
              ctrl_node->ivx += u;
              ctrl_node->ivy += v;
              ctrl_node->ivz += w;
          }

          vortex_node1 = vortex_node2;
          vortex_node2 = vortex_node2->next;
      }

      // Update velocity
      /*if(ctrl_node->z<0.01){      //introduce fade for induced vel. to try and limit end effects
        ctrl_node->vx += (ctrl_node->z/0.01)*ctrl_node->ivx;
        ctrl_node->vy += (ctrl_node->z/0.01)*ctrl_node->ivy;
        ctrl_node->vz += (ctrl_node->z/0.01)*ctrl_node->ivz;
      }else if(ctrl_node->z>0.03){
        ctrl_node->vx += ((0.04-ctrl_node->z)/0.01)*ctrl_node->ivx;
        ctrl_node->vy += ((0.04-ctrl_node->z)/0.01)*ctrl_node->ivy;
        ctrl_node->vz += ((0.04-ctrl_node->z)/0.01)*ctrl_node->ivz;
      }else{
        ctrl_node->vx += ctrl_node->ivx;
        ctrl_node->vy += ctrl_node->ivy;
        ctrl_node->vz += ctrl_node->ivz;
      }*/
      
      ctrl_node->vx += ctrl_node->ivx;
      ctrl_node->vy += ctrl_node->ivy;
      ctrl_node->vz += ctrl_node->ivz;

      // Update position
      ctrl_node->x += ctrl_node->vx * dt;
      ctrl_node->y += ctrl_node->vy * dt;
      ctrl_node->z += ctrl_node->vz * dt;

      ctrl_node = ctrl_node->next;
    }
  }