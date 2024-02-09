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
#include <array>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <vector>
//#include <armadillo>
#include <utility>
#include <map>
#include <tuple>
#include "/Users/alexmann/Library/Mobile Documents/com~apple~CloudDocs/Documents/School/Carleton/Masters/program_working_dir/inviscid_model/spline.h"
#include "/Users/alexmann/Library/Mobile Documents/com~apple~CloudDocs/Documents/School/Carleton/Masters/program_working_dir/inviscid_model/filament.h"


int main() {

  //ask user if they would like to 1) start a new run or 2) evaluate a velocity field for an old result:
  std::cout << "Would you like to 1) start a new run or 2) evaluate a velocity field for an old result?\n(Enter 1 or 2): ";
  int choice;
  std::cin >> choice;
  while (choice != 1 && choice != 2) {
    std::cout << "Invalid choice, please enter 1 or 2: ";
    std::cin >> choice;
  }

  std::ifstream configFile;   //define here so can be used in both cases
  std::string full_res_path;  //define here so can be used globablly 

  if(choice == 2){    //evaluate velocity field for old result ONLY
    //folder and name for results to evaluate (if not same folder):
    std::string res_folder = "RESULTS/5f_5kh_1b_jan/";
    full_res_path = res_folder + "5f_5kh_1b_timestep_40.vtu";
    std::string full_config_path = res_folder + "config.txt";
    //check if res_folder variable has been defined:
    if (res_folder.empty()) {
      configFile.open("config.txt");
      std::cout << "res_folder variable not defined, config in current directory opened..." << std::endl;
    }else{
      std::cout << "Config path: " << full_config_path << "\n";
      configFile.open(full_config_path);
      if (!configFile) {
        std::cerr << "Error opening config file(2)" << std::endl;
        return 1;
      }
      std::cout << "Config " << res_folder << "config.txt opened...\n";
    }

    std::cout << "Evaluating velocity field for result: " + full_res_path + "\n";
  }else{
    std::cout << "Starting new run...\n";
    configFile.open("config.txt");
    if (!configFile) {
      std::cerr << "Error opening config file(4)" << std::endl;
      return 1;
    }
    std::cout << "Config opened...\n";
  }

  //read config file in either case:
  std::string bg_vel_file; //= "/Users/alexmann/Documents/School/Carleton/Masters/program_working_dir/inviscid_model/vel_prof.txt" ;
  std::vector<double> amplitudes(3,0.0);  //= {0.02, 0.0, 0.0};
  std::vector<double> wavelengths(3,0.0);  //= {0.004, 1.0, 1.0};
  std::vector<double> phases(3,0.0);   //= {0.0, 0.0, 0.0};
  std::vector<int> directions(1,0);  //= {0};
  double KH_circulation, fil_circulation, core_radius, max_sl, min_sl, max_range, dz_i, dt, domain_width, braid_circ, braid_spacing1, braid_spacing2;
  int num_timesteps, num_billows, num_filaments, KH_num_steps, num_braid_fils;
  bool fil_type;
  bool initialize_fr_file;
  std::string initialize_file_path;
  
  //read setup parameters from config file
  std::string line;
  int line_count = 1;
  while(std::getline(configFile, line) && line_count < 20){
    
    if (line.empty() || line[0] == '#'){
      continue;
    }
    
    std::stringstream lineStream(line);

    switch(line_count){
      case 1:
        //lineStream >> bg_vel_file; //original but need to include space
        bg_vel_file = line;
        break;
      case 2:
        lineStream >> amplitudes[0] >> amplitudes[1] >> amplitudes[2];
        break;
      case 3:
        lineStream >> wavelengths[0] >> wavelengths[1] >> wavelengths[2];
        break;
      case 4:
        lineStream >> phases[0] >> phases[1] >> phases[2];
        break;
      case 5:
        lineStream >> directions[0];
        break;
      case 6:
        lineStream >> KH_circulation;
        break;
      case 7:
        lineStream >> core_radius >> max_range ;
        break;
      case 8:
        lineStream >> min_sl >> max_sl >> dz_i;
        break;
      case 9:
        lineStream >> domain_width;
        break;
      case 10:
        lineStream >> dt;
        break;
      case 11:
        lineStream >> num_timesteps;
        break;
      case 12:
        lineStream >> num_billows;
        break;
      case 13:
        lineStream >> num_filaments;
        break;
      case 14:
        lineStream >> fil_type;
        break;
      case 15:
        lineStream >> KH_num_steps;
        break;
      case 16:
        lineStream >> initialize_fr_file;
        break;
      case 17:
        initialize_file_path = line;
        break;
      case 18:
        lineStream >> braid_circ >> braid_spacing1 >> braid_spacing2;
        break;
      case 19:
        lineStream >> num_braid_fils;
        break;
    }

    //gotta get the node definitions inside here too to use the end points for filament initialization
    
    line_count++;

  }  

  std::vector<std::vector<double>> fil_initial_coords(num_filaments+num_braid_fils, std::vector<double>(2));

  //rewrite of initial filament-location-reading to allow for comments:
  int k=0;
  while(std::getline(configFile, line) && k<(num_filaments+num_braid_fils)){
    if (line.empty() || line[0] == '#'){
      continue;
    }

    std::stringstream lineStream(line);
    while(!lineStream.eof() && k<num_filaments){
      if(fil_type){
        //straight lines, just x, y for each fil, z_1 taken to be 0 later in implementation and z_2=domain_width
        lineStream >> fil_initial_coords[k][0];
        lineStream >> fil_initial_coords[k][1];
      }else{
        //in this case, circular filaments, x y of circular center, diameter taken as 'domain_width' and z assumed=0
        lineStream >> fil_initial_coords[k][0];
        lineStream >> fil_initial_coords[k][1];
      }
      k++;
    }
  }
/*
  std::cout << "Fil coords:\n";
  for(int i=0; i<num_filaments; i++){
    std::cout << fil_initial_coords[i][0] << " " << fil_initial_coords[i][1] <<"\n";
  }
*/
  std::cout << "Config file read finished...\n";

  configFile.close();

  SineBackgroundGenerator background_vel = SineBackgroundGenerator(bg_vel_file,
    amplitudes, wavelengths, phases, directions) ;

  //int num_gen = std::round(num_timesteps/25);
  //make sure initial increment is as close to dz_i as possible while still perfectly fitting domain width:
  int num_nodes_k = std::round(domain_width/dz_i);
  dz_i = domain_width/num_nodes_k;
  
  //create vector for all initial nodes on filament(s)
  std::vector<std::vector<std::vector<double>>> fil_node_coords(num_filaments + num_braid_fils);   //if fix end points, calc number of elements, can initialize better.
  //taking two endpoints and initializing filament(s) with evenly spaced nodes 
  //...(get_fil_coords only written for z-oriented filaments, get_fil_coords2 more generalized):
  if(choice==1){      //only need to do this if starting a new run
    for(int l=0; l<num_filaments; l++){
      std::cout << "filament " << l << "...\n";

      std::vector<std::vector<double>> temp_coords = get_fil_coords(fil_type, fil_initial_coords[l][0], fil_initial_coords[l][1], domain_width, dz_i);
      fil_node_coords[l] = temp_coords;

    }

    for(int l=num_filaments; l<(num_filaments+num_braid_fils); l++){
      std::cout << "braid filament " << l << "...\n";

      std::vector<std::vector<double>> temp_coords = get_fil_coords2(fil_type, fil_initial_coords[0][0]+0.002, fil_initial_coords[0][1]+0.01+fil_initial_coords[l][0], fil_initial_coords[l][1], fil_initial_coords[0][0]+0.033, fil_initial_coords[0][1]-0.01+fil_initial_coords[l][0], fil_initial_coords[l][1], dz_i);
      fil_node_coords[l] = temp_coords;
    }
  }
  
  int timestep_start = 0;
  //int num_billows = 7; //number of billows to initialize and have in domain at all times - should alter to make generically compatible with initialization from file
  int num_braid_pairs = std::trunc(domain_width/braid_spacing2);  //number of braid pairs within domain, dictated by domain width and pairing spacing
  
  //get z-coordinate correct for center of first vortex pair:
  double z0 = (domain_width/2-(num_braid_pairs-1)*braid_spacing2/2);  //z-coordinate of first vortex pair
  for(int m=num_filaments; m<(num_filaments+num_braid_fils); m++){
      //adjust z-coordinate to of each filament to be at correct location for first of vortex pair
      for(int n=0; n<fil_node_coords[m].size(); n++){
        fil_node_coords[m][n][2] = fil_node_coords[m][n][2]+z0;
      }
  }

  //create vector of Filaments for storage once created
  std::vector<std::unique_ptr<Filament>> filaments;
  std::cout << "num_billows: " << num_billows << "\nnum_filaments: " << num_filaments << "\nnum_braid_pairs: " << num_braid_pairs << "\nnum_braid_fils: " << num_braid_fils << "\n= " << num_billows*num_filaments + num_braid_pairs*2*num_braid_fils*num_billows << " filaments.\n";
  filaments.reserve(num_billows*num_filaments + num_braid_pairs*2*num_braid_fils*num_billows);   //2 vortices per braid pair, num_braid_fil filaments per vortex

  std::ofstream log_file;
  log_file.open("log_file.txt", std::ios_base::app);

  if(choice==2){
    std::cout << "Reading velocity field from file...(" + full_res_path + ")\n";
    read_vtk(full_res_path, filaments, fil_circulation, core_radius, max_sl, min_sl, max_range, background_vel);
    std::cout << "Back from read_vtk...\n";
    output_velocity_field("test_vel_field_40.vtu", filaments, background_vel, 0.0, 0.1, -0.02, 0.02, 0.01, 0.03, 40, 30, 40);
    return 0;
  }

  if(initialize_fr_file){
    //if initialize from file, initialize filaments vector from previous data (STILL NEED TO GET TIMESTEP)
    read_vtk(initialize_file_path, filaments, fil_circulation, core_radius, max_sl, min_sl, max_range, background_vel);
    //read timestep from end of file name:
    std::string timestep_str = initialize_file_path.substr(initialize_file_path.find_last_of("_")+1, initialize_file_path.find_last_of(".")-initialize_file_path.find_last_of("_")-1);
    int timestep_start = std::stoi(timestep_str);
  }else{
    //initialize with multiple filaments (then will replace those at end of domain with new ones to develop naturually-ish)
    for(int p=0; p<num_billows; p++){
      //loop through filaments and shift x coordinates by 0.035*p (avg. spacing between billows)
      //except then changed to mimic real billow locations at a random timestep from viscous simulation (hence why dif. #s for each)
      if(p==1){
        for(int n=0; n<(num_filaments+num_braid_fils); n++){
          for(int o=0; o<fil_node_coords[n].size(); o++){
            fil_node_coords[n][o][0] += 0.0373;
            fil_node_coords[n][o][1] += 0.0003;
          }
        }
      }else if(p==2){
        for(int n=0; n<(num_filaments+num_braid_fils); n++){
          for(int o=0; o<fil_node_coords[n].size(); o++){
            fil_node_coords[n][o][0] += 0.034;
            fil_node_coords[n][o][1] += 0.00015;
          }
        }
      }else if(p==3){
        for(int n=0; n<(num_filaments+num_braid_fils); n++){
          for(int o=0; o<fil_node_coords[n].size(); o++){
            fil_node_coords[n][o][0] += 0.037;
            fil_node_coords[n][o][1] -= 0.00085;
          }
        }
      }else if(p==4){
        for(int n=0; n<(num_filaments+num_braid_fils); n++){
          for(int o=0; o<fil_node_coords[n].size(); o++){
            fil_node_coords[n][o][0] += 0.02938;
            fil_node_coords[n][o][1] += 0.0052;
          }
        }
      }else if(p==5){
        for(int n=0; n<(num_filaments+num_braid_fils); n++){
          for(int o=0; o<fil_node_coords[n].size(); o++){
            fil_node_coords[n][o][0] += 0.04165;
            fil_node_coords[n][o][1] += 0.002;
          }
        }
      }else if(p==6){
        for(int n=0; n<(num_filaments+num_braid_fils); n++){
          for(int o=0; o<fil_node_coords[n].size(); o++){
            fil_node_coords[n][o][0] += 0.016;
            fil_node_coords[n][o][1] -= 0.0061;
          }
        }
      }

      //KH FILAMENTS:

      //fil_circ based on KH_circulation and number of filaments
      fil_circulation = KH_circulation/num_filaments;

      //create filament(s) and add to vector
      for(int m=0; m<num_filaments; m++){
        auto temp_f = std::make_unique<Filament>(fil_circulation, core_radius, max_sl, min_sl, max_range, fil_node_coords[m], background_vel, fil_type);
        int fil_i = p*num_filaments+m; //can use as index if push_back is the issue   ***  NO LONGER ACCURATE WITH BRAIDS  ***
        filaments.push_back(std::move(temp_f));
        int current_index = filaments.size()-1;
        log_file << "Filament " << fil_i << " added to filaments at index " << current_index << ".\n";
      }

      //BRAID VORTICES:
      //***OG LOCATION OF z0***

      //create braid vortices and add to vector                     ****  For cleanliness could make a function to do this elsewhere  ****
      for(int i=0; i<num_braid_pairs; i++){
        std::cout << "Braid vortex " << i << "...\n";
        std::cout << "Z: " << fil_node_coords[num_filaments][0][2] << ".\n";

        //if z coord is outside domain_width, error return:
        if(fil_node_coords[num_filaments][0][2] > domain_width){
          std::cout << "Error: z-coordinate of braid vortex outside domain width.\n";
          return 1;
        }

        for(int m=num_filaments; m<(num_filaments+num_braid_fils); m++){
          //adjust z-coordinate to of each filament to be at correct location for first of vortex pair
          for(int n=0; n<fil_node_coords[m].size(); n++){
            fil_node_coords[m][n][2] = fil_node_coords[m][n][2]-braid_spacing1/2;
          }
          //first vortex of pair:
          auto temp_f = std::make_unique<Filament>((-1)*braid_circ, core_radius, max_sl, min_sl, max_range, fil_node_coords[m], background_vel, fil_type);
          int fil_i = p*num_filaments+m; //can use as index if push_back is the issue  ***  NO LONGER ACCURATE WITH BRAIDS  ***
          filaments.push_back(std::move(temp_f));
          int current_index = filaments.size()-1;
          log_file << "Braid filament " << fil_i << " (1/2) added to filaments at index " << current_index << ".\n";

          //adjust z-coordinate to of each filament to be at correct location for second of vortex pair
          //***  doing this here means order of filaments will alternate between pairs for when there is more than one fil/vortex  ***
          //***  could do outside here but it would require a second loop through filaments  ***
          for(int n=0; n<fil_node_coords[m].size(); n++){
            fil_node_coords[m][n][2] = fil_node_coords[m][n][2]+braid_spacing1;
          }

          //second vortex of pair:
          auto temp_f2 = std::make_unique<Filament>(braid_circ, core_radius, max_sl, min_sl, max_range, fil_node_coords[m], background_vel, fil_type);
          fil_i = p*num_filaments+m; //can use as index if push_back is the issue  ***  NO LONGER ACCURATE WITH BRAIDS  ***
          filaments.push_back(std::move(temp_f2));
          current_index = filaments.size()-1;
          log_file << "Braid filament " << fil_i << " (2/2) added to filaments at index " << current_index << ".\n";

          //return z-coordinate to of each filament to center of pair, and increment to next position:
          for(int n=0; n<fil_node_coords[m].size(); n++){
            fil_node_coords[m][n][2] = fil_node_coords[m][n][2]-braid_spacing1/2+braid_spacing2;
          }

        }

      }

      //return braid vortices center to z-coordinate of first vortex pair for after next billow:
      for(int m=num_filaments; m<(num_filaments+num_braid_fils); m++){
        for(int n=0; n<fil_node_coords[m].size(); n++){
          fil_node_coords[m][n][2] = fil_node_coords[m][n][2]-braid_spacing2*num_braid_pairs;
        }
      }

    }

    //return filament coords to original positions (for transient billow introduction)
    for(int n=0; n<(num_filaments+num_braid_fils); n++){
      for(int o=0; o<fil_node_coords[n].size(); o++){
        //fil_node_coords[n][o][0] -= (0.04*(num_billows-1));
        fil_node_coords[n][o][0] -= (0.0373+0.034+0.037+0.02938+0.04165+0.016);
        fil_node_coords[n][o][1] -= (0.0003+0.00015-0.00085+0.0052+0.002-0.0061);
      }
    }
  }
  

  log_file.close();

  //counter variable to keep track of billow to be replaced
  int last_billow_i = num_billows-1;
  /*
  in order to properly keep track of everything, all filaments belonging to KH billows will be first in
  'filaments'.  Then the braid filaments will be added to the end of the vector.
  */

  //timestep_start (if initialize from file, 0 otherwise)
  for (int j = timestep_start; j < num_timesteps; j++) {
    //convert num_billows & num_filaments to string for filename:
    std::string fn_template = std::to_string(num_filaments) + "f_" + std::to_string(num_billows) + "kh_" + std::to_string(num_braid_fils) + "b_lesscirc_timestep_" + std::to_string(j);
    
    if(j != 0 && j%KH_num_steps == 0){
      int fil_start_i = last_billow_i*(num_filaments+num_braid_fils*num_braid_pairs*2);
      for(int m=0; m<num_filaments; m++){
            //    -----------------  need to figure out last billow for initialized result  ----------------- 
        //making new KH filament to replace old one: 
        std::cout << "Making temp_f...\n";
        auto temp_f = std::make_unique<Filament>(fil_circulation, core_radius, max_sl, min_sl, max_range, fil_node_coords[m], background_vel, fil_type);
        std::cout << "Made temp_f, moving to filaments... ";
        filaments[fil_start_i+m] = std::move(temp_f);
        std::cout << "Moved successfully\n";

        log_file.open("log_file.txt", std::ios_base::app);

        log_file << "Filament at index " << last_billow_i+m << " replaced at timestep " << j << "...\n";

        log_file.close();
      } 

      //same thing for braid filaments:
      //create braid vortices and add to vector                     ****  For cleanliness could make a function to do this elsewhere  ****
      for(int i=0; i<num_braid_pairs; i++){

        //if z coord is outside domain_width, error return:
        if(fil_node_coords[num_filaments][0][2] > domain_width){
          std::cout << "Error: z-coordinate of braid vortex outside domain width.\n";
          return 1;
        }

        for(int m=num_filaments; m<(num_filaments+num_braid_fils); m++){
          //adjust z-coordinate to of each filament to be at correct location for first of vortex pair
          for(int n=0; n<fil_node_coords[m].size(); n++){
            fil_node_coords[m][n][2] = fil_node_coords[m][n][2]-braid_spacing1/2;
          }
          //first vortex of pair:
          auto temp_f = std::make_unique<Filament>((-1)*braid_circ, core_radius, max_sl, min_sl, max_range, fil_node_coords[m], background_vel, fil_type);
          filaments[fil_start_i + (2*m - num_filaments + i*num_braid_fils*2)] = std::move(temp_f);

          //adjust z-coordinate to of each filament to be at correct location for second of vortex pair
          //***  doing this here means order of filaments will alternate between pairs for when there is more than one fil/vortex  ***
          //***  could do outside here but it would require a second loop through filaments  ***
          for(int n=0; n<fil_node_coords[m].size(); n++){
            fil_node_coords[m][n][2] = fil_node_coords[m][n][2]+braid_spacing1;
          }

          //second vortex of pair:
          auto temp_f2 = std::make_unique<Filament>(braid_circ, core_radius, max_sl, min_sl, max_range, fil_node_coords[m], background_vel, fil_type);
          filaments[fil_start_i + (2*m - num_filaments + i*num_braid_fils*2) + 1] = std::move(temp_f2);

          //return z-coordinate to of each filament to center of pair, and increment to next position:
          for(int n=0; n<fil_node_coords[m].size(); n++){
            fil_node_coords[m][n][2] = fil_node_coords[m][n][2]-braid_spacing1/2+braid_spacing2;
          }
        }
      }

      //return braid vortices center to z-coordinate of first vortex pair for after next billow:
      for(int m=num_filaments; m<(num_filaments+num_braid_fils); m++){
        for(int n=0; n<fil_node_coords[m].size(); n++){
          fil_node_coords[m][n][2] = fil_node_coords[m][n][2]-braid_spacing2*num_braid_pairs;
        }
      }


      //billow index for next time:
      if(last_billow_i == 0){
        last_billow_i = num_billows-1;
      }else{
        last_billow_i--;
      }
    }
    
    //output every fourth timestep, step each one
    if(j%5 == 0){
      save_multiple_fil_timestep(fn_template, filaments);
    }

    //if(j%10 ==0){
    //  output_velocity_field(fn_template, filaments, background_vel, 0.0, 0.15, -0.02, 0.02, 0.01, 0.03, 120, 60, 40);
    //}

    step(filaments, dt);

  }
  
  return 0;
  
}