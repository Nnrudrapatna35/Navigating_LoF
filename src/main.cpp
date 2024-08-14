/** ------ Libraries ---------------------------------------------------------*/
#include <chrono>
#include <iostream>
#include <vector>
#include <cmath>
#include "boost/multi_array.hpp"

/** ------ Project-specific header files -------------------------------------*/
#include "GlobalConfigurations.hpp"
#include "SimpleFunctions.hpp"
#include "MemoryAllocations.hpp"
#include "CTimeDependentHjbSolver.hpp"
#include "CTimeGrid.hpp"
#include "CFoodDensity.hpp"
#include "WritetoFile.hpp"
#include "CTimeDependentTracer.hpp"
#include "CSVReader.hpp"
#include "CTerrain.hpp"

/* Example list
 * ExampleRiskReward - two patch example with a single stage and variable risk
 *                     premium
 * ExampleRealWorld - environment inspired by data for Samango monkeys
 * ExampleMultiStage - two patch environment with food depletion and multiple stages
 */

/* Impact of food depletion on patch selection */
std::shared_ptr<CTimeDependentHjbSolver> ExampleMultiStage(const double aSpottingRate,
                                                  const double aMaxKillRate,
                                                  const bool aAllowModeSwitches, 
                                                  const objective_t aObjective,
                                                  const int aNPhysical,
                                                  const int aNe,
                                                  const int aNStages,
                                                  const double aT) {
  // Grid parameters
  const int n_physical = aNPhysical;
  const int n_energy = aNe;
  const double physical_min = 0;
  const double physical_max = 0.6;
  const double energy_min = 0;
  const double energy_max = 3;
  const double time_max = aT;
  const int n_stages = aNStages;
  const bool full_grid = true;

  // Food parameters
  const double food_radius = 0.05;
  const double harvest_rate = 0.01;
  const bool food_depletion = true;

  // Path parameters
  const std::vector<double> source_point {0.1, 0.3};
  const double initial_energy = 1;
  const bool return_home = false;

  // Environmental parameters
  const double running_cost_const = 0.08;
  const double speed_const = 5;

  // Set function parameters for both modes
  std::function<double(double, double, double, double)> cost_mode_1 =
    std::bind(varConstant4D, std::placeholders::_1, std::placeholders::_2,
              std::placeholders::_3, std::placeholders::_4, running_cost_const);
  std::function<double(double, double, double, double)> cost_mode_2 =
    std::bind(varConstant4D, std::placeholders::_1, std::placeholders::_2,
              std::placeholders::_3, std::placeholders::_4, running_cost_const);

  std::function<double(double, double)> food_density_fnc = &twoPeaksClose;

  std::function<double(double, double, double, double)> speed_mode_1 =
    std::bind(varConstant4D, std::placeholders::_1, std::placeholders::_2,
              std::placeholders::_3, std::placeholders::_4, speed_const);
  std::function<double(double, double, double, double)> speed_mode_2 =
    std::bind(varConstant4D, std::placeholders::_1, std::placeholders::_2,
              std::placeholders::_3, std::placeholders::_4, speed_const);

  std::function<double(double, double, double, double)> spotting_rate =
    std::bind(varConstant4D, std::placeholders::_1, std::placeholders::_2,
              std::placeholders::_3, std::placeholders::_4, aSpottingRate);
  std::function<double(double, double, double, double)> kill_rate =
    std::bind(&rectangle_kill, std::placeholders::_1, std::placeholders::_2,
              std::placeholders::_3, std::placeholders::_4, aMaxKillRate);
  std::function<double(double, double, double, double)> give_up_rate =
    &rectangle_give_up;

  std::function<double(double, double)> food_kernel = &constant_kernel;
  
  std::function<double(double, double)> home_base = &square_lower_left;
  std::function<double(double, double)> obstacle = &none;

  std::shared_ptr<CTerrain> terrain = std::make_shared<CTerrain>(
    CTerrain(spotting_rate, kill_rate, give_up_rate, home_base, obstacle,
             n_physical, physical_min, physical_max));

  double theta;
  double energy_threshold;
  objective_t solver_objective = aObjective;
  std::function<double(double)> utility_fnc;
  if (aObjective == SQUAREROOT) {
    theta = 0.5;
    utility_fnc = std::bind(power_utility, std::placeholders::_1, theta, energy_max);
    solver_objective = UTILITY;
  } else if (aObjective == SIGMOID) {
    theta = 0.1*energy_max;
    energy_threshold = 0.5*energy_max;
    utility_fnc = std::bind(sigmoid_utility_normcdf, std::placeholders::_1,
                            theta, energy_threshold, energy_max);
    solver_objective = UTILITY;
  } else if (aObjective == QUADRATIC) {
    theta = 2.0;
    utility_fnc = std::bind(power_utility, std::placeholders::_1, theta, energy_max);
    solver_objective = UTILITY;
  } else if (aObjective == LINEAR) {
    theta = 1.0;
    utility_fnc = std::bind(power_utility, std::placeholders::_1, theta, energy_max);
    solver_objective = UTILITY;
  } else if (aObjective == UTILITY) {
    theta = 1.0;
    utility_fnc = std::bind(power_utility, std::placeholders::_1, theta, energy_max);
  }

  std::shared_ptr<CFoodDensity> food_density = std::make_shared<CFoodDensity>(
  CFoodDensity(food_density_fnc, food_kernel, utility_fnc, terrain,
                food_radius, harvest_rate, food_depletion,
                n_physical, physical_min, physical_max));
  
  // Build solver
  CTimeDependentHjbSolver solver = CTimeDependentHjbSolver(cost_mode_1, cost_mode_2,
                                    speed_mode_1, speed_mode_2, food_density,
                                    terrain, energy_threshold, n_physical,
                                    n_energy, physical_min, physical_max,
                                    energy_min, energy_max, time_max,
                                    source_point, initial_energy,
                                    solver_objective, n_stages, full_grid, 
                                    return_home, aAllowModeSwitches);
  // return a pointer
  return std::make_shared<CTimeDependentHjbSolver>(solver);
}

/* ExampleRiskReward: Trade off between food available and risk of predation */
std::shared_ptr<CTimeDependentHjbSolver> ExampleRiskReward(const double aSpeed,
                                                   const bool aFoodDepletion,
                                                   const objective_t aObjective,
                                                   const int aNPhysical,
                                                   const int aNe,
                                                   const int aNStages,
                                                   const bool aIsBaselineShift,
                                                   const double aRateShift,
                                                   const bool aAllowModeSwitches, 
                                                   const bool aFullGrid) {
  // Set grid parameters
  const double physical_min = 0;
  const double physical_max = 1;
  const double energy_min = 0;
  const double energy_max = 10.0;
  const double time_max = 1;
  const double food_radius = 0.01;
  const double harvest_rate = 0.001;
  const bool food_depletion = aFoodDepletion;
  const bool return_home = true;

  // Set tracer parameters
  const std::vector<double> source_point {0.41,0.5};
  const double initial_energy = 0.7 * energy_max;

  const double running_cost_const = 3.0;
  const double speed_const = aSpeed;
  const double kill_rate_const = 4.0;
  const double give_up_rate_const = 1.0;

  // Function-based parameters for both modes
  std::function<double(double, double, double, double)> cost_mode_1 =
      std::bind(varConstant4D, std::placeholders::_1, std::placeholders::_2,
                std::placeholders::_3, std::placeholders::_4, running_cost_const);
  std::function<double(double, double, double, double)> cost_mode_2 =
      std::bind(varConstant4D, std::placeholders::_1, std::placeholders::_2,
                std::placeholders::_3, std::placeholders::_4, running_cost_const);

  std::function<double(double, double, double, double)> spotting_rate;
  std::function<double(double, double)> food_density_fnc;
  if (aIsBaselineShift == false) {
    food_density_fnc = &circleFood3Cts;
    spotting_rate = std::bind(varCirclePred3Cts, std::placeholders::_1, 
                              std::placeholders::_2, std::placeholders::_3, 
                              std::placeholders::_4, aRateShift);
  } else {
    food_density_fnc = &circleFood4Cts;
    spotting_rate = std::bind(varCirclePred4Cts, std::placeholders::_1, 
                              std::placeholders::_2, std::placeholders::_3, 
                              std::placeholders::_4, aRateShift);
  }

  std::function<double(double, double, double, double)> speed_mode_1 =
      std::bind(varConstant4D, std::placeholders::_1, std::placeholders::_2,
                std::placeholders::_3, std::placeholders::_4, speed_const);
  std::function<double(double, double, double, double)> speed_mode_2 =
      std::bind(varConstant4D, std::placeholders::_1, std::placeholders::_2,
                std::placeholders::_3, std::placeholders::_4, speed_const);

  std::function<double(double, double, double, double)> give_up_rate =
      std::bind(varConstant4D, std::placeholders::_1, std::placeholders::_2,
                std::placeholders::_3, std::placeholders::_4,
                give_up_rate_const);
  std::function<double(double, double, double, double)> kill_rate =
      std::bind(varConstant4D, std::placeholders::_1, std::placeholders::_2,
                std::placeholders::_3, std::placeholders::_4, kill_rate_const);

  std::function<double(double, double)> food_kernel = &constant_kernel;

  double theta;
  double energy_threshold;
  objective_t solver_objective = aObjective;
  std::function<double(double)> utility_fnc;
  if (aObjective == SQUAREROOT) {
    theta = 0.5;
    utility_fnc = std::bind(power_utility, std::placeholders::_1, theta, energy_max);
    solver_objective = UTILITY;
  } else if (aObjective == SIGMOID) {
    theta = 1.0;
    energy_threshold = 3.0;
    utility_fnc = std::bind(sigmoid_utility_normcdf,std::placeholders::_1,
                            theta, energy_threshold,energy_max);
    solver_objective = UTILITY;
  } else if (aObjective == QUADRATIC) {
    theta = 2.0;
    utility_fnc = std::bind(power_utility, std::placeholders::_1, theta, energy_max);
    solver_objective = UTILITY;
  } else if (aObjective == LINEAR) {
    theta = 1.0;
    utility_fnc = std::bind(power_utility, std::placeholders::_1, theta, energy_max);
    solver_objective = UTILITY;
  } else if (aObjective == UTILITY) {
    theta = 2.0;
    energy_threshold = 3.0;
    utility_fnc = std::bind(sigmoid_utility_normcdf,std::placeholders::_1,
                            theta, energy_threshold,energy_max);
  }

  std::function<double(double, double)> home_base = &risk_reward_hb;
  std::function<double(double, double)> obstacle = &none;

  std::shared_ptr<CTerrain> terrain = std::make_shared<CTerrain>(
    CTerrain(spotting_rate, kill_rate, give_up_rate, home_base, obstacle,
             aNPhysical, physical_min, physical_max));

  std::shared_ptr<CFoodDensity> food_density = std::make_shared<CFoodDensity>(
    CFoodDensity(food_density_fnc, food_kernel, utility_fnc, terrain,
                 food_radius, harvest_rate, food_depletion,
                 aNPhysical, physical_min, physical_max));

  // Build solver
  CTimeDependentHjbSolver solver = CTimeDependentHjbSolver(cost_mode_1, cost_mode_2,
                                    speed_mode_1, speed_mode_2, food_density,
                                    terrain, energy_threshold, aNPhysical,
                                    aNe, physical_min, physical_max,
                                    energy_min, energy_max, time_max,
                                    source_point, initial_energy,
                                    solver_objective, aNStages, aFullGrid, 
                                    return_home, aAllowModeSwitches);

  // Return a pointer to a solver
  return std::make_shared<CTimeDependentHjbSolver>(solver);
}

/* ExampleRealWorld: An example inspired by real world data */
std::shared_ptr<CTimeDependentHjbSolver> ExampleRealWorld(const bool aFoodDepletion,
                                                  const bool aAllowModeSwitches,
                                                  const objective_t aObjective,
                                                  const int aNPhysical,
                                                  const int aNe,
                                                  const int aNStages) {
  // Set grid parameters
  const int n_physical = aNPhysical;
  const int n_energy = aNe;
  const double physical_min = 0;
  const double physical_max = 1;
  const double energy_min = 0;
  const double energy_max = 1400;
  const double time_max = 6;
  const int n_stages = aNStages;
  const double food_radius = 0.05; 
  const double harvest_rate = 1e-6;

  // Extra parameters
  const std::vector<double> source_point {0.2,0.7};
  const double initial_energy = energy_max/2;

  const bool full_grid = true;
  const bool return_home = false;

  // Set function parameters for both modes
  std::shared_ptr<memory::array2D_t<double>> cost_mode_1 = std::make_shared<memory::array2D_t<double>>(memory::allocateArray2D<double>(n_physical, n_physical));
  std::shared_ptr<memory::array2D_t<double>> cost_mode_2 = std::make_shared<memory::array2D_t<double>>(memory::allocateArray2D<double>(n_physical, n_physical));
  std::shared_ptr<memory::array2D_t<double>> food_density_array = read_csv("visualization/imread/Psi_201.csv");
  std::shared_ptr<memory::array2D_t<double>> spotting_rate = read_csv("visualization/imread/MuS_201.csv");
  std::shared_ptr<memory::array2D_t<double>> speed_mode_1 = std::make_shared<memory::array2D_t<double>>(memory::allocateArray2D<double>(n_physical, n_physical));
  std::shared_ptr<memory::array2D_t<double>> speed_mode_2 = std::make_shared<memory::array2D_t<double>>(memory::allocateArray2D<double>(n_physical, n_physical));
  std::shared_ptr<memory::array2D_t<double>> give_up_rate = read_csv("visualization/imread/MuG_201.csv");
  std::shared_ptr<memory::array2D_t<double>> kill_rate = std::make_shared<memory::array2D_t<double>>(memory::allocateArray2D<double>(n_physical, n_physical));
  std::function<double(double, double)> food_kernel = &constant_kernel;
  std::shared_ptr<memory::array2D_t<double>> home_base = std::make_shared<memory::array2D_t<double>>(memory::allocateArray2D<double>(n_physical, n_physical));
  std::shared_ptr<memory::array2D_t<double>> obstacle = read_csv("visualization/imread/Obstacle_201.csv");
  std::shared_ptr<memory::array2D_t<double>> sleeping_sites = read_csv("visualization/imread/SleepingSites_201.csv");

  // Set environment parameters
  const double running_cost_const = 14.0;
  const double speed_const = 0.7;
  const double kill_rate_max = 5.0;
  const double spotting_rate_max = 0.12;
  const double spotting_rate_baseline = 0.01;
  const double food_density_max = 2.2e8;
  const double food_density_baseline = 2e6;
  const double give_up_rate_max = 10.0;
  const double give_up_rate_baseline = 1.0;

  for(int i = 0; i < n_physical; ++i){
    for(int j = 0; j < n_physical; ++j){
      double x = physical_min + (double) i / (n_physical-1) * (physical_max-physical_min);
      double y = physical_min + (double) j / (n_physical-1) * (physical_max-physical_min);

      (*cost_mode_1)[i][j] = running_cost_const;
      (*cost_mode_2)[i][j] = running_cost_const;
      (*food_density_array)[i][j] = food_density_max*(*food_density_array)[i][j]
                                    + food_density_baseline; 
      (*spotting_rate)[i][j] = spotting_rate_max*(*spotting_rate)[i][j]
                               + spotting_rate_baseline;
      (*speed_mode_1)[i][j] = speed_const;
      (*speed_mode_2)[i][j] = speed_const;
      (*give_up_rate)[i][j] = give_up_rate_max*(*give_up_rate)[i][j] 
                              + give_up_rate_baseline;
      (*kill_rate)[i][j] = kill_rate_max;
      (*home_base)[i][j] = (&none)(x,y);
    }
  }

  double theta;
  double energy_threshold;
  objective_t solver_objective = aObjective;
  std::function<double(double)> utility_fnc;
  if (aObjective == SQUAREROOT) {
    theta = 0.5;
    utility_fnc = std::bind(power_utility, std::placeholders::_1, theta, energy_max);
    solver_objective = UTILITY;
  } else if (aObjective == SIGMOID) {
    theta = 0.1*energy_max;
    energy_threshold = 0.1*energy_max;
    utility_fnc = std::bind(sigmoid_utility_normcdf, std::placeholders::_1,
                            theta, energy_threshold, energy_max);
    solver_objective = UTILITY;
  } else if (aObjective == QUADRATIC) {
    theta = 2.0;
    utility_fnc = std::bind(power_utility, std::placeholders::_1, theta, energy_max);
    solver_objective = UTILITY;
  } else if (aObjective == LINEAR) {
    theta = 1.0;
    utility_fnc = std::bind(power_utility, std::placeholders::_1, theta, energy_max);
    solver_objective = UTILITY;
  } else if (aObjective == UTILITY) {
    theta = 1.0;
    utility_fnc = std::bind(power_utility, std::placeholders::_1, theta, energy_max);
  }

  std::shared_ptr<CTerrain> terrain = std::make_shared<CTerrain>(
    CTerrain(spotting_rate, kill_rate, give_up_rate, home_base, obstacle,
             n_physical, physical_min, physical_max, sleeping_sites));

  std::shared_ptr<CFoodDensity> food_density = std::make_shared<CFoodDensity>(
    CFoodDensity(food_density_array, food_kernel, utility_fnc, terrain,
                 food_radius, harvest_rate, aFoodDepletion,
                 n_physical, physical_min, physical_max));

  // Build solver
  CTimeDependentHjbSolver solver = CTimeDependentHjbSolver(cost_mode_1, cost_mode_2,
                                    speed_mode_1, speed_mode_2, food_density,
                                    terrain, energy_threshold, n_physical,
                                    n_energy, physical_min, physical_max,
                                    energy_min, energy_max, time_max,
                                    source_point, initial_energy, 
                                    solver_objective, n_stages, full_grid, 
                                    return_home, aAllowModeSwitches);
  // Return a pointer
  return std::make_shared<CTimeDependentHjbSolver>(solver);
}


/** ------ Namespaces --------------------------------------------------------*/
int main(int argc, char* argv[]) {
  /** Select example based on command line argument */
  if (argc < 3) {
    std::cout << "Insufficient arguments given, terminating" << std::endl;
    assert(false);
  }

  std::string arg1 = std::string(argv[1]);

  std::string arg2 = std::string(argv[2]);
  objective_t objective;
  std::string objective_name;
  if (arg2 == "E") {
    objective = ENERGY;
    objective_name = "_Energy";
  } else if (arg2 == "S") {
    objective = SURVIVAL;
    objective_name = "_Survival";
  } else if (arg2 == "U") {
    objective = UTILITY;
    objective_name = "_Utility";
  } else if (arg2 == "M") {
    objective = SIGMOID;
    objective_name = "_Sigmoid";
  } else if (arg2 == "R") {
    objective = SQUAREROOT;
    objective_name = "_SquareRoot";
  } else if (arg2 == "Q") {
    objective = QUADRATIC;
    objective_name = "_Quadratic";
  } else if (arg2 == "L") {
    objective = LINEAR;
    objective_name = "_Linear";
  } else {
    std::cout << "Invalid objective entered" << std::endl;
    assert(false);
  }

  // Default stencil is 5 pt
  stencil_t stencil = FIVE_POINT;
  std::string stencil_name = "";
  bool path_tracing = false;
  if (argc == 4) {
    const std::string arg3 = std::string(argv[3]);
    std::cout << arg3 << std::endl;
    if (arg3 == "9") {
      stencil_name = "_9pt";
      stencil = NINE_POINT;
    } else if (arg3 == "5") {
      stencil = FIVE_POINT;
    } else if (arg3 == "P"){
      path_tracing = true;
    } else {
      std::cout << "Invalid stencil type entered" << std::endl;
      assert(false);
    }
  }

   if (argc == 5) {
    const std::string arg3 = std::string(argv[3]);
    const std::string arg4 = std::string(argv[4]);
    if (arg3 == "9") {
      stencil_name = "_9pt";
      stencil = NINE_POINT;
    } else if (arg3 == "5") {
      stencil = FIVE_POINT;
    } else {
      std::cout << "Invalid argument entered" << std::endl;
      assert(false);
    }
    if (arg4 == "P"){
      path_tracing = true;
    } else {
      std::cout << "Invalid argument entered" << std::endl;
      assert(false);
    }
  }

  if (argc > 5) {
    std::cout << "Too many arguments given, terminating" << std::endl;
    assert(false);
  }

  std::shared_ptr<CTimeDependentHjbSolver> timedepsolver;
  std::string filename = arg1 + objective_name + stencil_name;

  /* Decide if we enable predation in the path tracer */
  bool allow_predation = true;
  std::string params = "";

  /** Start timer */
  auto t1 = std::chrono::high_resolution_clock::now();

  if (arg1 == "ExampleRealWorld") {
    const int n_physical = 201;
    const int n_energy = 71;
    const int n_stages = 1;
    const bool food_depletion = true;

    timedepsolver = ExampleRealWorld(food_depletion, allow_predation, objective, 
                                     n_physical, n_energy, n_stages);
  }
  else if (arg1 == "ExampleMultiStage") {
    const int n_physical = 61;//121;
    const int n_energy = 101;//301;
    const double spotting_rate = 3;
    const double kill_rate = 3;
    int N = 5;
    double T = 1;

    allow_predation = false;
    for (int realization = 0; realization < 6; ++realization) {
      timedepsolver = ExampleMultiStage(spotting_rate, kill_rate, allow_predation, 
                                objective, n_physical, n_energy, N, T);
      params = "_N_" + std::to_string(N) + "_mus_" + std::to_string((int)spotting_rate) 
              + "_muk_" + std::to_string((int)kill_rate)
              + "_realization_" + std::to_string(realization);
      filename = arg1 + objective_name + stencil_name + params;

      timedepsolver->solveMultiStage(filename, stencil);

      allow_predation = true;
    }
    
    return 0;
  }
  else if (arg1 == "ExampleRiskReward") {
    /* This conditional creates the paths necessary for the threshold plot*/
    const int n_physical = 101;
    const int n_energy = 101;
    const int n_stages = 1;
    const int n_initial_energies = 20;

    const double max_energy = 10.0;
    const double speed = 2.0;

    allow_predation = false;
    const bool is_final_stage = true;
    const bool classify_path = true;
    const bool use_baseline_shift = false;
    const bool read_from_file = path_tracing;
    const bool read_from_file_inner = false;
    const bool full_grid = true;

    bool food_depletion;
    std::string include_food_depletion;
    std::cout << "Would you like to run this example with food depletion? (Y/N)" 
                 << std::endl;
    std::cin >> include_food_depletion;
    if (include_food_depletion == "Y") {
      food_depletion = true;
      include_food_depletion = "";
    } else if (include_food_depletion == "N") {
      food_depletion = false;
      include_food_depletion = "_NoDepletion";
    } else {
      std::cout << "Invalid response. Terminating." << std::endl;
      assert(false);
    }

    for (int rate_scaled = 1; rate_scaled <= 30; rate_scaled += 1) {
      double rate_shift = rate_scaled / 10.0;
      timedepsolver = ExampleRiskReward(speed, food_depletion,objective, 
                                        n_physical, n_energy, n_stages, 
                                        use_baseline_shift, rate_shift, 
                                        allow_predation, full_grid); 
      std::string param = "";
      filename = arg1 + objective_name + stencil_name + param;\
      
      if (!read_from_file) {
        const int current_stage = 0;
        timedepsolver->solveHJB(current_stage, filename, stencil);
        timedepsolver->writeResults(filename, current_stage);
      } else {
        timedepsolver->readValuesFromFile(filename);
      }

      // Trace paths for a range of initial energies
      std::string path_folder = "RR_PatchSelection/";
      std::cout << "Tracing paths for " << n_initial_energies 
                << " initial energy values..." << std::endl;
      for (int k = 0; k < n_initial_energies; ++k){
        const double initial_energy = (k+1)*(max_energy/n_initial_energies);
        if (use_baseline_shift) {
          param = "_BaselineShift_" + std::to_string(rate_scaled) 
                  + "_e0_" + std::to_string(k+1) + include_food_depletion;
        } else {
          param = "_UpperRateShift_" + std::to_string(rate_scaled) 
                  + "_e0_" + std::to_string(k+1) + include_food_depletion;
        }

        filename = path_folder + arg1 + objective_name + stencil_name + param;

        timedepsolver->setInitialEnergy(initial_energy);
        timedepsolver->resetFoodDensity();
        timedepsolver->traceHJB(is_final_stage, filename, read_from_file_inner, 
                                classify_path);
      }
      std::cout << "Done." << std::endl;
    }

    return 0;
  }
  else if (arg1 == "ExampleRiskRewardPT") {
    /* This conditional creates the paths necessary for trajectory comparisons */
    const int n_physical = 201;
    const int n_energy = 151;
    const int n_stages = 1;
    const int n_initial_energies = 20;
    int n_paths;
    std::string path_folder;
    std::string realizations;

    std::cout << "Would you like to run many realizations? (Y/N)" 
              << std::endl;
    std::cin >> realizations;
    if (realizations == "Y") {
      n_paths = 100;
      allow_predation = true;
      path_folder = "RR_FoodDepletion/";
    } else if (realizations == "N") {
      n_paths = 1;
      allow_predation = false;
      path_folder = "RR_PathDeformation/";
    } else {
      std::cout << "Invalid response. Terminating." << std::endl;
      assert(false);
    }

    const double max_energy = 10.0;
    const double speed = 2.0;
    const int rate_scaled = 20;
    const double rate_shift = rate_scaled / 10.0; // Convert to the physical range

    const bool is_final_stage = true;
    const bool use_baseline_shift = false;
    const bool read_from_file = path_tracing;
    const bool read_from_file_inner = false;
    const bool full_grid = true;

    bool food_depletion;
    std::string include_food_depletion;
    std::cout << "Would you like to run this example with food depletion? (Y/N)" 
              << std::endl;
    std::cin >> include_food_depletion;
    if (include_food_depletion == "Y") {
      food_depletion = true;
      include_food_depletion = "";
    } else if (include_food_depletion == "N") {
      food_depletion = false;
      include_food_depletion = "_NoDepletion";
    } else {
      std::cout << "Invalid response. Terminating." << std::endl;
      assert(false);
    }

    timedepsolver = ExampleRiskReward(speed, food_depletion, objective, 
                                      n_physical, n_energy, n_stages, 
                                      use_baseline_shift, rate_shift, 
                                      allow_predation, full_grid); 
    
    std::string param = "_RateScaled_" + std::to_string(rate_scaled);

    filename = "ExampleRiskReward" + objective_name+ stencil_name + param;
    if (!read_from_file) {
      std::cout << "Solving for value function ..." << std::endl;
      const int current_stage = 0;
      timedepsolver->solveHJB(current_stage, filename, stencil);
      timedepsolver->writeResults(filename, current_stage);
      std::cout << " Done." << std::endl;
    } else {
      std::cout << "Reading value function..." << std::endl;
      timedepsolver->readValuesFromFile(filename);
      std::cout << " Done." << std::endl;
    }

    /* Trace paths for a range of initial energies */
    std::cout << "Tracing paths for " << n_initial_energies 
              << " initial energy values..." << std::endl;
    for (int k = 0; k < n_initial_energies; ++k){
      const double initial_energy = (k+1)*(max_energy/n_initial_energies);
      if (use_baseline_shift) {
        param = "_BaselineShift_" + std::to_string(rate_scaled) 
                + "_e0_" + std::to_string(k+1) + include_food_depletion;
      } else {
        param = "_UpperRateShift_" + std::to_string(rate_scaled) 
                + "_e0_" + std::to_string(k+1) + include_food_depletion;
      }

      filename = path_folder + "ExampleRiskReward" + objective_name + stencil_name + param;

      if (food_depletion) {
        timedepsolver->resetFoodDensity();
      }
      timedepsolver->setInitialEnergy(initial_energy);
      timedepsolver->traceNPaths(n_paths, filename, read_from_file_inner);
    }
    std::cout << "Done." << std::endl;

    return 0;
  }
  else {
    std::cout << "Invalid argument entered" << std::endl;
    assert(false);
  }

  if (path_tracing) {
    std::cout << "Reading in  value function ..." << std::endl;
    timedepsolver->readValuesFromFile(filename);
    std::cout << "Done." << std::endl;
    const bool read_from_file = false;

    while (path_tracing) {

      int initial_mode;
      std::cout << "Starting mode? (1/2)" << std::endl;
      std::cin >> initial_mode;
      if ((initial_mode != 1) && (initial_mode != 2)) {
        std::cout << "Invalid mode entered." << std::endl;
        assert(false);
      }
      timedepsolver->setInitialMode(initial_mode);
      std::string mode_label = "";
      if (initial_mode == 2) {
        mode_label = "_Mode2";
      }

      bool food_depletion;
      std::string include_food_depletion;
      std::cout << "Would you like to run this example with food depletion? (Y/N)" 
                << std::endl;
      std::cin >> include_food_depletion;
      if (include_food_depletion == "Y") {
        food_depletion = true;
        include_food_depletion = "";
      } else if (include_food_depletion == "N") {
        food_depletion = false;
        include_food_depletion = "_NoDepletion";
      } else {
        std::cout << "Invalid response. Terminating." << std::endl;
        assert(false);
      }

      std::string allow_mode_switches;
      std::cout << "Would you like to allow mode switches? (Y/N)" 
                << std::endl;
      std::cin >> allow_mode_switches;
      if (allow_mode_switches == "Y") {
        timedepsolver->setModeSwitches(true);
        allow_mode_switches = "";
      } else if (allow_mode_switches == "N") {
        timedepsolver->setModeSwitches(false);
        allow_mode_switches = "_NoModeSwitches";
      } else {
        std::cout << "Invalid response. Terminating." << std::endl;
        assert(false);
      }

      std::string specify_energy;
      std::cout << "Would you like to specify the starting energy? (Y/N)" 
                << std::endl;
      std::cin >> specify_energy;
      double initial_energy;
      if (specify_energy == "Y") {
        std::cout << "Enter starting energy:" 
                  << std::endl;
        std::cin >> initial_energy;
        timedepsolver->setInitialEnergy(initial_energy);
      } else if (specify_energy == "N") {
        initial_energy = timedepsolver->getInitialEnergy();
      } else {
        std::cout << "Invalid response. Terminating." << std::endl;
        assert(false);
      }

      std::string pathFolder = "RW_Paths/";
      std::string pathParams = ""; 
      std::string pathFilename = pathFolder + arg1 + objective_name + stencil_name
                               + mode_label+ include_food_depletion + allow_mode_switches;

      int num_paths;
      std::cout << "How many paths?" << std::endl;
      std::cin >> num_paths;
      if (num_paths == 1) {
        std::string specify_spotting_times;
        std::cout << "Would you like to specify spotting time (Y/N)?" << std::endl;
        std::cin >> specify_spotting_times;

        if (specify_spotting_times == "Y") {
          double spotting_time;
          std::cout << "Enter spotting time." << std::endl;
          std::cin >> spotting_time;
          pathFilename = pathFilename + "_SpotTime_" + std::to_string((int)(spotting_time*10));
          timedepsolver->setSpottingTimes(std::vector({spotting_time}));
        } else if (specify_spotting_times == "N") {
          std::cout << "Ok." << std::endl;
        } else {
          std::cout << "Invalid response. Terminating." << std::endl;
          assert(false);
        }

        const bool finalStage = true;
        timedepsolver->traceHJB(finalStage, pathFilename, read_from_file);
      } else {
        bool starting_locations;
        std::string many_starting_locations;
        std::cout << "Would you like to include many starting locations (Y/N)?" << std::endl;
        std::cin >> many_starting_locations;

        if (many_starting_locations == "Y") {
          std::cout << " Tracing paths..." << std::endl;
          auto tstart = std::chrono::high_resolution_clock::now();
          timedepsolver->traceNDifferentPaths(num_paths, pathFilename, 
                                              read_from_file, initial_mode);
          /** Stop timer */
          auto tend = std::chrono::high_resolution_clock::now();
          std::cout << "Complete, took "
                    << std::chrono::duration_cast<std::chrono::milliseconds>(tend-tstart).count()
                    << " milliseconds" << std::endl;

        } else if (many_starting_locations == "N") {
          std::cout << " Tracing paths..." << std::endl;
          auto tstart = std::chrono::high_resolution_clock::now();
          timedepsolver->traceNPaths(num_paths, pathFilename);
          /** Stop timer */
          auto tend = std::chrono::high_resolution_clock::now();
          std::cout << "Complete, took "
                    << std::chrono::duration_cast<std::chrono::milliseconds>(tend-tstart).count()
                    << " milliseconds" << std::endl;
        } else {
          std::cout << "Invalid response. Terminating." << std::endl;
          assert(false);
        }
      }

      std::string more_path_tracing;
      std::cout << "Done. Would you like to trace more paths? (Y/N)" << std::endl;
      std::cin >> more_path_tracing;
      if (more_path_tracing == "Y") {
        path_tracing = true;
      } else if (more_path_tracing == "N") {
        path_tracing = false;
      } else {
        std::cout << "Invalid response. Terminating." << std::endl;
        assert(false);
      }
    }
  } else {
    timedepsolver->solveMultiStage(filename, stencil);
  }

  /** Stop timer */
  auto t2 = std::chrono::high_resolution_clock::now();
  std::cout << "Complete, took "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
            << " milliseconds" << std::endl;
  return 0;
}
