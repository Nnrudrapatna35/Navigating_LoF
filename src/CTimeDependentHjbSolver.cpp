/*
* ==============================================================================
*
*  Copyright (C) 2024 Marissa Gee Nicolas Gonzalez Granda Sunay Joshi 
                      Nagaprasad Rudrapatna Anne Somalwar 
*
*  This program is free software: you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation, either version 3 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*
* ------------------------------------------------------------------------------
* File: CTimeDependentHjbSolver.cpp
*
* Authors: REU 2022 Landscape of Fear Group
*   (based on code by Elliot Cartee, Qianli Song, Lexiao Lai)
*     (based on code by Marc Aurèle Gilles and Zachary Clawson)
*
* Description: This class is for solving the time-dependent HJB equation.
* It assumes that fGrid contains the speed and cost functions.
* It assumes a 2D regular grid with equal spacing in both directions.
* The target is assumed to be single points.
* Updates from both 5-point and 9-point stencils are possible,
*     depends on setting in GlobalConfigurations.hpp
* It depends on the CTimeGrid class.
* (See also CTimeDependentHjbSolver.hpp)
*
* ==============================================================================
*/
/* ------ Libraries ---------------------------------------------------------*/
#include <cmath>
#include <string>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <boost/progress.hpp>
#include <iomanip>



/* ------ Project-specific header files -------------------------------------*/
#include "CTimeDependentHjbSolver.hpp"
#include "WritetoFile.hpp"

/* ------ Namespaces --------------------------------------------------------*/
using namespace memory;
using namespace std;

/*==============================================================================
  Constructor
==============================================================================*/

// Function-based constructor
CTimeDependentHjbSolver::CTimeDependentHjbSolver(function<double(double, double, double, double)> aCost1,
                          function<double(double, double, double, double)> aCost2,
                          function<double(double, double, double, double)> aSpeed1,
                          function<double(double, double, double, double)> aSpeed2,
                          shared_ptr<CFoodDensity> aFoodDensity,
                          shared_ptr<CTerrain> aTerrain,
                          const double aEnergyThreshold,
                          const int aN, const int aNe,
                          const double aPhysMin, const double aPhysMax,
                          const double aEnergyMin, const double aEnergyMax,
                          const double aTimeMax, const vector<double> aSource,
                          const double aInitialEnergy,
                          const objective_t aObjective, const int aNStages,
                          const bool aFullGrid, const bool aReturnHome, 
                          const bool aAllowModeSwitches) {

  assert(aPhysMax >= aPhysMin);
  assert(aEnergyMax >= aEnergyMin);
  assert(aEnergyMax >= aEnergyThreshold);
  assert(aEnergyThreshold >= aEnergyMin);
  assert(aTimeMax > 0.0);

  
  fObjective = aObjective;
  fNStages = aNStages;
  fNSurvivedStages = 0;
  fTimeMax = aTimeMax;
  fEnergyThreshold = aEnergyThreshold;
  fFullGrid = aFullGrid;
  fReturnHome = aReturnHome;
  fFoodDensity = aFoodDensity;
  fTerrain = aTerrain;

  // Computing dt CFL
  const double h = (aPhysMax - aPhysMin) / (aN-1);
  const double de = (aEnergyMax - aEnergyMin) / (aNe-1);
  double maxVal = -1.0 * INF;
  double temp = 0.0;
  const double stageLength = aTimeMax/aNStages;
  // Print length of horizon *per stage*
  cout << "S: " << stageLength << endl;

  for (int i = 0; i < aN; i++) {
    for (int j = 0; j < aN; j++) {
      // for (int k = 0; k < aNe; k++) {
      // In general we allow these values to depend on e, however we haven't 
      // used those examples, so to save time only iterate through x and y

      const double x = aPhysMin + (double)i * h;
      const double y = aPhysMin + (double)j * h;
      // const double e = aEnergyMin + (double)k * de;

      // Current values; plug in 0 for fourth argument!
      const double f1 = aSpeed1(x,y,0,0);
      const double f2 = aSpeed2(x,y,0,0);

      const double K1 = aCost1(x,y,0,0);
      const double K2 = aCost2(x,y,0,0);

      const double F = fFoodDensity->getFoodAccumulated(i,j);

      const double cond1 = f1 + (h * abs(F - K1)) / (SQRT2 * de);
      const double cond2 = f2 + (h * K2) / (SQRT2 * de);

      temp = max(cond1, cond2);
      maxVal = max(temp, maxVal);
      // }
    }
  }

  const double dt = h / (SQRT2 * maxVal); // Time spacing
  const int nt = ceil((stageLength) / dt + 1); // Compute number of time gridpoints *per stage*!

  cout << "nt: " << nt << endl;
  cout << "dt: " << fTimeMax/(nt-1) << endl;

  // Build grids; pass in time *per stage* stageLength!
  if (fFullGrid) {
    fGrid1 = make_shared<CTimeGrid> (CTimeGrid(aCost1, aSpeed1, aN, aNe, nt,
                                      aPhysMin, aPhysMax, aEnergyMin,
                                      aEnergyMax, stageLength));
    fGrid2 = make_shared<CTimeGrid> (CTimeGrid(aCost2, aSpeed2, aN, aNe, nt,
                                      aPhysMin, aPhysMax, aEnergyMin, 
                                      aEnergyMax, stageLength));
  } else {
    assert(aNStages == 1);
    fGrid1 = make_shared<CTimeGrid> (CTimeGrid(aCost1, aSpeed1, aN, aNe, nt,
                                      aPhysMin, aPhysMax, aEnergyMin,
                                      aEnergyMax, aTimeMax, fFullGrid, nt-2));
    fGrid2 = make_shared<CTimeGrid> (CTimeGrid(aCost2, aSpeed2, aN, aNe, nt,
                                      aPhysMin, aPhysMax, aEnergyMin, 
                                      aEnergyMax, aTimeMax, fFullGrid, nt-2));
  }

  // Build tracer
  // Note that InitialMode is always 1 to start! Also, pass in aObjective
  fTracer = make_shared<CTimeDependentTracer> (CTimeDependentTracer(fGrid1, fGrid2,
                                                fFoodDensity, fTerrain, aSource, 
                                                1, aInitialEnergy,
                                                aAllowModeSwitches));


  // Array of ReturnTimes
  fReturnTimes = make_shared<array2D_t<double>>(allocateArray2D<double>(aN,aN));

  // Call CEikonalSolver.march
  CEikonalSolver eikonalSolver = CEikonalSolver(aSpeed1, aTerrain, aN, aPhysMin, aPhysMax);
  fReturnTimes = eikonalSolver.march();
}

// an alternate function-based constructor that allows us to specify the number of timesteps
CTimeDependentHjbSolver::CTimeDependentHjbSolver(function<double(double, double, double, double)> aCost1,
                          function<double(double, double, double, double)> aCost2,
                          function<double(double, double, double, double)> aSpeed1,
                          function<double(double, double, double, double)> aSpeed2,
                          shared_ptr<CFoodDensity> aFoodDensity,
                          shared_ptr<CTerrain> aTerrain,
                          const double aEnergyThreshold,
                          const int aNt, const int aN, const int aNe,
                          const double aPhysMin, const double aPhysMax,
                          const double aEnergyMin, const double aEnergyMax,
                          const double aTimeMax, const vector<double> aSource,
                          const double aInitialEnergy,
                          const objective_t aObjective, const int aNStages,
                          const bool aFullGrid, const bool aReturnHome, 
                          const bool aAllowModeSwitches) {

  assert(aPhysMax >= aPhysMin);
  assert(aEnergyMax >= aEnergyMin);
  assert(aEnergyMax >= aEnergyThreshold);
  assert(aEnergyThreshold >= aEnergyMin);
  assert(aTimeMax > 0.0);

  
  fObjective = aObjective;
  fNStages = aNStages;
  fNSurvivedStages = 0;
  fTimeMax = aTimeMax;
  fEnergyThreshold = aEnergyThreshold;
  fFullGrid = aFullGrid;
  fReturnHome = aReturnHome;
  fFoodDensity = aFoodDensity;
  fTerrain = aTerrain;

  // Computing dt CFL
  const double h = (aPhysMax - aPhysMin) / (aN-1);
  const double de = (aEnergyMax - aEnergyMin) / (aNe-1);
  double maxVal = -1.0 * INF;
  double temp = 0.0;
  const double stageLength = aTimeMax/aNStages;
  // Print length of horizon *per stage*
  cout << "S: " << stageLength << endl;

  for (int i = 0; i < aN; i++) {
    for (int j = 0; j < aN; j++) {
      // for (int k = 0; k < aNe; k++) {
      // In general we allow these values to depend on e, however we haven't 
      // used those examples, so to save time only iterate through x and y

      const double x = aPhysMin + (double)i * h;
      const double y = aPhysMin + (double)j * h;
      // const double e = aEnergyMin + (double)k * de;

      // Current values; plug in 0 for fourth argument!
      const double f1 = aSpeed1(x,y,0,0);
      const double f2 = aSpeed2(x,y,0,0);

      const double K1 = aCost1(x,y,0,0);
      const double K2 = aCost2(x,y,0,0);

      const double F = fFoodDensity->getFoodAccumulated(i,j);

      const double cond1 = f1 + (h * abs(F - K1)) / (SQRT2 * de);
      const double cond2 = f2 + (h * K2) / (SQRT2 * de);

      temp = max(cond1, cond2);
      maxVal = max(temp, maxVal);
      // }
    }
  }

  const double dt = h / (SQRT2 * maxVal); // Time spacing
  int nt = ceil((aTimeMax/aNStages) / dt + 1); // Compute number of time gridpoints *per stage*!

  assert(aNt >= nt); //ensure that our chose nt satisfies the CFL
  nt = aNt;

  cout << "nt: " << nt << endl;
  cout << "dt: " << fTimeMax/(nt-1) << endl;

  // Build grids; pass in time *per stage* stageLength!
  if (fFullGrid) {
    fGrid1 = make_shared<CTimeGrid> (CTimeGrid(aCost1, aSpeed1, aN, aNe, nt,
                                      aPhysMin, aPhysMax, aEnergyMin,
                                      aEnergyMax, stageLength));
    fGrid2 = make_shared<CTimeGrid> (CTimeGrid(aCost2, aSpeed2, aN, aNe, nt,
                                      aPhysMin, aPhysMax, aEnergyMin, 
                                      aEnergyMax, stageLength));
  } else {
    assert(aNStages == 1);
    fGrid1 = make_shared<CTimeGrid> (CTimeGrid(aCost1, aSpeed1, aN, aNe, nt,
                                      aPhysMin, aPhysMax, aEnergyMin,
                                      aEnergyMax, aTimeMax, fFullGrid, nt-2));
    fGrid2 = make_shared<CTimeGrid> (CTimeGrid(aCost2, aSpeed2, aN, aNe, nt,
                                      aPhysMin, aPhysMax, aEnergyMin, 
                                      aEnergyMax, aTimeMax, fFullGrid, nt-2));
  }

  // Build tracer
  // Note that InitialMode is always 1 to start! Also, pass in aObjective
  fTracer = make_shared<CTimeDependentTracer> (CTimeDependentTracer(fGrid1, fGrid2,
                                                fFoodDensity, fTerrain, aSource, 
                                                1, aInitialEnergy,
                                                aAllowModeSwitches));


  // Array of ReturnTimes
  fReturnTimes = make_shared<array2D_t<double>>(allocateArray2D<double>(aN,aN));

  // Call CEikonalSolver.march
  CEikonalSolver eikonalSolver = CEikonalSolver(aSpeed1, aTerrain, aN, aPhysMin, aPhysMax);
  fReturnTimes = eikonalSolver.march();
}


// Array-based constructor
CTimeDependentHjbSolver::CTimeDependentHjbSolver(shared_ptr<array2D_t<double>> aCost1,
                          shared_ptr<array2D_t<double>> aCost2,
                          shared_ptr<array2D_t<double>> aSpeed1,
                          shared_ptr<array2D_t<double>> aSpeed2,
                          shared_ptr<CFoodDensity> aFoodDensity,
                          shared_ptr<CTerrain> aTerrain,
                          const double aEnergyThreshold,
                          const int aN, const int aNe,
                          const double aPhysMin, const double aPhysMax,
                          const double aEnergyMin, const double aEnergyMax,
                          const double aTimeMax, const vector<double> aSource,
                          const double aInitialEnergy,
                          const objective_t aObjective, const int aNStages,
                          const bool aFullGrid, const bool aReturnHome, 
                          const bool aAllowModeSwitches) { 

  assert(aPhysMax >= aPhysMin);
  assert(aEnergyMax >= aEnergyMin);
  assert(aEnergyMax >= aEnergyThreshold);
  assert(aEnergyThreshold >= aEnergyMin);
  assert(aTimeMax > 0.0);

  // Print length of horizon *per stage*
  cout << "T: " << aTimeMax / aNStages << endl;
  fObjective = aObjective;
  fNStages = aNStages;
  fNSurvivedStages = 0;
  fTimeMax = aTimeMax;
  fReturnHome = aReturnHome;
  fFoodDensity = aFoodDensity;
  fTerrain = aTerrain;
  fEnergyThreshold = aEnergyThreshold;
  fFullGrid = aFullGrid;

  // Computing dt CFL
  const double h = (aPhysMax - aPhysMin) / (aN-1);
  const double de = (aEnergyMax - aEnergyMin) / (aNe-1);
  double maxVal = -1.0 * INF;
  double temp = 0.0;
  const double stageLength = aTimeMax/aNStages;
  // Print length of horizon *per stage*
  cout << "S: " << stageLength << endl;

  for (int i = 0; i < aN; i++) {
    for (int j = 0; j < aN; j++) {
      // for (int k = 0; k < aNe; k++) {
      // Current values;
      const double f1 = (*aSpeed1)[i][j];
      const double f2 = (*aSpeed2)[i][j];

      const double K1 = (*aCost1)[i][j];
      const double K2 = (*aCost2)[i][j];

      const double F = fFoodDensity->getFoodAccumulated(i,j);

      const double cond1 = f1 + (h * abs(F - K1)) / (SQRT2 * de);
      const double cond2 = f2 + (h * K2) / (SQRT2 * de);

      temp = max(cond1, cond2);
      maxVal = max(temp, maxVal);
      // }
    }
  }

  double dt = h / (SQRT2 * maxVal); // Time spacing
  int nt = ceil((stageLength) / dt + 1); // Compute number of time gridpoints *per stage*!

  cout << "nt: " << nt << endl;
  cout << "dt: " << fTimeMax/(nt-1) << endl;

  // Build grids; pass in time *per stage* stageLength!
  if (fFullGrid) {
    fGrid1 = make_shared<CTimeGrid> (CTimeGrid(aCost1, aSpeed1, aN, aNe, nt,
                                      aPhysMin, aPhysMax, aEnergyMin,
                                      aEnergyMax, stageLength));
    fGrid2 = make_shared<CTimeGrid> (CTimeGrid(aCost2, aSpeed2, aN, aNe, nt,
                                      aPhysMin, aPhysMax, aEnergyMin, 
                                      aEnergyMax, stageLength));
  } else {
    assert(aNStages == 1);
    fGrid1 = make_shared<CTimeGrid> (CTimeGrid(aCost1, aSpeed1, aN, aNe, nt,
                                      aPhysMin, aPhysMax, aEnergyMin,
                                      aEnergyMax, aTimeMax, fFullGrid, nt-2));
    fGrid2 = make_shared<CTimeGrid> (CTimeGrid(aCost2, aSpeed2, aN, aNe, nt,
                                      aPhysMin, aPhysMax, aEnergyMin, 
                                      aEnergyMax, aTimeMax, fFullGrid, nt-2));
  }

  // Build tracer
  // Note that InitialMode is always 1 to start! Also, pass in aObjective
  fTracer = make_shared<CTimeDependentTracer> (CTimeDependentTracer(fGrid1, fGrid2,
                                                fFoodDensity, fTerrain, aSource, 
                                                1, aInitialEnergy,
                                                aAllowModeSwitches));


  // Array of ReturnTimes
  fReturnTimes = make_shared<array2D_t<double>>(allocateArray2D<double>(aN,aN));

  // Call CEikonalSolver.march
  CEikonalSolver eikonalSolver = CEikonalSolver(aSpeed1, aTerrain, aN, aPhysMin, aPhysMax);
  fReturnTimes = eikonalSolver.march();
}


/*==============================================================================
  Main function for solving time-dependent HJB equation
==============================================================================*/
void CTimeDependentHjbSolver::solveHJB(const int aCurrStage, 
                                       const string aFilename, 
                                       const stencil_t aStencilType) {
  /* Get grid dimensions using fGrid1 */
  const int nx = fGrid1->getGridSizeX();
  const int ny = fGrid1->getGridSizeY();
  const int ne = fGrid1->getGridSizeE();
  const int nt = fGrid1->getGridSizeT();

  /* Compute terminal condition */
  if (fReturnHome) {
    expandingGoHomeTerminal(aCurrStage);
  } else {
    cout << "Vanilla terminal" << endl;
    vanillaTerminal();
  }
  // cout << "Receding horizon terminal condition." << endl;
  // recedingHorizonTerminal(aCurrStage, aStencilType);

  if (!fFullGrid) {
    fGrid1->writeSliceToFile(aFilename + "_Stage" + to_string(aCurrStage) + "_Mode1", nt-1);
    fGrid2->writeSliceToFile(aFilename + "_Stage" + to_string(aCurrStage) + "_Mode2", nt-1);
  }

  /* Enforce value = 0 if energy = 0. */
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      if (fFullGrid) {
        for (int l = 0; l < nt; ++l) {
          fGrid1->setValue(i,j,0,l,0);
          fGrid2->setValue(i,j,0,l,0);
        }
      } else {
        fGrid1->setValue(i,j,0,nt-1,0);
        fGrid1->setValue(i,j,0,nt-2,0);
        fGrid2->setValue(i,j,0,nt-1,0);
        fGrid2->setValue(i,j,0,nt-2,0);
      }
    }
  }

  if (aCurrStage == 0) {
    // If this is the first stage write the initial food information to file
    fFoodDensity->writeFoodGridToFile(aFilename);
  }

  /* Evolve value function backward in time */
  /* *Start* at k = 1 to avoid overwriting boundary condition! */
  cout << "Solving backwards in time." << endl;
  boost::progress_display show_progress(nt-1);
  for (int l = nt - 1; l > 0; --l) {
    for (int i = 0; i < nx; ++i) {
      for (int j = 0; j < ny; ++j) {
        for (int k = 1; k < ne; ++k) {
          if (fTerrain->isObstacle(i,j)) {
            fGrid1->setValue(i,j,k,l-1,0);
            fGrid2->setValue(i,j,k,l-1,0);
          } else {
            double u_new1;
            double u_new2;
            if (aStencilType == FIVE_POINT) {
              u_new1 = update_5pt(i,j,k,l,1);
              u_new2 = update_5pt(i,j,k,l,2);
            } else if (aStencilType == NINE_POINT){
              u_new1 = update_9pt(i,j,k,l,1);
              u_new2 = update_9pt(i,j,k,l,2);
            } else {
              cout << "Invalid stencil type" << endl;
              assert(false);
            }

            /* Set values for each grid */
            fGrid1->setValue(i,j,k,l-1,u_new1);
            fGrid2->setValue(i,j,k,l-1,u_new2);
          }
        }
      }
    }
    if (!fFullGrid) {
      fGrid1->writeSliceToFile(aFilename + "_Stage" + to_string(aCurrStage) + "_Mode1", l-1);
      fGrid2->writeSliceToFile(aFilename + "_Stage" + to_string(aCurrStage) + "_Mode2", l-1);

      fGrid1->advanceSlicesBackward();
      fGrid2->advanceSlicesBackward();
    }
    ++show_progress;
  }
  std::cout << std::endl;

  return;
}
// Solve multistage
void CTimeDependentHjbSolver::solveMultiStage(const string aFilename,
                                              const stencil_t aStencilType) {
  // Keep track of nSurvivedStages
  fNSurvivedStages = 0;

  // Iterate through stages
  for (int currStage = 0; currStage < fNStages; ++currStage) {
    cout << "Stage: " << currStage << endl;

    // Check if we are in the last stage
    const bool isFinalStage = (currStage == fNStages-1);

    // Solve the HJB equation for this stage
    solveHJB(currStage, aFilename, aStencilType);

    // // Trace one realization of the animal's path without food depletion
    // traceHJB(isFinalStage, aFilename + "_Stage" + to_string(currStage) + "_NoDepletion", false);

    // Trace one realization of the animal's path
    traceHJB(isFinalStage, aFilename + "_Stage" + to_string(currStage));

    // Update food accumulated
    fFoodDensity->updateFoodAccumulated();

    // get mode list, energy list
    Path currentPath = fTracer->getPath();
    vector<int> currentModeList = fTracer->getModeList();
    vector<double> currentEnergyList = fTracer->getEnergyList();

    // set initial conditions of tracer
    const double currentX = (currentPath.back()).x;
    const double currentY = (currentPath.back()).y;

    fTracer->setSource(currentX, currentY);
    fTracer->setInitialMode(currentModeList.back());
    fTracer->setInitialEnergy(currentEnergyList.back());

    // concatenate to the big lists
    if (currStage == 0) {
      fBigPath.insert(fBigPath.end(), currentPath.begin(), currentPath.end());
      fBigModeList.insert(fBigModeList.end(), currentModeList.begin(), currentModeList.end());
      fBigEnergyList.insert(fBigEnergyList.end(), currentEnergyList.begin(), currentEnergyList.end());
    } else {
      fBigPath.insert(fBigPath.end(), ++currentPath.begin(), currentPath.end());
      fBigModeList.insert(fBigModeList.end(), ++currentModeList.begin(), currentModeList.end());
      fBigEnergyList.insert(fBigEnergyList.end(), ++currentEnergyList.begin(), currentEnergyList.end());
    }

    writeResults(aFilename, currStage);

    // exit loop if currentModeList.back() is NOT 1 or 2!
    if (!(currentModeList.back()==1 || currentModeList.back()==2)) {
      break;
    }
    fNSurvivedStages++; // Increment
  }

  // update fNSurvivedStages
  cout << "The animal survived " << fNSurvivedStages << " stages!" << endl;

}

/*==============================================================================
  Path tracing functions
==============================================================================*/
/* Trace optimal path and write result to file */
void CTimeDependentHjbSolver::traceHJB(const bool aIsFinalStage,
                                       const string aFilename,
                                       const bool aReadFromFile,
                                       const bool aClassifyPath) {
  if (fFullGrid) {
    if (aReadFromFile) {
      fGrid1->readValuesFromFile(aFilename + "_Stage0_Mode1Value");
      fGrid2->readValuesFromFile(aFilename + "_Stage0_Mode2Value");
    }
    fTracer->tracePath(aIsFinalStage); 
  } else {
    if (aReadFromFile) {
      fTracer->tracePath(aIsFinalStage, aFilename + "_Stage0");
    } else {
      fTracer->tracePath(aIsFinalStage, aFilename);
    }
  }
   
  fTracer->printPathToFile(aFilename);
  if (aClassifyPath) {
    fTracer->classifyPath(aFilename);
  }
}

/* Trace many optimal paths and write results to file */
void CTimeDependentHjbSolver::traceNPaths(const int aNumPaths,
                                          const string aFilename,
                                          const bool aReadFromFile,
                                          const int aInitialMode) {
  if (fFullGrid) {
    fTracer->tracePathsParallel(aFilename, aNumPaths, aInitialMode); 
  } else {
    if (aReadFromFile) {
      fTracer->tracePathsParallel(aFilename + "_Stage0", aNumPaths, aInitialMode);
    } else {
      fTracer->tracePathsParallel(aFilename, aNumPaths, aInitialMode);
    }
  }
}

void CTimeDependentHjbSolver::resetFoodDensity() {
  fFoodDensity->resetFoodDensity();
}

/* Trace many optimal paths and write results to file */
void CTimeDependentHjbSolver::traceNDifferentPaths(const int aNumPaths,
                                                   const string aFilename,
                                                   const bool aReadFromFile,
                                                   const int aInitialMode) {
  if (fFullGrid) {
    fTracer->traceDifferentPathsParallel(aFilename, aNumPaths, aInitialMode); 
  } else {
    cout << "Need to update filepaths for this to work." << endl;
    assert(false);
    if (aReadFromFile) {
      fTracer->traceDifferentPathsParallel(aFilename + "_Stage0", aNumPaths,
                                           aInitialMode);
    } else {
      fTracer->traceDifferentPathsParallel(aFilename, aNumPaths, aInitialMode);
    }
  }
}

/*==============================================================================
  Terminal condition functions
==============================================================================*/
/* Need not return to the home base */
void CTimeDependentHjbSolver::vanillaTerminal() {
  /* Get grid dimensions */
  const int nx = fGrid1->getGridSizeX();
  const int ny = fGrid1->getGridSizeY();
  const int ne = fGrid1->getGridSizeE();
  const int nt = fGrid1->getGridSizeT();

  /* Enforce terminal conditions on value function */
  for (int i = 0; i < nx; ++i) {
    for(int j = 0; j < ny; ++j) {
      for (int k = 0; k < ne; ++k) {
        double value = 0;
        const double e = fGrid1->eGridToPhysical(k);

        /* Compute terminal reward */
        if (fTerrain->isObstacle(i,j)){
          value = 0;
        } else if (fObjective == SURVIVAL && e > fEnergyThreshold) {
          value = 1;
        } else if (fObjective == ENERGY && e > 0) {
          value = e;
        } else if (fObjective == UTILITY) {
          value = fFoodDensity->getUtility(e);
        } else {
          value = 0;
        }

        fGrid1->setValue(i,j,k,nt-1,value);
        fGrid2->setValue(i,j,k,nt-1,value);
      }
    }
  }
}

/* Must be able to reach home base in final stage */
void CTimeDependentHjbSolver::goHomeTerminal() {
  /* Get grid dimensions */
  const int nx = fGrid1->getGridSizeX();
  const int ny = fGrid1->getGridSizeY();
  const int ne = fGrid1->getGridSizeE();
  const int nt = fGrid1->getGridSizeT();

  /* Enforce terminal conditions on value function */
  for (int i = 0; i < nx; ++i) {
    for(int j = 0; j < ny; ++j) {
      for (int k = 0; k < ne; ++k) {
        double value = 0;
        const double e = fGrid1->eGridToPhysical(k);

        /* Only compute reward if animal is in the home base */
        if (fTerrain->isHomeBase(i,j) == false){
          value = 0;
        } else if (fTerrain->isObstacle(i,j)){
          value = 0;
        } else if (fObjective == SURVIVAL && e > fEnergyThreshold) {
          value = 1;          
        } else if (fObjective == ENERGY) {
          value = e;
        } else if (fObjective == UTILITY) {
          value = fFoodDensity->getUtility(e);
        } else {
          value = 0;
        }

        fGrid1->setValue(i,j,k,nt-1,value);
        fGrid2->setValue(i,j,k,nt-1,value);
      }
    }
  }
}

// Must be able to reach home base at the end of each stage
void CTimeDependentHjbSolver::expandingGoHomeTerminal(const int aCurrStage) {
  /* get grid dimensions */
  const int nx = fGrid1->getGridSizeX();
  const int ny = fGrid1->getGridSizeY();
  const int ne = fGrid1->getGridSizeE();
  const int nt = fGrid1->getGridSizeT();

  /* enforce terminal conditions on value function */
  /* NOTE: we are starting at *k=0* */
  /* all other values are 0 because grids are zeroed out in constructor */
  for (int i = 0; i < nx; ++i) {
    for(int j = 0; j < ny; ++j) {
      for (int k = 0; k < ne; ++k) {
        double value = 0;
        const double e = fGrid1->eGridToPhysical(k);

        /* Only compute reward if animal can reach the home base */
        if ((*fReturnTimes)[i][j] > (double)(fNStages - 1 - aCurrStage) * fTimeMax / fNStages ) {
          value = 0;
        } else if (fObjective == SURVIVAL && e > fEnergyThreshold) {
          value = 1;
        } else if (fObjective == ENERGY && e > 0) {
          value = e;
        } else if (fObjective == UTILITY) {
          value = fFoodDensity->getUtility(e);
        } else {
          value = 0;
        }

        fGrid1->setValue(i,j,k,nt-1,value);
        fGrid2->setValue(i,j,k,nt-1,value);
      }
    }
  }
}

void CTimeDependentHjbSolver::recedingHorizonTerminal(const int aCurrStage, 
                                                      const stencil_t aStencilType) {
  /* get grid dimensions */
  const int nx = fGrid1->getGridSizeX();
  const int ny = fGrid1->getGridSizeY();
  const int ne = fGrid1->getGridSizeE();
  const int nt = fGrid1->getGridSizeT();

  // cout << "fTime is " << fTime << endl;
  const float dt = fTimeMax/(fNStages*(nt-1));
  // cout << "dt is " << dt << endl;
  const float timeRemaining = (float)(fNStages - 1 - aCurrStage) * fTimeMax / fNStages;
  const int ntRemaining = round(timeRemaining/dt) + 1;
  cout << "Time remaining " << timeRemaining << endl;
  // cout << "nt remaining is " << ntRemaining << endl;
  cout << "Numerical time remaining " << dt*(ntRemaining-1) << endl;

  shared_ptr<array4D_t<float>> tempValueArray1 =
    make_shared<array4D_t<float>>(allocateArray4D<float>(nx, ny, ne, ntRemaining));
  shared_ptr<array4D_t<float>> tempValueArray2 =
    make_shared<array4D_t<float>>(allocateArray4D<float>(nx, ny, ne, ntRemaining));

  // Terminal condition for our terminal condition
  for (int i = 0; i < nx; ++i) {
    for(int j = 0; j < ny; ++j) {
      for (int k = 0; k < ne; ++k) {
        float value = 0;
        float e = fGrid1->eGridToPhysical(k);

        /* Only compute reward if animal is in the home base */
        // if (fTerrain->isHomeBase(i,j) == false){
        //   value = 0;
        // } else 
        if (fTerrain->isObstacle(i,j)){
          value = 0;
        } else if (fObjective == SURVIVAL && e > fEnergyThreshold) {
          value = 1;          
        } else if (fObjective == ENERGY) {
          value = e;
        } else if (fObjective == UTILITY) {
          value = fFoodDensity->getUtility(e);
        } else {
          value = 0;
        }

        (*tempValueArray1)[i][j][k][ntRemaining-1] = value;
        (*tempValueArray2)[i][j][k][ntRemaining-1] = value;

      }
    }
  }

  /* Enforce value = 0 if energy = 0. */
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int l = 0; l < ntRemaining; ++l) {
        (*tempValueArray1)[i][j][0][l] = 0;
        (*tempValueArray2)[i][j][0][l] = 0;
      }
    }
  }

  /* enforce terminal conditions on value function */
  /* NOTE: we are starting at *k=0* */
  /* all other values are 0 because grids are zeroed out in constructor */
  cout << "Computing terminal condition." << endl;
  boost::progress_display show_progress(ntRemaining-1);
  for (int l = ntRemaining - 1; l > 0; --l) {
    for (int i = 0; i < nx; ++i) {
      for(int j = 0; j < ny; ++j) {
        const float foodAccumulated = fFoodDensity->getFoodAccumulated(i,j);

        for (int k = 1; k < ne; ++k) {
          if (fTerrain->isObstacle(i,j)) {
            (*tempValueArray1)[i][j][k][l-1] = 0;
            (*tempValueArray2)[i][j][k][l-1] = 0;
          } else {
            float u_new1;
            float u_new2;
            /* Note that update_5pt takes mode and foodAccumulated */
            if (aStencilType == FIVE_POINT) {
              u_new1 = update_5pt_terminal(tempValueArray1, tempValueArray2, i,j,k,l,1,foodAccumulated);
              u_new2 = update_5pt_terminal(tempValueArray1, tempValueArray2, i,j,k,l,2,0); // Zero food is accumulated in mode 2
            } else if (aStencilType == NINE_POINT){
              assert(false);
              // u_new1 = update_9pt(i,j,k,l,1,foodAccumulated);
              // u_new2 = update_9pt(i,j,k,l,2,0); // Zero food is accumulated in mode 2
            } else {
              cout << "Invalid stencil type" << endl;
              assert(false);
            }

            /* Set values for each grid */
            (*tempValueArray1)[i][j][k][l-1] = u_new1;
            (*tempValueArray2)[i][j][k][l-1] = u_new2;
          }
        }
      }
    }
    ++show_progress;
  }

  for (int i = 0; i < nx; ++i) {
    for(int j = 0; j < ny; ++j) {
      for (int k = 0; k < ne; ++k) {
        const float value1 = (*tempValueArray1)[i][j][k][0];
        const float value2 = (*tempValueArray2)[i][j][k][0];
        fGrid1->setValue(i,j,k,nt-1,value1);
        fGrid2->setValue(i,j,k,nt-1,value2);
      }
    }
  }
}

/*==============================================================================
  Five-point stencil HJB update
==============================================================================*/
double CTimeDependentHjbSolver::update_5pt(const int aI, const int aJ,
                                           const int aK, const int aL,
                                           const int aMode) const
{

  /* Grid sizes */
  const int nx = fGrid1->getGridSizeX();
  const int ny = fGrid1->getGridSizeY();
  const int ne = fGrid1->getGridSizeE();

  const double h = fGrid1->getH();
  const double he = fGrid1->getE(); /* NOTE he = energy spacing!... */

  // assign current grid based on mode
  shared_ptr<CTimeGrid> currentGrid;
  if (aMode == 1) {
    currentGrid = fGrid1;
  } else if (aMode == 2) {
    currentGrid = fGrid2;
  }

  // code that is common to both modes
  const float u = currentGrid->getValue(aI,aJ,aK,aL); /* moved INTO if statement */

  /* Calculate directional derivatives of value function */
  double dplusx, dminusx, dplusy, dminusy, dpluse, dminuse;
  if (aI < nx-1) {
    dplusx  = (u-currentGrid->getValue(aI+1,aJ,aK,aL))/h;
  } else {
    dplusx = 0;
  }
  if (aI > 0) {
    dminusx = (u-currentGrid->getValue(aI-1,aJ,aK,aL))/h;
  } else {
    dminusx = 0;
  }
  if (aJ < ny-1) {
    dplusy  = (u-currentGrid->getValue(aI,aJ+1,aK,aL))/h;
  } else {
    dplusy = 0;
  }
  if (aJ > 0) {
    dminusy = (u-currentGrid->getValue(aI,aJ-1,aK,aL))/h;
  } else {
    dminusy = 0;
  }
  /* NOTE also need e derivative. Divide by he... */
  if (aK < ne-1) {
    dpluse  = (u-currentGrid->getValue(aI,aJ,aK+1,aL))/he;
  } else {
    dpluse = 0;
  }
  if (aK > 0) {
    dminuse = (u-currentGrid->getValue(aI,aJ,aK-1,aL))/he;
  } else {
    dminuse = 0;
  }

  /* Calculate partial derivatives using upwind discretization */
  double u_x, u_y, u_e;
  if (dplusx <= dminusx && dplusx < 0) {
    u_x = -dplusx;
    assert(u_x > 0);
  } else if (dminusx <= dplusx && dminusx < 0) {
    u_x = dminusx;
    assert(u_x < 0);
  } else {
    u_x = 0;
  }
  if (dplusy <= dminusy && dplusy < 0) {
    u_y = -dplusy;
    assert(u_y > 0);
  } else if (dminusy <= dplusy && dminusy < 0) {
    u_y = dminusy;
    assert(u_y < 0);
  } else {
    u_y = 0;
  }

  // need these parameters for energy upwind
  const double K = currentGrid->getCost(aI,aJ,aK,aL);
  double F;
  if (aMode == 1) {
    F = fFoodDensity->getFoodAccumulated(aI,aJ); // we can do this because we pass in zero for mode 2! ...
  } else {
    F = 0;
  }

  /* NOTE changed conditionals... */
  /* ADD does this not change for max vs min?... */
  if (F - K > 0) {
    u_e = -dpluse;
  } else if (F - K < 0) {
    u_e = dminuse;
  } else {
    u_e = 0;
  }

  const double g = currentGrid->getSpeed(aI,aJ,aK,aL);
  const double dt = currentGrid->getDt();

  const float u1 = fGrid1->getValue(aI,aJ,aK,aL);
  const float u2 = fGrid2->getValue(aI,aJ,aK,aL);

  // retrieve these parameters from the grids that they correspond to
  const double mu_s = fTerrain->getSpottingRate(aI,aJ,aK,aL);
  const double mu_k = fTerrain->getKillRate(aI,aJ,aK,aL);
  const double mu_g = fTerrain->getGiveUpRate(aI,aJ,aK,aL);

  double val = 0; // updated value to be returned

  if (aMode == 1) {
    val = u - mu_s*u*dt + mu_s*u2*dt + (F-K)*u_e*dt
          + g*sqrt(u_x*u_x + u_y*u_y)*dt;   // mode 1 HJB update
  } else if (aMode == 2) {
    val = u + mu_g*u1*dt - (mu_g+mu_k)*u*dt - K*u_e*dt
          + g*sqrt(u_x*u_x + u_y*u_y)*dt; // mode 2 HJB update
  }

  return val;
}

/*==============================================================================
  Nine-point stencil HJB update
==============================================================================*/
double CTimeDependentHjbSolver::update_9pt(const int aI, const int aJ,
                                          const int aK, const int aL,
                                          const int aMode) const {
  /* Grid sizes */
  const int nx = fGrid1->getGridSizeX();
  const int ny = fGrid1->getGridSizeY();
  const int ne = fGrid1->getGridSizeE();

  /* Get local variables */
  const double x = fGrid1->xGridToPhysical(aI);
  const double y = fGrid1->yGridToPhysical(aJ);
  const double e = fGrid1->eGridToPhysical(aK);

  const double h = fGrid1->getH();
  const double he = fGrid1->getE(); /* NOTE he = energy spacing!... */

  // assign current grid based on mode
  shared_ptr<CTimeGrid> currentGrid;
  if (aMode == 1) {
    currentGrid = fGrid1;
  }
  else if (aMode == 2) {
    currentGrid = fGrid2;
  }

  const float u = currentGrid->getValue(aI,aJ,aK,aL); /* moved INTO if statement */

  /* Calculate upwind e derivative. */
  double dpluse, dminuse, u_e;
  if (aK < ne-1) {
    dpluse  = (u-currentGrid->getValue(aI,aJ,aK+1,aL))/he;
  } else {
    dpluse = 0;
  }
  if (aK > 0) {
    dminuse = (u-currentGrid->getValue(aI,aJ,aK-1,aL))/he;
  } else {
    dminuse = 0;
  }

  // need these parameters for energy upwind
  const double K = currentGrid->getCost(aI,aJ,aK,aL);
  double F;
  if (aMode == 1) {
    F = fFoodDensity->getFoodAccumulated(aI,aJ);
  } else {
    F = 0;
  }
  /* NOTE changed conditionals... */
  /* ADD does this not change for max vs min?... */
  if (F - K > 0) {
    u_e = -dpluse;
  } else if (F - K < 0) {
    u_e = dminuse;
  } else {
    u_e = 0;
  }

  /* Array encoding nine point stencil */
  const int stencil[8][2] = {{1,0}, {1,1}, {0,1}, {-1,1}, {-1,0}, {-1,-1}, {0,-1}, {1,-1}};

  double u_x, u_y;
  /** Calculate upwind discretization of gradient (8 neighbors)
   *  (See Diagram below) */
  /*
     d[3]    d[2]   d[1]
          o----o----o
          |\   |   /|
          | \  |  / |
          |  \ | /  |
     d[4] o----u----o d[0]
          |  / | \  |
          | /  |  \ |
          |/   |   \|
          o----o----o
     d[5]   d[6]   d[7]
                             */

  /* Calculate the difference in value with all 8 neighbors */
  double dneighbors[8];
  for (int it = 0; it < 8; ++it) {
    const int i1 = stencil[it][0];
    const int j1 = stencil[it][1];
    if (aI +i1 >= 0 && aI + i1 <= nx-1 && aJ + j1 >= 0 && aJ + j1 <= ny-1) {
      dneighbors[it]  = (u - currentGrid->getValue(aI + stencil[it][0], aJ + stencil[it][1], aK, aL))/h;
    } else {
      dneighbors[it] = 0;
    }
  }

  /** i_direct, i_diagonal = index of the smallest direct/diagonal neighbor
   *   i_direct ∈ {0, 2, 4, 6}; i_diagonal ∈ {1, 3, 5, 7} */
  int i_direct = 0, i_diagonal = 1;
  double max_direct = dneighbors[0], max_diagonal = dneighbors[1];
  for (int i = 1; i < 4; ++i) {
    if (dneighbors[i*2] < max_direct) {
      i_direct = i*2;
      max_direct = dneighbors[i*2];
    }
    if (dneighbors[i*2 + 1] < max_diagonal) {
      i_diagonal = i*2 + 1;
      max_diagonal = dneighbors[i*2 + 1];
    }
  }

  /* i_diff indicates if the direct and diagonal neighbors chosen are adjacent */
  /* i_diff also helps encode (+)(-) relations between u_x and u_y in different simplices */
  int i_diff = i_diagonal - i_direct;
  if (i_diff == 7) {
    i_diff = -1;
  }

  /* If i_diff equals 1 or -1, then the direct and diagonal neighbors chosen are adjacent */
  if (i_diff * i_diff == 1) {
    double grad_u[2];
    computeGradient(i_diff, i_direct, max_diagonal, max_direct, grad_u);
    u_x = grad_u[0];
    u_y = grad_u[1];
  } else {
    /* The case where the direct and diagonal neighbors chosen are not adjacent */
    double nom_dir[2], nom_diag[2];
    int diff_dir, diff_diag;

    /* Diagonal neighbor nominates a simplex */
    if (dneighbors[(i_diagonal+1)%8]<=dneighbors[i_diagonal-1]) {
      diff_diag = -1;
      computeGradient(-1, (i_diagonal+1)%8, max_diagonal, dneighbors[(i_diagonal+1)%8], nom_diag);
    } else {
      diff_diag = 1;
      computeGradient(1, i_diagonal-1, max_diagonal, dneighbors[i_diagonal-1], nom_diag);
    }

    /* Direct neighbor nominates a simplex */
    if (dneighbors[i_direct+1]<=dneighbors[(i_direct-1)%8]) {
      diff_dir = 1;
      computeGradient(1, i_direct, dneighbors[i_direct+1], max_direct, nom_dir);
    } else {
      diff_dir = -1;
      computeGradient(-1, i_direct, dneighbors[(i_direct-1)%8], max_direct, nom_dir);
    }

    /* Compare the two candidate simplices nominated by direct and diagonal neighbors */
    if (nom_dir[0] * nom_dir[0] + nom_dir[1] * nom_dir[1] >= nom_diag[0] * nom_diag[0] + nom_diag[1] * nom_diag[1]) {
      u_x = nom_dir[0];
      u_y = nom_dir[1];
      i_diagonal = i_direct + diff_dir;
    } else {
      u_x = nom_diag[0];
      u_y = nom_diag[1];
      i_direct = i_diagonal - diff_diag;
    }
  }

  const double g = currentGrid->getSpeed(aI,aJ,aK,aL);
  const double dt = currentGrid->getDt();

  const float u1 = fGrid1->getValue(aI,aJ,aK,aL);
  const float u2 = fGrid2->getValue(aI,aJ,aK,aL);

  // retrieve these parameters from the grids that they correspond to
  const double phi = fTerrain->getSpottingRate(aI,aJ,aK,aL);
  const double mu_k = fTerrain->getKillRate(aI,aJ,aK,aL);
  const double mu_g = fTerrain->getGiveUpRate(aI,aJ,aK,aL);

  double u_new;
  /* Check upwinding condition */
  if (((i_direct == 0 || i_direct == 4) && abs(u_y) > abs(u_x)) || ((i_direct == 2 || i_direct == 6) && abs(u_y) < abs(u_x))) {
    /* Upwinding condition is not satisfied, use smaller of one-sided updates */
    double u_new1;
    double u_new2;
    if (aMode == 1) {
      u_new1 = u - phi*u*dt + phi*u2*dt + (F-K)*u_e*dt - g*dneighbors[i_direct]*dt;
      u_new2 = u + mu_g*u1*dt - (mu_g+mu_k)*u*dt - K*u_e*dt - g*dneighbors[i_diagonal]*1/sqrt(2)*dt;
    } else if (aMode == 2) {
      u_new1 = u - phi*u*dt + phi*u2*dt + (F-K)*u_e*dt - g*dneighbors[i_direct]*dt;
      u_new2 = u + mu_g*u1*dt - (mu_g+mu_k)*u*dt - K*u_e*dt - g*dneighbors[i_diagonal]*1/sqrt(2)*dt;
    }

    if (u_new1 >= u_new2){
      u_new = u_new1;
    } else {
      u_new = u_new2;
    }
  } else {
    /* Upwinding condition is satisfied, use Euler step */
    if (aMode == 1) {
      u_new = u - phi*u*dt + phi*u2*dt + (F-K)*u_e*dt
            + g*sqrt(u_x*u_x + u_y*u_y)*dt;   // mode 1 HJB update
    } else if (aMode == 2) {
      u_new = u + mu_g*u1*dt - (mu_g+mu_k)*u*dt - K*u_e*dt
            + g*sqrt(u_x*u_x + u_y*u_y)*dt; // mode 2 HJB update
    }
  }

  return u_new;
}

/*==============================================================================
  Helper function for computing gradient using nine-point stencil
================================================================================

               d[3]    d[2]   d[1]
                    o----o----o
                    |\ 3 | 2 /|
                    | \  |  / |
                    | 4\ | /1 |
               d[4] o----u----o d[0]
                    | 5/ | \8 |
                    | /  |  \ |
                    |/ 6 | 7 \|
                    o----o----o
               d[5]   d[6]   d[7]

        Compute the gradients in 8 simplices as follows:

        Simplex 1:  u_x = - d[0];
                    u_y = - d[1] - u_x;

        Simplex 2:  u_y = - d[2];
                    u_x = - d[1] - u_y;

        Simplex 3:  u_y = - d[2];
                    u_x = d[3] + u_y;

        Simplex 4:  u_x = d[4];
                    u_y = - d[3] + u_x;

        Simplex 5:  u_x = d[4];
                    u_y = d[5] - u_x;

        Simplex 6:  u_y = d[6];
                    u_x = d[5] - u_y;

        Simplex 7:  u_y = d[6];
                    u_x = - d[7] + u_y;

        Simplex 8:  u_x = - d[0];
                    u_y = d[7] + u_x;
                                            */
void CTimeDependentHjbSolver::computeGradient(const int i_diff, const int i_direct,
                                              const double v_diag,
                                              const double v_dir,
                                              double (&u)[2]) const {
  const int div = (i_direct+2) / 4;
  if (v_dir > 0) {
    u[0] = 0;
    u[1] = 0;
  }
  else if (v_diag >= 0) {
    if (i_direct == 2 || i_direct == 6) {
      u[1] = pow(-1, div) * v_dir;
      u[0] = 0;
    } else {
      u[0] = pow(-1, div+1) * v_dir;
      u[1] = 0;
    }
  } else {
    const double u_xy = v_diag;
    if (i_direct == 2 || i_direct == 6) {
      u[1] = pow(-1, div) * v_dir;
      if (abs(v_diag) < abs(v_dir)) {
        u[0] = 0;
      } else {
        u[0] = i_diff * pow(-1, div + 1) * (u_xy + pow(-1, div + 1)*u[1]);
      }
    } else {
      u[0] = pow(-1, div+1) * v_dir;
      if (abs(v_diag) < abs(v_dir)) {
        u[1] = 0;
      } else {
        u[1] = i_diff * pow(-1, div + 1) * (u_xy + pow(-1, div)*u[0]);
      }
    }
  }
}


float CTimeDependentHjbSolver::update_5pt_terminal(shared_ptr<array4D_t<float>> aValue1,
                                                   shared_ptr<array4D_t<float>> aValue2,
                                                   const int aI, const int aJ,
                                                   const int aK, const int aL,
                                                   const int aMode,
                                                   const float aFoodAccumulated) const
{

  /* Grid sizes */
  const int nx = fGrid1->getGridSizeX();
  const int ny = fGrid1->getGridSizeY();
  const int ne = fGrid1->getGridSizeE();

  /* Get local variables */
  const float x = fGrid1->xGridToPhysical(aI);
  const float y = fGrid1->yGridToPhysical(aJ);
  const float e = fGrid1->eGridToPhysical(aK);

  const float h = fGrid1->getH();
  const float he = fGrid1->getE(); /* NOTE he = energy spacing!... */

  // assign current grid based on mode
  shared_ptr<array4D_t<float>> currentValueArray;
  shared_ptr<CTimeGrid> currentGrid;
  if (aMode == 1) {
    currentValueArray = aValue1;
    currentGrid = fGrid1;
  }
  else if (aMode == 2) {
    currentValueArray = aValue2;
    currentGrid = fGrid2;
  }

  // code that is common to both modes
  const float u = (*currentValueArray)[aI][aJ][aK][aL]; /* moved INTO if statement */

  /* Calculate directional derivatives of value function */
  float dplusx, dminusx, dplusy, dminusy, dpluse, dminuse;
  if (aI < nx-1) {
    dplusx  = (u-(*currentValueArray)[aI+1][aJ][aK][aL])/h;
  } else {
    dplusx = 0;
  }
  if (aI > 0) {
    dminusx = (u-(*currentValueArray)[aI-1][aJ][aK][aL])/h;
  } else {
    dminusx = 0;
  }
  if (aJ < ny-1) {
    dplusy  = (u-(*currentValueArray)[aI][aJ+1][aK][aL])/h;
  } else {
    dplusy = 0;
  }
  if (aJ > 0) {
    dminusy = (u-(*currentValueArray)[aI][aJ-1][aK][aL])/h;
  } else {
    dminusy = 0;
  }
  /* NOTE also need e derivative. Divide by he... */
  if (aK < ne-1) {
    dpluse  = (u-(*currentValueArray)[aI][aJ][aK+1][aL])/he;
  } else {
    dpluse = 0;
  }
  if (aK > 0) {
    dminuse = (u-(*currentValueArray)[aI][aJ][aK-1][aL])/he;
  } else {
    dminuse = 0;
  }

  /* Calculate partial derivatives using upwind discretization */
  /* NOTE changing from min to max scheme! Did we do this right? */
  float u_x, u_y, u_e;
  if (dplusx <= dminusx && dplusx < 0) {
    u_x = -dplusx;
    assert(u_x > 0);
  } else if (dminusx <= dplusx && dminusx < 0) {
    u_x = dminusx;
    assert(u_x < 0);
  } else {
    u_x = 0;
  }
  if (dplusy <= dminusy && dplusy < 0) {
    u_y = -dplusy;
    assert(u_y > 0);
  } else if (dminusy <= dplusy && dminusy < 0) {
    u_y = dminusy;
    assert(u_y < 0);
  } else {
    u_y = 0;
  }

  // need these parameters for energy upwind
  const float K = currentGrid->getCost(aI,aJ,aK,aL);
  const float psi = aFoodAccumulated; // we can do this because we pass in zero for mode 2! ...

  /* NOTE changed conditionals... */
  /* ADD does this not change for max vs min?... */
  if (psi - K > 0) {
    u_e = -dpluse;
  } else if (psi - K < 0) {
    u_e = dminuse;
  } else {
    u_e = 0;
  }

  const float g = currentGrid->getSpeed(aI,aJ,aK,aL);
  const float dt = currentGrid->getDt();

  const float u1 = (*aValue1)[aI][aJ][aK][aL];
  const float u2 = (*aValue2)[aI][aJ][aK][aL];

  // retrieve these parameters from the grids that they correspond to
  const float phi = fTerrain->getSpottingRate(aI,aJ,aK,aL);
  const float mu_k = fTerrain->getKillRate(aI,aJ,aK,aL);
  const float mu_g = fTerrain->getGiveUpRate(aI,aJ,aK,aL);

  float val = 0; // updated value to be returned

  if (aMode == 1) {
    val = u - phi*u*dt + phi*u2*dt + (psi-K)*u_e*dt
          + g*sqrt(u_x*u_x + u_y*u_y)*dt;   // mode 1 HJB update
  } else if (aMode == 2) {
    val = u + mu_g*u1*dt - (mu_g+mu_k)*u*dt - K*u_e*dt
          + g*sqrt(u_x*u_x + u_y*u_y)*dt; // mode 2 HJB update
  }

  return val;
}

/*==============================================================================
  Read/write to file
==============================================================================*/
void CTimeDependentHjbSolver::writeResults(const string aFilename, 
                                           const int aCurrentStage) const {
  // writing BigPath, Modes, Energy to file
  if ((aCurrentStage == fNStages - 1) || !(fBigModeList.back()==1 || fBigModeList.back()==2)) {
    // Either we made it to the final stage or we died early
    // Either way, write the end of simulation results to file
    const int path_length = fBigPath.size();
    ofstream out("output/" + aFilename + "_OptimalTrajectory", ios::binary);
    for (int i = 0; i < path_length; ++i) {
      out.write((char*) &fBigPath[i].x, sizeof(double));
      out.write((char*) &fBigPath[i].y, sizeof(double));
    }

    ofstream outSteps("output/" + aFilename + "_TotalSteps", ios::binary);
    outSteps.write((char*) &path_length, sizeof(int));

    io::writeVectorToFile<double>(aFilename + "_EnergyList", fBigEnergyList);

    io::writeVectorToFile<int>(aFilename + "_ModeList", fBigModeList);

    // writing total number of stages to file
    ofstream outNStages("output/" + aFilename + "_NStages", ios::binary);
    outNStages.write((char*) &fNStages, sizeof(int));

    // write number of survived stages
    ofstream outNSurvivedStages("output/" + aFilename + "_NSurvivedStages", ios::binary);
    outNSurvivedStages.write((char*) &fNSurvivedStages, sizeof(int));

    // writing total number of stages to file
    ofstream outTime("output/" + aFilename + "_MaxTime", ios::binary);
    outTime.write((char*) &fTimeMax, sizeof(double));

    ofstream outEnergyThreshold("output/" + aFilename + "_EnergyThreshold", ios::binary);
    outEnergyThreshold.write((char*) &fEnergyThreshold, sizeof(int));

    fTerrain->writeTerrainToFile(aFilename);
    fFoodDensity->writeFoodParamsToFile(aFilename);
  }

  // Write stage dependent output to file
  fFoodDensity->writeFoodGridToFile(aFilename + "_Stage" + to_string(aCurrentStage));
  fGrid1->writeGridToFile(aFilename + "_Stage" + to_string(aCurrentStage) + "_Mode1");
  fGrid2->writeGridToFile(aFilename + "_Stage" + to_string(aCurrentStage) + "_Mode2");

}

void CTimeDependentHjbSolver::readValuesFromFile(const string aFilename) {
  if(fFullGrid) {
    fGrid1->readValuesFromFile(aFilename + "_Stage0_Mode1Value");
    fGrid2->readValuesFromFile(aFilename + "_Stage0_Mode2Value");
  } else {
    /* Read in first slice */
    fGrid1->setCurrentSlice(0);
    fGrid2->setCurrentSlice(0);
    fGrid1->readValuesFromFile(aFilename + "_Stage0_Mode1Value", 0);
    fGrid2->readValuesFromFile(aFilename + "_Stage0_Mode2Value", 0);
    fGrid1->setFilename(aFilename + "_Stage0_Mode1Value");
    fGrid2->setFilename(aFilename + "_Stage0_Mode2Value");
  }
}
