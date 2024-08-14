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
* File: CTimeDependentTracer.cpp
*
* Authors: Marissa Gee Nicolas Gonzalez Granda Sunay Joshi 
*          Nagaprasad Rudrapatna Anne Somalwar 
*          (based on code by Elliot Cartee, Qianli Song, Lexiao Lai)
*             (based on code by Marc Aurèle Gilles and Zachary Clawson)
*
* Description: This is the trajectory tracer class for time-dependent problems.
* Given the solution to the time-dependent HJB equation, this class computes the
* optimal trajectory.
* It assumes a 2D regular grid with equal spacing in both directions.
* It uses trilinear interpolation of the HJB solution.
* The target and source are assumed to be single points.
* This implementation computes a grid search over a uniform grid of directions.
* It depends on the CTimeGrid class.
* (see also CTimeDependentTracer.hpp)
*
* ==============================================================================
*/

/** ------ Libraries ---------------------------------------------------------*/
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <random>
#include "omp.h"

/** ------ Project-specific header files -------------------------------------*/
#include "CTimeDependentTracer.hpp"
#include "MemoryAllocations.hpp"

/** ------ Namespaces --------------------------------------------------------*/
using namespace std;
using namespace memory;

/*==============================================================================
  Constructor
==============================================================================*/

// Finite horizon constructor
CTimeDependentTracer::CTimeDependentTracer(shared_ptr<CTimeGrid> aGrid1,
                                           shared_ptr<CTimeGrid> aGrid2,
                                           shared_ptr<CFoodDensity> aFoodDensity,
                                           shared_ptr<CTerrain> aTerrain,
                                           const vector<double> aSource,
                                           const int aInitialMode,
                                           const double aInitialEnergy, 
                                           const bool aAllowModeSwitches) {

  fGrid1 = aGrid1;
  fGrid2 = aGrid2;
  fFoodDensity = aFoodDensity;
  fTerrain = aTerrain;
  fSource = aSource;
  fInitialMode = aInitialMode;
  fInitialEnergy = aInitialEnergy;
  fInitialTime = 0;
  fAllowModeSwitches = aAllowModeSwitches;
  fSpottingTimes = {};
  fDt = fGrid1->getDt();

  const double h = fGrid1->getH();
  const double alpha = fFoodDensity->getHarvestRate();
  const double r = fFoodDensity->getFoodRadius();
  const double BMax = fFoodDensity->getFoodKernel(0);
  const double BMin = fFoodDensity->getFoodKernel(r);
  // fTimeFactor = min(1.0, ceil(alpha * fDt/pow(h,2) * 1.0/2.0));
  fTimeFactor = max(100.0, ceil(alpha * fDt/pow(h,2) * BMax/BMin));
  cout << "Time factor is " << fTimeFactor << endl;
}

/*==============================================================================
  Compute optimal path
==============================================================================*/
void CTimeDependentTracer::tracePath(const bool aIsFinalStage,
                                     const string aFilename) {
  shared_ptr<array2D_t<double>> foodDensityPtr = fFoodDensity->getFoodDensityPtr();
  tracePath(fPath, fModeList, fEnergyList, foodDensityPtr, fSource[0], fSource[1],
            aIsFinalStage, aFilename);
}


void CTimeDependentTracer::tracePath(Path& aPath, vector<int>& aModeList,
                                     vector<double>& aEnergyList,
                                     shared_ptr<array2D_t<double>> aFoodDensity,
                                     const double aInitialX,
                                     const double aInitialY,
                                     const bool aIsFinalStage,
                                     const string aFilename) {
  // Reset fPath, fModeList, fEnergyList
  aPath.clear();
  aModeList.clear();
  aEnergyList.clear();

  /* Pseudo-timestep for path computation; same in modes 1, 2 */
  const double tau = fDt / fTimeFactor; /* Usually run with fTimeFactor == 1 */
  const int nt = fGrid1->getGridSizeT();
  const int maxSteps = (nt - 1) * fTimeFactor; // Maximum number of timesteps for finite horizon

  /* Add source point to path */
  double currentX = aInitialX;
  double currentY = aInitialY;
  aPath.push_back(Point2d(currentX, currentY));

  /* Set initial energy */
  double currentE = fInitialEnergy; // Current energy, initialized to initial energy
  aEnergyList.push_back(currentE);
  cout << "Initial Energy: " << currentE << endl;

  /* Set initial time */
  double currentT = fInitialTime;

  /* Set initial mode */
  int currentMode = fInitialMode; // Start in fInitialMode, not necessarily 1!
  aModeList.push_back(currentMode);

  int numSteps = 1;
  shared_ptr<CTimeGrid> currentGrid; // If we are solving finite horizon
  if (currentMode == 1) {
      currentGrid = fGrid1;
  } else if (currentMode == 2) {
      currentGrid = fGrid2;
  }

  /* Read in slices if not using full grid */
  int currentSlice = 0;
  if (!currentGrid->isFullGrid()) {
    currentGrid->setCurrentSlice(currentSlice);
    currentGrid->readValuesFromFile(currentSlice);
    currentGrid->readValuesFromFile(currentSlice + 1);
  }


  /* Create random number generator */
  random_device rand_dev;
  mt19937 randomGenerator; // Set the seed randomly
  if (fSpottingTimes.size() > 0) {
    randomGenerator.seed(123456);
  } else {
    randomGenerator.seed(rand_dev());
  }


  /** ======== Main Loop =====================================================*/
  /* Loop until (possibly infinite) horizon */
  /* Loop based on mode + positive energy */

  const double minE = currentGrid->getMinE();
  const double maxE = currentGrid->getMaxE();
  const double de = currentGrid->getE();
  const double minX = currentGrid->getMinX();
  const double maxX = currentGrid->getMaxX();
  const double minY = currentGrid->getMinY();
  const double maxY = currentGrid->getMaxY();
  const double maxT = currentGrid->getMaxT();

  vector<double> realizedFood;
  vector<double> predictedFood;

  bool aboveThreshold = true;

  while (numSteps < maxSteps && currentE > 0) {
    /* Assume we stay in currentMode for tau seconds */    
    double speed, cost, muS, muK, muG;

    speed = currentGrid->getSpeedPhysical(currentX, currentY, currentE, currentT);
    cost = currentGrid->getCostPhysical(currentX, currentY, currentE, currentT);
    muS = fTerrain->getSpottingRatePhysical(currentX, currentY);
    muK = fTerrain->getKillRatePhysical(currentX, currentY);
    muG = fTerrain->getGiveUpRatePhysical(currentX, currentY);

    // Update newE based on whether currentMode = 1,2
    double newX, newY, newE;
    const double newT = currentT + tau;
    if (newT > maxT) {
      cout << setprecision(24) << "currentT " << currentT << endl;
      cout << "newT " << newT << endl;
      cout << "tau " << tau << endl;
      cout << "numSteps " << numSteps << endl;
      cout << "maxSteps " << maxSteps << endl;
    }

    // Read in next slice if needed
    if (!currentGrid->isFullGrid()) {
      const int newSlice = floor(newT/fDt);
      if ((newSlice != currentSlice) && (newSlice < nt - 1)) {
        currentGrid->setCurrentSlice(newSlice);
        currentGrid->readValuesFromFile(newSlice);
        currentGrid->readValuesFromFile(newSlice + 1);
        currentSlice = newSlice;
      } else if (newSlice == nt - 1) {
        currentGrid->setCurrentSlice(currentSlice);
        currentGrid->readValuesFromFile(currentSlice);
        currentGrid->readValuesFromFile(currentSlice + 1);
        currentSlice = newSlice;
      }
    } 

    if (currentMode == 1) {
      if (fFoodDensity->getFoodDepletion()) {
        
        // Check how much space we have in our stomach
        const double foodSpace = maxE - currentE;
        double foodFactor = 1;
        if (foodSpace < de) {
          // We might fill up before we finish eating
          const double possibleFood = tau * fFoodDensity->calcFoodAccumulated(currentX, currentY);
          if (possibleFood > foodSpace + tau*cost) {
            // If we will run out of space, scale food eaten to space available 
            foodFactor = (foodSpace + tau*cost)/possibleFood;
          }
        }

        // For the sake of comparison, we may print depleteFood versus tau*calcFoodAccumulated (units of energy)
        // cout << "calcFoodAcc value: " << tau * fFoodDensity->calcFoodAccumulated(currentX, currentY) << endl;
        double totalFood = fFoodDensity->depleteFood(currentX, currentY, 
                                                     tau, foodFactor,
                                                     aFoodDensity);
        // cout << "depleteFood value: " << totalFood << endl;
        newE = currentE + totalFood - tau * cost;
      } else {
        double totalFood = fFoodDensity->calcFoodAccumulated(currentX, currentY);
        newE = currentE + tau * (totalFood - cost);
      }
    } else if (currentMode == 2) {
      newE = currentE - tau*cost;
    } else {
      newE = 0;
    }
    if (newE >= maxE) {
      assert(currentE >= maxE - de);
    }
    newE = max(minE, min(maxE, newE)); // Round to within bounds

    /* Note that we're MAXIMIZING, so initialize maxValue to -infty */
    const double stayInPlaceValue = currentGrid->getValuePhysical(currentX, 
                                                                  currentY, 
                                                                  newE, newT);
    double maxValue = stayInPlaceValue;
    double bestTheta = -1;

    /* Grid search over directions */
    for (int i = 0; i < N_THETA; ++i) {
      const double theta = 2 * PI * (i) / N_THETA;
      newX = currentX + cos(theta) * tau * speed;
      newY = currentY + sin(theta) * tau * speed;

      // Find tempValue
      double tempValue;

      // If within bounds, use getValue
      if (minX <= newX && newX <= maxX && minY <= newY && newY <= maxY) {
        tempValue = currentGrid->getValuePhysical(newX, newY, newE, newT); 
      } else {
        // animal went out of bounds
        tempValue = NEG_INF;
      }

      /* Update if bigger, since we are maximizing */
      if ((tempValue > maxValue) && (tempValue > stayInPlaceValue + EPSILON)) {
        bestTheta = theta;
        maxValue = tempValue;
      }
    }

    // Check maxValue in the LAST step of the LAST stage
    if (numSteps == maxSteps-1) {
      cout << "MAX VAL IN LAST STEP: " << maxValue << endl;
    }

    // Find coordinates of the best step, other than the stay-in-place control
    // (Not a const, because we might change it to stay-in-place)
    double bestX = currentX + cos(bestTheta) * tau * speed;
    double bestY = currentY + sin(bestTheta) * tau * speed;

    // Check if stay-in-place won
    if (bestTheta == -1) {
      bestX = currentX;
      bestY = currentY;
    }

    // Compute next mode
    int newMode = currentMode;
    if (fSpottingTimes.size() > 0) {
      if (currentMode == 1) {
        newMode = 1;
        for (int n = 0; n < fSpottingTimes.size(); ++n) {
          if (fabs(currentT - fSpottingTimes[n]) < tau) {
            newMode = 2;
            cout << "Spotted at " << currentT << endl;
          }
        }
      } else if (currentMode == 2){
        exponential_distribution<double> distributionGaveUp(muG);
        const double gaveUpTime = distributionGaveUp(randomGenerator);

        // If nothing happened in tau, stay in mode 2
        if (gaveUpTime > tau) {
          newMode = 2;
        } else {
          newMode = 1;
        }
      }
    } else {
      if ((currentMode == 1) && (fAllowModeSwitches)) {
        exponential_distribution<double> distributionSpotted(muS);
        const double spottedTime = distributionSpotted(randomGenerator);

        // If nothing happened in tau, stay in mode 1
        if (spottedTime > tau) {
          newMode = 1;
        } else {
          newMode = 2;
        }
      } else if ((currentMode == 2) && (fAllowModeSwitches)) {
        exponential_distribution<double> distributionKilled(muK);
        exponential_distribution<double> distributionGaveUp(muG);

        double killedTime = distributionKilled(randomGenerator);
        const double gaveUpTime = distributionGaveUp(randomGenerator);

        double minTime = min(killedTime, gaveUpTime);

        // If nothing happened in tau, stay in mode 2
        if (minTime > tau) {
          newMode = 2;
        } else if (minTime == killedTime) {
          newMode = 3;
        } else if (minTime == gaveUpTime) {
          newMode = 1;
        }
      }
    }

    // If you are killed or died, you land at 0 energy
    if (newMode == 3) {
      newE = 0;
    }
    if (newE == 0 && newMode != 3) {
      // no energy, but you weren't killed/died, then you're exhausted; set mode 4
      newMode = 4;
    }

    if (newMode == 1) {
      currentGrid = fGrid1;
    } else if (newMode == 2) {
      currentGrid = fGrid2;
    }

    // Update currentE, currentMode
    currentX = bestX;
    currentY = bestY;
    currentE = newE;
    currentT = newT;
    currentMode = newMode;
    
    // Update path
    aPath.push_back(Point2d(bestX, bestY));
    aModeList.push_back(newMode); // Used to be currentMode; this is the same thing, but less clever
    aEnergyList.push_back(newE);

    numSteps++;
  }

  if (aModeList.back() == 1 || aModeList.back() == 2) {
    cout << "Survived this stage!" << endl;
  } else if (aModeList.back() == 3) {
    cout << "Killed by predator :(" << endl;
  } else if (aModeList.back() == 4) {
    cout << "Died of exhaustion :(" << endl;
  }
}

void CTimeDependentTracer::tracePathsParallel(const string aFilename,
                                              const int aNumPaths,
                                              const int aInitialMode) {
  const bool isFinalStage = true;
  const int nPhysical = fGrid1->getGridSizeX();
  setInitialMode(aInitialMode);

  if (aNumPaths == 1) {
    tracePath(isFinalStage, aFilename);
    printPathToFile(aFilename);
  } else {
    #pragma omp parallel for shared(isFinalStage, nPhysical)
    for (int k = 0; k < aNumPaths; ++k) {
      Path path;
      vector<int> modeList;
      vector<double> energyList;
      shared_ptr<array2D_t<double>> foodDensity = make_shared<array2D_t<double>>(
        allocateArray2D<double>(nPhysical, nPhysical));

      // Set food density based on undepleted values
      for (int i = 0; i < nPhysical; ++i) {
        for (int j = 0; j < nPhysical; ++j) {
          (*foodDensity)[i][j] = fFoodDensity->getInitialFoodDensity(i,j);
        }
      }

      tracePath(path, modeList, energyList, foodDensity, fSource[0], fSource[1],
                isFinalStage, aFilename);
      printPathToFile(aFilename + "_Path_" + to_string(k+1), path, modeList, energyList);
    }
  }
}

void CTimeDependentTracer::traceDifferentPathsParallel(const string aFilename,
                                                       const int aNumPaths,
                                                       const int aInitialMode) {
  const bool isFinalStage = true;
  const bool classifyEachPath = false; 
  const int nPhysical = fGrid1->getGridSizeX();
  const double maxX = fGrid1->getMaxX();
  const double maxY = fGrid1->getMaxY();

  setInitialMode(aInitialMode);

  std::default_random_engine generator;
  std::uniform_real_distribution<double> xDistribution(0.0,1.0);
  std::uniform_real_distribution<double> yDistribution(0.0,1.0);

  #pragma omp parallel for shared(isFinalStage, nPhysical)
  for (int k = 0; k < aNumPaths; ++k) {
    Path path;
    vector<int> modeList;
    vector<double> energyList;
    shared_ptr<array2D_t<double>> foodDensity = make_shared<array2D_t<double>>(
      allocateArray2D<double>(nPhysical, nPhysical));

    // Set food density based on undepleted values
    for (int i = 0; i < nPhysical; ++i) {
      for (int j = 0; j < nPhysical; ++j) {
        (*foodDensity)[i][j] = fFoodDensity->getInitialFoodDensity(i,j);
      }
    }

    // Set random starting location
    double initialX = xDistribution(generator);
    double initialY = yDistribution(generator);
    
    bool obstacle = fTerrain->isObstaclePhysicalConservative(initialX, initialY);
    bool sleepingSite;
    if (fInitialMode == 1) {
      sleepingSite = fTerrain->isSleepingSitePhysicalConservative(initialX, initialY);
    } else {
      sleepingSite = true;
    }
    while(obstacle == true || sleepingSite == false) {
      initialX = xDistribution(generator);
      initialY = yDistribution(generator);

      obstacle = fTerrain->isObstaclePhysicalConservative(initialX, initialY);
      if (fInitialMode == 1) {
        sleepingSite = fTerrain->isSleepingSitePhysicalConservative(initialX, initialY);
      } 
    }

    tracePath(path, modeList, energyList, foodDensity, initialX, initialY, 
              isFinalStage, aFilename);
    printPathToFile(aFilename + "_Path_" + to_string(k+1), path, modeList, energyList);
  }
}


/*==============================================================================
  See which patch we go to
==============================================================================*/
patch_t CTimeDependentTracer::classifyPath(const string aFilename) {
  bool visitedTopPatch = false;
  bool visitedBottomPatch = false;

  for (int i = 0; i < fPath.size(); ++i) {
    const double x = fPath[i].x;
    const double y = fPath[i].y;

    if (!(fTerrain->isHomeBasePhysical(x,y))) {
      if (y > 0.5) {
        visitedTopPatch = true;
      } else if (y < 0.5) {
        visitedBottomPatch = true;
      }
    }
  }

  patch_t pathType = HOME;
  if (visitedBottomPatch && visitedTopPatch) {
    pathType = BOTH;
  } else if (visitedBottomPatch) {
    pathType = LOWER;
  } else if (visitedTopPatch) {
    pathType = UPPER;
  }

  ofstream outPathType("output/" + aFilename + "_PathType", ios::binary);
  outPathType.write((char*) &pathType, sizeof(int));

  return pathType;
}


/*==============================================================================
  Compute optimal path and write to file
==============================================================================*/
void CTimeDependentTracer::printPathToFile(const string aFilename) {
  printPathToFile(aFilename, fPath, fModeList, fEnergyList);
}

void CTimeDependentTracer::printPathToFile(const string aFilename, Path& aPath,
                                           vector<int>& aModeList,
                                           vector<double>& aEnergyList) {

  // Write path to file
  ofstream out("output/" + aFilename + "_Path", ios::binary);
  for (int i = 0; i < aPath.size(); ++i) {
    out.write((char*) &aPath[i].x, sizeof(double));
    out.write((char*) &aPath[i].y, sizeof(double));
  }

  // Write modes
  io::writeVectorToFile<int>(aFilename + "_Modes", aModeList);

  // Write energies
  io::writeVectorToFile<double>(aFilename + "_Energy", aEnergyList);

  // Write number of steps to a file
  int path_length = aPath.size();

  cout << "Steps taken: ";
  cout << path_length << endl;


  ofstream outSteps("output/" + aFilename + "_Steps", ios::binary);
  outSteps.write((char*) &path_length, sizeof(int));

  ofstream outTimeFactor("output/" + aFilename + "_TimeFactor", ios::binary);
  outTimeFactor.write((char*) &fTimeFactor, sizeof(int));
}
