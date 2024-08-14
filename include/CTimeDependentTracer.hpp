/*
* ==============================================================================
*
*  Copyright (C) 2019  Marc Aurèle Gilles
*  Copyright (C) 2019  Elliot Cartee, Qianli Song, Lexiao Lai
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
* File: CTimeDependentTracer.hpp
*
* Authors: REU 2022 Landscape of Fear Group
*   (based on code by Elliot Cartee, Qianli Song, Lexiao Lai)
*     (based on code by Marc Aurèle Gilles and Zachary Clawson)
*
* Description: This is the trajectory tracer class for time-dependent problems.
* Given the solution to the time-dependent HJB equation, this class computes the
* optimal trajectory.
* It assumes a 2D regular grid with equal spacing in both directions.
* It uses trilinear interpolation of the HJB solution.
* The target and source are assumed to be single points.
* This implementation computes a grid search over a uniform grid of directions.
* It depends on the CTimeGrid class.
*
* ==============================================================================
*/

#ifndef TIMEDEPENDENTTRACER_HPP
#define TIMEDEPENDENTTRACER_HPP

/** ------ Libraries ---------------------------------------------------------*/
#include <cmath>
#include <vector>
#include <string>

/** ------ Project-specific header files -------------------------------------*/
#include "GlobalConfigurations.hpp"
#include "CTimeGrid.hpp"
#include "CFoodDensity.hpp"
#include "WritetoFile.hpp"

/** ------ Data structures and typedefs --------------------------------------*/
/* Struct for a single Gridpoint */
struct Point2d {
  double x, y;
  Point2d(double a, double b) {
    x = a;
    y = b;
  }
};
/* Typedef for path objects */
typedef std::vector<Point2d> Path;

/** ------ Main Class definition ---------------------------------------------*/
class CTimeDependentTracer {
  public:
    /* ------ Constructors ---------------------------------------------------*/
    CTimeDependentTracer() = default;

    /**
     * The CTimeDependentTracer Constructor
     * This is the main constructor.
     * @param aGrid1 CTimeGrid object containing the value function of the PDE in mode 1.
     * @param aGrid2 CTimeGrid object containing the value function of the PDE in mode 2.
     * @param aFoodDensity CFoodDensity object.
     * @param aSource physical coordinates of the source.
     * @param aInitialEnergy initial energy of animal.
     */

    // constructor for finite horizon problem
    CTimeDependentTracer(std::shared_ptr<CTimeGrid> aGrid1,
                         std::shared_ptr<CTimeGrid> aGrid2,
                         std::shared_ptr<CFoodDensity> aFoodDensity,
                         std::shared_ptr<CTerrain> aTerrain,
                         const std::vector<double> aSource,
                         const int aInitialMode,
                         const double aInitialEnergy,
                         const bool aAllowModeSwitches = true);

    /*-- Main ----------------------------------------------------------------*/
    /**
    * A path printing function.
    * Computes the optimal path and prints it to file specified by aFilename
    * @param aFilename a string containing the name of the file.
    */
    void printPathToFile(const std::string aFilename);

    /**
    * The helper function for actually computing the optimal path
    */
    void tracePath(const bool aIsFinalStage, 
                   const std::string aFilename = "");

    void tracePathsParallel(const std::string aFilename, const int aNumPaths,
                            const int aInitialMode = 1);

    void traceDifferentPathsParallel(const std::string aFilename, 
                                     const int aNumPaths,
                                     const int aInitialMode = 1);

    // setters
    void setSource(double aX, double aY);
    void setInitialMode(int aInitialMode);
    void setInitialEnergy(double aInitialEnergy);
    void setInitialTime(double aInitialTime);
    void setModeSwitches(bool aAllowModeSwitches);
    void setSpottingTimes(std::vector<double> aSpottingTimes);

    double getInitialEnergy() const;
    Path getPath() const;
    std::vector<int> getModeList() const;
    std::vector<double> getEnergyList() const;

    int tGridToPhysical(const double aTPhysical) const;

    // Classify which patch we visit in risk-reward example
    patch_t classifyPath(const std::string aFilename); 

  private:
    /* The number of tracer steps per PDE timestep.
       tau = dt/fTimeFactor */
    int fTimeFactor;

    // delta t
    double fDt;

    /* Source X and Y coordinates */
    std::vector<double> fSource;

    /* Specified Spotting Times */
    std::vector<double> fSpottingTimes;

    // initial mode
    int fInitialMode;

    /* Initial energy */
    double fInitialEnergy;

    /* An initial time, must be compatible with horizon length and step size */
    double fInitialTime;

    /* Whether we allow mode switches to occur when path tracing */
    bool fAllowModeSwitches;

    /* Timegrids containing value functions of the PDE if finite horizon */
    std::shared_ptr<CTimeGrid> fGrid1;
    std::shared_ptr<CTimeGrid> fGrid2;

    // FoodDensity object with current food density
    std::shared_ptr<CFoodDensity> fFoodDensity;

    std::shared_ptr<CTerrain> fTerrain;

    // path, modelist, energylist objects
    Path fPath;
    std::vector<int> fModeList;
    std::vector<double> fEnergyList;

    void tracePath(Path& aPath, std::vector<int>& aModeList,
                   std::vector<double>& aEnergyList,
                   std::shared_ptr<memory::array2D_t<double>> aFoodDensity,
                   const double aInitialX, const double aInitialY,
                   const bool aIsFinalStage,
                   const std::string aFilename = "");

    void printPathToFile(const std::string aFilename, Path& aPath,
                         std::vector<int>& aModeList,
                         std::vector<double>& aEnergyList);
};

// setters
inline void CTimeDependentTracer::setSource(double aX, double aY) {
  fSource = {aX, aY};
}

inline void CTimeDependentTracer::setInitialMode(int aInitialMode) {
  fInitialMode = aInitialMode;
}

inline void CTimeDependentTracer::setInitialEnergy(double aInitialEnergy) {
  fInitialEnergy = aInitialEnergy;
}

inline void CTimeDependentTracer::setInitialTime(double aInitialTime) {
  fInitialTime = aInitialTime;
}

inline void CTimeDependentTracer::setModeSwitches(bool aAllowModeSwitches) {
  fAllowModeSwitches = aAllowModeSwitches;
}

inline void CTimeDependentTracer::setSpottingTimes(std::vector<double> aSpottingTimes) {
  fSpottingTimes = aSpottingTimes;
}

// getters
inline double CTimeDependentTracer::getInitialEnergy() const {
  return fInitialEnergy;
}

inline Path CTimeDependentTracer::getPath() const {
  return fPath;
}

inline std::vector<int> CTimeDependentTracer::getModeList() const {
  return fModeList;
}

inline std::vector<double> CTimeDependentTracer::getEnergyList() const {
  return fEnergyList;
}

inline int CTimeDependentTracer::tGridToPhysical(const double aTPhysical) const {
  const double tau = fDt/fTimeFactor;
  const int l = floor(aTPhysical/tau);
  return l;
}

#endif
