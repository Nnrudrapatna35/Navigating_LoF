/*
* ==============================================================================
*
*  Copyright (C) 2024 Marissa Gee Nicolas Gonzalez Granda Sunay Joshi 
*                     Nagaprasad Rudrapatna Anne Somalwar 
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
* File: CTimeDependentHjbSolver.hpp
*
* Authors: Marissa Gee Nicolas Gonzalez Granda Sunay Joshi 
*          Nagaprasad Rudrapatna Anne Somalwar 
*          (based on code by Elliot Cartee, Qianli Song, Lexiao Lai)
*             (based on code by Marc Aurèle Gilles and Zachary Clawson)
*
* Description: This class is for solving the time-dependent HJB equation.
* It assumes that fGrid contains the speed and cost functions.
* It assumes a 2D regular grid with equal spacing in both directions.
* The target is assumed to be single points.
* Updates from both 5-point and 9-point stencils are possible,
*     depends on setting in GlobalConfigurations.hpp
* It depends on the CTimeGrid class.
*
* ==============================================================================
*/

#ifndef CTIMEDEPENDENTHJBSOLVER_HPP
#define CTIMEDEPENDENTHJBSOLVER_HPP

/** ------ Libraries ---------------------------------------------------------*/
#include <cmath>
#include <string>

/** ------ Project-specific header files -------------------------------------*/
#include "CTimeGrid.hpp"
#include "CFoodDensity.hpp"
#include "GlobalConfigurations.hpp"
#include "CTimeDependentTracer.hpp"
#include "CEikonalSolver.hpp"
#include "MemoryAllocations.hpp"
#include "CTerrain.hpp"

class CTimeDependentHjbSolver
{
  public:
    /** Default constructor */
    CTimeDependentHjbSolver() = default;

    /** Constructor
    * @param aGrid1  shared pointer to a 4D TimeGrid representing value function in mode 1
    * @param aGrid2  shared pointer to a 4D TimeGrid representing value function in mode 2
    */

    // function-based constructor
    CTimeDependentHjbSolver(std::function<double(double, double, double, double)> aCost1,
                            std::function<double(double, double, double, double)> aCost2,
                            std::function<double(double, double, double, double)> aSpeed1,
                            std::function<double(double, double, double, double)> aSpeed2,
                            std::shared_ptr<CFoodDensity> aFoodDensity,
                            std::shared_ptr<CTerrain> aTerrain,
                            const double aEnergyThreshold, 
                            const int aN, const int aNe,
                            const double aPhysMin, const double aPhysMax,
                            const double aEnergyMin, const double aEnergyMax,
                            const double aTimeMax,
                            const std::vector<double> aSource,
                            const double aInitialEnergy,
                            const objective_t aObjective, const int aNStages,
                            const bool aFullGrid = true, 
                            const bool aReturnHome = true, 
                            const bool aAllowModeSwitches = true);

    // alternate function-based constructor
    CTimeDependentHjbSolver(std::function<double(double, double, double, double)> aCost1,
                            std::function<double(double, double, double, double)> aCost2,
                            std::function<double(double, double, double, double)> aSpeed1,
                            std::function<double(double, double, double, double)> aSpeed2,
                            std::shared_ptr<CFoodDensity> aFoodDensity,
                            std::shared_ptr<CTerrain> aTerrain,
                            const double aEnergyThreshold, const int aNt,
                            const int aN, const int aNe,
                            const double aPhysMin, const double aPhysMax,
                            const double aEnergyMin, const double aEnergyMax,
                            const double aTimeMax,
                            const std::vector<double> aSource,
                            const double aInitialEnergy,
                            const objective_t aObjective, const int aNStages,
                            const bool aFullGrid = true, 
                            const bool aReturnHome = true, 
                            const bool aAllowModeSwitches = true);

    // array-based constructor
    CTimeDependentHjbSolver(std::shared_ptr<memory::array2D_t<double>> aCost1,
                            std::shared_ptr<memory::array2D_t<double>> aCost2,
                            std::shared_ptr<memory::array2D_t<double>> aSpeed1,
                            std::shared_ptr<memory::array2D_t<double>> aSpeed2,
                            std::shared_ptr<CFoodDensity> aFoodDensity,
                            std::shared_ptr<CTerrain> aTerrain,
                            const double aEnergyThreshold,
                            const int aN, const int aNe,
                            const double aPhysMin, const double aPhysMax,
                            const double aEnergyMin, const double aEnergyMax,
                            const double aTimeMax,
                            const std::vector<double> aSource,
                            const double aInitialEnergy,
                            const objective_t aObjective, const int aNStages,
                            const bool aFullGrid = true, 
                            const bool aReturnHome = true, 
                            const bool aAllowModeSwitches = true);

    /** This is the main function for solving the HJB equation */
    void solveHJB(const int aCurrStage, const std::string aFilename, 
                  const stencil_t aStencilType);

    // this traces optimal path + prints to file
    void traceHJB(const bool aIsFinalStage, const std::string aFilename,
                  const bool aReadFromFile = false,
                  const bool aClassifyPath = false);


    void traceNPaths(const int aNumPaths, const std::string aFilename, 
                     const bool aReadFromFile = false,
                     const int aInitialMode = 1);

    void traceNDifferentPaths(const int aNumPaths, const std::string aFilename, 
                              const bool aReadFromFile = false,
                              const int aInitialMode = 1);

    // this undoes food depletion, helpful for tracing many paths in a row
    void resetFoodDensity();

    // Set the initial energy for path tracing
    void setInitialEnergy(const double aInitialEnergy);
    void setInitialMode(const int aInitialMode);
    void setModeSwitches(const bool aAllowModeSwitches);
    void setSpottingTimes(std::vector<double> aSpottingTimes);

    // Get the initial energy for path tracing
    double getInitialEnergy() const;

    // this solves multistage
    void solveMultiStage(const std::string aFilename,
                         const stencil_t aStencilType = FIVE_POINT);

    // this solves receding horizon
    void solveRecedingHorizon(const std::string aFilename,
                              const stencil_t aStencilType = FIVE_POINT);

    /** Function for writing results to file */
    void writeResults(const std::string aFilename, const int aCurrentStage) const;

    void readValuesFromFile(const std::string aFilename);

    void setFilename(const std::string aFilename);


  private:
    /* Shared pointers to *two* TimeGrid objects which contain
     * - The approximation u[i,j,k] of the value function (for all time l) (for each mode)
     * - Parameters of the model
     */
    std::shared_ptr<CTimeGrid> fGrid1;
    std::shared_ptr<CTimeGrid> fGrid2;
    std::shared_ptr<CFoodDensity> fFoodDensity; // added FoodDensity object
    std::shared_ptr<CTerrain> fTerrain; // added FoodDensity object

    objective_t fObjective; //*ENERGY MAX*
    int fNStages;
    int fNSurvivedStages; // added
    double fRadius;
    double fEnergyThreshold;
    bool fFullGrid;
    bool fReturnHome;

    double fTimeMax; // added

    // pointer to tracer objects
    std::shared_ptr<CTimeDependentTracer> fTracer;

    // pointer to 2d array of return times (for terminal conditions)
    std::shared_ptr<memory::array2D_t<double>> fReturnTimes;

    // big path, modelist, energylist objects
    // NOT initialized in constructor?...
    Path fBigPath;
    std::vector<int> fBigModeList;
    std::vector<double> fBigEnergyList;


    /* Helper function that computes HJB update using 5-point stencil */
    double update_5pt(const int aI, const int aJ, const int aK, const int aL, 
                      const int aMode) const;

    float update_5pt_terminal(std::shared_ptr<memory::array4D_t<float>> aValue1,
                              std::shared_ptr<memory::array4D_t<float>> aValue2,
                              const int aI, const int aJ, const int aK, 
                              const int aL, const int aMode,
                              const float aFoodAccumulated) const;
                    
    /* Terminal condition functions */
    void vanillaTerminal();
    void goHomeTerminal();
    void expandingGoHomeTerminal(const int aCurrStage); // added
    void recedingHorizonTerminal(const int aCurrStage, const stencil_t aStencilType);

    /* Helper function that computes HJB update using 9-point stencil */
    double update_9pt(const int aI, const int aJ, const int aK, const int aL, 
                      const int aMode) const;

    /* ADD figure out how to change below... */
    /* Helper function for evaluating gradient on 9-point stencil */
    void computeGradient(const int i_diff, const int i_direct, const double v_diag,
                         const double v_dir, double (&u)[2]) const;
};

/* ============================================================================
*    Inline function definitions
* ============================================================================*/

inline void CTimeDependentHjbSolver::setInitialEnergy(const double aInitialEnergy) {
  fTracer->setInitialEnergy(aInitialEnergy);
}

inline void CTimeDependentHjbSolver::setInitialMode(const int aInitialMode) {
  fTracer->setInitialMode(aInitialMode);
}

inline void CTimeDependentHjbSolver::setModeSwitches(const bool aAllowModeSwitches) {
  fTracer->setModeSwitches(aAllowModeSwitches);
}

inline void CTimeDependentHjbSolver::setSpottingTimes(std::vector<double> aSpottingTimes) {
  fTracer->setSpottingTimes(aSpottingTimes);
}

inline void CTimeDependentHjbSolver::setFilename(const std::string aFilename) {
  fGrid1->setFilename(aFilename);
  fGrid2->setFilename(aFilename);
}

inline double CTimeDependentHjbSolver::getInitialEnergy() const {
  return fTracer->getInitialEnergy();
}

#endif
