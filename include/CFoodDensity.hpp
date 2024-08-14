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
* File: CFoodDensity.hpp
*
* Authors: Marissa Gee Nicolas Gonzalez Granda Sunay Joshi 
*          Nagaprasad Rudrapatna Anne Somalwar 
*
* Description: This class defines a FoodDensity object.
*
* ==============================================================================
*/

#ifndef CFOODDENSITY_HPP
#define CFOODDENSITY_HPP

/* ----- Libraries ----------------------------------------------------------*/
#include <cmath>
#include <string>
#include <iostream>
#include "boost/multi_array.hpp"

/* ------ Project-specific header files -------------------------------------*/
#include "MemoryAllocations.hpp"
#include "SimpleFunctions.hpp"
#include "CTerrain.hpp"

class CFoodDensity
{
  private:
    /* Properties of FoodDensity object */
    /* Grid spacing: assuming to be the same in the two spatial variables */
    double fH;

    /* Physical bounds on the grid */
    double fMinX;
    double fMinY;
    double fMaxX;
    double fMaxY;

    /* Food density function
        This has to be given as an input.
        Computes food density at physical location (x,y,e,t).
        Signature is fFoodDensity(x,y,e,t).
    std::function<double(double, double, double, double)> fFoodDensity; */

    /* Food density, stored as pointer to a 2d array (x,y) */
    std::shared_ptr<memory::array2D_t<double>> fFoodDensity;

    /* Initial food density used in the constructor*/
    std::shared_ptr<memory::array2D_t<double>> fInitialFoodDensity;

    /* Food kernel function to be used in foodAccumulated function */
    std::function<double(double, double)> fFoodKernel;

    /* Utility function that determines value of energy */
    std::function<double(double)> fUtilityFunction;

    /* Radius parameter for fFoodKernel */
    double fRadius;

    /* Proportion of food available which is eaten */
    double fHarvestRate;

    /* Whether or not to allow the food to deplete */
    bool fFoodDepletion;

    std::shared_ptr<memory::array2D_t<double>> fFoodAccumulated;
    std::shared_ptr<memory::array2D_t<double>> fInitialFoodAccumulated; 

    std::shared_ptr<CTerrain> fTerrain;

  public:
    /* ========================================================================
    *    Constructors
    * ========================================================================*/
    /* Default constructor */
    CFoodDensity() = default;

    /*
     * This is the main CFoodDensity constructor
     * It initializes the grid and computes the grid spacing
     * for a square grid with aN points in both spatial dimensions.
     * @param aFoodDensity either a function with two real inputs or a shared pointer to a 2d array of real numbers. The food density.
     * @param aFoodKernel a function with two real inputs. The food kernel function.
     * @param aTerrain a class containing information about home base and obstacles
     * @param aRadius a real number. The radius of the kernel function inside which food density is updated.
     * @param aN an integer. The number of spatial grid points. Defaults to 101.
     * @param aPhysMin a real number. The lower bound on the physical grid coordinates. Defaults to 0.
     * @param aPhysMax a real number. The upper bound on the physical grid coordinates. Defaults to 1.
     * The domain is
     * [aPhysMin,aPhysMax] x [aPhysMin,aPhysMax].
     */

    /* Function-based constructor */
    CFoodDensity(std::function<double(double, double)> aFoodDensity,
                 std::function<double(double, double)> aFoodKernel,
                 std::function<double(double)> aUtilityFunction,
                 std::shared_ptr<CTerrain> aTerrain,
                 const double aRadius, const double aHarvestRate, 
                 const bool aFoodDepletion, const int aN = 101, 
                 const double aPhysMin = 0, const double aPhysMax = 1);

    /* Function-based constructor */
    CFoodDensity(std::function<double(double, double)> aFoodDensity,
                 std::function<double(double, double)> aFoodKernel,
                 std::function<double(double, double)> aFoodAccumulated,
                 std::function<double(double)> aUtilityFunction,
                 std::shared_ptr<CTerrain> aTerrain,
                 const double aRadius, const double aHarvestRate, 
                 const bool aFoodDepletion, const int aN = 101, 
                 const double aPhysMin = 0, const double aPhysMax = 1);

    /* Array-based constructor (for aFoodDensity) */
    CFoodDensity(std::shared_ptr<memory::array2D_t<double>> aFoodDensity,
                 std::function<double(double, double)> aFoodKernel,
                 std::function<double(double)> aUtilityFunction,
                 std::shared_ptr<CTerrain> aTerrain,
                 const double aRadius, const double aHarvestRate, 
                 const bool aFoodDepletion, const int aN = 101, 
                 const double aPhysMin = 0, const double aPhysMax = 1);


    void updateFoodAccumulated();

    /* ========================================================================
    *    Setters
    *=========================================================================*/

    /* Set fFoodDensity to the inputted food density */
    void setFoodDensity(std::shared_ptr<memory::array2D_t<double>> aFoodDensity);

    /* Reset fFoodDensity to the food density inputted to the constructor*/
    void resetFoodDensity();

    /* Set the food density array at position (aI,aJ) to aFoodValue */
    void setFoodDensityPoint(const int aI, const int aJ, const double aFoodValue);

    void setFoodDepletion(const bool aFoodDepletion);

    /* ========================================================================
    *    Getters
    * ========================================================================*/
    /* Grid-coordinate function */
    double getFoodDensity(const int aI, const int aJ) const;

    /* Physical-coordinate function */
    double getFoodDensityPhysical(const double aX, const double aY) const;

    /* Grid-coordinate function */
    double getInitialFoodDensity(const int aI, const int aJ) const;

    /* Get pointer to food density grid*/
    std::shared_ptr<memory::array2D_t<double>> getFoodDensityPtr() const;

    /* Calculates the value of the food kernel function for a given distance and radius */
    double getFoodKernel(const double aDist) const;

    /* Calculates the utility of a given energy level */
    double getUtility(const double aE) const;

    /* Return the food accumulated at a grid point */
    double getFoodAccumulated(const int aI, const int aJ) const;

    /* Uses bilinear interpolation to evaluate the food accumulated */
    double getFoodAccumulatedPhysical(const double aX, const double aY) const;

    /* Harvest rate */
    double getHarvestRate() const;

    /* Food Radius */
    double getFoodRadius() const;

    /* Food Depletion */
    bool getFoodDepletion() const;

    /* Grid parameters */
    int getGridSizeX() const;
    int getGridSizeY() const;
    double getH() const;
    double getMinX() const;
    double getMinY() const;
    double getMaxX() const;
    double getMaxY() const;

    /* ========================================================================
    *    Other
    * ========================================================================*/
    /* Mapping back and forth between grid and physical */
    double xGridToPhysical(const int aI) const;
    double yGridToPhysical(const int aJ) const;
    int xPhysicalToGrid(const double aX) const;
    int yPhysicalToGrid(const double aY) const;

    /*
     * A I/O member which prints the value and cost grids to file.
     * @param aFilename a string which contains the prefix to the names
     *      of the files to which the grids will be printed.
     * The (cost/value) grids will be printed to
     *      files called "aFilename"+(Cost/Value)
     */
     void writeFoodGridToFile(const std::string aFilename) const;
     void writeFoodParamsToFile(const std::string aFilename) const;

    /* Calculate amount of food accumulated at physical position (aX,aY) */
    double calcFoodAccumulated(const double aX, const double aY);

    /* Deplete the food density at (aX,aY) after eating at (aX,aY) for aTau units of time
     * aTau should be dt / TimeFactor when called in CTimeDependentTracer.cpp */
    double depleteFood(const double aX, const double aY, const double aTau, 
                       const double aFoodFactor);
    double depleteFood(const double aX, const double aY, 
                       const double aTau, const double aFoodFactor,
                       std::shared_ptr<memory::array2D_t<double>> aFoodDensityArray);
};
/* ============================================================================
*    Inline function definitions
* ============================================================================*/

/* ------ Inline definition of setter function ------------------------------*/
inline void CFoodDensity::setFoodDensityPoint(const int aI, const int aJ, 
                                              const double aFoodValue) {
  (*fFoodDensity)[aI][aJ] = aFoodValue;
}

/* ------ Inline definition of getters --------------------------------------*/
inline double CFoodDensity::getFoodDensity(const int aI, const int aJ) const {
  double psi = 0;
  if (!(fTerrain->isObstacle(aI,aJ))) {
    psi = (*fFoodDensity)[aI][aJ];
  }
  return psi;
}

inline void CFoodDensity::setFoodDepletion(const bool aFoodDepletion) {
  fFoodDepletion = aFoodDepletion;
}

inline double CFoodDensity::getInitialFoodDensity(const int aI, const int aJ) const {
  double psi = 0;
  if (!(fTerrain->isObstacle(aI,aJ))) {
    psi = (*fInitialFoodDensity)[aI][aJ];
  }
  return psi;
}

inline double CFoodDensity::getFoodDensityPhysical(const double aX, const double aY) const {
  std::cout << "Not implemented." << std::endl;
  assert(false);
  return 0;
}

inline std::shared_ptr<memory::array2D_t<double>> CFoodDensity::getFoodDensityPtr() const {
  return fFoodDensity;
}

inline double CFoodDensity::getFoodKernel(const double aDist) const {
  return fFoodKernel(aDist, fRadius+EPSILON);
}

inline double CFoodDensity::getUtility(const double aE) const {
    return fUtilityFunction(aE);
}

inline double CFoodDensity::getFoodAccumulated(const int aI, const int aJ) const {
  return (*fFoodAccumulated)[aI][aJ];
}

inline double CFoodDensity::getHarvestRate() const{
  return fHarvestRate;
}

inline double CFoodDensity::getFoodRadius() const{
  return fRadius;
}

inline bool CFoodDensity::getFoodDepletion() const{
  return fFoodDepletion;
}

inline int CFoodDensity::getGridSizeX() const{
  return fFoodDensity->shape()[0];
}
inline int CFoodDensity::getGridSizeY() const{
  return fFoodDensity->shape()[1];
}
inline double CFoodDensity::getH() const{
  return fH;
}
inline double CFoodDensity::getMinX() const{
  return fMinX;
}
inline double CFoodDensity::getMinY() const{
  return fMinY;
}
inline double CFoodDensity::getMaxX() const{
  return fMaxX;
}
inline double CFoodDensity::getMaxY() const{
  return fMaxY;
}

/* ------ Inline definition of conversion functions -------------------------*/
inline double CFoodDensity::xGridToPhysical(const int aI) const {
  const double intermediate = aI*fH;
  // std::cout << "X intermediate: " << intermediate << std::endl;
  return fMinX + intermediate;
}

inline double CFoodDensity::yGridToPhysical(const int aJ) const {
  const double intermediate = aJ*fH;
  // std::cout << "Y intermediate: " << intermediate << std::endl;
  return fMinY + intermediate;
}

inline int CFoodDensity::xPhysicalToGrid(const double aX) const {
  return std::round((aX - fMinX)/fH);
}

inline int CFoodDensity::yPhysicalToGrid(const double aY) const {
  return std::round((aY - fMinY)/fH);
}

#endif
