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
* File: CTerrain.hpp
*
* Authors: Marissa Gee Nicolas Gonzalez Granda Sunay Joshi 
*          Nagaprasad Rudrapatna Anne Somalwar 
*
* Description: This class defines a Terrain object.
*
* ==============================================================================
*/

#ifndef CTERRAIN_HPP
#define CTERRAIN_HPP

/* ----- Libraries ----------------------------------------------------------*/
#include <cmath>
#include <string>
#include "boost/multi_array.hpp"

/* ------ Project-specific header files -------------------------------------*/
#include "MemoryAllocations.hpp"
#include "SimpleFunctions.hpp"

class CTerrain
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

    /* Food density, stored as pointer to a 2d array (x,y) */
    std::shared_ptr<memory::array2D_t<double>> fSpottingRate;
    std::shared_ptr<memory::array2D_t<double>> fKillRate;
    std::shared_ptr<memory::array2D_t<double>> fGiveUpRate;


    std::shared_ptr<memory::array2D_t<double>> fHomeBase;
    std::shared_ptr<memory::array2D_t<double>> fObstacle;
    std::shared_ptr<memory::array2D_t<double>> fSleepingSite;

  public:
    /* ========================================================================
    *    Constructors
    * ========================================================================*/
    /* Default constructor */
    CTerrain() = default;

    /*
     * This is the main CFoodDensity constructor
     * It initializes the grid and computes the grid spacing
     * for a square grid with aN points in both spatial dimensions.
     * @param aFoodDensity either a function with two real inputs or a shared pointer to a 2d array of real numbers. The food density.
     * @param aFoodKernel a function with two real inputs. The food kernel function.
     * @param aHomeBase either a function with two real inputs or a shared pointer to a 2d array of real numbers. The home base.
     * @param aObstacle either a function with two real inputs or a shared pointer to a 2d array of real numbers. The obstacle.
     * @param aRadius a real number. The radius of the kernel function inside which food density is updated.
     * @param aN an integer. The number of spatial grid points. Defaults to 101.
     * @param aPhysMin a real number. The lower bound on the physical grid coordinates. Defaults to 0.
     * @param aPhysMax a real number. The upper bound on the physical grid coordinates. Defaults to 1.
     * The domain is
     * [aPhysMin,aPhysMax] x [aPhysMin,aPhysMax].
     */

    /* Function-based constructor */
    CTerrain(std::function<double(double, double, double, double)> aSpottingRate,
             std::function<double(double, double, double, double)> aKillRate,
             std::function<double(double, double, double, double)> aGiveUpRate,
             std::function<double(double, double)> aHomeBase,
             std::function<double(double, double)> aObstacle,
             const int aN = 101, const double aPhysMin = 0,
             const double aPhysMax = 1);

    /* Array-based constructor (for aFoodDensity, aHomeBase) */
    CTerrain(std::shared_ptr<memory::array2D_t<double>> aSpottingRate,
             std::shared_ptr<memory::array2D_t<double>> aKillRate,
             std::shared_ptr<memory::array2D_t<double>> aGiveUpRate,
             std::shared_ptr<memory::array2D_t<double>> aHomeBase,
             std::shared_ptr<memory::array2D_t<double>> aObstacle,
             const int aN = 101, const double aPhysMin = 0,
             const double aPhysMax = 1, 
             std::shared_ptr<memory::array2D_t<double>> aSleepingSite = nullptr);

    void writeTerrainToFile(const std::string aFilename) const;
    /* ========================================================================
    *    Setters
    *=========================================================================*/

    /* ========================================================================
    *    Getters
    * ========================================================================*/
    /* Grid-coordinate function */
    double getHomeBase(const int aI, const int aJ) const;
    double getSpottingRate(const int aI, const int aJ, const int aK = 0, const int aL = 0) const;
    double getKillRate(const int aI, const int aJ, const int aK = 0, const int aL = 0) const;
    double getGiveUpRate(const int aI, const int aJ, const int aK = 0, const int aL = 0) const;

    bool isHomeBase(const int aI, const int aJ) const;
    bool isObstacle(const int aI, const int aJ) const;
    bool isSleepingSite(const int aI, const int aJ) const;
    bool isRegularTerrain(const int aI, const int aJ) const;

    /* Physical-coordinate function */
    double getSpottingRatePhysical(const double aX, const double aY, const double aE = 0, const double aT = 0) const;
    double getKillRatePhysical(const double aX, const double aY, const double aE = 0, const double aT = 0) const;
    double getGiveUpRatePhysical(const double aX, const double aY, const double aE = 0, const double aT = 0) const;

    bool isHomeBasePhysical(const double aX, const double aY) const;
    bool isObstaclePhysical(const double aX, const double aY) const;
    bool isRegularTerrainPhysical(const double aX, const double aY) const;
    bool isObstaclePhysicalConservative(const double aX, const double aY) const;
    bool isSleepingSitePhysicalConservative(const double aX, const double aY) const;

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
};
/* ============================================================================
*    Inline function definitions
* ============================================================================*/
inline double CTerrain::getHomeBase(const int aI, const int aJ) const {
  return (*fHomeBase)[aI][aJ];
}

inline double CTerrain::getSpottingRate(const int aI, const int aJ, const int aK, const int aL) const {
  return (*fSpottingRate)[aI][aJ]*getHomeBase(aI, aJ);
}

inline double CTerrain::getKillRate(const int aI, const int aJ, const int aK, const int aL) const {
  return (*fKillRate)[aI][aJ]*getHomeBase(aI, aJ);
}

inline double CTerrain::getGiveUpRate(const int aI, const int aJ, const int aK, const int aL) const {
  double mug = (*fGiveUpRate)[aI][aJ];
  double scalingFactor = 10;
  if (mug == 0) {
    mug = mug + scalingFactor*(1-getHomeBase(aI, aJ));
  } else {
    mug = mug*(scalingFactor - (scalingFactor - 1 )*getHomeBase(aI, aJ));
  }
  return mug;
}

inline bool CTerrain::isHomeBase(const int aI, const int aJ) const {
  bool value = false;
  if ((*fHomeBase)[aI][aJ] < 1) {
    value = true;
  }
  return value;
}

inline bool CTerrain::isObstacle(const int aI, const int aJ) const {
  bool value = false;
  if ((*fObstacle)[aI][aJ] == 0) {
    value = true;
  }
  return value;
}

inline bool CTerrain::isSleepingSite(const int aI, const int aJ) const {
  bool value = false;
  if ((*fSleepingSite)[aI][aJ] == 0) {
    value = true;
  }
  return value;
}

inline bool CTerrain::isRegularTerrain(const int aI, const int aJ) const {
  bool value = true;
  if (isObstacle(aI, aJ) || isHomeBase(aI, aJ)) {
    value = false;
  }
  return value;
}

inline bool CTerrain::isHomeBasePhysical(const double aX, const double aY) const {
  const int i = xPhysicalToGrid(aX);
  const int j = yPhysicalToGrid(aY);
  bool value = false;

  if ((*fHomeBase)[i][j] < 1) {
    value = true;
  }
  return value;
}

inline bool CTerrain::isObstaclePhysical(const double aX, const double aY) const {
  const int i = xPhysicalToGrid(aX);
  const int j = yPhysicalToGrid(aY);
  bool value = false;

  if ((*fObstacle)[i][j] == 0) {
    value = true;
  }
  return value;
}

inline bool CTerrain::isRegularTerrainPhysical(const double aX, const double aY) const {
  bool value = true;

  if (isObstaclePhysical(aX, aY) || isHomeBasePhysical(aX, aY)) {
    value = false;
  }
  return value;
}

inline int CTerrain::getGridSizeX() const{
  return fHomeBase->shape()[0];
}
inline int CTerrain::getGridSizeY() const{
  return fHomeBase->shape()[1];
}
inline double CTerrain::getH() const{
  return fH;
}
inline double CTerrain::getMinX() const{
  return fMinX;
}
inline double CTerrain::getMinY() const{
  return fMinY;
}
inline double CTerrain::getMaxX() const{
  return fMaxX;
}
inline double CTerrain::getMaxY() const{
  return fMaxY;
}

/* ------ Inline definition of conversion functions -------------------------*/
inline double CTerrain::xGridToPhysical(const int aI) const {
  return fMinX + (double)aI * fH;
}

inline double CTerrain::yGridToPhysical(const int aJ) const {
  return fMinY + (double)aJ * fH;
}

inline int CTerrain::xPhysicalToGrid(const double aX) const {
  return std::round((aX - fMinX)/fH);
}

inline int CTerrain::yPhysicalToGrid(const double aY) const {
  return std::round((aY - fMinY)/fH);
}

#endif
