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
* File: CTerrain.cpp
*
* Authors: Marissa Gee Nicolas Gonzalez Granda Sunay Joshi Nagaprasad Rudrapatna 
*          Anne Somalwar 
*
* Description: This implements a FoodDensity object.
* Only the grid values are stored as an array.
* The speed and cost functions are stored only as function pointers, which are treated as inputs.
*
* (See also CFoodDensity.hpp)
*
* ==============================================================================
*/
/* ------ Libraries ---------------------------------------------------------*/
#include <string>
#include "boost/multi_array.hpp"

/* ------ Project-specific header files -------------------------------------*/
#include "CTerrain.hpp"
#include "WritetoFile.hpp"

/* ------ Namespaces --------------------------------------------------------*/
using namespace memory;
using namespace std;

/*==============================================================================
  Constructors
================================================================================*/

/* Function-based constructor */
CTerrain::CTerrain(function<double(double, double, double, double)> aSpottingRate,
                   function<double(double, double, double, double)> aKillRate,
                   function<double(double, double, double, double)> aGiveUpRate,
                   function<double(double, double)> aHomeBase,
                   function<double(double, double)> aObstacle,
                   const int aN, const double aPhysMin,
                   const double aPhysMax) {

  assert(aPhysMax >= aPhysMin);

  // Assign these member variables *before* calling x/yGridToPhysical below
  fMinX = aPhysMin;
  fMinY = aPhysMin;

  fMaxX = aPhysMax;
  fMaxY = aPhysMax;

  fH = (aPhysMax - aPhysMin)/(aN - 1);

  /* Initialize and fill in fFoodDensity, fHomeBase, fObstacle */
  fSpottingRate = make_shared<array2D_t<double>>(allocateArray2D<double>(aN,aN));
  fKillRate = make_shared<array2D_t<double>>(allocateArray2D<double>(aN,aN));
  fGiveUpRate = make_shared<array2D_t<double>>(allocateArray2D<double>(aN,aN));
  fHomeBase = make_shared<array2D_t<double>>(allocateArray2D<double>(aN,aN));
  fObstacle = make_shared<array2D_t<double>>(allocateArray2D<double>(aN,aN));
  fSleepingSite = make_shared<array2D_t<double>>(allocateArray2D<double>(aN,aN));

  for (int i = 0; i < aN; i++) {
    for (int j = 0; j < aN; j++) {
      double x = xGridToPhysical(i);
      double y = yGridToPhysical(j);
      (*fSpottingRate)[i][j] = aSpottingRate(x,y,0,0);
      (*fKillRate)[i][j] = aKillRate(x,y,0,0);
      (*fGiveUpRate)[i][j] = aGiveUpRate(x,y,0,0);
      (*fHomeBase)[i][j] = aHomeBase(x,y);
      (*fObstacle)[i][j] = aObstacle(x,y);
      (*fSleepingSite)[i][j] = 1;
    }
  }
}

/* Array-based constructor */
CTerrain::CTerrain(shared_ptr<array2D_t<double>> aSpottingRate,
                   shared_ptr<array2D_t<double>> aKillRate,
                   shared_ptr<array2D_t<double>> aGiveUpRate,
                   shared_ptr<array2D_t<double>> aHomeBase,
                   shared_ptr<array2D_t<double>> aObstacle,
                   const int aN, const double aPhysMin,
                   const double aPhysMax, 
                   shared_ptr<array2D_t<double>> aSleepingSite) {

  assert(aPhysMax >= aPhysMin);

  // Assign these member variables *before* calling x/yGridToPhysical below
  fMinX = aPhysMin;
  fMinY = aPhysMin;

  fMaxX = aPhysMax;
  fMaxY = aPhysMax;

  fH = (aPhysMax - aPhysMin)/(aN - 1);

  fHomeBase = aHomeBase;
  fObstacle = aObstacle;
  fSpottingRate = aSpottingRate;
  fKillRate = aKillRate;
  fGiveUpRate = aGiveUpRate;

  if (aSleepingSite) {
    fSleepingSite = aSleepingSite;
  } else {
    for (int i = 0; i < aN; i++) {
      for (int j = 0; j < aN; j++) {
        (*fSleepingSite)[i][j] = 1;
      }
    }
  }
}


void CTerrain::writeTerrainToFile(const string aFilename) const {
  const int nx = getGridSizeX();
  const int ny = getGridSizeY();

  shared_ptr<array2D_t<double>> spottingRate = make_shared<array2D_t<double>>(allocateArray2D<double>(nx,ny));
  shared_ptr<array2D_t<double>> killRate = make_shared<array2D_t<double>>(allocateArray2D<double>(nx,ny));
  shared_ptr<array2D_t<double>> giveUpRate = make_shared<array2D_t<double>>(allocateArray2D<double>(nx,ny));

  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      (*spottingRate)[i][j] = getSpottingRate(i,j);
      (*killRate)[i][j] = getKillRate(i,j);
      (*giveUpRate)[i][j] = getGiveUpRate(i,j);
    }
  }


  io::writeToFile2D<double>(aFilename+"_PredatorDensity", *spottingRate);
  io::writeToFile2D<double>(aFilename+"_KillRate", *killRate);
  io::writeToFile2D<double>(aFilename+"_GiveUpRate", *giveUpRate);

  io::writeToFile2D<double>(aFilename+"_HomeBase", *fHomeBase);
  io::writeToFile2D<double>(aFilename+"_Obstacle", *fObstacle);

  io::writeToFile2D<double>(aFilename+"_SleepingSite", *fSleepingSite);
}

/*==============================================================================
  Returns spotting rate at the physical location (aX,aY,aE,aT)
      using bilinear interpolation
==============================================================================*/

double CTerrain::getSpottingRatePhysical(const double aX, const double aY, const double aE, const double aT) const {

  /* Grid sizes  */
  const int nx = getGridSizeX();
  const int ny = getGridSizeY();

  /* Find grid coordinates of "lower-left corner" */
  /* Set k, l to 0, since functions do not depend on energy and time */
  const int i = floor((aX - fMinX) / fH);
  const int j = floor((aY - fMinY) / fH);
  const int k = 0;
  const int l = 0;

  /* Find interpolation coefficients */
  const double r_x = (aX - i*fH) / fH;
  const double r_y = (aY - j*fH) / fH;

  /* Find gridpoints for neighbors (used to handle boundary issues) */
  int i1 = i + 1;
  int j1 = j + 1;
  if (i == nx-1) { /* Check if point is on right side of grid */
    i1 = i;
  }
  if (j == ny-1) { /* Check if point is on top side of grid */
    j1 = j;
  }

  /* Interpolate along y-axis */
  const double mus1 = (1-r_y) * getSpottingRate(i,j,k,l) + r_y * getSpottingRate(i,j1,k,l);
  const double mus2 = (1-r_y) * getSpottingRate(i1,j,k,l) + r_y * getSpottingRate(i1,j1,k,l);

  /* Interpolate along x-axis */
  const double mus3 = (1-r_x) * mus1 + r_x * mus2;

  // Determine if we are inside the homebase
  int notHomeBase = 1;
  if (isHomeBasePhysical(aX, aY)) {
  // if (false) { // for Zeno
    notHomeBase = 0;
  }

  return mus3 * notHomeBase;
}

/*==============================================================================
  Returns kill rate at the physical location (aX,aY,aE,aT)
      using bilinear interpolation
==============================================================================*/

double CTerrain::getKillRatePhysical(const double aX, const double aY, const double aE, const double aT) const {

  /* Grid sizes  */
  const int nx = getGridSizeX();
  const int ny = getGridSizeY();

  /* Find grid coordinates of "lower-left corner" */
  /* Set k, l to 0, since functions do not depend on energy and time */
  const int i = floor((aX - fMinX) / fH);
  const int j = floor((aY - fMinY) / fH);
  const int k = 0;
  const int l = 0;

  /* Find interpolation coefficients */
  const double r_x = (aX - i*fH) / fH;
  const double r_y = (aY - j*fH) / fH;

  /* Find gridpoints for neighbors (used to handle boundary issues) */
  int i1 = i + 1;
  int j1 = j + 1;
  if (i == nx-1) { /* Check if point is on right side of grid */
    i1 = i;
  }
  if (j == ny-1) { /* Check if point is on top side of grid */
    j1 = j;
  }

  // Bilinear interpolation

  /* Interpolate along y-axis */
  const double muk1 = (1-r_y) * getKillRate(i,j,k,l) + r_y * getKillRate(i,j1,k,l);
  const double muk2 = (1-r_y) * getKillRate(i1,j,k,l) + r_y * getKillRate(i1,j1,k,l);

  /* Interpolate along x-axis */
  const double muk3 = (1-r_x) * muk1 + r_x * muk2;

  // Determine if we are inside the homebase
  int notHomeBase = 1;
  if (isHomeBasePhysical(aX, aY)) {
  // if (false) { // for Zeno
    notHomeBase = 0;
  }

  return muk3 * notHomeBase;
}

/*==============================================================================
  Returns giveup rate at the physical location (aX,aY,aE,aT)
      using bilinear interpolation
==============================================================================*/

double CTerrain::getGiveUpRatePhysical(const double aX, const double aY, const double aE, const double aT) const {

  /* Grid sizes  */
  const int nx = getGridSizeX();
  const int ny = getGridSizeY();

  /* Find grid coordinates of "lower-left corner" */
  /* Set k, l to 0, since functions do not depend on energy and time */
  const int i = floor((aX - fMinX) / fH);
  const int j = floor((aY - fMinY) / fH);
  const int k = 0;
  const int l = 0;

  /* Find interpolation coefficients */
  const double r_x = (aX - i*fH) / fH;
  const double r_y = (aY - j*fH) / fH;

  /* Find gridpoints for neighbors (used to handle boundary issues) */
  int i1 = i + 1;
  int j1 = j + 1;
  if (i == nx-1) { /* Check if point is on right side of grid */
    i1 = i;
  }
  if (j == ny-1) { /* Check if point is on top side of grid */
    j1 = j;
  }

  // Bilinear interpolation

  /* Interpolate along y-axis */
  const double mug1 = (1-r_y) * getGiveUpRate(i,j,k,l) + r_y * getGiveUpRate(i,j1,k,l);
  const double mug2 = (1-r_y) * getGiveUpRate(i1,j,k,l) + r_y * getGiveUpRate(i1,j1,k,l);

  /* Interpolate along x-axis */
  const double mug3 = (1-r_x) * mug1 + r_x * mug2;

  // Determine scaling depending on if we are inside homebase
  int scalingFactor = 1;
  if (isHomeBasePhysical(aX, aY)) {
  // if (false) { // for Zeno
    scalingFactor = 3;
  }

  // Multiply by scaling factor
  return mug3 * scalingFactor;
}

/*==============================================================================
  Returns a conservative estimate of whether the point (aX, aY) is inside of an 
      obstacle by checking all adjacent grid points
==============================================================================*/

bool CTerrain::isObstaclePhysicalConservative(const double aX, const double aY) const {
  /* Grid sizes  */
  const int nx = getGridSizeX();
  const int ny = getGridSizeY();

  /* Find grid coordinates of "lower-left corner" */
  const int i = floor((aX - fMinX) / fH);
  const int j = floor((aY - fMinY) / fH);

  /* Find gridpoints for neighbors (used to handle boundary issues) */
  int i1 = i + 1;
  int j1 = j + 1;
  if (i == nx-1) { /* Check if point is on right side of grid */
    i1 = i;
  }
  if (j == ny-1) { /* Check if point is on top side of grid */
    j1 = j;
  }

  bool value = false;
  if (((*fObstacle)[i][j] == 0) || ((*fObstacle)[i1][j] == 0) 
      || ((*fObstacle)[i][j1] == 0) || ((*fObstacle)[i1][j1] == 0)) {
    value = true;
  }
  return value;
}

bool CTerrain::isSleepingSitePhysicalConservative(const double aX, const double aY) const {
  /* Grid sizes  */
  const int nx = getGridSizeX();
  const int ny = getGridSizeY();

  /* Find grid coordinates of "lower-left corner" */
  const int i = floor((aX - fMinX) / fH);
  const int j = floor((aY - fMinY) / fH);

  /* Find gridpoints for neighbors (used to handle boundary issues) */
  int i1 = i + 1;
  int j1 = j + 1;
  if (i == nx-1) { /* Check if point is on right side of grid */
    i1 = i;
  }
  if (j == ny-1) { /* Check if point is on top side of grid */
    j1 = j;
  }

  bool value = false;
  if (((*fSleepingSite)[i][j] == 1) && ((*fSleepingSite)[i1][j] == 1) 
      && ((*fSleepingSite)[i][j1] == 1) && ((*fSleepingSite)[i1][j1] == 1)) {
    value = true;
  }
  return value;
}

