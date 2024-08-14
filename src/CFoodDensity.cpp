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
* File: CFoodDensity.cpp
*
* Authors: Marissa Gee Nicolas Gonzalez Granda Sunay Joshi Nagaprasad Rudrapatna 
*          Anne Somalwar 
*
* Description: This class stores and provides access to environment properties
* that depend only on an animal's position.
*
* (See also CTerrain.hpp)
*
* ==============================================================================
*/
/* ------ Libraries ---------------------------------------------------------*/
#include <string>
#include "boost/multi_array.hpp"

/* ------ Project-specific header files -------------------------------------*/
#include "CFoodDensity.hpp"
#include "WritetoFile.hpp"

/* ------ Namespaces --------------------------------------------------------*/
using namespace memory;
using namespace std;

/*==============================================================================
  Constructors
================================================================================*/

/* Function-based constructor */
CFoodDensity::CFoodDensity(function<double(double, double)> aFoodDensity,
                           function<double(double, double)> aFoodKernel,
                           function<double(double)> aUtilityFunction,
                           shared_ptr<CTerrain> aTerrain,
                           const double aRadius, const double aHarvestRate,
                           const bool aFoodDepletion,const int aN, 
                           const double aPhysMin, const double aPhysMax) {

  assert(aPhysMax >= aPhysMin);
  assert(aRadius > 0.0);

  // Assign these member variables *before* calling x/yGridToPhysical below
  fMinX = aPhysMin;
  fMinY = aPhysMin;
  fMaxX = aPhysMax;
  fMaxY = aPhysMax;

  fH = (aPhysMax - aPhysMin)/(aN-1);

  fFoodKernel = aFoodKernel;
  fUtilityFunction = aUtilityFunction;
  fRadius = aRadius;
  fHarvestRate = aHarvestRate;
  fTerrain = aTerrain;
  fFoodDepletion = aFoodDepletion;

  /* Initialize and fill in fFoodDensity, fObstacle */
  fFoodDensity = make_shared<array2D_t<double>>(allocateArray2D<double>(aN,aN));
  fInitialFoodDensity = make_shared<array2D_t<double>>(allocateArray2D<double>(aN,aN));
  fFoodAccumulated = make_shared<array2D_t<double>>(allocateArray2D<double>(aN,aN));
  for (int i = 0; i < aN; i++) {
    for (int j = 0; j < aN; j++) {
      double x = xGridToPhysical(i);
      double y = yGridToPhysical(j);
      double homebase = fTerrain->getHomeBase(i,j);
      (*fInitialFoodDensity)[i][j] = aFoodDensity(x,y)*homebase;
      (*fFoodDensity)[i][j] = aFoodDensity(x,y)*homebase;
    }
  }

  // Only calculate food accumulated after we have set the food density everywhere
  for (int i = 0; i < aN; i++) {
    for (int j = 0; j < aN; j++) {
      const double x = xGridToPhysical(i);
      const double y = yGridToPhysical(j);
      const double food_accumulated = calcFoodAccumulated(x,y);
      (*fFoodAccumulated)[i][j] = food_accumulated;
    }
  }
}

/* Function-based constructor with food accumulated argument */
CFoodDensity::CFoodDensity(function<double(double, double)> aFoodDensity,
                           function<double(double, double)> aFoodKernel,
                           function<double(double, double)> aFoodAccumulated,
                           function<double(double)> aUtilityFunction,
                           shared_ptr<CTerrain> aTerrain,
                           const double aRadius, const double aHarvestRate,
                           const bool aFoodDepletion, const int aN, 
                           const double aPhysMin, const double aPhysMax) {

  assert(aPhysMax >= aPhysMin);
  assert(aRadius > 0.0);

  // Assign these member variables *before* calling x/yGridToPhysical below
  fMinX = aPhysMin;
  fMinY = aPhysMin;
  fMaxX = aPhysMax;
  fMaxY = aPhysMax;

  fH = (aPhysMax - aPhysMin)/(aN-1);

  fFoodKernel = aFoodKernel;
  fUtilityFunction = aUtilityFunction;
  fRadius = aRadius;
  fHarvestRate = aHarvestRate;
  fTerrain = aTerrain;

  /* Initialize and fill in fFoodDensity */
  fFoodDensity = make_shared<array2D_t<double>>(allocateArray2D<double>(aN,aN));
  fInitialFoodDensity = make_shared<array2D_t<double>>(allocateArray2D<double>(aN,aN));
  fFoodAccumulated = make_shared<array2D_t<double>>(allocateArray2D<double>(aN,aN));
  for (int i = 0; i < aN; i++) {
    for (int j = 0; j < aN; j++) {
      double x = xGridToPhysical(i);
      double y = yGridToPhysical(j);
      double homebase = fTerrain->getHomeBase(i,j);
      (*fInitialFoodDensity)[i][j] = aFoodDensity(x,y)*homebase;
      (*fFoodDensity)[i][j] = aFoodDensity(x,y)*homebase;
    }
  }

  // Only calculate food accumulated after we have set the food density everywhere
  for (int i = 0; i < aN; i++) {
    for (int j = 0; j < aN; j++) {
      const double x = xGridToPhysical(i);
      const double y = yGridToPhysical(j);
      const double food_accumulated = aFoodAccumulated(x,y);
      (*fFoodAccumulated)[i][j] = food_accumulated;
    }
  }
}

/* Array-based constructor */
CFoodDensity::CFoodDensity(shared_ptr<array2D_t<double>> aFoodDensity,
                           function<double(double, double)> aFoodKernel,
                           function<double(double)> aUtilityFunction,
                           shared_ptr<CTerrain> aTerrain,
                           const double aRadius, const double aHarvestRate, 
                           const bool aFoodDepletion, const int aN, 
                           const double aPhysMin, const double aPhysMax) {

  assert(aPhysMax >= aPhysMin);
  assert(aRadius > 0.0);
  fMinX = aPhysMin;
  fMinY = aPhysMin;
  fMaxX = aPhysMax;
  fMaxY = aPhysMax;

  fH = (aPhysMax - aPhysMin)/(aN-1);

  /* assign fFoodDensity, fTerrain, fFoodKernel, fRadius */
  fTerrain = aTerrain;
  fFoodKernel = aFoodKernel;
  fUtilityFunction = aUtilityFunction;
  fRadius = aRadius;
  fHarvestRate = aHarvestRate;

  fFoodDensity = make_shared<array2D_t<double>>(allocateArray2D<double>(aN,aN));
  fFoodAccumulated = make_shared<array2D_t<double>>(allocateArray2D<double>(aN,aN));
  fInitialFoodDensity = aFoodDensity;
  for (int i = 0; i < aN; i++) {
    for (int j = 0; j < aN; j++) {
      const double homebase = fTerrain->getHomeBase(i,j);
      (*fFoodDensity)[i][j] = (*aFoodDensity)[i][j]*homebase;
    }
  }

  for (int i = 0; i < aN; i++) {
    for (int j = 0; j < aN; j++) {
      const double x = xGridToPhysical(i);
      const double y = yGridToPhysical(j);
      const double food_accumulated = calcFoodAccumulated(x,y);
      (*fFoodAccumulated)[i][j] = food_accumulated;
    }
  }
}


/*==============================================================================
  Modifying food density
==============================================================================*/
void CFoodDensity::setFoodDensity(shared_ptr<array2D_t<double>> aFoodDensity) {
  const int nx = getGridSizeX();
  const int ny = getGridSizeY();
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      setFoodDensityPoint(i, j, (*aFoodDensity)[i][j]);
    }
  }
}

void CFoodDensity::resetFoodDensity() {
  setFoodDensity(fInitialFoodDensity);
}

/*==============================================================================
  Write grid to file
==============================================================================*/
void CFoodDensity::writeFoodGridToFile(const string aFilename) const {

  /* Obtain grid sizes */
  int nx = getGridSizeX();
  int ny = getGridSizeY();

  /* Initialize and fill foodDensityArray */
  array2D_t<double> foodDensityArray = allocateArray2D<double>(nx, ny);
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      foodDensityArray[i][j] = getFoodDensity(i,j);
    }
  }

  io::writeToFile2D<double>(aFilename + "_FoodDensity", foodDensityArray);
  io::writeToFile2D<double>(aFilename + "_FoodAccumulated", *fFoodAccumulated);
}

void CFoodDensity::writeFoodParamsToFile(const string aFilename) const {
  // writing harvest rate to file
  ofstream outHarvestRate("output/" + aFilename + "_HarvestRate", ios::binary);
  outHarvestRate.write((char*) &fHarvestRate, sizeof(double));

  // writing radius to file
  ofstream outRadius("output/" + aFilename + "_Radius", ios::binary);
  outRadius.write((char*) &fRadius, sizeof(double));
}

/*==============================================================================
  Food accumulation and depletion
==============================================================================*/

/* Update entire food accumulated grid*/
void CFoodDensity::updateFoodAccumulated() {
  const int nx = getGridSizeX();
  const int ny = getGridSizeY();

  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      const double x = xGridToPhysical(i);
      const double y = yGridToPhysical(j);
      (*fFoodAccumulated)[i][j] = calcFoodAccumulated(x,y);
    }
  }
}

/* Food accumulation when animal is at (aX,aY)*/
double CFoodDensity::calcFoodAccumulated(const double aX, const double aY) {

  double food = 0.0;
  double normalization = 0.0;

  /* Obtain grid sizes */
  const int nx = getGridSizeX();
  const int ny = getGridSizeY();

  /* Obtain grid spacing */
  const double h = getH();

  /* Current gridpoint */
  int curr_i = xPhysicalToGrid(aX);
  int curr_j = yPhysicalToGrid(aY);

  /* If fRadius < h, then the animal accumulates all the food at the current point */
  if (fRadius < h) {
    return getFoodDensityPhysical(aX, aY);
  }

  /* Sidelength of square of interest */
  int num = ceil(fRadius / h);

  /* Boundaries of square
   * Need 0 <= i <= nx-2, 0 <= j <= ny-2 */
  int left_idx = max(curr_i - num - 1, 0);
  int right_idx = min(curr_i + num + 1, nx-1);
  int bottom_idx = max(curr_j - num - 1, 0);
  int top_idx = min(curr_j + num + 1, ny-1);

  /* Loop over square of sidelength num centered on (curr_i, curr_j) */
  int counter = 0;
  for (int i = left_idx; i <= right_idx; ++i) {
    for (int j = bottom_idx; j <= top_idx; ++j) {
      const double p_x = xGridToPhysical(i);
      const double p_y = yGridToPhysical(j);
      const double dist = sqrt(pow((aX-p_x), 2) + pow((aY-p_y), 2));

      // Only consider point in square if the distance from the current point is at most fRadius
      if (dist <= fRadius + EPSILON) {
        ++counter;
        double psi_point;
        if (getFoodDepletion()) {
          psi_point = getFoodDensity(i,j);
        } else {
          psi_point = getInitialFoodDensity(i,j);
        }
        food += getFoodKernel(dist) * getHarvestRate() * psi_point * h * h;
        if (!fTerrain->isObstacle(i,j)){
          normalization += getFoodKernel(dist) * h * h;
        }
      } 
    }
  }

  if (normalization == 0) {
    if (food == 0) {
      normalization = 1;
    } else {
      cout << "Error calculating food accumulated." << endl;
      assert(false);
    }
  }
  return food/normalization;
}

double CFoodDensity::getFoodAccumulatedPhysical(const double aX, const double aY) const {
  /* Grid sizes  */
  const int nx = getGridSizeX();
  const int ny = getGridSizeY();

  /* Find grid coordinates of "lower-left corner" */
  const int i = floor((aX - fMinX) / fH);
  const int j = floor((aY - fMinY) / fH);

  /* Find interpolation coefficients */
  double r_x = (aX - i*fH) / fH;
  if ((r_x < 0) || (r_x > 1)) {
    if (abs(r_x) < EPSILON_FLOAT) {
      r_x = 0;
    } else if (abs(1-r_x) < EPSILON_FLOAT) {
      r_x = 1;
    } else {
      cout << "rx is " << r_x << endl;
      assert(false);
    }
  }
  assert((r_x >= 0) && (r_x <= 1));
  double r_y = (aY - j*fH) / fH;
  if ((r_y < 0) || (r_y > 1)) {
    if (abs(r_y) < EPSILON_FLOAT) {
      r_y = 0;
    } else if (abs(1-r_y) < EPSILON_FLOAT) {
      r_y = 1;
    } else {
      cout << "ry is " <<  r_y << endl;
      assert(false);
    }
  }
  assert((r_y >= 0) && (r_y <= 1));

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
  const double F1 = (1-r_y) * getFoodAccumulated(i,j) 
                    + r_y * getFoodAccumulated(i,j1);
  const double F2 = (1-r_y) * getFoodAccumulated(i1,j) 
                    + r_y * getFoodAccumulated(i1,j1);

  /* Interpolate along x-axis */
  return (1-r_x) * F1 + r_x * F2;
}

/* Food depletion after eating at (aX,aY) for aTau units of time
 * Returns totalFood, i.e. units of *energy* (not energy/time!) */
double CFoodDensity::depleteFood(const double aX, const double aY, 
                                 const double aTau, const double aFoodFactor) {
  return depleteFood(aX, aY, aTau, aFoodFactor, fFoodDensity);
}

double CFoodDensity::depleteFood(const double aX, const double aY, 
                                 const double aTau, const double aFoodFactor,
                                 shared_ptr<array2D_t<double>> aFoodDensityArray) {

  /* Should be equivalent to calcFoodAccumulated if aTau = dt / TimeFactor */
  double totalFood = 0.0; // total food depleted (units = energy)

  /* Obtain grid sizes */
  const int nx = getGridSizeX();
  const int ny = getGridSizeY();

  /* Obtain grid spacing */
  const double h = getH();

  /* Current gridpoint */
  int curr_i = xPhysicalToGrid(aX);
  int curr_j = yPhysicalToGrid(aY);

  /* If fRadius < h, then the animal depletes all the food at the current point */
  if (fRadius < h) {
    totalFood = getFoodDensity(curr_i, curr_j) * aTau;
    // cout << "I just ate this much food! : " << totalFood << endl;

    return totalFood;
  }

  /* Sidelength of square of interest */
  int num = ceil(fRadius / h);

  /* Boundaries of square
   * Need 0 <= i <= nx-2, 0 <= j <= ny-2 */
  int left_idx = max(curr_i - num - 1, 0);
  int right_idx = min(curr_i + num + 1, nx-1);
  int bottom_idx = max(curr_j - num - 1, 0);
  int top_idx = min(curr_j + num + 1, ny-1);

  double normalization = 0;
  for (int i = left_idx; i <= right_idx; ++i) {
    for (int j = bottom_idx; j <= top_idx; ++j) {
      const double x = xGridToPhysical(i);
      const double y = yGridToPhysical(j);
      const double dist = sqrt(pow(x - aX, 2) + pow(y - aY, 2));

      // if ((dist <= fRadius + EPSILON) && (fTerrain->isRegularTerrain(i,j))) {
      if ((dist <= fRadius + EPSILON) && (!fTerrain->isObstacle(i,j))) {
        normalization += getFoodKernel(dist) * h * h;
      }
    }
  }


  if (normalization == 0) {
    // No forageable terrain
    totalFood = 0;
  } else {
    // Deplete food
    for (int i = left_idx; i <= right_idx; ++i) {
      for (int j = bottom_idx; j <= top_idx; ++j) {
        const double x = xGridToPhysical(i);
        const double y = yGridToPhysical(j);
        const double dist = sqrt(pow(x - aX, 2) + pow(y - aY, 2));

        // Only consider point in square if the distance from the current point is at most fRadius
        if (dist <= fRadius + EPSILON) {
          double psi_old = 0;
          if (!fTerrain->isObstacle(i,j)) {
            psi_old = (*aFoodDensityArray)[i][j];
          }
          // Deplete food from the entire pack of animals
          const double psi_change = psi_old * (getFoodKernel(dist)/normalization) 
                                    * getHarvestRate() * aTau * aFoodFactor;
          const double psi_new = psi_old - psi_change;
          (*aFoodDensityArray)[i][j] = psi_new;

          // Food depleted / area = sum of all changes in food density
          totalFood += psi_change;
        }
      }
    }
  }

  // Need to account for area (h^2) and divide by number of animals
  totalFood = totalFood * h * h;
  return totalFood;
}
