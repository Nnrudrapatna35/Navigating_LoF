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
* File: CTimeGrid.cpp
*
* Authors: REU 2022 Landscape of Fear Group
*   (based on code by Elliot Cartee, Qianli Song, Lexiao Lai)
*     (based on code by Marc Aurèle Gilles and Zachary Clawson)
*
* Description: This class is used to handle both the logical and physical
* representation of a 3D regular grid with fixed spacing in all 3 dimensions.
* It is used by CFMM, CTimeDependentHjbSolver, CTimeDependentTracer, and
* CMovingObserver. This class is an interface for between the underlying
* data structures and the rest of the code.
* Only the grid values are stored as an array. The speed and cost
* functions are stored only as function pointers, which are treated as inputs.
*
* (See also CTimeGrid.hpp)
*
* ==============================================================================
*/
/* ------ Libraries ---------------------------------------------------------*/
#include <string>
#include <iomanip>
#include "boost/multi_array.hpp"

/* ------ Project-specific header files -------------------------------------*/
#include "CTimeGrid.hpp"
#include "WritetoFile.hpp"

/* ------ Namespaces --------------------------------------------------------*/
using namespace memory;
using namespace std;

/*==============================================================================
  Constructor
================================================================================*/

// Function-based constructor
CTimeGrid::CTimeGrid(function<double(double, double, double, double)> aCost,
                     function<double(double, double, double, double)> aSpeed,
                     const int aN, const int aNe, const int aNt,
                     const double aPhysMin, const double aPhysMax,
                     const double aEnergyMin, const double aEnergyMax,
                     const double aTime,const bool aFullGrid, const int aCurrentL) {

  assert(aPhysMax >= aPhysMin);
  assert(aEnergyMax >= aEnergyMin);
  assert(aTime > 0.0);

  fMinX = aPhysMin;
  fMinY = aPhysMin;
  fMinE = aEnergyMin;
  fMinT = 0.0;
  fMaxX = aPhysMax;
  fMaxY = aPhysMax;
  fMaxE = aEnergyMax;
  fMaxT = aTime;

  fFullGrid = aFullGrid;
  fCurrentL = aCurrentL;

  fH = (aPhysMax - aPhysMin)/(aN-1);
  fE = (aEnergyMax - aEnergyMin)/(aNe-1);
  fDt = aTime / (aNt-1);

  /* Initialize and fill all arrays (fValues, fCost, fPredatorDensity, fSpeed, fMuG, fMuK, fHomeBase) */
  if (fFullGrid) {
    fValues = make_shared<array4D_t<float>>(allocateArray4D<float>(aNt,aN,aN,aNe));
  } else {
    fValuesCurrent = make_shared<array3D_t<float>>(allocateArray3D<float>(aN,aN,aNe));
    fValuesFuture = make_shared<array3D_t<float>>(allocateArray3D<float>(aN,aN,aNe));
  }
  fCost = make_shared<array2D_t<double>>(allocateArray2D<double>(aN,aN));
  fSpeed = make_shared<array2D_t<double>>(allocateArray2D<double>(aN,aN));

  // conversions to physical
  for(int i = 0; i < aN; i++) {
    double x = xGridToPhysical(i);
    for(int j = 0; j < aN; j++) {
      double y = yGridToPhysical(j);

      // While we have generally allowed these to depend on energy and time,
      // none of our examples currently make use of that dependence, and having
      // an unnecessary 4D for loop here slows things down significantly, so for
      // now we will assume no energy or time dependence
      (*fCost)[i][j] = aCost(x,y,0,0);
      (*fSpeed)[i][j] = aSpeed(x,y,0,0);
    }
  }
}

// Array-based constructor
CTimeGrid::CTimeGrid(shared_ptr<array2D_t<double>> aCost,
                     shared_ptr<array2D_t<double>> aSpeed,
                     const int aN, const int aNe, const int aNt,
                     const double aPhysMin, const double aPhysMax,
                     const double aEnergyMin, const double aEnergyMax,
                     const double aTime, const bool aFullGrid, const int aCurrentL) {

  assert(aPhysMax >= aPhysMin);
  assert(aEnergyMax >= aEnergyMin);
  assert(aTime > 0.0);
  fMinX = aPhysMin;
  fMinY = aPhysMin;
  fMinE = aEnergyMin;
  fMinT = 0.0;
  fMaxX = aPhysMax;
  fMaxY = aPhysMax;
  fMaxE = aEnergyMax;
  fMaxT = aTime;

  fFullGrid = aFullGrid;
  fCurrentL = aCurrentL;

  fH = (aPhysMax - aPhysMin)/(aN-1);
  fE = (aEnergyMax - aEnergyMin)/(aNe-1);
  fDt = aTime / (aNt-1);

  fCost = aCost;
  fSpeed = aSpeed;

  /* Initialize and fill fValues array */
  if (fFullGrid) {
    fValues = make_shared<array4D_t<float>>(allocateArray4D<float>(aNt,aN,aN,aNe));
  } else {
    fValuesCurrent = make_shared<array3D_t<float>>(allocateArray3D<float>(aN,aN,aNe));
    fValuesFuture = make_shared<array3D_t<float>>(allocateArray3D<float>(aN,aN,aNe));
  }
}

/*==============================================================================
  Write grid to file
==============================================================================*/
void CTimeGrid::writeGridToFile(const string aFilename) const {
  const int nx = getGridSizeX();
  const int ny = getGridSizeY();
  const int ne = getGridSizeE();
  const int nt = getGridSizeT();

  vector<double> params;
  params.push_back(nx);
  params.push_back(ny);
  params.push_back(ne);
  params.push_back(nt);
  params.push_back(fMinX);
  params.push_back(fMaxX);
  params.push_back(fDt);
  params.push_back(fH);
  params.push_back(fMinE);
  params.push_back(fMaxE);

  // Writing nx,ny,ne,nt to file, so they can be read in in Matlab
  io::writeVectorToFile<double>(aFilename + "Parameters", params);
  io::writeToFile2D<double>(aFilename + "Cost",  *fCost);
  io::writeToFile2D<double>(aFilename + "Speed",  *fSpeed);

  if (fFullGrid) {
    io::writeToFile4D<float>(aFilename + "Value", *fValues);
  }
}

void CTimeGrid::writeSliceToFile(const string aFilename, const int aSlice) const {
  assert(!fFullGrid);
  if (aSlice == fCurrentL) {
    io::writeToFile3D<float>(aFilename + "Value_Slice" + to_string(fCurrentL), 
                             *fValuesCurrent);
  } else if (aSlice == fCurrentL + 1) {
    io::writeToFile3D<float>(aFilename + "Value_Slice" + to_string(fCurrentL+1), 
                             *fValuesFuture);
  } else {
    cout << "Invalid slice requested." << endl;
    assert(false);
  }
}

/*==============================================================================
  Returns speed at the physical location (aX,aY,aE,aT)
      using bilinear interpolation
==============================================================================*/
double CTimeGrid::getSpeedPhysical(const double aX, const double aY, const double aE,
                        const double aT) const {

  /* Grid sizes  */
  const int nx = getGridSizeX();
  const int ny = getGridSizeY();
  const int ne = getGridSizeE();
  const int nt = getGridSizeT();

  /* Find grid coordinates of "lower-left corner" */
  const int i = floor((aX - fMinX) / fH);
  const int j = floor((aY - fMinY) / fH);
  const int k = floor((aE - fMinE) / fE);
  const int l = floor((aT - fMinT) / fDt);

  /* Find interpolation coefficients */
  double r_x = (aX - i*fH) / fH;
  if ((r_x < 0) || (r_x > 1)) {
    if (abs(r_x) < EPSILON_FLOAT) {
      r_x = 0;
    } else if (abs(1-r_x) < EPSILON_FLOAT) {
      r_x = 1;
    } else {
      cout << "r_x is " << r_x << endl;
      assert(false);
    }
  }
  double r_y = (aY - j*fH) / fH;
  if ((r_y < 0) || (r_y > 1)) {
    if (abs(r_y) < EPSILON_FLOAT) {
      r_y = 0;
    } else if (abs(1-r_y) < EPSILON_FLOAT) {
      r_y = 1;
    } else {
      cout << "r_y is " << r_y << endl;
      assert(false);
    }
  }

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
  const double f1 = (1-r_y) * getSpeed(i,j,k,l) + r_y * getSpeed(i,j1,k,l);
  const double f2 = (1-r_y) * getSpeed(i1,j,k,l) + r_y * getSpeed(i1,j1,k,l);

  /* Interpolate along x-axis */
  return (1-r_x) * f1 + r_x * f2;
}

/*==============================================================================
  Returns cost at the physical location (aX,aY,aE,aT)
      using bilinear interpolation
==============================================================================*/
double CTimeGrid::getCostPhysical(const double aX, const double aY, const double aE,
                        const double aT) const {

  /* Grid sizes  */
  const int nx = getGridSizeX();
  const int ny = getGridSizeY();
  const int ne = getGridSizeE();
  const int nt = getGridSizeT();

  /* Find grid coordinates of "lower-left corner" */
  const int i = floor((aX - fMinX) / fH);
  const int j = floor((aY - fMinY) / fH);
  const int k = floor((aE - fMinE) / fE);
  const int l = floor((aT - fMinT) / fDt);

  /* Find interpolation coefficients */
  double r_x = (aX - i*fH) / fH;
  if ((r_x < 0) || (r_x > 1)) {
    if (abs(r_x) < EPSILON_FLOAT) {
      r_x = 0;
    } else if (abs(1-r_x) < EPSILON_FLOAT) {
      r_x = 1;
    } else {
      cout << "r_x is " << r_x << endl;
      assert(false);
    }
  }
  double r_y = (aY - j*fH) / fH;
  if ((r_y < 0) || (r_y > 1)) {
    if (abs(r_y) < EPSILON_FLOAT) {
      r_y = 0;
    } else if (abs(1-r_y) < EPSILON_FLOAT) {
      r_y = 1;
    } else {
      cout << "r_y is " << r_y << endl;
      assert(false);
    }
  }

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
  const double K1 = (1-r_y) * getCost(i,j,k,l) + r_y * getCost(i,j1,k,l);
  const double K2 = (1-r_y) * getCost(i1,j,k,l) + r_y * getCost(i1,j1,k,l);

  /* Interpolate along x-axis */
  return (1-r_x) * K1 + r_x * K2;
}


/*==============================================================================
  Returns value function at the physical location (aX,aY,aE,aT)
      using quadrilinear interpolation
==============================================================================*/
double CTimeGrid::getValuePhysical(const double aX, const double aY, const double aE,
                                   const double aT) const {

  // if (fFullGrid == false) {
  //   assert(floor((aT - fMinT) / fDt) >= fCurrentL);
  //   assert(ceil((aT - fMinT) / fDt) <= fCurrentL + 1);
  // }
  /* Grid sizes  */
  const int nx = getGridSizeX();
  const int ny = getGridSizeY();
  const int ne = getGridSizeE();
  const int nt = getGridSizeT();

  /* Find grid coordinates of "lower-left corner" */
  const int i = floor((aX - fMinX) / fH);
  const int j = floor((aY - fMinY) / fH);
  const int k = floor((aE - fMinE) / fE);
  int l = floor((aT - fMinT) / fDt);
  if (!fFullGrid) {
    if (((aT - fMinT) >= (fCurrentL*fDt - EPSILON)) && 
        ((aT - fMinT) <= ((fCurrentL  + 1)*fDt + EPSILON))) {
      l = fCurrentL;
    } else if (aT >= fMaxT && aT <= fMaxT + EPSILON_FLOAT) {
      assert(fCurrentL == nt - 1 || fCurrentL + 1 == nt - 1); 
      l = nt-1;
    } else {
      cout << setprecision(24) << aT << endl;
      cout << fCurrentL*fDt << endl;
      cout << (fCurrentL+1)*fDt << endl;
      cout << "Invalid time requested." << endl;
      assert(false);
    }
  }

  /* Find interpolation coefficients */
  double r_x = (aX - xGridToPhysical(i)) / fH;
  if ((r_x < 0) || (r_x > 1)) {
    if (abs(r_x) < EPSILON_FLOAT) {
      r_x = 0;
    } else if (abs(1-r_x) < EPSILON_FLOAT) {
      r_x = 1;
    } else {
      cout << "r_x is " << r_x << endl;
      assert(false);
    }
  }
  assert((r_x >= 0) && (r_x <= 1));
  double r_y = (aY - yGridToPhysical(j)) / fH;
  if ((r_y < 0) || (r_y > 1)) {
    if (abs(r_y) < EPSILON_FLOAT) {
      r_y = 0;
    } else if (abs(1-r_y) < EPSILON_FLOAT) {
      r_y = 1;
    } else {
      cout << "r_y is " << r_y << endl;
      assert(false);
    }
  }
  assert((r_y >= 0) && (r_y <= 1));
  double r_e = (aE - eGridToPhysical(k)) / fE;
  if ((r_e < 0) || (r_e > 1)) {
    if (abs(r_e) < EPSILON_FLOAT) {
      r_e = 0;
    } else if (abs(1-r_e) < EPSILON_FLOAT) {
      r_e = 1;
    } else {
      cout << "r_e is " << r_e << endl;
      assert(false);
    }
  }
  assert((r_e >= 0) && (r_e <= 1));
  double r_t = (aT - tGridToPhysical(l)) / fDt;
  if ((r_t < 0) || (r_t > 1)) {
    if (abs(r_t) < EPSILON_FLOAT) {
      r_t = 0;
    } else if (abs(1-r_t) < EPSILON_FLOAT) {
      r_t = 1;
    } else {
      cout << "r_t is " << r_t << endl;
      assert(false);
    }
  }
  assert((r_t >= 0) && (r_t <= 1));

  /* Find gridpoints for neighbors (used to handle boundary issues) */
  int i1 = i + 1;
  int j1 = j + 1;
  int k1 = k + 1;
  int l1 = l + 1;
  if (i == nx-1) { /* Check if point is on right side of grid */
    i1 = i;
  }
  if (j == ny-1) { /* Check if point is on top side of grid */
    j1 = j;
  }
  if (k == ne-1) { /* Check if point is at final energy level */
    k1 = k;
  }
  if (l == nt-1) { /* Check if point is at final timestep */
    l1 = l;
  }

  // Quadrilinear interpolation

  /* Interpolate along t-axis */
  const double u1 = (1-r_t)*getValue(i,j,k,l)   + (r_t)*getValue(i,j,k,l1);
  const double u2 = (1-r_t)*getValue(i1,j,k,l)  + (r_t)*getValue(i1,j,k,l1);
  const double u3 = (1-r_t)*getValue(i,j1,k,l)  + (r_t)*getValue(i,j1,k,l1);
  const double u4 = (1-r_t)*getValue(i,j,k1,l)  + (r_t)*getValue(i,j,k1,l1);
  const double u5 = (1-r_t)*getValue(i1,j1,k,l) + (r_t)*getValue(i1,j1,k,l1);
  const double u6 = (1-r_t)*getValue(i,j1,k1,l) + (r_t)*getValue(i,j1,k1,l1);
  const double u7 = (1-r_t)*getValue(i1,j,k1,l) + (r_t)*getValue(i1,j,k1,l1);
  const double u8 = (1-r_t)*getValue(i1,j1,k1,l) + (r_t)*getValue(i1,j1,k1,l1);


  /* Interpolate along e-axis */
  const double u9 = (1-r_e)*u1 + (r_e)*u4;
  const double u10 = (1-r_e)*u2 + (r_e)*u7;
  const double u11 = (1-r_e)*u3 + (r_e)*u6;
  const double u12 = (1-r_e)*u5 + (r_e)*u8;

  /* Interpolate along y-axis */
  const double u13 = (1-r_y)*u9 + (r_y)*u11;
  const double u14 = (1-r_y)*u10 + (r_y)*u12;


  /* Interpolate along x-axis */
  return (1-r_x)*u13 + (r_x)*u14;
}

void CTimeGrid::advanceSlicesBackward() {
  const int nx = getGridSizeX();
  const int ny = getGridSizeY();
  const int ne = getGridSizeE();

  fValuesFuture = fValuesCurrent;
  fValuesCurrent = make_shared<array3D_t<float>>(allocateArray3D<float>(nx,ny,ne));
  fCurrentL = fCurrentL - 1;
}

void CTimeGrid::advanceSlicesForward() {
  const int nx = getGridSizeX();
  const int ny = getGridSizeY();
  const int ne = getGridSizeE();

  fValuesCurrent = fValuesFuture;
  fValuesFuture = make_shared<array3D_t<float>>(allocateArray3D<float>(nx,ny,ne));
  fCurrentL = fCurrentL + 1;
}

void CTimeGrid::readValuesFromFile(const string aFilename, const int aSlice) {
  if (fFullGrid) {
    io::readFromFile4D<float>(aFilename, *fValues);
  } else {
    if (aSlice == fCurrentL) {
      io::readFromFile3D<float>(aFilename + "_Slice" + to_string(aSlice), *fValuesCurrent);
    } else if (aSlice == fCurrentL + 1) {
      io::readFromFile3D<float>(aFilename + "_Slice" + to_string(aSlice), *fValuesFuture);
    } else {
      cout << "Invalid slice requested." << endl;
      assert(false);
    }
  }
}

void CTimeGrid::readValuesFromFile(const int aSlice) {
  assert(!fFullGrid);
  if (aSlice == fCurrentL) {
    io::readFromFile3D<float>(fFilename + "_Slice" + to_string(aSlice), *fValuesCurrent);
  } else if (aSlice == fCurrentL + 1) {
    io::readFromFile3D<float>(fFilename + "_Slice" + to_string(aSlice), *fValuesFuture);
  } else {
    cout << "Invalid slice requested." << endl;
    assert(false);
  }
}