/*
* ==============================================================================
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
* File: CTimeGrid.hpp
*
* Authors: Marissa Gee Nicolas Gonzalez Granda Sunay Joshi 
*          Nagaprasad Rudrapatna Anne Somalwar 
*          (based on code by Elliot Cartee, Qianli Song, Lexiao Lai)
*             (based on code by Marc Aurèle Gilles and Zachary Clawson)
*
* Description: This class is used to handle both the logical and physical
* representation of a 3D regular grid with fixed spacing in all 3 dimensions.
* It is used by CFMM, CTimeDependentHjbSolver, CTimeDependentTracer, and
* CMovingObserver. This class is an interface for between the underlying
* data structures and the rest of the code.
* Only the grid values are stored as an array. The speed and cost
* functions are stored only as function pointers, which are treated as inputs.
*
* ==============================================================================
*/

#ifndef CTIMEGRID_HPP
#define CTIMEGRID_HPP

/** ----- Libraries ----------------------------------------------------------*/
#include <cmath>
#include <string>
#include <iostream>
#include "boost/multi_array.hpp"

/** ------ Project-specific header files -------------------------------------*/
#include "MemoryAllocations.hpp"
#include "SimpleFunctions.hpp"

class CTimeGrid
{
  private:
    /** Properties of grid */
    /* Grid spacing: assuming to be the same in the two spatial variables */
    double fH;
    /* Energy spacing */
    double fE;
    /* Length of timestep */
    double fDt;

    std::string fFilename;

    /* Grid of value function values. Stored as a pointer to a 4D boost::multi_array */
    /* (i,j,k,l) correspond to (x, y, energy, time) */
    std::shared_ptr<memory::array4D_t<float>> fValues;

    std::shared_ptr<memory::array3D_t<float>> fValuesCurrent;
    std::shared_ptr<memory::array3D_t<float>> fValuesFuture;

    int fCurrentL;

    /* Physical bounds on the grid */
    double fMinX;
    double fMinY;
    double fMinE;
    double fMinT;
    double fMaxX;
    double fMaxY;
    double fMaxE;
    double fMaxT;

    bool fFullGrid;

    /* Rate of being killed/predator giving up, stored as a function pointer.
    std::function<double(double, double, double, double)> fMuK;
    std::function<double(double, double, double, double)> fMuG; */

    // mu_k, mu_g stored as pointers to 2d arrays of (x,y)
    // std::shared_ptr<memory::array2D_t<double>> fMuK;
    // std::shared_ptr<memory::array2D_t<double>> fMuG;

    /* Speed function, stored as a function pointer.
        This has to be given as an input.
        Computes speed f at physical location (x,y,e,t).
        Signature is fSpeed(x,y,e,t).
    std::function<double(double, double, double, double)> fSpeed; */

    // speed stored as a pointer to a 2d array of (x,y)
    std::shared_ptr<memory::array2D_t<double>> fSpeed;

    /* Cost function, stored as a function pointer.
         This has to be given as an input.
         Computes cost K at physical location (x,y,e,t).
         Signature is fCost(x,y,e,t).
    std::function<double(double, double, double, double)> fCost; */

    // cost stored as a pointer to a 2d array of (x,y)
    std::shared_ptr<memory::array2D_t<double>> fCost;

    

    // predator density as a pointer to a 2d array of (x,y)
    // std::shared_ptr<memory::array2D_t<double>> fPredatorDensity;

  public:
    /* ========================================================================
    *    Constructors
    * ========================================================================*/
    /** Default constructor */
    CTimeGrid() = default;

    /*
     * This is the main CTimeGrid constructor.
     * It initializes the grid and computes the grid spacing
     * for a square grid with aN points in both spatial dimensions,
     * aNe points in the energy dimension,
     * and aNt points in the time dimension.
     * @param aCost cost function pointer.
     * @param aSpeed speed function pointer.
     * @param aEnergyThreshold energy threshold for terminal condition.
     * @param aN an integer. The number of spatial grid points. Defaults to 101.
     * @param aNe an integer. The number of energy grid points. Defaults to 101.
     * @param aNt an integer. The number of time steps. Defaults to 1001.
     * @param aPhysMin a real number. The lower bound on the physical grid coordinates. Defaults to 0.
     * @param aPhysMax a real number. The upper bound on the physical grid coordinates. Defaults to 1.
     * @param aEnergyMin a real number. The lower bound on the energy grid. Defaults to 0.
     * @param aEnergyMax a real number. The upper bound on the energy grid. Defaults to 1.
     * @param aTime a real number. The length of the time interval [0, aTime]. Defaults to 1.
     * The domain is
     * [aPhysMin,aPhysMax] x [aPhysMin,aPhysMax] x [aEnergyMin, aEnergyMax] x [0,aTime].
     */

    // Function-based constructor
    CTimeGrid(std::function<double(double, double, double, double)> aCost,
              std::function<double(double, double, double, double)> aSpeed,
              const int aN = 101, const int aNe = 101, const int aNt = 1001,
              const double aPhysMin = 0, const double aPhysMax = 1,
              const double aEnergyMin = 0, const double aEnergyMax = 1,
              const double aTime = 1, const bool aFullGrid = true, 
              const int aCurrentL = 0);

    // Array-based constructor
    CTimeGrid(std::shared_ptr<memory::array2D_t<double>> aCost,
              std::shared_ptr<memory::array2D_t<double>> aSpeed,
              const int aN = 101, const int aNe = 101, const int aNt = 1001,
              const double aPhysMin = 0, const double aPhysMax = 1,
              const double aEnergyMin = 0, const double aEnergyMax = 1,
              const double aTime = 1, const bool aFullGrid = true, 
              const int aCurrentL = 0);

    /* ========================================================================
    *    Setters
    *=========================================================================*/
    /* Set the value function array at position (aI,aJ,aK,aL) to aValue */
    void setValue(const int aI, const int aJ, const int aK, const int aL,
                  const double aValue);

    void setCurrentSlice(const int aSlice);

    void setFilename(const std::string aFilename);

    /* ========================================================================
    *    Getters
    * ========================================================================*/
    /* Grid-coordinate functions */
    double getValue(const int aI, const int aJ, const int aK, const int aL) const;
    double getSpeed(const int aI, const int aJ, const int aK, const int aL) const;
    double getCost(const int aI, const int aJ, const int aK, const int aL) const;
  
    /* Physical-coordinate functions */
    double getValuePhysical(const double aX, const double aY, const double aE,
                            const double aT) const;
    double getSpeedPhysical(const double aX, const double aY, const double aE,
                            const double aT) const;
    double getCostPhysical(const double aX, const double aY, const double aE,
                            const double aT) const;
 
    /* Grid parameters */
    int getGridSizeX() const;
    int getGridSizeY() const;
    int getGridSizeE() const;
    int getGridSizeT() const;
    double getH() const;
    double getE() const;
    double getDt() const;
    double getMinX() const;
    double getMinY() const;
    double getMinE() const;
    double getMinT() const;
    double getMaxX() const;
    double getMaxY() const;
    double getMaxE() const;
    double getMaxT() const;
    // double getEnergyThreshold() const;

    /* ========================================================================
    *    Other
    * ========================================================================*/
    /* Mapping back and forth between grid and physical */
    double xGridToPhysical(const int aI) const;
    double yGridToPhysical(const int aJ) const;
    double eGridToPhysical(const int aK) const;
    double tGridToPhysical(const int aL) const;
    int xPhysicalToGrid(const double aX) const;
    int yPhysicalToGrid(const double aY) const;
    int ePhysicalToGrid(const double aE) const;
    int tPhysicalToGrid(const double aT) const;

    /*
     * A I/O member which prints the value and cost grids to file.
     * @param aFilename a string which contains the prefix to the names
     *      of the files to which the grids will be printed.
     * The (cost/value) grids will be printed to
     *      files called "aFilename"+(Cost/Value)
     */
    void writeGridToFile(const std::string aFilename) const;
    void writeSliceToFile(const std::string aFilename, const int aSlice) const;

    void advanceSlicesBackward();
    void advanceSlicesForward();

    void readValuesFromFile(const std::string aFilename, const int aSlice = 0);   
    void readValuesFromFile(const int aSlice);
    bool isFullGrid() const;

};
/* ============================================================================
*    Inline function definitions
* ============================================================================*/
/* ------ Inline definition of setter function ------------------------------*/
inline void CTimeGrid::setValue(const int aI, const int aJ, const int aK,
                                const int aL, const double aValue) {
  if (fFullGrid) {
    (*fValues)[aL][aI][aJ][aK] = aValue;
  } else {
    if (aL == fCurrentL) {
      (*fValuesCurrent)[aI][aJ][aK] = aValue;
    } else if (aL == fCurrentL + 1) {
      (*fValuesFuture)[aI][aJ][aK] = aValue;
    } else {
      std::cout << "Invalid time slice entered." << std::endl;
      std::cout << "Entered " << aL << " when current L is " << fCurrentL << std::endl;
      assert(false);
    }
  }
}

inline void CTimeGrid::setCurrentSlice(const int aSlice) {
  fCurrentL = aSlice;
}

inline void CTimeGrid::setFilename(const std::string aFilename) {
  fFilename = aFilename;
}


/* ------ Inline definition of getters --------------------------------------*/
/* Get from grid indices */
inline double CTimeGrid::getValue(const int aI, const int aJ,
                                  const int aK, const int aL) const {
  double returnValue;
  if (fFullGrid) {
    returnValue = (*fValues)[aL][aI][aJ][aK];
  } else {
    if (aL == fCurrentL) {
      returnValue = (*fValuesCurrent)[aI][aJ][aK];
    } else if (aL == fCurrentL + 1) {
      returnValue = (*fValuesFuture)[aI][aJ][aK];
    } else {
      std::cout << "Invalid time slice requested." << std::endl;
      std::cout << "Requested " << aL << " when current L is " << fCurrentL << std::endl;
      assert(false);
    }
  }
  return returnValue;
}

inline double CTimeGrid::getSpeed(const int aI, const int aJ,
                                  const int aK, const int aL) const {
  return (*fSpeed)[aI][aJ];
}

inline double CTimeGrid::getCost(const int aI, const int aJ,
                                 const int aK, const int aL) const {
  return (*fCost)[aI][aJ];
}

/* Get grid sizes (and rates mu) */
inline int CTimeGrid::getGridSizeX() const{
  return round((fMaxX - fMinX)/fH + 1);
}
inline int CTimeGrid::getGridSizeY() const{
  return round((fMaxY - fMinY)/fH + 1);
}
inline int CTimeGrid::getGridSizeE() const{
  return round((fMaxE - fMinE)/fE + 1);
}
inline int CTimeGrid::getGridSizeT() const{
  return round((fMaxT - fMinT)/fDt + 1);
}
inline double CTimeGrid::getH() const{
  return fH;
}
inline double CTimeGrid::getE() const{
  return fE;
}
inline double CTimeGrid::getDt() const{
  return fDt;
}
inline double CTimeGrid::getMinX() const{
  return fMinX;
}
inline double CTimeGrid::getMinY() const{
  return fMinY;
}
inline double CTimeGrid::getMinE() const{
  return fMinE;
}
inline double CTimeGrid::getMinT() const{
  return fMinT;
}
inline double CTimeGrid::getMaxX() const{
  return fMaxX;
}
inline double CTimeGrid::getMaxY() const{
  return fMaxY;
}
inline double CTimeGrid::getMaxE() const{
  return fMaxE;
}
inline double CTimeGrid::getMaxT() const{
  return fMaxT;
}

inline bool CTimeGrid::isFullGrid() const{
  return fFullGrid;
}

/** ------ Inline definition of conversion functions -------------------------*/
inline double CTimeGrid::xGridToPhysical(const int aI) const {
  return fMinX + (double)aI * fH;
}

inline double CTimeGrid::yGridToPhysical(const int aJ) const {
  return fMinY + (double)aJ * fH;
}

inline double CTimeGrid::eGridToPhysical(const int aK) const {
  return fMinE + (double)aK * fE;
}

inline double CTimeGrid::tGridToPhysical(const int aL) const {
  return fMinT + (double)aL * fDt;
}

inline int CTimeGrid::xPhysicalToGrid(const double aX) const {
  return std::round((aX - fMinX)/fH);
}

inline int CTimeGrid::yPhysicalToGrid(const double aY) const {
  return std::round((aY - fMinY)/fH);
}

inline int CTimeGrid::ePhysicalToGrid(const double aE) const {
  return std::round((aE - fMinE)/fE);
}

inline int CTimeGrid::tPhysicalToGrid(const double aT) const {
  return std::round((aT - fMinT)/fDt);
}

#endif
