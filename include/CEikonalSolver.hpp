/*
* ==============================================================================
*  Copyright (C) 2024 Marissa Gee Nicolas Gonzalez Granda Sunay Joshi 
*                     Nagaprasad Rudrapatna Anne Somalwar 
*  Copyright (C) 2021 Marissa Gee
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
*
* File: CEikonalSolver.hpp
*
* Author: Marissa Gee
*   (based on code by Elliot Cartee, Marc Aurèle Gilles, and Zachary Clawson)
*
* Description: This is a class for solving stationary PDEs
*              using the Fast Marching Method.
* It works on a 2D regular grid and a 5 point stencil.
* It allows for different spacing in the horizontal and vertical directions.
* The target set is determined by aDomain.
* The solution to the PDE is stored in the fPrimary grid.
* The FMM implementation uses the boost::heap and CGrid classes.
* It contains an augmented grid to compute the evolution of the boundary
*
* ==============================================================================
*/

#ifndef CEIKONALSOLVER_HPP
#define CEIKONALSOLVER_HPP

/** ------ Libraries ---------------------------------------------------------*/
#include <boost/heap/binomial_heap.hpp>

/** ------ Project-specific header files -------------------------------------*/
#include "GlobalConfigurations.hpp"
#include "MemoryAllocations.hpp"
#include "CTerrain.hpp"

class CEikonalSolver
{
  public:
    /** ========================================================================
      Constructors
    ==========================================================================*/
    /** Default Constructor */
    CEikonalSolver() = default;

    /**
     * The Eikonal CEikonalSolver Constructor (RepairToYou)
     * @param aPrimary pointer to CGrid object containing the primary grid
     *           on which the Eikonal equation is solved
     * @param aDomain  pointer to array defining domain in which to solve
     *           Eikonal equation
     * @param aBoundarySet  pointer to array defining location of the goal set,
     *           assumes value of zero at boundary
     */

     // array-based constructor
    CEikonalSolver(const std::shared_ptr<memory::array2D_t<double>> aSpeed,
                   const std::shared_ptr<CTerrain> aTerrain,
                   const int aN, const double aPhysMin, const double aPhysMax);

     // function-based constructor
    CEikonalSolver(const std::function<double(double, double, double, double)> aSpeed,
                   const std::shared_ptr<CTerrain> aTerrain,
                   const int aN, const double aPhysMin, const double aPhysMax);



    /** ========================================================================
      Compute functions
    ==========================================================================*/
    /**
     * This function computes the PDE solution using the Fast Marching Method.
     *
     * The solution is stored in fPrimary (accessible by getValue/setValue).
     */
     std::shared_ptr<memory::array2D_t<double>> march();

  private:
    /** A GridPoint Class to be used by the heap. */
    class CHeapGP {
      public:
        int fI;
        int fJ;
        double fValue;
        CHeapGP(int aI, int aJ, double aValue): fI(aI), fJ(aJ), fValue(aValue) {};
    };

    /** Struct comparison which is required by boost::heap */
    struct compare_CHeapGP
    {
      bool operator()(const CHeapGP& aPoint1, const CHeapGP& aPoint2) const
      {
        return aPoint1.fValue > aPoint2.fValue;
      }
    };

    double fH;
    int fN;

    /* Physical bounds on the grid */
    double fMinX;
    double fMinY;
    double fMaxX;
    double fMaxY;

    /** Typedef Heap types to make it easier to read. */
    typedef boost::heap::binomial_heap<CHeapGP, boost::heap::compare<compare_CHeapGP> > CFMMHeap_t;
    typedef typename boost::heap::binomial_heap<CHeapGP, boost::heap::compare<compare_CHeapGP> >::handle_type handle_t;

    /** A shared pointer to a CGrid for the solution of the Eikonal equation */
    std::shared_ptr<memory::array2D_t<double>> fPrimary;

    /** Array representing Domain */
    std::shared_ptr<memory::array2D_t<bool>> fDomain;
    /** Array representing the locations the of the goal set */
    std::shared_ptr<memory::array2D_t<bool>> fBoundarySet;
    /** Array representing the value function at the boundary set*/
    std::shared_ptr<memory::array2D_t<double>> fBoundaryValue;
    /** Array representing the value function at the boundary set*/
    std::shared_ptr<memory::array2D_t<double>> fSpeed;


    /** The status (far/considered/accepted) of each grid point. */
    std::shared_ptr<memory::array2D_t<status_t>> fStatus;
    /** Backpointers for the heap */
    std::shared_ptr<memory::array2D_t<handle_t>> fHeapPointers;

    /**
     * FMM Initialization.
     * This function sets the status of all nodes
     *      and adds the border to the heap.
     * @param CFMMHeap_t& aCFMMHeap Boost::heap passed by value
     *      it is initialized by this function.
     */
    void initialize_cfmm(CFMMHeap_t& aCFMMHeap);

    /**
    * A helper function which updates neighbors and adds them to the heap.
    * @param aCFMMHeap boost heap object
    * @param aCurrent_i int x logical coordinate of grid point
    * @param aCurrent_j int y logical coordinate of grid point
    */
    void update_neighbors(CFMMHeap_t& aCFMMheap,
                          const int aCurrent_i, const int aCurrent_j);

    /**
     * A helper function to figure out if grid point (aI,aJ) is in the domain.
     * @param aI int x logical coordinate
     * @param aJ int y logical coordinate
     */
    bool in_domain(const int aI, const int aJ) const;

    /** Returns the value of the smaller of the two horizontal neighbors.
    * @param aI int x logical coordinate
    * @param aJ int y logical coordinate
    */
    int smaller_h_neighbor(const int aI, const int aJ) const;

    /** Returns the value of the smaller of the two vertical neighbors.
    * @param aI int x logical coordinate
    * @param aJ int y logical coordinate
    */
    int smaller_v_neighbor(const int aI, const int aJ) const;

    /** Computes the Eikonal update.
    * @param aI int x logical coordinate
    * @param aJ int y logical coordinate
    * @param aVerticalNeighbor int direction of smaller vertical neighbor
    * @param aHorizontalNeighbor int direction of smaller horizontal neighbor
    */
    double compute_from_two_neighbors(const int aI, const int aJ,
                                      const int aVerticalNeighbor,
                                      const int aHorizontalNeighbor);

    /** Update a gridpoint. Compute the update, assign, and return if updated.
    * @param aI int x logical coordinate
    * @param aJ int y logical coordinate
    */
    bool update_gp(const int aI, const int aJ);
};

/** ============================================================================
  Inline function definitions
==============================================================================*/

/* Compute if gridpoint [aI,aJ] is inside the domain */
inline
bool CEikonalSolver::in_domain(const int aI, const int aJ) const {
  if (aI >= 0 && aI < fN && aJ >=0 && aJ < fN) {
    return (*fDomain)[aI][aJ];
  }
  return false;
}

#endif
