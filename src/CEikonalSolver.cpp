/*
* ==============================================================================
*
*  Copyright (C) 2021  Marissa Gee
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
* File: CEikonalSolver.cpp
*
* Author: Marissa Gee
*   (based on code by Elliot Cartee, Marc Aurèle Gilles, and Zachary Clawson)
**
* Description: This is a class for solving stationary Eikonal equations using
* the Fast Marching Method.
* It works on a 2D regular grid and a 5 point stencil.
* It allows for different spacing in the horizontal and vertical directions.
* The target set is determined by aDomain.
* The solution to the Eikonal is stored in the fPrimary grid.
* The FMM implementation uses the boost::heap and CGrid classes.
* It contains an augmented grid to compute the evolution of the boundary
*
* (See also CEikonalSolver.hpp)
*
* ==============================================================================
*/
/** ------ Libraries ---------------------------------------------------------*/
#include <cmath>
#include <iostream>
#include <cassert>

/** ------ Project-specific header files -------------------------------------*/
#include "CEikonalSolver.hpp"

/** ------ Namespaces --------------------------------------------------------*/
using namespace std;
using namespace memory;

/** ============================================================================
  Eikonal Constructor (RepairToYou)
==============================================================================*/

/* Array-based constructor */
CEikonalSolver::CEikonalSolver(const shared_ptr<array2D_t<double>> aSpeed,
                               const shared_ptr<CTerrain> aTerrain,
                               const int aN, const double aPhysMin, const double aPhysMax)
                               : fN(aN), fMinX(aPhysMin), fMinY(aPhysMin),
                                fMaxX(aPhysMax), fMaxY(aPhysMax) {
  /** Initialize properties and allocate arrays. */

  /** In this scenario the value function at the boundary set is always equal to
  *   zero
  */
  const int nx = fN;
  const int ny = fN;
  fH = (fMaxX - fMinX)/(fN-1);
  fBoundaryValue = make_shared<array2D_t<double>>(allocateArray2D<double>(nx, ny));
  fPrimary = make_shared<array2D_t<double>>(allocateArray2D<double>(nx, ny));
  fSpeed = make_shared<array2D_t<double>>(allocateArray2D<double>(nx, ny));
  fDomain = make_shared<array2D_t<bool>>(allocateArray2D<bool>(nx, ny));
  fBoundarySet = make_shared<array2D_t<bool>>(allocateArray2D<bool>(nx, ny));

  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      (*fBoundaryValue)[i][j] = 0;
      (*fDomain)[i][j] = true;
      if (aTerrain->isHomeBase(i,j)) {
        (*fBoundarySet)[i][j] = true;
      } else {
        (*fBoundarySet)[i][j] = false;
      }

      if (aTerrain->isObstacle(i,j)) {
        (*fSpeed)[i][j] = 0;
      } else {
        (*fSpeed)[i][j] = (*aSpeed)[i][j];
      }
    }
  }
}


/* Function-based constructor */
CEikonalSolver::CEikonalSolver(const function<double(double, double, double, double)> aSpeed,
                               const shared_ptr<CTerrain> aTerrain,
                               const int aN, const double aPhysMin, const double aPhysMax)
                               : fN(aN), fMinX(aPhysMin), fMinY(aPhysMin),
                                fMaxX(aPhysMax), fMaxY(aPhysMax) {
  /** Initialize properties and allocate arrays. */

  /** In this scenario the value function at the boundary set is always equal to
  *   zero
  */
  const int nx = fN;
  const int ny = fN;
  fH = (fMaxX - fMinX)/(fN-1);
  fBoundaryValue = make_shared<array2D_t<double>>(allocateArray2D<double>(nx, ny));
  fPrimary = make_shared<array2D_t<double>>(allocateArray2D<double>(nx, ny));
  fSpeed = make_shared<array2D_t<double>>(allocateArray2D<double>(nx, ny));
  fDomain = make_shared<array2D_t<bool>>(allocateArray2D<bool>(nx, ny));
  fBoundarySet = make_shared<array2D_t<bool>>(allocateArray2D<bool>(nx, ny));

  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      // compute x, y
      // NOTE (aN-1)!
      double x = aPhysMin + (double)i / (aN-1) * (aPhysMax - aPhysMin);
      double y = aPhysMin + (double)j / (aN-1) * (aPhysMax - aPhysMin);

      (*fBoundaryValue)[i][j] = 0;
      (*fDomain)[i][j] = true;

      // function call
      if (aTerrain->isHomeBasePhysical(x,y)) {
        (*fBoundarySet)[i][j] = true;
      } else {
        (*fBoundarySet)[i][j] = false;
      }

      // function call
      if (aTerrain->isObstaclePhysical(x,y)) {
        (*fSpeed)[i][j] = 0;
      } else {
        (*fSpeed)[i][j] = aSpeed(x,y,0,0); // function call, with DUMMY values for last two args
      }
    }
  }
}


/** ============================================================================
  Main compute functions
==============================================================================*/
std::shared_ptr<memory::array2D_t<double>> CEikonalSolver::march() {
  /** Initialize heap */
  CFMMHeap_t CFMMHeap;
  initialize_cfmm(CFMMHeap);

  /** Perform Djikstra's algorithm to update all grid points */
  while (!CFMMHeap.empty()) {
    CHeapGP current_GP = CFMMHeap.top();
    CFMMHeap.pop();
    (*fStatus)[current_GP.fI][current_GP.fJ] = ACCEPTED;
    update_neighbors(CFMMHeap, current_GP.fI, current_GP.fJ);
  }

  shared_ptr<array2D_t<double>> returnArray = make_shared<array2D_t<double>>(allocateArray2D<double>(fN,fN));

  for (int i = 0; i < fN; i++) {
    for (int j = 0; j < fN; j++) {
      (*returnArray)[i][j] = double((*fPrimary)[i][j]);
    }
  }
  return returnArray;
}

void CEikonalSolver::initialize_cfmm(CFMMHeap_t& aCFMMHeap) {
  const int nx = fN;
  const int ny = fN;
  fStatus = make_shared<array2D_t<status_t>>(allocateArray2D<status_t>(nx,ny));
  fHeapPointers = make_shared<array2D_t<handle_t>>(allocateArray2D<handle_t>(nx,ny));

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      if (!in_domain(i,j)) {
        /**
        * All points outside the domain are set to ACCEPTED and have infinite cost.
        */
        (*fPrimary)[i][j] = INF;
        (*fStatus)[i][j] = ACCEPTED;
      } else {
        /** If speed is 0, then value has to be infinity there. */
        const double f = double((*fSpeed)[i][j]);
        if (pow(f,2) == 0) {
          (*fPrimary)[i][j] = INF;
          (*fStatus)[i][j] = ACCEPTED;
        } else {
          /** Initialize gridpoints at the goal set */
          if ((*fBoundarySet)[i][j] == true) {
            const double tempVal = (*fBoundaryValue)[i][j];
            (*fPrimary)[i][j] = tempVal;
            (*fHeapPointers)[i][j] = aCFMMHeap.push(CHeapGP(i, j, tempVal));
            (*fStatus)[i][j] = CONSIDERED;
          } else {
            (*fPrimary)[i][j] = LARGE_NUMBER;
            (*fStatus)[i][j] = FAR;
          }
        }
      }
    }
  }
}


/** ============================================================================
  Helper functions
==============================================================================*/
void CEikonalSolver::update_neighbors(CFMMHeap_t& aCFMMHeap, const int aCurrent_i,
                           const int aCurrent_j) {
  /* Five-point Stencil */
  constexpr int stencil[4][2] = {{1,0}, {-1,0}, {0,1}, {0,-1}};

  /* Iterate over neighbors. */
  for (int k = 0; k < 4; ++ k) {
    const int nbr_i = aCurrent_i + stencil[k][0];
    const int nbr_j = aCurrent_j + stencil[k][1];

    if (in_domain(nbr_i,nbr_j)) {
      /* If accepted nothing to do, otherwise, try to update */
      if ((*fStatus)[nbr_i][nbr_j] != ACCEPTED) {
        const bool gp_was_updated = update_gp(nbr_i, nbr_j);

        /* If it was already considered and was updated then update heap */
        if ((*fStatus)[nbr_i][nbr_j] == CONSIDERED) {
          if (gp_was_updated) {
            CHeapGP gp = CHeapGP(nbr_i, nbr_j, (*fPrimary)[nbr_i][nbr_j]);
            aCFMMHeap.update((*fHeapPointers)[nbr_i][nbr_j], gp);
          }
        } else {
          /* Else add to heap */
          (*fStatus)[nbr_i][nbr_j] = CONSIDERED;
          CHeapGP gp = CHeapGP(nbr_i, nbr_j, (*fPrimary)[nbr_i][nbr_j]);
          (*fHeapPointers)[nbr_i][nbr_j] = aCFMMHeap.push(gp);
        }
      }
    }
  }
}

int CEikonalSolver::smaller_h_neighbor(const int aI, const int aJ) const {
  int bestI;
  if (aI == 0) {
    bestI = 1;
  } else if (aI == fN - 1) {
    bestI = -1;
  } else if ((*fPrimary)[aI-1][aJ] < (*fPrimary)[aI+1][aJ]) {
    bestI = -1;
  } else {
    bestI = 1;
  }

  return bestI;
}

int CEikonalSolver::smaller_v_neighbor(const int aI, const int aJ) const {
  int bestJ;
  if (aJ == 0) {
    bestJ = 1;
  } else if (aJ == fN - 1) {
    bestJ =  -1;
  } else if ((*fPrimary)[aI][aJ-1] < (*fPrimary)[aI][aJ+1]) {
    bestJ = -1;
  } else {
    bestJ = 1;
  }

  return bestJ;
}

/** ============================================================================
  Update functions
==============================================================================*/
double CEikonalSolver::compute_from_two_neighbors(const int aI, const int aJ,
                                      const int aVerticalNeighbor,
                                      const int aHorizontalNeighbor) {
  /** Allocate return value */
  double new_value;

  /** Get stepsizes */
  const double dx = fH;
  const double dy = fH;

  /** Get shorter names for local variables */
  const double f = double((*fSpeed)[aI][aJ]);
  const double K = 1;

  const double u_h = (*fPrimary)[aI + aHorizontalNeighbor][aJ];
  const double u_v = (*fPrimary)[aI][aJ + aVerticalNeighbor];


  /** Variables containing possible updates to primary variable */
  /** Signature is {quad1, quad2, horizontal, vertical} */
  std::vector<double> updates = {LARGE_NUMBER, LARGE_NUMBER, LARGE_NUMBER, LARGE_NUMBER};
  std::vector<bool> upwinding = {false, false, false, false};

  /** Compute discriminant */
  double A, B, C;
  A = pow(f/dx,2) + pow(f/dy,2);
  B = -2*u_h*pow(f/dx,2) - 2*u_v*pow(f/dy,2);
  C = pow(f*u_h/dx,2) + pow(f*u_v/dy,2) - pow(K,2);
  const double discriminant = B*B - 4*A*C;

  if (discriminant >= 0) {
    /** Compute two-sided update */
    if (A == 0) {
      /** Degenerate case: Quadratic update is actually just a linear update */
      if (B == 0) {
        updates[0] = INF;
        updates[1] = INF;
      } else {
        updates[0] = -C / B;
        updates[1] = updates[0];
      }
    } else {
      /** Quadratic formula */
      updates[0] = (-B + sqrt(discriminant)) / (2*A);
      updates[1] = (-B - sqrt(discriminant)) / (2*A);
    }

    /** Check upwinding condition for both two-sided updates */
    if (updates[0] >= max(u_h, u_v)) {
      upwinding[0] = true;
    }
    if (updates[1] >= max(u_h, u_v)) {
      upwinding[1] = true;
    }
  }

  /** Compute one sided updates. */
  updates[2] = u_h + dx * K / f;
  updates[3] = u_v + dy * K / f;

  /** Check upwinding condition for one-sided updates. */
  upwinding[2] = (updates[2] >= u_h);
  upwinding[3] = (updates[3] >= u_v);

  /** Find smallest upwind update */
  int min_upwind = -1;
  double min_update = INF;
  for (int i = 0; i < 4; i++) {
    if (upwinding[i] && (updates[i] < min_update)) {
      min_upwind = i;
      min_update = updates[i];
    }
  }
  if (min_upwind == -1) {
    /** No upwind updates */
    assert(false);
    new_value = LARGE_NUMBER;
  }

  /** Use smallest upwind update */
  new_value = updates[min_upwind];
  return new_value;
}

/*-----------------------------------------------------------------------------/
/-- Update gridpoint
/-----------------------------------------------------------------------------*/
bool CEikonalSolver::update_gp(const int aI, const int aJ) {
  /** Compute smaller vertical and horizontal neighbors. */
  const int hNeighbor = smaller_h_neighbor(aI,aJ);
  const int vNeighbor = smaller_v_neighbor(aI,aJ);

  /** Compute update based on smaller horizontal and vertical neighbors */
  const double new_value = compute_from_two_neighbors(aI, aJ, vNeighbor, hNeighbor);
  if (new_value < (*fPrimary)[aI][aJ]) {
    (*fPrimary)[aI][aJ] = new_value;
    return true;
  } else {
    return false;
  }
}
