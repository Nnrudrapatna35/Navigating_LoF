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
*  along with this program. If not, see <http://www.gnu.org/licenses/>.
*
* ------------------------------------------------------------------------------
* File: GlobalConfigurations.hpp
*
* Authors: Marissa Gee Nicolas Gonzalez Granda Sunay Joshi 
*          Nagaprasad Rudrapatna Anne Somalwar 
*          (based on code by Elliot Cartee, Qianli Song, Lexiao Lai)
*             (based on code by Marc Aurèle Gilles and Zachary Clawson)
*
* Description: This file contains the definitions of various numerical constants
* and types. It also contains all global parameters for all of the numerics.
*
* ==============================================================================
*/

#ifndef GLOBAL_CONFIG_HPP
#define GLOBAL_CONFIG_HPP

/*-- STL ---------------------------------------------------------------------*/
#include <cassert>
#include <functional>
#include <limits>

/*-----------------------------------------------------------------------------/
/-- Numerical Constants
/-----------------------------------------------------------------------------*/
constexpr double INF = std::numeric_limits<double>::max();
constexpr int INT_INF = std::numeric_limits<int>::max();
constexpr double NEG_INF = std::numeric_limits<double>::lowest();
constexpr double PI = 3.141592653589793;
constexpr double SQRT2 = 1.41421356237309504880168872420969807;
constexpr double LARGE_NUMBER = 1.0e32;
constexpr double LARGE_NEGATIVE_NUMBER = -1.0e32;
constexpr double EPSILON = 1e-16;
constexpr float EPSILON_FLOAT = 1e-8;

constexpr double IMPROVEMENT_THRESHOLD = 1e-7;
constexpr double GIVE_UP_THRESHOLD = 1e-7;


/*-----------------------------------------------------------------------------/
/-- Gridpoint status for FMM
/-----------------------------------------------------------------------------*/
enum status_t {FAR, CONSIDERED, ACCEPTED};
enum objective_t {SURVIVAL, ENERGY, UTILITY, SQUAREROOT, 
                  SIGMOID, QUADRATIC, LINEAR};
enum stencil_t {FIVE_POINT, NINE_POINT};
enum life_status_t : char {ALIVE, DEAD, LOWENERGY, INVALID};
enum patch_t {HOME, UPPER, LOWER, BOTH};

/* ------ Shorter version of function<float(float)> ------------------------*/
typedef float (*Function1D)(float);

/*-----------------------------------------------------------------------------/
/-- Global parameters & flags
/-----------------------------------------------------------------------------*/
/* Terminal conditions at all non-target points */
constexpr double TERMINAL_VALUE = 1;

/* Tolerance in outer loop iterations */
constexpr double TOLERANCE = 1e-6;

/* ------ Parameters for trajectory tracing ---------------------------------*/
/* N_THETA is the number of directions to search in.
 * Should be divisible by 4 so that the 4 grid direction are searched. */
  const int N_THETA = 360;

#endif
