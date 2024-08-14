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
* File: SimpleFunctions.hpp
*
* Authors: Marissa Gee Nicolas Gonzalez Granda Sunay Joshi 
*          Nagaprasad Rudrapatna Anne Somalwar 
*          (based on code by Elliot Cartee, Qianli Song, Lexiao Lai)
*             (based on code by Marc Aurèle Gilles and Zachary Clawson)
*
* Description: This file contains inline definitions of many small/simple
* functions used for speed or pointwise observability.
*
* ==============================================================================
*/

#ifndef SIMPLEFUNCTIONS_HPP
#define SIMPLEFUNCTIONS_HPP

#include <cmath>

#include "GlobalConfigurations.hpp"

/*****************************************************************************
 * Constant Functions
*****************************************************************************/
inline double varConstant(double aX, double aC) {
  return aC;
}

inline double varConstant2D(double aX, double aY, double aC) {
  return aC;
}

inline double varConstant3D(double aX, double aY, double aE, double aC) {
  return aC;
}

inline double varConstant4D(double aX, double aY, double aE, double aT, double aC) {
  return aC;
}


/*****************************************************************************
 * Basic Gaussian Functions
*****************************************************************************/
inline double varGaussian1D(double aX, double aY, double aMuY, double aSigma, 
                            double aScale, double aOffset) {
  return aOffset + aScale*exp(-((aY-aMuY)*(aY-aMuY)) / (2*aSigma*aSigma));
}

inline double varGaussian2D(double aX, double aY, double aMuX, double aMuY,
                            double aSigma, double aScale, double aOffset) {
  return aOffset + aScale*exp(-((aX-aMuX)*(aX-aMuX) + (aY-aMuY)*(aY-aMuY)) / (2*aSigma*aSigma));
}

inline double varEllipseGaussian2D(double aX, double aY, double aMuX, double aMuY,
                                   double aSigmaX, double aSigmaY, 
                                   double aScale, double aOffset) {
  return aOffset + aScale*exp(-((aX-aMuX)*(aX-aMuX)/(2*aSigmaX*aSigmaX) + (aY-aMuY)*(aY-aMuY)/(2*aSigmaY*aSigmaY)) );
}

/*****************************************************************************
 * Bump Functions
*****************************************************************************/

inline double square_region_center_bump(double x, double y) {
  const double centerX = 0.5;
  const double centerY = 0.5;
  const double widthX = 0.05;
  const double widthY = 0.05;
  const double scale = 0.4;
  double value = 0;
  if (fabs(x-centerX) >= widthX) {
    value = 0;
  } else if (fabs(y-centerY) >= widthY) {
    value = 0;
  } else {
    value = exp(-scale/(1 - pow((x - centerX)/widthX, 2))) 
            * exp(-scale/(1 - pow((y - centerY)/widthY, 2)));
  }

  return 1 - value/(exp(-scale)*exp(-scale));
}


/*****************************************************************************
 * Kernel Functions
*****************************************************************************/
/* Single-variable gaussian which integrates to 1 */
inline double small_gauss_kernel(double aX, double aR) {
  double sigma = aR/4;
  if (fabs(aX) > aR) {
    return 0;
  }
  return exp(-aX*aX / (2*pow(sigma, 2))) / (sqrt(2*PI)*sigma);
}

inline double constant_kernel(double aX, double aR) {
  if (aX > aR || aX < 0) {
    return 0;
  }
  return 1/(PI*aR*aR);
}

inline double small_gauss_kernel(double aX, double aR) {
  if (aX > aR || aX < 0) {
    return 0;
  }
  return 1/(PI*aR*aR) * exp(1)/(exp(1)-1) * exp(-aX*aX/(aR*aR));
}

/****************************************************************************
 ************************ Food Density Functions ****************************
 ****************************************************************************/

inline double twoPeaksClose(double aX, double aY) {
  const double baseline = 0;

  const double muY1 = 0.1;
  const double muY2 = 0.5;
  const double muX1 = 0.4;
  const double muX2 = 0.4;
  const double sigma = 0.035;
  const double offset = 0;
  const double scale1 = 500;
  const double scale2 = scale1;
  const double gaussian1 = varGaussian2D(aX, aY, muX1, muY1, sigma, scale1, offset);
  const double gaussian2 = varGaussian2D(aX, aY, muX2, muY2, sigma, scale2, offset);

  return baseline + gaussian1 + gaussian2;
}

/* Risk-reward Food Densities*/
inline double circleFood3Cts(double aX, double aY) {
  const double baseline = 0.0;

  const double muX = 0.5;
  const double muY1 = 0.1; // Smaller food density
  const double muY2 = 0.9; // Bigger food density
  const double sigma = 0.05; // Correct, or 0.1?
  const double offset = 0;
  const double scale1 = 6.0 / 0.001;
  const double scale2 = 10.0 / 0.001;
  const double gaussian1 = varGaussian2D(aX, aY, muX, muY1, sigma, scale1, offset);
  const double gaussian2 = varGaussian2D(aX, aY, muX, muY2, sigma, scale2, offset);

  return baseline + gaussian1 + gaussian2;
}

inline double circleFood4Cts(double aX, double aY) {
  const double baseline = 0.0;

  const double muX = 0.5;
  const double muY1 = 0.1; // Smaller food density
  const double muY2 = 0.9; // Bigger food density
  const double sigma = 0.05; // Correct, or 0.1?
  const double offset = 0;
  const double scale1 = 7.0 / 0.001;
  const double scale2 = 9.0 / 0.001;
  const double gaussian1 = varGaussian2D(aX, aY, muX, muY1, sigma, scale1, offset);
  const double gaussian2 = varGaussian2D(aX, aY, muX, muY2, sigma, scale2, offset);

  return baseline + gaussian1 + gaussian2;
}

/****************************************************************************
 ************************ Spotting Rate Functions ***************************
 ****************************************************************************/
/* Writing a new function to be used as BIG predator density*/

/* Risk-reward spotting rates */
inline double varCirclePred3Cts(double aX, double aY, double aE, double aT, double aUpperSpottingRate) {
  const double baseline = 0.5; // Replaced with constant baseline

  const double muX = 0.5;
  const double muY1 = 0.3; // Less dangerous
  const double muY2 = 0.7; // More dangerous
  const double sigma = 0.1; // Correct? Change back to 0.05?
  const double offset = 0;
  const double scale1 = 0.5;
  const double scale2 = scale1 + aUpperSpottingRate; // Slowly increase this number! Range from 0.5 to 1.3?
  const double gaussian1 = varGaussian2D(aX, aY, muX, muY1, sigma, scale1, offset);
  const double gaussian2 = varGaussian2D(aX, aY, muX, muY2, sigma, scale2, offset);

  return baseline + gaussian1 + gaussian2;
}

inline double varCirclePred4Cts(double aX, double aY, double aE, double aT, double aBaselineShift) {
  const double baseline = 0.1 + aBaselineShift; // Shift the baseline of 0.1

  const double muX = 0.5;
  const double muY1 = 0.3; // Less dangerous
  const double muY2 = 0.7; // More dangerous
  const double sigma = 0.1; 
  const double offset = 0;
  const double scale1 = 0.5;
  const double scale2 = 0.5 + aBaselineShift; // Shift the value of 0.5
  const double gaussian1 = varGaussian2D(aX, aY, muX, muY1, sigma, scale1, offset);
  const double gaussian2 = varGaussian2D(aX, aY, muX, muY2, sigma, scale2, offset);

  return baseline + gaussian1 + gaussian2;
}


/****************************************************************************
 ******************** Give Up and Kill Rate Functions ***********************
 ****************************************************************************/
inline double rectangle_give_up(double aX, double aY, double aE, double aT) {
  const double centerX = 0.4;
  const double centerY = 0.1;
  const double widthX = 0.15;
  const double widthY = 0.15;
  const double scale = 1;
  const double baseline = 2;
  const double height = 3;
  double value = 0;
  if (fabs(aX-centerX) >= widthX) {
    value = 0;
  } else if (fabs(aY-centerY) >= widthY) {
    value = 0;
  } else {
    value = exp(-scale/(1 - pow((aX - centerX)/widthX, 2))) 
            * exp(-scale/(1 - pow((aY - centerY)/widthY, 2)));
  }

  return baseline + height*value/(exp(-scale)*exp(-scale));
}

inline double rectangle_kill(double aX, double aY, double aE, double aT, double aHeight) {
  const double centerX = 0.4;
  const double centerY = 0.1;
  const double widthX = 0.15;
  const double widthY = 0.15;
  const double scale = 1;
  const double protection = 0.2;
  double value = 0;
  if (fabs(aX-centerX) >= widthX) {
    value = 0;
  } else if (fabs(aY-centerY) >= widthY) {
    value = 0;
  } else {
    value = exp(-scale/(1 - pow((aX - centerX)/widthX, 2))) 
            * exp(-scale/(1 - pow((aY - centerY)/widthY, 2)));
  }

  return aHeight - (protection*aHeight)*value/(exp(-scale)*exp(-scale));
}


/****************************************************************************
 *************************** Utility Functions ******************************
 ****************************************************************************/
inline double power_utility(double aE, double aTheta, double aMaxEnergy) {
  return pow((aE/aMaxEnergy), aTheta);
}

inline double threshold_utility(double aE, double aTheta) {
  double value;
  if (aE > aTheta) {
    value = 1;
  } else {
    value = 0;
  }
  return value;
}

inline double sigmoid_utility_atan(double aE, double aTheta, double aEThreshold,
                                   double aMaxEnergy) {
  return (atan(aTheta*(aE - aEThreshold)) + atan(aTheta*aEThreshold))
         / (atan(aTheta*(aMaxEnergy - aEThreshold)) + atan(aTheta*aEThreshold));
}

inline double sigmoid_utility_normcdf(double aE, double aTheta, double aEThreshold,
                                      double aMaxEnergy) {
  return (erf((aE - aEThreshold)/(sqrt(2)*aTheta)) - erf((-aEThreshold)/(sqrt(2)*aTheta)))
         / (erf((aMaxEnergy - aEThreshold)/(sqrt(2)*aTheta)) - erf((-aEThreshold)/(sqrt(2)*aTheta)));
}


/****************************************************************************
 ********************* Home Base & Obstacle Functions ***********************
 ****************************************************************************/
inline double none(double x, double y) {
  return 1;
}

inline double bar_region(double aX, double aY) {
  double value;
  if (aY <= 0.1) {
    value = 0;
  } else {
    value = 1;
  }
  return value;
}

inline double risk_reward_hb(double aX, double aY) {
	if (aX >= 0.4 && aX <= 0.5 && aY >= 0.45 && aY <= 0.55) {
    return 0;
  } else {
    return 1;
  }
}

inline double square_region(double x, double y) {
	if (x >= 0.55 && x <= 0.65 && y >= 0.55 && y <= 0.65) {
    return 0;
  } else {
    return 1;
  }
}

inline double square_center(double x, double y) {
	if (x >= 0.45 && x <= 0.55 && y >= 0.45 && y <= 0.55) {
    return 0;
  } else {
    return 1;
  }
}

inline double square_lower_left(double x, double y) {
	if (x >= 0.05 && x <= 0.15 && y >= 0.25 && y <= 0.35) {
    return 0;
  } else {
    return 1;
  }
}

inline double circle_region(double x, double y) {
	if (0.25*0.25 - ((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5)) > 0) {
		return 0;
	} else {
		return 1;
	}
}

inline double annulus_region(double x, double y) {
	if (0.25*0.25 - ((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5)) > 0) {
		if (0.125*0.125 - ((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5)) < 0) {
			return 0;
		} else {
			return 1;
		}
	} else {
		return 1;
	}
}

#endif
