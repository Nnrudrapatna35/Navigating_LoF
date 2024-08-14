/*
* ==============================================================================
*
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
* Authors: REU 2022 Landscape of Fear Group
*   (based on code by Elliot Cartee, Qianli Song, Lexiao Lai)
*     (based on code by Marc Aurèle Gilles and Zachary Clawson)
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
inline double babyGauss(double aX, double aR) {
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

inline double gausskernel(double aX, double aR) {
  if (aX > aR || aX < 0) {
    return 0;
  }
  return 1/(PI*aR*aR) * exp(1)/(exp(1)-1) * exp(-aX*aX/(aR*aR));
}

/****************************************************************************
 ************************ Food Density Functions ****************************
 ****************************************************************************/
inline double gaussian2D(double aX, double aY) {
  return 100 + 750*exp(-((aX-0.8)*(aX-0.8) + (aY-0.8)*(aY-0.8)) / (2*0.09*0.09));
  // return 2.5 + 5*exp(-((aX-0.2)*(aX-0.2) + (aY-0.2)*(aY-0.2)) / (2*0.07*0.07));
}

inline double twoPeaks(double aX, double aY) {
  const double baseline = 0;

  const double muY = 0.8;
  const double muX1 = 0.3;
  const double muX2 = 0.7;
  const double sigma = 0.1;
  const double offset = 0;
  const double scale1 = 4900;
  const double scale2 = 5000;
  const double gaussian1 = varGaussian2D(aX, aY, muX1, muY, sigma, scale1, offset);
  const double gaussian2 = varGaussian2D(aX, aY, muX2, muY, sigma, scale2, offset);

  return baseline + gaussian1 + gaussian2;
}

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

inline double threePeaks(double aX, double aY) {

  return  0.15*exp(-((aX-0.4)*(aX-0.4) + (aY-0.8)*(aY-0.8)) / (2*0.07*0.07))+
    0.55*exp(-((aX-0.2)*(aX-0.2) + (aY-0.4)*(aY-0.4)) / (2*0.07*0.07))+
    0.75*exp(-((aX-0.7)*(aX-0.7) + (aY-0.2)*(aY-0.2)) / (2*0.11*0.11)) + 0.1;
}

/* a VERY PARTICULAR gaussian */
inline double fourPeaks(double aX, double aY) {
  return  5*exp(-((aX-0.2)*(aX-0.2) + (aY-0.2)*(aY-0.2)) / (2*0.07*0.07))+
    7*exp(-((aX-0.2)*(aX-0.2) + (aY-0.8)*(aY-0.8)) / (2*0.11*0.11))+
    9*exp(-((aX-0.8)*(aX-0.8) + (aY-0.2)*(aY-0.2)) / (2*0.14*0.14))+
    3*exp(-((aX-0.8)*(aX-0.8) + (aY-0.8)*(aY-0.8)) / (2*0.09*0.09))+2.5;
}

inline double fourPeaks2(double aX, double aY) {
  // return  50*exp(-((aX-0.2)*(aX-0.2) + (aY-0.2)*(aY-0.2)) / (2*0.07*0.07))+
  //   70*exp(-((aX-0.2)*(aX-0.2) + (aY-0.8)*(aY-0.8)) / (2*0.11*0.11))+
  //   90*exp(-((aX-0.8)*(aX-0.8) + (aY-0.2)*(aY-0.2)) / (2*0.14*0.14))+
  //   30*exp(-((aX-0.8)*(aX-0.8) + (aY-0.8)*(aY-0.8)) / (2*0.09*0.09))+10;

  // return  10*exp(-((aX-0.2)*(aX-0.2) + (aY-0.2)*(aY-0.2)) / (2*0.07*0.07))+
  //   10*exp(-((aX-0.2)*(aX-0.2) + (aY-0.8)*(aY-0.8)) / (2*0.11*0.11))+
  //   10*exp(-((aX-0.8)*(aX-0.8) + (aY-0.2)*(aY-0.2)) / (2*0.14*0.14))+
  //   30*exp(-((aX-0.8)*(aX-0.8) + (aY-0.8)*(aY-0.8)) / (2*0.09*0.09))+1;

  // return  0.15*exp(-((aX-0.2)*(aX-0.2) + (aY-0.2)*(aY-0.2)) / (2*0.07*0.07))+
  //   0.15*exp(-((aX-0.2)*(aX-0.2) + (aY-0.8)*(aY-0.8)) / (2*0.11*0.11))+
  //   0.15*exp(-((aX-0.8)*(aX-0.8) + (aY-0.2)*(aY-0.2)) / (2*0.14*0.14))+
  //   0.75*exp(-((aX-0.8)*(aX-0.8) + (aY-0.8)*(aY-0.8)) / (2*0.09*0.09))+0.01;

  // return  0.4*exp(-((aX-0.2)*(aX-0.2) + (aY-0.2)*(aY-0.2)) / (2*0.07*0.07))+
  //   0.4*exp(-((aX-0.2)*(aX-0.2) + (aY-0.8)*(aY-0.8)) / (2*0.11*0.11))+
  //   0.4*exp(-((aX-0.8)*(aX-0.8) + (aY-0.2)*(aY-0.2)) / (2*0.14*0.14))+
  //   2*exp(-((aX-0.8)*(aX-0.8) + (aY-0.8)*(aY-0.8)) / (2*0.09*0.09))+0.3;

  // return  0.4*exp(-((aX-0.2)*(aX-0.2) + (aY-0.2)*(aY-0.2)) / (2*0.07*0.07))+
  //   0.4*exp(-((aX-0.2)*(aX-0.2) + (aY-0.8)*(aY-0.8)) / (2*0.11*0.11))+
  //   0.4*exp(-((aX-0.8)*(aX-0.8) + (aY-0.2)*(aY-0.2)) / (2*0.14*0.14))+
  //   0.7*exp(-((aX-0.8)*(aX-0.8) + (aY-0.8)*(aY-0.8)) / (2*0.09*0.09))+0.3;

  return  2*exp(-((aX-0.2)*(aX-0.2) + (aY-0.2)*(aY-0.2)) / (2*0.07*0.07))+
    2*exp(-((aX-0.2)*(aX-0.2) + (aY-0.8)*(aY-0.8)) / (2*0.11*0.11))+
    2*exp(-((aX-0.8)*(aX-0.8) + (aY-0.2)*(aY-0.2)) / (2*0.14*0.14))+
    3.9*exp(-((aX-0.8)*(aX-0.8) + (aY-0.8)*(aY-0.8)) / (2*0.09*0.09))+0.1;

}

inline double fourPeaks3(double aX, double aY) {

  // return  10*exp(-((aX-0.2)*(aX-0.2) + (aY-0.2)*(aY-0.2)) / (2*0.07*0.07))+
  //   10*exp(-((aX-0.2)*(aX-0.2) + (aY-0.8)*(aY-0.8)) / (2*0.11*0.11))+
  //   22*exp(-((aX-0.8)*(aX-0.8) + (aY-0.2)*(aY-0.2)) / (2*0.14*0.14))+
  //   17*exp(-((aX-0.8)*(aX-0.8) + (aY-0.8)*(aY-0.8)) / (2*0.09*0.09)) + 8;

  return  5*exp(-((aX-0.2)*(aX-0.2) + (aY-0.2)*(aY-0.2)) / (2*0.07*0.07))+
    5*exp(-((aX-0.2)*(aX-0.2) + (aY-0.8)*(aY-0.8)) / (2*0.11*0.11))+
    17*exp(-((aX-0.8)*(aX-0.8) + (aY-0.2)*(aY-0.2)) / (2*0.14*0.14))+
    10*exp(-((aX-0.8)*(aX-0.8) + (aY-0.8)*(aY-0.8)) / (2*0.09*0.09)) + 8;
}

inline double fourPeaks4(double aX, double aY) {

  // Zero-out food in a region slightly larger than the homebase
  if (aX >= 0.4 && aX <= 0.6 && aY >= 0.4 && aY <= 0.6) {
    return 0;
  }

  double psi_max = 42000000 / 20.0; // Maximum value of food density
  double baseline_food = 0; // Baseline amount of food

  // One peak is low, two are medium, and the top-right is high
  return baseline_food +
    psi_max * 0.5 * exp(-((aX-0.1)*(aX-0.1) + (aY-0.1)*(aY-0.1)) / (2*0.04*0.04)) +
    psi_max * 0.75 * exp(-((aX-0.1)*(aX-0.1) + (aY-0.9)*(aY-0.9)) / (2*0.07*0.07)) +
    psi_max * 0.75 * exp(-((aX-0.9)*(aX-0.9) + (aY-0.1)*(aY-0.1)) / (2*0.07*0.07)) +
    psi_max * exp(-((aX-0.9)*(aX-0.9) + (aY-0.9)*(aY-0.9)) / (2*0.10*0.10));

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

/* 1D Example Food Densities*/
inline double barFood(double aX, double aY) {
  double value;
  if (aY >= 0.9) {
    value = 17500;
  } else {
    value = 0;
  }
  return value;
}

inline double barFoodAcc(double aX, double aY) {
  double value;
  if (aY >= 0.9) {
    value = 17.5;
  } else {
    value = 0;
  }
  return value;
}

inline double breadcrumbs(double aX, double aY) {
  const double baseline = 0; // Make this linear in y, 0.3 at the bottom and 0.4 at the top
  const double sigma = 0.1;
  const double offset = 0;
  const double smallScale = 600;
  const double bigScale = 1000;

  double theta = 160;
  double muX = 0.9 + 0.8*(cos(theta*PI/180));
  double muY = 0.1 + 0.8*(sin(theta*PI/180));
  const double gaussian1 = varGaussian2D(aX, aY, muX, muY, sigma, smallScale, offset);
  theta = 135;
  muX = 0.9 + 0.8*(cos(theta*PI/180));
  muY = 0.1 + 0.8*(sin(theta*PI/180));
  const double gaussian2 = varGaussian2D(aX, aY, muX, muY, sigma, smallScale, offset);
  theta = 110;
  muX = 0.9 + 0.8*(cos(theta*PI/180));
  muY = 0.1 + 0.8*(sin(theta*PI/180));
  const double gaussian3 = varGaussian2D(aX, aY, muX, muY, sigma, smallScale, offset);
  theta = 90;
  muX = 0.9 + 0.8*(cos(theta*PI/180));
  muY = 0.1 + 0.8*(sin(theta*PI/180));
  const double gaussian4 = varGaussian2D(aX, aY, muX, muY, sigma, bigScale, offset);

  return baseline + gaussian1 + gaussian2 + gaussian3 + gaussian4;
}

/* Multistage Example Food Densities*/
inline double foodGrid3(double aX, double aY) {
  const double baseline = 0; // Make this linear in y, 0.3 at the bottom and 0.4 at the top
  const double sigma = 0.05;
  const double offset = 0;

  double muX = 0.2;
  double muY = 0.8;
  double scale = 40;
  const double gaussian1 = varGaussian2D(aX, aY, muX, muY, sigma, scale, offset);

  muX = 0.8;
  scale = 60;
  const double gaussian2 = varGaussian2D(aX, aY, muX, muY, sigma, scale, offset);

  muY = 0.2;
  scale = 50;
  const double gaussian3 = varGaussian2D(aX, aY, muX, muY, sigma, scale, offset);

  return baseline + gaussian1 + gaussian2 + gaussian3;
}

inline double foodGrid6(double aX, double aY) {
  const double baseline = 0; // Make this linear in y, 0.3 at the bottom and 0.4 at the top
  const double sigma = 0.05;
  const double offset = 0;

  double muX = 0.1;
  double muY = 0.4;
  double scale = 20;
  const double gaussian1 = varGaussian2D(aX, aY, muX, muY, sigma, scale, offset);

  muY = 0.7;
  scale = 40;
  const double gaussian2 = varGaussian2D(aX, aY, muX, muY, sigma, scale, offset);

  muX = 0.4;
  muY = 0.1;
  scale = 30;
  const double gaussian3 = varGaussian2D(aX, aY, muX, muY, sigma, scale, offset);

  muX = 0.7;
  scale = 15;
  const double gaussian4 = varGaussian2D(aX, aY, muX, muY, sigma, scale, offset);

  muX = 0.55;
  muY = 0.55;
  scale = 30;
  const double gaussian5 = varGaussian2D(aX, aY, muX, muY, sigma, scale, offset);

  muX = 0.85;
  muY = 0.85;
  scale = 60;
  const double gaussian6 = varGaussian2D(aX, aY, muX, muY, sigma, scale, offset);

  return baseline + gaussian1 + gaussian2 + gaussian3 + gaussian4 + gaussian5 + gaussian6;

}

inline double foodGrid8(double aX, double aY) {
  const double baseline = 0; // Make this linear in y, 0.3 at the bottom and 0.4 at the top
  const double sigma = 0.05;
  const double offset = 0;

  double muX = 0.15;
  double muY = 0.5;
  double scale = 40;
  const double gaussian1 = varGaussian2D(aX, aY, muX, muY, sigma, scale, offset);

  muY = 0.85;
  scale = 40;
  const double gaussian2 = varGaussian2D(aX, aY, muX, muY, sigma, scale, offset);

  muX = 0.5;
  muY = 0.15;
  scale = 40;
  const double gaussian3 = varGaussian2D(aX, aY, muX, muY, sigma, scale, offset);

  muY = 0.5;
  scale = 50;
  const double gaussian4 = varGaussian2D(aX, aY, muX, muY, sigma, scale, offset);

  muY = 0.85;
  scale = 50;
  const double gaussian5 = varGaussian2D(aX, aY, muX, muY, sigma, scale, offset);

  muX = 0.85;
  muY = 0.15;
  scale = 50;
  const double gaussian6 = varGaussian2D(aX, aY, muX, muY, sigma, scale, offset);

  muY = 0.5;
  scale = 40;
  const double gaussian7 = varGaussian2D(aX, aY, muX, muY, sigma, scale, offset);

  muY = 0.85;
  scale = 60;
  const double gaussian8 = varGaussian2D(aX, aY, muX, muY, sigma, scale, offset);

  return baseline + gaussian1 + gaussian2 + gaussian3 + gaussian4 + gaussian5
                  + gaussian6 + gaussian7 + gaussian8;

}

inline double randomPeaks(double aX, double aY) {
  const int numPeaks = 10;
  double muX[numPeaks] = {0.45,0.68,0.65,0.74,0.78,0.63,0.81,0.90,0.52,0.50};
  double muY[numPeaks] = {0.70,0.86,0.34,0.59,0.22,0.77,0.91,0.63,0.90,0.50};
  double sigma[numPeaks] = {0.75,0.26,0.51,0.70,0.89,0.34,0.78,0.55,0.69,10.0};
  double scale[numPeaks] = {0.59,0.55,0.47,0.75,0.86,0.66,0.99,0.81,0.69,0.10};

  double value = 0;
  const double offset = 0;

  for (int i = 0; i < numPeaks; ++i) {
    value += varGaussian2D(aX, aY, muX[i], muY[i], (1.0/10.0)*sigma[i], 100.0*scale[i], offset);
  }

  return value;
}

inline double gerbilPeaks(double aX, double aY) {
  const int numPeaks = 10;
  double muX[numPeaks] = {0.8444, 0.3045, 0.7805, 0.6753, 0.1067, 0.3022, 0.3868, 0.9160, 0.1012, 0.4624};
  double muY[numPeaks] = {0.9577, 0.3007, 0.6761, 0.2891, 0.5718, 0.6951, 0.0680, 0.2548, 0.2240, 0.9678};
  double sigma = 0.05;
  double scale = 400;

  double value = 0;
  const double offset = 0;

  for (int i = 0; i < numPeaks; ++i) {
    value += varGaussian2D(aX, aY, muX[i], muY[i], sigma, scale, offset);
  }

  return value;
}

/****************************************************************************
 ************************ Spotting Rate Functions ***************************
 ****************************************************************************/
/* Writing a new function to be used as BIG predator density*/
inline double fourPeaksPred(double aX, double aY, double aE, double aT) {
  //return 100*(aX+aY+0.2);
  return 100*exp(-((aX-0.25)*(aX-0.25) + (aY-0.25)*(aY-0.25)) / (2*0.02*0.02))+
  140*exp(-((aX-0.25)*(aX-0.25) + (aY-0.75)*(aY-0.75)) / (2*0.025*0.025))+
  180*exp(-((aX-0.75)*(aX-0.75) + (aY-0.25)*(aY-0.25)) / (2*0.03*0.03))+
  150*exp(-((aX-0.75)*(aX-0.75) + (aY-0.75)*(aY-0.75)) / (2*0.035*0.035))+25;
}

inline double fourPeaksPred2(double aX, double aY, double aE, double aT) {
  return 30*exp(-((aX-0.25)*(aX-0.25) + (aY-0.25)*(aY-0.25)) / (2*0.02*0.02))+
  40*exp(-((aX-0.25)*(aX-0.25) + (aY-0.75)*(aY-0.75)) / (2*0.025*0.025))+
  60*exp(-((aX-0.75)*(aX-0.75) + (aY-0.25)*(aY-0.25)) / (2*0.03*0.03))+
  50*exp(-((aX-0.75)*(aX-0.75) + (aY-0.75)*(aY-0.75)) / (2*0.035*0.035))+10;
}

inline double fourPeaksPred3(double aX, double aY, double aE, double aT) {

  double baseline_rate = 0; // 0.05;

  return baseline_rate +
  1.0 * exp(-((aX-0.25)*(aX-0.25) + (aY-0.25)*(aY-0.25)) / (2*0.02*0.02)) +
  0.2 * exp(-((aX-0.25)*(aX-0.25) + (aY-0.75)*(aY-0.75)) / (2*0.035*0.035)) +
  0.2 * exp(-((aX-0.75)*(aX-0.75) + (aY-0.25)*(aY-0.25)) / (2*0.035*0.035)) +
  0.5 * exp(-((aX-0.75)*(aX-0.75) + (aY-0.75)*(aY-0.75)) / (2*0.05*0.05));

}

inline double fourPeaksPred4(double aX, double aY, double aE, double aT) {

  double baseline_rate = 0;

  if (aX >= 0.2 && aX <= 0.3 && aY >= 0.2 && aY <= 0.3) {
    return baseline_rate + 1.0;
  }

  if (aX >= 0.2 && aX <= 0.3 && aY >= 0.7 && aY <= 0.8) {
    return baseline_rate + 0.2;
  }

  if (aX >= 0.7 && aX <= 0.8 && aY >= 0.2 && aY <= 0.3) {
    return baseline_rate + 0.2;
  }

  if (aX >= 0.7 && aX <= 0.8 && aY >= 0.7 && aY <= 0.8) {
    return baseline_rate + 0.5;
  }

  return baseline_rate;

}

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

inline double ellipsePredCts(double aX, double aY, double aE, double aT) {
  const double baseline = 0.3; // Make this linear in y, 0.3 at the bottom and 0.4 at the top

  const double muX = 0.5;
  const double muY1 = 0.3;
  const double muY2 = 0.7;
  const double sigmaX = 0.2;
  const double sigmaY = 0.05;
  const double offset = 0;
  const double scale1 = 0.5;
  const double scale2 = 1.0;
  const double gaussian1 = varEllipseGaussian2D(aX, aY, muX, muY1, sigmaX, sigmaY, scale1, offset);
  const double gaussian2 = varEllipseGaussian2D(aX, aY, muX, muY2, sigmaX, sigmaY, scale2, offset);

  return baseline + gaussian1 + gaussian2;
}

inline double gerbilPredator(double aX, double aY, double aE, double aT) {
  const int numPeaks = 2;
  double muX[numPeaks] = {0.95, 0.05};
  double muY[numPeaks] = {0.95, 0.05};
  double sigma = 0.3;
  double scale = 0.5;

  double value = 0;
  const double offset = 0.5;

  for (int i = 0; i < numPeaks; ++i) {
    value += varGaussian2D(aX, aY, muX[i], muY[i], sigma, scale, offset);
  }

  return value;
}

inline double gerbilGiveUp(double aX, double aY, double aE, double aT) {
  const int numPeaks = 2;
  double muX[numPeaks] = {0.95, 0.05};
  double muY[numPeaks] = {0.95, 0.05};
  double sigma = 0.3;
  double scale = 0.5;

  double value = 0;
  const double offset = -1.0;

  for (int i = 0; i < numPeaks; ++i) {
    value += varGaussian2D(aX, aY, muX[i], muY[i], sigma, scale, offset);
  }

  return -1*value;
}

/****************************************************************************
 ******************** Give Up and Kill Rate Functions ***********************
 ****************************************************************************/
/* muG is inverse of predator density*/
inline double invFourPeaksPred(double aX, double aY, double aE, double aT) {
  return 2500 / (100*exp(-((aX-0.25)*(aX-0.25) + (aY-0.25)*(aY-0.25)) / (2*0.02*0.02))+
  140*exp(-((aX-0.25)*(aX-0.25) + (aY-0.75)*(aY-0.75)) / (2*0.025*0.025))+
  180*exp(-((aX-0.75)*(aX-0.75) + (aY-0.25)*(aY-0.25)) / (2*0.03*0.03))+
  150*exp(-((aX-0.75)*(aX-0.75) + (aY-0.75)*(aY-0.75)) / (2*0.035*0.035))+25);
}

inline double invFourPeaksPred2(double aX, double aY, double aE, double aT) {
  return 1000 / (100*exp(-((aX-0.25)*(aX-0.25) + (aY-0.25)*(aY-0.25)) / (2*0.02*0.02))+
  140*exp(-((aX-0.25)*(aX-0.25) + (aY-0.75)*(aY-0.75)) / (2*0.025*0.025))+
  180*exp(-((aX-0.75)*(aX-0.75) + (aY-0.25)*(aY-0.25)) / (2*0.03*0.03))+
  150*exp(-((aX-0.75)*(aX-0.75) + (aY-0.75)*(aY-0.75)) / (2*0.035*0.035))+25); // testing with numerator
}

inline double invFourPeaksPred3(double aX, double aY, double aE, double aT) {
  return 500 / (100*exp(-((aX-0.25)*(aX-0.25) + (aY-0.25)*(aY-0.25)) / (2*0.02*0.02))+
  140*exp(-((aX-0.25)*(aX-0.25) + (aY-0.75)*(aY-0.75)) / (2*0.025*0.025))+
  180*exp(-((aX-0.75)*(aX-0.75) + (aY-0.25)*(aY-0.25)) / (2*0.03*0.03))+
  150*exp(-((aX-0.75)*(aX-0.75) + (aY-0.75)*(aY-0.75)) / (2*0.035*0.035))+25);
}

inline double invFourPeaksPred4(double aX, double aY, double aE, double aT) {
  return 167 / (30*exp(-((aX-0.25)*(aX-0.25) + (aY-0.25)*(aY-0.25)) / (2*0.02*0.02))+
  40*exp(-((aX-0.25)*(aX-0.25) + (aY-0.75)*(aY-0.75)) / (2*0.025*0.025))+
  60*exp(-((aX-0.75)*(aX-0.75) + (aY-0.25)*(aY-0.25)) / (2*0.03*0.03))+
  50*exp(-((aX-0.75)*(aX-0.75) + (aY-0.75)*(aY-0.75)) / (2*0.035*0.035))+10);
}

inline double twoEllipseGaussians(double aX, double aY, double aE, double aT) {
  const double baseline = 0.5; // Make this linear in y, 0.3 at the bottom and 0.4 at the top

  const double muX1 = 0.4;
  const double muX2 = 0.55;
  const double muY1 = 0.8;
  const double muY2 = 0.2;
  const double sigmaX1 = 0.3;
  const double sigmaX2 = 0.2;
  const double sigmaY = 0.1;
  const double offset = 0;
  const double scale = 2.0;
  const double gaussian1 = varEllipseGaussian2D(aX, aY, muX1, muY1, sigmaX1, 
                                                sigmaY, scale, offset);
  const double gaussian2 = varEllipseGaussian2D(aX, aY, muX2, muY2, sigmaX2, 
                                                sigmaY, scale, offset);

  return baseline + gaussian1 + gaussian2;
}

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

inline double square_region_real(double x, double y) {
	if (x >= 0.86 && x <= 0.88 && y >= 0.86 && y <= 0.88) {
    return 0;
  } else {
    return 1;
  }
}

// this is TWO-PART homebase for REAL example7
inline double square_region_real2(double x, double y) {
	if (x >= 0.55 && x <= 0.57 && y >= 0.55 && y <= 0.57) {
    return 0;
  }
  else if (x >= 0.34 && x <= 0.36 && y >= 0.72 && y <= 0.74) {
    return 0;
  }
  else {
    return 1;
  }
}

inline double square_region_real3(double x, double y) {
	if (x >= 0.65 && x <= 0.67 && y >= 0.55 && y <= 0.57) {
    return 0;
  }
  else if (x >= 0.83 && x <= 0.85 && y >= 0.49 && y <= 0.51) {
    return 0;
  }
  else {
    return 1;
  }
}

inline double square_region_double(double x, double y) {
	if (x >= 0.55 && x <= 0.65 && y >= 0.55 && y <= 0.65) {
    return 0;
  }
  // the second homebase!
  else if (x >= 0.9 && x <= 1 && y >= 0.9 && y <= 1) {
    return 0;
  }
  else {
    return 1;
  }
}

inline double square_region2(double x, double y) {
	if (x >= 0.75 && x <= 0.85 && y >= 0.25 && y <= 0.35) {
    return 0;
  } else {
    return 1;
  }
}

inline double square_region3(double x, double y) {
	if (x >= 0.6 && x <= 1 && y >= 0 && y <= 0.5) {
    return 0;
  } else {
    return 1;
  }
}

inline double square_region5(double x, double y) {
	if (x >= 0.0 && x <= 0.2 && y >= 0.0 && y <= 0.2) {
    return 0;
  } else {
    return 1;
  }
}

inline double square_region6(double x, double y) {
	if (x >= 0.45 && x <= 0.55 && y >= 0.15 && y <= 0.25) {
    return 0;
  } else {
    return 1;
  }
}

inline double rectangle_region(double x, double y) {
  if (x >= 0 && x <= .25 && y >= 0 && y <= 0.5) {
    return 0;
  } else {
    return 1;
  }
}

inline double rectangle_region2(double x, double y) {
  // h = 1 / 32 for coarsest grid (33x33x5xnt)
  // 8 / 32 = 0.25
  // 16 / 32 = 0.5
  // 12 / 32 = 0.375
  // 22 / 32 = 0.6875
  if (x >= 0 && x <= .375 && y >= 0 && y <= 0.6875) {
    return 0;
  } else {
    return 1;
  }
}

inline double rectangle_region4(double x, double y) {
  // h = 1 / 32 for coarsest grid (33x33x5xnt)
  // 8 / 32 = 0.25
  // 16 / 32 = 0.5
  // 12 / 32 = 0.375
  // 22 / 32 = 0.6875
  if (x >= 0 && x <= .1 && y >= 0 && y <= 0.1) {
    return 0;
  } else {
    return 1;
  }
}

inline double rectangle_region5(double x, double y) {
  // h = 1 / 32 for coarsest grid (33x33x5xnt)
  // 8 / 32 = 0.25
  // 16 / 32 = 0.5
  // 12 / 32 = 0.375
  // 22 / 32 = 0.6875
  if (x >= 0.21 && x <= .35 && y >= 0.21 && y <= 0.35) {
    return 0;
  } else {
    return 1;
  }
}

inline double rectangle_region3(double x, double y) {
  // h = 1 / 32 for coarsest grid (33x33x5xnt)
  // 8 / 32 = 0.25
  // 16 / 32 = 0.5
  // 12 / 32 = 0.375
  // 22 / 32 = 0.6875
  if (x >= 0 && x <= 0.25 && y >= 0 && y <= 0.5) {
    return 0;
  } else {
    return 1;
  }
}

/* A continuous version of rectange_region3. Varies continuously from 0 to aC.
  Changed the y bound and added arguments. */
inline double continuous_rectangle_region3(double aX, double aY, double aE, double aT, double aC) {
  // h = 1 / 32 = 0.03125 for coarsest grid (33x33x5xnt)
  // 2 / 32 = 0.0625
  // 4 / 32 = 0.125
  // 8 / 32 = 0.25
  // 16 / 32 = 0.5
  // 12 / 32 = 0.375
  // 22 / 32 = 0.6875
  if (aX >= 0 && aX <= 0.25 && aY >= 0 && aY <= 0.5) {
    return 0;
  }
  else if (aX > 0.25 && aX <= 0.3125 && aY > 0.5 && aY <= 0.625) {
    double epsilon = 0.0625 / 0.25;
    double factor = std::max( aX/0.25, aY/0.5 );
    return (factor-1) * aC / epsilon;
  }
  else {
    return aC;
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
