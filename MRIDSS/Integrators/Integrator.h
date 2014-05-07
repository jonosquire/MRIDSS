//
//  Integrator.h
//  TwoDFluid
//
//  Created by Jonathan Squire on 4/24/14.
//  Copyright (c) 2014 J Squire. All rights reserved.
//

#ifndef INTEGRATOR_H_
#define INTEGRATOR_H_

// Abstract integrator class for Direct Statistical Simulation

#include "../General_Definitions.h"

typedef std::complex<double> dcmplx;
typedef Eigen::Matrix<dcmplx,Eigen::Dynamic,Eigen::Dynamic> dcmplxMat;
// USE COLUMN_MAJOR FORMAT FOR MATRIX. THE ONLY PLACE WHERE THIS MIGHT CAUSE A PROBLEM IS IN FFTW, SO LONG AS I'M CAREFUL SHOULD BE EASY
// Reason for this choice is to facilitate transformation over from Matlab
typedef Eigen::Array<dcmplx, Eigen::Dynamic, 1> dcmplxVec;// Use array since no matrix multiplications
typedef Eigen::Array<double, Eigen::Dynamic, 1> doubVec; // For use in linear part

// Data is an array of Eigen complex matrices
class Integrator {
public:
    virtual ~Integrator() {};
    virtual int Step(double t, dcmplxVec *MF, dcmplxMat *C) = 0;
};


#endif  // INTEGRATOR_H_

