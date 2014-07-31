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
#include "../Models/Model.h"


// Data is an array of Eigen complex matrices
class Integrator {
public:
    virtual ~Integrator() {};
    virtual int Step(double t, dcmplxVec *MF, dcmplxMat *C) = 0;
    // Reinitialize linear operators (these are stored in integrator)
    virtual void Reinitialize_linear_Ops(double t) = 0;
};


#endif  // INTEGRATOR_H_

