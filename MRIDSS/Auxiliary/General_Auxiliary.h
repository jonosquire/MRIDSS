//
//  General_Auxiliary.h
//  MRIDSS
//
//  Created by Jonathan Squire on 5/5/14.
//  Copyright (c) 2014 Jonathan Squire. All rights reserved.
//

#ifndef __MRIDSS__General_Auxiliary__
#define __MRIDSS__General_Auxiliary__

#include "../General_Definitions.h"

#include "../Models/Model.h"
#include "../Integrators/Integrator.h"
#include "Input_parameters.h"
#include "fftwPlans.h"

void ShearingBox_Remap(Model* mod, dcmplxMat* Ckl);

// Check aspects of the solution
void CheckSolution(dcmplxVec* MF, dcmplxMat* Ckl, fftwPlans& fft);


// Save full data set.
// Currently this is just used so a simulation can be continued upon finish
void Save_Full_Data_for_Restart(Model* fluidEqs, Inputs& SP, MPIdata& mpi, double t, dcmplxVec* MF, dcmplxMat* Ckl);
// Load Data set
void Load_Full_Data_for_Restart(Model* fluidEqs, Inputs& SP, MPIdata& mpi, dcmplxVec* MF, dcmplxMat* Ckl);


#endif /* defined(__MRIDSS__General_Auxiliary__) */
