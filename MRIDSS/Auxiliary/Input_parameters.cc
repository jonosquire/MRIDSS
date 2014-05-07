//
//  Input_parameters.cc
//  MRIDSS
//
//  Created by Jonathan Squire on 5/7/14.
//  Copyright (c) 2014 Jonathan Squire. All rights reserved.
//

#include "Input_parameters.h"

// Class for handling the input parameters and storing the data.
// This should be a convenient way to input data and have the methods required to process this in a basic way


void Inputs::initialize() {
    // Time variables
    nsteps = t_final/dt;
    timevar_save_nsteps = timvar_save_interval/dt;
    fullsol_save_nsteps = fullsol_save_interval/dt;
    
    // Shearing box
    num_before_remap = 2*nsteps; // If no remap
    if (remapQ) {
        TSB = L[1]/L[0]/q;
        // Test for an integer number of steps
        if (TSB/dt-floor(TSB/dt) != 0.0) {
            std::cout << "Warning, non-integer number of steps before remapping: " << TSB/dt <<std::endl;
        }
        num_before_remap = TSB/dt;
    }
}