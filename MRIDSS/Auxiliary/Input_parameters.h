//
//  Input_parameters.h
//  MRIDSS
//
//  Created by Jonathan Squire on 5/7/14.
//  Copyright (c) 2014 Jonathan Squire. All rights reserved.
//

#ifndef __MRIDSS__Input_parameters__
#define __MRIDSS__Input_parameters__

#include "../General_Definitions.h"

// Class for handling the input parameters and storing the data.
// This should be a convenient way to input data and have the methods required to process this in a basic way

class Inputs {
public:
    // Constructor
    Inputs(double nu, double eta, double q,
           double noise,
           int nx, int ny, int nz, double lx,double ly, double lz,
           double t_initial, double t_final, double dt,
           double timvar_save_interval, double fullsol_save_interval,
           bool remapQ, bool QuasiLinearQ) :
    nu(nu), eta(eta), q(q), f_noise(noise),
    NXY{nx,ny}, NZ(nz), L{lx,ly,lz},
    t_initial(t_initial), t_final(t_final), dt(dt),
    timvar_save_interval(timvar_save_interval), fullsol_save_interval(fullsol_save_interval),
    remapQ(remapQ), QuasiLinearQ(QuasiLinearQ)
    { initialize();
    };
    
    ///////////////////////////////////////////////
    /////           USER INPUTS            ////////
    
    // Grid size
    int NXY[2]; // X and Y grid
    int NZ ; // Z grid
    // Box dimensions
    const double L[3];
    
    // Time domain
    const double dt;
    const double t_final; // Final time
    const double t_initial; // Initial time
    // Saving
    const double timvar_save_interval;   // Interval for saving energy, AM etc.
    const double fullsol_save_interval;   // Interval for saving full solution
    
    // Dissipation
    const double nu;   // Viscosity
    const double eta;   // Resisitivity - could ignore for hydrodynamic
    // Shear
    const double q;  // Shear rate
    // Noise
    const double f_noise;  // Driving noise (un-normalized)
    // Shearing box
    bool remapQ;   // Whether to remap or not
    // Include quasi-linear feedback
    bool QuasiLinearQ;
    
    ///////////////////////////////////////////////

    
    ///////////////////////////////////////////////
    /////      Derived quantities          ////////
    int nsteps; // Total number of steps
    int timevar_save_nsteps; // Steps before saving timevar
    int fullsol_save_nsteps; // Steps before full save
    // SB parameters
    double TSB; // Time before remap = LY/LX/q
    int num_before_remap;  // Number of steps before remap
    

    ////////////////////////////////////////////////
    //////              FUNCTIONS             //////
    
    void initialize();

};


#endif /* defined(__MRIDSS__Input_parameters__) */
