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

#include "../Auxiliary/MPIdata.h"

// Class for handling the input parameters and storing the data.
// This should be a convenient way to input data and have the methods required to process this in a basic way

// Really more like a struct

class Inputs {
public:
    // Constructor
    Inputs(const MPIdata& mpi,const std::string& input_file_name);
    

    
    ///////////////////////////////////////////////
    /////       USER INPUTS (through file)   //////
    
    // Grid size
    int NXY[2]; // X and Y grid
    int NZ ; // Z grid
    // Box dimensions
    double L[3];
    
    // Time domain
    double dt; // Negative for variable time-step 
    double CFL; // CFL number, might be ignored (for a fixed step integrator)
    double t_final; // Final time
    double t_initial; // Initial time
    double t_start; // Start time (different when loading from file)
    int i_start; // Start step (different from 0 when loading from file)

    // Saving
    double timvar_save_interval;   // Interval for saving energy, AM etc.
    double fullsol_save_interval;   // Interval for saving full solution
    
    // Mean magnetic fields
    double B0z;
    // Dissipation
    double nu;   // Viscosity
    double eta;   // Resisitivity - could ignore for hydrodynamic
    // Shear
    double q;  // Shear rate
    double omega; // Rotation
    // Noise
    double f_noise;  // Driving noise (un-normalized)
    double noise_range_low; // Low driving noise cutoff (default 0)
    double noise_range_high; // High driving noise cutoff (default inf)
    double v_noise_mult, b_noise_mult; // Multiply individual noise values by these amounts
    
    // Shearing box
    bool remapQ;   // Whether to remap or not
    // Include quasi-linear feedback
    bool QuasiLinearQ;
    
    // Time variables to save
    bool energy_save_Q;
    bool AM_save_Q;
    bool dissipation_save_Q;
    bool reynolds_save_Q;
    bool mean_field_save_Q;
    
    // Initial conditions
    double initial_By;
    
    // Start using saved data
    bool start_from_saved_Q;
    
    // Check that equations are correct in main.cc
    std::string equations_to_use;
    
    
    ///////////////////////////////////////////////

    
    ///////////////////////////////////////////////
    /////      Derived quantities          ////////
    int nsteps; // Total number of steps
    int timevar_save_nsteps; // Steps before saving timevar
    int fullsol_save_nsteps; // Steps before full save
    // SB parameters
    double TSB; // Time before remap = LY/LX/q
    int num_before_remap;  // Number of steps before remap
    // Times to use in saving
    double fullsave_t;
    double timevar_t;
    
    // Directory for data storage
    std::string simulation_dir;
    
    
    
    
    

    ////////////////////////////////////////////////
    //////              FUNCTIONS             //////
    
    void initialize_();
    // Read from text file
    template <class T>
    T Read_From_Input_File_(const std::string& varstr, const std::string& full_file, T default_value);
    // Find input file .DSSinput
    std::string Find_Input_File(const std::string& directory,const std::string& input_name);
    
private:
    int mpi_node_; // Mpi ode for internal usage

};


#endif /* defined(__MRIDSS__Input_parameters__) */
