//
//  TimeVariables.cc
//  MRIDSS
//
//  Created by Jonathan Squire on 6/25/14.
//  Copyright (c) 2014 Jonathan Squire. All rights reserved.
//

#include "TimeVariables.h"






/////////////////////////////////////////////////////
//                                                 //
//                 TimeVariables                   //
//                                                 //
/////////////////////////////////////////////////////

// Current position
int TimeVariables::curr_pos_ = 0;

// Constructor
TimeVariables::TimeVariables(Inputs& SP, int ENwidth, int numMF, int num_reynolds, int mpinode) :
width_(ENwidth), num_reynolds_saves_(num_reynolds),
mpi_node_(mpinode),
simulation_dir_(SP.simulation_dir), MFlen_(0), numMF_(numMF),
en_save_Q_(SP.energy_save_Q), AM_save_Q_(SP.AM_save_Q), diss_save_Q_(SP.dissipation_save_Q),
rey_save_Q_(SP.reynolds_save_Q),
MF_save_Q_(SP.mean_field_save_Q),
clk_diff_(0.0)// Timing
{
    ////////////////////////////////
    //    THIS NEEDS TIDYING!!!!!!!!
    
    mpi_node_bool_ =  (mpi_node_ ==0);
    // Main variables
    // Store only on 1-processor - presumbably this is better
    // TODO: Only allocate necessary memory?
    if (mpi_node_bool_) {
        energy_ = new double[width_];
        angular_momentum_ = new double[width_];
        dissipation_ = new double[width_];
        reynolds_ = new double[num_reynolds_saves_];
    }
    
    // Saving to disk - TODO improve to save at regular intervals?
    fname_energy_ = simulation_dir_ + "energy.dat";
    fname_angular_momentum_ = simulation_dir_ + "angular_momentum.dat";
    fname_dissipation_ = simulation_dir_ + "dissipation.dat";
    // Reynolds stress
    fname_reynolds_ = simulation_dir_ + "reynolds_stress.dat";
    // Mean fields
    fname_mean_fields_ = simulation_dir_ + "mean_fields.dat";
    // Time
    fname_time_ = simulation_dir_ + "time_vec.dat";
    
    std::ios_base::openmode o_mode;
    open_append = std::ios::binary | std::ios::app; // Keep since used regularly
    if (SP.start_from_saved_Q) {
        o_mode = open_append;
    } else {
        o_mode = std::ios::binary | std::ios::trunc;
    }
    
    // REMOVE APPEND - OTHERWISE ITS ANNOYING!
    if (mpi_node_bool_) {
        
        
        if (en_save_Q_) {
            fileS_energy_.open(fname_energy_.c_str(), o_mode);
            if (!fileS_energy_.is_open() ) {
                std::cout << "WARNING: " << fname_energy_ <<  " file unable to be opened" << std::endl;
            }
            fileS_energy_.close();
        }
        if (AM_save_Q_) {
            fileS_angular_momentum_.open(fname_angular_momentum_.c_str(), o_mode);
            if (!fileS_angular_momentum_.is_open() ) {
                std::cout << "WARNING: " << fname_angular_momentum_ <<  " file unable to be opened" << std::endl;
            }
            fileS_angular_momentum_.close();
        }
        if (diss_save_Q_) {
            fileS_dissipation_.open(fname_dissipation_.c_str(), o_mode);
            if (!fileS_dissipation_.is_open() ) {
                std::cout << "WARNING: " << fname_dissipation_ <<  " file unable to be opened" << std::endl;
            }
            fileS_dissipation_.close();
        }
        // Reynolds stress
        if (rey_save_Q_) {
            fileS_reynolds_.open(fname_reynolds_.c_str(), o_mode);
            if (!fileS_reynolds_.is_open() ) {
                std::cout << "WARNING: " << fname_reynolds_ <<  " file unable to be opened" << std::endl;
            }
            fileS_reynolds_.close();
        }
        
        // Mean fields
        if (MF_save_Q_) {
            fileS_mean_fields_.open(fname_mean_fields_.c_str(), o_mode);
            if (!fileS_mean_fields_.is_open() ) {
                std::cout << "WARNING: " << fname_mean_fields_ <<  " file unable to be opened" << std::endl;
            }
            fileS_mean_fields_.close();
        }
        //        Time data
        fileS_time_.open(fname_time_.c_str(), o_mode);
        if (!fileS_time_.is_open() ) {
            std::cout << "WARNING: " << fname_time_ <<  " file unable to be opened" << std::endl;
        }
        fileS_time_.close();

    }
    
    
}

// Destructor
TimeVariables::~TimeVariables() {
    if (mpi_node_bool_) {

        delete[] energy_;
        delete[] angular_momentum_;
        delete[] dissipation_;
        delete[] reynolds_;
    }
    
    
//    fileS_energy_.close();
//    fileS_angular_momentum_.close();
//    fileS_dissipation_.close();
//    fileS_reynolds_.close();
//    fileS_mean_fields_.close();
//    fileS_time_.close();

    
}


// Save the data to disk
void TimeVariables::Save_Data(double t){
    // Run only on root process
    if (mpi_node_bool_){
        
        // Open all the files
        fileS_energy_.open(fname_energy_.c_str(), open_append);
        fileS_angular_momentum_.open(fname_angular_momentum_.c_str(), open_append);
        fileS_dissipation_.open(fname_dissipation_.c_str(), open_append);
        fileS_reynolds_.open(fname_reynolds_.c_str(), open_append);
        fileS_time_.open(fname_time_.c_str(), open_append);
        
        
        if (en_save_Q_ ) {
            // Must write each individually since don't know that it is continuous in memory
            fileS_energy_.write( (char*) energy_, sizeof(energy_)*width_);
            // Read into matlab as a double (4*nsteps) array
        }
        if (AM_save_Q_ ) {
            // Must write each individually since don't know that it is continuous in memory
            fileS_angular_momentum_.write( (char*) angular_momentum_, sizeof(energy_)*width_);
            // Read into matlab as a double (4*nsteps) array
        }
        if (diss_save_Q_) {
            // Must write each individually since don't know that it is continuous in memory
               fileS_dissipation_.write( (char*) dissipation_, sizeof(energy_)*width_);
            // Read into matlab as a double (4*nsteps) array
        }
        if (rey_save_Q_ ) {
            // Must write each individually since don't know that it is continuous in memory
                fileS_reynolds_.write( (char*) reynolds_, sizeof(reynolds_)*num_reynolds_saves_);

            // Read into matlab as a double (num_reynolds_saves_*nsteps) array
        }
        // Write
        fileS_time_.write( (char*) &t, sizeof(double) );
        
        // Close all the files
        fileS_energy_.close();
        fileS_angular_momentum_.close();
        fileS_dissipation_.close();
        fileS_reynolds_.close();
        fileS_time_.close();
        
    }
    
    
    
}


// Save Mean field data
void TimeVariables::Save_Mean_Fields(dcmplxVec *MF, fftwPlans& fft){
    // Could certainly be improved - just dumps the data into fname_mean_fields
    // Start timing
    clk_start_ = clock();
    
    // IN MATLAB
    // Saves MF[i] sequentially from 0 to numMF. Data is double
    
    if ( MF_save_Q_ && mpi_node_bool_ ){
        
        fileS_mean_fields_.open(fname_mean_fields_.c_str(),open_append);
        // Initialize data
        if (MFlen_ == 0 ) {
            MFlen_ = MF[0].size();
            MFdata_c_ = dcmplxVec( MFlen_);
            MFdata_d_ = doubVec (MFlen_);
        }
        
        // Save each MF variable sequentially
        for (int i=0; i<numMF_; ++i) {
            // Take fft - could be done faster with c_to_r fft
            MFdata_c_ = MF[i];
            fft.back_1D(MFdata_c_.data());
            MFdata_d_ = MFdata_c_.real()/MFlen_;
            //Save
            fileS_mean_fields_.write( (char*) MFdata_d_.data(), sizeof(double)*MFlen_);
        }
        
        fileS_mean_fields_.close();
    }
    // Timing
    clk_diff_ += (clock() - clk_start_ )/ (double)CLOCKS_PER_SEC;
}
//                                                 //
/////////////////////////////////////////////////////
