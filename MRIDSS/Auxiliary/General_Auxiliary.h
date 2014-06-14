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

void ShearingBox_Remap(Model* mod, dcmplxMat* Ckl);

// FFTW plan storage (and execution) class
class fftwPlans {
    // Stores 1-D and 2-D transforms for the mean fields and Ckl
public:
    // Constructor
    fftwPlans(): plans_calculated_(0) {};
    // Destructor
    ~fftwPlans();
    // No Copy constructor, but should only have one instance anyway
    // Setup plans
    void calculatePlans( int NZ );
    
    /////////////////////////////////////////////////////
    // Functions for exectuting the various FFTs, given pointers to data
    
    // Ckl matrix FFTs
    void for_2D_dim1(dcmplx* Cklin){
        fftw_execute_dft(t2D_for_dim1_,CAST_T0_FFTW(Cklin),CAST_T0_FFTW(Cklin));
    };
    void for_2D_dim2(dcmplx* Cklin){
        fftw_execute_dft(t2D_for_dim2_,CAST_T0_FFTW(Cklin),CAST_T0_FFTW(Cklin));
    };
    void back_2D_dim1(dcmplx* Cklin){
        fftw_execute_dft(t2D_back_dim1_,CAST_T0_FFTW(Cklin),CAST_T0_FFTW(Cklin));
    };
    void back_2D_dim2(dcmplx* Cklin){
        fftw_execute_dft(t2D_back_dim2_,CAST_T0_FFTW(Cklin),CAST_T0_FFTW(Cklin));
    };

    // 1-D MF transforms
    void for_1D(dcmplx* MFin){
        fftw_execute_dft(t1D_for_,CAST_T0_FFTW(MFin),CAST_T0_FFTW(MFin));
    };
    void back_1D(dcmplx* MFin){
        fftw_execute_dft(t1D_back_,CAST_T0_FFTW(MFin),CAST_T0_FFTW(MFin));
    };

    /////////////////////////////////////////////////////

    
private:
    // Transforms of 2-D data
    // Need both dimensions (sceond is used for 2-D transform in Reynolds stress)
    fftw_plan t2D_for_dim1_, t2D_for_dim2_;  // Forward 2-D transform
    fftw_plan t2D_back_dim1_, t2D_back_dim2_;  // Backward 2-D transform
    
    // Transforms of 1-D MF data
    fftw_plan t1D_for_, t1D_back_;  // Forward 1-D transform
    
    bool plans_calculated_; // 1 if calculatePlans has been run
    
};



//  Timevar storage, i.e., enregy, angular momentum etc.
class TimeVariables {
public:
    // Constructor
    TimeVariables(Inputs SP,  int width, int mpinode); // nsteps is number of saves, width is number to save
    // Destructor
    ~TimeVariables();
    
    // Increase time save number (could improve later)
    void increase_savenum() {curr_pos_++;};
    // Returns pointer to current array
    double* current_energy() {return energy_[curr_pos_];};
    double* energy_at_n(int n) {return energy_[n];};
    double* current_AM() {return angular_momentum_[curr_pos_];};
    double* current_diss() {return dissipation_[curr_pos_];};
    
    int number_of_saves() {return nsteps_;};
    
    // MPI node
    int tv_mpi_node() {return mpi_node_;};
    
    // Save data - Change to each time-step?
    void Save_Data();
    
    // Save mean fields - could improve
    void Save_Mean_Fields(dcmplxVec *MFdata, int numMF, fftwPlans& fft);
    
    
private:
    // Basic data
    const int nsteps_;  // Number of saves (length of array)
    const int width_;  // Number of saves for each step in each
    
    // Actual storage arrays
    double ** energy_;
    double ** angular_momentum_;
    double ** dissipation_;
    // Saving these variables - file streams
    std::ofstream fileS_energy_;
    std::ofstream fileS_angular_momentum_;
    std::ofstream fileS_dissipation_;
    // Saving these variables - file name
    std::string fname_energy_;
    std::string fname_angular_momentum_;
    std::string fname_dissipation_;
    
    // Mean fields
    std::ofstream fileS_mean_fields_;
    std::string fname_mean_fields_;
    dcmplxVec MFdata_c_; // Space for taking fft and converting to real
    doubVec MFdata_d_;
    long MFlen_; // Set to zero once and change at first save
    
    // Directory
    std::string simulation_dir_;
    
    // Current position in the time dimension
    static int curr_pos_;

    
    // MPI node on which the data is stored
    int mpi_node_;
};

#endif /* defined(__MRIDSS__General_Auxiliary__) */
