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
    void calculatePlans( int MF_size,  int Ckl_size);
    
    /////////////////////////////////////////////////////
    // Functions for exectuting the various FFTs, given pointers to data
    
    // Ckl matrix FFTs
    void for_C2D_d1(dcmplx* Cklin){
        fftw_execute_dft(C2D_for_dim1_,CAST_T0_FFTW(Cklin),CAST_T0_FFTW(Cklin));
    };
    void for_C2D_d2(dcmplx* Cklin){
        fftw_execute_dft(C2D_for_dim2_,CAST_T0_FFTW(Cklin),CAST_T0_FFTW(Cklin));
    };
    void back_C2D_d1(dcmplx* Cklin){
        fftw_execute_dft(C2D_back_dim1_,CAST_T0_FFTW(Cklin),CAST_T0_FFTW(Cklin));
    };
    void back_C2D_d2(dcmplx* Cklin){
        fftw_execute_dft(C2D_back_dim2_,CAST_T0_FFTW(Cklin),CAST_T0_FFTW(Cklin));
    };

    // 1-D MF transforms
    void for_MF1D(dcmplx* MFin,dcmplx* MFout){
        fftw_execute_dft(MF1D_for_,CAST_T0_FFTW(MFin),CAST_T0_FFTW(MFout));
    };
    void back_MF1D(dcmplx* MFin,dcmplx* MFout){
        fftw_execute_dft(MF1D_back_,CAST_T0_FFTW(MFin),CAST_T0_FFTW(MFout));
    };
    void for_IP_MF1D(dcmplx* MFin){
        fftw_execute_dft(MF1D_IP_for_,CAST_T0_FFTW(MFin),CAST_T0_FFTW(MFin));
    };
    void back_IP_MF1D(dcmplx* MFin){
        fftw_execute_dft(MF1D_IP_back_,CAST_T0_FFTW(MFin),CAST_T0_FFTW(MFin));
    };

    // 1-D*NZ MF transforms
    void for_MF2D(dcmplx* mMFin,dcmplx* mMFout){
        fftw_execute_dft(MF2D_for_,CAST_T0_FFTW(mMFin),CAST_T0_FFTW(mMFout));
    };
    void back_MF2D(dcmplx* mMFin,dcmplx* mMFout){
        fftw_execute_dft(MF2D_back_,CAST_T0_FFTW(mMFin),CAST_T0_FFTW(mMFout));
    };
    void for_IP_MF2D(dcmplx* mMFin){
        fftw_execute_dft(MF2D_IP_for_,CAST_T0_FFTW(mMFin),CAST_T0_FFTW(mMFin));
    };
    void back_IP_MF2D(dcmplx* mMFin){
        fftw_execute_dft(MF2D_IP_back_,CAST_T0_FFTW(mMFin),CAST_T0_FFTW(mMFin));
    };
    /////////////////////////////////////////////////////

    
private:
    // Transforms of Ckl for Reynolds stress
    fftw_plan C2D_for_dim1_, C2D_for_dim2_;  // Forward 2-D transform
    fftw_plan C2D_back_dim1_, C2D_back_dim2_;  // Backward 2-D transform
    
    // Transforms of 1-D MF data
    fftw_plan MF1D_for_, MF1D_IP_for_;  // Forward 1-D transform
    fftw_plan MF1D_back_,MF1D_IP_back_; // Backward 1-D transform
    
    // 1-D transforms of matrices of MF-like quantities
    fftw_plan MF2D_for_, MF2D_IP_for_; // Forward transform, don't need dim1 and dim2, since it's NZ 1-D transforms
    fftw_plan MF2D_back_, MF2D_IP_back_; // Backward transform
    
    bool plans_calculated_; // 1 if calculatePlans has been run
    
};



//  Timevar storage, i.e., enregy, angular momentum etc.
class TimeVariables {
public:
    // Constructor
    TimeVariables(int nsteps, int width, int mpinode); // nsteps is number of saves, width is number to save
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
    
private:
    // Basic data
    const int nsteps_;  // Number of saves (length of array)
    const int width_;  // Number of saves for each step in each
    
    // Actual storage arrays
    double ** energy_;
    double ** angular_momentum_;
    double ** dissipation_;
    
    // Current position in the time dimension
    static int curr_pos_;
    
    // MPI node on which the data is stored
    int mpi_node_;
};

#endif /* defined(__MRIDSS__General_Auxiliary__) */
