//
//  General_Auxiliary.cc
//  MRIDSS
//
//  Created by Jonathan Squire on 5/5/14.
//  Copyright (c) 2014 Jonathan Squire. All rights reserved.
//

#include "General_Auxiliary.h"


// Remapping for the shearing box
// This is done in the opposite way to my matlab code -- rather than moving C, just change kx. This then requires nothing complicated and no MPI communication
void ShearingBox_Remap(Model* mod, dcmplxMat* Ckl){
// FRIEND TO THE MODEL CLASS - THIS CHANGES THE kx_ array;
    dcmplx kxLfac(0,-1/(2*PI/mod->L_[0]) );
    dcmplx kyLfac(0,-1/(2*PI/mod->L_[1]) );
    
    dcmplx * kxp = mod->kx_;// Pointer to kxp data
    for (int i=0; i<mod->Cdimxy(); ++i) {
        // Find new kx vector
        int k_i = i + mod->index_for_k_array(); // Index in kx_, ky_ for each process
        kxp[k_i] = kxp[k_i]*kxLfac + mod->ky_[k_i]*kyLfac;
        
        // If kx_ is too large zero out Ckl and set kx negative
        double overlap = kxp[k_i].real()+1 - mod->Nxy_[0]/2;
        if (overlap>0) {
            Ckl[i].setZero();
            kxp[k_i] = -mod->Nxy_[0]/2+overlap;
        }
        // Un-normalize again
        kxp[k_i] = kxp[k_i]/kxLfac;
    }
}



/////////////////////////////////////////////////////
//                                                 //
//                 fftw_Plans                      //
//                                                 //
/////////////////////////////////////////////////////

// fftw plans constructor - creates the plans
void fftwPlans::calculatePlans( int MF_size,  int Ckl_size) {
    
    
    ////////////////////////////////////////////////
    //      CREATE TEMPORARY EIGEN OBJECTS FOR PLANS
    //  I'm relatively sure it doesn't create fftw problems if you later delete these...
    
    
    //  THIS MIGHT BE CAUSING THE BUG!!! IF THE PLAN FUNCTIONS REFERENCE THE POINTER
    // TO THE DATA (WHICH THEY OBVIOUSLY DO!), AND THEN THAT DATA GOES OUT OF SCOPE
    // (WHICH IT DOES), THEN IF THE fftw_destroy_plan FUNCTION REFERENCES THE
    // POINTER DATA IT MIGHT NOT WORK.
    // BUT WHY WOULD THIS CAUSE A BUG WHEN a IS DESTROYED??
    
    
    
    // Ckl temp
    dcmplxMat  Ckl_tmp(Ckl_size,Ckl_size);
    dcmplx * Ckl_p = Ckl_tmp.data();
    
    // MF vector temps - need in and out separately
    dcmplxVec MF_tmp1( MF_size );
    dcmplx * MF_p1 = MF_tmp1.data();
    dcmplxVec MF_tmp2( MF_size );
    dcmplx * MF_p2 = MF_tmp2.data();
    // MF matrix temps - need in and out separately
    dcmplxMat  MFmat_tmp1(MF_size,MF_size);
    dcmplx * MFmat_p1 = MFmat_tmp1.data();
    dcmplxMat  MFmat_tmp2(MF_size,MF_size);
    dcmplx * MFmat_p2 = MFmat_tmp2.data();
    
    ////////////////////////////////////////////////
    //     2D PLANS - full data for Reynold's stress
    //
    // Need fft(fft( Ckl )')' (in Matlab notation), while standard 2-D fft is
    // fft(fft( Ckl ).').', hence the separate calculation for each dimension.
    //
    // These transforms can be in-place

    // Something funky going on - memory leak of some sort
    
    C2D_back_dim2_ = fftw_plan_many_dft(1, &Ckl_size, Ckl_size,
                                        CAST_T0_FFTW(Ckl_p), NULL, Ckl_size, 1,
                                        CAST_T0_FFTW(Ckl_p), NULL, Ckl_size, 1,
                                        FFTW_FORWARD, MY_FFTWPLAN);
    C2D_for_dim1_ = fftw_plan_many_dft(1, &Ckl_size, Ckl_size,
                                       CAST_T0_FFTW(Ckl_p), NULL, 1, Ckl_size,
                                       CAST_T0_FFTW(Ckl_p), NULL, 1, Ckl_size,
                                       FFTW_FORWARD, MY_FFTWPLAN);
    C2D_for_dim2_ = fftw_plan_many_dft(1, &Ckl_size, Ckl_size,
                                       CAST_T0_FFTW(Ckl_p), NULL, Ckl_size, 1,
                                       CAST_T0_FFTW(Ckl_p), NULL, Ckl_size, 1,
                                       FFTW_FORWARD, MY_FFTWPLAN);
    C2D_back_dim1_ = fftw_plan_many_dft(1, &Ckl_size, Ckl_size,
                                        CAST_T0_FFTW(Ckl_p), NULL, 1, Ckl_size,
                                        CAST_T0_FFTW(Ckl_p), NULL, 1, Ckl_size,
                                        FFTW_FORWARD, MY_FFTWPLAN);

    
    ////////////////////////////////////////////////
    //      MF plans
    //  Need both a 1-D plan and a 1-D*NZ plan
    //
    // Useful to have both in-place and out of place plans available
   
    // True 1-D plans
    MF1D_for_ = fftw_plan_dft_1d(MF_size, CAST_T0_FFTW(MF_p1), CAST_T0_FFTW(MF_p2),
                                 FFTW_FORWARD, MY_FFTWPLAN);
    MF1D_back_ = fftw_plan_dft_1d(MF_size, CAST_T0_FFTW(MF_p1), CAST_T0_FFTW(MF_p2),
                                  FFTW_BACKWARD, MY_FFTWPLAN);
    MF1D_IP_for_ = fftw_plan_dft_1d(MF_size, CAST_T0_FFTW(MF_p1), CAST_T0_FFTW(MF_p1),
                                 FFTW_FORWARD, MY_FFTWPLAN);
    MF1D_IP_back_ = fftw_plan_dft_1d(MF_size, CAST_T0_FFTW(MF_p1), CAST_T0_FFTW(MF_p1),
                                  FFTW_BACKWARD, MY_FFTWPLAN);
    
    // 1-D*NZ plans
    MF2D_for_ = fftw_plan_many_dft(1, &MF_size, MF_size,
                                   CAST_T0_FFTW(MFmat_p1), NULL, 1, MF_size,
                                   CAST_T0_FFTW(MFmat_p2), NULL, 1, MF_size,
                                   FFTW_FORWARD, MY_FFTWPLAN);
    MF2D_back_ = fftw_plan_many_dft(1, &MF_size, MF_size,
                                   CAST_T0_FFTW(MFmat_p1), NULL, 1, MF_size,
                                   CAST_T0_FFTW(MFmat_p2), NULL, 1, MF_size,
                                   FFTW_BACKWARD, MY_FFTWPLAN);
    MF2D_IP_for_ = fftw_plan_many_dft(1, &MF_size, MF_size,
                                   CAST_T0_FFTW(MFmat_p1), NULL, 1, MF_size,
                                   CAST_T0_FFTW(MFmat_p1), NULL, 1, MF_size,
                                   FFTW_FORWARD, MY_FFTWPLAN);
    MF2D_IP_back_ = fftw_plan_many_dft(1, &MF_size, MF_size,
                                    CAST_T0_FFTW(MFmat_p1), NULL, 1, MF_size,
                                    CAST_T0_FFTW(MFmat_p1), NULL, 1, MF_size,
                                    FFTW_BACKWARD, MY_FFTWPLAN);
    
    // Flag to use in destructor
    plans_calculated_ = 1;
    
}

// Destructor
fftwPlans::~fftwPlans() {
    ///////////////////////////////////////////////////////////////////////
    // This deletes data that has been created in fftwPlans::calculatePlans

   // std::cout << "Destroying " << std::endl;
    if (plans_calculated_) {
        // Ckl plans
        fftw_destroy_plan(C2D_for_dim1_);
        fftw_destroy_plan(C2D_for_dim2_);
        fftw_destroy_plan(C2D_back_dim1_);
        fftw_destroy_plan(C2D_back_dim2_);
        
        // MF plans
        fftw_destroy_plan(MF1D_for_);
        fftw_destroy_plan(MF1D_back_);
        fftw_destroy_plan(MF1D_IP_for_);
        fftw_destroy_plan(MF1D_IP_back_);
        fftw_destroy_plan(MF2D_for_);
        fftw_destroy_plan(MF2D_back_);
        fftw_destroy_plan(MF2D_IP_for_);
        fftw_destroy_plan(MF2D_IP_back_);
    }
    
    fftw_cleanup();
    
}

//                                                 //
/////////////////////////////////////////////////////





/////////////////////////////////////////////////////
//                                                 //
//                 TimeVariables                   //
//                                                 //
/////////////////////////////////////////////////////

// Current position
int TimeVariables::curr_pos_ = 0;

// Constructor
TimeVariables::TimeVariables(int nsteps, int width, int mpinode) :
nsteps_(nsteps), width_(width), mpi_node_(mpinode)
{
    // Main variables
    // Store only on 1-processor - presumbably this is better
    if (mpi_node_ == 0) {
        energy_ = new double*[nsteps];
        angular_momentum_ = new double*[nsteps];
        dissipation_ = new double*[nsteps];
        for (int i=0; i<nsteps; ++i) {
            energy_[i] = new double[width_];
            angular_momentum_[i] = new double[width_];
            dissipation_[i] = new double[width_];
        }
    }
}

// Destructor
TimeVariables::~TimeVariables() {
    if (mpi_node_==0) {
        for (int i=0; i<nsteps_; ++i) {
            delete[] energy_[i];
            delete[] angular_momentum_[i];
            delete[] dissipation_[i];
        }
    }
    delete[] energy_;
    delete[] angular_momentum_;
    delete[] dissipation_;
    
}

//                                                 //
/////////////////////////////////////////////////////

