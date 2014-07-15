//
//  fftwPlans.cc
//  MRIDSS
//
//  Created by Jonathan Squire on 6/25/14.
//  Copyright (c) 2014 Jonathan Squire. All rights reserved.
//

#include "fftwPlans.h"

/////////////////////////////////////////////////////
//                                                 //
//                 fftw_Plans                      //
//                                                 //
/////////////////////////////////////////////////////

// fftw plans constructor - creates the plans
void fftwPlans::calculatePlans( int NZ ) {
    
    
    ////////////////////////////////////////////////
    //      CREATE TEMPORARY EIGEN OBJECTS FOR PLANS
    //  I'm relatively sure it doesn't create fftw problems if you later delete these...
    
    // Ckl temp
    dcmplxMat  Ckl_tmp(NZ,NZ);
    dcmplx * Ckl_p = Ckl_tmp.data();
    
    // Vector temps
    dcmplxVec MF_tmp1( NZ );
    dcmplx * MF_p1 = MF_tmp1.data();
    

    ////////////////////////////////////////////////
    //     2D PLANS -  Reynold's stress
    //
    // Need fft(fft( Ckl )')' (in Matlab notation), while standard 2-D fft is
    // fft(fft( Ckl ).').', hence the separate calculation for each dimension.
    //
    // All transforms can be in-place
    
    
    
    t2D_for_dim1_ = fftw_plan_many_dft(1, &NZ, NZ,
                                       CAST_T0_FFTW(Ckl_p), NULL, 1, NZ,
                                       CAST_T0_FFTW(Ckl_p), NULL, 1, NZ,
                                       FFTW_FORWARD, MY_FFTWPLAN);
    t2D_for_dim2_ = fftw_plan_many_dft(1, &NZ, NZ,
                                       CAST_T0_FFTW(Ckl_p), NULL, NZ, 1,
                                       CAST_T0_FFTW(Ckl_p), NULL, NZ, 1,
                                       FFTW_FORWARD, MY_FFTWPLAN);
    t2D_back_dim1_ = fftw_plan_many_dft(1, &NZ, NZ,
                                        CAST_T0_FFTW(Ckl_p), NULL, 1, NZ,
                                        CAST_T0_FFTW(Ckl_p), NULL, 1, NZ,
                                        FFTW_BACKWARD, MY_FFTWPLAN);
    t2D_back_dim2_ = fftw_plan_many_dft(1, &NZ, NZ,
                                        CAST_T0_FFTW(Ckl_p), NULL, NZ, 1,
                                        CAST_T0_FFTW(Ckl_p), NULL, NZ, 1,
                                        FFTW_BACKWARD, MY_FFTWPLAN);
    
    
    ////////////////////////////////////////////////
    //   1D Plans
    //  Need both a 1-D plan and a 1-D*NZ plan
    //
    // Only in place transforms have been useful, could delete others
    
    // True 1-D plans
    t1D_for_ = fftw_plan_dft_1d(NZ, CAST_T0_FFTW(MF_p1), CAST_T0_FFTW(MF_p1),
                                FFTW_FORWARD, MY_FFTWPLAN);
    t1D_back_ = fftw_plan_dft_1d(NZ, CAST_T0_FFTW(MF_p1), CAST_T0_FFTW(MF_p1),
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
        fftw_destroy_plan(t2D_for_dim1_);
        fftw_destroy_plan(t2D_for_dim2_);
        fftw_destroy_plan(t2D_back_dim1_);
        fftw_destroy_plan(t2D_back_dim2_);
        
        // MF plans
        fftw_destroy_plan(t1D_for_);
        fftw_destroy_plan(t1D_back_);
    }
    
    fftw_cleanup();
    
}

// Full 2-D transform fft(fft(a)')')

void fftwPlans::for_2DFull(dcmplxMat& Cin){
    for_2D_dim1(Cin.data());
    Cin = Cin.conjugate();
    for_2D_dim2(Cin.data());
    Cin = Cin.conjugate();
}

void fftwPlans::back_2DFull(dcmplxMat& Cin){
    double nzm2 = 1.0/(Cin.size() );
    back_2D_dim1(Cin.data());
    Cin = Cin.conjugate();
    back_2D_dim2(Cin.data());
    Cin = Cin.conjugate()*nzm2;
}

// Full backwards transform without normalizing by NZ^2
void fftwPlans::back_2DFull_noNZ(dcmplxMat& Cin){
    back_2D_dim1(Cin.data());
    Cin = Cin.conjugate();
    back_2D_dim2(Cin.data());
    Cin = Cin.conjugate();
}
//                                                 //
/////////////////////////////////////////////////////

