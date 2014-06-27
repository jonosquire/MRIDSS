//
//  fftwPlans.h
//  MRIDSS
//
//  Created by Jonathan Squire on 6/25/14.
//  Copyright (c) 2014 Jonathan Squire. All rights reserved.
//

#ifndef __MRIDSS__fftwPlans__
#define __MRIDSS__fftwPlans__

#include "../General_Definitions.h"


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


#endif /* defined(__MRIDSS__fftwPlans__) */
