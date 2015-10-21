//
//  intRK2CN.h
//  MRIDSS
//
//  Created by Jonathan Squire on 6/5/14.
//  Copyright (c) 2014 Jonathan Squire. All rights reserved.
//

#ifndef __MRIDSS__intRK2CN__
#define __MRIDSS__intRK2CN__

#include "Integrator.h"


class Model;

class RK2CN : public Integrator {
public:
    // Constructor and destructor
    RK2CN(double t0, Inputs& SP, Model &model);
    ~RK2CN();
    // Main step
    double Step(double t, dcmplxVec* MF, dcmplxMat * Ckl);
    // Reinitialize linear operators
    void Reinitialize_linear_Ops(double t);
    // Average time step
    double mean_time_step() const {return dt_mean_/step_count_;};
private:
    const int dim_Ckl_array_;   // dimension of Ckl array of matrices (Nx*Ny)
    const int size_Ckl_;  //size of individual matrices (4*Nz)
    const int num_MF_;   //  Number of mean fields
    const int size_MF_;  // size of Mean fields (Nz)
    
    // Timestepping
    double dt_;      // timestep
    bool variable_dt_;
    double CFLnum_, dtmax_,dt_mean_;
    int step_count_;
    
    // Model
    Model &model_;    // functor to evaluate f(x,t)
    
    // SPACE FOR INTEGRATOR EVALUATION
    dcmplxMat *Ckl_rhs_, *Ckl_th2_;     // temporary space to hold rhs of CKL equation
    dcmplxVec *MF_rhs_, *MF_th2_;     // Temporary space to hold rhs of MF equation
    // Space to hold linear (diffusion) part of operator
    // Double is much quicker than dcmplx if possible
    doubVec * linop_Ckl_;
    doubVec * linop_Ckl_old_; // Storage from previous step
    
    doubVec * linop_MF_; // THIS ASSUMES linop_MF IS CONSTANT IN TIME
    // Coefficients for CN equation
    doubVec * linop_MF_linCo_dt_, * linop_MF_NLCo_dt_;// For the second RK step
    doubVec * linop_MF_linCo_dto2_, * linop_MF_NLCo_dto2_; // For the first RK step (Euler CN with dt/2)
    
    
};
#endif /* defined(__MRIDSS__intRK2CN__) */
