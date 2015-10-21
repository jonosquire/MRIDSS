//
//  RK3CN.h
//  QL_DNS
//
//  Created by Jonathan Squire on 9/17/14.
//  Copyright (c) 2014 J Squire. All rights reserved.
//

#ifndef __QL_DNS__RK3CN__
#define __QL_DNS__RK3CN__

#include "Integrator.h"

class Model;

class RK3CN : public Integrator {
public:
    RK3CN(double t0, Inputs& SP, Model & model);
    ~RK3CN();
    // Main step
    double Step(double t, dcmplxVec* MF, dcmplxMat * Ckl);
    // Reinitialize linear operators
    void Reinitialize_linear_Ops(double t);
    // Average time step
    double mean_time_step() const {return dt_mean_/step_count_;};
private:
    // Time-step
    double dt_, lin_dt_;
    double CFLnum_, dtmax_,dt_mean_;
    int step_count_;
    bool variable_dt_;
    
    const int dim_Ckl_array_;   // dimension of Ckl array of matrices (Nx*Ny)
    const int size_Ckl_;  //size of individual matrices (4*Nz)
    const int num_MF_;   //  Number of mean fields
    const int size_MF_;  // size of Mean fields (Nz)
    
    // Integrator parameters
    const double a0_,a1_,a2_,b1_,b2_,p1_,p2_;
    double linC_, NLC1_, NLC2_; // Coefficients for each step
    
    
    // Model
    Model& model_;
    
    // SPACE FOR INTEGRATOR EVALUATION
    dcmplxMat *Ckl_rhs_, *Ckl_rhs2_;     // temporary space to hold rhs of CKL equation
    dcmplxVec *MF_rhs_, *MF_rhs2_;     // Temporary space to hold rhs of MF equation
    // Space to hold linear (diffusion) part of operator
    // Double is much quicker than dcmplx if possible
    doubVec *linop_Ckl_;
    doubVec *linop_Ckl_old_; // Storage from previous step
    
    doubVec *linop_MF_; // THIS ASSUMES linop_MF IS CONSTANT IN TIME
    doubVec denomMF_; // Used for mean denominator
    //
};




#endif /* defined(__QL_DNS__RK3CN__) */


