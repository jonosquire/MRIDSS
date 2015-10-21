//
//  intEuler.h
//  TwoDFluid
//
//  Created by Jonathan Squire on 28/04/14.
//  Copyright (c) 2014 J Squire. All rights reserved.
//

#ifndef __MRIDSS__intEuler__
#define __MRIDSS__intEuler__


#include "Integrator.h"

class Model;

class Euler : public Integrator {
public:
    Euler(double t0, Inputs& SP, Model &model);
    ~Euler();
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
    double CFLnum_, dtmax_,dt_mean_;
    int step_count_;
    
    Model &model_;    // functor to evaluate f(x,t)
    
    // SPACE FOR INTEGRATOR EVALUATION
    dcmplxMat *Ckl_new_;     // temporary space to hold new C solution
    dcmplxVec *MF_new_;     // Temporary space to hold new MF solution
    // Space to hold linear (diffusion) part of operator
    doubVec * linop_Ckl_;
    doubVec * linop_MF_;  //Double is much quicker than dcmplx if possible
};


#endif /* defined(__TwoDFluid__intEuler__) */
