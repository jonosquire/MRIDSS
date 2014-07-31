//
//  intEulerCN.h
//  MRIDSS
//
//  Created by Jonathan Squire on 28/04/14.
//  Copyright (c) 2014 J Squire. All rights reserved.
//

#ifndef __MRIDSS__intEulerCN__
#define __MRIDSS__intEulerCN__

#include "Integrator.h"

class Model;

class EulerCN : public Integrator {
public:
    // Constructor and destructor
    EulerCN(double t0, double dt, Model &model);
    ~EulerCN();
    // Main step
    int Step(double t, dcmplxVec* MF, dcmplxMat * Ckl);
    // Reinitialize linear operators
    void Reinitialize_linear_Ops(double t);
private:
    const int dim_Ckl_array_;   // dimension of Ckl array of matrices (Nx*Ny)
    const int size_Ckl_;  //size of individual matrices (4*Nz)
    const int num_MF_;   //  Number of mean fields
    const int size_MF_;  // size of Mean fields (Nz)
    
    const double dt_;      // timestep
    Model &model_;    // functor to evaluate f(x,t)
    
    // SPACE FOR INTEGRATOR EVALUATION
    dcmplxMat *Ckl_new_;     // temporary space to hold new C solution
    dcmplxVec *MF_new_;     // Temporary space to hold new MF solution
    // Space to hold linear (diffusion) part of operator
    // Double is much quicker than dcmplx if possible
    doubVec * linop_Ckl_;
    doubVec * linop_Ckl_old_; // Storage from previous step
    
    doubVec * linop_MF_; // THIS ASSUMES linop_MF IS CONSTANT IN TIME
    doubVec * linop_MF_linCo_;
    doubVec * linop_MF_NLCo_;// Coefficients for CN equation

    
};


#endif /* defined(__MRIDSS__intEulerCN__) */
