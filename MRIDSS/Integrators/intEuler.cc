//
//  intEuler.cc
//  TwoDFluid
//
//  Created by Jonathan Squire on 28/04/14.
//  Copyright (c) 2014 J Squire. All rights reserved.
//

#include "intEuler.h"

Euler::Euler(double t0,Inputs& SP, Model &model) :
dim_Ckl_array_(model.Cdimxy()), size_Ckl_(model.Cdimz()),
num_MF_(model.num_MFs()), size_MF_(model.MFdimz()),
model_(model),
step_count_(0), dt_mean_(0.0)
{
    if (SP.dt<0)
        std::cout << "Variable time-step not supported by Euler integrator!"<< std::endl;
    dt_ = abs(SP.dt); // This integrator requires a time-step to be specified in inputs!!
    
    
    // Assigning necessary temporary space
    MF_new_ =new dcmplxVec[num_MF_];
    for (int i=0; i<num_MF_; i++)
        MF_new_[i]= dcmplxVec(size_MF_);
    
    Ckl_new_ =new dcmplxMat[dim_Ckl_array_];
    for (int i=0; i<dim_Ckl_array_; i++) // MPI -- Split
        Ckl_new_[i]= dcmplxMat(size_Ckl_,size_Ckl_);
    
    
    
    // Space for diffusion operators
    linop_Ckl_ = new doubVec[dim_Ckl_array_];
    for (int i=0; i<dim_Ckl_array_; i++) // MPI -- Split
        linop_Ckl_[i]= doubVec::Zero(size_Ckl_);//Diagonal since diffusion
    
    // MF linear operator is constant in time
    linop_MF_ = new doubVec[num_MF_];
    for (int i=0; i<num_MF_; i++) // MPI -- Split
        linop_MF_[i]= doubVec(size_MF_);//Diagonal since diffusion
    model_.linearOPs_Init(t0,linop_MF_,linop_Ckl_);
    
    
}

Euler::~Euler() {
    delete[] linop_Ckl_;
    delete[] linop_MF_;
    delete[] Ckl_new_;
    delete[] MF_new_;

}

// Reinitialize linear operators
void Euler::Reinitialize_linear_Ops(double t){
    
    for (int i=0; i<num_MF_; i++) // MPI -- Split
        linop_MF_[i]= doubVec(size_MF_);//Diagonal since diffusion
    model_.linearOPs_Init(t,linop_MF_,linop_Ckl_);
}


double Euler::Step(double t, dcmplxVec* MF, dcmplxMat * Ckl) {
    
    dt_mean_ += dt_;
    ++step_count_;
    
    model_.rhs(t, 0, MF, Ckl, MF_new_, Ckl_new_,linop_Ckl_);
    
    for (int i=0; i<dim_Ckl_array_; i++) {
        Ckl[i].noalias() += dt_*(Ckl_new_[i]+
                                 (linop_Ckl_[i].matrix().asDiagonal()*Ckl[i]+Ckl[i]*linop_Ckl_[i].matrix().asDiagonal()) );
    }
    
    for (int i=0; i<num_MF_; i++) {
        MF[i] += dt_*(MF_new_[i] + linop_MF_[i]*MF[i]);
    }
    
    return dt_;
}
