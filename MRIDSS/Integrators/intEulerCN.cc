//
//  intEulerCN.cc
//  MRIDSS
//
//  Created by Jonathan Squire on 28/04/14.
//  Copyright (c) 2014 J Squire. All rights reserved.
//

#include "intEulerCN.h"

EulerCN::EulerCN(double t0, Inputs& SP, Model &model) :
dim_Ckl_array_(model.Cdimxy()), size_Ckl_(model.Cdimz()),
num_MF_(model.num_MFs()), size_MF_(model.MFdimz()),
model_(model),
step_count_(0), dt_mean_(0.0)
{
    if (SP.dt<0)
        std::cout << "Variable time-step not supported by EulerCN integrator!"<< std::endl;
    dt_ = fabs(SP.dt); // This integrator requires a time-step to be specified in inputs!!
    variable_dt_ = 0;
    
    // Assigning necessary temporary space
    MF_new_ =new dcmplxVec[num_MF_];
    for (int i=0; i<num_MF_; i++)
        MF_new_[i]= dcmplxVec(size_MF_);
    
    Ckl_new_ =new dcmplxMat[dim_Ckl_array_];
    for (int i=0; i<dim_Ckl_array_; i++) // MPI -- Split
        Ckl_new_[i]= dcmplxMat(size_Ckl_,size_Ckl_);
    
    //////////////////////////////////////
    ////////   LINEAR OPERATORS //////////
    //////////////////////////////////////
  
    // Space for diffusion operators
    linop_Ckl_ = new doubVec[dim_Ckl_array_];
    linop_Ckl_old_ = new doubVec[dim_Ckl_array_];
    for (int i=0; i<dim_Ckl_array_; i++){
        linop_Ckl_[i]= doubVec::Zero(size_Ckl_);//Diagonal since diffusion
        linop_Ckl_old_[i]= doubVec(size_Ckl_);
    }
    
    /////////////////////////////////////
    // INITIALIZING LINEAR OPERATORS - MF is constant in time
    linop_MF_ = new doubVec[num_MF_];
    linop_MF_linCo_ = new doubVec[num_MF_];
    linop_MF_NLCo_ = new doubVec[num_MF_];
    for (int i=0; i<num_MF_; i++){
        linop_MF_[i]= doubVec(size_MF_);
        linop_MF_linCo_[i]= doubVec(size_MF_);
        linop_MF_NLCo_[i]= doubVec(size_MF_);
    }
    model_.linearOPs_Init(t0, linop_MF_, linop_Ckl_old_);
    for (int i=0; i<num_MF_; i++){
        linop_MF_NLCo_[i] = dt_/(1.0-dt_/2*linop_MF_[i]);
        linop_MF_linCo_[i] = (1.0+dt_/2*linop_MF_[i])/(1.0-dt_/2*linop_MF_[i]);
    }
    
    
    
}

EulerCN::~EulerCN() {
    delete[] linop_Ckl_;
    delete[] linop_Ckl_old_;
    delete[] linop_MF_;
    delete[] linop_MF_linCo_;
    delete[] linop_MF_NLCo_;
    delete[] Ckl_new_;
    delete[] MF_new_;
    
}

// Reinitialize linear operators
void EulerCN::Reinitialize_linear_Ops(double t){
    
    model_.linearOPs_Init(t, linop_MF_, linop_Ckl_old_);
    for (int i=0; i<num_MF_; i++){
        linop_MF_NLCo_[i] = dt_/(1.0-dt_/2*linop_MF_[i]);
        linop_MF_linCo_[i] = (1.0+dt_/2*linop_MF_[i])/(1.0-dt_/2*linop_MF_[i]);
    }
    
}


double EulerCN::Step(double t, dcmplxVec* MF, dcmplxMat * Ckl) {
    
    
    dt_mean_ += dt_;
    ++step_count_;
    
    
    model_.rhs(t, dt_, MF, Ckl, MF_new_, Ckl_new_,linop_Ckl_);
    
    dcmplx * dataC, * dataC_new;
    double * data_linopC, *data_linopC_old;
    double denom;// Temporary storage for loop
    for (int i=0; i<dim_Ckl_array_; i++) {
        // Get pointers to various data sets - saves calcuating full mat for lin_op
        dataC = Ckl[i].data();
        dataC_new = Ckl_new_[i].data();
        data_linopC = linop_Ckl_[i].data();
        data_linopC_old = linop_Ckl_old_[i].data();
        // NB: This method is very nearly as fast as using a pure pointer array (ie. native C++). The data() method seems to have very little overhead as expected
        
        // Replace Ckl_new with dt/(1-dt*L)*Ckl_new
        // and Replace Ckl with (1+dt*L)/(1-dt*L)*Ckl
        for (int j=0; j<size_Ckl_; ++j) {
            for (int k=0; k<size_Ckl_; ++k) {
                denom = 1/(1-dt_/2*(data_linopC[k]+data_linopC[j]));
                // Column/Row major order makes no difference here
                dataC[k+size_Ckl_*j]=(1+dt_/2*(data_linopC_old[k]+data_linopC_old[j]))*denom*dataC[k+size_Ckl_*j];
                dataC_new[k+size_Ckl_*j]=dt_*denom*dataC_new[k+size_Ckl_*j];
            }
        }
        
        // Euler Crank Nicholson integrator step
        Ckl[i] += Ckl_new_[i];
    }
    
    for (int i=0; i<num_MF_; i++) {
        MF[i] = linop_MF_NLCo_[i]*MF_new_[i] + linop_MF_linCo_[i]*MF[i];
    }
    
    // Move variables around
    for (int i=0; i<dim_Ckl_array_; i++)
        linop_Ckl_old_[i] = linop_Ckl_[i];
                       
                       
    return dt_;
}



