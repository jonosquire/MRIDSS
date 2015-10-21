//
//  RK3CN.cpp
//  QL_DNS
//
//  Created by Jonathan Squire on 9/17/14.
//  Copyright (c) 2014 J Squire. All rights reserved.
//

// Third order RK time integrator as in Lundblahd et al. 1992
// Checked error(dt) using energy and it seems to be perfectly 3rd order in time
// Not sure about where the linear operators should be evaluated though, need to work this out!!
#include "intRK3CN.h"


RK3CN::RK3CN(double t0, Inputs& SP, Model &model) :
dim_Ckl_array_(model.Cdimxy()), size_Ckl_(model.Cdimz()),
num_MF_(model.num_MFs()), size_MF_(model.MFdimz()),
model_(model),
a0_(8.0/15.0),a1_(5.0/12.0),a2_(3.0/4.0),b1_(-17.0/60.0),b2_(-5.0/12.0),
p1_(8.0/15.0),p2_(2.0/3.0),// Integrator parameters
step_count_(0), dt_mean_(0.0)
{
    if (SP.dt < 0 ){
        variable_dt_ = 1;
        CFLnum_ = SP.CFL;
        dtmax_ = -SP.dt;
    } else {
        dt_ = SP.dt;
        variable_dt_ = 0;
    }
    
    
    // Assigning necessary temporary space
    MF_rhs_ =new dcmplxVec[num_MF_];
    MF_rhs2_ =new dcmplxVec[num_MF_];
    for (int i=0; i<num_MF_; i++){
        MF_rhs_[i]= dcmplxVec(size_MF_);
        MF_rhs2_[i]= dcmplxVec(size_MF_);
    }
    
    // C matrices
    Ckl_rhs_ =new dcmplxMat[dim_Ckl_array_];
    Ckl_rhs2_ = new dcmplxMat[dim_Ckl_array_];
    for (int i=0; i<dim_Ckl_array_; i++) {// Split across processors
        Ckl_rhs_[i]= dcmplxMat(size_Ckl_,size_Ckl_);
        Ckl_rhs2_[i]= dcmplxMat(size_Ckl_,size_Ckl_);
    }

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
    // Assume MF operator constant in time 
    linop_MF_ = new doubVec[num_MF_];
    for (int i=0; i<num_MF_; i++){
        linop_MF_[i]= doubVec(size_MF_);
    }
    denomMF_ = doubVec(size_MF_);
    model_.linearOPs_Init(t0, linop_MF_, linop_Ckl_old_);
    
}

RK3CN::~RK3CN() {
    // Linear operators
    // C
    delete[] linop_Ckl_;
    delete[] linop_Ckl_old_;
    // MF
    delete[] linop_MF_;
    
    // NL parts
    delete[] Ckl_rhs_;
    delete[] MF_rhs_;
    delete[] Ckl_rhs2_;
    delete[] MF_rhs2_;
}


// Reinitialize linear operators
void RK3CN::Reinitialize_linear_Ops(double t){
    linop_MF_ = new doubVec[num_MF_];
    for (int i=0; i<num_MF_; i++){
        linop_MF_[i]= doubVec(size_MF_);
    }
    model_.linearOPs_Init(t, linop_MF_, linop_Ckl_old_);
}


///////////////////////////////////////////
//   STEP
// Step forward in time
double RK3CN::Step(double t, dcmplxVec* MF, dcmplxMat * Ckl){
    
    if (variable_dt_) {
        dt_ = CFLnum_/model_.Calculate_CFL(MF, Ckl);
        if (dt_ > dtmax_)
            dt_ = dtmax_;
        
    }
    dt_mean_ += dt_;
    ++step_count_;
    
    
    // Temporary pointer variables
    dcmplx * dataC, * dataC_rhs, *dataC_rhs2;
    double * data_linopC, *data_linopC_old;
    double denom;
    
    // Still not SURE linear part is correct on this, should check...
    //////////////////////////////////////////////////////
    //   STEP 1
    linC_ =dt_/2*(a0_);
    NLC1_ = a0_*dt_;
    lin_dt_ = dt_*p1_;
    
    model_.rhs(t, lin_dt_,  MF, Ckl, MF_rhs_, Ckl_rhs_,linop_Ckl_);
    
    // Linear fluctuations
    for (int i=0; i<dim_Ckl_array_; i++) {// Loop over kx ky
        // Get pointers to various data sets - saves calcuating full mat for lin_op
        dataC = Ckl[i].data();  // This is not re-written
        dataC_rhs = Ckl_rhs_[i].data();
        data_linopC = linop_Ckl_[i].data(); // At time t+dt/2
        data_linopC_old = linop_Ckl_old_[i].data(); //At time t
        
        for (int j=0; j<size_Ckl_; ++j) {
            for (int k=0; k<size_Ckl_; ++k) {
                denom = 1.0/(1.0-linC_*(data_linopC[k]+data_linopC[j]));
                // Column/Row major order makes no difference here
                dataC[k+size_Ckl_*j] = NLC1_*denom*dataC_rhs[k+size_Ckl_*j] +
                    (1+linC_*(data_linopC_old[k]+data_linopC_old[j]))*denom*dataC[k+size_Ckl_*j];
            }
        }
    }
    
    // Mean fields
    for (int i=0; i<num_MF_; ++i) {
        denomMF_ = 1.0/(1.0-linC_*linop_MF_[i]);
        MF[i] = NLC1_*denomMF_*MF_rhs_[i] + (1+linC_*linop_MF_[i])*denomMF_*MF[i];
    }

    
    // Move linear operator forward variables around
    for (int i=0; i<dim_Ckl_array_; i++)
        linop_Ckl_old_[i] = linop_Ckl_[i];
    //    END - STEP 1
    //////////////////////////////////////////////////////
    
    
    
    //////////////////////////////////////////////////////
    //   STEP 2 - Updates sol
    linC_ =dt_/2*(a1_+b1_);
    NLC1_ = a1_*dt_;
    NLC2_ = b1_*dt_;
    lin_dt_ = dt_*(p2_-p1_);
    
    model_.rhs(t+dt_*p1_, lin_dt_,  MF, Ckl, MF_rhs2_, Ckl_rhs2_,linop_Ckl_);
    
    // Linear fluctuations
    for (int i=0; i<dim_Ckl_array_; i++) {// Loop over kx ky
        // Get pointers to various data sets - saves calcuating full mat for lin_op
        dataC = Ckl[i].data();  // This is not re-written
        dataC_rhs = Ckl_rhs_[i].data();
        dataC_rhs2 = Ckl_rhs2_[i].data();
        data_linopC = linop_Ckl_[i].data(); // At time t+dt/2
        data_linopC_old = linop_Ckl_old_[i].data(); //At time t
        
        for (int j=0; j<size_Ckl_; ++j) {
            for (int k=0; k<size_Ckl_; ++k) {
                denom = 1.0/(1.0-linC_*(data_linopC[k]+data_linopC[j]));
                // Column/Row major order makes no difference here
                dataC[k+size_Ckl_*j] = NLC1_*denom*dataC_rhs2[k+size_Ckl_*j] + NLC2_*denom*dataC_rhs[k+size_Ckl_*j]  +  (1+linC_*(data_linopC_old[k]+data_linopC_old[j]))*denom*dataC[k+size_Ckl_*j];
            }
        }
    }
    
    // Mean fields
    for (int i=0; i<num_MF_; ++i) {
        denomMF_ = 1.0/(1.0-linC_*linop_MF_[i]);
        MF[i] = NLC1_*denomMF_*MF_rhs2_[i] +  NLC2_*denomMF_*MF_rhs_[i] + (1+linC_*linop_MF_[i])*denomMF_*MF[i];
    }
    
    
    // Move linear operator forward variables around
    for (int i=0; i<dim_Ckl_array_; i++)
        linop_Ckl_old_[i] = linop_Ckl_[i];

    //    END - STEP 2
    //////////////////////////////////////////////////////
    
    
    
    
    //////////////////////////////////////////////////////
    //   STEP 3 - Updates sol
    linC_ =dt_/2*(a2_+b2_);
    NLC1_ = a2_*dt_;
    NLC2_ = b2_*dt_;
    lin_dt_ = dt_*(1.0-p2_);
    
    model_.rhs(t+dt_*p2_, lin_dt_,  MF, Ckl, MF_rhs_, Ckl_rhs_,linop_Ckl_);
    
    // Linear fluctuations
    for (int i=0; i<dim_Ckl_array_; i++) {// Loop over kx ky
        // Get pointers to various data sets - saves calcuating full mat for lin_op
        dataC = Ckl[i].data();  // This is not re-written
        dataC_rhs = Ckl_rhs_[i].data();
        dataC_rhs2 = Ckl_rhs2_[i].data();
        data_linopC = linop_Ckl_[i].data(); // At time t+dt/2
        data_linopC_old = linop_Ckl_old_[i].data(); //At time t
        
        for (int j=0; j<size_Ckl_; ++j) {
            for (int k=0; k<size_Ckl_; ++k) {
                denom = 1.0/(1.0-linC_*(data_linopC[k]+data_linopC[j]));
                // Column/Row major order makes no difference here
                dataC[k+size_Ckl_*j] = NLC1_*denom*dataC_rhs[k+size_Ckl_*j] + NLC2_*denom*dataC_rhs2[k+size_Ckl_*j]  +  (1+linC_*(data_linopC_old[k]+data_linopC_old[j]))*denom*dataC[k+size_Ckl_*j];
            }
        }
    }
    
    // Mean fields
    for (int i=0; i<num_MF_; ++i) {
        denomMF_ = 1.0/(1.0-linC_*linop_MF_[i]);
        MF[i] = NLC1_*denomMF_*MF_rhs_[i] +  NLC2_*denomMF_*MF_rhs2_[i] + (1+linC_*linop_MF_[i])*denomMF_*MF[i];
    }
    
    
    // Move linear operator forward variables around
    for (int i=0; i<dim_Ckl_array_; i++)
        linop_Ckl_old_[i] = linop_Ckl_[i];
    
    //    END - STEP 3
    //////////////////////////////////////////////////////
    
    return dt_;
    
    
    
}