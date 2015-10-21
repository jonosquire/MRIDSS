//
//  intRK2CN.cc
//  MRIDSS
//
//  Created by Jonathan Squire on 6/5/14.
//  Copyright (c) 2014 Jonathan Squire. All rights reserved.
//

#include "intRK2CN.h"


RK2CN::RK2CN(double t0, Inputs& SP, Model &model) :
dim_Ckl_array_(model.Cdimxy()), size_Ckl_(model.Cdimz()),
num_MF_(model.num_MFs()), size_MF_(model.MFdimz()),
model_(model),
step_count_(0), dt_mean_(0.0)
{
    if (SP.dt<0)
        std::cout << "Variable time-step not supported by RK2CN integrator!"<< std::endl;
    dt_ = fabs(SP.dt); // This integrator requires a time-step to be specified in inputs!!
    variable_dt_ = 0;
    
    // Assigning necessary temporary space
    MF_rhs_ =new dcmplxVec[num_MF_];
    MF_th2_ =new dcmplxVec[num_MF_];
    for (int i=0; i<num_MF_; i++){
        MF_rhs_[i]= dcmplxVec(size_MF_);
        MF_th2_[i]= dcmplxVec(size_MF_);
    }
    
    // C matrices
    Ckl_rhs_ =new dcmplxMat[dim_Ckl_array_];
    Ckl_th2_ = new dcmplxMat[dim_Ckl_array_];
    for (int i=0; i<dim_Ckl_array_; i++) {// Split across processors
        Ckl_rhs_[i]= dcmplxMat(size_Ckl_,size_Ckl_);
        Ckl_th2_[i]= dcmplxMat(size_Ckl_,size_Ckl_);
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
    
    /////////////////////////////////////
    // INITIALIZING LINEAR OPERATORS - MF is constant in time
    linop_MF_ = new doubVec[num_MF_];
    linop_MF_linCo_dt_ = new doubVec[num_MF_];
    linop_MF_NLCo_dt_ = new doubVec[num_MF_];
    linop_MF_linCo_dto2_ = new doubVec[num_MF_];
    linop_MF_NLCo_dto2_ = new doubVec[num_MF_];
    for (int i=0; i<num_MF_; i++){
        linop_MF_[i]= doubVec(size_MF_);
        linop_MF_linCo_dt_[i]= doubVec(size_MF_);
        linop_MF_NLCo_dt_[i]= doubVec(size_MF_);
        linop_MF_linCo_dto2_[i]= doubVec(size_MF_);
        linop_MF_NLCo_dto2_[i]= doubVec(size_MF_);
    }
    model_.linearOPs_Init(t0, linop_MF_, linop_Ckl_old_);
    for (int i=0; i<num_MF_; i++){
        linop_MF_NLCo_dt_[i] = dt_/(1.0-dt_/2*linop_MF_[i]);
        linop_MF_linCo_dt_[i] = (1.0+dt_/2*linop_MF_[i])/(1.0-dt_/2*linop_MF_[i]);
        linop_MF_NLCo_dto2_[i] = dt_/2/(1.0-dt_/4*linop_MF_[i]);
        linop_MF_linCo_dto2_[i] = (1.0+dt_/4*linop_MF_[i])/(1.0-dt_/4*linop_MF_[i]);
    }
    delete[] linop_MF_; // No longer required for any reason
    
    
    
}

RK2CN::~RK2CN() {
    // Linear operators
    // C
    delete[] linop_Ckl_;
    delete[] linop_Ckl_old_;
    // MF
    delete[] linop_MF_linCo_dt_;
    delete[] linop_MF_NLCo_dt_;
    delete[] linop_MF_linCo_dto2_;
    delete[] linop_MF_NLCo_dto2_;
    
    delete[] Ckl_rhs_;
    delete[] MF_rhs_;
    delete[] Ckl_th2_;
    delete[] MF_th2_;

    
}

// Reinitialize linear operators
void RK2CN::Reinitialize_linear_Ops(double t){
    linop_MF_ = new doubVec[num_MF_];
    for (int i=0; i<num_MF_; i++){
        linop_MF_[i]= doubVec(size_MF_);
    }
    model_.linearOPs_Init(t, linop_MF_, linop_Ckl_old_);
    for (int i=0; i<num_MF_; i++){
        linop_MF_NLCo_dt_[i] = dt_/(1.0-dt_/2*linop_MF_[i]);
        linop_MF_linCo_dt_[i] = (1.0+dt_/2*linop_MF_[i])/(1.0-dt_/2*linop_MF_[i]);
        linop_MF_NLCo_dto2_[i] = dt_/2/(1.0-dt_/4*linop_MF_[i]);
        linop_MF_linCo_dto2_[i] = (1.0+dt_/4*linop_MF_[i])/(1.0-dt_/4*linop_MF_[i]);
    }
    delete[] linop_MF_;
}

double RK2CN::Step(double t, dcmplxVec* MF, dcmplxMat * Ckl) {
    
    dt_mean_ += dt_;
    ++step_count_;
    
    // Temporary pointer variables
    dcmplx * dataC, * dataC_rhs, *dataC_th2;
    double * data_linopC, *data_linopC_old;
    double denom, dt_tmp;// Temporary storage for loop
    
    //////////////////////////////////////////////////////////
    //                  FIRST STEP                          //
    dt_tmp = dt_/2;
    model_.rhs(t, dt_tmp, MF, Ckl, MF_rhs_, Ckl_rhs_,linop_Ckl_);
    
    // WARNING: SECOND SUBSTEP OF THIS INTEGRATOR REQUIRES INTEGRATOR TO WORK IF in AND out OBJECTS ARE THE SAME. THIS CAN BE A LITTLE TRICKY. HAVE TO BE CAREFUL!!
    // MF_rhs and Ckl_rhs are MF(t+h/2) and Ckl(t+h/2) after this step - saves storage
    for (int i=0; i<dim_Ckl_array_; i++) {
        // Get pointers to various data sets - saves calcuating full mat for lin_op
        dataC = Ckl[i].data();  // This is not re-written
        dataC_rhs = Ckl_rhs_[i].data();
        dataC_th2 = Ckl_th2_[i].data();
        data_linopC = linop_Ckl_[i].data(); // At time t+dt/2
        data_linopC_old = linop_Ckl_old_[i].data(); //At time t
        // NB: This method is very nearly as fast as using a pure pointer array (ie. native C++). The data() method seems to have very little overhead as expected
        
        // Replace Ckl_rhs with (1+dt/4*L)/(1-dt/4*L)*Ckl + dt/2/(1-dt/4*L)*Ckl_rhs
        for (int j=0; j<size_Ckl_; ++j) {
            for (int k=0; k<size_Ckl_; ++k) {
                denom = 1/(1-dt_tmp/2*(data_linopC[k]+data_linopC[j]));
                // Column/Row major order makes no difference here
                dataC_th2[k+size_Ckl_*j]=(1+dt_tmp/2*(data_linopC_old[k]+data_linopC_old[j]))
                        *denom*dataC[k+size_Ckl_*j]  +
                    dt_tmp*denom*dataC_rhs[k+size_Ckl_*j];
            }
        }
    }
    
    for (int i=0; i<num_MF_; i++) {
        MF_th2_[i] = linop_MF_NLCo_dto2_[i]*MF_rhs_[i] + linop_MF_linCo_dto2_[i]*MF[i];
    }
    
    
    //////////////////////////////////////////////////////////
    //                  SECOND STEP                          //
    dt_tmp = dt_;
    model_.rhs(t+dt_/2, dt_/2, MF_th2_, Ckl_th2_, MF_rhs_, Ckl_rhs_,linop_Ckl_);
    
    // NOTE data_linopC_old is still at t
    for (int i=0; i<dim_Ckl_array_; i++) {
        // Get pointers to various data sets - saves calcuating full mat for lin_op
        dataC = Ckl[i].data();
        dataC_rhs = Ckl_rhs_[i].data();
        data_linopC = linop_Ckl_[i].data(); //Evaluated at t+dt
        data_linopC_old = linop_Ckl_old_[i].data(); // At time t
        // NB: This method is very nearly as fast as using a pure pointer array (ie. native C++). The data() method seems to have very little overhead as expected
        
        // Replace Ckl_ with (1+dt/2*L)/(1-dt/2*L)*Ckl + dt/(1-dt/2*L)*Ckl_rhs
        for (int j=0; j<size_Ckl_; ++j) {
            for (int k=0; k<size_Ckl_; ++k) {
                denom = 1/(1-dt_tmp/2*(data_linopC[k]+data_linopC[j]));
                // Column/Row major order makes no difference here
                dataC[k+size_Ckl_*j]=(1+dt_tmp/2*(data_linopC_old[k]+data_linopC_old[j]))
                        *denom*dataC[k+size_Ckl_*j]  +
                    dt_tmp*denom*dataC_rhs[k+size_Ckl_*j];
            }
        }
    }
    
    for (int i=0; i<num_MF_; i++) {
        MF[i] = linop_MF_NLCo_dt_[i]*MF_rhs_[i] + linop_MF_linCo_dt_[i]*MF[i];
    }
    
    // Move variables around
    for (int i=0; i<dim_Ckl_array_; i++)
        linop_Ckl_old_[i] = linop_Ckl_[i];
    

    
    return dt_;
}

