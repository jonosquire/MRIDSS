//
//  MHD_BQlin.h
//  MRIDSS
//
//  Created by Jonathan Squire on 4/25/14.
//  Copyright (c) 2014 J Squire. All rights reserved.
//

#ifndef __MRIDSS__MHD_BQlin__
#define __MRIDSS__MHD_BQlin__


#include "model.h"
#include "../Auxiliary/Initialization_routines.h"
#include "../Auxiliary/Input_parameters.h"

//



// Model class for S3T/CE2 shearing box MHD model
// Basic Quasi-linear MHD
// Derived from Model (model.h)
class MHD_BQlin : public Model {
public:
    MHD_BQlin(const Inputs& sp);
    ~MHD_BQlin();
    
    
    // Equations themselves
    void rhs(double t, double dt_lin,
            const dcmplxVec *MFin, const dcmplxMat *Cin,
            dcmplxVec *MFout, dcmplxMat *Cout,
            doubVec * linop_Ckl);
    // These are split into linear operator (diffusion) and nonlinear part
    // MF Diffusion operator separate - only calculate in constructor
    void linearOPs_Init(double t0,
                        doubVec * linop_MF, doubVec * linop_Ckl_old);
    
    // Dimension of the model
    int Cdimxy() const { return nxy_full_; };// x y dimension
    int Cdimz() const { return nz_Cfull_; };// z matrix dimension
    int MFdimz() const { return NZ_; };  // Size of MF vectors in z
    int num_MFs() const { return num_MF_; };  // Number of mean fields
    
    
    //  AUXILIARY FUNCTIONS
    void Calc_Energy(double* energy, double t,const dcmplxVec *MFin, const dcmplxMat *Cin ); // Energy
    void Calc_AngularMomentum(double* AM, double t,const dcmplxVec *MFin, const dcmplxMat *Cin ); // Angular momentum
    void Calc_Dissipation(double* AM, double t,const dcmplxVec *MFin, const dcmplxMat *Cin ); // Dissipation

private:
    
    // 2 Mean fields and 4 fluctuating fields in this model
    const int num_MF_ = 2;
    const int num_fluct_ = 4;
    int nz_Cfull_;
    
    // INPUT PARAMETERS
    const double nu_;
    const double eta_;// viscosity & resistivity
    const double q_;
    const double f_noise_; // driving noise

    // Temporary variables
    doubVec lapFtmp_, lap2tmp_;
    double kxtmp_,kytmp_;
    // Qkl temporary
    doubVec Qkl_tmp_;
    
};



#endif /* defined(__TwoDFluid__MHD_BQlin__) */
