//
//  Model.h
//  TwoDFluid
//
//  Created by Jonathan Squire on 4/25/14.
//  Copyright (c) 2014 J Squire. All rights reserved.
//

#ifndef MRIDSS_Model_h
#define MRIDSS_Model_h

#include "../General_Definitions.h"
#include "../Auxiliary/MPIdata.h"
#include "../Auxiliary/TimeVariables.h"

// Abstract model class - container for equations of motion for MRIDSS
// Ckl data is a pointer array (size dimxy) of Eigen matrices (size dimz^2)
// Mean field data is a (num_MF) pointer array of Eigen matrices of size NZ


class Model {
public:
    
    Model(int NZ,const int NXY[2] ,const double * L) :
        L_(L) {
        // These dimensions do not depend on model (i.e., number of variables, MFs etc.)
        NZ_ = NZ;
        Nxy_[0] = NXY[0]; // Nyquist (and zero) frequncy included in NXY[0]
        Nxy_[1]= NXY[1]/2; // ONLY USE HALF OF THE TRANSFORM FOR THE Y DIMENSION DUE TO FULL Ckl BEING REAL!
        nxy_full_= (Nxy_[0]-1)*Nxy_[1] ;// Leave out Nyquist frequency in kx
        // NB: Currently evolving the (0,0) component of Ckl. This is uneccessary but more convenient (especially with parallelization), will just stay zero anyway.
    }
    virtual ~Model() {};
    
    // fx = f(x,t)
    // Returns nonlinear part of of RHS and a linear operator (linop_Ckl, linop_MF), so
    //  as to allow for semi-implicit integrators.
    // The linear operators are evaluated at t + dt_lin
    virtual void rhs(double t, double dt_lin,
                    dcmplxVec *MFin, dcmplxMat *Cin,
                    dcmplxVec *MFout, dcmplxMat *Cout,
                    doubVec * linop_Ckl) = 0;
    
    // Initialization of linear operators - for mean fields linop is constant anyway
    virtual void linearOPs_Init(double t0, doubVec * linop_MF, doubVec * linop_Ckl_old) = 0;
    
    // number of states (size of x)
    virtual int Cdimxy_full() const = 0; // Full C size in x,y
    virtual int Cdimz() const = 0;  // In z (size of eigen matrix)
    virtual int MFdimz() const = 0;  // Size of MF vectors in z
    virtual int num_MFs() const = 0;  // Number of mean fields
    // Box size
    virtual double box_length(int index) const =0; // Size of box
    
    // MPI related
    virtual int Cdimxy() const = 0; // Size of C array on given MPI process
    virtual int index_for_k_array() const =0;
    
    // Return kx - necessary for remapping
    dcmplx* kx_pointer() const {return kx_;};
    dcmplx* ky_pointer() const {return ky_;};
    
    // Remapping procedure
    friend void ShearingBox_Remap(Model* model,dcmplxMat* Ckl);
    
    //////////////////////////////////////////////////////////////////
    //////  AUXILIARY FUNCTIONS OPERATING ON SOLUTION   //////////////
    virtual void Calc_Energy_AM_Diss(TimeVariables& tv, double t,const dcmplxVec *MFin, const dcmplxMat *Cin ) = 0;
    //////////////////////////////////////////////////////////////////

protected:
    
    // DIMENSIONS - these do not depend on model
    int Nxy_[2];           // Dimensions (Nx, Ny/2)
    int NZ_;
    int nxy_full_;      // Ckl Nxy[0]*Nxy[1] total dimension
    const double * L_;       // Box size (Lx,Ly,Lz)
    
    // K VECTORS DATA
    dcmplx *kx_,*ky_;
    dcmplxVec kz_;
    doubVec kz2_; // kz^2
    
};



#endif  // MODEL_H_
