//
//  Constant_Damping.h
//  MRIDSS
//
//  Created by Jonathan Squire on 4/25/14.
//  Copyright (c) 2014 J Squire. All rights reserved.
//

#ifndef __MRIDSS__Constant_Damping__
#define __MRIDSS__Constant_Damping__


#include "Model.h"

#include "../Auxiliary/Initialization_routines.h"
#include "../Auxiliary/General_Auxiliary.h"
#include "../Auxiliary/Input_parameters.h"
#include "../Auxiliary/fftwPlans.h"
#include "../Auxiliary/TimeVariables.h"

//



// Model class for S3T/CE2 shearing box MHD model
// Basic Quasi-linear MHD
// Derived from Model (model.h)
class Constant_Damping : public Model {
public:
    Constant_Damping(const Inputs& sp, MPIdata& mpi, fftwPlans& fft) ;
    ~Constant_Damping();
    
    // Model name - for comparison wiht input file as a check
    const std::string equations_name;
    
    // Equations themselves
    void rhs(double t, double dt_lin,
             dcmplxVec *MFin, dcmplxMat *Cin,
             dcmplxVec *MFout, dcmplxMat *Cout,
             doubVec * linop_Ckl);
    // These are split into linear operator (diffusion) and nonlinear part
    // MF Diffusion operator separate - only calculate in constructor
    void linearOPs_Init(double t0,
                        doubVec * linop_MF, doubVec * linop_Ckl_old);
    
    // Dimension of the model
    int Cdimxy_full() const { return nxy_full_; };// x y dimension, all processes
    int Cdimz() const { return nz_Cfull_; };// z matrix dimension
    int MFdimz() const { return NZ_; };  // Size of MF vectors in z
    int num_MFs() const { return num_MF_; };  // Number of mean fields
    int num_linFs() const { return num_fluct_;}; // Number of fluctuating fields
    // MPI related
    int Cdimxy() const { return mpi_.nxy(); };// x y dimension
    int index_for_k_array() const { return mpi_.minxy_i(); }; // Index for each processor in k_ arrays
    // Box dimensions
    double box_length(int index) const { return L_[index];};
    
    
    //  AUXILIARY FUNCTIONS
    //  Calculate energy, angular momentum and dissipation
    void Calc_Energy_AM_Diss(TimeVariables& tv, double t,const dcmplxVec *MFin, const dcmplxMat *Cin );
    int num_Reynolds_saves(){return 2;};
    
private:
    
    // 2 Mean fields and 4 fluctuating fields in this model
    const int num_MF_;
    const int num_fluct_;
    int nz_Cfull_;
    
    // INPUT PARAMETERS
    const double nu_;
    const double eta_;// viscosity & resistivity
    const double q_;
    const double f_noise_; // driving noise
    const double QL_YN_; // Turn on quasi-linear feedback
    
    // MPI data
    MPIdata& mpi_; // Reference to MPI data
    fftwPlans& fft_;

    
    // Dealias setup
    int delaiasBnds_[2]; // Stores the bounds of the arrays for dealias
    
    // Useful arrays
    int totalN2_;
    double mult_noise_fac_; // Factor to multiply noise to get values consistent with previous numbers
    
    // Pre-calculate ilap2 related quantities to save computation (mainly for fft matrices)
    int* ky_index_; // Stores index in ky for a given full nxy index
    doubVec* lap2_, *ilap2_; // Laplacian and inverse

    void set_QL_YN(bool){mpi_.print1("DOES NOTHING!");};
    
    
    
    ////////////////////////////////////////////////////
    //               TEMPORARY VARIABLES              //
    doubVec lapFtmp_, lap2tmp_; // Laplacians - nice to still have lap2
    doubVec ilapFtmp_, ilap2tmp_; // Inverse Laplacians
    
    dcmplx kxctmp_, kyctmp_;
    double kxtmp_,kytmp_;
    // Qkl temporary
    doubVec Qkl_tmp_;
    
    doubVec Aop_tmp_; // A matrix, use asDiagonal
    
    
    
    //////////////////////////////////////////////////////
    //         PRIVATE FUNCTIONS FOR INITIALIZATION     //
    
    // Create and store arrays of lap2 to save computation of ffts
    void Define_Lap2_Arrays_(void);
    
};



#endif /* defined(__TwoDFluid__Constant_Damping__) */
