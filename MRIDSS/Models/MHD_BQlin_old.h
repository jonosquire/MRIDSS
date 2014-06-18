//
//  MHD_BQlin_old.h
//  MRIDSS
//
//  Created by Jonathan Squire on 4/25/14.
//  Copyright (c) 2014 J Squire. All rights reserved.
//

#ifndef __MRIDSS__MHD_BQlin_old__
#define __MRIDSS__MHD_BQlin_old__


#include "Model.h"

#include "../Auxiliary/Initialization_routines.h"
#include "../Auxiliary/General_Auxiliary.h"
#include "../Auxiliary/Input_parameters.h"

//



// Model class for S3T/CE2 shearing box MHD model
// Basic Quasi-linear MHD
// Derived from Model (model.h)
class MHD_BQlin_old : public Model {
public:
    MHD_BQlin_old(const Inputs& sp, MPIdata& mpi, fftwPlans& fft) ;
    ~MHD_BQlin_old();
    
    
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
    // MPI related
    int Cdimxy() const { return mpi_.nxy(); };// x y dimension
    int index_for_k_array() const { return mpi_.minxy_i(); }; // Index for each processor in k_ arrays
    
    
    // Dealiasing
    void dealias(dcmplx *arr); // 2-D version - TODO tidy up
    void dealias(dcmplxVec& vec); // 1-D version
    
    //  AUXILIARY FUNCTIONS
    //  Calculate energy, angular momentum and dissipation
    void Calc_Energy_AM_Diss(TimeVariables& tv, double t,const dcmplxVec *MFin, const dcmplxMat *Cin );

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
    
    // FFTW plans
    fftwPlans& fft_;
    
    // Dealias setup
    int delaiasBnds_[2]; // Stores the bounds of the arrays for dealias
    
    // Useful arrays to save computation
    dcmplxMat fft_identity_; // ifft of the identity matrix
    // Pre-calculate ilap2 related quantities to save computation (mainly for fft matrices)
    int* ky_index_; // Stores index in ky for a given full nxy index
    doubVec* lap2_, *ilap2_; // Laplacian and inverse
    dcmplxMat* fft_ilap2_,*fft_kzilap2_,*fft_kz2ilap2_; // ifft of ilap2, kz*ilap2 and kz^2*ilap2
    
    
    // Reynolds stress - store both complex and double for fft and to pass less data around with MPI
    // complex
    dcmplxVec bzux_m_uzbx_c_; // bz*ux-uz*bx
    dcmplxVec bzuy_m_uzby_c_; // bz*uy - uz*by
    // double
    doubVec bzux_m_uzbx_d_, bzux_m_uzbx_drec_; // bz*ux-uz*bx
    doubVec bzuy_m_uzby_d_, bzuy_m_uzby_drec_; // bz*uy - uz*by

    
    ////////////////////////////////////////////////////
    //               TEMPORARY VARIABLES              //
    doubVec lapFtmp_, lap2tmp_; // Laplacians - nice to still have lap2
    doubVec ilapFtmp_, ilap2tmp_; // Inverse Laplacians
    
    dcmplx kxctmp_, kyctmp_;
    double kxtmp_,kytmp_;
    // Qkl temporary
    doubVec Qkl_tmp_;
    // Operator matrix
    dcmplxMat Aop_tmp_;
    dcmplxMat Aop_block_tmp_,Aop_block_tmp2_; // Blocks to use as temporaries
    // Real versions of B and derivatives
    dcmplxVec rBy_tmp_;
    dcmplxVec rDzBy_tmp_;
    dcmplxVec rDzzBy_tmp_;
    // Reynolds stresses
    dcmplxMat reynolds_mat_tmp_; // Temporary matrix storage for fft
    Eigen::Matrix<dcmplx, Eigen::Dynamic, 1> rey_mkxky_tmp_,rey_kz_tmp_,rey_mkxkz_tmp_,rey_mky_tmp_; // Convenient to store vectors for converting between u, zeta etc. and u uy uz...

    
    
    //////////////////////////////////////////////////////
    //         PRIVATE FUNCTIONS FOR INITIALIZATION     //
    
    // Create and store arrays of lap2 to save computation of ffts
    void Define_Lap2_Arrays_(void);
    // fft of the identity matrix
    void Set_fft_identity_(void);
    
};



#endif /* defined(__TwoDFluid__MHD_BQlin_old__) */