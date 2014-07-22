//
//  MHD_BQlin.h
//  MRIDSS
//
//  Created by Jonathan Squire on 4/25/14.
//  Copyright (c) 2014 J Squire. All rights reserved.
//

#ifndef __MRIDSS__MHD_BQlin__
#define __MRIDSS__MHD_BQlin__


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
class MHD_BQlin : public Model {
public:
    MHD_BQlin(const Inputs& sp, MPIdata& mpi, fftwPlans& fft) ;
    ~MHD_BQlin();
    
    
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
    // Box dimensions
    double box_length(int index) const { return L_[index];};
    
    
    // Dealiasing
    void dealias(dcmplx *arr); // 2-D version - TODO tidy up
    void dealias(dcmplxVec& vec); // 1-D version
    void dealias(doubVec& vec);
    
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
    
    // Sizes for driving and energy
    long totalN2_;
    double mult_noise_fac_; // Factor to multiply noise to get values consistent with previous numbers
    
    // turn off driving of ky=0, kx,kz != 0 modes (i.e., non-shearing waves)
    bool dont_drive_ky0_modes_Q_;
    
    // Reynolds stress - store both complex and double for fft and to pass less data around with MPI
    // complex
    dcmplxVec bzux_m_uzbx_c_; // bz*ux-uz*bx
    dcmplxVec bzuy_m_uzby_c_; // bz*uy - uz*by
    // double
    doubVec bzux_m_uzbx_d_, bzux_m_uzbx_drec_; // bz*ux-uz*bx
    doubVec bzuy_m_uzby_d_, bzuy_m_uzby_drec_; // bz*uy - uz*by
    
    double * reynolds_save_tmp_; // Saving the various contributions to mean-field dynamo

    
    ////////////////////////////////////////////////////
    //               TEMPORARY VARIABLES              //
    doubVec lapFtmp_, lap2tmp_; // Laplacians - nice to still have lap2
    doubVec ilapFtmp_, ilap2tmp_; // Inverse Laplacians
    
    dcmplx kxctmp_, kyctmp_;
    double kxtmp_,kytmp_;
    // Qkl temporary
    doubVec Qkl_tmp_;
    // Operator matrices - put into submatrices rather than storing full matrix
    // Some are zero and some are diagonal (eigen .asDiagonal)
    dcmplxVecM Aop_v11_, Aop_v12_, Aop_v21_, Aop_v43_;
    dcmplxMat Aop_m13_, Aop_m14_, Aop_m23_, Aop_m24_, Aop_m31_, Aop_m41_, Aop_m42_;
    dcmplxMat Aop_block_tmp_;

    // Ckl_in submatrices - reset after each column for memory
    dcmplxMat C1_,C2_,C3_,C4_;
    
    // Real versions of B and derivatives
    dcmplxVec rBy_tmp_;
    dcmplxVec rDzBy_tmp_;
    dcmplxVec rDzzBy_tmp_;
    // Reynolds stresses
    dcmplxMat reynolds_mat_tmp_; // Temporary matrix storage for fft
    dcmplxVecM rey_mkxky_tmp_,rey_kz_tmp_,rey_mkxkz_tmp_,rey_mky_tmp_; // Convenient to store vectors for converting between u, zeta etc. and u uy uz...

    //////////////////////////////////////////////////////

    //    BLOCK MATRIX MULTIPLICATION                   //
    void Block_Matrix_Mult_(int column, dcmplxMat& Ckl_in_i,dcmplxMat& Ckl_out_i);
    
    
    //////////////////////////////////////////////////////
    //         PRIVATE FUNCTIONS FOR INITIALIZATION     //
    
    // Create and store arrays of lap2 to save computation of ffts
    void Define_Lap2_Arrays_(void);
    // fft of the identity matrix
    void Set_fft_identity_(void);
    
};



#endif /* defined(__TwoDFluid__MHD_BQlin__) */
