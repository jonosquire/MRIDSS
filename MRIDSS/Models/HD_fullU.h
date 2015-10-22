//
//  HD_fullU.h
//  MRIDSS
//
//  Created by Jonathan Squire on 4/25/14.
//  Copyright (c) 2014 J Squire. All rights reserved.
//

#ifndef __MRIDSS__HD_fullU__
#define __MRIDSS__HD_fullU__


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
class HD_fullU : public Model {
public:
    HD_fullU(const Inputs& sp, MPIdata& mpi, fftwPlans& fft) ;
    ~HD_fullU();
    
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
    
    // Quasi-linear state
    void set_QL_YN(bool QL){QL_YN_ = QL;};
    
    // Dealiasing
    void dealias(dcmplxMat &inMat);
    void dealias(dcmplxVec& vec); // 1-D version
    void dealias(doubVec& vec); // 1-D version
    
    //  AUXILIARY FUNCTIONS
    //  Calculate energy, angular momentum and dissipation
    void Calc_Energy_AM_Diss(TimeVariables& tv, double t,const dcmplxVec *MFin, const dcmplxMat *Cin );
    int num_Reynolds_saves(){return MFdimz()*num_MFs();}; //Change if Calc_Energy_AM_Diss is modified!
    
    //////////////////////////
    // CFL number
    double Calculate_CFL(const dcmplxVec *MFin,const dcmplxMat* Ckl);
    double kmax;

private:
    
    
    // 2 Mean fields and 4 fluctuating fields in this model
    const int num_MF_;
    const int num_fluct_;
    int nz_Cfull_;
    
    // INPUT PARAMETERS
    const double nu_;
    const double q_;
    const double Omega_;
    const double f_noise_; // driving noise
    bool QL_YN_; // Turn on quasi-linear feedback
    
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
    // Create and store arrays of lap2
    void Define_Lap2_Arrays_(void);
    
    // Sizes for driving and energy
    long totalN2_;
    double mult_noise_fac_; // Factor to multiply noise to get values consistent with previous numbers
    double fft2Dfac_, fft1Dfac_, fftFac_Reynolds_;

    
    // turn off driving of ky=0, kx,kz != 0 modes (i.e., non-shearing waves)
    bool dont_drive_ky0_modes_Q_;
    // Noise properties
    double noise_range_[2]; // High/low cutoff
    Eigen::Array<bool,Eigen::Dynamic,1> drive_condition_; // Store driving in z as you loop  though x,y
    void print_noise_range_();
   
    
    ////////////////////////////////////////////////////
    //    TEMPORARY VARIABLES - SAME ACROSS MODELS    //
    doubVec lapFtmp_, lap2tmp_; // Laplacians - nice to still have lap2
    doubVec ilapFtmp_, ilap2tmp_; // Inverse Laplacians
    
    dcmplx kxctmp_, kyctmp_;
    double kxtmp_,kytmp_;
    // Qkl temporary
    doubVec Qkl_tmp_;
    
    
    // Real versions of MFs and derivatives
    dcmplxVecM Uy_, dzUy_, dzdzUy_;
    dcmplxVecM Ux_, dzUx_, dzdzUx_, dzdzdzUx_;
    
    // Reynolds stresses
    dcmplxMat reynolds_mat_tmp_; // Temporary matrix storage for fft
    
    // store both complex and double for fft and to pass less data around with MPI
    dcmplxVec ux_rey_stress_c_, uy_rey_stress_c_;// U reynolds stresses
    // double
    doubVec ux_rey_stress_d_, uy_rey_stress_d_;// U reynolds stresses
    doubVec reynolds_stress_MPI_send_, reynolds_stress_MPI_receive_; // mpi buffers
    
    double * reynolds_save_tmp_; // Saving the various contributions to mean-field dynamo


    // Automatically generated temporary variables - class definition
    dcmplxVecM rey_TdzTiLap2TkxbTky, rey_TiLap2TkxbTkyP2, rey_TiLap2TkxbP2Tky, rey_TiLap2TkxbTky, rey_TdzTiLap2Tkxb, rey_TdzTiLap2Tky, rey_TdzP2TiLap2, rey_TiLap2Tky, rey_TdzTiLap2, rey_Tdz;

    //////////////////////////////////////////////////////
    
    
    /////////////////////////////////////////////////
    //  AUTO GENERATED VARIABLES

    // Automatically generated temporary variables - class definition
    dcmplxVecM T2TdzP2TiLap2Tkxb, T2TdzTiLap2Tky, T2TdzTOm, TdzTiLap2Tkxb, TdzTkxbPLUSTmdzTiLap2TkxbP3, TdzTqPLUSTm2TdzTOm, TiLap2TkxbP2Tky, TiLap2Tky, Tm2TdzTiLap2TkxbP2Tky, Tm2TiLap2TkxbTkyP2, TmdzTiLap2Tkxb, TmiLap2TkxbP2Tky;
    
    dcmplx Tkxb, Tky, Tm2TkxbTkyTq, Tmkxb, Tmky;
    
    dcmplxMat Ctmp_1_, Ctmp_2_, Ctmp_3_, Ctmp_4_, Ctmp_5_, Ctmp_6_;
    
    dcmplxMat C11_, C12_, C21_, C22_;

    
    
    /////////////////////////////////////////////////

};



#endif /* defined__MRIDSS_HF_fullU__) */
