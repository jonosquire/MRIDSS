//
//  Model_AutoGen_template.cc
//  MRIDSS
//
//  Created by Jonathan Squire on 4/25/14.
//  Copyright (c) 2014 J Squire. All rights reserved.
//

#include "Model_AutoGen_template.h"

// Model class for S3T/CE2 shearing box MHD model
// Template (not in c++ sense) for the equations to be copied into from automatically generated Mathematica file
//
// Derived from Model (model.h)


///   DOES NOT WORK RIGHT NOW!!!!!!!!


// Constructor for Model_AutoGen_template
Model_AutoGen_template::Model_AutoGen_template(const Inputs& sp, MPIdata& mpi, fftwPlans& fft) :
num_MF_(2), num_fluct_(4),
q_(sp.q), f_noise_(sp.f_noise), nu_(sp.nu), eta_(sp.eta),
QL_YN_(sp.QuasiLinearQ),
Model(sp.NZ, sp.NXY , sp.L), // Dimensions - stored in base
mpi_(mpi), // MPI data
fft_(fft) // FFT data
{
    // Dimensions that depend on model
    nz_Cfull_=num_fluct_*NZ_;
    // Calcuate kx and ky arrays - single dimesnsion NX*NY(/2) to simplify parallelization
    kx_ = new dcmplx[ nxy_full_ ]; // Entire kx_ and ky_ arrays are stored on all processors
    ky_ = new dcmplx[ nxy_full_ ];
    ky_index_ = new int[ nxy_full_ ];
    define_kxy_array(kx_, ky_,ky_index_, Nxy_, L_);
    
    // kz array - Eigen vector
    kz_ = dcmplxVec( NZ_ );
    kz2_ = doubVec( NZ_ );
    define_kz_array(kz_, NZ_, L_);
    kz2_ = (kz_*kz_).real(); // kz2_ is a double Eigen array (rather than dcmplx). But keeps its negative sign! i.e., use in the same way as kz*kz

    
    
    // Setup MPI
    mpi_.Split_NXY_Grid( nxy_full_ ); // Work out MPI splitting
    
    // Set fftw plans
    fft_.calculatePlans( NZ_ );
    
    // Delalias setup
    delaiasBnds_[0] = NZ_/3+1; // Delete array including these elements
    delaiasBnds_[1] = NZ_-NZ_/3-1; // PUT BACK TO NZ_!!!!

    // Sizes for
    totalN2_ = Nxy_[0]*Nxy_[1]*2*NZ_;
    totalN2_ = totalN2_*totalN2_;
    mult_noise_fac_ = 1.0/(16*32*32); // Defined so that consistent with (16, 32, 32) results
    mult_noise_fac_ = mult_noise_fac_*mult_noise_fac_;
    // Factor for inverse transform in 2-D 1/NZ^2
    fft2Dfac_=1.0/NZ_/NZ_;
    fft1Dfac_ = 1.0/NZ_;
    
    // turn off driving of ky=0, kx,kz != 0 modes (i.e., non-shearing waves)
    dont_drive_ky0_modes_Q_ = 0;
    
    ////////////////////////////////////////////////////
    // Useful arrays to save computation
    Define_Lap2_Arrays_(); // Lap2- Allocates data to lap2_,ilap2_
    
    // Reynolds stresses
    bzux_m_uzbx_c_ = dcmplxVec::Zero(NZ_);
    bzuy_m_uzby_c_ = dcmplxVec::Zero(NZ_);
    bzux_m_uzbx_d_ = doubVec::Zero(NZ_);
    bzuy_m_uzby_d_ = doubVec::Zero(NZ_);
    // Receive buffers for MPI, for some reason MPI::IN_PLACE is not giving the correct result!
    bzux_m_uzbx_drec_ = doubVec::Zero(NZ_);
    bzuy_m_uzby_drec_ = doubVec::Zero(NZ_);
    
    reynolds_save_tmp_ = new double[5];
    
    
    ////////////////////////////////////////////////////
    //                                                //
    //               TEMPORARY ARRAYS                 //
    //                                                //
    //  These are used during evaluation of the various parts of the model
    // There should never be any need to keep their value over more than 1 time-step
    
    //////////////////////////////////////////////////
    // These arrays are used for all sets of equations
    lapFtmp_ = doubVec( NZ_ ); //For (time-dependent) k^2
    ilapFtmp_ = doubVec( NZ_ ); // Inverse
    lap2tmp_ = doubVec( NZ_ ); // For ky^2+kz^2 - could be pre-assigned
    ilap2tmp_ = doubVec( NZ_ ); // Inverse
    // Qkl
    Qkl_tmp_ = doubVec( nz_Cfull_ );
    //////////////////////////////////////////////////


    // Real versions of B
    By_ = dcmplxVecM::Zero(NZ_);
    dzBy_ = dcmplxVecM::Zero(NZ_);
    dzdzBy_ = dcmplxVecM::Zero(NZ_);
    
    // Reynolds stresses - this need to be changed to automatic also
    reynolds_mat_tmp_ = dcmplxMat(NZ_,NZ_);
    // converting between u, zeta... and u, uy... for the Reynolds stresses
    rey_mkxky_tmp_ = dcmplxVecM( NZ_ );
    rey_kz_tmp_ = dcmplxVecM( NZ_ );
    rey_mkxkz_tmp_ = dcmplxVecM( NZ_ );
    rey_mky_tmp_ = dcmplxVecM( NZ_ );
    
    
    /////////////////////////////////////////////////
    //  AUTO GENERATED VARIABLES
    
    // Automatically generated temporary variables - class constructor
    T2Tdz = dcmplxVecM(NZ_);
    T2TdzTiLap2TkxbP2Tky = dcmplxVecM(NZ_);
    T2TdzTiLap2Tky = dcmplxVecM(NZ_);
    T2TiLap2TkxbTkyP2 = dcmplxVecM(NZ_);
    TdzP2TiLap2TkxbPLUSTmiLap2TkxbTkyP2 = dcmplxVecM(NZ_);
    TdzTiLap2Tkxb = dcmplxVecM(NZ_);
    TdzTqPLUSTm2Tdz = dcmplxVecM(NZ_);
    TiLap2Tky = dcmplxVecM(NZ_);
    TmdzTiLap2Tkxb = dcmplxVecM(NZ_);
    TmdzTq = dcmplxVecM(NZ_);
    TmiLap2Tky = dcmplxVecM(NZ_);
    
    Ctmp_1_ = dcmplxMat(NZ_,NZ_);
    Ctmp_2_ = dcmplxMat(NZ_,NZ_);
    Ctmp_3_ = dcmplxMat(NZ_,NZ_);
    
    C11_ = dcmplxMat(NZ_,NZ_);
    C12_ = dcmplxMat(NZ_,NZ_);
    C13_ = dcmplxMat(NZ_,NZ_);
    C14_ = dcmplxMat(NZ_,NZ_);
    C21_ = dcmplxMat(NZ_,NZ_);
    C22_ = dcmplxMat(NZ_,NZ_);
    C23_ = dcmplxMat(NZ_,NZ_);
    C24_ = dcmplxMat(NZ_,NZ_);
    C31_ = dcmplxMat(NZ_,NZ_);
    C32_ = dcmplxMat(NZ_,NZ_);
    C33_ = dcmplxMat(NZ_,NZ_);
    C34_ = dcmplxMat(NZ_,NZ_);
    C41_ = dcmplxMat(NZ_,NZ_);
    C42_ = dcmplxMat(NZ_,NZ_);
    C43_ = dcmplxMat(NZ_,NZ_);
    C44_ = dcmplxMat(NZ_,NZ_);
    /////////////////////////////////////////////////


}


// Destructor
Model_AutoGen_template::~Model_AutoGen_template(){
    // Here I am destroying Base members from the derived destructor - is this bad?
    delete[] kx_;
    delete[] ky_;
    
    delete[] lap2_;
    delete[] ilap2_;
    // Data arrays
    delete[] reynolds_save_tmp_;
}


// Main EulerFluid function, fx = f(t,x) called at ecah time-step
void Model_AutoGen_template::rhs(double t, double dt_lin,
                   dcmplxVec *MFin, dcmplxMat *Ckl_in,
                   dcmplxVec *MFout, dcmplxMat *Ckl_out,
                   doubVec * linop_Ckl) {
    
    // If desired, the function should work fine with Ckl_in = Ckl_out. This is
    // useful for integration schemes that do not require re-use of Ckl_in
    
    // Calculate MFs in real space - would be good to also include this in auto generation
    By_ = MFin[1].matrix()*fft1Dfac_; // fft back doesn't include normalization
    fft_.back_1D(By_.data());
    dzBy_ = (kz_*MFin[1]).matrix()*fft1Dfac_;
    fft_.back_1D(dzBy_.data());
    dzdzBy_ = (kz2_*MFin[1]).matrix()*fft1Dfac_;
    fft_.back_1D(dzdzBy_.data());
    
    // Reynolds stresses -- added to at each step in the loop
    bzux_m_uzbx_d_.setZero();
    bzuy_m_uzby_d_.setZero();
    
    
    /////////////////////////////////////
    //   ALL OF THIS IS PARALLELIZED
    for (int i=0; i<Cdimxy();  ++i){
        
        // Full Loop containing all the main work on Ckl
        
        ///////////////////////////////////////
        ///// MAIN CKL EQUATIONS
        
        int k_i = i + index_for_k_array(); // k index
        int ind_ky = ky_index_[k_i];
        
        // Form Laplacians using time-dependent kx
        kyctmp_ = ky_[k_i];
        kytmp_ = kyctmp_.imag(); // It is useful to have both complex and double versions
        kxctmp_ = kx_[k_i] + q_*t*kyctmp_;
        kxtmp_ = kxctmp_.imag();
        
        lap2tmp_ = lap2_[ind_ky];
        lapFtmp_ = -kxtmp_*kxtmp_ + lap2tmp_;
        ilapFtmp_ = 1/lapFtmp_;
        ilap2tmp_ = ilap2_[ind_ky];;
        
        ////////////////////////////////////////
        //              FORM Qkl              //
        
        bool Qkl_Zero_Q = kytmp_==0 && kxtmp_==0; // Whether Qkl  should be zero for this mode
        if (dont_drive_ky0_modes_Q_ && kytmp_==0)
            Qkl_Zero_Q = 1;  // Option in the code, set to zero for all ky=0 modes if set
        
        if (Qkl_Zero_Q) {
            Qkl_tmp_.setZero(); // Don't drive the mean components! (really could leave mean cmpts out entirely)
            ilapFtmp_(0) = 1; // Avoid infinities (lap2 already done, since stored)
        }
        else {
            lapFtmp_ = lap2tmp_/lapFtmp_; // lapFtmp_ is no longer lapF!!!
            dealias(lapFtmp_);
            dealias(lap2tmp_);
            // CANNOT USE AFTER THIS!
            if (kytmp_==0.0) { // Don't drive the ky=kz=0 component (variables not invertible)
                lapFtmp_(0)=0;
            }
            
            Qkl_tmp_ << lapFtmp_, lap2tmp_, lapFtmp_, lap2tmp_;
            Qkl_tmp_ = (f_noise_*f_noise_*totalN2_*mult_noise_fac_)*Qkl_tmp_.abs();

            
        }
        ////////////////////////////////////////
        
        
        /////////////////////////////////////////
        ///    INSERT VARIABLE DEFS HERE
        
        
        //Assign automatically generated variables in equations
        // Submatrix variable definition (automatic)
        C11_ = Ckl_in[i].block( 0, 0, NZ_, NZ_);
        C12_ = Ckl_in[i].block( 0, NZ_, NZ_, NZ_);
        C13_ = Ckl_in[i].block( 0, 2*NZ_, NZ_, NZ_);
        C14_ = Ckl_in[i].block( 0, 3*NZ_, NZ_, NZ_);
        C21_ = Ckl_in[i].block( NZ_, 0, NZ_, NZ_);
        C22_ = Ckl_in[i].block( NZ_, NZ_, NZ_, NZ_);
        C23_ = Ckl_in[i].block( NZ_, 2*NZ_, NZ_, NZ_);
        C24_ = Ckl_in[i].block( NZ_, 3*NZ_, NZ_, NZ_);
        C31_ = Ckl_in[i].block( 2*NZ_, 0, NZ_, NZ_);
        C32_ = Ckl_in[i].block( 2*NZ_, NZ_, NZ_, NZ_);
        C33_ = Ckl_in[i].block( 2*NZ_, 2*NZ_, NZ_, NZ_);
        C34_ = Ckl_in[i].block( 2*NZ_, 3*NZ_, NZ_, NZ_);
        C41_ = Ckl_in[i].block( 3*NZ_, 0, NZ_, NZ_);
        C42_ = Ckl_in[i].block( 3*NZ_, NZ_, NZ_, NZ_);
        C43_ = Ckl_in[i].block( 3*NZ_, 2*NZ_, NZ_, NZ_);
        C44_ = Ckl_in[i].block( 3*NZ_, 3*NZ_, NZ_, NZ_);
        
        // Scalar variable definition (automatic)
        Tky = kyctmp_;
        Tm2TkxbTkyTq = (-2.*kxctmp_*kyctmp_)*q_;
        Tmkxb = (-1.*kxctmp_);
        
        // Vector variable definition (automatic)
        T2Tdz = 2.*kz_.matrix();
        T2TdzTiLap2TkxbP2Tky = (2.*kxctmp_*kxctmp_*kyctmp_)*(ilap2tmp_*kz_).matrix();
        T2TdzTiLap2Tky = (2.*kyctmp_)*(ilap2tmp_*kz_).matrix();
        T2TiLap2TkxbTkyP2 = (2.*kxctmp_*kyctmp_*kyctmp_)*ilap2tmp_.matrix();
        TdzP2TiLap2TkxbPLUSTmiLap2TkxbTkyP2 = (-1.*kxctmp_*kyctmp_*kyctmp_)*ilap2tmp_.matrix() + (ilap2tmp_*kz_*kz_).matrix()*kxctmp_;
        TdzTiLap2Tkxb = (ilap2tmp_*kz_).matrix()*kxctmp_;
        TdzTqPLUSTm2Tdz = -2.*kz_.matrix() + kz_.matrix()*q_;
        TiLap2Tky = ilap2tmp_.matrix()*kyctmp_;
        TmdzTiLap2Tkxb = (-1.*kxctmp_)*(ilap2tmp_*kz_).matrix();
        TmdzTq = (-1.*q_)*kz_.matrix();
        TmiLap2Tky = (-1.*kyctmp_)*ilap2tmp_.matrix();
        
        
        /////////////////////////////////////////
        ///    INSERT EQUATIONS HERE
//        std::cout << "kx = " << kxtmp_ << ", ky = " << kytmp_ << std::endl;
//        std::cout << Ctmp_4_ << std::endl << std::endl;
        
        ////////////////////////////////////////////////////////
        ///       AUTOMATICALLY GENERATED EQUATIONS         /////
        ///       see GenerateC++Equations.nb in MMA        /////
        
        Ctmp_1_ = (fft2Dfac_*Tky)*C31_;
        fft_.back_2DFull( Ctmp_1_ );
        Ctmp_1_ = By_.asDiagonal()*Ctmp_1_;
        Ctmp_2_ = (fft2Dfac_*T2TdzTiLap2TkxbP2Tky).asDiagonal()*C31_ + (fft2Dfac_*T2TiLap2TkxbTkyP2).asDiagonal()*C41_;
        fft_.back_2DFull( Ctmp_2_ );
        Ctmp_2_ = dzBy_.asDiagonal()*Ctmp_2_;
        fft_.for_2DFull( Ctmp_1_ );
        fft_.for_2DFull( Ctmp_2_ );
        Ctmp_2_ = Ctmp_1_ + ilapFtmp_.matrix().asDiagonal()*(Ctmp_2_ + T2Tdz.asDiagonal()*C21_ + Tm2TkxbTkyTq*C11_);
        dealias(Ctmp_2_);
        Ckl_out[i].block( 0, 0, NZ_, NZ_) = Ctmp_2_;
        
        
        
        Ctmp_1_ = (fft2Dfac_*Tky)*C32_;
        fft_.back_2DFull( Ctmp_1_ );
        Ctmp_1_ = By_.asDiagonal()*Ctmp_1_;
        Ctmp_2_ = (fft2Dfac_*T2TdzTiLap2TkxbP2Tky).asDiagonal()*C32_ + (fft2Dfac_*T2TiLap2TkxbTkyP2).asDiagonal()*C42_;
        fft_.back_2DFull( Ctmp_2_ );
        Ctmp_2_ = dzBy_.asDiagonal()*Ctmp_2_;
        fft_.for_2DFull( Ctmp_1_ );
        fft_.for_2DFull( Ctmp_2_ );
        Ctmp_2_ = Ctmp_1_ + ilapFtmp_.matrix().asDiagonal()*(Ctmp_2_ + T2Tdz.asDiagonal()*C22_ + Tm2TkxbTkyTq*C12_);
        dealias(Ctmp_2_);
        Ckl_out[i].block( 0, NZ_, NZ_, NZ_) = Ctmp_2_;
        
        
        
        Ctmp_1_ = (fft2Dfac_*Tky)*C33_;
        fft_.back_2DFull( Ctmp_1_ );
        Ctmp_1_ = By_.asDiagonal()*Ctmp_1_;
        Ctmp_2_ = (fft2Dfac_*T2TdzTiLap2TkxbP2Tky).asDiagonal()*C33_ + (fft2Dfac_*T2TiLap2TkxbTkyP2).asDiagonal()*C43_;
        fft_.back_2DFull( Ctmp_2_ );
        Ctmp_2_ = dzBy_.asDiagonal()*Ctmp_2_;
        fft_.for_2DFull( Ctmp_1_ );
        fft_.for_2DFull( Ctmp_2_ );
        Ctmp_2_ = Ctmp_1_ + ilapFtmp_.matrix().asDiagonal()*(Ctmp_2_ + T2Tdz.asDiagonal()*C23_ + Tm2TkxbTkyTq*C13_);
        dealias(Ctmp_2_);
        Ckl_out[i].block( 0, 2*NZ_, NZ_, NZ_) = Ctmp_2_;
        
        
        
        Ctmp_1_ = (fft2Dfac_*Tky)*C34_;
        fft_.back_2DFull( Ctmp_1_ );
        Ctmp_1_ = By_.asDiagonal()*Ctmp_1_;
        Ctmp_2_ = (fft2Dfac_*T2TdzTiLap2TkxbP2Tky).asDiagonal()*C34_ + (fft2Dfac_*T2TiLap2TkxbTkyP2).asDiagonal()*C44_;
        fft_.back_2DFull( Ctmp_2_ );
        Ctmp_2_ = dzBy_.asDiagonal()*Ctmp_2_;
        fft_.for_2DFull( Ctmp_1_ );
        fft_.for_2DFull( Ctmp_2_ );
        Ctmp_2_ = Ctmp_1_ + ilapFtmp_.matrix().asDiagonal()*(Ctmp_2_ + T2Tdz.asDiagonal()*C24_ + Tm2TkxbTkyTq*C14_);
        dealias(Ctmp_2_);
        Ckl_out[i].block( 0, 3*NZ_, NZ_, NZ_) = Ctmp_2_;
        
        
        
        Ctmp_1_ = (fft2Dfac_*Tky)*C41_;
        fft_.back_2DFull( Ctmp_1_ );
        Ctmp_1_ = By_.asDiagonal()*Ctmp_1_;
        Ctmp_2_ = (fft2Dfac_*TmdzTiLap2Tkxb).asDiagonal()*C31_ + (fft2Dfac_*TmiLap2Tky).asDiagonal()*C41_;
        fft_.back_2DFull( Ctmp_2_ );
        Ctmp_2_ = dzdzBy_.asDiagonal()*Ctmp_2_;
        Ctmp_3_ = (fft2Dfac_*Tmkxb)*C31_;
        fft_.back_2DFull( Ctmp_3_ );
        Ctmp_3_ = dzBy_.asDiagonal()*Ctmp_3_;
        Ctmp_1_ = Ctmp_1_ + Ctmp_2_ + Ctmp_3_;
        fft_.for_2DFull( Ctmp_1_ );
        Ctmp_3_ = Ctmp_1_ + TdzTqPLUSTm2Tdz.asDiagonal()*C11_;
        dealias(Ctmp_3_);
        Ckl_out[i].block( NZ_, 0, NZ_, NZ_) = Ctmp_3_;
        
        
        
        Ctmp_1_ = (fft2Dfac_*Tky)*C42_;
        fft_.back_2DFull( Ctmp_1_ );
        Ctmp_1_ = By_.asDiagonal()*Ctmp_1_;
        Ctmp_2_ = (fft2Dfac_*TmdzTiLap2Tkxb).asDiagonal()*C32_ + (fft2Dfac_*TmiLap2Tky).asDiagonal()*C42_;
        fft_.back_2DFull( Ctmp_2_ );
        Ctmp_2_ = dzdzBy_.asDiagonal()*Ctmp_2_;
        Ctmp_3_ = (fft2Dfac_*Tmkxb)*C32_;
        fft_.back_2DFull( Ctmp_3_ );
        Ctmp_3_ = dzBy_.asDiagonal()*Ctmp_3_;
        Ctmp_1_ = Ctmp_1_ + Ctmp_2_ + Ctmp_3_;
        fft_.for_2DFull( Ctmp_1_ );
        Ctmp_3_ = Ctmp_1_ + TdzTqPLUSTm2Tdz.asDiagonal()*C12_;
        dealias(Ctmp_3_);
        Ckl_out[i].block( NZ_, NZ_, NZ_, NZ_) = Ctmp_3_;
        
        
        
        Ctmp_1_ = (fft2Dfac_*Tky)*C43_;
        fft_.back_2DFull( Ctmp_1_ );
        Ctmp_1_ = By_.asDiagonal()*Ctmp_1_;
        Ctmp_2_ = (fft2Dfac_*TmdzTiLap2Tkxb).asDiagonal()*C33_ + (fft2Dfac_*TmiLap2Tky).asDiagonal()*C43_;
        fft_.back_2DFull( Ctmp_2_ );
        Ctmp_2_ = dzdzBy_.asDiagonal()*Ctmp_2_;
        Ctmp_3_ = (fft2Dfac_*Tmkxb)*C33_;
        fft_.back_2DFull( Ctmp_3_ );
        Ctmp_3_ = dzBy_.asDiagonal()*Ctmp_3_;
        Ctmp_1_ = Ctmp_1_ + Ctmp_2_ + Ctmp_3_;
        fft_.for_2DFull( Ctmp_1_ );
        Ctmp_3_ = Ctmp_1_ + TdzTqPLUSTm2Tdz.asDiagonal()*C13_;
        dealias(Ctmp_3_);
        Ckl_out[i].block( NZ_, 2*NZ_, NZ_, NZ_) = Ctmp_3_;
        
        
        
        Ctmp_1_ = (fft2Dfac_*Tky)*C44_;
        fft_.back_2DFull( Ctmp_1_ );
        Ctmp_1_ = By_.asDiagonal()*Ctmp_1_;
        Ctmp_2_ = (fft2Dfac_*TmdzTiLap2Tkxb).asDiagonal()*C34_ + (fft2Dfac_*TmiLap2Tky).asDiagonal()*C44_;
        fft_.back_2DFull( Ctmp_2_ );
        Ctmp_2_ = dzdzBy_.asDiagonal()*Ctmp_2_;
        Ctmp_3_ = (fft2Dfac_*Tmkxb)*C34_;
        fft_.back_2DFull( Ctmp_3_ );
        Ctmp_3_ = dzBy_.asDiagonal()*Ctmp_3_;
        Ctmp_1_ = Ctmp_1_ + Ctmp_2_ + Ctmp_3_;
        fft_.for_2DFull( Ctmp_1_ );
        Ctmp_3_ = Ctmp_1_ + TdzTqPLUSTm2Tdz.asDiagonal()*C14_;
        dealias(Ctmp_3_);
        Ckl_out[i].block( NZ_, 3*NZ_, NZ_, NZ_) = Ctmp_3_;
        
        
        
        Ctmp_1_ = (fft2Dfac_*Tky)*C11_;
        fft_.back_2DFull( Ctmp_1_ );
        Ctmp_1_ = By_.asDiagonal()*Ctmp_1_;
        fft_.for_2DFull( Ctmp_1_ );
        dealias(Ctmp_1_);
        Ckl_out[i].block( 2*NZ_, 0, NZ_, NZ_) = Ctmp_1_;
        
        
        
        Ctmp_1_ = (fft2Dfac_*Tky)*C12_;
        fft_.back_2DFull( Ctmp_1_ );
        Ctmp_1_ = By_.asDiagonal()*Ctmp_1_;
        fft_.for_2DFull( Ctmp_1_ );
        dealias(Ctmp_1_);
        Ckl_out[i].block( 2*NZ_, NZ_, NZ_, NZ_) = Ctmp_1_;
        
        
        
        Ctmp_1_ = (fft2Dfac_*Tky)*C13_;
        fft_.back_2DFull( Ctmp_1_ );
        Ctmp_1_ = By_.asDiagonal()*Ctmp_1_;
        fft_.for_2DFull( Ctmp_1_ );
        dealias(Ctmp_1_);
        Ckl_out[i].block( 2*NZ_, 2*NZ_, NZ_, NZ_) = Ctmp_1_;
        
        
        
        Ctmp_1_ = (fft2Dfac_*Tky)*C14_;
        fft_.back_2DFull( Ctmp_1_ );
        Ctmp_1_ = By_.asDiagonal()*Ctmp_1_;
        fft_.for_2DFull( Ctmp_1_ );
        dealias(Ctmp_1_);
        Ckl_out[i].block( 2*NZ_, 3*NZ_, NZ_, NZ_) = Ctmp_1_;
        
        
        
        Ctmp_1_ = (fft2Dfac_*Tky)*C21_;
        fft_.back_2DFull( Ctmp_1_ );
        Ctmp_1_ = By_.asDiagonal()*Ctmp_1_;
        Ctmp_2_ = (fft2Dfac_*T2TdzTiLap2Tky).asDiagonal()*C21_ + (fft2Dfac_*TdzP2TiLap2TkxbPLUSTmiLap2TkxbTkyP2).asDiagonal()*C11_;
        fft_.back_2DFull( Ctmp_2_ );
        Ctmp_2_ = dzBy_.asDiagonal()*Ctmp_2_;
        Ctmp_3_ = (fft2Dfac_*TdzTiLap2Tkxb).asDiagonal()*C11_ + (fft2Dfac_*TiLap2Tky).asDiagonal()*C21_;
        fft_.back_2DFull( Ctmp_3_ );
        Ctmp_3_ = dzdzBy_.asDiagonal()*Ctmp_3_;
        Ctmp_1_ = Ctmp_1_ + Ctmp_2_ + Ctmp_3_;
        fft_.for_2DFull( Ctmp_1_ );
        Ctmp_3_ = Ctmp_1_ + TmdzTq.asDiagonal()*C31_;
        dealias(Ctmp_3_);
        Ckl_out[i].block( 3*NZ_, 0, NZ_, NZ_) = Ctmp_3_;
        
        
        
        Ctmp_1_ = (fft2Dfac_*Tky)*C22_;
        fft_.back_2DFull( Ctmp_1_ );
        Ctmp_1_ = By_.asDiagonal()*Ctmp_1_;
        Ctmp_2_ = (fft2Dfac_*T2TdzTiLap2Tky).asDiagonal()*C22_ + (fft2Dfac_*TdzP2TiLap2TkxbPLUSTmiLap2TkxbTkyP2).asDiagonal()*C12_;
        fft_.back_2DFull( Ctmp_2_ );
        Ctmp_2_ = dzBy_.asDiagonal()*Ctmp_2_;
        Ctmp_3_ = (fft2Dfac_*TdzTiLap2Tkxb).asDiagonal()*C12_ + (fft2Dfac_*TiLap2Tky).asDiagonal()*C22_;
        fft_.back_2DFull( Ctmp_3_ );
        Ctmp_3_ = dzdzBy_.asDiagonal()*Ctmp_3_;
        Ctmp_1_ = Ctmp_1_ + Ctmp_2_ + Ctmp_3_;
        fft_.for_2DFull( Ctmp_1_ );
        Ctmp_3_ = Ctmp_1_ + TmdzTq.asDiagonal()*C32_;
        dealias(Ctmp_3_);
        Ckl_out[i].block( 3*NZ_, NZ_, NZ_, NZ_) = Ctmp_3_;
        
        
        
        Ctmp_1_ = (fft2Dfac_*Tky)*C23_;
        fft_.back_2DFull( Ctmp_1_ );
        Ctmp_1_ = By_.asDiagonal()*Ctmp_1_;
        Ctmp_2_ = (fft2Dfac_*T2TdzTiLap2Tky).asDiagonal()*C23_ + (fft2Dfac_*TdzP2TiLap2TkxbPLUSTmiLap2TkxbTkyP2).asDiagonal()*C13_;
        fft_.back_2DFull( Ctmp_2_ );
        Ctmp_2_ = dzBy_.asDiagonal()*Ctmp_2_;
        Ctmp_3_ = (fft2Dfac_*TdzTiLap2Tkxb).asDiagonal()*C13_ + (fft2Dfac_*TiLap2Tky).asDiagonal()*C23_;
        fft_.back_2DFull( Ctmp_3_ );
        Ctmp_3_ = dzdzBy_.asDiagonal()*Ctmp_3_;
        Ctmp_1_ = Ctmp_1_ + Ctmp_2_ + Ctmp_3_;
        fft_.for_2DFull( Ctmp_1_ );
        Ctmp_3_ = Ctmp_1_ + TmdzTq.asDiagonal()*C33_;
        dealias(Ctmp_3_);
        Ckl_out[i].block( 3*NZ_, 2*NZ_, NZ_, NZ_) = Ctmp_3_;
        
        
        
        Ctmp_1_ = (fft2Dfac_*Tky)*C24_;
        fft_.back_2DFull( Ctmp_1_ );
        Ctmp_1_ = By_.asDiagonal()*Ctmp_1_;
        Ctmp_2_ = (fft2Dfac_*T2TdzTiLap2Tky).asDiagonal()*C24_ + (fft2Dfac_*TdzP2TiLap2TkxbPLUSTmiLap2TkxbTkyP2).asDiagonal()*C14_;
        fft_.back_2DFull( Ctmp_2_ );
        Ctmp_2_ = dzBy_.asDiagonal()*Ctmp_2_;
        Ctmp_3_ = (fft2Dfac_*TdzTiLap2Tkxb).asDiagonal()*C14_ + (fft2Dfac_*TiLap2Tky).asDiagonal()*C24_;
        fft_.back_2DFull( Ctmp_3_ );
        Ctmp_3_ = dzdzBy_.asDiagonal()*Ctmp_3_;
        Ctmp_1_ = Ctmp_1_ + Ctmp_2_ + Ctmp_3_;
        fft_.for_2DFull( Ctmp_1_ );
        Ctmp_3_ = Ctmp_1_ + TmdzTq.asDiagonal()*C34_;
        dealias(Ctmp_3_);
        Ckl_out[i].block( 3*NZ_, 3*NZ_, NZ_, NZ_) = Ctmp_3_;
        ///     AUTOMATICALLY GENERATED EQUATIONS - END       ///
        ////////////////////////////////////////////////////////

        
        /////                                ////
        /////////////////////////////////////////
        

        //Since I almost always want the Reynolds stress in a linear calculation better to still calculate and only change the MF feedback
        ///////////////////////////////////////////////////////////
        //////                                               //////
        //////             REYNOLDS STRESSES                 //////
        //////                                               //////
        
//        ubops={[id zs zs zs], [-kxt*ky*ilap2 K.kz.*ilap2 zs zs ],...
//            [-kxt*K.kz.*ilap2 -ky*ilap2 zs zs ],...
//            [zs zs id zs], [zs zs -kxt*ky*ilap2 K.kz.*ilap2],...
//            [zs zs -kxt*K.kz.*ilap2 -ky*ilap2]};
//        bzuxmuzbx(:,xxx,yyy)=real(diag(ifft(ifft(  ubops{6}*ckl*(ubops{1}')  )')')-...
//             diag(ifft(ifft(  ubops{3}*ckl*(ubops{4}')  )')'));
//        bzuymuzby(:,xxx,yyy)=real(diag(ifft(ifft(  ubops{6}*ckl*(ubops{2}')  )')')-...
//             diag(ifft(ifft(  ubops{3}*ckl*(ubops{5}')  )')'));
        // Faster to do it without using ubops, since that involves lots of multiplying zeros
        
        // mult_fac for sum - since no -ve ky values are used, count everything but ky=0 twice
        double mult_fac = 2.0;
        if (kytmp_== 0.0 )
            mult_fac = 1.0; // Only count ky=0 mode once
        
        double ftfac = mult_fac*fft2Dfac_; // NZ^2 factor for ifft is included here - more efficient since scalar
        // These are Eigen matrices instead, just for asDiagonal
        rey_mkxky_tmp_ = (-kyctmp_*kxctmp_ )*ilap2tmp_.cast<dcmplx>().matrix();
        rey_kz_tmp_ = (kz_*ilap2tmp_).matrix();
        rey_mkxkz_tmp_ = ((-kxctmp_ )*kz_*ilap2tmp_).matrix();
        rey_mky_tmp_ = (-kyctmp_ )*ilap2tmp_.cast<dcmplx>().matrix();
        
        // CAN PROBABLY REDUCE THE COMPUTATION BY HALF BY USING Ckl SYMMETRY!
        // bz*ux - uz*bx
        // bz*ux
        reynolds_mat_tmp_ = rey_mkxkz_tmp_.asDiagonal()*Ckl_in[i].block( 2*NZ_, 0, NZ_, NZ_) +rey_mky_tmp_.asDiagonal()*Ckl_in[i].block(3*NZ_, 0,NZ_, NZ_);
        // - uz*bx
        reynolds_mat_tmp_ -= rey_mkxkz_tmp_.asDiagonal()*Ckl_in[i].block(0,2*NZ_, NZ_, NZ_) +
             rey_mky_tmp_.asDiagonal()*Ckl_in[i].block(NZ_, 2*NZ_,NZ_, NZ_) ;
        
        // fft(fft( a )')'
        fft_.back_2DFull(reynolds_mat_tmp_);

        // Keep a running sum
        bzux_m_uzbx_d_ += ftfac*(reynolds_mat_tmp_.diagonal().array().real());
        
        // bz*uy - uz*by
        // bz*uy
        reynolds_mat_tmp_ = rey_mkxkz_tmp_.asDiagonal()*Ckl_in[i].block(0,2*NZ_, NZ_, NZ_)*rey_mkxky_tmp_.conjugate().asDiagonal();
        reynolds_mat_tmp_ += rey_mkxkz_tmp_.asDiagonal()*Ckl_in[i].block(NZ_,2*NZ_, NZ_, NZ_)*rey_kz_tmp_.conjugate().asDiagonal();
        reynolds_mat_tmp_ += rey_mky_tmp_.asDiagonal()*Ckl_in[i].block(0,3*NZ_, NZ_, NZ_)*rey_mkxky_tmp_.conjugate().asDiagonal();
        reynolds_mat_tmp_ += rey_mky_tmp_.asDiagonal()*Ckl_in[i].block(NZ_,3*NZ_, NZ_, NZ_)*rey_kz_tmp_.conjugate().asDiagonal();
        // -uz*by
        reynolds_mat_tmp_ -= rey_mkxkz_tmp_.asDiagonal()*Ckl_in[i].block( 2*NZ_,0, NZ_, NZ_)*rey_mkxky_tmp_.conjugate().asDiagonal();
        reynolds_mat_tmp_ -= rey_mkxkz_tmp_.asDiagonal()*Ckl_in[i].block(3*NZ_, 0, NZ_, NZ_)*rey_kz_tmp_.conjugate().asDiagonal();
        reynolds_mat_tmp_ -= rey_mky_tmp_.asDiagonal()*Ckl_in[i].block(2*NZ_, NZ_, NZ_, NZ_)*rey_mkxky_tmp_.conjugate().asDiagonal();
        reynolds_mat_tmp_ -= rey_mky_tmp_.asDiagonal()*Ckl_in[i].block(3*NZ_, NZ_, NZ_, NZ_)*rey_kz_tmp_.conjugate().asDiagonal();
        // fft(fft( a )')'
        fft_.back_2DFull(reynolds_mat_tmp_);
        
        // Keep a running sum
        bzuy_m_uzby_d_ -= ftfac*(reynolds_mat_tmp_.diagonal().array().real());
        
        
        
        //////                                               //////
        ///////////////////////////////////////////////////////////
    
        
        // Add to adjoint
        Ckl_out[i] += Ckl_out[i].adjoint().eval();
        
//        std::cout << "kx = " << kxtmp_ << ", ky = " << kytmp_ << std::endl;
//        std::cout << Ckl_out[i] << std::endl << std::endl;

        // Add driving noise
        Ckl_out[i] += Qkl_tmp_.cast<dcmplx>().matrix().asDiagonal();
        
        
        ////////////////////////////////////////
        //////   LINEAR PART
        // Need to re-evaluate laplacian, since different time.
        kxtmp_=kxtmp_ + q_*dt_lin*kytmp_;
        // Redefine lapFtmp and lap2tmp, they are not lapF and lap2 so cannot use after this!!!
        lapFtmp_ = nu_*((-kxtmp_*kxtmp_-kytmp_*kytmp_)+ kz2_);
        lap2tmp_ = (eta_/nu_)*lapFtmp_;// Saves some calculation
        linop_Ckl[i] << lapFtmp_, lapFtmp_, lap2tmp_, lap2tmp_;
        ////////////////////////////////////////


        
    }
    // Since I almost always want the Reynolds stress in a linear calculation better to still calculate and only change the MF feedback
    //////////////////////////////////////
    // Reynolds stress
    // Sum accross all processes
        
    // FOR SOME REASON MPI::IN_PLACE IS GIVING INCORRECT RESULTS!
//        mpi_.SumAllReduce_IP_double(bzux_m_uzbx_d_.data(), NZ_);
//        mpi_.SumAllReduce_IP_double(bzuy_m_uzby_d_.data(), NZ_);
    mpi_.SumAllReduce_double(bzux_m_uzbx_d_.data(), bzux_m_uzbx_drec_.data(), NZ_);
    mpi_.SumAllReduce_double(bzuy_m_uzby_d_.data(), bzuy_m_uzby_drec_.data(), NZ_);
    
    // Convert to complex - didn't seem worth the effort of fftw real transforms as time spent here (and in MF transforms) will be very minimal
    
    bzux_m_uzbx_c_ = bzux_m_uzbx_drec_.cast<dcmplx>();
    bzuy_m_uzby_c_ = bzuy_m_uzby_drec_.cast<dcmplx>();
    
    // Calculate stresses in Fourier space
    fft_.for_1D(bzux_m_uzbx_c_.data());
    fft_.for_1D(bzuy_m_uzby_c_.data());
    
    double ftfac = 1.0/(Nxy_[0]*2*Nxy_[1]);ftfac=ftfac*ftfac;
    bzux_m_uzbx_c_ = ftfac*kz_*bzux_m_uzbx_c_;
    bzuy_m_uzby_c_ = ftfac*kz_*bzuy_m_uzby_c_;
    dealias(bzux_m_uzbx_c_);
    dealias(bzuy_m_uzby_c_);
    //////////////////////////////////////


    
    //////////////////////////////////////
    ////   MEAN FIELDS    ////////////////
    // In some integrators MFout will be the same as MFin.
    // Because of this, important to calculate MFout[1] first, since it depends on MFin[1], and this would be updated otherwise
    // Changed this, should be safer now
    if (QL_YN_) {
        MFout[1] = -q_*MFin[0] + bzuy_m_uzby_c_;
        MFout[0] = bzux_m_uzbx_c_;
    } else { // Still calculate Reynolds stress in linear calculation
        MFout[1] = -q_*MFin[0];
        MFout[0].setZero();
    }
    
    
    

}


//////////////////////////////////////////////////////////
///   Initialization of the linear operators        //////

void Model_AutoGen_template::linearOPs_Init(double t0, doubVec * linop_MF, doubVec * linop_Ckl_old) {
    // Constructs initial (t=0) linear operators, since MF operator is constant this completely defines it for entire integration
    // Mean fields
    if (QL_YN_) {
        for (int i=0; i<num_MF_;  ++i){
            linop_MF[i] = eta_*kz2_;
        };
    } else { // Turn off mean dissipation when no QL effects
        for (int i=0; i<num_MF_;  ++i){
            linop_MF[i].setZero();
        };
    }
    
    
    // Ckl

    for (int i=0; i<Cdimxy();  ++i){
        int k_i = i + index_for_k_array(); // k index
        kytmp_=ky_[k_i].imag();
        kxtmp_=kx_[k_i].imag() + q_*t0*kytmp_;
        lapFtmp_ = -kxtmp_*kxtmp_-kytmp_*kytmp_+kz2_;
        linop_Ckl_old[i] << nu_*lapFtmp_, nu_*lapFtmp_, eta_*lapFtmp_, eta_*lapFtmp_;
    }
    
}


//////////////////////////////////////////////////////////
////                DEALIASING                      //////


// Applies dealiasing to 2-D matrix of size NZ_ in z Fourier space
void Model_AutoGen_template::dealias(dcmplxMat &inMat) {
    dcmplx * arr = inMat.data();
    // Could also use Block here, might be faster
    for (int j=delaiasBnds_[0]; j<delaiasBnds_[1]+1; ++j){
        for (int i=0; i<NZ_; ++i) {
            arr[i+NZ_*j]=0; // Column major, not that it makes any difference
        }
    }
    
    for (int j=0; j<NZ_; ++j) {
        for (int i=delaiasBnds_[0]; i<delaiasBnds_[1]+1; ++i){
            arr[i+NZ_*j]=0; // Column major, not that it makes any difference
        }
    }
}

//1-D version - takes a dcmplxVec input
void Model_AutoGen_template::dealias(dcmplxVec& vec) {
    for (int j=delaiasBnds_[0]; j<delaiasBnds_[1]+1; ++j){
        vec(j) = 0;
    }
}

//1-D version - takes a doubVec input
void Model_AutoGen_template::dealias(doubVec& vec) {
    for (int j=delaiasBnds_[0]; j<delaiasBnds_[1]+1; ++j){
        vec(j) = 0;
    }
}
//////////////////////////////////////////////////////////







//////////////////////////////////////////////////////////
//////                                              //////
//////          ENERGY, ANGULAR MOMENTUM etc.       //////
//////                                              //////
//////////////////////////////////////////////////////////

void Model_AutoGen_template::Calc_Energy_AM_Diss(TimeVariables& tv, double t, const dcmplxVec *MFin, const dcmplxMat *Cin ) {
    // Energy, angular momentum and dissipation of the solution MFin and Cin
    // TimeVariables class stores info and data about energy, AM etc.
    // t is time
    
    // OUTPUT: energy[1] and [2] contain U and B mean field energies (energy[1]=0 for this)
    // energy[3] and [4] contain u and b fluctuating energies.
    
    // Energy
    double energy_u=0, energy_b=0;
    double energy_u_f=0, energy_b_f=0; // Receive buffer for MPI_Reduce
    // Angular momentum
    double AM_u = 0, AM_b = 0;
    double AM_u_f=0, AM_b_f=0; // Receive buffer for MPI_Reduce

    int mult_fac;// Factor to account for only doing half fft sum.
    /////////////////////////////////////
    //   ALL OF THIS IS PARALLELIZED
    for (int i=0; i<Cdimxy();  ++i){
        
        int k_i = i + index_for_k_array(); // k index
        int ind_ky = ky_index_[k_i];
        // Form Laplacians using time-dependent kx
        kyctmp_ = ky_[k_i];
        kytmp_ = kyctmp_.imag(); // It is useful to have both complex and double versions
        kxctmp_ = kx_[k_i] + q_*t*kyctmp_;
        kxtmp_ = kxctmp_.imag();
        
        ilap2tmp_ = ilap2_[ind_ky];
        lapFtmp_ = -kxtmp_*kxtmp_ + lap2_[ind_ky];
        
        mult_fac = 2;
        if (kytmp_== 0.0 )
            mult_fac = 1; // Only count ky=0 mode once

        if (tv.energy_save_Q()){
            //////////////////////////////////////
            //     ENERGY
            // Use Qkl_tmp_ for Mkl to save memory
            lapFtmp_ = lapFtmp_*ilap2tmp_;
            Qkl_tmp_ << lapFtmp_, -ilap2tmp_,lapFtmp_, -ilap2tmp_;
            
            // Energy = trace(Mkl*Ckl), Mkl is diagonal
            Qkl_tmp_ = mult_fac*Qkl_tmp_ * Cin[i].real().diagonal().array();
            
            int num_u_b = num_fluct_/2; // Keep it clear where numbers are coming from.
            energy_u += Qkl_tmp_.head(NZ_*num_u_b).sum();
            energy_b += Qkl_tmp_.tail(NZ_*num_u_b).sum();
            //
            //////////////////////////////////////
        }
        if (tv.AngMom_save_Q()){
            //////////////////////////////////////
            //   Angular momentum
            
            // These are Eigen matrices instead, just for asDiagonal
            rey_mkxky_tmp_ = (-kyctmp_*kxctmp_ )*ilap2tmp_.cast<dcmplx>().matrix();
            rey_kz_tmp_ = (kz_*ilap2tmp_).matrix();
            
            // uy*ux - steal Reynolds stress variable for this
            bzux_m_uzbx_d_ = (rey_mkxky_tmp_.array()*Cin[i].block( 0, 0, NZ_, NZ_).diagonal().array()).real() + (rey_kz_tmp_.array()*Cin[i].block(NZ_, 0,NZ_, NZ_).diagonal().array()).real();
            AM_u += mult_fac*bzux_m_uzbx_d_.sum();
            // by*bx - steal Reynolds stress variable for this
            bzux_m_uzbx_d_ = (rey_mkxky_tmp_.array()*Cin[i].block( 2*NZ_,2*NZ_, NZ_, NZ_).diagonal().array()).real() +   (rey_kz_tmp_.array()*Cin[i].block(3*NZ_, 2*NZ_,NZ_, NZ_).diagonal().array()).real();
            AM_b += mult_fac*bzux_m_uzbx_d_.sum();
            //
            //////////////////////////////////////
        }
        
        
    }
    if (tv.energy_save_Q()){
        // Put the energy on processor 0
        mpi_.SumReduce_doub(&energy_u,&energy_u_f,1);
        mpi_.SumReduce_doub(&energy_b,&energy_b_f,1);
    }
    if (tv.AngMom_save_Q()){
        // Put the angular momentum on processor 0
        mpi_.SumReduce_doub(&AM_u,&AM_u_f,1);
        mpi_.SumReduce_doub(&AM_b,&AM_b_f,1);
    }
//    energy_u_f = energy_u; energy_b_f = energy_b;   // If in-place version is desired
//    mpi_.SumReduce_IP_doub(&energy_u_f,1);
//    mpi_.SumReduce_IP_doub(&energy_b_f,1);
    
    // Currently, TimeVariables is set to save on root process, may want to generalize this at some point (SumReduce_doub is also)
    if (mpi_.my_n_v() == 0) {
        ////////////////////////////////////////////
        ///// All this is only on processor 0  /////
        double divfac=1.0/totalN2_;
        energy_u_f = energy_u_f*divfac;
        energy_b_f = energy_b_f*divfac;
        AM_u_f = AM_u_f*divfac;
        AM_b_f = AM_b_f*divfac;
        
        
        ///////////////////////////////////////
        ///       MEAN FIELDS            //////
        // Only need to calculate on one processor
        double energy_MU=0, energy_MB=0;
        double AM_MU =0, AM_MB=0;
        
        energy_MB = (MFin[0].abs2().sum() + MFin[1].abs2().sum())/(NZ_*NZ_);
        AM_MB = (MFin[0]*MFin[1].conjugate()).real().sum()/(NZ_*NZ_);
        
        ///////////////////////////////////////
        //////         OUTPUT            //////
        
        // Energy
        double* en_point = tv.current_energy();
        en_point[0] = energy_MU/2;
        en_point[1] = energy_MB/2;
        en_point[2] = energy_u_f/2;
        en_point[3] = energy_b_f/2;
        
        // Angular momentum
        double* AM_point = tv.current_AM();
        AM_point[0] = AM_MU/2;
        AM_point[1] = AM_MB/2;
        AM_point[2] = AM_u_f/2;
        AM_point[3] = AM_b_f/2;
        
        
        
        if (tv.reynolds_save_Q()) {
            //////////////////////////////////////
            //   Reynolds stress
            // Saves quantities to do with the reynolds stress and dynamo
            // 1) Shear contribution: Re( -q Bx(k0) By(k0))/By(k0)
            // 2) y emf: Re( bzuy_m_uzby_c_*By(k0) )/By(k0)
            // 3) y dissipation: eta*k0^2*By(k0)
            // 4) x emf: Re( bzux_m_uzbx_c_*Bx(k0) )/Bx(k0)
            // 5) x dissipation: eta*k0^2*Bx(k0)
            // Have assumed k0 to be lowest kz! i.e., driving largest dynamo possible in the box
            // There may be slight errors here from saving using the updated values of MFin, presumably this is a small effect, especially at high resolution (low dt) and in steady state.
            double* rey_point = tv.current_reynolds();
            rey_point[0] = real(-q_*MFin[0](1)*conj(MFin[1](1)))/abs(MFin[1](1) );
            rey_point[1] = real(bzuy_m_uzby_c_(1)*conj(MFin[1](1)))/abs(MFin[1](1) );
            rey_point[2] = eta_ *kz2_(1)*abs(MFin[1](1) );
            rey_point[3] = real(bzux_m_uzbx_c_(1)*conj(MFin[0](1)))/abs(MFin[0](1) );
            rey_point[4] = eta_ * kz2_(1)*abs(MFin[0](1) );
//            rey_point[0] = real(bzuy_m_uzby_c_(1)*conj(MFin[1](1))) ;
//            rey_point[1] = real(bzux_m_uzbx_c_(1)*conj(MFin[1](1))) ;
            //
            //////////////////////////////////////
        }

        ///// All this is only on processor 0  /////
        ////////////////////////////////////////////
    }
    
    // Save the data
    tv.Save_Data();
    
}




//////////////////////////////////////////////////////
//         PRIVATE FUNCTIONS FOR INITIALIZATION     //

// Create and store arrays of lap2 to save computation of ffts
void Model_AutoGen_template::Define_Lap2_Arrays_(void){
    // Lap2 quantities are stored on all processors, even though I could save some memory here. Easy to change later if a problem, the amount of memory may actually be significant...
    
    // Maximum ky index
    int number_of_ky = ky_index_[Cdimxy_full()-1]+1;
    // Quick error check
    if (Cdimxy_full()%number_of_ky!=0)
        std::cout << "Warning, something funny going on in Define_Lap2_Arrays_" << std::endl;
    // Assign data to arrays
    lap2_ = new doubVec[ number_of_ky ];
    ilap2_ = new doubVec[ number_of_ky ];
    
    for (int i=0; i<number_of_ky;  ++i){
        
        // Form Laplacian
        kytmp_=ky_[i].imag();
        lap2tmp_ = -kytmp_*kytmp_ + kz2_;
        ilap2tmp_ = 1/lap2tmp_;
        
        // Fix up zero bits
        if (kytmp_==0 ) {
            lap2tmp_(0)=0;
            // Avoid infinities
            ilap2tmp_(0)=1;
        }
        
        // Assign to lap2_ and ilap2_
        lap2_[i] = lap2tmp_;
        ilap2_[i] = ilap2tmp_;
        
    }
    
}

//////////////////////////////////////////////////////







