//
//  MHD_FullUBQlin.cc
//  MRIDSS
//
//  Created by Jonathan Squire on 4/25/14.
//  Copyright (c) 2014 J Squire. All rights reserved.
//

#include "MHD_FullUBQlin.h"

// Model class for S3T/CE2 shearing box MHD model
// Template (not in c++ sense) for the equations to be copied into from automatically generated Mathematica file
//
// Derived from Model (model.h)

// Constructor for MHD_FullUBQlin
MHD_FullUBQlin::MHD_FullUBQlin(const Inputs& sp, MPIdata& mpi, fftwPlans& fft) :
equations_name("MHD_FullUBQlin"),
num_MF_(4), num_fluct_(4),
q_(sp.q), Omega_(sp.omega), // Rotation
f_noise_(sp.f_noise), nu_(sp.nu), eta_(sp.eta),B0z_(sp.B0z),// Mean B0z
QL_YN_(sp.QuasiLinearQ),
dont_drive_ky0_modes_Q_ (0),// turn off driving of ky=0, kx,kz != 0 modes (i.e., non-shearing waves)
drive_only_velocity_fluctuations_(0), // If True, only drive velocity field.
velocity_noise_mult_(sp.v_noise_mult),magnetic_noise_mult_(sp.b_noise_mult), //Multiply each noise by these numbers
Model(sp.NZ, sp.NXY , sp.L), // Dimensions - stored in base
mpi_(mpi), // MPI data
fft_(fft) // FFT data
{
    // Check that model is that specified in input file
    if (sp.equations_to_use != equations_name) {
        std::stringstream error_str;
        error_str << "Model name, " << equations_name << ", does not match that specified in input file, " << sp.equations_to_use << "!!" << std::endl;
        mpi.print1( error_str.str() );
        ABORT;
    }
    
    // Dimensions that depend on model
    nz_Cfull_=num_fluct_*NZ_;
    // Calcuate kx and ky arrays - single dimension NX*NY(/2) to simplify parallelization
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
    // Factor for xy transform (in Reynolds stress)
    fftFac_Reynolds_ = 1.0/(Nxy_[0]*2*Nxy_[1]);
    fftFac_Reynolds_=fftFac_Reynolds_*fftFac_Reynolds_;
    
    
    ////////////////////////////////////////////////////
    // Useful arrays to save computation
    Define_Lap2_Arrays_(); // Lap2- Allocates data to lap2_,ilap2_
    
    
    // CFL calculation
    kmax = 0;
    double dealiasfac[3] = {2.0, 2.0,3.0};
    int N[3] = {Nxy_[0],2*Nxy_[1],NZ_};
    for (int i=0; i<3; ++i) {
        if (kmax < 2*PI/box_length(i)*(N[i]/dealiasfac[i]))
            kmax = 2*PI/box_length(i)*(N[i]/dealiasfac[i]);
    }
    
    // Noise cutoff - compare to laplacian so squared
    noise_range_[0] = sp.noise_range_low*sp.noise_range_low;
    noise_range_[1] = sp.noise_range_high*sp.noise_range_high;
    drive_condition_ = Eigen::Array<bool,Eigen::Dynamic,1>(NZ_);
    print_noise_range_();
    
    // Reynolds stresses
    bzux_m_uzbx_c_ = dcmplxVec::Zero(NZ_); // Bx (take dz after)
    bzuy_m_uzby_c_ = dcmplxVec::Zero(NZ_); // By (take dz after)
    ux_rey_stress_c_ = dcmplxVec::Zero(NZ_); // Ux
    uy_rey_stress_c_ = dcmplxVec::Zero(NZ_); // Uy
    bzux_m_uzbx_d_ = doubVec::Zero(NZ_);// Bx (take dz after)
    bzuy_m_uzby_d_ = doubVec::Zero(NZ_); // By (take dz after)
    ux_rey_stress_d_ = doubVec::Zero(NZ_); // Ux
    uy_rey_stress_d_ = doubVec::Zero(NZ_); // Uy

    // Send/receive buffers for MPI, for some reason MPI::IN_PLACE is not giving the correct result!
    reynolds_stress_MPI_send_ = doubVec::Zero(num_MF_*NZ_);
    reynolds_stress_MPI_receive_ = doubVec::Zero(num_MF_*NZ_);
    
    reynolds_save_tmp_ = new double[5];
    
    
    // Random settings, print to avoid mistakes
    std::stringstream prnt;
    prnt << "Rotation set to " << Omega_ <<"\n";
    mpi_.print1(prnt.str());
    if (dont_drive_ky0_modes_Q_)
        mpi_.print1("ky=0 modes are excluded from this calculation!\n");
    if (drive_only_velocity_fluctuations_)
        mpi_.print1("No Magnetic driving noise!\n");
    
    
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


    // Real versions of B - need 3rd order Bx derivatives for full B Quasi-linear
    By_ = dcmplxVecM::Zero(NZ_);
    dzBy_ = dcmplxVecM::Zero(NZ_);
    dzdzBy_ = dcmplxVecM::Zero(NZ_);
    Bx_ = dcmplxVecM::Zero(NZ_);
    dzBx_ = dcmplxVecM::Zero(NZ_);
    dzdzBx_ = dcmplxVecM::Zero(NZ_);
    dzdzdzBx_ = dcmplxVecM::Zero(NZ_);
    // Real versions of U - need 3rd order Ux derivatives for full B Quasi-linear
    Uy_ = dcmplxVecM::Zero(NZ_);
    dzUy_ = dcmplxVecM::Zero(NZ_);
    dzdzUy_ = dcmplxVecM::Zero(NZ_);
    Ux_ = dcmplxVecM::Zero(NZ_);
    dzUx_ = dcmplxVecM::Zero(NZ_);
    dzdzUx_ = dcmplxVecM::Zero(NZ_);
    dzdzdzUx_ = dcmplxVecM::Zero(NZ_);

    
    // Reynolds stresses - this need to be changed to automatic also
    reynolds_mat_tmp_ = dcmplxMat(NZ_,NZ_);
    // Automatically generated temporary variables - class constructor
    rey_TdzTiLap2TkxbTky = dcmplxVecM(NZ_);
    rey_TiLap2TkxbTkyP2 = dcmplxVecM(NZ_);
    rey_TiLap2TkxbP2Tky = dcmplxVecM(NZ_);
    rey_TiLap2TkxbTky = dcmplxVecM(NZ_);
    rey_TdzTiLap2Tkxb = dcmplxVecM(NZ_);
    rey_TdzTiLap2Tky = dcmplxVecM(NZ_);
    rey_TdzP2TiLap2 = dcmplxVecM(NZ_);
    rey_TiLap2Tky = dcmplxVecM(NZ_);
    rey_TdzTiLap2 = dcmplxVecM(NZ_);
    rey_Tdz = dcmplxVecM(NZ_);
    
    
    /////////////////////////////////////////////////
    //  AUTO GENERATED VARIABLES
    
    
    // Automatically generated temporary variables - class constructor
    T2Tdz = dcmplxVecM(NZ_);
    T2TdzP2TiLap2Tkxb = dcmplxVecM(NZ_);
    T2TdzTiLap2TkxbP2Tky = dcmplxVecM(NZ_);
    T2TdzTiLap2Tky = dcmplxVecM(NZ_);
    T2TiLap2TkxbTkyP2 = dcmplxVecM(NZ_);
    TB0zTdz = dcmplxVecM(NZ_);
    TdzP2TiLap2TkxbPLUSTmiLap2TkxbTkyP2 = dcmplxVecM(NZ_);
    TdzTiLap2Tkxb = dcmplxVecM(NZ_);
    TdzTiLap2TkxbP3PLUSTmdzTkxb = dcmplxVecM(NZ_);
    TdzTkxbPLUSTmdzTiLap2TkxbP3 = dcmplxVecM(NZ_);
    TdzTqPLUSTm2Tdz = dcmplxVecM(NZ_);
    TiLap2TkxbP2Tky = dcmplxVecM(NZ_);
    TiLap2TkxbTkyP2PLUSTmdzP2TiLap2Tkxb = dcmplxVecM(NZ_);
    TiLap2Tky = dcmplxVecM(NZ_);
    Tm2TdzP2TiLap2Tkxb = dcmplxVecM(NZ_);
    Tm2TdzTiLap2TkxbP2Tky = dcmplxVecM(NZ_);
    Tm2TdzTiLap2Tky = dcmplxVecM(NZ_);
    Tm2TiLap2TkxbTkyP2 = dcmplxVecM(NZ_);
    TmdzTiLap2Tkxb = dcmplxVecM(NZ_);
    TmdzTq = dcmplxVecM(NZ_);
    TmiLap2TkxbP2Tky = dcmplxVecM(NZ_);
    TmiLap2Tky = dcmplxVecM(NZ_);
    
    Ctmp_1_ = dcmplxMat(NZ_,NZ_);
    Ctmp_2_ = dcmplxMat(NZ_,NZ_);
    Ctmp_3_ = dcmplxMat(NZ_,NZ_);
    Ctmp_4_ = dcmplxMat(NZ_,NZ_);
    Ctmp_5_ = dcmplxMat(NZ_,NZ_);
    Ctmp_6_ = dcmplxMat(NZ_,NZ_);
    Ctmp_7_ = dcmplxMat(NZ_,NZ_);
    Ctmp_8_ = dcmplxMat(NZ_,NZ_);
    Ctmp_9_ = dcmplxMat(NZ_,NZ_);
    Ctmp_10_ = dcmplxMat(NZ_,NZ_);
    Ctmp_11_ = dcmplxMat(NZ_,NZ_);
    Ctmp_12_ = dcmplxMat(NZ_,NZ_);
    
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
    /////////////////////////////////////////////////////////////


}


// Destructor
MHD_FullUBQlin::~MHD_FullUBQlin(){
    // Here I am destroying Base members from the derived destructor - is this bad?
    delete[] kx_;
    delete[] ky_;
    
    delete[] lap2_;
    delete[] ilap2_;
    // Data arrays
    delete[] reynolds_save_tmp_;
}


// Main EulerFluid function, fx = f(t,x) called at ecah time-step
void MHD_FullUBQlin::rhs(double t, double dt_lin,
                   dcmplxVec *MFin, dcmplxMat *Ckl_in,
                   dcmplxVec *MFout, dcmplxMat *Ckl_out,
                   doubVec * linop_Ckl) {
    
    // If desired, the function should work fine with Ckl_in = Ckl_out. This is
    // useful for integration schemes that do not require re-use of Ckl_in
    
    /////////////////////////////////////////////////////
    // Calculate MFs in real space
    // MF vector is (Bx, By, Ux, Uy) for compatibility with other models
    {
    // By
    By_ = MFin[1].matrix()*fft1Dfac_; // fft back doesn't include normalization
    fft_.back_1D(By_.data());
    dzBy_ = (kz_*MFin[1]).matrix()*fft1Dfac_;
    fft_.back_1D(dzBy_.data());
    dzdzBy_ = (kz2_*MFin[1]).matrix()*fft1Dfac_;
    fft_.back_1D(dzdzBy_.data());
    // Bx
    Bx_ = MFin[0].matrix()*fft1Dfac_; // fft back doesn't include normalization
    fft_.back_1D(Bx_.data());
    dzBx_ = (kz_*MFin[0]).matrix()*fft1Dfac_;
    fft_.back_1D(dzBx_.data());
    dzdzBx_ = (kz2_*MFin[0]).matrix()*fft1Dfac_;
    fft_.back_1D(dzdzBx_.data());
    dzdzdzBx_ = (kz_*kz2_*MFin[0]).matrix()*fft1Dfac_;
    fft_.back_1D(dzdzdzBx_.data());
    // Uy
    Uy_ = MFin[3].matrix()*fft1Dfac_; // fft back doesn't include normalization
    fft_.back_1D(Uy_.data());
    dzUy_ = (kz_*MFin[3]).matrix()*fft1Dfac_;
    fft_.back_1D(dzUy_.data());
    dzdzUy_ = (kz2_*MFin[3]).matrix()*fft1Dfac_;
    fft_.back_1D(dzdzUy_.data());
    // Ux
    Ux_ = MFin[2].matrix()*fft1Dfac_; // fft back doesn't include normalization
    fft_.back_1D(Ux_.data());
    dzUx_ = (kz_*MFin[2]).matrix()*fft1Dfac_;
    fft_.back_1D(dzUx_.data());
    dzdzUx_ = (kz2_*MFin[2]).matrix()*fft1Dfac_;
    fft_.back_1D(dzdzUx_.data());
    dzdzdzUx_ = (kz_*kz2_*MFin[2]).matrix()*fft1Dfac_;
    fft_.back_1D(dzdzdzUx_.data());
    }
    /////////////////////////////////////////////////////

    
    // Reynolds stresses -- added to at each step in the loop
    bzux_m_uzbx_d_.setZero();
    bzuy_m_uzby_d_.setZero();
    ux_rey_stress_d_.setZero();
    uy_rey_stress_d_.setZero();
    
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
        ilap2tmp_ = ilap2_[ind_ky];
        
        drive_condition_ = lapFtmp_> -noise_range_[1] && lapFtmp_< -noise_range_[0];
        
        ////////////////////////////////////////
        //              FORM Qkl              //
        {
        bool Qkl_Zero_Q = kytmp_==0 && kxtmp_==0; // Whether Qkl  should be zero for this mode
        if (dont_drive_ky0_modes_Q_ && kytmp_==0)
            Qkl_Zero_Q = 1;  // Option in the code, set to zero for all ky=0 modes if set
        
        if (Qkl_Zero_Q) {
            Qkl_tmp_.setZero(); // Don't drive the mean components! (really could leave mean cmpts out entirely)
            ilapFtmp_(0) = 1; // Avoid infinities (lap2 already done, since stored)
        }
        else {
            lapFtmp_ = drive_condition_.cast<double>()*lap2tmp_/lapFtmp_; // lapFtmp_ is no longer lapF!!!
            lap2tmp_ *= drive_condition_.cast<double>();
            dealias(lapFtmp_);
            dealias(lap2tmp_);
            // CANNOT USE AFTER THIS!
            if (kytmp_==0.0) { // Don't drive the ky=kz=0 component (variables not invertible)
                lapFtmp_(0)=0;
            }
            
            Qkl_tmp_ << velocity_noise_mult_*lapFtmp_, velocity_noise_mult_*lap2tmp_, magnetic_noise_mult_*lapFtmp_, magnetic_noise_mult_*lap2tmp_;
            
            Qkl_tmp_ = (f_noise_*f_noise_*totalN2_*mult_noise_fac_)*Qkl_tmp_.abs();
            
            if( drive_only_velocity_fluctuations_ ) // Set magnetic fluctuations to zero
                Qkl_tmp_.segment(2*NZ_, 2*NZ_).setZero();

        }
        }
        ////////////////////////////////////////
        
        
        /////////////////////////////////////////
        ///    INSERT VARIABLE DEFS HERE
        {
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
        Tkxb = (kxctmp_)*fft2Dfac_;
        Tky = (kyctmp_)*fft2Dfac_;
        Tm2TkxbTkyTq = ((-2.*(kxctmp_*kyctmp_))*q_)*fft2Dfac_;
        Tmkxb = ((-1.*kxctmp_))*fft2Dfac_;
        Tmky = ((-1.*kyctmp_))*fft2Dfac_;
        
        // Vector variable definition (automatic)
        T2Tdz = ((2.*Omega_)*kz_.matrix())*fft2Dfac_;
        T2TdzP2TiLap2Tkxb = ((2.*kxctmp_)*(ilap2tmp_*kz2_).matrix())*fft2Dfac_;
        T2TdzTiLap2TkxbP2Tky = ((2.*(kyctmp_*pow(kxctmp_,2)))*(ilap2tmp_*kz_).matrix())*fft2Dfac_;
        T2TdzTiLap2Tky = ((2.*kyctmp_)*(ilap2tmp_*kz_).matrix())*fft2Dfac_;
        T2TiLap2TkxbTkyP2 = ((2.*(kxctmp_*pow(kyctmp_,2)))*ilap2tmp_.matrix())*fft2Dfac_;
        TB0zTdz = (B0z_*kz_.matrix())*fft2Dfac_;
        TdzP2TiLap2TkxbPLUSTmiLap2TkxbTkyP2 = ((-1.*(kxctmp_*pow(kyctmp_,2)))*ilap2tmp_.matrix() + (ilap2tmp_*kz2_).matrix()*kxctmp_)*fft2Dfac_;
        TdzTiLap2Tkxb = ((ilap2tmp_*kz_).matrix()*kxctmp_)*fft2Dfac_;
        TdzTiLap2TkxbP3PLUSTmdzTkxb = ((-1.*kxctmp_)*kz_.matrix() + (ilap2tmp_*kz_).matrix()*pow(kxctmp_,3))*fft2Dfac_;
        TdzTkxbPLUSTmdzTiLap2TkxbP3 = ((-1.*pow(kxctmp_,3))*(ilap2tmp_*kz_).matrix() + kxctmp_*kz_.matrix())*fft2Dfac_;
        TdzTqPLUSTm2Tdz = ((-2.*Omega_)*kz_.matrix() + kz_.matrix()*q_)*fft2Dfac_;
        TiLap2TkxbP2Tky = (ilap2tmp_.matrix()*(kyctmp_*pow(kxctmp_,2)))*fft2Dfac_;
        TiLap2TkxbTkyP2PLUSTmdzP2TiLap2Tkxb = ((-1.*kxctmp_)*(ilap2tmp_*kz2_).matrix() + ilap2tmp_.matrix()*(kxctmp_*pow(kyctmp_,2)))*fft2Dfac_;
        TiLap2Tky = (ilap2tmp_.matrix()*kyctmp_)*fft2Dfac_;
        Tm2TdzP2TiLap2Tkxb = ((-2.*kxctmp_)*(ilap2tmp_*kz2_).matrix())*fft2Dfac_;
        Tm2TdzTiLap2TkxbP2Tky = ((-2.*(kyctmp_*pow(kxctmp_,2)))*(ilap2tmp_*kz_).matrix())*fft2Dfac_;
        Tm2TdzTiLap2Tky = ((-2.*kyctmp_)*(ilap2tmp_*kz_).matrix())*fft2Dfac_;
        Tm2TiLap2TkxbTkyP2 = ((-2.*(kxctmp_*pow(kyctmp_,2)))*ilap2tmp_.matrix())*fft2Dfac_;
        TmdzTiLap2Tkxb = ((-1.*kxctmp_)*(ilap2tmp_*kz_).matrix())*fft2Dfac_;
        TmdzTq = ((-1.*q_)*kz_.matrix())*fft2Dfac_;
        TmiLap2TkxbP2Tky = ((-1.*(kyctmp_*pow(kxctmp_,2)))*ilap2tmp_.matrix())*fft2Dfac_;
        TmiLap2Tky = ((-1.*kyctmp_)*ilap2tmp_.matrix())*fft2Dfac_;
        }
        
        
        
        ////////////////////////////////////////////////////////
        ///       AUTOMATICALLY GENERATED EQUATIONS         /////
        ///       see GenerateC++Equations.nb in MMA        /////
        {
            ////////////////////////////////////////////////////////
            ///       AUTOMATICALLY GENERATED EQUATIONS         /////
            ///       see GenerateC++Equations.nb in MMA        /////
            
            Ctmp_1_ = Tkxb*C31_;
            fft_.back_2DFull( Ctmp_1_ );
            Ctmp_1_ = Bx_.asDiagonal()*Ctmp_1_;
            Ctmp_2_ = Tky*C31_;
            fft_.back_2DFull( Ctmp_2_ );
            Ctmp_2_ = By_.asDiagonal()*Ctmp_2_;
            Ctmp_3_ = Tmkxb*C11_;
            fft_.back_2DFull( Ctmp_3_ );
            Ctmp_3_ = Ux_.asDiagonal()*Ctmp_3_;
            Ctmp_4_ = Tmky*C11_;
            fft_.back_2DFull( Ctmp_4_ );
            Ctmp_4_ = Uy_.asDiagonal()*Ctmp_4_;
            Ctmp_5_ = T2TdzP2TiLap2Tkxb.asDiagonal()*C11_ + T2TdzTiLap2Tky.asDiagonal()*C21_;
            fft_.back_2DFull( Ctmp_5_ );
            Ctmp_5_ = dzdzUx_.asDiagonal()*Ctmp_5_;
            Ctmp_6_ = T2TdzTiLap2TkxbP2Tky.asDiagonal()*C31_ + T2TiLap2TkxbTkyP2.asDiagonal()*C41_;
            fft_.back_2DFull( Ctmp_6_ );
            Ctmp_6_ = dzBy_.asDiagonal()*Ctmp_6_;
            Ctmp_7_ = TdzTiLap2Tkxb.asDiagonal()*C11_ + TiLap2Tky.asDiagonal()*C21_;
            fft_.back_2DFull( Ctmp_7_ );
            Ctmp_7_ = dzdzdzUx_.asDiagonal()*Ctmp_7_;
            Ctmp_8_ = Tm2TdzP2TiLap2Tkxb.asDiagonal()*C31_ + Tm2TdzTiLap2Tky.asDiagonal()*C41_;
            fft_.back_2DFull( Ctmp_8_ );
            Ctmp_8_ = dzdzBx_.asDiagonal()*Ctmp_8_;
            Ctmp_9_ = Tm2TdzTiLap2TkxbP2Tky.asDiagonal()*C11_ + Tm2TiLap2TkxbTkyP2.asDiagonal()*C21_;
            fft_.back_2DFull( Ctmp_9_ );
            Ctmp_9_ = dzUy_.asDiagonal()*Ctmp_9_;
            Ctmp_10_ = TmdzTiLap2Tkxb.asDiagonal()*C31_ + TmiLap2Tky.asDiagonal()*C41_;
            fft_.back_2DFull( Ctmp_10_ );
            Ctmp_10_ = dzdzdzBx_.asDiagonal()*Ctmp_10_;
            Ctmp_11_ = TdzTiLap2TkxbP3PLUSTmdzTkxb.asDiagonal()*C31_ + TiLap2TkxbP2Tky.asDiagonal()*C41_ + Tmky*C41_;
            fft_.back_2DFull( Ctmp_11_ );
            Ctmp_11_ = dzBx_.asDiagonal()*Ctmp_11_;
            Ctmp_12_ = TdzTkxbPLUSTmdzTiLap2TkxbP3.asDiagonal()*C11_ + Tky*C21_ + TmiLap2TkxbP2Tky.asDiagonal()*C21_;
            fft_.back_2DFull( Ctmp_12_ );
            Ctmp_12_ = dzUx_.asDiagonal()*Ctmp_12_;
            Ctmp_1_ = Ctmp_1_ + Ctmp_2_ + Ctmp_3_ + Ctmp_4_;
            fft_.for_2DFull( Ctmp_1_ );
            Ctmp_5_ = Ctmp_10_ + Ctmp_11_ + Ctmp_12_ + Ctmp_5_ + Ctmp_6_ + Ctmp_7_ + Ctmp_8_ + Ctmp_9_;
            fft_.for_2DFull( Ctmp_5_ );
            Ctmp_12_ = Ctmp_1_ + ilapFtmp_.matrix().asDiagonal()*(Ctmp_5_ + (T2Tdz/fft2Dfac_).asDiagonal()*C21_ + (Tm2TkxbTkyTq/fft2Dfac_)*C11_) + (TB0zTdz/fft2Dfac_).asDiagonal()*C31_;
            dealias(Ctmp_12_);
            Ckl_out[i].block( 0, 0, NZ_, NZ_) = Ctmp_12_;
            
            
            
            Ctmp_1_ = Tkxb*C32_;
            fft_.back_2DFull( Ctmp_1_ );
            Ctmp_1_ = Bx_.asDiagonal()*Ctmp_1_;
            Ctmp_2_ = Tky*C32_;
            fft_.back_2DFull( Ctmp_2_ );
            Ctmp_2_ = By_.asDiagonal()*Ctmp_2_;
            Ctmp_3_ = Tmkxb*C12_;
            fft_.back_2DFull( Ctmp_3_ );
            Ctmp_3_ = Ux_.asDiagonal()*Ctmp_3_;
            Ctmp_4_ = Tmky*C12_;
            fft_.back_2DFull( Ctmp_4_ );
            Ctmp_4_ = Uy_.asDiagonal()*Ctmp_4_;
            Ctmp_5_ = T2TdzP2TiLap2Tkxb.asDiagonal()*C12_ + T2TdzTiLap2Tky.asDiagonal()*C22_;
            fft_.back_2DFull( Ctmp_5_ );
            Ctmp_5_ = dzdzUx_.asDiagonal()*Ctmp_5_;
            Ctmp_6_ = T2TdzTiLap2TkxbP2Tky.asDiagonal()*C32_ + T2TiLap2TkxbTkyP2.asDiagonal()*C42_;
            fft_.back_2DFull( Ctmp_6_ );
            Ctmp_6_ = dzBy_.asDiagonal()*Ctmp_6_;
            Ctmp_7_ = TdzTiLap2Tkxb.asDiagonal()*C12_ + TiLap2Tky.asDiagonal()*C22_;
            fft_.back_2DFull( Ctmp_7_ );
            Ctmp_7_ = dzdzdzUx_.asDiagonal()*Ctmp_7_;
            Ctmp_8_ = Tm2TdzP2TiLap2Tkxb.asDiagonal()*C32_ + Tm2TdzTiLap2Tky.asDiagonal()*C42_;
            fft_.back_2DFull( Ctmp_8_ );
            Ctmp_8_ = dzdzBx_.asDiagonal()*Ctmp_8_;
            Ctmp_9_ = Tm2TdzTiLap2TkxbP2Tky.asDiagonal()*C12_ + Tm2TiLap2TkxbTkyP2.asDiagonal()*C22_;
            fft_.back_2DFull( Ctmp_9_ );
            Ctmp_9_ = dzUy_.asDiagonal()*Ctmp_9_;
            Ctmp_10_ = TmdzTiLap2Tkxb.asDiagonal()*C32_ + TmiLap2Tky.asDiagonal()*C42_;
            fft_.back_2DFull( Ctmp_10_ );
            Ctmp_10_ = dzdzdzBx_.asDiagonal()*Ctmp_10_;
            Ctmp_11_ = TdzTiLap2TkxbP3PLUSTmdzTkxb.asDiagonal()*C32_ + TiLap2TkxbP2Tky.asDiagonal()*C42_ + Tmky*C42_;
            fft_.back_2DFull( Ctmp_11_ );
            Ctmp_11_ = dzBx_.asDiagonal()*Ctmp_11_;
            Ctmp_12_ = TdzTkxbPLUSTmdzTiLap2TkxbP3.asDiagonal()*C12_ + Tky*C22_ + TmiLap2TkxbP2Tky.asDiagonal()*C22_;
            fft_.back_2DFull( Ctmp_12_ );
            Ctmp_12_ = dzUx_.asDiagonal()*Ctmp_12_;
            Ctmp_1_ = Ctmp_1_ + Ctmp_2_ + Ctmp_3_ + Ctmp_4_;
            fft_.for_2DFull( Ctmp_1_ );
            Ctmp_5_ = Ctmp_10_ + Ctmp_11_ + Ctmp_12_ + Ctmp_5_ + Ctmp_6_ + Ctmp_7_ + Ctmp_8_ + Ctmp_9_;
            fft_.for_2DFull( Ctmp_5_ );
            Ctmp_12_ = Ctmp_1_ + ilapFtmp_.matrix().asDiagonal()*(Ctmp_5_ + (T2Tdz/fft2Dfac_).asDiagonal()*C22_ + (Tm2TkxbTkyTq/fft2Dfac_)*C12_) + (TB0zTdz/fft2Dfac_).asDiagonal()*C32_;
            dealias(Ctmp_12_);
            Ckl_out[i].block( 0, NZ_, NZ_, NZ_) = Ctmp_12_;
            
            
            
            Ctmp_1_ = Tkxb*C33_;
            fft_.back_2DFull( Ctmp_1_ );
            Ctmp_1_ = Bx_.asDiagonal()*Ctmp_1_;
            Ctmp_2_ = Tky*C33_;
            fft_.back_2DFull( Ctmp_2_ );
            Ctmp_2_ = By_.asDiagonal()*Ctmp_2_;
            Ctmp_3_ = Tmkxb*C13_;
            fft_.back_2DFull( Ctmp_3_ );
            Ctmp_3_ = Ux_.asDiagonal()*Ctmp_3_;
            Ctmp_4_ = Tmky*C13_;
            fft_.back_2DFull( Ctmp_4_ );
            Ctmp_4_ = Uy_.asDiagonal()*Ctmp_4_;
            Ctmp_5_ = T2TdzP2TiLap2Tkxb.asDiagonal()*C13_ + T2TdzTiLap2Tky.asDiagonal()*C23_;
            fft_.back_2DFull( Ctmp_5_ );
            Ctmp_5_ = dzdzUx_.asDiagonal()*Ctmp_5_;
            Ctmp_6_ = T2TdzTiLap2TkxbP2Tky.asDiagonal()*C33_ + T2TiLap2TkxbTkyP2.asDiagonal()*C43_;
            fft_.back_2DFull( Ctmp_6_ );
            Ctmp_6_ = dzBy_.asDiagonal()*Ctmp_6_;
            Ctmp_7_ = TdzTiLap2Tkxb.asDiagonal()*C13_ + TiLap2Tky.asDiagonal()*C23_;
            fft_.back_2DFull( Ctmp_7_ );
            Ctmp_7_ = dzdzdzUx_.asDiagonal()*Ctmp_7_;
            Ctmp_8_ = Tm2TdzP2TiLap2Tkxb.asDiagonal()*C33_ + Tm2TdzTiLap2Tky.asDiagonal()*C43_;
            fft_.back_2DFull( Ctmp_8_ );
            Ctmp_8_ = dzdzBx_.asDiagonal()*Ctmp_8_;
            Ctmp_9_ = Tm2TdzTiLap2TkxbP2Tky.asDiagonal()*C13_ + Tm2TiLap2TkxbTkyP2.asDiagonal()*C23_;
            fft_.back_2DFull( Ctmp_9_ );
            Ctmp_9_ = dzUy_.asDiagonal()*Ctmp_9_;
            Ctmp_10_ = TmdzTiLap2Tkxb.asDiagonal()*C33_ + TmiLap2Tky.asDiagonal()*C43_;
            fft_.back_2DFull( Ctmp_10_ );
            Ctmp_10_ = dzdzdzBx_.asDiagonal()*Ctmp_10_;
            Ctmp_11_ = TdzTiLap2TkxbP3PLUSTmdzTkxb.asDiagonal()*C33_ + TiLap2TkxbP2Tky.asDiagonal()*C43_ + Tmky*C43_;
            fft_.back_2DFull( Ctmp_11_ );
            Ctmp_11_ = dzBx_.asDiagonal()*Ctmp_11_;
            Ctmp_12_ = TdzTkxbPLUSTmdzTiLap2TkxbP3.asDiagonal()*C13_ + Tky*C23_ + TmiLap2TkxbP2Tky.asDiagonal()*C23_;
            fft_.back_2DFull( Ctmp_12_ );
            Ctmp_12_ = dzUx_.asDiagonal()*Ctmp_12_;
            Ctmp_1_ = Ctmp_1_ + Ctmp_2_ + Ctmp_3_ + Ctmp_4_;
            fft_.for_2DFull( Ctmp_1_ );
            Ctmp_5_ = Ctmp_10_ + Ctmp_11_ + Ctmp_12_ + Ctmp_5_ + Ctmp_6_ + Ctmp_7_ + Ctmp_8_ + Ctmp_9_;
            fft_.for_2DFull( Ctmp_5_ );
            Ctmp_12_ = Ctmp_1_ + ilapFtmp_.matrix().asDiagonal()*(Ctmp_5_ + (T2Tdz/fft2Dfac_).asDiagonal()*C23_ + (Tm2TkxbTkyTq/fft2Dfac_)*C13_) + (TB0zTdz/fft2Dfac_).asDiagonal()*C33_;
            dealias(Ctmp_12_);
            Ckl_out[i].block( 0, 2*NZ_, NZ_, NZ_) = Ctmp_12_;
            
            
            
            Ctmp_1_ = Tkxb*C34_;
            fft_.back_2DFull( Ctmp_1_ );
            Ctmp_1_ = Bx_.asDiagonal()*Ctmp_1_;
            Ctmp_2_ = Tky*C34_;
            fft_.back_2DFull( Ctmp_2_ );
            Ctmp_2_ = By_.asDiagonal()*Ctmp_2_;
            Ctmp_3_ = Tmkxb*C14_;
            fft_.back_2DFull( Ctmp_3_ );
            Ctmp_3_ = Ux_.asDiagonal()*Ctmp_3_;
            Ctmp_4_ = Tmky*C14_;
            fft_.back_2DFull( Ctmp_4_ );
            Ctmp_4_ = Uy_.asDiagonal()*Ctmp_4_;
            Ctmp_5_ = T2TdzP2TiLap2Tkxb.asDiagonal()*C14_ + T2TdzTiLap2Tky.asDiagonal()*C24_;
            fft_.back_2DFull( Ctmp_5_ );
            Ctmp_5_ = dzdzUx_.asDiagonal()*Ctmp_5_;
            Ctmp_6_ = T2TdzTiLap2TkxbP2Tky.asDiagonal()*C34_ + T2TiLap2TkxbTkyP2.asDiagonal()*C44_;
            fft_.back_2DFull( Ctmp_6_ );
            Ctmp_6_ = dzBy_.asDiagonal()*Ctmp_6_;
            Ctmp_7_ = TdzTiLap2Tkxb.asDiagonal()*C14_ + TiLap2Tky.asDiagonal()*C24_;
            fft_.back_2DFull( Ctmp_7_ );
            Ctmp_7_ = dzdzdzUx_.asDiagonal()*Ctmp_7_;
            Ctmp_8_ = Tm2TdzP2TiLap2Tkxb.asDiagonal()*C34_ + Tm2TdzTiLap2Tky.asDiagonal()*C44_;
            fft_.back_2DFull( Ctmp_8_ );
            Ctmp_8_ = dzdzBx_.asDiagonal()*Ctmp_8_;
            Ctmp_9_ = Tm2TdzTiLap2TkxbP2Tky.asDiagonal()*C14_ + Tm2TiLap2TkxbTkyP2.asDiagonal()*C24_;
            fft_.back_2DFull( Ctmp_9_ );
            Ctmp_9_ = dzUy_.asDiagonal()*Ctmp_9_;
            Ctmp_10_ = TmdzTiLap2Tkxb.asDiagonal()*C34_ + TmiLap2Tky.asDiagonal()*C44_;
            fft_.back_2DFull( Ctmp_10_ );
            Ctmp_10_ = dzdzdzBx_.asDiagonal()*Ctmp_10_;
            Ctmp_11_ = TdzTiLap2TkxbP3PLUSTmdzTkxb.asDiagonal()*C34_ + TiLap2TkxbP2Tky.asDiagonal()*C44_ + Tmky*C44_;
            fft_.back_2DFull( Ctmp_11_ );
            Ctmp_11_ = dzBx_.asDiagonal()*Ctmp_11_;
            Ctmp_12_ = TdzTkxbPLUSTmdzTiLap2TkxbP3.asDiagonal()*C14_ + Tky*C24_ + TmiLap2TkxbP2Tky.asDiagonal()*C24_;
            fft_.back_2DFull( Ctmp_12_ );
            Ctmp_12_ = dzUx_.asDiagonal()*Ctmp_12_;
            Ctmp_1_ = Ctmp_1_ + Ctmp_2_ + Ctmp_3_ + Ctmp_4_;
            fft_.for_2DFull( Ctmp_1_ );
            Ctmp_5_ = Ctmp_10_ + Ctmp_11_ + Ctmp_12_ + Ctmp_5_ + Ctmp_6_ + Ctmp_7_ + Ctmp_8_ + Ctmp_9_;
            fft_.for_2DFull( Ctmp_5_ );
            Ctmp_12_ = Ctmp_1_ + ilapFtmp_.matrix().asDiagonal()*(Ctmp_5_ + (T2Tdz/fft2Dfac_).asDiagonal()*C24_ + (Tm2TkxbTkyTq/fft2Dfac_)*C14_) + (TB0zTdz/fft2Dfac_).asDiagonal()*C34_;
            dealias(Ctmp_12_);
            Ckl_out[i].block( 0, 3*NZ_, NZ_, NZ_) = Ctmp_12_;
            
            
            
            Ctmp_1_ = Tkxb*C41_;
            fft_.back_2DFull( Ctmp_1_ );
            Ctmp_1_ = Bx_.asDiagonal()*Ctmp_1_;
            Ctmp_2_ = Tky*C41_;
            fft_.back_2DFull( Ctmp_2_ );
            Ctmp_2_ = By_.asDiagonal()*Ctmp_2_;
            Ctmp_3_ = Tmkxb*C21_;
            fft_.back_2DFull( Ctmp_3_ );
            Ctmp_3_ = Ux_.asDiagonal()*Ctmp_3_;
            Ctmp_4_ = Tmky*C21_;
            fft_.back_2DFull( Ctmp_4_ );
            Ctmp_4_ = Uy_.asDiagonal()*Ctmp_4_;
            Ctmp_5_ = TdzTiLap2Tkxb.asDiagonal()*C11_ + TiLap2Tky.asDiagonal()*C21_;
            fft_.back_2DFull( Ctmp_5_ );
            Ctmp_5_ = dzdzUy_.asDiagonal()*Ctmp_5_;
            Ctmp_6_ = TdzTiLap2Tkxb.asDiagonal()*C41_ + TmiLap2TkxbP2Tky.asDiagonal()*C31_;
            fft_.back_2DFull( Ctmp_6_ );
            Ctmp_6_ = dzBx_.asDiagonal()*Ctmp_6_;
            Ctmp_7_ = TiLap2TkxbP2Tky.asDiagonal()*C11_ + TmdzTiLap2Tkxb.asDiagonal()*C21_;
            fft_.back_2DFull( Ctmp_7_ );
            Ctmp_7_ = dzUx_.asDiagonal()*Ctmp_7_;
            Ctmp_8_ = TmdzTiLap2Tkxb.asDiagonal()*C31_ + TmiLap2Tky.asDiagonal()*C41_;
            fft_.back_2DFull( Ctmp_8_ );
            Ctmp_8_ = dzdzBy_.asDiagonal()*Ctmp_8_;
            Ctmp_9_ = Tkxb*C11_;
            fft_.back_2DFull( Ctmp_9_ );
            Ctmp_9_ = dzUy_.asDiagonal()*Ctmp_9_;
            Ctmp_10_ = Tmkxb*C31_;
            fft_.back_2DFull( Ctmp_10_ );
            Ctmp_10_ = dzBy_.asDiagonal()*Ctmp_10_;
            Ctmp_1_ = Ctmp_1_ + Ctmp_10_ + Ctmp_2_ + Ctmp_3_ + Ctmp_4_ + Ctmp_5_ + Ctmp_6_ + Ctmp_7_ + Ctmp_8_ + Ctmp_9_;
            fft_.for_2DFull( Ctmp_1_ );
            Ctmp_10_ = Ctmp_1_ + (TB0zTdz/fft2Dfac_).asDiagonal()*C41_ + (TdzTqPLUSTm2Tdz/fft2Dfac_).asDiagonal()*C11_;
            dealias(Ctmp_10_);
            Ckl_out[i].block( NZ_, 0, NZ_, NZ_) = Ctmp_10_;
            
            
            
            Ctmp_1_ = Tkxb*C42_;
            fft_.back_2DFull( Ctmp_1_ );
            Ctmp_1_ = Bx_.asDiagonal()*Ctmp_1_;
            Ctmp_2_ = Tky*C42_;
            fft_.back_2DFull( Ctmp_2_ );
            Ctmp_2_ = By_.asDiagonal()*Ctmp_2_;
            Ctmp_3_ = Tmkxb*C22_;
            fft_.back_2DFull( Ctmp_3_ );
            Ctmp_3_ = Ux_.asDiagonal()*Ctmp_3_;
            Ctmp_4_ = Tmky*C22_;
            fft_.back_2DFull( Ctmp_4_ );
            Ctmp_4_ = Uy_.asDiagonal()*Ctmp_4_;
            Ctmp_5_ = TdzTiLap2Tkxb.asDiagonal()*C12_ + TiLap2Tky.asDiagonal()*C22_;
            fft_.back_2DFull( Ctmp_5_ );
            Ctmp_5_ = dzdzUy_.asDiagonal()*Ctmp_5_;
            Ctmp_6_ = TdzTiLap2Tkxb.asDiagonal()*C42_ + TmiLap2TkxbP2Tky.asDiagonal()*C32_;
            fft_.back_2DFull( Ctmp_6_ );
            Ctmp_6_ = dzBx_.asDiagonal()*Ctmp_6_;
            Ctmp_7_ = TiLap2TkxbP2Tky.asDiagonal()*C12_ + TmdzTiLap2Tkxb.asDiagonal()*C22_;
            fft_.back_2DFull( Ctmp_7_ );
            Ctmp_7_ = dzUx_.asDiagonal()*Ctmp_7_;
            Ctmp_8_ = TmdzTiLap2Tkxb.asDiagonal()*C32_ + TmiLap2Tky.asDiagonal()*C42_;
            fft_.back_2DFull( Ctmp_8_ );
            Ctmp_8_ = dzdzBy_.asDiagonal()*Ctmp_8_;
            Ctmp_9_ = Tkxb*C12_;
            fft_.back_2DFull( Ctmp_9_ );
            Ctmp_9_ = dzUy_.asDiagonal()*Ctmp_9_;
            Ctmp_10_ = Tmkxb*C32_;
            fft_.back_2DFull( Ctmp_10_ );
            Ctmp_10_ = dzBy_.asDiagonal()*Ctmp_10_;
            Ctmp_1_ = Ctmp_1_ + Ctmp_10_ + Ctmp_2_ + Ctmp_3_ + Ctmp_4_ + Ctmp_5_ + Ctmp_6_ + Ctmp_7_ + Ctmp_8_ + Ctmp_9_;
            fft_.for_2DFull( Ctmp_1_ );
            Ctmp_10_ = Ctmp_1_ + (TB0zTdz/fft2Dfac_).asDiagonal()*C42_ + (TdzTqPLUSTm2Tdz/fft2Dfac_).asDiagonal()*C12_;
            dealias(Ctmp_10_);
            Ckl_out[i].block( NZ_, NZ_, NZ_, NZ_) = Ctmp_10_;
            
            
            
            Ctmp_1_ = Tkxb*C43_;
            fft_.back_2DFull( Ctmp_1_ );
            Ctmp_1_ = Bx_.asDiagonal()*Ctmp_1_;
            Ctmp_2_ = Tky*C43_;
            fft_.back_2DFull( Ctmp_2_ );
            Ctmp_2_ = By_.asDiagonal()*Ctmp_2_;
            Ctmp_3_ = Tmkxb*C23_;
            fft_.back_2DFull( Ctmp_3_ );
            Ctmp_3_ = Ux_.asDiagonal()*Ctmp_3_;
            Ctmp_4_ = Tmky*C23_;
            fft_.back_2DFull( Ctmp_4_ );
            Ctmp_4_ = Uy_.asDiagonal()*Ctmp_4_;
            Ctmp_5_ = TdzTiLap2Tkxb.asDiagonal()*C13_ + TiLap2Tky.asDiagonal()*C23_;
            fft_.back_2DFull( Ctmp_5_ );
            Ctmp_5_ = dzdzUy_.asDiagonal()*Ctmp_5_;
            Ctmp_6_ = TdzTiLap2Tkxb.asDiagonal()*C43_ + TmiLap2TkxbP2Tky.asDiagonal()*C33_;
            fft_.back_2DFull( Ctmp_6_ );
            Ctmp_6_ = dzBx_.asDiagonal()*Ctmp_6_;
            Ctmp_7_ = TiLap2TkxbP2Tky.asDiagonal()*C13_ + TmdzTiLap2Tkxb.asDiagonal()*C23_;
            fft_.back_2DFull( Ctmp_7_ );
            Ctmp_7_ = dzUx_.asDiagonal()*Ctmp_7_;
            Ctmp_8_ = TmdzTiLap2Tkxb.asDiagonal()*C33_ + TmiLap2Tky.asDiagonal()*C43_;
            fft_.back_2DFull( Ctmp_8_ );
            Ctmp_8_ = dzdzBy_.asDiagonal()*Ctmp_8_;
            Ctmp_9_ = Tkxb*C13_;
            fft_.back_2DFull( Ctmp_9_ );
            Ctmp_9_ = dzUy_.asDiagonal()*Ctmp_9_;
            Ctmp_10_ = Tmkxb*C33_;
            fft_.back_2DFull( Ctmp_10_ );
            Ctmp_10_ = dzBy_.asDiagonal()*Ctmp_10_;
            Ctmp_1_ = Ctmp_1_ + Ctmp_10_ + Ctmp_2_ + Ctmp_3_ + Ctmp_4_ + Ctmp_5_ + Ctmp_6_ + Ctmp_7_ + Ctmp_8_ + Ctmp_9_;
            fft_.for_2DFull( Ctmp_1_ );
            Ctmp_10_ = Ctmp_1_ + (TB0zTdz/fft2Dfac_).asDiagonal()*C43_ + (TdzTqPLUSTm2Tdz/fft2Dfac_).asDiagonal()*C13_;
            dealias(Ctmp_10_);
            Ckl_out[i].block( NZ_, 2*NZ_, NZ_, NZ_) = Ctmp_10_;
            
            
            
            Ctmp_1_ = Tkxb*C44_;
            fft_.back_2DFull( Ctmp_1_ );
            Ctmp_1_ = Bx_.asDiagonal()*Ctmp_1_;
            Ctmp_2_ = Tky*C44_;
            fft_.back_2DFull( Ctmp_2_ );
            Ctmp_2_ = By_.asDiagonal()*Ctmp_2_;
            Ctmp_3_ = Tmkxb*C24_;
            fft_.back_2DFull( Ctmp_3_ );
            Ctmp_3_ = Ux_.asDiagonal()*Ctmp_3_;
            Ctmp_4_ = Tmky*C24_;
            fft_.back_2DFull( Ctmp_4_ );
            Ctmp_4_ = Uy_.asDiagonal()*Ctmp_4_;
            Ctmp_5_ = TdzTiLap2Tkxb.asDiagonal()*C14_ + TiLap2Tky.asDiagonal()*C24_;
            fft_.back_2DFull( Ctmp_5_ );
            Ctmp_5_ = dzdzUy_.asDiagonal()*Ctmp_5_;
            Ctmp_6_ = TdzTiLap2Tkxb.asDiagonal()*C44_ + TmiLap2TkxbP2Tky.asDiagonal()*C34_;
            fft_.back_2DFull( Ctmp_6_ );
            Ctmp_6_ = dzBx_.asDiagonal()*Ctmp_6_;
            Ctmp_7_ = TiLap2TkxbP2Tky.asDiagonal()*C14_ + TmdzTiLap2Tkxb.asDiagonal()*C24_;
            fft_.back_2DFull( Ctmp_7_ );
            Ctmp_7_ = dzUx_.asDiagonal()*Ctmp_7_;
            Ctmp_8_ = TmdzTiLap2Tkxb.asDiagonal()*C34_ + TmiLap2Tky.asDiagonal()*C44_;
            fft_.back_2DFull( Ctmp_8_ );
            Ctmp_8_ = dzdzBy_.asDiagonal()*Ctmp_8_;
            Ctmp_9_ = Tkxb*C14_;
            fft_.back_2DFull( Ctmp_9_ );
            Ctmp_9_ = dzUy_.asDiagonal()*Ctmp_9_;
            Ctmp_10_ = Tmkxb*C34_;
            fft_.back_2DFull( Ctmp_10_ );
            Ctmp_10_ = dzBy_.asDiagonal()*Ctmp_10_;
            Ctmp_1_ = Ctmp_1_ + Ctmp_10_ + Ctmp_2_ + Ctmp_3_ + Ctmp_4_ + Ctmp_5_ + Ctmp_6_ + Ctmp_7_ + Ctmp_8_ + Ctmp_9_;
            fft_.for_2DFull( Ctmp_1_ );
            Ctmp_10_ = Ctmp_1_ + (TB0zTdz/fft2Dfac_).asDiagonal()*C44_ + (TdzTqPLUSTm2Tdz/fft2Dfac_).asDiagonal()*C14_;
            dealias(Ctmp_10_);
            Ckl_out[i].block( NZ_, 3*NZ_, NZ_, NZ_) = Ctmp_10_;
            
            
            
            Ctmp_1_ = Tkxb*C11_;
            fft_.back_2DFull( Ctmp_1_ );
            Ctmp_1_ = Bx_.asDiagonal()*Ctmp_1_;
            Ctmp_2_ = Tky*C11_;
            fft_.back_2DFull( Ctmp_2_ );
            Ctmp_2_ = By_.asDiagonal()*Ctmp_2_;
            Ctmp_3_ = Tmkxb*C31_;
            fft_.back_2DFull( Ctmp_3_ );
            Ctmp_3_ = Ux_.asDiagonal()*Ctmp_3_;
            Ctmp_4_ = Tmky*C31_;
            fft_.back_2DFull( Ctmp_4_ );
            Ctmp_4_ = Uy_.asDiagonal()*Ctmp_4_;
            Ctmp_5_ = TdzTiLap2Tkxb.asDiagonal()*C11_ + TiLap2Tky.asDiagonal()*C21_;
            fft_.back_2DFull( Ctmp_5_ );
            Ctmp_5_ = dzBx_.asDiagonal()*Ctmp_5_;
            Ctmp_6_ = TmdzTiLap2Tkxb.asDiagonal()*C31_ + TmiLap2Tky.asDiagonal()*C41_;
            fft_.back_2DFull( Ctmp_6_ );
            Ctmp_6_ = dzUx_.asDiagonal()*Ctmp_6_;
            Ctmp_1_ = Ctmp_1_ + Ctmp_2_ + Ctmp_3_ + Ctmp_4_ + Ctmp_5_ + Ctmp_6_;
            fft_.for_2DFull( Ctmp_1_ );
            Ctmp_6_ = Ctmp_1_ + (TB0zTdz/fft2Dfac_).asDiagonal()*C11_;
            dealias(Ctmp_6_);
            Ckl_out[i].block( 2*NZ_, 0, NZ_, NZ_) = Ctmp_6_;
            
            
            
            Ctmp_1_ = Tkxb*C12_;
            fft_.back_2DFull( Ctmp_1_ );
            Ctmp_1_ = Bx_.asDiagonal()*Ctmp_1_;
            Ctmp_2_ = Tky*C12_;
            fft_.back_2DFull( Ctmp_2_ );
            Ctmp_2_ = By_.asDiagonal()*Ctmp_2_;
            Ctmp_3_ = Tmkxb*C32_;
            fft_.back_2DFull( Ctmp_3_ );
            Ctmp_3_ = Ux_.asDiagonal()*Ctmp_3_;
            Ctmp_4_ = Tmky*C32_;
            fft_.back_2DFull( Ctmp_4_ );
            Ctmp_4_ = Uy_.asDiagonal()*Ctmp_4_;
            Ctmp_5_ = TdzTiLap2Tkxb.asDiagonal()*C12_ + TiLap2Tky.asDiagonal()*C22_;
            fft_.back_2DFull( Ctmp_5_ );
            Ctmp_5_ = dzBx_.asDiagonal()*Ctmp_5_;
            Ctmp_6_ = TmdzTiLap2Tkxb.asDiagonal()*C32_ + TmiLap2Tky.asDiagonal()*C42_;
            fft_.back_2DFull( Ctmp_6_ );
            Ctmp_6_ = dzUx_.asDiagonal()*Ctmp_6_;
            Ctmp_1_ = Ctmp_1_ + Ctmp_2_ + Ctmp_3_ + Ctmp_4_ + Ctmp_5_ + Ctmp_6_;
            fft_.for_2DFull( Ctmp_1_ );
            Ctmp_6_ = Ctmp_1_ + (TB0zTdz/fft2Dfac_).asDiagonal()*C12_;
            dealias(Ctmp_6_);
            Ckl_out[i].block( 2*NZ_, NZ_, NZ_, NZ_) = Ctmp_6_;
            
            
            
            Ctmp_1_ = Tkxb*C13_;
            fft_.back_2DFull( Ctmp_1_ );
            Ctmp_1_ = Bx_.asDiagonal()*Ctmp_1_;
            Ctmp_2_ = Tky*C13_;
            fft_.back_2DFull( Ctmp_2_ );
            Ctmp_2_ = By_.asDiagonal()*Ctmp_2_;
            Ctmp_3_ = Tmkxb*C33_;
            fft_.back_2DFull( Ctmp_3_ );
            Ctmp_3_ = Ux_.asDiagonal()*Ctmp_3_;
            Ctmp_4_ = Tmky*C33_;
            fft_.back_2DFull( Ctmp_4_ );
            Ctmp_4_ = Uy_.asDiagonal()*Ctmp_4_;
            Ctmp_5_ = TdzTiLap2Tkxb.asDiagonal()*C13_ + TiLap2Tky.asDiagonal()*C23_;
            fft_.back_2DFull( Ctmp_5_ );
            Ctmp_5_ = dzBx_.asDiagonal()*Ctmp_5_;
            Ctmp_6_ = TmdzTiLap2Tkxb.asDiagonal()*C33_ + TmiLap2Tky.asDiagonal()*C43_;
            fft_.back_2DFull( Ctmp_6_ );
            Ctmp_6_ = dzUx_.asDiagonal()*Ctmp_6_;
            Ctmp_1_ = Ctmp_1_ + Ctmp_2_ + Ctmp_3_ + Ctmp_4_ + Ctmp_5_ + Ctmp_6_;
            fft_.for_2DFull( Ctmp_1_ );
            Ctmp_6_ = Ctmp_1_ + (TB0zTdz/fft2Dfac_).asDiagonal()*C13_;
            dealias(Ctmp_6_);
            Ckl_out[i].block( 2*NZ_, 2*NZ_, NZ_, NZ_) = Ctmp_6_;
            
            
            
            Ctmp_1_ = Tkxb*C14_;
            fft_.back_2DFull( Ctmp_1_ );
            Ctmp_1_ = Bx_.asDiagonal()*Ctmp_1_;
            Ctmp_2_ = Tky*C14_;
            fft_.back_2DFull( Ctmp_2_ );
            Ctmp_2_ = By_.asDiagonal()*Ctmp_2_;
            Ctmp_3_ = Tmkxb*C34_;
            fft_.back_2DFull( Ctmp_3_ );
            Ctmp_3_ = Ux_.asDiagonal()*Ctmp_3_;
            Ctmp_4_ = Tmky*C34_;
            fft_.back_2DFull( Ctmp_4_ );
            Ctmp_4_ = Uy_.asDiagonal()*Ctmp_4_;
            Ctmp_5_ = TdzTiLap2Tkxb.asDiagonal()*C14_ + TiLap2Tky.asDiagonal()*C24_;
            fft_.back_2DFull( Ctmp_5_ );
            Ctmp_5_ = dzBx_.asDiagonal()*Ctmp_5_;
            Ctmp_6_ = TmdzTiLap2Tkxb.asDiagonal()*C34_ + TmiLap2Tky.asDiagonal()*C44_;
            fft_.back_2DFull( Ctmp_6_ );
            Ctmp_6_ = dzUx_.asDiagonal()*Ctmp_6_;
            Ctmp_1_ = Ctmp_1_ + Ctmp_2_ + Ctmp_3_ + Ctmp_4_ + Ctmp_5_ + Ctmp_6_;
            fft_.for_2DFull( Ctmp_1_ );
            Ctmp_6_ = Ctmp_1_ + (TB0zTdz/fft2Dfac_).asDiagonal()*C14_;
            dealias(Ctmp_6_);
            Ckl_out[i].block( 2*NZ_, 3*NZ_, NZ_, NZ_) = Ctmp_6_;
            
            
            
            Ctmp_1_ = Tkxb*C21_;
            fft_.back_2DFull( Ctmp_1_ );
            Ctmp_1_ = Bx_.asDiagonal()*Ctmp_1_;
            Ctmp_2_ = Tky*C21_;
            fft_.back_2DFull( Ctmp_2_ );
            Ctmp_2_ = By_.asDiagonal()*Ctmp_2_;
            Ctmp_3_ = Tmkxb*C41_;
            fft_.back_2DFull( Ctmp_3_ );
            Ctmp_3_ = Ux_.asDiagonal()*Ctmp_3_;
            Ctmp_4_ = Tmky*C41_;
            fft_.back_2DFull( Ctmp_4_ );
            Ctmp_4_ = Uy_.asDiagonal()*Ctmp_4_;
            Ctmp_5_ = T2TdzTiLap2Tky.asDiagonal()*C21_ + TdzP2TiLap2TkxbPLUSTmiLap2TkxbTkyP2.asDiagonal()*C11_;
            fft_.back_2DFull( Ctmp_5_ );
            Ctmp_5_ = dzBy_.asDiagonal()*Ctmp_5_;
            Ctmp_6_ = TdzTiLap2Tkxb.asDiagonal()*C11_ + TiLap2Tky.asDiagonal()*C21_;
            fft_.back_2DFull( Ctmp_6_ );
            Ctmp_6_ = dzdzBy_.asDiagonal()*Ctmp_6_;
            Ctmp_7_ = TdzTiLap2Tkxb.asDiagonal()*C21_ + TmiLap2TkxbP2Tky.asDiagonal()*C11_;
            fft_.back_2DFull( Ctmp_7_ );
            Ctmp_7_ = dzBx_.asDiagonal()*Ctmp_7_;
            Ctmp_8_ = TiLap2TkxbP2Tky.asDiagonal()*C31_ + TmdzTiLap2Tkxb.asDiagonal()*C41_;
            fft_.back_2DFull( Ctmp_8_ );
            Ctmp_8_ = dzUx_.asDiagonal()*Ctmp_8_;
            Ctmp_9_ = TiLap2TkxbTkyP2PLUSTmdzP2TiLap2Tkxb.asDiagonal()*C31_ + Tm2TdzTiLap2Tky.asDiagonal()*C41_;
            fft_.back_2DFull( Ctmp_9_ );
            Ctmp_9_ = dzUy_.asDiagonal()*Ctmp_9_;
            Ctmp_10_ = TmdzTiLap2Tkxb.asDiagonal()*C31_ + TmiLap2Tky.asDiagonal()*C41_;
            fft_.back_2DFull( Ctmp_10_ );
            Ctmp_10_ = dzdzUy_.asDiagonal()*Ctmp_10_;
            Ctmp_1_ = Ctmp_1_ + Ctmp_10_ + Ctmp_2_ + Ctmp_3_ + Ctmp_4_ + Ctmp_5_ + Ctmp_6_ + Ctmp_7_ + Ctmp_8_ + Ctmp_9_;
            fft_.for_2DFull( Ctmp_1_ );
            Ctmp_10_ = Ctmp_1_ + (TB0zTdz/fft2Dfac_).asDiagonal()*C21_ + (TmdzTq/fft2Dfac_).asDiagonal()*C31_;
            dealias(Ctmp_10_);
            Ckl_out[i].block( 3*NZ_, 0, NZ_, NZ_) = Ctmp_10_;
            
            
            
            Ctmp_1_ = Tkxb*C22_;
            fft_.back_2DFull( Ctmp_1_ );
            Ctmp_1_ = Bx_.asDiagonal()*Ctmp_1_;
            Ctmp_2_ = Tky*C22_;
            fft_.back_2DFull( Ctmp_2_ );
            Ctmp_2_ = By_.asDiagonal()*Ctmp_2_;
            Ctmp_3_ = Tmkxb*C42_;
            fft_.back_2DFull( Ctmp_3_ );
            Ctmp_3_ = Ux_.asDiagonal()*Ctmp_3_;
            Ctmp_4_ = Tmky*C42_;
            fft_.back_2DFull( Ctmp_4_ );
            Ctmp_4_ = Uy_.asDiagonal()*Ctmp_4_;
            Ctmp_5_ = T2TdzTiLap2Tky.asDiagonal()*C22_ + TdzP2TiLap2TkxbPLUSTmiLap2TkxbTkyP2.asDiagonal()*C12_;
            fft_.back_2DFull( Ctmp_5_ );
            Ctmp_5_ = dzBy_.asDiagonal()*Ctmp_5_;
            Ctmp_6_ = TdzTiLap2Tkxb.asDiagonal()*C12_ + TiLap2Tky.asDiagonal()*C22_;
            fft_.back_2DFull( Ctmp_6_ );
            Ctmp_6_ = dzdzBy_.asDiagonal()*Ctmp_6_;
            Ctmp_7_ = TdzTiLap2Tkxb.asDiagonal()*C22_ + TmiLap2TkxbP2Tky.asDiagonal()*C12_;
            fft_.back_2DFull( Ctmp_7_ );
            Ctmp_7_ = dzBx_.asDiagonal()*Ctmp_7_;
            Ctmp_8_ = TiLap2TkxbP2Tky.asDiagonal()*C32_ + TmdzTiLap2Tkxb.asDiagonal()*C42_;
            fft_.back_2DFull( Ctmp_8_ );
            Ctmp_8_ = dzUx_.asDiagonal()*Ctmp_8_;
            Ctmp_9_ = TiLap2TkxbTkyP2PLUSTmdzP2TiLap2Tkxb.asDiagonal()*C32_ + Tm2TdzTiLap2Tky.asDiagonal()*C42_;
            fft_.back_2DFull( Ctmp_9_ );
            Ctmp_9_ = dzUy_.asDiagonal()*Ctmp_9_;
            Ctmp_10_ = TmdzTiLap2Tkxb.asDiagonal()*C32_ + TmiLap2Tky.asDiagonal()*C42_;
            fft_.back_2DFull( Ctmp_10_ );
            Ctmp_10_ = dzdzUy_.asDiagonal()*Ctmp_10_;
            Ctmp_1_ = Ctmp_1_ + Ctmp_10_ + Ctmp_2_ + Ctmp_3_ + Ctmp_4_ + Ctmp_5_ + Ctmp_6_ + Ctmp_7_ + Ctmp_8_ + Ctmp_9_;
            fft_.for_2DFull( Ctmp_1_ );
            Ctmp_10_ = Ctmp_1_ + (TB0zTdz/fft2Dfac_).asDiagonal()*C22_ + (TmdzTq/fft2Dfac_).asDiagonal()*C32_;
            dealias(Ctmp_10_);
            Ckl_out[i].block( 3*NZ_, NZ_, NZ_, NZ_) = Ctmp_10_;
            
            
            
            Ctmp_1_ = Tkxb*C23_;
            fft_.back_2DFull( Ctmp_1_ );
            Ctmp_1_ = Bx_.asDiagonal()*Ctmp_1_;
            Ctmp_2_ = Tky*C23_;
            fft_.back_2DFull( Ctmp_2_ );
            Ctmp_2_ = By_.asDiagonal()*Ctmp_2_;
            Ctmp_3_ = Tmkxb*C43_;
            fft_.back_2DFull( Ctmp_3_ );
            Ctmp_3_ = Ux_.asDiagonal()*Ctmp_3_;
            Ctmp_4_ = Tmky*C43_;
            fft_.back_2DFull( Ctmp_4_ );
            Ctmp_4_ = Uy_.asDiagonal()*Ctmp_4_;
            Ctmp_5_ = T2TdzTiLap2Tky.asDiagonal()*C23_ + TdzP2TiLap2TkxbPLUSTmiLap2TkxbTkyP2.asDiagonal()*C13_;
            fft_.back_2DFull( Ctmp_5_ );
            Ctmp_5_ = dzBy_.asDiagonal()*Ctmp_5_;
            Ctmp_6_ = TdzTiLap2Tkxb.asDiagonal()*C13_ + TiLap2Tky.asDiagonal()*C23_;
            fft_.back_2DFull( Ctmp_6_ );
            Ctmp_6_ = dzdzBy_.asDiagonal()*Ctmp_6_;
            Ctmp_7_ = TdzTiLap2Tkxb.asDiagonal()*C23_ + TmiLap2TkxbP2Tky.asDiagonal()*C13_;
            fft_.back_2DFull( Ctmp_7_ );
            Ctmp_7_ = dzBx_.asDiagonal()*Ctmp_7_;
            Ctmp_8_ = TiLap2TkxbP2Tky.asDiagonal()*C33_ + TmdzTiLap2Tkxb.asDiagonal()*C43_;
            fft_.back_2DFull( Ctmp_8_ );
            Ctmp_8_ = dzUx_.asDiagonal()*Ctmp_8_;
            Ctmp_9_ = TiLap2TkxbTkyP2PLUSTmdzP2TiLap2Tkxb.asDiagonal()*C33_ + Tm2TdzTiLap2Tky.asDiagonal()*C43_;
            fft_.back_2DFull( Ctmp_9_ );
            Ctmp_9_ = dzUy_.asDiagonal()*Ctmp_9_;
            Ctmp_10_ = TmdzTiLap2Tkxb.asDiagonal()*C33_ + TmiLap2Tky.asDiagonal()*C43_;
            fft_.back_2DFull( Ctmp_10_ );
            Ctmp_10_ = dzdzUy_.asDiagonal()*Ctmp_10_;
            Ctmp_1_ = Ctmp_1_ + Ctmp_10_ + Ctmp_2_ + Ctmp_3_ + Ctmp_4_ + Ctmp_5_ + Ctmp_6_ + Ctmp_7_ + Ctmp_8_ + Ctmp_9_;
            fft_.for_2DFull( Ctmp_1_ );
            Ctmp_10_ = Ctmp_1_ + (TB0zTdz/fft2Dfac_).asDiagonal()*C23_ + (TmdzTq/fft2Dfac_).asDiagonal()*C33_;
            dealias(Ctmp_10_);
            Ckl_out[i].block( 3*NZ_, 2*NZ_, NZ_, NZ_) = Ctmp_10_;
            
            
            
            Ctmp_1_ = Tkxb*C24_;
            fft_.back_2DFull( Ctmp_1_ );
            Ctmp_1_ = Bx_.asDiagonal()*Ctmp_1_;
            Ctmp_2_ = Tky*C24_;
            fft_.back_2DFull( Ctmp_2_ );
            Ctmp_2_ = By_.asDiagonal()*Ctmp_2_;
            Ctmp_3_ = Tmkxb*C44_;
            fft_.back_2DFull( Ctmp_3_ );
            Ctmp_3_ = Ux_.asDiagonal()*Ctmp_3_;
            Ctmp_4_ = Tmky*C44_;
            fft_.back_2DFull( Ctmp_4_ );
            Ctmp_4_ = Uy_.asDiagonal()*Ctmp_4_;
            Ctmp_5_ = T2TdzTiLap2Tky.asDiagonal()*C24_ + TdzP2TiLap2TkxbPLUSTmiLap2TkxbTkyP2.asDiagonal()*C14_;
            fft_.back_2DFull( Ctmp_5_ );
            Ctmp_5_ = dzBy_.asDiagonal()*Ctmp_5_;
            Ctmp_6_ = TdzTiLap2Tkxb.asDiagonal()*C14_ + TiLap2Tky.asDiagonal()*C24_;
            fft_.back_2DFull( Ctmp_6_ );
            Ctmp_6_ = dzdzBy_.asDiagonal()*Ctmp_6_;
            Ctmp_7_ = TdzTiLap2Tkxb.asDiagonal()*C24_ + TmiLap2TkxbP2Tky.asDiagonal()*C14_;
            fft_.back_2DFull( Ctmp_7_ );
            Ctmp_7_ = dzBx_.asDiagonal()*Ctmp_7_;
            Ctmp_8_ = TiLap2TkxbP2Tky.asDiagonal()*C34_ + TmdzTiLap2Tkxb.asDiagonal()*C44_;
            fft_.back_2DFull( Ctmp_8_ );
            Ctmp_8_ = dzUx_.asDiagonal()*Ctmp_8_;
            Ctmp_9_ = TiLap2TkxbTkyP2PLUSTmdzP2TiLap2Tkxb.asDiagonal()*C34_ + Tm2TdzTiLap2Tky.asDiagonal()*C44_;
            fft_.back_2DFull( Ctmp_9_ );
            Ctmp_9_ = dzUy_.asDiagonal()*Ctmp_9_;
            Ctmp_10_ = TmdzTiLap2Tkxb.asDiagonal()*C34_ + TmiLap2Tky.asDiagonal()*C44_;
            fft_.back_2DFull( Ctmp_10_ );
            Ctmp_10_ = dzdzUy_.asDiagonal()*Ctmp_10_;
            Ctmp_1_ = Ctmp_1_ + Ctmp_10_ + Ctmp_2_ + Ctmp_3_ + Ctmp_4_ + Ctmp_5_ + Ctmp_6_ + Ctmp_7_ + Ctmp_8_ + Ctmp_9_;
            fft_.for_2DFull( Ctmp_1_ );
            Ctmp_10_ = Ctmp_1_ + (TB0zTdz/fft2Dfac_).asDiagonal()*C24_ + (TmdzTq/fft2Dfac_).asDiagonal()*C34_;
            dealias(Ctmp_10_);
            Ckl_out[i].block( 3*NZ_, 3*NZ_, NZ_, NZ_) = Ctmp_10_;
            ///     AUTOMATICALLY GENERATED EQUATIONS - END       ///
            ////////////////////////////////////////////////////////
        }
        // Legacy code used to start QL evolution at a specific time. 
//        else // THIS ASSUMES QL_YN_=0 IMPLIES ZERO MEAN FIELD!!!
//        {
//            // Evolves C in zero mean field - much faster when possible!
//            Coperator_with_zero_mean_fields(1, C11_, C21_, C31_, Ckl_out[i]);
//            Coperator_with_zero_mean_fields(2, C12_, C22_, C32_, Ckl_out[i]);
//            Coperator_with_zero_mean_fields(3, C13_, C23_, C33_, Ckl_out[i]);
//            Coperator_with_zero_mean_fields(4, C14_, C24_, C34_, Ckl_out[i]);
//
//        }
        ///     AUTOMATICALLY GENERATED EQUATIONS - END       ///
        ////////////////////////////////////////////////////////
        
        
    
        // Add to adjoint
        Ckl_out[i] += Ckl_out[i].adjoint().eval();
        
        // Add driving noise
        Ckl_out[i] += Qkl_tmp_.cast<dcmplx>().matrix().asDiagonal();
        
        /////                                ////
        /////////////////////////////////////////
        

        //Since I almost always want the Reynolds stress in a linear calculation better to still calculate and only change the MF feedback when QL_YN_=0
        ///////////////////////////////////////////////////////////
        //////                                               //////
        //////             REYNOLDS STRESSES                 //////
        //////                                               //////
        
        {
        // mult_fac for sum - since no -ve ky values are used, count everything but ky=0 twice
        double mult_fac = 2.0;
        if (kytmp_== 0.0 )
            mult_fac = 1.0; // Only count ky=0 mode once
        
        double ftfac = mult_fac*fft2Dfac_; // NZ^2 factor for ifft is included here - more efficient since scalar

        //Assign automatically generated variables in equations
        // Vector variable definition Reynolds stress (automatic)
        rey_TdzTiLap2TkxbTky = (ilap2tmp_.cast<dcmplx>()*kz_).matrix()*(kxctmp_*kyctmp_);
        rey_TiLap2TkxbTkyP2 = ilap2tmp_.cast<dcmplx>().matrix()*(kxctmp_*pow(kyctmp_,2));
        rey_TiLap2TkxbP2Tky = ilap2tmp_.cast<dcmplx>().matrix()*(kyctmp_*pow(kxctmp_,2));
        rey_TiLap2TkxbTky = ilap2tmp_.cast<dcmplx>().matrix()*(kxctmp_*kyctmp_);
        rey_TdzTiLap2Tkxb = (ilap2tmp_.cast<dcmplx>()*kz_).matrix()*kxctmp_;
        rey_TdzTiLap2Tky = (ilap2tmp_.cast<dcmplx>()*kz_).matrix()*kyctmp_;
        rey_TdzP2TiLap2 = (ilap2tmp_.cast<dcmplx>()*kz_.pow(2)).matrix();
        rey_TiLap2Tky = ilap2tmp_.cast<dcmplx>().matrix()*kyctmp_;
        rey_TdzTiLap2 = (ilap2tmp_.cast<dcmplx>()*kz_).matrix();
        rey_Tdz = kz_.matrix();

        
        // bz*ux - uz*bx
        // Bx reynolds stress
        reynolds_mat_tmp_ = C31_*rey_TdzTiLap2Tkxb.asDiagonal()-C32_*rey_TiLap2Tky.asDiagonal()-rey_TdzTiLap2Tkxb.asDiagonal()*C31_-rey_TiLap2Tky.asDiagonal()*C41_;
        // fft(fft( a )')'
        fft_.back_2DFull(reynolds_mat_tmp_);
        // Keep a running sum
        bzux_m_uzbx_d_ += ftfac*(reynolds_mat_tmp_.diagonal().array().real());
        
        
        // bz*uy - uz*by
        // By reynolds stress
        reynolds_mat_tmp_ = rey_TdzTiLap2.asDiagonal()*C41_*rey_TdzTiLap2Tkxb.asDiagonal()+rey_TdzTiLap2Tkxb.asDiagonal()*C31_*rey_TiLap2TkxbTky.asDiagonal()+rey_TdzTiLap2Tkxb.asDiagonal()*C32_*rey_TdzTiLap2.asDiagonal()+rey_TiLap2TkxbTky.asDiagonal()*C32_*rey_TiLap2Tky.asDiagonal()+rey_TiLap2Tky.asDiagonal()*C41_*rey_TiLap2TkxbTky.asDiagonal()+rey_TiLap2Tky.asDiagonal()*C42_*rey_TdzTiLap2.asDiagonal()-rey_TdzTiLap2.asDiagonal()*C42_*rey_TiLap2Tky.asDiagonal()-rey_TiLap2TkxbTky.asDiagonal()*C31_*rey_TdzTiLap2Tkxb.asDiagonal();
        // fft(fft( a )')'
        fft_.back_2DFull(reynolds_mat_tmp_);
        // Keep a running sum
        bzuy_m_uzby_d_ += ftfac*(reynolds_mat_tmp_.diagonal().array().real());
        
        // Ux reynolds stress
        reynolds_mat_tmp_ = -kxctmp_*C11_ + kxctmp_*C33_ + kyctmp_*C11_*rey_TiLap2TkxbTky.asDiagonal() + kyctmp_*C12_*rey_TdzTiLap2.asDiagonal() - kyctmp_*C33_*rey_TiLap2TkxbTky.asDiagonal() - kyctmp_*C34_*rey_TdzTiLap2.asDiagonal() + rey_Tdz.asDiagonal()*C11_*rey_TdzTiLap2Tkxb.asDiagonal() - rey_Tdz.asDiagonal()*C12_*rey_TiLap2Tky.asDiagonal() - rey_Tdz.asDiagonal()*C33_*rey_TdzTiLap2Tkxb.asDiagonal() + rey_Tdz.asDiagonal()*C34_*rey_TiLap2Tky.asDiagonal();
        // fft(fft( a )')'
        fft_.back_2DFull(reynolds_mat_tmp_);
        // Keep a running sum
        ux_rey_stress_d_ += ftfac*(reynolds_mat_tmp_.diagonal().array().real());
        
        // Uy reynolds stress
        reynolds_mat_tmp_ = rey_TdzP2TiLap2.asDiagonal()*C21_*rey_TdzTiLap2Tkxb.asDiagonal() - rey_TdzP2TiLap2.asDiagonal()*C22_*rey_TiLap2Tky.asDiagonal() - rey_TdzP2TiLap2.asDiagonal()*C43_*rey_TdzTiLap2Tkxb.asDiagonal() + rey_TdzP2TiLap2.asDiagonal()*C44_*rey_TiLap2Tky.asDiagonal() - rey_TdzTiLap2Tkxb.asDiagonal()*C21_ + rey_TdzTiLap2Tkxb.asDiagonal()*C43_ - rey_TdzTiLap2TkxbTky.asDiagonal()*C11_*rey_TdzTiLap2Tkxb.asDiagonal() + rey_TdzTiLap2TkxbTky.asDiagonal()*C12_*rey_TiLap2Tky.asDiagonal() + rey_TdzTiLap2TkxbTky.asDiagonal()*C33_*rey_TdzTiLap2Tkxb.asDiagonal() - rey_TdzTiLap2TkxbTky.asDiagonal()*C34_*rey_TiLap2Tky.asDiagonal() + rey_TdzTiLap2Tky.asDiagonal()*C21_*rey_TiLap2TkxbTky.asDiagonal() + rey_TdzTiLap2Tky.asDiagonal()*C22_*rey_TdzTiLap2.asDiagonal() - rey_TdzTiLap2Tky.asDiagonal()*C43_*rey_TiLap2TkxbTky.asDiagonal() - rey_TdzTiLap2Tky.asDiagonal()*C44_*rey_TdzTiLap2.asDiagonal() + rey_TiLap2TkxbP2Tky.asDiagonal()*C11_ - rey_TiLap2TkxbP2Tky.asDiagonal()*C33_ - rey_TiLap2TkxbTkyP2.asDiagonal()*C11_*rey_TiLap2TkxbTky.asDiagonal() - rey_TiLap2TkxbTkyP2.asDiagonal()*C12_*rey_TdzTiLap2.asDiagonal() + rey_TiLap2TkxbTkyP2.asDiagonal()*C33_*rey_TiLap2TkxbTky.asDiagonal() + rey_TiLap2TkxbTkyP2.asDiagonal()*C34_*rey_TdzTiLap2.asDiagonal();
        // fft(fft( a )')'
        fft_.back_2DFull(reynolds_mat_tmp_);
        // Keep a running sum
        uy_rey_stress_d_ += ftfac*(reynolds_mat_tmp_.diagonal().array().real());
        }
        
        //////                                               //////
        ///////////////////////////////////////////////////////////
    
        
        
        
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
    {
    reynolds_stress_MPI_send_ << bzux_m_uzbx_d_, bzuy_m_uzby_d_, ux_rey_stress_d_, uy_rey_stress_d_;// Put into a single variable to reduce mpi calls
    
    // FOR SOME REASON MPI::IN_PLACE IS GIVING INCORRECT RESULTS!
    mpi_.SumAllReduce_double(reynolds_stress_MPI_send_.data(), reynolds_stress_MPI_receive_.data(), num_MFs()*NZ_);
    // Extract from ..._receive_ variable and convert to complex
    // (didn't seem worth the effort of fftw real transforms as time spent here (and in MF transforms) will be very minimal)
    bzux_m_uzbx_c_ = reynolds_stress_MPI_receive_.segment(0,NZ_).cast<dcmplx>();
    bzuy_m_uzby_c_ = reynolds_stress_MPI_receive_.segment(NZ_,NZ_).cast<dcmplx>();
    ux_rey_stress_c_ = reynolds_stress_MPI_receive_.segment(2*NZ_,NZ_).cast<dcmplx>();
    uy_rey_stress_c_ = reynolds_stress_MPI_receive_.segment(3*NZ_,NZ_).cast<dcmplx>();

    
    // Calculate stresses in Fourier space
    fft_.for_1D(bzux_m_uzbx_c_.data());
    fft_.for_1D(bzuy_m_uzby_c_.data());
    fft_.for_1D(ux_rey_stress_c_.data());
    fft_.for_1D(uy_rey_stress_c_.data());

    
    bzux_m_uzbx_c_ = fftFac_Reynolds_*kz_*bzux_m_uzbx_c_;
    bzuy_m_uzby_c_ = fftFac_Reynolds_*kz_*bzuy_m_uzby_c_;
    ux_rey_stress_c_ *= fftFac_Reynolds_;
    uy_rey_stress_c_ *= fftFac_Reynolds_;

    dealias(bzux_m_uzbx_c_);
    dealias(bzuy_m_uzby_c_);
    dealias(ux_rey_stress_c_);
    dealias(uy_rey_stress_c_);
    }
    //////////////////////////////////////

    
    //////////////////////////////////////
    ////   MEAN FIELDS    ////////////////
    if (QL_YN_) {
        // B update
        MFout[0] = B0z_*kz_*MFin[2] + bzux_m_uzbx_c_;
        MFout[1] = -q_*MFin[0] + B0z_*kz_*MFin[3] + bzuy_m_uzby_c_;
        // U update
        MFout[2] = 2*Omega_*MFin[3] + B0z_*kz_*MFin[0] + ux_rey_stress_c_;
        MFout[3] = (q_-2.0*Omega_)*MFin[2] + B0z_*kz_*MFin[1] + uy_rey_stress_c_;
    } else { // Still calculate Reynolds stress in linear calculation
        // B update
        MFout[1].setZero();
        MFout[0].setZero();
        // U update
        MFout[2].setZero();
        MFout[3].setZero();
        
    }
    
    
    

}


//////////////////////////////////////////////////////////
///   Initialization of the linear operators        //////

void MHD_FullUBQlin::linearOPs_Init(double t0, doubVec * linop_MF, doubVec * linop_Ckl_old) {
    // Constructs initial (t=0) linear operators, since MF operator is constant this completely defines it for entire integration
    // Mean fields (Bx, By, Ux, Uy)
    if (QL_YN_) {
        linop_MF[0] = eta_*kz2_;
        linop_MF[1] = eta_*kz2_;
        linop_MF[2] = nu_*kz2_;
        linop_MF[3] = nu_*kz2_;

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
void MHD_FullUBQlin::dealias(dcmplxMat &inMat) {
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
void MHD_FullUBQlin::dealias(dcmplxVec& vec) {
    for (int j=delaiasBnds_[0]; j<delaiasBnds_[1]+1; ++j){
        vec(j) = 0;
    }
}

//1-D version - takes a doubVec input
void MHD_FullUBQlin::dealias(doubVec& vec) {
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

void MHD_FullUBQlin::Calc_Energy_AM_Diss(TimeVariables& tv, double t, const dcmplxVec *MFin, const dcmplxMat *Cin ) {
    // Energy, angular momentum and dissipation of the solution MFin and Cin
    // TimeVariables class stores info and data about energy, AM etc.
    // t is time
    
    tv.start_timing();
    
    // OUTPUT: energy[1] and [2] contain U and B mean field energies (energy[1]=0 for this)
    // energy[3] and [4] contain u and b fluctuating energies.
    
    // Energy
    double energy_u=0, energy_b=0;
    // Angular momentum
    double AM_u = 0, AM_b = 0;
    // Dissipation
    double diss_u=0,diss_b=0;
    
    // MPI buffers - this method passes around some unecessary data but will be very minimal
    const int num_to_mpi = 6;

    
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
            lap2tmp_ = lapFtmp_*ilap2tmp_; // lap2tmp_ just a convenient storage
            Qkl_tmp_ << lap2tmp_, -ilap2tmp_,lap2tmp_, -ilap2tmp_;
            
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
            rey_TiLap2TkxbTky = (-kyctmp_*kxctmp_ )*ilap2tmp_.cast<dcmplx>().matrix();
            rey_TdzTiLap2 = (kz_*ilap2tmp_).matrix();
            
            // uy*ux - steal Reynolds stress variable for this
            bzux_m_uzbx_d_ = (rey_TiLap2TkxbTky.array()*Cin[i].block( 0, 0, NZ_, NZ_).diagonal().array()).real() + (rey_TdzTiLap2.array()*Cin[i].block(NZ_, 0,NZ_, NZ_).diagonal().array()).real();
            AM_u += mult_fac*bzux_m_uzbx_d_.sum();
            // by*bx - steal Reynolds stress variable for this
            bzux_m_uzbx_d_ = (rey_TiLap2TkxbTky.array()*Cin[i].block( 2*NZ_,2*NZ_, NZ_, NZ_).diagonal().array()).real() +   (rey_TdzTiLap2.array()*Cin[i].block(3*NZ_, 2*NZ_,NZ_, NZ_).diagonal().array()).real();
            AM_b += mult_fac*bzux_m_uzbx_d_.sum();
            //
            //////////////////////////////////////
        }
        if (tv.dissip_save_Q()){
            //////////////////////////////////////
            //     DISSIPATION
            // Use Qkl_tmp_ for Mkl to save memory
            // CHANGES ilap2tmp_ SO PUT LAST OUT OF DIAGNOSTICS!!
            ilap2tmp_ *= lapFtmp_; // lap2tmp_ just a convenient storage
            lap2tmp_ = -lapFtmp_*ilap2tmp_;
            Qkl_tmp_ << nu_*lap2tmp_, nu_*ilap2tmp_,eta_*lap2tmp_, eta_*ilap2tmp_;
            
            // Energy = trace(Mkl*Ckl), Mkl is diagonal
            Qkl_tmp_ = mult_fac*Qkl_tmp_ * Cin[i].real().diagonal().array();
            
            int num_u_b = num_fluct_/2; // Keep it clear where numbers are coming from.
            diss_u += Qkl_tmp_.head(NZ_*num_u_b).sum();
            diss_b += Qkl_tmp_.tail(NZ_*num_u_b).sum();
            //
            //////////////////////////////////////
        }
        
        
    }
    double mpi_send_buff[num_to_mpi] = {energy_u,energy_b,AM_u,AM_b,diss_u,diss_b};
    double mpi_receive_buff[num_to_mpi];

    // Put the everything on processor 0
    mpi_.SumReduce_doub(mpi_send_buff,mpi_receive_buff,num_to_mpi);
//    mpi_.SumReduce_IP_doub(&energy_u_f,1); // Is this working?

    
    // Currently, TimeVariables is set to save on root process, may want to generalize this at some point (SumReduce_doub is also)
    if (mpi_.my_n_v() == 0) {
        ////////////////////////////////////////////
        ///// All this is only on processor 0  /////
        double divfac=box_length(2)/totalN2_;
        energy_u = mpi_receive_buff[0]*divfac;
        energy_b = mpi_receive_buff[1]*divfac;
        AM_u = mpi_receive_buff[2]*divfac;
        AM_b = mpi_receive_buff[3]*divfac;
        diss_u = mpi_receive_buff[4]*divfac;
        diss_b = mpi_receive_buff[5]*divfac;
        
        
        ///////////////////////////////////////
        ///       MEAN FIELDS            //////
        // Only need to calculate on one processor
        double energy_MU=0, energy_MB=0;
        double AM_MU =0, AM_MB=0;
        double diss_MU =0, diss_MB =0;
        
        energy_MB = box_length(2)*(MFin[0].abs2().sum() + MFin[1].abs2().sum())/(NZ_*NZ_);
        energy_MU = box_length(2)*(MFin[2].abs2().sum() + MFin[3].abs2().sum())/(NZ_*NZ_);

        AM_MB = box_length(2)*(MFin[0]*MFin[1].conjugate()).real().sum()/(NZ_*NZ_);
        AM_MU = box_length(2)*(MFin[2]*MFin[3].conjugate()).real().sum()/(NZ_*NZ_);
        
        diss_MB = box_length(2)*(-eta_/(NZ_*NZ_))*((MFin[0].abs2()*kz2_).sum() + (MFin[1].abs2()*kz2_).sum() );
        diss_MU = box_length(2)*(-nu_/(NZ_*NZ_))*((MFin[2].abs2()*kz2_).sum() + (MFin[3].abs2()*kz2_).sum());
        
        ///////////////////////////////////////
        //////         OUTPUT            //////
        
        // Energy
        double* en_point = tv.current_energy();
        en_point[0] = energy_MU/2;
        en_point[1] = energy_MB/2;
        en_point[2] = energy_u/2;
        en_point[3] = energy_b/2;
        
        // Angular momentum
        double* AM_point = tv.current_AM();
        AM_point[0] = AM_MU;
        AM_point[1] = AM_MB;
        AM_point[2] = AM_u;
        AM_point[3] = AM_b;
        
        // Energy
        double* diss_point = tv.current_diss();
        diss_point[0] = diss_MU;
        diss_point[1] = diss_MB;
        diss_point[2] = diss_u;
        diss_point[3] = diss_b;
        
        if (tv.reynolds_save_Q()) {
            //////////////////////////////////////
            //   Reynolds stress
            // Saves quantities to do with the reynolds stress and dynamo
            // 1) Shear contribution: Re( -q Bx By)/|By|
            // 2) y emf: Re( bzuy_m_uzby_c_*By )/|By|
            // 3) y dissipation: eta*k^2*By
            // 4) x emf: Re( bzux_m_uzbx_c_*Bx )/|Bx|
            // 5) x dissipation: eta*k^2*Bx
            // Have assumed k0 to be lowest kz! i.e., driving largest dynamo possible in the box
            // There may be slight errors here from saving using the updated values of MFin, presumably this is a small effect, especially at high resolution (low dt) and in steady state.
            double* rey_point = tv.current_reynolds();
            rey_point[0] = -q_*(MFin[0]*MFin[1].conjugate()).real().sum()/sqrt(MFin[1].abs2().sum());
            rey_point[1] = (bzuy_m_uzby_c_*MFin[1].conjugate()).real().sum()/sqrt(MFin[1].abs2().sum());
            rey_point[2] = -eta_ *sqrt((kz2_*MFin[1]).abs2().sum());
            rey_point[3] = (bzux_m_uzbx_c_*MFin[0].conjugate()).real().sum()/sqrt(MFin[0].abs2().sum());
            rey_point[4] = -eta_ *sqrt((kz2_*MFin[0]).abs2().sum());
//            rey_point[0] = real(bzuy_m_uzby_c_(1)*conj(MFin[1](1))) ;
//            rey_point[1] = real(bzux_m_uzbx_c_(1)*conj(MFin[1](1))) ;
            //
            //////////////////////////////////////
        }

        ///// All this is only on processor 0  /////
        ////////////////////////////////////////////
    }
    
    // Save the data
    tv.Save_Data(t);
    
    tv.finish_timing();
    
}




//////////////////////////////////////////////////////
//         PRIVATE FUNCTIONS FOR INITIALIZATION     //

// Create and store arrays of lap2 to save computation of ffts
void MHD_FullUBQlin::Define_Lap2_Arrays_(void){
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


void MHD_FullUBQlin::print_noise_range_(){
    // Print out total number of driven modes
    int tot_count = 0, driv_count = 0;
    int tot_count_all = 0, driv_count_all = 0;// MPI reduced versions
    for (int i=0; i<Cdimxy(); ++i) {
        int k_i = i + index_for_k_array(); // k index
        int ind_ky = ky_index_[k_i];
         // Form Laplacians using time-dependent kx
        kytmp_ = ky_[k_i].imag();
        kxtmp_ = kx_[k_i].imag();
        lap2tmp_ = lap2_[ind_ky];
        lapFtmp_ = -kxtmp_*kxtmp_ + lap2tmp_;
        dealias(lapFtmp_);
        
        drive_condition_ = lapFtmp_ != 0.0;
        tot_count += drive_condition_.cast<int>().sum();;
        drive_condition_ = lapFtmp_> -noise_range_[1] && lapFtmp_< -noise_range_[0];
        driv_count += drive_condition_.cast<int>().sum();
    }
    mpi_.SumReduce_int(&tot_count, &tot_count_all, 1);
    mpi_.SumReduce_int(&driv_count, &driv_count_all, 1);
    std::stringstream noise_output;
    noise_output << "Driving " << driv_count_all << " out of a total of " << tot_count_all << " modes (k = " << sqrt(noise_range_[0]) << " to " << sqrt(noise_range_[1]) << ")\n";
    mpi_.print1(noise_output.str());
    
}
//////////////////////////////////////////////////////


//////////////////////////
// CFL number
double MHD_FullUBQlin::Calculate_CFL(const dcmplxVec *MFin,const dcmplxMat* Ckl)  {
    // Returns CFL/dt to calculate dt - in this case CFL/dt = kmax By + q
    
    // Mean fields have been previously calculated, so may as well use them
    double Bxmax = sqrt(Bx_.array().abs2().maxCoeff());
    double Bymax = sqrt(By_.array().abs2().maxCoeff());
    double Uxmax = sqrt(Ux_.array().abs2().maxCoeff());
    double Uymax = sqrt(Uy_.array().abs2().maxCoeff());
    // CFL
    return kmax*(Bymax + Bxmax + Uxmax + Uymax) + q_;
}



///////////////////////////////////////////////////////////
//      OPERATOR WTIH B=0 -  FOR SPEED WHEN POSSIBLE     //

// Evolves C in zero mean field - much faster when possible!
void MHD_FullUBQlin::Coperator_with_zero_mean_fields(int row, dcmplxMat& C1i, dcmplxMat& C2i, dcmplxMat& C3i, dcmplxMat& Cklout){

    Cklout.block( 0, (row-1)*NZ_, NZ_, NZ_) = ilapFtmp_.matrix().asDiagonal()*((T2Tdz/fft2Dfac_).asDiagonal()*C2i + (Tm2TkxbTkyTq/fft2Dfac_)*C1i);
    Cklout.block( NZ_, (row-1)*NZ_, NZ_, NZ_) = (TdzTqPLUSTm2Tdz/fft2Dfac_).asDiagonal()*C1i;
    Cklout.block( 2*NZ_, (row-1)*NZ_, NZ_, NZ_).setZero();
    Cklout.block( 3*NZ_, (row-1)*NZ_, NZ_, NZ_) = (TmdzTq/fft2Dfac_).asDiagonal()*C3i;
    
}




