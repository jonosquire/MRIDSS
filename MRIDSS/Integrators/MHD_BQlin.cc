//
//  MHD_BQlin.cc
//  MRIDSS
//
//  Created by Jonathan Squire on 4/25/14.
//  Copyright (c) 2014 J Squire. All rights reserved.
//

#include "MHD_BQlin.h"


// Model class for S3T/CE2 shearing box MHD model
// Basic Quasi-linear MHD
//
// Derived from Model (model.h)


// Constructor for MHD_BQlin
MHD_BQlin::MHD_BQlin(const Inputs& sp, MPIdata& mpi, fftwPlans& fft) :
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
    delaiasBnds_[1] = NZ_-NZ_/3-1;

    
    ////////////////////////////////////////////////////
    // Useful arrays to save computation
    // ifft of identity and 1/lap2 etc.
    
    Set_fft_identity_();// fft of identity
    Define_Lap2_Arrays_(); // Lap2 and ffts - Allocates data to lap2_,ilap2_,fft_ilap2_,fft_kzilap2_,fft_kz2ilap2_
    
    
    // Reynolds stresses
    bzux_m_uzbx_c_ = dcmplxVec::Zero(NZ_);
    bzuy_m_uzby_c_ = dcmplxVec::Zero(NZ_);
    bzux_m_uzbx_d_ = doubVec::Zero(NZ_);
    bzuy_m_uzby_d_ = doubVec::Zero(NZ_);
    // Receive buffers for MPI, for some reason MPI::IN_PLACE is not giving the correct result!
    bzux_m_uzbx_drec_ = doubVec::Zero(NZ_);
    bzuy_m_uzby_drec_ = doubVec::Zero(NZ_);
    
    
    ////////////////////////////////////////////////////
    //                                                //
    //               TEMPORARY ARRAYS                 //
    //                                                //
    //  These are used during evaluation of the various parts of the model
    // There should never be any need to keep their value over more than 1 time-step
    
    lapFtmp_ = doubVec( NZ_ ); //For (time-dependent) k^2
    ilapFtmp_ = doubVec( NZ_ ); // Inverse
    lap2tmp_ = doubVec( NZ_ ); // For ky^2+kz^2 - could be pre-assigned
    ilap2tmp_ = doubVec( NZ_ ); // Inverse
    // Qkl
    Qkl_tmp_ = doubVec( nz_Cfull_ );
    
    // A operator matrix - lots of submatrices
    dcmplxVec Aop_v11_, Aop_v12_, Aop_v21_, Aop_v43_;
    dcmplxMat Aop_m13_, Aop_m14_, Aop_m23_, Aop_m24_, Aop_m31_, Aop_m41_, Aop_m42_;
    // Zero - A22 A32 A33 A34 A44
    // Vectors
    Aop_v11_ = dcmplxVec(NZ_);
    Aop_v12_ = dcmplxVec(NZ_);
    Aop_v21_ = dcmplxVec(NZ_);
    Aop_v43_ = dcmplxVec(NZ_);
    // Matrices
    Aop_m13_= dcmplxMat(NZ_,NZ_);
    Aop_m14_= dcmplxMat(NZ_,NZ_);
    Aop_m23_= dcmplxMat(NZ_,NZ_);
    Aop_m24_= dcmplxMat(NZ_,NZ_);
    Aop_m31_= dcmplxMat(NZ_,NZ_);
    Aop_m41_= dcmplxMat(NZ_,NZ_);
    Aop_m42_= dcmplxMat(NZ_,NZ_);
    Aop_block_tmp_= dcmplxMat(NZ_,NZ_); // Useful temporary

    
    // Real versions of B
    rBy_tmp_ = dcmplxVec::Zero(NZ_);
    rDzBy_tmp_ = dcmplxVec::Zero(NZ_);
    rDzzBy_tmp_ = dcmplxVec::Zero(NZ_);
    
    // Reynolds stresses
    reynolds_mat_tmp_ = dcmplxMat(NZ_,NZ_);
    // converting between u, zeta... and u, uy... for the Reynolds stresses
    rey_mkxky_tmp_ = Eigen::Matrix<dcmplx, Eigen::Dynamic, 1>( NZ_ );
    rey_kz_tmp_ = Eigen::Matrix<dcmplx, Eigen::Dynamic, 1>( NZ_ );
    rey_mkxkz_tmp_ = Eigen::Matrix<dcmplx, Eigen::Dynamic, 1>( NZ_ );
    rey_mky_tmp_ = Eigen::Matrix<dcmplx, Eigen::Dynamic, 1>( NZ_ );

}


// Destructor
MHD_BQlin::~MHD_BQlin(){
    // Here I am destroying Base members from the derived destructor - is this bad?
    delete[] kx_;
    delete[] ky_;
    
    // useful fft arrays
    delete[] ky_index_;
    delete[] lap2_;
    delete[] ilap2_;
    delete[] fft_ilap2_;
    delete[] fft_kzilap2_;
    delete[] fft_kz2ilap2_;
    
}


// Main EulerFluid function, fx = f(t,x) called at ecah time-step
void MHD_BQlin::rhs(double t, double dt_lin,
                   dcmplxVec *MFin, dcmplxMat *Ckl_in,
                   dcmplxVec *MFout, dcmplxMat *Ckl_out,
                   doubVec * linop_Ckl) {
    
    // If desired, the function should work fine with Ckl_in = Ckl_out. This is
    // useful for integration schemes that do not require re-use of Ckl_in
    
    // Calculate MFs in real space - for this model, only need MF[2], ie., By
    rBy_tmp_ = MFin[1]*(1.0/NZ_); // fft back doesn't include normalization
    fft_.back_1D(rBy_tmp_.data());
    rDzBy_tmp_ = kz_*MFin[1]*(1.0/NZ_);
    fft_.back_1D(rDzBy_tmp_.data());
    rDzzBy_tmp_ = kz2_*MFin[1]*(1.0/NZ_);
    fft_.back_1D(rDzzBy_tmp_.data());
    
    // Reynolds stresses -- added to at each step in the loop
    bzux_m_uzbx_d_.setZero();
    bzuy_m_uzby_d_.setZero();
    
    
    /////////////////////////////////////
    //   ALL OF THIS MUST BE PARALLELIZED
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
        if (kytmp_==0 && kxtmp_==0) {
            Qkl_tmp_.setZero(); // Don't drive the mean components! (really could leave mean cmpts out entirely)
            ilapFtmp_(0) = 1; // Avoid infinities (lap2 already done, since stored)
        }
        else {
            lapFtmp_ = lap2tmp_/lapFtmp_; // lapFtmp_ is no longer lapF!!!
            if (kytmp_==0) { // Don't drive the ky=kz=0 component (variables not invertible)
                lapFtmp_(0)=0;
            }
            Qkl_tmp_ << lapFtmp_, lap2tmp_, lapFtmp_, lap2tmp_;
            Qkl_tmp_ = (f_noise_*f_noise_)*Qkl_tmp_.abs();
        }
        ////////////////////////////////////////
        
        
        ////////////////////////////////////////
        //                                    //
        //      MAIN WORK - FORM OPMAT        //
        //                                    //
        
        // To avoid problems, use kxctmp and kyctmp
        
        //  ROW 1
        
        // u u Component
        // -2*q*ilapF*kxt*ky
        Aop_v11_= ((-2*q_*kxctmp_*kyctmp_)*ilapFtmp_.cast<dcmplx>()).matrix().asDiagonal();
        
        // u zeta Component
        // 2*ilapF.*K.kz
        Aop_v12_ = (2*ilapFtmp_*kz_).matrix().asDiagonal();
        
        // u b Component
        // dealiasmat.*( fft(ky*realB.*idft) + 2*kxt^2*ky*ilapF*fft(realDzB.*rilap2kz) )
        Aop_m13_ = (fft_identity_.array().colwise()*(kyctmp_*rBy_tmp_)).matrix();
        fft_.for_2D_dim1(Aop_m13_.data());
        Aop_block_tmp_ = (fft_kzilap2_[ind_ky].array().colwise()*rDzBy_tmp_).matrix();
        fft_.for_2D_dim1(Aop_block_tmp_.data());
        Aop_m13_ +=
            ((2.0*kxctmp_*kxctmp_*kyctmp_)*ilapFtmp_.cast<dcmplx>()).matrix().asDiagonal()*
                Aop_block_tmp_;
        dealias(Aop_m13_.data()); // Could also do delaliasing in assignment
        // Final assignment done
        
        // u eta component
        // dealiasmat.*( 2*kxt*ky^2*ilapF*fft( realDzB.*rilap2 )  )
        Aop_m14_ = (fft_ilap2_[ind_ky].array().colwise()*rDzBy_tmp_).matrix();
        fft_.for_2D_dim1(Aop_m14_.data());
        Aop_m14_ =
            ((2.0*kxctmp_*kyctmp_*kyctmp_)*ilapFtmp_.cast<dcmplx>()).matrix().asDiagonal()*
                Aop_m14_;
        dealias(Aop_m14_.data());
        // Final assignment done
        
        
        // ROW 2
        
        // zeta u
        // (q-2)*K.kz
        Aop_v21_ = ((q_-2.0)*kz_).matrix().asDiagonal();
        
        // zeta zeta
        // 0
        
        // zeta b
        // dealiasmat.*( fft(-kxt*realDzzB.*rilap2kz-kxt*realDzB.*idft) )
        Aop_m23_ = (  fft_kzilap2_[ind_ky].array().colwise()*((-kxctmp_)*rDzzBy_tmp_) +
            fft_identity_.array().colwise()*(-kxctmp_*rDzBy_tmp_)  ).matrix();
        fft_.for_2D_dim1(Aop_m23_.data());
        dealias(Aop_m23_.data());
        // Final assignment done
        
        // zeta eta
        // dealiasmat.*( fft( ky*realB.*idft - ky*realDzzB.*rilap2 ))
        Aop_m24_ = (  fft_ilap2_[ind_ky].array().colwise()*((-kyctmp_)*rDzzBy_tmp_) +
            fft_identity_.array().colwise()*(kyctmp_*rBy_tmp_)  ).matrix();
        fft_.for_2D_dim1(Aop_m24_.data());
        dealias(Aop_m24_.data());
        // Final assignment done

        
        // ROW 3

        // b u
        // dealiasmat.*( fft(ky*realB.*idft) )
        Aop_m31_ = (  fft_identity_.array().colwise()*(kyctmp_*rBy_tmp_)  ).matrix();
        fft_.for_2D_dim1(Aop_m31_.data());
        dealias(Aop_m31_.data());
        // Final assignment done
        
        // b zeta
        // 0
        
        // b b
        // 0 (except linear)
        
        // b eta
        //0
        
        
        // ROW 4
        
        // eta u
        // dealiasmat.*(fft( 2*kxt*realDzB.*ifft(ilap2.*K.kz.^2) + kxt*realDzzB.*rilap2kz - kxt*realDzB.*idft  ))
        Aop_m41_ =
            (  fft_kz2ilap2_[ind_ky].array().colwise()*((2.0*kxctmp_)*rDzBy_tmp_) +
            fft_kzilap2_[ind_ky].array().colwise()*(kxctmp_*rDzzBy_tmp_)  +
            fft_identity_.array().colwise()*((-kxctmp_)*rDzBy_tmp_)  ).matrix();
        fft_.for_2D_dim1(Aop_m41_.data());
        dealias(Aop_m41_.data());
        
        // eta zeta
        // dealiasmat.*(fft(  ky*realB.*idft + 2*ky*realDzB.*rilap2kz + ky*realDzzB.*rilap2  ))
        Aop_m42_ =
            (  fft_ilap2_[ind_ky].array().colwise()*(kyctmp_*rDzzBy_tmp_) +
             fft_kzilap2_[ind_ky].array().colwise()*((2.0*kyctmp_)*rDzBy_tmp_)  +
             fft_identity_.array().colwise()*(kyctmp_*rBy_tmp_)  ).matrix();
        fft_.for_2D_dim1(Aop_m42_.data());
        dealias(Aop_m42_.data());
        
        
        // eta b
        // -q*K.kz
        Aop_v43_ = ((-q_)*kz_).matrix().asDiagonal();
        
        // eta eta
        // 0
        
        //                                    //
        ////////////////////////////////////////

    
    if (QL_YN_) {
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
        
        double ftfac = mult_fac/(NZ_*NZ_); // NZ^2 factor for ifft is included here - more efficient since scalar
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
        fft_.back_2D_dim1(reynolds_mat_tmp_.data());
        reynolds_mat_tmp_ = reynolds_mat_tmp_.conjugate();
        fft_.back_2D_dim2(reynolds_mat_tmp_.data());
        reynolds_mat_tmp_ = reynolds_mat_tmp_.conjugate(); // Not necessary, but better not to be too confusing
        
//        std::cout << "kx = " << kxtmp_ << ", ky = " << kytmp_ << std::endl;
//        std::cout << reynolds_mat_tmp_ << std::endl;
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
        fft_.back_2D_dim1(reynolds_mat_tmp_.data());
        reynolds_mat_tmp_ = reynolds_mat_tmp_.conjugate();
        fft_.back_2D_dim2(reynolds_mat_tmp_.data());
        reynolds_mat_tmp_ = reynolds_mat_tmp_.conjugate(); // Not necessary, but better not to be too confusing
        // Keep a running sum
        bzuy_m_uzby_d_ -= ftfac*(reynolds_mat_tmp_.diagonal().array().real());
        
        //////                                               //////
        ///////////////////////////////////////////////////////////
    }
        
        /////////////////////////////////////////
        //    Multiply to get RHS              //
        
        // Should experiment with putting this before Reynolds stress, may help stability
        
        // With the block method of matrix multiplication I can no longer use Ckl_out = Ckl_in
        // to save memory in the integrator!
        
        
        // A*C matrix product - rather complicated due to all the submatrices
        
        
//        Ckl_out[i] += Ckl_out[i].adjoint().eval(); // Nearly 2* faster than above version
//        Ckl_out[i] += Qkl_tmp_.cast<dcmplx>().matrix().asDiagonal();
        Ckl_out[i].setZero();
        
        //        // PRINTING FOR DEBUG
        //        for (int i=0; i<Cdimz(); ++i) {
        //            for (int j=0; j<Cdimz(); ++j){
        //                if (imag(Aop_tmp_(i,j))<0.00001 && imag(Aop_tmp_(i,j))>-0.00001) {
        //                    Aop_tmp_(i,j).imag(0);
        //                }
        //                if (real(Aop_tmp_(i,j))<0.00001 && real(Aop_tmp_(i,j))>-0.00001) {
        //                    Aop_tmp_(i,j).real(0);
        //                }
        //            }
        //        }
        //        std::cout << "kx = " << kxtmp_ << ", ky = " << kytmp_ << std::endl;
        //        std::cout << Aop_tmp_ << std::endl;
        
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
    if (QL_YN_) {
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
    } else { // Set Reynolds stress to zero
        bzux_m_uzbx_c_.setZero();
        bzuy_m_uzby_c_.setZero();
    }
//    std::cout << "kx = " << kxtmp_-q_*dt_lin*kytmp_ << ", ky = " << kytmp_ << std::endl;
//    std::cout << bzux_m_uzbx_c_.transpose() << std::endl;
//    std::cout << bzuy_m_uzby_c_.transpose() << std::endl;
//    
    //////////////////////////////////////
    ////   MEAN FIELDS    ////////////////
    // In some integrators MFout will be the same as MFin.
    // Because of this, important to calculate MFout[1] first, since it depends on MFin[1], and this would be updated otherwise
    MFout[1] = -q_*MFin[0] + bzuy_m_uzby_c_;
    MFout[0] = bzux_m_uzbx_c_;
    

}


//////////////////////////////////////////////////////////
///   Initialization of the linear operators        //////

void MHD_BQlin::linearOPs_Init(double t0, doubVec * linop_MF, doubVec * linop_Ckl_old) {
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
void MHD_BQlin::dealias(dcmplx *arr) {
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
void MHD_BQlin::dealias(dcmplxVec& vec) {
    for (int j=delaiasBnds_[0]; j<delaiasBnds_[1]+1; ++j){
        vec(j) = 0;
    }
}


//////////////////////////////////////////////////////////
//////                                              //////
//////          ENERGY, ANGULAR MOMENTUM etc.       //////
//////                                              //////
//////////////////////////////////////////////////////////

void MHD_BQlin::Calc_Energy_AM_Diss(TimeVariables& tv, double t, const dcmplxVec *MFin, const dcmplxMat *Cin ) {
    // Energy, angular momentum and dissipation of the solution MFin and Cin
    // TimeVariables class stores info and data about energy, AM etc.
    // t is time
    
    // OUTPUT: energy[1] and [2] contain U and B mean field energies (energy[1]=0 for this)
    // energy[3] and [4] contain u and b fluctuating energies.
    
    double energy_u=0, energy_b=0;
    double energy_u_f=0, energy_b_f=0; // Receive buffer for MPI_Reduce

    int mult_fac;// Factor to account for only doing half fft sum.
    /////////////////////////////////////
    //   ALL OF THIS IS PARALLELIZED
    for (int i=0; i<Cdimxy();  ++i){
        
        int k_i = i + index_for_k_array(); // k index
        int ind_ky = ky_index_[k_i];
        // Form Laplacians using time-dependent kx
        kytmp_=ky_[k_i].imag();
        kxtmp_=kx_[k_i].imag() + q_*t*kytmp_;
        ilap2tmp_ = ilap2_[ind_ky];
        lapFtmp_ = -kxtmp_*kxtmp_ + lap2_[ind_ky];
        
        mult_fac = 2;
        if (kytmp_== 0.0 )
            mult_fac = 1; // Only count ky=0 mode once

        // Use Qkl_tmp_ for Mkl to save memory
        lapFtmp_ = lapFtmp_*ilap2tmp_;
        Qkl_tmp_ << lapFtmp_, -ilap2tmp_,lapFtmp_, -ilap2tmp_;

        
                // Energy = trace(Mkl*Ckl), Mkl is diagonal
        Qkl_tmp_ = mult_fac*Qkl_tmp_ * Cin[i].real().diagonal().array();
        
        
//        std::cout << "kx = " << kx_[i] << ", ky = " <<ky_[i] << std::endl;
//        std::cout << Qkl_tmp_.sum() << std::endl;
        int num_u_b = num_fluct_/2; // Keep it clear where numbers are coming from.
        energy_u += Qkl_tmp_.head(NZ_*num_u_b).sum();
        energy_b += Qkl_tmp_.tail(NZ_*num_u_b).sum();
        
//        std::cout << Qkl_tmp_ << std::endl;
//        std::cout << "a" << std::endl;

    }
    // Put the energy on processor 0
    mpi_.SumReduce_doub(&energy_u,&energy_u_f,1);
    mpi_.SumReduce_doub(&energy_b,&energy_b_f,1);
//    energy_u_f = energy_u; energy_b_f = energy_b;   // If in-place version is desired
//    mpi_.SumReduce_IP_doub(&energy_u_f,1);
//    mpi_.SumReduce_IP_doub(&energy_b_f,1);
    
    if (mpi_.my_n_v() == tv.tv_mpi_node()) {
        ////////////////////////////////////////////
        ///// All this is only on processor 0  /////
        double divfac=1.0/(Nxy_[0]*Nxy_[1]*2*NZ_);
        energy_u_f = energy_u_f*divfac*divfac;
        energy_b_f = energy_b_f*divfac*divfac;
        
        
        ///////////////////////////////////////
        ///       MEAN FIELDS            //////
        // Only need to calculate on one processor
        double energy_MU=0, energy_MB=0;
        
        energy_MB = (MFin[0].abs2().sum() + MFin[1].abs2().sum())/(NZ_*NZ_);
        
        ///////////////////////////////////////
        //////         OUTPUT            //////
        double* en_point = tv.current_energy();
        en_point[0] = energy_MU/2;
        en_point[1] = energy_MB/2;
        en_point[2] = energy_u_f/2;
        en_point[3] = energy_b_f/2;
        
        ///// All this is only on processor 0  /////
        ////////////////////////////////////////////
    }
    
    // Increase the number to save
    tv.increase_savenum();
    
}







//////////////////////////////////////////////////////
//         PRIVATE FUNCTIONS FOR INITIALIZATION     //

// fft of the identity matrix
void MHD_BQlin::Set_fft_identity_(void) {
    fft_identity_ = dcmplxMat::Identity(NZ_,NZ_);
    fft_.back_2D_dim1(fft_identity_.data());
    fft_identity_ /= NZ_;
}

// Create and store arrays of lap2 to save computation of ffts
void MHD_BQlin::Define_Lap2_Arrays_(void){
    // Lap2 quantities are stored on all processors, even though I could save some memory here. Easy to change later if a problem, the amount of memory may actually be significant...
    
    // Maximum ky index
    int number_of_ky = ky_index_[Cdimxy_full()-1]+1;
    // Quick error check
    if (Cdimxy_full()%number_of_ky!=0)
        std::cout << "Warning, something funny going on in Define_Lap2_Arrays_" << std::endl;
    // Assign data to arrays
    lap2_ = new doubVec[ number_of_ky ];
    ilap2_ = new doubVec[ number_of_ky ];
    fft_ilap2_ = new dcmplxMat[ number_of_ky ];
    fft_kzilap2_ = new dcmplxMat[ number_of_ky ];
    fft_kz2ilap2_ = new dcmplxMat[ number_of_ky ];

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
        
        // TAKE FFTs
        fft_ilap2_[i] = (ilap2tmp_/NZ_).cast<dcmplx>().matrix().asDiagonal();
        fft_.back_2D_dim1(fft_ilap2_[i].data());// NZ factor is added in when vector
        
        fft_kzilap2_[i] = (kz_*(ilap2tmp_/NZ_).cast<dcmplx>()).matrix().asDiagonal();
        fft_.back_2D_dim1(fft_kzilap2_[i].data());
        
        fft_kz2ilap2_[i] = (kz2_*ilap2tmp_/NZ_).cast<dcmplx>().matrix().asDiagonal();
        fft_.back_2D_dim1(fft_kz2ilap2_[i].data());
        
    }

}

//////////////////////////////////////////////////////


