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
    fft_.calculatePlans( NZ_,  nz_Cfull_);
    
    // Delalias setup
    delaiasBnds_[0] = NZ_/3+1; // Delete array including these elements
    delaiasBnds_[1] = NZ_-NZ_/3-1;

    
    ////////////////////////////////////////////////////
    // Useful arrays to save computation
    // ifft of identity and 1/lap2 etc.
    
    Set_fft_identity_();// fft of identity
    Define_Lap2_Arrays_(); // Lap2 and ffts - Allocates data to lap2_,ilap2_,fft_ilap2_,fft_kzilap2_,fft_kz2ilap2_
    
    
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
    
    // A operator matrix
    Aop_tmp_= dcmplxMat(nz_Cfull_,nz_Cfull_);
    Aop_block_tmp_= dcmplxMat(NZ_,NZ_);
    Aop_block_tmp2_= dcmplxMat(NZ_,NZ_);

    
    // Real versions of B
    rBy_tmp_ = dcmplxVec::Zero(NZ_);
    rDzBy_tmp_ = dcmplxVec::Zero(NZ_);
    rDzzBy_tmp_ = dcmplxVec::Zero(NZ_);
    
    // fft of full C matrix
    ifft_Ckl_tmp_ = dcmplxMat(NZ_,NZ_);
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
    fft_.back_MF1D(MFin[1].data(), rBy_tmp_.data());
    rBy_tmp_ /= NZ_; // fft back doesn't include normalization
    rDzBy_tmp_ = kz_*MFin[1];
    fft_.back_IP_MF1D(rDzBy_tmp_.data());
    rDzBy_tmp_ /= NZ_;
    rDzzBy_tmp_ = kz2_*MFin[1];
    fft_.back_IP_MF1D(rDzzBy_tmp_.data());
    rDzzBy_tmp_ /= NZ_;
    
    // Operator (A) matrix -- to be filled up
    Aop_tmp_.setZero();
    
    
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
            Qkl_tmp_ = Qkl_tmp_.abs();
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
        Aop_tmp_.block(0, 0, NZ_, NZ_) += (-2*q_*kxctmp_*kyctmp_*ilapFtmp_).cast<dcmplx>().matrix().asDiagonal();
        
        // u zeta Component
        // 2*ilapF.*K.kz
        Aop_tmp_.block(0, NZ_, NZ_, NZ_) += (2*ilapFtmp_*kz_).cast<dcmplx>().matrix().asDiagonal();
        
        // u b Component
        // dealiasmat.*( fft(ky*realB.*idft) + 2*kxt^2*ky*ilapF*fft(realDzB.*rilap2kz) )
        Aop_block_tmp_ = (fft_identity_.array().colwise()*(kyctmp_*rBy_tmp_)).matrix();
        fft_.for_IP_MF2D(Aop_block_tmp_.data());
        Aop_block_tmp2_ = (fft_kzilap2_[ind_ky].array().colwise()*rDzBy_tmp_).matrix();
        fft_.for_IP_MF2D(Aop_block_tmp2_.data());
        Aop_block_tmp_ +=
            ((2.0*kxctmp_*kxctmp_*kyctmp_)*ilapFtmp_.cast<dcmplx>()).matrix().asDiagonal()*
                Aop_block_tmp2_;
        dealias(Aop_block_tmp_.data()); // Could also do delaliasing in assignment
        // Final assignment
        Aop_tmp_.block(0, 2*NZ_, NZ_, NZ_) = Aop_block_tmp_;
        
        // u eta component
        // dealiasmat.*( 2*kxt*ky^2*ilapF*fft( realDzB.*rilap2 )  )
        Aop_block_tmp_ = (fft_ilap2_[ind_ky].array().colwise()*rDzBy_tmp_).matrix();
        fft_.for_IP_MF2D(Aop_block_tmp_.data());
        Aop_block_tmp_ =
            ((2.0*kxctmp_*kyctmp_*kyctmp_)*ilapFtmp_.cast<dcmplx>()).matrix().asDiagonal()*
                Aop_block_tmp_;
        dealias(Aop_block_tmp_.data());
        // Final assignment
        Aop_tmp_.block(0, 3*NZ_, NZ_, NZ_) = Aop_block_tmp_;
        
        
        // ROW 2
        
        // zeta u
        // (q-2)*K.kz
        Aop_tmp_.block(NZ_, 0, NZ_, NZ_) += ((q_-2.0)*kz_).matrix().asDiagonal();
        
        // zeta zeta
        // 0
        
        // zeta b
        // dealiasmat.*( fft(-kxt*realDzzB.*rilap2kz-kxt*realDzB.*idft) )
        Aop_block_tmp_ = (  fft_kzilap2_[ind_ky].array().colwise()*((-kxctmp_)*rDzzBy_tmp_) +
            fft_identity_.array().colwise()*(-kxctmp_*rDzBy_tmp_)  ).matrix();
        fft_.for_IP_MF2D(Aop_block_tmp_.data());
        dealias(Aop_block_tmp_.data());
        // Final assignment
        Aop_tmp_.block(NZ_, 2*NZ_, NZ_, NZ_) = Aop_block_tmp_;
        
        // zeta eta
        // dealiasmat.*( fft( ky*realB.*idft - ky*realDzzB.*rilap2 ))
        Aop_block_tmp_ = (  fft_ilap2_[ind_ky].array().colwise()*((-kyctmp_)*rDzzBy_tmp_) +
            fft_identity_.array().colwise()*(kyctmp_*rBy_tmp_)  ).matrix();
        fft_.for_IP_MF2D(Aop_block_tmp_.data());
        dealias(Aop_block_tmp_.data());
        // Final assignment
        Aop_tmp_.block(NZ_, 3*NZ_, NZ_, NZ_) = Aop_block_tmp_;

        
        // ROW 3

        // b u
        // dealiasmat.*( fft(ky*realB.*idft) )
        Aop_block_tmp_ = (  fft_identity_.array().colwise()*(kyctmp_*rBy_tmp_)  ).matrix();
        fft_.for_IP_MF2D(Aop_block_tmp_.data());
        dealias(Aop_block_tmp_.data());
        // Final assignment
        Aop_tmp_.block(2*NZ_, 0, NZ_, NZ_) = Aop_block_tmp_;
        
        // b zeta
        // 0
        
        // b b
        // 0 (except linear)
        
        // b eta
        //0
        
        
        // ROW 4
        
        // eta u
        // dealiasmat.*(fft( 2*kxt*realDzB.*ifft(ilap2.*K.kz.^2) + kxt*realDzzB.*rilap2kz - kxt*realDzB.*idft  ))
        Aop_block_tmp_ =
            (  fft_kz2ilap2_[ind_ky].array().colwise()*((2.0*kxctmp_)*rDzBy_tmp_) +
            fft_kzilap2_[ind_ky].array().colwise()*(kxctmp_*rDzzBy_tmp_)  +
            fft_identity_.array().colwise()*((-kxctmp_)*rDzBy_tmp_)  ).matrix();
        fft_.for_IP_MF2D(Aop_block_tmp_.data());
        dealias(Aop_block_tmp_.data());
        // Final assignment
        Aop_tmp_.block(3*NZ_, 0, NZ_, NZ_) = Aop_block_tmp_;
        
        
        // eta zeta
        // dealiasmat.*(fft(  ky*realB.*idft + 2*ky*realDzB.*rilap2kz + ky*realDzzB.*rilap2  ))
        Aop_block_tmp_ =
            (  fft_ilap2_[ind_ky].array().colwise()*(kyctmp_*rDzzBy_tmp_) +
             fft_kzilap2_[ind_ky].array().colwise()*((2.0*kyctmp_)*rDzBy_tmp_)  +
             fft_identity_.array().colwise()*(kyctmp_*rBy_tmp_)  ).matrix();
        fft_.for_IP_MF2D(Aop_block_tmp_.data());
        dealias(Aop_block_tmp_.data());
        // Final assignment
        Aop_tmp_.block(3*NZ_, NZ_, NZ_, NZ_) = Aop_block_tmp_;
        
        // eta b
        // -q*K.kz
        Aop_tmp_.block(3*NZ_, 2*NZ_, NZ_, NZ_) += ((-q_)*kz_).matrix().asDiagonal();
        
        // eta eta
        // 0
        
        //                                    //
        ////////////////////////////////////////
        
        
        /////////////////////////////////////////
        //    Multiply to get RHS              //
        
        // In some integrators Ckl_out will be the same as Ckl_in.
        // MUST REMEMBER THIS
        Ckl_out[i] = Aop_tmp_*Ckl_in[i] + Ckl_in[i]*Aop_tmp_.adjoint();
        Ckl_out[i] += Qkl_tmp_.cast<dcmplx>().matrix().asDiagonal();
        

        
        ////////////////////////////////////////
        //////   LINEAR PART
        // Need to re-evaluate laplacian, since different time.
        kxtmp_=kxtmp_ + q_*dt_lin*kytmp_;
        lapFtmp_ = (-kxtmp_*kxtmp_-kytmp_*kytmp_)+ kz2_;
        linop_Ckl[i] << nu_*lapFtmp_, nu_*lapFtmp_, eta_*lapFtmp_, eta_*lapFtmp_;
    }
    

    //////////////////////////////////////
    ////   MEAN FIELDS    ////////////////
    for (int i=0; i<num_MF_; ++i){
        MFout[i] = MFin[i];
    }

}


//////////////////////////////////////////////////////////
///   Initialization of the linear operators        //////

void MHD_BQlin::linearOPs_Init(double t0, doubVec * linop_MF, doubVec * linop_Ckl_old) {
    // Constructs initial (t=0) linear operators, since MF operator is constant this completely defines it for entire integration
    // Mean fields
    for (int i=0; i<num_MF_;  ++i){
        linop_MF[i] = eta_*kz2_;
    };
    
    // Ckl
    double kxtmp,kytmp; // Loop through x,y array

    for (int i=0; i<Cdimxy();  ++i){
        int k_i = i + index_for_k_array(); // k index
        kytmp=ky_[k_i].imag();
        kxtmp=kx_[k_i].imag() + q_*t0*kytmp;
        lapFtmp_ = -kxtmp*kxtmp-kytmp*kytmp+kz2_;
        linop_Ckl_old[i] << nu_*lapFtmp_, nu_*lapFtmp_, eta_*lapFtmp_, eta_*lapFtmp_;
    }
    
}


//////////////////////////////////////////////////////////
////                DEALIASING                      //////
// Applies dealiasing to matrix of size NZ_ in z Fourier space
void MHD_BQlin::dealias(dcmplx *arr) {
    // Could also use Block here, might be faster
    dcmplx zc(0,0);
    for (int i=delaiasBnds_[0]; i<=delaiasBnds_[1]; i++) {
        for (int j=delaiasBnds_[0]; j<=delaiasBnds_[1]; j++)
            arr[i+NZ_*j]=zc; // Column major, not that it makes any difference
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
        // Form Laplacians using time-dependent kx
        kytmp_=ky_[k_i].imag();
        kxtmp_=kx_[k_i].imag() + q_*t*kytmp_;
        lap2tmp_ = -kytmp_*kytmp_ + kz2_;
        lapFtmp_ = -kxtmp_*kxtmp_ + lap2tmp_;
        
        mult_fac = 2;
        if (kytmp_== 0.0 ) {
            lap2tmp_(0) = 1.0; // Shouldn't be energy in ky=kz=0 anyway
            mult_fac = 1; // Only count ky=0 mode once
        }
        
        // Use Qkl_tmp_ to save memory
        lap2tmp_ = (1/lap2tmp_).abs();
        lapFtmp_ = (lapFtmp_*lap2tmp_).abs();
        Qkl_tmp_ << lapFtmp_, lap2tmp_,lapFtmp_,lap2tmp_;
        
        // Energy = trace(Mkl*Ckl), Mkl is diagonal
        Qkl_tmp_ = mult_fac*Qkl_tmp_ * Cin[i].real().diagonal().array();
        
        int num_u_b = num_fluct_/2; // Keep it clear where numbers are coming from.
        energy_u += Qkl_tmp_.head(NZ_*num_u_b).sum();
        energy_b += Qkl_tmp_.tail(NZ_*num_u_b).sum();
        
    }
    // Put the energy on processor 0
    mpi_.SumReduce_doub(&energy_u,&energy_u_f,1);
    mpi_.SumReduce_doub(&energy_b,&energy_b_f,1);
    
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
    fft_.back_IP_MF2D(fft_identity_.data());
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
        fft_.back_IP_MF2D(fft_ilap2_[i].data());// NZ factor is added in when vector
        
        fft_kzilap2_[i] = (kz_*(ilap2tmp_/NZ_).cast<dcmplx>()).matrix().asDiagonal();
        fft_.back_IP_MF2D(fft_kzilap2_[i].data());
        
        fft_kz2ilap2_[i] = (kz2_*ilap2tmp_/NZ_).cast<dcmplx>().matrix().asDiagonal();
        fft_.back_IP_MF2D(fft_kz2ilap2_[i].data());
        
    }

}

//////////////////////////////////////////////////////


