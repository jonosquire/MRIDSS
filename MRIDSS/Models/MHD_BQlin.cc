//
//  MHD_BQlin.cc
//  MRIDSS
//
//  Created by Jonathan Squire on 4/25/14.
//  Copyright (c) 2014 J Squire. All rights reserved.
//

#include "MHD_BQlin.h"


#define fftplan FFTW_ESTIMATE

// Model class for S3T/CE2 shearing box MHD model
// Basic Quasi-linear MHD
//
// Derived from Model (model.h)


// Constructor for MHD_BQlin
MHD_BQlin::MHD_BQlin(const Inputs& sp) :
q_(sp.q), f_noise_(sp.f_noise), nu_(sp.nu), eta_(sp.eta),
Model(sp.NZ, sp.NXY , sp.L) // Dimensions - stored in base
{
    // Dimensions that depend on model
    nz_Cfull_=num_fluct_*NZ_;
    // Calcuate kx and ky arrays - single dimesnsion NX*NY(/2) to simplify parallelization
    kx_ = new dcmplx[ nxy_full_ ];
    ky_ = new dcmplx[ nxy_full_ ];
    define_kxy_array(kx_, ky_, Nxy_, L_);
    // kz array - Eigen vector
    kz_ = dcmplxVec( NZ_ );
    kz2_ = doubVec( NZ_ );
    define_kz_array(kz_, NZ_, L_);
    kz2_ = (kz_*kz_).real(); // kz2_ is a double Eigen array (rather than dcmplx)
    
    // Temporary arrays to use in evaluation
    lapFtmp_ = doubVec( NZ_ ); //For (time-dependent) k^2
    lap2tmp_ = doubVec( NZ_ ); // For ky^2+kz^2 - could be pre-assigned
    // Qkl
    Qkl_tmp_ = doubVec( nz_Cfull_ );
}


// Destructor
MHD_BQlin::~MHD_BQlin(){
    // Here I am destroying Base members from the derived destructor - is this bad?
    delete[] kx_;
    delete[] ky_;
    
    fftw_destroy_plan(plan_2D_for);
    fftw_destroy_plan(plan_2D_back);
}


// Main EulerFluid function, fx = f(t,x) called at ecah time-step
void MHD_BQlin::rhs(double t, double dt_lin,
                   const dcmplxVec *MFin, const dcmplxMat *Ckl_in,
                   dcmplxVec *MFout, dcmplxMat *Ckl_out,
                   doubVec * linop_Ckl) {
    
    
    double r = 0.5; // Constant damping parameter
    
    /////////////////////////////////////
    //   ALL OF THIS MUST BE PARALLELIZED
    for (int i=0; i<nxy_full_;  ++i){
        
        // Full Loop containing all the main work on Ckl
        
        ///////////////////////////////////////
        ///// MAIN CKL EQUATIONS
        
        // Form Laplacians using time-dependent kx
        kytmp_=ky_[i].imag();
        kxtmp_=kx_[i].imag() + q_*t*kytmp_;
        lap2tmp_ = -kytmp_*kytmp_ + kz2_;
        lapFtmp_ = -kxtmp_*kxtmp_ + lap2tmp_;
        
        // Qkl
        if (kytmp_==0 && kxtmp_==0) {
            Qkl_tmp_.setZero(); // Don't drive the mean components! (really could leave this out entirely)
        }
        else {
            lapFtmp_ = lap2tmp_/lapFtmp_;
            if (kytmp_==0) { // Don't drive the ky=kz=0 component (variables not invertible)
                lapFtmp_(0)=0;
                lap2tmp_(0)=0;
            }
            Qkl_tmp_ << lapFtmp_, lap2tmp_, lapFtmp_, lap2tmp_;
            Qkl_tmp_ = Qkl_tmp_.abs();
        }
        
        Ckl_out[i] = -r*Ckl_in[i]-Ckl_in[i]*r;
        for (int j=0; j<nz_Cfull_; ++j) {
            // Check speed of this - eigen accessor rather fast with -O3
            Ckl_out[i](j,j) = Ckl_out[i](j,j).real() + Qkl_tmp_(j);
        }
        

        
        ////////////////////////////////////////
        //////   LINEAR PART
        // Need to re-evaluate laplacian, since different time.
        kxtmp_=kxtmp_ + q_*dt_lin*kytmp_;
        lapFtmp_ = (-kxtmp_*kxtmp_-kytmp_*kytmp_)+ kz2_;
        linop_Ckl[i] << nu_*lapFtmp_, nu_*lapFtmp_, eta_*lapFtmp_, eta_*lapFtmp_;
    }
    
    
    
//    for (int i=0; i<nxy_full_;  ++i) {
//        fftw_execute_dft(plan_2D_for,
//                         CAST_T0_FFTW(Ckl_out[i].data()) ,
//                         CAST_T0_FFTW(Ckl_out[i].data()));
//    }
    
    

    //////////////////////////////////////
    ////   MEAN FIELDS    ////////////////
    for (int i=0; i<num_MF_; ++i)
        MFout[i] = MFin[i];
    

};


//////////////////////////////////////////////////////////
///   Initialization of the linear operators        //////

void MHD_BQlin::linearOPs_Init(double t0, doubVec * linop_MF, doubVec * linop_Ckl_old) {
    // Constructs initial (t=0) linear operators, since MF operator is constant this completely defines it for entire integration
    // Mean fields
    for (int i=0; i<num_MF_;  ++i){
        for (int j=0; j<NZ_; ++j)
            linop_MF[i](j) = -0.0*eta_*kz_[j].imag()*kz_[j].imag();
    };
    
    // Ckl
    double kxtmp,kytmp; // Loop through x,y array

    for (int i=0; i<nxy_full_;  ++i){
        kytmp=ky_[i].imag();
        kxtmp=kx_[i].imag() + q_*t0*kytmp;
        lapFtmp_ = -kxtmp*kxtmp-kytmp*kytmp+kz2_;
        linop_Ckl_old[i] << nu_*lapFtmp_, nu_*lapFtmp_, eta_*lapFtmp_, eta_*lapFtmp_;
    }
    
}



//////////////////////////////////////////////////////////
//////                                              //////
//////          ENERGY, ANGULAR MOMENTUM etc.       //////
//////                                              //////
//////////////////////////////////////////////////////////

void MHD_BQlin::Calc_Energy(double* energy, double t, const dcmplxVec *MFin, const dcmplxMat *Cin ) {
    // Energy of the solution MFin and Cin
    // energy should be pointer array of LENGTH 4
    // t is time
    
    // OUTPUT: energy[1] and [2] contain U and B mean field energies (energy[1]=0 for this)
    // energy[3] and [4] contain u and b fluctuating energies.
    
    double energy_u=0, energy_b=0;
    int mult_fac;// Factor to account for only doing half fft sum.
    /////////////////////////////////////
    //   ALL OF THIS MUST BE PARALLELIZED
    for (int i=0; i<nxy_full_;  ++i){
        
        // Form Laplacians using time-dependent kx
        kytmp_=ky_[i].imag();
        kxtmp_=kx_[i].imag() + q_*t*kytmp_;
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
        //  Requires MPI_Reduce !!!

//        REMOVE!!!
//        std::cout << "(kx,ky)=(" << kx_[i].imag() << "," << kytmp_ <<")"<<std::endl;
//        std::cout << mult_fac*Qkl_tmp_.matrix().asDiagonal()*Cin[i].real() << std::endl;
        
        // Energy = trace(Mkl*Ckl), Mkl is diagonal
        Qkl_tmp_ = mult_fac*Qkl_tmp_ * Cin[i].real().diagonal().array();
        
        int num_u_b = num_fluct_/2; // Keep it clear where numbers are coming from.
        energy_u += Qkl_tmp_.head(NZ_*num_u_b).sum();
        energy_b += Qkl_tmp_.tail(NZ_*num_u_b).sum();
        
        

    }
    
    // MPI_Reduce(.....)  - put energy on one of the processes
    double divfac=1.0/(Nxy_[0]*Nxy_[1]*2*NZ_);
    energy_u = energy_u*divfac*divfac;
    energy_b = energy_b*divfac*divfac;
    ///////////////////////////////////////
    ///       MEAN FIELDS            //////
    double energy_MU=0, energy_MB=0;
    
    energy_MB = (MFin[0].abs2().sum() + MFin[1].abs2().sum())/(NZ_*NZ_);
    
    ///////////////////////////////////////
    //////         OUTPUT            //////
    energy[0] = energy_MU;
    energy[1] = energy_MB;
    energy[2] = energy_u;
    energy[3] = energy_b;

}

void MHD_BQlin::Calc_AngularMomentum(double* AM, double t,const dcmplxVec *MFin, const dcmplxMat *Cin ){} // Angular momentum
void MHD_BQlin::Calc_Dissipation(double* AM, double t,const dcmplxVec *MFin, const dcmplxMat *Cin ) {} // Dissipation
