//
//  Constant_Damping.cc
//  MRIDSS
//
//  Created by Jonathan Squire on 4/25/14.
//  Copyright (c) 2014 J Squire. All rights reserved.
//

#include "Constant_Damping.h"


// Model class for S3T/CE2 shearing box MHD model
// Operator matrix is simply -rI, i.e., constant damping. This is do check
// driving noise
//
// Derived from Model (model.h)


// Constructor for Constant_Damping
Constant_Damping::Constant_Damping(const Inputs& sp, MPIdata& mpi, fftwPlans& fft) :
num_MF_(2), num_fluct_(4),
q_(sp.q), f_noise_(sp.f_noise), nu_(sp.nu), eta_(sp.eta),
QL_YN_(sp.QuasiLinearQ),
Model(sp.NZ, sp.NXY , sp.L), // Dimensions - stored in base
mpi_(mpi) // MPI data
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
    
    
    // Delalias setup
    delaiasBnds_[0] = NZ_/3+1; // Delete array including these elements
    delaiasBnds_[1] = NZ_-NZ_/3-1;
    
    
    ////////////////////////////////////////////////////
    // Useful arrays to save computation
    // ifft of identity and 1/lap2 etc.
    Define_Lap2_Arrays_(); // Lap2 and ffts - Allocates data to lap2_,ilap2_,fft_ilap2_,fft_kzilap2_,fft_kz2ilap2_
    
    totalN_ = Nxy_[0]*Nxy_[1]*2*NZ_;

    mult_noise_fac_ = 1.0/(16*32*32); // Defined so that consistent with (16, 32, 32) results
    
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
    // A operator
    Aop_tmp_ = doubVec( nz_Cfull_ );
    
    
}


// Destructor
Constant_Damping::~Constant_Damping(){
    // Here I am destroying Base members from the derived destructor - is this bad?
    delete[] kx_;
    delete[] ky_;
    
    // useful fft arrays
    delete[] ky_index_;
    delete[] lap2_;
    delete[] ilap2_;
    
}


// Main EulerFluid function, fx = f(t,x) called at ecah time-step
void Constant_Damping::rhs(double t, double dt_lin,
                    dcmplxVec *MFin, dcmplxMat *Ckl_in,
                    dcmplxVec *MFout, dcmplxMat *Ckl_out,
                    doubVec * linop_Ckl) {
    
    // If desired, the function should work fine with Ckl_in = Ckl_out. This is
    // useful for integration schemes that do not require re-use of Ckl_in
    
    double r = 0.5;
    Aop_tmp_.setConstant(-r);
    
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
            // CANNOT USE AFTER THIS!
            if (kytmp_==0.0) { // Don't drive the ky=kz=0 component (variables not invertible)
                lapFtmp_(0)=0;
        }
//            
//            /////////////   DEBUGING ////////////////
//            /////      TO DELETE     ////////////////
//            double kxmax = 2*PI/L_[0]*7, kymax = 2*PI/L_[1]*15;
//            if (kxtmp_ - q_*t*kytmp_ > kxmax || kytmp_ > kymax) {
//                lapFtmp_.setZero();
//                lap2tmp_.setZero();
//            }
//            double kzmax = 2*PI/L_[2]*15;
//            for (int i=0; i<NZ_; ++i) {
//                if (abs(kz_(i)) > kzmax) {
//                    lapFtmp_(i) = 0;
//                    lap2tmp_(i) = 0;
//                }
//            }

            
            Qkl_tmp_ << lapFtmp_, lap2tmp_, lapFtmp_, lap2tmp_;
            Qkl_tmp_ = (f_noise_*f_noise_*totalN_*mult_noise_fac_)*Qkl_tmp_.abs();
        }
        ////////////////////////////////////////
        
        
        
        
        
        /////////////////////////////////////////
        //    Multiply to get RHS              //
        
        // Should experiment with putting this before Reynolds stress, may help stability
        Ckl_out[i] = Aop_tmp_.matrix().asDiagonal()*Ckl_in[i];
        
        
        // Add to adjoint
        Ckl_out[i] += Ckl_out[i].adjoint().eval();
        
        // Add driving noise
        Ckl_out[i] += Qkl_tmp_.cast<dcmplx>().matrix().asDiagonal();
        
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

    //    std::cout << "kx = " << kxtmp_-q_*dt_lin*kytmp_ << ", ky = " << kytmp_ << std::endl;
    //    std::cout << bzux_m_uzbx_c_.transpose() << std::endl;
    //    std::cout << bzuy_m_uzby_c_.transpose() << std::endl;
    //
    //////////////////////////////////////
    ////   MEAN FIELDS    ////////////////
 
    MFout[1].setZero();
    MFout[0].setZero();
    
}


//////////////////////////////////////////////////////////
///   Initialization of the linear operators        //////

void Constant_Damping::linearOPs_Init(double t0, doubVec * linop_MF, doubVec * linop_Ckl_old) {
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
//////                                              //////
//////          ENERGY, ANGULAR MOMENTUM etc.       //////
//////                                              //////
//////////////////////////////////////////////////////////

void Constant_Damping::Calc_Energy_AM_Diss(TimeVariables& tv, double t, const dcmplxVec *MFin, const dcmplxMat *Cin ) {
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
            
            int num_u_b = num_fluct_/2; // Split into u and b parts
            energy_u += Qkl_tmp_.head(NZ_*num_u_b).sum();
            energy_b += Qkl_tmp_.tail(NZ_*num_u_b).sum();
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
    
    if (mpi_.my_n_v() == 0) {
        ////////////////////////////////////////////
        ///// All this is only on processor 0  /////
        double divfac=1.0/totalN_;
        divfac = divfac*divfac;
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
            // Just save the previously assigned Reynolds stress to compare By with Bx stress
            double* rey_point = tv.current_reynolds();
            rey_point[0] = 0;
            rey_point[1] = 0;
            //
            //////////////////////////////////////
        }
        
        ///// All this is only on processor 0  /////
        ////////////////////////////////////////////
    }
    
    // Increase the number to save
    tv.Save_Data();
    
}








//////////////////////////////////////////////////////
//         PRIVATE FUNCTIONS FOR INITIALIZATION     //



// Create and store arrays of lap2 to save computation of ffts
void Constant_Damping::Define_Lap2_Arrays_(void){
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


