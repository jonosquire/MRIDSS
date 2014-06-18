//
//  General_Auxiliary.cc
//  MRIDSS
//
//  Created by Jonathan Squire on 5/5/14.
//  Copyright (c) 2014 Jonathan Squire. All rights reserved.
//

#include "General_Auxiliary.h"


// Remapping for the shearing box
// This is done in the opposite way to my matlab code -- rather than moving C, just change kx. This then requires nothing complicated and no MPI communication
void ShearingBox_Remap(Model* mod, dcmplxMat* Ckl){
// FRIEND TO THE MODEL CLASS - THIS CHANGES THE kx_ array;
    dcmplx kxLfac(0,-1/(2*PI/mod->L_[0]) );
    dcmplx kyLfac(0,-1/(2*PI/mod->L_[1]) );
    
    dcmplx * kxp = mod->kx_;// Pointer to kxp data
    for (int i=0; i<mod->Cdimxy(); ++i) {
        // Find new kx vector
        int k_i = i + mod->index_for_k_array(); // Index in kx_, ky_ for each process
        kxp[k_i] = kxp[k_i]*kxLfac + mod->ky_[k_i]*kyLfac;
        
        // If kx_ is too large zero out Ckl and set kx negative
        double overlap = kxp[k_i].real()+1 - mod->Nxy_[0]/2;
        if (overlap>0) {
            Ckl[i].setZero();
            kxp[k_i] = -mod->Nxy_[0]/2+overlap;
        }
        // Un-normalize again
        kxp[k_i] = kxp[k_i]/kxLfac;
    }
}



//////////////////////////////////////////////////////////
//              CHECK SOLUTION IS OK                    //
void CheckSolution(dcmplxVec* MF, dcmplxMat* Ckl){
    // FRIEND TO THE MODEL CLASS
    // Checks important aspects of solution e.g., hasn't gotten too large
    
    // Mean field is real
    if (MF[1].real().abs().sum()/MF[1].imag().abs().sum() < 1e10 ){
        std::cout << "<<<<<<<<<<<<<<>>>>>>>>>>>>>>" << std::endl;
        std::cout << "ERROR: Mean field has developed an imaginary part!!" << std::endl;
        std::cout << "<<<<<<<<<<<<<<>>>>>>>>>>>>>>" << std::endl;
        ABORT;
    }
    
    // Stability
    if (MF[1].abs().maxCoeff()>1e20 || MF[0].abs().maxCoeff()>1e20) {
        std::cout << "<<<<<<<<<<<<<<>>>>>>>>>>>>>>" << std::endl;
        std::cout << "ERROR: Solution is very large, probably unstable!!" << std::endl;
        std::cout << "<<<<<<<<<<<<<<>>>>>>>>>>>>>>" << std::endl;
        ABORT;
    }
    if (!std::isfinite(MF[1].abs().sum()) || !std::isfinite(MF[1].abs().sum())) {
        std::cout << "<<<<<<<<<<<<<<>>>>>>>>>>>>>>" << std::endl;
        std::cout << "ERROR: Solution contains NaN or Inf!!" << std::endl;
        std::cout << "<<<<<<<<<<<<<<>>>>>>>>>>>>>>" << std::endl;
        ABORT;
    }
    
}

//////////////////////////////////////////////////////////


/////////////////////////////////////////////////////
//                                                 //
//                 fftw_Plans                      //
//                                                 //
/////////////////////////////////////////////////////

// fftw plans constructor - creates the plans
void fftwPlans::calculatePlans( int NZ ) {
    
    
    ////////////////////////////////////////////////
    //      CREATE TEMPORARY EIGEN OBJECTS FOR PLANS
    //  I'm relatively sure it doesn't create fftw problems if you later delete these...
    
    // Ckl temp
    dcmplxMat  Ckl_tmp(NZ,NZ);
    dcmplx * Ckl_p = Ckl_tmp.data();
    
    // Vector temps
    dcmplxVec MF_tmp1( NZ );
    dcmplx * MF_p1 = MF_tmp1.data();

    
    ////////////////////////////////////////////////
    //     2D PLANS - full data for Reynold's stress
    //
    // Need fft(fft( Ckl )')' (in Matlab notation), while standard 2-D fft is
    // fft(fft( Ckl ).').', hence the separate calculation for each dimension.
    //
    // All transforms can be in-place

    
    
    t2D_for_dim1_ = fftw_plan_many_dft(1, &NZ, NZ,
                                       CAST_T0_FFTW(Ckl_p), NULL, 1, NZ,
                                       CAST_T0_FFTW(Ckl_p), NULL, 1, NZ,
                                       FFTW_FORWARD, MY_FFTWPLAN);
    t2D_for_dim2_ = fftw_plan_many_dft(1, &NZ, NZ,
                                       CAST_T0_FFTW(Ckl_p), NULL, NZ, 1,
                                       CAST_T0_FFTW(Ckl_p), NULL, NZ, 1,
                                       FFTW_FORWARD, MY_FFTWPLAN);
    t2D_back_dim1_ = fftw_plan_many_dft(1, &NZ, NZ,
                                        CAST_T0_FFTW(Ckl_p), NULL, 1, NZ,
                                        CAST_T0_FFTW(Ckl_p), NULL, 1, NZ,
                                        FFTW_BACKWARD, MY_FFTWPLAN);
    t2D_back_dim2_ = fftw_plan_many_dft(1, &NZ, NZ,
                                        CAST_T0_FFTW(Ckl_p), NULL, NZ, 1,
                                        CAST_T0_FFTW(Ckl_p), NULL, NZ, 1,
                                        FFTW_BACKWARD, MY_FFTWPLAN);

    
    ////////////////////////////////////////////////
    //   1D Plans
    //  Need both a 1-D plan and a 1-D*NZ plan
    //
    // Only in place transforms have been useful, could delete others
   
    // True 1-D plans
    t1D_for_ = fftw_plan_dft_1d(NZ, CAST_T0_FFTW(MF_p1), CAST_T0_FFTW(MF_p1),
                                 FFTW_FORWARD, MY_FFTWPLAN);
    t1D_back_ = fftw_plan_dft_1d(NZ, CAST_T0_FFTW(MF_p1), CAST_T0_FFTW(MF_p1),
                                  FFTW_BACKWARD, MY_FFTWPLAN);
    
    
    // Flag to use in destructor
    plans_calculated_ = 1;
    
}

// Destructor
fftwPlans::~fftwPlans() {
    ///////////////////////////////////////////////////////////////////////
    // This deletes data that has been created in fftwPlans::calculatePlans

   // std::cout << "Destroying " << std::endl;
    if (plans_calculated_) {
        // Ckl plans
        fftw_destroy_plan(t2D_for_dim1_);
        fftw_destroy_plan(t2D_for_dim2_);
        fftw_destroy_plan(t2D_back_dim1_);
        fftw_destroy_plan(t2D_back_dim2_);
        
        // MF plans
        fftw_destroy_plan(t1D_for_);
        fftw_destroy_plan(t1D_back_);
    }
    
    fftw_cleanup();
    
}

//                                                 //
/////////////////////////////////////////////////////





/////////////////////////////////////////////////////
//                                                 //
//                 TimeVariables                   //
//                                                 //
/////////////////////////////////////////////////////

// Current position
int TimeVariables::curr_pos_ = 0;

// Constructor
TimeVariables::TimeVariables(Inputs SP, int width, int numMF, int mpinode) :
nsteps_(SP.nsteps/SP.timevar_save_nsteps + 1), width_(width), mpi_node_(mpinode),
simulation_dir_(SP.simulation_dir), MFlen_(0), numMF_(numMF),
en_save_Q_(SP.energy_save_Q), AM_save_Q_(SP.AM_save_Q), diss_save_Q_(SP.dissipation_save_Q),
rey_save_Q_(SP.reynolds_save_Q),
MF_save_Q_(SP.mean_field_save_Q)
{
    ////////////////////////////////
    //    THIS NEEDS TIDYING!!!!!!!!
    
    // Main variables
    // Store only on 1-processor - presumbably this is better
    // TODO: Only allocate necessary memory?
    if (mpi_node_ == 0) {
        energy_ = new double*[nsteps_];
        angular_momentum_ = new double*[nsteps_];
        dissipation_ = new double*[nsteps_];
        reynolds_ = new double*[nsteps_];
        for (int i=0; i<nsteps_; ++i) {
            energy_[i] = new double[width_];
            angular_momentum_[i] = new double[width_];
            dissipation_[i] = new double[width_];
            reynolds_[i] = new double[numMF_];
        }
    }
    
    // Saving to disk - TODO improve to save at regular intervals?
    fname_energy_ = simulation_dir_ + "energy.dat";
    fname_angular_momentum_ = simulation_dir_ + "angular_momentum.dat";
    fname_dissipation_ = simulation_dir_ + "dissipation.dat";
    // Reynolds stress
    fname_reynolds_ = simulation_dir_ + "reynolds_stress.dat";
    // Mean fields
    fname_mean_fields_ = simulation_dir_ + "mean_fields.dat";
    
    if (mpi_node_ == 0) {
        if (en_save_Q_) {
            fileS_energy_.open(fname_energy_.c_str(), std::ios::binary);
            if (!fileS_energy_.is_open() ) {
                std::cout << "WARNING: " << fname_energy_ <<  " file unable to be opened" << std::endl;
            }
        }
        if (AM_save_Q_) {
            fileS_angular_momentum_.open(fname_angular_momentum_.c_str(), std::ios::binary);
            if (!fileS_angular_momentum_.is_open() ) {
                std::cout << "WARNING: " << fname_angular_momentum_ <<  " file unable to be opened" << std::endl;
            }
        }
        if (diss_save_Q_) {
            fileS_dissipation_.open(fname_dissipation_.c_str(), std::ios::binary);
            if (!fileS_dissipation_.is_open() ) {
                std::cout << "WARNING: " << fname_dissipation_ <<  " file unable to be opened" << std::endl;
            }
        }
        // Reynolds stress
        if (rey_save_Q_) {
            fileS_reynolds_.open(fname_reynolds_.c_str(), std::ios::binary);
            if (!fileS_reynolds_.is_open() ) {
                std::cout << "WARNING: " << fname_reynolds_ <<  " file unable to be opened" << std::endl;
            }
        }
        // Mean fields
        if (MF_save_Q_) {
            fileS_mean_fields_.open(fname_mean_fields_.c_str(), std::ios::binary);
            if (!fileS_mean_fields_.is_open() ) {
                std::cout << "WARNING: " << fname_mean_fields_ <<  " file unable to be opened" << std::endl;
            }
        }
        
        ////////////////////////////////////////
        //        Save time data
        std::ofstream fileS_time;
        fileS_time.open((simulation_dir_+"time_vec.dat").c_str(),std::ios::binary);
        double *t_vec = new double[nsteps_];
        double dttmp = SP.timvar_save_interval;
        for (int i=0; i<nsteps_; ++i) {
            t_vec[i] = i*dttmp;
        }
        // Write
        fileS_time.write( (char*) t_vec, sizeof(t_vec)*nsteps_);
        // Clean-up
        delete [] t_vec;
        fileS_time.close();
        ////////////////////////////////////////
    }
    
    
}

// Destructor
TimeVariables::~TimeVariables() {
    if (mpi_node_==0) {
        for (int i=0; i<nsteps_; ++i) {
            delete[] energy_[i];
            delete[] angular_momentum_[i];
            delete[] dissipation_[i];
            delete[] reynolds_[i];
        }
    }
    delete[] energy_;
    delete[] angular_momentum_;
    delete[] dissipation_;
    delete[] reynolds_;
    
    fileS_energy_.close();
    fileS_angular_momentum_.close();
    fileS_dissipation_.close();
    fileS_reynolds_.close();
    fileS_mean_fields_.close();
    
}


// Save the data to disk
void TimeVariables::Save_Data(){
    // This can be run on all processes, but file will only be open on processor 0
    
    if (fileS_energy_.is_open()) {
        // Must write each individually since don't know that it is continuous in memory
        for (int i=0; i<nsteps_; ++i) {
            fileS_energy_.write( (char*) energy_[i], sizeof(energy_)*width_);
        }
        // Read into matlab as a double (4*nsteps) array
    }
    if (fileS_angular_momentum_.is_open()) {
        // Must write each individually since don't know that it is continuous in memory
        for (int i=0; i<nsteps_; ++i) {
            fileS_angular_momentum_.write( (char*) angular_momentum_[i], sizeof(energy_)*width_);
        }
        // Read into matlab as a double (4*nsteps) array
    }
    if (fileS_dissipation_.is_open()) {
        // Must write each individually since don't know that it is continuous in memory
        for (int i=0; i<nsteps_; ++i) {
            fileS_dissipation_.write( (char*) dissipation_[i], sizeof(energy_)*width_);
        }
        // Read into matlab as a double (4*nsteps) array
    }
    if (fileS_reynolds_.is_open()) {
        // Must write each individually since don't know that it is continuous in memory
        for (int i=0; i<nsteps_; ++i) {
            fileS_reynolds_.write( (char*) reynolds_[i], sizeof(reynolds_)*numMF_);
        }
        // Read into matlab as a double (4*nsteps) array
    }


    
}


// Save Mean field data
void TimeVariables::Save_Mean_Fields(dcmplxVec *MF, fftwPlans& fft){
    // Could certainly be improved - just dumps the data into fname_mean_fields
    
    // IN MATLAB
    // Saves MF[i] sequentially from 0 to numMF. Data is double
    
    // Not very fast I'm sure, but won't be important
    // is_open() will only evaluate to true on proc 1
    if (fileS_mean_fields_.is_open()) {
        // Initialize data
        if (MFlen_ == 0 ) {
            MFlen_ = MF[0].size();
            MFdata_c_ = dcmplxVec( MFlen_);
            MFdata_d_ = doubVec (MFlen_);
        }
        
        // Save each MF variable sequentially
        for (int i=0; i<numMF_; ++i) {
            // Take fft - could be done faster with c_to_r fft
            MFdata_c_ = MF[i];
            fft.back_1D(MFdata_c_.data());
            MFdata_d_ = MFdata_c_.real()/MFlen_;
            //Save
            fileS_mean_fields_.write( (char*) MFdata_d_.data(), sizeof(double)*MFlen_);
        }
        
        
    }
}
//                                                 //
/////////////////////////////////////////////////////

