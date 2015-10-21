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

// DONT USE THIS!!!
void ShearingBox_Remap(double qt, Model* mod, dcmplxMat* Ckl){
    //FRIEND TO THE MODEL CLASS -THIS CHANGES THE kx_ array
    double kxLfac=2*PI/mod->box_length(0);
    int nx = mod->Nxy_[0]-1;
    // THIS IS NOT PERFECT!! IS ASYMMETRICAL AS IT IS
    dcmplx * kxp = mod->kx_pointer();// Pointer to kx data
    dcmplx * kyp = mod->ky_pointer();// Pointer to ky data
    int k_i;
    double kxt;
    int numk;
    
    // Go through and fix everything
    // Since this is only called every qt, can be more than one "nx" out of range
    for (int i=0; i<mod->Cdimxy(); ++i) {
        // Find new kx vector
        // This will become different on each processor even though each is storing the whole array. This is a little strange but probably inconsequential (aside from slight memory usage).
        k_i = i + mod->index_for_k_array(); // Index in kx_, ky_ for each process
        kxt = kxp[k_i].imag() + qt*kyp[k_i].imag();
        if (kxt > nx/2.0*kxLfac) {
            
            //std::cout << "Remapped kx " << kxp[k_i].imag()+ qt*kyp[k_i].imag() << " ky " << kyp[k_i].imag() << std::endl;
            // Put kx back in correct range
            kxp[k_i] = kxp[k_i] - dcmplx(0,nx*kxLfac);
            // zero out Ckl
            Ckl[i].setZero();
        }
    }
    
}


// Remapping for the shearing box
// Using the Lithwick method of continuously remapping wavenumbers when they become too large
void ShearingBox_Continuous_Remap(double qt, Model* mod, dcmplxMat* Ckl){
    // FRIEND TO THE MODEL CLASS -THIS CHANGES THE kx_ array
    double kxLfac=2*PI/mod->box_length(0);
    int nx = mod->Nxy_[0]-1;
    
    dcmplx * kxp = mod->kx_pointer();// Pointer to kx data
    dcmplx * kyp = mod->ky_pointer();// Pointer to ky data
    int k_i;
    double kxt;
    for (int i=0; i<mod->Cdimxy(); ++i) {
        // Find new kx vector
        // This will become different on each processor even though each is storing the whole array. This is a little strange but probably inconsequential (aside from slight memory usage).
        k_i = i + mod->index_for_k_array(); // Index in kx_, ky_ for each process
        kxt = kxp[k_i].imag() + qt*kyp[k_i].imag();
        if (kxt > nx/2.0*kxLfac) {
            //std::cout << "Remapped kx " << kxp[k_i].imag()+ qt*kyp[k_i].imag() << " ky " << kyp[k_i].imag() << std::endl;
            // Put kx back in correct range
            kxp[k_i] = kxp[k_i] - dcmplx(0,nx*kxLfac);
            // zero out Ckl
            Ckl[i].setZero();
        }
    }
    
    
}



//////////////////////////////////////////////////////////
//              CHECK SOLUTION IS OK                    //
void CheckSolution(dcmplxVec* MF, dcmplxMat* Ckl, fftwPlans& fft){
    // Checks important aspects of solution e.g., hasn't gotten too large
    
    // Take an in place transform
    fft.back_1D(MF[1].data());
    MF[1] /= MF[1].size();
    
//    std::cout << MF[1].transpose() << std::endl;
    
    if (MF[1].real().abs().sum()/MF[1].imag().abs().sum() < 1e10 &&  MF[1].real().abs().sum()>1e-2){
        std::cout << "<<<<<<<<<<<<<<>>>>>>>>>>>>>>" << std::endl;
        std::cout << "ERROR: Mean field has developed an imaginary part!!" << std::endl;
        std::cout << "<<<<<<<<<<<<<<>>>>>>>>>>>>>>" << std::endl;
        ABORT;
    }
    // Transform back
    fft.for_1D(MF[1].data());
    
    // Stability
    if (MF[1].abs().maxCoeff()>1e20 || MF[0].abs().maxCoeff()>1e20) {
        std::cout << "<<<<<<<<<<<<<<>>>>>>>>>>>>>>" << std::endl;
        std::cout << "ERROR: Solution is very large, probably unstable!!" << std::endl;
        std::cout << "<<<<<<<<<<<<<<>>>>>>>>>>>>>>" << std::endl;
        ABORT;
    }
    
    //  THIS ISN'T WORKING FOR SOME REASON!!
    if (!std::isfinite(MF[1].abs().sum()) || !std::isfinite(MF[1].abs().sum())) {
        std::cout << "<<<<<<<<<<<<<<>>>>>>>>>>>>>>" << std::endl;
        std::cout << "ERROR: Solution contains NaN or Inf!!" << std::endl;
        std::cout << "<<<<<<<<<<<<<<>>>>>>>>>>>>>>" << std::endl;
        ABORT;
    }
    
}

//////////////////////////////////////////////////////////







//////////////////////////////////////////////
//   Save final state 

void Save_Full_Data_for_Restart(Model* fluidEqs, Inputs& SP, MPIdata& mpi, double t, dcmplxVec* MF, dcmplxMat* Ckl){
    // Saves, number of processors (int) and index in k array, then i (int) and t (double) (for kx), then MFs (dcmplx), then all of Ckl (dcmplx)
    
    std::string state_dir(SP.simulation_dir+"FINAL_STATE/");
    if (mpi.my_n_v() == 0) {
        tinydir_dir tmp_dir; // Check if directory exists
        int dir_status = tinydir_open(&tmp_dir, state_dir.c_str());
        tinydir_close(&tmp_dir);
        // Make directory if there isn't one
        if ( dir_status == -1 ){
            int status = mkdir(state_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
            if (status == -1) {
                std::cout << "Warning: Failed to create simulation directory: " << state_dir << std::endl;
            }
        }
        
    }
    mpi.Barrier();
    
    std::stringstream final_state_name;
    final_state_name << SP.simulation_dir+"FINAL_STATE/" << "FINAL_CKL_STATE_Proc" << mpi.my_n_v() << ".dat";
    std::ofstream final_state_file(final_state_name.str().c_str(), std::ios::binary);
    if (!final_state_file.is_open() ) {
        std::cout << "WARNING: " << final_state_name.str() <<  " file unable to be opened" << std::endl;
    }
    
    /////////////////////////////////////////
    // Write Data
    
    // Number of processors, index (on k array) for each processor
    final_state_file.write( (char *) mpi.total_n_p(), sizeof( int ) );
    int minxy = mpi.minxy_i(), maxxy = mpi.maxxy_i();
    final_state_file.write( (char *) &minxy, sizeof( int ) );
    final_state_file.write( (char *) &maxxy, sizeof( int ) );
    // i, just for simplicity
    final_state_file.write( (char *) &(SP.nsteps), sizeof( int ) );
    // t, just for simplicity
    final_state_file.write( (char *) &t, sizeof( double) );
    // Save k_x - this is necessary because it gets shifted around by Remapping
    final_state_file.write( (char *) fluidEqs->kx_pointer(), sizeof( dcmplx)*(fluidEqs->Cdimxy_full()) );
    // Mean fields
    for (int i=0; i<fluidEqs->num_MFs(); ++i) {
        final_state_file.write( (char *) MF[i].data() , sizeof(dcmplx)*(fluidEqs->MFdimz()) );
    }
    // Ckl - change to single before saving
    Ckl_storage_mat Ckl_single(fluidEqs->Cdimz(), fluidEqs->Cdimz());
    for (int i=0; i<fluidEqs->Cdimxy(); ++i) {
        Ckl_single = Ckl[i].cast<Ckl_storage>();
        final_state_file.write( (char *) Ckl_single.data() , sizeof(Ckl_storage)*(fluidEqs->Cdimz())*(fluidEqs->Cdimz()) );
    }
    final_state_file.close();
}
//
///////////////////////////////////////////////




//////////////////////////////////////////////
//   Load final state

void Load_Full_Data_for_Restart(Model* fluidEqs, Inputs& SP, MPIdata& mpi, dcmplxVec* MF, dcmplxMat* Ckl) {
    // Loads, number of processors (int) and index in k array, then i (int) and t (double) (for kx), then MFs (dcmplx), then all of Ckl (dcmplx)
    // Puts all these in the right place so it can start where it left off
    
    // Check directory exists
    std::string state_dir(SP.simulation_dir+"FINAL_STATE/");
    
    
    // Check saved directory exists
    tinydir_dir tmp_dir;
    int dir_status =tinydir_open(&tmp_dir, state_dir.c_str());
    if ( dir_status == -1 ){
        mpi.print1("Chosen to start from saved state but no FINAL_STATE folder found!\n");
        mpi.print1("Proceding evaluation as specified in input\n");
    } else {
        tinydir_close(&tmp_dir);
            
        std::stringstream final_state_name;
        final_state_name << state_dir << "FINAL_CKL_STATE_Proc" << mpi.my_n_v() << ".dat";
        std::ifstream final_state_file(final_state_name.str().c_str(), std::ios::binary);
        if (!final_state_file.is_open() ) {
            std::cout << "WARNING: " << final_state_name.str() <<  " file unable to be opened" << std::endl;
            ABORT;
        }
        
        /////////////////////////////////////////
        // Read Data
        
        // MPI stuff
        int read_np=-1, read_minxy=-1,read_maxxy=-1; // Read these to check
        // Number of processors, index (on k array) for each processor
        final_state_file.read( (char *) &read_np, sizeof( int ) );
        final_state_file.read( (char *) &read_minxy, sizeof( int ) );
        final_state_file.read( (char *) &read_maxxy, sizeof( int ) );
        
        if (read_np != mpi.total_n_v()) {
            std::cout << "Proc " << mpi.my_n_v() << " ERROR: Number of processors read from FINAL_STATE does not match those in current simulation!" << std::endl;
            ABORT;
        }
        if (read_minxy != mpi.minxy_i() || read_maxxy != mpi.maxxy_i()) {
            std::cout << "Proc " << mpi.my_n_v() << " ERROR: kx, ky range does not match that read from FINAL_STATE!!" << std::endl;
            ABORT;
        }
        
        
        // Time stuff
        int i_curr;
        double t_init;
        // i, just for simplicity
        final_state_file.read( (char *) &i_curr, sizeof( int ) );
        // t, just for simplicity
        final_state_file.read( (char *) &t_init, sizeof( double) );
        // Changes t_start and i_start in SP
        SP.t_start = t_init;
        SP.i_start = i_curr;
    
        
        // kx array - this is necessary because it gets shifted around by Remapping
        dcmplx* kx_point = fluidEqs->kx_pointer();
        final_state_file.read( (char *) kx_point, sizeof( dcmplx)*(fluidEqs->Cdimxy_full()) );
        
        
        // Mean fields
        for (int i=0; i<fluidEqs->num_MFs(); ++i) {
            final_state_file.read( (char *) MF[i].data() , sizeof(dcmplx)*(fluidEqs->MFdimz()) );
        }
        // Ckl
        Ckl_storage_mat Ckl_single(fluidEqs->Cdimz(), fluidEqs->Cdimz());
        for (int i=0; i<fluidEqs->Cdimxy(); ++i) {
            final_state_file.read( (char *) Ckl_single.data() , sizeof(Ckl_storage)*(fluidEqs->Cdimz())*(fluidEqs->Cdimz()) );
            Ckl[i] = Ckl_single.cast<dcmplx>();
            
        }
        final_state_file.close();
        
        std::stringstream printstr;
        printstr << "Loading full solution from disk at t = " << t_init << "\n\n";
        mpi.print1(printstr.str());
    }
}
//
///////////////////////////////////////////////

