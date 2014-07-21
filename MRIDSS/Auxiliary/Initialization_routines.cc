//
//  Initialization_routines.cc
//  MRIDSS
//
//  Created by Jonathan Squire on 1/05/14.
//  Copyright (c) 2014 Jonathan Squire. All rights reserved.
//

#include "Initialization_routines.h"



// Initialize the memory for the solutions, MF and Cin
void Initialize_Solutions(dcmplxVec *MF, dcmplxMat* Ckl,Model* fluidEqs) {
    // Take dimensinos from Model
    int MFdimz = fluidEqs->MFdimz();
    int nzC = fluidEqs->Cdimz();
    for (int i=0; i<fluidEqs->num_MFs(); ++i)
        MF[i] = dcmplxVec( MFdimz );
    // Split across processes
    for (int i=0; i < fluidEqs->Cdimxy(); i++)
        Ckl[i]= dcmplxMat(nzC,nzC);
}


// Initial conditions - outputs in Fourier space
void InitialConditions(dcmplxVec * MF, dcmplxMat* Ckl, Inputs SP, Model* equations, MPIdata& mpi, fftwPlans& fft){
//    // Generate random distribution
//    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
//    std::default_random_engine generator(seed);
//    std::normal_distribution<double> distribution(0.0,1.0);
//    auto normal_dist = std::bind ( distribution, generator ); // Asign simple call
    
    int NXY = mpi.nxy();
    int numMF = equations->num_MFs();
    
    /////////////////////////////////
    //      EVERYTHING ZERO        //
    
    // Assign to all entries of Ckl
    int nz = Ckl[0].rows();
    for (int i=0; i<NXY; ++i) {
        for (int j=0; j<nz; ++j) {
            for (int k=0; k<nz; ++k)
                Ckl[i](j,k)=dcmplx(0,0);
        }
        Ckl[i] += Ckl[i].adjoint().eval();
        
    }

//    // Assign to mean fields
//    for (int i=0; i<numMF; ++i) {
//        for (int j=0; j<nz/4; ++j) {
//            MF[i](j) = dcmplx(0,0);
//        }
//    }
//    //                             //
//    /////////////////////////////////
    
//    /////////////////////////////////
//    //   LOAD FROM (MATLAB) FILE   //
//    
//    // DOES NOT WORK WITH MPI - NEED TO UPGRADE
//    std::string init_filename("/Users/jsquire/Documents/MRIDSS/MRIDSS/Data/SW_test/initial_conditions.dat");
//    std::ifstream init_file(init_filename, std::ios::in | std::ios::binary);
//    if (init_file.is_open()){
//        // Read into Ckl
//        long nzf = Ckl[0].rows();
//        double *realbuf = new double[nzf*nzf];
//        double *imagbuf = new double[nzf*nzf];
//        for (int i=0; i<NXY; ++i) {
//            init_file.read((char*)realbuf, sizeof(realbuf)*nzf*nzf);
//            init_file.read((char*)imagbuf, sizeof(imagbuf)*nzf*nzf);
//            
//            dcmplx* Cp = Ckl[i].data();
//            for (int i=0; i<nzf*nzf; i++) {
//                Cp[i] = dcmplx(realbuf[i],imagbuf[i]);
//            }
////            std::cout << Ckl[i] << std::endl;
////            std::cout << std::endl;
//
//        }
//        init_file.close();
//        // End if
//        delete[] realbuf;
//        delete[] imagbuf;
//    }
//    else {
//        std::cout << "ERROR: Initial condition failed to open" << std::endl;
//    }

    // Assign to mean fields
    doubVec zg = doubVec::LinSpaced( MF[0].size(), 0, 2*PI*(1.0-1.0/MF[0].size()) );
    
    // Decide what MF initial conditions based on SP.initial_By
    if (SP.initial_By > 0.0) { // Lowest Kz mode in the box, amplitude from SP.initial_By
        MF[0].setZero();
        for (int i=0; i<MF[0].size(); ++i) {
            MF[1](i) = (dcmplx) SP.initial_By*cos( zg(i) );
        }
    } else { // Sets to random in Bx and By, amplitude from SP.initial_By
        // Sets to random in Bx and By
        double mult_fac[2] = {-SP.initial_By, -0.1*SP.initial_By};
        for (int i=0; i<numMF; ++i) {
            MF[i].real().setRandom();
            MF[i].imag().setZero();
            MF[i] *= mult_fac[i];
        }

    }

    // Take the Fourier transform
    for (int i=0; i<numMF; ++i) {
        fft.for_1D(MF[i].data());
    }


    
}

// INITIALIZING K - THESE ARE BETTER HERE TO REUSE BETWEEN DIFFERENT MODELS

// Defining kx and ky variables in the constructor for EulerFluid
// kx and ky are vectors of complex, k2 is 2-D array of doubles
void define_kxy_array(dcmplx* kx,dcmplx* ky, int* ky_index, int* Nxy, const double* L){
    // Define kx and ky arrays - length is Nx*Ny
    
    // Leave out the (zero,zero) frequncy 

    // Include only positive ky values!
    for (int i=0; i<Nxy[0]/2; ++i) {
        for (int j=0; j<Nxy[1]; ++j) {
            kx[j+Nxy[1]*i]=dcmplx(0,i*2*PI/L[0]);
            ky[j+Nxy[1]*i]=dcmplx(0,j*2*PI/L[1]);
            ky_index[j+Nxy[1]*i]=j;
        }
    }
    // Leave out the Nyquist frequncy in kx
    for (int i=Nxy[0]/2+1; i<Nxy[0]; ++i) {
        for (int j=0; j<Nxy[1]; ++j) {
            kx[j+Nxy[1]*(i-1)]=dcmplx(0,(-Nxy[0]+i)*2*PI/L[0]);
            ky[j+Nxy[1]*(i-1)]=dcmplx(0,j*2*PI/L[1]);
            ky_index[j+Nxy[1]*(i-1)]=j;
        }
    }
}


// Define 1-D kz array
void define_kz_array(dcmplxVec& kz, int NZ, const double* L) {
    // kz is length NZ pointer array
    for (int i=0; i<NZ/2; ++i)
        kz(i) = dcmplx(0,i*2*PI/L[2]);
    for (int i=NZ/2; i<NZ; ++i)
        kz(i) = dcmplx(0,(-NZ+i)*2*PI/L[2]);
    
}









