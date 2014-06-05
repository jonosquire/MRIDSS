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
void InitialConditions(dcmplxVec * MF, dcmplxMat* Ckl, int numMF, int NXY){
    // Generate random distribution
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::normal_distribution<double> distribution(0.0,1.0);
//    auto normal_dist = std::bind ( distribution, generator ); // Asign simple call
    
    // Assign to all entries of Ckl
    int nz = Ckl[0].rows();
    for (int i=0; i<NXY; ++i) {
        for (int j=0; j<nz; ++j) {
            for (int k=0; k<nz; ++k)
                Ckl[i](j,k)=dcmplx(0,0);
        }
    }
    
    // Assign to mean fields
    for (int i=0; i<numMF; ++i) {
        for (int j=0; j<nz/4; ++j) {
            MF[i](j) = dcmplx(0,0);
        }
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
// Another version
void define_kxy_array(dcmplx* kx,dcmplx* ky, int* Nxy, const double* L){
    // Define kx and ky arrays - length is Nx*Ny
    
    // Leave out the (zero,zero) frequncy
    
    // Include only positive ky values!
    for (int i=0; i<Nxy[0]/2; ++i) {
        for (int j=0; j<Nxy[1]; ++j) {
            kx[j+Nxy[1]*i]=dcmplx(0,i*2*PI/L[0]);
            ky[j+Nxy[1]*i]=dcmplx(0,j*2*PI/L[1]);
        }
    }
    // Leave out the Nyquist frequncy in kx
    for (int i=Nxy[0]/2+1; i<Nxy[0]; ++i) {
        for (int j=0; j<Nxy[1]; ++j) {
            kx[j+Nxy[1]*(i-1)]=dcmplx(0,(-Nxy[0]+i)*2*PI/L[0]);
            ky[j+Nxy[1]*(i-1)]=dcmplx(0,j*2*PI/L[1]);
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









