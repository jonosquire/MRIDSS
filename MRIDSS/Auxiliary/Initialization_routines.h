//
//  Initialization_routines.h
//  MRIDSS
//
//  Created by Jonathan Squire on 1/05/14.
//  Copyright (c) 2014 Jonathan Squire. All rights reserved.
//

#ifndef __MRIDSS__Initialization_routines__
#define __MRIDSS__Initialization_routines__

#include "../General_Definitions.h"

#include "../Models/Model.h"
#include "../Integrators/Integrator.h"


// Assigning FFTW plans in the Model Constructor
void createFFTplans(dcmplx* Ckl0, Model* eq);
// Initial conditions
void InitialConditions(dcmplxVec * MF, dcmplxMat* Ckl, int numMF, int NXY);
// Initialize solutions
void Initialize_Solutions(dcmplxVec *MF, dcmplxMat* Ckl,
                          int num_MFs, int MFdimz, int nxy_full, int nz_Cfull);

// Initializing k
void define_kxy_array(dcmplx* kx,dcmplx* ky, int* Nxy, const double* L);
void define_kz_array(dcmplxVec& kz, int NZ, const double* L);

#endif /* defined(__MRIDSS__Initialization_routines__) */
