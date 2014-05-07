//
//  main.cpp
//  MRIDSS
//
//  Created by Jonathan Squire on 1/05/14.
//  Copyright (c) 2014 Jonathan Squire. All rights reserved.
//

#include "General_Definitions.h"

// Models
#include "Models/MHD_BQlin.h"
// Integrators
#include "Integrators/intEuler.h"
#include "Integrators/intEulerCN.h"
// Auxiliary
#include "Auxiliary/Input_parameters.h"

int main(int argc, char ** argv)
{
    ////////////////////////////////////////////////////////////////
    /////////                                           ////////////
    /////////                INPUTS                     ////////////
    /////////                                           ////////////
    ////////////////////////////////////////////////////////////////
    
    Inputs SP(
        // Viscosity and resistivity
        0.001, 0.001, // (nu,eta)
        // Shear parameter
        1.5,     // q
        // Quasi-linear forcing
        1,      // f_noise
        // Grid resolution
        4, 4, 4,   // (NX, NY, NZ)
        // Box dimensions
        2*PI, 2*PI, 2*PI,     //  (LX,LY,LZ)
        // Time interval, time-step
        0, 10,      //  t_initial, t_final
        0.1,        // dt
        // Time intervals for saving
        1,      // timevar_save_interval -- for mean fields, energy, AM etc.
        10,     // fullsol_save_interval -- saving full solution
        // Remap (Y or N)
        1,       //remapQ
        // Quasi-linear feedback (Y or N)
        1       // QuasiLinearQ
                    );
    MPI_Init(&argc,&argv);

    // Construct "model" object
    Model* fluidEqs = new MHD_BQlin(SP);
    const int nxy_full = fluidEqs->Cdimxy();   // number of matrices
    const int nz_Cfull = fluidEqs->Cdimz();   // matrix size
    const int num_MFs = fluidEqs->num_MFs();   // Number of mean fields
    const int MFdimz = fluidEqs->MFdimz();   // Dimension of mean fields (NZ)
    
    // Initialize solution
    dcmplxVec *MF = new dcmplxVec[num_MFs];
    dcmplxMat *Ckl =new dcmplxMat[nxy_full];
    Initialize_Solutions(MF, Ckl,
                              num_MFs, MFdimz, nxy_full, nz_Cfull);
    
    
    //  Energy, angular momentum, etc.
//    double **energy = new double*[nsteps];
//    for (int i=0; i<nsteps; ++i) {
//        energy[i] = new double[4];
//    }
    
    
    // Create the FFT plans (should be done before assigning the initial conditions)
    createFFTplans(Ckl[0].data(), fluidEqs);
    
    // Intitial conditions 
    InitialConditions(MF, Ckl,  num_MFs, nxy_full);

    // Set up integration scheme
    double t=SP.t_initial;
    Integrator *integrator = new EulerCN(t, SP.dt, *fluidEqs);
    
    // Main loop
    for (int i = 1; i < SP.nsteps+1; i++) {
        integrator->Step(t, MF, Ckl);
        t = i * SP.dt;
        
        
        // Calculate energy, AM, etc.
        if (i%SP.timevar_save_nsteps==0){
            std::cout << "Done step " << i << std::endl;
        }
        
        if (i%SP.num_before_remap == 0) {
            ShearingBox_Remap(fluidEqs, Ckl);
            std::cout << "<<<< Remapping >>>>" << std::endl;
        }
        
        
        // Full save of all data
        if (i%SP.fullsol_save_nsteps ==0) {
            std::cout << "Saving data, step " << i << std::endl;
        };
        
    }

//    for (int i=0; i<nxy_full; ++i) {
//        std::cout << "Array " << i << std::endl;
//        std::cout << Ckl[i] << std::endl;
//    }
//    for (int i=0; i<num_MFs; ++i) {
//        std::cout << "MF " << i << std::endl;
//        std::cout << MF[i] << std::endl;
//    }
    
    MPI_Finalize();
    return 0;
}

