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
#include "Models/MHD_fullBQlin.h"
#include "Models/MHD_FullUBQlin.h"
//#include "Models/Constant_Damping.h"
// Integrators
#include "Integrators/intEuler.h"
#include "Integrators/intEulerCN.h"
#include "Integrators/intRK2CN.h"
// Auxiliary
#include "Auxiliary/Input_parameters.h"
#include "Auxiliary/MPIdata.h"
#include "Auxiliary/TimeVariables.h"


int main(int argc, char ** argv)
{
    MPIdata mpi; // Storage of MPI data
#ifdef USE_MPI_FLAG
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, mpi.total_n_p());
    MPI_Comm_rank(MPI_COMM_WORLD, mpi.my_n_p());
    MPI_Comm_size(MPI_COMM_WORLD, mpi.comm_size_p());
#endif
    
    ////////////////////////////////////////////////////////////////
    /////////                                           ////////////
    /////////                INPUTS                     ////////////
    /////////                                           ////////////
    ////////////////////////////////////////////////////////////////
    
    Inputs SP( mpi   );
    
    ////////////////////////////////////////////////////////////////

    // Construct "model" object
    fftwPlans fft;
    Model* fluidEqs = new MHD_FullUBQlin(SP, mpi, fft);
    const int num_MFs = fluidEqs->num_MFs();
    
    // Initialize solution
    dcmplxVec *MF = new dcmplxVec[ num_MFs ];
    dcmplxMat *Ckl =new dcmplxMat[ fluidEqs->Cdimxy() ]; // NB: Cdimxy is size on each proc
    Initialize_Solutions(MF, Ckl, fluidEqs);

    
    // Intitial conditions
    InitialConditions(MF, Ckl, SP, fluidEqs,  mpi, fft); //can be changed by reading from FINAL_STATE

    
    // Load from FINAL_STATE if specified in input file
    if (SP.start_from_saved_Q) {
        Load_Full_Data_for_Restart(fluidEqs, SP, mpi, MF, Ckl);
    }

    // Set up time variables -- energy, engular momentum etc.
    TimeVariables time_vars(SP, 4, fluidEqs->num_MFs(), 5, mpi.my_n_v());
    
    double t=SP.t_start ; // Initial time
    bool QL_effects_are_on=0;

    // Set up integration scheme
    Integrator *integrator = new RK2CN(t, SP.dt, *fluidEqs);
    
    clock_t start = clock();
    double diff;
    
    // Time variables at zero before we start
    if (!SP.start_from_saved_Q) {
        fluidEqs->Calc_Energy_AM_Diss(time_vars, t, MF, Ckl );
        time_vars.Save_Mean_Fields(MF,  fft);
    }
    
    // Main loop
    for (int i = SP.i_start+1; i < SP.nsteps+1; i++) {
        integrator->Step(t, MF, Ckl);
        
        if (SP.remapQ) // Remap at every step now
            ShearingBox_Continuous_Remap(SP.q*t,fluidEqs, Ckl);
        
        // Calculate energy, AM, etc.
        if (i%SP.timevar_save_nsteps==0){
            std::stringstream done_step;
            done_step << "Done step " << i << std::endl;
            mpi.print1( done_step.str() );
            
            fluidEqs->Calc_Energy_AM_Diss(time_vars, t, MF, Ckl );
            time_vars.Save_Mean_Fields(MF,  fft);
            CheckSolution(MF, Ckl, fft);
        }
        
        
        // Full save of all data
        if (i%SP.fullsol_save_nsteps ==0) {
        //    std::cout << "Saving data, step " << i << std::endl;
        };
        
        t = t + SP.dt;
        
        if (t > 100.0 && !QL_effects_are_on) {
            fluidEqs->set_QL_YN(1);
            // Reinitialize integrator to recalculate MF linear operator
            integrator->Reinitialize_linear_Ops(t);

            QL_effects_are_on=1;
        }
        
    }

    
    diff = (clock() - start ) / (double)CLOCKS_PER_SEC;
    std::stringstream time_str;
    time_str <<"Full time "<< diff << "s" << std::endl;
    mpi.print1( time_str.str() );

    
    if (SP.fullsol_save_interval !=0) {
        Save_Full_Data_for_Restart(fluidEqs,SP,mpi,t,MF,Ckl);
        mpi.print1("Saved Full Data to disk\n");
    }
    
    delete[] MF;
    delete[] Ckl;
    delete integrator;
    delete fluidEqs;
#ifdef USE_MPI_FLAG
    MPI_Finalize();
#endif
    return 0;
}

