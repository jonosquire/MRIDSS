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
#include "Models/MHD_fullBQlin_noMaxwell.h"
#include "Models/HD_fullU.h"
#include "Models/Constant_Damping.h"
// Integrators
#include "Integrators/intEuler.h"
#include "Integrators/intEulerCN.h"
#include "Integrators/intRK2CN.h"
#include "Integrators/intRK3CN.h"

// Auxiliary
#include "Auxiliary/Input_parameters.h"
#include "Auxiliary/MPIdata.h"
#include "Auxiliary/TimeVariables.h"
#include "Auxiliary/InterfaceOutput.h"


int main(int argc, char *argv[])
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
    std::string input_file_name;
    if (argc == 2) {
        input_file_name = argv[1]; // If an argument is passed, this specifies input file
    } else {
       input_file_name = "null";
    }
    
    Inputs SP( mpi , input_file_name);
    
    ////////////////////////////////////////////////////////////////

    // Construct "model" object
    fftwPlans fft;
    Model* fluidEqs;
    std::stringstream printstr;
    printstr << "Using model: " << SP.equations_to_use << std::endl;;;
    if (SP.equations_to_use == "MHD_BQlin") {
        fluidEqs = new MHD_BQlin(SP, mpi, fft);
    } else if (SP.equations_to_use == "MHD_fullBQlin") {
        fluidEqs = new MHD_fullBQlin(SP, mpi, fft);
    } else if (SP.equations_to_use == "MHD_FullUBQlin") {
        fluidEqs = new MHD_FullUBQlin(SP, mpi, fft);
    } else if (SP.equations_to_use == "MHD_fullBQlin_noMaxwell") {
        fluidEqs = new MHD_fullBQlin_noMaxwell(SP, mpi, fft);
    } else if (SP.equations_to_use == "HD_fullU") {
        fluidEqs = new HD_fullU(SP, mpi, fft);
    } else {
        std::cout << "ERROR: no matching model found!" << std::endl;
        ABORT;
        fluidEqs = new MHD_BQlin(SP, mpi, fft); // Initialize to shut up compiler!
    }
    mpi.print1(printstr.str());

    
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

    double t=SP.t_start ; // Initial time
    double dt; // Time-step, set by integrator (possibly from input file)
    int step_since_TV = 0, step_since_dump = 0; // Counting steps since timevar/dump
    
    // Set up time variables -- energy, engular momentum etc.
    TimeVariables time_vars(SP, fluidEqs->num_linFs(), fluidEqs->num_MFs(), fluidEqs->num_Reynolds_saves(), mpi.my_n_v());
    
    bool QL_effects_are_on=SP.QuasiLinearQ;

    // Set up integration scheme
    Integrator *integrator = new RK3CN(t, SP, *fluidEqs);
    
    // Time variables at zero before we start
    if (!SP.start_from_saved_Q) {
        fluidEqs->Calc_Energy_AM_Diss(time_vars, t, MF, Ckl );
        time_vars.Save_Mean_Fields(MF,  fft);
    }
    
    // Interface output and timing
    InterfaceOutput *interface = new InterfaceOutput(&mpi, &time_vars, integrator, &SP);
    
    // Main loop
    while ( t < SP.t_final && !interface->StopSimulationQ()){
        
        dt = integrator->Step(t, MF, Ckl);
        
        // Update times
        t += dt;
        SP.timevar_t += dt;
        SP.fullsave_t += dt;
        ++step_since_TV;++step_since_dump; // Steps since save
        
        if (SP.remapQ) // Remap at every step now
            ShearingBox_Continuous_Remap(SP.q*t,fluidEqs, Ckl);
        
        
        // Calculate energy, AM, etc.
        if (SP.timevar_t - SP.timvar_save_interval > -1e-8){
            fluidEqs->Calc_Energy_AM_Diss(time_vars, t, MF, Ckl );
            time_vars.Save_Mean_Fields(MF,  fft);
            SP.timevar_t -= SP.timvar_save_interval;
            step_since_TV = 0;
            // Check solution
            CheckSolution(MF, Ckl, fft);
        }
        
//        if (t > 100.0 && !QL_effects_are_on) {
//            fluidEqs->set_QL_YN(1);
//            // Reinitialize integrator to recalculate MF linear operator
//            integrator->Reinitialize_linear_Ops(t);
//
//            QL_effects_are_on=1;
//        }
        
        // Interface
        interface->CheckAll(t);
        if (interface->StopSimulationQ()){
            step_since_TV = -1;
        }
        
    }
    // Final energy and mean fields
    if (step_since_TV != 0){
        fluidEqs->Calc_Energy_AM_Diss(time_vars, t, MF, Ckl );
        time_vars.Save_Mean_Fields(MF,  fft);
    }
    
    
    if (SP.fullsol_save_interval !=0 || interface->StopSimulationQ()) {
        Save_Full_Data_for_Restart(fluidEqs,SP,mpi,t,MF,Ckl);
        mpi.print1("Saved Full Data to disk\n\n");
    }
    
    // Timing
    std::stringstream time_str;
    time_str << "Finished calculation at  t = " << t << " with average time-step " << integrator->mean_time_step() << std::endl << "Full time "<< interface->ElapsedTime() << "s:" << " Time variables: " << time_vars.TVtime() << "s" << std::endl;
    mpi.print1( time_str.str() );
    
    delete interface;
    delete[] MF;
    delete[] Ckl;
    delete integrator;
    delete fluidEqs;
#ifdef USE_MPI_FLAG
    MPI_Finalize();
#endif
    return 0;
}

