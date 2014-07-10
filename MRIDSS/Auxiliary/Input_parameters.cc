//
//  Input_parameters.cc
//  MRIDSS
//
//  Created by Jonathan Squire on 5/7/14.
//  Copyright (c) 2014 Jonathan Squire. All rights reserved.
//

#include "Input_parameters.h"

// Class for handling the input parameters and storing the data.
// This should be a convenient way to input data and have the methods required to process this in a basic way

Inputs::Inputs(const MPIdata& mpi): mpi_node_(mpi.my_n_v()) {
    // Constructor for Inputs class (really a struct)
    
    // Find input file - this is a file in base directory with extension .DSSinput
    std::string file_name = Find_Input_File(CURR_BASE_DIR);
    
    simulation_dir = DATA_DIR+file_name+"/";
    if (mpi_node_ == 0) {
        // Make a new directory in /Data named the same as the input file
        // Check if it exists - create if it does not
        tinydir_dir tmp_dir;
        int dir_status =tinydir_open(&tmp_dir, simulation_dir.c_str());
        if ( dir_status == -1 ){
            std::cout << "Creating directory: " << simulation_dir << std::endl;
            int status = mkdir(simulation_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
            if (status == -1) {
                std::cout << "Warning: Failed to create simulation directory: " << simulation_dir << std::endl;
            }
        }
        tinydir_close(&tmp_dir);
    }
    
    // Reads data from input file

    std::ifstream input_file( (CURR_BASE_DIR + file_name + ".DSSinput").c_str() ); // Full input file
    
    if (input_file.is_open()){
        // Convert file into a string
        std::string fullfile ( std::istreambuf_iterator<char>( input_file ),
                         (std::istreambuf_iterator<char>()) );
        if (mpi_node_ == 0) {
            // Copy file into simulation_dir
            std::ofstream copy_of_input_file( (simulation_dir+file_name+".DSSinput").c_str());
            if (copy_of_input_file.is_open()) {
                copy_of_input_file << fullfile;
            } else {
                std::cout << "Warning: Unable to create copy of input file in " << simulation_dir << std::endl;
            }
            copy_of_input_file.close();
        }
        
        ////////////////////////////////////////////
        // Find each variable from input
        
        
        ////////////////////////////////////////////
        // Read_From_Input_File seems to cause a memory
        // leak when it is reading a long number (e.g., 0.06666666667)
        // However, it is read correctly and the leak is
        // tiny, so I guess will just leave this.
        ////////////////////////////////////////////
        
        // Find each variable in the input file
        // Grid - nx, ny nz, lx,ly, lz
        NXY[0] = Read_From_Input_File_<int>("nx_",fullfile,16);
        NXY[1] = Read_From_Input_File_<int>("ny_",fullfile,16);
        NZ = Read_From_Input_File_<int>("nz_",fullfile,16);
        L[0] = Read_From_Input_File_<double>("lx_",fullfile,2*PI);
        L[1] = Read_From_Input_File_<double>("ly_",fullfile,2*PI);
        L[2] = Read_From_Input_File_<double>("lz_",fullfile,2*PI);
        
        // Time domain variables: dt, t_final, t_initial
        dt = Read_From_Input_File_<double>("dt_",fullfile,0.33333333333333);
        t_final = Read_From_Input_File_<double>("t_final_",fullfile,10.0);
        t_initial = Read_From_Input_File_<double>("t_initial_",fullfile,0.0);
        // Saving: tv_save, full_save
        timvar_save_interval = Read_From_Input_File_<double>("tv_save_",fullfile,dt);
        fullsol_save_interval = Read_From_Input_File_<double>("full_save_",fullfile,0.0);
        if (fullsol_save_interval<t_final && fullsol_save_interval != 0) {
            mpi.print1("Warning: intermediate saves not yet supported!\nSetting fullsol_save_interval to t_final\n");
            fullsol_save_interval = t_final + 1;
        }
        
        // Other physical parameters: nu, eta, q, f_noise
        nu = Read_From_Input_File_<double>("nu_",fullfile,0.01);
        eta = Read_From_Input_File_<double>("eta_",fullfile,0.01);
        q = Read_From_Input_File_<double>("q_",fullfile, 1.5);
        f_noise = Read_From_Input_File_<double>("f_noise_",fullfile, 0);
        
        // Boolean flags: remapQ, QuasiLinearQ
        remapQ = Read_From_Input_File_<bool>("Remap?", fullfile, 1);
        QuasiLinearQ = Read_From_Input_File_<bool>("QuasiLinear?", fullfile, 1);
        
        // Saving time variables
        energy_save_Q = Read_From_Input_File_<bool>("save_energy?", fullfile, 0);
        AM_save_Q = Read_From_Input_File_<bool>("save_angular_mom?", fullfile, 0);
        dissipation_save_Q = Read_From_Input_File_<bool>("save_dissipation?", fullfile, 0);
        reynolds_save_Q = Read_From_Input_File_<bool>("save_reynolds_stress?", fullfile, 0);
        
        mean_field_save_Q = Read_From_Input_File_<bool>("mean_field_save?", fullfile, 1);
        
        // Initial conditions
        // Initial By, multiplying cos(k_z(1) z). If set to negative, uses random noise (multiplied by specified number)
        initial_By = Read_From_Input_File_<double>("initial_By_",fullfile, 1e-7);
        
        // Restarting simulation
        start_from_saved_Q = Read_From_Input_File_<bool>("start_from_saved?", fullfile, 0);
        // If True, takes step, time, and data from FINAL_CKL_STATE_Proc#
        // This is loaded after other initialization, with t_final taken from input

    } else {
        std::cout << "Error: Input file not found!" << std::endl;
        ABORT;
    }
    
    // May want to resave the input file here, once I have a better saving mechanism
    
    // Work out various derived quantities
    initialize_();
}

void Inputs::initialize_() {
    // Time variables
    nsteps = t_final/dt;
    t_start = t_initial;
    i_start = 0;
    timevar_save_nsteps = round(timvar_save_interval/dt);
    fullsol_save_nsteps = round(fullsol_save_interval/dt);
    
    // Shearing box
    num_before_remap = 2*nsteps; // If no remap
    if (remapQ) {
        TSB = L[1]/L[0]/q;
        // Test for an integer number of steps
        if (fabs(TSB/dt-round(TSB/dt)) > 1e-10) {
            if (mpi_node_== 0){
                std::cout << "Warning, non-integer number of steps before remapping: " << TSB/dt <<std::endl;
            }
        }
        num_before_remap = round(TSB/dt);
    }
}



template <class T>
T Inputs::Read_From_Input_File_(const std::string& varstr, const std::string& full_file,
                                T default_value){
    // Reads the variable that matches varstr from the text file (in full_file)
    // should be in format
    //      varstr = val
    // Where val matches the chosen template type
    // default value is used if nothing is found
    
    // Find position of variable
    std::size_t pos = full_file.find(varstr);
    
    if (pos != std::string::npos) {
        
        // Take the part of the string containing this (truncate at 100 characters)
        pos = full_file.find_first_of("=",pos);
        std::string str_end = full_file.substr(pos+1,pos+100);
        // Read into variable
        T out;
        std::istringstream str_end_SS(str_end);
        str_end_SS >> out;
        
        return out;
    } else {
        if (mpi_node_==0) {
            std::cout << "Warning: Nonexistent variable " << varstr << " in input file. ";
            std::cout << "Using default value " << default_value << std::endl;
        }

        return default_value;
    }
}

// Find file in directory that has extensino .DSSinput
std::string Inputs::Find_Input_File(const std::string& directory){
    tinydir_dir dir;
    tinydir_open(&dir, directory.c_str());
    
    std::string inputfile;// Input file name string
    unsigned found_input = 0; // Flag for error checking
    while (dir.has_next)
    {
        tinydir_file file;
        tinydir_readfile(&dir, &file);
        
        std::string filename(file.name);
        
        std::size_t in_ext_pos = filename.find("DSSinput");// Does it contain DSSinput?
        if (in_ext_pos != std::string::npos) {
            // Found file
            inputfile = filename.substr(0,in_ext_pos-1);
            found_input+=1;
        }
        tinydir_next(&dir);
    }
    tinydir_close(&dir);
    
    if (found_input ==1) {
        return inputfile;
    } else if (found_input > 1) {
        if (mpi_node_ ==0) {
            std::cout <<"Warning: Found multiple input files, using first: " << inputfile << ".DSSinput" << std::endl;
        }
        return inputfile;
    } else {
        if (mpi_node_ ==0) {
            std::cout <<"Error: Found no input file (extension .DSSinput)" << std::endl;
        }
        ABORT;
        return 0;
    }

}