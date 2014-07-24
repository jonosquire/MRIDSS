//
//  General_Definitions.h
//  MRIDSS
//
//  Created by Jonathan Squire on 5/6/14.
//  Copyright (c) 2014 Jonathan Squire. All rights reserved.
//

#ifndef MRIDSS_General_Definitions_h
#define MRIDSS_General_Definitions_h


//////////////////////////////////////////////////
/////              MPI FLAG               ////////
// This is better as a preprocessor directive since it lets you debug in serial more easily

#define USE_MPI_FLAG

//////////////////////////////////////////////////

// Handle for general definitions - all programs should #include "General_Definitions.h"
// Contains necessary libraries, typedefs, global variable declarations etc.

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <sstream>
#include <iomanip>
#include <complex>
//#include <random>
#include <ctime>
#include <cmath>
#include <sys/types.h>
#include <sys/stat.h>


// External libraries
#ifdef USE_MPI_FLAG
    #include <mpi.h>
#endif
#include <fftw3.h>
#define EIGEN_NO_DEBUG  // Define to turn off eigen range checking
#include "../Eigen/Dense"
#include "../External_headers/tinydir.h"



//////////////////////////////////////////////////
// Convenient FFTW definitions

// FFTW planning mechanism
#define MY_FFTWPLAN FFTW_EXHAUSTIVE
// Reinterpret dcmplx for use in fftw routines
#define CAST_T0_FFTW(a) reinterpret_cast<fftw_complex*>(&(a)[0])

//////////////////////////////////////////////////


/////////////////////////////////////////////////////
//  Aborting - useful to have MPI abort. (obviously would be better with execptions)
#ifdef USE_MPI_FLAG
    #define ABORT MPI_Abort(MPI_COMM_WORLD,1);
#else
    #define ABORT abort();
#endif

//////////////////////////////////////////////////
// Eigen matrix and array typedefs

typedef std::complex<double> dcmplx;
typedef Eigen::Matrix<dcmplx,Eigen::Dynamic,Eigen::Dynamic> dcmplxMat;
// USE COLUMN_MAJOR FORMAT FOR MATRIX. THE ONLY PLACE WHERE THIS MIGHT CAUSE A PROBLEM IS IN FFTW, SO LONG AS I'M CAREFUL SHOULD BE EASY
// Reason for this choice is to facilitate transformation over from Matlab
typedef Eigen::Array<dcmplx, Eigen::Dynamic, 1> dcmplxVec;
typedef Eigen::Matrix<dcmplx, Eigen::Dynamic, 1> dcmplxVecM; // vec for eigen matrix object - used for .asDiagonal objects
typedef Eigen::Array<double, Eigen::Dynamic, 1> doubVec;// For use in linear part

// Saving Ckl
typedef std::complex<float> Ckl_storage; // For saving Ckl (kept general for backwards compatibility)
typedef Eigen::Matrix<Ckl_storage, Eigen::Dynamic, Eigen::Dynamic> Ckl_storage_mat;

//////////////////////////////////////////////////


const double PI=atan(1)*4;

const std::string CURR_BASE_DIR = "/Users/jsquire/Documents/MRIDSS/MRIDSS/";
const std::string DATA_DIR = "/Users/jsquire/Documents/MRIDSS/MRIDSS/Data/";
// Data directory is defined from root so that /p/mridss can be used on the cluster

#endif
