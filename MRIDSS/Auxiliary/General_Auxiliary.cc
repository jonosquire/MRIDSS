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
// FRIEND TO THE MODEL CLASS
    dcmplx kxLfac(0,-1/(2*PI/mod->L_[0]) );
    dcmplx kyLfac(0,-1/(2*PI/mod->L_[1]) );
    
    dcmplx * kxp = mod->kx_;// Pointer to kxp data
    for (int i=0; i<mod->nxy_full_; ++i) {
        // Find new kx vector
        kxp[i] = kxp[i]*kxLfac + mod->ky_[i]*kyLfac;
        
        // If kx_ is too large zero out Ckl and set kx negative
        double overlap = kxp[i].real()+1 - mod->Nxy_[0]/2;
        if (overlap>0) {
            Ckl[i].setZero();
            kxp[i] = -mod->Nxy_[0]/2+overlap;
        }
        
        // Un-normalize again
        kxp[i] = kxp[i]/kxLfac;
    }
};