//
//  MPIfunctions.h
//  MRIDSS
//
//  Created by Jonathan Squire on 5/7/14.
//  Copyright (c) 2014 Jonathan Squire. All rights reserved.
//

#ifndef __MRIDSS__MPIfunctions__
#define __MRIDSS__MPIfunctions__

#include "../General_Definitions.h"

#endif /* defined(__MRIDSS__MPIfunctions__) */

class MPIdata {
public:
    
    // General MPI data
    int total_nodes;  // Number of processors
    int my_node;     // My node
    int communicator_size; // Something - might want later
    
    // Data for each node in nxy_full array
    // This is x,y dimensions data of Ckl, stored as a single pointer array
    int my_nxyi_min, my_nxyi_max;
    
    
};
