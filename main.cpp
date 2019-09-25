//
//  main.cpp
//  Subgradient
//
//  Created by Rui Shibasaki on 25/04/17.
//  Copyright © 2017 Rui Shibasaki. All rights reserved.
//

#include <iostream>
#include <iomanip>
#include "subgradient_g.hpp"
#include "mcnd.hpp"

int main(int argc, const char * argv[]) {
    
    
    
    std::string instance(argv[1]);
    MCND  mcnd_data(instance);
    int primalsize = mcnd_data.narcs * mcnd_data.ndemands + mcnd_data.narcs;
    int dualsize = mcnd_data.nnodes * mcnd_data.ndemands;

    subgrad MCNDsubgrad(primalsize, dualsize, "subgradient.par", instance);

    std::cout<<instance<<" ";

#ifndef WIN32
    // start time measurement
    double t0;
    struct tms timearr; clock_t tres;
    tres = times(&timearr);
    t0 = timearr.tms_utime;
#endif
    
    std::cout<<("\nSTART SOLVING...\n");
    // invoke volume algorithm
    if (MCNDsubgrad.solve(mcnd_data) < 0) {
        std::cout <<"solve failed..."<<std::endl;//printf("solve failed...\n");
    } else {
        
#ifndef WIN32
        // end time measurement
        tres = times(&timearr);
        double t = (timearr.tms_utime-t0)/100.;
#endif
        std::cout<<"Lower bound: "<<std::setprecision(15)<<MCNDsubgrad.lstar<<" stat: ";
        std::cout<<t<<" "; //printf(" Total Time: %f secs\n", t);
        std::cout<<MCNDsubgrad.iter_<<std::endl;  //printf(" Iteration finale: %d\n", volp.iter() );
        
    }
    return 0;
}


