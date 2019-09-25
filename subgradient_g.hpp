//
//  subgradient.hpp
//  Subgradient
//
//  Created by Rui Shibasaki on 25/04/17.
//  Copyright © 2017 Rui Shibasaki. All rights reserved.
//

#ifndef subgradient_g_hpp
#define subgradient_g_hpp

#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <sstream>
#include <cmath>
#include <limits>
#include <fstream>

#ifndef WIN32
#include <sys/times.h>
#include <mach/mach.h>
#endif


class user_Provide;

//CLASS SUBGRADIENT METHOD PARAMETERS---------
class subgrad_parms {
public:
    unsigned int printflag;
    unsigned int printinvl;
    unsigned int greentestinvl;
    unsigned int yellowtestinvl;
    unsigned int redtestinvl;
    double lambdainit;
    double lambdamin;
    double greenfact;
    double yellowfact;
    double redfact;
    
    unsigned int maxsgriters;
    unsigned int  ascent_first_check;
    unsigned int ascent_check_invl;
    double primal_abs_precision;
    double minimum_rel_ascent;
    double maxtime;
    
    subgrad_parms(){}
    ~subgrad_parms() {}
};

//MAIN CLASS SUBGRADIENT METHOD---------
class subgrad{
public:    
    // ATTRIBUTES -----------
    int psize;
    int dsize;
    int iter_;
    unsigned int contgreen, contred;
    
    double viol;
    double lcost;
    double pcost;
    double lstar;
    double ub;
    double lambda;
    
    std::vector<double> x;
    std::vector<double> g;


    std::vector<double> rc; // reduced costs
    std::vector<double> u; // dual vector
    subgrad_parms parm;
    //std::ofstream outfile;
    
#ifndef WIN32
    struct task_basic_info t_info; //mem usage
    mach_msg_type_number_t t_info_count;
    double maxMem;
#endif
    
    //INITIALIZATIONS, CONSTRUCTORS AND DESTRUCTOR -----------

    void read_parms(const char* filename);

    subgrad(int primalsize, int dualsize, const char* filename, std::string instance);
    
    inline ~subgrad() {x.clear(); rc.clear();u.clear(); }
    
    // MAIN FUNCTIONS -----------------------------------
    
    int solve(user_Provide & hooks);
    double compute_stepsize();
    int take_step();
    int update_lambda(unsigned int itertype);
    
    // AUXILIARY FUNCTIONS-------------------------------
    double compute_norm(const std::vector<double> & v)const; // ||.||^2
    double vecDotProduct(const std::vector<double> & v1, const std::vector<double> & v2) const;
    int vecScalarProd(std::vector<double> & v1, double a);
    double sg_abs(const double &);
    double sg_max(const double &,const double &);


    int sumVec(std::vector<double>& w,const std::vector<double> & v1, const std::vector<double> & v2);
    double dotProduct(const std::vector<double> & v1, const std::vector<double> & v2);

    double readjust_target(double oldtarget);
    void print_info();
    
    
};


/** The user hooks should be overridden by the user to provide the
 problem specific routines for the volume algorithm. The user
 should derive a class ...
 
 for all hooks: return value of -1 means that subgrad should quit
 */
class user_Provide {
public:
    virtual ~user_Provide() {}
public:
    /** compute reduced costs
     @param u (IN) the dual variables
     @param rc (OUT) the reduced cost with respect to the dual values
     */
    virtual int compute_rc(const std::vector<double>& u, std::vector<double>& rc) = 0;
    
    /** Solve the subproblem for the subgradient step.
     @param u (IN) the dual variables
     @param rc (IN) the reduced cost with respect to the dual values
     @param g (OUT) the subgradient (constraints violation)
     @param lcost (OUT) the lagrangean cost with respect to the dual values
     @param x (OUT) the primal result of solving the subproblem
     @param pcost (OUT) the primal objective value of <code>x</code>
     */
    virtual int solve_subproblem(const std::vector<double>& u, const std::vector<double>& rc,
                                 double& lcost,std::vector<double> &g, std::vector<double>& x,
                                 double& pcost) = 0;
    
    /** Solve an heuristic for UB computing.
     @param x (IN) current solution variables
     @param ub (OUT) the upper bound set
     */
    virtual int set_UB(double& ub, std::vector<double>& x) = 0;
    
    /** Solve an heuristic for UB computing.
    @param u (OUT) starting dual vector
    */
    virtual int warm_start(std::vector<double>& u) = 0;

};



#endif /* subgradient_hpp */
