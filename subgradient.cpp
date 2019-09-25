//
//  subgradient.cpp
//  Subgradient
//
//  Created by Rui Shibasaki on 25/04/17.
//  Copyright © 2017 Rui Shibasaki. All rights reserved.
//


#include "subgradient.hpp"


// AUXILIARY FUNCTIONS-------------------------------

void subgrad::print_info(){
    
    std::cout<<std::setprecision(15)<<iter_<<" L*= "<<lstar<<" target= "<<ub<<" L= "<<lcost<<" iter g,y,r = {"<<contgreen<<" , "
    <<contyellow<<" , "<<contred<<"}"<<std::endl;
}

double subgrad::sg_abs(const double & f)const{
    if (f<0.000) {
        return -f;
    }else return f;
}

double subgrad::sg_max(const double & a,const double & b)const {
    if (a<b) {
        return b;
    }else return a;
}



double subgrad::compute_norm(const std::vector<double> & v) const{
    double sum=0;
    std::vector<double>::const_iterator it;
    for(it=v.begin();it!=v.end();++it){
        sum+= (*it)*(*it);
    }
    return sum;
}

double subgrad::vecDotProduct(const std::vector<double> & v1, const std::vector<double> & v2) const{
    double sum=0;
    if (v1.size()!=v2.size()) {
        std::cout<<"ERROR: (subgrad::vecDotProduct) vectors with different sizes"<<std::endl;
        return -1;
    }
    std::vector<double>::const_iterator it1=v1.begin();
    std::vector<double>::const_iterator it2=v2.begin();
    for(int p=0;p<v1.size();++p){
        sum+= (*(it1))*(*(it2));
        ++it1;++it2;
    }
    return sum;
}

int subgrad::vecScalarProd(std::vector<double> & v1, double a){
    
    std::vector<double>::iterator it=v1.begin();
    for(;it!=v1.end();++it){
        *it *= a;
    }
    
    return 1;
}


int subgrad::sumVec(std::vector<double>& w,const std::vector<double> & v1, const std::vector<double> & v2, double scalar){
    if (v1.size()!=v2.size() || v1.size()!= w.size()) {
        std::cout<<"ERROR: (subgrad::sumVec) vectors with different sizes"<<std::endl;
        return 0;
    }
    std::vector<double>::iterator itw=w.begin();
    std::vector<double>::const_iterator it1=v1.begin();
    std::vector<double>::const_iterator it2=v2.begin();

    for(int p=0;p<v1.size();++p){
        (*(itw)) = (*(it1))+scalar*(*(it2));
        ++itw;++it1;++it2;
    }

    return 1;
}


double subgrad::compute_teta(const double & ddotg) const{
    double teta=0;
    
    if (ddotg < -1e-10) {
        teta = compute_norm(g)/compute_norm(d);
    }else teta =0.00;
    
    return  std::sqrt(teta);
}


double subgrad::readjust_target(double oldtarget){
    double target = ub;
    if (lstar >= target - sg_abs(target) * 0.5) {
        if (sg_abs(lstar) < 10.0) {
            target = 10.0;
        } else {
            target += 0.025 * sg_abs(target);
            target = sg_max(target, lstar + 0.05 * sg_abs(lstar));
        }
        if (target != oldtarget && (parm.printflag & 2)) {
            printf("     **** readjusting target!!! new target = %f *****\n",
                   target);
        }
    }
    return target;
}



// MAIN FUNCTIONS -----------------------------------

double subgrad::compute_stepsize(){
    double stepsize=0;
    double gdotd = vecDotProduct(g, d);

    if(sg_abs(gdotd) > 1e-10)
        stepsize = lambda*(ub - lstar)/gdotd;
    else stepsize = lambda*(ub - lstar)/compute_norm(g);
    return stepsize;
}


int subgrad::take_step(const double & teta){
    std::vector<double>::iterator itu = u.begin();
    std::vector<double>::const_iterator itd = d.begin();
    double stepsize=0;
    if(sumVec(d, g, d,teta)){
        stepsize = compute_stepsize();
        for (int p=0; p<dsize; ++p){
            (*(itu)) += stepsize*(*(itd));
            ++itu;++itd;
        }
        return 0;
    }else return -1;
    
}

int subgrad::update_lambda(unsigned int itertype){
    
    switch (itertype) {
        case 1:
            if (contgreen>=parm.greentestinvl) {
                lambda *= parm.greenfact;
                contgreen=0;
            }
            break;
            
        case 2:
            if (contyellow>=parm.yellowtestinvl) {
                lambda *= parm.yellowfact;
                contyellow=0;
            }
            break;
            
        case 3:
            if ((contred>=parm.redtestinvl) && (lambda>parm.lambdamin)){
                lambda *= parm.redfact;
                contred=0;
            }
            break;
            
        default:
            break;
    }
    
    return 0;
}




int subgrad::solve(user_Provide & hooks){
    
#ifndef WIN32
    // start time measurement
    double t0;
    struct tms timearr; clock_t tres;
    tres = times(&timearr);
    t0 = timearr.tms_utime;
#endif

    int retval = 0;
    retval = hooks.warm_start(u); // compute reduced costs
    if (retval < 0)  return -1;
    retval = hooks.compute_rc(u, rc); // compute reduced costs
    if (retval < 0)  return -1;
    retval = hooks.solve_subproblem(u, rc, lcost, g, x, pcost);
    if (retval < 0)  return -1;
    double target = -__DBL_MAX__/2;
    target = readjust_target(target);
    
    iter_ = 0;
    unsigned int itertype=0;
    double teta = 0.0;
    double ddotg=0.0;

    double * lcost_sequence = new double[parm.ascent_check_invl];
    const int ascent_first_check = std::max(parm.ascent_first_check, parm.ascent_check_invl);
    
    for (iter_ = 1; iter_ <= parm.maxsgriters; ++iter_) {  // main iteration
        
        
        teta = compute_teta(ddotg);
        retval = take_step(teta); // take a dual step
        if (retval < 0)  break;
        /*for (int p=0; p<dsize; ++p) {
            std::cout <<u[p]<<std::endl;
        }*/
        
        retval = hooks.compute_rc(u, rc); // compute reduced costs
        if (retval < 0)  break;
        
        retval = hooks.solve_subproblem(u, rc, lcost, g, x, pcost);  // solve relaxed problem
        if (retval < 0)  break;
       
        retval = hooks.set_UB(ub, x);  // set the UB target
        if (retval < 0)  break;
        
        if(ub>0.0 && ub>pcost){
            ub = pcost;
        }
        
        target = readjust_target(target);
        ub = target;
        ddotg = vecDotProduct(g, d);
        
        if (lcost > lstar) {
            lstar = lcost;
            if(ddotg<0){
                contyellow++;
                contgreen=0;
                contred=0;
                itertype = 2;
            }
            else{
                contgreen++;
                contred=0;
                contyellow=0;
                itertype = 1;
            }
        }else{
            contred++;
            contgreen=0;
            contyellow=0;
            itertype = 3;
        }
        
        retval = update_lambda(itertype);
        if (retval < 0)  break;
        
        // test terminating criteria
        const bool primal_feas = (std::sqrt(compute_norm(g)) < parm.primal_abs_precision);
        
        // test optimality
        if (primal_feas){
            if (parm.printflag) printf(" small violation \n");
            break;
        }
        
        // test for non-improvement
        const int k = iter_ % parm.ascent_check_invl;
        if (iter_ > ascent_first_check) {
            if (lcost - lcost_sequence[k] < sg_abs(lcost_sequence[k]) * parm.minimum_rel_ascent){
                if (parm.printflag) printf(" small improvement \n");
                break;
            }
        }
        lcost_sequence[k] = lcost;
        
#ifndef WIN32
        // end time measurement
        tres = times(&timearr);
        double t = (timearr.tms_utime-t0)/100.;
        
        if(iter_%parm.printinvl == 0){
            //print_info();
            //outfile<<std::setprecision(15)<<lstar<<" "<<t<<std::endl;
        }
        
        if (t > parm.maxtime) {
            break;
        }
        
        //memcheck
       /* t_info_count = TASK_BASIC_INFO_COUNT;
        if(KERN_SUCCESS != task_info(mach_task_self(), TASK_BASIC_INFO, (task_info_t)&t_info,&t_info_count)){
            std::cout<<"PROBLEM WITH MEM USAGE CHECK"<<std::endl;
        }
        
        if (t_info.resident_size > maxMem) {
            maxMem = t_info.resident_size;
        }*/

#endif
    }
    delete[] lcost_sequence;
    
    return retval;
}



//INITIALIZATIONS, CONSTRUCTORS AND DESTRUCTOR -----------

subgrad::subgrad(int primalsize, int dualsize, const char* filename, std::string instance){
    psize = primalsize;
    dsize = dualsize;
    read_parms(filename);
    ub=10;
    viol=0;
    lcost=0;
    pcost =0;
    lstar = 0;
    contred=0;
    contgreen=0;
    contyellow=0;
    g.resize(dsize, 0.0);
    d.resize(dsize, 0.0);
    x.resize(psize, 0.0);
    rc.resize(psize, 0.0); // reduced costs
    u.resize(dsize, 0.0); // dual vector
    lambda = parm.lambdainit;
    const char *file = (instance.substr(instance.size()-14,1)+instance.substr(instance.size()-10,6)).c_str();
    //outfile.open(file);
}


void subgrad::read_parms(const char *filename){
    char s[100];
    FILE * infile = fopen(filename, "r");
    if (!infile) {
        printf("Failure to open ufl datafile: %s\n ", filename);
        std::abort();
    }
    
    while (fgets(s, 100, infile)) {
        unsigned long len = strlen(s) - 1;
        if (s[len] == '\n')
            s[len] = 0;
        std::string ss(s);
        
        if (ss.find("printflag") == 0) {
            unsigned long i = ss.find("=");
            parm.printflag = atoi(&s[i+1]);
            
        } else if (ss.find("printinvl") == 0) {
            unsigned long i = ss.find("=");
            parm.printinvl = atoi(&s[i+1]);
            
        } else if (ss.find("maxsgriters") == 0) {
            unsigned long i = ss.find("=");
            parm.maxsgriters = atoi(&s[i+1]);
            
        } else if (ss.find("greentestinvl") == 0) {
            unsigned long i = ss.find("=");
            parm.greentestinvl = atoi(&s[i+1]);
            
        } else if (ss.find("yellowtestinvl") == 0) {
            unsigned long i = ss.find("=");
            parm.yellowtestinvl = atoi(&s[i+1]);
            
        } else if (ss.find("redtestinvl") == 0) {
            unsigned long i = ss.find("=");
            parm.redtestinvl = atoi(&s[i+1]);
            
        } else if (ss.find("greenfact") == 0) {
            unsigned long i = ss.find("=");
            parm.greenfact = atof(&s[i+1]);
            
        } else if (ss.find("yellowfact") == 0) {
            unsigned long i = ss.find("=");
            parm.yellowfact = atof(&s[i+1]);
            
        } else if (ss.find("redfact") == 0) {
            unsigned long i = ss.find("=");
            parm.redfact = atof(&s[i+1]);
            
        } else if (ss.find("lambdainit") == 0) {
            unsigned long i = ss.find("=");
            parm.lambdainit = atof(&s[i+1]);
            
        } else if (ss.find("lambdamin") == 0) {
            unsigned long i = ss.find("=");
            parm.lambdamin = atof(&s[i+1]);
            
        } else if (ss.find("ascent_first_check") == 0) {
            unsigned long i = ss.find("=");
            parm.ascent_first_check = atoi(&s[i+1]);
            
        } else if (ss.find("ascent_check_invl") == 0) {
            unsigned long i = ss.find("=");
            parm.ascent_check_invl = atoi(&s[i+1]);
            
        } else if (ss.find("primal_abs_precision") == 0) {
            unsigned long i = ss.find("=");
            parm.primal_abs_precision = atof(&s[i+1]);
            
        } else if (ss.find("minimum_rel_ascent") == 0) {
            unsigned long i = ss.find("=");
            parm.minimum_rel_ascent = atof(&s[i+1]);
            
        } else if (ss.find("maxtime") == 0) {
            unsigned long i = ss.find("=");
            parm.maxtime = atof(&s[i+1]);
        }
        
    }
    fclose(infile);
}



