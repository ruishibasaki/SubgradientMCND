//
//  mcnd.hpp
//  Subgradient
//
//  Created by Rui Shibasaki on 25/04/17.
//  Copyright © 2017 Rui Shibasaki. All rights reserved.
//

#ifndef mcnd_hpp
#define mcnd_hpp

#include <vector>
#include <list>
#include <queue>
#include <sstream>
#include <cstdio>
#include <iomanip>
#include <string>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <limits>



#include "subgradient_g.hpp"



struct Pair{
    int node, pos;
    Pair &operator=(const Pair &);
    inline Pair(int n,int p):node(n), pos(p) {};
    inline Pair():node(0),pos(-1){};
    inline ~Pair(){};
};



struct Arc{
    int i, j;
    double capa;
    double f;
    double c;
};

struct Demand{
    int O, D;
    int quantity;
};

struct HeapCell{
    int k;
    double rc_;
    
    inline HeapCell(){};
    inline HeapCell(int dk, double rc):k(dk), rc_(rc){}
    HeapCell(const HeapCell & copy){
        this->k = copy.k;
        this->rc_ = copy.rc_;
    }
    
    HeapCell& operator=(HeapCell other){
        this->k = other.k;
        this->rc_= other.rc_;
        return *this;
    }
};

class comp{
    
public:
    bool operator()(const HeapCell & x, const HeapCell & y)const{
        return (x.rc_>y.rc_);
    }
};


class MCND: public user_Provide{
public:
    // for all hooks: return value of -1 means that volume should quit
    // compute reduced costs
    int compute_rc(const std::vector<double> & u, std::vector<double> & rc);
    // solve relaxed problem
    int solve_subproblem(const std::vector<double> & u, const std::vector<double> & rc,
                         double& lcost, std::vector<double> &g, std::vector<double> & x,
                         double& pcost);
    
    int set_UB(double& ub, std::vector<double>& x);
    
    int warm_start(std::vector<double>& u);
    
    int preprocess();
    
    double dijkstra(int source,std::vector<double> & potentials, std::vector<Pair> & preced, const std::vector<double>& cij);
    
    double heuristic(const std::vector<double>& u, const std::vector<double>& x, std::vector<double>& primal);
    
    
public:
    //VOL_dvector dcos; // cost for opening facilities
    std::vector<Arc> arcs;
    std::vector<Demand> d_k;
    std::vector<std::vector<double> >b_ijk;
    int ndemands, narcs, nnodes; // number of demands, number of arcs, number of nodes
    std::list<HeapCell> heap;
    std::vector<std::vector<bool> >rchble;
    std::vector<std::vector<Pair> > neighbors;
    
public:
    MCND() { }
    MCND(std::string fname);
    
    virtual ~MCND(){
        arcs.clear();
        d_k.clear();
        heap.clear();
        for(int k=0;k<ndemands;++k){
            b_ijk[k].clear();
            //rchble[k].clear();
        }
        b_ijk.clear();
       // rchble.clear();
        for (int n=0; n<nnodes; ++n) {
            neighbors[n].clear();
        }
        neighbors.clear();
    }
};



#endif /* mcnd_hpp */
