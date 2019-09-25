//
//  mcnd.cpp
//  Subgradient
//
//  Created by Rui Shibasaki on 25/04/17.
//  Copyright © 2017 Rui Shibasaki. All rights reserved.
//

#include "mcnd.hpp"


int
MCND::compute_rc(const std::vector<double>& u, std::vector<double>& rc)
{
    for(int k=0; k<ndemands; ++k){
        for(int e=0; e<narcs; ++e){
            rc[narcs + k*narcs + e] = arcs[e].c + u[k*nnodes + arcs[e].i-1] - u[k*nnodes + arcs[e].j-1];
        }
    }
    for(int e =0; e<narcs; ++e)
        rc[e] = arcs[e].f;
    
    
    return 0;
}



// IN: dual vector u
// OUT: primal solution to the Lagrangian subproblem (x)
//      optimal value of Lagrangian subproblem (lcost)
//      g = difference between the rhs and lhs when substituting
//                  x into the relaxed constraints (v)
//      objective value of x substituted into the original problem (pcost)
//      xrc
// return value: -1 (volume should quit) 0 ow

int
MCND::solve_subproblem(const std::vector<double>& u, const std::vector<double>& rc,
                       double& lcost, std::vector<double> &g,
                       std::vector<double>& x, double& pcost)
{
    double kpsack;
    double cost_e;
    double fillUp;
    HeapCell cell(0,0.0);
    for(int p=0; p<x.size();++p)
        x[p] =0.0;
    
    lcost =0.0;
    for(int e=0; e<narcs; ++e){
        heap.clear();
        kpsack =0.0;
        fillUp =0.0;
        //get reduced cost for each commodity in arc e
        for(int k=0;k<ndemands;++k){
            if(rc[narcs + k*narcs + e]< 0.0){
                //if(rchble[d_k[k].O-1][arcs[e].i-1]&&
                  // rchble[arcs[e].j-1][d_k[k].D-1]){
                    cell.k = k;
                    cell.rc_ =rc[narcs + k*narcs + e];
                    heap.push_back(cell);
                //}
            }
        }
        
        heap.sort(comp());
        //std::stable_sort(heap.begin(), heap.end(), comp());
        
        
        //solve knap
        while(heap.size()>0 && (fillUp < arcs[e].capa)){
            x[narcs + heap.back().k*narcs + e] = std::min((arcs[e].capa - fillUp),  b_ijk[heap.back().k][e]);
            fillUp += x[narcs +  heap.back().k*narcs + e];
            kpsack += heap.back().rc_ * x[narcs +  heap.back().k*narcs + e];
            heap.pop_back();
        }
        
        //check total reduced cost to set y
        cost_e = kpsack + arcs[e].f;
        if( cost_e < 0.0){
            x[e] = 1.0;	//calcul of y_ij
            lcost += cost_e;
        }
        
    }
    
    //calcul of pcost and x_ij
    pcost= 0;
    for(int e=0; e<narcs;++e){
        pcost += arcs[e].f * x[e];
        for(int k=0; k<ndemands; ++k){
            if(x[e] == 0.0)
                x[narcs + k*narcs + e] = 0.0;
            pcost += arcs[e].c * x[narcs + k*narcs + e];
        }
    }
    
    //calcul of g
    for(int k=0; k<ndemands; ++k){
        lcost += d_k[k].quantity * ( u[k*nnodes + d_k[k].D-1] - u[k*nnodes + d_k[k].O-1]);
        for(int i=0; i<nnodes; ++i){
            
            g[k*nnodes + i] = 0.0;
            
            for(int e =0; e<narcs; ++e){
                if(i == arcs[e].i-1)
                    g[k*nnodes + i] += x[narcs + k*narcs+e];
                else if(i == arcs[e].j-1)
                    g[k*nnodes + i] -= x[narcs + k*narcs+e];
            }
            
            if(i == d_k[k].O-1)
                g[k*nnodes + i] -= d_k[k].quantity;
            else if( i == d_k[k].D-1)
                g[k*nnodes + i] += d_k[k].quantity;
            
        }
    }
    return 0;
}



//--------------------- DIJKSTRA-------------//
Pair &Pair::operator = (const Pair &source){
    node = source.node;
    pos = source.pos;
    return *this;
}

void remove_i(int p, std::list<int> &N){
    
    std::list<int>::iterator it = N.begin();
    std::advance(it,p);
    N.erase(it);
}

Pair arg_min(const std::list<int> & nodes, const std::vector<double> & costs){
    Pair k;
    k.node = nodes.front();
    k.pos = 0;
    std::list<int>::const_iterator it = nodes.begin();
    for (int i=1; i<nodes.size(); ++i){
        if(costs[k.node-1]>costs[*(++it)-1]){
            k.node = *(it);
            k.pos = i;
        }
    }
    return k;
    
}

double MCND::dijkstra(int source, std::vector<double> & potentials, std::vector<Pair> & preced, const std::vector<double>& cij){
    std::list<int> nodes(nnodes);
    std::list<int>::iterator it;
    double pathcost=0.0;
    Pair i;
    
    it = nodes.begin();
    for (int i=0; i<nnodes; ++i) {
        *(it++) = i+1;
        potentials[i] = __DBL_MAX__;
        preced[i].node = 0;
        preced[i].pos = -1;
    }
    
    potentials[source] = 0.0;

    while (nodes.size()>=1) {
        i=arg_min(nodes, potentials);
        remove_i(i.pos, nodes);
        for (int n=0;n<neighbors[i.node-1].size(); ++n) {
            int j = neighbors[i.node-1][n].node;
            int a = neighbors[i.node-1][n].pos;
            if (potentials[j-1]>potentials[i.node-1]+cij[a]){
                potentials[j-1]=potentials[i.node-1]+cij[a];
                preced[j-1].node=i.node;
                preced[j-1].pos=a;
            }
        }
    }
    
    return pathcost;
}

//----------------------------------------------------

double MCND::heuristic(const std::vector<double>& u,const std::vector<double>& x, std::vector<double>& primal){
    
    for(int p=0; p<x.size();++p)
        primal[p] =0.0;

    
    //sort commodities in a decreasing order cijk*(Dk - Ok)
    HeapCell cell(0,0.0);
    heap.clear();
    for(int k=0;k<ndemands;++k){
        cell.k = k;
        cell.rc_ =d_k[k].quantity * ( u[k*nnodes + d_k[k].D-1] - u[k*nnodes + d_k[k].O-1]);
        heap.push_back(cell);
    }
    heap.sort(comp());
    
    
    //main loop
    std::vector<double> wij(narcs,0.0);
    std::vector<double> uij(narcs);
    for (int a=0; a<narcs; ++a) {
        wij[a] = arcs[a].c*(1-x[a]);
        uij[a] = arcs[a].capa;
    }
    
    std::vector<Pair> preced(nnodes);
    std::vector<double> potentials(nnodes,0);
    std::vector<int> path;
    

    int dmand=0;
    double epsP=0.0;
   
    while(!heap.empty()) {
        dmand = heap.front().k;
        heap.pop_front();
        epsP = d_k[dmand].quantity;
        double filled=0.0;
       
        while(filled<d_k[dmand].quantity){
           
            dijkstra(d_k[dmand].O-1,potentials, preced, wij);//find path
            
            //path retrival
            int i,j, arcij;
            j=d_k[dmand].D;
            while (j != d_k[dmand].O){
                i = preced[j-1].node;
                arcij = preced[j-1].pos;
                
                if (uij[arcij]<=1e-10) { //check if impossible path
                    return -1; //no feasible solution found
                }
                
                if(epsP > uij[arcij])
                    epsP = uij[arcij]; //update max flow in path
                
                path.push_back(arcij);
                //std::cout<<j<<" "<<arcij<<std::endl;
                j = i;
            }//std::cout<<j<<std::endl;
            
            //push flow epsP and update capacities
            while(!path.empty()){
                arcij = path.back();
                path.pop_back();
                primal[arcij]=1.0;
                primal[narcs + dmand*narcs+arcij] += epsP;
                uij[arcij]-=epsP;
                
                if(uij[arcij]<=1e-10)//arc is saturated
                    wij[arcij] = __DBL_MAX__; //block
            }
            
            filled+=epsP;
            epsP = d_k[dmand].quantity-filled;
        }
    }
    
    double solvalue=0.0;
    for (int a=0; a<narcs; ++a){
        solvalue+= primal[a]*arcs[a].f;
        for (int k=0; k<ndemands; ++k) {
            solvalue+=primal[narcs + dmand*narcs+a]*arcs[a].c;
            //if(primal[narcs + k*narcs+a]>0.0 || primal[a]>0.0)
              //  std::cout<<a+1<<primal[a]<<" "<<primal[narcs + dmand*narcs+a]<<std::endl;
        }
    }
    return solvalue;
}


//--------------------- HOT START-------------//

int MCND::warm_start(std::vector<double>& u){
    
    //auxiliary vectors
    std::vector<double>rc(narcs);
    std::vector<Pair> preced(nnodes);
    std::vector<double> costs(nnodes,0);
    
    //initialization of auxiliary vectors
    for (int a=0; a<narcs; ++a) {
        //get pondered arc costs
        rc[a] = arcs[a].c + arcs[a].f/arcs[a].capa;
    }
    
    //preprocess(neighbors);
    
    for(int k=0; k<ndemands; ++k){
        
        //reinitialize
        dijkstra(d_k[k].O-1,costs, preced, rc);
        
        for (int n=0; n<nnodes; ++n) {
            u[k*nnodes + n]=costs[n];
        }
    }

    return 0;
}



//------------------------------------------

int MCND::set_UB(double& ub, std::vector<double>& x){
    return 0;
}


int MCND::preprocess(){
   
    rchble.resize(nnodes);
    std::queue<int> fifo;
    int s;
    for(int n=0;n<nnodes;++n){
        rchble[n].resize(nnodes,false);
        
        //BFS - Breadth First Search
        fifo.push(n);
        rchble[n][n] = true;
        while(!fifo.empty()){
            s = fifo.front();
            fifo.pop();
            //std::cout<<"acces "<<s<<std::endl;
            for (int v=0;v<neighbors[s].size(); ++v){
                int j = neighbors[s][v].node-1;
                if(!rchble[n][j]){
                    fifo.push(j);
                    rchble[n][j] = true;
                }
            }
        }
    }
    
    for(int n=0;n<nnodes;++n)
         for(int nn=0;nn<nnodes;++nn)
             if(rchble[n][nn]==false)
                std::cout<<"NOR REACH "<<n+1<<" "<<nn+1<<std::endl;
    std::cout<<"cabo "<<std::endl;
    return 0;
}



MCND::MCND(std::string fname) {
    
    std::ifstream file;
    file.open(fname.c_str());
    if (!file.is_open()) {
        std::cout<<"Failure to open ufl datafile: %s\n "<<fname;
        abort();
    }
    

    std::string s;
    std::istringstream ss;
    
    getline(file,s);
    getline(file,s);
    ss.str(s);

    // read number of locations and number of customers
    
    ss>>nnodes;
    ss>>narcs;
    ss>>ndemands;
    
    ss.clear();
    
    d_k.resize(ndemands);
    arcs.resize(narcs);
    neighbors.resize(nnodes);
    Pair neighb;
    
    for(int i=0; i<narcs; ++i){
        getline(file,s);
        ss.str(s);
        ss>>arcs[i].i;
        ss>>arcs[i].j;
        ss>>arcs[i].c;
        ss>>arcs[i].capa;
        ss>>arcs[i].f;
        neighb.node=arcs[i].j;
        neighb.pos=i;
        neighbors[arcs[i].i-1].push_back(neighb);
        ss.clear();
    }
    
    for(int i=0; i<ndemands; ++i){
        getline(file,s);
        ss.str(s);
        ss>>d_k[i].O >> d_k[i].D >>d_k[i].quantity;
        ss.clear();
    }
    
    
    b_ijk.resize(ndemands);
    //calculate variable bounds
    for(int k=0;k<ndemands;++k){
        b_ijk[k].resize(narcs);
        for(int e=0; e<narcs; ++e){
            b_ijk[k][e] = d_k[k].quantity;
            if(d_k[k].quantity > arcs[e].capa)
                b_ijk[k][e] = arcs[e].capa;
        }
    }
    
    file.close();
    
}
