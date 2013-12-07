#include "../include/ModelAwas.h"
#include<Eigen/Core>
#include <iostream>
#include <cmath>

ModelAwas::ModelAwas()
{
    stiffness = 10000.0;
    inertia = (0.27),
    Wx = (2000),
    Wu0 = (0.1),
    Wu1 = (0.0),
    Wterminal = (1100),
    Wlim = (10),
    dT = (0.0005),
    Tdes = (0.0),
    mu = (0.05),
    alpha = (0.01),
    n = (10),
    window = (5);

}

ModelAwas::~ModelAwas()
{
    //dtor
}



/* --- Evolution du système Exponentielle -----------*/
ModelAwas::State_t ModelAwas::evolutionT1(const State_t& s,const Control_t& C) const{
State_t S(1);
S(0) = s(0);
return S;
}
/* --- Evolution du système Sinusoidale -------------- */
ModelAwas::State_t ModelAwas::evolutionT2(const State_t& s,const Control_t& C) const{
State_t S(2);
S(0) = s(1);
S(1) = -s(0);
return S;
}

/* --- Evolution du modèle Awas ------------------------- */
ModelAwas::State_t ModelAwas::evolutionRK4(const State_t& s, const Control_t& c) const{
    State_t S(4);
    S(0) = s(1);
    S(1) = - stiffness/inertia*s(2)*s(2)*sin(s(0)) - c(0);
    S(2) = s(3) ;
    S(3) = c(1);
    return S;
}


/* --- Initialisation Du Model AwAS : vecteur de taille 4, theta dtheta r dr --------- */
ModelAwas::State_t ModelAwas::init(const double& theta, const double& r) const{
State_t state(4);
state(0) = theta;
state(2) = r;
state(1) = 0.0; state(3) = 0.0;
return state;
}

  void ModelAwas::torqueWanted(const double& t) {
  ModelAwas::Tdes = t;
  }

/* --- DISPLAY FUNCTIONS -------------- */
void ModelAwas::display14(const State_t& state) const{
std::cout << state(0) <<std::endl;
std::cout << state(1) <<std::endl;
std::cout << state(2) <<std::endl;
std::cout << state(3) <<std::endl;
return;
}

void ModelAwas::display12(const State_t& state) const{
std::cout << state(0) <<std::endl;
std::cout << state(1) <<std::endl;
return;
}


void ModelAwas::display22(const State_dx& sdx) const{
std::cout << sdx(0,0)<<" ; "<< sdx(0,1) << std::endl;
std::cout << sdx(1,0)<<" ; "<< sdx(1,1) << std::endl;
return;
}

void ModelAwas::display44(const State_dx& sdx) const{
std::cout << sdx(0,0)<<" ; "<< sdx(0,1) << ";" << sdx(0,2) << " ; " << sdx(0,3) <<std::endl;
std::cout << sdx(1,0)<<" ; "<< sdx(1,1) << ";" << sdx(1,2) << " ; " << sdx(1,3) <<std::endl;
std::cout << sdx(2,0)<<" ; "<< sdx(2,1) << ";" << sdx(2,2) << " ; " << sdx(2,3) <<std::endl;
std::cout << sdx(3,0)<<" ; "<< sdx(3,1) << ";" << sdx(3,2) << " ; " << sdx(3,3) <<std::endl;

return;
}

void ModelAwas::display42(const State_du& sdx) const{
std::cout << sdx(0,0)<<" ; "<< sdx(0,1) <<std::endl;
std::cout << sdx(1,0)<<" ; "<< sdx(1,1) <<std::endl;
std::cout << sdx(2,0)<<" ; "<< sdx(2,1) <<std::endl;
std::cout << sdx(3,0)<<" ; "<< sdx(3,1) <<std::endl;

return;
}
  /* --- DERIVATION FUNCTIONS ----------------------------------------------- */
  ModelAwas::State_t ModelAwas::  evolution(State_t& state, const Control_t& control)
  {
    const double theta = state(0)+dT*state(1);
    const double r = state(2)+dT*state(3);
    const double rp = state(3)+dT*control(1);
    const double thetap = state(1)+dT*(-stiffness/inertia*state(2)*state(2)*state(0)-control(0));

    State_t state_new(4);
    state_new(0) = theta;
    state_new(1) = thetap;
    state_new(2) = r;
    state_new(3) = rp;

    return state_new;
  }

  ModelAwas::State_dx ModelAwas::  evolution_dx   (const State_t& state) const
  {
    State_dx state_dx(4,4);
    state_dx(0,0) = 1;
    state_dx(0,1) = dT;
    state_dx(1,0) = -dT*stiffness/inertia*state(2)*state(2);
    state_dx(1,1) = 1;
    state_dx(1,2) = -2*dT*stiffness/inertia*state(0)*state(2);
    state_dx(2,2) = 1;
    state_dx(2,3) = dT;
    state_dx(3,3) = 1;

    state_dx(0,2) = 0 ;
    state_dx(0,3) = 0;
    state_dx(1,3) = 0;
    state_dx(2,0) = 0;
    state_dx(2,1) = 0;
    state_dx(3,0) = 0;
    state_dx(3,1) = 0;
    state_dx(3,2) = 0;

    return state_dx;
  }

  ModelAwas::State_du ModelAwas::  evolution_du   () const
  {
    State_du state_du(4,2);
    state_du(1,0) = -dT;
    state_du(3,1) = dT;

    state_du(0,0) = 0;
    state_du(0,1) = 0;
    state_du(1,1) = 0;
    state_du(2,0) = 0;
    state_du(2,1) = 0;
    state_du(3,0) = 0;

    return state_du;
  }


    /* --- COST FUNCTIONS ------------*/
  ModelAwas::Cost_t ModelAwas::  instCost     (const State_t& state, const Control_t& control) const
  {
   Cost_t cost;
   cost = (stiffness*state(2)*state(2)*sin(state(0))-Tdes) * (stiffness*state(2)*state(2)*sin(state(0))-Tdes);
      return cost;
  }

  ModelAwas::Cost_dx ModelAwas::instCost_dx (const State_t& state) const
{ Cost_dx cost_dx(4);

double r = state(2);
double theta = state(0);

double dLr = -4.0*stiffness*r*sin(theta) *(-stiffness*sin(theta)*r*r + Tdes);
double dLtheta = -2.0*stiffness*r*r*cos(theta)*(- stiffness*sin(theta)*r*r + Tdes);

// Fonction barrière
double DL =0;
if(r<0.06){
    DL = -1.0/0.01*1.0/((r-0.05)/0.01)/((r-0.05)/0.01) -2/0.01*((r-0.05)/0.01) -1/0.01;}
else{
        if(r>0.14){
    DL = 1.0/0.01*1.0/((-r+0.15)/0.01)/((-r+0.15)/0.01) +2/0.01*((-r+0.15)/0.01)+1/0.01;
}
    }

cost_dx (0) = dLtheta;
cost_dx (1) = 0.0;
cost_dx (2) = dLr + Wlim*DL;
cost_dx (3) = 0.0;

return Wx *cost_dx;
}

ModelAwas::Cost_dxx ModelAwas::instCost_dxx (const State_t& state) const
{Cost_dxx cost_dxx(4,4);

double r = state(2);
double theta = state(0);

double dLrr = 8.0*stiffness*stiffness*r*r*sin(theta)*sin(theta) - 4.0*stiffness*sin(theta)*(-stiffness*sin(theta)*r*r + Tdes);
double dLthetatheta = 2.0*stiffness*stiffness*r*r*r*r*cos(theta)*cos(theta) + 2.0* stiffness*r*r*sin(theta)*(-stiffness*sin(theta)*r*r + Tdes);
double dLrtheta = 4.0*stiffness*stiffness*r*r*r*cos(theta)*sin(theta) - 4.0*stiffness*r*cos(theta)*(-stiffness*sin(theta)*r*r + Tdes);

// Fonction barrière
    double DL =0;
    if(r<0.06){
        DL = 2.0/0.01/0.01*1.0/((r-0.05)/0.01)/((r-0.05)/0.01)/((r-0.05)/0.01) -2/0.01/0.01;
        }
    else{
        if(r>0.14){
            DL = -2.0/0.01/0.01*1.0/((-r+0.15)/0.01)/((-r+0.15)/0.01)/((-r+0.15)/0.01) -2/0.01/0.01;
}
    }


cost_dxx(0,0) = dLthetatheta;
cost_dxx(0,1) = 0.0;
cost_dxx(0,2) = dLrtheta;
cost_dxx(0,3) = 0.0;

cost_dxx(1,0) = 0.0;
cost_dxx(1,1) = 0.0;
cost_dxx(1,2) = 0.0;
cost_dxx(1,3) = 0.0;


cost_dxx(2,0) = dLrtheta;
cost_dxx(2,1) = 0.0;
cost_dxx(2,2) = dLrr + Wlim*DL;
cost_dxx(2,3) = 0.0;


cost_dxx(3,0) = 0.0;
cost_dxx(3,1) = 0.0;
cost_dxx(3,2) = 0.0;
cost_dxx(3,3) = 0.0;

return Wx*cost_dxx;
}

ModelAwas::Cost_du ModelAwas::instCost_du (const Control_t& control) const
{Cost_du cost_du(2);
cost_du(0) = Wu0 * 2 * control(0);
cost_du(1) = Wu1 * 2 * control(1);
return cost_du;
}

ModelAwas::Cost_duu ModelAwas::instCost_duu (const Control_t& control) const
{Cost_duu cost_duu(2,2);
cost_duu(0,0) = Wu0 * 2 ;
cost_duu(1,1) = Wu1 * 2 ;
cost_duu(1,0) = 0.0;
cost_duu(0,1) = 0.0;
return cost_duu;
}

ModelAwas::Cost_dx ModelAwas::termCost_dx (const State_t& state) const
{Cost_dx cost_dx = Wterminal * instCost_dx(state);
return cost_dx;
}

ModelAwas::Cost_dxx ModelAwas::termCost_dxx (const State_t& state) const
{Cost_dxx cost_dxx = Wterminal * instCost_dxx(state);
return cost_dxx;
}

/* ---- Functions for ReInitialization ------------------------*/
 ModelAwas::Cost_dx ModelAwas::initInstCost_dx (const State_t& state) const{
 Cost_dx cost_dx(4);

double r = state(2);
double theta = state(0);

double dLtheta = -2.0*stiffness*r*r*cos(theta)*(- stiffness*sin(theta)*r*r + Tdes);

cost_dx (0) = dLtheta;
cost_dx (1) = 0.0;
cost_dx (2) = 0.0;
cost_dx (3) = 0.0;

return cost_dx;
 }

 ModelAwas::Cost_dxx ModelAwas::initInstCost_dxx (const State_t& state) const
{Cost_dxx cost_dxx(4,4);

double r = state(2);
double theta = state(0);

double dLthetatheta = 2.0*stiffness*stiffness*r*r*r*r*cos(theta)*cos(theta) + 2.0* stiffness*r*r*sin(theta)*(-stiffness*sin(theta)*r*r + Tdes);

cost_dxx(0,0) = dLthetatheta;
cost_dxx(0,1) = 0.0;
cost_dxx(0,2) = 0.0;
cost_dxx(0,3) = 0.0;

cost_dxx(1,0) = 0.0;
cost_dxx(1,1) = 0.0;
cost_dxx(1,2) = 0.0;
cost_dxx(1,3) = 0.0;


cost_dxx(2,0) = 0.0;
cost_dxx(2,1) = 0.0;
cost_dxx(2,2) = 0.0;
cost_dxx(2,3) = 0.0;


cost_dxx(3,0) = 0.0;
cost_dxx(3,1) = 0.0;
cost_dxx(3,2) = 0.0;
cost_dxx(3,3) = 0.0;

return cost_dxx;
}

 ModelAwas::State_dx ModelAwas::  initEvolution_dx   (const State_t& state) const  {
    State_dx state_dx(4,4);
    state_dx(0,0) = 1;
    state_dx(0,1) = dT;
    state_dx(1,0) = -dT*stiffness/inertia*state(2)*state(2);
    state_dx(1,1) = 1;
    state_dx(1,2) = 0;
    state_dx(2,2) = 1;
    state_dx(2,3) = 0;
    state_dx(3,3) = 1;

    state_dx(0,2) = 0 ;
    state_dx(0,3) = 0;
    state_dx(1,3) = 0;
    state_dx(2,0) = 0;
    state_dx(2,1) = 0;
    state_dx(3,0) = 0;
    state_dx(3,1) = 0;
    state_dx(3,2) = 0;

    return state_dx;
  }

  ModelAwas::State_du ModelAwas::  initEvolution_du   () const
  {
    State_du state_du(4,2);
    state_du(1,0) = -dT;
    state_du(3,1) = 0;

    state_du(0,0) = 0;
    state_du(0,1) = 0;
    state_du(1,1) = 0;
    state_du(2,0) = 0;
    state_du(2,1) = 0;
    state_du(3,0) = 0;

    return state_du;
  }


