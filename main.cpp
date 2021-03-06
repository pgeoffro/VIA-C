#include <iostream>
#include <fstream>
#include <cmath>
#include <list>
#include <Eigen/Core>

#include "include/ModelAwas.h"
#include "include/integratorRK4.h"
#include "include/ilqr.h"
#include "include/Robot3R.h"
#include "include/MathTools.h"
#include "include/ddpls.h"

using namespace std;


int main()
{
    typedef ModelAwas::State_t State_t;
    typedef ModelAwas::Control_t Control_t;
    typedef ilqr<ModelAwas>::StateList_t StateList_t;
    typedef ilqr<ModelAwas>::ControlList_t ControlList_t;

    integratorRK4<ModelAwas> robot;

 //   ModelAwas model;
/* --- Test de RK4 ------------------*/
/*
//  Test RK4 avec Exp
  ModelAwas::State_t state(1);
   state(0) = 1;

 // Test RK4 avec Sin et Cos
    ModelAwas::State_t state(2);
    state(0) = 1;
    ModelAwas::Control_t control(1);

  for (int i=0; i<314/2; i++){
////   robot.integrate(state,control);
    State_t k1 = 0.01 * model.evolutionT2 (state, control);
    State_t k2 = 0.01 * model.evolutionT2 (state + k1*(0.01/2.0), control);
    State_t k3 = 0.01 * model.evolutionT2 (state + k2*(0.01/2.0), control);
    State_t k4 = 0.01 * model.evolutionT2 (state + k3*0.01/2.0, control);

    state = state + 1.0/6.0 *(k1 + 2*k2 + 2*k3 + k4);
    cout << state(0) << " VS " << cos((i+1)*0.01) << " difference "<< state(0) - cos((i+1)*0.01) <<endl;
 }
*/


/* ---- tests des fonctions de ModelAwas */
/*State_t state = model.init(0.05,0.1);
cout << "state :"<< endl;
model.display14(state);

Control_t control(2);
state = model.evolution(state,control);
cout << "state � t+dt, Control nul :"<< endl;
model.display14(state);

control(0) = 0.3;
control(1) = -0.05;
state = model.evolution(state,control);
cout << "state � t+2*dt Control non nul:"<< endl;
model.display14(state);

ModelAwas::State_dx state_dx(4,4);
state_dx = model.evolution_dx(state);
cout << "state dx : "<< endl;
model.display44(state_dx);

ModelAwas::State_du state_du(2,4);
state_du = model.evolution_du();
cout << "state du : "<< endl;
model.display42(state_du);

ModelAwas::Cost_t cost = model.instCost(state,control);
cout << "cost : " << cost << endl;

ModelAwas::Cost_dx cost_dx = model.instCost_dx(state,control);
cout << "cost dx : " << endl;
model.display14(cost_dx);

state(2) = 0.056;
cost_dx = model.instCost_dx(state,control);
cout << "cost dx hors des bornes (r=0.056 m): " << endl;
model.display14(cost_dx);

state(2) = 0.143;
cost_dx = model.instCost_dx(state,control);
cout << "cost dx hors des bornes (r = 0.143): " << endl;
model.display14(cost_dx);

state(2) = 0.08;
ModelAwas::Cost_dxx cost_dxx = model.instCost_dxx(state,control);
cout << "cost dxx  (r = 0.08): " << endl;
model.display44(cost_dxx);

state(2) = 0.056;
cost_dxx = model.instCost_dxx(state,control);
cout << "cost dxx  (r = 0.056): " << endl;
model.display44(cost_dxx);

state(2) = 0.143;
cost_dxx = model.instCost_dxx(state,control);
cout << "cost dxx  (r = 0.143): " << endl;
model.display44(cost_dxx);

ModelAwas::Cost_du cost_du = model.instCost_du(control);
cout << "cost du: " << endl;
model.display12(cost_du);

ModelAwas::Cost_duu cost_duu = model.instCost_duu(control);
cout << "cost duu: " << endl;
model.display22(cost_duu);

state(2) = 0.08;
cost_dx = model.termCost_dx(state,control);
cout << "terminal cost dxx  (r = 0.08): " << endl;
model.display14(cost_dx);

cost_dxx = model.termCost_dxx(state,control);
cout << "terminal cost dxx  (r = 0.08): " << endl;
model.display44(cost_dxx);
*/

/* --- Test de l'int�grateur RK4 -------*/
/*   ModelAwas::State_t state(4);
    state(2) = 0.1;
    state(0) = 0.0;
    state(1) = 0.0;
    state(3) = 0.0;
    ModelAwas::Control_t control(2);
    control(1) = 0; control(0) = 0;
//    integratorRK4<ModelAwas> robot;
  for(int i = 0; i<1000;i++){
 state = robot.integrate(state,control);
    cout << state << endl<<endl;
    }
*/




/* --- Tests de ilqr ------------- */
/*
ilqr<ModelAwas> algo;
algo.init();
//algo.displayVx();
State_t s=model.init(0.05,0.1);
algo.initState(s);
algo.Loops();
//algo.displayV(algo.controlList);
algo.displayV(algo.stateList);
cout<< endl;
algo.displayV(algo.controlList);
*/


/* --- Mise en place de l'algorithme complet -------------*/
/*ilqr<ModelAwas> algo;
State_t s = model.init(0.05,0.1);
algo.completeAlgo(s);

//algo.displayV(algo.torqueList);

ofstream fichier("RTorque.txt", ios::out | ios::trunc);
for(std::list<double>::iterator i=algo.torqueList.begin();i!=algo.torqueList.end();i++){
fichier<<*i;
fichier<<endl;}

ofstream fichier2("HTorque.txt", ios::out | ios::trunc);
for(std::list<double>::iterator i=algo.torqueWanted.begin();i!=algo.torqueWanted.end();i++){
fichier2<<*i;
fichier2<<endl;}
*/

/* --- Tests des fonctions de Robot3R --------------- */
/*
Robot3R r3r;
Robot3R::State_t S(3);
S(0) = 0; S(1) = 0.0; S(2) = 3.14159/2.0;
cout<<S<<endl;
cout<<"xeff : "<<r3r.xeff(S)<<" and yeff : "<<r3r.yeff(S)<<endl;

Robot3R::Control_t C(3);
C(0) = 0; C(1) = 0; C(2) = 0;
cout<<"Control : "<<endl<<C<<endl;
cout<<"Cost_du : "<<endl<<r3r.instCost_du(C)<<endl;
cout<<"Cost_duu : "<<endl<<r3r.instCost_duu(C)<<endl;

cout<<"Cost : "<<endl<<r3r.instCost(S,C)<<endl;
cout<<"Cost_dx : "<<endl<<r3r.instCost_dx(S)<<endl;
cout<<"Cost_dxx : "<<endl<<r3r.instCost_dxx(S)<<endl;

cout<<"TermCost_dx : "<<endl<<r3r.termCost_dx(S)<<endl;
cout<<"TermCost_dxx : "<<endl<<r3r.termCost_dxx(S)<<endl;
*/

/* --- DDP with Robot3R ------------*/
/*Robot3R model;
Robot3R::State_t S(3);
S(0) = 3; S(1) = -0.6; S(2) = 0.0;

integratorRK4<Robot3R> rk4;
ilqr<Robot3R> algo;
algo.completeAlgo(S);

//algo.displayV(algo.stateList);

ofstream fichier1("q1.txt", ios::out | ios::trunc);

for(std::list<Robot3R::State_t>::iterator i=algo.stateHist.begin();i!=algo.stateHist.end();i++){
fichier1<<*i;
fichier1<<endl;
}
*/


/* --- Tests des outils math�matiques ---------------- */
/*
MathTools MT;
/*MathTools::Matrix I = MT.Id(4);
MathTools::Vector A(3); A(0) = 12; A(1)= 6; A(2) = -4;
MathTools::Vector B(3); B(0) = -51; B(1) = 167; B(2) = 24;

/*MathTools::Matrix M(3,3);
M(0,0) = 12; M(0,1) = -51; M(0,2) = 4;
M(1,0) = 6; M(1,1) = 167; M(1,2) = -68;
M(2,0) = -4; M(2,1) = 24; M(2,2) = -41;*/

/*MathTools::Matrix M(3,4);
M(0,0) = 0; M(0,1) = 4; M(0,2) = 3; M(0,3) = 3;
M(1,0) = -1; M(1,1) = -3; M(1,2) = -2; M(1,3) = -3;
M(2,0) = 0; M(2,1) = 0; M(2,2) = 5; M(2,3) = 3;*/

/*
MathTools::Matrix M(4,3);
M(0,0) = 2; M(0,1) = 4; M(0,2) = 0;
M(1,0) = -1; M(1,1) = 0; M(1,2) = 0;
M(2,0) = 3; M(2,1) = 3; M(2,2) = 5;
M(3,0) = 2; M(3,1) = 3; M(3,2) = 3;*/

/*MathTools::Matrix M(10,6);
M(0,0) = -1.492 ;
M(0,1) = -1.351;
M(0,2) = -0.6755;
M(1,0) = -2.4678 ;
M(1,1) = -0.7374 ;
M(1,2) = -0.7374 ;
M(2,3) = 0.01;
M(3,4) = 0.01;
M(4,5) = 0.01;
//MT.Proj(A,A,3);
//MathTools::Matrix Q = MT.GramSchmidtQ(M,3,3);
MathTools::Matrix U = MT.HouseholderQ(M,10,6);
MathTools::Matrix R = MT.GramSchmidtR(M,U,10,6);

//std::cout<<MT.Proj(B,A,3)<<std::endl<<std::endl;
std::cout<<M<<std::endl<<std::endl;
//std::cout<<Q<<std::endl<<std::endl;
std::cout<<"U : "<<std::endl<<U<<std::endl<<std::endl;
std::cout<<"R :"<<std::endl<<R<<std::endl;
std::cout<<"M :"<<std::endl<<U*R<<std::endl;
*/

/* --- Tests des fonctions LS du Robot R3 --------- */
/*MathTools MT;
Robot3R model;
Robot3R::State_t S(3);
S(0) = 3; S(1) = -0.6; S(2) = 0.0;
Robot3R::Control_t C(3);
C(0) = 1; C(1) = 2; C(2) = -1;
//std::cout<<model.instCostLS(S,C)<<std::endl;
//std::cout<<model.instCostLS_dx(S,C)<<std::endl;
//std::cout<<model.instCostLS_du(S,C)<<std::endl;

Robot3R::CostLS_du Rurx = model.instCostLSrurx(S,C);
MathTools::Matrix Q = MT.HouseholderQ(Rurx,5,6);
MathTools::Matrix R = MT.GramSchmidtR(Rurx,Q,5,6);

std::cout<<"Rurx : "<<Rurx<<std::endl;
std::cout<<"Q : "<<Q<<std::endl;
std::cout<<"R : "<<R<<std::endl;
*/

/* --- DDP LS ---------- */
MathTools MT;
Robot3R model;
Robot3R::State_t S(3);
S(0) = 3; S(1) = 1.5; S(2) = 0.0;
Robot3R::Control_t C(3);
C(0) = 1; C(1) = 2; C(2) = -1;

ddpls<Robot3R> algo;
algo.init();
algo.initState(S);
//algo.initR();
//algo.backwardLoop();
//algo.forwardLoop(0);
//algo.Loops();
algo.completeAlgo(S);
ofstream fichier1("ddpls.txt", ios::out | ios::trunc);

for(std::list<Robot3R::State_t>::iterator i=algo.stateHist.begin();i!=algo.stateHist.end();i++){
fichier1<<*i;
fichier1<<endl;
}

ilqr<Robot3R> algo2;
algo2.completeAlgo(S);

//algo.displayV(algo.stateList);

ofstream fichier2("q1.txt", ios::out | ios::trunc);

for(std::list<Robot3R::State_t>::iterator i=algo2.stateHist.begin();i!=algo2.stateHist.end();i++){
fichier2<<*i;
fichier2<<endl;
}
//std::cout<<"G : "<<algo.G<<std::endl;

return 0;
}
