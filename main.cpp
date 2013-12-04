#include <iostream>
#include <fstream>
#include <cmath>
#include <list>
#include <Eigen/Core>
#include <io.h>


#include "include/ModelAwas.h"
#include "include/integratorRK4.h"
#include "include/ilqr.h"

#define SLEEP_LGTH 2  // sleep time in seconds
#define NPOINTS    50 // length of array

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__TOS_WIN__)
 #include <conio.h>   //for getch(), needed in wait_for_key()
 #include <windows.h> //for Sleep()
 void sleep(int i) { Sleep(i*1000); }
#endif

using namespace std;


int main()
{
    typedef ModelAwas::State_t State_t;
    typedef ModelAwas::Control_t Control_t;
    typedef ilqr<ModelAwas>::StateList_t StateList_t;
    typedef ilqr<ModelAwas>::ControlList_t ControlList_t;

    integratorRK4<ModelAwas> robot;

    ModelAwas model;
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
    control(1) = 8; control(0) = 3;
    integratorRK4<ModelAwas> robot;
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
ilqr<ModelAwas> algo;
State_t s = model.init(0.05,0.1);
algo.completeAlgo(s);

algo.displayV(algo.torqueList);

ofstream fichier("test.txt", ios::out | ios::trunc);
for(std::list<double>::iterator i=algo.torqueList.begin();i!=algo.torqueList.end();i++){
fichier<<*i;
fichier<<endl;}


/*
Gnuplot::set_GNUPlotPath("C:/Program Files (x86)/gnuplot/bin");

system("cd C:/Program Files (x86)/gnuplot/bin/ && gnuplot.exe macro.txt ");
system("pause");

try
    {
        Gnuplot g1("lines");



    }
    catch (GnuplotException ge)
    {
        cout << ge.what() << endl;
    }
*/
return 0;
}
