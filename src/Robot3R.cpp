#include "../include/Robot3R.h"

#include<Eigen/Core>
#include<Eigen/LU>
#include<cmath>

Robot3R::Robot3R()
{
    L = 1.0;
    dT = 0.001;

    Wterminal = 100;
    Winst = 1;
    Wu = 0.01;

    xdes = 2.0;
    ydes = -1.3;

    mu = 0.05;
    alpha = 0.1;
}

Robot3R::~Robot3R()
{
    //dtor
}

   /* --- Position of effector ---------------------------*/

double Robot3R::xeff(const State_t& state)const{
return L*(cos(state(0)) + cos(state(0)+state(1)) + cos(state(0)+state(1)+state(2)));
}

double Robot3R::dxeff_q0(const State_t& state)const{
return -L*(sin(state(0)) + sin(state(0)+state(1)) + sin(state(0)+state(1)+state(2)));
}

double Robot3R::dxeff_q1(const State_t& state)const{
return -L*(sin(state(0)+state(1)) + sin(state(0)+state(1)+state(2)));
}

double Robot3R::dxeff_q2(const State_t& state)const{
return -L*(sin(state(0)+state(1)+state(2)));
}

double Robot3R::ddxeff_q0(const State_t& state)const{
return -L*(cos(state(0)) + cos(state(0)+state(1)) + cos(state(0)+state(1)+state(2)));
}

double Robot3R::ddxeff_q1(const State_t& state)const{
return -L*(cos(state(0)+state(1)) + cos(state(0)+state(1)+state(2)));
}

double Robot3R::ddxeff_q2(const State_t& state)const{
return -L*(cos(state(0)+state(1)+state(2)));
}

double Robot3R::ddxeff_q0q1(const State_t& state)const{
return -L*(cos(state(0)+state(1)) + cos(state(0)+state(1)+state(2)));
}

double Robot3R::ddxeff_q1q2(const State_t& state)const{
return -L*(cos(state(0)+state(1)+state(2)));
}

double Robot3R::ddxeff_q0q2(const State_t& state)const{
return -L*(cos(state(0)+state(1)+state(2)));
}


/* --- yeff --------------- */

double Robot3R::yeff(const State_t& state)const{
return L*(sin(state(0)) + sin(state(0)+state(1)) + sin(state(0)+state(1)+state(2)));
}

double Robot3R::dyeff_q0(const State_t& state)const{
return L*(cos(state(0)) + cos(state(0)+state(1)) + cos(state(0)+state(1)+state(2)));
}

double Robot3R::dyeff_q1(const State_t& state)const{
return L*(cos(state(0)+state(1)) + cos(state(0)+state(1)+state(2)));
}

double Robot3R::dyeff_q2(const State_t& state)const{
return L*(cos(state(0)+state(1)+state(2)));
}

double Robot3R::ddyeff_q0(const State_t& state)const{
return -L*(sin(state(0)) + sin(state(0)+state(1)) + sin(state(0)+state(1)+state(2)));
}

double Robot3R::ddyeff_q1(const State_t& state)const{
return -L*(sin(state(0)+state(1)) + sin(state(0)+state(1)+state(2)));
}

double Robot3R::ddyeff_q2(const State_t& state)const{
return -L*(sin(state(0)+state(1)+state(2)));
}

double Robot3R::ddyeff_q0q1(const State_t& state)const{
return -L*(sin(state(0)+state(1)) + sin(state(0)+state(1)+state(2)));
}

double Robot3R::ddyeff_q0q2(const State_t& state)const{
return -L*(sin(state(0)+state(1)+state(2)));
}

double Robot3R::ddyeff_q1q2(const State_t& state)const{
return -L*(sin(state(0)+state(1)+state(2)));
}

/* --- Evolution Functions ------------- */

Robot3R::State_t Robot3R::evolution(State_t& state, const Control_t& control){
    State_t S(3);
    S(0) = state(0) + dT*control(0);
    S(1) = state(1) + dT*control(1);
    S(2) = state(2) + dT*control(2);
    return S;
}

Robot3R::State_t Robot3R::evolutionRK4(const State_t& state, const Control_t& control)const{
    State_t S(3);
    S(0) = control(0);
    S(1) = control(1);
    S(2) = control(2);
    return S;
}

 Robot3R::State_dx Robot3R::evolution_dx   (const State_t& state) const{
    State_dx S(3,3);
    S(0,0) = 1; S(1,1) = 1; S(2,2) = 1;
    S(0,1) = 0; S(1,0) = 0;
    S(1,2) = 0; S(2,1) = 0;
    S(2,0) = 0; S(0,2) = 0;

    return S;
    }

 Robot3R::State_du Robot3R::evolution_du   () const{
     State_du S(3,3);
     S(0,0) = dT; S(1,1) = dT; S(2,2) = dT;
    S(0,1) = 0; S(1,0) = 0;
    S(1,2) = 0; S(2,1) = 0;
    S(2,0) = 0; S(0,2) = 0;
     return S;
 }

/* --- Cost functions --------------- */
Robot3R::Cost_t Robot3R::instCost(const State_t& state, const Control_t& control) const{
Cost_t C;
C = 0.0;
C = Winst*((xeff(state) - xdes)*(xeff(state) - xdes) + (yeff(state) - ydes)*(yeff(state) - ydes)) + Wu*(control(0)*control(0) + control(1)*control(1) + control(2)*control(2));
return C;
}

Robot3R::Cost_dx Robot3R::instCost_dx(const State_t& state) const{
    Cost_dx C(3);
    C(0) = 2*Winst*(xeff(state)-xdes)*dxeff_q0(state) + 2*Winst*(yeff(state)-ydes)*dyeff_q0(state);
    C(1) = 2*Winst*(xeff(state)-xdes)*dxeff_q1(state) + 2*Winst*(yeff(state)-ydes)*dyeff_q1(state);
    C(2) = 2*Winst*(xeff(state)-xdes)*dxeff_q2(state) + 2*Winst*(yeff(state)-ydes)*dyeff_q2(state);
    return C/2;
}

Robot3R::Cost_dxx Robot3R::instCost_dxx(const State_t& state) const{
    Cost_dxx C(3,3);
 // Normal
    C(0,0) = 2*Winst*(dxeff_q0(state)*dxeff_q0(state) + ddxeff_q0(state)*(xeff(state)-xdes) + dyeff_q0(state)*dyeff_q0(state) + ddyeff_q0(state)*(yeff(state)-ydes));
    C(1,1) = 2*Winst*(dxeff_q1(state)*dxeff_q1(state) + ddxeff_q1(state)*(xeff(state)-xdes) + dyeff_q1(state)*dyeff_q1(state) + ddyeff_q1(state)*(yeff(state)-ydes));
    C(2,2) = 2*Winst*(dxeff_q2(state)*dxeff_q2(state) + ddxeff_q2(state)*(xeff(state)-xdes) + dyeff_q2(state)*dyeff_q2(state) + ddyeff_q2(state)*(yeff(state)-ydes));

    C(0,1) = 2*Winst*(dxeff_q0(state)*dxeff_q1(state) + (xeff(state)-xdes)*ddxeff_q0q1(state) + dyeff_q0(state)*dyeff_q1(state) + ddyeff_q0q1(state)*(yeff(state)-ydes));
    C(1,0) = C(0,1);

    C(1,2) = 2*Winst*(dxeff_q2(state)*dxeff_q1(state) + (xeff(state)-xdes)*ddxeff_q1q2(state) + dyeff_q2(state)*dyeff_q1(state) + ddyeff_q1q2(state)*(yeff(state)-ydes));
    C(2,1) = C(1,2);

    C(0,2) = 2*Winst*(dxeff_q0(state)*dxeff_q2(state) + (xeff(state)-xdes)*ddxeff_q0q2(state) + dyeff_q0(state)*dyeff_q2(state) + ddyeff_q0q2(state)*(yeff(state)-ydes));
    C(2,0) = C(0,2);

    // Hess = J'J
 /*   Cost_dx J = instCost_dx(state);
   C = J*J.transpose(); */

    return C/2;
}

Robot3R::Cost_du Robot3R::instCost_du(const Control_t& control) const{
    Cost_du C;
    C = 2*Wu*control;
    return C/2;
}

Robot3R::Cost_duu Robot3R::instCost_duu(const Control_t& control) const{
    Cost_duu C(3,3);
    C(0,0)=2*Wu; C(1,1)= 2*Wu; C(2,2) = 2*Wu;
    C(1,0) = 0; C(0,1) = 0; C(2,0) = 0; C(0,2) = 0; C(2,1) = 0; C(1,2) = 0;

     // Hess = J'J
/*         Cost_du J = instCost_du(control);
         C = J*J.transpose();
*/
    return C/2;
}

Robot3R::Cost_dx Robot3R::termCost_dx(const State_t& state) const{
    Cost_dx C;
    C = Wterminal/Winst*instCost_dx(state);
    return C;
}

Robot3R::Cost_dxx Robot3R::termCost_dxx(const State_t& state) const{
    Cost_dxx C;
    C = Wterminal/Winst*instCost_dxx(state);
    return C;
}

/* --- FUNSTIONS FOR LEAST SQUARE DDP ----------------------- */

/* --- Cost --------------------*/

Robot3R::CostLS Robot3R::instCostLS(const State_t& state, const Control_t& control) const{
CostLS C(5);
C(0) = Winst*((xeff(state) - xdes));
C(1) = Winst*(yeff(state) - ydes);
C(2) = Wu *control(0);
C(3) = Wu *control(1);
C(4) = Wu *control(2);

return C;
}

Robot3R::CostLS_dx Robot3R::instCostLS_dx(const State_t& state, const Control_t& control) const{
CostLS_dx C(3,5);
C(0,0) = Winst*dxeff_q0(state); C(1,0) = Winst*dxeff_q1(state); C(2,0) = Winst*dxeff_q2(state);
C(0,1) = Winst*dyeff_q0(state); C(1,1) = Winst*dyeff_q1(state); C(2,1) = Winst*dyeff_q2(state);
C(0,2) = 0.0; C(1,2) = 0.0; C(2,2) = 0.0;
C(0,3) = 0.0; C(1,3) = 0.0; C(2,3) = 0.0;
C(0,4) = 0.0; C(1,4) = 0.0; C(2,4) = 0.0;

return C.transpose();
}

Robot3R::CostLS_du Robot3R::instCostLS_du(const State_t& state, const Control_t& control) const{
CostLS_du C(3,5);
C(0,0) = 0.0; C(1,0) = 0.0; C(2,0) = 0.0;
C(0,1) = 0.0; C(1,1) = 0.0; C(2,1) = 0.0;
C(0,2) = Wu; C(1,2) = 0.0; C(2,2) = 0.0;
C(0,3) = 0.0; C(1,3) = Wu; C(2,3) = 0.0;
C(0,4) = 0.0; C(1,4) = 0.0; C(2,4) = Wu;

return C.transpose();
}

Robot3R::CostLS_du Robot3R::instCostLSrurx(const State_t& state, const Control_t& control) const{
CostLS_du C(5,6);

CostLS_du Cu = instCostLS_du(state,control);
CostLS_dx Cx = instCostLS_dx(state,control);

for(int i=0;i<3;i++){
    for (int j=0;j<5;j++){
        C(j,i) = Cx(j,i);
        C(j,i+3) = Cu(j,i);
    }
}
return C;
}


Robot3R::CostLS Robot3R::termCostLS(const State_t& state) const{
CostLS C(5);
C(0) = Wterminal*((xeff(state) - xdes));
C(1) = Wterminal*(yeff(state) - ydes);
C(2) = 0.0;
C(3) = 0.0;
C(4) = 0.0;

return C;
}

Robot3R::CostLS_dx Robot3R::termCostLS_dx(const State_t& state) const{
CostLS_dx C(3,5);
C(0,0) = Wterminal*dxeff_q0(state); C(1,0) = Wterminal*dxeff_q1(state); C(2,0) = Wterminal*dxeff_q2(state);
C(0,1) = Wterminal*dyeff_q0(state); C(1,1) = Wterminal*dyeff_q1(state); C(2,1) = Wterminal*dyeff_q2(state);
C(0,2) = 0.0; C(1,2) = 0.0; C(2,2) = 0.0;
C(0,3) = 0.0; C(1,3) = 0.0; C(2,3) = 0.0;
C(0,4) = 0.0; C(1,4) = 0.0; C(2,4) = 0.0;

return C.transpose();
}

