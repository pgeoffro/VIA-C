#include "../include/ilqr.h"
#include "../include/ModelAwas.h"
#include "../include/integratorRK4.h"

#include<Eigen/Core>
#include<Eigen/LU>
#include<cmath>

template<typename Model_t>
  ilqr<Model_t>::
  ilqr (void)
     {
         nbPreviewSteps = 10;
         NTotal = 50;
  }

/* --- Display values of list<VectorXd> ----------- */
template<typename Model_t> void ilqr<Model_t>::displayV(std::list<VectorXd> L){
for (std::list<VectorXd>::iterator iter=L.begin(); iter!=L.end();++iter){
        std::cout<<*iter<<std::endl;
        std::cout<<std::endl;
}
}

/* --- Display values of list<double> ----------- */
template<typename Model_t> void ilqr<Model_t>::displayV(std::list<double> L){
for (std::list<double>::iterator iter=L.begin(); iter!=L.end();++iter){
        std::cout<<*iter<<std::endl;
        std::cout<<std::endl;
}
}

/* --- The complete DDP algo - the only one used in main with displayV----------------*/
template<typename Model_t> void ilqr<Model_t>::completeAlgo(const State_t& state){
    init();
    initState(state);
    for(int i = 0; i<NTotal;i++){
    Loops();
    computeControl();
    switching();
    extend();
    initLoops();
    }
    }

/* --- Create Lists with correct sizes ----------- */
template<typename Model_t> void
    ilqr<Model_t>::init(){

    assert( model.n>0);

    controlList.resize(model.n-1);
    stateList.resize(model.n);
    gainList.resize(model.n-1);
    vxList.resize(model.n);
    vxxList.resize(model.n);
    openLoopList.resize(model.n-1);
    stateHist.resize(0);
    torqueList.resize(0);

    isInit = true;

    }

/* --- Initialzation of StateList with the same previous state everywhere & ControlList with null command ---- */
template<typename Model_t> void
    ilqr<Model_t>::initState(const State_t& state){
        StateList_t::iterator iterS = stateList.begin();
        *iterS = state;
        State_t s = state;

        Control_t C(2); C(0) = 0.0; C(1) = 0.0;
        for(ControlList_t::iterator iterC = controlList.begin() ; iterC!=controlList.end() ; iterC++ ){
            *iterC = C;
            iterS ++;
             //s = model.evolution(s,C);
            *iterS = model.evolution(s,C);
        }

    }

/* --- Use the computed control to calculate new states (RK4) --------------*/
    template<typename Model_t>
  void ilqr<Model_t>::
  computeControl ()
  {
    assert(isInit);
    integratorRK4<Model_t> robot;

    std::list<VectorXd>::iterator iterState = stateList.begin();
    State_t state = *iterState;

    for (ControlList_t::iterator iterControl=controlList.begin(); iterControl!=controlList.end();++iterControl){
        State_t state = robot.integrate(*iterState,*iterControl);
        iterState++;
        *iterState = state;
    }
    };

/* --- Erase states and controls computed for the window and stock them in historical lists  ----------------*/
template<typename Model_t> void ilqr<Model_t>::switching(){
    double Torque;
    State_t S;
    for (int i = 0; i< model.window;i++){
            S = *stateList.begin();
            stateHist.insert(stateHist.end(),*stateList.begin());
            stateList.erase(stateList.begin());
            controlList.erase(controlList.begin());
            Torque = model.stiffness * S(2) * S(2) *sin(S(0));
            torqueList.insert(torqueList.end(),Torque);
    }
}

/* --- Resize state and control Lists after switching function ----------------*/
template<typename Model_t> void ilqr<Model_t>::extend(){
    StateList_t::iterator iterState=stateList.begin();
    State_t S = *iterState;

    controlList.resize(model.n-1);
    stateList.resize(model.n);

    initState(S);

   /* ControlList_t::iterator iterControl=controlList.end(); iterControl--;
    StateList_t::iterator iterState=stateList.end(); iterState--;
    State_t S = *iterState;
    integratorRK4<Model_t> robot;

    Control_t C(2);
    C(0)= 0.0; C(1) =0.0;

    //controlList.insert(controlList.end(),model.window,*iterControl);

    controlList.insert(controlList.end(),model.window,C);

    for (int i = 0; i< model.window;i++){
    S = robot.integrate(S,C); //*iterControl);
    stateList.insert(stateList.end(),S);
    }*/
}

/* --- Backward Loop ----------- */
template<typename Model_t> void ilqr<Model_t>::backwardLoop () {
    StateList_t::iterator iterState = stateList.end(); iterState--;
    VxList::iterator iterVx = vxList.end();  iterVx--;
    VxxList::iterator iterVxx = vxxList.end(); iterVxx--;
    OpenLoopTermList::iterator iterOpenLoop = openLoopList.end();
    GainList::iterator iterGain = gainList.end();

    State_t state = *iterState;

	// Final State
	Cost_dx vx = model.termCost_dx(state);
    *iterVx= vx; //iterVx--;
	Cost_dxx vxx = model.termCost_dxx(state);
	*iterVxx = vxx; //iterVxx--;

/*
std::cout<<"state : "<<std::endl;
std::cout<<state<<std::endl;

std::cout<<"vx : "<<std::endl;
std::cout<<vx<<std::endl;

std::cout<<"vxx : "<<std::endl;
std::cout<<vxx<<std::endl;
*/
for (ControlList_t::reverse_iterator iterControl=controlList.rbegin(); iterControl!=controlList.rend();++iterControl){
    iterState --;
    iterVx--;
    iterVxx --;
    iterOpenLoop--;
    iterGain--;

    Control_t control = *iterControl;

    State_dx fx = model.evolution_dx(state);
	State_du fu = model. evolution_du();

	Cost_dx Lx = model.instCost_dx(state);
	Cost_dxx Lxx = model.instCost_dxx(state);
	Cost_du Lu = model.instCost_du(control);
	Cost_duu Luu = model.instCost_duu(control);

    for (int i=0;i<4;i++){
	  vxx(i,i) = vxx(i,i)+ model.mu;}

	VectorXd Qx = Lx + fx.transpose()*vx;
	VectorXd Qu = Lu + fu.transpose()*vx;
	MatrixXd Qxx = Lxx + fx.transpose()*vxx*fx;
	MatrixXd Quu = Luu + fu.transpose()*vxx*fu;
	MatrixXd Qux = fu.transpose()*vxx*fx;

	MatrixXd QuuInv = Quu.inverse();
	Quu(0,1) =  Quu(0,1);
	Quu(1,0) =  Quu(1,0);

    //std::cout<<*iterGain<<std::endl<<std::endl;

    *iterOpenLoop =-(QuuInv*Qu);
	*iterGain = -(QuuInv*Qux);
	*iterVx = Qx - (Qux.transpose()*QuuInv*Qu);
	*iterVxx = Qxx - (Qux.transpose()*QuuInv*Qux);

   // std::cout<<*iterGain<<std::endl<<std::endl;

    vx = *iterVx;
    vxx = *iterVxx;



}
}

/* --- Forward Loop ----------- */
template<typename Model_t> void ilqr<Model_t>::forwardLoop (const double alpha) {
    StateList_t::iterator iterState = stateList.begin();
    OpenLoopTermList::iterator iterOpenLoop = openLoopList.begin();
    GainList::iterator iterGain = gainList.begin();

    State_t state = *iterState;
    State_t state_new = *iterState;

for (ControlList_t::iterator iterControl=controlList.begin(); iterControl!=controlList.end();++iterControl){
    iterState++;
    *iterControl = *iterControl + alpha*(*iterOpenLoop) + *iterGain*(state_new - state) ;
       // std::cout<<*iterControl<<std::endl<<std::endl;
    state = *iterState;
    *iterState = model.evolution(state_new,*iterControl);
    state_new = *iterState;
//    Uti = U(:,i) + alpha*k(:,i) + K(:,4*i-3:4*i)*(Xt(:,i) - X(:,i));

    iterOpenLoop++;
    iterGain++;
}

}


/* --- Back- and For-ward Loops computed nbPreviousSteps times  ----------- */
template<typename Model_t> void ilqr<Model_t>::Loops(){
for(int i=0;i<nbPreviewSteps;i++){
backwardLoop();
forwardLoop(model.alpha);
}
}

/* --- Initialization before new iteration --------------- */
template<typename Model_t> void ilqr<Model_t>::initLoops(){
for(int i=0;i<15;i++){
initBackwardLoop();
//std::cout<< "1 : \n";
//displayV(openLoopList);
forwardLoop(0.1);
//std::cout<< "2 : \n";
//displayV(stateList);
}

}

/* --- Backward Loop ----------- */
template<typename Model_t> void ilqr<Model_t>::initBackwardLoop () {
    StateList_t::iterator iterState = stateList.end(); iterState--;
    VxList::iterator iterVx = vxList.end();  iterVx--;
    VxxList::iterator iterVxx = vxxList.end(); iterVxx--;
    OpenLoopTermList::iterator iterOpenLoop = openLoopList.end();
    GainList::iterator iterGain = gainList.end();

    State_t state = *iterState;

	// Final State
	Cost_dx vx = 1000.0*model.initInstCost_dx(state);
    *iterVx= vx; //iterVx--;
   //std::cout<<"vx : "<<std::endl<<vx<<std::endl;
	Cost_dxx vxx = 1000.0 *model.initInstCost_dxx(state);
	*iterVxx = vxx; //iterVxx--;
    //std::cout<<"vxx : "<<std::endl<<vxx<<std::endl;

for (ControlList_t::reverse_iterator iterControl=controlList.rbegin(); iterControl!=controlList.rend();++iterControl){
    iterState --;
    iterVx--;
    iterVxx --;
    iterOpenLoop--;
    iterGain--;

    Control_t control = *iterControl;

    State_dx fx = model.initEvolution_dx(state); //std::cout<<"fx : "<<std::endl<<fx<<std::endl;
	State_du fu = model. initEvolution_du(); //std::cout<<"fu : "<<std::endl<<fu<<std::endl;

	Cost_dx Lx = model.initInstCost_dx(state); //std::cout<<"Lx : "<<std::endl<<Lx<<std::endl;
	Cost_dxx Lxx = model.initInstCost_dxx(state); //std::cout<<"Lxx : "<<std::endl<<Lxx<<std::endl;

    for (int i=0;i<4;i++){
	  vxx(i,i) = vxx(i,i)+ model.mu;}

	VectorXd Qx = Lx + fx.transpose()*vx;
	VectorXd Qu =  fu.transpose()*vx;
	MatrixXd Qxx = Lxx + fx.transpose()*vxx*fx;
	MatrixXd Quu = fu.transpose()*vxx*fu;
	MatrixXd Qux = fu.transpose()*vxx*fx;

    double QInv =  1.0/Quu(0,0);
	MatrixXd QuuInv(2,2); //= Quu.inverse(); std::cout<<"quu-1 : "<<std::endl<<QInv<<std::endl;
	QuuInv(0,1) =  0.0;
	QuuInv(1,0) =  0.0;
	QuuInv(0,0) = QInv;
	QuuInv(1,1) = 1.0;

    *iterOpenLoop =-(QuuInv*Qu);
	*iterGain = -(QuuInv*Qux);
	*iterVx = Qx - (Qux.transpose()*QuuInv*Qu);
	*iterVxx = Qxx - (Qux.transpose()*QuuInv*Qux);

    vx = *iterVx; //std::cout<<"vx : "<<std::endl<<vx<<std::endl;
    vxx = *iterVxx;



}
}


template class ilqr<ModelAwas>;







