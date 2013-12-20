//#ifndef "INTEGRATORRK4_H"
//#define "INTERGRATORRK4_H"

#include "../include/integratorRK4.h"
#include "../include/ModelAwas.h"
#include "../include/Robot3R.h"

template<typename Model_t> integratorRK4<Model_t>:: integratorRK4(void)
{

}

template<typename Model_t> integratorRK4<Model_t>::~integratorRK4()
{
    //dtor
}

  template<typename Model_t>
  typename Model_t::State_t integratorRK4<Model_t>::
  integrate (const State_t& state, const Control_t& control) const
  {
    double dt;
    dt = model.dT/40;
    State_t state_new = state;
    //std::cout<<state_new<<std::endl<<std::endl;
    for (int i = 0;i<40;i++){
    State_t k1 = dt*model.evolutionRK4(state_new, control);
    State_t k2 = dt*model.evolutionRK4 (state_new + k1*(dt/2), control);
    State_t k3 = dt*model.evolutionRK4 (state_new+ k2*(dt/2), control);
    State_t k4 = dt*model.evolutionRK4 (state_new + k3*dt/2, control);

    state_new = state_new + 1.0/6.0 *(k1 + 2*k2 + 2*k3 + k4);}
    return state_new;
  }



template class integratorRK4<ModelAwas>;

template class integratorRK4<Robot3R>;

//#endif
