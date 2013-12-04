#include "../include/integratorRK4.h"

template<typename Model_t> integratorRK4<Model_t>::    integratorRK4(void)
{

}

template<typename Model_t> integratorRK4<Model_t>::~integratorRK4()
{
    //dtor
}


  template<typename Model_t>
  typename Model_t::State_t& integratorRK4<Model_t>::
  integrate (const State_t state, const Control_t& control)
  {
    double dt = Model_t::dT/100;
    State_t k1 = dt*evolutionT1 (state, control);
    State_t k2 = dt*evolutionT1 (state + k1*(dt/2), control);
    State_t k3 = dt*evolutionT1 (state + k2*(dt/2), control);
    State_t k4 = dt*evolutionT1 (state + k3*dt/2, control);

    State_t state_new = state + 1.0/6.0 *(k1 + 2*k2 + 2*k3 + k4);
    return state_new;
  }
