#ifndef INTEGRATORRK4_H
#define INTEGRATORRK4_H


template<typename Model_t>
class integratorRK4
{
    public: /* -- Types -- */
    typedef typename Model_t::State_t State_t;
    typedef typename Model_t::Control_t Control_t;

    public :
         Model_t model;

    public:
        integratorRK4(void);
        virtual ~integratorRK4();


    public :
        State_t integrate (const State_t& state, const Control_t& control) const;


    protected:
    private:
};


//#include"../src/integratorRK4.t.cpp"

#endif // INTEGRATORRK4_H
