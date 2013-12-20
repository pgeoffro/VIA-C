#ifndef ROBOT3R_H
#define ROBOT3R_H

#include <iostream>
#include<Eigen/Core>

class Robot3R
{
    public :
    typedef Eigen::VectorXd State_t;
    typedef Eigen::VectorXd Control_t;
    typedef double Cost_t;

    typedef Eigen::MatrixXd State_dx;
    typedef Eigen::MatrixXd State_du;
    typedef Eigen::VectorXd Cost_dx;
    typedef Eigen::VectorXd Cost_du;
    typedef Eigen::MatrixXd Cost_dxx;
    typedef Eigen::MatrixXd Cost_duu;

    public :
        double L;
        double dT;

        double Wterminal;
        double Winst;
        double Wu;

        double xdes;
        double ydes;

        double mu;
        double alpha;

    public:
        Robot3R();
        virtual ~Robot3R();

        /* --- Position of effector ---------------------------*/

    virtual double xeff(const State_t& state) const;
    virtual double dxeff_q0(const State_t& state) const;
    virtual double dxeff_q1(const State_t& state) const;
    virtual double dxeff_q2(const State_t& state) const;

    virtual double ddxeff_q0(const State_t& state) const;
    virtual double ddxeff_q1(const State_t& state) const;
    virtual double ddxeff_q2(const State_t& state) const;

    virtual double ddxeff_q0q1(const State_t& state) const;
    virtual double ddxeff_q1q2(const State_t& state) const;
    virtual double ddxeff_q0q2(const State_t& state) const;



    virtual double yeff(const State_t& state) const;
    virtual double dyeff_q0(const State_t& state) const;
    virtual double dyeff_q1(const State_t& state) const;
    virtual double dyeff_q2(const State_t& state) const;

    virtual double ddyeff_q0(const State_t& state) const;
    virtual double ddyeff_q1(const State_t& state) const;
    virtual double ddyeff_q2(const State_t& state) const;

    virtual double ddyeff_q0q1(const State_t& state) const;
    virtual double ddyeff_q1q2(const State_t& state) const;
    virtual double ddyeff_q0q2(const State_t& state) const;



    /* --- Evolution Functions ------------- */
    virtual State_t evolution(State_t& state, const Control_t& control) ;
    virtual State_t evolutionRK4(const State_t& state, const Control_t& control) const;

    virtual State_dx evolution_dx   (const State_t& state) const;
    virtual State_du evolution_du   () const;

    // Cost Functions
    virtual Cost_t instCost(const State_t& state, const Control_t& control) const;
    virtual Cost_dx instCost_dx(const State_t& state) const;
    virtual Cost_dxx instCost_dxx(const State_t& state) const;
    virtual Cost_du instCost_du(const Control_t& control) const;
    virtual Cost_duu instCost_duu(const Control_t& control) const;
    virtual Cost_dx termCost_dx(const State_t& state) const;
    virtual Cost_dxx termCost_dxx(const State_t& state) const;


    protected:
    private:
};

#endif // ROBOT3R_H
