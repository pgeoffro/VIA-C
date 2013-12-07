#ifndef MODELAWAS_H
#define MODELAWAS_H

#include <iostream>
#include<Eigen/Core>

class ModelAwas
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


    public:
    double stiffness;
    double inertia;
    double Wx;
    double Wu0;
    double Wu1;
    double Wterminal;
    double Wlim;
    double dT;
    double Tdes;
    double mu;
    double alpha;
    int n; //size of the vector during optimization
    double window;

  public:
    ModelAwas();

    virtual ~ModelAwas();
    // tests de RK4
    virtual State_t evolutionT1(const State_t& s, const Control_t& c) const;
    virtual State_t evolutionT2(const State_t& s, const Control_t& c) const;
    // The real one for Model Awas and RK4
    virtual State_t evolutionRK4(const State_t& s, const Control_t& c) const;


    // AwAS Model
    virtual State_t init(const double& theta, const double& r) const;
    virtual void torqueWanted(const double& t) ;


   // Display Functions
    virtual void display14(const State_t& state) const;
    virtual void display12(const State_t& state) const;
    virtual void display22(const State_dx& state) const;
    virtual void display44(const State_dx& sdx) const;
    virtual void display42(const State_du& sdx) const;

    // Derivation Functions
    virtual State_t evolution(State_t& state, const Control_t& control) ;
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

    // Re-Init functions
    virtual Cost_dx initInstCost_dx(const State_t& state) const;
    virtual Cost_dxx initInstCost_dxx(const State_t& state) const;
    virtual State_dx initEvolution_dx   (const State_t& state) const;
    virtual State_du initEvolution_du   () const;

    protected:
    private:
};

#endif // MODELAWAS_H
