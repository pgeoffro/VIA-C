#ifndef ILQR_H
#define ILQR_H

#include<Eigen/Core>
#include<list>

template<typename Model_t>
class ilqr
{
      public: /* -- Types -- */
    typedef typename Model_t::State_t State_t;
    typedef typename Model_t::Control_t Control_t;
    typedef typename Model_t::State_dx State_dx;
    typedef typename Model_t::State_du State_du;
    typedef typename Model_t::Cost_t Cost_t;
    typedef typename Model_t::Cost_dx Cost_dx;
    typedef typename Model_t::Cost_dxx Cost_dxx;
    typedef typename Model_t::Cost_du Cost_du;
    typedef typename Model_t::Cost_duu Cost_duu;
    typedef typename Eigen::MatrixXd MatrixXd;
    typedef typename Eigen::VectorXd VectorXd;

    typedef std::list<VectorXd> StateList_t;
    typedef std::list<VectorXd> ControlList_t;
    typedef std::list<VectorXd> VxList ;
    typedef std::list<MatrixXd> VxxList ;
    typedef std::list<VectorXd> OpenLoopTermList ;
    typedef std::list<MatrixXd> GainList ;
    typedef std::list<double> TorqueList ;

  public: //protected: /* Parameters */
    Model_t model;
    int nbPreviewSteps;
    double timeStep;
    bool isInit;
    int compteur;
    int NTotal;
    int n; //size of the vector during optimization
    int window;

  public: //protected: /* Intermediate results */
    StateList_t stateList;
    StateList_t stateHist;
    ControlList_t controlList;
    VxList vxList;
    VxxList vxxList;
    GainList gainList;
    OpenLoopTermList openLoopList;
    Cost_t costI;
    double tCost;
    TorqueList torqueList;
    TorqueList torqueWanted;

    public:
        ilqr();
        //virtual ~ilqr();

    public:
    void displayV(std::list<VectorXd> L);
    void displayV(std::list<double> L);
    void init ();
    void initState(const State_t& state);
    void initTorque(const int i);
    void computeControl();
    void switching();
    void extend();

    public:
    void backwardLoop ();
    void forwardLoop(const double alpha);
    void Loops();
   /* void initLoops(); */
   /* void initBackwardLoop(); */
    void completeAlgo(const State_t& state);
    //void computeControl (const State_t& state);

    protected:
    private:
};


#endif // ILQR_H
