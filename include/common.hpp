#ifndef COMMON_H
#define COMMON_H

#include <casadi/casadi.hpp>

typedef struct{
    double control;
}costWeights;

typedef struct{
    int phaseLength;
    double time;
    costWeights _costWeights;
    double solverVerbosity;
    /*
    ipopt's available linear solvers = {"ma27", "ma57", "ma77", "ma86", "ma97", "pardiso", "wsmp", "mumps"};
    used for it's interior point method
    */
    std::string ipoptLinearSolver;
}Settings;

typedef struct {
    casadi::DM state = casadi::DM::zeros(4,1);
}State;

typedef struct{
    double l;
    double m1;
    double m2;
    double g;
}cart_pole_model;

typedef struct{
    double d_min;
    double d_max;
    double u_min;
    double u_max;
}constraint_value;
#endif