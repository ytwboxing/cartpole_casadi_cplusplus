#pragma once

#include "DirectTranscriptSolver.hpp"

using casadi::DM;
using casadi::MX;
using casadi::Function;
using casadi::Opti;

class DirectCollocationSolver: public DirectTransSolver{
public:
    DirectCollocationSolver(const cart_pole_model& _model, const constraint_value& _constraint);
    ~DirectCollocationSolver(){ };
    bool setupProblemColloc(const Settings& _settings);
    bool solveColloc(const State& initialState, const State& finalState);
    void getSolutionColloc(casadi::DM& state, casadi::DM& control);
private:
    double CubicHermitePoly(const double& x0, const double& xf, const double& v0, const double& vf, 
                    const double& totalTime, const double& nowTime);
    /*
    Simpson quadrature, also known as Simpson’s rule for integration, works by approximating the 
    integrand of the integral as a piecewise quadratic function
    see https://epubs.siam.org/doi/pdf/10.1137/16M1062569
    */
    double SimpsonCollocation(const double& x0, const double& x_mid, const double& xf, const double& totalTime);
    Function systemDynamics;
    Function getSystemDynamics();
    void setOptColloc();
};