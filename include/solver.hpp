#ifndef SOLVER_H
#define SOLVER_H

#include <common.hpp>
#include <vector>
#include <string>
#include <cassert>
#include <memory>

class Solver{
public:
    Solver(const cart_pole_model& _model, const constraint_value& _constraint);
    ~Solver() { };
    bool setupProblem(const Settings& _settings);
    bool solve(const State& initialState, const State& finalState);
    void getSolution(casadi::DM& state, casadi::DM& control);

private:
    enum SolverState{NOT_INITIALIZED=0, PROBLEM_SET, PROBLEM_SOLVED};
    SolverState solverState;
    Settings settings;
    casadi::Function integratorDynamics;
    casadi::Function accelerationConsisitencyConstraint;
    casadi::MX initialStateParameters;
    casadi::MX finalStateParameters;
    casadi::MX X, A, U, T;
    casadi::MX minCartHorizonPos, maxCartHorizonPos, minU, maxU;
    
    cart_pole_model model;

    casadi::Opti opti;
    std::unique_ptr<casadi::OptiSol> solution;
    casadi::Function getIntegratorDynamics();
    casadi::Function getAccelerationConsistencyConstraintFunction();
    void setupOpti();
    void setParametersValue(const State& initialState, const State& finalState);
};

#endif