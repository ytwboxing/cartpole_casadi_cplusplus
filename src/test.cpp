#include "solver.hpp"

int main(){
    Settings settings;
    settings.phaseLength = 10;
    settings._costWeights.control = 1;
    settings.solverVerbosity = 1;
    settings.ipoptLinearSolver = "mumps";
    
    State initialState, finalState;
    initialState.state.zeros(4);
    finalState.state = {1, 3.14, 0, 0};
    
    cart_pole_model model = {0.5, 1, 0.3, 9.81};
    constraint_value constraint = {0, 2, -20, 20};

    Solver solver(model, constraint);
    solver.setupProblem(settings);
    bool ok = solver.solve(initialState, finalState);
    
    casadi::DM sol_state, sol_control;
    if(ok) solver.getSolution(sol_state, sol_control);

    std::cout << sol_state << "\n\n" << sol_control << std::endl;
    return 0;
}