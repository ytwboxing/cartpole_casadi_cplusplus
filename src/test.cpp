#include "DirectTranscriptSolver.hpp"
#include "DirectCollocationSolver.hpp"

int main(){
    Settings settings;
    settings.phaseLength = 10;
    settings.time = 2;
    settings._costWeights.control = 1;
    settings.solverVerbosity = 1;
    settings.ipoptLinearSolver = "mumps";
    
    State initialState, finalState;
    initialState.state.zeros(4);
    finalState.state = {1, 3.14, 0, 0};
    
    cart_pole_model model = {0.5, 1, 0.3, 9.81};
    constraint_value constraint = {0, 2, -100, 100};

    DirectTransSolver solver(model, constraint);
    solver.setTranscriptMethod(TranscriptMethod::SINGLE_SHOOTING);
    solver.setupProblem(settings);
    bool ok = solver.solve(initialState, finalState);
    
    casadi::DM sol_state, sol_control;
    if(ok) solver.getSolution(sol_state, sol_control);
    std::cout<<"single shooting transcript\n";
    std::cout << "state:\n" << sol_state << "\ncontrol:\n" << sol_control << "\n\n";

    solver.setTranscriptMethod(TranscriptMethod::MULTIPLE_SHOOTING);
    solver.setupProblem(settings);
    ok = solver.solve(initialState, finalState);
    
    if(ok) solver.getSolution(sol_state, sol_control);
    std::cout<<"multiple shooting transcript\n";
    std::cout << "state:\n" << sol_state << "\ncontrol:\n" << sol_control << "\n\n";

    Settings co_settings;
    co_settings.phaseLength = 10;
    co_settings.time = 2;
    co_settings._costWeights.control = 1;
    co_settings.solverVerbosity = 1;
    co_settings.ipoptLinearSolver = "mumps";
    DirectCollocationSolver co_solver(model, constraint);
    co_solver.setupProblemColloc(co_settings);
    bool co_ok = co_solver.solveColloc(initialState, finalState);
    casadi::DM co_sol_state, co_sol_control;
    if(co_ok) co_solver.getSolutionColloc(co_sol_state, co_sol_control);
    std::cout<<"collocation\n";
    std::cout << "state:\n" << co_sol_state << "\ncontrol:\n" << co_sol_control << std::endl;

    return 0;
}