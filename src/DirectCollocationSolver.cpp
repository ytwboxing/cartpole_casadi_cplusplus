#include "DirectCollocationSolver.hpp"
using S1 = casadi::Slice;

DirectCollocationSolver::DirectCollocationSolver(const cart_pole_model& _model, const constraint_value& _constraint)
    :DirectTransSolver(_model, _constraint){
        solution = nullptr;
        minCartHorizonPos = opti.parameter();
        maxCartHorizonPos = opti.parameter();
        minU = opti.parameter();
        maxU = opti.parameter();
        T = opti.parameter();
        initialStateParameters = opti.parameter(4);
        finalStateParameters = opti.parameter(4);
        opti.set_value(minCartHorizonPos, _constraint.d_min);
        opti.set_value(maxCartHorizonPos, _constraint.d_max);
        opti.set_value(minU, _constraint.u_min);
        opti.set_value(maxU, _constraint.u_max);
    }
double DirectCollocationSolver::CubicHermitePoly(const double& x0, const double& xf, const double& v0, const double& vf, 
                    const double& totalTime, const double& nowTime)
    {
        double a, b, c, d, xt;
        double t = nowTime;
        a = x0;
        b = v0;
        c = -( (3 * (x0 - xf) + totalTime * (2 * v0 + vf)) ) / std::pow(totalTime, 2);
        d = ( 2 * (x0 - xf) + totalTime * (v0 + vf) ) / std::pow(totalTime, 3);

        xt = a + b * t + c * std::pow(t, 2) + d * std::pow(t, 3);

        return xt;
    }
double DirectCollocationSolver::SimpsonCollocation(const double& x0, const double& x_mid, const double& xf, const double& totalTime)
    {
        double xt = totalTime * (x0 + 4 * x_mid + xf) / 6;

        return xt;
    }
Function DirectCollocationSolver::getSystemDynamics(){
    casadi::MX X = casadi::MX::sym("x", 4);
    casadi::MX U = casadi::MX::sym("u", 1);
    casadi::MX dt = casadi::MX::sym("dt");
    MX A(2, 1);
    casadi::MX currentPosition = X(casadi::Slice(0,2));
    casadi::MX currentVelocity = X(casadi::Slice(2,4));

    casadi::MX temp1 = model.l * model.m2 * sin(currentPosition(1)) * pow(currentVelocity(1), 2);
    casadi::MX temp2 = model.m2 * model.g * cos(currentPosition(1)) * sin(currentPosition(1));
    casadi::MX temp3 = model.m1 + model.m2 * (1 - pow(cos(currentPosition(1)), 2));
    A(0) = (temp1 + U + temp2) / temp3;

    temp1 = model.l * model.m2 * cos(currentPosition(1)) * sin(currentPosition(1)) * pow(currentVelocity(1), 2);
    temp2 = (model.m1 + model.m2) * model.g * sin(currentPosition(1));
    temp3 = model.l * model.m1 + model.l * model.m2 * (1 - pow(cos(currentPosition(1)), 2));
    A(1) = (temp1 + U * cos(currentPosition(1)) + temp2) / temp3;
    casadi::MX dynamics = casadi::MX::vertcat({X(S1(2,4)), A});

    return casadi::Function("dynamics", {X, U, dt}, {dynamics});
}
void DirectCollocationSolver::setOptColloc(){
    casadi_int phaseLength = static_cast<casadi_int> (settings.phaseLength);
    opti.set_value(T, settings.time);
    casadi_int N = 2 * phaseLength;
    X = opti.variable(4, N + 1);
    A = opti.variable(2, N + 1);
    U = opti.variable(1, N + 1);

    opti.subject_to(X(S1(), 0) == initialStateParameters);
    opti.subject_to(X(S1(), N) == finalStateParameters);

    casadi::MX dT = T / phaseLength;
    MX f_curr, f_mid, f_next;
    MX costFunction = 0;
    costWeights w = settings._costWeights;
    for(casadi_int k = 0; k < N; ++k){
        //samples: x(k), x(k+1), x(k+2),,,,,,
        if(k % 2 == 0 && k + 2 <= N){
            f_curr = casadi::MX::vertcat(systemDynamics({X(S1(), k), U(S1(), k), 0}));
            f_next = casadi::MX::vertcat(systemDynamics({X(S1(), k + 2), U(S1(), k + 2), 0}));
            f_mid = casadi::MX::vertcat(systemDynamics({X(S1(), k + 1), U(S1(), k + 1), 0}));
            opti.subject_to(X(S1(), k + 2) == X(S1(), k) + dT * (f_curr + 4 * f_mid + f_next) / 6);
            
            costFunction += w.control * dT / 6 * ( pow(U(0,k),2) + 4 * pow(U(0,k+1),2) + pow(U(0,k+2),2));
        }
        /*
        collocation points, because states uses cubic hermite polynomial, collocation points are midium points at 
        each interval
        x(k+0.5), x(k+1+0.5),,,,,
        */
        else if(k % 2 != 0){
            opti.subject_to(X(S1(), k) == 0.5 * (X(S1(), k - 1) + X(S1(), k + 1)) + dT * (f_curr - f_next) / 8);
            opti.subject_to(U(S1(), k) == 0.5 * (U(S1(), k - 1) + U(S1(), k + 1)));

        }
        opti.subject_to(minCartHorizonPos <= X(0, k) <= maxCartHorizonPos);
        opti.subject_to(minU <= U(0, k) <= maxU);
    }   
    opti.subject_to(minCartHorizonPos <= X(0, N) <= maxCartHorizonPos);
    opti.subject_to(minU <= U(0, N) <= maxU);
    opti.minimize(costFunction);
}
bool DirectCollocationSolver::setupProblemColloc(const Settings& _settings){
    settings = _settings;
    systemDynamics = getSystemDynamics();
    //accelerationConsisitencyConstraint = getAccelerationConsistencyConstraintFunction();

    setOptColloc();

    casadi::Dict casadiOptions;
    casadi::Dict ipoptOptions;

    casadiOptions["expand"] = true;//Replace MX with SX expressions in problem formulation, speed up
    unsigned long solverVerbosity = settings.solverVerbosity;
    if (solverVerbosity) {
        casadi_int ipoptVerbosity = static_cast<long long>(solverVerbosity - 1);
        ipoptOptions["print_level"] = ipoptVerbosity;
        casadiOptions["print_time"] = true;
        casadiOptions["bound_consistency"] = false;
    } else {
        ipoptOptions["print_level"] = 0;
        casadiOptions["print_time"] = false;
        //casadiOptions["bound_consistency"] = false;
        //ipoptOptions["fixed_variable_treatment"] = "make_constraint";
    }
    ipoptOptions["linear_solver"] = settings.ipoptLinearSolver;

    opti.solver("ipopt", casadiOptions, ipoptOptions);

    solverState = SolverState::PROBLEM_SET;

    return true;
}
bool DirectCollocationSolver::solveColloc(const State& initialState, const State& finalState){
    if(solverState == SolverState::NOT_INITIALIZED){
        throw std::runtime_error("problem not initialized");
        return false;
    }
    setParametersValue(initialState, finalState);
    casadi_int npoints = 2 * static_cast<casadi_int> (settings.phaseLength);
    casadi::DM initPos = casadi::DM::zeros(2,1);
    casadi::DM finalPos = casadi::DM::zeros(2,1);
    initPos = initialState.state(casadi::Slice(0,2));
    finalPos = finalState.state(casadi::Slice(0,2));
    casadi::DM interpolatedPosition(2, 1);
    casadi::DM linSpacePoints = casadi::DM::linspace(0, 1, npoints + 1);
    for(casadi_int k = 0; k < npoints + 1; ++k){
        /*
        initial guess, pos using linear interpolatation from init to final, 
        vel and control all set zero
        */
        interpolatedPosition = initPos + linSpacePoints(k) * (finalPos - initPos);
        opti.set_initial(X(casadi::Slice(0,2), k), interpolatedPosition);
        opti.set_initial(X(casadi::Slice(2,4), k), 0);
    }
    for(casadi_int k = 0; k < npoints; ++k){
        opti.set_initial(U(casadi::Slice(), k), 0);
    }
    
    solverState = SolverState::PROBLEM_SET;
    
    try{
        solution = std::make_unique<casadi::OptiSol>(opti.solve());
    }catch(std::exception &e){
        opti.debug().show_infeasibilities(1e-5);
        std::cerr << "error while solving the optimization" << std::endl;
        std::cerr << "Details:\n " << e.what() << std::endl;
        return false;
    }
    solverState = SolverState::PROBLEM_SOLVED;
    std::cout << "\nsolve success\n\n";
    return true;
}
void DirectCollocationSolver::getSolutionColloc(casadi::DM& state, casadi::DM& control){
    DM state_all = solution -> value(X);
    DM control_all = solution -> value(U);
    int a = 0;
    state = DM::zeros(4, settings.phaseLength + 1);
    control = DM::zeros(1, settings.phaseLength + 1);
    for(int i = 0; i < 2 * settings.phaseLength + 1; ++i){
        if(i % 2 == 0){
            //std::cout <<"i\n" << i << "\n"<<control_all(S1(), i) << "\n\n";
            state(S1(),a) = DM::vertcat({state_all(S1(), i)});
            control(S1(),a) = control_all(S1(), i);
            a++;
        }
    }
}
