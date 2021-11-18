#include "DirectTranscriptSolver.hpp"

DirectTransSolver::DirectTransSolver(const cart_pole_model& _model, const constraint_value& _constraint)
 {
    model = _model;
    solverState = SolverState::NOT_INITIALIZED;
    transMethod = TranscriptMethod::SINGLE_SHOOTING;
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

 casadi::Function DirectTransSolver::getIntegratorDynamics()
 {
     casadi::MX x = casadi::MX::sym("x", 4);
     casadi::MX a = casadi::MX::sym("a", 2);
     casadi::MX dT = casadi::MX::sym("dt");

     casadi::MX p = x(casadi::Slice(0,2));
     casadi::MX v = x(casadi::Slice(2,4));
     casadi::MX rhs = casadi::MX::vertcat({p + dT * (v + 0.5 * dT * a),
                                           v + dT * a});
    
    return casadi::Function("Integrator", {x,a,dT}, {rhs});
 }

 casadi::Function DirectTransSolver::getAccelerationConsistencyConstraintFunction()
 {
    casadi::MX X = casadi::MX::sym("x", 4);
    casadi::MX U = casadi::MX::sym("u", 1);
    casadi::MX A = casadi::MX::sym("a", 2);

    casadi::MX currentPosition = X(casadi::Slice(0,2));
    casadi::MX currentVelocity = X(casadi::Slice(2,4));

    casadi::MX temp1 = model.l * model.m2 * sin(currentPosition(1)) * pow(currentVelocity(1), 2);
    casadi::MX temp2 = model.m2 * model.g * cos(currentPosition(1)) * sin(currentPosition(1));
    casadi::MX temp3 = model.m1 + model.m2 * (1 - pow(cos(currentPosition(1)), 2));
    casadi::MX constraint_cart = A(0) - (temp1 + U + temp2) / temp3;

    temp1 = model.l * model.m2 * cos(currentPosition(1)) * sin(currentPosition(1)) * pow(currentVelocity(1), 2);
    temp2 = (model.m1 + model.m2) * model.g * sin(currentPosition(1));
    temp3 = model.l * model.m1 + model.l * model.m2 * (1 - pow(cos(currentPosition(1)), 2));
    casadi::MX constraint_pole = A(1) + (temp1 + U * cos(currentPosition(1)) + temp2) / temp3;

    casadi::MX constarint = casadi::MX::vertcat({constraint_cart, constraint_pole});
    return casadi::Function("accelerationConsistency", {X, U, A}, {constarint});

 }

 void DirectTransSolver::setupOpti()
 {
    using S1 = casadi::Slice;
    casadi_int phaseLength, N;
    opti.set_value(T, settings.time);
    if(transMethod == TranscriptMethod::SINGLE_SHOOTING){
        phaseLength = static_cast<casadi_int> (settings.phaseLength);
        N = phaseLength;
    }
    else if(transMethod == TranscriptMethod::MULTIPLE_SHOOTING){
        phaseLength = static_cast<casadi_int> (settings.phaseLength / 2);
        N = phaseLength * 2;
    }
    X = opti.variable(4, N + 1);
    A = opti.variable(2, N);
    // casadi::MX X(4, N + 1);// = opti.variable(4, N + 1);
    // casadi::MX A(2, N);// = opti.variable(2, N);
    U = opti.variable(1, N);

    opti.subject_to(X(S1(), 0) == initialStateParameters);
    opti.subject_to(X(S1(), N) == finalStateParameters);
    if(transMethod == TranscriptMethod::SINGLE_SHOOTING){
        casadi::MX dT = T / phaseLength;
        for(casadi_int k = 0; k < N; ++k){
            opti.subject_to(X(S1(), k + 1) == casadi::MX::vertcat(integratorDynamics({X(S1(), k), A(S1(), k), dT})));
            opti.subject_to(minCartHorizonPos <= X(0, k) <= maxCartHorizonPos);
            opti.subject_to(minU <= U(0, k) <= maxU);
            opti.subject_to(casadi::MX::vertcat(accelerationConsisitencyConstraint({X(S1(), k), U(S1(), k), A(S1(), k)})) == casadi::MX::zeros(2, 1));
        } 
    }
    else if(transMethod == TranscriptMethod::MULTIPLE_SHOOTING){
        casadi::MX dT = T / (phaseLength * 2);
        for(casadi_int phase = 0;phase<2;++phase){
            for(casadi_int k = phaseLength*phase;k < phaseLength*(phase+1);++k){
                if(k != (phaseLength-1)){
                    opti.subject_to(X(S1(), k + 1) == casadi::MX::vertcat(integratorDynamics({X(S1(), k), A(S1(), k), dT})));
                }
                opti.subject_to(minCartHorizonPos <= X(0, k) <= maxCartHorizonPos);
                opti.subject_to(minU <= U(0, k) <= maxU);

                opti.subject_to(casadi::MX::vertcat(accelerationConsisitencyConstraint({X(S1(), k), U(S1(), k), A(S1(), k)})) == casadi::MX::zeros(2, 1));
            }
    }
    opti.subject_to( X(S1(), phaseLength) 
        == casadi::MX::vertcat(integratorDynamics({X(S1(), phaseLength-1), A(S1(), phaseLength-1), dT}))
        - casadi::DM::vertcat({0.2, 0.4, 0.01, 0.01}) 
        );
    }
    costWeights w = settings._costWeights;
    casadi::MX costFunction = w.control * casadi::MX::sumsqr(U(0, S1()));
    opti.minimize(costFunction);
 }

bool DirectTransSolver::setupProblem(const Settings& _settings)
{
    settings = _settings;
    integratorDynamics = getIntegratorDynamics();
    accelerationConsisitencyConstraint = getAccelerationConsistencyConstraintFunction();

    setupOpti();

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

void DirectTransSolver::setParametersValue(const State& initialState, const State& finalState)
{
    opti.set_value(initialStateParameters, initialState.state);
    opti.set_value(finalStateParameters, finalState.state);
}

bool DirectTransSolver::solve(const State& initialState, const State& finalState)
{
    if(solverState == SolverState::NOT_INITIALIZED){
        throw std::runtime_error("problem not initialized");
        return false;
    }
    setParametersValue(initialState, finalState);
    casadi_int npoints = static_cast<casadi_int> (settings.phaseLength);

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

void DirectTransSolver::getSolution(casadi::DM& state, casadi::DM& control){
    casadi::DM state_all = solution -> value(X);
    casadi::DM control_all = solution -> value(U);
    state = casadi::DM::zeros(4, settings.phaseLength + 1);
    control = casadi::DM::zeros(1, settings.phaseLength + 1);
    state = state_all;
    control = control_all;
}