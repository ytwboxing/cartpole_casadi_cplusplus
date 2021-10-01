#include <casadi/casadi.hpp>
#include <iostream>

int main(){
    casadi::Opti opti;
    casadi::MX x,y;
    std::unique_ptr<casadi::OptiSol> sol;
    casadi::DM D_sol;

    x = opti.variable();
    y = opti.variable();

    opti.minimize(pow(x, 3)+pow(y, 2));
    opti.subject_to(0 <= x <= 0.4);
    opti.subject_to(0.7 <= y <= 0.9);

    casadi::Dict casadiOptions;
    casadi::Dict ipoptOptions;
    casadiOptions["expand"] = true;//Replace MX with SX expressions in problem formulation, speed up
    ipoptOptions["print_level"] = 0;
    ipoptOptions["linear_solver"] = "mumps";
    opti.solver("ipopt", casadiOptions, ipoptOptions);
    opti.set_initial(x, 0);
    opti.set_initial(y, 0);
    sol = std::make_unique<casadi::OptiSol>(opti.solve());

    D_sol = casadi::DM::vertcat({sol -> value(x), sol -> value(y)});
    std::cout << "x: " << D_sol(0) << "\n" << "y: " << D_sol(1) << std::endl;

    return 0;
}


