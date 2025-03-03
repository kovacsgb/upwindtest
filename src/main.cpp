#include <iostream>
#include "Grid.hpp"
#include "Equation.hpp"
#include "Solver.hpp"
#include "Flux.hpp"
#include "Solverstep.hpp"
#include "Quantity.hpp"
#include <cmath>
#include <fstream>


Quantity<3>* f(double x, std::vector<double> params)
{
    Euler* q = new Euler();
    double A = params[0];
    double mu = params[1];
    double sigma = params[2];
    (*q)[0] = A * exp(-pow((x - mu) / sigma, 2))+A*0.1;
    (*q)[1] = (*q)[0];
    (*q)[2] = 3/2 * std::pow((*q)[0], 1.4)+0.5*(*q)[1]*(*q)[1]/(*q)[0];

    return q;
}

Transport* f_transport_gauss(double x, std::vector<double> params)
{
    Transport* q = new Transport();
    double A = params[0];
    double mu = params[1];
    double sigma = params[2];
    (*q)[0] = A * exp(-pow((x - mu) / sigma, 2));

    return q;
}
Transport* f_transport_discont(double x, std::vector<double> params)
{
    Transport* q = new Transport();
    double A = params[0];
    double mu = params[1];
    (*q)[0] = A ? x< mu :0;

    return q;
}

Euler* f_Sod_shock_tube(double x, std::vector<double> params)
{
    Euler *q = new Euler();
    double rhoL = params[0];
    double uL = params[1];
    double pL = params[2];
    double rhoR = params[3];
    double uR = params[4];
    double pR = params[5];
    double x0 = params[6];

    if(x <= x0)
    {
        (*q)[0] = rhoL;
        (*q)[1] = rhoL * uL;
        (*q)[2] = 3/2*pL + 0.5*rhoL*uL*uL;
    }
    else
    {
        (*q)[0] = rhoR;
        (*q)[1] = rhoR * uR;
        (*q)[2] = 3/2*pR + 0.5*rhoR*uR*uR;
    }

    return q;
}


int main() {

    std::cout << "Test Quantity" << std::endl;
    std::cout << "----------------" << std::endl;
    Euler q1 = {1,2,3};
    std::unique_ptr<Quantity<3>> q2 = q1.clone();
    std::unique_ptr<Quantity<3>> q3 = q1.clone();
    std::cout << q1[0] << " " << (*q2)[0] << " " << (*q3)[0] << std::endl;
    std::cout << q1[1] << " " << (*q2)[1] << " " << (*q3)[1] << std::endl;
    std::cout << q1[2] << " " << (*q2)[2] << " " << (*q3)[2] << std::endl;

    std::cout << "Test Grid" << std::endl;
    std::cout << "----------------" << std::endl;
    Grid<3> grid1 = Grid<3>::GenerateFromBorders(100, 0.0, 1.0);
    grid1.setupY(f, {1.0, 0.5, 0.1});
    Grid<3> grid2 = grid1;
    Grid<3> grid3 = grid1;
    std::cout << grid1.getX(0) << " " << grid2.getX(0) << " " << grid3.getX(0) << std::endl;

    Grid<3> grid = Grid<3>::GenerateFromBorders(100, 0.0, 1.0);
    std::vector<double> params = {1.0, 0.5, 0.1};
    grid.setupY(f, params);

    std::cout << (Euler{1,1,1} / Euler{2,2,2})[0];

    std::ofstream output{"out.txt"};

    for (int i = 0; i < grid.totalsize(); i++)
    {
        //std::cout << "Hello, World!" << std::endl;
        Euler euler = grid.getY<Euler>(i);
        //std::cout << "Hello, World!" << std::endl;
        output << grid.getX(i) << " " << euler.getrho() << " " << euler.getu() << " " << euler.getp() << std::endl;
    }
    
    EulerEquation equation;
    CFL<Euler,3,1> testCFL{0.5, equation};
    std::cout << testCFL(grid, 0.01) << std::endl;
    std::cout << "Hello, World!" << std::endl;

    std::cout << "First test: Upwind with transport" << std::endl;
    std::cout << "---------------------------------" << std::endl;
    
    Grid<1> grid_transport = Grid<1>::GenerateFromBorders(500, 0.0, 5.0);
    std::vector<double> params_transport = {1.0, 0.5, 0.1};
    grid_transport.setupY(f_transport_gauss, params_transport);
    grid_transport.setupY(f_transport_discont, {1.,0.5});
    TransportEquation transport{0.5};
    grid_transport.updateBoundary(transport);


    for (int i=0; i<grid_transport.totalsize(); i++)
    {
        Transport transport = grid_transport.getY<Transport>(i);
        std::cout << grid_transport.getX(i) << " " << transport[0] << std::endl;
    }

    //Upwind<Transport,1,1> upwind{transport};
    vonLeer<Transport,1,1> upwind{transport};
    CFL<Transport,1,1>CFL_transport{0.25, transport};
    explicitStep<Transport,1,1> step_transport{transport, upwind};

    std::cout << "CFL is called" << std::endl;
    std::cout << CFL_transport(grid_transport, 0.01) << std::endl;
    std::cout << "Start solver" << std::endl;
    Solver<Transport, explicitStep<Transport,1,1>,CFL<Transport,1,1>,1> solver_transport{grid_transport, grid_transport.getDx(), CFL_transport, step_transport};


    solver_transport.setT0(0.0);
    solver_transport.setT1(200*CFL_transport(grid_transport, grid_transport.getDx()));
    solver_transport.solve();

    Upwind<Transport,1,1> upwind_transport{transport};
    explicitStep<Transport,1,1> step_transport2{transport, upwind_transport};
    Solver<Transport, explicitStep<Transport,1,1>,CFL<Transport,1,1>,1> solver_transport2{grid_transport, grid_transport.getDx(), CFL_transport, step_transport2};
    solver_transport2.setT0(0.0);
    solver_transport2.setT1(200*CFL_transport(grid_transport, grid_transport.getDx()));
    solver_transport2.solve();


    std::ofstream output_transport{"out_transport.txt"};

    auto out_grid=solver_transport.getGrid();
    auto out_grid2=solver_transport2.getGrid();
    for (int i = 0; i < out_grid.totalsize(); i++)
    {
        Transport transport = out_grid.getY<Transport>(i);
        Transport transport_old = grid_transport.getY<Transport>(i);
        Transport transport2 = out_grid2.getY<Transport>(i);
        output_transport << out_grid.getX(i) << " " << transport[0]  << " " << transport_old[0] << " " << transport2[0] << std::endl;
    }

    std::cout << "Second test: Burgers equation" << std::endl;
    std::cout << "---------------------------------" << std::endl;
    Grid<1> grid_burgers = Grid<1>::GenerateFromBorders(500, 0.0, 5.0);
    grid_burgers.setupY(f_transport_discont, {1,0.5});
    BurgersEquation burgers;
    grid_burgers.updateBoundary(burgers);

    for (int i = 0; i < grid_burgers.totalsize(); i++)
    {
        Transport transport = grid_burgers.getY<Transport>(i);
        std::cout << grid_burgers.getX(i) << " " << transport[0] << std::endl;
    }

    Upwind<Transport,1,1> upwind_burgers{burgers};
    vonLeer<Transport,1,1> vonLeer_burgers{burgers};
    CFL<Transport,1,1> CFL_burgers{0.25, burgers};
    explicitStep<Transport,1,1> step_burgers{burgers, upwind_burgers};
    explicitStep<Transport,1,1> step_burgers_vonLeer{burgers, vonLeer_burgers};

    std::cout << "CFL is called" << std::endl;
    std::cout << CFL_burgers(grid_burgers, 0.01) << std::endl;
    std::cout << "Start solver" << std::endl;
    Solver<Transport, explicitStep<Transport,1,1>,CFL<Transport,1,1>,1> solver_burgers{grid_burgers, grid_burgers.getDx(), CFL_burgers, step_burgers};
    Solver<Transport, explicitStep<Transport,1,1>,CFL<Transport,1,1>,1> solver_burgers_vonLeer{grid_burgers, grid_burgers.getDx(), CFL_burgers, step_burgers_vonLeer};
    solver_burgers.setT0(0.0);
    solver_burgers.setT1(500 * CFL_burgers(grid_burgers, grid_burgers.getDx()));
    solver_burgers.solve();
    solver_burgers_vonLeer.setT0(0.0);
    solver_burgers_vonLeer.setT1(500 * CFL_burgers(grid_burgers, grid_burgers.getDx()));
    solver_burgers_vonLeer.solve();


    std::ofstream output_burgers{"out_burgers.txt"};

    auto out_grid_burgers = solver_burgers.getGrid();
    auto out_grid_burgers_vonLeer = solver_burgers_vonLeer.getGrid();
    for (int i = 0; i < out_grid_burgers.totalsize(); i++)
    {
        Transport transport = out_grid_burgers.getY<Transport>(i);
        Transport transport_old = grid_burgers.getY<Transport>(i); 
        Transport transport_vonLeer = out_grid_burgers_vonLeer.getY<Transport>(i);
        output_burgers << out_grid_burgers.getX(i) << " " << transport[0] << " " << transport_old[0]  << " " << transport_vonLeer[0]<< std::endl;
    }








    std::cout<< "Third test: Sod shock tube" << std::endl;
    std::cout << "---------------------------------" << std::endl;

    Grid<3> grid_sod = Grid<3>::GenerateFromBorders(200, -5, 5);
    std::vector<double> params_sod = {1.0, 0.0, 1.0, 0.125, 0.0, 0.1, 0.0};
    grid_sod.setupY(f_Sod_shock_tube, params_sod);
    //grid_sod.setupY(f, {1.0, 0.0,0.2});
    EulerEquation equation_sod;
    grid_sod.updateBoundary(equation_sod);


    std::cout << "Dump starting position" << std::endl;
    std::ofstream sod_start_out{"sod_start.txt"};
    for(int i=0; i<grid_sod.totalsize(); i++)
    {
        Euler euler = grid_sod.getY<Euler>(i);
        sod_start_out << grid_sod.getX(i) << " " << euler.getrho() << " " << euler.getu() << " " << euler.getp() << std::endl;
    }

    std::cout << "Set up solver" << std::endl;
    vonLeer<Euler,1,3> upwind_sod(equation_sod);
    FluxWithViscosity flux_with_visc(upwind_sod,4);
    CFL<Euler,3,1> CFL_sod{0.1, equation_sod};
    //explicitStep<Euler,1,3> step_sod{equation_sod, upwind_sod};
    explicitStep<Euler,1,3> step_sod{equation_sod, flux_with_visc};
    //Runge2Step<Euler,1,3> step_sod{equation_sod, flux_with_visc,grid_sod};
    Solver<Euler, explicitStep<Euler,1,3>,CFL<Euler,3,1>,3> solver_sod{grid_sod, grid_sod.getDx(), CFL_sod, step_sod};
    //Solver<Euler, Runge2Step<Euler,1,3>,CFL<Euler,3,1>,3> solver_sod{grid_sod, grid_sod.getDx(), CFL_sod, step_sod};

    solver_sod.setT0(0.0);
    solver_sod.setT1(1.7);
    solver_sod.solve();

    std::ofstream output_sod{"out_sod.txt"};

    auto out_grid_sod=solver_sod.getGrid();
    for (int i = 0; i < out_grid_sod.totalsize(); i++)
    {
        Euler euler = out_grid_sod.getY<Euler>(i);
        output_sod << out_grid_sod.getX(i) << " " << euler.getrho() << " " << euler.getu() << " " << euler.getp() << std::endl;
    }



    return 0;
}