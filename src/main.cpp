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
    (*q)[0] = A * exp(-pow((x - mu) / sigma, 2));
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
    //grid_transport.setupY(f_transport_discont, {1.,0.5});
    TransportEquation transport{0.5};
    grid_transport.updateBoundary(transport);


    for (int i=0; i<grid_transport.totalsize(); i++)
    {
        Transport transport = grid_transport.getY<Transport>(i);
        std::cout << grid_transport.getX(i) << " " << transport[0] << std::endl;
    }

    Upwind<Transport,1,1> upwind{transport};
    CFL<Transport,1,1>CFL_transport{0.25, transport};
    explicitStep<Transport,1,1> step_transport{transport, upwind};

    std::cout << "CFL is called" << std::endl;
    std::cout << CFL_transport(grid_transport, 0.01) << std::endl;
    std::cout << "Start solver" << std::endl;
    Solver<Transport, explicitStep<Transport,1,1>,CFL<Transport,1,1>,1> solver_transport{grid_transport, grid_transport.getDx(), CFL_transport, step_transport};


    solver_transport.setT0(0.0);
    solver_transport.setT1(100*CFL_transport(grid_transport, grid_transport.getDx()));
    solver_transport.solve();

    std::ofstream output_transport{"out_transport.txt"};

    auto out_grid=solver_transport.getGrid();
    for (int i = 0; i < out_grid.totalsize(); i++)
    {
        Transport transport = out_grid.getY<Transport>(i);
        Transport transport_old = grid_transport.getY<Transport>(i);
        output_transport << out_grid.getX(i) << " " << transport[0]  << " " << transport_old[0] << std::endl;
    }


    return 0;
}