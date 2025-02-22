#include <iostream>
#include "Grid.hpp"
#include "Quantity.hpp"
#include <cmath>


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

int main() {

    Grid<3> grid = Grid<3>::GenerateFromBorders(100, 0.0, 1.0);
    std::vector<double> params = {1.0, 0.5, 0.1};
    grid.setupY(f, params);

    for (int i = 0; i < grid.size(); i++)
    {
        //std::cout << "Hello, World!" << std::endl;
        Euler euler = grid.getY<Euler>(i);
        //std::cout << "Hello, World!" << std::endl;
        std::cout << grid.getX(i) << " " << euler.getrho() << " " << euler.getu() << " " << euler.getp() << std::endl;
    }
    std::cout << "Hello, World!" << std::endl;
    return 0;
}