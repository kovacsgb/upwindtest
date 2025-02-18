#ifndef GRID_HPP
#define GRID_HPP


#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>


struct Grid
{
    private:
        int N;
        double dx;
        std::vector<double> x;
    public:
        Grid() = default;
        Grid(int N,double x0, double dx);
        static Grid GenerateFromBorders(int N, double x0, double x1);


};

#endif //GRID_HPP
