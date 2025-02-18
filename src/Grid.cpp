#include "Grid.hpp"

Grid::Grid(int N, double x0, double dx) : N(N), dx(dx)
{
    x.resize(N);
    for (int i = 0; i < N; i++)
    {
        x[i] = x0 + i * dx;
    }
}

Grid Grid::GenerateFromBorders(int N, double x0, double x1)
{
    double dx = (x1 - x0) / (N - 1);
    return Grid(N, x0, dx);
}
