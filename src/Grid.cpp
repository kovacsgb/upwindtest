#include "Grid.hpp"

template <int M>
Grid<M>::Grid(int N, double x0, double dx) : N(N), dx(dx)
{
    x.resize(N);
    for (int i = 0; i < N; i++)
    {
        x[i] = x0 + i * dx;
    }
}

template <int M>
Grid<M> Grid<M>::GenerateFromBorders(int N, double x0, double x1)
{
    double dx = (x1 - x0) / (N - 1);
    return Grid<M>(N, x0, dx);
}


template <int M>
double &Grid<M>::getX(int i)
{
    return x[i];
}

template <int M>
template <typename T>
T& Grid<M>::getY(int i)
{
    return *dynamic_cast<T*>(y[i].get());
}

template <int M>
int Grid<M>::size()
{
    return N;
}

template <int M>
double Grid<M>::getDx()
{
    return dx;
}  

template <int M>
void Grid<M>::setupY(std::function<Quantity<M>(double, std::vector<double>)> f, std::vector<double> params)
{
    y.resize(N);
    std::transform(x.begin(), x.end(), y.begin(), [&](double xi) {
        return std::make_unique<Quantity<M>>(f(xi, params));
    });
}