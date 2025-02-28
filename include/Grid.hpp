#ifndef GRID_HPP
#define GRID_HPP


#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <memory>
#include <functional>

#include "Quantity.hpp"
#include "Equation.hpp"

template<typename T,int N> class EquationBase;

template<int M>
struct Grid
{
    private:
        int N;
        int NGhost;
        double dx;
        std::vector<double> x;
        std::vector<std::unique_ptr<Quantity<M>>> y;
    public:
        Grid() = default;
        Grid(int N, int NGhost, double x0, double dx);

        Grid(const Grid& grid_) : N(grid_.N), NGhost(grid_.NGhost), dx(grid_.dx), x(grid_.x)
        {
            //std::cerr << "Grid copy contructor called" << std::endl;
            //std::cerr << "Data is:" << N << " " << NGhost << " " << dx << std::endl;
            y.resize(N+2*NGhost);
            std::transform(grid_.y.begin(), grid_.y.end(), y.begin(), [&](const auto& i) {
                return std::move(i->clone());
            });

            size_t j=0;
           /* for (const auto& i : grid_.y)
            {
                //std::cerr << j++ << std::endl;
                std::unique_ptr<Quantity<M>> q = i->clone();
                //std::cerr << "Quantity cloned" << std::endl;
                y.push_back(std::move(q));
            }*/
        }

        Grid& operator=(const Grid& grid_)
        {
            //std::cerr << "Grid copy = called" << std::endl;
            N = grid_.N;
            NGhost = grid_.NGhost;
            dx = grid_.dx;
            x = grid_.x;
            y.resize(N+2*NGhost);
            size_t j=0;
            std::transform(grid_.y.begin(), grid_.y.end(), y.begin(), [&](const auto& i) {
                return std::move(i->clone());
            });
            return *this;
        }

        static Grid GenerateFromBorders(int N, double x0, double x1);
        double& getX(int i);
        template<typename T>
        T& getY(int i);
        template<typename T>
        T& getY(int i) const;

    
        int size();
        int totalsize() { return N+2*NGhost; }
        int startid();
        double getDx();
        auto begin() { return x.begin(); }
        auto end() { return x.end(); }
        void setupY(std::function<Quantity<M>*(double,std::vector<double>)> f, std::vector<double> params);

        template<typename T, int L>
        void updateBoundary(EquationBase<T,L>& eq)
        {
            eq.template updateBoundary<M>(*this);
        }
};

template <int M>
Grid<M>::Grid(int N, int NGhost, double x0, double dx) : N(N), NGhost(NGhost), dx(dx)
{
    x.resize(N+2*NGhost);
    for (int i = 0; i < N+2*NGhost; i++)
    {
        x[i] = x0 + i * dx;
    }
}

template <int M>
Grid<M> Grid<M>::GenerateFromBorders(int N, double x0, double x1)
{
    double dx = (x1 - x0) / (N - 1);
    return Grid<M>(N,2, x0-2*dx, dx);
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
template<typename T>
T& Grid<M>::getY(int i) const
{
    return *dynamic_cast<T*>(y[i].get());
}
template <int M>
int Grid<M>::size()
{
    return N+NGhost;
}


template <int M>
int Grid<M>::startid()
{
    return NGhost;
}

template <int M>
double Grid<M>::getDx()
{
    return dx;
}  

template <int M>
void Grid<M>::setupY(std::function<Quantity<M>*(double, std::vector<double>)> f, std::vector<double> params)
{
    y.resize(N+2*NGhost);
    std::transform(x.begin(), x.end(), y.begin(), [&](double xi) {
        return std::unique_ptr<Quantity<M>>(f(xi, params));
    });
}
#endif //GRID_HPP
