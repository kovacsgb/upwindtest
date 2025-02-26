#include "Grid.hpp"
#include "Quantity.hpp"

template <typename T,int N>
class Solverstep
{
    public:

    virtual T operator()(Grid<N>& grid_old, double dx, double d, int i) = 0;
};


template <typename T,int N, typename F>
class explicitStep : public Solverstep<T,N>
{
    public:

    explicitStep(double cfl) : cfl(cfl) {}

    T operator()(Grid<N>& grid_old, double dx, double dt, int i) override
    {
        T old = grid_old.getY<T>(i);
        T new_val = old + dt
    }

    private:

};