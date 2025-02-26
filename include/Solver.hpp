#include "Grid.hpp"
//#include "Solverstep.hpp"
#include "Equation.hpp"
#include <algorithm>

template <typename T, int N, int M>
class CFL
{
    private:
    double cfl;
    EquationBase<T,M>& equation;

public:
    CFL(double cfl, EquationBase<T,M>& equation) : cfl(cfl), equation(equation) {}

    double operator()(Grid<N>& grid, double dx)
    {
        auto max_c_s = std::max_element(grid.begin(), grid.end(), [&](const auto& a, const auto& b) {
            return equation.getSoundSpeed(grid.template getY<T>(a)) < equation.getSoundSpeed(grid.template getY<T>(b));
        });

        double max_sound_speed = equation.getSoundSpeed(grid.template getY<T>(*max_c_s));
        return cfl * dx / max_sound_speed;
    }



    double getCFL()
    {
        return cfl;
    }
    void setCFL(double cfl)
    {
        this->cfl = cfl;
    }
};


template <typename SolvStep, typename CFL_calc,int N>
class Solver
{
private:
    Grid<N> grid_old;
    Grid<N> grid_new;
    double dx;
    CFL_calc CFL;
    double t0;
    double t1;
    SolvStep step;
public:
    Solver(Grid<N> grid, double dx, CFL_calc CFL, SolvStep step_) : grid_old(grid), dx(dx), CFL(CFL) , step(step_)
    {
        grid_new = grid_old;
    }

    void solve()
    {
        double dt = CFL(grid_old, dx);
        t0 = 0;
        t1 = 0;
        while (t1 < 1)
        {
            for (int i = 0; i < grid_old.size(); i++)
            {
                grid_new.setY(i, step(grid_old, dx, dt, i));
            }
            grid_old = grid_new;
            t1 += dt;
        }
    }

};