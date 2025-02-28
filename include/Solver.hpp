#ifndef SOLVER_HPP
#define SOLVER_HPP 

#include "Grid.hpp"
//#include "Solverstep.hpp"
#include "Equation.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>
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
        /*
        auto max_c_s = std::max_element(grid.begin(), grid.end(), [&](const auto& a, const auto& b) {
            return equation.getSoundSpeed(grid.template getY<T>(a)) < equation.getSoundSpeed(grid.template getY<T>(b));
        });*/


        //double max_sound_speed = equation.getSoundSpeed(grid.template getY<T>(*max_c_s));
        
        double max_sound_speed = 0;
        for (int i = grid.startid(); i < grid.size(); i++)
        {

            const T& q = grid.template getY<T>(i);
            std::cerr << q[0] << std::endl;
            max_sound_speed = std::max(max_sound_speed, equation.getSoundSpeed(q));
        }
        
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


template <typename T, typename SolvStep, typename CFL_calc,int N>
class Solver
{
private:
    Grid<N> grid_old;
    Grid<N> grid_new;
    double dx;
    CFL_calc& CFL;
    double t0;
    double t1;
    SolvStep& step;
public:
    Solver(const Grid<N>& grid, double dx, CFL_calc& CFL, SolvStep& step_) :
     grid_old(grid), grid_new(grid), dx(dx), CFL(CFL) , t0(0), t1(1), step(step_)
    //Solver(Grid<N> grid, double dx, CFL_calc CFL, SolvStep step_)
     {
        std::cerr << "Solver constructor" << std::endl;
        /*grid_old = grid;
        grid_new = grid;
        this->dx = dx;
        this->CFL = CFL;
        this->step = step_;
        t0=0;
        t1=1;
*/
        //this->grid_new  = grid_old;
    }

    void solve()
    {
        

        double dt = CFL(grid_old, dx);
        double t=t0;
        std::ofstream output_transport{"out_ts.txt"};
        while (t < t1)
        {
            double dt = CFL(grid_old, dx);
            std::cerr << "Time: " << t0 << t1 << std::endl;
            for (int i = grid_old.startid(); i < grid_old.size(); i++)
            {
                grid_new.template getY<T>(i) = step(grid_old, dx, dt, i);
            }
            

            
            for (int i = 0; i < grid_new.totalsize(); i++)
            {
                T transport = grid_new.template getY<T>(i);
                output_transport << grid_new.getX(i) << " " << transport[0]  << std::endl;
            }

            output_transport << std::endl;

            //grid_new.updateBoundary(step.getEquation());
            grid_old = grid_new;
            t += dt;
        }
    }

    void setT0(double t0)
    {
        this->t0 = t0;
    }
    void setT1(double t1)
    {
        this->t1 = t1;
    }
    double getT0()
    {
        return t0;
    }
    double getT1()
    {
        return t1;
    }
    Grid<N>& getGrid()
    {
        return grid_old;
    }
    void setGrid(Grid<N>& grid)
    {
        grid_old = grid;
    }
    double getDx()
    {
        return dx;
    }
    void setDx(double dx)
    {
        this->dx = dx;
    }

};

#endif // SOLVER_HPP