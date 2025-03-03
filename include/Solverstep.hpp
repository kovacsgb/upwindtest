#ifndef SOLVERSTEP_HPP
#define SOLVERSTEP_HPP

#include "Grid.hpp"
#include "Quantity.hpp"
#include "Equation.hpp"
#include "Flux.hpp"

template <typename T,int N, int M>
class Solverstep
{
    protected:
    Flux<T,N,M> &Fluxterm;
    EquationBase<T,N,M>& equation;
    public:
    Solverstep(Flux<T,N,M>& Fluxterm, EquationBase<T,N,M>& equation) : Fluxterm(Fluxterm), equation(equation) {}
    virtual void operator()(Grid<M>& grid_old,Grid<M>& grid_new, double dx, double dt) = 0;
    virtual EquationBase<T,N,M>& getEquation() { return equation; }
};


template <typename T,int N, int M>
class explicitStep : public Solverstep<T,N,M>
{

    public:

    explicitStep(EquationBase<T,N,M>& equation, Flux<T,N,M>& Fluxterm) : Solverstep<T,N,M>{Fluxterm,equation} {}

    void operator()(Grid<M>& grid_old, Grid<M>& grid_new, double dx, double dt) override
    {
        for (int i = grid_old.startid(); i < grid_old.size(); i++)
        {
            T old = grid_old.template getY<T>(i);
            T new_val = old -this->Fluxterm(grid_old, i) *  (dt /dx )-this->equation.SourceTerms(old);
            grid_new.template getY<T>(i) = new_val;
        }
        
        
        /*if (i == 2)
        {
            std::cerr << "Flux_term:" << this->Fluxterm(grid_old, i)[0] << std::endl;
            std::cerr << "Old value:" << old[0] << std::endl;  
            std::cerr << "New value:" << new_val[0] << std::endl;
            std::cerr << "dt/dx" << dt/dx << std::endl;
            
        }*/
        //return new_val;
    }

};


template <typename T,int N, int M>
class Runge2Step : public Solverstep<T,N,M>
{
    private:
    Grid<M> grid_middle;

    public:
    Runge2Step(EquationBase<T,N,M>& equation, Flux<T,N,M>& Fluxterm,const Grid<M>& template_grid) : Solverstep<T,N,M>{Fluxterm,equation}, grid_middle(template_grid) {}

    void operator()(Grid<M>& grid_old,Grid<M>& grid_new, double dx, double dt) override
    {
        for (int i = grid_old.startid(); i < grid_old.size(); i++)
        {
            T old = grid_old.template getY<T>(i);
            T middle = old - this->Fluxterm(grid_old, i) * 0.5 - (dt /dx )-this->equation.SourceTerms(old);
            grid_middle.template getY<T>(i) = middle;
        }

        grid_middle.updateBoundary(this->equation);

        for (int i = grid_old.startid(); i < grid_old.size(); i++)
        {
            T old = grid_old.template getY<T>(i);
            T middle = grid_middle.template getY<T>(i);
            T new_val = old - this->Fluxterm(grid_middle, i) * (dt /dx )-this->equation.SourceTerms(middle);
            grid_new.template getY<T>(i) = new_val;
        }

    }
};

#endif // SOLVERSTEP_HPP