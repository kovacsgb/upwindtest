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
    EquationBase<T,M>& equation;
    public:
    Solverstep(Flux<T,N,M>& Fluxterm, EquationBase<T,M>& equation) : Fluxterm(Fluxterm), equation(equation) {}
    virtual T operator()(Grid<N>& grid_old, double dx, double dt, int i) = 0;
};


template <typename T,int N, int M>
class explicitStep : public Solverstep<T,N,M>
{
    public:

    explicitStep(EquationBase<T,M>& equation, Flux<T,N,M>& Fluxterm) : Solverstep<T,N,M>{Fluxterm,equation} {}

    T operator()(Grid<N>& grid_old, double dx, double dt, int i) override
    {
        T old = grid_old.template getY<T>(i);
        T new_val = old +this->Fluxterm(grid_old, i) *  (dt /dx )+this->equation.SourceTerms(old);
        return new_val;
    }

};

#endif // SOLVERSTEP_HPP