#ifndef EQUATION_HPP
#define EQUATION_HPP
#include "Quantity.hpp"
#include "Grid.hpp"
#include <cmath>
#include <tuple>
#include <array>

template <int N> struct Grid;

template <typename T, int N,int M>
class EquationBase {

    public:
    virtual std::array<std::tuple<T,T>,N> AdvectionTerms(T inputval) = 0;
    virtual std::array<T,N> AdvectionCoeffs(T inputval) = 0;
    virtual double getAdvSpeed(const T inputval) = 0;
    virtual std::array<std::tuple<T,T>,N> DiffusionTerms(T inputval) = 0;
    
    virtual void updateBoundary(Grid<M>& grid)
    {
        grid.template getY<T>(1) = grid.template getY<T>(2);
        grid.template getY<T>(0) = grid.template getY<T>(1);
        grid.template getY<T>(grid.totalsize()-2) = grid.template getY<T>(grid.totalsize()-3);
        grid.template getY<T>(grid.totalsize()-1) = grid.template getY<T>(grid.totalsize()-2);
    };
    virtual double getSoundSpeed(const T inputval) = 0;
    virtual T SourceTerms(T inputval) = 0;
    //virtual T PreviousStep(T inputval) = 0;

};


class EulerEquation : public EquationBase<Euler,1,3>
{
    public:
    std::array<std::tuple<Euler,Euler>,1> AdvectionTerms(Euler inputval) override;
    std::array<Euler,1> AdvectionCoeffs(Euler inputval) override;
    std::array<std::tuple<Euler,Euler>,1> DiffusionTerms(Euler inputval) override;
    Euler SourceTerms(Euler inputval) override;
    double getSoundSpeed(const Euler inputval) override;
    double getAdvSpeed(const Euler inputval) override;
};

class TransportEquation : public EquationBase<Transport,1,1>
{
    private:
    double v;

    public:
    std::array<std::tuple<Transport,Transport>,1> AdvectionTerms(Transport inputval) override;
    std::array<std::tuple<Transport,Transport>,1> DiffusionTerms(Transport inputval) override;
    std::array<Transport,1> AdvectionCoeffs(Transport inputval) override;
    Transport SourceTerms(Transport inputval) override;
    double getSoundSpeed(const Transport inputval) override;
    double getAdvSpeed(const Transport inputval) override;
    TransportEquation(double v) : v(v) {}
    double getv() { return v; }
    void setv(double v) { this->v = v; }

};


class BurgersEquation : public EquationBase<Transport,1,1>
{
    public:
    std::array<std::tuple<Transport,Transport>,1> AdvectionTerms(Transport inputval) override;
    std::array<Transport,1> AdvectionCoeffs(Transport inputval) override;
    std::array<std::tuple<Transport,Transport>,1> DiffusionTerms(Transport inputval) override;
    Transport SourceTerms(Transport inputval) override;
    double getSoundSpeed(const Transport inputval) override;
    double getAdvSpeed(const Transport inputval) override;
};

/*

class EulerWithDiffusion : public EquationBase<Euler,1,3>
{
    private:
    std::vector<double> PI;
    double xi=0.1;
    void calc_PI(Grid<3>& grid);
    public:

    EulerWithDiffusion(int N) : PI(N,0) {}

    std::array<std::tuple<Euler,Euler>,1> AdvectionTerms(Euler inputval) override;
    std::array<Euler,1> AdvectionCoeffs(Euler inputval) override;
    std::array<std::tuple<Euler,Euler>,1> DiffusionTerms(Euler inputval) override;

    Euler SourceTerms(Euler inputval) override;
    double getSoundSpeed(const Euler inputval) override;
    double getAdvSpeed(const Euler inputval) override;
    double get_xi() { return xi; }
    void set_xi(double xi) { this->xi = xi; }
};
*/
#endif
