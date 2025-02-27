#ifndef EQUATION_HPP
#define EQUATION_HPP
#include "Quantity.hpp"
#include <cmath>
#include <tuple>
#include <array>


template <typename T, int N>
class EquationBase {

    public:
    virtual std::array<std::tuple<T,T>,N> AdvectionTerms(T inputval) = 0;
    virtual std::array<T,N> AdvectionCoeffs(T inputval) = 0;
    virtual std::array<std::tuple<T,T>,N> DiffusionTerms(T inputval) = 0;
    virtual double getSoundSpeed(T inputval) = 0;
    virtual T SourceTerms(T inputval) = 0;
    //virtual T PreviousStep(T inputval) = 0;

};


class EulerEquation : public EquationBase<Euler,1>
{
    public:
    std::array<std::tuple<Euler,Euler>,1> AdvectionTerms(Euler inputval) override;
    std::array<Euler,1> AdvectionCoeffs(Euler inputval) override;
    std::array<std::tuple<Euler,Euler>,1> DiffusionTerms(Euler inputval) override;
    Euler SourceTerms(Euler inputval) override;
    double getSoundSpeed(Euler inputval) override;
};

class TransportEquation : public EquationBase<Transport,1>
{
    private:
    double v;

    public:
    std::array<std::tuple<Transport,Transport>,1> AdvectionTerms(Transport inputval) override;
    std::array<std::tuple<Transport,Transport>,1> DiffusionTerms(Transport inputval) override;
    std::array<Transport,1> AdvectionCoeffs(Transport inputval) override;
    Transport SourceTerms(Transport inputval) override;
    double getSoundSpeed(Transport inputval) override;

    TransportEquation(double v) : v(v) {}
    double getv() { return v; }
    void setv(double v) { this->v = v; }

};

#endif
