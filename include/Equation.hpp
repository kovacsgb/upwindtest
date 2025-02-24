#include "Quantity.hpp"
#include <cmath>
#include <tuple>
#include <array>


template <typename T, int N>
class EquationBase {

    public:
    virtual std::array<std::tuple<T,T>,N> AdvectionTerms(T inputval) = 0;
    virtual std::array<std::tuple<T,T>,N> DiffusionTerms(T inputval) = 0;
    virtual T SourceTerms(T inputval) = 0;
    virtual T PreviousStep(T inputval) = 0;

};


class EulerEquation : public EquationBase<Euler,1>
{
    public:
    std::array<std::tuple<Euler,Euler>,1> AdvectionTerms(Euler inputval) override;
    std::array<std::tuple<Euler,Euler>,1> DiffusionTerms(Euler inputval) override;
    Euler SourceTerms(Euler inputval) override;
};