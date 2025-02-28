#include "Equation.hpp"
#include "Quantity.hpp"

std::array<std::tuple<Transport,Transport>,1> TransportEquation::AdvectionTerms(Transport inputval)
{
    std::array<std::tuple<Transport,Transport>,1> result;
    result[0] = std::make_tuple(this->getv(), inputval);
    return result;
}

std::array<std::tuple<Transport,Transport>,1> TransportEquation::DiffusionTerms(Transport inputval)
{
    std::array<std::tuple<Transport,Transport>,1> result;
    result[0] = std::make_tuple(0.0, inputval);
    return result;
}

std::array<Transport,1> TransportEquation::AdvectionCoeffs(Transport inputval)
{
    return std::array<Transport,1>{this->getv()};
}


Transport TransportEquation::SourceTerms(Transport inputval)
{
    return 0.0;
}

double TransportEquation::getSoundSpeed(const Transport inputval)
{
    return this->getv();
}