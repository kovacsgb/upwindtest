#include "Equation.hpp"
#include "Quantity.hpp"

std::array<std::tuple<Euler,Euler>,1> EulerEquation::AdvectionTerms(Euler inputval)
{
    Euler terms;
    terms[0]=inputval[1];
    terms[1]=2/3.*(inputval[1]*inputval[1]/inputval[0])+inputval[2];
    terms[3]=5./3.*inputval[1]*inputval[2]/inputval[0]-1/3.*inputval[1]*inputval[1]/inputval[0];

    Euler coeff{1,1,1};
    return {std::make_tuple(coeff,terms)};
}

std::array<std::tuple<Euler,Euler>,1> EulerEquation::DiffusionTerms(Euler inputval)
{
    Euler terms;
    terms[0]=0;
    terms[1]=0;
    terms[2]=0;

    Euler coeff{0,0,0};
    return {std::make_tuple(coeff,terms)};
}

Euler EulerEquation::SourceTerms(Euler inputval)
{
    Euler terms;
    terms[0]=0;
    terms[1]=0;
    terms[2]=0;

    return terms;
}

std::array<Euler,1> EulerEquation::AdvectionCoeffs(Euler inputval)
{
    Euler coeff{1,1,1};

    return std::array<Euler,1>{coeff};
}

double EulerEquation::getSoundSpeed(Euler inputval)
{
    return sqrt(1.4*inputval.getp()/inputval.getrho());
}

/*
std::array<Euler,1> EulerEquation::AdvectionCoeffs(Euler inputval)
{
    Euler coeff{1,1,1};

    return std::array<Euler,1>{coeff};

}*/