#include "Equation.hpp"
#include <cmath>

std::array<std::tuple<Transport, Transport>, 1> BurgersEquation::AdvectionTerms(Transport inputval) {
    // Define the advection terms for Burgers' equation
    Transport term = inputval * inputval / 2.0;
    return {std::make_tuple(Transport{1}, term)};
}

std::array<Transport, 1> BurgersEquation::AdvectionCoeffs(Transport inputval) {
    // Define the advection coefficients for Burgers' equation
    return {1};
}

std::array<std::tuple<Transport, Transport>, 1> BurgersEquation::DiffusionTerms(Transport inputval) {
    // Define the diffusion terms for Burgers' equation
    Transport term = 0.0; // Assuming no diffusion for simplicity
    return {std::make_tuple(term, term)};
}

Transport BurgersEquation::SourceTerms(Transport inputval) {
    // Define the source terms for Burgers' equation
    return 0.0; // Assuming no source term for simplicity
}

double BurgersEquation::getSoundSpeed(const Transport inputval) {
    // Define the sound speed for Burgers' equation
    return std::abs(inputval[0]);
}

double BurgersEquation::getAdvSpeed(const Transport inputval) {
    // Define the advection speed for Burgers' equation
    return inputval[0];
}