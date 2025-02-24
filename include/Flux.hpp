#include "Quantity.hpp"



class Flux {
public:
    Flux() = default;
    virtual ~Flux() = default;

    virtual double operator()(const double& x, const double& t) const = 0;
};