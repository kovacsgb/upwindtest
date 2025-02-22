#ifndef QUANTITY_HPP
#define QUANTITY_HPP

#include <array>

template <int N>
class Quantity
{
protected:
    std::array<double, N> data;
public:
    Quantity() = default;
    Quantity(Quantity const &) = default;
    Quantity(Quantity &&) = default;
    Quantity &operator=(Quantity const &) = default;
    Quantity &operator=(Quantity &&) = default;
    virtual ~Quantity() = default;

    double &operator[](int i)
    {
        return data[i];
    }
};

class Euler : public Quantity<3>
/*
Euler quantities for grid:
0 - density
1 - momentum
2 - total energy 
*/
{
public:
    Euler() = default;
    Euler(Euler const &) = default;
    Euler(Euler &&) = default;
    Euler &operator=(Euler const &) = default;
    Euler &operator=(Euler &&) = default;
    virtual ~Euler() = default;
    double getrho() { return data[0]; }
    double getu() { return data[1] / data[0]; }
    double getp() { return (data[2] - 0.5 * data[1] * data[1] / data[0]) * (2./3.); }
};

#endif // QUANTITY_HPP