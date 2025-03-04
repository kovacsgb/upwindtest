#ifndef QUANTITY_HPP
#define QUANTITY_HPP

#include <array>
#include <memory>
#include <iostream>
template <int N>
class Quantity
{
protected:
    std::array<double, N> data;
public:
    Quantity() = default;
    Quantity(Quantity const & q_) : data(q_.data) {}
    Quantity(Quantity &&) = default;
    Quantity(std::array<double, N> data) : data(data) {}
    Quantity(double data) {
        for (int i = 0; i < N; ++i)
        {
            this->data[i] = data;
        }
    }
    Quantity &operator=(Quantity const &) = default;
    Quantity &operator=(Quantity &&) = default;
    virtual ~Quantity() = default;
    

    virtual std::unique_ptr<Quantity<N>> clone()
    {
        return std::make_unique<Quantity<N>>(this->data);
    }
    
    double &operator[](int i)
    {
        return data[i];
    }

    double &operator[](int i) const
    {
        const double &a = data[i];
        return const_cast<double &>(a);
    }



    // Addition operator
    Quantity<N> operator+(const Quantity<N>& other) const
    {
        Quantity<N> result;
        for (int i = 0; i < N; ++i)
        {
            result.data[i] = this->data[i] + other.data[i];
        }
        return result;
    }

    // Subtraction operator
    Quantity<N> operator-(const Quantity<N>& other) const
    {
        Quantity<N> result;
        for (int i = 0; i < N; ++i)
        {
            result.data[i] = this->data[i] - other.data[i];
        }
        return result;
    }

    Quantity<N> operator*(const Quantity<N>& other) const
    {
        Quantity<N> result;
        for (int i = 0; i < N; ++i)
        {
            result.data[i] = this->data[i] * other.data[i];
        }
        return result;
    }

    Quantity<N> operator/(const Quantity<N>& other) const
    {
        Quantity<N> result;
        for (int i = 0; i < N; ++i)
        {
            result.data[i] = this->data[i] / other.data[i];
        }
        return result;
    }

    Quantity<N> operator*(const double& other) const
    {
        Quantity<N> result;
        for (int i = 0; i < N; ++i)
        {
            result.data[i] = this->data[i] * other;
        }
        return result;
    }

    Quantity<N> operator/(const double& other) const
    {
        Quantity<N> result;
        for (int i = 0; i < N; ++i)
        {
            result.data[i] = this->data[i] / other;
        }
        return result;
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
    virtual ~Euler()  = default;
    Euler(double r, double p, double e) : Quantity<3>() {
        data[0] = r;
        data[1] = p;
        data[2] = e;
    }
    Euler(std::array<double, 3> data) : Quantity<3>(data) {}
    Euler(double r) : Quantity<3>() {
        data[0] = r;
        data[1] = r;
        data[2] = r;
    }


    std::unique_ptr<Quantity<3>> clone() override
    {
        return std::make_unique<Euler>(this->data);
    }

    double getrho() const { return data[0]; }
    double getu() const  { return data[1] / data[0]; }
    double getp() const { return (data[2] - 0.5 * data[1] * data[1] / data[0]) * (2./5.); }
    
    Euler operator+(const Euler& other) const
    {
        Euler result;
        for (int i = 0; i < 3; ++i)
        {
            result.data[i] = this->data[i] + other.data[i];
        }
        return result;
    }

    Euler operator-(const Euler& other) const
    {
        Euler result;
        for (int i = 0; i < 3; ++i)
        {
            result.data[i] = this->data[i] - other.data[i];
        }
        return result;
    }

    Euler operator*(const Euler& other) const
    {
        Euler result;
        for (int i = 0; i < 3; ++i)
        {
            result.data[i] = this->data[i] * other.data[i];
        }
        return result;
    }

    Euler operator/(const Euler& other) const
    {
        Euler result;
        for (int i = 0; i < 3; ++i)
        {
            result.data[i] = this->data[i] / other.data[i];
        }
        return result;
    }

    Euler operator*(const double& other) const
    {
        Euler result;
        for (int i = 0; i < 3; ++i)
        {
            result.data[i] = this->data[i] * other;
        }
        return result;
    }

    Euler operator/(const double& other) const
    {
        Euler result;
        for (int i = 0; i < 3; ++i)
        {
            result.data[i] = this->data[i] / other;
        }
        return result;
    }

};



class Transport : public Quantity<1>
{
public:
    Transport() = default;
    Transport(Transport const & Tr) : Quantity<1>(Tr.data) {};
    Transport(Transport &&) = default;
    Transport &operator=(Transport const &) = default;
    Transport &operator=(Transport &&) = default;
    virtual ~Transport() = default;

    Transport(std::array<double, 1> data) : Quantity<1>(data) {}

    Transport(double r) : Quantity<1>() {
        data[0] = r;
    }
    double getrho() const { return data[0]; }

    std::unique_ptr<Quantity<1>> clone() override
    {
       // std::cerr << "Transport clone" << std::endl;
        return std::make_unique<Transport>(this->data[0]);
    }

    Transport operator+(const Transport& other) const
    {
        Transport result;
        result.data[0] = this->data[0] + other.data[0];
        return result;
    }

    Transport operator-(const Transport& other) const
    {
        Transport result;
        result.data[0] = this->data[0] - other.data[0];
        return result;
    }

    Transport operator*(const Transport& other) const
    {
        Transport result;
        result.data[0] = this->data[0] * other.data[0];
        return result;
    }

    Transport operator/(const Transport& other) const
    {
        Transport result;
        result.data[0] = this->data[0] / other.data[0];
        return result;
    }

    Transport operator*(const double& other) const
    {
        Transport result;
        result.data[0] = this->data[0] * other;
        return result;
    }

    Transport operator/(const double& other) const
    {
        Transport result;
        result.data[0] = this->data[0] / other;
        return result;
    }

};

#endif // QUANTITY_HPP