#ifndef FLUX_HPP
#define FLUX_HPP

#include "Quantity.hpp"
#include "Grid.hpp"
#include "Equation.hpp"
#include <memory>
#include <algorithm>
#include <numeric>
#include <array>

template <typename T, int N, int M>
class Flux {

    protected:
    EquationBase<T,N>& equations;
    virtual std::array<T,N> F_imhalf(Grid<M>& grid,const std::array<T,N>& coeffs, int i)   = 0;
    virtual std::array<T,N> F_iphalf(Grid<M>& grid,const std::array<T,N>& coeffs, int i)   = 0;

public:
    Flux() = default;
    virtual ~Flux() = default;

    Flux(EquationBase<T,N>& equations) : equations(equations) {}

    virtual T operator()(Grid<M>& grid, int i)
    {
        auto coeffs = equations.AdvectionCoeffs(grid.template getY<T>(i));
        auto F_imhalf = this->F_imhalf( grid,coeffs, i);
        auto F_iphalf = this->F_iphalf( grid,coeffs, i);

        std::array<T,N> F_j; 
        size_t j=0;
        std::transform(F_iphalf.begin(), F_iphalf.end(), F_imhalf.begin(), F_j.begin(), [&](const auto& a, const auto& b) {
            return a * coeffs[j] - b * coeffs[j++];
        });
        return std::accumulate(F_j.begin(), F_j.end(), T{0});
    }
};

template <typename T, int N, int M>
class Upwind : public Flux<T,N,M>
{

    private:
    std::array<T,N> F_imhalf(Grid<M>& grid,const std::array<T,N>& coeffs, int i)   override
    {
        std::array<T,N> F_imhalf;
        size_t j=0;
        auto left = this->equations.AdvectionTerms(grid.template getY<T>(i-1));
        auto right = this->equations.AdvectionTerms(grid.template getY<T>(i));

        for (const auto& coeff : coeffs)
        {
  
            for(auto k=0; k< M; k++)
            {
                const double& c = coeff[k];
                if (c > 0)
                {
                    F_imhalf[j][k] =  std::get<1>(right[j])[k];
                }
                else
                {
                    F_imhalf[j][k] = std::get<1>(left[j])[k];
                }
            }
            j+=1;
        }
        return F_imhalf;
    }


    std::array<T,N> F_iphalf(Grid<M>& grid,const std::array<T,N>& coeffs, int i) override
    {
        std::array<T,N> F_iphalf;
        size_t j=0;
        auto left = this->equations.AdvectionTerms(grid.template getY<T>(i));
        auto right = this->equations.AdvectionTerms(grid.template getY<T>(i+1));

        for (const auto& coeff : coeffs)
        {
            for(auto k=0;k<M;k++)
            {
                const double& c= coeff[k];
                if (c > 0)
                {
                    F_iphalf[j][k] =  std::get<1>(right[j])[k];
                }
                else
                {
                    F_iphalf[j][k] = std::get<1>(left[j])[k];
                }
            }
            j+=1;
        }
        return F_iphalf;
    }



    public:
    Upwind(EquationBase<T,N>& equations) : Flux<T,N,M>(equations) {}

};
#endif