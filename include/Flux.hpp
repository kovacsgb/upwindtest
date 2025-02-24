#include "Quantity.hpp"
#include "Grid.hpp"
#include "Equation.hpp"
#include <memory>

template <typename T, int N>
class Flux {

    protected:
    std::shared_ptr<EquationBase<T,N>> equations;

public:
    Flux() = default;
    virtual ~Flux() = default;

    virtual T operator()(Grid<N> grid, int i) const = 0;
};

template<typename T, int N>
class Upwind : public Flux<T,N> {
    public:
    virtual T calcFlux(Grid <N> grid, int iphalf) = 0;
    virtual T operator()(Grid<N> grid, int i)
    {
        T F_iphalf = calcFlux(grid, i);
        T F_imhalf = calcFlux(grid, i);

        
    }

};

template <typename T, int N>
class BasicUpwind : public Upwind<T,N> {
};
