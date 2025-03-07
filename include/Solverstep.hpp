#ifndef SOLVERSTEP_HPP
#define SOLVERSTEP_HPP

#include "Grid.hpp"
#include "Quantity.hpp"
#include "Equation.hpp"
#include "Flux.hpp"
#include <fstream>
#include <string>


template <int N>
struct CheckIfThree
{
    static constexpr bool value = false;
};

// Specialization for N = 3
template <>
struct CheckIfThree<3>
{
    static constexpr bool value = true;
};

template <int N>
bool isThree()
{
    return CheckIfThree<N>::value;
}

template <int M,class T>
 T apply_visc(Grid<M> grid_old, int i)
{
    return T{0};
}

template <>
 Euler apply_visc(Grid<3> grid_old, int i)
{

    Euler old_i = grid_old.template getY<Euler>(i);
    Euler old_ip1 = grid_old.template getY<Euler>(i+1);
    Euler old_ip2 = grid_old.template getY<Euler>(i+1);
    Euler old_im1 = grid_old.template getY<Euler>(i-1);
    Euler old_im2 = grid_old.template getY<Euler>(i-1);
/*
    Euler vis_val{0};
    double rho_iphalf = 0.5*(old_ip1.getrho()+old_i.getrho());
    double rho_imhalf = 0.5*(old_i.getrho()+old_im1.getrho());

    double deltaU_iphalf=0.5*(old_ip2.getu()-old_i.getu());
    double deltaU_imhalf=0.5*(old_i.getu()-old_im2.getu());
    double k=1;

    double PI_iphalf = k*rho_iphalf*deltaU_iphalf*deltaU_iphalf;
    double PI_imhalf = k*rho_imhalf*deltaU_imhalf*deltaU_iphalf;


    double PI_i = k*old_i.getrho()*0.25*(old_ip1.getu()-old_im1.getu())*(old_ip1.getu()-old_im1.getu());
    double E_part = PI_i*(old_ip1.getu()-old_im1.getu());

        vis_val =(old_i.getu()-old_im1.getu()) > 0 ? Euler{0,PI_iphalf-PI_imhalf,E_part} : Euler(0);

        if(vis_val[1] !=0) std::cerr << vis_val[1] <<" "<< vis_val[2] << std::endl;*/
    Euler Val_ip1 = (old_ip1-old_i);
    Euler Val_i = (old_i - old_im1);

    auto absval = [](Euler val){return Euler{std::abs(val[0]),
        std::abs(val[1]),
        std::abs(val[2])};};

    Euler abs_ip1 = absval(Val_ip1);
    Euler abs_i = absval(Val_i);

    Euler vis_val = Euler{0,1,1}* (abs_ip1 * Val_ip1 - abs_i * Val_i); 
    
    //if (Val_ip1[1] <= 0) vis_val = Euler{0,0,0}; 
    return vis_val;
}



template <typename T,int N, int M>
class Solverstep
{
    protected:
    Flux<T,N,M> &Fluxterm;
    EquationBase<T,N,M>& equation;
    public:
    Solverstep(Flux<T,N,M>& Fluxterm, EquationBase<T,N,M>& equation) : Fluxterm(Fluxterm), equation(equation) {}
    virtual void operator()(Grid<M>& grid_old,Grid<M>& grid_new, double dx, double dt) = 0;
    virtual EquationBase<T,N,M>& getEquation() { return equation; }
};

template <typename T,int N, int M>
class explicitStep : public Solverstep<T,N,M>
{
    private:
    std::ofstream debug_out{"out_fluxes"};

    double xi=1;
    public:

    explicitStep(EquationBase<T,N,M>& equation, Flux<T,N,M>& Fluxterm) : Solverstep<T,N,M>{Fluxterm,equation} {}

    void operator()(Grid<M>& grid_old, Grid<M>& grid_new, double dx, double dt) override
    {
        for (int i = grid_old.startid(); i < grid_old.size(); i++)
        {
            T old = grid_old.template getY<T>(i);
            T fluxterm = this->Fluxterm(grid_old, i);
            /*debug_out << grid_old.getX(i);
            for(auto j=0;j<M;j++)
            {
                debug_out << " " << fluxterm[j];
            }*/
            debug_out << std::endl;
            T new_val = old - fluxterm *  (dt /dx )-this->equation.SourceTerms(old);
            grid_new.template getY<T>(i) = new_val;
        }
        

        //artificial viscosity
        if(isThree<M>())
        {
            grid_old = grid_new;
            grid_old.updateBoundary(this->equation);
            for (int i = grid_old.startid(); i < grid_old.size(); i++)
            {

                T old_i = grid_old.template getY<T>(i);
                T vis_val = apply_visc<M,T>(grid_old,i);

                T new_val = old_i +  vis_val  * (xi*dt/dx);
                grid_new.template getY<T>(i) = new_val;
            }
        }
        


        debug_out << std::endl;
        
        /*if (i == 2)
        {
            std::cerr << "Flux_term:" << this->Fluxterm(grid_old, i)[0] << std::endl;
            std::cerr << "Old value:" << old[0] << std::endl;  
            std::cerr << "New value:" << new_val[0] << std::endl;
            std::cerr << "dt/dx" << dt/dx << std::endl;
            
        }*/
        //return new_val;
    }
    void setXi(double new_xi) {
        xi = new_xi;
    }

    double getXi() const {
        return xi;
    }

};


template <typename T,int N, int M>
class Runge2Step : public Solverstep<T,N,M>
{
    private:
    Grid<M> grid_middle;

    public:
    Runge2Step(EquationBase<T,N,M>& equation, Flux<T,N,M>& Fluxterm,const Grid<M>& template_grid) : Solverstep<T,N,M>{Fluxterm,equation}, grid_middle(template_grid) {}

    void operator()(Grid<M>& grid_old,Grid<M>& grid_new, double dx, double dt) override
    {
        for (int i = grid_old.startid(); i < grid_old.size(); i++)
        {
            T old = grid_old.template getY<T>(i);
            T middle = old - this->Fluxterm(grid_old, i) * 0.5 - (dt /dx )-this->equation.SourceTerms(old);
            grid_middle.template getY<T>(i) = middle;
        }

        if(isThree<M>())
        {
            grid_new = grid_middle;
            grid_new.updateBoundary(this->equation);
            for (int i = grid_new.startid(); i < grid_new.size(); i++)
            {

                T old_i = grid_new.template getY<T>(i);
                T vis_val = apply_visc<M,T>(grid_new,i);

                T new_val = old_i +  vis_val  * (0.5*dt/dx);
                grid_middle.template getY<T>(i) = new_val;
            }
        }


        grid_middle.updateBoundary(this->equation);

        for (int i = grid_old.startid(); i < grid_old.size(); i++)
        {
            T old = grid_old.template getY<T>(i);
            T middle = grid_middle.template getY<T>(i);
            T new_val = old - this->Fluxterm(grid_middle, i) * (dt /dx )-this->equation.SourceTerms(middle);
            grid_new.template getY<T>(i) = new_val;
        }

        if(isThree<M>())
        {
            grid_old = grid_new;
            grid_old.updateBoundary(this->equation);
            for (int i = grid_old.startid(); i < grid_old.size(); i++)
            {

                T old_i = grid_old.template getY<T>(i);
                T vis_val = apply_visc<M,T>(grid_old,i);

                T new_val = old_i +  vis_val  * (0.5*dt/dx);
                grid_new.template getY<T>(i) = new_val;
            }
        }
        
    }
};

#endif // SOLVERSTEP_HPP