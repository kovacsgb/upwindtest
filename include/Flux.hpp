#ifndef FLUX_HPP
#define FLUX_HPP

#include "Quantity.hpp"
#include "Grid.hpp"
#include "Equation.hpp"
#include <memory>
#include <algorithm>
#include <numeric>
#include <array>
#include <iostream>

template<typename t, int N, int M>
class DonorCell;


template <typename T, int N, int M>
class Flux {

    protected:
    EquationBase<T,N,M>& equations;
    virtual std::array<T,N> F_imhalf(Grid<M>& grid,const std::array<T,N>& coeffs, int i)   = 0;
    virtual std::array<T,N> F_iphalf(Grid<M>& grid,const std::array<T,N>& coeffs, int i)   = 0;

public:
    Flux() = default;
    virtual ~Flux() = default;

    Flux(EquationBase<T,N,M>& equations) : equations(equations) {}

    virtual T operator()(Grid<M>& grid, int i)
    {
        auto coeffs = equations.AdvectionCoeffs(grid.template getY<T>(i));
        auto F_imhalf = this->F_imhalf( grid,coeffs, i);
        auto F_iphalf = this->F_iphalf( grid,coeffs, i);

        std::array<T,N> F_j; 
        size_t j=0;
        std::transform(F_iphalf.begin(), F_iphalf.end(), F_imhalf.begin(), F_j.begin(), [&](const auto& a, const auto& b) {
            /*if(i==2)
            {std::cerr << a[0] << " " << b[0] << std::endl;
            std::cerr << coeffs[j][0] << std::endl;
            std::cerr << (a*coeffs[j]-b*coeffs[j])[0] << std::endl;
            std::cerr << "----------------" << std::endl;}*/
            auto c=a * coeffs[j] - b * coeffs[j];
            j++;
           /* if(i==2)
                {std::cerr << c[0] << std::endl;}*/
            return c;
        });
        return std::accumulate(F_j.begin(), F_j.end(), T{0});
    } 
    friend class DonorCell<T,N,M>;
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
        auto velocity = 0.5*(this->equations.getAdvSpeed(grid.template getY<T>(i))+this->equations.getAdvSpeed(grid.template getY<T>(i-1)));
        for (const auto& coeff : coeffs)
        {
  
            for(auto k=0; k< M; k++)
            {
                const double& c = coeff[k];
                if (velocity < 0)
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
        auto velocity = 0.5*(this->equations.getAdvSpeed(grid.template getY<T>(i))+this->equations.getAdvSpeed(grid.template getY<T>(i+1)));
        for (const auto& coeff : coeffs)
        {
            for(auto k=0;k<M;k++)
            {
                const double& c= coeff[k];
                if (velocity < 0)
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
    Upwind(EquationBase<T,N,M>& equations) : Flux<T,N,M>(equations) {}
    friend class DonorCell<T,N,M>;

};

template <typename T, int N, int M>
class vonLeer : public Flux<T,N,M>
{
    private:

    T delta_q(T q_im1, T q_i, T q_ip1)
    {
        auto delta_q_iphalf = q_ip1 - q_i;
        auto delta_q_imhalf = q_i - q_im1;
        auto delta_factor=delta_q_iphalf * delta_q_imhalf;

        T delta_q=delta_factor/(delta_q_iphalf + delta_q_imhalf);
        for(auto k=0; k<M; k++)
        {
            if(delta_factor[k] <= 0)
            delta_q[k] = 0;
        }
        return delta_q;

    }

    std::array<T,N> shockvell_calc(Grid<M>& grid, int i)
    {
        auto v_i=this->equations.getAdvSpeed(grid.template getY<T>(i));
        auto v_ip1=this->equations.getAdvSpeed(grid.template getY<T>(i+1));

        return std::array<T,N>{0.5*(v_i+v_ip1)};
    }

    std::array<T,N> F_iphalf(Grid<M>& grid,const std::array<T,N>& coeffs, int i)   override
    {
        std::array<T,N> F_iphalf;
        std::array<T,N> shockvel;
        
        T delta_q_left=delta_q(grid.template getY<T>(i-1),grid.template getY<T>(i),grid.template getY<T>(i+1));
        T delta_q_right=delta_q(grid.template getY<T>(i),grid.template getY<T>(i+1),grid.template getY<T>(i+2));

        auto right = this->equations.AdvectionTerms(grid.template getY<T>(i+1)-delta_q_right);
        auto left = this->equations.AdvectionTerms(grid.template getY<T>(i)+delta_q_left);
/*
        auto F_ip1=(this->equations.AdvectionTerms(grid.template getY<T>(i+1)));
        auto F_i=(this->equations.AdvectionTerms(grid.template getY<T>(i)));
        auto q_diff=grid.template getY<T>(i+1)-grid.template getY<T>(i);
        */
        auto v_i=this->equations.getAdvSpeed(grid.template getY<T>(i));
        auto v_ip1=this->equations.getAdvSpeed(grid.template getY<T>(i+1));
        auto qfac = v_i * v_ip1;
/*
        for (auto l=0; l<N;l++)
        {
            for (auto k=0; k<M; k++)
            {
                auto F_diff= std::get<1>(F_ip1[l])-std::get<1>(F_i[l]);
                if (std::abs(q_diff[k]) >=1e-5)
                {
                    shockvel[l][k] = F_diff[k]/q_diff[k];
                }
                else
                {
                    shockvel[l][k] =  this->equations.getAdvSpeed(grid.template getY<T>(i));
                }
            }
        }
*/
        shockvel = shockvell_calc(grid,i);
        size_t j=0;
        for (const auto& coeff : coeffs)
        {
            
            for(auto k=0; k< M; k++)
            {
                const double& c = shockvel[j][k];
                //std::cerr<< qfac[k] << std::endl;
                if(true)
                {
                    if (c < 0 )//&& v_ip1 < 0)
                    {
                        F_iphalf[j][k] =  std::get<1>(right[j])[k];
                    }
                    else //if(c > 0 && v_i > 0)
                    {
                        F_iphalf[j][k] = std::get<1>(left[j])[k];
                    }
                }
                else
                {
                    F_iphalf[j][k] = 0;
                }
            }
            j+=1;
        }
        return F_iphalf;
    }

    std::array<T,N> F_imhalf(Grid<M>& grid,const std::array<T,N>& coeffs, int i) override
    {
        std::array<T,N> F_imhalf = F_iphalf(grid,coeffs,i-1);
        return F_imhalf;
    }   

    public:
    vonLeer(EquationBase<T,N,M>& equations) : Flux<T,N,M>(equations) {}
    friend class DonorCell<T,N,M>;
};

template <typename T, int N, int M>
class Central : public Flux<T,N,M>
{
    std::array<T,N> F_iphalf(Grid<M>& grid,const std::array<T,N>& coeffs, int i)   override
    {
        size_t j=0;
        std::array<T,N> Fp1half;
        for (const auto& Val :this->equations.AdvectionTerms(grid.template getY<T>(i+1)))
        {
            Fp1half[j] = std::get<1>(Val)*0.5;
            j++; 
        }
        return Fp1half;
    }
    std::array<T,N> F_imhalf(Grid<M>& grid,const std::array<T,N>& coeffs, int i) override
    {
        return this->F_iphalf(grid,coeffs,i-2);
    }
    public:
    Central(EquationBase<T,N,M>& equations) : Flux<T,N,M>(equations) {}

    friend class DonorCell<T,N,M>;
};

template<typename T, int N, int M>
class DonorCell : public Flux<T,N,M>
{
    private:
    Central<T,N,M> CS;
    Flux<T,N,M>& BS;
    std::array<T,N> F_iphalf(Grid<M>& grid,const std::array<T,N>& coeffs, int i)
    {
        auto v_i = this->equations.getAdvSpeed(grid.template getY<T>(i));
        auto v_ip1 =  this->equations.getAdvSpeed(grid.template getY<T>(i+1));
        auto v_iphalf = 0.5*(v_i + v_ip1);
        //auto cSound = 0.5*( this->equations.getSoundSpeed(grid.template getY<T>(i))+
        //this->equations.getSoundSpeed(grid.template getY<T>(i+1)));
        auto cSound = this->equations.getSoundSpeed(grid.template getY<T>(i));

        double theta = v_i/cSound;
        theta = (theta < 0.2) ? 0.2 : ((theta > 1) ? 1 : theta);
        auto CS_arr=CS.F_iphalf(grid,coeffs,i);
        auto BS_arr=BS.F_iphalf(grid,coeffs,i);

        std::array<T,N> Final_F;
        for (auto j=0;j<N;j++)
        {
            Final_F[j] = T{1-theta} * CS_arr[j] + T{theta} * BS_arr[j];
        }

        return Final_F;

    }
    std::array<T,N> F_imhalf(Grid<M>& grid,const std::array<T,N>& coeffs, int i)
    {
        auto v_i =  this->equations.getAdvSpeed(grid.template getY<T>(i));
        auto v_im1 =  this->equations.getAdvSpeed(grid.template getY<T>(i-1));
        auto v_imhalf = 0.5*(v_i + v_im1);
        //auto cSound = 0.5*( this->equations.getSoundSpeed(grid.template getY<T>(i));//)+
        //this->equations.getSoundSpeed(grid.template getY<T>(i-1)));
        auto cSound =this->equations.getSoundSpeed(grid.template getY<T>(i));

        double theta = v_i/cSound;
        theta = (theta < 0.2) ? 0.2 : ((theta > 1) ? 1 : theta);

        auto CS_arr=CS.F_imhalf(grid,coeffs,i);
        auto BS_arr=BS.F_imhalf(grid,coeffs,i);

        std::array<T,N> Final_F;
        for (auto j=0;j<N;j++)
        {
            Final_F[j] = T{1-theta} * CS_arr[j] + T{theta} * BS_arr[j];
        }

        return Final_F;
    }
    public:
    DonorCell(EquationBase<T,N,M>& equations,Flux<T,N,M>& bs) : Flux<T,N,M>(equations), CS(equations), BS(bs)  {}
};


class FluxWithViscosity : public Flux<Euler,1,3>
{
    private:
    Flux<Euler,1,3>& Calc;
    double xi;

    std::array<Euler,1> F_iphalf(Grid<3>& grid,const std::array<Euler,1>& coeffs, int i)
    {
        return std::array<Euler,1>{0};
    }

    std::array<Euler,1> F_imhalf(Grid<3>& grid,const std::array<Euler,1>& coeffs, int i)
    {
        return std::array<Euler,1>{0};
    }


    public:
        FluxWithViscosity(Flux<Euler,1,3>& Calc,double xi=0.1) :Flux<Euler,1,3>(equations), Calc(Calc), xi(xi) {}

        double getxi() {return this->xi;}
        void setxi(double xx){this->xi = xx;}

        Euler operator()(Grid<3>& grid, int i) override
        {
            auto F_j = Calc(grid,i);
            /*double v_ip2=grid.getY<Euler>(i+1).getu();
            double v_ip1=grid.getY<Euler>(i+1).getu();
            double v_im2=grid.getY<Euler>(i-2).getu();
            double v_i= grid.getY<Euler>(i).getu();
            double v_im1= grid.getY<Euler>(i-1).getu();
            double rho_i=grid.getY<Euler>(i).getrho();
            double rho_im1=grid.getY<Euler>(i-1).getrho();
            double rho_im2=grid.getY<Euler>(i-2).getrho();
            double rho_ip1=grid.getY<Euler>(i-2).getrho();

            double v_iphalf=0.5*(v_i+v_ip1);
            double v_imhalf=0.5*(v_im1+v_i);


            double  PI_i, PI_im1,PI_ip1;

            if(v_ip1-v_im1<0)
            {
                PI_i =1* xi * rho_i * (v_ip1-v_im1)*(v_ip1-v_im1);
            }
            else
            {
                PI_i=0;
            }
            if(v_im1-v_im2<0)
            {
                PI_im1=1 * xi *rho_im1* (v_i-v_im2)*(v_i-v_im2);
            }
            else
            {
                PI_im1=0;
            }
            if(v_ip2-v_i<0)
            {
                PI_ip1=0.25 * xi *rho_ip1* (v_ip2-v_i)*(v_ip2-v_i);
            }
            else
            {
                PI_ip1=0;
            }

            double PI_iphalf=PI_i;//0.5*(PI_i+PI_ip1);
            double PI_imhalf=PI_im1;//0.5*(PI_i+PI_im1);

            //F_j[0] -= 
            F_j[1] -=  PI_i; // PI_iphalf - PI_imhalf;
            F_j[2] -=   PI_i*v_i;//v_iphalf* PI_iphalf - v_imhalf * PI_imhalf;

            return F_j;*/



        }
};


#endif