#include "Sod_exact.hpp"
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <tuple>
#include <algorithm>

#include <functional>

double P3_func(double P3, std::vector<double> params)
{
    double gamma = params[0];
    double P1 = params[1];
    double rho1 = params[2];
    double P5 = params[3];
    double rho5 = params[4];

    double Gamma = (gamma-1)/(gamma+1);
    double Beta = (gamma-1)/2/gamma;

    double u4 = (P3-P5)*std::sqrt((1-Gamma)/(rho5*(P3+Gamma*P5)));
    double u3 = (std::pow(P1,Beta)-std::pow(P3,Beta))*std::sqrt((1-Gamma*Gamma)*std::pow(P1,1./gamma)
    /(Gamma*Gamma*rho1));

    return u3-u4;

}

double secant_method(std::function<double(double,std::vector<double>)> f,double left, double right, std::vector<double> params)
{
    double tol = 1e-6;
    int max_iter = 100;
    double f_left = f(left, params);
    double f_right = f(right, params);

    if (std::abs(f_left) < tol) return left;
    if (std::abs(f_right) < tol) return right;

    for (int i = 0; i < max_iter; ++i) {
        double mid = right - f_right * (right - left) / (f_right - f_left);
        double f_mid = f(mid, params);

        if (std::abs(f_mid) < tol) return mid;

        left = right;
        f_left = f_right;
        right = mid;
        f_right = f_mid;
    }

    throw std::runtime_error("Secant method did not converge");
}


std::vector<std::array<double,3>> SOD_exact(std::vector<double> X_arr,std::vector<double> params,double t, double gamma)
{
    double rhol = params[0];
    double vl = params[1];
    double Pl = params[2];
    double rhor = params[3];
    double vr = params[4];
    double Pr = params[5];
    double x_center = params[6];

    //sound speeds
    double cs_l = std::sqrt(gamma*Pl/rhol);
    double cs_r = std::sqrt(gamma*Pr/rhor);

    // Helper const
    double Gamma = (gamma-1)/(gamma+1);
    double Beta = (gamma-1)/2/gamma;
    //Calculate iteratively P3

    auto cs_calc = [&gamma](double p,double rho){return std::sqrt(gamma*p/rho);};


    double P3 = secant_method(P3_func,Pl,Pr,{gamma,Pl,rhol,Pr,rhor});


    double rho3 = rhol*std::pow(P3/Pl,1/gamma);
    
    double v3 = vr + (P3-Pr)/std::sqrt(0.5*rhor*((gamma+1)*P3+(gamma-1)*Pr));

    double c3 = cs_calc(P3,rho3);

    double P4 = P3;

    double rho4 = rhor* (P4+ Gamma * Pr)/(Pr+ Gamma* P4);

    double v4 = v3;

    double c4 = cs_calc(P4,rho4);

    //rarefaction fan functions

    auto v2   = [=](double x){return 2./(gamma+1)* (cs_l+ (x-x_center)/t);};
    auto P2   = [=](double x){return Pl * std::pow(1- 0.5* (gamma-1)* v2(x)/cs_l,1./Beta);};
    auto rho2 = [=](double x){return rhol* std::pow(1- 0.5*(gamma-1) * v2(x)/cs_l,2./(gamma-1));};


    //calculate shock velocities

    double c12 = vl-cs_l;
    double c23 = v3-c3;
    double c34 = (rho4*v4-rho3*v3)/(rho4-rho3);
    double c45 = vr+cs_r*std::sqrt( 1+(gamma+1)/(2*gamma)*(P4/Pr-1) );//cs_r+(rhor*vr - rho4*v4)/(rhor-rho4);


    //contact positions

    double x_12 = x_center + c12*t;
    double x_23 = x_center + c23 *t;
    double x_34 = x_center + c34*t;
    double x_45 = x_center + c45*t;
    std::cout << "CS=" << c12 << " " << c23 << " " << c34 << " " << c45 << std::endl;
    std::cout << "Borders are:" << x_12 << " " << x_23 << " " << x_34 << " " << x_45 << std::endl; 
    std::vector<std::array<double,3>> output{X_arr.size(),{0,0,0}};

    std::transform(X_arr.begin(),X_arr.end(),output.begin(),[&](auto x)
    {
        double rho,v,p;
        rho = rhol;
        p = Pl;
        v = vl;
        if (x < x_45)
        {
            if(x > x_12)
            {
                if (x>x_23)
                {
                    if (x >= x_34)
                    {
                        rho = rho4;
                        v = v4;
                        p = P4;
                    }
                    else
                    {
                        rho = rho3;
                        p = P3;
                        v = v3;
                    }
                }
                else
                {
                    rho = rho2(x);
                    p = P2(x);
                    v = v2(x);
                }
            }
        }
        else
        {
            p = Pr;
            rho = rhor;
            v = vr;
        }
        return std::array<double,3> {rho,v,p};
    });

    return output;

}