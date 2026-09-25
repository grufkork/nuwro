#include "generatormt.h"
#include "jednostki.h"
#include "kinsolver.h"
#include <cassert>
#include "particle.h"
#include "qel_sigma.h"
#include <iostream>
#include <fstream>
#include <iterator>
#include <sstream>
#include <cmath>
#include "params.h"
#include "event1.h"
#include "kinematics.h"
#include "pdg.h"
#include "nucleus.h"
#include "nucleus_data.h"
#include "dis/LeptonMass.h"
#include <cstdlib>

#define LOCALKF localkf_O

#include "rpa_2013.h"
#include "hyperon_interaction.h"

#include "ab_initio_event.h"
#include "ab_initio_data/lfg_R00.h"
#include "ab_initio_data/lfg_R0z.h"
#include "ab_initio_data/lfg_Rxx.h"
#include "ab_initio_data/lfg_Rxy.h"
#include "ab_initio_data/lfg_Rzz.h"


int Stupid_2D_to_1D(int n, int m, int N){
    return n*N + m;
}

////////////////////////////////////////////////////////////////////////
double ab_initio_event(params &p, event &e, nucleus &t, bool nc)
    ////////////////////////////////////////////////////////////////////////
{


    // Print log
    // std::cout << "QEL event: " << (nc ? "NC" : "CC") << std::endl;
    // std::cout << "  Neutrino: " << e.in[0] << std::endl;
    // std::cout << "VALLS" << std::endl;
    // std::cout << "  Nucleon: " << e.in[1] << std::endl;
    // std::exit(0);
    e.weight=0;

    particle lepton_in=e.in[0];    // Incoming neutrino (or electron)
    particle nucleon_in=e.in[1];    // initial nucleon
    particle lepton_out;        // Outgoing lepton
    particle nucleon_out;            // Outgoing nucleon

    nucleon_out.r=nucleon_in.r; // 4-position? at time of event
    lepton_out.r=nucleon_in.r;

    int kind=0; // 0 - cc //  1 - nc proton // 2 - nc neutron
    // std::cout << " -- " << std::endl;
    if(nc)
    {
        lepton_out=lepton_in;
        nucleon_out=nucleon_in;
        kind=(nucleon_in.pdg==pdg_proton?1:2);
    }
    else if((lepton_in.pdg>0 && nucleon_in.pdg==PDG::pdg_proton) ||( lepton_in.pdg<0 && nucleon_in.pdg==PDG::pdg_neutron))
    {
        // std::cout << "Invalid CC event: " << lepton_in << " " << nucleon_in << std::endl;
        return 0;
    }
    else
    {
        lepton_out.pdg=lepton_in.pdg-(lepton_in.pdg>0 ? 1 :-1);
        lepton_out.set_mass(PDG::mass(lepton_out.pdg));
        nucleon_out.pdg=(lepton_in.pdg>0 ? PDG::pdg_proton : PDG::pdg_neutron);//zmiana JN
        switch(p.qel_rpa)
        {
            case 2:
            case 3:
                nucleon_in.set_mass(t.Mf());
                nucleon_out.set_mass(t.Mf());
                break;
            default:
                nucleon_out.set_mass(PDG::mass(nucleon_out.pdg));
                break;
        }
    }

    if (lepton_in.pdg==11 && nucleon_in.pdg==pdg_proton)
        kind=10;
    if (lepton_in.pdg==11 && nucleon_in.pdg==pdg_neutron)
        kind=11;//will be used when FF are selected in the file ff.cc


    vect aa;
    aa = vect (nucleon_in);
    double dlu=sqrt( (aa + vect(lepton_in))*(aa+vect(lepton_in)) );
    if (dlu<= (lepton_out.mass() + nucleon_out.mass() ) )
    {
        // std::cout << "Event kinematically forbidden: " << lepton_in << " " << nucleon_in << std::endl;
        e.weight=0;
        return 0;
    }

    double m_l = lepton_out.mass();

    double coswidth = 2.0;
    // double cosval = coswidth * frandom();
    double costheta = coswidth * frandom() - 1.0;


    double eps = lepton_in.t;
    double k = eps; // Because neutrino mass is approx zero

    double eps_prim_min = m_l;
    double eps_prim_max = eps;

    double w_min = eps-eps_prim_max;
    double w_max = eps-eps_prim_min;

    // double wwidth = 1200.0;
    double wwidth = w_max - w_min;

    double w = wwidth * frandom() + w_min;
    double eps_prim = eps - (w);
    double k_prim = sqrt(eps_prim*eps_prim- m_l*m_l);
    // double k_prim = sqrt(eps_prim*eps_prim - m_l*m_l);
    double q = sqrt(-(costheta * 2.0 * k * k_prim - k*k - k_prim*k_prim));

    double xfade_point = 400.0;
    double xfade_width = 0.0;
    double crossfade = 1.0 - max(0.0, min(1.0, (q - xfade_point) / xfade_width));
    // crossfade = 0.0;

    if (crossfade == 0.0){
        e.weight = 0.0;
        return 0.0;
    }

    if ( w > q) {
        std::cout << "draw fail" << std::endl;
        e.weight = 0.0;
        return 0.0;
    }

    double xsec = 0;
    double omega = w;
    bool is_anti = lepton_in.pdg < 0;

    xsec = calc_xsec(costheta, q, omega, k, k_prim, eps, eps_prim, m_l, is_anti) * wwidth * coswidth;// * 2.0 * Pi;// * sin(acos(costheta));


    // std::cout << q_final << std::endl;

    double momentum_perpendicular = sqrt(1-costheta*costheta) * k_prim;
    double azimuth = 2.0 * M_PI * frandom();

    lepton_out.t = eps_prim;
    lepton_out.x = momentum_perpendicular * cos(azimuth);
    lepton_out.y = momentum_perpendicular * sin(azimuth);
    lepton_out.z = k_prim * costheta;

    // xsec *= 10e-3; // To GeV?

    if (xsec < 0.0){
        e.weight = 0.0;
        return 0.0;
    }


    e.temp.push_back(lepton_out);
    e.temp.push_back(nucleon_out);
    e.out.push_back(lepton_out);
    e.out.push_back(nucleon_out);
    e.weight=xsec/cm2 * crossfade;
    // std::cout << "accepted" << std::endl;



    return e.weight; // *cm2;
}

double calc_xsec(double costheta, double q, double omega, double k, double k_prim, double eps, double eps_prim, double m_l, bool is_anti){
    double spacing = 4.0;
    int gridsteps = 300;
    double PI = 3.14159265359;





    // double Q_sq = -q2;
    double Q_sq = q*q - omega*omega;

    // std::cout << "q: " << q << " omega: " << omega << " k: " << k << " k_prim: " << k_prim << " eps: " << eps << " eps': " << eps_prim << std::endl;


    // q = k - k'
    // q^2 = kk - 2kk' + k'k'
    // kk' = |k||k'|cos(theta)
    // cos(theta) = (kk + k'k' - q^2) / (2|k||k'|)
    // double temp = (k*k + k_prim * k_prim - q*q) / (2.0 * k * k_prim);
    // -> costheta * 2 * k * k_prim - k * k - k_prim * k_prim = q^2
    // double costheta = temp;




    //  w/y
    //  ^  2  3 
    //  |  
    //  |  0  1
    //   -----> q/x

    int table_x, table_y;
    double xfrac, yfrac;

    table_x = (int)floor(q/spacing);
    table_y = (int)floor(omega/spacing);
    // std::cout << table_x << " " << table_y << std::endl;
    // std::cout << q << " " << omega << std::endl;
    xfrac = q/spacing - table_x;
    yfrac = omega /spacing - table_y;

    // List of pointers to curves
    double *curves[] = {lfg_R00, lfg_R0z, lfg_Rzz, lfg_Rxx, lfg_Rxy};
    double resolved_vals[] = {0.0, 0.0, 0.0, 0.0, 0.0};

    for(int i = 0; i < 5; i++){
        double val__, valx_, val_y, valxy;

        val__ = curves[i][Stupid_2D_to_1D(table_x + 0, table_y + 0, gridsteps)];
        valx_ = curves[i][Stupid_2D_to_1D(table_x + 1, table_y + 0, gridsteps)];
        val_y = curves[i][Stupid_2D_to_1D(table_x + 0, table_y + 1, gridsteps)];
        valxy = curves[i][Stupid_2D_to_1D(table_x + 1, table_y + 1, gridsteps)];

        double valA, valB;

        valA = (1.0 - xfrac) * val__ + xfrac * valx_;
        valB = (1.0 - xfrac) * val_y + xfrac * valxy;

        double final_val = (1.0 - yfrac) * valA + yfrac * valB;
        resolved_vals[i] = final_val;
    }

    double v00, v0z, vzz, vxx, vxy;
    v00 = 2.0 * eps * eps_prim * ( 1 + k_prim/eps_prim * costheta );
    v0z = omega / q * (m_l*m_l + v00) + m_l*m_l / q * (eps + eps_prim);
    vzz = omega * omega / (q*q) * (m_l * m_l + v00) + m_l*m_l / (q * q) * (m_l * m_l + 2.0 * omega * (eps + eps_prim) + q*q);
    vxx = Q_sq + Q_sq * (m_l * m_l + v00) / (2.0 * q * q) - m_l*m_l / (q * q) * (m_l*m_l / 2.0 + omega * (eps + eps_prim));
    vxy = Q_sq * (eps + eps_prim) / q - m_l * m_l * omega / q;

    double vxy_sign = -1.0; // If neutrino, I think it should be positive for antineutrinos
    if (is_anti){
        vxy_sign *= -1.0;
    }

    double coefficients[] = {v00, -v0z, vzz, vxx, vxy * vxy_sign};

    double xsec = 0.0;
    for(int i = 0; i < 5; i++){
        xsec += coefficients[i] * resolved_vals[i];
        // xsec += resolved_vals[i];
    }
    // std::cout << "  " << xsec << std::endl;

    // std::cout << q << " " << omega << " " << xsec << std::endl;

    xsec *= G * G / (8. * PI * PI) * k_prim / eps;

    return xsec;
}
