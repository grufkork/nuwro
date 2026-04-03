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
#include "lfg_data/lfg_R00.h"
#include "lfg_data/lfg_R0z.h"
#include "lfg_data/lfg_Rxx.h"
#include "lfg_data/lfg_Rxy.h"
#include "lfg_data/lfg_Rzz.h"


int Stupid_2D_to_1D(int n, int m, int N){
    return n*N + m;
}

////////////////////////////////////////////////////////////////////////
double ab_initio_event(params &p, event &e, nucleus &t, bool nc)
    ////////////////////////////////////////////////////////////////////////
{
    double spacing = 10.0;
    int gridsteps = 120;
    double fermi_constant = 1.1663787e-11;
    double PI = 3.14159265359;


    // Print log
    // std::cout << "QEL event: " << (nc ? "NC" : "CC") << std::endl;
    // std::cout << "  Neutrino: " << e.in[0] << std::endl;
    // std::cout << "VALLS" << std::endl;
    // std::cout << "  Nucleon: " << e.in[1] << std::endl;
    // std::exit(0);
    e.weight=0;

    particle nu=e.in[0];    // Incoming neutrino (or electron)
    particle N0=e.in[1];    // initial nucleon
    particle lepton;        // Outgoing lepton
    particle N1;            // Outgoing nucleon

    N1.r=N0.r;
    lepton.r=N0.r;

    int kind=0; // 0 - cc //  1 - nc proton // 2 - nc neutron
    if(nc)
    {
        lepton=nu;
        N1=N0;
        kind=(N0.pdg==pdg_proton?1:2);
    }
    else if((nu.pdg>0 && N0.pdg==PDG::pdg_proton) ||( nu.pdg<0 && N0.pdg==PDG::pdg_neutron))
    {
        return 0;
    }
    else
    {
        lepton.pdg=nu.pdg-(nu.pdg>0 ? 1 :-1);
        lepton.set_mass(PDG::mass(lepton.pdg));
        N1.pdg=(nu.pdg>0 ? PDG::pdg_proton : PDG::pdg_neutron);//zmiana JN
        switch(p.qel_rpa)
        {
            case 2:
            case 3:
                N0.set_mass(t.Mf());
                N1.set_mass(t.Mf());
                break;
            default:
                N1.set_mass(PDG::mass(N1.pdg));
                break;
        }
    }

    if (nu.pdg==11 && N0.pdg==pdg_proton)
        kind=10;
    if (nu.pdg==11 && N0.pdg==pdg_neutron)
        kind=11;//will be used when FF are selected in the file ff.cc

    double _E_bind=0; //binding energy

    if(t.A() > 1)
    {
        switch(p.nucleus_target)
        {
            case 0: _E_bind = 0; break; // free nucleon
            case 1: _E_bind = p.nucleus_E_b; break; // GFG
            case 2: _E_bind = t.Ef(N0) + p.kaskada_w; break; // LFG
            case 3: _E_bind = bodek_binding_energy(N0, t.p, t.n); break; // Bodek-Ritchie
            case 4: _E_bind = binen (N0.p(), p.nucleus_p, p.nucleus_n); break; // effective SF
            case 5: _E_bind = deuter_binen (N0.p()); break; // deuterium
            case 6: _E_bind = 0; break; // effective potential
            default: _E_bind= 0;
        }
    }

    vect aa;
    aa = vect (N0);
    aa.t-=_E_bind;
    double dlu=sqrt( (aa + vect(nu))*(aa+vect(nu)) );
    if (dlu<= (lepton.mass() + N1.mass() ) )
    {
        e.weight=0;
        return 0;
    }

    // cross section (is 0 until the reaction occurs)
    double xsec = 0;
    double q2,jakobian;

    q2 = czarek_kinematics2(_E_bind, nu, N0, lepton, N1,jakobian);


    vect nu4 = nu;
    nu4.boost (-N0.v ());  // go to nucleon rest frame
    vect lepton4 = lepton;
    lepton4.boost(-N0.v());

    double Enu0 = nu4.t;   // neutrino energy in target frame

    double m_l = lepton.mass();

    // double omega = w;
    double eps = nu4.t;
    double eps_prim = lepton4.t;//eps - omega;
    double omega = eps-eps_prim;
    double Q_sq = -q2;
    double q = sqrt(Q_sq + omega * omega);

    double k = eps;
    double k_prim = sqrt(eps_prim * eps_prim - m_l*m_l);


    double temp = (k*k + k_prim * k_prim - q*q) / (2.0 * k * k_prim);
    double costheta = std::acos(temp);
    //  w/y
    //  ^  2  3 
    //  |  
    //  |  0  1
    //   -----> q/x

    // for(double iterval = 0.0; iterval < 300.0; iterval += 1.0){
        int table_x, table_y;
        double xfrac, yfrac;

        // q = 300;
        // omega = iterval;


        // std::cout << "q: " << q << " w: " << w << std::endl;

        table_x = (int)floor(q/spacing);
        table_y = (int)floor(omega/spacing);
        xfrac = q/spacing - table_x;
        yfrac = omega/spacing - table_y;

        // std::cout << "table_x: " << table_x << " table_y: " << table_y << std::endl;

        // List of pointers to curves
        double *curves[] = {lfg_R00, lfg_R0z, lfg_Rzz, lfg_Rxx, lfg_Rxy};
        double resolved_vals[] = {0.0, 0.0, 0.0, 0.0, 0.0};

        for(int i = 0; i < 5; i++){
            double val__, valx_, val_y, valxy;

            val__ = curves[i][Stupid_2D_to_1D(table_x + 0, table_y + 0, gridsteps)];
            valx_ = curves[i][Stupid_2D_to_1D(table_x + 1, table_y + 0, gridsteps)];
            val_y = curves[i][Stupid_2D_to_1D(table_x + 0, table_y + 1, gridsteps)];
            valxy = curves[i][Stupid_2D_to_1D(table_x + 1, table_y + 1, gridsteps)];

            // std::cout << val__ << " " << valx_ << " " << val_y << " " << valxy << std::endl;

            double valA, valB;

            valA = (1.0 - xfrac) * val__ + xfrac * valx_;
            valB = (1.0 - xfrac) * val_y + xfrac * valxy;

            double final_val = (1.0 - yfrac) * valA + yfrac * valB;
            resolved_vals[i] = final_val;
            // std::cout << final_val << std::endl;
        }

        // std::cout << q << "  " << omega << " - " << resolved_vals[0] << std::endl;
    // }



    double v00, v0z, vzz, vxx, vxy;
    v00 = 2.0 * eps * eps_prim * ( 1 + k_prim/eps_prim * costheta );
    v0z = omega / q * (m_l*m_l + v00) + m_l*m_l / q * (eps + eps_prim);
    vzz = omega * omega / (q*q) * (m_l * m_l + v00) + m_l*m_l / (q * q) * (m_l * m_l + 2.0 * omega * (eps + eps_prim) + q*q);
    vxx = Q_sq + Q_sq * (m_l * m_l + v00) / (2.0 * q * q) - m_l*m_l / (q * q) * (m_l*m_l / 2.0 + omega * (eps + eps_prim));
    vxy = Q_sq * (eps + eps_prim) / q - m_l * m_l * omega / q;

    double vxy_sign = -1.0; // If neutrino, I think it should be positive for antineutrinos
    if (nu.pdg < 0){
        vxy_sign = 1.0;
    }

    double coefficients[] = {v00, -v0z, vzz, vxx, vxy * vxy_sign};

    xsec = 0.0;
    for(int i = 0; i < 5; i++){
        // std::cout << resolved_vals[i] << " ";
        xsec += coefficients[i] * 1.0;//resolved_vals[i];
    }
    // std::cout << "  " << xsec << std::endl;

    xsec *= fermi_constant * fermi_constant / (8. * PI * PI) * k_prim / eps;


    // xsec = jakobian * qel_sigma(Enu0, q2, kind, nu.pdg<0, lepton.mass(), N0.mass());
    xsec *= jakobian;

    xsec *= 10e-3; // To MeV?

    // std::cout << "yippee" << std::endl;
    // std::exit(0);



    e.temp.push_back(lepton);
    e.temp.push_back(N1);
    e.out.push_back(lepton);
    e.out.push_back(N1);
    e.weight=xsec/cm2;

    return e.weight*cm2;
}
