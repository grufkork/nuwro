#include <cassert>          
#include "params.h"
#include "particle.h"
#include "event1.h"
#include "nucleus.h"

double ab_initio_event(params &p, event &e, nucleus &t, bool nc);
double calc_xsec(double costheta, double q, double omega, double k, double k_prim, double eps, double eps_prim, double m_l, bool is_anti, double _E_bind);
int Stupid_2D_to_1D(int n, int m, int N);
