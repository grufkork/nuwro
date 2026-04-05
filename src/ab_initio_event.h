#include <cassert>          
#include "params.h"
#include "particle.h"
#include "event1.h"
#include "nucleus.h"

double ab_initio_event(params &p, event &e, nucleus &t, bool nc);
double calc_xsec(double q, double omega, double eps, double m_l, bool is_anti);
int Stupid_2D_to_1D(int n, int m, int N);
