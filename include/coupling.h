#ifndef COUPLING_H
#define COUPLING_H

#include <iostream>
#include <iomanip>
#include <complex>

using namespace btas;
using std::complex;

namespace coupling {

   void J1J2_1D(bool,double,TArray<double,2> &);

   void J1J2_2D(bool,double,TArray<double,2> &);

   void rand(TArray<double,2> &);

};

#endif
