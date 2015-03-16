#ifndef GLOBAL_H
#define GLOBAL_H

#include <iostream>
#include <fstream>
#include <vector>
#include <complex>

#include <btas/common/blas_cxx_interface.h>
#include <btas/common/TVector.h>
#include <btas/DENSE/TArray.h>

using std::ostream;
using std::vector;
using std::complex;

using namespace btas;

class Random;
class SpinQuantum;

namespace global {

   //!random number generator
   extern Random RN;

   //!lenght of the lattice
   extern int L;

   //!bond dimension of the tensors
   extern int D;

   //!physical dimension
   extern int d;

   //physical indices
   extern Qshapes<SpinQuantum> qp;

   extern SpinQuantum qt;

   //!initializer
   void init(const SpinQuantum &,int,int,int);

   template<typename T>
      T rgen();

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
