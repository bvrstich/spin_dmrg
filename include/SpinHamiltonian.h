#ifndef SPINHAMILTONIAN_H
#define SPINHAMILTONIAN_H

#include <iostream>
#include <iomanip>

#include "MPSblas.h"

using namespace btas;
using namespace mpsxx;

class SpinQuantum;

namespace SpinHamiltonian {

   //some functions which initialize an MPO to a certian Hamiltonian
   MPO<SpinQuantum> ising(int,int,double,double);

   MPO<SpinQuantum> XY(int,int,double,double);

   MPO<SpinQuantum> heisenberg(int,int,double,double,double);

   MPO<SpinQuantum> raise(int,int);

   MPO<SpinQuantum> lower(int,int);

   MPO<SpinQuantum> Sz(int,int);

}

#endif
