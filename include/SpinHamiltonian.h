#ifndef SPINHAMILTONIAN_H
#define SPINHAMILTONIAN_H

#include <iostream>
#include <iomanip>

#include "MPSblas.h"

using namespace btas;
using namespace mpsxx;

class SpinQuantum;

namespace SpinHamiltonian {

   //some functions which initialize an MPO to a certain Hamiltonian
   MPO<SpinQuantum> ising(double,double);

   MPO<SpinQuantum> heisenberg(double,double,double);

   MPO<SpinQuantum> heisenberg(const DArray<2> &,double);

   std::vector<double> local_spin(const MPS<SpinQuantum> &);

   //general functions to make MPO construction easier
   void insert_id(QSDArray<4> &O,int row,int column,double value);
   void insert_Sz(QSDArray<4> &O,int row,int column,double value);
   void insert_Sp(QSDArray<4> &O,int row,int column,double value);
   void insert_Sm(QSDArray<4> &O,int row,int column,double value);

}

#endif
