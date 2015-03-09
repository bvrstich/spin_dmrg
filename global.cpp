#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <complex>

using std::cout;
using std::endl;
using std::vector;
using std::complex;
using std::ofstream;

#include "include.h"

namespace global{

   int D;

   int L;

   int d;

   Random RN;

   Qshapes<SpinQuantum> qp;

   /**
    * @param D_in virtual dimension of the trial
    * @param d_in physical dimension
    * @param L_in input length of the chain
    */
   void init(int D_in,int d_in,int L_in){

      L = L_in;
      d = d_in;
      D = D_in;

      //set physical quantum numbers for spinquantum
      qp.clear();

      int m = -d + 1;

      while(m < d){

         qp.push_back(SpinQuantum(m));

         m += 2;

      }

   }

   //!function which generates random complex numbers uniformly on a square of side 2 [(-1,1):(-1,1)]
   template<>
      complex<double> rgen(){ 

         return complex<double>(2.0*RN() - 1.0,2.0*RN() - 1.0); 

      }

   //!function which generates uniform random numbers between [-1:1]
   template<>
      double rgen(){ 

         return 2.0*RN() - 1.0;

      }

}
