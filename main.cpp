#include <iostream>
#include <iomanip>
using namespace std;

#include "include.h"

/**
 * simple dmrg algorithm based on Naoki Nakatini's example dmrg code from the BTAS library.
 * I added my own MPO's and MPS/MPO functionalities from MPSblas.h
 **/
int main(int argc, char* argv[]){

   cout.precision(16);

   int L =  1000;
   int D = 81;

   int d = 3;

   global::init(D,d,L);

   //set MPO to the Heisenberg model
   mpsxx::MPO<Quantum> mpo = SpinHamiltonian::heisenberg(1.0,1.0,0.0);

   //and canonicalize it
   compress(true,mpo,mpsxx::Left,0);
   compress(true,mpo,mpsxx::Right,0);
   compress(true,mpo,mpsxx::Left,0);
   compress(true,mpo,mpsxx::Right,0);

   //initialize the mps structure
   mpsxx::MPS<Quantum> mps = mpsxx::create<Quantum>(L,Quantum(0),global::qp,1,global::rgen<double>); 

   //and canonicalize it
   compress(true,mps,mpsxx::Left,0);
   compress(true,mps,mpsxx::Right,D);

   normalize(mps);

   compress(true,mps,mpsxx::Left,0);
   compress(false,mps,mpsxx::Right,0);

   //two-site
   double energy = algorithm::dmrg(mpo,mps, algorithm::TWOSITE);

   cout.precision(16);
   cout << "\tGround state energy (two-site) = " << setw(20) << fixed << energy << endl << endl;

   //one-site
   energy = dmrg(mpo,mps, algorithm::ONESITE);

   cout.precision(16);
   cout << "\tGround state energy (one-site) = " << setw(20) << fixed << energy << endl << endl;

   return 0;

}
