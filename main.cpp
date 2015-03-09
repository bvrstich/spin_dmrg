#include <iostream>
#include <iomanip>
using namespace std;

#include "include.h"

int main(int argc, char* argv[]){

   cout.precision(16);

  //
  // define working space for 20 sites chain
  //

   int L =  100;
   int D = 200;

   int d = 2;

   global::init(D,d,L);

   //set MPO to the Heisenberg model
   mpsxx::MPO<Quantum> mpo = SpinHamiltonian::heisenberg(L,d,1.0,1.0,0.0);

   //initialize the mps structure
   mpsxx::MPS<Quantum> mps = mpsxx::create<Quantum>(L,Quantum(0),global::qp,10,global::rgen<double>); 

   //and canonicalize it
   compress(mps,mpsxx::Left,0);
   compress(mps,mpsxx::Right,D);

   double energy = algorithm::dmrg(mpo,mps, algorithm::TWOSITE);
/*
   cout.precision(16);
   cout << "\tGround state energy (two-site) = " << setw(20) << fixed << energy << endl << endl;

   cout << "\tCalling DMRG program ( one-site algorithm) " << endl;

   energy = dmrg(sites, ONESITE, M);
   cout.precision(16);
   cout << "\tGround state energy (one-site) = " << setw(20) << fixed << energy << endl << endl;
*/
   return 0;

}
