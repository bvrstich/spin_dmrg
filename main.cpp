#include <iostream>
#include <iomanip>
using namespace std;

#include "include.h"

int main(int argc, char* argv[]){

   cout.precision(16);

  //
  // define working space for 20 sites chain
  //

   int L =  10;
   int D = 10;

   int d = 2;

   global::init(D,d,L);

   //set MPO to the Heisenberg model
   mpsxx::MPO<Quantum> mpo = SpinHamiltonian::heisenberg(1.0,1.0,0.0);

   //and canonicalize it
   compress(true,mpo,mpsxx::Left,0);
   compress(true,mpo,mpsxx::Right,0);
   compress(true,mpo,mpsxx::Left,0);
   compress(true,mpo,mpsxx::Right,0);

   cout << mpo << endl;

   //initialize the mps structure
   mpsxx::MPS<Quantum> mps = mpsxx::create<Quantum>(L,Quantum(0),global::qp,10,global::rgen<double>); 

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
