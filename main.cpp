/**
 * @mainpage 
 * This is a simple dmrg algorithm based on Naoki Nakatini's example dmrg code from the BTAS library.
 * I added my own MPO's and MPS/MPO functionalities from MPSblas.h, to compress, canonicalize and contracts MPS and MPO's.
 * It uses the two and one-site algorithm to optimize general spin systems. There are some predefined MPO's in one and two dimensions, 
 * and the possibility to construct a general Hamiltonian: H = \sum_ij J_ij S_i . S_j
 * @author Brecht Verstichel
 * @date 11-03-2015
 */
#include <iostream>
#include <iomanip>
using namespace std;

#include "include.h"

int main(int argc, char* argv[]){

   cout.precision(16);

   int L =  10;
   int D = 256;

   double J2 = 0.0;

   int d = 2;

   global::init(Quantum(2),D,d,L);

   //set the coupling matrix:
   DArray<2> J(L,L);
   coupling::J1J2_1D(false,J2,J);

   //set MPO to the Heisenberg model
   mpsxx::MPO<Quantum> mpo = SpinHamiltonian::heisenberg(J,0.0);

   //and canonicalize it
   compress(true,mpo,mpsxx::Left,0);
   compress(true,mpo,mpsxx::Right,0);
   compress(true,mpo,mpsxx::Left,0);
   compress(true,mpo,mpsxx::Right,0);

   //initialize the mps structure
   mpsxx::MPS<Quantum> mps = mpsxx::create<Quantum>(L,global::qt,global::qp,1,global::rgen<double>); 

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

   std::vector<double> spin_exp = SpinHamiltonian::local_spin(mps);

   cout << endl;
   cout << "Expectation value of local Sz on site i:" << endl;
   cout << endl;

   double total_spin = 0.0;

   for(int i = 0;i < L;++i){

      cout << i << "\t" << spin_exp[i] << endl;

      total_spin += spin_exp[i];

   }

   cout << endl;
   cout << "Expectation value of total Sz:\t" << total_spin << endl;

   return 0;

}
