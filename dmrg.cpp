#include <iostream>
#include <iomanip>
#include <vector>
using namespace std;

#include "include.h"

using namespace btas;

namespace algorithm {

   /**
    * imlementation of the single-site algorithm for dmrg update
    * @param forward or backward sweep
    * @param site , location of the current tensor
    * @param mpo containing the hamiltonian
    * @param mps current state 
    * @param LO left renormalized operator
    * @param RO right renormalized operator
    */
   template<class Q>
      double optimize_onesite(bool forward,int site, const mpsxx::MPO<Q> &mpo, mpsxx::MPS<Q> &mps,
            
            std::vector< QSDArray<3> > &LO, std::vector< QSDArray<3> > &RO ) {

      boost::function<void(const QSDArray<3>&, QSDArray<3>&)>
         f_contract = boost::bind(ComputeSigmaVector, mpo[site], LO[site], RO[site], _1, _2);

      QSDArray<3> diag(mps[site].q(), mps[site].qshape());

      ComputeDiagonal(mpo[site], LO[site], RO[site], diag);

      double energy = davidson::diagonalize(f_contract, diag, mps[site]);

      if(forward) {

         QSDArray<3> U;
         Canonicalize(true, mps[site] , U, global::D);

         QSDArray<3> tmp3;
         ComputeGuess(true, U, mps[site], mps[site+1], tmp3);

         //move the updated tensors
         mps[site] = std::move(U);
         mps[site+1] = std::move(tmp3);

         //make new operator
         Renormalize (true, mpo[site],LO[site], mps[site], mps[site], LO[site+1]);

      }
      else {

         QSDArray<3> VT;
         Canonicalize(false, mps[site], VT, global::D);

         QSDArray<3> tmp3;
         ComputeGuess(false, VT, mps[site], mps[site-1], tmp3);

         //move the updated tensors
         mps[site] = std::move(VT);
         mps[site-1] = std::move(tmp3);

         Renormalize (false, mpo[site], RO[site], mps[site], mps[site], RO[site-1]);

      }

      return energy;

   }

/*
   double optimize_twosite(bool forward, MpSite& sysdot, MpSite& envdot, int M) {

      QSDArray<4> wfnc;
      QSDArray<4> diag;
      boost::function<void(const QSDArray<4>&, QSDArray<4>&)> f_contract;
      if(forward) {
         QSDgemm(NoTrans, NoTrans, 1.0, sysdot.wfnc, envdot.rmps, 1.0, wfnc);
         f_contract = boost::bind(ComputeSigmaVector, sysdot.mpo, envdot.mpo, sysdot.lopr, envdot.ropr, _1, _2);
         diag.resize(wfnc.q(), wfnc.qshape());
         ComputeDiagonal(sysdot.mpo, envdot.mpo, sysdot.lopr, envdot.ropr, diag);
      }
      else {
         QSDgemm(NoTrans, NoTrans, 1.0, envdot.lmps, sysdot.wfnc, 1.0, wfnc);
         f_contract = boost::bind(ComputeSigmaVector, envdot.mpo, sysdot.mpo, envdot.lopr, sysdot.ropr, _1, _2);
         diag.resize(wfnc.q(), wfnc.qshape());
         ComputeDiagonal(envdot.mpo, sysdot.mpo, envdot.lopr, sysdot.ropr, diag);
      }

      double energy = davidson::diagonalize(f_contract, diag, wfnc);

      if(forward) {
         Canonicalize(1,        wfnc, sysdot.lmps, envdot.wfnc, M);
         envdot.lopr.clear();
         Renormalize (1, sysdot.mpo,  sysdot.lopr, sysdot.lmps, sysdot.lmps, envdot.lopr);
      }
      else {
         Canonicalize(0,        wfnc, sysdot.rmps, envdot.wfnc, M);
         envdot.ropr.clear();
         Renormalize (0, sysdot.mpo,  sysdot.ropr, sysdot.rmps, sysdot.rmps, envdot.ropr);
      }

      return energy;
   }

  */
   /**
    * one sweep of the DMRG alorithm,
    * @param mpo object containing the interactions
    * @param mps initial MPS object
    * @param algo flag for two or one site update algorithm
    */
   template<class Q>
      double dmrg_sweep(const mpsxx::MPO<Q> &mpo,mpsxx::MPS<Q> &mps,std::vector< QSDArray<3> > &LO,
            
            std::vector< QSDArray<3> > &RO, DMRG_ALGORITHM algo){

         double emin = 1.0e8;
            
         // forward sweep
         cout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
         cout << "\t\t\tFORWARD SWEEP" << endl;
         cout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

         //for(int i = 0; i < global::L-1; ++i){
int i = 0;
            // diagonalize
            double eswp;

          //  if(algo == ONESITE)
               eswp = optimize_onesite(true,i, mpo, mps,LO,RO );
               /*
           // else
               eswp = optimize_twosite(true, sites[i], sites[i+1], M);

            if(eswp < emin)
               emin = eswp;

            // print result
            cout.precision(16);
            cout << "\t\t\tEnergy = " << setw(24) << fixed << eswp << endl;

         }

         // backward sweep
         cout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
         cout << "\t\t\tBACKWARD SWEEP" << endl;
         cout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

         for(int i = L-1; i > 0; --i) {

            // diagonalize
            double eswp;

            if(algo == ONESITE)
               eswp = optimize_onesite(false, sites[i], sites[i-1], M);
            else
               eswp = optimize_twosite(false, sites[i], sites[i-1], M);

            if(eswp < emin)
               emin = eswp;

            // print result
            cout.precision(16);
            cout << "\t\t\tEnergy = " << setw(24) << fixed << eswp << endl;

         }
*/
         return emin;

      }

   /**
    * implementation of sweeping algorithm for dmrg optimization
    * @param mpo object containing the interactions
    * @param mps initial MPS object
    * @param algo flag for two or one site update algorithm
    */
   template<class Q>
      double dmrg(const mpsxx::MPO<Q> &mpo,mpsxx::MPS<Q> &mps, DMRG_ALGORITHM algo){

         double esav = 1.0e8;

         //left and right renormalized operators
         std::vector< QSDArray<3> > LO(global::L);
         std::vector< QSDArray<3> > RO(global::L);

         //initialize the right renormalized operator
         algorithm::init_ro(mpo,mps,LO,RO);

         for(int iter = 0; iter < 100; ++iter){

            cout << "\t====================================================================================================" << endl;
            cout << "\t\tSWEEP ITERATION [ " << setw(4) << iter << " ] "   << endl;
            cout << "\t====================================================================================================" << endl;

            double eswp = dmrg_sweep(mpo,mps,LO,RO,algo);

            double edif = eswp - esav;

            cout << "\t====================================================================================================" << endl;
            cout << "\t\tSWEEP ITERATION [ " << setw(4) << iter << " ] FINISHED" << endl;

            cout.precision(16);

            cout << "\t\t\tSweep Energy = " << setw(24) << fixed << eswp << " ( delta E = ";
            cout.precision(2);
            cout << setw(8) << scientific << edif << " ) " << endl;
            cout << "\t====================================================================================================" << endl;
            cout << endl;

            esav = eswp;

            if(fabs(edif) < 1.0e-8)
               break;
         }

         return esav;

      }

   /**
    * initialize the right renormalized operator using a right canonical input mps
    * @param mpo object containing the interactions
    * @param mps initial MPS object
    * @param RO on output the right renormalized operator
    */
   template<class Q> 
      void init_ro(const mpsxx::MPO<Q> &mpo,const mpsxx::MPS<Q> &mps,std::vector< QSDArray<3> > &LO,std::vector< QSDArray<3> > &RO){

         Qshapes<Q> qz(1, Q::zero());
         Dshapes dz(qz.size(), 1);

         TVector<Qshapes<Quantum>, 3> qshape = make_array( qz, qz, qz);
         TVector<Dshapes,          3> dshape = make_array( dz, dz, dz);

         RO[global::L-1].resize(Q::zero(), qshape, dshape);
         RO[global::L-1] = 1.0;

         for(int i = global::L-1; i > 0; --i) {

            RO[i-1].clear();
            Renormalize (0, mpo[i], RO[i], mps[i], mps[i], RO[i-1]);

         }

         LO[0].resize(Q::zero(), qshape, dshape);
         LO[0] = 1.0;

      }

   //forward declarations
   template double dmrg<Quantum>(const mpsxx::MPO<Quantum> &mpo,mpsxx::MPS<Quantum> &mps, DMRG_ALGORITHM algo);

   template double dmrg_sweep<Quantum>(const mpsxx::MPO<Quantum> &mpo,mpsxx::MPS<Quantum> &mps,
         
         std::vector< QSDArray<3> > &LO,std::vector< QSDArray<3> > &RO, DMRG_ALGORITHM algo);

   template void init_ro<Quantum>(const mpsxx::MPO<Quantum> &,const mpsxx::MPS<Quantum> &,std::vector< QSDArray<3> > &LO,std::vector< QSDArray<3> > &RO);

   template double optimize_onesite<Quantum>(bool forward,int site, const mpsxx::MPO<Quantum> &, mpsxx::MPS<Quantum> &,
         
         std::vector< QSDArray<3> > &, std::vector< QSDArray<3> > & );


}
