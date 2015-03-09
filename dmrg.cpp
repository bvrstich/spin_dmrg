#include <iostream>
#include <iomanip>
#include <vector>
using namespace std;

#include "include.h"

using namespace btas;

double prototype::optimize_onesite(bool forward, MpSite& sysdot, MpSite& envdot, int M)
{
   boost::function<void(const QSDArray<3>&, QSDArray<3>&)>
      f_contract = boost::bind(ComputeSigmaVector, sysdot.mpo, sysdot.lopr, sysdot.ropr, _1, _2);
   QSDArray<3> diag(sysdot.wfnc.q(), sysdot.wfnc.qshape());
   ComputeDiagonal(sysdot.mpo, sysdot.lopr, sysdot.ropr, diag);
   double energy = davidson::diagonalize(f_contract, diag, sysdot.wfnc);

   if(forward) {
      Canonicalize(1, sysdot.wfnc, sysdot.lmps, M);
      ComputeGuess(1, sysdot.lmps, sysdot.wfnc, envdot.rmps, envdot.wfnc);
      envdot.lopr.clear();
      Renormalize (1, sysdot.mpo,  sysdot.lopr, sysdot.lmps, sysdot.lmps, envdot.lopr);
   }
   else {
      Canonicalize(0, sysdot.wfnc, sysdot.rmps, M);
      ComputeGuess(0, sysdot.rmps, sysdot.wfnc, envdot.lmps, envdot.wfnc);
      envdot.ropr.clear();
      Renormalize (0, sysdot.mpo,  sysdot.ropr, sysdot.rmps, sysdot.rmps, envdot.ropr);
   }

   return energy;
}

double prototype::optimize_twosite(bool forward, MpSite& sysdot, MpSite& envdot, int M)
{
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

double prototype::dmrg_sweep(MpStorages& sites, DMRG_ALGORITHM algo, int M)
{
   int    L    = sites.size();
   double emin = 1.0e8;
   // fowrad sweep
   cout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
   cout << "\t\t\tFORWARD SWEEP" << endl;
   cout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
   for(int i = 0; i < L-1; ++i) {
      // diagonalize
      double eswp;
      if(algo == ONESITE) eswp = optimize_onesite(1, sites[i], sites[i+1], M);
      else                eswp = optimize_twosite(1, sites[i], sites[i+1], M);
      if(eswp < emin) emin = eswp;
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
      if(algo == ONESITE) eswp = optimize_onesite(0, sites[i], sites[i-1], M);
      else                eswp = optimize_twosite(0, sites[i], sites[i-1], M);
      if(eswp < emin) emin = eswp;
      // print result
      cout.precision(16);
      cout << "\t\t\tEnergy = " << setw(24) << fixed << eswp << endl;
   }
   return emin;
}

double prototype::dmrg(MpStorages& sites, DMRG_ALGORITHM algo, int M)
{

   double esav = 1.0e8;
   for(int iter = 0; iter < 100; ++iter) {
      cout << "\t====================================================================================================" << endl;
      cout << "\t\tSWEEP ITERATION [ " << setw(4) << iter << " ] "   << endl;
      cout << "\t====================================================================================================" << endl;
      double eswp = dmrg_sweep(sites, algo, M);
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
      if(fabs(edif) < 1.0e-8) break;
   }

   return esav;
}
