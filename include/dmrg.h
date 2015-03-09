#ifndef _PROTOTYPE_DMRG_H
#define _PROTOTYPE_DMRG_H 1

#include "MPSblas.h"

namespace algorithm {

   enum DMRG_ALGORITHM { ONESITE, TWOSITE };
/*
   double optimize_onesite(bool forward, MpSite& sysdot, MpSite& envdot, int M = 0);
   double optimize_twosite(bool forward, MpSite& sysdot, MpSite& envdot, int M = 0);
 */
   template<class Q>
      double dmrg_sweep(const mpsxx::MPO<Q> &,mpsxx::MPS<Q> &, DMRG_ALGORITHM algo);

   template<class Q>
      double dmrg(const mpsxx::MPO<Q> &,mpsxx::MPS<Q> &, DMRG_ALGORITHM algo);

   template<class Q> 
      void init_ro(const mpsxx::MPO<Q> &,const mpsxx::MPS<Q> &,std::vector< QSDArray<3> > &RO,std::vector< QSDArray<3> > &LO);

};

#endif // _PROTOTYPE_DMRG_H
