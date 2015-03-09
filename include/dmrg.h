#ifndef _PROTOTYPE_DMRG_H
#define _PROTOTYPE_DMRG_H 1

#include "MPSblas.h"

namespace algorithm {

   enum DMRG_ALGORITHM { ONESITE, TWOSITE };
/*
   double optimize_onesite(bool forward, MpSite& sysdot, MpSite& envdot, int M = 0);
   double optimize_twosite(bool forward, MpSite& sysdot, MpSite& envdot, int M = 0);

   double dmrg_sweep(MpStorages& sites, DMRG_ALGORITHM algo, int M = 0);
  */
   template<class Q>
      double dmrg(mpsxx::MPO<Q> &,mpsxx::MPS<Q> &, DMRG_ALGORITHM algo);

};

#endif // _PROTOTYPE_DMRG_H
