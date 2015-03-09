#ifndef _PROTOTYPE_DMRG_H
#define _PROTOTYPE_DMRG_H 1

#include "mpsite.h"

#include "MPSblas.h"


namespace prototype {

   enum DMRG_ALGORITHM { ONESITE, TWOSITE };

   double optimize_onesite(bool forward, MpSite& sysdot, MpSite& envdot, int M = 0);
   double optimize_twosite(bool forward, MpSite& sysdot, MpSite& envdot, int M = 0);

   double dmrg_sweep(MpStorages& sites, DMRG_ALGORITHM algo, int M = 0);

   double dmrg(MpStorages& sites, DMRG_ALGORITHM algo, int M = 0);

};

#endif // _PROTOTYPE_DMRG_H
