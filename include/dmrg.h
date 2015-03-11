#ifndef _PROTOTYPE_DMRG_H
#define _PROTOTYPE_DMRG_H 1

#include "MPSblas.h"

namespace algorithm {

   enum DMRG_ALGORITHM { ONESITE, TWOSITE };

   template<class Q>
      double dmrg_sweep(const mpsxx::MPO<Q> &,mpsxx::MPS<Q> &,std::vector< QSDArray<3> > &LO,std::vector< QSDArray<3> > &RO,DMRG_ALGORITHM algo);

   template<class Q>
      double dmrg(const mpsxx::MPO<Q> &,mpsxx::MPS<Q> &, DMRG_ALGORITHM algo);

   template<class Q> 
      void init_ro(const mpsxx::MPO<Q> &,const mpsxx::MPS<Q> &,std::vector< QSDArray<3> > &LO,std::vector< QSDArray<3> > &RO);
  
   template<class Q>
      double optimize_onesite(bool forward,int site, const mpsxx::MPO<Q> &mpo, mpsxx::MPS<Q> &mps,
            
            std::vector< QSDArray<3> > &LO, std::vector< QSDArray<3> > &RO );

   template<class Q>
      double optimize_twosite(bool forward,int site, const mpsxx::MPO<Q> &mpo, mpsxx::MPS<Q> &mps,
            
            std::vector< QSDArray<3> > &LO, std::vector< QSDArray<3> > &RO );


};

#endif // _PROTOTYPE_DMRG_H
