#include <iostream>
#include <iomanip>

using std::cout;
using std::endl;
using std::ostream;

#include "include.h"

namespace SpinHamiltonian {


   /**
    * initialize the MPO<SpinQuantum> to represent a nearest-neighbour ising Hamiltonian on a lattice of size L and with coupling constant J
    * @param L length of the chain
    * @param d local dimension: i.e. defines the size of the local spins
    * @param J coupling constant
    * @param B magnetic fieldstrength in z direction
    */
   MPO<SpinQuantum> ising(int L,int d,double J,double B){

      double sz = 0.5 * (d - 1.0);//size of local spin

      MPO<SpinQuantum> mpo(L);

      Qshapes<SpinQuantum> qz;
      qz.push_back(SpinQuantum::zero());

      //incoming
      Qshapes<SpinQuantum> qi;
      qi.push_back(SpinQuantum::zero());//I has spin 0
      qi.push_back(SpinQuantum::zero());//Sz has spin 0
      qi.push_back(SpinQuantum::zero());//B has spin 0

      //outgoing
      Qshapes<SpinQuantum> qo;
      qo.push_back(SpinQuantum::zero());//I has spin 0
      qo.push_back(SpinQuantum::zero());//Sz has spin 0
      qo.push_back(SpinQuantum::zero());//B has spin 0

      TVector<Qshapes<SpinQuantum>,4> qshape = make_array(qz,global::qp,-global::qp,qo);

      //initialize the quantumnumbers of the MPO<SpinQuantum>
      mpo[0].resize(SpinQuantum::zero(),qshape);

      qshape = make_array(qi,global::qp,-global::qp,qo);

      for(int i = 1;i < L-1;++i)
         mpo[i].resize(SpinQuantum::zero(),qshape);

      qshape = make_array(qi,global::qp,-global::qp,-qz);

      mpo[L-1].resize(SpinQuantum::zero(),qshape);

      //first the identity
      double mz = -sz;

      for(int m = 0;m < d;++m){

         // set block elements
         DArray<4> I_op(1, 1, 1, 1);//identity
         I_op = 1.0;

         DArray<4> Sz_op(1, 1, 1, 1);//identity
         Sz_op = mz;

         DArray<4> B_op(1, 1, 1, 1);//identity
         B_op = -B*mz;

         DArray<4> J_op(1, 1, 1, 1);//identity
         J_op = J*mz;

         mpo[0].insert(shape(0,m,m,0),I_op);
         mpo[0].insert(shape(0,m,m,1),Sz_op);
         mpo[0].insert(shape(0,m,m,2),B_op);

         for(int i = 1;i < L - 1;++i){

            mpo[i].insert(shape(0,m,m,0),I_op);
            mpo[i].insert(shape(0,m,m,1),Sz_op);
            mpo[i].insert(shape(0,m,m,2),B_op);
            mpo[i].insert(shape(1,m,m,2),J_op);
            mpo[i].insert(shape(2,m,m,2),I_op);

         }

         mpo[L-1].insert(shape(0,m,m,0),B_op);
         mpo[L-1].insert(shape(1,m,m,0),J_op);
         mpo[L-1].insert(shape(2,m,m,0),I_op);

         mz += 1.0;

      }

      //merge everything together
      TVector<Qshapes<SpinQuantum>,1> qmerge;
      TVector<Dshapes,1> dmerge;

      qmerge[0] = mpo[0].qshape(3);
      dmerge[0] = mpo[0].dshape(3);

      QSTmergeInfo<1> info(qmerge,dmerge);

      QSDArray<4> tmp;
      QSTmerge(mpo[0],info,tmp);

      mpo[0] = tmp;

      for(int i = 1;i < L - 1;++i){

         //first merge the row
         qmerge[0] = mpo[i].qshape(0);
         dmerge[0] = mpo[i].dshape(0);

         info.reset(qmerge,dmerge);

         tmp.clear();

         QSTmerge(info,mpo[i],tmp);

         //then merge the column
         qmerge[0] = tmp.qshape(3);
         dmerge[0] = tmp.dshape(3);

         info.reset(qmerge,dmerge);

         mpo[i].clear();

         QSTmerge(tmp,info,mpo[i]);

      }

      //only merge row for i = L - 1
      qmerge[0] = mpo[L - 1].qshape(0);
      dmerge[0] = mpo[L - 1].dshape(0);

      info.reset(qmerge,dmerge);

      tmp.clear();

      QSTmerge(info,mpo[L - 1],tmp);

      mpo[L - 1] = tmp;


      return mpo;

   }

   /**
    * initialize the MPO<SpinQuantum> to represent a Sz operator
    * @param L length of the chain
    * @param d local dimension: i.e. defines the size of the local spins
    */
   MPO<SpinQuantum> Sz(int L,int d){

      double sz = 0.5 * (d - 1.0);//size of local spin

      MPO<SpinQuantum> mpo(L);

      Qshapes<SpinQuantum> qz;
      qz.push_back(SpinQuantum::zero());//Sz has spin 0

      //incoming
      Qshapes<SpinQuantum> qi;
      qi.push_back(SpinQuantum::zero());//Sz has spin 0
      qi.push_back(SpinQuantum::zero());//I has spin 0

      //outgoing
      Qshapes<SpinQuantum> qo;
      qo.push_back(SpinQuantum::zero());//Sz has spin 0
      qo.push_back(SpinQuantum::zero());//I has spin 0

      TVector<Qshapes<SpinQuantum>,4> qshape = make_array(qz,global::qp,-global::qp,qo);

      //initialize the quantumnumbers of the MPO<SpinQuantum>
      mpo[0].resize(SpinQuantum::zero(),qshape);

      qshape = make_array(qi,global::qp,-global::qp,qo);

      for(int i = 1;i < L-1;++i)
         mpo[i].resize(SpinQuantum::zero(),qshape);

      qshape = make_array(qi,global::qp,-global::qp,qz);

      mpo[L-1].resize(SpinQuantum::zero(),qshape);

      double mz = -sz;

      for(int m = 0;m < d;++m){

         // set block elements
         DArray<4> I_op(1, 1, 1, 1);//identity
         I_op = 1.0;

         DArray<4> Sz_op(1, 1, 1, 1);//Sz
         Sz_op = mz;

         mpo[0].insert(shape(0,m,m,0),I_op);
         mpo[0].insert(shape(0,m,m,1),Sz_op);

         for(int i = 1;i < L - 1;++i){

            mpo[i].insert(shape(0,m,m,0),I_op);
            mpo[i].insert(shape(0,m,m,1),Sz_op);
            mpo[i].insert(shape(1,m,m,1),I_op);

         }

         mpo[L-1].insert(shape(0,m,m,0),Sz_op);
         mpo[L-1].insert(shape(1,m,m,0),I_op);

         mz += 1.0;

      }

      //merge everything together
      TVector<Qshapes<SpinQuantum>,1> qmerge;
      TVector<Dshapes,1> dmerge;

      qmerge[0] = mpo[0].qshape(3);
      dmerge[0] = mpo[0].dshape(3);

      QSTmergeInfo<1> info(qmerge,dmerge);

      QSDArray<4> tmp;
      QSTmerge(mpo[0],info,tmp);

      mpo[0] = tmp;

      for(int i = 1;i < L - 1;++i){

         //first merge the row
         qmerge[0] = mpo[i].qshape(0);
         dmerge[0] = mpo[i].dshape(0);

         info.reset(qmerge,dmerge);

         tmp.clear();

         QSTmerge(info,mpo[i],tmp);

         //then merge the column
         qmerge[0] = tmp.qshape(3);
         dmerge[0] = tmp.dshape(3);

         info.reset(qmerge,dmerge);

         mpo[i].clear();

         QSTmerge(tmp,info,mpo[i]);

      }

      //only merge row for i = L - 1
      qmerge[0] = mpo[L - 1].qshape(0);
      dmerge[0] = mpo[L - 1].dshape(0);

      info.reset(qmerge,dmerge);

      tmp.clear();

      QSTmerge(info,mpo[L - 1],tmp);

      mpo[L - 1] = tmp;

      return mpo;

   }

   /**
    * initialize the MPO<SpinQuantum> to represent a spin raising operator
    * @param L length of the chain
    * @param d local dimension: i.e. defines the size of the local spins
    */
   MPO<SpinQuantum> raise(int L,int d){

      double sz = 0.5 * (d - 1.0);//size of local spin

      MPO<SpinQuantum> mpo(L);

      Qshapes<SpinQuantum> qz;
      qz.push_back(SpinQuantum::zero());

      Qshapes<SpinQuantum> qt;
      qt.push_back(SpinQuantum(2));

      //incoming
      Qshapes<SpinQuantum> qi;
      qi.push_back(SpinQuantum::zero());//I has spin 0
      qi.push_back(SpinQuantum(2));//S+ has spin 2

      //outgoing
      Qshapes<SpinQuantum> qo;
      qo.push_back(SpinQuantum::zero());//I has spin 0
      qo.push_back(SpinQuantum(-2));//S+ has spin 2

      TVector<Qshapes<SpinQuantum>,4> qshape = make_array(qz,global::qp,-global::qp,qo);

      //initialize the quantumnumbers of the MPO<SpinQuantum>
      mpo[0].resize(SpinQuantum::zero(),qshape);

      qshape = make_array(qi,global::qp,-global::qp,qo);

      for(int i = 1;i < L-1;++i)
         mpo[i].resize(SpinQuantum::zero(),qshape);

      qshape = make_array(qi,global::qp,-global::qp,-qt);

      mpo[L-1].resize(SpinQuantum::zero(),qshape);

      //first the identity
      double mz = -sz;

      for(int m = 0;m < d;++m){

         // set block elements
         DArray<4> I_op(1, 1, 1, 1);//identity
         I_op = 1.0;

         mpo[0].insert(shape(0,m,m,0),I_op);

         for(int i = 1;i < L - 1;++i){

            mpo[i].insert(shape(0,m,m,0),I_op);
            mpo[i].insert(shape(1,m,m,1),I_op);

         }

         mpo[L-1].insert(shape(1,m,m,0),I_op);

         mz += 1.0;

      }

      //then the raising operator
      mz = -sz;

      for(int m = 0; m < d - 1; ++m) {

         // set block elements
         DArray<4> Sp(1, 1, 1, 1);
         Sp = std::sqrt( (sz - mz) * (sz + mz + 1.0) );

         // insert blocks
         mpo[0].insert(shape(0,m+1,m,1),Sp);

         for(int i = 1; i < L-1; ++i) 
            mpo[i].insert(shape(0,m+1,m,1),Sp);

         mpo[L-1].insert(shape(0,m+1,m,0),Sp);

         mz  += 1.0;

      }

      return mpo;

   }

   /**
    * initialize the MPO<SpinQuantum> to represent a spin lowering operator
    * @param L length of the chain
    * @param d local dimension: i.e. defines the size of the local spins
    */
   MPO<SpinQuantum> lower(int L,int d){

      double sz = 0.5 * (d - 1.0);//size of local spin

      MPO<SpinQuantum> mpo(L);

      Qshapes<SpinQuantum> qz;
      qz.push_back(SpinQuantum::zero());

      Qshapes<SpinQuantum> qt;
      qt.push_back(SpinQuantum(-2));

      //incoming
      Qshapes<SpinQuantum> qi;
      qi.push_back(SpinQuantum::zero());//I has spin 0
      qi.push_back(SpinQuantum(-2));//S+ has spin 2

      //outgoing
      Qshapes<SpinQuantum> qo;
      qo.push_back(SpinQuantum::zero());//I has spin 0
      qo.push_back(SpinQuantum(2));//S+ has spin 2

      TVector<Qshapes<SpinQuantum>,4> qshape = make_array(qz,global::qp,-global::qp,qo);

      //initialize the quantumnumbers of the MPO<SpinQuantum>
      mpo[0].resize(SpinQuantum::zero(),qshape);

      qshape = make_array(qi,global::qp,-global::qp,qo);

      for(int i = 1;i < L-1;++i)
         mpo[i].resize(SpinQuantum::zero(),qshape);

      qshape = make_array(qi,global::qp,-global::qp,-qt);

      mpo[L-1].resize(SpinQuantum::zero(),qshape);

      //first the identity
      double mz = -sz;

      for(int m = 0;m < d;++m){

         // set block elements
         DArray<4> I_op(1, 1, 1, 1);//identity
         I_op = 1.0;

         mpo[0].insert(shape(0,m,m,0),I_op);

         for(int i = 1;i < L - 1;++i){

            mpo[i].insert(shape(0,m,m,0),I_op);
            mpo[i].insert(shape(1,m,m,1),I_op);

         }

         mpo[L-1].insert(shape(1,m,m,0),I_op);

         mz += 1.0;

      }

      //then the lowering operator
      mz = sz;

      for(int m = d - 1;m > 0;--m) {

         // set block elements
         DArray<4> Sm(1, 1, 1, 1);
         Sm = std::sqrt( (sz + mz) * (sz - mz + 1.0) );

         // insert blocks
         mpo[0].insert(shape(0,m-1,m,1),Sm);

         for(int i = 1; i < L-1; ++i) 
            mpo[i].insert(shape(0,m-1,m,1),Sm);

         mpo[L-1].insert(shape(0,m-1,m,0),Sm);

         mz  -= 1.0;

      }

      return mpo;

   }

   /**
    * initialize the MPO<SpinQuantum> to represent a nearest-neighbour ising Hamiltonian on a lattice of size L and with coupling constant J
    * @param L length of the chain
    * @param d local dimension: i.e. defines the size of the local spins
    * @param J coupling constant between S+S-
    * @param B magnetic fieldstrength in z direction
    */
   MPO<SpinQuantum> XY(int L,int d,double J,double B){

      double sz = 0.5 * (d - 1.0);//size of local spin

      MPO<SpinQuantum> mpo(L);

      Qshapes<SpinQuantum> qz;
      qz.push_back(SpinQuantum::zero());

      //incoming
      Qshapes<SpinQuantum> qi;
      qi.push_back(SpinQuantum::zero());//I has spin 0
      qi.push_back(SpinQuantum(2));//S- has spin -2
      qi.push_back(SpinQuantum(-2));//S+ has spin +2
      qi.push_back(SpinQuantum::zero());//B has spin 0

      //outgoing
      Qshapes<SpinQuantum> qo;
      qo.push_back(SpinQuantum::zero());//I has spin 0
      qo.push_back(SpinQuantum(-2));//S+ has spin 2
      qo.push_back(SpinQuantum(2));//S- has spin 2
      qo.push_back(SpinQuantum::zero());//B has spin 0

      TVector<Qshapes<SpinQuantum>,4> qshape = make_array(qz,global::qp,-global::qp,qo);

      //initialize the quantumnumbers of the MPO<SpinQuantum>
      mpo[0].resize(SpinQuantum::zero(),qshape);

      qshape = make_array(qi,global::qp,-global::qp,qo);

      for(int i = 1;i < L-1;++i)
         mpo[i].resize(SpinQuantum::zero(),qshape);

      qshape = make_array(qi,global::qp,-global::qp,qz);

      mpo[L-1].resize(SpinQuantum::zero(),qshape);

      //first the identity and local term
      double mz = -sz;

      for(int m = 0;m < d;++m){

         // set block elements
         DArray<4> I_op(1, 1, 1, 1);//identity
         I_op = 1.0;

         DArray<4> B_op(1, 1, 1, 1);//identity
         B_op = -B*mz;

         mpo[0].insert(shape(0,m,m,0),I_op);
         mpo[0].insert(shape(0,m,m,3),B_op);

         for(int i = 1;i < L - 1;++i){

            mpo[i].insert(shape(0,m,m,0),I_op);
            mpo[i].insert(shape(0,m,m,3),B_op);
            mpo[i].insert(shape(3,m,m,3),I_op);

         }

         mpo[L-1].insert(shape(0,m,m,0),B_op);
         mpo[L-1].insert(shape(3,m,m,0),I_op);

         mz += 1.0;

      }

      //then the raising operator
      mz = -sz;

      for(int m = 0; m < d - 1; ++m) {

         // set block elements
         DArray<4> Sp(1, 1, 1, 1);
         Sp = std::sqrt( J*(sz - mz) * (sz + mz + 1.0) );

         // insert blocks
         mpo[0].insert(shape(0,m+1,m,1),Sp);

         for(int i = 1;i < L - 1;++i){

            mpo[i].insert(shape(0,m+1,m,1),Sp);
            mpo[i].insert(shape(2,m+1,m,3),Sp);

         }

         mpo[L-1].insert(shape(2,m+1,m,0),Sp);

         mz  += 1.0;

      }

      //then the lowering operator
      mz = sz;

      for(int m = d-1; m > 0;--m) {

         // set block elements
         DArray<4> Sm(1, 1, 1, 1);
         Sm = std::sqrt( J*(sz + mz) * (sz - mz + 1.0) );

         // insert blocks
         mpo[0].insert(shape(0,m-1,m,2),Sm);

         for(int i = 1;i < L - 1;++i){

            mpo[i].insert(shape(0,m-1,m,2),Sm);
            mpo[i].insert(shape(1,m-1,m,3),Sm);

         }

         mpo[L-1].insert(shape(1,m-1,m,0),Sm);

         mz -= 1.0;

      }


      //merge everything together
      TVector<Qshapes<SpinQuantum>,1> qmerge;
      TVector<Dshapes,1> dmerge;

      qmerge[0] = mpo[0].qshape(3);
      dmerge[0] = mpo[0].dshape(3);

      QSTmergeInfo<1> info(qmerge,dmerge);

      QSDArray<4> tmp;
      QSTmerge(mpo[0],info,tmp);

      mpo[0] = tmp;

      for(int i = 1;i < L - 1;++i){

         //first merge the row
         qmerge[0] = mpo[i].qshape(0);
         dmerge[0] = mpo[i].dshape(0);

         info.reset(qmerge,dmerge);

         tmp.clear();

         QSTmerge(info,mpo[i],tmp);

         //then merge the column
         qmerge[0] = tmp.qshape(3);
         dmerge[0] = tmp.dshape(3);

         info.reset(qmerge,dmerge);

         mpo[i].clear();

         QSTmerge(tmp,info,mpo[i]);

      }

      //only merge row for i = L - 1
      qmerge[0] = mpo[L - 1].qshape(0);
      dmerge[0] = mpo[L - 1].dshape(0);

      info.reset(qmerge,dmerge);

      tmp.clear();

      QSTmerge(info,mpo[L - 1],tmp);

      mpo[L - 1] = tmp;

      return mpo;

   }

   /**
    * initialize the MPO<SpinQuantum> to represent a nearest-neighbour anisotropic Heisenberg Hamiltonian on a lattice of size L 
    * and with coupling constant Jz for the z spins and Jxy for XY
    * inside a magnetic field B which defines the z direction.
    * nearest neighbour interaction is repulsive, i.e. Jz and Jxy > 0.0!
    * @param L length of the chain
    * @param d local dimension: i.e. defines the size of the local spins
    * @param Jz coupling constant > 0 
    * @param Jxy coupling constant > 0
    * @param B magnetic fieldstrength
    */
   MPO<SpinQuantum> heisenberg(int L,int d,double Jz,double Jxy,double B){

      double sz = 0.5 * (d - 1.0);//size of local spin

      MPO<SpinQuantum> mpo(L);

      Qshapes<SpinQuantum> qz;
      qz.push_back(SpinQuantum::zero());

      //incoming
      Qshapes<SpinQuantum> qi;
      qi.push_back(SpinQuantum::zero());//I has spin 0
      qi.push_back(SpinQuantum(2));//S- has spin -2
      qi.push_back(SpinQuantum(-2));//S+ has spin +2
      qi.push_back(SpinQuantum::zero());//Sz has spin 0
      qi.push_back(SpinQuantum::zero());//B has spin 0

      //outgoing
      Qshapes<SpinQuantum> qo;
      qo.push_back(SpinQuantum::zero());//I has spin 0
      qo.push_back(SpinQuantum(-2));//S+ has spin 2
      qo.push_back(SpinQuantum(2));//S- has spin 2
      qo.push_back(SpinQuantum::zero());//Sz has spin 0
      qo.push_back(SpinQuantum::zero());//B has spin 0

      TVector<Qshapes<SpinQuantum>,4> qshape = make_array(qz,global::qp,-global::qp,qo);

      //initialize the quantumnumbers of the MPO<SpinQuantum>
      mpo[0].resize(SpinQuantum::zero(),qshape);

      qshape = make_array(qi,global::qp,-global::qp,qo);

      for(int i = 1;i < L-1;++i)
         mpo[i].resize(SpinQuantum::zero(),qshape);

      qshape = make_array(qi,global::qp,-global::qp,-qz);

      mpo[L-1].resize(SpinQuantum::zero(),qshape);

      //first the elements diagonal in the physical indices (B,Sz,I)
      double mz = -sz;

      for(int m = 0;m < d;++m){

         // set block elements
         DArray<4> I_op(1, 1, 1, 1);//identity
         I_op = 1.0;

         DArray<4> Sz_op(1, 1, 1, 1);//identity
         Sz_op = std::sqrt(Jz) * mz;

         DArray<4> B_op(1, 1, 1, 1);//identity
         B_op = -B*mz;

         mpo[0].insert(shape(0,m,m,0),I_op);
         mpo[0].insert(shape(0,m,m,3),Sz_op);
         mpo[0].insert(shape(0,m,m,4),B_op);

         for(int i = 1;i < L - 1;++i){

            mpo[i].insert(shape(0,m,m,0),I_op);
            mpo[i].insert(shape(0,m,m,3),Sz_op);
            mpo[i].insert(shape(0,m,m,4),B_op);
            mpo[i].insert(shape(3,m,m,4),Sz_op);
            mpo[i].insert(shape(4,m,m,4),I_op);

         }

         mpo[L-1].insert(shape(0,m,m,0),B_op);
         mpo[L-1].insert(shape(3,m,m,0),Sz_op);
         mpo[L-1].insert(shape(4,m,m,0),I_op);

         mz += 1.0;

      }

      //then the raising operator
      mz = -sz;

      for(int m = 0; m < d - 1; ++m) {

         // set block elements
         DArray<4> Sp(1, 1, 1, 1);
         Sp = std::sqrt( Jxy * (sz - mz) * (sz + mz + 1.0) );

         // insert blocks
         mpo[0].insert(shape(0,m+1,m,1),Sp);

         for(int i = 1;i < L - 1;++i){

            mpo[i].insert(shape(0,m+1,m,1),Sp);
            mpo[i].insert(shape(2,m+1,m,4),Sp);

         }

         mpo[L-1].insert(shape(2,m+1,m,0),Sp);

         mz  += 1.0;

      }

      //finally the lowering operator
      mz = sz;

      for(int m = d - 1;m > 0;--m) {

         // set block elements
         DArray<4> Sm(1, 1, 1, 1);
         Sm = std::sqrt( Jxy * (sz + mz) * (sz - mz + 1.0) );

         // insert blocks
         mpo[0].insert(shape(0,m-1,m,2),Sm);

         for(int i = 1;i < L - 1;++i){

            mpo[i].insert(shape(0,m-1,m,2),Sm);
            mpo[i].insert(shape(1,m-1,m,4),Sm);

         }

         mpo[L-1].insert(shape(1,m-1,m,0),Sm);

         mz -= 1.0;

      }

      //merge everything together
      TVector<Qshapes<SpinQuantum>,1> qmerge;
      TVector<Dshapes,1> dmerge;

      qmerge[0] = mpo[0].qshape(3);
      dmerge[0] = mpo[0].dshape(3);

      QSTmergeInfo<1> info(qmerge,dmerge);

      QSDArray<4> tmp;
      QSTmerge(mpo[0],info,tmp);

      mpo[0] = tmp;

      for(int i = 1;i < L - 1;++i){

         //first merge the row
         qmerge[0] = mpo[i].qshape(0);
         dmerge[0] = mpo[i].dshape(0);

         info.reset(qmerge,dmerge);

         tmp.clear();

         QSTmerge(info,mpo[i],tmp);

         //then merge the column
         qmerge[0] = tmp.qshape(3);
         dmerge[0] = tmp.dshape(3);

         info.reset(qmerge,dmerge);

         mpo[i].clear();

         QSTmerge(tmp,info,mpo[i]);

      }

      //only merge row for i = L - 1
      qmerge[0] = mpo[L - 1].qshape(0);
      dmerge[0] = mpo[L - 1].dshape(0);

      info.reset(qmerge,dmerge);

      tmp.clear();

      QSTmerge(info,mpo[L - 1],tmp);

      mpo[L - 1] = tmp;

      return mpo;

   }

}
