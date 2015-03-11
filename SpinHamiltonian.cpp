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
   MPO<SpinQuantum> ising(double J,double B){

      MPO<SpinQuantum> mpo(global::L);

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

      for(int i = 1;i < global::L-1;++i)
         mpo[i].resize(SpinQuantum::zero(),qshape);

      qshape = make_array(qi,global::qp,-global::qp,-qz);

      mpo[global::L-1].resize(SpinQuantum::zero(),qshape);

      //first site
      insert_id(mpo[0],0,0,1.0);//I
      insert_Sz(mpo[0],0,1,1.0);//Sz_i Sz_i+1
      insert_Sz(mpo[0],0,2,-B);//B

      //middle sites
      for(int i = 1;i < global::L - 1;++i){

         //for next site
         insert_id(mpo[i],0,0,1.0);//I
         insert_Sz(mpo[i],0,1,1.0);//Sz_i Sz_i+1
         insert_Sz(mpo[i],0,2,-B);//B

         //close down previous sites
         insert_Sz(mpo[i],1,2,J);//Sz_i Sz_i+1
         insert_id(mpo[i],2,2,1.0);//id

      }

      //last site
     
      //close down previous sites
      insert_Sz(mpo[global::L-1],0,0,-B);//B
      insert_Sz(mpo[global::L-1],1,0,J);//Sz_i Sz_i+1
      insert_id(mpo[global::L-1],2,0,1.0);//id

      //merge everything together
      TVector<Qshapes<SpinQuantum>,1> qmerge;
      TVector<Dshapes,1> dmerge;

      qmerge[0] = mpo[0].qshape(3);
      dmerge[0] = mpo[0].dshape(3);

      QSTmergeInfo<1> info(qmerge,dmerge);

      QSDArray<4> tmp;
      QSTmerge(mpo[0],info,tmp);

      mpo[0] = tmp;

      for(int i = 1;i < global::L - 1;++i){

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
      qmerge[0] = mpo[global::L - 1].qshape(0);
      dmerge[0] = mpo[global::L - 1].dshape(0);

      info.reset(qmerge,dmerge);

      tmp.clear();

      QSTmerge(info,mpo[global::L - 1],tmp);

      mpo[global::L - 1] = tmp;

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
   MPO<SpinQuantum> heisenberg(double Jxy,double Jz,double B){

      MPO<SpinQuantum> mpo(global::L);

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

      for(int i = 1;i < global::L-1;++i)
         mpo[i].resize(SpinQuantum::zero(),qshape);

      qshape = make_array(qi,global::qp,-global::qp,-qz);

      mpo[global::L-1].resize(SpinQuantum::zero(),qshape);

      //fill the mpo!

      //fist site
      insert_id(mpo[0],0,0,1.0);//I
      insert_Sp(mpo[0],0,1,1.0);//S+_i S-_i+1
      insert_Sm(mpo[0],0,2,1.0);//S-_i S+_i+1
      insert_Sz(mpo[0],0,3,1.0);//Sz_i Sz_i+1
      insert_Sz(mpo[0],0,4,-B);//B

      //middle sites
      for(int i = 1;i < global::L - 1;++i){

         //for next site
         insert_id(mpo[i],0,0,1.0);//I
         insert_Sp(mpo[i],0,1,1.0);//S+_i S-_i+1
         insert_Sm(mpo[i],0,2,1.0);//S-_i S+_i+1
         insert_Sz(mpo[i],0,3,1.0);//Sz_i Sz_i+1
         insert_Sz(mpo[i],0,4,-B);//B
 
         //close down previous sites
         insert_Sm(mpo[i],1,4,0.5 * Jxy);
         insert_Sp(mpo[i],2,4,0.5 * Jxy);
         insert_Sz(mpo[i],3,4,Jz);
         insert_id(mpo[i],4,4,1.0);//id

      }

      //close down stuff on last site
      insert_Sz(mpo[global::L-1],0,0,-B);//B
      insert_Sm(mpo[global::L-1],1,0,0.5 * Jxy);
      insert_Sp(mpo[global::L-1],2,0,0.5 * Jxy);
      insert_Sz(mpo[global::L-1],3,0,Jz);
      insert_id(mpo[global::L-1],4,0,1.0);//id

      //merge everything together
      TVector<Qshapes<SpinQuantum>,1> qmerge;
      TVector<Dshapes,1> dmerge;

      qmerge[0] = mpo[0].qshape(3);
      dmerge[0] = mpo[0].dshape(3);

      QSTmergeInfo<1> info(qmerge,dmerge);

      QSDArray<4> tmp;
      QSTmerge(mpo[0],info,tmp);

      mpo[0] = tmp;

      for(int i = 1;i < global::L - 1;++i){

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
      qmerge[0] = mpo[global::L - 1].qshape(0);
      dmerge[0] = mpo[global::L - 1].dshape(0);

      info.reset(qmerge,dmerge);

      tmp.clear();

      QSTmerge(info,mpo[global::L - 1],tmp);

      mpo[global::L - 1] = tmp;

      return mpo;

   }

   /**
    * insert identity operator in mpo O
    */
   void insert_id(QSDArray<4> &O,int row,int col,double value){

      DArray<4> Ip(1,1,1,1);
      Ip = value;

      for(int m = 0;m < global::d;++m)
         O.insert(shape(row,m,m,col),Ip);

   }

   /**
    * insert Sz operator in mpo O
    */
   void insert_Sz(QSDArray<4> &O,int row,int col,double value){

      DArray<4> Sz_op(1, 1, 1, 1);

      double sz = 0.5 * (global::d - 1.0);//size of local spin

      double mz = -sz;

      for(int m = 0;m < global::d;++m){

         Sz_op = mz * value;
         O.insert(shape(row,m,m,col),Sz_op);

         mz += 1.0;

      }

   }

   /**
    * insert Sp operator in mpo O
    */
   void insert_Sp(QSDArray<4> &O,int row,int col,double value){

      DArray<4> Sp_op(1, 1, 1, 1);

      double sz = 0.5 * (global::d - 1.0);//size of local spin

      double mz = -sz;

      for(int m = 0;m < global::d - 1;++m){

         Sp_op = value * std::sqrt( (sz - mz) * (sz + mz + 1.0) );

         O.insert(shape(row,m+1,m,col),Sp_op);

         mz += 1.0;

      }

   }

   /**
    * insert Sm operator in mpo O
    */
   void insert_Sm(QSDArray<4> &O,int row,int col,double value){

      DArray<4> Sm_op(1, 1, 1, 1);

      double sz = 0.5 * (global::d - 1.0);//size of local spin

      double mz = -sz + 1;

      for(int m = 1;m < global::d;++m){

         Sm_op = value * std::sqrt( (sz + mz) * (sz - mz + 1.0) );

         O.insert(shape(row,m-1,m,col),Sm_op);

         mz += 1.0;

      }

   }

}
