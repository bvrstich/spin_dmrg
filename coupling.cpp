#include <iostream>
#include <iomanip>

using std::cout;
using std::endl;
using std::ostream;

#include "include.h"

namespace coupling {

   /**
    * fill the TArray<double,double,2> object J with the coupling matrix for a 1D J1J2 model
    * @param pbc flag for periodic or open boundary conditions
    * @param J2 the J2 interaction strength as a function of J1
    * @param J the DArray object, will be filled with correct numbers on output
    */
   void J1J2_1D(bool pbc,double J2,TArray<double,2> &J){

      int L = J.shape(0);

      J = 0.0;

      //nn-interaction
      for(int i = 0;i < L - 1;++i){

         J(i,i + 1) = 1.0;
         J(i + 1,i) = 1.0;

      }

      //nnn-interaction
      for(int i = 0;i < L - 2;++i){

         J(i,i + 2) = J2;
         J(i + 2,i) = J2;

      }

      if(pbc){

         //nn
         J(0,L-1) = 1.0;
         J(L-1,0) = 1.0;

         //nnn
         J(0,L-2) = J2;
         J(L-2,0) = J2;

         J(1,L-1) = J2;
         J(L-1,1) = J2;

      }

   }


   /**
    * set the coupling matrix for the J1J2 model with J2 given as argument in units of J1
    * @param pbc periodic boundary conditions or not
    * @param J2 coupling strength
    * @param J coupling matrix, TArray<double,2> so matrix 
    */
   void J1J2_2D(bool pbc,double J2,TArray<double,2> &J){

      J = 0.0;

      int L = sqrt(J.shape(0));

      if(pbc){

         for (int row=0; row< L; row++)
            for (int col=0; col< L; col++){

               int number = row + L * col;

               int neighbour1 = (row + 1       )%L + L * col;
               int neighbour2 = (row - 1 + L)%L + L * col;
               int neighbour3 = row                   + L * ((col - 1 + L)%L);
               int neighbour4 = row                   + L * ((col + 1       )%L);

               J(number,neighbour1) = 1.0;
               J(number,neighbour2) = 1.0;
               J(number,neighbour3) = 1.0;
               J(number,neighbour4) = 1.0;

               int neighbour5 = (row + 1       )%L + L * ((col + 1       )%L);
               int neighbour6 = (row + 1       )%L + L * ((col - 1 + L)%L);
               int neighbour7 = (row - 1 + L)%L + L * ((col - 1 + L)%L);
               int neighbour8 = (row - 1 + L)%L + L * ((col + 1       )%L);

               J(number,neighbour5) = J2;
               J(number,neighbour6) = J2;
               J(number,neighbour7) = J2;
               J(number,neighbour8) = J2;

            }


      }
      else{//obc

         for (int row=0; row< L; row++)
            for (int col=0; col< L; col++){

               int number = row + L * col;

               int neighbour1 = (row + 1       )%L + L * col;
               int neighbour2 = (row - 1 + L)%L + L * col;
               int neighbour3 = row                   + L * ((col - 1 + L)%L);
               int neighbour4 = row                   + L * ((col + 1       )%L);

               J(number,neighbour1) = 1.0;
               J(number,neighbour2) = 1.0;
               J(number,neighbour3) = 1.0;
               J(number,neighbour4) = 1.0;

               int neighbour5 = (row + 1       )%L + L * ((col + 1       )%L);
               int neighbour6 = (row + 1       )%L + L * ((col - 1 + L)%L);
               int neighbour7 = (row - 1 + L)%L + L * ((col - 1 + L)%L);
               int neighbour8 = (row - 1 + L)%L + L * ((col + 1       )%L);

               J(number,neighbour5) = J2;
               J(number,neighbour6) = J2;
               J(number,neighbour7) = J2;
               J(number,neighbour8) = J2;

            }

         //no coupling from row 0 -> L-1
         for(int col1 = 0;col1 < L;++col1)
            for(int col2 = 0;col2 < L;++col2){

               int num1 = L * col1;//row 0
               int num2 = L-1 + L * col2;//row L - 1

               J(num1,num2) = 0.0;
               J(num2,num1) = 0.0;

            }

         //no coupling from col 0 -> L-1
         for(int row1 = 0;row1 < L;++row1)
            for(int row2 = 0;row2 < L;++row2){

               int num1 = row1;//col 0
               int num2 = row2 + L * (L-1);

               J(num1,num2) = 0.0;
               J(num2,num1) = 0.0;

            }

      }

   }

}
