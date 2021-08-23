#include<iostream>
#include"Engine/cuArray.hpp"

#ifndef __DERIVCU_HPP__
#define __DERIVCU_HPP__

template<typename Tprec>
__global__
void transposeArrayJ(ArrayDev3D<Tprec> a, ArrayDev3D<Tprec> b);


template<typename Tprec>
__global__
void transposeArrayK(ArrayDev3D<Tprec> a, ArrayDev3D<Tprec> b);


template<typename Tprec>
__global__
void getC4RHS(ArrayDev3D<Tprec> f,  ArrayDev3D<Tprec> r, Tprec delta);


template<typename Tprec>
__global__
void cuDerivCS(ArrayDev3D<Tprec> f, 
             ArrayDev3D<Tprec> r,
             ArrayDev3D<Tprec> gam
            );

//Deriv functions
namespace cu{
namespace cs{

template<typename Tprec>
void dFdX(
         ArrayDev3D<Tprec>  &fprima, 
         ArrayDev3D<Tprec>  &f, 
         Tprec               delta
         );
 

template<typename Tprec>
void dFdY(
         ArrayDev3D<Tprec>  &fprima, 
         ArrayDev3D<Tprec>  &f, 
         Tprec               delta
         );

template<typename Tprec>
void dFdZ(
          ArrayDev3D<Tprec>  &fprima, 
          ArrayDev3D<Tprec>  &f, 
          Tprec               delta
          );
//end namespaces
}
}

#include"Schemes/derivCompact.cu"

#endif
