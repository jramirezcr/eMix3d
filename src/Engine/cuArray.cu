#include<cuda.h>

//todo 

template<typename Tprec>
__host__ 
ArrayDev<Tprec>::ArrayDev(int _dimsize){

  dimsize = _dimsize > 0? _dimsize: 1;

  allocate();

  isCopy = false;

}

template<typename Tprec>
ArrayDev<Tprec>::ArrayDev(ArrayDev<Tprec> &_orig){

  *this = _orig;
  isCopy = true;

}


template<typename Tprec>
ArrayDev<Tprec>::~ArrayDev(){
     if(!isCopy){
       cudaFree(d_data);
     }

}

    
template<typename Tprec>
__host__
void ArrayDev<Tprec>::resize(int _dimsize){


      dimsize = _dimsize > 0? _dimsize: 1;
      allocate();
      isCopy = false;
}


template<typename Tprec>
__device__ 
Tprec& ArrayDev<Tprec>::operator()(int i) 
{

      return  d_data[i];
}


/*
template<typename Tprec>
__host__
ArrayDev<Tprec>& ArrayDev<Tprec>::operator=(ArrayDev<Tprec>& array){

     if(this != &array){
        dimsize = array.getDim();
        allocate();

        copyDevToDev(array.devPtr());
     }

   return *this;
}
*/

template<typename Tprec>
__host__ __device__ 
int ArrayDev<Tprec>::getDim() const { 
    return dimsize;
} 


template<typename Tprec>
__host__ 
Tprec* ArrayDev<Tprec>::devPtr(){
    return d_data; 
} 

template<typename Tprec>
__host__ void ArrayDev<Tprec>::allocate(){

  cudaError_t retult = 
  cudaMalloc((void**)&d_data, dimsize*sizeof(Tprec));
}


template<typename Tprec>
__host__ 
void ArrayDev<Tprec>::copyToHost(Tprec* hostPtr) {

   cudaMemcpy(hostPtr, d_data, dimsize*sizeof(Tprec), 
              cudaMemcpyDeviceToHost);
   
}

template<typename Tprec>
__host__ 
void ArrayDev<Tprec>::copyFromHost(Tprec* hostPtr){

   cudaMemcpy(d_data, hostPtr, dimsize*sizeof(Tprec), 
              cudaMemcpyHostToDevice);
   
}

template<typename Tprec>
__host__ 
void ArrayDev<Tprec>::copyDevToDev(Tprec* devPtr){

   cudaMemcpy(d_data, devPtr, dimsize*sizeof(Tprec), 
              cudaMemcpyDeviceToDevice);
   
}

template<typename Tprec>
__device__ 
Tprec& ArrayDev3D<Tprec>::operator()(int i, int j, int k) 
{
      return  ArrayDev<Tprec>::d_data[nx*ny*k + nx*j + i];
}


