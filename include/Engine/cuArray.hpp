#ifndef _CUARRAYDEV_HPP
#define _CUARRAYDEV_HPP

#include<cuda.h>

//todo 

template<typename Tprec>
class ArrayDev{

    bool isCopy;

protected:

public:

    Tprec *d_data;
    int dimsize;

public:

    __host__ 
    explicit ArrayDev(): d_data(0), dimsize(0){}

    __host__ 
    explicit ArrayDev(int _dimsize);

    ArrayDev(ArrayDev<Tprec>&);

    virtual ~ArrayDev();


    //Services
     __host__
     void resize(int _dimsize);         //Change size array

    __host__ __device__ 
    int getDim() const;           //Get array dim 

    __host__ 
    Tprec* devPtr();                    
    //Return dev pointer to device data

    __host__ 
    void copyToHost(Tprec* hostPtr);
    // Copy array data to host array

    __host__ 
    void copyFromHost(Tprec* hostPtr);  
   // Copy data from a host array

    __host__
    void copyDevToDev(Tprec* devPtr);

    //Overloaded 
    //    __host__ 
    //ArrayDev<Tprec>& operator=(ArrayDev<Tprec>& );

    __device__ 
    Tprec& operator()(int i); //Allows notation syle A(i) 

private:

    __host__ 
    void allocate();

    __host__ 
    void deallocate();
};

template<typename Tprec>
class ArrayDev3D
: public ArrayDev<Tprec>{

private:

   int nx;
   int ny;
   int nz;

public:
   explicit ArrayDev3D(): ArrayDev<Tprec>(){}
   explicit ArrayDev3D(int i, int j, int k)
   :ArrayDev<Tprec>(){

      nx = i > 0 ? i : 1;
      ny = j > 0 ? j : 1;
      nz = k > 0 ? k : 1;

      ArrayDev<Tprec>::resize(nx*ny*nz);
   }

   void resize(int i, int j, int k){

      nx = i > 0 ? i : 1;
      ny = j > 0 ? j : 1;
      nz = k > 0 ? k : 1;

      ArrayDev<Tprec>::resize(nx*ny*nz);

   }


   virtual ~ArrayDev3D(){};

    __device__ 
    Tprec& operator()(int i, int j, int k); 
    //Allows notation syle A(i,j,k) 

    __host__ __device__
    int getDim(int _dim){
    if(_dim==1)     {return nx;}
    else if(_dim==2){return ny;}
    else if(_dim==3){return nz;}
    else            {return 0 ;}
    }

};

#include"Engine/cuArray.cu"

#endif
