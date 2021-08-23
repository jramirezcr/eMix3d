/*

    Writen by Jorge Ramirez Cruz
    jramirezcr@iingen.unam.mx
    jramirezcr@gmail.com
    Instituto de Ingenieria, UNAM 

*/

#ifndef __ARRAY_H
#define __ARRAY_H

#include<iostream>
#include<cstdlib>

//TODO: resizing method
//TODO: 

template<int N_rank>
struct selector{
   enum{ value = N_rank};
};

template <typename Tprec, int N_rank>
class Array 
{

private:
   Tprec* data;
   int dim[N_rank];

   void store();

public:

   //Constructors

   Array()
   :data(0)
   {}

   Array(int lenght1)
   {
      dim[0] = lenght1; 
      if(1 != N_rank){std::cout<< "ARRAY WRONG Rank\n";exit(0);}
      store();
   }

   Array(int lenght1, int lenght2)
   {
      dim[0] = lenght1; 
      dim[1] = lenght2; 

      if(2 != N_rank){std::cout<< "ARRAY WRONG Rank\n";exit(0);}
      store();
   }

   Array(int lenght1, int lenght2, int lenght3)
   {
      dim[0] = lenght1; 
      dim[1] = lenght2; 
      dim[2] = lenght3; 

      if(3 != N_rank){std::cout<< "ARRAY WRONG Rank\n";exit(0);}
      store();
   }

   void resize(int lenght1)
   {
      dim[0] = lenght1; 
      if(1 != N_rank){std::cout<< "ARRAY WRONG Rank\n";exit(0);}
      store();
   }

   void resize(int lenght1, int lenght2)
   {
      dim[0] = lenght1; 
      dim[1] = lenght2; 

      if(2 != N_rank){std::cout<< "ARRAY WRONG Rank\n";exit(0);}
      store();
   }

   void resize(int lenght1, int lenght2, int lenght3)
   {
      dim[0] = lenght1; 
      dim[1] = lenght2; 
      dim[2] = lenght3; 

      if(3 != N_rank){std::cout<< "ARRAY WRONG Rank\n";exit(0);}
      store();

   }

   
   Array(const Array<Tprec,1>& other){

      dim[0] = other.dim[0]; 
      store();

      for(int i = 0; i < dim[0]; i++){
        data[i] = other.data[i]; 
      }
      
   }

   Array(const Array<Tprec,2>& other){

      dim[0] = other.dim[0]; 
      dim[1] = other.dim[1]; 
      store();

      for(int j = 0; j < dim[1]; j++){
         for(int i = 0; i < dim[0]; i++){
            data[dim[0]*j + i] = other.data[dim[0]*j + i]; 
         }
      }
      
   }


   Array(const Array<Tprec,3>& other){

      dim[0] = other.dim[0]; 
      dim[1] = other.dim[1]; 
      dim[2] = other.dim[2]; 
      store();

      for(int k = 0; k < dim[2]; k++){
         for(int j = 0; j < dim[1]; j++){
            for(int i = 0; i < dim[0]; i++){
               data[dim[0]*dim[1]*k + dim[0]*j + i] = 
               other.data[dim[0]*dim[1]*k + dim[0]*j + i]; 
            }
         }
      }
      
   }


   ~Array(){delete[] data;}



// Acces element data overloading

   Tprec& operator()(int i, 
                    selector<N_rank> dimType = selector<1>()){

       return data[i];
   }

   Tprec& operator()(int i, int j, 
                    selector<N_rank> dimType = selector<2>()){

       return data[dim[0]*j + i];
   }

   Tprec& operator()(int i, int j, int k,
                    selector<N_rank> dimType = selector<3>()){

       return data[dim[1]*dim[0]*k + dim[0]*j + i];
   }


   Array<Tprec, N_rank>& operator=(const Array<Tprec, 1>& other){
      if(this != &other) {
         for(int i = 0; i < dim[0]; i++){
           data[i] = other.data[i]; 
         }
      }

      return *this;
   }

   Array<Tprec, N_rank>& operator=(const Array<Tprec, 2>& other){
      if(this != &other) {
         for(int j = 0; j < dim[1]; j++){
            for(int i = 0; i < dim[0]; i++){
               data[dim[0]*j + i] = other.data[dim[0]*j + i]; 
            }
         }
      }

      return *this;
   }

   Array<Tprec, N_rank>& operator=(const Array<Tprec, 3>& other){
      if(this != &other) {
         for(int k = 0; k < dim[2]; k++){
            for(int j = 0; j < dim[1]; j++){
               for(int i = 0; i < dim[0]; i++){
                  data[dim[0]*dim[1]*k + dim[0]*j + i] = 
                  other.data[dim[0]*dim[1]*k + dim[0]*j + i]; 
               }
            }
         }
      }

      return *this;
   }

   void constantAsignation(Tprec value, selector<1>){

         for(int i = 0; i < dim[0]; i++){
           data[i] = value; 
         }
   }

   void constantAsignation(Tprec value, selector<2>){
         for(int j = 0; j < dim[1]; j++){
            for(int i = 0; i < dim[0]; i++){
               data[dim[0]*j + i] = value; 
            }
         }
   }

   void constantAsignation(Tprec value, selector<3>){
      for(int k = 0; k < dim[2]; k++){
         for(int j = 0; j < dim[1]; j++){
            for(int i = 0; i < dim[0]; i++){
               data[dim[0]*dim[1]*k + dim[0]*j + i] = value;
            }
         }
      }
   }

   Array<Tprec, N_rank>& operator=(Tprec value){
      constantAsignation(value, selector<N_rank>());
      return *this;
   }



//-----------------------------------------------------------------

   Array<Tprec, N_rank> operator*(const Array<Tprec, 1>& other){
      Array<Tprec, N_rank> result(other.dim[0]);  
      for(int i = 0; i < dim[0]; i++){
        result.data[i] = data[i]*other.data[i]; 
      }
      return result;
   }


   Array<Tprec, N_rank> operator*(const Array<Tprec, 2>& other){
      Array<Tprec, N_rank> result(dim[0], dim[1]);  

         for(int j = 0; j < dim[1]; j++){
            for(int i = 0; i < dim[0]; i++){
               result.data[dim[0]*j + i] = 
                      data[dim[0]*j + i]*other.data[dim[0]*j + i]; 
            }
         }


      return result;
   }


   Array<Tprec, N_rank> operator*(const Array<Tprec, 3>& other){

     Array<Tprec, N_rank> result(dim[0], dim[1], dim[2]);  

     for(int k = 0; k < dim[2]; k++){
        for(int j = 0; j < dim[1]; j++){
           for(int i = 0; i < dim[0]; i++){
              result.data[dim[0]*dim[1]*k + dim[0]*j + i] = 
                     data[dim[0]*dim[1]*k + dim[0]*j + i]*
                     other.data[dim[0]*dim[1]*k + dim[0]*j + i]; 
           }
        }
     }

      return result;
   }

   Array<Tprec, N_rank> operator*(const Tprec value);
   Array<Tprec, N_rank> operator/(const Tprec value);
   Array<Tprec, N_rank> operator+(const Tprec value);
   Array<Tprec, N_rank> operator-(const Tprec value);

   //Sublayer multiply overloading
   Array<Tprec, N_rank> constantMul(const Tprec value, selector<1>);
   Array<Tprec, N_rank> constantMul(const Tprec value, selector<2>);
   Array<Tprec, N_rank> constantMul(const Tprec value, selector<3>);

   //Sublayer divisor overloading
   Array<Tprec, N_rank> constantDiv(const Tprec value, selector<1>);
   Array<Tprec, N_rank> constantDiv(const Tprec value, selector<2>);
   Array<Tprec, N_rank> constantDiv(const Tprec value, selector<3>);

   //Sublayer Plus overloading
   Array<Tprec, N_rank> constantAdd(const Tprec value, selector<1>);
   Array<Tprec, N_rank> constantAdd(const Tprec value, selector<2>);
   Array<Tprec, N_rank> constantAdd(const Tprec value, selector<3>);

   //Sublayer Plus overloading
   Array<Tprec, N_rank> constantSub(const Tprec value, selector<1>);
   Array<Tprec, N_rank> constantSub(const Tprec value, selector<2>);
   Array<Tprec, N_rank> constantSub(const Tprec value, selector<3>);

//----------------------------------------------------------------
//-----------------------------------------------------------------

   Array<Tprec, N_rank> operator/(const Array<Tprec, 1>& other){
      Array<Tprec, N_rank> result(other.dim[0]);  
      for(int i = 0; i < dim[0]; i++){
        result.data[i] = data[i]/other.data[i]; 
      }
      return result;
   }


   Array<Tprec, N_rank> operator/(const Array<Tprec, 2>& other){
      Array<Tprec, N_rank> result(dim[0], dim[1]);  

         for(int j = 0; j < dim[1]; j++){
            for(int i = 0; i < dim[0]; i++){
               result.data[dim[0]*j + i] = 
                      data[dim[0]*j + i]/other.data[dim[0]*j + i]; 
            }
         }


      return result;
   }

   Array<Tprec, N_rank> operator/(const Array<Tprec, 3>& other){
      Array<Tprec, N_rank> result(dim[0], dim[1], dim[2]);  

     for(int k = 0; k < dim[2]; k++){
        for(int j = 0; j < dim[1]; j++){
           for(int i = 0; i < dim[0]; i++){
              result.data[dim[0]*dim[1]*k + dim[0]*j + i] = 
                     data[dim[0]*dim[1]*k + dim[0]*j + i]/
                     other.data[dim[0]*dim[1]*k + dim[0]*j + i]; 
           }
        }
     }

      return result;
   }

//----------------------------------------------------------------
//-----------------------------------------------------------------

   Array<Tprec, N_rank> operator-(const Array<Tprec, 1>& other){
      Array<Tprec, N_rank> result(other.dim[0]);  
      for(int i = 0; i < dim[0]; i++){
        result.data[i] = data[i]-other.data[i]; 
      }
      return result;
   }


   Array<Tprec, N_rank> operator-(const Array<Tprec, 2>& other){
      Array<Tprec, N_rank> result(dim[0], dim[1]);  

         for(int j = 0; j < dim[1]; j++){
            for(int i = 0; i < dim[0]; i++){
               result.data[dim[0]*j + i] = 
                      data[dim[0]*j + i]-other.data[dim[0]*j + i]; 
            }
         }


      return result;
   }

   Array<Tprec, N_rank> operator-(const Array<Tprec, 3>& other){
      Array<Tprec, N_rank> result(dim[0], dim[1], dim[2]);  

     for(int k = 0; k < dim[2]; k++){
        for(int j = 0; j < dim[1]; j++){
           for(int i = 0; i < dim[0]; i++){
              result.data[dim[0]*dim[1]*k + dim[0]*j + i] = 
                     data[dim[0]*dim[1]*k + dim[0]*j + i]-
                     other.data[dim[0]*dim[1]*k + dim[0]*j + i]; 
           }
        }
     }

      return result;
   }

//----------------------------------------------------------------
//-----------------------------------------------------------------

   Array<Tprec, N_rank> operator+(const Array<Tprec, 1>& other){
      Array<Tprec, N_rank> result(other.dim[0]);  
      for(int i = 0; i < dim[0]; i++){
        result.data[i] = data[i]+other.data[i]; 
      }
      return result;
   }


   Array<Tprec, N_rank> operator+(const Array<Tprec, 2>& other){
      Array<Tprec, N_rank> result(dim[0], dim[1]);  

         for(int j = 0; j < dim[1]; j++){
            for(int i = 0; i < dim[0]; i++){
               result.data[dim[0]*j + i] = 
                      data[dim[0]*j + i]+other.data[dim[0]*j + i]; 
            }
         }


      return result;
   }

   Array<Tprec, N_rank> operator+(const Array<Tprec, 3>& other){
      Array<Tprec, N_rank> result(dim[0], dim[1], dim[2]);  

     for(int k = 0; k < dim[2]; k++){
        for(int j = 0; j < dim[1]; j++){
           for(int i = 0; i < dim[0]; i++){
              result.data[dim[0]*dim[1]*k + dim[0]*j + i] = 
                     data[dim[0]*dim[1]*k + dim[0]*j + i]+
                     other.data[dim[0]*dim[1]*k + dim[0]*j + i]; 
           }
        }
     }

      return result;
   }


   int getDim(int _dim){

       //TODO: check dim for each dim type

       return dim[_dim - 1];
   }
};

template <typename Tprec, int N_rank>
void Array<Tprec,N_rank>::store()
{
   int offset = 1;
   for(int i = 0; i < N_rank; i++){
      offset *= dim[i]; 
   }  

   data = new Tprec[offset];
}


//---------------------------------------------------------------
// Overloading type A*value
//--------------------------------------------------------------

template <typename Tprec, int N_rank>
Array<Tprec, N_rank> Array<Tprec, N_rank>::operator*(const Tprec value)
{
   return constantMul(value, selector<N_rank>());
}


template <typename Tprec, int N_rank> Array<Tprec, N_rank> 
Array<Tprec, N_rank>::constantMul(const Tprec value, selector<1>){
   Array<Tprec, N_rank> result(dim[0]);  
   for(int i = 0; i < dim[0]; i++){
     result.data[i] = data[i]*value; 
   }
   return result;
}



template <typename Tprec, int N_rank> Array<Tprec, N_rank> 
Array<Tprec, N_rank>::constantMul(const Tprec value, selector<2>){
   Array<Tprec, N_rank> result(dim[0], dim[1]);  

      for(int j = 0; j < dim[1]; j++){
         for(int i = 0; i < dim[0]; i++){
            result.data[dim[0]*j + i] = 
                   data[dim[0]*j + i]*value; 
         }
      }


   return result;
}


template <typename Tprec, int N_rank> Array<Tprec, N_rank> 
Array<Tprec, N_rank>::constantMul(const Tprec value, selector<3>){

  Array<Tprec, N_rank> result(dim[0], dim[1], dim[2]);  
 
  for(int k = 0; k < dim[2]; k++){
     for(int j = 0; j < dim[1]; j++){
        for(int i = 0; i < dim[0]; i++){
           result.data[dim[0]*dim[1]*k + dim[0]*j + i] = 
                  data[dim[0]*dim[1]*k + dim[0]*j + i]*value;
        }
     }
   }
   return result;
}


//---------------------------------------------------------------
// Overloading type A/value
//--------------------------------------------------------------

template <typename Tprec, int N_rank>
Array<Tprec, N_rank> Array<Tprec, N_rank>::operator/(const Tprec value)
{
   return constantDiv(value, selector<N_rank>());
}


template <typename Tprec, int N_rank> Array<Tprec, N_rank> 
Array<Tprec, N_rank>::constantDiv(const Tprec value, selector<1>){
   Array<Tprec, N_rank> result(dim[0]);  
   for(int i = 0; i < dim[0]; i++){
     result.data[i] = data[i]/value; 
   }
   return result;
}



template <typename Tprec, int N_rank> Array<Tprec, N_rank> 
Array<Tprec, N_rank>::constantDiv(const Tprec value, selector<2>){
   Array<Tprec, N_rank> result(dim[0], dim[1]);  

      for(int j = 0; j < dim[1]; j++){
         for(int i = 0; i < dim[0]; i++){
            result.data[dim[0]*j + i] = 
                   data[dim[0]*j + i]/value; 
         }
      }


   return result;
}


template <typename Tprec, int N_rank> Array<Tprec, N_rank> 
Array<Tprec, N_rank>::constantDiv(const Tprec value, selector<3>){

  Array<Tprec, N_rank> result(dim[0], dim[1], dim[2]);  
 
  for(int k = 0; k < dim[2]; k++){
     for(int j = 0; j < dim[1]; j++){
        for(int i = 0; i < dim[0]; i++){
           result.data[dim[0]*dim[1]*k + dim[0]*j + i] = 
                  data[dim[0]*dim[1]*k + dim[0]*j + i]/value;
        }
     }
  }

   return result;

}

//---------------------------------------------------------------
// Overloading type A + value
//--------------------------------------------------------------

template <typename Tprec, int N_rank>
Array<Tprec, N_rank> Array<Tprec, N_rank>::operator+(const Tprec value)
{
   return constantAdd(value, selector<N_rank>());
}


template <typename Tprec, int N_rank> Array<Tprec, N_rank> 
Array<Tprec, N_rank>::constantAdd(const Tprec value, selector<1>){
   Array<Tprec, N_rank> result(dim[0]);  
   for(int i = 0; i < dim[0]; i++){
     result.data[i] = data[i]+value; 
   }
   return result;
}



template <typename Tprec, int N_rank> Array<Tprec, N_rank> 
Array<Tprec, N_rank>::constantAdd(const Tprec value, selector<2>){
   Array<Tprec, N_rank> result(dim[0], dim[1]);  

      for(int j = 0; j < dim[1]; j++){
         for(int i = 0; i < dim[0]; i++){
            result.data[dim[0]*j + i] = 
                   data[dim[0]*j + i]+value; 
         }
      }


   return result;
}


template <typename Tprec, int N_rank> Array<Tprec, N_rank> 
Array<Tprec, N_rank>::constantAdd(const Tprec value, selector<3>){

  Array<Tprec, N_rank> result(dim[0], dim[1], dim[2]);  
 
  for(int k = 0; k < dim[2]; k++){
     for(int j = 0; j < dim[1]; j++){
        for(int i = 0; i < dim[0]; i++){
           result.data[dim[0]*dim[1]*k + dim[0]*j + i] = 
                  data[dim[0]*dim[1]*k + dim[0]*j + i]+value;
        }
     }
  }

   return result;
}

//---------------------------------------------------------------
// Overloading type A - value
//--------------------------------------------------------------

template <typename Tprec, int N_rank>
Array<Tprec, N_rank> Array<Tprec, N_rank>::operator-(const Tprec value)
{
   return constantSub(value, selector<N_rank>());
}


template <typename Tprec, int N_rank> Array<Tprec, N_rank> 
Array<Tprec, N_rank>::constantSub(const Tprec value, selector<1>){
   Array<Tprec, N_rank> result(dim[0]);  
   for(int i = 0; i < dim[0]; i++){
     result.data[i] = data[i] - value; 
   }
   return result;
}



template <typename Tprec, int N_rank> Array<Tprec, N_rank> 
Array<Tprec, N_rank>::constantSub(const Tprec value, selector<2>){
   Array<Tprec, N_rank> result(dim[0], dim[1]);  

      for(int j = 0; j < dim[1]; j++){
         for(int i = 0; i < dim[0]; i++){
            result.data[dim[0]*j + i] = 
                   data[dim[0]*j + i] - value; 
         }
      }


   return result;
}


template <typename Tprec, int N_rank> Array<Tprec, N_rank> 
Array<Tprec, N_rank>::constantSub(const Tprec value, selector<3>){

  Array<Tprec, N_rank> result(dim[0], dim[1], dim[2]);  
 
  for(int k = 0; k < dim[2]; k++){
     for(int j = 0; j < dim[1]; j++){
        for(int i = 0; i < dim[0]; i++){
           result.data[dim[0]*dim[1]*k + dim[0]*j + i] = 
                  data[dim[0]*dim[1]*k + dim[0]*j + i] - value;
        }
     }
  }

   return result;

}

#endif

