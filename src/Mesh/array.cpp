template<typename Tprec>
Array2d<Tprec>::Array2d(){

   nx = 1; 
   ny = 1; 

   data = new Tprec[nx*ny];

}


template<typename Tprec>
Array2d<Tprec>::Array2d(int _nx, int _ny){

   nx = _nx > 0 ? _nx : 1; 
   ny = _ny > 0 ? _ny : 1; 

   data = new Tprec[nx*ny];

}


template<typename Tprec>
Array2d<Tprec>::Array2d(const Array2d<Tprec>& _array){

    if(this != &_array){
        nx = _array.nx;
        ny = _array.ny;
        resize(nx,ny);
        int id = 0;
        for(int i = 0; i < nx; i++){
           for(int j = 0; j < ny; j++){
              id = i + j*nx;
              data[id] = _array.data[id];
           }
        }
    }

}
template<typename Tprec>
void Array2d<Tprec>::resize(int _nx, int _ny){
   
    delete [] data;

    nx = _nx > 0 ? _nx : 1; 
    ny = _ny > 0 ? _ny : 1; 

    data = new Tprec[nx*ny]; 
    
}


template<typename Tprec>
int Array2d<Tprec>::getDim(int sizeDim) const{
    if(1 == sizeDim){
        return nx;
    }
    else if(2 == sizeDim){
        return ny;
    }
    else{
        return -1;
    }
}



template<typename Tprec>
inline Tprec &Array2d<Tprec>::operator()(int i, int j) const{
    return data[i + nx*j];
}


// Field members


template<typename Tprec>
Field1d<Tprec>::Field1d(){

   nx = 1; 

   data = new Tprec[nx];

}


template<typename Tprec>
Field1d<Tprec>::Field1d(int _nx){

   nx = _nx > 0 ? _nx : 1; 

   data = new Tprec[nx];

}

template<typename Tprec>
void Field1d<Tprec>::resize(int _nx){
   
    delete [] data;

    nx = _nx > 0 ? _nx : 1; 

    data = new Tprec[nx]; 
    
}


template<typename Tprec>
Field1d<Tprec>::Field1d(const Field1d<Tprec>& _array){

    if(this != &_array){
        nx = _array.nx;
        resize(nx);
        for(int i = 0; i < nx; i++){
           data[i] = _array.data[i];
        }
    }

}

template<typename Tprec>
int Field1d<Tprec>::getDim(int sizeDim) const{
    if(1 == sizeDim){
        return nx;
    }
    else{
        return -1;
    }
}


template<typename Tprec>
Tprec& Field1d<Tprec>::operator()(int i) const{
    return data[i];
}



