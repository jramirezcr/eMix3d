#ifndef _CUARMMesh_HPP
#define _CUARMMesh_HPP
#include <iostream>


template<typename Tprec>
class Mesh{

protected:

    virtual void make(){};
    
public:
    Mesh(){}
    ~Mesh(){};
};


template<typename ArrayType, typename Tprec>
class Mesh2d : public Mesh<Tprec>{


protected:
    int nx;                // #cells Axis-x
    int ny;                // #cells Axis-y

    Tprec Lx;              // Axis X distace
    Tprec Ly;              // Axis Y distance

    Tprec dx;
    Tprec dy;

    


public:
   Mesh2d(Tprec _lx, Tprec _ly, int _nx, int _ny); // _lx, _ly - longitud 
   ~Mesh2d(){};

   inline int getNx()const {return nx;}
   inline int getNy()const {return ny;}

   inline Tprec getDeltaX( ) const {return dx;}
   inline Tprec getDeltaY( ) const {return dy;}

   inline Tprec getLx( ) const {return Lx;}
   inline Tprec getLy( ) const {return Ly;}

   inline void setNx(int _nx){nx = _nx > 0   ? _nx : 1;}
   inline void setNy(int _ny){ny = _ny > 0   ? _ny : 1;}
   inline void setLx(Tprec _Lx){Lx = _Lx > 0.0   ? _Lx : 1.0;}
   inline void setLy(Tprec _Ly){Ly = _Ly > 0.0   ? _Ly : 1.0;}

   void solve(Tprec _centerx, Tprec _centery);

   Mesh2d* child;
   Mesh2d* sibil;
   int     level;
   int     id;

    ArrayType x;
    ArrayType y;

};

template<typename ArrayType, typename Tprec>
Mesh2d<ArrayType,Tprec>::Mesh2d(Tprec _lx, Tprec _ly, int _nx, int _ny)
:x(_nx,_ny),y(_nx,_ny)
{
   setLx(_lx);
   setLy(_ly);
   setNx(_nx);
   setNy(_ny);

   dx = Lx / (float)(nx - 1);
   dy = Ly / (float)(ny - 1);
}

template<typename ArrayType, typename Tprec>
void Mesh2d<ArrayType,Tprec>::solve(Tprec _centerx, 
                                    Tprec _centery){

    for(int j = 0; j < ny; j++){
        for(int i = 0; i < nx; i++){
           x(i,j) = dx*(Tprec)i + _centerx; 
           y(i,j) = dy*(Tprec)j + _centery ;
        }
    }

}

//#include "Mesh/cuMesh.cpp"

#endif
