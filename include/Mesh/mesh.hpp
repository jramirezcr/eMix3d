#ifndef _Mesh_HPP
#define _Mesh_HPP
#include <iostream>

#include "Engine/array.hpp"
#include "Engine/cuArray.hpp"


template<typename ArrayType, typename Tprec>
class Mesh2d{


protected:
    int nx;                // #cells Axis-x
    int ny;                // #cells Axis-y

    Tprec Lx;              // Axis X distace
    Tprec Ly;              // Axis Y distance

    Tprec dx;
    Tprec dy;

    ArrayType xx;
    ArrayType yy;

    virtual void make();

public:
   Mesh2d():x(xx), y(yy){};
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

   const ArrayType& x;
   const ArrayType& y;
};


//-------3D

template<typename ArrayType, typename Tprec>
class Mesh3d{


protected:

    int nx;                // #cells Axis-x
    int ny;                // #cells Axis-y
    int nz;


    Tprec dx;
    Tprec dy;
    Tprec dz;

    Tprec Lx;
    Tprec Ly;
    Tprec Lz;


   ArrayType xx;
   ArrayType yy;
   ArrayType zz;

public:

   Mesh3d():x(xx), y(yy), z(zz){};
   Mesh3d(Tprec, Tprec, Tprec, int, int, int);

   ~Mesh3d(){};

   inline int getNx()const {return nx;}
   inline int getNy()const {return ny;}
   inline int getNz()const {return nz;}

   inline Tprec getDeltaX( ) const {return dx;}
   inline Tprec getDeltaY( ) const {return dy;}
   inline Tprec getDeltaZ( ) const {return dz;}

   inline Tprec getLx( ) const {return Lx;}
   inline Tprec getLy( ) const {return Ly;}
   inline Tprec getLz( ) const {return Lz;}

   inline void setNx(int _nx){nx = _nx > 0   ? _nx : 1;}
   inline void setNy(int _ny){ny = _ny > 0   ? _ny : 1;}
   inline void setNz(int _nz){nz = _nz > 0   ? _nz : 1;}

   inline void setLx(Tprec _lx){Lx = _lx != 0.0   ? _lx : 1.0;}
   inline void setLy(Tprec _ly){Ly = _ly != 0.0   ? _ly : 1.0;}
   inline void setLz(Tprec _lz){Lz = _lz != 0.0   ? _lz : 1.0;}

   const ArrayType& x;
   const ArrayType& y;
   const ArrayType& z;

};

template<typename ArrayType, typename Tprec>
Mesh3d<ArrayType,Tprec>::Mesh3d(Tprec _lx, Tprec _ly, Tprec _lz,  
                                int _nx, int _ny, int _nz)

:x(xx), y(yy), z(zz)
{


   setLx(_lx);
   setLy(_lx);
   setLz(_lx);

   setNx(_nx);
   setNy(_nx);
   setNz(_nx);

}


//#include "Mesh/mesh.cpp"

#endif



