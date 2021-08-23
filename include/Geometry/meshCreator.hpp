#ifndef _MESHCREATOR_H
#define _MESHCREATOR_H

template<typename Tprec>
void regularMesh(
                Tprec *xx, Tprec *yy, Tprec *zz,
                Tprec Lx, Tprec Ly, Tprec Lz,
                Tprec Ox, Tprec Oy, Tprec Oz,
                unsigned int Nx, unsigned int Ny, unsigned int  Nz
                );
                
#include "meshCreator.cpp"
#endif
