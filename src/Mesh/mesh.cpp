template<typename ArrayType, typename Tprec>
Mesh2d<ArrayType,Tprec>::Mesh2d(Tprec _lx, Tprec _ly, int _nx, int _ny)
: x(xx), y(yy)
{

   nx = _nx > 0 ? _nx : 1; 
   ny = _ny > 0 ? _ny : 1; 

   Lx = _lx > 0.0 ? _lx : 1.0; 
   Ly = _ly > 0.0 ? _ly : 1.0; 


   make();
}

template<typename ArrayType, typename Tprec>
Mesh2d<ArrayType,Tprec>::Mesh2d(Mesh2d &meshToCopy){
    nx = meshToCopy.nx;
    ny = meshToCopy.ny;

    xx.resize(nx,ny);
    yy.resize(nx,ny);

   for(int j = 0; j < ny; j++){
       for(int i = 0; i < nx; i++){
           xx(i,j) = meshToCopy.xx(i,j);
           yy(i,j) = meshToCopy.yy(i,j);
       } 
   }

}


template<typename ArrayType, typename Tprec>
void Mesh2d<ArrayType,Tprec>::make(){

   xx.resize(nx,ny);
   yy.resize(nx,ny);

   for(int i = 0; i < nx; i++){
       for(int j = 0; j < ny; j++){
           xx(i,j) = Lx*((Tprec)i)/ ((Tprec)nx -1.0);
           yy(i,j) = Ly*((Tprec)j)/ ((Tprec)ny -1.0);
       } 
   }

   dx = xx(1,0) - xx(0,0);
   dy = yy(0,1) - yy(0,0);
}


//-------3D

template<typename ArrayType, typename Tprec>
Mesh3d<ArrayType,Tprec>::Mesh3d(
                                Tprec _lx,  Tprec _ly, Tprec _lz, 
                                int    _nx, int   _ny, int   _nz
                                )
: Mesh2d<ArrayType,Tprec>(_lx, _ly, _nx, _ny),
  z(zz)
{

   nz = _nz > 0   ? _nz : 1; 
   Lz = _lz > 0.0 ? _lz : 1.0; 

   make();
}


template<typename ArrayType, typename Tprec>
Mesh3d<ArrayType,Tprec>::Mesh3d(Mesh3d &meshToCopy){


 
    this->nx = 10;//meshToCopy.nx;
    this->ny = meshToCopy.ny;
    nz = meshToCopy.nz;

    this->xx.resize(this->getNx(), this->getNy(), nz);
    this->yy.resize(this->getNx(), this->getNy(), nz);
    zz.resize(this->getNx(), this->getNy(), nz);

    this->xx = meshToCopy.xx;
    this->yy = meshToCopy.yy;
    zz = meshToCopy.zz;

}


template<typename ArrayType, typename Tprec>
void Mesh3d<ArrayType,Tprec>::make(){


   this->xx.resize(this->getNx(), this->getNy(), nz);
   this->yy.resize(this->getNx(), this->getNy(), nz);
   zz.resize(this->getNx(), this->getNy(), nz);

   for(int i = 0; i < this->nx; i++){
       for(int j = 0; j < this->ny; j++){
           for(int k = 0; k < nz; k++){
               this->xx(i,j,k) = 
               this->Lx*((Tprec)i)/ ((Tprec)this->nx - 1.0);

               this->yy(i,j,k) = 
               this->Ly*((Tprec)j)/ ((Tprec)this->ny - 1.0);

               zz(i,j,k) = Lz*((Tprec)k)/ ((Tprec)nz - 1.0);
           } 
       }
   }

   this->dx = this->xx(1,0,0) - this->xx(0,0,0);
   this->dy = this->yy(0,1,0) - this->yy(0,0,0);
   dz = zz(0,0,1) - zz(0,0,0);

}

/*
*/
