
template<typename Tprec>
void regularMesh(
                Tprec *xx, Tprec *yy, Tprec *zz,
                Tprec Lx, Tprec Ly, Tprec Lz,
                Tprec Ox, Tprec Oy, Tprec Oz,
                unsigned int Nx, unsigned int Ny, unsigned int  Nz
                ){
   Tprec dx = Lx / (Tprec)(Nx - 1);
   Tprec dy = Ly / (Tprec)(Ny - 1);
   Tprec dz = Lz / (Tprec)(Nz - 1);
   unsigned int id = 0;

   for(unsigned int k = 0; k <= Nz - 1; k++){
      for(unsigned int j = 0; j <= Ny - 1; j++){
         for(unsigned int i = 0; i <= Nx - 1; i++){
            id = i + Nx*j + Nx*Ny*k;
            xx[id] = Lx*((Tprec)i- 1.0) / ((Tprec)Nx - 1.0) + Ox;
            yy[id] = Ly*((Tprec)j- 1.0) / ((Tprec)Ny - 1.0) + Oy;
            zz[id] = Lz*((Tprec)k- 1.0) / ((Tprec)Nz - 1.0) + Oz;
         }
      }
   }


}
