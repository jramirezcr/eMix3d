#ifndef _stlIO_H_
#define _stlIO_H_

#include <fstream>
#include <iostream>

using namespace std;

template<typename Tprec>
struct Vector3d{
   Tprec *x;
   Tprec *y;
   Tprec *z;
};

template<typename Tprec>
struct triangleSTL{
   Tprec vrtx1[3];
   Tprec vrtx2[3];
   Tprec vrtx3[3];
   Tprec n[3];
};

template<typename Tprec>
class StlIO {

public:
   StlIO(const char _name[]); 
   void getNameSTL() const;
   unsigned int getNumTrngl() const;
   void read();
   void getVrtx(triangleSTL<Tprec>*);
   triangleSTL<Tprec> *gmtry;
   Tprec getMaxAxisX() const;
   void rotateAxisZ(float alpha);

private:
   char *file_name, header[81];
   std::ifstream fileSTL;
   unsigned int num_facets;
};

#include"stlIO.cpp"

#endif
