
#ifndef VECTOR_H
#define VECTOR_H

#include <iostream>
#include <cmath>
#include <cstdlib>

template<typename Tprec>
class Vector
{
private:
  int dim;      // Dimension del vector
  Tprec *data;  // Puntero para almancenar los datos
  
public:
  Vector();              // Constructor por omision
  Vector(int);           // Constructor con argumentos
  Vector(int, Tprec);    // Constructor con argumentos (sobrecarga)
  Vector(const Vector<Tprec>&); // Constructor copia
  ~Vector();             // Destructor
  
//--- FUNCIONES MIEMBRO DE LA CLASE
  int getDim() const {return dim; }                  // Dimension
  void resize(int);                                  // Redimensionar
  Tprec calcNormL2() const;                          // Norma L-2 
  Tprec calcNormMax() const;                         // Norma infty
  Tprec& operator()(int i) const { return data[i]; } // v(i)

  Vector<Tprec>& operator=(const Vector<Tprec>&);                  
  // asignacion
  Vector<Tprec>& operator+=(const Vector<Tprec>&);                
  // v += v2
  Vector<Tprec>& operator-=(const Vector<Tprec>&);                
  // v -= v2

//--- FUNCIONES AMIGAS DE LA CLASE
  template<typename U>
  Vector<U> friend operator+(const Vector<U>&, const Vector<U>&); 
  // v = v1+v2
  template<typename U>
  Vector<U> friend operator-(const Vector<U>&, const Vector<U>&); 
  // v = v1-v2
  template<typename U>
  Vector<U> friend operator-(const Vector<U>&, const Vector<U>&); 
  template<typename U>
  U friend dot(const Vector<U>&, const Vector<U>&);       
  // dot product
  template<typename U>
  U friend operator*(const Vector<U>&, const Vector<U>&);  
  // scalar-vec 
  template<typename U>
  Vector<U> friend operator*(const U, const Vector<U>&);  
  // scalar-vec 
  template<typename U>
  Vector<U> friend operator*(const Vector<U>&, const U);  
  // vec-scalar

  //  friend Vector operator*(const Matriz&, const Vector&);  // vec-scalar

  template<typename U>
  friend std::ostream& operator<<(std::ostream&, const Vector<U>&); 
  // output
};

#include "../../src/Mesh/Vector.cpp"

#endif // VECTOR_H
