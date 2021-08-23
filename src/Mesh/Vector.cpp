#include "Vector.hpp"
//#include "Matriz.hpp"
//
//--- CONSTRUCTORES
//
template<typename Tprec>
Vector<Tprec>::Vector()
{
  dim = 1;
  data = new Tprec[1];
}

template<typename Tprec>
Vector<Tprec>::Vector(int n)
{
  dim = n;
  data = new Tprec[dim];
}

template<typename Tprec>
Vector<Tprec>::Vector(int n, Tprec a)
{
  dim = n;
  data = new Tprec[dim];
  for (int i = 0; i < dim; i++)  data[i] = a;
}

template<typename Tprec>
Vector<Tprec>::Vector(const Vector<Tprec>& otro) 
{
  dim = otro.dim;
  data = new Tprec[dim];
  for(int i = 0; i < dim; ++i)
    data[i] = otro.data[i];
}
//
//--- DESTRUCTOR
//

template<typename Tprec>
Vector<Tprec>::~Vector()
{
  if ( dim > 0 )
    delete [] data;
}
//
//--- FUNCIONES MIEMBRO DE LA CLASE
//

template<typename Tprec>
void Vector<Tprec>::resize(int n)
{
  if ( dim > 0 )
    delete [] data;
  dim = n;
  data = new Tprec[dim];
}

template<typename Tprec>
Tprec Vector<Tprec>::calcNormL2() const
{
  Tprec norm = 0;
  for(int i = 0; i < dim; i++)
    norm += data[i] * data[i];
  return sqrt(norm);
}


template<typename Tprec>
Tprec Vector<Tprec>::calcNormMax() const 
{
  Tprec norm = 0;
  for (int i = 0; i < dim; i++) 
    if (norm < std::abs(data[i]) ) 
      norm = std::abs(data[i]);
  return norm;
}


template<typename Tprec>
Vector<Tprec>& Vector<Tprec>::operator=(const Vector<Tprec>& otro){
  if (this != &otro) {
    delete [] data;
    dim = otro.dim;
    data = new Tprec[dim];
    for(int i = 0; i < dim; ++i)
      data[i] = otro.data[i];
  }
  return *this;
}


template<typename Tprec>
Vector<Tprec>& Vector<Tprec>::operator+=(const Vector<Tprec>& v) 
{
  if (dim != v.dim ) {
    std::cout << "\n bad vector sizes \n";
    exit(1);
  }
  for (int i = 0; i < dim; i++) 
    data[i] += v(i);

  return *this;
}


template<typename Tprec>
Vector<Tprec>& Vector<Tprec>::operator-=(const Vector<Tprec>& v) 
{
  if (dim != v.dim ) {
    std::cout << "\n bad vector sizes \n";
  }
  for (int i = 0; i < dim; i++) 
    data[i] -= v(i);
  
  return *this;
}
//
//--- FUNCIONES AMIGAS DE LA CLASE
//

template<typename Tprec>
Vector<Tprec> operator+(const Vector<Tprec>& v1, const Vector<Tprec>& v2) 
{
  if (v1.dim != v2.dim ) {
    std::cout << "\n bad vector sizes \n";
    exit(1);
  }          
  Vector<Tprec> sum = v1; 
  sum += v2;
  return sum;
}

template<typename Tprec>
Vector<Tprec> operator-(const Vector<Tprec>& v1, const Vector<Tprec>& v2) 
{ 
  if (v1.dim != v2.dim ) {
    std::cout << "\n bad vector sizes \n";
    exit(1);
  } 
  Vector<Tprec> sum = v1; 
  sum -= v2;
  return sum;
}


template<typename Tprec>
Tprec dot(const Vector<Tprec>& v1, const Vector<Tprec>& v2)
{
  if (v1.dim != v2.dim ) {
    std::cout << "\n bad vector sizes \n";
    exit(1);
  } 
  Tprec dot = 0;
  for(int i = 0; i < v1.dim; i++)
    dot += v1(i) * v2(i);

  return dot;
}


template<typename Tprec>
Tprec operator*(const Vector<Tprec>& v1, const Vector<Tprec>& v2)
{
  if (v1.dim != v2.dim ) {
    std::cout << "\n bad vector sizes \n";
    exit(1);
  } 
  Tprec dot = 0;
  for(int i = 0; i < v1.dim; i++)
    dot += v1(i) * v2(i);

  return dot;
}

template<typename Tprec>
Vector<Tprec> operator*(const Tprec scalar, const Vector<Tprec>& v) 
{
  int n = v.getDim();
  Vector<Tprec> tm(n);
  for (int i = 0; i < n; i++) 
    tm(i) = scalar * v(i); 
  
  return tm;
}


template<typename Tprec>
Vector<Tprec> operator*(const Vector<Tprec>& v, const Tprec scalar) 
{
  return scalar * v; 
}



template<typename Tprec>
std::ostream& operator<<(std::ostream& s, const Vector<Tprec>& v)
{
  int dim = v.getDim();
  for(int i = 0; i < v.getDim(); ++i) {
      s << v(i) << std::endl;
  }
  return s;
}


