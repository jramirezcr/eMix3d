#ifndef _Array_HPP
#define _Array_HPP

#include "Eigen/Dense"

template<typename Tprec>
class Array2d{

    Tprec *data;
    int nx;
    int ny;

public:
    Array2d();
    Array2d(int, int);
    Array2d(const Array2d&);
    ~Array2d(){delete [] data;}

    void resize(int, int); // Dynamic resize Array
    int getDim(int) const; // Get dimension x(1) y(2) 

    //Overloaded 
    Tprec& operator()(int, int) const;


};

template<typename Tprec>
class Field1d{

    Tprec *data;
    int nx;

public:
    Field1d();
    Field1d(int);
    Field1d(const Field1d&);
    ~Field1d(){delete [] data;}

    virtual void resize(int); // Dynamic resize Array
    virtual int getDim(int) const; // Get dimension x(1)

    //Overloaded 
    virtual Tprec& operator()(int) const;
};


#include "Mesh/array.cpp"

#endif
