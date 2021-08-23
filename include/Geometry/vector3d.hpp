/************************************************
Written by Jorge Ramirez Cruz
fga_jorgeram@hotmail.com
************************************************/

#ifndef _VECTOR3D_H
#define _VECTOR3D_H

//External Includes
//#include<iostream>
#include<math.h>


namespace vecMath{

template<typename real>
class Vector3d
{

public:

    real x; //Value X-Axis
    real y; //Value Y-Axis
    real z; //Value Z-Axis

private:

    real garb;

public:

    //Default Constructor
    Vector3d(): x(0.0), y(0.0), z(0.0) {}

    //Explicit Constructor
    Vector3d(const real x_garb, const real y_garb, const real z_garb)
    : x(x_garb), y(y_garb), z(z_garb){}

    //Change the direction
    void invert()
    {
        x = -x;
        y = -y;
        z = -z;
    }

    //Magnitude of vector
    real magnitude() const
    {
    return sqrt(x*x + y*y + z*z);
    }

    //Magnitude squared of vector
    real squareMagnitude() const
    {
    return x*x + y*y + z*z;
    }


    //normalize the vector, magnitude = 1.0
    void normalize()
    {
        real l = magnitude();
        if (l > 0.0)
        {
           x *= ((real)1)/l;
           y *= ((real)1)/l;
           z *= ((real)1)/l;
        }
    }



    // Multiplies OPETATIONS
    //Multiplies the vector by a scalar
    void operator *= (const real value)
    {
        x *= value; 
        y *= value; 
        z *= value; 
    }

    //Return a COPY of the vector multiply by a scalar
    Vector3d operator * ( const real value)
    {
    	return Vector3d(value*x, value*y, value*z);
    }

    // ADD OPETATIONS
    //Overload vector sum add the original  vector
    void operator += (const Vector3d& v)
    {
       x += v.x;
       y += v.y;
       z += v.z;
    }

    //Overload vector sum add a copy of the vector
    Vector3d operator + (const Vector3d& v)
    {
         return Vector3d(x + v.x, y + v.y, z + v.z);
    }

    // SUBSTRACT OPETATIONS
    //Overload vector sum add the original  vector
    void operator -= (const Vector3d& v)
    {
       x -= v.x;
       y -= v.y;
       z -= v.z;
    }

    //Overload vector sum add a copy of the vector
    Vector3d operator - (const Vector3d& v)
    {
    return Vector3d(x - v.x, y - v.y, z - v.z);
    }


    //Saxpy A = A + c*B
    void addScaledVector(const Vector3d & v, const real scalar)
    { 
        x += v.x * scalar;
        y += v.y * scalar;
        z += v.z * scalar;
    }

    //Vector multiply
    Vector3d componentProduct(const Vector3d &vector) const
    {
        return Vector3d(x * vector.x, y * vector.y, z * vector.z);
    }

    void componentProductUpdate(const Vector3d &vector)
    {
        x *= vector.x;
        y *= vector.y;
        z *= vector.z;
    }

    //Scalar dot product
    real scalarProduct(const Vector3d &vector) const
    {
        return x*vector.x + y*vector.y + z*vector.z;
    }

    real operator *(const Vector3d &vector) const
    {
    return x*vector.x + y*vector.y + z*vector.z;
    }

    //Return the Vector Product
    Vector3d vectorProduct(const Vector3d &vector) const
    {
       return Vector3d( y*vector.z-z*vector.y ,
                        z*vector.x-x*vector.z ,
                        x*vector.y-y*vector.x);
    }

    void operator %= (const Vector3d &vector)
    {
        *this = vectorProduct(vector);
    }

    Vector3d operator % (const Vector3d &vector) const
    {
        return Vector3d( y*vector.z-z*vector.y,
                         z*vector.x-x*vector.z,
                         x*vector.y-y*vector.x );
}
};
//End class Vector3d


}

#endif
