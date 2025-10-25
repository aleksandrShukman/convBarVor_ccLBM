#ifndef VECTOR3D_H
#define VECTOR3D_H

#ifndef ASSERT
#define ASSERT

#include "assert.h"
#include <math.h>

class Vector3D
{
public:

	// Constructor
	Vector3D()
	{
		coordinates_[0]=0;
		coordinates_[1]=0;
		coordinates_[2]=0;
	}

	Vector3D(double x, double y, double z)
	{
		coordinates_[0]=x;
		coordinates_[1]=y;
		coordinates_[2]=z;
	}

	Vector3D(const Vector3D& vec)
	{
			*this=vec;
	}
	// Assignment operator
	Vector3D& operator=(const Vector3D& rhs)
	{
		coordinates_[0]=rhs.getX();
		coordinates_[1]=rhs.getY();
		coordinates_[2]=rhs.getZ();
		return *this;
	}

	// Members

		// Get
		inline double getX() const { return coordinates_[0]; }
		inline double getY() const { return coordinates_[1]; }
		inline double getZ() const { return coordinates_[2]; }
		inline const double* getComponentArray() const { return coordinates_; }

		inline double getComponent(unsigned short i)
		{
			ASSERT(i<3);
			return coordinates_[i];
		}

		// Set
		inline void setVector3D(double x,double y,double z) { coordinates_[0]=x; coordinates_[1]=y; coordinates_[2]=z; }
		inline void setComponent(unsigned char comp,double value) { coordinates_[comp]=value; }

		// Operators
		inline double operator[] (unsigned short index)
		{
			ASSERT(index<3);
			return coordinates_[index];
		}

		inline Vector3D operator+=(const Vector3D& n)			// Adds a vector to a vector component by component (with assignment)
		{
			const double* temp=n.getComponentArray();
			for (unsigned short i=0; i<3;i++) coordinates_[i]+=temp[i];
			return (*this);
		}

		inline Vector3D operator+(const Vector3D& n1) const		// Adds a vector to a vector component by component
		{
			Vector3D result(coordinates_[0]+n1.getX(),coordinates_[1]+n1.getY(),coordinates_[2]+n1.getZ());
			return result;
		}

		inline Vector3D operator-=(const Vector3D& n)			// Substracts a vector to a vector component by component (with assignment)
		{
			const double* temp=n.getComponentArray();
			for (unsigned short i=0; i<3;i++) coordinates_[i]-=temp[i];
			return (*this);
		}

		inline Vector3D operator-(const Vector3D& n1) const		// Substracts a vector to a vector component by component
		{
			Vector3D result;
			result.setVector3D(this->getX()-n1.getX(),this->getY()-n1.getY(),this->getZ()-n1.getZ());
			return result;
		}

		inline Vector3D operator*=(const double& s)				// Multiplies a vector with a scalar (by components) (with assignment)
		{
			for (unsigned short i=0; i<3;i++) coordinates_[i]*=s;
			return (*this);
		}

		inline Vector3D operator*(const double& s)				// Multiplies a vector with a scalar (by components)
		{
			Vector3D result(coordinates_[0]*s,coordinates_[1]*s,coordinates_[2]*s);
			return result;
		}

		inline Vector3D operator/=(const double& s)				// Devides a vector by a scalar (by components) (with assignment)
		{
			for (unsigned short i=0; i<3;i++) coordinates_[i]/=s;
			return (*this);
		}

		inline Vector3D operator/(const double& s)				// Devides a vector by a scalar (by components)
		{
			Vector3D result(coordinates_[0]/s,coordinates_[1]/s,coordinates_[2]/s);
			return result;
		}

		inline double absoluteValue() const					// Calcs absolute value of a vector
		{
			double val=sqrt(coordinates_[0]*coordinates_[0]+coordinates_[1]*coordinates_[1]+coordinates_[2]*coordinates_[2]);
			return val;
		}

		inline Vector3D crossProduct(Vector3D &n1, Vector3D &n2)	// Calcs cross product of two vectors
		{
			coordinates_[0]=n1.getY()*n2.getZ()-n1.getZ()*n2.getY();
			coordinates_[1]=n1.getZ()*n2.getX()-n1.getX()*n2.getZ();
			coordinates_[2]=n1.getX()*n2.getY()-n1.getY()*n2.getX();
			return (*this);
		}

		inline double dotProduct(Vector3D &n1, Vector3D &n2)			//Calculates the scalar product of two vectors
		{
			return (n1.getX()*n2.getX()+n1.getY()*n2.getY()+n1.getZ()*n2.getZ());
		}


private:

	double coordinates_[3];

};

#endif
#endif
