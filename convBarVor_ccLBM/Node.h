#pragma once

#include "Control.h"

using namespace std;

class Node
{
public:
	//!constructors
	Node() {} //!standard

	Node(double x, double y, double z, int indexI, int indexJ, int indexK, int level) : x_(x), y_(y), z_(z), indexI_(indexI), indexJ_(indexJ), indexK_(indexK), level_(level) //!init all data
	{

	}

	//!destructor
	virtual ~Node()
	{	

	}

	//!get and set methods
	virtual inline int    getIndexI() const { return indexI_; }
	virtual inline int    getIndexJ() const { return indexJ_; }
	virtual inline int    getIndexK() const { return indexK_; }
	virtual inline int    getLevel() const { return level_; }
	virtual inline double getX() const { return x_; }
	virtual inline double getY() const { return y_; }
	virtual inline double getZ() const { return z_; }


private:

	//!coordinates
	double x_;
	double y_;
	double z_;

protected:

	//!level
	int level_;

	//!node indices
	int indexI_;
	int indexJ_;
	int indexK_;

};