#pragma once

#include "Voxel.h"

using namespace std;

class InterfaceVoxel : public Voxel
{
public:

	//!constructor
	InterfaceVoxel() : Voxel() { tag_ = "gridcoupling"; }

	InterfaceVoxel(double x, double y, double z, int indexI, int indexJ, int indexK, int level) : Voxel(x, y, z, indexI, indexJ, indexK, level)
	{
		tag_ = "gridcoupling";
	}

	//!destructor
	virtual ~InterfaceVoxel()
	{

	}

	//!members
	virtual void explode(bool& alternatingL1, bool& alternatingL0);
	virtual void explodeLinear(bool& alternatingL1, bool& alternatingL0, Control& bc, string& interfaceType);
	virtual void coalesce(bool& alternatingL1, bool& alternatingL0, string& interfaceType);
	virtual void calcGhostDistribution(bool& alternatingL1, double* ghostArr);

	virtual inline void transport(bool& distr) //!transport step (push step)
	{
		for (int dir = 0; dir < 18; dir++)
		{
			if (neighbors_[dir] != NULL)
			{
				neighbors_[dir]->setDistribution(distr, dir, distr_[!distr][dir]);
			}
		}
		distr_[distr][18] = distr_[!distr][18];
	}

	virtual inline void transportPull(bool& distr) //!transport step (pull step)
	{
		int counterDir;
		bool negDistr = !distr;

		for (int dir = 0; dir < 18; dir++)
		{
			if (neighbors_[dir] != NULL)
			{
				switch (dir)
				{
				case 0:  counterDir = 2;  break;
				case 1:  counterDir = 3;  break;
				case 2:  counterDir = 0;  break;
				case 3:  counterDir = 1;  break;
				case 4:  counterDir = 5;  break;
				case 5:  counterDir = 4;  break;
				case 6:  counterDir = 11; break;
				case 7:  counterDir = 10; break;
				case 8:  counterDir = 13; break;
				case 9:  counterDir = 12; break;
				case 10: counterDir = 7;  break;
				case 11: counterDir = 6;  break;
				case 12: counterDir = 9;  break;
				case 13: counterDir = 8;  break;
				case 14: counterDir = 15; break;
				case 15: counterDir = 14; break;
				case 16: counterDir = 17; break;
				case 17: counterDir = 16; break;
				}
				distr_[distr][dir] = neighbors_[counterDir]->getDistribution(negDistr, dir);
			}
		}
		distr_[distr][18] = distr_[negDistr][18];
	}

private:


};