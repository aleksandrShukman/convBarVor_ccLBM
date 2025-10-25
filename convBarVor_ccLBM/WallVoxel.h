#pragma once

#include "Voxel.h"

using namespace std;

class WallVoxel: public Voxel
{

public:

	//!constructor
	WallVoxel() : Voxel(){ tag_ = "wall"; } //!standard

	WallVoxel(double x, double y, double z, long long indexI, long long indexJ, long long indexK, long long level) : Voxel(x, y, z, indexI, indexJ, indexK, level)
	{
		tag_ = "wall";
	}

	//!destructor
	virtual ~WallVoxel() 
	{

	}

	virtual inline void halfWayBounceBack(bool& distr, Control& bc)
	{
		for (int dir = 0; dir < 18; dir++) 
		{
			if (neighbors_[dir] != NULL)
			{
				//!calc opposite direction for bounce back
				int counterDir = 18;
				switch (dir)
				{
					case 0:  counterDir = 2;  break; //!east to west
					case 1:  counterDir = 3;  break; //!north to south
					case 2:  counterDir = 0;  break; //!west to east
					case 3:  counterDir = 1;  break; //!south to north
					case 4:  counterDir = 5;  break; //!top to bottom
					case 5:  counterDir = 4;  break; //!bottom to top
					case 6:  counterDir = 11; break; //!top east to bottom west
					case 7:  counterDir = 10; break; //!bottom east to top west
					case 8:  counterDir = 13; break; //!top north to bottom south
					case 9:  counterDir = 12; break; //!bottom north to top south
					case 10: counterDir = 7;  break; //!top west to bottom east
					case 11: counterDir = 6;  break; //!bottom west to top east
					case 12: counterDir = 9;  break; //!top south to bottom north
					case 13: counterDir = 8;  break; //!bottom south to top north
					case 14: counterDir = 15; break; //!north west to south east
					case 15: counterDir = 14; break; //!south east to north west
					case 16: counterDir = 17; break; //!north east to south west
					case 17: counterDir = 16; break; //!south west to north east
				}
				//!perform bounce-back
				neighbors_[dir]->setDistribution(distr, dir, distr_[distr][counterDir]);
			}
		}
	} //!end bounceBack

	virtual inline void bouzidiBounceBack()
	{

	}

private:

};