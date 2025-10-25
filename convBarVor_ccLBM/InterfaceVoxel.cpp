#include "InterfaceVoxel.h"
#include "Vector3D.H"
#include <vector>
#include <map>
#include <string>


void InterfaceVoxel::explode(bool& alternatingL1, bool& alternatingL0)
{
	//!loop over all discrete velocity directions
	for (int dir = 0; dir < 19; dir++) //!explode pops in all lattice directions since coalesce() overrides unused pops anyway...
	{
		double explodedState = distr_[alternatingL0][dir];
		for (int par = 0; par < 8; par++) this->getPartner(par)->setDistribution(alternatingL1, dir, explodedState);
	}
}

void InterfaceVoxel::coalesce(bool& alternatingL1, bool& alternatingL0, string& interface)
{
	vector<int> directionsBottom = { 4, 12, 8, 10, 6 };
	vector<int> directionsTop = { 5, 13, 9, 7, 11 };
	vector<int> directionsLeft = { 1, 16, 8, 14, 9 };
	vector<int> directionsRight = { 3, 15, 12, 17, 13 };
	vector<int> directionsInlet = { 0, 6,  7,  15, 16 };
	vector<int> directionsOutlet = { 2, 10, 11, 14, 17 };

	vector<int> directionsBottomLeft = { 4, 1, 16, 12, 8, 14, 9, 10, 6 };
	vector<int> directionsBottomRight = { 4, 3, 15, 12, 8, 17, 13, 10, 6 };
	vector<int> directionsTopLeft = { 5, 1, 16, 8, 14, 13, 9, 7, 11 };
	vector<int> directionsTopRight = { 5, 3, 15, 12, 17, 13, 9, 7, 11 };
	vector<int> directionsInletBottom = { 0, 6, 7, 15, 16, 4, 12, 8, 10 };
	vector<int> directionsInletTop = { 0, 6,  7,  15, 16, 5, 13, 9, 11 };
	vector<int> directionsInletLeft = { 0, 6, 7, 15, 1, 16, 8, 14, 9 };
	vector<int> directionsInletRight = { 0, 6, 7, 16, 3, 15, 12, 17, 13 };
	vector<int> directionsOutletLeft = { 2, 10, 11, 14, 17, 1, 16, 8, 9 };
	vector<int> directionsOutletRight = { 2, 10, 11, 14, 17, 3, 15, 12, 13 };
	vector<int> directionsOutletBottom = { 2, 10, 11, 14, 17, 4, 12, 8, 6 };
	vector<int> directionsOutletTop = { 2, 10, 11, 14, 17, 5, 13, 9, 7 };

	vector<int> directionsInletBottomLeft = { 0, 6, 7, 15, 16, 4, 12, 8, 10, 1, 14, 9 };
	vector<int> directionsInletBottomRight = { 0, 6, 7, 15, 16, 4, 12, 8, 10, 3, 17, 13 };
	vector<int> directionsInletTopLeft = { 0, 6, 7, 15, 16, 5, 13, 9, 11, 1, 8, 14 };
	vector<int> directionsInletTopRight = { 0, 6, 7, 15, 16, 5, 13, 9, 11, 3, 12, 17 };
	vector<int> directionsOutletBottomLeft = { 2, 10, 11, 14, 17, 4, 12, 8, 6, 1, 16, 9 };
	vector<int> directionsOutletBottomRight = { 2, 10, 11, 14, 17, 4, 12, 8, 6, 3, 15, 13 };
	vector<int> directionsOutletTopLeft = { 2, 10, 11, 14, 17, 5, 13, 7, 1, 16, 8, 9 };
	vector<int> directionsOutletTopRight = { 2, 10, 11, 14, 5, 13, 9, 7, 3, 15, 12, 17 };

	//!create a map of interface strings and corresponding directions
	map<string, vector<int>> dirMap{ {"interfaceBottom", directionsBottom},  {"interfaceTop", directionsTop},  {"interfaceLeft", directionsLeft},
									 {"interfaceRight", directionsRight},  {"interfaceInlet", directionsInlet},  {"interfaceOutlet", directionsOutlet},
									 {"interfaceBottomLeft", directionsBottomLeft},  {"interfaceBottomRight", directionsBottomRight},  {"interfaceTopLeft", directionsTopLeft},
									 {"interfaceTopRight", directionsTopRight},  {"interfaceInletBottom", directionsInletBottom},  {"interfaceInletTop", directionsInletTop},
									 {"interfaceInletLeft", directionsInletLeft},  {"interfaceInletRight", directionsInletRight},  {"interfaceOutletLeft", directionsOutletLeft},
									 {"interfaceOutletRight", directionsOutletRight},  {"interfaceOutletBottom", directionsOutletBottom},  {"interfaceOutletTop", directionsOutletTop},
									 {"interfaceInletBottomLeft", directionsInletBottomLeft},  {"interfaceInletBottomRight", directionsInletBottomRight},  {"interfaceInletTopLeft", directionsInletTopLeft},
									 {"interfaceInletTopRight", directionsInletTopRight},  {"interfaceOutletBottomLeft", directionsOutletBottomLeft},  {"interfaceOutletBottomRight", directionsOutletBottomRight},
									 {"interfaceOutletTopLeft", directionsOutletTopLeft},  {"interfaceOutletTopRight", directionsOutletTopRight} };

	for (int dir = 0; dir <= dirMap[interface].size() - 1; dir++)
	{
		double coalescedState = 0.;
		for (int par = 0; par < 8; par++)
		{
			coalescedState += this->getPartner(par)->getDistribution(alternatingL1, dirMap[interface][dir]);
		}
		coalescedState /= 8.;
		distr_[alternatingL0][dirMap[interface][dir]] = coalescedState;
	}
}

void InterfaceVoxel::calcGhostDistribution(bool& alternatingL1, double* ghostArr)
{
	for (int dir = 0; dir < 19; dir++)
	{
		double coalescedState = 0.;
		for (int par = 0; par < 8; par++)
		{
			coalescedState += this->getPartner(par)->getDistribution(alternatingL1, dir);
		}
		coalescedState /= 8.;
		ghostArr[dir] = coalescedState;
	}
}


void InterfaceVoxel::explodeLinear(bool& alternatingL1, bool& alternatingL0, Control& bc, string& interface)
{
	//!spacing on coarse level
	double spacing = bc.getSpacing() / pow(2., level_);

	double* xsiX = bc.getXsiX();
	double* xsiY = bc.getXsiY();
	double* xsiZ = bc.getXsiZ();
	double xsiSquared;

	double posVectorCoarse[3] = { 0. };
	posVectorCoarse[0] = this->getX();
	posVectorCoarse[1] = this->getY();
	posVectorCoarse[2] = this->getZ();

	double distrTemp1 = 0., distrTemp2 = 0., distrTemp3 = 0., distrTemp4 = 0.;
	double distrGradient[3] = { 0. };
	double posVectorFine[3] = { 0. };
	double posVectorDiff[3] = { 0. };
	double dotProductXsi = 0., dotProductPos = 0.;
	double explodedState = 0.;

	if ((interface == "interfaceBottom") || (interface == "interfaceTop"))
	{
		for (int dir = 0; dir < 19; dir++)
		{
			//!post-collision distributions in direction dir of coarse interface neighbors for linear interpolation in explode procedure
			distrTemp1 = this->getNeighbor(0)->getDistribution(alternatingL0, dir);
			distrTemp2 = this->getNeighbor(2)->getDistribution(alternatingL0, dir);
			distrTemp3 = this->getNeighbor(1)->getDistribution(alternatingL0, dir);
			distrTemp4 = this->getNeighbor(3)->getDistribution(alternatingL0, dir);

			distrGradient[0] = (distrTemp1 - distrTemp2) / (2. * spacing);
			distrGradient[1] = (distrTemp3 - distrTemp4) / (2. * spacing);
			distrGradient[2] = 0.;

			dotProductXsi = xsiX[dir] * distrGradient[0] + xsiY[dir] * distrGradient[1] + xsiZ[dir] * distrGradient[2];
			xsiSquared = pow(xsiX[dir], 2.) + pow(xsiY[dir], 2.) + pow(xsiZ[dir], 2.);

			if (dir != 18)
			{
				distrGradient[0] -= (xsiX[dir] * dotProductXsi) / xsiSquared;
				distrGradient[1] -= (xsiY[dir] * dotProductXsi) / xsiSquared;
				distrGradient[2] -= (xsiZ[dir] * dotProductXsi) / xsiSquared;
			}
			else {}

			for (int par = 0; par < 8; par++)
			{
				//!underlying fine voxels
				posVectorFine[0] = this->getPartner(par)->getX();
				posVectorFine[1] = this->getPartner(par)->getY();
				posVectorFine[2] = this->getPartner(par)->getZ();

				posVectorDiff[0] = posVectorFine[0] - posVectorCoarse[0];
				posVectorDiff[1] = posVectorFine[1] - posVectorCoarse[1];
				posVectorDiff[2] = posVectorFine[2] - posVectorCoarse[2];

				dotProductPos = posVectorDiff[0] * distrGradient[0] + posVectorDiff[1] * distrGradient[1] + posVectorDiff[2] * distrGradient[2];

				explodedState = distr_[alternatingL0][dir] + dotProductPos;
				this->getPartner(par)->setDistribution(alternatingL1, dir, explodedState);
			}
		}
	}

	else if ((interface == "interfaceLeft") || (interface == "interfaceRight"))
	{
		for (int dir = 0; dir < 19; dir++)
		{
			//!post-collision distributions in direction dir of coarse interface neighbors for linear interpolation in explode procedure
			distrTemp1 = this->getNeighbor(0)->getDistribution(alternatingL0, dir);
			distrTemp2 = this->getNeighbor(2)->getDistribution(alternatingL0, dir);
			distrTemp3 = this->getNeighbor(4)->getDistribution(alternatingL0, dir);
			distrTemp4 = this->getNeighbor(5)->getDistribution(alternatingL0, dir);

			distrGradient[0] = (distrTemp1 - distrTemp2) / (2. * spacing);
			distrGradient[1] = 0.;
			distrGradient[2] = (distrTemp3 - distrTemp4) / (2. * spacing);

			dotProductXsi = xsiX[dir] * distrGradient[0] + xsiY[dir] * distrGradient[1] + xsiZ[dir] * distrGradient[2];
			xsiSquared = pow(xsiX[dir], 2.) + pow(xsiY[dir], 2.) + pow(xsiZ[dir], 2.);

			if (dir != 18)
			{
				distrGradient[0] -= (xsiX[dir] * dotProductXsi) / xsiSquared;
				distrGradient[1] -= (xsiY[dir] * dotProductXsi) / xsiSquared;
				distrGradient[2] -= (xsiZ[dir] * dotProductXsi) / xsiSquared;
			}
			else {}

			for (int par = 0; par < 8; par++)
			{
				//!underlying fine voxels
				posVectorFine[0] = this->getPartner(par)->getX();
				posVectorFine[1] = this->getPartner(par)->getY();
				posVectorFine[2] = this->getPartner(par)->getZ();

				posVectorDiff[0] = posVectorFine[0] - posVectorCoarse[0];
				posVectorDiff[1] = posVectorFine[1] - posVectorCoarse[1];
				posVectorDiff[2] = posVectorFine[2] - posVectorCoarse[2];

				dotProductPos = posVectorDiff[0] * distrGradient[0] + posVectorDiff[1] * distrGradient[1] + posVectorDiff[2] * distrGradient[2];

				explodedState = distr_[alternatingL0][dir] + dotProductPos;
				this->getPartner(par)->setDistribution(alternatingL1, dir, explodedState);
			}
		}
	}

	else if ((interface == "interfaceInlet") || (interface == "interfaceOutlet"))
	{
		for (int dir = 0; dir < 19; dir++)
		{
			//!post-collision distributions in direction dir of coarse interface neighbors for linear interpolation in explode procedure
			distrTemp1 = this->getNeighbor(1)->getDistribution(alternatingL0, dir);
			distrTemp2 = this->getNeighbor(3)->getDistribution(alternatingL0, dir);
			distrTemp3 = this->getNeighbor(4)->getDistribution(alternatingL0, dir);
			distrTemp4 = this->getNeighbor(5)->getDistribution(alternatingL0, dir);

			distrGradient[0] = 0.;
			distrGradient[1] = (distrTemp1 - distrTemp2) / (2. * spacing);
			distrGradient[2] = (distrTemp3 - distrTemp4) / (2. * spacing);

			dotProductXsi = xsiX[dir] * distrGradient[0] + xsiY[dir] * distrGradient[1] + xsiZ[dir] * distrGradient[2];
			xsiSquared = pow(xsiX[dir], 2.) + pow(xsiY[dir], 2.) + pow(xsiZ[dir], 2.);

			if (dir != 18)
			{
				distrGradient[0] -= (xsiX[dir] * dotProductXsi) / xsiSquared;
				distrGradient[1] -= (xsiY[dir] * dotProductXsi) / xsiSquared;
				distrGradient[2] -= (xsiZ[dir] * dotProductXsi) / xsiSquared;
			}
			else {}

			for (int par = 0; par < 8; par++)
			{
				//!underlying fine voxels
				posVectorFine[0] = this->getPartner(par)->getX();
				posVectorFine[1] = this->getPartner(par)->getY();
				posVectorFine[2] = this->getPartner(par)->getZ();

				posVectorDiff[0] = posVectorFine[0] - posVectorCoarse[0];
				posVectorDiff[1] = posVectorFine[1] - posVectorCoarse[1];
				posVectorDiff[2] = posVectorFine[2] - posVectorCoarse[2];

				dotProductPos = posVectorDiff[0] * distrGradient[0] + posVectorDiff[1] * distrGradient[1] + posVectorDiff[2] * distrGradient[2];

				explodedState = distr_[alternatingL0][dir] + dotProductPos;
				this->getPartner(par)->setDistribution(alternatingL1, dir, explodedState);
			}
		}
	}
}

