#pragma once

#include "Voxel.h"
#include "FluidVoxel.h"
#include "BoundaryVoxel.h"
#include "WallVoxel.h"
#include "Control.h"
#include "InterfaceVoxel.h"
#include "Node.h"
#include <omp.h>
#include <string>
#include <fstream>
#include <iomanip>

using namespace std;

class Lattice
{

public:

	//!constructors

	Lattice(Control& bc) : bc_(bc) //!copy reference
	{
		//!create lattice data structure
		int maxLevel = bc.getLevelMax(); //!maximum level
		int numberOfLevels = maxLevel + 1; //!number of levels

		lattice_ = new Voxel****[numberOfLevels]; //!create structure for every refinement level

		for (int level = 0; level <= maxLevel; level++)
		{
			int imax = bc.getNumberOfVoxelsX() * pow(2., level) - 1; //!index of last voxel in x direction in ->this level
			int jmax = bc.getNumberOfVoxelsY() * pow(2., level) - 1; //!index of last voxel in y direction in ->this level
			int kmax = bc.getNumberOfVoxelsZ() * pow(2., level) - 1; //!index of last voxel in z direction in ->this level

			int numberOfVoxelsX = imax + 1; //!total number of voxels in x direction in ->this level
			int numberOfVoxelsY = jmax + 1; //!total number of voxels in y direction in ->this level
			int numberOfVoxelsZ = kmax + 1; //!total number of voxels in z direction in ->this level

			lattice_[level] = new Voxel***[numberOfVoxelsX]; //!create arrays in x direction for every level

			for(int i = 0; i <= imax; i++)
			{
				lattice_[level][i] = new Voxel**[numberOfVoxelsY]; //!create arrays in y direction for every level

				for(int j = 0; j <= jmax; j++)
				{
					lattice_[level][i][j] = new Voxel*[numberOfVoxelsZ]; //!create arrays in z direction for every level
				} //!end j
			} //!end i
		} //!end level

		//!create mesh data structure
		mesh_ = new Node**** [numberOfLevels]; //!create structure for every refinement level

		for (int level = 0; level <= maxLevel; level++)
		{
			int imax = (bc.getNumberOfNodesX() - 1) * pow(2., level); //!index of last node in x direction in ->this level
			int jmax = (bc.getNumberOfNodesY() - 1) * pow(2., level); //!index of last node in y direction in ->this level
			int kmax = (bc.getNumberOfNodesZ() - 1) * pow(2., level); //!index of last node in z direction in ->this level

			int numberOfNodesX = imax + 1; //!total number of nodes in x direction in ->this level
			int numberOfNodesY = jmax + 1; //!total number of nodes in y direction in ->this level
			int numberOfNodesZ = kmax + 1; //!total number of nodes in z direction in ->this level

			mesh_[level] = new Node*** [numberOfNodesX]; //!create arrays in x direction for every level

			for (int i = 0; i <= imax; i++)
			{
				mesh_[level][i] = new Node** [numberOfNodesY]; //!create arrays in y direction for every level

				for (int j = 0; j <= jmax; j++)
				{
					mesh_[level][i][j] = new Node* [numberOfNodesZ]; //!create arrays in z direction for every level
				} //!end j
			} //!end i
		} //!end level

		createLattice(); //!create lattice
		createMesh(); //!create mesh
		setConnectivity(); //!define connectivity

	} //!end constructor


	//!destructor
	~Lattice() //!release all dynamically allocated memory
	{
		if (lattice_ != NULL)
		{
			int maxLevel = bc_.getLevelMax(); //!maximum level
			int numberOfLevels = maxLevel + 1; //!number of levels

			for (int level = 0; level <= maxLevel; level++)
			{
				if (lattice_[level])
				{
					int imax = (bc_.getNumberOfVoxelsX() - 1) * pow(2., level); //!index of last node in x direction in ->this level
					int jmax = (bc_.getNumberOfVoxelsY() - 1) * pow(2., level); //!index of last node in y direction in ->this level
					int kmax = (bc_.getNumberOfVoxelsZ() - 1) * pow(2., level); //!index of last node in z direction in ->this level

					for(int i = 0; i <= imax; i++)
					{
						if (lattice_[level][i])
						{
							for(int j = 0; j <= jmax; j++)
							{
								if (lattice_[level][i][j])
								{
									for(int k = 0; k <= kmax; k++)
									{
										if (lattice_[level][i][j][k]) delete lattice_[level][i][j][k];
									}
									delete[] lattice_[level][i][j]; //!delete arrays in z direction for every level
								}
							} //!end j
							delete[] lattice_[level][i]; //!delete arrays in y direction for every level
						}
					} //!end i
					delete[] lattice_[level]; //!delete arrays in x direction for every level
				}
			} //!end level
			delete[] lattice_; //!delete levels
		}

		if (mesh_ != NULL)
		{
			int maxLevel = bc_.getLevelMax(); //!maximum level
			int numberOfLevels = maxLevel + 1; //!number of levels

			for (int level = 0; level <= maxLevel; level++)
			{
				if (mesh_[level])
				{
					int imax = (bc_.getNumberOfNodesX() - 1) * pow(2., level); //!index of last node in x direction in ->this level
					int jmax = (bc_.getNumberOfNodesY() - 1) * pow(2., level); //!index of last node in y direction in ->this level
					int kmax = (bc_.getNumberOfNodesZ() - 1) * pow(2., level); //!index of last node in z direction in ->this level

					for (int i = 0; i <= imax; i++)
					{
						if (mesh_[level][i])
						{
							for (int j = 0; j <= jmax; j++)
							{
								if (mesh_[level][i][j])
								{
									for (int k = 0; k <= kmax; k++)
									{
										if (mesh_[level][i][j][k]) delete mesh_[level][i][j][k];
									}
									delete[] mesh_[level][i][j]; //!delete arrays in z direction for every level
								}
							} //!end j
							delete[] mesh_[level][i]; //!delete arrays in y direction for every level
						}
					} //!end i
					delete[] mesh_[level]; //!delete arrays in x direction for every level
				}
			} //!end level
			delete[] mesh_; //!delete levels
		}
	} //!end denstructor

	//!member
	virtual void createLattice();
	virtual void createMesh();
	virtual void setConnectivity();
	virtual void solve(bool &alternating);
	virtual void solve(bool& alternatingL0, bool& alternatingL1);
	virtual void timeSteppingSRT(int& level,bool &alternating);
	virtual void timeSteppingMRT_GSbasis(int& level, bool& alternating);
	virtual void timeSteppingMRT_RMbasis(int& level, bool& alternating);
	virtual void timeSteppingRR(int& level, bool& alternating);
	virtual void timeSteppingHRR(int& level, bool& alternating);
	virtual void nestedTimeSteppingSRT(bool &alternatingL1,bool &alternatingL0, Control &bc_);
	virtual void nestedTimeSteppingMRT_GSbasis(bool& alternatingL1, bool& alternatingL0, Control& bc_);
	virtual void nestedTimeSteppingMRT_RMbasis(bool& alternatingL1, bool& alternatingL0, Control& bc_);
	virtual void nestedTimeSteppingRR(bool& alternatingL1, bool& alternatingL0, Control& bc_);
	virtual void nestedTimeSteppingHRR(bool& alternatingL1, bool& alternatingL0, Control& bc_);
	virtual void writeResultsVTK(string filename, bool &alternating, bool& alternatingL1);
	virtual void writeResultsAtVoxel(string filename, bool& alternatingL0, bool& alternatingL1, Control& bc_);
	virtual void initialization(Control& bc);

	virtual inline double getDeviation() { return deviation_; }
	virtual inline double getDistribution(bool& alternating, int& dir, int& level, int& i, int& j, int& k) { if (dir < 19) return lattice_[level][i][j][k]->getDistribution(alternating, dir); else return NULL; };
	virtual double readValue(string filename, string keyword); //!filename: name of initialization file; keyword: string to be read from file; return value: (double) to corresponding keyword

private:
	Control &bc_; //!reference to instance of control class
	Voxel***** lattice_; //!pointer for creating lattice data structure consisting of voxels
	Node***** mesh_; //!pointer for creating mesh data structure consisting of nodes
};
