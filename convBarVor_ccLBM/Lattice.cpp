#define _USE_MATH_DEFINES

#include "Lattice.h"
#include "math.h"

void Lattice::createLattice()
{
	//!create voxels in lattice structure
	//!lattice topology: bottom and top layer (+/- z normal) of lattice are periodic boundaries
	//!left and right boundaries (+/- y normal) are periodic boundaries
	//!front and back boundaries (+/- x normal) are periodic boundaries
	//!remaining voxels are fluid

	int refLayer, refLayerX;

	//!create voxels in array
	for (int level = 0; level <= bc_.getLevelMax(); level++)
	{
		int imax = bc_.getNumberOfVoxelsX() * pow(2., level) - 1; //!index of last voxel in x direction in ->this level
		int jmax = bc_.getNumberOfVoxelsY() * pow(2., level) - 1; //!index of last voxel in y direction in ->this level
		int kmax = bc_.getNumberOfVoxelsZ() * pow(2., level) - 1; //!index of last voxel in z direction in ->this level

		refLayer = bc_.getRefinementLayer() * pow(2, level);
		refLayerX = bc_.getRefinementLayerX() * pow(2, level);

		for (int k = 0; k <= kmax; k++)
		{
			for (int j = 0; j <= jmax; j++)
			{
				for (int i = 0; i <= imax; i++)
				{
					double xCoord = ((2. * i + 1) / pow(2., level + 1)) * bc_.getSpacing();
					double yCoord = ((2. * j + 1) / pow(2., level + 1)) * bc_.getSpacing();
					double zCoord = ((2. * k + 1) / pow(2., level + 1)) * bc_.getSpacing();

					//!level 0 with refinement
					if ((level == 0) && (bc_.getLevelMax() > 0))
					{
						//!interface voxels front and rear interface
						if (((i == refLayerX - 1 || i == imax - refLayer + 1)) && ((j >= refLayer - 1) && (j <= jmax - refLayer + 1)) && ((k >= refLayer - 1) && (k <= kmax - refLayer + 1))) lattice_[level][i][j][k] = new InterfaceVoxel(xCoord, yCoord, zCoord, i, j, k, level); //!if voxel is at interface --> interface voxel

						//!gridcoupling nodes bottom and top
						else if (((k == refLayer - 1) || (k == kmax - refLayer + 1)) && ((j >= refLayer - 1) && (j <= jmax - refLayer + 1)) && ((i >= refLayerX - 1) && (i <= imax - refLayer + 1))) lattice_[level][i][j][k] = new InterfaceVoxel(xCoord, yCoord, zCoord, i, j, k, level); //!if voxel is on interface --> interface voxel

						//!gridcoupling nodes left and right
						else if (((j == refLayer - 1) || (j == jmax - refLayer + 1)) && ((i >= refLayerX - 1) && (i <= imax - refLayer + 1)) && ((k >= refLayer - 1) && (k <= kmax - refLayer + 1))) lattice_[level][i][j][k] = new InterfaceVoxel(xCoord, yCoord, zCoord, i, j, k, level); //!if voxel is on interface --> interface voxel

						//!standard fluid voxels
						else lattice_[level][i][j][k] = new FluidVoxel(xCoord, yCoord, zCoord, i, j, k, level);   //!all other voxels are fluid
					}

					//!level 0 without refinement
					else if ((level == 0) && (bc_.getLevelMax() == 0))
					{
						//!standard fluid voxels
						lattice_[level][i][j][k] = new FluidVoxel(xCoord, yCoord, zCoord, i, j, k, level);	//!all other voxels are fluid
					}

					//!level 1
					else if (level == 1)
					{
						//!interface voxels front and rear interface
						if ((i == refLayerX - 1 || i == refLayerX - 2 || i == imax - refLayer + 1 || i == imax - refLayer + 2) && ((j >= refLayer - 2) && (j <= jmax - refLayer + 2)) && ((k >= refLayer - 2) && (k <= kmax - refLayer + 2))) lattice_[level][i][j][k] = new InterfaceVoxel(xCoord, yCoord, zCoord, i, j, k, level); //!if voxel is at interface --> interface voxel

						//!interface voxels left and right interface
						else if ((j == refLayer - 1 || j == refLayer - 2 || j == jmax - refLayer + 1 || j == jmax - refLayer + 2) && ((i >= refLayerX - 2) && (i <= imax - refLayer + 2)) && ((k >= refLayer - 2) && (k <= kmax - refLayer + 2))) lattice_[level][i][j][k] = new InterfaceVoxel(xCoord, yCoord, zCoord, i, j, k, level); //!if voxel is at interface --> interface voxel

						//!interface voxels bottom and top interface
						else if ((k == refLayer - 1 || k == refLayer - 2 || k == kmax - refLayer + 1 || k == kmax - refLayer + 2) && ((i >= refLayerX - 2) && (i <= imax - refLayer + 2)) && ((j >= refLayer - 2) && (j <= jmax - refLayer + 2))) lattice_[level][i][j][k] = new InterfaceVoxel(xCoord, yCoord, zCoord, i, j, k, level); //!if voxel is at interface --> interface voxel

						//!standard fluid voxels
						else lattice_[level][i][j][k] = new FluidVoxel(xCoord, yCoord, zCoord, i, j, k, level); //!all other voxels are fluid
					}
				} //!end k
			} //!end j
		} //!end i
	} //!end level
} //!end createLattice


void Lattice::createMesh()
{
	//!create nodes in mesh structure
	//!mesh topology: bottom and top layer (+/- z normal) of mesh are periodic boundaries
	//!left and right boundaries (+/- y normal) are periodic boundaries
	//!front and back boundaries (+/- x normal) are periodic boundaries
	//!remaining nodes are fluid

	//!create nodes in array
	for (int level = 0; level <= bc_.getLevelMax(); level++)
	{
		int imax = (bc_.getNumberOfNodesX() - 1) * pow(2, level); //!index of last node in x direction in ->this level
		int jmax = (bc_.getNumberOfNodesY() - 1) * pow(2, level); //!index of last node in y direction in ->this level
		int kmax = (bc_.getNumberOfNodesZ() - 1) * pow(2, level); //!index of last node in z direction in ->this level

		for (int k = 0; k <= kmax; k++)
		{
			for (int j = 0; j <= jmax; j++)
			{
				for (int i = 0; i <= imax; i++)
				{
					double xCoord = i * (bc_.getSpacing() / pow(2., level));
					double yCoord = j * (bc_.getSpacing() / pow(2., level));
					double zCoord = k * (bc_.getSpacing() / pow(2., level));

					mesh_[level][i][j][k] = new Node(xCoord, yCoord, zCoord, i, j, k, level);
				} //!end k
			} //!end j
		} //!end i
	} //!end level
} //!end createMesh


void Lattice::setConnectivity()
{
	//!create connectivity, i.e., set neighbor pointer of node to the according neighbor node
	//!for wall nodes the neighbor pointer is set to NULL if neighbor is missing
	//!neighbor pointer for periodic nodes are set to periodic neighbor to achieve periodicity condition

	//!set connectivity for fluid nodes
	for (int level = 0; level <= bc_.getLevelMax(); level++)
	{
		int imax = bc_.getNumberOfVoxelsX() * pow(2., level) - 1; //!index of last voxel in x direction in ->this level
		int jmax = bc_.getNumberOfVoxelsY() * pow(2., level) - 1; //!index of last voxel in y direction in ->this level
		int kmax = bc_.getNumberOfVoxelsZ() * pow(2., level) - 1; //!index of last voxel in z direction in ->this level
		int refLayer = bc_.getRefinementLayer() * pow(2, level);
		int refLayerX = bc_.getRefinementLayerX() * pow(2, level);

		for (int i = 0; i <= imax; i++) //!all nodes in x direction except inlet and outlet boundaries
		{
			for (int j = 0; j <= jmax; j++) //!all nodes except left and right wall
			{
				for (int k = 0; k <= kmax; k++) //!all nodes in z direction except bottom and top wall
				{
					for (int dir = 0; dir < 18; dir++)
					{
						int iNeighbor, jNeighbor, kNeighbor; //!node indices of next neighbor in this direction
						switch (dir) //!calc indices of neighbor
						{
						case 0:  iNeighbor = i + 1; jNeighbor = j;		kNeighbor = k;     break;
						case 1:  iNeighbor = i;		jNeighbor = j + 1;	kNeighbor = k;	   break;
						case 2:	 iNeighbor = i - 1; jNeighbor = j;		kNeighbor = k;	   break;
						case 3:  iNeighbor = i;		jNeighbor = j - 1;	kNeighbor = k;     break;
						case 4:  iNeighbor = i;     jNeighbor = j;  	kNeighbor = k + 1; break;
						case 5:  iNeighbor = i;		jNeighbor = j;	    kNeighbor = k - 1; break;
						case 6:  iNeighbor = i + 1; jNeighbor = j;	    kNeighbor = k + 1; break;
						case 7:  iNeighbor = i + 1; jNeighbor = j;	    kNeighbor = k - 1; break;
						case 8:  iNeighbor = i;		jNeighbor = j + 1;	kNeighbor = k + 1; break;
						case 9:	 iNeighbor = i;		jNeighbor = j + 1;	kNeighbor = k - 1; break;
						case 10: iNeighbor = i - 1; jNeighbor = j;	    kNeighbor = k + 1; break;
						case 11: iNeighbor = i - 1; jNeighbor = j;	    kNeighbor = k - 1; break;
						case 12: iNeighbor = i;		jNeighbor = j - 1;	kNeighbor = k + 1; break;
						case 13: iNeighbor = i;		jNeighbor = j - 1;	kNeighbor = k - 1; break;
						case 14: iNeighbor = i - 1; jNeighbor = j + 1;	kNeighbor = k;     break;
						case 15: iNeighbor = i + 1; jNeighbor = j - 1;	kNeighbor = k;     break;
						case 16: iNeighbor = i + 1; jNeighbor = j + 1;	kNeighbor = k;     break;
						case 17: iNeighbor = i - 1; jNeighbor = j - 1;	kNeighbor = k;     break;
						} //!end switch

						//!consider periodic boundaries
						if (iNeighbor == imax + 1) iNeighbor = 0; else if (iNeighbor < 0) iNeighbor = imax; //!set neighbor to back or front at periodic boundaries
						if (jNeighbor == jmax + 1) jNeighbor = 0; else if (jNeighbor < 0) jNeighbor = jmax; //!set neighbor to left or right at periodic boundaries
						if (kNeighbor == kmax + 1) kNeighbor = 0; else if (kNeighbor < 0) kNeighbor = kmax; //!set neighbor to top or bottom at periodic boundaries

						Voxel* neighbor = lattice_[level][iNeighbor][jNeighbor][kNeighbor];
						//!check if neighbor was found
						if (neighbor == NULL) cerr << " Error: could not find neighbor at node i j k" << i << " " << j << " " << k << ". Direction is " << dir;
						lattice_[level][i][j][k]->setNeighbor(dir, neighbor);	//!set pointer to neighbor
					} //!end dir
					if ((level == 0) && (bc_.getLevelMax() > 0) && ((i <= refLayerX - 1) || (i >= imax - refLayer + 1) || (j <= refLayer - 1) || (j >= jmax - refLayer + 1) || (k <= refLayer - 1) || (k >= kmax - refLayer + 1)))
					{
						int n = 2;
						int ii = i * n;
						int jj = j * n;
						int kk = k * n;
						lattice_[level][i][j][k]->setPartner(lattice_[level + 1][ii][jj][kk], 0);
						lattice_[level][i][j][k]->setPartner(lattice_[level + 1][ii + 1][jj][kk], 1);
						lattice_[level][i][j][k]->setPartner(lattice_[level + 1][ii][jj + 1][kk], 2);
						lattice_[level][i][j][k]->setPartner(lattice_[level + 1][ii][jj][kk + 1], 3);
						lattice_[level][i][j][k]->setPartner(lattice_[level + 1][ii + 1][jj + 1][kk], 4);
						lattice_[level][i][j][k]->setPartner(lattice_[level + 1][ii + 1][jj][kk + 1], 5);
						lattice_[level][i][j][k]->setPartner(lattice_[level + 1][ii][jj + 1][kk + 1], 6);
						lattice_[level][i][j][k]->setPartner(lattice_[level + 1][ii + 1][jj + 1][kk + 1], 7);
					}
					else
					{
						Voxel* nullPtr = NULL;
						lattice_[level][i][j][k]->setPartner(nullPtr, 0);
						lattice_[level][i][j][k]->setPartner(nullPtr, 1);
						lattice_[level][i][j][k]->setPartner(nullPtr, 2);
						lattice_[level][i][j][k]->setPartner(nullPtr, 3);
						lattice_[level][i][j][k]->setPartner(nullPtr, 4);
						lattice_[level][i][j][k]->setPartner(nullPtr, 5);
						lattice_[level][i][j][k]->setPartner(nullPtr, 6);
						lattice_[level][i][j][k]->setPartner(nullPtr, 7);
					}
				} //!end k
			} //!end j
		} //!end i
	} //!end level
} //!end set connectivity



void Lattice::initialization(Control& bc) //!initialize all distribution values according to specified density and zero velocity
{
	int maxLevel = bc_.getLevelMax(); //!maximum level

	//!initialize populations based on zero or analytical maximum velocity equilibrium distribution
	for (int level = 0; level <= maxLevel; level++)	//!loop over all levels
	{
		double colFreq = (pow(bc_.getSpeedOfSound(), 2.) * (bc_.getTimeStep() / pow(2., level))) / (bc_.getKinematicViscosity() + 0.5 * pow(bc_.getSpeedOfSound(), 2.) * (bc_.getTimeStep() / pow(2., level))); //!calc dimensionless collision frequency
		double equil[19] = { 0. };
		double timeStep = bc_.getTimeStep() / pow(2., level);
		double dens = bc_.getDensity();
		double velo[3] = { 0 };

		int imax = bc_.getNumberOfVoxelsX() * pow(2., level) - 1; //!index of last voxel in x direction in ->this level
		int jmax = bc_.getNumberOfVoxelsY() * pow(2., level) - 1; //!index of last voxel in y direction in ->this level
		int kmax = bc_.getNumberOfVoxelsZ() * pow(2., level) - 1; //!index of last voxel in z direction in ->this level

		velo[0] = bc_.getInitVelo()[0];
		velo[1] = bc_.getInitVelo()[1];
		velo[2] = bc_.getInitVelo()[2];
		double Rc = 0.06;
		double xCenter = -6. * Rc;
		double yCenter = 0.;
		double cs = bc_.getSpeedOfSound();

		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all voxels in y direction
			{
				for (int k = 0; k <= kmax; k++) //!all voxels in z direction
				{
					lattice_[level][i][j][k]->setTimeStep(timeStep);
					double xNode = lattice_[level][i][j][k]->getX() - bc_.getChannelSizeX_() / 2.;
					double yNode = lattice_[level][i][j][k]->getY() - bc_.getChannelSizeY_() / 2.;

					//!initialize convected vortex
					//double epsilon = 0.15 * velo[0];
					double epsilon = 0.15 * cs;
					double expexpRho = exp(-(pow(xNode - xCenter, 2.) + pow(yNode - yCenter, 2.)) / (Rc * Rc));
					double expexpVelo = exp(-(pow(xNode - xCenter, 2.) + pow(yNode - yCenter, 2.)) / (2. * Rc * Rc));
					double rhoNode = dens * exp(-((epsilon * epsilon) / (2. * cs * cs)) * expexpRho);
					double uNode = velo[0] - epsilon * ((yNode - yCenter) / Rc) * expexpVelo;
					double vNode = epsilon * ((xNode - xCenter) / Rc) * expexpVelo;
					double veloNode[3] = { uNode, vNode, 0. };

                    if (bc.getOrderOfEquilibrium() == 2)
                    {
                        lattice_[level][i][j][k]->calcEquilibrium2_GHbasis(equil, bc, rhoNode, veloNode);
                    }
                    else if (bc.getOrderOfEquilibrium() == 3)
                    {
                        lattice_[level][i][j][k]->calcEquilibrium3_GHbasis(equil, bc, rhoNode, veloNode);
                    }
                    else if (bc.getOrderOfEquilibrium() == 4)
                    {
                        lattice_[level][i][j][k]->calcEquilibrium4_GHbasis(equil, bc, rhoNode, veloNode);
                    }
                    else
                    {
                        cerr << "Wrong order of equilibrium!" << flush;
                        exit(1);
                    }

					for (int dir = 0; dir <= 18; dir++) //!initialize all directions
					{
						lattice_[level][i][j][k]->setDistribution(true, dir, equil[dir]); //!set both distribution arrays
						lattice_[level][i][j][k]->setDistribution(false, dir, equil[dir]);
						lattice_[level][i][j][k]->setDistributionPreCol(dir, equil[dir]);
						lattice_[level][i][j][k]->setCubicMachCorrection(dir, 0.);
						lattice_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, 0.);
						lattice_[level][i][j][k]->setCollisionFrequency(colFreq); //!set collision frequency
					} //!end dir

					//!initialize strain rate and deviatoric stress tensor
					for (int x = 0; x < 3; x++)
					{
						for (int y = 0; y < 3; y++)
						{
							lattice_[level][i][j][k]->setPi1Tensor(x, y, 0.);
							lattice_[level][i][j][k]->setStrainRateTensor(x, y, 0.);
							lattice_[level][i][j][k]->setStrainRateTensorFD(x, y, 0.);
						}
					}
				} //!end k
			} //!end j
		} //!end i
	} //!end level
} //!end initialize

void Lattice::timeSteppingSRT(int& level, bool& alternating) //!perform a collide stream step for this level
{
	//!the boolean variable alternating is used to differentiate between the two distribution arrays (0/1)
	bool negAlternating = !alternating; //!accordingly by negAlternating is opposite array of distribution
	int imax = bc_.getNumberOfVoxelsX() - 1; //!index of last voxel in x direction in ->this level
	int jmax = bc_.getNumberOfVoxelsY() - 1; //!index of last voxel in y direction in ->this level
	int kmax = bc_.getNumberOfVoxelsZ() - 1; //!index of last voxel in z direction in ->this level

	//!safe pre-collision distribution on all voxels
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all voxels in y direction
			{
				for (int k = 0; k <= kmax; k++) //!all voxels in z direction
				{
					for (int dir = 0; dir < 19; dir++)
					{
						double distrTemp = lattice_[level][i][j][k]->getDistribution(alternating, dir); //!initialize temporary distribution
						lattice_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
					}
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

	//!collide and stream
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel;

		//!do collision and streaming for all fluid voxels
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all voxels in y direction execpt walls
			{
				for (int k = 0; k <= kmax; k++) //!all voxels in z direction execpt walls
				{
					lattice_[level][i][j][k]->collideSRT(bc_, alternating); //!perform collide at voxel
					lattice_[level][i][j][k]->transport(negAlternating); //!transport the distribution from alternating (just performed collide) to neighbors (write into negAlternating array)
				} //!end k
			} //!end j
		} //!end i
	} //!end omp
} //!end timeSteppingSRT


void Lattice::timeSteppingMRT_GSbasis(int& level, bool& alternating) //!perform a collide stream step for this level
{
	//!the boolean variable alternating is used to differentiate between the two distribution arrays (0/1)
	bool negAlternating = !alternating; //!accordingly by negAlternating is opposite array of distribution
	int imax = bc_.getNumberOfVoxelsX() - 1; //!index of last voxel in x direction in ->this level
	int jmax = bc_.getNumberOfVoxelsY() - 1; //!index of last voxel in y direction in ->this level
	int kmax = bc_.getNumberOfVoxelsZ() - 1; //!index of last voxel in z direction in ->this level

	//!safe pre-collision distribution on all voxels
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all voxels in y direction
			{
				for (int k = 0; k <= kmax; k++) //!all voxels in z direction
				{
					for (int dir = 0; dir < 19; dir++)
					{
						double distrTemp = lattice_[level][i][j][k]->getDistribution(alternating, dir); //!initialize temporary distribution
						lattice_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
					}
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

	//!collide and stream
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel;

		//!do collision and streaming for all fluid voxels
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all voxels in y direction execpt walls
			{
				for (int k = 0; k <= kmax; k++) //!all voxels in z direction execpt walls
				{
					lattice_[level][i][j][k]->collideMRT_GSbasis(bc_, alternating); //!perform collide at voxel
					lattice_[level][i][j][k]->transport(negAlternating); //!transport the distribution from alternating (just performed collide) to neighbors (write into negAlternating array)
				} //!end k
			} //!end j
		} //!end i
	} //!end omp
} //!end timeSteppingMRT_GSbasis


void Lattice::timeSteppingMRT_RMbasis(int& level, bool& alternating) //!perform a collide stream step for this level
{
	//!the boolean variable alternating is used to differentiate between the two distribution arrays (0/1)
	bool negAlternating = !alternating; //!accordingly by negAlternating is opposite array of distribution
	int imax = bc_.getNumberOfVoxelsX() - 1; //!index of last voxel in x direction in ->this level
	int jmax = bc_.getNumberOfVoxelsY() - 1; //!index of last voxel in y direction in ->this level
	int kmax = bc_.getNumberOfVoxelsZ() - 1; //!index of last voxel in z direction in ->this level

	//!safe pre-collision distribution on all voxels
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all voxels in y direction
			{
				for (int k = 0; k <= kmax; k++) //!all voxels in z direction
				{
					for (int dir = 0; dir < 19; dir++)
					{
						double distrTemp = lattice_[level][i][j][k]->getDistribution(alternating, dir); //!initialize temporary distribution
						lattice_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
					}
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

	//!collide and stream
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel;

		//!do collision and streaming for all fluid voxels
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all voxels in y direction execpt walls
			{
				for (int k = 0; k <= kmax; k++) //!all voxels in z direction execpt walls
				{
					lattice_[level][i][j][k]->collideMRT_RMbasis(bc_, alternating); //!perform collide at voxel
					lattice_[level][i][j][k]->transport(negAlternating); //!transport the distribution from alternating (just performed collide) to neighbors (write into negAlternating array)
				} //!end k
			} //!end j
		} //!end i
	} //!end omp
} //!end timeSteppingMRT_RMbasis


void Lattice::timeSteppingRR(int& level, bool& alternating) //!perform a collide stream step for this level
{
	//!the boolean variable alternating is used to differentiate between the two distribution arrays (0/1)
	bool negAlternating = !alternating; //!accordingly by negAlternating is opposite array of distribution
	int imax = bc_.getNumberOfVoxelsX() - 1; //!index of last voxel in x direction in ->this level
	int jmax = bc_.getNumberOfVoxelsY() - 1; //!index of last voxel in y direction in ->this level
	int kmax = bc_.getNumberOfVoxelsZ() - 1; //!index of last voxel in z direction in ->this level

	//!safe pre-collision distribution on all voxels
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all voxels in y direction
			{
				for (int k = 0; k <= kmax; k++) //!all voxels in z direction
				{
					for (int dir = 0; dir < 19; dir++)
					{
						double distrTemp = lattice_[level][i][j][k]->getDistribution(alternating, dir); //!initialize temporary distribution
						lattice_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
					}
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

	//!collide and stream
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel;

		//!do collision and streaming for all fluid voxels
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all voxels in y direction execpt walls
			{
				for (int k = 0; k <= kmax; k++) //!all voxels in z direction execpt walls
				{
					lattice_[level][i][j][k]->collideRR(bc_, alternating); //!perform collide at voxel
					lattice_[level][i][j][k]->transport(negAlternating); //!transport the distribution from alternating (just performed collide) to neighbors (write into negAlternating array)
				} //!end k
			} //!end j
		} //!end i
	} //!end omp
} //!end timeSteppingRR


void Lattice::timeSteppingHRR(int& level, bool& alternating) //!perform a collide stream step for this level
{
	//!the boolean variable alternating is used to differentiate between the two distribution arrays (0/1)
	bool negAlternating = !alternating; //!accordingly by negAlternating is opposite array of distribution
	int imax = bc_.getNumberOfVoxelsX() - 1; //!index of last voxel in x direction in ->this level
	int jmax = bc_.getNumberOfVoxelsY() - 1; //!index of last voxel in y direction in ->this level
	int kmax = bc_.getNumberOfVoxelsZ() - 1; //!index of last voxel in z direction in ->this level

	//!safe pre-collision distribution on all voxels
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all voxels in y direction
			{
				for (int k = 0; k <= kmax; k++) //!all voxels in z direction
				{
					for (int dir = 0; dir < 19; dir++)
					{
						double distrTemp = lattice_[level][i][j][k]->getDistribution(alternating, dir); //!initialize temporary distribution
						lattice_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
					}
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

	//!collide and stream
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel;

		//!do collision and streaming for all fluid voxels
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all voxels in y direction execpt walls
			{
				for (int k = 0; k <= kmax; k++) //!all voxels in z direction execpt walls
				{
					lattice_[level][i][j][k]->collideHRR(bc_, alternating); //!perform collide at voxel
					lattice_[level][i][j][k]->transport(negAlternating); //!transport the distribution from alternating (just performed collide) to neighbors (write into negAlternating array)
				} //!end k
			} //!end j
		} //!end i
	} //!end omp
} //!end timeSteppingHRR


void Lattice::solve(bool& alternatingL0)
{
	int level = 0;
	string collisionModel = bc_.getCollisionModel();

	if (collisionModel == "SRT")
	{
		this->timeSteppingSRT(level, alternatingL0);
	}
	else if (collisionModel == "MRT-GS")
	{
		this->timeSteppingMRT_GSbasis(level, alternatingL0);
	}
	else if (collisionModel == "MRT-RM")
	{
		this->timeSteppingMRT_RMbasis(level, alternatingL0);
	}
	else if (collisionModel == "RR")
	{
		this->timeSteppingRR(level, alternatingL0);
	}
	else if (collisionModel == "HRR")
	{
		this->timeSteppingHRR(level, alternatingL0);
	}
} //!end solve

void Lattice::solve(bool& alternatingL0, bool& alternatingL1)
{
	string collisionModel = bc_.getCollisionModel();

	if (collisionModel == "SRT")
	{
		this->nestedTimeSteppingSRT(alternatingL1, alternatingL0, bc_);
	}
	else if (collisionModel == "MRT-GS")
	{
		this->nestedTimeSteppingMRT_GSbasis(alternatingL1, alternatingL0, bc_);
	}
	else if (collisionModel == "MRT-RM")
	{
		this->nestedTimeSteppingMRT_RMbasis(alternatingL1, alternatingL0, bc_);
	}
	else if (collisionModel == "RR")
	{
		this->nestedTimeSteppingRR(alternatingL1, alternatingL0, bc_);
	}
	else if (collisionModel == "HRR")
	{
		this->nestedTimeSteppingHRR(alternatingL1, alternatingL0, bc_);
	}
} //!end solve


void Lattice::writeResultsVTK(string filename, bool& alternatingL0, bool& alternatingL1)
{
	int refLayer = bc_.getRefinementLayer();	//!number of refined cells on Level 0
	int refLayerX = bc_.getRefinementLayerX();	//!number of refined cells on Level 0

	filename.append(".vtk");
	ofstream paraview(filename.c_str()); //!open file stream to write

	if (!paraview) //!if stream is corrupted, quit
	{
		cerr << "An error occured while writing paraview file\n";
		return;
	}
	else
	{
		//!write header of paraview file
		paraview << "# vtk DataFile Version 2.0" << "\n";
		paraview << "LBM results for grid refinement benchmark in channel flow\n";
		paraview << "ASCII" << "\n";
		paraview << "DATASET UNSTRUCTURED_GRID" << "\n";

		paraview << "POINTS ";
		int pntPos = static_cast<int>(paraview.tellp()); //!remember stream position to write number of nodes to the file (this will get important later)
		paraview << "              " << " float" << "\n";
		paraview << setiosflags(ios::left | ios::showpoint | ios::scientific) << setprecision(8);

		int nodeCntrBot = 0; //!number of fine nodes in bottom region
		int nodeCntrTop = 0; //!number of fine nodes in top region
		int nodeCntrLft = 0; //!number of fine nodes in left region
		int nodeCntrRgt = 0; //!number of fine nodes in right region
		int nodeCntrInl = 0; //!number of fine nodes in inlet region
		int nodeCntrOut = 0; //!number of fine nodes in outlet region
		int nodeCntrMid = 0; //!number of coarse nodes in mid region
		int pntCntr     = 0;

		if (bc_.getLevelMax() > 0)
		{
			//------------------------------------------//
			//        coordinates in *.vtk-file:		//
			//------------------------------------------//

			//---------------------//
			// bottom part level 1 //
			//---------------------//

			int level = 1;
			int imax = (bc_.getNumberOfNodesX() - 1) * pow(2, level); //!index of last node in x direction in ->this level
			int jmax = (bc_.getNumberOfNodesY() - 1) * pow(2, level); //!index of last node in y direction in ->this level
			int kmax = (bc_.getNumberOfNodesZ() - 1) * pow(2, level); //!index of last node in z direction in ->this level
			refLayer = bc_.getRefinementLayer() * pow(2, level);
            refLayerX = bc_.getRefinementLayerX() * pow(2, level);

			for (int k = 0; k <= refLayer - 2; k++) //!all nodes in z direction in refined bottom region
			{
				for (int j = 0; j <= jmax; j++) //!all nodes in y direction
				{
					for (int i = 0; i <= imax; i++) //!all nodes in x direction
					{
						if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
						{
							paraview << mesh_[level][i][j][k]->getX() << " " << mesh_[level][i][j][k]->getY() << " " << mesh_[level][i][j][k]->getZ() << "\n"; //!write coordinates
							nodeCntrBot++; //!count nodes
						}
					} //!end i
				} //!end j
			} //!end k

			//------------------//
			// top part level 1 //
			//------------------//

			for (int k = kmax - refLayer + 2; k <= kmax; k++) //!all nodes in z direction in refined top region
			{
				for (int j = 0; j <= jmax; j++) //!all nodes in y direction
				{
					for (int i = 0; i <= imax; i++) //!all nodes in x direction
					{
						if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
						{
							paraview << mesh_[level][i][j][k]->getX() << " " << mesh_[level][i][j][k]->getY() << " " << mesh_[level][i][j][k]->getZ() << "\n"; //!write coordinates
							nodeCntrTop++; //!count nodes
						}
					} //!end i
				} //!end j
			} //!end k

			//-------------------//
			// left part level 1 //
			//-------------------//

			for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all nodes in z direction
			{
				for (int j = 0; j <= refLayer - 2; j++) //!all nodes in y direction in refined left region
				{
					for (int i = 0; i <= imax; i++) //!all nodes in x direction
					{
						if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
						{
							paraview << mesh_[level][i][j][k]->getX() << " " << mesh_[level][i][j][k]->getY() << " " << mesh_[level][i][j][k]->getZ() << "\n"; //!write coordinates
							nodeCntrLft++; //!count nodes
						}
					} //!end i
				} //!end j
			} //!end k

			//--------------------//
			// right part level 1 //
			//--------------------//

			for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all nodes in z direction except those considered before
			{
				for (int j = jmax - refLayer + 2; j <= jmax; j++) //!all nodes in y direction in refined right region
				{
					for (int i = 0; i <= imax; i++) //!all nodes in x direction
					{
						if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
						{
							paraview << mesh_[level][i][j][k]->getX() << " " << mesh_[level][i][j][k]->getY() << " " << mesh_[level][i][j][k]->getZ() << "\n"; //!write coordinates
							nodeCntrRgt++; //!count nodes
						}
					} //!end i
				} //!end j
			} //!end k

			//--------------------//
			// inlet part level 1 //
			//--------------------//

			for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all nodes in z direction except those considered before
			{
				for (int j = refLayer - 2; j <= jmax - refLayer + 2; j++) //!all nodes in y direction except those considered before
				{
					for (int i = 0; i <= refLayerX - 2; i++) //!all nodes in x direction in refined inlet region
					{
						if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
						{
							paraview << mesh_[level][i][j][k]->getX() << " " << mesh_[level][i][j][k]->getY() << " " << mesh_[level][i][j][k]->getZ() << "\n"; //!write coordinates
							nodeCntrInl++; //!count nodes
						}
					} //!end i
				} //!end j
			} //!end k

			//---------------------//
			// outlet part level 1 //
			//---------------------//

			for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all nodes in z direction except those considered before
			{
				for (int j = refLayer - 2; j <= jmax - refLayer + 2; j++) //!all nodes in y direction except those considered before
				{
					for (int i = imax - refLayer + 2; i <= imax; i++) //!all nodes in x direction in refined outlet region
					{
						if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
						{
							paraview << mesh_[level][i][j][k]->getX() << " " << mesh_[level][i][j][k]->getY() << " " << mesh_[level][i][j][k]->getZ() << "\n"; //!write coordinates
							nodeCntrOut++; //!count nodes
						}
					} //!end i
				} //!end j
			} //!end k

			//------------------//
			// mid part level 0 //
			//------------------//

			level = 0;
			imax = bc_.getNumberOfNodesX() - 1; //!index of last node in x direction in ->this level
			jmax = bc_.getNumberOfNodesY() - 1; //!index of last node in y direction in ->this level
			kmax = bc_.getNumberOfNodesZ() - 1; //!index of last node in z direction in ->this level
			refLayer = bc_.getRefinementLayer();
			refLayerX = bc_.getRefinementLayerX() * pow(2, level);

			for (int k = refLayer - 1; k <= kmax - refLayer + 1; k++) //!all nodes in z direction in unrefined mid region
			{
				for (int j = refLayer - 1; j <= jmax - refLayer + 1; j++) //!all nodes in y direction
				{
					for (int i = refLayerX - 1; i <= imax - refLayer + 1; i++) //!all nodes in x direction
					{
						if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
						{
							paraview << mesh_[level][i][j][k]->getX() << " " << mesh_[level][i][j][k]->getY() << " " << mesh_[level][i][j][k]->getZ() << "\n"; //!write coordinates
							nodeCntrMid++; //!count nodes
						}
					} //!end i
				} //!end j
			} //!end k

			//!assemble and assign index
			paraview << "\n";
			int nodeCntrTot = nodeCntrBot + nodeCntrTop + nodeCntrLft + nodeCntrRgt + nodeCntrInl + nodeCntrOut; //!total number of fine nodes
			pntCntr = nodeCntrTot + nodeCntrMid; //!store numberOfPoints

			paraview << "CELLS "; //!write cells
			int cellPos = static_cast<int>(paraview.tellp()); //!remember cell position
			int cellCntr = 0;
			paraview << "         " << " " << "         " << "\n";

			//------------------------------------------//
			//		  connectivity in *.vtk-file		//
			//------------------------------------------//

			//---------------------//
			// bottom part level 1 //
			//---------------------//

			level = 1;
			imax = (bc_.getNumberOfNodesX() - 1) * pow(2, level); //!index of last node in x direction in ->this level
			jmax = (bc_.getNumberOfNodesY() - 1) * pow(2, level); //!index of last node in y direction in ->this level
			kmax = (bc_.getNumberOfNodesZ() - 1) * pow(2, level); //!index of last node in z direction in ->this level
			refLayer = bc_.getRefinementLayer() * pow(2, level);
			refLayerX = bc_.getRefinementLayerX() * pow(2, level);

			int icmax = bc_.getNumberOfVoxelsX() * pow(2, level) - 1; //!index of last voxel in x direction in ->this level
			int jcmax = bc_.getNumberOfVoxelsY() * pow(2, level) - 1; //!index of last voxel in y direction in ->this level
			int kcmax = bc_.getNumberOfVoxelsZ() * pow(2, level) - 1; //!index of last voxel in z direction in ->this level

			int nodeCntrX = imax + 1; //!number of nodes in x direction in this level
			int nodeCntrY = jmax + 1; //!number of nodes in y direction in this level
			int nodeCntrXY = nodeCntrX * nodeCntrY; //!number of nodes per xy plane in this level
			int nodeCntrRef = refLayer - 1;
			int nodeCntrRefX = refLayerX - 1;
			int nodeCntrZ = 0;
			int nodeCntrJ = 0;
			int nodeCntrI = 0;

			for (int k = 0; k <= refLayer - 3; k++) //!all cells in z direction in refined bottom region
			{
				for (int j = 0; j <= jcmax; j++) //!all cells in y direction
				{
					for (int i = 0; i <= icmax; i++) //!all cells in x direction
					{
						if (lattice_[level][i][j][k] != NULL) //!write only existing voxels
						{
							//!calc index of node 0 (i.e. line in paraview file) from i,j and k
							int firstNodeIndex = i + j * nodeCntrX + k * nodeCntrXY;

							//!write connectivity
							paraview << "8 "; //!write number of points in cell (i.e. number of variables to read in this line)
							paraview << firstNodeIndex << " " << firstNodeIndex + 1 << " " << firstNodeIndex + nodeCntrX << " " << firstNodeIndex + nodeCntrX + 1 << " "; //!nodes 0 to 3
							paraview << firstNodeIndex + nodeCntrXY << " " << firstNodeIndex + nodeCntrXY + 1 << " "; //!nodes 4 and 5
							paraview << firstNodeIndex + nodeCntrXY + nodeCntrX << " " << firstNodeIndex + nodeCntrXY + nodeCntrX + 1 << "\n"; //!nodes 6 and 7
							cellCntr++; //!count cells
						}
					} //!end i
				} //!end j
			} //!end k

			//------------------//
			// top part level 1 //
			//------------------//

			nodeCntrZ = 0;

			for (int k = kcmax - refLayer + 3; k <= kcmax; k++) //!all cells in z direction in refined top region
			{
				for (int j = 0; j <= jcmax; j++) //!all cells in y direction
				{
					for (int i = 0; i <= icmax; i++) //!all cells in x direction
					{
						if (lattice_[level][i][j][k] != NULL) //!write only existing nodes
						{
							//!calc index of node 0 (i.e. line in paraview file) from i,j and k
							int firstNodeIndex = i + j * nodeCntrX + nodeCntrZ * nodeCntrXY + nodeCntrBot;

							//!write connectivity
							paraview << "8 "; //!write number of points in cell (i.e. number of variables to read in this line)
							paraview << firstNodeIndex << " " << firstNodeIndex + 1 << " " << firstNodeIndex + nodeCntrX << " " << firstNodeIndex + nodeCntrX + 1 << " "; //!nodes 0 to 3
							paraview << firstNodeIndex + nodeCntrXY << " " << firstNodeIndex + nodeCntrXY + 1 << " "; //!nodes 4 and 5
							paraview << firstNodeIndex + nodeCntrXY + nodeCntrX << " " << firstNodeIndex + nodeCntrXY + nodeCntrX + 1 << "\n"; //!nodes 6 and 7
							cellCntr++; //!count cells

						}
					} //!end i
				} //!end j
				nodeCntrZ++;
			} //!end k

			//-------------------//
			// left part level 1 //
			//-------------------//

			nodeCntrZ = 0;

			for (int k = refLayer - 2; k <= kcmax - refLayer + 2; k++) //!all cells in z direction in refined left region except cells already considered before
			{
				for (int j = 0; j <= refLayer - 3; j++) //!all cells in y direction in refined left region
				{
					for (int i = 0; i <= icmax; i++) //!all cells in x direction
					{
						if (lattice_[level][i][j][k] != NULL) //!write only existing nodes
						{
							//!calc index of node 0 (i.e. line in paraview file) from i,j and k
							int firstNodeIndex = i + j * nodeCntrX + nodeCntrZ * nodeCntrRef * nodeCntrX + nodeCntrBot + nodeCntrTop;

							//!write connectivity
							paraview << "8 "; //!write number of points in cell (i.e. number of variables to read in this line)
							paraview << firstNodeIndex << " " << firstNodeIndex + 1 << " " << firstNodeIndex + nodeCntrX << " " << firstNodeIndex + nodeCntrX + 1 << " "; //!nodes 0 to 3
							paraview << firstNodeIndex + nodeCntrRef * nodeCntrX << " " << firstNodeIndex + nodeCntrRef * nodeCntrX + 1 << " "; //!nodes 4 and 5
							paraview << firstNodeIndex + nodeCntrRef * nodeCntrX + nodeCntrX << " " << firstNodeIndex + nodeCntrRef * nodeCntrX + nodeCntrX + 1 << "\n"; //!nodes 6 and 7
							cellCntr++; //!count cells
						}
					} //!end i
				} //!end j
				nodeCntrZ++;
			} //!end k

			//--------------------//
			// right part level 1 //
			//--------------------//

			nodeCntrZ = 0;

			for (int k = refLayer - 2; k <= kcmax - refLayer + 2; k++) //!all cells in z direction in refined right region except nodes already considered before
			{
				for (int j = jcmax - refLayer + 3; j <= jcmax; j++) //!all cells in y direction in refined right region
				{
					for (int i = 0; i <= icmax; i++) //!all cells in x direction
					{
						if (lattice_[level][i][j][k] != NULL) //!write only existing nodes
						{
							//!calc index of node 0 (i.e. line in paraview file) from i,j and k
							int firstNodeIndex = i + nodeCntrJ * nodeCntrX + nodeCntrZ * nodeCntrRef * nodeCntrX + nodeCntrBot + nodeCntrTop + nodeCntrLft;

							//!write connectivity
							paraview << "8 "; //!write number of points in cell (i.e. number of variables to read in this line)
							paraview << firstNodeIndex << " " << firstNodeIndex + 1 << " " << firstNodeIndex + nodeCntrX << " " << firstNodeIndex + nodeCntrX + 1 << " "; //!nodes 0 to 3
							paraview << firstNodeIndex + nodeCntrRef * nodeCntrX << " " << firstNodeIndex + nodeCntrRef * nodeCntrX + 1 << " "; //!nodes 4 and 5
							paraview << firstNodeIndex + nodeCntrRef * nodeCntrX + nodeCntrX << " " << firstNodeIndex + nodeCntrRef * nodeCntrX + nodeCntrX + 1 << "\n"; //!nodes 6 and 7
							cellCntr++; //!count cells
						}
					} //!end i
					nodeCntrJ++;
				} //!end j
				nodeCntrZ++;
				nodeCntrJ = 0;
			} //!end k

			//--------------------//
			// inlet part level 1 //
			//--------------------//

			nodeCntrZ = 0;
			nodeCntrJ = 0;

			for (int k = refLayer - 2; k <= kcmax - refLayer + 2; k++) //!all cells in z direction in refined inlet region except nodes already considered before
			{
				for (int j = refLayer - 2; j <= jcmax - refLayer + 2; j++) //!all cells in y direction in refined inlet region except nodes already considered before
				{
					for (int i = 0; i <= refLayerX - 3; i++) //!all cells in x direction in refined inlet region
					{
						if (lattice_[level][i][j][k] != NULL) //!write only existing nodes
						{
							//!calc index of node 0 (i.e. line in paraview file) from i,j and k
							int firstNodeIndex = i + nodeCntrJ * nodeCntrRefX + nodeCntrZ * nodeCntrRefX * (nodeCntrY - 2 * (nodeCntrRef - 1)) + nodeCntrBot + nodeCntrTop + nodeCntrLft + nodeCntrRgt;

							//!write connectivity
							paraview << "8 "; //!write number of points in cell (i.e. number of variables to read in this line)
							paraview << firstNodeIndex << " " << firstNodeIndex + 1 << " " << firstNodeIndex + nodeCntrRefX << " " << firstNodeIndex + nodeCntrRefX + 1 << " "; //!nodes 0 to 3
							paraview << firstNodeIndex + nodeCntrRefX * (nodeCntrY - 2 * (nodeCntrRef - 1)) << " " << firstNodeIndex + nodeCntrRefX * (nodeCntrY - 2 * (nodeCntrRef - 1)) + 1 << " "; //!nodes 4 and 5
							paraview << firstNodeIndex + nodeCntrRefX * (nodeCntrY - 2 * (nodeCntrRef - 1)) + nodeCntrRefX << " " << firstNodeIndex + nodeCntrRefX * (nodeCntrY - 2 * (nodeCntrRef - 1)) + nodeCntrRefX + 1 << "\n"; //!nodes 6 and 7
							cellCntr++; //!count cells
						}
					} //!end i
					nodeCntrJ++;
				} //!end j
				nodeCntrZ++;
				nodeCntrJ = 0;
			} //!end k

			//---------------------//
			// outlet part level 1 //
			//---------------------//

			nodeCntrZ = 0;
			nodeCntrJ = 0;

			for (int k = refLayer - 2; k <= kcmax - refLayer + 2; k++) //!all cells in z direction in refined outlet region except nodes already considered before
			{
				for (int j = refLayer - 2; j <= jcmax - refLayer + 2; j++) //!all cells in y direction in refined outlet region except nodes already considered before
				{
					for (int i = icmax - refLayer + 3; i <= icmax; i++) //!all cells in x direction in refined outlet region
					{
						if (lattice_[level][i][j][k] != NULL) //!write only existing nodes
						{
							//!calc index of node 0 (i.e. line in paraview file) from i,j and k
							int firstNodeIndex = nodeCntrI + nodeCntrJ * nodeCntrRef + nodeCntrZ * nodeCntrRef * (nodeCntrY - 2 * (nodeCntrRef - 1)) + nodeCntrBot + nodeCntrTop + nodeCntrLft + nodeCntrRgt + nodeCntrInl;

							//!write connectivity
							paraview << "8 "; //!write number of points in cell (i.e. number of variables to read in this line)
							paraview << firstNodeIndex << " " << firstNodeIndex + 1 << " " << firstNodeIndex + nodeCntrRef << " " << firstNodeIndex + nodeCntrRef + 1 << " "; //!nodes 0 to 3
							paraview << firstNodeIndex + nodeCntrRef * (nodeCntrY - 2 * (nodeCntrRef - 1)) << " " << firstNodeIndex + nodeCntrRef * (nodeCntrY - 2 * (nodeCntrRef - 1)) + 1 << " "; //!nodes 4 and 5
							paraview << firstNodeIndex + nodeCntrRef * (nodeCntrY - 2 * (nodeCntrRef - 1)) + nodeCntrRef << " " << firstNodeIndex + nodeCntrRef * (nodeCntrY - 2 * (nodeCntrRef - 1)) + nodeCntrRef + 1 << "\n"; //!nodes 6 and 7
							cellCntr++; //!count cells
						}
						nodeCntrI++;
					} //!end i
					nodeCntrJ++;
					nodeCntrI = 0;
				} //!end j
				nodeCntrZ++;
				nodeCntrJ = 0;
			} //!end k

			//------------------//
			// mid part level 0 //
			//------------------//

			level = 0;
			imax = bc_.getNumberOfNodesX() - 1; //!index of last node in x direction in ->this level
			jmax = bc_.getNumberOfNodesY() - 1; //!index of last node in y direction in ->this level
			kmax = bc_.getNumberOfNodesZ() - 1; //!index of last node in z direction in ->this level
			refLayer = bc_.getRefinementLayer();
			refLayerX = bc_.getRefinementLayerX() * pow(2, level);

			icmax = bc_.getNumberOfVoxelsX() - 1; //!index of last voxel in x direction in ->this level
			jcmax = bc_.getNumberOfVoxelsY() - 1; //!index of last voxel in y direction in ->this level
			kcmax = bc_.getNumberOfVoxelsZ() - 1; //!index of last voxel in z direction in ->this level

			nodeCntrRef = refLayer - 1;
			nodeCntrRefX = refLayerX - 1;

			nodeCntrX = imax + 1; //!number of nodes in x direction in ->this level
			nodeCntrY = jmax + 1; //!number of nodes in y direction in ->this level
			nodeCntrXY = (nodeCntrX - (nodeCntrRefX + nodeCntrRef)) * (nodeCntrY - 2 * nodeCntrRef); //!nodes per xy plane in this level

			paraview << "\n";

			nodeCntrZ = 0;
			nodeCntrZ = 0;
			nodeCntrJ = 0;
			nodeCntrI = 0;

			for (int k = refLayer - 1; k <= kcmax - refLayer + 1; k++) //!all cells in z direction in unrefined mid region
			{
				for (int j = refLayer - 1; j <= jcmax - refLayer + 1; j++) //!all cells in y direction
				{
					for (int i = refLayerX - 1; i <= icmax - refLayer + 1; i++) //!all cells in x direction
					{
						if (lattice_[level][i][j][k] != NULL) //!write only existing nodes
						{
							//!calc index of node 0 (i.e. line in paraview file) from i,j and k
							int firstNodeIndex = nodeCntrI + nodeCntrJ * (nodeCntrX - (nodeCntrRefX + nodeCntrRef)) + nodeCntrZ * nodeCntrXY + nodeCntrTot;

							//!write connectivity
							paraview << "8 "; //!write number of points in cell (i.e. number of variables to read in this line)
							paraview << firstNodeIndex << " " << firstNodeIndex + 1 << " " << firstNodeIndex + nodeCntrX - (nodeCntrRefX + nodeCntrRef) << " " << firstNodeIndex + nodeCntrX - (nodeCntrRefX + nodeCntrRef) + 1 << " "; //!nodes 0 to 3
							paraview << firstNodeIndex + nodeCntrXY << " " << firstNodeIndex + nodeCntrXY + 1 << " "; //!nodes 4 and 5
							paraview << firstNodeIndex + nodeCntrXY + nodeCntrX - (nodeCntrRefX + nodeCntrRef) << " " << firstNodeIndex + nodeCntrXY + nodeCntrX - (nodeCntrRefX + nodeCntrRef) + 1 << "\n"; //!nodes 6 and 7
							cellCntr++; //!count cell
							nodeCntrI++;
						}
					} //!end i
					nodeCntrJ++;
					nodeCntrI = 0;
				} //!end j
				nodeCntrZ++;
				nodeCntrJ = 0;
			} //!end k

			paraview << "\n";
			paraview << "CELL_TYPES " << cellCntr << "\n"; //!write cell types to file

			for (int i = 0; i < cellCntr; i++)		//!write cell types to paraview file
			{
				paraview << "11" << "\n";
			}
			paraview << "\n";

			//-------------------------//
			// solutions in *.vtk-file //
			//-------------------------//

			//->> COLLISION FREQUENCY <<-//

			paraview << "CELL_DATA " << cellCntr << "\n";
			paraview << "SCALARS Omega float 1" << "\n";
			paraview << "LOOKUP_TABLE default" << "\n";

			//---------------------//
			// bottom part level 1 //
			//---------------------//

			level = 1;
			imax = bc_.getNumberOfVoxelsX() * pow(2, level) - 1; //!index of last voxel in x direction in ->this level
			jmax = bc_.getNumberOfVoxelsY() * pow(2, level) - 1; //!index of last voxel in y direction in ->this level
			kmax = bc_.getNumberOfVoxelsZ() * pow(2, level) - 1; //!index of last voxel in z direction in ->this level
			refLayer = bc_.getRefinementLayer() * pow(2, level);
			refLayerX = bc_.getRefinementLayerX() * pow(2, level);

			for (int k = 0; k <= refLayer - 3; k++) //!all voxels in z direction in refined bottom region
			{
				for (int j = 0; j <= jmax; j++) //!all voxels in y direction
				{
					for (int i = 0; i <= imax; i++) //!all voxels in x direction
					{
						if (lattice_[level][i][j][k] != NULL) //!write only existing voxels
						{
							paraview << lattice_[level][i][j][k]->getCollisionFrequency() << "\n";
						}
					} //!end i
				} //!end j
			} //!end k

			//------------------//
			// top part level 1 //
			//------------------//

			for (int k = kmax - refLayer + 3; k <= kmax; k++) //!all voxels in z direction in refined top region
			{
				for (int j = 0; j <= jmax; j++) //!all voxels in y direction
				{
					for (int i = 0; i <= imax; i++) //!all voxels in x direction
					{
						if (lattice_[level][i][j][k] != NULL) //!write only existing voxels
						{
							paraview << lattice_[level][i][j][k]->getCollisionFrequency() << "\n";
						}
					} //!end i
				} //!end j
			} //!end k

			//-------------------//
			// left part level 1 //
			//-------------------//

			for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction in refined top region
			{
				for (int j = 0; j <= refLayer - 3; j++) //!all voxels in y direction
				{
					for (int i = 0; i <= imax; i++) //!all voxels in x direction
					{
						if (lattice_[level][i][j][k] != NULL) //!write only existing voxels
						{
							paraview << lattice_[level][i][j][k]->getCollisionFrequency() << "\n";
						}
					} //!end i
				} //!end j
			} //!end k

			//--------------------//
			// right part level 1 //
			//--------------------//

			for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction in refined right region
			{
				for (int j = jmax - refLayer + 3; j <= jmax; j++) //!all voxels in y direction
				{
					for (int i = 0; i <= imax; i++) //!all voxels in x direction
					{
						if (lattice_[level][i][j][k] != NULL) //!write only existing voxels
						{
							paraview << lattice_[level][i][j][k]->getCollisionFrequency() << "\n";
						}
					} //!end i
				} //!end j
			} //!end k

			//--------------------//
			// inlet part level 1 //
			//--------------------//

			for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction in refined inlet region
			{
				for (int j = refLayer - 2; j <= jmax - refLayer + 2; j++) //!all voxels in y direction
				{
					for (int i = 0; i <= refLayerX - 3; i++) //!all voxels in x direction
					{
						if (lattice_[level][i][j][k] != NULL) //!write only existing voxels
						{
							paraview << lattice_[level][i][j][k]->getCollisionFrequency() << "\n";
						}
					} //!end i
				} //!end j
			} //!end k

			//---------------------//
			// outlet part level 1 //
			//---------------------//

			for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction in refined outlet region
			{
				for (int j = refLayer - 2; j <= jmax - refLayer + 2; j++) //!all voxels in y direction
				{
					for (int i = imax - refLayer + 3; i <= imax; i++) //!all voxels in x direction
					{
						if (lattice_[level][i][j][k] != NULL) //!write only existing voxels
						{
							paraview << lattice_[level][i][j][k]->getCollisionFrequency() << "\n";
						}
					} //!end i
				} //!end j
			} //!end k

			//------------------//
			// mid part level 0 //
			//------------------//

			level = 0;
			imax = bc_.getNumberOfVoxelsX() - 1; //!index of last node in x direction in ->this level
			jmax = bc_.getNumberOfVoxelsY() - 1; //!index of last node in y direction in ->this level
			kmax = bc_.getNumberOfVoxelsZ() - 1; //!index of last node in z direction in ->this level
			refLayer = bc_.getRefinementLayer();
			refLayerX = bc_.getRefinementLayerX() * pow(2, level);

			for (int k = refLayer - 1; k <= kmax - refLayer + 1; k++) //!all voxels in z direction in unrefined mid region
			{
				for (int j = refLayer - 1; j <= jmax - refLayer + 1; j++) //!all voxels in y direction
				{
					for (int i = refLayerX - 1; i <= imax - refLayer + 1; i++) //!all voxels in x direction
					{
						if (lattice_[level][i][j][k] != NULL) //!write only existing voxels
						{
							paraview << lattice_[level][i][j][k]->getCollisionFrequency() << "\n";
						}
					} //!end i
				} //!end j
			} //!end k

			//->> PRESSURE <<-//

			paraview << "SCALARS Pressure float 1" << "\n";
			paraview << "LOOKUP_TABLE default" << "\n";
			double speedOfSoundSq = pow(bc_.getSpeedOfSound(), 2.); //!square of speed of sound

			//---------------------//
			// bottom part level 1 //
			//---------------------//

			level = 1;
			imax = bc_.getNumberOfVoxelsX() * pow(2, level) - 1; //!index of last voxel in x direction in ->this level
			jmax = bc_.getNumberOfVoxelsY() * pow(2, level) - 1; //!index of last voxel in y direction in ->this level
			kmax = bc_.getNumberOfVoxelsZ() * pow(2, level) - 1; //!index of last voxel in z direction in ->this level
			refLayer = bc_.getRefinementLayer() * pow(2, level);
			refLayerX = bc_.getRefinementLayerX() * pow(2, level);

			for (int k = 0; k <= refLayer - 3; k++) //!all voxels in z direction in refined bottom region
			{
				for (int j = 0; j <= jmax; j++) //!all voxels in y direction
				{
					for (int i = 0; i <= imax; i++) //!all voxels in x direction
					{
						if (lattice_[level][i][j][k] != NULL) //!write only existing voxels
						{
							paraview << lattice_[level][i][j][k]->calcDensity(bc_, alternatingL1) * speedOfSoundSq << "\n"; //!for fluid nodes calculate pressure
						}
					} //!end i
				} //!end j
			} //!end k

			//------------------//
			// top part level 1 //
			//------------------//

			for (int k = kmax - refLayer + 3; k <= kmax; k++) //!all voxels in z direction in refined top region
			{
				for (int j = 0; j <= jmax; j++) //!all voxels in y direction
				{
					for (int i = 0; i <= imax; i++) //!all voxels in x direction
					{
						if (lattice_[level][i][j][k] != NULL) //!write only existing voxels
						{
							paraview << lattice_[level][i][j][k]->calcDensity(bc_, alternatingL1) * speedOfSoundSq << "\n"; //!for fluid nodes calculate pressure
						}
					} //!end i
				} //!end j
			} //!end k

			//-------------------//
			// left part level 1 //
			//-------------------//

			for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction in refined top region
			{
				for (int j = 0; j <= refLayer - 3; j++) //!all voxels in y direction
				{
					for (int i = 0; i <= imax; i++) //!all voxels in x direction
					{
						if (lattice_[level][i][j][k] != NULL) //!write only existing voxels
						{
							paraview << lattice_[level][i][j][k]->calcDensity(bc_, alternatingL1) * speedOfSoundSq << "\n"; //!for fluid nodes calculate pressure
						}
					} //!end i
				} //!end j
			} //!end k

			//--------------------//
			// right part level 1 //
			//--------------------//

			for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction in refined right region
			{
				for (int j = jmax - refLayer + 3; j <= jmax; j++) //!all voxels in y direction
				{
					for (int i = 0; i <= imax; i++) //!all voxels in x direction
					{
						if (lattice_[level][i][j][k] != NULL) //!write only existing voxels
						{
							paraview << lattice_[level][i][j][k]->calcDensity(bc_, alternatingL1) * speedOfSoundSq << "\n"; //!for fluid nodes calculate pressure
						}
					} //!end i
				} //!end j
			} //!end k

			//--------------------//
			// inlet part level 1 //
			//--------------------//

			for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction in refined inlet region
			{
				for (int j = refLayer - 2; j <= jmax - refLayer + 2; j++) //!all voxels in y direction
				{
					for (int i = 0; i <= refLayerX - 3; i++) //!all voxels in x direction
					{
						if (lattice_[level][i][j][k] != NULL) //!write only existing voxels
						{
							paraview << lattice_[level][i][j][k]->calcDensity(bc_, alternatingL1) * speedOfSoundSq << "\n"; //!for fluid nodes calculate pressure
						}
					} //!end i
				} //!end j
			} //!end k

			//---------------------//
			// outlet part level 1 //
			//---------------------//

			for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction in refined outlet region
			{
				for (int j = refLayer - 2; j <= jmax - refLayer + 2; j++) //!all voxels in y direction
				{
					for (int i = imax - refLayer + 3; i <= imax; i++) //!all voxels in x direction
					{
						if (lattice_[level][i][j][k] != NULL) //!write only existing voxels
						{
							paraview << lattice_[level][i][j][k]->calcDensity(bc_, alternatingL1) * speedOfSoundSq << "\n"; //!for fluid nodes calculate pressure
						}
					} //!end i
				} //!end j
			} //!end k

			//------------------//
			// mid part level 0 //
			//------------------//

			level = 0;
			imax = bc_.getNumberOfVoxelsX() - 1; //!index of last node in x direction in ->this level
			jmax = bc_.getNumberOfVoxelsY() - 1; //!index of last node in y direction in ->this level
			kmax = bc_.getNumberOfVoxelsZ() - 1; //!index of last node in z direction in ->this level
			refLayer = bc_.getRefinementLayer();
			refLayerX = bc_.getRefinementLayerX() * pow(2, level);

			for (int k = refLayer - 1; k <= kmax - refLayer + 1; k++) //!all voxels in z direction in unrefined mid region
			{
				for (int j = refLayer - 1; j <= jmax - refLayer + 1; j++) //!all voxels in y direction
				{
					for (int i = refLayerX - 1; i <= imax - refLayer + 1; i++) //!all voxels in x direction
					{
						if (lattice_[level][i][j][k] != NULL) //!write only existing voxels
						{
							paraview << lattice_[level][i][j][k]->calcDensity(bc_, alternatingL0) * speedOfSoundSq << "\n"; //!for fluid and boundary voxels calc pressure
						}
					}//end i
				}//end j
			}//end k


			//->> DENSITY <<-//

			paraview << "SCALARS Density float 1" << "\n";
			paraview << "LOOKUP_TABLE default" << "\n";

			//---------------------//
			// bottom part level 1 //
			//---------------------//

			level = 1;
			imax = bc_.getNumberOfVoxelsX() * pow(2, level) - 1; //!index of last voxel in x direction in ->this level
			jmax = bc_.getNumberOfVoxelsY() * pow(2, level) - 1; //!index of last voxel in y direction in ->this level
			kmax = bc_.getNumberOfVoxelsZ() * pow(2, level) - 1; //!index of last voxel in z direction in ->this level
			refLayer = bc_.getRefinementLayer() * pow(2, level);
			refLayerX = bc_.getRefinementLayerX() * pow(2, level);

			for (int k = 0; k <= refLayer - 3; k++) //!all voxels in z direction in refined bottom region
			{
				for (int j = 0; j <= jmax; j++) //!all voxels in y direction
				{
					for (int i = 0; i <= imax; i++) //!all voxels in x direction
					{
						if (lattice_[level][i][j][k] != NULL) //!write only existing voxels
						{
							paraview << lattice_[level][i][j][k]->calcDensity(bc_, alternatingL1) << "\n"; //!for fluid nodes calculate density
						}
					} //!end i
				} //!end j
			} //!end k

			//------------------//
			// top part level 1 //
			//------------------//

			for (int k = kmax - refLayer + 3; k <= kmax; k++) //!all voxels in z direction in refined top region
			{
				for (int j = 0; j <= jmax; j++) //!all voxels in y direction
				{
					for (int i = 0; i <= imax; i++) //!all voxels in x direction
					{
						if (lattice_[level][i][j][k] != NULL) //!write only existing voxels
						{

							paraview << lattice_[level][i][j][k]->calcDensity(bc_, alternatingL1) << "\n"; //!for fluid nodes calculate density
						}
					} //!end i
				} //!end j
			} //!end k

			//-------------------//
			// left part level 1 //
			//-------------------//

			for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction in refined top region
			{
				for (int j = 0; j <= refLayer - 3; j++) //!all voxels in y direction
				{
					for (int i = 0; i <= imax; i++) //!all voxels in x direction
					{
						if (lattice_[level][i][j][k] != NULL) //!write only existing voxels
						{

							paraview << lattice_[level][i][j][k]->calcDensity(bc_, alternatingL1) << "\n"; //!for fluid nodes calculate density
						}
					} //!end i
				} //!end j
			} //!end k

			//--------------------//
			// right part level 1 //
			//--------------------//

			for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction in refined right region
			{
				for (int j = jmax - refLayer + 3; j <= jmax; j++) //!all voxels in y direction
				{
					for (int i = 0; i <= imax; i++) //!all voxels in x direction
					{
						if (lattice_[level][i][j][k] != NULL) //!write only existing voxels
						{
							paraview << lattice_[level][i][j][k]->calcDensity(bc_, alternatingL1) << "\n"; //!for fluid nodes calculate density
						}
					} //!end i
				} //!end j
			} //!end k

			//--------------------//
			// inlet part level 1 //
			//--------------------//

			for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction in refined inlet region
			{
				for (int j = refLayer - 2; j <= jmax - refLayer + 2; j++) //!all voxels in y direction
				{
					for (int i = 0; i <= refLayerX - 3; i++) //!all voxels in x direction
					{
						if (lattice_[level][i][j][k] != NULL) //!write only existing voxels
						{
							paraview << lattice_[level][i][j][k]->calcDensity(bc_, alternatingL1) << "\n"; //!for fluid nodes calculate density
						}
					} //!end i
				} //!end j
			} //!end k

			//---------------------//
			// outlet part level 1 //
			//---------------------//

			for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction in refined outlet region
			{
				for (int j = refLayer - 2; j <= jmax - refLayer + 2; j++) //!all voxels in y direction
				{
					for (int i = imax - refLayer + 3; i <= imax; i++) //!all voxels in x direction
					{
						if (lattice_[level][i][j][k] != NULL) //!write only existing voxels
						{
							paraview << lattice_[level][i][j][k]->calcDensity(bc_, alternatingL1) << "\n"; //!for fluid nodes calculate density

						}
					} //!end i
				} //!end j
			} //!end k

			//------------------//
			// mid part level 0 //
			//------------------//

			level = 0;
			imax = bc_.getNumberOfVoxelsX() - 1; //!index of last node in x direction in ->this level
			jmax = bc_.getNumberOfVoxelsY() - 1; //!index of last node in y direction in ->this level
			kmax = bc_.getNumberOfVoxelsZ() - 1; //!index of last node in z direction in ->this level
			refLayer = bc_.getRefinementLayer();
			refLayerX = bc_.getRefinementLayerX() * pow(2, level);

			for (int k = refLayer - 1; k <= kmax - refLayer + 1; k++) //!all voxels in z direction in unrefined mid region
			{
				for (int j = refLayer - 1; j <= jmax - refLayer + 1; j++) //!all voxels in y direction
				{
					for (int i = refLayerX - 1; i <= imax - refLayer + 1; i++) //!all voxels in x direction
					{
						if (lattice_[level][i][j][k] != NULL) //!write only existing voxels
						{
							paraview << lattice_[level][i][j][k]->calcDensity(bc_, alternatingL0) << "\n"; //!for fluid and boundary voxels calc density
						}
					}//end i
				}//end j
			}//end k


			//->> VELOCITY <<-//

			paraview << "\nVECTORS Velocity float " << "\n";

			//---------------------//
			// bottom part level 1 //
			//---------------------//

			level = 1;
			imax = bc_.getNumberOfVoxelsX() * pow(2, level) - 1; //!index of last voxel in x direction in ->this level
			jmax = bc_.getNumberOfVoxelsY() * pow(2, level) - 1; //!index of last voxel in y direction in ->this level
			kmax = bc_.getNumberOfVoxelsZ() * pow(2, level) - 1; //!index of last voxel in z direction in ->this level
			refLayer = bc_.getRefinementLayer() * pow(2, level);
			refLayerX = bc_.getRefinementLayerX() * pow(2, level);

			for (int k = 0; k <= refLayer - 3; k++) //!all voxels in z direction in refined bottom region
			{
				for (int j = 0; j <= jmax; j++) //!all voxels in y direction
				{
					for (int i = 0; i <= imax; i++) //!all voxels in x direction
					{
						if (lattice_[level][i][j][k] != NULL) //!write only existing voxels
						{
							double velo[3] = { 0. };
							double dens = lattice_[level][i][j][k]->calcDensity(bc_, alternatingL1);
							lattice_[level][i][j][k]->calcVelocity(bc_, alternatingL1, dens, velo);
							paraview << velo[0] << " " << velo[1] << " " << velo[2] << " " << "\n"; //!write velocity to file
						}
					} //!end i
				} //!end j
			} //!end k

			//------------------//
			// top part level 1 //
			//------------------//

			for (int k = kmax - refLayer + 3; k <= kmax; k++) //!all voxels in z direction in refined top region
			{
				for (int j = 0; j <= jmax; j++) //!all voxels in y direction
				{
					for (int i = 0; i <= imax; i++) //!all voxels in x direction
					{
						if (lattice_[level][i][j][k] != NULL) //!write only existing voxels
						{
							double velo[3] = { 0. };
							double dens = lattice_[level][i][j][k]->calcDensity(bc_, alternatingL1);
							lattice_[level][i][j][k]->calcVelocity(bc_, alternatingL1, dens, velo);
							paraview << velo[0] << " " << velo[1] << " " << velo[2] << " " << "\n"; //!write velocity to file
						}
					} //!end i
				} //!end j
			} //!end k

			//-------------------//
			// left part level 1 //
			//-------------------//

			for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction in refined top region
			{
				for (int j = 0; j <= refLayer - 3; j++) //!all voxels in y direction
				{
					for (int i = 0; i <= imax; i++) //!all voxels in x direction
					{
						if (lattice_[level][i][j][k] != NULL) //!write only existing voxels
						{
							double velo[3] = { 0. };
							double dens = lattice_[level][i][j][k]->calcDensity(bc_, alternatingL1);
							lattice_[level][i][j][k]->calcVelocity(bc_, alternatingL1, dens, velo);
							paraview << velo[0] << " " << velo[1] << " " << velo[2] << " " << "\n"; //!write velocity to file
						}
					} //!end i
				} //!end j
			} //!end k

			//--------------------//
			// right part level 1 //
			//--------------------//

			for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction in refined right region
			{
				for (int j = jmax - refLayer + 3; j <= jmax; j++) //!all voxels in y direction
				{
					for (int i = 0; i <= imax; i++) //!all voxels in x direction
					{
						if (lattice_[level][i][j][k] != NULL) //!write only existing voxels
						{
							double velo[3] = { 0. };
							double dens = lattice_[level][i][j][k]->calcDensity(bc_, alternatingL1);
							lattice_[level][i][j][k]->calcVelocity(bc_, alternatingL1, dens, velo);
							paraview << velo[0] << " " << velo[1] << " " << velo[2] << " " << "\n"; //!write velocity to file
						}
					} //!end i
				} //!end j
			} //!end k

			//--------------------//
			// inlet part level 1 //
			//--------------------//

			for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction in refined inlet region
			{
				for (int j = refLayer - 2; j <= jmax - refLayer + 2; j++) //!all voxels in y direction
				{
					for (int i = 0; i <= refLayerX - 3; i++) //!all voxels in x direction
					{
						if (lattice_[level][i][j][k] != NULL) //!write only existing voxels
						{
							double velo[3] = { 0. };
							double dens = lattice_[level][i][j][k]->calcDensity(bc_, alternatingL1);
							lattice_[level][i][j][k]->calcVelocity(bc_, alternatingL1, dens, velo);
							paraview << velo[0] << " " << velo[1] << " " << velo[2] << " " << "\n"; //!write velocity to file
						}
					} //!end i
				} //!end j
			} //!end k

			//---------------------//
			// outlet part level 1 //
			//---------------------//

			for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction in refined outlet region
			{
				for (int j = refLayer - 2; j <= jmax - refLayer + 2; j++) //!all voxels in y direction
				{
					for (int i = imax - refLayer + 3; i <= imax; i++) //!all voxels in x direction
					{
						if (lattice_[level][i][j][k] != NULL) //!write only existing voxels
						{
							double velo[3] = { 0. };
							double dens = lattice_[level][i][j][k]->calcDensity(bc_, alternatingL1);
							lattice_[level][i][j][k]->calcVelocity(bc_, alternatingL1, dens, velo);
							paraview << velo[0] << " " << velo[1] << " " << velo[2] << " " << "\n"; //!write velocity to file
						}
					} //!end i
				} //!end j
			} //!end k

			//------------------//
			// mid part level 0 //
			//------------------//

			level = 0;
			imax = bc_.getNumberOfVoxelsX() - 1; //!index of last voxel in x direction in ->this level
			jmax = bc_.getNumberOfVoxelsY() - 1; //!index of last voxel in y direction in ->this level
			kmax = bc_.getNumberOfVoxelsZ() - 1; //!index of last voxel in z direction in ->this level
			refLayer = bc_.getRefinementLayer();
			refLayerX = bc_.getRefinementLayerX() * pow(2, level);

			for (int k = refLayer - 1; k <= kmax - refLayer + 1; k++) //!all voxels in z direction in unrefined mid region
			{
				for (int j = refLayer - 1; j <= jmax - refLayer + 1; j++) //!all voxels in y direction
				{
					for (int i = refLayerX - 1; i <= imax - refLayer + 1; i++) //!all voxels in x direction
					{
						if (lattice_[level][i][j][k] != NULL) //!write voxels existing nodes
						{
							double velo[3] = { 0. };
							double dens = lattice_[level][i][j][k]->calcDensity(bc_, alternatingL0);
							lattice_[level][i][j][k]->calcVelocity(bc_, alternatingL0, dens, velo);
							paraview << velo[0] << " " << velo[1] << " " << velo[2] << " " << endl; //!write velocity to file
						}
					} //!end i
				} //!end j
			} //!end k

			//-------//
			// final //
			//-------//

			//!update node and cell numbers
			paraview.seekp(pntPos); //!go to nodePosition of the stream pointer
			paraview << pntCntr;	//!insert correct node number
			paraview.seekp(cellPos); //!go to cellPosition of the stream pointer
			paraview << cellCntr << " " << 9 * cellCntr; //!insert correct cell number and number of indices
			paraview.close(); //!close file
		}
		else
		{
			//!write all node coordinates to file
			int level = 0;
			int imax = bc_.getNumberOfNodesX() - 1; //!index of last node in x direction in ->this level
			int jmax = bc_.getNumberOfNodesY() - 1; //!index of last node in y direction in ->this level
			int kmax = bc_.getNumberOfNodesZ() - 1; //!index of last node in z direction in ->this level

			int nodeCntrX = imax + 1; //!nodes in x direction in this level
			int nodeCntrY = jmax + 1; //!nodes in y direction in this level
			int nodeCntrXY = nodeCntrX * nodeCntrY; //!nodes per xy plane in this level

			int index = 0;

			for (int k = 0; k <= kmax; k++) //!all nodes in z direction
			{
				for (int j = 0; j <= jmax; j++) //!all nodes in y direction
				{
					for (int i = 0; i <= imax; i++) //!all nodes in x direction
					{
						if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
						{
							paraview << mesh_[level][i][j][k]->getX() << " " << mesh_[level][i][j][k]->getY() << " " << mesh_[level][i][j][k]->getZ() << "\n"; //!write coordinates
							index++; //!count nodes
						}
					} //!end i
				} //!end j
			} //!end k

			paraview << "\n";

			pntCntr = index;	//!store numberOfPoints

			paraview << "CELLS "; //!write cells
			int cellPos = static_cast<int>(paraview.tellp()); //!remember cell position
			int cellCntr = 0;
			paraview << "         " << " " << "         " << "\n";

			int icmax = bc_.getNumberOfVoxelsX() - 1; //!index of last voxel in x direction in ->this level
			int jcmax = bc_.getNumberOfVoxelsY() - 1; //!index of last voxel in y direction in ->this level
			int kcmax = bc_.getNumberOfVoxelsZ() - 1; //!index of last voxel in z direction in ->this level

			//!write cell connectivity and write file; needed for visualization

			for (int k = 0; k <= kcmax; k++) //!all voxels in z direction
			{
				for (int j = 0; j <= jcmax; j++) //!all voxels in y direction
				{
					for (int i = 0; i <= icmax; i++) //!all voxels in x direction
					{
						if (lattice_[level][i][j][k] != NULL) //!write only existing voxels
						{
							//!write indices of nodes, who create the cell to file
							//!paraview cell connectivity (for type vtk_voxel (=11)) is:
							/*
							  6--------7
							 /|	      /|
							4-|------5 |
							| 2------|-3   z  y
							|/		 |/	    |/
							0--------1		-->x
							*/
							//!calc index of node 0 (i.e. line in paraview file) from i,j and k
							int nodeIndex0 = i + j * nodeCntrX + k * nodeCntrXY;

							//!write connectivity
							paraview << "8 "; //!write number of points in cell (i.e. number of variables to read in this line)
							paraview << nodeIndex0 << " " << nodeIndex0 + 1 << " " << nodeIndex0 + nodeCntrX << " " << nodeIndex0 + nodeCntrX + 1 << " "; //!nodes 0 to 3
							paraview << nodeIndex0 + nodeCntrXY << " " << nodeIndex0 + nodeCntrXY + 1 << " "; //!nodes 4 and 5
							paraview << nodeIndex0 + nodeCntrXY + nodeCntrX << " " << nodeIndex0 + nodeCntrXY + nodeCntrX + 1 << "\n"; //!nodes 6 and 7
							cellCntr++; //!count cell
						}
					} //!end i
				} //!end j
			} //!end k

			paraview << "\n";
			paraview << "CELL_TYPES " << cellCntr << "\n"; //!write cell types to file

			for (int i = 0; i < cellCntr; i++)		//!write cell types to paraview file
			{
				paraview << "11" << "\n";
			}
			paraview << "\n";

			//!write all node data to file

			//!first of all: pressure
			paraview << "CELL_DATA " << cellCntr << "\n";
			paraview << "SCALARS Pressure float 1" << "\n";
			paraview << "LOOKUP_TABLE default" << "\n";

			double speedOfSoundSq = pow(bc_.getSpeedOfSound(), 2.); //!square of speed of sound

			for (int k = 0; k <= kcmax; k++) //!all voxels in z direction
			{
				for (int j = 0; j <= jcmax; j++) //!all voxels in y direction
				{
					for (int i = 0; i <= icmax; i++) //!all voxels in x direction
					{
						if (lattice_[level][i][j][k] != NULL) //!write only existing voxels
						{
							paraview << lattice_[level][i][j][k]->calcDensity(bc_, alternatingL0) * speedOfSoundSq << "\n"; //!calculate pressure
						}
					} //!end i
				} //!end j
			} //!end k

			paraview << "SCALARS Density float 1" << "\n";
			paraview << "LOOKUP_TABLE default" << "\n";

			for (int k = 0; k <= kcmax; k++) //!all voxels in z direction
			{
				for (int j = 0; j <= jcmax; j++) //!all voxels in y direction
				{
					for (int i = 0; i <= icmax; i++) //!all voxels in x direction
					{
						if (lattice_[level][i][j][k] != NULL) //!write only existing voxels
						{
							paraview << lattice_[level][i][j][k]->calcDensity(bc_, alternatingL0) << "\n"; //!calculate density
						}
					} //!end i
				} //!end j
			} //!end k

			paraview << "SCALARS Omega float 1" << "\n";
			paraview << "LOOKUP_TABLE default" << "\n";

			for (int k = 0; k <= kcmax; k++) //!all voxels in z direction
			{
				for (int j = 0; j <= jcmax; j++) //!all voxels in y direction
				{
					for (int i = 0; i <= icmax; i++) //!all voxels in x direction
					{
						if (lattice_[level][i][j][k] != NULL) //!write only existing voxels
						{
							paraview << lattice_[level][i][j][k]->getCollisionFrequency() << "\n";
						}
					} //!end i
				} //!end j
			} //!end k

		//!next: write velocity to file
			paraview << "\nVECTORS Velocity float " << "\n";

			for (int k = 0; k <= kcmax; k++) //!all voxels in z direction
			{
				for (int j = 0; j <= jcmax; j++) //!all voxels in y direction
				{
					for (int i = 0; i <= icmax; i++) //!all voxels in x direction
					{
						if (lattice_[level][i][j][k] != NULL) //!write only existing voxels
						{
							double velo[3] = { 0. };
							double dens = lattice_[level][i][j][k]->calcDensity(bc_, alternatingL0);
							lattice_[level][i][j][k]->calcVelocity(bc_, alternatingL0, dens, velo);
							paraview << velo[0] << " " << velo[1] << " " << velo[2] << " " << "\n"; //!write velocity to file
						}
					} //!end i
				} //!end j
			} //!end k

		 //!update node and cell numbers
			paraview.seekp(pntPos);	//!go to nodePosition of the stream pointer
			paraview << pntCntr;	//!insert correct node number
			paraview.seekp(cellPos); //!go to cellPosition of the stream pointer
			paraview << cellCntr << " " << 9 * cellCntr;	//!insert correct cell number and number of indices
			paraview.close(); //!close file
		} //!end else (write to file)
	}
} //!end write results


//!ALGORITHM USED BY ROHDE
void Lattice::nestedTimeSteppingSRT(bool& alternatingL1, bool& alternatingL0, Control& bc_) //!perform nested time step with sequence of collide stream step for all level
{
	int refLayer, refLayerX; //!number of refined cells
	bool negAlternatingL1 = !alternatingL1;
	bool negAlternatingL0 = !alternatingL0;

	//!STRING VARIABLES FOR DIFFERENCIATING BETWEEN VARIOUS INTERFACE VOXELS -> NEEDED FOR EXPLODE AND COALESCE METHOD
	//!'core' interface voxels
	string interfaceBottom = "interfaceBottom";
	string interfaceTop = "interfaceTop";
	string interfaceLeft = "interfaceLeft";
	string interfaceRight = "interfaceRight";
	string interfaceInlet = "interfaceInlet";
	string interfaceOutlet = "interfaceOutlet";

	//!'edge' interface voxels
	string interfaceBottomLeft = "interfaceBottomLeft";
	string interfaceBottomRight = "interfaceBottomRight";
	string interfaceTopLeft = "interfaceTopLeft";
	string interfaceTopRight = "interfaceTopRight";
	string interfaceInletBottom = "interfaceInletBottom";
	string interfaceInletTop = "interfaceInletTop";
	string interfaceOutletBottom = "interfaceOutletBottom";
	string interfaceOutletTop = "interfaceOutletTop";
	string interfaceInletLeft = "interfaceInletLeft";
	string interfaceInletRight = "interfaceInletRight";
	string interfaceOutletLeft = "interfaceOutletLeft";
	string interfaceOutletRight = "interfaceOutletRight";

	//!'corner' interface voxels
	string interfaceInletBottomLeft = "interfaceInletBottomLeft";
	string interfaceInletBottomRight = "interfaceInletBottomRight";
	string interfaceInletTopLeft = "interfaceInletTopLeft";
	string interfaceInletTopRight = "interfaceInletTopRight";
	string interfaceOutletBottomLeft = "interfaceOutletBottomLeft";
	string interfaceOutletBottomRight = "interfaceOutletBottomRight";
	string interfaceOutletTopLeft = "interfaceOutletTopLeft";
	string interfaceOutletTopRight = "interfaceOutletTopRight";

	//--------------------------------------------------//
	//--------------------------------------------------//
	//             collision and transport				//
	//--------------------------------------------------//
	//--------------------------------------------------//

	//!current distribution (t=0) level 0: alternating
	//!current distribution (t=0) level 1: alternating

	int level = 0;
	int imax = bc_.getNumberOfVoxelsX() - 1; //!index of last voxel in x direction in ->this level
	int jmax = bc_.getNumberOfVoxelsY() - 1; //!index of last voxel in y direction in ->this level
	int kmax = bc_.getNumberOfVoxelsZ() - 1; //!index of last voxel in z direction in ->this level
	refLayer = bc_.getRefinementLayer(); //!number of refined coarse voxels
	refLayerX = bc_.getRefinementLayerX(); //!number of refined coarse voxels

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do outer loop in parallel
		for (int i = refLayerX - 1; i <= imax - refLayer + 1; i++) //!all nodes in x direction
		{
			for (int j = refLayer - 1; j <= jmax - refLayer + 1; j++) //!all nodes in y direction
			{
				for (int k = refLayer - 1; k <= kmax - refLayer + 1; k++) //!nodes in z direction from coarse lower interface to upper interface
				{
					for (int dir = 0; dir < 19; dir++)
					{
						double distrTemp = lattice_[level][i][j][k]->getDistribution(alternatingL0, dir); //!initialize temporary distribution
						lattice_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
					}
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

	//!collide in all coarse voxels
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do outer loop in parallel
		for (int i = refLayerX - 1; i <= imax - refLayer + 1; i++) //!all nodes in x direction
		{
			for (int j = refLayer - 1; j <= jmax - refLayer + 1; j++) //!all nodes in y direction
			{
				for (int k = refLayer - 1; k <= kmax - refLayer + 1; k++) //!nodes in z direction from coarse lower interface to upper interface
				{
					lattice_[level][i][j][k]->collideSRT(bc_, alternatingL0); //!perform collide in voxel with recursive-regularized collision operator (RR)
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

	level = 1;
	refLayer = bc_.getRefinementLayer() * pow(2, level);
	refLayerX = bc_.getRefinementLayerX() * pow(2, level);
	imax = bc_.getNumberOfVoxelsX() * pow(2., level) - 1; //!index of last voxel in x direction in ->this level
	jmax = bc_.getNumberOfVoxelsY() * pow(2., level) - 1; //!index of last voxel in y direction in ->this level
	kmax = bc_.getNumberOfVoxelsZ() * pow(2., level) - 1; //!index of last voxel in z direction in ->this level

	//!collide in inner fine voxels (no collision in fine interface voxels!)
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all voxels in y direction
			{
				//!bottom refined region
				for (int k = 0; k <= refLayer - 3; k++) //!all fine inner voxels in refined bottom region
				{
					lattice_[level][i][j][k]->collideSRT(bc_, alternatingL1); //!perform collide in voxel with recursive-regularized collision operator (RR)
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all voxels in y direction
			{
				//!top refined region
				for (int k = kmax - refLayer + 3; k <= kmax; k++) //!all fine inner voxels in refined top region
				{
					lattice_[level][i][j][k]->collideSRT(bc_, alternatingL1); //!perform collide in voxel with recursive-regularized collision operator (RR)
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction
			{
				//!left refined region
				for (int j = 0; j <= refLayer - 3; j++) //!all fine inner voxels in refined left region
				{
					lattice_[level][i][j][k]->collideSRT(bc_, alternatingL1); //!perform collide in voxel with recursive-regularized collision operator (RR)
				} //!end j
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction
			{
				//!right refined region
				for (int j = jmax - refLayer + 3; j <= jmax; j++) //!all fine inner voxels in refined right region
				{
					lattice_[level][i][j][k]->collideSRT(bc_, alternatingL1); //!perform collide in voxel with recursive-regularized collision operator (RR)
				} //!end j
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer - 2; j <= jmax - refLayer + 2; j++) //!all voxels in y direction
		{
			for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction
			{
				//!inlet refined region
				for (int i = 0; i <= refLayerX - 3; i++) //!all fine inner voxels in refined bottom region
				{
					lattice_[level][i][j][k]->collideSRT(bc_, alternatingL1); //!perform collide in voxel with recursive-regularized collision operator (RR)
				} //!end i
			} //!end k
		} //!end j
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer - 2; j <= jmax - refLayer + 2; j++) //!all voxels in y direction
		{
			for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction
			{
				//!outlet refined region
				for (int i = imax - refLayer + 3; i <= imax; i++) //!all fine inner voxels in refined top region
				{
					lattice_[level][i][j][k]->collideSRT(bc_, alternatingL1); //!perform collide in voxel with recursive-regularized collision operator (RR)
				} //!end i
			} //!end k
		} //!end j
	} //!end omp

	level = 0;
	imax = bc_.getNumberOfVoxelsX() - 1; //!index of last voxel in x direction in ->this level
	jmax = bc_.getNumberOfVoxelsY() - 1; //!index of last voxel in y direction in ->this level
	kmax = bc_.getNumberOfVoxelsZ() - 1; //!index of last voxel in z direction in ->this level
	refLayer = bc_.getRefinementLayer(); //!number of refined coarse voxels
	refLayerX = bc_.getRefinementLayerX() * pow(2, level);

	//!explode states in coarse interface voxels
	//!INTERFACE EDGES AND CORNERS INCLUDED (ALL DISTRIBUTIONS ARE BEING EXPLODED!)
	if (bc_.getExplosionOrder() == 1)
	{
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int i = refLayerX - 1; i <= imax - refLayer + 1; i++) //!all coarse interface voxels in x direction
			{
				for (int j = refLayer - 1; j <= jmax - refLayer + 1; j++) //!all coarse interface voxels in y direction
				{
					//!explode states in coarse bottom interface voxels
					int k = refLayer - 1;
					lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

					//!explode states in coarse top interface voxels
					k = kmax - refLayer + 1;
					lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int i = refLayerX - 1; i <= imax - refLayer + 1; i++) //!all coarse interface voxels in x direction
			{
				for (int k = refLayer; k <= kmax - refLayer; k++) //!all coarse interface voxels in z direction
				{
					//!explode states in coarse right interface voxels
					int j = refLayer - 1;
					lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

					//!explode states in coarse left interface voxels
					j = jmax - refLayer + 1;
					lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int j = refLayer; j <= jmax - refLayer; j++) //!all coarse interface voxels in y direction
			{
				for (int k = refLayer; k <= kmax - refLayer; k++) //!all coarse interface voxels in z direction
				{
					//!explode states in coarse inlet interface voxels
					int i = refLayerX - 1;
					lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

					//!explode states in coarse outlet interface voxels
					i = imax - refLayer + 1;
					lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);
				} //!end k
			} //!end i
		} //!end omp
	} //!end if

	else if (bc_.getExplosionOrder() == 2)
	{
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int i = refLayerX; i <= imax - refLayer; i++) //!all coarse interface voxels in x direction
			{
				for (int j = refLayer; j <= jmax - refLayer; j++) //!all coarse interface voxels in y direction
				{
					//!explode states in coarse bottom interface voxels
					int k = refLayer - 1;
					lattice_[level][i][j][k]->explodeLinear(alternatingL1, alternatingL0, bc_, interfaceBottom);

					//!explode states in coarse top interface voxels
					k = kmax - refLayer + 1;
					lattice_[level][i][j][k]->explodeLinear(alternatingL1, alternatingL0, bc_, interfaceTop);
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int i = refLayerX; i <= imax - refLayer; i++) //!all coarse interface voxels in x direction
			{
				for (int k = refLayer; k <= kmax - refLayer; k++) //!all coarse interface voxels in z direction
				{
					//!explode states in coarse right interface voxels
					int j = refLayer - 1;
					lattice_[level][i][j][k]->explodeLinear(alternatingL1, alternatingL0, bc_, interfaceLeft);

					//!explode states in coarse left interface voxels
					j = jmax - refLayer + 1;
					lattice_[level][i][j][k]->explodeLinear(alternatingL1, alternatingL0, bc_, interfaceRight);
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int j = refLayer; j <= jmax - refLayer; j++) //!all coarse interface voxels in y direction
			{
				for (int k = refLayer; k <= kmax - refLayer; k++) //!all coarse interface voxels in z direction
				{
					//!explode states in coarse inlet interface voxels
					int i = refLayerX - 1;
					lattice_[level][i][j][k]->explodeLinear(alternatingL1, alternatingL0, bc_, interfaceInlet);

					//!explode states in coarse outlet interface voxels
					i = imax - refLayer + 1;
					lattice_[level][i][j][k]->explodeLinear(alternatingL1, alternatingL0, bc_, interfaceOutlet);
				} //!end k
			} //!end i
		} //!end omp

		/*COARSE INTERFACE EDGES*/
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel;

			//!explode coarse states on interface edge voxels along x direction edges
			for (int i = refLayerX - 1; i <= imax - refLayer + 1; i++)
			{
				//!bottom left edge
				int j = refLayer - 1;
				int k = refLayer - 1;
				lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

				//!bottom right edge
				j = jmax - refLayer + 1;
				k = refLayer - 1;
				lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

				//!top left edge
				j = refLayer - 1;
				k = kmax - refLayer + 1;
				lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

				//!top right edge
				j = jmax - refLayer + 1;
				k = kmax - refLayer + 1;
				lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel;

			//!explode coarse states on interface edge voxels along y direction edges
			for (int j = refLayer; j <= jmax - refLayer; j++)
			{
				//!inlet bottom edge
				int i = refLayerX - 1;
				int k = refLayer - 1;
				lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

				//!inlet top edge
				i = refLayerX - 1;
				k = kmax - refLayer + 1;
				lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

				//!outlet bottom edge
				i = imax - refLayer + 1;
				k = refLayer - 1;
				lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

				//!outlet top edge
				i = imax - refLayer + 1;
				k = kmax - refLayer + 1;
				lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);
			} //!end j
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel;

			//!explode coarse states on interface edge voxels along z direction edges
			for (int k = refLayer; k <= kmax - refLayer; k++)
			{
				//!inlet left edge
				int i = refLayerX - 1;
				int j = refLayer - 1;
				lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

				//!inlet right edge
				i = refLayerX - 1;
				j = jmax - refLayer + 1;
				lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

				//!outlet left edge
				i = imax - refLayer + 1;
				j = refLayer - 1;
				lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

				//!outlet right edge
				i = imax - refLayer + 1;
				j = jmax - refLayer + 1;
				lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);
			} //!end k
		} //!end omp
	} //!end else if

	//!transport in all coarse voxels
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX - 1; i <= imax - refLayer + 1; i++) //!all voxels in x direction
		{
			for (int j = refLayer - 1; j <= jmax - refLayer + 1; j++) //!all voxels in y direction
			{
				for (int k = refLayer - 1; k <= kmax - refLayer + 1; k++)
				{
					lattice_[level][i][j][k]->transport(negAlternatingL0);
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

	//!transport in all fine voxels
	level = 1;
	refLayer = bc_.getRefinementLayer() * pow(2, level);
	refLayerX = bc_.getRefinementLayerX() * pow(2, level);
	imax = bc_.getNumberOfVoxelsX() * pow(2., level) - 1; //!index of last voxel in x direction in ->this level
	jmax = bc_.getNumberOfVoxelsY() * pow(2., level) - 1; //!index of last voxel in y direction in ->this level
	kmax = bc_.getNumberOfVoxelsZ() * pow(2., level) - 1; //!index of last voxel in z direction in ->this level

	//!transport in all fine voxels
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all voxels in y direction except walls
			{
				//!bottom region
				for (int k = 0; k <= refLayer - 1; k++) //!all fine voxels in refined bottom region except walls
				{
					lattice_[level][i][j][k]->transport(negAlternatingL1);
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all voxels in y direction except walls
			{
				//!top region
				for (int k = kmax - refLayer + 1; k <= kmax; k++) //!all fine voxels in refined top region except walls
				{
					lattice_[level][i][j][k]->transport(negAlternatingL1);
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int k = refLayer; k <= kmax - refLayer; k++) //!all voxels in z direction except walls and voxels considered before
			{
				//!left region
				for (int j = 0; j <= refLayer - 1; j++) //!all fine voxels in refined left region
				{
					lattice_[level][i][j][k]->transport(negAlternatingL1);
				} //!end k
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int k = refLayer; k <= kmax - refLayer; k++) //!all voxels in z direction except walls and voxels considered before
			{
				//!right region
				for (int j = jmax - refLayer + 1; j <= jmax; j++) //!all fine voxels in refined right region
				{
					lattice_[level][i][j][k]->transport(negAlternatingL1);
				} //!end j
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer; j <= jmax - refLayer; j++) //!all voxels in y direction except except walls and voxels considered before
		{
			for (int k = refLayer; k <= kmax - refLayer; k++) //!all voxels in z direction except walls and voxels considered before
			{
				//!inlet region
				for (int i = 0; i <= refLayerX - 1; i++) //!all fine voxels in refined inlet region except boundary
				{
					lattice_[level][i][j][k]->transport(negAlternatingL1);
				} //!end i
			} //!end k
		} //!end j
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer; j <= jmax - refLayer; j++) //!all voxels in y direction except except walls and voxels considered before
		{
			for (int k = refLayer; k <= kmax - refLayer; k++) //!all voxels in z direction except walls and voxels considered before
			{
				//!outlet region
				for (int i = imax - refLayer + 1; i <= imax; i++) //!all fine voxels in refined outlet region except boundary
				{
					lattice_[level][i][j][k]->transport(negAlternatingL1);
				} //!end i
			} //!end k
		} //!end j
	} //!end omp

	//!collide in inner fine voxels (no collision in fine interface voxels!)
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all voxels in y direction
			{
				//!bottom refined region
				for (int k = 0; k <= refLayer - 3; k++) //!all fine inner voxels in refined bottom region
				{
					lattice_[level][i][j][k]->collideSRT(bc_, negAlternatingL1); //!perform collide in voxel with recursive-regularized collision operator (RR)
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all voxels in y direction
			{
				//!top refined region
				for (int k = kmax - refLayer + 3; k <= kmax; k++) //!all fine inner voxels in refined top region
				{
					lattice_[level][i][j][k]->collideSRT(bc_, negAlternatingL1); //!perform collide in voxel with recursive-regularized collision operator (RR)
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction
			{
				//!left refined region
				for (int j = 0; j <= refLayer - 3; j++) //!all fine inner voxels in refined left region
				{
					lattice_[level][i][j][k]->collideSRT(bc_, negAlternatingL1); //!perform collide in voxel with recursive-regularized collision operator (RR)
				} //!end j
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction
			{
				//!right refined region
				for (int j = jmax - refLayer + 3; j <= jmax; j++) //!all fine inner voxels in refined right region
				{
					lattice_[level][i][j][k]->collideSRT(bc_, negAlternatingL1); //!perform collide in voxel with recursive-regularized collision operator (RR)
				} //!end j
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer - 2; j <= jmax - refLayer + 2; j++) //!all voxels in y direction
		{
			for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction
			{
				//!inlet refined region
				for (int i = 0; i <= refLayerX - 3; i++) //!all fine inner voxels in refined bottom region
				{
					lattice_[level][i][j][k]->collideSRT(bc_, negAlternatingL1); //!perform collide in voxel with recursive-regularized collision operator (RR)
				} //!end i
			} //!end k
		} //!end j
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer - 2; j <= jmax - refLayer + 2; j++) //!all voxels in y direction
		{
			for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction
			{
				//!outlet refined region
				for (int i = imax - refLayer + 3; i <= imax; i++) //!all fine inner voxels in refined top region
				{
					lattice_[level][i][j][k]->collideSRT(bc_, negAlternatingL1); //!perform collide in voxel with recursive-regularized collision operator (RR)
				} //!end i
			} //!end k
		} //!end j
	} //!end omp

//!transport in all fine voxels
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all voxels in y direction except walls
			{
				//!bottom region
				for (int k = 0; k <= refLayer - 1; k++) //!all fine voxels in refined bottom region except walls
				{
					lattice_[level][i][j][k]->transport(alternatingL1);
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all voxels in y direction except walls
			{
				//!top region
				for (int k = kmax - refLayer + 1; k <= kmax; k++) //!all fine voxels in refined top region except walls
				{
					lattice_[level][i][j][k]->transport(alternatingL1);
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int k = refLayer; k <= kmax - refLayer; k++) //!all voxels in z direction except walls and voxels considered before
			{
				//!left region
				for (int j = 0; j <= refLayer - 1; j++) //!all fine voxels in refined left region
				{
					lattice_[level][i][j][k]->transport(alternatingL1);
				} //!end k
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int k = refLayer; k <= kmax - refLayer; k++) //!all voxels in z direction except walls and voxels considered before
			{
				//!right region
				for (int j = jmax - refLayer + 1; j <= jmax; j++) //!all fine voxels in refined right region
				{
					lattice_[level][i][j][k]->transport(alternatingL1);
				} //!end j
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer; j <= jmax - refLayer; j++) //!all voxels in y direction except except walls and voxels considered before
		{
			for (int k = refLayer; k <= kmax - refLayer; k++) //!all voxels in z direction except walls and voxels considered before
			{
				//!inlet region
				for (int i = 0; i <= refLayerX - 1; i++) //!all fine voxels in refined inlet region except boundary
				{
					lattice_[level][i][j][k]->transport(alternatingL1);
				} //!end i
			} //!end k
		} //!end j
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer; j <= jmax - refLayer; j++) //!all voxels in y direction except except walls and voxels considered before
		{
			for (int k = refLayer; k <= kmax - refLayer; k++) //!all voxels in z direction except walls and voxels considered before
			{
				//!outlet region
				for (int i = imax - refLayer + 1; i <= imax; i++) //!all fine voxels in refined outlet region except boundary
				{
					lattice_[level][i][j][k]->transport(alternatingL1);
				} //!end i
			} //!end k
		} //!end j
	} //!end omp

	//!level 0
	level = 0;
	refLayer = bc_.getRefinementLayer();
	refLayerX = bc_.getRefinementLayerX() * pow(2, level);
	imax = bc_.getNumberOfVoxelsX() - 1; //!index of last voxel in x direction in ->this level
	jmax = bc_.getNumberOfVoxelsY() - 1; //!index of last voxel in y direction in ->this level
	kmax = bc_.getNumberOfVoxelsZ() - 1; //!index of last voxel in z direction in ->this level

	//!coalesce states in coarse interface voxels
	//!ATTENTION: COARSE EDGE AND CORNER VOXELS NEED SEPARATE TREATMENT
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do outer loop in parallel
		for (int i = refLayerX; i <= imax - refLayer; i++) //!all voxels in x direction
		{
			for (int j = refLayer; j <= jmax - refLayer; j++) //!all voxels in y direction
			{
				//!coalesce states in coarse bottom interface voxels
				int k = refLayer - 1;
				lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceBottom);

				//!coalesce states in coarse top interface voxels
				k = kmax - refLayer + 1;
				lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceTop);
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do outer loop in parallel
		for (int i = refLayerX; i <= imax - refLayer; i++) //!all coarse interface voxels in x direction
		{
			for (int k = refLayer; k <= kmax - refLayer; k++) //!all coarse interface voxels in z direction
			{
				//!coalesce states in coarse right interface voxels
				int j = refLayer - 1;
				lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceLeft);

				//!coalesce states in coarse left interface voxels
				j = jmax - refLayer + 1;
				lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceRight);
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do outer loop in parallel
		for (int j = refLayer; j <= jmax - refLayer; j++) //!all coarse interface voxels in y direction
		{
			for (int k = refLayer; k <= kmax - refLayer; k++) //!all coarse interface voxels in z direction
			{
				//!coalesce states in coarse inlet interface voxels
				int i = refLayerX - 1;
				lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceInlet);

				//!coalesce states in coarse outlet interface voxels
				i = imax - refLayer + 1;
				lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceOutlet);
			} //!end k
		} //!end i
	} //!end omp

	/*COARSE INTERFACE EDGES (NO CORNERS!)*/
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel;
		//!coalesce fine states in coarse interface edge voxels along x direction edges
		for (int i = refLayerX; i <= imax - refLayer; i++)
		{
			//!bottom left edge
			int j = refLayer - 1;
			int k = refLayer - 1;
			lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceBottomLeft);

			//!bottom right edge
			j = jmax - refLayer + 1;
			k = refLayer - 1;
			lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceBottomRight);

			//!top left edge
			j = refLayer - 1;
			k = kmax - refLayer + 1;
			lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceTopLeft);

			//!top right edge
			j = jmax - refLayer + 1;
			k = kmax - refLayer + 1;
			lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceTopRight);
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel;
		//!coalesce fine states in coarse interface edge voxels along y direction edges
		for (int j = refLayer; j <= jmax - refLayer; j++)
		{
			//!inlet bottom edge
			int i = refLayerX - 1;
			int k = refLayer - 1;
			lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceInletBottom);

			//!inlet top edge
			i = refLayerX - 1;
			k = kmax - refLayer + 1;
			lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceInletTop);

			//!outlet bottom edge
			i = imax - refLayer + 1;
			k = refLayer - 1;
			lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceOutletBottom);

			//!outlet top edge
			i = imax - refLayer + 1;
			k = kmax - refLayer + 1;
			lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceOutletTop);
		} //!end j
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel;
		//!coalesce fine states in coarse interface edge voxels along z direction edges
		for (int k = refLayer; k <= kmax - refLayer; k++)
		{
			//!inlet left edge
			int i = refLayerX - 1;
			int j = refLayer - 1;
			lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceInletLeft);

			//!inlet right edge
			i = refLayerX - 1;
			j = jmax - refLayer + 1;
			lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceInletRight);

			//!outlet left edge
			i = imax - refLayer + 1;
			j = refLayer - 1;
			lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceOutletLeft);

			//!outlet right edge
			i = imax - refLayer + 1;
			j = jmax - refLayer + 1;
			lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceOutletRight);
		} //!end k
	} //!end omp

	/*COARSE INTERFACE CORNERS*/
	//!inlet bottom left corner
	int i = refLayerX - 1;
	int j = refLayer - 1;
	int k = refLayer - 1;
	lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceInletBottomLeft);

	//!inlet bottom right corner
	i = refLayerX - 1;
	j = jmax - refLayer + 1;
	k = refLayer - 1;
	lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceInletBottomRight);

	//!inlet top left corner
	i = refLayerX - 1;
	j = refLayer - 1;
	k = kmax - refLayer + 1;
	lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceInletTopLeft);

	//!inlet top right corner
	i = refLayerX - 1;
	j = jmax - refLayer + 1;
	k = kmax - refLayer + 1;
	lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceInletTopRight);

	//!outlet bottom left corner
	i = imax - refLayer + 1;
	j = refLayer - 1;
	k = refLayer - 1;
	lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceOutletBottomLeft);

	//!outlet bottom right corner
	i = imax - refLayer + 1;
	j = jmax - refLayer + 1;
	k = refLayer - 1;
	lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceOutletBottomRight);

	//!outlet top left corner
	i = imax - refLayer + 1;
	j = refLayer - 1;
	k = kmax - refLayer + 1;
	lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceOutletTopLeft);

	//!outlet top right corner
	i = imax - refLayer + 1;
	j = jmax - refLayer + 1;
	k = kmax - refLayer + 1;
	lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceOutletTopRight);
} //!end timeStepNestedSRT


void Lattice::nestedTimeSteppingMRT_GSbasis(bool& alternatingL1, bool& alternatingL0, Control& bc_) //!perform nested time step with sequence of collide stream step for all level
{
	int refLayer, refLayerX; //!number of refined cells
	bool negAlternatingL1 = !alternatingL1;
	bool negAlternatingL0 = !alternatingL0;

	//!STRING VARIABLES FOR DIFFERENCIATING BETWEEN VARIOUS INTERFACE VOXELS -> NEEDED FOR EXPLODE AND COALESCE METHOD
	//!'core' interface voxels
	string interfaceBottom = "interfaceBottom";
	string interfaceTop = "interfaceTop";
	string interfaceLeft = "interfaceLeft";
	string interfaceRight = "interfaceRight";
	string interfaceInlet = "interfaceInlet";
	string interfaceOutlet = "interfaceOutlet";

	//!'edge' interface voxels
	string interfaceBottomLeft = "interfaceBottomLeft";
	string interfaceBottomRight = "interfaceBottomRight";
	string interfaceTopLeft = "interfaceTopLeft";
	string interfaceTopRight = "interfaceTopRight";
	string interfaceInletBottom = "interfaceInletBottom";
	string interfaceInletTop = "interfaceInletTop";
	string interfaceOutletBottom = "interfaceOutletBottom";
	string interfaceOutletTop = "interfaceOutletTop";
	string interfaceInletLeft = "interfaceInletLeft";
	string interfaceInletRight = "interfaceInletRight";
	string interfaceOutletLeft = "interfaceOutletLeft";
	string interfaceOutletRight = "interfaceOutletRight";

	//!'corner' interface voxels
	string interfaceInletBottomLeft = "interfaceInletBottomLeft";
	string interfaceInletBottomRight = "interfaceInletBottomRight";
	string interfaceInletTopLeft = "interfaceInletTopLeft";
	string interfaceInletTopRight = "interfaceInletTopRight";
	string interfaceOutletBottomLeft = "interfaceOutletBottomLeft";
	string interfaceOutletBottomRight = "interfaceOutletBottomRight";
	string interfaceOutletTopLeft = "interfaceOutletTopLeft";
	string interfaceOutletTopRight = "interfaceOutletTopRight";

	//--------------------------------------------------//
	//--------------------------------------------------//
	//             collision and transport				//
	//--------------------------------------------------//
	//--------------------------------------------------//

	//!current distribution (t=0) level 0: alternating
	//!current distribution (t=0) level 1: alternating

	int level = 0;
	int imax = bc_.getNumberOfVoxelsX() - 1; //!index of last voxel in x direction in ->this level
	int jmax = bc_.getNumberOfVoxelsY() - 1; //!index of last voxel in y direction in ->this level
	int kmax = bc_.getNumberOfVoxelsZ() - 1; //!index of last voxel in z direction in ->this level
	refLayer = bc_.getRefinementLayer(); //!number of refined coarse voxels
	refLayerX = bc_.getRefinementLayerX(); //!number of refined coarse voxels

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do outer loop in parallel
		for (int i = refLayerX - 1; i <= imax - refLayer + 1; i++) //!all nodes in x direction
		{
			for (int j = refLayer - 1; j <= jmax - refLayer + 1; j++) //!all nodes in y direction
			{
				for (int k = refLayer - 1; k <= kmax - refLayer + 1; k++) //!nodes in z direction from coarse lower interface to upper interface
				{
					for (int dir = 0; dir < 19; dir++)
					{
						double distrTemp = lattice_[level][i][j][k]->getDistribution(alternatingL0, dir); //!initialize temporary distribution
						lattice_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
					}
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

	//!collide in all coarse voxels
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do outer loop in parallel
		for (int i = refLayerX - 1; i <= imax - refLayer + 1; i++) //!all nodes in x direction
		{
			for (int j = refLayer - 1; j <= jmax - refLayer + 1; j++) //!all nodes in y direction
			{
				for (int k = refLayer - 1; k <= kmax - refLayer + 1; k++) //!nodes in z direction from coarse lower interface to upper interface
				{
					lattice_[level][i][j][k]->collideMRT_GSbasis(bc_, alternatingL0); //!perform collide in voxel with recursive-regularized collision operator (RR)
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

	level = 1;
	refLayer = bc_.getRefinementLayer() * pow(2, level);
	refLayerX = bc_.getRefinementLayerX() * pow(2, level);
	imax = bc_.getNumberOfVoxelsX() * pow(2., level) - 1; //!index of last voxel in x direction in ->this level
	jmax = bc_.getNumberOfVoxelsY() * pow(2., level) - 1; //!index of last voxel in y direction in ->this level
	kmax = bc_.getNumberOfVoxelsZ() * pow(2., level) - 1; //!index of last voxel in z direction in ->this level

	//!collide in inner fine voxels (no collision in fine interface voxels!)
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all voxels in y direction
			{
				//!bottom refined region
				for (int k = 0; k <= refLayer - 3; k++) //!all fine inner voxels in refined bottom region
				{
					lattice_[level][i][j][k]->collideMRT_GSbasis(bc_, alternatingL1); //!perform collide in voxel with recursive-regularized collision operator (RR)
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all voxels in y direction
			{
				//!top refined region
				for (int k = kmax - refLayer + 3; k <= kmax; k++) //!all fine inner voxels in refined top region
				{
					lattice_[level][i][j][k]->collideMRT_GSbasis(bc_, alternatingL1); //!perform collide in voxel with recursive-regularized collision operator (RR)
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction
			{
				//!left refined region
				for (int j = 0; j <= refLayer - 3; j++) //!all fine inner voxels in refined left region
				{
					lattice_[level][i][j][k]->collideMRT_GSbasis(bc_, alternatingL1); //!perform collide in voxel with recursive-regularized collision operator (RR)
				} //!end j
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction
			{
				//!right refined region
				for (int j = jmax - refLayer + 3; j <= jmax; j++) //!all fine inner voxels in refined right region
				{
					lattice_[level][i][j][k]->collideMRT_GSbasis(bc_, alternatingL1); //!perform collide in voxel with recursive-regularized collision operator (RR)
				} //!end j
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer - 2; j <= jmax - refLayer + 2; j++) //!all voxels in y direction
		{
			for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction
			{
				//!inlet refined region
				for (int i = 0; i <= refLayerX - 3; i++) //!all fine inner voxels in refined bottom region
				{
					lattice_[level][i][j][k]->collideMRT_GSbasis(bc_, alternatingL1); //!perform collide in voxel with recursive-regularized collision operator (RR)
				} //!end i
			} //!end k
		} //!end j
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer - 2; j <= jmax - refLayer + 2; j++) //!all voxels in y direction
		{
			for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction
			{
				//!outlet refined region
				for (int i = imax - refLayer + 3; i <= imax; i++) //!all fine inner voxels in refined top region
				{
					lattice_[level][i][j][k]->collideMRT_GSbasis(bc_, alternatingL1); //!perform collide in voxel with recursive-regularized collision operator (RR)
				} //!end i
			} //!end k
		} //!end j
	} //!end omp

	level = 0;
	imax = bc_.getNumberOfVoxelsX() - 1; //!index of last voxel in x direction in ->this level
	jmax = bc_.getNumberOfVoxelsY() - 1; //!index of last voxel in y direction in ->this level
	kmax = bc_.getNumberOfVoxelsZ() - 1; //!index of last voxel in z direction in ->this level
	refLayer = bc_.getRefinementLayer(); //!number of refined coarse voxels
	refLayerX = bc_.getRefinementLayerX() * pow(2, level);

	//!explode states in coarse interface voxels
	//!INTERFACE EDGES AND CORNERS INCLUDED (ALL DISTRIBUTIONS ARE BEING EXPLODED!)
	if (bc_.getExplosionOrder() == 1)
	{
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int i = refLayerX - 1; i <= imax - refLayer + 1; i++) //!all coarse interface voxels in x direction
			{
				for (int j = refLayer - 1; j <= jmax - refLayer + 1; j++) //!all coarse interface voxels in y direction
				{
					//!explode states in coarse bottom interface voxels
					int k = refLayer - 1;
					lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

					//!explode states in coarse top interface voxels
					k = kmax - refLayer + 1;
					lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int i = refLayerX - 1; i <= imax - refLayer + 1; i++) //!all coarse interface voxels in x direction
			{
				for (int k = refLayer; k <= kmax - refLayer; k++) //!all coarse interface voxels in z direction
				{
					//!explode states in coarse right interface voxels
					int j = refLayer - 1;
					lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

					//!explode states in coarse left interface voxels
					j = jmax - refLayer + 1;
					lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int j = refLayer; j <= jmax - refLayer; j++) //!all coarse interface voxels in y direction
			{
				for (int k = refLayer; k <= kmax - refLayer; k++) //!all coarse interface voxels in z direction
				{
					//!explode states in coarse inlet interface voxels
					int i = refLayerX - 1;
					lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

					//!explode states in coarse outlet interface voxels
					i = imax - refLayer + 1;
					lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);
				} //!end k
			} //!end i
		} //!end omp
	} //!end if

	else if (bc_.getExplosionOrder() == 2)
	{
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int i = refLayerX; i <= imax - refLayer; i++) //!all coarse interface voxels in x direction
			{
				for (int j = refLayer; j <= jmax - refLayer; j++) //!all coarse interface voxels in y direction
				{
					//!explode states in coarse bottom interface voxels
					int k = refLayer - 1;
					lattice_[level][i][j][k]->explodeLinear(alternatingL1, alternatingL0, bc_, interfaceBottom);

					//!explode states in coarse top interface voxels
					k = kmax - refLayer + 1;
					lattice_[level][i][j][k]->explodeLinear(alternatingL1, alternatingL0, bc_, interfaceTop);
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int i = refLayerX; i <= imax - refLayer; i++) //!all coarse interface voxels in x direction
			{
				for (int k = refLayer; k <= kmax - refLayer; k++) //!all coarse interface voxels in z direction
				{
					//!explode states in coarse right interface voxels
					int j = refLayer - 1;
					lattice_[level][i][j][k]->explodeLinear(alternatingL1, alternatingL0, bc_, interfaceLeft);

					//!explode states in coarse left interface voxels
					j = jmax - refLayer + 1;
					lattice_[level][i][j][k]->explodeLinear(alternatingL1, alternatingL0, bc_, interfaceRight);
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int j = refLayer; j <= jmax - refLayer; j++) //!all coarse interface voxels in y direction
			{
				for (int k = refLayer; k <= kmax - refLayer; k++) //!all coarse interface voxels in z direction
				{
					//!explode states in coarse inlet interface voxels
					int i = refLayerX - 1;
					lattice_[level][i][j][k]->explodeLinear(alternatingL1, alternatingL0, bc_, interfaceInlet);

					//!explode states in coarse outlet interface voxels
					i = imax - refLayer + 1;
					lattice_[level][i][j][k]->explodeLinear(alternatingL1, alternatingL0, bc_, interfaceOutlet);
				} //!end k
			} //!end i
		} //!end omp

		/*COARSE INTERFACE EDGES*/
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel;

			//!explode coarse states on interface edge voxels along x direction edges
			for (int i = refLayerX - 1; i <= imax - refLayer + 1; i++)
			{
				//!bottom left edge
				int j = refLayer - 1;
				int k = refLayer - 1;
				lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

				//!bottom right edge
				j = jmax - refLayer + 1;
				k = refLayer - 1;
				lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

				//!top left edge
				j = refLayer - 1;
				k = kmax - refLayer + 1;
				lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

				//!top right edge
				j = jmax - refLayer + 1;
				k = kmax - refLayer + 1;
				lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel;

			//!explode coarse states on interface edge voxels along y direction edges
			for (int j = refLayer; j <= jmax - refLayer; j++)
			{
				//!inlet bottom edge
				int i = refLayerX - 1;
				int k = refLayer - 1;
				lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

				//!inlet top edge
				i = refLayerX - 1;
				k = kmax - refLayer + 1;
				lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

				//!outlet bottom edge
				i = imax - refLayer + 1;
				k = refLayer - 1;
				lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

				//!outlet top edge
				i = imax - refLayer + 1;
				k = kmax - refLayer + 1;
				lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);
			} //!end j
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel;

			//!explode coarse states on interface edge voxels along z direction edges
			for (int k = refLayer; k <= kmax - refLayer; k++)
			{
				//!inlet left edge
				int i = refLayerX - 1;
				int j = refLayer - 1;
				lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

				//!inlet right edge
				i = refLayerX - 1;
				j = jmax - refLayer + 1;
				lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

				//!outlet left edge
				i = imax - refLayer + 1;
				j = refLayer - 1;
				lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

				//!outlet right edge
				i = imax - refLayer + 1;
				j = jmax - refLayer + 1;
				lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);
			} //!end k
		} //!end omp
	} //!end else if

	//!transport in all coarse voxels
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX - 1; i <= imax - refLayer + 1; i++) //!all voxels in x direction
		{
			for (int j = refLayer - 1; j <= jmax - refLayer + 1; j++) //!all voxels in y direction
			{
				for (int k = refLayer - 1; k <= kmax - refLayer + 1; k++)
				{
					lattice_[level][i][j][k]->transport(negAlternatingL0);
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

	//!transport in all fine voxels
	level = 1;
	refLayer = bc_.getRefinementLayer() * pow(2, level);
	refLayerX = bc_.getRefinementLayerX() * pow(2, level);
	imax = bc_.getNumberOfVoxelsX() * pow(2., level) - 1; //!index of last voxel in x direction in ->this level
	jmax = bc_.getNumberOfVoxelsY() * pow(2., level) - 1; //!index of last voxel in y direction in ->this level
	kmax = bc_.getNumberOfVoxelsZ() * pow(2., level) - 1; //!index of last voxel in z direction in ->this level

	//!transport in all fine voxels
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all voxels in y direction except walls
			{
				//!bottom region
				for (int k = 0; k <= refLayer - 1; k++) //!all fine voxels in refined bottom region except walls
				{
					lattice_[level][i][j][k]->transport(negAlternatingL1);
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all voxels in y direction except walls
			{
				//!top region
				for (int k = kmax - refLayer + 1; k <= kmax; k++) //!all fine voxels in refined top region except walls
				{
					lattice_[level][i][j][k]->transport(negAlternatingL1);
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int k = refLayer; k <= kmax - refLayer; k++) //!all voxels in z direction except walls and voxels considered before
			{
				//!left region
				for (int j = 0; j <= refLayer - 1; j++) //!all fine voxels in refined left region
				{
					lattice_[level][i][j][k]->transport(negAlternatingL1);
				} //!end k
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int k = refLayer; k <= kmax - refLayer; k++) //!all voxels in z direction except walls and voxels considered before
			{
				//!right region
				for (int j = jmax - refLayer + 1; j <= jmax; j++) //!all fine voxels in refined right region
				{
					lattice_[level][i][j][k]->transport(negAlternatingL1);
				} //!end j
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer; j <= jmax - refLayer; j++) //!all voxels in y direction except except walls and voxels considered before
		{
			for (int k = refLayer; k <= kmax - refLayer; k++) //!all voxels in z direction except walls and voxels considered before
			{
				//!inlet region
				for (int i = 0; i <= refLayerX - 1; i++) //!all fine voxels in refined inlet region except boundary
				{
					lattice_[level][i][j][k]->transport(negAlternatingL1);
				} //!end i
			} //!end k
		} //!end j
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer; j <= jmax - refLayer; j++) //!all voxels in y direction except except walls and voxels considered before
		{
			for (int k = refLayer; k <= kmax - refLayer; k++) //!all voxels in z direction except walls and voxels considered before
			{
				//!outlet region
				for (int i = imax - refLayer + 1; i <= imax; i++) //!all fine voxels in refined outlet region except boundary
				{
					lattice_[level][i][j][k]->transport(negAlternatingL1);
				} //!end i
			} //!end k
		} //!end j
	} //!end omp

	//!collide in inner fine voxels (no collision in fine interface voxels!)
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all voxels in y direction
			{
				//!bottom refined region
				for (int k = 0; k <= refLayer - 3; k++) //!all fine inner voxels in refined bottom region
				{
					lattice_[level][i][j][k]->collideMRT_GSbasis(bc_, negAlternatingL1); //!perform collide in voxel with recursive-regularized collision operator (RR)
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all voxels in y direction
			{
				//!top refined region
				for (int k = kmax - refLayer + 3; k <= kmax; k++) //!all fine inner voxels in refined top region
				{
					lattice_[level][i][j][k]->collideMRT_GSbasis(bc_, negAlternatingL1); //!perform collide in voxel with recursive-regularized collision operator (RR)
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction
			{
				//!left refined region
				for (int j = 0; j <= refLayer - 3; j++) //!all fine inner voxels in refined left region
				{
					lattice_[level][i][j][k]->collideMRT_GSbasis(bc_, negAlternatingL1); //!perform collide in voxel with recursive-regularized collision operator (RR)
				} //!end j
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction
			{
				//!right refined region
				for (int j = jmax - refLayer + 3; j <= jmax; j++) //!all fine inner voxels in refined right region
				{
					lattice_[level][i][j][k]->collideMRT_GSbasis(bc_, negAlternatingL1); //!perform collide in voxel with recursive-regularized collision operator (RR)
				} //!end j
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer - 2; j <= jmax - refLayer + 2; j++) //!all voxels in y direction
		{
			for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction
			{
				//!inlet refined region
				for (int i = 0; i <= refLayerX - 3; i++) //!all fine inner voxels in refined bottom region
				{
					lattice_[level][i][j][k]->collideMRT_GSbasis(bc_, negAlternatingL1); //!perform collide in voxel with recursive-regularized collision operator (RR)
				} //!end i
			} //!end k
		} //!end j
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer - 2; j <= jmax - refLayer + 2; j++) //!all voxels in y direction
		{
			for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction
			{
				//!outlet refined region
				for (int i = imax - refLayer + 3; i <= imax; i++) //!all fine inner voxels in refined top region
				{
					lattice_[level][i][j][k]->collideMRT_GSbasis(bc_, negAlternatingL1); //!perform collide in voxel with recursive-regularized collision operator (RR)
				} //!end i
			} //!end k
		} //!end j
	} //!end omp

//!transport in all fine voxels
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all voxels in y direction except walls
			{
				//!bottom region
				for (int k = 0; k <= refLayer - 1; k++) //!all fine voxels in refined bottom region except walls
				{
					lattice_[level][i][j][k]->transport(alternatingL1);
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all voxels in y direction except walls
			{
				//!top region
				for (int k = kmax - refLayer + 1; k <= kmax; k++) //!all fine voxels in refined top region except walls
				{
					lattice_[level][i][j][k]->transport(alternatingL1);
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int k = refLayer; k <= kmax - refLayer; k++) //!all voxels in z direction except walls and voxels considered before
			{
				//!left region
				for (int j = 0; j <= refLayer - 1; j++) //!all fine voxels in refined left region
				{
					lattice_[level][i][j][k]->transport(alternatingL1);
				} //!end k
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int k = refLayer; k <= kmax - refLayer; k++) //!all voxels in z direction except walls and voxels considered before
			{
				//!right region
				for (int j = jmax - refLayer + 1; j <= jmax; j++) //!all fine voxels in refined right region
				{
					lattice_[level][i][j][k]->transport(alternatingL1);
				} //!end j
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer; j <= jmax - refLayer; j++) //!all voxels in y direction except except walls and voxels considered before
		{
			for (int k = refLayer; k <= kmax - refLayer; k++) //!all voxels in z direction except walls and voxels considered before
			{
				//!inlet region
				for (int i = 0; i <= refLayerX - 1; i++) //!all fine voxels in refined inlet region except boundary
				{
					lattice_[level][i][j][k]->transport(alternatingL1);
				} //!end i
			} //!end k
		} //!end j
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer; j <= jmax - refLayer; j++) //!all voxels in y direction except except walls and voxels considered before
		{
			for (int k = refLayer; k <= kmax - refLayer; k++) //!all voxels in z direction except walls and voxels considered before
			{
				//!outlet region
				for (int i = imax - refLayer + 1; i <= imax; i++) //!all fine voxels in refined outlet region except boundary
				{
					lattice_[level][i][j][k]->transport(alternatingL1);
				} //!end i
			} //!end k
		} //!end j
	} //!end omp

	//!level 0
	level = 0;
	refLayer = bc_.getRefinementLayer();
	refLayerX = bc_.getRefinementLayerX() * pow(2, level);
	imax = bc_.getNumberOfVoxelsX() - 1; //!index of last voxel in x direction in ->this level
	jmax = bc_.getNumberOfVoxelsY() - 1; //!index of last voxel in y direction in ->this level
	kmax = bc_.getNumberOfVoxelsZ() - 1; //!index of last voxel in z direction in ->this level

	//!coalesce states in coarse interface voxels
	//!ATTENTION: COARSE EDGE AND CORNER VOXELS NEED SEPARATE TREATMENT
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do outer loop in parallel
		for (int i = refLayerX; i <= imax - refLayer; i++) //!all voxels in x direction
		{
			for (int j = refLayer; j <= jmax - refLayer; j++) //!all voxels in y direction
			{
				//!coalesce states in coarse bottom interface voxels
				int k = refLayer - 1;
				lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceBottom);

				//!coalesce states in coarse top interface voxels
				k = kmax - refLayer + 1;
				lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceTop);
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do outer loop in parallel
		for (int i = refLayerX; i <= imax - refLayer; i++) //!all coarse interface voxels in x direction
		{
			for (int k = refLayer; k <= kmax - refLayer; k++) //!all coarse interface voxels in z direction
			{
				//!coalesce states in coarse right interface voxels
				int j = refLayer - 1;
				lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceLeft);

				//!coalesce states in coarse left interface voxels
				j = jmax - refLayer + 1;
				lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceRight);
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do outer loop in parallel
		for (int j = refLayer; j <= jmax - refLayer; j++) //!all coarse interface voxels in y direction
		{
			for (int k = refLayer; k <= kmax - refLayer; k++) //!all coarse interface voxels in z direction
			{
				//!coalesce states in coarse inlet interface voxels
				int i = refLayerX - 1;
				lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceInlet);

				//!coalesce states in coarse outlet interface voxels
				i = imax - refLayer + 1;
				lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceOutlet);
			} //!end k
		} //!end i
	} //!end omp

	/*COARSE INTERFACE EDGES (NO CORNERS!)*/
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel;
		//!coalesce fine states in coarse interface edge voxels along x direction edges
		for (int i = refLayerX; i <= imax - refLayer; i++)
		{
			//!bottom left edge
			int j = refLayer - 1;
			int k = refLayer - 1;
			lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceBottomLeft);

			//!bottom right edge
			j = jmax - refLayer + 1;
			k = refLayer - 1;
			lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceBottomRight);

			//!top left edge
			j = refLayer - 1;
			k = kmax - refLayer + 1;
			lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceTopLeft);

			//!top right edge
			j = jmax - refLayer + 1;
			k = kmax - refLayer + 1;
			lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceTopRight);
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel;
		//!coalesce fine states in coarse interface edge voxels along y direction edges
		for (int j = refLayer; j <= jmax - refLayer; j++)
		{
			//!inlet bottom edge
			int i = refLayerX - 1;
			int k = refLayer - 1;
			lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceInletBottom);

			//!inlet top edge
			i = refLayerX - 1;
			k = kmax - refLayer + 1;
			lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceInletTop);

			//!outlet bottom edge
			i = imax - refLayer + 1;
			k = refLayer - 1;
			lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceOutletBottom);

			//!outlet top edge
			i = imax - refLayer + 1;
			k = kmax - refLayer + 1;
			lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceOutletTop);
		} //!end j
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel;
		//!coalesce fine states in coarse interface edge voxels along z direction edges
		for (int k = refLayer; k <= kmax - refLayer; k++)
		{
			//!inlet left edge
			int i = refLayerX - 1;
			int j = refLayer - 1;
			lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceInletLeft);

			//!inlet right edge
			i = refLayerX - 1;
			j = jmax - refLayer + 1;
			lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceInletRight);

			//!outlet left edge
			i = imax - refLayer + 1;
			j = refLayer - 1;
			lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceOutletLeft);

			//!outlet right edge
			i = imax - refLayer + 1;
			j = jmax - refLayer + 1;
			lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceOutletRight);
		} //!end k
	} //!end omp

	/*COARSE INTERFACE CORNERS*/
	//!inlet bottom left corner
	int i = refLayerX - 1;
	int j = refLayer - 1;
	int k = refLayer - 1;
	lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceInletBottomLeft);

	//!inlet bottom right corner
	i = refLayerX - 1;
	j = jmax - refLayer + 1;
	k = refLayer - 1;
	lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceInletBottomRight);

	//!inlet top left corner
	i = refLayerX - 1;
	j = refLayer - 1;
	k = kmax - refLayer + 1;
	lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceInletTopLeft);

	//!inlet top right corner
	i = refLayerX - 1;
	j = jmax - refLayer + 1;
	k = kmax - refLayer + 1;
	lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceInletTopRight);

	//!outlet bottom left corner
	i = imax - refLayer + 1;
	j = refLayer - 1;
	k = refLayer - 1;
	lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceOutletBottomLeft);

	//!outlet bottom right corner
	i = imax - refLayer + 1;
	j = jmax - refLayer + 1;
	k = refLayer - 1;
	lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceOutletBottomRight);

	//!outlet top left corner
	i = imax - refLayer + 1;
	j = refLayer - 1;
	k = kmax - refLayer + 1;
	lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceOutletTopLeft);

	//!outlet top right corner
	i = imax - refLayer + 1;
	j = jmax - refLayer + 1;
	k = kmax - refLayer + 1;
	lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceOutletTopRight);
} //!end timeStepNestedMRT_GSbasis


void Lattice::nestedTimeSteppingMRT_RMbasis(bool& alternatingL1, bool& alternatingL0, Control& bc_) //!perform nested time step with sequence of collide stream step for all level
{
	int refLayer, refLayerX; //!number of refined cells
	bool negAlternatingL1 = !alternatingL1;
	bool negAlternatingL0 = !alternatingL0;

	//!STRING VARIABLES FOR DIFFERENCIATING BETWEEN VARIOUS INTERFACE VOXELS -> NEEDED FOR EXPLODE AND COALESCE METHOD
	//!'core' interface voxels
	string interfaceBottom = "interfaceBottom";
	string interfaceTop = "interfaceTop";
	string interfaceLeft = "interfaceLeft";
	string interfaceRight = "interfaceRight";
	string interfaceInlet = "interfaceInlet";
	string interfaceOutlet = "interfaceOutlet";

	//!'edge' interface voxels
	string interfaceBottomLeft = "interfaceBottomLeft";
	string interfaceBottomRight = "interfaceBottomRight";
	string interfaceTopLeft = "interfaceTopLeft";
	string interfaceTopRight = "interfaceTopRight";
	string interfaceInletBottom = "interfaceInletBottom";
	string interfaceInletTop = "interfaceInletTop";
	string interfaceOutletBottom = "interfaceOutletBottom";
	string interfaceOutletTop = "interfaceOutletTop";
	string interfaceInletLeft = "interfaceInletLeft";
	string interfaceInletRight = "interfaceInletRight";
	string interfaceOutletLeft = "interfaceOutletLeft";
	string interfaceOutletRight = "interfaceOutletRight";

	//!'corner' interface voxels
	string interfaceInletBottomLeft = "interfaceInletBottomLeft";
	string interfaceInletBottomRight = "interfaceInletBottomRight";
	string interfaceInletTopLeft = "interfaceInletTopLeft";
	string interfaceInletTopRight = "interfaceInletTopRight";
	string interfaceOutletBottomLeft = "interfaceOutletBottomLeft";
	string interfaceOutletBottomRight = "interfaceOutletBottomRight";
	string interfaceOutletTopLeft = "interfaceOutletTopLeft";
	string interfaceOutletTopRight = "interfaceOutletTopRight";

	//--------------------------------------------------//
	//--------------------------------------------------//
	//             collision and transport				//
	//--------------------------------------------------//
	//--------------------------------------------------//

	//!current distribution (t=0) level 0: alternating
	//!current distribution (t=0) level 1: alternating

	int level = 0;
	int imax = bc_.getNumberOfVoxelsX() - 1; //!index of last voxel in x direction in ->this level
	int jmax = bc_.getNumberOfVoxelsY() - 1; //!index of last voxel in y direction in ->this level
	int kmax = bc_.getNumberOfVoxelsZ() - 1; //!index of last voxel in z direction in ->this level
	refLayer = bc_.getRefinementLayer(); //!number of refined coarse voxels
	refLayerX = bc_.getRefinementLayerX(); //!number of refined coarse voxels

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do outer loop in parallel
		for (int i = refLayerX - 1; i <= imax - refLayer + 1; i++) //!all nodes in x direction
		{
			for (int j = refLayer - 1; j <= jmax - refLayer + 1; j++) //!all nodes in y direction
			{
				for (int k = refLayer - 1; k <= kmax - refLayer + 1; k++) //!nodes in z direction from coarse lower interface to upper interface
				{
					for (int dir = 0; dir < 19; dir++)
					{
						double distrTemp = lattice_[level][i][j][k]->getDistribution(alternatingL0, dir); //!initialize temporary distribution
						lattice_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
					}
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

	//!collide in all coarse voxels
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do outer loop in parallel
		for (int i = refLayerX - 1; i <= imax - refLayer + 1; i++) //!all nodes in x direction
		{
			for (int j = refLayer - 1; j <= jmax - refLayer + 1; j++) //!all nodes in y direction
			{
				for (int k = refLayer - 1; k <= kmax - refLayer + 1; k++) //!nodes in z direction from coarse lower interface to upper interface
				{
					lattice_[level][i][j][k]->collideMRT_RMbasis(bc_, alternatingL0); //!perform collide in voxel with recursive-regularized collision operator (RR)
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

	level = 1;
	refLayer = bc_.getRefinementLayer() * pow(2, level);
	refLayerX = bc_.getRefinementLayerX() * pow(2, level);
	imax = bc_.getNumberOfVoxelsX() * pow(2., level) - 1; //!index of last voxel in x direction in ->this level
	jmax = bc_.getNumberOfVoxelsY() * pow(2., level) - 1; //!index of last voxel in y direction in ->this level
	kmax = bc_.getNumberOfVoxelsZ() * pow(2., level) - 1; //!index of last voxel in z direction in ->this level

	//!collide in inner fine voxels (no collision in fine interface voxels!)
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all voxels in y direction
			{
				//!bottom refined region
				for (int k = 0; k <= refLayer - 3; k++) //!all fine inner voxels in refined bottom region
				{
					lattice_[level][i][j][k]->collideMRT_RMbasis(bc_, alternatingL1); //!perform collide in voxel with recursive-regularized collision operator (RR)
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all voxels in y direction
			{
				//!top refined region
				for (int k = kmax - refLayer + 3; k <= kmax; k++) //!all fine inner voxels in refined top region
				{
					lattice_[level][i][j][k]->collideMRT_RMbasis(bc_, alternatingL1); //!perform collide in voxel with recursive-regularized collision operator (RR)
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction
			{
				//!left refined region
				for (int j = 0; j <= refLayer - 3; j++) //!all fine inner voxels in refined left region
				{
					lattice_[level][i][j][k]->collideMRT_RMbasis(bc_, alternatingL1); //!perform collide in voxel with recursive-regularized collision operator (RR)
				} //!end j
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction
			{
				//!right refined region
				for (int j = jmax - refLayer + 3; j <= jmax; j++) //!all fine inner voxels in refined right region
				{
					lattice_[level][i][j][k]->collideMRT_RMbasis(bc_, alternatingL1); //!perform collide in voxel with recursive-regularized collision operator (RR)
				} //!end j
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer - 2; j <= jmax - refLayer + 2; j++) //!all voxels in y direction
		{
			for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction
			{
				//!inlet refined region
				for (int i = 0; i <= refLayerX - 3; i++) //!all fine inner voxels in refined bottom region
				{
					lattice_[level][i][j][k]->collideMRT_RMbasis(bc_, alternatingL1); //!perform collide in voxel with recursive-regularized collision operator (RR)
				} //!end i
			} //!end k
		} //!end j
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer - 2; j <= jmax - refLayer + 2; j++) //!all voxels in y direction
		{
			for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction
			{
				//!outlet refined region
				for (int i = imax - refLayer + 3; i <= imax; i++) //!all fine inner voxels in refined top region
				{
					lattice_[level][i][j][k]->collideMRT_RMbasis(bc_, alternatingL1); //!perform collide in voxel with recursive-regularized collision operator (RR)
				} //!end i
			} //!end k
		} //!end j
	} //!end omp

	level = 0;
	imax = bc_.getNumberOfVoxelsX() - 1; //!index of last voxel in x direction in ->this level
	jmax = bc_.getNumberOfVoxelsY() - 1; //!index of last voxel in y direction in ->this level
	kmax = bc_.getNumberOfVoxelsZ() - 1; //!index of last voxel in z direction in ->this level
	refLayer = bc_.getRefinementLayer(); //!number of refined coarse voxels
	refLayerX = bc_.getRefinementLayerX() * pow(2, level);

	//!explode states in coarse interface voxels
	//!INTERFACE EDGES AND CORNERS INCLUDED (ALL DISTRIBUTIONS ARE BEING EXPLODED!)
	if (bc_.getExplosionOrder() == 1)
	{
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int i = refLayerX - 1; i <= imax - refLayer + 1; i++) //!all coarse interface voxels in x direction
			{
				for (int j = refLayer - 1; j <= jmax - refLayer + 1; j++) //!all coarse interface voxels in y direction
				{
					//!explode states in coarse bottom interface voxels
					int k = refLayer - 1;
					lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

					//!explode states in coarse top interface voxels
					k = kmax - refLayer + 1;
					lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int i = refLayerX - 1; i <= imax - refLayer + 1; i++) //!all coarse interface voxels in x direction
			{
				for (int k = refLayer; k <= kmax - refLayer; k++) //!all coarse interface voxels in z direction
				{
					//!explode states in coarse right interface voxels
					int j = refLayer - 1;
					lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

					//!explode states in coarse left interface voxels
					j = jmax - refLayer + 1;
					lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int j = refLayer; j <= jmax - refLayer; j++) //!all coarse interface voxels in y direction
			{
				for (int k = refLayer; k <= kmax - refLayer; k++) //!all coarse interface voxels in z direction
				{
					//!explode states in coarse inlet interface voxels
					int i = refLayerX - 1;
					lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

					//!explode states in coarse outlet interface voxels
					i = imax - refLayer + 1;
					lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);
				} //!end k
			} //!end i
		} //!end omp
	} //!end if

	else if (bc_.getExplosionOrder() == 2)
	{
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int i = refLayerX; i <= imax - refLayer; i++) //!all coarse interface voxels in x direction
			{
				for (int j = refLayer; j <= jmax - refLayer; j++) //!all coarse interface voxels in y direction
				{
					//!explode states in coarse bottom interface voxels
					int k = refLayer - 1;
					lattice_[level][i][j][k]->explodeLinear(alternatingL1, alternatingL0, bc_, interfaceBottom);

					//!explode states in coarse top interface voxels
					k = kmax - refLayer + 1;
					lattice_[level][i][j][k]->explodeLinear(alternatingL1, alternatingL0, bc_, interfaceTop);
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int i = refLayerX; i <= imax - refLayer; i++) //!all coarse interface voxels in x direction
			{
				for (int k = refLayer; k <= kmax - refLayer; k++) //!all coarse interface voxels in z direction
				{
					//!explode states in coarse right interface voxels
					int j = refLayer - 1;
					lattice_[level][i][j][k]->explodeLinear(alternatingL1, alternatingL0, bc_, interfaceLeft);

					//!explode states in coarse left interface voxels
					j = jmax - refLayer + 1;
					lattice_[level][i][j][k]->explodeLinear(alternatingL1, alternatingL0, bc_, interfaceRight);
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int j = refLayer; j <= jmax - refLayer; j++) //!all coarse interface voxels in y direction
			{
				for (int k = refLayer; k <= kmax - refLayer; k++) //!all coarse interface voxels in z direction
				{
					//!explode states in coarse inlet interface voxels
					int i = refLayerX - 1;
					lattice_[level][i][j][k]->explodeLinear(alternatingL1, alternatingL0, bc_, interfaceInlet);

					//!explode states in coarse outlet interface voxels
					i = imax - refLayer + 1;
					lattice_[level][i][j][k]->explodeLinear(alternatingL1, alternatingL0, bc_, interfaceOutlet);
				} //!end k
			} //!end i
		} //!end omp

		/*COARSE INTERFACE EDGES*/
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel;

			//!explode coarse states on interface edge voxels along x direction edges
			for (int i = refLayerX - 1; i <= imax - refLayer + 1; i++)
			{
				//!bottom left edge
				int j = refLayer - 1;
				int k = refLayer - 1;
				lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

				//!bottom right edge
				j = jmax - refLayer + 1;
				k = refLayer - 1;
				lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

				//!top left edge
				j = refLayer - 1;
				k = kmax - refLayer + 1;
				lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

				//!top right edge
				j = jmax - refLayer + 1;
				k = kmax - refLayer + 1;
				lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel;

			//!explode coarse states on interface edge voxels along y direction edges
			for (int j = refLayer; j <= jmax - refLayer; j++)
			{
				//!inlet bottom edge
				int i = refLayerX - 1;
				int k = refLayer - 1;
				lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

				//!inlet top edge
				i = refLayerX - 1;
				k = kmax - refLayer + 1;
				lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

				//!outlet bottom edge
				i = imax - refLayer + 1;
				k = refLayer - 1;
				lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

				//!outlet top edge
				i = imax - refLayer + 1;
				k = kmax - refLayer + 1;
				lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);
			} //!end j
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel;

			//!explode coarse states on interface edge voxels along z direction edges
			for (int k = refLayer; k <= kmax - refLayer; k++)
			{
				//!inlet left edge
				int i = refLayerX - 1;
				int j = refLayer - 1;
				lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

				//!inlet right edge
				i = refLayerX - 1;
				j = jmax - refLayer + 1;
				lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

				//!outlet left edge
				i = imax - refLayer + 1;
				j = refLayer - 1;
				lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

				//!outlet right edge
				i = imax - refLayer + 1;
				j = jmax - refLayer + 1;
				lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);
			} //!end k
		} //!end omp
	} //!end else if

	//!transport in all coarse voxels
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX - 1; i <= imax - refLayer + 1; i++) //!all voxels in x direction
		{
			for (int j = refLayer - 1; j <= jmax - refLayer + 1; j++) //!all voxels in y direction
			{
				for (int k = refLayer - 1; k <= kmax - refLayer + 1; k++)
				{
					lattice_[level][i][j][k]->transport(negAlternatingL0);
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

	//!transport in all fine voxels
	level = 1;
	refLayer = bc_.getRefinementLayer() * pow(2, level);
	refLayerX = bc_.getRefinementLayerX() * pow(2, level);
	imax = bc_.getNumberOfVoxelsX() * pow(2., level) - 1; //!index of last voxel in x direction in ->this level
	jmax = bc_.getNumberOfVoxelsY() * pow(2., level) - 1; //!index of last voxel in y direction in ->this level
	kmax = bc_.getNumberOfVoxelsZ() * pow(2., level) - 1; //!index of last voxel in z direction in ->this level
	int jmid = jmax / 2;
	int kmid = kmax / 2;

	//!transport in all fine voxels
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all voxels in y direction except walls
			{
				//!bottom region
				for (int k = 0; k <= refLayer - 1; k++) //!all fine voxels in refined bottom region except walls
				{
					lattice_[level][i][j][k]->transport(negAlternatingL1);
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all voxels in y direction except walls
			{
				//!top region
				for (int k = kmax - refLayer + 1; k <= kmax; k++) //!all fine voxels in refined top region except walls
				{
					lattice_[level][i][j][k]->transport(negAlternatingL1);
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int k = refLayer; k <= kmax - refLayer; k++) //!all voxels in z direction except walls and voxels considered before
			{
				//!left region
				for (int j = 0; j <= refLayer - 1; j++) //!all fine voxels in refined left region
				{
					lattice_[level][i][j][k]->transport(negAlternatingL1);
				} //!end k
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int k = refLayer; k <= kmax - refLayer; k++) //!all voxels in z direction except walls and voxels considered before
			{
				//!right region
				for (int j = jmax - refLayer + 1; j <= jmax; j++) //!all fine voxels in refined right region
				{
					lattice_[level][i][j][k]->transport(negAlternatingL1);
				} //!end j
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer; j <= jmax - refLayer; j++) //!all voxels in y direction except except walls and voxels considered before
		{
			for (int k = refLayer; k <= kmax - refLayer; k++) //!all voxels in z direction except walls and voxels considered before
			{
				//!inlet region
				for (int i = 0; i <= refLayerX - 1; i++) //!all fine voxels in refined inlet region except boundary
				{
					lattice_[level][i][j][k]->transport(negAlternatingL1);
				} //!end i
			} //!end k
		} //!end j
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer; j <= jmax - refLayer; j++) //!all voxels in y direction except except walls and voxels considered before
		{
			for (int k = refLayer; k <= kmax - refLayer; k++) //!all voxels in z direction except walls and voxels considered before
			{
				//!outlet region
				for (int i = imax - refLayer + 1; i <= imax; i++) //!all fine voxels in refined outlet region except boundary
				{
					lattice_[level][i][j][k]->transport(negAlternatingL1);
				} //!end i
			} //!end k
		} //!end j
	} //!end omp

	//!collide in inner fine voxels (no collision in fine interface voxels!)
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all voxels in y direction
			{
				//!bottom refined region
				for (int k = 0; k <= refLayer - 3; k++) //!all fine inner voxels in refined bottom region
				{
					lattice_[level][i][j][k]->collideMRT_RMbasis(bc_, negAlternatingL1); //!perform collide in voxel with recursive-regularized collision operator (RR)
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all voxels in y direction
			{
				//!top refined region
				for (int k = kmax - refLayer + 3; k <= kmax; k++) //!all fine inner voxels in refined top region
				{
					lattice_[level][i][j][k]->collideMRT_RMbasis(bc_, negAlternatingL1); //!perform collide in voxel with recursive-regularized collision operator (RR)
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction
			{
				//!left refined region
				for (int j = 0; j <= refLayer - 3; j++) //!all fine inner voxels in refined left region
				{
					lattice_[level][i][j][k]->collideMRT_RMbasis(bc_, negAlternatingL1); //!perform collide in voxel with recursive-regularized collision operator (RR)
				} //!end j
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction
			{
				//!right refined region
				for (int j = jmax - refLayer + 3; j <= jmax; j++) //!all fine inner voxels in refined right region
				{
					lattice_[level][i][j][k]->collideMRT_RMbasis(bc_, negAlternatingL1); //!perform collide in voxel with recursive-regularized collision operator (RR)
				} //!end j
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer - 2; j <= jmax - refLayer + 2; j++) //!all voxels in y direction
		{
			for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction
			{
				//!inlet refined region
				for (int i = 0; i <= refLayerX - 3; i++) //!all fine inner voxels in refined bottom region
				{
					lattice_[level][i][j][k]->collideMRT_RMbasis(bc_, negAlternatingL1); //!perform collide in voxel with recursive-regularized collision operator (RR)
				} //!end i
			} //!end k
		} //!end j
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer - 2; j <= jmax - refLayer + 2; j++) //!all voxels in y direction
		{
			for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction
			{
				//!outlet refined region
				for (int i = imax - refLayer + 3; i <= imax; i++) //!all fine inner voxels in refined top region
				{
					lattice_[level][i][j][k]->collideMRT_RMbasis(bc_, negAlternatingL1); //!perform collide in voxel with recursive-regularized collision operator (RR)
				} //!end i
			} //!end k
		} //!end j
	} //!end omp

//!transport in all fine voxels
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all voxels in y direction except walls
			{
				//!bottom region
				for (int k = 0; k <= refLayer - 1; k++) //!all fine voxels in refined bottom region except walls
				{
					lattice_[level][i][j][k]->transport(alternatingL1);
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all voxels in y direction except walls
			{
				//!top region
				for (int k = kmax - refLayer + 1; k <= kmax; k++) //!all fine voxels in refined top region except walls
				{
					lattice_[level][i][j][k]->transport(alternatingL1);
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int k = refLayer; k <= kmax - refLayer; k++) //!all voxels in z direction except walls and voxels considered before
			{
				//!left region
				for (int j = 0; j <= refLayer - 1; j++) //!all fine voxels in refined left region
				{
					lattice_[level][i][j][k]->transport(alternatingL1);
				} //!end k
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int k = refLayer; k <= kmax - refLayer; k++) //!all voxels in z direction except walls and voxels considered before
			{
				//!right region
				for (int j = jmax - refLayer + 1; j <= jmax; j++) //!all fine voxels in refined right region
				{
					lattice_[level][i][j][k]->transport(alternatingL1);
				} //!end j
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer; j <= jmax - refLayer; j++) //!all voxels in y direction except except walls and voxels considered before
		{
			for (int k = refLayer; k <= kmax - refLayer; k++) //!all voxels in z direction except walls and voxels considered before
			{
				//!inlet region
				for (int i = 0; i <= refLayerX - 1; i++) //!all fine voxels in refined inlet region except boundary
				{
					lattice_[level][i][j][k]->transport(alternatingL1);
				} //!end i
			} //!end k
		} //!end j
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer; j <= jmax - refLayer; j++) //!all voxels in y direction except except walls and voxels considered before
		{
			for (int k = refLayer; k <= kmax - refLayer; k++) //!all voxels in z direction except walls and voxels considered before
			{
				//!outlet region
				for (int i = imax - refLayer + 1; i <= imax; i++) //!all fine voxels in refined outlet region except boundary
				{
					lattice_[level][i][j][k]->transport(alternatingL1);
				} //!end i
			} //!end k
		} //!end j
	} //!end omp

	//!level 0
	level = 0;
	refLayer = bc_.getRefinementLayer();
	refLayerX = bc_.getRefinementLayerX() * pow(2, level);
	imax = bc_.getNumberOfVoxelsX() - 1; //!index of last voxel in x direction in ->this level
	jmax = bc_.getNumberOfVoxelsY() - 1; //!index of last voxel in y direction in ->this level
	kmax = bc_.getNumberOfVoxelsZ() - 1; //!index of last voxel in z direction in ->this level

	//!coalesce states in coarse interface voxels
	//!ATTENTION: COARSE EDGE AND CORNER VOXELS NEED SEPARATE TREATMENT
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do outer loop in parallel
		for (int i = refLayerX; i <= imax - refLayer; i++) //!all voxels in x direction
		{
			for (int j = refLayer; j <= jmax - refLayer; j++) //!all voxels in y direction
			{
				//!coalesce states in coarse bottom interface voxels
				int k = refLayer - 1;
				lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceBottom);

				//!coalesce states in coarse top interface voxels
				k = kmax - refLayer + 1;
				lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceTop);
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do outer loop in parallel
		for (int i = refLayerX; i <= imax - refLayer; i++) //!all coarse interface voxels in x direction
		{
			for (int k = refLayer; k <= kmax - refLayer; k++) //!all coarse interface voxels in z direction
			{
				//!coalesce states in coarse right interface voxels
				int j = refLayer - 1;
				lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceLeft);

				//!coalesce states in coarse left interface voxels
				j = jmax - refLayer + 1;
				lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceRight);
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do outer loop in parallel
		for (int j = refLayer; j <= jmax - refLayer; j++) //!all coarse interface voxels in y direction
		{
			for (int k = refLayer; k <= kmax - refLayer; k++) //!all coarse interface voxels in z direction
			{
				//!coalesce states in coarse inlet interface voxels
				int i = refLayerX - 1;
				lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceInlet);

				//!coalesce states in coarse outlet interface voxels
				i = imax - refLayer + 1;
				lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceOutlet);
			} //!end k
		} //!end i
	} //!end omp

	/*COARSE INTERFACE EDGES (NO CORNERS!)*/
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel;
		//!coalesce fine states in coarse interface edge voxels along x direction edges
		for (int i = refLayerX; i <= imax - refLayer; i++)
		{
			//!bottom left edge
			int j = refLayer - 1;
			int k = refLayer - 1;
			lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceBottomLeft);

			//!bottom right edge
			j = jmax - refLayer + 1;
			k = refLayer - 1;
			lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceBottomRight);

			//!top left edge
			j = refLayer - 1;
			k = kmax - refLayer + 1;
			lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceTopLeft);

			//!top right edge
			j = jmax - refLayer + 1;
			k = kmax - refLayer + 1;
			lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceTopRight);
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel;
		//!coalesce fine states in coarse interface edge voxels along y direction edges
		for (int j = refLayer; j <= jmax - refLayer; j++)
		{
			//!inlet bottom edge
			int i = refLayerX - 1;
			int k = refLayer - 1;
			lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceInletBottom);

			//!inlet top edge
			i = refLayerX - 1;
			k = kmax - refLayer + 1;
			lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceInletTop);

			//!outlet bottom edge
			i = imax - refLayer + 1;
			k = refLayer - 1;
			lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceOutletBottom);

			//!outlet top edge
			i = imax - refLayer + 1;
			k = kmax - refLayer + 1;
			lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceOutletTop);
		} //!end j
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel;
		//!coalesce fine states in coarse interface edge voxels along z direction edges
		for (int k = refLayer; k <= kmax - refLayer; k++)
		{
			//!inlet left edge
			int i = refLayerX - 1;
			int j = refLayer - 1;
			lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceInletLeft);

			//!inlet right edge
			i = refLayerX - 1;
			j = jmax - refLayer + 1;
			lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceInletRight);

			//!outlet left edge
			i = imax - refLayer + 1;
			j = refLayer - 1;
			lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceOutletLeft);

			//!outlet right edge
			i = imax - refLayer + 1;
			j = jmax - refLayer + 1;
			lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceOutletRight);
		} //!end k
	} //!end omp

	/*COARSE INTERFACE CORNERS*/
	//!inlet bottom left corner
	int i = refLayerX - 1;
	int j = refLayer - 1;
	int k = refLayer - 1;
	lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceInletBottomLeft);

	//!inlet bottom right corner
	i = refLayerX - 1;
	j = jmax - refLayer + 1;
	k = refLayer - 1;
	lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceInletBottomRight);

	//!inlet top left corner
	i = refLayerX - 1;
	j = refLayer - 1;
	k = kmax - refLayer + 1;
	lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceInletTopLeft);

	//!inlet top right corner
	i = refLayerX - 1;
	j = jmax - refLayer + 1;
	k = kmax - refLayer + 1;
	lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceInletTopRight);

	//!outlet bottom left corner
	i = imax - refLayer + 1;
	j = refLayer - 1;
	k = refLayer - 1;
	lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceOutletBottomLeft);

	//!outlet bottom right corner
	i = imax - refLayer + 1;
	j = jmax - refLayer + 1;
	k = refLayer - 1;
	lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceOutletBottomRight);

	//!outlet top left corner
	i = imax - refLayer + 1;
	j = refLayer - 1;
	k = kmax - refLayer + 1;
	lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceOutletTopLeft);

	//!outlet top right corner
	i = imax - refLayer + 1;
	j = jmax - refLayer + 1;
	k = kmax - refLayer + 1;
	lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceOutletTopRight);
} //!end timeStepNestedMRT_RMbasis


void Lattice::nestedTimeSteppingRR(bool& alternatingL1, bool& alternatingL0, Control& bc_) //!perform nested time step with sequence of collide stream step for all level
{
	if (bc_.getCubicMachCorrection() == "yes")
	{
		int refLayer, refLayerX; //!number of refined cells
		bool negAlternatingL1 = !alternatingL1;
		bool negAlternatingL0 = !alternatingL0;

		//!STRING VARIABLES FOR DIFFERENCIATING BETWEEN VARIOUS INTERFACE VOXELS -> NEEDED FOR EXPLODE AND COALESCE METHOD
		//!'core' interface voxels
		string interfaceBottom = "interfaceBottom";
		string interfaceTop = "interfaceTop";
		string interfaceLeft = "interfaceLeft";
		string interfaceRight = "interfaceRight";
		string interfaceInlet = "interfaceInlet";
		string interfaceOutlet = "interfaceOutlet";

		//!'edge' interface voxels
		string interfaceBottomLeft = "interfaceBottomLeft";
		string interfaceBottomRight = "interfaceBottomRight";
		string interfaceTopLeft = "interfaceTopLeft";
		string interfaceTopRight = "interfaceTopRight";
		string interfaceInletBottom = "interfaceInletBottom";
		string interfaceInletTop = "interfaceInletTop";
		string interfaceOutletBottom = "interfaceOutletBottom";
		string interfaceOutletTop = "interfaceOutletTop";
		string interfaceInletLeft = "interfaceInletLeft";
		string interfaceInletRight = "interfaceInletRight";
		string interfaceOutletLeft = "interfaceOutletLeft";
		string interfaceOutletRight = "interfaceOutletRight";

		//!'corner' interface voxels
		string interfaceInletBottomLeft = "interfaceInletBottomLeft";
		string interfaceInletBottomRight = "interfaceInletBottomRight";
		string interfaceInletTopLeft = "interfaceInletTopLeft";
		string interfaceInletTopRight = "interfaceInletTopRight";
		string interfaceOutletBottomLeft = "interfaceOutletBottomLeft";
		string interfaceOutletBottomRight = "interfaceOutletBottomRight";
		string interfaceOutletTopLeft = "interfaceOutletTopLeft";
		string interfaceOutletTopRight = "interfaceOutletTopRight";

		//--------------------------------------------------//
		//--------------------------------------------------//
		//             collision and transport				//
		//--------------------------------------------------//
		//--------------------------------------------------//

		//!current distribution (t=0) level 0: alternating
		//!current distribution (t=0) level 1: alternating

		int level = 1;
		refLayer = bc_.getRefinementLayer() * pow(2, level);
		refLayerX = bc_.getRefinementLayerX() * pow(2, level);
		int imax = bc_.getNumberOfVoxelsX() * pow(2., level) - 1; //!index of last voxel in x direction in ->this level
		int jmax = bc_.getNumberOfVoxelsY() * pow(2., level) - 1; //!index of last voxel in y direction in ->this level
		int kmax = bc_.getNumberOfVoxelsZ() * pow(2., level) - 1; //!index of last voxel in z direction in ->this level

		//!safe pre-collide distributions for HRR collision operator on fine inner voxels
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!all voxels in x direction
			{
				for (int j = 0; j <= jmax; j++) //!all voxels in y direction
				{
					//!bottom refined region
					for (int k = 0; k <= refLayer - 3; k++) //!all fine inner voxels in refined bottom region
					{
						for (int dir = 0; dir < 19; dir++)
						{
							double distrTemp = lattice_[level][i][j][k]->getDistribution(alternatingL1, dir); //!initialize temporary distribution
							lattice_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
							double psiTemp = lattice_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
							lattice_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
						}
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!all voxels in x direction
			{
				for (int j = 0; j <= jmax; j++) //!all voxels in y direction
				{
					//!top refined region
					for (int k = kmax - refLayer + 3; k <= kmax; k++) //!all fine inner voxels in refined top region
					{
						for (int dir = 0; dir < 19; dir++)
						{
							double distrTemp = lattice_[level][i][j][k]->getDistribution(alternatingL1, dir); //!initialize temporary distribution
							lattice_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
							double psiTemp = lattice_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
							lattice_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
						}
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!all voxels in x direction
			{
				for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction
				{
					//!left refined region
					for (int j = 0; j <= refLayer - 3; j++) //!all fine inner voxels in refined left region
					{
						for (int dir = 0; dir < 19; dir++)
						{
							double distrTemp = lattice_[level][i][j][k]->getDistribution(alternatingL1, dir); //!initialize temporary distribution
							lattice_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
							double psiTemp = lattice_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
							lattice_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
						}
					} //!end j
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!all voxels in x direction
			{
				for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction
				{
					//!right refined region
					for (int j = jmax - refLayer + 3; j <= jmax; j++) //!all fine inner voxels in refined right region
					{
						for (int dir = 0; dir < 19; dir++)
						{
							double distrTemp = lattice_[level][i][j][k]->getDistribution(alternatingL1, dir); //!initialize temporary distribution
							lattice_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
							double psiTemp = lattice_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
							lattice_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
						}
					} //!end j
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer - 2; j <= jmax - refLayer + 2; j++) //!all voxels in y direction
			{
				for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction
				{
					//!inlet refined region
					for (int i = 0; i <= refLayerX - 3; i++) //!all fine inner voxels in refined bottom region
					{
						for (int dir = 0; dir < 19; dir++)
						{
							double distrTemp = lattice_[level][i][j][k]->getDistribution(alternatingL1, dir); //!initialize temporary distribution
							lattice_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
							double psiTemp = lattice_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
							lattice_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
						}
					} //!end i
				} //!end k
			} //!end j
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer - 2; j <= jmax - refLayer + 2; j++) //!all voxels in y direction
			{
				for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction
				{
					//!outlet refined region
					for (int i = imax - refLayer + 3; i <= imax; i++) //!all fine inner voxels in refined top region
					{
						for (int dir = 0; dir < 19; dir++)
						{
							double distrTemp = lattice_[level][i][j][k]->getDistribution(alternatingL1, dir); //!initialize temporary distribution
							lattice_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
							double psiTemp = lattice_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
							lattice_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
						}
					} //!end i
				} //!end k
			} //!end j
		} //!end omp

		level = 0;
		imax = bc_.getNumberOfVoxelsX() - 1; //!index of last voxel in x direction in ->this level
		jmax = bc_.getNumberOfVoxelsY() - 1; //!index of last voxel in y direction in ->this level
		kmax = bc_.getNumberOfVoxelsZ() - 1; //!index of last voxel in z direction in ->this level
		refLayer = bc_.getRefinementLayer(); //!number of refined coarse voxels
		refLayerX = bc_.getRefinementLayerX() * pow(2, level);

		//!safe pre-collide distributions for HRR collision operator on coarse voxels
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int i = refLayerX - 1; i <= imax - refLayer + 1; i++) //!all nodes in x direction
			{
				for (int j = refLayer - 1; j <= jmax - refLayer + 1; j++) //!all nodes in y direction
				{
					for (int k = refLayer - 1; k <= kmax - refLayer + 1; k++) //!nodes in z direction from coarse lower interface to upper interface
					{
						for (int dir = 0; dir < 19; dir++)
						{
							double distrTemp = lattice_[level][i][j][k]->getDistribution(alternatingL0, dir); //!initialize temporary distribution
							lattice_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
							double psiTemp = lattice_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
							lattice_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
						}
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

		//!collide in all coarse voxels
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int i = refLayerX - 1; i <= imax - refLayer + 1; i++) //!all nodes in x direction
			{
				for (int j = refLayer - 1; j <= jmax - refLayer + 1; j++) //!all nodes in y direction
				{
					for (int k = refLayer - 1; k <= kmax - refLayer + 1; k++) //!nodes in z direction from coarse lower interface to upper interface
					{
						lattice_[level][i][j][k]->collideRR(bc_, alternatingL0); //!perform collide in voxel with hybrid recursive-regularized collision operator (HRR)
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

		//!explode states in coarse interface voxels
		//!INTERFACE EDGES AND CORNERS INCLUDED (ALL DISTRIBUTIONS ARE BEING EXPLODED!)
		if (bc_.getExplosionOrder() == 1)
		{
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
			{
#pragma omp for //!do outer loop in parallel
				for (int i = refLayerX - 1; i <= imax - refLayer + 1; i++) //!all coarse interface voxels in x direction
				{
					for (int j = refLayer - 1; j <= jmax - refLayer + 1; j++) //!all coarse interface voxels in y direction
					{
						//!explode states in coarse bottom interface voxels
						int k = refLayer - 1;
						lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

						//!explode states in coarse top interface voxels
						k = kmax - refLayer + 1;
						lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);
					} //!end j
				} //!end i
			} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
			{
#pragma omp for //!do outer loop in parallel
				for (int i = refLayerX - 1; i <= imax - refLayer + 1; i++) //!all coarse interface voxels in x direction
				{
					for (int k = refLayer; k <= kmax - refLayer; k++) //!all coarse interface voxels in z direction
					{
						//!explode states in coarse right interface voxels
						int j = refLayer - 1;
						lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

						//!explode states in coarse left interface voxels
						j = jmax - refLayer + 1;
						lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);
					} //!end k
				} //!end i
			} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
			{
#pragma omp for //!do outer loop in parallel
				for (int j = refLayer; j <= jmax - refLayer; j++) //!all coarse interface voxels in y direction
				{
					for (int k = refLayer; k <= kmax - refLayer; k++) //!all coarse interface voxels in z direction
					{
						//!explode states in coarse inlet interface voxels
						int i = refLayerX - 1;
						lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

						//!explode states in coarse outlet interface voxels
						i = imax - refLayer + 1;
						lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);
					} //!end k
				} //!end i
			} //!end omp
		} //!end if

		else if (bc_.getExplosionOrder() == 2)
		{
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
			{
#pragma omp for //!do outer loop in parallel
				for (int i = refLayerX; i <= imax - refLayer; i++) //!all coarse interface voxels in x direction
				{
					for (int j = refLayer; j <= jmax - refLayer; j++) //!all coarse interface voxels in y direction
					{
						//!explode states in coarse bottom interface voxels
						int k = refLayer - 1;
						lattice_[level][i][j][k]->explodeLinear(alternatingL1, alternatingL0, bc_, interfaceBottom);

						//!explode states in coarse top interface voxels
						k = kmax - refLayer + 1;
						lattice_[level][i][j][k]->explodeLinear(alternatingL1, alternatingL0, bc_, interfaceTop);
					} //!end j
				} //!end i
			} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
			{
#pragma omp for //!do outer loop in parallel
				for (int i = refLayerX; i <= imax - refLayer; i++) //!all coarse interface voxels in x direction
				{
					for (int k = refLayer; k <= kmax - refLayer; k++) //!all coarse interface voxels in z direction
					{
						//!explode states in coarse right interface voxels
						int j = refLayer - 1;
						lattice_[level][i][j][k]->explodeLinear(alternatingL1, alternatingL0, bc_, interfaceLeft);

						//!explode states in coarse left interface voxels
						j = jmax - refLayer + 1;
						lattice_[level][i][j][k]->explodeLinear(alternatingL1, alternatingL0, bc_, interfaceRight);
					} //!end k
				} //!end i
			} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
			{
#pragma omp for //!do outer loop in parallel
				for (int j = refLayer; j <= jmax - refLayer; j++) //!all coarse interface voxels in y direction
				{
					for (int k = refLayer; k <= kmax - refLayer; k++) //!all coarse interface voxels in z direction
					{
						//!explode states in coarse inlet interface voxels
						int i = refLayerX - 1;
						lattice_[level][i][j][k]->explodeLinear(alternatingL1, alternatingL0, bc_, interfaceInlet);

						//!explode states in coarse outlet interface voxels
						i = imax - refLayer + 1;
						lattice_[level][i][j][k]->explodeLinear(alternatingL1, alternatingL0, bc_, interfaceOutlet);
					} //!end k
				} //!end i
			} //!end omp

			/*COARSE INTERFACE EDGES*/
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
			{
#pragma omp for //!do most outer loop in parallel;

				//!explode coarse states on interface edge voxels along x direction edges
				for (int i = refLayerX - 1; i <= imax - refLayer + 1; i++)
				{
					//!bottom left edge
					int j = refLayer - 1;
					int k = refLayer - 1;
					lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

					//!bottom right edge
					j = jmax - refLayer + 1;
					k = refLayer - 1;
					lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

					//!top left edge
					j = refLayer - 1;
					k = kmax - refLayer + 1;
					lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

					//!top right edge
					j = jmax - refLayer + 1;
					k = kmax - refLayer + 1;
					lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);
				} //!end i
			} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
			{
#pragma omp for //!do most outer loop in parallel;

				//!explode coarse states on interface edge voxels along y direction edges
				for (int j = refLayer; j <= jmax - refLayer; j++)
				{
					//!inlet bottom edge
					int i = refLayerX - 1;
					int k = refLayer - 1;
					lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

					//!inlet top edge
					i = refLayerX - 1;
					k = kmax - refLayer + 1;
					lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

					//!outlet bottom edge
					i = imax - refLayer + 1;
					k = refLayer - 1;
					lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

					//!outlet top edge
					i = imax - refLayer + 1;
					k = kmax - refLayer + 1;
					lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);
				} //!end j
			} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
			{
#pragma omp for //!do most outer loop in parallel;

				//!explode coarse states on interface edge voxels along z direction edges
				for (int k = refLayer; k <= kmax - refLayer; k++)
				{
					//!inlet left edge
					int i = refLayerX - 1;
					int j = refLayer - 1;
					lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

					//!inlet right edge
					i = refLayerX - 1;
					j = jmax - refLayer + 1;
					lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

					//!outlet left edge
					i = imax - refLayer + 1;
					j = refLayer - 1;
					lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

					//!outlet right edge
					i = imax - refLayer + 1;
					j = jmax - refLayer + 1;
					lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);
				} //!end k
			} //!end omp
		} //!end else if

		//!transport in all coarse voxels
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX - 1; i <= imax - refLayer + 1; i++) //!all voxels in x direction
			{
				for (int j = refLayer - 1; j <= jmax - refLayer + 1; j++) //!all voxels in y direction
				{
					for (int k = refLayer - 1; k <= kmax - refLayer + 1; k++)
					{
						lattice_[level][i][j][k]->transport(negAlternatingL0);
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

		level = 1;
		refLayer = bc_.getRefinementLayer() * pow(2, level);
		refLayerX = bc_.getRefinementLayerX() * pow(2, level);
		imax = bc_.getNumberOfVoxelsX() * pow(2., level) - 1; //!index of last voxel in x direction in ->this level
		jmax = bc_.getNumberOfVoxelsY() * pow(2., level) - 1; //!index of last voxel in y direction in ->this level
		kmax = bc_.getNumberOfVoxelsZ() * pow(2., level) - 1; //!index of last voxel in z direction in ->this level

		//!safe pre-collide distributions for HRR collision operator on fine first-layer voxels after explode
		//!for inner fine voxels states have already been saved...
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX - 2; i <= imax - refLayer + 2; i++) //!all voxels in x direction
			{
				for (int j = refLayer - 2; j <= jmax - refLayer + 2; j++) //!all voxels in y direction
				{
					//!bottom refined region
					int k = refLayer - 2;
					for (int dir = 0; dir < 19; dir++)
					{
						double distrTemp = lattice_[level][i][j][k]->getDistribution(alternatingL1, dir); //!initialize temporary distribution
						lattice_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
						double psiTemp = lattice_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
						lattice_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
					}
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX - 2; i <= imax - refLayer + 2; i++) //!all voxels in x direction
			{
				for (int j = refLayer - 2; j <= jmax - refLayer + 2; j++) //!all voxels in y direction
				{
					//!top refined region
					int k = kmax - refLayer + 2;
					for (int dir = 0; dir < 19; dir++)
					{
						double distrTemp = lattice_[level][i][j][k]->getDistribution(alternatingL1, dir); //!initialize temporary distribution
						lattice_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
						double psiTemp = lattice_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
						lattice_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
					}
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX - 2; i <= imax - refLayer + 2; i++) //!all voxels in x direction
			{
				for (int k = refLayer - 1; k <= kmax - refLayer + 1; k++) //!all voxels in z direction
				{
					//!left refined region
					int j = refLayer - 2;
					for (int dir = 0; dir < 19; dir++)
					{
						double distrTemp = lattice_[level][i][j][k]->getDistribution(alternatingL1, dir); //!initialize temporary distribution
						lattice_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
						double psiTemp = lattice_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
						lattice_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
					}
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX - 2; i <= imax - refLayer + 2; i++) //!all voxels in x direction
			{
				for (int k = refLayer - 1; k <= kmax - refLayer + 1; k++) //!all voxels in z direction
				{
					//!right refined region
					int j = jmax - refLayer + 2;
					for (int dir = 0; dir < 19; dir++)
					{
						double distrTemp = lattice_[level][i][j][k]->getDistribution(alternatingL1, dir); //!initialize temporary distribution
						lattice_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
						double psiTemp = lattice_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
						lattice_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
					}
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer - 1; j <= jmax - refLayer + 1; j++) //!all voxels in y direction
			{
				for (int k = refLayer - 1; k <= kmax - refLayer + 1; k++) //!all voxels in z direction
				{
					//!inlet refined region
					int i = refLayerX - 2;
					for (int dir = 0; dir < 19; dir++)
					{
						double distrTemp = lattice_[level][i][j][k]->getDistribution(alternatingL1, dir); //!initialize temporary distribution
						lattice_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
						double psiTemp = lattice_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
						lattice_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
					}
				} //!end k
			} //!end j
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer - 1; j <= jmax - refLayer + 1; j++) //!all voxels in y direction
			{
				for (int k = refLayer - 1; k <= kmax - refLayer + 1; k++) //!all voxels in z direction
				{
					//!inlet refined region
					int i = imax - refLayer + 2;
					for (int dir = 0; dir < 19; dir++)
					{
						double distrTemp = lattice_[level][i][j][k]->getDistribution(alternatingL1, dir); //!initialize temporary distribution
						lattice_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
						double psiTemp = lattice_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
						lattice_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
					}
				} //!end k
			} //!end j
		} //!end omp

		//!collide in inner fine voxels (no collision in fine interface voxels!)
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!all voxels in x direction
			{
				for (int j = 0; j <= jmax; j++) //!all voxels in y direction
				{
					//!bottom refined region
					for (int k = 0; k <= refLayer - 3; k++) //!all fine inner voxels in refined bottom region
					{
						lattice_[level][i][j][k]->collideRR(bc_, alternatingL1); //!perform collide in voxel with hybrid recursive-regularized collision operator (HRR)
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!all voxels in x direction
			{
				for (int j = 0; j <= jmax; j++) //!all voxels in y direction
				{
					//!top refined region
					for (int k = kmax - refLayer + 3; k <= kmax; k++) //!all fine inner voxels in refined top region
					{
						lattice_[level][i][j][k]->collideRR(bc_, alternatingL1); //!perform collide in voxel with hybrid recursive-regularized collision operator (HRR)
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!all voxels in x direction
			{
				for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction
				{
					//!left refined region
					for (int j = 0; j <= refLayer - 3; j++) //!all fine inner voxels in refined left region
					{
						lattice_[level][i][j][k]->collideRR(bc_, alternatingL1); //!perform collide in voxel with hybrid recursive-regularized collision operator (HRR)
					} //!end j
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!all voxels in x direction
			{
				for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction
				{
					//!right refined region
					for (int j = jmax - refLayer + 3; j <= jmax; j++) //!all fine inner voxels in refined right region
					{
						lattice_[level][i][j][k]->collideRR(bc_, alternatingL1); //!perform collide in voxel with hybrid recursive-regularized collision operator (HRR)
					} //!end j
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer - 2; j <= jmax - refLayer + 2; j++) //!all voxels in y direction
			{
				for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction
				{
					//!inlet refined region
					for (int i = 0; i <= refLayerX - 3; i++) //!all fine inner voxels in refined bottom region
					{
						lattice_[level][i][j][k]->collideRR(bc_, alternatingL1); //!perform collide in voxel with hybrid recursive-regularized collision operator (HRR)
					} //!end i
				} //!end k
			} //!end j
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer - 2; j <= jmax - refLayer + 2; j++) //!all voxels in y direction
			{
				for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction
				{
					//!outlet refined region
					for (int i = imax - refLayer + 3; i <= imax; i++) //!all fine inner voxels in refined top region
					{
						lattice_[level][i][j][k]->collideRR(bc_, alternatingL1); //!perform collide in voxel with hybrid recursive-regularized collision operator (HRR)
					} //!end i
				} //!end k
			} //!end j
		} //!end omp

		//!transport in all fine voxels
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!all voxels in x direction
			{
				for (int j = 0; j <= jmax; j++) //!all voxels in y direction except walls
				{
					//!bottom region
					for (int k = 0; k <= refLayer - 1; k++) //!all fine voxels in refined bottom region except walls
					{
						lattice_[level][i][j][k]->transport(negAlternatingL1);
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!all voxels in x direction
			{
				for (int j = 0; j <= jmax; j++) //!all voxels in y direction except walls
				{
					//!top region
					for (int k = kmax - refLayer + 1; k <= kmax; k++) //!all fine voxels in refined top region except walls
					{
						lattice_[level][i][j][k]->transport(negAlternatingL1);
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!all voxels in x direction
			{
				for (int k = refLayer; k <= kmax - refLayer; k++) //!all voxels in z direction except walls and voxels considered before
				{
					//!left region
					for (int j = 0; j <= refLayer - 1; j++) //!all fine voxels in refined left region
					{
						lattice_[level][i][j][k]->transport(negAlternatingL1);
					} //!end k
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!all voxels in x direction
			{
				for (int k = refLayer; k <= kmax - refLayer; k++) //!all voxels in z direction except walls and voxels considered before
				{
					//!right region
					for (int j = jmax - refLayer + 1; j <= jmax; j++) //!all fine voxels in refined right region
					{
						lattice_[level][i][j][k]->transport(negAlternatingL1);
					} //!end j
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer; j <= jmax - refLayer; j++) //!all voxels in y direction except except walls and voxels considered before
			{
				for (int k = refLayer; k <= kmax - refLayer; k++) //!all voxels in z direction except walls and voxels considered before
				{
					//!inlet region
					for (int i = 0; i <= refLayerX - 1; i++) //!all fine voxels in refined inlet region except boundary
					{
						lattice_[level][i][j][k]->transport(negAlternatingL1);
					} //!end i
				} //!end k
			} //!end j
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer; j <= jmax - refLayer; j++) //!all voxels in y direction except except walls and voxels considered before
			{
				for (int k = refLayer; k <= kmax - refLayer; k++) //!all voxels in z direction except walls and voxels considered before
				{
					//!outlet region
					for (int i = imax - refLayer + 1; i <= imax; i++) //!all fine voxels in refined outlet region except boundary
					{
						lattice_[level][i][j][k]->transport(negAlternatingL1);
					} //!end i
				} //!end k
			} //!end j
		} //!end omp

		//!safe pre-collide distributions for HRR collision operator on fine inner and first-layer voxels
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!all voxels in x direction
			{
				for (int j = 0; j <= jmax; j++) //!all voxels in y direction
				{
					//!bottom refined region
					for (int k = 0; k <= refLayer - 2; k++) //!all fine inner voxels in refined bottom region
					{
						for (int dir = 0; dir < 19; dir++)
						{
							double distrTemp = lattice_[level][i][j][k]->getDistribution(negAlternatingL1, dir); //!initialize temporary distribution
							lattice_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
							double psiTemp = lattice_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
							lattice_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
						}
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!all voxels in x direction
			{
				for (int j = 0; j <= jmax; j++) //!all voxels in y direction
				{
					//!top refined region
					for (int k = kmax - refLayer + 2; k <= kmax; k++) //!all fine inner voxels in refined top region
					{
						for (int dir = 0; dir < 19; dir++)
						{
							double distrTemp = lattice_[level][i][j][k]->getDistribution(negAlternatingL1, dir); //!initialize temporary distribution
							lattice_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
							double psiTemp = lattice_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
							lattice_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
						}
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!all voxels in x direction
			{
				for (int k = refLayer - 1; k <= kmax - refLayer + 1; k++) //!all voxels in z direction
				{
					//!left refined region
					for (int j = 0; j <= refLayer - 2; j++) //!all fine inner voxels in refined left region
					{
						for (int dir = 0; dir < 19; dir++)
						{
							double distrTemp = lattice_[level][i][j][k]->getDistribution(negAlternatingL1, dir); //!initialize temporary distribution
							lattice_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
							double psiTemp = lattice_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
							lattice_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
						}
					} //!end j
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!all voxels in x direction
			{
				for (int k = refLayer - 1; k <= kmax - refLayer + 1; k++) //!all voxels in z direction
				{
					//!right refined region
					for (int j = jmax - refLayer + 2; j <= jmax; j++) //!all fine inner voxels in refined right region
					{
						for (int dir = 0; dir < 19; dir++)
						{
							double distrTemp = lattice_[level][i][j][k]->getDistribution(negAlternatingL1, dir); //!initialize temporary distribution
							lattice_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
							double psiTemp = lattice_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
							lattice_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
						}
					} //!end j
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer - 1; j <= jmax - refLayer + 1; j++) //!all voxels in y direction
			{
				for (int k = refLayer - 1; k <= kmax - refLayer + 1; k++) //!all voxels in z direction
				{
					//!inlet refined region
					for (int i = 0; i <= refLayerX - 2; i++) //!all fine inner voxels in refined bottom region
					{
						for (int dir = 0; dir < 19; dir++)
						{
							double distrTemp = lattice_[level][i][j][k]->getDistribution(negAlternatingL1, dir); //!initialize temporary distribution
							lattice_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
							double psiTemp = lattice_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
							lattice_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
						}
					} //!end i
				} //!end k
			} //!end j
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer - 1; j <= jmax - refLayer + 1; j++) //!all voxels in y direction
			{
				for (int k = refLayer - 1; k <= kmax - refLayer + 1; k++) //!all voxels in z direction
				{
					//!outlet refined region
					for (int i = imax - refLayer + 2; i <= imax; i++) //!all fine inner voxels in refined top region
					{
						for (int dir = 0; dir < 19; dir++)
						{
							double distrTemp = lattice_[level][i][j][k]->getDistribution(negAlternatingL1, dir); //!initialize temporary distribution
							lattice_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
							double psiTemp = lattice_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
							lattice_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
						}
					} //!end i
				} //!end k
			} //!end j
		} //!end omp

		//!collide in inner fine voxels (no collision in fine interface voxels!)
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!all voxels in x direction
			{
				for (int j = 0; j <= jmax; j++) //!all voxels in y direction
				{
					//!bottom refined region
					for (int k = 0; k <= refLayer - 3; k++) //!all fine inner voxels in refined bottom region
					{
						lattice_[level][i][j][k]->collideRR(bc_, negAlternatingL1); //!perform collide in voxel with hybrid recursive-regularized collision operator (HRR)
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!all voxels in x direction
			{
				for (int j = 0; j <= jmax; j++) //!all voxels in y direction
				{
					//!top refined region
					for (int k = kmax - refLayer + 3; k <= kmax; k++) //!all fine inner voxels in refined top region
					{
						lattice_[level][i][j][k]->collideRR(bc_, negAlternatingL1); //!perform collide in voxel with hybrid recursive-regularized collision operator (HRR)
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!all voxels in x direction
			{
				for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction
				{
					//!left refined region
					for (int j = 0; j <= refLayer - 3; j++) //!all fine inner voxels in refined left region
					{
						lattice_[level][i][j][k]->collideRR(bc_, negAlternatingL1); //!perform collide in voxel with hybrid recursive-regularized collision operator (HRR)
					} //!end j
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!all voxels in x direction
			{
				for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction
				{
					//!right refined region
					for (int j = jmax - refLayer + 3; j <= jmax; j++) //!all fine inner voxels in refined right region
					{
						lattice_[level][i][j][k]->collideRR(bc_, negAlternatingL1); //!perform collide in voxel with hybrid recursive-regularized collision operator (HRR)
					} //!end j
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer - 2; j <= jmax - refLayer + 2; j++) //!all voxels in y direction
			{
				for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction
				{
					//!inlet refined region
					for (int i = 0; i <= refLayerX - 3; i++) //!all fine inner voxels in refined bottom region
					{
						lattice_[level][i][j][k]->collideRR(bc_, negAlternatingL1); //!perform collide in voxel with hybrid recursive-regularized collision operator (HRR)
					} //!end i
				} //!end k
			} //!end j
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer - 2; j <= jmax - refLayer + 2; j++) //!all voxels in y direction
			{
				for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction
				{
					//!outlet refined region
					for (int i = imax - refLayer + 3; i <= imax; i++) //!all fine inner voxels in refined top region
					{
						lattice_[level][i][j][k]->collideRR(bc_, negAlternatingL1); //!perform collide in voxel with hybrid recursive-regularized collision operator (HRR)
					} //!end i
				} //!end k
			} //!end j
		} //!end omp

		//!transport in all fine voxels
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!all voxels in x direction
			{
				for (int j = 0; j <= jmax; j++) //!all voxels in y direction except walls
				{
					//!bottom region
					for (int k = 0; k <= refLayer - 1; k++) //!all fine voxels in refined bottom region except walls
					{
						lattice_[level][i][j][k]->transport(alternatingL1);
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!all voxels in x direction
			{
				for (int j = 0; j <= jmax; j++) //!all voxels in y direction except walls
				{
					//!top region
					for (int k = kmax - refLayer + 1; k <= kmax; k++) //!all fine voxels in refined top region except walls
					{
						lattice_[level][i][j][k]->transport(alternatingL1);
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!all voxels in x direction
			{
				for (int k = refLayer; k <= kmax - refLayer; k++) //!all voxels in z direction except walls and voxels considered before
				{
					//!left region
					for (int j = 0; j <= refLayer - 1; j++) //!all fine voxels in refined left region
					{
						lattice_[level][i][j][k]->transport(alternatingL1);
					} //!end k
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!all voxels in x direction
			{
				for (int k = refLayer; k <= kmax - refLayer; k++) //!all voxels in z direction except walls and voxels considered before
				{
					//!right region
					for (int j = jmax - refLayer + 1; j <= jmax; j++) //!all fine voxels in refined right region
					{
						lattice_[level][i][j][k]->transport(alternatingL1);
					} //!end j
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer; j <= jmax - refLayer; j++) //!all voxels in y direction except except walls and voxels considered before
			{
				for (int k = refLayer; k <= kmax - refLayer; k++) //!all voxels in z direction except walls and voxels considered before
				{
					//!inlet region
					for (int i = 0; i <= refLayerX - 1; i++) //!all fine voxels in refined inlet region except boundary
					{
						lattice_[level][i][j][k]->transport(alternatingL1);
					} //!end i
				} //!end k
			} //!end j
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer; j <= jmax - refLayer; j++) //!all voxels in y direction except except walls and voxels considered before
			{
				for (int k = refLayer; k <= kmax - refLayer; k++) //!all voxels in z direction except walls and voxels considered before
				{
					//!outlet region
					for (int i = imax - refLayer + 1; i <= imax; i++) //!all fine voxels in refined outlet region except boundary
					{
						lattice_[level][i][j][k]->transport(alternatingL1);
					} //!end i
				} //!end k
			} //!end j
		} //!end omp

		//!level 0
		level = 0;
		refLayer = bc_.getRefinementLayer();
		refLayerX = bc_.getRefinementLayerX() * pow(2, level);
		imax = bc_.getNumberOfVoxelsX() - 1; //!index of last voxel in x direction in ->this level
		jmax = bc_.getNumberOfVoxelsY() - 1; //!index of last voxel in y direction in ->this level
		kmax = bc_.getNumberOfVoxelsZ() - 1; //!index of last voxel in z direction in ->this level

		//!coalesce states in coarse interface voxels
		//!ATTENTION: COARSE EDGE AND CORNER VOXELS NEED SEPARATE TREATMENT
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int i = refLayerX; i <= imax - refLayer; i++) //!all voxels in x direction
			{
				for (int j = refLayer; j <= jmax - refLayer; j++) //!all voxels in y direction
				{
					//!coalesce states in coarse bottom interface voxels
					int k = refLayer - 1;
					lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceBottom);

					//!coalesce states in coarse top interface voxels
					k = kmax - refLayer + 1;
					lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceTop);
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int i = refLayerX; i <= imax - refLayer; i++) //!all coarse interface voxels in x direction
			{
				for (int k = refLayer; k <= kmax - refLayer; k++) //!all coarse interface voxels in z direction
				{
					//!coalesce states in coarse right interface voxels
					int j = refLayer - 1;
					lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceLeft);

					//!coalesce states in coarse left interface voxels
					j = jmax - refLayer + 1;
					lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceRight);
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int j = refLayer; j <= jmax - refLayer; j++) //!all coarse interface voxels in y direction
			{
				for (int k = refLayer; k <= kmax - refLayer; k++) //!all coarse interface voxels in z direction
				{
					//!coalesce states in coarse inlet interface voxels
					int i = refLayerX - 1;
					lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceInlet);

					//!coalesce states in coarse outlet interface voxels
					i = imax - refLayer + 1;
					lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceOutlet);
				} //!end k
			} //!end i
		} //!end omp

		/*COARSE INTERFACE EDGES (NO CORNERS!)*/
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel;
			//!coalesce fine states in coarse interface edge voxels along x direction edges
			for (int i = refLayerX; i <= imax - refLayer; i++)
			{
				//!bottom left edge
				int j = refLayer - 1;
				int k = refLayer - 1;
				lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceBottomLeft);

				//!bottom right edge
				j = jmax - refLayer + 1;
				k = refLayer - 1;
				lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceBottomRight);

				//!top left edge
				j = refLayer - 1;
				k = kmax - refLayer + 1;
				lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceTopLeft);

				//!top right edge
				j = jmax - refLayer + 1;
				k = kmax - refLayer + 1;
				lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceTopRight);
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel;
			//!coalesce fine states in coarse interface edge voxels along y direction edges
			for (int j = refLayer; j <= jmax - refLayer; j++)
			{
				//!inlet bottom edge
				int i = refLayerX - 1;
				int k = refLayer - 1;
				lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceInletBottom);

				//!inlet top edge
				i = refLayerX - 1;
				k = kmax - refLayer + 1;
				lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceInletTop);

				//!outlet bottom edge
				i = imax - refLayer + 1;
				k = refLayer - 1;
				lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceOutletBottom);

				//!outlet top edge
				i = imax - refLayer + 1;
				k = kmax - refLayer + 1;
				lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceOutletTop);
			} //!end j
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel;
			//!coalesce fine states in coarse interface edge voxels along z direction edges
			for (int k = refLayer; k <= kmax - refLayer; k++)
			{
				//!inlet left edge
				int i = refLayerX - 1;
				int j = refLayer - 1;
				lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceInletLeft);

				//!inlet right edge
				i = refLayerX - 1;
				j = jmax - refLayer + 1;
				lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceInletRight);

				//!outlet left edge
				i = imax - refLayer + 1;
				j = refLayer - 1;
				lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceOutletLeft);

				//!outlet right edge
				i = imax - refLayer + 1;
				j = jmax - refLayer + 1;
				lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceOutletRight);
			} //!end k
		} //!end omp

		/*COARSE INTERFACE CORNERS*/
		//!inlet bottom left corner
		int i = refLayerX - 1;
		int j = refLayer - 1;
		int k = refLayer - 1;
		lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceInletBottomLeft);

		//!inlet bottom right corner
		i = refLayerX - 1;
		j = jmax - refLayer + 1;
		k = refLayer - 1;
		lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceInletBottomRight);

		//!inlet top left corner
		i = refLayerX - 1;
		j = refLayer - 1;
		k = kmax - refLayer + 1;
		lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceInletTopLeft);

		//!inlet top right corner
		i = refLayerX - 1;
		j = jmax - refLayer + 1;
		k = kmax - refLayer + 1;
		lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceInletTopRight);

		//!outlet bottom left corner
		i = imax - refLayer + 1;
		j = refLayer - 1;
		k = refLayer - 1;
		lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceOutletBottomLeft);

		//!outlet bottom right corner
		i = imax - refLayer + 1;
		j = jmax - refLayer + 1;
		k = refLayer - 1;
		lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceOutletBottomRight);

		//!outlet top left corner
		i = imax - refLayer + 1;
		j = refLayer - 1;
		k = kmax - refLayer + 1;
		lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceOutletTopLeft);

		//!outlet top right corner
		i = imax - refLayer + 1;
		j = jmax - refLayer + 1;
		k = kmax - refLayer + 1;
		lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceOutletTopRight);
	}
	else
	{
		int refLayer, refLayerX; //!number of refined cells
		bool negAlternatingL1 = !alternatingL1;
		bool negAlternatingL0 = !alternatingL0;

		//!STRING VARIABLES FOR DIFFERENCIATING BETWEEN VARIOUS INTERFACE VOXELS -> NEEDED FOR EXPLODE AND COALESCE METHOD
		//!'core' interface voxels
		string interfaceBottom = "interfaceBottom";
		string interfaceTop = "interfaceTop";
		string interfaceLeft = "interfaceLeft";
		string interfaceRight = "interfaceRight";
		string interfaceInlet = "interfaceInlet";
		string interfaceOutlet = "interfaceOutlet";

		//!'edge' interface voxels
		string interfaceBottomLeft = "interfaceBottomLeft";
		string interfaceBottomRight = "interfaceBottomRight";
		string interfaceTopLeft = "interfaceTopLeft";
		string interfaceTopRight = "interfaceTopRight";
		string interfaceInletBottom = "interfaceInletBottom";
		string interfaceInletTop = "interfaceInletTop";
		string interfaceOutletBottom = "interfaceOutletBottom";
		string interfaceOutletTop = "interfaceOutletTop";
		string interfaceInletLeft = "interfaceInletLeft";
		string interfaceInletRight = "interfaceInletRight";
		string interfaceOutletLeft = "interfaceOutletLeft";
		string interfaceOutletRight = "interfaceOutletRight";

		//!'corner' interface voxels
		string interfaceInletBottomLeft = "interfaceInletBottomLeft";
		string interfaceInletBottomRight = "interfaceInletBottomRight";
		string interfaceInletTopLeft = "interfaceInletTopLeft";
		string interfaceInletTopRight = "interfaceInletTopRight";
		string interfaceOutletBottomLeft = "interfaceOutletBottomLeft";
		string interfaceOutletBottomRight = "interfaceOutletBottomRight";
		string interfaceOutletTopLeft = "interfaceOutletTopLeft";
		string interfaceOutletTopRight = "interfaceOutletTopRight";

		//--------------------------------------------------//
		//--------------------------------------------------//
		//             collision and transport				//
		//--------------------------------------------------//
		//--------------------------------------------------//

		//!current distribution (t=0) level 0: alternating
		//!current distribution (t=0) level 1: alternating

		int level = 0;
		int imax = bc_.getNumberOfVoxelsX() - 1; //!index of last voxel in x direction in ->this level
		int jmax = bc_.getNumberOfVoxelsY() - 1; //!index of last voxel in y direction in ->this level
		int kmax = bc_.getNumberOfVoxelsZ() - 1; //!index of last voxel in z direction in ->this level
		refLayer = bc_.getRefinementLayer(); //!number of refined coarse voxels
		refLayerX = bc_.getRefinementLayerX(); //!number of refined coarse voxels

	#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
	#pragma omp for //!do outer loop in parallel
			for (int i = refLayerX - 1; i <= imax - refLayer + 1; i++) //!all nodes in x direction
			{
				for (int j = refLayer - 1; j <= jmax - refLayer + 1; j++) //!all nodes in y direction
				{
					for (int k = refLayer - 1; k <= kmax - refLayer + 1; k++) //!nodes in z direction from coarse lower interface to upper interface
					{
						for (int dir = 0; dir < 19; dir++)
						{
							double distrTemp = lattice_[level][i][j][k]->getDistribution(alternatingL0, dir); //!initialize temporary distribution
							lattice_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
						}
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

		//!collide in all coarse voxels
	#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
	#pragma omp for //!do outer loop in parallel
			for (int i = refLayerX - 1; i <= imax - refLayer + 1; i++) //!all nodes in x direction
			{
				for (int j = refLayer - 1; j <= jmax - refLayer + 1; j++) //!all nodes in y direction
				{
					for (int k = refLayer - 1; k <= kmax - refLayer + 1; k++) //!nodes in z direction from coarse lower interface to upper interface
					{
						lattice_[level][i][j][k]->collideRR(bc_, alternatingL0); //!perform collide in voxel with recursive-regularized collision operator (RR)
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

		level = 1;
		refLayer = bc_.getRefinementLayer() * pow(2, level);
		refLayerX = bc_.getRefinementLayerX() * pow(2, level);
		imax = bc_.getNumberOfVoxelsX() * pow(2., level) - 1; //!index of last voxel in x direction in ->this level
		jmax = bc_.getNumberOfVoxelsY() * pow(2., level) - 1; //!index of last voxel in y direction in ->this level
		kmax = bc_.getNumberOfVoxelsZ() * pow(2., level) - 1; //!index of last voxel in z direction in ->this level

		//!collide in inner fine voxels (no collision in fine interface voxels!)
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!all voxels in x direction
			{
				for (int j = 0; j <= jmax; j++) //!all voxels in y direction
				{
					//!bottom refined region
					for (int k = 0; k <= refLayer - 3; k++) //!all fine inner voxels in refined bottom region
					{
						lattice_[level][i][j][k]->collideRR(bc_, alternatingL1); //!perform collide in voxel with recursive-regularized collision operator (RR)
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!all voxels in x direction
			{
				for (int j = 0; j <= jmax; j++) //!all voxels in y direction
				{
					//!top refined region
					for (int k = kmax - refLayer + 3; k <= kmax; k++) //!all fine inner voxels in refined top region
					{
						lattice_[level][i][j][k]->collideRR(bc_, alternatingL1); //!perform collide in voxel with recursive-regularized collision operator (RR)
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!all voxels in x direction
			{
				for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction
				{
					//!left refined region
					for (int j = 0; j <= refLayer - 3; j++) //!all fine inner voxels in refined left region
					{
						lattice_[level][i][j][k]->collideRR(bc_, alternatingL1); //!perform collide in voxel with recursive-regularized collision operator (RR)
					} //!end j
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!all voxels in x direction
			{
				for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction
				{
					//!right refined region
					for (int j = jmax - refLayer + 3; j <= jmax; j++) //!all fine inner voxels in refined right region
					{
						lattice_[level][i][j][k]->collideRR(bc_, alternatingL1); //!perform collide in voxel with recursive-regularized collision operator (RR)
					} //!end j
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer - 2; j <= jmax - refLayer + 2; j++) //!all voxels in y direction
			{
				for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction
				{
					//!inlet refined region
					for (int i = 0; i <= refLayerX - 3; i++) //!all fine inner voxels in refined bottom region
					{
						lattice_[level][i][j][k]->collideRR(bc_, alternatingL1); //!perform collide in voxel with recursive-regularized collision operator (RR)
					} //!end i
				} //!end k
			} //!end j
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer - 2; j <= jmax - refLayer + 2; j++) //!all voxels in y direction
			{
				for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction
				{
					//!outlet refined region
					for (int i = imax - refLayer + 3; i <= imax; i++) //!all fine inner voxels in refined top region
					{
						lattice_[level][i][j][k]->collideRR(bc_, alternatingL1); //!perform collide in voxel with recursive-regularized collision operator (RR)
					} //!end i
				} //!end k
			} //!end j
		} //!end omp

		level = 0;
		imax = bc_.getNumberOfVoxelsX() - 1; //!index of last voxel in x direction in ->this level
		jmax = bc_.getNumberOfVoxelsY() - 1; //!index of last voxel in y direction in ->this level
		kmax = bc_.getNumberOfVoxelsZ() - 1; //!index of last voxel in z direction in ->this level
		refLayer = bc_.getRefinementLayer(); //!number of refined coarse voxels
		refLayerX = bc_.getRefinementLayerX() * pow(2, level);

		//!explode states in coarse interface voxels
		//!INTERFACE EDGES AND CORNERS INCLUDED (ALL DISTRIBUTIONS ARE BEING EXPLODED!)
		if (bc_.getExplosionOrder() == 1)
		{
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
			{
#pragma omp for //!do outer loop in parallel
				for (int i = refLayerX - 1; i <= imax - refLayer + 1; i++) //!all coarse interface voxels in x direction
				{
					for (int j = refLayer - 1; j <= jmax - refLayer + 1; j++) //!all coarse interface voxels in y direction
					{
						//!explode states in coarse bottom interface voxels
						int k = refLayer - 1;
						lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

						//!explode states in coarse top interface voxels
						k = kmax - refLayer + 1;
						lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);
					} //!end j
				} //!end i
			} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
			{
#pragma omp for //!do outer loop in parallel
				for (int i = refLayerX - 1; i <= imax - refLayer + 1; i++) //!all coarse interface voxels in x direction
				{
					for (int k = refLayer; k <= kmax - refLayer; k++) //!all coarse interface voxels in z direction
					{
						//!explode states in coarse right interface voxels
						int j = refLayer - 1;
						lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

						//!explode states in coarse left interface voxels
						j = jmax - refLayer + 1;
						lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);
					} //!end k
				} //!end i
			} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
			{
#pragma omp for //!do outer loop in parallel
				for (int j = refLayer; j <= jmax - refLayer; j++) //!all coarse interface voxels in y direction
				{
					for (int k = refLayer; k <= kmax - refLayer; k++) //!all coarse interface voxels in z direction
					{
						//!explode states in coarse inlet interface voxels
						int i = refLayerX - 1;
						lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

						//!explode states in coarse outlet interface voxels
						i = imax - refLayer + 1;
						lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);
					} //!end k
				} //!end i
			} //!end omp
		} //!end if

		else if (bc_.getExplosionOrder() == 2)
		{
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
			{
#pragma omp for //!do outer loop in parallel
				for (int i = refLayerX; i <= imax - refLayer; i++) //!all coarse interface voxels in x direction
				{
					for (int j = refLayer; j <= jmax - refLayer; j++) //!all coarse interface voxels in y direction
					{
						//!explode states in coarse bottom interface voxels
						int k = refLayer - 1;
						lattice_[level][i][j][k]->explodeLinear(alternatingL1, alternatingL0, bc_, interfaceBottom);

						//!explode states in coarse top interface voxels
						k = kmax - refLayer + 1;
						lattice_[level][i][j][k]->explodeLinear(alternatingL1, alternatingL0, bc_, interfaceTop);
					} //!end j
				} //!end i
			} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
			{
#pragma omp for //!do outer loop in parallel
				for (int i = refLayerX; i <= imax - refLayer; i++) //!all coarse interface voxels in x direction
				{
					for (int k = refLayer; k <= kmax - refLayer; k++) //!all coarse interface voxels in z direction
					{
						//!explode states in coarse right interface voxels
						int j = refLayer - 1;
						lattice_[level][i][j][k]->explodeLinear(alternatingL1, alternatingL0, bc_, interfaceLeft);

						//!explode states in coarse left interface voxels
						j = jmax - refLayer + 1;
						lattice_[level][i][j][k]->explodeLinear(alternatingL1, alternatingL0, bc_, interfaceRight);
					} //!end k
				} //!end i
			} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
			{
#pragma omp for //!do outer loop in parallel
				for (int j = refLayer; j <= jmax - refLayer; j++) //!all coarse interface voxels in y direction
				{
					for (int k = refLayer; k <= kmax - refLayer; k++) //!all coarse interface voxels in z direction
					{
						//!explode states in coarse inlet interface voxels
						int i = refLayerX - 1;
						lattice_[level][i][j][k]->explodeLinear(alternatingL1, alternatingL0, bc_, interfaceInlet);

						//!explode states in coarse outlet interface voxels
						i = imax - refLayer + 1;
						lattice_[level][i][j][k]->explodeLinear(alternatingL1, alternatingL0, bc_, interfaceOutlet);
					} //!end k
				} //!end i
			} //!end omp

			/*COARSE INTERFACE EDGES*/
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
			{
#pragma omp for //!do most outer loop in parallel;

				//!explode coarse states on interface edge voxels along x direction edges
				for (int i = refLayerX - 1; i <= imax - refLayer + 1; i++)
				{
					//!bottom left edge
					int j = refLayer - 1;
					int k = refLayer - 1;
					lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

					//!bottom right edge
					j = jmax - refLayer + 1;
					k = refLayer - 1;
					lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

					//!top left edge
					j = refLayer - 1;
					k = kmax - refLayer + 1;
					lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

					//!top right edge
					j = jmax - refLayer + 1;
					k = kmax - refLayer + 1;
					lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);
				} //!end i
			} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
			{
#pragma omp for //!do most outer loop in parallel;

				//!explode coarse states on interface edge voxels along y direction edges
				for (int j = refLayer; j <= jmax - refLayer; j++)
				{
					//!inlet bottom edge
					int i = refLayerX - 1;
					int k = refLayer - 1;
					lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

					//!inlet top edge
					i = refLayerX - 1;
					k = kmax - refLayer + 1;
					lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

					//!outlet bottom edge
					i = imax - refLayer + 1;
					k = refLayer - 1;
					lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

					//!outlet top edge
					i = imax - refLayer + 1;
					k = kmax - refLayer + 1;
					lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);
				} //!end j
			} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
			{
#pragma omp for //!do most outer loop in parallel;

				//!explode coarse states on interface edge voxels along z direction edges
				for (int k = refLayer; k <= kmax - refLayer; k++)
				{
					//!inlet left edge
					int i = refLayerX - 1;
					int j = refLayer - 1;
					lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

					//!inlet right edge
					i = refLayerX - 1;
					j = jmax - refLayer + 1;
					lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

					//!outlet left edge
					i = imax - refLayer + 1;
					j = refLayer - 1;
					lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

					//!outlet right edge
					i = imax - refLayer + 1;
					j = jmax - refLayer + 1;
					lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);
				} //!end k
			} //!end omp
		} //!end else if

		//!transport in all coarse voxels
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX - 1; i <= imax - refLayer + 1; i++) //!all voxels in x direction
			{
				for (int j = refLayer - 1; j <= jmax - refLayer + 1; j++) //!all voxels in y direction
				{
					for (int k = refLayer - 1; k <= kmax - refLayer + 1; k++)
					{
						lattice_[level][i][j][k]->transport(negAlternatingL0);
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

		//!transport in all fine voxels
		level = 1;
		refLayer = bc_.getRefinementLayer() * pow(2, level);
		refLayerX = bc_.getRefinementLayerX() * pow(2, level);
		imax = bc_.getNumberOfVoxelsX() * pow(2., level) - 1; //!index of last voxel in x direction in ->this level
		jmax = bc_.getNumberOfVoxelsY() * pow(2., level) - 1; //!index of last voxel in y direction in ->this level
		kmax = bc_.getNumberOfVoxelsZ() * pow(2., level) - 1; //!index of last voxel in z direction in ->this level

		//!transport in all fine voxels
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!all voxels in x direction
			{
				for (int j = 0; j <= jmax; j++) //!all voxels in y direction except walls
				{
					//!bottom region
					for (int k = 0; k <= refLayer - 1; k++) //!all fine voxels in refined bottom region except walls
					{
						lattice_[level][i][j][k]->transport(negAlternatingL1);
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!all voxels in x direction
			{
				for (int j = 0; j <= jmax; j++) //!all voxels in y direction except walls
				{
					//!top region
					for (int k = kmax - refLayer + 1; k <= kmax; k++) //!all fine voxels in refined top region except walls
					{
						lattice_[level][i][j][k]->transport(negAlternatingL1);
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!all voxels in x direction
			{
				for (int k = refLayer; k <= kmax - refLayer; k++) //!all voxels in z direction except walls and voxels considered before
				{
					//!left region
					for (int j = 0; j <= refLayer - 1; j++) //!all fine voxels in refined left region
					{
						lattice_[level][i][j][k]->transport(negAlternatingL1);
					} //!end k
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!all voxels in x direction
			{
				for (int k = refLayer; k <= kmax - refLayer; k++) //!all voxels in z direction except walls and voxels considered before
				{
					//!right region
					for (int j = jmax - refLayer + 1; j <= jmax; j++) //!all fine voxels in refined right region
					{
						lattice_[level][i][j][k]->transport(negAlternatingL1);
					} //!end j
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer; j <= jmax - refLayer; j++) //!all voxels in y direction except except walls and voxels considered before
			{
				for (int k = refLayer; k <= kmax - refLayer; k++) //!all voxels in z direction except walls and voxels considered before
				{
					//!inlet region
					for (int i = 0; i <= refLayerX - 1; i++) //!all fine voxels in refined inlet region except boundary
					{
						lattice_[level][i][j][k]->transport(negAlternatingL1);
					} //!end i
				} //!end k
			} //!end j
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer; j <= jmax - refLayer; j++) //!all voxels in y direction except except walls and voxels considered before
			{
				for (int k = refLayer; k <= kmax - refLayer; k++) //!all voxels in z direction except walls and voxels considered before
				{
					//!outlet region
					for (int i = imax - refLayer + 1; i <= imax; i++) //!all fine voxels in refined outlet region except boundary
					{
						lattice_[level][i][j][k]->transport(negAlternatingL1);
					} //!end i
				} //!end k
			} //!end j
		} //!end omp

		//!collide in inner fine voxels (no collision in fine interface voxels!)
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!all voxels in x direction
			{
				for (int j = 0; j <= jmax; j++) //!all voxels in y direction
				{
					//!bottom refined region
					for (int k = 0; k <= refLayer - 3; k++) //!all fine inner voxels in refined bottom region
					{
						lattice_[level][i][j][k]->collideRR(bc_, negAlternatingL1); //!perform collide in voxel with recursive-regularized collision operator (RR)
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!all voxels in x direction
			{
				for (int j = 0; j <= jmax; j++) //!all voxels in y direction
				{
					//!top refined region
					for (int k = kmax - refLayer + 3; k <= kmax; k++) //!all fine inner voxels in refined top region
					{
						lattice_[level][i][j][k]->collideRR(bc_, negAlternatingL1); //!perform collide in voxel with recursive-regularized collision operator (RR)
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!all voxels in x direction
			{
				for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction
				{
					//!left refined region
					for (int j = 0; j <= refLayer - 3; j++) //!all fine inner voxels in refined left region
					{
						lattice_[level][i][j][k]->collideRR(bc_, negAlternatingL1); //!perform collide in voxel with recursive-regularized collision operator (RR)
					} //!end j
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!all voxels in x direction
			{
				for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction
				{
					//!right refined region
					for (int j = jmax - refLayer + 3; j <= jmax; j++) //!all fine inner voxels in refined right region
					{
						lattice_[level][i][j][k]->collideRR(bc_, negAlternatingL1); //!perform collide in voxel with recursive-regularized collision operator (RR)
					} //!end j
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer - 2; j <= jmax - refLayer + 2; j++) //!all voxels in y direction
			{
				for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction
				{
					//!inlet refined region
					for (int i = 0; i <= refLayerX - 3; i++) //!all fine inner voxels in refined bottom region
					{
						lattice_[level][i][j][k]->collideRR(bc_, negAlternatingL1); //!perform collide in voxel with recursive-regularized collision operator (RR)
					} //!end i
				} //!end k
			} //!end j
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer - 2; j <= jmax - refLayer + 2; j++) //!all voxels in y direction
			{
				for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction
				{
					//!outlet refined region
					for (int i = imax - refLayer + 3; i <= imax; i++) //!all fine inner voxels in refined top region
					{
						lattice_[level][i][j][k]->collideRR(bc_, negAlternatingL1); //!perform collide in voxel with recursive-regularized collision operator (RR)
					} //!end i
				} //!end k
			} //!end j
		} //!end omp

	//!transport in all fine voxels
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!all voxels in x direction
			{
				for (int j = 0; j <= jmax; j++) //!all voxels in y direction except walls
				{
					//!bottom region
					for (int k = 0; k <= refLayer - 1; k++) //!all fine voxels in refined bottom region except walls
					{
						lattice_[level][i][j][k]->transport(alternatingL1);
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!all voxels in x direction
			{
				for (int j = 0; j <= jmax; j++) //!all voxels in y direction except walls
				{
					//!top region
					for (int k = kmax - refLayer + 1; k <= kmax; k++) //!all fine voxels in refined top region except walls
					{
						lattice_[level][i][j][k]->transport(alternatingL1);
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!all voxels in x direction
			{
				for (int k = refLayer; k <= kmax - refLayer; k++) //!all voxels in z direction except walls and voxels considered before
				{
					//!left region
					for (int j = 0; j <= refLayer - 1; j++) //!all fine voxels in refined left region
					{
						lattice_[level][i][j][k]->transport(alternatingL1);
					} //!end k
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!all voxels in x direction
			{
				for (int k = refLayer; k <= kmax - refLayer; k++) //!all voxels in z direction except walls and voxels considered before
				{
					//!right region
					for (int j = jmax - refLayer + 1; j <= jmax; j++) //!all fine voxels in refined right region
					{
						lattice_[level][i][j][k]->transport(alternatingL1);
					} //!end j
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer; j <= jmax - refLayer; j++) //!all voxels in y direction except except walls and voxels considered before
			{
				for (int k = refLayer; k <= kmax - refLayer; k++) //!all voxels in z direction except walls and voxels considered before
				{
					//!inlet region
					for (int i = 0; i <= refLayerX - 1; i++) //!all fine voxels in refined inlet region except boundary
					{
						lattice_[level][i][j][k]->transport(alternatingL1);
					} //!end i
				} //!end k
			} //!end j
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer; j <= jmax - refLayer; j++) //!all voxels in y direction except except walls and voxels considered before
			{
				for (int k = refLayer; k <= kmax - refLayer; k++) //!all voxels in z direction except walls and voxels considered before
				{
					//!outlet region
					for (int i = imax - refLayer + 1; i <= imax; i++) //!all fine voxels in refined outlet region except boundary
					{
						lattice_[level][i][j][k]->transport(alternatingL1);
					} //!end i
				} //!end k
			} //!end j
		} //!end omp

		//!level 0
		level = 0;
		refLayer = bc_.getRefinementLayer();
		refLayerX = bc_.getRefinementLayerX() * pow(2, level);
		imax = bc_.getNumberOfVoxelsX() - 1; //!index of last voxel in x direction in ->this level
		jmax = bc_.getNumberOfVoxelsY() - 1; //!index of last voxel in y direction in ->this level
		kmax = bc_.getNumberOfVoxelsZ() - 1; //!index of last voxel in z direction in ->this level

		//!coalesce states in coarse interface voxels
		//!ATTENTION: COARSE EDGE AND CORNER VOXELS NEED SEPARATE TREATMENT
	#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
	#pragma omp for //!do outer loop in parallel
			for (int i = refLayerX; i <= imax - refLayer; i++) //!all voxels in x direction
			{
				for (int j = refLayer; j <= jmax - refLayer; j++) //!all voxels in y direction
				{
					//!coalesce states in coarse bottom interface voxels
					int k = refLayer - 1;
					lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceBottom);

					//!coalesce states in coarse top interface voxels
					k = kmax - refLayer + 1;
					lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceTop);
				} //!end j
			} //!end i
		} //!end omp

	#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
	#pragma omp for //!do outer loop in parallel
			for (int i = refLayerX; i <= imax - refLayer; i++) //!all coarse interface voxels in x direction
			{
				for (int k = refLayer; k <= kmax - refLayer; k++) //!all coarse interface voxels in z direction
				{
					//!coalesce states in coarse right interface voxels
					int j = refLayer - 1;
					lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceLeft);

					//!coalesce states in coarse left interface voxels
					j = jmax - refLayer + 1;
					lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceRight);
				} //!end k
			} //!end i
		} //!end omp

	#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
	#pragma omp for //!do outer loop in parallel
			for (int j = refLayer; j <= jmax - refLayer; j++) //!all coarse interface voxels in y direction
			{
				for (int k = refLayer; k <= kmax - refLayer; k++) //!all coarse interface voxels in z direction
				{
					//!coalesce states in coarse inlet interface voxels
					int i = refLayerX - 1;
					lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceInlet);

					//!coalesce states in coarse outlet interface voxels
					i = imax - refLayer + 1;
					lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceOutlet);
				} //!end k
			} //!end i
		} //!end omp

		/*COARSE INTERFACE EDGES (NO CORNERS!)*/
	#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
	#pragma omp for //!do most outer loop in parallel;
			//!coalesce fine states in coarse interface edge voxels along x direction edges
			for (int i = refLayerX; i <= imax - refLayer; i++)
			{
				//!bottom left edge
				int j = refLayer - 1;
				int k = refLayer - 1;
				lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceBottomLeft);

				//!bottom right edge
				j = jmax - refLayer + 1;
				k = refLayer - 1;
				lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceBottomRight);

				//!top left edge
				j = refLayer - 1;
				k = kmax - refLayer + 1;
				lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceTopLeft);

				//!top right edge
				j = jmax - refLayer + 1;
				k = kmax - refLayer + 1;
				lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceTopRight);
			} //!end i
		} //!end omp

	#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
	#pragma omp for //!do most outer loop in parallel;
			//!coalesce fine states in coarse interface edge voxels along y direction edges
			for (int j = refLayer; j <= jmax - refLayer; j++)
			{
				//!inlet bottom edge
				int i = refLayerX - 1;
				int k = refLayer - 1;
				lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceInletBottom);

				//!inlet top edge
				i = refLayerX - 1;
				k = kmax - refLayer + 1;
				lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceInletTop);

				//!outlet bottom edge
				i = imax - refLayer + 1;
				k = refLayer - 1;
				lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceOutletBottom);

				//!outlet top edge
				i = imax - refLayer + 1;
				k = kmax - refLayer + 1;
				lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceOutletTop);
			} //!end j
		} //!end omp

	#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
	#pragma omp for //!do most outer loop in parallel;
			//!coalesce fine states in coarse interface edge voxels along z direction edges
			for (int k = refLayer; k <= kmax - refLayer; k++)
			{
				//!inlet left edge
				int i = refLayerX - 1;
				int j = refLayer - 1;
				lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceInletLeft);

				//!inlet right edge
				i = refLayerX - 1;
				j = jmax - refLayer + 1;
				lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceInletRight);

				//!outlet left edge
				i = imax - refLayer + 1;
				j = refLayer - 1;
				lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceOutletLeft);

				//!outlet right edge
				i = imax - refLayer + 1;
				j = jmax - refLayer + 1;
				lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceOutletRight);
			} //!end k
		} //!end omp

		/*COARSE INTERFACE CORNERS*/
		//!inlet bottom left corner
		int i = refLayerX - 1;
		int j = refLayer - 1;
		int k = refLayer - 1;
		lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceInletBottomLeft);

		//!inlet bottom right corner
		i = refLayerX - 1;
		j = jmax - refLayer + 1;
		k = refLayer - 1;
		lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceInletBottomRight);

		//!inlet top left corner
		i = refLayerX - 1;
		j = refLayer - 1;
		k = kmax - refLayer + 1;
		lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceInletTopLeft);

		//!inlet top right corner
		i = refLayerX - 1;
		j = jmax - refLayer + 1;
		k = kmax - refLayer + 1;
		lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceInletTopRight);

		//!outlet bottom left corner
		i = imax - refLayer + 1;
		j = refLayer - 1;
		k = refLayer - 1;
		lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceOutletBottomLeft);

		//!outlet bottom right corner
		i = imax - refLayer + 1;
		j = jmax - refLayer + 1;
		k = refLayer - 1;
		lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceOutletBottomRight);

		//!outlet top left corner
		i = imax - refLayer + 1;
		j = refLayer - 1;
		k = kmax - refLayer + 1;
		lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceOutletTopLeft);

		//!outlet top right corner
		i = imax - refLayer + 1;
		j = jmax - refLayer + 1;
		k = kmax - refLayer + 1;
		lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceOutletTopRight);
	}
} //!end timeStepNestedRR

//!nested time stepping algorithm
void Lattice::nestedTimeSteppingHRR(bool& alternatingL1, bool& alternatingL0, Control& bc_) //!perform nested time step with sequence of collide stream step for all level
{
	int refLayer, refLayerX; //!number of refined cells
	bool negAlternatingL1 = !alternatingL1;
	bool negAlternatingL0 = !alternatingL0;

	//!STRING VARIABLES FOR DIFFERENCIATING BETWEEN VARIOUS INTERFACE VOXELS -> NEEDED FOR EXPLODE AND COALESCE METHOD
	//!'core' interface voxels
	string interfaceBottom = "interfaceBottom";
	string interfaceTop = "interfaceTop";
	string interfaceLeft = "interfaceLeft";
	string interfaceRight = "interfaceRight";
	string interfaceInlet = "interfaceInlet";
	string interfaceOutlet = "interfaceOutlet";

	//!'edge' interface voxels
	string interfaceBottomLeft = "interfaceBottomLeft";
	string interfaceBottomRight = "interfaceBottomRight";
	string interfaceTopLeft = "interfaceTopLeft";
	string interfaceTopRight = "interfaceTopRight";
	string interfaceInletBottom = "interfaceInletBottom";
	string interfaceInletTop = "interfaceInletTop";
	string interfaceOutletBottom = "interfaceOutletBottom";
	string interfaceOutletTop = "interfaceOutletTop";
	string interfaceInletLeft = "interfaceInletLeft";
	string interfaceInletRight = "interfaceInletRight";
	string interfaceOutletLeft = "interfaceOutletLeft";
	string interfaceOutletRight = "interfaceOutletRight";

	//!'corner' interface voxels
	string interfaceInletBottomLeft = "interfaceInletBottomLeft";
	string interfaceInletBottomRight = "interfaceInletBottomRight";
	string interfaceInletTopLeft = "interfaceInletTopLeft";
	string interfaceInletTopRight = "interfaceInletTopRight";
	string interfaceOutletBottomLeft = "interfaceOutletBottomLeft";
	string interfaceOutletBottomRight = "interfaceOutletBottomRight";
	string interfaceOutletTopLeft = "interfaceOutletTopLeft";
	string interfaceOutletTopRight = "interfaceOutletTopRight";

	//--------------------------------------------------//
	//--------------------------------------------------//
	//             collision and transport				//
	//--------------------------------------------------//
	//--------------------------------------------------//

	//!current distribution (t=0) level 0: alternating
	//!current distribution (t=0) level 1: alternating

	int level = 1;
	refLayer = bc_.getRefinementLayer() * pow(2, level);
	refLayerX = bc_.getRefinementLayerX() * pow(2, level);
	int imax = bc_.getNumberOfVoxelsX() * pow(2., level) - 1; //!index of last voxel in x direction in ->this level
	int jmax = bc_.getNumberOfVoxelsY() * pow(2., level) - 1; //!index of last voxel in y direction in ->this level
	int kmax = bc_.getNumberOfVoxelsZ() * pow(2., level) - 1; //!index of last voxel in z direction in ->this level

	//!safe pre-collide distributions for HRR collision operator on fine inner voxels
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all voxels in y direction
			{
				//!bottom refined region
				for (int k = 0; k <= refLayer - 3; k++) //!all fine inner voxels in refined bottom region
				{
					for (int dir = 0; dir < 19; dir++)
					{
						double distrTemp = lattice_[level][i][j][k]->getDistribution(alternatingL1, dir); //!initialize temporary distribution
						lattice_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
						double psiTemp = lattice_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
						lattice_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
					}
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all voxels in y direction
			{
				//!top refined region
				for (int k = kmax - refLayer + 3; k <= kmax; k++) //!all fine inner voxels in refined top region
				{
					for (int dir = 0; dir < 19; dir++)
					{
						double distrTemp = lattice_[level][i][j][k]->getDistribution(alternatingL1, dir); //!initialize temporary distribution
						lattice_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
						double psiTemp = lattice_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
						lattice_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
					}
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction
			{
				//!left refined region
				for (int j = 0; j <= refLayer - 3; j++) //!all fine inner voxels in refined left region
				{
					for (int dir = 0; dir < 19; dir++)
					{
						double distrTemp = lattice_[level][i][j][k]->getDistribution(alternatingL1, dir); //!initialize temporary distribution
						lattice_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
						double psiTemp = lattice_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
						lattice_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
					}
				} //!end j
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction
			{
				//!right refined region
				for (int j = jmax - refLayer + 3; j <= jmax; j++) //!all fine inner voxels in refined right region
				{
					for (int dir = 0; dir < 19; dir++)
					{
						double distrTemp = lattice_[level][i][j][k]->getDistribution(alternatingL1, dir); //!initialize temporary distribution
						lattice_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
						double psiTemp = lattice_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
						lattice_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
					}
				} //!end j
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer - 2; j <= jmax - refLayer + 2; j++) //!all voxels in y direction
		{
			for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction
			{
				//!inlet refined region
				for (int i = 0; i <= refLayerX - 3; i++) //!all fine inner voxels in refined bottom region
				{
					for (int dir = 0; dir < 19; dir++)
					{
						double distrTemp = lattice_[level][i][j][k]->getDistribution(alternatingL1, dir); //!initialize temporary distribution
						lattice_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
						double psiTemp = lattice_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
						lattice_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
					}
				} //!end i
			} //!end k
		} //!end j
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer - 2; j <= jmax - refLayer + 2; j++) //!all voxels in y direction
		{
			for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction
			{
				//!outlet refined region
				for (int i = imax - refLayer + 3; i <= imax; i++) //!all fine inner voxels in refined top region
				{
					for (int dir = 0; dir < 19; dir++)
					{
						double distrTemp = lattice_[level][i][j][k]->getDistribution(alternatingL1, dir); //!initialize temporary distribution
						lattice_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
						double psiTemp = lattice_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
						lattice_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
					}
				} //!end i
			} //!end k
		} //!end j
	} //!end omp

	level = 0;
	imax = bc_.getNumberOfVoxelsX() - 1; //!index of last voxel in x direction in ->this level
	jmax = bc_.getNumberOfVoxelsY() - 1; //!index of last voxel in y direction in ->this level
	kmax = bc_.getNumberOfVoxelsZ() - 1; //!index of last voxel in z direction in ->this level
	refLayer = bc_.getRefinementLayer(); //!number of refined coarse voxels
	refLayerX = bc_.getRefinementLayerX() * pow(2, level);

	//!safe pre-collide distributions for HRR collision operator on coarse voxels
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do outer loop in parallel
		for (int i = refLayerX - 1; i <= imax - refLayer + 1; i++) //!all nodes in x direction
		{
			for (int j = refLayer - 1; j <= jmax - refLayer + 1; j++) //!all nodes in y direction
			{
				for (int k = refLayer - 1; k <= kmax - refLayer + 1; k++) //!nodes in z direction from coarse lower interface to upper interface
				{
					for (int dir = 0; dir < 19; dir++)
					{
						double distrTemp = lattice_[level][i][j][k]->getDistribution(alternatingL0, dir); //!initialize temporary distribution
						lattice_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
						double psiTemp = lattice_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
						lattice_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
					}
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

	//!collide in all coarse voxels
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do outer loop in parallel
		for (int i = refLayerX - 1; i <= imax - refLayer + 1; i++) //!all nodes in x direction
		{
			for (int j = refLayer - 1; j <= jmax - refLayer + 1; j++) //!all nodes in y direction
			{
				for (int k = refLayer - 1; k <= kmax - refLayer + 1; k++) //!nodes in z direction from coarse lower interface to upper interface
				{
					lattice_[level][i][j][k]->collideHRR(bc_, alternatingL0); //!perform collide in voxel with hybrid recursive-regularized collision operator (HRR)
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

	//!explode states in coarse interface voxels
	//!INTERFACE EDGES AND CORNERS INCLUDED (ALL DISTRIBUTIONS ARE BEING EXPLODED!)
	if (bc_.getExplosionOrder() == 1)
	{
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int i = refLayerX - 1; i <= imax - refLayer + 1; i++) //!all coarse interface voxels in x direction
			{
				for (int j = refLayer - 1; j <= jmax - refLayer + 1; j++) //!all coarse interface voxels in y direction
				{
					//!explode states in coarse bottom interface voxels
					int k = refLayer - 1;
					lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

					//!explode states in coarse top interface voxels
					k = kmax - refLayer + 1;
					lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int i = refLayerX - 1; i <= imax - refLayer + 1; i++) //!all coarse interface voxels in x direction
			{
				for (int k = refLayer; k <= kmax - refLayer; k++) //!all coarse interface voxels in z direction
				{
					//!explode states in coarse right interface voxels
					int j = refLayer - 1;
					lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

					//!explode states in coarse left interface voxels
					j = jmax - refLayer + 1;
					lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int j = refLayer; j <= jmax - refLayer; j++) //!all coarse interface voxels in y direction
			{
				for (int k = refLayer; k <= kmax - refLayer; k++) //!all coarse interface voxels in z direction
				{
					//!explode states in coarse inlet interface voxels
					int i = refLayerX - 1;
					lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

					//!explode states in coarse outlet interface voxels
					i = imax - refLayer + 1;
					lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);
				} //!end k
			} //!end i
		} //!end omp
	} //!end if

	else if (bc_.getExplosionOrder() == 2)
	{
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int i = refLayerX; i <= imax - refLayer; i++) //!all coarse interface voxels in x direction
			{
				for (int j = refLayer; j <= jmax - refLayer; j++) //!all coarse interface voxels in y direction
				{
					//!explode states in coarse bottom interface voxels
					int k = refLayer - 1;
					lattice_[level][i][j][k]->explodeLinear(alternatingL1, alternatingL0, bc_, interfaceBottom);

					//!explode states in coarse top interface voxels
					k = kmax - refLayer + 1;
					lattice_[level][i][j][k]->explodeLinear(alternatingL1, alternatingL0, bc_, interfaceTop);
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int i = refLayerX; i <= imax - refLayer; i++) //!all coarse interface voxels in x direction
			{
				for (int k = refLayer; k <= kmax - refLayer; k++) //!all coarse interface voxels in z direction
				{
					//!explode states in coarse right interface voxels
					int j = refLayer - 1;
					lattice_[level][i][j][k]->explodeLinear(alternatingL1, alternatingL0, bc_, interfaceLeft);

					//!explode states in coarse left interface voxels
					j = jmax - refLayer + 1;
					lattice_[level][i][j][k]->explodeLinear(alternatingL1, alternatingL0, bc_, interfaceRight);
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int j = refLayer; j <= jmax - refLayer; j++) //!all coarse interface voxels in y direction
			{
				for (int k = refLayer; k <= kmax - refLayer; k++) //!all coarse interface voxels in z direction
				{
					//!explode states in coarse inlet interface voxels
					int i = refLayerX - 1;
					lattice_[level][i][j][k]->explodeLinear(alternatingL1, alternatingL0, bc_, interfaceInlet);

					//!explode states in coarse outlet interface voxels
					i = imax - refLayer + 1;
					lattice_[level][i][j][k]->explodeLinear(alternatingL1, alternatingL0, bc_, interfaceOutlet);
				} //!end k
			} //!end i
		} //!end omp

		/*COARSE INTERFACE EDGES*/
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel;

			//!explode coarse states on interface edge voxels along x direction edges
			for (int i = refLayerX - 1; i <= imax - refLayer + 1; i++)
			{
				//!bottom left edge
				int j = refLayer - 1;
				int k = refLayer - 1;
				lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

				//!bottom right edge
				j = jmax - refLayer + 1;
				k = refLayer - 1;
				lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

				//!top left edge
				j = refLayer - 1;
				k = kmax - refLayer + 1;
				lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

				//!top right edge
				j = jmax - refLayer + 1;
				k = kmax - refLayer + 1;
				lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel;

			//!explode coarse states on interface edge voxels along y direction edges
			for (int j = refLayer; j <= jmax - refLayer; j++)
			{
				//!inlet bottom edge
				int i = refLayerX - 1;
				int k = refLayer - 1;
				lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

				//!inlet top edge
				i = refLayerX - 1;
				k = kmax - refLayer + 1;
				lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

				//!outlet bottom edge
				i = imax - refLayer + 1;
				k = refLayer - 1;
				lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

				//!outlet top edge
				i = imax - refLayer + 1;
				k = kmax - refLayer + 1;
				lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);
			} //!end j
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel;

			//!explode coarse states on interface edge voxels along z direction edges
			for (int k = refLayer; k <= kmax - refLayer; k++)
			{
				//!inlet left edge
				int i = refLayerX - 1;
				int j = refLayer - 1;
				lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

				//!inlet right edge
				i = refLayerX - 1;
				j = jmax - refLayer + 1;
				lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

				//!outlet left edge
				i = imax - refLayer + 1;
				j = refLayer - 1;
				lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);

				//!outlet right edge
				i = imax - refLayer + 1;
				j = jmax - refLayer + 1;
				lattice_[level][i][j][k]->explode(alternatingL1, alternatingL0);
			} //!end k
		} //!end omp
	} //!end else if

	//!transport in all coarse voxels
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX - 1; i <= imax - refLayer + 1; i++) //!all voxels in x direction
		{
			for (int j = refLayer - 1; j <= jmax - refLayer + 1; j++) //!all voxels in y direction
			{
				for (int k = refLayer - 1; k <= kmax - refLayer + 1; k++)
				{
					lattice_[level][i][j][k]->transport(negAlternatingL0);
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

	level = 1;
	refLayer = bc_.getRefinementLayer() * pow(2, level);
	refLayerX = bc_.getRefinementLayerX() * pow(2, level);
	imax = bc_.getNumberOfVoxelsX() * pow(2., level) - 1; //!index of last voxel in x direction in ->this level
	jmax = bc_.getNumberOfVoxelsY() * pow(2., level) - 1; //!index of last voxel in y direction in ->this level
	kmax = bc_.getNumberOfVoxelsZ() * pow(2., level) - 1; //!index of last voxel in z direction in ->this level

	//!safe pre-collide distributions for HRR collision operator on fine first-layer voxels after explode
	//!for inner fine voxels states have already been saved...
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX - 2; i <= imax - refLayer + 2; i++) //!all voxels in x direction
		{
			for (int j = refLayer - 2; j <= jmax - refLayer + 2; j++) //!all voxels in y direction
			{
				//!bottom refined region
				int k = refLayer - 2;
				for (int dir = 0; dir < 19; dir++)
				{
					double distrTemp = lattice_[level][i][j][k]->getDistribution(alternatingL1, dir); //!initialize temporary distribution
					lattice_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
					double psiTemp = lattice_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
					lattice_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
				}
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX - 2; i <= imax - refLayer + 2; i++) //!all voxels in x direction
		{
			for (int j = refLayer - 2; j <= jmax - refLayer + 2; j++) //!all voxels in y direction
			{
				//!top refined region
				int k = kmax - refLayer + 2;
				for (int dir = 0; dir < 19; dir++)
				{
					double distrTemp = lattice_[level][i][j][k]->getDistribution(alternatingL1, dir); //!initialize temporary distribution
					lattice_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
					double psiTemp = lattice_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
					lattice_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
				}
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX - 2; i <= imax - refLayer + 2; i++) //!all voxels in x direction
		{
			for (int k = refLayer - 1; k <= kmax - refLayer + 1; k++) //!all voxels in z direction
			{
				//!left refined region
				int j = refLayer - 2;
				for (int dir = 0; dir < 19; dir++)
				{
					double distrTemp = lattice_[level][i][j][k]->getDistribution(alternatingL1, dir); //!initialize temporary distribution
					lattice_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
					double psiTemp = lattice_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
					lattice_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
				}
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX - 2; i <= imax - refLayer + 2; i++) //!all voxels in x direction
		{
			for (int k = refLayer - 1; k <= kmax - refLayer + 1; k++) //!all voxels in z direction
			{
				//!right refined region
				int j = jmax - refLayer + 2;
				for (int dir = 0; dir < 19; dir++)
				{
					double distrTemp = lattice_[level][i][j][k]->getDistribution(alternatingL1, dir); //!initialize temporary distribution
					lattice_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
					double psiTemp = lattice_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
					lattice_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
				}
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer - 1; j <= jmax - refLayer + 1; j++) //!all voxels in y direction
		{
			for (int k = refLayer - 1; k <= kmax - refLayer + 1; k++) //!all voxels in z direction
			{
				//!inlet refined region
				int i = refLayerX - 2;
				for (int dir = 0; dir < 19; dir++)
				{
					double distrTemp = lattice_[level][i][j][k]->getDistribution(alternatingL1, dir); //!initialize temporary distribution
					lattice_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
					double psiTemp = lattice_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
					lattice_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
				}
			} //!end k
		} //!end j
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer - 1; j <= jmax - refLayer + 1; j++) //!all voxels in y direction
		{
			for (int k = refLayer - 1; k <= kmax - refLayer + 1; k++) //!all voxels in z direction
			{
				//!inlet refined region
				int i = imax - refLayer + 2;
				for (int dir = 0; dir < 19; dir++)
				{
					double distrTemp = lattice_[level][i][j][k]->getDistribution(alternatingL1, dir); //!initialize temporary distribution
					lattice_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
					double psiTemp = lattice_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
					lattice_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
				}
			} //!end k
		} //!end j
	} //!end omp


	//!collide in inner fine voxels (no collision in fine interface voxels!)
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all voxels in y direction
			{
				//!bottom refined region
				for (int k = 0; k <= refLayer - 3; k++) //!all fine inner voxels in refined bottom region
				{
					lattice_[level][i][j][k]->collideHRR(bc_, alternatingL1); //!perform collide in voxel with hybrid recursive-regularized collision operator (HRR)
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all voxels in y direction
			{
				//!top refined region
				for (int k = kmax - refLayer + 3; k <= kmax; k++) //!all fine inner voxels in refined top region
				{
					lattice_[level][i][j][k]->collideHRR(bc_, alternatingL1); //!perform collide in voxel with hybrid recursive-regularized collision operator (HRR)
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction
			{
				//!left refined region
				for (int j = 0; j <= refLayer - 3; j++) //!all fine inner voxels in refined left region
				{
					lattice_[level][i][j][k]->collideHRR(bc_, alternatingL1); //!perform collide in voxel with hybrid recursive-regularized collision operator (HRR)
				} //!end j
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction
			{
				//!right refined region
				for (int j = jmax - refLayer + 3; j <= jmax; j++) //!all fine inner voxels in refined right region
				{
					lattice_[level][i][j][k]->collideHRR(bc_, alternatingL1); //!perform collide in voxel with hybrid recursive-regularized collision operator (HRR)
				} //!end j
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer - 2; j <= jmax - refLayer + 2; j++) //!all voxels in y direction
		{
			for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction
			{
				//!inlet refined region
				for (int i = 0; i <= refLayerX - 3; i++) //!all fine inner voxels in refined bottom region
				{
					lattice_[level][i][j][k]->collideHRR(bc_, alternatingL1); //!perform collide in voxel with hybrid recursive-regularized collision operator (HRR)
				} //!end i
			} //!end k
		} //!end j
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer - 2; j <= jmax - refLayer + 2; j++) //!all voxels in y direction
		{
			for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction
			{
				//!outlet refined region
				for (int i = imax - refLayer + 3; i <= imax; i++) //!all fine inner voxels in refined top region
				{
					lattice_[level][i][j][k]->collideHRR(bc_, alternatingL1); //!perform collide in voxel with hybrid recursive-regularized collision operator (HRR)
				} //!end i
			} //!end k
		} //!end j
	} //!end omp

	//!transport in all fine voxels
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all voxels in y direction except walls
			{
				//!bottom region
				for (int k = 0; k <= refLayer - 1; k++) //!all fine voxels in refined bottom region except walls
				{
					lattice_[level][i][j][k]->transport(negAlternatingL1);
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all voxels in y direction except walls
			{
				//!top region
				for (int k = kmax - refLayer + 1; k <= kmax; k++) //!all fine voxels in refined top region except walls
				{
					lattice_[level][i][j][k]->transport(negAlternatingL1);
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int k = refLayer; k <= kmax - refLayer; k++) //!all voxels in z direction except walls and voxels considered before
			{
				//!left region
				for (int j = 0; j <= refLayer - 1; j++) //!all fine voxels in refined left region
				{
					lattice_[level][i][j][k]->transport(negAlternatingL1);
				} //!end k
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int k = refLayer; k <= kmax - refLayer; k++) //!all voxels in z direction except walls and voxels considered before
			{
				//!right region
				for (int j = jmax - refLayer + 1; j <= jmax; j++) //!all fine voxels in refined right region
				{
					lattice_[level][i][j][k]->transport(negAlternatingL1);
				} //!end j
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer; j <= jmax - refLayer; j++) //!all voxels in y direction except except walls and voxels considered before
		{
			for (int k = refLayer; k <= kmax - refLayer; k++) //!all voxels in z direction except walls and voxels considered before
			{
				//!inlet region
				for (int i = 0; i <= refLayerX - 1; i++) //!all fine voxels in refined inlet region except boundary
				{
					lattice_[level][i][j][k]->transport(negAlternatingL1);
				} //!end i
			} //!end k
		} //!end j
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer; j <= jmax - refLayer; j++) //!all voxels in y direction except except walls and voxels considered before
		{
			for (int k = refLayer; k <= kmax - refLayer; k++) //!all voxels in z direction except walls and voxels considered before
			{
				//!outlet region
				for (int i = imax - refLayer + 1; i <= imax; i++) //!all fine voxels in refined outlet region except boundary
				{
					lattice_[level][i][j][k]->transport(negAlternatingL1);
				} //!end i
			} //!end k
		} //!end j
	} //!end omp

	//!safe pre-collide distributions for HRR collision operator on fine inner and first-layer voxels
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all voxels in y direction
			{
				//!bottom refined region
				for (int k = 0; k <= refLayer - 2; k++) //!all fine inner voxels in refined bottom region
				{
					for (int dir = 0; dir < 19; dir++)
					{
						double distrTemp = lattice_[level][i][j][k]->getDistribution(negAlternatingL1, dir); //!initialize temporary distribution
						lattice_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
						double psiTemp = lattice_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
						lattice_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
					}
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all voxels in y direction
			{
				//!top refined region
				for (int k = kmax - refLayer + 2; k <= kmax; k++) //!all fine inner voxels in refined top region
				{
					for (int dir = 0; dir < 19; dir++)
					{
						double distrTemp = lattice_[level][i][j][k]->getDistribution(negAlternatingL1, dir); //!initialize temporary distribution
						lattice_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
						double psiTemp = lattice_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
						lattice_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
					}
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int k = refLayer - 1; k <= kmax - refLayer + 1; k++) //!all voxels in z direction
			{
				//!left refined region
				for (int j = 0; j <= refLayer - 2; j++) //!all fine inner voxels in refined left region
				{
					for (int dir = 0; dir < 19; dir++)
					{
						double distrTemp = lattice_[level][i][j][k]->getDistribution(negAlternatingL1, dir); //!initialize temporary distribution
						lattice_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
						double psiTemp = lattice_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
						lattice_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
					}
				} //!end j
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int k = refLayer - 1; k <= kmax - refLayer + 1; k++) //!all voxels in z direction
			{
				//!right refined region
				for (int j = jmax - refLayer + 2; j <= jmax; j++) //!all fine inner voxels in refined right region
				{
					for (int dir = 0; dir < 19; dir++)
					{
						double distrTemp = lattice_[level][i][j][k]->getDistribution(negAlternatingL1, dir); //!initialize temporary distribution
						lattice_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
						double psiTemp = lattice_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
						lattice_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
					}
				} //!end j
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer - 1; j <= jmax - refLayer + 1; j++) //!all voxels in y direction
		{
			for (int k = refLayer - 1; k <= kmax - refLayer + 1; k++) //!all voxels in z direction
			{
				//!inlet refined region
				for (int i = 0; i <= refLayerX - 2; i++) //!all fine inner voxels in refined bottom region
				{
					for (int dir = 0; dir < 19; dir++)
					{
						double distrTemp = lattice_[level][i][j][k]->getDistribution(negAlternatingL1, dir); //!initialize temporary distribution
						lattice_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
						double psiTemp = lattice_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
						lattice_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
					}
				} //!end i
			} //!end k
		} //!end j
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer - 1; j <= jmax - refLayer + 1; j++) //!all voxels in y direction
		{
			for (int k = refLayer - 1; k <= kmax - refLayer + 1; k++) //!all voxels in z direction
			{
				//!outlet refined region
				for (int i = imax - refLayer + 2; i <= imax; i++) //!all fine inner voxels in refined top region
				{
					for (int dir = 0; dir < 19; dir++)
					{
						double distrTemp = lattice_[level][i][j][k]->getDistribution(negAlternatingL1, dir); //!initialize temporary distribution
						lattice_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
						double psiTemp = lattice_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
						lattice_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
					}
				} //!end i
			} //!end k
		} //!end j
	} //!end omp

		//!collide in inner fine voxels (no collision in fine interface voxels!)
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all voxels in y direction
			{
				//!bottom refined region
				for (int k = 0; k <= refLayer - 3; k++) //!all fine inner voxels in refined bottom region
				{
					lattice_[level][i][j][k]->collideHRR(bc_, negAlternatingL1); //!perform collide in voxel with hybrid recursive-regularized collision operator (HRR)
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all voxels in y direction
			{
				//!top refined region
				for (int k = kmax - refLayer + 3; k <= kmax; k++) //!all fine inner voxels in refined top region
				{
					lattice_[level][i][j][k]->collideHRR(bc_, negAlternatingL1); //!perform collide in voxel with hybrid recursive-regularized collision operator (HRR)
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction
			{
				//!left refined region
				for (int j = 0; j <= refLayer - 3; j++) //!all fine inner voxels in refined left region
				{
					lattice_[level][i][j][k]->collideHRR(bc_, negAlternatingL1); //!perform collide in voxel with hybrid recursive-regularized collision operator (HRR)
				} //!end j
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction
			{
				//!right refined region
				for (int j = jmax - refLayer + 3; j <= jmax; j++) //!all fine inner voxels in refined right region
				{
					lattice_[level][i][j][k]->collideHRR(bc_, negAlternatingL1); //!perform collide in voxel with hybrid recursive-regularized collision operator (HRR)
				} //!end j
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer - 2; j <= jmax - refLayer + 2; j++) //!all voxels in y direction
		{
			for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction
			{
				//!inlet refined region
				for (int i = 0; i <= refLayerX - 3; i++) //!all fine inner voxels in refined bottom region
				{
					lattice_[level][i][j][k]->collideHRR(bc_, negAlternatingL1); //!perform collide in voxel with hybrid recursive-regularized collision operator (HRR)
				} //!end i
			} //!end k
		} //!end j
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer - 2; j <= jmax - refLayer + 2; j++) //!all voxels in y direction
		{
			for (int k = refLayer - 2; k <= kmax - refLayer + 2; k++) //!all voxels in z direction
			{
				//!outlet refined region
				for (int i = imax - refLayer + 3; i <= imax; i++) //!all fine inner voxels in refined top region
				{
					lattice_[level][i][j][k]->collideHRR(bc_, negAlternatingL1); //!perform collide in voxel with hybrid recursive-regularized collision operator (HRR)
				} //!end i
			} //!end k
		} //!end j
	} //!end omp

	//!transport in all fine voxels
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all voxels in y direction except walls
			{
				//!bottom region
				for (int k = 0; k <= refLayer - 1; k++) //!all fine voxels in refined bottom region except walls
				{
					lattice_[level][i][j][k]->transport(alternatingL1);
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all voxels in y direction except walls
			{
				//!top region
				for (int k = kmax - refLayer + 1; k <= kmax; k++) //!all fine voxels in refined top region except walls
				{
					lattice_[level][i][j][k]->transport(alternatingL1);
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int k = refLayer; k <= kmax - refLayer; k++) //!all voxels in z direction except walls and voxels considered before
			{
				//!left region
				for (int j = 0; j <= refLayer - 1; j++) //!all fine voxels in refined left region
				{
					lattice_[level][i][j][k]->transport(alternatingL1);
				} //!end k
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all voxels in x direction
		{
			for (int k = refLayer; k <= kmax - refLayer; k++) //!all voxels in z direction except walls and voxels considered before
			{
				//!right region
				for (int j = jmax - refLayer + 1; j <= jmax; j++) //!all fine voxels in refined right region
				{
					lattice_[level][i][j][k]->transport(alternatingL1);
				} //!end j
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer; j <= jmax - refLayer; j++) //!all voxels in y direction except except walls and voxels considered before
		{
			for (int k = refLayer; k <= kmax - refLayer; k++) //!all voxels in z direction except walls and voxels considered before
			{
				//!inlet region
				for (int i = 0; i <= refLayerX - 1; i++) //!all fine voxels in refined inlet region except boundary
				{
					lattice_[level][i][j][k]->transport(alternatingL1);
				} //!end i
			} //!end k
		} //!end j
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer; j <= jmax - refLayer; j++) //!all voxels in y direction except except walls and voxels considered before
		{
			for (int k = refLayer; k <= kmax - refLayer; k++) //!all voxels in z direction except walls and voxels considered before
			{
				//!outlet region
				for (int i = imax - refLayer + 1; i <= imax; i++) //!all fine voxels in refined outlet region except boundary
				{
					lattice_[level][i][j][k]->transport(alternatingL1);
				} //!end i
			} //!end k
		} //!end j
	} //!end omp

	//!level 0
	level = 0;
	refLayer = bc_.getRefinementLayer();
	refLayerX = bc_.getRefinementLayerX() * pow(2, level);
	imax = bc_.getNumberOfVoxelsX() - 1; //!index of last voxel in x direction in ->this level
	jmax = bc_.getNumberOfVoxelsY() - 1; //!index of last voxel in y direction in ->this level
	kmax = bc_.getNumberOfVoxelsZ() - 1; //!index of last voxel in z direction in ->this level

	//!coalesce states in coarse interface voxels
	//!ATTENTION: COARSE EDGE AND CORNER VOXELS NEED SEPARATE TREATMENT
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do outer loop in parallel
		for (int i = refLayerX; i <= imax - refLayer; i++) //!all voxels in x direction
		{
			for (int j = refLayer; j <= jmax - refLayer; j++) //!all voxels in y direction
			{
				//!coalesce states in coarse bottom interface voxels
				int k = refLayer - 1;
				lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceBottom);

				//!coalesce states in coarse top interface voxels
				k = kmax - refLayer + 1;
				lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceTop);
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do outer loop in parallel
		for (int i = refLayerX; i <= imax - refLayer; i++) //!all coarse interface voxels in x direction
		{
			for (int k = refLayer; k <= kmax - refLayer; k++) //!all coarse interface voxels in z direction
			{
				//!coalesce states in coarse right interface voxels
				int j = refLayer - 1;
				lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceLeft);

				//!coalesce states in coarse left interface voxels
				j = jmax - refLayer + 1;
				lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceRight);
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do outer loop in parallel
		for (int j = refLayer; j <= jmax - refLayer; j++) //!all coarse interface voxels in y direction
		{
			for (int k = refLayer; k <= kmax - refLayer; k++) //!all coarse interface voxels in z direction
			{
				//!coalesce states in coarse inlet interface voxels
				int i = refLayerX - 1;
				lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceInlet);

				//!coalesce states in coarse outlet interface voxels
				i = imax - refLayer + 1;
				lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceOutlet);
			} //!end k
		} //!end i
	} //!end omp

	/*COARSE INTERFACE EDGES (NO CORNERS!)*/
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel;
		//!coalesce fine states in coarse interface edge voxels along x direction edges
		for (int i = refLayerX; i <= imax - refLayer; i++)
		{
			//!bottom left edge
			int j = refLayer - 1;
			int k = refLayer - 1;
			lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceBottomLeft);

			//!bottom right edge
			j = jmax - refLayer + 1;
			k = refLayer - 1;
			lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceBottomRight);

			//!top left edge
			j = refLayer - 1;
			k = kmax - refLayer + 1;
			lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceTopLeft);

			//!top right edge
			j = jmax - refLayer + 1;
			k = kmax - refLayer + 1;
			lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceTopRight);
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel;
		//!coalesce fine states in coarse interface edge voxels along y direction edges
		for (int j = refLayer; j <= jmax - refLayer; j++)
		{
			//!inlet bottom edge
			int i = refLayerX - 1;
			int k = refLayer - 1;
			lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceInletBottom);

			//!inlet top edge
			i = refLayerX - 1;
			k = kmax - refLayer + 1;
			lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceInletTop);

			//!outlet bottom edge
			i = imax - refLayer + 1;
			k = refLayer - 1;
			lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceOutletBottom);

			//!outlet top edge
			i = imax - refLayer + 1;
			k = kmax - refLayer + 1;
			lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceOutletTop);
		} //!end j
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel;
		//!coalesce fine states in coarse interface edge voxels along z direction edges
		for (int k = refLayer; k <= kmax - refLayer; k++)
		{
			//!inlet left edge
			int i = refLayerX - 1;
			int j = refLayer - 1;
			lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceInletLeft);

			//!inlet right edge
			i = refLayerX - 1;
			j = jmax - refLayer + 1;
			lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceInletRight);

			//!outlet left edge
			i = imax - refLayer + 1;
			j = refLayer - 1;
			lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceOutletLeft);

			//!outlet right edge
			i = imax - refLayer + 1;
			j = jmax - refLayer + 1;
			lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceOutletRight);
		} //!end k
	} //!end omp

	/*COARSE INTERFACE CORNERS*/
	//!inlet bottom left corner
	int i = refLayerX - 1;
	int j = refLayer - 1;
	int k = refLayer - 1;
	lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceInletBottomLeft);

	//!inlet bottom right corner
	i = refLayerX - 1;
	j = jmax - refLayer + 1;
	k = refLayer - 1;
	lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceInletBottomRight);

	//!inlet top left corner
	i = refLayerX - 1;
	j = refLayer - 1;
	k = kmax - refLayer + 1;
	lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceInletTopLeft);

	//!inlet top right corner
	i = refLayerX - 1;
	j = jmax - refLayer + 1;
	k = kmax - refLayer + 1;
	lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceInletTopRight);

	//!outlet bottom left corner
	i = imax - refLayer + 1;
	j = refLayer - 1;
	k = refLayer - 1;
	lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceOutletBottomLeft);

	//!outlet bottom right corner
	i = imax - refLayer + 1;
	j = jmax - refLayer + 1;
	k = refLayer - 1;
	lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceOutletBottomRight);

	//!outlet top left corner
	i = imax - refLayer + 1;
	j = refLayer - 1;
	k = kmax - refLayer + 1;
	lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceOutletTopLeft);

	//!outlet top right corner
	i = imax - refLayer + 1;
	j = jmax - refLayer + 1;
	k = kmax - refLayer + 1;
	lattice_[level][i][j][k]->coalesce(alternatingL1, negAlternatingL0, interfaceOutletTopRight);
} //!end timeStepNestedHRR


void Lattice::writeResultsAtVoxel(string filename, bool& alternatingL0, bool& alternatingL1, Control& bc_)
{
    fstream results;
    results.open(filename, std::fstream::in | std::fstream::out | std::fstream::app);

    //!results are written after every coarse cycle, so sampling frequency depends on coarse time step
    int level0 = 0;
    int level1 = 1;
    int MPnum = bc_.getMPnum();
    double timestep = bc_.getTimeStep();
    double radius = bc_.getMPradius();
    double z = bc_.getMPZ();
    double* phi = new double[MPnum];

    int** MPindex = new int*[MPnum];
    for (int index = 0; index < MPnum; index++)
    {
        MPindex[index] = new int[3];
        //!calculate node indices for each monitor point from coordinates
        phi[index] = index * (M_PI/(0.5 * MPnum)) + M_PI/MPnum;
        double x = radius * cos(phi[index]) + bc_.getChannelSizeX_()/2.;
        double y = radius * sin(phi[index]) + bc_.getChannelSizeY_()/2.;

        if(x < bc_.getRefinementLayerX()*bc_.getSpacing())
        {
            //!monitor points on the fine grid
            MPindex[index][2] = (((z / bc_.getSpacing()) * pow(2., level1 + 1.)) - 1.) / 2.;
            MPindex[index][1] = (((y / bc_.getSpacing()) * pow(2., level1 + 1.)) - 1.) / 2.;
            MPindex[index][0] = (((x / bc_.getSpacing()) * pow(2., level1 + 1.)) - 1.) / 2.;
        }

        else if(x >= bc_.getRefinementLayerX()*bc_.getSpacing())
        {
            //!monitor points on the coarse grid
            MPindex[index][2] = (((z / bc_.getSpacing()) * pow(2., level0 + 1.)) - 1.) / 2.;
            MPindex[index][1] = (((y / bc_.getSpacing()) * pow(2., level0 + 1.)) - 1.) / 2.;
            MPindex[index][0] = (((x / bc_.getSpacing()) * pow(2., level0 + 1.)) - 1.) / 2.;
        }
        else
        {
            cerr << "Something seems to be wrong with the monitor point coordinates." << flush;
            exit(1);
        }
    }

    //!find nodes
    Voxel** nodes = new Voxel*[MPnum];

    for(int index = 0; index < MPnum; index++)
    {
        double x = radius * cos(phi[index]) + bc_.getChannelSizeX_()/2.;

        if(x < bc_.getRefinementLayerX()*bc_.getSpacing())
        {
            nodes[index] = lattice_[level1][MPindex[index][0]][MPindex[index][1]][MPindex[index][2]];
        }
        else if(x >= bc_.getRefinementLayerX()*bc_.getSpacing())
        {
            nodes[index] = lattice_[level0][MPindex[index][0]][MPindex[index][1]][MPindex[index][2]];
        }
    }

    //!if file does not exist, create new file
    //!write header
    if ((bc_.getTime_() == 0) && results)
    {
        results << "time step in level 0 [s] = " << timestep << endl;
        results << "number of time steps" << "\t" << "(x, y, z)" << "\t"<< "density [kgm-3]" << "\t" << "pressure [Pa]" << "\t" << "velocity x [ms-1]" << "\t" << "velocity y [ms-1]" << "\t" << "velocity z [ms-1]" << endl << endl;
        results.close();

    }
    else if (bc_.getTime_() > 0)
    {
        results << bc_.getTime_() << "\t";

        for(int index = 0; index < MPnum; index++)
        {
            double x = radius * cos(phi[index]) + bc_.getChannelSizeX_()/2;

            //!write monitor points on fine grid
            if ((bc_.getTime_() > 0) && (x < bc_.getRefinementLayerX()*bc_.getSpacing()))
            {
                double pressure = nodes[index]->calcPressure(alternatingL1, bc_);
                double density = nodes[index]->calcDensity(bc_, alternatingL1);
                double velocity[3] = { 0. };
                nodes[index]->calcVelocity(bc_, alternatingL1, density, velocity);

                results << fixed << setprecision(12) << "(" << nodes[index]->getX() << ", " << nodes[index]->getY() << ", " << nodes[index]->getZ() << ")" << "\t" << density << " " << pressure << " " << velocity[0] << " " << velocity[1] << " " << velocity[2] << "\t";
            }
            //!write monitor points on coarse grid
            else if ((bc_.getTime_() > 0) && (x >= bc_.getRefinementLayerX()*bc_.getSpacing()))
            {
                double pressure = nodes[index]->calcPressure(alternatingL0, bc_);
                double density = nodes[index]->calcDensity(bc_, alternatingL0);
                double velocity[3] = { 0. };
                nodes[index]->calcVelocity(bc_, alternatingL0, density, velocity);

                results << fixed << setprecision(12) << "(" << nodes[index]->getX() << ", " << nodes[index]->getY() << ", " << nodes[index]->getZ() << ")" << "\t" << density << " " << pressure << " " << velocity[0] << " " << velocity[1] << " " << velocity[2] << "\t";
            }
        }
        results << endl;
        results.close();
    }
    else
    {
        cerr << "File seems to be corrupt..." << flush;
        exit(1);
    }

    for(int index = 0; index < MPnum; index++)
    {
        delete[] MPindex[index];
    }
    delete[] MPindex;//!delete arrays

    delete[] nodes;//!delete arrays
    delete[] phi;

} //!end write results


double Lattice::readValue(string filename, string keyword) //!filename: name of initialization file; keyword: string to be read from file; return value: (double) to corresponding keyword
{
	ifstream file;
	string s;
	file.open(filename.c_str());
	istringstream sin;
	string op;
	double number = 0.;

	if (!file)
	{
		cerr << "Error reading initialization file " << filename << "!";
		exit(1);
	}

	while (file.good())
	{
		op = " ";
		getline(file, s);
		sin.clear();
		sin.str(s);
		sin >> op;
		if (op == keyword)
		{
			sin >> number;
			return number;
		}
	}
	cerr << "No Value for keyword \"" + keyword + "\" found in " << filename << "!";
	exit(1);
	return number;
}
