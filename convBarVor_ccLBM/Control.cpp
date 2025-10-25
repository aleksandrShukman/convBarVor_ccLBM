#define _USE_MATH_DEFINES
#include "Control.h"
#include "math.h"
#include <algorithm>
#include <array>

Control::Control(string filename) : timeStepMax_(static_cast<int>(readValue(filename, "Max-Time-Step"))),
									writeInterval_(static_cast<int>(readValue(filename, "WriteInterval"))),
									density_(static_cast<double>(readValue(filename, "Density"))),
									numberOfThreads_(static_cast<int>(readValue(filename, "Threads"))),
									kinematicViscosity_(static_cast<double>(readValue(filename, "KinematicViscosity"))),
									caseName_(readFileName(filename, "Save-as")),
									explosionOrder_(static_cast<int>(readValue(filename, "Explosion-Order"))),
									equilOrder_(static_cast<int>(readValue(filename, "Order-Of-Equilibrium"))),
									collisionModel_(readString(filename, "Collision-Model")),
									cubicMachCorrection_(readString(filename, "Cubic-Mach-Correction")),
									hybridPar_(static_cast<double>(readValue(filename, "Hybrid-Parameter"))),
									spacing_(static_cast<double>(readValue(filename, "Spacing"))),
                                    MPradius_(static_cast<double>(readValue(filename, "MP-Radius"))),
                                    MPZ_(static_cast<double>(readValue(filename, "MP-Z"))),
                                    MPnum_(static_cast<int>(readValue(filename, "Number-Of-MP"))),
                                    channelSizeX_(static_cast<double>(readValue(filename, "Domain-Size-X"))),
									channelSizeY_(static_cast<double>(readValue(filename, "Domain-Size-Y"))),
									channelSizeZ_(static_cast<double>(readValue(filename, "Domain-Size-Z"))),
									refLayer_(static_cast<int>(readValue(filename, "RefinementLayerYZ"))),
                                    refLayerX_(static_cast<int>(readValue(filename, "RefinementLayerX"))),
									levelMax_(static_cast<int>(readValue(filename, "MaxLevel"))),
									machNumber_(static_cast<double>(readValue(filename, "MACH-Number")))
		{
            //!initialization velocity
			readVector(filename, "Init-Velocity", initVelocity_);

			//!calculate number of voxels in cartesian system
			double epsilon  = 1e-8; //!geometrical tolerance
			numberOfVoxelsX_ = ceil(channelSizeX_ / spacing_ + epsilon) - 1;
			numberOfVoxelsY_ = ceil(channelSizeY_ / spacing_ + epsilon) - 1;
			numberOfVoxelsZ_ = ceil(channelSizeZ_ / spacing_ + epsilon) - 1;

			//!calculate number of nodes in cartesian system
			numberOfNodesX_ = numberOfVoxelsX_ + 1;
			numberOfNodesY_ = numberOfVoxelsY_ + 1;
			numberOfNodesZ_ = numberOfVoxelsZ_ + 1;

			//!convected barotropic vortex
            timeStep_ = machNumber_ * spacing_ / (sqrt(3.) * initVelocity_[0]);
            molecularVelocity_ = spacing_ / timeStep_;
            speedOfSound_ = molecularVelocity_ / sqrt(3.);
            collisionFrequency_ = (pow(speedOfSound_, 2) * timeStep_) / (kinematicViscosity_ + 0.5 * pow(speedOfSound_, 2) * timeStep_);
			reynoldsNumber_ = (initVelocity_[0] * 0.06) / kinematicViscosity_; //!Reynolds number based on vortex diameter --> initialization() in Lattice.cpp

			//!lattice parameters (D3Q19)
			xsiX_[0] = molecularVelocity_;
			xsiX_[1] = 0.;
			xsiX_[2] = -molecularVelocity_;
			xsiX_[3] = 0.;
			xsiX_[4] = 0.;
			xsiX_[5] = 0.;
			xsiX_[6] = molecularVelocity_;
			xsiX_[7] = molecularVelocity_;
			xsiX_[8] = 0.;
			xsiX_[9] = 0.;
			xsiX_[10] = -molecularVelocity_;
			xsiX_[11] = -molecularVelocity_;
			xsiX_[12] = 0.;
			xsiX_[13] = 0.;
			xsiX_[14] = -molecularVelocity_;
			xsiX_[15] = molecularVelocity_;
			xsiX_[16] = molecularVelocity_;
			xsiX_[17] = -molecularVelocity_;
			xsiX_[18] = 0.;

			xsiY_[0] = 0.;
			xsiY_[1] = molecularVelocity_;
			xsiY_[2] = 0.;
			xsiY_[3] = -molecularVelocity_;
			xsiY_[4] = 0.;
			xsiY_[5] = 0.;
			xsiY_[6] = 0.;
			xsiY_[7] = 0.;
			xsiY_[8] = molecularVelocity_;
			xsiY_[9] = molecularVelocity_;
			xsiY_[10] = 0.;
			xsiY_[11] = 0.;
			xsiY_[12] = -molecularVelocity_;
			xsiY_[13] = -molecularVelocity_;
			xsiY_[14] = molecularVelocity_;
			xsiY_[15] = -molecularVelocity_;
			xsiY_[16] = molecularVelocity_;
			xsiY_[17] = -molecularVelocity_;
			xsiY_[18] = 0.;

			xsiZ_[0] = 0.;
			xsiZ_[1] = 0.;
			xsiZ_[2] = 0.;
			xsiZ_[3] = 0.;
			xsiZ_[4] = molecularVelocity_;
			xsiZ_[5] = -molecularVelocity_;
			xsiZ_[6] = molecularVelocity_;
			xsiZ_[7] = -molecularVelocity_;
			xsiZ_[8] = molecularVelocity_;
			xsiZ_[9] = -molecularVelocity_;
			xsiZ_[10] = molecularVelocity_;
			xsiZ_[11] = -molecularVelocity_;
			xsiZ_[12] = molecularVelocity_;
			xsiZ_[13] = -molecularVelocity_;
			xsiZ_[14] = 0.;
			xsiZ_[15] = 0.;
			xsiZ_[16] = 0.;
			xsiZ_[17] = 0.;
			xsiZ_[18] = 0.;

			//!set weight factors
			for (unsigned dir = 0; dir < 19; dir++)
			{
				if (dir < 6) weightFactors_[dir]	   = 1. / 18.;
				else if (dir < 18) weightFactors_[dir] = 1. / 36.;
				else  weightFactors_[dir]              = 1. / 3.;
			}
		}


double Control::readValue(string filename,string keyword) //!filename: name of control file; keyword: string to be read from file; return value: (double) to corresponding keyword
{
	ifstream file;
	string s;
	file.open(filename.c_str());
	istringstream sin;
	string op;
	double number = 0.;

	if (!file)
	{
		cerr << "Error reading boundary conditions file " << filename;
		exit(1);
	}

	while (file.good())
	{
		op=" ";
		getline(file, s);
		sin.clear();
		sin.str(s);
		sin >> op ;
		if (op == keyword)
		{
			sin >> number;
			return number;
		}
	}
	cerr << "No Value for keyword \"" + keyword + "\" found in " << filename;
	exit(1);
	return number;
}

string Control::readString(string filename,string keyword)			//!filename: name of control file; keyword: string to be read from file; return value: (string) to corresponding keyword
{
	ifstream file;
	string s;
	file.open(filename.c_str());
	istringstream sin;
	string op;
	string theString;

	if (!file)
	{
		cerr << "Error reading boundary conditions file " << filename; exit(1);
	}

	while (file.good())
	{
		op = " ";
		getline(file, s);
		sin.clear();
		sin.str(s);
		sin >> op ;
		if (op == keyword)
		{
			if (keyword == "file-name")
			{
				string temp = s.substr(keyword.length());
				size_t found = temp.find_first_not_of(" ");
				theString = temp.substr(found);
			}
			else { sin >> theString; }
			return theString;
		}
	}
	cerr << "No Value for keyword \"" << keyword << "\" found in " << filename;
	exit(1);
	return theString;
}

string Control::readFileName(string filename,string keyword)		//!filename: name of control file; keyword: string to be read from file; return value: (string) filname (a string including posibility for whitespaces)
{
	ifstream file;
	string s;
	file.open(filename.c_str());
	istringstream sin;
	string op;
	string theString;

	if (!file)
	{
		cerr << "Error reading boundary conditions file " << filename; exit(1);
	}

	while (file.good())
	{
		op = " ";
		getline(file, s);
		sin.clear();
		sin.str(s);
		sin >> op ;
		if (op == keyword)
		{
			string temp = s.substr(keyword.length());
			size_t found = temp.find_first_not_of(" ");
			theString = temp.substr(found);
			return theString;
		}
	}
	cerr << "No Value for keyword \"" << keyword << "\" found in " << filename;
	exit(1);
	return theString;
}

void Control::readVector(string filename,string keyword,double* vector)	//!filename: name of control file; keyword: string to be read from file; return value a vector (3D!) to corresponding keyword --> the 3 components are written to array
{
	ifstream file;
	string s;
	file.open(filename.c_str());
	istringstream sin;
	string op;
	double x, y, z;

	if (!file)
	{
		cerr << "Error reading boundary conditions file " << filename; exit(1);
	}

	while (file.good())
	{
		op = " ";
		getline(file, s);
		sin.clear();
		sin.str(s);
		sin >> op ;
		if (op == keyword)
		{
			sin >> x >> y >> z;
			vector[0] = x;
			vector[1] = y;
			vector[2] = z;
			return;
		}
	}
	cerr << "No Value for keyword \"" << keyword << "\" found in " << filename;
	exit(1);
	return;
}
