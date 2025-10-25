#pragma once

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>

using namespace std;

class Control
{
	public:

	//!constructor
	Control() {}

	Control(string filename);

	//!destructor
	~Control() {}

	//!members: read functions for parameters from control file
	double readValue(string filename,string keyword); //!filename: name of control file; keyword: string to be read from file; return value: (double) to corresponding keyword
	string readString(string filename,string keyword); //!filename: name of control file; keyword: string to be read from file; return value: (string) to corresponding keyword
	string readFileName(string filename,string keyword); //!filename: name of control file; keyword: string to be read from file; return value: (string) filname (a string including posibility for whitespaces)
	void readVector(string filename,string keyword,double* vector);	//!filename: name of control file; keyword: string to be read from file; return value a vector (3D!) to corresponding keyword --> the 3 components are written to array

	//!get and set members
	inline int			  getNumberOfVoxelsX() { return numberOfVoxelsX_; }
	inline int		      getNumberOfVoxelsY() { return numberOfVoxelsY_; }
	inline int			  getNumberOfVoxelsZ() { return numberOfVoxelsZ_; }
	inline int			  getNumberOfNodesX() { return numberOfNodesX_; }
	inline int		      getNumberOfNodesY() { return numberOfNodesY_; }
	inline int			  getNumberOfNodesZ() { return numberOfNodesZ_; }
	inline int            getLevelMax() { return levelMax_; }
	inline int			  getTimeStepMax() { return timeStepMax_; }
	inline int			  getWriteInterval() { return writeInterval_; }
	inline int			  getNumberOfThreads() { return numberOfThreads_; }
	inline int            getExplosionOrder() { return explosionOrder_; }
	inline int		      getRefinementLayer() { return refLayer_; }
    inline int		      getRefinementLayerX() { return refLayerX_; }
    inline int            getMPnum() {return MPnum_;}
    inline int            getOrderOfEquilibrium() { return equilOrder_; }
    inline void           setTime_(double time) { time_ = time; }
    inline bool           getTurbulenceModelling() { return turbulenceModelling_; }
	inline double		  getSpacing() {return spacing_;}
	inline double		  getConvCriterion() { return convCriterion_; }
	inline double		  getHybridParameter() { return hybridPar_; }
	inline double		  getDensity() {return density_;}
    inline double		  getMPradius() { return MPradius_; }
    inline double		  getMPZ() { return MPZ_; }
	inline double		  getMachNumber() {return machNumber_;}
	inline double		  getReynoldsNumber() { return reynoldsNumber_; }
	inline double		  getTime_() { return time_; }
	inline double		  getMolecularVelocity() {return molecularVelocity_;}
	inline double		  getSpeedOfSound() {return speedOfSound_;}
	inline double		  getKinematicViscosity() {return kinematicViscosity_;}
	inline double         getSmagorinsky() {return smagorinskyConstant_;}
	inline double         getTimeStep() {return timeStep_;}
	inline double		  getChannelSizeX_() {return channelSizeX_;}
	inline double		  getChannelSizeY_() {return channelSizeY_;}
	inline double		  getChannelSizeZ_() {return channelSizeZ_;}
    inline double*        getInitVelo() { return initVelocity_; }
	inline double*        getXsiX() {return xsiX_;} //!unsafe: reference to local variable is passed; watch out when using
	inline double*        getXsiY() {return xsiY_;} //!unsafe: reference to local variable is passed; watch out when using
	inline double*        getXsiZ() {return xsiZ_;} //!unsafe: reference to local variable is passed; watch out when using
	inline double*        getWeightFactors() {return weightFactors_;} //!unsafe: reference to loacl variable is passed; watch out when using
	virtual inline double getCollisionFrequency() { return collisionFrequency_; }
	inline string         getCaseName() { return caseName_; }
	inline string		  getCollisionModel() { return collisionModel_; }
	inline string		  getCubicMachCorrection() { return cubicMachCorrection_; }

	private:

	//!grid variables
	int    levelMax_; //!max. level of mesh
	int    numberOfVoxelsX_, numberOfVoxelsY_, numberOfVoxelsZ_; //!number of voxels in x, y and z direction
	int    numberOfNodesX_, numberOfNodesY_, numberOfNodesZ_; //!number of nodes in x, y and z direction
	int	   refLayer_; //!number of refined coarse voxels at boundaries
    int	   refLayerX_; //!number of refined coarse cells at boundaries
	double spacing_; //!mesh spacing in level 0
	double convCriterion_; //!convergence criterion
	double channelSizeX_, channelSizeY_, channelSizeZ_; //!channel dimensions in x, y and z direction

	//!hybrid recursive-regularization
	double hybridPar_; //!hybridization parameter

	//!order of accuracy of explosion
	int explosionOrder_;
	int equilOrder_;

	//!lattice variables
	double xsiX_[19], xsiY_[19], xsiZ_[19]; //!components of molecular velocity vetor (D3Q19 lattice)
	double weightFactors_[19]; //!lattice weight factors (D3Q19)

	//!physical variables
    int MPnum_;                     //!number of monitor points
	double density_;				//!fluid density
	double machNumber_;				//!global Mach number
	double molecularVelocity_;		//!molecular velocity (delta_x/delta_t)
	double speedOfSound_;			//!speed of sound
	double kinematicViscosity_;		//!kinematic velocity
    double initVelocity_[3];		//!prescribed initialization velocity
    double MPradius_;               //!radial distance from center axis to monitor point position
    double MPZ_;                    //!axial position of monitor points
	double timeStep_;				//!time step in level 0
	double collisionFrequency_;	    //!dimensionless collision frequency on level 0
	double reynoldsNumber_;		    //!reynolds number

	//!turbulence parameters (needed later)
	bool   turbulenceModelling_;    //!if true LES is performed
	double smagorinskyConstant_;	//!smagorinsky constant

	//!control parameter
	int    timeStepMax_; //!number of time steps to run(in level 0)
	int    writeInterval_; //!number of time steps between to consecutive result files
	int    numberOfThreads_; //!number of cpus for SMP
    double time_; //!physical time
	string caseName_; //!absolute path to the file to be saved
	string collisionModel_;
	string cubicMachCorrection_;
};
