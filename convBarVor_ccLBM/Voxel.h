#pragma once

#include "Control.h"
#include "Vector3D.H"

using namespace std;

class Voxel
{
public:
	//!constructors
	Voxel() {} //!standard

	Voxel(double x, double y, double z, int indexI, int indexJ, int indexK, int level) : x_(x), y_(y), z_(z), indexI_(indexI), indexJ_(indexJ), indexK_(indexK), level_(level) //!init all data
	{
		for (long long dir = 0; dir < 19; dir++)
		{
			distr_[false][dir] = 0;
			distr_[true][dir] = 0;

			distrPreCol_[dir] = 0;
		}
		turbVisc_ = 0;
		isCorner_ = false;
	}

	//!destructor
	virtual ~Voxel()
	{

	}

	string tag_;

	//!get and set methods
	virtual inline int    getIndexI() const { return indexI_; }
	virtual inline int    getIndexJ() const { return indexJ_; }
	virtual inline int    getIndexK() const { return indexK_; }
	virtual inline int    getLevel() const { return level_; }
    virtual inline int    getNormalDir() { return normalDir_; }
    virtual inline bool   getIsCorner() { return isCorner_; }
    virtual inline double getTimeStep() const { return timeStep_; }
	virtual inline double getX() const { return x_; }
	virtual inline double getY() const { return y_; }
	virtual inline double getZ() const { return z_; }
	virtual inline double getCollisionFrequency() const { return colFreq_; }
	virtual inline double getTurbulentViscosity() const { return turbVisc_; }
	virtual inline double getDistribution(bool& distr, int& dir) { if (dir < 19) return distr_[distr][dir]; else return NULL; };
	virtual inline double getDistributionPreCol(int& dir) { if (dir < 19) return distrPreCol_[dir]; else return NULL; };
	virtual inline double getEquilibriumDistribution(int& dir, double* equil) { return equil[dir]; };
	virtual inline double getRecursiveRegularizedOffEquilibrium(int& dir, double* nEquilRR) { return nEquilRR[dir]; };
	virtual inline double getPi1Tensor(int i, int j) { return Pi1Tensor_[i][j]; }
	virtual inline double getCubicMachCorrection(int& dir) { if (dir < 19) return psi_[dir]; else return NULL; };
	virtual inline double getCubicMachCorrectionPreCol(int& dir) { if (dir < 19) return psiPreCol_[dir]; else return NULL; };
	virtual inline double getStrainRateTensor(int i, int j) { return strainRateTensor_[i][j]; }
	virtual inline double getStrainRateTensorFD(int i, int j) { return strainRateTensorFD_[i][j]; }
    virtual inline string getTag() { return tag_; }
    virtual inline Voxel* getNeighbor(int dir) { if (dir < 18) return neighbors_[dir]; else return NULL; }
	virtual inline Voxel* getPartner(int idx) { return partner_[idx]; }
	virtual inline Voxel* getNormalNeighbor() { return normalNeighbor_; }
	virtual inline Voxel* getTangentNeighbor1(int dir) { return tangentNeighbor1_[dir]; }
	virtual inline Voxel* getTangentNeighbor2(int dir) { return tangentNeighbor2_[dir]; }
	virtual inline Voxel* getNormalNeighborSnd() { return normalNeighborSnd_; }

	virtual inline void	setDistributionPreCol(int& dir, double distrVal) { if (dir < 19) distrPreCol_[dir] = distrVal; }
	virtual inline void setPi1Tensor(int i, int j, double value) { Pi1Tensor_[i][j] = value; }
	virtual inline void	setCubicMachCorrection(int& dir, double psiVal) { if (dir < 19) psi_[dir] = psiVal; }
	virtual inline void	setCubicMachCorrectionPreCol(int& dir, double psiVal) { if (dir < 19) psiPreCol_[dir] = psiVal; }
	virtual inline void setStrainRateTensor(int i, int j, double value) { strainRateTensor_[i][j] = value; }
	virtual inline void setStrainRateTensorFD(int i, int j, double value) { strainRateTensorFD_[i][j] = value; }
	virtual inline void setDistribution(bool alternating, int& dir, double distrVal) { if (dir < 19) distr_[alternating][dir] = distrVal; }
	virtual inline void setCollisionFrequency(double& colFreq) { colFreq_ = colFreq; }
	virtual inline void setTurbulentViscosity(double& turbVisc) { turbVisc_ = turbVisc; }
	virtual inline void setVolume(double& volume) { volume_ = volume; }
	virtual inline void setLevel(int& level) { level_ = level; }
	virtual inline void setNeighbor(int& dir, Voxel*& neighbor) { if (dir < 18) neighbors_[dir] = neighbor; }
	virtual inline void setPartner(Voxel*& partner, int idx) { partner_[idx] = partner; }
	virtual inline void setNormalNeighbor(Voxel*& normalNeighbor) { normalNeighbor_ = normalNeighbor; }
	virtual inline void setNormalNeighborSnd(Voxel*& normalNeighborSnd) { normalNeighborSnd_ = normalNeighborSnd; }
	virtual inline void setNormalDir(int& dir) { normalDir_ = dir; }
	virtual inline void setTangentNeighbor1(Voxel*& tangentNeighbor1, int dir) { tangentNeighbor1_[dir] = tangentNeighbor1; }
	virtual inline void setTangentNeighbor2(Voxel*& tangentNeighbor2, int dir) { tangentNeighbor2_[dir] = tangentNeighbor2; }
	virtual inline void setIsCorner(bool value) { isCorner_ = value; }
	virtual inline void setTimeStep(double& timeStep) { timeStep_ = timeStep; }

	//!members
	virtual void collideSRT(Control& bc, bool& alternating);
	virtual void collideMRT_GSbasis(Control& bc, bool& alternating);
	virtual void collideMRT_RMbasis(Control& bc, bool& alternating);
	virtual void collideRR(Control& bc, bool& distr);
	virtual void collideHRR(Control& bc, bool& distr);
	virtual void calcVelocity(Control& bc, bool& distr, double& dens, double* velo);
	virtual void calcVelocityPreCol(Control& bc, double& dens, double* velo);
	virtual void calcVelocity(Control& bc, double& dens, double* velo, double* distrArr, double* psiArr);
    virtual void calcEquilibrium2_GHbasis(double* equil, Control& bc, double& dens, double* velo);
    virtual void calcEquilibrium3_GHbasis(double* equil, Control& bc, double& dens, double* velo);
    virtual void calcEquilibrium4_GHbasis(double* equil, Control& bc, double& dens, double* velo);
    virtual void calcRRoffEquilibrium2_GHbasis(bool& alternating, double* nEquilRR, double* equil, Control& bc, double* velo);
    virtual void calcRRoffEquilibrium3_GHbasis(bool& alternating, double* nEquilRR, double* equil, Control& bc, double* velo);
    virtual void calcRRoffEquilibrium4_GHbasis(bool& alternating, double* nEquilRR, double* equil, Control& bc, double* velo);
    virtual void calcHRRoffEquilibrium2_GHbasis(bool& alternating, double* nEquilHRR, double* equil, Control& bc, double* velo);
    virtual void calcHRRoffEquilibrium3_GHbasis(bool& alternating, double* nEquilHRR, double* equil, Control& bc, double* velo);
    virtual void calcHRRoffEquilibrium4_GHbasis(bool& alternating, double* nEquilHRR, double* equil, Control& bc, double* velo);
	virtual void calcPi1Tensor(bool& alternating, double* equil, Control& bc, double* velo);
	virtual void transformMomentsToDistributions_RMbasis(double* distributions, double* moments, Control& bC);
	virtual void calcAllMoments_GSbasis(Control& bc, double* moments, bool& alternating);
	virtual void calcAllMoments_GSbasis(Control& bc, double* distr, double* moments, bool& alternating);
	virtual void calcAllMoments_RMbasis(Control& bc, double* moments, double* distribution) const;
	virtual void calcEquilibrium4_RMbasis(double* distributions, double density, double* velo, Control& bc);
	virtual void calcEquilibrium4_RMbasis(double* distributions, double density, Vector3D* velo, Control& bc);
	virtual void calcAllMomentsPreCol_GSbasis(Control& bc, double* moments);
	virtual void calcEquilibriumMoments2_GSbasis(double* equilibriumMoments, double& density, double& momentumX, double& momentumY, double& momentumZ);
	virtual void transformMomentsToDistributions_GSbasis(double* distributions, double* moments, Control& bc);
	virtual void calcStrainRateTensorFDandCubicMachCorrection(Control& bc);
	virtual double calcDensity(Control& bc, bool& alternating);
	virtual double calcDensityPreCol(Control& bc);
	virtual double calcDensity(Control& bc, double* distrArr, double* psiArr);
	virtual double calcPressure(bool& alternating, Control& bc);
	virtual inline void transport(bool& /*distr*/) { cerr << "\nTransport from voxel invoked." << flush; exit(1); } //!salvation implementation
	virtual inline void transportPull(bool& /*distr*/) { cerr << "\nTransport from voxel invoked." << flush; exit(1); } //!salvation implementation

	//!bounce-back methods
	virtual inline void halfWayBounceBack(bool&/*alternating*/, Control&/*bc*/) {}; //!virtual declaration for wallnode

	//!new members from Interface
	virtual inline void explode(bool&/*distributionFine*/, bool& /*distributionCoarse*/) { cerr << "\nExplode from voxel invoked." << flush; exit(1); };
	virtual inline void coalesce(bool&/*distributionFine*/, bool& /*distributionCoarse*/, string& /*interface*/) { cerr << "\nCoalescence from voxel invoked." << flush; exit(1); };
	virtual inline void explodeLinear(bool&/*distributionCoarse*/, bool& /*distributionFine*/, Control& /*bc*/, string& /*interface*/) {};
	virtual inline void calcGhostDistribution(bool& /*alternatingL1*/, double* /*ghostArr*/) {};


private:

	//!coordinates
	double x_;
	double y_;
	double z_;

protected:

	//!level
	int level_;

	//!voxel indices
	int indexI_;
	int indexJ_;
	int indexK_;

	//!distributions
	double distr_[2][19];

	//!deviatoric stress tensor
	double Pi1Tensor_[3][3];

	//!cubic Mach correction terms
	double psi_[19];

	//!cubic Mach correction terms
	double psiPreCol_[19];

	//!strain-rate tensor (second order central finite difference scheme)
	double strainRateTensor_[3][3];

	//!strain-rate tensor (second order central finite difference scheme)
	double strainRateTensorFD_[3][3];

	//!pre-collision distribution
	double distrPreCol_[19];

	//!relaxation time
	double colFreq_;

	//!eddy viscosity
	double turbVisc_;

	//!time step
	double timeStep_;

	//!array of pointers to adjacent nodes --> create connectivity D3Q19
	Voxel* neighbors_[18];

	//!partner voxels
	Voxel* partner_[8];

	double volume_;

	//!first and second normal-direction neighbors of boundary nodes for LODI pressure BC
	Voxel* normalNeighbor_;
	Voxel* normalNeighborSnd_;
	Voxel* tangentNeighbor1_[2];
	Voxel* tangentNeighbor2_[2];
	int normalDir_;
	bool isCorner_;

};
