/* D3Q19 LBM SIMULATION OF CONVECTED BAROPTROPIC VORTEX */
//!SEE "Schukmann, A.; Haas, V.; Schneider, A. Spurious Aeroacoustic Emissions in Lattice Boltzmann Simulations on Non-Uniform
//!Grids. Fluids 2025, 10, 31. https://doi.org/10.3390/fluids10020031"
//!AND "Schukmann, A.; Schneider, A.; Haas, V.; Böhle, M. Analysis of Hierarchical Grid Refinement Techniques for the Lattice
//!Boltzmann Method by Numerical Experiments. Fluids 2023, 8, 103. https://doi.org/10.3390/fluids8030103 " FOR FURTHER DETAILS

#include <sstream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <chrono>
#include <ctime>
#include <vector>
#include "Control.h"
#include "Lattice.h"
#include "Node.h"
#include "Voxel.h"

using namespace std;
using namespace chrono;

ostream&
display(ostream& os, nanoseconds ns)
{
    typedef duration<int, ratio<86400>> days;
    char fill = os.fill();
    os.fill('0');
    auto d = duration_cast<days>(ns);
    ns -= d;
    auto h = duration_cast<hours>(ns);
    ns -= h;
    auto m = duration_cast<minutes>(ns);
    ns -= m;
    auto s = duration_cast<seconds>(ns);
    os << setw(2) << d.count() << "d:"
        << setw(2) << h.count() << "h:"
        << setw(2) << m.count() << "m:"
        << setw(2) << s.count() << 's' << endl << endl;
    os.fill(fill);
    return os;
};

int main(int argc, char* argv[])
{
    steady_clock::time_point beginTotal = steady_clock::now();
    steady_clock::time_point endTotal;
    steady_clock::time_point beginCalc;
    steady_clock::time_point endCalc;
    nanoseconds calcTime;
    string inputFile = "/absolute-path-to-BC-file/BoundaryConditions.txt"; //!name of input file; could also be replaced by argv for dynamic file names
    string initTimeStepString;
    Control boundaryConditions(inputFile); //!create instance of control class and read parameters from file
    Lattice simulationCase(boundaryConditions); //!create instance of lattice
    int time = 0;
    boundaryConditions.setTime_(static_cast<double>(time));
    int initTime = 0;
    int writeInterval = boundaryConditions.getWriteInterval(); //!get write interval
    int maxTimeStep = boundaryConditions.getTimeStepMax();
    int maxLevel = boundaryConditions.getLevelMax();
    bool alternatingL0 = true; //!variable for switching between distribution arrays
    bool negAlternatingL0 = !alternatingL0;
    bool alternatingL1 = true;

    //!write to console
    cout << "LB convected barotropic vortex solver" << endl;
    cout << "initializing from uniform zero-velocity equilibrium distribution" << endl << endl;

    if (maxLevel == 0) //!single level
    {
        cout << "all parameters are being displayed in SI units" << endl << endl;
        cout << "        dimensionless collision frequency = " << boundaryConditions.getCollisionFrequency() << endl;
        cout << "                           time step size = " << boundaryConditions.getTimeStep() << endl;
        cout << "                          REYNOLDS number = " << boundaryConditions.getReynoldsNumber() << endl;
        cout << "                      convective velocity = " << boundaryConditions.getInitVelo()[0] << endl;
        cout << "                      kinematic viscosity = " << boundaryConditions.getKinematicViscosity() << endl;
        cout << "                              MACH number = " << boundaryConditions.getMachNumber() << endl;
        cout << "                     number of time steps = " << boundaryConditions.getTimeStepMax() << endl;
        cout << "                          collision model = " << boundaryConditions.getCollisionModel() << endl << endl;

        //!initialize nodes
        simulationCase.initialization(boundaryConditions);

        //!set calculation clock
        beginCalc = steady_clock::now();

        //!print current time to console
        auto startTime = system_clock::now();
        time_t cur_time = system_clock::to_time_t(startTime);

        cout << "current date and time is: " << ctime(&cur_time) << endl;

        while (time < maxTimeStep) //!loop over all time steps
        {
            simulationCase.solve(alternatingL0); //!solve case for this time step

            time++; //!increment time
            boundaryConditions.setTime_(static_cast<double>(time));

            if (time - initTime == 1)
            {
                endCalc = steady_clock::now();
                calcTime = duration_cast<nanoseconds>((endCalc - beginCalc) * (maxTimeStep - initTime));

                cout << "estimated total computation time: ";
                display(cout, calcTime);
            } //!end if
            else if ((time % writeInterval) == 0) //!write results every write interval
            {
                //!write result file
                string file = boundaryConditions.getCaseName(); //!get case name
                ostringstream timeStep; //!create string stream for conversion from integer to string
                timeStep << time; //!convert time (integer) to string
                file.append("_" + timeStep.str()); //!add current time step to file name
                simulationCase.writeResultsVTK(file, negAlternatingL0, alternatingL1); //!write result file

                //!write status to console
                cout << "solved " << fixed << setprecision(2) << (static_cast<double>(time) - static_cast<double>(initTime)) / (maxTimeStep - static_cast<double>(initTime)) * 100 << "% of the case" << endl;

                if (time - initTime == writeInterval)
                {
                    endCalc = steady_clock::now();
                    calcTime = duration_cast<nanoseconds>((endCalc - beginCalc) * ((maxTimeStep - initTime) / writeInterval));

                    cout << endl << "corrected computation time estimation (including time for file creation): ";
                    display(cout, calcTime);
                } //!end if
            } //!end else if

            alternatingL0 = !alternatingL0; //!switch distribution array
            negAlternatingL0 = !negAlternatingL0;
        } //!end time
    } //!end if
    else //!multi-level
    {
        cout << "all parameters are being displayed in SI units" << endl << endl;
        cout << "result file(s): " << boundaryConditions.getCaseName() << "_t.vtk" << endl << endl;
        cout << "dimensionless collision frequency on coarse grid = " << boundaryConditions.getCollisionFrequency() << endl;
        cout << "  dimensionless collision frequency on fine grid = " << 1. / ((2. / boundaryConditions.getCollisionFrequency()) - 0.5) << endl;
        cout << "                   time step size on coarse grid = " << boundaryConditions.getTimeStep() << endl;
        cout << "                     time step size on fine grid = " << boundaryConditions.getTimeStep() / 2. << endl;
        cout << "                          spacing on coarse grid = " << boundaryConditions.getSpacing() << endl;
        cout << "                            spacing on fine grid = " << boundaryConditions.getSpacing() / 2. << endl;
        cout << "                                 REYNOLDS number = " << boundaryConditions.getReynoldsNumber() << endl;
        cout << "                             convective velocity = " << boundaryConditions.getInitVelo()[0] << endl;
        cout << "                             kinematic viscosity = " << boundaryConditions.getKinematicViscosity() << endl;
        cout << "                                         density = " << boundaryConditions.getDensity() << endl;
        cout << "                                     MACH number = " << boundaryConditions.getMachNumber() << endl;
        cout << "                             order of explosion  = " << boundaryConditions.getExplosionOrder() << endl;
        cout << "                            order of equilibrium = " << boundaryConditions.getOrderOfEquilibrium() << endl;
        cout << "                                 collision model = " << boundaryConditions.getCollisionModel() << endl;
        if(boundaryConditions.getCollisionModel() == "HRR")
        {
            cout << "                         hybridization parameter = " << boundaryConditions.getHybridParameter() << endl;
        }
        cout << "                                   domain size x = " << boundaryConditions.getChannelSizeX_() << endl;
        cout << "                                   domain size y = " << boundaryConditions.getChannelSizeY_() << endl;
        cout << "                                   domain size z = " << boundaryConditions.getChannelSizeZ_() << endl;
        cout << "                                cubic correction = " << boundaryConditions.getCubicMachCorrection() << endl;
        cout << "                            number of time steps = " << boundaryConditions.getTimeStepMax() << endl;
        cout << "                                 write intervall = " << boundaryConditions.getWriteInterval() << endl;
        cout << "refinement layers at Y and Z boundaries and back = " << boundaryConditions.getRefinementLayer() << endl;
        cout << "                    refinement layers at X front = " << boundaryConditions.getRefinementLayerX() << endl << endl;

        //!initialize nodes
        simulationCase.initialization(boundaryConditions);

        //!set calculation clock
        beginCalc = steady_clock::now();

        //!print current time to console
        auto startTime = std::chrono::system_clock::now();
        time_t cur_time = system_clock::to_time_t(startTime);

        cout << "current date and time is: ";
        cout << ctime(&cur_time) << endl;

        while (time < maxTimeStep) //!loop over all time steps
        {
            if (time == 0)
            {
                //!write initialized fields
                string file = boundaryConditions.getCaseName(); //!get case name
                ostringstream timeStep; //!create string stream for conversion from integer to string
                timeStep << time; //!convert time (integer) to string
                file.append("_" + timeStep.str()); //!add current time step to file name
                simulationCase.writeResultsVTK(file, alternatingL0, alternatingL1); //!write result file
            }

            simulationCase.solve(alternatingL0, alternatingL1); //!solve case for this time step

            //!write result files at nodes
            string fileMP = boundaryConditions.getCaseName().append("_MONITORPOINTS.dat"); //!Write monitor points
            simulationCase.writeResultsAtVoxel(fileMP, negAlternatingL0, alternatingL1, boundaryConditions); //!write result file

            time++; //!increment time
            boundaryConditions.setTime_(static_cast<double>(time));

            if (time - initTime == 1)
            {
                endCalc = steady_clock::now();
                calcTime = duration_cast<nanoseconds>((endCalc - beginCalc) * (maxTimeStep - initTime));

                cout << "estimated total computation time: ";
                display(cout, calcTime);
            } //!end if
            else if ((time % writeInterval) == 0) //!write results every write interval
            {
                //!write result file
                string file = boundaryConditions.getCaseName(); //!get case name
                ostringstream timeStep; //!create string stream for conversion from integer to string
                timeStep << time; //!convert time (integer) to string
                file.append("_" + timeStep.str()); //!add current time step to file name
                simulationCase.writeResultsVTK(file, negAlternatingL0, alternatingL1); //!write result file

                //!write status to console
                cout << "solved " << fixed << setprecision(2) << (static_cast<double>(time) - static_cast<double>(initTime)) / (maxTimeStep - static_cast<double>(initTime)) * 100 << "% of the case" << endl;

                if (time - initTime == writeInterval)
                {
                    endCalc = steady_clock::now();
                    calcTime = duration_cast<nanoseconds>((endCalc - beginCalc) * ((maxTimeStep - initTime) / writeInterval));

                    cout << endl << "corrected computation time estimation (including time for file creation): ";
                    display(cout, calcTime);
                } //!end if
            } //!end else if
            alternatingL0 = !alternatingL0; //!switch distribution array
            negAlternatingL0 = !negAlternatingL0;
        } //!end time
    } //!end else

    endTotal = steady_clock::now();
    cout << endl << "given maximum number of time steps or convergence criterion reached" << endl << endl;
    cout << "total computation time: ";
    display(cout, duration_cast<nanoseconds>(endTotal - beginTotal));
} //!end main
