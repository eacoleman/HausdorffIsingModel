/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * IsingModel.cpp                                                              *
 * Author: Evan Coleman and Brad Marston, 2016                                 *
 *                                                                             *
 * Class definitions for Ising model simulation. Key characteristics:          *
 *  - Generates lattice of a given Hausdorff dimension and depth               *
 *  - Multithreaded computation of partition function                          *
 *  - Multithreaded Monte Carlo steps                                          *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "interface/IsingModel.h"

// Constructors/destructors implemented simply
// (because of number of options)
IsingModel::~IsingModel() {};


/* (void) setNumThreads
 *    | How many threads to use at a time.
 *  I | (int) number of threads
 */
void IsingModel::setNumThreads(const int num) {
    if(nThreads < 1) return;
    nThreads = num;
    hasBeenSetup=false;
}


/* (void) setNumMCSteps 
 *    | How many MC steps to perform 
 *  I | (int) number of steps 
 */
void IsingModel::setNumMCSteps(const int num) {
    if(num < 1) return;
    nMCSteps = num;
    hasBeenSetup=false;
}


/* (void) setLatticeDepth
 *    | How many steps to simulate into the fractal lattice 
 *  I | (int) depth 
 */
void IsingModel::setLatticeDepth(const int num) {
    latticeDepth=num;
    hasBeenSetup=false;
}


/* (void) setInteractionSigma
 *    | The distance coupling exponent on the J_ij 
 *  I | (int) number of spins 
 */
void IsingModel::setInteractionSigma(const double sig) {
    interactionSigma=sig;
    hasBeenSetup=false;
}


/* (void) setHausdorffDimension
 *    | Set the lattice Hausdorff dimension.
 *  I | (double) dimension to use 
 */
void IsingModel::setHausdorffDimension(const double dim) {
    if(dim <= 0) return;
    hausdorffDim=dim;
    hasBeenSetup=false;
}


/* (void) setHausdorffMethod
 *    | Set the lattice Hausdorff scaling method.
 *  I | (double) method to use:
 *    |         - SCALING   = modify separation
 *    |         - SPLITTING = modify # divisions
 */
void IsingModel::setHausdorffMethod(char* const hmtd) {
    hausdorffMethod=hmtd;
    hasBeenSetup=false;
}


/* (void) setCouplingConsts
 *    | Set the values of H,J in the hamiltonian 
 *  I | (double) value of H, magnetic field coupling 
 *    | (double) value of J, neighbor couplings
 */
void IsingModel::setCouplingConsts(const double tH, const double tJ) {
    H=tH;
    J=tJ;
    hasBeenSetup=false;
}


/* (void) setTemperature
 *    | Set the temperature of the system
 *  I | (double) value of k_B * T (>0) to use 
 */
void IsingModel::setTemperature(const double tkbT) {
    if (tkbT < 0) return;
    kbT=tkbT;
    hasBeenSetup=false;
}


/* (vector<int>) getSpinArray 
 *    | Returns an array of the spins (+1,-1, or 0)
 */
const std::vector<int> IsingModel::getSpinArray() {
    std::vector<int> spins(0);
    for(const auto &it : spinArray) {
        spins.push_back(it.active ? it.S : 0);
    }
    return spins;
}


/* (vector<int>) getLatticeDimensions 
 *    | Returns an array of the number of spins along
 *    | each lattice edge
 */
const std::vector<int> IsingModel::getLatticeDimensions() {
    return latticeDimensions;
}


/* (int) getMagnetization() 
 *    | Returns the magnetization of the lattice 
 */
const int IsingModel::getMagnetization() {
    int mag=0;
    for(const auto &it : spinArray) {
       mag += it.S*it.active;
    }
    return mag;
}


/* (double) getFreeEnergy 
 *    | Get the free energy of the system (no multithread)
 *  I | (int) single spin to flip 
 */
const double IsingModel::getFreeEnergy(const int flip) {
    std::vector<int> tflips;
    tflips.push_back(flip);
    return getFreeEnergy(tflips);
}


/* (double) getFreeEnergy 
 *    | Get the free energy of the system (no multithread)
 *  I | (vector<int> (default: empty)) array of spin indices to flip 
 */
const double IsingModel::getFreeEnergy(const std::vector<int>& flips) {
    double energy=0;
    for(int i=0; i<nLatticePoints; i++) {
        spin s=spinArray.at(i);
        if (!s.active) continue;
        bool spinFlip=
            (std::find(flips.begin(),flips.end(),i)!=flips.end() ? -1 : 1); 

        energy += h()*s.S*spinFlip;

        for(int j=0; j<nLatticePoints; j++) {
            spin ts=spinArray.at(j);
            if(!ts.active) continue;
            if(i==j)       continue;
            bool tspinFlip=
                (std::find(flips.begin(),flips.end(),j)!=flips.end() ? -1 : 1); 

            energy += K()*pow(getDistanceSq(s,ts),interactionSigma/2)
                *s.S*ts.S*spinFlip*tspinFlip;
        }
    }
    return energy;
}


/* (double) computePartitionFunction 
 *    | Get the partition function of the system (no multithread)
 *  I | (int (default: 0)) index to start the trace at 
 *    | (vector<int> (default: empty)) array of spin indices to flip 
 */
const double IsingModel::computePartitionFunction(const int start, const std::vector<int>& flips) {
    double Z=0;
    
    std::vector<int> newflips = flips; 
    newflips.push_back(start);

    // add e^-BH with S_i = +1 and S_i = -1
    Z+=exp(getFreeEnergy(flips));
    Z+=exp(getFreeEnergy(newflips));

    // branch into computations with history 
    // S_i = +1 and S_i = -1
    Z+=computePartitionFunction(start+1,newflips);
    Z+=computePartitionFunction(start+1,flips);

    return Z; 
}


/* (void) addSpins 
 *    | Recursively add spins to the spinArray 
 *  I | (int) depth to build to
 *    | (double) vector of coordinates to start fractal at
 *    | (double) vector of coordinates to end fractal at
 */
void IsingModel::addSpins(const int depth, const std::vector<double>& x0, const std::vector<double>& x1) {

    // Each iteration will add the spins between the endpoints
    // TODO: More development. Current template is slow. 

    // i is the dimension we are building spins in 
    // j is the index of the hausdorff slice we are in
    // k is an index to construct other coordinates of the sites
    for(int i=0; i <  latticeDimensions.size(); i++) {
    for(int k=0; k <  latticeDimensions.size(); k++) {
    for(int j=1; j <= hausdorffSlices; j++) {
        double delta   = abs(x1.at(i)-x0.at(i));
        double slices  = hausdorffSlices*hausdorffScale*delta;
        double spacing = (delta - slices)/(hausdorffSlices-1);
        double newstart = x0.at(i) + (j-1)*hausdorffScale*delta + j*spacing; 
        double newend   = newstart + hausdorffScale*delta;

        // j is the index of new lattice points on slice i
        spin ts;
        ts.active=1;
        ts.S=1;
        ts.coords.push_back(newstart);
        
        spinArray.push_back(ts);
        
        if(depth == 0) {
            return;
        } else {
            addSpins(depth-1,newstart,newend);
        }
    }}}

}


/* (void) setup 
 *    | Prepare the class object for simulation 
 */
void IsingModel::setup() {
    
    // Calculate the lattice dimensions from the input
    // Hausdorff dimension
    // (See my notes for derivation)
    if(hausdorffMethod=="SCALING") {
        hausdorffScale=pow(hausdorffSlices,-1/hausdorffDim);
    }

    // Check that the dimensions are reasoable
    if(hausdorffScale > 1/hausdorffSlices) {
        std::cout<<"ERROR: Invalid Hausdorff scaling"<<std::endl;
        exit(EXIT_FAILURE); 
    }
    
    // Keep track of the number of site coordinates along each axis
    for(int i=0; i<ceil(hausdorffDim); i++) {
        latticeDimensions.push_back(0);
    }

    // Keep track of the number of lattice points,
    // this is the length of the spinArray
    nLatticePoints=1;
    for(int i : latticeDimensions) nLatticePoints *= i;

    // Generate the lattice array
    std::vector<double> x0;
    std::vector<double> x1;
    for(int i=0; i<latticeDimensions.size(); i++) {
        x0.push_back(0);
        x1.push_back(1);
    }

    addSpins(latticeDepth,x0,x1);

    hasBeenSetup=true;
}


/* (void) runMonteCarlo 
 *    | Run the Monte Carlo simulation (spin-flipping) 
 */
void IsingModel::runMonteCarlo() {
    if(!hasBeenSetup) {
        std::cout<<"ERROR: Object has not been setup!"<<std::endl;
        exit(EXIT_FAILURE); 
    }

    // TODO: Multithread
    for(int i=0; i < nMCSteps; i++) {
        monteCarloStep();
    }
}


/* (void) monteCarloStep 
 *    | Perform one run over the lattice, attempting spin-flips 
 */
void IsingModel::monteCarloStep() {

    // loop over spins
    for(int i=0; i < nSpins; i++) {
        double tE = getFreeEnergy(i);
        bool spinFlip=false;

        if(tE-freeEnergy<0) {
            spinFlip = true;
        } else {
            // check probability to flip spin 
        }
            
        if(spinFlip) {
            spinArray.at(i).S=-spinArray.at(i).S;
            // TODO: keep track of spin flips?
        }
    }

}


/* (void) simpleMonteCarloStep 
 *    | Perform one run over the lattice, attempting spin-flips 
 *    | SIMPLE: Nearest-neighbors only!
 */
void IsingModel::simpleMonteCarloStep() {

}
