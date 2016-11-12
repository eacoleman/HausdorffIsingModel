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
IsingModel::IsingModel() {}
IsingModel::~IsingModel() {}


/* (void) setNumThreads
 *    | How many threads to use at a time.
 *  I | (int) number of threads
 */
void IsingModel::setNumThreads(int num) {
    if(nThreads < 1) return;
    nThreads = num;
}


/* (void) setNumSpins
 *    | How many spins to use.
 *  I | (int) number of spins 
 */
void IsingModel::setNumSpins(int num) {
    nSpins=num;
}


/* (void) setHausdorffDimension
 *    | Set the lattice Hausdorff dimension.
 *  I | (double) dimension to use 
 */
void IsingModel::setHausdorffDimension(double dim) {
    if(dim <= 0) return;
    hausdorffDim=dim;
}


/* (void) setCouplingConsts
 *    | Set the values of H,J in the hamiltonian 
 *  I | (double) value of H, magnetic field coupling 
 *    | (double) value of J, neighbor couplings
 */
void IsingModel::setCouplingConsts(double tH, double tJ) {
    H=tH;
    J=tJ;
}


/* (void) setTemperature
 *    | Set the temperature of the system
 *  I | (double) value of k_B * T (>0) to use 
 */
void IsingModel::setTemperature(double tkbT) {
    if (tkbT < 0) return;
    kbT=tkbT;
}


/* (vector<int>) getSpinArray 
 *    | Returns an array of the spins (+1,-1, or 0)
 */
std::vector<int> IsingModel::getSpinArray() {
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
std::vector<int> IsingModel::getLatticeDimensions() {
    return latticeDimensions;
}


/* (int) getMagnetization() 
 *    | Returns the magnetization of the lattice 
 */
int IsingModel::getMagnetization() {
    int mag=0;
    for(const auto &it : spinArray) {
       mag += it.S*it.active;
    }
    return mag;
}


/* (double) getFreeEnergy 
 *    | Get the free energy of the system (no multithread)
 *  I | (vector<int> (default: empty)) array of spin indices to flip 
 */
double IsingModel::getFreeEnergy(const std::vector<int>& flips) {
    double energy=0;
    for(int i=0; i<nSpins; i++) {
        spin s=spinArray.at(i);
        if (!s.active) continue;
        bool spinFlip=
            (std::find(flips.begin(),flips.end(),i)!=flips.end() ? -1 : 1); 

        energy += h()*s.S*spinFlip;

        for(int j=0; j<nSpins; j++) {
            spin ts=spinArray.at(j);
            if(!ts.active) continue;
            if(i==j)       continue;
            bool tspinFlips=
                (std::find(flips.begin(),flips.end(),j)!=flips.end() ? -1 : 1); 

            energy += K()*pow(getDistanceSq(s,ts),interactionSigma/2)
                *s.S*ts.S*spinFlip*tspinFlip;
        }
    }
    return energy;
}


/* (double) getFreeEnergy 
 *    | Get the partition function of the system (no multithread)
 *  I | (int (default: 0)) index to start the trace at 
 *    | (vector<int> (default: empty)) array of spin indices to flip 
 */
double IsingModel::computePartitionFunction(int start, std::vector<int> flips) {
    double Z=0;
    
    std::vector<int> newflips = flips; 
    newflips.push_back(start);

    // add e^-BH with S_i = +1 and S_i = -1
    Z+=exp(getFreeEnergy(spinArray,flips));
    Z+=exp(getFreeEnergy(spinArray,newflips));

    // branch into computations with history 
    // S_i = +1 and S_i = -1
    Z+=computePartitionFunction(start+1,newflips);
    Z+=computePartitionFunction(start+1,flips);

    return Z; 
}


void IsingModel::setup() {

}


void IsingModel::runMonteCarlo() {

}
